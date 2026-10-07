"""
PocketMapper: map and compare binding pockets across protein structures.

`PocketMapper.search()` is the whole workflow. It resolves its settings, then runs the steps in
`pocketmapper.steps` in order, each handing the next its files under `results_dir`:

1. `configure_workflow` -> job file over arguments; `resolve_settings` -> Settings, job_settings.json.
2. `steps.parse` -> query_records.json, target_records.json, cache_dirs.json; failed_entries.json
   started afresh.
3. `steps.fetch` -> structures, a bundled Foldseek database and PISA interfaces into the cache.
4. `steps.align` -> foldseek or the local sequence aligner, per `aligner` -> alignment.tsv. Against
   a PDB Foldseek database, its hits are appended to the target records as pisa records.
5. `steps.pockets` -> pockets.json. The Pocket shape itself is declared in `pockets/pocket.py`.
6. `steps.compare` -> pocket_comparison.tsv.
7. `steps.superpose` -> the top align_count targets per query, superposed into aligned_structures/.

Every step but parse drops the records it fails, into failed_entries.json, so each records file
holds only usable records. Each step also runs on its own; see `pocketmapper.commands`.

This is the only module that knows about `Settings`; components are handed the individual values
they need, so none of them has to build one to be usable on its own.

Author: Lachlan Ellingboe
"""

import json
import logging
import os
from dataclasses import asdict
from dataclasses import fields

from pocketmapper.constants import DEFAULT_ALIGN_COUNT
from pocketmapper.constants import DEFAULT_ALIGN_STRUCT_METHOD
from pocketmapper.constants import DEFAULT_ALIGNER
from pocketmapper.constants import DEFAULT_CACHE_DIR
from pocketmapper.constants import DEFAULT_DELETE_TMP
from pocketmapper.constants import DEFAULT_PISA_SOURCE
from pocketmapper.constants import DEFAULT_POCKET_METHOD
from pocketmapper.constants import DEFAULT_VERBOSITY
from pocketmapper.exceptions import PocketMapperError
from pocketmapper.lib import delete_temp_dir
from pocketmapper.lib import empty_temp_dir
from pocketmapper.lib import log_to_file
from pocketmapper.lib import make_dir
from pocketmapper.records import CACHE_MANIFEST_KEYS
from pocketmapper.records import fsdb_record
from pocketmapper.records import read_cache_manifest
from pocketmapper.records import read_records
from pocketmapper.settings import Settings
from pocketmapper.settings import check_fsdb_align_struct_method
from pocketmapper.settings import check_fsdb_aligner
from pocketmapper.settings import require_foldseek
from pocketmapper.settings import resolve_align_struct_method
from pocketmapper.settings import resolve_aligner
from pocketmapper.settings import resolve_delete_tmp
from pocketmapper.settings import resolve_paths
from pocketmapper.settings import resolve_pisa_source
from pocketmapper.settings import resolve_threads
from pocketmapper.steps.align import align_chains
from pocketmapper.steps.compare import compare_aligned_pockets
from pocketmapper.steps.fetch import fetch_inputs
from pocketmapper.steps.parse import parse_inputs
from pocketmapper.steps.pockets import build_pockets
from pocketmapper.steps.superpose import superpose_top_targets

logger = logging.getLogger(__name__)


class PocketMapper:
    """
    The pipeline. `search()` runs the whole workflow; see the module docstring for its steps.

    Holds no state between calls.
    """

    def search(
        self,
        query=None,
        target=None,
        job_file=None,
        cache_dir=DEFAULT_CACHE_DIR,
        results_dir=None,
        verbosity=DEFAULT_VERBOSITY,
        threads=None,
        aligner=DEFAULT_ALIGNER,
        align_count=DEFAULT_ALIGN_COUNT,
        align_struct_method=DEFAULT_ALIGN_STRUCT_METHOD,
        query_pocket_method=DEFAULT_POCKET_METHOD,
        target_pocket_method=DEFAULT_POCKET_METHOD,
        delete_tmp=DEFAULT_DELETE_TMP,
        pisa_source=DEFAULT_PISA_SOURCE,
        pdb_dir=None,
        alphafold_dir=None,
        pocket_dir=None,
        foldseek_preprocessed_structure_dir=None,
        temp_dir=None,
        aligned_structure_dir=None,
        alignment_path=None,
        pocket_comparison_path=None,
        query_records_path=None,
        target_records_path=None,
        pockets_path=None,
        failed_entries_path=None,
        job_settings_path=None,
        log_path=None,
        fsdb_dir=None,
    ):
        """
        Run the full PocketMapper search workflow.

        For the length of the call, sets the `pocketmapper` logger's level and adds a handler writing
        to `log_path`; both are restored on return. Empties `temp_dir` on the way in and, unless
        `delete_tmp` is 0, deletes it on the way out. Each directory is created only when something
        is first written into it. Besides the results, writes the steps' hand-off files into
        `results_dir` -- the records files, cache_dirs.json, pockets.json and failed_entries.json --
        so any step can be rerun there on its own.

        Args:
            query (str, optional): Query identifier, string or path to a list. Required here or in
                `job_file`, not both.
            target (str, optional): Target structure identifier, string or path to a list. As `query`.
            job_file (str, optional): Path to a JSON job file of Settings field name -> value. Every
                value it sets wins over the matching argument.
            cache_dir (str, optional): Directory to cache intermediate structures.
                Defaults to DEFAULT_CACHE_DIR.
            results_dir (str, optional): Directory to output results to.
                Defaults to pocketmapper_results_<YYMMDD_HHMMSS>.
            verbosity (int, optional): Control logging level. Defaults to DEFAULT_VERBOSITY.
            threads (int, optional): Cap on the cores Foldseek uses. Defaults to one per available core.
            aligner (str, optional): Chain aligner -- 'foldseek', which needs the foldseek binary, or
                'seq' for the local BLOSUM62 sequence aligner. Defaults to DEFAULT_ALIGNER.
            align_count (int, optional): Number of top targets to superpose onto each query.
                Defaults to DEFAULT_ALIGN_COUNT.
            align_struct_method (str, optional): Which transform superposes a target onto its query --
                'foldseek' for Foldseek's whole-chain fit, 'pocket' for the fit of the two pockets
                on their overlapping residues, or 'auto' (the default) for 'foldseek' with
                aligner 'foldseek' and 'pocket' with 'seq', which produces no chain transform at all.
            query_pocket_method (str, optional): Force a pocket method for every query entry --
                'pisa', 'passthrough', 'vdw', 'whole_chain' or 'foldseek_db' -- or 'auto' (the
                default) to infer it per entry from the input string.
            target_pocket_method (str, optional): As `query_pocket_method`, for the target side.
            delete_tmp (int, optional): 1 deletes temp_dir at the end of the run; 0 keeps it for
                inspection. Defaults to DEFAULT_DELETE_TMP.
            pisa_source (str, optional): Where PISA interfaces are fetched from -- 'ftp' for the EBI
                FTP server or 'api' for the paced PDBe API. Both serve the same data into the same
                cache. Defaults to DEFAULT_PISA_SOURCE.
            pdb_dir (str, optional): Cache of fetched PDB structures.
                Defaults to <cache_dir>/pdb_structures.
            alphafold_dir (str, optional): Cache of fetched AlphaFold structures.
                Defaults to <cache_dir>/alphafold_structures.
            pocket_dir (str, optional): Cache of parsed pockets. Defaults to <cache_dir>/pockets.
            foldseek_preprocessed_structure_dir (str, optional): Cache of the single-chain structures
                Foldseek is given. Defaults to <cache_dir>/foldseek_preprocessed_structures.
            temp_dir (str, optional): Per-run scratch, wiped before use and deleted at the end of the
                run. Defaults to <results_dir>/tmp.
            aligned_structure_dir (str, optional): Where the superposed structures for the top hits
                are written. Defaults to <results_dir>/aligned_structures.
            alignment_path (str, optional): Where the alignment table is written.
                Defaults to <results_dir>/alignment.tsv.
            pocket_comparison_path (str, optional): Where the pocket comparison table is written.
                Defaults to <results_dir>/pocket_comparison.tsv.
            query_records_path (str, optional): Where the query records are written.
                Defaults to <results_dir>/query_records.json.
            target_records_path (str, optional): Where the target records are written.
                Defaults to <results_dir>/target_records.json.
            pockets_path (str, optional): Where the pockets are written.
                Defaults to <results_dir>/pockets.json.
            failed_entries_path (str, optional): Where the entries dropped along the way are written.
                Defaults to <results_dir>/failed_entries.json.
            job_settings_path (str, optional): Where this run's resolved settings are dumped.
                Defaults to <results_dir>/job_settings.json.
            log_path (str, optional): Where the run log is written.
                Defaults to <results_dir>/info.log.
            fsdb_dir (str, optional): Cache of bundled Foldseek databases.
                Defaults to <cache_dir>/fsdb.

        Returns:
            dict: The resolved Settings as a dictionary. Results are written to `results_dir` --
                read pocket_comparison.tsv and alignment.tsv from the paths it names.

        Raises:
            PocketMapperError: On a bad job file, a query/target given both ways or neither way,
                or any pipeline failure.
        """
        # The Settings fields as passed to this call. Kept beside the signature it mirrors, so the
        # two cannot drift; `job_file` is absent because Settings has no such field.
        arguments = {
            "query": query,
            "target": target,
            "cache_dir": cache_dir,
            "results_dir": results_dir,
            "aligner": aligner,
            "verbosity": verbosity,
            "threads": threads,
            "align_count": align_count,
            "align_struct_method": align_struct_method,
            "query_pocket_method": query_pocket_method,
            "target_pocket_method": target_pocket_method,
            "delete_tmp": delete_tmp,
            "pisa_source": pisa_source,
            "pdb_dir": pdb_dir,
            "alphafold_dir": alphafold_dir,
            "pocket_dir": pocket_dir,
            "foldseek_preprocessed_structure_dir": foldseek_preprocessed_structure_dir,
            "temp_dir": temp_dir,
            "aligned_structure_dir": aligned_structure_dir,
            "alignment_path": alignment_path,
            "pocket_comparison_path": pocket_comparison_path,
            "query_records_path": query_records_path,
            "target_records_path": target_records_path,
            "pockets_path": pockets_path,
            "failed_entries_path": failed_entries_path,
            "job_settings_path": job_settings_path,
            "log_path": log_path,
            "fsdb_dir": fsdb_dir,
        }

        values = self.configure_workflow(job_file, arguments)

        # The only directory made up front: log_to_file opens its file handler immediately.
        # Every other directory is made by whatever first writes into it.
        make_dir(os.path.dirname(values["log_path"]), {"stage": "Configuring Settings"})
        with log_to_file(values["log_path"], values["verbosity"]):
            settings = self.resolve_settings(values)
            cache_dirs = {key: getattr(settings, key) for key in CACHE_MANIFEST_KEYS}
            parse_inputs(
                settings.query,
                settings.target,
                settings.query_pocket_method,
                settings.target_pocket_method,
                cache_dirs,
                settings.results_dir,
                settings.query_records_path,
                settings.target_records_path,
                settings.failed_entries_path,
            )
            # Checked before the database is downloaded
            if fsdb_record(read_records(settings.target_records_path)) is not None:
                check_fsdb_aligner(settings.aligner)
                check_fsdb_align_struct_method(settings.align_struct_method)
            # The manifest is the parse step's absolute copy of cache_dirs
            cache_dirs = read_cache_manifest(settings.results_dir)

            fetch_inputs(
                settings.query_records_path,
                settings.target_records_path,
                settings.query_records_path,
                settings.target_records_path,
                settings.failed_entries_path,
                cache_dirs["pocket_dir"],
                settings.pisa_source,
                settings.threads,
                settings.temp_dir,
            )
            align_chains(
                settings.query_records_path,
                settings.target_records_path,
                settings.query_records_path,
                settings.target_records_path,
                settings.alignment_path,
                settings.failed_entries_path,
                settings.aligner,
                settings.threads,
                settings.verbosity,
                settings.temp_dir,
                cache_dirs,
                settings.pisa_source,
            )
            # The fetch step already tried every PISA entry, and only successes are cached
            build_pockets(
                [settings.query_records_path, settings.target_records_path],
                settings.pockets_path,
                settings.failed_entries_path,
                cache_dirs["pocket_dir"],
                settings.pisa_source,
                download_pisa=False,
            )
            compare_aligned_pockets(
                settings.query_records_path,
                settings.target_records_path,
                settings.alignment_path,
                settings.pockets_path,
                settings.pocket_comparison_path,
                settings.results_dir,
            )
            superpose_top_targets(
                settings.query_records_path,
                settings.target_records_path,
                settings.pocket_comparison_path,
                settings.alignment_path,
                settings.aligned_structure_dir,
                settings.align_struct_method,
                settings.align_count,
                settings.threads,
            )
            delete_temp_dir(settings.temp_dir, settings.delete_tmp, [settings.cache_dir, settings.results_dir])

            logger.info("PocketMapper search completed successfully.", extra={"stage": "End"})

        return asdict(settings)

    def configure_workflow(self, job_file, arguments):
        """
        Layer the job file over the arguments and fill in every path left unset.

        Args:
            job_file (str or None): Path to a JSON job file, or None for none. Any value it sets wins
                over `arguments`.
            arguments (dict): Settings field name -> value from `search()`.

        Returns:
            dict: Settings field name -> value, every path set but nothing else checked.

        Raises:
            PocketMapperError: If the job file is missing, unreadable or names an unknown setting, or if
                query or target is given both ways or neither way.
        """
        log_extra = {"stage": "Configuring Settings"}

        # 1. The job file, if any
        job = {}
        if job_file is not None:
            if not os.path.isfile(job_file):
                logger.critical(f"Job file not found: {job_file}", extra=log_extra)
                raise PocketMapperError(f"Job file not found: {job_file}")
            try:
                with open(job_file) as f:
                    job = json.load(f)
            except Exception as e:
                logger.critical(f"Error reading job file: {job_file}. Is it in JSON format?", extra=log_extra)
                raise PocketMapperError(f"Error reading job file: {job_file}. Is it in JSON format?") from e
            if not isinstance(job, dict):
                msg = f'Job file {job_file} must hold a JSON object of {{"option": value}}.'
                logger.critical(msg, extra=log_extra)
                raise PocketMapperError(msg)
            unknown = sorted(set(job) - {f.name for f in fields(Settings)})
            if unknown:
                msg = f"Unknown setting(s) in {job_file}: {', '.join(unknown)}"
                logger.critical(msg, extra=log_extra)
                raise PocketMapperError(msg)

        # 2. The job file wins over the arguments, except that query and target must come from exactly
        # one of the two: a silently discarded positional would search something other than what the
        # command line shows.
        for key in ("query", "target"):
            if job.get(key) is not None and arguments[key] is not None:
                msg = f"{key} is set both in the job file and as an argument; give it only once."
                logger.critical(msg, extra=log_extra)
                raise PocketMapperError(msg)
        values = {**arguments, **job}
        for key in ("query", "target"):
            if values[key] is None:
                msg = f"No {key} given; pass it as an argument or set it in the job file."
                logger.critical(msg, extra=log_extra)
                raise PocketMapperError(msg)

        # 3. Paths left unset by both
        return resolve_paths(values)

    def resolve_settings(self, values):
        """
        Check and resolve every setting, and write job_settings.json.

        Empties `temp_dir`. Probes the foldseek binary when that aligner is
        selected. Logs, so run it with the run's log open.

        Args:
            values (dict): Settings field name -> value, from `configure_workflow`.

        Returns:
            Settings: The resolved configuration. Also written to `job_settings_path`.

        Raises:
            PocketMapperError: If a setting has an unknown value, or foldseek is selected but cannot run.
        """
        log_extra = {"stage": "Configuring Settings"}

        # 4a. Empty this run's scratch space, so a rerun into the same results_dir does not inherit
        # the last run's structures. After the log is open, so the warning for a temp_dir that cannot
        # be emptied reaches info.log.
        empty_temp_dir(values["temp_dir"], [values["cache_dir"], values["results_dir"]])

        # 4b. Validate the aligner and, for foldseek, probe the binary. Must come before anything is
        # fetched, so a missing binary fails without wasted downloads, and before the settings are
        # dumped below, so job_settings.json records the normalised value.
        values["aligner"] = resolve_aligner(values["aligner"])
        if values["aligner"] == "foldseek":
            require_foldseek(
                "The foldseek aligner was selected",
                "Or pass --aligner seq to use the built-in BLOSUM62 sequence aligner.",
            )

        # 4c. Reads the aligner just validated, so it must follow it. After the log is open, so the
        # 'auto' resolution is visible.
        values["align_struct_method"] = resolve_align_struct_method(values["align_struct_method"], values["aligner"])

        # 4d. Same reasoning as 4b/4c: after the log is open so the resolution is visible, and
        # before the settings are logged and dumped, so job_settings.json records a concrete count.
        values["threads"] = resolve_threads(values["threads"])

        # 4e. Before the settings are dumped, so job_settings.json records a checked value.
        values["delete_tmp"] = resolve_delete_tmp(values["delete_tmp"])

        # 4f. Before the settings are dumped, so job_settings.json records the normalised value.
        values["pisa_source"] = resolve_pisa_source(values["pisa_source"])

        settings = Settings(**values)
        logger.info(f"Settings: {json.dumps(asdict(settings), indent=4)}", extra=log_extra)

        # 5. Output dump
        try:
            os.makedirs(os.path.dirname(settings.job_settings_path), exist_ok=True)
            with open(settings.job_settings_path, "w") as f:
                json.dump(asdict(settings), f, indent=4)
            logger.info(f"Settings successfully dumped to {settings.job_settings_path}", extra=log_extra)
        except Exception as e:
            logger.error(f"Failed to dump settings to {settings.job_settings_path}: {e}", extra=log_extra)
        return settings
