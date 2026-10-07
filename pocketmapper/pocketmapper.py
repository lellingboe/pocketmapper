"""
PocketMapper: map and compare binding pockets across protein structures.

`PocketMapper.search()` is the whole workflow. It resolves its settings, then runs the steps in
`pocketmapper.steps` in order. Each step parses the query and target entries again from the same
settings; what they hand each other is the cache and the files under `results_dir`:

1. `configure_workflow` -> arguments over job file; `resolve_settings` -> Settings, job_settings.json.
2. `steps.parse` -> every entry checked before anything is fetched; failed_entries.json started
   afresh.
3. `steps.fetch_structures` -> structures and a bundled Foldseek database into the cache.
4. `steps.align` -> foldseek or the local sequence aligner, per `aligner` -> alignment.tsv.
5. `steps.pockets` -> PISA interfaces into the cache, then pockets.json. Against a PDB Foldseek
   database, also a pisa pocket per interface of each hit. The Pocket shape itself is declared in
   `pockets/pocket.py`.
6. `steps.compare` -> pocket_comparison.tsv.
7. `steps.superpose` -> the top align_count targets per query, superposed into aligned_structures/.

Each step skips what it cannot use and lists it, with the reason, in failed_entries.json. `search`
calls each step's entry function with its resolved Settings as the job file, so each step also runs
on its own, the same way; see `pocketmapper.steps`.

This is the only module that builds a `Settings`. Steps take a job dict keyed by its field names,
and components are handed the individual values they need, so none of them has to build one to be
usable on its own.

Author: Lachlan Ellingboe
"""

import json
import logging
from dataclasses import asdict
from dataclasses import fields

from pocketmapper.entries import fsdb_record
from pocketmapper.lib import log_to_file
from pocketmapper.lib import run_scope
from pocketmapper.lib import temp_dir_scope
from pocketmapper.settings import Settings
from pocketmapper.settings import check_fsdb_align_struct_method
from pocketmapper.settings import check_fsdb_aligner
from pocketmapper.settings import dump_settings
from pocketmapper.settings import layer_settings
from pocketmapper.settings import require_foldseek
from pocketmapper.settings import require_setting
from pocketmapper.settings import resolve_align_struct_method
from pocketmapper.settings import resolve_aligner
from pocketmapper.settings import resolve_delete_tmp
from pocketmapper.settings import resolve_paths
from pocketmapper.settings import resolve_pisa_source
from pocketmapper.settings import resolve_threads
from pocketmapper.steps.align import align
from pocketmapper.steps.compare import compare
from pocketmapper.steps.fetch_structures import fetch_structures
from pocketmapper.steps.parse import parse
from pocketmapper.steps.parse import parse_job_entries
from pocketmapper.steps.pockets import pockets
from pocketmapper.steps.superpose import superpose

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
        cache_dir=None,
        results_dir=None,
        work_dir=None,
        verbosity=None,
        threads=None,
        aligner=None,
        align_count=None,
        align_struct_method=None,
        query_pocket_method=None,
        target_pocket_method=None,
        delete_tmp=None,
        pisa_source=None,
        pdb_dir=None,
        alphafold_dir=None,
        pocket_dir=None,
        foldseek_preprocessed_structure_dir=None,
        temp_dir=None,
        aligned_structure_dir=None,
        alignment_path=None,
        pocket_comparison_path=None,
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
        is first written into it. Besides the results, writes pockets.json, failed_entries.json and
        job_settings.json into `results_dir`, so any step can be rerun there on its own, with
        job_settings.json as its job file.

        Args:
            query (str, optional): Query identifier, string or path to a list. Required here or in
                `job_file`, not both.
            target (str, optional): Target structure identifier, string or path to a list. As `query`.
            job_file (str or dict, optional): JSON job file of job key -> value, or the same
                already loaded. Any argument given overrides it; a key Settings lacks is ignored.
            cache_dir (str, optional): Directory to cache intermediate structures.
                Defaults to DEFAULT_CACHE_DIR.
            results_dir (str, optional): Directory to output results to.
                Defaults to pocketmapper_results_<YYMMDD_HHMMSS>.
            work_dir (str, optional): Directory that entries and relative paths resolve against.
                Defaults to the working directory.
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
            pockets_path (str, optional): Where the pockets are written.
                Defaults to <results_dir>/pockets.json.
            failed_entries_path (str, optional): Where the entries dropped along the way are written.
                Defaults to <results_dir>/failed_entries.json.
            job_settings_path (str, optional): Where this run's resolved settings are dumped. Never
                read from a job file. Defaults to <results_dir>/job_settings.json.
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
            "work_dir": work_dir,
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
            "pockets_path": pockets_path,
            "failed_entries_path": failed_entries_path,
            "job_settings_path": job_settings_path,
            "log_path": log_path,
            "fsdb_dir": fsdb_dir,
        }

        values = self.configure_workflow(job_file, arguments)

        # The steps' own scopes are no-ops inside this one, so only job_settings.json is written
        with run_scope("search"), log_to_file(values["log_path"], values["verbosity"]):
            settings = self.resolve_settings(values)
            # Every step gets the whole run's settings; their own log and temp scopes are no-ops
            # inside these, so the log is written once and temp_dir lives for the whole run
            job = asdict(settings)
            with temp_dir_scope(settings.temp_dir, settings.delete_tmp, [settings.cache_dir, settings.results_dir]):
                parse(job_file=job)
                # Checked before the database is downloaded
                if fsdb_record(parse_job_entries(job, "search")["target"]) is not None:
                    check_fsdb_aligner(settings.aligner)
                    check_fsdb_align_struct_method(settings.align_struct_method)
                fetch_structures(job_file=job)
                align(job_file=job)
                pockets(job_file=job)
                compare(job_file=job)
                superpose(job_file=job)

            logger.info("PocketMapper search completed successfully.", extra={"stage": "End"})

        return job

    def configure_workflow(self, job_file, arguments):
        """
        Layer the arguments over the job file and fill in every path left unset.

        Args:
            job_file (str, dict or None): JSON job file, the same already loaded, or None for none.
                Any argument given overrides it.
            arguments (dict): Settings field name -> value from `search()`, None for unset.

        Returns:
            dict: Job key -> value, every path set but nothing else checked.

        Raises:
            PocketMapperError: If the job file is missing, unreadable or names an unknown setting, or if
                query or target is given both ways or neither way.
        """
        values = layer_settings(job_file, arguments)
        for key in ("query", "target"):
            require_setting(values, key)
        return resolve_paths(values, "search")

    def resolve_settings(self, values):
        """
        Check and resolve every setting, and write job_settings.json.

        Probes the foldseek binary when that aligner is selected. Logs, so run it with the run's log
        open.

        Args:
            values (dict): Job key -> value, from `configure_workflow`. Keys Settings lacks, which
                only some steps take, are dropped.

        Returns:
            Settings: The resolved configuration. Also written to `job_settings_path`, without that
                field, as a job file any step can take.

        Raises:
            PocketMapperError: If a setting has an unknown value, foldseek is selected but cannot run, or
                the settings cannot be written.
        """
        log_extra = {"stage": "Configuring Settings"}

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

        settings = Settings(**{field.name: values[field.name] for field in fields(Settings)})
        logger.info(f"Settings: {json.dumps(asdict(settings), indent=4)}", extra=log_extra)

        # 5. Output dump
        dump_settings(asdict(settings))
        return settings
