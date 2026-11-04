"""
PocketMapper: map and compare binding pockets across protein structures.

`PocketMapper.search()` is the whole workflow. It resolves its settings, then runs the steps in order:

1. `lib.settings.layer_command` -> arguments over job file; `resolve_settings` -> every setting checked,
   job_settings.json.
2. `parse.parse` -> every entry checked before anything is fetched; entries.tsv written afresh.
3. `fetch.fetch` -> structures and a bundled Foldseek database into the cache.
4. `align.align` -> foldseek or the local sequence aligner, per `aligner` -> alignment.tsv.
5. `pockets.pockets` -> PISA interfaces into the cache, then pockets.json. Against a PDB Foldseek
   database, also a pisa pocket per interface of each hit. The Pocket shape itself is declared in
   `pockets/pocket.py`.
6. `compare.compare` -> pocket_comparison.tsv.
7. `superpose.superpose` -> the top align_count targets per query, superposed into aligned_structures/.

Each step skips what it cannot use and marks it, with the reason, in entries.tsv.

Steps: each step is a package, and `pocketmapper.<step>.<step>` holds its two layers. The entry function,
named after the step's command (`parse`, `fetch_structures`, `align`, `pockets`, `compare`, `superpose`),
layers its settings, logs to `log_path` for the length of the call, validates its settings, manages its
scratch directory and calls the core function, so each step also runs on its own. The core function
takes explicit paths, values and records as arguments, and keeps no state between calls.

Settings: an argument given beats the job file (a JSON path, or the same already loaded as a dict), which
beats `lib.options.SETTING_DEFAULTS`. Only `search` and `parse` take `cache_dir` and `results_dir` and
derive the other paths from them; each writes its settings, every path resolved, as a job file for the
steps after it (job_settings.json, `<results_dir>/parse_settings.json`). Every other step takes each path
it uses from its arguments or job file, and an argument moves only the path it names.

Entries: `parse` writes entries.tsv, a row per query and target entry (`lib.entries`). Every later step
reads its records from it, skipping the entries an earlier step dropped, and marks the ones it drops
itself. The cache and the files under `results_dir` are all the steps hand on. A step clears its own
marks when it starts, so every step is idempotent on rerun.

Settings travel as a job dict, job key -> value, `command_keys("search")` for search. Entry functions take
one as `job_file`; core functions and components take only the individual values they need, so none of
them needs a job dict to be usable on its own.

Author: Lachlan Ellingboe
"""

import json
import logging

from pocketmapper.align.align import align
from pocketmapper.compare.compare import compare
from pocketmapper.fetch.fetch import fetch_structures
from pocketmapper.lib.entries import entry_sides, fsdb_record, read_entries
from pocketmapper.lib.foldseek import require_foldseek
from pocketmapper.lib.logs import log_to_file
from pocketmapper.lib.options import command_keys
from pocketmapper.lib.scopes import run_scope, temp_dir_scope
from pocketmapper.lib.settings import (
    check_fsdb_align_struct_method,
    check_fsdb_aligner,
    check_offline_refresh,
    dump_settings,
    layer_command,
    require_setting,
    resolve_auto_align_struct_method,
    resolve_paths,
    validate_options,
)
from pocketmapper.parse.parse import parse
from pocketmapper.pockets.pockets import pockets
from pocketmapper.superpose.superpose import superpose

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
        local_root=None,
        verbosity=None,
        threads=None,
        aligner=None,
        align_count=None,
        align_struct_method=None,
        query_pocket_method=None,
        target_pocket_method=None,
        delete_tmp=None,
        pisa_source=None,
        pisa_bsa_threshold=None,
        offline=None,
        refresh_structures=None,
        pdb_dir=None,
        alphafold_dir=None,
        pocket_dir=None,
        foldseek_preprocessed_structure_dir=None,
        temp_dir=None,
        aligned_structure_dir=None,
        alignment_path=None,
        pocket_comparison_path=None,
        pockets_path=None,
        entries_path=None,
        job_settings_path=None,
        log_path=None,
        fsdb_dir=None,
    ):
        """
        Run the full PocketMapper search workflow.

        For the length of the call, sets the `pocketmapper` logger's level and adds a handler writing
        to `log_path`; both are restored on return. Empties `temp_dir` on the way in and, unless
        `delete_tmp` is 0, deletes it on a clean return; a raised error leaves it in place. Each
        directory is created only when something is first written into it. Besides the results,
        writes pockets.json, entries.tsv and job_settings.json into `results_dir`, so any
        step can be rerun there on its own, with job_settings.json as its job file.

        Args:
            query (str, optional): Query identifier, string or path to a list. Required here or in
                `job_file`, not both.
            target (str, optional): Target structure identifier, string or path to a list. As `query`.
            job_file (str or dict, optional): JSON job file of job key -> value, or the same
                already loaded. Any argument given overrides it. A key search does not take is rejected
                unless it is one of `lib.options.STEP_ONLY_KEYS`, which only some steps take, and is
                left out of job_settings.json.
            cache_dir (str, optional): Directory to cache intermediate structures.
                Defaults to DEFAULT_CACHE_DIR.
            results_dir (str, optional): Directory to output results to.
                Defaults to pocketmapper_results_<YYMMDD_HHMMSS>.
            local_root (str, optional): Directory that relative local structure files and Foldseek databases
                resolve against. Defaults to the working directory.
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
                'pisa', 'passthrough', 'vdw' or 'whole_chain' -- or 'auto' (the default) to infer it
                per entry from the input string.
            target_pocket_method (str, optional): As `query_pocket_method`, for the target side, which
                also takes 'foldseek_db'.
            delete_tmp (int, optional): 1 deletes temp_dir at the end of a successful run; 0 keeps
                it for inspection. Defaults to DEFAULT_DELETE_TMP.
            pisa_source (str, optional): Where PISA interfaces are fetched from -- 'ftp' for the EBI
                FTP server or 'api' for the paced PDBe API. Both serve the same data into the same
                cache. Defaults to DEFAULT_PISA_SOURCE.
            pisa_bsa_threshold (float, optional): Build pisa pockets from buried surface: a residue
                joins when the partner chain buries more than this percent of its accessible surface,
                at least 0 and below 100. Defaults to None, which builds them from interface bonds.
            offline (int, optional): 1 downloads nothing -- no structures, PISA interfaces or Foldseek
                databases -- and skips and marks whatever is not cached; 0 downloads what the cache
                lacks. Defaults to DEFAULT_OFFLINE.
            refresh_structures (int, optional): 1 re-downloads the PDB and AlphaFold structures in the
                fetch step even when they are cached; 0 keeps them. Cannot be 1 with `offline` 1.
                Defaults to DEFAULT_REFRESH_STRUCTURES.
            pdb_dir (str, optional): Cache of fetched PDB structures.
                Defaults to <cache_dir>/pdb_structures.
            alphafold_dir (str, optional): Cache of fetched AlphaFold structures.
                Defaults to <cache_dir>/alphafold_structures.
            pocket_dir (str, optional): Cache of PISA data. Defaults to <cache_dir>/pockets.
            foldseek_preprocessed_structure_dir (str, optional): Cache of the single-chain structures
                Foldseek is given. Defaults to <cache_dir>/foldseek_preprocessed_structures.
            temp_dir (str, optional): Per-run scratch, wiped before use and deleted at the end of a
                successful run. Defaults to <results_dir>/tmp.
            aligned_structure_dir (str, optional): Where the superposed structures for the top hits
                are written. Defaults to <results_dir>/aligned_structures.
            alignment_path (str, optional): Where the alignment table is written.
                Defaults to <results_dir>/alignment.tsv.
            pocket_comparison_path (str, optional): Where the pocket comparison table is written.
                Defaults to <results_dir>/pocket_comparison.tsv.
            pockets_path (str, optional): Where the pockets are written.
                Defaults to <results_dir>/pockets.json.
            entries_path (str, optional): Where the entries table is written, each entry dropped along
                the way marked with the reason. Defaults to <results_dir>/entries.tsv.
            job_settings_path (str, optional): Where this run's resolved settings are dumped. Never
                read from a job file. Defaults to <results_dir>/job_settings.json.
            log_path (str, optional): Where the run log is written.
                Defaults to <results_dir>/info.log.
            fsdb_dir (str, optional): Cache of bundled Foldseek databases.
                Defaults to <cache_dir>/fsdb.

        Returns:
            dict: The resolved settings, job key -> value, as written to job_settings.json plus
                `job_settings_path`. Results are written to `results_dir` -- read pocket_comparison.tsv
                and alignment.tsv from the paths it names.

        Raises:
            PocketMapperError: On a bad job file, a query/target given both ways or neither way, a
                setting with a bad value, offline and refresh_structures both 1, or any pipeline
                failure.
        """
        # First, while the locals are only the parameters
        values = layer_command("search", job_file, locals())
        for key in ("query", "target"):
            require_setting(values, key)
        values = resolve_paths(values, "search")

        # The steps' own scopes are no-ops inside this one, so only job_settings.json is written
        with run_scope("search"), log_to_file(values["log_path"], values["verbosity"]):
            # Every step gets the whole run's settings; their own log and temp scopes are no-ops
            # inside these, so the log is written once and temp_dir lives for the whole run
            job = self.resolve_settings(values)
            with temp_dir_scope(job["temp_dir"], job["delete_tmp"]):
                table = parse(job_file=job)
                logger.info(
                    f"Entries:\n{table.to_string(index=False)}", extra={"stage": "Determine Query/Target Types"}
                )
                # Checked before the database is downloaded
                if fsdb_record(entry_sides(read_entries(job["entries_path"]))["target"]) is not None:
                    check_fsdb_aligner(job["aligner"])
                    check_fsdb_align_struct_method(job["align_struct_method"])
                fetch_structures(job_file=job)
                align(job_file=job)
                pockets(job_file=job)
                compare(job_file=job)
                superpose(job_file=job)

            logger.info("PocketMapper search completed successfully.", extra={"stage": "End"})

        return job

    def resolve_settings(self, values):
        """
        Check and resolve every setting, and write job_settings.json.

        Probes the foldseek binary when that aligner is selected. Logs, so run it with the run's log
        open.

        Args:
            values (dict): Job key -> value, every path resolved. Keys search does not take, which
                only some steps take, are dropped.

        Returns:
            dict: Job key -> value for each of `command_keys("search")`, checked and normalised. Also
                written to `job_settings_path`, without that key, as a job file any step can take.

        Raises:
            PocketMapperError: If a setting has an unknown value, offline and refresh_structures are
                both 1, foldseek is selected but cannot run, or the settings cannot be written.
        """
        log_extra = {"stage": "Configuring Settings"}

        # Everything is checked before anything is fetched, so a bad value fails without wasted downloads,
        # and before the settings are dumped below, so job_settings.json records checked, normalised
        # values. After the log is open, so a failure and the `threads` and 'auto' resolutions reach it.
        values = validate_options(values, "search")
        check_offline_refresh(values["offline"], values["refresh_structures"])
        # Reads the aligner just validated. job_settings.json records a concrete method.
        values["align_struct_method"] = resolve_auto_align_struct_method(
            values["align_struct_method"], values["aligner"]
        )
        # Last: the one check that runs a program
        if values["aligner"] == "foldseek":
            require_foldseek(
                "The foldseek aligner was selected",
                "Or pass --aligner seq to use the built-in BLOSUM62 sequence aligner.",
            )

        job = {key: values[key] for key in command_keys("search")}
        logger.info(f"Settings: {json.dumps(job, indent=4)}", extra=log_extra)
        dump_settings(job)
        return job
