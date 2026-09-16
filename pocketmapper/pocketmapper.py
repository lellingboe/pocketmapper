"""
PocketMapper: map and compare binding pockets across protein structures.

`search()` is the only public method; everything else on the class is an internal step of it. The
command line lives in `cli.py`, which is the only module that knows about argv or exit codes.

`search()` is the whole workflow, top to bottom:

1. `configure_workflow` -> job file over arguments -> Settings, directories, job_settings.json, logging.
2. `configure_query_target` -> QTProcessor -> one DataFrame of QTRecords per side.
3. `fetch_missing_structures` (or `fetch_missing_fsdb`) -> mmCIF into pdb_dir / alphafold_dir.
4. `alignment` -> foldseek or the local sequence aligner, per `aligner` -> alignment.tsv.
5. `get_pockets` -> a `retrieve_*_pockets` builder per pocket method, merged into one
   pocket_id -> Pocket dict. The Pocket shape itself is declared in `pocket.py`.
6. `compare_pockets_based_on_alignment` -> pocket_comparison.compare_pockets -> pocket_comparison.tsv.
7. `align_structs` -> superposes the top align_count targets per query into aligned_structures/.

This is the only module that knows about `Settings`; components are handed the individual values
they need, so none of them has to build one to be usable on its own.

Author: Lachlan Ellingboe
"""

import json
import logging
import logging.config
import os
import shutil
from dataclasses import asdict
from dataclasses import dataclass
from dataclasses import fields
from datetime import datetime

import pandas as pd

from pocketmapper.constants import ALIGN_STRUCT_METHODS
from pocketmapper.constants import ALIGNERS
from pocketmapper.constants import DEFAULT_ALIGN_COUNT
from pocketmapper.constants import DEFAULT_ALIGN_STRUCT_METHOD
from pocketmapper.constants import DEFAULT_ALIGNER
from pocketmapper.constants import DEFAULT_CACHE_DIR
from pocketmapper.constants import DEFAULT_DELETE_TMP
from pocketmapper.constants import DEFAULT_VERBOSITY
from pocketmapper.constants import FOLDSEEK_FORMAT_OUTPUT
from pocketmapper.constants import FOLDSEEK_INSTALL_HINT
from pocketmapper.constants import LOG_FORMAT
from pocketmapper.downloads.pisa_downloader import PisaDownloader
from pocketmapper.downloads.structure_downloader import StructureDownloader
from pocketmapper.exceptions import PocketMapperError
from pocketmapper.foldseek import bundled_foldseek_dbs
from pocketmapper.foldseek import check_foldseek
from pocketmapper.foldseek import run_foldseek
from pocketmapper.lib import StageFilter
from pocketmapper.lib import is_within
from pocketmapper.lib import jsonify_dict
from pocketmapper.lib import parse_foldseek_pdb_entry_name
from pocketmapper.lib import split_chain_info
from pocketmapper.pisa_parser import PisaParser
from pocketmapper.pocket_calculator import PocketCalculator
from pocketmapper.pocket_comparison import compare_pockets
from pocketmapper.pocket_parser import parse_pocket_from_struct
from pocketmapper.qt_processor import QTProcessor
from pocketmapper.sequence_aligner import SequenceAligner
from pocketmapper.structure_aligner import StructureAligner
from pocketmapper.structure_preprocessor import StructurePreprocessor


@dataclass
class Settings:
    """
    Fully resolved PocketMapper run configuration.

    A record of what a run actually used, built once at the end of `configure_workflow` and dumped
    to job_settings.json. It has no defaults: each field arrives from the job file, the arguments to
    `search()` or run-time resolution, in that priority order. The pocket methods are the only
    optional fields, because None means "infer the method from each entry".
    """

    query: str
    target: str
    cache_dir: str
    results_dir: str
    query_pocket_method: str | None
    target_pocket_method: str | None
    # "foldseek" or "seq" (the local BLOSUM62 sequence aligner).
    aligner: str
    align_count: int
    # "pocket" or "foldseek"; "auto" is resolved to one of those before a Settings is built.
    align_struct_method: str
    verbosity: int
    threads: int
    # temp_dir holds the per-run inputs actually handed to the aligner, and delete_tmp removes it
    # on the way out. False keeps it: a run that produced no rows is diagnosed from what it was
    # given, which is gone by the time anyone looks.
    delete_tmp: bool
    pdb_dir: str
    alphafold_dir: str
    pocket_dir: str
    foldseek_preprocessed_structure_dir: str
    temp_dir: str
    aligned_structure_dir: str
    alignment_path: str
    pocket_comparison_path: str
    job_settings_path: str
    log_path: str
    fsdb_dir: str


def resolve_paths(values):
    """
    Fill in `results_dir` and any derived path left unset.

    Paths already set -- via the job file or an argument -- are left untouched. `results_dir`
    defaults to a timestamped name; the derived paths default to locations under `cache_dir` or
    `results_dir`.

    Args:
        values (dict): Settings field name -> value, with `cache_dir` set.

    Returns:
        dict: A copy of `values` with every path set.
    """
    values = dict(values)
    if values["results_dir"] is None:
        values["results_dir"] = f"pocketmapper_results_{datetime.now().strftime('%y%m%d_%H%M%S')}"
    cache_dir = values["cache_dir"]
    results_dir = values["results_dir"]
    derived = {
        "pdb_dir": os.path.join(cache_dir, "pdb_structures"),
        "alphafold_dir": os.path.join(cache_dir, "alphafold_structures"),
        "pocket_dir": os.path.join(cache_dir, "pockets"),
        "foldseek_preprocessed_structure_dir": os.path.join(cache_dir, "foldseek_preprocessed_structures"),
        "temp_dir": os.path.join(results_dir, "tmp"),
        "aligned_structure_dir": os.path.join(results_dir, "aligned_structures"),
        "alignment_path": os.path.join(results_dir, "alignment.tsv"),
        "pocket_comparison_path": os.path.join(results_dir, "pocket_comparison.tsv"),
        "job_settings_path": os.path.join(results_dir, "job_settings.json"),
        "log_path": os.path.join(results_dir, "info.log"),
        "fsdb_dir": os.path.join(cache_dir, "fsdb"),
    }
    for key, path in derived.items():
        if values[key] is None:
            values[key] = path
    return values


class PocketMapper:
    """
    The pipeline.

    `search()` is the only public method and runs the whole workflow; see the module docstring for
    its steps. Every other method is an internal step and carries a leading underscore.

    Usable as a library -- nothing here needs a terminal -- but note `search()` has global side
    effects: it reconfigures the *root* logger via `logging.config.dictConfig`, and deletes its
    temporary directories on the way out.
    """

    def __init__(self):
        """
        Install a CRITICAL-only root handler and initialise the Foldseek-database flags.

        `configure_workflow` can fail on a bad settings file before `configure_logging` has run, so
        the root logger needs a handler that formats with `LOG_FORMAT` from construction. The
        handler is built by hand rather than through `basicConfig(format=...)` so `StageFilter` can
        be attached to it; without the filter a record carrying no `stage` fails to format.
        """
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter(LOG_FORMAT))
        handler.addFilter(StageFilter())
        logging.basicConfig(level=logging.CRITICAL, handlers=[handler])

        self.fsdb_target = False
        self.fsdb_pdb_target = False

    def configure_logging(self, verbosity, log_path):
        """
        Configure logging level and handlers based on user settings.

        Sets both a console stream handler and a file handler (to `log_path` in `settings`),
        adjusting the verbosity depending on the user's input.

        Args:
            verbosity (int): The verbosity level (4=DEBUG, 3=INFO, 2=WARNING, else ERROR).
            log_path (str): The path to the log file.
        """
        # Set log level based on verbosity setting (default to INFO if not set)
        log_level = None
        if verbosity == 4:
            log_level = "DEBUG"
        elif verbosity == 3:
            log_level = "INFO"
        elif verbosity == 2:
            log_level = "WARNING"
        else:
            log_level = "ERROR"

        log_config = {
            "version": 1,
            "formatters": {
                "standard": {"format": LOG_FORMAT},
            },
            # On the handlers rather than on the root logger: a handler filter also sees records
            # propagating up from third-party loggers, which never pass a "stage".
            "filters": {
                "stage": {"()": "pocketmapper.lib.StageFilter"},
            },
            "handlers": {
                "console": {
                    "class": "logging.StreamHandler",
                    "level": log_level,
                    "formatter": "standard",
                    "filters": ["stage"],
                    "stream": "ext://sys.stdout",
                },
                "file": {
                    "class": "logging.FileHandler",
                    "level": log_level,
                    "formatter": "standard",
                    "filters": ["stage"],
                    "filename": log_path,
                },
            },
            "root": {
                "handlers": ["console", "file"],
                "level": log_level,
                "propagate": True,
            },
        }
        logging.config.dictConfig(log_config)

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
        query_pocket_method=None,
        target_pocket_method=None,
        delete_tmp=DEFAULT_DELETE_TMP,
        pdb_dir=None,
        alphafold_dir=None,
        pocket_dir=None,
        foldseek_preprocessed_structure_dir=None,
        temp_dir=None,
        aligned_structure_dir=None,
        alignment_path=None,
        pocket_comparison_path=None,
        job_settings_path=None,
        log_path=None,
        fsdb_dir=None,
    ):
        """
        Orchestrate and run the full PocketMapper search workflow.

        Every value set in `job_file` wins over the matching argument here. `query` and `target`
        are required from exactly one of the two.

        Args:
            query (str, optional): Query identifier, string or path to a list.
            target (str, optional): Target structure identifier, string or path to a list.
            job_file (str, optional): Path to a JSON job file of Settings field name -> value.
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
                'pisa', 'passthrough', 'vdw', 'whole_chain' or 'foldseek_db'. Left unset, it is
                inferred per entry from the input string.
            target_pocket_method (str, optional): As `query_pocket_method`, for the target side.
            delete_tmp (bool, optional): Delete temp_dir at the end of the run. Defaults to
                DEFAULT_DELETE_TMP; False keeps it for inspection.
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
        # The Settings fields as passed to this call; the job file is layered on top. `job_file` is
        # deliberately absent: Settings has no such field. Keeping the dict here, next to the
        # signature it mirrors, is what keeps a new option from being added to one and not the other.
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
            "pdb_dir": pdb_dir,
            "alphafold_dir": alphafold_dir,
            "pocket_dir": pocket_dir,
            "foldseek_preprocessed_structure_dir": foldseek_preprocessed_structure_dir,
            "temp_dir": temp_dir,
            "aligned_structure_dir": aligned_structure_dir,
            "alignment_path": alignment_path,
            "pocket_comparison_path": pocket_comparison_path,
            "job_settings_path": job_settings_path,
            "log_path": log_path,
            "fsdb_dir": fsdb_dir,
        }

        self.settings = self.configure_workflow(job_file, arguments)
        self.query_df, self.target_df = (
            self.configure_query_target()
        )  # parses the query and target inputs to determine their types and sets up the relevant data structures for each entry

        self.query_df = self.fetch_missing_structures("query", self.query_df)
        if self.fsdb_target:
            self.fetch_missing_fsdb(self.target_df, self.foldseek_tmp_dir)  # Fetch any missing foldseek databases
        else:
            self.target_df = self.fetch_missing_structures("target", self.target_df)  # Fetch any missing structures

        self.alignment()  # Align the query and target structures using either local sequence alignment or foldseek based on the settings
        pockets = self.get_pockets()  # Adds seq_pos and ca-coords to the pocket info dict
        self.compare_pockets_based_on_alignment(pockets)
        self.align_structs()
        self.delete_tmp()

        logging.info("PocketMapper search completed successfully.", extra={"stage": "End"})

        return asdict(self.settings)

    def configure_workflow(self, job_file, arguments):
        """
        Build the fully resolved `Settings` for this run.

        Layers an optional JSON job file over the arguments passed to `search()` -- any value the job
        file sets wins -- then resolves the unset paths and run-dependent values, creates the
        directories and writes job_settings.json.

        Args:
            job_file (str or None): Path to a JSON job file, or None for none.
            arguments (dict): Settings field name -> value from `search()`.

        Returns:
            Settings: The resolved configuration. Also written to `job_settings_path`.

        Raises:
            PocketMapperError: If the job file is missing, unreadable or names an unknown setting, if
                query or target is given both ways or neither way, or if a directory cannot be made.
        """
        log_extra = {"stage": "Configuring Settings"}

        # 1. The job file, if any
        job = {}
        if job_file is not None:
            if not os.path.isfile(job_file):
                logging.critical(f"Job file not found: {job_file}", extra=log_extra)
                raise PocketMapperError(f"Job file not found: {job_file}")
            try:
                with open(job_file) as f:
                    job = json.load(f)
            except Exception as e:
                logging.critical(f"Error reading job file: {job_file}. Is it in JSON format?", extra=log_extra)
                raise PocketMapperError(f"Error reading job file: {job_file}. Is it in JSON format?") from e
            if not isinstance(job, dict):
                msg = f'Job file {job_file} must hold a JSON object of {{"option": value}}.'
                logging.critical(msg, extra=log_extra)
                raise PocketMapperError(msg)
            unknown = sorted(set(job) - {f.name for f in fields(Settings)})
            if unknown:
                msg = f"Unknown setting(s) in {job_file}: {', '.join(unknown)}"
                logging.critical(msg, extra=log_extra)
                raise PocketMapperError(msg)

        # 2. The job file wins over the arguments, except that query and target must come from exactly
        # one of the two: a silently discarded positional would search something other than what the
        # command line shows.
        for key in ("query", "target"):
            if job.get(key) is not None and arguments[key] is not None:
                msg = f"{key} is set both in the job file and as an argument; give it only once."
                logging.critical(msg, extra=log_extra)
                raise PocketMapperError(msg)
        values = {**arguments, **job}
        for key in ("query", "target"):
            if values[key] is None:
                msg = f"No {key} given; pass it as an argument or set it in the job file."
                logging.critical(msg, extra=log_extra)
                raise PocketMapperError(msg)

        # 3. Paths left unset by both
        values = resolve_paths(values)

        # Ensure all necessary directories exist before proceeding, creating them if needed.
        # results_dir is listed in its own right: configure_logging opens a file handler under it
        # immediately below, and every other path here is settable away from it.
        dirs_to_create = [
            "results_dir",
            "pdb_dir",
            "alphafold_dir",
            "pocket_dir",
            "foldseek_preprocessed_structure_dir",
            "aligned_structure_dir",
        ]
        for dir_key in dirs_to_create:
            path = values[dir_key]
            try:
                os.makedirs(path, exist_ok=True)
            except OSError as e:
                logging.critical(f"Error creating directory {path}", extra=log_extra)
                raise PocketMapperError(f"Error creating directory {path}") from e

        self.configure_logging(values["verbosity"], values["log_path"])

        # 4a. Lay out this run's scratch space. Must come after configure_logging: the root logger is
        # still at CRITICAL before it, so the warning for a temp_dir that cannot be wiped would be
        # swallowed. Nothing reads the scratch directories until step 4 of the pipeline.
        self.configure_temp_dir(values["temp_dir"], values["cache_dir"], values["results_dir"])

        # 4b. Validate the aligner and, for foldseek, probe the binary. Must come before anything is
        # fetched, so a missing binary fails without wasted downloads, and before the settings are
        # dumped below, so job_settings.json records the normalised value.
        values["aligner"] = self.resolve_aligner(values["aligner"])

        # 4c. Reads the aligner resolve_aligner just validated, so it must follow it. After
        # configure_logging, so the 'auto' resolution is visible.
        values["align_struct_method"] = self.resolve_align_struct_method(
            values["align_struct_method"], values["aligner"]
        )

        # 4d. Same reasoning as 4b/4c: after configure_logging so the resolution is visible, and
        # before the settings are logged and dumped, so job_settings.json records a concrete count.
        values["threads"] = self.resolve_threads(values["threads"])

        settings = Settings(**values)
        logging.info(f"Settings: {json.dumps(asdict(settings), indent=4)}", extra=log_extra)

        # 5. Output dump
        try:
            os.makedirs(os.path.dirname(settings.job_settings_path), exist_ok=True)
            with open(settings.job_settings_path, "w") as f:
                json.dump(asdict(settings), f, indent=4)
            logging.info(f"Settings successfully dumped to {settings.job_settings_path}", extra=log_extra)
        except Exception as e:
            logging.error(f"Failed to dump settings to {settings.job_settings_path}: {e}", extra=log_extra)
        return settings

    def configure_temp_dir(self, temp_dir, cache_dir, results_dir):
        """
        Empty this run's scratch directory and lay out the subdirectories under it.

        The wipe is guarded by `is_within` and the creation is not: `temp_dir` is settable by both the
        job file and the command line, so a mistyped one costs a stray directory rather than its
        contents, but the run still needs somewhere to put its scratch. Foldseek's own scratch
        subdirectory is not created here -- Foldseek creates it itself.

        Args:
            temp_dir (str): This run's scratch directory.
            cache_dir (str): The cache root; `temp_dir` is only emptied when under it or `results_dir`.
            results_dir (str): The results root.

        Returns:
            None: Sets `query_tmp_dir`, `target_tmp_dir` and `foldseek_tmp_dir` on the instance.

        Raises:
            PocketMapperError: If a scratch directory cannot be created.
        """
        log_extra = {"stage": "Configuring Workflow"}

        # Subdirectories rather than Settings fields: every Settings field is reachable from the
        # command line, and these are placed by --temp_dir alone.
        self.query_tmp_dir = os.path.join(temp_dir, "query_structures")
        self.target_tmp_dir = os.path.join(temp_dir, "target_structures")
        self.foldseek_tmp_dir = os.path.join(temp_dir, "foldseek_tmp")

        # A rerun into the same results_dir would otherwise hand Foldseek's createdb whatever the
        # previous run left behind.
        if is_within(temp_dir, [cache_dir, results_dir]):
            shutil.rmtree(temp_dir, ignore_errors=True)
        else:
            logging.warning(
                f"Reusing temp_dir {temp_dir} without emptying it: it is outside cache_dir "
                "and results_dir. Empty it yourself if a previous run left anything there.",
                extra=log_extra,
            )

        for path in (temp_dir, self.query_tmp_dir, self.target_tmp_dir):
            try:
                os.makedirs(path, exist_ok=True)
            except OSError as e:
                logging.critical(f"Error creating directory {path}", extra=log_extra)
                raise PocketMapperError(f"Error creating directory {path}") from e

    def resolve_aligner(self, aligner):
        """
        Validate the `aligner` setting and check that its binary can run.

        "seq" needs nothing external and never probes for foldseek. "foldseek" is an optional external
        binary, so it is checked here, before any structure is fetched, rather than failing partway
        through the run at the first `run_foldseek` call.

        Args:
            aligner (str): The requested aligner, "foldseek" or "seq", in any case.

        Returns:
            str: The aligner, lowercased.

        Raises:
            PocketMapperError: If the value is not one of ALIGNERS, or foldseek cannot be run.
        """
        log_extra = {"stage": "Configuring Settings"}

        # The CLI always hands over a str, but a job file can hold anything at all.
        normalised = aligner.lower() if isinstance(aligner, str) else aligner
        if normalised not in ALIGNERS:
            msg = f"Unknown aligner {aligner!r}. Choose one of: {', '.join(ALIGNERS)}."
            logging.critical(msg, extra=log_extra)
            raise PocketMapperError(msg)

        if normalised == "foldseek" and not check_foldseek():
            msg = (
                "The foldseek aligner was selected but 'foldseek' could not be run; it is either not on "
                f"PATH or not executable (run with --verbosity 4 for the reason). {FOLDSEEK_INSTALL_HINT} "
                "Or pass --aligner seq to use the built-in BLOSUM62 sequence aligner."
            )
            logging.critical(msg, extra=log_extra)
            raise PocketMapperError(msg)

        return normalised

    def resolve_align_struct_method(self, align_struct_method, aligner):
        """
        Turn the tri-value `align_struct_method` setting into "pocket" or "foldseek".

        "foldseek" uses Foldseek's whole-chain transform from alignment.tsv; "pocket" uses the
        superposition of the two pockets on their overlapping residues, which `compare_pockets`
        already writes to pocket_comparison.tsv. "auto" picks whichever the run can actually do:
        the local BLOSUM62 aligner writes "-" for the chain transform, so it has only the pocket one.

        An explicit "foldseek" with the "seq" aligner is an error rather than a silent switch to
        "pocket": better to fail before anything is downloaded than to hand back a method the user did
        not ask for.

        Args:
            align_struct_method (str): The requested method.
            aligner (str): This run's aligner, already validated.

        Returns:
            str: "pocket" or "foldseek".

        Raises:
            PocketMapperError: If the value is not one of ALIGN_STRUCT_METHODS, or "foldseek" was
                asked for with the "seq" aligner.
        """
        log_extra = {"stage": "Configuring Settings"}

        # The CLI always hands over a str, but a job file can hold anything at all.
        method = align_struct_method.lower() if isinstance(align_struct_method, str) else align_struct_method
        if method not in ALIGN_STRUCT_METHODS:
            msg = (
                f"Unknown align_struct_method {align_struct_method!r}. "
                f"Choose one of: {', '.join(ALIGN_STRUCT_METHODS)}."
            )
            logging.critical(msg, extra=log_extra)
            raise PocketMapperError(msg)

        if method == "auto":
            method = "foldseek" if aligner == "foldseek" else "pocket"
            logging.info(
                f"align_struct_method 'auto' resolved to '{method}' (aligner '{aligner}' is in use)",
                extra=log_extra,
            )
        elif method == "foldseek" and aligner != "foldseek":
            msg = (
                "align_struct_method 'foldseek' needs Foldseek's whole-chain transform, but this run "
                "uses the local BLOSUM62 aligner, which does not produce one. Use "
                "--align_struct_method pocket, or --aligner foldseek."
            )
            logging.critical(msg, extra=log_extra)
            raise PocketMapperError(msg)

        return method

    def resolve_threads(self, threads):
        """
        Turn an unset `threads` setting into a concrete core count.

        Unset means one per available core, which is also what Foldseek does when given no
        `--threads`, so the default changes nothing about how Foldseek runs. Resolving it here
        rather than leaving it None means job_settings.json records the number the run used.

        Args:
            threads (int or None): The requested core count.

        Returns:
            int: A positive core count.

        Raises:
            PocketMapperError: If `threads` is not a positive integer.
        """
        log_extra = {"stage": "Configuring Settings"}

        if threads is None:
            # os.cpu_count(), not os.process_cpu_count() (3.13+) or os.sched_getaffinity (Linux
            # only): the floor is 3.10 and this has to work on macOS. None on an exotic platform.
            threads = os.cpu_count() or 1
            logging.info(f"threads unset, using one per available core ({threads})", extra=log_extra)
            return threads

        # argparse guarantees an int, but a job file can hold anything at all -- and a bool is
        # an int to isinstance, while --threads is not a flag.
        if isinstance(threads, bool) or not isinstance(threads, int) or threads < 1:
            msg = f"threads must be a positive integer, got {threads!r}."
            logging.critical(msg, extra=log_extra)
            raise PocketMapperError(msg)

        return threads

    def configure_query_target(self):
        """
        Parse the query and target inputs into record DataFrames.

        Uses `QTProcessor` to parse, validate and resolve each side's structure type and pocket method.
        Also rejects `align_struct_method="pocket"` against a Foldseek-DB target here, before anything is
        fetched -- see the check itself for why one early rejection covers both kinds of database.

        Returns:
            tuple: (query_df, target_df). Note these are returned, not stored on the instance; the caller
                assigns them.

        Raises:
            PocketMapperError: If either side is unusable, or on the align_struct_method rejection above.
        """
        log_extra = {"stage": "Determine Query/Target Types"}

        qtprocessor = QTProcessor(
            pdb_dir=self.settings.pdb_dir,
            alphafold_dir=self.settings.alphafold_dir,
            foldseek_preprocessed_structure_dir=self.settings.foldseek_preprocessed_structure_dir,
            fsdb_dir=self.settings.fsdb_dir,
        )
        q_df = qtprocessor.process_qt_cmdline_input(
            qt_input=self.settings.query,
            name="query",
            pocket_method=self.settings.query_pocket_method,
        )
        t_df = qtprocessor.process_qt_cmdline_input(
            qt_input=self.settings.target,
            name="target",
            pocket_method=self.settings.target_pocket_method,
        )

        errors = []
        if len(q_df) < 1:
            logging.critical("No valid query entries after processing", extra=log_extra)
            errors.append("no valid query entries")
        if len(t_df) < 1:
            logging.critical("No valid target entries after processing", extra=log_extra)
            errors.append("no valid target entries")
        if errors:
            raise PocketMapperError("; ".join(errors))

        if t_df.loc[0, "struct_type"] == "foldseek_db":
            if self.settings.aligner == "foldseek":
                self.fsdb_target = True
                # Neither kind of Foldseek DB can be superposed on its pocket, and which kind this is
                # is not known until expand_fsdb_pdb_targets has read the hit names -- so reject both
                # here, before anything is fetched. Unreachable from "auto", which resolves to
                # "foldseek" whenever the foldseek aligner is in use.
                if self.settings.align_struct_method == "pocket":
                    msg = (
                        "align_struct_method 'pocket' is not available against a Foldseek database "
                        "target. A human_domains-style hit has no coordinates to superpose at all, and "
                        "a PDB database's structures are assemblies while its pockets come from the "
                        "wwPDB asymmetric unit, so a pocket fit would be applied in the wrong frame. "
                        "Use --align_struct_method foldseek."
                    )
                    logging.critical(msg, extra=log_extra)
                    raise PocketMapperError(msg)
            else:
                msg = "A Foldseek database target requires --aligner foldseek."
                logging.critical(msg, extra=log_extra)
                raise PocketMapperError(msg)
        return q_df, t_df

    def fetch_missing_structures(self, name, qt_df):
        """
        Download the reference structures a side's records need.

        Uses `StructureDownloader`, which fetches each record to its own `struct_path`; records whose
        fetch fails are flagged rather than dropped, so the reason survives into the results. A side
        with nothing left is an error.

        Args:
            name (str): Which side this is, e.g. "query" or "target". Used in logging.
            qt_df (pandas.DataFrame): That side's records.

        Returns:
            pandas.DataFrame: The records with `success` and `failure_reason` updated.

        Raises:
            PocketMapperError: If no structure for this side could be fetched.
        """
        log_extra = {"stage": "Downloading Structures"}
        logging.debug(f"{name.capitalize()} data before fetching structures: \n{qt_df.head()}", extra=log_extra)
        structure_downloader = StructureDownloader()
        unique_records = qt_df.drop_duplicates(subset="struct_info").to_dict(orient="records")
        results = structure_downloader.download_missing_structures(unique_records)
        logging.debug(f"Structure fetcher results: {results}", extra=log_extra)

        # Update the dataframe with success/failure information
        qt_df["success"] = qt_df["struct_info"].map(results).fillna(False)
        qt_df.loc[~qt_df["success"], "failure_reason"] = "structure_not_found"

        # Logging results of structure fetching and updating query and target data with success/failure info
        logging.info(
            f"{sum(results.values())}/{len(results)} {name} required structures available",
            extra=log_extra,
        )
        if len(qt_df.query("success == False")) > 0:
            logging.warning(
                f"Missing structures for {name}(s): {', '.join(qt_df.loc[~qt_df['success'], 'pocket_id'].unique().tolist())}",
                extra=log_extra,
            )

        # Verifying sufficient structures were found to continue
        if qt_df["success"].sum() < 1:
            logging.critical(f"Insufficient {name} structures after fetching", extra=log_extra)
            raise PocketMapperError(f"Insufficient {name} structures after fetching. No valid {name} entries remain.")
        return qt_df

    def fetch_missing_fsdb(self, qt_df, tmp_dir):
        """
        Download a bundled Foldseek database if it is not already on disk.

        Args:
            qt_df (pandas.DataFrame): The target records; the database is named by row 0.
            tmp_dir (str): Scratch directory for `foldseek databases`.

        Returns:
            None: The database is written to the record's `struct_path`.

        Raises:
            PocketMapperError: If the destination cannot be created, or -- from `run_foldseek` --
                if the download fails.
        """
        log_extra = {"stage": "Fetching Missing Foldseek Database"}
        fsdb_name = qt_df.loc[0, "struct_info"].upper()
        fsdb_path = qt_df.loc[0, "struct_path"]
        if not os.path.exists(fsdb_path):
            logging.info(f"Fetching bundled Foldseek database '{fsdb_name}' to {fsdb_path}", extra=log_extra)
            try:
                os.makedirs(os.path.dirname(fsdb_path), exist_ok=True)
            except Exception as e:
                msg = f"Failed to create the directory for Foldseek database '{fsdb_name}' at {fsdb_path}: {e}"
                logging.critical(msg, extra=log_extra)
                raise PocketMapperError(msg) from e
            run_foldseek(
                ["databases", fsdb_name, fsdb_path, tmp_dir, "--threads", str(self.settings.threads)],
                log_extra,
            )
            logging.info(f"Successfully fetched Foldseek database '{fsdb_name}'", extra=log_extra)

    def alignment(self):
        """
        Coordinate structural alignment routes bridging query and target proteins.

        Dispatches to `foldseek_alignment()` for the "foldseek" aligner, else to
        `local_alignment()`. Requisites like `foldseek_preprocessing()`
        precede foldseek routines.

        Returns:
            None
        """
        log_extra = {"stage": "Alignment"}
        if self.settings.aligner == "foldseek":
            logging.info("Preprocessing structures for Foldseek...", extra=log_extra)
            self.foldseek_preprocessing()
            logging.info("Running Foldseek easy-search...", extra=log_extra)
            self.foldseek_alignment()
        else:
            logging.info("Running local pairwise aligner...", extra=log_extra)
            self.local_alignment()

    def foldseek_preprocessing(self):
        """
        Split complex PDB structures into single chain mmCIF files for Foldseek processing.
        Copies structure to the query/target directories, ensuring that the Foldseek preprocessed structure directory is populated.

        Instantiates a `StructurePreprocessor` mapping items sourced from the
        local datastores towards isolated files targeting specific sequences.

        Returns:
            None
        """
        log_extra = {"stage": "Preprocessing Structures"}

        structure_preprocessor = StructurePreprocessor()
        qtdf_dir_iter = [(self.query_df, self.query_tmp_dir)]
        if not self.fsdb_target:
            qtdf_dir_iter.append((self.target_df, self.target_tmp_dir))

        for df, search_dir in qtdf_dir_iter:
            records = df.drop_duplicates(subset=["preprocess_name", "chain_info"]).to_dict(orient="records")
            logging.debug(f"Records to preprocess: {json.dumps(records, indent=4)}", extra=log_extra)

            structure_preprocessor.set_output_directory(self.settings.foldseek_preprocessed_structure_dir)
            structure_preprocessor.update_cache()
            results = structure_preprocessor.preprocess_records(records=records, search_dir=search_dir)
            logging.debug(f"Preprocessing results: {json.dumps(results, indent=4)}", extra=log_extra)

            # Updating success and failure cols based on preprocessing results
            df.set_index("pocket_id", inplace=True)
            for index, success in results.items():
                if not success:
                    df.loc[index, "success"] = False
                    df.loc[index, "failure_reason"] = "structure_preprocessing_failed"
            df.reset_index(inplace=True)

        logging.info("Finished preprocessing structures", extra=log_extra)
        logging.debug(f"Query data after preprocessing: \n{self.query_df.head()}", extra=log_extra)
        logging.debug(f"Target data after preprocessing: \n{self.target_df.head()}", extra=log_extra)

    def foldseek_alignment(self):
        """
        Build the query (and, unless the target is a database, target) Foldseek DB, then search them.

        `foldseek easy-search` writes its raw tabular matches to `self.settings.alignment_path`.
        Every invocation goes through `foldseek.run_foldseek`, which logs and reports the failures.

        Returns:
            None
        """
        log_extra = {"stage": "Foldseek Alignment"}
        logging.info("Running Foldseek alignment...", extra=log_extra)

        # Setting up paths for foldseek databases
        self.query_db_path = os.path.join(self.query_tmp_dir, "query_db")
        run_foldseek(
            ["createdb", self.query_tmp_dir, self.query_db_path, "--threads", str(self.settings.threads)],
            log_extra,
        )

        if self.fsdb_target:
            self.target_db_path = self.target_df.loc[0, "struct_path"]
            logging.debug(f"Targeting bundled human_domains Foldseek DB at {self.target_db_path}", extra=log_extra)
        else:
            self.target_db_path = os.path.join(self.target_tmp_dir, "target_db")
            run_foldseek(
                ["createdb", self.target_tmp_dir, self.target_db_path, "--threads", str(self.settings.threads)],
                log_extra,
            )

        query_target_align_cmd = [
            "easy-search",
            self.query_db_path,
            self.target_db_path,
            self.settings.alignment_path,
            self.foldseek_tmp_dir,
            "--format-output",
            FOLDSEEK_FORMAT_OUTPUT,
            "--format-mode",
            "4",
            "-e",
            "0.001",
            "--file-include",
            r".*\.cif\.gz",
            "--max-seqs",
            "5000",
            "--threads",
            str(self.settings.threads),
            "-v",  # verbosity
            str(
                min(3, self.settings.verbosity)
            ),  # cap foldseek verbosity at 3 (info level) since it can be very verbose at higher levels and we already have our own logging verbosity control
        ]
        run_foldseek(query_target_align_cmd, log_extra)
        logging.debug("Foldseek alignment completed successfully", extra=log_extra)

    def local_alignment(self):
        """
        Execute traditional sequence-level alignment across inputs.

        Uses the `SequenceAligner` to locally match structured frames, resolving pairwise dependencies
        within internal dataframe representations. Outputs directly to the configuration TSV alignment file.

        Returns:
            None
        """
        log_extra = {"stage": "Local Alignment"}
        logging.info("Running local sequence alignments...", extra=log_extra)

        # Run the sequence aligner on the same preprocessed structures as foldseek uses
        aligner = SequenceAligner()
        query_records = (
            self.query_df.query("success").drop_duplicates(subset="preprocess_name").to_dict(orient="records")
        )
        target_records = (
            self.target_df.query("success").drop_duplicates(subset="preprocess_name").to_dict(orient="records")
        )
        alignment = aligner.align_records(
            query_records,
            target_records,
        )
        alignment.to_csv(self.settings.alignment_path, index=False, sep="\t")

    def dump_pockets(self, pockets, filename):
        """
        Write a pocket collection to the pocket cache directory as JSON, for inspection only.

        Nothing reads these files back -- pockets are recomputed every run -- so they are debug
        artefacts rather than a cache. `asdict` gives the nested shape declared in `pocket.py`:
        residues sit under a `residues` key rather than alongside the metadata. Written compact for
        the same reason: on a bundled-`pdb` Foldseek run `pisa_pockets.json` holds thousands of
        pockets, and nothing reads it to justify the whitespace.

        A None survives into the file as `null` rather than raising here. Three of the four builders
        store whatever their producer returned, including None for a missing structure or chain, and
        the failure for those belongs where it already is -- in `compare_pockets`.

        Args:
            pockets (dict): pocket_id -> Pocket (or None).
            filename (str): Basename to write under `pocket_dir`.

        Returns:
            None: Writes a file.
        """
        serialisable = {pid: asdict(pocket) if pocket is not None else None for pid, pocket in pockets.items()}
        with open(os.path.join(self.settings.pocket_dir, filename), "w") as f:
            json.dump(serialisable, f)

    def download_pisa_interfaces(self, pdb_list):
        """
        Download PISA summaries, assemblies and interfaces for `pdb_list` into the pocket cache.

        Both callers -- `expand_fsdb_pdb_targets` and `retrieve_pisa_pockets` -- go through here, so
        the two share one set of directories. Let them compute their own and the second re-downloads
        everything the first already fetched, behind PisaDownloader's per-entry sleep.

        Args:
            pdb_list (list): PDB IDs to fetch interfaces for.

        Returns:
            str: The interface directory -- the only one of the three anything reads back.
        """
        pisa_response_dir = os.path.join(self.settings.pocket_dir, "pisa_responses")
        interface_dir = os.path.join(pisa_response_dir, "interfaces")
        PisaDownloader().download_missing_interfaces(
            pdb_list=pdb_list,
            summary_dir=os.path.join(pisa_response_dir, "summaries"),
            asm_dir=os.path.join(pisa_response_dir, "assemblies"),
            interface_dir=interface_dir,
            error_path=os.path.join(pisa_response_dir, "errors.json"),
        )
        return interface_dir

    def select_pocket_records(self, pocket_method, label, dedup_subset=None):
        """
        Select the successfully fetched records of one pocket method, from both sides at once.

        Args:
            pocket_method (str): The `pocket_method` value to match.
            label (str): Human-readable name of the method, used in the log lines.
            dedup_subset (list, optional): Columns to drop duplicate records on. Only `pisa` passes
                one: its pocket is fully determined by structure and chain, and the Foldseek-DB
                expansion generates the same pair once per hit. The other methods must not dedup on
                structure and chain -- two passthrough pockets on one chain differ only in
                `residue_info`. Defaults to None, meaning keep every record.

        Returns:
            pandas.DataFrame: The matching records, possibly empty.
        """
        log_extra = {"stage": f"Retrieving {label} Pockets"}
        logging.info(f"Checking for {label} pockets...", extra=log_extra)

        qt_df = pd.concat([self.query_df, self.target_df], ignore_index=True).query(
            f"success and pocket_method == '{pocket_method}'"
        )
        if dedup_subset is not None:
            qt_df = qt_df.drop_duplicates(subset=dedup_subset)

        if len(qt_df) == 0:
            logging.info(f"No {label} pockets to retrieve", extra=log_extra)
        else:
            logging.info(f"{len(qt_df)} {label} pockets to retrieve", extra=log_extra)
        return qt_df

    def get_pockets(self):
        """
        Build every pocket in the run, dispatching each record to its pocket method.

        Fans out over the `builders` table rather than over hand-written calls: one row per
        `pocket_method`, naming the log label, the dedup columns and the builder, with the key
        doubling as the dump filename stem. Every builder returns the same Pocket shape (see
        `pocket.py`), so nothing downstream needs to know which method produced a given pocket.

        Returns:
            dict: pocket_id -> Pocket, across both sides.
        """
        log_extra = {"stage": "Getting Pockets"}
        logging.info("Starting pocket retrieval...", extra=log_extra)

        # Turns PDB Foldseek-database hits into ordinary pisa target records, so the retrieval below
        # picks them up like any other pisa entry. No-op for every other kind of target.
        self.expand_fsdb_pdb_targets()

        # pocket_method -> (log label, drop_duplicates subset, builder). Insertion order is the
        # merge order. A new pocket method is one row here plus its builder, and nothing else.
        builders = {
            "pisa": ("PISA", ["struct_info", "chain_info"], self.retrieve_pisa_pockets),
            "passthrough": ("passthrough", None, self.retrieve_passthrough_pockets),
            "vdw": ("VDW", None, self.retrieve_vdw_pockets),
            "whole_chain": ("whole chain", None, self.retrieve_whole_chain_pockets),
        }

        pockets = {}
        for pocket_method, (label, dedup_subset, builder) in builders.items():
            qt_df = self.select_pocket_records(pocket_method, label, dedup_subset)
            if len(qt_df) == 0:
                continue
            method_pockets = builder(qt_df)
            self.dump_pockets(method_pockets, f"{pocket_method}_pockets.json")
            logging.debug(f"Extracted {label} pockets: {method_pockets}", extra=log_extra)
            pockets |= method_pockets

        logging.debug(f"Combined pockets: {pockets}", extra=log_extra)
        return pockets

    def expand_fsdb_pdb_targets(self):
        """
        Turn the hits of a PDB Foldseek-database search into ordinary PISA target records.

        A Foldseek-database target has no per-chain records of its own, so `compare_pockets` normally
        synthesises a whole-chain pseudo-pocket for each hit and leaves every `pocket_2_*` column empty.
        The PDB database is built from real PDB entries, though, so its hits have real PISA interfaces:
        this reads the hit names out of the alignment table, resolves each to a PDB ID and chain, asks
        PISA which chains that chain touches, and appends one `pisa` record per interface to
        `self.target_df`. `retrieve_pisa_pockets` then handles them like any other pisa entry.

        The generated records carry the Foldseek entry name as their `preprocess_name` rather than the
        one `QTProcessor` derives, because that is the key alignments are stored under -- it is what
        joins these pockets back to their alignment rows and to their Foldseek transforms.

        Runs after `alignment`, so `alignment.tsv` exists. Does nothing unless the target is a Foldseek
        database whose entries are named in the PDB style; a database of anything else (human_domains)
        keeps the synthesised whole-chain pockets.

        Returns:
            None: appends to `self.target_df` and sets `self.fsdb_pdb_target`.
        """
        if not self.fsdb_target:
            return

        log_extra = {"stage": "Expanding Foldseek DB Targets"}

        alignment_df = pd.read_csv(self.settings.alignment_path, sep="\t", engine="c")
        hits = {}  # foldseek entry name -> (pdb_id, chain_id)
        for hit_name in alignment_df["target"].unique().tolist():
            resolved = parse_foldseek_pdb_entry_name(hit_name)
            if resolved is not None:
                hits[hit_name] = resolved
        if not hits:
            logging.info(
                "Foldseek database target is not a PDB database; keeping whole-chain target pockets",
                extra=log_extra,
            )
            return
        self.fsdb_pdb_target = True

        pdb_list = sorted({pdb_id for pdb_id, _ in hits.values()})
        logging.info(
            f"Retrieving PISA interfaces for {len(pdb_list)} PDB entries behind {len(hits)} Foldseek hits",
            extra=log_extra,
        )

        # Shares retrieve_pisa_pockets' cache, so its own call is a no-op for everything fetched here.
        interface_dir = self.download_pisa_interfaces(pdb_list)

        # Building one record per interface the hit chain takes part in
        parser = PisaParser()
        qtprocessor = QTProcessor(
            pdb_dir=self.settings.pdb_dir,
            alphafold_dir=self.settings.alphafold_dir,
            foldseek_preprocessed_structure_dir=self.settings.foldseek_preprocessed_structure_dir,
            fsdb_dir=self.settings.fsdb_dir,
        )
        records = []
        for hit_name, (pdb_id, chain_id) in hits.items():
            for partner in parser.get_interface_partners(pdb_id, chain_id, interface_dir):
                record = qtprocessor.parse_individual_qt(f"{pdb_id}:{chain_id}_{partner}", pocket_method="pisa")
                if record is None:
                    continue
                # The alignment is keyed by the Foldseek entry name, not by the name QTProcessor derives.
                # Nothing preprocesses these structures, so the preprocessing paths are meaningless here.
                record.preprocess_name = hit_name
                record.preprocess_path = None
                record.preprocess_path_gz = None
                records.append(asdict(record))
        if not records:
            logging.warning("No PISA interfaces found for any Foldseek hit", extra=log_extra)
            return

        # Fetching structures last, so only entries that actually produced a pocket are downloaded.
        # Hits whose structure can't be fetched come back marked success=False and are dropped here;
        # fetch_missing_structures raises only if not one of them could be fetched, which would leave
        # nothing to compare against at all.
        target_df = self.fetch_missing_structures("foldseek hit", pd.DataFrame(records))
        target_df = target_df.query("success")

        logging.info(
            f"Added {len(target_df)} PISA target pockets from {target_df['preprocess_name'].nunique()} Foldseek hits",
            extra=log_extra,
        )
        # Row 0 stays the database record itself -- foldseek_alignment and align_structs read its struct_path.
        self.target_df = pd.concat([self.target_df, target_df], ignore_index=True)

    def retrieve_pisa_pockets(self, pisa_df):
        """
        Build a Pocket per record from the PDBe PISA interface it names.

        Downloads (or reuses from cache) the PISA files for every PDB entry in `pisa_df`, takes each
        record's pocket to be the residues its own chain contributes to the named interface, then
        reads the CA coordinates for those residues out of the structure. `PisaParser` skips records
        whose entry or interface cannot be resolved, so the result may be smaller than `pisa_df`.

        Args:
            pisa_df (pandas.DataFrame): Records with `pocket_method == "pisa"`, as returned by
                `select_pocket_records`.

        Returns:
            dict: pocket_id -> Pocket.
        """
        log_extra = {"stage": "Retrieving PISA Pockets"}

        pisa_pdb_list = pisa_df["struct_info"].unique().tolist()
        logging.debug(f"PDBs for which to retrieve PISA pockets: {pisa_pdb_list}", extra=log_extra)
        interface_dir = self.download_pisa_interfaces(pisa_pdb_list)

        parser = PisaParser()
        pisa_pockets = parser.get_pockets_from_records(records=pisa_df.to_dict(orient="records"), in_dir=interface_dir)
        logging.debug(f"PISA pockets before coordinates: {pisa_pockets}", extra=log_extra)

        # PisaParser gives residue ids but no geometry; this second pass fills in seq_pos and the CA
        # coordinates on the same Pocket, which is what the comparison and superposition need.
        for _, row in pisa_df.iterrows():
            if row["pocket_id"] in pisa_pockets:
                domain_chain, _ = split_chain_info(row["chain_info"])
                pisa_pockets[row["pocket_id"]] = parse_pocket_from_struct(
                    struct=row["struct_path"],
                    chain_id=domain_chain,
                    pocket_residues=[int(x) for x in pisa_pockets[row["pocket_id"]].res_auth_ids],
                    pocket=pisa_pockets[row["pocket_id"]],
                )
        return pisa_pockets

    def retrieve_passthrough_pockets(self, pt_df):
        """
        Build a Pocket from the residue ids the entry names outright.

        A passthrough entry carries its own residue list ("P24941:A:160,161"), so there is nothing to
        compute: the ids are read from `residue_info` and looked up in the structure. They are sorted
        numerically rather than kept in the order they were typed -- `res_auth_ids` is the order the
        comparison pairs residues against the other pocket in, and every other pocket method produces
        it ascending.

        An entry naming a residue the chain cannot supply is skipped rather than compared, so every
        pocket returned here has a `residues` entry for every id in its `res_auth_ids`.

        Args:
            pt_df (pandas.DataFrame): Records with `pocket_method == "passthrough"`, as returned by
                `select_pocket_records`.

        Returns:
            dict: pocket_id -> Pocket, possibly smaller than `pt_df`.
        """
        log_extra = {"stage": "Retrieving passthrough Pockets"}

        passthrough_pockets = {}
        for _, row in pt_df.iterrows():
            domain_chain, _ = split_chain_info(row["chain_info"])
            pocket = parse_pocket_from_struct(
                struct=row["struct_path"],
                chain_id=domain_chain,
                pocket_residues=sorted(int(x) for x in row["residue_info"].split(",")),
            )
            # A missing structure or chain gives None back, which would fail later with an opaque
            # TypeError inside compare_pockets.
            if pocket is None:
                logging.warning(
                    f"Could not parse chain {domain_chain} of {row['struct_info']} for {row['pocket_id']}, "
                    "skipping this entry",
                    extra=log_extra,
                )
                continue
            # residues holds exactly the requested ids the chain walk reached with CA coordinates, so
            # anything else left in res_auth_ids is an id the comparison would raise a KeyError on.
            unusable = [res_id for res_id in pocket.res_auth_ids if res_id not in pocket.residues]
            if unusable:
                logging.warning(
                    f"Residue(s) {','.join(unusable)} of {row['pocket_id']} are not in chain {domain_chain} "
                    "or have no CA atom, skipping this entry",
                    extra=log_extra,
                )
                continue
            passthrough_pockets[row["pocket_id"]] = pocket
        return passthrough_pockets

    def retrieve_vdw_pockets(self, vdw_df):
        """
        Build a Pocket from the van der Waals contacts between the entry's two chains.

        A chain pair on a structure PISA cannot serve -- a local file or an AlphaFold model --
        resolves to this method. The pocket is the residues of the first chain whose atoms come
        within van der Waals contact of the second; see `PocketCalculator.pocket_overlap`.

        Args:
            vdw_df (pandas.DataFrame): Records with `pocket_method == "vdw"`, as returned by
                `select_pocket_records`.

        Returns:
            dict: pocket_id -> Pocket.
        """
        vdw_pockets = {}
        pc = PocketCalculator()
        for _, row in vdw_df.iterrows():
            domain_chain, motif_chain = split_chain_info(row["chain_info"])
            vdw_pockets[row["pocket_id"]] = pc.pocket_overlap(
                structure=row["struct_path"],
                domain_chain=domain_chain,
                motif_chain=motif_chain,
            )
        return vdw_pockets

    def retrieve_whole_chain_pockets(self, wc_df):
        """
        Build the "pocket" for an open search: every CA-bearing residue of the chain.

        An entry that names a structure but no pocket ("4Q5J:B", or "4Q5J" for the default chain) asks
        whether the query pocket resembles anything on that chain at all. Unlike the whole-chain pseudo-pocket
        `compare_pockets` synthesises for Foldseek-database hits, this is a real Pocket parsed from the
        structure, so it carries residue codes and CA coordinates and can be superposed.

        Args:
            wc_df (pandas.DataFrame): Records with `pocket_method == "whole_chain"`, as returned by
                `select_pocket_records`.

        Returns:
            dict: pocket_id -> Pocket.
        """
        log_extra = {"stage": "Retrieving whole chain Pockets"}

        whole_chain_pockets = {}
        for _, row in wc_df.iterrows():
            domain_chain, _ = split_chain_info(row["chain_info"])
            pocket = parse_pocket_from_struct(
                struct=row["struct_path"],
                chain_id=domain_chain,
                pocket_residues=None,  # None means the whole chain
            )
            # A missing structure or chain gives None back. Storing that would fail later with an opaque
            # TypeError inside compare_pockets, so drop the entry and say which one it was -- an unreadable
            # chain is much more likely here, where the chain can come from the default rather than the user.
            if pocket is None:
                logging.warning(
                    f"Could not parse chain {row['chain_info']} of {row['struct_info']} for {row['pocket_id']}, "
                    "skipping this entry",
                    extra=log_extra,
                )
                continue
            whole_chain_pockets[row["pocket_id"]] = pocket
        return whole_chain_pockets

    def compare_pockets_based_on_alignment(self, pockets):
        """
        Merge sequence or fold structure alignments traversing spatial pocket metrics internally to evaluate correlations.

        Processes raw mapping dictionaries tracking unresolved structural aliases. Utilizes a referenced BLOSUM
        scorecard for computing structural compatibility values output into a final result dataframe file.

        Args:
            pockets (dict): Populated pool of pocket residue locations collected via initial configurations.

        Returns:
            None
        """
        log_extra = {"stage": "Comparing Pockets Based on Alignment"}

        logging.info("Reading alignment results...", extra=log_extra)
        alignment_df = pd.read_csv(self.settings.alignment_path, sep="\t", engine="c")
        blosum_path = os.path.join(os.path.dirname(__file__), "blosum62.bla")

        logging.info(f"{len(alignment_df)} alignment pairs to compare", extra=log_extra)
        logging.debug(f"Alignment pairs: \n{alignment_df.head()}", extra=log_extra)

        preproc_to_ids = {}
        for _, row in pd.concat([self.query_df, self.target_df], ignore_index=True).iterrows():
            # Keyed by preprocess_name: one chain can carry several pockets (several pocket_ids), so
            # membership must be tested on the key, not on the pocket_id, or each chain keeps only its last.
            if row["preprocess_name"] in preproc_to_ids:
                if row["pocket_id"] not in preproc_to_ids[row["preprocess_name"]]:
                    preproc_to_ids[row["preprocess_name"]].append(row["pocket_id"])
            else:
                preproc_to_ids[row["preprocess_name"]] = [row["pocket_id"]]
        logging.debug(f"Preprocessed name to pocket ID mapping: {preproc_to_ids}", extra=log_extra)

        # A PDB Foldseek database has real PISA pockets for its hits (see expand_fsdb_pdb_targets);
        # every other database has no target records at all, so its pockets must be synthesised.
        synthesise_target_pockets = self.fsdb_target and not self.fsdb_pdb_target

        # Read from target_df rather than target_db_path so this does not depend on a step-4 side
        # effect. Row 0 is still the single database record here: expand_fsdb_pdb_targets only rewrites
        # target_df on the PDB path, which synthesise_target_pockets excludes.
        offset_table_path = None
        if synthesise_target_pockets:
            # Matched on the resolved path, so naming the bundled DB by its path and naming it
            # "human_domains" give the same answer.
            offset_paths = {
                db["db_path"]: db["offset_path"] for db in bundled_foldseek_dbs(self.settings.fsdb_dir).values()
            }
            offset_table_path = offset_paths.get(self.target_df.loc[0, "struct_path"])
            if offset_table_path is None:
                logging.info(
                    "Foldseek database ships no offset table; target residue ids will be positions "
                    "within each database entry rather than UniProt coordinates",
                    extra=log_extra,
                )
            elif not os.path.isfile(offset_table_path):
                # Packaged with the database (pyproject.toml's human_domains/*), so absence means a
                # broken install. Falling back would emit positional ids under a column that the
                # bundled database is documented to report in UniProt coordinates.
                msg = f"Bundled human-domains offset table is missing from the installation: {offset_table_path}"
                logging.critical(msg, extra=log_extra)
                raise PocketMapperError(msg)
            else:
                logging.info(
                    f"Mapping target residue ids to UniProt coordinates using {offset_table_path}", extra=log_extra
                )

        pockets_df, unknown_alias, incorrect_mapping = compare_pockets(
            alignment_df,
            pockets,
            preproc_to_ids=preproc_to_ids,
            blosum_path=blosum_path,
            synthesise_target_pockets=synthesise_target_pockets,
            offset_table_path=offset_table_path,
        )

        # Logging cases where a residue was given a single cahr name unfamiliar to pocketmapper
        if len(unknown_alias) > 0:
            unknown_alias_path = os.path.join(self.settings.results_dir, "unknown_ids.json")
            logging.warning("Unknown Foldseek Alias, see unknown_alias.json in results directory", extra=log_extra)
            with open(unknown_alias_path, "w") as f:
                json.dump(jsonify_dict(dict(unknown_alias)), f)

        # logging cases where foldseek mapping had low sequence identity to the parsed structure
        if len(incorrect_mapping) > 0:
            incorrect_mapping_path = os.path.join(self.settings.results_dir, "incorrect_mapping.json")
            logging.warning("Foldseek mapping with low sequence identity to parsed structure", extra=log_extra)
            with open(incorrect_mapping_path, "w") as f:
                json.dump(jsonify_dict(dict(incorrect_mapping)), f)

        # Writing pocket comparison results to output file
        output_path = self.settings.pocket_comparison_path
        pockets_df.to_csv(output_path, index=False, sep="\t")
        logging.info(f"Pocket comparison results saved to {output_path}", extra=log_extra)

    def align_structs(self):
        """
        Perform structural superposition of target structures against the query reference frame.

        Unpacks this run's settings and target shape for `StructureAligner.align_structs`, which takes
        the top `align_count` targets of each query from the pocket comparison results and writes them
        superposed onto it into `aligned_structure_dir`.

        Returns:
            None
        """
        # A Foldseek-database target has no structures of its own, so the aligner rebuilds them from
        # the database. Which records it gets says how to read a target id: the PDB database's hits
        # were expanded into real pockets by expand_fsdb_pdb_targets, while for any other database the
        # ids are database entry names and target_df holds only the database record itself.
        fsdb_path = self.target_df.loc[0, "struct_path"] if self.fsdb_target else None
        if self.fsdb_target and not self.fsdb_pdb_target:
            target_records = []
        else:
            target_records = self.target_df.to_dict(orient="records")

        aligner = StructureAligner()
        aligner.align_structs(
            query_records=self.query_df.to_dict(orient="records"),
            target_records=target_records,
            pocket_comparison=self.settings.pocket_comparison_path,
            out_dir=self.settings.aligned_structure_dir,
            method=self.settings.align_struct_method,  # already "pocket" or "foldseek"
            align_count=self.settings.align_count,
            alignment=self.settings.alignment_path,
            threads=self.settings.threads,
            fsdb_path=fsdb_path,
        )

    def delete_tmp(self):
        """
        Delete this run's scratch directory, unless `delete_tmp` says otherwise.

        A temp_dir resolving outside cache_dir and results_dir is left alone and warned about: it is
        settable by both the settings file and the command line, so a mistyped `--temp_dir` would
        otherwise hand an unrelated directory to `shutil.rmtree`.

        Returns:
            None
        """
        log_extra = {"stage": "Cleaning Up"}

        path = self.settings.temp_dir

        # Named rather than announced: a run kept for inspection should say where its inputs are.
        if not self.settings.delete_tmp:
            logging.info(f"delete_tmp is False; keeping {path}", extra=log_extra)
            return

        if not is_within(path, [self.settings.cache_dir, self.settings.results_dir]):
            logging.warning(
                f"Not deleting temp_dir {path}: it is outside cache_dir and results_dir. "
                "Remove it yourself if that was intended.",
                extra=log_extra,
            )
            return
        shutil.rmtree(path)
