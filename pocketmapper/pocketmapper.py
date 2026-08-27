"""
PocketMapper: A tool for mapping and analyzing protein pockets.

Author: Lachlan Ellingboe

"""

from dataclasses import asdict, dataclass, field, replace
import fire
import logging
import logging.config
import json
import subprocess
import pandas as pd
import os
import sys
from datetime import datetime
import shutil
from pocketmapper.lib import jsonify_dict, parse_foldseek_pdb_entry_name, safe_filename
from pocketmapper.exceptions import PocketMapperError
from pocketmapper.pisa_downloader import PisaDownloader
from pocketmapper.pisa_parser import PisaParser
from pocketmapper.pocket_comparison import compare_pockets, parse_pocket_transform
from pocketmapper.sequence_aligner import SequenceAligner
from pocketmapper.pocket_calculator import PocketCalculator
from pocketmapper.qt_processor import QTProcessor
from pocketmapper.structure_aligner import StructureAligner
from pocketmapper.structure_fetcher import StructureFetcher
from pocketmapper.structure_preprocessor import StructurePreprocessor
from pocketmapper.lib_struct import parse_pocket_from_struct
from pocketmapper.constants import (
    ALIGN_STRUCT_METHODS,
    FOLDSEEK_FORMAT_OUTPUT,
    FOLDSEEK_INSTALL_HINT,
    HELP_MESSAGE,
)


@dataclass
class Settings:
    """
    Fully resolved PocketMapper run configuration.

    Built by layering three sources in priority order (lowest to highest):
    dataclass defaults -> JSON settings file -> explicit CLI arguments.
    Derived paths (structure_dir, alignment_path, etc.) are filled in
    afterward by resolve_paths(), unless already set by the settings file.
    """

    query: str | None = None
    target: str | None = None
    cache_dir: str = "pocketmapper_cache"
    results_dir: str = field(default_factory=lambda: f"pocketmapper_results_{datetime.now().strftime('%y%m%d_%H%M%S')}")
    query_pocket_method: str | None = None
    target_pocket_method: str | None = None
    # Tri-state: None (the default) means "auto" -- use Foldseek when the binary is on PATH and
    # fall back to the local aligner when it is not. _resolve_foldseek() turns this into a concrete
    # bool before anything else reads it, so the rest of the pipeline only ever sees True/False.
    foldseek: bool | None = None
    align_count: int = 10
    # Which transform superposes a target onto its query in step 7: "foldseek" (Foldseek's whole-chain
    # fit) or "pocket" (the fit of the two pockets on their overlapping residues). The default "auto"
    # is collapsed to one of those by _resolve_align_struct_method(), so nothing downstream sees it.
    align_struct_method: str = "auto"
    verbosity: int = 3

    # Derived paths -- left unset (None) until resolve_paths() fills them in,
    # unless explicitly provided via the settings file.
    structure_dir: str | None = None
    pocket_dir: str | None = None
    foldseek_tmp_dir: str | None = None
    foldseek_preprocessed_structure_dir: str | None = None
    query_dir: str | None = None
    target_dir: str | None = None
    aligned_structure_dir: str | None = None
    alignment_path: str | None = None
    pocket_comparison_path: str | None = None
    job_settings_path: str | None = None
    log_path: str | None = None
    fsdb_dir: str | None = None

    def resolve_paths(self):
        """
        Return a copy of these settings with any unset derived paths filled
        in from cache_dir/results_dir. Paths already set (e.g. via the
        settings file) are left untouched.
        """
        derived = {
            "structure_dir": os.path.join(self.cache_dir, "ref_structures"),
            "pocket_dir": os.path.join(self.cache_dir, "pockets"),
            "foldseek_tmp_dir": os.path.join(self.cache_dir, "foldseek_tmp"),
            "foldseek_preprocessed_structure_dir": os.path.join(self.cache_dir, "foldseek_preprocessed_structures"),
            "query_dir": os.path.join(self.results_dir, "query_structures"),
            "target_dir": os.path.join(self.results_dir, "target_structures"),
            "aligned_structure_dir": os.path.join(self.results_dir, "aligned_structures"),
            "alignment_path": os.path.join(self.results_dir, "alignment.tsv"),
            "pocket_comparison_path": os.path.join(self.results_dir, "pocket_comparison.tsv"),
            "job_settings_path": os.path.join(self.results_dir, "job_settings.json"),
            "log_path": os.path.join(self.results_dir, "info.log"),
            "fsdb_dir": os.path.join(self.cache_dir, "fsdb"),
        }
        unset = {key: path_val for key, path_val in derived.items() if getattr(self, key) is None}
        return replace(self, **unset)


class PocketMapper:
    def __init__(self):
        """
        Initialize the PocketMapper instance.

        Sets up the default logging configuration, formatting strings, and points to
        bundled structural databases if applicable.
        """
        self._log_extra = {"stage": "init"}
        self._log_fmt = "%(levelname)s: %(stage)s - %(msg)s"
        logging.getLogger(__name__)
        logging.basicConfig(level=logging.CRITICAL, format=self._log_fmt)

        self.fsdb_target = False
        self._fsdb_pdb_target = False
        # Set for real by _resolve_foldseek(); read by _configure_query_target to explain *why*
        # a Foldseek-DB target was rejected. True here so a caller that skips search() is not
        # told the binary is missing when nothing has looked for it.
        self._foldseek_available = True

    def _configure_logging(self, settings):
        """
        Configure logging level and handlers based on user settings.

        Sets both a console stream handler and a file handler (to `log_path` in `settings`),
        adjusting the verbosity depending on the user's input.

        Args:
            settings (Settings): Configuration object with `verbosity` and `log_path` attributes.
                verbosity: 4=DEBUG, 3=INFO, 2=WARNING, else ERROR.
        """
        self._log_extra.update({"stage": "Configuring Logging"})

        # Set log level based on verbosity setting (default to INFO if not set)
        log_level = None
        if settings.verbosity == 4:
            log_level = "DEBUG"
        elif settings.verbosity == 3:
            log_level = "INFO"
        elif settings.verbosity == 2:
            log_level = "WARNING"
        else:
            log_level = "ERROR"

        log_config = {
            "version": 1,
            "formatters": {
                "standard": {"format": self._log_fmt},
            },
            "handlers": {
                "console": {
                    "class": "logging.StreamHandler",
                    "level": log_level,
                    "formatter": "standard",
                    "stream": "ext://sys.stdout",
                },
                "file": {
                    "class": "logging.FileHandler",
                    "level": log_level,
                    "formatter": "standard",
                    "filename": settings.log_path,
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
        query=None,  # settings passed to configure
        target=None,
        settings=None,
        cache_dir=None,
        results_dir=None,
        verbosity=None,
        help=None,
        foldseek=None,
        align_count=None,
        align_struct_method=None,
        query_pocket_method=None,
        target_pocket_method=None,
    ):
        """
        Orchestrate and run the full PocketMapper search workflow.

        Args:
            query (str): Target query identifier, string or path to a list.
            target (str): Target structure identifier, string or path to a list.
            settings (str, optional): Path to a JSON settings file.
            cache_dir (str, optional): Directory to cache intermediate structures.
            results_dir (str, optional): Directory to output results to.
            verbosity (int, optional): Control logging level.
            help (bool, optional): Output the help message and exit.
            foldseek (bool, optional): Use foldseek for structure alignment instead of local sequence
                alignment. Left unset, foldseek is used when the binary is on PATH and the local
                aligner is used with a warning when it is not. True makes foldseek a hard
                requirement -- a missing binary is an error; False always uses the local aligner.
            align_count (int, optional): Number of top targets to superpose onto each query.
            align_struct_method (str, optional): Which transform superposes a target onto its query --
                'foldseek' for Foldseek's whole-chain fit, 'pocket' for the fit of the two pockets
                on their overlapping residues, or 'auto' (the default) for 'foldseek' when Foldseek is
                in use and 'pocket' with the local aligner, which produces no chain transform at all.
        """
        self._log_extra = {
            "stage": "Starting Search"
        }  # dict needed for logging extra info, can be updated throughout the process to indicate the current stage in logs

        # Storing input parameters
        self._query = query
        self._target = target
        self._settings_file = settings
        self._cache_dir = cache_dir
        self._results_dir = results_dir
        self._verbosity = verbosity
        self._help = help
        self._foldseek = foldseek
        self._align_count = align_count
        self._align_struct_method = align_struct_method
        self._query_pocket_method = query_pocket_method
        self._target_pocket_method = target_pocket_method

        self._check_help_search()  # Checks if help flag is set and if so prints the help message and exits
        self._settings = self._configure_workflow()  # configures the settings which have already been read
        self._query_df, self._target_df = (
            self._configure_query_target()
        )  # parses the query and target inputs to determine their types and sets up the relevant data structures for each entry

        self._query_df = self._fetch_missing_structures(
            "query", self._query_df, self._settings.structure_dir
        )  # Fetch any missing structures
        if self.fsdb_target:
            self._fetch_missing_fsdb(
                self._target_df, self._settings.foldseek_tmp_dir
            )  # Fetch any missing foldseek databases
        else:
            self._target_df = self._fetch_missing_structures(
                "target", self._target_df, self._settings.structure_dir
            )  # Fetch any missing structures

        self._alignment()  # Align the query and target structures using either local sequence alignment or foldseek based on the settings
        pockets = self._get_pockets()  # Adds seq_pos and ca-coords to the pocket info dict
        self._compare_pockets_based_on_alignment(pockets)
        self._align_structs()
        self._delete_tmp()

        logging.info("PocketMapper search completed successfully.", extra={"stage": "End"})

    def _check_help_search(self):
        """
        Display help information for the PocketMapper tool and exit the program.

        If the 'self._help' parameter is provided and evaluates to True, this method prints a
        help message describing the usage, options, and features of the PocketMapper package,
        then terminates execution.

        Returns:
            None: Process exits if `self._help` is True.
        """

        if self._help:
            print(HELP_MESSAGE)
            exit()

    def _configure_workflow(self):
        """
        Set up and evaluate configuration settings for the workflow.

        This process:
        1. Sets base defaults for caching and results.
        2. Overrides these defaults via an optional JSON settings file.
        3. Prioritizes CLI arguments (passed via the `search` method) over file settings.
        4. Computes all derivative working directories needed during pipeline execution.
        5. Saves final configuration payload locally, setting up necessary directories.

        Returns:
            None
        """
        self._log_extra.update({"stage": "Configuring Settings"})

        # 1. Base defaults
        settings = Settings()

        # 2. Populate settings from the settings file if provided
        if self._settings_file is not None:
            if not os.path.isfile(self._settings_file):
                logging.critical(f"Settings file not found: {self._settings_file}", extra=self._log_extra)
                raise PocketMapperError(f"Settings file not found: {self._settings_file}")
            try:
                with open(self._settings_file) as f:
                    settings_data_from_file = json.load(f)
                settings = replace(settings, **settings_data_from_file)
            except TypeError as e:
                logging.critical(f"Unknown setting(s) in {self._settings_file}: {e}", extra=self._log_extra)
                raise PocketMapperError(f"Unknown setting(s) in {self._settings_file}: {e}") from e
            except Exception as e:
                logging.critical(
                    f"Error reading settings file: {self._settings_file}. Is it in JSON format?", extra=self._log_extra
                )
                raise PocketMapperError(
                    f"Error reading settings file: {self._settings_file}. Is it in JSON format?"
                ) from e

        # 3. Override settings with explicit command-line arguments
        cli_overrides = {
            "query": self._query,
            "target": self._target,
            "cache_dir": self._cache_dir,
            "results_dir": self._results_dir,
            "foldseek": self._foldseek,
            "verbosity": self._verbosity,
            "align_count": self._align_count,
            "align_struct_method": self._align_struct_method,
            "query_pocket_method": self._query_pocket_method,
            "target_pocket_method": self._target_pocket_method,
        }
        cli_overrides = {key: value for key, value in cli_overrides.items() if value is not None}
        settings = replace(settings, **cli_overrides)

        # 4. Computed paths (only fills in paths not already set via the settings file)
        settings = settings.resolve_paths()

        # Ensure all necessary directories exist before proceeding, creating them if needed
        dirs_to_create = [
            "structure_dir",
            "query_dir",
            "target_dir",
            "pocket_dir",
            "foldseek_preprocessed_structure_dir",
            "aligned_structure_dir",
        ]
        for dir_key in dirs_to_create:
            path = getattr(settings, dir_key)
            try:
                os.makedirs(path, exist_ok=True)
            except OSError as e:
                logging.critical(f"Error creating directory {path}", extra=self._log_extra)
                raise PocketMapperError(f"Error creating directory {path}") from e

        self._configure_logging(settings)

        # 4b. Resolve the tri-state foldseek setting into a concrete bool. Must come after
        # _configure_logging (the root logger is still at CRITICAL before it, so the fallback
        # warning would be swallowed) and before the settings are logged and dumped below, so
        # job_settings.json records what the run actually did.
        settings = self._resolve_foldseek(settings)

        # 4c. Same reasoning, and it reads the bool _resolve_foldseek just settled, so it must follow it.
        settings = self._resolve_align_struct_method(settings)

        logging.info(f"Settings: {json.dumps(asdict(settings), indent=4)}", extra=self._log_extra)

        # 5. Output dump
        try:
            os.makedirs(os.path.dirname(settings.job_settings_path), exist_ok=True)
            with open(settings.job_settings_path, "w") as f:
                json.dump(asdict(settings), f, indent=4)
            logging.info(f"Settings successfully dumped to {settings.job_settings_path}", extra=self._log_extra)
        except Exception as e:
            logging.error(f"Failed to dump settings to {settings.job_settings_path}: {e}", extra=self._log_extra)
        return settings

    def _resolve_foldseek(self, settings):
        """
        Turn the tri-state `foldseek` setting into a concrete bool.

        Foldseek is an optional external binary, so the default (None, "auto") is resolved against
        what is actually installed: foldseek when it is on PATH, the local BLOSUM62 aligner with a
        warning when it is not. An explicit True is a hard requirement and errors instead of falling
        back; an explicit False always means the local aligner and never probes for the binary.

        Called before any structure is fetched, so an unmet requirement fails without wasted
        downloads rather than as a raw FileNotFoundError from the first `foldseek` subprocess call.

        Args:
            settings (Settings): Settings whose `foldseek` field may still be None.

        Returns:
            Settings: A copy with `foldseek` set to True or False.

        Raises:
            PocketMapperError: If foldseek was explicitly requested but is not installed.
        """
        # Local stage dict: _configure_logging leaves self._log_extra reading "Configuring Logging".
        stage = {"stage": "Configuring Settings"}

        if settings.foldseek is False:
            return settings

        self._foldseek_available = shutil.which("foldseek") is not None
        if self._foldseek_available:
            return replace(settings, foldseek=True)

        if settings.foldseek is True:
            msg = (
                "Foldseek alignment was requested but 'foldseek' was not found on PATH. "
                f"{FOLDSEEK_INSTALL_HINT} Alternatively, set --foldseek False to use the local "
                "BLOSUM62 sequence aligner."
            )
            logging.critical(msg, extra=stage)
            raise PocketMapperError(msg)

        logging.warning(
            "'foldseek' not found on PATH; falling back to the local BLOSUM62 sequence aligner. "
            "The local aligner produces no whole-chain transform, so aligned_structures/*.pdb are "
            f"superposed on the pocket instead (see --align_struct_method). {FOLDSEEK_INSTALL_HINT}",
            extra=stage,
        )
        return replace(settings, foldseek=False)

    def _resolve_align_struct_method(self, settings):
        """
        Turn the tri-value `align_struct_method` setting into "pocket" or "foldseek".

        "foldseek" uses Foldseek's whole-chain transform from alignment.tsv; "pocket" uses the
        superposition of the two pockets on their overlapping residues, which `compare_pockets`
        already writes to pocket_comparison.tsv. The default "auto" picks whichever the run can
        actually do: the local BLOSUM62 aligner writes "-" for the chain transform, so it has only
        the pocket one.

        An explicit "foldseek" without the binary is an error rather than a silent switch to "pocket" --
        the same call as an unmet `--foldseek True`, and for the same reason: better to fail before
        anything is downloaded than to hand back a method the user did not ask for.

        Called from _configure_workflow after _resolve_foldseek, whose resolved bool it reads, and
        before the settings are logged and dumped, so job_settings.json records what the run did.

        Args:
            settings (Settings): Settings whose `foldseek` is already a concrete bool.

        Returns:
            Settings: A copy with `align_struct_method` set to "pocket" or "foldseek".

        Raises:
            PocketMapperError: If the value is not one of ALIGN_STRUCT_METHODS, or "foldseek" was
                asked for on the local-aligner path.
        """
        stage = {"stage": "Configuring Settings"}

        method = settings.align_struct_method
        # fire hands over whatever was typed, and a settings file can hold anything at all.
        method = method.lower() if isinstance(method, str) else method
        if method not in ALIGN_STRUCT_METHODS:
            msg = (
                f"Unknown align_struct_method {settings.align_struct_method!r}. "
                f"Choose one of: {', '.join(ALIGN_STRUCT_METHODS)}."
            )
            logging.critical(msg, extra=stage)
            raise PocketMapperError(msg)

        if method == "auto":
            method = "foldseek" if settings.foldseek else "pocket"
            logging.info(
                f"align_struct_method 'auto' resolved to '{method}' "
                f"({'foldseek' if settings.foldseek else 'the local aligner'} is in use)",
                extra=stage,
            )
        elif method == "foldseek" and not settings.foldseek:
            msg = (
                "align_struct_method 'foldseek' needs Foldseek's whole-chain transform, but this run "
                "uses the local BLOSUM62 aligner, which does not produce one. Use "
                "--align_struct_method pocket, or enable foldseek."
            )
            logging.critical(msg, extra=stage)
            raise PocketMapperError(msg)

        return replace(settings, align_struct_method=method)

    def _configure_query_target(self):
        """
        Determine and process input data (formats and types) for the query and target constraints.

        Uses the `QTProcessor` class to load, parse, and validate query vs. target identities
        and requested pocket methodologies (e.g. "pisa", "passthrough"). Updates local dataframes.

        Returns:
            None: Instantiates `self._query_df` and `self._target_df`.
        """
        self._log_extra.update({"stage": "Determine Query/Target Types"})

        qtprocessor = QTProcessor(
            structure_dir=self._settings.structure_dir,
            foldseek_preprocessed_structure_dir=self._settings.foldseek_preprocessed_structure_dir,
            fsdb_dir=self._settings.fsdb_dir,
        )
        q_df = qtprocessor.process_qt_cmdline_input(
            qt_input=self._settings.query,
            name="query",
            pocket_method=self._settings.query_pocket_method,
        )
        t_df = qtprocessor.process_qt_cmdline_input(
            qt_input=self._settings.target,
            name="target",
            pocket_method=self._settings.target_pocket_method,
        )

        errors = []
        if len(q_df) < 1:
            logging.critical("No valid query entries after processing", extra=self._log_extra)
            errors.append("no valid query entries")
        if len(t_df) < 1:
            logging.critical("No valid target entries after processing", extra=self._log_extra)
            errors.append("no valid target entries")
        if errors:
            raise PocketMapperError("; ".join(errors))

        if t_df.loc[0, "struct_type"] == "foldseek_db":
            if self._settings.foldseek:
                self.fsdb_target = True
                # Neither kind of Foldseek DB can be superposed on its pocket, and which kind this is
                # is not known until _expand_fsdb_pdb_targets has read the hit names -- so reject both
                # here, before anything is fetched. Unreachable from "auto": a DB target forces
                # foldseek on, and auto resolves to "foldseek" whenever it is on.
                if self._settings.align_struct_method == "pocket":
                    msg = (
                        "align_struct_method 'pocket' is not available against a Foldseek database "
                        "target. A human_domains-style hit has no coordinates to superpose at all, and "
                        "a PDB database's structures are assemblies while its pockets come from the "
                        "wwPDB asymmetric unit, so a pocket fit would be applied in the wrong frame. "
                        "Use --align_struct_method foldseek."
                    )
                    logging.critical(msg, extra=self._log_extra)
                    raise PocketMapperError(msg)
            else:
                # foldseek is already resolved to a concrete bool here, so False means either the
                # binary is missing or the user turned it off -- say which, since the fixes differ.
                if not self._foldseek_available:
                    msg = (
                        "A Foldseek database was specified as the target, which requires the "
                        f"'foldseek' binary, but it was not found on PATH. {FOLDSEEK_INSTALL_HINT}"
                    )
                else:
                    msg = (
                        "Foldseek database specified as target but foldseek is not enabled. "
                        "Remove --foldseek False to use it."
                    )
                logging.critical(msg, extra=self._log_extra)
                raise PocketMapperError(msg)
        return q_df, t_df

    def _fetch_missing_structures(self, name, qt_df, out_dir):
        """
        Identify and download required structure files.

        Iterates over `self._query_df` and `self._target_df` determining if structure
        components are locally accessible, utilizing `StructureFetcher`.
        Failed retrieval flags will trigger program termination if zero data persists.

        Returns:
            None
        """
        # Downloading structures
        self._log_extra.update({"stage": "Fetching Missing Structures"})
        logging.info("Starting", extra=self._log_extra)
        structure_fetcher = StructureFetcher()

        logging.debug(f"{name.capitalize()} data before fetching structures: \n{qt_df.head()}", extra=self._log_extra)

        # Get list of unique structures to fetch based on struct_info and struct_type
        unique_records = qt_df.drop_duplicates(subset="struct_info").to_dict(orient="records")

        # Update structure fetcher and fetch structures
        structure_fetcher.set_output_directory(out_dir)
        structure_fetcher.update_cache()
        results = structure_fetcher.fetch_structures(unique_records)
        logging.debug(f"Structure fetcher results: {results}", extra=self._log_extra)

        # Update the dataframe with success/failure information
        qt_df["success"] = qt_df["struct_info"].map(results).fillna(False)
        qt_df.loc[~qt_df["success"], "failure_reason"] = "structure_not_found"

        # Logging results of structure fetching and updating query and target data with success/failure info
        logging.info(
            f"{sum(results.values())}/{len(results)} {name} required structures available",
            extra=self._log_extra,
        )
        if len(qt_df.query("success == False")) > 0:
            logging.warning(
                f"Missing structures for {name}(s): {', '.join(qt_df.loc[~qt_df['success'], 'pocket_id'].unique().tolist())}",
                extra=self._log_extra,
            )

        # Verifying sufficient structures were found to continue
        if qt_df["success"].sum() < 1:
            logging.critical(f"Insufficient {name} structures after fetching", extra=self._log_extra)
            raise PocketMapperError(f"Insufficient {name} structures after fetching. No valid {name} entries remain.")
        return qt_df

    def _fetch_missing_fsdb(self, qt_df, tmp_dir):
        self._log_extra.update({"stage": "Fetching Missing Foldseek Database"})
        fsdb_name = qt_df.loc[0, "struct_info"].upper()
        fsdb_path = qt_df.loc[0, "struct_path"]
        if not os.path.exists(fsdb_path):
            logging.info(f"Fetching bundled Foldseek database '{fsdb_name}' to {fsdb_path}", extra=self._log_extra)
            try:
                os.makedirs(os.path.dirname(fsdb_path), exist_ok=True)
                cmd = ["foldseek", "databases", fsdb_name, fsdb_path, tmp_dir]
                logging.debug(f"Running command: {' '.join([str(x) for x in cmd])}", extra=self._log_extra)
                subprocess.run(cmd, check=True)
                logging.info(f"Successfully fetched Foldseek database '{fsdb_name}'", extra=self._log_extra)
            except Exception as e:
                logging.critical(
                    f"Failed to fetch Foldseek database '{fsdb_name}' to {fsdb_path}: {e}", extra=self._log_extra
                )
                raise PocketMapperError(f"Failed to fetch Foldseek database '{fsdb_name}' to {fsdb_path}: {e}") from e

    def _alignment(self):
        """
        Coordinate structural alignment routes bridging query and target proteins.

        Dispatches to `_foldseek_alignment()` if the foldseek flag is set, else rolls
        back to `_local_alignment()`. Requisites like `_foldseek_preprocessing()`
        precede foldseek routines.

        Returns:
            None
        """
        log_extra = {"stage": "Alignment"}
        if self._settings.foldseek:
            logging.info("Preprocessing structures for Foldseek...", extra=log_extra)
            self._foldseek_preprocessing()
            logging.info("Running Foldseek easy-search...", extra=log_extra)
            self._foldseek_alignment()
        else:
            logging.info("Running local pairwise aligner...", extra=log_extra)
            self._local_alignment()

    def _foldseek_preprocessing(self):
        """
        Split complex PDB structures into single chain mmCIF files for Foldseek processing.
        Copies structure to the query/target directories, ensuring that the Foldseek preprocessed structure directory is populated.

        Instantiates a `StructurePreprocessor` mapping items sourced from the
        local datastores towards isolated files targeting specific sequences.

        Returns:
            None
        """
        stage = {"stage": "Preprocessing Structures"}

        structure_preprocessor = StructurePreprocessor()
        qtdf_dir_iter = [(self._query_df, self._settings.query_dir)]
        if not self.fsdb_target:
            qtdf_dir_iter.append((self._target_df, self._settings.target_dir))

        for df, search_dir in qtdf_dir_iter:
            records = df.drop_duplicates(subset=["preprocess_name", "chain_info"]).to_dict(orient="records")
            logging.debug(f"Records to preprocess: {json.dumps(records, indent=4)}", extra=stage)

            structure_preprocessor.set_output_directory(self._settings.foldseek_preprocessed_structure_dir)
            structure_preprocessor.update_cache()
            results = structure_preprocessor.preprocess_records(records=records, search_dir=search_dir)
            logging.debug(f"Preprocessing results: {json.dumps(results, indent=4)}", extra=stage)

            # Updating success and failure cols based on preprocessing results
            df.set_index("pocket_id", inplace=True)
            for index, success in results.items():
                if not success:
                    df.loc[index, "success"] = False
                    df.loc[index, "failure_reason"] = "structure_preprocessing_failed"
            df.reset_index(inplace=True)

        logging.info("Finished preprocessing structures", extra=stage)
        logging.debug(f"Query data after preprocessing: \n{self._query_df.head()}", extra=stage)
        logging.debug(f"Target data after preprocessing: \n{self._target_df.head()}", extra=stage)

    def _foldseek_alignment(self):
        """
        Execute foldseek sub-commands via bash interfacing.

        Triggers `foldseek easy-search`, pushing input query targets directly against formatted targets
        using `self._settings.alignment_path` to store raw tabular matches.

        Returns:
            None
        """
        stage = {"stage": "Foldseek Alignment"}
        logging.info("Running Foldseek alignment...", extra=stage)

        # Setting up paths for foldseek databases
        self._query_db_path = os.path.join(self._settings.query_dir, "query_db")
        query_db_cmd = [
            "foldseek",
            "createdb",
            self._settings.query_dir,
            self._query_db_path,
        ]
        logging.debug(
            f"Running Foldseek createdb for query with command: {' '.join([str(x) for x in query_db_cmd])}", extra=stage
        )
        subprocess.run(query_db_cmd, check=True)

        if self.fsdb_target:
            self._target_db_path = self._target_df.loc[0, "struct_path"]
            logging.debug(f"Targeting bundled human_domains Foldseek DB at {self._target_db_path}", extra=stage)
        else:
            self._target_db_path = os.path.join(self._settings.target_dir, "target_db")
            target_db_cmd = [
                "foldseek",
                "createdb",
                self._settings.target_dir,
                self._target_db_path,
            ]
            logging.debug(
                f"Running Foldseek createdb for target with command: {' '.join([str(x) for x in target_db_cmd])}",
                extra=stage,
            )
            subprocess.run(target_db_cmd, check=True)

        query_target_align_cmd = [
            "foldseek",
            "easy-search",
            self._query_db_path,
            self._target_db_path,
            self._settings.alignment_path,
            self._settings.foldseek_tmp_dir,
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
            "-v",  # verbosity
            str(
                min(3, self._settings.verbosity)
            ),  # cap foldseek verbosity at 3 (info level) since it can be very verbose at higher levels and we already have our own logging verbosity control
        ]
        logging.debug(
            f"Running Foldseek with command: {' '.join([str(x) for x in query_target_align_cmd])}", extra=stage
        )
        subprocess.run(query_target_align_cmd, check=True)
        logging.debug("Foldseek alignment completed successfully", extra=stage)

    def _local_alignment(self):
        """
        Execute traditional sequence-level alignment across inputs.

        Uses the `SequenceAligner` to locally match structured frames, resolving pairwise dependencies
        within internal dataframe representations. Outputs directly to the configuration TSV alignment file.

        Returns:
            None
        """
        stage = {"stage": "Local Alignment"}
        logging.info("Running local sequence alignments...", extra=stage)

        # Run the sequence aligner on the same preprocessed structures as foldseek uses
        aligner = SequenceAligner()
        query_records = (
            self._query_df.query("success").drop_duplicates(subset="preprocess_name").to_dict(orient="records")
        )
        target_records = (
            self._target_df.query("success").drop_duplicates(subset="preprocess_name").to_dict(orient="records")
        )
        alignment = aligner.align_records(
            query_records,
            target_records,
        )
        alignment.to_csv(self._settings.alignment_path, index=False, sep="\t")

    def _get_pockets(self):
        """
        Aggregate pocket coordinate arrays based on configured mapping logic (PISA, Passthrough, VDW, Whole Chain).

        Executes targeted extraction requests across the configured pocket methodologies. Combines all derived
        pocket residues/points into a standard composite mapping object.

        Returns:
            dict: An aggregated collection of pocket records.
        """
        stage = {"stage": "Getting Pockets"}
        logging.info("Starting pocket retrieval...", extra=stage)

        # Turns PDB Foldseek-database hits into ordinary pisa target records, so the retrieval below
        # picks them up like any other pisa entry. No-op for every other kind of target.
        self._expand_fsdb_pdb_targets()

        pisa_pockets = self._retrieve_pisa_pockets()
        passthrough_pockets = self._retrieve_passthrough_pockets()
        vdw_pockets = self._retrieve_vdw_pockets()
        whole_chain_pockets = self._retrieve_whole_chain_pockets()

        pockets = pisa_pockets | passthrough_pockets | vdw_pockets | whole_chain_pockets
        logging.debug(f"Combined pockets: {pockets}", extra=stage)
        return pockets

    def _expand_fsdb_pdb_targets(self):
        """
        Turn the hits of a PDB Foldseek-database search into ordinary PISA target records.

        A Foldseek-database target has no per-chain records of its own, so `compare_pockets` normally
        synthesises a whole-chain pseudo-pocket for each hit and leaves every `pocket_2_*` column empty.
        The PDB database is built from real PDB entries, though, so its hits have real PISA interfaces:
        this reads the hit names out of the alignment table, resolves each to a PDB ID and chain, asks
        PISA which chains that chain touches, and appends one `pisa` record per interface to
        `self._target_df`. `_retrieve_pisa_pockets` then handles them like any other pisa entry.

        The generated records carry the Foldseek entry name as their `preprocess_name` rather than the
        one `QTProcessor` derives, because that is the key alignments are stored under -- it is what
        joins these pockets back to their alignment rows and to their Foldseek transforms.

        Runs after `_alignment`, so `alignment.tsv` exists. Does nothing unless the target is a Foldseek
        database whose entries are named in the PDB style; a database of anything else (human_domains)
        keeps the synthesised whole-chain pockets.

        Returns:
            None: appends to `self._target_df` and sets `self._fsdb_pdb_target`.
        """
        if not self.fsdb_target:
            return

        stage = {"stage": "Expanding Foldseek DB Targets"}

        alignment_df = pd.read_csv(self._settings.alignment_path, sep="\t", engine="c")
        hits = {}  # foldseek entry name -> (pdb_id, chain_id)
        for hit_name in alignment_df["target"].unique().tolist():
            resolved = parse_foldseek_pdb_entry_name(hit_name)
            if resolved is not None:
                hits[hit_name] = resolved
        if not hits:
            logging.info(
                "Foldseek database target is not a PDB database; keeping whole-chain target pockets",
                extra=stage,
            )
            return
        self._fsdb_pdb_target = True

        pdb_list = sorted({pdb_id for pdb_id, _ in hits.values()})
        logging.info(
            f"Retrieving PISA interfaces for {len(pdb_list)} PDB entries behind {len(hits)} Foldseek hits",
            extra=stage,
        )

        # Same directories _retrieve_pisa_pockets uses, so the two share one cache and its own call is a
        # no-op for everything downloaded here.
        pisa_response_dir = os.path.join(self._settings.pocket_dir, "pisa_responses")
        interface_dir = os.path.join(pisa_response_dir, "interfaces")
        downloader = PisaDownloader()
        downloader.get_interfaces(
            pdb_list=pdb_list,
            summary_dir=os.path.join(pisa_response_dir, "summaries"),
            asm_dir=os.path.join(pisa_response_dir, "assemblies"),
            interface_dir=interface_dir,
        )

        # Building one record per interface the hit chain takes part in
        parser = PisaParser()
        qtprocessor = QTProcessor(
            structure_dir=self._settings.structure_dir,
            foldseek_preprocessed_structure_dir=self._settings.foldseek_preprocessed_structure_dir,
            fsdb_dir=self._settings.fsdb_dir,
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
            logging.warning("No PISA interfaces found for any Foldseek hit", extra=stage)
            return

        # Fetching structures last, so only entries that actually produced a pocket are downloaded.
        # Hits whose structure can't be fetched come back marked success=False and are dropped here;
        # _fetch_missing_structures raises only if not one of them could be fetched, which would leave
        # nothing to compare against at all.
        target_df = self._fetch_missing_structures("foldseek hit", pd.DataFrame(records), self._settings.structure_dir)
        target_df = target_df.query("success")

        logging.info(
            f"Added {len(target_df)} PISA target pockets from {target_df['preprocess_name'].nunique()} Foldseek hits",
            extra=stage,
        )
        # Row 0 stays the database record itself -- _foldseek_alignment and _align_structs read its struct_path.
        self._target_df = pd.concat([self._target_df, target_df], ignore_index=True)

    def _retrieve_pisa_pockets(self):
        """
        Request, parse, and translate remote pocket mapping endpoints through the PDBe PISA service.

        Downloads corresponding assembly and interface data over REST points. Extracts residue maps
        associating directly to structured indices previously resolved, persisting local copies to the cache.

        Returns:
            dict: Translated pocket coordinates indexed by `pocket_id`.
        """
        stage = {"stage": "Retrieving PISA Pockets"}
        logging.info("Checking for PISA pockets...", extra=stage)

        # Selecting relevant records from quert and taget dataframes
        pisa_df = (
            pd.concat([self._query_df, self._target_df], ignore_index=True)
            .query("success and pocket_method == 'pisa'")
            .drop_duplicates(subset=["struct_info", "chain_info"])
        )
        if len(pisa_df) == 0:
            logging.info("No PISA pockets to retrieve", extra=stage)
            return {}
        else:
            logging.info(f"{len(pisa_df)} PISA pockets to retrieve", extra=stage)

        pisa_response_dir = os.path.join(self._settings.pocket_dir, "pisa_responses")
        pisa_pdb_list = pisa_df["struct_info"].unique().tolist()
        logging.debug(f"PDBs for which to retrieve PISA pockets: {pisa_pdb_list}", extra=self._log_extra)
        downloader = PisaDownloader()
        downloader.get_interfaces(
            pdb_list=pisa_pdb_list,
            summary_dir=os.path.join(pisa_response_dir, "summaries"),
            asm_dir=os.path.join(pisa_response_dir, "assemblies"),
            interface_dir=os.path.join(pisa_response_dir, "interfaces"),
        )

        parser = PisaParser()
        pisa_pockets = parser.get_pockets_from_records(
            records=pisa_df.to_dict(orient="records"),
            in_dir=os.path.join(pisa_response_dir, "interfaces"),
        )
        logging.debug(f"Extracted PISA pockets: {pisa_pockets}", extra=self._log_extra)

        for _, row in pisa_df.iterrows():
            if row["pocket_id"] in pisa_pockets:
                pisa_pockets[row["pocket_id"]] = parse_pocket_from_struct(
                    struct=row["struct_path"],
                    chain_id=row["chain_info"].split("_")[0],
                    pocket_residues=[int(x) for x in pisa_pockets[row["pocket_id"]]["res_auth_ids"]],
                    pocket=pisa_pockets[row["pocket_id"]],
                )
        logging.debug(f"Extracted PISA pockets with coords: {pisa_pockets}", extra=self._log_extra)

        with open(os.path.join(self._settings.pocket_dir, "pisa_pockets.json"), "w") as f:
            json.dump(pisa_pockets, f)

        return pisa_pockets

    def _retrieve_passthrough_pockets(self):
        """
        Convert manually defined list-based user parameters directly into atomic coordinate graphs.

        Target data frames assigned a 'passthrough' target type have explicitly identified residues parsed
        statically based on existing indices within pre-existing mmcif structures.

        Returns:
            dict: Translated pocket coordinates indexed by `pocket_id`.
        """
        stage = {"stage": "Passthrough Pocket Calculation"}
        logging.info("Checking for passthrough pockets...", extra=stage)

        pt_df = pd.concat([self._query_df, self._target_df], ignore_index=True).query(
            "success and pocket_method == 'passthrough'"
        )
        if len(pt_df) == 0:
            logging.info("No passthrough pockets to retrieve", extra=stage)
            return {}
        else:
            logging.info(f"{len(pt_df)} passthrough pockets to retrieve", extra=stage)

        # for each pocket in query and target parse pocket info from the structure and store in a dict
        passthrough_pockets = {}
        for _, row in pt_df.iterrows():
            pocket_residues = [int(x) for x in row["residue_info"].split(",")]
            pocket = parse_pocket_from_struct(
                struct=row["struct_path"],
                chain_id=row["chain_info"].split("_")[0],
                pocket_residues=pocket_residues,
            )
            passthrough_pockets[row["pocket_id"]] = pocket
        with open(os.path.join(self._settings.pocket_dir, "passthrough_pockets.json"), "w") as f:
            json.dump(passthrough_pockets, f, indent=4)
        logging.debug(f"Extracted passthrough pockets: {passthrough_pockets}", extra=self._log_extra)

        return passthrough_pockets

    def _retrieve_whole_chain_pockets(self):
        """
        Build the "pocket" for an open search: every CA-bearing residue of the chain.

        An entry that names a structure but no pocket ("4Q5J:B", or "4Q5J" for the default chain) asks
        whether the query pocket resembles anything on that chain at all. Unlike the whole-chain pseudo-pocket
        `compare_pockets` synthesises for Foldseek-database hits, this is a real pocket dict parsed from the
        structure, so it carries residue codes and CA coordinates and can be superposed.

        Returns:
            dict: Translated pocket coordinates indexed by `pocket_id`.
        """
        stage = {"stage": "Whole Chain Pocket Calculation"}
        logging.info("Checking for whole chain pockets...", extra=stage)

        wc_df = pd.concat([self._query_df, self._target_df], ignore_index=True).query(
            "success and pocket_method == 'whole_chain'"
        )
        if len(wc_df) == 0:
            logging.info("No whole chain pockets to retrieve", extra=stage)
            return {}
        else:
            logging.info(f"{len(wc_df)} whole chain pockets to retrieve", extra=stage)

        whole_chain_pockets = {}
        for _, row in wc_df.iterrows():
            pocket = parse_pocket_from_struct(
                struct=row["struct_path"],
                chain_id=row["chain_info"].split("_")[0],
                pocket_residues=None,  # None means the whole chain
            )
            # A missing structure or chain gives None back. Storing that would fail later with an opaque
            # TypeError inside compare_pockets, so drop the entry and say which one it was -- an unreadable
            # chain is much more likely here, where the chain can come from the default rather than the user.
            if pocket is None:
                logging.warning(
                    f"Could not parse chain {row['chain_info']} of {row['struct_info']} for {row['pocket_id']}, "
                    "skipping this entry",
                    extra=stage,
                )
                continue
            whole_chain_pockets[row["pocket_id"]] = pocket
        with open(os.path.join(self._settings.pocket_dir, "whole_chain_pockets.json"), "w") as f:
            json.dump(whole_chain_pockets, f, indent=4)
        logging.debug(f"Extracted whole chain pockets: {whole_chain_pockets}", extra=self._log_extra)

        return whole_chain_pockets

    def _retrieve_vdw_pockets(self):
        """
        Evaluate pocket clusters structurally using Van-der-Waals (VDW) interaction overlapping metrics.

        Passes valid structures referencing 'vdw' methods into `PocketCalculator` implementations to
        synthesize interactive residue domains from a protein-motif complex.

        Returns:
            dict: Translated pocket coordinates indexed by `pocket_id`.
        """
        stage = {"stage": "VDW Pocket Calculation"}
        logging.info("Checking for VDW pockets...", extra=stage)
        vdw_df = pd.concat([self._query_df, self._target_df], ignore_index=True).query(
            "success and pocket_method == 'vdw'"
        )
        if len(vdw_df) == 0:
            logging.info("No VDW pockets to retrieve", extra=stage)
            return {}
        else:
            logging.info(f"{len(vdw_df)} VDW pockets to retrieve", extra=stage)

        vdw_pockets = {}
        pc = PocketCalculator()
        for _, row in vdw_df.iterrows():

            pocket = pc.pocket_overlap(
                structure=row["struct_path"],
                domain_chain=row["chain_info"].split("_")[0],
                motif_chain=row["chain_info"].split("_")[1],
            )
            vdw_pockets[row["pocket_id"]] = pocket
        with open(os.path.join(self._settings.pocket_dir, "vdw_pockets.json"), "w") as f:
            json.dump(vdw_pockets, f, indent=4)
        logging.debug(f"Extracted VDW pockets: {vdw_pockets}", extra=self._log_extra)

        return vdw_pockets

    def _compare_pockets_based_on_alignment(self, pockets):
        """
        Merge sequence or fold structure alignments traversing spatial pocket metrics internally to evaluate correlations.

        Processes raw mapping dictionaries tracking unresolved structural aliases. Utilizes a referenced BLOSUM
        scorecard for computing structural compatibility values output into a final result dataframe file.

        Args:
            pockets (dict): Populated pool of pocket residue locations collected via initial configurations.

        Returns:
            None
        """
        stage = {"stage": "Comparing Pockets Based on Alignment"}

        logging.info("Reading alignment results...", extra=stage)
        alignment_df = pd.read_csv(self._settings.alignment_path, sep="\t", engine="c")
        blosum_path = os.path.join(os.path.dirname(__file__), "blosum62.bla")

        logging.info(f"{len(alignment_df)} alignment pairs to compare", extra=stage)
        logging.debug(f"Alignment pairs: \n{alignment_df.head()}", extra=stage)

        preproc_to_ids = {}
        for _, row in pd.concat([self._query_df, self._target_df], ignore_index=True).iterrows():
            # Keyed by preprocess_name: one chain can carry several pockets (several pocket_ids), so
            # membership must be tested on the key, not on the pocket_id, or each chain keeps only its last.
            if row["preprocess_name"] in preproc_to_ids:
                if row["pocket_id"] not in preproc_to_ids[row["preprocess_name"]]:
                    preproc_to_ids[row["preprocess_name"]].append(row["pocket_id"])
            else:
                preproc_to_ids[row["preprocess_name"]] = [row["pocket_id"]]
        logging.debug(f"Preprocessed name to pocket ID mapping: {preproc_to_ids}", extra=stage)

        pockets_df, unknown_alias, incorrect_mapping = compare_pockets(
            alignment_df,
            pockets,
            preproc_to_ids=preproc_to_ids,
            blosum_path=blosum_path,
            # A PDB Foldseek database has real PISA pockets for its hits (see _expand_fsdb_pdb_targets);
            # every other database has no target records at all, so its pockets must be synthesised.
            synthesise_target_pockets=self.fsdb_target and not self._fsdb_pdb_target,
        )

        # Logging cases where a residue was given a single cahr name unfamiliar to pocketmapper
        if len(unknown_alias) > 0:
            unknown_alias_path = os.path.join(self._settings.results_dir, "unknown_ids.json")
            logging.warning("Unknown Foldseek Alias, see unknown_alias.json in results directory", extra=stage)
            with open(unknown_alias_path, "w") as f:
                json.dump(jsonify_dict(dict(unknown_alias)), f)

        # logging cases where foldseek mapping had low sequence identity to the parsed structure
        if len(incorrect_mapping) > 0:
            incorrect_mapping_path = os.path.join(self._settings.results_dir, "incorrect_mapping.json")
            logging.warning("Foldseek mapping with low sequence identity to parsed structure", extra=stage)
            with open(incorrect_mapping_path, "w") as f:
                json.dump(jsonify_dict(dict(incorrect_mapping)), f)

        # Writing pocket comparison results to output file
        output_path = self._settings.pocket_comparison_path
        pockets_df.to_csv(output_path, index=False, sep="\t")
        logging.info(f"Pocket comparison results saved to {output_path}", extra=stage)

    def _align_structs(self):
        """
        Perform structural superposition of target structures against the query reference frame.

        Takes the top N targets (as defined by `align_count`) from the pocket comparison results
        and performs a structural alignment using the `StructureAligner` class.
        The aligned structures will be saved to the target directory for downstream analysis.

        Returns:
            None
        """
        stage = {"stage": "Structural Alignment"}
        if self._settings.align_count <= 0:
            logging.info("No Aligned Structures to Process", extra=stage)
            return

        method = self._settings.align_struct_method  # already "pocket" or "foldseek"
        logging.info(f"Performing structural alignment of target structures on the {method}...", extra=stage)

        # Pre-loading
        aligner = StructureAligner()
        pocket_comparison_df = pd.read_csv(self._settings.pocket_comparison_path, sep="\t", engine="c")
        alignment_df = None
        pocket_transform_df = None
        if method == "foldseek":
            alignment_df = pd.read_csv(
                self._settings.alignment_path,
                sep="\t",
                engine="c",
                index_col=["query", "target"],
            )
        else:
            # (pocket_1, pocket_2) is unique -- compare_pockets' existing_calcs scores each pair once.
            pocket_transform_df = pocket_comparison_df.dropna(subset=["p2_to_p1_u", "p2_to_p1_t"]).set_index(
                ["pocket_1", "pocket_2"]
            )[["p2_to_p1_u", "p2_to_p1_t"]]

        # For each query structure, align the top N target structures
        qt_id_map = {}
        unique_target_ids = set()
        for record in self._query_df.to_dict(orient="records"):
            query_id = record["pocket_id"]
            logging.debug(f"Processing query {query_id} for structural alignment", extra=stage)

            # Select the top N target structures based on pocket comparison metrics. Targets sharing no
            # pocket residues with the query are excluded: there is no common set of residues to superpose
            # on, and their overlap metrics are empty so they would sort arbitrarily.
            #
            # A whole-chain target -- an open search, or a Foldseek-DB hit -- has no jaccard_index, so it
            # sorts to the end and is ranked by the secondary key, min_overlap_similarity, instead.
            candidates = pocket_comparison_df.query(f"pocket_1 == '{query_id}' and overlap_count > 0")
            overlapping_count = len(candidates)
            if method == "pocket":
                # _superpose fits nothing below three overlapping residues, so those targets have no
                # transform. Drop them here rather than when writing, or they would eat align_count
                # slots and the run would quietly produce fewer structures than asked for.
                candidates = candidates.dropna(subset=["p2_to_p1_u", "p2_to_p1_t"])

            target_ids = (
                candidates.sort_values(by=["jaccard_index", "min_overlap_similarity"], ascending=False)
                .head(self._settings.align_count)
                .loc[:, "pocket_2"]
                .to_list()
            )
            if not target_ids:
                if overlapping_count:
                    logging.info(
                        f"No target overlaps the pocket of query {query_id} by the three residues a "
                        "superposition needs; skipping its structural alignment",
                        extra=stage,
                    )
                else:
                    logging.info(
                        f"No target overlaps the pocket of query {query_id}; skipping its structural alignment",
                        extra=stage,
                    )
                continue
            logging.debug(f"Top target IDs for query {query_id}: {target_ids}", extra=stage)
            qt_id_map[query_id] = target_ids
            unique_target_ids.update(target_ids)

        if not unique_target_ids:
            logging.info("No query/target pair shares pocket residues, nothing to superpose", extra=stage)
            return

        self._query_df = self._query_df.set_index("pocket_id")
        if self.fsdb_target is False:
            target_record_df = self._target_df.set_index("pocket_id")
        else:  # If the target is a Foldseek database we need to make pdb structures from required entries
            source_db_path = self._target_df.loc[0, "struct_path"]
            logging.debug(f"Using Foldseek database at {source_db_path} for structural alignment", extra=stage)

            # With a PDB database the target IDs are pocket IDs ("4Q5J:B_F"), not database entry names, so
            # map them back through the records built by _expand_fsdb_pdb_targets. One pocket ID can come
            # from more than one entry (the same chain in two assemblies) -- keep the first, or the lookup
            # below returns duplicate rows and the structure gets superposed twice.
            if self._fsdb_pdb_target:
                id_to_entry = (
                    self._target_df.dropna(subset=["preprocess_name"])
                    .drop_duplicates(subset="pocket_id", keep="first")
                    .set_index("pocket_id")["preprocess_name"]
                )
                target_entry_names = {target_id: id_to_entry[target_id] for target_id in unique_target_ids}
            else:
                target_entry_names = {target_id: target_id for target_id in unique_target_ids}

            # Get chain IDs corresponding to the required entries from the Foldseek database lookup file
            source_db_lookup_path = source_db_path + ".lookup"
            source_db_lookup_df = pd.read_csv(
                source_db_lookup_path, sep="\t", header=None, names=["chain_id", "name", "struct_id"]
            )
            source_db_lookup_df = source_db_lookup_df.set_index("name")
            chain_ids = source_db_lookup_df.loc[list(target_entry_names.values()), "chain_id"].tolist()

            # Make directory for subdb
            subdb_dir = os.path.join(self._settings.aligned_structure_dir, "fsdb")
            os.makedirs(subdb_dir, exist_ok=True)

            # Create a file listing the required chain IDs for the subdb creation
            subdb_chain_id_path = os.path.join(subdb_dir, "required_chain_ids.txt")
            with open(subdb_chain_id_path, "w") as f:
                for target_id in chain_ids:
                    f.write(f"{target_id}\n")

            # Create the subdb using foldseek's createsubdb command
            subdb_path = os.path.join(subdb_dir, "subdb")
            subdb_command = [
                "foldseek",
                "createsubdb",
                subdb_chain_id_path,
                source_db_path,
                subdb_path,
            ]
            subprocess.run(subdb_command, check=True)

            # Create a directory for extracted structures
            subdb_struct_dir = os.path.join(self._settings.aligned_structure_dir, "fsdb_structures")
            os.makedirs(subdb_struct_dir, exist_ok=True)

            # Convert the subdb to PDB format using foldseek's convert2pdb command
            convert2pdb_command = [
                "foldseek",
                "convert2pdb",
                "--pdb-output-mode",
                "1",
                subdb_path,
                subdb_struct_dir,
            ]
            subprocess.run(convert2pdb_command, check=True)

            # Make record df for the target records based on the unique target IDs and the subdb structure
            # directory. chain_info stays None: each extracted structure holds exactly the one chain of its
            # database entry, which foldseek_transform takes as the domain chain.
            target_record_df = pd.DataFrame(
                {"pocket_id": list(target_entry_names.keys()), "preprocess_name": list(target_entry_names.values())}
            )
            target_record_df["chain_info"] = None
            target_record_df["struct_path"] = target_record_df["preprocess_name"].apply(
                lambda x: os.path.join(subdb_struct_dir, f"{x}.pdb")
            )
            target_record_df = target_record_df.set_index("pocket_id")

        for query_id, target_ids in qt_id_map.items():
            # pocket_id is the index of _query_df, so it is not in the row dict -- put it back, since
            # foldseek_transform reads it for the COMPND metadata.
            query_record = self._query_df.loc[query_id].to_dict()
            query_record["pocket_id"] = query_id
            logging.debug(f"Query record for '{query_id}': {json.dumps(query_record, indent=4)}", extra=stage)
            # Fetch the corresponding target records. A pocket_2 need not be a target: when a query and a
            # target share a chain they share a preprocess_name, so compare_pockets pairs every pocket on
            # that chain with every other and some rows come back with a query-only pocket_id in pocket_2.
            # Those have no target structure to superpose, so drop them rather than let .loc raise.
            known_target_ids = [target_id for target_id in target_ids if target_id in target_record_df.index]
            missing_target_ids = [target_id for target_id in target_ids if target_id not in target_record_df.index]
            if missing_target_ids:
                logging.debug(
                    f"Skipping non-target pocket(s) {missing_target_ids} when superposing onto '{query_id}'",
                    extra=stage,
                )
            top_target_records = target_record_df.loc[known_target_ids].reset_index().to_dict(orient="records")
            logging.debug(
                f"Top target records for query '{query_id}': {json.dumps(top_target_records, indent=4)}", extra=stage
            )

            # The query is the reference frame every target is superposed onto, so it must lead the list.
            aln_records = [query_record] + top_target_records
            if len(aln_records) > 1:
                out_path = os.path.join(self._settings.aligned_structure_dir, f"{safe_filename(query_id)}.pdb")
                if method == "foldseek":
                    aligner.foldseek_transform(
                        aln_records=aln_records,
                        alignment_df=alignment_df,
                        out_path=out_path,
                    )
                else:
                    # transforms is positional, not keyed by pocket_id: a self-comparison gives the
                    # query and a target the same pocket_id, so a dict would collide.
                    transforms = [None]  # the query is the reference frame, placed untransformed
                    for record in top_target_records:
                        try:
                            row = pocket_transform_df.loc[(query_id, record["pocket_id"])]
                        except KeyError:
                            transforms.append(None)
                            logging.warning(
                                f"No pocket superposition for {query_id} against {record['pocket_id']}",
                                extra=stage,
                            )
                            continue
                        transforms.append(parse_pocket_transform(row["p2_to_p1_u"], row["p2_to_p1_t"]))
                    aligner.transform(aln_records=aln_records, transforms=transforms, out_path=out_path)

    def _delete_tmp(self):
        """
        Cleanup temporary cache directories holding extracted domains following completion cycles.

        Erases dynamically built sub-directories holding processed query logic, removing intermediate
        uncompressed sequences unless flagged for extended review formats via `human_domains`.

        Returns:
            None
        """
        tmp_dirs = [
            "query_dir",
            "target_dir",
        ]
        if self._settings.foldseek:
            tmp_dirs.append("foldseek_tmp_dir")

        # TODO this is unsafe
        for dir in tmp_dirs:
            shutil.rmtree(getattr(self._settings, dir))


def main():
    try:
        fire.Fire(PocketMapper())
    except PocketMapperError:
        # Already logged with full stage context at the raise site.
        sys.exit(1)


if __name__ == "__main__":
    main()
