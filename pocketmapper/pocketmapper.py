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
from pocketmapper.lib import jsonify_dict, safe_filename
from pocketmapper.exceptions import PocketMapperError
from pocketmapper.pisa_downloader import PisaDownloader
from pocketmapper.pisa_parser import PisaParser
from pocketmapper.pocket_comparison import compare_pockets
from pocketmapper.sequence_aligner import SequenceAligner
from pocketmapper.pocket_calculator import PocketCalculator
from pocketmapper.qt_processor import QTProcessor
from pocketmapper.structure_aligner import StructureAligner
from pocketmapper.structure_fetcher import StructureFetcher
from pocketmapper.structure_preprocessor import StructurePreprocessor
from pocketmapper.lib_struct import parse_pocket_from_struct
from pocketmapper.constants import HELP_MESSAGE


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
    foldseek: bool = False
    align_count: int = 10
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
            foldseek (bool, optional): Use foldseek for structure alignment instead of local sequence alignment.
            align_count (int, optional): Number of top targets to superpose onto each query.
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
            else:
                logging.critical(
                    "Foldseek database specified as target but foldseek is not enabled. Please set --foldseek True.",
                    extra=self._log_extra,
                )
                raise PocketMapperError(
                    "Foldseek database specified as target but foldseek is not enabled. Please set --foldseek True."
                )
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
            "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,lddt,qaln,taln,u,t,qseq,tseq",
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
        Aggregate pocket coordinate arrays based on configured mapping logic (PISA, Passthrough, VDW).

        Executes targeted extraction requests across the configured pocket methodologies. Combines all derived
        pocket residues/points into a standard composite mapping object.

        Returns:
            dict: An aggregated collection of pocket records.
        """
        stage = {"stage": "Getting Pockets"}
        logging.info("Starting pocket retrieval...", extra=stage)

        pisa_pockets = self._retrieve_pisa_pockets()
        passthrough_pockets = self._retrieve_passthrough_pockets()
        vdw_pockets = self._retrieve_vdw_pockets()

        pockets = pisa_pockets | passthrough_pockets | vdw_pockets
        logging.debug(f"Combined pockets: {pockets}", extra=stage)
        return pockets

    # TODO If target is PDB foldseek database, a set of target queries needs to be generated to be used as pockets
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

        # TODO PDB: get list of pdbs from alignment.tsv and extend the query list
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

        # TODO PDB: need a parser that can handle just 1 chain ID
        parser = PisaParser()
        pisa_pockets = parser.get_pockets_from_records(
            records=pisa_df.to_dict(orient="records"),
            in_dir=os.path.join(pisa_response_dir, "interfaces"),
        )
        logging.debug(f"Extracted PISA pockets: {pisa_pockets}", extra=self._log_extra)

        # TODO PDB: Should have all the info neeed to rerun earlier parts of the pipeline and set up the target df
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
            alignment_df, pockets, preproc_to_ids=preproc_to_ids, blosum_path=blosum_path, alphafold=self.fsdb_target
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
        else:
            logging.info("Performing structural alignment of target structures...", extra=stage)

        # Pre-loading
        aligner = StructureAligner()
        alignment_df = pd.read_csv(
            self._settings.alignment_path,
            sep="\t",
            engine="c",
            index_col=["query", "target"],
        )
        pocket_comparison_df = pd.read_csv(self._settings.pocket_comparison_path, sep="\t", engine="c")

        # For each query structure, align the top N target structures
        qt_id_map = {}
        unique_target_ids = set()
        for record in self._query_df.to_dict(orient="records"):
            query_id = record["pocket_id"]
            logging.debug(f"Processing query {query_id} for structural alignment", extra=stage)

            # Select the top N target structures based on pocket comparison metrics. Targets sharing no
            # pocket residues with the query are excluded: there is no common set of residues to superpose
            # on, and their overlap metrics are empty so they would sort arbitrarily.
            target_ids = (
                pocket_comparison_df.query(f"pocket_1 == '{query_id}' and overlap_count > 0")
                .sort_values(by=["pocket_1_pct_overlap", "min_overlap_similarity"], ascending=False)
                .head(self._settings.align_count)
                .loc[:, "pocket_2"]
                .to_list()
            )
            if not target_ids:
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

            # Get chain IDs corresponding to the unique target IDs from the Foldseek database lookup file
            source_db_lookup_path = source_db_path + ".lookup"
            source_db_lookup_df = pd.read_csv(
                source_db_lookup_path, sep="\t", header=None, names=["chain_id", "name", "struct_id"]
            )
            source_db_lookup_df = source_db_lookup_df.set_index("name")
            chain_ids = source_db_lookup_df.loc[list(unique_target_ids), "chain_id"].tolist()

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

            # Make record df for the target records based on the unique target IDs and the subdb structure directory
            target_record_df = pd.DataFrame({"preprocess_name": list(unique_target_ids)})
            target_record_df["chain_info"] = None
            target_record_df["pocket_id"] = target_record_df["preprocess_name"]
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
            # Fetch the corresponding target records
            top_target_records = target_record_df.loc[target_ids].reset_index().to_dict(orient="records")
            logging.debug(
                f"Top target records for query '{query_id}': {json.dumps(top_target_records, indent=4)}", extra=stage
            )

            # The query is the reference frame every target is superposed onto, so it must lead the list.
            aln_records = [query_record] + top_target_records
            if len(aln_records) > 1:
                aligner.foldseek_transform(
                    aln_records=aln_records,
                    alignment_df=alignment_df,
                    out_path=os.path.join(self._settings.aligned_structure_dir, f"{safe_filename(query_id)}.pdb"),
                )

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
