"""
PocketMapper: A tool for mapping and analyzing protein pockets.

Author: Lachlan Ellingboe

"""

from importlib.resources import files
import fire
import logging
import logging.config
import json
import subprocess
import pandas as pd
import os
from datetime import datetime
import shutil
from pocketmapper import lib
from pocketmapper.pisa_downloader import PisaDownloader
from pocketmapper.pisa_parser import PisaParser
from pocketmapper.sequence_aligner import SequenceAligner
from pocketmapper.pocket_calculator import PocketCalculator
from pocketmapper.qt_processor import QTProcessor
from pocketmapper.structure_aligner import StructureAligner
from pocketmapper.structure_fetcher import StructureFetcher
from pocketmapper.structure_preprocessor import StructurePreprocessor
from pocketmapper.lib_struct import parse_pocket_from_struct
from pocketmapper import human_domains


class PocketMapper:
    def __init__(self):
        """
        Initialize the PocketMapper instance.

        Sets up the default logging configuration, formatting strings, and points to
        bundled structural databases if applicable.
        """
        self._log_extra = {"stage": "init"}
        self._human_domains_db_path = files(human_domains).joinpath("human_v3_20260531")
        self._log_fmt = "%(levelname)s: %(stage)s - %(msg)s"
        logging.getLogger(__name__)
        logging.basicConfig(level=logging.CRITICAL, format=self._log_fmt)

    def _configure_logging(self, settings):
        """
        Configure logging level and handlers based on user settings.

        Sets both a console stream handler and a file handler (to `log_path` in `settings`),
        adjusting the verbosity depending on the user's input.

        Args:
            settings (dict): Configuration dictionary containing keys:
                - "verbosity" (int): The verbosity level (4=DEBUG, 3=INFO, 2=WARNING, else ERROR).
                - "log_path" (str): Path to write the info.log file.
        """
        self._log_extra.update({"stage": "Configuring Logging"})

        # Set log level based on verbosity setting (default to INFO if not set)
        log_level = None
        if settings["verbosity"] == 4:
            log_level = "DEBUG"
        elif settings["verbosity"] == 3:
            log_level = "INFO"
        elif settings["verbosity"] == 2:
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
                    "filename": settings["log_path"],
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
        align_struct=None,
        align_count=None,
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
            align_struct (bool, optional): Align target structures after pocket comparison.
            align_count (int, optional): Number of top alignments to consider for pocket comparison.
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
        self._align_struct = align_struct
        self._align_count = align_count

        self._check_help_search()  # Checks if help flag is set and if so prints the help message and exits
        self._configure_workflow()  # configures the settings which have already been read
        self._configure_query_target()  # parses the query and target inputs to determine their types and sets up the relevant data structures for each entry
        self._fetch_missing_structures()  # Fetch any missing structures
        self._alignment()  # Align the query and target structures using either local sequence alignment or foldseek based on the settings
        pockets = self._get_pockets()  # Adds seq_pos and ca-coords to the pocket info dict

        # TODO Remove this hack after preserving the method
        if False:  # hack to make ATP pocket search working
            pc = PocketCalculator()
            pockets = pc.atp_pocket_overlap(
                r"/Users/lellingboe/Work/data/kinase_edit/atp_pocket/pocketmapper_cache/divided_structs/3BU5_A_B.cif.gz",
                "A",
                "3BU5_A_B",
            )
            with open(
                r"/Users/lellingboe/Work/data/kinase_edit/atp_pocket/pocketmapper_cache/pockets/3BU5_A_B_atp_pocket.json",
                "w",
            ) as f:
                json.dump(pockets, f, indent=4)

        self._compare_pockets_based_on_alignment(pockets)
        self._align_structures()
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
            help_msg = """
        PocketMapper - A tool for mapping and analyzing protein pockets.

        Usage:
            pocketmapper search [OPTIONS]

        Primary options:
            --query QUERY            Query identifier or path. Accepts:
                        - 'PDB_CHAIN_CHAIN' (e.g., 1ABC_A_B)
                        - path to a file listing PDB_CHAIN_CHAIN entries (each line)
            --target TARGET          Target identifier or path. Accepts:
                        - 'PDB_CHAIN_CHAIN' (e.g., 2XYZ_C_D)
                        - path to a file listing PDB_CHAIN_CHAIN entries (each line)
                        - special foldseek DB alias 'human_domains' to use the bundled Foldseek DB
            --settings FILE          Path to a JSON file of {"ARG": "VALUE", ...} (overridden by explicit CLI args)
            --cache_dir DIR          Directory for caching files
            --results_dir DIR        Directory for writing results
            --verbosity LEVEL        Set verbosity level (4=DEBUG, 3=INFO, 2=WARNING, else ERROR)
            --foldseek BOOL          Whether to use foldseek for structure alignment instead of local sequence alignment
            --align_struct BOOL      Whether to align target structures after pocket comparison
            --help                   Show this help message and exit

        Advanced Options (set via settings JSON):
            structure_dir            Directory to store downloaded/available structures
            pocket_dir               Directory to store calculated pockets
            foldseek_preprocessed_struct_dir       Directory for preprocessed/divided structures
            query_dir                Temporary directory for query divided structures
            target_dir               Temporary directory for target divided structures
            alignment_path           Path to write alignment TSV
            pocket_comparison_path   Path to write pocket comparison TSV
            foldseek                 Use Foldseek for alignment (bool). If true and target == 'ted', uses bundled DB.
            pisa_pockets             Retrieve pockets via PISA (bool)
            structure                If set, treat inputs as raw structure files (bool)

        Description:
            Orchestrates fetching/preprocessing of structures, runs local or Foldseek alignments,
            fetches pockets (PISA), extracts atom coordinates from mmCIF files, compares pockets
            using alignments and scoring, and writes results to the results directory.

        Examples:
            # Single pair using local alignment and default settings
            pocketmapper search --query 1ABC_A_B --target 2XYZ_C_D --results_dir ./out

            # Batch mode using files with one PDB_CHAIN_CHAIN per line
            pocketmapper search --query queries.txt --target targets.txt --settings config.json

            # Use Foldseek (set foldseek True). When using the built-in human_domains DB:
            pocketmapper search --query 1ABC_A_B --target human_domains --foldseek True --results_dir ./out_fs

            # Override cache and set verbosity to debug
            pocketmapper search --query 1ABC_A_B --target 2XYZ_C_D --cache_dir /tmp/cache --verbosity 4

        Notes:
            - Query/target inputs are interpreted either as single PDB_CHAIN_CHAIN strings or as file paths.
            - Boolean settings can be provided on the command line (e.g., --foldseek True).
            - Use a settings JSON to persist complex configurations; CLI options override settings file values.

        For more information, see the project README or the github repository.
                """
            print(help_msg)
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
        now = datetime.now().strftime("%y%m%d_%H%M%S")
        self._settings = {
            "query": None,
            "target": None,
            "cache_dir": "pocketmapper_cache",
            "results_dir": f"pocketmapper_results_{now}",
            "query_pocket_method": None,
            "target_pocket_method": None,
            "foldseek": False,
            "align_struct": False,
            "align_count": 10,
            "verbosity": 3,  # default to info level logging
        }

        # 2. Populate settings from the settings file if provided
        if self._settings_file is not None:
            if not os.path.isfile(self._settings_file):
                logging.critical(f"Settings file not found: {self._settings_file}", extra=self._log_extra)
                exit(1)
            try:
                with open(self._settings_file) as f:
                    settings_data_from_file = json.load(f)
                    self._settings.update(settings_data_from_file)
            except Exception:
                logging.critical(
                    f"Error reading settings file: {self._settings_file}. Is it in JSON format?", extra=self._log_extra
                )
                exit(1)

        # 3. Override settings with explicit command-line arguments
        # Using a mapping of settings keys to the instance attributes
        cli_args_mapping = {
            "query": self._query,
            "target": self._target,
            "cache_dir": self._cache_dir,
            "results_dir": self._results_dir,
            "foldseek": self._foldseek,
            "verbosity": self._verbosity,
            "align_struct": self._align_struct,
            "align_count": self._align_count,
        }
        for key, value in cli_args_mapping.items():
            if value is not None:
                self._settings[key] = value

        # 4. Computed paths
        cache_dir = self._settings["cache_dir"]
        results_dir = self._settings["results_dir"]

        derived_paths = {
            # Default directories stemming from cache_dir
            "structure_dir": os.path.join(cache_dir, "ref_structures"),
            "pocket_dir": os.path.join(cache_dir, "pockets"),
            "foldseek_tmp_dir": os.path.join(cache_dir, "foldseek_tmp"),
            "foldseek_preprocessed_structure_dir": os.path.join(cache_dir, "foldseek_preprocessed_structures"),
            # Default directories and files stemming from results_dir
            "query_dir": os.path.join(results_dir, "query_structures"),
            "target_dir": os.path.join(results_dir, "target_structures"),
            "alignment_path": os.path.join(results_dir, "alignment.tsv"),
            "pocket_comparison_path": os.path.join(results_dir, "pocket_comparison.tsv"),
            "job_settings_path": os.path.join(results_dir, "job_settings.json"),
            "log_path": os.path.join(results_dir, "info.log"),
        }

        # Only overwrite paths if they were not explicitly provided via the settings file
        for key, path_val in derived_paths.items():
            if key not in self._settings:
                self._settings[key] = path_val

        # Ensure all necessary directories exist before proceeding, creating them if needed
        dirs_to_create = [
            "structure_dir",
            "query_dir",
            "target_dir",
            "pocket_dir",
            "foldseek_preprocessed_structure_dir",
        ]
        for dir_key in dirs_to_create:
            path = self._settings[dir_key]
            try:
                os.makedirs(path, exist_ok=True)
            except OSError:
                logging.critical(f"Error creating directory {path}", extra=self._log_extra)
                exit(1)

        self._configure_logging(self._settings)
        logging.debug(f"Settings: {json.dumps(self._settings, indent=4)}", extra=self._log_extra)

        # 5. Output dump
        try:
            os.makedirs(os.path.dirname(self._settings["job_settings_path"]), exist_ok=True)
            with open(self._settings["job_settings_path"], "w") as f:
                json.dump(self._settings, f, indent=4)
            logging.info(
                f"Settings successfully dumped to {self._settings['job_settings_path']}", extra=self._log_extra
            )
        except Exception as e:
            logging.error(
                f"Failed to dump settings to {self._settings['job_settings_path']}: {e}", extra=self._log_extra
            )

    def _configure_query_target(self):
        """
        Determine and process input data (formats and types) for the query and target constraints.

        Uses the `QTProcessor` class to load, parse, and validate query vs. target identities
        and requested pocket methodologies (e.g. "pisa", "passthrough"). Updates local dataframes.

        Returns:
            None: Instantiates `self._query_df` and `self._target_df`.
        """
        self._log_extra.update({"stage": "Determine Query/Target Types"})

        qtprocessor = QTProcessor(settings=self._settings)
        self._query_df, self._target_df = qtprocessor.process_qt_cmdline_input()

    def _fetch_missing_structures(self):
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

        for name, df in [("query", self._query_df), ("target", self._target_df)]:
            logging.debug(f"{name.capitalize()} data before fetching structures: \n{df.head()}", extra=self._log_extra)

            # Get list of unique structures to fetch based on struct_info and struct_type
            unique_records = df.drop_duplicates(subset="struct_info").to_dict(orient="records")

            # Update structure fetcher and fetch structures
            structure_fetcher.set_output_directory(self._settings["structure_dir"])
            structure_fetcher.update_cache()
            results = structure_fetcher.fetch_structures(unique_records)
            logging.debug(f"Structure fetcher results: {results}", extra=self._log_extra)

            # Update the dataframe with success/failure information
            df["success"] = df["struct_info"].map(results).fillna(False)
            df.loc[~df["success"], "failure_reason"] = "structure_not_found"

            # Logging results of structure fetching and updating query and target data with success/failure info
            logging.info(
                f"{sum(results.values())}/{len(results)} {name} required structures available",
                extra=self._log_extra,
            )
            if len(df.query("success == False")) > 0:
                logging.warning(
                    f"Missing structures for {name}(s): {', '.join(df.loc[~df['success'], 'pocket_id'].unique().tolist())}",
                    extra=self._log_extra,
                )

        # Verifying sufficient structures were found to continue
        exit_flag = False
        if self._query_df["success"].sum() < 1:
            logging.critical("Insufficient query structures after fetching", extra=self._log_extra)
            exit_flag = True
        if self._target_df["success"].sum() < 1:
            logging.critical("Insufficient target structures after fetching", extra=self._log_extra)
            exit_flag = True
        if exit_flag:
            exit(1)

    def _alignment(self):
        """
        Coordinate structural alignment routes bridging query and target proteins.

        Dispatches to `_run_foldseek()` if global flag is true, else rolls
        back to `_local_alignment()`. Requisites like `_preprocess_structures()`
        precede foldseek routines.

        Returns:
            None
        """
        self._log_extra.update({"stage": "Alignment"})
        if self._settings["foldseek"]:
            # if not os.path.isfile(self._settings["alignment_path"]):
            logging.info("Preprocessing structures for Foldseek...", extra=self._log_extra)
            self._preprocess_structures()
            logging.info("Running Foldseek easy-search...", extra=self._log_extra)
            self._run_foldseek()
        else:
            logging.info("Running local pairwise aligner...", extra=self._log_extra)
            self._local_alignment()

    def _preprocess_structures(self):
        """
        Split complex PDB structures into single chain mmCIF files for Foldseek processing.

        Instantiates a `StructurePreprocessor` mapping items sourced from the
        local datastores towards isolated files targeting specific sequences.

        Returns:
            None
        """
        stage = {"stage": "Preprocessing Structures"}

        structure_preprocessor = StructurePreprocessor(
            source_dir=self._settings["structure_dir"],
            out_dir=self._settings["foldseek_preprocessed_structure_dir"],
        )
        qt_iter = zip(
            [self._query_df, self._target_df],
            [self._settings["query_dir"], self._settings["target_dir"]],
        )

        for df, search_dir in qt_iter:
            records = df.drop_duplicates(subset=["preprocess_name", "chain_info"]).to_dict(orient="records")
            logging.debug(f"Records to preprocess: {records}", extra=stage)
            results = structure_preprocessor.preprocess_records(records=records, search_dir=search_dir)
            logging.debug(f"Preprocessing results: {results}", extra=stage)

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

    def _run_foldseek(self):
        """
        Execute foldseek sub-commands via bash interfacing.

        Triggers `foldseek easy-search`, pushing input query targets directly against formatted targets
        using `self._settings["alignment_path"]` to store raw tabular matches.

        Returns:
            None
        """
        stage = {"stage": "Foldseek Alignment"}
        if self._settings["target"] == "human_domains":
            self._settings["target_dir"] = self._human_domains_db_path
        cmd = [
            "foldseek",
            "easy-search",
            self._settings["query_dir"],
            self._settings["target_dir"],
            self._settings["alignment_path"],
            self._settings["foldseek_tmp_dir"],
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
                min(3, self._settings["verbosity"])
            ),  # cap foldseek verbosity at 3 (info level) since it can be very verbose at higher levels and we already have our own logging verbosity control
        ]
        logging.debug(f"Running Foldseek with command: {' '.join([str(x) for x in cmd])}", extra=stage)
        subprocess.run(cmd, check=True)
        logging.debug("Foldseek alignment completed successfully", extra=stage)

    def _local_alignment(self):
        """
        Execute traditional sequence-level alignment across inputs.

        Uses the `SequenceAligner` to locally match structured frames, resolving pairwise dependencies
        within internal dataframe representations. Outputs directly to the configuration TSV alignment file.

        Returns:
            None
        """
        # TODO fix local alignment with work with query/target data
        stage = {"stage": "Local Alignment"}
        logging.info("Running local sequence alignments...", extra=stage)
        aligner = SequenceAligner()
        alignment = aligner.align_df(
            self._query_df, self._target_df, self._settings["foldseek_preprocessed_structure_dir"]
        )
        alignment.to_csv(self._settings["alignment_path"], index=False, sep="\t")

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

        pisa_response_dir = os.path.join(self._settings["pocket_dir"], "pisa_responses")
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
                    struct=os.path.join(self._settings["structure_dir"], f"{row['struct_info']}.cif.gz"),
                    chain_id=row["chain_info"].split("_")[0],
                    pocket_residues=[int(x) for x in pisa_pockets[row["pocket_id"]]["res_auth_ids"]],
                    pocket=pisa_pockets[row["pocket_id"]],
                )
        logging.debug(f"Extracted PISA pockets with coords: {pisa_pockets}", extra=self._log_extra)

        with open(os.path.join(self._settings["pocket_dir"], "pisa_pockets.json"), "w") as f:
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
                struct=os.path.join(self._settings["structure_dir"], f"{row['struct_info']}.cif.gz"),
                chain_id=row["chain_info"].split("_")[0],
                pocket_residues=pocket_residues,
            )
            passthrough_pockets[row["pocket_id"]] = pocket
        with open(os.path.join(self._settings["pocket_dir"], "passthrough_pockets.json"), "w") as f:
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
            match row["struct_type"]:
                case "local_file":
                    struct_path = row["struct_info"]
                case _:
                    struct_path = os.path.join(self._settings["structure_dir"], f"{row['struct_info']}.cif.gz")

            pocket = pc.pocket_overlap(
                structure=struct_path,
                domain_chain=row["chain_info"].split("_")[0],
                motif_chain=row["chain_info"].split("_")[1],
            )
            vdw_pockets[row["pocket_id"]] = pocket
        with open(os.path.join(self._settings["pocket_dir"], "vdw_pockets.json"), "w") as f:
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
        alignment_df = pd.read_csv(self._settings["alignment_path"], sep="\t", engine="c")
        blosum_path = os.path.join(os.path.dirname(__file__), "blosum62.bla")

        logging.info(f"{len(alignment_df)} alignment pairs to compare", extra=stage)
        logging.debug(f"Alignment pairs: \n{alignment_df.head()}", extra=stage)

        preproc_to_ids = {}
        for _, row in pd.concat([self._query_df, self._target_df], ignore_index=True).iterrows():
            if row["pocket_id"] in preproc_to_ids:
                if row["pocket_id"] not in preproc_to_ids[row["preprocess_name"]]:
                    preproc_to_ids[row["preprocess_name"]].append(row["pocket_id"])
            else:
                preproc_to_ids[row["preprocess_name"]] = [row["pocket_id"]]
        logging.debug(f"Preprocessed name to pocket ID mapping: {preproc_to_ids}", extra=stage)

        alphafold = (
            self._settings["target"] == "human_domains"
        )  # TODO this is actually if we should trat the target as having no defnied pocket, not if we are searching alphafold
        pockets_df, unknown_alias, incorrect_mapping = lib.compare_pockets(
            alignment_df, pockets, preproc_to_ids=preproc_to_ids, blosum_path=blosum_path, alphafold=alphafold
        )

        # Logging cases where a residue was given a single cahr name unfamiliar to pocketmapper
        if len(unknown_alias) > 0:
            unknown_alias_path = os.path.join(self._settings["results_dir"], "unknown_ids.json")
            logging.warning("Unknown Foldseek Alias, see unknown_alias.json in results directory", extra=stage)
            with open(unknown_alias_path, "w") as f:
                json.dump(lib.jsonify_dict(dict(unknown_alias)), f)

        # logging cases where foldseek mapping had low sequence identity to the parsed structure
        if len(incorrect_mapping) > 0:
            incorrect_mapping_path = os.path.join(self._settings["results_dir"], "incorrect_mapping.json")
            logging.warning("Foldseek mapping with low sequence identity to parsed structure", extra=stage)
            with open(incorrect_mapping_path, "w") as f:
                json.dump(lib.jsonify_dict(dict(incorrect_mapping)), f)

        # Writing pocket comparison results to output file
        output_path = self._settings["pocket_comparison_path"]
        pockets_df.to_csv(output_path, index=False, sep="\t")
        logging.info(f"Pocket comparison results saved to {output_path}", extra=stage)

    def _align_structures(self):
        """
        Perform structural superposition of target structures against the query reference frame.

        If the `align_struct` flag is set, this method will take the top N alignments (as defined by `align_count`)
        from the alignment results and perform a structural alignment using the `StructureAligner` class.
        The aligned structures will be saved to the target directory for downstream analysis.

        Returns:
            None
        """
        stage = {"stage": "Structural Alignment"}
        if not self._settings["align_struct"]:
            logging.info("Structural alignment skipped based on settings", extra=stage)
        else:
            logging.info("Performing structural alignment of target structures...", extra=stage)

        # Pre-loading
        aligner = StructureAligner()
        # alignment_df = pd.read_csv(self._settings["alignment_path"], sep="\t", engine="c")
        pocket_comparison_df = pd.read_csv(self._settings["pocket_comparison_path"], sep="\t", engine="c")

        # For each query structure, align the top N target structures
        for _, row in self._query_df.iterrows():
            query_id = row["pocket_id"]
            logging.debug(f"Processing query {query_id} for structural alignment", extra=stage)
            top_target_ids = (
                pocket_comparison_df.query(f"pocket_1 == '{row['pocket_id']}'")
                .sort_values(by=["min_pct_overlap", "min_overlap_similarity"], ascending=False)
                .head(self._settings["align_count"])
                .loc[:, "pocket_2"]
                .to_list()
            )
            logging.debug(f"Top target IDs for query {query_id}: {top_target_ids}", extra=stage)
            # query_data = {
            #    "name": row["preprocess_name"],
            # }
            aln_ids = [query_id] + top_target_ids
            logging.debug(f"Aligning structures for query {row['pocket_id']}: {aln_ids}", extra=stage)
            if len(aln_ids) > 1:
                aligner.foldseek_transform(
                    aln_ids,
                    self._settings["alignment_path"],
                    self._settings["foldseek_preprocessed_structure_dir"],
                    self._settings["structure_dir"],
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
        ]
        if self._target != "human_domains":  # We don't want to delete the human domains foldseek db if we used it
            tmp_dirs.append("target_dir")
        if self._settings.get("foldseek"):
            tmp_dirs.append("foldseek_tmp_dir")

        # TODO this is unsafe
        for dir in tmp_dirs:
            shutil.rmtree(self._settings[dir])


def main():
    fire.Fire(PocketMapper())


if __name__ == "__main__":
    main()
