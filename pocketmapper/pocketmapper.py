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
from pocketmapper import pisa
from pocketmapper.sequence_aligner import SequenceAligner
from pocketmapper.pocket_calculator import PocketCalculator
from pocketmapper.qt_processor import QTProcessor
from pocketmapper.structure_fetcher import StructureFetcher
from pocketmapper.structure_preprocessor import StructurePreprocessor
from pocketmapper.lib_struct import parse_pocket_from_struct
from pocketmapper import human_domains


class PocketMapper:
    def __init__(self):
        self._log_extra = {"stage": "init"}
        self.human_domains_db_path = files(human_domains).joinpath("human_260310")
        self.log_fmt = "%(levelname)s: %(stage)s - %(msg)s"
        logging.getLogger(__name__)
        logging.basicConfig(level=logging.CRITICAL, format=self.log_fmt)

    def _configure_logging(self, settings):
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
                "standard": {"format": self.log_fmt},
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

    # TODO implement caching option
    def search(
        self,
        query=None,  # settings passed to configure
        target=None,
        settings=None,
        cache_dir=None,
        results_dir=None,
        verbosity=None,  # set verbosity level (see pocketmapper search --help for details)
        help=None,  # help option
        foldseek=None,  # whether to use foldseek for alignment (if false, uses local sequence alignment)
        align_struct=None,  # whether to align structures after pocket comparison
    ):
        """
        Orchestrate and run the full PocketMapper search workflow.
        See pocketmapper search --help for details.
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

        # Main try-except block to catch unhandled exceptions
        self._check_help_search()  # checks if help flag is set and if so prints the help message and exits
        self._configure_workflow()  # configures the settings which have already been read
        self._configure_query_target()  # parses the query and target inputs to determine their types and sets up the relevant data structures for each entry
        self._fetch_missing_structures()  # Fetch any missing structures

        # TODO put preprocessing in alignment call
        if self._settings["foldseek"]:
            self._preprocess_structures()
        self._alignment()  # Align the query and target structures using either local sequence alignment or foldseek based on the settings
        pockets = self._get_pockets()  # Adds seq_pos and cacoords to the pocket info dict

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
        self._delete_tmp()

        logging.info("PocketMapper search completed successfully.", extra={"stage": "End"})

    def _check_help_search(self):
        """
        Displays help information for the PocketMapper tool and exits the program.

        If the 'help' parameter is provided and evaluates to True, this method prints a help message
        describing the usage, options, and features of the PocketMapper package, then terminates execution.

        Parameters:
            help (bool): If True, triggers the display of the help message.

        Usage:
            Call this method when the user requests help (e.g., via a command-line flag).
        """

        if self._help:
            help_msg = """
        PocketMapper - A tool for mapping and analyzing protein pockets.

        Usage:
            pocketmapper search [OPTIONS]

        Primary options (passed to PocketMapper.search):
            --query QUERY            Query identifier or path. Accepts:
                        - 'PDB_CHAIN_CHAIN' (e.g., 1ABC_A_B)
                        - path to a file listing PDB_CHAIN_CHAIN entries (each line)
            --target TARGET          Target identifier or path. Accepts:
                        - 'PDB_CHAIN_CHAIN' (e.g., 2XYZ_C_D)
                        - path to a file listing PDB_CHAIN_CHAIN entries (each line)
                        - special foldseek DB alias 'human_domains' to use the bundled Foldseek DB
            --settings FILE          Path to a JSON settings file (overridden by explicit CLI args)
            --dump_settings FILE     Path to dump the finalized JSON configuration to
            --cache_dir DIR          Directory for caching intermediate files (overrides settings.cache_dir)
            --results_dir DIR        Directory for writing results (overrides settings.results_dir)
            --verbose                Enable more detailed (info) logging
            --debug                  Enable debug-level logging
            --help                   Show this help message and exit

        Relevant settings (can be placed in settings JSON or passed as CLI kwargs):
            cache_dir                Base cache directory (default: pocketmapper_cache)
            results_dir              Results directory (default: pocketmapper_results_<timestamp>)
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

            # Use Foldseek (set foldseek true). When using the built-in human_domains DB:
            pocketmapper search --query 1ABC_A_B --target human_domains --foldseek True --results_dir ./out_fs

            # Override cache and enable debug logging
            pocketmapper search --query 1ABC_A_B --target 2XYZ_C_D --cache_dir /tmp/cache --debug

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
        Reads and configures settings for the PocketMapper workflow.

        1) Sets base defaults
        2) Overrides with provided JSON settings file
        3) Overrides with explicitly passed CLI arguments
        4) Computes derived paths based on the finalized cache/results directories
        5) Dumps configuration to a JSON file in the results directory for record-keeping
        """
        self._log_extra.update({"stage": "Configuring Settings"})

        # 1. Base defaults
        now = datetime.now().strftime("%y%m%d_%H%M%S")
        self._settings = {
            "cache_dir": "pocketmapper_cache",
            "results_dir": f"pocketmapper_results_{now}",
            "query_pocket_method": None,
            "target_pocket_method": None,
            "query": None,
            "target": None,
            "foldseek": False,
            "align_struct": False,
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
            "align_struct": self._align_struct,
            "verbosity": self._verbosity,
        }
        for key, value in cli_args_mapping.items():
            if value is not None:
                self._settings[key] = value

        # Update self properties from settings so they're accessible everywhere consistently
        self._query = self._settings.get("query")
        self._target = self._settings.get("target")

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
        Determining the inputs types for query and target based on the format of the provided values.
        """
        self._log_extra.update({"stage": "Determine Query/Target Types"})

        qtprocessor = QTProcessor(
            query=self._query,
            target=self._target,
            query_pocket_method=self._settings["query_pocket_method"],
            target_pocket_method=self._settings["target_pocket_method"],
        )
        self._query_data, self._target_data = qtprocessor.main()
        logging.debug(f"Query data after processing: \n{self._query_data.head()}", extra=self._log_extra)
        logging.debug(f"Target data after processing: \n{self._target_data.head()}", extra=self._log_extra)

    def _fetch_missing_structures(self):
        """
        1) Downloads structures for the PDB entries in self._pdb_df['interaction_pdb'].

        2) Verifies that structures were found for the query and target entries (if required based on their types)
           and logs the results. If no structures are found for either query or target when required, logs a critical
           error and exits.

        TODO Update to download alphafold structures if the input is a uniprot id
        """
        # Downloading structures
        self._log_extra.update({"stage": "Fetching Missing Structures"})
        logging.info("Starting", extra=self._log_extra)

        structure_fetcher = StructureFetcher(out_dir=self._settings["structure_dir"])
        records = pd.concat([self._query_data, self._target_data]).to_dict(orient="records")
        results = structure_fetcher.get_structures(records)
        logging.debug(f"Structure fetcher results: {results}", extra=self._log_extra)

        success_list = [pdb for pdb, success in results.items() if success]
        failed_list = [pdb for pdb, success in results.items() if not success]

        # Logging results of structure fetching and updating query and target data with success/failure info
        logging.info(
            f"Finished checking for structures. Successfully found structures for {len(success_list)}/{len(results)} PDBs.",
            extra=self._log_extra,
        )
        if len(failed_list) > 0:
            logging.warning(
                f"Failed to find structures for the following PDBs: {', '.join(failed_list)}", extra=self._log_extra
            )

        # updating success column to false if structure fetching failed and adding failure reason
        self._query_data["success"] = self._query_data["struct_info"].map(results).fillna(False)
        self._target_data["success"] = self._target_data["struct_info"].map(results).fillna(False)
        self._query_data.loc[~self._query_data["success"], "failure_reason"] = "structure_not_found"
        self._target_data.loc[~self._target_data["success"], "failure_reason"] = "structure_not_found"

        # Logging the query and target data after structure fetching to verify the updates
        logging.debug(f"Query data after structure fetching: \n{self._query_data.head()}", extra=self._log_extra)
        logging.debug(f"Target data after structure fetching: \n{self._target_data.head()}", extra=self._log_extra)

        # Verifying enough structures were found to continue
        exit_flag = False
        if self._query_data["success"].sum() < 1:
            logging.critical("Insufficient query structures after fetching", extra=self._log_extra)
            exit_flag = True
        if self._target_data["success"].sum() < 1:
            logging.critical("Insufficient target structures after fetching", extra=self._log_extra)
            exit_flag = True
        if exit_flag:
            exit(1)

    def _preprocess_structures(self):
        """
        Divides PDB structures into relevant domains and chains
        """
        self._log_extra.update({"stage": "Preprocess Structures"})

        structure_preprocessor = StructurePreprocessor(
            source_dir=self._settings["structure_dir"], out_dir=self._settings["foldseek_preprocessed_structure_dir"]
        )
        qt_iter = zip(
            [self._query_data, self._target_data], [self._settings["query_dir"], self._settings["target_dir"]]
        )

        for df, search_dir in qt_iter:
            records = df.drop_duplicates(subset=["struct_info", "chain_info"]).to_dict(orient="records")
            logging.debug(f"Records to preprocess: {records}", extra=self._log_extra)
            results = structure_preprocessor.preprocess_records(records=records, search_dir=search_dir)
            logging.debug(f"Preprocessing results: {results}", extra=self._log_extra)

            df.set_index("pocket_id", inplace=True)
            for index, success in results.items():
                if not success:
                    df.loc[index, "success"] = False
                    df.loc[index, "failure_reason"] = "structure_preprocessing_failed"
            df.reset_index(inplace=True)

        logging.info("Finished preprocessing structures", extra=self._log_extra)
        logging.debug(f"Query data after preprocessing: \n{self._query_data.head()}", extra=self._log_extra)
        logging.debug(f"Target data after preprocessing: \n{self._target_data.head()}", extra=self._log_extra)

    def _get_pockets(self):
        """
        Retrieves pockets for the query and target structures based on the specified pocket methods.
        Currently only supports PISA, but can be extended in the future to support other methods.

        For PISA pockets, retrieves the relevant PISA interfaces based on the query and target pocket info and
        extracts the pocket residues and their coordinates from the mmCIF files. The extracted pocket info is then
        stored in a dictionary for later comparison.
        """
        self._log_extra.update({"stage": "Get Pockets"})
        pisa_pockets = self._retrieve_pisa_pockets()
        passthrough_pockets = self.retrieve_passthrough_pockets()
        pockets = pisa_pockets | passthrough_pockets
        logging.debug(f"Combined pockets: {pockets}", extra=self._log_extra)
        return pockets

    def _retrieve_pisa_pockets(self):
        """
        Determines the PISA interfaces from query and targets,
        Retrives the pisa data from the PISA API,
        Parses coordinates for the pockets from the mmCIF files,
        Stores the pocket info in a dict for later comparison
        """

        self._log_extra.update({"stage": "Pisa Pocket Calculation"})
        all_pisa_df = pd.concat([self._query_data, self._target_data], ignore_index=True).query(
            "success and pocket_method == 'pisa'"
        )

        pisa_response_dir = os.path.join(self._settings["pocket_dir"], "pisa_responses")
        # Dispatch to relevant pocket retrieval/calculation method based on the specified pocket method for the query and target (currently only PISA, but can be extended in the future)
        pisa_pdb_list = all_pisa_df["struct_info"].unique().tolist()
        logging.debug(f"PDBs for which to retrieve PISA pockets: {pisa_pdb_list}", extra=self._log_extra)
        logging.info("Retrieving PISA pockets...", extra=self._log_extra)
        downloader = pisa.PisaDownloader()
        downloader.get_interfaces(
            pdb_list=pisa_pdb_list,
            summary_dir=os.path.join(pisa_response_dir, "summaries"),
            asm_dir=os.path.join(pisa_response_dir, "assemblies"),
            interface_dir=os.path.join(pisa_response_dir, "interfaces"),
        )

        # for each interface in query and target extract the relevant info from the parsed pockets
        pisa_pocketids = all_pisa_df["sanitized_pocket_id"].unique().tolist()
        logging.debug(f"Pocket IDs for which to retrieve PISA pockets: {pisa_pocketids}", extra=self._log_extra)
        pisa_pockets = lib.get_pisa_pockets(
            pocket_id_arr=pisa_pocketids,
            in_dir=os.path.join(pisa_response_dir, "interfaces"),
        )
        logging.warning(f"pockets example: {list(pisa_pockets.items())[:2]}", extra=self._log_extra)
        logging.debug(f"Extracted PISA pockets: {pisa_pockets}", extra=self._log_extra)
        for _, row in all_pisa_df.iterrows():
            if row["sanitized_pocket_id"] in pisa_pockets:
                pisa_pockets[row["sanitized_pocket_id"]] = parse_pocket_from_struct(
                    struct=os.path.join(self._settings["structure_dir"], f"{row['struct_info']}.cif.gz"),
                    chain_id=row["chain_info"].split("_")[0],
                    pocket_residues=[int(x) for x in pisa_pockets[row["sanitized_pocket_id"]]["res_auth_ids"]],
                    pocket=pisa_pockets[row["sanitized_pocket_id"]],
                )

        with open(os.path.join(self._settings["pocket_dir"], "pisa_pockets.json"), "w") as f:
            json.dump(pisa_pockets, f)

        return pisa_pockets

    def retrieve_passthrough_pockets(self):
        self._log_extra.update({"stage": "Passthrough Pocket Calculation"})
        all_passthrough_data = pd.concat([self._query_data, self._target_data], ignore_index=True)
        all_passthrough_data = all_passthrough_data.query(
            "success and pocket_method == 'passthrough' and struct_type == 'pdb'"
        )

        # for each pocket in query and target parse pocket info from the structure and store in a dict
        passthrough_pockets = {}
        for _, row in all_passthrough_data.iterrows():
            pocket_residues = [int(x) for x in row["residue_info"].split(",")]
            pocket = parse_pocket_from_struct(
                struct=os.path.join(self._settings["structure_dir"], f"{row['struct_info']}.cif.gz"),
                chain_id=row["chain_info"].split("_")[0],
                pocket_residues=pocket_residues,
            )
            passthrough_pockets[row["sanitized_pocket_id"]] = pocket
        with open(os.path.join(self._settings["pocket_dir"], "passthrough_pockets.json"), "w") as f:
            json.dump(passthrough_pockets, f, indent=4)

        return passthrough_pockets

    def _alignment(self):
        self._log_extra.update({"stage": "Alignment"})
        if self._settings["foldseek"]:
            # if not os.path.isfile(self._settings["alignment_path"]):
            logging.info("Running Foldseek easy-search...", extra=self._log_extra)
            self._run_foldseek()
        else:
            logging.info("Running local alignments...", extra=self._log_extra)
            self._local_alignment()

    def _run_foldseek(self):
        """ """
        if self._target == "human_domains":
            self._settings["target_dir"] = self.human_domains_db_path
        self._log_extra.update({"stage": "Foldseek Alignment"})
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
            r".*_[0-9A-Za-z]\.cif\.gz",
            "--max-seqs",
            "2500",
            "-v",  # verbosity
            "2",  # str(min(3, self._settings["verbosity"])),  # cap foldseek verbosity at 3 (info level) since it can be very verbose at higher levels and we already have our own logging verbosity control
        ]
        logging.debug(f"Running Foldseek with command: {' '.join([str(x) for x in cmd])}", extra=self._log_extra)
        subprocess.run(cmd, check=True)
        logging.debug("Foldseek alignment completed successfully", extra=self._log_extra)

    def _local_alignment(self):
        # TODO fix local calignment with work with query/target data
        aligner = SequenceAligner()
        alignment = aligner.align_df(self._pdb_df, self._settings["foldseek_preprocessed_structure_dir"])
        print(alignment)
        alignment.to_csv(self._settings["alignment_path"], index=False, sep="\t")

    def _compare_pockets_based_on_alignment(self, pockets):
        self._log_extra.update({"stage": "Pocket Comparison"})
        logging.info("Comparing pockets...", extra=self._log_extra)
        alignment_df = pd.read_csv(self._settings["alignment_path"], sep="\t", engine="c")
        blosum_path = os.path.join(os.path.dirname(__file__), "blosum62.bla")

        alphafold = (
            self._target == "human_domains"
        )  # TODO this is actually if we should trat the target as having no defnied pocket, not if we are searching alphafold
        pockets_df, unknown_alias, incorrect_mapping = lib.compare_pockets(
            alignment_df, pockets, blosum_path=blosum_path, alphafold=alphafold
        )

        if len(unknown_alias) > 0:
            unknown_alias_path = os.path.join(self._settings["results_dir"], "unknown_ids.json")
            logging.warning(
                "Unknown Foldseek Alias, see unknown_alias.json in results directory", extra=self._log_extra
            )
            with open(unknown_alias_path, "w") as f:
                json.dump(lib.jsonify_dict(dict(unknown_alias)), f)

        if len(incorrect_mapping) > 0:
            incorrect_mapping_path = os.path.join(self._settings["results_dir"], "incorrect_mapping.json")
            logging.warning("Foldseek mapping with low sequence identity to parsed structure", extra=self._log_extra)
            with open(incorrect_mapping_path, "w") as f:
                json.dump(lib.jsonify_dict(dict(incorrect_mapping)), f)

        # Map sanitized pocket ids back to original pocket ids for output
        sanitized_to_id_map = {}
        for i, row in self._query_data.iterrows():
            sanitized_to_id_map[row["sanitized_pocket_id"]] = row["pocket_id"]
        pockets_df["pocket_1"] = pockets_df["pocket_1"].map(sanitized_to_id_map)

        # Writing pocket comparison results to output file
        output_path = self._settings["pocket_comparison_path"]
        pockets_df.to_csv(output_path, index=False, sep="\t")
        with open(os.path.join(self._settings["results_dir"], "sanitized_to_id_map.json"), "w") as f:
            json.dump(sanitized_to_id_map, f, indent=2)
        logging.info(f"Pocket comparison results saved to {output_path}", extra=self._log_extra)

    def _delete_tmp(self):
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
