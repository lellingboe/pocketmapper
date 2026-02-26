"""
PocketMapper: A tool for mapping and analyzing protein pockets.

Author: Lachlan Ellingboe

"""

from importlib.resources import files
import fire
import logging
import json
import subprocess
import pandas as pd
import os
from datetime import datetime
import shutil
import gemmi
from pocketmapper import lib
from pocketmapper import pisa
from pocketmapper.sequence_aligner import SequenceAligner
from pocketmapper.pocket_calculator import PocketCalculator
from pocketmapper.qt_processor import QTProcessor
from pocketmapper import human_domains


class PocketMapper:
    def __init__(self):
        self._settings = {}
        self._stage = {"stage": "init"}
        self._requires_structures = ["pdb_chain_chain", "file"]

    # TODO implement caching option
    def search(
        self,
        query=None,  # settings passed to configure
        target=None,
        settings=None,
        cache_dir=None,
        results_dir=None,
        verbose=False,  # makes log file more verbose
        debug=False,  # make log file even more verbose
        help=False,  # help option
        query_pocket_method=None,  # method to calculate pockets, default is PISA
        target_pocket_method=None,  # method to calculate pockets, default is PISA
        foldseek=None,  # whether to use foldseek for alignment (if false, uses local sequence alignment)
        align_struct=None,  # whether to align structures after pocket comparison
    ):
        """
        Orchestrate and run the full PocketMapper search workflow.
        See pocketmapper search --help for details.
        """
        self._stage = {
            "stage": "Starting Search"
        }  # dict needed for logging extra info, can be updated throughout the process to indicate the current stage in logs

        # Storing input parameters
        self._query = query
        self._target = target
        self._settings_file = settings
        self._cache_dir = cache_dir
        self._results_dir = results_dir
        self._verbose = verbose
        self._debug = debug
        self._help = help
        self._foldseek = foldseek
        self._query_pocket_method = query_pocket_method
        self._target_pocket_method = target_pocket_method
        self._align_struct = align_struct

        self.human_domains_db_path = files(human_domains).joinpath("human")

        # Main try-except block to catch unhandled exceptions
        try:
            # Setting up things
            self._search_help()
            self._setup_logging()
            self._configure()  # configures the settings which have already been read

            # Preprocessing/Downloading required data
            self._setup_query_target()
            self._prepare_directories()
            self._fetch_pdb_structures()
            if self._foldseek:
                self._preprocess_structures()

            # Alignment
            self._alignment()

            # Pockets
            pockets = self._get_pockets()  # Adds seq_pos and cacoords to the pocket info dict
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

            self._compare_pockets_and_save(pockets)
            self._delete_tmp()

            logging.info("PocketMapper search completed successfully.", extra={"stage": "End"})

        # Unhandled exception stops the process and logs the error
        except Exception as e:
            logging.exception(str(e), extra=self._stage)
            exit(1)

    def align(self):
        """
        Docstring for align

        To be implemented, allow you to aligne 2 structures with foldseek or local alignment without doing the full pocketmapper search workflow
        """
        pass

    def _search_help(self):
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
            pisa_dir                 Directory for PISA related files
            divided_struct_dir       Directory for preprocessed/divided structures
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

    def _setup_logging(self):
        """
        Sets up logging configuration for the PocketMapper workflow.
        The logging level is determined based on the 'debug' and 'verbose' flags provided during initialization.
        Logs are output to both the console and a file named 'info.log' in the current directory.
        The log format includes the log level, stage, and message.

        """
        self._stage = {"stage": "Logging Setup"}  # Updating stage for logging context

        # Determining log level based on debug and verbose flags
        if self._debug:
            log_level = logging.DEBUG
        elif self._verbose:
            log_level = logging.INFO
        else:
            log_level = logging.WARNING

        # Configuring logging to console
        fmt = "%(levelname)s: %(stage)s - %(msg)s"
        self.logger = logging.getLogger("pocketmapper")
        logging.basicConfig(level=log_level, format=fmt)

        # Configuring logging to file
        # TODO make log file name include timestamp and maybe some info about the run (e.g., query and target)
        # TODO allow use to specify log file name and location via settings
        # TODO use output directory
        formatter = logging.Formatter(fmt)
        fh = logging.FileHandler("info.log")
        fh.setLevel(log_level)
        fh.setFormatter(formatter)
        self.logger.addHandler(fh)

        self.logger.debug(
            "Level set to DEBUG", extra=self._stage
        )  # example of how to use the logger with the stage info

    def _configure(self):
        """
        Reads and configures settings for the PocketMapper workflow.

        1) Reads setting from provided JSON
        2) Overrides settings with any provided as explicit keyword arguments
        3) Sets default values for any settings not provided in either of the above steps
        """
        self._stage.update({"stage": "Configuring Settings"})

        # Populates settings from the settings file if provided
        if self._settings_file is not None:
            if not os.path.isfile(self._settings_file):
                logging.critical(f"Settings file not found: {self._settings_file}", extra=self._stage)
                exit(1)
            try:
                with open(self._settings_file) as f:
                    settings_data_from_file = json.load(f)
            except Exception:
                logging.exception(
                    f"Error reading settings file: {self._settings_file}. Is it in JSON format?", extra=self._stage
                )
                exit(1)
            finally:
                self._settings.update(settings_data_from_file)

        # Override settings_file with any provided command-line arguments
        cdm_line_args = {
            "query": self._query,
            "target": self._target,
            "cache_dir": self._cache_dir,
            "results_dir": self._results_dir,
            "foldseek": self._foldseek,
            "query_pocket_method": self._query_pocket_method,
            "target_pocket_method": self._target_pocket_method,
            "align_struct": self._align_struct,
        }
        for key, value in cdm_line_args.items():
            if value is not None:
                self._settings[key] = value

        # Using defaults setting for any settings not provided by either of the above steps
        cache_dir = self._settings.get("cache_dir", "pocketmapper_cache")
        now = datetime.now().strftime("%y%m%d_%H%M%S")
        results_dir = self._settings.get("results_dir", f"pocketmapper_results_{now}")
        defaults = {
            # Default directories stemming from cache_dir
            "cache_dir": cache_dir,
            "structure_dir": os.path.join(cache_dir, "pdb_structures"),
            "pocket_dir": os.path.join(cache_dir, "pockets"),
            "foldseek_tmp_dir": os.path.join(cache_dir, "foldseek_tmp"),
            "pisa_dir": os.path.join(cache_dir, "pisa_pockets"),
            "divided_struct_dir": os.path.join(cache_dir, "divided_structs"),
            # Default directories and files stemming from results_dir
            "results_dir": results_dir,
            "query_dir": os.path.join(results_dir, "query_structures"),
            "target_dir": os.path.join(results_dir, "target_structures"),
            "alignment_path": os.path.join(results_dir, "alignment.tsv"),
            "pocket_comparison_path": os.path.join(results_dir, "pocket_comparison.tsv"),
            "foldseek": False,
            "query_pocket_method": None,
            "target_pocket_method": None,
            "align_struct": False,
        }
        for key, value in defaults.items():
            if key not in self._settings:
                self._settings[key] = value

        logging.debug("Internal settings:", extra=self._stage)
        for k, v in self._settings.items():
            logging.debug(f"{k}: {v}", extra=self._stage)

    def _setup_query_target(self):
        """
        Determining the inputs types for query and target based on the format of the provided values.
        """
        self._stage.update({"stage": "Determine Query/Target Types"})

        qtprocessor = QTProcessor(
            query=self._query,
            target=self._target,
            query_pocket_method=self._settings["query_pocket_method"],
            target_pocket_method=self._settings["target_pocket_method"],
        )
        self._query_data, self._target_data = qtprocessor.main()

    def _prepare_directories(self):
        self._stage.update({"stage": "Directory Preparation"})
        dirs_to_create = [
            "structure_dir",
            "query_dir",
            "target_dir",
            "pocket_dir",
            "divided_struct_dir",
        ]

        for dir_key in dirs_to_create:
            path = self._settings[dir_key]
            try:
                os.makedirs(path, exist_ok=True)
            except OSError:
                logging.critical(f"Error creating directory {path}", extra=self._stage)
                exit(1)

    def _fetch_pdb_structures(self):
        """
        1) Downloads structures for the PDB entries in self._pdb_df['interaction_pdb'].

        2) Verifies that structures were found for the query and target entries (if required based on their types)
           and logs the results. If no structures are found for either query or target when required, logs a critical
           error and exits.

        TODO Update to download alphafold structures if the input is a uniprot id
        """
        # Downloading structures
        self._stage.update({"stage": "Fetching Structures"})
        logging.info("Checking for mmCIF structures...", extra=self._stage)

        # List of all unique PDBs
        query_pdbs = set(self._query_data.query("struct_type == 'pdb'")["struct_info"])
        target_pdbs = set(self._target_data.query("struct_type == 'pdb'")["struct_info"])
        all_pdbs = list(query_pdbs.union(target_pdbs))
        success_map = lib.get_mmcifs(
            pdb_list=all_pdbs,
            out_dir=self._settings["structure_dir"],
        )
        failed_list = [pdb for pdb, success in success_map.items() if not success]

        # Logging results of structure fetching and updating query and target data with success/failure info
        logging.info(
            f"Finished checking for structures. Successfully found structures for {len(all_pdbs) - len(failed_list)}/{len(all_pdbs)} PDBs.",
            extra=self._stage,
        )
        if len(failed_list) > 0:
            logging.warning(
                f"Failed to find structures for the following PDBs: {', '.join(failed_list)}", extra=self._stage
            )

        # updating success column to false if structure fetching failed and adding failure reason
        self._query_data["success"] = ~self._query_data["struct_info"].isin(failed_list)
        self._target_data["success"] = ~self._target_data["struct_info"].isin(failed_list)
        self._query_data.loc[~self._query_data["success"], "failure_reason"] = "structure_not_found"
        self._target_data.loc[~self._target_data["success"], "failure_reason"] = "structure_not_found"

        # Logging the query and target data after structure fetching to verify the updates
        self.logger.debug(f"Query data after structure fetching: \n{self._query_data.head()}", extra=self._stage)
        self.logger.debug(f"Target data after structure fetching: \n{self._target_data.head()}", extra=self._stage)

        # Verifying enough structures were found to continue
        exit_flag = False
        if self._query_data["success"].sum() < 1:
            logging.critical("Insufficient query after fetching structures found", extra=self._stage)
            exit_flag = True
        if self._target_data["success"].sum() < 1:
            logging.critical("Insufficient targets after fetching structures found", extra=self._stage)
            exit_flag = True
        if exit_flag:
            exit(1)

    def _preprocess_structures(self):
        """
        Divides PDB structures into relevant domains and chains
        """
        self._stage.update({"stage": "Preprocess Structures"})
        logging.info("Dividing mmCIF structures...", extra=self._stage)
        query_divided_map = lib.pdb_preprocessing_gemmi(
            df=self._query_data.query("success and struct_type == 'pdb'"),
            ref_dir=self._settings["structure_dir"],
            cache_dir=self._settings["divided_struct_dir"],
            out_dir=self._settings["query_dir"],
        )
        for index, success in query_divided_map.items():
            if not success:
                self._query_data.loc[index, "success"] = False
                self._query_data.loc[index, "failure_reason"] = "structure_preprocessing_failed"
        target_divided_map = lib.pdb_preprocessing_gemmi(
            df=self._target_data.query("success and struct_type == 'pdb'"),
            ref_dir=self._settings["structure_dir"],
            cache_dir=self._settings["divided_struct_dir"],
            out_dir=self._settings["target_dir"],
        )
        for index, success in target_divided_map.items():
            if not success:
                self._target_data.loc[index, "success"] = False
                self._target_data.loc[index, "failure_reason"] = "structure_preprocessing_failed"
        logging.info("Finished dividing structures", extra=self._stage)
        logging.debug(f"Query data after structure preprocessing: \n{self._query_data.head()}", extra=self._stage)
        logging.debug(f"Target data after structure preprocessing: \n{self._target_data.head()}", extra=self._stage)

    def _get_pockets(self):
        """
        Retrieves pockets for the query and target structures based on the specified pocket methods.
        Currently only supports PISA, but can be extended in the future to support other methods.

        For PISA pockets, retrieves the relevant PISA interfaces based on the query and target pocket info and
        extracts the pocket residues and their coordinates from the mmCIF files. The extracted pocket info is then
        stored in a dictionary for later comparison.
        """
        pockets = {}
        pisa_pockets = self._retrieve_pisa_pockets()
        pockets.update(pisa_pockets)
        passthrough_pockets = self.retrieve_passthrough_pockets()
        pockets.update(passthrough_pockets)
        return pockets

    def _retrieve_pisa_pockets(self):
        self._stage.update({"stage": "Pisa Pocket Calculation"})
        all_pisa_data = pd.concat([self._query_data, self._target_data], ignore_index=True)
        all_pisa_data = all_pisa_data.query("success and pocket_method == 'pisa' and struct_type == 'pdb'")

        # Dispatch to relevant pocket retrieval/calculation method based on the specified pocket method for the query and target (currently only PISA, but can be extended in the future)
        pisa_pdb_list = all_pisa_data["struct_info"].unique().tolist()
        logging.debug(f"PDBs for which to retrieve PISA pockets: {pisa_pdb_list}", extra=self._stage)
        logging.info("Retrieving PISA pockets...", extra=self._stage)
        downloader = pisa.PisaDownloader()
        downloader.get_interfaces(
            pdb_list=pisa_pdb_list,
            summary_dir=os.path.join(self._settings["pisa_dir"], "summaries"),
            asm_dir=os.path.join(self._settings["pisa_dir"], "assemblies"),
            interface_dir=os.path.join(self._settings["pisa_dir"], "interfaces"),
        )

        # for each interface in query and target extract the relevant info from the parsed pockets
        query_pisa_pocketid_list = all_pisa_data["pocket_id"].unique().tolist()
        all_pisa_pocket_ids = list(set(query_pisa_pocketid_list))
        logging.debug(f"Pocket IDs for which to retrieve PISA pockets: {all_pisa_pocket_ids}", extra=self._stage)
        pisa_pockets = lib.get_pisa_pockets(
            pocket_id_arr=all_pisa_pocket_ids,
            in_dir=os.path.join(self._settings["pisa_dir"], "interfaces"),
        )
        pisa_pockets = self._add_atom_coords_to_pisa_pockets(
            pisa_pockets
        )  # Adds seq_pos and cacoords to the pocket info dict

        with open(os.path.join(self._settings["pisa_dir"], "pisa_pockets.json"), "w") as f:
            json.dump(pisa_pockets, f)

        return pisa_pockets

    def _add_atom_coords_to_pisa_pockets(self, pockets):
        """
        Adds atom coords to pisa pockets
        """
        self._stage.update({"stage": "Getting atom coords"})
        for pocket_id, pocket in pockets.items():
            try:
                struct_path = os.path.join(self._settings["structure_dir"], f"{pocket_id.split(':')[0]}.cif.gz")
                st = gemmi.read_structure(struct_path, format=gemmi.CoorFormat.Mmcif)
                domain_chain = st[0][pocket_id.split(":")[1].split("_")[0]]
                pocket_keys = pocket.keys()
                seq_pos = 0
                for res in domain_chain:
                    # mapping ca_seq position
                    ca_atom = res.get_ca()
                    res_id = str(res.seqid.num)

                    # If the residue has a CA atom specified, save the info to the pocket and
                    if ca_atom is not None:
                        if res_id in pocket_keys:
                            pockets[pocket_id][res_id]["seq_pos"] = seq_pos
                            pockets[pocket_id][res_id]["ca_coords"] = list(res.get_ca().pos)
                        seq_pos += 1
                    else:
                        if res_id in pocket_keys:
                            msg = (
                                f"Pocket residue {res_id} in {pocket_id} "
                                "does not have CA coords and will be excluded from the comparison"
                            )
                            logging.warning(
                                msg,
                                extra=self._stage,
                            )
                            pockets[pocket_id][res_id]["seq_pos"] = -1  # Removes it from later comparison
                pockets[pocket_id]["has_coords"] = True

            except Exception:
                logging.warning(f"Error getting coords for {pocket_id}", extra=self._stage)
                pockets[pocket_id]["has_coords"] = False
        return pockets

    def retrieve_passthrough_pockets(self):
        self._stage.update({"stage": "Passthrough Pocket Calculation"})
        all_passthrough_data = pd.concat([self._query_data, self._target_data], ignore_index=True)
        all_passthrough_data = all_passthrough_data.query(
            "success and pocket_method == 'passthrough' and struct_type == 'pdb'"
        )
        passthrough_pockets = lib.passthrough_pockets(all_passthrough_data, self._settings["structure_dir"])

        with open(os.path.join(self._settings["pocket_dir"], "passthrough_pockets.json"), "w") as f:
            json.dump(passthrough_pockets, f, indent=4)

        return passthrough_pockets

    def _alignment(self):
        self._stage.update({"stage": "Alignment"})
        if self._settings["foldseek"]:
            logging.info("Running Foldseek easy-search...", extra=self._stage)
            self._run_foldseek()
        else:
            logging.info("Running local alignments...", extra=self._stage)
            self._local_alignment()

    def _run_foldseek(self):
        """ """
        if self._target == "human_domains":
            self._settings["target_dir"] = self.human_domains_db_path
        self._stage.update({"stage": "Foldseek Alignment"})
        cmd = [
            "foldseek",
            "easy-search",
            self._settings["query_dir"],
            self._settings["target_dir"],
            self._settings["alignment_path"],
            self._settings["foldseek_tmp_dir"],
            "--format-output",
            "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,lddt,qaln,taln,u,t",
            "--format-mode",
            "4",
            "-e",
            "0.001",
            "--file-include",
            r"[0-9A-Z]{4}_[0-9A-Za-z]\.cif\.gz",
            "--max-seqs",
            "2500",
            "-v",  # verbosity
            "2",
        ]
        self.logger.debug(f"Running Foldseek with command: {' '.join([str(x) for x in cmd])}", extra=self._stage)
        subprocess.run(cmd, check=True)
        self.logger.debug("Foldseek alignment completed successfully", extra=self._stage)

    def _local_alignment(self):
        aligner = SequenceAligner()
        alignment = aligner.align_df(self._pdb_df, self._settings["divided_struct_dir"])
        print(alignment)
        alignment.to_csv(self._settings["alignment_path"], index=False, sep="\t")

    def _compare_pockets_and_save(self, pockets):
        self._stage.update({"stage": "Pocket Comparison"})
        logging.info("Comparing pockets...", extra=self._stage)
        alignment_df = pd.read_csv(self._settings["alignment_path"], sep="\t", engine="c")
        blosum_path = os.path.join(os.path.dirname(__file__), "blosum62.bla")

        alphafold = (
            self._target == "human_domains"
        )  # TODO this is actually if we should trat the target as having no defnied pocket, not if we are searching alphafold
        pockets_df, unknown_alias = lib.compare_pockets(
            alignment_df, pockets, blosum_path=blosum_path, alphafold=alphafold
        )

        if len(unknown_alias) > 0:
            unknown_alias_path = os.path.join(self._settings["results_dir"], "unknown_ids.json")
            logging.warning("Unknown Foldseek Alias, see unknown_alias.json in results directory", extra=self._stage)
            with open(unknown_alias_path, "w") as f:
                json.dump(lib.jsonify_dict(dict(unknown_alias)), f)

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
        logging.info(f"Pocket comparison results saved to {output_path}", extra=self._stage)

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
