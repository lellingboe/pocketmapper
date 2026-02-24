"""
Functions
given target and query structures...
fetch the files (optional)
calculate pockets
store calculated pockets
make domain only and domain-motif structures
Align these structures
"""

import fire
import logging
import json
import subprocess
import pandas as pd
import os
import re
from datetime import datetime
import shutil
import gemmi
from importlib.resources import files
from pocketmapper import lib
from pocketmapper import pisa
from pocketmapper.sequence_aligner import SequenceAligner
from pocketmapper.pocket_calculator import PocketCalculator
from pocketmapper import human_domains


class PocketMapper:
    def __init__(self):
        self._settings = {}
        self._stage = {"stage": "init"}
        self._requires_structures = ["pdb_chain_chain", "file"]
        self._pdb_df = None

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
        help=None,  # help option
        foldseek=None,
        pocket_method=None,  # method to calculate pockets, default is PISA
        align_struct=None,  # whether to align structures after pocket comparison
    ):
        """
        Orchestrate and run the full PocketMapper search workflow.
        See pocketmapper search --help for details.
        """
        self._stage = {
            "stage": "Start"
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

        # Main try-except block to catch unhandled exceptions
        try:
            # Setting up things
            self._search_help()
            self._setup_logging()
            self._configure(  # configures the settings which have already been read
                cache_dir=cache_dir,
                results_dir=results_dir,
                query=query,
                target=target,
                foldseek=foldseek,
            )
            self._setup_query_target()

            # Preparing structures for later
            self._prepare_directories()
            self._fetch_and_verify_structures()
            self._divide_structures()

            self._alignment()

            if True:  # hack to make ATP pocket search working
                pockets = self._retrieve_pockets()
                pockets = self._get_atom_coords_from_cif(pockets)  # Adds seq_pos and cacoords to the pocket info dict
            else:
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

        if self._help:
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

    def _configure(self, **kwargs):
        """
        Reads and configures settings for the PocketMapper workflow.

        1) Reads setting from provided JSON
        2) Overrides settings with any provided as explicit keyword arguments
        3) Sets default values for any settings not provided in either of the above steps
        """
        self._stage.update({"stage": "Configuring Settings"})

        # Populates settings from the settings file if provided
        if self._settings_file is not None:
            # Checking that the supplied settings file exists
            if not os.path.isfile(self._settings_file):
                logging.critical(f"Settings file not found: {self._settings_file}", extra=self._stage)
                exit(1)
            # Try to read the settings file and update the settings dict, catching JSON parsing errors
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
        for key, value in kwargs.items():
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
            # Other settings
            "foldseek": False,
            "pisa_pockets": True,
            "structure": False,
        }
        for key, value in defaults.items():
            if key not in self._settings:
                self._settings[key] = value

        logging.debug(f"\nInternal settings dictionary:\n{self._settings}", extra=self._stage)

    def _setup_query_target(self):
        """
        Determining the inputs types for query and target based on the format of the provided values.
        """
        self._stage.update({"stage": "Determine Query/Target Types"})

        # What is determined in this stage
        self._query_type = None  # {file, pdb_chain_chain}
        self._target_type = None  # {file, foldseek_db, pdb_chain_chain}
        self._pdb_df = (
            None  # DataFrame to hold all the relevant info about the query and target PDBs, chains, motifs, etc.
        )

        query = self._settings.get("query")
        if os.path.isfile(query):
            self._query_type = "file"
        else:
            self._query_type = "pdb_chain_chain"

        target = self._settings.get("target")
        if os.path.isfile(target):
            if ".txt" in target:
                self._target_type = "file"
        elif target.lower() == "human_domains":
            self._target_type = "foldseek_db"
            self._settings["target_dir"] = files(human_domains).joinpath(
                "human"
            )  # Overriding target dir to point to the db
        else:
            self._target_type = "pdb_chain_chain"

        self.logger.debug(f"Determined query type: {self._query_type}", extra=self._stage)
        self.logger.debug(f"Determined target type: {self._target_type}", extra=self._stage)

        # Checking a target and query is specified
        self._stage.update({"stage": "Input Validation"})
        if not self._settings.get("query"):
            raise ValueError("No query specified. Use --query")
        if not self._settings.get("target"):
            raise ValueError("No target specified. Use --target")

        # Reading the data
        self._stage.update({"stage": "Reading Input Data"})
        query_data = self._parse_inputs(self._query, self._query_type, type="query")
        if type(query_data) is pd.DataFrame:
            self.logger.debug(f"query_data:\n{query_data.head()}", extra=self._stage)
        else:
            self.logger.debug(f"query_data:\n{query_data}", extra=self._stage)

        target_data = self._parse_inputs(self._target, self._target_type, type="target")
        if type(target_data) is pd.DataFrame:
            self.logger.debug(f"target_data:\n{target_data.head()}", extra=self._stage)
        else:
            self.logger.debug(f"target_data:\n{target_data}", extra=self._stage)

        """
        query_df = self._make_tq_df("query", "query_file", type="query")
        target_df = self._make_tq_df("target", "target_file", type="target")
        """
        pdb_data = [query_data]
        if type(target_data) is pd.DataFrame:
            pdb_data.append(target_data)
        self._pdb_df = pd.concat(pdb_data, ignore_index=True)

        logging.debug(f"Combined Query/Target DataFrame:\n{self._pdb_df.head(10)}", extra=self._stage)

    def _parse_inputs(self, value, value_type, **kwargs):
        pdb_chain_chain_re = re.compile(r"[A-Za-z0-9]{4}_[A-Za-z0-9]_[A-Za-z0-9]")
        if value_type == "file":
            pdbs = []
            domains = []
            motifs = []
            valid_lines = 0
            with open(value) as f:
                for i, line in enumerate(f.readlines()):
                    if line[0] == "#":
                        continue
                    line = line.strip()
                    if pdb_chain_chain_re.match(line):
                        pdb, domain, motif = line.split("_")
                        pdbs.append(pdb)
                        domains.append(domain)
                        motifs.append(motif)
                        valid_lines += 1
                    else:
                        self.logger.warning(
                            f"Line {i + 1} in file '{value}' does not match PDB_CHAIN_CHAIN format: {line}",
                            extra=self._stage,
                        )
            if valid_lines == 0:
                self.logger.critical(f"No valid lines found in file '{value}'.", extra=self._stage)
                exit(1)
            df = self.pdb_info_to_df(
                pdbs=pdbs,
                domains=domains,
                motifs=motifs,
                **kwargs,
            )
            return df
        elif value_type == "pdb_chain_chain":
            if not pdb_chain_chain_re.match(value):
                self.logger.critical(f"'{value}' does not match required format 'PDB_CHAIN_CHAIN'.", extra=self._stage)
                exit(1)
            pdb, domain, motif = value.split("_")
            df = self.pdb_info_to_df(
                pdbs=[pdb],
                domains=[domain],
                motifs=[motif],
                **kwargs,
            )
            return df
        elif value_type == "foldseek_db":
            # For foldseek db, the value is the path to the db
            return value
        else:
            self.logger.critical(f"Unknown value_type '{value_type}'", extra=self._stage)
            exit(1)

    def pdb_info_to_df(self, pdbs, domains, motifs, **kwargs):
        """
        Docstring for pdb_info_to_df

        :param self: Description
        :param pdbs: Description
        :param domains: Description
        :param motifs: Description
        :param kwargs: Description
        """
        data = {
            "interaction_pdb": [pdb.upper() for pdb in pdbs],
            "domain_chain": domains,
            "motif_chain": motifs,
        }

        df = pd.DataFrame.from_dict(data)
        df["pdb_domain"] = df.apply(lambda x: x.interaction_pdb + "_" + x.domain_chain, axis=1)
        df["pdb_domain_motif"] = df.apply(lambda x: x.pdb_domain + "_" + x.motif_chain, axis=1)
        for k, v in kwargs.items():
            df[k] = v
        return df

    def _prepare_directories(self):
        self._stage.update({"stage": "Directory Preparation"})
        dirs_to_create = [
            "structure_dir",
            "query_dir",
            "pocket_dir",
            "pisa_dir",
            "divided_struct_dir",
            "results_dir",
        ]
        if self._target_type != "foldseek_db":
            dirs_to_create.append("target_dir")

        for dir_key in dirs_to_create:
            path = self._settings[dir_key]
            try:
                os.makedirs(path, exist_ok=True)
            except OSError as e:
                raise OSError(f"Error creating directory {path}: {e}")

    def _fetch_and_verify_structures(self):
        """
        1) Downloads structures for the PDB entries in self._pdb_df['interaction_pdb'].

        2) Verifies that structures were found for the query and target entries (if required based on their types)
           and logs the results. If no structures are found for either query or target when required, logs a critical
           error and exits.
        """
        # Downloading structures
        self._stage.update({"stage": "Fetch Structures"})
        logging.info("Checking for mmCIF structures...", extra=self._stage)
        found_map = lib.get_mmcifs(
            pdb_list=self._pdb_df["interaction_pdb"].unique(),
            out_dir=self._settings["structure_dir"],
        )
        self._pdb_df["structure_found"] = self._pdb_df["interaction_pdb"].map(found_map)

        # Verifying enough structures were found to continue
        self._stage.update({"stage": "Verify Structures"})
        # Query
        if self._query_type in self._requires_structures:
            if not self._pdb_df.query("structure_found and type == 'query'").empty:
                logging.info("Query structures found.", extra=self._stage)
            else:
                logging.critical("No query structures found locally or via download", extra=self._stage)
                exit(1)
        # Target
        if self._target_type in self._requires_structures:
            if not self._pdb_df.query("structure_found and type == 'target'").empty:
                logging.info("Target structures found.", extra=self._stage)
            else:
                logging.critical("No target structures found locally or via download", extra=self._stage)
                exit(1)

    def _divide_structures(self):
        self._stage.update({"stage": "Preprocess Structures"})
        logging.info("Dividing mmCIF structures...", extra=self._stage)
        divided_map = lib.pdb_preprocessing_gemmi(
            df=self._pdb_df.query("structure_found"),
            ref_dir=self._settings["structure_dir"],
            cache_dir=self._settings["divided_struct_dir"],
            query_dir=self._settings["query_dir"],
            target_dir=self._settings["target_dir"],
        )
        self._pdb_df["divided_struct"] = self._pdb_df["pdb_domain_motif"].map(divided_map).fillna(False)

    def _retrieve_pockets(self):
        self._stage.update({"stage": "Pocket Calculation"})

        # WRITING PISA POCKETS
        logging.info("Retrieving PISA pockets...", extra=self._stage)
        downloader = pisa.PisaDownloader()
        downloader.get_interfaces(
            pdb_list=self._pdb_df.query("divided_struct")["interaction_pdb"].str.lower().unique(),
            summary_dir=os.path.join(self._settings["pisa_dir"], "summaries"),
            asm_dir=os.path.join(self._settings["pisa_dir"], "assemblies"),
            interface_dir=os.path.join(self._settings["pisa_dir"], "interfaces"),
        )
        pisa_pockets = lib.get_pisa_pockets(
            df=self._pdb_df.query("divided_struct"),
            in_dir=os.path.join(self._settings["pisa_dir"], "interfaces"),
            out_dir=self._settings["pocket_dir"],
        )
        with open(os.path.join(self._settings["pisa_dir"], "all_pockets_1.json"), "w") as f:
            json.dump(pisa_pockets, f)

        return pisa_pockets

    def _get_atom_coords_from_cif(self, pockets):
        self._stage.update({"stage": "Getting atom coords"})
        for pocket_id, pocket in pockets.items():
            try:
                struct_path = os.path.join(self._settings["divided_struct_dir"], f"{pocket_id[:-2]}.cif.gz")
                st = gemmi.read_structure(struct_path, format=gemmi.CoorFormat.Mmcif)
                domain_chain = st[0][pocket_id[5]]
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

        with open(os.path.join(self._settings["pisa_dir"], "all_pockets_2.json"), "w") as f:
            json.dump(pockets, f)
        return pockets

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
            "3",
        ]
        subprocess.run(cmd, check=True)

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

        alphafold = self._target_type == "foldseek_db"
        pockets_df, unknown_alias = lib.compare_pockets(
            alignment_df, pockets, blosum_path=blosum_path, alphafold=alphafold
        )

        if len(unknown_alias) > 0:
            unknown_alias_path = os.path.join(self._settings["results_dir"], "unknown_ids.json")
            logging.warning("Unknown Foldseek Alias, see unknown_alias.json in results directory", extra=self._stage)
            with open(unknown_alias_path, "w") as f:
                json.dump(lib.jsonify_dict(dict(unknown_alias)), f)

        output_path = self._settings["pocket_comparison_path"]
        pockets_df.to_csv(output_path, index=False, sep="\t")
        logging.info(f"Pocket comparison results saved to {output_path}", extra=self._stage)

    def _delete_tmp(self):
        tmp_dirs = [
            "query_dir",
        ]
        if self._target_type != "foldseek_db":
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
