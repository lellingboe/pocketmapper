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
import lib
import subprocess
import pandas as pd
import os
import re
from datetime import datetime
import pisa
import shutil
import gemmi
from importlib.resources import files


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
        **kwargs,
    ):
        """
        Main orchestration method to run the pocket mapping workflow.
        """
        self._stage = {"stage": "Start"}

        # Storing input parameters
        self._query = query
        self._target = target
        self._settings_file = settings
        self._cache_dir = cache_dir
        self._results_dir = results_dir
        self._verbose = verbose
        self._debug = debug
        self._help = help
        self._uncaught_args = kwargs

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
            )
            self._determine_query_target_types()
            self._check_query_target_set()
            self._read_query_target_data()

            # Preparing structures for later
            self._prepare_directories()
            self._fetch_and_verify_structures()
            self._divide_structures()

            self._run_foldseek()

            pockets = self._retrieve_pockets()
            pockets = self._get_atom_coords_from_cif(pockets)

            self._compare_pockets_and_save(pockets)
            self._delete_tmp()

            logging.info("PocketMapper search completed successfully.", extra={"stage": "End"})

        # Unhandled exception stops the process and logs the error
        except Exception as e:
            logging.exception(str(e), extra=self._stage)
            exit(1)

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

Options:
    --query QUERY            Specify the query in 'PDB_CHAIN_CHAIN' format (e.g., 1ABC_A_B)
    --target TARGET          Specify the target in 'PDB_CHAIN_CHAIN' format (e.g., 2XYZ_C_D)
    --query_file FILE        Path to a TSV file listing query structures
    --target_file FILE       Path to a TSV file listing target structures
    --settings FILE          Path to a JSON settings file
    --cache_dir DIR          Directory for caching intermediate files (default: pocketmapper_cache)
    --results_dir DIR        Directory for output results (default: pocketmapper_results_<timestamp>)
    --verbose                Enable more detailed logging
    --debug                  Enable debug-level logging
    --help                   Show this help message and exit

Description:
    PocketMapper orchestrates the workflow for fetching structures, calculating and comparing protein pockets,
    and aligning structures using Foldseek. It supports both single structure and batch (file-based) modes.

Examples:
    pocketmapper search --query 1ABC_A_B --target 2XYZ_C_D
    pocketmapper search --query_file queries.tsv --target_file targets.tsv --settings config.json

For more information, visit: https://github.com/yourorg/pocketmapper
                """

        if self._help:
            print(help_msg)
            exit()

    def _setup_logging(self):
        self._stage = {"stage": "Logging Setup"}
        if self._debug:
            log_level = logging.DEBUG
        elif self._verbose:
            log_level = logging.INFO
        else:
            log_level = logging.WARNING
        fmt = "%(levelname)s: %(stage)s - %(msg)s"
        self.logger = logging.getLogger("pocketmapper")
        logging.basicConfig(level=log_level, format=fmt)

        # Writing to file
        formatter = logging.Formatter(fmt)
        fh = logging.FileHandler("test.log")
        fh.setLevel(log_level)
        fh.setFormatter(formatter)
        self.logger.addHandler(fh)

        self.logger.debug("Level set to DEBUG", extra=self._stage)

    def _configure(self, **kwargs):
        """
        Configures settings for the PocketMapper workflow.
        """
        self._stage.update({"stage": "Configuring Settings"})

        # Unrecognised arguments
        if len(self._uncaught_args) > 0:
            logging.critical(f"Unrecognised args: {list(self._uncaught_args.keys())}", extra=self._stage)
            exit(1)

        # Populates settings from the settings file if provided
        if self._settings_file is not None:
            # Checking that the supplied settings file exists
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
        for key, value in kwargs.items():
            if value is not None:
                self._settings[key] = value

        # Default settings
        cache_dir = self._settings.get("cache_dir", "pocketmapper_cache")
        now = datetime.now().strftime("%y%m%d_%H%M%S")
        results_dir = self._settings.get("results_dir", f"pocketmapper_results_{now}")
        defaults = {
            "cache_dir": cache_dir,
            "structure_dir": os.path.join(cache_dir, "pdb_structures"),
            "pocket_dir": os.path.join(cache_dir, "pockets"),
            "foldseek_tmp_dir": os.path.join(cache_dir, "foldseek_tmp"),
            "pisa_dir": os.path.join(cache_dir, "pisa_pockets"),
            "divided_struct_dir": os.path.join(cache_dir, "divided_structs"),
            "results_dir": results_dir,
            "query_dir": os.path.join(results_dir, "query_structures"),
            "target_dir": os.path.join(results_dir, "target_structures"),
            "foldseek_path": os.path.join(results_dir, "foldseek_results.tsv"),
            "pocket_comparison_path": os.path.join(results_dir, "pocket_comparison.tsv"),
            "foldseek": True,
            "pisa_pockets": True,
            "structure": False,
        }
        for key, value in defaults.items():
            if key not in self._settings:
                self._settings[key] = value

        logging.debug(f"\n{self._settings}", extra=self._stage)

    def _determine_query_target_types(self):
        self._stage.update({"stage": "Determine Query/Target Types"})

        # possible types for both query and target:
        #   pdb_chain_chain
        #   file (pdb_chain_chain per line)
        #   alphafold3 structure (implement later)

        # possible types for target only:
        #   foldseek database
        #

        self._query_type = None
        self._target_type = None

        query = self._settings.get("query")
        if os.path.isfile(query):
            self._query_type = "file"
        else:
            self._query_type = "pdb_chain_chain"

        target = self._settings.get("target")
        if os.path.isfile(target):
            if ".txt" in target:
                self._target_type = "file"
        elif target.lower() == "ted":
            self._target_type = "foldseek_db"
            self._settings["target_dir"] = files("TED_fs").joinpath("TED")  # Overriding target dir to point to the db
        else:
            self._target_type = "pdb_chain_chain"

        self.logger.debug(f"Determined query type: {self._query_type}", extra=self._stage)
        self.logger.debug(f"Determined target type: {self._target_type}", extra=self._stage)

    def _check_query_target_set(self):
        # Checking a target and query is specified
        self._stage.update({"stage": "Input Validation"})
        if not self._settings.get("query"):
            raise ValueError("No query specified. Use --query")
        if not self._settings.get("target"):
            raise ValueError("No target specified. Use --target")

        # Checking single pdb inputs

        """input_re = re.compile(r"[A-Za-z0-9]{4}_[A-Za-z0-9]_[A-Za-z0-9]")
        for key in ["query", "target"]:
            value = self._settings.get(key)
            if value and not input_re.match(value):
                self.logger.critical(
                    f"{key.capitalize()} '{value}' does not match required format 'PDB_CHAIN_CHAIN'.", extra=self._stage
                )
                exit(1)"""

    def _read_query_target_data(self):
        self._stage.update({"stage": "Reading Input Data"})

        self._query_data = self._parse_inputs(self._query, self._query_type, type="query")
        if type(self._query_data) is pd.DataFrame:
            self.logger.debug(f"query_data:\n{self._query_data.head()}", extra=self._stage)
        else:
            self.logger.debug(f"query_data:\n{self._query_data}", extra=self._stage)

        self._target_data = self._parse_inputs(self._target, self._target_type, type="target")
        if type(self._target_data) is pd.DataFrame:
            self.logger.debug(f"target_data:\n{self._target_data.head()}", extra=self._stage)
        else:
            self.logger.debug(f"target_data:\n{self._target_data}", extra=self._stage)

        """
        query_df = self._make_tq_df("query", "query_file", type="query")
        target_df = self._make_tq_df("target", "target_file", type="target")
        """
        pdb_data = [self._query_data]
        if type(self._target_data) is pd.DataFrame:
            pdb_data.append(self._target_data)
        self._pdb_df = pd.concat(pdb_data, ignore_index=True)

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
            # For foldseek db, we just return the path to the db
            return value
        else:
            self.logger.critical(f"Unknown value_type '{value_type}'", extra=self._stage)
            exit(1)

    def pdb_info_to_df(self, pdbs, domains, motifs, **kwargs):
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
        self._stage.update({"stage": "Fetch Structures"})
        logging.info("Checking for mmCIF structures...", extra=self._stage)
        found_map = lib.get_mmcifs(
            pdb_list=self._pdb_df["interaction_pdb"].unique(),
            out_dir=self._settings["structure_dir"],
        )
        self._pdb_df["structure_found"] = self._pdb_df["interaction_pdb"].map(found_map)

        self._stage.update({"stage": "Verify Structures"})
        if (
            not self._pdb_df.query("structure_found and type == 'query'").empty
            and self._query_type in self._requires_structures
        ):
            logging.info("Query structures found.", extra=self._stage)
        else:
            logging.critical("No query structures found locally or via download", extra=self._stage)
            exit(1)

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

        """self._stage.update({"stage": "Verify Preprocessing"})
        if self._pdb_df.query("divided_struct and type == 'query'").empty:
            logging.critical("No query structure could be preprocessed", extra=self._stage)
            exit(1)
        if self._pdb_df.query("divided_struct and type == 'target'").empty:
            logging.critical("No target structure could be preprocessed", extra=self._stage)
            exit(1)"""

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

    def _run_foldseek(self):
        self._stage.update({"stage": "Foldseek Alignment"})
        logging.info("Running Foldseek easy-search...", extra=self._stage)
        cmd = [
            "foldseek",
            "easy-search",
            self._settings["query_dir"],
            self._settings["target_dir"],
            self._settings["foldseek_path"],
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
            "-v",
            "0",
        ]
        subprocess.run(cmd, check=True)

    def _compare_pockets_and_save(self, pockets):
        self._stage.update({"stage": "Pocket Comparison"})
        logging.info("Comparing pockets...", extra=self._stage)
        alignment_df = pd.read_csv(self._settings["foldseek_path"], sep="\t", engine="c")
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
            "foldseek_tmp_dir",
            "query_dir",
            # "target_dir",
        ]
        # TODO this is unsafe
        for dir in tmp_dirs:
            shutil.rmtree(self._settings[dir])


def main():
    fire.Fire(PocketMapper())


if __name__ == "__main__":
    main()
