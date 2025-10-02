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


class PocketMapper:
    def __init__(self):
        self._settings = {}
        self._stage = {"stage": "init"}
        self._all_df = None

    # TODO implement caching option
    def search(
        self,
        query=None,  # settings passed to configure
        target=None,
        query_file=None,
        target_file=None,
        settings=None,
        cache_dir=None,
        results_dir=None,
        verbose=False,  # settings passed to logging
        debug=False,
        help=None,  # help option
        **kwargs,
    ):
        """
        Main orchestration method to run the pocket mapping workflow.
        """

        try:
            # Setting up things
            self._help(help)
            self._setup_logging(debug, verbose)
            self._configure(  # configures the settings which have already been read
                settings_file=settings,
                cache_dir=cache_dir,
                results_dir=results_dir,
                query=query,
                target=target,
                query_file=query_file,
                target_file=target_file,
                uncaught_args=kwargs,
            )
            self._validate_inputs()
            self._prepare_directories()

            # Preparing structures for later
            self._prepare_dataframes()
            self._fetch_and_verify_structures()
            self._divide_structures()

            pockets = self._retrieve_pockets()
            pockets = self._get_atom_coords_from_cif(pockets)

            self._run_foldseek()

            self._compare_pockets_and_save(pockets)
            self._delete_tmp()

            logging.info("PocketMapper search completed successfully.", extra={"stage": "End"})

        # Unhandled exception stops the process and logs the error
        except Exception as e:
            logging.exception(str(e), extra=self._stage)
            exit(1)

    def _help(self, help):
        if help:
            print("Displaying help!")
            exit()

    def _setup_logging(self, debug, verbose):
        self._stage.update({"stage": "Logging Setup"})

        if debug:
            log_level = logging.DEBUG
        elif verbose:
            log_level = logging.INFO
        else:
            log_level = logging.WARNING
        fmt = "%(levelname)s: %(stage)s - %(msg)s"
        logging.basicConfig(level=log_level, format=fmt, force=True)

    def _configure(self, settings_file, uncaught_args, **kwargs):
        self._stage.update({"stage": "Configuring Settings"})

        # Unrecognised arguments
        if len(uncaught_args) > 0:
            logging.critical(f"Unrecognised args: {list(uncaught_args.keys())}", extra=self._stage)
            exit(1)

        # get all info from settings
        if settings_file:
            try:
                with open(settings_file) as f:
                    job_data = json.load(f)
            except FileNotFoundError:
                logging.critical("No settings file found at specified location", extra=self._stage)
                exit(1)
            except Exception:
                logging.exception("Unexpected error reading settings file", extra=self._stage)
                exit(1)
            finally:
                self._settings.update(job_data)

        # Override settings_file with any provided command-line arguments
        for key, value in kwargs.items():
            if value is not None:
                self._settings[key] = value

        # Defult settings
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

    def _validate_inputs(self):
        # Checking a target and query is specified
        self._stage.update({"stage": "Input Validation"})
        if not self._settings.get("query") and not self._settings.get("query_file"):
            raise ValueError("No query specified. Use --query or --query_file.")
        if not self._settings.get("target") and not self._settings.get("target_file"):
            raise ValueError("No target specified. Use --target or --target_file.")

        # Checking single pdb inputs
        input_re = re.compile(r"[A-Za-z0-9]{4}_[A-Za-z0-9]_[A-Za-z0-9]")
        for key in ["query", "target"]:
            value = self._settings.get(key)
            if value and not input_re.match(value):
                raise ValueError(f"{key.capitalize()} '{value}' does not match required format 'PDB_CHAIN_CHAIN'.")

    def _prepare_directories(self):
        self._stage.update({"stage": "Directory Preparation"})
        dirs_to_create = [
            "structure_dir",
            "query_dir",
            "target_dir",
            "pocket_dir",
            "pisa_dir",
            "divided_struct_dir",
            "results_dir",
        ]
        for dir_key in dirs_to_create:
            path = self._settings[dir_key]
            try:
                os.makedirs(path, exist_ok=True)
            except OSError as e:
                raise OSError(f"Error creating directory {path}: {e}")

    def _prepare_dataframes(self):
        self._stage.update({"stage": "Data Preparation"})
        query_df = self._make_tq_df("query", "query_file", type="query")
        target_df = self._make_tq_df("target", "target_file", type="target")
        self._all_df = pd.concat([query_df, target_df], ignore_index=True)

    def _fetch_and_verify_structures(self):
        self._stage.update({"stage": "Fetch Structures"})
        logging.info("Checking for mmCIF structures...", extra=self._stage)
        found_map = lib.get_mmcifs(
            pdb_list=self._all_df["interaction_pdb"].unique(),
            out_dir=self._settings["structure_dir"],
        )
        self._all_df["structure_found"] = self._all_df["interaction_pdb"].map(found_map)

        self._stage.update({"stage": "Verify Structures"})
        if not self._all_df.query("structure_found and type == 'query'").empty:
            logging.info("Query structures found.", extra=self._stage)
        else:
            logging.critical("No query structures found locally or via download", extra=self._stage)
            exit(1)

        if not self._all_df.query("structure_found and type == 'target'").empty:
            logging.info("Target structures found.", extra=self._stage)
        else:
            logging.critical("No target structures found locally or via download", extra=self._stage)
            exit(1)

    def _divide_structures(self):
        self._stage.update({"stage": "Preprocess Structures"})
        logging.info("Dividing mmCIF structures...", extra=self._stage)
        divided_map = lib.pdb_preprocessing_gemmi(
            df=self._all_df.query("structure_found"),
            ref_dir=self._settings["structure_dir"],
            cache_dir=self._settings["divided_struct_dir"],
            query_dir=self._settings["query_dir"],
            target_dir=self._settings["target_dir"],
        )
        self._all_df["divided_struct"] = self._all_df["pdb_domain_motif"].map(divided_map).fillna(False)

        self._stage.update({"stage": "Verify Preprocessing"})
        if self._all_df.query("divided_struct and type == 'query'").empty:
            logging.critical("No query structure could be preprocessed", extra=self._stage)
            exit(1)
        if self._all_df.query("divided_struct and type == 'target'").empty:
            logging.critical("No target structure could be preprocessed", extra=self._stage)
            exit(1)

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
            "--exhaustive-search",
        ]
        subprocess.run(cmd, check=True)

    def _retrieve_pockets(self):
        self._stage.update({"stage": "Pocket Calculation"})

        # WRITING PISA POCKETS
        logging.info("Retrieving PISA pockets...", extra=self._stage)
        downloader = pisa.PisaDownloader()
        downloader.get_interfaces(
            pdb_list=self._all_df.query("divided_struct")["interaction_pdb"].str.lower().unique(),
            summary_dir=os.path.join(self._settings["pisa_dir"], "summaries"),
            asm_dir=os.path.join(self._settings["pisa_dir"], "assemblies"),
            interface_dir=os.path.join(self._settings["pisa_dir"], "interfaces"),
        )
        pisa_pockets = lib.get_pisa_pockets(
            df=self._all_df.query("divided_struct"),
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
                for res in domain_chain:
                    res_id = str(res.seqid.num)
                    if res_id in pocket_keys:
                        pockets[pocket_id][res_id]["ca_coords"] = list(res.get_ca().pos)
                pockets[pocket_id]["has_coords"] = True
            except Exception:
                logging.warning(f"Error is getting coords for {pocket_id}", extra=self._stage)
                pockets[pocket_id]["has_coords"] = False

        with open(os.path.join(self._settings["pisa_dir"], "all_pockets_2.json"), "w") as f:
            json.dump(pockets, f)
        return pockets

    def _compare_pockets_and_save(self, pockets):
        self._stage.update({"stage": "Pocket Comparison"})
        logging.info("Comparing pockets...", extra=self._stage)
        alignment_df = pd.read_csv(self._settings["foldseek_path"], sep="\t", engine="c")
        blosum_path = os.path.join(os.path.dirname(__file__), "blosum62.bla")

        pockets_df = lib.compare_pockets(alignment_df, pockets, blosum_path=blosum_path)

        output_path = self._settings["pocket_comparison_path"]
        pockets_df.to_csv(output_path, index=False, sep="\t")
        logging.info(f"Pocket comparison results saved to {output_path}", extra=self._stage)

    def _make_tq_df(self, single, file, **kwargs):
        if self._settings.get(single) is not None:
            pdb, domain, motif = self._settings.get(single).split("_")
            try:
                df = pd.DataFrame.from_dict(
                    {
                        0: {
                            "interaction_pdb": pdb.upper(),
                            "domain_chain": domain,
                            "motif_chain": motif,
                        }
                    },
                    orient="index",
                )
            except Exception:
                logging.critical(f"Error with parsing {self._settings.get(single)}", extra=self._stage)
                exit()
        else:
            df = pd.read_csv(self._settings[file], sep="\t", index_col=False)
            df.interaction_pdb = df.interaction_pdb.str.upper()
            logging.debug("\n" + str(df.head(5)), extra={"stage": f"{file}"})

        df["pdb_domain"] = df.apply(lambda x: x.interaction_pdb + "_" + x.domain_chain, axis=1)
        df["pdb_domain_motif"] = df.apply(lambda x: x.pdb_domain + "_" + x.motif_chain, axis=1)
        for k, v in kwargs.items():
            df[k] = v
        return df

    def _delete_tmp(self):
        tmp_dirs = [
            "foldseek_tmp_dir",
            "query_dir",
            "target_dir",
        ]
        for dir in tmp_dirs:
            shutil.rmtree(self._settings[dir])


def main():
    fire.Fire(PocketMapper())


if __name__ == "__main__":
    main()
