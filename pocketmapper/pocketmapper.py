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


class PocketMapper:
    def __init__(self):
        self._settings = {}
        self._stage = {"stage": "init"}

    # TODO implement caching option
    def search(
        self,
        query=None,
        target=None,
        query_file=None,
        target_file=None,
        verbose=False,
        debug=False,
        settings=None,
        cache_dir=None,
        results_dir=None,
    ):
        """
        Main orchestration method to run the pocket mapping workflow.
        """
        try:
            self._setup_logging(debug, verbose)
            self._configure(
                settings,
                cache_dir=cache_dir,
                results_dir=results_dir,
                query=query,
                target=target,
                query_file=query_file,
                target_file=target_file,
            )
            self._validate_inputs()
            self._prepare_directories()

            # all_df tracks structures and failures through the workflow
            all_df = self._prepare_dataframes()
            all_df = self._fetch_and_verify_structures(all_df)
            all_df = self._preprocess_structures(all_df)

            self._run_foldseek()

            pockets = self._calculate_and_retrieve_pockets(all_df)
            self._compare_pockets_and_save(pockets)

            logging.info("PocketMapper search completed successfully.", extra={"stage": "End"})

        except Exception as e:
            logging.exception(str(e), extra=self._stage)
            exit(1)

    def _setup_logging(self, debug, verbose):
        self._stage.update({"stage": "Logging Setup"})
        log_level = logging.WARNING
        if debug:
            log_level = logging.DEBUG
        elif verbose:
            log_level = logging.INFO
        fmt = "%(levelname)s: %(stage)s - %(msg)s"
        logging.basicConfig(level=log_level, format=fmt, force=True)

    def _configure(self, settings_file, cache_dir, results_dir, **kwargs):
        # Defult settings
        if cache_dir is None:
            cache_dir = "pocketmapper_cache"
        if results_dir is None:
            now = datetime.now().strftime("%y%m%d_%H%M%S")
            results_dir = f"pocketmapper_results_{now}"
        self._settings.update(
            {
                "structure_dir": os.path.join(cache_dir, "pdb_structures"),
                "pocket_dir": os.path.join(cache_dir, "pockets"),
                "foldseek_tmp_dir": os.path.join(cache_dir, "foldseek_tmp"),
                "query_dir": os.path.join(results_dir, "query_structures"),
                "target_dir": os.path.join(results_dir, "target_structures"),
                "foldseek_path": os.path.join(results_dir, "foldseek_results.tsv"),
                "pocket_comparison_path": os.path.join(results_dir, "pocket_comparison.tsv"),
                "foldseek": True,
                "pisa_pockets": True,
                "pisa_dir": "/Users/lellingboe/Work/data/PISA/interfaces",
                "structure": False,
            }
        )

        # Override defaults with settings file if provided
        self._stage.update({"stage": "Configuration"})
        if settings_file:
            self._read_settings(settings_file)

        # Override settings with any provided command-line arguments
        for key, value in kwargs.items():
            if value is not None:
                self._settings[key] = value

    def _read_settings(self, settings):
        try:
            with open(settings) as f:
                job_data = json.load(f)
        except FileNotFoundError:
            logging.critical("No settings file found at specified location", extra=self._stage)
            exit()
        except Exception:
            logging.exception("Error reading settings file", extra=self._stage)
            exit()
        self._settings.update(job_data)
        logging.info(job_data, extra=self._stage)

    def _validate_inputs(self):
        self._stage.update({"stage": "Input Validation"})
        if not self._settings.get("query") and not self._settings.get("query_file"):
            raise ValueError("No query specified. Use --query or --query_file.")
        if not self._settings.get("target") and not self._settings.get("target_file"):
            raise ValueError("No target specified. Use --target or --target_file.")

        input_re = re.compile(r"[A-Za-z0-9]{4}_[A-Za-z0-9]_[A-Za-z0-9]")
        for key in ["query", "target"]:
            value = self._settings.get(key)
            if value and not input_re.match(value):
                raise ValueError(f"{key.capitalize()} '{value}' does not match required format 'PDB_CHAIN_CHAIN'.")

    def _prepare_directories(self):
        self._stage.update({"stage": "Directory Preparation"})
        dirs_to_create = ["structure_dir", "query_dir", "target_dir", "pocket_dir"]
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
        return pd.concat([query_df, target_df], ignore_index=True)

    def _fetch_and_verify_structures(self, df):
        self._stage.update({"stage": "Fetch Structures"})
        logging.info("Checking for mmCIF structures...", extra=self._stage)
        found_map = lib.get_mmcifs(
            pdb_list=df["interaction_pdb"].unique(),
            out_dir=self._settings["structure_dir"],
        )
        df["structure_found"] = df["interaction_pdb"].map(found_map)

        self._stage.update({"stage": "Verify Structures"})
        if not df.query("structure_found and type == 'query'").empty:
            logging.info("Query structures found.", extra=self._stage)
        else:
            raise FileNotFoundError("No query structures found locally or via download.")

        if not df.query("structure_found and type == 'target'").empty:
            logging.info("Target structures found.", extra=self._stage)
        else:
            raise FileNotFoundError("No target structures found locally or via download.")
        return df

    def _preprocess_structures(self, df):
        self._stage.update({"stage": "Preprocess Structures"})
        logging.info("Dividing mmCIF structures...", extra=self._stage)
        divided_map = lib.pdb_preprocessing_gemmi(
            df=df.query("structure_found"),
            ref_dir=self._settings["structure_dir"],
            query_dir=self._settings["query_dir"],
            target_dir=self._settings["target_dir"],
        )
        df["divided_struct"] = df["pdb_domain_motif"].map(divided_map).fillna(False)

        self._stage.update({"stage": "Verify Preprocessing"})
        if df.query("divided_struct and type == 'query'").empty:
            raise ValueError("Failed to preprocess any query structures.")
        if df.query("divided_struct and type == 'target'").empty:
            raise ValueError("Failed to preprocess any target structures.")
        return df

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
            r"[0-9A-Z]{4}_[0-9A-Za-z]\.cif",
            "--exhaustive-search",
            "-v",
            "2",
        ]
        subprocess.run(cmd, check=True)

    def _calculate_and_retrieve_pockets(self, df):
        self._stage.update({"stage": "Pocket Calculation"})
        logging.info("Calculating/retrieving pockets...", extra=self._stage)
        pockets, problem_atoms, problem_residues = lib.calculate_pockets(
            df=df.query("divided_struct"),
            target_dir=self._settings["target_dir"],
            query_dir=self._settings["query_dir"],
            pocket_dir=self._settings["pocket_dir"],
        )
        if self._settings["pisa_pockets"]:
            logging.info("Retrieving PISA pockets...", extra=self._stage)
            pisa_pockets = lib.get_pisa_pocket(
                df=df.query("divided_struct"),
                pocket_dir=self._settings["pisa_dir"],
            )
            print(pisa_pockets.keys())

        if problem_atoms:
            logging.warning(f"Atoms with no VdW radii: {problem_atoms}")
        if problem_residues:
            logging.warning(f"Residues with no single AA code: {problem_residues}")
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


def main():
    fire.Fire(PocketMapper)


if __name__ == "__main__":
    main()
