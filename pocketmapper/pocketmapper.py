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


class PocketMapper:
    def __init__(self):
        # Defaults
        self._settings = {
            "structure_dir": "pdb",
            "query_dir": "query",
            "target_dir": "target",
            "pocket_dir": ".",
            "foldseek_path": "foldseek_alignment.tab",
            "foldseek_tmp_dir": "foldseek_temp",
            "pocket_comparison_path": "pocket_comparison.tab",
            "foldseek": True,
            "structure": False,
        }
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
        structure=None,
        foldseek=None,
        af2_target=None,
        **kwargs,
    ):

        # CREATING LOGGER
        logging.getLogger()
        if debug:
            log_level = logging.DEBUG
        elif verbose:
            log_level = logging.INFO
        else:
            log_level = logging.WARNING
        fmt = "%(levelname)s: %(stage)s - %(msg)s"
        logging.basicConfig(level=log_level, format=fmt)

        # VALIDATING COMMAND LINE OPTIONS

        # Checking for unrecognized command line options
        self._stage.update({"stage": "Checking Commandline Options"})
        if kwargs:
            msg = f"Unrecognised inputs: {",".join(kwargs.keys())}"
            logging.critical(msg, extra=self._stage)
            exit()

        # Reading the options file
        self._stage.update({"stage": "Reading Settings File"})
        if settings is not None:
            self._reading_options(settings=settings)

        # If value is specified by the command line use it
        if query is not None:
            self._settings["query"] = query
        if query_file is not None:
            self._settings["query_file"] = query_file
        if target is not None:
            self._settings["target"] = target
        if target_file is not None:
            self._settings["target_file"] = target_file
        if structure is not None:
            self._settings["structure"] = structure
        if foldseek is not None:
            self._settings["foldseek"] = foldseek

        # Checking prescence of query and target
        self._stage.update({"stage": "checking inputs"})
        if self._settings.get("query") is None and self._settings.get("query_file") is None:
            logging.critical(self._settings, extra=self._stage)
            logging.critical("No query specified", extra=self._stage)
            exit()
        if self._settings.get("target") is None and self._settings.get("target_file") is None:
            logging.critical("No target specified", extra=self._stage)
            exit()

        # Checking query/target format
        input_re = re.compile(r"[A-Za-z1-9]{4}_[A-Za-z1-9]_[A-Za-z1-9]")
        if self._settings.get("query") is not None:
            if not input_re.match(self._settings.get("query")):
                msg = f"Query '{self._settings.get("query")}' does not match required format"
                logging.critical(msg, extra=self._stage)
                exit()
        if self._settings.get("target") is not None:
            if not input_re.match(self._settings.get("target")):
                msg = f"Target '{self._settings.get("target")}' does not match required format"
                logging.critical(msg, extra=self._stage)
                exit()

        # If any of the specified directories don't exist, make them
        self._stage.update({"stage": "Checking/Creating Directories"})
        dirs = [
            "structure_dir",
            "query_dir",
            "target_dir",
            "pocket_dir",
        ]
        for dir in dirs:
            try:
                if not os.path.exists(self._settings[dir]):
                    os.mkdir(self._settings[dir])
            except Exception:
                logging.critical(
                    f"Error creating directory {self._settings[dir]}",
                    extra=self._stage,
                )
                exit()

        # Format queries
        self._stage.update({"stage": "Formatting Query"})
        query_df = self._make_tq_df("query", "query_file", type="query")
        target_df = self._make_tq_df("target", "target_file", type="target")
        all_df = pd.concat([query_df, target_df], ignore_index=True)
        status = {}

        # Downloading missing PDB files
        self._stage.update({"stage": "Checking mmCIF structures"})
        logging.info("", extra=self._stage)
        status["structure_found"] = lib.get_mmcifs(
            pdb_list=all_df["interaction_pdb"].unique(),
            out_dir=self._settings["structure_dir"],
        )
        logging.debug(status, extra={"stage": "Available MMCif Structures"})
        # Updating query_df with a column indicating if the structure is available
        all_df["structure_found"] = all_df["interaction_pdb"].map(status["structure_found"])
        logging.debug(all_df.head(), extra={"stage": "Download Results in DataFrame"})

        self._stage.update({"stage": "Verifying downloaded/cached structures"})
        # Ensuring both a target and query structure were found
        if len(all_df.query("structure_found & type == 'query'")) < 1:
            logging.critical("No query structures were found locally or able to be downloaded", extra=self._stage)
            exit()
        if len(all_df.query("structure_found & type == 'target'")) < 1:
            logging.critical("No target structures were found locally or able to be downloaded", extra=self._stage)
            exit()

        # Dividing the structure files into
        self._stage.update({"stage": "Dividing mmCIF structures"})
        logging.info("", extra=self._stage)
        status["divided_struct"] = lib.pdb_preprocessing(
            df=all_df.query("structure_found"),
            ref_dir=self._settings["structure_dir"],
            query_dir=self._settings["query_dir"],
            target_dir=self._settings["target_dir"],
        )

        # Updating query_df with a column indicating if the minimal interraction structures are available
        all_df["divided_struct"] = all_df["pdb_domain_motif"].map(status["divided_struct"]).fillna(False)

        # Ensuring both a target and query structure were found
        self._stage.update({"stage": "Verifying divided structures"})
        if len(all_df.query("divided_struct & type == 'query'")) < 1:
            logging.critical("No query structures are able to be divided", extra=self._stage)
            exit()
        if len(all_df.query("divided_struct & type == 'target'")) < 1:
            logging.critical("No target structures  are able to be divided", extra=self._stage)
            exit()

        # Running foldseek
        subprocess.run(
            [
                "foldseek",
                "easy-search",
                self._settings["query_dir"],  # query folder of structure
                self._settings["target_dir"],  # target folder of structures
                self._settings["foldseek_path"],  # output file
                self._settings["foldseek_tmp_dir"],  # temp folder
                "--format-output",
                "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,lddt,qaln,taln,u,t",
                "--format-mode",  # BLAST with headers
                "4",
                "-e",  # e-value threshold
                "0.001",
                "--file-include",
                r"[0-9A-Z]{4}_[0-9A-Za-z]\.cif",
                "--exhaustive-search",
            ]
        )

        # Retrieving/calculating pockets
        self._stage.update({"stage": "Calculating/retrieving Pockets"})
        logging.info("", extra=self._stage)
        pockets, problem_atoms, problem_residues = lib.calculate_pockets(
            df=all_df.query("divided_struct"),
            target_dir=self._settings["target_dir"],
            query_dir=self._settings["query_dir"],
            pocket_dir=self._settings["pocket_dir"],
        )
        if len(problem_atoms) > 0:
            logging.warning(f"Atoms with no VdW radii: {problem_atoms}")
        if len(problem_residues) > 0:
            logging.warning(
                f"Residues with no single AA code: {problem_residues}",
                extra=self._stage,
            )

        print("Comparing Pockets")
        p_c_path = self._settings["pocket_comparison_path"]
        blosum_path = os.path.join(os.path.dirname(__file__), "blosum62.bla")
        alignment_df = pd.read_csv(self._settings["foldseek_path"], sep="\t", engine="c")
        pockets_df = lib.compare_pockets(
            alignment_df, pockets, blosum_path=blosum_path
        )  # , alphafold=ALPHAFOLD, alphafold_dir=ALPHAFOLD_DIR)
        pockets_df.to_csv(p_c_path, index=False, sep="\t")

        # Saving the query_df for reference
        query_df.to_csv(r"/Users/lellingboe/Work/data/pocketmapper/test/query.csv")

        ##############################
        logging.critical("Success!", extra={"stage": "End"})
        exit()

    def _reading_options(self, settings):
        try:
            with open(settings) as f:
                job_data = json.load(f)
        except Exception:
            logging.critical("Could not read options file", extra=self._stage)
            exit()
        self._settings.update(job_data)
        logging.info(job_data, extra=self._stage)

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
