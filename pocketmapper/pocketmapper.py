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


class PocketMapper:
    def __init__(
        self, job_file=None, log_file=None, verbose=False, debug=False
    ):
        # Setting up logger
        logging.getLogger()

        if debug:
            log_level = logging.DEBUG
        elif verbose:
            log_level = logging.INFO
        else:
            log_level = logging.WARNING

        if log_file:
            logging.basicConfig(filename=log_file, level=log_level)
        else:
            logging.basicConfig(level=log_level)  # Gets default streamHandler

    def log_test(self):
        logging.debug("debug works!")
        logging.info("info works!")
        logging.warning("warning works!")

    # TODO implement caching option
    def search(self, file, caching=False):
        """
        Requires a job file with entries
            queries
            structure_dir
        """

        # Reading the job file
        logging.info("Reading Job File")
        with open(file) as f:
            job_data = json.load(f)
        self.job_data = job_data
        logging.debug(self.job_data)

        # If any of the specified directories don't exist, make them
        for dir in ["structure_dir", "domain_dir", "motif_dir", "pocket_dir"]:
            if not os.path.exists(job_data[dir]):
                os.mkdir(job_data[dir])

        # Formatting queries and downloading missing PDB files
        query_df = pd.DataFrame.from_dict(job_data["queries"], orient="index")
        query_df["pdb_domain"] = query_df.apply(
            lambda x: x.interaction_pdb + "_" + x.domain_chain, axis=1
        )
        query_df["pdb_domain_motif"] = query_df.apply(
            lambda x: x.pdb_domain + "_" + x.motif_chain, axis=1
        )
        status = {"input": job_data["queries"]}

        # Downloading missing PDB files
        logging.info("Checking mmCIF structures...")
        status["structure_found"] = lib.get_mmcifs(
            pdb_list=query_df["interaction_pdb"].unique(),
            out_dir=job_data["structure_dir"],
        )
        # Updating query_df with a column indicating if the structure is available
        query_df["structure_found"] = query_df["interaction_pdb"].map(
            status["structure_found"]
        )

        # Dividing the structure files into
        logging.info("Dividing mmCIF structures...")
        status["divided_struct"] = lib.pdb_preprocessing(
            queries=query_df.query("structure_found")
            .iloc[:, :3]
            .itertuples(index=False),
            ref_dir=job_data["structure_dir"],
            domain_dir=job_data["domain_dir"],
            motif_dir=job_data["motif_dir"],
        )
        # Updating query_df with a column indicating if the minimal interraction structures are available
        query_df["divided_struct"] = (
            query_df["pdb_domain_motif"]
            .map(status["divided_struct"])
            .fillna(False)
        )

        # regex of all interactionPdb_domainChain string to match the files for foldseek comparison
        fs_file_regex = (
            "("
            + ")|(".join(
                query_df.query("divided_struct")["pdb_domain"].unique()
            )
            + ")"
        )

        # Running foldseek
        subprocess.run(
            [
                "foldseek",
                "easy-search",
                job_data["domain_dir"],  # query folder of structure
                job_data["domain_dir"],  # target folder of structures
                job_data["foldseek_path"],  # output file
                job_data["foldseek_tmp_path"],  # temp folder
                "--format-output",
                "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,lddt,qaln,taln,u,t",
                "--format-mode",  # BLAST with headers
                "4",
                "-e",  # e-value threshold
                "0.001",
                "--file-include",
                fs_file_regex,
                "--exhaustive-search",
            ]
        )

        # Retrieving/calculating pockets
        print("Getting Pockets")
        pockets, problem_atoms, problem_residues = lib.calculate_pockets(
            queries=query_df.query("divided_struct")
            .iloc[:, :3]
            .itertuples(index=False),
            motif_dir=job_data["motif_dir"],
            pocket_dir=job_data["pocket_dir"],
        )
        if len(problem_atoms) > 0:
            logging.warning(f"Atoms with no VdW radii: {problem_atoms}")
        if len(problem_residues) > 0:
            logging.warning(
                f"Residues with no single AA code: {problem_residues}"
            )

        print("Comparing Pockets")
        p_c_path = job_data["pocket_comparison_path"]
        if os.path.exists(p_c_path):
            pockets_df = pd.read_csv(p_c_path, sep="\t")
        else:
            blosum_path = os.path.join(
                os.path.dirname(__file__), "blosum62.bla"
            )
            alignment_df = pd.read_csv(
                job_data["foldseek_path"], sep="\t", engine="c"
            )
            pockets_df = lib.compare_pockets(
                alignment_df, pockets, blosum_path=blosum_path
            )  # , alphafold=ALPHAFOLD, alphafold_dir=ALPHAFOLD_DIR)
            pockets_df.to_csv(p_c_path, index=False, sep="\t")

        # Saving the query_df for reference
        query_df.to_csv(
            r"/Users/lellingboe/Work/data/pocketmapper/test/query.csv"
        )


def main():
    fire.Fire(PocketMapper)


if __name__ == "__main__":
    main()
