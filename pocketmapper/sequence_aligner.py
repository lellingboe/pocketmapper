"""
Code related to local sequence alignment
"""

from Bio import Align
from Bio.Align import substitution_matrices
import gemmi
from itertools import product
import pandas as pd


class SequenceAligner:
    def __init__(self):
        self.single_aa_code = {
            "CYS": "C",
            "ASP": "D",
            "SER": "S",
            "GLN": "Q",
            "LYS": "K",
            "ILE": "I",
            "PRO": "P",
            "THR": "T",
            "PHE": "F",
            "ASN": "N",
            "GLY": "G",
            "HIS": "H",
            "LEU": "L",
            "ARG": "R",
            "TRP": "W",
            "ALA": "A",
            "VAL": "V",
            "GLU": "E",
            "TYR": "Y",
            "MET": "M",
            "SEP": "S",  # phosphoserine
            "TPO": "T",  # phosphotheonine
            "PTR": "Y",  # phosphotyrosine
            "MSE": "M",  # selenomethionine
        }

    def _replaceNonCommonResidues(self, peptide):
        processed_peptide = list(peptide)
        common_aas = list("ACDEFGHIKLMNPQRSTVWY")

        for i in range(0, len(peptide)):
            if processed_peptide[i] not in common_aas:
                processed_peptide[i] = "X"

        return "".join(processed_peptide)

    def _align_seqs(self, peptide1, peptide2, aligner=None):

        if aligner is None:
            aligner = Align.PairwiseAligner()
            aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")

        peptide1 = self._replaceNonCommonResidues(peptide1)
        peptide2 = self._replaceNonCommonResidues(peptide2)
        alignments = aligner.align(peptide1, peptide2)

        peptide1_aligned = [peptide1[i] if i != -1 else "-" for i in alignments[0].indices[0]]
        peptide2_aligned = [peptide2[i] if i != -1 else "-" for i in alignments[0].indices[1]]

        return [peptide1_aligned, peptide2_aligned]

    def align_records(self, query_records, target_records):
        """
        TODO docstring
        specs for df?
        """

        # building a mapping of preprocess_name to sequence for all queries and targets (to avoid redundant structure parsing and sequence extraction during alignment) - this assumes that preprocess_name is unique across queries and targets, which should be the case if they are in the format "P12345_A" or "1ABC_A"
        # Reads the full reference structure rather than the per-chain files under
        # foldseek_preprocessed_structure_dir: that preprocessing only runs on the Foldseek branch, and
        # we select the chain ourselves below, so the split copies buy us nothing here.
        name_to_seq = {}
        for record in query_records + target_records:
            name = record["preprocess_name"]
            if name in name_to_seq:
                continue  # skip if already processed
            path = record["struct_path"]
            st = gemmi.read_structure(path)  # format inferred from the extension, so local .pdb inputs work too
            st.setup_entities()
            aln_chain = record["chain_info"][0]  # e.g., "A" from "A_B" or "A"
            seq = "".join(
                [self.single_aa_code.get(res.name, "X") for res in st[0][aln_chain].get_polymer() if "CA" in res]
            )
            name_to_seq[name] = seq

        # Performing pairwise alignment
        aligner = Align.PairwiseAligner()
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        result_rows = []
        for q_record, t_record in product(query_records, target_records):
            query = q_record["preprocess_name"]
            target = t_record["preprocess_name"]
            qseq = name_to_seq[query]
            tseq = name_to_seq[target]
            qaln, taln = self._align_seqs(qseq, tseq, aligner)
            aln_len = len(qaln)

            identity = 0
            mismatch = 0
            gapopen = 0
            a_prev = "X"
            b_prev = "X"
            for a, b in zip(qaln, taln):
                if a == b and a != "-":
                    identity += 1
                elif a != b and a != "-" and b != "-":
                    mismatch += 1
                if a == "-" and a_prev != "-":
                    gapopen += 1
                if b == "-" and b_prev != "-":
                    gapopen += 1
                a_prev = a
                b_prev = b
            qend = len([x for x in qaln if x != "-"])
            tend = len([x for x in taln if x != "-"])

            result = {
                "query": query,
                "target": target,
                "fident": identity / aln_len,
                "alnlen": aln_len,
                "mismatch": mismatch / aln_len,
                "gapopen": gapopen,
                "qstart": 1,
                "qend": qend,
                "tstart": 1,
                "tend": tend,
                "evalue": "-",
                "lddt": "-",
                "qaln": "".join(qaln),
                "taln": "".join(taln),
                "u": "-",
                "t": "-",
                "qseq": qseq,
                "tseq": tseq,
            }

            result_rows.append(result)
        return pd.DataFrame(result_rows)
