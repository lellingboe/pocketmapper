from Bio import Align
from Bio.Align import substitution_matrices
import os
import gemmi
from itertools import combinations
import pandas as pd


class LocalAligner:
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
            "SEP": "S",  # phospho
            "TPO": "T",  # phospho
            "PTR": "Y",  # phospho
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

    def align_df(self, df, struct_dir):

        cols = "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,lddt,qaln,taln,u,t".split(
            ","
        )

        aligner = Align.PairwiseAligner()
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        df = df.query("divided_struct")
        name_to_seq = {}
        for index, row in df.iterrows():
            path = os.path.join(struct_dir, f"{row['pdb_domain']}.cif.gz")
            st = gemmi.read_structure(path, format=gemmi.CoorFormat.Mmcif)
            st.setup_entities()
            seq = "".join([self.single_aa_code.get(x.name, "X") for x in st[0][row["domain_chain"]].get_polymer()])
            name_to_seq[row["pdb_domain"]] = seq

        result_rows = []
        # TODO make this query and target specific
        # queries = df.query("type == 'query' & divided_struct")
        # targets = df.query("type == 'target' & divided_struct")
        for q, t in combinations(name_to_seq.keys(), 2):
            seq1 = name_to_seq[q]
            seq2 = name_to_seq[t]
            qaln, taln = self._align_seqs(seq1, seq2, aligner)
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
                "query": q,
                "target": t,
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
            }

            result_rows.append(result)
            # Here you would typically store or process the alignment results
            # For example, you could create a new DataFrame with the results

        return pd.DataFrame(result_rows, columns=cols)
