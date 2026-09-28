"""
Local pairwise sequence alignment, the `--aligner seq` alternative to Foldseek.

Biopython's PairwiseAligner over BLOSUM62 stands in for the structural aligner and produces the
same alignment table. It has no structural information to offer, so it writes "-" for the `u` and
`t` transform columns.
"""

from itertools import product

import gemmi
import pandas as pd
from Bio import Align
from Bio.Align import substitution_matrices

from pocketmapper.constants import ALIGNMENT_COLUMNS
from pocketmapper.lib import one_letter_code
from pocketmapper.lib import split_chain_info


class SequenceAligner:
    """
    Produces the alignment table from sequence alone, without Foldseek.
    """

    def replaceNonCommonResidues(self, peptide):
        """
        Map anything outside the 20 standard amino acids to "X".

        Args:
            peptide (str): A single-letter sequence.

        Returns:
            str: The sequence with non-standard codes replaced by "X".
        """
        # BLOSUM62 would raise rather than score a code outside its alphabet
        processed_peptide = list(peptide)
        common_aas = list("ACDEFGHIKLMNPQRSTVWY")

        for i in range(0, len(peptide)):
            if processed_peptide[i] not in common_aas:
                processed_peptide[i] = "X"

        return "".join(processed_peptide)

    def align_seqs(self, peptide1, peptide2, aligner):
        """
        Align two sequences and expand the result to gapped strings.

        Args:
            peptide1 (str): Query sequence.
            peptide2 (str): Target sequence.
            aligner (Bio.Align.PairwiseAligner): Configured aligner; only the top alignment is used.

        Returns:
            list: [query_aligned, target_aligned], each a list of characters with "-" for gaps.
        """
        peptide1 = self.replaceNonCommonResidues(peptide1)
        peptide2 = self.replaceNonCommonResidues(peptide2)
        alignments = aligner.align(peptide1, peptide2)

        peptide1_aligned = [peptide1[i] if i != -1 else "-" for i in alignments[0].indices[0]]
        peptide2_aligned = [peptide2[i] if i != -1 else "-" for i in alignments[0].indices[1]]

        return [peptide1_aligned, peptide2_aligned]

    def align_records(self, query_records, target_records):
        """
        Align every query against every target and build the alignment table.

        Args:
            query_records (list): QTRecord dicts for the query side.
            target_records (list): QTRecord dicts for the target side.

        Returns:
            pandas.DataFrame: One row per query/target pair, columns in `ALIGNMENT_COLUMNS` order.
        """

        # One sequence per preprocess_name, so a chain on both sides is parsed once. Read from the full
        # reference structure, selecting the chain here.
        name_to_seq = {}
        for record in query_records + target_records:
            name = record["preprocess_name"]
            if name in name_to_seq:
                continue  # skip if already processed
            path = record["struct_path"]
            st = gemmi.read_structure(path)  # format inferred from the extension, so local .pdb inputs work too
            st.setup_entities()
            aln_chain, _ = split_chain_info(record["chain_info"])
            seq = "".join([one_letter_code(res.name) for res in st[0][aln_chain].get_polymer() if "CA" in res])
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
            qaln, taln = self.align_seqs(qseq, tseq, aligner)
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
        # Pinned to the shared column order rather than left to dict order
        return pd.DataFrame(result_rows, columns=ALIGNMENT_COLUMNS)
