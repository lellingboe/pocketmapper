"""
Local pairwise sequence alignment, the fallback for when Foldseek is unavailable.

Biopython's PairwiseAligner over BLOSUM62 stands in for the structural aligner and produces the
same alignment table, so the rest of the pipeline cannot tell the two apart. It has no structural
information to offer, so it writes "-" for the `u` and `t` transform columns -- which is why
structural superposition by `align_struct_method="foldseek"` yields the query alone on this path.
"""

from Bio import Align
from Bio.Align import substitution_matrices
import gemmi
from itertools import product
import pandas as pd

from pocketmapper.constants import ALIGNMENT_COLUMNS, SINGLE_AA_CODE


class SequenceAligner:
    """
    Produces the alignment table from sequence alone, without Foldseek.
    """

    def _replaceNonCommonResidues(self, peptide):
        """
        Map anything outside the 20 standard amino acids to "X".

        BLOSUM62 is only defined over the standard alphabet, so an unmapped residue code would raise
        rather than score.

        Args:
            peptide (str): A single-letter sequence.

        Returns:
            str: The sequence with non-standard codes replaced by "X".
        """
        processed_peptide = list(peptide)
        common_aas = list("ACDEFGHIKLMNPQRSTVWY")

        for i in range(0, len(peptide)):
            if processed_peptide[i] not in common_aas:
                processed_peptide[i] = "X"

        return "".join(processed_peptide)

    def _align_seqs(self, peptide1, peptide2, aligner):
        """
        Align two sequences and expand the result to gapped strings.

        Args:
            peptide1 (str): Query sequence.
            peptide2 (str): Target sequence.
            aligner (Bio.Align.PairwiseAligner): Configured aligner; only the top alignment is used.

        Returns:
            list: [query_aligned, target_aligned], each a list of characters with "-" for gaps.
        """
        peptide1 = self._replaceNonCommonResidues(peptide1)
        peptide2 = self._replaceNonCommonResidues(peptide2)
        alignments = aligner.align(peptide1, peptide2)

        peptide1_aligned = [peptide1[i] if i != -1 else "-" for i in alignments[0].indices[0]]
        peptide2_aligned = [peptide2[i] if i != -1 else "-" for i in alignments[0].indices[1]]

        return [peptide1_aligned, peptide2_aligned]

    def align_records(self, query_records, target_records):
        """
        Align every query against every target and build the alignment table.

        Sequences are extracted once per `preprocess_name` and reused across pairs, so a chain appearing
        on both sides is parsed only once.

        The returned columns are pinned to `constants.ALIGNMENT_COLUMNS`. That order is a positional
        contract shared with `_foldseek_alignment` and with `pocket_comparison`, which unpacks each row
        into an `AlignmentRow` by position -- see the note above the constant.

        Args:
            query_records (list): QTRecord dicts for the query side.
            target_records (list): QTRecord dicts for the target side.

        Returns:
            pandas.DataFrame: One row per query/target pair, columns in `ALIGNMENT_COLUMNS` order.
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
            seq = "".join([SINGLE_AA_CODE.get(res.name, "X") for res in st[0][aln_chain].get_polymer() if "CA" in res])
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
        # Columns are pinned to the shared contract rather than left to dict order: pocket_comparison
        # unpacks these rows positionally, so the order matters as much as the names.
        return pd.DataFrame(result_rows, columns=ALIGNMENT_COLUMNS)
