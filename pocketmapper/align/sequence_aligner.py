"""
Local pairwise sequence alignment, the `--aligner seq` alternative to Foldseek.

Biopython's PairwiseAligner, in local mode over BLOSUM62 with affine gap scores, stands in for the
structural aligner and produces the same alignment table. As in Foldseek's table, `qaln`/`taln` hold only
the aligned region, `qstart`..`qend` and `tstart`..`tend` locate it (1-based, inclusive), and `mismatch` is a
count. It has no structural information to offer, so it writes "-" for the `evalue`, `lddt`, `u` and `t`
columns.
"""

import logging
from itertools import product

import pandas as pd
from Bio import Align
from Bio.Align import substitution_matrices

from pocketmapper.lib.constants import ALIGNMENT_COLUMNS
from pocketmapper.lib.exceptions import PocketMapperError
from pocketmapper.lib.residues import ca_residues, one_letter_code, split_chain_info
from pocketmapper.lib.structure_manager import load_structure

logger = logging.getLogger(__name__)

# Gap scores in Biopython's convention: a gap of length n scores OPEN_GAP_SCORE + (n - 1) * EXTEND_GAP_SCORE
OPEN_GAP_SCORE = -11
EXTEND_GAP_SCORE = -1


def build_aligner():
    """
    Build the local BLOSUM62 aligner with affine gap scores.

    Mode and gap scores are set explicitly, because Biopython's defaults (global mode, linear gap score -1) vary
    across releases and would open two gaps rather than accept a poor mismatch.

    Returns:
        Bio.Align.PairwiseAligner: The configured aligner.
    """
    return Align.PairwiseAligner(
        mode="local",
        substitution_matrix=substitution_matrices.load("BLOSUM62"),
        open_gap_score=OPEN_GAP_SCORE,
        extend_gap_score=EXTEND_GAP_SCORE,
    )


class SequenceAligner:
    """
    Produces the alignment table from sequence alone, without Foldseek.
    """

    def __init__(self, manager):
        """
        Keep the structure locator.

        Args:
            manager (StructureManager): Locates each record's reference structure.
        """
        self.manager = manager

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
        Align two sequences and summarise the top alignment in Foldseek's alignment columns.

        Args:
            peptide1 (str): Query sequence.
            peptide2 (str): Target sequence.
            aligner (Bio.Align.PairwiseAligner): A local-mode aligner, as from `build_aligner`; only the top
                alignment is used.

        Returns:
            dict or None: `fident`, `alnlen`, `mismatch`, `gapopen`, `qstart`, `qend`, `tstart`, `tend`,
                `qaln` and `taln`, or None when no residue pair scores above zero, so there is no local
                alignment.
        """
        peptide1 = self.replaceNonCommonResidues(peptide1)
        peptide2 = self.replaceNonCommonResidues(peptide2)
        # Empty when no residue pair scores above zero
        alignment = next(iter(aligner.align(peptide1, peptide2)), None)
        if alignment is None:
            return None

        # In local mode the indices cover only the aligned region; -1 marks a gap
        qaln = "".join(peptide1[i] if i != -1 else "-" for i in alignment.indices[0])
        taln = "".join(peptide2[i] if i != -1 else "-" for i in alignment.indices[1])
        # The first and last coordinates are the region's 0-based start and exclusive end on each sequence,
        # which are Foldseek's 1-based start and inclusive end once the start is shifted by one
        (qstart, tstart), (qend, tend) = alignment.coordinates[:, 0] + 1, alignment.coordinates[:, -1]

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

        return {
            "fident": identity / len(qaln),
            "alnlen": len(qaln),
            "mismatch": mismatch,
            "gapopen": gapopen,
            "qstart": int(qstart),
            "qend": int(qend),
            "tstart": int(tstart),
            "tend": int(tend),
            "qaln": qaln,
            "taln": taln,
        }

    def align_records(self, query_records, target_records):
        """
        Align every query against every target and build the alignment table.

        Args:
            query_records (list): QTRecord dicts for the query side.
            target_records (list): QTRecord dicts for the target side.

        Returns:
            pandas.DataFrame: One row per query/target pair, columns in `ALIGNMENT_COLUMNS` order.

        Raises:
            PocketMapperError: If a record's structure file is missing or cannot be read.
        """
        log_extra = {"stage": "Local Alignment"}

        # One sequence per preprocess_name, so a chain on both sides is parsed once. Read from the full
        # reference structure, selecting the chain here.
        name_to_seq = {}
        for record in query_records + target_records:
            name = record["preprocess_name"]
            if name in name_to_seq:
                continue  # skip if already processed
            path = self.manager.record_path(record)
            # Format inferred from the extension, so local .pdb inputs work too; ca_residues reads the
            # entity types load_structure sets
            st = load_structure(path, log_extra)
            if st is None:
                msg = f"Cannot align {record['pocket_id']}: structure file {path} is missing or unreadable"
                logger.critical(msg, extra=log_extra)
                raise PocketMapperError(msg)
            aln_chain, _ = split_chain_info(record["chain_info"])
            # The same walk as a pocket's seq_pos, so the alignment positions match it
            name_to_seq[name] = "".join(one_letter_code(res.name) for res in ca_residues(st[0][aln_chain]))

        aligner = build_aligner()
        result_rows = []
        for q_record, t_record in product(query_records, target_records):
            query = q_record["preprocess_name"]
            target = t_record["preprocess_name"]
            qseq = name_to_seq[query]
            tseq = name_to_seq[target]
            columns = self.align_seqs(qseq, tseq, aligner)
            # No row, as Foldseek gives none for a pair it cannot align
            if columns is None:
                logger.warning(
                    f"No local alignment between {query} and {target}: no residue pair scores above zero",
                    extra=log_extra,
                )
                continue

            result_rows.append(
                {
                    "query": query,
                    "target": target,
                    **columns,
                    "evalue": "-",
                    "lddt": "-",
                    "u": "-",
                    "t": "-",
                    "qseq": qseq,
                    "tseq": tseq,
                }
            )
        # Pinned to the shared column order rather than left to dict order
        return pd.DataFrame(result_rows, columns=ALIGNMENT_COLUMNS)
