"""
Building of aligned multi-structure PDB files.

Takes one query structure plus one or more matched targets and writes them back out as a single
human-readable PDB, as glue around gemmi: read mmCIF, apply per-target rotations and translations,
rename chains so each structure stays distinguishable, and emit a combined PDB with COMPND
metadata naming each model.
"""

import logging
import string
from itertools import count, permutations
import gemmi
import numpy as np


class StructureAligner:
    """
    Superposes structures onto a reference and writes them as one PDB.

    Two entry points, matching the two transform sources the `align_struct_method` setting selects:
    `transform` applies transforms it is handed ("pocket"), and `foldseek_transform` sources them from
    a Foldseek alignment table first ("foldseek").
    """

    def __init__(self):
        """
        Initialise the logger. The aligner holds no other state.
        """
        self.logger = logging.getLogger(__name__)
        self._log_extra = {"stage": "StructureAligner"}
        logging.debug("Initialized", extra=self._log_extra)

    def _char_gen(self):
        """
        Yield short, PDB-friendly chain identifiers for aligned output.

        The output structure may hold more chains than any one input, so names are generated rather than
        reused. "0" is skipped throughout: it is reserved for the domain chain, which
        `_apply_transformation` renames so it can be found consistently across models.

        Single characters are emitted first, then two-character pairs once that space is exhausted.

        Yields:
            str: A unique chain label suitable for writing into the final PDB.
        """
        nice_chars = string.digits + string.ascii_letters
        for x in nice_chars[1:]:  # 0 is reserved for domain names
            yield x
        for x, y in permutations(nice_chars, 2):
            yield (x + y)

    def _apply_transformation(self, structs, domain_chains, motif_chains, us, ts):
        """
        Apply rigid-body transforms to a set of structures and merge them.

        The first structure is the reference frame. The domain chain of each input is renamed to "0" in
        the output so it can be recognised consistently; any motif chain gets a generated name.

        Args:
            structs (list): gemmi.Structures to align and merge.
            domain_chains (list): Domain chain id to preserve from each input structure.
            motif_chains (list): Motif chain id per input, or None where there is none.
            us (list): 3x3 rotation matrices, one per structure, in gemmi's LEFT-multiplying convention.
            ts (list): Translation vectors of length 3, one per structure.

        Returns:
            gemmi.Structure: A merged structure holding one model per input.
        """
        # Align everything to the first struct
        ref_st = gemmi.Structure()
        chain_names = self._char_gen()

        for i, cn, st, dc, mc, u, t in zip(count(1), chain_names, structs, domain_chains, motif_chains, us, ts):
            ref_st.add_model(gemmi.Model(i))

            # Apply transformation
            trans = gemmi.Transform(gemmi.Mat33(u), gemmi.Vec3(*t))
            st[0].transform_pos_and_adp(trans)

            # Add chains to the reference structure
            st[0][dc].name = "0"
            ref_st[-1].add_chain(st[0]["0"])
            if mc is not None:
                st[0][mc].name = cn
                ref_st[-1].add_chain(st[0][cn])

        ref_st.setup_entities()
        return ref_st

    def transform(self, aln_records, transforms, out_path):
        """
        Build an aligned multi-structure PDB from ready-made rigid-body transforms.

        The general entry point: it knows nothing about where a transform came from, only how to apply it.
        `foldseek_transform` sources them from Foldseek's whole-chain alignment; `_align_structs` sources
        them from the pocket superposition in pocket_comparison.tsv.

        The first record is the reference frame and is always placed untransformed, so `transforms[0]` is
        ignored. Every other entry is applied to its record's whole structure.

        `transforms` is positional rather than keyed by `pocket_id` on purpose: a query can be compared
        against itself, so the reference and a target may carry the same `pocket_id`.

        Args:
            aln_records (list): Ordered records, reference first. Fields read: `pocket_id`,
                `struct_path`, `chain_info`.
            transforms (list): Parallel to `aln_records`. Entry 0 is ignored; entry i is either a
                `(u, t)` pair -- a 3x3 rotation in gemmi's LEFT-multiplying convention plus a 3-vector
                translation -- or None to drop that record from the output. A caller passing None is
                expected to have logged why.
            out_path (str): Destination path for the aligned PDB file.

        Returns:
            None: Writes the merged PDB to `out_path`.
        """
        structs = []
        domain_chains = []
        motif_chains = []
        us = []
        ts = []
        kept_records = []
        dropped = []

        for i, record in enumerate(aln_records):
            try:
                if i == 0:
                    struct_u = np.eye(3)
                    struct_t = np.zeros(3)
                else:
                    if transforms[i] is None:
                        dropped.append(record["pocket_id"])
                        continue
                    struct_u, struct_t = transforms[i]

                struct = gemmi.read_structure(record["struct_path"])

                if record["chain_info"] is None:
                    domain_chain = 0  # first chain
                    motif_chain = None
                else:
                    chains = record["chain_info"].split("_")
                    domain_chain = chains[0]
                    if len(chains) > 1:
                        motif_chain = chains[1]
                    else:
                        motif_chain = None

                # If everything has been successful add it things to be processed
                structs.append(struct)
                us.append(struct_u)
                ts.append(struct_t)
                domain_chains.append(domain_chain)
                motif_chains.append(motif_chain)
                kept_records.append(record)

            except Exception as e:
                dropped.append(record["pocket_id"])
                self.logger.error(f"Problem processing {record['pocket_id']}: {e}", extra={"stage": "StructureAligner"})

        if dropped:
            self.logger.warning(
                f"Not superposing {dropped}; they are absent from {out_path}",
                extra={"stage": "StructureAligner"},
            )
        if not kept_records:
            self.logger.error(f"No structure could be placed, not writing {out_path}", extra=self._log_extra)
            return

        aligned_struct = self._apply_transformation(structs, domain_chains, motif_chains, us, ts)
        self._write_aligned(kept_records, aligned_struct, out_path)

    def _write_aligned(self, kept_records, aligned_struct, out_path):
        """
        Write a merged structure out as a PDB with a COMPND header naming each model.

        Takes the records that actually made it into `aligned_struct`, not the records the caller started
        with: a record dropped for want of a transform used to keep its COMPND entry, so the header named
        models the file did not contain.

        The chain labels must come from a fresh `_char_gen()` consumed in the same order
        `_apply_transformation` consumed its own, or the header and the coordinates disagree.

        Args:
            kept_records (list): Records present in `aligned_struct`, in model order.
            aligned_struct (gemmi.Structure): The merged structure to write.
            out_path (str): Destination path for the PDB file.

        Returns:
            None: Writes the PDB to `out_path`.
        """
        pdb_str = aligned_struct.make_pdb_string()

        model_nums = (str(x) for x in count(1))
        model_names = [record["pocket_id"] for record in kept_records]
        chain_names = self._char_gen()
        header = ""
        for model_num, model_name, chain_name in zip(model_nums, model_names, chain_names):
            line_nums = (str(x) for x in count(1))
            header += f"""
COMPND {next(line_nums).zfill(3)} MOL_ID: {model_num};
COMPND {next(line_nums).zfill(3)} MOLECULE: {model_name[:70]};
COMPND {next(line_nums).zfill(3)} CHAIN: {chain_name};
"""

        with open(out_path, "w") as f:
            f.write(header)
            f.write(pdb_str)

    def foldseek_transform(self, aln_records, alignment_df, out_path):
        """
        Build an aligned multi-structure PDB from Foldseek-style alignment results.

        The first record in `aln_records` is the reference; every other must have a row in `alignment_df`.
        The stored Foldseek transform strings are already in the LEFT-multiplying convention `transform`
        wants, so they are parsed and passed straight through.

        The local BLOSUM62 aligner writes "-" for `u` and `t`, so every target is dropped here and the
        output holds the query alone -- use the pocket transforms with `transform` instead (see the
        `align_struct_method` setting).

        Args:
            aln_records (list): Ordered alignment records, reference first. Fields read: `pocket_id`,
                `preprocess_name`, `struct_path`, `chain_info`.
            alignment_df (pandas.DataFrame): Alignment table carrying the Foldseek transforms, indexed by
                (query, target) `preprocess_name`.
            out_path (str): Destination path for the aligned PDB file.

        Returns:
            None: Writes the merged PDB to `out_path`.
        """
        query_preprocess_name = aln_records[0]["preprocess_name"]

        transforms = [None]  # the reference is placed untransformed
        for record in aln_records[1:]:
            try:
                row = alignment_df.loc[query_preprocess_name, record["preprocess_name"]]
                struct_u = np.array([float(x) for x in row["u"].split(",")]).reshape((3, 3))
                struct_t = np.array([float(x) for x in row["t"].split(",")])
                transforms.append((struct_u, struct_t))
            except Exception as e:
                transforms.append(None)
                self.logger.error(
                    f"Problem processing {record['pocket_id']}: {e}", extra={"stage": "foldseek_transform"}
                )

        self.transform(aln_records, transforms, out_path)
