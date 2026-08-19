"""
Utilities for building aligned multi-structure models.

The alignment workflow in PocketMapper needs a lightweight way to take one query
structure plus one or more matched target structures and write them back out as a
single, human-readable PDB file. This module provides that glue around `gemmi`:

- read mmCIF structures from disk
- apply per-target rotation and translation matrices
- rename chains so each aligned structure remains distinguishable in the output
- emit a combined PDB with simple COMPND metadata

The public entry point is :class:`StructureAligner`, which is currently used by
the PocketMapper structural-alignment step after pocket comparison.
"""

import logging
import string
from itertools import count, permutations
import gemmi
import numpy as np


class StructureAligner:
    def __init__(self):
        self.logger = logging.getLogger(__name__)
        self._log_extra = {"stage": "StructureAligner"}
        logging.debug("Initialized", extra=self._log_extra)

    def _char_gen(self):
        """
        Yield short, PDB-friendly chain identifiers for aligned output.

        The output structure may contain more chains than the original inputs.
        To avoid collisions with the preserved domain chain name ``0``, this
        generator emits:

        - single-character identifiers starting at ``1`` and then ``A``/``B``...
        - two-character identifiers when the single-character space is exhausted

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

        The first structure is treated as the reference frame. Each structure is
        transformed in place using the supplied rotation matrix and translation
        vector, then copied into a new combined ``gemmi.Structure``. The domain
        chain from each input is renamed to ``0`` in the output so it can be
        recognized consistently, while any motif chain is given a generated name.

        Args:
            structs (list[gemmi.Structure]): Structures to align and merge.
            domain_chains (list[str]): Domain chain IDs to preserve from each input structure.
            motif_chains (list[str | None]): Optional motif chain IDs to preserve from each input structure.
            us (list[numpy.ndarray]): 3x3 rotation matrices, one per structure.
            ts (list[numpy.ndarray]): Translation vectors of length 3, one per structure.

        Returns:
            gemmi.Structure: A merged aligned structure containing one model per input.
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

    def foldseek_transform(self, aln_records, alignment_df, out_path):
        """
        Build an aligned multi-structure PDB from Foldseek-style alignment results.

        The first record in ``aln_records`` is treated as the reference structure.
        Every subsequent record must correspond to a row in ``alignment_df`` whose
        index is the reference ``preprocess_name`` and whose columns include the
        target ``preprocess_name``. The stored Foldseek transform strings are
        parsed into a 3x3 rotation matrix ``u`` and a 3-vector translation ``t``.

        Each structure is read from ``record["struct_path"]`` using mmCIF parsing,
        transformed into the reference frame, and written to ``out_path`` as a
        single PDB containing one model per input record.

        Expected record fields:
            - ``pocket_id``: stable identifier used in output metadata
            - ``preprocess_name``: alignment key used to look up ``u`` and ``t``
            - ``struct_path``: path to the mmCIF or mmCIF.GZ structure file
            - ``chain_info``: domain chain or ``domain_motif`` pair

        Args:
            aln_records (list[dict]): Ordered alignment records, with the reference first.
            alignment_df (pandas.DataFrame): Alignment table containing Foldseek transforms.
            out_path (str): Destination path for the aligned PDB file.

        Returns:
            None
        """

        target_preprocess_name = aln_records[0]["preprocess_name"]
        structs = []
        domain_chains = []
        motif_chains = []
        us = []
        ts = []
        for i, record in enumerate(aln_records):
            try:
                # Load the structure
                struct_path = record["struct_path"]
                struct = gemmi.read_structure(struct_path)  # , format=gemmi.CoorFormat.Mmcif)

                # If the structure is not the target, get the transformation matrices from the alignment dataframe
                if i > 0:
                    row = alignment_df.loc[target_preprocess_name, record["preprocess_name"]]
                    struct_u = np.array([float(x) for x in row["u"].split(",")]).reshape((3, 3))
                    struct_t = np.array([float(x) for x in row["t"].split(",")])
                else:
                    struct_u = np.eye(3)
                    struct_t = np.zeros(3)

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

            except Exception as e:
                self.logger.error(
                    f"Problem processing {record['pocket_id']}: {e}", extra={"stage": "foldseek_transform"}
                )

        aligned_struct = self._apply_transformation(structs, domain_chains, motif_chains, us, ts)
        pdb_str = aligned_struct.make_pdb_string()

        model_nums = (str(x) for x in count(1))
        model_names = [record["pocket_id"] for record in aln_records]
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
