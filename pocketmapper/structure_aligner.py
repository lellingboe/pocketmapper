"""
Code related to creating an aligned structure

StructureAligner:
- align: method which takes paths to files and returns aligned gemmi structures
- align_objects: align multiple gemmi structures given transformation matrices
"""

import string
from itertools import count, permutations
import gemmi
import pandas as pd
import os
import numpy as np


class StructureAligner:
    def __init__(self):
        pass

    def _char_gen(self):
        nice_chars = string.digits + string.ascii_letters
        for x in nice_chars[1:]:  # Skipping 0 for domain names
            yield x
        for x, y in permutations(nice_chars, 2):
            yield (x + y)

    def _apply_transformation(self, structs, domain_chains, motif_chains, us, ts):
        """
        Function which alignes a set of structure given rotation matrices and translation vectors

        #length n lists
        -structs: list of gemmi structures
        -domain_chains: list of chain ids for the domain chains in each structure
        -motif_chains: list of chain ids for the motif chains in each structure

        # length n-1 lists
        -us: list of rotation matrices (numpy arrays) to apply to each structure (except the first)
        -ts: list of translation vectors (numpy arrays) to apply to each structure (except the first)
        """
        # Align everything to the first struct
        ref_st, ref_dc, ref_mc = structs[0], domain_chains[0], motif_chains[0]
        ref_st[0][ref_dc].name = "0"
        chain_names = self._char_gen()
        ref_st[0][ref_mc].name = next(chain_names)

        for i, cn, st, dc, mc, u, t in zip(
            count(2), chain_names, structs[1:], domain_chains[1:], motif_chains[1:], us, ts
        ):
            ref_st.add_model(gemmi.Model(i))

            # Apply transformation
            trans = gemmi.Transform(gemmi.Mat33(u), gemmi.Vec3(*t))
            st[0].transform_pos_and_adp(trans)

            # Add chains to the reference structure
            st[0][dc].name = "0"
            st[0][mc].name = cn
            ref_st[-1].add_chain(st[0]["0"])
            ref_st[-1].add_chain(st[0][cn])

        ref_st.setup_entities()
        return ref_st

    def foldseek_transform(self, struct_names, alignment_path, struct_dir, out_path):
        alignments = pd.read_csv(alignment_path, sep="\t", engine="c")
        alignments = alignments.set_index(["query", "target"])

        target = struct_names[0]
        structs = []
        domain_chains = []
        motif_chains = []
        us = []
        ts = []
        for struct_name in struct_names:
            try:
                struct_path = os.path.join(struct_dir, f"{struct_name}.cif.gz")
                if os.path.exists(struct_path):
                    struct = gemmi.read_structure(struct_path, format=gemmi.CoorFormat.Mmcif)

                row_id = (target[:-2], struct_name[:-2])
                if struct_name != target and row_id in alignments.index:
                    row = alignments.loc[row_id]
                    struct_u = np.array([float(x) for x in row["u"].split(",")]).reshape((3, 3))
                    struct_t = np.array([float(x) for x in row["t"].split(",")])

                # If everything has been successful add it things to be processed
                structs.append(struct)
                domain_chains.append(struct_name.split("_")[1])
                motif_chains.append(struct_name.split("_")[2])
                if struct_name != target:
                    us.append(struct_u)
                    ts.append(struct_t)
            except Exception as e:
                print(f"Problem processing {struct_name}: {e}")

        aligned_struct = self._apply_transformation(structs, domain_chains, motif_chains, us, ts)
        pdb_str = aligned_struct.make_pdb_string()
        aligned_struct.write_pdb(out_path)

        model_nums = (str(x) for x in count(1))
        header = ""
        for model_num, pdb_name, chain_name in zip(model_nums, struct_names, self._char_gen()):
            header += f"""
COMPND{next(model_nums).zfill(4)} MOL_ID: {model_num};\n
COMPND{next(model_nums).zfill(4)} MOLECULE: {pdb_name};\n
COMPND{next(model_nums).zfill(4)} CHAIN: {chain_name};\n
"""

        with open(out_path, "w") as f:
            f.write(header)
            f.write(pdb_str)
