"""
Code related to creating an aligned structure

StructureAligner:
- align: method which takes paths to files and returns aligned gemmi structures
- align_objects: align multiple gemmi structures given transformation matrices
"""

import logging
import string
from itertools import count, permutations
import gemmi
import numpy as np


class StructureAligner:
    def __init__(self):
        pass

    def _char_gen(self):
        nice_chars = string.digits + string.ascii_letters
        for x in nice_chars[1:]:  # 0 is reserved for domain names
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
        -us: list of rotation matrices (numpy arrays) to apply to each structure (except the first)
        -ts: list of translation vectors (numpy arrays) to apply to each structure (except the first)
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
        Aligns a set of structures based on the transformation matrices in the alignment dataframe
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
                struct = gemmi.read_structure(struct_path, format=gemmi.CoorFormat.Mmcif)

                # If the structure is not the target, get the transformation matrices from the alignment dataframe
                if i > 0:
                    row = alignment_df.loc[target_preprocess_name, record["preprocess_name"]]
                    struct_u = np.array([float(x) for x in row["u"].split(",")]).reshape((3, 3))
                    struct_t = np.array([float(x) for x in row["t"].split(",")])
                else:
                    struct_u = np.eye(3)
                    struct_t = np.zeros(3)

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
                logging.error(f"Problem processing {record['pocket_id']}: {e}")

        aligned_struct = self._apply_transformation(structs, domain_chains, motif_chains, us, ts)
        pdb_str = aligned_struct.make_pdb_string()

        model_nums = (str(x) for x in count(1))
        model_names = [record["pocket_id"].replace(":", "_").replace(",", "_") for record in aln_records]
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
