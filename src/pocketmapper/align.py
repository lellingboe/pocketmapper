"""
Function for aligning a community of pockets
    Cached and uncached version

Inputs:
    Path to structure folder

"""

import string
from itertools import cycle
import gemmi
import numpy as np


def char_gen():
    nice_chars = string.digits + string.ascii_letters
    nice_chars = nice_chars.replace("q", "")  # 9 is reserved for the domain chain
    for x in cycle(nice_chars):
        yield (x)


class Aligner:
    def __init__(self):
        pass

    def align_objects(self, structs, domain_chains, motif_chains, us, ts):
        # Align everything to the first struct
        ref_st, ref_dc = structs[0], domain_chains[0]
        ref_st[0][ref_dc].name = "0"

        for i, st, dc, mc, u, t in zip(range(len(us)), structs[1:], domain_chains[1:], motif_chains[1:], us, ts):
            print(i, st, dc, mc, u, t)
            ref_st.add_model(gemmi.Model(2))

            trans = gemmi.Transform(gemmi.Mat33(u), gemmi.Vec3(*t))
            st[0].transform_pos_and_adp(trans)

            st[0][dc].name = "0"
            ref_st[-1].add_chain(st[0]["0"])
            ref_st[-1].add_chain(st[0][mc], unique_name=True)
        ref_st.setup_entities()
        return ref_st

    def test(self):
        struct_path = r"/Users/lellingboe/Work/data/pocketmapper/example/pocketmapper_cache/pdb_structures/2w2h.cif.gz"
        st = gemmi.read_structure(struct_path, format=gemmi.CoorFormat.Mmcif)
        model = st[0]
        model.add_chain(model["B"], unique_name=True)
        del model["B"]
        del model["D"]
        del model["S"]
        test_out = r"/Users/lellingboe/Work/data/pocketmapper/align_test/test.cif"
        groups = gemmi.MmcifOutputGroups(False, atoms=True, group_pdb=True)
        st.make_mmcif_document(groups).write_file(test_out)


if __name__ == "__main__":
    aligner = Aligner()
    st1 = gemmi.read_structure(
        r"/Users/lellingboe/Work/data/pocketmapper/example/pocketmapper_cache/divided_structs/4M1Z_A_C.cif.gz",
        format=gemmi.CoorFormat.Mmcif,
    )
    st2 = gemmi.read_structure(
        r"/Users/lellingboe/Work/data/pocketmapper/example/pocketmapper_cache/divided_structs/2W2M_A_E.cif.gz",
        format=gemmi.CoorFormat.Mmcif,
    )
    st1.setup_entities()
    st2.setup_entities()

    structs = [st1, st2]
    domain_chains = ["A", "A"]
    motif_chains = ["C", "E"]

    u = np.array([float(x) for x in "-0.516,-0.856,0.024,0.486,-0.270,0.831,-0.705,0.441,0.555".split(",")]).reshape(
        (3, 3)
    )
    t = np.array([float(x) for x in "-9.822,-25.162,-28.802".split(",")])
    us = [u]
    ts = [t]

    test_st = aligner.align_objects(structs, domain_chains, motif_chains, us, ts)

    # groups = gemmi.MmcifOutputGroups(False, atoms=True, group_pdb=True)
    test_st.make_mmcif_document().write_file(r"/Users/lellingboe/Work/data/pocketmapper/align_test/test.cif")
