import gzip
import hashlib
import shutil
import os
import logging
from copy import deepcopy
from Bio.SVDSuperimposer import SVDSuperimposer
from Bio.PDB import MMCIFParser
from glob import glob
from itertools import product
from tqdm import tqdm
import json
from collections import defaultdict
import pandas as pd
import gemmi
from numpy import array
from numpy import linalg as LA
from pocketmapper import pisa_downloader

from pocketmapper.constants import SINGLE_AA_CODE, VDW_RADII


def pdb_preprocessing_gemmi(df, ref_dir, cache_dir, out_dir):
    """
    Docstring for pdb_preprocessing_gemmi

    :param df: Description
    :param ref_dir: directory for reference pdb files to be divided
    :param cache_dir: directory for divided pdbs to be cached
    :param out_dir: directory to be used with foldseek
    """
    status_dict = {}
    stage = {"stage": "Dividing structures"}

    for i, row in tqdm(df.iterrows()):
        pdb = row["struct_info"]
        chain_info = row["chain_info"]  # e.g. A_B or A
        chains = chain_info.split("_")  # [A, B] or [A]
        try:
            # Ensuring divided structure is in the cache directory
            while len(chains) > 0:
                pdb_chains = pdb + "_" + "_".join(chains)
                cache_path = os.path.join(cache_dir, f"{pdb_chains}.cif")
                cache_path_gz = cache_path + ".gz"

                if not os.path.exists(cache_path_gz):
                    ref_path = os.path.join(ref_dir, f"{pdb}.cif.gz")
                    st = gemmi.read_structure(ref_path, format=gemmi.CoorFormat.Mmcif)

                    # Taking first model and deleting the rest
                    del st[1:]
                    model = st[0]

                    # verify structure contains all interaction chains
                    model_chains = set([chain.name for chain in model])
                    if not set(chains).issubset(model_chains):
                        msg = f"Preprocessing: {pdb} does not contain all interaction chains {chains}"
                        logging.warning(
                            msg,
                            extra=stage,
                        )
                        status_dict[i] = False
                        chains = []  # to skip to next pdb
                        continue

                    # Detaching all non interaction chains
                    for chain_id in model_chains:
                        if chain_id not in chains:
                            del model[chain_id]

                    # Output the domain and motif pdb file
                    groups = gemmi.MmcifOutputGroups(False, atoms=True, group_pdb=True)
                    st.make_mmcif_document(groups).write_file(cache_path)
                    with open(cache_path, "rb") as f_in:
                        with gzip.open(cache_path_gz, "wb") as f_out:
                            shutil.copyfileobj(f_in, f_out)
                    os.remove(cache_path)

                out_path_gz = os.path.join(out_dir, f"{pdb_chains}.cif.gz")
                shutil.copyfile(cache_path_gz, out_path_gz)  # copying to foldseek directory

                status_dict[i] = True
                chains = chains[:-1]  # e.g. A_B -> A

        except Exception as e:
            logging.warning(f"Could not divide {pdb} with chain info {chain_info}", extra=stage)
            logging.debug("Exception info", exc_info=e, extra=stage)
            status_dict[i] = False

    return status_dict


def calculate_pockets(df, target_dir, query_dir, pocket_dir):
    """Takes in a path to a pdb file"""
    parser = MMCIFParser()
    pocket_cache = glob(pocket_dir + "/*.json")

    all_problem_atoms = defaultdict(lambda: 0)
    all_problem_residues = defaultdict(lambda: 0)
    pocket_dict = {}
    for i, row in tqdm(df.iterrows()):
        pocket_path = os.path.join(pocket_dir, f"{row.pdb_domain_motif}.json")
        if pocket_path in pocket_cache:  # If cache exists, just load that
            with open(pocket_path, "r") as f:
                pocket = json.load(f)
        else:
            if row.type == "query":
                tmp_dir = query_dir
            if row.type == "target":
                tmp_dir = target_dir

            # Load the structure
            try:
                structure = parser.get_structure(
                    row.pdb_domain_motif,
                    os.path.join(tmp_dir, f"{row.pdb_domain_motif}.cif"),
                )
            except Exception:
                logging.exception(
                    f"Error parsing structure {row.pdb_domain_motif}",
                    extra={"stage": "Calculating Pockets"},
                )
                continue

            # Calculate the pocket from that structure
            try:
                pocket, problem_atoms, problem_residues = pocket_overlap(structure, row.domain_chain, row.motif_chain)
            except Exception:
                logging.exception(
                    f"Error calculating pocket {row.pdb_domain_motif}",
                    extra={"stage": "Calculating Pockets"},
                )
                continue

            # Update problem cases
            for atom in problem_atoms:
                all_problem_atoms[atom] += 1
            for res in problem_residues:
                all_problem_residues[res] += 1

            with open(pocket_path, "w") as f:
                json.dump(pocket, f)
        pocket_dict[row.pdb_domain_motif] = pocket

    return pocket_dict, all_problem_atoms, all_problem_residues


def jsonify_dict(item):
    """
    Recursively looks for sets in a dictionary and turns then into lists
    This allows dicts with sets to become JSON serializeable
    """
    if isinstance(item, set):
        return list(item)
    elif isinstance(item, dict):
        return {str(k): jsonify_dict(v) for k, v in item.items()}
    else:
        return item


def safe_filename(name, max_len=80):
    """
    Build a filesystem-safe filename stem from a pocket_id.

    pocket_ids can embed long comma-separated residue lists (e.g. for
    passthrough/VDW queries), which can exceed OS filename length limits.
    Names longer than max_len are truncated and given an md5 suffix to keep
    them unique.
    """
    safe_name = name.replace(":", "_").replace(",", "_")
    if len(safe_name) <= max_len:
        return safe_name
    name_hash = hashlib.md5(name.encode()).hexdigest()
    return f"{safe_name[:max_len]}_{name_hash}"


# reimplement with scipy.spatial.distance.cdist
def pocket_overlap(structure, domain_chain, motif_chain):
    """
    structure: Biopython model
    chain1, chain2 : Strings -> Chain IDs
    """

    model = structure[0]

    pocket_res_ids = dict()
    motif_res_ids = dict()
    full_interaction = dict()

    problem_atoms = set()
    problem_residues = set()

    # Filter out hetatoms
    domain_residues = [x for x in model[domain_chain].get_residues() if x.id[0] != "W"]  # removing water molecules
    motif_residues = [x for x in model[motif_chain].get_residues()]

    for res1, res2 in product(domain_residues, motif_residues):
        # atom ordering per residue: ['N', 'CA', 'C', 'O', 'CB', R1, R1, ...]
        if res1.get_resname() == "GLY":
            backbone1 = [0, 2, 3]
        else:
            backbone1 = [0, 1, 2, 3]
        if res2.get_resname() == "GLY":
            backbone2 = [0, 2, 3]
        else:
            backbone2 = [0, 1, 2, 3]

        for (pos1, atom1), (pos2, atom2) in product(enumerate(res1.get_atoms()), enumerate(res2.get_atoms())):
            distance = atom1 - atom2
            if distance > 5:
                continue

            # Skipping pocket residues not in the standard 20
            if atom1.parent.resname not in SINGLE_AA_CODE:
                problem_residues.add(res1.resname)
                continue

            # VDW Radii
            try:
                vdw1 = VDW_RADII[atom1.id[0]]
            except KeyError:
                problem_atoms.add(atom1.id)
                continue
            try:
                vdw2 = VDW_RADII[atom2.id[0]]
            except KeyError:
                problem_atoms.add(atom2.id)
                continue

            vdw_range = vdw1 + vdw2
            overlap = vdw_range - distance
            if overlap > -0.4:
                (full_interaction.setdefault(res1.id[1], dict()).setdefault(res2.id[1], set())).add(
                    (pos1 not in backbone1, pos2 not in backbone2)
                )

                pocket_res_ids.setdefault(res1.id[1], False)
                if pos1 not in backbone1:
                    pocket_res_ids[res1.id[1]] = True
                motif_res_ids.setdefault(res2.id[1], False)
                if pos2 not in backbone1:
                    motif_res_ids[res2.id[1]] = True

    if len(problem_atoms) > 0:
        logging.warning(
            f"No vdw radius for {list(problem_atoms)} in {structure.id}",
            extra={"stage": "Calculating Pocket"},
        )
    if len(problem_residues) > 0:
        logging.warning(
            f"No single AA code for {problem_residues}: {structure.id}",
            extra={"stage": "Calculating Pocket"},
        )

    # Dict for mapping residue id to sequence position
    res_id_to_pos = {}
    res_pos_coords = {}
    seq = []
    for i, res in enumerate(domain_residues):
        atoms = list(res.get_atoms())
        if len(atoms) > 1 and atoms[1].id == "CA":
            res_id_to_pos[res.id[1]] = i
            res_pos_coords[i] = atoms[1].coord.tolist()
        seq.append(SINGLE_AA_CODE.get(res.get_resname(), "X"))
    seq = "".join(seq)

    # mapping pocket ids to sequence position for foldseek
    if pocket_res_ids:
        pocket_res_pos = {res_id_to_pos[k]: v for k, v in pocket_res_ids.items() if k in res_id_to_pos}

    pocket = jsonify_dict(
        {
            "pocket_exists": len(pocket_res_ids) > 0,
            "pocket_res_ids": pocket_res_ids,
            "pocket_res_pos": pocket_res_pos,
            "res_id_to_pos": res_id_to_pos,
            "pocket_to_motif_sidechain_overlap": full_interaction,
            "res_pos_coords": res_pos_coords,
            "seq": seq,
        }
    )

    return pocket, problem_atoms, problem_residues


"""
Foldseek output format used for comparison:
0 query
1 target
2 fident
3 alnlen
4 mismatch
5 gapopen
6 qstart
7 qend
8 tstart
9 tend
10 evalue
11 lddt
12 qaln
13 taln
14 u
15 t
16 qseq
17 tseq
"""


def compare_pockets(
    alignment_df,
    pocket_dict,
    preproc_to_ids,
    blosum_path=r"/home/data/motif_aligner/blosum62.bla",
    alphafold=False,
):
    """
    Compare two pockets based on foldseek alignment
    """

    # stage = {"stage": "Pocket Comparison"}
    blosum_similarity_matrix = read_blast_similarity_matrix(blosum_path)

    unknown_ids = defaultdict(lambda: defaultdict(set))  # for saving tri-code ids which are unknown
    incorrect_mapping = defaultdict(dict)  # for saving cases where foldseek mapping doesn't match pocketmapper sequence

    # Setting up vars for use later
    existing_calcs = set()
    output_rows = []
    sup = SVDSuperimposer()

    # TODO divide this into common things for each pocket and a cross-comparison
    for row in tqdm(alignment_df.itertuples(index=False)):
        domain_1 = row[0]
        domain_2 = row[1]
        try:

            # MAPPING SEQUENCE POSITION TO ALIGNMENT POSITION #####
            p1_seq_to_aln = {}
            i = 0
            for j, res in enumerate(row[12]):
                if res != "-":
                    p1_seq_to_aln[i] = j
                    i += 1
            # p1_aln_to_seq = {v: k for k, v in p1_seq_to_aln.items()}

            p2_seq_to_aln = {}
            i = 0
            for j, res in enumerate(row[13]):
                if res != "-":
                    p2_seq_to_aln[i] = j
                    i += 1
            # p2_aln_to_seq = {v: k for k, v in p2_seq_to_aln.items()}

            # GETTING POCKETS WHICH CORRESPOND TO THE FOLDSEEK NAME
            pockets_1 = {}
            for pocket_id in preproc_to_ids.get(domain_1):
                if pocket_id in pocket_dict:
                    pockets_1[pocket_id] = pocket_dict[pocket_id]
            if len(pockets_1) == 0:
                continue

            if alphafold:
                pockets_2 = {
                    domain_2: {
                        "res_auth_ids": [str(k) for k in range(row[9])],
                        "id_pos_codes_match": True,
                        "pocket_exists": True,
                        "has_coords": False,
                        "ca_sequence": row[17],
                    }
                }
                pockets_2[domain_2].update({str(k): {"seq_pos": k} for k in range(row[9])})
            else:
                pockets_2 = {}
                for pocket_id in preproc_to_ids.get(domain_2):
                    if pocket_id in pocket_dict:
                        pockets_2[pocket_id] = pocket_dict[pocket_id]
                if len(pockets_2) == 0:
                    continue

            # Iterating through aligned pairs
            for pocket_id_1, pocket_id_2 in product(pockets_1.keys(), pockets_2.keys()):

                # No self comparisons
                # if interaction_1 == interaction_2:
                #    continue

                # Checking for A-B comparison if B-A has already been calculated
                if (pocket_id_1, pocket_id_2) in existing_calcs:
                    continue
                existing_calcs.add((pocket_id_1, pocket_id_2))

                if pocket_id_1 in incorrect_mapping or pocket_id_2 in incorrect_mapping:
                    continue

                # Starting the output list
                output = {
                    "pocket_1": pocket_id_1,
                    "pocket_2": pocket_id_2,
                    "evalue": row[10],
                    "lddt": row[11],
                }

                # Getting the pockets
                p1 = deepcopy(pockets_1.get(pocket_id_1))
                p2 = deepcopy(pockets_2.get(pocket_id_2))
                if not p1["pocket_exists"] or not p2["pocket_exists"]:
                    continue

                p1_seq_identity = sum(map(str.__eq__, row[16], p1["ca_sequence"])) / len(p1["ca_sequence"])
                if p1_seq_identity < 0.8:
                    incorrect_mapping[pocket_id_1] = {
                        "p1_seq_identity": p1_seq_identity,
                        "p1_seq": p1["ca_sequence"],
                        "fs_seq": row[16],
                    }

                p2_seq_identity = sum(map(str.__eq__, row[17], p2["ca_sequence"])) / len(p2["ca_sequence"])
                if p2_seq_identity < 0.8:
                    incorrect_mapping[pocket_id_2] = {
                        "p2_seq_identity": p2_seq_identity,
                        "p2_seq": p2["ca_sequence"],
                        "fs_seq": row[17],
                    }

                # Calculated Metrics
                # pocket 1
                p1_adj = 1 - row[6]  # row[6] is qstart
                p1_fs_len = row[7] - row[6] + 1  # row[7] is qend
                p1_in_aln_region_count = 0
                p1["foldseek_pos"] = []
                for res in p1["res_auth_ids"]:
                    fs_start_adj_pos = int(p1[res]["seq_pos"]) + p1_adj
                    in_aln_region = -1 < fs_start_adj_pos < p1_fs_len
                    p1[res]["in_fs_aln_region"] = in_aln_region
                    if in_aln_region:
                        p1_in_aln_region_count += 1
                        fs_res_pos = p1_seq_to_aln[fs_start_adj_pos]
                        p1["foldseek_pos"].append(fs_res_pos)
                        p1[res]["fs_pos"] = fs_res_pos
                        p1[res]["fs_res_code"] = row[12][fs_res_pos]

                        # Checking fs single res code and pocektmapper single res codes match
                        if p1[res]["fs_res_code"] != p1[res]["res_code_single"]:
                            unknown_ids[p1[res]["fs_res_code"]][p1[res]["res_code"]].add(
                                ",".join([pocket_id_2, pocket_id_1, res])
                            )

                output["pocket_1_res_ids"] = ",".join(p1["res_auth_ids"])
                output["pocket_1_len"] = len(p1["res_auth_ids"])
                output["pocket_1_seq"] = "".join([p1[x]["res_code_single"] for x in p1["res_auth_ids"]])
                output["pocket_1_pct_aln"] = p1_in_aln_region_count / len(p1["res_auth_ids"])

                # pocket 2
                p2_adj = 1 - row[8]
                p2_fs_len = row[9] - row[8] + 1
                p2_in_aln_region_count = 0
                p2["foldseek_pos"] = []
                for res in p2["res_auth_ids"]:
                    fs_start_adj_pos = int(p2[res]["seq_pos"]) + p2_adj
                    in_aln_region = -1 < fs_start_adj_pos < p2_fs_len
                    p2[res]["in_fs_aln_region"] = in_aln_region
                    if in_aln_region:
                        p2_in_aln_region_count += 1
                        fs_res_pos = p2_seq_to_aln[fs_start_adj_pos]
                        p2["foldseek_pos"].append(fs_res_pos)
                        p2[res]["fs_pos"] = fs_res_pos
                        p2[res]["fs_res_code"] = row[13][fs_res_pos]

                        # Checking fs single res code and pocektmapper single res codes match

                        if not alphafold and (p2[res]["fs_res_code"] != p2[res]["res_code_single"]):
                            debug_id = ",".join([pocket_id_1, pocket_id_2, res])
                            unknown_ids[p2[res]["fs_res_code"]][p2[res]["res_code"]].add(debug_id)

                if not alphafold:
                    output["pocket_2_res_ids"] = ",".join(p2["res_auth_ids"])
                    output["pocket_2_len"] = len(p2["res_auth_ids"])
                    output["pocket_2_seq"] = "".join([p2[x]["res_code_single"] for x in p2["res_auth_ids"]])
                    output["pocket_2_pct_aln"] = p2_in_aln_region_count / len(p2["res_auth_ids"])

                # OVERLAP
                overlapping_residues = [x for x in p1["foldseek_pos"] if x in p2["foldseek_pos"]]
                output["overlap_count"] = len(overlapping_residues)
                if len(overlapping_residues) == 0:  # If no overlapping resides, end here
                    output_rows.append(output)
                    continue

                # overlap ids
                p1_overlap_ids = []
                for res in p1["res_auth_ids"]:
                    if p1[res].get("fs_pos", -1) in overlapping_residues:
                        p1_overlap_ids.append(res)
                output["pocket_1_overlap_ids"] = ",".join(p1_overlap_ids)
                p2_overlap_ids = []
                for res in p2["res_auth_ids"]:
                    if p2[res].get("fs_pos", -1) in overlapping_residues:
                        p2_overlap_ids.append(res)
                output["pocket_2_overlap_ids"] = ",".join(p2_overlap_ids)

                # percent overlap
                p1_pct_overlap = len(overlapping_residues) / len(p1["res_auth_ids"])
                output["pocket_1_pct_overlap"] = p1_pct_overlap

                if not alphafold:
                    p2_pct_overlap = len(overlapping_residues) / len(p2["res_auth_ids"])
                    output["pocket_2_pct_overlap"] = p2_pct_overlap
                    output["min_pct_overlap"] = min(p1_pct_overlap, p2_pct_overlap)
                    output["max_pct_overlap"] = max(p1_pct_overlap, p2_pct_overlap)

                #####
                # Aligned residues as a sequence
                p1_aln_seq = "".join([row[12][x] for x in overlapping_residues])
                p2_aln_seq = "".join([row[13][x] for x in overlapping_residues])
                output["pocket_1_seq_overlap"] = p1_aln_seq
                output["pocket_2_seq_overlap"] = p2_aln_seq

                # Identity
                overlap_identity = sum(map(str.__eq__, p1_aln_seq, p2_aln_seq)) / len(overlapping_residues)
                output["overlap_identity"] = overlap_identity

                # BLOSUM62 similarity
                similarity = binary_similarity(p1_aln_seq, p2_aln_seq, blosum_similarity_matrix)
                output["overlap_similarity_binary"] = similarity

                similarity_1_2 = full_similarity(p1_aln_seq, p2_aln_seq, blosum_similarity_matrix)
                similarity_2_1 = full_similarity(p2_aln_seq, p1_aln_seq, blosum_similarity_matrix)
                output["overlap_similarity_1_2"] = similarity_1_2
                output["overlap_similarity_2_1"] = similarity_2_1
                output["min_overlap_similarity"] = min(similarity_1_2, similarity_2_1)
                output["max_overlap_similarity"] = max(similarity_1_2, similarity_2_1)

                # RMSD dings
                if p1["has_coords"] and p2["has_coords"]:
                    x = array([p1[str(x)]["ca_coords"] for x in p1_overlap_ids])
                    y = array([p2[str(x)]["ca_coords"] for x in p2_overlap_ids])

                    if len(overlapping_residues) > 2:
                        sup.set(x, y)
                        sup.run()
                        u, t = sup.get_rotran()
                        output["p2_to_p1_u"] = u.flatten().tolist()
                        output["p2_to_p1_t"] = t.tolist()

                        # TODO do this with matrix algebra isntead of doing it twice
                        sup.set(y, x)
                        sup.run()
                        u, t = sup.get_rotran()
                        output["p1_to_p2_u"] = u.flatten().tolist()
                        output["p1_to_p2_t"] = t.tolist()
                        output["rmsd"] = sup.get_rms()

                        x_on_y = sup.get_transformed()
                        ca_dists = LA.norm(x_on_y - y, axis=1)
                        output["ca_dists"] = ",".join([str(round(x, 3)) for x in ca_dists])

                output_rows.append(output)

        except KeyError:
            logging.warning(
                f"Uncontrolled KeyError calculating {domain_1} and {domain_2}",
                extra={"stage": "Pocket Comparison"},
            )
            raise
        except Exception:
            logging.exception(
                f"Uncontrolled error calculating {domain_1} and {domain_2}",
                extra={"stage": "Pocket Comparison"},
            )
            raise

    return pd.DataFrame.from_dict(output_rows), unknown_ids, incorrect_mapping


def binary_similarity(seqA, seqB, similarity_matrix):
    """
    Similarity of A->B and B->A
    score = 1 if full similarity score > 0 else score = 0
    normalized to the length of the seqence
    """
    seqA = seqA.replace("U", "X").upper()
    seqB = seqB.replace("U", "X").upper()

    similarity = [(similarity_matrix[x][y] > 0) for x, y in zip(seqA, seqB)]  # True or False
    similarity_score = sum(similarity) / len(similarity)
    return similarity_score


def full_similarity(seqA, seqB, similarity_matrix):
    """
    Similarity of A->B
    blosum scores normalized by
    normalized to the length of the seqence
    """
    seqA = seqA.replace("U", "X").upper()
    seqB = seqB.replace("U", "X").upper()

    similarity = [similarity_matrix[x][y] for x, y in zip(seqA, seqB)]
    similarity_max = [similarity_matrix[x][x] for x in seqA]
    similarity_normalized = [x / y for x, y in zip(similarity, similarity_max)]
    return sum(similarity_normalized) / len(similarity_normalized)


def read_blast_similarity_matrix(similarity_matrix_path, delimiter=" "):
    similarity_matrix = {}
    file_content = open(similarity_matrix_path).read().strip().split("\n")
    header = None
    idx_to_aa = None
    max_score = 0
    row_counter = 0
    for line in file_content:
        # Skip comment lines
        if line[0] == "#":
            continue

        # parsing the first lines into a dict
        if header is None:
            if delimiter == " ":
                header = line.strip().split()  # splits on whitespace and discards empty results
            else:
                header = line.strip().split(delimiter)
            idx_to_aa = dict(list(zip(list(range(0, len(header))), header)))
        else:
            if delimiter == " ":
                fields = line.strip().split()
            else:
                fields = line.strip().split(delimiter)
            from_aa = idx_to_aa[row_counter]
            row_counter += 1
            similarity_matrix[from_aa] = {}
            for idx, score in enumerate(fields):
                to_aa = idx_to_aa[idx]
                similarity_matrix[from_aa][to_aa] = float(score)
                if to_aa not in similarity_matrix:
                    similarity_matrix[to_aa] = {}
                similarity_matrix[to_aa][from_aa] = float(score)
                if float(score) > max_score:
                    max_score = float(score)
            similarity_matrix[from_aa]["-"] = -4

    # Add and construct the "-" entry
    similarity_matrix["-"] = {}
    for aa in similarity_matrix:
        if aa == "-":
            similarity_matrix["-"][aa] = -1
        else:
            similarity_matrix["-"][aa] = -4
    return similarity_matrix


def download_pisa_info(pdb_list, summary_dir, assembly_dir, interface_dir):
    downloader = pisa_downloader.PisaDownloader()
    downloader.get_interfaces(pdb_list, summary_dir, assembly_dir, interface_dir)
