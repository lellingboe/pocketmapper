import gzip
import os
import logging
from urllib.request import urlcleanup, urlretrieve
from Bio.SVDSuperimposer import SVDSuperimposer
from Bio.PDB import MMCIFParser, MMCIFIO
from glob import glob
from concurrent.futures import ThreadPoolExecutor
from itertools import repeat, product
from tqdm import tqdm
import json
from collections import defaultdict
import numpy as np
import pandas as pd

# TODO keep phospho information
SINGLE_AA_CODE = {
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
}
VDW_RADII = {"C": 1.88, "N": 1.64, "O": 1.46, "S": 1.77, "P": 1.87, "H": 1.0}
# https://www.cgl.ucsf.edu/chimerax/docs/user/commands/clashes.html


def get_mmcif(pdb_code, out_dir, cache):
    pdb_code = pdb_code.lower()
    out_fname = os.path.join(out_dir, f"{pdb_code}.cif")
    if not (out_fname in cache):

        url = f"https://files.wwpdb.org/pub/pdb/data/structures/divided/mmCIF/{pdb_code[1:3]}/{pdb_code}.cif.gz"
        gz_fname = os.path.join(out_dir, f"{pdb_code}.temp.gz")
        try:
            urlcleanup()
            urlretrieve(url, gz_fname)
        except OSError:
            logging.warning(f"get_mmcif: Could not download {pdb_code}")
            return (pdb_code, False)
        except Exception as e:
            logging.warning(f"Atypical issue when downloading {pdb_code}")
        else:
            with gzip.open(gz_fname, "rb") as gz:
                with open(out_fname, "wb") as out:
                    out.writelines(gz)
            os.remove(gz_fname)
    return (pdb_code, True)


def get_mmcifs(pdb_list, out_dir):
    cache = glob(os.path.join(out_dir, "*.cif"))
    with ThreadPoolExecutor(max_workers=100) as e:
        result = e.map(
            get_mmcif,
            pdb_list,
            repeat(out_dir),
            repeat(cache),
        )
    return {x.upper():y  for x,y in result}

# TODO thread pool executor version
def pdb_preprocessing(queries, ref_dir, domain_dir, motif_dir):
    """
    queries: a lit of tuple of the form (pdb_id, domain_chains, motif_chains)
        all tuple elements are strings

    writes out .pdb files to self.pdb_directory
    """
    status_dict = {}
    parser = MMCIFParser(QUIET=True)
    io = MMCIFIO()

    if not os.path.exists(domain_dir):
        os.mkdir(domain_dir)
    if not os.path.exists(motif_dir):
        os.mkdir(motif_dir)

    motif_cache = glob(os.path.join(motif_dir, "*.cif"))
    domain_cache = glob(os.path.join(domain_dir, "*.cif"))

    for pdb_id, motif_chain, domain_chains in tqdm(queries):
        name = f"{pdb_id}_{domain_chains}_{motif_chain}"
        try:
            ref_path = os.path.join(ref_dir,f"{pdb_id}.cif")
            domain_out = os.path.join(domain_dir,f"{pdb_id}_{domain_chains}.cif")
            motif_out = os.path.join(motif_dir,f"{name}.cif")

            interaction_chains = list(domain_chains + motif_chain)

            if (motif_out in motif_cache) and (domain_out in domain_cache):
                status_dict[f"{name}"] = True
            else:
                structure = parser.get_structure(pdb_id, ref_path)

                # Taking first model and detaching the rest
                model_gen = structure.get_models()
                model = next(model_gen)
                for dup_model in list(model_gen):
                    structure.detach_child(dup_model.id)

                # verify structure contains all interaction chains
                model_chains = {x.id for x in model.get_chains()}
                if not set(interaction_chains).issubset(model_chains):
                    logging.warning(
                        f"Preprocessing: {pdb_id} does not contain all interaction chains {interaction_chains}"
                    )
                    status_dict[f"{name}"] = False
                    continue
                else:
                    status_dict[f"{name}"] = True

                # Detaching all non interaction chains
                for chain in list(model.get_chains()):
                    if chain.id not in interaction_chains:
                        model.detach_child(chain.id)

                # Output the domain and motif pdb file
                io.set_structure(structure)
                io.save(motif_out)

                # output the domain pdb file
                model.detach_child(motif_chain)
                io.set_structure(structure)
                io.save(domain_out)
        except Exception as e:
            logging.warning(f"Could not divide {pdb_id}")
            status_dict[f"{name}"] = False


    return status_dict


def calculate_pockets(queries, motif_dir, pocket_dir):
    """Takes in a path to a pdb file"""
    parser = MMCIFParser()

    if not os.path.exists(pocket_dir):
        os.mkdir(pocket_dir)
    pocket_cache = glob(pocket_dir + "/*.json")

    all_problem_atoms = defaultdict(lambda:0)
    all_problem_residues = defaultdict(lambda:0)
    pocket_dict = {}
    for pdb_id, motif_chain, domain_chains in tqdm(queries):
        pocket_name = f"{pdb_id}_{domain_chains}_{motif_chain}"
        pocket_path = os.path.join(pocket_dir,f"{pocket_name}.json")
        if pocket_path in pocket_cache:  # If cache exists, just load that
            with open(pocket_path, "r") as f:
                pocket = json.load(f)
        else:
            # Load the structure
            try:
                structure = parser.get_structure(
                    pocket_name, os.path.join(motif_dir,f"{pocket_name}.cif")
                )
            except:
                logging.exception(f"Error parsing structure {pocket_name}")
                continue

            # Calculate the pocket from that structure
            try:
                pocket, problem_atoms, problem_residues = pocket_overlap(structure, domain_chains, motif_chain)
            except:
                logging.exception(f"Error calculating pocket {pocket_name}")
                continue

            # Update problem cases
            for atom in problem_atoms:
                all_problem_atoms[atom] +=1
            for res in problem_residues:
                all_problem_residues[res] +=1

            with open(pocket_path, "w") as f:
                json.dump(pocket, f)
        pocket_dict[pocket_name] = pocket

    return pocket_dict, all_problem_atoms, all_problem_residues


def sets_to_lists(item):
    """
    Recursively looks for sets in a dictionary and turns then into lists
    This allows dicts with sets to become JSON serializeable
    """
    if isinstance(item, set):
        return list(item)
    elif isinstance(item, dict):
        return {k: sets_to_lists(v) for k, v in item.items()}
    else:
        return item

# reimplement with scipy.spatial.distance.cdist
def pocket_overlap(structure, domain_chain, motif_chain):
    """
    structure1, structure2: Biopython models
    chain1, chain2 : Strings -> Chain IDs
    """

    model = structure[0]

    pocket_res_ids = dict()
    motif_res_ids = dict()
    full_interaction = dict()

    problem_atoms = set()
    problem_residues = set()

    # Filter out hetatoms
    domain_residues = [x for x in model[domain_chain].get_residues() if x.id[0] == " "]
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

        for (pos1, atom1), (pos2, atom2) in product(
            enumerate(res1.get_atoms()), enumerate(res2.get_atoms())
        ):
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
            except KeyError as e:
                problem_atoms.add(atom1.id)
                continue
            try:
                vdw2 = VDW_RADII[atom2.id[0]]
            except KeyError as e:
                problem_atoms.add(atom2.id)
                continue

            vdw_range = vdw1 + vdw2
            overlap = vdw_range - distance
            if overlap > -0.4:
                (
                    full_interaction.setdefault(res1.id[1], dict()).setdefault(
                        res2.id[1], set()
                    )
                ).add((pos1 not in backbone1, pos2 not in backbone2))

                pocket_res_ids.setdefault(res1.id[1], False)
                if pos1 not in backbone1:
                    pocket_res_ids[res1.id[1]] = True
                motif_res_ids.setdefault(res2.id[1], False)
                if pos2 not in backbone1:
                    motif_res_ids[res2.id[1]] = True

    if len(problem_atoms) > 0:
        logging.warning(f"No vdw radius for {list(problem_atoms)} in {structure.id}")
    if len(problem_residues) > 0:
        logging.warning(f"No single AA code for {problem_residues}: {structure.id}")

    # Dict for mapping residue id to sequence position
    res_id_to_pos = {}
    res_pos_coords = {}
    for i, res in enumerate(domain_residues):
        atoms = list(res.get_atoms())
        if len(atoms) > 1 and atoms[1].id == "CA":
            res_id_to_pos[res.id[1]] = i
            res_pos_coords[i] = atoms[1].coord.tolist()


    # mapping pocket ids to sequence position for foldseek
    pocket_res_pos = {res_id_to_pos[k]: v for k, v in pocket_res_ids.items()}

    full_interaction = sets_to_lists(full_interaction)  # sets are not JSON serializable
    pocket = {
        "pocket_exists": len(pocket_res_ids) > 0,
        "pocket_res_ids": pocket_res_ids,
        "pocket_res_pos": pocket_res_pos,
        "res_id_to_pos": res_id_to_pos,
        "pocket_to_motif_sidechain_overlap": full_interaction,
        "res_pos_coords": res_pos_coords,
    }


    return pocket, problem_atoms, problem_residues

def compare_pockets(
    alignment_df,
    pocket_dict,
    blosum_path=r"/home/data/motif_aligner/blosum62.bla",
    alphafold=False,
    alphafold_dir=None,
):
    """
    Compare two pockets based on foldseek alignment

    output = [i1, i2, ]
    """

    blosum_similarity_matrix = read_blast_similarity_matrix(blosum_path)

    domain_pocket_dict = defaultdict(dict)
    for k, v in pocket_dict.items():
        domain_pocket_dict[k[:6]][k[7]] = v

    # Setting up vars for use later
    existing_calcs = set()
    output_rows = []
    sup = SVDSuperimposer()

    for row in tqdm(alignment_df.itertuples(index=False)):

        domain_1 = row[0]
        domain_2 = row[1]
        try:
            # % non-gaps to gaps each way, similarity

            p1_seq_to_aln = {}
            i = 0
            for j, res in enumerate(row[12]):
                if res != "-":
                    p1_seq_to_aln[i] = j
                    i += 1
            p1_aln_to_seq = {v: k for k, v in p1_seq_to_aln.items()}

            p2_seq_to_aln = {}
            i = 0
            for j, res in enumerate(row[13]):
                if res != "-":
                    p2_seq_to_aln[i] = j
                    i += 1
            p2_aln_to_seq = {v: k for k, v in p2_seq_to_aln.items()}
        
            # Check that both pockets are loaded:
            pockets_1 = domain_pocket_dict.get(domain_1)
            if not pockets_1:
                continue
            if alphafold:
                pockets_2 = {
                    'A': {
                        'pocket_exists':True,
                        'pocket_res_pos':dict(zip(range(row[8] - 1, row[9]), repeat(True))),
                    }
                }
            else:
                pockets_2 = domain_pocket_dict.get(domain_2)
                if not pockets_2:
                    continue

            # Iterating through aligned pairs
            for motif_1, motif_2 in product(pockets_1.keys(), pockets_2.keys()):
                # Defining interaction names
                interaction_1 = domain_1 + "_" + motif_1
                if alphafold:
                    interaction_2 = domain_2
                else:
                    interaction_2 = domain_2 + "_" + motif_2

                # No self comparisons
                if interaction_1 == interaction_2:
                    continue

                # Checking for A-B comparison if B-A has already been calculated
                if (interaction_2, interaction_1) in existing_calcs:
                    continue
                existing_calcs.add((interaction_1, interaction_2))

                # Starting the output list
                output = {
                    "pocket_1": domain_1 + "_" + motif_1,
                    "pocket_2": domain_2 + "_" + motif_2,
                    "evalue": row[10],
                    "lddt": row[11],
                }

                # Getting the pockets
                p1 = pockets_1.get(motif_1)
                p2 = pockets_2.get(motif_2)
                if not p1["pocket_exists"] or not p2["pocket_exists"]:
                    continue

                ## Calculated Metrics

                # subtract start of alignment       
                p1_adj = 1 - row[6]
                p2_adj = 1 - row[8]


                p1_adjusted_start = {
                    (int(pos) + p1_adj): sidechain
                    for pos, sidechain in p1["pocket_res_pos"].items()
                }
                p2_adjusted_start = {
                    (int(pos) + p2_adj): sidechain
                    for pos, sidechain in p2["pocket_res_pos"].items()
                }
                p1_in_aln_region = {  # only the indices that are in the alignment region
                    k: v
                    for k, v in p1_adjusted_start.items()
                    if (
                        -1 < k and k <= (row[7] - row[6])
                    )  # row[7] - row[6] is the length of the aligned region
                }
                p2_in_aln_region = {
                    k: v
                    for k, v in p2_adjusted_start.items()
                    if (-1 < k and k <= (row[9] - row[8]))
                }
                p1_percent_in_aln = len(p1_in_aln_region) / len(p1_adjusted_start)
                p2_percent_in_aln = len(p2_in_aln_region) / len(p2_adjusted_start)
                
                output["pocket_1_len"] = len(p1["pocket_res_pos"])
                output["pocket_2_len"] = len(p2["pocket_res_pos"])
                output["pocket_1_pct_aln"] = p1_percent_in_aln
                output["pocket_2_pct_aln"] = p2_percent_in_aln

                """
                Finished Here
                1) need to make sure AF residue is actually aligned to reference pocket
                """

                # if either pocket isn't in the aligned region, end here
                if min(p1_percent_in_aln, p2_percent_in_aln) == 0:
                    output_rows.append(output)
                    continue

                # map sequence position to alignment position
                p1_aln_pos = {p1_seq_to_aln[k]: v for k, v in p1_in_aln_region.items()}
                p2_aln_pos = {p2_seq_to_aln[k]: v for k, v in p2_in_aln_region.items()}

                # count the overlapping residues
                overlapping_residues = [x for x in p1_aln_pos if x in p2_aln_pos]
                output["overlap_res"] = overlapping_residues
                output["overlap_count"] = len(overlapping_residues)

                # mapping common coords back to pocket 1 and 2 position coords
                p1_overlap_pos = [
                    p1_aln_to_seq[x] - p1_adj for x in overlapping_residues
                ]
                p2_overlap_pos = [
                    p2_aln_to_seq[x] - p2_adj for x in overlapping_residues
                ]
                
                x = np.array([p1["res_pos_coords"][str(x)] for x in p1_overlap_pos])
                y = np.array([p2["res_pos_coords"][str(x)] for x in p2_overlap_pos])

                if len(overlapping_residues) > 2:
                    sup.set(x, y)
                    sup.run()
                    u, t = sup.get_rotran()
                    output["p2_to_p1_u"] = u.flatten().tolist()
                    output["p2_to_p1_t"] = t.tolist()
                    sup.set(y, x)
                    sup.run()
                    u, t = sup.get_rotran()
                    output["p1_to_p2_u"] = u.flatten().tolist()
                    output["p1_to_p2_t"] = t.tolist()
                    output["rmsd"] = sup.get_rms()
                    output["rmsd_score"] = -np.log2(sup.get_rms()) + 1
                

                # If no overlapping resides, end here
                if len(overlapping_residues) == 0:
                    output_rows.append(output)
                    continue

                # Percentage overlap
                p1_pct_overlap = len(overlapping_residues) / len(p1_adjusted_start)
                p2_pct_overlap = len(overlapping_residues) / len(p2_adjusted_start)
                output["pocket_1_pct_overlap"] = p1_pct_overlap
                output["pocket_2_pct_overlap"] = p2_pct_overlap
                output["min_pct_overlap"] = min(p1_pct_overlap, p2_pct_overlap)
                output["max_pct_overlap"] = max(p1_pct_overlap, p2_pct_overlap)

                # Aligned residues as a sequence
                p1_aln_seq = "".join([row[12][x] for x in overlapping_residues])
                p2_aln_seq = "".join([row[13][x] for x in overlapping_residues])
                output["pocket_1_aln_seq"] = p1_aln_seq
                output["pocket_2_aln_seq"] = p2_aln_seq

                # Identity
                overlap_identity = sum(map(str.__eq__, p1_aln_seq, p2_aln_seq)) / len(
                    overlapping_residues
                )
                output["overlap_identity"] = overlap_identity

                # BLOSUM62 similarity
                similarity = binary_similarity(
                    p1_aln_seq, p2_aln_seq, blosum_similarity_matrix
                )
                output["overlap_similarity_binary"] = similarity

                similarity_1_2 = full_similarity(
                    p1_aln_seq, p2_aln_seq, blosum_similarity_matrix
                )
                similarity_2_1 = full_similarity(
                    p2_aln_seq, p1_aln_seq, blosum_similarity_matrix
                )
                output["overlap_similarity_1_2"] = similarity_1_2
                output["overlap_similarity_2_1"] = similarity_2_1
                output["min_overlap_similarity"] = min(similarity_1_2, similarity_2_1)
                output["max_overlap_similarity"] = max(similarity_1_2, similarity_2_1)

                # Sidechain interactor conservation p1
                p1_sidechain_pos = [x for x in overlapping_residues if p1_aln_pos[x]]
                p1_seq_p1_sc_contact = "".join([row[12][x] for x in p1_sidechain_pos])
                p2_seq_p1_sc_contact = "".join([row[13][x] for x in p1_sidechain_pos])
                output["p1_seq_p1_sc_contact"] = p1_seq_p1_sc_contact
                output["p2_seq_p1_sc_contact"] = p2_seq_p1_sc_contact
                p1_sc_contact_count = len(p1_sidechain_pos)
                output["p1_sc_contact_count"] = p1_sc_contact_count
                if p1_sc_contact_count > 0:
                    p1_sc_similarity = full_similarity(
                        p1_seq_p1_sc_contact,
                        p2_seq_p1_sc_contact,
                        blosum_similarity_matrix,
                    )
                    output["p1_sc_similarity"] = p1_sc_similarity

                # Sidechain interactor conservation p2
                p2_sidechain_pos = [x for x in overlapping_residues if p2_aln_pos[x]]
                p1_seq_p2_sc_contact = "".join([row[12][x] for x in p2_sidechain_pos])
                p2_seq_p2_sc_contact = "".join([row[13][x] for x in p2_sidechain_pos])
                output["p1_seq_p2_sc_contact"] = p1_seq_p2_sc_contact
                output["p2_seq_p2_sc_contact"] = p2_seq_p2_sc_contact
                p2_sc_contact_count = len(p2_sidechain_pos)
                output["p2_sc_binder_count"] = p2_sc_contact_count
                if p2_sc_contact_count > 0:
                    p2_sc_similarity = full_similarity(
                        p2_seq_p2_sc_contact,
                        p1_seq_p2_sc_contact,
                        blosum_similarity_matrix,
                    )
                    output["p2_sc_similarity"] = p2_sc_similarity

                output["min_sc_similarity"] = min(p1_sc_similarity, p2_sc_similarity)

                output_rows.append(output)
        except:
            output_rows.append(output)
            logging.exception(
                f"Pocket Comparison Failure: {domain_1}_{motif_1}, {domain_2}_{motif_2}"
            )

    return pd.DataFrame.from_dict(output_rows)


def binary_similarity(seqA, seqB, similarity_matrix):
    """
    Similarity of A->B and B->A
    score = 1 if full similarity score > 0 else score = 0
    normalized to the length of the seqence
    """
    seqA = seqA.replace("U", "X").upper()
    seqB = seqB.replace("U", "X").upper()

    similarity = [
        (similarity_matrix[x][y] > 0) for x, y in zip(seqA, seqB)
    ]  # True or False
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
        if header == None:
            if delimiter == " ":
                header = (
                    line.strip().split()
                )  # splits on whitespace and discards empty results
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

