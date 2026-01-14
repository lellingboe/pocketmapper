import gzip
import shutil
import os
import logging
from urllib.request import urlcleanup, urlretrieve
from copy import deepcopy

from Bio.SVDSuperimposer import SVDSuperimposer
from Bio.PDB import MMCIFParser
from glob import glob
from concurrent.futures import ThreadPoolExecutor
from itertools import repeat, product
from tqdm import tqdm
import json
from collections import defaultdict
import pandas as pd
import gemmi
import pisa
from numpy import array
from numpy import linalg as LA

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
    "PTR": "Y",  # phospho
    "MSE": "M",  # selenomethionine
}
TRIPLE_AA_CODE = defaultdict(list)
for k, v in SINGLE_AA_CODE.items():
    TRIPLE_AA_CODE[v].append(k)
VDW_RADII = {"C": 1.88, "N": 1.64, "O": 1.46, "S": 1.77, "P": 1.87, "H": 1.0}
# https://www.cgl.ucsf.edu/chimerax/docs/user/commands/clashes.html


def get_mmcif(pdb_code, out_dir, cache):
    stage = {"stage": "Downloading PDB File"}
    pdb_code = pdb_code.lower()
    out_fname = os.path.join(out_dir, f"{pdb_code}.cif.gz")
    if not (out_fname in cache):
        url = f"https://files.wwpdb.org/pub/pdb/data/structures/divided/mmCIF/{pdb_code[1:3]}/{pdb_code}.cif.gz"
        try:
            urlcleanup()
            urlretrieve(url, out_fname)
        except OSError:
            logging.warning(f"get_mmcif: Could not download {pdb_code}", extra=stage)
            return (pdb_code, False)
        except Exception:
            logging.warning(f"Atypical issue when downloading {pdb_code}", extra=stage)
            return (pdb_code, False)
        # else:
        #    with gzip.open(gz_fname, "rb") as gz:
        #        with open(out_fname, "wb") as out:
        #            out.writelines(gz)
        #    os.remove(gz_fname)
    return (pdb_code, True)


def get_mmcifs(pdb_list, out_dir):
    cache = glob(os.path.join(out_dir, "*.cif.gz"))
    with ThreadPoolExecutor(max_workers=100) as e:
        result = e.map(
            get_mmcif,
            pdb_list,
            repeat(out_dir),
            repeat(cache),
        )
    return {x.upper(): y for x, y in result}


def pdb_preprocessing_gemmi(df, ref_dir, cache_dir, target_dir, query_dir):
    """
    queries: a lit of tuple of the form (pdb_id, domain_chains, motif_chains)
        all tuple elements are strings

    writes out .pdb files to self.pdb_directory
    """
    status_dict = {}
    stage = {"stage": "Dividing structures"}

    for i, row in tqdm(df.iterrows()):
        try:
            cache_domain_path = os.path.join(cache_dir, f"{row.pdb_domain}.cif")
            cache_domain_path_gz = cache_domain_path + ".gz"
            cache_motif_path = os.path.join(cache_dir, f"{row.pdb_domain_motif}.cif")
            cache_motif_path_gz = cache_motif_path + ".gz"

            # Ensuring divided structure is in the cache directory
            if not os.path.exists(cache_domain_path_gz):
                ref_path = os.path.join(ref_dir, f"{row.interaction_pdb}.cif.gz")
                st = gemmi.read_structure(ref_path, format=gemmi.CoorFormat.Mmcif)

                # Taking first model and deleting the rest
                del st[1:]
                model = st[0]

                # verify structure contains all interaction chains
                interaction_chains = list(row.domain_chain + row.motif_chain)
                model_chains = set([chain.name for chain in model])
                if not set(interaction_chains).issubset(model_chains):
                    msg = f"Preprocessing: {row.interaction_pdb}"
                    "does not contain all interaction chains {interaction_chains}"
                    logging.warning(
                        msg,
                        extra=stage,
                    )
                    status_dict[f"{row.pdb_domain_motif}"] = False
                    continue

                # Detaching all non interaction chains
                for chain_id in model_chains:
                    if chain_id not in interaction_chains:
                        del model[chain_id]

                # Output the domain and motif pdb file
                groups = gemmi.MmcifOutputGroups(False, atoms=True, group_pdb=True)
                st.make_mmcif_document(groups).write_file(cache_motif_path)
                with open(cache_motif_path, "rb") as f_in:
                    with gzip.open(cache_motif_path_gz, "wb") as f_out:
                        shutil.copyfileobj(f_in, f_out)
                os.remove(cache_motif_path)

                # output the domain pdb file
                del model[row.motif_chain]
                st.make_mmcif_document(groups).write_file(cache_domain_path)
                with open(cache_domain_path, "rb") as f_in:
                    with gzip.open(cache_domain_path_gz, "wb") as f_out:
                        shutil.copyfileobj(f_in, f_out)
                os.remove(cache_domain_path)

                status_dict[f"{row.pdb_domain_motif}"] = True

            # Copying divides structures into forders for foldseek
            if row.type == "query":
                out_dir = query_dir
            elif row.type == "target":
                out_dir = target_dir
            else:
                logging.warning(f"row with {row.interaction_pdb} is neither target nor query", extra=stage)
                continue

            domain_out = os.path.join(out_dir, f"{row.pdb_domain}.cif.gz")
            motif_out = os.path.join(out_dir, f"{row.pdb_domain_motif}.cif.gz")
            shutil.copyfile(cache_domain_path_gz, domain_out)
            shutil.copyfile(cache_motif_path_gz, motif_out)
            status_dict[f"{row.pdb_domain_motif}"] = True

        except Exception as e:
            logging.warning(f"Could not divide {row.interaction_pdb}", extra=stage)
            logging.debug("Exception info", exc_info=e, extra="stage")
            status_dict[f"{row.pdb_domain_motif}"] = False

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


def get_pisa_pockets(df, in_dir, out_dir):
    stage = {"stage": "Calculating Pockets"}
    """Takes in a path to a pdb file"""
    bond_types = ["hydrogen_bonds", "salt_bridges", "disulfide_bonds", "covalent_bonds", "other_bonds"]
    pockets = {}
    for i, row in tqdm(df.iterrows()):
        pdb_id = row.interaction_pdb.lower()

        # Loading PISA pocket file
        in_path = os.path.join(in_dir, f"{pdb_id}.json")
        if not os.path.exists(in_path):
            logging.warning(f"Could not load PISA data for {pdb_id}", extra=stage)
            continue
        with open(in_path, "r") as f:
            pisa_data = json.load(f)

        # Extracting the relevant interface
        interface_chains = "".join(sorted(row.domain_chain + row.motif_chain))
        if interface_chains not in pisa_data:
            logging.warning(f"No PISA data for {pdb_id} interface {interface_chains}", extra=stage)
            continue
        pisa_data = pisa_data[interface_chains]

        # Checking the interfaces features 2 molecules
        if not len(pisa_data["molecules"]) == 2:
            logging.warning(f"More than two molecules in {pdb_id} interface {interface_chains}", extra=stage)
            continue

        # Getting the molecule id for the domain chain
        pocket_mol_id = None
        for mol in pisa_data["molecules"]:
            if mol["chain_id"] == row.domain_chain:
                pocket_mol_id = mol["molecule_id"]
                break
        if pocket_mol_id is None:
            logging.warning(f"Could not find domain chain in {pdb_id} interface {interface_chains}", extra=stage)
            continue

        # Making output pocket
        pocket = {
            "res_auth_ids": set(),
            "id_pos_codes_match": True,
        }

        # Getting the pocket residues
        for bond_type in bond_types:
            bonds_dict = pisa_data[bond_type]
            res_auth_ids = bonds_dict[f"atom_site_{pocket_mol_id}_seq_nums"]
            pocket["res_auth_ids"].update(set(res_auth_ids))
            for i, res_auth_id in enumerate(res_auth_ids):
                if res_auth_id not in pocket:  # initializing dict
                    pocket[res_auth_id] = {}
                res_dict = {}
                res_dict["res_code"] = bonds_dict[f"atom_site_{pocket_mol_id}_residues"][i]
                res_dict["res_code_single"] = SINGLE_AA_CODE.get(res_dict["res_code"], "X")
                res_dict["uniprot_pos"] = bonds_dict[f"atom_site_{pocket_mol_id}_unp_nums"][i]

                pocket[res_auth_id] = res_dict

        sorted_res_auth_ids = [str(x) for x in sorted([int(x) for x in pocket["res_auth_ids"]])]
        pocket["res_auth_ids"] = sorted_res_auth_ids
        if len(pocket["res_auth_ids"]) > 0:
            pocket["pocket_exists"] = True

        # Making JSON serializable
        pocket = jsonify_dict(pocket)
        pockets[row.pdb_domain_motif] = pocket

    return pockets


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


def compare_pockets(
    alignment_df,
    pocket_dict,
    blosum_path=r"/home/data/motif_aligner/blosum62.bla",
    alphafold=False,
):
    """
    Compare two pockets based on foldseek alignment
    """

    stage = {"stage": "Pocket Comparison"}
    blosum_similarity_matrix = read_blast_similarity_matrix(blosum_path)

    unknown_ids = defaultdict(lambda: defaultdict(set))

    # {pdb_domain_motif: info} -> {pdb_domain: {motif1: info, motif2:info}}
    domain_pocket_dict = defaultdict(dict)
    for k, v in pocket_dict.items():
        domain_pocket_dict[k[:6]][k[7]] = v

    # Setting up vars for use later
    existing_calcs = set()
    output_rows = []
    sup = SVDSuperimposer()

    # TODO divide this into common things for each pocket and a cross-comparison
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
            # p1_aln_to_seq = {v: k for k, v in p1_seq_to_aln.items()}

            p2_seq_to_aln = {}
            i = 0
            for j, res in enumerate(row[13]):
                if res != "-":
                    p2_seq_to_aln[i] = j
                    i += 1
            # p2_aln_to_seq = {v: k for k, v in p2_seq_to_aln.items()}

            # Check that both pockets are loaded:
            pockets_1 = domain_pocket_dict.get(domain_1)
            if not pockets_1:
                logging.debug(f"domain:{domain_1}", extra=stage)
                continue
            if alphafold:
                pockets_2 = {
                    "A": {
                        "res_auth_ids": [str(k) for k in range(row[9])],
                        "id_pos_codes_match": True,
                        "pocket_exists": True,
                        "has_coords": False,
                    }
                }
                pockets_2["A"].update({str(k): {"seq_pos": k} for k in range(row[9])})
            else:
                pockets_2 = domain_pocket_dict.get(domain_2)
                if not pockets_2:
                    logging.debug(f"domain:{domain_2}", extra=stage)
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
                # if interaction_1 == interaction_2:
                #    continue

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
                p1 = deepcopy(pockets_1.get(motif_1))
                p2 = deepcopy(pockets_2.get(motif_2))
                if not p1["pocket_exists"] or not p2["pocket_exists"]:
                    continue

                # Calculated Metrics
                # pocket 1
                p1_adj = 1 - row[6]
                p1_fs_len = row[7] - row[6] + 1
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
                                ",".join([interaction_2, interaction_1, res])
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
                            debug_id = ",".join([interaction_1, interaction_2, res])
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
                f"Uncontrolled KeyError calculating {domain_1}_{motif_1} and {domain_2}_{motif_2}",
                extra={"stage": "Pocket Comparison"},
            )
        except Exception:
            logging.exception(
                f"Uncontrolled error calculating {domain_1}_{motif_1} and {domain_2}_{motif_2}",
                extra={"stage": "Pocket Comparison"},
            )
            exit(1)

    return pd.DataFrame.from_dict(output_rows), unknown_ids


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
    downloader = pisa.PisaDownloader()
    downloader.get_interfaces(pdb_list, summary_dir, assembly_dir, interface_dir)
