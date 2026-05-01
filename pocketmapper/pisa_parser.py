import os
import json
import logging
from pocketmapper.lib import jsonify_dict
from pocketmapper.constants import SINGLE_AA_CODE


class PisaParser:
    def __init__(self):
        logging.getLogger(__name__)

    def get_pockets_from_records(self, records, in_dir):
        stage = {"stage": "Calculating Pockets"}
        """Takes in a path to a pdb file"""
        bond_types = ["hydrogen_bonds", "salt_bridges", "disulfide_bonds", "covalent_bonds", "other_bonds"]
        pockets = {}
        for record in records:
            pdb_id = record["struct_info"]

            # Loading PISA pocket file
            in_path = os.path.join(in_dir, f"{pdb_id}.json")
            if not os.path.exists(in_path):
                logging.warning(f"Could not load PISA data for {pdb_id}", extra=stage)
                continue
            with open(in_path, "r") as f:
                pisa_data = json.load(f)

            # Extracting the relevant interface
            interface_chains = "".join(sorted(record["chain_info"].split("_")))
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
                if mol["chain_id"] == record["chain_info"][0]:  # Assuming first chain is the domain chain
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
            pockets[record["pocket_id"]] = pocket

        return pockets
