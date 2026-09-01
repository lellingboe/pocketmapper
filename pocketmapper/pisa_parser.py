"""
Reading of cached PISA interface data into pocket residue sets.

Consumes the JSON files `PisaDownloader` writes and turns one interface into the Pocket the rest of
the pipeline expects. This is the `pisa` pocket method, available for PDB entries only --
AlphaFold models and local files have no PISA data.
"""

import os
import json
import logging
from pocketmapper.constants import SINGLE_AA_CODE
from pocketmapper.pocket import Pocket, PocketResidue


class PisaParser:
    """
    Turns cached PISA interfaces into Pockets.
    """

    def __init__(self):
        """
        Initialise the parser. State lives in the cache directory, not on the instance.
        """
        logging.getLogger(__name__)

    def _load_interfaces(self, pdb_id, in_dir):
        """
        Load the cached interface file for a PDB entry, or None if there isn't one.

        PisaDownloader writes these files under a lower-cased PDB code while callers hold whatever case
        the user typed, so the name is tried verbatim first and then lower-cased. Without that an
        upper-cased input finds nothing on a case-sensitive filesystem.

        Args:
            pdb_id (str): PDB entry to load, in any case.
            in_dir (str): Directory of parsed interface files.

        Returns:
            dict: The entry's interfaces keyed by sorted chain pair, or None if not cached.
        """
        for candidate in (pdb_id, pdb_id.lower()):
            in_path = os.path.join(in_dir, f"{candidate}.json")
            if os.path.exists(in_path):
                with open(in_path, "r") as f:
                    return json.load(f)
        return None

    def get_interface_partners(self, pdb_id, chain_id, in_dir):
        """
        List the chains a given chain shares a PISA interface with.

        Interface files are keyed by the two chain ids of the interface, sorted and concatenated
        (e.g. "BF"), so a chain's partners are the other half of every key it appears in. A
        homodimer key ("BB") yields the chain itself.

        Args:
            pdb_id (str): PDB entry the chain belongs to.
            chain_id (str): Chain whose partners are wanted.
            in_dir (str): Directory of parsed interface files written by PisaDownloader.

        Returns:
            list[str]: Partner chain ids, in file order. Empty if there is no data for this entry
                or the chain takes part in no interface.
        """
        stage = {"stage": "Calculating Pockets"}
        pisa_data = self._load_interfaces(pdb_id, in_dir)
        if pisa_data is None:
            logging.debug(f"Could not load PISA data for {pdb_id}", extra=stage)
            return []

        partners = []
        for interface_chains in pisa_data:
            if chain_id not in interface_chains:
                continue
            partner = interface_chains.replace(chain_id, "", 1)
            if partner and partner not in partners:
                partners.append(partner)
        if not partners:
            logging.debug(f"No PISA interface involving chain {chain_id} of {pdb_id}", extra=stage)
        return partners

    def get_pockets_from_records(self, records, in_dir):
        """
        Build a Pocket per record from its PISA interface.

        The pocket is the set of residues on the record's own chain that take part in any bond of the
        interface, across all five bond types. Records whose entry, interface or domain chain cannot be
        resolved are logged and skipped, so the result may be smaller than `records`.

        Only the interface with exactly two molecules is usable, since the pocket is defined against a
        single partner.

        Args:
            records (list): QTRecord dicts; reads `struct_info`, `chain_info` and `pocket_id`.
            in_dir (str): Directory of parsed interface files written by PisaDownloader.

        Returns:
            dict: pocket_id -> Pocket, carrying `res_auth_ids` and one residue per interface residue.
                Lacks the CA data `pocket_parser.parse_pocket_from_struct` adds later, so `ca_sequence`,
                `has_coords` and every residue's `seq_pos`/`ca_coords` are still at their defaults.
        """
        stage = {"stage": "Calculating Pockets"}
        bond_types = ["hydrogen_bonds", "salt_bridges", "disulfide_bonds", "covalent_bonds", "other_bonds"]
        pockets = {}
        for record in records:
            pdb_id = record["struct_info"]

            # Loading PISA pocket file
            pisa_data = self._load_interfaces(pdb_id, in_dir)
            if pisa_data is None:
                logging.warning(f"Could not load PISA data for {pdb_id}", extra=stage)
                continue

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
            pocket = Pocket()

            # Getting the pocket residues. Residues are keyed as strings to match the author seqids
            # parse_pocket_from_struct will later use when it adds coordinates to this same Pocket.
            all_res_auth_ids = set()
            for bond_type in bond_types:
                bonds_dict = pisa_data[bond_type]
                res_auth_ids = bonds_dict[f"atom_site_{pocket_mol_id}_seq_nums"]
                all_res_auth_ids.update(res_auth_ids)
                for i, res_auth_id in enumerate(res_auth_ids):
                    res_code = bonds_dict[f"atom_site_{pocket_mol_id}_residues"][i]
                    pocket.residues[str(res_auth_id)] = PocketResidue(
                        res_code=res_code,
                        res_code_single=SINGLE_AA_CODE.get(res_code, "X"),
                        uniprot_pos=bonds_dict[f"atom_site_{pocket_mol_id}_unp_nums"][i],
                    )

            # Sorted numerically: the bond types are walked in list order, so insertion order into
            # residues is arbitrary, and res_auth_ids is the ordering the comparison relies on.
            pocket.res_auth_ids = [str(x) for x in sorted(int(x) for x in all_res_auth_ids)]
            pocket.pocket_exists = len(pocket.res_auth_ids) > 0

            pockets[record["pocket_id"]] = pocket

        return pockets
