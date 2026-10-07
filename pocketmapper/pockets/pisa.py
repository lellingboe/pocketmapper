"""
Reading of cached PISA interface data into pocket residue sets.

Consumes the JSON files `PisaDownloader` writes and turns one interface into the Pocket the rest of
the pipeline expects. This is the `pisa` pocket method, available for PDB entries only --
AlphaFold models and local files have no PISA data.

`pisa_pockets` is the method's builder: download, parse, then add coordinates from the structure.
`pisa_cache_layout` owns the cache layout under `pocket_dir/pisa/`.
"""

import json
import logging
import os

from tqdm import tqdm

from pocketmapper.constants import DEFAULT_PISA_SOURCE
from pocketmapper.downloads.pisa_downloader import PisaDownloader
from pocketmapper.lib import one_letter_code
from pocketmapper.lib import split_chain_info
from pocketmapper.pockets.pocket import Pocket
from pocketmapper.pockets.pocket import PocketResidue
from pocketmapper.pockets.structure import parse_pocket_from_struct

logger = logging.getLogger(__name__)


class PisaParser:
    """
    Turns cached PISA interfaces into Pockets.
    """

    def load_interfaces(self, pdb_id, in_dir):
        """
        Load the cached interface file for a PDB entry, or None if there isn't one.

        Args:
            pdb_id (str): PDB entry to load, in any case. Tried verbatim, then lower-cased.
            in_dir (str): Directory of parsed interface files.

        Returns:
            dict: The entry's interfaces keyed by sorted, comma-joined chain pair ("A,B-2"), or None
                if not cached.
        """
        # Cached files are named by the lower-cased code; pdb_id may be in any case.
        for candidate in (pdb_id, pdb_id.lower()):
            in_path = os.path.join(in_dir, f"{candidate}.json")
            if os.path.exists(in_path):
                with open(in_path, "r") as f:
                    return json.load(f)
        return None

    def get_interface_partners(self, pdb_id, chain_id, in_dir):
        """
        List the chains a given chain shares a PISA interface with.

        Args:
            pdb_id (str): PDB entry the chain belongs to.
            chain_id (str): Chain whose partners are wanted.
            in_dir (str): Directory of parsed interface files.

        Returns:
            list[str]: Partner chain ids, in file order; a homodimer yields `chain_id` itself. Empty if
                there is no data for this entry or the chain takes part in no interface.
        """
        log_extra = {"stage": "Calculating Pockets"}
        pisa_data = self.load_interfaces(pdb_id, in_dir)
        if pisa_data is None:
            logger.debug(f"Could not load PISA data for {pdb_id}", extra=log_extra)
            return []

        # Keys are the interface's two chain ids, sorted and comma-joined ("A,B-2"); a partner is the
        # other half of every key the chain appears in.
        partners = []
        for interface_chains in pisa_data:
            chains = interface_chains.split(",")
            if chain_id not in chains:
                continue
            chains.remove(chain_id)
            partner = chains[0]
            if partner not in partners:
                partners.append(partner)
        if not partners:
            logger.debug(f"No PISA interface involving chain {chain_id} of {pdb_id}", extra=log_extra)
        return partners

    def get_pockets_from_records(self, records, in_dir):
        """
        Build a Pocket per record from the residues its chain contributes to any bond of its interface.

        Args:
            records (list): QTRecord dicts; reads `struct_info`, `chain_info` and `pocket_id`.
            in_dir (str): Directory of parsed interface files.

        Returns:
            dict: pocket_id -> Pocket with residue codes but no coordinates: `ca_sequence`,
                `has_coords` and every residue's `seq_pos` and `ca_coords` are left at their defaults.
                A record whose entry, interface or domain chain cannot be resolved, or whose interface
                does not have exactly two molecules, is skipped with a warning.
        """
        log_extra = {"stage": "Calculating Pockets"}
        bond_types = ["hydrogen_bonds", "salt_bridges", "disulfide_bonds", "covalent_bonds", "other_bonds"]
        pockets = {}
        for record in records:
            pdb_id = record["struct_info"]

            # Load the entry's interfaces
            pisa_data = self.load_interfaces(pdb_id, in_dir)
            if pisa_data is None:
                logger.warning(f"Could not load PISA data for {pdb_id}", extra=log_extra)
                continue

            # Extract the interface. A record with no partner chain is skipped, not left to raise in sorted().
            domain_chain, motif_chain = split_chain_info(record["chain_info"])
            if motif_chain is None:
                logger.warning(f"No partner chain in chain_info '{record['chain_info']}' for {pdb_id}", extra=log_extra)
                continue
            interface_chains = ",".join(sorted([domain_chain, motif_chain]))
            if interface_chains not in pisa_data:
                logger.warning(f"No PISA data for {pdb_id} interface {interface_chains}", extra=log_extra)
                continue
            pisa_data = pisa_data[interface_chains]

            # A pocket is defined against a single partner
            if not len(pisa_data["molecules"]) == 2:
                logger.warning(f"More than two molecules in {pdb_id} interface {interface_chains}", extra=log_extra)
                continue

            # Find the molecule id of the domain chain
            pocket_mol_id = None
            for mol in pisa_data["molecules"]:
                if mol["chain_id"] == domain_chain:
                    pocket_mol_id = mol["molecule_id"]
                    break
            if pocket_mol_id is None:
                logger.warning(f"Could not find domain chain in {pdb_id} interface {interface_chains}", extra=log_extra)
                continue

            pocket = Pocket()

            # Collect the domain chain's residues across all bond types, keyed by author seqid as a string
            all_res_auth_ids = set()
            for bond_type in bond_types:
                bonds_dict = pisa_data[bond_type]
                res_auth_ids = bonds_dict[f"atom_site_{pocket_mol_id}_seq_nums"]
                all_res_auth_ids.update(res_auth_ids)
                for i, res_auth_id in enumerate(res_auth_ids):
                    res_code = bonds_dict[f"atom_site_{pocket_mol_id}_residues"][i]
                    pocket.residues[str(res_auth_id)] = PocketResidue(
                        res_code=res_code,
                        res_code_single=one_letter_code(res_code),
                        uniprot_pos=bonds_dict[f"atom_site_{pocket_mol_id}_unp_nums"][i],
                    )

            # Sorted numerically: walking the bond types leaves the residues in arbitrary order
            pocket.res_auth_ids = [str(x) for x in sorted(int(x) for x in all_res_auth_ids)]
            pocket.pocket_exists = len(pocket.res_auth_ids) > 0

            pockets[record["pocket_id"]] = pocket

        return pockets


def pisa_cache_layout(pocket_dir):
    """
    The PISA cache's locations under a pocket cache directory.

    Args:
        pocket_dir (str): Pocket cache directory.

    Returns:
        dict: "summary_dir", "asm_dir", "interface_dir" and "error_path", each under `pocket_dir/pisa/`.
            `interface_dir` holds the per-entry interface files `PisaParser` reads.
    """
    pisa_dir = os.path.join(pocket_dir, "pisa")
    return {
        "summary_dir": os.path.join(pisa_dir, "summaries"),
        "asm_dir": os.path.join(pisa_dir, "assemblies"),
        "interface_dir": os.path.join(pisa_dir, "interface_pairs"),
        "error_path": os.path.join(pisa_dir, "errors.json"),
    }


def download_pisa_interfaces(pdb_list, pocket_dir, pisa_source):
    """
    Populate the PISA interface cache under `pocket_dir/pisa/` for a list of PDB entries.

    When any entry fails, `pisa/errors.json` is overwritten with the failures.

    Args:
        pdb_list (list): PDB codes, in any case.
        pocket_dir (str): Pocket cache directory.
        pisa_source (str): Where assembly interfaces are fetched from: "ftp" or "api".

    Returns:
        str: The directory of per-entry interface files, for `PisaParser`.
    """
    layout = pisa_cache_layout(pocket_dir)
    PisaDownloader(source=pisa_source).download_missing_interfaces(pdb_list=pdb_list, **layout)
    return layout["interface_dir"]


def pisa_pockets(records, pocket_dir, pisa_source=DEFAULT_PISA_SOURCE, download=True):
    """
    Build a Pocket per record from the PDBe PISA interface it names, with coordinates from its structure.

    Unless `download` is False, downloads any PISA files not already cached under `pocket_dir/pisa/`.
    With the "api" source that is one paced request per assembly, so an uncached list can take a long
    time. When any entry fails, `pisa/errors.json` is overwritten with the failures.

    Args:
        records (list): QTRecord dicts with `pocket_method == "pisa"`.
        pocket_dir (str): Pocket cache directory.
        pisa_source (str): Where assembly interfaces are fetched from: "ftp" or "api".
            Defaults to DEFAULT_PISA_SOURCE.
        download (bool): False reads only what is already cached. Defaults to True.

    Returns:
        dict: pocket_id -> Pocket. A record whose interface cannot be resolved is skipped with a
            warning; one whose structure or chain cannot be read maps to None.
    """
    log_extra = {"stage": "Retrieving PISA Pockets"}

    pdb_list = list(dict.fromkeys(record["struct_info"] for record in records))
    logger.debug(f"PDBs for which to retrieve PISA pockets: {pdb_list}", extra=log_extra)
    if download:
        interface_dir = download_pisa_interfaces(pdb_list, pocket_dir, pisa_source)
    else:
        interface_dir = pisa_cache_layout(pocket_dir)["interface_dir"]

    pockets = PisaParser().get_pockets_from_records(records=records, in_dir=interface_dir)
    logger.debug(f"PISA pockets before coordinates: {pockets}", extra=log_extra)

    # PisaParser gives residue ids only; add seq_pos and CA coordinates to the same Pocket
    for record in tqdm(records):
        if record["pocket_id"] in pockets:
            domain_chain, _ = split_chain_info(record["chain_info"])
            pockets[record["pocket_id"]] = parse_pocket_from_struct(
                struct=record["struct_path"],
                chain_id=domain_chain,
                pocket_residues=[int(x) for x in pockets[record["pocket_id"]].res_auth_ids],
                pocket=pockets[record["pocket_id"]],
            )
    return pockets
