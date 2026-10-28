"""
Reading of cached PISA interface data into pocket residue sets.

Consumes the JSON files `PisaDownloader` writes (`load_interfaces`) and turns one interface into the
Pocket the rest of the pipeline expects (`pisa_pockets_from_records`). This is the `pisa` pocket method,
available for PDB entries only -- AlphaFold models and local files have no PISA data. A pocket is the
domain chain's residues in any bond of the interface (`bonded_residues`) or, given a threshold, those
whose accessible surface the partner buries past it (`buried_residues`).

`pisa_pockets` is the method's builder: download, parse, then add coordinates from the structure.
`download_pisa_interfaces` places the PISA cache at `pocket_dir/pisa/`; `PisaDownloader` names what is
inside it.
"""

import json
import logging
import os

from pocketmapper.lib.constants import DEFAULT_PISA_SOURCE
from pocketmapper.lib.residues import one_letter_code, split_chain_info
from pocketmapper.pockets.pisa_downloader import PisaDownloader
from pocketmapper.pockets.pocket import Pocket, PocketResidue
from pocketmapper.pockets.structure import parse_pocket_from_struct, records_by_structure

logger = logging.getLogger(__name__)


def load_interfaces(pdb_id, in_dir):
    """
    Load the cached interface file for a PDB entry, or None if there isn't one.

    Args:
        pdb_id (str): PDB entry to load, in any case. Tried verbatim, then lower-cased.
        in_dir (str): Directory of parsed interface files.

    Returns:
        dict: The entry's interfaces keyed by sorted, comma-joined chain pair ("A,B-2"), or None if not
            cached.
    """
    # Cached files are named by the lower-cased code; pdb_id may be in any case.
    for candidate in (pdb_id, pdb_id.lower()):
        in_path = os.path.join(in_dir, f"{candidate}.json")
        if os.path.exists(in_path):
            with open(in_path, "r") as f:
                return json.load(f)
    return None


def interface_partners(interfaces, chain_id):
    """
    List the chains a given chain shares a PISA interface with.

    Args:
        interfaces (dict): One entry's interfaces, as `load_interfaces` returns them.
        chain_id (str): Chain whose partners are wanted.

    Returns:
        list: Partner chain ids, in file order; a homodimer yields `chain_id` itself. Empty if the chain
            takes part in no interface.
    """
    # Keys are the interface's two chain ids, sorted and comma-joined ("A,B-2"); a partner is the other
    # half of every key the chain appears in.
    partners = []
    for interface_chains in interfaces:
        chains = interface_chains.split(",")
        if chain_id not in chains:
            continue
        chains.remove(chain_id)
        partner = chains[0]
        if partner not in partners:
            partners.append(partner)
    return partners


def pisa_pockets_from_records(records, in_dir, bsa_threshold=None):
    """
    Build a Pocket per record from the residues its chain contributes to its interface.

    Each entry's interface file is read once, however many records name it.

    Args:
        records (list): QTRecord dicts; reads `struct_info`, `chain_info` and `pocket_id`.
        in_dir (str): Directory of parsed interface files.
        bsa_threshold (float or None, optional): Percent of a residue's accessible surface the partner
            must bury, strictly exceeded, for the residue to join (`buried_residues`). Defaults to None,
            which takes the residues in any bond instead (`bonded_residues`).

    Returns:
        dict: pocket_id -> Pocket with residue codes but no coordinates: `ca_sequence`, `has_coords` and
            every residue's `seq_pos` and `ca_coords` are left at their defaults. A record whose entry,
            interface or domain chain cannot be resolved, or whose interface does not have exactly two
            molecules, is left out with a warning.
    """
    log_extra = {"stage": "Retrieving PISA Pockets"}
    groups = {}
    for record in records:
        groups.setdefault(record["struct_info"], []).append(record)

    pockets = {}
    for pdb_id, group in groups.items():
        entry = load_interfaces(pdb_id, in_dir)
        if entry is None:
            logger.warning(f"Could not load PISA data for {pdb_id}, skipping its {len(group)} entries", extra=log_extra)
            continue

        for record in group:
            # A record with no partner chain is skipped, not left to raise in sorted()
            domain_chain, motif_chain = split_chain_info(record["chain_info"])
            if motif_chain is None:
                logger.warning(f"No partner chain in chain_info '{record['chain_info']}' for {pdb_id}", extra=log_extra)
                continue
            interface_chains = ",".join(sorted([domain_chain, motif_chain]))
            if interface_chains not in entry:
                logger.warning(f"No PISA data for {pdb_id} interface {interface_chains}", extra=log_extra)
                continue
            interface = entry[interface_chains]

            # A pocket is defined against a single partner
            if len(interface["molecules"]) != 2:
                logger.warning(
                    f"{pdb_id} interface {interface_chains} does not have exactly two molecules", extra=log_extra
                )
                continue

            pocket_mol = next((mol for mol in interface["molecules"] if mol["chain_id"] == domain_chain), None)
            if pocket_mol is None:
                logger.warning(f"Could not find domain chain in {pdb_id} interface {interface_chains}", extra=log_extra)
                continue

            pocket = Pocket()
            if bsa_threshold is None:
                pocket.residues = bonded_residues(interface, pocket_mol["molecule_id"])
            else:
                pocket.residues = buried_residues(pocket_mol, bsa_threshold)

            # Sorted numerically: neither rule leaves the residues in ascending order
            pocket.res_auth_ids = [str(x) for x in sorted(int(x) for x in pocket.residues)]
            pocket.pocket_exists = len(pocket.res_auth_ids) > 0

            pockets[record["pocket_id"]] = pocket

    return pockets


def bonded_residues(interface, mol_id):
    """
    Collect a molecule's residues across every bond type of a PISA interface.

    Args:
        interface (dict): One interface of a cached interface file.
        mol_id (str): `molecule_id` of the molecule whose residues are wanted.

    Returns:
        dict: Author seqid as a string -> PocketResidue with `res_code`, `res_code_single` and
            `uniprot_pos`, in bond order.
    """
    bond_types = ["hydrogen_bonds", "salt_bridges", "disulfide_bonds", "covalent_bonds", "other_bonds"]
    residues = {}
    for bond_type in bond_types:
        bonds_dict = interface[bond_type]
        res_auth_ids = bonds_dict[f"atom_site_{mol_id}_seq_nums"]
        for i, res_auth_id in enumerate(res_auth_ids):
            res_code = bonds_dict[f"atom_site_{mol_id}_residues"][i]
            residues[str(res_auth_id)] = PocketResidue(
                res_code=res_code,
                res_code_single=one_letter_code(res_code),
                uniprot_pos=bonds_dict[f"atom_site_{mol_id}_unp_nums"][i],
            )
    return residues


def buried_residues(molecule, bsa_threshold):
    """
    Collect a molecule's residues whose accessible surface its partner buries past a threshold.

    A residue joins when `100 * bsa / asa` is strictly greater than `bsa_threshold`; one with no
    accessible surface never joins.

    Args:
        molecule (dict): One entry of an interface's `molecules`.
        bsa_threshold (float): Percent buried a residue must exceed.

    Returns:
        dict: Author seqid as a string -> PocketResidue with `res_code` and `res_code_single`, in
            molecule order. `uniprot_pos` is None: the molecule lists carry no UniProt numbering.
    """
    residues = {}
    seen = set()
    rows = zip(
        molecule["residue_seq_ids"],
        molecule["residue_label_comp_ids"],
        molecule["accessible_surface_areas"],
        molecule["buried_surface_areas"],
    )
    for res_auth_id, res_code, asa, bsa in rows:
        # A seq id listed twice is one residue in alternative forms (microheterogeneity); its first row
        # decides, as the first conformer does wherever a chain is walked.
        if res_auth_id in seen:
            continue
        seen.add(res_auth_id)
        if asa > 0 and 100 * bsa / asa > bsa_threshold:
            residues[str(res_auth_id)] = PocketResidue(
                res_code=res_code,
                res_code_single=one_letter_code(res_code),
                uniprot_pos=None,
            )
    return residues


def download_pisa_interfaces(pdb_list, pocket_dir, pisa_source, download=True):
    """
    Populate the PISA interface cache under `pocket_dir/pisa/` for a list of PDB entries.

    When any entry fails, `pisa/errors.json` is overwritten with the failures. With `download` False,
    nothing is fetched or written, and an entry not cached is only counted in a warning.

    Args:
        pdb_list (list): PDB codes, in any case.
        pocket_dir (str): Pocket cache directory.
        pisa_source (str): Where assembly interfaces are fetched from: "ftp" or "api".
        download (bool, optional): Whether entries not cached may be downloaded. Defaults to True.

    Returns:
        str: The directory of per-entry interface files, for `load_interfaces`.
    """
    downloader = PisaDownloader(os.path.join(pocket_dir, "pisa"), source=pisa_source, download=download)
    downloader.download_missing_interfaces(pdb_list)
    return downloader.interface_dir


def pisa_pockets(records, pocket_dir, pisa_source=DEFAULT_PISA_SOURCE, bsa_threshold=None, download=True):
    """
    Build a Pocket per record from the PDBe PISA interface it names, with coordinates from its structure.

    Downloads any PISA files not already cached under `pocket_dir/pisa/`, unless `download` is False.
    With the "api" source that is one paced request per assembly, so an uncached list can take a long
    time. When any entry fails, `pisa/errors.json` is overwritten with the failures.

    Args:
        records (list): QTRecord dicts with `pocket_method == "pisa"`.
        pocket_dir (str): Pocket cache directory.
        pisa_source (str, optional): Where assembly interfaces are fetched from: "ftp" or "api".
            Defaults to DEFAULT_PISA_SOURCE.
        bsa_threshold (float or None, optional): Percent of a residue's accessible surface the partner
            must bury, strictly exceeded, for the residue to join. Not validated here. Defaults to None,
            which builds the pocket from interface bonds instead.
        download (bool, optional): Whether PISA entries not cached may be downloaded. Defaults to True.

    Returns:
        dict: pocket_id -> Pocket, for built pockets only. A record whose interface cannot be resolved
            (see `pisa_pockets_from_records`; with `download` False, also one whose entry is not cached),
            or whose structure or chain cannot be read, is left out with a warning. One with no residue
            under the rule keeps `pocket_exists` False.
    """
    log_extra = {"stage": "Retrieving PISA Pockets"}

    pdb_list = list(dict.fromkeys(record["struct_info"] for record in records))
    logger.debug(f"PDBs for which to retrieve PISA pockets: {pdb_list}", extra=log_extra)
    interface_dir = download_pisa_interfaces(pdb_list, pocket_dir, pisa_source, download=download)

    parsed = pisa_pockets_from_records(records, interface_dir, bsa_threshold=bsa_threshold)
    logger.debug(f"{len(parsed)} PISA pockets parsed from interfaces", extra=log_extra)

    # The interfaces give residue ids only; seq_pos and CA coordinates are added to the same Pocket
    pockets = {}
    parsed_records = [record for record in records if record["pocket_id"] in parsed]
    for st, group in records_by_structure(parsed_records, log_extra):
        for record in group:
            domain_chain, _ = split_chain_info(record["chain_info"])
            parsed_pocket = parsed[record["pocket_id"]]
            residue_ids = [int(x) for x in parsed_pocket.res_auth_ids]
            pocket = None if st is None else parse_pocket_from_struct(st, domain_chain, residue_ids, parsed_pocket)
            if pocket is None:
                logger.warning(
                    f"Could not read chain {domain_chain} of {record['struct_info']} for {record['pocket_id']}, "
                    "skipping this entry",
                    extra=log_extra,
                )
                continue
            pockets[record["pocket_id"]] = pocket
    return pockets
