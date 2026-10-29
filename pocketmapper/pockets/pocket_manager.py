"""
Pockets through one interface: building them, the PISA cache, Foldseek hit expansion and the pocket files.

`POCKET_BUILDERS` is the whole table of pocket methods. Every builder takes `(records, structures, **options)`,
where `structures` is the manager's `StructureManager`, and returns pocket_id -> Pocket for the pockets it
built, leaving out (and logging) the records it could not build. `PocketManager` sets each method's options
once, so `build` needs to know nothing about any one method. Built pockets are not cached: every `build`
builds them afresh.

The pocket files are module-level functions and need no manager: `write_pockets_file`/`read_pockets_file`
for the versioned pockets file, which also maps each aligned chain to the pockets on it, and
`write_pockets_table` for a table of every pocket.
"""

import json
import logging
import os
import re
from dataclasses import asdict

import pandas as pd

from pocketmapper.lib.constants import DEFAULT_PISA_SOURCE
from pocketmapper.lib.download import DOWNLOAD_WORKERS
from pocketmapper.lib.entries import failure
from pocketmapper.lib.exceptions import PocketMapperError
from pocketmapper.lib.foldseek import parse_foldseek_pdb_entry_name
from pocketmapper.lib.residues import split_chain_info
from pocketmapper.lib.structure_manager import StructureManager
from pocketmapper.parse.qt_processor import QTProcessor
from pocketmapper.pockets.pisa import interface_partners, load_interfaces, pisa_pockets
from pocketmapper.pockets.pisa_downloader import PisaDownloader
from pocketmapper.pockets.pocket import Pocket, PocketResidue
from pocketmapper.pockets.structure import passthrough_pockets, whole_chain_pockets
from pocketmapper.pockets.vdw import vdw_pockets

logger = logging.getLogger(__name__)

# The pockets file format `write_pockets_file` writes and `read_pockets_file` accepts
POCKETS_FILE_VERSION = 2

# pocket_method -> (log label, dedup keys, builder), in merge order.
# Only pisa dedups: its pocket is fully determined by structure and chain. The others must not --
# two passthrough pockets on one chain differ only in `residue_info`.
POCKET_BUILDERS = {
    "pisa": ("PISA", ("struct_info", "chain_info"), pisa_pockets),
    "passthrough": ("Passthrough", None, passthrough_pockets),
    "vdw": ("VdW", None, vdw_pockets),
    "whole_chain": ("Whole Chain", None, whole_chain_pockets),
}


class PocketManager:
    """
    Builds pockets and expands Foldseek hits into pockets, under one download policy.

    Holds one `PisaDownloader`, caching under `pocket_dir/pisa/` and shared by the pisa builder and
    `expand_pdb_hits`, and one `StructureManager` for the hits' structures. `download` governs both, and
    `max_workers` sets the width of both download pools. Builders locate structures with that
    manager's `record_path` and read them with `lib.structure_manager.load_structure`, which never downloads.
    """

    def __init__(
        self,
        pocket_dir,
        pdb_dir=None,
        alphafold_dir=None,
        fsdb_dir=None,
        local_root=None,
        pisa_source=DEFAULT_PISA_SOURCE,
        bsa_threshold=None,
        download=True,
        max_workers=DOWNLOAD_WORKERS,
    ):
        """
        Set up the PISA cache, the structure fetcher and each pocket method's options.

        Writes nothing.

        Args:
            pocket_dir (str): Parent of the PISA cache, `pocket_dir/pisa/`.
            pdb_dir (str or None, optional): Cache directory of PDB structures. Needed to locate a `pdb`
                record's structure, so by `build` for such records and by `expand_pdb_hits`. Defaults to None.
            alphafold_dir (str or None, optional): Cache directory of AlphaFold models. Needed to locate an
                `alphafold` record's structure. Defaults to None.
            fsdb_dir (str or None, optional): Directory of downloaded Foldseek databases. Needed to locate a
                bundled Foldseek database. Defaults to None.
            local_root (str or None, optional): Directory a relative local structure file resolves against.
                Defaults to None, the current directory.
            pisa_source (str, optional): Where PISA assembly interfaces are fetched from: "ftp" or "api".
                Defaults to DEFAULT_PISA_SOURCE.
            bsa_threshold (float or None, optional): Percent of a residue's accessible surface the partner
                chain must bury, strictly exceeded, for the residue to join a pisa pocket. Not validated
                here. Defaults to None, which builds pisa pockets from interface bonds.
            download (bool, optional): False sends no request: only cached PISA data and structures are
                used. Defaults to True.
            max_workers (int, optional): Width of the PISA and structure download pools. Defaults to
                DOWNLOAD_WORKERS.

        Raises:
            PocketMapperError: If `pisa_source` is not one of PISA_SOURCES.
        """
        self.pisa = PisaDownloader(
            os.path.join(pocket_dir, "pisa"), source=pisa_source, max_workers=max_workers, download=download
        )
        self.structures = StructureManager(
            pdb_dir, alphafold_dir, fsdb_dir, local_root, download=download, max_workers=max_workers
        )
        # pocket_method -> keyword arguments of its builder; a method absent here takes none
        self.builder_options = {"pisa": {"downloader": self.pisa, "bsa_threshold": bsa_threshold}}

    def build(self, records):
        """
        Build a Pocket for each record, by the builder of its `pocket_method`.

        Downloads any PISA data not already cached, unless downloads are off. PISA records sharing a
        structure and chain are built once, under the first one's pocket_id. Methods are merged in
        `POCKET_BUILDERS` order, so a later method's pocket replaces an earlier one's of the same pocket_id.

        Args:
            records (list): QTRecord dicts; reads `pocket_id`, `pocket_method`, `struct_info`,
                `struct_type`, `chain_info` and `residue_info`. Query and target records may be mixed.
                Used as given: filter out any whose structure is unavailable before calling.

        Returns:
            dict: pocket_id -> Pocket, for built pockets only. A record whose pocket cannot be built is left
                out and logged by its builder; one built but empty is kept, with `pocket_exists=False`.
        """
        log_extra = {"stage": "Getting Pockets"}
        logger.info("Starting pocket retrieval...", extra=log_extra)

        pockets = {}
        for pocket_method, (label, dedup_keys, builder) in POCKET_BUILDERS.items():
            method_log_extra = {"stage": f"Retrieving {label} Pockets"}
            logger.info(f"Checking for {label} pockets...", extra=method_log_extra)

            method_records = [record for record in records if record["pocket_method"] == pocket_method]
            if dedup_keys is not None:
                unique = {}
                for record in method_records:
                    unique.setdefault(tuple(record[k] for k in dedup_keys), record)
                method_records = list(unique.values())

            if not method_records:
                logger.info(f"No {label} pockets to retrieve", extra=method_log_extra)
                continue
            logger.info(f"{len(method_records)} {label} pockets to retrieve", extra=method_log_extra)

            method_pockets = builder(method_records, self.structures, **self.builder_options.get(pocket_method, {}))
            logger.debug(f"Built {len(method_pockets)} {label} pockets", extra=method_log_extra)
            pockets |= method_pockets

        logger.debug(f"Built {len(pockets)} pockets in all", extra=log_extra)
        return pockets

    def expand_pdb_hits(self, hit_names):
        """
        Turn the hits of a PDB Foldseek-database search into pisa records, one per PISA interface.

        Downloads the PISA data and the structures those records need, unless downloads are off. Each
        record's `preprocess_name` is the hit name rather than the one `QTProcessor` derives.

        Args:
            hit_names (list): The database entry names the search hit.

        Returns:
            tuple: (records, chains, failures). `records` holds those whose structure is available.
                `chains` maps every hit name to the pocket_ids of its records, including those whose
                structure could not be fetched; it is empty for a hit that is not PDB-style or has no
                interface the pisa pocket method can spell. `failures` has one `structure_not_found`
                failure per such record, as `lib.entries.failure` builds it, with no side.

        Raises:
            PocketMapperError: If `pdb_dir` was not given, or if interfaces were found but none of their
                structures could be fetched.
        """
        log_extra = {"stage": "Expanding Foldseek DB Targets"}
        qtprocessor = QTProcessor()

        chains = {hit_name: [] for hit_name in hit_names}
        hits = {}  # foldseek entry name -> (pdb_id, chain_id)
        for hit_name in hit_names:
            resolved = parse_foldseek_pdb_entry_name(hit_name)
            if resolved is not None:
                hits[hit_name] = resolved

        pdb_list = sorted({pdb_id for pdb_id, _ in hits.values()})
        logger.info(
            f"Retrieving PISA interfaces for {len(pdb_list)} PDB entries behind {len(hits)} Foldseek hits",
            extra=log_extra,
        )
        self.pisa.download_missing_interfaces(pdb_list)

        # Each entry's interfaces are read once, however many hits share it
        hit_chains = {}  # pdb_id -> the chain ids hit on it, as dict keys in hit order
        for pdb_id, chain_id in hits.values():
            hit_chains.setdefault(pdb_id, {})[chain_id] = None
        partners = {}  # (pdb_id, chain_id) -> partner chain ids
        for pdb_id, chain_ids in hit_chains.items():
            interfaces = load_interfaces(pdb_id, self.pisa.interface_dir)
            if interfaces is None:
                logger.debug(f"Could not load PISA data for {pdb_id}", extra=log_extra)
                continue
            for chain_id in chain_ids:
                partners[(pdb_id, chain_id)] = interface_partners(interfaces, chain_id)
                if not partners[(pdb_id, chain_id)]:
                    logger.debug(f"No PISA interface involving chain {chain_id} of {pdb_id}", extra=log_extra)

        # Building one record per interface the hit chain takes part in.
        # PISA stores chain pairs the input grammar cannot spell (multi-character chain ids); those are
        # counted and skipped here rather than each rejected with a warning.
        pisa_pattern = qtprocessor.pocket_methods["pisa"][0]
        unspellable = 0
        records = []
        for hit_name, (pdb_id, chain_id) in hits.items():
            for partner in partners.get((pdb_id, chain_id), []):
                if not re.match(pisa_pattern, f"{chain_id}_{partner}"):
                    unspellable += 1
                    continue
                record, reason = qtprocessor.parse_individual_qt(f"{pdb_id}:{chain_id}_{partner}", pocket_method="pisa")
                if record is None:
                    logger.warning(f"{reason}; skipping this interface", extra=log_extra)
                    continue
                # The alignment is keyed by the Foldseek entry name
                record.preprocess_name = hit_name
                records.append(asdict(record))
                chains[hit_name].append(record.pocket_id)
        if unspellable:
            logger.info(
                f"Skipped {unspellable} PISA interfaces whose chain ids the pisa pocket method cannot express",
                extra=log_extra,
            )
        if not records:
            logger.warning("No PISA interfaces found for any Foldseek hit", extra=log_extra)
            return [], chains, []

        # Fetched last, so only entries with an interface are downloaded. Failed fetches are dropped.
        found = self.structures.fetch(records)
        logger.info(
            f"{sum(found.values())}/{len(found)} foldseek hit required structures available",
            extra={"stage": "Downloading Structures"},
        )
        if not any(found.values()):
            logger.critical("Insufficient foldseek hit structures after fetching", extra=log_extra)
            raise PocketMapperError(
                "Insufficient foldseek hit structures after fetching. No valid foldseek hit entries remain."
            )
        failures = [
            failure(None, record, "structure_not_found") for record in records if not found[record["struct_info"]]
        ]
        records = [record for record in records if found[record["struct_info"]]]

        logger.info(
            f"Added {len(records)} PISA target pockets from {len({r['preprocess_name'] for r in records})} Foldseek hits",
            extra=log_extra,
        )
        return records, chains, failures


def write_pockets_file(pockets, chains, path):
    """
    Write a pockets file: the pockets plus the map from aligned chain to the pockets on it.

    Creates the file's directory if needed.

    Args:
        pockets (dict): pocket_id -> Pocket, or None for an entry given and not built.
        chains (dict): preprocess_name -> the pocket_ids on that chain; empty for a chain with none.
        path (str): File to write.

    Returns:
        None
    """
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as f:
        json.dump({"version": POCKETS_FILE_VERSION, "pockets": serialise_pockets(pockets), "chains": chains}, f)


def read_pockets_file(path):
    """
    Read a pockets file written by `write_pockets_file`.

    Only the shape and version are checked. A hand-written file with `res_auth_ids` out of order or a
    wrong `seq_pos` reads fine and gives silently wrong overlap.

    Args:
        path (str): File to read.

    Returns:
        tuple: (pockets, chains), as `write_pockets_file` takes them; pockets rebuilt as by
            `deserialise_pockets`.

    Raises:
        PocketMapperError: If the file is not JSON, is not a pockets file, or is of another version.
    """
    try:
        with open(path) as f:
            serialised = json.load(f)
    except ValueError as e:
        msg = f"Could not read the pockets file {path}: {e}"
        logger.critical(msg)
        raise PocketMapperError(msg) from e
    version = serialised.get("version") if isinstance(serialised, dict) else None
    if version != POCKETS_FILE_VERSION:
        msg = f"{path} is not a version {POCKETS_FILE_VERSION} pockets file (version {version}); rerun pockets"
        logger.critical(msg)
        raise PocketMapperError(msg)
    try:
        return deserialise_pockets(serialised["pockets"]), dict(serialised["chains"])
    except (TypeError, AttributeError, KeyError, ValueError) as e:
        msg = f"Could not read the pockets file {path}: {e}"
        logger.critical(msg)
        raise PocketMapperError(msg) from e


def serialise_pockets(pockets):
    """
    Turn a pocket collection into JSON-ready dicts.

    Args:
        pockets (dict): pocket_id -> Pocket, or None.

    Returns:
        dict: pocket_id -> the Pocket's fields, or None.
    """
    return {pid: asdict(pocket) if pocket is not None else None for pid, pocket in pockets.items()}


def deserialise_pockets(serialised):
    """
    Rebuild a pocket collection from `serialise_pockets` output.

    Args:
        serialised (dict): pocket_id -> a Pocket's fields, or None.

    Returns:
        dict: pocket_id -> Pocket, or None. Residues are rebuilt as PocketResidues, keeping
            `res_auth_ids` in order and every None field.
    """
    pockets = {}
    for pid, fields in serialised.items():
        if fields is None:
            pockets[pid] = None
            continue
        residues = {res_id: PocketResidue(**residue) for res_id, residue in fields["residues"].items()}
        pockets[pid] = Pocket(**(fields | {"residues": residues}))
    return pockets


def write_pockets_table(records, pockets, path):
    """
    Write a table of every pocket: pocket_id, chain, res_auth_ids and method.

    Creates the table's directory if needed.

    Args:
        records (list): The QTRecord dicts the pockets were built for; a pocket_id listed twice (on both
            sides, or as a given record and a Foldseek hit) is written once, from its first record.
        pockets (dict): pocket_id -> Pocket, or None for one not built, whose res_auth_ids are left
            empty.
        path (str): The table to write, tab-separated.

    Returns:
        None
    """
    log_extra = {"stage": "Getting Pockets"}
    rows = {}
    for record in records:
        pocket = pockets.get(record["pocket_id"])
        rows.setdefault(
            record["pocket_id"],
            {
                "pocket_id": record["pocket_id"],
                "chain": split_chain_info(record["chain_info"])[0],
                "res_auth_ids": ",".join(map(str, pocket.res_auth_ids)) if pocket is not None else "",
                "method": record["pocket_method"],
            },
        )
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    pd.DataFrame(list(rows.values()), columns=["pocket_id", "chain", "res_auth_ids", "method"]).to_csv(
        path, sep="\t", index=False
    )
    logger.info(f"Pockets listed in {path}", extra=log_extra)
