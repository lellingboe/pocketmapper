"""
Pockets through one interface: building them, the PISA cache and Foldseek hit expansion.

`POCKET_BUILDERS` is the whole table of pocket methods. Every builder takes `(records, structures)`, plus
the keyword options in `builder_options` for a method that has any, where `structures` is the manager's
`StructureManager`, and returns pocket_id -> Pocket for the pockets it built, leaving out (and logging) the
records it could not build. `PocketManager` sets each method's options once, so `build` needs to know nothing
about any one method. Built pockets are not cached: every `build` builds them afresh.
"""

import logging
import os
import re
from dataclasses import asdict

from pocketmapper.lib.constants import DEFAULT_PISA_SOURCE
from pocketmapper.lib.download import DOWNLOAD_WORKERS
from pocketmapper.lib.entries import failure
from pocketmapper.lib.exceptions import PocketMapperError
from pocketmapper.lib.foldseek import parse_foldseek_pdb_entry_name
from pocketmapper.lib.structure_manager import StructureManager
from pocketmapper.parse.qt_processor import QTProcessor
from pocketmapper.pockets.pisa import interface_partners, load_interfaces, pisa_pockets
from pocketmapper.pockets.pisa_downloader import PisaDownloader
from pocketmapper.pockets.structure import passthrough_pockets, whole_chain_pockets
from pocketmapper.pockets.vdw import vdw_pockets

logger = logging.getLogger(__name__)

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
