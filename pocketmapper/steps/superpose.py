"""
Step 6: superpose the top targets of each query onto it, one PDB per query.
"""

import logging

from pocketmapper.records import fsdb_record
from pocketmapper.records import read_records
from pocketmapper.records import synthesise_target_pockets
from pocketmapper.structure_aligner import StructureAligner

logger = logging.getLogger(__name__)


def superpose_top_targets(
    query_records,
    target_records,
    pocket_comparison,
    alignment,
    aligned_structure_dir,
    align_struct_method,
    align_count,
    threads,
):
    """
    Superpose the top `align_count` targets of each query onto it, into `aligned_structure_dir`.

    Against a Foldseek database, the target structures are rebuilt out of it under
    `aligned_structure_dir`.

    Args:
        query_records (str): The query records file.
        target_records (str): The target records file.
        pocket_comparison (str): The pocket comparison table.
        alignment (str): The alignment table; its transforms are read with the "foldseek" method.
        aligned_structure_dir (str): Where the superposed structures are written.
        align_struct_method (str): "pocket" or "foldseek", already resolved.
        align_count (int): Most targets to superpose onto each query; 0 writes nothing.
        threads (int): Thread count for rebuilding structures out of a Foldseek database.

    Returns:
        None

    Raises:
        PocketMapperError: If a records file cannot be read, or Foldseek fails rebuilding structures.
    """
    query = read_records(query_records)
    target = read_records(target_records)

    # With a Foldseek database, a hit's transform fits the database's own structure, so that is
    # what is superposed. Given no target records, the aligner reads target ids as entry names: only
    # a PDB database's hits have records of their own.
    database = fsdb_record(target)
    fsdb_path = database["struct_path"] if database is not None else None
    if synthesise_target_pockets(target):
        target = []

    StructureAligner().align_structs(
        query_records=query,
        target_records=target,
        pocket_comparison=pocket_comparison,
        out_dir=aligned_structure_dir,
        method=align_struct_method,
        align_count=align_count,
        alignment=alignment,
        threads=threads,
        fsdb_path=fsdb_path,
    )
