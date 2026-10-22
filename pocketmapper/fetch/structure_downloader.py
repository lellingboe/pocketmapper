"""
Retrieval of mmCIF structures from the wwPDB and AlphaFold.

`STRUCTURE_SOURCES` maps each downloadable `struct_type` to its URL builder and the transform its response
needs; a new source is one row. Structures land as gzipped mmCIF at each record's `struct_path`, which
doubles as the on-disk cache between runs.
"""

import logging

from pocketmapper.lib.download import DOWNLOAD_WORKERS, DownloadJob, download_many
from pocketmapper.lib.files import gzip_file

logger = logging.getLogger(__name__)


def pdb_url(pdb_code):
    """
    Build the wwPDB address of a PDB entry's divided mmCIF, which is served already gzipped.

    Args:
        pdb_code (str): The 4-character PDB code, in any case.

    Returns:
        str: The download URL.
    """
    code = pdb_code.lower()
    return f"https://files.wwpdb.org/pub/pdb/data/structures/divided/mmCIF/{code[1:3]}/{code}.cif.gz"


def alphafold_url(uniprot_acc):
    """
    Build the AlphaFold DB address of a UniProt accession's v6 model, which is served as plain mmCIF.

    Args:
        uniprot_acc (str): The UniProt accession.

    Returns:
        str: The download URL.
    """
    return f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_acc}-F1-model_v6.cif"


# struct_type -> (URL builder, transform). AlphaFold is gzipped on the way in to match the cache.
STRUCTURE_SOURCES = {"pdb": (pdb_url, None), "alphafold": (alphafold_url, gzip_file)}
# Types already on disk, or fetched elsewhere
NOT_DOWNLOADED = ("local_file", "foldseek_db")


def download_structures(records, max_workers=DOWNLOAD_WORKERS):
    """
    Concurrently fetch the structures a list of records names.

    Writes each structure to its record's `struct_path`, unless a file is already there. Records whose
    type is in `NOT_DOWNLOADED` report True without fetching; any other type warns and reports
    False.

    Args:
        records (list of dict): Records, each with `struct_type`, `struct_info` and, for the downloading
            types, `struct_path`.
        max_workers (int, optional): Downloads to run concurrently. Defaults to DOWNLOAD_WORKERS.

    Returns:
        dict: struct_info -> whether its structure is available.
    """
    log_extra = {"stage": "Downloading Structures"}
    found = {}
    jobs = {}
    for record in records:
        struct_type, struct_info = record["struct_type"], record["struct_info"]
        if struct_type in STRUCTURE_SOURCES:
            url_for, transform = STRUCTURE_SOURCES[struct_type]
            jobs[struct_info] = DownloadJob(url_for(struct_info), record["struct_path"], transform)
        elif struct_type in NOT_DOWNLOADED:
            found[struct_info] = True
        else:
            logger.warning(f"Unknown structure type {struct_type} for struct_info {struct_info}", extra=log_extra)
            found[struct_info] = False
    found.update(download_many(jobs, max_workers=max_workers, log_extra=log_extra))
    return found
