"""
Retrieval of mmCIF structures from the wwPDB and AlphaFold.

Downloads are concurrent -- the pool width is a constructor argument -- and land as gzipped mmCIF at
each record's `struct_path`, which doubles as the on-disk cache between runs.
"""

import logging
import os
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import as_completed

from tqdm import tqdm

from pocketmapper.downloads.lib_download import download_file
from pocketmapper.lib import gzip_file

logger = logging.getLogger(__name__)

# Default pool width. Fixed rather than tied to a thread count: a worker waits on a socket, not on a core.
DOWNLOAD_WORKERS = 8


class StructureDownloader:
    """
    Downloads PDB and AlphaFold structures to the paths their records name.

    Each record carries its own destination, so there is no shared output directory and no call
    order to observe. A destination whose parent directory does not exist is reported as a failed
    download rather than created; the caller owns the directory.

    How many downloads run at once is fixed at construction, in `max_workers`.
    """

    def __init__(self, max_workers=DOWNLOAD_WORKERS, max_retries=5, base_delay=0.25, max_delay=30.0):
        """
        Store the pool width and the retry budget shared by every download this instance makes.

        Args:
            max_workers (int): Downloads to run concurrently. Defaults to DOWNLOAD_WORKERS.
            max_retries (int): Attempts per download before giving up. Defaults to 5.
            base_delay (float): Seconds to wait after the first failed attempt. Defaults to 0.25.
            max_delay (float): Ceiling on the doubling backoff delay. Defaults to 30.0.
        """
        self.max_workers = max_workers
        self.max_retries = max_retries
        self.base_delay = base_delay
        self.max_delay = max_delay
        self.log_extra = {"stage": "Downloading Structures"}

    def download_missing_structures(self, records):
        """
        Concurrently fetch multiple structures based on query records.

        Args:
            records (list of dict): A list of dictionaries, where each dict has:
                - "struct_type" (str): Type of the structure (e.g., "alphafold", "pdb").
                - "struct_info" (str): Identifier for the structure (e.g., "P12345", "1ABC").
                - "struct_path" (str): Destination file path, read for the downloading types only.

        Returns:
            dict: A mapping of the structure identifier to a boolean status indicating
                  whether the fetch was successful (True) or not (False).
        """
        with ThreadPoolExecutor(max_workers=self.max_workers) as e:
            futures = [e.submit(self.download_missing_structure, record) for record in records]
            results = [f.result() for f in tqdm(as_completed(futures), total=len(futures))]
        return dict(results)

    def download_missing_structure(self, record):
        """
        Fetch a single structure, dispatching on its type.

        Local-file and Foldseek-database records are assumed already on disk and report success without
        downloading anything.

        Args:
            record (dict): Describes the structure to fetch; reads `struct_type`, `struct_info` and,
                for the downloading types, `struct_path`.

        Returns:
            tuple: (struct_info, succeeded).
        """
        match record["struct_type"]:
            case "alphafold":
                return self.download_alphafold(record["struct_info"], record["struct_path"])
            case "pdb":
                return self.download_mmcif(record["struct_info"], record["struct_path"])
            case "local_file":
                return (record["struct_info"], True)
            case "foldseek_db":
                # We assume that the foldseek db is already downloaded and available at the specified path, so we just check if the file exists
                return (record["struct_info"], True)
            case _:
                logger.warning(
                    f"Unknown structure type {record['struct_type']} for struct_info {record['struct_info']}",
                    extra=self.log_extra,
                )
                return (record["struct_info"], False)

    def download_alphafold(self, uniprot_acc, out_fpath, version="v6"):
        """
        Download an AlphaFold model in mmCIF format and compress it to gzip.

        Args:
            uniprot_acc (str): The UniProt accession number for the target structure.
            out_fpath (str): Destination path for the gzipped mmCIF. An existing file here is
                treated as cached and nothing is fetched.
            version (str, optional): The AlphaFold database version. Defaults to "v6".

        Returns:
            tuple: A pair containing the Uniprot accession (str) and a boolean
                   indicating success (True) or failure (False).
        """
        if not os.path.exists(out_fpath):
            url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_acc}-F1-model_{version}.cif"
            # AlphaFold serves plain mmCIF; the cache holds it gzipped, so it is compressed on the
            # way in rather than stored twice.
            if not download_file(
                url,
                out_fpath,
                transform=gzip_file,
                max_retries=self.max_retries,
                base_delay=self.base_delay,
                max_delay=self.max_delay,
                log_extra=self.log_extra,
            ):
                return (uniprot_acc, False)
        return (uniprot_acc, True)

    def download_mmcif(self, pdb_code, out_fpath):
        """
        Download a PDB structure, which the wwPDB already serves as gzipped mmCIF.

        Args:
            pdb_code (str): The 4-character PDB code for the target structure.
            out_fpath (str): Destination path for the gzipped mmCIF. An existing file here is
                treated as cached and nothing is fetched.

        Returns:
            tuple: A pair containing the PDB code (str) and a boolean indicating
                   success (True) or failure (False).
        """
        pdb_code_lowered = pdb_code.lower()
        if not os.path.exists(out_fpath):
            url = f"https://files.wwpdb.org/pub/pdb/data/structures/divided/mmCIF/{pdb_code_lowered[1:3]}/{pdb_code_lowered}.cif.gz"
            if not download_file(
                url,
                out_fpath,
                max_retries=self.max_retries,
                base_delay=self.base_delay,
                max_delay=self.max_delay,
                log_extra=self.log_extra,
            ):
                return (pdb_code, False)
        return (pdb_code, True)
