"""
Retrieval of mmCIF structures from the wwPDB and AlphaFold.

Downloads are concurrent and land as gzipped mmCIF in the output directory, which doubles as the
on-disk cache between runs.
"""

import logging
import os
from concurrent.futures import ThreadPoolExecutor

from pocketmapper.downloads.lib_download import download_file
from pocketmapper.lib import gzip_file


class StructureDownloader:
    """
    Downloads PDB and AlphaFold structures into a cached output directory.

    Call in order -- `set_output_directory()`, then `update_cache()`, then `fetch_structures()`.
    The ordering is required and nothing enforces it.
    """

    def __init__(self, max_retries=5, base_delay=0.25, max_delay=30.0):
        """
        Initialise with no output directory; `set_output_directory` supplies it later.

        Args:
            max_retries (int): Attempts per download before giving up. Defaults to 5.
            base_delay (float): Seconds to wait after the first failed attempt. Defaults to 0.25.
            max_delay (float): Ceiling on the doubling backoff delay. Defaults to 30.0.
        """
        self.max_retries = max_retries
        self.base_delay = base_delay
        self.max_delay = max_delay
        self.out_dir = None
        self.cache = None
        self.logger = logging.getLogger(__name__)
        self.log_extra = {"stage": "StructureFetcher"}
        self.logger.debug(f"Initialized with cache: {self.cache}", extra=self.log_extra)

    def set_output_directory(self, out_dir):
        """
        Set or update the target output directory, creating it if it does not exist.

        Args:
            out_dir (str): The new output directory path.
        """
        self.out_dir = out_dir
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

    def update_cache(self):
        """
        Update the internal cache of files present in the output directory.

        Holds bare filenames, so the callers below must test `os.path.basename(...)` for membership --
        testing the joined output path can never match and silently re-fetches every structure on every
        run. Left read-only during `fetch_structures`, which reads it from 100 threads at once.

        A `<name>.cif.gz.part` or `<name>.cif.gz.raw.part` left by an interrupted download lands here
        too, but neither can match a `.cif.gz` lookup, so they are inert.
        """
        self.cache = set(os.listdir(self.out_dir))

    def fetch_structures(self, records):
        """
        Concurrently fetch multiple structures based on query records.

        Args:
            records (list of dict): A list of dictionaries, where each dict has:
                - "struct_type" (str): Type of the structure (e.g., "alphafold", "pdb").
                - "struct_info" (str): Identifier for the structure (e.g., "P12345", "1ABC").

        Returns:
            dict: A mapping of the structure identifier to a boolean status indicating
                  whether the fetch was successful (True) or not (False).
        """
        with ThreadPoolExecutor(max_workers=100) as e:
            results = e.map(self.fetch_structure, records)
        collected_result = {query: result for query, result in results}
        return collected_result

    def fetch_structure(self, record):
        """
        Fetch a single structure, dispatching on its type.

        Local-file and Foldseek-database records are assumed already on disk and report success without
        downloading anything.

        Args:
            record (dict): Describes the structure to fetch; reads `struct_type` and `struct_info`.

        Returns:
            tuple: (struct_info, succeeded).
        """
        stage = {"stage": "Fetching Structure"}
        match record["struct_type"]:
            case "alphafold":
                return self.fetch_alphafold(record["struct_info"])
            case "pdb":
                return self.fetch_mmcif(record["struct_info"])
            case "local_file":
                return (record["struct_info"], True)
            case "foldseek_db":
                # We assume that the foldseek db is already downloaded and available at the specified path, so we just check if the file exists
                return (record["struct_info"], True)
            case _:
                self.logger.warning(
                    f"Unknown structure type {record['struct_type']} for struct_info {record['struct_info']}",
                    extra=stage,
                )
                return (record["struct_info"], False)

    def fetch_alphafold(self, uniprot_acc, version="v6"):
        """
        Download an AlphaFold model in mmCIF format and compress it to gzip.

        Args:
            uniprot_acc (str): The UniProt accession number for the target structure.
            version (str, optional): The AlphaFold database version. Defaults to "v6".

        Returns:
            tuple: A pair containing the Uniprot accession (str) and a boolean
                   indicating success (True) or failure (False).
        """
        stage = {"stage": "Downloading AlphaFold File"}
        out_fpath = os.path.join(self.out_dir, f"{uniprot_acc}.cif.gz")
        if os.path.basename(out_fpath) not in self.cache:
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
                log_extra=stage,
            ):
                return (uniprot_acc, False)
        return (uniprot_acc, True)

    def fetch_mmcif(self, pdb_code):
        """
        Download a PDB structure, which the wwPDB already serves as gzipped mmCIF.

        Args:
            pdb_code (str): The 4-character PDB code for the target structure.

        Returns:
            tuple: A pair containing the PDB code (str) and a boolean indicating
                   success (True) or failure (False).
        """
        stage = {"stage": "Downloading PDB File"}
        out_fpath = os.path.join(self.out_dir, f"{pdb_code}.cif.gz")
        pdb_code_lowered = pdb_code.lower()
        if os.path.basename(out_fpath) not in self.cache:
            url = f"https://files.wwpdb.org/pub/pdb/data/structures/divided/mmCIF/{pdb_code_lowered[1:3]}/{pdb_code_lowered}.cif.gz"
            if not download_file(
                url,
                out_fpath,
                max_retries=self.max_retries,
                base_delay=self.base_delay,
                max_delay=self.max_delay,
                log_extra=stage,
            ):
                return (pdb_code, False)
        return (pdb_code, True)
