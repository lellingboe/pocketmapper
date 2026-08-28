"""
Retrieval of mmCIF structures from the wwPDB and AlphaFold.

Downloads are concurrent and land as gzipped mmCIF in the output directory, which doubles as the
on-disk cache between runs.
"""

import os
import logging
from concurrent.futures import ThreadPoolExecutor
from urllib.request import urlcleanup, urlretrieve
import shutil
import gzip


def _discard_partial(fpath):
    """
    Delete a partially written download, ignoring the case where it never got created.

    Args:
        fpath (str): Path to remove.

    Returns:
        None
    """
    if os.path.exists(fpath):
        os.remove(fpath)


class StructureFetcher:
    """
    Downloads PDB and AlphaFold structures into a cached output directory.

    Call in order -- `set_output_directory()`, then `update_cache()`, then `fetch_structures()`.
    The ordering is required and nothing enforces it.
    """

    def __init__(self):
        """
        Initialise with no output directory; `set_output_directory` supplies it later.
        """
        self.out_dir = None
        self.cache = None
        self.logger = logging.getLogger(__name__)
        self._log_extra = {"stage": "StructureFetcher"}
        self.logger.debug(f"Initialized with cache: {self.cache}", extra=self._log_extra)

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

        A `<name>.cif.gz.part` left by an interrupted download lands here too, but can never match a
        `.cif.gz` lookup, so it is inert.
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
        temp_fpath = os.path.join(self.out_dir, f"{uniprot_acc}.cif")
        part_fpath = f"{out_fpath}.part"
        if os.path.basename(out_fpath) not in self.cache:
            url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_acc}-F1-model_{version}.cif"
            try:
                urlcleanup()
                urlretrieve(url, temp_fpath)
            except OSError:
                _discard_partial(temp_fpath)
                self.logger.warning(f"OSError when downloading {uniprot_acc}", extra=stage)
                return (uniprot_acc, False)
            except Exception:
                _discard_partial(temp_fpath)
                self.logger.warning(f"Atypical error when downloading {uniprot_acc}", extra=stage)
                return (uniprot_acc, False)

            # Compressing the downloaded cif file to gz format and removing the original cif file to save
            # space. The gzip goes to a .part first and is moved into place only once complete: the cache
            # trusts any .cif.gz it finds, so a truncated one under the real name would be served for good.
            with open(temp_fpath, "rb") as f_in:
                with gzip.open(part_fpath, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.replace(part_fpath, out_fpath)
            os.remove(temp_fpath)

        return (uniprot_acc, True)

    def fetch_mmcif(self, pdb_code):
        """
        Download a PDB structure in mmCIF format and compress it to gzip.

        Args:
            pdb_code (str): The 4-character PDB code for the target structure.

        Returns:
            tuple: A pair containing the PDB code (str) and a boolean indicating
                   success (True) or failure (False).
        """
        stage = {"stage": "Downloading PDB File"}
        out_fpath = os.path.join(self.out_dir, f"{pdb_code}.cif.gz")
        part_fpath = f"{out_fpath}.part"
        pdb_code_lowered = pdb_code.lower()
        if os.path.basename(out_fpath) not in self.cache:
            url = f"https://files.wwpdb.org/pub/pdb/data/structures/divided/mmCIF/{pdb_code_lowered[1:3]}/{pdb_code_lowered}.cif.gz"
            self.logger.debug(f"Attempting to download {pdb_code} from {url}", extra=stage)
            try:
                urlcleanup()
                urlretrieve(url, part_fpath)
            except OSError:
                _discard_partial(part_fpath)
                self.logger.warning(f"Unable to download {pdb_code}", extra=stage)
                return (pdb_code, False)
            except Exception:
                _discard_partial(part_fpath)
                self.logger.warning(f"Atypical error when downloading {pdb_code}", extra=stage)
                return (pdb_code, False)
            # Moved into place only once the download is complete -- see fetch_alphafold.
            os.replace(part_fpath, out_fpath)
        return (pdb_code, True)
