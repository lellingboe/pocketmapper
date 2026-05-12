import os
import logging
from concurrent.futures import ThreadPoolExecutor
from urllib.request import urlcleanup, urlretrieve
import shutil
import gzip


class StructureFetcher:
    def __init__(self, out_dir):
        """
        Initialize the StructureFetcher with a target output directory.

        Args:
            out_dir (str): Directory where the downloaded structured files (e.g., mmCIF) will be saved and cached.
        """
        self.out_dir = out_dir
        self.cache = os.listdir(out_dir)
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
        """
        self.cache = os.listdir(self.out_dir)

    def get_structures(self, records):
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

    def fetch_structure(self, query):
        """
        Fetch a single structure based on its type and identifier.

        Args:
            query (dict): A dictionary describing the structure to fetch containing
                "struct_type" and "struct_info".

        Returns:
            tuple: A pair containing the structure identifier (str) and a boolean
                   indicating success (True) or failure (False).
        """
        stage = {"stage": "Fetching Structure"}
        match query["struct_type"]:
            case "alphafold":
                return self.fetch_alphafold(query["struct_info"])
            case "pdb":
                return self.fetch_mmcif(query["struct_info"])
            case "local_file":
                return (query["struct_info"], True)
            case "foldseek_db":
                # We assume that the foldseek db is already downloaded and available at the specified path, so we just check if the file exists
                return (query["struct_info"], True)
            case _:
                self.logger.warning(
                    f"Unknown structure type {query['struct_type']} for struct_info {query['struct_info']}", extra=stage
                )
                return (query["struct_info"], False)

    def fetch_alphafold(self, uniprot_acc, version="v6"):
        """
        Download an AlphaFold model in mmCIF format and compress it to gzip.

        Args:
            uniprot_acc (str): The UniProt accession number for the target structure.
            version (str, optional): The AlphaFold database version. Defaults to "v6".

        Returns:
            tuple: A pair containing the UniProt accession number (str) and a boolean
                   indicating success (True) or failure (False).
        """
        stage = {"stage": "Downloading AlphaFold File"}
        out_fname = f"{uniprot_acc}.cif.gz"
        temp_fname = f"{uniprot_acc}.cif"
        if not (out_fname in self.cache):
            url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_acc}-F1-model_{version}.cif"
            try:
                urlcleanup()
                urlretrieve(url, os.path.join(self.out_dir, temp_fname))
            except OSError:
                self.logger.warning(f"OSError when downloading {uniprot_acc}", extra=stage)
                return (uniprot_acc, False)
            except Exception:
                self.logger.warning(f"Atypical error when downloading {uniprot_acc}", extra=stage)
                return (uniprot_acc, False)

            # Compressing the downloaded cif file to gz format and removing the original cif file to save space
            with open(os.path.join(self.out_dir, temp_fname), "rb") as f_in:
                with gzip.open(os.path.join(self.out_dir, out_fname), "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(os.path.join(self.out_dir, temp_fname))

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
        out_fname = f"{pdb_code}.cif.gz"
        pdb_code_lowered = pdb_code.lower()
        if not (out_fname in self.cache):
            url = f"https://files.wwpdb.org/pub/pdb/data/structures/divided/mmCIF/{pdb_code_lowered[1:3]}/{pdb_code_lowered}.cif.gz"
            self.logger.debug(f"Attempting to download {pdb_code} from {url}", extra=stage)
            try:
                urlcleanup()
                urlretrieve(url, os.path.join(self.out_dir, out_fname))
            except OSError:
                self.logger.warning(f"Unable to download {pdb_code}", extra=stage)
                return (pdb_code, False)
            except Exception:
                self.logger.warning(f"Atypical error when downloading {pdb_code}", extra=stage)
                return (pdb_code, False)
        return (pdb_code, True)
