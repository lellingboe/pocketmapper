import os
import logging
from concurrent.futures import ThreadPoolExecutor
from urllib.request import urlcleanup, urlretrieve


class StructureFetcher:
    def __init__(self, out_dir):
        """
        :param out_dir: Directory where the mmCIF files will be saved.
        """
        self.out_dir = out_dir
        self.cache = os.listdir(out_dir)
        self.logger = logging.getLogger(__name__)
        self._log_extra = {"stage": "StructureFetcher"}
        self.logger.debug(f"Initialized with cache: {self.cache}", extra=self._log_extra)

    def _set_output_directory(self, out_dir):
        self.out_dir = out_dir
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

    def update_cache(self):
        self.cache = os.listdir(self.out_dir)

    def get_structures(self, records):
        """
        queries: list of records with keys "struct_type" and "struct_info"
        [
            {"struct_type": "alphafold", "struct_info": "P12345"},
            {"struct_type": "pdb", "struct_info": "1ABC"},
            ...
        ]
        returns:
        {
            "P12345": True,  # True if successfully fetched, False otherwise
            "1ABC": False,
            ...
        }
        """
        with ThreadPoolExecutor(max_workers=100) as e:
            results = e.map(self.get_structure, records)
        collected_result = {query: result for query, result in results}
        return collected_result

    def get_structure(self, query):
        stage = {"stage": "Fetching Structure"}
        match query["struct_type"]:
            case "alphafold":
                return self.get_alphafold(query["struct_info"])
            case "pdb":
                return self.get_mmcif(query["struct_info"])
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

    def get_alphafold(self, uniprot_acc):
        stage = {"stage": "Downloading AlphaFold File"}
        out_fname = f"{uniprot_acc}.cif"
        if not (out_fname in self.cache):
            url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_acc}-F1-model_v6.cif"
            try:
                urlcleanup()
                urlretrieve(url, os.path.join(self.out_dir, out_fname))
            except OSError:
                self.logger.warning(f"OSError when downloading {uniprot_acc}", extra=stage)
                return (uniprot_acc, False)
            except Exception:
                self.logger.warning(f"Atypical error when downloading {uniprot_acc}", extra=stage)
                return (uniprot_acc, False)
        return (uniprot_acc, True)

    def get_mmcif(self, pdb_code):
        stage = {"stage": "Downloading PDB File"}
        out_fname = f"{pdb_code}.cif.gz"
        pdb_code_lowered = pdb_code.lower()
        if not (out_fname in self.cache):
            url = f"https://files.wwpdb.org/pub/pdb/data/structures/divided/mmCIF/{pdb_code_lowered[1:3]}/{pdb_code_lowered}.cif.gz"
            try:
                urlcleanup()
                urlretrieve(url, os.path.join(self.out_dir, out_fname))
            except OSError:
                self.logger.warning(f"OSError when downloading {pdb_code}", extra=stage)
                return (pdb_code, False)
            except Exception:
                self.logger.warning(f"Atypical error when downloading {pdb_code}", extra=stage)
                return (pdb_code, False)
        return (pdb_code, True)
