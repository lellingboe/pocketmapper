"""
Reduction of fetched structures to the single chain the alignment needs.

Foldseek indexes one structure per chain, so each record's reference structure is split down to
its alignment chain before the search directory is built. Parsing and writing are gemmi
throughout; the single-chain copies are written as gzipped mmCIF.
"""

import logging
import os
import shutil

import gemmi
from tqdm import tqdm

from pocketmapper.lib import gzip_file
from pocketmapper.lib import split_chain_info

logger = logging.getLogger(__name__)


class StructurePreprocessor:
    """
    Splits reference structures into single-chain copies for Foldseek to index.

    Each record's cached copy is named by its `preprocess_name`, so there is no call order to observe.
    """

    def __init__(
        self,
    ):
        """
        Initialise the logging stage. The preprocessor holds no other state.
        """
        self.log_extra = {"stage": "Preprocessing Structures"}
        logger.debug("Initialized")

    def preprocess_records(self, records, cache_dir, search_dir):
        """
        Split each record's reference structure down to its single alignment chain.

        Writes each single-chain copy to `<cache_dir>/<preprocess_name>.cif.gz`, creating the directory,
        unless a file is already there, and copies it into `search_dir` under the same name.

        Args:
            records (list): QTRecord dicts carrying `struct_path` and `preprocess_name`.
            cache_dir (str): Directory the single-chain copies are cached in.
            search_dir (str): Directory Foldseek will read the single-chain structures from.

        Returns:
            dict: pocket_id -> whether preprocessing succeeded. Foldseek-database records count as
                succeeded untouched.
        """
        status_dict = {}

        for record in tqdm(records):
            if record["struct_type"] == "foldseek_db":
                status_dict[record["pocket_id"]] = True  # foldseek db records are already preprocessed
                continue

            struct_info = record["struct_info"]
            chain_info = record["chain_info"]  # e.g., A_B or A
            chain, _ = split_chain_info(chain_info)
            # Ensuring divided structure is in the cache directory, e.g. <cache>/P12345_A_<md5>.cif.gz
            out_path_gz = os.path.join(cache_dir, f"{record['preprocess_name']}.cif.gz")
            out_path = out_path_gz.removesuffix(".gz")  # scratch copy, deleted once gzipped

            if not os.path.exists(out_path_gz):
                ref_path = record["struct_path"]  # e.g., /path/to/alphafold_dir/P12345.cif.gz
                st = gemmi.read_structure(ref_path)

                # Taking first model and deleting the rest
                del st[1:]
                model = st[0]

                # verify structure contains all interaction chains
                model_chains = set([chain.name for chain in model])
                if chain not in model_chains:
                    msg = f"Preprocessing: {struct_info} does not contain chain '{chain}' specified in chain_info '{chain_info}'"
                    logger.warning(
                        msg,
                        extra=self.log_extra,
                    )
                    status_dict[record["pocket_id"]] = False
                    # Falling through would write, and cache, an empty structure marked as a success
                    continue

                # Detaching all non interaction chains
                for chain_id in model_chains:
                    if chain_id != chain:
                        del model[chain_id]

                # Output the domain and motif pdb file
                os.makedirs(cache_dir, exist_ok=True)
                groups = gemmi.MmcifOutputGroups(False, atoms=True, group_pdb=True)
                st.make_mmcif_document(groups).write_file(out_path)
                # Through a .part: the cache trusts any .cif.gz it finds, so a truncated one would stick
                part_path_gz = f"{out_path_gz}.part"
                gzip_file(out_path, part_path_gz)
                os.replace(part_path_gz, out_path_gz)
                os.remove(out_path)

            search_path = os.path.join(search_dir, os.path.basename(out_path_gz))
            shutil.copyfile(out_path_gz, search_path)

            status_dict[record["pocket_id"]] = True

        return status_dict
