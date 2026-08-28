"""
Reduction of fetched structures to the single chain the alignment needs.

Foldseek indexes one structure per chain, so each record's reference structure is split down to
its alignment chain before the search directory is built. Parsing and writing are gemmi
throughout; the single-chain copies are written as gzipped mmCIF.
"""

import gzip
import logging
import os
import shutil
import gemmi
from tqdm import tqdm


class StructurePreprocessor:
    """
    Splits reference structures into single-chain copies for Foldseek to index.

    Call in order -- `set_output_directory()`, then `update_cache()`, then `preprocess_records()`.
    The ordering is required: the cache is read from the output directory, and preprocessing needs
    both. Nothing enforces it.
    """

    def __init__(
        self,
    ):
        """
        Initialise with no output directory; `set_output_directory` supplies it later.
        """
        self.logger = logging.getLogger(__name__)
        self._log_extra = {"stage": "Structure Preprocessor"}
        self.logger.debug("Initialized", extra=self._log_extra)

        self.out_dir = None
        self.cache = None

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

        Note this cache never hits: it holds `os.listdir` basenames while the callers below test a full
        output path for membership, so every structure is re-fetched on every run. Fix the comparison
        rather than working around it if this becomes the bottleneck.
        """
        self.cache = os.listdir(self.out_dir)

    def preprocess_records(self, records, search_dir):
        """
        Split each record's reference structure down to its single alignment chain.

        The single-chain copy is cached under the output directory set by `set_output_directory()` and
        then copied into `search_dir` for Foldseek to index. Records already in a Foldseek database are
        passed through untouched.

        Args:
            records (list): QTRecord dicts carrying `struct_path` and the `preprocess_*` paths.
            search_dir (str): Directory Foldseek will read the single-chain structures from.

        Returns:
            dict: pocket_id -> whether preprocessing succeeded.
        """
        status_dict = {}
        stage = {"stage": "Dividing structures"}

        for record in tqdm(records):
            if record["struct_type"] == "foldseek_db":
                status_dict[record["pocket_id"]] = True  # foldseek db records are already preprocessed
                continue
            elif record["success"] is False:
                continue

            struct_info = record["struct_info"]
            chain_info = record["chain_info"]  # e.g., A_B or A
            chain = chain_info[0]  # e.g., "A"
            # Ensuring divided structure is in the cache directory
            out_path = record[
                "preprocess_path"
            ]  # e.g., /path/to/foldseek_preprocessed_structure_dir/P12345_A_<md5>.cif
            out_path_gz = record[
                "preprocess_path_gz"
            ]  # e.g., /path/to/foldseek_preprocessed_structure_dir/P12345_A_<md5>.cif.gz

            if not (out_path_gz in self.cache):
                ref_path = record["struct_path"]  # e.g., /path/to/structure_dir/P12345.cif.gz
                st = gemmi.read_structure(ref_path, format=gemmi.CoorFormat.Mmcif)

                # Taking first model and deleting the rest
                del st[1:]
                model = st[0]

                # verify structure contains all interaction chains
                model_chains = set([chain.name for chain in model])
                if chain not in model_chains:
                    msg = f"Preprocessing: {struct_info} does not contain chain '{chain}' specified in chain_info '{chain_info}'"
                    logging.warning(
                        msg,
                        extra=stage,
                    )
                    status_dict[record["pocket_id"]] = False

                # Detaching all non interaction chains
                for chain_id in model_chains:
                    if chain_id != chain:
                        del model[chain_id]

                # Output the domain and motif pdb file
                groups = gemmi.MmcifOutputGroups(False, atoms=True, group_pdb=True)
                st.make_mmcif_document(groups).write_file(out_path)
                with open(out_path, "rb") as f_in:
                    with gzip.open(out_path_gz, "wb") as f_out:
                        shutil.copyfileobj(f_in, f_out)
                os.remove(out_path)

            search_path = os.path.join(search_dir, f"{record['preprocess_name']}.cif.gz")
            shutil.copyfile(out_path_gz, search_path)  # copying to foldseek directory

            status_dict[record["pocket_id"]] = True

        return status_dict
