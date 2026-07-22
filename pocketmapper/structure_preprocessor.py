import gzip
import logging
import os
import shutil
import gemmi
from tqdm import tqdm


class StructurePreprocessor:
    def __init__(
        self,
    ):
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
        """
        self.cache = os.listdir(self.out_dir)

    def preprocess_records(self, records, search_dir):
        """
        Docstring for pdb_preprocessing_gemmi

        :param df: Description
        :param ref_dir: directory for reference pdb files to be divided
        :param cache_dir: directory for divided pdbs to be cached
        :param out_dir: directory to be used with foldseek
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

            # except Exception as e:
            #    logging.warning(f"Could not divide {struct_info} with chain info {chain_info}", extra=stage)
            #    logging.debug("Exception info", exc_info=e, extra=stage)
            #    status_dict[record["pocket_id"]] = False

        return status_dict
