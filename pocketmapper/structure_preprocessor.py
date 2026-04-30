import gzip
import logging
import os
import shutil
import gemmi
from tqdm import tqdm


class StructurePreprocessor:
    def __init__(self, source_dir, out_dir):
        self.logger = logging.getLogger(__name__)
        self._log_extra = {"stage": "Structure Preprocessor"}
        self.logger.debug("Initialized", extra=self._log_extra)

        self.source_dir = source_dir
        self.out_dir = out_dir
        self.cache = os.listdir(out_dir)

    def update_cache(self):
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
            out_struct_name = struct_info + "_" + chain
            out_path = os.path.join(self.out_dir, f"{out_struct_name}.cif")
            out_path_gz = out_path + ".gz"

            if not (out_path_gz in self.cache):
                match record["struct_type"]:
                    case "alphafold":
                        ref_path = os.path.join(self.source_dir, f"{struct_info}.cif.gz")
                    case "pdb":
                        ref_path = os.path.join(self.source_dir, f"{struct_info}.cif.gz")
                    case "local_file":
                        # TODO Implement local file handling
                        raise NotImplementedError("Local file handling not implemented yet")
                    case _:
                        self.logger.warning(
                            f"Unknown structure type {record['struct_type']} for struct_info {struct_info}",
                            extra=stage,
                        )

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
                    return status_dict

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

            search_path = os.path.join(search_dir, f"{out_struct_name}.cif.gz")
            shutil.copyfile(out_path_gz, search_path)  # copying to foldseek directory

            status_dict[record["pocket_id"]] = True

            # except Exception as e:
            #    logging.warning(f"Could not divide {struct_info} with chain info {chain_info}", extra=stage)
            #    logging.debug("Exception info", exc_info=e, extra=stage)
            #    status_dict[record["pocket_id"]] = False

        return status_dict
