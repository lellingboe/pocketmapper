"""
Code for processing query/target pairs and orchestrating the workflow
"""

# TODO Folder input - iterate through files in folder with correct format

import hashlib
import logging
import os
import re
import pandas as pd
import json


class QTProcessor:
    """
    Class for dealing with the procesing of query and target data, including determining types, processing input files, and preparing data for downstream analysis.
    """

    def __init__(self, settings):
        # Initialize logger
        self.logger = logging.getLogger(__name__)
        self._log_extra = {"stage": "QTProcessor Initialization"}
        self.logger.debug(
            "Started",
            extra=self._log_extra,
        )

        self.settings = settings

        # TODO regexes for validating output when pocket_method is specified
        # Structure type regex patterns
        self.pdb_regex = r"[a-zA-Z0-9]{4}$"
        self.uniprot_regex = r"([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$"  # https://www.uniprot.org/help/accession_numbers

        # Pocket method regex patterns
        self.pisa_regex = r"^[A-Za-z0-9]_[A-Za-z0-9]$"  # pattern like "A_B"
        self.passthrough_regex = r"^[A-Za-z0-9](\:(\d+\,?)*)?$"  # pattern like "A(:1,2,3)"
        self.vdw_regex = r"^[A-Za-z0-9](_[A-Za-z0-9])?(\:(\d+\,?)*)?$"  # pattern like "A_B:1,2,3"

        self._bundled_foldseek_dbs = {"human_domains"}
        self._available_pocket_methods = {"pisa", "passthrough", "vdw"}

    def process_qt_cmdline_input(self):
        """
        qt_input: input from parsing query/target on command line
        """
        dfs = []
        for name in ["query", "target"]:
            self._log_extra.update({"stage": f"Processing {name}"})
            self.logger.debug("Processing {name}", extra=self._log_extra)

            qt_input = self.settings[name]
            pocket_method = self.settings[f"{name}_pocket_method"]
            self.processing_target = name == "Target"

            # Check that query and target are specified
            if isinstance(qt_input, type(None)):
                self.logger.critical(f"{name} input is required. Exiting.", extra=self._log_extra)
                exit(1)

            data = []
            if os.path.isfile(qt_input):  # if it's a file, process each line as a separate query/target
                try:
                    with open(qt_input) as f:
                        for line in f.readlines():
                            data.append(self.parse_individual_qt(line.strip(), pocket_method=pocket_method))
                except Exception as e:
                    self.logger.critical(f"Problem reading the file {qt_input}: {e}", extra=self._log_extra)
                    exit(1)
            else:
                data.append(self.parse_individual_qt(qt_input, pocket_method=pocket_method))

            data = [
                x for x in data if x is not None
            ]  # removing any None entries that may have been added due to errors
            dfs.append(pd.DataFrame.from_dict(data))
        return dfs[0], dfs[1]  # return query and target dataframes

    def parse_individual_qt(self, qt, pocket_method):
        """
        Parses input of the form "struct_info:chain_info:residues" and returns a dictionary with
        structured information about the structure and pocket. If pocket_method is not provided,
        it will attempt to determine a default pocket method
        """
        data = {
            "pocket_id": qt,
            "struct_info": None,
            "chain_info": None,
            "residue_info": None,
            "struct_type": None,
            "struct_path": None,
            "preprocess_name": None,
            "preprocess_path": None,
            "pocket_method": None,
            "success": True,
            "failure_reason": "",
        }

        # Bundled foldseek databases have a special format and are treated differently
        if qt in self._bundled_foldseek_dbs:
            data["struct_info"] = qt
            data["struct_type"] = "foldseek_db"
            data["struct_path"] = "foldseek_db"
            return data

        # Unpack the input string into its components
        categories = ["struct_info", "chain_info", "residue_info"]
        for x, y in zip(categories, qt.split(":")):
            data[x] = y

        # If chain info is not provided, default to chain A for targets
        if data["chain_info"] is None:
            self.logger.warning(f"No specified chain for {qt}, skipping", extra=self._log_extra)
            return None

        # determining structure info
        data["struct_type"] = self.determine_struct_type(data["struct_info"])
        if data["struct_type"] is None:
            self.logger.warning(f"Could not determine structure type for {qt}", extra=self._log_extra)
            return None
        data["struct_path"] = self.determine_ref_struct_path(data["struct_info"], data["struct_type"])

        # Generate a unique name for the structure based on its components and a hash of the name
        input_fname = os.path.basename(data["struct_info"]).split(".")[0]
        name = input_fname + "_" + data["chain_info"][0]  # e.g., "P12345_A" or "1ABC_A"
        name_md5 = hashlib.md5(name.encode()).hexdigest()
        data["preprocess_name"] = name + "_" + name_md5
        data["preprocess_path"] = os.path.join(
            self.settings["foldseek_preprocessed_structure_dir"], data["preprocess_name"] + ".cif"
        )
        data["preprocess_path_gz"] = data["preprocess_path"] + ".gz"

        if pocket_method is not None:
            data["pocket_method"] = pocket_method
        else:
            data["pocket_method"] = self.determine_pocket_method(qt, data["struct_type"])
        if data["pocket_method"] is None:
            self.logger.warning(f"Could not determine pocket method for {qt}", extra=self._log_extra)
            return None

        self.logger.debug(f"Processed {qt} into structured data: {json.dumps(data, indent=4)}", extra=self._log_extra)
        return data

    def determine_struct_type(self, struct_str):
        """
        If struct_str is a file:
            return "local_file"
        If struct_str matches a PDB ID regex pattern:
            return "pdb"
        If struct_str matches a Uniprot ID regex pattern:
            return "alphafold"
        else:
            log critical error and exit
        """
        if re.match(self.pdb_regex, struct_str):
            return "pdb"
        elif re.match(self.uniprot_regex, struct_str):
            return "alphafold"
        elif os.path.isfile(struct_str):
            return "local_file"
        elif os.isdir(struct_str):
            logging.critical(f"Directory input is not currently supported: {struct_str}", extra=self._log_extra)
            exit(1)
        else:
            logging.warning(f"Could not determine structure type for {struct_str}", extra=self._log_extra)
            return None

    def determine_ref_struct_path(self, struct_info, struct_type):
        """
        Determine the path to the structure file based on its type and identifier.

        Args:
            struct_info (str): Identifier for the structure (e.g., "P12345", "1ABC").
            struct_type (str): Type of the structure ("alphafold", "pdb", "local_file").

        Returns:
            str: Path to the structure file.
        """
        match struct_type:
            case "alphafold":
                return os.path.join(self.settings["structure_dir"], f"{struct_info}.cif.gz")
            case "pdb":
                return os.path.join(self.settings["structure_dir"], f"{struct_info}.cif.gz")
            case "local_file":
                return struct_info
            case _:
                logging.critical(
                    f"Unknown structure type {struct_type} for struct_info {struct_info}", extra=self._log_extra
                )
                exit(1)

    def determine_pocket_method(self, qt_str, struct_type):
        """
        Daterming pocket method based on regex pattern matching
        returns one of {"pisa", "foldseek", "contact", "other"}
        """
        pocket_info_str = qt_str.split(":", 1)[1]  # Assuming pocket info is always after the first ":"
        self.logger.debug(
            f"Determining pocket method for {pocket_info_str} using regex patterns", extra=self._log_extra
        )
        match struct_type:
            case "alphafold":
                if re.match(self.passthrough_regex, pocket_info_str):
                    return "passthrough"
                else:
                    return None
            case "pdb":
                if re.match(self.pisa_regex, pocket_info_str):
                    return "pisa"
                elif re.match(self.passthrough_regex, pocket_info_str):
                    return "passthrough"
                elif re.match(self.vdw_regex, pocket_info_str):
                    return "vdw"
                else:
                    return None
            case "local_file":
                if re.match(self.passthrough_regex, pocket_info_str):
                    return "passthrough"
                elif re.match(self.vdw_regex, pocket_info_str):
                    return "vdw"
                else:
                    return None
            case _:
                return None
