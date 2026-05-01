"""
Code for processing query/target pairs and orchestrating the workflow
"""

# TODO Folder input - iterate through files in folder with correct format

import logging
import os
import re
import pandas as pd


class QTProcessor:
    """
    Class for dealing with the procesing of query and target data, including determining types, processing input files, and preparing data for downstream analysis.
    """

    def __init__(self, query, target, query_pocket_method, target_pocket_method):
        self._available_pocket_methods = {"pisa", "passthrough", "vdw"}

        # Store the input parameters
        self._query = query
        self._target = target
        self._query_pocket_method = query_pocket_method
        self._target_pocket_method = target_pocket_method

        # Initialize logger
        self.logger = logging.getLogger(__name__)
        self._log_extra = {"stage": "QTProcessor Initialization"}
        self.logger.debug(
            f"\n  query: {query}\n  target: {target}\n  query_pocket_method: {query_pocket_method}\n  target_pocket_method: {target_pocket_method}",
            extra=self._log_extra,
        )

        # Check that query and target are specified
        if isinstance(query, type(None)) or isinstance(target, type(None)):
            self.logger.critical(
                f"Query and/or target not specified:\nquery: {query}\ntarget: {target}", extra=self._log_extra
            )
            exit(1)

        # TODO regexes for validating output when pocket_method is specified
        # setting regex patterns for determining default structure types
        self.pdb_regex = r"[a-zA-Z0-9]{4}$"
        self.uniprot_regex = r"([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$"  # https://www.uniprot.org/help/accession_numbers
        # setting regex patterns for determining default pocket methods
        self.pisa_regex = r"^[A-Za-z0-9]_[A-Za-z0-9]$"  # pattern like "A_B"
        self.passthrough_regex = r"^[A-Za-z0-9](\:(\d+\,?)*)?$"  # pattern like "A(:1,2,3)"
        self.vdw_regex = r"^[A-Za-z0-9](_[A-Za-z0-9])?(\:(\d+\,?)*)?$"  # pattern like "A_B:1,2,3"

    def process_qt(self):
        """
        Main method to process the query and target data. This includes determining the types of
        query and target, processing the input data accordingly, and preparing it for downstream analysis.
        """
        self._log_extra.update({"stage": "Main QT Processing"})
        self.logger.debug("Starting main", extra=self._log_extra)

        # Processing query/target
        processed_query = self.process_qt_cmdline_input(self._query, pocket_method=self._query_pocket_method)
        processed_target = self.process_qt_cmdline_input(self._target, pocket_method=self._target_pocket_method)

        # Logging the processed data and determined pocket methods for debugging purposes
        self.logger.debug(f"Processed query: \n{processed_query}", extra=self._log_extra)
        self.logger.debug(f"Processed target: \n{processed_target}", extra=self._log_extra)

        return processed_query, processed_target

    def process_qt_cmdline_input(self, qt_input, pocket_method=None):
        """
        qt_input: input from parsing query/target on command line
        """
        self._log_extra.update({"stage": "Processing Query/Target"})
        self.logger.debug("Starting process_qt_cmdline_input", extra=self._log_extra)

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

        data = [x for x in data if x is not None]  # removing any None entries that may have been added due to errors
        return pd.DataFrame.from_dict(data)

    def parse_individual_qt(self, qt, pocket_method=None):
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
            "preprocess_name": None,
            "pocket_method": None,
            "success": True,
            "failure_reason": "",
        }

        # Bundled foldseek databases have a special format and are treated differently
        if qt == "human_domains":
            data["struct_info"] = "human_domains"
            data["struct_type"] = "foldseek_db"
            return data

        # Unpack the input string into its components
        categories = ["struct_info", "chain_info", "residue_info"]
        for x, y in zip(categories, qt.split(":")):
            data[x] = y
        if data["chain_info"] is None:
            self.logger.warning(f"Chain info not specified in {qt}", extra=self._log_extra)
            return None

        input_fname = os.path.basename(data["struct_info"]).split(".")[0]
        data["preprocess_name"] = input_fname + "_" + data["chain_info"][0]  # e.g., "P12345_A" or "1ABC_A"

        # determining structure info
        data["struct_type"] = self.determine_struct_type(data["struct_info"])
        if data["struct_type"] is None:
            self.logger.warning(f"Could not determine structure type for {qt}", extra=self._log_extra)
            return None

        if pocket_method is not None:
            data["pocket_method"] = pocket_method
        else:
            data["pocket_method"] = self.determine_pocket_method(qt, data["struct_type"])
        if data["pocket_method"] is None:
            self.logger.warning(f"Could not determine pocket method for {qt}", extra=self._log_extra)
            return None

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
        if os.path.isfile(struct_str):
            return "local_file"
        elif re.match(self.pdb_regex, struct_str):
            return "pdb"
        elif re.match(self.uniprot_regex, struct_str):
            return "alphafold"
        else:
            return None

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
