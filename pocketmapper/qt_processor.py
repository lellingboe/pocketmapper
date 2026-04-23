"""
Code for processing query/target pairs and orchestrating the workflow
"""

import logging
import os
import pandas as pd
import re


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

    def main(self):
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

        return pd.DataFrame.from_dict(data)

    def parse_individual_qt(self, qt, pocket_method=None):
        """
        Parses input of the form "struct:pocket_info" and returns a dictionary with
        structured information about the structure and pocket. If pocket_method is not provided,
        it will attempt to determine a default pocket method
        """
        data = {
            "pocket_id": qt,
            "struct_info": None,
            "chain_info": None,
            "residue_info": None,
            "struct_type": None,
            "pocket_method": None,
            "success": True,
            "failure_reason": None,
            "sanitized_pocket_id": qt.replace(":", "_").replace(",", "_"),  # For use in file names, etc.
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

        # determining structure info
        data["struct_type"] = self.determine_struct_type(data["struct_info"])
        if pocket_method is not None:
            data["pocket_method"] = pocket_method
        else:
            data["pocket_method"] = self.default_pocket_method(qt)
        return data

    def determine_struct_type(self, struct_str):
        """
        Daterming structure type based on regex pattern matching
        returns one of {"local_file", "pdb", "alphafold"}
        """
        pdb_regex = r"[a-zA-Z0-9]{4}$"
        uniprot_regex = r"([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$"  # https://www.uniprot.org/help/accession_numbers
        if os.path.isfile(struct_str):
            return "local_file"
        elif re.match(pdb_regex, struct_str):
            return "pdb"
        elif re.match(uniprot_regex, struct_str):
            return "alphafold"
        else:
            self.logger.critical(f"Could not determine structure type for {struct_str}", extra=self._log_extra)
            exit(1)

    def default_pocket_method(self, qt_str):
        """
        Daterming pocket method based on regex pattern matching
        returns one of {"pisa", "foldseek", "contact", "other"}
        """
        pocket_info_str = qt_str.split(":", 1)[1]  # Assuming pocket info is always after the first ":"
        pisa_regex = r"^[A-Za-z0-9]_[A-Za-z0-9]$"  # pattern like "A_B"
        chain_chain_regex = r"^[A-Za-z0-9](_[A-Za-z0-9])?\:(\d+\,?)+$"  # pattern like "A_B:1,2,3 or A:1,2,3"
        self.logger.debug(
            f"Determining pocket method for {pocket_info_str} using regex patterns", extra=self._log_extra
        )
        if re.match(pisa_regex, pocket_info_str):
            return "pisa"
        elif re.match(chain_chain_regex, pocket_info_str):
            return "passthrough"
        else:
            self.logger.critical(f"Could not determine pocket method for {pocket_info_str}", extra=self._log_extra)
            exit(1)
