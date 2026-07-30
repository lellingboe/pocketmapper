"""
Code related to downloading and parsing PISA interfaces
"""

import os
import logging
from urllib.request import urlcleanup, urlretrieve
import pandas as pd
from time import sleep
from tqdm import tqdm
import json
from collections import defaultdict
from glob import glob
from pocketmapper.exceptions import PocketMapperError


class PisaDownloader:
    def __init__(self, max_retries=5, base_delay=0.25, max_delay=30.0):
        self.logger = logging.getLogger(__name__)
        self._stage = {}
        self.max_retries = max_retries
        self.base_delay = base_delay
        self.max_delay = max_delay

    def _fetch_with_backoff(self, url, out_fname):
        """
        Download url to out_fname, retrying on failure with exponential backoff.

        Delay starts at base_delay and doubles after each failed attempt, capped
        at max_delay. Returns True on success, False once max_retries is exhausted.
        """
        delay = self.base_delay
        for attempt in range(1, self.max_retries + 1):
            try:
                urlcleanup()
                urlretrieve(url, out_fname)
                return True
            except Exception as e:
                if attempt == self.max_retries:
                    logging.warning(f"Giving up on {url} after {self.max_retries} attempts: {e}", extra=self._stage)
                    return False
                logging.debug(
                    f"Attempt {attempt}/{self.max_retries} failed for {url} ({e}); retrying in {delay:.2f}s",
                    extra=self._stage,
                )
                sleep(delay)
                delay = min(delay * 2, self.max_delay)
        return False

    def get_interfaces(self, pdb_list, summary_dir, asm_dir, interface_dir):
        for dir in [summary_dir, asm_dir, interface_dir]:
            os.makedirs(dir, exist_ok=True)

        existing_files = glob(r"*.json", root_dir=interface_dir)
        missing_pdbs = [x.lower() for x in pdb_list if f"{x.lower()}.json" not in existing_files]

        self._stage = {"stage": "Checking cache for interfaces"}
        if len(missing_pdbs) > 0:
            local_count = len(pdb_list) - len(missing_pdbs)
            logging.info(f"{local_count}/{len(pdb_list)} interfaces found locally", extra=self._stage)

            found_pdbs = self.get_summaries(missing_pdbs, summary_dir)
            assembly_dict = self.parse_summaries(
                found_pdbs, summary_dir
            )  # Dictionary with pdb_code as key and list of assemblies as value
            self.get_assemblies(assembly_dict, asm_dir)
            self.parse_assemblies(assembly_dict, asm_dir, interface_dir)
        else:
            self.logger.info("All interfaces found locally", extra=self._stage)

    def get_summaries(self, pdb_codes, summary_dir):
        self._stage = {"stage": "Downloading PISA summaries"}
        print("Downloading summaries")
        problems = []
        valid = []
        for pdb_code in tqdm(pdb_codes):
            out_fname = os.path.join(summary_dir, f"{pdb_code}.json")
            if os.path.exists(out_fname):
                valid.append(pdb_code)
            else:
                url = f"https://www.ebi.ac.uk/pdbe/api/v2/pdb/entry/summary/{pdb_code}"
                if self._fetch_with_backoff(url, out_fname):
                    valid.append(pdb_code)
                else:
                    problems.append(pdb_code)
                sleep(self.base_delay)
        pd.Series(problems).to_csv(os.path.join(summary_dir, "_Failed.txt"), header=False, index=False)
        return valid

    def parse_summaries(self, pdb_codes, summary_dir):
        self._stage = {"stage": "Parsing PISA summaries"}
        print("Parsing summaries")
        asm_dict = defaultdict(list)
        for pdb_code in tqdm(pdb_codes):
            try:
                fname = os.path.join(summary_dir, f"{pdb_code}.json")
                with open(fname) as f:
                    data = json.load(f)
                logging.debug(f"Summary data for {pdb_code}: {data}", extra=self._stage)
                if len(data[pdb_code]) != 1:
                    logging.critical(f"More than one entry in summary for {pdb_code}", extra=self._stage)
                    continue
                for assembly in data[pdb_code][0]["assemblies"]:
                    asm_dict[pdb_code].append(assembly["assembly_id"])
            except Exception as e:
                logging.exception(f"Issue parsing summary for {pdb_code}", extra=self._stage)
                raise PocketMapperError(f"Issue parsing summary for {pdb_code}: {e}") from e
        return asm_dict

    def get_assemblies(self, asm_dict, asm_dir):
        self._stage = {"stage": "Downloading PISA assemblies"}
        print("Downloading assemblies")
        problems = []
        for pdb_code, assemblies in tqdm(asm_dict.items()):
            for asm in assemblies:
                out_fname = os.path.join(asm_dir, f"{pdb_code}_{asm}.json")
                if not os.path.exists(out_fname):
                    url = f"https://www.ebi.ac.uk/pdbe/api/pisa/interfaces/{pdb_code}/{asm}"
                    if not self._fetch_with_backoff(url, out_fname):
                        problems.append(f"{pdb_code}_{asm}")
                    sleep(self.base_delay)
        pd.Series(problems).to_csv(os.path.join(asm_dir, "_Failed.txt"), header=False, index=False)

    def parse_assemblies(self, asm_dict, asm_dir, interface_dir):
        self._stage = {"stage": "Parsing PISA assemblies"}
        print("Parsing assemblies")
        for pdb_code, assemblies in tqdm(asm_dict.items()):

            all_interfaces = {}
            for asm in assemblies:
                asm_fname = os.path.join(asm_dir, f"{pdb_code}_{asm}.json")
                # Opening assembly file
                try:
                    with open(asm_fname) as f:
                        data = json.load(f)
                except FileNotFoundError:
                    continue
                if len(data.keys()) != 1:
                    logging.critical(f"More than one entry in assembly for {pdb_code}_{asm}", extra=self._stage)
                    continue

                if "PISA" in data:
                    data = data["PISA"]
                    pdb_code = data["pdb_id"]
                else:
                    data = data[pdb_code]

                try:
                    for interface in data["assembly"]["interfaces"]:
                        # Checking if interface is between two molecules
                        if len(interface["molecules"]) != 2:
                            logging.critical(
                                f"More than one molecule in {pdb_code}, {interface["interface_id"]}", extra=self._stage
                            )
                            continue
                        # Checking the chain meet the expected specifications
                        chain_ids = []
                        for molecule in interface["molecules"]:
                            chain_ids.append(molecule["chain_id"])
                        if not all(len(c) == 1 for c in chain_ids):
                            continue
                        entry_name = "".join(sorted(chain_ids))
                        all_interfaces[entry_name] = interface
                except Exception:
                    logging.exception(f"Issue parsing assembly for {pdb_code}_{asm}", extra=self._stage)
                    continue

            interface_fname = os.path.join(interface_dir, f"{pdb_code}.json")
            with open(interface_fname, "w") as out_f:
                json.dump(all_interfaces, out_f)
