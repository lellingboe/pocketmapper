"""
Download and flatten PDBe PISA interface data into a per-entry cache.

Three stages, each cached on disk so a rerun costs nothing: entry summaries give the assembly ids,
one request per assembly gives its interfaces, and those are flattened into a single
`<pdb_code>.json` per entry keyed by sorted chain pair -- the shape `PisaParser` reads.

Every request is spaced to stay within the PDBe API's tolerance, which is what makes the first run
over a large hit list slow. The spacing grows if the API starts refusing requests and is not lowered
again, so a rate-limited run slows down and stays slow rather than re-provoking the API.
"""

import json
import logging
import os
from collections import defaultdict
from glob import glob

import pandas as pd
from tqdm import tqdm

from pocketmapper.downloads.lib_download import download_api
from pocketmapper.exceptions import PocketMapperError


class PisaDownloader:
    """
    Fetches PISA interfaces from the PDBe API into a local cache.

    `get_interfaces` is the entry point; the remaining methods are its stages and are separately
    usable. Failed downloads are recorded in a `_Failed.txt` beside the files they belong to rather
    than raising, so one dead entry does not abort a large batch. Responses are written through a
    `.part` file, so an interrupted run cannot leave a truncated response for a later run to trust.
    """

    def __init__(self, max_retries=5, base_delay=0.25, max_delay=30.0):
        """
        Configure the retry and rate-limiting behaviour shared by every request.

        Args:
            max_retries (int): Attempts per URL before giving up. Defaults to 5.
            base_delay (float): Seconds between requests, before any backoff. Defaults to 0.25.
            max_delay (float): Ceiling on the doubling backoff delay. Defaults to 30.0.
        """
        self.logger = logging.getLogger(__name__)
        self.stage = {"stage": "PisaDownloader"}
        self.max_retries = max_retries
        self.base_delay = base_delay
        self.max_delay = max_delay

    def fetch_with_backoff(self, url, out_fname):
        """
        Download `url` to `out_fname` under this instance's pacing and retry settings.

        Args:
            url (str): Address to fetch.
            out_fname (str): Path to write the response to.

        Returns:
            bool: True on success, False once the attempts are exhausted.
        """
        return download_api(
            url,
            out_fname,
            max_retries=self.max_retries,
            base_delay=self.base_delay,
            max_delay=self.max_delay,
            log_extra=self.stage,
        )

    def get_interfaces(self, pdb_list, summary_dir, asm_dir, interface_dir):
        """
        Populate the interface cache for a list of PDB entries.

        Entries already cached in `interface_dir` are skipped, so only the missing ones cost requests.
        The remainder are taken through all three stages -- summaries, assemblies, then flattening.

        Args:
            pdb_list (list): PDB codes, in any case.
            summary_dir (str): Cache directory for entry summaries.
            asm_dir (str): Cache directory for per-assembly interface responses.
            interface_dir (str): Output directory for the flattened per-entry files.

        Returns:
            None: Writes one `<pdb_code>.json` per entry into `interface_dir`.
        """
        for dir in [summary_dir, asm_dir, interface_dir]:
            os.makedirs(dir, exist_ok=True)

        existing_files = glob(r"*.json", root_dir=interface_dir)
        missing_pdbs = [x.lower() for x in pdb_list if f"{x.lower()}.json" not in existing_files]

        self.stage = {"stage": "Checking cache for interfaces"}
        if len(missing_pdbs) > 0:
            local_count = len(pdb_list) - len(missing_pdbs)
            logging.info(f"{local_count}/{len(pdb_list)} interfaces found locally", extra=self.stage)

            found_pdbs = self.get_summaries(missing_pdbs, summary_dir)
            assembly_dict = self.parse_summaries(
                found_pdbs, summary_dir
            )  # Dictionary with pdb_code as key and list of assemblies as value
            self.get_assemblies(assembly_dict, asm_dir)
            self.parse_assemblies(assembly_dict, asm_dir, interface_dir)
        else:
            self.logger.info("All interfaces found locally", extra=self.stage)

    def get_summaries(self, pdb_codes, summary_dir):
        """
        Download the PDBe entry summary for each code, which names its assemblies.

        Args:
            pdb_codes (list): Lower-cased PDB codes to fetch.
            summary_dir (str): Directory to cache summaries in; already-present files are not refetched.

        Returns:
            list: The codes whose summary is now on disk. Failures are logged and written to
                `_Failed.txt` in `summary_dir` rather than raising.
        """
        self.stage = {"stage": "Downloading PISA summaries"}
        print("Downloading summaries")
        problems = []
        valid = []
        for pdb_code in tqdm(pdb_codes):
            out_fname = os.path.join(summary_dir, f"{pdb_code}.json")
            if os.path.exists(out_fname):
                valid.append(pdb_code)
            else:
                url = f"https://www.ebi.ac.uk/pdbe/api/v2/pdb/entry/summary/{pdb_code}"
                if self.fetch_with_backoff(url, out_fname):
                    valid.append(pdb_code)
                else:
                    problems.append(pdb_code)
        pd.Series(problems).to_csv(os.path.join(summary_dir, "_Failed.txt"), header=False, index=False)
        return valid

    def parse_summaries(self, pdb_codes, summary_dir):
        """
        Read cached summaries and collect each entry's assembly ids.

        Args:
            pdb_codes (list): Codes whose summaries are cached.
            summary_dir (str): Directory holding them.

        Returns:
            collections.defaultdict: pdb_code -> list of assembly ids.

        Raises:
            PocketMapperError: If a cached summary cannot be parsed.
        """
        self.stage = {"stage": "Parsing PISA summaries"}
        print("Parsing summaries")
        asm_dict = defaultdict(list)
        for pdb_code in tqdm(pdb_codes):
            try:
                fname = os.path.join(summary_dir, f"{pdb_code}.json")
                with open(fname) as f:
                    data = json.load(f)
                logging.debug(f"Summary data for {pdb_code}: {data}", extra=self.stage)
                if len(data[pdb_code]) != 1:
                    logging.critical(f"More than one entry in summary for {pdb_code}", extra=self.stage)
                    continue
                for assembly in data[pdb_code][0]["assemblies"]:
                    asm_dict[pdb_code].append(assembly["assembly_id"])
            except Exception as e:
                logging.exception(f"Issue parsing summary for {pdb_code}", extra=self.stage)
                raise PocketMapperError(f"Issue parsing summary for {pdb_code}: {e}") from e
        return asm_dict

    def get_assemblies(self, asm_dict, asm_dir):
        """
        Download the PISA interfaces for every assembly of every entry.

        One request per assembly, paced by the shared downloader; already-cached assemblies are skipped.

        Args:
            asm_dict (dict): pdb_code -> list of assembly ids, as returned by `parse_summaries`.
            asm_dir (str): Directory to cache the responses in.

        Returns:
            None: Failures are collected into `_Failed.txt` in `asm_dir` rather than raising.
        """
        self.stage = {"stage": "Downloading PISA assemblies"}
        print("Downloading assemblies")
        problems = []
        for pdb_code, assemblies in tqdm(asm_dict.items()):
            for asm in assemblies:
                out_fname = os.path.join(asm_dir, f"{pdb_code}_{asm}.json")
                if not os.path.exists(out_fname):
                    url = f"https://www.ebi.ac.uk/pdbe/api/pisa/interfaces/{pdb_code}/{asm}"
                    if not self.fetch_with_backoff(url, out_fname):
                        problems.append(f"{pdb_code}_{asm}")
        pd.Series(problems).to_csv(os.path.join(asm_dir, "_Failed.txt"), header=False, index=False)

    def parse_assemblies(self, asm_dict, asm_dir, interface_dir):
        """
        Flatten each entry's assemblies into one interface file keyed by chain pair.

        Interfaces from every assembly of an entry are merged into a single dict keyed by the two chain
        ids sorted and concatenated (e.g. "BF"), which is how `PisaParser` looks them up. Where two
        assemblies describe the same chain pair the later one wins.

        Only two-molecule interfaces with single-character chain ids are kept -- a pocket is defined
        against exactly one partner chain, and multi-character ids do not survive the concatenated key.

        Args:
            asm_dict (dict): pdb_code -> list of assembly ids.
            asm_dir (str): Directory holding the cached assembly responses.
            interface_dir (str): Output directory for the per-entry files.

        Returns:
            None: Writes one `<pdb_code>.json` per entry into `interface_dir`.
        """
        self.stage = {"stage": "Parsing PISA assemblies"}
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
                    logging.critical(f"More than one entry in assembly for {pdb_code}_{asm}", extra=self.stage)
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
                                f"More than one molecule in {pdb_code}, {interface['interface_id']}", extra=self.stage
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
                    logging.exception(f"Issue parsing assembly for {pdb_code}_{asm}", extra=self.stage)
                    continue

            interface_fname = os.path.join(interface_dir, f"{pdb_code}.json")
            with open(interface_fname, "w") as out_f:
                json.dump(all_interfaces, out_f)
