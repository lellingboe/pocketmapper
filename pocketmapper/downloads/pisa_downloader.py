"""
Download and flatten PDBe PISA interface data into a per-entry cache.

Three stages, each cached on disk so a rerun costs nothing: entry summaries give the assembly ids,
one request per assembly gives its interfaces, and those are flattened into a single
`<pdb_code>.json` per entry keyed by sorted chain pair -- the shape `PisaParser` reads.

Every request is spaced to stay within the PDBe API's tolerance, which is what makes the first run
over a large hit list slow. The spacing grows if the API starts refusing requests and is not lowered
again, so a rate-limited run slows down and stays slow rather than re-provoking the API.

No stage raises on a bad entry. Each returns the ids it could not handle and the entry point
collects them into one JSON report, so the entries PISA lacks cost a line in that file rather than
the rest of the batch.
"""

import json
import logging
import os
from collections import defaultdict
from glob import glob

from tqdm import tqdm

from pocketmapper.downloads.lib_download import download_api
from pocketmapper.exceptions import PocketMapperError


class PisaDownloader:
    """
    Fetches PISA interfaces from the PDBe API into a local cache.

    `download_missing_interfaces` is the entry point; the remaining methods are its stages and are
    separately usable, though each expects its directory to exist already. Every stage returns what
    it could not handle instead of raising, so one dead entry does not abort a large batch.
    Responses are written through a `.part` file, so an interrupted run cannot leave a truncated
    response for a later run to trust.
    """

    def __init__(self, max_retries=5, base_delay=0.25, max_delay=30.0):
        """
        Configure the retry and rate-limiting behaviour shared by every request.

        Args:
            max_retries (int): Attempts per URL before giving up. Defaults to 5.
            base_delay (float): Seconds between requests, before any backoff. Defaults to 0.25.
            max_delay (float): Ceiling on the doubling backoff delay. Defaults to 30.0.
        """
        self.log_extra = {"stage": "Downloading PISA Interfaces"}
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
            log_extra=self.log_extra,
        )

    def download_missing_interfaces(self, pdb_list, summary_dir, asm_dir, interface_dir, error_path):
        """
        Populate the interface cache for a list of PDB entries.

        Entries already cached in `interface_dir` are skipped, so only the missing ones cost requests.
        The remainder are taken through all four stages -- summaries downloaded then parsed for the
        assembly ids, assemblies downloaded then flattened. The three directories are created if they
        do not exist.

        Args:
            pdb_list (list): PDB codes, in any case.
            summary_dir (str): Cache directory for entry summaries.
            asm_dir (str): Cache directory for per-assembly interface responses.
            interface_dir (str): Output directory for the flattened per-entry files.
            error_path (str): Path for the failure report, written only if something failed.

        Returns:
            dict: Stage name -> the ids that stage could not handle, keyed `summary_downloading`,
                `summary_parsing`, `assembly_downloading` and `assembly_parsing`. A stage with
                nothing to report is absent, so an empty dict means a clean run. Writes one
                `<pdb_code>.json` per entry into `interface_dir`.
        """

        existing_files = glob(r"*.json", root_dir=interface_dir)
        missing_pdbs = [x.lower() for x in pdb_list if f"{x.lower()}.json" not in existing_files]
        local_count = len(pdb_list) - len(missing_pdbs)
        logging.info(f"{local_count}/{len(pdb_list)} interfaces found locally", extra=self.log_extra)

        all_failures = {}
        if len(missing_pdbs) > 0:
            os.makedirs(summary_dir, exist_ok=True)
            os.makedirs(asm_dir, exist_ok=True)
            os.makedirs(interface_dir, exist_ok=True)

            success, sd_failures = self.download_missing_summaries(missing_pdbs, summary_dir)
            if sd_failures:
                all_failures["summary_downloading"] = sd_failures
            assembly_dict, sp_failures = self.parse_summaries(success, summary_dir)
            if sp_failures:
                all_failures["summary_parsing"] = sp_failures

            success, ad_failures = self.download_missing_assemblies(assembly_dict, asm_dir)
            if ad_failures:
                all_failures["assembly_downloading"] = ad_failures
            ap_failures = self.parse_assemblies(success, asm_dir, interface_dir)
            if ap_failures:
                all_failures["assembly_parsing"] = ap_failures

            if all_failures:
                with open(error_path, "w") as f:
                    json.dump(all_failures, f)
                logging.warning(
                    f"See {error_path} for details of {sum(len(v) for v in all_failures.values())} PISA failures",
                    extra=self.log_extra,
                )

        return all_failures

    def download_missing_summaries(self, pdb_codes, summary_dir):
        """
        Download the PDBe entry summary for each code, which names its assemblies.

        Args:
            pdb_codes (list): Lower-cased PDB codes to fetch.
            summary_dir (str): Directory to cache summaries in; already-present files are not refetched.

        Returns:
            tuple: (list of codes whose summary is now on disk, list of codes that could not be
                fetched). A failure is logged and collected rather than raised.
        """
        logging.info(f"Downloading {len(pdb_codes)} PISA summaries", extra=self.log_extra)
        success = []
        failure = []
        for pdb_code in tqdm(pdb_codes):
            out_fname = os.path.join(summary_dir, f"{pdb_code}.json")
            if os.path.exists(out_fname):
                success.append(pdb_code)
            else:
                url = f"https://www.ebi.ac.uk/pdbe/api/v2/pdb/entry/summary/{pdb_code}"
                if self.fetch_with_backoff(url, out_fname):
                    success.append(pdb_code)
                else:
                    failure.append(pdb_code)
        if len(failure) > 0:
            logging.warning(f"Failed to download {len(failure)} summaries", extra=self.log_extra)
        return success, failure

    def parse_summaries(self, pdb_codes, summary_dir):
        """
        Read cached summaries and collect each entry's assembly ids.

        Args:
            pdb_codes (list): Codes whose summaries are cached.
            summary_dir (str): Directory holding them.

        Returns:
            tuple: (defaultdict of pdb_code -> list of assembly ids, list of codes whose summary
                could not be read). A summary naming anything other than exactly one entry counts as
                a failure; neither case raises.
        """
        logging.info(f"Parsing {len(pdb_codes)} PISA summaries", extra=self.log_extra)
        asm_dict = defaultdict(list)
        failure = []
        for pdb_code in tqdm(pdb_codes):
            try:
                fname = os.path.join(summary_dir, f"{pdb_code}.json")
                with open(fname) as f:
                    data = json.load(f)
                logging.debug(f"Summary data for {pdb_code}: {data}", extra=self.log_extra)
                if len(data[pdb_code]) != 1:
                    raise PocketMapperError(f"More than one entry in summary for {pdb_code}")
                for assembly in data[pdb_code][0]["assemblies"]:
                    asm_dict[pdb_code].append(assembly["assembly_id"])
            except Exception as e:
                logging.debug(f"Issue parsing summary for {pdb_code}: ({e})", extra=self.log_extra)
                failure.append(pdb_code)
        if len(failure) > 0:
            logging.warning(f"Failed to parse {len(failure)} PISA summaries", extra=self.log_extra)
        return asm_dict, failure

    def download_missing_assemblies(self, asm_dict, asm_dir):
        """
        Download the PISA interfaces for every assembly of every entry.

        One request per assembly, paced by the shared downloader. An assembly already on disk counts
        as a success without costing a request.

        Args:
            asm_dict (dict): pdb_code -> list of assembly ids, as returned by `parse_summaries`.
            asm_dir (str): Directory to cache the responses in.

        Returns:
            tuple: (defaultdict of pdb_code -> list of the assembly ids now on disk, list of
                `<pdb_code>_<assembly_id>` that could not be fetched).
        """
        logging.info(f"Downloading {sum(len(v) for v in asm_dict.values())} PISA assemblies", extra=self.log_extra)
        failure = []
        success = defaultdict(list)
        for pdb_code, assemblies in tqdm(asm_dict.items()):
            for asm in assemblies:
                out_fname = os.path.join(asm_dir, f"{pdb_code}_{asm}.json")
                url = f"https://www.ebi.ac.uk/pdbe/api/pisa/interfaces/{pdb_code}/{asm}"
                if os.path.exists(out_fname) or self.fetch_with_backoff(url, out_fname):
                    success[pdb_code].append(asm)
                else:
                    failure.append(f"{pdb_code}_{asm}")
        if len(failure) > 0:
            logging.warning(f"Failed to download {len(failure)} assemblies", extra=self.log_extra)
        return success, failure

    def parse_assemblies(self, asm_dict, asm_dir, interface_dir):
        """
        Flatten each entry's assemblies into one interface file keyed by chain pair.

        Interfaces from every assembly of an entry are merged into a single dict keyed by the two chain
        ids sorted and concatenated (e.g. "BF"), which is how `PisaParser` looks them up. Where two
        assemblies describe the same chain pair the later one wins, and an entry left with no usable
        interface gets no file at all.

        Only two-molecule interfaces with single-character chain ids are kept -- a pocket is defined
        against exactly one partner chain, and multi-character ids do not survive the concatenated key.
        Skipping one is routine, so it is reported alongside the genuine failures rather than raised.

        Args:
            asm_dict (dict): pdb_code -> list of assembly ids.
            asm_dir (str): Directory holding the cached assembly responses.
            interface_dir (str): Output directory for the per-entry files.

        Returns:
            list: What did not make it into a file -- `<pdb_code>_<assembly_id>` for an assembly that
                could not be read, `<pdb_code>_<assembly_id>_parse_error` for one that raised, and
                `<pdb_code>_<assembly_id>_<interface_id>` for each skipped interface. Writes one
                `<pdb_code>.json` per entry that has at least one usable interface.
        """
        logging.info(
            f"Parsing {sum(len(v) for v in asm_dict.values())} PISA assemblies into per-interface files",
            extra=self.log_extra,
        )
        failure = []
        for pdb_code, assemblies in tqdm(asm_dict.items()):
            pdb_interfaces = {}
            for asm_id in assemblies:
                try:
                    # Opening assembly file
                    asm_fname = os.path.join(asm_dir, f"{pdb_code}_{asm_id}.json")
                    with open(asm_fname) as f:
                        data = json.load(f)

                    if len(data.keys()) != 1:
                        logging.debug(f"More than one entry in assembly for {pdb_code}_{asm_id}", extra=self.log_extra)
                        failure.append(f"{pdb_code}_{asm_id}")
                        continue

                    if "PISA" in data:
                        data = data["PISA"]
                    else:
                        data = data[pdb_code]

                    for interface in data["assembly"]["interfaces"]:
                        interface_id = interface["interface_id"]
                        # Checking if interface is between two molecules
                        if len(interface["molecules"]) != 2:
                            logging.debug(
                                f"More than one molecule in {pdb_code}, {interface['interface_id']}",
                                extra=self.log_extra,
                            )
                            failure.append(f"{pdb_code}_{asm_id}_{interface_id}")
                            continue
                        # Checking the chain ids are single character
                        # TODO - Find a way to handle multi character chain ids
                        chain_ids = []
                        for molecule in interface["molecules"]:
                            chain_ids.append(molecule["chain_id"])
                        if not all(len(c) == 1 for c in chain_ids):
                            logging.debug(
                                f"Multi-character chain ids in {pdb_code}, {interface['interface_id']}",
                                extra=self.log_extra,
                            )
                            failure.append(f"{pdb_code}_{asm_id}_{interface_id}")
                            continue
                        entry_name = "".join(sorted(chain_ids))
                        pdb_interfaces[entry_name] = interface
                except Exception as e:
                    logging.debug(
                        f"Unexpected error parsing assembly for {pdb_code}_{asm_id}: ({e})", extra=self.log_extra
                    )
                    failure.append(f"{pdb_code}_{asm_id}_parse_error")

            if len(pdb_interfaces) > 0:
                interface_fname = os.path.join(interface_dir, f"{pdb_code}.json")
                with open(interface_fname, "w") as out_f:
                    json.dump(pdb_interfaces, out_f)
        if len(failure) > 0:
            logging.warning(f"Failed to parse {len(failure)} PISA assemblies", extra=self.log_extra)
        return failure
