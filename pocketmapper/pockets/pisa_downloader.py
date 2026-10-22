"""
Download and flatten PDBe PISA interface data into a per-entry cache.

Three stages, each cached on disk so a rerun costs nothing: entry summaries, fetched in batches,
give the assembly ids, one file per assembly gives its interfaces, and those are flattened into a single
`<pdb_code>.json` per entry keyed by sorted chain pair -- the shape `PisaParser` reads.

`PisaDownloader` names the cache's contents under the `pisa_dir` it is given: `summaries/`,
`assemblies/`, `interface_pairs/` and `errors.json`.

The per-assembly files come from one of two sources serving identical JSON, so both share one cache:
the EBI FTP server, fetched concurrently, or the PDBe API, one request per assembly.

Every API request is spaced to stay within the PDBe API's tolerance, which is what makes the first
API run over a large hit list slow. The spacing grows if the API starts refusing requests and is not
lowered again, so a rate-limited run slows down and stays slow rather than re-provoking the API.

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

from pocketmapper.lib.constants import DEFAULT_PISA_SOURCE, PISA_SOURCES
from pocketmapper.lib.download import DOWNLOAD_WORKERS, DownloadJob, download_api, download_file, download_many
from pocketmapper.lib.exceptions import PocketMapperError

logger = logging.getLogger(__name__)


def ftp_assembly_url(pdb_code, asm):
    """
    The EBI FTP address of one assembly's interfaces.

    Args:
        pdb_code (str): Lower-cased PDB code.
        asm (str): Assembly id.

    Returns:
        str: The URL.
    """
    return (
        "https://ftp.ebi.ac.uk/pub/databases/msd/pdb-assemblies-analysis/split/"
        f"{pdb_code[1:3]}/{pdb_code}_assembly{asm}_interfaces.json"
    )


def api_assembly_url(pdb_code, asm):
    """
    The PDBe API address of one assembly's interfaces.

    Args:
        pdb_code (str): Lower-cased PDB code.
        asm (str): Assembly id.

    Returns:
        str: The URL.
    """
    return f"https://www.ebi.ac.uk/pdbe/api/pisa/interfaces/{pdb_code}/{asm}"


# Keys match lib.constants.PISA_SOURCES
ASSEMBLY_SOURCES = {"ftp": (ftp_assembly_url, download_file), "api": (api_assembly_url, download_api)}


class PisaDownloader:
    """
    Fetches PISA interfaces from the EBI FTP server or the PDBe API into a local cache.

    `download_missing_interfaces` is the entry point; the remaining methods are its stages and are
    separately usable. Only `parse_assemblies` expects its output directory, `interface_dir`, to exist
    already. Every stage returns what it could not handle instead of raising, so one dead entry does not
    abort a large batch. Responses are written through a `.part` file, so an interrupted run cannot leave
    a truncated response for a later run to trust.
    """

    def __init__(
        self,
        pisa_dir,
        source=DEFAULT_PISA_SOURCE,
        summary_batch_size=50,
        max_workers=DOWNLOAD_WORKERS,
    ):
        """
        Lay out the cache under `pisa_dir` and configure where assembly interfaces come from.

        Args:
            pisa_dir (str): Cache directory. Holds `summaries/`, `assemblies/`, `interface_pairs/` (the
                per-entry files `PisaParser` reads) and the failure report `errors.json`.
            source (str, optional): Where assembly interfaces are fetched from: "ftp" or "api". Summaries
                always come from the API. Defaults to DEFAULT_PISA_SOURCE.
            summary_batch_size (int, optional): Entries per summary request. The API fails with a 500
                somewhere above 800, so keep it well below that. Defaults to 50.
            max_workers (int, optional): Concurrent FTP downloads. Defaults to DOWNLOAD_WORKERS.

        Raises:
            PocketMapperError: If `source` is not one of PISA_SOURCES.
        """
        self.log_extra = {"stage": "Downloading PISA Interfaces"}
        if source not in PISA_SOURCES:
            msg = f"Unknown PISA source {source!r}. Choose one of: {', '.join(PISA_SOURCES)}."
            logger.critical(msg, extra=self.log_extra)
            raise PocketMapperError(msg)
        self.source = source
        self.summary_batch_size = summary_batch_size
        self.max_workers = max_workers
        self.summary_dir = os.path.join(pisa_dir, "summaries")
        self.asm_dir = os.path.join(pisa_dir, "assemblies")
        self.interface_dir = os.path.join(pisa_dir, "interface_pairs")
        self.error_path = os.path.join(pisa_dir, "errors.json")

    def download_missing_interfaces(self, pdb_list):
        """
        Populate the interface cache for a list of PDB entries.

        Entries already in `interface_dir` cost no requests. For the rest, caches the summaries and
        assembly responses, writes one `<pdb_code>.json` per entry into `interface_dir`, and writes
        `error_path` if anything failed.

        Args:
            pdb_list (list): PDB codes, in any case.

        Returns:
            dict: Stage name -> the ids that stage could not handle, keyed `summary_downloading`,
                `summary_parsing`, `assembly_downloading` and `assembly_parsing`. A stage with
                nothing to report is absent, so an empty dict means a clean run.
        """

        existing_files = glob(r"*.json", root_dir=self.interface_dir)
        missing_pdbs = [x.lower() for x in pdb_list if f"{x.lower()}.json" not in existing_files]
        local_count = len(pdb_list) - len(missing_pdbs)
        logger.info(f"{local_count}/{len(pdb_list)} interfaces found locally", extra=self.log_extra)

        all_failures = {}
        if len(missing_pdbs) > 0:
            # The only directory written with `open` before anything has been downloaded into it; it also
            # makes `pisa_dir` for `error_path`. The download directories are made by their downloads.
            os.makedirs(self.interface_dir, exist_ok=True)

            success, sd_failures = self.download_missing_summaries(missing_pdbs)
            if sd_failures:
                all_failures["summary_downloading"] = sd_failures
            assembly_dict, sp_failures = self.parse_summaries(success)
            if sp_failures:
                all_failures["summary_parsing"] = sp_failures

            success, ad_failures = self.download_missing_assemblies(assembly_dict)
            if ad_failures:
                all_failures["assembly_downloading"] = ad_failures
            ap_failures = self.parse_assemblies(success)
            if ap_failures:
                all_failures["assembly_parsing"] = ap_failures

            if all_failures:
                with open(self.error_path, "w") as f:
                    json.dump(all_failures, f)
                logger.warning(
                    f"See {self.error_path} for details of {sum(len(v) for v in all_failures.values())} PISA failures",
                    extra=self.log_extra,
                )

        return all_failures

    def download_missing_summaries(self, pdb_codes):
        """
        Download the PDBe entry summary for each code into `summary_dir`, which names its assemblies.

        Codes already there are not refetched. The rest are requested `summary_batch_size` at a time.
        Each batch response is written to `_batch.json` in `summary_dir`, which creates the directory,
        split into one `<pdb_code>.json` per entry, and deleted. The API silently omits an entry it has
        no summary for, so a code absent from its batch's response counts as a failure, as does every
        code in a batch whose request failed.

        Args:
            pdb_codes (list): Lower-cased PDB codes to fetch.

        Returns:
            tuple: (list of codes whose summary is now on disk, list of codes that could not be
                fetched). A failure is logged and collected rather than raised.
        """
        logger.info(f"Downloading {len(pdb_codes)} PISA summaries", extra=self.log_extra)
        success = []
        failure = []
        to_fetch = []
        for pdb_code in dict.fromkeys(pdb_codes):
            if os.path.exists(os.path.join(self.summary_dir, f"{pdb_code}.json")):
                success.append(pdb_code)
            else:
                to_fetch.append(pdb_code)

        url = "https://www.ebi.ac.uk/pdbe/api/v2/pdb/entry/summary/"
        batch_fname = os.path.join(self.summary_dir, "_batch.json")
        for start in tqdm(range(0, len(to_fetch), self.summary_batch_size)):
            end = start + self.summary_batch_size
            batch = to_fetch[start:end]
            try:
                if not download_api(url, batch_fname, log_extra=self.log_extra, data=",".join(batch).encode()):
                    failure.extend(batch)
                    continue
                with open(batch_fname) as f:
                    data = json.load(f)
                for pdb_code in batch:
                    if pdb_code not in data:
                        failure.append(pdb_code)
                        continue
                    out_fname = os.path.join(self.summary_dir, f"{pdb_code}.json")
                    with open(f"{out_fname}.part", "w") as f:
                        json.dump({pdb_code: data[pdb_code]}, f)
                    os.replace(f"{out_fname}.part", out_fname)
                    success.append(pdb_code)
            finally:
                if os.path.exists(batch_fname):
                    os.remove(batch_fname)
        if len(failure) > 0:
            logger.warning(f"Failed to download {len(failure)} summaries", extra=self.log_extra)
        return success, failure

    def parse_summaries(self, pdb_codes):
        """
        Read cached summaries from `summary_dir` and collect each entry's assembly ids.

        Args:
            pdb_codes (list): Codes whose summaries are cached.

        Returns:
            tuple: (defaultdict of pdb_code -> list of assembly ids, list of codes whose summary
                could not be read). A summary naming anything other than exactly one entry counts as
                a failure; neither case raises.
        """
        logger.info(f"Parsing {len(pdb_codes)} PISA summaries", extra=self.log_extra)
        asm_dict = defaultdict(list)
        failure = []
        for pdb_code in tqdm(pdb_codes):
            try:
                fname = os.path.join(self.summary_dir, f"{pdb_code}.json")
                with open(fname) as f:
                    data = json.load(f)
                logger.debug(f"Summary data for {pdb_code}: {data}", extra=self.log_extra)
                if len(data[pdb_code]) != 1:
                    raise PocketMapperError(f"More than one entry in summary for {pdb_code}")
                for assembly in data[pdb_code][0]["assemblies"]:
                    asm_dict[pdb_code].append(assembly["assembly_id"])
            except Exception as e:
                logger.debug(f"Issue parsing summary for {pdb_code}: ({e})", extra=self.log_extra)
                failure.append(pdb_code)
        if len(failure) > 0:
            logger.warning(f"Failed to parse {len(failure)} PISA summaries", extra=self.log_extra)
        return asm_dict, failure

    def download_missing_assemblies(self, asm_dict):
        """
        Download the PISA interfaces for every assembly of every entry into `asm_dir`, from this source.

        Writes one file per assembly. An assembly already there counts as a success without costing a
        request, whichever source wrote it. FTP files are fetched `max_workers` at a time, and one the
        server does not have (404) is not retried; API requests are paced, one at a time.

        Args:
            asm_dict (dict): pdb_code -> list of assembly ids, as returned by `parse_summaries`.

        Returns:
            tuple: (defaultdict of pdb_code -> list of the assembly ids now on disk, list of
                `<pdb_code>_<assembly_id>` that could not be fetched), both in `asm_dict` order.
        """
        logger.info(
            f"Downloading {sum(len(v) for v in asm_dict.values())} PISA assemblies from {self.source}",
            extra=self.log_extra,
        )
        url_for, fetch = ASSEMBLY_SOURCES[self.source]
        jobs = {
            (pdb_code, asm): DownloadJob(url_for(pdb_code, asm), os.path.join(self.asm_dir, f"{pdb_code}_{asm}.json"))
            for pdb_code, assemblies in asm_dict.items()
            for asm in assemblies
        }
        fetched = download_many(jobs, fetch, self.max_workers, self.log_extra)

        failure = []
        success = defaultdict(list)
        for (pdb_code, asm), ok in fetched.items():
            if ok:
                success[pdb_code].append(asm)
            else:
                failure.append(f"{pdb_code}_{asm}")
        if len(failure) > 0:
            logger.warning(f"Failed to download {len(failure)} assemblies", extra=self.log_extra)
        return success, failure

    def parse_assemblies(self, asm_dict):
        """
        Flatten each entry's assemblies in `asm_dir` into one interface file keyed by chain pair.

        Writes `<pdb_code>.json` into `interface_dir`, which must already exist, for each entry with at
        least one usable interface: every assembly's interfaces merged into one dict keyed by the two
        chain ids sorted and joined by a comma (e.g. "A,B-2"), a later assembly winning on a repeated pair.

        Args:
            asm_dict (dict): pdb_code -> list of assembly ids.

        Returns:
            list: What did not make it into a file -- `<pdb_code>_<assembly_id>` for an assembly that
                could not be read, `<pdb_code>_<assembly_id>_parse_error` for one that raised, and
                `<pdb_code>_<assembly_id>_<interface_id>` for each interface skipped for not having
                exactly two molecules.
        """
        logger.info(
            f"Parsing {sum(len(v) for v in asm_dict.values())} PISA assemblies into per-interface files",
            extra=self.log_extra,
        )
        failure = []
        for pdb_code, assemblies in tqdm(asm_dict.items()):
            pdb_interfaces = {}
            for asm_id in assemblies:
                try:
                    # Opening assembly file
                    asm_fname = os.path.join(self.asm_dir, f"{pdb_code}_{asm_id}.json")
                    with open(asm_fname) as f:
                        data = json.load(f)

                    if len(data.keys()) != 1:
                        logger.debug(f"More than one entry in assembly for {pdb_code}_{asm_id}", extra=self.log_extra)
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
                            logger.debug(
                                f"More than one molecule in {pdb_code}, {interface['interface_id']}",
                                extra=self.log_extra,
                            )
                            failure.append(f"{pdb_code}_{asm_id}_{interface_id}")
                            continue
                        # Comma-separated, so chain ids of any length stay distinguishable
                        chain_ids = [molecule["chain_id"] for molecule in interface["molecules"]]
                        entry_name = ",".join(sorted(chain_ids))
                        pdb_interfaces[entry_name] = interface
                except Exception as e:
                    logger.debug(
                        f"Unexpected error parsing assembly for {pdb_code}_{asm_id}: ({e})", extra=self.log_extra
                    )
                    failure.append(f"{pdb_code}_{asm_id}_parse_error")

            if len(pdb_interfaces) > 0:
                interface_fname = os.path.join(self.interface_dir, f"{pdb_code}.json")
                with open(interface_fname, "w") as out_f:
                    json.dump(pdb_interfaces, out_f)
        if len(failure) > 0:
            logger.warning(f"Failed to parse {len(failure)} PISA assemblies", extra=self.log_extra)
        return failure
