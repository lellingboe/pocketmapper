"""
Download and flatten PDBe PISA interface data into a per-entry cache.

Three stages, each cached on disk so a rerun costs nothing: entry summaries, fetched in batches,
give the assembly ids, one file per assembly gives its interfaces, and those are flattened into a single
`<pdb_code>.json` per entry, `{"assemblies": [...], "interfaces": {...}}`. `"interfaces"` is keyed by sorted
chain pair and is what `pisa.load_interfaces` returns; `"assemblies"` lists the assemblies it was built from,
so an entry with no usable interface is cached too, and one built from only some of its assemblies is
completed by a later run.

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

With downloads off, the entry point uses only what `interface_pairs/` already holds, complete or not, and
sends no request.
"""

import json
import logging
import os
from collections import defaultdict
from glob import glob

from tqdm import tqdm

from pocketmapper.lib.constants import DEFAULT_PISA_SOURCE, PISA_SOURCES
from pocketmapper.lib.download import (
    DOWNLOAD_WORKERS,
    DownloadJob,
    check_json,
    download_api,
    download_file,
    download_many,
)
from pocketmapper.lib.exceptions import PocketMapperError
from pocketmapper.lib.files import part_path

logger = logging.getLogger(__name__)

# What a malformed or unexpected JSON file can raise while it is read and indexed
PARSE_ERRORS = (OSError, ValueError, KeyError, TypeError, IndexError, AttributeError)


def write_json(data, path):
    """
    Write `data` as JSON to `path` through a scratch `<path>.<hex>.part`, creating its directory.

    The file appears at `path` only once it is complete. The scratch name is unique per call, so two
    writers of one `path` never share it, and it is removed whether or not the write succeeds.

    Args:
        data (dict): JSON-serialisable data.
        path (str): Destination file.
    """
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    part_fpath = part_path(path)
    try:
        with open(part_fpath, "w") as f:
            json.dump(data, f)
        os.replace(part_fpath, path)
    finally:
        if os.path.exists(part_fpath):
            os.remove(part_fpath)


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
    separately usable, and each makes the directories it writes into. Every stage returns what it could not
    handle instead of raising, so one dead entry does not abort a large batch. Every file, downloaded or
    written here, is written through a uniquely named `.part` file, so an interrupted run cannot leave a
    truncated file for a later run to trust, and two runs sharing the cache cannot write each other's. An
    entry's interface file records the assemblies it was built from (`is_complete`), so an online run
    rebuilds one that lacks an assembly its summary names. With `download` False,
    `download_missing_interfaces` runs no stage, so it sends no request; the download stages called directly
    still fetch.
    """

    def __init__(
        self,
        pisa_dir,
        source=DEFAULT_PISA_SOURCE,
        summary_batch_size=50,
        max_workers=DOWNLOAD_WORKERS,
        download=True,
    ):
        """
        Lay out the cache under `pisa_dir` and configure where assembly interfaces come from.

        Args:
            pisa_dir (str): Cache directory. Holds `summaries/`, `assemblies/`, `interface_pairs/` (the
                per-entry files `pisa.load_interfaces` reads) and the failure report `errors.json`.
            source (str, optional): Where assembly interfaces are fetched from: "ftp" or "api". Summaries
                always come from the API. Defaults to DEFAULT_PISA_SOURCE.
            summary_batch_size (int, optional): Entries per summary request. The API fails with a 500
                somewhere above 800, so keep it well below that. Defaults to 50.
            max_workers (int, optional): Concurrent FTP downloads. Defaults to DOWNLOAD_WORKERS.
            download (bool, optional): Whether `download_missing_interfaces` may fetch and rebuild entries
                missing from `interface_pairs/` or incomplete there. Defaults to True.

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
        self.download = download
        self.summary_dir = os.path.join(pisa_dir, "summaries")
        self.asm_dir = os.path.join(pisa_dir, "assemblies")
        self.interface_dir = os.path.join(pisa_dir, "interface_pairs")
        self.error_path = os.path.join(pisa_dir, "errors.json")

    def download_missing_interfaces(self, pdb_list):
        """
        Populate the interface cache for a list of PDB entries.

        Entries complete in `interface_dir` (`is_complete`) cost no requests. For the rest, caches the
        summaries and assembly responses, writes one `<pdb_code>.json` per entry it could read a summary for
        into `interface_dir`, and writes `error_path` if anything failed. An entry rebuilt from cached
        summary and assemblies, such as one whose file is in the older flat shape, costs no request either.

        With `download` False, an entry counts as cached if its `<pdb_code>.json` exists, complete or not, and
        the rest skip every stage: one warning gives their count, and nothing is written. That includes an
        entry whose summaries and assemblies are cached but which has no `<pdb_code>.json`; it is not
        reparsed.

        Args:
            pdb_list (list): PDB codes, in any case.

        Returns:
            dict: Stage name -> the ids that stage could not handle, keyed `summary_downloading`,
                `summary_parsing`, `assembly_downloading` and `assembly_parsing`. A stage with
                nothing to report is absent, so an empty dict means a clean run. Always empty with
                `download` False: the entries skipped are counted only in the warning.
        """
        pdb_codes = list(dict.fromkeys(x.lower() for x in pdb_list))
        if self.download:
            missing_pdbs = [x for x in pdb_codes if not self.is_complete(x)]
            found = "complete"
        else:
            existing_files = set(glob("*.json", root_dir=self.interface_dir))
            missing_pdbs = [x for x in pdb_codes if f"{x}.json" not in existing_files]
            found = "found"
        local_count = len(pdb_codes) - len(missing_pdbs)
        logger.info(f"{local_count}/{len(pdb_codes)} PISA interface files {found} locally", extra=self.log_extra)

        all_failures = {}
        if len(missing_pdbs) > 0 and not self.download:
            logger.warning(
                f"Downloads are off: skipped {len(missing_pdbs)} PISA entries not cached locally",
                extra=self.log_extra,
            )
        elif len(missing_pdbs) > 0:
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
                write_json(all_failures, self.error_path)
                logger.warning(
                    f"See {self.error_path} for details of {sum(len(v) for v in all_failures.values())} PISA failures",
                    extra=self.log_extra,
                )

        return all_failures

    def is_complete(self, pdb_code):
        """
        Whether an entry's interface file is built from every assembly its cached summary names.

        False if `interface_dir` has no `<pdb_code>.json`, or one that is not in the shape `parse_assemblies`
        writes (the older flat dict of chain pairs included). True if the file's `"assemblies"` holds every id
        `summary_assemblies` gives, and also if no summary is cached or it cannot be read, so that checking
        never costs a request.

        Args:
            pdb_code (str): Lower-cased PDB code.

        Returns:
            bool: Whether the file can be used as is.
        """
        try:
            with open(os.path.join(self.interface_dir, f"{pdb_code}.json")) as f:
                data = json.load(f)
            used = set(data["assemblies"])
            if not isinstance(data["interfaces"], dict):
                return False
        except PARSE_ERRORS:
            return False
        try:
            listed = self.summary_assemblies(pdb_code)
        except PARSE_ERRORS:
            return True
        return used.issuperset(listed)

    def download_missing_summaries(self, pdb_codes):
        """
        Download the PDBe entry summary for each code into `summary_dir`, which names its assemblies.

        Codes already there are not refetched. The rest are requested `summary_batch_size` at a time.
        Each batch response is written to a scratch file in `summary_dir` named uniquely per batch
        (`_batch.json.<hex>.part`), which creates the directory, split into one `<pdb_code>.json` per
        entry, and deleted. The API silently omits an entry it has no summary for, so a code absent from
        its batch's response counts as a failure, as does every code in a batch whose request failed or
        whose response is not JSON.

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
        for start in tqdm(range(0, len(to_fetch), self.summary_batch_size)):
            end = start + self.summary_batch_size
            batch = to_fetch[start:end]
            batch_fname = part_path(os.path.join(self.summary_dir, "_batch.json"))
            try:
                if not download_api(
                    url, batch_fname, log_extra=self.log_extra, data=",".join(batch).encode(), check=check_json
                ):
                    failure.extend(batch)
                    continue
                with open(batch_fname) as f:
                    data = json.load(f)
                for pdb_code in batch:
                    if pdb_code not in data:
                        failure.append(pdb_code)
                        continue
                    write_json({pdb_code: data[pdb_code]}, os.path.join(self.summary_dir, f"{pdb_code}.json"))
                    success.append(pdb_code)
            finally:
                if os.path.exists(batch_fname):
                    os.remove(batch_fname)
        if len(failure) > 0:
            logger.warning(f"Failed to download {len(failure)} summaries", extra=self.log_extra)
        return success, failure

    def summary_assemblies(self, pdb_code):
        """
        Read an entry's cached summary from `summary_dir` and list its assembly ids.

        Args:
            pdb_code (str): Lower-cased PDB code.

        Returns:
            list: Assembly ids, in summary order; empty for an entry with no assemblies.

        Raises:
            OSError: If the summary is not cached or cannot be read.
            ValueError: If it is not JSON, or names anything other than exactly one entry.
            KeyError: If it lacks an expected key; TypeError, IndexError or AttributeError for another
                unexpected shape (all in PARSE_ERRORS).
        """
        with open(os.path.join(self.summary_dir, f"{pdb_code}.json")) as f:
            data = json.load(f)
        if len(data[pdb_code]) != 1:
            raise ValueError(f"expected one entry in summary for {pdb_code}, found {len(data[pdb_code])}")
        return [assembly["assembly_id"] for assembly in data[pdb_code][0]["assemblies"]]

    def parse_summaries(self, pdb_codes):
        """
        Read cached summaries from `summary_dir` and collect each entry's assembly ids.

        Args:
            pdb_codes (list): Codes whose summaries are cached.

        Returns:
            tuple: (defaultdict of pdb_code -> list of assembly ids, list of codes whose summary
                could not be read). Every code read is a key, with an empty list if its summary names no
                assembly. A summary naming anything other than exactly one entry counts as a failure;
                neither case raises.
        """
        logger.info(f"Parsing {len(pdb_codes)} PISA summaries", extra=self.log_extra)
        asm_dict = defaultdict(list)
        failure = []
        for pdb_code in tqdm(pdb_codes):
            try:
                asm_dict[pdb_code] = self.summary_assemblies(pdb_code)
            except PARSE_ERRORS as e:
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
        server does not have (404) is not retried; API requests are paced, one at a time. A response that
        is not JSON is not retried and not written, and counts as a failure.

        Args:
            asm_dict (dict): pdb_code -> list of assembly ids, as returned by `parse_summaries`.

        Returns:
            tuple: (defaultdict of pdb_code -> list of the assembly ids now on disk, list of
                `<pdb_code>_<assembly_id>` that could not be fetched), both in `asm_dict` order. Every
                entry of `asm_dict` is a key of the first, with an empty list if none of its assemblies is
                on disk.
        """
        logger.info(
            f"Downloading {sum(len(v) for v in asm_dict.values())} PISA assemblies from {self.source}",
            extra=self.log_extra,
        )
        url_for, fetch = ASSEMBLY_SOURCES[self.source]
        jobs = {
            (pdb_code, asm): DownloadJob(
                url_for(pdb_code, asm), os.path.join(self.asm_dir, f"{pdb_code}_{asm}.json"), check=check_json
            )
            for pdb_code, assemblies in asm_dict.items()
            for asm in assemblies
        }
        fetched = download_many(jobs, fetch, self.max_workers, self.log_extra)

        failure = []
        success = defaultdict(list, {pdb_code: [] for pdb_code in asm_dict})
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

        Writes `<pdb_code>.json` into `interface_dir` for every entry, replacing any file there, as
        `{"assemblies": [...], "interfaces": {...}}`. `"interfaces"` merges every used assembly's interfaces
        into one dict keyed by the two chain ids sorted and joined by a comma (e.g. "A,B-2"), a later
        assembly winning on a repeated pair; it is empty if no assembly has a usable interface.
        `"assemblies"` lists, in `asm_dict` order, the assemblies read and parsed without error, including
        one whose interfaces were all skipped. An assembly that could not be read or raised contributes no
        interface and is not listed, so the entry is not complete (`is_complete`).

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
            used = []
            for asm_id in assemblies:
                try:
                    asm_fname = os.path.join(self.asm_dir, f"{pdb_code}_{asm_id}.json")
                    with open(asm_fname) as f:
                        data = json.load(f)

                    if len(data.keys()) != 1:
                        logger.debug(
                            f"Expected one entry in assembly for {pdb_code}_{asm_id}, found {len(data.keys())}",
                            extra=self.log_extra,
                        )
                        failure.append(f"{pdb_code}_{asm_id}")
                        continue

                    if "PISA" in data:
                        data = data["PISA"]
                    else:
                        data = data[pdb_code]

                    # Merged only once the whole assembly has parsed, so an unlisted one contributes nothing
                    asm_interfaces = {}
                    for interface in data["assembly"]["interfaces"]:
                        interface_id = interface["interface_id"]
                        if len(interface["molecules"]) != 2:
                            logger.debug(
                                f"Expected two molecules in {pdb_code}_{asm_id} interface {interface_id}, "
                                f"found {len(interface['molecules'])}",
                                extra=self.log_extra,
                            )
                            failure.append(f"{pdb_code}_{asm_id}_{interface_id}")
                            continue
                        # Comma-separated, so chain ids of any length stay distinguishable
                        chain_ids = [molecule["chain_id"] for molecule in interface["molecules"]]
                        entry_name = ",".join(sorted(chain_ids))
                        asm_interfaces[entry_name] = interface
                    pdb_interfaces.update(asm_interfaces)
                    used.append(asm_id)
                except PARSE_ERRORS as e:
                    logger.debug(
                        f"Unexpected error parsing assembly for {pdb_code}_{asm_id}: ({e})", extra=self.log_extra
                    )
                    failure.append(f"{pdb_code}_{asm_id}_parse_error")

            entry = {"assemblies": used, "interfaces": pdb_interfaces}
            write_json(entry, os.path.join(self.interface_dir, f"{pdb_code}.json"))
        if len(failure) > 0:
            logger.warning(f"Failed to parse {len(failure)} PISA assemblies", extra=self.log_extra)
        return failure
