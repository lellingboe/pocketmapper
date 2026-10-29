"""
Finding, fetching and loading the reference structures PocketMapper works on.

`STRUCTURE_SOURCES` maps each downloadable `struct_type` to its URL builder and the transform its response
needs; a new source is one row. Structures land as gzipped mmCIF at each record's `struct_path`, which
doubles as the on-disk cache between runs. A local file is read where it is and never copied.

`StructureManager.path` works out where a structure or Foldseek database lives from its `struct_type` and
`struct_info`, and `record_path` does the same for a record. For now `fetch` still reads each record's own
`struct_path`, so it needs no directory; `path`, `record_path` and `get` need the directory of the type they
resolve. `load_structure` never downloads.
"""

import logging
import os

import gemmi

from pocketmapper.lib.download import DOWNLOAD_WORKERS, DownloadJob, download_many
from pocketmapper.lib.exceptions import PocketMapperError
from pocketmapper.lib.files import gzip_file
from pocketmapper.lib.foldseek import BUNDLED_FSDB_NAMES, bundled_foldseek_dbs

logger = logging.getLogger(__name__)


def pdb_url(pdb_code):
    """
    Build the wwPDB address of a PDB entry's divided mmCIF, which is served already gzipped.

    Args:
        pdb_code (str): The 4-character PDB code, in any case.

    Returns:
        str: The download URL.
    """
    code = pdb_code.lower()
    return f"https://files.wwpdb.org/pub/pdb/data/structures/divided/mmCIF/{code[1:3]}/{code}.cif.gz"


def alphafold_url(uniprot_acc):
    """
    Build the AlphaFold DB address of a UniProt accession's v6 model, which is served as plain mmCIF.

    Args:
        uniprot_acc (str): The UniProt accession.

    Returns:
        str: The download URL.
    """
    return f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_acc}-F1-model_v6.cif"


# struct_type -> (URL builder, transform). AlphaFold is gzipped on the way in to match the cache.
STRUCTURE_SOURCES = {"pdb": (pdb_url, None), "alphafold": (alphafold_url, gzip_file)}
# Types already on disk, or fetched elsewhere
NOT_DOWNLOADED = ("local_file", "foldseek_db")


class StructureManager:
    """
    The cache layout, downloads and loading of reference structures.
    """

    def __init__(
        self,
        pdb_dir=None,
        alphafold_dir=None,
        fsdb_dir=None,
        local_root=None,
        download=True,
        refresh=False,
        max_workers=DOWNLOAD_WORKERS,
    ):
        """
        Set up the cache directories and the download policy.

        Args:
            pdb_dir (str or None, optional): Cache directory of PDB structures, needed only to resolve a
                `pdb` structure (`path`, `record_path`, `get`). Defaults to None.
            alphafold_dir (str or None, optional): Cache directory of AlphaFold models, needed only to
                resolve an `alphafold` structure. Defaults to None.
            fsdb_dir (str or None, optional): Cache directory of downloaded Foldseek databases, needed only
                to resolve a bundled `foldseek_db` name. Defaults to None.
            local_root (str or None, optional): Directory that a relative local file or user Foldseek
                database path resolves against, read when `path` is called. Defaults to None, the current
                directory.
            download (bool, optional): False sends no request: a structure is available only if its file
                exists. Defaults to True.
            refresh (bool, optional): True downloads PDB and AlphaFold structures even when a cached file
                exists. Defaults to False.
            max_workers (int, optional): Downloads to run concurrently. Defaults to DOWNLOAD_WORKERS.

        Raises:
            PocketMapperError: If `refresh` is True and `download` is False.
        """
        if refresh and not download:
            msg = "Structures cannot be refreshed with downloads turned off"
            logger.critical(msg, extra={"stage": "Downloading Structures"})
            raise PocketMapperError(msg)
        self.dirs = {"pdb": pdb_dir, "alphafold": alphafold_dir}
        self.fsdb_dir = fsdb_dir
        self.local_root = local_root
        self.download = download
        self.refresh = refresh
        self.max_workers = max_workers

    def path(self, struct_type, struct_info):
        """
        Build the path a structure or Foldseek database is cached at, or read from.

        Args:
            struct_type (str): "pdb", "alphafold", "local_file" or "foldseek_db".
            struct_info (str): The PDB code or UniProt accession, in any case; the local file's path; or
                a bundled Foldseek database name or a user database path.

        Returns:
            str: `<dir>/<struct_info upper-cased>.cif.gz` for pdb and alphafold; for a local file or a user
                database, `struct_info` joined onto `local_root` (an absolute path comes back unchanged);
                for a bundled database name, its path from `lib.foldseek.bundled_foldseek_dbs`.

        Raises:
            PocketMapperError: If the type is unknown, or the directory it needs was not given.
        """
        log_extra = {"stage": "Locating Structures"}
        if struct_type == "local_file":
            return self.local_path(struct_info)
        if struct_type == "foldseek_db":
            if struct_info not in BUNDLED_FSDB_NAMES:
                return self.local_path(struct_info)
            if self.fsdb_dir is None:
                msg = f"The path of bundled Foldseek database {struct_info} needs fsdb_dir"
                logger.critical(msg, extra=log_extra)
                raise PocketMapperError(msg)
            return bundled_foldseek_dbs(self.fsdb_dir)[struct_info]
        if struct_type not in self.dirs:
            msg = f"Unknown structure type {struct_type} for struct_info {struct_info}"
            logger.critical(msg, extra=log_extra)
            raise PocketMapperError(msg)
        struct_dir = self.dirs[struct_type]
        if struct_dir is None:
            msg = f"The path of {struct_type} structure {struct_info} needs {struct_type}_dir"
            logger.critical(msg, extra=log_extra)
            raise PocketMapperError(msg)
        return os.path.join(struct_dir, f"{struct_info.upper()}.cif.gz")

    def local_path(self, path):
        """
        Resolve a user-given path against `local_root`, or the current directory if that is None.

        Args:
            path (str): The path as the user wrote it.

        Returns:
            str: `path` unchanged if absolute, otherwise joined onto the root.
        """
        root = os.getcwd() if self.local_root is None else self.local_root
        return os.path.join(root, path)

    def record_path(self, record):
        """
        Build the path of the structure or Foldseek database a record names.

        Args:
            record (dict): A record with `struct_type` and `struct_info`.

        Returns:
            str: The path, as `path` gives it.

        Raises:
            PocketMapperError: As `path`.
        """
        return self.path(record["struct_type"], record["struct_info"])

    def fetch(self, records):
        """
        Concurrently fetch the structures a list of records names.

        Writes each PDB and AlphaFold structure to its record's `struct_path`, unless a file is already
        there and `refresh` is off. With `download` off, nothing is fetched. Records sharing a
        `struct_info` give one download, from the first of them.

        Args:
            records (list of dict): Records, each with `struct_type`, `struct_info` and, for every type but
                `foldseek_db`, `struct_path`.

        Returns:
            dict: struct_info -> whether its structure is available. A local file is available if it
                exists, a Foldseek database always, and an unknown type never (with a warning).
        """
        log_extra = {"stage": "Downloading Structures"}
        found = {}
        jobs = {}
        for record in records:
            struct_type, struct_info = record["struct_type"], record["struct_info"]
            if struct_info in found or struct_info in jobs:
                continue
            if struct_type in STRUCTURE_SOURCES and self.download:
                url_for, transform = STRUCTURE_SOURCES[struct_type]
                jobs[struct_info] = DownloadJob(url_for(struct_info), record["struct_path"], transform)
            elif struct_type in STRUCTURE_SOURCES or struct_type == "local_file":
                found[struct_info] = os.path.exists(record["struct_path"])
            elif struct_type in NOT_DOWNLOADED:
                found[struct_info] = True
            else:
                logger.warning(f"Unknown structure type {struct_type} for struct_info {struct_info}", extra=log_extra)
                found[struct_info] = False
        found.update(download_many(jobs, max_workers=self.max_workers, log_extra=log_extra, overwrite=self.refresh))
        return found

    def get(self, struct_type, struct_info):
        """
        Fetch a structure if needed and load it.

        Args:
            struct_type (str): "pdb", "alphafold" or "local_file".
            struct_info (str): The PDB code or UniProt accession, in any case, or the local file's path.

        Returns:
            gemmi.Structure: The structure with its entity types set, or None (with a warning) if it is not
                available or cannot be read.

        Raises:
            PocketMapperError: If the type is unknown, or the cache directory it needs was not given.
        """
        struct_path = self.path(struct_type, struct_info)
        if struct_type != "local_file":
            struct_info = struct_info.upper()  # as parse writes it, so the URL builders get it upper-cased
        self.fetch([{"struct_type": struct_type, "struct_info": struct_info, "struct_path": struct_path}])
        return load_structure(struct_path, {"stage": "Loading Structures"})


def load_structure(struct, log_extra=None, entities=True):
    """
    Return a parsed structure with its entity types set, reading it first when given a path.

    Never downloads. `setup_entities()` is called on every structure returned, one passed in included,
    which modifies that structure in place. `lib.residues.ca_residues` needs the entity types: without
    them, a local mmCIF with no entity data gives empty pockets.

    Args:
        struct (gemmi.Structure or str): A parsed structure, or a path gemmi will read, its format
            inferred from the extension.
        log_extra (dict or None, optional): `extra` for the warnings. Defaults to None.
        entities (bool, optional): False skips `setup_entities()`, leaving the structure as read. Defaults
            to True.

    Returns:
        gemmi.Structure: The structure, or None (with a warning) if the file does not exist or cannot be
            read, or the structure has no model or its first model has no chain.
    """
    if isinstance(struct, gemmi.Structure):
        st = struct
    else:
        if not os.path.exists(struct):
            logger.warning(f"Structure file {struct} does not exist.", extra=log_extra)
            return None
        try:
            st = gemmi.read_structure(struct)
        except (RuntimeError, ValueError, IndexError, OSError) as e:
            logger.warning(f"Structure file {struct} cannot be read: {e}", extra=log_extra)
            return None
    if len(st) == 0:
        logger.warning(f"Structure {st.name} has no model.", extra=log_extra)
        return None
    if len(st[0]) == 0:
        logger.warning(f"The first model of structure {st.name} has no chain.", extra=log_extra)
        return None
    if entities:
        st.setup_entities()
    return st
