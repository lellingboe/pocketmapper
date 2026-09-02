"""
Parsing of query and target input strings into structured records.

Input grammar is `struct_info[:chain_info[:residue_info]]`, and either side may instead be a file
holding one such string per line -- README's "Input format" table documents the forms.
`determine_struct_type` and `determine_pocket_method` implement them, against the regexes defined
in `QTProcessor.__init__`.

The original input string is kept verbatim as `pocket_id`, which is the identifier used throughout
the results. Orchestration lives in `pocketmapper.py`; this module only parses.
"""

# TODO Folder input - iterate through files in folder with correct format

from dataclasses import asdict, dataclass
import hashlib
from importlib.resources import files
import logging
import os
import re
import pandas as pd
import json
from pocketmapper.constants import DEFAULT_CHAIN
from pocketmapper.exceptions import PocketMapperError

# The bundled human-domains Foldseek DB, named here and nowhere else -- bump it on a DB refresh.
BUNDLED_HUMAN_DOMAINS_DB = "human_v3_20260901"

# The UniProt coordinates of every entry in the DB above, shipped beside it and keyed by the same
# entry names. `pocketmapper._compare_pockets_based_on_alignment` resolves it and
# `pocket_comparison._synthesise_target_pocket` applies it, so the two files are now coupled: refresh
# the table with the DB, or a hit on an entry the table has lost aborts the run.
BUNDLED_HUMAN_DOMAINS_OFFSET_TABLE = "offset_table.tsv"


def _bundled_human_domains_path(filename):
    """
    Resolve a file shipped in the package's `human_domains` directory.

    Anchored on the `pocketmapper` package, not on `human_domains` itself: that directory holds only
    Foldseek DB files and has no __init__.py, so passing it to `files()` resolves a namespace package,
    which importlib.resources only learned to handle in 3.12. Going through the parent regular package
    works on every supported version.

    Args:
        filename (str): Name of the file within the `human_domains` directory.

    Returns:
        str: Absolute path to that file.
    """
    return str(files("pocketmapper").joinpath("human_domains", filename))


def bundled_human_domains_offset_table(struct_path):
    """
    The UniProt offset table shipped beside the bundled human-domains DB, or None for any other DB.

    Only the bundled DB ships a table, so a user-supplied Foldseek database keeps the 0-indexed
    positional target residue ids it has always had rather than being renumbered against a table that
    does not describe it. Matched on the resolved path rather than on the input string, so naming the
    bundled DB by its path is treated the same as naming it "human_domains".

    Args:
        struct_path (str): The Foldseek DB path in use, from the target record's `struct_path`.

    Returns:
        str | None: Path to the table, or None when `struct_path` is not the bundled DB.
    """
    if struct_path != _bundled_human_domains_path(BUNDLED_HUMAN_DOMAINS_DB):
        return None
    return _bundled_human_domains_path(BUNDLED_HUMAN_DOMAINS_OFFSET_TABLE)


@dataclass
class QTRecord:
    """
    A single parsed query/target entry.

    Holds the raw input string as `pocket_id` plus everything derived from it -- structure location,
    preprocessing paths and pocket method. `success` and `failure_reason` let a record survive a failed
    fetch so the reason can be reported alongside the ones that worked.
    """

    pocket_id: str
    struct_info: str | None = None
    chain_info: str | None = None
    residue_info: str | None = None
    struct_type: str | None = None
    struct_path: str | None = None
    preprocess_name: str | None = None
    preprocess_path: str | None = None
    preprocess_path_gz: str | None = None
    pocket_method: str | None = None
    success: bool = True
    failure_reason: str = ""


class QTProcessor:
    """
    Parses query and target input into `QTRecord` DataFrames.

    Handles both sides identically; `process_qt_cmdline_input` is the entry point and is called once
    per side.
    """

    def __init__(self, structure_dir, foldseek_preprocessed_structure_dir, fsdb_dir):
        """
        Store the directories that record paths are resolved against, and compile the input regexes.

        Args:
            structure_dir (str): Directory fetched reference structures are written to; where
                `pdb`/`alphafold` records get their `struct_path`.
            foldseek_preprocessed_structure_dir (str): Directory the Foldseek preprocessing step
                writes to; where records get their `preprocess_path`.
            fsdb_dir (str): Directory holding downloaded Foldseek databases, used to locate the
                bundled `pdb` database.
        """
        # Initialize logger
        self.logger = logging.getLogger(__name__)
        self._log_extra = {"stage": "QTProcessor Initialization"}
        self.logger.debug(
            "Started",
            extra=self._log_extra,
        )

        self.structure_dir = structure_dir
        self.foldseek_preprocessed_structure_dir = foldseek_preprocessed_structure_dir

        # TODO regexes for validating output when pocket_method is specified
        # Structure type regex patterns
        self.pdb_regex = r"^[a-zA-Z0-9]{4}$"
        self.uniprot_regex = r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$"  # https://www.uniprot.org/help/accession_numbers

        # Pocket method regex patterns.
        # whole_chain is checked before the others: a bare chain also matches the passthrough and vdw
        # patterns, so it has to win or an open entry would be read as an empty pocket.
        self.whole_chain_regex = r"^[A-Za-z0-9]?\:?$"  # pattern like "A", "A:", or nothing at all
        self.pisa_regex = r"^[A-Za-z0-9]_[A-Za-z0-9]$"  # pattern like "A_B"
        self.passthrough_regex = r"^[A-Za-z0-9]\:(\d+\,?)+$"  # pattern like "A:1,2,3"
        self.vdw_regex = r"^[A-Za-z0-9](_[A-Za-z0-9])?(\:(\d+\,?)*)?$"  # pattern like "A_B:1,2,3"

        self._bundled_foldseek_dbs = {
            "human_domains": _bundled_human_domains_path(BUNDLED_HUMAN_DOMAINS_DB),
            "pdb": os.path.join(fsdb_dir, "pdb"),
        }

    def process_qt_cmdline_input(self, qt_input, name, pocket_method=None):
        """
        Parse one side of the comparison -- a query or a target -- into a DataFrame of `QTRecord`s.

        Call it once per side; `name` only labels the side in log messages and errors.

        Args:
            qt_input (str): A query or target string ("struct_info:chain_info:residue_info"), or a
                path to a file holding one such string per line.
            name (str): Which side this input is, e.g. "query" or "target". Used in logging and
                error messages.
            pocket_method (str | None): Pocket method to force for every entry, or None to infer
                it from each input string.

        Returns:
            pandas.DataFrame: the parsed records for this side.
        """
        self._log_extra.update({"stage": f"Processing {name}"})
        self.logger.debug(f"Processing {name}", extra=self._log_extra)

        # Check that the input is specified
        if isinstance(qt_input, type(None)):
            self.logger.critical(f"{name} input is required. Exiting.", extra=self._log_extra)
            raise PocketMapperError(f"{name} input is required.")

        records = []
        if pocket_method != "foldseek_db" and os.path.isfile(
            qt_input
        ):  # if it's a file, process each line as a separate query/target
            try:
                with open(qt_input) as f:
                    for line in f.readlines():
                        records.append(self.parse_individual_qt(line.strip(), pocket_method=pocket_method))
            except Exception as e:
                self.logger.critical(f"Problem reading the file {qt_input}: {e}", extra=self._log_extra)
                raise PocketMapperError(f"Problem reading the file {qt_input}: {e}") from e
        else:
            records.append(self.parse_individual_qt(qt_input, pocket_method=pocket_method))

        records = [
            r for r in records if r is not None
        ]  # removing any None entries that may have been added due to errors
        return pd.DataFrame([asdict(r) for r in records])

    def parse_individual_qt(self, qt, pocket_method):
        """
        Parse one input string into a `QTRecord`.

        Also computes `preprocess_name` -- `<basename>_<chain><md5>` -- which is the key alignments are
        stored under, while pockets are keyed by `pocket_id`. One `preprocess_name` can serve several
        pocket_ids, since the same chain can carry more than one pocket.

        Args:
            qt (str): One input entry, "struct_info:chain_info:residue_info", or the name of a Foldseek
                database.
            pocket_method (str | None): Pocket method to force, or None to infer it from the string.

        Returns:
            QTRecord: The parsed record, or None if the structure type or pocket method could not be
                determined -- both are logged as warnings so one bad line does not abort the batch.
        """
        # Foldseek databases have a special format and are treated differently
        if qt in self._bundled_foldseek_dbs or pocket_method == "foldseek_db":
            return QTRecord(
                pocket_id=qt,
                struct_info=qt,
                struct_type="foldseek_db",
                struct_path=self._bundled_foldseek_dbs.get(
                    qt, qt
                ),  # Use the bundled path if available, otherwise use the input as is
            )

        # Unpack the input string into its components
        parts = qt.split(":")
        struct_info = parts[0] if len(parts) > 0 else None
        chain_info = parts[1] if len(parts) > 1 else None
        residue_info = parts[2] if len(parts) > 2 else None

        # An entry that names no chain is an open search over DEFAULT_CHAIN. A chain is still required
        # downstream -- preprocess_name bakes it in, and the pocket methods all index by it -- so fill
        # it in here rather than carrying a None through the pipeline.
        if not chain_info:
            chain_info = DEFAULT_CHAIN

        # determining structure info
        struct_type = self.determine_struct_type(struct_info)
        if struct_type is None:
            self.logger.warning(f"Could not determine structure type for {qt}", extra=self._log_extra)
            return None
        struct_path = self.determine_ref_struct_path(struct_info, struct_type)

        # Generate a unique name for the structure based on its components and a hash of the name
        input_fname = os.path.basename(struct_info).split(".")[0]
        name = input_fname + "_" + chain_info[0]  # e.g., "P12345_A" or "1ABC_A"
        name_md5 = hashlib.md5(name.encode()).hexdigest()
        preprocess_name = name + "_" + name_md5
        preprocess_path = os.path.join(self.foldseek_preprocessed_structure_dir, preprocess_name + ".cif")
        preprocess_path_gz = preprocess_path + ".gz"

        resolved_pocket_method = (
            pocket_method if pocket_method is not None else self.determine_pocket_method(qt, struct_type)
        )
        if resolved_pocket_method is None:
            self.logger.warning(f"Could not determine pocket method for {qt}", extra=self._log_extra)
            return None

        record = QTRecord(
            pocket_id=qt,
            struct_info=struct_info,
            chain_info=chain_info,
            residue_info=residue_info,
            struct_type=struct_type,
            struct_path=struct_path,
            preprocess_name=preprocess_name,
            preprocess_path=preprocess_path,
            preprocess_path_gz=preprocess_path_gz,
            pocket_method=resolved_pocket_method,
        )
        self.logger.debug(
            f"Processed {qt} into structured data: {json.dumps(asdict(record), indent=4)}", extra=self._log_extra
        )
        return record

    def determine_struct_type(self, struct_str):
        """
        Classify a structure identifier as "pdb", "alphafold" or "local_file".

        The regexes are tried before the filesystem check, so an identifier that also happens to name a
        file in the working directory is still read as an accession.

        Args:
            struct_str (str): The `struct_info` portion of an input entry.

        Returns:
            str: One of "pdb", "alphafold", "local_file", or None if nothing matched (logged as a
                warning).

        Raises:
            PocketMapperError: If `struct_str` names a directory, which is not supported.
        """
        if re.match(self.pdb_regex, struct_str):
            return "pdb"
        elif re.match(self.uniprot_regex, struct_str):
            return "alphafold"
        elif os.path.isfile(struct_str):
            return "local_file"
        elif os.path.isdir(struct_str):
            logging.critical(f"Directory input is not currently supported: {struct_str}", extra=self._log_extra)
            raise PocketMapperError(f"Directory input is not currently supported: {struct_str}")
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
                return os.path.join(self.structure_dir, f"{struct_info}.cif.gz")
            case "pdb":
                return os.path.join(self.structure_dir, f"{struct_info}.cif.gz")
            case "local_file":
                return struct_info
            case _:
                logging.critical(
                    f"Unknown structure type {struct_type} for struct_info {struct_info}", extra=self._log_extra
                )
                raise PocketMapperError(f"Unknown structure type {struct_type} for struct_info {struct_info}")

    def determine_pocket_method(self, qt_str, struct_type):
        """
        Determine the pocket method from the entry's pocket info and structure type.

        An entry that names no pocket -- a bare chain, or no chain at all -- is an open search:
        "whole_chain", meaning every CA-bearing residue of the chain is treated as the pocket.

        Which methods are reachable depends on the structure type: PISA is PDB-only, so a local file or
        AlphaFold model with a chain pair resolves to "vdw" instead.

        Args:
            qt_str (str): The full input entry; everything after the first ":" is the pocket info.
            struct_type (str): As returned by `determine_struct_type`.

        Returns:
            str: One of "whole_chain", "pisa", "passthrough", "vdw", or None if no pattern matched.
        """
        # An entry may be a bare structure ("4Q5J"), in which case there is no pocket info at all.
        pocket_info_str = qt_str.split(":", 1)[1] if ":" in qt_str else ""
        self.logger.debug(
            f"Determining pocket method for {pocket_info_str} using regex patterns", extra=self._log_extra
        )
        match struct_type:
            case "alphafold":
                if re.match(self.whole_chain_regex, pocket_info_str):
                    return "whole_chain"
                elif re.match(self.passthrough_regex, pocket_info_str):
                    return "passthrough"
                else:
                    return None
            case "pdb":
                if re.match(self.whole_chain_regex, pocket_info_str):
                    return "whole_chain"
                elif re.match(self.pisa_regex, pocket_info_str):
                    return "pisa"
                elif re.match(self.passthrough_regex, pocket_info_str):
                    return "passthrough"
                elif re.match(self.vdw_regex, pocket_info_str):
                    return "vdw"
                else:
                    return None
            case "local_file":
                if re.match(self.whole_chain_regex, pocket_info_str):
                    return "whole_chain"
                elif re.match(self.passthrough_regex, pocket_info_str):
                    return "passthrough"
                elif re.match(self.vdw_regex, pocket_info_str):
                    return "vdw"
                else:
                    return None
            case _:
                return None
