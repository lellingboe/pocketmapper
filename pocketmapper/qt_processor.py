"""
Parsing of query and target input strings into structured records.

Input grammar is `struct_info[:chain_info[:residue_info]]`, and either side may instead be a file
holding one such string per line -- README's "Input format" table documents the forms.
`determine_struct_type` and `determine_pocket_method` implement them, against the regexes defined
in `QTProcessor.__init__`. A caller may force a pocket method instead of the default "auto";
`validate_pocket_method` holds a forced one to the same patterns, so no record leaves here without
the chains and residues its method reads.

The original input string is kept verbatim as `pocket_id`, which is the identifier used throughout
the results. Orchestration lives in `pocketmapper.py`; this module only parses.
"""

# TODO Folder input - iterate through files in folder with correct format

import hashlib
import json
import logging
import os
import re
from dataclasses import asdict
from dataclasses import dataclass

import pandas as pd

from pocketmapper.constants import DEFAULT_CHAIN
from pocketmapper.constants import DEFAULT_POCKET_METHOD
from pocketmapper.exceptions import PocketMapperError
from pocketmapper.foldseek import bundled_foldseek_dbs
from pocketmapper.lib import split_chain_info


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

    def __init__(self, pdb_dir, alphafold_dir, foldseek_preprocessed_structure_dir, fsdb_dir):
        """
        Store the directories that record paths are resolved against, and compile the input regexes.

        Args:
            pdb_dir (str): Directory fetched PDB structures are written to; where `pdb` records get
                their `struct_path`.
            alphafold_dir (str): Directory fetched AlphaFold structures are written to; where
                `alphafold` records get their `struct_path`.
            foldseek_preprocessed_structure_dir (str): Directory the Foldseek preprocessing step
                writes to; where records get their `preprocess_path`.
            fsdb_dir (str): Directory holding downloaded Foldseek databases, used to locate the
                bundled `pdb` database.
        """
        # Held on the instance, not built per function: process_qt_cmdline_input names the side
        # being processed and the determine_* helpers it drives all log under that name.
        self.log_extra = {"stage": "Processing Inputs"}
        logging.debug("Started")

        self.pdb_dir = pdb_dir
        self.alphafold_dir = alphafold_dir
        self.foldseek_preprocessed_structure_dir = foldseek_preprocessed_structure_dir

        # Structure type regex patterns
        self.pdb_regex = r"^[a-zA-Z0-9]{4}$"
        self.uniprot_regex = r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$"  # https://www.uniprot.org/help/accession_numbers

        # Pocket method -> (the pocket-info pattern that spells it, what an entry must carry to use
        # it). Every pattern spells a chain as a single character, which is what lib.split_chain_info
        # relies on. One table for both directions: determine_pocket_method picks the first method
        # whose pattern the entry matches, and validate_pocket_method checks an entry against the
        # pattern of the method it was given, so an inferred method and a forced one cannot disagree.
        self.pocket_methods = {
            "whole_chain": (r"^[A-Za-z0-9]?\:?$", "a single chain and no residue list, e.g. '4Q5J:B'"),
            "pisa": (r"^[A-Za-z0-9]_[A-Za-z0-9]$", "a chain pair, e.g. '4Q5J:B_F'"),
            "passthrough": (r"^[A-Za-z0-9]\:(\d+\,?)+$", "a chain and a residue list, e.g. '4Q5J:A:10,11,12'"),
            # The partner chain is required, though vdw is tried last and a bare chain would be
            # claimed by whole_chain long before reaching it. That made the pattern's optional
            # partner unreachable when inferring, and wrong when validating: a forced vdw on "4Q5J:A"
            # matched, then reached gemmi as chain None. A trailing residue list still matches and is
            # still ignored -- the contacts are what define the pocket.
            "vdw": (r"^[A-Za-z0-9]_[A-Za-z0-9](\:(\d+\,?)*)?$", "a chain pair, e.g. '4Q5J:A_B'"),
        }

        # struct_type -> the pocket methods it supports, in the order they are tried. whole_chain is
        # first everywhere: a bare chain also matches the passthrough and vdw patterns, so it has to
        # win or an open entry would be read as an empty pocket. PISA is PDB-only, and an AlphaFold
        # model is a single chain, so neither reaches a method that needs a partner chain.
        self.struct_type_pocket_methods = {
            "alphafold": ("whole_chain", "passthrough"),
            "pdb": ("whole_chain", "pisa", "passthrough", "vdw"),
            "local_file": ("whole_chain", "passthrough", "vdw"),
        }

        self.bundled_foldseek_dbs = bundled_foldseek_dbs(fsdb_dir)

    def accepted_pocket_methods(self):
        """
        The pocket method values a caller may pass, in the order they are reported.

        Wider than the keys of `pocket_methods` by "auto", which infers the method from each entry,
        and by "foldseek_db", which names a whole database rather than a way of deriving a pocket from
        a structure. Neither has a pocket-info pattern of its own.

        Returns:
            tuple: The accepted `pocket_method` values.
        """
        return ("auto",) + tuple(self.pocket_methods) + ("foldseek_db",)

    def process_qt_cmdline_input(self, qt_input, name, pocket_method=DEFAULT_POCKET_METHOD):
        """
        Parse one side of the comparison -- a query or a target -- into a DataFrame of `QTRecord`s.

        Call it once per side; `name` only labels the side in log messages and errors.

        Args:
            qt_input (str): A query or target string ("struct_info:chain_info:residue_info"), or a
                path to a file holding one such string per line.
            name (str): Which side this input is, e.g. "query" or "target". Used in logging and
                error messages.
            pocket_method (str, optional): Pocket method to force for every entry, or "auto" to
                infer it from each input string. Defaults to DEFAULT_POCKET_METHOD.

        Raises:
            PocketMapperError: If `qt_input` is None, `pocket_method` is not one of
                `accepted_pocket_methods`, or the input file cannot be read.

        Returns:
            pandas.DataFrame: the parsed records for this side.
        """
        self.log_extra.update({"stage": f"Processing {name}"})
        logging.debug(f"Processing {name}", extra=self.log_extra)

        # Check that the input is specified
        if isinstance(qt_input, type(None)):
            logging.critical(f"{name} input is required. Exiting.", extra=self.log_extra)
            raise PocketMapperError(f"{name} input is required.")

        # The method applies to every entry, so an unrecognised one is a setting to correct
        # rather than an entry to skip -- raise once, before anything is parsed or fetched.
        if pocket_method not in self.accepted_pocket_methods():
            msg = (
                f"Unknown {name} pocket method {pocket_method!r}. "
                f"Choose one of: {', '.join(self.accepted_pocket_methods())}."
            )
            logging.critical(msg, extra=self.log_extra)
            raise PocketMapperError(msg)

        records = []
        if pocket_method != "foldseek_db" and os.path.isfile(
            qt_input
        ):  # if it's a file, process each line as a separate query/target
            try:
                with open(qt_input) as f:
                    for line in f.readlines():
                        records.append(self.parse_individual_qt(line.strip(), pocket_method=pocket_method))
            except Exception as e:
                logging.critical(f"Problem reading the file {qt_input}: {e}", extra=self.log_extra)
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
            pocket_method (str): Pocket method to force, or "auto" to infer it from the string.

        Returns:
            QTRecord: The parsed record, or None if the structure type or pocket method could not be
                determined -- both are logged as warnings so one bad line does not abort the batch.
        """
        # Foldseek databases have a special format and are treated differently
        if qt in self.bundled_foldseek_dbs or pocket_method == "foldseek_db":
            db = self.bundled_foldseek_dbs.get(qt)
            return QTRecord(
                pocket_id=qt,
                struct_info=qt,
                struct_type="foldseek_db",
                struct_path=db["db_path"] if db else qt,  # Bundled path if available, otherwise the input as is
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
            logging.warning(f"Could not determine structure type for {qt}", extra=self.log_extra)
            return None
        struct_path = self.determine_ref_struct_path(struct_info, struct_type)

        # Generate a unique name for the structure based on its components and a hash of the name
        input_fname = os.path.basename(struct_info).split(".")[0]
        domain_chain, _ = split_chain_info(chain_info)
        name = input_fname + "_" + domain_chain  # e.g., "P12345_A" or "1ABC_A"
        name_md5 = hashlib.md5(name.encode()).hexdigest()
        preprocess_name = name + "_" + name_md5
        preprocess_path = os.path.join(self.foldseek_preprocessed_structure_dir, preprocess_name + ".cif")
        preprocess_path_gz = preprocess_path + ".gz"

        resolved_pocket_method = (
            pocket_method if pocket_method != "auto" else self.determine_pocket_method(qt, struct_type)
        )
        if resolved_pocket_method is None:
            logging.warning(f"Could not determine pocket method for {qt}", extra=self.log_extra)
            return None

        # Run unconditionally rather than only for a forced method: for an inferred one it is a
        # tautology, since the method was chosen by the pattern it is checked against, and running it
        # either way makes "a record carries what its pocket method needs" hold for every record.
        if not self.validate_pocket_method(qt, resolved_pocket_method, struct_type):
            return None

        # The residue list IS the pocket on the passthrough path, so it is normalised here rather than
        # where the pocket is built -- before any structure is fetched. The pattern above has already
        # established that there is a list; what is left is the ids it holds.
        if resolved_pocket_method == "passthrough":
            residue_info = self.parse_residue_info(qt, residue_info)
            if residue_info is None:
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
        logging.debug(
            f"Processed {qt} into structured data: {json.dumps(asdict(record), indent=4)}", extra=self.log_extra
        )
        return record

    def parse_residue_info(self, qt, residue_info):
        """
        Normalise a passthrough entry's residue list.

        Rejects a list that cannot name residues at all -- absent, or holding anything but positive
        integers -- and collapses repeats. A repeat would otherwise reach `Pocket.res_auth_ids` twice
        and pair the two sides of a comparison off by one, with no error. Both checks still earn their
        place next to the passthrough pattern, which requires a list of digits: the absent case because
        this method is usable on its own, and the integer case because "0" is digits and is not a
        residue id.

        Args:
            qt (str): The whole input entry, named in the log messages.
            residue_info (str | None): The entry's `residue_info` portion.

        Returns:
            str: The comma-joined residue ids, repeats dropped and the typed order kept, or None if
                the list is unusable -- logged as a warning, so one bad entry does not abort a batch.
        """
        if not residue_info:
            logging.warning(
                f"No residue ids in {qt}, which the passthrough pocket method requires", extra=self.log_extra
            )
            return None

        res_ids = []
        duplicates = []
        for res_id in residue_info.split(","):
            # isdecimal rather than isdigit: int() accepts every decimal digit but not every digit,
            # so isdigit would let a superscript through to a ValueError further down.
            if not res_id.isdecimal() or int(res_id) < 1:
                logging.warning(
                    f"Residue id '{res_id}' in {qt} is not a positive integer; skipping this entry",
                    extra=self.log_extra,
                )
                return None
            res_id = str(int(res_id))  # Canonical, so "07" and "7" are recognised as the same residue
            if res_id in res_ids:
                if res_id not in duplicates:  # An id repeated three times is still one message
                    duplicates.append(res_id)
            else:
                res_ids.append(res_id)

        if duplicates:
            logging.warning(
                f"Residue id(s) {','.join(duplicates)} listed more than once in {qt}; using each one once",
                extra=self.log_extra,
            )
        return ",".join(res_ids)

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
            logging.critical(f"Directory input is not currently supported: {struct_str}", extra=self.log_extra)
            raise PocketMapperError(f"Directory input is not currently supported: {struct_str}")
        else:
            logging.warning(f"Could not determine structure type for {struct_str}", extra=self.log_extra)
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
                return os.path.join(self.alphafold_dir, f"{struct_info}.cif.gz")
            case "pdb":
                return os.path.join(self.pdb_dir, f"{struct_info}.cif.gz")
            case "local_file":
                return struct_info
            case _:
                logging.critical(
                    f"Unknown structure type {struct_type} for struct_info {struct_info}", extra=self.log_extra
                )
                raise PocketMapperError(f"Unknown structure type {struct_type} for struct_info {struct_info}")

    def pocket_info(self, qt_str):
        """
        The pocket portion of an input entry -- everything after the first ":".

        Args:
            qt_str (str): One input entry, "struct_info:chain_info:residue_info".

        Returns:
            str: The chain and residue fields as typed, or "" for a bare structure ("4Q5J"), which
                names no pocket at all.
        """
        return qt_str.split(":", 1)[1] if ":" in qt_str else ""

    def determine_pocket_method(self, qt_str, struct_type):
        """
        Determine the pocket method from the entry's pocket info and structure type.

        An entry that names no pocket -- a bare chain, or no chain at all -- is an open search:
        "whole_chain", meaning every CA-bearing residue of the chain is treated as the pocket.

        Which methods are reachable depends on the structure type, as `struct_type_pocket_methods`
        lays out: PISA is PDB-only, so a local file with a chain pair resolves to "vdw" instead, and
        an AlphaFold model, being a single chain, reaches neither "pisa" nor "vdw".

        Args:
            qt_str (str): The full input entry; everything after the first ":" is the pocket info.
            struct_type (str): As returned by `determine_struct_type`.

        Returns:
            str: One of "whole_chain", "pisa", "passthrough", "vdw", or None if no pattern matched.
        """
        pocket_info_str = self.pocket_info(qt_str)
        logging.debug(f"Determining pocket method for {pocket_info_str} using regex patterns", extra=self.log_extra)
        for method in self.struct_type_pocket_methods.get(struct_type, ()):
            if re.match(self.pocket_methods[method][0], pocket_info_str):
                return method
        return None

    def validate_pocket_method(self, qt_str, pocket_method, struct_type):
        """
        Check that an entry can supply what its pocket method needs.

        Two ways it cannot: the method is unavailable for this kind of structure -- PISA needs the
        PDB's interface data, and a chain pair means nothing on a single-chain AlphaFold model -- or
        the entry does not spell the chains and residues the method reads.

        Args:
            qt_str (str): The full input entry, named in the log messages.
            pocket_method (str): The method resolved for it, inferred or forced.
            struct_type (str): As returned by `determine_struct_type`.

        Returns:
            bool: True if the entry is usable, False if it is not -- logged as a warning, so one bad
                entry does not abort a batch.
        """
        supported = self.struct_type_pocket_methods.get(struct_type, ())
        if pocket_method not in supported:
            logging.warning(
                f"The {pocket_method} pocket method is not available for the {struct_type} entry {qt_str}; "
                f"{struct_type} entries support: {', '.join(supported)}. Skipping this entry",
                extra=self.log_extra,
            )
            return False

        pattern, needs = self.pocket_methods[pocket_method]
        pocket_info_str = self.pocket_info(qt_str)
        if not re.match(pattern, pocket_info_str):
            logging.warning(
                f"'{pocket_info_str}' in {qt_str} is not what the {pocket_method} pocket method reads; "
                f"it needs {needs}. Skipping this entry",
                extra=self.log_extra,
            )
            return False
        return True
