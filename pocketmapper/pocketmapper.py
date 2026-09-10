"""
PocketMapper: map and compare binding pockets across protein structures.

`search()` is the only public method; everything else on the class is an internal step of it. The
command line lives in `cli.py`, which is the only module that knows about argv or exit codes.

`search()` is the whole workflow, top to bottom:

1. `configure_workflow` -> Settings, directories, job_settings.json, logging.
2. `configure_query_target` -> QTProcessor -> one DataFrame of QTRecords per side.
3. `fetch_missing_structures` (or `fetch_missing_fsdb`) -> mmCIF into structure_dir.
4. `alignment` -> foldseek or the local aligner -> alignment.tsv.
5. `get_pockets` -> a `retrieve_*_pockets` builder per pocket method, merged into one
   pocket_id -> Pocket dict. The Pocket shape itself is declared in `pocket.py`.
6. `compare_pockets_based_on_alignment` -> pocket_comparison.compare_pockets -> pocket_comparison.tsv.
7. `align_structs` -> superposes the top align_count targets per query into aligned_structures/.

This is the only module that knows about `Settings`; components are handed the individual values
they need, so none of them has to build one to be usable on its own.

Author: Lachlan Ellingboe
"""

import json
import logging
import logging.config
import os
import shutil
from dataclasses import asdict
from dataclasses import dataclass
from dataclasses import field
from dataclasses import replace
from datetime import datetime

import pandas as pd

from pocketmapper.constants import ALIGN_STRUCT_METHODS
from pocketmapper.constants import FOLDSEEK_FORMAT_OUTPUT
from pocketmapper.constants import FOLDSEEK_INSTALL_HINT
from pocketmapper.downloads.pisa_downloader import PisaDownloader
from pocketmapper.downloads.structure_fetcher import StructureFetcher
from pocketmapper.exceptions import PocketMapperError
from pocketmapper.foldseek import bundled_foldseek_dbs
from pocketmapper.foldseek import check_foldseek
from pocketmapper.foldseek import run_foldseek
from pocketmapper.lib import is_within
from pocketmapper.lib import jsonify_dict
from pocketmapper.lib import parse_foldseek_pdb_entry_name
from pocketmapper.lib import safe_filename
from pocketmapper.lib import split_chain_info
from pocketmapper.pisa_parser import PisaParser
from pocketmapper.pocket_calculator import PocketCalculator
from pocketmapper.pocket_comparison import compare_pockets
from pocketmapper.pocket_comparison import parse_pocket_transform
from pocketmapper.pocket_parser import parse_pocket_from_struct
from pocketmapper.qt_processor import QTProcessor
from pocketmapper.sequence_aligner import SequenceAligner
from pocketmapper.structure_aligner import StructureAligner
from pocketmapper.structure_preprocessor import StructurePreprocessor


@dataclass
class Settings:
    """
    Fully resolved PocketMapper run configuration.

    Built by layering three sources in priority order (lowest to highest):
    dataclass defaults -> JSON settings file -> explicit CLI arguments.
    Derived paths (structure_dir, alignment_path, etc.) are filled in
    afterward by resolve_paths(), unless already set by the settings file.
    """

    query: str | None = None
    target: str | None = None
    cache_dir: str = "pocketmapper_cache"
    results_dir: str = field(default_factory=lambda: f"pocketmapper_results_{datetime.now().strftime('%y%m%d_%H%M%S')}")
    query_pocket_method: str | None = None
    target_pocket_method: str | None = None
    # Tri-state: None (the default) means "auto" -- use Foldseek when the binary is on PATH and
    # fall back to the local aligner when it is not. resolve_foldseek() turns this into a concrete
    # bool before anything else reads it, so the rest of the pipeline only ever sees True/False.
    foldseek: bool | None = None
    align_count: int = 10
    # Which transform superposes a target onto its query in step 7: "foldseek" (Foldseek's whole-chain
    # fit) or "pocket" (the fit of the two pockets on their overlapping residues). The default "auto"
    # is collapsed to one of those by resolve_align_struct_method(), so nothing downstream sees it.
    align_struct_method: str = "auto"
    verbosity: int = 3
    # query_dir, target_dir and foldseek_tmp_dir hold the per-run inputs actually handed to the
    # aligner and the pocket parser, and delete_tmp removes them on the way out. False keeps them:
    # a run that produced no rows is diagnosed from what it was given, which is gone by the time
    # anyone looks.
    delete_tmp: bool = True

    # Derived paths -- left unset (None) until resolve_paths() fills them in, unless explicitly
    # provided via the settings file or the matching command-line option.
    structure_dir: str | None = None
    pocket_dir: str | None = None
    foldseek_tmp_dir: str | None = None
    foldseek_preprocessed_structure_dir: str | None = None
    query_dir: str | None = None
    target_dir: str | None = None
    aligned_structure_dir: str | None = None
    alignment_path: str | None = None
    pocket_comparison_path: str | None = None
    job_settings_path: str | None = None
    log_path: str | None = None
    fsdb_dir: str | None = None

    def resolve_paths(self):
        """
        Return a copy of these settings with any unset derived paths filled in.

        Paths already set -- via the settings file or a command-line option -- are left untouched.
        Always call this on a Settings you built yourself; `search()` does it for you. Skipping it
        leaves the derived paths None and yields an opaque `TypeError: expected str, bytes or
        os.PathLike object, not NoneType` from inside `os.path.join`.

        Returns:
            Settings: A new instance with the derived paths resolved against cache_dir/results_dir.
        """
        derived = {
            "structure_dir": os.path.join(self.cache_dir, "ref_structures"),
            "pocket_dir": os.path.join(self.cache_dir, "pockets"),
            "foldseek_tmp_dir": os.path.join(self.cache_dir, "foldseek_tmp"),
            "foldseek_preprocessed_structure_dir": os.path.join(self.cache_dir, "foldseek_preprocessed_structures"),
            "query_dir": os.path.join(self.results_dir, "query_structures"),
            "target_dir": os.path.join(self.results_dir, "target_structures"),
            "aligned_structure_dir": os.path.join(self.results_dir, "aligned_structures"),
            "alignment_path": os.path.join(self.results_dir, "alignment.tsv"),
            "pocket_comparison_path": os.path.join(self.results_dir, "pocket_comparison.tsv"),
            "job_settings_path": os.path.join(self.results_dir, "job_settings.json"),
            "log_path": os.path.join(self.results_dir, "info.log"),
            "fsdb_dir": os.path.join(self.cache_dir, "fsdb"),
        }
        unset = {key: path_val for key, path_val in derived.items() if getattr(self, key) is None}
        return replace(self, **unset)


class PocketMapper:
    """
    The pipeline.

    `search()` is the only public method and runs the whole workflow; see the module docstring for
    its steps. Every other method is an internal step and carries a leading underscore.

    Usable as a library -- nothing here needs a terminal -- but note `search()` has global side
    effects: it reconfigures the *root* logger via `logging.config.dictConfig`, and deletes its
    temporary directories on the way out.
    """

    def __init__(self):
        """
        Initialise logging defaults and the Foldseek-database flags.

        Sets the root log format `"%(levelname)s: %(stage)s - %(msg)s"`, which is why every log call in
        the package must pass `extra={"stage": ...}` -- a record without it fails to format. Nothing
        enforces that.
        """
        self.log_extra = {"stage": "init"}
        self.log_fmt = "%(levelname)s: %(stage)s - %(msg)s"
        logging.getLogger(__name__)
        logging.basicConfig(level=logging.CRITICAL, format=self.log_fmt)

        self.fsdb_target = False
        self.fsdb_pdb_target = False
        # Set for real by resolve_foldseek(); read by configure_query_target to explain *why*
        # a Foldseek-DB target was rejected. True here so a caller that skips search() is not
        # told the binary is missing when nothing has looked for it.
        self.foldseek_available = True

    def configure_logging(self, verbosity, log_path):
        """
        Configure logging level and handlers based on user settings.

        Sets both a console stream handler and a file handler (to `log_path` in `settings`),
        adjusting the verbosity depending on the user's input.

        Args:
            verbosity (int): The verbosity level (4=DEBUG, 3=INFO, 2=WARNING, else ERROR).
            log_path (str): The path to the log file.
        """
        self.log_extra.update({"stage": "Configuring Logging"})

        # Set log level based on verbosity setting (default to INFO if not set)
        log_level = None
        if verbosity == 4:
            log_level = "DEBUG"
        elif verbosity == 3:
            log_level = "INFO"
        elif verbosity == 2:
            log_level = "WARNING"
        else:
            log_level = "ERROR"

        log_config = {
            "version": 1,
            "formatters": {
                "standard": {"format": self.log_fmt},
            },
            "handlers": {
                "console": {
                    "class": "logging.StreamHandler",
                    "level": log_level,
                    "formatter": "standard",
                    "stream": "ext://sys.stdout",
                },
                "file": {
                    "class": "logging.FileHandler",
                    "level": log_level,
                    "formatter": "standard",
                    "filename": log_path,
                },
            },
            "root": {
                "handlers": ["console", "file"],
                "level": log_level,
                "propagate": True,
            },
        }
        logging.config.dictConfig(log_config)

    def search(
        self,
        query=None,  # settings passed to configure
        target=None,
        settings=None,
        cache_dir=None,
        results_dir=None,
        verbosity=None,
        foldseek=None,
        align_count=None,
        align_struct_method=None,
        query_pocket_method=None,
        target_pocket_method=None,
        delete_tmp=None,
        structure_dir=None,
        pocket_dir=None,
        foldseek_tmp_dir=None,
        foldseek_preprocessed_structure_dir=None,
        query_dir=None,
        target_dir=None,
        aligned_structure_dir=None,
        alignment_path=None,
        pocket_comparison_path=None,
        job_settings_path=None,
        log_path=None,
        fsdb_dir=None,
    ):
        """
        Orchestrate and run the full PocketMapper search workflow.

        Args:
            query (str): Target query identifier, string or path to a list.
            target (str): Target structure identifier, string or path to a list.
            settings (str, optional): Path to a JSON settings file.
            cache_dir (str, optional): Directory to cache intermediate structures.
            results_dir (str, optional): Directory to output results to.
            verbosity (int, optional): Control logging level.
            foldseek (bool, optional): Use foldseek for structure alignment instead of local sequence
                alignment. Left unset, foldseek is used when the binary is on PATH and the local
                aligner is used with a warning when it is not. True makes foldseek a hard
                requirement -- a missing binary is an error; False always uses the local aligner.
            align_count (int, optional): Number of top targets to superpose onto each query.
            align_struct_method (str, optional): Which transform superposes a target onto its query --
                'foldseek' for Foldseek's whole-chain fit, 'pocket' for the fit of the two pockets
                on their overlapping residues, or 'auto' (the default) for 'foldseek' when Foldseek is
                in use and 'pocket' with the local aligner, which produces no chain transform at all.
            query_pocket_method (str, optional): Force a pocket method for every query entry --
                'pisa', 'passthrough', 'vdw', 'whole_chain' or 'foldseek_db'. Left unset, it is
                inferred per entry from the input string.
            target_pocket_method (str, optional): As `query_pocket_method`, for the target side.
            delete_tmp (bool, optional): Delete query_dir, target_dir and -- when Foldseek did the
                aligning -- foldseek_tmp_dir at the end of the run. Defaults to True; False keeps
                them for inspection.
            structure_dir (str, optional): Cache of fetched reference structures.
                Defaults to <cache_dir>/ref_structures.
            pocket_dir (str, optional): Cache of parsed pockets. Defaults to <cache_dir>/pockets.
            foldseek_tmp_dir (str, optional): Foldseek's scratch directory, deleted after a Foldseek
                run. Defaults to <cache_dir>/foldseek_tmp.
            foldseek_preprocessed_structure_dir (str, optional): Cache of the single-chain structures
                Foldseek is given. Defaults to <cache_dir>/foldseek_preprocessed_structures.
            query_dir (str, optional): Per-run query structures, deleted at the end of the run.
                Defaults to <results_dir>/query_structures.
            target_dir (str, optional): Per-run target structures, deleted at the end of the run.
                Defaults to <results_dir>/target_structures.
            aligned_structure_dir (str, optional): Where the superposed structures for the top hits
                are written. Defaults to <results_dir>/aligned_structures.
            alignment_path (str, optional): Where the alignment table is written.
                Defaults to <results_dir>/alignment.tsv.
            pocket_comparison_path (str, optional): Where the pocket comparison table is written.
                Defaults to <results_dir>/pocket_comparison.tsv.
            job_settings_path (str, optional): Where this run's resolved settings are dumped.
                Defaults to <results_dir>/job_settings.json.
            log_path (str, optional): Where the run log is written.
                Defaults to <results_dir>/info.log.
            fsdb_dir (str, optional): Cache of bundled Foldseek databases.
                Defaults to <cache_dir>/fsdb.

        Returns:
            None: Results are written to `results_dir` -- read pocket_comparison.tsv and
                alignment.tsv from there.
        """
        self.log_extra = {
            "stage": "Starting Search"
        }  # dict needed for logging extra info, can be updated throughout the process to indicate the current stage in logs

        # The Settings fields this call is overriding. `settings` is deliberately absent: it names the
        # JSON file the overrides sit on top of, and Settings has no such field. Keeping the dict here,
        # next to the signature it mirrors, is what keeps a new option from being added to one and not
        # the other.
        cli_overrides = {
            "query": query,
            "target": target,
            "cache_dir": cache_dir,
            "results_dir": results_dir,
            "foldseek": foldseek,
            "verbosity": verbosity,
            "align_count": align_count,
            "align_struct_method": align_struct_method,
            "query_pocket_method": query_pocket_method,
            "target_pocket_method": target_pocket_method,
            "delete_tmp": delete_tmp,
            "structure_dir": structure_dir,
            "pocket_dir": pocket_dir,
            "foldseek_tmp_dir": foldseek_tmp_dir,
            "foldseek_preprocessed_structure_dir": foldseek_preprocessed_structure_dir,
            "query_dir": query_dir,
            "target_dir": target_dir,
            "aligned_structure_dir": aligned_structure_dir,
            "alignment_path": alignment_path,
            "pocket_comparison_path": pocket_comparison_path,
            "job_settings_path": job_settings_path,
            "log_path": log_path,
            "fsdb_dir": fsdb_dir,
        }

        self.settings = self.configure_workflow(settings, cli_overrides)
        self.query_df, self.target_df = (
            self.configure_query_target()
        )  # parses the query and target inputs to determine their types and sets up the relevant data structures for each entry

        self.query_df = self.fetch_missing_structures(
            "query", self.query_df, self.settings.structure_dir
        )  # Fetch any missing structures
        if self.fsdb_target:
            self.fetch_missing_fsdb(
                self.target_df, self.settings.foldseek_tmp_dir
            )  # Fetch any missing foldseek databases
        else:
            self.target_df = self.fetch_missing_structures(
                "target", self.target_df, self.settings.structure_dir
            )  # Fetch any missing structures

        self.alignment()  # Align the query and target structures using either local sequence alignment or foldseek based on the settings
        pockets = self.get_pockets()  # Adds seq_pos and ca-coords to the pocket info dict
        self.compare_pockets_based_on_alignment(pockets)
        self.align_structs()
        self.delete_tmp()

        logging.info("PocketMapper search completed successfully.", extra={"stage": "End"})

    def configure_workflow(self, settings_file, cli_overrides):
        """
        Build the fully resolved `Settings` for this run.

        Layers three sources in priority order -- dataclass defaults, then an optional JSON settings file,
        then the arguments passed to `search()` -- then resolves the derived paths, creates the
        directories and writes job_settings.json.

        Args:
            settings_file (str or None): Path to a JSON settings file, or None for none.
            cli_overrides (dict): Settings field name -> value from `search()`. A None value means
                "not supplied" and is dropped, which is what leaves the settings file in charge of
                that field; anything else wins over the file.

        Returns:
            Settings: The resolved configuration. Also written to `job_settings_path`.
        """
        self.log_extra.update({"stage": "Configuring Settings"})

        # 1. Base defaults
        settings = Settings()

        # 2. Populate settings from the settings file if provided
        if settings_file is not None:
            if not os.path.isfile(settings_file):
                logging.critical(f"Settings file not found: {settings_file}", extra=self.log_extra)
                raise PocketMapperError(f"Settings file not found: {settings_file}")
            try:
                with open(settings_file) as f:
                    settings_data_from_file = json.load(f)
                settings = replace(settings, **settings_data_from_file)
            except TypeError as e:
                logging.critical(f"Unknown setting(s) in {settings_file}: {e}", extra=self.log_extra)
                raise PocketMapperError(f"Unknown setting(s) in {settings_file}: {e}") from e
            except Exception as e:
                logging.critical(
                    f"Error reading settings file: {settings_file}. Is it in JSON format?", extra=self.log_extra
                )
                raise PocketMapperError(f"Error reading settings file: {settings_file}. Is it in JSON format?") from e

        # 3. Override settings with the arguments explicitly passed to search()
        supplied = {key: value for key, value in cli_overrides.items() if value is not None}
        settings = replace(settings, **supplied)

        # 4. Computed paths (only fills in paths not already set by the settings file or an argument)
        settings = settings.resolve_paths()

        # Ensure all necessary directories exist before proceeding, creating them if needed
        dirs_to_create = [
            "structure_dir",
            "query_dir",
            "target_dir",
            "pocket_dir",
            "foldseek_preprocessed_structure_dir",
            "aligned_structure_dir",
        ]
        for dir_key in dirs_to_create:
            path = getattr(settings, dir_key)
            try:
                os.makedirs(path, exist_ok=True)
            except OSError as e:
                logging.critical(f"Error creating directory {path}", extra=self.log_extra)
                raise PocketMapperError(f"Error creating directory {path}") from e

        self.configure_logging(settings.verbosity, settings.log_path)

        # 4b. Resolve the tri-state foldseek setting into a concrete bool. Must come after
        # configure_logging (the root logger is still at CRITICAL before it, so the fallback
        # warning would be swallowed) and before the settings are logged and dumped below, so
        # job_settings.json records what the run actually did.
        settings = self.resolve_foldseek(settings)

        # 4c. Same reasoning, and it reads the bool resolve_foldseek just settled, so it must follow it.
        settings = self.resolve_align_struct_method(settings)

        logging.info(f"Settings: {json.dumps(asdict(settings), indent=4)}", extra=self.log_extra)

        # 5. Output dump
        try:
            os.makedirs(os.path.dirname(settings.job_settings_path), exist_ok=True)
            with open(settings.job_settings_path, "w") as f:
                json.dump(asdict(settings), f, indent=4)
            logging.info(f"Settings successfully dumped to {settings.job_settings_path}", extra=self.log_extra)
        except Exception as e:
            logging.error(f"Failed to dump settings to {settings.job_settings_path}: {e}", extra=self.log_extra)
        return settings

    def resolve_foldseek(self, settings):
        """
        Turn the tri-state `foldseek` setting into a concrete bool.

        Foldseek is an optional external binary, so the default (None, "auto") is resolved against
        what is actually installed: foldseek when it is on PATH, the local BLOSUM62 aligner with a
        warning when it is not. An explicit True is a hard requirement and errors instead of falling
        back; an explicit False always means the local aligner and never probes for the binary.

        Called before any structure is fetched, so an unmet requirement fails without wasted
        downloads rather than partway through the run at the first `run_foldseek` call.

        Args:
            settings (Settings): Settings whose `foldseek` field may still be None.

        Returns:
            Settings: A copy with `foldseek` set to True or False.

        Raises:
            PocketMapperError: If foldseek was explicitly requested but is not installed.
        """
        # Local stage dict: configure_logging leaves self.log_extra reading "Configuring Logging".
        stage = {"stage": "Configuring Settings"}

        if settings.foldseek is False:
            return settings

        self.foldseek_available = check_foldseek()
        if self.foldseek_available:
            return replace(settings, foldseek=True)

        if settings.foldseek is True:
            msg = (
                "Foldseek alignment was requested but 'foldseek' was not found on PATH. "
                f"{FOLDSEEK_INSTALL_HINT} Alternatively, set --foldseek False to use the local "
                "BLOSUM62 sequence aligner."
            )
            logging.critical(msg, extra=stage)
            raise PocketMapperError(msg)

        logging.warning(
            "'foldseek' not found on PATH; falling back to the local BLOSUM62 sequence aligner. "
            "The local aligner produces no whole-chain transform, so aligned_structures/*.pdb are "
            f"superposed on the pocket instead (see --align_struct_method). {FOLDSEEK_INSTALL_HINT}",
            extra=stage,
        )
        return replace(settings, foldseek=False)

    def resolve_align_struct_method(self, settings):
        """
        Turn the tri-value `align_struct_method` setting into "pocket" or "foldseek".

        "foldseek" uses Foldseek's whole-chain transform from alignment.tsv; "pocket" uses the
        superposition of the two pockets on their overlapping residues, which `compare_pockets`
        already writes to pocket_comparison.tsv. The default "auto" picks whichever the run can
        actually do: the local BLOSUM62 aligner writes "-" for the chain transform, so it has only
        the pocket one.

        An explicit "foldseek" without the binary is an error rather than a silent switch to "pocket" --
        the same call as an unmet `--foldseek True`, and for the same reason: better to fail before
        anything is downloaded than to hand back a method the user did not ask for.

        Called from configure_workflow after resolve_foldseek, whose resolved bool it reads, and
        before the settings are logged and dumped, so job_settings.json records what the run did.

        Args:
            settings (Settings): Settings whose `foldseek` is already a concrete bool.

        Returns:
            Settings: A copy with `align_struct_method` set to "pocket" or "foldseek".

        Raises:
            PocketMapperError: If the value is not one of ALIGN_STRUCT_METHODS, or "foldseek" was
                asked for on the local-aligner path.
        """
        stage = {"stage": "Configuring Settings"}

        method = settings.align_struct_method
        # The CLI always hands over a str, but a settings file can hold anything at all.
        method = method.lower() if isinstance(method, str) else method
        if method not in ALIGN_STRUCT_METHODS:
            msg = (
                f"Unknown align_struct_method {settings.align_struct_method!r}. "
                f"Choose one of: {', '.join(ALIGN_STRUCT_METHODS)}."
            )
            logging.critical(msg, extra=stage)
            raise PocketMapperError(msg)

        if method == "auto":
            method = "foldseek" if settings.foldseek else "pocket"
            logging.info(
                f"align_struct_method 'auto' resolved to '{method}' "
                f"({'foldseek' if settings.foldseek else 'the local aligner'} is in use)",
                extra=stage,
            )
        elif method == "foldseek" and not settings.foldseek:
            msg = (
                "align_struct_method 'foldseek' needs Foldseek's whole-chain transform, but this run "
                "uses the local BLOSUM62 aligner, which does not produce one. Use "
                "--align_struct_method pocket, or enable foldseek."
            )
            logging.critical(msg, extra=stage)
            raise PocketMapperError(msg)

        return replace(settings, align_struct_method=method)

    def configure_query_target(self):
        """
        Parse the query and target inputs into record DataFrames.

        Uses `QTProcessor` to parse, validate and resolve each side's structure type and pocket method.
        Also rejects `align_struct_method="pocket"` against a Foldseek-DB target here, before anything is
        fetched -- see the check itself for why one early rejection covers both kinds of database.

        Returns:
            tuple: (query_df, target_df). Note these are returned, not stored on the instance; the caller
                assigns them.

        Raises:
            PocketMapperError: If either side is unusable, or on the align_struct_method rejection above.
        """
        self.log_extra.update({"stage": "Determine Query/Target Types"})

        qtprocessor = QTProcessor(
            structure_dir=self.settings.structure_dir,
            foldseek_preprocessed_structure_dir=self.settings.foldseek_preprocessed_structure_dir,
            fsdb_dir=self.settings.fsdb_dir,
        )
        q_df = qtprocessor.process_qt_cmdline_input(
            qt_input=self.settings.query,
            name="query",
            pocket_method=self.settings.query_pocket_method,
        )
        t_df = qtprocessor.process_qt_cmdline_input(
            qt_input=self.settings.target,
            name="target",
            pocket_method=self.settings.target_pocket_method,
        )

        errors = []
        if len(q_df) < 1:
            logging.critical("No valid query entries after processing", extra=self.log_extra)
            errors.append("no valid query entries")
        if len(t_df) < 1:
            logging.critical("No valid target entries after processing", extra=self.log_extra)
            errors.append("no valid target entries")
        if errors:
            raise PocketMapperError("; ".join(errors))

        if t_df.loc[0, "struct_type"] == "foldseek_db":
            if self.settings.foldseek:
                self.fsdb_target = True
                # Neither kind of Foldseek DB can be superposed on its pocket, and which kind this is
                # is not known until expand_fsdb_pdb_targets has read the hit names -- so reject both
                # here, before anything is fetched. Unreachable from "auto": a DB target forces
                # foldseek on, and auto resolves to "foldseek" whenever it is on.
                if self.settings.align_struct_method == "pocket":
                    msg = (
                        "align_struct_method 'pocket' is not available against a Foldseek database "
                        "target. A human_domains-style hit has no coordinates to superpose at all, and "
                        "a PDB database's structures are assemblies while its pockets come from the "
                        "wwPDB asymmetric unit, so a pocket fit would be applied in the wrong frame. "
                        "Use --align_struct_method foldseek."
                    )
                    logging.critical(msg, extra=self.log_extra)
                    raise PocketMapperError(msg)
            else:
                # foldseek is already resolved to a concrete bool here, so False means either the
                # binary is missing or the user turned it off -- say which, since the fixes differ.
                if not self.foldseek_available:
                    msg = (
                        "A Foldseek database was specified as the target, which requires the "
                        f"'foldseek' binary, but it was not found on PATH. {FOLDSEEK_INSTALL_HINT}"
                    )
                else:
                    msg = (
                        "Foldseek database specified as target but foldseek is not enabled. "
                        "Remove --foldseek False to use it."
                    )
                logging.critical(msg, extra=self.log_extra)
                raise PocketMapperError(msg)
        return q_df, t_df

    def fetch_missing_structures(self, name, qt_df, out_dir):
        """
        Download the reference structures a side's records need.

        Uses `StructureFetcher`; records whose fetch fails are flagged rather than dropped, so the reason
        survives into the results. A side with nothing left is an error.

        Args:
            name (str): Which side this is, e.g. "query" or "target". Used in logging.
            qt_df (pandas.DataFrame): That side's records.
            out_dir (str): Directory to fetch into.

        Returns:
            pandas.DataFrame: The records with `success` and `failure_reason` updated.

        Raises:
            PocketMapperError: If no structure for this side could be fetched.
        """
        # Downloading structures
        self.log_extra.update({"stage": "Fetching Missing Structures"})
        logging.info("Starting", extra=self.log_extra)
        structure_fetcher = StructureFetcher()

        logging.debug(f"{name.capitalize()} data before fetching structures: \n{qt_df.head()}", extra=self.log_extra)

        # Get list of unique structures to fetch based on struct_info and struct_type
        unique_records = qt_df.drop_duplicates(subset="struct_info").to_dict(orient="records")

        # Update structure fetcher and fetch structures
        structure_fetcher.set_output_directory(out_dir)
        structure_fetcher.update_cache()
        results = structure_fetcher.fetch_structures(unique_records)
        logging.debug(f"Structure fetcher results: {results}", extra=self.log_extra)

        # Update the dataframe with success/failure information
        qt_df["success"] = qt_df["struct_info"].map(results).fillna(False)
        qt_df.loc[~qt_df["success"], "failure_reason"] = "structure_not_found"

        # Logging results of structure fetching and updating query and target data with success/failure info
        logging.info(
            f"{sum(results.values())}/{len(results)} {name} required structures available",
            extra=self.log_extra,
        )
        if len(qt_df.query("success == False")) > 0:
            logging.warning(
                f"Missing structures for {name}(s): {', '.join(qt_df.loc[~qt_df['success'], 'pocket_id'].unique().tolist())}",
                extra=self.log_extra,
            )

        # Verifying sufficient structures were found to continue
        if qt_df["success"].sum() < 1:
            logging.critical(f"Insufficient {name} structures after fetching", extra=self.log_extra)
            raise PocketMapperError(f"Insufficient {name} structures after fetching. No valid {name} entries remain.")
        return qt_df

    def fetch_missing_fsdb(self, qt_df, tmp_dir):
        """
        Download a bundled Foldseek database if it is not already on disk.

        Args:
            qt_df (pandas.DataFrame): The target records; the database is named by row 0.
            tmp_dir (str): Scratch directory for `foldseek databases`.

        Returns:
            None: The database is written to the record's `struct_path`.

        Raises:
            PocketMapperError: If the destination cannot be created, or -- from `run_foldseek` --
                if the download fails.
        """
        self.log_extra.update({"stage": "Fetching Missing Foldseek Database"})
        fsdb_name = qt_df.loc[0, "struct_info"].upper()
        fsdb_path = qt_df.loc[0, "struct_path"]
        if not os.path.exists(fsdb_path):
            logging.info(f"Fetching bundled Foldseek database '{fsdb_name}' to {fsdb_path}", extra=self.log_extra)
            try:
                os.makedirs(os.path.dirname(fsdb_path), exist_ok=True)
            except Exception as e:
                msg = f"Failed to create the directory for Foldseek database '{fsdb_name}' at {fsdb_path}: {e}"
                logging.critical(msg, extra=self.log_extra)
                raise PocketMapperError(msg) from e
            run_foldseek(["databases", fsdb_name, fsdb_path, tmp_dir], self.log_extra)
            logging.info(f"Successfully fetched Foldseek database '{fsdb_name}'", extra=self.log_extra)

    def alignment(self):
        """
        Coordinate structural alignment routes bridging query and target proteins.

        Dispatches to `foldseek_alignment()` if the foldseek flag is set, else rolls
        back to `local_alignment()`. Requisites like `foldseek_preprocessing()`
        precede foldseek routines.

        Returns:
            None
        """
        log_extra = {"stage": "Alignment"}
        if self.settings.foldseek:
            logging.info("Preprocessing structures for Foldseek...", extra=log_extra)
            self.foldseek_preprocessing()
            logging.info("Running Foldseek easy-search...", extra=log_extra)
            self.foldseek_alignment()
        else:
            logging.info("Running local pairwise aligner...", extra=log_extra)
            self.local_alignment()

    def foldseek_preprocessing(self):
        """
        Split complex PDB structures into single chain mmCIF files for Foldseek processing.
        Copies structure to the query/target directories, ensuring that the Foldseek preprocessed structure directory is populated.

        Instantiates a `StructurePreprocessor` mapping items sourced from the
        local datastores towards isolated files targeting specific sequences.

        Returns:
            None
        """
        stage = {"stage": "Preprocessing Structures"}

        structure_preprocessor = StructurePreprocessor()
        qtdf_dir_iter = [(self.query_df, self.settings.query_dir)]
        if not self.fsdb_target:
            qtdf_dir_iter.append((self.target_df, self.settings.target_dir))

        for df, search_dir in qtdf_dir_iter:
            records = df.drop_duplicates(subset=["preprocess_name", "chain_info"]).to_dict(orient="records")
            logging.debug(f"Records to preprocess: {json.dumps(records, indent=4)}", extra=stage)

            structure_preprocessor.set_output_directory(self.settings.foldseek_preprocessed_structure_dir)
            structure_preprocessor.update_cache()
            results = structure_preprocessor.preprocess_records(records=records, search_dir=search_dir)
            logging.debug(f"Preprocessing results: {json.dumps(results, indent=4)}", extra=stage)

            # Updating success and failure cols based on preprocessing results
            df.set_index("pocket_id", inplace=True)
            for index, success in results.items():
                if not success:
                    df.loc[index, "success"] = False
                    df.loc[index, "failure_reason"] = "structure_preprocessing_failed"
            df.reset_index(inplace=True)

        logging.info("Finished preprocessing structures", extra=stage)
        logging.debug(f"Query data after preprocessing: \n{self.query_df.head()}", extra=stage)
        logging.debug(f"Target data after preprocessing: \n{self.target_df.head()}", extra=stage)

    def foldseek_alignment(self):
        """
        Build the query (and, unless the target is a database, target) Foldseek DB, then search them.

        `foldseek easy-search` writes its raw tabular matches to `self.settings.alignment_path`.
        Every invocation goes through `foldseek.run_foldseek`, which logs and reports the failures.

        Returns:
            None
        """
        stage = {"stage": "Foldseek Alignment"}
        logging.info("Running Foldseek alignment...", extra=stage)

        # Setting up paths for foldseek databases
        self.query_db_path = os.path.join(self.settings.query_dir, "query_db")
        run_foldseek(["createdb", self.settings.query_dir, self.query_db_path], stage)

        if self.fsdb_target:
            self.target_db_path = self.target_df.loc[0, "struct_path"]
            logging.debug(f"Targeting bundled human_domains Foldseek DB at {self.target_db_path}", extra=stage)
        else:
            self.target_db_path = os.path.join(self.settings.target_dir, "target_db")
            run_foldseek(["createdb", self.settings.target_dir, self.target_db_path], stage)

        query_target_align_cmd = [
            "easy-search",
            self.query_db_path,
            self.target_db_path,
            self.settings.alignment_path,
            self.settings.foldseek_tmp_dir,
            "--format-output",
            FOLDSEEK_FORMAT_OUTPUT,
            "--format-mode",
            "4",
            "-e",
            "0.001",
            "--file-include",
            r".*\.cif\.gz",
            "--max-seqs",
            "5000",
            "-v",  # verbosity
            str(
                min(3, self.settings.verbosity)
            ),  # cap foldseek verbosity at 3 (info level) since it can be very verbose at higher levels and we already have our own logging verbosity control
        ]
        run_foldseek(query_target_align_cmd, stage)
        logging.debug("Foldseek alignment completed successfully", extra=stage)

    def local_alignment(self):
        """
        Execute traditional sequence-level alignment across inputs.

        Uses the `SequenceAligner` to locally match structured frames, resolving pairwise dependencies
        within internal dataframe representations. Outputs directly to the configuration TSV alignment file.

        Returns:
            None
        """
        stage = {"stage": "Local Alignment"}
        logging.info("Running local sequence alignments...", extra=stage)

        # Run the sequence aligner on the same preprocessed structures as foldseek uses
        aligner = SequenceAligner()
        query_records = (
            self.query_df.query("success").drop_duplicates(subset="preprocess_name").to_dict(orient="records")
        )
        target_records = (
            self.target_df.query("success").drop_duplicates(subset="preprocess_name").to_dict(orient="records")
        )
        alignment = aligner.align_records(
            query_records,
            target_records,
        )
        alignment.to_csv(self.settings.alignment_path, index=False, sep="\t")

    def dump_pockets(self, pockets, filename):
        """
        Write a pocket collection to the pocket cache directory as JSON, for inspection only.

        Nothing reads these files back -- pockets are recomputed every run -- so they are debug
        artefacts rather than a cache. `asdict` gives the nested shape declared in `pocket.py`:
        residues sit under a `residues` key rather than alongside the metadata. Written compact for
        the same reason: on a bundled-`pdb` Foldseek run `pisa_pockets.json` holds thousands of
        pockets, and nothing reads it to justify the whitespace.

        A None survives into the file as `null` rather than raising here. Three of the four builders
        store whatever their producer returned, including None for a missing structure or chain, and
        the failure for those belongs where it already is -- in `compare_pockets`.

        Args:
            pockets (dict): pocket_id -> Pocket (or None).
            filename (str): Basename to write under `pocket_dir`.

        Returns:
            None: Writes a file.
        """
        serialisable = {pid: asdict(pocket) if pocket is not None else None for pid, pocket in pockets.items()}
        with open(os.path.join(self.settings.pocket_dir, filename), "w") as f:
            json.dump(serialisable, f)

    def download_pisa_interfaces(self, pdb_list):
        """
        Download PISA summaries, assemblies and interfaces for `pdb_list` into the pocket cache.

        Both callers -- `expand_fsdb_pdb_targets` and `retrieve_pisa_pockets` -- go through here, so
        the two share one set of directories. Let them compute their own and the second re-downloads
        everything the first already fetched, behind PisaDownloader's per-entry sleep.

        Args:
            pdb_list (list): PDB IDs to fetch interfaces for.

        Returns:
            str: The interface directory -- the only one of the three anything reads back.
        """
        pisa_response_dir = os.path.join(self.settings.pocket_dir, "pisa_responses")
        interface_dir = os.path.join(pisa_response_dir, "interfaces")
        PisaDownloader().get_interfaces(
            pdb_list=pdb_list,
            summary_dir=os.path.join(pisa_response_dir, "summaries"),
            asm_dir=os.path.join(pisa_response_dir, "assemblies"),
            interface_dir=interface_dir,
        )
        return interface_dir

    def select_pocket_records(self, pocket_method, label, dedup_subset=None):
        """
        Select the successfully fetched records of one pocket method, from both sides at once.

        Args:
            pocket_method (str): The `pocket_method` value to match.
            label (str): Human-readable name of the method, used in the log lines.
            dedup_subset (list, optional): Columns to drop duplicate records on. Only `pisa` passes
                one: its pocket is fully determined by structure and chain, and the Foldseek-DB
                expansion generates the same pair once per hit. The other methods must not dedup on
                structure and chain -- two passthrough pockets on one chain differ only in
                `residue_info`. Defaults to None, meaning keep every record.

        Returns:
            pandas.DataFrame: The matching records, possibly empty.
        """
        stage = {"stage": f"Retrieving {label} Pockets"}
        logging.info(f"Checking for {label} pockets...", extra=stage)

        qt_df = pd.concat([self.query_df, self.target_df], ignore_index=True).query(
            f"success and pocket_method == '{pocket_method}'"
        )
        if dedup_subset is not None:
            qt_df = qt_df.drop_duplicates(subset=dedup_subset)

        if len(qt_df) == 0:
            logging.info(f"No {label} pockets to retrieve", extra=stage)
        else:
            logging.info(f"{len(qt_df)} {label} pockets to retrieve", extra=stage)
        return qt_df

    def get_pockets(self):
        """
        Build every pocket in the run, dispatching each record to its pocket method.

        Fans out over the `builders` table rather than over hand-written calls: one row per
        `pocket_method`, naming the log label, the dedup columns and the builder, with the key
        doubling as the dump filename stem. Every builder returns the same Pocket shape (see
        `pocket.py`), so nothing downstream needs to know which method produced a given pocket.

        Returns:
            dict: pocket_id -> Pocket, across both sides.
        """
        stage = {"stage": "Getting Pockets"}
        logging.info("Starting pocket retrieval...", extra=stage)

        # Turns PDB Foldseek-database hits into ordinary pisa target records, so the retrieval below
        # picks them up like any other pisa entry. No-op for every other kind of target.
        self.expand_fsdb_pdb_targets()

        # pocket_method -> (log label, drop_duplicates subset, builder). Insertion order is the
        # merge order. A new pocket method is one row here plus its builder, and nothing else.
        builders = {
            "pisa": ("PISA", ["struct_info", "chain_info"], self.retrieve_pisa_pockets),
            "passthrough": ("passthrough", None, self.retrieve_passthrough_pockets),
            "vdw": ("VDW", None, self.retrieve_vdw_pockets),
            "whole_chain": ("whole chain", None, self.retrieve_whole_chain_pockets),
        }

        pockets = {}
        for pocket_method, (label, dedup_subset, builder) in builders.items():
            qt_df = self.select_pocket_records(pocket_method, label, dedup_subset)
            if len(qt_df) == 0:
                continue
            method_pockets = builder(qt_df)
            self.dump_pockets(method_pockets, f"{pocket_method}_pockets.json")
            logging.debug(f"Extracted {label} pockets: {method_pockets}", extra=stage)
            pockets |= method_pockets

        logging.debug(f"Combined pockets: {pockets}", extra=stage)
        return pockets

    def expand_fsdb_pdb_targets(self):
        """
        Turn the hits of a PDB Foldseek-database search into ordinary PISA target records.

        A Foldseek-database target has no per-chain records of its own, so `compare_pockets` normally
        synthesises a whole-chain pseudo-pocket for each hit and leaves every `pocket_2_*` column empty.
        The PDB database is built from real PDB entries, though, so its hits have real PISA interfaces:
        this reads the hit names out of the alignment table, resolves each to a PDB ID and chain, asks
        PISA which chains that chain touches, and appends one `pisa` record per interface to
        `self.target_df`. `retrieve_pisa_pockets` then handles them like any other pisa entry.

        The generated records carry the Foldseek entry name as their `preprocess_name` rather than the
        one `QTProcessor` derives, because that is the key alignments are stored under -- it is what
        joins these pockets back to their alignment rows and to their Foldseek transforms.

        Runs after `alignment`, so `alignment.tsv` exists. Does nothing unless the target is a Foldseek
        database whose entries are named in the PDB style; a database of anything else (human_domains)
        keeps the synthesised whole-chain pockets.

        Returns:
            None: appends to `self.target_df` and sets `self.fsdb_pdb_target`.
        """
        if not self.fsdb_target:
            return

        stage = {"stage": "Expanding Foldseek DB Targets"}

        alignment_df = pd.read_csv(self.settings.alignment_path, sep="\t", engine="c")
        hits = {}  # foldseek entry name -> (pdb_id, chain_id)
        for hit_name in alignment_df["target"].unique().tolist():
            resolved = parse_foldseek_pdb_entry_name(hit_name)
            if resolved is not None:
                hits[hit_name] = resolved
        if not hits:
            logging.info(
                "Foldseek database target is not a PDB database; keeping whole-chain target pockets",
                extra=stage,
            )
            return
        self.fsdb_pdb_target = True

        pdb_list = sorted({pdb_id for pdb_id, _ in hits.values()})
        logging.info(
            f"Retrieving PISA interfaces for {len(pdb_list)} PDB entries behind {len(hits)} Foldseek hits",
            extra=stage,
        )

        # Shares retrieve_pisa_pockets' cache, so its own call is a no-op for everything fetched here.
        interface_dir = self.download_pisa_interfaces(pdb_list)

        # Building one record per interface the hit chain takes part in
        parser = PisaParser()
        qtprocessor = QTProcessor(
            structure_dir=self.settings.structure_dir,
            foldseek_preprocessed_structure_dir=self.settings.foldseek_preprocessed_structure_dir,
            fsdb_dir=self.settings.fsdb_dir,
        )
        records = []
        for hit_name, (pdb_id, chain_id) in hits.items():
            for partner in parser.get_interface_partners(pdb_id, chain_id, interface_dir):
                record = qtprocessor.parse_individual_qt(f"{pdb_id}:{chain_id}_{partner}", pocket_method="pisa")
                if record is None:
                    continue
                # The alignment is keyed by the Foldseek entry name, not by the name QTProcessor derives.
                # Nothing preprocesses these structures, so the preprocessing paths are meaningless here.
                record.preprocess_name = hit_name
                record.preprocess_path = None
                record.preprocess_path_gz = None
                records.append(asdict(record))
        if not records:
            logging.warning("No PISA interfaces found for any Foldseek hit", extra=stage)
            return

        # Fetching structures last, so only entries that actually produced a pocket are downloaded.
        # Hits whose structure can't be fetched come back marked success=False and are dropped here;
        # fetch_missing_structures raises only if not one of them could be fetched, which would leave
        # nothing to compare against at all.
        target_df = self.fetch_missing_structures("foldseek hit", pd.DataFrame(records), self.settings.structure_dir)
        target_df = target_df.query("success")

        logging.info(
            f"Added {len(target_df)} PISA target pockets from {target_df['preprocess_name'].nunique()} Foldseek hits",
            extra=stage,
        )
        # Row 0 stays the database record itself -- foldseek_alignment and align_structs read its struct_path.
        self.target_df = pd.concat([self.target_df, target_df], ignore_index=True)

    def retrieve_pisa_pockets(self, pisa_df):
        """
        Build a Pocket per record from the PDBe PISA interface it names.

        Downloads (or reuses from cache) the PISA files for every PDB entry in `pisa_df`, takes each
        record's pocket to be the residues its own chain contributes to the named interface, then
        reads the CA coordinates for those residues out of the structure. `PisaParser` skips records
        whose entry or interface cannot be resolved, so the result may be smaller than `pisa_df`.

        Args:
            pisa_df (pandas.DataFrame): Records with `pocket_method == "pisa"`, as returned by
                `select_pocket_records`.

        Returns:
            dict: pocket_id -> Pocket.
        """
        stage = {"stage": "Retrieving PISA Pockets"}

        pisa_pdb_list = pisa_df["struct_info"].unique().tolist()
        logging.debug(f"PDBs for which to retrieve PISA pockets: {pisa_pdb_list}", extra=stage)
        interface_dir = self.download_pisa_interfaces(pisa_pdb_list)

        parser = PisaParser()
        pisa_pockets = parser.get_pockets_from_records(records=pisa_df.to_dict(orient="records"), in_dir=interface_dir)
        logging.debug(f"PISA pockets before coordinates: {pisa_pockets}", extra=stage)

        # PisaParser gives residue ids but no geometry; this second pass fills in seq_pos and the CA
        # coordinates on the same Pocket, which is what the comparison and superposition need.
        for _, row in pisa_df.iterrows():
            if row["pocket_id"] in pisa_pockets:
                domain_chain, _ = split_chain_info(row["chain_info"])
                pisa_pockets[row["pocket_id"]] = parse_pocket_from_struct(
                    struct=row["struct_path"],
                    chain_id=domain_chain,
                    pocket_residues=[int(x) for x in pisa_pockets[row["pocket_id"]].res_auth_ids],
                    pocket=pisa_pockets[row["pocket_id"]],
                )
        return pisa_pockets

    def retrieve_passthrough_pockets(self, pt_df):
        """
        Build a Pocket from the residue ids the entry names outright.

        A passthrough entry carries its own residue list ("P24941:A:160,161"), so there is nothing to
        compute: the ids are read from `residue_info` and looked up in the structure.

        Args:
            pt_df (pandas.DataFrame): Records with `pocket_method == "passthrough"`, as returned by
                `select_pocket_records`.

        Returns:
            dict: pocket_id -> Pocket.
        """
        passthrough_pockets = {}
        for _, row in pt_df.iterrows():
            domain_chain, _ = split_chain_info(row["chain_info"])
            passthrough_pockets[row["pocket_id"]] = parse_pocket_from_struct(
                struct=row["struct_path"],
                chain_id=domain_chain,
                pocket_residues=[int(x) for x in row["residue_info"].split(",")],
            )
        return passthrough_pockets

    def retrieve_vdw_pockets(self, vdw_df):
        """
        Build a Pocket from the van der Waals contacts between the entry's two chains.

        A chain pair on a structure PISA cannot serve -- a local file or an AlphaFold model --
        resolves to this method. The pocket is the residues of the first chain whose atoms come
        within van der Waals contact of the second; see `PocketCalculator.pocket_overlap`.

        Args:
            vdw_df (pandas.DataFrame): Records with `pocket_method == "vdw"`, as returned by
                `select_pocket_records`.

        Returns:
            dict: pocket_id -> Pocket.
        """
        vdw_pockets = {}
        pc = PocketCalculator()
        for _, row in vdw_df.iterrows():
            domain_chain, motif_chain = split_chain_info(row["chain_info"])
            vdw_pockets[row["pocket_id"]] = pc.pocket_overlap(
                structure=row["struct_path"],
                domain_chain=domain_chain,
                motif_chain=motif_chain,
            )
        return vdw_pockets

    def retrieve_whole_chain_pockets(self, wc_df):
        """
        Build the "pocket" for an open search: every CA-bearing residue of the chain.

        An entry that names a structure but no pocket ("4Q5J:B", or "4Q5J" for the default chain) asks
        whether the query pocket resembles anything on that chain at all. Unlike the whole-chain pseudo-pocket
        `compare_pockets` synthesises for Foldseek-database hits, this is a real Pocket parsed from the
        structure, so it carries residue codes and CA coordinates and can be superposed.

        Args:
            wc_df (pandas.DataFrame): Records with `pocket_method == "whole_chain"`, as returned by
                `select_pocket_records`.

        Returns:
            dict: pocket_id -> Pocket.
        """
        stage = {"stage": "Retrieving whole chain Pockets"}

        whole_chain_pockets = {}
        for _, row in wc_df.iterrows():
            domain_chain, _ = split_chain_info(row["chain_info"])
            pocket = parse_pocket_from_struct(
                struct=row["struct_path"],
                chain_id=domain_chain,
                pocket_residues=None,  # None means the whole chain
            )
            # A missing structure or chain gives None back. Storing that would fail later with an opaque
            # TypeError inside compare_pockets, so drop the entry and say which one it was -- an unreadable
            # chain is much more likely here, where the chain can come from the default rather than the user.
            if pocket is None:
                logging.warning(
                    f"Could not parse chain {row['chain_info']} of {row['struct_info']} for {row['pocket_id']}, "
                    "skipping this entry",
                    extra=stage,
                )
                continue
            whole_chain_pockets[row["pocket_id"]] = pocket
        return whole_chain_pockets

    def compare_pockets_based_on_alignment(self, pockets):
        """
        Merge sequence or fold structure alignments traversing spatial pocket metrics internally to evaluate correlations.

        Processes raw mapping dictionaries tracking unresolved structural aliases. Utilizes a referenced BLOSUM
        scorecard for computing structural compatibility values output into a final result dataframe file.

        Args:
            pockets (dict): Populated pool of pocket residue locations collected via initial configurations.

        Returns:
            None
        """
        stage = {"stage": "Comparing Pockets Based on Alignment"}

        logging.info("Reading alignment results...", extra=stage)
        alignment_df = pd.read_csv(self.settings.alignment_path, sep="\t", engine="c")
        blosum_path = os.path.join(os.path.dirname(__file__), "blosum62.bla")

        logging.info(f"{len(alignment_df)} alignment pairs to compare", extra=stage)
        logging.debug(f"Alignment pairs: \n{alignment_df.head()}", extra=stage)

        preproc_to_ids = {}
        for _, row in pd.concat([self.query_df, self.target_df], ignore_index=True).iterrows():
            # Keyed by preprocess_name: one chain can carry several pockets (several pocket_ids), so
            # membership must be tested on the key, not on the pocket_id, or each chain keeps only its last.
            if row["preprocess_name"] in preproc_to_ids:
                if row["pocket_id"] not in preproc_to_ids[row["preprocess_name"]]:
                    preproc_to_ids[row["preprocess_name"]].append(row["pocket_id"])
            else:
                preproc_to_ids[row["preprocess_name"]] = [row["pocket_id"]]
        logging.debug(f"Preprocessed name to pocket ID mapping: {preproc_to_ids}", extra=stage)

        # A PDB Foldseek database has real PISA pockets for its hits (see expand_fsdb_pdb_targets);
        # every other database has no target records at all, so its pockets must be synthesised.
        synthesise_target_pockets = self.fsdb_target and not self.fsdb_pdb_target

        # Read from target_df rather than target_db_path so this does not depend on a step-4 side
        # effect. Row 0 is still the single database record here: expand_fsdb_pdb_targets only rewrites
        # target_df on the PDB path, which synthesise_target_pockets excludes.
        offset_table_path = None
        if synthesise_target_pockets:
            # Matched on the resolved path, so naming the bundled DB by its path and naming it
            # "human_domains" give the same answer.
            offset_paths = {
                db["db_path"]: db["offset_path"] for db in bundled_foldseek_dbs(self.settings.fsdb_dir).values()
            }
            offset_table_path = offset_paths.get(self.target_df.loc[0, "struct_path"])
            if offset_table_path is None:
                logging.info(
                    "Foldseek database ships no offset table; target residue ids will be positions "
                    "within each database entry rather than UniProt coordinates",
                    extra=stage,
                )
            elif not os.path.isfile(offset_table_path):
                # Packaged with the database (pyproject.toml's human_domains/*), so absence means a
                # broken install. Falling back would emit positional ids under a column that the
                # bundled database is documented to report in UniProt coordinates.
                msg = f"Bundled human-domains offset table is missing from the installation: {offset_table_path}"
                logging.critical(msg, extra=stage)
                raise PocketMapperError(msg)
            else:
                logging.info(
                    f"Mapping target residue ids to UniProt coordinates using {offset_table_path}", extra=stage
                )

        pockets_df, unknown_alias, incorrect_mapping = compare_pockets(
            alignment_df,
            pockets,
            preproc_to_ids=preproc_to_ids,
            blosum_path=blosum_path,
            synthesise_target_pockets=synthesise_target_pockets,
            offset_table_path=offset_table_path,
        )

        # Logging cases where a residue was given a single cahr name unfamiliar to pocketmapper
        if len(unknown_alias) > 0:
            unknown_alias_path = os.path.join(self.settings.results_dir, "unknown_ids.json")
            logging.warning("Unknown Foldseek Alias, see unknown_alias.json in results directory", extra=stage)
            with open(unknown_alias_path, "w") as f:
                json.dump(jsonify_dict(dict(unknown_alias)), f)

        # logging cases where foldseek mapping had low sequence identity to the parsed structure
        if len(incorrect_mapping) > 0:
            incorrect_mapping_path = os.path.join(self.settings.results_dir, "incorrect_mapping.json")
            logging.warning("Foldseek mapping with low sequence identity to parsed structure", extra=stage)
            with open(incorrect_mapping_path, "w") as f:
                json.dump(jsonify_dict(dict(incorrect_mapping)), f)

        # Writing pocket comparison results to output file
        output_path = self.settings.pocket_comparison_path
        pockets_df.to_csv(output_path, index=False, sep="\t")
        logging.info(f"Pocket comparison results saved to {output_path}", extra=stage)

    def align_structs(self):
        """
        Perform structural superposition of target structures against the query reference frame.

        Takes the top N targets (as defined by `align_count`) from the pocket comparison results
        and performs a structural alignment using the `StructureAligner` class.
        The aligned structures will be saved to the target directory for downstream analysis.

        Returns:
            None
        """
        stage = {"stage": "Structural Alignment"}
        if self.settings.align_count <= 0:
            logging.info("No Aligned Structures to Process", extra=stage)
            return

        method = self.settings.align_struct_method  # already "pocket" or "foldseek"
        logging.info(f"Performing structural alignment of target structures on the {method}...", extra=stage)

        # Pre-loading
        aligner = StructureAligner()
        pocket_comparison_df = pd.read_csv(self.settings.pocket_comparison_path, sep="\t", engine="c")
        alignment_df = None
        pocket_transform_df = None
        if method == "foldseek":
            alignment_df = pd.read_csv(
                self.settings.alignment_path,
                sep="\t",
                engine="c",
                index_col=["query", "target"],
            )
        else:
            # (pocket_1, pocket_2) is unique -- compare_pockets' existing_calcs scores each pair once.
            pocket_transform_df = pocket_comparison_df.dropna(subset=["p2_to_p1_u", "p2_to_p1_t"]).set_index(
                ["pocket_1", "pocket_2"]
            )[["p2_to_p1_u", "p2_to_p1_t"]]

        # For each query structure, align the top N target structures
        qt_id_map = {}
        unique_target_ids = set()
        for record in self.query_df.to_dict(orient="records"):
            query_id = record["pocket_id"]
            logging.debug(f"Processing query {query_id} for structural alignment", extra=stage)

            # Select the top N target structures based on pocket comparison metrics. Targets sharing no
            # pocket residues with the query are excluded: there is no common set of residues to superpose
            # on, and their overlap metrics are empty so they would sort arbitrarily.
            #
            # A whole-chain target -- an open search, or a Foldseek-DB hit -- has no jaccard_index, so it
            # sorts to the end and is ranked by the secondary key, min_overlap_similarity, instead.
            candidates = pocket_comparison_df.query(f"pocket_1 == '{query_id}' and overlap_count > 0")
            overlapping_count = len(candidates)
            if method == "pocket":
                # superpose fits nothing below three overlapping residues, so those targets have no
                # transform. Drop them here rather than when writing, or they would eat align_count
                # slots and the run would quietly produce fewer structures than asked for.
                candidates = candidates.dropna(subset=["p2_to_p1_u", "p2_to_p1_t"])

            target_ids = (
                candidates.sort_values(by=["jaccard_index", "min_overlap_similarity"], ascending=False)
                .head(self.settings.align_count)
                .loc[:, "pocket_2"]
                .to_list()
            )
            if not target_ids:
                if overlapping_count:
                    logging.info(
                        f"No target overlaps the pocket of query {query_id} by the three residues a "
                        "superposition needs; skipping its structural alignment",
                        extra=stage,
                    )
                else:
                    logging.info(
                        f"No target overlaps the pocket of query {query_id}; skipping its structural alignment",
                        extra=stage,
                    )
                continue
            logging.debug(f"Top target IDs for query {query_id}: {target_ids}", extra=stage)
            qt_id_map[query_id] = target_ids
            unique_target_ids.update(target_ids)

        if not unique_target_ids:
            logging.info("No query/target pair shares pocket residues, nothing to superpose", extra=stage)
            return

        self.query_df = self.query_df.set_index("pocket_id")
        if self.fsdb_target is False:
            target_record_df = self.target_df.set_index("pocket_id")
        else:  # If the target is a Foldseek database we need to make pdb structures from required entries
            source_db_path = self.target_df.loc[0, "struct_path"]
            logging.debug(f"Using Foldseek database at {source_db_path} for structural alignment", extra=stage)

            # With a PDB database the target IDs are pocket IDs ("4Q5J:B_F"), not database entry names, so
            # map them back through the records built by expand_fsdb_pdb_targets. One pocket ID can come
            # from more than one entry (the same chain in two assemblies) -- keep the first, or the lookup
            # below returns duplicate rows and the structure gets superposed twice.
            if self.fsdb_pdb_target:
                id_to_entry = (
                    self.target_df.dropna(subset=["preprocess_name"])
                    .drop_duplicates(subset="pocket_id", keep="first")
                    .set_index("pocket_id")["preprocess_name"]
                )
                target_entry_names = {target_id: id_to_entry[target_id] for target_id in unique_target_ids}
            else:
                target_entry_names = {target_id: target_id for target_id in unique_target_ids}

            # Get chain IDs corresponding to the required entries from the Foldseek database lookup file
            source_db_lookup_path = source_db_path + ".lookup"
            source_db_lookup_df = pd.read_csv(
                source_db_lookup_path, sep="\t", header=None, names=["chain_id", "name", "struct_id"]
            )
            source_db_lookup_df = source_db_lookup_df.set_index("name")
            chain_ids = source_db_lookup_df.loc[list(target_entry_names.values()), "chain_id"].tolist()

            # Make directory for subdb
            subdb_dir = os.path.join(self.settings.aligned_structure_dir, "fsdb")
            os.makedirs(subdb_dir, exist_ok=True)

            # Create a file listing the required chain IDs for the subdb creation
            subdb_chain_id_path = os.path.join(subdb_dir, "required_chain_ids.txt")
            with open(subdb_chain_id_path, "w") as f:
                for target_id in chain_ids:
                    f.write(f"{target_id}\n")

            # Create the subdb using foldseek's createsubdb command
            subdb_path = os.path.join(subdb_dir, "subdb")
            run_foldseek(["createsubdb", subdb_chain_id_path, source_db_path, subdb_path], stage)

            # Create a directory for extracted structures
            subdb_struct_dir = os.path.join(self.settings.aligned_structure_dir, "fsdb_structures")
            os.makedirs(subdb_struct_dir, exist_ok=True)

            # Convert the subdb to PDB format using foldseek's convert2pdb command
            run_foldseek(["convert2pdb", "--pdb-output-mode", "1", subdb_path, subdb_struct_dir], stage)

            # Make record df for the target records based on the unique target IDs and the subdb structure
            # directory. chain_info stays None: each extracted structure holds exactly the one chain of its
            # database entry, which foldseek_transform takes as the domain chain.
            target_record_df = pd.DataFrame(
                {"pocket_id": list(target_entry_names.keys()), "preprocess_name": list(target_entry_names.values())}
            )
            target_record_df["chain_info"] = None
            target_record_df["struct_path"] = target_record_df["preprocess_name"].apply(
                lambda x: os.path.join(subdb_struct_dir, f"{x}.pdb")
            )
            target_record_df = target_record_df.set_index("pocket_id")

        for query_id, target_ids in qt_id_map.items():
            # pocket_id is the index of query_df, so it is not in the row dict -- put it back, since
            # foldseek_transform reads it for the COMPND metadata.
            query_record = self.query_df.loc[query_id].to_dict()
            query_record["pocket_id"] = query_id
            logging.debug(f"Query record for '{query_id}': {json.dumps(query_record, indent=4)}", extra=stage)
            # Fetch the corresponding target records. A pocket_2 need not be a target: when a query and a
            # target share a chain they share a preprocess_name, so compare_pockets pairs every pocket on
            # that chain with every other and some rows come back with a query-only pocket_id in pocket_2.
            # Those have no target structure to superpose, so drop them rather than let .loc raise.
            known_target_ids = [target_id for target_id in target_ids if target_id in target_record_df.index]
            missing_target_ids = [target_id for target_id in target_ids if target_id not in target_record_df.index]
            if missing_target_ids:
                logging.debug(
                    f"Skipping non-target pocket(s) {missing_target_ids} when superposing onto '{query_id}'",
                    extra=stage,
                )
            top_target_records = target_record_df.loc[known_target_ids].reset_index().to_dict(orient="records")
            logging.debug(
                f"Top target records for query '{query_id}': {json.dumps(top_target_records, indent=4)}", extra=stage
            )

            # The query is the reference frame every target is superposed onto, so it must lead the list.
            aln_records = [query_record] + top_target_records
            if len(aln_records) > 1:
                out_path = os.path.join(self.settings.aligned_structure_dir, f"{safe_filename(query_id)}.pdb")
                if method == "foldseek":
                    aligner.foldseek_transform(
                        aln_records=aln_records,
                        alignment_df=alignment_df,
                        out_path=out_path,
                    )
                else:
                    # transforms is positional, not keyed by pocket_id: a self-comparison gives the
                    # query and a target the same pocket_id, so a dict would collide.
                    transforms = [None]  # the query is the reference frame, placed untransformed
                    for record in top_target_records:
                        try:
                            row = pocket_transform_df.loc[(query_id, record["pocket_id"])]
                        except KeyError:
                            transforms.append(None)
                            logging.warning(
                                f"No pocket superposition for {query_id} against {record['pocket_id']}",
                                extra=stage,
                            )
                            continue
                        transforms.append(parse_pocket_transform(row["p2_to_p1_u"], row["p2_to_p1_t"]))
                    aligner.transform(aln_records=aln_records, transforms=transforms, out_path=out_path)

    def delete_tmp(self):
        """
        Delete this run's scratch directories, unless `delete_tmp` says otherwise.

        Removes query_dir and target_dir, plus foldseek_tmp_dir when Foldseek did the aligning.
        Anything resolving outside cache_dir and results_dir is left alone and warned about: these
        three are settable by both the settings file and the command line, so a mistyped
        `--query_dir` would otherwise hand an unrelated directory to `shutil.rmtree`.

        Returns:
            None
        """
        self.log_extra.update({"stage": "Cleaning Up"})

        tmp_dirs = [
            "query_dir",
            "target_dir",
        ]
        if self.settings.foldseek:
            tmp_dirs.append("foldseek_tmp_dir")

        # Named rather than counted: which directories survive depends on the aligner, so a run kept
        # for inspection should say where its inputs actually are.
        if not self.settings.delete_tmp:
            kept = ", ".join(getattr(self.settings, dir_key) for dir_key in tmp_dirs)
            logging.info(f"delete_tmp is False; keeping {kept}", extra=self.log_extra)
            return

        roots = [self.settings.cache_dir, self.settings.results_dir]
        for dir_key in tmp_dirs:
            path = getattr(self.settings, dir_key)
            if not is_within(path, roots):
                logging.warning(
                    f"Not deleting {dir_key} {path}: it is outside cache_dir and results_dir. "
                    "Remove it yourself if that was intended.",
                    extra=self.log_extra,
                )
                continue
            shutil.rmtree(path)
