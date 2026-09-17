"""
Building of aligned multi-structure PDB files.

Takes one query structure plus one or more matched targets and writes them back out as a single
human-readable PDB, as glue around gemmi: read mmCIF, apply per-target rotations and translations,
rename chains so each structure stays distinguishable, and emit a combined PDB with COMPND
metadata naming each model.

`StructureAligner.align_structs` drives the whole step from a finished pocket comparison, so a caller
holding one of those needs nothing else from the pipeline to produce the aligned structures.
"""

import json
import logging
import os
import string
from itertools import count
from itertools import permutations

import gemmi
import numpy as np
import pandas as pd

from pocketmapper.constants import RESOLVED_ALIGN_STRUCT_METHODS
from pocketmapper.exceptions import PocketMapperError
from pocketmapper.foldseek import extract_fsdb_structures
from pocketmapper.lib import safe_filename
from pocketmapper.lib import split_chain_info
from pocketmapper.pocket_comparison import parse_pocket_transform

logger = logging.getLogger(__name__)


def as_dataframe(table, **read_csv_kwargs):
    """
    Accept a result table either as a DataFrame or as the path of the TSV holding it.

    Args:
        table (pandas.DataFrame or str): The table itself, or a path to read it from.
        **read_csv_kwargs: Extra arguments for `pandas.read_csv`, used only when reading a path.

    Returns:
        pandas.DataFrame: `table` unchanged, or the parsed contents of the file it named.
    """
    if isinstance(table, pd.DataFrame):
        return table
    return pd.read_csv(table, sep="\t", engine="c", **read_csv_kwargs)


class StructureAligner:
    """
    Superposes structures onto a reference and writes them as one PDB.

    `align_structs` is the whole step: it picks the best targets for each query out of a pocket
    comparison table, finds a transform for each and writes one PDB per query.

    Below it sit the two transform sources the `align_struct_method` setting selects: `transform`
    applies transforms it is handed ("pocket"), and `foldseek_transform` sources them from a Foldseek
    alignment table first ("foldseek"). Either is usable on its own when the caller has already chosen
    which structures to superpose.
    """

    def __init__(self):
        """
        Initialise the logging stage. The aligner holds no other state.
        """
        self.log_extra = {"stage": "Structural Alignment"}
        logger.debug("Initialized")

    def char_gen(self):
        """
        Yield short, PDB-friendly chain identifiers for aligned output.

        The output structure may hold more chains than any one input, so names are generated rather than
        reused. "0" is skipped throughout: it is reserved for the domain chain, which
        `apply_transformation` renames so it can be found consistently across models.

        Single characters are emitted first, then two-character pairs once that space is exhausted.

        Yields:
            str: A unique chain label suitable for writing into the final PDB.
        """
        nice_chars = string.digits + string.ascii_letters
        for x in nice_chars[1:]:  # 0 is reserved for domain names
            yield x
        for x, y in permutations(nice_chars, 2):
            yield (x + y)

    def apply_transformation(self, structs, domain_chains, motif_chains, us, ts):
        """
        Apply rigid-body transforms to a set of structures and merge them.

        The first structure is the reference frame. The domain chain of each input is renamed to "0" in
        the output so it can be recognised consistently; any motif chain gets a generated name.

        Args:
            structs (list): gemmi.Structures to align and merge.
            domain_chains (list): Domain chain id to preserve from each input structure.
            motif_chains (list): Motif chain id per input, or None where there is none.
            us (list): 3x3 rotation matrices, one per structure, in gemmi's LEFT-multiplying convention.
            ts (list): Translation vectors of length 3, one per structure.

        Returns:
            gemmi.Structure: A merged structure holding one model per input.
        """
        # Align everything to the first struct
        ref_st = gemmi.Structure()
        chain_names = self.char_gen()

        for i, cn, st, dc, mc, u, t in zip(count(1), chain_names, structs, domain_chains, motif_chains, us, ts):
            ref_st.add_model(gemmi.Model(i))

            # Apply transformation
            trans = gemmi.Transform(gemmi.Mat33(u), gemmi.Vec3(*t))
            st[0].transform_pos_and_adp(trans)

            # Add chains to the reference structure
            st[0][dc].name = "0"
            ref_st[-1].add_chain(st[0]["0"])
            if mc is not None:
                st[0][mc].name = cn
                ref_st[-1].add_chain(st[0][cn])

        ref_st.setup_entities()
        return ref_st

    def transform(self, aln_records, transforms, out_path):
        """
        Build an aligned multi-structure PDB from ready-made rigid-body transforms.

        The general entry point: it knows nothing about where a transform came from, only how to apply it.
        `foldseek_transform` sources them from Foldseek's whole-chain alignment; `align_structs` sources
        them from the pocket superposition in a pocket comparison table.

        The first record is the reference frame and is always placed untransformed, so `transforms[0]` is
        ignored. Every other entry is applied to its record's whole structure.

        `transforms` is positional rather than keyed by `pocket_id` on purpose: a query can be compared
        against itself, so the reference and a target may carry the same `pocket_id`.

        Args:
            aln_records (list): Ordered records, reference first. Fields read: `pocket_id`,
                `struct_path`, `chain_info`.
            transforms (list): Parallel to `aln_records`. Entry 0 is ignored; entry i is either a
                `(u, t)` pair -- a 3x3 rotation in gemmi's LEFT-multiplying convention plus a 3-vector
                translation -- or None to drop that record from the output. A caller passing None is
                expected to have logged why.
            out_path (str): Destination path for the aligned PDB file.

        Returns:
            str: `out_path`, once the merged PDB is written; None when no structure could be placed
                and nothing was written.
        """
        structs = []
        domain_chains = []
        motif_chains = []
        us = []
        ts = []
        kept_records = []
        dropped = []

        for i, record in enumerate(aln_records):
            try:
                if i == 0:
                    struct_u = np.eye(3)
                    struct_t = np.zeros(3)
                else:
                    if transforms[i] is None:
                        dropped.append(record["pocket_id"])
                        continue
                    struct_u, struct_t = transforms[i]

                struct = gemmi.read_structure(record["struct_path"])

                if record["chain_info"] is None:
                    domain_chain = 0  # first chain
                    motif_chain = None
                else:
                    domain_chain, motif_chain = split_chain_info(record["chain_info"])

                # If everything has been successful add it things to be processed
                structs.append(struct)
                us.append(struct_u)
                ts.append(struct_t)
                domain_chains.append(domain_chain)
                motif_chains.append(motif_chain)
                kept_records.append(record)

            except Exception as e:
                dropped.append(record["pocket_id"])
                logger.error(f"Problem processing {record['pocket_id']}: {e}", extra=self.log_extra)

        if dropped:
            logger.warning(
                f"Not superposing {dropped}; they are absent from {out_path}",
                extra=self.log_extra,
            )
        if not kept_records:
            logger.error(f"No structure could be placed, not writing {out_path}", extra=self.log_extra)
            return None

        aligned_struct = self.apply_transformation(structs, domain_chains, motif_chains, us, ts)
        self.write_aligned(kept_records, aligned_struct, out_path)
        return out_path

    def write_aligned(self, kept_records, aligned_struct, out_path):
        """
        Write a merged structure out as a PDB with a COMPND header naming each model.

        Takes the records that actually made it into `aligned_struct`, not the records the caller started
        with: a record dropped for want of a transform used to keep its COMPND entry, so the header named
        models the file did not contain.

        The chain labels must come from a fresh `char_gen()` consumed in the same order
        `apply_transformation` consumed its own, or the header and the coordinates disagree.

        Args:
            kept_records (list): Records present in `aligned_struct`, in model order.
            aligned_struct (gemmi.Structure): The merged structure to write.
            out_path (str): Destination path for the PDB file.

        Returns:
            None: Writes the PDB to `out_path`.
        """
        pdb_str = aligned_struct.make_pdb_string()

        model_nums = (str(x) for x in count(1))
        model_names = [record["pocket_id"] for record in kept_records]
        chain_names = self.char_gen()
        header = ""
        for model_num, model_name, chain_name in zip(model_nums, model_names, chain_names):
            line_nums = (str(x) for x in count(1))
            header += f"""
COMPND {next(line_nums).zfill(3)} MOL_ID: {model_num};
COMPND {next(line_nums).zfill(3)} MOLECULE: {model_name[:70]};
COMPND {next(line_nums).zfill(3)} CHAIN: {chain_name};
"""

        with open(out_path, "w") as f:
            f.write(header)
            f.write(pdb_str)

    def foldseek_transform(self, aln_records, alignment_df, out_path):
        """
        Build an aligned multi-structure PDB from Foldseek-style alignment results.

        The first record in `aln_records` is the reference; every other must have a row in `alignment_df`.
        The stored Foldseek transform strings are already in the LEFT-multiplying convention `transform`
        wants, so they are parsed and passed straight through.

        The local BLOSUM62 aligner writes "-" for `u` and `t`, so every target is dropped here and the
        output holds the query alone -- use the pocket transforms with `transform` instead (see the
        `align_struct_method` setting).

        Args:
            aln_records (list): Ordered alignment records, reference first. Fields read: `pocket_id`,
                `preprocess_name`, `struct_path`, `chain_info`.
            alignment_df (pandas.DataFrame): Alignment table carrying the Foldseek transforms, indexed by
                (query, target) `preprocess_name`.
            out_path (str): Destination path for the aligned PDB file.

        Returns:
            str: `out_path`, once the merged PDB is written; None when no structure could be placed
                and nothing was written.
        """
        query_preprocess_name = aln_records[0]["preprocess_name"]

        transforms = [None]  # the reference is placed untransformed
        for record in aln_records[1:]:
            try:
                row = alignment_df.loc[query_preprocess_name, record["preprocess_name"]]
                struct_u = np.array([float(x) for x in row["u"].split(",")]).reshape((3, 3))
                struct_t = np.array([float(x) for x in row["t"].split(",")])
                transforms.append((struct_u, struct_t))
            except Exception as e:
                transforms.append(None)
                logger.error(f"Problem processing {record['pocket_id']}: {e}", extra=self.log_extra)

        return self.transform(aln_records, transforms, out_path)

    def select_targets(self, query_ids, pocket_comparison_df, align_count, method):
        """
        Pick the targets to superpose onto each query, best first.

        Targets sharing no pocket residues with the query are excluded: there is no common set of
        residues to superpose on, and their overlap metrics are empty so they would sort arbitrarily.
        Ranking is by `jaccard_index` then `min_overlap_similarity`; a whole-chain target -- an open
        search, or a Foldseek-DB hit -- has no jaccard_index, so it sorts to the end and is ranked by
        the secondary key instead.

        Args:
            query_ids (list): Query `pocket_id`s to select for.
            pocket_comparison_df (pandas.DataFrame): Pocket comparison table.
            align_count (int): Most targets to keep per query.
            method (str): "pocket" or "foldseek"; "pocket" additionally drops targets with no
                superposition.

        Returns:
            dict: Query `pocket_id` -> its target `pocket_id`s, best first. Queries with no usable
                target are absent.
        """
        qt_id_map = {}
        for query_id in query_ids:
            logger.debug(f"Processing query {query_id} for structural alignment", extra=self.log_extra)
            candidates = pocket_comparison_df[
                (pocket_comparison_df["pocket_1"] == query_id) & (pocket_comparison_df["overlap_count"] > 0)
            ]
            overlapping_count = len(candidates)
            if method == "pocket":
                # superpose fits nothing below three overlapping residues, so those targets have no
                # transform. Drop them here rather than when writing, or they would eat align_count
                # slots and the run would quietly produce fewer structures than asked for.
                candidates = candidates.dropna(subset=["p2_to_p1_u", "p2_to_p1_t"])

            target_ids = (
                candidates.sort_values(by=["jaccard_index", "min_overlap_similarity"], ascending=False)
                .head(align_count)
                .loc[:, "pocket_2"]
                .to_list()
            )
            if not target_ids:
                if overlapping_count:
                    logger.info(
                        f"No target overlaps the pocket of query {query_id} by the three residues a "
                        "superposition needs; skipping its structural alignment",
                        extra=self.log_extra,
                    )
                else:
                    logger.info(
                        f"No target overlaps the pocket of query {query_id}; skipping its structural alignment",
                        extra=self.log_extra,
                    )
                continue
            logger.debug(f"Top target IDs for query {query_id}: {target_ids}", extra=self.log_extra)
            qt_id_map[query_id] = target_ids
        return qt_id_map

    def fsdb_target_records(self, fsdb_path, target_records, target_ids, out_dir, threads):
        """
        Build target records by rebuilding the needed structures out of a Foldseek database.

        How a target id names a database entry depends on the database, and the two cases are told
        apart by whether `target_records` holds anything at all -- never per id, or one column would
        mix entries resolved two different ways:

        - With records (a PDB database, whose hits were expanded into real pockets), an id is a
          `pocket_id` and its record's `preprocess_name` is the entry name. One pocket id can come from
          more than one entry -- the same chain in two assemblies -- so the first is kept; an id no
          record covers is dropped, since the caller filters those out again by structure.
        - Without them, each id is itself an entry name.

        Args:
            fsdb_path (str): Path to the Foldseek database.
            target_records (list): Target records, or an empty list when the ids are entry names.
            target_ids (list): Target ids needing a structure.
            out_dir (str): Directory to extract the structures under.
            threads (int): Thread count for the extraction.

        Returns:
            dict: Target `pocket_id` -> a record carrying `preprocess_name`, `struct_path` and a None
                `chain_info`.
        """
        if target_records:
            named = [record for record in target_records if not pd.isna(record.get("preprocess_name"))]
            entry_of_id = self.records_by_id(named)
            entry_names = {
                target_id: entry_of_id[target_id]["preprocess_name"]
                for target_id in target_ids
                if target_id in entry_of_id
            }
        else:
            entry_names = {target_id: target_id for target_id in target_ids}

        logger.debug(f"Using Foldseek database at {fsdb_path} for structural alignment", extra=self.log_extra)
        struct_paths = extract_fsdb_structures(
            fsdb_path,
            list(dict.fromkeys(entry_names.values())),
            out_dir,
            threads,
            self.log_extra,
        )

        # chain_info stays None: each extracted structure holds exactly the one chain of its database
        # entry, which `transform` takes as the domain chain.
        return {
            target_id: {
                "pocket_id": target_id,
                "preprocess_name": entry_name,
                "chain_info": None,
                "struct_path": struct_paths[entry_name],
            }
            for target_id, entry_name in entry_names.items()
        }

    def pocket_transforms(self, query_id, target_records, pocket_transform_df):
        """
        Look up the pocket superposition of each target against one query.

        Args:
            query_id (str): The query's `pocket_id`.
            target_records (list): Target records, in the order they will be superposed.
            pocket_transform_df (pandas.DataFrame): `p2_to_p1_u`/`p2_to_p1_t`, indexed by
                (`pocket_1`, `pocket_2`).

        Returns:
            list: One entry longer than `target_records` -- a leading None for the reference frame,
                then a `(u, t)` pair per target, or None where the pair has no superposition.
        """
        # Positional rather than keyed by pocket_id: a self-comparison gives the query and a target the
        # same pocket_id, so a dict would collide.
        transforms = [None]  # the query is the reference frame, placed untransformed
        for record in target_records:
            try:
                row = pocket_transform_df.loc[(query_id, record["pocket_id"])]
            except KeyError:
                transforms.append(None)
                logger.warning(
                    f"No pocket superposition for {query_id} against {record['pocket_id']}",
                    extra=self.log_extra,
                )
                continue
            transforms.append(parse_pocket_transform(row["p2_to_p1_u"], row["p2_to_p1_t"]))
        return transforms

    def records_by_id(self, records, id_field="pocket_id"):
        """
        Index records by an id field, keeping the first record carrying each id.

        Args:
            records (list): The records to index.
            id_field (str, optional): Field to key on. Defaults to "pocket_id".

        Returns:
            dict: Id -> the first record carrying it.
        """
        by_id = {}
        for record in records:
            by_id.setdefault(record[id_field], record)
        return by_id

    def warn_unknown_ids(self, side, wanted_ids, known_ids):
        """
        Log any requested pocket id that nothing in this run knows about.

        Args:
            side (str): "query" or "target", for the message.
            wanted_ids (list): The ids the caller asked to limit to.
            known_ids (container): The ids actually available.

        Returns:
            None
        """
        unknown = [wanted_id for wanted_id in wanted_ids if wanted_id not in known_ids]
        if unknown:
            logger.warning(f"No {side} in this run matches {unknown}", extra=self.log_extra)

    def align_structs(
        self,
        query_records,
        target_records,
        pocket_comparison,
        out_dir,
        method,
        align_count=10,
        alignment=None,
        threads=1,
        fsdb_path=None,
        query_ids=None,
        target_ids=None,
        overwrite=True,
    ):
        """
        Superpose the best targets of a finished pocket comparison onto each query.

        Writes one PDB per query into `out_dir`, named by `lib.safe_filename` of the query's
        `pocket_id` -- so a file is identified by the `MOLECULE` records inside it, not by its name.
        A query whose targets all turn out to be unusable is skipped rather than written empty.

        Args:
            query_records (list): Query records. Fields read: `pocket_id`, `struct_path`,
                `chain_info`, and `preprocess_name` for the "foldseek" method.
            target_records (list): Target records, same fields. Empty is allowed only alongside
                `fsdb_path`, where it means the target ids are database entry names.
            pocket_comparison (pandas.DataFrame or str): The pocket comparison table, or the path of
                the TSV holding it.
            out_dir (str): Directory to write the aligned PDBs into. Created if missing.
            method (str): Which transform to superpose with -- "pocket", read from the comparison
                table, or "foldseek", read from `alignment`.
            align_count (int, optional): Most targets to superpose onto each query. Defaults to 10;
                zero or less writes nothing.
            alignment (pandas.DataFrame or str, optional): The Foldseek alignment table, or the path
                of the TSV holding it. Required by the "foldseek" method and ignored otherwise.
            threads (int, optional): Thread count for rebuilding structures out of a Foldseek
                database. Defaults to 1 and is ignored without `fsdb_path`.
            fsdb_path (str, optional): Foldseek database to rebuild the target structures from, for a
                run whose targets were a database rather than structures of their own.
            query_ids (list, optional): Only superpose onto these queries. Defaults to None, all of
                them.
            target_ids (list, optional): Only consider these targets. Defaults to None, all of them.
            overwrite (bool, optional): Defaults to True. False keeps any PDB already in `out_dir` and
                skips that query entirely, so its targets cost nothing either.

        Returns:
            dict: Query `pocket_id` -> the path written for it. Queries that produced no file are
                absent, so an empty dict means nothing was written.

        Raises:
            PocketMapperError: If `method` is neither "pocket" nor "foldseek".
        """
        if method not in RESOLVED_ALIGN_STRUCT_METHODS:
            msg = (
                f"Unknown align_struct_method {method!r}. "
                f"Choose one of: {', '.join(RESOLVED_ALIGN_STRUCT_METHODS)}."
            )
            logger.critical(msg, extra=self.log_extra)
            raise PocketMapperError(msg)
        if align_count <= 0:
            logger.info("No Aligned Structures to Process", extra=self.log_extra)
            return {}

        logger.info(f"Performing structural alignment of target structures on the {method}...", extra=self.log_extra)
        os.makedirs(out_dir, exist_ok=True)

        # Pre-loading
        pocket_comparison_df = as_dataframe(pocket_comparison)
        alignment_df = None
        pocket_transform_df = None
        if method == "foldseek":
            alignment_df = as_dataframe(alignment)
            if "query" in alignment_df.columns:
                alignment_df = alignment_df.set_index(["query", "target"])
        else:
            # (pocket_1, pocket_2) is unique -- compare_pockets' existing_calcs scores each pair once.
            pocket_transform_df = pocket_comparison_df.dropna(subset=["p2_to_p1_u", "p2_to_p1_t"]).set_index(
                ["pocket_1", "pocket_2"]
            )[["p2_to_p1_u", "p2_to_p1_t"]]

        # A record list can hold the same pocket twice -- the same entry given twice on the command
        # line. Keying by pocket_id keeps the first of each, so a query is not written twice and a
        # target is not superposed twice into the same file.
        query_by_id = self.records_by_id(query_records)
        if query_ids is not None:
            self.warn_unknown_ids("query", query_ids, query_by_id)
            query_by_id = {query_id: query_by_id[query_id] for query_id in query_ids if query_id in query_by_id}
        if target_ids is not None:
            self.warn_unknown_ids("target", target_ids, set(pocket_comparison_df["pocket_2"]))
            pocket_comparison_df = pocket_comparison_df[pocket_comparison_df["pocket_2"].isin(list(target_ids))]

        out_paths = {}
        if not overwrite:
            for query_id in list(query_by_id):
                out_path = os.path.join(out_dir, f"{safe_filename(query_id)}.pdb")
                if os.path.exists(out_path):
                    logger.info(
                        f"Keeping the existing {out_path}; not superposing onto {query_id} again",
                        extra=self.log_extra,
                    )
                    out_paths[query_id] = out_path
                    del query_by_id[query_id]

        qt_id_map = self.select_targets(list(query_by_id), pocket_comparison_df, align_count, method)
        unique_target_ids = list(dict.fromkeys(t for targets in qt_id_map.values() for t in targets))
        if not unique_target_ids:
            logger.info("No query/target pair shares pocket residues, nothing to superpose", extra=self.log_extra)
            return out_paths

        if fsdb_path is None:
            target_by_id = self.records_by_id(target_records)
        else:
            target_by_id = self.fsdb_target_records(fsdb_path, target_records, unique_target_ids, out_dir, threads)

        for query_id, query_target_ids in qt_id_map.items():
            query_record = query_by_id[query_id]
            logger.debug(f"Query record for '{query_id}': {json.dumps(query_record, indent=4)}", extra=self.log_extra)

            # A pocket_2 need not be a target: when a query and a target share a chain they share a
            # preprocess_name, so compare_pockets pairs every pocket on that chain with every other and
            # some rows come back with a query-only pocket_id in pocket_2. Those have no target
            # structure to superpose, so drop them.
            missing_target_ids = [t for t in query_target_ids if t not in target_by_id]
            if missing_target_ids:
                logger.debug(
                    f"Skipping non-target pocket(s) {missing_target_ids} when superposing onto '{query_id}'",
                    extra=self.log_extra,
                )
            top_target_records = [target_by_id[t] for t in query_target_ids if t in target_by_id]
            logger.debug(
                f"Top target records for query '{query_id}': {json.dumps(top_target_records, indent=4)}",
                extra=self.log_extra,
            )
            if not top_target_records:
                continue

            # The query is the reference frame every target is superposed onto, so it must lead the list.
            aln_records = [query_record] + top_target_records
            out_path = os.path.join(out_dir, f"{safe_filename(query_id)}.pdb")
            if method == "foldseek":
                written = self.foldseek_transform(
                    aln_records=aln_records,
                    alignment_df=alignment_df,
                    out_path=out_path,
                )
            else:
                transforms = self.pocket_transforms(query_id, top_target_records, pocket_transform_df)
                written = self.transform(aln_records=aln_records, transforms=transforms, out_path=out_path)
            if written is not None:
                out_paths[query_id] = written

        return out_paths
