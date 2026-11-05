# Outputs

What a run writes, what it caches, and what every column of `pocket_comparison.tsv` means.

## Result files

Everything is written under `--results_dir`:

| File | Contents |
| --- | --- |
| `pocket_comparison.tsv` | The main result: one row per pocket pair. [Columns below](#pocket_comparisontsv-columns). |
| `alignment.tsv` | The chain alignments the comparison is built on, from Foldseek or the local aligner. |
| `aligned_structures/*.pdb` | The top `--align_count` targets superposed onto each query. Named after the query, so identify a file by its `MOLECULE` records rather than its filename. |
| `job_settings.json` | The fully resolved settings for the run. `parse` run on its own writes `parse_settings.json` instead. |
| `info.log` | The run log. |
| `pockets.json` | Every entry's pocket, keyed by pocket id (`null` where it could not be built), and `chains`: each aligned chain's name -> the pocket ids on it. Each pocket is an object of metadata fields plus a `residues` map keyed by author residue id (`100`, or `100A` with an insertion code). Versioned; a file from an older version is rejected, so rerun `pockets`. |
| `entries.tsv` | A row per query and target entry: `side`, `pocket_id` (the entry as typed), how it parsed (`struct_info`, `chain_info`, `residue_info`, `struct_type`, `preprocess_name`, `pocket_method`), `source` (the input or file it came from), and, once a step has dropped it, `failed_step`, `failed_reason` (`invalid_entry`, `structure_not_found`, `structure_preprocessing_failed`, `pocket_not_built`) and `failed_detail`. The table holds no paths: a step finds each structure from its `struct_type` and `struct_info`, in its own cache dirs. A Foldseek-database hit is no entry, so one dropped is only logged. |
| `unknown_ids.json` | Pocket residues whose one-letter code in the alignment differs from the structure's. Only written if any were seen, beside `pocket_comparison.tsv`. |
| `incorrect_mapping.json` | Pockets dropped because their own sequence differed from the aligner's for that chain (case ignored; an `X` on either side matches) — typically assembly vs asymmetric unit numbering. Each record gives `side`, `reason`, `first_mismatch` (0-based), `pocket_seq` and `aligner_seq`. Only written if any were dropped, beside `pocket_comparison.tsv`. |

## Cache

Under `--cache_dir` and surviving between runs: the downloaded mmCIF files (`pdb_structures/` and
`alphafold_structures/`), the single-chain copies Foldseek is given
(`foldseek_preprocessed_structures/`), the downloaded `pdb` Foldseek database (`fsdb/`) and the raw PISA
responses (`pockets/pisa/`). Pockets are not cached: every run builds them afresh into `pockets.json`.

## `pocket_comparison.tsv` columns

The `query_*` columns describe the query pocket, the `target_*` columns the target. The table always
has all of these columns: a comparison that stops early — no overlapping residues, no coordinates,
fewer than three residues to superpose, or an open target — leaves the remaining fields **empty**
rather than dropping them.

### Identity

| Column | Meaning |
| --- | --- |
| `query`, `target` | The two pocket ids, in the input-entry form (`4Q5J:B_F`). |
| `evalue` | Alignment E-value for the underlying chain pair. `-` with the local aligner. |
| `lddt` | Alignment lDDT reported by Foldseek. `-` with the local aligner. |

### Pocket description

| Column | Meaning |
| --- | --- |
| `query_res_ids`, `target_res_ids` | The pocket's author residue ids, comma-separated: the residue number with its insertion code, if any (`100A`). |
| `query_len`, `target_len` | Number of residues in the pocket. |
| `query_seq`, `target_seq` | The pocket's residues as single-letter codes, in `res_ids` order. |
| `query_pct_aln`, `target_pct_aln` | Fraction of the pocket's residues that fall inside the aligned region at all. A low value means the pocket sits largely outside the alignment. |

### Overlap

| Column | Meaning |
| --- | --- |
| `overlap_count` | How many alignment positions both pockets occupy. The headline number. |
| `query_overlap_ids` | The query residues in the overlap, as author residue ids, in alignment order. |
| `target_overlap_ids` | The target residues they map onto, in alignment order and paired with `query_overlap_ids` by index. UniProt numbers for a `human_domains` target; see [Databases](input.md#databases). |
| `jaccard_index` | `overlap_count` divided by the union of the two pockets. Empty for an open or database target. |

### Overlap scoring (BLOSUM62)

| Column | Meaning |
| --- | --- |
| `query_seq_overlap`, `target_seq_overlap` | The two overlap sequences, aligned position for position. |
| `overlap_identity` | Fraction of overlap positions where the two residues are identical. |
| `overlap_similarity_binary` | Fraction of overlap positions with a positive BLOSUM62 score — how much of the overlap is conservatively substituted, ignoring how strongly. Positions where the query residue is X (any unknown residue) are skipped; empty when every position is skipped. |
| `overlap_similarity_1_2` | Mean BLOSUM62 score over the overlap, normalised per position by the query residue's self-score: in -1 to 1, and 1.0 for an identical pair. Positions where the query residue is X are skipped; empty when every position is skipped. |
| `overlap_similarity_2_1` | The same, normalised against the target residue instead, so positions where the target residue is X are skipped. The two differ because self-scores differ between residues. |
| `min_overlap_similarity`, `max_overlap_similarity` | The smaller and larger of the two directional scores. An empty side is ignored. |

### Superposition

Populated only when both pockets have coordinates and at least three residues overlap.

| Column | Meaning |
| --- | --- |
| `target_to_query_u`, `target_to_query_t` | Rotation (nine values) and translation (three) that put the target pocket onto the query pocket, fitted on the overlapping CA atoms. Comma-separated at six decimals (`alignment.tsv`'s `u` and `t` have three). |
| `query_to_target_u`, `query_to_target_t` | The same fit in the other direction. |
| `rmsd` | RMSD of the overlapping CA atoms after superposition, in Å. |
| `ca_dists` | Per-residue CA distance after superposition, comma-separated, in overlap order. |

The rotation matrices are stored in Biopython's right-multiplying convention
(`coords · u + t`). If you feed one to gemmi, which left-multiplies, transpose it first — or use
`pocketmapper.compare.pocket_comparison.parse_pocket_transform`, which does the conversion for you.
