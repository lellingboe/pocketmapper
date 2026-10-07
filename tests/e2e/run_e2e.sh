#!/usr/bin/env bash
#
# End-to-end test cases for PocketMapper.
#
# Each case runs the real CLI against real remote services (wwPDB, AlphaFold,
# PDBe PISA) and asserts that the run exits cleanly and writes the expected
# result files. There are no mocks -- these are smoke tests for the whole
# pipeline, not unit tests.
#
# Run `./run_e2e.sh --help` for usage.

# Must be executed, never sourced: `set -u` below -- and every `exit` further
# down -- would otherwise apply to the calling shell. Under Terminal.app that
# surfaces as "-bash: HISTTIMEFORMAT: unbound variable" at each prompt (its
# per-prompt history hook reads that unset variable) and --list/--help kill the
# login shell outright, leaving a dead window.
if [ "${BASH_SOURCE[0]}" != "$0" ]; then
    echo "run_e2e.sh must be run, not sourced -- use: ./run_e2e.sh $*" >&2
    return 2
fi

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FIXTURES_DIR="$SCRIPT_DIR/fixtures"

# ---------------------------------------------------------------------------
# Test cases
#
#   name | tags | expect | description | args
#
# tags     space-separated; used by --tag and to gate cases needing extra
#          resources (see `needs-*` handling below).
# expect   space-separated tokens, all asserted. At most one of:
#          `rows` = run must succeed AND pocket_comparison.tsv must hold at
#          least one data row; `ok` = run must succeed and write the file, but
#          zero hits is a legitimate outcome for that pair; `fail` = the run
#          must exit non-zero (a rejected input or option combination);
#          `queries=N` = run must succeed AND the `query` column of
#          pocket_comparison.tsv must hold at least N distinct values, for a
#          case where losing one query's rows is the failure. With none of
#          them, the run must succeed. Plus any of: `files=A,B` = each path
#          exists under the case's results dir; `failed=REASON@POCKET_ID` =
#          failed_entries.json lists that entry with that reason;
#          `same=DIR` = pocket_comparison.tsv and alignment.tsv match, once
#          sorted, those in the results dir's subdirectory DIR.
# args     passed to `pocketmapper search` verbatim (word-split on spaces).
#          Or, when the first word is a command (parse, fetch_structures,
#          align, pockets, compare, superpose, search), a chain of commands
#          separated by ` ; `, run in order: every one before the last must
#          succeed, and the last one's exit status is the run's. A `cd DIR`
#          segment sets the working directory of the commands after it. Each
#          command gets --verbosity
#          and, unless it sets them, --results_dir and (parse, search)
#          --cache_dir.
#          @PDB_FSDB@ expands to $POCKETMAPPER_PDB_FSDB, @CACHE@ to the shared
#          cache dir and @OUT@ to this case's own results dir -- the last two
#          let a case aim a path option somewhere real without hardcoding a
#          machine-specific path. Foldseek is the CLI's default, so a case
#          running search or align is assumed to need the binary and is
#          skipped when it is missing; a case opts out with an explicit
#          `--aligner seq`, which exercises the local BLOSUM62 aligner and
#          still runs without the binary.
#
# Cases are grouped by what they exercise, and each group is named by its
# prefix and numbered within itself: test_core_* structure-vs-structure pairs,
# test_open_* open whole-chain targets, test_domains_* human_domains DB
# targets, test_fsdb_* the larger Foldseek DB targets, test_local_* the local
# aligner, test_invalid_* input that must be rejected or skipped rather than
# compared, test_settings_* how a run is configured rather than what it
# computes, test_steps_* the step commands run on their own. The prefix tracks
# the group, not the tag column -- test_open_*, test_invalid_*, the
# test_settings_* and most of the test_local_* and test_steps_* cases are
# tagged 'core' as well.
#
# Append to a group and nothing else moves; inserting mid-group still renumbers
# that group's tail, so CI pins (see .github/workflows/test_and_deploy.yml) name
# cases at the head of a group. Blank lines separate the groups and are skipped
# by the runner; a '#' comment inside the block would NOT be, so keep
# annotations up here.
# ---------------------------------------------------------------------------
read -r -d '' CASES <<'EOF'
test_core_1|core|rows|PISA interface pair (PDB vs PDB)|4Q5J:A_E 4Q5J:B_F
test_core_2|core|rows|Batch file vs batch file, PISA interfaces both sides|pdb_pisa_in.txt pdb_pisa_in.txt
test_core_3|core|rows|Chain-ID case sensitivity (4DX9 a_b vs A_B)|4DX9.txt 4DX9.txt
test_core_4|core|ok|PISA interface vs single-residue AlphaFold pocket (human CDK2)|4Q5J:B_F P24941:A:160
test_core_5|core|ok|PISA interface vs single-residue AlphaFold pocket (mouse ortholog)|4Q5J:B_F P97377:A:160
test_core_6|core|rows|AlphaFold passthrough vs AlphaFold passthrough|P06493:A:160,161,162,163,164,165 P24941:A:160,161,162,163,164,165
test_core_7|core|rows|Two pockets on one query chain (pisa + passthrough)|multi_pocket_chain.txt 4Q5J:B_F
test_core_8|core|rows|Superposing on the pocket rather than the chain, with foldseek|4Q5J:A_E 4Q5J:B_F --align_struct_method pocket
test_core_9|core|queries=2|Same-named local files in different directories|same_name.txt same_name.txt
test_core_10|core|rows|Microheterogeneous residue (4Z0Y:A 252 is HS8 and HIS)|4Z0Y:A_E 4Z0Y:C_G

test_open_1|core|rows|PISA interface vs an open whole-chain target|4Q5J:A_E 4Q5J:B
test_open_2|core|rows|PISA interface vs a bare structure, chain defaulting to A|4Q5J:B_F 4Q5J

test_domains_1|human_domains|rows|Single PISA interface vs human domains|4Q5J:B_F human_domains
test_domains_2|human_domains|rows|Mixed batch file (PDB, local mmCIF, AlphaFold) vs human domains|testfile.txt human_domains
test_domains_3|human_domains|rows|AlphaFold passthrough residues vs human domains|P06493:A:160,161,162,163,164,165 human_domains
test_domains_4|human_domains|rows|Large CDK2 pocket residue list vs human domains|1B38:A:8,9,10,11,12,13,14,15,16,17,18,19,20,30,31,32,33,34,35,47,48,49,50,51,52,53,54,55,56,57,58,59,61,62,63,64,65,66,67,68,69,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,143,144,145,146,147,148,149 human_domains
test_domains_5|human_domains|rows|Large kinase pocket residue list vs human domains|4WB5:A:47,48,49,50,51,52,53,54,55,56,57,58,59,69,70,71,72,73,74,87,88,89,90,91,92,93,94,95,96,97,98,99,101,102,103,104,105,106,107,108,109,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,182,183,184,185,186,187,188 human_domains

test_fsdb_1|needs-pdb-fsdb slow|rows|PISA interface vs a local Foldseek PDB database|4Q5J:B_F @PDB_FSDB@ --target_pocket_method foldseek_db
test_fsdb_2|needs-pdb-download huge|rows|PISA interface vs the bundled full-PDB Foldseek database|4Q5J:A_E pdb

test_local_1|core local|rows|Local BLOSUM62 sequence alignment, no Foldseek (same pair as test_core_1)|4Q5J:A_E 4Q5J:B_F --aligner seq
test_local_2|local|rows|Local aligner over mixed input types (PDB, local mmCIF, AlphaFold)|testfile.txt testfile.txt --aligner seq
test_local_3|core local|rows|Open whole-chain target on the local aligner|4Q5J:A_E 4Q5J:B --aligner seq
test_local_4|core local|rows|Explicit pocket superposition on the local aligner|4Q5J:A_E 4Q5J:B_F --aligner seq --align_struct_method pocket
test_local_5|core local|fail|align_struct_method foldseek rejected on the local aligner|4Q5J:A_E 4Q5J:B_F --aligner seq --align_struct_method foldseek
test_local_6|core local|fail|Unknown align_struct_method rejected|4Q5J:A_E 4Q5J:B_F --aligner seq --align_struct_method bogus
test_local_7|core local|queries=2|Same-named local files on the local aligner|same_name.txt same_name.txt --aligner seq

test_invalid_1|core|rows failed=pocket_not_built@4Q5J:A:9999|Passthrough residue id absent from the chain is skipped|invalid_residues.txt 4Q5J:B_F --aligner seq
test_invalid_2|core|rows|Duplicated passthrough residue ids collapsed|4Q5J:A:1101,1101,1104 4Q5J:B_F --aligner seq
test_invalid_3|core|fail|Forced passthrough with no residue list rejected|4Q5J:A 4Q5J:B_F --aligner seq --query_pocket_method passthrough
test_invalid_4|core|fail|Unknown forced pocket method rejected|4Q5J:A_E 4Q5J:B_F --aligner seq --query_pocket_method psia
test_invalid_5|core|fail|Forced pisa with no partner chain rejected|4Q5J:A 4Q5J:B_F --aligner seq --query_pocket_method pisa
test_invalid_6|core|fail|Forced pisa on an AlphaFold entry rejected|P24941:A_B 4Q5J:B_F --aligner seq --query_pocket_method pisa
test_invalid_7|core|fail|Forced vdw with no partner chain rejected|4Q5J:A 4Q5J:B_F --aligner seq --query_pocket_method vdw
test_invalid_8|core|rows failed=invalid_entry@4Q5J:A|Forced method skips only the entries that cannot use it|forced_pisa_mixed.txt 4Q5J:B_F --aligner seq --query_pocket_method pisa
test_invalid_9|core|fail|Unknown aligner rejected|4Q5J:A_E 4Q5J:B_F --aligner bogus

test_settings_1|core settings|rows|Path options set on the command line|4Q5J:A_E 4Q5J:B_F --aligner seq --pdb_dir @CACHE@/pdb_structures --alphafold_dir @CACHE@/alphafold_structures --pocket_dir @CACHE@/pockets --alignment_path @OUT@/custom_alignment.tsv --aligned_structure_dir @OUT@/custom_aligned --job_settings_path @OUT@/custom_settings.json --log_path @OUT@/custom.log --temp_dir @OUT@/custom_temp
test_settings_2|core settings|rows|CLI arguments win over the job file|--job_file job_file.json 4Q5J:A_E 4Q5J:B_F --aligner seq --align_count 5
test_settings_3|core settings|rows|Temp directories kept with --delete_tmp 0|4Q5J:A_E 4Q5J:B_F --aligner seq --delete_tmp 0
test_settings_4|core settings|rows|Explicit --threads accepted and honoured|4Q5J:A_E 4Q5J:B_F --aligner seq --threads 2
test_settings_5|core settings|rows|Job file supplies query and target|--job_file job_file_qt.json --aligner seq
test_settings_6|core settings|fail|Query given both positionally and in the job file|--job_file job_file_qt.json 4Q5J:A_E 4Q5J:B_F --aligner seq
test_settings_7|core settings|rows|Explicit auto pocket method infers as the default does|4Q5J:A_E 4Q5J:B_F --aligner seq --query_pocket_method auto --target_pocket_method auto
test_settings_8|core settings|fail|delete_tmp other than 1 or 0 rejected|4Q5J:A_E 4Q5J:B_F --aligner seq --delete_tmp 2

test_steps_1|core local|files=query_records.json,target_records.json,cache_dirs.json,failed_entries.json,parse_settings.json|parse writes the records files and the cache manifest|parse 4Q5J:A_E 4Q5J:B_F
test_steps_2|core local|files=query_records.json,target_records.json|fetch_structures after parse|parse 4Q5J:A_E 4Q5J:B_F ; fetch_structures
test_steps_3|core local|files=alignment.tsv|align after fetch_structures|parse 4Q5J:A_E 4Q5J:B_F ; fetch_structures ; align --aligner seq
test_steps_4|core local|files=pockets.json|pockets after align|parse 4Q5J:A_E 4Q5J:B_F ; fetch_structures ; align --aligner seq ; pockets @OUT@/query_records.json @OUT@/target_records.json
test_steps_5|core local|rows|compare after pockets|parse 4Q5J:A_E 4Q5J:B_F ; fetch_structures ; align --aligner seq ; pockets @OUT@/query_records.json @OUT@/target_records.json ; compare
test_steps_6|core local|rows files=aligned_structures|superpose after align --aligner seq, the aligner not restated|parse 4Q5J:A_E 4Q5J:B_F ; fetch_structures ; align --aligner seq ; pockets @OUT@/query_records.json @OUT@/target_records.json ; compare ; superpose
test_steps_7|core|rows same=search|Chained commands give what search gives (PISA batch vs itself)|search pdb_pisa_in.txt pdb_pisa_in.txt --results_dir @OUT@/search ; parse pdb_pisa_in.txt pdb_pisa_in.txt ; fetch_structures ; align ; pockets @OUT@/query_records.json @OUT@/target_records.json ; compare ; superpose
test_steps_8|human_domains|rows same=search|Chained commands give what search gives (mixed batch vs human domains), run from the results dir after parse|search testfile.txt human_domains --results_dir @OUT@/search ; parse testfile.txt human_domains ; cd @OUT@ ; fetch_structures ; align ; pockets @OUT@/query_records.json @OUT@/target_records.json ; compare ; superpose
test_steps_9|needs-pdb-fsdb slow|rows same=rerun|align rerun on its own output against a PDB Foldseek database gives the same rows|parse 4Q5J:B_F @PDB_FSDB@ --target_pocket_method foldseek_db ; fetch_structures ; align ; pockets @OUT@/query_records.json @OUT@/target_records.json ; compare ; align --query_records_path @OUT@/rerun/query_records.json --target_records_path @OUT@/rerun/target_records.json --alignment_path @OUT@/rerun/alignment.tsv ; pockets @OUT@/rerun/query_records.json @OUT@/rerun/target_records.json --alignment @OUT@/rerun/alignment.tsv --pockets_path @OUT@/rerun/pockets.json ; compare --target_records @OUT@/rerun/target_records.json --alignment @OUT@/rerun/alignment.tsv --pockets @OUT@/rerun/pockets.json --pocket_comparison_path @OUT@/rerun/pocket_comparison.tsv
test_steps_10|core local|fail|compare rejects an alignment naming chains pockets never saw|parse 4Q5J:A_E 4Q5J:B_F ; fetch_structures ; align --aligner seq ; pockets @OUT@/query_records.json ; compare
test_steps_11|core local|fail|parse rejects a Foldseek database beside a structure target|parse 4Q5J:A_E fsdb_mixed_target.txt
test_steps_12|core local|fail failed=invalid_entry@human_domains|parse rejects a Foldseek database query|parse human_domains 4Q5J:B_F
test_steps_14|core local|rows|parse takes query and target from a job file; pockets defaults to both records files|parse --job_file job_file_qt.json ; fetch_structures ; align --aligner seq ; pockets ; compare
test_steps_15|core local|rows files=compare_settings.json|A step reads its inputs from another run's job_settings.json; an argument overrides it|search 4Q5J:A_E 4Q5J:B_F --aligner seq --align_count 0 --results_dir @OUT@/search ; compare --job_file @OUT@/search/job_settings.json --pocket_comparison_path @OUT@/pocket_comparison.tsv
EOF

# ---------------------------------------------------------------------------
# Defaults (all overridable)
# ---------------------------------------------------------------------------
OUT_DIR="${POCKETMAPPER_E2E_OUT:-$PWD/e2e_results}"
CACHE_DIR="${POCKETMAPPER_E2E_CACHE:-}"
PDB_FSDB="${POCKETMAPPER_PDB_FSDB:-}"
VERBOSITY="${POCKETMAPPER_E2E_VERBOSITY:-4}"
POCKETMAPPER_BIN="${POCKETMAPPER_BIN:-pocketmapper}"
KEEP=0
LIST=0
DRY_RUN=0
TAG_FILTER=""
SELECTED=""

usage() {
    cat <<USAGE
End-to-end test cases for PocketMapper.

Usage: $(basename "$0") [OPTIONS] [TEST_NAME...]

With no TEST_NAME, runs every case except those tagged 'huge' or whose
required resources are unavailable (those are reported as SKIP).

Options:
  -o, --out-dir DIR     Where each case writes its results, one subdirectory
                        per case. Default: \$PWD/e2e_results
                        (env: POCKETMAPPER_E2E_OUT)
  -c, --cache-dir DIR   Shared PocketMapper cache. Reused across cases and
                        across runs, so a warm cache makes reruns much faster.
                        Default: <out-dir>/pocketmapper_cache
                        (env: POCKETMAPPER_E2E_CACHE)
  -t, --tag TAG         Only run cases carrying TAG (e.g. core, slow,
                        human_domains).
  -k, --keep            Keep any existing results instead of clearing each
                        case's directory before it runs.
  -n, --dry-run         Print the commands that would run, then exit.
  -l, --list            List the available cases and exit.
  -v, --verbosity N     PocketMapper verbosity (4=DEBUG .. 1=ERROR). Default: 4
  -h, --help            Show this message.

Environment:
  POCKETMAPPER_BIN        pocketmapper executable to test. Default: pocketmapper
  POCKETMAPPER_PDB_FSDB   Path to a prebuilt Foldseek PDB database. Required
                          for test_fsdb_1, which is skipped when unset.

Notes:
  * Foldseek is the CLI default, so most cases need the 'foldseek' binary on
    PATH and are skipped without it. Cases tagged 'local' pass
    '--aligner seq' to use the built-in BLOSUM62 aligner and still run.
  * Cases hit wwPDB, AlphaFold and PDBe PISA, so they need network access.
  * test_fsdb_2 downloads the full PDB Foldseek database (2GB download, 7Gb unzipped) and is
    therefore excluded unless named explicitly.

Examples:
  $(basename "$0") -o /tmp/pm_e2e             # everything, into /tmp/pm_e2e
  $(basename "$0") -t core                    # quick cases only
  $(basename "$0") test_core_1 test_fsdb_1   # two specific cases
USAGE
}

while [ $# -gt 0 ]; do
    case "$1" in
        -o|--out-dir)    OUT_DIR="$2"; shift 2 ;;
        -c|--cache-dir)  CACHE_DIR="$2"; shift 2 ;;
        -t|--tag)        TAG_FILTER="$2"; shift 2 ;;
        -v|--verbosity)  VERBOSITY="$2"; shift 2 ;;
        -k|--keep)       KEEP=1; shift ;;
        -n|--dry-run)    DRY_RUN=1; shift ;;
        -l|--list)       LIST=1; shift ;;
        -h|--help)       usage; exit 0 ;;
        -*)              echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
        *)               SELECTED="$SELECTED $1"; shift ;;
    esac
done

if [ "$LIST" -eq 1 ]; then
    printf '%-16s %-26s %s\n' "NAME" "TAGS" "DESCRIPTION"
    while IFS='|' read -r name tags expect desc args; do
        [ -z "$name" ] && continue
        printf '%-16s %-26s %s\n' "$name" "$tags" "$desc"
    done <<< "$CASES"
    exit 0
fi

# Resolve to an absolute path: cases run with the fixtures directory as their
# working directory (testfile.txt refers to 4Q5J.cif.gz relatively), so a
# relative --out-dir would otherwise land inside the repo.
mkdir -p "$OUT_DIR" || { echo "Cannot create out-dir: $OUT_DIR" >&2; exit 2; }
OUT_DIR="$(cd "$OUT_DIR" && pwd)"
[ -z "$CACHE_DIR" ] && CACHE_DIR="$OUT_DIR/pocketmapper_cache"
mkdir -p "$CACHE_DIR" || { echo "Cannot create cache-dir: $CACHE_DIR" >&2; exit 2; }
CACHE_DIR="$(cd "$CACHE_DIR" && pwd)"

# ---------------------------------------------------------------------------
# Preflight
# ---------------------------------------------------------------------------
if ! command -v "$POCKETMAPPER_BIN" >/dev/null 2>&1; then
    echo "ERROR: '$POCKETMAPPER_BIN' not found on PATH. Install with 'pip install -e .'" >&2
    exit 2
fi
HAVE_FOLDSEEK=1
if ! command -v foldseek >/dev/null 2>&1; then
    HAVE_FOLDSEEK=0
    echo "WARNING: 'foldseek' not found on PATH; only the '--aligner seq' cases will run." >&2
fi

echo "pocketmapper : $(command -v "$POCKETMAPPER_BIN")"
echo "out-dir      : $OUT_DIR"
echo "cache-dir    : $CACHE_DIR"
echo

PASS=0; FAIL=0; SKIP=0
FAILED_NAMES=""

selected() {
    [ -z "$SELECTED" ] && return 1
    for s in $SELECTED; do [ "$s" = "$1" ] && return 0; done
    return 1
}

has_tag() {
    for t in $1; do [ "$t" = "$2" ] && return 0; done
    return 1
}

COMMANDS="parse fetch_structures align pockets compare superpose search"

is_command() {
    for c in $COMMANDS; do [ "$c" = "$1" ] && return 0; done
    return 1
}

# The command each ';'-separated segment of a case's args runs; `search` for a
# case that names none.
segment_commands() {
    local first=1 word
    is_command "${1%% *}" || { echo search; return; }
    # shellcheck disable=SC2086  # deliberate word-splitting of the args field
    for word in $1; do
        if [ "$word" = ";" ]; then first=1; continue; fi
        [ "$first" -eq 1 ] && echo "$word"
        first=0
    done
}

while IFS='|' read -r name tags expect desc args; do
    [ -z "$name" ] && continue

    explicit=0
    if [ -n "$SELECTED" ]; then
        selected "$name" || continue
        explicit=1
    fi
    if [ -n "$TAG_FILTER" ] && ! has_tag "$tags" "$TAG_FILTER"; then
        continue
    fi

    # --- gating -----------------------------------------------------------
    # Foldseek is needed by a search or align segment, unless the case uses the
    # local aligner.
    skip_reason=""
    uses_foldseek=0
    for cmd in $(segment_commands "$args"); do
        case "$cmd" in search|align) uses_foldseek=1 ;; esac
    done
    case "$args" in
        *"--aligner seq"*) uses_foldseek=0 ;;
    esac
    if [ "$uses_foldseek" -eq 1 ] && [ "$HAVE_FOLDSEEK" -eq 0 ]; then
        skip_reason="foldseek not installed"
    elif has_tag "$tags" "needs-pdb-fsdb" && [ -z "$PDB_FSDB" ]; then
        skip_reason="POCKETMAPPER_PDB_FSDB not set"
    elif has_tag "$tags" "needs-pdb-fsdb" && [ ! -e "$PDB_FSDB" ]; then
        skip_reason="POCKETMAPPER_PDB_FSDB=$PDB_FSDB does not exist"
    elif has_tag "$tags" "huge" && [ "$explicit" -eq 0 ]; then
        skip_reason="tagged 'huge'; name it explicitly to run"
    fi

    if [ -n "$skip_reason" ]; then
        printf 'SKIP  %-16s %s (%s)\n' "$name" "$desc" "$skip_reason"
        SKIP=$((SKIP + 1))
        continue
    fi

    # --- build the commands -----------------------------------------------
    case_out="$OUT_DIR/$name"
    resolved_args="${args//@PDB_FSDB@/$PDB_FSDB}"
    resolved_args="${resolved_args//@CACHE@/$CACHE_DIR}"
    resolved_args="${resolved_args//@OUT@/$case_out}"
    # A case that names no command is one search
    is_command "${args%% *}" || resolved_args="search $resolved_args"

    if [ "$KEEP" -eq 0 ] && [ -d "$case_out" ] && [ "$DRY_RUN" -eq 0 ]; then
        # Scoped to this case's own directory; never touches OUT_DIR itself.
        rm -rf "$case_out"
    fi

    # One line per segment: the working directory, a tab, then the command.
    # `cd DIR` segments only move the working directory of those after them.
    segments=""
    cwd="$FIXTURES_DIR"
    segment=""
    # shellcheck disable=SC2086  # deliberate word-splitting of the args field
    for word in $resolved_args ";"; do
        if [ "$word" != ";" ]; then
            segment="$segment $word"
            continue
        fi
        set -- $segment
        segment=""
        [ $# -eq 0 ] && continue
        if [ "$1" = "cd" ]; then
            cwd="$2"
            continue
        fi
        extra="--verbosity $VERBOSITY"
        case " $* " in *" --results_dir "*) ;; *) extra="$extra --results_dir $case_out" ;; esac
        case "$1" in
            search|parse)
                case " $* " in *" --cache_dir "*) ;; *) extra="$extra --cache_dir $CACHE_DIR" ;; esac ;;
        esac
        segments="$segments$cwd	$* $extra
"
    done

    if [ "$DRY_RUN" -eq 1 ]; then
        while IFS='	' read -r dir cmd; do
            [ -z "$cmd" ] && continue
            printf 'DRY   %-16s (cd %s && %s %s)\n' "$name" "$dir" "$POCKETMAPPER_BIN" "$cmd"
        done <<< "$segments"
        continue
    fi

    printf 'RUN   %-16s %s\n' "$name" "$desc"
    log="$OUT_DIR/$name.log"
    : > "$log"
    started=$(date +%s)
    # Every segment but the last must succeed; the last one's status is the case's.
    status=0
    early_failure=0
    remaining=$(printf '%s' "$segments" | grep -c .)
    while IFS='	' read -r dir cmd; do
        [ -z "$cmd" ] && continue
        remaining=$((remaining - 1))
        echo "### (cd $dir && pocketmapper $cmd)" >> "$log"
        # shellcheck disable=SC2086  # deliberate word-splitting of the command
        ( mkdir -p "$dir" && cd "$dir" && "$POCKETMAPPER_BIN" $cmd ) >> "$log" 2>&1 < /dev/null
        status=$?
        if [ "$status" -ne 0 ]; then
            [ "$remaining" -gt 0 ] && early_failure=1
            break
        fi
    done <<< "$segments"
    elapsed=$(( $(date +%s) - started ))

    # --- assertions -------------------------------------------------------
    comparison="$case_out/pocket_comparison.tsv"
    problem=""
    # The one token about pocket_comparison.tsv or the exit status, if any
    outcome=""
    for token in $expect; do
        case "$token" in
            rows|ok|fail|queries=*) outcome="$token" ;;
        esac
    done
    if [ "$early_failure" -eq 1 ]; then
        problem="a command before the last exited $status"
    elif [ "$outcome" = "fail" ]; then
        # A rejected input or option combination
        [ "$status" -eq 0 ] && problem="exit=0 (expected a rejection)"
    elif [ "$status" -ne 0 ]; then
        problem="exit=$status"
    fi
    rows=-1
    if [ -z "$problem" ] && [ -n "$outcome" ] && [ "$outcome" != "fail" ]; then
        if [ ! -f "$comparison" ]; then
            problem="no pocket_comparison.tsv"
        else
            rows=$(( $(wc -l < "$comparison") - 1 ))
            [ "$rows" -lt 0 ] && rows=0
            # Distinct values of the column headed `query`
            queries=$(awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == "query") c = i; next }
                                   c { print $c }' "$comparison" | sort -u | wc -l | tr -d ' ')
            if [ "$outcome" = "rows" ] && [ "$rows" -lt 1 ]; then
                problem="0 comparison rows (expected >=1)"
            elif [ "${outcome#queries=}" != "$outcome" ] && [ "$queries" -lt "${outcome#queries=}" ]; then
                problem="$queries distinct queries (expected >=${outcome#queries=})"
            fi
        fi
    fi
    if [ -z "$problem" ]; then
        for token in $expect; do
            case "$token" in
                files=*)
                    for f in $(echo "${token#files=}" | tr ',' ' '); do
                        [ -e "$case_out/$f" ] || { problem="no $f"; break; }
                    done ;;
                failed=*)
                    reason="${token#failed=}"; pocket_id="${reason#*@}"; reason="${reason%%@*}"
                    python3 - "$case_out/failed_entries.json" "$reason" "$pocket_id" <<'PY' || problem="failed_entries.json has no $reason entry for $pocket_id"
import json, sys
path, reason, pocket_id = sys.argv[1:]
entries = json.load(open(path))
sys.exit(0 if any(e["reason"] == reason and e["pocket_id"] == pocket_id for e in entries) else 1)
PY
                    ;;
                same=*)
                    other="$case_out/${token#same=}"
                    for f in pocket_comparison.tsv alignment.tsv; do
                        if [ ! -f "$case_out/$f" ] || [ ! -f "$other/$f" ]; then
                            problem="$f missing from $case_out or $other"; break
                        fi
                        # Row order is not stable between runs
                        if ! cmp -s <(sort "$case_out/$f") <(sort "$other/$f"); then
                            problem="$f differs from $other/$f"; break
                        fi
                    done ;;
            esac
            [ -n "$problem" ] && break
        done
    fi

    if [ -n "$problem" ]; then
        printf '  FAIL  %s after %ds -- see %s\n' "$problem" "$elapsed" "$log"
        FAIL=$((FAIL + 1)); FAILED_NAMES="$FAILED_NAMES $name"
    elif [ "$outcome" = "fail" ]; then
        printf '  PASS  rejected as expected in %ds\n' "$elapsed"
        PASS=$((PASS + 1))
    elif [ "$rows" -ge 0 ]; then
        printf '  PASS  %d comparison rows in %ds\n' "$rows" "$elapsed"
        PASS=$((PASS + 1))
    else
        printf '  PASS  in %ds\n' "$elapsed"
        PASS=$((PASS + 1))
    fi
done <<< "$CASES"

[ "$DRY_RUN" -eq 1 ] && exit 0

echo
echo "passed: $PASS  failed: $FAIL  skipped: $SKIP"
if [ "$FAIL" -gt 0 ]; then
    echo "failed cases:$FAILED_NAMES"
    exit 1
fi
exit 0
