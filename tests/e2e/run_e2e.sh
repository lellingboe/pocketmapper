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
# expect   `rows` = run must succeed AND pocket_comparison.tsv must hold at
#          least one data row; `ok` = run must succeed and write the file, but
#          zero hits is a legitimate outcome for that pair.
# args     passed to `pocketmapper search` verbatim (word-split on spaces).
#          @PDB_FSDB@ expands to $POCKETMAPPER_PDB_FSDB. Foldseek is the
#          CLI's default, so a case is assumed to need the binary and is
#          skipped when it is missing; a case opts out with an explicit
#          `--foldseek False`, which exercises the local BLOSUM62 aligner
#          and still runs without the binary.
#
# Cases are grouped by what they exercise -- structure-vs-structure pairs, then
# human_domains DB targets, then the larger Foldseek DB targets, then the local
# aligner -- and numbered in that order. Blank lines separate the groups and are
# skipped by the runner; a '#' comment inside the block would NOT be, so keep
# annotations up here.
# ---------------------------------------------------------------------------
read -r -d '' CASES <<'EOF'
test_1|core|rows|PISA interface pair (PDB vs PDB)|4Q5J:A_E 4Q5J:B_F --foldseek
test_2|core|rows|Batch file vs batch file, PISA interfaces both sides|pdb_pisa_in.txt pdb_pisa_in.txt --foldseek
test_3|core|rows|Chain-ID case sensitivity (4DX9 a_b vs A_B)|4DX9.txt 4DX9.txt --foldseek
test_4|core|ok|PISA interface vs single-residue AlphaFold pocket (human CDK2)|4Q5J:B_F P24941:A:160 --foldseek
test_5|core|ok|PISA interface vs single-residue AlphaFold pocket (mouse ortholog)|4Q5J:B_F P97377:A:160 --foldseek
test_6|core|rows|AlphaFold passthrough vs AlphaFold passthrough|P06493:A:160,161,162,163,164,165 P24941:A:160,161,162,163,164,165 --foldseek
test_7|core|rows|Two pockets on one query chain (pisa + passthrough)|multi_pocket_chain.txt 4Q5J:B_F --foldseek

test_8|human_domains|rows|Single PISA interface vs human domains|4Q5J:B_F human_domains --foldseek
test_9|human_domains|rows|Mixed batch file (PDB, local mmCIF, AlphaFold) vs human domains|testfile.txt human_domains --foldseek
test_10|human_domains|rows|AlphaFold passthrough residues vs human domains|P06493:A:160,161,162,163,164,165 human_domains --foldseek
test_11|human_domains|rows|Large CDK2 pocket residue list vs human domains|--query 1B38:A:8,9,10,11,12,13,14,15,16,17,18,19,20,30,31,32,33,34,35,47,48,49,50,51,52,53,54,55,56,57,58,59,61,62,63,64,65,66,67,68,69,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,143,144,145,146,147,148,149 --target human_domains --foldseek
test_12|human_domains|rows|Large kinase pocket residue list vs human domains|--query 4WB5:A:47,48,49,50,51,52,53,54,55,56,57,58,59,69,70,71,72,73,74,87,88,89,90,91,92,93,94,95,96,97,98,99,101,102,103,104,105,106,107,108,109,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,182,183,184,185,186,187,188 --target human_domains --foldseek

test_13|needs-pdb-fsdb slow|rows|PISA interface vs a local Foldseek PDB database|4Q5J:B_F @PDB_FSDB@ --target_pocket_method foldseek_db --foldseek
test_14|needs-pdb-download huge|rows|PISA interface vs the bundled full-PDB Foldseek database|4Q5J:A_E pdb --foldseek

test_15|core local|rows|Local BLOSUM62 sequence alignment, no Foldseek (same pair as test_1)|4Q5J:A_E 4Q5J:B_F --foldseek False
test_16|local|rows|Local aligner over mixed input types (PDB, local mmCIF, AlphaFold)|testfile.txt testfile.txt --foldseek False
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
                          for test_13, which is skipped when unset.

Notes:
  * Foldseek is the CLI default, so most cases need the 'foldseek' binary on
    PATH and are skipped without it. Cases tagged 'local' pass
    '--foldseek False' to force the built-in BLOSUM62 aligner and still run.
  * Cases hit wwPDB, AlphaFold and PDBe PISA, so they need network access.
  * test_14 downloads the full PDB Foldseek database (2GB download, 7Gb unzipped) and is
    therefore excluded unless named explicitly.

Examples:
  $(basename "$0") -o /tmp/pm_e2e             # everything, into /tmp/pm_e2e
  $(basename "$0") -t core                    # quick cases only
  $(basename "$0") test_1 test_15             # two specific cases
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
    printf '%-9s %-26s %s\n' "NAME" "TAGS" "DESCRIPTION"
    while IFS='|' read -r name tags expect desc args; do
        [ -z "$name" ] && continue
        printf '%-9s %-26s %s\n' "$name" "$tags" "$desc"
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
    echo "WARNING: 'foldseek' not found on PATH; only the '--foldseek False' cases will run." >&2
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
    skip_reason=""
    case "$args" in
        *"--foldseek False"*) uses_foldseek=0 ;;
        *)                    uses_foldseek=1 ;;
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
        printf 'SKIP  %-9s %s (%s)\n' "$name" "$desc" "$skip_reason"
        SKIP=$((SKIP + 1))
        continue
    fi

    # --- build the command ------------------------------------------------
    case_out="$OUT_DIR/$name"
    resolved_args="${args//@PDB_FSDB@/$PDB_FSDB}"

    if [ "$KEEP" -eq 0 ] && [ -d "$case_out" ]; then
        # Scoped to this case's own directory; never touches OUT_DIR itself.
        rm -rf "$case_out"
    fi

    # shellcheck disable=SC2086  # deliberate word-splitting of the args field
    set -- $resolved_args \
        --verbosity "$VERBOSITY" \
        --cache_dir "$CACHE_DIR" \
        --results_dir "$case_out"

    if [ "$DRY_RUN" -eq 1 ]; then
        printf 'DRY   %-9s (cd %s && %s search %s)\n' "$name" "$FIXTURES_DIR" "$POCKETMAPPER_BIN" "$*"
        continue
    fi

    printf 'RUN   %-9s %s\n' "$name" "$desc"
    log="$OUT_DIR/$name.log"
    started=$(date +%s)
    ( cd "$FIXTURES_DIR" && "$POCKETMAPPER_BIN" search "$@" ) > "$log" 2>&1
    status=$?
    elapsed=$(( $(date +%s) - started ))

    # --- assertions -------------------------------------------------------
    comparison="$case_out/pocket_comparison.tsv"
    if [ "$status" -ne 0 ]; then
        printf '  FAIL  exit=%d after %ds -- see %s\n' "$status" "$elapsed" "$log"
        FAIL=$((FAIL + 1)); FAILED_NAMES="$FAILED_NAMES $name"
    elif [ ! -f "$comparison" ]; then
        printf '  FAIL  no pocket_comparison.tsv after %ds -- see %s\n' "$elapsed" "$log"
        FAIL=$((FAIL + 1)); FAILED_NAMES="$FAILED_NAMES $name"
    else
        rows=$(( $(wc -l < "$comparison") - 1 ))
        [ "$rows" -lt 0 ] && rows=0
        if [ "$expect" = "rows" ] && [ "$rows" -lt 1 ]; then
            printf '  FAIL  0 comparison rows after %ds (expected >=1) -- see %s\n' "$elapsed" "$log"
            FAIL=$((FAIL + 1)); FAILED_NAMES="$FAILED_NAMES $name"
        else
            printf '  PASS  %d comparison rows in %ds\n' "$rows" "$elapsed"
            PASS=$((PASS + 1))
        fi
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
