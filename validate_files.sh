#!/bin/bash
# ============================================================================
# Validate input ROOT files for a given .list and tree name.
#
# Opens each file, gets the requested tree, and reads through ALL its entries
# in the small subset of branches the analysis touches. Files that crash, error,
# or return zero entries are flagged. Produces:
#   <list>_clean.list  — only files that fully passed
#   <list>_bad.list    — files that failed (with reason)
#
# Each file is checked in its OWN process so a SEGV in one doesn't take down
# the validator. Runs N at a time in parallel.
#
# Usage:
#   ./validate_files.sh                                   # reads from config.json
#   ./validate_files.sh <list_file> <tree_name>           # explicit
#   ./validate_files.sh <list_file> <tree_name> [N=16]    # custom concurrency
# ============================================================================

set -uo pipefail

LIST="${1:-}"
TREE="${2:-}"
N_CONCURRENT="${3:-16}"

# If no args: pull source list and tree name from config.json
if [[ -z "$LIST" || -z "$TREE" ]]; then
    CONFIG="config.json"
    [[ -f "$CONFIG" ]] || {
        echo "ERROR: no args given and config.json not found in CWD."
        echo "Usage: $0 [<list_file> <tree_name> [N_concurrent=16]]"
        exit 1
    }
    command -v python3 >/dev/null || { echo "ERROR: python3 needed to parse $CONFIG"; exit 1; }
    LIST=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['input']['source'])" "$CONFIG")
    TREE=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['input']['tree_name'])" "$CONFIG")
    echo "(auto-detected from $CONFIG: list=$LIST, tree=$TREE)"
fi

[[ -f "$LIST" ]] || { echo "ERROR: list not found: $LIST"; exit 1; }
[[ "$LIST" == *.list ]] || { echo "ERROR: source must be a .list file (got: $LIST)"; exit 1; }
command -v root >/dev/null || { echo "ERROR: root not in PATH"; exit 1; }

CLEAN="${LIST%.list}_clean.list"
BAD="${LIST%.list}_bad.list"
TMP=$(mktemp -d -t fat_validate.XXXXXX)
trap 'rm -rf "$TMP"' EXIT

# ---------- per-file checker (run in subshell so SEGV is contained) ----------
check_one() {
    local path="$1"
    local tree="$2"
    local stamp="$3"   # writes $stamp.ok or $stamp.bad on completion

    # ROOT in batch mode, exit code 0 = ok, anything else = bad
    root -l -b -q -e "
{
    auto f = TFile::Open(\"$path\");
    if (!f || f->IsZombie()) { gApplication->Terminate(11); }
    auto t = (TTree*)f->Get(\"$tree\");
    if (!t) { gApplication->Terminate(12); }
    Long64_t n = t->GetEntries();
    if (n <= 0) { gApplication->Terminate(13); }
    // Read every entry — this is the only way to catch corrupt baskets.
    for (Long64_t i = 0; i < n; ++i) {
        if (t->GetEntry(i) <= 0) { gApplication->Terminate(14); }
    }
    f->Close();
    gApplication->Terminate(0);
}
" >/dev/null 2>&1
    local rc=$?
    if (( rc == 0 )); then
        echo "$path" > "$stamp.ok"
    else
        echo "$path  rc=$rc" > "$stamp.bad"
    fi
}
export -f check_one

# ---------- enumerate files (handle quoted-comma-separated format) ----------
mapfile -t LINES < <(grep -vE '^\s*(#|//|$)' "$LIST")
TOTAL=${#LINES[@]}
echo "Validating $TOTAL files from $LIST against tree '$TREE'..."
echo "(parallel: $N_CONCURRENT)"

idx=0
for line in "${LINES[@]}"; do
    # Strip quotes/commas/whitespace
    path=$(echo "$line" | sed -E 's/^[[:space:]]*"?//;s/"?[[:space:],;]*$//')
    [[ -n "$path" ]] || continue

    while (( $(jobs -rp | wc -l) >= N_CONCURRENT )); do
        wait -n 2>/dev/null || true
    done
    stamp="$TMP/$(printf "%05d" $idx)"
    echo "$line" > "$stamp.line"
    check_one "$path" "$TREE" "$stamp" &
    ((idx++))
done
wait

# ---------- collate results, preserving original line format ----------
> "$CLEAN"
> "$BAD"
n_ok=0; n_bad=0
for stamp_line in "$TMP"/*.line; do
    base="${stamp_line%.line}"
    line=$(cat "$stamp_line")
    if [[ -f "$base.ok" ]]; then
        echo "$line" >> "$CLEAN"
        ((n_ok++))
    else
        reason=$(cat "$base.bad" 2>/dev/null || echo "(no rc)")
        echo "$line   # $reason" >> "$BAD"
        ((n_bad++))
    fi
done

echo
echo "================================================================"
echo "Validated $TOTAL files: $n_ok OK, $n_bad BAD"
echo "  Clean list:  $CLEAN"
echo "  Bad list:    $BAD"
echo "================================================================"

if (( n_bad > 0 )); then
    echo
    echo "Bad files:"
    cat "$BAD"
fi
