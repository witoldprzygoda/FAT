#!/bin/bash
# ============================================================================
# trigger_scan parallel launcher — split input list, run N trigger_scan jobs,
# hadd outputs.
#
# Mirrors the structure of the main run_parallel.sh. Reads a flat scanner
# config (input_source, tree_name, trig_pt3, trig_pt2, output_file), splits
# input_source into N chunks, generates a temporary chunk-config per chunk
# pointing at the chunk-list and a per-chunk output ROOT, runs trigger_scan
# on each in parallel (throttled), then hadd-merges the per-chunk outputs
# into the original output_file.
#
# Does NOT modify the project sources, the input config, or the input list.
# All temp files live in $TMPDIR (printed below); cleaned up on success only.
#
# Usage (from research/calibration/):
#   ./run_parallel_scan.sh [config_epem.json] [N_parts=32] [N_concurrent=N_parts]
#
# Examples:
#   ./run_parallel_scan.sh                            # config_epem.json default
#   ./run_parallel_scan.sh config_epep.json
#   ./run_parallel_scan.sh config_emem.json 32 16
# ============================================================================

set -uo pipefail

CONFIG="${1:-config_epem.json}"
N_PARTS="${2:-32}"
N_CONCURRENT="${3:-$N_PARTS}"
N_RETRIES="${4:-2}"
LAUNCH_STAGGER="${5:-1}"

# ---------- sanity ----------
[[ -f "$CONFIG"        ]] || { echo "ERROR: config not found: $CONFIG"; exit 1; }
[[ -x ./trigger_scan   ]] || { echo "ERROR: ./trigger_scan missing or not executable (run 'make' first)"; exit 1; }
command -v hadd    >/dev/null || { echo "ERROR: hadd (ROOT) not in PATH"; exit 1; }
command -v python3 >/dev/null || { echo "ERROR: python3 not in PATH"; exit 1; }
command -v split   >/dev/null || { echo "ERROR: split (coreutils) not in PATH"; exit 1; }

# Pull flat config fields. Bail out cleanly on missing keys.
SOURCE=$(python3 -c "import json,sys; d=json.load(open(sys.argv[1])); print(d['input_source'])" "$CONFIG")
OUTPUT=$(python3 -c "import json,sys; d=json.load(open(sys.argv[1])); print(d['output_file'])"  "$CONFIG")
TREE=$(  python3 -c "import json,sys; d=json.load(open(sys.argv[1])); print(d['tree_name'])"    "$CONFIG")

[[ "$SOURCE" == *.list ]] || { echo "ERROR: input_source must be a .list file (got: $SOURCE)"; exit 1; }
[[ -f "$SOURCE"        ]] || { echo "ERROR: input list not found: $SOURCE"; exit 1; }

# Pick work-dir base — same convention as run_parallel.sh (FAT_TMPDIR > /mnt/scratch
# > $TMPDIR > /tmp). Per-chunk outputs are tiny (one TTree entry per file in
# the chunk list), so disk requirements are negligible compared to the main
# launcher.
if [[ -n "${FAT_TMPDIR:-}" ]]; then
    WORK_BASE="$FAT_TMPDIR"
elif [[ -d /mnt/scratch && -w /mnt/scratch ]]; then
    WORK_BASE="/mnt/scratch/fat_tmp"
    mkdir -p "$WORK_BASE" 2>/dev/null
else
    WORK_BASE="${TMPDIR:-/tmp}"
fi
[[ -d "$WORK_BASE" && -w "$WORK_BASE" ]] || {
    echo "ERROR: work-dir base does not exist or is not writable: $WORK_BASE"
    exit 1
}
TMPDIR=$(mktemp -d -p "$WORK_BASE" fat_scan_parallel.XXXXXX)

cat <<EOF
================================================================
trigger_scan parallel launcher
  Config:         $CONFIG
  Input list:     $SOURCE
  Tree:           $TREE
  Final output:   $OUTPUT
  Parts:          $N_PARTS
  Concurrent:     $N_CONCURRENT
  Work dir:       $TMPDIR
================================================================
EOF

# Strip blanks/comments before splitting.
grep -vE '^\s*(#|//|$)' "$SOURCE" > "$TMPDIR/all.list"
TOTAL=$(wc -l < "$TMPDIR/all.list")
echo "Source has $TOTAL non-comment lines."

if (( N_PARTS > TOTAL )); then
    echo "WARNING: requested $N_PARTS parts but only $TOTAL files. Capping to $TOTAL."
    N_PARTS=$TOTAL
fi

split -n "l/$N_PARTS" -d -a 3 --additional-suffix=.list \
      "$TMPDIR/all.list" "$TMPDIR/part_"

n_empty=0
for f in "$TMPDIR"/part_*.list; do
    if [[ ! -s "$f" ]]; then
        rm -f "$f"
        ((n_empty++))
    fi
done
if (( n_empty > 0 )); then
    echo "Dropped $n_empty empty chunk(s) (split-rounding artefact)."
fi

# Generate one chunk-config + output filename per chunk
for f in "$TMPDIR"/part_*.list; do
    idx=$(basename "$f" | sed -E 's/part_([0-9]+)\.list/\1/')
    cfg="$TMPDIR/cfg_$idx.json"
    out="$TMPDIR/output_part_$idx.root"
    python3 - "$CONFIG" "$f" "$out" "$cfg" <<'PYEOF'
import json, sys
src_cfg, src_list, out_root, dst_cfg = sys.argv[1:5]
with open(src_cfg) as fh:
    c = json.load(fh)
c['input_source'] = src_list
c['output_file']  = out_root
with open(dst_cfg, 'w') as fh:
    json.dump(c, fh, indent=2)
PYEOF
done

CFG_COUNT=$(ls "$TMPDIR"/cfg_*.json | wc -l)
echo "Generated $CFG_COUNT chunked configs. Launching..."
echo "(Per-job stdout/stderr in $TMPDIR/cfg_NNN.log)"
echo

# ---------- run jobs with controlled concurrency ----------
START=$SECONDS
FIRST_JOB=1
for cfg in "$TMPDIR"/cfg_*.json; do
    while (( $(jobs -rp | wc -l) >= N_CONCURRENT )); do
        wait -n 2>/dev/null || true
    done
    if (( FIRST_JOB == 0 && LAUNCH_STAGGER > 0 )); then
        sleep "$LAUNCH_STAGGER"
    fi
    FIRST_JOB=0
    (
        base_log="${cfg%.json}.log"
        success=0
        for try in $(seq 0 $N_RETRIES); do
            if (( try == 0 )); then
                log="$base_log"
            else
                sleep $((5 + RANDOM % 10))
                log="${cfg%.json}_retry${try}.log"
                printf '  [retry %d] %s\n' "$try" "$(basename "$cfg")"
            fi
            if ./trigger_scan "$cfg" > "$log" 2>&1; then
                success=1
                break
            fi
        done
        if (( success )); then
            if (( try == 0 )); then
                printf '  [OK]   %s\n' "$(basename "$cfg")"
            else
                printf '  [OK after %d retries] %s\n' "$try" "$(basename "$cfg")"
            fi
        else
            printf '  [FAIL after %d tries] %s  (see %s and *_retry*.log)\n' \
                $((N_RETRIES + 1)) "$(basename "$cfg")" "$base_log"
        fi
    ) &
done
wait

ELAPSED=$((SECONDS - START))
echo
echo "All scan jobs finished in ${ELAPSED}s."

# ---------- verify and merge ----------
shopt -s nullglob
OUTPUTS=("$TMPDIR"/output_part_*.root)
shopt -u nullglob

if (( ${#OUTPUTS[@]} == 0 )); then
    echo "ERROR: no output_part_*.root produced — aborting (tmp kept at $TMPDIR)"
    exit 2
fi
if (( ${#OUTPUTS[@]} != CFG_COUNT )); then
    echo "WARNING: expected $CFG_COUNT outputs, found ${#OUTPUTS[@]} — some jobs failed."
    echo "Inspect logs in $TMPDIR. Continuing with hadd of what we have."
fi

echo "Merging ${#OUTPUTS[@]} files into $OUTPUT ..."
hadd -f "$OUTPUT" "${OUTPUTS[@]}"
HADD_RC=$?

if (( HADD_RC != 0 )); then
    echo "ERROR: hadd failed (rc=$HADD_RC) — tmp kept at $TMPDIR"
    exit 3
fi

if (( ${#OUTPUTS[@]} == CFG_COUNT )); then
    rm -rf "$TMPDIR"
else
    echo
    echo "NOTE: ${#OUTPUTS[@]}/$CFG_COUNT chunks succeeded — tmp kept at $TMPDIR"
    echo "      Failed chunks: inspect cfg_NNN.log and rerun those configs:"
    for cfg in "$TMPDIR"/cfg_*.json; do
        idx=$(basename "$cfg" .json | sed 's/cfg_//')
        out="$TMPDIR/output_part_${idx}.root"
        [[ -f "$out" ]] || echo "        ./trigger_scan $cfg"
    done
fi

echo "================================================================"
echo "Done. Output: $OUTPUT"
echo "================================================================"
