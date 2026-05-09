#!/bin/bash
# ============================================================================
# FAT parallel launcher — split input list, run N ana jobs, hadd outputs.
#
# Reads config_<channel>.json (or any other .json), splits the .list source
# into N chunks, generates a temporary config per chunk with its own output
# filename, launches up to N_CONCURRENT ana jobs in parallel, then merges
# all part outputs into the original output filename via hadd.
#
# Does NOT modify the project sources, the input config, or the input list.
# All temp files live in $TMPDIR (printed below); cleaned up on success only.
#
# Usage:
#   ./run_parallel.sh [config_epem.json] [N_parts=32] [N_concurrent=N_parts]
#
# Examples:
#   ./run_parallel.sh                              # config_epem.json default
#   ./run_parallel.sh config_epep.json 16 8        # 16 parts, 8 concurrent
#   ./run_parallel.sh config_emem_cal.json 64 16   # cal-mode 64 chunks
# ============================================================================

set -uo pipefail

CONFIG="${1:-config_epem.json}"
N_PARTS="${2:-32}"
N_CONCURRENT="${3:-$N_PARTS}"
N_RETRIES="${4:-2}"          # per-chunk retries on failure
LAUNCH_STAGGER="${5:-3}"     # seconds between consecutive job starts (0 = no stagger)

# ---------- sanity ----------
[[ -f "$CONFIG"   ]] || { echo "ERROR: config not found: $CONFIG"; exit 1; }
[[ -x ./ana      ]] || { echo "ERROR: ./ana missing or not executable (run 'make' first)"; exit 1; }
command -v hadd    >/dev/null || { echo "ERROR: hadd (ROOT) not in PATH"; exit 1; }
command -v python3 >/dev/null || { echo "ERROR: python3 not in PATH"; exit 1; }
command -v split   >/dev/null || { echo "ERROR: split (coreutils) not in PATH"; exit 1; }

SOURCE=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['input']['source'])" "$CONFIG")
OUTPUT=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['output']['filename'])" "$CONFIG")

[[ "$SOURCE" == *.list ]] || { echo "ERROR: source must be a .list file (got: $SOURCE)"; exit 1; }
[[ -f "$SOURCE"        ]] || { echo "ERROR: source list not found: $SOURCE"; exit 1; }

# Pick work-dir base.
# Order: FAT_TMPDIR env var > /mnt/scratch/fat_tmp (auto-create if /mnt/scratch
# is a writable mount) > $TMPDIR > /tmp. Each chunk holds ~3-4 intermediate
# .root files of ~600-800 MB while running, so /tmp is usually too small.
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
TMPDIR=$(mktemp -d -p "$WORK_BASE" fat_parallel.XXXXXX)

# Warn if free space looks tight — rough rule: ~1 GB per concurrent job + 0.5 GB per chunk
FREE_GB=$(df -BG --output=avail "$TMPDIR" | tail -1 | tr -d ' G')
NEEDED_GB=$(( N_CONCURRENT * 4 + N_PARTS / 2 ))
if (( FREE_GB < NEEDED_GB )); then
    echo "================================================================"
    echo "WARNING: only ${FREE_GB} GB free in $WORK_BASE — may run out of space."
    echo "         Estimated need: ~${NEEDED_GB} GB (heuristic)."
    echo "         If runs fail with 'No space left on device', use:"
    echo "           FAT_TMPDIR=/path/to/large/disk ./run_parallel.sh ..."
    echo "================================================================"
fi

cat <<EOF
================================================================
FAT parallel launcher
  Config:         $CONFIG
  Source list:    $SOURCE
  Output:         $OUTPUT
  Parts:          $N_PARTS
  Concurrent:     $N_CONCURRENT
  Work dir:       $TMPDIR
================================================================
EOF

# Strip blank lines and comments before splitting (so each part has only real entries)
grep -vE '^\s*(#|//|$)' "$SOURCE" > "$TMPDIR/all.list"
TOTAL=$(wc -l < "$TMPDIR/all.list")
echo "Source has $TOTAL non-comment lines."

if (( N_PARTS > TOTAL )); then
    echo "WARNING: requested $N_PARTS parts but only $TOTAL files. Capping to $TOTAL."
    N_PARTS=$TOTAL
fi

# split -n l/N — N nearly-equal chunks, line-aligned
# Use --additional-suffix=.list so ana recognises each chunk as a file list.
split -n "l/$N_PARTS" -d -a 3 --additional-suffix=.list \
      "$TMPDIR/all.list" "$TMPDIR/part_"

# Drop empty chunks. `split -n l/N` can leave some chunks empty when N is close
# to the line count, due to byte-rounded line boundaries. Empty chunks would
# fail downstream in ana ("no files in list") and pollute the OK/FAIL summary.
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

# Generate one config + output filename per chunk
for f in "$TMPDIR"/part_*.list; do
    idx=$(basename "$f" | sed -E 's/part_([0-9]+)\.list/\1/')
    cfg="$TMPDIR/cfg_$idx.json"
    out="$TMPDIR/output_part_$idx.root"
    python3 - "$CONFIG" "$f" "$out" "$cfg" <<'PYEOF'
import json, sys
src_cfg, src_list, out_root, dst_cfg = sys.argv[1:5]
with open(src_cfg) as fh:
    c = json.load(fh)
c['input']['source']    = src_list
c['output']['filename'] = out_root
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
N_FINISHED=0

FIRST_JOB=1
for cfg in "$TMPDIR"/cfg_*.json; do
    # Throttle: wait until a slot frees up
    while (( $(jobs -rp | wc -l) >= N_CONCURRENT )); do
        wait -n 2>/dev/null || true
    done
    # Stagger launches (avoid all jobs hitting chain->GetEntries() simultaneously)
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
                # Backoff between retries — most ROOT/HDD glitches resolve after a pause.
                sleep $((5 + RANDOM % 10))
                log="${cfg%.json}_retry${try}.log"
                printf '  [retry %d] %s\n' "$try" "$(basename "$cfg")"
            fi
            if ./ana "$cfg" > "$log" 2>&1; then
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

# Wait for the final wave
wait

ELAPSED=$((SECONDS - START))
echo
echo "All jobs finished in ${ELAPSED}s ($(printf '%.1f' "$(echo "$ELAPSED/60" | bc -l)") min)."

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

# Cleanup only when ALL jobs succeeded (no missing outputs)
if (( ${#OUTPUTS[@]} == CFG_COUNT )); then
    rm -rf "$TMPDIR"
else
    echo
    echo "NOTE: ${#OUTPUTS[@]}/$CFG_COUNT chunks succeeded — tmp kept at $TMPDIR"
    echo "      Failed chunks: inspect cfg_NNN.log and rerun those configs:"
    for cfg in "$TMPDIR"/cfg_*.json; do
        idx=$(basename "$cfg" .json | sed 's/cfg_//')
        out="$TMPDIR/output_part_${idx}.root"
        [[ -f "$out" ]] || echo "        ./ana $cfg"
    done
fi

echo "================================================================"
echo "Done. Output: $OUTPUT"
echo "================================================================"
