#!/bin/bash
# ============================================================================
# Live monitor for run_parallel.sh — auto-finds the latest fat_parallel.* dir
# in /tmp and displays per-chunk status (done / latest ETA line).
#
# Usage:
#   ./monitor_parallel.sh             # auto-detect latest tmp dir, refresh 5s
#   ./monitor_parallel.sh /tmp/...    # explicit tmp dir
#   ./monitor_parallel.sh /tmp/... 2  # custom refresh interval (s)
# ============================================================================

set -uo pipefail

TMP="${1:-}"
INTERVAL="${2:-5}"

if [[ -z "$TMP" ]]; then
    # Search known launcher work-dir locations and pick the most recent.
    # Honors FAT_TMPDIR env var if set (matches run_parallel.sh convention).
    candidates=()
    [[ -n "${FAT_TMPDIR:-}" ]] && candidates+=("$FAT_TMPDIR/fat_parallel."*)
    candidates+=(/mnt/scratch/fat_tmp/fat_parallel.*  \
                 /tmp/fat_parallel.*                  \
                 "$HOME/tmp"/fat_parallel.*)

    TMP=$(ls -dt "${candidates[@]}" 2>/dev/null | head -1)
    [[ -d "$TMP" ]] || {
        echo "No fat_parallel.* dir found in any known location."
        echo "Searched: \$FAT_TMPDIR, /mnt/scratch/fat_tmp, /tmp, \$HOME/tmp"
        echo "Pass the path explicitly: $0 /path/to/fat_parallel.XXXXXX"
        exit 1
    }
fi

[[ -d "$TMP" ]] || { echo "Not a directory: $TMP"; exit 1; }

# A chunk is "done" when its log contains the final "Analysis Complete!" marker
# (printed by ana right before exit). NOT just when output_part_NNN.root exists —
# that file is created up-front by RECREATE and stays tiny until the end.
is_done() {
    local cfg="$1"
    # Check the most recent log for this chunk (any of base.log or _retryN.log)
    local log
    log=$(ls -t "${cfg%.json}"*.log 2>/dev/null | head -1)
    [[ -n "$log" ]] && grep -q "Analysis Complete" "$log" 2>/dev/null
}

while true; do
    TOTAL=$(ls "$TMP"/cfg_*.json 2>/dev/null | wc -l)

    DONE=0; RUNNING=0; WAITING=0
    running_lines=""; recent_done=""

    for cfg in "$TMP"/cfg_*.json; do
        idx=$(basename "$cfg" .json | sed 's/cfg_//')
        log=$(ls -t "${cfg%.json}"*.log 2>/dev/null | head -1)

        if is_done "$cfg"; then
            DONE=$((DONE + 1))
            out="$TMP/output_part_${idx}.root"
            sz="?"
            [[ -f "$out" ]] && sz=$(stat -c %s "$out" | numfmt --to=iec)
            recent_done+=$(printf '  cfg_%s  [done] %s\n' "$idx" "$sz")$'\n'
        elif [[ -n "$log" ]]; then
            RUNNING=$((RUNNING + 1))
            line=$(tail -c 400 "$log" 2>/dev/null | tr '\r' '\n' \
                   | grep -oE '[0-9.]+ *% \| ETA: [0-9:-]+' | tail -1)
            running_lines+=$(printf '  cfg_%s  %s\n' "$idx" "${line:-(starting…)}")$'\n'
        else
            WAITING=$((WAITING + 1))
        fi
    done

    clear
    printf '== monitor: %s ==\n' "$TMP"
    printf '== %d done / %d running / %d waiting    (total %d, %s) ==\n\n' \
           "$DONE" "$RUNNING" "$WAITING" "$TOTAL" "$(date '+%H:%M:%S')"

    if (( RUNNING > 0 )); then
        printf 'Running:\n'
        printf '%s' "$running_lines"
        printf '\n'
    fi

    if (( DONE > 0 )); then
        # Show only the last 10 completed chunks (avoids scroll spam at high N)
        N_SHOW=10
        printf 'Last %d done (of %d):\n' "$N_SHOW" "$DONE"
        printf '%s' "$recent_done" | tail -n "$N_SHOW"
    fi

    if (( DONE >= TOTAL && TOTAL > 0 )); then
        echo
        echo "All chunks done. (Launcher should be hadd-ing now.)"
        break
    fi

    sleep "$INTERVAL"
done
