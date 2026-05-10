#!/bin/bash
# ============================================================================
# Live monitor for run_parallel_scan.sh — auto-finds the latest
# fat_scan_parallel.* dir and displays per-chunk status.
#
# Usage:
#   ./monitor_parallel_scan.sh             # auto-detect, refresh 5s
#   ./monitor_parallel_scan.sh /tmp/...    # explicit tmp dir
#   ./monitor_parallel_scan.sh /tmp/... 2  # custom refresh interval (s)
# ============================================================================

set -uo pipefail

TMP="${1:-}"
INTERVAL="${2:-5}"

if [[ -z "$TMP" ]]; then
    candidates=()
    [[ -n "${FAT_TMPDIR:-}" ]] && candidates+=("$FAT_TMPDIR/fat_scan_parallel."*)
    candidates+=(/mnt/scratch/fat_tmp/fat_scan_parallel.*  \
                 /tmp/fat_scan_parallel.*                  \
                 "$HOME/tmp"/fat_scan_parallel.*)

    TMP=$(ls -dt "${candidates[@]}" 2>/dev/null | head -1)
    [[ -d "$TMP" ]] || {
        echo "No fat_scan_parallel.* dir found in any known location."
        echo "Searched: \$FAT_TMPDIR, /mnt/scratch/fat_tmp, /tmp, \$HOME/tmp"
        echo "Pass the path explicitly: $0 /path/to/fat_scan_parallel.XXXXXX"
        exit 1
    }
fi

[[ -d "$TMP" ]] || { echo "Not a directory: $TMP"; exit 1; }

# A chunk is "done" when its log contains the "Scan Complete!" marker
# (printed by trigger_scan right before exit). NOT just when output_part_NNN.root
# exists — RECREATE creates the file early.
is_done() {
    local cfg="$1"
    local log
    log=$(ls -t "${cfg%.json}"*.log 2>/dev/null | head -1)
    [[ -n "$log" ]] && grep -q "Scan Complete" "$log" 2>/dev/null
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
            # trigger_scan prints "[NN/TOTAL] N=… PT3=… PT2=… <path>" lines —
            # show the latest one as the progress indicator.
            line=$(tail -c 800 "$log" 2>/dev/null \
                   | grep -oE '\[ *[0-9]+/ *[0-9]+\][^[:cntrl:]]*' | tail -1)
            running_lines+=$(printf '  cfg_%s  %s\n' "$idx" "${line:-(starting…)}")$'\n'
        else
            WAITING=$((WAITING + 1))
        fi
    done

    clear
    printf '== monitor (scan): %s ==\n' "$TMP"
    printf '== %d done / %d running / %d waiting    (total %d, %s) ==\n\n' \
           "$DONE" "$RUNNING" "$WAITING" "$TOTAL" "$(date '+%H:%M:%S')"

    if (( RUNNING > 0 )); then
        printf 'Running:\n'
        printf '%s' "$running_lines"
        printf '\n'
    fi

    if (( DONE > 0 )); then
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
