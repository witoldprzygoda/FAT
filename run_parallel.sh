#!/bin/bash
# ============================================================================
# FAT parallel launcher — split input list, run N ana jobs, hadd outputs.
#
# Reads config.json (or any other .json), splits the .list source into N
# chunks, generates a temporary config per chunk with its own output filename,
# launches up to N_CONCURRENT ana jobs in parallel, then merges all part
# outputs into the original output filename via hadd.
#
# Does NOT modify the project sources, the input config, or the input list.
# All temp files live in $TMPDIR (printed below); cleaned up on success only.
#
# Entry-range split: when the .list has exactly ONE line naming a plain file
# and N_parts > 1, the entries of that file are split instead of the list.
# N = entries of the configured tree (TChain), range [S, E) with
# S = input.start_event and E = min(S + input.max_events, N) (E = N when
# max_events <= 0), L = E - S. Part k gets start_event = S + floor(k L / P),
# max_events = floor((k+1) L / P) - floor(k L / P). Before the merge the
# ranges that the jobs report in their logs must cover [S, E) without gaps or
# overlaps. No entry split (the file runs as ONE part, with a warning) when
# trigger.bias_correction is on (the PT3 weight is computed per file over the
# entries a job reads) or when the line is not a plain file (pattern, URL).
#
# Merge: before hadd, an empty ntuple of a part whose columns differ from the
# parts with entries is rewritten with their columns (hadd takes the layout of
# the first file); after hadd every ntuple must have those columns and the sum
# of the part entries.
#
# Usage:
#   ./run_parallel.sh [config.json] [N_parts=32] [N_concurrent=min(N_parts,15)]
#                     [N_retries=2] [launch_stagger_s=3]
#   N_parts 1..999 (three-digit part numbers keep the hadd order).
#
# Examples:
#   ./run_parallel.sh                       # 32 parts, 15 concurrent
#   ./run_parallel.sh config.json 16 8      # 16 parts, max 8 concurrent
#   ./run_parallel.sh config.json 64 16     # 64 chunks, 16 at a time
#   ./run_parallel.sh config_sim_lep.json 15 15   # one file, 15 entry ranges
#
# Exit codes: 0 merged output written; 1 bad arguments, config, list or input
# (nothing launched); 2 part outputs missing; 3 hadd, its ntuple layout check
# or rename failed; 4 one or more parts failed; 5 entry ranges of the parts do
# not cover the input range (split mode); 129 / 130 / 143 launcher stopped by
# SIGHUP / SIGINT / SIGTERM (running parts are killed). No merge after any
# failure.
# ============================================================================

set -uo pipefail

CONFIG="${1:-config.json}"
N_PARTS="${2:-32}"
N_CONCURRENT="${3:-}"        # default min(N_PARTS, 15)
N_RETRIES="${4:-2}"          # per-chunk retries on failure
LAUNCH_STAGGER="${5:-3}"     # seconds between consecutive job starts (0 = no stagger)

# ---------- arguments (checked before anything is created) ----------
is_uint() { [[ "$1" =~ ^[0-9]+$ ]]; }
is_uint "$N_PARTS" && (( 10#$N_PARTS >= 1 && 10#$N_PARTS <= 999 )) || {
    echo "ERROR: N_parts must be an integer from 1 to 999 (got '$N_PARTS'). Nothing launched."; exit 1; }
N_PARTS=$(( 10#$N_PARTS ))
if [[ -z "$N_CONCURRENT" ]]; then
    N_CONCURRENT=$(( N_PARTS < 15 ? N_PARTS : 15 ))
fi
is_uint "$N_CONCURRENT" && (( 10#$N_CONCURRENT >= 1 )) || {
    echo "ERROR: N_concurrent must be a positive integer (got '$N_CONCURRENT'). Nothing launched."; exit 1; }
N_CONCURRENT=$(( 10#$N_CONCURRENT ))
is_uint "$N_RETRIES" || {
    echo "ERROR: N_retries must be an integer >= 0 (got '$N_RETRIES'). Nothing launched."; exit 1; }
N_RETRIES=$(( 10#$N_RETRIES ))
is_uint "$LAUNCH_STAGGER" || {
    echo "ERROR: launch stagger must be an integer number of seconds >= 0 (got '$LAUNCH_STAGGER'). Nothing launched."; exit 1; }
LAUNCH_STAGGER=$(( 10#$LAUNCH_STAGGER ))

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

# The merged output is written only at the very end: check its directory now.
OUTDIR=$(dirname "$OUTPUT")
[[ -d "$OUTDIR" && -w "$OUTDIR" ]] || {
    echo "ERROR: output directory missing or not writable: $OUTDIR"
    exit 1
}

# Called on every failure: a file left at the output path by an earlier run
# must not be taken for the result of this run.
warn_stale_output() {
    if [[ -e "$OUTPUT" ]]; then
        echo "WARNING: $OUTPUT exists from an EARLIER run" \
             "($(date -r "$OUTPUT" '+%F %T')) — it is NOT the result of this run."
    fi
}

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
  Parts:          $N_PARTS requested
  Concurrent:     $N_CONCURRENT
  Work dir:       $TMPDIR
================================================================
EOF

# Strip blank lines and comments before splitting (so each part has only real entries)
grep -vE '^\s*(#|//|$)' "$SOURCE" > "$TMPDIR/all.list"
TOTAL=$(wc -l < "$TMPDIR/all.list")
echo "Source has $TOTAL non-comment lines."

SPLIT_ENTRIES=0
if (( TOTAL == 1 && N_PARTS > 1 )); then
    SPLIT_ENTRIES=1
    # File path as NTupleReader::openFromList reads it: text between the first
    # two quotes (chain->Add("...")), otherwise the line without trailing ' ;\t'.
    IN_FILE=$(python3 -c '
import sys
line = open(sys.argv[1]).readline().strip()
q = line.find("\"")
if q >= 0:
    e = line.find("\"", q + 1)
    line = line[q + 1:e] if e > q else ""
print(line.rstrip(" \t;"))' "$TMPDIR/all.list")
    BIAS_CORR=$(python3 -c "
import json, sys
v = json.load(open(sys.argv[1])).get('trigger', {}).get('bias_correction', True)
print(1 if (v is True or v == 1 or str(v).lower() == 'true') else 0)" "$CONFIG")
    if (( BIAS_CORR )); then
        echo "WARNING: trigger.bias_correction is on: ana computes the PT3 weight of a file over the entries"
        echo "         one job reads, so an entry split would change trigger_corr. No entry split:"
        echo "         the single file runs as ONE part."
        SPLIT_ENTRIES=0
    elif [[ ! -f "$IN_FILE" ]]; then
        echo "WARNING: the line of the list is not a plain file ($IN_FILE): a pattern, URL or"
        echo "         file/tree name cannot be split by entries. No entry split: the list runs as ONE part."
        SPLIT_ENTRIES=0
    fi
fi

if (( SPLIT_ENTRIES )); then
    # ---------- entry-range split of a single input file ----------
    TREE=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['input'].get('tree_name', 'PPip_ID'))" "$CONFIG")
    read -r CFG_START CFG_MAX < <(python3 -c "
import json, sys
i = json.load(open(sys.argv[1]))['input']
print(int(float(i.get('start_event', 0))), int(float(i.get('max_events', -1))))" "$CONFIG")
    if (( CFG_START < 0 )); then
        echo "ERROR: invalid input.start_event $CFG_START (must be >= 0). Nothing launched."
        rm -rf "$TMPDIR"
        warn_stale_output
        exit 1
    fi

    # Entries of the configured tree, counted as ana counts them (TChain, highest cycle)
    COUNT_LOG="$TMPDIR/entries.log"
    # (root -e returns the value of the last expression, not a status: the
    # FAT_ENTRIES line and the absence of ROOT errors are checked instead.)
    root -l -b -q -e "TChain c(\"$TREE\"); c.Add(\"$IN_FILE\"); printf(\"FAT_ENTRIES %lld\\n\", c.GetEntries());" \
         > "$COUNT_LOG" 2>&1
    N_ENTRIES=$(grep -oE '^FAT_ENTRIES -?[0-9]+' "$COUNT_LOG" | awk '{print $2}')
    if [[ -z "$N_ENTRIES" ]] || (( N_ENTRIES <= 0 )) \
       || grep -qE '(Error|SysError|Fatal) in <' "$COUNT_LOG"; then
        echo "ERROR: cannot count the entries of tree '$TREE' in $IN_FILE (wrong tree name or file?):"
        sed 's/^/    /' "$COUNT_LOG"
        echo "Nothing launched. Tmp kept at $TMPDIR"
        warn_stale_output
        exit 1
    fi

    RANGE_START=$CFG_START
    RANGE_END=$N_ENTRIES
    if (( CFG_MAX > 0 && CFG_START + CFG_MAX < N_ENTRIES )); then
        RANGE_END=$(( CFG_START + CFG_MAX ))
    fi
    RANGE_LEN=$(( RANGE_END - RANGE_START ))
    if (( RANGE_LEN <= 0 )); then
        echo "ERROR: empty entry range [$RANGE_START, $RANGE_END) of $N_ENTRIES entries" \
             "(input.start_event $CFG_START, input.max_events $CFG_MAX). Nothing launched."
        rm -rf "$TMPDIR"
        warn_stale_output
        exit 1
    fi
    if (( N_PARTS > RANGE_LEN )); then
        echo "WARNING: requested $N_PARTS parts but only $RANGE_LEN entries. Capping to $RANGE_LEN."
        N_PARTS=$RANGE_LEN
    fi
    echo "Parts used:     $N_PARTS"
    echo "Entry-range split: tree $TREE, $N_ENTRIES entries in $IN_FILE;" \
         "range [$RANGE_START, $RANGE_END) ($RANGE_LEN entries) in $N_PARTS parts:"

    : > "$TMPDIR/ranges.txt"
    for (( k = 0; k < N_PARTS; k++ )); do
        idx=$(printf '%03d' "$k")
        a=$(( RANGE_START + k * RANGE_LEN / N_PARTS ))
        b=$(( RANGE_START + (k + 1) * RANGE_LEN / N_PARTS ))
        cfg="$TMPDIR/cfg_$idx.json"
        out="$TMPDIR/output_part_$idx.root"
        echo "$idx $a $b" >> "$TMPDIR/ranges.txt"
        printf '  part %s: entries [%d, %d)  (%d entries)\n' "$idx" "$a" "$b" $(( b - a ))
        python3 - "$CONFIG" "$TMPDIR/all.list" "$out" "$cfg" "$a" $(( b - a )) <<'PYEOF'
import json, sys
src_cfg, src_list, out_root, dst_cfg = sys.argv[1:5]
start, count = int(sys.argv[5]), int(sys.argv[6])
with open(src_cfg) as fh:
    c = json.load(fh)
c['input']['source']      = src_list
c['input']['start_event'] = start
c['input']['max_events']  = count
c['output']['filename']   = out_root
with open(dst_cfg, 'w') as fh:
    json.dump(c, fh, indent=2)
PYEOF
    done
else

if (( N_PARTS > TOTAL )); then
    echo "WARNING: requested $N_PARTS parts but only $TOTAL files. Capping to $TOTAL."
    N_PARTS=$TOTAL
fi
echo "Parts used:     $N_PARTS (empty chunks are dropped below)"

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

fi   # SPLIT_ENTRIES

CFG_COUNT=$(ls "$TMPDIR"/cfg_*.json | wc -l)
echo "Generated $CFG_COUNT chunked configs. Launching..."
echo "(Per-job stdout/stderr in $TMPDIR/cfg_NNN.log)"
echo

# ---------- stop the parts when the launcher is stopped ----------
# Without this, killing the launcher (Ctrl+C, kill, closed ssh session) would
# leave the ana jobs running as orphans. The part subshells are killed first
# (so that none starts a retry), then their ana processes. Nothing is merged.
on_signal() {
    local name=$1 code=$2 pids kids="" p
    trap - INT TERM HUP
    echo
    echo "INTERRUPTED by SIG$name: stopping the running parts ..."
    pids=$(jobs -p)
    for p in $pids; do kids+=" $(pgrep -P "$p" 2>/dev/null | tr '\n' ' ')"; done
    [[ -n "$pids" ]] && kill -TERM $pids 2>/dev/null
    [[ -n "${kids// /}" ]] && kill -TERM $kids 2>/dev/null
    wait 2>/dev/null
    echo "No merge: $OUTPUT not written. Tmp kept at $TMPDIR"
    warn_stale_output
    exit "$code"
}
trap 'on_signal INT 130' INT
trap 'on_signal TERM 143' TERM
trap 'on_signal HUP 129' HUP

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
            ./ana "$cfg" > "$log" 2>&1
            rc=$?
            if (( rc == 0 )); then
                success=1
                break
            fi
            # A retry cannot help: 2 = input variables missing in the tree,
            # 3 = exception in the event loop, 130 = interrupted (Ctrl+C).
            if (( rc == 2 || rc == 3 || rc == 130 )); then
                break
            fi
        done
        if (( success )); then
            # Positive record of success, required for the merge below; it
            # holds the log of the successful try (read by the range check).
            printf '%s\n' "$log" > "${cfg%.json}.ok"
            if (( try == 0 )); then
                printf '  [OK]   %s\n' "$(basename "$cfg")"
            else
                printf '  [OK after %d retries] %s\n' "$try" "$(basename "$cfg")"
            fi
        else
            printf '  [FAIL after %d tries, exit code %d] %s  (see %s and *_retry*.log)\n' \
                $((try + 1)) "$rc" "$(basename "$cfg")" "$base_log"
            # Informational marker; the merge below requires cfg_NNN.ok.
            touch "${cfg%.json}.failed"
        fi
    ) &
done

# Wait for the final wave
wait

ELAPSED=$((SECONDS - START))
echo
echo "All jobs finished in ${ELAPSED}s ($((ELAPSED / 60)) min $((ELAPSED % 60)) s)."

# ---------- verify and merge ----------
# A part is good only if one of its tries returned zero (marker cfg_NNN.ok).
# A part without the marker has failed: all its tries returned non-zero, or
# its job was killed. Any failed part: name it, do NOT merge, keep the work
# dir with the logs.
shopt -s nullglob
FAILED=()
for cfg in "$TMPDIR"/cfg_*.json; do
    [[ -e "${cfg%.json}.ok" ]] || FAILED+=("$cfg")
done
OUTPUTS=("$TMPDIR"/output_part_*.root)
shopt -u nullglob

if (( ${#FAILED[@]} > 0 )); then
    echo "ERROR: ${#FAILED[@]} of $CFG_COUNT parts failed — no hadd, $OUTPUT not written."
    for cfg in "${FAILED[@]}"; do
        idx=$(basename "$cfg" .json | sed 's/cfg_//')
        echo "    part $idx: logs $TMPDIR/cfg_${idx}.log (and cfg_${idx}_retry*.log);" \
             "rerun with: ./ana $TMPDIR/cfg_${idx}.json"
    done
    warn_stale_output
    echo "Tmp kept at $TMPDIR"
    exit 4
fi

# Split mode: the entry ranges the jobs report ("Processing events A to B")
# must be the planned ones and cover [RANGE_START, RANGE_END) without gaps or
# overlaps; every job must report that it processed all of its entries.
if (( SPLIT_ENTRIES )); then
    RANGE_ERR=0
    expect=$RANGE_START
    while read -r idx a b; do
        okfile="$TMPDIR/cfg_$idx.ok"
        log=$(head -1 "$okfile")
        line=$(grep -aoE 'Processing events [0-9]+ to [0-9]+' "$log" | tail -1)
        ra=$(awk '{print $3}' <<< "$line")
        rb=$(awk '{print $5}' <<< "$line")
        done_n=$(grep -aoE 'Events processed: [0-9]+' "$log" | tail -1 | awk '{print $3}')
        if [[ -z "$ra" || -z "$rb" || -z "$done_n" ]]; then
            echo "ERROR: part $idx: no 'Processing events' / 'Events processed' line in $log"
            RANGE_ERR=1
            continue
        fi
        if (( ra != a || rb != b )); then
            echo "ERROR: part $idx processed entries [$ra, $rb), planned [$a, $b)"
            RANGE_ERR=1
        fi
        if (( ra != expect )); then
            echo "ERROR: part $idx starts at entry $ra, expected $expect (gap or overlap)"
            RANGE_ERR=1
        fi
        if (( done_n != rb - ra )); then
            echo "ERROR: part $idx processed $done_n entries, its range has $(( rb - ra ))"
            RANGE_ERR=1
        fi
        expect=$rb
    done < <(sort -n -k2 "$TMPDIR/ranges.txt")
    if (( expect != RANGE_END )); then
        echo "ERROR: the parts end at entry $expect, the range ends at $RANGE_END"
        RANGE_ERR=1
    fi
    if (( RANGE_ERR )); then
        echo "ERROR: entry ranges of the parts do not cover [$RANGE_START, $RANGE_END) — no hadd, $OUTPUT not written."
        warn_stale_output
        echo "Tmp kept at $TMPDIR"
        exit 5
    fi
    echo "Entry ranges checked: the $CFG_COUNT parts cover [$RANGE_START, $RANGE_END) without gaps or overlaps."
fi

if (( ${#OUTPUTS[@]} == 0 )); then
    echo "ERROR: no output_part_*.root produced — aborting (tmp kept at $TMPDIR)"
    warn_stale_output
    exit 2
fi
if (( ${#OUTPUTS[@]} != CFG_COUNT )); then
    echo "ERROR: expected $CFG_COUNT outputs, found ${#OUTPUTS[@]} — no hadd (tmp kept at $TMPDIR)"
    warn_stale_output
    exit 2
fi

# ---------- ntuple layouts ----------
# hadd takes the layout of every tree from the first file that has it. A part
# with no row of an ntuple holds an empty TNtuple with only the columns set so
# far (trigger_corr), and as the first file it would make hadd drop every
# other column of that ntuple without an error. Before hadd such empty
# ntuples get the columns of the parts with rows; after hadd every ntuple must
# have those columns and the sum of the part entries.
printf '%s\n' "${OUTPUTS[@]}" > "$TMPDIR/hadd.list"
cat > "$TMPDIR/fat_merge_layout.C" <<'CEOF'
// fat_merge_layout.C -- written by run_parallel.sh (named macro, interpreted).
// fat_merge_layout(list, merged): list = part outputs in hadd order, one per line.
//   merged == "": prepare. An EMPTY TNtuple of a part whose columns differ from those of the parts with
//                 entries is rewritten (empty, with their columns). Parts with entries must agree.
//   merged != "": check. Every tree of the merged file has the columns of the parts with entries and the
//                 sum of the part entries.
// Prints FAT_LAYOUT_OK / FAT_MERGE_OK on success; FAT_LAYOUT_* / FAT_MERGE_ERROR lines otherwise.
#include <TFile.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TKey.h>
#include <TClass.h>
#include <TObjArray.h>
#include <TDirectory.h>
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <vector>

struct FatTreeInfo { Long64_t entries; std::vector<std::string> cols; std::string cls, title; };

void fatCollect(TDirectory* d, const std::string& path, std::map<std::string, FatTreeInfo>& m) {
  std::set<std::string> names;
  TIter it(d->GetListOfKeys());
  TKey* k;
  while ((k = (TKey*)it())) names.insert(k->GetName());
  for (const std::string& name : names) {
    TKey* key = d->GetKey(name.c_str());                 // highest cycle
    TClass* c = key ? TClass::GetClass(key->GetClassName()) : nullptr;
    if (!c) continue;
    if (c->InheritsFrom(TDirectory::Class())) {
      fatCollect((TDirectory*)key->ReadObj(), path + name + "/", m);
      continue;
    }
    if (!c->InheritsFrom(TTree::Class())) continue;
    TTree* t = (TTree*)key->ReadObj();
    FatTreeInfo info;
    info.entries = t->GetEntries();
    info.cls = key->GetClassName();
    info.title = t->GetTitle();
    TObjArray* leaves = t->GetListOfLeaves();
    for (int i = 0; i < leaves->GetEntries(); ++i) info.cols.push_back(leaves->At(i)->GetName());
    m[path + name] = info;
    delete t;
  }
}

std::string fatJoin(const std::vector<std::string>& v, const char* sep) {
  std::string s;
  for (size_t i = 0; i < v.size(); ++i) { if (i) s += sep; s += v[i]; }
  return s;
}

void fat_merge_layout(const char* listfn, const char* merged) {
  std::vector<std::string> parts;
  std::ifstream in(listfn);
  std::string line;
  while (std::getline(in, line)) if (!line.empty()) parts.push_back(line);
  std::vector<std::map<std::string, FatTreeInfo> > info(parts.size());
  std::set<std::string> trees;
  for (size_t p = 0; p < parts.size(); ++p) {
    TFile* f = TFile::Open(parts[p].c_str(), "READ");
    if (!f || f->IsZombie()) { printf("FAT_LAYOUT_ERROR cannot open %s\n", parts[p].c_str()); return; }
    fatCollect(f, "", info[p]);
    f->Close();
    delete f;
    for (std::map<std::string, FatTreeInfo>::const_iterator it = info[p].begin(); it != info[p].end(); ++it) trees.insert(it->first);
  }
  // reference columns per tree: those of the first part with entries
  std::map<std::string, std::vector<std::string> > ref;
  std::map<std::string, Long64_t> sum;
  bool ok = true;
  for (const std::string& tr : trees) {
    sum[tr] = 0;
    bool haveRef = false;
    for (size_t p = 0; p < parts.size(); ++p) {
      std::map<std::string, FatTreeInfo>::const_iterator it = info[p].find(tr);
      if (it == info[p].end()) continue;
      sum[tr] += it->second.entries;
      if (it->second.entries == 0) continue;
      if (!haveRef) { ref[tr] = it->second.cols; haveRef = true; continue; }
      if (it->second.cols != ref[tr]) {
        printf("FAT_LAYOUT_ERROR %s: part %s has %zu columns with entries, an earlier part %zu (columns differ)\n",
               tr.c_str(), parts[p].c_str(), it->second.cols.size(), ref[tr].size());
        ok = false;
      }
    }
    if (!haveRef) {                                     // no part has rows: take the first part's columns
      for (size_t p = 0; p < parts.size() && !haveRef; ++p) {
        std::map<std::string, FatTreeInfo>::const_iterator it = info[p].find(tr);
        if (it != info[p].end()) { ref[tr] = it->second.cols; haveRef = true; }
      }
    }
  }
  if (!ok) return;
  if (std::string(merged).empty()) {
    for (size_t p = 0; p < parts.size(); ++p) {
      for (std::map<std::string, FatTreeInfo>::const_iterator it = info[p].begin(); it != info[p].end(); ++it) {
        const FatTreeInfo& ti = it->second;
        if (ti.entries != 0 || ti.cols == ref[it->first]) continue;
        if (ti.cls != "TNtuple") {
          printf("FAT_LAYOUT_ERROR %s in %s: empty %s with other columns, cannot be rewritten\n", it->first.c_str(), parts[p].c_str(), ti.cls.c_str());
          return;
        }
        TFile f(parts[p].c_str(), "UPDATE");
        TDirectory* d = &f;
        std::string name = it->first;
        const size_t sl = name.rfind('/');
        if (sl != std::string::npos) { d = f.GetDirectory(name.substr(0, sl).c_str()); name = name.substr(sl + 1); }
        if (!d) { printf("FAT_LAYOUT_ERROR %s: directory missing in %s\n", it->first.c_str(), parts[p].c_str()); return; }
        d->cd();
        d->Delete((name + ";*").c_str());
        TNtuple* nt = new TNtuple(name.c_str(), ti.title.c_str(), fatJoin(ref[it->first], ":").c_str());
        nt->Write();
        f.Close();
        printf("FAT_LAYOUT_FIXED %s in %s: empty, %zu -> %zu columns\n", it->first.c_str(), parts[p].c_str(), ti.cols.size(), ref[it->first].size());
      }
    }
    printf("FAT_LAYOUT_OK %zu trees in %zu parts\n", trees.size(), parts.size());
    return;
  }
  std::map<std::string, FatTreeInfo> mi;
  TFile* fm = TFile::Open(merged, "READ");
  if (!fm || fm->IsZombie()) { printf("FAT_MERGE_ERROR cannot open %s\n", merged); return; }
  fatCollect(fm, "", mi);
  fm->Close();
  delete fm;
  for (const std::string& tr : trees) {
    std::map<std::string, FatTreeInfo>::const_iterator it = mi.find(tr);
    if (it == mi.end()) { printf("FAT_MERGE_ERROR %s missing in the merged file\n", tr.c_str()); ok = false; continue; }
    if (it->second.entries != sum[tr]) {
      printf("FAT_MERGE_ERROR %s: %lld entries merged, the parts have %lld\n", tr.c_str(), it->second.entries, sum[tr]);
      ok = false;
    }
    if (it->second.cols != ref[tr]) {
      printf("FAT_MERGE_ERROR %s: %zu columns merged, the parts have %zu\n", tr.c_str(), it->second.cols.size(), ref[tr].size());
      ok = false;
    }
    printf("FAT_MERGE_TREE %s entries %lld columns %zu\n", tr.c_str(), it->second.entries, it->second.cols.size());
  }
  if (ok) printf("FAT_MERGE_OK %zu trees\n", trees.size());
}
CEOF
root -l -b -q "$TMPDIR/fat_merge_layout.C(\"$TMPDIR/hadd.list\", \"\")" > "$TMPDIR/layout.log" 2>&1
if ! grep -q '^FAT_LAYOUT_OK' "$TMPDIR/layout.log"; then
    echo "ERROR: the ntuple layouts of the parts cannot be merged safely — no hadd, $OUTPUT not written:"
    sed 's/^/    /' "$TMPDIR/layout.log"
    warn_stale_output
    echo "Tmp kept at $TMPDIR"
    exit 3
fi
grep '^FAT_LAYOUT_FIXED' "$TMPDIR/layout.log" | sed 's/^FAT_LAYOUT_FIXED/Layout fixed before hadd:/'

# Merge into a temporary name next to the target and rename on success: a
# failed hadd must not leave a truncated file at the output path.
MERGE_TMP="${OUTPUT%.root}.merging.root"
echo "Merging ${#OUTPUTS[@]} files into $OUTPUT ..."
hadd -f "$MERGE_TMP" "${OUTPUTS[@]}"
HADD_RC=$?

if (( HADD_RC != 0 )); then
    rm -f "$MERGE_TMP"
    echo "ERROR: hadd failed (rc=$HADD_RC) — $OUTPUT not written, tmp kept at $TMPDIR"
    warn_stale_output
    exit 3
fi
root -l -b -q "$TMPDIR/fat_merge_layout.C(\"$TMPDIR/hadd.list\", \"$MERGE_TMP\")" > "$TMPDIR/merge_check.log" 2>&1
if ! grep -q '^FAT_MERGE_OK' "$TMPDIR/merge_check.log"; then
    rm -f "$MERGE_TMP"
    echo "ERROR: the merged ntuples do not have the columns or the entries of the parts — $OUTPUT not written:"
    sed 's/^/    /' "$TMPDIR/merge_check.log"
    warn_stale_output
    echo "Tmp kept at $TMPDIR"
    exit 3
fi
echo "Merged ntuples checked (columns and entries):"
grep '^FAT_MERGE_TREE' "$TMPDIR/merge_check.log" | sed 's/^FAT_MERGE_TREE/   /'
if ! mv -f "$MERGE_TMP" "$OUTPUT"; then
    echo "ERROR: could not rename $MERGE_TMP to $OUTPUT — tmp kept at $TMPDIR"
    warn_stale_output
    exit 3
fi

# All parts succeeded and were merged. The job logs go away with the work
# dir: show first what ROOT or ana warned about in them (errors of a try that
# was repeated are listed too).
LOG_NOTES=$(grep -ohE $'(Warning|Error|SysError) in <[^\r]*|Warning: [^\r]*' \
                 "$TMPDIR"/cfg_*.log 2>/dev/null | sort | uniq -c)
if [[ -n "$LOG_NOTES" ]]; then
    echo "Messages found in the job logs (count, text):"
    echo "$LOG_NOTES" | head -30
fi

# Remove the work dir
rm -rf "$TMPDIR"

echo "================================================================"
echo "Done. Output: $OUTPUT"
echo "================================================================"
