#!/bin/bash
# isbest_emulation_run.sh — run plots/isbest_emulation_betanew.C for every PAT
# file of a list and the three patterns (pm, pp, mm) in parallel, then hadd
# per pattern.  Usage:
#   plots/isbest_emulation_run.sh [list] [jobs] [R_pt3] [R_pt2] [pool_depth] [patterns] [tag] [vertex_z_min]
# Defaults: list = /hdd3/przygoda/HADES/files/exp/gen4/lep_betanew.list,
#           jobs = 8, R_pt3 = 2, R_pt2 = 10, pool_depth = 20, patterns = "pm pp mm",
#           tag = vtx, vertex_z_min = -500 (FAT event cut; -1e9 = none).
# Outputs:  plots/output/isbest_emu_<tag>/<pat>_<day>.root (per file) and
#           plots/output/isbest_emu_<tag>_<pat>.root (hadd over the list).
set -euo pipefail
cd "$(dirname "$0")/.."
LIST=${1:-/hdd3/przygoda/HADES/files/exp/gen4/lep_betanew.list}
JOBS=${2:-8}
R3=${3:-2}
R2=${4:-10}
POOL=${5:-20}
PATS=${6:-"pm pp mm"}
TAG=${7:-vtx}
VTX=${8:--500}
OUT=plots/output/isbest_emu_${TAG}
mkdir -p "$OUT"
# compile once (avoids concurrent ACLiC builds); the dummy call returns after the pattern check
root -l -b -q 'plots/isbest_emulation_betanew.C+("/dev/null","xx","/dev/null")' > "$OUT/compile.log" 2>&1 || true
test -f plots/isbest_emulation_betanew_C.so || { echo "compile failed, see $OUT/compile.log"; exit 1; }
# stale per-day outputs of the requested patterns are removed so that hadd sums only this run
for pat in $PATS; do find "$OUT" -maxdepth 1 -type f \( -name "${pat}_*.root" -o -name "${pat}_*.log" \) -delete; done
run_one() {
    local f=$1 pat=$2 r3=$3 r2=$4 pool=$5 out=$6 vtx=$7
    local day; day=$(basename "$f" .root); day=${day#lep_}
    local seed=$((10#$day))
    root -l -b -q "plots/isbest_emulation_betanew.C+(\"$f\",\"$pat\",\"$out/${pat}_${day}.root\",$r3,$r2,$pool,-1,$seed,\"/hdd2/przygoda/hades/pp45/exp/pp45_pid_cuts_exp_betanew_v1.root\",3000000,$vtx)" \
        > "$out/${pat}_${day}.log" 2>&1 \
        && echo "done $pat $day" || { mv -f "$out/${pat}_${day}.root" "$out/${pat}_${day}.root.failed" 2>/dev/null; echo "FAILED $pat $day (see $out/${pat}_${day}.log)"; }
}
export -f run_one
grep -v '^\s*$' "$LIST" | while read -r f; do for pat in $PATS; do echo "$f $pat"; done; done \
    | xargs -P "$JOBS" -n 2 bash -c 'run_one "$0" "$1" '"$R3 $R2 $POOL $OUT $VTX"
for pat in $PATS; do
    inputs=$(grep -v '^\s*$' "$LIST" | while read -r f; do day=$(basename "$f" .root); day=${day#lep_}; echo "$OUT/${pat}_${day}.root"; done | xargs ls 2>/dev/null)
    hadd -f "plots/output/isbest_emu_${TAG}_${pat}.root" $inputs > "$OUT/hadd_${pat}.log" 2>&1 && echo "hadd $pat ok ($(echo $inputs | wc -w) files)" || echo "hadd $pat FAILED"
done
