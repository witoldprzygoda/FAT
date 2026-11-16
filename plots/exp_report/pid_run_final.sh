#!/bin/bash
# pid_run_final.sh -- PID chapter: run pid_final_fill.C on all 29 day files of the betanew PAT production
# (6 processes in parallel) and merge the outputs with hadd.
# Inputs : /hdd3/przygoda/HADES/files/exp/gen4/lep_betanew.list (read-only), the compact scan file and the
#          data polygons (paths below).
# Outputs: $OUT/final/pid_final_<day>.root, $OUT/pid_final.root (merged), $OUT/final/log_<day>.txt
# Run    : bash /home/przygoda/HADES/FAT/FAT/plots/exp_report/pid_run_final.sh
set -e
MAC=/home/przygoda/HADES/FAT/FAT/plots/exp_report/pid_final_fill.C
OUT=/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/pid
BUILD=/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/aclic_pid
LIST=/hdd3/przygoda/HADES/files/exp/gen4/lep_betanew.list
COMPACT=/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/scan_v2/scan_compact.root
CUTS=/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/pid/inputs/pp45_pid_cuts_exp_betanew_v1.root
mkdir -p $OUT/final $BUILD
# compile once (a test on 10 rows), then run in parallel on the compiled library
root -l -b -q -e "gSystem->SetBuildDir(\"$BUILD\", kTRUE)" "$MAC+(\"$(head -1 $LIST)\", \"$COMPACT\", \"$CUTS\", \"$OUT/final/test.root\", 10)" > $OUT/final/log_compile.txt 2>&1
rm -f $OUT/final/test.root
grep -v '^#' $LIST | xargs -P 6 -I{} bash -c '
  f={}; d=$(basename $f .root); d=${d#lep_}
  root -l -b -q -e "gSystem->SetBuildDir(\"'$BUILD'\", kTRUE)" "'$MAC'+(\"$f\", \"'$COMPACT'\", \"'$CUTS'\", \"'$OUT'/final/pid_final_$d.root\")" > '$OUT'/final/log_$d.txt 2>&1'
hadd -f $OUT/pid_final.root $OUT/final/pid_final_0*.root > $OUT/final/log_hadd.txt 2>&1
grep -h "PipPimEpEm_ID:" $OUT/final/log_0*.txt
