#!/usr/bin/env python3
# ============================================================================
# wp2_lumi.py  --  WP2: luminosity per day and per period for the ANALYSED runs
# ============================================================================
# Reads the luminosity CSVs of the START-free chain "N" (copies in the scratch
# directory final/remote/, A3 note) and writes
#   docs/latex/tables_all/wp2_lumi_perday.txt   one line per day (CSV, '#' comments)
#   docs/latex/tables_all/wp2_lumi_perperiod.txt per period / early / late / all
# Analysed runs = runs with at least one merged DST in the FAT input
# (lumi_coverage_perrun.csv, n_dst_merged > 0): 17 557 runs, 6079.2 nb^-1.
# No partial-run scaling (spec D7: L = 6079.2).
# Usage: python3 -I plots/all_report/wp2_lumi.py
# ============================================================================
import csv, math, os, sys

REMOTE = "/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/remote"
TABDIR = "/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all"
RHO = 1.0973   # rho(sector 2) of the luminosity chain (per-run CSV header, FACTORS line)

def rows(path):
    with open(path) as f:
        lines = [l for l in f if not l.startswith("#")]
    return list(csv.DictReader(lines))

cov = rows(os.path.join(REMOTE, "lumi_coverage_perrun.csv"))
run = rows(os.path.join(REMOTE, "lumi_startfree_perrun.csv"))
day = rows(os.path.join(REMOTE, "lumi_startfree_perday.csv"))
per = rows(os.path.join(REMOTE, "sf_luminosity_perperiod.csv"))

# stat of each run (L6_stat) from the per-run file, keyed by dst_stem
stat_of = {r["dst_stem"]: float(r["L6_stat"]) for r in run}

D = {}
for r in cov:
    d = r["day"]
    e = D.setdefault(d, dict(nrun=0, nan=0, Lall=0.0, Lan=0.0, s2an=0.0, Lpart=0.0))
    L = float(r["L6_SF"])
    e["nrun"] += 1
    e["Lall"] += L
    if int(r["n_dst_merged"]) > 0:
        e["nan"] += 1
        e["Lan"] += L
        e["s2an"] += stat_of[r["dst_stem"]] ** 2
        e["Lpart"] += L * int(r["n_dst_merged"]) / int(r["n_dst_list"])

def period_of(d):
    d = int(d)
    if 34 <= d <= 46: return "I"
    if 47 <= d <= 48: return "IIA"
    if 49 <= d <= 53: return "IIB"
    if 59 <= d <= 67: return "III"
    return "?"

dd = {r["day"]: r for r in day}
os.makedirs(TABDIR, exist_ok=True)
out = open(os.path.join(TABDIR, "wp2_lumi_perday.txt"), "w")
out.write("# wp2_lumi.py: luminosity per day, START-free chain N (lumi_startfree_perday.csv, lumi_coverage_perrun.csv)\n")
out.write("# L6_all = all runs of the day; L6_an = analysed runs (>= 1 merged DST in the FAT input); L6_an_part = with the merged-DST fraction\n")
out.write("# S2fac6 = sector-2 PT2 factor already INSIDE L (days 034-044): sector 2 replaced by rho x mean of the other five sectors, rho = 1.0973;\n")
out.write("#   eps2_PT2 = sector-2 elastic efficiency implied by S2fac6: ((5 + rho)/S2fac6 - 5)/rho (1 = no deficit)\n")
out.write("# w6 = target-window fraction removed from L (per period); LF2, LF3 = PT2/PT3 DAQ live fractions; LF3/LF2 x 64/63\n")
out.write("day,period,nrun,nrun_an,L6_all,L6_an,L6_an_stat,L6_an_part,S2fac6,eps2_PT2,w6,LF2,LF3,LF3oLF2x64o63,early\n")
for d in sorted(D):
    e = D[d]; r = dd[d]
    s2 = float(r["S2fac6"])
    out.write("%s,%s,%d,%d,%.4f,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%d\n" % (
        d, period_of(d), e["nrun"], e["nan"], e["Lall"], e["Lan"], math.sqrt(e["s2an"]), e["Lpart"], s2, ((5.0 + RHO) / s2 - 5.0) / RHO,
        float(r["w6"]), float(r["LF2"]), float(r["LF3"]), float(r["LF3_over_LF2"]) * 64.0 / 63.0, 1 if int(d) <= 44 else 0))
out.close()

# per period, early (034-044), late (045-067), all
groups = [("I", lambda d: period_of(d) == "I"), ("IIA", lambda d: period_of(d) == "IIA"),
          ("IIB", lambda d: period_of(d) == "IIB"), ("III", lambda d: period_of(d) == "III"),
          ("early_034_044", lambda d: int(d) <= 44), ("late_045_067", lambda d: int(d) >= 45),
          ("I_late_045_046", lambda d: 45 <= int(d) <= 46), ("all", lambda d: True)]
pp = {r["period"]: r for r in per}
out = open(os.path.join(TABDIR, "wp2_lumi_perperiod.txt"), "w")
out.write("# wp2_lumi.py: luminosity per group of days [nb^-1], START-free chain N\n")
out.write("# L6_an = analysed runs (sum = 6079.2 = alldefs::kL); stat = Poisson of the elastic count (per-run L6_stat in quadrature)\n")
out.write("# unc_sec, unc_hades = per-period 'Sectors (+ sector-2 correction)' and total HADES syst of sf_luminosity_perperiod.csv\n")
out.write("#   (all runs; NOT independent between periods; given only for the four periods and 'all')\n")
out.write("# LwS2 = L-weighted mean of S2fac6; Lw_w6 = L-weighted w6; LwLF = L-weighted LF3/LF2 x 64/63\n")
out.write("group,ndays,nrun_an,L6_all,L6_an,stat,unc_sec,unc_hades,LwS2,Lw_w6,LwLF\n")
for name, sel in groups:
    ds = [d for d in sorted(D) if sel(d)]
    La = sum(D[d]["Lan"] for d in ds); Lall = sum(D[d]["Lall"] for d in ds)
    st = math.sqrt(sum(D[d]["s2an"] for d in ds))
    s2 = sum(D[d]["Lan"] * float(dd[d]["S2fac6"]) for d in ds) / La
    w6 = sum(D[d]["Lan"] * float(dd[d]["w6"]) for d in ds) / La
    lf = sum(D[d]["Lan"] * float(dd[d]["LF3_over_LF2"]) * 64.0 / 63.0 for d in ds) / La
    key = {"I": "I", "IIA": "IIA", "IIB": "IIB", "III": "III", "all": "all"}.get(name)
    us = float(pp[key]["L6_unc"]) if key else float("nan")
    uh = float(pp[key]["L6_hades"]) if key else float("nan")
    out.write("%s,%d,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.5f,%.5f,%.5f\n" % (
        name, len(ds), sum(D[d]["nan"] for d in ds), Lall, La, st, us, uh, s2, w6, lf))
out.close()
print(open(os.path.join(TABDIR, "wp2_lumi_perperiod.txt")).read())
