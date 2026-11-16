#!/usr/bin/env python3
# ============================================================================
# v5_figE2_prep.py -- input table for v5_figE2.C (corrected drawing of WP5 Fig. E2: legend no longer hides the S5 data)
# Numbers are taken from WP5's own output (results.txt 'bin' lines: N, k CB, fitted background, fitted eta parts) and
# from the v5 dump of hists_v1.root (event-group variance of N, CB replica variance, raw templates at SMASH x PDG BR).
# Writes /hdd3/przygoda/HADES/FAT_out/all_report/v5/v5_E2_input.txt. Run: python3 plots/all_report/v5_figE2_prep.py
# ============================================================================
import numpy as np
D = "/hdd3/przygoda/HADES/FAT_out/all_report/v5/"
NB2RAW = 6079.2 * 63.0 / 64.0
rows = {}
kE2 = None
for line in open("/hdd3/przygoda/HADES/FAT_out/all_report/wp5/results.txt"):
    p = line.split()
    if p and p[0] == "values" and p[1] == "FMAIN": kE2 = float(p[13])
    if p and p[0] == "bin" and p[1].endswith("E2"):
        rows[(p[1], int(p[2]))] = [float(x) for x in p[3:]]
def rd(fn):
    d = {}
    for line in open(fn):
        q = line.split(); d[q[0]] = np.array([float(x) for x in q[1:]])
    return d
def rep(fn):
    d = {}
    for line in open(fn):
        q = line.split(); d.setdefault(q[0], []).append([float(x) for x in q[2:]])
    return {k: np.array(v) for k, v in d.items()}
dat = rd(D + "data_S0.txt"); cbr = rep(D + "cbrep_S0.txt"); sim = rd(D + "sim_S0.txt")
out = open(D + "v5_E2_input.txt", "w")
out.write("# slice xbin lo hi D varD C2 C34 C5 Geta smash   (D = N - kCB - background; raw rows per bin)\n")
for j in range(1, 6):
    h = "F10_S%dE2" % j
    vcb = cbr["CBrep_" + h].var(axis=0)          # population variance over the 200 replicas (as WP5)
    eg = dat["Nerr_" + h] ** 2
    for x in range(1, 71):
        lo, hi, N, CB, kCB, bg, bgs, bgm, c2, c3, c4, c5, ge, mu, infit = rows[(h, x)]
        Dv = N - kCB - bg
        varD = eg[x - 1] + kE2 ** 2 * vcb[x - 1]
        sm = sum(sim[c + "_" + h][x - 1] for c in ("C2", "C3", "C4", "C5", "Geta")) * NB2RAW
        out.write("%d %d %.5f %.5f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n" % (j - 1, x, lo, hi, Dv, varD, c2, c3 + c4, c5, ge, sm))
out.close()
print("k_E2 =", kE2, "rows:", len(rows))
