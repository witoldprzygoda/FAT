#!/usr/bin/env python3
# ============================================================================
# wp2_csec.py  --  WP2: correction C_sec for the sector-2 deficit of days 034-045
# ============================================================================
# Inputs (written by wp2_lumi.py, wp2_sector.C, wp2_simfs.C):
#   docs/latex/tables_all/wp2_lumi_perday.txt   L6 of the analysed runs per day
#   docs/latex/tables_all/wp2_sector_eps.txt    eps per sector pattern (lines "tm|S0|0|N|q|<cat>,eps,err,share")
#   docs/latex/tables_all/wp2_simfs.txt         simulated shares (lines "fs|<group>|<window>|<cat>,f,err,W")
# Correction (plain words): in the deficit days a row with a track in sector 2 survives with probability eps; rows
# without such a track are not affected (measured). If a share f of the true signal has a sector-2 track and the
# deficit days carry L_e of the total L, the observed yield is Y (1 - f (1 - eps) L_e/L), so
#   C_sec = 1 / (1 - f (1 - eps) L_e / L)                         (exact; nominal)
#   C_sec' = 1 + f (1/eps - 1) L_e / L                            (the spec's first-order form, printed for comparison)
# f from the simulation (spec), eps from data (double ratio against the rows without a sector-2 track).
# Variants: loss summed over sector patterns B0..B3 or over the number of sector-2 tracks (signal mix differs from the
# inclusive rows); f from the data of the late days (N - CB share in the window); deficit span 034-044 (spec) instead of
# 034-045. Per period: only period I has deficit days.
# Usage: python3 -I plots/all_report/wp2_csec.py
# ============================================================================
import math

TAB = "/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all"

Lday = {}
for line in open(TAB + "/wp2_lumi_perday.txt"):
    if line.startswith("#") or line.startswith("day"): continue
    t = line.strip().split(",")
    Lday[int(t[0])] = float(t[5])
L = sum(Lday.values())
L_I = sum(v for d, v in Lday.items() if 34 <= d <= 46)
span = {0: (34, 45), 1: (34, 44)}
Le = {q: sum(v for d, v in Lday.items() if span[q][0] <= d <= span[q][1]) for q in span}

eps = {}      # (q, cat) -> (eps, err)
fdata = {}    # window -> late N-CB share with a sector-2 track
for line in open(TAB + "/wp2_sector_eps.txt"):
    if line.startswith("tm|S0|0|N|"):
        key, a, b, c = line.strip().rsplit(",", 3)   # category names may contain commas
        p = key.split("|")
        v = [float(a), float(b), float(c)]
        eps[(int(p[4]), p[5])] = (v[0], v[1])
    elif line.startswith("S0|") and "|NmCB|045-067" in line:
        key, vals = line.strip().split(",", 1)
        w = int(key.split("|")[1])
        fdata[w] = float(vals.split(",")[-1])

fs = {}
for line in open(TAB + "/wp2_simfs.txt"):
    if not line.startswith("fs|"): continue
    key, a, b, c = line.strip().rsplit(",", 3)
    p = key.split("|")
    fs[(p[1], p[2], p[3])] = (float(a), float(b))

CAT_S2 = "uses s2"
PATS = ["B0 s2 only", "B1 s2 + 1 other", "B2 s2 + 2 others", "B3 s2 + 3 others"]
PATS_SIM = ["B0", "B1", "B2", "B3"]
N2 = ["one s2 track", ">=2 s2 tracks"]
signals = [("bump C2+C3 in B", "C2+C3", "B", 1), ("peak C4+C5 in P", "C4+C5", "P", 2), ("C5 in P", "C5", "P", 2),
           ("C4 in P", "C4", "P", 2), ("C2 in B", "C2", "B", 1), ("C3 in B", "C3", "B", 1),
           ("T12 in CR-U", "T12", "CR-U", 3), ("T12 all M4", "T12", "all M4", 0)]

def C_exact(loss, le, l): return 1.0 / (1.0 - loss * le / l)
def C_spec(f, e, le, l): return 1.0 + f * (1.0 / e - 1.0) * le / l

out = open(TAB + "/wp2_csec.txt", "w")
summ = open(TAB + "/wp2_csec_summary.txt", "w")
out.write("# wp2_csec.py: C_sec = 1/(1 - f (1 - eps) L_e/L); f = simulated share with a sector-2 track (S0, V3, mode A),\n")
out.write("# eps = data double ratio (S0, N, all M4) of rows with a sector-2 track, deficit days vs late days\n")
out.write("# L = %.3f; nominal span 034-045 (L_e %.3f, L_e/L %.4f); spec span 034-044 (L_e %.3f, L_e/L %.4f); L_I = %.3f\n"
          % (L, Le[0], Le[0] / L, Le[1], Le[1] / L, L_I))
out.write("# eps(uses s2): span 034-045 vs 046-067 = %.4f +- %.4f; span 034-044 vs 047-067 = %.4f +- %.4f\n\n"
          % (eps[(0, CAT_S2)][0], eps[(0, CAT_S2)][1], eps[(1, CAT_S2)][0], eps[(1, CAT_S2)][1]))
summ.write("# wp2_csec.py summary: signal, C_sec (all periods), stat, syst, total, C_sec(period I), total error (I)\n")
for name, g, w, wi in signals:
    f, ef = fs[(g, w, CAT_S2)]
    e0, ee0 = eps[(0, CAT_S2)]
    loss1 = f * (1 - e0)
    # stat of the loss: f and eps errors
    sl = math.hypot(ef * (1 - e0), f * ee0)
    # pattern decomposition
    loss2 = sum(fs[(g, w, ps)][0] * (1 - eps[(0, pd)][0]) for ps, pd in zip(PATS_SIM, PATS))
    loss3 = fs[(g, w, "one s2 track")][0] * (1 - eps[(0, N2[0])][0]) + fs[(g, w, ">=2 s2 tracks")][0] * (1 - eps[(0, N2[1])][0])
    # data f variant
    loss4 = fdata.get(wi, f) * (1 - e0)
    # spec span 034-044
    e1 = eps[(1, CAT_S2)][0]
    loss5 = f * (1 - e1)
    C1 = C_exact(loss1, Le[0], L)
    Cs = C_spec(f, e0, Le[0], L)
    C2 = C_exact(loss2, Le[0], L); C3 = C_exact(loss3, Le[0], L); C4 = C_exact(loss4, Le[0], L)
    C5 = C_exact(loss5, Le[1], L)
    stat = C1 * C1 * sl * Le[0] / L
    syst_dec = max(abs(C2 - C1), abs(C3 - C1))
    syst_f = abs(C4 - C1)
    syst_span = abs(C5 - C1)
    syst = math.sqrt(syst_dec ** 2 + syst_f ** 2 + syst_span ** 2)
    tot = math.hypot(stat, syst)
    # period I
    CI = C_exact(loss1, Le[0], L_I)
    totI = tot * (CI - 1) / (C1 - 1) if C1 > 1 else 0.0
    out.write("[%s]\n" % name)
    out.write("  f(sim, sector 2) = %.4f +- %.4f   (f in the other sectors: %s)\n" % (f, ef, " ".join(
        "%.3f" % fs[(g, w, "uses s%d" % s)][0] for s in (0, 1, 3, 4, 5))))
    out.write("  f(data, late days, N - CB in the window) = %s\n" % ("%.4f" % fdata[wi] if wi in fdata else "-"))
    out.write("  pattern shares (sim) B0..B3: %s ; one / >=2 sector-2 tracks: %.4f / %.4f\n" % (
        " ".join("%.4f" % fs[(g, w, ps)][0] for ps in PATS_SIM), fs[(g, w, "one s2 track")][0], fs[(g, w, ">=2 s2 tracks")][0]))
    out.write("  loss per unit L_e/L: simple f(1-eps) %.5f +- %.5f ; by pattern %.5f ; by number of s2 tracks %.5f ; f from data %.5f\n" % (
        loss1, sl, loss2, loss3, loss4))
    out.write("  C_sec nominal (exact, 034-045)      = %.5f +- %.5f (stat) +- %.5f (syst) = +- %.5f\n" % (C1, stat, syst, tot))
    out.write("      spec first-order form          = %.5f\n" % Cs)
    out.write("      by pattern / by n(s2) / f data = %.5f / %.5f / %.5f\n" % (C2, C3, C4))
    out.write("      spec span 034-044 (L_e %.1f)  = %.5f\n" % (Le[1], C5))
    out.write("      syst parts: decomposition %.5f, f data-sim %.5f, span %.5f\n" % (syst_dec, syst_f, syst_span))
    out.write("  period I only (L_I %.1f):  C_sec(I) = %.5f +- %.5f ; periods IIA, IIB, III: 1\n\n" % (L_I, CI, totI))
    summ.write("%s,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n" % (name, C1, stat, syst, tot, CI, totI))
# ratio peak/bump
out.write("[ratios] C_sec(peak C4+C5)/C_sec(bump C2+C3) and C5/C2 (enter DR and the golden ratio): see the values above;\n")
out.close(); summ.close()
print(open(TAB + "/wp2_csec.txt").read())
