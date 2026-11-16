#!/usr/bin/env python3
# ============================================================================
# wp2_tableN1.py  --  WP2: Table N1 (normalisation items: L, C_sec, f_win) from the WP2 tables
# ============================================================================
# Inputs (docs/latex/tables_all/): wp2_lumi_perperiod.txt, wp2_csec.txt, wp2_csec_summary.txt, wp2_win.txt,
#   wp2_win_summary.txt, wp2_vzfit.txt, wp2_wincomp.txt
# Outputs: docs/latex/tables_all/wp2_N1.txt (all numbers, machine lines "n1|key,value,...") and wp2_N1.tex (booktabs)
# Recommended f_win for the eta yields (a_eta, sigma(pp -> eta X)): N - CB share of the material in the windows B + P
# (S0, L-weighted over the three empty-target days, normaliser ch[-200,-170)). Uncertainties:
#   stat   = empty and full counts and both normalisers (event-group variances)
#   norm   = flux-ratio range: lower = vz-profile fit (wp2_vzfit.C), upper = positivity bound b <= F/E in [-170,-160)
#   gas    = one-sided (lower): on-axis share g of the in-cell empty rows x the in-cell part of f_win (wp2_wincomp.C)
#   comp   = symmetric: |f(B + P) - f(0.40-0.70)| (the eta share of the material yield may follow the whole fit range)
# Usage: python3 -I plots/all_report/wp2_tableN1.py
# ============================================================================
import math
import re

TAB = "/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all"
DAYS = ["034", "047", "065"]

# ---------------------------------------------------------------- luminosity
Lp = {}
for line in open(TAB + "/wp2_lumi_perperiod.txt"):
    if line.startswith("#") or line.startswith("group"):
        continue
    t = line.strip().split(",")
    Lp[t[0]] = (float(t[4]), float(t[3]), int(t[2]), float(t[5]))   # L6_an, L6_all, nrun_an, stat
Lw = [Lp["I"][0], Lp["IIA"][0] + Lp["IIB"][0], Lp["III"][0]]
Ltot = sum(Lw)
w = [x / Ltot for x in Lw]

# ---------------------------------------------------------------- C_sec
csec_txt = open(TAB + "/wp2_csec.txt").read()
m = re.search(r"L_e ([\d.]+), L_e/L ([\d.]+)\); spec span 034-044 \(L_e ([\d.]+)", csec_txt)
Le45, Le45f, Le44 = float(m.group(1)), float(m.group(2)), float(m.group(3))
m = re.search(r"span 034-045 vs 046-067 = ([\d.]+) \+- ([\d.]+); span 034-044 vs 047-067 = ([\d.]+) \+- ([\d.]+)", csec_txt)
eps45, eeps45, eps44, eeps44 = map(float, m.groups())
fsim = {}
for blk in csec_txt.split("\n[")[1:]:
    name = blk.split("]")[0]
    mm = re.search(r"f\(sim, sector 2\) = ([\d.]+) \+- ([\d.]+)", blk)
    if mm:
        fsim[name] = (float(mm.group(1)), float(mm.group(2)))
csec = {}
for line in open(TAB + "/wp2_csec_summary.txt"):
    if line.startswith("#"):
        continue
    t = line.strip().split(",")
    csec[t[0]] = [float(x) for x in t[1:]]   # C, stat, syst, tot, C_I, tot_I

# ---------------------------------------------------------------- f_win
summ = {}
for line in open(TAB + "/wp2_win_summary.txt"):
    if line.startswith("#"):
        continue
    key, rest = line.strip().split(",", 1)
    summ[key] = [float(x) for x in rest.split(",")]
# per-day f and stat (nominal normaliser) from wp2_win.txt
perday = {}
for line in open(TAB + "/wp2_win.txt"):
    if not line.startswith("  ") or "|" not in line or len(line) < 40:
        continue
    sel, win, qn, day = line[2:7].strip(), line[8:21].strip(), line[22:27].strip(), line[28:31]
    nums = re.findall(r"(-?[\d.]+) \+- +(-?[\d.]+)", line.split("|", 1)[1])
    if nums:
        perday[(sel, win, qn, day)] = [(float(a), float(b)) for a, b in nums]
# b per day: nominal (wincomp b_ch), b_max (wincomp), b_fit (vzfit nominal)
bch, bmax, bfit, ebfit = {}, {}, {}, {}
for line in open(TAB + "/wp2_wincomp.txt"):
    mm = re.match(r"\s+(\d{3})\s+b_ch ([\d.]+) \+- ([\d.]+)\s+b_max ([\d.]+) \+- ([\d.]+)", line)
    if mm:
        bch[mm.group(1)] = float(mm.group(2))
        bmax[mm.group(1)] = float(mm.group(4))
gas = {}
for line in open(TAB + "/wp2_wincomp.txt"):
    if line.startswith("wc|"):
        k, rest = line.strip()[3:].split(",", 1)
        v = [float(x) for x in rest.split(",")]
        gas[k] = v   # total, entrance, cell, exit, gas bound
    mm = re.search(r"g = .* = ([\d.]+)", line)
    if mm:
        g_share = float(mm.group(1))
for line in open(TAB + "/wp2_vzfit.txt"):
    mm = re.match(r"\s+nominal\s+(\d{3}) .*b\s+([\d.]+) \+-\s+([\d.]+)", line)
    if mm:
        bfit[mm.group(1)] = float(mm.group(2))
        ebfit[mm.group(1)] = float(mm.group(3))
bnorm = {}
for line in open(TAB + "/wp2_win.txt"):
    mm = re.match(r"\s+(\d{3}) \| ch ([\d.]+) ([\d.]+) ([\d.]+) \+- ([\d.]+)", line)
    if mm:
        bnorm[mm.group(1)] = (float(mm.group(4)), float(mm.group(5)))


def fwin(sel, win, qn="N-CB"):
    """combined f [%] for the nominal normaliser, stat, and the fit / positivity variants"""
    v = summ["%s|%s|%s" % (sel, win, qn)]
    f, st = v[0], v[1]
    fd = [perday[(sel, win, qn, d)][0] for d in DAYS]   # (f, err) with normaliser 0
    ffit = sum(w[i] * fd[i][0] * bfit[d] / bch[d] for i, d in enumerate(DAYS))
    fmax = sum(w[i] * fd[i][0] * bmax[d] / bch[d] for i, d in enumerate(DAYS))
    return dict(f=f, stat=st, fit=ffit, max=fmax, eenear=v[2], eefar=v[3], chw=v[4], chfar=v[5], days=fd)


out = []
P = out.append
P("# Table N1 (WP2): normalisation items for absolute quantities; numbers from the WP2 tables (see wp2_tableN1.py)")
P("# L = START-free elastic luminosity (chain N) of the analysed runs. It ALREADY contains: the PT2 DAQ live fraction and beam")
P("#   duty factor, the sector-2 PT2 (elastic) trigger deficit of days 034-044 (S_2, sector 2 replaced by rho x mean of the")
P("#   other sectors), and (1 - w6): target-window elastic events removed (LH2 only). It does NOT contain: the PT3 word")
P("#   factor 64/63, the PT3 live fraction (one-sided +7.3/-0 %, not applied), the sector-2 deficit of the PT3 channel")
P("#   rows (C_sec below), the material (empty-target) contribution of the PT3 channel rows (f_win below).")
P("")
P("[luminosity] analysed runs, nb^-1 (stat = Poisson of the elastic count)")
for k in ["I", "IIA", "IIB", "III", "all"]:
    P("  %-4s L = %9.3f +- %5.3f   (all runs %9.3f, analysed runs %d)" % (k, Lp[k][0], Lp[k][3], Lp[k][1], Lp[k][2]))
P("  deficit span 034-045: L_e = %.3f (L_e/L = %.4f); 034-044: %.3f; clean days 046-067: %.3f; 045-067: %.3f"
  % (Le45, Le45f, Le44, Lp["all"][0] - Le45, Lp["late_045_067"][0]))
P("  correlated L uncertainty: 5.56 % HADES syst (6.9 % total with sigma_int and the Argonne normalisation)")
P("n1|L,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f" % (Lp["all"][0], Lp["I"][0], Lp["IIA"][0], Lp["IIB"][0], Lp["III"][0], Le45,
                                             Lp["all"][0] - Le45))
P("")
P("[C_sec] sector-2 PT3 deficit: C_sec = 1/(1 - f (1 - eps) L_e/L); eps from data (S0 rows with a sector-2 track, double ratio")
P("  against rows without one, days 034-045 vs 046-067), f = simulated share with a sector-2 track (S0, V3)")
P("  eps = %.4f +- %.4f (034-045 vs 046-067); %.4f +- %.4f (spec span 034-044 vs 047-067)" % (eps45, eeps45, eps44, eeps44))
for k, lab in [("bump C2+C3 in B", "bump C2+C3 in B"), ("peak C4+C5 in P", "peak C4+C5 in P"), ("C2 in B", "C2 in B"),
               ("C3 in B", "C3 in B"), ("C4 in P", "C4 in P"), ("C5 in P", "C5 in P"), ("T12 in CR-U", "T12 in CR-U"),
               ("T12 all M4", "T12 all M4")]:
    c = csec[k]
    P("  %-18s f = %.4f +- %.4f   C_sec = %.4f +- %.4f (stat) +- %.4f (syst)   period I only: %.4f +- %.4f"
      % (lab, fsim[k][0], fsim[k][1], c[0], c[1], c[2], c[4], c[5]))
    P("n1|csec|%s,%.5f,%.5f,%.5f,%.5f,%.5f" % (lab, c[0], c[1], c[2], c[4], c[5]))
cb, cp = csec["bump C2+C3 in B"], csec["peak C4+C5 in P"]
rat = cp[0] / cb[0]
P("  ratio peak/bump C_sec = %.4f (enters DR); C5/C2 = %.4f (enters the golden ratio): both 1 within 0.2 %%"
  % (rat, csec["C5 in P"][0] / csec["C2 in B"][0]))
P("")
P("[f_win] material (empty-target) share of the full-target rows, N - CB, L-weighted over days 034 (I), 047 (IIA+IIB), 065 (III) [%]")
P("  normaliser b = U_full/U_empty, channel rows in [-200,-170) mm. Per day b: " +
  ", ".join("%s %.3f +- %.3f (fit %.3f +- %.3f, bound %.3f)" % (d, bnorm[d][0], bnorm[d][1], bfit[d], ebfit[d], bmax[d]) for d in DAYS))
P("  window      f_win +- stat | fit-b   bound-b | EpEm near / far (excluded) | per day 034 047 065 (stat)")
rows = {}
for win in ["all M4", "B", "P", "B+P", "0.40-0.70", "CR-U"]:
    r = fwin("S0", win)
    rows[win] = r
    P("  %-10s %6.2f +- %4.2f | %6.2f  %6.2f | %6.2f / %6.2f | %s" % (
        win, r["f"], r["stat"], r["fit"], r["max"], r["eenear"], r["eefar"],
        "  ".join("%5.2f +- %4.2f" % x for x in r["days"])))
    P("n1|fwin|S0|%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f" % (win, r["f"], r["stat"], r["fit"], r["max"], r["eenear"], r["eefar"],
                                                              r["days"][0][0], r["days"][2][0]))
P("  composition by vertex region (B+P): entrance %.2f, in-cell %.2f, exit %.2f; on-axis share of in-cell empty rows g = %.3f"
  % (gas["B+P"][1], gas["B+P"][2], gas["B+P"][3], g_share))
P("  selection variants, B+P (N - CB): " + "; ".join(
    "%s %.2f +- %.2f" % (s, summ["%s|B+P|N-CB" % s][0], summ["%s|B+P|N-CB" % s][1]) for s in ["S0", "pre", "S0b1", "vz150", "vz170", "dz40"]))
P("  selection variants, all M4 (N - CB): " + "; ".join(
    "%s %.2f +- %.2f" % (s, summ["%s|all M4|N-CB" % s][0], summ["%s|all M4|N-CB" % s][1]) for s in ["S0", "pre", "S0b1", "vz150", "vz170", "dz40"]))
P("")
# recommended value
r = rows["B+P"]
f0, st = r["f"], r["stat"]
up_norm, dn_norm = r["max"] - f0, f0 - r["fit"]
dn_gas = gas["B+P"][4]
comp = abs(f0 - rows["0.40-0.70"]["f"])
up = math.sqrt(st ** 2 + up_norm ** 2 + comp ** 2)
dn = math.sqrt(st ** 2 + dn_norm ** 2 + dn_gas ** 2 + comp ** 2)
P("[recommended] f_win for the eta yields (S0, isBest >= 0): %.1f %%  +%.1f / -%.1f  (stat %.1f, normaliser +%.1f/-%.1f, gas -%.1f, composition %.1f)"
  % (f0, up, dn, st, up_norm, dn_norm, dn_gas, comp))
one = 1 - f0 / 100
P("  1 - f_win = %.3f +%.3f / -%.3f   (relative +%.1f / -%.1f %%)" % (one, dn / 100, up / 100, 100 * dn / 100 / one, 100 * up / 100 / one))
cs = 0.5 * (cb[0] + cp[0])
ecs = max(cb[3], cp[3])
K = cs * one
eKu = K * math.hypot(ecs / cs, dn / 100 / one)
eKd = K * math.hypot(ecs / cs, up / 100 / one)
P("  C_sec (eta, mean of bump and peak) = %.4f +- %.4f;  C_sec x (1 - f_win) = %.3f +%.3f / -%.3f" % (cs, ecs, K, eKu, eKd))
P("  sigma(pp -> eta X) = a_eta x 1.0665 mb x C_sec x (1 - f_win) (spec WP8); a_eta and its bias treatment are WP5/WP8 items")
P("n1|recommended,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.5f,%.5f,%.4f,%.4f,%.4f" % (f0, up, dn, st, up_norm, dn_norm, dn_gas, comp, cs, ecs, K, eKu, eKd))
fB, fP = rows["B"], rows["P"]
rr = (1 - fP["f"] / 100) / (1 - fB["f"] / 100)
err = rr * math.hypot(fP["stat"] / 100 / (1 - fP["f"] / 100), fB["stat"] / 100 / (1 - fB["f"] / 100))
P("  ratios (DR, golden ratio): (1 - f_P)/(1 - f_B) = %.3f +- %.3f (stat) -> no correction; quote +-%.1f %% as a systematic" % (rr, err, 100 * max(abs(rr - 1), err)))
P("n1|ratioPB,%.4f,%.4f" % (rr, err))
P("")
P("[per period] factors for the WP7 stability test: K_p = C_sec,p x (1 - f_win,p); f_win,p = B+P value of the empty-target day of the period")
cI = 0.5 * (cb[4] + cp[4])
ecI = max(cb[5], cp[5])
fd = r["days"]
for p, (lab, Lval, c, ec, fi) in enumerate([("I", Lp["I"][0], cI, ecI, fd[0]), ("IIA", Lp["IIA"][0], 1.0, 0.0, fd[1]),
                                              ("IIB", Lp["IIB"][0], 1.0, 0.0, fd[1]), ("III", Lp["III"][0], 1.0, 0.0, fd[2])]):
    Kp = c * (1 - fi[0] / 100)
    eK = Kp * math.hypot(ec / c, fi[1] / 100 / (1 - fi[0] / 100))
    P("  %-4s L = %8.3f  C_sec = %.4f +- %.4f  f_win = %5.2f +- %4.2f %%  K = %.4f +- %.4f" % (lab, Lval, c, ec, fi[0], fi[1], Kp, eK))
    P("n1|period|%s,%.3f,%.5f,%.5f,%.4f,%.4f,%.5f,%.5f" % (lab, Lval, c, ec, fi[0], fi[1], Kp, eK))
P("  (IIA and IIB share the empty-target day 047, so their f_win are fully correlated; the f_win differences between the days")
P("   are within about 2 sigma of the per-day statistics: test with and without per-period f_win)")
P("")
P("[fall-backs] C_sec determinable: yes. f_win determinable: yes (empty-target runs exist: 7 + 16 + 104 runs, 1270 DSTs);")
P("  the spec fall-back range 0.7-2.1 % (elastic w6 x A-dependence) is NOT used: the measured share is about 8x larger.")
txt = "\n".join(out) + "\n"
open(TAB + "/wp2_N1.txt", "w").write(txt)
print(txt)

# ---------------------------------------------------------------- LaTeX snippet
T = []
T.append(r"% Table N1 (WP2), generated by plots/all_report/wp2_tableN1.py from docs/latex/tables_all/wp2_*.txt")
T.append(r"\begin{table}[htbp]\centering\small")
T.append(r"\caption{Normalisation items for absolute quantities (WP2). $L$ = START-free elastic luminosity of the analysed runs; it already"
         r" contains the sector-2 PT2 deficit removal and the elastic target-window removal ($w_6$). $C_\mathrm{sec}$ corrects the"
         r" sector-2 PT3 trigger deficit of days 034--045 for the channel rows; $f_\mathrm{win}$ is the material (empty-target) share"
         r" of $N-\mathrm{CB}$ (S0, isBest $\geq 0$), measured in the empty-target runs of days 034, 047 and 065 and combined with"
         r" the luminosity of the period each day represents.}")
T.append(r"\label{tab:N1}")
T.append(r"\begin{tabular}{llr}")
T.append(r"\toprule")
T.append(r"item & selection / window & value \\")
T.append(r"\midrule")
T.append(r"$L$ [nb$^{-1}$] & all / I / IIA / IIB / III & %.1f / %.1f / %.1f / %.1f / %.1f \\" % (
    Lp["all"][0], Lp["I"][0], Lp["IIA"][0], Lp["IIB"][0], Lp["III"][0]))
T.append(r"$L$ [nb$^{-1}$] & deficit days 034--045 / clean 046--067 & %.1f / %.1f \\" % (Le45, Lp["all"][0] - Le45))
T.append(r"$\varepsilon_{s2}$ & rows with a sector-2 track, 034--045 & $%.4f \pm %.4f$ \\" % (eps45, eeps45))
T.append(r"$f_{s2}$ (sim) & bump C2+C3 / peak C4+C5 & %.3f / %.3f \\" % (fsim["bump C2+C3 in B"][0], fsim["peak C4+C5 in P"][0]))
T.append(r"$C_\mathrm{sec}$ & bump C2+C3 in B & $%.4f \pm %.4f$ \\" % (cb[0], cb[3]))
T.append(r"$C_\mathrm{sec}$ & peak C4+C5 in P & $%.4f \pm %.4f$ \\" % (cp[0], cp[3]))
T.append(r"$C_\mathrm{sec}$ & period I only (bump / peak) & %.4f / %.4f \\" % (cb[4], cp[4]))
T.append(r"\midrule")
for win, lab in [("all M4", r"all $M_4$"), ("B", "bump B"), ("P", "peak P"), ("B+P", "B + P"), ("0.40-0.70", r"0.40--0.70"), ("CR-U", "CR-U")]:
    rw = rows[win]
    T.append(r"$f_\mathrm{win}$ [\%%] & %s & $%.1f \pm %.1f$ (fit-$b$ %.1f, bound %.1f) \\" % (lab, rw["f"], rw["stat"], rw["fit"], rw["max"]))
T.append(r"$f_\mathrm{win}$ [\%%] & recommended for $\eta$ (B + P) & $%.1f^{+%.1f}_{-%.1f}$ \\" % (f0, up, dn))
T.append(r"$C_\mathrm{sec}\,(1-f_\mathrm{win})$ & $\sigma(pp\to\eta X)$ & $%.3f^{+%.3f}_{-%.3f}$ \\" % (K, eKu, eKd))
T.append(r"$(1-f_P)/(1-f_B)$ & ratios (DR, golden) & $%.3f \pm %.3f$ (not applied) \\" % (rr, err))
T.append(r"\bottomrule")
T.append(r"\end{tabular}")
T.append(r"\end{table}")
open(TAB + "/wp2_N1.tex", "w").write("\n".join(T) + "\n")
