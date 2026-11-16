#!/usr/bin/env python3
# ============================================================================
# x17p26_tables.py  --  X17 limit with PDG 2026 (analyst ruling D25, 15 November 2026). Prefix x17p26. Written 2026-11-15.
# Writes NEW tables (suffix _pdg26) next to the frozen ones; no existing file is changed.
# Inputs (read only):
#   /hdd3/przygoda/HADES/FAT_out/all_report/x17p26/x17p26_state.txt  (x17p26_limit.C: frozen code rerun with ONLY BR_dir changed:
#        r -> BR with 2.76e-4 and the BR_dir part of the X17 yield nuisance 0.11/2.68 -> 0.16/2.76 (sX1 0.2467 -> 0.2501);
#        normalisations (2), (3) with the WP5 eta scales rescaled to PDG 2026: a_eta / 1.00436, a_P / 1.0224)
#   docs/latex/tables_all/x17d_scan.txt, x17c_expected.txt, x17c_expected_asimov.txt, x17d_D1.txt, x17d_checks.txt (frozen, PDG 2020)
# Conversion of the numbers that were NOT rerun (toys, FC, template replicas, folding, counting, TRolke):
#   normalisation (1): x 2.76/2.68 = 1.029851 x k, k = rerun/frozen ratio at that mass (observed limit for data toys, FC and replicas;
#                      Asimov median for the expected toys; 1.0011-1.0014 everywhere);
#   normalisation (2), counting check, TRolke: x 1.00436 (a_eta of WP5 fitted with PDG 2020 BR(eta -> pi+pi-pi0) in C2/C3; v8);
#   normalisation (3): x 1.0224 (a_P = Y_peak / SMASH peak C4 + C5; v8).
# Run: python3 -I plots/all_report/x17p26_tables.py
# ============================================================================
import math

TAB = "/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all/"
ST = "/hdd3/przygoda/HADES/FAT_out/all_report/x17p26/x17p26_state.txt"
FBR = 2.76 / 2.68
F2, F3 = 1.00436, 1.0224

# ---------------------------------------------------------------- rerun state (PDG 2026)
S = {}
for l in open(ST):
    if l.startswith("#") or not l.strip(): continue
    p = [float(x) for x in l.split()]
    m = round(p[0], 2)
    S[m] = dict(r20=p[1], r26=p[2], ul=p[3:6], rh=p[6:9], re=p[9:12], c2=p[12], e1=p[13:18], e2=p[18:23], e3=p[23:28], ul20=p[28:31],
                c220=p[31], e120=p[32:37])
def kobs(m): return S[m]["r26"] / S[m]["r20"]
def kasi(m): return S[m]["e1"][2] / (S[m]["e120"][2] * FBR)
def F1obs(m): return FBR * kobs(m)
def F1asi(m): return FBR * kasi(m)

# ---------------------------------------------------------------- frozen tables (PDG 2020)
scan = {}
for l in open(TAB + "x17d_scan.txt"):
    if l.startswith("#") or not l.strip(): continue
    p = [float(x) for x in l.replace("|", " ").split()]
    scan[round(p[0], 2)] = p
asim = {}
for l in open(TAB + "x17c_expected_asimov.txt"):
    if l.startswith("#") or not l.strip(): continue
    p = [float(x) for x in l.replace("|", " ").split()]
    asim[round(p[0], 2)] = p
exp = {}   # (m, norm, method) -> [n, q0..q4]
for l in open(TAB + "x17c_expected.txt"):
    p = l.split()
    if len(p) >= 9 and p[2] in ("Asimov", "toys"):
        exp[(round(float(p[0]), 2), int(p[1]), p[2])] = [p[3]] + [float(x) for x in p[4:9]]
tcl = {}; mcs = {}; fc = []
for l in open(TAB + "x17d_checks.txt"):
    s = l.split()
    if l.startswith("  m ") and "UL asym" in l:
        tcl[(round(float(s[1]), 2), int(s[3].rstrip(":")))] = (float(s[6]), float(s[8]), float(s[10]))
    elif l.startswith("MC m"):
        m, n = round(float(s[2]), 2), int(s[4].rstrip(":"))
        mcs[(m, n)] = dict(nom=float(s[7].rstrip(";")), cbmean=float(s[13]), cbrms=float(s[15].rstrip(",")), cbr=float(s[18].rstrip(";")),
                           c4mean=float(s[24]), c4rms=float(s[26].rstrip(",")), c4r=float(s[29].rstrip(";")), comb=float(s[32]))
    elif l.startswith("  0.") and "/" in l:
        fc.append([float(x) for x in s[:7]])
    elif l.startswith("FC 90"):
        fcint = (float(l.split("[")[1].split(",")[0]), float(l.split(",")[1].split("]")[0]))

out = []
def W(fn, lines):
    with open(TAB + fn, "w") as f: f.write("\n".join(lines) + "\n")
    out.append(fn)

HDR = ["# PDG 2026 (analyst ruling D25, 15 November 2026): BR(eta -> pi+pi-e+e-) = (2.76 +- 0.16) e-4 (S = 1.5) instead of the frozen PDG 2020",
       "# (2.68 +- 0.11) e-4. Normalisation (1): frozen code rerun with ONLY BR_dir changed (r -> BR x 2.76e-4; BR_dir part of the X17 yield",
       "# nuisance 4.10 % -> 5.80 %, sX1 0.2467 -> 0.2501; x17p26_limit.C); r_UL changes by +0.11..+0.14 %. (2): x 1.00436, (3): x 1.0224 (WP5 eta",
       "# scales a_eta, a_P fitted with PDG 2020 BRs; rerun with the rescaled scales). Toy-based numbers: frozen toys x 2.76/2.68 x k (k = rerun ratio).",
       "# Written by plots/all_report/x17p26_tables.py from the frozen PDG 2020 table of the same name without _pdg26."]

# ---------------------------------------------------------------- scan (38 masses)
L = HDR + ["# x17d_scan_pdg26.txt -- OBSERVED 90 % CLs upper limits on BR(eta -> pi+pi-X17) BR(X17 -> e+e-), PDG 2026, same columns as x17d_scan.txt.",
           "# expA = Asimov median of the OBSERVED r = 0 fit x 2.76/2.68 x k (information); Z, p0, N_X rows, N5 copied (independent of the PDG edition).",
           "# m[MeV]   R_X   UL(1)      UL(2)      UL(3)     | r_hat(1)    err      | Z(1)   p0(1)   | expA(1) med | N_X UL(1) rows | UL(1) cos2 | N5_hat | r_UL(1)"]
for m in sorted(S):
    p, s = scan[m], S[m]
    L.append("%6.2f %7.4f %11.4e %11.4e %11.4e | %+11.4e %10.4e | %5.2f %7.4f | %11.4e | %8.1f | %11.4e | %7.1f | %.5f" % (
        m, p[1], s["ul"][0], s["ul"][1], s["ul"][2], s["rh"][0], s["re"][0], p[7], p[8], p[9] * F1obs(m), p[10], s["c2"], p[12], s["r26"]))
W("x17d_scan_pdg26.txt", L)

# selection of tab:x17scan (report_all), LaTeX rows
sel = [10, 12, 15, 16.7, 16.85, 17.0, 18, 20, 22, 25, 30, 35, 40]
L = ["% x17d_scan_pdg26.tex -- rows of Table tab:x17scan with PDG 2026 (x17p26_tables.py). UL (1), r_hat: BR_dir = 2.76e-4 (rerun with the PDG 2026",
     "% BR_dir uncertainty in s_X); (2) x 1.00436, (3) x 1.0224; frozen median: Asimov median of the frozen background-only set, rerun; cos2: rerun.",
     "\\begin{tabular}{@{}r r r r r r r r@{}}", "\\toprule",
     "$m_X$ [\\MeVcc] & UL (1) & UL (2) & UL (3) & $\\hat r\\times2.76\\times10^{-4}$ & local $Z$ & frozen median & UL (1), $1+\\cos^2$ \\\\", "\\midrule"]
for i, m in enumerate(sel):
    m = round(float(m), 2); s = S[m]; z = scan[m][7]
    zs = "\\textbf{%.2f}" % z if abs(m - 22) < 1e-6 else ("0" if z == 0 else "%.2f" % z)
    mm = ("%g" % m)
    if i == 0:
        L.append("%s & $%.2f\\times10^{-5}$ & $%.2f\\times10^{-5}$ & $%.2f\\times10^{-5}$ & $(%+.2f\\pm%.2f)\\times10^{-5}$ & %s & $%.2f\\times10^{-5}$ & $%.2f\\times10^{-5}$ \\\\" % (
            mm, 1e5 * s["ul"][0], 1e5 * s["ul"][1], 1e5 * s["ul"][2], 1e5 * s["rh"][0], 1e5 * s["re"][0], zs, 1e5 * s["e1"][2], 1e5 * s["c2"]))
    else:
        L.append("%s & %.2f & %.2f & %.2f & $%+.2f\\pm%.2f$ & %s & %.2f & %.2f \\\\" % (
            mm, 1e5 * s["ul"][0], 1e5 * s["ul"][1], 1e5 * s["ul"][2], 1e5 * s["rh"][0], 1e5 * s["re"][0], zs, 1e5 * s["e1"][2], 1e5 * s["c2"]))
L += ["\\bottomrule", "\\end{tabular}"]
W("x17d_scan_pdg26.tex", L)

# ---------------------------------------------------------------- frozen expectation: Asimov per mass
L = HDR + ["# x17c_expected_asimov_pdg26.txt -- EXPECTED 90 % CLs limits (Asimov, frozen background-only set of the data side bands), PDG 2026.",
           "# (1) rerun with BR_dir = 2.76e-4 and its PDG 2026 uncertainty in sX; (2), (3) rerun with a_eta / 1.00436, a_P / 1.0224. Same columns as x17c_expected_asimov.txt.",
           "# m[MeV]   R_X    (1): -2s        -1s        median     +1s        +2s     | (2) median  | (3) median | X17 rows at (1) median | s_X2D [rows/BR, a=1]"]
for m in sorted(S):
    a, s = asim[m], S[m]
    L.append("%6.2f %7.4f  %11.4e %11.4e %11.4e %11.4e %11.4e | %11.4e | %11.4e | %8.1f | %.4e" % (m, a[1], *s["e1"], s["e2"][2], s["e3"][2], a[9], a[10]))
W("x17c_expected_asimov_pdg26.txt", L)

# ---------------------------------------------------------------- frozen expectation: Asimov + toys at the quoted masses
L = HDR + ["# x17c_expected_pdg26.txt -- expected 90 % CLs limits of the FROZEN procedure in PDG 2026 units; same layout as x17c_expected.txt.",
           "# Asimov rows: rerun (x17p26_limit.C). Toy rows: frozen toy quantiles x 2.76/2.68 x k_Asimov(m) in (1), x 1.00436 in (2), x 1.0224 in (3).",
           "m      norm method   ntoys        -2s        -1s     median        +1s        +2s   toy/Asimov median"]
expc = {}
for m in (16.7, 16.85, 17.0, 10.0, 25.0, 40.0):
    for n in (1, 2, 3):
        e = S[m]["e%d" % n]
        t = exp[(m, n, "toys")]
        f = F1asi(m) if n == 1 else (F2 if n == 2 else F3)
        tq = [x * f for x in t[1:]]
        expc[(m, n)] = (e, tq, t[0])
        L.append("%6.2f %4d %-7s %6s %11.4e %11.4e %11.4e %11.4e %11.4e" % (m, n, "Asimov", "-", *e))
        L.append("%6.2f %4d %-7s %6s %11.4e %11.4e %11.4e %11.4e %11.4e   %.3f" % (m, n, "toys", t[0], *tq, tq[2] / e[2]))
W("x17c_expected_pdg26.txt", L)
L = ["% x17c_expected_pdg26.tex -- x17c_expected.tex with PDG 2026 (x17p26_tables.py): Asimov rerun with BR_dir = 2.76e-4 and its PDG 2026",
     "% uncertainty; toys x 2.76/2.68 x k_Asimov(m) (1.0011-1.0013, rerun ratio); (2) x 1.00436, (3) x 1.0224. Expected 90 % CLs UL x 1e5.",
     "\\begin{tabular}{r ccc cccr cc}", "\\toprule",
     " & \\multicolumn{3}{c}{(1) Asimov} & \\multicolumn{4}{c}{(1) background-only toys} & (2) & (3) \\\\",
     "$m_X$ [MeV/$c^2$] & median & $\\pm1\\sigma$ & $\\pm2\\sigma$ & median & $\\pm1\\sigma$ & $\\pm2\\sigma$ & toys & median & median \\\\", "\\midrule"]
for m in (16.7, 16.85, 17.0, 10.0, 25.0, 40.0):
    e, tq, nt = expc[(m, 1)]
    L.append("%.2f & %.2f & [%.2f, %.2f] & [%.2f, %.2f] & %.2f & [%.2f, %.2f] & [%.2f, %.2f] & %s & %.2f & %.2f \\\\" % (
        m, 1e5 * e[2], 1e5 * e[1], 1e5 * e[3], 1e5 * e[0], 1e5 * e[4], 1e5 * tq[2], 1e5 * tq[1], 1e5 * tq[3], 1e5 * tq[0], 1e5 * tq[4], nt,
        1e5 * expc[(m, 2)][0][2], 1e5 * expc[(m, 3)][0][2]))
L += ["\\bottomrule", "\\end{tabular}"]
W("x17c_expected_pdg26.tex", L)

# ---------------------------------------------------------------- Table D1
D1 = {}
for l in open(TAB + "x17d_D1.txt"):
    if l.startswith("#") or not l.strip() or not l[0].isspace() and not l[0].isdigit(): continue
    s = l.split()
    try:
        m, n = round(float(s[0]), 2), int(s[1])
    except ValueError:
        continue
    D1[(m, n)] = s
L = HDR + ["# x17d_D1_pdg26.txt -- Table D1 in PDG 2026 units (same columns as x17d_D1.txt).",
           "# UL_asym, r_hat: rerun; UL_toy, MC_RMS: frozen x factor; exp: frozen toys x factor (x17c_expected_pdg26.txt); N_X_UL rows unchanged;",
           "# sBRBR column dropped (superseded: use x17_fold_pdg26.txt).",
           "m      norm    UL_asym     UL_toy  toy/as      exp median [16%, 84%]     r_hat(BR)            Z     p0   N_X_UL    MC_RMS    UL_cos2"]
d1rows = []
for m in (16.7, 16.85, 17.0):
    s = S[m]
    for n in (1, 2, 3):
        a20, t20, ratio = tcl[(m, n)]
        f = F1obs(m) if n == 1 else (F2 if n == 2 else F3)
        ulA = s["ul"][n - 1]; ulT = t20 * f
        e, tq, nt = expc[(m, n)]
        mc = mcs[(m, n)]["comb"] * f
        z, p0, nx = float(D1[(m, n)][9]), float(D1[(m, n)][10]), float(D1[(m, n)][11])
        c2 = s["c2"] if n == 1 else float("nan")
        d1rows.append((m, n, ulA, ulT, tq, s["rh"][n - 1], s["re"][n - 1], z, nx))
        L.append("%6.2f %4d %11.4e %11.4e %7.3f %11.4e [%.3e, %.3e] %+.3e+-%.3e %6.2f %6.3f %8.1f %9.3e %11.4e" % (
            m, n, ulA, ulT, ulT / ulA, tq[2], tq[1], tq[3], s["rh"][n - 1], s["re"][n - 1], z, p0, nx, mc, c2))
cnt_br, cnt_err = 1.9234e-05, 1.4996e-04   # frozen counting check (x17d_D1.txt: +1.92e-05 +- 1.50e-04)
L.append("counting X17c: excess +153.9 +- 1198.2 rows, BR = excess / (a_eta s1) = %+.2e +- %.2e (a_eta / 1.00436; frozen %+.2e +- %.2e)" % (cnt_br * F2, cnt_err * F2, cnt_br, cnt_err))
L.append("Feldman-Cousins 90 %% CL interval at 17 MeV/c2 (normalisation 1): [%.3e, %.3e] (frozen %.4e x 2.76/2.68 x %.5f)" % (fcint[0], fcint[1] * F1obs(17.0), fcint[1], kobs(17.0)))
L.append("simulation expectation of Part I (comparison only, no BR_dir; SMASH background with PDG 2026 changes it by < 0.5 %): 1.9e-5 / 5.7e-6 (unchanged)")
W("x17d_D1_pdg26.txt", L)
L = ["% x17d_D1_pdg26.tex -- Table D1 with PDG 2026 (x17p26_tables.py). Observed 90% CLs upper limits on BR(eta->pi+pi-X17)BR(X17->e+e-), frozen procedure;",
     "% (1) rerun with BR_dir = 2.76e-4 (central value and PDG 2026 uncertainty in s_X); (2) x 1.00436, (3) x 1.0224; toy columns converted.",
     "\\begin{tabular}{llcccccc}", "\\toprule",
     "$m_X$ [MeV/$c^2$] & norm. & observed UL & toy-CL$_s$ UL & expected (median, 68\\%) & $\\hat{r}\\times2.76\\cdot10^{-4}$ & local $Z$ & $N_X^{\\rm UL}$ \\\\", "\\midrule"]
for (m, n, ulA, ulT, tq, rh, re, z, nx) in d1rows:
    L.append("%.2f & (%d) & $%.2f\\cdot10^{-5}$ & $%.2f\\cdot10^{-5}$ & $%.2f\\,[%.2f, %.2f]\\cdot10^{-5}$ & $(%+.2f\\pm%.2f)\\cdot10^{-5}$ & %.2f & %.0f \\\\" % (
        m, n, 1e5 * ulA, 1e5 * ulT, 1e5 * tq[2], 1e5 * tq[1], 1e5 * tq[3], 1e5 * rh, 1e5 * re, z, nx))
    if n == 3 and m < 16.99: L.append("\\addlinespace")
L += ["\\bottomrule", "\\end{tabular}"]
W("x17d_D1_pdg26.tex", L)

# ---------------------------------------------------------------- checks (toy CLs, FC, replicas) in PDG 2026 units
L = HDR + ["# x17d_checks_pdg26.txt -- x17d_checks.txt in PDG 2026 units (factors as in the header; nothing rerun except the asymptotic UL)."]
L.append("toy CLs vs asymptotic (data):")
for m in (16.7, 16.85, 17.0):
    for n in (1, 2, 3):
        a20, t20, ratio = tcl[(m, n)]
        f = F1obs(m) if n == 1 else (F2 if n == 2 else F3)
        L.append("  m %.2f norm %d: UL asym %.4e toy %.4e ratio %.4f   (frozen: %.4e / %.4e)" % (m, n, S[m]["ul"][n - 1], t20 * f, t20 * f / S[m]["ul"][n - 1], a20, t20))
L.append("FC at 17 MeV/c2 (norm 1): r, BR (PDG 2026 = r x 2.76e-4 x %.5f), t_obs, t68, t90, t95, p" % kobs(17.0))
for p in fc:
    L.append("  %.4f %.4e %.4f %.4f %.4f %.4f %.4f" % (p[0], p[0] * 2.76e-4 * kobs(17.0), p[2], p[3], p[4], p[5], p[6]))
L.append("FC 90 %% CL interval: [%.4e, %.4e]" % (fcint[0], fcint[1] * F1obs(17.0)))
for m in (16.7, 16.85, 17.0):
    for n in (1, 2, 3):
        c = mcs[(m, n)]
        f = F1obs(m) if n == 1 else (F2 if n == 2 else F3)
        L.append("MC m %.2f norm %d: nominal UL %.4e; CB replicas 200: UL mean %.4e RMS %.3e, r_hat RMS %.3e; C4 replicas 100: UL mean %.4e RMS %.3e, r_hat RMS %.3e; combined RMS %.3e (%.1f %%)" % (
            m, n, S[m]["ul"][n - 1], c["cbmean"] * f, c["cbrms"] * f, c["cbr"] * (FBR if n == 1 else f), c["c4mean"] * f, c["c4rms"] * f, c["c4r"] * (FBR if n == 1 else f),
            c["comb"] * f, 100 * c["comb"] * f / S[m]["ul"][n - 1]))
W("x17d_checks_pdg26.txt", L)

# ---------------------------------------------------------------- TRolke (no BR_dir; a_eta x 1/1.00436)
L = ["# trolke_x17_pdg26.txt -- trolke_x17.txt with PDG 2026: BR = N_X / (a_eta s1), a_eta(WP5) = 0.3823 / 1.00436 (PDG 2026 BR(eta -> pi+pi-pi0)",
     "# in the C2/C3 normalisation of the bump; v8); BR_dir does not enter. N_X intervals unchanged. (x17p26_tables.py)"]
for l in open(TAB + "trolke_x17.txt"):
    s = l.strip()
    if not s: continue
    hi = float(s.split("[")[-1].split(",")[1].rstrip("]"))
    L.append(s.split("-> BR")[0] + "-> BR x BR in [0.00e+00, %.2e]   (frozen PDG 2020: %.2e)" % (hi * F2, hi))
W("trolke_x17_pdg26.txt", L)

# ---------------------------------------------------------------- nuisance table
L = ["% x17c_nuisances_pdg26.tex -- x17c_nuisances.tex with the PDG 2026 BR_dir uncertainty in the X17 yield nuisance (1) (x17p26_tables.py).",
     "% Only the row 'X17 yield (1)' changes: BR_dir 0.11/2.68 = 4.1 % -> 0.16/2.76 = 5.8 %, total 24.7 % -> 25.0 %."]
for l in open(TAB + "x17c_nuisances.tex"):
    l = l.rstrip("\n")
    if l.startswith("%"): continue
    if l.startswith("X17 yield (1)"):
        l = "X17 yield (1), log-normal & $\\pm25.0$\\,\\% & kinematics, $N_5$ bias, direct shape, peak width, BR$_\\mathrm{dir}$ (PDG 2026: 5.8\\,\\%) \\\\"
    L.append(l)
W("x17c_nuisances_pdg26.tex", L)

# ---------------------------------------------------------------- master conversion list
m = 17.0
s17 = S[17.0]
L = ["# x17_pdg26_conversion.txt -- X17 numbers: frozen PDG 2020 value -> PDG 2026 value (x17p26_tables.py; see x17p26_limit.C).",
     "# factors: (1) x %.6f (2.76/2.68) x k, k(17.0) observed %.5f / Asimov %.5f (BR_dir uncertainty 4.1 -> 5.8 %% in s_X); (2) x %.5f; (3) x %.4f" % (FBR, kobs(17.0), kasi(17.0), F2, F3),
     "r_UL (1) at 16.70/16.85/17.00: %.5f/%.5f/%.5f -> %.5f/%.5f/%.5f" % (S[16.7]["r20"], S[16.85]["r20"], S[17.0]["r20"], S[16.7]["r26"], S[16.85]["r26"], S[17.0]["r26"]),
     "UL (1) 16.70/16.85/17.00: %.4e/%.4e/%.4e -> %.4e/%.4e/%.4e" % (S[16.7]["ul20"][0], S[16.85]["ul20"][0], S[17.0]["ul20"][0], S[16.7]["ul"][0], S[16.85]["ul"][0], S[17.0]["ul"][0]),
     "UL (2) 16.70/16.85/17.00: %.4e/%.4e/%.4e -> %.4e/%.4e/%.4e" % (S[16.7]["ul20"][1], S[16.85]["ul20"][1], S[17.0]["ul20"][1], S[16.7]["ul"][1], S[16.85]["ul"][1], S[17.0]["ul"][1]),
     "UL (3) 16.70/16.85/17.00: %.4e/%.4e/%.4e -> %.4e/%.4e/%.4e" % (S[16.7]["ul20"][2], S[16.85]["ul20"][2], S[17.0]["ul20"][2], S[16.7]["ul"][2], S[16.85]["ul"][2], S[17.0]["ul"][2]),
     "UL (1) 1+cos2 16.70/16.85/17.00: %.4e/%.4e/%.4e -> %.4e/%.4e/%.4e" % (S[16.7]["c220"], S[16.85]["c220"], S[17.0]["c220"], S[16.7]["c2"], S[16.85]["c2"], S[17.0]["c2"]),
     "frozen Asimov (1) at 17.0: %s -> %s" % (" ".join("%.4e" % x for x in s17["e120"]), " ".join("%.4e" % x for x in s17["e1"])),
     "frozen toys (1) at 17.0 [2.3, 15.9, 50, 84.1, 97.7 %%]: %s -> %s" % (" ".join("%.4e" % x for x in exp[(17.0, 1, "toys")][1:]), " ".join("%.4e" % x for x in expc[(17.0, 1)][1])),
     "toy CLs on the data (1) 16.70/16.85/17.00: %.4e/%.4e/%.4e -> %.4e/%.4e/%.4e" % (tcl[(16.7, 1)][1], tcl[(16.85, 1)][1], tcl[(17.0, 1)][1],
                                                                                tcl[(16.7, 1)][1] * F1obs(16.7), tcl[(16.85, 1)][1] * F1obs(16.85), tcl[(17.0, 1)][1] * F1obs(17.0)),
     "FC upper edge 17.0: %.4e -> %.4e" % (fcint[1], fcint[1] * F1obs(17.0)),
     "MC combined RMS (1) 17.0: %.3e -> %.3e; CB %.3e -> %.3e; C4 %.3e -> %.3e; C4 replica mean %.4e -> %.4e" % (
         mcs[(17.0, 1)]["comb"], mcs[(17.0, 1)]["comb"] * F1obs(17.0), mcs[(17.0, 1)]["cbrms"], mcs[(17.0, 1)]["cbrms"] * F1obs(17.0), mcs[(17.0, 1)]["c4rms"],
         mcs[(17.0, 1)]["c4rms"] * F1obs(17.0), mcs[(17.0, 1)]["c4mean"], mcs[(17.0, 1)]["c4mean"] * F1obs(17.0)),
     "r_hat (1) 17.0: %+.4e +- %.4e -> %+.4e +- %.4e" % (scan[17.0][5], scan[17.0][6], s17["rh"][0], s17["re"][0]),
     "whole mass range (1): %.3e - %.3e -> %.3e - %.3e" % (min(S[x]["ul20"][0] for x in S), max(S[x]["ul20"][0] for x in S), min(S[x]["ul"][0] for x in S), max(S[x]["ul"][0] for x in S)),
     "counting X17c: %+.3e +- %.3e -> %+.3e +- %.3e" % (cnt_br, cnt_err, cnt_br * F2, cnt_err * F2)]
W("x17_pdg26_conversion.txt", L)
print("x17p26_tables: wrote", ", ".join(out))
print("k_obs:", {m: round(kobs(m), 5) for m in (16.7, 16.85, 17.0, 10.0, 25.0, 40.0)})
print("k_asi:", {m: round(kasi(m), 5) for m in (16.7, 16.85, 17.0, 10.0, 25.0, 40.0)})
