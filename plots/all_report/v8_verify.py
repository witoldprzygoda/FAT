#!/usr/bin/env python3
# ============================================================================================================================
# v8_verify.py -- independent verifier of the WP8 physics quantities (prefix v8), 2026-11-14.
# Own arithmetic from the PRIMARY inputs (WP5 fit parameters, WP6 estimator yields and SMASH window sums, WP7 family deltas);
# no WP8 code is used. Outputs: docs/latex/tables_all/v8_verification.txt, v8_tableP2_corr.tex, v8_literature.tex;
# /hdd3/.../all_report/v8/v8_verify.log (stdout).  Run: python3 -I plots/all_report/v8_verify.py
# External numbers typed here, with their (verified) source:
#   PDG 2020 (used by the simulation, sim_eta_physics.h kBRppee/kBRppg, A2 note): eta->pi+pi-pi0 22.92+-0.28 %, eta->pi+pi-gamma
#     4.22+-0.08 %, eta->pi+pi-e+e- 2.68+-0.11 e-4, pi0->e+e-gamma 1.174+-0.035 % (P.A. Zyla et al., PTEP 2020, 083C01)
#   PDG 2026 (current edition, listing created 2026-06-01; F. Takahashi et al. (PDG), Int. J. Mod. Phys. A 41, 2630011 (2026)):
#     23.02+-0.25 %, 4.28+-0.07 %, 2.76+-0.16 e-4 (S=1.5), 1.174+-0.035 %  (pdg.lbl.gov/2026/listings/rpp2026-list-eta.pdf, -pi-zero.pdf)
#   SMASH eta per collision of the simulated sample: N_eta = 2.364e8 in 1e10 events, sigma = 2.364e8 x 45.239 mb / 1e10 = 1.0694 mb
#     (A2 note, sect. 4; the sim study's 1.069 mb); 1.0665 +- 0.0016 mb is a 200-file sub-sample estimate.
#   C_sec(peak C4+C5 in P) = 1.0116 +- 0.0029 (WP2 Table N1), f_eta = 3 +4/-3 % (AMENDMENT 1).
# ============================================================================================================================
import math, os, re

AR = "/hdd3/przygoda/HADES/FAT_out/all_report"
TAB = "/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all"
L = 6079.2
YF = 64.0 / 63.0
OUT = []
def P(s=""):
    print(s); OUT.append(s)

def lines(path):
    with open(path) as f:
        return f.read().splitlines()

# ---------------------------------------------------------------- primary inputs
w5 = lines(os.path.join(AR, "wp5/results.txt"))
w6 = lines(os.path.join(AR, "wp6/results.txt"))
w7 = lines(os.path.join(AR, "wp7/syst.txt"))
def d5(block, name):
    for l in w5:
        a = l.split()
        if len(a) > 3 and a[0] == "derived" and a[1] == block and a[2] == name:
            d = {"v": float(a[3])}
            for i in range(4, len(a) - 1, 2):
                try: d[a[i]] = float(a[i + 1])
                except ValueError: pass
            return d
    raise KeyError(block + " " + name)
def e6(prefix):
    for l in w6:
        if l.startswith("est " + prefix):
            m = re.search(r"\)?\s([-\d.e+]+) stat ([-\d.e+]+) N ([-\d.e+]+) CB ([-\d.e+]+) A ([-\d.e+]+) F ([-\d.e+]+) mc ([-\d.e+]+)", l)
            return dict(v=float(m.group(1)), stat=float(m.group(2)), mc=float(m.group(7)))
    raise KeyError(prefix)
sims = {}
for l in w6:
    if l.startswith("simsmash"):
        a = l.split()[1:]
        sims = {a[i]: float(a[i + 1]) for i in range(0, len(a), 2)}
syst7 = {}
for l in w7:
    if l.startswith("syst "):
        a = l.split(); v = [float(x) for x in a[2:]]
        syst7[a[1]] = dict(v=v[0], stat=v[1], fam=v[2:14], syst=v[14], mc=v[18], nup=v[19], ndn=v[20])

a_eta = d5("FMAIN", "a_eta"); rho = d5("FMAIN", "rho_d"); kap = d5("FMAIN", "kappa")
al2 = d5("FGOLD", "alpha_2"); al5 = d5("FGOLD", "alpha_5")
YbF = d5("FMAIN", "Y_bump"); YpF = d5("FMAIN", "Y_peak")
Yp = e6("PSB_Ypeak "); Yb = e6("B4_Ybump_all_extrap"); YbE2 = e6("B4_Ybump_E2")
YpE1 = e6("PSB_Ypeak_E1"); YpE2 = e6("PSB_Ypeak_E2")

P("# v8_verification.txt -- plots/all_report/v8_verify.py (independent verifier v8, 2026-11-14)")
P("# S0, isBest >= 0, L = 6079.2 nb^-1; data yields = raw rows x 64/63; SMASH window sums at L (no 63/64) from wp6 'simsmash'.")
P("")
P("## A. Reproduction with own arithmetic (WP8 value in brackets)")
BR0 = dict(ppp0=0.2292, p0D=0.01174, ppg=0.0422, ppee=2.68e-4)          # PDG 2020 = simulation
BRE0 = dict(ppp0=0.0028, p0D=0.00035, ppg=0.0008, ppee=0.11e-4)
BR1 = dict(ppp0=0.2302, p0D=0.01174, ppg=0.0428, ppee=2.76e-4)          # PDG 2026
BRE1 = dict(ppp0=0.0025, p0D=0.00035, ppg=0.0007, ppee=0.16e-4)
SIG0, SIG1 = 1.0665, 2.364e8 * 45.239 / 1e10
CS0, CS1 = 1.0110, 1.0116
FE = 0.97

BR_G = BR0["ppee"] * al5["v"] / al2["v"]
BR_A = BR0["ppee"] * rho["v"]
peakS, bumpS = sims["peakP_all"], sims["bumpB_all"]
aP = Yp["v"] / peakS
aB = Yb["v"] / bumpS
aBE2 = YbE2["v"] / sims["bumpB_E2"]
DR = aP / aB
DRF = (YpF["v"] / YbF["v"]) / (peakS / bumpS)
sig = aP * SIG0 * CS0 * FE
sigF = a_eta["v"] * SIG0 * CS0 * FE
P("BR_G = 2.68e-4 x alpha_5/alpha_2 = 2.68e-4 x %.5f/%.5f = %.5e   [2.6546e-4]" % (al5["v"], al2["v"], BR_G))
P("BR_A = 2.68e-4 x rho_d = %.5e   [3.5247e-4]" % BR_A)
P("a_P = Y_peak/(C4+C5 SMASH in P) = %.2f/%.2f = %.6f   [0.418249]" % (Yp["v"], peakS, aP))
P("a_B(all) = Y_bump(all)/(C2+C3 SMASH in B) = %.1f/%.1f = %.6f ; a_B(E2) = %.6f [0.311677]" % (Yb["v"], bumpS, aB, aBE2))
P("DR = a_P/a_B(all) = %.5f   [1.39743];  DR_F = (%.1f/%.1f)/(%.1f/%.1f) = %.5f   [1.16896]" % (DR, YpF["v"], YbF["v"], peakS, bumpS, DRF))
P("sigma = a_P x 1.0665 x 1.0110 x 0.97 = %.5f mb   [0.43744];  F-main: %.5f mb [0.39982]" % (sig, sigF))
P("sigma/sigma_SMASH = %.4f [0.410];  visible peak %.5f nb [0.83695], bump %.4f nb [2.696]" % (sig / SIG0, Yp["v"] / L, Yb["v"] / L))
P("sigma x BR: ppee %.1f nb [117.2], ppg %.0f nb [18460], ppp0 (a_B E2) %.0f nb [74714]; F-main a rho: %.1f nb [140.9]" % (
    sig * 1e6 * BR0["ppee"], sig * 1e6 * BR0["ppg"], aBE2 * SIG0 * 1e6 * BR0["ppp0"] * CS0 * FE, a_eta["v"] * rho["v"] * SIG0 * 1e6 * BR0["ppee"] * CS0 * FE))
# units check of the side-band yield: raw limit 5008.5 x 64/63
P("units: Y_peak raw = %.1f (WP6 iteration limit 5008.5); a_P from raw rows / (SMASH x 63/64) = %.6f" % (Yp["v"] / YF, (Yp["v"] / YF) / (peakS * 63 / 64)))
P("L_sim factor: 6079.2/221.05 x 63/64 = %.4f (WP0/WP6: 27.0718)" % (L / 221.05 * 63 / 64))

# ---------------------------------------------------------------- B. N rows (normalisation), own propagation
def peak_br_err(B, E):
    c4 = sims["C4P"] * B["ppg"] / BR0["ppg"]; c5 = sims["C5P"] * B["ppee"] / BR0["ppee"]
    s4, s5 = c4 / (c4 + c5), c5 / (c4 + c5)
    return math.hypot(s4 * E["ppg"] / B["ppg"], s5 * E["ppee"] / B["ppee"]), s4
eb0, s40 = peak_br_err(BR0, BRE0)
eL = 0.069
eC_up = math.hypot(0.03 / FE, 0.0035 / CS0); eC_dn = math.hypot(0.04 / FE, 0.0035 / CS0)
nup = math.sqrt(eL ** 2 + eC_up ** 2 + eb0 ** 2); ndn = math.sqrt(eL ** 2 + eC_dn ** 2 + eb0 ** 2)
P("")
P("## B. Normalisation (N) of sigma, own propagation")
P("peak BR mix (PDG 2020, C4 share %.4f): %.2f %%  [WP7 2.3 %%];  C_sec(1-f_eta): +%.2f/-%.2f %%;  L 6.9 %%" % (s40, 100 * eb0, 100 * eC_up, 100 * eC_dn))
P("N(sigma) = +%.2f/-%.2f %%  [WP7 +7.915/-8.357 %%]" % (100 * nup, 100 * ndn))
ebump0 = math.hypot(BRE0["ppp0"] / BR0["ppp0"], BRE0["p0D"] / BR0["p0D"])
P("bump/C2 BRs (ppp0 (+) pi0 Dalitz): %.2f %% [WP7 3.22 %%]; DR N = %.2f %% [3.97 %%]" % (100 * ebump0, 100 * math.hypot(ebump0, eb0)))
st, sy, mc = syst7["sigma_eta"]["stat"] / syst7["sigma_eta"]["v"], None, syst7["sigma_eta"]["mc"] / syst7["sigma_eta"]["v"]
sy0 = 0.16008 + 0  # WP7/WP8 total syst (S12 = 2.5 % pion vertex)
tot0 = math.sqrt(st ** 2 + sy0 ** 2 + mc ** 2 + ndn ** 2)
P("total (stat %.2f, syst %.2f, MC %.2f, N-dn %.2f %%) = %.2f %% -> %.4f mb  [WP8 20.7 %%, 0.091 mb]" % (100 * st, 100 * sy0, 100 * mc, 100 * ndn, 100 * tot0, tot0 * sig))
P("live time: 1/0.932 - 1 = %.2f %% (spec +7.3 %%); with the L-weighted LF3/LF2 x 64/63 = 0.9334 (5675.7/6080.5): %.2f %%" % (100 * (1 / 0.932 - 1), 100 * (6080.5 / 5675.7 - 1)))

# ---------------------------------------------------------------- C. corrections (confirmed factor errors)
P("")
P("## C. Corrected factors (v8): PDG 2026 BRs, SMASH sigma of the simulated sample, C_sec of the peak")
f4, f5 = BR1["ppg"] / BR0["ppg"], BR1["ppee"] / BR0["ppee"]; f23 = BR1["ppp0"] / BR0["ppp0"]
peak1 = sims["C4P"] * f4 + sims["C5P"] * f5; bump1 = bumpS * f23
aP1 = Yp["v"] / peak1
sig1 = aP1 * SIG1 * CS1 * FE
DR1 = aP1 / (Yb["v"] / bump1)
BRG1 = BR_G * f23; BRA1 = BR_A * f23
P("class factors PDG2026/PDG2020: C2,C3 x %.5f; C4 x %.5f; C5 x %.5f" % (f23, f4, f5))
P("SMASH peak in P: %.2f -> %.2f rows (x %.5f); bump in B x %.5f" % (peakS, peak1, peak1 / peakS, f23))
P("a_P: %.5f -> %.5f ; sigma_SMASH(sample) %.4f mb (was 1.0665); C_sec(peak) 1.0116 (was 1.0110 mean)" % (aP, aP1, SIG1))
P("sigma(pp -> eta X) = %.5f x %.4f x 1.0116 x 0.97 = %.4f mb   (WP8 0.4374; change %+.2f %%)" % (aP1, SIG1, sig1, 100 * (sig1 / sig - 1)))
P("   of which: PDG BRs %+.2f %%, SMASH sample sigma %+.2f %%, C_sec(peak) %+.2f %%" % (100 * (peakS / peak1 - 1), 100 * (SIG1 / SIG0 - 1), 100 * (CS1 / CS0 - 1)))
P("sigma/sigma_SMASH(sample) = %.4f (unchanged by the SMASH factor)" % (sig1 / SIG1))
P("DR: %.4f -> %.4f ; DR_F: %.4f -> %.4f" % (DR, DR1, DRF, DRF * f23 / (peak1 / peakS)))
P("BR_G: %.4e -> %.4e ; BR_A: %.4e -> %.4e (scale with BR(eta->pi+pi-pi0) of the C2 normalisation; the 2.68e-4 factor cancels)" % (BR_G, BRG1, BR_A, BRA1))
P("PDG reference for the BR comparison: 2.76 +- 0.16 e-4 (PDG 2026; 2.68 +- 0.11 e-4 in PDG 2020-2024)")
eb1, s41 = peak_br_err(BR1, BRE1)
ebump1 = math.hypot(BRE1["ppp0"] / BR1["ppp0"], BRE1["p0D"] / BR1["p0D"])
nup1 = math.sqrt(eL ** 2 + eC_up ** 2 + eb1 ** 2); ndn1 = math.sqrt(eL ** 2 + eC_dn ** 2 + eb1 ** 2)
P("N with PDG 2026: peak BR mix %.2f %% (C4 share %.4f); N(sigma) +%.2f/-%.2f %%; bump/C2 BRs %.2f %%; DR N %.2f %%" % (
    100 * eb1, s41, 100 * nup1, 100 * ndn1, 100 * ebump1, 100 * math.hypot(ebump1, eb1)))

# ---------------------------------------------------------------- D. missing systematic: peak composition (conversion scale of C4)
P("")
P("## D. Missing systematic of sigma (side bands): composition of the peak (photon-conversion scale of C4)")
c4, c5 = sims["C4P"] * f4, sims["C5P"] * f5
up = (c4 + c5) / (0.8 * c4 + c5); dn = (c4 + c5) / (1.2 * c4 + c5)
s13 = (abs(up - 1) + abs(dn - 1)) / 2
P("a = a_P (C4 + C5)/(kappa C4 + C5); kappa = 1 +- 0.2 (the analysis' own prior, WP3/WP5 fall-back): %+.2f %% / %+.2f %% -> +-%.2f %%" % (100 * (up - 1), 100 * (dn - 1), 100 * s13))
# data indication from the per-class side-band yields (E1 < 10, E2 >= 10 MeV), PDG BR(ppee) fixed (rho = 1)
c4E2 = 45006.3 * 3.39316e-6 * L * f4; c5E2 = 285.822 * 2.43037e-3 * L * f5     # Table P1 (S0, V3) x L, PDG 2026
c4E1 = c4 - c4E2; c5E1 = c5 - c5E2
y1, y2 = YpE1["v"], YpE2["v"]; e1 = math.hypot(YpE1["stat"], YpE1["mc"]); e2 = math.hypot(YpE2["stat"], YpE2["mc"])
det = c4E1 * c5E2 - c4E2 * c5E1
A = (c4E1 * y2 - c4E2 * y1) / det          # a (direct scale, rho = 1)
U = (y1 * c5E2 - y2 * c5E1) / det          # a kappa
dA = math.hypot(c4E1 * e2, c4E2 * e1) / det
dU = math.hypot(c5E2 * e1, c5E1 * e2) / det
P("SMASH peak by class (rows, PDG 2026): C4 E1 %.0f E2 %.0f ; C5 E1 %.0f E2 %.0f ; data E1 %.0f +- %.0f, E2 %.0f +- %.0f (stat (+) MC)" % (
    c4E1, c4E2, c5E1, c5E2, y1, e1, y2, e2))
P("solve with rho = 1: a = %.3f +- %.3f, a kappa = %.3f +- %.3f -> kappa = %.2f ; sigma would be %.3f +- %.3f mb (vs a_P route %.3f)" % (
    A, dA, U, dU, U / A, A * SIG1 * CS1 * FE, dA * SIG1 * CS1 * FE, sig1))
P("(indication only: E1 is the class with the unmodelled component, k_E1 = 7.6; the F-main kappa is 0.83 +- 0.30 with the prior, 0.33 +- 0.52 free)")
fam = syst7["sigma_eta"]["fam"]
s12_old = fam[11]
fam_new = [x / syst7["sigma_eta"]["v"] for x in fam]
s12_new = max(s12_old / syst7["sigma_eta"]["v"], s13)
fam_new[11] = s12_new
sy1 = math.sqrt(sum(x * x for x in fam_new))
sy1q = math.sqrt(sum(x * x for x in fam_new[:11]) + (s12_old / syst7["sigma_eta"]["v"]) ** 2 + s13 ** 2)
tot1 = math.sqrt(st ** 2 + sy1 ** 2 + mc ** 2 + ndn1 ** 2)
P("syst(sigma) with the composition term in S12 (spec max rule): %.2f %% (was %.2f %%); as a separate family: %.2f %%" % (100 * sy1, 100 * sy0, 100 * sy1q))

# ---------------------------------------------------------------- E. corrected headline
P("")
P("## E. Corrected headline (v8)")
P("sigma(pp -> eta X) = %.3f +- %.3f (stat) +- %.3f (syst) +- %.3f (MC) +%.3f/-%.3f (norm.) mb; total %.3f mb (%.1f %%)" % (
    sig1, st * sig1, sy1 * sig1, mc * sig1, nup1 * sig1, ndn1 * sig1, tot1 * sig1, 100 * tot1))
P("   data/SMASH = %.3f +- %.3f;  sigma x BR(eta -> pi+pi-e+e-, PDG 2026) = %.1f nb;  sigma x BR(pi+pi-gamma) = %.0f nb" % (
    sig1 / SIG1, tot1 * sig1 / SIG1, sig1 * 1e6 * BR1["ppee"], sig1 * 1e6 * BR1["ppg"]))
P("   one-sided (not in the total, all would RAISE sigma): live time +7.1..+7.3 %; START-based track-sorter loss +1.3..+2.8 % (report_exp);")
P("   cor_ok rule +2.2 +- 1.0 % (data drop cor_ok == 0 rows, the simulation has none: checked, 0 of 230469 std / 7441490 lep rows);")
P("   event-vertex efficiency of vz > -500 (data keeps 0.565 of the e+e- rows in P, about 0.68 of N - CB; simulation >= 0.96):")
P("   NOT MEASURED for the eta signal, up to x 1.4 (N - CB rate) or x 1.7 (row rate); (pT, y) production model +19/-7 % (WP8).")
P("BR_G = %.2f e-4 (vs PDG 2026 2.76 +- 0.16); BR_A = %.2f e-4; DR = %.3f; DR_F = %.3f" % (BRG1 * 1e4, BRA1 * 1e4, DR1, DRF * f23 / (peak1 / peakS)))

# ---------------------------------------------------------------- F. literature comparison
P("")
P("## F. Literature (verified on line, 2026-11-14)")
lit_incl = [(1.14, 0.20, "HADES p+p 3.5 GeV inclusive, PYTHIA extrapolation"), (0.93, 0.14, "same, UrQMD extrapolation")]
for v, e, t in lit_incl:
    z = (v - sig1) / math.hypot(e, tot1 * sig1)
    P("  %-52s %.2f +- %.2f mb ; this/lit = %.2f ; difference %.1f sigma" % (t, v, e, sig1 / v, z))
P("  [G. Agakishiev et al. (HADES), Eur. Phys. J. A 48 (2012) 64, arXiv:1112.3607, Table 1; T = 3.5 GeV, sqrt(s) = 3.18 GeV]")
P("  preliminary, same beam time (eta -> pi+pi-pi0, pi0 -> gamma gamma, ECAL): 1.15 mb (+-1 sigma window) / 1.72 mb (+-2 sigma), no errors")
P("  [A. Strach, talk, Workshop at 1 GeV scale: From mesons to axions, Jagiellonian University, Krakow, 20 Sep 2024, slide 13] -> ratio %.2f / %.2f" % (sig1 / 1.15, sig1 / 1.72))
with open(os.path.join(TAB, "v8_verification.txt"), "w") as f:
    f.write("\n".join(OUT) + "\n")

# ---------------------------------------------------------------- TeX snippets
tex = []
tex.append("% v8_tableP2_corr.tex -- plots/all_report/v8_verify.py (verifier v8, 2026-11-14): corrected sigma(pp -> eta X) and literature")
tex.append("\\begin{table}[H]\n\\centering\n\\caption{Inclusive $\\eta$ production cross section after the independent verification (v8). "
           "Peak side bands, S0, $\\mathrm{isBest}\\geq0$, $L=6079.2$\\,nb$^{-1}$; PDG 2026 branching ratios; "
           "$\\sigma_{\\mathrm{SMASH}}$ of the simulated sample $1.069$\\,mb; $C_{\\mathrm{sec}}(1-f_\\eta)=1.0116\\times0.97$. "
           "Systematic error includes the composition of the peak (conversion scale of C4, $\\kappa=1\\pm0.2$). One-sided items are listed apart; "
           "all of them would raise $\\sigma$.}\n\\label{tab:v8sigma}\n\\small\n\\begin{tabular}{@{}p{0.37\\textwidth} p{0.59\\textwidth}@{}}\n\\toprule")
tex.append("$\\sigma(pp\\to\\eta X)$ & $%.3f\\pm%.3f\\,(\\mathrm{stat})\\pm%.3f\\,(\\mathrm{syst})\\pm%.3f\\,(\\mathrm{MC})^{+%.3f}_{-%.3f}\\,(\\mathrm{norm.})$\\,mb \\\\" % (
    sig1, st * sig1, sy1 * sig1, mc * sig1, nup1 * sig1, ndn1 * sig1))
tex.append("$\\sigma/\\sigma_{\\mathrm{SMASH}}$ & $%.3f\\pm%.3f$ \\\\" % (sig1 / SIG1, tot1 * sig1 / SIG1))
tex.append("\\midrule\nlive time (PT3 live fraction) & $+7.1$ to $+7.3$\\,\\% \\\\\nSTART-based track-sorter loss & $+1.3$ to $+2.8$\\,\\% \\\\")
tex.append("\\texttt{cor\\_ok} rule (absent in the simulation) & $+2.2\\pm1.0$\\,\\% \\\\")
tex.append("event-vertex efficiency of $v_z>-500$\\,mm & not measured for the $\\eta$; data keep 0.565 of the rows in P, simulation $\\geq0.96$ \\\\")
tex.append("$(p_T,y)$ production model & $+19/-7$\\,\\% (WP8) \\\\")
tex.append("\\midrule\n\\multicolumn{2}{@{}l}{\\emph{Literature}} \\\\")
tex.append("HADES, $p+p$ 3.5\\,GeV, inclusive (EPJA 48 (2012) 64) & $1.14\\pm0.20$\\,mb (PYTHIA), $0.93\\pm0.14$\\,mb (UrQMD) \\\\")
tex.append("HADES, same beam time, preliminary talk (Krak\\'ow 2024) & 1.15\\,mb ($\\pm1\\sigma$ window), 1.72\\,mb ($\\pm2\\sigma$); no errors, unpublished \\\\")
tex.append("HADES, $pp\\to pp\\eta$ 3.5\\,GeV (exclusive) & $0.211\\pm0.047$\\,mb (prelim., IJMPA 26 (2011) 694); $0.137\\pm0.010$\\,mb (thesis) \\\\")
tex.append("$pp\\to pp\\eta$ models, $\\sqrt{s}=3.46$\\,GeV (arXiv:2501.14952) & 0.081 / 0.116 / 0.137\\,mb \\\\")
tex.append("\\bottomrule\n\\end{tabular}\n\\end{table}")
with open(os.path.join(TAB, "v8_tableP2_corr.tex"), "w") as f:
    f.write("\n".join(tex) + "\n")
print("wrote", os.path.join(TAB, "v8_tableP2_corr.tex"))
