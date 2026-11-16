#!/usr/bin/env python3
# ============================================================================
# fix_numbers.py  --  review fix (prefix fix): the numbers changed or added by the review of report_all, recomputed from
# the published inputs of the packages (no refit). Every input is quoted with its source; output tables_all/fix_numbers.txt.
#   1. composition of the eta peak from the data (reviewer 1, finding 2): per-M(ee)-class peak side bands (WP6/v8) and the
#      X17 2D fit (x17d N5, N4); direct-decay and conversion routes; composition systematic of sigma(pp -> eta X);
#   2. sigma(pp -> eta X): new systematic and total, data/SMASH, alternative values (routes), known-sign items;
#   3. consistent PDG-2026 rescaling of the cross-check values (F-main, edge, P_main r = 1 / 0.89, tag) and SMASH visible;
#   4. a_B(E2), a_B(all), DR in the kappa conventions; sigma x BR(eta -> pi+pi-pi0);
#   5. displaced-window closure: uncertainty of the rms (non-overlapping subsets, n_eff); E2 sigma_tr;
#   6. S8 widened to the observed range of g (linear extrapolation of the WP7 "g free" variant);
#   7. alternative totals (WP7 rules / without R1 / literal Barlow) with the new composition term;
#   8. fold z values and calibrated per-class edge-step errors.
# Run: python3 -I plots/all_report/fix_numbers.py   (from /home/przygoda/HADES/FAT/FAT; < 1 s)
# ============================================================================
import math
import numpy as np

OUT = "/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all/fix_numbers.txt"
F = open(OUT, "w")
def P(*a):
    s = " ".join(str(x) for x in a); print(s); F.write(s + "\n")
def q(*v): return math.sqrt(sum(x * x for x in v))

# ---------------------------------------------------------------- inputs (sources in brackets)
SIG_SMASH = 1.0694          # mb, simulated sample (v8)
CSEC_P, CSEC_B, CSEC_M = 1.0116, 1.0104, 1.0110   # WP2 (peak, bump, mean)
ONE_FETA = 0.97             # Amendment 1
K_P = CSEC_P * ONE_FETA
PDG26 = dict(C2=1.00436, C3=1.00436, C4=1.01422, C5=1.02985)   # v8 rescaling factors PDG 2020 -> 2026
YP, YP_STAT, YP_MC = 5087.97, 431.6, 245.0                     # peak side bands (WP6)
SM_P20 = 12164.94                                              # SMASH C4 + C5 in P, PDG 2020, rows at L (v8)
C4P_nb, C5P_nb = 0.954863, 1.04621                             # SMASH in P [nb], PDG 2020 (v8)
s4_20 = C4P_nb / (C4P_nb + C5P_nb)
C4P26, C5P26 = C4P_nb * PDG26["C4"], C5P_nb * PDG26["C5"]
s4 = C4P26 / (C4P26 + C5P26)
aP20 = YP / SM_P20
aP = aP20 * (C4P_nb + C5P_nb) / (C4P26 + C5P26)               # PDG 2026
P("# fix_numbers.txt -- review fixes of report_all (fix_numbers.py); PDG 2026 unless stated")
P("a_P (PDG 2020) = %.5f ; a_P (PDG 2026) = %.5f ; C4 share of the simulated peak in P: %.4f (2020), %.4f (2026)" % (aP20, aP, s4_20, s4))
sig_adopt = aP * SIG_SMASH * K_P
P("sigma_eta adopted (v8) = a_P x 1.0694 x 1.0116 x 0.97 = %.4f mb" % sig_adopt)

# ---------------------------------------------------------------- 1. composition from the data
P("\n## 1. composition of the eta peak (data)")
# (a) per-M(ee)-class peak side bands (WP6), SMASH shapes per class (v8; rows, these are PDG 2026 numbers: 4946 + 942 = 5888 = 5804.8 x 1.0142)
YE1, YE1s, YE1m = 2380., 291., 164.
YE2, YE2s, YE2m = 2780., 326., 179.
c4 = (4946., 942.); c5 = (2201., 4349.)
# E2 displaced-window closure term (WP6: rms 1.572 in E2 alone): sigma_tr(E2) = sigma_stat(E2) sqrt(rms^2 - 1)
str_E2 = YE2s * math.sqrt(1.572 ** 2 - 1)
P("per-class side-band yields: E1 %.0f +- %.0f (stat) +- %.0f (MC); E2 %.0f +- %.0f +- %.0f; E2 closure term sigma_tr(E2) = %.0f rows (rms 1.572)" % (YE1, YE1s, YE1m, YE2, YE2s, YE2m, str_E2))
det = c4[0] * c5[1] - c5[0] * c4[1]
def solve(y1, y2):   # y = K C4 + A C5 per class -> (A = a, K = a kappa)
    A = (c4[0] * y2 - c4[1] * y1) / det
    Kc = (c5[1] * y1 - c5[0] * y2) / det
    return A, Kc
A, Kc = solve(YE1, YE2)
def errAK(s1, s2):
    sA = math.sqrt((c4[0] * s2) ** 2 + (c4[1] * s1) ** 2) / det
    sK = math.sqrt((c5[1] * s1) ** 2 + (c5[0] * s2) ** 2) / det
    cov = (-(c4[1] * c5[1]) * s1 ** 2 - (c4[0] * c5[0]) * s2 ** 2) / det ** 2
    return sA, sK, cov
s1, s2 = q(YE1s, YE1m), q(YE2s, YE2m)
sA0, sK0, _ = errAK(s1, s2)
s2t = q(YE2s, YE2m, str_E2)
sA, sK, cAK = errAK(s1, s2t)
kap_sb = Kc / A
skap_sb = kap_sb * math.sqrt((sK / Kc) ** 2 + (sA / A) ** 2 - 2 * cAK / (A * Kc))
P("side bands per class (rho = 1): direct route a = %.3f +- %.3f (stat+MC; v8 0.59 +- 0.10) +- -> %.3f with the E2 closure term; conversion route a kappa = %.3f +- %.3f; kappa = %.2f +- %.2f"
  % (A, sA0, sA, Kc, sK, kap_sb, skap_sb))
P("per-class a_P: E1 %.3f, E2 %.3f (SMASH E2 share %.3f, data %.3f)" % (YE1 / (c4[0] + c5[0]), YE2 / (c4[1] + c5[1]), (c4[1] + c5[1]) / (sum(c4) + sum(c5)), YE2 / (YE1 + YE2)))
# (b) X17 2D fit (x17d, r = 0, all cells): N5 = 3659 +- 503, N4 = 1648 +- 383; references at a_eta = 1 (PDG 2020): 7458.8, 6656.3
N5, sN5, N4, sN4 = 3659., 503., 1648., 383.
N5r, N4r = 7458.8 * PDG26["C5"], 6656.3 * PDG26["C4"]
a5, sa5 = N5 / N5r, sN5 / N5r
a4, sa4 = N4 / N4r, sN4 / N4r
P("X17 2D fit: N5/N5ref = %.3f +- %.3f, N4/N4ref = %.3f +- %.3f (PDG 2026; PDG 2020: %.3f, %.3f); kappa = %.2f" % (a5, sa5, a4, sa4, N5 / 7458.8, N4 / 6656.3, a4 / a5))
P("F-main E1+E2: kappa = 0.504 +- 0.162 (WP5); F-main E2: a_eta rho_d = 0.3823 x 1.3152 = %.3f (PDG 2020) -> %.3f (2026)" % (0.38228 * 1.3152, 0.38228 * 1.3152 / PDG26["C5"]))
w1, w2 = 1 / sA ** 2, 1 / sa5 ** 2
a_dir = (A * w1 + a5 * w2) / (w1 + w2); sa_dir = 1 / math.sqrt(w1 + w2)
w1c, w2c = 1 / sK ** 2, 1 / sa4 ** 2
a_conv = (Kc * w1c + a4 * w2c) / (w1c + w2c); sa_conv = 1 / math.sqrt(w1c + w2c)
P("direct-decay route (weighted mean of side bands and X17 fit): a = %.3f +- %.3f ; conversion route: a = %.3f +- %.3f" % (a_dir, sa_dir, a_conv, sa_conv))
comp = 0.5 * (a_dir - a_conv) / aP
P("COMPOSITION SYSTEMATIC of sigma_eta = half the difference of the two routes / a_P = (%.3f - %.3f)/2/%.4f = +-%.1f %%" % (a_dir, a_conv, aP, 100 * comp))
P("  for comparison: X17 fit alone +-%.1f %%; side bands alone +-%.1f %%; v8 prior kappa = 1 +- 0.2: +10.5/-8.7 %% (+-9.6 %%)"
  % (50 * (a5 - a4) / aP, 50 * (A - Kc) / aP))
for kap in (0.37, 0.46, 0.50, 0.51, 0.8, 1.2):
    up = 1 / (kap * s4 + 1 - s4); dn = 1 / (s4 + (1 - s4) / kap)
    P("  kappa_data = %.2f: a = a_P x %.3f if C5 is right (direct reading), a_P x %.3f if C4 is right (conversion reading)" % (kap, up, dn))

# ---------------------------------------------------------------- 2. sigma_eta
P("\n## 2. sigma(pp -> eta X)")
S_WP7 = 16.01          # WP7/WP8 total syst of sigma (per cent) incl. S12 = 2.5 and S8 = 0.34
S8_new = 1.0           # S8 widened to g = 0.47-2.6 (section 6)
syst = math.sqrt(S_WP7 ** 2 - 2.5 ** 2 - 0.34 ** 2 + S8_new ** 2 + (100 * comp) ** 2)
stat = 8.5; mc = 5.9    # WP7 Table S1, sigma_eta column (MC includes the SMASH denominator)
Nup, Ndn = 8.20, 8.64  # v8 (PDG 2026)
tot_up = math.sqrt(stat ** 2 + syst ** 2 + mc ** 2 + Nup ** 2); tot_dn = math.sqrt(stat ** 2 + syst ** 2 + mc ** 2 + Ndn ** 2)
tot = 0.5 * (tot_up + tot_dn)
P("syst: sqrt(16.01^2 - 2.5^2 - 0.34^2 + %.1f^2 (S8) + %.1f^2 (composition, S12 family maximum)) = %.1f %%" % (S8_new, 100 * comp, syst))
P("sigma = %.3f +- %.3f (stat) +- %.3f (syst) +- %.3f (MC) +%.3f/-%.3f (norm.) mb; total %.1f %% (+%.1f/-%.1f)"
  % (sig_adopt, sig_adopt * stat / 100, sig_adopt * syst / 100, sig_adopt * mc / 100, sig_adopt * Nup / 100, sig_adopt * Ndn / 100, tot, tot_up, tot_dn))
rS = sig_adopt / SIG_SMASH
P("sigma/sigma_SMASH = %.3f +- %.3f (C_sec(1 - f_eta) included)" % (rS, rS * tot / 100))
P("syst with the spec's prior kappa = 1 +- 0.2 (v8): 18.5 %; without any composition term: 16.0 %")
for nm, a, sa in (("direct route (side bands per class, rho = 1)", A, sA), ("direct route (X17 fit N5)", a5, sa5), ("direct route, weighted", a_dir, sa_dir),
                  ("conversion route (side bands)", Kc, sK), ("conversion route (X17 fit N4)", a4, sa4), ("conversion route, weighted", a_conv, sa_conv)):
    P("  sigma via the %-44s = %.3f +- %.3f (stat) mb" % (nm, a * SIG_SMASH * K_P, sa * SIG_SMASH * K_P))
P("known-sign items (central value of sigma after the item; not applied):")
items = [("r = 0.89 (data peak width 5.44 vs 6.09 MeV)", 1 - 0.108), ("f_eta = 12 % (upper end of the nucleon-count expectation)", 0.88 / 0.97),
         ("y-only reweighting", 1 - 0.080), ("QED gamma* angle of C2", 1 - 0.001), ("live time", 1.071), ("live time (upper)", 1.073),
         ("START-based sorter", 1.013), ("sorter (upper)", 1.028), ("cor_ok rule", 1.022), ("pT-only reweighting", 1.210),
         ("event-vertex efficiency at the N - CB rate", 1.4), ("event-vertex efficiency at the raw-row rate", 1.7),
         ("composition, direct route (weighted)", a_dir / aP), ("composition, conversion route (weighted)", a_conv / aP)]
for nm, f in items:
    P("  %-60s x %.3f -> %.3f mb" % (nm, f, sig_adopt * f))
up_lo = 1.4 * 1.071 * 1.013 * 1.022; up_hi = 1.7 * 1.073 * 1.028 * 1.022
P("  vertex x live time x sorter x cor_ok: x %.2f-%.2f -> %.2f-%.2f mb; with the direct route in addition: %.2f-%.2f mb; with pT also: up to %.2f mb"
  % (up_lo, up_hi, sig_adopt * up_lo, sig_adopt * up_hi, sig_adopt * up_lo * a_dir / aP, sig_adopt * up_hi * A / aP, sig_adopt * up_hi * A / aP * 1.21))
for nm, v, e in (("HADES 3.5 GeV PYTHIA", 1.14, 0.20), ("HADES 3.5 GeV UrQMD", 0.93, 0.14)):
    P("  %s %.2f +- %.2f mb: ratio %.2f; naive distance with the new total %.1f sigma (not a significance: one-sided items missing)"
      % (nm, v, e, sig_adopt / v, (v - sig_adopt) / math.sqrt(e ** 2 + (sig_adopt * tot / 100) ** 2)))

# ---------------------------------------------------------------- 3. cross-checks, PDG 2026, 1.0694 mb
P("\n## 3. cross-check values of sigma, all with PDG 2026, sigma_SMASH = 1.0694 mb and the C_sec of their window")
xc = [("F-main a_eta (bump-dominated, C2/C3 rescaled, mean C_sec)", 0.38228 / PDG26["C2"], CSEC_M),
      ("f_eta(edge) = 1141.9/2915.7 (C2+C3)", (1141.9 / 2915.7) / PDG26["C2"], CSEC_B),
      ("f_eta(tag) = 0.454 (C2+C3, emulated tag)", 0.454 / PDG26["C2"], CSEC_B),
      ("P_main S_data/S_sim, r = 1 (4388/10881)", (4388. / 10881.) / ((C4P26 + C5P26) / (C4P_nb + C5P_nb)), CSEC_P),
      ("P_main S_data/S_sim, r = 0.89 (4388/12098)", (4388. / 12098.) / ((C4P26 + C5P26) / (C4P_nb + C5P_nb)), CSEC_P),
      ("peak side bands a_P (adopted)", aP, CSEC_P)]
for nm, a, cs in xc:
    P("  %-60s scale %.4f -> sigma %.3f mb (WP8 value with PDG 2020, 1.0665: %.3f)" % (nm, a, a * SIG_SMASH * cs * ONE_FETA,
      (a * (PDG26["C2"] if "bump" in nm or "edge" in nm or "tag" in nm else (C4P26 + C5P26) / (C4P_nb + C5P_nb) if "P_main" in nm else 1)) * 1.0665 * CSEC_M * ONE_FETA))
P("SMASH visible: C4 + C5 in P %.3f nb (PDG 2020 2.001), C2 + C3 in B %.3f nb (PDG 2020 9.008)" % (C4P26 + C5P26, 9.00759 * PDG26["C2"]))

# ---------------------------------------------------------------- 4. a_B and DR conventions
P("\n## 4. bump scales and DR")
YB_E2, YB_all = 7822.17, 16389.3
SB_E2, SB_all = 25097.0, 54758.9        # SMASH C2 + C3 (kappa = 1) in B, rows (v8, PDG 2020)
aBE2_20, aBall_20 = YB_E2 / SB_E2, YB_all / SB_all
aBE2, aBall = aBE2_20 / PDG26["C2"], aBall_20 / PDG26["C2"]
P("a_B(E2) = Y_bump(E2)/(C2+C3)_E2 = %.4f (2020) %.4f (2026); a_B(all) = Y_bump(all)/(C2+C3)_all = %.4f (2020) %.4f (2026)" % (aBE2_20, aBE2, aBall_20, aBall))
P("  (Y_bump(all) = Y_bump(E2)/f_E2 with f_E2(kappa = 0.832) = 0.4773; the SMASH denominators use kappa = 1)")
DR20, DR26 = aP20 / aBall_20, aP / aBall
P("DR = a_P / a_B(all) = %.4f (2020) %.4f (2026)  [the adopted definition]" % (DR20, DR26))
fE2_k1 = SB_E2 / SB_all
P("kappa = 1 throughout (f_E2 = %.4f): a_B(all) = a_B(E2): DR = %.3f (2020) %.3f (2026)" % (fE2_k1, aP20 / aBE2_20, aP / aBE2))
C2E2, C3E2 = 2870 * 0.001150, 244442 * 0.00000339       # nb (Table P1, PDG 2020)
rk = (C2E2 + 0.832 * C3E2) / (C2E2 + C3E2)
aPk20 = aP20 / (0.832 * s4_20 + 1 - s4_20); aPk26 = aP / (0.832 * s4 + 1 - s4)
P("kappa = 0.832 throughout (numerator and denominators): a_B = %.4f, a_P = %.4f: DR = %.3f (2020) %.3f (2026)" % (aBE2_20 / rk, aPk20, aPk20 / (aBE2_20 / rk), aPk26 / (aBE2 / rk)))
for nm, a in (("a_B(all)", aBall_20), ("a_B(E2)", aBE2_20)):
    P("sigma x BR(eta -> pi+pi-pi0) from %s = a x 1.0694 mb x 0.2292 x C_sec(bump) x 0.97 = %.1f ub (PDG-edition independent)" % (nm, a * SIG_SMASH * 0.2292 * CSEC_B * ONE_FETA * 1e3))

# ---------------------------------------------------------------- 5. displaced windows
P("\n## 5. displaced-window closure of the peak side bands")
pulls = [0.97, 0.12, -0.49, -0.38, 0.84, -0.78, 1.94, 1.10, -1.17, 2.29, 0.64, 0.08, 1.53]
pu = np.array(pulls); rms = math.sqrt((pu ** 2).mean())
P("13 windows (10 MeV steps, 37.8 MeV wide): rms %.3f, mean %+.3f" % (rms, pu.mean()))
sst = 431.6
for off in range(4):
    sub = pu[off::4]; r_ = math.sqrt((sub ** 2).mean())
    P("  non-overlapping subset %d (%s): rms %.3f -> sigma_tr = %.0f rows" % (off, ", ".join("%+.2f" % x for x in sub), r_, sst * math.sqrt(max(0, r_ ** 2 - 1))))
neff = 13 * 10 / 37.8
ms = (pu ** 2).mean(); sms = ms * math.sqrt(2 / neff)
lo, hi = max(0., ms - sms), ms + sms
P("n_eff = %.1f: rms^2 = %.2f +- %.2f -> rms %.2f (+%.2f/-%.2f); sigma_tr = %.0f +%.0f/-%.0f rows; mean pull %+.2f +- %.2f"
  % (neff, ms, sms, rms, math.sqrt(hi) - rms, rms - math.sqrt(lo), sst * math.sqrt(ms - 1), sst * math.sqrt(hi - 1) - sst * math.sqrt(ms - 1), sst * math.sqrt(ms - 1), pu.mean(), rms / math.sqrt(neff)))
s10_hi = 100 * sst * math.sqrt(hi - 1) / YP
P("  with the upper value sigma_tr = %.0f rows (%.1f %%): Y_peak syst %.1f %% instead of 11.4 %%" % (sst * math.sqrt(hi - 1), s10_hi, math.sqrt(11.44 ** 2 - 4.7 ** 2 + s10_hi ** 2)))
P("E2 class: rms 1.572 -> sigma_tr(E2) = 326 x sqrt(1.572^2 - 1) = %.0f rows (%.0f %% of 2780)" % (str_E2, 100 * str_E2 / YE2))

# ---------------------------------------------------------------- 6. S8 widened
P("\n## 6. family S8 widened to the observed range of g (0.47-2.6), linear in g from the WP7 'g free' variant (g = 1.717 -> 2.56)")
dg = 2.56 - 1.71677
for nm, val, dfree, wp7 in (("Y_bump", 16389., 851.6, 4.55), ("Y_peak (and sigma_eta)", 5088., -35.31, 0.34), ("a_eta", 0.38228, -0.02347, 5.72), ("Y4", 1847., 92.28, 3.42)):
    d047 = dfree / dg * (0.47 - 1.71677)
    P("  %-24s g free %+.2f %%; extrapolated to g = 0.47: %+.1f %%; S8 WP7 %.2f %% -> widened %.1f %%" % (nm, 100 * dfree / val, 100 * d047 / val, wp7, max(wp7, abs(100 * d047 / val))))

# ---------------------------------------------------------------- 7. alternative totals
P("\n## 7. totals of the systematic uncertainty (per cent)")
P("Y_peak: WP7 rules (MC(Delta) subtraction + R1) 11.4; without R1 15.2; literal Barlow with stat only 14.7 (WP7 syst.txt)")
for nm, base in (("WP7 rules", 16.01), ("without R1", 18.6), ("literal Barlow (stat only)", 17.7)):
    P("sigma_eta, %-28s: without composition %.1f ; with the data composition term %.1f" % (nm, base, math.sqrt(base ** 2 - 2.5 ** 2 + (100 * comp) ** 2)))

# ---------------------------------------------------------------- 8. small items
P("\n## 8. small items")
for nm, v0, e0, v1, e1 in (("a_eta", 0.356, 0.0195, 0.3946, 0.0204),):
    P("fold z (printed errors, independent halves): %s %.2f" % (nm, (v0 - v1) / q(e0, e1)))
P("folds (WP7 wp7_tests.txt): a_eta z -1.37, Y_bump(F-main) -1.22, Y_peak(F-main) -1.55, DR_F -0.91, Y_peak (side bands) -1.43")
P("edge step per class, statistical error x 1.2 (Part II calibration): E1 0.024 -> %.3f, E2 0.036 -> %.3f" % (0.024 * 1.2, 0.036 * 1.2))
P("X17: r_UL = 3.06e-5/2.68e-4 = %.4f ; with PDG 2026 (2.76e-4): %.2e ; with the C4-replica offset (6 %%): %.2e (PDG 2026 %.2e)" % (3.0636e-5 / 2.68e-4, 3.0636e-5 * 2.76 / 2.68, 3.25e-5, 3.25e-5 * 2.76 / 2.68))
F.close()
