#!/usr/bin/env python3
# ============================================================================================================
# fix2_numbers.py -- final review of the four reports (report_exp, report_exp_short, report_all, report_all_short),
# 2026-11-14. Recomputes the numbers that the two final reviewers found wrong or incomplete. Inputs are published
# numbers only (no data are read): notes_all/ve_verify_note.md, ve_data_note.md, ve_sim_cause_note.md,
# tables_all/vtx_ana*.txt, vtx_numbers.txt, fix_numbers.txt, report_exp.tex (period normalisation +-5 %).
#   1. eps_vtx(eta) with the statistical error of the edge step E x 1.2 (Part II calibration), and the combinations
#      S(smooth like-sign CB) + T + E, S(nominal CB) + T + E, S(smooth) + E, T + E (bump);
#   2. the efficiency needed for 0.93 / 1.14 mb and its distance from every measure;
#   3. the normalisation error of sigma(pp -> eta X) with the +-5 % period normalisation of Part II added in quadrature,
#      the new totals, sigma/sigma_SMASH, and the normalisation error of sigma_vis(peak) = Y_peak/L;
#   4. the direct-decay route with the one-sided items, separately for its weighted value and its per-class side-band value;
#   5. the like-sign / unlike-sign ratio of the vertex loss; the crossing-point distances from the LH2 cell;
#   6. the close-pair reading (D32): inputs with their status.
# Output: docs/latex/tables_all/fix2_numbers.txt.   Run: python3 -I plots/all_report/fix2_numbers.py
# ============================================================================================================
import math

OUT = "/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all/fix2_numbers.txt"
L = []
def P(s=""):
    print(s); L.append(s)

def comb(ms):
    w = [1 / e ** 2 for _, e in ms]
    m = sum(v * wi for (v, _), wi in zip(ms, w)) / sum(w)
    e = 1 / math.sqrt(sum(w))
    chi2 = sum((v - m) ** 2 * wi for (v, _), wi in zip(ms, w))
    return m, e, chi2

def eps(a, ea, b, eb):   # a/(a+b), independent errors
    s = a + b
    return a / s, math.sqrt(b * b * ea * ea + a * a * eb * eb) / (s * s)

P("# fix2_numbers.txt -- plots/all_report/fix2_numbers.py (final review, 2026-11-14); inputs: verifier ve_verify, ve_data, ve_sim_cause,")
P("# tables_all/vtx_*.txt, vtx_numbers.txt, fix_numbers.txt, Part II (period normalisation +-5 %)")
P()
P("## 1. eps_vtx(eta) = X_V/(X_V + X_Nv) with the edge-step statistical error x 1.2 (Part II calibration)")
CAL = 1.2
eS_nom = eps(4657, 397, 902, 411)            # stat_g (event-group) errors, verifier
eS_smo = eps(4670, 370, 263, 324)            # smooth like-sign CB, mean of 12 variants, plain fit errors
eS_N = eps(4646, 370, 165, 324)              # unlike-sign rows only, plain fit errors
eT = eps(3297, 110, -4, 77)                  # ECAL-tagged eta
eE0 = eps(1241, 76, 21, 80)                  # edge step, fit errors (as printed before)
eE = eps(1241, 76 * CAL, 21, 80 * CAL)       # edge step, fit errors x 1.2
E_mod = 0.028                                # rms over 18 model variations applied to both classes
eE_tot = math.hypot(eE[1], E_mod)
P("E(V) = 1241 +- %.0f, E(Nv) = 21 +- %.0f (fit errors 76 / 80 x %.1f)" % (76 * CAL, 80 * CAL, CAL))
P("S nominal CB      : %.3f +- %.3f (stat_g); model 0.060 -> total %.3f" % (eS_nom[0], eS_nom[1], math.hypot(eS_nom[1], 0.060)))
P("S smooth CB       : %.3f +- %.3f" % eS_smo)
P("S N only (no CB)  : %.3f +- %.3f" % eS_N)
P("T (ECAL-tagged)   : %.3f +- %.3f" % eT)
P("E (edge step)     : %.3f +- %.3f (stat, x1.2) +- %.3f (model) = %.3f total   [before: +- %.3f stat, %.3f total]"
  % (eE[0], eE[1], E_mod, eE_tot, eE0[1], math.hypot(eE0[1], E_mod)))
m, e, c = comb([eS_smo, eT, (eE[0], eE_tot)])
P("combined S(smooth)+T+E : %.4f +- %.4f, chi2 %.2f/2; one-sided 95 %% CL lower bound %.4f   [before, E without x1.2: 0.9934 +- 0.0209]" % (m, e, c, m - 1.645 * e))
m2, e2, c2 = comb([(eS_nom[0], math.hypot(eS_nom[1], 0.060)), eT, (eE[0], eE_tot)])
P("combined S(nominal CB, stat_g (+) model)+T+E: %.4f +- %.4f, chi2 %.2f/2" % (m2, e2, c2))
m3, e3, c3 = comb([eS_smo, (eE[0], eE_tot)])
P("combined S(smooth)+E (no ECAL): %.4f +- %.4f, chi2 %.2f/1" % (m3, e3, c3))
m4, e4, c4 = comb([eT, (eE[0], eE_tot)])
P("bump scale, T + E: %.4f +- %.4f (chi2 %.2f/1); simulation C2+C3 0.9986 -> C_vtx(bump) = %.3f +- %.3f" % (m4, e4, c4, 0.9986 / m4, 0.9986 / m4 * e4 / m4))
P("adopted (rounded, unchanged): eps_vtx = 0.994 +- 0.021; C_vtx = 0.99970/0.994 = %.4f +- %.4f; 95 %% CL: C_vtx < %.4f"
  % (0.99970 / 0.994, 0.99970 / 0.994 * 0.021 / 0.994, 0.99970 / (0.994 - 1.645 * 0.021)))
P("V contains 'pre' (vz > -160): S(pre) 4671 +- 367 in the verifier's code; the extra rows -500 < vz <= -160 hold S = -19 +- 133;")
P("   S(V) = 4657 against S(pre) + S(slice) = %.0f; T(V) = 3297 (verifier) against T(pre) = 3314 +- 105 (Part II code, ve_data V: 3301)" % (4671 - 19))
P()
P("## 2. efficiency needed for 0.93 / 1.14 mb (sigma scales as eps_sim/eps_data) and its distance from every measure")
s0 = 0.429
for target in (0.93, 1.14):
    need = round(s0 * 0.99970 / target, 2)
    P("sigma = %.2f mb needs eps = %.2f: T %.1f sigma, E %.1f sigma (stat x1.2 (+) model; %.1f with stat only), S(smooth) %.1f sigma, S(nominal) %.1f sigma, combined %.0f sigma"
      % (target, need, (eT[0] - need) / eT[1], (eE[0] - need) / eE_tot, (eE[0] - need) / eE[1], (eS_smo[0] - need) / eS_smo[1],
         (eS_nom[0] - need) / math.hypot(eS_nom[1], 0.060), (m - need) / e))
P()
P("## 3. normalisation of sigma(pp -> eta X) with the period normalisation of Part II (+-5 %, steps of the per-day rates)")
stat, mc = 8.5, 5.9
comp = (0.497 - 0.237) / 2 / 0.40909 * 100
syst = math.sqrt(16.01 ** 2 - 2.5 ** 2 - 0.34 ** 2 + 1.0 ** 2 + comp ** 2)
Nup0, Ndn0 = 8.20, 8.64                       # L 6.9 (+) C_sec(1 - f_eta) +3.1/-4.1 (+) PDG peak BR mix 3.15 (v8)
PER = 5.0
Nup, Ndn = math.hypot(Nup0, PER), math.hypot(Ndn0, PER)
tu0 = math.sqrt(stat ** 2 + syst ** 2 + mc ** 2 + Nup0 ** 2); td0 = math.sqrt(stat ** 2 + syst ** 2 + mc ** 2 + Ndn0 ** 2)
tu = math.sqrt(stat ** 2 + syst ** 2 + mc ** 2 + Nup ** 2); td = math.sqrt(stat ** 2 + syst ** 2 + mc ** 2 + Ndn ** 2)
P("syst %.2f %% (composition %.2f %%); N before +%.2f/-%.2f %% -> with the period term +%.2f/-%.2f %%" % (syst, comp, Nup0, Ndn0, Nup, Ndn))
P("sigma = 0.429 +- %.3f (stat) +- %.3f (syst) +- %.3f (MC) +%.3f/-%.3f (norm.) mb" % (s0 * stat / 100, s0 * syst / 100, s0 * mc / 100, s0 * Nup / 100, s0 * Ndn / 100))
P("total: before %.2f %% (+%.2f/-%.2f); now %.2f %% (+%.2f/-%.2f)" % (0.5 * (tu0 + td0), tu0, td0, 0.5 * (tu + td), tu, td))
rS = s0 / 1.0694
P("sigma/sigma_SMASH = %.3f +- %.3f (before +- %.3f)" % (rS, rS * 0.5 * (tu + td) / 100, rS * 0.5 * (tu0 + td0) / 100))
Lr = 6.9
P("sigma_vis(peak) = Y_peak/L = 0.837 visible rows per nb^-1: normalisation L 6.9 %% -> +-%.3f; with the period term sqrt(6.9^2 + 5^2) = %.2f %% -> +-%.3f"
  % (0.837 * Lr / 100, math.hypot(Lr, PER), 0.837 * math.hypot(Lr, PER) / 100))
P("sigma x BR (eta -> pi+pi-e+e-) 118.5 nb and (eta -> pi+pi-gamma) 18.4 ub keep +-38 %% (now %.1f %%)" % (0.5 * (tu + td)))
P("X17 sigma x BR x BR: refolded with s_rel = %.3f in fix2_x17_fold.py (tables_all/fix2_x17_fold.txt)" % (0.5 * (tu + td) / 100))
P()
P("## 4. sigma with the one-sided items (live time x sorter x cor_ok), direct-decay route split by its two values")
lt = (1.071, 1.073); so = (1.013, 1.028); co = 1.022
plo = lt[0] * so[0] * co; phi = lt[1] * so[1] * co; v95 = 0.99970 / (0.994 - 1.645 * 0.021)
P("product %.4f - %.4f; C_vtx 95 %% bound x%.4f" % (plo, phi, v95))
for nm, v in (("adopted", 0.429), ("direct route, weighted", 0.522), ("direct route, per-class side bands", 0.621)):
    P("  %-36s %.3f mb -> one-sided items %.3f - %.3f; + pT (x1.21) %.3f - %.3f; + C_vtx 95 %% %.3f - %.3f"
      % (nm, v, v * plo, v * phi, v * plo * 1.21, v * phi * 1.21, v * plo * v95, v * phi * v95))
P("known-sign items that lower sigma: r_res = 0.89 x0.892, f_eta = 12 % x0.907, y-only reweighting x0.920 -> -8 to -11 %")
P("one-sided items that raise it: live time x sorter x cor_ok +%.1f to +%.1f %%; C_vtx at most +4.2 %% more" % (100 * (plo - 1), 100 * (phi - 1)))
P("HADES 3.5 GeV UrQMD 0.93 +- 0.14 mb against 0.847 mb (side bands + one-sided + pT): %.1f sigma (error of the published value only)" % ((0.93 - 0.847) / 0.14))
P()
P("## 5. vertex loss of like-sign and unlike-sign rows; position of the crossing point")
P("rows removed by vz > -500: e+e- 32.3 %%, e+e+ 70.5 %%, e-e- 75.9 %% -> ratio of the loss fractions %.2f / %.2f" % (70.5 / 32.3, 75.9 / 32.3))
P("   (other definitions: loss odds %.1f / %.1f)" % ((0.705 / 0.295) / (0.323 / 0.677), (0.759 / 0.241) / (0.323 / 0.677)))
cell = (-148.4, -98.6); ctr = 0.5 * (cell[0] + cell[1])
for nm, z in (("mode", -197.5), ("mean", -207.1)):
    P("crossing z of Nv, %s %.1f mm: %.1f cm in front of the upstream end of the LH2 cell (%.1f mm), %.1f cm from its centre (%.1f mm)"
      % (nm, z, (cell[0] - z) / 10, cell[0], (ctr - z) / 10, ctr))
P("no-vertex share of the e+e- rows by day span (ve_sim_cause): days 034-045 0.31, 046-053 0.37, 059-067 0.275")
P("target box (ve_data): pion pairs crossing at -170 <= z < -85 mm, r < 8 mm: V 2 182 631 rows, Nv 3339 -> %.2f %% without a vertex" % (100 * 3339 / (2182631 + 3339)))
P("resolution estimate (ve_sim_cause): rows of target events with at most one own track in the finder window 0.0219 (data), 0.0009 (sim);")
P("   conditional vertex loss of such a row 0.33-0.38 (simulation) -> data loss %.4f (upper estimate), simulation %.4f" % (0.0219 * 0.35, 0.0009 * 0.35))
P()
P("## 6. close-pair reading (D32), inputs and their status")
P("direct-decay route a = 0.497 +- 0.059 (stat; weights: X17 fit N5 %.0f %%, side bands per class %.0f %%)"
  % (100 * (1 / 0.065 ** 2) / (1 / 0.065 ** 2 + 1 / 0.139 ** 2), 100 * (1 / 0.139 ** 2) / (1 / 0.065 ** 2 + 1 / 0.139 ** 2)))
P("   systematics of its main input N5 (x17ab): half spread over the background and CB variants 10.9 %; three background classes +22 %;")
P("   on the simulation-study toys of the specification N5 came out +55 % (normalisation (1) failed the slope test there)")
P("conversion route a = 0.237 +- 0.049 (stat)")
P("bump: a_B(all) 0.299, a_B(E2) 0.312 NOT DETERMINED (syst 157-165 %%); f_eta(edge) 0.392 +- 0.050 (determined): %.1f sigma below the direct route"
  % ((0.497 - 0.392) / math.hypot(0.059, 0.050)))
P("f_eta(tag) = 0.454 +- 0.043 (pi0 Dalitz pairs of the same bump decay, M(ee) < 0.10; emulated simulated tag, UNVERIFIED): %.1f sigma from the direct route"
  % ((0.497 - 0.454) / math.hypot(0.059, 0.043)))
P("F-main rho_d x a_eta = 1.3152 x 0.3823 = 0.503 (PDG 2020), 0.488 (PDG 2026): cross-check (closure failed)")
P("E2 share of the peak: data 0.546 +- 0.053 against SMASH 0.425 (C4/C5 shapes of the side-band fit; WP8: 0.423): +%.1f sigma" % ((0.546 - 0.425) / 0.053))
open(OUT, "w").write("\n".join(L) + "\n")
