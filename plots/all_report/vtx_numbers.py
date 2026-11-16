#!/usr/bin/env python3
# ============================================================================================================
# vtx_numbers.py -- integration of the event-vertex efficiency study into report_all (2026-11-14).
# Recomputes, from the published numbers of the verifier ve_verify (notes_all/ve_verify_note.md, tables_all/vtx_*.txt)
# and of the review fixes (tables_all/fix_numbers.txt, fix_x17_fold.txt), every number of report_all that depended on the
# event-vertex efficiency:
#   - the combination of the three data measures S (smooth like-sign CB), T, E of eps_vtx(eta) (inverse-variance mean);
#   - C_vtx = eps_sim / eps_data and its one-sided 95 % CL bound;
#   - the efficiency that would be needed for 0.93 / 1.14 mb and its distance from each measure;
#   - the sigma(pp -> eta X) ranges "with the one-sided items", now WITHOUT the vertex item (the earlier x1.4-1.7 is gone),
#     and with C_vtx at its 95 % CL bound;
#   - the X17 sigma x BR x BR limit with C_vtx (the folded limit scales linearly with sigma, fix_x17_fold.txt);
#   - the close-pair reading (analyst's interpretation): routes, bump scale, E2 shares.
# Output: docs/latex/tables_all/vtx_numbers.txt.   Run: python3 -I plots/all_report/vtx_numbers.py
# ============================================================================================================
import math

OUT = "/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all/vtx_numbers.txt"
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

P("# vtx_numbers.txt -- plots/all_report/vtx_numbers.py; inputs: ve_verify (vtx_ana*.txt, vtx_sim.txt), fix_numbers.txt, fix_x17_fold.txt")
P("## 1. eps_vtx(eta, data) = X_V / (X_V + X_Nv), classes V (vz > -500) and Nv (vz = -1000); base PT3, isBest >= 0, cor_ok")
eS_nom = eps(4657, 397, 902, 411)            # stat_g errors (vtx_ana.txt A/B)
eS_smo = eps(4670, 370, 263, 324)            # smooth like-sign CB, mean of 12 variants (vtx_ana3.txt A)
eS_N = eps(4646, 370, 165, 324)              # N only, no CB (vtx_ana2.txt A)
eT = eps(3297, 110, -4, 77)                  # ECAL-tagged (vtx_ana.txt G)
eE = eps(1241, 76, 21, 80)                   # edge step (vtx_ana.txt H)
eE_tot = math.hypot(eE[1], 0.028)            # + model rms over 18 variations applied to both classes
P("S nominal CB      : %.3f +- %.3f" % eS_nom)
P("S smooth CB       : %.3f +- %.3f   (rms over 12 smooth variants of S_Nv: 30 rows)" % eS_smo)
P("S N only (no CB)  : %.3f +- %.3f" % eS_N)
P("T (ECAL-tagged)   : %.3f +- %.3f" % eT)
P("E (edge step)     : %.3f +- %.3f (stat) +- 0.028 (model) = %.3f total" % (eE[0], eE[1], eE_tot))
m, e, c = comb([eS_smo, eT, (eE[0], eE_tot)])
P("combined S(smooth)+T+E : %.4f +- %.4f, chi2 %.2f/2;  one-sided 95 %% CL lower bound %.3f" % (m, e, c, m - 1.645 * e))
m2, e2, c2 = comb([(eS_nom[0], math.hypot(eS_nom[1], 0.060)), eT, (eE[0], eE_tot)])   # nominal S: stat_g (+) model 0.060 (verifier)
P("combined S(nominal, stat (+) model)+T+E: %.4f +- %.4f, chi2 %.2f/2" % (m2, e2, c2))
m3, e3, c3 = comb([eS_smo, (eE[0], eE_tot)])
P("combined S(smooth)+E (no ECAL): %.4f +- %.4f, chi2 %.2f/1" % (m3, e3, c3))
eps_d, eeps_d = 0.994, 0.021            # adopted (verifier, rounded)
eps_sim = 0.99970                        # C4+C5, all M4 (vtx_sim.txt)
P()
P("## 2. correction factor (not applied)")
C = eps_sim / eps_d
P("C_vtx = eps_sim/eps_data = %.5f/%.3f = %.4f +- %.4f (relative %.1f %%)" % (eps_sim, eps_d, C, C * eeps_d / eps_d, 100 * eeps_d / eps_d))
lo95 = eps_d - 1.645 * eeps_d
P("95 %% CL: eps_data > %.4f -> C_vtx < %.4f (+%.1f %%)" % (lo95, eps_sim / lo95, 100 * (eps_sim / lo95 - 1)))
P("bump: T and E combined 0.999 +- 0.022 (verifier); simulation C2+C3 0.9986 +- 0.0007 (all M4), 0.9981 +- 0.0009 (bump window)")
P("      C_vtx(bump) = %.3f -> 1.00 +- 0.02; ratios peak/bump (DR, BR_G, BR_A): C_vtx(peak)/C_vtx(bump) = 1.00 +- 0.03" % (0.9986 / 0.999))
s0 = 0.429
P("sigma with C_vtx applied: %.4f mb; at the 95 %% CL bound: %.4f mb" % (s0 * C, s0 * eps_sim / lo95))
P()
P("## 3. what eps would give the published values (sigma scales as eps_sim/eps_data)")
for target in (0.93, 1.14):
    need = round(s0 * eps_sim / target, 2)   # 0.46 / 0.38 as in the verifier note
    P("sigma = %.2f mb needs eps_data = %.3f;  distance: T %.1f sigma, E %.1f sigma, S(smooth) %.1f sigma, S(nominal) %.1f sigma, combined %.0f sigma"
      % (target, need, (eT[0] - need) / eT[1], (eE[0] - need) / eE_tot, (eS_smo[0] - need) / eS_smo[1],
         (eS_nom[0] - need) / math.hypot(eS_nom[1], 0.060), (m - need) / e))
P("row survival (not an eta efficiency): all unlike-sign rows 0.6773, cor_ok 0.6787, peak window 0.565")
P()
P("## 4. sigma(pp -> eta X) with the one-sided items, WITHOUT the vertex item (it was x1.4-1.7 before the measurement)")
lt = (1.071, 1.073); so = (1.013, 1.028); co = 1.022
plo = lt[0] * so[0] * co; phi = lt[1] * so[1] * co
P("live time x sorter x cor_ok = %.4f - %.4f" % (plo, phi))
P("adopted 0.429 mb x that = %.3f - %.3f mb" % (s0 * plo, s0 * phi))
dw, dsb = 0.522, 0.621      # direct-decay route: weighted, side bands per class (fix_numbers.txt)
P("direct-decay route (weighted 0.522, side bands 0.621) x that = %.3f - %.3f mb  [lower: weighted x low product; upper: side bands x high product, as in fix_numbers.txt]"
  % (dw * plo, dsb * phi))
P("with the pT-only production item (x1.21) in addition, upper end: %.3f mb" % (dsb * phi * 1.21))
P("with C_vtx at its 95 %% CL bound (x%.3f) in addition: %.3f - %.3f mb (adopted), %.3f - %.3f mb (direct route)"
  % (eps_sim / lo95, s0 * plo * eps_sim / lo95, s0 * phi * eps_sim / lo95, dw * plo * eps_sim / lo95, dsb * phi * eps_sim / lo95))
P("before (fix_numbers.txt): vertex x live time x sorter x cor_ok x1.55-1.92 -> 0.67-0.82 mb; with the direct route 0.81-1.19 mb; with pT up to 1.44 mb")
P()
P("## 5. X17 sigma x BR x BR (normalisation (1), frozen PDG 2020, sigma uncertainty 37.9 % folded; fix_x17_fold.txt)")
ul_f, ul_p = 15.01, 13.14
P("folded 15.01 nb; with C_vtx applied x%.4f: %.2f nb; at the 95 %% CL bound: %.2f nb; plain product 13.14 -> %.2f nb"
  % (C, ul_f * C, ul_f * eps_sim / lo95, ul_p * C))
P("before: 21-26 nb 'if the eta lose their vertex at the N - CB / raw-row rate' (x1.4-1.7): excluded by the measurement")
P()
P("## 6. close-pair reading (analyst's interpretation, to be tested; numbers from fix_numbers.txt, WP5-WP8)")
P("direct-decay route a = 0.497 +- 0.059 (side bands per class 0.592 +- 0.139, X17 fit 0.476 +- 0.065); F-main rho_d x a_eta = 1.3152 x 0.3823 = %.3f (PDG 2020; 0.488 with PDG 2026)" % (1.3152 * 0.3823))
P("conversion route a = 0.237 +- 0.049; bump scale a_B(all) 0.299, a_B(E2) 0.312, f_eta(edge) 0.392 -> 0.30-0.39")
P("ratio direct / conversion = %.2f +- %.2f; direct / bump (0.30-0.39) = %.2f-%.2f" % (0.497 / 0.237, 0.497 / 0.237 * math.hypot(0.059 / 0.497, 0.049 / 0.237), 0.497 / 0.39, 0.497 / 0.30))
P("E2 share of the peak (M(ee) >= 10 MeV): data 0.546 +- 0.053 (WP6) against 0.423 (WP8) / 0.425 (C4/C5 shapes) in SMASH: +%.1f sigma" % ((0.546 - 0.423) / 0.053))
P("HGeant C2 gamma* angle isotropic: simulated C2 efficiency 11-12 %% too high -> bump scale x1.09-1.12 at most: %.2f-%.2f" % (0.299 * 1.09, 0.312 * 1.12))
P("sigma by the direct-decay route 0.52 +- 0.06 mb (stat); with the one-sided items %.2f-%.2f mb" % (dw * plo, dw * phi))
open(OUT, "w").write("\n".join(L) + "\n")
