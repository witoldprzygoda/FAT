#!/usr/bin/env python3
# ============================================================================
# peaktail_summary.py  --  package peaktail, task 4: the effect of the tail scale t on the other Part III quantities and
#                          the collection of all numbers (peaktail_numbers.json)
# ============================================================================
# Written 2026-11-15 (prefix peaktail). Reads the json of peaktail_sim / _ufit / _line / _rad / _prop and peaktail_bedge.txt.
# Computes in addition (simulated rows of peaktail_dump.C, S0, mode A, w x V3):
#   - the peak side bands (WP6 estimator: Y = [(N-CB)_P - R'(N-CB)_SB]/(1 - R' f), R' = 1.1859) and a_P, sigma(peak) for
#     t = t_data +- 1 sigma and the QED tail: only the leakage f = (C4+C5)(SB)/(C4+C5)(P) and the P fraction change;
#   - the double ratio DR = a_P / f_eta(u-edge);
#   - a proxy of the ECAL tag (a Gaussian content with a pol2, as the frozen tag fit): the frozen Gaussian (0.5439, 6.09 MeV)
#     + pol2 on 0.49-0.62 fitted to the simulated peak C4 + C5 in M4' with the tail scaled by t;
#   - the X17 fit, normalisation (1) (X17 rows = r N5 R_X): the change of R_X = X17-like / direct-decay content in the X2D
#     cells (M4 0.5225-0.61, M(ee) 2-62 MeV), with C5 rows of M(ee) 15.5-24.5 MeV standing for the X17 line.
# Output: /hdd3/.../peaktail/peaktail_numbers.json, peaktail_summary.log
# Run: python3 -I plots/all_report/peaktail_summary.py
# ============================================================================
import json
import math
import sys

import numpy as np

sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import edgeu_core as C
import peaktail_core as P

LOG = open(P.OUT + "peaktail_summary.log", "w")
def L(*a):
    s = " ".join(str(x) for x in a)
    print(s); LOG.write(s + "\n"); LOG.flush()
J = {k: json.load(open(P.OUT + "peaktail_%s.json" % k)) for k in ("sim", "ufit", "line", "rad", "prop")}
BE = {}
for l in open(P.OUT + "peaktail_bedge.txt"):
    if l.startswith("PT "):
        x = l.split(); tag = x[1]
        vals = dict(zip(x[2::2], [float(v) for v in x[3::2]]))
        BE.setdefault(tag, vals)        # first occurrence = WP6 weights (nominal bedge)
OUTN = {}
td, sd = J["prop"]["t_comb"]["t_data"], J["prop"]["t_comb"]["s_data"]
tq5 = J["prop"]["t_qed"]["t_C5"]
S = P.load_sim(); w = P.sim_weight(S); gS = S["g"]
us, mpp, m4p = C.u_sim(S)
M4r, Jac = P.m4_and_jacobian(S)
dl = C.S_NOM * (Jac[0] * (S["ep_p"] - S["ep_ptrue"]) + Jac[1] * (S["em_p"] - S["em_ptrue"])).astype(np.float64) * 1e-3
flag = (dl < -0.005)
RR = np.load(P.OUT + "peaktail_rad_rows.npz")
pk = ((gS == 2) | (gS == 3)) & (w > 0)
def wt_t(t4, t5, rng=None):
    out = w.copy()
    for gg, t in ((2, t4), (3, t5)):
        m = (gS == gg) & (w > 0)
        if rng is not None:
            m &= (mpp >= rng[0]) & (mpp < rng[1])
        N = w[m].sum(); Nt = (w * flag)[m].sum(); Nc = N - Nt
        cc = (N - t * Nt) / Nc
        out[m] = w[m] * ((~flag[m]) * cc + t * flag[m])
    return out
# ------------------------------------------------------------------ side bands, a_P, sigma(peak), DR
PW = (0.53364, 0.55416); SBW = ((0.5250, 0.53364), (0.55416, 0.5628)); Rp = 1.1859
def win(wt, x, sel):
    inP = (x >= PW[0]) & (x < PW[1]); inSB = ((x >= SBW[0][0]) & (x < SBW[0][1])) | ((x >= SBW[1][0]) & (x < SBW[1][1]))
    return wt[sel & inP].sum(), wt[sel & inSB].sum(), wt[sel].sum()
P1, SB1, N1 = win(w, m4p, pk); f1 = SB1 / P1
def sbq(Pt, SBt, Nt):
    ft = SBt / Pt; yr = (1 - Rp * f1) / (1 - Rp * ft); aPr = yr * P1 / Pt
    return dict(f=ft, Y_rel=yr, aP_rel=aPr, sigma_rel=aPr, Pfrac=Pt / Nt)
SBR = {}
for tag, t4, t5 in (("t_data", td, td), ("t_data+1s", td + sd, td + sd), ("t_data-1s", td - sd, td - sd), ("t0", 0, 0), ("t2", 2, 2), ("qed", 1, tq5)):
    SBR[tag] = sbq(*win(wt_t(t4, t5), m4p, pk))
# QED toy rows
w4 = np.where(gS == 2, w, 0.0)
Pa, SBa, Na = win(w4, m4p, pk & (gS == 2))
PW_r = (RR["m4p"] >= PW[0]) & (RR["m4p"] < PW[1])
SB_r = ((RR["m4p"] >= SBW[0][0]) & (RR["m4p"] < SBW[0][1])) | ((RR["m4p"] >= SBW[1][0]) & (RR["m4p"] < SBW[1][1]))
SBR["qed_rows"] = sbq(Pa + RR["w"][PW_r].sum(), SBa + RR["w"][SB_r].sum(), Na + RR["w"].sum())
L("PEAK SIDE BANDS (leakage f = SB/P of the simulated peak; Y_peak, a_P, sigma(peak) relative to the nominal; nominal f = %.4f)" % f1)
for k, v in SBR.items():
    L("   %-10s f %.4f P/all %.4f  Y_peak x %.4f  a_P = sigma(peak) x %.4f (%+.1f %%)" % (k, v["f"], v["Pfrac"], v["Y_rel"], v["aP_rel"], 100 * (v["aP_rel"] - 1)))
slope_aP = (SBR["t_data+1s"]["aP_rel"] - SBR["t_data-1s"]["aP_rel"]) / (2 * sd)
OUTN["sidebands"] = SBR; OUTN["aP_slope_per_t"] = slope_aP
L("   slope of a_P (and sigma from the peak) around t_data: %+.3f per unit t; +-1 sigma(t): %+.1f / %+.1f %%" %
  (slope_aP, 100 * (SBR["t_data+1s"]["aP_rel"] / SBR["t_data"]["aP_rel"] - 1), 100 * (SBR["t_data-1s"]["aP_rel"] / SBR["t_data"]["aP_rel"] - 1)))
FU = J["prop"]["f_u"]
DR = {}
for k in ("t_data", "t_data+1s", "t_data-1s", "t0", "t2", "qed"):
    DR[k] = SBR[k]["aP_rel"] / (FU[k]["f"] / FU["nominal"]["f"])
L("   double ratio DR_u = a_P / f(u-edge), relative to the nominal 0.98: " + ", ".join("%s x %.3f" % (k, v) for k, v in DR.items()))
OUTN["DR_rel"] = DR
# ------------------------------------------------------------------ ECAL-tag proxy: Gaussian content with pol2
xe = np.arange(0.49, 0.62 + 1e-9, 0.0025); xc = 0.5 * (xe[1:] + xe[:-1])
from scipy.special import erf
def gfrac(lo, hi, mu, s):
    return 0.5 * (erf((hi - mu) / (math.sqrt(2) * s)) - erf((lo - mu) / (math.sqrt(2) * s)))
Gc = gfrac(xe[:-1], xe[1:], P.MPEAK, P.SPEAK_S)
def gcontent(wt):
    h = np.histogram(m4p[pk], bins=xe, weights=wt[pk])[0]
    tt = (xc - 0.555) / 0.065
    F = np.vstack([np.ones_like(tt), tt, tt * tt, Gc]).T
    p, *_ = np.linalg.lstsq(F, h, rcond=None)   # equal weights: in the data the background dominates the bin errors (S/B 0.1-0.2)
    return p[3], h.sum()
A1, H1 = gcontent(w)
TAG = {}
for tag, t4, t5 in (("t0", 0, 0), ("t_data-1s", td - sd, td - sd), ("t_data", td, td), ("t_data+1s", td + sd, td + sd), ("t2", 2, 2), ("qed", 1, tq5)):
    A_, H_ = gcontent(wt_t(t4, t5))
    TAG[tag] = dict(content_rel=A_ / A1, ftag_rel=A1 / A_)
L("\nECAL-TAG PROXY (Gaussian 0.5439 / 6.09 MeV + pol2 on 0.49-0.62 fitted to the simulated peak; nominal content %.3f of the rows in the range)" % (A1 / H1))
for k, v in TAG.items():
    L("   %-10s Gaussian content x %.4f -> f(tag) = T_data / T_sim x %.4f (%+.1f %%)" % (k, v["content_rel"], v["ftag_rel"], 100 * (v["ftag_rel"] - 1)))
OUTN["tag_proxy"] = TAG
# ------------------------------------------------------------------ X17 normalisation (1)
x2d = (m4p >= 0.5225) & (m4p < 0.61) & (S["mee"] >= 0.002) & (S["mee"] < 0.062)
c5 = (gS == 3) & (w > 0)
lineX = c5 & (S["mee"] >= 0.0155) & (S["mee"] < 0.0245)
band = (m4p >= 0.5325) & (m4p < 0.5575)
def rx(wt):
    return (wt[lineX & band & x2d].sum() / wt[lineX].sum()) / (wt[c5 & x2d].sum() / wt[c5].sum())
rx0 = rx(w)
X17 = {}
for tag, t5 in (("t0", 0.0), ("t_data-1s", td - sd), ("t_data+1s", td + sd), ("t2", 2.0), ("qed", tq5)):
    X17[tag] = rx(wt_t(1, t5)) / rx0
# QED toy rows: X17-like rows and C5 in X2D
mr = RR["rows"]; m4r = RR["m4p"]; meer = RR["mee"]; wr = RR["w"]
x2dr = (m4r >= 0.5225) & (m4r < 0.61) & (meer >= 0.002) & (meer < 0.062)
lineXr = (S["mee"][mr] >= 0.0155) & (S["mee"][mr] < 0.0245)
bandr = (m4r >= 0.5325) & (m4r < 0.5575) & (meer >= 0.0155) & (meer < 0.0245)
rxq = (wr[lineXr & bandr & x2dr].sum() / wr[lineXr].sum()) / (wr[x2dr].sum() / wr.sum())
X17["qed_rows"] = rxq / rx0
OUTN["x17_RX_rel"] = X17
L("\nX17, normalisation (1): R_X = (X17-like rows in the peak band of the line) / (direct decays in the X2D cells), relative:")
L("   " + ", ".join("%s x %.4f" % (k, v) for k, v in X17.items()) + "   (the limit on BR scales with 1/R_X)")
# ------------------------------------------------------------------ collect
OUTN["t_data"] = dict(t=td, s=sd); OUTN["t_qed"] = J["prop"]["t_qed"]; OUTN["t_comb"] = J["prop"]["t_comb"]
OUTN["f_u"] = {k: dict(f=v["f"], rel=v["rel"]) for k, v in FU.items()}
OUTN["f_u_slope"] = J["prop"]["f_u_slope_rel_per_t"]
nb = BE["nominal"]["fc_rec"]
OUTN["f_bedge"] = {k: dict(f=v["fc_rec"], rel=v["fc_rec"] / nb - 1, f_all=v["f_all"], bumponly=v["bumponly"]) for k, v in BE.items()}
OUTN["f_bedge_slope"] = (BE["t2"]["fc_rec"] - BE["t0"]["fc_rec"]) / 2 / nb
OUTN["syst_u"] = abs(OUTN["f_u_slope"]) * sd
OUTN["syst_bedge"] = abs(OUTN["f_bedge_slope"]) * sd
OUTN["checker"] = J["ufit"]["checker"]; OUTN["decomp"] = J["ufit"]["decomp"]
OUTN["t_u"] = J["prop"]["t_u"]; OUTN["t_ctrl"] = J["prop"]["t_ctrl"]; OUTN["t_u_class"] = J["prop"]["t_u_class"]
L("\nf(u-edge) slope %+.4f / unit t, systematic +-%.2f %%; f(bedge) slope %+.4f / unit t, systematic +-%.2f %%" %
  (OUTN["f_u_slope"], 100 * OUTN["syst_u"], OUTN["f_bedge_slope"], 100 * OUTN["syst_bedge"]))
for k in ("t_data", "t_data+1s", "t_data-1s", "t0", "t2", "qed"):
    L("   %-10s f(u) %.4f (%+.2f %%)   f(bedge) %.4f (%+.2f %%)   difference of the relative shifts %+.2f %%" %
      (k, FU[k]["f"], 100 * FU[k]["rel"], BE[k]["fc_rec"], 100 * (BE[k]["fc_rec"] / nb - 1), 100 * (FU[k]["rel"] - (BE[k]["fc_rec"] / nb - 1))))
json.dump(OUTN, open(P.OUT + "peaktail_numbers.json", "w"), indent=1)
L("\nwritten: peaktail_numbers.json, peaktail_summary.log")
