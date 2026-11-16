#!/usr/bin/env python3
# ============================================================================
# vedgeu_checks.py  --  INDEPENDENT VERIFIER of the package edgeu: further checks with own code
# ============================================================================
# Written 2026-11-15 (prefix vedgeu). Stat of the q = 0 alignment and of pol3; E1/E2 classes and the kappa split;
# four-slice estimator in u; periods; the e-e- like-sign step (class shares, period pattern, shape).
# Run: python3 -I plots/all_report/vedgeu_checks.py > /hdd3/przygoda/HADES/FAT_out/all_report/vedgeu/vedgeu_checks_stdout.txt
# ============================================================================
import json
import math
import sys

import numpy as np
from scipy.stats import chi2 as C2

sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import vedgeu_core as V

OUT = {}
D, S = V.load()
s0d = (D["bits"] & 1) != 0; s0s = (S["bits"] & 1) != 0
ud = V.u_data(D); bD = V.ubin(ud); sD = V.slice_index(D["pp"])
us, mpp = V.u_sim(S); bS = V.ubin(us); sS = V.slice_index(mpp)
w = np.where(s0s, S["w"], 0.0); g = S["c"]; cf = np.where(S["ee"] < 0.010, V.CF1, V.CF2)
T = V.step(0.0067)
def dH(mask, b=bD, s=sD):
    return [V.hist(b[mask & (D["k"] == k)], s[mask & (D["k"] == k)]) for k in range(3)]
def sT(gg, mask=None, b=bS, s=sS, c2=True):
    m = np.isin(g, gg) & (np.ones(len(g), bool) if mask is None else mask)
    ww = w * np.where((g == 0) & c2, cf, 1.0)
    return V.hist(b[m], s[m], ww[m])          # per slice
H = dH(s0d)
N, CB, CBC = V.spectrum(*H)
BUMPs = sT([0, 1]); PEAKs = sT([2, 3])
R0 = V.gls_fit(N, CB, CBC, T, fixed=V.AP * PEAKs.sum(0))
hs0 = V.proj(R0, BUMPs.sum(0))
print("nominal: h %.2f +- %.2f, rel %.2f %%; f %.4f" % (R0.h, R0.eh, 100 * R0.eh / R0.h, R0.h / hs0))
# ---------------------------------------------------------------- stat with E(m, 0) and with pol3
ud0 = V.u_data(D, useq=False); b0 = V.ubin(ud0)
us0, _ = V.u_sim(S, useq=False); bs0 = V.ubin(us0)
N0, CB0, CBC0 = V.spectrum(*dH(s0d, b=b0))
Rq = V.gls_fit(N0, CB0, CBC0, T, fixed=V.AP * sT([2, 3], b=bs0).sum(0))
hq = V.proj(Rq, sT([0, 1], b=bs0).sum(0))
R3 = V.gls_fit(N, CB, CBC, T, fixed=V.AP * PEAKs.sum(0), order=3)
OUT["stat_rel"] = dict(nominal=R0.eh / R0.h, q0=Rq.eh / Rq.h, pol3=R3.eh / R3.h, pol3_over_pol2=(R3.eh / V.proj(R3, BUMPs.sum(0))) / (R0.eh / hs0))
print("relative stat (fit): nominal %.2f %%, E(m, 0) %.2f %% (f %.4f), pol3 %.2f %% (error ratio pol3/pol2 on f %.2f)" %
      (100 * R0.eh / R0.h, 100 * Rq.eh / Rq.h, Rq.h / hq, 100 * R3.eh / R3.h, OUT["stat_rel"]["pol3_over_pol2"]))
# ---------------------------------------------------------------- M(ee) classes (CB = CB(all) x like-sign class share in -100..+100 MeV, per slice)
xm = (V.X >= -0.100) & (V.X < 0.100)
def class_spec(cm):
    Hc = dH(s0d & cm)
    Nc = Hc[0].sum(0); CBc = np.zeros(V.NB); COVc = np.zeros((V.NB, V.NB))
    for i in range(4):
        cb, cv = V.cb_one(H[1][i], H[2][i])
        den = (H[1][i] + H[2][i])[xm].sum()
        sh = (Hc[1][i] + Hc[2][i])[xm].sum() / den
        CBc += cb * sh; COVc += cv * sh * sh
    return Nc, CBc, COVc
CL = {}
for lab, cmd, cms in [("E1", D["ee"] < 0.010, S["ee"] < 0.010), ("E2", D["ee"] >= 0.010, S["ee"] >= 0.010)]:
    Nc, CBc, COVc = class_spec(cmd)
    Rc = V.gls_fit(Nc, CBc, COVc, T, fixed=V.AP * sT([2, 3], mask=cms).sum(0))
    c2c = V.proj(Rc, sT([0], mask=cms).sum(0)); c3 = V.proj(Rc, sT([1], mask=cms).sum(0))
    CL[lab] = dict(h=Rc.h, eh=Rc.eh, c2c=c2c, c3=c3, f=Rc.h / (c2c + c3))
    print("class %s: h %.1f +- %.1f, C2c %.1f, C3 %.1f, f %.4f +- %.4f (fit)" % (lab, Rc.h, Rc.eh, c2c, c3, Rc.h / (c2c + c3), Rc.eh / (c2c + c3)))
e1, e2 = CL["E1"], CL["E2"]
rho = e1["h"] / e2["h"]
kap = (rho * e2["c2c"] - e1["c2c"]) / (e1["c3"] - rho * e2["c3"]); aD = e2["h"] / (e2["c2c"] + kap * e2["c3"])
# stat of kappa by error propagation of the two independent class fits
def ksol(h1, h2):
    r = h1 / h2
    k = (r * e2["c2c"] - e1["c2c"]) / (e1["c3"] - r * e2["c3"])
    return k, h2 / (e2["c2c"] + k * e2["c3"])
dk1 = (ksol(e1["h"] + e1["eh"], e2["h"])[0] - kap); dk2 = (ksol(e1["h"], e2["h"] + e2["eh"])[0] - kap)
OUT["classes"] = CL; OUT["kappa"] = dict(kappa=kap, a=aD, stat=math.hypot(dk1, dk2))
print("E1/E2 split: a_D %.3f, kappa %.3f +- %.3f (stat, linear propagation)" % (aD, kap, math.hypot(dk1, dk2)))
# ---------------------------------------------------------------- four-slice estimator in u
hd4, hs4 = 0.0, 0.0
for i in range(4):
    cb, cv = V.cb_one(H[1][i], H[2][i])
    Ri = V.gls_fit(H[0][i].astype(float), cb, cv, T, fixed=V.AP * PEAKs[i])
    hd4 += Ri.h; hs4 += V.proj(Ri, BUMPs[i])
    print("   slice %d: h_data %.1f +- %.1f, h_sim %.1f, f %.4f" % (i, Ri.h, Ri.eh, V.proj(Ri, BUMPs[i]), Ri.h / V.proj(Ri, BUMPs[i])))
OUT["four_slice"] = dict(hd=hd4, hs=hs4, f=hd4 / hs4)
print("four-slice estimator: sum h_data %.1f, sum h_sim %.1f, f %.4f" % (hd4, hs4, hd4 / hs4))
# ---------------------------------------------------------------- periods
per = np.full(len(D["fid"]), -1)
per[(D["fid"] >= 34) & (D["fid"] <= 46)] = 0; per[(D["fid"] >= 47) & (D["fid"] <= 48)] = 1
per[(D["fid"] >= 49) & (D["fid"] <= 53)] = 2; per[(D["fid"] >= 59) & (D["fid"] <= 67)] = 3
LP = [2809.453, 372.610, 1078.022, 1819.110]; CS = [1.0241, 1.0, 1.0, 1.0]
fp, ep = [], []
for p in range(4):
    Np_, CBp, CBCp = V.spectrum(*dH(s0d & (per == p)))
    Rp = V.gls_fit(Np_, CBp, CBCp, T, fixed=V.AP * PEAKs.sum(0) * LP[p] / 6079.2)
    hsp = V.proj(Rp, BUMPs.sum(0)) * LP[p] / 6079.2
    fp.append(Rp.h / hsp * CS[p]); ep.append(Rp.eh / hsp * CS[p])
fp, ep = np.array(fp), np.array(ep)
wm = np.sum(fp / ep ** 2) / np.sum(1 / ep ** 2); ch = float(np.sum((fp - wm) ** 2 / ep ** 2))
OUT["periods"] = dict(f=fp.tolist(), e=ep.tolist(), chi2=ch, p=float(1 - C2.cdf(ch, 3)))
print("periods (x C_sec): " + ", ".join("%.3f +- %.3f" % (a, b) for a, b in zip(fp, ep)) + "; chi2 %.2f/3, p %.3f" % (ch, 1 - C2.cdf(ch, 3)))
# ---------------------------------------------------------------- the e-e- step
mmT = H[2].sum(0); ppT = H[1].sum(0)
lin = lambda Hh: (V.proj(R0, Hh), math.sqrt(float(R0.Lh @ np.diag(np.maximum(Hh[R0.idx], 1.0)) @ R0.Lh)))
for lab, cm in [("E1", D["ee"] < 0.010), ("E2", D["ee"] >= 0.010)]:
    Hm = dH(s0d & cm)[2].sum(0); Hp = dH(s0d & cm)[1].sum(0)
    fr = Hm[R0.idx].sum() / mmT[R0.idx].sum()
    hm, em = lin(Hm)
    print("e-e- class %s: %.1f %% of the e-e- rows in -45..+45; step %.1f +- %.1f; expected if the step were uniform in the e-e- rows: %.1f" %
          (lab, 100 * fr, hm, em, fr * lin(mmT)[0]))
    OUT["mm_" + lab] = dict(frac=fr, step=hm, err=em, expect_uniform=fr * lin(mmT)[0])
hm_p, em_p = [], []
for p in range(4):
    a, b = lin(dH(s0d & (per == p))[2].sum(0)); hm_p.append(a); em_p.append(b)
hm_p, em_p = np.array(hm_p), np.array(em_p)
tot = hm_p.sum(); Lw = np.array(LP) / sum(LP)
# share of the e-e- rows per period (the step relative to the e-e- rows)
nmm = np.array([dH(s0d & (per == p))[2].sum(0)[R0.idx].sum() for p in range(4)])
for lab, wts in [("L", Lw), ("e-e- rows", nmm / nmm.sum())]:
    exp_ = tot * wts
    chp = float(np.sum((hm_p - exp_) ** 2 / em_p ** 2))
    print("e-e- step per period %s; expected proportional to %s: %s; chi2 %.2f/3, p %.3f" %
          (", ".join("%.1f +- %.1f" % (a, b) for a, b in zip(hm_p, em_p)), lab, ", ".join("%.1f" % x for x in exp_), chp, 1 - C2.cdf(chp, 3)))
    OUT["mm_periods_" + lab] = dict(step=hm_p.tolist(), err=em_p.tolist(), expected=exp_.tolist(), chi2=chp, p=float(1 - C2.cdf(chp, 3)))
# shape of the e-e- spectrum near u = 0 (rows per 2.5 MeV) with a smooth pol2 from the outer bins
print("e-e- rows per 2.5 MeV bin, -45..+45 MeV:")
print("  " + " ".join("%d" % v for v in mmT[R0.idx]))
print("e+e+ rows:")
print("  " + " ".join("%d" % v for v in ppT[R0.idx]))
json.dump(OUT, open(V.DIR + "vedgeu_checks.json", "w"), indent=1, default=float)
