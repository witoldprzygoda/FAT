#!/usr/bin/env python3
# ============================================================================
# vedgeu_toys.py  --  INDEPENDENT VERIFIER of the package edgeu: pseudo-data with a mis-modelled background
# ============================================================================
# Written 2026-11-15 (prefix vedgeu). Own truth construction (not edgeu_toys.py), from the WP4 recipe as described:
#   truth = sum over slices of the data CB + c B(u) [1 + b X1(u) + q X2(u)] + a (C2c + C3) + a G_eta(u) + a_P (C4 + C5)
#   B = SMASH non-eta background (G_D + G_C) in u, smoothed; X1, X2 = mean of x and x^2 of the background rows per u bin,
#   x = (M4' - 0.55)/0.1, smoothed by a weighted polynomial (own choice, not a running mean); G_eta smoothed.
#   Two own smoothers: (i) exp(pol4) fitted by a scaled-Poisson likelihood on -100..+100 MeV, (ii) Gaussian kernel
#   (12.5 MeV) on the full grid (edge-corrected by the kernel normalisation).
#   Toys: Poisson total N, Poisson like-sign per slice (data counts), CB rebuilt per slice, nominal estimator (own GLS),
#   h_sim = the truth's own bump template (AA). 1000 toys per truth.
# Run: python3 -I plots/all_report/vedgeu_toys.py > /hdd3/przygoda/HADES/FAT_out/all_report/vedgeu/vedgeu_toys_stdout.txt
# ============================================================================
import json
import math
import sys
import time

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import vedgeu_core as V

T0 = time.time()
D, S = V.load()
ud = V.u_data(D); bD = V.ubin(ud); sD = V.slice_index(D["pp"]); s0d = (D["bits"] & 1) != 0
Hk = [V.hist(bD[s0d & (D["k"] == k)], sD[s0d & (D["k"] == k)]) for k in range(3)]
LSp, LSm = Hk[1], Hk[2]
CBs = [V.cb_one(LSp[i], LSm[i])[0] for i in range(4)]
us, mpp = V.u_sim(S); bS = V.ubin(us); sS = V.slice_index(mpp)
s0s = (S["bits"] & 1) != 0; w = np.where(s0s, S["w"], 0.0); g = S["c"]
cf = np.where(S["ee"] < 0.010, V.CF1, V.CF2)
m4p = V.m4p_sim(S); xx = (m4p - 0.55) / 0.1
def T_(gg, ww):
    m = np.isin(g, gg)
    return V.hist(bS[m], sS[m], ww[m]).sum(0)
BUMP = T_([0], w * cf) + T_([1], w); PEAK = T_([2, 3], w); GETA = T_([4], w)
BKG = T_([5, 6], w); BX1 = T_([5, 6], w * xx); BX2 = T_([5, 6], w * xx * xx)
G490 = T_([5, 6], w * np.exp(-0.5 * ((m4p - 0.490) / 0.015) ** 2)); G505 = T_([5, 6], w * np.exp(-0.5 * ((m4p - 0.505) / 0.008) ** 2))
X = V.X
SIG = 0.00670; TT = V.step(SIG)

def smooth_exp(y, order=4, lo=-0.100, hi=0.100):
    m = (X >= lo) & (X < hi); t = (X[m] - 0.5 * (lo + hi)) / (0.5 * (hi - lo)); yy = np.maximum(y[m], 0)
    sc = yy.sum() / max((yy ** 2).sum(), 1e-9) * yy.sum()   # scaled Poisson: effective-count scaling (own choice)
    def nll(p):
        mu = np.exp(np.clip(np.polyval(p, t), -50, 50))
        return float(np.sum(mu - yy * np.log(mu + 1e-300)))
    p = np.polyfit(t, np.log(np.maximum(yy, 1e-3)), order)
    p = minimize(nll, p, method="Nelder-Mead", options=dict(maxiter=40000, xatol=1e-10, fatol=1e-10)).x
    p = minimize(nll, p, method="BFGS").x
    out = y.copy(); out[m] = np.exp(np.polyval(p, t))
    return out
def smooth_kernel(y, s=0.0125):
    K = np.exp(-0.5 * ((X[:, None] - X[None, :]) / s) ** 2)
    return (K @ y) / K.sum(1)
def smooth_ratio(num, den, order=3, lo=-0.100, hi=0.100):
    r = np.where(den > 0, num / np.where(den > 0, den, 1), 0.)
    m = (X >= lo) & (X < hi) & (den > 0)
    p = np.polyfit(X[m], r[m], order, w=np.sqrt(den[m]))
    return np.polyval(p, X)
R1, R2 = smooth_ratio(BX1, BKG), smooth_ratio(BX2, BKG)
R490, R505 = smooth_ratio(G490, BKG, 2), smooth_ratio(G505, BKG, 2)
# the narrow structures are not smooth in u: keep their exact (raw-ratio) shape, kernel-smoothed lightly
R490k = smooth_kernel(np.where(BKG > 0, G490 / np.where(BKG > 0, BKG, 1), 0), 0.005)
R505k = smooth_kernel(np.where(BKG > 0, G505 / np.where(BKG > 0, BKG, 1), 0), 0.005)
SM = {"exp4": (smooth_exp(BKG, 4), smooth_exp(GETA, 4)), "kern": (smooth_kernel(BKG), smooth_kernel(GETA))}

def estimator(Ntot, pp, mm, Vfix=None):
    CB = np.zeros(V.NB); COV = np.zeros((V.NB, V.NB))
    for i in range(4):
        c, cv = V.cb_one(pp[i], mm[i]); CB += c; COV += cv
    return V.gls_fit(Ntot, CB, COV, TT, fixed=V.AP * PEAK, Vfix=Vfix)
Rdat = estimator(Hk[0].sum(0), LSp, LSm)
HS = V.proj(Rdat, BUMP)
print("data: h %.2f, h_sim %.2f, f %.4f" % (Rdat.h, HS, Rdat.h / HS))

def truth(sm="exp4", c=1.35, b=-0.15, q=0.0, a=0.42, A=0.0, Rg=None):
    B, Ge = SM[sm]
    bkg = c * B * (1 + b * R1 + q * R2 + (A * Rg if Rg is not None else 0.0))
    return sum(CBs) + bkg + a * BUMP + a * Ge + V.AP * PEAK
SC = [("a", dict()), ("q-", dict(q=-0.173)), ("q--", dict(q=-0.346)), ("q+", dict(q=0.173)),
      ("narrow", dict(A=0.05, Rg=R505k)), ("wide", dict(A=0.10, Rg=R490k)),
      ("a_kern", dict(sm="kern")), ("q-_kern", dict(sm="kern", q=-0.173))]
rng = np.random.Generator(np.random.PCG64(0x70E5))
NT = 1000
OUT = {}
for tag, kw in SC:
    mu = truth(**kw)
    Ra = estimator(mu, LSp, LSm)
    fas = Ra.h / V.proj(Ra, BUMP)
    if tag in ("a", "q-", "a_kern", "q-_kern"):
        fs, pulls = [], []
        for it in range(NT):
            Nt = rng.poisson(mu).astype(float)
            pp = [rng.poisson(LSp[i]).astype(float) for i in range(4)]; mm = [rng.poisson(LSm[i]).astype(float) for i in range(4)]
            R = estimator(Nt, pp, mm)
            hs = V.proj(R, BUMP); fs.append(R.h / hs); pulls.append((R.h / hs - 0.42) / (R.eh / hs))
        fs = np.array(fs); pulls = np.array(pulls)
        OUT[tag] = dict(asimov=fas, bias_asimov=100 * (fas / 0.42 - 1), toy_mean=float(fs.mean()), toy_err=float(fs.std(ddof=1) / math.sqrt(NT)),
                        bias_toy=100 * (fs.mean() / 0.42 - 1), pull_sd=float(pulls.std(ddof=1)), pull_mean=float(pulls.mean()))
        print("%-8s Asimov f %.4f (%+.2f %%); toys mean %.4f +- %.4f (%+.2f %%), pull mean %.2f sd %.3f" %
              (tag, fas, OUT[tag]["bias_asimov"], fs.mean(), OUT[tag]["toy_err"], OUT[tag]["bias_toy"], pulls.mean(), pulls.std(ddof=1)))
    else:
        OUT[tag] = dict(asimov=fas, bias_asimov=100 * (fas / 0.42 - 1))
        print("%-8s Asimov f %.4f (%+.2f %%)" % (tag, fas, OUT[tag]["bias_asimov"]))
for tag in OUT:
    ref = "a_kern" if tag.endswith("_kern") else "a"
    OUT[tag]["shift_vs_a_asimov"] = OUT[tag]["bias_asimov"] - OUT[ref]["bias_asimov"]
    if "bias_toy" in OUT[tag]:
        OUT[tag]["shift_vs_a_toy"] = OUT[tag]["bias_toy"] - OUT[ref]["bias_toy"]
    print("   %-8s shift vs (a): Asimov %+.2f %%%s" % (tag, OUT[tag]["shift_vs_a_asimov"], ("; toys %+.2f %%" % OUT[tag]["shift_vs_a_toy"]) if "shift_vs_a_toy" in OUT[tag] else ""))
# the q-term in M4 (B4-like): how much of c q B X2 is absorbed by pol2 + step? fraction of the q-term in the fit range that the step takes
dq = 1.35 * (-0.173) * SM["exp4"][0] * R2
print("q-term: change of the truth in -45..+45 MeV: %.0f rows (%.1f %% of N there); step projection %.1f rows (%.2f %% of a h_sim)" %
      (dq[Rdat.idx].sum(), 100 * dq[Rdat.idx].sum() / Hk[0].sum(0)[Rdat.idx].sum(), V.proj(Rdat, dq), 100 * V.proj(Rdat, dq) / (0.42 * HS)))
# the smooth-background step projections (U11-like): background and G_eta parts with both smoothers, and raw
for sm in ("exp4", "kern"):
    B, Ge = SM[sm]
    print("smoother %s: step projection of 1.35 B (1 - 0.15 X1) %.1f, of 0.42 G_eta %.1f -> bias %+.2f %%" %
          (sm, V.proj(Rdat, 1.35 * B * (1 - 0.15 * R1)), V.proj(Rdat, 0.42 * Ge), 100 * (V.proj(Rdat, 1.35 * B * (1 - 0.15 * R1)) + V.proj(Rdat, 0.42 * Ge)) / (0.42 * HS)))
for o in (2, 3, 5, 6):
    B = smooth_exp(BKG, o); Ge = smooth_exp(GETA, o)
    print("exp(pol%d): step projection of B %.1f, of G_eta %.1f" % (o, V.proj(Rdat, B), V.proj(Rdat, Ge)))
print("raw: B %.1f, G_eta %.1f" % (V.proj(Rdat, BKG), V.proj(Rdat, GETA)))
json.dump(OUT, open(V.DIR + "vedgeu_toys.json", "w"), indent=1)
print("done in %.0f s" % (time.time() - T0))
