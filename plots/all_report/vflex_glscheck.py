#!/usr/bin/env python3
# vflex_glscheck.py -- independent verification of package flex: is the hybrid GLS chi2 a valid model-comparison statistic?
#   The hybrid GLS covariance W(theta) contains terms that grow with the scales (BB-lite diagonal ~ s^2, replica covariance of
#   the smoothed part ~ s s'). (1) evaluate the F8 hybrid GLS at the parameter points of the nested tied fits; (2) the same
#   with W frozen at the F8 fit point (one common metric); (3) minimise the hybrid GLS of F8 directly over theta.
# Run: python3 -B plots/all_report/vflex_glscheck.py
import sys, json
import numpy as np
from scipy import optimize, linalg
sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import vflex_core as V

D = V.load(); L = V.Layout(D, (1,))
M = V.FlexModel(["C2", "C3", "C4", "C5", "Geta", "GD_om", "GD_nom", "GC"])
T, Vt = V.templates(D, L, M)
TR = [V.templates(D, L, M, rep=r)[0] for r in range(V.NREPS)]
R8 = V.fit(M, L, T, Vt, egcov=L.eg)
ties = json.load(open(V.DIR + "/vflex_ties.json"))


def W_of(th):
    mu, J, vbb, k, _ = M.eval(th, L, T, Vt)
    W = V.data_cov(L.y, mu, L.eg, floor=1.0) + k * k * L.ccb + np.diag(vbb)
    B = np.array([M.eval(th, L, Tr, Vt)[4] for Tr in TR]); B = B - B.mean(0)
    return W + B.T @ B / (len(TR) - 1), mu


def gls(th, W=None):
    Wt, mu = W_of(th)
    if W is None:
        W = Wt
    r = L.y - mu
    return r @ linalg.solve(W, r, assume_a='pos') + ((th[-1] - 1) / 0.05) ** 2


W8, _ = W_of(R8['th'])
print(f"F8 fit point: hybrid GLS {gls(R8['th']):.2f}; frozen-W {gls(R8['th'], W8):.2f}; log det W {np.linalg.slogdet(W8)[1]:.2f}")
out = {}
for tag, v in ties.items():
    if "th_full" not in v:
        continue
    th = np.array(v["th_full"])
    Wt, _ = W_of(th)
    out[tag] = dict(hyb=gls(th), frozen=gls(th, W8), dlogdet=np.linalg.slogdet(Wt)[1] - np.linalg.slogdet(W8)[1])
    print(f"[{tag:28s}] hybrid GLS {out[tag]['hyb']:8.2f} | with W frozen at the F8 point {out[tag]['frozen']:8.2f} (d {out[tag]['frozen'] - gls(R8['th'], W8):+7.2f}) | d log det W {out[tag]['dlogdet']:+6.2f}")
# direct minimisation of the hybrid GLS over the F8 parameters
r = optimize.minimize(gls, R8['th'], method="Powell", options=dict(maxiter=20000, xtol=1e-6, ftol=1e-9))
print(f"minimum of the F8 hybrid GLS over theta: {r.fun:.2f} at " + " ".join(f"{n}={x:.3f}" for n, x in zip(M.names, r.x)))
out["min_hybrid"] = dict(val=r.fun, th=r.x.tolist())
json.dump(out, open(V.DIR + "/vflex_glscheck.json", "w"), indent=1, default=float)
