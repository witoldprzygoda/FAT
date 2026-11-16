#!/usr/bin/env python3
# ============================================================================
# edgeu_bgsmooth.py  --  package edgeu: goodness of the smoothers of the simulated non-eta background in u (for the E7-like term)
# ============================================================================
# Written 2026-11-15 (prefix edgeu). The step projection of a smooth background (the part that pol2 cannot absorb) depends on
# the smoother. This script fits exp(pol n), n = 2..6, on u in [-150, 150) and [-100, 100) MeV to the simulated non-eta
# background (G_D + G_C, S0, mode A, V3, 0.30 <= M(pipi)' < 0.40), with the event-group variance, and reports the chi2 of each
# smoother next to its step projection (estimator of the nominal data fit). Output: /hdd3/.../edgeu/edgeu_bgsmooth.json.
# Run: python3 -I plots/all_report/edgeu_bgsmooth.py
# ============================================================================
import json
import math
import sys

import numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2 as chi2dist

sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import edgeu_core as C

ANA = json.load(open(C.OUT + "edgeu_ana.json"))
D = C.load_data(); S = C.load_sim()
g = C.GRID_NOM
SL = C.slices_of(0.30, 0.40)
ud = C.u_data(D)
Hd = C.hist_data(D, ud, (D["bits"] & 1) != 0)
sp = C.data_spec(Hd, SL)
u, mpp, m4p = C.u_sim(S)
w = C.sim_weight(S, c2corr=True)
Hs = C.hist_sim(S, u, mpp, w, g)
Vs = C.hist_sim_var(S, u, mpp, w, g)
def comp(H, names):
    return sum(H[C.G[n], a:b].sum((0, 1)) for n in names for a, b in SL)
peak = comp(Hs, ["C4", "C5"])
T = C.step_gauss(g, ANA["reso"]["sigma_core"])
R0 = C.fit_step(sp, T, g, -0.045, 0.045, 2, fixed=C.AP * peak)
hs = C.h_lin(R0, comp(Hs, ["C2", "C3"]))
out = {}
for lab, names in (("nonEta", ["GDom", "GDnom", "GCom", "GCnom"]), ("Geta", ["Geta"])):
    y = comp(Hs, names); v = comp(Vs, names)
    for lo, hi in ((-0.150, 0.150), (-0.100, 0.100)):
        m = (g.x >= lo) & (g.x < hi)
        xx = (g.x[m] - 0.5 * (lo + hi)) / (0.5 * (hi - lo)); yy = np.maximum(y[m], 0); vv = np.maximum(v[m], 1e-6)
        kfac = vv.sum() / max(yy.sum(), 1e-9)   # scaled Poisson: effective counts y / k
        def nll(p):
            mu = np.exp(np.clip(np.polyval(p[::-1], xx), -50, 50))
            return float(np.sum(mu - yy * np.log(mu + 1e-300)))
        p = np.array([math.log(max(yy.mean(), 1e-6))])
        for o in range(1, 7):
            p = minimize(nll, np.r_[p, 0.0], method="BFGS").x
            if o < 2:
                continue
            sm = np.zeros(g.nb); sm[m] = np.exp(np.polyval(p[::-1], xx))
            chi2 = float(np.sum((yy - sm[m]) ** 2 / vv)); ndf = int(m.sum() - (o + 1))
            K = C.h_lin(R0, sm)
            out["%s_%g_%g_pol%d" % (lab, lo, hi, o)] = dict(chi2=chi2, ndf=ndf, p=float(1 - chi2dist.cdf(chi2, ndf)), K=K, K_rel=K / hs)
            print("%-7s %+.3f..%+.3f exp(pol%d): chi2 %.1f/%d (p %.3f), step projection %.2f rows (%.4f of h_sim)" % (lab, lo, hi, o, chi2, ndf, 1 - chi2dist.cdf(chi2, ndf), K, K / hs))
json.dump(out, open(C.OUT + "edgeu_bgsmooth.json", "w"), indent=1)
