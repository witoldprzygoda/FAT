#!/usr/bin/env python3
# vflex_estim.py -- independent verification of package flex: estimator dependence of the nominal F8 scales
#   (a) the flex quasi-likelihood (BB-lite variance in the weights) -- reproduced in vflex_fit.py
#   (b) Poisson maximum likelihood with the templates taken as exact (no template-noise treatment)
#   (c) Barlow-Beeston profile likelihood with Gaussian template nuisances (one per raw template and bin), profiled
#       analytically (= additive BB-lite nuisance on the template sum, variance sum_raw s_t^2 var(T_t,i))
#   plus the per-window balance of data and model, and the same quantities for RF pseudo-data (is the overshoot of the
#   model total a property of the estimator with noisy templates?)
# Run: python3 -B plots/all_report/vflex_estim.py
import sys, json
import numpy as np
from scipy import optimize
sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import vflex_core as V

D = V.load(); L = V.Layout(D, (1,))
M = V.FlexModel(["C2", "C3", "C4", "C5", "Geta", "GD_om", "GD_nom", "GC"]); T, Vt = V.templates(D, L, M)
raw = ~M.bg
y = L.y


def bbl_nll(t, y=y, cb=L.cb, T=T, Vt=Vt, kc=1.0):
    """profile -2 lnL with one Gaussian nuisance per raw template and bin (T_t,i measures tau_t,i with variance var_t,i);
    profiled analytically: mu = [(m0 - v) + sqrt((m0 - v)^2 + 4 v y)] / 2, penalty v (1 - y/mu)^2, v = sum s_t^2 var_t,i
    (equivalent to the additive Barlow-Beeston-lite nuisance delta_i ~ N(0, v_i) on the template sum)"""
    s, k = t[:8], t[-1]
    m0 = k * cb + T @ s
    v = (Vt * s ** 2).sum(1)
    d = m0 - v
    mu = np.where(v > 0, 0.5 * (d + np.sqrt(d * d + 4 * v * y)), m0)
    mu = np.maximum(mu, 1e-9)
    pen = np.where(v > 0, v * (1 - y / mu) ** 2, 0.0)
    nll = 2 * (mu - y * np.log(mu)).sum() + pen.sum() + ((k - kc) / 0.05) ** 2
    return nll, mu, mu - m0


def pois_nll(t, y=y, cb=L.cb, T=T, kc=1.0):
    mu = np.maximum(t[-1] * cb + T @ t[:8], 1e-9)
    return 2 * (mu - y * np.log(mu)).sum() + ((t[-1] - kc) / 0.05) ** 2


def mini(f, x0):
    best = None
    for meth in ("L-BFGS-B", "Powell", "BFGS"):
        r = optimize.minimize(f, x0 if best is None else best.x, method=meth, options=dict(maxiter=100000, **({"ftol": 1e-15, "gtol": 1e-10} if meth == "L-BFGS-B" else {})))
        if best is None or r.fun < best.fun - 1e-9:
            best = r
    return best


R = V.fit(M, L, T, Vt, egcov=L.eg)
th_q = R['th']
rp = mini(pois_nll, th_q)
rb = mini(lambda t: bbl_nll(t)[0], th_q)
# BB-lite errors: inverse Hessian (numerical) of the profile -2 lnL / 2
def hess(f, x, h=1e-4):
    n = len(x); H = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            e_i = np.zeros(n); e_j = np.zeros(n); e_i[i] = h; e_j[j] = h
            H[i, j] = H[j, i] = (f(x + e_i + e_j) - f(x + e_i - e_j) - f(x - e_i + e_j) + f(x - e_i - e_j)) / (4 * h * h)
    return H
Hb = hess(lambda t: bbl_nll(t)[0] / 2, rb.x)
eb = np.sqrt(np.diag(np.linalg.inv(Hb)))
out = {"names": M.names, "quasi": th_q.tolist(), "poisson": rp.x.tolist(), "bblite": rb.x.tolist(), "bblite_err_hess": eb.tolist()}
print(f"{'scale':9s} {'quasi(flex)':>11s} {'Poisson ML':>11s} {'BB prof L':>10s} {'+- (Hess)':>9s}")
for i, n in enumerate(M.names):
    print(f"{n:9s} {th_q[i]:11.4f} {rp.x[i]:11.4f} {rb.x[i]:10.4f} {eb[i]:9.4f}")
mub = bbl_nll(rb.x)[1]
xc = 0.5 * (L.lo + L.hi)
for wn, (a, b) in dict(fit=(0.40, 0.70), low=(0.40, 0.455), B=(0.455, 0.5225), P=(0.5336, 0.5542), CRU=(0.565, 0.70)).items():
    s = (xc > a) & (xc < b)
    print(f"window {wn:4s}: N {y[s].sum():8.0f}  quasi model {R['mu'][s].sum():9.1f} ({100 * (R['mu'][s].sum() / y[s].sum() - 1):+.2f} %)  "
          f"BB-lite mu {mub[s].sum():9.1f} ({100 * (mub[s].sum() / y[s].sum() - 1):+.2f} %)")
out["win_quasi_over_data_fit"] = float(R['mu'].sum() / y.sum() - 1)

# RF pseudo-data: does the quasi-likelihood fit also overshoot the total? (truth = F8 fit of the data)
TR = [V.templates(D, L, M, rep=r)[0] for r in range(V.NREPS)]
muT = M.eval(th_q, L, T, Vt)[0]
over = []; overP = []; sP = (xc > 0.5336) & (xc < 0.5542)
for t in range(200):
    rng = np.random.Generator(np.random.PCG64(777 + t))
    yt = rng.poisson(muT).astype(float)
    cons = {M.np - 1: (th_q[-1] + 0.05 * rng.standard_normal(), 0.05)}
    Rt = V.fit(M, L, TR[t % 100], Vt, y=yt, cb=L.cbrep[:, t % 200], start=th_q, cons=cons, egcov=None, sandwich=False)
    over.append(Rt['mu'].sum() / yt.sum() - 1); overP.append(Rt['mu'][sP].sum() / yt[sP].sum() - 1)
over = np.array(over); overP = np.array(overP)
print(f"RF toys (200): model total / data total - 1 = {100 * over.mean():+.2f} +- {100 * over.std(ddof=1):.2f} % (data: {100 * out['win_quasi_over_data_fit']:+.2f} %); "
      f"peak window P: {100 * overP.mean():+.2f} +- {100 * overP.std(ddof=1):.2f} %")
out["rf_over_mean"] = float(over.mean()); out["rf_over_rms"] = float(over.std(ddof=1))
out["rf_overP_mean"] = float(overP.mean()); out["rf_overP_rms"] = float(overP.std(ddof=1))
json.dump(out, open(V.DIR + "/vflex_estim.json", "w"), indent=1)
