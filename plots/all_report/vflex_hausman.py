#!/usr/bin/env python3
# vflex_hausman.py -- independent verification of package flex: is the estimator dependence of the F8 scales larger than
# a correct model with the real template noise would give? (Hausman-type test)
#   For the data and for RF pseudo-data (truth = F8 fit of the data; bootstrap replica templates as in the bias test) the
#   F8 scales are estimated twice: flex quasi-likelihood (BB-lite variance in the weights) and Poisson ML with the templates
#   taken as exact. The data difference is compared with the toy distribution of the difference.
#   Also: the Gaussian Barlow-Beeston profile likelihood with all scales bounded to [0, 10] (physical region).
# Run: python3 -B plots/all_report/vflex_hausman.py [ntoy]
import sys, json
import numpy as np
from scipy import optimize
from multiprocessing import Pool
sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import vflex_core as V

ntoy = int(sys.argv[1]) if len(sys.argv) > 1 else 300
D = V.load(); L = V.Layout(D, (1,))
M = V.FlexModel(["C2", "C3", "C4", "C5", "Geta", "GD_om", "GD_nom", "GC"]); T, Vt = V.templates(D, L, M)
TR = [V.templates(D, L, M, rep=r)[0] for r in range(V.NREPS)]


def pois_fit(y, cb, Tm, kc, x0):
    def f(t):
        mu = np.maximum(t[-1] * cb + Tm @ t[:8], 1e-9)
        return 2 * (mu - y * np.log(mu)).sum() + ((t[-1] - kc) / 0.05) ** 2
    def g(t):
        mu = np.maximum(t[-1] * cb + Tm @ t[:8], 1e-9)
        w = 2 * (1 - y / mu)
        return np.append(Tm.T @ w, cb @ w + 2 * (t[-1] - kc) / 0.05 ** 2)
    r = optimize.minimize(f, x0, jac=g, method="L-BFGS-B", options=dict(maxiter=50000, ftol=1e-16, gtol=1e-9))
    return r.x


def bb_fit(y, cb, Tm, kc, x0, bounds):
    def f(t):
        s, k = t[:8], t[-1]
        m0 = k * cb + Tm @ s; v = (Vt * s ** 2).sum(1); d = m0 - v
        mu = np.maximum(np.where(v > 0, 0.5 * (d + np.sqrt(d * d + 4 * v * y)), m0), 1e-9)
        return 2 * (mu - y * np.log(mu)).sum() + np.where(v > 0, v * (1 - y / mu) ** 2, 0).sum() + ((k - kc) / 0.05) ** 2
    r = optimize.minimize(f, x0, method="L-BFGS-B", bounds=bounds, options=dict(maxiter=50000, ftol=1e-16, gtol=1e-9))
    return r.x, r.fun


Rd = V.fit(M, L, T, Vt, egcov=L.eg)
thq = Rd['th']
thp = pois_fit(L.y, L.cb, T, 1.0, thq)
bnds = [(0, 10)] * 8 + [(0.5, 1.5)]
thb, fb = bb_fit(L.y, L.cb, T, 1.0, np.clip(thq, 0.01, 9), bnds)
muT = M.eval(thq, L, T, Vt)[0]


def one(t):
    rng = np.random.Generator(np.random.PCG64(4242 + t))
    y = rng.poisson(muT).astype(float)
    kc = thq[-1] + 0.05 * rng.standard_normal()
    cb = L.cbrep[:, t % 200]; Tm = TR[t % 100]
    R = V.fit(M, L, Tm, Vt, y=y, cb=cb, start=thq, cons={M.np - 1: (kc, 0.05)}, egcov=None, sandwich=False)
    p = pois_fit(y, cb, Tm, kc, thq)
    return R['th'], p, R['ok']


if __name__ == "__main__":
    with Pool(6) as pool:
        res = pool.map(one, range(ntoy), chunksize=10)
    Q = np.array([r[0] for r in res]); Pp = np.array([r[1] for r in res]); ok = np.array([r[2] for r in res])
    dif = (Pp - Q)[ok]
    dd = thp - thq
    print(f"data: Poisson-ML minus quasi; RF toys ({ok.sum()}): mean and rms of the same difference; z = (data - toy mean)/toy rms")
    out = {"names": M.names, "quasi": thq.tolist(), "poisson": thp.tolist(), "bb_bounded": thb.tolist(), "diff_data": dd.tolist(),
           "diff_toy_mean": dif.mean(0).tolist(), "diff_toy_rms": dif.std(0, ddof=1).tolist(), "toy_poisson_mean": Pp[ok].mean(0).tolist()}
    C = np.cov(dif.T, ddof=1)
    for k, n in enumerate(M.names):
        z = (dd[k] - dif[:, k].mean()) / dif[:, k].std(ddof=1)
        print(f"{n:9s} quasi {thq[k]:7.4f}  PoissonML {thp[k]:7.4f}  BB-profile[0,10] {thb[k]:7.4f} | data diff {dd[k]:+7.4f}  toy {dif[:, k].mean():+7.4f} +- {dif[:, k].std(ddof=1):6.4f}  z {z:+6.2f} | toy PoissonML mean {Pp[ok, k].mean():7.4f} (truth {thq[k]:.4f})")
    e = dd - dif.mean(0)
    chi2 = e @ np.linalg.solve(C, e)
    print(f"joint chi2 of the data difference (9 scales, toy covariance): {chi2:.1f}; toys exceeding it: "
          f"{np.mean([(x - dif.mean(0)) @ np.linalg.solve(C, x - dif.mean(0)) >= chi2 for x in dif]):.3f}")
    out["hausman_chi2"] = float(chi2)
    json.dump(out, open(V.DIR + "/vflex_hausman.json", "w"), indent=1)
