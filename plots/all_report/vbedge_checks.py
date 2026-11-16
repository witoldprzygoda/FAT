"""vbedge_checks.py -- INDEPENDENT VERIFIER of the package bedge: tests of the plain-words conclusions.

Written 2026-11-15 (prefix vbedge). Own code (vbedge_core.py, arrays of vbedge_fill.C).
  A. per-slice data/simulation pattern (rule R1 of the Dalitz-plot variants): chi2 with data errors only (as bedge)
     and with the template (MC) errors of the per-slice simulated steps added (own replicas)
  B. control slices 0.40-0.41 and 0.42-0.46 (fixed edge 0.515, half range 45 MeV)
  C. E7: step projection of a smooth non-eta background, own smoothing (scaled-Poisson ML exp(pol n) on 0.40-0.70,
     n = 2, 3 (bedge nominal), 4, 5; each order started from the previous solution); SMASH-like (kappa_B = 1) and stage-A-like (c_j (1 + b x)(G_D + 2.93 G_C))
  D. the noise of the kappa = 0.48 and C2-corrected readings and the E1/E2 difference
Output: stdout and /hdd3/przygoda/HADES/FAT_out/all_report/vbedge/vbedge_checks.json
Run: python3 -I plots/all_report/vbedge_checks.py
"""
import json
import sys

import numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2 as chi2d

sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
from vbedge_core import *  # noqa

d = load()
N, S, R = d["N"], d["S"], d["R"]
OUTJ = {}


def simy(arr, cls, pa, pbb):
    return arr[:, pa:pbb].sum(axis=(0, 1)) if cls < 0 else arr[cls, pa:pbb].sum(0)


# ------------------------------------------------------------------------------------------- A. pattern test
print("A. per-slice h_data / h_sim (all M(ee)); chi2 against a constant")
cfgnames = {0: "KLOE-2 (nominal)", 2: "phase space (100)", 3: "weights off (0)"}
hd, ehd = [], []
for lo, hi in NOMSET:
    y, cb, cov, u = data_spec(N, -1, pb(lo), pb(hi))
    h, eh, *_ = fit(y, cb, cov, Cfg(lo, hi)); hd.append(h); ehd.append(eh)
hd, ehd = np.array(hd), np.array(ehd)
OUTJ["pattern"] = {}
for cf, nm in cfgnames.items():
    hs = np.array([fit(simy(S[cf, 2], -1, pb(lo), pb(hi)), None, None, Cfg(lo, hi))[0] for lo, hi in NOMSET])
    hr = np.array([[fit(simy(R[r, cf, 2], -1, pb(lo), pb(hi)), None, None, Cfg(lo, hi))[0] for lo, hi in NOMSET] for r in range(NRS)])
    mc = hr.std(0, ddof=1)
    rat = hd / hs
    e_d = ehd / hs
    e_t = np.sqrt(e_d ** 2 + (rat * mc / hs) ** 2)
    def c2(e):
        w = 1 / e ** 2; m = (rat * w).sum() / w.sum(); return ((rat - m) ** 2 * w).sum(), m
    cd, md = c2(e_d); ct, mt = c2(e_t)
    print(f"   {nm:18s} ratios {' '.join('%.3f' % x for x in rat)}; data err {' '.join('%.3f' % x for x in e_d)}; MC rel {' '.join('%.3f' % x for x in mc / hs)}")
    print(f"   {'':18s} chi2 data-only {cd:6.1f}/3 (p {chi2d.sf(cd, 3):.1e});  with MC {ct:6.1f}/3 (p {chi2d.sf(ct, 3):.1e})")
    OUTJ["pattern"][nm] = dict(rat=rat.tolist(), e_d=e_d.tolist(), mc_rel=(mc / hs).tolist(), chi2_d=cd, chi2_t=ct)
# the same with the MC noise of the DIFFERENCE pattern (variant vs nominal, same events): does the data prefer KLOE?
# direct likelihood comparison: chi2 of data vs a f x h_sim (free f), errors data (+) MC of that variant
print()

# ------------------------------------------------------------------------------------------- B. control slices
print("B. control slices (not in E)")
for lo, hi, fe in ((0.40, 0.41, None), (0.42, 0.46, 0.515)):
    for cls in (-1, 0, 1):
        c = Cfg(lo, hi, hr=0.045 if fe else None, fixed_edge=fe)
        y, cb, cov, u = data_spec(N, cls, pb(lo), pb(hi))
        h, eh, *_ = fit(y, cb, cov, c)
        hs = fit(simy(S[0, 2], cls, pb(lo), pb(hi)), None, None, c)[0]
        print(f"   {lo:.2f}-{hi:.2f} cls {cls:2d} edge {c.edge:.4f}: data {h:7.1f} +- {eh:5.1f}, sim {hs:7.1f}")
print()

# ------------------------------------------------------------------------------------------- C. E7
print("C. E7: step projection of a smooth non-eta background (data covariance of each slice)")
xs = X
sel = np.where((X > 0.40) & (X < 0.70))[0]
GLx, GLw = np.polynomial.legendre.leggauss(8)


def smooth(c, order, p0=None):
    """scaled-Poisson ML fit of exp(pol(order)) in x = (M4 - 0.55)/0.1 to bin contents c on the 0.40-0.70 bins;
    returns the bin integrals on the full grid (0 outside 0.40-0.70)."""
    cc = c[sel]
    edges_lo = M4LO + sel * M4W
    xq = np.array([[(e + 0.5 * M4W + 0.5 * M4W * g - 0.55) / 0.1 for g in GLx] for e in edges_lo])
    wq = 0.5 * M4W * GLw
    def model(p):
        P = np.zeros_like(xq)
        for k in range(order + 1):
            P += p[k] * xq ** k
        return (np.exp(P) * wq).sum(1)
    def nll(p):
        m = np.maximum(model(p), 1e-300)
        return float((m - cc * np.log(m)).sum())
    if p0 is None:
        p0 = np.zeros(order + 1); p0[0] = np.log(cc.sum() / (0.30))
    best = None
    for meth in ("Nelder-Mead", "BFGS", "Powell", "BFGS"):
        r = minimize(nll, p0 if best is None else best.x, method=meth, options=dict(maxiter=50000))
        if best is None or r.fun < best.fun:
            best = r
    out = np.zeros(NM4); out[sel] = model(best.x)
    return out, best.x


res = {}
for k in range(3):
    cls = k - 1
    tot = {}
    for lo, hi in NOMSET:
        pa, pbb = pb(lo), pb(hi)
        c = Cfg(lo, hi)
        y, cb, cov, u = data_spec(N, cls, pa, pbb)
        V = fit(y, cb, cov, c)[4]
        sl_j = {0.30: 0.9816, 0.34: 0.9844, 0.37: 0.9505, 0.39: 1.0473}[lo]
        for gname, gi in (("GD", 4), ("GC", 5)):
            yr = simy(S[0, gi], cls, pa, pbb)
            px = None
            for order in (2, 3, 4, 5):
                ys, px = smooth(yr, order, None if px is None else np.append(px, 0.0))
                tilt = (X - 0.55) / 0.1
                h = fit(ys, None, None, c, Vfix=V)[0]
                hx = fit(ys * tilt, None, None, c, Vfix=V)[0]
                tot.setdefault((gname, order), [0, 0]); tot[(gname, order)][0] += h; tot[(gname, order)][1] += hx
                tot.setdefault((gname, order, "cj"), [0, 0]); tot[(gname, order, "cj")][0] += sl_j * h; tot[(gname, order, "cj")][1] += sl_j * hx
            hraw = fit(yr, None, None, c, Vfix=V)[0]
            tot.setdefault((gname, "raw"), 0); tot[(gname, "raw")] += hraw
        B = sum(fit(simy(S[0, gi], cls, pa, pbb), None, None, c, Vfix=V)[0] for gi in (0, 1))
        tot["B"] = tot.get("B", 0) + B
    res[k] = tot
    b, kB = -0.0639, 2.932
    for order in (2, 3, 4, 5):
        Ksim = tot[("GD", order)][0] + tot[("GC", order)][0]
        KA = (tot[("GD", order, "cj")][0] + b * tot[("GD", order, "cj")][1]) + kB * (tot[("GC", order, "cj")][0] + b * tot[("GC", order, "cj")][1])
        f = 0.39162
        print(f"   class {k} exp(pol{order}): G_D {tot[('GD', order)][0]:7.1f}  G_C {tot[('GC', order)][0]:7.1f}  K_sim {Ksim:7.1f}  K_A {KA:7.1f}"
              f"  B {tot['B']:.1f}  bias(SMASH-like data) {(Ksim - f*Ksim)/tot['B']:+.4f}  bias(stage-A data) {(KA - f*Ksim)/tot['B']:+.4f}")
    print(f"   class {k} raw: G_D {tot[('GD', 'raw')]:.1f} G_C {tot[('GC', 'raw')]:.1f}")
OUTJ["E7"] = {str(k): {str(kk): v for kk, v in t.items()} for k, t in res.items()}
print()

# ------------------------------------------------------------------------------------------- D. E1 vs E2
a = json.load(open(f"{OUT}/vbedge_ana.json"))
print("D. E1 vs E2 (own numbers): f_E1 %.4f, f_E2 %.4f" % (a["q0"]["f_E1"], a["q0"]["f_E2"]))
d12 = a["q0"]["f_E2"] - a["q0"]["f_E1"]
e = (a["stat_boot"]["f_E1"] ** 2 + a["stat_boot"]["f_E2"] ** 2 + a["mc"]["f_E1"] ** 2 + a["mc"]["f_E2"] ** 2) ** .5
print("   difference %.4f +- %.4f (stat (+) MC, classes treated as independent) = %.2f sigma" % (d12, e, d12 / e))
json.dump(OUTJ, open(f"{OUT}/vbedge_checks.json", "w"), indent=1, default=float)
