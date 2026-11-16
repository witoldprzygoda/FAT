"""vbedge_checks2.py -- INDEPENDENT VERIFIER of bedge: folds, periods, raw-rest MC share, smoother convergence.

Written 2026-11-15 (prefix vbedge). Own code; arrays of vbedge_fill.C (dSub.bin: fold0 fold1 perI perIIA perIIB perIII).
Run: python3 -I plots/all_report/vbedge_checks2.py
"""
import sys
import numpy as np
from scipy.stats import chi2 as chi2d
sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
from vbedge_core import *  # noqa

d = load()
N, S, R = d["N"], d["S"], d["R"]
Sub = np.fromfile(f"{OUT}/dSub.bin").reshape(6, NS, NCL, NPP, NM4)


def Edata(Nx, cls=-1, sset=NOMSET):
    E = v = 0
    for lo, hi in sset:
        y, cb, cov, u = data_spec(Nx, cls, pb(lo), pb(hi))
        h, eh, *_ = fit(y, cb, cov, Cfg(lo, hi)); E += h; v += eh * eh
    return E, v ** .5


def Esim(arr, cls=-1, sset=NOMSET):
    return sum(fit((arr[:, pb(lo):pb(hi)].sum(axis=(0, 1)) if cls < 0 else arr[cls, pb(lo):pb(hi)].sum(0)), None, None, Cfg(lo, hi))[0] for lo, hi in sset)


Es0 = Esim(S[0, 2])
tot_us = N[0].sum()
print("folds (simulation x row fraction of the us rows):")
fv = []
for k in range(2):
    fr = Sub[k, 0].sum() / tot_us
    E, e = Edata(Sub[k])
    fv.append((E / (Es0 * fr), e / (Es0 * fr)))
    print("   fold %d: row fraction %.4f, f = %.4f +- %.4f" % (k, fr, *fv[-1]))
z = (fv[1][0] - fv[0][0]) / (fv[0][1] ** 2 + fv[1][1] ** 2) ** .5
print("   z = %.2f (fit errors), %.2f (x 1.036)" % (z, z / 1.036))
Lp = [2809.453, 372.610, 1078.022, 1819.110]
Cs = [1.0241, 1.0, 1.0, 1.0]
for usec in (0, 1):
    xs = []
    for p in range(4):
        E, e = Edata(Sub[2 + p])
        sc = Es0 * Lp[p] / 6079.2
        k = Cs[p] if usec else 1.0
        xs.append((E / sc * k, e / sc * k * 1.036))
    w = np.array([1 / e ** 2 for x, e in xs]); x = np.array([x for x, e in xs])
    m = (w * x).sum() / w.sum(); ch = (w * (x - m) ** 2).sum()
    print("periods %s: %s; chi2 %.2f/3, p = %.3f" % ("x C_sec" if usec else "raw    ", " ".join("%.3f(%.0f)" % (a, 1000 * b) for a, b in xs), ch, chi2d.sf(ch, 3)))

# raw rest (T12 - C2 - C3) share of the MC variance, own replicas (data covariance of the nominal fits)
Vd = []
for lo, hi in NOMSET:
    y, cb, cov, u = data_spec(N, -1, pb(lo), pb(hi)); Vd.append(fit(y, cb, cov, Cfg(lo, hi))[4])
def lin(arr):
    return sum(fit(arr[:, pb(lo):pb(hi)].sum(axis=(0, 1)), None, None, Cfg(lo, hi), Vfix=V)[0] for (lo, hi), V in zip(NOMSET, Vd))
rest = np.array([lin(R[r, 0, 2]) - lin(R[r, 0, 0]) - lin(R[r, 0, 1]) for r in range(NRS)])
es = np.array([Esim(R[r, 0, 2]) for r in range(NRS)])
print("raw rest of the simulated step: nominal %.1f; replica rms %.1f; E_sim replica rms %.1f -> share of variance %.2f" % (
    lin(S[0, 2]) - lin(S[0, 0]) - lin(S[0, 1]), rest.std(ddof=1), es.std(ddof=1), (rest.std(ddof=1) / es.std(ddof=1)) ** 2))
