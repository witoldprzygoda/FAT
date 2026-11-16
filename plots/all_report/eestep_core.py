#!/usr/bin/env python3
# ============================================================================
# eestep_core.py  --  package eestep: loaders of the eestep_dump.C / eestep_simdump.C tables and the step estimator
# ============================================================================
# Written 2026-11-15 (prefix eestep). The step estimator is the one of edgeu (edgeu_core.fit_step, pol2 + h T(u) on
# -45 <= u < +45 MeV/c2, 2.5 MeV/c2 bins, T = smeared step with the simulated core width of edgeu_ana.json, variance =
# model (Poisson), the edgeu like-sign null test). u = M4 - E(M(pipi), M(ee)) (frozen edge function, data scale).
# ============================================================================
import json, math, sys
import numpy as np
sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import edgeu_core as C

OUT = "/hdd3/przygoda/HADES/FAT_out/all_report/eestep/"
G = C.GRID_NOM
SIG_U = json.load(open(C.OUT + "edgeu_ana.json"))["reso"]["sigma_core"]
T = C.step_gauss(G, SIG_U)

def load(table, layout="eestep_layout.txt"):
    lines = open(OUT + layout).read().split("\n")
    i = 0; tabs = {}
    while i < len(lines):
        l = lines[i]
        if l and not l.startswith(" "):
            nm, n, nc = l.split(); n, nc = int(n), int(nc)
            cols = [lines[i + 1 + j].strip() for j in range(nc)]
            tabs[nm] = (n, cols); i += nc + 1
        else:
            i += 1
    n, cols = tabs[table]
    D = {c: np.fromfile(OUT + "%s_%s.f32" % (table, c), dtype=np.float32) for c in cols}
    for c in cols:
        assert D[c].size == n, (table, c, D[c].size, n)
    D["n"] = n
    return D

def u_of(D):
    e = C.edge(D["mpipi"].astype(np.float64), D["mee"].astype(np.float64))
    u = D["m4"].astype(np.float64) - e
    u[~(D["mee"] < C.MPI0)] = np.nan
    return u

def in_slices(mpipi, lo=0.30, hi=0.40):
    pb = C.ppbin_of(np.asarray(mpipi, dtype=np.float32))
    return (pb >= C.PB(lo)) & (pb < C.PB(hi))

def hist_u(u, mask, weights=None):
    b = G.bin(u)
    ok = mask & (b >= 0)
    return np.bincount(b[ok], weights=None if weights is None else weights[ok], minlength=G.nb).astype(np.float64)

def step(y, lo=-0.045, hi=0.045, order=2):
    """edgeu like-sign null-test estimator on a histogram y (fine grid)"""
    R = C.fit_step(C.Spec(y, np.zeros(G.nb), np.zeros((G.nb, G.nb))), T, G, lo, hi, order)
    return R

def step_rows(u, mask, weights=None):
    y = hist_u(u, mask, weights)
    R = step(y)
    n = float(y[R.idx].sum())
    near = float(y[(G.x > -0.0151) & (G.x < 0.0)].sum() / 6.0)   # rows per bin in -15..0 MeV/c2
    return dict(h=R.h, eh=R.eh, chi2=R.chi2, ndf=R.ndf, rows=n, near=near, share=R.h / near if near > 0 else np.nan, R=R, y=y)

def fmt(r):
    return "h = %6.1f +- %5.1f  rows %6.0f  rows/bin(-15..0) %6.1f  share %5.3f +- %5.3f  chi2 %.1f/%d" % (
        r["h"], r["eh"], r["rows"], r["near"], r["share"], r["eh"] / r["near"] if r["near"] > 0 else np.nan, r["chi2"], r["ndf"])

def period_of_day(day):
    return C.period_of_day(np.asarray(day).astype(np.int16))
