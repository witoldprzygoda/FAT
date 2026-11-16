#!/usr/bin/env python3
"""vpt_core.py -- INDEPENDENT VERIFIER of the package peaktail (the low-mass tail of the eta peak): own shared layer.

Written 2026-11-15 (prefix vpt). No peaktail_* code or output and no edgeu_core / bedge code is imported: the edge
function, u, the CB rule, the step template and the GLS fit are re-written here from the written definitions
(defs.h Edge, edgeu note Sect. 3, ScanEngine CB rule). Inputs: the columns of vpt_dump.C in
/hdd3/przygoda/HADES/FAT_out/all_report/vpt/.
"""
import math
import numpy as np
from scipy.special import erfc, erf

OUT = "/hdd3/przygoda/HADES/FAT_out/all_report/vpt/"
META, MPI0, ME_E, MPI = 0.547862, 0.1349768, 0.000510999, 0.13957
ESC, EPP = 0.99274, 0.9970
S_NOM, SPP = 0.99615, 0.9952
AP = 0.418249
CF_E1, CF_E2 = 0.879, 0.8825
SIG_STEP = 0.006701875966571207
PEAKM = 0.5439


def lay():
    d = {}
    for l in open(OUT + "layout.txt"):
        k, v = l.split(); d[k] = int(v)
    return d


def col(n, dt, N):
    a = np.fromfile(OUT + n, dtype=dt); assert a.size == N, (n, a.size, N); return a


def load_data():
    N = lay()["data"]
    D = {k: col("d_%s" % f, dt, N) for k, f, dt in [("k", "k.i8", np.int8), ("grp", "grp.i32", np.int32), ("m4", "m4.f32", np.float32),
                                                     ("mpipi", "mpipi.f32", np.float32), ("mee", "mee.f32", np.float32), ("bits", "bits.u8", np.uint8)]}
    return D


def load_sim():
    N = lay()["sim"]
    spec = [("tree", "tree.i8", np.int8), ("grp", "grp.i32", np.int32), ("g", "g.i8", np.int8), ("bits", "bits.u8", np.uint8),
            ("m4", "m4.f32", np.float32), ("m4t", "m4t.f32", np.float32), ("pp", "pp.f32", np.float32), ("ppt", "ppt.f32", np.float32),
            ("ee", "ee.f32", np.float32), ("eet", "eet.f32", np.float32), ("w", "w.f64", np.float64), ("f108", "f108.f32", np.float32)]
    S = {k: col("s_" + f, dt, N) for k, f, dt in spec}
    for n in ("ep", "em", "pip", "pim"):
        for q in ("p", "pt", "pr", "th", "tht", "ph"):
            S[n + "_" + q] = col("s_%s_%s.f32" % (n, q), np.float32, N)
    for n in ("ep", "em"):
        S[n + "_vr"] = col("s_%s_vr.f32" % n, np.float32, N); S[n + "_vz"] = col("s_%s_vz.f32" % n, np.float32, N)
        S[n + "_proc"] = col("s_%s_proc.i32" % n, np.int32, N); S[n + "_par"] = col("s_%s_par.i32" % n, np.int32, N)
    return S


def m4max(m, q):
    m = np.asarray(m, float); q = np.asarray(q, float)
    E0 = (META ** 2 + MPI0 ** 2 - m * m) / (2 * META)
    p2 = E0 * E0 - MPI0 ** 2
    with np.errstate(invalid="ignore"):
        a = META ** 2 - META * (E0 - np.sqrt(np.where(p2 >= 0, p2, np.nan))) * (1 - q * q / MPI0 ** 2)
        return np.where((p2 >= 0) & (a > 0), np.sqrt(a), np.nan)


def edge(m, q):
    return ESC * m4max(np.asarray(m, float) / EPP, q)


def sim_m4p(S, s=S_NOM, r=1.0, m4=None):
    m4 = S["m4"].astype(float) if m4 is None else m4
    mt = S["m4t"].astype(float)
    return np.where(mt > 0, s * (mt + r * (m4 - mt)), s * m4)


def sim_weight(S, sel=1, c2=True):
    w = S["w"] * S["f108"]
    ok = ((S["bits"] & sel) != 0) & ((S["bits"] & 8) != 0)
    w = np.where(ok, w, 0.0)
    if c2:
        isc2 = S["g"] == 0
        w = np.where(isc2 & (S["ee"] < 0.010), w * CF_E1, w)
        w = np.where(isc2 & ~(S["ee"] < 0.010), w * CF_E2, w)
    return w


# ------------------------------------------------------------------ four-vectors (GeV), for substitution studies
def p4(p, th, ph, m):
    p = p.astype(float) / 1000.0; th = np.radians(th.astype(float)); ph = np.radians(ph.astype(float))
    px, py, pz = p * np.sin(th) * np.cos(ph), p * np.sin(th) * np.sin(ph), p * np.cos(th)
    return np.sqrt(p * p + m * m), px, py, pz


def minv(*vs):
    E = sum(v[0] for v in vs); x = sum(v[1] for v in vs); y = sum(v[2] for v in vs); z = sum(v[3] for v in vs)
    return np.sqrt(np.maximum(E * E - x * x - y * y - z * z, 0))


# ------------------------------------------------------------------ grid, CB, step, fit
class Grid:
    def __init__(self, lo=-0.200, bw=0.0025, nb=140, reb=4):
        self.lo, self.bw, self.nb, self.reb = lo, bw, nb, reb
        self.e = lo + bw * np.arange(nb + 1); self.x = lo + bw * (np.arange(nb) + 0.5)

    def h(self, v, w=None):
        return np.histogram(v, bins=self.e, weights=w)[0].astype(float)


def cb_cells(pp, mm, g, mode=0):
    """2 sqrt(N++ N--) in cells of reb fine bins, linearly interpolated between cell centres (flat beyond the outer
    centres); returns CB per fine bin and its covariance (cell variance N++ + N--)"""
    nw = g.nb // g.reb
    P = pp[:nw * g.reb].reshape(nw, g.reb).sum(1); M = mm[:nw * g.reb].reshape(nw, g.reb).sum(1)
    cw = 2 * np.sqrt(P * M) if mode == 0 else P + M
    A = np.zeros((g.nb, nw))
    cc = g.lo + g.reb * g.bw * (np.arange(nw) + 0.5)
    for i, x in enumerate(g.x):
        t = (x - cc[0]) / (g.reb * g.bw)
        if t <= 0:
            A[i, 0] = 1.0 / g.reb
        elif t >= nw - 1:
            A[i, nw - 1] = 1.0 / g.reb
        else:
            j = int(t); f = t - j; A[i, j] += (1 - f) / g.reb; A[i, j + 1] += f / g.reb
    return A @ cw, A @ np.diag(P + M) @ A.T


def step_t(g, s, u0=0.0):
    def I(x):
        return 0.5 * (x * erfc(x / (math.sqrt(2) * s)) - math.sqrt(2 / math.pi) * s * np.exp(-x * x / (2 * s * s)))
    return (I(g.e[1:] - u0) - I(g.e[:-1] - u0)) / g.bw


def gls(y, F, V):
    Vi = np.linalg.inv(V); C = np.linalg.inv(F.T @ Vi @ F); L = C @ F.T @ Vi; p = L @ y; r = y - F @ p
    return p, C, L, float(r @ Vi @ r)


class Fit:
    pass


def fitstep(sig, cb, cbcov, T, g, lo=-0.045, hi=0.045, order=2, fixed=None, extra=None, Vfix=None, nit=6):
    """N - CB - fixed = pol(order) [+ extra columns] + h T; V = cov(CB) + diag(max(model of N, 1)) iterated"""
    i = np.where((g.x >= lo - 1e-12) & (g.x < hi - 1e-12))[0]
    t = (g.x[i] - 0.5 * (lo + hi)) / (0.5 * (hi - lo))
    cols = [t ** j for j in range(order + 1)] + ([c[i] for c in extra] if extra else []) + [T[i]]
    F = np.vstack(cols).T
    fx = np.zeros(len(i)) if fixed is None else fixed[i]
    y = sig[i] - fx
    V0 = cbcov[np.ix_(i, i)]
    if Vfix is not None:
        V = Vfix; p, C, L, chi2 = gls(y, F, V)
    else:
        V = V0 + np.diag(np.maximum(y + fx + cb[i], 1.0)); p, C, L, chi2 = gls(y, F, V)
        for _ in range(nit):
            h0 = p[-1]
            V = V0 + np.diag(np.maximum(F @ p + fx + cb[i], 1.0)); p, C, L, chi2 = gls(y, F, V)
            if abs(p[-1] - h0) < 1e-4 * math.sqrt(C[-1, -1]):
                break
    R = Fit(); R.i, R.p, R.C, R.L, R.chi2, R.V, R.F = i, p, C, L, chi2, V, F
    R.h, R.eh, R.ndf = p[-1], math.sqrt(C[-1, -1]), len(i) - F.shape[1]
    R.y = y; R.model = F @ p
    return R


def hlin(R, y):
    return float(R.L[-1] @ y[R.i])
