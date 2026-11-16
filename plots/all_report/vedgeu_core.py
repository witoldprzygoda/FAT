#!/usr/bin/env python3
# ============================================================================
# vedgeu_core.py  --  INDEPENDENT VERIFIER of the package edgeu: own columns, edge, u, CB, step fit, replicas
# ============================================================================
# Written 2026-11-15 (prefix vedgeu). Does NOT import edgeu_core.py or read any edgeu column file. Reads the columns
# of vedgeu_fill.C in /hdd3/przygoda/HADES/FAT_out/all_report/vedgeu/ (layout in its header).
# Everything below is written from the descriptions in report_all / defs.h / the edgeu note, not copied:
#   edge   E(m, q) = 0.99274 sqrt(m_eta^2 - m_eta (E_pi0 - p_pi0)(1 - q^2/m_pi0^2)), m' = m / 0.9970
#   u      data: M4 - E(M(pipi), M(ee)); sim: 0.99615 [m4_t + r (m4 - m4_t)] - E(0.9952 M(pipi), M(ee))
#   CB     2 sqrt(N++ N--) in 10 MeV cells of u, linear interpolation between cell centres (flat outside), per slice
#   fit    GLS of N - CB - a_P K on pol2(t) + h T(u), V = cov(CB) + diag(max(model of N, 1)), iterated
#   boot   Poisson(1) per (tree, event group), numpy PCG64 with own seeds (data: 200, sim: 100)
# ============================================================================
import math

import numpy as np
from scipy.special import erfc

DIR = "/hdd3/przygoda/HADES/FAT_out/all_report/vedgeu/"
MET, MP0 = 0.547862, 0.1349768
ESC, EPP = 0.99274, 0.9970
S_N, R_N, SPP_N = 0.99615, 1.0, 0.9952
CF1, CF2 = 0.879, 0.8825
AP = 0.418249


def m4max(m, q):
    """kinematic upper limit of M(pi+pi-e+e-) in eta -> pi+pi-pi0, pi0 -> gamma e+e-(q), on the true scale"""
    m = np.asarray(m, float); q = np.asarray(q, float)
    E0 = (MET ** 2 + MP0 ** 2 - m ** 2) / (2 * MET)
    p2 = E0 ** 2 - MP0 ** 2
    with np.errstate(invalid="ignore"):
        p = np.sqrt(np.where(p2 >= 0, p2, np.nan))
        egam = (E0 - p) * (MP0 ** 2 - q ** 2) / (2 * MP0 ** 2)     # softest photon energy in the eta frame
        a = MET ** 2 - 2 * MET * egam
        return np.where(a > 0, np.sqrt(np.where(a > 0, a, np.nan)), np.nan)


def E(m, q, esc=ESC, epp=EPP):
    return esc * m4max(np.asarray(m, float) / epp, q)


def load():
    def c(n, dt):
        return np.fromfile(DIR + n, dtype=dt)
    D = dict(k=c("d_k.i1", np.int8), grp=c("d_grp.i4", np.int32), fid=c("d_fid.i2", np.int16), m4=c("d_m4.f4", np.float32),
             pp=c("d_pp.f4", np.float32), ee=c("d_ee.f4", np.float32), bits=c("d_bits.u1", np.uint8))
    S = dict(tree=c("s_tree.i1", np.int8), grp=c("s_grp.i4", np.int32), c=c("s_c.i1", np.int8), om=c("s_om.i1", np.int8),
             bits=c("s_bits.u1", np.uint8), m4=c("s_m4.f4", np.float32), m4t=c("s_m4t.f4", np.float32), pp=c("s_pp.f4", np.float32),
             ppt=c("s_ppt.f4", np.float32), ee=c("s_ee.f4", np.float32), eet=c("s_eet.f4", np.float32), w=c("s_w.f8", np.float64),
             w0=c("s_w0.f8", np.float64))
    return D, S


def u_data(D, esc=ESC, epp=EPP, useq=True):
    q = D["ee"].astype(float) if useq else np.zeros(len(D["ee"]))
    u = D["m4"].astype(float) - E(D["pp"].astype(float), q, esc, epp)
    u[~(D["ee"] < MP0)] = np.nan
    return u


def m4p_sim(S, s=S_N, r=R_N):
    m4 = S["m4"].astype(float); mt = S["m4t"].astype(float)
    return np.where(mt > 0, s * (mt + r * (m4 - mt)), s * m4)


def u_sim(S, s=S_N, r=R_N, spp=SPP_N, useq=True):
    mpp = spp * S["pp"].astype(float)
    q = S["ee"].astype(float) if useq else np.zeros(len(S["ee"]))
    u = m4p_sim(S, s, r) - E(mpp, q)
    u[~(S["ee"] < MP0)] = np.nan
    return u, mpp


# ---------------------------------------------------------------- u grid and histograms
LO, BW, NB = -0.200, 0.0025, 140
X = LO + BW * (np.arange(NB) + 0.5)
EDG = LO + BW * np.arange(NB + 1)
SLICES = [(0.30, 0.34), (0.34, 0.37), (0.37, 0.39), (0.39, 0.40)]


def ubin(u, lo=LO, bw=BW, nb=NB):
    b = np.full(len(u), -1, dtype=np.int64)
    ok = np.isfinite(u)
    bb = np.floor((u[ok] - lo) / bw).astype(np.int64)
    bb[(bb < 0) | (bb >= nb)] = -1
    b[ok] = bb
    return b


def slice_index(mpp, slices=SLICES):
    s = np.full(len(mpp), -1, dtype=np.int64)
    m = np.asarray(mpp, dtype=np.float32)
    for i, (a, b) in enumerate(slices):
        s[(m >= np.float32(a)) & (m < np.float32(b))] = i
    return s


def hist(b, s, w=None, ns=4, nb=NB):
    ok = (b >= 0) & (s >= 0)
    flat = s[ok] * nb + b[ok]
    return np.bincount(flat, weights=None if w is None else w[ok], minlength=ns * nb).reshape(ns, nb)


def hist_reps(b, s, W, ns=4, nb=NB):
    """W: (nrep, nrow) weights -> (nrep, ns, nb)"""
    ok = (b >= 0) & (s >= 0)
    flat = s[ok] * nb + b[ok]
    out = np.empty((W.shape[0], ns * nb))
    for r in range(W.shape[0]):
        out[r] = np.bincount(flat, weights=W[r, ok], minlength=ns * nb)
    return out.reshape(W.shape[0], ns, nb)


# ---------------------------------------------------------------- CB
def cb_one(pp, mm, reb=4, lo=LO, bw=BW, mode="geo", k=1.0):
    nb = len(pp); nc = nb // reb
    P = pp[:nc * reb].reshape(nc, reb).sum(1); M = mm[:nc * reb].reshape(nc, reb).sum(1)
    v = k * (2 * np.sqrt(P * M) if mode == "geo" else P + M)
    var = k * k * (P + M)
    cc = lo + reb * bw * (np.arange(nc) + 0.5)
    x = lo + bw * (np.arange(nb) + 0.5)
    A = np.zeros((nb, nc))
    for i in range(nb):
        if x[i] <= cc[0]:
            A[i, 0] = 1
        elif x[i] >= cc[-1]:
            A[i, -1] = 1
        else:
            j = int(np.searchsorted(cc, x[i]) - 1)
            f = (x[i] - cc[j]) / (cc[j + 1] - cc[j])
            A[i, j] = 1 - f; A[i, j + 1] = f
    A /= reb
    return A @ v, A @ np.diag(var) @ A.T


def spectrum(Hus, Hpp, Hmm, mode="geo", reb=4, k=1.0):
    """sum over slices of N, CB (per slice), cov(CB)"""
    N = Hus.sum(0)
    CB = np.zeros(Hus.shape[1]); COV = np.zeros((Hus.shape[1], Hus.shape[1]))
    for i in range(Hus.shape[0]):
        cb, cv = cb_one(Hpp[i], Hmm[i], reb=reb, mode=mode, k=k)
        CB += cb; COV += cv
    return N, CB, COV


# ---------------------------------------------------------------- step template and fit
def step(sig, edges=EDG, u0=0.0):
    """bin average of 0.5 erfc((u - u0) / (sqrt(2) sig)) by Simpson integration (own numerics, 41 points per bin)"""
    out = np.empty(len(edges) - 1)
    for i in range(len(edges) - 1):
        xx = np.linspace(edges[i], edges[i + 1], 41) - u0
        yy = 0.5 * erfc(xx / (math.sqrt(2) * sig))
        out[i] = (yy[0] + yy[-1] + 4 * yy[1:-1:2].sum() + 2 * yy[2:-1:2].sum()) / (3 * 40)
    return out


class Fit:
    pass


def gls_fit(N, CB, CBCOV, T, fixed=None, lo=-0.045, hi=0.045, order=2, x=X, Vfix=None, niter=8, extra=None):
    idx = np.where((x > lo - 1e-9) & (x < hi - 1e-9))[0]
    t = x[idx] / (0.5 * (hi - lo)) - (0.5 * (hi + lo)) / (0.5 * (hi - lo))
    cols = [t ** j for j in range(order + 1)]
    if extra is not None:
        cols += [e[idx] for e in extra]
    cols.append(T[idx])
    F = np.array(cols).T
    fx = np.zeros(len(idx)) if fixed is None else fixed[idx]
    y = (N - CB)[idx] - fx
    V0 = CBCOV[np.ix_(idx, idx)]
    def solve(V):
        Vi = np.linalg.inv(V)
        Cm = np.linalg.inv(F.T @ Vi @ F)
        Lm = Cm @ F.T @ Vi
        p = Lm @ y
        r = y - F @ p
        return p, Cm, Lm, float(r @ Vi @ r), Vi
    if Vfix is not None:
        V = Vfix
        p, Cm, Lm, chi2, Vi = solve(V)
    else:
        V = V0 + np.diag(np.maximum(N[idx], 1.0))
        p, Cm, Lm, chi2, Vi = solve(V)
        for _ in range(niter):
            V = V0 + np.diag(np.maximum(F @ p + fx + CB[idx], 1.0))
            p, Cm, Lm, chi2, Vi = solve(V)
    F0 = F[:, :-1]
    Ci = np.linalg.inv(F0.T @ Vi @ F0); p0 = Ci @ F0.T @ Vi @ y; r0 = y - F0 @ p0
    R = Fit()
    R.idx, R.p, R.C, R.L, R.V, R.F, R.y = idx, p, Cm, Lm, V, F, y
    R.h, R.eh, R.chi2, R.chi20, R.ndf = float(p[-1]), float(math.sqrt(Cm[-1, -1])), chi2, float(r0 @ Vi @ r0), len(idx) - F.shape[1]
    R.Lh = Lm[-1]
    R.smooth = F[:, :order + 1] @ p[:order + 1]
    return R


def proj(R, y):
    return float(R.Lh @ np.asarray(y)[R.idx])


# ---------------------------------------------------------------- Poisson(1) replica weights per (tree, group)
def group_weights(tree, grp, n, seed):
    key = tree.astype(np.int64) * (1 << 32) + grp.astype(np.int64)
    uq, inv = np.unique(key, return_inverse=True)
    rng = np.random.Generator(np.random.PCG64(seed))
    k = rng.poisson(1.0, size=(n, uq.size)).astype(np.uint8)
    return k[:, inv]


def sd(a):
    a = np.asarray(a, float)
    return float(a.std(ddof=1))
