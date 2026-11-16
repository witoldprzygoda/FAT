#!/usr/bin/env python3
# ============================================================================
# veestep_core.py  --  INDEPENDENT verification of eestep (prefix veestep): shared layer
# ============================================================================
# Written 2026-11-15 by the independent verifier. Own code: no eestep_* or edgeu_* module is imported. Uses only
#   - the row columns of veestep_dump.C (/hdd3/przygoda/HADES/FAT_out/all_report/veestep/d_*),
#   - the edgeu simulation column files (sim_*, simls_*; read with own loaders),
#   - the frozen definitions of defs.h (edge function, scales, slices), re-implemented here.
# Estimator (own choice, not the edgeu GLS): spectra in 2.5 MeV bins of u on -45..+45 MeV; model pol2(u) + h S(u),
# S(u) = 0.5 erfc(u / (sqrt2 sigma)), sigma = 6.7 MeV; h = step height in rows per 2.5 MeV bin below the edge.
#   like-sign alone: Poisson maximum likelihood (exact for small counts);
#   N - CB and simulation: weighted least squares with variances fixed from a smoothed first pass (linear estimator);
#   errors: 300 Poisson(1) event-group bootstrap replicas of the data (groups = adjacent rows of one event, the compact
#   grp column, per sample and day), for every quantity that mixes samples.
# CB (own construction of the D5 rule): per M(pipi) slice (0.30-0.34-0.37-0.39-0.40) and per M(ee) band, 10 MeV cells of u
# aligned at u = 0: CB(cell) = 2 sqrt(N++ N--); inside a cell the CB is shared over the 2.5 MeV bins in proportion to
# (N++ + N--) of the bin (keeps the like-sign shape inside the cell).
# ============================================================================
import math

import numpy as np
from scipy.optimize import minimize
from scipy.special import erfc

OUT = "/hdd3/przygoda/HADES/FAT_out/all_report/veestep/"
EDG = "/hdd3/przygoda/HADES/FAT_out/all_report/edgeu/"
MEta, MPi0 = 0.547862, 0.1349768
ESC, EPP = 0.99274, 0.9970          # defs.h kEdgeScale, kEdgePiPi
S_M4, S_PP = 0.99615, 0.9952         # D8 transfer of the simulation
CF_E1, CF_E2 = 0.879, 0.8825         # C2 gamma* factors (D35 / WP8)
A_P = 0.418249                       # peak scale (WP6/WP7)
SIG = 0.0067
SLICES = [(0.30, 0.34), (0.34, 0.37), (0.37, 0.39), (0.39, 0.40)]
BW = 0.0025
LO, HI = -0.045, 0.045
EDGES = np.round(np.arange(-0.050, 0.050 + 1e-9, BW), 6)   # 40 bins, cells of 4 bins aligned at 0
NB = len(EDGES) - 1
XC = 0.5 * (EDGES[:-1] + EDGES[1:])
FIT = (XC > LO) & (XC < HI)          # 36 bins


def m4max(m, q):
    """kinematic limit on the true scale: sqrt(m_eta^2 - m_eta (E0 - p0)(1 - q^2/m_pi0^2))"""
    m = np.asarray(m, np.float64); q = np.asarray(q, np.float64)
    E0 = (MEta ** 2 + MPi0 ** 2 - m * m) / (2 * MEta)
    p2 = E0 * E0 - MPi0 ** 2
    with np.errstate(invalid="ignore"):
        a = MEta ** 2 - MEta * (E0 - np.sqrt(np.where(p2 >= 0, p2, np.nan))) * (1 - q * q / MPi0 ** 2)
        return np.where((p2 >= 0) & (a > 0), np.sqrt(a), np.nan)


def edge(m, q):
    """data-scale edge E(m, q) of defs.h"""
    return ESC * m4max(np.asarray(m, np.float64) / EPP, q)


def col(path, dt):
    return np.fromfile(path, dtype=dt)


def load_data():
    n = int(open(OUT + "veestep_layout.txt").read().split()[1])
    spec = dict(k="i8", bits="u8", day="i16", ent="i32", grp="i32", nrow="i32", nlep="i32", s1="i32", s2="i32", sr="i32",
                m4="f32", mpipi="f32", mee="f32", oa="f32", p1="f32", p2="f32", t1="f32", t2="f32", f1="f32", f2="f32",
                ppip="f32", ppim="f32", vz="f32", lpos="i16", lneg="i16", hpos="i16", hneg="i16", tmult="i16", nus="i16",
                ntwin="i16", nls="i16", nflip="i16", joinok="u8")
    dmap = {"i8": np.int8, "u8": np.uint8, "i16": np.int16, "i32": np.int32, "f32": np.float32}
    D = {}
    for c, t in spec.items():
        a = col(OUT + "d_%s.%s" % (c, t), dmap[t])
        assert a.size == n, (c, a.size, n)
        D[c] = a
    D["n"] = n
    with np.errstate(invalid="ignore"):
        e = edge(D["mpipi"], D["mee"])
        u = D["m4"].astype(np.float64) - e
        u[~(D["mee"] < MPi0)] = np.nan
    D["u"] = u
    D["E"] = e
    D["S0"] = (D["bits"] & 1) != 0
    D["gkey"] = (D["k"].astype(np.int64) << 40) | (D["day"].astype(np.int64) << 32) | (D["grp"].astype(np.int64) & 0xFFFFFFFF)
    return D


def load_sim_us():
    """edgeu sim e+e- columns (own loader)"""
    lay = {}
    for l in open(EDG + "dump_layout.txt"):
        if not l.startswith("#"):
            a, b = l.split(); lay[a] = int(b)
    n = lay["sim"]
    S = {}
    for c, t in [("tree", np.int8), ("grp", np.int32), ("g", np.int8), ("cls", np.int8), ("lcls", np.int8), ("pcls", np.int8),
                 ("bits", np.uint8), ("m4", np.float32), ("m4true", np.float32), ("mpipi", np.float32), ("mee", np.float32),
                 ("w", np.float64), ("f108", np.float32)]:
        a = col(EDG + "sim_%s.%s" % (c, {np.int8: "i8", np.int32: "i32", np.uint8: "u8", np.float32: "f32", np.float64: "f64"}[t]), t)
        assert a.size == n
        S[c] = a
    m4 = S["m4"].astype(np.float64); mt = S["m4true"].astype(np.float64)
    S["m4p"] = np.where(mt > 0, S_M4 * (mt + (m4 - mt)), S_M4 * m4)
    S["mpp"] = S_PP * S["mpipi"].astype(np.float64)
    with np.errstate(invalid="ignore"):
        S["u"] = S["m4p"] - edge(S["mpp"], S["mee"])
        S["u"][~(S["mee"] < MPi0)] = np.nan
    ok = ((S["bits"] & 1) != 0) & ((S["bits"] & 8) != 0)
    w = S["w"] * S["f108"]
    c2 = S["g"] == 0
    w = np.where(c2 & (S["mee"] < 0.010), w * CF_E1, w)
    w = np.where(c2 & ~(S["mee"] < 0.010), w * CF_E2, w)
    S["wt"] = np.where(ok, w, 0.0)
    return S


def load_sim_ls():
    n = int(open(EDG + "dumpls_layout.txt").read().split()[1])
    S = {}
    for c, t, ext in [("tree", np.int8, "i8"), ("grp", np.int32, "i32"), ("tmpl", np.int8, "i8"), ("cls", np.int8, "i8"),
                      ("lcls", np.int8, "i8"), ("pcls", np.int8, "i8"), ("flags", np.int32, "i32"), ("ghost", np.int32, "i32"),
                      ("bits", np.uint8, "u8"), ("id1", np.int32, "i32"), ("id2", np.int32, "i32"), ("par1", np.int32, "i32"),
                      ("par2", np.int32, "i32"), ("gpar1", np.int32, "i32"), ("gpar2", np.int32, "i32"), ("proc1", np.int32, "i32"),
                      ("proc2", np.int32, "i32"), ("m4", np.float32, "f32"), ("m4true", np.float32, "f32"),
                      ("mpipi", np.float32, "f32"), ("mpipitrue", np.float32, "f32"), ("mee", np.float32, "f32"),
                      ("meetrue", np.float32, "f32"), ("w", np.float64, "f64"), ("f108", np.float32, "f32")]:
        a = col(EDG + "simls_%s.%s" % (c, ext), t)
        assert a.size == n, c
        S[c] = a
    m4 = S["m4"].astype(np.float64); mt = S["m4true"].astype(np.float64)
    S["m4p"] = np.where(mt > 0, S_M4 * (mt + (m4 - mt)), S_M4 * m4)
    S["mpp"] = S_PP * S["mpipi"].astype(np.float64)
    with np.errstate(invalid="ignore"):
        S["u"] = S["m4p"] - edge(S["mpp"], S["mee"])
        S["u"][~(S["mee"] < MPi0)] = np.nan
        S["utrue"] = mt - m4max(S["mpipitrue"], S["meetrue"])
    return S


def hist(u, w=None):
    h, _ = np.histogram(u, bins=EDGES, weights=w)
    return h.astype(np.float64)


def step_t(sigma=SIG, u0=0.0):
    """bin average of 0.5 erfc((u - u0)/(sqrt2 sigma)) over the 2.5 MeV bins"""
    def Gi(x):
        return 0.5 * (x * erfc(x / (math.sqrt(2) * sigma)) - math.sqrt(2 / math.pi) * sigma * np.exp(-x * x / (2 * sigma * sigma)))
    a = EDGES[:-1] - u0; b = EDGES[1:] - u0
    return (Gi(b) - Gi(a)) / BW


def design(sigma=SIG, u0=0.0, order=2):
    t = XC[FIT] / 0.045
    cols = [t ** j for j in range(order + 1)] + [step_t(sigma, u0)[FIT]]
    return np.vstack(cols).T


def fit_pois(y, sigma=SIG, u0=0.0, order=2, nostep=False):
    """Poisson ML fit of counts y (fine grid) with pol + h S; returns h, error, -2lnL, fitted model"""
    F = design(sigma, u0, order)
    if nostep:
        F = F[:, :-1]
    yy = y[FIT]
    p0 = np.linalg.lstsq(F, yy, rcond=None)[0]

    def nll(p):
        mu = F @ p
        if np.any(mu <= 0):
            return 1e30
        return float(np.sum(mu - yy * np.log(mu)))

    def grad(p):
        mu = F @ p
        return F.T @ (1 - yy / np.maximum(mu, 1e-12))
    r = minimize(nll, p0, jac=grad, method="BFGS", options=dict(gtol=1e-8, maxiter=2000))
    p = r.x
    mu = F @ p
    H = F.T @ (F * (yy / np.maximum(mu, 1e-12) ** 2)[:, None])
    C = np.linalg.pinv(H)
    m2 = 2 * float(np.sum(mu - yy + np.where(yy > 0, yy * np.log(np.where(yy > 0, yy, 1) / mu), 0)))
    if nostep:
        return dict(h=0.0, eh=0.0, chi2=m2, p=p, mu=mu)
    return dict(h=float(p[-1]), eh=float(math.sqrt(C[-1, -1])), chi2=m2, p=p, mu=mu, ndf=int(FIT.sum() - F.shape[1]))


class Lin:
    """linear least-squares step estimator with fixed variances var (fine grid, FIT bins)"""
    def __init__(self, var, sigma=SIG, u0=0.0, order=2):
        F = design(sigma, u0, order)
        Wt = 1.0 / np.maximum(var[FIT], 1.0)
        A = F.T @ (F * Wt[:, None])
        self.C = np.linalg.inv(A)
        self.L = self.C @ (F * Wt[:, None]).T
        self.F = F
        self.Wt = Wt

    def h(self, y):
        return float(self.L[-1] @ y[FIT])

    def p(self, y):
        return self.L @ y[FIT]

    def chi2(self, y):
        p = self.p(y)
        r = y[FIT] - self.F @ p
        return float(np.sum(r * r * self.Wt))


def cb_spec(pp, mm):
    """CB on the fine grid from like-sign fine histograms of ONE slice: 2 sqrt(N++ N--) per 10 MeV cell, shared in proportion
    to (N++ + N--) of the bins"""
    nc = NB // 4
    P = pp.reshape(nc, 4); M = mm.reshape(nc, 4)
    cb = 2 * np.sqrt(P.sum(1) * M.sum(1))
    s = P + M
    tot = s.sum(1)
    share = np.where(tot[:, None] > 0, s / np.where(tot[:, None] > 0, tot[:, None], 1), 0.25)
    return (cb[:, None] * share).reshape(NB)


def slice_masks(D, base):
    m = D["mpipi"]
    return [base & (m >= a) & (m < b) for (a, b) in SLICES]


def spectra(D, base, w=None):
    """N, N++, N--, CB summed over the slices (CB per slice); w: row weights or None"""
    N = np.zeros(NB); PP = np.zeros(NB); MM = np.zeros(NB); CB = np.zeros(NB)
    for sm in slice_masks(D, base):
        ww = None if w is None else w[sm]
        k = D["k"][sm]; u = D["u"][sm]
        wk = [None if ww is None else ww[k == j] for j in range(3)]
        n = hist(u[k == 0], wk[0]); p = hist(u[k == 1], wk[1]); q = hist(u[k == 2], wk[2])
        N += n; PP += p; MM += q; CB += cb_spec(p, q)
    return N, PP, MM, CB


def boot_weights(D, rows, nrep=300, seed=0x7EE57E9):
    key = D["gkey"][rows]
    uq, inv = np.unique(key, return_inverse=True)
    rng = np.random.default_rng(seed)
    K = rng.poisson(1.0, size=(nrep, uq.size)).astype(np.float32)
    return K[:, inv]


def sim_hist(S, mask, w=None):
    return hist(S["u"][mask], (S["wt"] if w is None else w)[mask])
