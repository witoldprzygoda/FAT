#!/usr/bin/env python3
# ============================================================================
# edgeu_core.py  --  package edgeu: shared layer (columns, edge function, u, spectra, CB rule, step fit, replicas)
# ============================================================================
# Written 2026-11-15 (prefix edgeu). The eta bump (eta -> pi+pi-pi0, pi0 -> e+e-gamma or gamma gamma with one photon
# converted, the other photon lost) measured in the distance to its kinematic edge,
#     u = M4 - E(M(pipi), M(ee)),
# E = the frozen edge function of defs.h (alldefs::Edge, Sect. 2 of report_all: E(m, q) = 0.99274 sqrt(m_eta^2 - m_eta
# (E_pi0 - p_pi0)(1 - q^2/m_pi0^2)), m' = m/0.9970). Every bump event has u <= 0 (up to the resolution), whatever its
# M(pipi) and M(ee): u puts the edges of all events on top of each other.
#
# Inputs: the column files of edgeu_dump.C in /hdd3/przygoda/HADES/FAT_out/all_report/edgeu/ (layout in dump_layout.txt).
# Conventions (Part III): data rows S0 (nominal) or pre, cor_ok == 1, isBest >= 0; N - CB in RAW rows per u bin, CB =
# 2 sqrt(N++ N--) in 10 MeV cells of u, interpolated to the fine bins (the ScanEngine rule, here applied to the variable
# u as scan_bump.C did for its aligned spectra), per M(pipi) slice; simulation: expected RAW data rows (w x RowFactor x
# 27.0718), mode A owned rows, V3 (code 108), M4' = 0.99615 [m4_t + r (m4 - m4_t)], M(pipi)' = 0.9952 M(pipi), M(ee) not
# scaled; C2 x c(E-class) for the HGeant gamma* angle (WP8 means, E1 0.879, E2 0.8825; D35).
# ============================================================================
import math
import os

import numpy as np
from scipy.special import erfc

OUT = "/hdd3/przygoda/HADES/FAT_out/all_report/edgeu/"
FIG = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_all/"
TAB = "/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all/"
NOTE = "/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/work_all/edgeu_results.md"

# ---------------------------------------------------------------- constants (defs.h, ScanEngine.h, Part III inputs)
ME, MPI0 = 0.547862, 0.1349768          # scan::kMeta, scan::kMpi0
EDGE_SCALE, EDGE_PIPI = 0.99274, 0.9970 # alldefs::kEdgeScale (eta peak 0.5439/m_eta), kEdgePiPi (K0S in the data)
S_NOM, R_NOM, SPP_NOM = 0.99615, 1.0, 0.9952   # D8 mass transfer of the simulation
S_ERR = 0.0009
L_DATA = 6079.2
PT3 = 64.0 / 63.0
SIM_TO_RAW = L_DATA / 221.05 / PT3      # 27.0718 (already applied in the dump)
AP, AP_STAT, AP_SYST, AP_MC = 0.418249, 0.035, 0.067, 0.025   # peak scale Y_peak / SMASH (C4 + C5) in P (WP6/WP7)
CF_E1, CF_E2, CF_ALL = 0.879, 0.8825, 0.8855                  # C2 gamma* angle factors (WP8 toy/closure means; bedge "rec")
CF_ERR = 0.024
GROUPS = ["C2", "C3", "C4", "C5", "Geta", "GDom", "GDnom", "GCom", "GCnom"]
G = {n: i for i, n in enumerate(GROUPS)}

# ---------------------------------------------------------------- grids
PP_LO, PP_W, NPP = 0.27, 0.005, 38      # M(pipi) 5 MeV bins 0.27-0.46 (data scale)
def ppbin_of(m):
    """index of the 5 MeV M(pipi) bin; edges rounded to 1e-3 as bedge (PPEdge)"""
    e = np.round((PP_LO + PP_W * np.arange(NPP + 1)) * 1000.) / 1000.
    b = np.searchsorted(e.astype(np.float32), m.astype(np.float32), side="right") - 1
    b[(b < 0) | (b >= NPP)] = -1
    return b
def PB(x):
    return int(round((x - PP_LO) / PP_W))
# M(ee) classes: 0 = < 10 MeV (E1), 1 = 10-20 MeV, 2 = >= 20 MeV (E2 = 1 + 2); coarse G3/F5g edges only (blinding D19)
NCLS = 3
def cls_of(mee):
    c = np.full(mee.shape, 2, dtype=np.int8)
    c[mee < 0.020] = 1
    c[mee < 0.010] = 0
    return c
CLS_SETS = {"all": [0, 1, 2], "E1": [0], "E2": [1, 2], "E20": [2]}

class Grid:
    """u grid: nb bins of width bw from lo; reb fine bins per CB cell (10 MeV)"""
    def __init__(self, lo=-0.200, bw=0.0025, nb=140, reb=4):
        self.lo, self.bw, self.nb, self.reb = lo, bw, nb, reb
        self.edges = lo + bw * np.arange(nb + 1)
        self.x = lo + bw * (np.arange(nb) + 0.5)
    def bin(self, u):
        uu = np.where(np.isfinite(u), u, self.lo - 1.0)
        b = np.floor((uu - self.lo) / self.bw).astype(np.int64)
        b[(b < 0) | (b >= self.nb)] = -1
        return b
GRID_NOM = Grid()

# ---------------------------------------------------------------- edge function (alldefs::Edge, vectorised)
def m4max_true(m, q):
    """kinematic limit on the TRUE scale (no scale factors), m = M(pipi), q = M(ee); NaN above the limit"""
    m = np.asarray(m, dtype=np.float64); q = np.asarray(q, dtype=np.float64)
    E0 = (ME * ME + MPI0 * MPI0 - m * m) / (2. * ME)
    p2 = E0 * E0 - MPI0 * MPI0
    with np.errstate(invalid="ignore"):
        a = ME * ME - ME * (E0 - np.sqrt(np.where(p2 >= 0, p2, np.nan))) * (1. - q * q / (MPI0 * MPI0))
        return np.where((p2 >= 0) & (a > 0), np.sqrt(a), np.nan)
def edge(m, q, scale=EDGE_SCALE, spp=EDGE_PIPI, useq=True):
    """E(m, q) on the data scale (defs.h Edge): scale x M4max(m/spp, q)"""
    qq = q if useq else np.zeros_like(np.asarray(q, dtype=np.float64))
    return scale * m4max_true(np.asarray(m, dtype=np.float64) / spp, qq)
def m4min_true(m):
    E0 = (ME * ME + MPI0 * MPI0 - m * m) / (2. * ME)
    p2 = E0 * E0 - MPI0 * MPI0
    with np.errstate(invalid="ignore"):
        a = ME * ME - ME * (E0 + np.sqrt(np.where(p2 >= 0, p2, np.nan)))
        return np.where((p2 >= 0) & (a > 0), np.sqrt(a), np.nan)

# ---------------------------------------------------------------- columns
def _col(name, dt, n):
    a = np.fromfile(OUT + name, dtype=dt)
    assert a.size == n, (name, a.size, n)
    return a
def layout():
    n = {}
    for l in open(OUT + "dump_layout.txt"):
        if l.startswith("#"):
            continue
        k, v = l.split()
        n[k] = int(v)
    return n
def load_data():
    n = layout()["data"]
    D = dict(k=_col("data_k.i8", np.int8, n), grp=_col("data_grp.i32", np.int32, n), m4=_col("data_m4.f32", np.float32, n),
             mpipi=_col("data_mpipi.f32", np.float32, n), mee=_col("data_mee.f32", np.float32, n),
             bits=_col("data_bits.u8", np.uint8, n), day=_col("data_day.i16", np.int16, n), fold=_col("data_fold.i8", np.int8, n))
    D["n"] = n
    return D
def load_sim():
    n = layout()["sim"]
    S = dict(tree=_col("sim_tree.i8", np.int8, n), grp=_col("sim_grp.i32", np.int32, n), g=_col("sim_g.i8", np.int8, n),
             cls=_col("sim_cls.i8", np.int8, n), lcls=_col("sim_lcls.i8", np.int8, n), pcls=_col("sim_pcls.i8", np.int8, n),
             bits=_col("sim_bits.u8", np.uint8, n), m4=_col("sim_m4.f32", np.float32, n), m4t=_col("sim_m4true.f32", np.float32, n),
             mpipi=_col("sim_mpipi.f32", np.float32, n), mpipit=_col("sim_mpipitrue.f32", np.float32, n),
             mee=_col("sim_mee.f32", np.float32, n), meet=_col("sim_meetrue.f32", np.float32, n), w=_col("sim_w.f64", np.float64, n))
    for c in (108, 116, 100, 107, 76):
        S["f%d" % c] = _col("sim_f%d.f32" % c, np.float32, n)
    S["n"] = n
    return S

# ---------------------------------------------------------------- periods (defs.h B)
def period_of_day(day):
    p = np.full(day.shape, -1, dtype=np.int8)
    p[(day >= 34) & (day <= 46)] = 0
    p[(day >= 47) & (day <= 48)] = 1
    p[(day >= 49) & (day <= 53)] = 2
    p[(day >= 59) & (day <= 67)] = 3
    return p
LP = [2809.453, 372.610, 1078.022, 1819.110]   # analysed L per period (bedge_ana.C, WP2/WP7)
CSEC_P = [1.0241, 1.0, 1.0, 1.0]

# ---------------------------------------------------------------- u of data and simulated rows
def u_data(D, scale=EDGE_SCALE, spp=EDGE_PIPI, useq=True, qcut=True, mshift=0.0):
    """u of every data row (NaN where the edge is not defined or M(ee) >= m_pi0 with qcut); mshift: pseudo edge
    E(M(pipi) + mshift, q) for the control region without a bump"""
    e = edge(D["mpipi"].astype(np.float64) + mshift, D["mee"].astype(np.float64), scale, spp, useq)
    u = D["m4"].astype(np.float64) - e
    if qcut:
        u[~(D["mee"] < MPI0)] = np.nan
    return u
def sim_m4p(S, s=S_NOM, r=R_NOM):
    m4 = S["m4"].astype(np.float64); mt = S["m4t"].astype(np.float64)
    return np.where(mt > 0, s * (mt + r * (m4 - mt)), s * m4)
CPP_SIM = 1.0018   # simulated M(pipi) scale reco/true (K0S, D8); replaced by the measured C2+C3 median in edgeu_ana.py
def sim_mppp(S, spp=SPP_NOM, rpp=None, cpp=None):
    """M(pipi)' = spp M(pipi); with rpp: spp [c m_t + rpp (M(pipi) - c m_t)] (resolution transfer of M(pipi), c = reco/true
    scale of the simulation)"""
    mp = S["mpipi"].astype(np.float64)
    if rpp is None:
        return spp * mp
    c = CPP_SIM if cpp is None else cpp
    mt = S["mpipit"].astype(np.float64)
    return np.where(mt > 0, spp * (c * mt + rpp * (mp - c * mt)), spp * mp)
def u_sim(S, s=S_NOM, r=R_NOM, spp=SPP_NOM, rpp=None, scale=EDGE_SCALE, epp=EDGE_PIPI, useq=True, qcut=True, cpp=None,
          mshift=0.0):
    m4p = sim_m4p(S, s, r); mpp = sim_mppp(S, spp, rpp, cpp)
    u = m4p - edge(mpp + mshift, S["mee"].astype(np.float64), scale, epp, useq)
    if qcut:
        u[~(S["mee"] < MPI0)] = np.nan
    return u, mpp, m4p

# ---------------------------------------------------------------- sim weights
def sim_weight(S, code=108, sel="S0", mode="A", best=False, c2corr=True):
    w = S["w"] * (S["f%d" % code] if code else 1.0)
    b = S["bits"]
    m = (b & (1 if sel == "S0" else 2)) != 0
    m &= (b & (8 if mode == "A" else 16)) != 0
    if best:
        m &= (b & 4) != 0
    w = np.where(m, w, 0.0)
    if c2corr:
        c2 = S["g"] == G["C2"]
        w = np.where(c2 & (S["mee"] < 0.010), w * CF_E1, w)
        w = np.where(c2 & ~(S["mee"] < 0.010), w * CF_E2, w)
    return w

# ---------------------------------------------------------------- histograms [.., pp bin, cls, u bin]
def hist_data(D, u, selmask, grid=GRID_NOM, weights=None):
    """H[k 3][pp NPP][cls NCLS][nb]; weights: None or (nrep x nrows) array -> H[rep][k][pp][cls][nb]"""
    b = grid.bin(u); pb = ppbin_of(D["mpipi"]); c = cls_of(D["mee"])
    ok = selmask & (b >= 0) & (pb >= 0)
    flat = ((D["k"][ok].astype(np.int64) * NPP + pb[ok]) * NCLS + c[ok]) * grid.nb + b[ok]
    size = 3 * NPP * NCLS * grid.nb
    if weights is None:
        return np.bincount(flat, minlength=size).astype(np.float64).reshape(3, NPP, NCLS, grid.nb)
    W = weights[:, ok]
    out = np.empty((W.shape[0], size))
    for r in range(W.shape[0]):
        out[r] = np.bincount(flat, weights=W[r], minlength=size)
    return out.reshape(W.shape[0], 3, NPP, NCLS, grid.nb)
def hist_sim(S, u, mpp, w, grid=GRID_NOM, groups=None, reps=None):
    """H[g 9][pp][cls][nb] of expected raw rows; reps: (rep weights array nrep x nrows) -> H[rep][g][pp][cls][nb]"""
    b = grid.bin(u); pb = ppbin_of(mpp.astype(np.float32)); c = cls_of(S["mee"])
    ok = (b >= 0) & (pb >= 0) & (w != 0)
    if groups is not None:
        ok &= np.isin(S["g"], groups)
    flat = ((S["g"][ok].astype(np.int64) * NPP + pb[ok]) * NCLS + c[ok]) * grid.nb + b[ok]
    size = 9 * NPP * NCLS * grid.nb
    if reps is None:
        return np.bincount(flat, weights=w[ok], minlength=size).reshape(9, NPP, NCLS, grid.nb)
    out = np.empty((reps.shape[0], size))
    wk = w[ok]
    R = reps[:, ok]
    for r in range(R.shape[0]):
        out[r] = np.bincount(flat, weights=wk * R[r], minlength=size)
    return out.reshape(reps.shape[0], 9, NPP, NCLS, grid.nb)
def hist_sim_var(S, u, mpp, w, grid=GRID_NOM):
    """event-group variance of the sim histograms (sum over groups of (sum w in the bin)^2)"""
    b = grid.bin(u); pb = ppbin_of(mpp.astype(np.float32)); c = cls_of(S["mee"])
    ok = (b >= 0) & (pb >= 0) & (w != 0)
    flat = ((S["g"][ok].astype(np.int64) * NPP + pb[ok]) * NCLS + c[ok]) * grid.nb + b[ok]
    key = S["tree"][ok].astype(np.int64) * (1 << 32) + S["grp"][ok].astype(np.int64)
    kk = key * (9 * NPP * NCLS * grid.nb) + flat
    uq, inv = np.unique(kk, return_inverse=True)
    s = np.bincount(inv, weights=w[ok])
    fl = uq % (9 * NPP * NCLS * grid.nb)
    return np.bincount(fl, weights=s * s, minlength=9 * NPP * NCLS * grid.nb).reshape(9, NPP, NCLS, grid.nb)

# ---------------------------------------------------------------- CB rule (ScanEngine::Spectra, applied to u)
def interp_matrix(grid):
    nw = grid.nb // grid.reb
    A = np.zeros((grid.nb, nw))
    for i in range(grid.nb):
        t = (grid.x[i] - (grid.lo + 0.5 * grid.reb * grid.bw)) / (grid.reb * grid.bw)
        if t <= 0:
            A[i, 0] += 1.0 / grid.reb
        elif t >= nw - 1:
            A[i, nw - 1] += 1.0 / grid.reb
        else:
            j0 = int(t); f = t - j0
            A[i, j0] += (1 - f) / grid.reb; A[i, j0 + 1] += f / grid.reb
    return A
_AM = {}
def cb_rule(pp, mm, grid, k=1.0, mode=0):
    """CB per fine bin and its covariance from like-sign counts (one slice): 2 sqrt(N++ N--) in cells (mode 0) or
    N++ + N-- (mode 1), x k; var(cell) = k^2 (N++ + N--) (ScanEngine)"""
    key = (grid.lo, grid.bw, grid.nb, grid.reb)
    if key not in _AM:
        _AM[key] = interp_matrix(grid)
    A = _AM[key]
    nw = grid.nb // grid.reb
    wp = pp[:nw * grid.reb].reshape(nw, grid.reb).sum(1); wm = mm[:nw * grid.reb].reshape(nw, grid.reb).sum(1)
    cbW = k * (2.0 * np.sqrt(wp * wm) if mode == 0 else (wp + wm))
    cbV = k * k * (wp + wm)
    return A @ cbW, (A * cbV) @ A.T

class Spec:
    """N, CB, cov(CB), N - CB on the fine grid"""
    def __init__(self, N, cb, cbcov, pp=None, mm=None):
        self.N, self.cb, self.cbcov, self.pp, self.mm = N, cb, cbcov, pp, mm
        self.sig = N - cb

def slices_of(pp_lo, pp_hi, edges=None):
    """slices of 5 MeV units for the CB: the bedge slices 0.30-0.34-0.37-0.39-0.40 inside [pp_lo, pp_hi), extended by
    one slice per 5 MeV unit outside them"""
    if edges is None:
        edges = [0.30, 0.34, 0.37, 0.39, 0.40]
    cuts = sorted(set([pp_lo, pp_hi] + [e for e in edges if pp_lo < e < pp_hi]))
    out = []
    for a, b in zip(cuts[:-1], cuts[1:]):
        if a < edges[0] or b > edges[-1]:
            x = a
            while x < b - 1e-9:
                out.append((PB(x), PB(min(b, x + PP_W)))); x += PP_W
        else:
            out.append((PB(a), PB(b)))
    return out

def data_spec(H, slices, cls="all", grid=GRID_NOM, cbk=1.0, cbmode=0, clsrule=0, share_rng=(-0.100, 0.100), cbslice=True):
    """N - CB summed over slices (CB per slice), classes 'all', 'E1', 'E2', 'E20'.
    H: [k][pp][cls][nb]. clsrule 0: CB(all classes) x like-sign class share over share_rng (per slice); 1: own rule."""
    cs = CLS_SETS[cls]
    N = np.zeros(grid.nb); CB = np.zeros(grid.nb); COV = np.zeros((grid.nb, grid.nb))
    PPs = np.zeros(grid.nb); MMs = np.zeros(grid.nb)
    sl = slices if cbslice else [(min(a for a, b in slices), max(b for a, b in slices))]
    xm = (grid.x >= share_rng[0]) & (grid.x < share_rng[1])
    for (a, b) in sl:
        n = H[0, a:b][:, cs].sum((0, 1)); p_all = H[1, a:b].sum((0, 1)); m_all = H[2, a:b].sum((0, 1))
        p_c = H[1, a:b][:, cs].sum((0, 1)); m_c = H[2, a:b][:, cs].sum((0, 1))
        if clsrule == 0 or cls == "all":
            if cbmode in (2, 3):
                # shape of one like-sign sample (2: e+e+, 3: e-e-) in cells, normalised to the nominal CB over share_rng
                cb0, cov0 = cb_rule(p_all, m_all, grid, cbk, 0)
                one = p_all if cbmode == 2 else m_all
                cbs, covs = cb_rule(one, one, grid, 1.0, 1)   # cells of N (mode 1 gives 2 N per cell)
                nrm = cb0[xm].sum() / cbs[xm].sum() if cbs[xm].sum() > 0 else 0.
                cb, cov = cbs * nrm, covs * nrm * nrm * 2.0   # var of 2 N per cell = 4 N = 2 x (cb_rule var 2 N)
            else:
                cb, cov = cb_rule(p_all, m_all, grid, cbk, cbmode)
            den = (p_all + m_all)[xm].sum()
            sh = (p_c + m_c)[xm].sum() / den if den > 0 else 0.
            cb, cov = cb * sh, cov * sh * sh
        else:
            cb, cov = cb_rule(p_c, m_c, grid, cbk, cbmode)
        N += n; CB += cb; COV += cov; PPs += p_c; MMs += m_c
    return Spec(N, CB, COV, PPs, MMs)

def sim_sum(Hs, groups, slices, cls="all"):
    cs = CLS_SETS[cls]
    y = 0.
    for g in groups:
        for (a, b) in slices:
            y = y + Hs[G[g], a:b][:, cs].sum((0, 1))
    return np.asarray(y, dtype=np.float64)

# ---------------------------------------------------------------- step templates
def step_gauss(grid, sigma, u0=0.0):
    """bin average of 0.5 erfc((u - u0)/(sqrt2 sigma))"""
    def Gi(x):
        return 0.5 * (x * erfc(x / (math.sqrt(2) * sigma)) - math.sqrt(2 / math.pi) * sigma * np.exp(-x * x / (2 * sigma * sigma)))
    a = grid.edges[:-1] - u0; b = grid.edges[1:] - u0
    return (Gi(b) - Gi(a)) / grid.bw
def step_emp(grid, delta, wd, u0=0.0):
    """bin average of P(delta > u - u0) for the weighted resolution sample (delta, wd)"""
    W = wd.sum()
    a = grid.edges[:-1] - u0
    out = np.zeros(grid.nb)
    for i in range(grid.nb):
        out[i] = (wd * np.clip(delta - a[i], 0, grid.bw)).sum() / (W * grid.bw)
    return out
def step_box(grid, sigma):
    """sharp step convolved with a box... (not used)"""
    return step_gauss(grid, sigma)

# ---------------------------------------------------------------- the fit (copy of the logic of scan::FitEdgeStep)
class FitRes:
    pass
def gls(y, F, V):
    Vi = np.linalg.inv(V)
    A = F.T @ Vi @ F
    C = np.linalg.inv(A)
    Lm = C @ F.T @ Vi
    p = Lm @ y
    r = y - F @ p
    return p, C, Lm, float(r @ Vi @ r)
def fit_step(spec, T, grid=GRID_NOM, lo=-0.045, hi=0.045, order=2, fixed=None, free_cols=None, Vfix=None, niter=6,
             own_var=None):
    """N - CB - fixed = sum_j p_j t^j [+ free_cols] + h T. fixed: fixed model part (e.g. a_P x peak template) on the
    fine grid; free_cols: list of extra template columns with free amplitudes. V = cov(CB) + diag(max(model of N, 1)),
    iterated (6 times, as scan::FitEdgeStep); Vfix: fixed covariance. own_var: for a simulated spectrum (no CB),
    var = max(model, 1) iterated (WP6 SimSpectra fit)."""
    idx = np.where((grid.x >= lo - 1e-12) & (grid.x < hi - 1e-12))[0]
    x = grid.x[idx]
    x0, dx = 0.5 * (lo + hi), 0.5 * (hi - lo)
    t = (x - x0) / dx
    cols = [t ** j for j in range(order + 1)]
    if free_cols is not None:
        cols += [c[idx] for c in free_cols]
    cols.append(T[idx])
    F = np.vstack(cols).T
    fx = np.zeros(len(idx)) if fixed is None else fixed[idx]
    if own_var is not None:
        y = own_var[idx] - fx
        ad = np.zeros(len(idx)); V0 = np.zeros((len(idx), len(idx)))
    else:
        y = spec.sig[idx] - fx
        ad = spec.cb[idx]; V0 = spec.cbcov[np.ix_(idx, idx)]
    if Vfix is not None:
        V = Vfix
        p, C, Lm, chi2 = gls(y, F, V)
    else:
        var = np.maximum(y + fx + ad, 1.0)
        V = V0 + np.diag(var)
        p, C, Lm, chi2 = gls(y, F, V)
        for it in range(niter):
            hold = p[-1]
            var = np.maximum(F @ p + fx + ad, 1.0)
            V = V0 + np.diag(var)
            p, C, Lm, chi2 = gls(y, F, V)
            if abs(p[-1] - hold) < 1e-4 * math.sqrt(C[-1, -1]):
                break
    F0 = F[:, :-1]
    p0, C0, L0, chi20 = gls(y, F0, V)
    R = FitRes()
    R.idx, R.x, R.F, R.V, R.p, R.C, R.L, R.y, R.fx = idx, x, F, V, p, C, Lm, y, fx
    R.h, R.eh = float(p[-1]), float(math.sqrt(C[-1, -1]))
    R.chi2, R.chi20, R.ndf = chi2, chi20, len(idx) - F.shape[1]
    R.Lh = Lm[-1]          # h = Lh . y (linear in the spectrum for the fixed V)
    R.model = F @ p; R.smooth = F[:, :order + 1] @ p[:order + 1]
    R.lo, R.hi, R.order = lo, hi, order
    return R
def h_lin(R, y_full):
    """h of another spectrum (fine grid) with the same linear estimator (fixed V of fit R)"""
    return float(R.Lh @ y_full[R.idx])

# ---------------------------------------------------------------- replica generators
_MASK = (1 << 64) - 1
def _mix(x):
    with np.errstate(over="ignore"):
        x = x + np.uint64(0x9E3779B97F4A7C15)
        x = (x ^ (x >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
        x = (x ^ (x >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
        return x ^ (x >> np.uint64(31))
_P1CDF = None
def _pcdf():
    global _P1CDF
    if _P1CDF is None:
        p = 0.36787944117144233; c = p; cdf = [c]
        for k in range(1, 31):
            p /= k; c += p; cdf.append(c)
        _P1CDF = np.array(cdf[:30])
    return _P1CDF
def group_pois(seed, tree, grp, n):
    """bedge::GroupPois (WP0 counter-based Poisson(1) replica weights), vectorised over groups: returns (n x ngroups) uint8"""
    tree = np.asarray(tree).astype(np.uint64); grp = np.asarray(grp).astype(np.int64).astype(np.uint32).astype(np.uint64)
    with np.errstate(over="ignore"):
        base = _mix(np.uint64(seed) ^ _mix((tree << np.uint64(32)) | grp))
        out = np.empty((n, base.size), dtype=np.uint8)
        cdf = _pcdf()
        for r in range(n):
            h = _mix(base + np.uint64(0x632BE59BD9B4E019) * np.uint64(r + 1))
            uu = (h >> np.uint64(11)).astype(np.float64) * (1.0 / 9007199254740992.0)
            out[r] = np.minimum(np.searchsorted(cdf, uu, side="right"), 30)
    return out
SEED_SIM = 0x5157A11
def sim_replica_weights(S, rows, n=100, seed=SEED_SIM):
    """(n x len(rows)) Poisson(1) weights of the WP0 generator for the selected sim rows (per event group)"""
    key = S["tree"][rows].astype(np.int64) * (1 << 32) + S["grp"][rows].astype(np.int64)
    uq, inv = np.unique(key, return_inverse=True)
    k = group_pois(seed, (uq >> 32), (uq & 0xFFFFFFFF), n)
    return k[:, inv]
SEED_DATA = 0xED6E0B00
def data_boot_weights(D, rows, n=200, seed=SEED_DATA):
    """(n x len(rows)) Poisson(1) weights per event group of each tree (us, pp, mm independent, as WP0/bedge), own seed"""
    key = D["k"][rows].astype(np.int64) * (1 << 32) + D["grp"][rows].astype(np.int64)
    uq, inv = np.unique(key, return_inverse=True)
    rng = np.random.default_rng(seed)
    k = rng.poisson(1.0, size=(n, uq.size)).astype(np.uint8)
    return k[:, inv]

# ---------------------------------------------------------------- small statistics helpers
def sd(a):
    a = np.asarray(a, dtype=np.float64)
    return float(np.sqrt(np.mean((a - a.mean()) ** 2) * len(a) / max(len(a) - 1, 1)))
def h68(a):
    a = np.sort(np.asarray(a, dtype=np.float64))
    return float(0.5 * (np.percentile(a, 84.135) - np.percentile(a, 15.865)))
