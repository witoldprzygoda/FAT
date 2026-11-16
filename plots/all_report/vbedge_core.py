"""vbedge_core.py -- INDEPENDENT VERIFIER of the package bedge: own edge-step fitter (numpy).

Written 2026-11-15 (prefix vbedge). Re-implemented from the written description of the frozen procedure
(ScanEngine FitEdgeStep / scan_bump.C as documented in the WP6 note): N - CB = pol2 [+ eta Gaussian] + h x step,
step = 0.5 erfc((m - E)/(sqrt(2) s)) averaged uniformly over the M(pipi) slice (21 points), E = 0.99277 M4max(m/0.9970),
s = 5.4 MeV, half range 45 MeV (35 MeV for slices from 0.39), GLS with the CB covariance plus an iterated
Poisson variance max(model + CB, 1). No code of bedge_* is imported; arrays come from vbedge_fill.C.
"""
import numpy as np
from math import erf, erfc, sqrt

OUT = "/hdd3/przygoda/HADES/FAT_out/all_report/vbedge"
NM4, NPP, NCL, NS = 240, 40, 3, 3
M4LO, M4W, PPLO, PPW = 0.35, 0.0025, 0.27, 0.005
NRD, NRS = 200, 100
META, MPI0 = 0.547862, 0.1349768
PEAKM, PEAKS = 0.5439, 0.0054
SCALE, SCALEPP = PEAKM / META, 0.9970
X = M4LO + (np.arange(NM4) + 0.5) * M4W        # bin centres


def load():
    d = {}
    d["N"] = np.fromfile(f"{OUT}/dN.bin").reshape(NS, NCL, NPP, NM4)
    d["Npre"] = np.fromfile(f"{OUT}/dNpre.bin").reshape(NS, NCL, NPP, NM4)
    d["B"] = np.fromfile(f"{OUT}/dB.bin").reshape(NRD, NS, NCL, NPP, NM4)
    d["S"] = np.fromfile(f"{OUT}/sS.bin").reshape(4, 6, NCL, NPP, NM4)     # cfg, grp (C2 C3 T12 otherEta G_D G_C)
    d["R"] = np.fromfile(f"{OUT}/sR.bin").reshape(NRS, 4, 3, NCL, NPP, NM4)  # rep, cfg, grp (C2 C3 T12)
    d["W"] = np.fromfile(f"{OUT}/sW.bin").reshape(2, NCL, 2)
    d["WR"] = np.fromfile(f"{OUT}/sWR.bin").reshape(NRS, 2, NCL, 2)
    return d


def pb(x):
    return int(round((x - PPLO) / PPW))


def m4max(mp):
    E0 = (META ** 2 + MPI0 ** 2 - mp ** 2) / (2 * META)
    p2 = E0 * E0 - MPI0 ** 2
    if p2 < 0:
        return -1.0
    m2 = META ** 2 - META * (E0 - sqrt(p2))
    return sqrt(m2) if m2 > 0 else -1.0


def edge(mp, scale=SCALE):
    e = m4max(mp / SCALEPP)
    return scale * e if e > 0 else -1.0


def step_tmpl(x, lo, hi, s, scale=SCALE, fixed_edge=None):
    q = 1.0 / (sqrt(2.0) * s)
    if fixed_edge is not None:
        return np.array([0.5 * erfc((xx - fixed_edge) * q) for xx in x])
    es = [edge(lo + (i + 0.5) * (hi - lo) / 21, scale) for i in range(21)]
    es = [e for e in es if e > 0]
    return np.array([np.mean([0.5 * erfc((xx - e) * q) for e in es]) for xx in x])


def gfrac(x1, x2, m=PEAKM, s=PEAKS):
    a = 1.0 / (sqrt(2.0) * s)
    return 0.5 * (erf((x2 - m) * a) - erf((x1 - m) * a))


# ------------------------------------------------------------------------------------------- CB
def cb_all(pp, mm, k=1.0, reb=4, mode=0):
    """CB per fine bin and its covariance (2 sqrt(N++ N--) in reb-bin cells, linear interpolation between cell centres)."""
    nw = NM4 // reb
    wp = pp[: nw * reb].reshape(nw, reb).sum(1)
    wm = mm[: nw * reb].reshape(nw, reb).sum(1)
    cw = k * (2 * np.sqrt(wp * wm) if mode == 0 else (wp + wm))
    vw = k * k * (wp + wm)
    t = (X - (M4LO + 0.5 * reb * M4W)) / (reb * M4W)
    j0 = np.clip(np.floor(t).astype(int), 0, nw - 1)
    f = t - np.floor(t)
    j1 = np.minimum(j0 + 1, nw - 1)
    low = t <= 0
    high = t >= nw - 1
    j0 = np.where(low, 0, np.where(high, nw - 1, j0))
    j1 = np.where(low, 0, np.where(high, nw - 1, j1))
    f = np.where(low | high, 0.0, f)
    A = np.zeros((NM4, nw))
    A[np.arange(NM4), j0] += (1 - f) / reb
    A[np.arange(NM4), j1] += f / reb
    cb = A @ cw
    cov = A @ np.diag(vw) @ A.T
    return cb, cov, A, cw, vw


def spread_frac(all_pp, all_mm, cl_pp, cl_mm):
    """class share of like-sign rows per 50 MeV band (bands from 0.30), per fine bin."""
    band = np.floor((X - 0.30) / 0.05 + 1e-9).astype(int)
    fr = np.zeros(NM4)
    for b in np.unique(band):
        sel = band == b
        la = (all_pp[sel] + all_mm[sel]).sum()
        ls = (cl_pp[sel] + cl_mm[sel]).sum()
        fr[sel] = ls / la if la > 0 else 0.0
    return fr


def data_spec(N, cls, pa, pbb, cbopt=None):
    """N: [3 samples][3 cls][40 pp][240]; cls -1 = all. Returns y = N - CB, cb, covCB."""
    cbopt = cbopt or {}
    u = N[0, :, pa:pbb].sum(axis=(0, 1)); p = N[1, :, pa:pbb].sum(axis=(0, 1)); m = N[2, :, pa:pbb].sum(axis=(0, 1))
    cb, cov, *_ = cb_all(p, m, **cbopt)
    if cls < 0:
        return u - cb, cb, cov, u
    uc = N[0, cls, pa:pbb].sum(0); pc = N[1, cls, pa:pbb].sum(0); mc = N[2, cls, pa:pbb].sum(0)
    fr = spread_frac(p, m, pc, mc)
    cbc = cb * fr
    covc = cov * np.outer(fr, fr)
    return uc - cbc, cbc, covc, uc


# ------------------------------------------------------------------------------------------- fit
class Cfg:
    def __init__(self, lo, hi, hr=None, order=2, s=PEAKS, scale=SCALE, fixed_edge=None, dlo=0.0, dhi=0.0):
        self.lo, self.hi = lo, hi
        self.hr = hr if hr is not None else (0.035 if lo >= 0.39 - 1e-9 else 0.045)
        self.order, self.s, self.scale, self.fixed_edge = order, s, scale, fixed_edge
        self.dlo, self.dhi = dlo, dhi
        self._prep()

    def _prep(self):
        e = self.fixed_edge if self.fixed_edge is not None else edge(0.5 * (self.lo + self.hi), self.scale)
        lo, hi = e - self.hr - self.dlo, e + self.hr + self.dhi
        if hi <= PEAKM - 3 * PEAKS:
            we = False
        elif hi < PEAKM - PEAKS:
            hi = PEAKM - 3 * PEAKS; we = False
        else:
            hi = max(hi, PEAKM + 3 * PEAKS); we = True
        lo = M4LO + M4W * np.floor((lo - M4LO) / M4W + 0.5)
        hi = M4LO + M4W * np.floor((hi - M4LO) / M4W + 0.5)
        self.idx = np.where((X >= lo - 1e-12) & (X < hi - 1e-12))[0]
        x = X[self.idx]
        x0, dx = 0.5 * (lo + hi), 0.5 * (hi - lo)
        t = (x - x0) / dx
        cols = [t ** j for j in range(self.order + 1)]
        if we:
            cols.append(np.array([gfrac(xx - 0.5 * M4W, xx + 0.5 * M4W) for xx in x]))
        st = step_tmpl(x, self.lo, self.hi, self.s, self.scale, self.fixed_edge)
        cols.append(st)
        self.F = np.column_stack(cols)
        self.st = st
        self.edge, self.rlo, self.rhi, self.we = e, lo, hi, we


def gls(y, V, F):
    L = np.linalg.cholesky(V)
    Fi = np.linalg.solve(L, F); yi = np.linalg.solve(L, y)
    C = np.linalg.inv(Fi.T @ Fi)
    p = C @ (Fi.T @ yi)
    r = yi - Fi @ p
    return p, C, float(r @ r)


def fit(y_full, cb_full, cov_full, cfg, Vfix=None):
    """returns h, eh, chi2, ndf, final V (npts x npts), par"""
    i = cfg.idx
    y = y_full[i]; cb = cb_full[i] if cb_full is not None else np.zeros(len(i))
    V0 = cov_full[np.ix_(i, i)] if cov_full is not None else np.zeros((len(i), len(i)))
    F = cfg.F
    if Vfix is not None:
        p, C, chi2 = gls(y, Vfix, F)
        return p[-1], sqrt(C[-1, -1]), chi2, len(i) - F.shape[1], Vfix, p
    var = np.maximum(y + cb, 1.0)
    p, C, chi2 = gls(y, V0 + np.diag(var), F)
    for _ in range(6):
        hold = p[-1]
        var = np.maximum(F @ p + cb, 1.0)
        p, C, chi2 = gls(y, V0 + np.diag(var), F)
        if abs(p[-1] - hold) < 1e-4 * sqrt(C[-1, -1]):
            break
    return p[-1], sqrt(C[-1, -1]), chi2, len(i) - F.shape[1], V0 + np.diag(var), p


NOMSET = [(0.30, 0.34), (0.34, 0.37), (0.37, 0.39), (0.39, 0.40)]
