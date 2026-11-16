#!/usr/bin/env python3
# ============================================================================
# peaktail_core.py  --  package peaktail: shared layer (columns of peaktail_dump.C, tail decomposition, helpers)
# ============================================================================
# Written 2026-11-15 (prefix peaktail). The question: how large is the low-mass tail of the eta peak (classes C4 = eta ->
# pi+pi-gamma with the photon converted, C5 = eta -> pi+pi-e+e- direct) in the data, and what does it do to the bump-edge
# measurements f_eta(u-edge) (edgeu) and f_eta(edge) (bedge).
# Conventions are those of edgeu_core.py (imported): data S0, cor_ok == 1, isBest >= 0; simulation mode A owned rows, w x V3
# (code 108), M4' = 0.99615 [m4_t + r (m4 - m4_t)], M(pipi)' = 0.9952 M(pipi), expected raw data rows.
# The "tail" of a simulated peak row is defined by its reconstruction residual d = M4' - 0.99615 m4_true (the generated
# four-body mass is a delta function: no radiative photon, no conversion deficit, see peaktail_sim.py):
#   core(d) = the mirror image of the high side, rho(2 mu - d) for d < mu (mu = mode of d);
#   tail(d) = rho(d) - core(d) on the low side (the asymmetric excess; zero above mu);
#   P_tail(d) = tail(d) / rho(d): the probability that a row with residual d belongs to the tail.
# A peak template with the tail scaled by t keeps the number of peak rows: weight x [ (1 - P_tail) c(t) + t P_tail ],
# c(t) = (N - t N_tail) / N_core (radiation moves rows from the core into the tail).
# ============================================================================
import math
import os
import sys

import numpy as np

sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import edgeu_core as EC

OUT = "/hdd3/przygoda/HADES/FAT_out/all_report/peaktail/"
FIG = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_all/"
NOTE = "/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/work_all/peaktail_results.md"
MPEAK = 0.5439          # data eta peak position (frozen, Part II)
SPEAK_D = 0.0054        # data eta peak width (Part II Gaussian)
SPEAK_S = 0.00609       # simulated eta peak width (WP6/bedge)
ME_ = 0.000511
MPI = 0.13957

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
             mpipi=_col("data_mpipi.f32", np.float32, n), mee=_col("data_mee.f32", np.float32, n), oa=_col("data_oa.f32", np.float32, n),
             mm4=_col("data_mm4.f32", np.float32, n), ep_p=_col("data_ep_p.f32", np.float32, n), em_p=_col("data_em_p.f32", np.float32, n),
             bits=_col("data_bits.u8", np.uint8, n), day=_col("data_day.i16", np.int16, n), fold=_col("data_fold.i8", np.int8, n))
    D["n"] = n
    return D
def load_sim():
    n = layout()["sim"]
    S = dict(tree=_col("sim_tree.i8", np.int8, n), grp=_col("sim_grp.i32", np.int32, n), g=_col("sim_g.i8", np.int8, n),
             cls=_col("sim_cls.i8", np.int8, n), lcls=_col("sim_lcls.i8", np.int8, n), pcls=_col("sim_pcls.i8", np.int8, n),
             bits=_col("sim_bits.u8", np.uint8, n), m4=_col("sim_m4.f32", np.float32, n), m4t=_col("sim_m4true.f32", np.float32, n),
             mpipi=_col("sim_mpipi.f32", np.float32, n), mpipit=_col("sim_mpipitrue.f32", np.float32, n),
             mee=_col("sim_mee.f32", np.float32, n), meet=_col("sim_meetrue.f32", np.float32, n), w=_col("sim_w.f64", np.float64, n),
             f108=_col("sim_f108.f32", np.float32, n), oa=_col("sim_oa.f32", np.float32, n), oat=_col("sim_oatrue.f32", np.float32, n))
    for a in ("ep", "em", "pip", "pim"):
        for c in ("p", "ptrue", "th", "thtrue", "ph"):
            S["%s_%s" % (a, c)] = _col("sim_%s_%s.f32" % (a, c), np.float32, n)
    for c in ("ep_vr", "em_vr", "ep_vz", "em_vz"):
        S[c] = _col("sim_%s.f32" % c, np.float32, n)
    S["n"] = n
    return S
def sim_weight(S, sel="S0", c2corr=True):
    """edgeu_core.sim_weight with the code 108 column of this dump"""
    w = S["w"] * S["f108"]
    b = S["bits"]
    m = ((b & (1 if sel == "S0" else 2)) != 0) & ((b & 8) != 0)
    w = np.where(m, w, 0.0)
    if c2corr:
        c2 = S["g"] == EC.G["C2"]
        w = np.where(c2 & (S["mee"] < 0.010), w * EC.CF_E1, w)
        w = np.where(c2 & ~(S["mee"] < 0.010), w * EC.CF_E2, w)
    return w

# ---------------------------------------------------------------- residual and the tail decomposition
def resid(S, s=EC.S_NOM, r=EC.R_NOM):
    """d = M4' - s m4_true (GeV); M4' = s [m4_t + r (m4 - m4_t)]"""
    m4 = S["m4"].astype(np.float64); mt = S["m4t"].astype(np.float64)
    return s * r * (m4 - mt)
class TailModel:
    """mirror decomposition of the residual distribution rho(d) of the weighted rows (d, w):
    histogram in bw bins on [-lo, +hi]; mode mu from a parabola around the maximum; core(d) = rho(2 mu - d) for d < mu;
    tail = max(rho - core, 0), smoothed; P_tail(d) by linear interpolation of tail/rho at the bin centres."""
    def __init__(self, d, w, bw=0.001, lo=-0.120, hi=0.060, smooth=2):
        self.e = np.arange(lo, hi + 1e-12, bw); self.x = 0.5 * (self.e[1:] + self.e[:-1])
        h, _ = np.histogram(d, bins=self.e, weights=w)
        if smooth:
            k = np.ones(2 * smooth + 1) / (2 * smooth + 1)
            hs = np.convolve(h, k, mode="same")
        else:
            hs = h
        i = int(np.argmax(hs))
        y0, y1, y2 = hs[i - 1], hs[i], hs[i + 1]
        den = (y0 - 2 * y1 + y2)
        self.mu = float(self.x[i] + (0.5 * (y0 - y2) / den * bw if den != 0 else 0.0))
        mir = np.interp(2 * self.mu - self.x, self.x, hs, left=0.0, right=0.0)
        core = np.where(self.x < self.mu, np.minimum(mir, hs), hs)
        tail = np.maximum(hs - core, 0.0)
        self.h, self.hs, self.core, self.tail = h, hs, core, tail
        self.ptail = np.where(hs > 0, tail / np.maximum(hs, 1e-300), 0.0)
        self.ftail = float(tail.sum() / hs.sum())
    def P(self, d):
        return np.interp(d, self.x, self.ptail, left=float(self.ptail[0]), right=0.0)
def tail_weights(w, ptail, t):
    """row weights with the tail scaled by t, the total kept: w [(1 - P) c + t P], c = (N - t Nt)/Nc"""
    N = w.sum(); Nt = (w * ptail).sum(); Nc = N - Nt
    c = (N - t * Nt) / Nc
    return w * ((1.0 - ptail) * c + t * ptail), Nt / N

# ---------------------------------------------------------------- 4-vectors and the linearised residual
def fourvec(p, th, ph, m):
    th = np.radians(th.astype(np.float64)); ph = np.radians(ph.astype(np.float64)); p = p.astype(np.float64)
    px, py, pz = p * np.sin(th) * np.cos(ph), p * np.sin(th) * np.sin(ph), p * np.cos(th)
    return np.sqrt(p * p + m * m), px, py, pz
def m4_and_jacobian(S):
    """M4 from the reconstructed (p, theta, phi) of the four particles and dM4/d|p_i| (i = ep, em, pip, pim)"""
    parts = [("ep", ME_), ("em", ME_), ("pip", MPI), ("pim", MPI)]
    V = [fourvec(S[a + "_p"] * 1e-3 if np.nanmedian(S[a + "_p"]) > 20 else S[a + "_p"], S[a + "_th"], S[a + "_ph"], m) for a, m in parts]
    E = sum(v[0] for v in V); PX = sum(v[1] for v in V); PY = sum(v[2] for v in V); PZ = sum(v[3] for v in V)
    M = np.sqrt(np.maximum(E * E - PX * PX - PY * PY - PZ * PZ, 0))
    J = []
    for (a, m), v in zip(parts, V):
        p = np.sqrt(v[1] ** 2 + v[2] ** 2 + v[3] ** 2)
        beta = p / v[0]
        nx, ny, nz = v[1] / p, v[2] / p, v[3] / p
        J.append((E * beta - (PX * nx + PY * ny + PZ * nz)) / M)
    return M, J

# ---------------------------------------------------------------- small helpers
def wq(a, w, qs):
    """weighted quantiles"""
    o = np.argsort(a); a = a[o]; w = w[o]
    c = np.cumsum(w); c = c / c[-1]
    return [float(np.interp(q, c, a)) for q in qs]
def frac(m, w):
    return float(w[m].sum() / w.sum()) if w.sum() > 0 else float("nan")
