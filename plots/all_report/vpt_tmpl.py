#!/usr/bin/env python3
"""vpt_tmpl.py -- verifier of peaktail: per-row quantities of the simulated peak rows (C4, C5) and the tail-scaled
templates K(t) = K_core (N - t N_tail)/N_core + t K_tail (own code; tail flags by SUBSTITUTION of true momenta).

Tail definitions (a row of C4/C5 is a tail row when):
  'sub5'  : the leptons' momentum deficit, computed by putting the true |p| of both leptons into the reconstructed
            four-vectors (reco angles kept), lowers M4' by more than 5 MeV/c2;
  'sub10' : the same, by more than 10 MeV/c2;
  'p95'   : a lepton has p_reco/p_true < 0.95 (no mass involved);
  'res15' : the residual d = M4' - 0.99615 m4_true is below -15 MeV/c2 (pure line-shape definition; resolution included).
"""
import sys
import numpy as np
sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import vpt_core as V

_cache = {}


def rows():
    if "R" in _cache:
        return _cache["R"]
    S = V.load_sim()
    w = V.sim_weight(S)                       # S0, mode A, w x V3, C2 x c (C2 only)
    m45 = ((S["g"] == 2) | (S["g"] == 3)) & (w > 0)
    mass = {"ep": V.ME_E, "em": V.ME_E, "pip": V.MPI, "pim": V.MPI}
    vr = {n: V.p4(S[n + "_p"], S[n + "_th"], S[n + "_ph"], mass[n]) for n in mass}
    vt = {n: V.p4(S[n + "_pt"], S[n + "_th"], S[n + "_ph"], mass[n]) for n in mass}
    m4r = V.minv(vr["ep"], vr["em"], vr["pip"], vr["pim"])
    m4l = V.minv(vt["ep"], vt["em"], vr["pip"], vr["pim"])
    dlep = V.S_NOM * (m4r - m4l)
    d = V.S_NOM * (S["m4"].astype(float) - S["m4t"].astype(float))
    rat = np.minimum(S["ep_p"] / S["ep_pt"], S["em_p"] / S["em_pt"]).astype(float)
    mpp = V.SPP * S["pp"].astype(float)
    m4p = V.sim_m4p(S)
    u = m4p - V.edge(mpp, S["ee"].astype(float))
    u[~(S["ee"] < V.MPI0)] = np.nan
    R = dict(S=S, w=w, m45=m45, dlep=dlep, d=d, rat=rat, mpp=mpp, m4p=m4p, u=u)
    R["flag"] = {"sub5": dlep < -0.005, "sub10": dlep < -0.010, "p95": rat < 0.95, "res15": d < -0.015}
    _cache["R"] = R
    return R


def in_pp(mpp, a, b):
    m = mpp.astype(np.float32)
    return (m >= np.float32(a)) & (m < np.float32(b))


def templates(x, edges, region, flag, w=None, cls_sel=None, extra=None):
    """per class (C4, C5): K_all, K_core, K_tail histograms of variable x (array per row) and the class sums in the
    region (all x).  Returns dict cls -> (Kall, Kcore, Ktail, N, Ncore, Ntail)."""
    R = rows(); S = R["S"]
    w = R["w"] if w is None else w
    out = {}
    for c, gi in (("C4", 2), ("C5", 3)):
        m = R["m45"] & (S["g"] == gi) & region
        if extra is not None:
            m &= extra
        tl = m & flag; co = m & ~flag
        N, Nt = w[m].sum(), w[tl].sum()
        Ka = np.histogram(x[m], bins=edges, weights=w[m])[0]
        Kc = np.histogram(x[co], bins=edges, weights=w[co])[0]
        Kt = np.histogram(x[tl], bins=edges, weights=w[tl])[0]
        out[c] = (Ka, Kc, Kt, N, N - Nt, Nt)
    return out


def K_of_t(T, tC4, tC5=None):
    tC5 = tC4 if tC5 is None else tC5
    K = 0
    for c, t in (("C4", tC4), ("C5", tC5)):
        Ka, Kc, Kt, N, Nc, Nt = T[c]
        K = K + Kc * (N - t * Nt) / Nc + t * Kt
    return K


def dK(T, which=("C4", "C5")):
    """dK/dt (K linear in t)"""
    out = 0
    for c in which:
        Ka, Kc, Kt, N, Nc, Nt = T[c]
        out = out + Kt - Kc * Nt / Nc
    return out


def group_boot(n, seed=20261115):
    """Poisson(1) weights per simulated event group (tree, grp) for the C4/C5 rows: (n x rows) on the full row set"""
    R = rows(); S = R["S"]
    key = S["tree"].astype(np.int64) * (1 << 32) + S["grp"].astype(np.int64)
    sel = R["m45"]
    uq, inv = np.unique(key[sel], return_inverse=True)
    rng = np.random.default_rng(seed)
    k = rng.poisson(1.0, size=(n, uq.size)).astype(np.float64)
    W = np.zeros((n, len(key)))
    W[:, sel] = k[:, inv]
    return W


def smooth_c4(K=400, seed=314, pool_mee=0.020, ppw=0.010):
    """Smoothed C4 entries: every C4 row i (S0, mode A) is replaced by K entries x = s m4_true,i + d5_k, where d5_k is the
    residual (M4' - s m4_true) of a C5 row drawn (prob. ~ weight) from the C5 rows with M(ee) < pool_mee in the same ppw
    M(pipi)' bin; u = x - E_i (row i's own reco edge); tail flags (sub5, p95, sub10, res15) taken from the C5 row.
    Returns dict with arrays (row index i, m4, u, mpp, w, flags)."""
    if ("smooth", K, seed, pool_mee, ppw) in _cache:
        return _cache[("smooth", K, seed, pool_mee, ppw)]
    R = rows(); S = R["S"]
    rng = np.random.default_rng(seed)
    c4 = np.where(R["m45"] & (S["g"] == 2))[0]
    c5 = R["m45"] & (S["g"] == 3) & (S["ee"] < pool_mee) & (S["m4t"] > 0)
    b = np.floor(R["mpp"] / ppw).astype(int)
    E = V.edge(R["mpp"], S["ee"].astype(float))
    out = {k: [] for k in ("i", "m4", "u", "mpp", "w", "sub5", "p95", "sub10", "res15")}
    for bb in np.unique(b[c4]):
        rows4 = c4[b[c4] == bb]
        pool = np.where(c5 & (b == bb))[0]
        if pool.size < 50:
            pool = np.where(c5 & (np.abs(b - bb) <= 1))[0]
        pw = R["w"][pool] / R["w"][pool].sum()
        for i in rows4:
            k = rng.choice(pool, size=K, p=pw)
            x = V.S_NOM * float(S["m4t"][i]) + R["d"][k]
            out["i"].append(np.full(K, i)); out["m4"].append(x)
            uu = x - E[i] if S["ee"][i] < V.MPI0 else np.full(K, np.nan)
            out["u"].append(uu); out["mpp"].append(np.full(K, R["mpp"][i])); out["w"].append(np.full(K, R["w"][i] / K))
            for f in ("sub5", "p95", "sub10", "res15"):
                out[f].append(R["flag"][f][k])
    o = {k: np.concatenate(v) for k, v in out.items()}
    _cache[("smooth", K, seed, pool_mee, ppw)] = o
    return o


def templates_smooth(var, edges, pp_lo, pp_hi, flag, wC5=None, wC4row=None, extra5=None, extra4=None, sm=None):
    """as templates(), but with C4 from smooth_c4(). var: 'u' or 'm4'. wC5: optional replacement weights (all rows) for C5;
    wC4row: optional per-row weights (all rows) applied to the C4 parent rows (bootstrap)."""
    R = rows(); S = R["S"]
    sm = smooth_c4() if sm is None else sm
    x5 = R["u"] if var == "u" else R["m4p"]
    reg5 = in_pp(R["mpp"], pp_lo, pp_hi) & R["m45"] & (S["g"] == 3)
    if var == "u":
        reg5 &= np.isfinite(R["u"])
    if extra5 is not None:
        reg5 &= extra5
    w5 = R["w"] if wC5 is None else wC5
    fl5 = R["flag"][flag]
    out = {}
    tl = reg5 & fl5; co = reg5 & ~fl5
    out["C5"] = (np.histogram(x5[reg5], edges, weights=w5[reg5])[0], np.histogram(x5[co], edges, weights=w5[co])[0],
                 np.histogram(x5[tl], edges, weights=w5[tl])[0], w5[reg5].sum(), w5[co].sum(), w5[tl].sum())
    x4 = sm[var]
    reg4 = in_pp(sm["mpp"], pp_lo, pp_hi)
    if var == "u":
        reg4 &= np.isfinite(x4)
    if extra4 is not None:
        reg4 &= extra4[sm["i"]]
    w4 = sm["w"] if wC4row is None else sm["w"] * wC4row[sm["i"]]
    f4 = sm[flag]
    tl = reg4 & f4; co = reg4 & ~f4
    out["C4"] = (np.histogram(x4[reg4], edges, weights=w4[reg4])[0], np.histogram(x4[co], edges, weights=w4[co])[0],
                 np.histogram(x4[tl], edges, weights=w4[tl])[0], w4[reg4].sum(), w4[co].sum(), w4[tl].sum())
    return out
