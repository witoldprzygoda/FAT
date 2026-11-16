#!/usr/bin/env python3
# ============================================================================
# vflex_core.py -- independent verification of package flex (prefix vflex): inputs, smoothing, models, estimator
# ============================================================================
# Written from the Part III / flex conventions as documented (not from the flex code paths):
#   inputs  : hists_v1.root F10 (S0, all, isBest >= 0) and wp5_egcov.root, dumped by vflex_dump.C (plain text);
#   template: SMASH sigma x PDG BR x V3 [nb per bin] x L / (64/63) = expected RAW data rows;
#             eta classes C2..C5, G_eta raw with BB-lite variance; non-eta G_D, G_C (and the omega/rho0 split) smoothed
#             by exp(pol n) in x = (M4 - 0.55)/0.1 with a scaled-Poisson likelihood (one scale V/C per histogram),
#             n by N_eff = C^2/V (>= 50: 3, >= 20: 2, >= 8: 1, else 0), parent shape scaled if N_eff < 20;
#             bin integrals by a 16-point Gauss-Legendre rule (flex: 8 points), maximisation by scipy trust-exact;
#   model   : mu = k CB + sum_t s_t f_t T_t, V = max(mu, 1) + k^2 var(CB) + sum_raw s_t^2 var(T_t);
#   estimator: the quasi-likelihood score equations J^T (y - mu) / V - P (theta - theta0) = 0 solved by a ROOT FINDER
#             (scipy.optimize.root, MINPACK hybrd/lm) -- not by Fisher scoring as in flex_lib.h;
#   errors  : sandwich with the event-group covariance scaled to the model, the CB replica covariance and the constraint;
#             MC = spread of the refits with the 100 WP0 bootstrap replicas (smoothing redone, BB variance of the full sample).
# ============================================================================
import os
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import optimize, linalg

DIR = "/hdd3/przygoda/HADES/FAT_out/all_report/vflex"
NB2RAW = 6079.2 * 63.0 / 64.0
NS, NX, NREPC, NREPS = 5, 70, 200, 100
COMPS = ["C2", "C3", "C4", "C5", "Geta", "GD", "Gconv", "Gpi0D", "Gvec", "Goth", "Gpi0D_om", "Gvec_om", "Goth_om",
         "Gpi0D_nom", "Gvec_nom", "Goth_nom", "Gconv_om", "Gconv_nom", "T3"]
# fit components: name -> (parts, smoothed, parent parts)
FCOMP = {
    "C2": (["C2"], False, None), "C3": (["C3"], False, None), "C4": (["C4"], False, None), "C5": (["C5"], False, None),
    "Geta": (["Geta"], False, None),
    "GD": (["GD"], True, None), "GC": (["Gconv"], True, None),
    "GD_om": (["Gpi0D_om", "Gvec_om", "Goth_om"], True, ["GD"]),
    "GD_nom": (["Gpi0D_nom", "Gvec_nom", "Goth_nom"], True, ["GD"]),
    "GC_om": (["Gconv_om"], True, ["Gconv"]), "GC_nom": (["Gconv_nom"], True, ["Gconv"]),
    "T3": (["T3"], True, None),
}


# ---------------------------------------------------------------------------- inputs
def load():
    cache = DIR + "/vflex_inputs.npz"
    if os.path.exists(cache):
        z = np.load(cache)
        return {k: z[k] for k in z.files}
    d = np.loadtxt(DIR + "/data_F10.txt")
    N = np.zeros((NS, 2, NX)); NE = np.zeros_like(N); CB = np.zeros_like(N)
    lo = np.zeros(NX); hi = np.zeros(NX)
    for r in d:
        j, e, x = int(r[0]), int(r[1]), int(r[2]) - 1
        N[j, e, x], NE[j, e, x], CB[j, e, x] = r[3], r[4], r[5]
        lo[x], hi[x] = r[6], r[7]
    q = np.loadtxt(DIR + "/cbrep_F10.txt")
    CBR = np.zeros((NS, 2, NX, NREPC))
    for r in q:
        CBR[int(r[0]), int(r[1]), int(r[2]) - 1] = r[3:]
    EG = np.loadtxt(DIR + "/egcov_F10.txt")
    SC = np.zeros((len(COMPS), NS, 2, NX)); SE = np.zeros_like(SC); SR = np.zeros((len(COMPS), NS, 2, NX, NREPS))
    with open(DIR + "/sim_F10.txt") as f:
        for line in f:
            t = line.split(); c = COMPS.index(t[0]); j, e, x = int(t[1]), int(t[2]), int(t[3]) - 1
            SC[c, j, e, x] = float(t[4]); SE[c, j, e, x] = float(t[5])
    with open(DIR + "/simrep_F10.txt") as f:
        for line in f:
            t = line.split(); c = COMPS.index(t[0]); j, e, x = int(t[1]), int(t[2]), int(t[3]) - 1
            SR[c, j, e, x] = np.array(t[4:], dtype=float)
    out = dict(N=N, NE=NE, CB=CB, lo=lo, hi=hi, CBR=CBR, EG=EG, SC=SC * NB2RAW, SV=(SE * NB2RAW) ** 2, SR=SR * NB2RAW)
    np.savez(cache, **out)
    return out


# ---------------------------------------------------------------------------- smoothing
GLX, GLW = leggauss(16)


class Smooth:
    def __init__(self):
        self.order, self.p, self.neff, self.scaled = -1, np.zeros(4), 0.0, False

    def bins(self, e):
        if self.order < 0:
            return np.zeros(len(e) - 1)
        h = 0.5 * (e[1:] - e[:-1]); m = 0.5 * (e[1:] + e[:-1])
        X = ((m[:, None] + h[:, None] * GLX[None, :]) - 0.55) / 0.1
        P = self.p
        return (h[:, None] * GLW[None, :] * np.exp(P[0] + X * (P[1] + X * (P[2] + X * P[3])))).sum(1)


def smooth_fit(c, v, e, parent=None):
    S = Smooth()
    C, V = c.sum(), v.sum()
    if not C > 0:
        return S
    S.neff = C * C / V if V > 0 else 0.0
    if S.neff < 20 and parent is not None and parent.order >= 0:
        P = Smooth(); P.order, P.p, P.neff, P.scaled = parent.order, parent.p.copy(), S.neff, True
        P.p[0] += np.log(C / parent.bins(e).sum())
        return P
    order = 3 if S.neff >= 50 else 2 if S.neff >= 20 else 1 if S.neff >= 8 else 0
    S.order = order
    sbar = V / C if V > 0 else 1.0
    nn = c / sbar
    h = 0.5 * (e[1:] - e[:-1]); m = 0.5 * (e[1:] + e[:-1])
    X = ((m[:, None] + h[:, None] * GLX[None, :]) - 0.55) / 0.1
    W = h[:, None] * GLW[None, :]
    PW = X[..., None] ** np.arange(order + 1)          # (n, q, k)

    def parts(p):
        ex = W * np.exp(PW @ p)                          # (n, q)
        f = ex.sum(1)
        df = np.einsum('nq,nqk->nk', ex, PW)
        d2f = np.einsum('nq,nqk,nql->nkl', ex, PW, PW)
        return f, df, d2f

    def fun(p):
        f, df, d2f = parts(p)
        mm = np.maximum(f / sbar, 1e-300)
        L = np.sum(np.where(nn > 0, nn * np.log(mm), 0.0)) - mm.sum()
        g = ((nn / mm - 1.0)[:, None] * df / sbar).sum(0)
        H = (-(nn / mm ** 2)[:, None, None] * df[:, :, None] * df[:, None, :] / sbar ** 2
             + (nn / mm - 1.0)[:, None, None] * d2f / sbar).sum(0)
        return -L, -g, -H

    p0 = np.zeros(order + 1); p0[0] = np.log(C / (e[-1] - e[0]))
    r = optimize.minimize(lambda p: fun(p)[0], p0, jac=lambda p: fun(p)[1], hess=lambda p: fun(p)[2],
                          method='trust-exact', options=dict(gtol=1e-10, maxiter=500))
    S.p = np.zeros(4); S.p[:order + 1] = r.x
    return S


# ---------------------------------------------------------------------------- fit bins (E2 or E1+E2), templates
class Layout:
    def __init__(self, D, classes=(1,), mlo=0.40, mhi=0.70):
        self.classes = list(classes)
        lo, hi = D['lo'], D['hi']
        xs = [x for x in range(NX) if lo[x] >= mlo - 1e-9 and hi[x] <= mhi + 1e-9]
        self.idx = [(e, j, x) for e in self.classes for j in range(NS) for x in xs]
        self.n = len(self.idx)
        self.e = np.array([t[0] for t in self.idx]); self.j = np.array([t[1] for t in self.idx]); self.x = np.array([t[2] for t in self.idx])
        self.lo = lo[self.x]; self.hi = hi[self.x]
        self.xc = (0.5 * (self.lo + self.hi) - 0.55) / 0.1
        self.edges = np.append(lo, hi[-1])
        self.y = np.array([D['N'][j, e, x] for (e, j, x) in self.idx])
        self.cb = np.array([D['CB'][j, e, x] for (e, j, x) in self.idx])
        cells = np.array([(2 * j + e) * NX + x for (e, j, x) in self.idx])
        self.eg = D['EG'][np.ix_(cells, cells)]
        R = np.array([D['CBR'][j, e, x] for (e, j, x) in self.idx])    # (n, 200)
        self.cbrep = R
        Rm = R - R.mean(1, keepdims=True)
        self.ccb = Rm @ Rm.T / R.shape[1]                                # 1/n as wp5
        self.vcb = np.diag(self.ccb).copy()
        self.xs = xs


def comp_hist(D, name, j, e, rep=-1):
    """content (70) and variance (70) of a fit component in histogram (j, e); rep >= 0: bootstrap replica content"""
    parts, sm, par = FCOMP[name]
    c = np.zeros(NX); v = np.zeros(NX)
    for p in parts:
        k = COMPS.index(p)
        c += D['SC'][k, j, e] if rep < 0 else D['SR'][k, j, e, :, rep]
        v += D['SV'][k, j, e]
    return c, v


def template(D, L, name, rep=-1, info=None):
    """per-bin template values and BB variance of a component for a layout (smoothed per fit histogram on 0.40-0.70)"""
    parts, sm, par = FCOMP[name]
    t = np.zeros(L.n); v = np.zeros(L.n)
    E = L.edges
    for e in L.classes:
        for j in range(NS):
            c, vv = comp_hist(D, name, j, e, rep)
            if sm:
                P = None
                if par is not None:
                    pc = np.zeros(NX); pv = np.zeros(NX)
                    for p in par:
                        k = COMPS.index(p)
                        pc += D['SC'][k, j, e] if rep < 0 else D['SR'][k, j, e, :, rep]
                        pv += D['SV'][k, j, e]
                    P = smooth_fit(pc, pv, E)
                S = smooth_fit(c, vv, E, P)
                if info is not None:
                    info[(name, j, e)] = (S.order, S.neff, S.scaled)
                vals = S.bins(E); var = np.zeros(NX)
            else:
                vals, var = c, vv
            sel = (L.e == e) & (L.j == j)
            t[sel] = vals[L.x[sel]]; v[sel] = var[L.x[sel]]
    return t, v


# ---------------------------------------------------------------------------- models
class FlexModel:
    """mu = k CB + sum_t s_t f_t T_t, f_t = 1 + b_j x for background terms if tilt (0 none, 1 common, 2 per slice)"""
    def __init__(self, comps, tilt=0, ksig=0.05, ties=None):
        self.comps, self.tilt, self.ksig = list(comps), tilt, ksig
        self.bg = np.array([FCOMP[c][1] for c in comps])
        self.names = ["s_" + c for c in comps] + (["b"] if tilt == 1 else []) + ([f"b_{j + 1}" for j in range(NS)] if tilt == 2 else []) + ["k_E2"]
        self.np = len(self.names)
        self.start = np.array([1.0] * len(comps) + [0.0] * (1 if tilt == 1 else 5 if tilt == 2 else 0) + [1.0])
        self.cons = {self.np - 1: (1.0, ksig)}
        self.nscale = len(comps)

    def eval(self, th, L, T, Vt):
        m = len(self.comps)
        s = th[:m]; k = th[-1]
        if self.tilt == 0:
            tl = np.ones(L.n)
        elif self.tilt == 1:
            tl = 1 + th[m] * L.xc
        else:
            tl = 1 + th[m + L.j] * L.xc
        F = np.where(self.bg[None, :], tl[:, None], 1.0)               # (n, m)
        mu = k * L.cb + (s[None, :] * F * T).sum(1)
        J = np.zeros((L.n, self.np))
        J[:, :m] = F * T
        bgsum = (s[None, :] * T * self.bg[None, :]).sum(1) * L.xc       # d mu / d b
        if self.tilt == 1:
            J[:, m] = bgsum
        elif self.tilt == 2:
            for j in range(NS):
                J[:, m + j] = np.where(L.j == j, bgsum, 0.0)
        J[:, -1] = L.cb
        vbb = ((s[None, :] * F) ** 2 * Vt).sum(1)
        bgpart = (s[None, :] * F * T * self.bg[None, :]).sum(1)
        return mu, J, vbb, k, bgpart


class P3Model:
    """Part III F-main: a_eta [C2 + kappa (C3 + C4) + rho_d C5 + g G_eta] + c_j (1 + b x) [G_D + kappa_B G_C] + k CB"""
    comps = ["C2", "C3", "C4", "C5", "Geta", "GD", "GC"]

    def __init__(self):
        self.names = ["a_eta", "kappa", "rho_d", "g"] + [f"c_{j + 1}" for j in range(NS)] + ["b", "kappa_B", "k_E2"]
        self.np = len(self.names)
        self.start = np.array([0.4, 1, 1, 1, 1, 1, 1, 1, 1, -0.07, 2.0, 1.0])
        self.cons = {1: (1.0, 0.2), 3: (1.0, 0.3), 11: (1.0, 0.05)}

    def eval(self, th, L, T, Vt):
        a, ka, rd, g = th[:4]; c = th[4 + L.j]; b = th[9]; kb = th[10]; k = th[11]
        C2, C3, C4, C5, Ge, GD, GC = T.T
        v2, v3, v4, v5, ve = Vt[:, 0], Vt[:, 1], Vt[:, 2], Vt[:, 3], Vt[:, 4]
        tl = 1 + b * L.xc
        eta = C2 + ka * (C3 + C4) + rd * C5 + g * Ge
        bg = GD + kb * GC
        mu = k * L.cb + a * eta + c * tl * bg
        J = np.zeros((L.n, self.np))
        J[:, 0] = eta; J[:, 1] = a * (C3 + C4); J[:, 2] = a * C5; J[:, 3] = a * Ge
        for j in range(NS):
            J[:, 4 + j] = np.where(L.j == j, tl * bg, 0.0)
        J[:, 9] = c * L.xc * bg; J[:, 10] = c * tl * GC; J[:, 11] = L.cb
        vbb = a ** 2 * (v2 + ka ** 2 * (v3 + v4) + rd ** 2 * v5 + g ** 2 * ve)
        return mu, J, vbb, k, c * tl * bg


def templates(D, L, model, rep=-1, rawbg=False, info=None):
    T = np.zeros((L.n, len(model.comps))); V = np.zeros_like(T)
    for a, c in enumerate(model.comps):
        if rawbg and FCOMP[c][1]:
            for e in L.classes:
                for j in range(NS):
                    cc, vv = comp_hist(D, c, j, e, rep)
                    sel = (L.e == e) & (L.j == j)
                    T[sel, a] = cc[L.x[sel]]; V[sel, a] = vv[L.x[sel]]
        else:
            t, v = template(D, L, c, rep, info)
            T[:, a] = t; V[:, a] = v          # raw: full-sample variance also for a replica; smoothed: 0
    return T, V


# ---------------------------------------------------------------------------- estimator
def score(th, model, L, y, T, Vt, cons, lam=1.0):
    mu, J, vbb, k, _ = model.eval(th, L, T, Vt)
    V = np.maximum(mu, 1.0) + k * k * L.vcb + lam * vbb
    g = J.T @ ((y - mu) / V)
    for i, (c0, s) in cons.items():
        g[i] -= (th[i] - c0) / s ** 2
    return g


def data_cov(y, mu, egcov, floor=0.0):
    """event-group covariance of N scaled to the model, eg_ab sqrt(mu_a mu_b / (y_a y_b)); empty bins: max(mu, floor) on
    the diagonal (fit: floor 0, hybrid GLS: floor 1, as flex_lib.h); egcov None: Poisson diag(mu)"""
    if egcov is None:
        return np.diag(np.maximum(mu, floor))
    sf = np.where(y > 0, np.sqrt(np.maximum(mu, 1e-9) / np.where(y > 0, y, 1)), 0.0)
    SN = egcov * np.outer(sf, sf)
    bad = y <= 0
    if bad.any():
        SN[bad, :] = 0; SN[:, bad] = 0
        SN[bad, bad] = np.maximum(mu[bad], floor)
    return SN


def fit(model, L, T, Vt, y=None, cb=None, start=None, cons=None, egcov=None, sandwich=True, fixed=None):
    """root of the quasi-likelihood score; returns dict with theta, mu, V, Vd, A^-1, sandwich parts"""
    import copy
    if y is None:
        y = L.y
    if cb is not None:
        L = copy.copy(L); L.cb = cb
    cons = model.cons if cons is None else cons
    th0 = model.start.copy() if start is None else np.array(start, float)
    free = np.ones(model.np, bool)
    if fixed:
        for i, v in fixed.items():
            free[i] = False; th0[i] = v

    def g(tf, lam=1.0):
        th = th0.copy(); th[free] = tf
        return score(th, model, L, y, T, Vt, cons, lam)[free]

    # scale the problem by the Fisher information at the start (the score has very different units per parameter)
    mu, J, vbb, k, _ = model.eval(th0, L, T, Vt)
    V = np.maximum(mu, 1.0) + k * k * L.vcb + vbb
    A0 = (J[:, free] / V[:, None]).T @ J[:, free]
    sc = 1.0 / np.sqrt(np.maximum(np.diag(A0), 1e-30))
    # 1) root of the score by a bounded trust-region least-squares solver (scipy least_squares 'trf' on the scaled score);
    #    the bounds exclude the spurious roots at |s| -> infinity (there the BB-lite variance makes every weight vanish)
    lo = np.where(np.arange(model.np) < getattr(model, 'nscale', model.np), -30.0, -5.0)[free]
    hi = -lo
    lo = np.minimum(lo, th0[free] - 1); hi = np.maximum(hi, th0[free] + 1)
    #    homotopy in the BB-lite variance (lambda = 0 -> 1): at lambda = 0 there is no root at infinity; each step starts
    #    from the previous root
    tf = th0[free].copy(); nfev = 0
    for lam in (0.0, 0.25, 0.5, 0.75, 1.0):
        r = optimize.least_squares(lambda z: g(z, lam) * sc, np.clip(tf, lo + 1e-9, hi - 1e-9), bounds=(lo, hi), method='trf',
                                   xtol=1e-15, ftol=1e-15, gtol=1e-15, x_scale=sc, max_nfev=4000)
        tf = r.x; nfev += r.nfev
    r.nfev = nfev
    # 2) polish / cross-check by the frozen-variance fixed-point iteration (V at the current point, exact weighted linear
    #    least squares of the linearised model; for the scale-only models this is exact)
    for it in range(200):
        th = th0.copy(); th[free] = tf
        mu, J, vbb, k, _ = model.eval(th, L, T, Vt)
        Vv = np.maximum(mu, 1.0) + k * k * L.vcb + vbb
        Jf = J[:, free]
        Pd = np.zeros(free.sum()); P0 = np.zeros(free.sum())
        for i, (c0, s_) in cons.items():
            if free[i]:
                q = list(np.where(free)[0]).index(i); Pd[q] = 1 / s_ ** 2; P0[q] = c0
        A = (Jf / Vv[:, None]).T @ Jf + np.diag(Pd)
        b = Jf.T @ ((y - mu) / Vv) - Pd * (tf - P0)
        d = linalg.solve(A, b, assume_a='pos')
        tf = tf + d
        if np.max(np.abs(d) * np.sqrt(np.diag(A))) < 1e-9:
            break
    ok = bool(np.max(np.abs(g(tf) * sc)) < 1e-6 and np.all(np.abs(tf) < 1e3))
    res_ls = r.x
    th = th0.copy(); th[free] = tf
    mu, J, vbb, k, bgp = model.eval(th, L, T, Vt)
    Vd = np.maximum(mu, 1.0) + k * k * L.vcb
    V = Vd + vbb
    Jf = J[:, free]
    P = np.zeros(free.sum())
    fidx = np.where(free)[0]
    for i, (c0, s) in cons.items():
        if free[i]:
            P[list(fidx).index(i)] = 1.0 / s ** 2
    A = (Jf / V[:, None]).T @ Jf + np.diag(P)
    Ai = linalg.inv(A)
    res = dict(ok=ok, th=th, mu=mu, V=V, Vd=Vd, free=free, Ai=Ai, k=k, bg=bgp, L=L, y=y, ls_th=res_ls, ls_nfev=r.nfev, fp_it=it, eg=egcov)
    res['pearson'] = np.sum((y - mu) ** 2 / V) + sum(((th[i] - c0) / s) ** 2 for i, (c0, s) in cons.items())
    if sandwich:
        U = Jf / V[:, None]
        SN = data_cov(y, mu, egcov)
        kk = k * np.ones(L.n)
        SC = L.ccb * np.outer(kk, kk)
        def emb(M):
            O = np.zeros((model.np, model.np)); O[np.ix_(free, free)] = M; return O
        Cn = emb(Ai @ (U.T @ SN @ U) @ Ai); Cc = emb(Ai @ (U.T @ SC @ U) @ Ai); Cp = emb(Ai @ np.diag(P) @ Ai)
        res.update(Cn=Cn, Ccb=Cc, Ccons=Cp, Cstat=Cn + Cc + Cp, SN=SN, SC=SC)
    return res


def hybrid_gls(model, L, T, Vt, R, reps_T):
    """data covariance (event groups scaled to the model + CB replicas) + BB diagonal + replica covariance of the smoothed
    part of the prediction at the fitted theta"""
    n = L.n
    y, mu = R['y'], R['mu']
    W = data_cov(y, mu, R['eg'], floor=1.0) + R['SC']
    Wd = W.copy()
    W = W + np.diag(R['V'] - R['Vd'])
    bgs = np.array([model.eval(R['th'], L, Tr, Vt)[4] for Tr in reps_T])
    if len(bgs) > 1:
        Bm = bgs - bgs.mean(0)
        W = W + Bm.T @ Bm / (len(bgs) - 1)
    r = y - mu
    cz = sum(((R['th'][i] - c0) / s) ** 2 for i, (c0, s) in model.cons.items())
    chi2 = r @ linalg.solve(W, r, assume_a='pos') + cz
    chi2d = r @ linalg.solve(Wd, r, assume_a='pos') + cz
    return chi2, chi2d


def replica_fits(model, L, T_reps, Vt, th0, y=None, cb=None, cons=None):
    out = []
    for Tr in T_reps:
        r = fit(model, L, Tr, Vt, y=y, cb=cb, start=th0, cons=cons, sandwich=False)
        out.append(r['th'] if r['ok'] else np.full(model.np, np.nan))
    return np.array(out)
