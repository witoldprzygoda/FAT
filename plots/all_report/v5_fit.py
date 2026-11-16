#!/usr/bin/env python3
# ============================================================================
# v5_fit.py  --  INDEPENDENT verification of WP5 (prefix v5), SPEC report_all
# ============================================================================
# Own code (numpy/scipy only); it does not include or call wp4_lib.h / wp5_lib.h / s5_blind_ana_fit.h.
# Inputs: the plain-text dump of hists_v1.root made by v5_dump.C (/hdd3/.../all_report/v5/*.txt).
# What it does:
#   1. own smoothing of the non-eta templates G_D, G_C: exp(pol3) in x = (M4 - 0.55)/0.1, scaled-Poisson likelihood with one
#      scale per histogram (sum var / sum content), bin integrals by 16-point Gauss-Legendre, minimised with scipy BFGS;
#      order rule by N_eff (>= 50: 3, 20-50: 2, 8-20: 1), F5g G_C with N_eff < 20 takes the F10 E2 parent shape;
#   2. F-main nominal (F10 E2, kappa = 1 +- 0.2, g = 1 +- 0.3, k_E2 = 1 +- 0.05) with three estimators:
#      E-IRLS : the spec estimator (weights V = max(mu,1) + k^2 var(CB) + sum coef^2 var(raw template)), solved as an
#               iterated weighted least-squares problem with scipy.optimize.least_squares (trust-region), V refreshed until
#               the fixed point is reached;
#      E-POIS : the plain Poisson likelihood (Baker-Cousins deviance + Gaussian constraints), no template/CB variance;
#      E-BBL  : Poisson likelihood with a Barlow-Beeston-lite Gaussian nuisance per bin (Conway), variance = template + CB;
#   3. errors of E-IRLS: Fisher (BB-lite), sandwich with the WP3 event-group covariance (WP3 file, not WP5's) and the CB
#      replica covariance; MC = refits with the 100 WP0 simulation replicas (own smoothing redone);
#   4. yields (raw rows x 64/63), DR, BR_A; Pearson / GLS chi2; pull runs; correlations; F-gold on F5g; s scan.
# Run:  python3 plots/all_report/v5_fit.py   (from /home/przygoda/HADES/FAT/FAT; about 2-4 min, < 2 GB)
# ============================================================================
import sys, os, math, json
import numpy as np
from scipy.optimize import least_squares, minimize
from numpy.polynomial.legendre import leggauss

D = "/hdd3/przygoda/HADES/FAT_out/all_report/v5/"
OUT = D
L = 6079.2
PT3 = 64.0 / 63.0
NB2RAW = L / PT3          # nb -> expected raw rows
BRDIR = 2.68e-4
NREP = 100
LOG = None   # opened in main() only (the module is also imported by v5_sscan.py)
def P(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    if LOG is not None: LOG.write(s + "\n"); LOG.flush()

# ---------------------------------------------------------------- binning
edges = [round(0.400 + 0.005 * i, 6) for i in range(25)] + [round(0.520 + 0.0025 * i, 6) for i in range(1, 21)] + \
        [round(0.570 + 0.005 * i, 6) for i in range(1, 27)]
edges = np.array(edges); assert len(edges) == 71
cen = 0.5 * (edges[1:] + edges[:-1])
xc = (cen - 0.55) / 0.1
glx, glw = leggauss(16)
# nodes per bin
hw = 0.5 * (edges[1:] - edges[:-1]); mid = 0.5 * (edges[1:] + edges[:-1])
NODEx = ((mid[:, None] + hw[:, None] * glx[None, :]) - 0.55) / 0.1      # 70 x 16, tilt variable at nodes
NODEw = hw[:, None] * glw[None, :]                                        # weights in GeV
F10 = ["F10_S%dE%d" % (j, e) for j in range(1, 6) for e in (1, 2)]
F5G = ["F5g_S%d" % j for j in range(1, 6)]
HN = F10 + F5G

# ---------------------------------------------------------------- read dumps
def read_rows(fn):
    d = {}
    for line in open(fn):
        p = line.split()
        d[p[0]] = np.array([float(x) for x in p[1:]])
    return d
def read_rep(fn):
    d = {}
    for line in open(fn):
        p = line.split()
        d.setdefault(p[0], []).append(np.array([float(x) for x in p[2:]]))
    return {k: np.array(v) for k, v in d.items()}

COMPS = ["C2", "C3", "C4", "C5", "Geta", "GD", "Gconv"]
SIG = ["C2", "C3", "C4", "C5", "Geta"]

def load_sim(sel):
    s = read_rows(D + "sim_%s.txt" % sel)
    T = {}
    for c in COMPS:
        for h in HN:
            T[(c, h)] = s[c + "_" + h] * NB2RAW
            T[(c, h, "v")] = (s[c + "err_" + h] * NB2RAW) ** 2
        T[(c, "win")] = s[c + "_win"] * NB2RAW
        T[(c, "tot")] = s[c + "_tot"] * NB2RAW
    return T, s
def load_simrep(sel):
    r = read_rep(D + "simrep_%s.txt" % sel)
    R = {}
    for c in COMPS:
        for h in HN + ["win", "tot"]:
            R[(c, h)] = r[c + "_" + h] * NB2RAW
    return R
def load_data(sel):
    d = read_rows(D + "data_%s.txt" % sel)
    rep = read_rep(D + "cbrep_%s.txt" % sel)
    X = {}
    for h in HN:
        X[("N", h)] = d["N_" + h]; X[("Nerr", h)] = d["Nerr_" + h]; X[("CB", h)] = d["CB_" + h]
        X[("CBrep", h)] = rep["CBrep_" + h]   # 200 x 70
    X["N_win"] = d["N_win"]; X["CB_win"] = d["CB_win"]; X["N_tot"] = d["N_tot"]; X["CB_tot"] = d["CB_tot"]
    return X

# ---------------------------------------------------------------- own smoothing
def smooth_fit(c, v, parent=None, force=None):
    """returns (order, params[4], neff, scaled) ; density exp(pol) per GeV in x"""
    C = c.sum(); V = v.sum()
    if not C > 0: return (-1, np.zeros(4), 0., False)
    neff = C * C / V if V > 0 else 0.
    if neff < 20 and parent is not None and parent[0] >= 0 and force is None:
        p = parent[1].copy()
        I = bin_integrals(parent).sum()
        p[0] += math.log(C / I)
        return (parent[0], p, neff, True)
    order = 3 if neff >= 50 else 2 if neff >= 20 else 1 if neff >= 8 else 0
    if force is not None: order = force
    s = V / C if V > 0 else 1.
    n = c / s
    npar = order + 1
    def F(p):
        pp = np.zeros(4); pp[:npar] = p
        poly = pp[0] + NODEx * (pp[1] + NODEx * (pp[2] + NODEx * pp[3]))
        return (NODEw * np.exp(poly)).sum(axis=1)
    def nll(p):
        m = np.maximum(F(p) / s, 1e-300)
        return -(np.where(n > 0, n * np.log(m), 0.) - m).sum()
    def grad(p):
        pp = np.zeros(4); pp[:npar] = p
        poly = pp[0] + NODEx * (pp[1] + NODEx * (pp[2] + NODEx * pp[3]))
        ex = NODEw * np.exp(poly)
        f = ex.sum(axis=1); m = np.maximum(f / s, 1e-300)
        g = np.zeros(npar)
        for k in range(npar):
            dk = (ex * NODEx ** k).sum(axis=1) / s
            g[k] = -((n / m - 1.) * dk).sum()
        return g
    p0 = np.zeros(npar); p0[0] = math.log(C / 0.3)
    best = None
    for meth in ("BFGS", "Nelder-Mead"):
        r = minimize(nll, p0 if best is None else best.x, jac=grad if meth == "BFGS" else None, method=meth,
                     options={"maxiter": 20000, "gtol": 1e-10} if meth == "BFGS" else {"maxiter": 20000, "xatol": 1e-10, "fatol": 1e-12})
        if best is None or r.fun < best.fun: best = r
    r = minimize(nll, best.x, jac=grad, method="BFGS", options={"maxiter": 20000, "gtol": 1e-12})
    if r.fun < best.fun: best = r
    p = np.zeros(4); p[:npar] = best.x
    return (order, p, neff, False)
def bin_integrals(sm):
    if sm[0] < 0: return np.zeros(70)
    p = sm[1]
    poly = p[0] + NODEx * (p[1] + NODEx * (p[2] + NODEx * p[3]))
    return (NODEw * np.exp(poly)).sum(axis=1)

def smooth_all(T, getc, force=None):
    """T: dict of raw templates; getc(c,h) -> content; variances from T[(c,h,'v')]. returns dict (c,h)->bin integrals, params"""
    S = {}; PR = {}
    for c in ("GD", "Gconv"):
        for h in F10:
            sm = smooth_fit(getc(c, h), T[(c, h, "v")], None, force); PR[(c, h)] = sm; S[(c, h)] = bin_integrals(sm)
        for j in range(1, 6):
            h = "F5g_S%d" % j; par = PR[(c, "F10_S%dE2" % j)]
            sm = smooth_fit(getc(c, h), T[(c, h, "v")], par, force); PR[(c, h)] = sm; S[(c, h)] = bin_integrals(sm)
    return S, PR

# ---------------------------------------------------------------- fit configuration
class Cfg:
    pass
def make_bins(X, T, S, hists, mlo=0.400, mhi=0.700, fix=None):
    """stack the selected histograms into flat arrays"""
    B = Cfg()
    sel = np.where((edges[:-1] >= mlo - 1e-9) & (edges[1:] <= mhi + 1e-9))[0]
    rows = []
    for h in hists:
        j = int(h.split("_S")[1][0]) - 1
        e = 0 if h.endswith("E1") else 1
        for i in sel: rows.append((h, j, e, i))
    B.rows = rows; B.n = len(rows)
    B.h = [r[0] for r in rows]; B.j = np.array([r[1] for r in rows]); B.e = np.array([r[2] for r in rows]); B.i = np.array([r[3] for r in rows])
    B.y = np.array([X[("N", r[0])][r[3]] for r in rows])
    B.cb = np.array([X[("CB", r[0])][r[3]] for r in rows])
    B.vcb = np.array([np.var(X[("CBrep", r[0])][:, r[3]], ddof=1) for r in rows])
    B.egv = np.array([X[("Nerr", r[0])][r[3]] ** 2 for r in rows])
    B.x = xc[B.i]
    B.gd = np.array([S[("GD", r[0])][r[3]] for r in rows]); B.gc = np.array([S[("Gconv", r[0])][r[3]] for r in rows])
    B.s = {c: np.array([T[(c, r[0])][r[3]] for r in rows]) for c in SIG}
    B.vs = {c: np.array([T[(c, r[0], "v")][r[3]] for r in rows]) for c in SIG}
    B.fix = np.zeros(B.n) if fix is None else fix
    return B
def swap_templates(B, Trep, Srep, Tvar):
    """replica templates (variances stay those of the full sample)"""
    B2 = Cfg(); B2.__dict__.update(B.__dict__)
    B2.gd = np.array([Srep[("GD", r[0])][r[3]] for r in B.rows]); B2.gc = np.array([Srep[("Gconv", r[0])][r[3]] for r in B.rows])
    B2.s = {c: np.array([Trep[(c, r[0])][r[3]] for r in B.rows]) for c in SIG}
    return B2

# parameter layouts
MAIN = ["a_eta", "kappa", "rho_d", "g", "c_1", "c_2", "c_3", "c_4", "c_5", "b", "kappa_B", "k_E2"]
def model_main(th, B, quad=False):
    a, ka, rho, g = th[0], th[1], th[2], th[3]
    c = th[4:9][B.j]; b = th[9]; kB = th[10]; k = th[11]
    tl = 1. + b * B.x
    bg = B.gd + kB * B.gc
    co = {"C2": a, "C3": a * ka, "C4": a * ka, "C5": a * rho, "Geta": a * g}
    mu = k * B.cb + c * tl * bg + sum(co[t] * B.s[t] for t in SIG) + B.fix
    vs = sum(co[t] ** 2 * B.vs[t] for t in SIG)
    return mu, vs, k
CONS_MAIN = [(1, 1.0, 0.2), (3, 1.0, 0.3), (11, 1.0, 0.05)]

GOLD = ["alpha_2", "alpha_5", "g", "c_1", "c_2", "c_3", "c_4", "c_5", "b", "k_E2"]
def make_model_gold(kB):
    def model_gold(th, B):
        a2, a5, g = th[0], th[1], th[2]
        c = th[3:8][B.j]; b = th[8]; k = th[9]
        bg = B.gd + kB * B.gc
        co = {"C2": a2, "C5": a5, "Geta": a2 * g}
        mu = k * B.cb + c * (1. + b * B.x) * bg + sum(co[t] * B.s[t] for t in co) + B.fix
        vs = sum(co[t] ** 2 * B.vs[t] for t in co)
        return mu, vs, k
    return model_gold
CONS_GOLD = [(2, 1.0, 0.3), (9, 1.0, 0.05)]

def V_of(B, th, model):
    mu, vs, k = model(th, B)
    return np.maximum(mu, 1.0) + k * k * B.vcb + vs, mu

# ---------------------------------------------------------------- estimator E-IRLS (fixed point via trust-region LSQ)
def fit_irls(B, model, cons, th0, tol=1e-10, maxouter=200):
    th = np.array(th0, float)
    V, _ = V_of(B, th, model)
    for outer in range(maxouter):
        sw = 1. / np.sqrt(V)
        def res(t):
            mu = model(t, B)[0]
            r = (B.y - mu) * sw
            pc = [(t[k] - v) / s for (k, v, s) in cons]
            return np.concatenate([r, pc])
        r = least_squares(res, th, method="trf", x_scale="jac", xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=5000)
        thn = r.x
        Vn, _ = V_of(B, thn, model)
        d = np.max(np.abs(thn - th) / np.maximum(1e-6, np.abs(th)))
        th = thn; V = Vn
        if d < tol and outer > 0: break
    return th, outer

def jac(B, model, th, eps=1e-6):
    mu0 = model(th, B)[0]
    J = np.zeros((B.n, len(th)))
    for k in range(len(th)):
        h = eps * max(1., abs(th[k]))
        tp = th.copy(); tm = th.copy(); tp[k] += h; tm[k] -= h
        J[:, k] = (model(tp, B)[0] - model(tm, B)[0]) / (2 * h)
    return J, mu0

def errors(B, model, cons, th, egfull=None, cbcov=None):
    V, mu = V_of(B, th, model)
    J, _ = jac(B, model, th)
    Pm = np.zeros((len(th), len(th)))
    for (k, v, s) in cons: Pm[k, k] += 1. / s ** 2
    A = J.T @ (J / V[:, None]) + Pm
    C = np.linalg.inv(A)
    U = J / V[:, None]
    # event-group covariance of N scaled to the model (sqrt(mu_a mu_b / (y_a y_b)))
    pos = B.y > 0
    sf = np.where(pos, np.sqrt(np.maximum(mu, 1e-9) / np.where(pos, B.y, 1.)), 0.)
    if egfull is not None: Sn = egfull * np.outer(sf, sf)
    else: Sn = np.diag(B.egv * sf * sf)
    Sn = Sn + np.diag(np.where(pos, 0., np.maximum(mu, 0.)))   # empty bins: Poisson variance of the model
    k = model(th, B)[2]
    Scb = (cbcov if cbcov is not None else np.diag(B.vcb)) * k * k
    Mn = U.T @ Sn @ U; Mc = U.T @ Scb @ U
    Cn = C @ Mn @ C; Ccb = C @ Mc @ C
    Ccons = []
    for (kk, v, s) in cons:
        Mp = np.zeros_like(Pm); Mp[kk, kk] = 1. / s ** 2
        Ccons.append(C @ Mp @ C)
    Cstat = Cn + Ccb + sum(Ccons)
    chi2P = ((B.y - mu) ** 2 / V).sum() + sum(((th[kk] - v) / s) ** 2 for (kk, v, s) in cons)
    m = np.maximum(mu, 1e-9)
    dev = 2 * (m - B.y + np.where(B.y > 0, B.y * np.log(np.where(B.y > 0, B.y, 1) / m), 0.)).sum() + sum(((th[kk] - v) / s) ** 2 for (kk, v, s) in cons)
    return dict(C=C, Cn=Cn, Ccb=Ccb, Ccons=Ccons, Cstat=Cstat, chi2P=chi2P, dev=dev, mu=mu, V=V, J=J)

# ---------------------------------------------------------------- estimator E-POIS and E-BBL
def fit_pois(B, model, cons, th0, bbl=False):
    # E-BBL: the nuisance variance (template + CB) is FIXED at th0 (the E-IRLS solution), as in HistFactory/Conway BB-lite;
    # a theta-dependent variance without its normalisation term would let the fit inflate the coefficients (runaway).
    mu0, vs0, k0 = model(np.array(th0, float), B)
    s2fix = vs0 + k0 * k0 * B.vcb
    def nll(t):
        mu, vs, k = model(t, B)
        if bbl:
            s2 = s2fix
            m = np.maximum(mu, 1e-9)
            # Conway BB-lite: profile delta
            bq = m + s2
            disc = np.maximum(bq * bq - 4 * s2 * (m - B.y), 0.)
            dl = 0.5 * (-bq + np.sqrt(disc))
            mt = np.maximum(m + dl, 1e-9)
            pen = np.where(s2 > 0, dl * dl / np.maximum(s2, 1e-300), 0.)
            val = 2 * (mt - B.y + np.where(B.y > 0, B.y * np.log(np.where(B.y > 0, B.y, 1) / mt), 0.)).sum() + pen.sum()
        else:
            m = np.maximum(mu, 1e-9)
            val = 2 * (m - B.y + np.where(B.y > 0, B.y * np.log(np.where(B.y > 0, B.y, 1) / m), 0.)).sum()
        return val + sum(((t[kk] - v) / s) ** 2 for (kk, v, s) in cons)
    th = np.array(th0, float)
    best = None
    for it in range(4):
        r = minimize(nll, th, method="BFGS", options={"maxiter": 50000, "gtol": 1e-7})
        r2 = minimize(nll, r.x, method="Nelder-Mead", options={"maxiter": 40000, "xatol": 1e-9, "fatol": 1e-9, "adaptive": True})
        th = r2.x if r2.fun < r.fun else r.x
        f = min(r.fun, r2.fun)
        if best is not None and abs(best - f) < 1e-7: break
        best = f
    # numerical Hessian -> covariance (2 x NLL scale: cov = 2 H^-1)
    n = len(th); H = np.zeros((n, n)); f0 = nll(th)
    hs = [1e-4 * max(1., abs(t)) for t in th]
    for a_ in range(n):
        for b_ in range(a_, n):
            ta = th.copy(); ta[a_] += hs[a_]; ta[b_] += hs[b_]; fpp = nll(ta)
            ta = th.copy(); ta[a_] += hs[a_]; ta[b_] -= hs[b_]; fpm = nll(ta)
            ta = th.copy(); ta[a_] -= hs[a_]; ta[b_] += hs[b_]; fmp = nll(ta)
            ta = th.copy(); ta[a_] -= hs[a_]; ta[b_] -= hs[b_]; fmm = nll(ta)
            H[a_, b_] = H[b_, a_] = (fpp - fpm - fmp + fmm) / (4 * hs[a_] * hs[b_])
    try: C = 2 * np.linalg.inv(H)
    except Exception: C = np.full((n, n), np.nan)
    return th, f0, C

# ---------------------------------------------------------------- derived quantities
def W(T, c, w):  # window sums (raw rows): win histogram bins 1 = B, 2 = P, 3 = SBlo, 4 = SBhi; tot
    if w == "tot": return T[(c, "tot")][0]
    return T[(c, "win")][{"B": 0, "P": 1}[w]]
def derived_main(th, T):
    a, ka, rho, g = th[0], th[1], th[2], th[3]
    d = {}
    d["a_eta"] = a; d["kappa"] = ka; d["rho_d"] = rho; d["g"] = g; d["kappa_B"] = th[10]
    d["Y2_B"] = PT3 * a * W(T, "C2", "B"); d["Y3_B"] = PT3 * a * ka * W(T, "C3", "B")
    d["Y_bump"] = d["Y2_B"] + d["Y3_B"]
    d["Y4_P"] = PT3 * a * ka * W(T, "C4", "P"); d["Y5_P"] = PT3 * a * rho * W(T, "C5", "P")
    d["Y_peak"] = d["Y4_P"] + d["Y5_P"]
    for c, f in (("C2", 1.), ("C3", ka), ("C4", ka), ("C5", rho)):
        d["Y%s_tot" % c[1]] = PT3 * a * f * W(T, c, "tot")
    sim = (W(T, "C4", "P") + W(T, "C5", "P")) / (W(T, "C2", "B") + W(T, "C3", "B"))
    d["DR"] = (d["Y_peak"] / d["Y_bump"]) / sim
    d["BR_A"] = BRDIR * rho
    d["sigvis_bump"] = d["Y_bump"] / L; d["sigvis_peak"] = d["Y_peak"] / L
    return d
def grad_der(fn, th, T, key, C):
    g = np.zeros(len(th))
    for k in range(len(th)):
        h = max(1e-7, 1e-4 * math.sqrt(max(C[k, k], 0.)))
        if C[k, k] <= 0: continue
        tp = th.copy(); tm = th.copy(); tp[k] += h; tm[k] -= h
        g[k] = (fn(tp, T)[key] - fn(tm, T)[key]) / (2 * h)
    return g

# ---------------------------------------------------------------- main
def main():
    global LOG, NREP
    NREP = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    LOG = open(OUT + "v5_fit.log", "w")
    P("v5_fit.py: independent verification of WP5 F-main / F-gold;", "NREP =", NREP)
    T, sraw = load_sim("S0")
    X = load_data("S0")
    P("units: NB2RAW = %.4f (L x 63/64), PT3 = %.6f" % (NB2RAW, PT3))
    P("data S0: N_win(P) = %.0f, CB_win(P) = %.2f, CB/N = %.4f; N_win(B) = %.0f" % (X["N_win"][1], X["CB_win"][1], X["CB_win"][1] / X["N_win"][1], X["N_win"][0]))
    # consistency: sum of F10 cells in B equals the win histogram (data N and sim C2)
    bB = np.where((edges[:-1] >= 0.455 - 1e-9) & (edges[1:] <= 0.5225 + 1e-9))[0]
    nB = sum(X[("N", h)][bB].sum() for h in F10)
    c2B = sum(T[("C2", h)][bB].sum() for h in F10)
    P("check: data N in B from F10 cells %.1f vs win %.1f; sim C2 in B from F10 %.4f vs win %.4f (raw rows)" % (nB, X["N_win"][0], c2B, W(T, "C2", "B")))
    # ---- own smoothing
    S, PR = smooth_all(T, lambda c, h: T[(c, h)])
    P("\n# own smoothing (order, p0..p3, neff, scaled) vs WP5 'smooth' lines")
    wp5s = {}
    for line in open("/hdd3/przygoda/HADES/FAT_out/all_report/wp5/results.txt"):
        if line.startswith("smooth "):
            p = line.split()
            wp5s[(p[1], "GD")] = [float(x) for x in p[5:9]]; wp5s[(p[1], "Gconv")] = [float(x) for x in p[16:20]]
    maxrel = 0
    for h in HN:
        for c in ("GD", "Gconv"):
            sm = PR[(c, h)]
            mine = bin_integrals(sm)
            pw = np.array(wp5s[(h, c)])
            o = sm[0]
            pw_poly = pw[0] + NODEx * (pw[1] + NODEx * (pw[2] + NODEx * pw[3]))
            theirs = (NODEw * np.exp(pw_poly)).sum(axis=1)
            rel = np.max(np.abs(mine - theirs) / np.maximum(theirs, 1e-3 * theirs.max()))
            maxrel = max(maxrel, rel)
            P("  %-10s %-5s order %d p = %s neff %.1f scaled %d | max rel. diff of bin integrals vs WP5 %.2e" % (h, c, o, " ".join("%.5f" % v for v in sm[1]), sm[2], sm[3], rel))
    P("  max relative difference of smoothed bin integrals (own vs WP5) = %.2e" % maxrel)

    # ---- F-main nominal (E2)
    hE2 = ["F10_S%dE2" % j for j in range(1, 6)]
    B = make_bins(X, T, S, hE2)
    P("\n# F-main nominal: %d bins (F10 E2)" % B.n)
    th0 = np.array([0.4, 1.0, 1.0, 1.0, 1, 1, 1, 1, 1, -0.07, 2.0, 1.0])
    thI, nout = fit_irls(B, model_main, CONS_MAIN, th0)
    # full event-group covariance (WP3 file) restricted to the E2 cells, CB replica covariance
    eg3 = np.loadtxt(D + "egcov_wp3_S0_F10.txt"); eg5 = np.loadtxt(D + "egcov_wp5_S0_F10.txt")
    P("event-group covariance: WP3 vs WP5 files max |diff| = %.3e (max element %.1f)" % (np.max(np.abs(eg3 - eg5)), np.max(np.abs(eg3))))
    diagN = np.concatenate([X[("Nerr", h)] ** 2 for h in F10])
    P("event-group covariance diagonal vs hists_v1 N errors^2: max rel diff %.2e" % (np.max(np.abs(np.diag(eg3) - diagN) / np.maximum(diagN, 1))))
    hidx = {h: k for k, h in enumerate(F10)}
    cells = np.array([hidx[B.h[q]] * 70 + B.i[q] for q in range(B.n)])
    egE2 = eg3[np.ix_(cells, cells)]
    cbrepE2 = np.concatenate([X[("CBrep", h)] for h in hE2], axis=1)   # 200 x 350
    cbcov = np.cov(cbrepE2, rowvar=False, ddof=1)
    allN = sum(X[("N", h)].sum() for h in F10); allV = eg3.sum()
    P("variance of the F10 total / N = %.4f" % (allV / allN))
    E = errors(B, model_main, CONS_MAIN, thI, egE2, cbcov)
    errN = np.sqrt(np.diag(E["Cn"])); errCB = np.sqrt(np.diag(E["Ccb"])); errS = np.sqrt(np.diag(E["Cstat"])); errF = np.sqrt(np.diag(E["C"]))
    # WP5 reference values
    ref = {}
    for line in open("/hdd3/przygoda/HADES/FAT_out/all_report/wp5/results.txt"):
        p = line.split()
        if len(p) > 2 and p[1] == "FMAIN" and p[0] in ("values", "err_N", "err_CB", "err_stat", "err_mc", "mc_shift"):
            ref[p[0]] = np.array([float(x) for x in p[2:]])
        if len(p) > 2 and p[0] == "derived" and p[1] == "FMAIN":
            ref["d_" + p[2]] = dict(v=float(p[3]), N=float(p[5]), CB=float(p[7]), stat=float(p[11]), mc=float(p[13]), shift=float(p[17]))
        if len(p) > 2 and p[0] == "gof" and p[1] == "FMAIN":
            ref["gof"] = {p[k]: float(p[k + 1]) for k in range(2, len(p) - 1, 2)}
    P("E-IRLS converged after %d outer iterations; Pearson chi2 = %.2f / %d; BC deviance = %.2f" % (nout, E["chi2P"], B.n - 12 + 3, E["dev"]))
    P("WP5: chi2P %.3f chi2BC %.3f" % (ref["gof"]["chi2P"], ref["gof"]["chi2BC"]))
    P("%-8s %12s %12s | %9s %9s | %9s %9s | %9s %9s | %9s" % ("param", "v5 IRLS", "WP5", "errN v5", "WP5", "errCB v5", "WP5", "stat v5", "WP5", "Fisher"))
    for k, nm in enumerate(MAIN):
        P("%-8s %12.6f %12.6f | %9.5f %9.5f | %9.5f %9.5f | %9.5f %9.5f | %9.5f" % (nm, thI[k], ref["values"][k], errN[k], ref["err_N"][k], errCB[k], ref["err_CB"][k], errS[k], ref["err_stat"][k], errF[k]))
    # correlations
    Cst = E["Cstat"]
    def corr(C, a, b): return C[a, b] / math.sqrt(C[a, a] * C[b, b])
    P("corr stat: (a_eta, b) = %.3f, (rho_d, kappa) = %.3f, (a_eta, g) = %.3f, (c_1, kappa_B) = %.3f" % (corr(Cst, 0, 9), corr(Cst, 2, 1), corr(Cst, 0, 3), corr(Cst, 4, 10)))
    dI = derived_main(thI, T)
    # ---- other estimators
    thP, fP, CP = fit_pois(B, model_main, CONS_MAIN, thI, bbl=False)
    thB, fB, CB_ = fit_pois(B, model_main, CONS_MAIN, thI, bbl=True)
    dP = derived_main(thP, T); dB = derived_main(thB, T)
    P("\n# estimator comparison (F-main nominal, E2): value (+- Hessian error for the likelihood fits)")
    P("%-8s %10s %10s %16s %16s" % ("param", "WP5", "E-IRLS", "E-POIS", "E-BBL"))
    for k, nm in enumerate(MAIN):
        P("%-8s %10.4f %10.4f %9.4f+-%5.4f %9.4f+-%5.4f" % (nm, ref["values"][k], thI[k], thP[k], math.sqrt(max(CP[k, k], 0)), thB[k], math.sqrt(max(CB_[k, k], 0))))
    P("E-POIS deviance %.2f, E-BBL profiled deviance %.2f" % (fP, fB))
    # ---- MC replicas (own smoothing redone per replica)
    R = load_simrep("S0")
    REPS = []
    for r in range(NREP):
        Tr = {}
        for c in COMPS:
            for h in HN: Tr[(c, h)] = R[(c, h)][r]; Tr[(c, h, "v")] = T[(c, h, "v")]
            Tr[(c, "win")] = R[(c, "win")][r]; Tr[(c, "tot")] = R[(c, "tot")][r]
        Sr, _ = smooth_all(Tr, lambda c, h: Tr[(c, h)])
        REPS.append((Tr, Sr))
    P("replica templates smoothed: %d" % len(REPS))
    rep = []; repd = []
    for r in range(NREP):
        Tr, Sr = REPS[r]
        Br = swap_templates(B, Tr, Sr, T)
        thr, _ = fit_irls(Br, model_main, CONS_MAIN, thI, tol=1e-8)
        rep.append(thr); repd.append(derived_main(thr, Tr))
        if r % 20 == 0: P("  replica %d: a_eta %.4f" % (r, thr[0]))
    rep = np.array(rep)
    mcE = rep.std(axis=0, ddof=1); mcS = rep.mean(axis=0) - thI
    P("\n# MC (replica) errors and shifts: v5 vs WP5")
    for k, nm in enumerate(MAIN):
        P("%-8s mc %9.5f (WP5 %9.5f)  shift %+9.5f (WP5 %+9.5f)  median %.4f" % (nm, mcE[k], ref["err_mc"][k], mcS[k], ref["mc_shift"][k], np.median(rep[:, k])))
    # derived quantities with errors
    keys = ["a_eta", "kappa", "rho_d", "g", "kappa_B", "Y2_B", "Y3_B", "Y_bump", "Y4_P", "Y5_P", "Y_peak", "Y2_tot", "Y3_tot", "Y4_tot", "Y5_tot", "DR", "BR_A", "sigvis_bump", "sigvis_peak"]
    P("\n# derived quantities: v5 value / N / CB / stat / MC vs WP5; E-POIS and E-BBL values")
    out = {}
    for key in keys:
        g = grad_der(derived_main, thI, T, key, E["C"])
        eN = math.sqrt(max(g @ E["Cn"] @ g, 0)); eCB = math.sqrt(max(g @ E["Ccb"] @ g, 0)); eS = math.sqrt(max(g @ E["Cstat"] @ g, 0))
        vals = np.array([d[key] for d in repd]); eMC = vals.std(ddof=1); sh = vals.mean() - dI[key]
        rf = ref.get("d_" + key, None)
        P("%-12s v5 %12.6g N %9.4g CB %9.4g stat %9.4g MC %9.4g shift %+9.4g | WP5 %12.6g N %9.4g CB %9.4g stat %9.4g MC %9.4g shift %+9.4g | POIS %12.6g BBL %12.6g" % (
            key, dI[key], eN, eCB, eS, eMC, sh, rf["v"] if rf else float("nan"), rf["N"] if rf else 0, rf["CB"] if rf else 0, rf["stat"] if rf else 0, rf["mc"] if rf else 0, rf["shift"] if rf else 0, dP[key], dB[key]))
        out[key] = dict(v=dI[key], N=eN, CB=eCB, stat=eS, mc=eMC, shift=sh, pois=dP[key], bbl=dB[key])
    # ---- GLS chi2 with MC covariance of the prediction at thI
    mus = []
    for r in range(NREP):
        Tr, Sr = REPS[r]
        mus.append(model_main(thI, swap_templates(B, Tr, Sr, T))[0])
    mus = np.array(mus); Cmc = np.cov(mus, rowvar=False, ddof=1)
    mu = E["mu"]; pos = B.y > 0; sf = np.where(pos, np.sqrt(np.maximum(mu, 1e-9) / np.where(pos, B.y, 1.)), 0.)
    Cd = egE2 * np.outer(sf, sf) + np.diag(np.where(pos, 0., np.maximum(mu, 1.))) + thI[11] ** 2 * cbcov
    r_ = B.y - mu
    pen = sum(((thI[kk] - v) / s) ** 2 for (kk, v, s) in CONS_MAIN)
    chiG = r_ @ np.linalg.solve(Cd + Cmc, r_) + pen
    chiGd = r_ @ np.linalg.solve(np.diag(np.diag(Cd)) + Cmc, r_) + pen
    P("\nGLS chi2 at the fit (data full cov + MC cov) = %.1f / 341 = %.3f (WP5 %.1f); diag data cov: %.1f (WP5 %.1f)" % (chiG, chiG / 341, ref["gof"]["chi2GLS"], chiGd, ref["gof"]["chi2GLSdiag"]))
    # Pearson by region
    pull = (B.y - mu) / np.sqrt(E["V"])
    for (lo, hi, nm) in ((0.400, 0.455, "0.400-0.455"), (0.455, 0.5225, "B"), (0.5225, 0.565, "0.5225-0.565"), (0.565, 0.575, "0.565-0.575"), (0.575, 0.700, "CR-U 0.575-0.700")):
        m = (edges[B.i] >= lo - 1e-9) & (edges[B.i + 1] <= hi + 1e-9)
        P("  Pearson in %-18s %7.1f for %3d bins" % (nm, (pull[m] ** 2).sum(), m.sum()))
    # runs in 0.52-0.56
    P("pull runs in 0.52-0.56 (same sign), per slice:")
    for j in range(5):
        m = (B.j == j) & (edges[B.i] >= 0.52 - 1e-9) & (edges[B.i + 1] <= 0.56 + 1e-9)
        p = pull[m]; L_ = 0; cur = 0; sg = 0
        for v in p:
            s_ = 1 if v > 0 else -1
            cur = cur + 1 if s_ == sg else 1; sg = s_; L_ = max(L_, cur)
        P("  S%d: n %d longest run %d; signs %s" % (j + 1, m.sum(), L_, "".join("+" if v > 0 else "-" for v in p)))
    # ---- F-gold
    hG = F5G
    # fixed term: a kappa (C3 + C4) from F-main in F5g
    fixG = np.concatenate([thI[0] * thI[1] * (T[("C3", h)] + T[("C4", h)]) for h in hG])
    BG = make_bins(X, T, S, hG, fix=fixG)
    BG.e = np.ones(BG.n, int)
    mg = make_model_gold(thI[10])
    thG0 = np.array([0.4, 0.4, 1.0, 1, 1, 1, 1, 1, -0.07, 1.0])
    thG, _ = fit_irls(BG, mg, CONS_GOLD, thG0)
    eg5g = np.loadtxt(D + "egcov_wp5_S0_F5g.txt")
    cbg = np.cov(np.concatenate([X[("CBrep", h)] for h in hG], axis=1), rowvar=False, ddof=1)
    EG = errors(BG, mg, CONS_GOLD, thG, eg5g, cbg)
    brg = BRDIR * thG[1] / thG[0]
    g_ = np.zeros(len(thG)); g_[0] = -brg / thG[0]; g_[1] = BRDIR / thG[0]
    eBR = math.sqrt(g_ @ EG["Cstat"] @ g_)
    thGp, fGp, CGp = fit_pois(BG, mg, CONS_GOLD, thG, bbl=False)
    thGb, fGb, CGb = fit_pois(BG, mg, CONS_GOLD, thG, bbl=True)
    P("\n# F-gold (F5g, M(ee) >= 20 MeV; kappa_B = %.4f and a kappa (C3 + C4) fixed from F-main)" % thI[10])
    for k, nm in enumerate(GOLD):
        P("  %-8s IRLS %.5f +- %.5f (stat)  POIS %.5f  BBL %.5f" % (nm, thG[k], math.sqrt(EG["Cstat"][k, k]), thGp[k], thGb[k]))
    P("  BR_G = %.4e +- %.2e (stat); POIS %.4e; BBL %.4e; Pearson chi2 %.1f/342" % (brg, eBR, BRDIR * thGp[1] / thGp[0], BRDIR * thGb[1] / thGb[0], EG["chi2P"]))
    # MC of F-gold
    brs = []
    for r in range(NREP):
        Tr, Sr = REPS[r]
        thr = rep[r]
        fixr = np.concatenate([thr[0] * thr[1] * (Tr[("C3", h)] + Tr[("C4", h)]) for h in hG])
        Bgr = swap_templates(BG, Tr, Sr, T); Bgr.fix = fixr
        tg, _ = fit_irls(Bgr, make_model_gold(thr[10]), CONS_GOLD, thG, tol=1e-8)
        brs.append(BRDIR * tg[1] / tg[0])
    brs = np.array(brs); q = np.percentile(brs, [15.87, 50, 84.13])
    P("  BR_G MC: rms %.3e, mc68 %.3e, mean shift %+.3e, median %.4e" % (brs.std(ddof=1), 0.5 * (q[2] - q[0]), brs.mean() - brg, q[1]))
    # ---- E-POIS MC spread (same replicas), to compare estimators at equal footing
    repP = []
    for r in range(min(NREP, 30)):
        Tr, Sr = REPS[r]
        tp, _, _ = fit_pois(swap_templates(B, Tr, Sr, T), model_main, CONS_MAIN, thP, bbl=False)
        repP.append(derived_main(tp, Tr))
    for key in ("a_eta", "Y_bump", "Y_peak", "rho_d"):
        v = np.array([d[key] for d in repP])
        P("E-POIS replicas (%d): %-7s rms %.4g, mean shift %+.4g (nominal %.5g)" % (len(repP), key, v.std(ddof=1), v.mean() - dP[key], dP[key]))
    # ---- per-slice a_eta (kappa, rho_d, g common)
    def model_slice(th, Bx):
        a = th[0:5][Bx.j]; ka, rho, g = th[5], th[6], th[7]
        c = th[8:13][Bx.j]; b = th[13]; kB = th[14]; k = th[15]
        bg = Bx.gd + kB * Bx.gc
        co = {"C2": a, "C3": a * ka, "C4": a * ka, "C5": a * rho, "Geta": a * g}
        mu = k * Bx.cb + c * (1. + b * Bx.x) * bg + sum(co[t] * Bx.s[t] for t in SIG)
        vs = sum(co[t] ** 2 * Bx.vs[t] for t in SIG)
        return mu, vs, k
    consS = [(5, 1.0, 0.2), (7, 1.0, 0.3), (15, 1.0, 0.05)]
    thS0 = np.concatenate([[thI[0]] * 5, thI[1:]])
    thS, _ = fit_irls(B, model_slice, consS, thS0)
    ES = errors(B, model_slice, consS, thS, egE2, cbcov)
    eS = np.sqrt(np.diag(ES["Cstat"]))[:5]
    w = 1 / eS ** 2; m_ = (w * thS[:5]).sum() / w.sum(); ch = (w * (thS[:5] - m_) ** 2).sum()
    from scipy.stats import chi2 as CHI2
    P("per-slice a_eta: %s; stat %s; chi2(stat, diag) %.2f/4 p = %.3f" % (" ".join("%.3f" % v for v in thS[:5]), " ".join("%.3f" % v for v in eS), ch, CHI2.sf(ch, 4)))
    Cs = ES["Cstat"][:5, :5]
    dvec = thS[:5] - m_
    # full-covariance consistency (stat): chi2 = min over common mean
    Ci = np.linalg.inv(Cs); one = np.ones(5); mm = (one @ Ci @ thS[:5]) / (one @ Ci @ one); dd = thS[:5] - mm
    chf = dd @ Ci @ dd
    P("per-slice a_eta full stat covariance: chi2 %.2f/4 p = %.3f (WP5: 7.57/4, p = 0.109)" % (chf, CHI2.sf(chf, 4)))
    # ---- s scan with fixed weights (V of the nominal fit)
    sv = read_rows(D + "sim_S0_svar.txt")
    Vfix = E["V"].copy()
    def fit_fixedV(Bx, th0):
        def res(t):
            return np.concatenate([(Bx.y - model_main(t, Bx)[0]) / np.sqrt(Vfix), [(t[kk] - v) / s_ for (kk, v, s_) in CONS_MAIN]])
        r = least_squares(res, th0, method="trf", x_scale="jac", xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=5000)
        return r.x, 2 * r.cost
    chis = {}
    for var in ("sDn", "nom", "sUp", "sK0", "r089", "r110"):
        if var == "nom": Tv = T
        else:
            Tv = {}
            for c in COMPS:
                for h in HN:
                    Tv[(c, h)] = sv[var + ":" + c + "_" + h] * NB2RAW; Tv[(c, h, "v")] = (sv[var + ":" + c + "err_" + h] * NB2RAW) ** 2
        Sv, _ = smooth_all(Tv, lambda c, h: Tv[(c, h)])
        Bv = make_bins(X, Tv, Sv, hE2)
        tv, cv = fit_fixedV(Bv, thI)
        chis[var] = cv
        P("s-scan %-5s fixed-weight chi2 %.2f  a_eta %.4f rho_d %.4f" % (var, cv, tv[0], tv[2]))
    sg = np.array([0.99525, 0.99615, 0.99705]); cg = np.array([chis["sDn"], chis["nom"], chis["sUp"]])
    A_, B_, C_ = np.polyfit(sg, cg, 2)
    shat = -B_ / (2 * A_); ssig = 1 / math.sqrt(A_)
    P("free s (whole fit, fixed weights): s_hat = %.5f +- %.5f -> %.2f sigma from 0.99615 (WP5 0.99653 +- 0.00019)" % (shat, ssig, (shat - 0.99615) / ssig))
    json.dump(dict(main_irls=list(thI), main_pois=list(thP), main_bbl=list(thB), derived=out, gold=list(thG), BR_G=brg, BR_G_stat=eBR,
                   BR_G_mc_rms=float(brs.std(ddof=1)), BR_G_mc68=float(0.5 * (q[2] - q[0])), chi2P=E["chi2P"], chi2GLS=chiG,
                   gold_pois=list(thGp), gold_bbl=list(thGb), mc=list(mcE), mcshift=list(mcS)), open(OUT + "v5_fit.json", "w"), indent=1)
    P("done")

if __name__ == "__main__":
    main()
