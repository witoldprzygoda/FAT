#!/usr/bin/env python3
# ============================================================================
# prevstyle_ana.py  --  package prevstyle of report_all: data and SMASH after the typical eta cuts (E1) drawn in the style
#                       of the earlier analysis (collaboration-meeting slides); figures, fits and numbers (prefix prevstyle)
# ============================================================================
# Written 2026-11-16. New file; nothing else is changed.
# Inputs : /hdd3/przygoda/HADES/FAT_out/all_report/prevstyle/prevstyle_data.json, prevstyle_sim.json (prevstyle_fill.C)
#          /hdd3/przygoda/HADES/FAT_out/all_report/e1/e1_data.json, e1_sim.json (package e1, read only: reproduction check)
# Outputs: docs/latex/figs_all/prevstyle_m4_E1.pdf/.png        (F1: M4 spectra, (a) SMASH scaled to the data integral,
#                                                               (b) SMASH at its cross sections with the measured scales)
#          docs/latex/figs_all/prevstyle_globalfit_E1.pdf/.png (F2: three-scale global fit, raw and smoothed templates)
#          docs/latex/figs_all/prevstyle_mee_E1.pdf/.png       (F3: M(ee) in the eta peak window against the side bands)
#          /hdd3/przygoda/HADES/FAT_out/all_report/prevstyle/prevstyle_numbers.json, prevstyle_table.txt, prevstyle_ana.log
# Conventions (package e1 / Part III):
#   selection pre + E1 (alpha_lab < 50 deg, alpha_eta > 140 deg, M(pipi) < 0.46 GeV/c2, OA(e+e-) > 2 deg); two variants:
#   "all" = isBest >= 0 (all stored combinations; NOMINAL, drawn first) and "best" = isBest == 1 (best combination only, as in
#   the earlier analysis; comparison only; data AND SMASH with isBest == 1);
#   data in raw counts: N, CB = 2 sqrt(N++ N--) from the like-sign data with the same cuts (k = 1; for isBest == 1 the naive
#   like-sign formula of the earlier analysis), N - CB; errors: N event groups (+) CB replicas;
#   SMASH in expected raw counts = nb x L x 63/64 (L = 6079.2 nb^-1), sigma x PDG BR x V3, mode A, M4' = 0.99615 M4,
#   M(pipi)' = 0.9952 M(pipi); C2 with the per-count QED gamma* weights (c2recon); errors: event groups of the simulation.
#   classes: C2 eta -> pi+pi-pi0(e+e-gamma), C3 eta -> pi+pi-pi0(gamma gamma), one photon converted, C4 eta -> pi+pi-gamma
#   with the photon converted, C5 direct eta -> pi+pi-e+e-; "rest" = other eta (G_eta) + non-eta (G_D + G_C).
#   measured scales: bump C2, C3 x f_eta = 0.426; peak C4, C5 x a_P = 0.418; G_eta x a_eta g = 0.656; G_D + G_C x c_b, c_b
#   from the data in 0.575-0.70 GeV/c2 (package e1: 1.450 with all eta x 0.40).
# Run: python3 -I plots/all_report/prevstyle_ana.py   (about 1 min)
# ============================================================================
import json
import math
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.optimize import minimize
from scipy.special import voigt_profile

OUT = "/hdd3/przygoda/HADES/FAT_out/all_report/prevstyle/"
E1DIR = "/hdd3/przygoda/HADES/FAT_out/all_report/e1/"
FIG = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_all/"
L = 6079.2
PT3 = 64. / 63.
NB2RAW = L / PT3
F_BUMP, A_P, G_ETA = 0.426, 0.418, 0.656
C_B_E1 = 1.450
FLEX = 0.46
SELS = ["all", "best"]          # nominal first
SELLAB = {"all": "isBest $\\geq$ 0: all stored combinations\n(nominal)",
          "best": "isBest == 1: best combination only,\nas in the earlier analysis (comparison)"}
SELSHORT = {"all": "isBest >= 0 (nominal)", "best": "isBest == 1 (earlier analysis, comparison)"}
PAL = dict(N="black", CB="red", sig="#1f3fbf", tot="black", C2="#ff6a00", C3="#d400d4", C4="#009900", C5="#00a5c8",
           rest="#8a8a8a", fit="#d62728", pre="#555555", X17="#cc0077", sb="#e07b00")
CLAB = {"C5": r"$\eta\to\pi^+\pi^-e^+e^-$",
        "C4": r"$\eta\to\pi^+\pi^-\gamma(e^+e^-)$",
        "C2": r"$\eta\to\pi^+\pi^-\pi^0(e^+e^-\gamma)$",
        "C3": r"$\eta\to\pi^+\pi^-\pi^0(\gamma\gamma)$, $\gamma$ converted",
        "rest": r"rest: non-$\eta$ + other $\eta$ decays"}
plt.rcParams.update({"font.size": 13, "axes.labelsize": 15, "axes.titlesize": 13, "legend.fontsize": 10.5, "xtick.labelsize": 12,
                     "ytick.labelsize": 12, "pdf.fonttype": 42})
LOG = open(OUT + "prevstyle_ana.log", "w")


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    LOG.write(s + "\n")
    LOG.flush()


DJ = json.load(open(OUT + "prevstyle_data.json"))
SJ = json.load(open(OUT + "prevstyle_sim.json"))
PARTS = {k: tuple(v) for k, v in DJ["meta"]["parts"].items()}
WIN = [w[0] for w in DJ["meta"]["windows"]]
WLOHI = {w[0]: (w[1], w[2]) for w in DJ["meta"]["windows"]}
iW = {w: i for i, w in enumerate(WIN)}
MEEWIN = DJ["meta"]["meewin"]           # SBlo, Wpk, SBhi, PB
NMEE = 32


def part(arr, tag):
    o, n = PARTS[tag]
    return np.array(arr[o:o + n], dtype=float)


def D(s, tag, what="sig"):
    d = DJ["sel"][s]
    N, NV = part(d["N"], tag), part(d["NV"], tag)
    CB, CBe = part(d["CB"], tag), part(d["CBrep"], tag)
    if what == "N":
        return N, NV
    if what == "CB":
        return CB, CBe ** 2
    return N - CB, NV + CBe ** 2


def S(s, comps, tag):
    """SMASH expected raw counts and variance (sum of components)"""
    if isinstance(comps, str):
        comps = [comps]
    v = 0.
    e = 0.
    for c in comps:
        x = SJ["sel"][s][c]
        v = v + part(x["s"], tag) * NB2RAW
        e = e + part(x["v"], tag) * NB2RAW ** 2
    return v, e


GROUPS = {"C2": ["C2"], "C3": ["C3"], "C4": ["C4"], "C5": ["C5"], "rest": ["Geta", "GD", "GC"], "Geta": ["Geta"],
          "nonEta": ["GD", "GC"], "peak": ["C4", "C5"], "bump": ["C2", "C3"], "T12": ["T12"]}

# ---------------------------------------------------------------- checks
CHK = []


def check(what, val, ref, tol):
    ok = abs(val - ref) <= tol
    CHK.append(f"{what:<78s} {val:15.7g} ref {ref:15.7g}  {'OK' if ok else 'DIFFERS'}")
    return ok


E1D = json.load(open(E1DIR + "e1_data.json"))
E1S = json.load(open(E1DIR + "e1_sim.json"))
EP = {k: tuple(v) for k, v in E1D["meta"]["parts"].items()}


def epart(arr, tag):
    o, n = EP[tag]
    return np.array(arr[o:o + n], dtype=float)


for key in ["N", "NV", "CB", "CBrep", "Npp", "Nmm"]:
    a = part(DJ["sel"]["all"][key], "b5")
    b = epart(E1D["sel"]["pre_E1"][key], "b5")
    check(f"data isBest>=0 vs e1 pre_E1: max |diff| of {key} in the 5 MeV bins", float(np.max(np.abs(a - b))), 0., 1e-6)
for key in ["N", "CB", "CBrep"]:
    a = part(DJ["sel"]["all"][key], "mee")[3 * NMEE:4 * NMEE]
    b = epart(E1D["sel"]["pre_E1"][key], "mee")
    check(f"data isBest>=0 vs e1 pre_E1: M(ee) in the peak band, max |diff| of {key}", float(np.max(np.abs(a - b))), 0., 1e-6)
check("data isBest>=0: N in the e1 X17 cell (PB, 15.5-24.5 MeV)", part(DJ["sel"]["all"]["N"], "xr")[6], epart(E1D["sel"]["pre_E1"]["N"], "x17")[0], 1e-9)
check("data isBest>=0: CB in the e1 X17 cell", part(DJ["sel"]["all"]["CB"], "xr")[6], epart(E1D["sel"]["pre_E1"]["CB"], "x17")[0], 1e-6)
check("data isBest>=0: N in P = 46 001 (e1 note)", D("all", "win", "N")[0][iW["P"]], 46001, 0.5)
check("data isBest>=0: CB in P = 5628 (e1 note)", D("all", "win", "CB")[0][iW["P"]], 5628, 0.6)
for mine, theirs in [("C2iso", "C2"), ("C3", "C3"), ("C4", "C4"), ("C5", "C5"), ("Geta", "Geta"), ("GD", "GD"), ("GC", "Gconv"), ("C6", "C6")]:
    a = part(SJ["sel"]["all"][mine]["s"], "b5")
    b = epart(E1S["sel"]["pre_E1"][theirs]["s"], "b5")
    check(f"SMASH isBest>=0 {mine} vs e1 pre_E1 {theirs}: max rel |diff| (5 MeV bins)", float(np.max(np.abs(a - b)) / np.max(np.abs(b))), 0., 1e-9)
t12, _ = S("all", "T12", "b5")
csum = sum(S("all", c, "b5")[0] for c in ["C2", "C3", "C4", "C5", "Geta", "GD", "GC"])
check("SMASH isBest>=0: T12 = C2(QED) + C3 + C4 + C5 + G_eta + G_D + G_C (max |diff|, counts)", float(np.max(np.abs(t12 - csum))), 0., 1e-4)
t12b, _ = S("best", "T12", "b5")
csumb = sum(S("best", c, "b5")[0] for c in ["C2", "C3", "C4", "C5", "Geta", "GD", "GC"])
check("SMASH isBest==1: T12 = sum of the classes (max |diff|, counts)", float(np.max(np.abs(t12b - csumb))), 0., 1e-4)
q = SJ["qed"]
check("QED join: C2 counts with a weight = c2recon list", q["nC2pre"], q["nC2ref"], 0)
check("QED join: missing weights", q["miss"], 0, 0)
check("QED join: M4 mismatches", q["m4bad"], 0, 0)
check("QED: weighted mean q_A of all pre C2 counts (c2recon all C2, S0: 0.8637)", q["meanq"], 0.8637, 0.002)
for s in SELS:
    for w in ["SBlo", "Wpk", "SBhi"]:
        qi = MEEWIN.index(w)
        mN = part(DJ["sel"][s]["N"], "mee")[qi * NMEE:(qi + 1) * NMEE].sum()
        check(f"data {s}: sum of the M(ee) cells = window N ({w})", mN, D(s, "win", "N")[0][iW[w]], 1e-9)
        mC = part(DJ["sel"][s]["CB"], "mee")[qi * NMEE:(qi + 1) * NMEE].sum()
        check(f"data {s}: sum of the M(ee) CB cells = window CB ({w})", mC, D(s, "win", "CB")[0][iW[w]], 1e-6 * max(1, mC))
for c in CHK:
    P(c)
NBAD = sum("DIFFERS" in c for c in CHK)
P(f"checks: {len(CHK)} done, {NBAD} differ")

# ---------------------------------------------------------------- binning
XL = 0.30 + 0.005 * np.arange(160)            # lower edges of the b5 bins 1..160
XC = XL + 0.0025
RNG = (XL >= 0.35 - 1e-9) & (XL < 0.90 - 1e-9)
M47 = (XC > 0.40) & (XC < 0.70)
MGF = (XC > 0.43) & (XC < 0.75)
FXL = 0.30 + 0.0025 * np.arange(320)
FXC = FXL + 0.00125


def b5(a):
    return a[1:161]


def fine(a):
    return a[1:321]


def chi2_model(dv, dvar, mv, mvar, mask, npar):
    var = dvar + mvar
    m = mask & (var > 0)
    return float(np.sum((dv[m] - mv[m]) ** 2 / var[m])), int(m.sum()) - npar


def gfit(dv, dvar, T, Tvar, mask, it=100):
    """chi2 fit d = sum_k a_k T_k, variance = dvar + sum_k a_k^2 Tvar_k (iterated); returns a, cov, chi2, ndf"""
    K = len(T)
    a = np.ones(K)
    m = mask & (dvar > 0)
    for _ in range(it):
        var = dvar + sum(a[k] ** 2 * Tvar[k] for k in range(K))
        w = 1. / var[m]
        A = np.array([[np.sum(T[i][m] * T[j][m] * w) for j in range(K)] for i in range(K)])
        b = np.array([np.sum(dv[m] * T[i][m] * w) for i in range(K)])
        an = np.linalg.solve(A, b)
        if np.max(np.abs(an - a)) < 1e-10:
            a = an
            break
        a = an
    var = dvar + sum(a[k] ** 2 * Tvar[k] for k in range(K))
    w = 1. / var[m]
    A = np.array([[np.sum(T[i][m] * T[j][m] * w) for j in range(K)] for i in range(K)])
    cov = np.linalg.inv(A)
    mod = sum(a[k] * T[k] for k in range(K))
    chi2 = float(np.sum((dv[m] - mod[m]) ** 2 * w))
    return a, cov, chi2, int(m.sum()) - K


# ---------------------------------------------------------------- smoothing functions (earlier analysis: Voigt / Weibull / A(x-B)^C exp(-D(x-B)^E))
def f_voigt(x, p):
    A, m, sg, gm = p
    return A * voigt_profile(x - m, abs(sg), abs(gm))


def f_weib(x, p):
    A, B, l, k = p
    u = np.clip((x - B) / l, 0, None)
    return np.where(x > B, A * (k / l) * u ** (k - 1) * np.exp(-u ** k), 0.)


def f_bkg(x, p):
    A, B, C, Dd, E = p
    u = np.clip(x - B, 0, None)
    return np.where(x > B, A * u ** C * np.exp(-Dd * u ** E), 0.)


SUB = (np.arange(8) + 0.5) / 8.


def bin_int(f, p, lo, w):
    """integral of the density f over bins [lo, lo + w) (8-point midpoint rule)"""
    xs = lo[:, None] + w * SUB[None, :]
    return f(xs, p).mean(axis=1) * w


def smooth_fit(f, p0, v, var, lo_fit, hi_fit):
    """weighted-Poisson deviance fit of the density f to the fine template v (counts per 2.5 MeV bin) in [lo_fit, hi_fit)"""
    m = (FXL >= lo_fit - 1e-9) & (FXL < hi_fit - 1e-9)
    wbar = float(np.sum(var[m]) / np.sum(v[m]))
    n = np.clip(v[m], 0, None) / wbar
    lo = FXL[m]

    def dev(p):
        mu = bin_int(f, p, lo, 0.0025) / wbar
        if not np.all(np.isfinite(mu)) or np.any(mu < 0):
            return 1e30
        mu = np.maximum(mu, 1e-12)
        t = mu - n
        nz = n > 0
        t[nz] += n[nz] * np.log(n[nz] / mu[nz])
        return 2. * float(np.sum(t))

    best = None
    p = np.array(p0, float)
    for meth in ["Nelder-Mead", "Powell", "Nelder-Mead"]:
        r = minimize(dev, p, method=meth, options=dict(maxiter=40000, maxfev=40000, xatol=1e-10, fatol=1e-10) if meth == "Nelder-Mead"
                     else dict(maxiter=40000, xtol=1e-10, ftol=1e-12))
        if best is None or r.fun < best.fun:
            best = r
        p = best.x
    ndf = int(m.sum()) - len(p0)
    return best.x, best.fun, ndf, wbar


# ---------------------------------------------------------------- per selection
R = {}
for s in SELS:
    r = {}
    dv, dvar = (b5(x) for x in D(s, "b5"))
    Nv, NVv = (b5(x) for x in D(s, "b5", "N"))
    CBv, CBvar = (b5(x) for x in D(s, "b5", "CB"))
    comp = {g: tuple(b5(x) for x in S(s, cs, "b5")) for g, cs in GROUPS.items()}
    compf = {g: tuple(fine(x) for x in S(s, cs, "fine")) for g, cs in GROUPS.items()}
    compw = {g: S(s, cs, "win") for g, cs in GROUPS.items()}
    dw, dwv = D(s, "win")
    Nw = D(s, "win", "N")[0]
    CBw, CBwv = D(s, "win", "CB")
    r.update(dv=dv, dvar=dvar, N=Nv, NV=NVv, CB=CBv, CBvar=CBvar, comp=comp, compf=compf)
    # integrals
    ia = iW["all"]
    r["data_all"] = (float(dw[ia]), float(math.sqrt(dwv[ia])))
    r["N_all"], r["CB_all"] = float(Nw[ia]), float(CBw[ia])
    r["T12_all"] = (float(compw["T12"][0][ia]), float(math.sqrt(compw["T12"][1][ia])))
    # (a) SMASH scaled to the data integral in 0.35-0.90
    sint = dw[ia] / compw["T12"][0][ia]
    sinte = sint * math.sqrt(dwv[ia] / dw[ia] ** 2 + compw["T12"][1][ia] / compw["T12"][0][ia] ** 2)
    r["sint"], r["sinte"] = float(sint), float(sinte)
    T, Tv = comp["T12"]
    r["chi2a"] = chi2_model(dv, dvar, sint * T, sint ** 2 * Tv, M47, 1)
    # (b) measured scales; c_b from 0.575-0.70
    ic = iW["CR"]
    etab_cr = F_BUMP * compw["bump"][0][ic] + A_P * compw["peak"][0][ic] + G_ETA * compw["Geta"][0][ic]
    etab_crv = F_BUMP ** 2 * compw["bump"][1][ic] + A_P ** 2 * compw["peak"][1][ic] + G_ETA ** 2 * compw["Geta"][1][ic]
    ne, nev = compw["nonEta"][0][ic], compw["nonEta"][1][ic]
    cb = (dw[ic] - etab_cr) / ne
    cbe = math.sqrt(dwv[ic] + etab_crv + cb ** 2 * nev) / ne
    eta40 = sum(compw[g][0][ic] for g in ["bump", "peak", "Geta"])
    cb40 = (dw[ic] - 0.40 * eta40) / ne
    r["cb"], r["cbe"], r["cb40"] = float(cb), float(cbe), float(cb40)
    mb = F_BUMP * comp["bump"][0] + A_P * comp["peak"][0] + G_ETA * comp["Geta"][0] + cb * comp["nonEta"][0]
    mbv = F_BUMP ** 2 * comp["bump"][1] + A_P ** 2 * comp["peak"][1] + G_ETA ** 2 * comp["Geta"][1] + cb ** 2 * comp["nonEta"][1]
    r["mb"], r["mbv"] = mb, mbv
    r["chi2b"] = chi2_model(dv, dvar, mb, mbv, M47, 1)
    # same with c_b = 1.450 (package e1) for the record
    mb145 = F_BUMP * comp["bump"][0] + A_P * comp["peak"][0] + G_ETA * comp["Geta"][0] + C_B_E1 * comp["nonEta"][0]
    r["chi2b145"] = chi2_model(dv, dvar, mb145, mbv, M47, 0)
    # region ratios data / model
    r["reg"] = {}
    for w in ["below", "B", "P", "Wpk", "r56", "r70", "fit1", "all"]:
        i = iW[w]
        ma = sint * compw["T12"][0][i]
        mbw = (F_BUMP * compw["bump"][0][i] + A_P * compw["peak"][0][i] + G_ETA * compw["Geta"][0][i] + cb * compw["nonEta"][0][i])
        r["reg"][w] = dict(data=float(dw[i]), err=float(math.sqrt(dwv[i])), N=float(Nw[i]), CB=float(CBw[i]), a=float(dw[i] / ma), b=float(dw[i] / mbw),
                           T12=float(compw["T12"][0][i]))
    # eta shares of SMASH as is in the windows
    r["share"] = {w: {g: float(compw[g][0][iW[w]] / compw["T12"][0][iW[w]]) for g in ["bump", "peak", "Geta", "nonEta"]} for w in ["B", "P", "Wpk"]}
    # ---- F2: three-scale global fit, raw templates (components = SMASH scaled to the data integral, as before the fit)
    Tk = [sint * comp["peak"][0], sint * comp["bump"][0], sint * comp["rest"][0]]
    Tkv = [sint ** 2 * comp["peak"][1], sint ** 2 * comp["bump"][1], sint ** 2 * comp["rest"][1]]
    a, cov, chi2, ndf = gfit(dv, dvar, Tk, Tkv, MGF)
    r["gf_raw"] = dict(a=a.tolist(), e=np.sqrt(np.diag(cov)).tolist(), cov=cov.tolist(), chi2=chi2, ndf=ndf,
                       eta=(a * sint).tolist(), etae=(np.sqrt(np.diag(cov)) * sint).tolist(),
                       rho=[float(cov[i, j] / math.sqrt(cov[i, i] * cov[j, j])) for i, j in [(0, 1), (0, 2), (1, 2)]])
    r["Tk"], r["Tkv"] = Tk, Tkv
    # data-error-only variant (no SMASH statistics in the variance)
    a0, cov0, chi20, ndf0 = gfit(dv, dvar, Tk, [0 * t for t in Tkv], MGF, it=1)
    r["gf_raw_dataonly"] = dict(a=a0.tolist(), e=np.sqrt(np.diag(cov0)).tolist(), chi2=chi20, ndf=ndf0, eta=(a0 * sint).tolist())
    # ---- smoothed templates (fit each SMASH template on the 2.5 MeV grid)
    sm = {}
    pk, pkv = compf["peak"]
    pp, dev, nd, wb = smooth_fit(f_voigt, [pk.sum(), 0.5437, 0.0045, 0.0010], pk, pkv, 0.46, 0.62)
    sm["peak"] = dict(f=f_voigt, p=pp, dev=dev, ndf=nd, wbar=wb, rng=(0.46, 0.62), name="Voigt")
    bm, bmv = compf["bump"]
    pb, dev, nd, wb = smooth_fit(f_weib, [bm.sum(), 0.40, 0.10, 5.0], bm, bmv, 0.40, 0.58)
    sm["bump"] = dict(f=f_weib, p=pb, dev=dev, ndf=nd, wbar=wb, rng=(0.40, 0.58), name="Weibull")
    rs, rsv = compf["rest"]
    # start: A (x-B)^C exp(-D (x-B)^E) with the mode near 0.62
    B0, C0, E0 = 0.425, 1.5, 1.0
    D0 = C0 / 0.2
    xs = np.linspace(B0, 1.3, 4000)
    g0 = f_bkg(xs, [1., B0, C0, D0, E0])
    A0 = rs.sum() / np.trapz(g0, xs)
    # nominal: smoothed inside the global-fit range 0.43-0.75 (the function cannot follow the sharp threshold at 0.43 that the
    # alpha_eta > 140 deg cut makes together with the tail above 0.75); variant "wide": 0.40-0.90 (sensitivity only)
    pr, dev, nd, wb = smooth_fit(f_bkg, [A0, B0, C0, D0, E0], rs, rsv, 0.43, 0.75)
    sm["rest"] = dict(f=f_bkg, p=pr, dev=dev, ndf=nd, wbar=wb, rng=(0.43, 0.75), name="A(x-B)^C exp(-D(x-B)^E)")
    prw, devw, ndw, wbw = smooth_fit(f_bkg, [A0, B0, C0, D0, E0], rs, rsv, 0.40, 0.90)
    sm["rest_wide"] = dict(f=f_bkg, p=prw, dev=devw, ndf=ndw, wbar=wbw, rng=(0.40, 0.90), name="A(x-B)^C exp(-D(x-B)^E), 0.40-0.90")
    for g in sm:
        sm[g]["b5"] = bin_int(sm[g]["f"], sm[g]["p"], XL, 0.005)
        sm[g]["fine"] = bin_int(sm[g]["f"], sm[g]["p"], FXL, 0.0025)
        raw = comp[g.replace("_wide", "")][0]
        sm[g]["int_gf"] = (float(sm[g]["b5"][MGF].sum()), float(raw[MGF].sum()))
    r["sm"] = sm
    Tks = [sint * sm["peak"]["b5"], sint * sm["bump"]["b5"], sint * sm["rest"]["b5"]]
    a2, cov2, chi22, ndf2 = gfit(dv, dvar, Tks, [0 * t for t in Tks], MGF, it=1)
    r["gf_sm"] = dict(a=a2.tolist(), e=np.sqrt(np.diag(cov2)).tolist(), cov=cov2.tolist(), chi2=chi22, ndf=ndf2,
                      eta=(a2 * sint).tolist(), etae=(np.sqrt(np.diag(cov2)) * sint).tolist(),
                      rho=[float(cov2[i, j] / math.sqrt(cov2[i, i] * cov2[j, j])) for i, j in [(0, 1), (0, 2), (1, 2)]])
    r["Tks"] = Tks
    # ---- F3: M(ee) in the peak window and the side bands
    def mee(tagw, src="data"):
        qi = MEEWIN.index(tagw)
        sl = slice(qi * NMEE, (qi + 1) * NMEE)
        if src == "data":
            v_, e_ = D(s, "mee")
            return v_[sl], e_[sl]
        v_, e_ = S(s, GROUPS.get(src, [src]), "mee")
        return v_[sl], e_[sl]

    def xr(tagw, k, src="data"):
        qi = MEEWIN.index(tagw)
        if src == "data":
            v_, e_ = D(s, "xr")
        else:
            v_, e_ = S(s, GROUPS.get(src, [src]), "xr")
        return float(v_[qi * 2 + k]), float(e_[qi * 2 + k])

    Wv, Wvar = mee("Wpk")
    Lv, Lvar = mee("SBlo")
    Hv, Hvar = mee("SBhi")
    SBv, SBvar = Lv + Hv, Lvar + Hvar
    iWp, iL, iH = iW["Wpk"], iW["SBlo"], iW["SBhi"]
    sb_tot = dw[iL] + dw[iH]
    sb_totv = dwv[iL] + dwv[iH]
    gr = r["gf_raw"]
    aa = np.array(gr["a"])
    # expected background in the window from the raw-template global fit: a2 s (C2 + C3) + a3 s (rest), exact M4' window sums
    bump_W, rest_W, peak_W = (sint * compw[g][0][iWp] for g in ["bump", "rest", "peak"])
    BW = aa[1] * bump_W + aa[2] * rest_W
    cv = np.array(gr["cov"])
    BWe = math.sqrt(bump_W ** 2 * cv[1, 1] + rest_W ** 2 * cv[2, 2] + 2 * bump_W * rest_W * cv[1, 2])
    kSB = BW / sb_tot
    kSBe = kSB * math.sqrt((BWe / BW) ** 2 + sb_totv / sb_tot ** 2)
    # alternative: smoothed-template fit, and the plain width ratio 27/24
    a_s = np.array(r["gf_sm"]["a"])
    Tksw = [sint * sm["peak"]["b5"], sint * sm["bump"]["b5"], sint * sm["rest_wide"]["b5"]]
    a3, cov3, chi23, ndf3 = gfit(dv, dvar, Tksw, [0 * t for t in Tksw], MGF, it=1)
    r["gf_sm_wide"] = dict(a=a3.tolist(), e=np.sqrt(np.diag(cov3)).tolist(), chi2=chi23, ndf=ndf3, eta=(a3 * sint).tolist(),
                           etae=(np.sqrt(np.diag(cov3)) * sint).tolist())
    sm_bump_W = sint * float(np.sum(sm["bump"]["fine"] * np.array([max(0., min(FXL[i] + 0.0025, 0.557) - max(FXL[i], 0.530)) / 0.0025 for i in range(320)])))
    sm_rest_W = sint * float(np.sum(sm["rest"]["fine"] * np.array([max(0., min(FXL[i] + 0.0025, 0.557) - max(FXL[i], 0.530)) / 0.0025 for i in range(320)])))
    kSB_sm = (a_s[1] * sm_bump_W + a_s[2] * sm_rest_W) / sb_tot
    kSB_w = 27. / 24.
    PW = aa[0] * peak_W
    r["sb"] = dict(BW=float(BW), BWe=float(BWe), SB=float(sb_tot), SBe=float(math.sqrt(sb_totv)), k=float(kSB), ke=float(kSBe),
                   k_sm=float(kSB_sm), k_w=kSB_w, W=float(dw[iWp]), We=float(math.sqrt(dwv[iWp])), PW=float(PW),
                   NW=float(Nw[iWp]), CBW=float(CBw[iWp]), NSB=float(Nw[iL] + Nw[iH]), CBSB=float(CBw[iL] + CBw[iH]),
                   peak_in_SB=float(aa[0] * sint * (compw["peak"][0][iL] + compw["peak"][0][iH])),
                   bump_W=float(bump_W), rest_W=float(rest_W), peak_W=float(peak_W))
    r["mee"] = dict(W=Wv, Wvar=Wvar, SB=SBv, SBvar=SBvar, NW=mee("Wpk")[0] + 0, )
    r["mee"]["NWraw"] = D(s, "mee", "N")[0][MEEWIN.index("Wpk") * NMEE:(MEEWIN.index("Wpk") + 1) * NMEE]
    r["mee"]["CBW"] = D(s, "mee", "CB")[0][MEEWIN.index("Wpk") * NMEE:(MEEWIN.index("Wpk") + 1) * NMEE]
    for g in ["C4", "C5", "peak", "C2", "C3", "rest", "C6"]:
        r["mee"][g] = mee("Wpk", g)
    # X17 windows (exact cells): k = 0 15.5-24.5, k = 1 12-22
    r["x17"] = {}
    for k, nm in [(0, "15.5-24.5"), (1, "12-22")]:
        w_, wv_ = xr("Wpk", k)
        l_, lv_ = xr("SBlo", k)
        h_, hv_ = xr("SBhi", k)
        sbx, sbxv = l_ + h_, lv_ + hv_
        diff = w_ - kSB * sbx
        diffe = math.sqrt(wv_ + kSB ** 2 * sbxv)
        pkx = aa[0] * sint * xr("Wpk", k, "peak")[0]
        c5x = aa[0] * sint * xr("Wpk", k, "C5")[0]
        c4x = aa[0] * sint * xr("Wpk", k, "C4")[0]
        c6x = xr("Wpk", k, "C6")[0]
        pkxv = (aa[0] * sint) ** 2 * xr("Wpk", k, "peak")[1]
        pkaP = A_P * xr("Wpk", k, "peak")[0]
        r["x17"][nm] = dict(W=w_, We=math.sqrt(wv_), SBs=kSB * sbx, SBse=kSB * math.sqrt(sbxv), diff=diff, diffe=diffe, peak=pkx, C5=c5x, C4=c4x,
                            excess=diff - pkx, C6share=c6x / S(s, "C6", "win")[0][iWp],
                            # systematic views of the same number: error of k, SMASH statistics of the eta peak, other k, other eta scale
                            e_k=sbx * kSBe, e_mc=math.sqrt(pkxv),
                            excess_kw=w_ - kSB_w * sbx - pkx, excess_ksm=w_ - kSB_sm * sbx - pkx,
                            peak_aP=pkaP, excess_aP=diff - pkaP, SBraw=sbx)
    # chi2 of (window - scaled side bands) against the SMASH eta peak (a1 s (C4 + C5)) in 0-50 MeV
    dd = Wv[:25] - kSB * SBv[:25]
    ddv = Wvar[:25] + kSB ** 2 * SBvar[:25]
    pkm, pkmv = r["mee"]["peak"]
    pm = aa[0] * sint * pkm[:25]
    pmv = (aa[0] * sint) ** 2 * pkmv[:25]
    r["mee_chi2"] = chi2_model(dd, ddv, pm, pmv, np.ones(25, bool), 0)
    r["mee_sum"] = dict(diff=float(dd.sum()), diffe=float(math.sqrt(ddv.sum())), peak=float(pm.sum()))
    R[s] = r

# ---------------------------------------------------------------- print numbers
P("\n==== F1: M4 after E1 ====")
for s in SELS:
    r = R[s]
    P(f"[{SELSHORT[s]}] data N 0.35-0.90 = {r['N_all']:.0f}, CB = {r['CB_all']:.1f}, N - CB = {r['data_all'][0]:.1f} +- {r['data_all'][1]:.1f}; "
      f"SMASH T12 = {r['T12_all'][0]:.1f} +- {r['T12_all'][1]:.1f} expected counts")
    P(f"  (a) integral scale s = {r['sint']:.4f} +- {r['sinte']:.4f}; chi2/ndf 0.40-0.70 = {r['chi2a'][0]:.1f}/{r['chi2a'][1]} = {r['chi2a'][0] / r['chi2a'][1]:.2f}")
    P(f"  (b) c_b = {r['cb']:.4f} +- {r['cbe']:.4f} (with all eta x 0.40 as in e1: {r['cb40']:.4f}); chi2/ndf = {r['chi2b'][0]:.1f}/{r['chi2b'][1]} = "
      f"{r['chi2b'][0] / r['chi2b'][1]:.2f}; with c_b = 1.450: {r['chi2b145'][0]:.1f}/{r['chi2b145'][1]}")
    for w, d in r["reg"].items():
        P(f"    {w:6s} data N-CB {d['data']:10.1f} +- {d['err']:7.1f}  (N {d['N']:9.0f}, CB {d['CB']:8.1f})  data/(a) {d['a']:.3f}  data/(b) {d['b']:.3f}")
    for w, d in r["share"].items():
        P(f"    SMASH as is, shares in {w}: bump {d['bump']:.3f}, peak {d['peak']:.3f}, other eta {d['Geta']:.3f}, non-eta {d['nonEta']:.3f}")
P("\n==== F2: three-scale global fit in 0.43-0.75 (5 MeV bins) ====")
NAMES = ["alpha1 (C4+C5)", "alpha2 (C2+C3)", "alpha3 (rest)"]
for s in SELS:
    r = R[s]
    for key, txt in [("gf_raw", "raw SMASH templates, variance = data + SMASH statistics"), ("gf_raw_dataonly", "raw templates, data errors only"),
                     ("gf_sm", "smoothed templates (Voigt / Weibull / A(x-B)^C exp(-D(x-B)^E) in 0.43-0.75), data errors"),
                     ("gf_sm_wide", "smoothed templates, rest smoothed in 0.40-0.90 (sensitivity), data errors")]:
        g = r[key]
        P(f"[{SELSHORT[s]}] {txt}: chi2/ndf = {g['chi2']:.1f}/{g['ndf']} = {g['chi2'] / g['ndf']:.2f}")
        for k in range(3):
            P(f"    {NAMES[k]} = {g['a'][k]:.4f} +- {g['e'][k]:.4f}  -> scale relative to SMASH = alpha x s = {g['eta'][k]:.4f}"
              + (f" +- {g['etae'][k]:.4f}" if "etae" in g else ""))
        if "rho" in g:
            P(f"    correlations (1,2) {g['rho'][0]:.3f}, (1,3) {g['rho'][1]:.3f}, (2,3) {g['rho'][2]:.3f}")
    for g in ["peak", "bump", "rest", "rest_wide"]:
        sm = r["sm"][g]
        P(f"    smoothing {g}: {sm['name']} p = {np.array2string(np.array(sm['p']), precision=5)}; deviance {sm['dev']:.1f}/{sm['ndf']} in "
          f"{sm['rng'][0]:.2f}-{sm['rng'][1]:.2f}; integral 0.43-0.75 smooth {sm['int_gf'][0]:.1f} vs raw {sm['int_gf'][1]:.1f} "
          f"({sm['int_gf'][0] / sm['int_gf'][1]:.4f}); mean weight per simulated count {sm['wbar']:.2f}")
P("\n==== F3: M(ee) in the window 530-557 against the side bands 518-530 + 557-569 ====")
for s in SELS:
    r = R[s]
    b = r["sb"]
    P(f"[{SELSHORT[s]}] window N {b['NW']:.0f}, CB {b['CBW']:.1f}, N - CB {b['W']:.1f} +- {b['We']:.1f}; side bands N {b['NSB']:.0f}, CB {b['CBSB']:.1f}, "
      f"N - CB {b['SB']:.1f} +- {b['SBe']:.1f}")
    P(f"    expected background in the window (global fit, raw templates): {b['BW']:.1f} +- {b['BWe']:.1f}; eta peak (a1 s (C4+C5)) {b['PW']:.1f}; "
      f"side-band scale k = {b['k']:.4f} +- {b['ke']:.4f} (smoothed-template fit {b['k_sm']:.4f}; width ratio 27/24 = {b['k_w']:.4f}); "
      f"eta peak inside the side bands {b['peak_in_SB']:.1f}")
    P(f"    window - k x side bands, 0-50 MeV: {r['mee_sum']['diff']:.1f} +- {r['mee_sum']['diffe']:.1f}; SMASH eta peak there {r['mee_sum']['peak']:.1f}; "
      f"chi2 (difference vs SMASH eta peak shape, 25 bins) {r['mee_chi2'][0]:.1f}/{r['mee_chi2'][1]}")
    for nm, d in r["x17"].items():
        P(f"    M(ee) {nm} MeV: window {d['W']:.1f} +- {d['We']:.1f}; scaled side bands {d['SBs']:.1f} +- {d['SBse']:.1f}; difference {d['diff']:.1f} +- "
          f"{d['diffe']:.1f}; SMASH eta peak {d['peak']:.1f} (C5 {d['C5']:.1f}, C4 {d['C4']:.1f}); difference - eta peak {d['excess']:.1f} +- {d['diffe']:.1f}; "
          f"X17 share of the window in this range {d['C6share']:.3f}")
        P(f"        same difference - eta peak: error of k {d['e_k']:.1f}, SMASH statistics of the eta peak {d['e_mc']:.1f}; with k = 27/24: "
          f"{d['excess_kw']:.1f}; with k of the smoothed-template fit: {d['excess_ksm']:.1f}; with the eta peak at a_P = {A_P} "
          f"(SMASH eta peak {d['peak_aP']:.1f}): {d['excess_aP']:.1f}")


# ---------------------------------------------------------------- drawing helpers
def steps(ax, y, color, lw=1.6, ls="-", label=None, rng=RNG, alpha=1.):
    xl = XL[rng]
    xe = np.r_[xl, xl[-1] + 0.005]
    yy = np.r_[y[rng], y[rng][-1]]
    ax.step(xe, yy, where="post", color=color, lw=lw, ls=ls, label=label, alpha=alpha)


def ratio_pad(ax, dv, dvar, mv, mvar, rng=RNG, color=None, ylim=(0.5, 1.5), mask_min=0.02):
    ok = rng & (mv > mask_min * np.max(mv[rng]))
    rr = dv[ok] / mv[ok]
    er = np.sqrt(dvar[ok]) / mv[ok]
    xl = XL[rng]
    okr = ok[rng]
    band = np.where(okr, np.sqrt(mvar[rng]) / np.where(mv[rng] > 0, mv[rng], 1), np.nan)
    ax.fill_between(np.r_[xl, xl[-1] + 0.005], np.r_[1 - band, (1 - band)[-1]], np.r_[1 + band, (1 + band)[-1]], step="post",
                    color="#bbbbbb", alpha=0.6, lw=0)
    ax.axhline(1, color="k", lw=0.8)
    ax.errorbar(XC[ok], rr, er, fmt="o", ms=3, color=color or PAL["sig"], lw=0.9)
    ax.set_ylim(*ylim)


def fmt_pm(v, e, d=3):
    return f"{v:.{d}f} $\\pm$ {e:.{d}f}"


# ================================================================ F1
from matplotlib.gridspec import GridSpecFromSubplotSpec

SELTITLE = {"all": r"isBest $\geq$ 0: all stored combinations (nominal)",
            "best": r"isBest == 1: best combination only, as in the earlier analysis (comparison)"}


def panel_grid(fig, top, bottom):
    outer = GridSpec(2, 2, figure=fig, hspace=0.16, wspace=0.17, left=0.065, right=0.985, top=top, bottom=bottom)
    cells = {}
    for row in range(2):
        for col in range(2):
            inner = GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[row, col], height_ratios=[3, 1.05], hspace=0.06)
            ax = fig.add_subplot(inner[0])
            axr = fig.add_subplot(inner[1], sharex=ax)
            cells[(row, col)] = (ax, axr)
    return cells


fig = plt.figure(figsize=(17, 16))
cells = panel_grid(fig, 0.862, 0.05)
ymax = 0
for s in SELS:
    ymax = max(ymax, np.max((R[s]["N"] + np.sqrt(R[s]["NV"]))[RNG]))
ymax *= 1.30
for row, s in enumerate(SELS):
    r = R[s]
    for col in range(2):
        ax, axr = cells[(row, col)]
        if col == 0:
            sc = {g: r["sint"] for g in ["C2", "C3", "C4", "C5", "rest"]}
            tot, totv = r["sint"] * r["comp"]["T12"][0], r["sint"] ** 2 * r["comp"]["T12"][1]
            restv = r["sint"] * r["comp"]["rest"][0]
        else:
            sc = {"C2": F_BUMP, "C3": F_BUMP, "C4": A_P, "C5": A_P}
            tot, totv = r["mb"], r["mbv"]
            restv = G_ETA * r["comp"]["Geta"][0] + r["cb"] * r["comp"]["nonEta"][0]
        m = RNG
        ax.errorbar(XC[m], r["N"][m], np.sqrt(r["NV"][m]), fmt="o", mfc="white", mec=PAL["N"], ecolor=PAL["N"], ms=4.2, lw=0.8, label="data $N$")
        ax.errorbar(XC[m], r["CB"][m], np.sqrt(r["CBvar"][m]), fmt="s", color=PAL["CB"], ms=3.2, lw=0.8, label=r"CB $=2\sqrt{N_{++}N_{--}}$ (data)")
        ax.errorbar(XC[m], r["dv"][m], np.sqrt(r["dvar"][m]), fmt="o", color=PAL["sig"], ms=4.2, lw=0.9, label="data $N-$CB")
        steps(ax, tot, PAL["tot"], lw=2.0, label="SMASH total")
        for g in ["C5", "C4", "C2", "C3"]:
            steps(ax, sc[g] * r["comp"][g][0], PAL[g], lw=1.6, label=CLAB[g])
        steps(ax, restv, PAL["rest"], lw=1.6, label=CLAB["rest"])
        ax.set_xlim(0.35, 0.90)
        ax.set_ylim(0, ymax)
        ax.set_ylabel(r"Counts per 5 MeV/$c^2$")
        plt.setp(ax.get_xticklabels(), visible=False)
        if col == 0:
            txt = (f"SMASH $\\times$ {r['sint']:.3f}, scaled to the data\nintegral in 0.35$-$0.90 GeV/$c^2$\n"
                   f"$\\chi^2$/ndf (0.40$-$0.70) = {r['chi2a'][0]:.0f}/{r['chi2a'][1]} = {r['chi2a'][0] / r['chi2a'][1]:.1f}")
        else:
            txt = (f"$\\eta\\to\\pi^+\\pi^-\\pi^0$ $\\times$ {F_BUMP}, peak decays $\\times$ {A_P},\n"
                   f"other $\\eta$ $\\times$ {G_ETA}, non-$\\eta$ $\\times$ $c_b$ = {r['cb']:.3f}\n"
                   f"($c_b$ from the data in 0.575$-$0.70 GeV/$c^2$)\n"
                   f"$\\chi^2$/ndf (0.40$-$0.70) = {r['chi2b'][0]:.0f}/{r['chi2b'][1]} = {r['chi2b'][0] / r['chi2b'][1]:.1f}")
        ax.text(0.985, 0.97, txt, transform=ax.transAxes, ha="right", va="top", fontsize=11.5,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#999999", alpha=0.92))
        ax.set_title(SELTITLE[s], loc="left", fontsize=12, fontweight="bold")
        if row == 0 and col == 0:
            hl, ll = ax.get_legend_handles_labels()
            order = [6, 8, 7, 0, 1, 2, 3, 4, 5] if len(hl) == 9 else list(range(len(hl)))
            fig.legend([hl[i] for i in order], [ll[i] for i in order], loc="upper center", bbox_to_anchor=(0.5, 0.972), ncol=5,
                       fontsize=11.5, frameon=False, columnspacing=1.6)
        ratio_pad(axr, r["dv"], r["dvar"], tot, totv)
        axr.set_ylabel("data / model", fontsize=13)
        if row == 1:
            axr.set_xlabel(r"$M(\pi^+\pi^-e^+e^-)$ [GeV/$c^2$]")
        else:
            plt.setp(axr.get_xticklabels(), visible=False)
        for xx in (0.40, 0.70):
            axr.axvline(xx, color="#666666", ls=":", lw=1)
fig.text(0.29, 0.897, "(a) SMASH scaled to the data integral, as in the earlier analysis", ha="center", fontsize=14.5, fontweight="bold")
fig.text(0.775, 0.897, "(b) SMASH at its cross sections with our measured scales", ha="center", fontsize=14.5, fontweight="bold")
fig.suptitle(r"$pp\to\pi^+\pi^-e^+e^-X$ at 4.5 GeV after the typical cuts ($\alpha_{lab}<50^\circ$, $\alpha_\eta>140^\circ$, "
             r"$M(\pi^+\pi^-)<0.46$ GeV/$c^2$, OA$(e^+e^-)>2^\circ$): data and SMASH in the style of the earlier analysis",
             fontsize=14.5, y=0.993)
for ext in ("pdf", "png"):
    fig.savefig(FIG + f"prevstyle_m4_E1.{ext}", dpi=110)
plt.close(fig)
P("written prevstyle_m4_E1")

# ================================================================ F2
fig = plt.figure(figsize=(17, 16))
cells = panel_grid(fig, 0.885, 0.05)
ymax2 = 0
for s in SELS:
    ymax2 = max(ymax2, np.max((R[s]["dv"] + np.sqrt(R[s]["dvar"]))[RNG]))
ymax2 *= 1.50
GFR = RNG & (XC > 0.43) & (XC < 0.75)      # model lines only over the fit range
for row, key in enumerate(["gf_raw", "gf_sm"]):
    for col, s in enumerate(SELS):
        r = R[s]
        g = r[key]
        a = np.array(g["a"])
        Tk = r["Tk"] if key == "gf_raw" else r["Tks"]
        Tkv = r["Tkv"] if key == "gf_raw" else [0 * t for t in Tk]
        ax, axr = cells[(row, col)]
        ax.axvspan(0.43, 0.75, color="#f2f2f2", zorder=0)
        m = RNG
        ax.errorbar(XC[m], r["dv"][m], np.sqrt(r["dvar"][m]), fmt="o", color=PAL["sig"], ms=4.2, lw=0.9, label="data $N-$CB", zorder=5)
        pre = sum(Tk)
        steps(ax, pre, PAL["pre"], lw=1.8, ls=":", label="sum of components before global fit\n(SMASH scaled to the data integral)", rng=GFR)
        fitm = sum(a[k] * Tk[k] for k in range(3))
        steps(ax, fitm, PAL["fit"], lw=2.2, label="global fit", rng=GFR)
        steps(ax, a[0] * Tk[0], PAL["C5"], lw=1.7, ls="--", label=r"$\alpha_1$ [$\eta\to\pi^+\pi^-e^+e^-$ + $\eta\to\pi^+\pi^-\gamma(e^+e^-)$]", rng=GFR)
        steps(ax, a[1] * Tk[1], PAL["C2"], lw=1.7, ls="--", label=r"$\alpha_2$ [both $\eta\to\pi^+\pi^-\pi^0$ decays]", rng=GFR)
        steps(ax, a[2] * Tk[2], PAL["rest"], lw=1.7, ls="--", label=r"$\alpha_3$ [rest: multipion background]", rng=GFR)
        ax.set_xlim(0.35, 0.90)
        ax.set_ylim(0, ymax2)
        ax.set_ylabel(r"Counts per 5 MeV/$c^2$")
        plt.setp(ax.get_xticklabels(), visible=False)
        sint = r["sint"]
        txt = (f"$\\alpha_1$ = {fmt_pm(g['a'][0], g['e'][0])}  ($\\eta$ peak: {g['eta'][0]:.3f} $\\pm$ {g['etae'][0]:.3f} of SMASH)\n"
               f"$\\alpha_2$ = {fmt_pm(g['a'][1], g['e'][1])}  ($\\eta\\to\\pi^+\\pi^-\\pi^0$: {g['eta'][1]:.3f} $\\pm$ {g['etae'][1]:.3f} of SMASH)\n"
               f"$\\alpha_3$ = {fmt_pm(g['a'][2], g['e'][2])}  (rest: {g['eta'][2]:.3f} $\\pm$ {g['etae'][2]:.3f} of SMASH)\n"
               f"$\\chi^2$/ndf = {g['chi2']:.0f}/{g['ndf']} = {g['chi2'] / g['ndf']:.2f}   (scale of SMASH = $\\alpha$ $\\times$ {sint:.3f})")
        ax.text(0.985, 0.97, txt, transform=ax.transAxes, ha="right", va="top", fontsize=11,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#999999", alpha=0.92), zorder=6)
        ax.set_title(SELTITLE[s] + "\n" + ["raw SMASH templates (errors: data and SMASH statistics)",
                                           "SMASH templates smoothed with Voigt / Weibull / $A(x-B)^C e^{-D(x-B)^E}$ (errors: data only)"][row],
                     loc="left", fontsize=12)
        if row == 0 and col == 0:
            hl, ll = ax.get_legend_handles_labels()
            order = [5, 0, 1, 2, 3, 4] if len(hl) == 6 else list(range(len(hl)))
            fig.legend([hl[i] for i in order], [ll[i] for i in order], loc="upper center", bbox_to_anchor=(0.5, 0.962), ncol=3,
                       fontsize=11.5, frameon=False, columnspacing=1.6)
        ratio_pad(axr, r["dv"], r["dvar"], fitm, sum(a[k] ** 2 * Tkv[k] for k in range(3)), rng=GFR, ylim=(0.6, 1.4))
        axr.axvspan(0.43, 0.75, color="#f2f2f2", zorder=0)
        axr.set_ylabel("data / fit", fontsize=13)
        if row == 1:
            axr.set_xlabel(r"$M(\pi^+\pi^-e^+e^-)$ [GeV/$c^2$]")
        else:
            plt.setp(axr.get_xticklabels(), visible=False)
fig.suptitle(r"Three-scale global fit $N(M)=\alpha_1[\eta\to\pi^+\pi^-e^+e^- + \eta\to\pi^+\pi^-\gamma] + \alpha_2[\eta\to\pi^+\pi^-\pi^0] + "
             r"\alpha_3[\mathrm{rest}]$ to $N-$CB in 0.43$-$0.75 GeV/$c^2$ (grey band), typical cuts", fontsize=14.5, y=0.993)
for ext in ("pdf", "png"):
    fig.savefig(FIG + f"prevstyle_globalfit_E1.{ext}", dpi=110)
plt.close(fig)
P("written prevstyle_globalfit_E1")

# ================================================================ F3
MX = 0.002 * np.arange(25)       # lower edges 0-48 MeV (GeV/c2)
MXC = MX + 0.001
fig = plt.figure(figsize=(17, 10.5))
gs = GridSpec(2, 2, figure=fig, height_ratios=[3, 1.6], hspace=0.07, wspace=0.17, left=0.065, right=0.985, top=0.90, bottom=0.075)
ymx = max(np.max(R[s]["mee"]["W"][:25] + np.sqrt(R[s]["mee"]["Wvar"][:25])) for s in SELS) * 1.15
dmin = min(np.min(R[s]["mee"]["W"][:25] - R[s]["sb"]["k"] * R[s]["mee"]["SB"][:25]
                  - np.sqrt(R[s]["mee"]["Wvar"][:25] + R[s]["sb"]["k"] ** 2 * R[s]["mee"]["SBvar"][:25])) for s in SELS)
dmax = max(np.max(R[s]["mee"]["W"][:25] - R[s]["sb"]["k"] * R[s]["mee"]["SB"][:25]
                  + np.sqrt(R[s]["mee"]["Wvar"][:25] + R[s]["sb"]["k"] ** 2 * R[s]["mee"]["SBvar"][:25])) for s in SELS)
X17SC = {}
for col, s in enumerate(SELS):
    r = R[s]
    me = r["mee"]
    k = r["sb"]["k"]
    a1s = r["gf_raw"]["a"][0] * r["sint"]
    ax = fig.add_subplot(gs[0, col])
    axd = fig.add_subplot(gs[1, col], sharex=ax)
    for a_ in (ax, axd):
        a_.axvspan(0.0155, 0.0245, color="#fff2a8", zorder=0)
        a_.axvline(0.012, color="#b8860b", ls="--", lw=1.2, zorder=1)
        a_.axvline(0.022, color="#b8860b", ls="--", lw=1.2, zorder=1)
    xe = np.r_[MX, 0.050]

    def st(a_, y, **kw):
        a_.step(xe, np.r_[y[:25], y[24]], where="post", **kw)

    sbs = k * me["SB"][:25]
    sbe = k * np.sqrt(me["SBvar"][:25])
    ax.fill_between(xe, np.r_[sbs - sbe, (sbs - sbe)[-1]], np.r_[sbs + sbe, (sbs + sbe)[-1]], step="post", color=PAL["sb"], alpha=0.25, lw=0)
    st(ax, sbs, color=PAL["sb"], lw=1.8, label=f"side bands 518$-$530 + 557$-$569 MeV/$c^2$, $N-$CB $\\times$ {k:.3f}\n(scaled to the expected background in the window)")
    pk = a1s * me["peak"][0][:25]
    st(ax, sbs + pk, color="black", lw=2.0, label=r"scaled side bands + SMASH $\eta$ peak ($\alpha_1$ fit scale)")
    st(ax, a1s * me["C5"][0], color=PAL["C5"], lw=1.6, ls="--", label=CLAB["C5"] + r" (SMASH, $\times$" + f"{a1s:.3f})")
    st(ax, a1s * me["C4"][0], color=PAL["C4"], lw=1.6, ls="--", label=CLAB["C4"] + r" (SMASH, $\times$" + f"{a1s:.3f})")
    c6 = me["C6"][0][:25]
    # arbitrary normalisation, the same in both pads: maximum = 0.6 x the largest window - side band difference
    x17sc = 0.6 * np.max(me["W"][:25] - sbs) / np.max(c6)
    X17SC[s] = x17sc
    st(ax, x17sc * c6, color=PAL["X17"], lw=1.8, ls="-.", label="X17 (17 MeV/$c^2$) from SMASH, arbitrary normalisation")
    ax.errorbar(MXC, me["W"][:25], np.sqrt(me["Wvar"][:25]), fmt="o", color=PAL["sig"], ms=5, lw=1.0, label="data $N-$CB in the $\\eta$ peak window 530$-$557 MeV/$c^2$", zorder=5)
    ax.set_xlim(0, 0.050)
    ax.set_ylim(0, ymx)
    ax.set_ylabel(r"Counts per 2 MeV/$c^2$")
    plt.setp(ax.get_xticklabels(), visible=False)
    ax.set_title(SELTITLE[s], loc="left", fontsize=12, fontweight="bold")
    b = r["sb"]
    x0 = r["x17"]["15.5-24.5"]
    x1 = r["x17"]["12-22"]
    txt = (f"window: $N-$CB = {b['W']:.0f} $\\pm$ {b['We']:.0f};  side bands: {b['SB']:.0f} $\\pm$ {b['SBe']:.0f}\n"
           f"expected background in the window (global fit) = {b['BW']:.0f} $\\pm$ {b['BWe']:.0f}\n"
           f"15.5$-$24.5 MeV/$c^2$: window $-$ side bands = {x0['diff']:.0f} $\\pm$ {x0['diffe']:.0f}, SMASH $\\eta$ = {x0['peak']:.0f}\n"
           f"12$-$22 MeV/$c^2$: window $-$ side bands = {x1['diff']:.0f} $\\pm$ {x1['diffe']:.0f}, SMASH $\\eta$ = {x1['peak']:.0f}")
    ax.text(0.985, 0.97, txt, transform=ax.transAxes, ha="right", va="top", fontsize=11,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#999999", alpha=0.92), zorder=6)
    if col == 0:
        h, l = ax.get_legend_handles_labels()
        h += [Patch(color="#fff2a8"), Line2D([], [], color="#b8860b", ls="--")]
        l += ["X17 window of this analysis, 15.5$-$24.5 MeV/$c^2$", "X17 window of the earlier analysis, 12$-$22 MeV/$c^2$"]
        ax.legend(h, l, loc="upper right", bbox_to_anchor=(0.985, 0.74), fontsize=10, frameon=True, framealpha=0.92)
    dd = me["W"][:25] - sbs
    dde = np.sqrt(me["Wvar"][:25] + k ** 2 * me["SBvar"][:25])
    axd.axhline(0, color="k", lw=0.8)
    st(axd, pk, color="black", lw=1.8, label=r"SMASH $\eta$ peak ($\eta\to\pi^+\pi^-e^+e^-$ + $\eta\to\pi^+\pi^-\gamma(e^+e^-)$, $\times$" + f"{a1s:.3f})")
    st(axd, x17sc * c6, color=PAL["X17"], lw=1.6, ls="-.", label="X17 shape (arbitrary normalisation)")
    axd.legend(loc="upper right", bbox_to_anchor=(0.985, 0.83), fontsize=10, frameon=True, framealpha=0.92)
    axd.errorbar(MXC, dd, dde, fmt="o", color=PAL["sig"], ms=5, lw=1.0, zorder=5)
    axd.set_ylim(dmin * 1.15 if dmin < 0 else -0.05 * dmax, dmax * 1.15)
    axd.set_ylabel("window $-$ side bands")
    axd.set_xlabel(r"$M(e^+e^-)$ [GeV/$c^2$]")
    axd.text(0.985, 0.94, f"$\\chi^2$ (difference vs SMASH $\\eta$ peak, 0$-$50 MeV/$c^2$) = {r['mee_chi2'][0]:.1f}/{r['mee_chi2'][1]}",
             transform=axd.transAxes, ha="right", va="top", fontsize=11, bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#999999", alpha=0.9))
fig.suptitle(r"$M(e^+e^-)$ in the $\eta$ peak window 530$-$557 MeV/$c^2$ against the side bands, typical cuts, data $N-$CB", fontsize=14.5)
for ext in ("pdf", "png"):
    fig.savefig(FIG + f"prevstyle_mee_E1.{ext}", dpi=110)
plt.close(fig)
P("written prevstyle_mee_E1")


# ---------------------------------------------------------------- numbers
def clean(o):
    if isinstance(o, dict):
        return {k: clean(v) for k, v in o.items() if k not in ("f",)}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, np.ndarray):
        return clean(o.tolist())
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    return o


NUM = {}
for s in SELS:
    r = R[s]
    NUM[s] = dict(sint=r["sint"], sinte=r["sinte"], chi2a=r["chi2a"], cb=r["cb"], cbe=r["cbe"], cb40=r["cb40"], chi2b=r["chi2b"],
                  chi2b145=r["chi2b145"], reg=r["reg"], share=r["share"], data_all=r["data_all"], N_all=r["N_all"], CB_all=r["CB_all"],
                  T12_all=r["T12_all"], gf_raw=r["gf_raw"], gf_raw_dataonly=r["gf_raw_dataonly"], gf_sm=r["gf_sm"], gf_sm_wide=r["gf_sm_wide"],
                  smooth={g: dict(name=r["sm"][g]["name"], p=r["sm"][g]["p"], dev=r["sm"][g]["dev"], ndf=r["sm"][g]["ndf"], wbar=r["sm"][g]["wbar"],
                                  int_gf=r["sm"][g]["int_gf"]) for g in r["sm"]},
                  sb=r["sb"], x17=r["x17"], mee_chi2=r["mee_chi2"], mee_sum=r["mee_sum"], x17_arbitrary_scale=X17SC[s])
NUM["checks"] = CHK
NUM["scales"] = dict(F_BUMP=F_BUMP, A_P=A_P, G_ETA=G_ETA, C_B_E1=C_B_E1, FLEX=FLEX)
json.dump(clean(NUM), open(OUT + "prevstyle_numbers.json", "w"), indent=1)
P("written prevstyle_numbers.json")
LOG.close()
