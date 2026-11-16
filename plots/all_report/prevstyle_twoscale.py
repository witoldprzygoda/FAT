#!/usr/bin/env python3
# ============================================================================
# prevstyle_twoscale.py  --  package prevstyle of report_all, addition: two free eta scales on a background fixed on the slope
#                            DEMONSTRATION in the style of the earlier analysis (not a measurement; no result of the report uses it)
# ============================================================================
# Written 2026-11-16. New file; prevstyle_ana.py and its outputs are only read.
# Inputs : /hdd3/przygoda/HADES/FAT_out/all_report/prevstyle/prevstyle_data.json, prevstyle_sim.json (prevstyle_fill.C)
#          /hdd3/przygoda/HADES/FAT_out/all_report/prevstyle/prevstyle_numbers.json (check: c_b of display (b))
# Outputs: docs/latex/figs_all/prevstyle_twoscale_E1.pdf/.png
#          /hdd3/przygoda/HADES/FAT_out/all_report/prevstyle/prevstyle_twoscale_numbers.json, prevstyle_twoscale.log
# Method (analyst's request, 2026-11-16):
#   selection pre + E1 (alpha_lab < 50 deg, alpha_eta > 140 deg, M(pipi) < 0.46 GeV/c2, OA(e+e-) > 2 deg); "all" = isBest >= 0
#   (all stored combinations; NOMINAL, top row) and "best" = isBest == 1 (best combination only, as in the earlier analysis;
#   comparison, bottom row); data N - CB with CB = 2 sqrt(N++ N--) from the like-sign data (k = 1).
#   (1) background fixed on the slope right of the eta peak: B = c_b (G_D + G_C) + 0.656 G_eta, c_b such that B + the eta bump and
#       peak terms match the data in the normalisation window 0.575-0.70 GeV/c2 (start: eta scales 0.426 / 0.418, as display (b));
#   (2) with B fixed, two free scales fitted to N - CB in the fit range 0.43-0.575 GeV/c2 (5 MeV/c2 bins):
#       N - CB = B + f_bump (C2 + C3) + a_peak (C4 + C5)
#       C2 = eta -> pi+pi-pi0, pi0 -> e+e-gamma, with its per-count QED gamma* weights (c2recon); C3 = eta -> pi+pi-pi0,
#       pi0 -> gamma gamma, one photon converted; C4 = eta -> pi+pi-gamma, photon converted; C5 = eta -> pi+pi-e+e-.
#       chi2 with the variance data (N event groups + CB replicas) + SMASH statistics of every term (scaled), iterated.
#   (3) c_b recomputed once with the fitted eta scales in the window, fit repeated (a second iteration is checked as well).
#   Sensitivity: normalisation window (0.575-0.65, 0.60-0.70, 0.575-0.75), fit range (0.40-0.575, 0.45-0.575, 0.47-0.575,
#   0.43-0.56), and the background shape (linear tilt of G_D + G_C by +-3 % per 0.1 GeV/c2, pivot at the window centre).
#   Errors of the scales: fit (data + SMASH statistics) and the propagated statistical error of c_b, quoted separately.
#   Window sums: the nominal window 0.575-0.70 uses the exact window cells of prevstyle_fill.C (as display (b)); the other windows
#   sum the 5 MeV/c2 bins, with the N and CB variances scaled by the ratio exact/bin-sum found in 0.575-0.70 (correlations between
#   bins from several combinations per event and from the 10 MeV/c2 CB cells).
# Run: python3 -I plots/all_report/prevstyle_twoscale.py   (a few seconds)
# ============================================================================
import json
import math

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

OUT = "/hdd3/przygoda/HADES/FAT_out/all_report/prevstyle/"
FIG = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_all/"
L = 6079.2
PT3 = 64. / 63.
NB2RAW = L / PT3
F_BUMP, A_P, G_ETA = 0.426, 0.418, 0.656       # dedicated measurements (bump edge, peak side bands) and other eta decays
FLEX, FLEXE = 0.46, 0.06                       # common eta scale of the flexible description
DISP, DISPE = 0.443, 0.043                     # one-scale display fit (c) after the typical cuts, package e1 (E1)
SELS = ["all", "best"]
SELSHORT = {"all": "isBest >= 0 (nominal)", "best": "isBest == 1 (earlier analysis, comparison)"}
SELTITLE = {"all": r"isBest $\geq$ 0: all stored combinations (nominal)",
            "best": r"isBest == 1: best combination only, as in the earlier analysis (comparison)"}
PAL = dict(sig="#1f3fbf", bump="#ff6a00", peak="#00a5c8", bkg="#7a7a7a", tot="black", norm="#ffe9a8", fitr="#e3ecff")
plt.rcParams.update({"font.size": 13, "axes.labelsize": 15, "axes.titlesize": 13, "legend.fontsize": 11, "xtick.labelsize": 12,
                     "ytick.labelsize": 12, "pdf.fonttype": 42})
LOG = open(OUT + "prevstyle_twoscale.log", "w")


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    LOG.write(s + "\n")
    LOG.flush()


DJ = json.load(open(OUT + "prevstyle_data.json"))
SJ = json.load(open(OUT + "prevstyle_sim.json"))
NUMREF = json.load(open(OUT + "prevstyle_numbers.json"))
PARTS = {k: tuple(v) for k, v in DJ["meta"]["parts"].items()}
WIN = [w[0] for w in DJ["meta"]["windows"]]
WLOHI = {w[0]: (w[1], w[2]) for w in DJ["meta"]["windows"]}
iW = {w: i for i, w in enumerate(WIN)}


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
    v = 0.
    e = 0.
    for c in comps:
        x = SJ["sel"][s][c]
        v = v + part(x["s"], tag) * NB2RAW
        e = e + part(x["v"], tag) * NB2RAW ** 2
    return v, e


GROUPS = {"bump": ["C2", "C3"], "peak": ["C4", "C5"], "Geta": ["Geta"], "nonEta": ["GD", "GC"], "T12": ["T12"]}
XL = 0.30 + 0.005 * np.arange(160)
XC = XL + 0.0025
RNG = (XL >= 0.35 - 1e-9) & (XL < 0.90 - 1e-9)
M47 = (XC > 0.40) & (XC < 0.70)


def b5(a):
    return a[1:161]


def bins_in(lo, hi):
    return (XL >= lo - 1e-9) & (XL + 0.005 <= hi + 1e-9)


CHK = []


def check(what, val, ref, tol):
    ok = abs(val - ref) <= tol
    CHK.append(f"{what:<84s} {val:15.7g} ref {ref:15.7g}  {'OK' if ok else 'DIFFERS'}")
    return ok


# ---------------------------------------------------------------- per selection inputs
IN = {}
for s in SELS:
    d = {}
    d["dv"], d["dvar"] = (b5(x) for x in D(s, "b5"))
    d["N"], d["NV"] = (b5(x) for x in D(s, "b5", "N"))
    d["CB"], d["CBvar"] = (b5(x) for x in D(s, "b5", "CB"))
    d["comp"] = {g: tuple(b5(x) for x in S(s, cs, "b5")) for g, cs in GROUPS.items()}
    d["compw"] = {g: S(s, cs, "win") for g, cs in GROUPS.items()}
    d["dw"], d["dwv"] = D(s, "win")
    # variance ratio exact window / bin sum in 0.575-0.70 (N: event groups, CB: replicas), used for the other windows
    m = bins_in(0.575, 0.70)
    ic = iW["CR"]
    nvw = D(s, "win", "N")[1][ic]
    cbvw = D(s, "win", "CB")[1][ic]
    d["rN"] = float(nvw / d["NV"][m].sum())
    d["rCB"] = float(cbvw / d["CBvar"][m].sum())
    d["rS"] = {g: float(d["compw"][g][1][ic] / d["comp"][g][1][m].sum()) if d["comp"][g][1][m].sum() > 0 else 1.0 for g in GROUPS}
    check(f"{s}: data N - CB, sum of 5 MeV bins = window cell (0.575-0.70)", float(d["dv"][m].sum()), float(d["dw"][ic]), 1e-3)
    for g in ["nonEta", "Geta", "bump", "peak"]:
        check(f"{s}: SMASH {g}, sum of 5 MeV bins = window cell (0.575-0.70)", float(d["comp"][g][0][m].sum()), float(d["compw"][g][0][ic]),
              1e-6 * max(1., float(d["compw"][g][0][ic])))
    IN[s] = d


def window_sums(s, lo, hi, tilt=None):
    """data and SMASH sums (values, variances) in the normalisation window [lo, hi)"""
    d = IN[s]
    if (abs(lo - 0.575) < 1e-9 and abs(hi - 0.70) < 1e-9) and tilt is None:
        ic = iW["CR"]
        return (float(d["dw"][ic]), float(d["dwv"][ic]), {g: (float(d["compw"][g][0][ic]), float(d["compw"][g][1][ic])) for g in GROUPS})
    m = bins_in(lo, hi)
    dvv = float(d["NV"][m].sum() * d["rN"] + d["CBvar"][m].sum() * d["rCB"])
    comps = {}
    for g in GROUPS:
        v, e = d["comp"][g]
        if g == "nonEta" and tilt is not None:
            v, e = v * tilt, e * tilt ** 2
        comps[g] = (float(v[m].sum()), float(e[m].sum() * d["rS"][g]))
    return float(d["dv"][m].sum()), dvv, comps


def cb_from_window(s, lo, hi, fB, aP, tilt=None):
    dw, dwv, c = window_sums(s, lo, hi, tilt)
    eta = fB * c["bump"][0] + aP * c["peak"][0] + G_ETA * c["Geta"][0]
    etav = fB ** 2 * c["bump"][1] + aP ** 2 * c["peak"][1] + G_ETA ** 2 * c["Geta"][1]
    ne, nev = c["nonEta"]
    cb = (dw - eta) / ne
    cbe = math.sqrt(dwv + etav + cb ** 2 * nev) / ne
    return cb, cbe, dict(data=dw, eta=eta, nonEta=ne, eta_share=eta / dw)


def fit2(s, cb, lo, hi, tilt=None, it=200):
    """chi2 fit of N - CB - B(c_b) = f_bump bump + a_peak peak in [lo, hi); variance data + SMASH statistics (iterated)"""
    d = IN[s]
    ne, nev = d["comp"]["nonEta"]
    if tilt is not None:
        ne, nev = ne * tilt, nev * tilt ** 2
    ge, gev = d["comp"]["Geta"]
    B = cb * ne + G_ETA * ge
    Bv = cb ** 2 * nev + G_ETA ** 2 * gev
    T = [d["comp"]["bump"][0], d["comp"]["peak"][0]]
    Tv = [d["comp"]["bump"][1], d["comp"]["peak"][1]]
    y = d["dv"] - B
    m = bins_in(lo, hi) & (d["dvar"] > 0)
    a = np.array([F_BUMP, A_P])
    for _ in range(it):
        var = d["dvar"] + Bv + a[0] ** 2 * Tv[0] + a[1] ** 2 * Tv[1]
        w = 1. / var[m]
        A = np.array([[np.sum(T[i][m] * T[j][m] * w) for j in range(2)] for i in range(2)])
        b = np.array([np.sum(y[m] * T[i][m] * w) for i in range(2)])
        an = np.linalg.solve(A, b)
        if np.max(np.abs(an - a)) < 1e-11:
            a = an
            break
        a = an
    var = d["dvar"] + Bv + a[0] ** 2 * Tv[0] + a[1] ** 2 * Tv[1]
    w = 1. / var[m]
    A = np.array([[np.sum(T[i][m] * T[j][m] * w) for j in range(2)] for i in range(2)])
    cov = np.linalg.inv(A)
    mod = B + a[0] * T[0] + a[1] * T[1]
    modv = Bv + a[0] ** 2 * Tv[0] + a[1] ** 2 * Tv[1]
    chi2 = float(np.sum((d["dv"][m] - mod[m]) ** 2 * w))
    return dict(a=a, cov=cov, e=np.sqrt(np.diag(cov)), rho=float(cov[0, 1] / math.sqrt(cov[0, 0] * cov[1, 1])), chi2=chi2,
                ndf=int(m.sum()) - 2, mod=mod, modv=modv, B=B, Bv=Bv, nbins=int(m.sum()))


def fit2_dataonly(s, cb, lo, hi):
    """same fit with the data errors only (no SMASH statistics in the variance): shows how much of the error is template noise"""
    d = IN[s]
    B = cb * d["comp"]["nonEta"][0] + G_ETA * d["comp"]["Geta"][0]
    T = [d["comp"]["bump"][0], d["comp"]["peak"][0]]
    m = bins_in(lo, hi) & (d["dvar"] > 0)
    w = 1. / d["dvar"][m]
    y = d["dv"] - B
    A = np.array([[np.sum(T[i][m] * T[j][m] * w) for j in range(2)] for i in range(2)])
    b = np.array([np.sum(y[m] * T[i][m] * w) for i in range(2)])
    a = np.linalg.solve(A, b)
    cov = np.linalg.inv(A)
    mod = B + a[0] * T[0] + a[1] * T[1]
    return dict(a=a, e=np.sqrt(np.diag(cov)), chi2=float(np.sum((d["dv"][m] - mod[m]) ** 2 * w)), ndf=int(m.sum()) - 2)


def run(s, nlo=0.575, nhi=0.70, flo=0.43, fhi=0.575, tilt_pct=None):
    """full procedure: c_b with the measured eta scales -> fit -> c_b with the fitted scales -> fit; propagated c_b error"""
    tilt = None
    if tilt_pct is not None:
        x0 = 0.5 * (nlo + nhi)
        tilt = 1. + tilt_pct / 100. * (XC - x0) / 0.1
    cb0, cbe0, w0 = cb_from_window(s, nlo, nhi, F_BUMP, A_P, tilt)
    f0 = fit2(s, cb0, flo, fhi, tilt)
    cb1, cbe1, w1 = cb_from_window(s, nlo, nhi, f0["a"][0], f0["a"][1], tilt)
    f1 = fit2(s, cb1, flo, fhi, tilt)
    cb2, _, _ = cb_from_window(s, nlo, nhi, f1["a"][0], f1["a"][1], tilt)
    f2 = fit2(s, cb2, flo, fhi, tilt)
    # propagated statistical error of c_b (window data and SMASH statistics): refit with c_b +- 1 sigma
    fp = fit2(s, cb1 + cbe1, flo, fhi, tilt)
    fm = fit2(s, cb1 - cbe1, flo, fhi, tilt)
    ecb = 0.5 * np.abs(fp["a"] - fm["a"])
    # sensitivity per 1 % of the fixed background in the fit range (c_b x 1.01, no renormalisation)
    f101 = fit2(s, cb1 * 1.01, flo, fhi, tilt)
    return dict(cb0=cb0, cbe0=cbe0, cb1=cb1, cbe1=cbe1, cb2=cb2, f0=f0, f1=f1, f2=f2, ecb=ecb, d1pc=(f101["a"] - f1["a"]),
                win=w1, nlo=nlo, nhi=nhi, flo=flo, fhi=fhi, tilt=tilt_pct, tiltarr=tilt)


def model_windows(s, r):
    """data / model in the report regions (exact window cells; model = fixed background + fitted eta terms)"""
    d = IN[s]
    out = {}
    fB, aP = r["f1"]["a"]
    cb = r["cb1"]
    for w in ["below", "B", "SBlo", "P", "Wpk", "SBhi", "r56", "r70", "CR"]:
        i = iW[w]
        cw = d["compw"]
        mod = cb * cw["nonEta"][0][i] + G_ETA * cw["Geta"][0][i] + fB * cw["bump"][0][i] + aP * cw["peak"][0][i]
        bkg = cb * cw["nonEta"][0][i] + G_ETA * cw["Geta"][0][i]
        out[w] = dict(data=float(d["dw"][i]), err=float(math.sqrt(d["dwv"][i])), model=float(mod), ratio=float(d["dw"][i] / mod),
                      bkg=float(bkg), bump=float(fB * cw["bump"][0][i]), peak=float(aP * cw["peak"][0][i]),
                      bkg_over_bump=float(bkg / (fB * cw["bump"][0][i])) if cw["bump"][0][i] > 0 else None,
                      bkg_over_peak=float(bkg / (aP * cw["peak"][0][i])) if cw["peak"][0][i] > 0 else None)
    return out


# ---------------------------------------------------------------- nominal and variants
VARS = [("nominal", dict()),
        ("norm 0.575-0.65", dict(nlo=0.575, nhi=0.65)),
        ("norm 0.60-0.70", dict(nlo=0.60, nhi=0.70)),
        ("norm 0.575-0.75", dict(nlo=0.575, nhi=0.75)),
        ("fit 0.40-0.575", dict(flo=0.40)),
        ("fit 0.45-0.575", dict(flo=0.45)),
        ("fit 0.47-0.575", dict(flo=0.47)),
        ("fit 0.43-0.56", dict(fhi=0.56)),
        ("tilt +3%/0.1", dict(tilt_pct=3.0)),
        ("tilt -3%/0.1", dict(tilt_pct=-3.0))]
R = {}
for s in SELS:
    R[s] = {}
    for name, kw in VARS:
        R[s][name] = run(s, **kw)
    R[s]["nominal"]["reg"] = model_windows(s, R[s]["nominal"])
    # chi2 in 0.40-0.70 of the nominal model (three adjusted numbers: c_b, f_bump, a_peak)
    r = R[s]["nominal"]
    d = IN[s]
    m = M47 & (d["dvar"] > 0)
    var = d["dvar"] + r["f1"]["modv"]
    r["chi2_4070"] = (float(np.sum((d["dv"][m] - r["f1"]["mod"][m]) ** 2 / var[m])), int(m.sum()) - 3)
    r["dataonly"] = fit2_dataonly(s, r["cb1"], 0.43, 0.575)
    check(f"{s}: start c_b (eta 0.426 / 0.418 in 0.575-0.70) = c_b of display (b)", r["cb0"], NUMREF[s]["cb"], 1e-9)
    check(f"{s}: its error = error of display (b)", r["cbe0"], NUMREF[s]["cbe"], 1e-9)

for c in CHK:
    P(c)
NBAD = sum("DIFFERS" in c for c in CHK)
P(f"checks: {len(CHK)} done, {NBAD} differ")
P(f"variance ratios exact/bin-sum in 0.575-0.70: " + ", ".join(f"{s}: N {IN[s]['rN']:.3f}, CB {IN[s]['rCB']:.3f}" for s in SELS))

P("\n==== two free eta scales on a background fixed on the slope (DEMONSTRATION, not a measurement) ====")
for s in SELS:
    P(f"\n[{SELSHORT[s]}]")
    for name, _ in VARS:
        r = R[s][name]
        f1 = r["f1"]
        P(f"  {name:18s} c_b start {r['cb0']:.4f} -> iterated {r['cb1']:.4f} +- {r['cbe1']:.4f} (2nd iteration {r['cb2']:.5f}); "
          f"f_bump = {f1['a'][0]:.3f} +- {f1['e'][0]:.3f} (fit) +- {r['ecb'][0]:.3f} (c_b); "
          f"a_peak = {f1['a'][1]:.3f} +- {f1['e'][1]:.3f} (fit) +- {r['ecb'][1]:.3f} (c_b); rho {f1['rho']:+.2f}; "
          f"chi2/ndf = {f1['chi2']:.1f}/{f1['ndf']} = {f1['chi2'] / f1['ndf']:.2f}")
        P(f"  {'':18s} before the c_b iteration: f_bump {r['f0']['a'][0]:.4f}, a_peak {r['f0']['a'][1]:.4f}; after the 2nd: {r['f2']['a'][0]:.4f}, "
          f"{r['f2']['a'][1]:.4f}; eta share of the window data {r['win']['eta_share']:.4f}; per +1 % of the fixed background: "
          f"d f_bump {r['d1pc'][0]:+.4f}, d a_peak {r['d1pc'][1]:+.4f}")
    r = R[s]["nominal"]
    P(f"  chi2 of the nominal model in 0.40-0.70 (60 bins, 3 adjusted numbers): {r['chi2_4070'][0]:.1f}/{r['chi2_4070'][1]}")
    P(f"  c_b iteration (nominal): {r['cb0']:.6f} -> {r['cb1']:.6f} (change {r['cb1'] - r['cb0']:+.6f}) -> {r['cb2']:.6f}; scales "
      f"{r['f0']['a'][0]:.5f}, {r['f0']['a'][1]:.5f} -> {r['f1']['a'][0]:.5f}, {r['f1']['a'][1]:.5f}")
    do = r["dataonly"]
    P(f"  nominal with data errors only: f_bump {do['a'][0]:.3f} +- {do['e'][0]:.3f}, a_peak {do['a'][1]:.3f} +- {do['e'][1]:.3f}, "
      f"chi2/ndf {do['chi2']:.1f}/{do['ndf']} = {do['chi2'] / do['ndf']:.2f}")
    tb = [R[s][n]["f1"]["a"][0] for n in ("tilt +3%/0.1", "tilt -3%/0.1")]
    tp = [R[s][n]["f1"]["a"][1] for n in ("tilt +3%/0.1", "tilt -3%/0.1")]
    P(f"  background tilt +-3 % per 0.1 GeV/c2: f_bump {min(tb):.3f}-{max(tb):.3f}, a_peak {min(tp):.3f}-{max(tp):.3f}; "
      f"tilt factor at 0.50 GeV/c2: {1 + 0.03 * (0.50 - 0.6375) / 0.1:.4f} / {1 - 0.03 * (0.50 - 0.6375) / 0.1:.4f}")
    for w, x in r["reg"].items():
        lo, hi = WLOHI[w]
        P(f"    {w:6s} {lo:.4f}-{hi:.4f}: data N-CB {x['data']:10.1f} +- {x['err']:6.1f}  model {x['model']:10.1f}  data/model {x['ratio']:.3f}  "
          f"(fixed bkg {x['bkg']:.1f}, bump {x['bump']:.1f}, peak {x['peak']:.1f}; bkg/bump {x['bkg_over_bump'] if x['bkg_over_bump'] is None else round(x['bkg_over_bump'], 2)}, "
          f"bkg/peak {x['bkg_over_peak'] if x['bkg_over_peak'] is None else round(x['bkg_over_peak'], 2)})")
    nb = [R[s][n]["f1"]["a"][0] for n, _ in VARS if n.startswith("norm") or n == "nominal"]
    fr = [R[s][n]["f1"]["a"][0] for n, _ in VARS if n.startswith("fit") or n == "nominal"]
    nbp = [R[s][n]["f1"]["a"][1] for n, _ in VARS if n.startswith("norm") or n == "nominal"]
    frp = [R[s][n]["f1"]["a"][1] for n, _ in VARS if n.startswith("fit") or n == "nominal"]
    P(f"  ranges: normalisation window -> f_bump {min(nb):.3f}-{max(nb):.3f}, a_peak {min(nbp):.3f}-{max(nbp):.3f}; "
      f"fit range -> f_bump {min(fr):.3f}-{max(fr):.3f}, a_peak {min(frp):.3f}-{max(frp):.3f}")


# ================================================================ figure
def steps(ax, y, color, lw=1.6, ls="-", label=None, rng=RNG, zorder=3):
    xl = XL[rng]
    xe = np.r_[xl, xl[-1] + 0.005]
    yy = np.r_[y[rng], y[rng][-1]]
    ax.step(xe, yy, where="post", color=color, lw=lw, ls=ls, label=label, zorder=zorder)


FH = 13.6                                   # figure height [in]; header positions in inches from the top
fig = plt.figure(figsize=(17, FH))
TOP, BOT = 1 - 2.72 / FH, 0.70 / FH
gl = GridSpec(2, 1, figure=fig, hspace=0.19, left=0.058, right=0.640, top=TOP, bottom=BOT)
gr = GridSpec(2, 1, figure=fig, hspace=0.19, left=0.795, right=0.985, top=TOP, bottom=BOT)
ymax = max(np.max((IN[s]["dv"] + np.sqrt(IN[s]["dvar"]))[RNG]) for s in SELS) * 1.68
VLAB = {"nominal": "nominal: window\n0.575$-$0.70, fit 0.43$-$0.575",
        "norm 0.575-0.65": "window 0.575$-$0.65", "norm 0.60-0.70": "window 0.60$-$0.70", "norm 0.575-0.75": "window 0.575$-$0.75",
        "fit 0.40-0.575": "fit range 0.40$-$0.575", "fit 0.45-0.575": "fit range 0.45$-$0.575", "fit 0.47-0.575": "fit range 0.47$-$0.575",
        "fit 0.43-0.56": "fit range 0.43$-$0.56",
        "tilt +3%/0.1": "background shape tilted\n+3 % per 0.1 GeV/$c^2$", "tilt -3%/0.1": "background shape tilted\n$-$3 % per 0.1 GeV/$c^2$"}
for row, s in enumerate(SELS):
    r = R[s]["nominal"]
    d = IN[s]
    f1 = r["f1"]
    fB, aP = f1["a"]
    inner = GridSpecFromSubplotSpec(2, 1, subplot_spec=gl[row, 0], height_ratios=[3, 1.1], hspace=0.06)
    ax = fig.add_subplot(inner[0])
    axr = fig.add_subplot(inner[1], sharex=ax)
    for a_ in (ax, axr):
        a_.axvspan(0.575, 0.70, color=PAL["norm"], zorder=0, lw=0)
        a_.axvspan(0.43, 0.575, color=PAL["fitr"], zorder=0, lw=0)
    m = RNG
    ax.errorbar(XC[m], d["dv"][m], np.sqrt(d["dvar"][m]), fmt="o", color=PAL["sig"], ms=4.2, lw=0.9, label="data $N-$CB (CB from the like-sign data)",
                zorder=5)
    B = f1["B"]
    steps(ax, B, PAL["bkg"], lw=2.0, ls="--", label=f"fixed background: non-$\\eta$ SMASH $\\times$ $c_b$ + other $\\eta$ decays $\\times$ {G_ETA}")
    steps(ax, fB * d["comp"]["bump"][0], PAL["bump"], lw=2.0, ls="--",
          label=r"fitted $\eta$ Dalitz bump: $\eta\to\pi^+\pi^-\pi^0$, $\pi^0\to e^+e^-\gamma$ or $\gamma\gamma$ ($\gamma$ converted)")
    steps(ax, aP * d["comp"]["peak"][0], PAL["peak"], lw=2.0, ls="--",
          label=r"fitted $\eta$ peak: $\eta\to\pi^+\pi^-e^+e^-$ + $\eta\to\pi^+\pi^-\gamma$ ($\gamma$ converted)")
    steps(ax, f1["mod"], PAL["tot"], lw=2.0, label="sum: background + bump + peak")
    ax.set_xlim(0.35, 0.90)
    ax.set_ylim(0, ymax)
    ax.set_ylabel(r"Counts per 5 MeV/$c^2$")
    plt.setp(ax.get_xticklabels(), visible=False)
    ax.set_title(SELTITLE[s], loc="left", fontsize=13, fontweight="bold")
    txt = (f"scales relative to SMASH:\n"
           f"$\\eta$ Dalitz bump  {fB:.3f} $\\pm$ {f1['e'][0]:.3f} (fit) $\\pm$ {r['ecb'][0]:.3f} ($c_b$)\n"
           f"$\\eta$ peak            {aP:.3f} $\\pm$ {f1['e'][1]:.3f} (fit) $\\pm$ {r['ecb'][1]:.3f} ($c_b$)\n"
           f"non-$\\eta$ $c_b$ = {r['cb1']:.3f} $\\pm$ {r['cbe1']:.3f} (from 0.575$-$0.70)\n"
           f"$\\chi^2$/ndf (0.43$-$0.575) = {f1['chi2']:.1f}/{f1['ndf']} = {f1['chi2'] / f1['ndf']:.2f}\n"
           f"our measurements: bump {F_BUMP} (kinematic\nedge), peak {A_P} (side bands)")
    ax.text(0.985, 0.975, txt, transform=ax.transAxes, ha="right", va="top", fontsize=11.5, linespacing=1.25,
            bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#999999", alpha=0.94), zorder=6)
    if row == 0:
        hl, ll = ax.get_legend_handles_labels()
        hl += [Patch(color=PAL["norm"]), Patch(color=PAL["fitr"])]
        ll += ["normalisation window of the background", "fit range of the two $\\eta$ scales"]
        order = [4, 0, 1, 2, 3, 5, 6] if len(hl) == 7 else list(range(len(hl)))
        fig.legend([hl[i] for i in order], [ll[i] for i in order], loc="upper left", bbox_to_anchor=(0.045, 1 - 1.17 / FH), ncol=2, fontsize=11.5,
                   frameon=False, handlelength=2.6, labelspacing=0.35, columnspacing=1.6)
    # ratio pad (bins where the model is above 2 % of its maximum)
    mv, mvv = f1["mod"], f1["modv"]
    ok = RNG & (mv > 0.02 * np.max(mv[RNG]))
    xb = XL[ok]
    band = np.sqrt(mvv[ok]) / mv[ok]
    axr.fill_between(np.r_[xb, xb[-1] + 0.005], np.r_[1 - band, (1 - band)[-1]], np.r_[1 + band, (1 + band)[-1]], step="post",
                     color="#9a9a9a", alpha=0.45, lw=0, zorder=1)
    axr.axhline(1, color="k", lw=0.8, zorder=2)
    axr.errorbar(XC[ok], d["dv"][ok] / mv[ok], np.sqrt(d["dvar"][ok]) / mv[ok], fmt="o", color=PAL["sig"], ms=3.2, lw=0.9, zorder=5)
    axr.set_ylim(0.75, 1.25)
    axr.set_ylabel("data / model", fontsize=13)
    if row == 1:
        axr.set_xlabel(r"$M(\pi^+\pi^-e^+e^-)$ [GeV/$c^2$]")
    else:
        plt.setp(axr.get_xticklabels(), visible=False)
    # sensitivity panel
    axs = fig.add_subplot(gr[row, 0])
    axs.axvspan(FLEX - FLEXE, FLEX + FLEXE, color="#dddddd", zorder=0, lw=0)
    axs.axvspan(DISP - DISPE, DISP + DISPE, facecolor="none", edgecolor="#8c8c8c", hatch="///", zorder=0, lw=0)
    axs.axvline(F_BUMP, color=PAL["bump"], ls=":", lw=2.2, zorder=1)
    axs.axvline(A_P, color=PAL["peak"], ls=":", lw=2.2, zorder=1)
    names = [n for n, _ in VARS]
    ys = np.arange(len(names))[::-1].astype(float)
    for y, n in zip(ys, names):
        rr = R[s][n]
        ff = rr["f1"]
        eb = math.hypot(ff["e"][0], rr["ecb"][0])
        ep = math.hypot(ff["e"][1], rr["ecb"][1])
        axs.errorbar(ff["a"][0], y + 0.16, xerr=eb, fmt="o", color=PAL["bump"], ms=6.5, lw=1.4, capsize=2.5, zorder=4)
        axs.errorbar(ff["a"][1], y - 0.16, xerr=ep, fmt="s", color=PAL["peak"], ms=6.0, lw=1.4, capsize=2.5, zorder=4)
        if n == "nominal":
            axs.axhspan(y - 0.47, y + 0.47, color="#e3ecff", zorder=0, lw=0)
    for ysep in (ys[0] - 0.5, ys[3] - 0.5, ys[7] - 0.5):
        axs.axhline(ysep, color="#bbbbbb", lw=0.8)
    axs.set_yticks(ys)
    axs.set_yticklabels([VLAB[n] for n in names], fontsize=11)
    axs.set_ylim(ys[-1] - 0.6, ys[0] + 0.6)
    axs.set_xlim(0.2, 1.0)
    axs.set_xlabel("scale relative to SMASH")
    axs.set_title("sensitivity of the two scales (GeV/$c^2$)", fontsize=12.5, loc="left")
    axs.tick_params(axis="y", length=0, pad=4)
    if row == 0:
        hs = [Line2D([], [], color=PAL["bump"], marker="o", ls="", ms=6.5), Line2D([], [], color=PAL["peak"], marker="s", ls="", ms=6),
              Line2D([], [], color=PAL["bump"], ls=":", lw=2.2), Line2D([], [], color=PAL["peak"], ls=":", lw=2.2),
              Patch(color="#dddddd"), Patch(facecolor="none", edgecolor="#8c8c8c", hatch="///")]
        hs = [hs[0], hs[1], hs[4], hs[2], hs[3], hs[5]]
        ls_ = ["this fit: $\\eta$ Dalitz bump", "this fit: $\\eta$ peak", f"flexible description: {FLEX} $\\pm$ {FLEXE}",
               f"bump, kinematic edge: {F_BUMP}", f"peak, side bands: {A_P}", f"one-scale display fit: {DISP:.2f} $\\pm$ {DISPE:.2f}"]
        fig.legend(hs, ls_, loc="upper left", bbox_to_anchor=(0.640, 1 - 1.17 / FH), ncol=2, fontsize=10.5, frameon=False, handlelength=1.9,
                   labelspacing=0.35, columnspacing=0.9)
fig.text(0.5, 1 - 0.12 / FH, "Demonstration in the style of the earlier analysis (not a measurement)", ha="center", va="top", fontsize=20,
         fontweight="bold", color="#9b1c1c")
fig.text(0.5, 1 - 0.57 / FH, r"Background (non-$\eta$ SMASH $\times$ $c_b$ + other $\eta$ decays $\times$ 0.656) fixed on the slope right of the $\eta$ "
         r"peak, 0.575$-$0.70 GeV/$c^2$; then two free $\eta$ scales fitted to $N-$CB in 0.43$-$0.575 GeV/$c^2$" "\n"
         r"typical cuts: $\alpha_{lab}<50^\circ$, $\alpha_\eta>140^\circ$, $M(\pi^+\pi^-)<0.46$ GeV/$c^2$, OA$(e^+e^-)>2^\circ$; "
         r"errors: data and SMASH statistics; no result of the report uses these scales",
         ha="center", va="top", fontsize=13, linespacing=1.4)
for ext in ("pdf", "png"):
    fig.savefig(FIG + f"prevstyle_twoscale_E1.{ext}", dpi=110)
plt.close(fig)
P("written prevstyle_twoscale_E1")


# ---------------------------------------------------------------- numbers
def clean(o):
    if isinstance(o, dict):
        return {k: clean(v) for k, v in o.items() if k not in ("mod", "modv", "B", "Bv", "tiltarr", "cov")}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, np.ndarray):
        return clean(o.tolist())
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    return o


NUM = {s: {n: R[s][n] for n, _ in VARS} for s in SELS}
for s in SELS:
    NUM[s]["nominal_reg"] = R[s]["nominal"]["reg"]
    NUM[s]["chi2_4070"] = R[s]["nominal"]["chi2_4070"]
NUM["checks"] = CHK
NUM["scales"] = dict(F_BUMP=F_BUMP, A_P=A_P, G_ETA=G_ETA, FLEX=FLEX, FLEXE=FLEXE, DISP=DISP, DISPE=DISPE)
NUM["note"] = "DEMONSTRATION in the style of the earlier analysis; not a measurement; no result of the report uses it"
json.dump(clean(NUM), open(OUT + "prevstyle_twoscale_numbers.json", "w"), indent=1)
P("written prevstyle_twoscale_numbers.json")
LOG.close()
