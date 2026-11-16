#!/usr/bin/env python3
# ============================================================================
# e1_figfix.py  --  COPY of e1_ana.py (built by the writer with scratchpad/final/e1/make_figfix.py) that redraws the e1 figures
# used in report_all with the verifier's corrections (e1v): X17 pair-momentum median from the upper bin edges (424, not 414 MeV/c),
# display fit drawn only over its range 0.40-0.70, en dash in the panel title, no ratio points where the model is almost empty,
# no vertical drop of the dashed X17 line, no line across the pre | S0 divider. Figures: figs_all/<name>_v2.pdf/.png;
# tables go to /hdd3/przygoda/HADES/FAT_out/all_report/e1/figfix/ (the e1 tables are not touched). Numbers are unchanged.
# Run: python3 -I plots/all_report/e1_figfix.py
# Original header of e1_ana.py follows.
# e1_ana.py  --  package e1 of report_all: data against SMASH after the "typical" eta cuts (E1-E4) and the X17-oriented
# OA < 6 deg cuts (X1-X4); figures, tables and numbers (prefix e1)
# ============================================================================
# Inputs: /hdd3/przygoda/HADES/FAT_out/all_report/e1/e1_data.json, e1_sim.json (plots/all_report/e1_fill.C)
#         /home/przygoda/HADES/FAT/FAT_simwt/plots/eta_sim/tables/rw_compare_enriched.txt (earlier study, read only)
# Conventions (Part III): data (N - CB) x 64/63 / L, L = 6079.2 nb^-1, CB = 2 sqrt(N++ N--) (k = 1, 10 MeV cells, spread rule);
#   SMASH sigma x PDG BR x V3 (mode A), M4' = 0.99615 M4, M(pipi)' = 0.9952 M(pipi); expected raw rows = nb x L x 63/64.
#   eta channels = C2 + C3 + C4 + C5 + G_eta; non-eta background = G_pi0D + G_conv + G_vec + G_oth.
#   (a) SMASH as is; (b) eta x 0.40 and the non-eta background x c (c from data in 0.575-0.70 GeV/c2, one factor per selection);
#   (c) two-parameter DISPLAY fit in 0.40-0.70 (5 MeV bins): data = a eta + c bkg (not a yield).
#   X17: C6 per unit BR; the 90 % CLs limit of S0 is N_X < 392 raw rows in the X2D M4 range (0.5225-0.61, all M(ee));
#   in a selection: N_X(sel, cell) = 392 x C6(sel, cell) / C6(S0, X2D).
# Outputs: docs/latex/figs_all/e1_*.pdf (+ .png), docs/latex/tables_all/e1_*.txt (+ .tex)
# Run: python3 -I plots/all_report/e1_ana.py   (about 1 min)
# ============================================================================
import json
import math
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

E1DIR = "/hdd3/przygoda/HADES/FAT_out/all_report/e1/"
FIG = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_all/"
TAB = "/hdd3/przygoda/HADES/FAT_out/all_report/e1/figfix/"
import os as _os
_os.makedirs(TAB, exist_ok=True)
OLDTAB = "/home/przygoda/HADES/FAT/FAT_simwt/plots/eta_sim/tables/rw_compare_enriched.txt"
L = 6079.2
PT3 = 64. / 63.
NB2RAW = L / PT3            # nb -> expected raw rows
RAW2NB = PT3 / L
LOLD = 6080.5
ETA_SCALE = 0.40
NX_UL = 392.                # 90 % CLs limit on X17 rows, S0, X2D M4 range, 17 MeV/c2 (x17d)
LREL = 0.069
PAL = dict(N="black", CB="red", NmCB="#1f3fbf", GD="#d0d0d0", GC="#707070", Geta="#c9a3f5", C2="#ff6a00", C3="#d400d4",
           C4="#009900", C5="#00cccc", X17="#cc0077", fit="#8c4b00")
plt.rcParams.update({"font.size": 13, "axes.labelsize": 15, "axes.titlesize": 14, "legend.fontsize": 11, "xtick.labelsize": 12,
                     "ytick.labelsize": 12, "pdf.fonttype": 42})

DJ = json.load(open(E1DIR + "e1_data.json"))
SJ = json.load(open(E1DIR + "e1_sim.json"))
PARTS = {k: tuple(v) for k, v in DJ["meta"]["parts"].items()}
WIN = [w[0] for w in DJ["meta"]["windows"]]
WTXT = {w[0]: (w[1], w[2], w[3]) for w in DJ["meta"]["windows"]}
SELS = list(DJ["sel"].keys())
assert SELS == list(SJ["sel"].keys())

LAB = {"pre": "pre", "pre_E1": "E1", "pre_E2": "E2", "pre_E3": "E3", "pre_E3oa": "E3oa", "pre_E4": "E4", "pre_E4oa": "E4oa",
       "pre_X1": "X1", "pre_X2": "X2", "pre_O6": "X4", "pre_X2old": "X2old"}
for s in SELS:
    if s.startswith("S0"):
        LAB[s] = "S0" if s == "S0" else ("X3" if s == "S0_O6" else "S0+" + s[3:])
CUTDEF = {"": "", "E1": r"$\alpha_{lab}<50^\circ$, $\alpha_\eta>140^\circ$, $M_{\pi\pi}<0.46$, OA$>2^\circ$",
          "E2": r"$\alpha_{lab}<50^\circ$, $\alpha_\eta>140^\circ$, $M_{\pi\pi}<0.46$",
          "E3": r"$\alpha_{lab}<50^\circ$, $\alpha_\eta>140^\circ$, $M_{\pi\pi}<0.42$",
          "E3oa": r"$\alpha_{lab}<50^\circ$, $\alpha_\eta>140^\circ$, $M_{\pi\pi}<0.42$, OA$>2^\circ$",
          "E4": r"$\alpha_{lab}<50^\circ$, $\alpha_\eta>140^\circ$", "E4oa": r"$\alpha_{lab}<50^\circ$, $\alpha_\eta>140^\circ$, OA$>2^\circ$",
          "X1": r"$\alpha_{lab}<50^\circ$, $\alpha_\eta>140^\circ$, $M_{\pi\pi}<0.46$, OA$<6^\circ$",
          "X2": r"$\alpha_{lab}<50^\circ$, $\alpha_\eta>140^\circ$, $M_{\pi\pi}<0.46$, $2^\circ<$OA$<6^\circ$",
          "O6": r"OA$<6^\circ$", "X2old": r"$\alpha_{lab}<50^\circ$, $\alpha_\eta>140^\circ$, $M_{\pi\pi}<0.42$, $2^\circ<$OA$<6^\circ$"}
CUTTXT = {"": "", "E1": "oa_lab < 50, oa_eta > 140, M(pipi) < 0.46, OA > 2", "E2": "oa_lab < 50, oa_eta > 140, M(pipi) < 0.46",
          "E3": "oa_lab < 50, oa_eta > 140, M(pipi) < 0.42", "E3oa": "oa_lab < 50, oa_eta > 140, M(pipi) < 0.42, OA > 2",
          "E4": "oa_lab < 50, oa_eta > 140", "E4oa": "oa_lab < 50, oa_eta > 140, OA > 2",
          "X1": "oa_lab < 50, oa_eta > 140, M(pipi) < 0.46, OA < 6", "X2": "oa_lab < 50, oa_eta > 140, M(pipi) < 0.46, 2 < OA < 6",
          "O6": "OA < 6", "X2old": "oa_lab < 50, oa_eta > 140, M(pipi) < 0.42, 2 < OA < 6"}


def base_cut(s):
    b = "S0" if s.startswith("S0") else "pre"
    c = s[len(b) + 1:] if len(s) > len(b) else ""
    return b, c


def sel_title(s):
    b, c = base_cut(s)
    bt = "pre" if b == "pre" else "S0"
    return f"{LAB[s]}: {bt}" + (" + " + CUTDEF[c] if c else "")


def sel_text(s):
    b, c = base_cut(s)
    return f"{LAB[s]} = {b}" + (" + " + CUTTXT[c] if c else "")


def part(arr, tag):
    o, n = PARTS[tag]
    return np.array(arr[o:o + n], dtype=float)


# ---------------------------------------------------------------- accessors
def D(s, tag, what="sig"):
    """data in part tag: 'N', 'CB', 'sig' (N - CB raw); returns (value, variance) arrays; CB variance = replica RMS^2"""
    d = DJ["sel"][s]
    N, NV = part(d["N"], tag), part(d["NV"], tag)
    CB, CBe = part(d["CB"], tag), part(d["CBrep"], tag)
    if what == "N":
        return N, NV
    if what == "CB":
        return CB, CBe ** 2
    return N - CB, NV + CBe ** 2


def Dnb(s, tag):
    v, e = D(s, tag)
    return v * RAW2NB, e * RAW2NB ** 2


def S(s, comps, tag):
    """simulation nb (sum of components), variance (event-group; components treated as independent, which they are row-wise)"""
    if isinstance(comps, str):
        comps = [comps]
    v = 0.
    e = 0.
    for c in comps:
        x = SJ["sel"][s][c]
        v = v + part(x["s"], tag)
        e = e + part(x["v"], tag)
    return v, e


ETA = ["C2", "C3", "C4", "C5", "Geta"]
BKG = ["Gpi0D", "Gconv", "Gvec", "Goth"]
iW = {w: i for i, w in enumerate(WIN)}


def wsel(arr, w):
    return arr[iW[w]]


# ---------------------------------------------------------------- checks of the totals
chk = []
def check(what, val, ref, tol):
    ok = abs(val - ref) <= tol
    chk.append(f"{what:<70s} {val:14.6g} ref {ref:14.6g}  {'OK' if ok else 'DIFFERS'}")
    return ok

T12w, _ = S("S0", "T12", "win")
sumw = sum(S("S0", c, "win")[0] for c in ETA + BKG)
check("S0: T12 = eta + non-eta in P (nb)", wsel(sumw, "P"), wsel(T12w, "P"), 1e-6)
check("S0: T12 in P = wp1_B1 7.9961 nb", wsel(T12w, "P"), 7.9961, 6e-5)
check("pre: T12 in P = wp1_B1 9.5901 nb", wsel(S("pre", "T12", "win")[0], "P"), 9.5901, 6e-5)
check("S0: data N in P = 55947", wsel(D("S0", "win", "N")[0], "P"), 55947, 0.5)
check("S0: data CB in P = 3279.0", wsel(D("S0", "win", "CB")[0], "P"), 3279.0, 0.06)
check("S0: data (N-CB) x 64/63/L in P = 8.8012 nb", wsel(Dnb("S0", "win")[0], "P"), 8.8012, 6e-5)
check("pre: data N in P = 77513", wsel(D("pre", "win", "N")[0], "P"), 77513, 0.5)
check("pre: data CB in P = 7948.8", wsel(D("pre", "win", "CB")[0], "P"), 7948.8, 0.06)
C6X2D_S0 = wsel(S("S0", "C6", "win")[0], "X2D")
check("S0: C6 per unit BR in the X2D M4 range = 6257.5 nb (x17d)", C6X2D_S0, 6257.5, 0.06)
check("pre: T12 in R5 = wp1_B1 108.4174 nb", wsel(S("pre", "T12", "win")[0], "R5"), 108.4174, 6e-4)

# ---------------------------------------------------------------- per-selection analysis
REG = [("below", "below the bump 0.35-0.455"), ("B", "bump window 0.455-0.5225"), ("P", "peak window 0.5336-0.5542"),
       ("r56", "0.56-0.70"), ("r70", "0.70-1.00")]
o5, n5 = PARTS["b5"]
X5 = 0.30 + 0.005 * (np.arange(1, 161) - 0.5)          # centres of b5 bins 1..160
FITM = (X5 > 0.40) & (X5 < 0.70)


def b5(arr):
    return arr[1:161]


def ratio(a, ea2, b, eb2):
    r = a / b
    return r, abs(r) * math.sqrt(ea2 / a ** 2 + eb2 / b ** 2) if a != 0 and b != 0 else float("nan")


def display_fit(dv, dvar, ev, evar, bv, bvar, mask):
    a, c = 0.4, 1.0
    for _ in range(60):
        var = dvar + a * a * evar + c * c * bvar
        m = mask & (var > 0)
        w = 1. / var[m]
        Saa = np.sum(ev[m] ** 2 * w); Scc = np.sum(bv[m] ** 2 * w); Sac = np.sum(ev[m] * bv[m] * w)
        Sya = np.sum(dv[m] * ev[m] * w); Syc = np.sum(dv[m] * bv[m] * w)
        det = Saa * Scc - Sac ** 2
        a, c = (Sya * Scc - Syc * Sac) / det, (Saa * Syc - Sac * Sya) / det
    var = dvar + a * a * evar + c * c * bvar
    m = mask & (var > 0)
    chi2 = float(np.sum((dv[m] - a * ev[m] - c * bv[m]) ** 2 / var[m]))
    return dict(a=a, ea=math.sqrt(Scc / det), c=c, ec=math.sqrt(Saa / det), rho=-Sac / det / math.sqrt(Scc / det * Saa / det),
                chi2=chi2, ndf=int(m.sum()) - 2)


def chi2_model(dv, dvar, mv, mvar, mask, npar):
    var = dvar + mvar
    m = mask & (var > 0)
    return float(np.sum((dv[m] - mv[m]) ** 2 / var[m])), int(m.sum()) - npar


def mismatch(dv, dvar, mv, mvar, mask, npar):
    """fractional model mismatch delta per bin such that sum r^2/(var + (delta m)^2) = ndf (0 if chi2/ndf <= 1 already)"""
    m = mask & ((dvar + mvar) > 0)
    r2 = (dv[m] - mv[m]) ** 2; var = dvar[m] + mvar[m]; mm = mv[m] ** 2; ndf = int(m.sum()) - npar
    f = lambda d: float(np.sum(r2 / (var + d * d * mm))) - ndf
    if f(0.) <= 0:
        return 0.
    lo, hi = 0., 1.
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if f(mid) > 0 else (lo, mid)
    return 0.5 * (lo + hi)


R = {}
for s in SELS:
    r = {}
    dw, dwv = Dnb(s, "win")
    ew, ewv = S(s, ETA, "win")
    bw, bwv = S(s, BKG, "win")
    tw, twv = S(s, "T12", "win")
    # (b) background factor from 0.575-0.70
    cr = iW["CR"]
    cb = (dw[cr] - ETA_SCALE * ew[cr]) / bw[cr]
    cbe = math.sqrt(dwv[cr] + ETA_SCALE ** 2 * ewv[cr] + cb ** 2 * bwv[cr]) / bw[cr]
    r["cb"], r["cbe"] = cb, cbe
    # region ratios
    r["reg"] = {}
    Nw, _ = D(s, "win", "N"); CBw, CBv = D(s, "win", "CB"); NVw = part(DJ["sel"][s]["NV"], "win")
    for w, _t in REG + [("R2", ""), ("post", ""), ("Fit", ""), ("all", ""), ("CR", ""), ("PB", "")]:
        i = iW[w]
        ma, mav = tw[i], twv[i]
        mb, mbv = ETA_SCALE * ew[i] + cb * bw[i], ETA_SCALE ** 2 * ewv[i] + cb ** 2 * bwv[i]
        r["reg"][w] = dict(N=Nw[i], CB=CBw[i], CBe=math.sqrt(CBv[i]), sig=Nw[i] - CBw[i], sigE=math.sqrt(NVw[i] + CBv[i]),
                           dnb=dw[i], dnbE=math.sqrt(dwv[i]), sim=ma, simE=math.sqrt(mav), eta=ew[i], bkg=bw[i],
                           ra=ratio(dw[i], dwv[i], ma, mav), rb=ratio(dw[i], dwv[i], mb, mbv), mb=mb)
    # b5 spectra in 0.40-0.70: chi2 of (a), (b) and the display fit (c)
    dv, dvar = (b5(x) for x in Dnb(s, "b5"))
    ev, evar = (b5(x) for x in S(s, ETA, "b5"))
    bv, bvar = (b5(x) for x in S(s, BKG, "b5"))
    tv, tvar = (b5(x) for x in S(s, "T12", "b5"))
    r["chi2a"] = chi2_model(dv, dvar, tv, tvar, FITM, 0)
    r["chi2b"] = chi2_model(dv, dvar, ETA_SCALE * ev + cb * bv, ETA_SCALE ** 2 * evar + cb ** 2 * bvar, FITM, 1)
    r["fit"] = display_fit(dv, dvar, ev, evar, bv, bvar, FITM)
    ft_ = r["fit"]
    r["delta"] = dict(a=mismatch(dv, dvar, tv, tvar, FITM, 0),
                      b=mismatch(dv, dvar, ETA_SCALE * ev + cb * bv, ETA_SCALE ** 2 * evar + cb ** 2 * bvar, FITM, 1),
                      c=mismatch(dv, dvar, ft_["a"] * ev + ft_["c"] * bv, ft_["a"] ** 2 * evar + ft_["c"] ** 2 * bvar, FITM, 2))
    # one common scale (data = k x SMASH) for reference
    m = FITM
    w_ = 1. / (dvar[m] + tvar[m])
    k1 = np.sum(dv[m] * tv[m] * w_) / np.sum(tv[m] ** 2 * w_)
    r["k1"] = k1
    r["chi2k1"] = chi2_model(dv, dvar, k1 * tv, k1 ** 2 * tvar, FITM, 1)
    # X17
    x17N, x17NV = D(s, "x17", "N")
    x17CB, x17CBv = D(s, "x17", "CB")
    c6x = S(s, "C6", "x17")[0][0]
    c6pb = wsel(S(s, "C6", "win")[0], "PB")
    r["x17"] = dict(N=x17N[0], NV=x17NV[0], CB=x17CB[0], CBe=math.sqrt(x17CBv[0]), c6=c6x, c6pb=c6pb,
                    NX=NX_UL * c6x / C6X2D_S0, NXpb=NX_UL * c6pb / C6X2D_S0)
    # model (b) in the X17 window (raw rows)
    ex = S(s, ETA, "x17")[0][0]; bx = S(s, BKG, "x17")[0][0]; c5x = S(s, "C5", "x17")[0][0]; c4x = S(s, "C4", "x17")[0][0]
    r["x17"]["model_b"] = (ETA_SCALE * ex + cb * bx) * NB2RAW
    r["x17"]["C5b"] = ETA_SCALE * c5x * NB2RAW
    r["x17"]["C4b"] = ETA_SCALE * c4x * NB2RAW
    # M(ee) shape in the peak band: data N / (CB + model (b)) in M(ee) 2-14, the X17 window and 26-62 MeV/c2
    Nm = D(s, "mee", "N")[0]; CBm = D(s, "mee", "CB")[0]
    mod = CBm + (ETA_SCALE * S(s, ETA, "mee")[0] + cb * S(s, BKG, "mee")[0]) * NB2RAW
    r["x17"]["rat_lo"] = float(np.sum(Nm[1:7]) / np.sum(mod[1:7]))
    r["x17"]["rat_hi"] = float(np.sum(Nm[13:31]) / np.sum(mod[13:31]))
    r["x17"]["rat_win"] = float(x17N[0] / (x17CB[0] + r["x17"]["model_b"]))
    r["x17"]["NX_share"] = float(r["x17"]["NX"] / (x17CB[0] + r["x17"]["model_b"]))
    R[s] = r
# relative X17 quantities
for s in SELS:
    x, x0 = R[s]["x17"], R["S0"]["x17"]
    x["eps_rel"] = x["c6"] / x0["c6"]
    x["eps_rel_pb"] = x["c6pb"] / x0["c6pb"]
    x["B_rel"] = x["N"] / x0["N"]
    x["proxy"] = math.sqrt(x["N"]) / x["eps_rel"] / math.sqrt(x0["N"])
    vs, v0 = x["N"] + x["CBe"] ** 2, x0["N"] + x0["CBe"] ** 2
    x["proxy_cb"] = math.sqrt(vs) / x["eps_rel"] / math.sqrt(v0)
    x["Z"] = x["NX"] / math.sqrt(vs)
    x["SB_gain"] = x["eps_rel"] / x["B_rel"]

# ---------------------------------------------------------------- X17 opening angle (S0, peak band)
def oa_hist(s, comp, tag):
    return S(s, comp, tag)[0]


oaC6 = oa_hist("S0", "C6", "oa"); oaC6x = oa_hist("S0", "C6", "oaX")
oaC5 = oa_hist("S0", "C5", "oa"); oaC5x = oa_hist("S0", "C5", "oaX")
oaD, oaDv = D("S0", "oa"); oaDx, oaDxv = D("S0", "oaX")


def frac(h, lo, hi):
    return float(np.sum(h[lo:hi]) / np.sum(h))


OAF = {}
for name, h in [("C6 peak band", oaC6), ("C6 X17 window", oaC6x), ("C5 peak band", oaC5), ("C5 X17 window", oaC5x),
                ("data N-CB peak band", oaD), ("data N-CB X17 window", oaDx)]:
    OAF[name] = dict(lt2=frac(h, 0, 2), lt6=frac(h, 0, 6), b26=frac(h, 2, 6), gt2=frac(h, 2, len(h)), gt6=frac(h, 6, len(h)),
                     mean=float(np.sum(h[:40] * (np.arange(40) + 0.5)) / np.sum(h[:40])), over40=frac(h, 40, 41))
peeC6 = S("S0", "C6", "pee")[0]; peeC5 = S("S0", "C5", "pee")[0]
PCEN = 20. * (np.arange(101) + 0.5)
p_thr = 2. * 17.0 / math.radians(6.)       # p_ee above which a symmetric 17 MeV pair has OA < 6 deg
p_thr2 = 2. * 17.0 / math.radians(2.)
PEDG = 20. * np.arange(102)
def _share_above(h, thr):
    lo = PEDG[:-1]; hi = PEDG[1:]
    full = np.sum(h[lo >= thr]); k = int(np.searchsorted(hi, thr, side='right'))
    part = h[k] * (hi[k] - thr) / (hi[k] - lo[k]) if lo[k] < thr < hi[k] else 0.
    return float((full + part) / np.sum(h))
peeF = dict(gt325=_share_above(peeC6, p_thr), median=float(np.interp(0.5, np.r_[0., np.cumsum(peeC6) / np.sum(peeC6)], PEDG)),
            median_old=float(np.interp(0.5, np.cumsum(peeC6) / np.sum(peeC6), PCEN)), gt325_old=float(np.sum(peeC6[PCEN > p_thr]) / np.sum(peeC6)))
print('FIGFIX pee C6: median (upper edges) %.1f, share above %.1f = %.4f; old median %.1f, old share %.4f' % (
      peeF['median'], p_thr, peeF['gt325'], peeF['median_old'], peeF['gt325_old']))

# ---------------------------------------------------------------- reproduction of the earlier study (old convention)
OLDREG = {"below": "below bump 0.35-0.455", "bump": "bump window 0.455-0.5225", "gap": "gap 0.5225-0.5336",
          "peak": "peak window 0.5336-0.5542", "post": "0.5542-0.56", "r56": "0.56-0.70", "r70": "0.70-1.00", "all": "all 0.35-1.00"}
oldref = {}
sec1 = False
for line in open(OLDTAB):
    if line.startswith("## 1."):
        sec1 = True; continue
    if line.startswith("## ") and sec1:
        break
    if not sec1:
        continue
    m = re.match(r"^(\S+)\s+(.+?)\s+([-\d.]+)\+-([\d.]+)\s+([\d.]+) -> ([\d.]+)\s+\|\s+([\d.]+) \+- ([\d.]+) \(([\d.]+)\)\s+([\d.]+) \+- ([\d.]+)", line)
    if m:
        oldref[(m.group(1), m.group(2).strip())] = dict(d=float(m.group(3)), de=float(m.group(4)), sb=float(m.group(5)), sa=float(m.group(6)),
                                                        rb=float(m.group(7)), ra=float(m.group(10)), rae=float(m.group(11)))
REPRO = []
for o in ["pre", "E1", "E2", "E3", "E3oa", "E4", "E4oa"]:
    for k, name in OLDREG.items():
        dn, dcb, dsig, dsigE = DJ["old"][o][k]
        sv, se = SJ["old"][o][k]
        d, de = dsig / LOLD, dsigE / LOLD
        ra = d / sv
        rae = ra * math.hypot(de / d, se / sv)
        ref = oldref.get((o, name))
        REPRO.append((o, name, d, de, sv, ra, rae, ref))

# ================================================================ TABLES
def fmt_ve(v, e, nd=3):
    return f"{v:.{nd}f} +- {e:.{nd}f}"


with open(TAB + "e1_repro.txt", "w") as f:
    f.write("# e1_repro.txt -- plots/all_report/e1_ana.py: reproduction of the earlier study tables/rw_compare_enriched.txt (section 1,\n")
    f.write("#   'after' = with the physics weights V3) in ITS convention: data (N - CB) / 6080.5 (no 64/63), CB on the 0.35-1.00 grid\n")
    f.write("#   (scan::FillMany, 10 MeV, reb 4); SMASH T1 + T2 at the UNSCALED M4 and the raw M(pipi), w x factor 108, mode A.\n")
    f.write("#   ref = the earlier table (printed digits); diff columns: this - ref.\n")
    f.write(f"{'sel':6s} {'region':28s} {'data N-CB [nb]':>18s} {'ref':>8s} {'sim [nb]':>9s} {'ref':>8s} {'data/sim':>16s} {'ref':>7s}  verdict\n")
    nbad = 0
    for o, name, d, de, sv, ra, rae, ref in REPRO:
        if ref:
            ok = abs(d - ref["d"]) < 6e-5 and abs(sv - ref["sa"]) < 6e-4 and abs(ra - ref["ra"]) < 6e-4
            nbad += 0 if ok else 1
            f.write(f"{o:6s} {name:28s} {d:9.4f}+-{de:7.4f} {ref['d']:8.4f} {sv:9.3f} {ref['sa']:8.3f} {ra:7.3f} +- {rae:5.3f} {ref['ra']:7.3f}  "
                    f"{'OK' if ok else 'DIFFERS'}\n")
        else:
            f.write(f"{o:6s} {name:28s} {d:9.4f}+-{de:7.4f} {'--':>8s} {sv:9.3f} {'--':>8s} {ra:7.3f} +- {rae:5.3f} {'--':>7s}\n")
    f.write(f"# rows that differ from the earlier table beyond the last printed digit: {nbad}\n")
    f.write("\n# checks of the Part III convention against Part III tables\n")
    for c in chk:
        f.write("# " + c + "\n")
NBAD_REPRO = nbad

ORDER_PRE = ["pre", "pre_E1", "pre_E2", "pre_E3", "pre_E3oa", "pre_E4", "pre_E4oa", "pre_X1", "pre_X2", "pre_X2old", "pre_O6"]
ORDER_S0 = ["S0", "S0_E1", "S0_E2", "S0_E3", "S0_E3oa", "S0_E4", "S0_E4oa", "S0_X1", "S0_X2", "S0_X2old", "S0_O6"]
ORDER = ORDER_PRE + ORDER_S0

with open(TAB + "e1_table.txt", "w") as f:
    f.write("# e1_table.txt -- plots/all_report/e1_ana.py (package e1): data against SMASH after the typical eta cuts and OA < 6 deg\n")
    f.write("# data: (N - CB) x 64/63 / L, L = 6079.2 nb^-1, CB = 2 sqrt(N++ N--) (k = 1, 10 MeV cells; replica errors); isBest >= 0\n")
    f.write("# SMASH: sigma x PDG BR x V3 (mode A), M4' = 0.99615 M4, M(pipi)' = 0.9952 M(pipi), T12 = eta + non-eta;\n")
    f.write("#   eta = C2 + C3 + C4 + C5 + G_eta, non-eta = G_pi0D + G_conv + G_vec + G_oth\n")
    f.write("# (a) SMASH as is; (b) eta x 0.40 and non-eta x c_b (c_b from the data in 0.575-0.70); (c) DISPLAY fit 0.40-0.70, 5 MeV bins,\n")
    f.write("#   data = a eta + c non-eta (NOT a yield; target-material rows and the CB at small OA are absorbed by c)\n")
    f.write("# errors: data stat (event groups) (+) CB replicas (+) SMASH MC; the 6.9 % of L is common and not included\n")
    f.write("# selections:\n")
    for s in ORDER:
        f.write(f"#   {sel_text(s)}\n")
    f.write("\n## 1. rows kept (raw rows): data N, CB (replica error), N - CB per region\n")
    f.write(f"{'sel':8s}" + "".join(f" | {t:>30s}" for _, t in REG) + "\n")
    for s in ORDER:
        f.write(f"{LAB[s]:8s}")
        for w, _ in REG:
            g = R[s]["reg"][w]
            f.write(f" | {g['N']:9.0f} {g['CB']:8.1f}+-{g['CBe']:5.1f} {g['sig']:8.0f}")
        f.write("\n")
    f.write("\n## 2. data / SMASH per region: (a) as is | (b) scaled\n")
    f.write(f"{'sel':8s}" + "".join(f" | {t[:22]:>27s}" for _, t in REG) + "\n")
    for s in ORDER:
        f.write(f"{LAB[s]:8s}")
        for w, _ in REG:
            g = R[s]["reg"][w]
            f.write(f" | {g['ra'][0]:5.3f}({1000 * g['ra'][1]:3.0f}) {g['rb'][0]:6.3f}({1000 * g['rb'][1]:3.0f})")
        f.write("\n")
    f.write("# numbers in brackets: error in units of 0.001\n")
    f.write("\n## 2b. eta share of the model in the windows: eta / (eta + non-eta), SMASH as is | scaled (b); last column: C5 x 0.4 / model (b) in P\n")
    f.write(f"{'sel':8s} | {'B as is':>8s} {'B (b)':>7s} | {'P as is':>8s} {'P (b)':>7s} | {'C5 share in P (b)':>17s}\n")
    for s in ORDER:
        g = R[s]["reg"]; cb_ = R[s]["cb"]
        def sh(w, f):
            e, b = g[w]["eta"], g[w]["bkg"]
            return f * e / (f * e + (cb_ if f != 1. else 1.) * b)
        c5p = wsel(S(s, "C5", "win")[0], "P") * ETA_SCALE / (ETA_SCALE * g["P"]["eta"] + cb_ * g["P"]["bkg"])
        f.write(f"{LAB[s]:8s} | {sh('B', 1.):8.3f} {sh('B', ETA_SCALE):7.3f} | {sh('P', 1.):8.3f} {sh('P', ETA_SCALE):7.3f} | {c5p:17.3f}\n")
    f.write("\n## 3. scales and goodness of fit in 0.40-0.70 (5 MeV bins)\n")
    f.write(f"{'sel':8s} | {'c_b (0.575-0.70)':>17s} | {'chi2/ndf (a)':>16s} | {'chi2/ndf (b)':>16s} | {'one scale k: k, chi2/ndf':>26s} | "
            f"{'fit (c): a_eta':>16s} {'c_bkg':>15s} {'rho':>6s} {'chi2/ndf':>15s} | {'delta (a) (b) (c)':>18s}\n")
    for s in ORDER:
        r = R[s]; ft = r["fit"]
        f.write(f"{LAB[s]:8s} | {r['cb']:7.3f} +- {r['cbe']:5.3f} | {r['chi2a'][0]:7.1f}/{r['chi2a'][1]:2d}={r['chi2a'][0] / r['chi2a'][1]:5.2f} | "
                f"{r['chi2b'][0]:7.1f}/{r['chi2b'][1]:2d}={r['chi2b'][0] / r['chi2b'][1]:5.2f} | {r['k1']:6.3f} {r['chi2k1'][0] / r['chi2k1'][1]:6.2f}"
                f"{'':12s} | {ft['a']:6.3f} +- {ft['ea']:5.3f} {ft['c']:6.3f} +- {ft['ec']:5.3f} {ft['rho']:6.2f} "
                f"{ft['chi2']:6.1f}/{ft['ndf']:2d}={ft['chi2'] / ft['ndf']:4.2f} | {r['delta']['a']:6.3f} {r['delta']['b']:5.3f} {r['delta']['c']:5.3f}\n")
    f.write("# delta = fractional model mismatch per 5 MeV bin that must be added (x model, in quadrature) to the errors for chi2/ndf = 1;\n")
    f.write("#   unlike chi2/ndf it does not fall when a cut removes rows (larger statistical errors); 0 = described within the errors\n")
    f.write("\n## 4. X17 (17 MeV/c2): window = eta peak band M4 0.5325-0.5575 and M(ee) 15.5-24.5 MeV/c2\n")
    f.write("# eps_rel = C6(sel, window) / C6(S0, window) (simulation); N_X = X17 rows in the window at the S0 limit, 392 x C6(sel, window) /\n")
    f.write("#   C6(S0, X2D M4 range) (= 6257.5 nb per unit BR); B = data N in the window; proxy = sqrt(B_sel) / eps_rel / sqrt(B_S0) (counting\n")
    f.write("#   proxy, NOT the frozen likelihood; smaller is better; proxy_cb adds var(CB)); Z = N_X / sqrt(N + var CB): an X17 at the present\n")
    f.write("#   S0 limit seen as a counting excess in this window; model (b) = eta x 0.40 + non-eta x c_b in the window (raw rows, no CB);\n")
    f.write("#   S/B gain = eps_rel / (B/B_S0): the change of signal over background (matters when the background systematics dominate)\n")
    f.write(f"{'sel':8s} | {'eps_rel':>7s} {'eps_rel(PB)':>11s} | {'N':>6s} {'CB':>12s} {'N-CB':>7s} {'model(b)':>8s} {'C5(b)':>6s} {'C4(b)':>6s} | "
            f"{'N_X':>6s} {'Z':>5s} | {'B/B_S0':>6s} {'proxy':>6s} {'proxy_cb':>8s} {'S/B gain':>8s}\n")
    for s in ORDER:
        x = R[s]["x17"]
        f.write(f"{LAB[s]:8s} | {x['eps_rel']:7.3f} {x['eps_rel_pb']:11.3f} | {x['N']:6.0f} {x['CB']:6.1f}+-{x['CBe']:4.1f} {x['N'] - x['CB']:7.0f} "
                f"{x['model_b']:8.0f} {x['C5b']:6.0f} {x['C4b']:6.0f} | {x['NX']:6.1f} {x['Z']:5.2f} | {x['B_rel']:6.3f} {x['proxy']:6.3f} {x['proxy_cb']:8.3f} {x['SB_gain']:8.3f}\n")
    f.write("\n## 4b. M(ee) shape in the peak band: data N / (CB + model (b)); NX_share = X17 rows at the S0 limit / (CB + model (b)) in the window\n")
    f.write(f"{'sel':8s} | {'2-14 MeV':>8s} {'15.5-24.5':>9s} {'26-62 MeV':>9s} | {'NX_share':>8s}\n")
    for s in ORDER:
        x = R[s]["x17"]
        f.write(f"{LAB[s]:8s} | {x['rat_lo']:8.3f} {x['rat_win']:9.3f} {x['rat_hi']:9.3f} | {x['NX_share']:8.3f}\n")
    f.write("\n## 5. opening angle of the simulated X17 (C6, 17.000 MeV/c2) and of the direct decays (C5), S0, peak band M4' 0.5325-0.5575\n")
    f.write(f"{'distribution':24s} {'OA<2':>6s} {'2-6':>6s} {'OA<6':>6s} {'OA>6':>6s} {'mean(0-40)':>10s} {'OA>=40':>7s}\n")
    for k, v in OAF.items():
        f.write(f"{k:24s} {v['lt2']:6.3f} {v['b26']:6.3f} {v['lt6']:6.3f} {v['gt6']:6.3f} {v['mean']:10.2f} {v['over40']:7.3f}\n")
    f.write(f"# X17 pair momentum (S0, peak band): median {peeF['median']:.0f} MeV/c; share above 2 m_X / (6 deg) = {p_thr:.0f} MeV/c "
            f"(the momentum a symmetric 17 MeV/c2 pair needs for OA < 6 deg): {peeF['gt325']:.3f}\n")

# ---------------------------------------------------------------- tex snippets
def tex_lab(s):
    return LAB[s].replace("+", "\\,+\\,")


with open(TAB + "e1_tableR.tex", "w") as f:
    f.write("% e1_tableR.tex -- plots/all_report/e1_ana.py: data/SMASH per M4 region, (a) as is / (b) eta x 0.40, non-eta x c_b\n")
    f.write("\\begin{tabular}{@{}l r r r r r r r r r r r@{}}\n\\toprule\n")
    f.write("& \\multicolumn{2}{c}{0.35--0.455} & \\multicolumn{2}{c}{B} & \\multicolumn{2}{c}{P} & \\multicolumn{2}{c}{0.56--0.70} & "
            "\\multicolumn{2}{c}{0.70--1.00} & \\\\\n")
    f.write("selection & (a) & (b) & (a) & (b) & (a) & (b) & (a) & (b) & (a) & (b) & $c_b$ \\\\\n\\midrule\n")
    for s in ORDER:
        if s == "S0":
            f.write("\\midrule\n")
        r = R[s]
        row = [tex_lab(s)]
        for w, _ in REG:
            g = r["reg"][w]
            row += [f"{g['ra'][0]:.2f}", f"{g['rb'][0]:.2f}"]
        row.append(f"{r['cb']:.3f}")
        f.write(" & ".join(row) + " \\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n")

with open(TAB + "e1_tableF.tex", "w") as f:
    f.write("% e1_tableF.tex -- plots/all_report/e1_ana.py: goodness of the description in 0.40-0.70 and the display-fit scales\n")
    f.write("\\begin{tabular}{@{}l r r r r r r r r r@{}}\n\\toprule\n")
    f.write("selection & $\\chi^2/n$ (a) & $\\delta$ (a) & $c_b$ & $\\chi^2/n$ (b) & $\\delta$ (b) & $a_\\eta$ (c) & $c$ (c) & $\\chi^2/n$ (c) & rows $N$ in 0.40--0.70 \\\\\n\\midrule\n")
    for s in ORDER:
        if s == "S0":
            f.write("\\midrule\n")
        r = R[s]; ft = r["fit"]
        f.write(f"{tex_lab(s)} & {r['chi2a'][0] / r['chi2a'][1]:.1f} & {r['delta']['a']:.3f} & {r['cb']:.3f} & {r['chi2b'][0] / r['chi2b'][1]:.1f} & {r['delta']['b']:.3f} & "
                f"${ft['a']:.3f}\\pm{ft['ea']:.3f}$ & ${ft['c']:.3f}\\pm{ft['ec']:.3f}$ & {ft['chi2'] / ft['ndf']:.1f} & {r['reg']['Fit']['N']:.0f} \\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n")

with open(TAB + "e1_tableX.tex", "w") as f:
    f.write("% e1_tableX.tex -- plots/all_report/e1_ana.py: X17 window (M4 0.5325-0.5575, M(ee) 15.5-24.5 MeV/c2), counting proxy\n")
    f.write("\\begin{tabular}{@{}l r r r r r r r r@{}}\n\\toprule\n")
    f.write("selection & $\\varepsilon_X/\\varepsilon_X(\\mathrm{S0})$ & $N$ & CB & $N-\\mathrm{CB}$ & $N_X$ at the S0 limit & $Z$ & proxy & S/B gain \\\\\n\\midrule\n")
    for s in ORDER:
        if s == "S0":
            f.write("\\midrule\n")
        x = R[s]["x17"]
        f.write(f"{tex_lab(s)} & {x['eps_rel']:.3f} & {x['N']:.0f} & ${x['CB']:.0f}\\pm{x['CBe']:.0f}$ & {x['N'] - x['CB']:.0f} & {x['NX']:.0f} & "
                f"{x['Z']:.2f} & {x['proxy']:.2f} & {x['SB_gain']:.2f} \\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n")

# ================================================================ FIGURES
def save(fig, name):
    assert not _os.path.exists(FIG + name + "_v2.pdf"), name + "_v2.pdf exists"
    fig.savefig(FIG + name + "_v2.pdf")
    fig.savefig(FIG + name + "_v2.png", dpi=110)
    plt.close(fig)


def stairs_stack(ax, edges, layers, labels=None):
    bottom = np.zeros(len(edges) - 1)
    for k, (h, col) in enumerate(layers):
        ax.stairs(bottom + h, edges, baseline=bottom, fill=True, color=col, lw=0, label=labels[k] if labels else None)
        bottom = bottom + h
    return bottom


F10E = np.round(np.concatenate([0.400 + 0.005 * np.arange(25), 0.520 + 0.0025 * np.arange(1, 21), 0.570 + 0.005 * np.arange(1, 27)]), 6)


def f10(vals, vvals, fine, fvar):
    """F10 binning from b5 (5 MeV bins, exact variance) and the fine (2.5 MeV) bins in 0.52-0.57"""
    out, outv = [], []
    for lo, hi in zip(F10E[:-1], F10E[1:]):
        if hi - lo > 0.004:
            x = int(round((lo - 0.30) / 0.005)) + 1
            out.append(vals[x]); outv.append(vvals[x])
        else:
            i = int(round((lo - 0.30) / 0.0025)) + 1
            out.append(fine[i]); outv.append(fvar[i])
    return np.array(out), np.array(outv)


def spectra(s, lay):
    """dict of arrays (nb per bin) for layout 'log' (b5 0.35-1.00) or 'lin' (F10)"""
    def get(fn):
        a5, v5 = fn("b5"); af, vf = fn("fine")
        if lay == "log":
            sl = slice(int(round((0.35 - 0.30) / 0.005)) + 1, int(round((1.00 - 0.30) / 0.005)) + 1)
            return a5[sl], v5[sl]
        return f10(a5, v5, af, vf)
    out = {}
    out["N"] = get(lambda t: tuple(x * RAW2NB ** (1 + k) for k, x in enumerate(D(s, t, "N"))))
    out["CB"] = get(lambda t: tuple(x * RAW2NB ** (1 + k) for k, x in enumerate(D(s, t, "CB"))))
    out["sig"] = get(lambda t: Dnb(s, t))
    for c in ["C2", "C3", "C4", "C5", "Geta", "Gpi0D", "Gconv", "Gvec", "Goth", "T12", "GD"]:
        out[c] = get(lambda t, c=c: S(s, c, t))
    out["eta"] = get(lambda t: S(s, ETA, t))
    out["bkg"] = get(lambda t: S(s, BKG, t))
    edges = np.round(0.35 + 0.005 * np.arange(131), 6) if lay == "log" else F10E
    return edges, out


def m4_figure(s):
    r = R[s]
    cb, ft = r["cb"], r["fit"]
    fig = plt.figure(figsize=(15, 12.5))
    gs = fig.add_gridspec(5, 2, height_ratios=[3, 1.05, 0.42, 3, 1.05], hspace=0.06, wspace=0.18, left=0.07, right=0.985, top=0.885, bottom=0.055)
    for row, lay in enumerate(["log", "lin"]):
        edges, H = spectra(s, lay)
        wdt = np.diff(edges)
        for col in range(2):
            ax = fig.add_subplot(gs[3 * row, col]); ar = fig.add_subplot(gs[3 * row + 1, col], sharex=ax)
            fe, fb = (1., 1.) if col == 0 else (ETA_SCALE, cb)
            layers = [(H["GD"][0] * fb / wdt, PAL["GD"]), (H["Gconv"][0] * fb / wdt, PAL["GC"]), (H["Geta"][0] * fe / wdt, PAL["Geta"]),
                      (H["C2"][0] * fe / wdt, PAL["C2"]), (H["C3"][0] * fe / wdt, PAL["C3"]), (H["C4"][0] * fe / wdt, PAL["C4"]),
                      (H["C5"][0] * fe / wdt, PAL["C5"])]
            tot = stairs_stack(ax, edges, layers)
            totv = (fe ** 2 * H["eta"][1] + fb ** 2 * H["bkg"][1])
            xc = 0.5 * (edges[1:] + edges[:-1])
            if col == 1:
                fitm = (ft["a"] * H["eta"][0] + ft["c"] * H["bkg"][0]) / wdt
                fr = (edges[:-1] >= 0.40 - 1e-9) & (edges[1:] <= 0.70 + 1e-9)
                fe_ = np.r_[edges[:-1][fr], edges[1:][fr][-1]]
                ax.stairs(fitm[fr], fe_, color=PAL["fit"], lw=1.8, ls="--", baseline=None)
            if lay == "log":
                ax.errorbar(xc, H["N"][0] / wdt, yerr=np.sqrt(H["N"][1]) / wdt, fmt="o", ms=2.6, color=PAL["N"], lw=0.8)
                ax.errorbar(xc, H["CB"][0] / wdt, yerr=np.sqrt(H["CB"][1]) / wdt, fmt="s", ms=2.2, color=PAL["CB"], lw=0.8)
            ax.errorbar(xc, H["sig"][0] / wdt, yerr=np.sqrt(H["sig"][1]) / wdt, fmt="o", ms=3.0, color=PAL["NmCB"], lw=0.9)
            for xw in (0.455, 0.5225, 0.53364, 0.55416):
                ax.axvline(xw, color="#999999", ls=":", lw=0.9)
            if lay == "log":
                ax.set_yscale("log")
                ymax = max(np.max(H["N"][0] / wdt), np.max(tot)) * 40
                pos = H["CB"][0][H["CB"][0] > 0] / wdt[H["CB"][0] > 0]
                ymin = max(np.min(pos) * 0.5 if len(pos) else 1e-2, ymax * 3e-5)
                ax.set_ylim(ymin, ymax)
                ax.set_xlim(0.35, 1.00)
            else:
                ymax = max(np.max(H["sig"][0] / wdt), np.max(tot)) * 1.32
                ax.set_ylim(0, ymax)
                ax.set_xlim(0.40, 0.70)
            ax.set_ylabel(r"d$\sigma$/d$M_4$ [nb/(GeV/$c^2$)]")
            plt.setp(ax.get_xticklabels(), visible=False)
            ttl = "(a) SMASH as is: $\\sigma\\times$BR$\\times$V3" if col == 0 else \
                f"(b) $\\eta\\times{ETA_SCALE:.2f}$, non-$\\eta\\times c_b$, $c_b={cb:.3f}$ (0.575\u20130.70)"
            ax.set_title(ttl, fontsize=13, loc="left")
            # text: region ratios
            g = r["reg"]
            key = "ra" if col == 0 else "rb"
            txt = (f"data/model: B {g['B'][key][0]:.2f}, P {g['P'][key][0]:.2f}, 0.56-0.70 {g['r56'][key][0]:.2f}, 0.70-1.00 {g['r70'][key][0]:.2f}\n"
                   + (f"$\\chi^2/n$(0.40-0.70) = {r['chi2a'][0]:.0f}/{r['chi2a'][1]}" if col == 0 else
                      f"$\\chi^2/n$ = {r['chi2b'][0]:.0f}/{r['chi2b'][1]}; display fit (dashed): $a_\\eta$={ft['a']:.2f}$\\pm${ft['ea']:.2f}, "
                      f"$c$={ft['c']:.3f}, $\\chi^2/n$={ft['chi2']:.0f}/{ft['ndf']}"))
            if lay == "log" or True:
                ax.text(0.015, 0.97, txt, transform=ax.transAxes, va="top", ha="left", fontsize=10.5,
                        bbox=dict(facecolor="white", alpha=0.8, edgecolor="none", pad=1.5))
            # ratio
            with np.errstate(divide="ignore", invalid="ignore"):
                rr = H["sig"][0] / (tot * wdt)
                re_ = np.sqrt(H["sig"][1]) / (tot * wdt)
                band = np.sqrt(totv) / (tot * wdt)
            ok = tot > 0
            ar.fill_between(edges, np.r_[1 - np.hypot(band, LREL), 1 - np.hypot(band[-1], LREL)],
                            np.r_[1 + np.hypot(band, LREL), 1 + np.hypot(band[-1], LREL)], step="post", color="#dddddd", lw=0)
            ar.axhline(1, color="k", ls=":", lw=1)
            ar.errorbar(xc[ok], rr[ok], yerr=re_[ok], fmt="o", ms=2.8, color=PAL["NmCB"], lw=0.8)
            if col == 1:
                with np.errstate(divide="ignore", invalid="ignore"):
                    ar.stairs((fitm / tot)[fr], fe_, color=PAL["fit"], lw=1.4, ls="--", baseline=None)
            ar.set_ylim(0.4, 1.9)
            ar.set_ylabel("data/model", fontsize=12)
            if row == 1:
                ar.set_xlabel(r"$M(\pi^+\pi^-e^+e^-)$ [GeV/$c^2$]")
            else:
                ar.set_xlabel(r"$M(\pi^+\pi^-e^+e^-)$ [GeV/$c^2$]", fontsize=12)
    handles = [Line2D([], [], marker="o", ls="", color=PAL["N"], ms=5, label="data $N$"),
               Line2D([], [], marker="s", ls="", color=PAL["CB"], ms=5, label="CB $=2\\sqrt{N_{++}N_{--}}$"),
               Line2D([], [], marker="o", ls="", color=PAL["NmCB"], ms=5, label="data $N-$CB"),
               Patch(color=PAL["GD"], label="$G_D$ ($\\pi^0$ Dalitz, vector, other)"), Patch(color=PAL["GC"], label="$G_C$ conversion pair"),
               Patch(color=PAL["Geta"], label="$G_\\eta$ other $\\eta$"), Patch(color=PAL["C2"], label="C2 bump $\\pi^0$ Dalitz"),
               Patch(color=PAL["C3"], label="C3 bump conversion"), Patch(color=PAL["C4"], label="C4 $\\eta\\to\\pi^+\\pi^-\\gamma$ conv."),
               Patch(color=PAL["C5"], label="C5 $\\eta\\to\\pi^+\\pi^-e^+e^-$"),
               Line2D([], [], color=PAL["fit"], ls="--", lw=1.8, label="(c) display fit $a_\\eta\\,\\eta+c\\,$non-$\\eta$"),
               Patch(color="#dddddd", label="ratio band: MC $\\oplus$ $L$ 6.9 %")]
    fig.legend(handles=handles, loc="upper center", ncol=6, fontsize=10.5, frameon=False, bbox_to_anchor=(0.5, 0.965))
    fig.suptitle(f"{sel_title(s)};  isBest $\\geq$ 0;  data $\\times$64/63, $L$ = 6079.2 nb$^{{-1}}$;  SMASH $M_4'=0.99615\\,M_4$",
                 fontsize=14, y=0.995)
    save(fig, "e1_m4_" + s)


for s in ["pre_E1", "pre_X2", "S0_E1", "S0_O6"]:
    m4_figure(s)

# ---------------------------------------------------------------- compact overview: linear 0.40-0.70, (b) scaled, selected cuts
def overview(sels, name, title):
    n = len(sels)
    ncol = 3
    nrow = int(math.ceil(n / ncol))
    fig = plt.figure(figsize=(16, 5.6 * nrow + 0.8))
    gs = fig.add_gridspec(3 * nrow, ncol, height_ratios=[3, 1, 0.42] * nrow, hspace=0.08, wspace=0.22, left=0.06, right=0.99,
                          top=1 - 0.9 / (5.6 * nrow + 0.8), bottom=0.02)
    for k, s in enumerate(sels):
        rr_, cc_ = divmod(k, ncol)
        ax = fig.add_subplot(gs[3 * rr_, cc_]); ar = fig.add_subplot(gs[3 * rr_ + 1, cc_], sharex=ax)
        edges, H = spectra(s, "lin")
        wdt = np.diff(edges); xc = 0.5 * (edges[1:] + edges[:-1])
        cb = R[s]["cb"]
        layers = [(H["GD"][0] * cb / wdt, PAL["GD"]), (H["Gconv"][0] * cb / wdt, PAL["GC"]), (H["Geta"][0] * ETA_SCALE / wdt, PAL["Geta"]),
                  (H["C2"][0] * ETA_SCALE / wdt, PAL["C2"]), (H["C3"][0] * ETA_SCALE / wdt, PAL["C3"]), (H["C4"][0] * ETA_SCALE / wdt, PAL["C4"]),
                  (H["C5"][0] * ETA_SCALE / wdt, PAL["C5"])]
        tot = stairs_stack(ax, edges, layers)
        ax.stairs(H["T12"][0] / wdt, edges, color="k", lw=1.2, ls=":")
        ax.errorbar(xc, H["sig"][0] / wdt, yerr=np.sqrt(H["sig"][1]) / wdt, fmt="o", ms=3, color=PAL["NmCB"], lw=0.9)
        ax.set_xlim(0.40, 0.70)
        ax.set_ylim(0, max(np.max(H["sig"][0] / wdt), np.max(H["T12"][0] / wdt), np.max(tot)) * 1.30)
        g = R[s]["reg"]
        ax.text(0.02, 0.97, f"{LAB[s]}" + (": " + CUTDEF[base_cut(s)[1]] if base_cut(s)[1] else "") +
                f"\n$c_b$ = {cb:.3f}; P: (a) {g['P']['ra'][0]:.2f}, (b) {g['P']['rb'][0]:.2f}; B: (a) {g['B']['ra'][0]:.2f}, (b) {g['B']['rb'][0]:.2f}",
                transform=ax.transAxes, va="top", fontsize=10.5)
        plt.setp(ax.get_xticklabels(), visible=False)
        if cc_ == 0:
            ax.set_ylabel(r"d$\sigma$/d$M_4$ [nb/(GeV/$c^2$)]")
            ar.set_ylabel("data/model", fontsize=12)
        okr = (tot > 0.02 * np.max(tot)) & (H["T12"][0] / wdt > 0.02 * np.max(H["T12"][0] / wdt))
        with np.errstate(divide="ignore", invalid="ignore"):
            ar.errorbar(xc[okr], (H["sig"][0] / (tot * wdt))[okr], yerr=(np.sqrt(H["sig"][1]) / (tot * wdt))[okr], fmt="o", ms=2.6,
                        color=PAL["NmCB"], lw=0.8, label="/ (b) scaled")
            ar.errorbar(xc[okr], (H["sig"][0] / H["T12"][0])[okr], yerr=(np.sqrt(H["sig"][1]) / H["T12"][0])[okr], fmt="o", ms=2.6,
                        mfc="white", color="#777777", lw=0.6, label="/ (a) as is")
        ar.axhline(1, color="k", ls=":", lw=1)
        ar.set_ylim(0.5, 1.8)
        ar.set_xlabel(r"$M(\pi^+\pi^-e^+e^-)$ [GeV/$c^2$]", fontsize=12)
        if k == 0:
            ar.legend(fontsize=9, loc="upper left", ncol=2, frameon=False)
    handles = [Line2D([], [], marker="o", ls="", color=PAL["NmCB"], ms=5, label="data $N-$CB"),
               Line2D([], [], color="k", ls=":", lw=1.2, label="SMASH as is (a)"),
               Patch(color=PAL["GD"], label="$G_D\\times c_b$"), Patch(color=PAL["GC"], label="$G_C\\times c_b$"),
               Patch(color=PAL["Geta"], label="$G_\\eta\\times0.4$"), Patch(color=PAL["C2"], label="C2$\\times0.4$"),
               Patch(color=PAL["C3"], label="C3$\\times0.4$"), Patch(color=PAL["C4"], label="C4$\\times0.4$"), Patch(color=PAL["C5"], label="C5$\\times0.4$")]
    fig.legend(handles=handles, loc="upper center", ncol=9, fontsize=10.5, frameon=False, bbox_to_anchor=(0.5, 1 - 0.35 / (5.6 * nrow + 0.8)))
    fig.suptitle(title, fontsize=14, y=0.999)
    save(fig, name)


overview(["pre", "pre_E1", "pre_E2", "pre_E3", "pre_E4", "pre_X1", "pre_X2", "pre_X2old", "pre_O6"], "e1_m4_overview_pre",
         "M$_4$ 0.40-0.70 after the typical cuts (base pre): data $N-$CB against SMASH, (b) $\\eta\\times0.4$, non-$\\eta\\times c_b$ "
         "(0.575-0.70); $\\times$64/63, $L$=6079.2 nb$^{-1}$")
overview(["S0", "S0_E1", "S0_E2", "S0_E3", "S0_E4", "S0_X1", "S0_X2", "S0_X2old", "S0_O6"], "e1_m4_overview_S0",
         "M$_4$ 0.40-0.70 after the typical cuts (base S0): data $N-$CB against SMASH, (b) $\\eta\\times0.4$, non-$\\eta\\times c_b$ "
         "(0.575-0.70); $\\times$64/63, $L$=6079.2 nb$^{-1}$")

# ---------------------------------------------------------------- region ratios and scales per selection
def summary_fig():
    fig, axs = plt.subplots(3, 1, figsize=(15, 13), sharex=True, gridspec_kw=dict(hspace=0.08, left=0.07, right=0.74, top=0.95, bottom=0.08))
    x = np.arange(len(ORDER))
    cols = {"below": "#888888", "B": PAL["C2"], "P": PAL["C4"], "r56": "#1f3fbf", "r70": "#000000"}
    off = {"below": -0.24, "B": -0.12, "P": 0.0, "r56": 0.12, "r70": 0.24}
    for k, key in enumerate(["ra", "rb"]):
        ax = axs[k]
        for w, t in REG:
            v = np.array([R[s]["reg"][w][key][0] for s in ORDER]); e = np.array([R[s]["reg"][w][key][1] for s in ORDER])
            ax.errorbar(x + off[w], v, yerr=e, fmt="o" if k == 0 else "s", ms=5, color=cols[w], label=t, lw=1)
        ax.axhline(1, color="k", ls=":", lw=1)
        ax.axvline(10.5, color="k", lw=0.8)
        ax.set_ylabel("data/SMASH (a) as is" if k == 0 else "data/model (b) scaled")
        ax.set_ylim(0.6, 1.7)
        if k == 0:
            ax.legend(loc="upper left", bbox_to_anchor=(1.005, 1.0), fontsize=11, title="$M_4$ region [GeV/$c^2$]")
    ax = axs[2]
    ax.errorbar(x - 0.1, [R[s]["cb"] for s in ORDER], yerr=[R[s]["cbe"] for s in ORDER], fmt="D", color="#555555", label="$c_b$ of (b) (0.575-0.70)")
    ax.errorbar(x + 0.1, [R[s]["fit"]["c"] for s in ORDER], yerr=[R[s]["fit"]["ec"] for s in ORDER], fmt="s", color="#999999", mfc="white",
                label="$c$ of the display fit (c)")
    ax.errorbar(x, [R[s]["fit"]["a"] for s in ORDER], yerr=[R[s]["fit"]["ea"] for s in ORDER], fmt="o", color=PAL["C3"], label="$a_\\eta$ of the display fit (c)")
    ax.axhline(0.40, color=PAL["C3"], ls="--", lw=1, label="Part III $\\eta$ scale 0.40")
    ax.axhline(1.0, color="k", ls=":", lw=1)
    ax.axvline(10.5, color="k", lw=0.8)
    ax.set_ylabel("scale factor")
    ax.set_ylim(0.0, 1.8)
    ax.legend(loc="upper left", bbox_to_anchor=(1.005, 1.0), fontsize=11)
    ax.set_xticks(x)
    ax.set_xticklabels([LAB[s] for s in ORDER], rotation=45, ha="right")
    fig.suptitle("Data/SMASH per $M_4$ region and the scale factors after each selection (left: base pre, right: base S0)", fontsize=14)
    save(fig, "e1_regions")


# summary_fig()  (unchanged: figs_all/e1_regions.pdf)

# ---------------------------------------------------------------- X17 opening angle study
def x17_oa_fig():
    fig, axs = plt.subplots(1, 3, figsize=(17, 5.6), gridspec_kw=dict(wspace=0.28, left=0.055, right=0.985, top=0.86, bottom=0.13))
    oe = np.arange(42) * 1.0; oe[-1] = 41.0
    ax = axs[0]
    def norm(h):
        return h / np.sum(h)
    ax.stairs(norm(oaC6)[:40], np.arange(41), color=PAL["X17"], lw=2.4, ls="--", label="X17 (C6, 17 MeV/$c^2$), peak band")
    ax.stairs(norm(oaC5x)[:40], np.arange(41), color=PAL["C5"], lw=2.0, label="C5 direct, $M_{ee}$ 15.5-24.5")
    ax.stairs(norm(oaC5)[:40], np.arange(41), color="#007777", lw=1.5, ls=":", label="C5 direct, all $M_{ee}$")
    hd = oaDx / np.sum(oaDx)
    ax.errorbar(np.arange(40) + 0.5, hd[:40], yerr=np.sqrt(oaDxv[:40]) / np.sum(oaDx), fmt="o", ms=3.5, color=PAL["NmCB"],
                label="data $N-$CB, $M_{ee}$ 15.5-24.5")
    for v in (2, 6):
        ax.axvline(v, color="k", ls=":", lw=1)
    ax.set_xlim(0, 40); ax.set_xlabel("OA($e^+e^-$) [deg]"); ax.set_ylabel("fraction per deg")
    ax.set_title("(a) S0, $\\eta$ peak band $M_4$ 0.5325-0.5575", loc="left", fontsize=13)
    ax.legend(fontsize=10)
    ax = axs[1]
    for h, col, ls, lab in [(oaC6, PAL["X17"], "--", "X17 (C6)"), (oaC5x, PAL["C5"], "-", "C5, $M_{ee}$ 15.5-24.5"),
                            (oaC5, "#007777", ":", "C5, all $M_{ee}$"), (oaDx, PAL["NmCB"], "-.", "data $N-$CB, $M_{ee}$ 15.5-24.5")]:
        cum = np.cumsum(h) / np.sum(h)
        ax.plot(np.arange(1, 42), cum, color=col, ls=ls, lw=2, label=lab)
    for v in (2, 6):
        ax.axvline(v, color="k", ls=":", lw=1)
    ax.set_xlim(0, 40); ax.set_ylim(0, 1.02)
    ax.set_xlabel("OA cut value [deg]"); ax.set_ylabel("fraction with OA < cut")
    f6 = OAF["C6 peak band"]
    ax.set_title(f"(b) X17: OA<6$^\\circ$ keeps {f6['lt6']:.2f}, 2-6$^\\circ$ keeps {f6['b26']:.2f}", loc="left", fontsize=13)
    ax.legend(fontsize=10, loc="lower right")
    ax = axs[2]
    ax.stairs(peeC6[:100] / np.sum(peeC6), 20. * np.arange(101), color=PAL["X17"], lw=2.2, ls="--", label="X17 (C6)")
    ax.stairs(peeC5[:100] / np.sum(peeC5), 20. * np.arange(101), color=PAL["C5"], lw=1.8, label="C5 direct (all $M_{ee}$)")
    ax.axvline(p_thr, color="k", ls=":", lw=1.2)
    ax.text(p_thr + 330, ax.get_ylim()[1] * 0.62, f"dotted: {p_thr:.0f} MeV/$c$, where a\nsymmetric 17 MeV/$c^2$ pair\nhas OA = 6$^\\circ$", fontsize=10.5, va="top")
    ax.set_xlim(0, 1400); ax.set_xlabel("$p(e^+e^-)$ [MeV/$c$]"); ax.set_ylabel("fraction per 20 MeV/$c$")
    ax.set_title(f"(c) pair momentum; X17 median {peeF['median']:.0f} MeV/$c$", loc="left", fontsize=13)
    ax.legend(fontsize=10)
    fig.suptitle("Opening angle of X17 $\\to e^+e^-$ in the HADES acceptance (SMASH, S0, isBest $\\geq$ 0): an X17 of 17 MeV/$c^2$ needs "
                 f"$p_{{ee}}>{p_thr:.0f}$ MeV/$c$ for OA < 6$^\\circ$", fontsize=13.5)
    save(fig, "e1_x17_oa")


x17_oa_fig()

# ---------------------------------------------------------------- M(ee) in the peak band with the X17 at the S0 limit
MEEE = 2. * np.arange(32)    # 0..62 MeV, the last cell (overflow) dropped


def mee_panel(ax, ar, s, show_leg=False):
    cb = R[s]["cb"]
    N, NV = D(s, "mee", "N"); CB, CBv = D(s, "mee", "CB")
    comps = {c: S(s, c, "mee")[0] * NB2RAW for c in ["GD", "Gconv", "Geta", "C2", "C3", "C4", "C5", "C6"]}
    n = 31
    w = 2.0
    layers = [(CB[:n] / w, PAL["CB"]), (comps["GD"][:n] * cb / w, PAL["GD"]), (comps["Gconv"][:n] * cb / w, PAL["GC"]),
              ((comps["Geta"] + comps["C2"] + comps["C3"])[:n] * ETA_SCALE / w, PAL["Geta"]),
              (comps["C4"][:n] * ETA_SCALE / w, PAL["C4"]), (comps["C5"][:n] * ETA_SCALE / w, PAL["C5"])]
    tot = stairs_stack(ax, MEEE, layers)
    xline = NX_UL * S(s, "C6", "mee")[0][:n] / C6X2D_S0 / w
    okm = tot > 0
    ax.stairs(np.where(okm, tot + xline, np.nan), MEEE, color=PAL["X17"], lw=2.2, ls="--", baseline=None)
    xc = MEEE[:-1] + 1
    ax.errorbar(xc, N[:n] / w, yerr=np.sqrt(NV[:n]) / w, fmt="o", ms=3.2, color=PAL["N"], lw=0.8)
    ax.set_yscale("log")
    ax.set_xlim(0, 62)
    ax.set_ylim(3., max(np.max(N[:n] / w), np.max(tot)) * 6)
    x = R[s]["x17"]
    ax.text(0.97, 0.96, f"{LAB[s]}\n$c_b$={cb:.3f}; X17 window: $N$={x['N']:.0f}\n$N_X$(limit)={x['NX']:.0f}, Z={x['Z']:.2f}, "
            f"$\\varepsilon_X/\\varepsilon_X$(S0)={x['eps_rel']:.2f}", transform=ax.transAxes, ha="right", va="top", fontsize=10,
            bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", pad=1.5))
    plt.setp(ax.get_xticklabels(), visible=False)
    with np.errstate(divide="ignore", invalid="ignore"):
        ar.errorbar(xc, N[:n] / w / tot, yerr=np.sqrt(NV[:n]) / w / tot, fmt="o", ms=2.8, color=PAL["N"], lw=0.8)
        ar.stairs(np.where(okm, (tot + xline) / tot, np.nan), MEEE, color=PAL["X17"], lw=2.0, ls="--", baseline=None)
    ar.axhline(1, color="k", ls=":", lw=1)
    ar.axvspan(15.5, 24.5, color="#fff2a8", alpha=0.6, lw=0)
    ar.set_ylim(0.6, 1.4)
    ar.set_xlabel("$M(e^+e^-)$ [MeV/$c^2$]", fontsize=12)


def mee_fig(sels, name, title):
    ncol = 4
    nrow = int(math.ceil(len(sels) / ncol))
    fig = plt.figure(figsize=(18, 5.4 * nrow + 0.9))
    gs = fig.add_gridspec(3 * nrow, ncol, height_ratios=[3, 1.1, 0.45] * nrow, hspace=0.08, wspace=0.2, left=0.05, right=0.99,
                          top=1 - 1.0 / (5.4 * nrow + 0.9), bottom=0.02)
    for k, s in enumerate(sels):
        rr_, cc_ = divmod(k, ncol)
        ax = fig.add_subplot(gs[3 * rr_, cc_]); ar = fig.add_subplot(gs[3 * rr_ + 1, cc_], sharex=ax)
        mee_panel(ax, ar, s)
        if cc_ == 0:
            ax.set_ylabel("rows per MeV/$c^2$")
            ar.set_ylabel("data/model", fontsize=12)
    handles = [Line2D([], [], marker="o", ls="", color=PAL["N"], ms=5, label="data $N$"), Patch(color=PAL["CB"], label="CB (like-sign data)"),
               Patch(color=PAL["GD"], label="$G_D\\times c_b$"), Patch(color=PAL["GC"], label="$G_C\\times c_b$ (conv. pair)"),
               Patch(color=PAL["Geta"], label="other $\\eta$: C2+C3+$G_\\eta$ ($\\times$0.4)"), Patch(color=PAL["C4"], label="C4 conv. ($\\times$0.4)"),
               Patch(color=PAL["C5"], label="C5 direct ($\\times$0.4)"),
               Line2D([], [], color=PAL["X17"], ls="--", lw=2.2, label="+ X17 at the S0 limit (392 rows in X2D)"),
               Patch(color="#fff2a8", label="X17 window 15.5-24.5")]
    fig.legend(handles=handles, loc="upper center", ncol=5, fontsize=10.5, frameon=False, bbox_to_anchor=(0.5, 1 - 0.3 / (5.4 * nrow + 0.9)))
    fig.suptitle(title, fontsize=14, y=0.999)
    save(fig, name)


mee_fig(ORDER_PRE, "e1_mee_pre",
        "$M(e^+e^-)$ in the $\\eta$ peak band ($M_4$ 0.5325-0.5575), base pre: data $N$ against CB + SMASH (b); dashed: + X17 at the S0 limit, "
        "scaled by the simulated X17 efficiency (raw rows)")
mee_fig(ORDER_S0, "e1_mee_S0",
        "$M(e^+e^-)$ in the $\\eta$ peak band ($M_4$ 0.5325-0.5575), base S0: data $N$ against CB + SMASH (b); dashed: + X17 at the S0 limit, "
        "scaled by the simulated X17 efficiency (raw rows)")

# ---------------------------------------------------------------- OA distributions in the peak band
def oa_fig(sels, name):
    fig = plt.figure(figsize=(17, 10.5))
    gs = fig.add_gridspec(2, len(sels), hspace=0.32, wspace=0.22, left=0.055, right=0.99, top=0.86, bottom=0.07)
    for j, s in enumerate(sels):
        cb = R[s]["cb"]
        for row, tag in enumerate(["oa", "oaX"]):
            ax = fig.add_subplot(gs[row, j])
            d, dv = D(s, tag)
            comps = {c: S(s, c, tag)[0] * NB2RAW for c in ["GD", "Gconv", "Geta", "C2", "C3", "C4", "C5", "C6"]}
            e = np.arange(41.)
            layers = [(comps["GD"][:40] * cb, PAL["GD"]), (comps["Gconv"][:40] * cb, PAL["GC"]),
                      ((comps["Geta"] + comps["C2"] + comps["C3"])[:40] * ETA_SCALE, PAL["Geta"]),
                      (comps["C4"][:40] * ETA_SCALE, PAL["C4"]), (comps["C5"][:40] * ETA_SCALE, PAL["C5"])]
            tot = stairs_stack(ax, e, layers)
            xl = NX_UL * S(s, "C6", tag)[0][:40] / C6X2D_S0
            ax.stairs(xl, e, color=PAL["X17"], lw=2.4, ls="--")
            ok = d[:40] > 0
            xc = np.arange(40) + 0.5
            ax.errorbar(xc[ok], d[:40][ok], yerr=np.sqrt(dv[:40][ok]), fmt="o", ms=3.4, color=PAL["NmCB"], lw=0.9)
            CBh = D(s, tag, "CB")[0]
            ax.stairs(CBh[:40], e, color=PAL["CB"], lw=1.4)
            ax.set_yscale("log")
            ax.set_xlim(0, 40)
            ax.set_ylim(max(0.5, np.min(xl[xl > 0]) * 0.3 if np.any(xl > 0) else 1), max(np.max(d[:40]), np.max(tot)) * 3)
            for v in (2, 6):
                ax.axvline(v, color="k", ls=":", lw=1)
            ax.set_xlabel("OA($e^+e^-$) [deg]")
            if j == 0:
                ax.set_ylabel("rows per deg")
            x = R[s]["x17"]
            win = "all $M_{ee}$" if tag == "oa" else "$M_{ee}$ 15.5-24.5 MeV/$c^2$ (X17 window)"
            ax.set_title(f"{LAB[s]}, peak band, {win}", fontsize=12, loc="left")
            xs = NX_UL * S(s, "C6", tag)[0] / C6X2D_S0
            ax.text(0.97, 0.96, f"X17 at the S0 limit: {np.sum(xs):.0f} rows\nof which OA<6$^\\circ$: {np.sum(xs[:6]):.0f}",
                    transform=ax.transAxes, ha="right", va="top", fontsize=10, bbox=dict(facecolor="white", alpha=0.85, edgecolor="none"))
    handles = [Line2D([], [], marker="o", ls="", color=PAL["NmCB"], ms=5, label="data $N-$CB"), Line2D([], [], color=PAL["CB"], lw=1.4, label="CB"),
               Patch(color=PAL["GD"], label="$G_D\\times c_b$"), Patch(color=PAL["GC"], label="$G_C\\times c_b$"),
               Patch(color=PAL["Geta"], label="other $\\eta$ $\\times$0.4"), Patch(color=PAL["C4"], label="C4 $\\times$0.4"),
               Patch(color=PAL["C5"], label="C5 $\\times$0.4"), Line2D([], [], color=PAL["X17"], ls="--", lw=2.4, label="X17 alone, at the S0 limit")]
    fig.legend(handles=handles, loc="upper center", ncol=8, fontsize=10.5, frameon=False, bbox_to_anchor=(0.5, 0.955))
    fig.suptitle("OA($e^+e^-$) in the $\\eta$ peak band ($M_4$ 0.5325-0.5575): data $N-$CB, SMASH (b) stacked, X17 at the present S0 limit",
                 fontsize=14, y=0.995)
    save(fig, name)


# oa_fig(["S0", "pre_E2", "S0_E2"], "e1_oa")  (unchanged: figs_all/e1_oa.pdf)

# ---------------------------------------------------------------- X17 sensitivity per selection
assert len(ORDER) == 22 and ORDER[10] == "pre_O6" and ORDER[11] == "S0", ORDER


def sens_fig():
    fig, axs = plt.subplots(2, 1, figsize=(16, 9), sharex=True, gridspec_kw=dict(hspace=0.08, left=0.07, right=0.67, top=0.93, bottom=0.12))
    x = np.arange(len(ORDER))
    ax = axs[0]
    def _seg(ax, y, fmt, **kw):
        y = np.asarray(y, dtype=float)
        ax.plot(x[:11], y[:11], fmt, **kw)
        kw.pop("label", None)
        ax.plot(x[11:], y[11:], fmt, **kw)
    _seg(ax, [R[s]["x17"]["eps_rel"] for s in ORDER], "o-", color=PAL["X17"], label="X17 efficiency / S0 (SMASH, C6)")
    _seg(ax, [R[s]["x17"]["B_rel"] for s in ORDER], "s-", color="k", label="data rows $N$ in the X17 window / S0")
    _seg(ax, [R[s]["x17"]["C5b"] / R["S0"]["x17"]["C5b"] for s in ORDER], "^-", color=PAL["C5"], label="direct C5 (SMASH) / S0")
    ax.axhline(1, color="k", ls=":", lw=1); ax.axvline(10.5, color="k", lw=0.8)
    ax.set_ylabel("relative to S0"); ax.set_ylim(0, 1.6)
    ax.legend(loc="upper left", bbox_to_anchor=(1.005, 1), fontsize=11)
    ax = axs[1]
    _seg(ax, [R[s]["x17"]["proxy"] for s in ORDER], "D-", color="#1f3fbf", label="counting proxy $\\sqrt{B}/(\\varepsilon/\\varepsilon_{S0})/\\sqrt{B_{S0}}$\n(smaller = better)")
    _seg(ax, [R[s]["x17"]["Z"] / R["S0"]["x17"]["Z"] for s in ORDER], "o--", color=PAL["X17"], mfc="white",
            label="$Z$ of an X17 at the S0 limit / $Z$(S0)\n(larger = better)")
    _seg(ax, [R[s]["x17"]["SB_gain"] for s in ORDER], "s:", color="#555555", mfc="white",
            label="S/B gain $=(\\varepsilon/\\varepsilon_{S0})/(B/B_{S0})$\n(larger = better if the background\nsystematics dominate)")
    ax.axhline(1, color="k", ls=":", lw=1); ax.axvline(10.5, color="k", lw=0.8)
    ax.set_ylabel("relative to S0"); ax.set_ylim(0.5, 1.5)
    ax.legend(loc="upper left", bbox_to_anchor=(1.005, 1), fontsize=11)
    ax.set_xticks(x); ax.set_xticklabels([LAB[s] for s in ORDER], rotation=45, ha="right")
    fig.suptitle("X17 window ($M_4$ 0.5325-0.5575, $M_{ee}$ 15.5-24.5 MeV/$c^2$): what each cut keeps (simple counting proxy, not the frozen fit)",
                 fontsize=13.5)
    save(fig, "e1_x17_sens")


sens_fig()

# ---------------------------------------------------------------- dump of the key numbers (for the note)
with open(TAB + "e1_numbers.txt", "w") as f:
    f.write("# e1_numbers.txt -- plots/all_report/e1_ana.py: key numbers (machine readable: sel key value)\n")
    for s in ORDER:
        r = R[s]
        for w, _ in REG + [("R2", ""), ("post", ""), ("Fit", ""), ("all", ""), ("CR", ""), ("PB", "")]:
            g = r["reg"][w]
            f.write(f"{s} {w}_N {g['N']:.0f}\n{s} {w}_CB {g['CB']:.2f} {g['CBe']:.2f}\n{s} {w}_dnb {g['dnb']:.5f} {g['dnbE']:.5f}\n"
                    f"{s} {w}_sim {g['sim']:.5f} {g['simE']:.5f}\n{s} {w}_eta {g['eta']:.5f}\n{s} {w}_bkg {g['bkg']:.5f}\n"
                    f"{s} {w}_ra {g['ra'][0]:.4f} {g['ra'][1]:.4f}\n{s} {w}_rb {g['rb'][0]:.4f} {g['rb'][1]:.4f}\n")
        f.write(f"{s} delta a={r['delta']['a']:.4f} b={r['delta']['b']:.4f} c={r['delta']['c']:.4f}\n")
        f.write(f"{s} cb {r['cb']:.4f} {r['cbe']:.4f}\n{s} chi2a {r['chi2a'][0]:.2f} {r['chi2a'][1]}\n{s} chi2b {r['chi2b'][0]:.2f} {r['chi2b'][1]}\n")
        f.write(f"{s} k1 {r['k1']:.4f} chi2 {r['chi2k1'][0]:.2f} {r['chi2k1'][1]}\n")
        ft = r["fit"]
        f.write(f"{s} fit a {ft['a']:.4f} {ft['ea']:.4f} c {ft['c']:.4f} {ft['ec']:.4f} rho {ft['rho']:.3f} chi2 {ft['chi2']:.2f} {ft['ndf']}\n")
        x = r["x17"]
        f.write(f"{s} x17 " + " ".join(f"{k}={v:.5g}" for k, v in x.items()) + "\n")
    for k, v in OAF.items():
        f.write(f"OA {k.replace(' ', '_')} " + " ".join(f"{a}={b:.4f}" for a, b in v.items()) + "\n")
    f.write(f"PEE C6 median={peeF['median']:.1f} share_above_{p_thr:.0f}={peeF['gt325']:.4f}\n")
    f.write(f"C6_X2D_S0 {C6X2D_S0:.3f}\nrepro_rows_differing {NBAD_REPRO}\n")

print("\n".join(chk))
print("repro rows differing:", NBAD_REPRO)
print(open(TAB + "e1_table.txt").read())
