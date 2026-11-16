#!/usr/bin/env python3
# e1v_ana.py -- independent verifier of package e1 (prefix e1v): numbers from the raw sums of e1v_verify.C
# Own implementation of the Part III conventions (no code of e1_ana.py / e1_fill.C / wp0_fill.C reused):
#   data nb = (N - CB) x 64/63 / 6079.2; CB = 2 sqrt(N++ N--) in 10 MeV cells of the 0.30-1.10 grid, linear interpolation
#   to the 2.5 MeV centres (constant beyond the first / last centre), windows with fractional fine bins; CB of the X17
#   window = sum over the peak-band fine bins of CB x f(M(ee) window | 50 MeV band) (like-sign share); k = 1
#   SMASH nb = sum w x RowFactor(108) / 221.05 (own mode A), M4' = 0.99615 M4; X17 C6 per unit BR (weight 2)
# Run: python3 -I plots/all_report/e1v_ana.py
import json
import math

import numpy as np

D = json.load(open("/hdd3/przygoda/HADES/FAT_out/all_report/e1v/e1v_data.json"))
S = json.load(open("/hdd3/przygoda/HADES/FAT_out/all_report/e1v/e1v_sim.json"))
NF, FLO, FW = 320, 0.30, 0.0025
OW = D["oWin"]
WIN = D["wins"]
iw = {w: i for i, w in enumerate(WIN)}
WLO = dict(B=0.4550, P=0.53364, r56=0.56, CR=0.575, PB=0.5325, X17=0.5325, X2D=0.5225, X2Dall=0.5225, Fit=0.40, R5=0.565)
WHI = dict(B=0.5225, P=0.55416, r56=0.70, CR=0.70, PB=0.5575, X17=0.5575, X2D=0.6100, X2Dall=0.6100, Fit=0.70, R5=0.70)
L, PT3, LOLD, LSIM = 6079.2, 64. / 63., 6080.5, 221.05
NB2RAW = L / PT3
ETA = 0.40
NX_UL = 392.

out = []


def P(s=""):
    print(s)
    out.append(s)


# ------------------------------------------------------------------ CB
def cb_model(pp, mm, lo_cell=0, hi_cell=80):
    """fine CB (320) and the interpolation matrix A[i, j] (fine i, wide cell j); cells lo_cell..hi_cell-1 only
    (old convention: 0.35-1.00 = cells 5..69, constant extrapolation beyond)"""
    wp = pp[:NF].reshape(80, 4).sum(1)
    wm = mm[:NF].reshape(80, 4).sum(1)
    cbW = 2. * np.sqrt(wp * wm)
    varW = wp + wm
    A = np.zeros((NF, 80))
    nw = hi_cell - lo_cell
    for i in range(NF):
        x = FLO + (i + 0.5) * FW
        t = (x - (FLO + lo_cell * 0.01 + 0.005)) / 0.01
        if t <= 0:
            k0 = k1 = 0; f = 0.
        elif t >= nw - 1:
            k0 = k1 = nw - 1; f = 0.
        else:
            k0 = int(t); k1 = k0 + 1; f = t - k0
        A[i, lo_cell + k0] += (1 - f) / 4.
        A[i, lo_cell + k1] += f / 4.
    return A @ cbW, A, varW


def frac_w(lo, hi):
    u = np.zeros(NF)
    for i in range(NF):
        a, b = FLO + i * FW, FLO + (i + 1) * FW
        o = min(b, hi) - max(a, lo)
        u[i] = max(0., min(1., o / FW))
    return u


def cb_sum(u, cb, A, varW):
    c = u @ A
    return float(u @ cb), float(np.sum(c * c * varW))


def arr(sel, key):
    return np.array(D["sel"][sel][key], dtype=float)


def data_sel(sel, old=False):
    us, usV, pp, mm = arr(sel, "us"), arr(sel, "usV"), arr(sel, "pp"), arr(sel, "mm")
    cb, A, varW = cb_model(pp, mm, 5, 70) if old else cb_model(pp, mm)
    r = {"us": us, "usV": usV, "cb": cb, "A": A, "varW": varW, "npass": D["sel"][sel]["npass"]}
    r["win"] = {}
    for w in ["B", "P", "r56", "CR", "Fit", "R5", "PB"]:
        N, NV = us[OW + iw[w]], usV[OW + iw[w]]
        C, CV = cb_sum(frac_w(WLO[w], WHI[w]), cb, A, varW)
        r["win"][w] = dict(N=N, NV=NV, CB=C, CBV=CV, sig=N - C, sigV=NV + CV)
    # X17 window: peak-band fine bins 93..102 x like-sign share of the M(ee) window in the band
    lsx = arr(sel, "lsx").reshape(16, 2)
    f_in = np.where(lsx.sum(1) > 0, lsx[:, 0] / np.maximum(lsx.sum(1), 1e-300), 0.)
    u = np.zeros(NF)
    for i in range(93, 103):
        u[i] = f_in[min(15, i // 20)]
    C, CV = cb_sum(u, cb, A, varW)
    N, NV = us[OW + iw["X17"]], usV[OW + iw["X17"]]
    r["win"]["X17"] = dict(N=N, NV=NV, CB=C, CBV=CV, sig=N - C, sigV=NV + CV)
    return r


def sim_sel(sel, comp):
    c = S["sel"][sel][comp]
    return np.array(c["s"]) / LSIM, np.array(c["v"]) / LSIM ** 2


def swin(sel, comp, w):
    s, v = sim_sel(sel, comp)
    return s[OW + iw[w]], v[OW + iw[w]]


SELS = ["pre", "E1", "E2", "X1", "S0", "S0+E1", "X3", "S0+X1"]
DS = {s: data_sel(s) for s in SELS}

# ------------------------------------------------------------------ 1. data
P("# e1v_ana.py -- independent recomputation (prefix e1v); data raw rows, CB k = 1, linear CB errors")
P("## 1. data: rows passing (us, pp, mm) and N, CB, N - CB in the windows (raw rows)")
P(f"{'sel':7s} {'us':>9s} {'pp':>7s} {'mm':>7s} | {'win':4s} {'N':>9s} {'CB':>16s} {'N-CB':>10s} {'nb (x64/63/L)':>14s}")
for s in SELS:
    d = DS[s]
    for w in ["B", "P", "r56", "CR", "X17"]:
        x = d["win"][w]
        nb = x["sig"] * PT3 / L
        P(f"{s:7s} {d['npass'][0]:9d} {d['npass'][1]:7d} {d['npass'][2]:7d} | {w:4s} {x['N']:9.0f} {x['CB']:9.1f} +- {math.sqrt(x['CBV']):5.1f} "
          f"{x['sig']:10.1f} {nb:14.4f}")

# ------------------------------------------------------------------ 2. SMASH as is, ratios (a)
P("\n## 2. SMASH as is (nb, T12 = eta + non-eta) and data/SMASH (a)")
reg = ["B", "P", "r56", "CR", "R5"]
RES = {}
for s in SELS:
    d = DS[s]
    RES[s] = {}
    for w in reg:
        e, ev = swin(s, "eta", w)
        b, bv = swin(s, "non", w)
        t, tv = e + b, ev + bv
        x = d["win"][w]
        dn, dnv = x["sig"] * PT3 / L, x["sigV"] * (PT3 / L) ** 2
        ra = dn / t
        rae = ra * math.sqrt(dnv / dn ** 2 + tv / t ** 2)
        RES[s][w] = dict(dn=dn, dnv=dnv, eta=e, etav=ev, non=b, nonv=bv, T=t, Tv=tv, ra=ra, rae=rae)
    P(f"{s:7s} " + "  ".join(f"{w}: data {RES[s][w]['dn']:8.4f} T12 {RES[s][w]['T']:8.4f} (a) {RES[s][w]['ra']:.3f}+-{RES[s][w]['rae']:.3f}" for w in ["B", "P", "r56"]))

# ------------------------------------------------------------------ 3. scaled (b): c_b from 0.575-0.70
P("\n## 3. scaled SMASH (b): eta x 0.40, non-eta x c_b, c_b = (data - 0.40 eta) / non-eta in 0.575-0.70")
for s in SELS:
    r = RES[s]["CR"]
    cbf = (r["dn"] - ETA * r["eta"]) / r["non"]
    cbe = math.sqrt(r["dnv"] + ETA ** 2 * r["etav"] + cbf ** 2 * r["nonv"]) / r["non"]
    RES[s]["cb"] = (cbf, cbe)
    rb = {w: RES[s][w]["dn"] / (ETA * RES[s][w]["eta"] + cbf * RES[s][w]["non"]) for w in ["B", "P", "r56"]}
    RES[s]["rb"] = rb
    P(f"{s:7s} c_b = {cbf:.4f} +- {cbe:.4f}   (b): B {rb['B']:.3f}  P {rb['P']:.3f}  0.56-0.70 {rb['r56']:.3f}")


# ------------------------------------------------------------------ 4. display fit (c) and goodness, 5 MeV bins 0.40-0.70
def b5(a):
    return a[:NF].reshape(160, 2).sum(1)


def data_b5(s):
    d = DS[s]
    N = b5(d["us"]); NV = b5(d["usV"])
    cb5 = b5(d["cb"])
    A5 = d["A"].reshape(160, 2, 80).sum(1)
    CV = np.array([np.sum(A5[k] ** 2 * d["varW"]) for k in range(160)])
    return (N - cb5) * PT3 / L, (NV + CV) * (PT3 / L) ** 2


def sim_b5(s, comp):
    v, vv = sim_sel(s, comp)
    return b5(v), b5(vv)


X5 = 0.30 + 0.005 * (np.arange(160) + 0.5)
M = (X5 > 0.40) & (X5 < 0.70)


def fit2(dv, dvar, ev, evar, bv, bvar):
    a, c = 0.4, 1.3
    for _ in range(100):
        var = dvar + a * a * evar + c * c * bvar
        w = 1. / var[M]
        X = np.vstack([ev[M], bv[M]]).T
        Fm = X.T @ (X * w[:, None])
        y = X.T @ (dv[M] * w)
        cov = np.linalg.inv(Fm)
        a, c = cov @ y
    var = dvar + a * a * evar + c * c * bvar
    chi2 = float(np.sum((dv[M] - a * ev[M] - c * bv[M]) ** 2 / var[M]))
    return a, math.sqrt(cov[0, 0]), c, math.sqrt(cov[1, 1]), cov[0, 1] / math.sqrt(cov[0, 0] * cov[1, 1]), chi2, int(M.sum()) - 2


def delta(dv, dvar, mv, mvar, npar):
    r2 = (dv[M] - mv[M]) ** 2; var = dvar[M] + mvar[M]; m2 = mv[M] ** 2; ndf = int(M.sum()) - npar
    f = lambda x: float(np.sum(r2 / (var + x * x * m2))) - ndf
    if f(0.) <= 0:
        return 0.
    lo, hi = 0., 2.
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if f(mid) > 0 else (lo, mid)
    return 0.5 * (lo + hi)


P("\n## 4. 0.40-0.70, 5 MeV bins: chi2/ndf and delta of (a) and (b); display fit (c) a_eta, c")
for s in SELS:
    dv, dvar = data_b5(s)
    ev, evar = sim_b5(s, "eta")
    bv, bvar = sim_b5(s, "non")
    cbf = RES[s]["cb"][0]
    ta, tav = ev + bv, evar + bvar
    tb, tbv = ETA * ev + cbf * bv, ETA ** 2 * evar + cbf ** 2 * bvar
    chia = float(np.sum((dv[M] - ta[M]) ** 2 / (dvar[M] + tav[M])))
    chib = float(np.sum((dv[M] - tb[M]) ** 2 / (dvar[M] + tbv[M])))
    a, ea, c, ec, rho, chic, ndfc = fit2(dv, dvar, ev, evar, bv, bvar)
    RES[s]["fit"] = (a, ea, c, ec, rho, chic, ndfc)
    da, db = delta(dv, dvar, ta, tav, 0), delta(dv, dvar, tb, tbv, 1)
    RES[s]["delta"] = (da, db)
    rows = float(np.sum(DS[s]["us"][40:160]))
    P(f"{s:7s} chi2/ndf (a) {chia / 60:6.1f}  delta(a) {da:.3f} | (b) {chib / 59:5.2f} delta(b) {db:.3f} | fit (c): a_eta {a:.3f} +- {ea:.3f}, "
      f"c {c:.3f} +- {ec:.3f}, rho {rho:+.2f}, chi2/ndf {chic / ndfc:.2f} | rows 0.40-0.70 {rows:.0f}")

# ------------------------------------------------------------------ 5. X17
P("\n## 5. X17 window (M4 0.5325-0.5575, M(ee) 15.5-24.5 MeV/c2): C6 efficiency relative to S0 and counting proxy")
c6x2d, _ = swin("S0", "C6", "X2D")
c6x2da, _ = swin("S0", "C6", "X2Dall")
P(f"C6(S0) per unit BR: X2D M4 range, M(ee) 2-62 MeV/c2: {c6x2d:.2f} nb; X2D M4 range, all M(ee): {c6x2da:.2f} nb (x17d: 6257.5)")
c6S0, _ = swin("S0", "C6", "X17")
c5S0, _ = swin("S0", "C5", "X17")
N0 = DS["S0"]["win"]["X17"]["N"]
P(f"{'sel':7s} {'C6 [nb/BR]':>11s} {'eps/eps_S0':>10s} {'C5 rel':>7s} {'N':>6s} {'CB':>12s} {'N-CB':>7s} {'N_X':>6s} {'Z':>5s} {'B/B_S0':>6s} {'proxy':>6s} {'S/B gain':>8s}")
for s in SELS:
    c6, _ = swin(s, "C6", "X17")
    c5, _ = swin(s, "C5", "X17")
    x = DS[s]["win"]["X17"]
    eps = c6 / c6S0
    NX = NX_UL * c6 / c6x2da
    Z = NX / math.sqrt(x["N"] + x["CBV"])
    Brel = x["N"] / N0
    proxy = math.sqrt(x["N"]) / eps / math.sqrt(N0)
    RES[s]["x17"] = dict(eps=eps, NX=NX, Z=Z, proxy=proxy, N=x["N"], CB=x["CB"])
    P(f"{s:7s} {c6:11.2f} {eps:10.3f} {c5 / c5S0:7.3f} {x['N']:6.0f} {x['CB']:6.1f} +- {math.sqrt(x['CBV']):3.1f} {x['sig']:7.0f} {NX:6.1f} {Z:5.2f} "
          f"{Brel:6.3f} {proxy:6.3f} {eps / Brel:8.3f}")
# raw-M(pipi) variant of the SMASH cut
for s, sr in [("E1", "E1raw"), ("X1", "X1raw")]:
    c6a, _ = swin(s, "C6", "X17"); c6b, _ = swin(sr, "C6", "X17")
    ta = swin(s, "eta", "r56")[0] + swin(s, "non", "r56")[0]; tb = swin(sr, "eta", "r56")[0] + swin(sr, "non", "r56")[0]
    pa = swin(s, "eta", "P")[0] + swin(s, "non", "P")[0]; pb = swin(sr, "eta", "P")[0] + swin(sr, "non", "P")[0]
    P(f"SMASH M(pipi) scale in the cut, {s}: C6 X17 window scaled/raw {c6a / c6b:.4f}; T12 0.56-0.70 scaled/raw {ta / tb:.4f}; T12 P scaled/raw {pa / pb:.4f}")

# ------------------------------------------------------------------ 6. C6 opening angle and pair momentum (S0)
P("\n## 6. simulated X17 (C6, S0): opening angle and pair momentum")
c6 = np.array(S["sel"]["S0"]["C6"]["s"])
oPB, oX = c6[D["oOaPB"]:D["oOaPB"] + 41], c6[D["oOaX"]:D["oOaX"] + 41]
for nm, h in [("peak band (all M(ee))", oPB), ("X17 window", oX)]:
    t = h.sum()
    P(f"{nm:22s}: OA<2 {h[:2].sum() / t:.3f}  2-6 {h[2:6].sum() / t:.3f}  OA<6 {h[:6].sum() / t:.3f}  OA>6 {h[6:].sum() / t:.3f}")
c5 = np.array(S["sel"]["S0"]["C5"]["s"])
h = c5[D["oOaX"]:D["oOaX"] + 41]
P(f"C5 X17 window          : OA<6 {h[:6].sum() / h.sum():.3f}")
d0 = DS["S0"]
# data N - CB in the X17 window per OA bin is not recomputed (only N): N with OA<6 share of the e+e- rows
hN = d0["us"][D["oOaX"]:D["oOaX"] + 41]
P(f"data N (not CB-subtracted), X17 window: OA<6 share {hN[:6].sum() / hN.sum():.3f}")
pPB = np.array(S["pPB"]); pX = np.array(S["pX"])
P(f"C6 pair momentum, peak band: median {np.median(pPB):.1f} MeV/c, share above 325 MeV/c {np.mean(pPB > 325.):.3f} (rows {len(pPB)})")
P(f"C6 pair momentum, X17 window: median {np.median(pX):.1f} MeV/c, share above 325 MeV/c {np.mean(pX > 325.):.3f}")
P(f"symmetric 17 MeV/c2 pair with OA = 6 deg: P = m / tan(3 deg) = {17. / math.tan(math.radians(3.)):.1f} MeV/c")

# ------------------------------------------------------------------ 7. Part III anchors and old convention
P("\n## 7. anchors")
x = DS["S0"]["win"]["P"]
P(f"S0 P: N {x['N']:.0f}, CB {x['CB']:.1f} +- {math.sqrt(x['CBV']):.1f}, (N-CB) x 64/63 / L = {x['sig'] * PT3 / L:.4f} nb, T12 {RES['S0']['P']['T']:.4f} nb")
x = DS["pre"]["win"]["P"]
P(f"pre P: N {x['N']:.0f}, CB {x['CB']:.1f} +- {math.sqrt(x['CBV']):.1f}, T12 {RES['pre']['P']['T']:.4f} nb")
P(f"pre R5 T12 {RES['pre']['R5']['T']:.4f} nb (wp1: 108.4174)")
# old convention (L 6080.5, no 64/63, CB cells 0.35-1.00, SMASH raw M4 and raw M(pipi))
for s, so, sim in [("pre", "pre", "pre_old"), ("E1", "E1", "E1_old")]:
    d = data_sel(s, old=True)
    line = []
    for w, lo, hi in [("bump", 0.455, 0.5225), ("peak", 0.53364, 0.55416), ("r56", 0.56, 0.70)]:
        key = {"bump": "B", "peak": "P", "r56": "r56"}[w]
        N = d["us"][OW + iw[key]]
        C, _ = cb_sum(frac_w(lo, hi), d["cb"], d["A"], d["varW"])
        dn = (N - C) / LOLD
        t = swin(sim, "eta", key)[0] + swin(sim, "non", key)[0]
        line.append(f"{w}: data {dn:.4f} nb, sim {t:.3f} nb, ratio {dn / t:.3f}")
    P(f"old convention {s:4s}: " + "; ".join(line))

open("/hdd3/przygoda/HADES/FAT_out/all_report/e1v/e1v_numbers.txt", "w").write("\n".join(out) + "\n")

# ------------------------------------------------------------------ 8. eta share and the scaled model in the X17 window
P("\n## 8. eta share of SMASH in B and P (as is | scaled b), C5 share of model (b) in P; data N / (CB + model b) in the X17 window")
for s in SELS:
    cbf = RES[s]["cb"][0]
    sh = []
    for w in ["B", "P"]:
        e, b = RES[s][w]["eta"], RES[s][w]["non"]
        sh.append(f"{w}: {e / (e + b):.3f} | {ETA * e / (ETA * e + cbf * b):.3f}")
    c5p = swin(s, "C5", "P")[0]
    c5sh = ETA * c5p / (ETA * RES[s]["P"]["eta"] + cbf * RES[s]["P"]["non"])
    ex, bx = swin(s, "eta", "X17")[0], swin(s, "non", "X17")[0]
    mb = (ETA * ex + cbf * bx) * NB2RAW
    x = DS[s]["win"]["X17"]
    P(f"{s:7s} {'  '.join(sh)}  C5(b)/model(b) in P {c5sh:.3f} | X17 window: model(b) {mb:.0f} rows, N/(CB+model) {x['N'] / (x['CB'] + mb):.3f}, "
      f"N_X/(CB+model) {RES[s]['x17']['NX'] / (x['CB'] + mb):.3f}")
open("/hdd3/przygoda/HADES/FAT_out/all_report/e1v/e1v_numbers.txt", "w").write("\n".join(out) + "\n")
