#!/usr/bin/env python3
# ============================================================================
# veestep_ana.py  --  INDEPENDENT verification of eestep (prefix veestep): data analysis
# ============================================================================
# Written 2026-11-15 by the independent verifier. Inputs: veestep_dump.C columns, edgeu sim columns (own loaders in
# veestep_core.py). Output: veestep_ana.log, veestep_ana.json in /hdd3/przygoda/HADES/FAT_out/all_report/veestep/.
# Sections
#   1  like-sign steps, all M(ee) (reproduction of the claims with an own estimator)
#   2  steps per M(ee) band: e+e+, e-e-, N, CB, N - CB (peak term fixed), SMASH x f, k_cross; bootstrap errors
#   3  opening-angle classes of the e-e- step; isBest; M(pipi) slices; control region 0.42-0.46
#   4  event content (PAT joins of veestep_dump.C: lpos_mult, e+e- entries of the event, twins, shared rings, n_lep)
#   5  charge flip: e-e- step at small M(ee), negative copies of the e+ (nflip), the effect on f
#   6  alignment: step in u, in M4 - E(m, 0), free position
#   7  peak window and X17 reach of the edge
# Run: python3 -I plots/all_report/veestep_ana.py
# ============================================================================
import json
import math
import sys

import numpy as np
import scipy.sparse as sp

sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import veestep_core as V

LOG = open(V.OUT + "veestep_ana.log", "w")
RES = {}


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    LOG.write(s + "\n")
    LOG.flush()


D = V.load_data()
S = V.load_sim_us()
okU = np.isfinite(D["u"]) & (np.abs(D["u"]) < 0.050)
inPP = (D["mpipi"] >= 0.30) & (D["mpipi"] < 0.40)
BASE = D["S0"] & okU & inPP
smU = (S["mpp"] >= 0.30) & (S["mpp"] < 0.40) & np.isfinite(S["u"])
BUMP = lambda m: V.sim_hist(S, smU & np.isin(S["g"], [0, 1]) & m)
PEAK = lambda m: V.sim_hist(S, smU & np.isin(S["g"], [2, 3]) & m)

# ---------------------------------------------------------------- bootstrap machinery (data)
rows = np.where(D["S0"] & okU & (D["mpipi"] >= 0.30) & (D["mpipi"] < 0.46))[0]
NREP = 300
K = V.boot_weights(D, rows, NREP)
pos = np.full(D["n"], -1, np.int64); pos[rows] = np.arange(rows.size)
sl_idx = np.full(D["n"], -1, np.int64)
for j, (a, b) in enumerate(V.SLICES):
    sl_idx[(D["mpipi"] >= a) & (D["mpipi"] < b)] = j
ubin = np.full(D["n"], -1, np.int64)
ub = np.floor((D["u"][okU] - V.EDGES[0]) / V.BW + 1e-9).astype(np.int64)
ubin[np.where(okU)[0]] = np.clip(ub, 0, V.NB - 1)


def rep_spectra(mask, reps=False):
    """nominal (and with reps=True replica) spectra N, PP, MM, CB (summed over slices, CB per slice) for mask & BASE"""
    m = mask & BASE
    r = np.where(m)[0]
    flat = (sl_idx[r] * 3 + D["k"][r].astype(np.int64)) * V.NB + ubin[r]
    A = sp.csr_matrix((np.ones(r.size), (np.arange(r.size), flat)), shape=(r.size, 4 * 3 * V.NB))
    H0 = np.asarray(A.sum(0)).reshape(4, 3, V.NB)
    if reps:
        HR = (sp.csr_matrix(K[:, pos[r]]) @ A).toarray().reshape(NREP, 4, 3, V.NB) if r.size else np.zeros((NREP, 4, 3, V.NB))
    else:
        HR = np.zeros((1, 4, 3, V.NB))

    def red(H):
        N = H[..., 0, :].sum(-2); PP = H[..., 1, :].sum(-2); MM = H[..., 2, :].sum(-2)
        nc = V.NB // 4
        Pc = H[..., 1, :].reshape(H.shape[:-2] + (nc, 4)); Mc = H[..., 2, :].reshape(H.shape[:-2] + (nc, 4))
        cbc = 2 * np.sqrt(Pc.sum(-1) * Mc.sum(-1))
        s = Pc + Mc; tot = s.sum(-1, keepdims=True)
        share = np.where(tot > 0, s / np.where(tot > 0, tot, 1), 0.25)
        CB = (cbc[..., None] * share).reshape(H.shape[:-2] + (V.NB,)).sum(-2)
        return N, PP, MM, CB
    return red(H0), red(HR)


def sd(a):
    a = np.asarray(a, float)
    return float(np.std(a, ddof=1))


# ---------------------------------------------------------------- nominal estimator (all M(ee))
(N0, PP0, MM0, CB0), _ = rep_spectra(np.ones(D["n"], bool))
r0 = V.fit_pois(N0)
var0 = np.maximum(np.interp(V.XC, V.XC[V.FIT], r0["mu"]), 1) + CB0
L0 = V.Lin(var0)
bump0, peak0 = BUMP(np.ones(S["u"].size, bool)), PEAK(np.ones(S["u"].size, bool))
hd0 = L0.h(N0 - CB0 - V.A_P * peak0); hs0 = L0.h(bump0)
F0 = hd0 / hs0
P("# veestep_ana.py: independent verification of eestep, 2026-11-15")
P("rows: data %d; S0, |u| < 50 MeV, 0.30 <= M(pipi) < 0.40: us %d pp %d mm %d" %
  (D["n"], (BASE & (D["k"] == 0)).sum(), (BASE & (D["k"] == 1)).sum(), (BASE & (D["k"] == 2)).sum()))
P("join checks: like-sign rows joined %d of %d" % (D["joinok"][D["k"] > 0].sum(), (D["k"] > 0).sum()))

# ---------------------------------------------------------------- 1
P("\n== 1. like-sign steps, all M(ee) (Poisson ML, pol2 + step, sigma 6.7 MeV, -45..45 MeV)")
for nm, y in [("e-e-", MM0), ("e+e+", PP0)]:
    r = V.fit_pois(y); r0n = V.fit_pois(y, nostep=True)
    P("  %s: rows %d, h = %.1f +- %.1f, -2lnL %.1f (ndf %d), delta(-2lnL) no step %.1f" % (nm, y[V.FIT].sum(), r["h"], r["eh"], r["chi2"], r["ndf"], r0n["chi2"] - r["chi2"]))
    RES["ls_all_" + nm] = [r["h"], r["eh"]]
P("  N alone: h = %.1f; CB: h = %.1f (linear estimator); peak term a_P K: %.1f" % (L0.h(N0), L0.h(CB0), L0.h(V.A_P * peak0)))
P("  f (own estimator) = h_data / h_sim = %.1f / %.1f = %.4f (nominal edgeu 0.4200)" % (hd0, hs0, F0))
RES["f_own"] = F0

# ---------------------------------------------------------------- 2
P("\n== 2. per M(ee) band (CB of each band from its own like-sign rows; errors: %d event-group bootstrap replicas)" % NREP)
bands = [(0, 10), (10, 30), (30, 50), (50, 70), (70, 90), (90, 110), (110, 135), (0, 50), (90, 135), (70, 135), (50, 135), (100, 135)]
P("  %-9s %13s %13s %9s %8s %14s %8s %14s %9s" % ("M(ee)", "e+e+", "e-e-", "N", "CB", "N-CB", "SMASHxf", "k_cross", "z(artef)"))
RES["bands"] = {}
for (a, b) in bands:
    me = D["mee"] * 1000
    m = (me >= a) & (me < b)
    msim = (S["mee"] * 1000 >= a) & (S["mee"] * 1000 < b)
    (N, PPs, MMs, CB), (NR, PPR, MMR, CBR) = rep_spectra(m, True)
    rp = V.fit_pois(PPs); rm = V.fit_pois(MMs)
    pk = V.A_P * PEAK(msim); bu = BUMP(msim)
    hN = L0.h(N - pk); hCB = L0.h(CB); hNC = hN - hCB
    hNR = np.array([L0.h(NR[i] - pk) for i in range(NREP)]); hCBR = np.array([L0.h(CBR[i]) for i in range(NREP)])
    hNCR = hNR - hCBR
    pred = F0 * L0.h(bu)
    # sim stat error of the prediction: event-group Poisson bootstrap of the sim rows
    kcr = (hN - pred) / hCB if abs(hCB) > 1e-9 else float("nan")
    kcrR = (hNR - pred) / np.where(np.abs(hCBR) > 1e-9, hCBR, np.nan)
    e_nc = sd(hNCR)
    z_art = (hNC - (pred - hCB)) / e_nc
    z_cb = (hNC - pred) / e_nc
    lab = "%d-%d" % (a, b)
    P("  %-9s %6.1f+-%5.1f %6.1f+-%5.1f %9.1f %8.1f %7.1f+-%5.1f %8.1f %7.2f+-%5.2f %9.2f   (z CB rule %.2f; boot e+e+ %.1f e-e- %.1f)" %
      (lab, rp["h"], rp["eh"], rm["h"], rm["eh"], hN, hCB, hNC, e_nc, pred, kcr, float(np.nanstd(kcrR, ddof=1)), z_art, z_cb,
       sd([L0.h(PPR[i]) for i in range(NREP)]), sd([L0.h(MMR[i]) for i in range(NREP)])))
    RES["bands"][lab] = dict(pp=[rp["h"], rp["eh"]], mm=[rm["h"], rm["eh"]], N=hN, eN=sd(hNR), CB=hCB, eCB=sd(hCBR), NmCB=[hNC, e_nc],
                             pred=pred, k=[kcr, float(np.nanstd(kcrR, ddof=1))], z_art=z_art, z_cb=z_cb,
                             ratio_pp_mm=rp["h"] / rm["h"] if rm["h"] else None, rows_mm=float(MMs[V.FIT].sum()),
                             rows_pp=float(PPs[V.FIT].sum()))
for lab in ["90-135", "70-135", "100-135"]:
    x = RES["bands"][lab]
    rr = x["pp"][0] / x["mm"][0]
    er = abs(rr) * math.hypot(x["pp"][1] / x["pp"][0], x["mm"][1] / x["mm"][0])
    P("  e+e+/e-e- step ratio at M(ee) %s: %.3f +- %.3f" % (lab, rr, er))
    RES["ratio_" + lab] = [rr, er]

# CB rule non-linearity: the CB response to a like-sign component with a ++/-- ratio different from the background
P("\n  CB-rule response to a cross-pair step with ++/-- ratio rc on a background with ratio rB (cells of the M(ee) >= 90 band):")
m = (D["mee"] * 1000 >= 90) & (D["mee"] * 1000 < 135)
(Nh, PPh, MMh, CBh), _ = rep_spectra(m)
rB = PPh[V.FIT].sum() / MMh[V.FIT].sum()
for rc in [0.08, 0.13, 0.30]:
    resp = math.sqrt(rB) + rc / math.sqrt(rB); true = 2 * math.sqrt(rc)
    P("    rB = %.3f, rc = %.2f: CB step = %.3f x (e-e- step), factorised partners 2 sqrt(rc) = %.3f: ratio %.3f" % (rB, rc, resp, true, resp / true))
RES["rB_90"] = rB

# ---------------------------------------------------------------- 3
P("\n== 3. e-e- step by opening angle, isBest, M(pipi) slice (all M(ee), S0, 0.30-0.40)")
oa = D["oa"]
for (a, b) in [(0, 10), (10, 20), (20, 35), (35, 60), (60, 180), (0, 20), (20, 60)]:
    m = (oa >= a) & (oa < b)
    (N, PPs, MMs, CB), _ = rep_spectra(m)
    rm = V.fit_pois(MMs); rp = V.fit_pois(PPs)
    P("  OA %3d-%3d deg: e-e- %6.1f +- %5.1f (rows %5d) | e+e+ %5.1f +- %4.1f | N - CB - peak %7.1f" %
      (a, b, rm["h"], rm["eh"], MMs[V.FIT].sum(), rp["h"], rp["eh"], L0.h(N - CB)))
    RES["oa_%d_%d" % (a, b)] = [rm["h"], rm["eh"]]
# OA of the e-e- step rows for M(ee) >= 90 only
for (a, b) in [(0, 20), (20, 35), (35, 60), (60, 180)]:
    m = (oa >= a) & (oa < b) & (D["mee"] >= 0.090)
    (N, PPs, MMs, CB), _ = rep_spectra(m)
    rm = V.fit_pois(MMs)
    P("  M(ee) >= 90, OA %3d-%3d: e-e- %6.1f +- %5.1f (rows %d)" % (a, b, rm["h"], rm["eh"], MMs[V.FIT].sum()))
for nm, m in [("isBest==1", (D["bits"] & 4) != 0), ("isBest!=1", (D["bits"] & 4) == 0)]:
    (N, PPs, MMs, CB), _ = rep_spectra(m)
    rm = V.fit_pois(MMs)
    P("  %s: e-e- %.1f +- %.1f (rows %d)" % (nm, rm["h"], rm["eh"], MMs[V.FIT].sum()))
for (a, b) in [(0.30, 0.34), (0.34, 0.37), (0.37, 0.40)]:
    m = (D["mpipi"] >= a) & (D["mpipi"] < b)
    (N, PPs, MMs, CB), _ = rep_spectra(m)
    rm = V.fit_pois(MMs); rm9 = V.fit_pois(rep_spectra(m & (D["mee"] >= 0.070))[0][2])
    P("  M(pipi) %.2f-%.2f: e-e- %.1f +- %.1f; M(ee) >= 70: %.1f +- %.1f" % (a, b, rm["h"], rm["eh"], rm9["h"], rm9["eh"]))
# control region 0.42-0.46 with a pseudo edge E(M(pipi) - 0.10, q)
with np.errstate(invalid="ignore"):
    up = D["m4"] - V.edge(D["mpipi"] - 0.10, D["mee"])
cr = D["S0"] & (D["mpipi"] >= 0.42) & (D["mpipi"] < 0.46) & (D["mee"] < V.MPi0) & (D["k"] == 2)
for nm, mm in [("all M(ee)", cr), ("M(ee) 50-70 (M4 threshold allows the window)", cr & (D["mee"] >= 0.050) & (D["mee"] < 0.070))]:
    rm = V.fit_pois(V.hist(up[mm]))
    P("  control M(pipi) 0.42-0.46, pseudo edge E(m - 0.10, q), %s: e-e- %.1f +- %.1f" % (nm, rm["h"], rm["eh"]))

# ---------------------------------------------------------------- 4
P("\n== 4. event content of the e-e- rows (all M(ee) and M(ee) >= 90; S0, 0.30-0.40)")
LS = D["k"] == 2
P("  join: e-e- rows with PAT values %d of %d; e+e- entries of the event found by key blocks (0 non-adjacent keys)" %
  ((D["joinok"] & LS).sum(), LS.sum()))
cats = [("all", np.ones(D["n"], bool)), ("lpos_mult == 0", D["lpos"] == 0), ("lpos_mult >= 1", D["lpos"] >= 1),
        ("no e+e- entry in event", D["nus"] == 0), ("e+e- entry in event", D["nus"] > 0),
        ("twin (same pi+ pi- and one e-)", D["ntwin"] > 0), ("e+e- entries but no twin", (D["nus"] > 0) & (D["ntwin"] == 0)),
        ("lneg_mult == 2", D["lneg"] == 2), ("lneg_mult >= 3", D["lneg"] >= 3), ("shared ring", D["sr"] != 0),
        ("n_lep == 2 (group)", D["nlep"] == 2), ("n_lep >= 3 (group)", D["nlep"] >= 3)]
for nm, c in cats:
    out = []
    for mq in [np.ones(D["n"], bool), D["mee"] >= 0.090]:
        (N, PPs, MMs, CB), _ = rep_spectra(c & mq)
        rm = V.fit_pois(MMs)
        out.append((rm["h"], rm["eh"], MMs[V.FIT].sum()))
    P("  %-34s all M(ee): %6.1f +- %5.1f (rows %5d) | M(ee) >= 90: %6.1f +- %5.1f (rows %5d)" % ((nm,) + out[0] + out[1]))
    RES["ev_" + nm] = out
sel = BASE & LS
P("  share of e-e- rows (S0, 0.30-0.40, |u| < 45 MeV) with lpos_mult == 0: %.3f; without e+e- entry: %.3f" %
  ((D["lpos"][sel & (np.abs(D["u"]) < 0.045)] == 0).mean(), (D["nus"][sel & (np.abs(D["u"]) < 0.045)] == 0).mean()))
# e+e+ rows: negative leptons in the event
LP = D["k"] == 1
(N, PPs, MMs, CB), _ = rep_spectra(D["lneg"] == 0)
rp = V.fit_pois(PPs)
P("  e+e+ rows with lneg_mult == 0: step %.1f +- %.1f (rows %d); share of e+e+ rows with lneg_mult == 0: %.3f" %
  (rp["h"], rp["eh"], PPs[V.FIT].sum(), (D["lneg"][BASE & LP] == 0).mean()))

# ---------------------------------------------------------------- 5
P("\n== 5. charge flip")
me = D["mee"]
(Nl, PPl, MMl, CBl), (NlR, PPlR, MMlR, CBlR) = rep_spectra(me < 0.050, True)
rm = V.fit_pois(MMl); rp = V.fit_pois(PPl)
msim = S["mee"] < 0.050
hNl = L0.h(Nl - V.A_P * PEAK(msim) - CBl)
rflip = rm["h"] / hNl
P("  M(ee) < 50: e-e- step %.1f +- %.1f, e+e+ step %.1f +- %.1f, bump step of N - CB - peak %.1f" % (rm["h"], rm["eh"], rp["h"], rp["eh"], hNl))
P("  flip rate per bump pair (e+ -> e-): r = %.4f +- %.4f; 90 %% CL upper limit (Gaussian, r >= 0 not imposed) %.4f" %
  (rflip, rm["eh"] / hNl, (rm["h"] + 1.2816 * rm["eh"]) / hNl))
# effect on f: a flipped pair leaves N (lost bump row) and adds an e-e- row (raises the CB)
rUL = (rm["h"] + 1.2816 * rm["eh"]) / hNl
dCB = math.sqrt(PPl[V.FIT].sum() / MMl[V.FIT].sum())   # d CB / d N-- for a small e-e- addition at the background ratio
P("  at r = %.4f: lost N rows %.2f %% of the bump; CB rise %.2f %% (dCB/dN-- = sqrt(N++/N--) = %.3f); f higher by about %.2f %%" %
  (rUL, 100 * rUL, 100 * rUL * dCB, dCB, 100 * rUL * (1 + dCB)))
RES["flip"] = dict(r=rflip, er=rm["eh"] / hNl, ul=rUL, df_ul=rUL * (1 + dCB))
# negative copies of the e+ of e+e- rows (nflip) and e-e- entries in the event
US = D["k"] == 0
for nm, c in [("e+e- rows with a negative copy of the e+ (nflip > 0)", D["nflip"] > 0), ("e+e- rows with e-e- entries in the event", D["nls"] > 0)]:
    (N, PPs, MMs, CB), _ = rep_spectra(c & (D["nls"] >= 0))
    r = V.fit_pois(N)
    P("  %s: step of N %.1f +- %.1f (rows %d) = %.4f of the step of all e+e- rows" % (nm, r["h"], r["eh"], N[V.FIT].sum(), r["h"] / L0.h(N0)))
P("  e+e- rows of the edge region with PAT join: %d; with nflip > 0: %d" % (((D["nls"] >= 0) & US).sum(), ((D["nflip"] > 0) & US).sum()))

# ---------------------------------------------------------------- 6
P("\n== 6. alignment of the e-e- step (M(ee) >= 70, S0, 0.30-0.40)")
mm = BASE & (D["k"] == 2) & (D["mee"] >= 0.070)
with np.errstate(invalid="ignore"):
    alt = {"u = M4 - E(m, q)": D["u"], "M4 - E(m, 0)": D["m4"] - V.edge(D["mpipi"], 0 * D["mee"]),
           "u - 10 MeV": D["u"] - 0.010, "u + 10 MeV": D["u"] + 0.010}
for nm, x in alt.items():
    sel = D["S0"] & inPP & (D["k"] == 2) & (D["mee"] >= 0.070) & (D["mee"] < V.MPi0)
    y = V.hist(x[sel])
    r = V.fit_pois(y); r0n = V.fit_pois(y, nostep=True)
    P("  %-18s step %6.1f +- %5.1f, delta(-2lnL) %.1f" % (nm, r["h"], r["eh"], r0n["chi2"] - r["chi2"]))
y = V.hist(D["u"][mm])
prof = []
for u0 in np.arange(-0.010, 0.0101, 0.0025):
    r = V.fit_pois(y, u0=u0)
    prof.append((u0 * 1000, r["chi2"], r["h"]))
c0 = min(p[1] for p in prof)
P("  free position: " + ", ".join("u0 %+.1f: d %.1f (h %.0f)" % (a, b - c0, h) for a, b, h in prof))
for sg in [0.004, 0.0067, 0.010, 0.013]:
    r = V.fit_pois(y, sigma=sg)
    P("  width %.1f MeV: h %.1f +- %.1f, -2lnL %.1f" % (sg * 1000, r["h"], r["eh"], r["chi2"]))

# ---------------------------------------------------------------- 7
P("\n== 7. reach of the edge into the peak window (P 0.53364-0.55416) and the X17 band")
for m in [0.28, 0.30, 0.35, 0.40]:
    qs = np.arange(0, 0.135, 0.0005)
    e = V.edge(np.full(qs.size, m), qs)
    q1 = qs[np.argmax(e >= 0.53364)] if np.any(e >= 0.53364) else float("nan")
    P("  M(pipi) %.2f: E(m, 0) %.5f; E >= 0.53364 for M(ee) >= %.1f MeV; E(m, 24.5 MeV) %.5f" % (m, e[0], q1 * 1000, V.edge(m, 0.0245)))
P("  lowest M(pipi) of the data rows 0.26 (2 m_pi = 0.2791): E(0.2791, 0.0245) = %.5f" % V.edge(0.2791, 0.0245))
json.dump(RES, open(V.OUT + "veestep_ana.json", "w"), indent=1, default=float)
P("done")
