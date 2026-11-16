#!/usr/bin/env python3
# ============================================================================
# edgeu_toys.py  --  package edgeu: closure of the u-edge estimator on pseudo-data with a known truth
# ============================================================================
# Written 2026-11-15 (prefix edgeu). The WP4 recipe (wp4_lib.h: truth = data CB + c x non-eta background x deformation +
# eta classes x scales; toys = Poisson N, Poisson like-sign, CB rebuilt) applied to the u spectra of the nominal estimator:
#   truth per M(pipi) slice s (0.30-0.34-0.37-0.39-0.40), expected raw rows per 2.5 MeV u bin:
#     mu_s = CB_s                                  the data CB of the slice (2 sqrt(N++ N--) of the data like-sign rows)
#          + c B(u) [1 + b x + q x^2] [1 + A G]    non-eta background G_D + G_C of SMASH summed over the slices (only the total
#                                                  N enters the fit), smoothed by exp(pol4) on -100..+100 MeV (pol3/5/6 and the
#                                                  raw spectrum as variants); the deformation factors are the running mean of
#                                                  the row factor in each u bin, x = (M4' - 0.55)/0.1 as D11,
#                                                  G = exp(-(M4' - m0)^2/2w^2); stage A: 0.98 (1 - 0.064 x)(G_D + 2.93 G_C)
#          + a [C2c + kappa C3] + a G_eta           bump (C2 with the gamma* factor) and other-eta rows (smoothed)
#          + a_Pt [kappa C4 + C5]                   eta peak, scale a_Pt (nominal 0.418 = the a_P the estimator subtracts)
#   toy: N_s = Poisson(mu_s), N++_s, N--_s = Poisson(data like-sign counts), CB rebuilt per slice by the rule; the
#   nominal estimator (pol2 + step, -45..+45 MeV, sigma 6.70 MeV, a_P x peak template fixed);
#     AA: h_sim of the full simulated sample (the truth's own templates: the estimator alone);
#     RF: h_sim and the peak template of the WP0 bootstrap replica r = toy mod 100 (template noise included).
#   Pull = (f - truth)/err, err = stat (fit, x 0.98 = bootstrap/fit of the data) for AA, (+) MC (replica rms) for RF.
# Linear predictions (expected shift = projection of the truth change with the estimator of the data fit) are printed next to
# the toy means. Outputs: /hdd3/.../edgeu/edgeu_toys.json, edgeu_toys.log.
# Run (from /home/przygoda/HADES/FAT/FAT; about 6 GB, a few minutes):  python3 -I plots/all_report/edgeu_toys.py [ntoy]
# ============================================================================
import json
import math
import sys
import time

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import edgeu_core as C

NTOY = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
T0 = time.time()
LOG = open(C.OUT + "edgeu_toys.log", "w")
def P(*a):
    s = " ".join(str(x) for x in a)
    print(s); LOG.write(s + "\n"); LOG.flush()
ANA = json.load(open(C.OUT + "edgeu_ana.json"))
SIG = ANA["reso"]["sigma_core"]
MC_REL = ANA["nom"]["mc"] / ANA["nom"]["f"]
BOOT_FIT = ANA["nom"]["statboot"] / ANA["nom"]["statfit"]
D = C.load_data(); S = C.load_sim()
P("# edgeu_toys.py, 2026-11-15: %d toys per truth; sigma %.2f MeV; MC rel %.4f; bootstrap/fit %.3f" % (NTOY, 1e3 * SIG, MC_REL, BOOT_FIT))
g = C.GRID_NOM
SL = C.slices_of(0.30, 0.40)
T = C.step_gauss(g, SIG)
LO, HI = -0.045, 0.045

# ---------------------------------------------------------------- data CB and like-sign counts per slice
ud = C.u_data(D)
Hd = C.hist_data(D, ud, (D["bits"] & 1) != 0)
LS = [(Hd[1, a:b].sum((0, 1)), Hd[2, a:b].sum((0, 1))) for a, b in SL]
CBs = [C.cb_rule(p, m, g)[0] for p, m in LS]
Nd = sum(Hd[0, a:b].sum((0, 1)) for a, b in SL)

# ---------------------------------------------------------------- simulated components per slice
def sim_hists(s=C.S_NOM, r=C.R_NOM, fac=None, groups=None):
    u, mpp, m4p = C.u_sim(S, s, r)
    w = C.sim_weight(S, c2corr=True)
    if fac is not None:
        w = w * fac(m4p, mpp)
    return C.hist_sim(S, u, mpp, w, g, groups=groups), m4p, mpp
Hs, M4P, MPP = sim_hists()
def comp(H, names, sl):
    a, b = sl
    return sum(H[C.G[n], a:b].sum((0, 1)) for n in names)
BUMP = sum(comp(Hs, ["C2", "C3"], sl) for sl in SL)
PEAK = sum(comp(Hs, ["C4", "C5"], sl) for sl in SL)
# non-eta background: the SUM over the slices (only the total N enters the fit; the CB stays per slice), smoothed by
# exp(pol n) on -100..+100 MeV (adequate chi2 there, edgeu_bgsmooth.py); deformations = smooth background x the running mean
# (5 bins) of the row-factor mean in each u bin
SM_LO, SM_HI = -0.100, 0.100
def smooth_exp(y, lo=SM_LO, hi=SM_HI, order=4):
    m = (g.x >= lo) & (g.x < hi)
    xx = (g.x[m] - 0.5 * (lo + hi)) / (0.5 * (hi - lo)); yy = np.maximum(y[m], 0)
    def nll(p):
        mu = np.exp(np.clip(np.polyval(p[::-1], xx), -50, 50))
        return float(np.sum(mu - yy * np.log(mu + 1e-300)))
    p = np.array([math.log(max(yy.mean(), 1e-6))])
    for o in range(1, order + 1):
        p = minimize(nll, np.r_[p, 0.0], method="BFGS").x
    out = np.zeros(g.nb); out[m] = np.exp(np.polyval(p[::-1], xx))
    out[~m] = y[~m]   # outside the smoothing range (never inside the fit range) the raw spectrum
    return out
def runmean(r, n=5):
    k = np.ones(n) / n
    return np.convolve(np.r_[np.full(n // 2, r[0]), r, np.full(n // 2, r[-1])], k, mode="valid")
FAC = {"x": lambda m4, mp: (m4 - 0.55) / 0.1, "x2": lambda m4, mp: ((m4 - 0.55) / 0.1) ** 2,
       "g490": lambda m4, mp: np.exp(-0.5 * ((m4 - 0.490) / 0.015) ** 2), "g505": lambda m4, mp: np.exp(-0.5 * ((m4 - 0.505) / 0.008) ** 2),
       "y": lambda m4, mp: (mp - 0.35) / 0.05}
HF = {k: sim_hists(fac=f, groups=[5, 6, 7, 8])[0] for k, f in FAC.items()}
NONETA = ["GDom", "GDnom", "GCom", "GCnom"]
PARTS = {"all": NONETA, "GD": ["GDom", "GDnom"], "GC": ["GCom", "GCnom"], "om": ["GDom", "GCom"], "nom": ["GDnom", "GCnom"]}
BGS = {}
raw_all = sum(comp(Hs, NONETA, sl) for sl in SL)
RATIO = {}
for k in FAC:
    rk = sum(comp(HF[k], NONETA, sl) for sl in SL)
    RATIO[k] = runmean(np.where(raw_all > 0, rk / np.where(raw_all > 0, raw_all, 1), 0.))
for order in (3, 4, 5, 6):
    for part, names in PARTS.items():
        raw = sum(comp(Hs, names, sl) for sl in SL)
        BGS[(part, order)] = smooth_exp(raw, order=order)
        BGS[(part, "raw")] = raw
    BGS[("eta", order)] = smooth_exp(sum(comp(Hs, ["Geta"], sl) for sl in SL), order=order)
BGS[("eta", "raw")] = sum(comp(Hs, ["Geta"], sl) for sl in SL)
BG = []
for sl in SL:
    d = {}
    d["C2"] = comp(Hs, ["C2"], sl); d["C3"] = comp(Hs, ["C3"], sl); d["C4"] = comp(Hs, ["C4"], sl); d["C5"] = comp(Hs, ["C5"], sl)
    BG.append(d)
# mass-transfer variant of the eta truth (s + 0.0009, r = 0.89)
HsM, _, _ = sim_hists(C.S_NOM + C.S_ERR, 0.89)
for i, sl in enumerate(SL):
    for n in ("C2", "C3", "C4", "C5"):
        BG[i][n + "_mass"] = comp(HsM, [n], sl)
# replicas of the bump and the peak (RF)
rows = np.where(S["g"] <= 3)[0]
K = C.sim_replica_weights(S, rows, 100)
u0, mpp0, _ = C.u_sim(S)
w0 = C.sim_weight(S, c2corr=True)
Ssub = {kk: S[kk][rows] for kk in ("g", "mee")}
HR = C.hist_sim(Ssub, u0[rows], mpp0[rows], w0[rows], g, reps=K)
BUMP_R = np.array([sum(comp(HR[r], ["C2", "C3"], sl) for sl in SL) for r in range(100)])
PEAK_R = np.array([sum(comp(HR[r], ["C4", "C5"], sl) for sl in SL) for r in range(100)])
del K, HR
P("   inputs ready (%.0f s)" % (time.time() - T0))

# ---------------------------------------------------------------- truth
def truth(sc):
    """per-slice expected rows: slice 0 carries the whole background (only the total N enters the fit)"""
    order = "raw" if sc.get("raw") else sc.get("order", 4)
    def bgpart(part):
        B = BGS[(part, order)]
        out = B.copy()
        for key, k in (("b", "x"), ("q", "x2"), ("ty", "y")):
            if sc.get(key, 0):
                out = out + sc[key] * B * RATIO[k]
        if sc.get("A", 0):
            out = out + sc["A"] * B * RATIO[sc["gk"]]
        return out
    if sc.get("stageA"):
        bkg = 0.98 * (1 + sc.get("b", 0) * RATIO["x"]) * (bgpart("GD") + 2.93 * bgpart("GC"))
    elif "fom" in sc:
        bkg = sc.get("c", 1.35) * (sc["fom"] * bgpart("om") + sc["fnom"] * bgpart("nom"))
    else:
        bkg = sc.get("c", 1.35) * bgpart("all")
    a, kap = sc.get("a", 0.42), sc.get("kap", 1.0)
    mus = []
    for i, sl in enumerate(SL):
        d = BG[i]
        mt = "_mass" if sc.get("mass") else ""
        eta = a * (d["C2" + mt] + kap * d["C3" + mt])
        pk = sc.get("aPt", C.AP) * (kap * d["C4" + mt] + d["C5" + mt])
        mu = CBs[i] + eta + pk
        if i == 0:
            mu = mu + bkg + a * BGS[("eta", order)]
        mus.append(np.maximum(mu, 0.))
    return mus
def estimator(Ns, pps, mms, peak=PEAK, Vfix=None):
    N = np.zeros(g.nb); CB = np.zeros(g.nb); COV = np.zeros((g.nb, g.nb))
    for n, p, m in zip(Ns, pps, mms):
        cb, cov = C.cb_rule(p, m, g)
        N += n; CB += cb; COV += cov
    sp = C.Spec(N, CB, COV)
    return C.fit_step(sp, T, g, LO, HI, 2, fixed=C.AP * peak, Vfix=Vfix)
# reference: the estimator on the data
R0 = estimator([Hd[0, a:b].sum((0, 1)) for a, b in SL], [p for p, m in LS], [m for p, m in LS])
HS0 = C.h_lin(R0, BUMP)
P("   data: h %.2f, h_sim %.2f, f %.4f (edgeu_ana: %.4f)" % (R0.h, HS0, R0.h / HS0, ANA["nom"]["f"]))

SC = [
    dict(tag="a", text="WP4 (a): non-eta x 1.35, tilt b = -0.15, a = 0.42 (smoothed SMASH background)", b=-0.15),
    dict(tag="c1", text="non-eta x 1, no tilt, a = 0.42", c=1.0),
    dict(tag="a1", text="(a) with a = 1", b=-0.15, a=1.0),
    dict(tag="null", text="(a) with a = 0 (no eta at all)", b=-0.15, a=0.0, aPt=0.0),
    dict(tag="q-", text="(a) + q = -0.173 (the B4-killing quadratic term)", b=-0.15, q=-0.173),
    dict(tag="q--", text="(a) + q = -0.346", b=-0.15, q=-0.346),
    dict(tag="q+", text="(a) + q = +0.173", b=-0.15, q=0.173),
    dict(tag="om", text="WP4 (e): omega/rho0-pion background x 0.3, the rest x 1.6", b=-0.15, c=1.0, fom=0.3, fnom=1.6),
    dict(tag="ty", text="(a) + M(pipi) tilt 1 + 0.3 (M(pipi) - 0.35)/0.05 of the non-eta background", b=-0.15, ty=0.3),
    dict(tag="wide", text="(a) x (1 + 0.10 Gauss(0.490, 15 MeV)) in M4", b=-0.15, A=0.10, gk="g490"),
    dict(tag="narrow", text="(a) x (1 + 0.05 Gauss(0.505, 8 MeV)) in M4", b=-0.15, A=0.05, gk="g505"),
    dict(tag="raw", text="(a) with the RAW SMASH background and G_eta (template noise in the truth)", b=-0.15, raw=True),
    dict(tag="mass", text="(a) with the eta truth at s + 0.0009, r = 0.89 (data-like resolution; templates nominal)", b=-0.15, mass=True),
    dict(tag="kap", text="(a) with conversions (C3, C4) x 0.48 in the truth", b=-0.15, kap=0.48),
    dict(tag="pk30", text="(a) with the true peak scale 0.30 (the estimator subtracts 0.418)", b=-0.15, aPt=0.30),
    dict(tag="pk55", text="(a) with the true peak scale 0.55", b=-0.15, aPt=0.55),
    dict(tag="datalike", text="data-like: a = 0.42, conversions x 0.48, non-eta x 1.16, no tilt", c=1.16, kap=0.48),
    dict(tag="sA", text="stage-A-like background: 0.98 (1 - 0.064 x)(G_D + 2.93 G_C), a = 0.42", stageA=True, b=-0.064),
    dict(tag="sm3", text="(a) with the background smoothed by exp(pol3)", b=-0.15, order=3),
    dict(tag="sm5", text="(a) with the background smoothed by exp(pol5)", b=-0.15, order=5),
    dict(tag="sm6", text="(a) with the background smoothed by exp(pol6)", b=-0.15, order=6),
]
rng = np.random.default_rng(0xED6E7015)
OUTS = {}
mu_a = None
for sc in SC:
    t1 = time.time()
    mus = truth(sc)
    a, kap = sc.get("a", 0.42), sc.get("kap", 1.0)
    # truth of the estimated quantity: f (= a) for kappa = 1; for kappa != 1 the kappa-0.48 reading closes at a
    lin2c = C.h_lin(R0, sum(BG[i]["C2"] for i in range(len(SL)))); lin3 = C.h_lin(R0, sum(BG[i]["C3"] for i in range(len(SL))))   # C2 already x gamma* factor
    if kap != 1.0:
        ftrue = a * (lin2c + kap * lin3) / (lin2c + lin3)   # what f measures when the conversions are scaled
    else:
        ftrue = a
    # Asimov (linear prediction): the estimator on the expected spectrum with expected like-sign
    Ra = estimator(mus, [p for p, m in LS], [m for p, m in LS], Vfix=None)
    f_asi = Ra.h / HS0
    if sc["tag"] == "a":
        mu_a = mus; f_a_asi = f_asi
    fa, fr, pa, pr = [], [], [], []
    for it in range(NTOY):
        Ns = [rng.poisson(m) for m in mus]
        pps = [rng.poisson(p) for p, m in LS]; mms = [rng.poisson(m) for p, m in LS]
        R = estimator(Ns, pps, mms)
        hs_aa = C.h_lin(R, BUMP)
        r = it % 100
        Rr = estimator(Ns, pps, mms, peak=PEAK_R[r], Vfix=R.V)
        hs_rf = C.h_lin(R, BUMP_R[r])
        f1 = R.h / hs_aa; f2 = Rr.h / hs_rf
        e1 = R.eh / hs_aa * BOOT_FIT
        fa.append(f1); fr.append(f2)
        if ftrue > 0:
            pa.append((f1 - ftrue) / e1); pr.append((f2 - ftrue) / math.hypot(e1, MC_REL * f2))
        else:
            pa.append(f1 / e1); pr.append(f2 / math.hypot(e1, MC_REL * 0.42))
    fa, fr, pa, pr = map(np.array, (fa, fr, pa, pr))
    o = dict(text=sc["text"], ftrue=ftrue, asimov=f_asi, mean_aa=float(fa.mean()), rms_aa=float(fa.std()), mean_rf=float(fr.mean()), rms_rf=float(fr.std()),
             pull_sd_aa=float(pa.std()), pull_mean_aa=float(pa.mean()), pull_sd_rf=float(pr.std()), pull_mean_rf=float(pr.mean()),
             cov_aa=float(np.mean(np.abs(pa) < 1)), cov_rf=float(np.mean(np.abs(pr) < 1)))
    if ftrue > 0:
        o["bias_aa"] = 100 * (o["mean_aa"] / ftrue - 1); o["bias_rf"] = 100 * (o["mean_rf"] / ftrue - 1); o["bias_asimov"] = 100 * (f_asi / ftrue - 1)
    else:
        o["bias_aa"] = o["mean_aa"]; o["bias_rf"] = o["mean_rf"]; o["bias_asimov"] = f_asi
    OUTS[sc["tag"]] = o
    P("   %-8s truth %.4f | Asimov %.4f (%+.2f %%) | AA mean %.4f (%+.2f %%) rms %.4f pull %.2f/%.2f cov %.2f | RF mean %.4f (%+.2f %%) rms %.4f pull %.2f/%.2f cov %.2f  [%s] %.0f s" %
      (sc["tag"], ftrue, f_asi, o["bias_asimov"], o["mean_aa"], o["bias_aa"], o["rms_aa"], o["pull_mean_aa"], o["pull_sd_aa"], o["cov_aa"],
       o["mean_rf"], o["bias_rf"], o["rms_rf"], o["pull_mean_rf"], o["pull_sd_rf"], o["cov_rf"], sc["text"], time.time() - t1))
# shifts relative to (a), Asimov and toys
for k, o in OUTS.items():
    if OUTS["a"]["ftrue"] > 0 and o["ftrue"] > 0:
        o["shift_vs_a_asimov"] = 100 * (o["asimov"] / o["ftrue"] - OUTS["a"]["asimov"] / OUTS["a"]["ftrue"])
        o["shift_vs_a_aa"] = o["bias_aa"] - OUTS["a"]["bias_aa"]
# free smearing on the Asimov of scenario (a): what the data fit's free sigma should give if the data resolution = simulation
sg = np.arange(0.0040, 0.0111, 0.0002); ch = []
for s_ in sg:
    Tq = C.step_gauss(g, s_)
    N = sum(mu_a); CB = sum(CBs)
    sp = C.Spec(N, CB, np.zeros((g.nb, g.nb)))
    R = C.fit_step(sp, Tq, g, LO, HI, 2, fixed=C.AP * PEAK, Vfix=np.diag(np.maximum(Nd[R0.idx], 1.0)))
    ch.append(R.chi2)
ch = np.array(ch); i = int(np.argmin(ch)); pf = np.polyfit(sg[max(0, i - 6):i + 7], ch[max(0, i - 6):i + 7], 2)
OUTS["asimov_free_sigma"] = dict(sigma=float(-pf[1] / (2 * pf[0])), err=float(1 / math.sqrt(pf[0])) if pf[0] > 0 else None)
P("   free smearing on the Asimov of (a) (data-sized errors): sigma = %.2f +- %.2f MeV (data: %.2f +- %.2f MeV)" %
  (1e3 * OUTS["asimov_free_sigma"]["sigma"], 1e3 * OUTS["asimov_free_sigma"]["err"], 1e3 * ANA["free_sigma"]["data"]["sigma"], 1e3 * ANA["free_sigma"]["data"]["err"]))
json.dump(OUTS, open(C.OUT + "edgeu_toys.json", "w"), indent=1)
P("done in %.0f s" % (time.time() - T0))
