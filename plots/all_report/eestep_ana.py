#!/usr/bin/env python3
# ============================================================================
# eestep_ana.py  --  package eestep: where do the e-e- rows of the step at u = 0 come from?
# ============================================================================
# Written 2026-11-15 (prefix eestep). Inputs (read-only):
#   - /hdd3/przygoda/HADES/FAT_out/all_report/eestep/ls_*.f32, usr_*.f32 (eestep_dump.C: data like-sign rows and e+e- rows
#     of the edge region with PAT track details and event links)
#   - /hdd3/przygoda/HADES/FAT_out/all_report/eestep/sus_*.f32, sls_*.f32 (eestep_simdump.C: simulated rows with truth and
#     like-sign twins); simls_pcls.i8 of edgeu_dumpls.C (pion-pair class of the simulated like-sign rows, same row order)
#   - the edgeu column files (edgeu_dump.C) for the f estimator per M(ee) band (data N, N++, N--; SMASH C2-C5)
# Output: eestep_ana.log, eestep_ana.json, eestep_figdata.npz in /hdd3/przygoda/HADES/FAT_out/all_report/eestep/.
# Estimators: the edgeu ones (edgeu_core.fit_step: pol2 + h T(u) on -45..+45 MeV/c2; like-sign rows with the Poisson model
# variance, N - CB with cov(CB) + model; h_sim = the same linear estimator applied to SMASH C2 + C3 with the C2 gamma*
# factors; peak term fixed at a_P x SMASH (C4 + C5)). S0, 0.30 <= M(pipi) < 0.40 (data scale), u = M4 - E(M(pipi), M(ee)).
# Run: python3 -I plots/all_report/eestep_ana.py   (numpy only, about 1 min)
# ============================================================================
import json, math, sys, time
import numpy as np
from scipy.stats import chi2 as chi2dist
sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import edgeu_core as C
import eestep_core as E

T0 = time.time()
OUT = E.OUT
LOG = open(OUT + "eestep_ana.log", "w")
def P(*a):
    s = " ".join(str(x) for x in a)
    print(s); LOG.write(s + "\n"); LOG.flush()
RES = {}
def hres(r):
    return dict(h=float(r["h"]), eh=float(r["eh"]), rows=float(r["rows"]), near=float(r["near"]), share=float(r["share"]),
                chi2=float(r["chi2"]), ndf=int(r["ndf"]))

# =====================================================================================================================
P("# eestep_ana.py, 2026-11-15: origin of the step at u = 0 in the like-sign e-e- rows")
L = E.load("ls"); uL = E.u_of(L)
U = E.load("usr"); uU = E.u_of(U)
baseL = (L["s0"] == 1) & E.in_slices(L["mpipi"])
baseU = (U["s0"] == 1) & E.in_slices(U["mpipi"])
MM = baseL & (L["k"] == 2); PPm = baseL & (L["k"] == 1)
P("rows: like-sign %d (pp %d, mm %d), e+e- edge region %d; S0, 0.30-0.40: pp %d, mm %d, e+e- %d" %
  (L["n"], (L["k"] == 1).sum(), (L["k"] == 2).sum(), U["n"], PPm.sum(), MM.sum(), baseU.sum()))
P("step template: Gaussian-smeared step, sigma = %.2f MeV/c2 (edgeu_ana.json)" % (1e3 * E.SIG_U))

# ---------------------------------------------------------------------------------------------------------------- 1
P("\n1. REPRODUCTION OF THE EDGEU NULL TESTS")
r_mm = E.step_rows(uL, MM); r_pp = E.step_rows(uL, PPm); r_us = E.step_rows(uU, baseU)
P("   e-e-: " + E.fmt(r_mm)); P("   e+e+: " + E.fmt(r_pp)); P("   e+e- (N, no CB): " + E.fmt(r_us))
RES["null"] = dict(mm=hres(r_mm), pp=hres(r_pp), us=hres(r_us))

# ---------------------------------------------------------------------------------------------------------------- 2
P("\n2. EVENT CONTENT OF THE e-e- ROWS (links to the e+e- rows of the same event; PAT event multiplicities)")
P("   lk_type: 0 = no e+e- PAT entry in the event; 2 = an e+e- entry with the same pi+, pi- and one of the two e- ('twin')")
splits = [
    ("no e+e- entry in the event (lk_type 0)", L["lk_type"] == 0),
    ("e+e- twin exists (lk_type 2)", L["lk_type"] == 2),
    ("  twin, odd e- shares the RICH ring of the twin e+", (L["lk_type"] == 2) & (L["lk_ring"] == 1)),
    ("  twin, no shared ring", (L["lk_type"] == 2) & (L["lk_ring"] == 0)),
    ("  twin, odd e- within 2 deg of the twin e+", (L["lk_type"] == 2) & (L["lk_ang"] < 2.0)),
    ("no positive lepton in the event (lpos_mult 0)", L["ev_lpos_mult"] == 0),
    ("positive lepton in the event (lpos_mult >= 1)", L["ev_lpos_mult"] >= 1),
    ("exactly 2 negative leptons (lneg_mult 2)", L["ev_lneg_mult"] == 2),
    ("3 or more negative leptons", L["ev_lneg_mult"] >= 3),
    ("isBest == 1", L["isBest"] == 1), ("isBest == 0", L["isBest"] == 0),
    ("the two e- share a RICH ring", L["shared_ring"] == 1), ("no shared ring", L["shared_ring"] == 0),
    ("one row in the event group (nrow 1)", L["nrow"] == 1), ("nrow >= 2", L["nrow"] >= 2),
    ("START hit (evtStart 1)", L["ev_evtStart"] == 1), ("no START", L["ev_evtStart"] != 1),
]
RES["content"] = {}
for lab, m in splits:
    r = E.step_rows(uL, MM & m)
    P("   %-55s %s" % (lab, E.fmt(r)))
    RES["content"][lab.strip()] = hres(r)
mreg = MM & (uL > -0.045) & (uL < 0.045)
P("   fractions of the e-e- rows in -45..+45 MeV/c2: lk_type 0 %.3f, lpos_mult 0 %.3f, isBest 1 %.3f, shared ring %.3f" %
  ((L["lk_type"][mreg] == 0).mean(), (L["ev_lpos_mult"][mreg] == 0).mean(), (L["isBest"][mreg] == 1).mean(), (L["shared_ring"][mreg] == 1).mean()))
P("   e+e- rows of the edge region with an e-e- twin (same pi+, pi-, e-; the e+ replaced by another negative candidate):")
for lab, m in [("no e-e- twin", U["lkmm_type"] != 2), ("e-e- twin", U["lkmm_type"] == 2), ("no e+e+ twin", U["lkpp_type"] != 2),
               ("e+e+ twin", U["lkpp_type"] == 2)]:
    r = E.step_rows(uU, baseU & m)
    P("   %-55s %s" % (lab, E.fmt(r)))
    RES["content"]["us " + lab] = hres(r)

# ---------------------------------------------------------------------------------------------------------------- 3
P("\n3. THE STEP PER M(ee) BAND: like-sign rows, e+e- rows, CB and N - CB against SMASH x f (edgeu estimator)")
D = C.load_data(); uD = C.u_data(D)
S = C.load_sim(); uS, mppS, _ = C.u_sim(S)
wS = C.sim_weight(S, 108, "S0", "A", c2corr=False)
s0D = (D["bits"] & 1) != 0
SL = C.slices_of(0.30, 0.40)
pbD = C.ppbin_of(D["mpipi"]); pbS = C.ppbin_of(mppS.astype(np.float32))
G = E.G; T = E.T
F_NOM, A_P = 0.4200, C.AP
def hk(k, mask):
    """per-slice histograms of data sample k"""
    out = []
    b = G.bin(uD)
    for (a, c) in SL:
        m = mask & s0D & (D["k"] == k) & (pbD >= a) & (pbD < c) & (b >= 0)
        out.append(np.bincount(b[m], minlength=G.nb).astype(np.float64))
    return out
def hsim(groups, mask, c2=True):
    b = G.bin(uS)
    m = mask & np.isin(S["g"], [C.G[g] for g in groups]) & (pbS >= SL[0][0]) & (pbS < SL[-1][1]) & (b >= 0) & (wS != 0)
    w = wS.copy()
    if c2:
        c2m = S["g"] == C.G["C2"]
        w = np.where(c2m & (S["mee"] < 0.010), w * C.CF_E1, w)
        w = np.where(c2m & ~(S["mee"] < 0.010), w * C.CF_E2, w)
    return np.bincount(b[m], weights=w[m], minlength=G.nb)
def band_fit(dmask, smask, cbmode=0):
    """edgeu nominal estimator for the rows in the masks; CB per slice from the like-sign rows of the same rows"""
    Nn = sum(hk(0, dmask)); pps = hk(1, dmask); mms = hk(2, dmask)
    CB = np.zeros(G.nb); COV = np.zeros((G.nb, G.nb))
    for p_, m_ in zip(pps, mms):
        cb, cov = C.cb_rule(p_, m_, G, 1.0, cbmode)
        CB += cb; COV += cov
    sp = C.Spec(Nn, CB, COV, sum(pps), sum(mms))
    K = hsim(["C4", "C5"], smask, c2=False)
    R = C.fit_step(sp, T, G, -0.045, 0.045, 2, fixed=A_P * K)
    ysim = hsim(["C2", "C3"], smask)
    hs = C.h_lin(R, ysim)
    hcb = C.h_lin(R, CB)
    hN = C.h_lin(R, Nn - A_P * K)   # step of N alone with the same linear estimator
    # like-sign step with their own Poisson fits
    rpp = E.step(sum(pps)); rmm = E.step(sum(mms))
    # separate errors of the N part and of the CB part of h (fixed V of the fit)
    Lh = R.Lh; var_model = np.maximum(R.model + R.fx + sp.cb[R.idx], 1.0)
    ehN = float(np.sqrt((Lh ** 2 * var_model).sum())); ehCB = float(np.sqrt(Lh @ COV[np.ix_(R.idx, R.idx)] @ Lh))
    return dict(hd=R.h, ehd=R.eh, ehN=ehN, ehCB=ehCB, hs=hs, f=R.h / hs if hs != 0 else np.nan, ef=R.eh / hs if hs != 0 else np.nan, hcb=hcb, hN=hN,
                chi2=R.chi2, ndf=R.ndf, hpp=rpp.h, ehpp=rpp.eh, hmm=rmm.h, ehmm=rmm.eh,
                npp=float(sum(pps)[rpp.idx].sum()), nmm=float(sum(mms)[rmm.idx].sum()), nN=float(Nn[R.idx].sum()), R=R, sp=sp, K=K, ysim=ysim)
allD = np.ones(D["n"], bool); allS = np.ones(S["n"], bool)
nomb = band_fit(allD, allS)
P("   all M(ee) (check against edgeu: h_data 1305.3, f 0.4200): h(N - CB) = %.1f +- %.1f, h_sim = %.1f, f = %.4f; h(N) = %.1f, "
  "h(CB) = %.1f; like-sign e+e+ %.1f +- %.1f, e-e- %.1f +- %.1f" % (nomb["hd"], nomb["ehd"], nomb["hs"], nomb["f"], nomb["hN"], nomb["hcb"],
                                                                   nomb["hpp"], nomb["ehpp"], nomb["hmm"], nomb["ehmm"]))
BANDS = [(0.0, 0.010), (0.010, 0.030), (0.030, 0.050), (0.050, 0.070), (0.070, 0.090), (0.090, 0.110), (0.110, C.MPI0)]
WIDE = [(0.0, 0.050), (0.050, 0.090), (0.090, C.MPI0), (0.070, C.MPI0)]
RES["bands"] = []
P("   %-13s %8s %8s %8s %8s | %8s %8s %8s | %8s %8s %8s %8s" % ("M(ee) MeV", "h(e+e+)", "+-", "h(e-e-)", "+-", "h(N)", "h(CB)", "h(N-CB)",
                                                           "+-", "f x hsim", "f_band", "rows e-e-"))
for (a, b) in BANDS + WIDE:
    md = (D["mee"] >= a) & (D["mee"] < b); ms = (S["mee"] >= a) & (S["mee"] < b)
    r = band_fit(md, ms)
    P("   %5.0f-%5.0f   %8.1f %8.1f %8.1f %8.1f | %8.1f %8.1f %8.1f | %8.1f %8.1f %8.3f %8.0f" %
      (1e3 * a, 1e3 * b, r["hpp"], r["ehpp"], r["hmm"], r["ehmm"], r["hN"], r["hcb"], r["hd"], r["ehd"], F_NOM * r["hs"], r["f"], r["nmm"]))
    RES["bands"].append(dict(lo=a, hi=b, **{k: float(v) for k, v in r.items() if k not in ("R", "sp", "K", "ysim")}))
# the test of the CB in the large-M(ee) band: N - CB against the SMASH bump x f
bb = [x for x in RES["bands"] if abs(x["lo"] - 0.070) < 1e-9 and abs(x["hi"] - C.MPI0) < 1e-9][0]
zz = (bb["hd"] - F_NOM * bb["hs"]) / bb["ehd"]
P("   M(ee) >= 70 MeV: h(N - CB) = %.1f +- %.1f against f x h_sim = %.1f (z = %.2f); h(N) = %.1f, h(CB) = %.1f; if the e-e- step were absent "
  "from the e+e- sample (artefact), N - CB would be about f x h_sim - h(CB) = %.1f" % (bb["hd"], bb["ehd"], F_NOM * bb["hs"], zz, bb["hN"], bb["hcb"], F_NOM * bb["hs"] - bb["hcb"]))
RES["cbtest70"] = dict(hd=bb["hd"], ehd=bb["ehd"], fhs=F_NOM * bb["hs"], z=zz, hN=bb["hN"], hcb=bb["hcb"])
# both readings in the two high bands, and the scale k_cross of the cross pairs in N relative to the CB rule:
#   N step = f h_sim + k_cross h(CB)  ->  k_cross = (h(N) - f h_sim) / h(CB)   (k_cross = 1: the CB rule; 0: like-sign artefact)
RES["readings"] = {}
for (a, b) in ((0.090, C.MPI0), (0.070, C.MPI0)):
    x = [y for y in RES["bands"] if abs(y["lo"] - a) < 1e-9 and abs(y["hi"] - b) < 1e-9][0]
    fh = F_NOM * x["hs"]
    z_phys = (x["hd"] - fh) / x["ehd"]
    z_art = (x["hd"] - (fh - x["hcb"])) / x["ehd"]
    k = (x["hN"] - fh) / x["hcb"]
    ek = math.sqrt(x["ehN"] ** 2 + (k * x["ehCB"]) ** 2) / x["hcb"]
    rat = x["hpp"] / x["hmm"]; erat = rat * math.sqrt((x["ehpp"] / max(abs(x["hpp"]), 1e-9)) ** 2 + (x["ehmm"] / x["hmm"]) ** 2)
    P("   M(ee) >= %.0f MeV: h(N) = %.1f +- %.1f, h(CB) = %.1f +- %.1f, h(N - CB) = %.1f +- %.1f, f x h_sim = %.1f;"
      " CB rule right (cross pairs also in e+e-): z = %+.2f; like-sign artefact (N - CB = f h_sim - h(CB) = %.1f): z = %+.2f;"
      " k_cross = %.2f +- %.2f; e+e+/e-e- step ratio %.2f +- %.2f" % (1e3 * a, x["hN"], x["ehN"], x["hcb"], x["ehCB"], x["hd"], x["ehd"], fh, z_phys,
                                                                    fh - x["hcb"], z_art, k, ek, rat, erat))
    RES["readings"]["%.0f" % (1e3 * a)] = dict(hN=x["hN"], ehN=x["ehN"], hcb=x["hcb"], ehCB=x["ehCB"], hd=x["hd"], ehd=x["ehd"], fhs=fh,
                                              z_phys=z_phys, z_art=z_art, k=k, ek=ek, ratio_pp_mm=rat, eratio=erat)
RES["kcross_f_shift"] = {}
for band in ("90", "70"):
    kk = RES["readings"][band]["k"]; ekk = RES["readings"][band]["ek"]
    P("   effect on f_eta(u-edge) of k_cross (M(ee) >= %s MeV value) instead of 1 (the CB step is %.1f rows = %.1f %% of h_data): %+.1f %% +- %.1f %%" %
      (band, nomb["hcb"], 100 * nomb["hcb"] / nomb["hd"], 100 * (1 - kk) * nomb["hcb"] / nomb["hd"], 100 * ekk * nomb["hcb"] / nomb["hd"]))
    RES["kcross_f_shift"][band] = dict(shift=(1 - kk) * nomb["hcb"] / nomb["hd"], err=ekk * nomb["hcb"] / nomb["hd"])
# opening-angle bands for the e-e- rows
P("   e-e- step per opening angle of the two leptons (deg):")
RES["oa"] = []
for (a, b) in [(0, 10), (10, 20), (20, 35), (35, 60), (60, 180)]:
    r = E.step_rows(uL, MM & (L["oa"] >= a) & (L["oa"] < b)); ru = E.step_rows(uU, baseU & (U["oa"] >= a) & (U["oa"] < b))
    P("      OA %3d-%3d: e-e- %s | e+e- N: h %.1f +- %.1f" % (a, b, E.fmt(r), ru["h"], ru["eh"]))
    RES["oa"].append(dict(lo=a, hi=b, mm=hres(r), us=hres(ru)))

# ---------------------------------------------------------------------------------------------------------------- 4
P("\n4. DOES THE e-e- STEP FOLLOW THE EDGE E(M(pipi), M(ee))? (e-e- rows, M(ee) >= 70 MeV)")
HI = MM & (L["mee"] >= 0.070)
def chi2_drop(var, mask, lo=-0.045, hi=0.045):
    y = E.hist_u(var, mask); R = E.step(y, lo, hi)
    return R.h, R.eh, R.chi20 - R.chi2, R.chi2, R.ndf
e_q = C.edge(L["mpipi"].astype(np.float64), L["mee"].astype(np.float64))
e_0 = C.edge(L["mpipi"].astype(np.float64), np.zeros(L["n"]))
var = {"u = M4 - E(m, q) (nominal)": L["m4"] - e_q, "M4 - E(m, 0) (no M(ee) term)": L["m4"] - e_0,
       "M4 - 0.535 (fixed M4)": L["m4"] - 0.535, "M4 - E(m, q) - 10 MeV": L["m4"] - e_q - 0.010, "M4 - E(m, q) + 10 MeV": L["m4"] - e_q + 0.010}
RES["align"] = {}
for k_, v in var.items():
    v = np.where(L["mee"] < C.MPI0, v, np.nan)
    h, eh, dchi, ch, nd = chi2_drop(v, HI)
    P("   %-32s step %.1f +- %.1f, delta chi2 (step vs none) %.1f, chi2 %.1f/%d" % (k_, h, eh, dchi, ch, nd))
    RES["align"][k_] = dict(h=h, eh=eh, dchi2=dchi)
# free position: scan u0
best = None; prof = []
for u0 in np.arange(-0.015, 0.0151, 0.0025):
    Ts = C.step_gauss(G, E.SIG_U, u0)
    y = E.hist_u(uL, HI)
    R = C.fit_step(C.Spec(y, np.zeros(G.nb), np.zeros((G.nb, G.nb))), Ts, G, -0.045, 0.045, 2)
    prof.append((float(u0), float(R.chi2), float(R.h)))
    if best is None or R.chi2 < best[1]:
        best = (u0, R.chi2, R.h)
P("   free edge position u0 (step at u0, same width): chi2 profile " + ", ".join("%+.1f: %.1f" % (1e3 * a, b) for a, b, _ in prof))
cmin = min(p[1] for p in prof)
inside = [p[0] for p in prof if p[1] <= cmin + 1.0]
P("   best u0 = %+.1f MeV/c2 (chi2 %.1f), 1-sigma range %+.1f .. %+.1f MeV/c2 (grid 2.5 MeV)" % (1e3 * best[0], best[1], 1e3 * min(inside), 1e3 * max(inside)))
RES["u0_profile"] = prof
# width scan
profw = []
for sg in [0.004, 0.0055, 0.0067, 0.008, 0.010, 0.013]:
    Ts = C.step_gauss(G, sg)
    y = E.hist_u(uL, HI)
    R = C.fit_step(C.Spec(y, np.zeros(G.nb), np.zeros((G.nb, G.nb))), Ts, G, -0.045, 0.045, 2)
    profw.append((sg, float(R.chi2), float(R.h)))
P("   step width scan (sigma MeV/c2: chi2, h): " + ", ".join("%.1f: %.1f, %.0f" % (1e3 * a, b, c) for a, b, c in profw))
RES["width_profile"] = profw
# per M(pipi) slice
P("   per M(pipi) slice (e-e-, M(ee) >= 70 MeV): the edge moves in M4 with M(pipi); the step stays at u = 0 if it is tied to E")
RES["slices"] = []
for (a, b) in [(0.30, 0.34), (0.34, 0.37), (0.37, 0.40)]:
    m = HI & E.in_slices(L["mpipi"], a, b)
    r = E.step_rows(uL, m)
    v0 = np.where(L["mee"] < C.MPI0, L["m4"] - 0.535, np.nan)
    h2, eh2, d2, _, _ = chi2_drop(v0, m)
    P("      %.2f-%.2f: in u %s | in M4 - 0.535: h %.1f +- %.1f" % (a, b, E.fmt(r), h2, eh2))
    RES["slices"].append(dict(lo=a, hi=b, u=hres(r), fixed=dict(h=h2, eh=eh2)))

# control region without eta kinematics (as the edgeu control test of N): M(pipi) 0.42-0.46, pseudo edge E(M(pipi) - 0.10, M(ee))
e_c = C.edge(L["mpipi"].astype(np.float64) - 0.10, L["mee"].astype(np.float64))
uc = np.where(L["mee"] < C.MPI0, L["m4"] - e_c, np.nan)
mc = (L["s0"] == 1) & (L["k"] == 2) & E.in_slices(L["mpipi"], 0.42, 0.46)
rc_all = E.step_rows(uc, mc); rc_hi = E.step_rows(uc, mc & (L["mee"] >= 0.070))
P("   control region M(pipi) 0.42-0.46 (no eta -> pi+pi-pi0 possible), pseudo edge E(M(pipi) - 0.10, M(ee)): e-e- all M(ee) %s" % E.fmt(rc_all))
P("                                                                                         e-e- M(ee) >= 70 MeV %s" % E.fmt(rc_hi))
RES["control"] = dict(all=hres(rc_all), hi=hres(rc_hi))

# ---------------------------------------------------------------------------------------------------------------- 5
P("\n5. TRACK AND ROW PROPERTIES OF THE STEP ROWS (e-e-, M(ee) >= 70 MeV; step share = h / rows per bin in -15..0 MeV/c2)")
P("   (a flipped Dalitz e+ would show up as a slot with bad track quality or a particular momentum / polar angle; the share is")
P("    flat in a variable when the step rows are distributed like the other e-e- rows)")
def scan(name, x, mask=HI, edges=None, q=(1 / 3, 2 / 3)):
    reg = mask & (uL > -0.045) & (uL < 0.045)
    if edges is None:
        edges = list(np.nanquantile(x[reg], q))
    e = [-np.inf] + list(edges) + [np.inf]
    parts = []; vals = []
    for a, b in zip(e[:-1], e[1:]):
        r = E.step_rows(uL, mask & (x >= a) & (x < b))
        parts.append("[%.3g,%.3g) %.2f+-%.2f" % (a, b, r["share"], r["eh"] / max(r["near"], 1)))
        vals.append((float(a), float(b), float(r["h"]), float(r["eh"]), float(r["share"]), float(r["eh"] / max(r["near"], 1))))
    sh = np.array([v[4] for v in vals]); es = np.array([v[5] for v in vals])
    w_ = 1 / es ** 2; mean = (sh * w_).sum() / w_.sum(); c2 = ((sh - mean) ** 2 * w_).sum()
    p = 1 - chi2dist.cdf(c2, len(sh) - 1)
    P("   %-26s %s   (flat: chi2 %.1f/%d, p %.2f)" % (name, " ; ".join(parts), c2, len(sh) - 1, p))
    RES["scan"][name] = dict(bins=vals, chi2=float(c2), ndf=len(sh) - 1, p=float(p))
RES["scan"] = {}
mx = np.maximum; mn = np.minimum
scan("p max (MeV/c)", mx(L["l1_pcor"], L["l2_pcor"])); scan("p min (MeV/c)", mn(L["l1_pcor"], L["l2_pcor"]))
scan("theta max (deg)", mx(L["l1_theta"], L["l2_theta"])); scan("theta min (deg)", mn(L["l1_theta"], L["l2_theta"]))
scan("RK chi2 max", mx(L["l1_rkchi2"], L["l2_rkchi2"])); scan("inner-seg chi2 max", mx(L["l1_mdcinnerchi2"], L["l2_mdcinnerchi2"]))
scan("outer-seg chi2 max", mx(L["l1_mdcouterchi2"], L["l2_mdcouterchi2"])); scan("META match max", mx(L["l1_metamatchquality"], L["l2_metamatchquality"]))
scan("RICH match max", mx(L["l1_richmatchquality"], L["l2_richmatchquality"])); scan("ring pads min", mn(L["l1_rich_padnum"], L["l2_rich_padnum"]))
scan("beta min", mn(L["l1_beta_new"], L["l2_beta_new"])); scan("systems (0 RPC,1 TOF) sum", L["l1_system"] + L["l2_system"], edges=[0.5, 1.5])
sec_same = (L["l1_sector"] == L["l2_sector"]).astype(np.float32)
scan("same sector (0/1)", sec_same, edges=[0.5])
scan("sector-2 lepton (0/1)", ((L["l1_sector"] == 2) | (L["l2_sector"] == 2)).astype(np.float32), edges=[0.5])
scan("opening angle (deg)", L["oa"]); scan("M(ee) (GeV/c2)", L["mee"])
scan("isBest", L["isBest"], edges=[0.5]); scan("day", L["day"], edges=[46.5, 52.5])
scan("pi- p (MeV/c)", L["pim_p"]); scan("pi+ p (MeV/c)", L["pip_p"]); scan("vertex z (mm)", L["vz"])
scan("hneg_mult", L["ev_hneg_mult"], edges=[1.5]); scan("hpos_mult", L["ev_hpos_mult"], edges=[1.5, 2.5])
# sector of each lepton: 6 sectors, share per sector (any lepton in sector s)
P("   sector content (rows with at least one e- in sector s), e-e- M(ee) >= 70 MeV:")
RES["sector"] = []
for s in range(6):
    r = E.step_rows(uL, HI & ((L["l1_sector"] == s) | (L["l2_sector"] == s)))
    P("      sector %d: %s" % (s, E.fmt(r)))
    RES["sector"].append(hres(r))

# ---------------------------------------------------------------------------------------------------------------- 6
P("\n6. THE CHARGE-FLIP READING: LIMIT FROM THE M(ee) BAND OF THE BUMP")
P("   A flipped Dalitz / conversion e+ keeps its momentum, so the e-e- row has the same M(ee) and u as the e+e- bump row.")
P("   92 % of the bump step in N is at M(ee) < 50 MeV; the flipped rows would be there too.")
lo50 = [x for x in RES["bands"] if x["lo"] == 0.0 and abs(x["hi"] - 0.050) < 1e-9][0]
rflip = lo50["hmm"] / lo50["hN"]; erflip = lo50["ehmm"] / lo50["hN"]
ul90 = (lo50["hmm"] + 1.2816 * lo50["ehmm"]) / lo50["hN"]
P("   M(ee) < 50 MeV: e-e- step %.1f +- %.1f against the e+e- step %.1f (N) -> flip rate r = %.4f +- %.4f, r < %.4f (90 %% CL, Gaussian)" %
  (lo50["hmm"], lo50["ehmm"], lo50["hN"], rflip, erflip, ul90))
# effect on f if every flipped row were a lost bump row and the CB took it in: migration (r) + CB part (sqrt(N++/N--) r h)
ppn = lo50["npp"]; mmn = lo50["nmm"]
cbpart = math.sqrt(ppn / mmn) * ul90 * lo50["hN"] / nomb["hd"]
P("   at the limit: lost bump rows %.1f %% of h(N), CB over-subtraction sqrt(N++/N--) x flipped rows = %.2f %% of h_data -> f higher by < %.1f %% (90 %% CL)" %
  (100 * ul90, 100 * cbpart, 100 * (ul90 + cbpart)))
RES["flip"] = dict(r=rflip, er=erflip, ul90=ul90, cbpart=cbpart, fbias_ul=ul90 + cbpart, h_mm_lo50=lo50["hmm"], eh_mm_lo50=lo50["ehmm"],
                   h_N_lo50=lo50["hN"])
twin = RES["content"]["us e-e- twin"]
P("   e+e- rows of the edge region that have an e-e- twin (wrong-sign copy of the e+): step %.1f +- %.1f (bump step in N %.1f): clone rate %.4f +- %.4f" %
  (twin["h"], twin["eh"], r_us["h"], twin["h"] / r_us["h"], twin["eh"] / r_us["h"]))
RES["clone"] = dict(h=twin["h"], eh=twin["eh"], rate=twin["h"] / r_us["h"], erate=twin["eh"] / r_us["h"])

# ---------------------------------------------------------------------------------------------------------------- 7
P("\n7. THE PEAK WINDOW AND THE X17 WINDOW")
PLO, PHI = 0.53364, 0.55416
# kinematic reach of the cross-photon pairs: they obey M4 <= E(M(pipi), M(ee)), so in P only if E >= 0.53364
qgrid = np.arange(0.0, C.MPI0, 0.0005)
for m_ in (0.28, 0.30, 0.35, 0.40):
    e_ = C.edge(np.full_like(qgrid, m_), qgrid)
    qm = qgrid[np.nanargmax(e_ >= PLO)] if np.any(e_ >= PLO) else np.nan
    P("   M(pipi) = %.2f: E(m, q) reaches the peak window (0.53364) for M(ee) >= %.1f MeV/c2; E(m, 24.5 MeV) = %.4f" % (m_, 1e3 * qm, C.edge(np.array([m_]), np.array([0.0245]))[0]))
emax_x17 = float(np.nanmax(C.edge(np.arange(0.27, 0.41, 0.001), np.full(140, 0.0245))))
P("   largest E for M(ee) <= 24.5 MeV/c2 (X17 band): %.4f GeV/c2 = %.1f MeV/c2 below the peak window (%.1f sigma of the u resolution)" %
  (emax_x17, 1e3 * (PLO - emax_x17), (PLO - emax_x17) / E.SIG_U))
# Inside a window the u distribution is cut by the window edges, so the step cannot be fitted there. Instead: the e-e- rows
# of the window that sit at or below their own edge (u <= 0; a cross-photon pair can only be there), times the cross-pair
# share at the edge measured in M(ee) >= 90 MeV (step / rows per bin), give the cross pairs of the window; the CB of the
# window changes by sqrt(N++/N--) per e-e- row (2 sqrt(N++ N--)); Y_peak = (N - CB)_P - R'(N - CB)_SB (WP6, R' about 1).
r90 = E.step_rows(uL, MM & (L["mee"] >= 0.090))
share90, eshare90 = r90["share"], r90["eh"] / r90["near"]
P("   cross-pair share of the e-e- rows just below the edge (M(ee) >= 90 MeV, 0.30-0.40): %.2f +- %.2f" % (share90, eshare90))
RES["peak"] = dict(share90=share90, eshare90=eshare90, emax_x17=emax_x17)
lsS0 = (L["s0"] == 1) & np.isfinite(uL)
for wn, a, b in (("SBlo", 0.5225, PLO), ("P", PLO, PHI), ("SBhi", PHI, 0.565)):
    inw = lsS0 & (L["m4"] >= a) & (L["m4"] < b)
    npp_w = ((L["k"] == 1) & (L["s0"] == 1) & (L["m4"] >= a) & (L["m4"] < b)).sum()
    nmm_w = ((L["k"] == 2) & (L["s0"] == 1) & (L["m4"] >= a) & (L["m4"] < b)).sum()
    bel = (inw & (L["k"] == 2) & (uL <= 0)).sum()
    belhi = (inw & (L["k"] == 2) & (uL <= 0) & (L["mee"] >= 0.090)).sum()
    xmm = share90 * bel
    dcb = math.sqrt(npp_w / nmm_w) * xmm * (1 + 0.13 * nmm_w / npp_w)   # e+e+ partners at 0.13 of the e-e- cross pairs
    P("   %-4s M4 %.5f-%.5f: N++ %d, N-- %d; e-e- rows with u <= 0: %d (%d with M(ee) >= 90 MeV) -> cross pairs about %.0f; "
      "CB carries about %.0f rows of them (%.2f %% of Y_peak 5088 x 63/64)" % (wn, a, b, npp_w, nmm_w, bel, belhi, xmm, dcb, 100 * dcb / (5088 * 63 / 64)))
    RES["peak"][wn] = dict(npp=int(npp_w), nmm=int(nmm_w), below=int(bel), below_hi=int(belhi), cross=float(xmm), dcb=float(dcb))
dP = RES["peak"]["P"]["dcb"]; dS = RES["peak"]["SBlo"]["dcb"] + RES["peak"]["SBhi"]["dcb"]
P("   if these like-sign rows had NO unlike-sign partners (artefact reading): the CB of P too high -> Y_peak lower by about %.0f rows (%.1f %%) "
  "from P alone; with the side bands at R' = 1 the net change is %+.0f rows (%+.1f %%, Y_peak higher); with partners (cross-photon pairs, "
  "the CB rule) the effect is zero" %
  (dP * 64 / 63, 100 * dP / (5088 * 63 / 64), (dS - dP) * 64 / 63, 100 * (dS - dP) / (5088 * 63 / 64)))
RES["peak"]["dY_P_only"] = float(-dP * 64 / 63); RES["peak"]["dY_with_SB"] = float((dS - dP) * 64 / 63)
x17 = (L["k"] == 2) & (L["s0"] == 1) & (L["m4"] >= PLO) & (L["m4"] < PHI) & (L["mee"] >= 0.0155) & (L["mee"] < 0.0245)
x17pp = (L["k"] == 1) & (L["s0"] == 1) & (L["m4"] >= PLO) & (L["m4"] < PHI) & (L["mee"] >= 0.0155) & (L["mee"] < 0.0245)
P("   X17P window (P x M(ee) 15.5-24.5 MeV, S0, M(pipi) 0.26-0.46 dumped): e-e- rows %d, e+e+ rows %d (D19 full window: 223 / 14); "
  "none of them can be a cross-photon pair below its edge (E <= %.4f)" % (x17.sum(), x17pp.sum(), emax_x17))
RES["x17"] = dict(mm=int(x17.sum()), pp=int(x17pp.sum()))

# ---------------------------------------------------------------------------------------------------------------- 8
P("\n8. SIMULATION (mechanism only; the simulation has no combinatorial background, no k or CB number is taken from it)")
Us = E.load("sus", "eestep_simlayout.txt"); Ls = E.load("sls", "eestep_simlayout.txt")
pcls = np.fromfile(C.OUT + "simls_pcls.i8", dtype=np.int8)
assert pcls.size == Ls["n"] and np.array_equal(np.fromfile(C.OUT + "simls_m4.f32", dtype=np.float32), Ls["m4"])
def usim(Dd):
    Ss = dict(m4=Dd["m4"], m4t=Dd["m4t"], mpipi=Dd["mpipi"], mpipit=Dd["mpipit"], mee=Dd["mee"])
    u_, mp_, _ = C.u_sim(Ss); return u_, mp_
uUs, mUs = usim(Us); uLs, mLs = usim(Ls)
wUs = Us["w"] * Us["f108"]; wLs = Ls["w"] * Ls["f108"]
bU = Us["bits"].astype(int); bL = Ls["bits"].astype(int)
sU = ((bU & 2) != 0) & ((bU & 16) != 0) & (mUs >= 0.30) & (mUs < 0.40) & np.isfinite(uUs)
sL = ((bL & 2) != 0) & ((bL & 16) != 0) & (mLs >= 0.30) & (mLs < 0.40) & np.isfinite(uLs)
RES["sim"] = {}
P("   (a) wrong-charge true leptons: like-sign rows whose two slots hold a true e+ and a true e- of the same source, per e+e- true pair")
for st, nm in ((0, "std"), (1, "lep"), (2, "eta")):
    tp = sU & (Us["set"] == st) & np.isin(Us["lcls"], [1, 2, 3, 4, 5, 6, 10])
    for k, kn in ((2, "e-e-"), (1, "e+e+")):
        lk = sL & (Ls["set"] == st) & (Ls["k"] == k)
        wrong = lk & (((Ls["l1_id"] == 2) & (Ls["l2_id"] == 3)) | ((Ls["l1_id"] == 3) & (Ls["l2_id"] == 2)))
        same = wrong & (((Ls["l1_par"] == Ls["l2_par"]) & np.isin(Ls["l1_par"], [1, 7, 17])) if st == 0 else ((Ls["l1_par"] == -1) & (Ls["l2_par"] == -1)))
        rate_n = same.sum() / max(tp.sum(), 1); rate_w = wLs[same].sum() / wUs[tp].sum()
        P("      %s %s: true pairs %d rows; same-source wrong-charge rows %d (weighted ratio %.5f, row ratio %.5f; with an e+e- clone %d)" %
          (nm, kn, tp.sum(), same.sum(), rate_w, rate_n, (same & (Ls["lk_type"] == 2)).sum()))
        RES["sim"]["flip_%s_%s" % (nm, "mm" if k == 2 else "pp")] = dict(true_pairs=int(tp.sum()), rows=int(same.sum()), rate_w=float(rate_w), rate_n=float(rate_n))
# momentum dependence (lep set, largest statistics): wrong-charge e+ in an e- slot vs e+ of true pairs
tp = sU & (Us["set"] == 1) & np.isin(Us["lcls"], [1, 2, 3, 4, 5, 6, 10]) & (Us["l1_id"] == 2)
lk = sL & (Ls["set"] == 1) & (Ls["k"] == 2)
pe = np.concatenate([Ls["l1_p"][lk & (Ls["l1_id"] == 2)], Ls["l2_p"][lk & (Ls["l2_id"] == 2)]])
pb_ = [0, 100, 200, 300, 500, 800, 5000]
hw, _ = np.histogram(pe, pb_); hd_, _ = np.histogram(Us["l1_p"][tp], pb_)
P("   (b) lep set, true e+ reconstructed as negative (any partner) per true e+ of e+e- true pairs, by reconstructed momentum [MeV/c]:")
P("      " + "; ".join("%d-%d: %d / %d = %.4f" % (pb_[i], pb_[i + 1], hw[i], hd_[i], hw[i] / max(hd_[i], 1)) for i in range(len(hw))))
RES["sim"]["flip_vs_p"] = dict(edges=pb_, wrong=hw.tolist(), ref=hd_.tolist())
# (c) cross-photon like-sign pairs of eta -> pi+pi-pi0 (std set): pions from the same eta (pcls 2), both leptons from pi0 photons
def anc(i):
    par, gpar = Ls["l%d_par" % i], Ls["l%d_gpar" % i]
    return (par == 7) | ((par == 1) & (gpar == 7))
P("   (c) std set: like-sign rows with the pions of one eta (pcls 2) and both leptons from pi0 photons (Dalitz or conversion):")
for k, kn in ((2, "e-e-"), (1, "e+e+")):
    lk = sL & (Ls["set"] == 0) & (Ls["k"] == k)
    cr = lk & (pcls == 2) & anc(1) & anc(2)
    near = cr & (uLs > -0.045) & (uLs < 0.045)
    y = E.hist_u(uLs, cr, wLs); y2 = E.hist_u(uLs, cr, wLs * wLs)
    R = E.step(y); eh = float(np.sqrt((R.Lh ** 2 * y2[R.idx]).sum()))
    mb = cr & (uLs > -0.015) & (uLs <= 0); ma = cr & (uLs > 0) & (uLs <= 0.015)
    bel = wLs[mb].sum() / 6; abo = wLs[ma].sum() / 6
    ebel = math.sqrt((wLs[mb] ** 2).sum()) / 6; eabo = math.sqrt((wLs[ma] ** 2).sum()) / 6
    P("      %s near the edge: %d rows in -15 < u <= 0 (%.1f +- %.1f per 2.5 MeV/c2), %d rows in 0 < u <= 15 (%.1f +- %.1f)" %
      (kn, mb.sum(), bel, ebel, ma.sum(), abo, eabo))
    P("      %s: %d rows of %d (weighted %.0f of %.0f expected raw rows in 0.30-0.40); in -45..45 MeV/c2 %d rows, median M(ee) %.0f MeV/c2, "
      "median OA %.0f deg; step %.0f +- %.0f (weighted); per 2.5 MeV/c2 below the edge %.1f, above %.1f" %
      (kn, cr.sum(), lk.sum(), wLs[cr].sum(), wLs[lk].sum(), near.sum(), 1e3 * np.median(Ls["mee"][near]) if near.sum() else np.nan,
       np.median(Ls["oa"][near]) if near.sum() else np.nan, R.h, eh, bel, abo))
    RES["sim"]["cross_%s" % ("mm" if k == 2 else "pp")] = dict(rows=int(cr.sum()), w=float(wLs[cr].sum()), w_all=float(wLs[lk].sum()), near=int(near.sum()),
        mee_med=float(np.median(Ls["mee"][near])) if near.sum() else None, oa_med=float(np.median(Ls["oa"][near])) if near.sum() else None,
        h=float(R.h), eh=eh, below=float(bel), above=float(abo), ebelow=float(ebel), eabove=float(eabo), nbelow=int(mb.sum()), nabove=int(ma.sum()))
cm = RES["sim"]["cross_mm"]; cp = RES["sim"]["cross_pp"]
P("      weighted e+e+/e-e- ratio of these pairs: %.3f (data like-sign step ratio in M(ee) >= 90 MeV: see Sect. 3)" % (cp["w"] / cm["w"]))
hs_nom = nomb["hs"]
P("      SMASH x f (f = %.3f) density of e-e- cross pairs just below the edge: %.1f rows per 2.5 MeV/c2, i.e. %.3f of h_sim x f "
  "(data: e-e- step %.1f +- %.1f = %.3f +- %.3f of h_data)" % (F_NOM, F_NOM * cm["below"], cm["below"] / hs_nom, r_mm["h"], r_mm["eh"],
                                                               r_mm["h"] / nomb["hd"], r_mm["eh"] / nomb["hd"]))
RES["sim"]["cross_below_fx"] = F_NOM * cm["below"]
d_ = cm["below"] - cm["above"]; ed_ = math.hypot(cm["ebelow"], cm["eabove"])
P("      SMASH x f: at most %.1f +- %.1f rows per 2.5 MeV/c2 (density below the edge), %.1f +- %.1f as below minus above; data e-e- step %.1f +- %.1f"
  " -> data / (SMASH x f) >= %.1f (upper SMASH value), z = %.1f" % (F_NOM * cm["below"], F_NOM * cm["ebelow"], F_NOM * d_, F_NOM * ed_, r_mm["h"], r_mm["eh"],
                                                                   r_mm["h"] / (F_NOM * cm["below"]), (r_mm["h"] - F_NOM * cm["below"]) / math.hypot(r_mm["eh"], F_NOM * cm["ebelow"])))
RES["sim"]["cross_cmp"] = dict(upper=F_NOM * cm["below"], eupper=F_NOM * cm["ebelow"], diff=F_NOM * d_, ediff=F_NOM * ed_,
                               ratio_min=r_mm["h"] / (F_NOM * cm["below"]), z=(r_mm["h"] - F_NOM * cm["below"]) / math.hypot(r_mm["eh"], F_NOM * cm["ebelow"]))

# ---------------------------------------------------------------------------------------------------------------- 9
P("\n9. FIGURE DATA")
fd = dict(x=G.x)
for lab, m in (("mm_all", MM), ("pp_all", PPm), ("mm_lo", MM & (L["mee"] < 0.050)), ("mm_mid", MM & (L["mee"] >= 0.050) & (L["mee"] < 0.090)),
               ("mm_hi", MM & (L["mee"] >= 0.090)), ("pp_hi", PPm & (L["mee"] >= 0.090)), ("pp_lo", PPm & (L["mee"] < 0.090))):
    r = E.step_rows(uL, m)
    fd["y_" + lab] = r["y"]; fd["m_" + lab] = r["R"].model; fd["s_" + lab] = r["R"].smooth; fd["i_" + lab] = r["R"].idx
fd["usN_lo"] = E.step_rows(uU, baseU & (U["mee"] < 0.050))["y"]; fd["usN_hi"] = E.step_rows(uU, baseU & (U["mee"] >= 0.090))["y"]
fd["oa_mm"] = L["oa"][mreg]; fd["u_mm"] = uL[mreg]; fd["mee_mm"] = L["mee"][mreg]
ureg = baseU & (uU > -0.045) & (uU < 0.045)
fd["oa_us"] = U["oa"][ureg]; fd["u_us"] = uU[ureg]; fd["mee_us"] = U["mee"][ureg]
crm = sL & (Ls["set"] == 0) & (Ls["k"] == 2) & (pcls == 2) & anc(1) & anc(2)
fd["sim_cross_u"] = uLs[crm]; fd["sim_cross_w"] = wLs[crm]; fd["sim_cross_mee"] = Ls["mee"][crm]
np.savez(OUT + "eestep_figdata.npz", **fd)
json.dump(RES, open(OUT + "eestep_ana.json", "w"), indent=1, default=float)
P("written eestep_ana.json, eestep_figdata.npz; %.0f s" % (time.time() - T0))
