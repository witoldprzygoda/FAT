#!/usr/bin/env python3
# ============================================================================
# c2recon_edgeu.py  --  package c2recon: the C2 gamma* factor that the u-edge estimator (edgeu) actually needs
# ============================================================================
# Written 2026-11-15 (prefix c2recon). New file; edgeu_core.py is imported and only read (no edgeu output is written).
#
# Steps:
#   1. rebuild the nominal edgeu fit with the edgeu code (edgeu_core: u, CB rule, step template, iterated GLS;
#      the configuration of edgeu_ana.py NOM) and check it against edgeu_ana.json (h_data, h_sim, C2 raw / corrected, C3);
#   2. join the simulated C2 rows of the edgeu columns with the per-row QED weights of c2recon_weights.py
#      (key: tree, event group, M4, true M4);
#   3. pass the REWEIGHTED C2 template through the same linear estimator (final data covariance) -> h(C2 QED); the
#      factor the estimator needs is c_u = h(C2 QED) / h(C2 isotropic), per M(ee) class and for all;
#   4. f_eta(u-edge) with the reweighted template, and with the applied class factors (0.879 / 0.8825) for comparison;
#   5. the statistical error of c_u (100 WP0 replicas of the simulation, plus 2000 row-level Poisson replicas),
#      the dependence on the density variant of the weights, and on the estimator variants (range, order, smearing,
#      E(m, 0), slices);
#   6. other definitions on the same rows: all M4, window B, -15...+5 MeV, far below, per softer-lepton momentum.
#
# Output: /hdd3/przygoda/HADES/FAT_out/all_report/c2recon/c2recon_edgeu.json, c2recon_edgeu.log, c2recon_rowcoef_edgeu.npz
# Run: python3 -I plots/all_report/c2recon_edgeu.py   (about 1 min, 3 GB)
# ============================================================================
import json
import math
import sys
import time

import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import edgeu_core as C

O = "/hdd3/przygoda/HADES/FAT_out/all_report/c2recon/"
T0 = time.time()
LOG = open(O + "c2recon_edgeu.log", "w")


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s); LOG.write(s + "\n"); LOG.flush()


RES = {}
D = C.load_data()
S = C.load_sim()
gS = S["g"]
P("# c2recon_edgeu.py, 2026-11-15: data rows %d, simulated rows %d (%.0f s)" % (D["n"], S["n"], time.time() - T0))

# ---- the simulated u resolution at the edge (copy of edgeu_ana.py section 0)
w_nom_raw = C.sim_weight(S, c2corr=False)
u_s0, mpp_s0, m4p_s0 = C.u_sim(S)
ut_s = C.EDGE_SCALE * (S["m4t"].astype(np.float64) - C.m4max_true(S["mpipit"].astype(np.float64), S["meet"].astype(np.float64)))


def reso_sample(u, mpp, w, near=-0.015):
    m = (gS <= 1) & np.isfinite(u) & np.isfinite(ut_s) & (w > 0) & (mpp >= 0.30) & (mpp < 0.40) & (ut_s > near)
    return (u - ut_s)[m], w[m]


def gauss_core(d, ww):
    h, e = np.histogram(d, bins=np.arange(-0.040, 0.0401, 0.002), weights=ww)
    h2, _ = np.histogram(d, bins=e, weights=ww * ww)
    x = 0.5 * (e[1:] + e[:-1])
    f = lambda x, a, mu, s: a * np.exp(-0.5 * (x - mu) ** 2 / s ** 2)
    p, cv = curve_fit(f, x, h, p0=[h.max(), 0, 0.007], sigma=np.sqrt(np.maximum(h2, 1.0)), absolute_sigma=True)
    return float(p[1]), float(abs(p[2])), float(math.sqrt(cv[2, 2]))


dS, wdS = reso_sample(u_s0, mpp_s0, w_nom_raw)
_, SIG_U, _ = gauss_core(dS, wdS)
P("   step template width (edgeu_ana.py rule): %.4f MeV (edgeu: 6.70)" % (1e3 * SIG_U))

# ---- the nominal estimator (edgeu_ana.py NOM)
SLEDGES = (0.30, 0.34, 0.37, 0.39, 0.40)
CF = np.array([C.CF_E1, C.CF_E2, C.CF_E2])


def build(lo=-0.045, hi=0.045, order=2, sigma=None, useq=True, sledges=SLEDGES, s=C.S_NOM, r=C.R_NOM, sigma_sim=None, cls="all"):
    g = C.Grid(-0.200, 0.0025, 140, 4)
    sig = SIG_U if sigma is None else sigma
    T = C.step_gauss(g, sig)
    ud = C.u_data(D, useq=useq)
    m = (D["bits"] & 1) != 0
    H = C.hist_data(D, ud, m, g)
    sl = C.slices_of(sledges[0], sledges[-1], list(sledges))
    sp = C.data_spec(H, sl, cls, g)
    us, mpp, _ = C.u_sim(S, s=s, r=r, useq=useq)
    Hs = C.hist_sim(S, us, mpp, w_nom_raw, g)
    def simpart(groups, c2=True, cls=tuple(C.CLS_SETS[cls])):
        y = np.zeros(g.nb)
        for gname in groups:
            gi = C.G[gname]
            for (a, b) in sl:
                for k in cls:
                    y = y + Hs[gi, a:b, k].sum(0) * (CF[k] if (gname == "C2" and c2) else 1.0)
        return y
    peakT = simpart(["C4", "C5"])
    R = C.fit_step(sp, T, g, lo, hi, order, fixed=C.AP * peakT)
    if sigma_sim is not None:   # own smearing for the simulation: a separate fit (not linear in the same way)
        Ts = C.step_gauss(g, sigma_sim)
    else:
        Ts = None
    return dict(g=g, sl=sl, sp=sp, R=R, Hs=Hs, simpart=simpart, us=us, mpp=mpp, Ts=Ts, lo=lo, hi=hi, order=order, cls=cls)


def row_coef(B):
    """linear coefficient of every simulated row in h_sim (rows outside the estimator: 0)"""
    g, R = B["g"], B["R"]
    b = g.bin(B["us"])
    pb = C.ppbin_of(B["mpp"].astype(np.float32))
    inslice = np.zeros(S["n"], bool)
    for (a, bb) in B["sl"]:
        inslice |= (pb >= a) & (pb < bb)
    L = np.zeros(g.nb); L[R.idx] = R.Lh
    incls = np.isin(C.cls_of(S["mee"]), C.CLS_SETS[B["cls"]])
    coef = np.where((b >= 0) & inslice & incls & np.isfinite(B["us"]), L[np.clip(b, 0, g.nb - 1)], 0.0)
    return coef


B0 = build()
R0 = B0["R"]
lin = {k: C.h_lin(R0, B0["simpart"]([k], c2=False)) for k in ("C2", "C3")}
lin["C2c"] = C.h_lin(R0, B0["simpart"](["C2"], c2=True))
hs0 = lin["C2c"] + lin["C3"]
P("   nominal: h_data = %.2f +- %.2f, chi2 %.2f/%d; h_sim = %.2f (C2 raw %.2f, C2 x 0.879/0.8825 %.2f, C3 %.2f); f = %.5f, uncorrected %.5f" %
  (R0.h, R0.eh, R0.chi2, R0.ndf, hs0, lin["C2"], lin["C2c"], lin["C3"], R0.h / hs0, R0.h / (lin["C2"] + lin["C3"])))
ref = json.load(open(C.OUT + "edgeu_ana.json"))["nom"]
P("   edgeu_ana.json: h_data %.2f, h_sim %.2f, C2 %.2f, C2c %.2f, C3 %.2f, f %.5f" % (ref["hd"], ref["hs"], ref["lin"]["C2"], ref["lin"]["C2c"], ref["lin"]["C3"], ref["f"]))
RES["nominal_check"] = dict(hd=R0.h, ehd=R0.eh, hs=hs0, C2=lin["C2"], C2c=lin["C2c"], C3=lin["C3"], f=R0.h / hs0,
                            ref=dict(hd=ref["hd"], hs=ref["hs"], C2=ref["lin"]["C2"], C2c=ref["lin"]["C2c"], C3=ref["lin"]["C3"], f=ref["f"]))

# ---- join the per-row weights
Wt = np.load(O + "c2recon_weights.npz")
key_w = {}
for i in range(len(Wt["tree"])):
    key_w[(int(Wt["tree"][i]), int(Wt["grp"][i]), float(np.float32(Wt["m4"][i])), float(np.float32(Wt["m4_true"][i])))] = i
c2rows = np.where(gS == 0)[0]
jw = np.full(S["n"], -1, np.int64)
for i in c2rows:
    jw[i] = key_w.get((int(S["tree"][i]), int(S["grp"][i]), float(S["m4"][i]), float(S["m4t"][i])), -1)
nC2used = int(((w_nom_raw > 0) & (gS == 0)).sum())
nJ = int(((jw >= 0) & (w_nom_raw > 0) & (gS == 0)).sum())
P("   C2 rows in the edgeu columns: %d (S0, mode A, w > 0: %d); joined with a weight: %d" % (len(c2rows), nC2used, nJ))
assert nJ == nC2used
QV = ("q_A", "q_B", "q_Ac", "q_flat", "q_A_kloe")
Q = {k: np.ones(S["n"]) for k in QV}
for k in QV:
    Q[k][jw >= 0] = Wt[k][jw[jw >= 0]]
isC2 = (gS == 0) & (w_nom_raw > 0)
w = w_nom_raw
mee = S["mee"]
E1 = mee < 0.010
P("   check: sum W over the joined C2 rows %.3f (c2recon_weights S0: %.3f)" % (w[isC2].sum(), (Wt["W108"] * (Wt["s0"] > 0)).sum()))


def factor(coef, q, m):
    num = (coef * w * q)[m].sum(); den = (coef * w)[m].sum()
    return num / den, num, den


coef0 = row_coef(B0)
# the reproduction through the row coefficients
P("   row coefficients: sum a W (C2) = %.3f (linear h of the C2 template %.3f)" % ((coef0 * w)[isC2].sum(), lin["C2"]))
out = {}
for cl, mc in (("all", isC2), ("E1", isC2 & E1), ("E2", isC2 & ~E1)):
    out[cl] = {}
    for k in QV:
        c, num, den = factor(coef0, Q[k], mc)
        out[cl][k] = dict(c=c, hq=num, hiso=den)
P("   FACTOR c_u = h(C2, QED weights) / h(C2, isotropic) of the u-edge estimator (density A = nominal):")
for cl in ("all", "E1", "E2"):
    P("      %-3s: %.4f  (B %.4f, coarse %.4f, flat %.4f, KLOE-X %.4f); h_iso %.2f -> h_QED %.2f" %
      (cl, out[cl]["q_A"]["c"], out[cl]["q_B"]["c"], out[cl]["q_Ac"]["c"], out[cl]["q_flat"]["c"], out[cl]["q_A_kloe"]["c"],
       out[cl]["q_A"]["hiso"], out[cl]["q_A"]["hq"]))
sh2 = out["E2"]["q_A"]["hiso"] / out["all"]["q_A"]["hiso"]
P("   E2 share of the isotropic C2 step: %.3f (verifier: 0.656); applied class factors give the effective %.4f" %
  (sh2, (C.CF_E1 * out["E1"]["q_A"]["hiso"] + C.CF_E2 * out["E2"]["q_A"]["hiso"]) / out["all"]["q_A"]["hiso"]))
RES["factor"] = out
RES["E2share"] = sh2
# f with the reweighted template
h2q = out["all"]["q_A"]["hq"]
f_new = R0.h / (h2q + lin["C3"])
f_app = R0.h / hs0
f_unc = R0.h / (lin["C2"] + lin["C3"])
P("   f_eta(u-edge): applied factors %.5f; event-by-event QED weights %.5f (%+.2f %%); uncorrected %.5f" % (f_app, f_new, 100 * (f_new / f_app - 1), f_unc))
P("   C2 share of the step (QED-weighted): %.4f" % (h2q / (h2q + lin["C3"])))
RES["f"] = dict(applied=f_app, qed=f_new, unc=f_unc, rel=f_new / f_app - 1, h2q=h2q, C3=lin["C3"], hd=R0.h)

# ---- statistical error of the factor
# (a) the WP0 replicas (Poisson(1) per event group; the edgeu MC replicas)
rows = np.where(isC2)[0]
BW = C.sim_replica_weights(S, rows, 100)
num = (BW * (coef0 * w * Q["q_A"])[rows]).sum(1); den = (BW * (coef0 * w)[rows]).sum(1)
cr = num / den
fr = R0.h / (num + lin["C3"])     # C3 fixed: the factor part of the MC error only
P("   factor MC error (100 WP0 replicas of the C2 rows): c_u mean %.4f, rms %.4f" % (cr.mean(), C.sd(cr)))
# (b) row-level Poisson replicas (2000, independent seed), per class
rng = np.random.default_rng(0xC2EC0DE)
keyg = S["tree"][rows].astype(np.int64) * (1 << 32) + S["grp"][rows].astype(np.int64)
uq, inv = np.unique(keyg, return_inverse=True)
K = rng.poisson(1.0, size=(2000, uq.size))[:, inv]
mE1 = E1[rows]
stat = {}
for cl, mm in (("all", np.ones(len(rows), bool)), ("E1", mE1), ("E2", ~mE1)):
    a = (coef0 * w * Q["q_A"])[rows] * mm; b = (coef0 * w)[rows] * mm
    cc = (K * a).sum(1) / (K * b).sum(1)
    stat[cl] = dict(rms=float(C.sd(cc)), h68=float(C.h68(cc)), mean=float(cc.mean()))
    P("   factor stat (2000 group-Poisson replicas) %-3s: rms %.4f, half 68 %% %.4f, mean %.4f" % (cl, stat[cl]["rms"], stat[cl]["h68"], stat[cl]["mean"]))
RES["stat_wp0"] = dict(mean=float(cr.mean()), rms=float(C.sd(cr)))
RES["stat"] = stat
# what part of the factor's error is in the applied-constant approach anyway: the MC of h(C2) itself
hiso_r = den
P("   for comparison: the replica rms of h(C2 iso) is %.1f rows (%.1f %%) -- the factor error is far smaller because numerator and denominator move together" %
  (C.sd(hiso_r), 100 * C.sd(hiso_r) / lin["C2"]))

# ---- estimator variants (the factor that each variant would need)
var = {}
for name, kw in [("range +-35", dict(lo=-0.035, hi=0.035)), ("range +-55", dict(lo=-0.055, hi=0.055)), ("pol3", dict(order=3)),
                 ("sigma 5.76 MeV", dict(sigma=0.00576)), ("sigma 7.10 MeV", dict(sigma=0.00710)), ("E(m,0)", dict(useq=False)),
                 ("r = 0.89", dict(r=0.89)), ("r = 1.10", dict(r=1.10)), ("upper end +30", dict(hi=0.030)), ("upper end +60", dict(hi=0.060))]:
    Bv = build(**kw)
    cv = row_coef(Bv)
    mm = isC2
    c_all = factor(cv, Q["q_A"], mm)[0]; c1 = factor(cv, Q["q_A"], mm & E1)[0]; c2_ = factor(cv, Q["q_A"], mm & ~E1)[0]
    var[name] = dict(all=c_all, E1=c1, E2=c2_)
    P("   variant %-15s: c_u all %.4f, E1 %.4f, E2 %.4f" % (name, c_all, c1, c2_))
RES["variants"] = var

# ---- other definitions on the same rows (weights W x 1; S0, mode A; the edgeu u)
uC = B0["us"]; mppC = B0["mpp"]
inpp = (mppC >= 0.30) & (mppC < 0.40)
m4p = m4p_s0
defs = {}
def wmean(q, m):
    mm = m & isC2
    return float((w * q)[mm].sum() / w[mm].sum()), float(w[mm].sum())
for name, m in [("all C2 rows (S0, all M(pipi), all M4)", np.ones(S["n"], bool)),
                ("0.30 <= M(pipi) < 0.40, all M4", inpp),
                ("window B (0.4550 <= M4 < 0.5225)", (m4p >= 0.4550) & (m4p < 0.5225)),
                ("edge window -15...+5 MeV in u (ecalacc)", np.isfinite(uC) & (uC >= -0.015) & (uC < 0.005)),
                ("edge window -15...+5 MeV, 0.30-0.40", np.isfinite(uC) & (uC >= -0.015) & (uC < 0.005) & inpp),
                ("fit range -45...+45 MeV, 0.30-0.40", np.isfinite(uC) & (uC >= -0.045) & (uC < 0.045) & inpp),
                ("far below u < -60 MeV", np.isfinite(uC) & (uC < -0.060)),
                ("u < -15 MeV, 0.30-0.40", np.isfinite(uC) & (uC < -0.015) & inpp)]:
    for cl, mc in (("all", np.ones(S["n"], bool)), ("E1", E1), ("E2", ~E1)):
        v, Wsum = wmean(Q["q_A"], m & mc)
        defs["%s | %s" % (name, cl)] = dict(c=v, W=Wsum)
    P("   definition %-44s: all %.4f, E1 %.4f, E2 %.4f" % (name, defs["%s | all" % name]["c"], defs["%s | E1" % name]["c"], defs["%s | E2" % name]["c"]))
RES["definitions"] = defs
# where the estimator puts its weight: share of positive / negative coefficient weight, by u region
cw = coef0 * w
for lo_, hi_ in [(-0.045, -0.015), (-0.015, 0.0), (0.0, 0.015), (0.015, 0.045)]:
    m = isC2 & np.isfinite(uC) & (uC >= lo_) & (uC < hi_)
    P("   estimator weight of the C2 rows in u %+.3f..%+.3f: sum a W = %8.1f (sum |a| W %.1f), <w> of these rows %.4f" %
      (lo_, hi_, cw[m].sum(), np.abs(cw[m]).sum(), (w * Q["q_A"])[m].sum() / w[m].sum()))
# ---- class fits (E1, E2 fitted separately, as the edgeu class check and the kappa split)
P("   CLASS FITS (each M(ee) class fitted on its own):")
cls_out = {}
for cl in ("E1", "E2"):
    Bc = build(cls=cl)
    Rc = Bc["R"]
    l2 = C.h_lin(Rc, Bc["simpart"](["C2"], c2=False)); l2c = C.h_lin(Rc, Bc["simpart"](["C2"], c2=True)); l3 = C.h_lin(Rc, Bc["simpart"](["C3"], c2=False))
    cc = row_coef(Bc)
    cq_, hq_, hi_ = factor(cc, Q["q_A"], isC2)
    cls_out[cl] = dict(hd=Rc.h, C2=l2, C2c=l2c, C3=l3, cq=cq_, C2q=hq_, f_app=Rc.h / (l2c + l3), f_qed=Rc.h / (hq_ + l3), f_unc=Rc.h / (l2 + l3))
    P("      %s: h_data %.2f; C2 %.2f (applied %.2f, QED %.2f: factor %.4f), C3 %.2f; f applied %.4f, QED %.4f (%+.2f %%), uncorrected %.4f" %
      (cl, Rc.h, l2, l2c, hq_, cq_, l3, cls_out[cl]["f_app"], cls_out[cl]["f_qed"], 100 * (cls_out[cl]["f_qed"] / cls_out[cl]["f_app"] - 1), cls_out[cl]["f_unc"]))
def kappa_solve(A1, B1, A2, B2, h1, h2):
    rho = h1 / h2
    kap = (rho * A2 - A1) / (B1 - rho * B2)
    return h2 / (A2 + kap * B2), kap
e1, e2 = cls_out["E1"], cls_out["E2"]
ka = kappa_solve(e1["C2c"], e1["C3"], e2["C2c"], e2["C3"], e1["hd"], e2["hd"])
kq = kappa_solve(e1["C2q"], e1["C3"], e2["C2q"], e2["C3"], e1["hd"], e2["hd"])
P("      E1/E2 split: applied factors a_D %.4f, kappa %.4f (edgeu: 0.444, 0.876); QED weights a_D %.4f, kappa %.4f" % (ka[0], ka[1], kq[0], kq[1]))
RES["classes"] = cls_out
RES["kappa"] = dict(applied=ka, qed=kq)
# kappa = 0.48 reading
ak_app = R0.h / (lin["C2c"] + 0.48 * lin["C3"]); ak_q = R0.h / (h2q + 0.48 * lin["C3"])
P("   kappa = 0.48 reading: applied %.4f, QED %.4f (%+.2f %%)" % (ak_app, ak_q, 100 * (ak_q / ak_app - 1)))
RES["k48"] = dict(applied=ak_app, qed=ak_q)

# ---- windows (visible yields): f x (C2 + C3) with the C2 template in each reading
ref_w = json.load(open(C.OUT + "edgeu_ana.json"))["windows"]
inB = (m4p >= 0.4550) & (m4p < 0.5225)
yl = {}
for lab, sel_ in [("B", inB), ("all", np.ones(S["n"], bool))]:
    m2 = isC2 & sel_
    W2 = float(w[m2].sum()); W2q = float((w * Q["q_A"])[m2].sum())
    W2c = ref_w[lab]["C2c"]; W3 = ref_w[lab]["C3"]
    y_app = f_app * (W2c + W3) * C.PT3; y_q = f_new * (W2q + W3) * C.PT3
    yl[lab] = dict(C2=W2, C2q=W2q, C2c=W2c, C3=W3, Y_app=y_app, Y_qed=y_q, cwin=W2q / W2)
    P("   window %s: SMASH C2 %.1f (ref %.1f), applied factors %.1f, QED %.1f (window factor %.4f); C3 %.1f; Y applied %.0f, QED %.0f (%+.2f %%)" %
      (lab, W2, ref_w[lab]["C2"], W2c, W2q, W2q / W2, W3, y_app, y_q, 100 * (y_q / y_app - 1)))
RES["yields"] = yl

# ---- the factor and the momentum of the softer lepton; sensitivity to the soft-lepton efficiency
import c2recon_weights as CW
RR = CW.load_rows()
pe = np.sqrt(RR["ep_px"] ** 2 + RR["ep_py"] ** 2 + RR["ep_pz"] ** 2) * 1e-3
pm = np.sqrt(RR["em_px"] ** 2 + RR["em_py"] ** 2 + RR["em_pz"] ** 2) * 1e-3
psoft_r = np.minimum(pe, pm)
psoft = np.full(S["n"], np.nan); psoft[jw >= 0] = psoft_r[jw[jw >= 0]]
pbins = [0.0, 0.06, 0.08, 0.10, 0.15, 0.20, 0.30, 0.50, 2.0]
vs = []
for a_, b_ in zip(pbins[:-1], pbins[1:]):
    m = isC2 & (psoft >= a_) & (psoft < b_)
    if m.sum() == 0:
        continue
    vs.append(dict(lo=a_, hi=b_, n=int(m.sum()), Wshare=float(w[m].sum() / w[isC2].sum()), c=float((w * Q["q_A"])[m].sum() / w[m].sum()),
                   coefshare=float((coef0 * w)[m].sum() / (coef0 * w)[isC2].sum())))
P("   factor vs momentum of the softer lepton (all S0 C2 rows): " + "; ".join("%.2f-%.2f GeV/c: %.3f (W share %.3f, step share %.3f)" %
                                                                            (x["lo"], x["hi"], x["c"], x["Wshare"], x["coefshare"]) for x in vs))
RES["vs_psoft"] = vs
P("   median softer-lepton momentum of the S0 C2 rows: %.3f GeV/c (weighted)" % float(np.median(psoft[isC2])))
sens = {}
for lab, pcut, fac in [("eff(p_soft < 0.10) x 0.8", 0.10, 0.8), ("eff(p_soft < 0.10) x 1.2", 0.10, 1.2), ("eff(p_soft < 0.15) x 0.8", 0.15, 0.8),
                       ("eff(p_soft < 0.15) x 1.2", 0.15, 1.2)]:
    sw_ = np.where(psoft < pcut, fac, 1.0)
    num = (coef0 * w * sw_ * Q["q_A"])[isC2].sum(); den = (coef0 * w * sw_)[isC2].sum()
    sens[lab] = float(num / den)
    P("   sensitivity: %-26s -> c_u %.4f (%+.4f)" % (lab, num / den, num / den - out["all"]["q_A"]["c"]))
RES["sens_lepton_eff"] = sens

# ---- paired statistical error of the variant differences (are the variant factors different beyond noise?)
for name, kw in [("pol3", dict(order=3)), ("range +-35", dict(lo=-0.035, hi=0.035)), ("E(m,0)", dict(useq=False))]:
    cv = row_coef(build(**kw))
    a1 = (coef0 * w * Q["q_A"])[rows]; b1 = (coef0 * w)[rows]; a2 = (cv * w * Q["q_A"])[rows]; b2 = (cv * w)[rows]
    d = (K * a2).sum(1) / (K * b2).sum(1) - (K * a1).sum(1) / (K * b1).sum(1)
    P("   variant %-11s: c_u - c_u(nominal) = %+.4f +- %.4f (paired replicas)" % (name, a2.sum() / b2.sum() - a1.sum() / b1.sum(), C.sd(d)))
    RES["variants"][name]["diff_err"] = float(C.sd(d))

np.savez_compressed(O + "c2recon_rowcoef_edgeu.npz", coef=coef0[isC2], w=w[isC2], q=Q["q_A"][isC2], u=uC[isC2], mee=mee[isC2],
                    mpp=mppC[isC2], tree=S["tree"][isC2], grp=S["grp"][isC2])
json.dump(RES, open(O + "c2recon_edgeu.json", "w"), indent=1, default=float)
P("done (%.0f s)" % (time.time() - T0))
