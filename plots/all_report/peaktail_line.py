#!/usr/bin/env python3
# ============================================================================
# peaktail_line.py  --  package peaktail, task 2: the tail of the eta peak measured in the DATA
# ============================================================================
# Written 2026-11-15 (prefix peaktail). Two data constraints on the tail scale t (t = 1: the SMASH tail; tail = simulated
# peak rows whose leptons' momentum deficit shifts M4 by more than 5 MeV/c2, see peaktail_core.py / peaktail_ufit.py):
#   A. the line shape of the peak in M4 where the bump cannot contribute: M(pi+pi-) >= 0.415 GeV/c2 (no eta -> pi+pi-pi0 above
#      0.4116 on the data scale), N - CB (2 sqrt(N++ N--) in 10 MeV cells, interpolated), fit pol2/pol3 + peak, the peak =
#      SMASH C4 + C5 in M4' of the same region with the tail scaled by t; amplitude fixed (a_P) or free; core width as SMASH
#      (r = 1) or data-like (r = 0.893, applied to the core rows only);
#   B. the u fit itself per M(e+e-) class (E1: conversions + low-mass direct decays; E2 and >= 20 MeV: direct decays only),
#      a_P fixed, t free (the region 0 < u < +45 MeV, where the tail and the smooth part compete).
# Errors: fit errors (GLS, V = cov(CB) + diag(model)), checked with 200 data bootstrap copies for the main choices.
# Outputs: /hdd3/.../peaktail/peaktail_line.json, peaktail_line.log, peaktail_linefig.npz
# Run: python3 -I plots/all_report/peaktail_line.py
# ============================================================================
import json
import math
import sys

import numpy as np

sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import edgeu_core as C
import peaktail_core as P

LOG = open(P.OUT + "peaktail_line.log", "w")
def L(*a):
    s = " ".join(str(x) for x in a)
    print(s); LOG.write(s + "\n"); LOG.flush()
RES = {}
D = P.load_data(); S = P.load_sim()
w = P.sim_weight(S); gS = S["g"]
us, mpp, m4p = C.u_sim(S)
M4r, J = P.m4_and_jacobian(S)
dl = C.S_NOM * (J[0] * (S["ep_p"] - S["ep_ptrue"]) + J[1] * (S["em_p"] - S["em_ptrue"])).astype(np.float64) * 1e-3
flag = (dl < -0.005).astype(np.float64)
pk = (gS == 2) | (gS == 3)
S0 = (D["bits"] & 1) != 0
m4d = D["m4"].astype(np.float64); mppd = D["mpipi"].astype(np.float64); meed = D["mee"].astype(np.float64)

# ------------------------------------------------------------------------------------------------ A. M4 line shape
G4 = C.Grid(0.40, 0.0025, 120, 4)           # M4 0.40-0.70, 2.5 MeV bins, CB cells 10 MeV
def data_m4(sel, pplo, pphi, slices, clsmask=None):
    """N - CB in M4 for pplo <= M(pipi) < pphi; CB per M(pipi) slice (2 sqrt(N++ N--) in 10 MeV cells, interpolated)"""
    N = np.zeros(G4.nb); CB = np.zeros(G4.nb); COV = np.zeros((G4.nb, G4.nb))
    b = G4.bin(m4d)
    for a, c in slices:
        m = sel & (mppd >= a) & (mppd < c) & (b >= 0)
        if clsmask is not None:
            m &= clsmask
        n = np.bincount(b[m & (D["k"] == 0)], minlength=G4.nb).astype(float)
        pp = np.bincount(b[m & (D["k"] == 1)], minlength=G4.nb).astype(float)
        mm = np.bincount(b[m & (D["k"] == 2)], minlength=G4.nb).astype(float)
        cb, cov = C.cb_rule(pp, mm, G4)
        N += n; CB += cb; COV += cov
    return C.Spec(N, CB, COV)
def sim_m4(wt, pplo, pphi, groups=(2, 3), r=1.0, rcore_only=True, clsmask=None):
    """simulated rows in M4' (expected raw rows); r: resolution transfer of the core rows (tail rows keep r = 1 if
    rcore_only) -- the tail of radiation is not a resolution effect"""
    m4 = S["m4"].astype(np.float64); mt = S["m4t"].astype(np.float64)
    rr = np.where(flag > 0, 1.0, r) if rcore_only else np.full(S["n"], r)
    x = np.where(mt > 0, C.S_NOM * (mt + rr * (m4 - mt)), C.S_NOM * m4)
    b = G4.bin(x)
    m = np.isin(gS, groups) & (b >= 0) & (mpp >= pplo) & (mpp < pphi) & (wt != 0)
    if clsmask is not None:
        m &= clsmask
    return np.bincount(b[m], weights=wt[m], minlength=G4.nb).astype(float)
def parts(pplo, pphi, r=1.0, clsmask=None):
    """K(t) = K0 + t KK in M4' with the per-class normalisation over the region"""
    K0 = np.zeros(G4.nb); KK = np.zeros(G4.nb)
    for gg in (2, 3):
        m = (gS == gg) & (mpp >= pplo) & (mpp < pphi) & (w > 0)
        if clsmask is not None:
            m &= clsmask
        N = w[m].sum(); Nt = (w * flag)[m].sum(); Nc = N - Nt
        if N <= 0:
            continue
        kc = sim_m4(w * (1 - flag), pplo, pphi, (gg,), r, clsmask=clsmask); kt = sim_m4(w * flag, pplo, pphi, (gg,), r, clsmask=clsmask)
        K0 += kc * N / Nc; KK += kt - kc * Nt / Nc
    return K0, KK
def lfit(spec, cols_fixed, cols_free, lo, hi, order, Vfix=None):
    """GLS: N - CB - fixed = pol(order) + sum free columns; V = cov(CB) + diag(max(model + CB, 1)) iterated"""
    idx = np.where((G4.x >= lo - 1e-12) & (G4.x < hi - 1e-12))[0]
    x = G4.x[idx]; x0, dx = 0.5 * (lo + hi), 0.5 * (hi - lo); tt = (x - x0) / dx
    F = np.vstack([tt ** j for j in range(order + 1)] + [c[idx] for c in cols_free]).T
    fx = np.zeros(len(idx)) if cols_fixed is None else cols_fixed[idx]
    y = spec.sig[idx] - fx; ad = spec.cb[idx]; V0 = spec.cbcov[np.ix_(idx, idx)]
    var = np.maximum(y + fx + ad, 1.0)
    for it in range(8):
        V = V0 + np.diag(var) if Vfix is None else Vfix
        p, Cm, Lm, chi2 = C.gls(y, F, V)
        var = np.maximum(F @ p + fx + ad, 1.0)
    return dict(p=p, C=Cm, chi2=chi2, ndf=len(idx) - F.shape[1], idx=idx, model=F @ p + fx, smooth=F[:, :order + 1] @ p[:order + 1], V=V)
L("A. LINE SHAPE OF THE PEAK WHERE THE BUMP CANNOT CONTRIBUTE (data S0, N - CB in M4; SMASH C4 + C5 with the tail scaled by t)")
regions = {"0.415-0.46": (0.415, 0.46, [(0.415, 0.43), (0.43, 0.445), (0.445, 0.46)]),
           "0.42-0.46": (0.42, 0.46, [(0.42, 0.44), (0.44, 0.46)]),
           "0.415-0.50": (0.415, 0.50, [(0.415, 0.43), (0.43, 0.445), (0.445, 0.46), (0.46, 0.48), (0.48, 0.50)])}
lineA = {}
for rn, (a, b, sl) in regions.items():
    sp = data_m4(S0, a, b, sl)
    for r in (1.0, 0.893):
        K0, KK = parts(a, b, r)
        tailfrac = float((KK + K0).sum())
        for lo, hi, order in ((0.49, 0.60, 2), (0.48, 0.62, 3), (0.50, 0.59, 2), (0.47, 0.62, 3)):
            key = "%s|r%.3f|%.2f-%.2f pol%d" % (rn, r, lo, hi, order)
            # (i) amplitude fixed at a_P, t free
            R1 = lfit(sp, C.AP * K0, [C.AP * KK], lo, hi, order)
            t1, e1 = float(R1["p"][order + 1]), float(math.sqrt(R1["C"][order + 1, order + 1]))
            # (ii) amplitude and t free
            R2 = lfit(sp, None, [K0, KK], lo, hi, order)
            a0, a1 = float(R2["p"][order + 1]), float(R2["p"][order + 2]); Cm = R2["C"][order + 1:order + 3, order + 1:order + 3]
            t2 = a1 / a0; e2 = abs(t2) * math.sqrt(max(Cm[1, 1] / a1 ** 2 + Cm[0, 0] / a0 ** 2 - 2 * Cm[0, 1] / (a0 * a1), 0)) if a1 != 0 else float("nan")
            # (iii) t = 1 and t = 0, amplitude free: chi2 difference
            R3 = lfit(sp, None, [K0 + KK], lo, hi, order); R4 = lfit(sp, None, [K0], lo, hi, order)
            lineA[key] = dict(t_fixA=t1, et_fixA=e1, chi2_fixA=R1["chi2"], ndf_fixA=R1["ndf"], amp=a0, eamp=float(math.sqrt(Cm[0, 0])), t_freeA=t2, et_freeA=e2,
                              chi2_t1=R3["chi2"], chi2_t0=R4["chi2"], amp_t1=float(R3["p"][order + 1]), peak_rows=float(C.AP * (K0 + KK).sum()))
            L("   %-38s a_P fixed: t = %5.2f +- %4.2f (chi2 %5.1f/%d) | amp. free: amp %.3f +- %.3f, t = %5.2f +- %4.2f | chi2(t=1) %5.1f chi2(t=0) %5.1f (amp %.3f)" %
              (key, t1, e1, R1["chi2"], R1["ndf"], a0, math.sqrt(Cm[0, 0]), t2, e2, R3["chi2"], R4["chi2"], R3["p"][order + 1]))
    # signal over background in the peak window, for the record
    K0, KK = parts(a, b, 1.0)
    iP = (G4.x >= 0.5336) & (G4.x < 0.5542)
    L("   %s: peak window P: N %.0f, CB %.0f, expected peak a_P x SMASH %.0f rows (S/B %.3f)" % (rn, sp.N[iP].sum(), sp.cb[iP].sum(),
      C.AP * (K0 + KK)[iP].sum(), C.AP * (K0 + KK)[iP].sum() / max(sp.sig[iP].sum() - C.AP * (K0 + KK)[iP].sum(), 1)))
RES["lineA"] = lineA
# tail share of the region in SMASH (rows with dl < -5 MeV) and the fraction below -15 MeV
for rn, (a, b, sl) in regions.items():
    m = pk & (mpp >= a) & (mpp < b) & (w > 0)
    d = P.resid(S)
    L("   SMASH in %s: tail share %.3f; d < -15 MeV %.3f; d < -25 MeV %.3f; sum %.1f expected raw rows" %
      (rn, (w * flag)[m].sum() / w[m].sum(), P.frac(d[m] < -0.015, w[m]), P.frac(d[m] < -0.025, w[m]), w[m].sum()))

# bootstrap of the main line-shape choice
rn0, lo0, hi0, or0 = "0.415-0.46", 0.49, 0.60, 2
a, b, sl = regions[rn0]
K0, KK = parts(a, b, 1.0)
sel_rows = np.where(S0 & (mppd >= a) & (mppd < b) & (G4.bin(m4d) >= 0))[0]
key = D["k"][sel_rows].astype(np.int64) * (1 << 32) + D["grp"][sel_rows].astype(np.int64)
uq, inv = np.unique(key, return_inverse=True)
rng = np.random.default_rng(0x7A11)
tb = []
for it in range(200):
    kk = rng.poisson(1.0, uq.size)[inv].astype(float)
    N = np.zeros(G4.nb); CB = np.zeros(G4.nb); COV = np.zeros((G4.nb, G4.nb))
    bb = G4.bin(m4d[sel_rows]); mp_ = mppd[sel_rows]; kd = D["k"][sel_rows]
    for aa, cc in sl:
        m = (mp_ >= aa) & (mp_ < cc)
        n = np.bincount(bb[m & (kd == 0)], weights=kk[m & (kd == 0)], minlength=G4.nb)
        pp = np.bincount(bb[m & (kd == 1)], weights=kk[m & (kd == 1)], minlength=G4.nb)
        mm = np.bincount(bb[m & (kd == 2)], weights=kk[m & (kd == 2)], minlength=G4.nb)
        cb, cov = C.cb_rule(pp, mm, G4); N += n; CB += cb; COV += cov
    R1 = lfit(C.Spec(N, CB, COV), C.AP * K0, [C.AP * KK], lo0, hi0, or0)
    tb.append(float(R1["p"][or0 + 1]))
RES["lineA_boot"] = dict(sd=C.sd(tb), mean=float(np.mean(tb)))
L("   bootstrap (200) of %s %.2f-%.2f pol%d, a_P fixed: t sd %.2f (fit %.2f), mean %.2f" % (rn0, lo0, hi0, or0, C.sd(tb),
  lineA["%s|r1.000|%.2f-%.2f pol%d" % (rn0, lo0, hi0, or0)]["et_fixA"], np.mean(tb)))
# figure data: the main choice
sp = data_m4(S0, a, b, sl)
R1 = lfit(sp, C.AP * K0, [C.AP * KK], lo0, hi0, or0)
R3 = lfit(sp, None, [K0 + KK], lo0, hi0, or0); R4 = lfit(sp, None, [K0], lo0, hi0, or0)
figA = dict(x=G4.x, sig=sp.sig, err=np.sqrt(np.maximum(sp.N, 1) + np.diag(sp.cbcov)), cb=sp.cb, K0=K0, KK=KK, idx=R1["idx"], m1=R1["model"], s1=R1["smooth"],
            p1=R1["p"], m3=R3["model"], s3=R3["smooth"], p3=R3["p"], m4_=R4["model"], s4=R4["smooth"], p4=R4["p"])

# ------------------------------------------------------------------------------------------------ B. the u fit per M(ee) class
L("\nB. THE u FIT PER M(e+e-) CLASS (edgeu estimator, class CB by the like-sign share rule; a_P fixed, t free; class peak templates)")
Ded = C.load_data(); ud = C.u_data(Ded); S0e = (Ded["bits"] & 1) != 0
g = C.GRID_NOM
H = C.hist_data(Ded, ud, S0e, g); SL = C.slices_of(0.30, 0.40)
import json as _j
A_ = _j.load(open(C.OUT + "edgeu_ana.json")); T = C.step_gauss(g, A_["reso"]["sigma_core"])
wc = P.sim_weight(S, c2corr=True)
pa, pb = C.PB(0.30), C.PB(0.40)
def utm(wt, groups, clsm):
    bb = g.bin(us); ppb = C.ppbin_of(mpp.astype(np.float32))
    m = np.isin(gS, groups) & (bb >= 0) & (ppb >= pa) & (ppb < pb) & (wt != 0) & clsm
    return np.bincount(bb[m], weights=wt[m], minlength=g.nb).astype(float)
inR = (mpp >= 0.30) & (mpp < 0.40)
CLS = {"all": np.ones(S["n"], bool), "E1": S["mee"] < 0.010, "E2": S["mee"] >= 0.010, "E20": S["mee"] >= 0.020}
lineB = {}
for cn, cm in CLS.items():
    sp = C.data_spec(H, SL, cn, g)
    ys = utm(wc, (0, 1), cm)
    K0 = np.zeros(g.nb); KK = np.zeros(g.nb)
    for gg in (2, 3):
        m = (gS == gg) & inR & (w > 0) & cm
        N = w[m].sum(); Nt = (w * flag)[m].sum(); Nc = N - Nt
        if N <= 0 or Nc <= 0:      # e.g. no simulated C4 row with M(ee) >= 20 MeV
            continue
        kc = utm(w * (1 - flag), (gg,), cm); kt = utm(w * flag, (gg,), cm)
        K0 += kc * N / Nc; KK += kt - kc * Nt / Nc
    Rn = C.fit_step(sp, T, g, fixed=C.AP * (K0 + KK)); fn = Rn.h / C.h_lin(Rn, ys)
    out = dict(f_nom=fn)
    for lab, lo, hi, order in (("-45..+45 pol2", -0.045, 0.045, 2), ("-45..+60 pol2", -0.045, 0.060, 2), ("-55..+70 pol3", -0.055, 0.070, 3)):
        R = C.fit_step(sp, T, g, lo, hi, order, fixed=C.AP * K0, free_cols=[C.AP * KK])
        it = order + 1
        out[lab] = dict(t=float(R.p[it]), et=float(math.sqrt(R.C[it, it])), f=R.h / C.h_lin(R, ys), chi2=R.chi2, ndf=R.ndf)
        L("   %-4s %-14s t = %5.2f +- %4.2f, f %.4f (nominal f of the class %.4f), chi2 %.1f/%d" % (cn, lab, R.p[it], math.sqrt(R.C[it, it]), R.h / C.h_lin(R, ys), fn, R.chi2, R.ndf))
    # projection of the tail on the class step per unit t
    out["dfdt_rel"] = C.h_lin(Rn, C.AP * KK) / Rn.h * -1.0
    lineB[cn] = out
    L("   %-4s df/dt (relative, fixed a_P) = %+.4f per unit t" % (cn, out["dfdt_rel"]))
RES["lineB"] = lineB
json.dump(RES, open(P.OUT + "peaktail_line.json", "w"), indent=1)
np.savez_compressed(P.OUT + "peaktail_linefig.npz", **figA)
L("\nwritten: peaktail_line.json, peaktail_line.log, peaktail_linefig.npz")
