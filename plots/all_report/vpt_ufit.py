#!/usr/bin/env python3
"""vpt_ufit.py -- verifier of peaktail: (1) tail scale t from the u fit (own code, own tail flags), (2) f_eta(u-edge)
versus t, (3) the checker's Gaussian peak term (chk_hsim.py logic re-coded) and its matched version.
Output: stdout (tee vpt_ufit.log) and vpt_ufit.json in /hdd3/przygoda/HADES/FAT_out/all_report/vpt/.
"""
import json
import math
import sys
import numpy as np
from scipy.special import erf
sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import vpt_core as V
import vpt_tmpl as TM
from vpt_nom import data_u, data_spec, SL, in_sl

g = V.Grid()
Z = np.load(V.OUT + "vpt_nom.npz")
N, CB, COV, T = Z["N"], Z["CB"], Z["COV"], Z["T"]
sig = N - CB
BUMP = Z["C2"] + Z["C3"]
R = TM.rows()
reg = TM.in_pp(R["mpp"], 0.30, 0.40) & np.isfinite(R["u"])
out = {}


def fit(fixed=None, extra=None, lo=-0.045, hi=0.045, order=2, s=sig, cb=CB, cov=COV):
    return V.fitstep(s, cb, cov, T, g, lo=lo, hi=hi, order=order, fixed=fixed, extra=extra)


def f_of(Rf, bump=BUMP):
    return Rf.h / V.hlin(Rf, bump)


# ------------------------------------------------------------------------------------------ nominal check
TT = {k: TM.templates(R["u"], g.e, reg, R["flag"][k]) for k in R["flag"]}
K1 = TM.K_of_t(TT["sub5"], 1.0)
R0 = fit(fixed=V.AP * K1)
print("nominal (K(1) from vpt_tmpl): h %.4f f %.5f chi2 %.2f/%d" % (R0.h, f_of(R0), R0.chi2, R0.ndf))
for k in TT:
    T4, T5 = TT[k]["C4"], TT[k]["C5"]
    print("  tail flag %-6s: C4 tail share %.3f, C5 %.3f (0.30-0.40, finite u)" % (k, T4[5] / T4[3], T5[5] / T5[3]))
# core / tail projection onto h per unit a_P (rows of step)
for k in ("sub5",):
    Ka, Kc, Kt, Nn, Nc, Nt = TT[k]["C4"]; Ka5, Kc5, Kt5, N5, Nc5, Nt5 = TT[k]["C5"]
    print("  projection on h per unit a_P (%s): core C4 %+.1f C5 %+.1f, tail C4 %+.1f C5 %+.1f, total %+.2f; "
          "no tail (core x N/Ncore): %+.1f" % (k, V.hlin(R0, Kc), V.hlin(R0, Kc5), V.hlin(R0, Kt), V.hlin(R0, Kt5),
                                               V.hlin(R0, Ka + Ka5), V.hlin(R0, TM.K_of_t(TT[k], 0.0))))

# ------------------------------------------------------------------------------------------ (1) t from the u fit
FITS = [(-0.045, 0.045, 2), (-0.045, 0.055, 2), (-0.045, 0.060, 2), (-0.045, 0.070, 3), (-0.055, 0.070, 3), (-0.040, 0.050, 2)]
NB = 200
Wb = TM.group_boot(NB)
print("\n(1) tail scale t from the u fit, a_P fixed (column a_P dK/dt free); stat = fit error; MC = sd over %d sim group bootstraps" % NB)
res1 = {}
for k in ("sub5", "sub10", "p95", "res15"):
    Tk = TT[k]
    Kb = TM.K_of_t(Tk, 1.0); dKk = TM.dK(Tk)
    # bootstrap templates
    Kr, dKr = [], []
    for r in range(NB):
        Tr = TM.templates(R["u"], g.e, reg, R["flag"][k], w=R["w"] * Wb[r])
        Kr.append(TM.K_of_t(Tr, 1.0)); dKr.append(TM.dK(Tr))
    for lo, hi, od in FITS:
        Rf = fit(fixed=V.AP * Kb, extra=[V.AP * dKk], lo=lo, hi=hi, order=od)
        c, ec = Rf.p[-2], math.sqrt(Rf.C[-2, -2])
        rho = Rf.C[-2, -1] / math.sqrt(Rf.C[-2, -2] * Rf.C[-1, -1])
        cr = []
        for r in range(NB):
            Rr = fit(fixed=V.AP * Kr[r], extra=[V.AP * dKr[r]], lo=lo, hi=hi, order=od)
            cr.append(Rr.p[-2])
        mc = float(np.std(cr, ddof=1))
        # amplitude free: chi2 profile over t
        tg = np.linspace(-2, 4, 121); ch = []
        for t in tg:
            ch.append(fit(extra=[TM.K_of_t(Tk, t)], lo=lo, hi=hi, order=od).chi2)
        ch = np.array(ch); j = int(np.argmin(ch))
        inside = tg[ch <= ch[j] + 1.0]
        res1[(k, lo, hi, od)] = dict(t=1 + c, stat=ec, mc=mc, rho=rho, chi2=Rf.chi2, ndf=Rf.ndf, tfree=float(tg[j]),
                                     tfree_lo=float(inside.min()), tfree_hi=float(inside.max()))
        print("  %-6s %+3.0f..%+3.0f pol%d: t = %5.2f +- %.2f (stat) +- %.2f (MC)  rho(t,h) %+.2f  chi2 %.1f/%d | amplitude free: "
              "t = %5.2f [%5.2f, %5.2f]" % (k, lo * 1e3, hi * 1e3, od, 1 + c, ec, mc, rho, Rf.chi2, Rf.ndf, tg[j], inside.min(), inside.max()))
out["ufit_t"] = {"%s_%g_%g_%d" % kk: v for kk, v in res1.items()}

# per M(ee) class (E1 < 10 MeV, E2 >= 10 MeV): own CB per class (own-class rule: like-sign rows of that class)
D = V.load_data(); ud = data_u(D); S0 = (D["bits"] & 1) != 0
print("\n    per M(ee) class (own-class CB; sim templates of the class), sub5, a_P fixed:")
for cn, cm_d, cm_s in (("E1", D["mee"] < 0.010, R["S"]["ee"] < 0.010), ("E2", ~(D["mee"] < 0.010), ~(R["S"]["ee"] < 0.010))):
    Nc, CBc, COVc = data_spec(D, ud, S0 & cm_d)
    Tc = TM.templates(R["u"], g.e, reg & cm_s, R["flag"]["sub5"])
    for lo, hi, od in ((-0.045, 0.045, 2), (-0.045, 0.060, 2)):
        Rf = fit(fixed=V.AP * TM.K_of_t(Tc, 1.0), extra=[V.AP * TM.dK(Tc)], lo=lo, hi=hi, order=od, s=Nc - CBc, cb=CBc, cov=COVc)
        print("     %s %+3.0f..%+3.0f pol%d: t = %5.2f +- %.2f (stat)  chi2 %.1f/%d" % (cn, lo * 1e3, hi * 1e3, od, 1 + Rf.p[-2], math.sqrt(Rf.C[-2, -2]), Rf.chi2, Rf.ndf))

# ------------------------------------------------------------------------------------------ (2) f_eta(u-edge) versus t
print("\n(2) f_eta(u-edge) with the fixed peak term a_P K(t) (nominal fit, iterated V; h_sim with each fit's estimator)")
f0 = f_of(R0)
fsl = {}
for k in ("sub5", "sub10", "p95", "res15"):
    vals = []
    for t in (0.0, 0.34, 0.5, 0.72, 1.0, 1.10, 1.5, 2.0):
        Rt = fit(fixed=V.AP * TM.K_of_t(TT[k], t))
        vals.append((t, f_of(Rt), Rt.chi2))
    sl = np.polyfit([v[0] for v in vals], [v[1] / f0 for v in vals], 1)[0]
    fsl[k] = sl
    print("  %-6s: " % k + "  ".join("t=%.2f f=%.4f (%+.1f%%, chi2 %.1f)" % (t, f, 100 * (f / f0 - 1), c2) for t, f, c2 in vals)
          + "   slope d ln f / dt = %+.4f" % sl)
    # per class
    for cl in ("C4", "C5"):
        Rt = fit(fixed=V.AP * (TM.K_of_t(TT[k], 0.0, 1.0) if cl == "C4" else TM.K_of_t(TT[k], 1.0, 0.0)))
        print("          %s alone at t = 0 (other class t = 1): %+.2f %%" % (cl, 100 * (f_of(Rt) / f0 - 1)))
out["f_slope"] = fsl
# template noise of the shift at t = 0
sh = []
for r in range(100):
    Tr = TM.templates(R["u"], g.e, reg, R["flag"]["sub5"], w=R["w"] * Wb[r])
    a = f_of(fit(fixed=V.AP * TM.K_of_t(Tr, 0.0))); b = f_of(fit(fixed=V.AP * TM.K_of_t(Tr, 1.0)))
    sh.append(a / b - 1)
print("  MC noise of the relative shift f(0)/f(1) - 1 (sub5): mean %+.4f sd %.4f" % (np.mean(sh), np.std(sh, ddof=1)))

# ------------------------------------------------------------------------------------------ (3) the checker's estimator
print("\n(3) the checker's Gaussian peak term, re-coded from chk_hsim.py's description")
m4d = D["m4"].astype(float); mp = D["mpipi"]
base = S0 & (D["k"] == 0) & in_sl(mp, 0.30, 0.40) & (D["mee"] < V.MPI0)
inP = base & (m4d >= 0.53364) & (m4d < 0.55416)
inSB = base & (((m4d >= 0.5225) & (m4d < 0.53364)) | ((m4d >= 0.55416) & (m4d < 0.56530)))
rr = 20.52 / 22.28
E = V.edge(mp.astype(float), D["mee"].astype(float))


def gauss_tmpl(cen, wts, s):
    K = np.array([(wts * 0.5 * (erf((g.e[j + 1] - cen) / (math.sqrt(2) * s)) - erf((g.e[j] - cen) / (math.sqrt(2) * s)))).sum() for j in range(g.nb)])
    return K / K.sum()


cen_d = V.PEAKM - np.concatenate([E[inP], E[inSB]]); w_d = np.concatenate([np.ones(inP.sum()), -rr * np.ones(inSB.sum())])
print("  data positions: P rows %d, SB rows %d x %.4f, net %.0f" % (inP.sum(), inSB.sum(), rr, inP.sum() - rr * inSB.sum()))
Kd54 = gauss_tmpl(cen_d, w_d, 0.0054)
Rk = fit(extra=[Kd54])
fchk = f_of(Rk)
print("  checker: h_data %.2f +- %.2f chi2 %.2f/%d amplitude %.0f +- %.0f h_sim %.2f f %.5f (%+.2f %%)" % (
    Rk.h, Rk.eh, Rk.chi2, Rk.ndf, Rk.p[-2], math.sqrt(Rk.C[-2, -2]), V.hlin(Rk, BUMP), fchk, 100 * (fchk / f0 - 1)))
LK = V.hlin(Rk, V.AP * K1)
print("  the checker's estimator applied to the SMASH peak a_P K(1): %+.1f rows (core %+.1f, tail %+.1f, sub5)" % (
    LK, V.hlin(Rk, V.AP * (TT["sub5"]["C4"][1] + TT["sub5"]["C5"][1])), V.hlin(Rk, V.AP * (TT["sub5"]["C4"][2] + TT["sub5"]["C5"][2]))))
fm = (Rk.h - LK) / V.hlin(Rk, BUMP)
print("  matched checker f = (h_data - L.(a_P K))/L.bump = %.5f (%+.2f %%)" % (fm, 100 * (fm / f0 - 1)))
for t in (0.0, 0.72, 2.0):
    LKt = V.hlin(Rk, V.AP * TM.K_of_t(TT["sub5"], t))
    print("     matched with K(t = %.2f): %.5f (%+.2f %%)" % (t, (Rk.h - LKt) / V.hlin(Rk, BUMP), 100 * ((Rk.h - LKt) / V.hlin(Rk, BUMP) / f0 - 1)))
# decomposition
Rs = R["S"]; m45 = R["m45"] & reg
cen_s = V.PEAKM - V.edge(R["mpp"][m45], Rs["ee"][m45].astype(float))   # SMASH rows at their own edge-based positions
Ks54 = gauss_tmpl(cen_s, R["w"][m45], 0.0054)
Ks61 = gauss_tmpl(cen_s, R["w"][m45], 0.00609)
Kd61 = gauss_tmpl(cen_d, w_d, 0.00609)
rows_ = [("SMASH peak shape, amplitude free", dict(extra=[K1])),
         ("SMASH no tail (sub5, t = 0), a_P fixed", dict(fixed=V.AP * TM.K_of_t(TT["sub5"], 0.0))),
         ("SMASH no tail (sub5, t = 0), amplitude free", dict(extra=[TM.K_of_t(TT["sub5"], 0.0)])),
         ("no peak term", dict()),
         ("Gaussian 5.4 MeV at SMASH positions, free", dict(extra=[Ks54])),
         ("Gaussian 6.09 MeV at SMASH positions, free", dict(extra=[Ks61])),
         ("Gaussian 6.09 MeV at data positions, free", dict(extra=[Kd61])),
         ("checker Gaussian, amplitude fixed at a_P x SMASH sum", dict(fixed=Kd54 * V.AP * K1.sum()))]
for name, kw in rows_:
    Rx = fit(**kw)
    print("   %-55s h %.1f  f %.4f (%+.1f %%)" % (name, Rx.h, f_of(Rx), 100 * (f_of(Rx) / f0 - 1)))
# where do the data positions differ from the SMASH ones
def wmean(x, w):
    return float((x * w).sum() / w.sum())
print("  peak position u_p = 0.5439 - E: data (P - SB) mean %.2f MeV, SMASH rows mean %.2f MeV" % (1e3 * wmean(cen_d, w_d), 1e3 * wmean(cen_s, R["w"][m45])))
# data bootstrap of the checker's shift
Wd = np.random.default_rng(7).poisson(1.0, size=(100, 1))
keyd = D["k"].astype(np.int64) * (1 << 32) + D["grp"].astype(np.int64)
sel = S0 & in_sl(mp, 0.30, 0.40)
uq, inv = np.unique(keyd[sel], return_inverse=True)
rng = np.random.default_rng(11)
shifts = []
idx_sel = np.where(sel)[0]
for r in range(60):
    kw = np.zeros(len(D["k"])); kw[idx_sel] = rng.poisson(1.0, uq.size)[inv]
    Nb = np.zeros(g.nb); CBb = np.zeros(g.nb); COVb = np.zeros((g.nb, g.nb))
    for a, b in SL:
        m = S0 & in_sl(mp, a, b) & np.isfinite(ud)
        n = g.h(ud[m & (D["k"] == 0)], kw[m & (D["k"] == 0)]); p = g.h(ud[m & (D["k"] == 1)], kw[m & (D["k"] == 1)])
        q = g.h(ud[m & (D["k"] == 2)], kw[m & (D["k"] == 2)])
        c, cv = V.cb_cells(p, q, g); Nb += n; CBb += c; COVb += cv
    sb = Nb - CBb
    wP = kw[inP]; wS = kw[inSB]
    Kb = gauss_tmpl(cen_d, np.concatenate([wP, -rr * wS]), 0.0054)
    a = f_of(V.fitstep(sb, CBb, COVb, T, g, extra=[Kb])); b = f_of(V.fitstep(sb, CBb, COVb, T, g, fixed=V.AP * K1))
    shifts.append(a - b)
print("  data bootstrap (60) of f(checker) - f(nominal): mean %+.4f sd %.4f" % (np.mean(shifts), np.std(shifts, ddof=1)))
out.update(dict(f0=f0, fchk=fchk, fmatched=fm, LK=LK, chk_shift_sd=float(np.std(shifts, ddof=1))))
json.dump(out, open(V.OUT + "vpt_ufit.json", "w"), indent=1)
