#!/usr/bin/env python3
# ============================================================================
# peaktail_prop.py  --  package peaktail, tasks 2-3: the data-constrained tail scale t and its propagation into
#                       f_eta(u-edge) with the edgeu estimator; arrays for the per-slice edge (bedge, peaktail_bedge.C)
# ============================================================================
# Written 2026-11-15 (prefix peaktail). Reads peaktail_ufit.json, peaktail_line.json, peaktail_rad.json and the columns.
#   A. t from the data: (i) the u fit (0.30-0.40, the region 0 < u < +45..+70 MeV), five fit choices; (ii) the M4 line shape
#      where the bump cannot contribute (M(pipi) 0.42-0.46), four fit choices; template noise of t from the WP0 replicas;
#      the peak composition of Part III (conversion route 0.237, direct route 0.497 instead of 0.418 for both) as a variant.
#   B. f_eta(u-edge) (edgeu estimator, a_P fixed, h_sim = the same linear estimator on C2 + C3) for: SMASH tail, t_data +- 1
#      sigma, no tail (t = 0), the checker's Gaussian, doubled tail, the QED radiative tail of C5 (peaktail_rad.py), the
#      composition; and each with the noise of the difference (simulation replicas).
#   C. bedge input arrays: C4 and C5 of the nominal bedge configuration split into core and tail rows, and C5 with the
#      radiation toy, in the bedge binning (M4' 0.35-0.95 / 2.5 MeV, M(pipi)' 0.27 + 5 MeV x 40, class E1 / E2 / no value).
# Outputs: /hdd3/.../peaktail/peaktail_prop.json, peaktail_prop.log, bedge_C45_{core,tail,rad}.bin, bedge_C45_check.txt
# Run: python3 -I plots/all_report/peaktail_prop.py
# ============================================================================
import json
import math
import struct
import sys

import numpy as np
from scipy.special import erf

sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import edgeu_core as C
import peaktail_core as P

LOG = open(P.OUT + "peaktail_prop.log", "w")
def L(*a):
    s = " ".join(str(x) for x in a)
    print(s); LOG.write(s + "\n"); LOG.flush()
RES = {}
UF = json.load(open(P.OUT + "peaktail_ufit.json")); LN = json.load(open(P.OUT + "peaktail_line.json")); RD = json.load(open(P.OUT + "peaktail_rad.json"))
A = json.load(open(C.OUT + "edgeu_ana.json")); FD = dict(np.load(C.OUT + "edgeu_figdata.npz"))
g = C.GRID_NOM; T = C.step_gauss(g, A["reso"]["sigma_core"])
D = C.load_data(); u = C.u_data(D); S0 = (D["bits"] & 1) != 0
H = C.hist_data(D, u, S0, g); SL = C.slices_of(0.30, 0.40); sp = C.data_spec(H, SL, "all", g)
ysim = FD["sim_C2"] + FD["sim_C3"]
S = P.load_sim(); w = P.sim_weight(S); gS = S["g"]
us, mpp, m4p = C.u_sim(S)
M4r, J = P.m4_and_jacobian(S)
dl = C.S_NOM * (J[0] * (S["ep_p"] - S["ep_ptrue"]) + J[1] * (S["em_p"] - S["em_ptrue"])).astype(np.float64) * 1e-3
flag = (dl < -0.005).astype(np.float64)
inR = (mpp >= 0.30) & (mpp < 0.40)
pa, pb = C.PB(0.30), C.PB(0.40)
bU = g.bin(us); ppbU = C.ppbin_of(mpp.astype(np.float32))
AP4, AP5 = 0.237, 0.497          # Part III peak composition (fix_numbers.txt): conversion route, direct route
def utm(wt, groups, extra=None):
    m = np.isin(gS, groups) & (bU >= 0) & (ppbU >= pa) & (ppbU < pb) & (wt != 0)
    if extra is not None:
        m &= extra
    return np.bincount(bU[m], weights=wt[m], minlength=g.nb).astype(float)
def parts_u(wt=None):
    """per class: K(t) = K0 + t KK (u, 0.30-0.40), normalisation per class over 0.30-0.40"""
    wt = w if wt is None else wt
    out = {}
    for gg in (2, 3):
        m = (gS == gg) & inR & (wt > 0)
        N = wt[m].sum(); Nt = (wt * flag)[m].sum(); Nc = N - Nt
        kc = utm(wt * (1 - flag), (gg,)); kt = utm(wt * flag, (gg,))
        out[gg] = (kc * N / Nc, kt - kc * Nt / Nc)
    return out
PU = parts_u()
def peak_term(t4, t5, a4=C.AP, a5=C.AP, pu=None):
    pu = PU if pu is None else pu
    return a4 * (pu[2][0] + t4 * pu[2][1]) + a5 * (pu[3][0] + t5 * pu[3][1])
def est(fixed=None, free=None, lo=-0.045, hi=0.045, order=2, spec=None):
    R = C.fit_step(sp if spec is None else spec, T, g, lo, hi, order, fixed=fixed, free_cols=free)
    return R, R.h / C.h_lin(R, ysim)
R0, f0 = est(fixed=peak_term(1, 1))
L("0. nominal reproduced: f = %.5f (edgeu 0.42004)" % f0)

# ================================================================================================ A. t from the data
L("\nA. THE TAIL SCALE t FROM THE DATA")
UCH = [("-45..+45 pol2", -0.045, 0.045, 2), ("-45..+55 pol2", -0.045, 0.055, 2), ("-45..+60 pol2", -0.045, 0.060, 2),
       ("-45..+70 pol3", -0.045, 0.070, 3), ("-55..+70 pol3", -0.055, 0.070, 3)]
rows_pk = np.where(((gS == 2) | (gS == 3)) & inR & (w > 0))[0]
RW = C.sim_replica_weights(S, rows_pk, 100)
PUR = []
for r in range(100):
    wr = w.copy(); wr[rows_pk] = w[rows_pk] * RW[r]
    PUR.append(parts_u(wr))
def tfit_u(lo, hi, order, a4=C.AP, a5=C.AP, pu=None, spec=None):
    pu = PU if pu is None else pu
    fixed = a4 * pu[2][0] + a5 * pu[3][0]; col = a4 * pu[2][1] + a5 * pu[3][1]
    R = C.fit_step(sp if spec is None else spec, T, g, lo, hi, order, fixed=fixed, free_cols=[col])
    return float(R.p[order + 1]), float(math.sqrt(R.C[order + 1, order + 1])), R
tu = {}
for ck, comp, a4, a5 in (("SMASH", "SMASH composition (a_P 0.418 both)", C.AP, C.AP), ("PartIII", "Part III composition (0.237 C4, 0.497 C5)", AP4, AP5)):
    for lab, lo, hi, order in UCH:
        t, et, R = tfit_u(lo, hi, order, a4, a5)
        # template noise: replica templates
        tr = [tfit_u(lo, hi, order, a4, a5, PUR[r])[0] for r in range(100)]
        tu["%s|%s" % (ck, lab)] = dict(t=t, stat=et, mc=C.sd(tr), mc_mean=float(np.mean(tr)), chi2=R.chi2, ndf=R.ndf)
        L("   u fit, %-42s %-14s t = %5.2f +- %4.2f (stat) +- %4.2f (MC templates; replica mean %.2f), chi2 %.1f/%d" %
          (comp, lab, t, et, C.sd(tr), np.mean(tr), R.chi2, R.ndf))
RES["t_u"] = tu
# per class with the Part III composition
L("   per M(ee) class (u fit, Part III composition):")
CLSM = {"E1": S["mee"] < 0.010, "E2": S["mee"] >= 0.010}
wc = P.sim_weight(S, c2corr=True)
tcls = {}
for cn, cm in CLSM.items():
    spc = C.data_spec(H, SL, cn, g)
    pu = {}
    for gg in (2, 3):
        m = (gS == gg) & inR & (w > 0) & cm
        N = w[m].sum(); Nt = (w * flag)[m].sum(); Nc = N - Nt
        kc = utm(w * (1 - flag), (gg,), cm); kt = utm(w * flag, (gg,), cm)
        pu[gg] = (kc * N / Nc, kt - kc * Nt / Nc)
    for comp, a4, a5 in (("SMASH", C.AP, C.AP), ("PartIII", AP4, AP5)):
        for lab, lo, hi, order in UCH[:3]:
            t, et, R = tfit_u(lo, hi, order, a4, a5, pu, spc)
            tcls["%s|%s|%s" % (cn, comp, lab)] = dict(t=t, stat=et, chi2=R.chi2, ndf=R.ndf)
            L("     %s %-7s %-14s t = %5.2f +- %4.2f (chi2 %.1f/%d)" % (cn, comp, lab, t, et, R.chi2, R.ndf))
RES["t_u_class"] = tcls
# control region (peaktail_line.py), 0.42-0.46, amplitude free; scale factor sqrt(chi2/ndf)
L("   M4 line shape where no bump can contribute (peaktail_line.py, M(pipi) 0.42-0.46, amplitude free; error x sqrt(chi2/ndf)):")
tc = {}
for key, v in LN["lineA"].items():
    if not key.startswith("0.42-0.46"):
        continue
    ndf = int(key.split("pol")[1]) + 1
    # ndf of the amplitude-free fit = ndf_fixA - 1
    sf = math.sqrt(max(v["chi2_t1"] / (v["ndf_fixA"] - 1), 1.0))
    tc[key] = dict(t=v["t_freeA"], stat=v["et_freeA"], stat_scaled=v["et_freeA"] * sf, sf=sf, amp=v["amp"], t_fixA=v["t_fixA"])
    L("     %-38s t = %5.2f +- %4.2f (x %.2f -> %4.2f); amplitude %.3f; with a_P fixed t = %.2f" % (key, v["t_freeA"], v["et_freeA"], sf, v["et_freeA"] * sf, v["amp"], v["t_fixA"]))
RES["t_ctrl"] = tc
# combination
uvals = [tu["SMASH|%s" % c[0]]["t"] for c in UCH]
u_main = tu["SMASH|-45..+60 pol2"]
t_u = float(np.mean(uvals)); s_u_stat = math.hypot(u_main["stat"], u_main["mc"]); s_u_ch = float(np.std(uvals))
cvals = [tc[k]["t"] for k in tc if "r1.000" in k]
c_main = tc["0.42-0.46|r1.000|0.49-0.60 pol2"]
c_r = tc["0.42-0.46|r0.893|0.49-0.60 pol2"]
t_c = float(np.mean(cvals)); s_c_stat = c_main["stat_scaled"]; s_c_ch = math.hypot(float(np.std(cvals)), c_r["t"] - c_main["t"])
s_u = math.hypot(s_u_stat, s_u_ch); s_c = math.hypot(s_c_stat, s_c_ch)
wu, wc_ = 1 / s_u ** 2, 1 / s_c ** 2
t_d = (wu * t_u + wc_ * t_c) / (wu + wc_); s_d = 1 / math.sqrt(wu + wc_)
uP = [tu["PartIII|%s" % c[0]]["t"] for c in UCH]
t_uP = float(np.mean(uP))
t_dP = (wu * t_uP + wc_ * t_c) / (wu + wc_)
comb = dict(t_u=t_u, s_u_stat=s_u_stat, s_u_choice=s_u_ch, s_u=s_u, t_c=t_c, s_c_stat=s_c_stat, s_c_choice=s_c_ch, s_c=s_c, t_data=t_d, s_data=s_d,
            t_u_PartIII=t_uP, t_data_PartIII=t_dP)
RES["t_comb"] = comb
L("   u fit: mean of the five choices t = %.2f +- %.2f (stat + MC of -45..+60 pol2) +- %.2f (spread of the choices) = +- %.2f" % (t_u, s_u_stat, s_u_ch, s_u))
L("   no-bump region: mean of the four choices t = %.2f +- %.2f (stat, scaled) +- %.2f (choices and core width r = 0.893) = +- %.2f" % (t_c, s_c_stat, s_c_ch, s_c))
L("   weighted mean: t_data = %.2f +- %.2f   (u fit with the Part III composition: %.2f -> combined %.2f)" % (t_d, s_d, t_uP, t_dP))
qed = RD["toy"]["fit 0.30-0.40|all"]["t_eff"]
m4c = (gS == 2) & inR & (w > 0); m5c = (gS == 3) & inR & (w > 0)
n4t = (w * flag)[m4c].sum(); n5t = (w * flag)[m5c].sum()
t_qed_comb = (n4t * 1.0 + n5t * qed) / (n4t + n5t)
RES["t_qed"] = dict(t_C5=qed, t_combined=t_qed_comb)
L("   QED expectation (SMASH material + the radiative tail of C5 that SMASH lacks): t_C5 = %.2f, C4 1 -> combined t = %.2f" % (qed, t_qed_comb))

# ================================================================================================ B. propagation into f(u-edge)
L("\nB. f_eta(u-edge) (edgeu estimator, a_P fixed unless stated, h_sim = the same linear estimator on C2 + C3)")
cases = {}
def case(tag, fixed=None, free=None, text="", pu_fn=None):
    R, f = est(fixed=fixed, free=free)
    # replica noise of the difference to the nominal (peak template replicas; the data part cancels for fixed terms)
    sd = None
    if pu_fn is not None:
        dd = []
        for r in range(100):
            pr = PUR[r]
            fx_v = pu_fn(pr); fx_n = peak_term(1, 1, pu=pr)
            hv = R0.h - C.h_lin(R0, fx_v - peak_term(1, 1)); hn = R0.h - C.h_lin(R0, fx_n - peak_term(1, 1))
            dd.append((hv - hn) / C.h_lin(R0, ysim))
        sd = C.sd(dd)
    cases[tag] = dict(f=f, rel=f / f0 - 1, h=R.h, chi2=R.chi2, ndf=R.ndf, mc_of_shift=sd, text=text)
    L("   %-16s f = %.4f (%+6.2f %%)%s  chi2 %.1f/%d  %s" % (tag, f, 100 * (f / f0 - 1), "" if sd is None else " +- %.4f (MC of the shift)" % sd, R.chi2, R.ndf, text))
case("nominal", peak_term(1, 1), text="SMASH tail (t = 1)")
case("t_data", peak_term(t_d, t_d), text="t = t_data = %.2f" % t_d, pu_fn=lambda pr: peak_term(t_d, t_d, pu=pr))
case("t_data+1s", peak_term(t_d + s_d, t_d + s_d), text="t = %.2f" % (t_d + s_d), pu_fn=lambda pr: peak_term(t_d + s_d, t_d + s_d, pu=pr))
case("t_data-1s", peak_term(t_d - s_d, t_d - s_d), text="t = %.2f" % (t_d - s_d), pu_fn=lambda pr: peak_term(t_d - s_d, t_d - s_d, pu=pr))
case("t0", peak_term(0, 0), text="no tail (t = 0, SMASH core)", pu_fn=lambda pr: peak_term(0, 0, pu=pr))
case("t2", peak_term(2, 2), text="doubled tail (t = 2)", pu_fn=lambda pr: peak_term(2, 2, pu=pr))
case("qed", peak_term(1, qed), text="QED radiative tail of C5 (t_C5 = %.2f, rows of the toy: scale form)" % qed, pu_fn=lambda pr: peak_term(1, qed, pu=pr))
# QED toy rows directly (shape of the toy, not only the scale)
RR = np.load(P.OUT + "peaktail_rad_rows.npz")
br = g.bin(RR["u"]); ppr = C.ppbin_of(RR["mpp"].astype(np.float32))
mr = (br >= 0) & (ppr >= pa) & (ppr < pb)
K5rad = np.bincount(br[mr], weights=RR["w"][mr], minlength=g.nb)
K4n = utm(w, (2,))
case("qed_rows", C.AP * (K4n + K5rad), text="QED toy rows of C5 in u (C4 as SMASH)")
case("comp", peak_term(1, 1, AP4, AP5), text="Part III composition 0.237 C4 + 0.497 C5, SMASH tail")
case("comp_t0", peak_term(0, 0, AP4, AP5), text="Part III composition, no tail")
case("comp_tdP", peak_term(t_dP, t_dP, AP4, AP5), text="Part III composition, t = %.2f" % t_dP)
case("comp_qed", AP4 * K4n + AP5 * K5rad, text="Part III composition, QED toy rows of C5")
# the checker's fit and its fair version (from peaktail_ufit.json)
cases["checker"] = dict(f=UF["checker"]["f"], rel=UF["checker"]["rel"], text="the checker's Gaussian (no tail, amplitude free)")
cases["checker_fair"] = dict(f=UF["decomp"]["chk_fairSMASH"]["f"], rel=UF["decomp"]["chk_fairSMASH"]["rel"], text="the checker's estimator with the SMASH peak in h_sim")
L("   checker          f = %.4f (%+6.2f %%)  the checker's Gaussian, amplitude free (no tail)" % (cases["checker"]["f"], 100 * cases["checker"]["rel"]))
L("   checker_fair     f = %.4f (%+6.2f %%)  the same estimator, SMASH bump + SMASH peak (with its tail) in h_sim" % (cases["checker_fair"]["f"], 100 * cases["checker_fair"]["rel"]))
RES["f_u"] = cases
slope = (cases["t2"]["f"] - cases["t0"]["f"]) / 2 / f0
RES["f_u_slope_rel_per_t"] = slope
L("   slope: d ln f / dt = %+.4f per unit t (linear: f(t) = f0 [1 + %.4f (t - 1)])" % (slope, slope))
syst_u = abs(slope) * s_d
RES["syst_u_rel"] = syst_u
L("   proposed systematic of f_eta(u-edge) from the tail: |slope| x sigma(t_data) = %.4f x %.2f = %.2f %%" % (abs(slope), s_d, 100 * syst_u))

# ================================================================================================ C. bedge arrays
L("\nC. BEDGE ARRAYS (C4, C5 of the nominal bedge configuration in the bedge binning)")
NM4, NPP, NCLS = 240, 40, 3
eD = np.round((0.27 + 0.005 * np.arange(NPP + 1)) * 1000.) / 1000.
def bedge_bins(m4x, mppx, meex):
    ib = np.where((m4x >= 0.35) & (m4x < 0.95), np.floor((m4x - 0.35) / 0.0025), -1).astype(np.int64)
    p = np.searchsorted(eD, mppx, side="right") - 1
    ok = (ib >= 0) & (ib < NM4) & (mppx >= eD[0]) & (mppx < eD[-1]) & (p >= 0) & (p < NPP)
    e = np.where(meex < 0.010, 0, np.where(meex >= 0.010, 1, 2))
    return ib, p, e, ok
wb = S["w"] * S["f108"]                                     # bedge nominal: w x V3, no C2 factor (applied in Quant)
own = ((S["bits"] & 1) != 0) & ((S["bits"] & 8) != 0)
m4b = np.where(S["m4t"] > 0, C.S_NOM * (S["m4t"].astype(np.float64) + C.R_NOM * (S["m4"].astype(np.float64) - S["m4t"].astype(np.float64))), C.S_NOM * S["m4"].astype(np.float64))
mppb = 0.9952 * S["mpipi"].astype(np.float64)
ib, pbb, eb, okb = bedge_bins(m4b, mppb, S["mee"].astype(np.float64))
def arr(wt, sel):
    A_ = np.zeros((2, NCLS, NPP, NM4))
    for k, gg in enumerate((2, 3)):
        m = sel & okb & (gS == gg) & (wt != 0)
        np.add.at(A_[k], (eb[m], pbb[m], ib[m]), wt[m])
    return A_
base = own & (wb != 0)
Acore = arr(wb * (1 - flag), base); Atail = arr(wb * flag, base)
# the radiation toy for C5 (rows, weights already w x V3 x selection; no C2 factor needed for C5)
ib_r, p_r, e_r, ok_r = bedge_bins(RR["m4p"], 0.9952 * S["mpipi"][RR["rows"]].astype(np.float64), RR["mee"])
Arad = np.zeros((2, NCLS, NPP, NM4))
Arad[0] = Acore[0] + Atail[0]
m = ok_r
np.add.at(Arad[1], (e_r[m], p_r[m], ib_r[m]), RR["w"][m])
def warr(fn, a):
    v = a.ravel().astype(np.float64)
    with open(P.OUT + fn, "wb") as fp:
        fp.write(struct.pack("<Q", v.size)); fp.write(v.tobytes())
warr("bedge_C45_core.bin", Acore); warr("bedge_C45_tail.bin", Atail); warr("bedge_C45_rad.bin", Arad)
# check against the bedge sim array (nominal configuration)
SS = np.fromfile("/hdd3/przygoda/HADES/FAT_out/all_report/bedge/sim_S.bin", dtype=np.float64, offset=8)
NSC, NG = 15, 9
SS = SS.reshape(NSC, NG, NCLS, NPP, NM4)
chk = []
for k, gg in enumerate((2, 3)):
    ref = SS[0, gg]; mine = Acore[k] + Atail[k]
    chk.append((float(np.max(np.abs(ref - mine))), float(ref.sum()), float(mine.sum())))
    L("   group %s: max |bedge sim_S - core - tail| = %.3e (sums %.4f vs %.4f)" % ("C4" if gg == 2 else "C5", chk[-1][0], chk[-1][1], chk[-1][2]))
# per-class normalisation over 0.30-0.40 (pp bins 6..26)
pa_b, pb_b = int(round((0.30 - 0.27) / 0.005)), int(round((0.40 - 0.27) / 0.005))
norm = []
for k in range(2):
    N = (Acore[k] + Atail[k])[:, pa_b:pb_b].sum(); Nt = Atail[k][:, pa_b:pb_b].sum()
    norm.append((N, Nt))
with open(P.OUT + "bedge_C45_check.txt", "w") as fo:
    fo.write("# group maxdiff sum_bedge sum_mine N(0.30-0.40) Ntail(0.30-0.40)\n")
    for k, gg in enumerate(("C4", "C5")):
        fo.write("%s %.6e %.6f %.6f %.6f %.6f\n" % (gg, chk[k][0], chk[k][1], chk[k][2], norm[k][0], norm[k][1]))
    fo.write("t_data %.6f s_data %.6f t_qed_C5 %.6f\n" % (t_d, s_d, qed))
RES["bedge_arrays"] = dict(check=chk, norm=norm)
json.dump(RES, open(P.OUT + "peaktail_prop.json", "w"), indent=1)
L("\nwritten: peaktail_prop.json, peaktail_prop.log, bedge_C45_{core,tail,rad}.bin, bedge_C45_check.txt")
