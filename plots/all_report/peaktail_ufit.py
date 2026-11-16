#!/usr/bin/env python3
# ============================================================================
# peaktail_ufit.py  --  package peaktail, tasks 2 (the u > 0 region) and 3 (propagation into f_eta(u-edge)):
#                       the edgeu fit (edgeu_core.fit_step, the nominal estimator of edgeu_ana.py) with peak terms of
#                       different low-mass tails
# ============================================================================
# Written 2026-11-15 (prefix peaktail). Uses edgeu_core (data columns, u, CB rule, step template, the fit) unchanged and the
# simulated rows of peaktail_dump.C for the peak template (identical rows and weights as edgeu_dump.C in 0.30-0.40; the
# nominal template K is reproduced to 1e-9 first).
# Tail of a simulated peak row: its leptons' momentum deficit shifts M4 by more than 5 MeV/c2 (linearised, dl < -5 MeV;
# peaktail_sim.py: the rest is symmetric). Template with the tail scaled by t, rows kept:
#     K(t) = K_core (N - t N_tail)/N_core + t K_tail.
# Steps:
#   0. reproduce edgeu (h 1305.2, h_sim 3107.3, f 0.4200) and the checker's Gaussian peak term (h 1217.3, h_sim 3171.0, f 0.3839)
#   1. decompose the checker's -8.6 %: shape, free amplitude, width, per-row peak positions, the same treatment of SMASH
#   2. f(t) for t = 0 ... 3 (a_P fixed; and a_P(t) re-derived from the side bands)
#   3. t free in the u fit (the data constraint from u > 0); with a_P fixed and free
#   4. noise: simulation replicas (WP0 generator) of the variant difference; data bootstrap of the free-t fit
# Outputs: /hdd3/.../peaktail/peaktail_ufit.json, peaktail_ufit.log, peaktail_ufitfig.npz
# Run: python3 -I plots/all_report/peaktail_ufit.py
# ============================================================================
import json
import math
import sys

import numpy as np
from scipy.special import erf

sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import edgeu_core as C
import peaktail_core as P

LOG = open(P.OUT + "peaktail_ufit.log", "w")
def L(*a):
    s = " ".join(str(x) for x in a)
    print(s); LOG.write(s + "\n"); LOG.flush()
RES = {}
A = json.load(open(C.OUT + "edgeu_ana.json"))
FD = dict(np.load(C.OUT + "edgeu_figdata.npz"))
g = C.GRID_NOM
SIG = A["reso"]["sigma_core"]
T = C.step_gauss(g, SIG)
# ------------------------------------------------------------------ data (edgeu columns, as edgeu_ana.py / chk_hsim.py)
D = C.load_data(); u = C.u_data(D); S0 = (D["bits"] & 1) != 0
H = C.hist_data(D, u, S0, g); SL = C.slices_of(0.30, 0.40); sp = C.data_spec(H, SL, "all", g)
ysim = FD["sim_C2"] + FD["sim_C3"]          # bump C2 (gamma* corrected) + C3, expected raw rows
# ------------------------------------------------------------------ simulated peak rows (peaktail dump)
S = P.load_sim(); w = P.sim_weight(S); gS = S["g"]
us, mpp, m4p = C.u_sim(S)
M4r, J = P.m4_and_jacobian(S)
dl = C.S_NOM * (J[0] * (S["ep_p"] - S["ep_ptrue"]) + J[1] * (S["em_p"] - S["em_ptrue"])).astype(np.float64) * 1e-3
TAILX = 0.005
flag = (dl < -TAILX).astype(np.float64)
pa, pb = C.PB(0.30), C.PB(0.40)
def utmpl(wt, groups=(2, 3), extra=None, r=None, s=None):
    """sum over 0.30-0.40 of the simulated rows of the groups in u (fine grid), weights wt"""
    uu = us if (r is None and s is None) else C.u_sim(S, s=(C.S_NOM if s is None else s), r=(C.R_NOM if r is None else r))[0]
    b = g.bin(uu); ppb = C.ppbin_of(mpp.astype(np.float32))
    m = np.isin(gS, groups) & (b >= 0) & (ppb >= pa) & (ppb < pb) & (wt != 0)
    if extra is not None:
        m &= extra
    return np.bincount(b[m], weights=wt[m], minlength=g.nb).astype(np.float64)
pk = (gS == 2) | (gS == 3)
inR = (mpp >= 0.30) & (mpp < 0.40)
K1 = utmpl(w)
L("0. peak template from the peaktail dump vs edgeu (figdata peakT): max |diff| = %.2e rows (sum %.3f vs %.3f)" %
  (np.max(np.abs(K1 - FD["peakT"])), K1.sum(), FD["peakT"].sum()))
RES["repro_template_maxdiff"] = float(np.max(np.abs(K1 - FD["peakT"])))
def wt_t(t, base=None, fl=None):
    """row weights of the peak classes with the tail scaled by t; normalisation per class over the u-fit M(pipi) range"""
    base = w if base is None else base
    fl = flag if fl is None else fl
    out = base.copy()
    for gg in (2, 3):
        m = (gS == gg) & inR & (base > 0)
        N = base[m].sum(); Nt = (base * fl)[m].sum(); Nc = N - Nt
        c = (N - t * Nt) / Nc
        out[m] = base[m] * ((1 - fl[m]) * c + t * fl[m])
    return out
def Kt(t, groups=(2, 3), **kw):
    return utmpl(wt_t(t), groups, **kw)
for gg, nm in ((2, "C4"), (3, "C5")):
    m = (gS == gg) & inR & (w > 0)
    L("   %s: tail share (dl < -5 MeV) %.4f of %.1f expected raw rows (0.30-0.40)" % (nm, (w * flag)[m].sum() / w[m].sum(), w[m].sum()))
    RES["tailshare_%s" % nm] = float((w * flag)[m].sum() / w[m].sum())

# ------------------------------------------------------------------ the estimator
def est(fixed=None, free=None, simextra=None, lo=-0.045, hi=0.045, order=2, spec=None):
    """edgeu nominal estimator: fit N - CB - fixed = pol2 + [free] + h T; h_sim = the same linear estimator on the bump
    (+ simextra: an extra simulated spectrum added to the bump, e.g. the SMASH peak for 'the same treatment of SMASH')"""
    R = C.fit_step(sp if spec is None else spec, T, g, lo, hi, order, fixed=fixed, free_cols=free)
    ys = ysim if simextra is None else ysim + simextra
    hs = C.h_lin(R, ys)
    return R, hs, R.h / hs
R0, hs0, f0 = est(fixed=C.AP * K1)
L("   nominal: h %.2f +- %.2f, chi2 %.2f/%d, h_sim %.2f, f %.5f   [edgeu 1305.2 +- 62.3, 27.3/32, 3107.3, 0.42004]" % (R0.h, R0.eh, R0.chi2, R0.ndf, hs0, f0))
RES["nominal"] = dict(h=R0.h, eh=R0.eh, chi2=R0.chi2, ndf=R0.ndf, hs=hs0, f=f0)
# projections of the template parts on the step (rows of h per unit of a_P), nominal estimator
Kc = utmpl(w * (1 - flag)); Ktl = utmpl(w * flag)
K4c = utmpl(w * (1 - flag), (2,)); K4t = utmpl(w * flag, (2,)); K5c = utmpl(w * (1 - flag), (3,)); K5t = utmpl(w * flag, (3,))
prj = dict(K=C.h_lin(R0, K1), core=C.h_lin(R0, Kc), tail=C.h_lin(R0, Ktl), C4core=C.h_lin(R0, K4c), C4tail=C.h_lin(R0, K4t),
           C5core=C.h_lin(R0, K5c), C5tail=C.h_lin(R0, K5t), Ksum=float(K1.sum()), tailsum=float(Ktl.sum()))
RES["projections"] = prj
L("   projections on h (rows of step per unit a_P, nominal estimator): K %.2f = core %.2f + tail %.2f; C4 core %.2f tail %.2f; C5 core %.2f tail %.2f"
  % (prj["K"], prj["core"], prj["tail"], prj["C4core"], prj["C4tail"], prj["C5core"], prj["C5tail"]))
L("   (rows of the template in the fine grid: K %.1f, tail %.1f)" % (prj["Ksum"], prj["tailsum"]))

# ------------------------------------------------------------------ the checker's Gaussian peak term (chk_hsim.py / edgeu_partII_datastep.py)
m4d = D["m4"].astype(np.float64); pbd = C.ppbin_of(D["mpipi"])
base = S0 & (D["k"] == 0) & (pbd >= C.PB(0.30)) & (pbd < C.PB(0.40)) & (D["mee"] < C.MPI0)
inP = base & (m4d >= 0.53364) & (m4d < 0.55416)
inSB = base & (((m4d >= 0.5225) & (m4d < 0.53364)) | ((m4d >= 0.55416) & (m4d < 0.56530)))
rsb = 20.52 / 22.28
Ed = C.edge(D["mpipi"], D["mee"])
def gauss_tmpl(cen, wr, sig, norm=True):
    K = np.array([(wr * 0.5 * (erf((g.edges[j + 1] - cen) / (math.sqrt(2) * sig)) - erf((g.edges[j] - cen) / (math.sqrt(2) * sig)))).sum() for j in range(g.nb)])
    return K / K.sum() if norm else K
wr_d = np.concatenate([np.ones(inP.sum()), -rsb * np.ones(inSB.sum())])
cen_d = P.MPEAK - np.concatenate([Ed[inP], Ed[inSB]])
KG = gauss_tmpl(cen_d, wr_d, P.SPEAK_D)
Rk, hsk, fk = est(free=[KG])
L("\n1. THE CHECKER'S FIT: data Gaussian peak term (5.4 MeV, per-row centres 0.5439 - E of the data P-window rows minus side bands), amplitude free:")
L("   h %.2f +- %.2f, chi2 %.2f/%d, amplitude %.0f +- %.0f, h_sim %.2f, f %.5f (%+.2f %%)   [checker 1217.25, 3170.95, 0.38388]" %
  (Rk.h, Rk.eh, Rk.chi2, Rk.ndf, Rk.p[3], math.sqrt(Rk.C[3, 3]), hsk, fk, 100 * (fk / f0 - 1)))
RES["checker"] = dict(h=Rk.h, eh=Rk.eh, chi2=Rk.chi2, ndf=Rk.ndf, amp=float(Rk.p[3]), eamp=float(math.sqrt(Rk.C[3, 3])), hs=hsk, f=fk, rel=fk / f0 - 1)
# decomposition
dec = {}
def rec(tag, R, hs, f, text, extra=None):
    dec[tag] = dict(h=R.h, eh=R.eh, chi2=R.chi2, ndf=R.ndf, hs=hs, f=f, rel=f / f0 - 1, text=text)
    if extra:
        dec[tag].update(extra)
    L("   %-10s h %7.1f +- %5.1f chi2 %5.1f/%d h_sim %7.1f f %.4f (%+5.1f %%)  %s" % (tag, R.h, R.eh, R.chi2, R.ndf, hs, f, 100 * (f / f0 - 1), text))
L("   decomposition (each line: h_data, chi2, h_sim of the same linear estimator, f, change vs nominal):")
R, hs, f = est(fixed=C.AP * K1); rec("nom", R, hs, f, "SMASH peak (C4 + C5), a_P fixed [edgeu nominal]")
R, hs, f = est(free=[K1]); rec("smashFree", R, hs, f, "SMASH peak shape, amplitude free (a_P = %.3f)" % (R.p[3]), dict(aP=float(R.p[3])))
R, hs, f = est(); rec("noPeak", R, hs, f, "no peak term (Part II data-only fit)")
Kt0 = Kt(0.0)
R, hs, f = est(fixed=C.AP * Kt0); rec("core_fix", R, hs, f, "SMASH peak without its tail (t = 0, rows kept), a_P fixed")
R, hs, f = est(free=[Kt0]); rec("core_free", R, hs, f, "SMASH peak without its tail, amplitude free (%.3f)" % R.p[3], dict(aP=float(R.p[3])))
# SMASH peak rows as plain Gaussians at their own centres (0.5439 - E of the simulated row, data scale), 5.4 MeV
mS = pk & inR & (w > 0) & np.isfinite(us)
cen_s = P.MPEAK - C.edge(mpp[mS], S["mee"][mS].astype(np.float64))
KGs = gauss_tmpl(cen_s, w[mS], P.SPEAK_D)
R, hs, f = est(free=[KGs]); rec("gaussSMASHpos", R, hs, f, "Gaussian 5.4 MeV at the SMASH peak rows' centres, amplitude free")
KGs61 = gauss_tmpl(cen_s, w[mS], P.SPEAK_S)
R, hs, f = est(free=[KGs61]); rec("gaussSMASH61", R, hs, f, "Gaussian 6.09 MeV (SMASH width) at the SMASH centres, amplitude free")
KG61 = gauss_tmpl(cen_d, wr_d, P.SPEAK_S)
R, hs, f = est(free=[KG61]); rec("gaussData61", R, hs, f, "data centres, Gaussian 6.09 MeV, amplitude free")
# amplitude of the data Gaussian fixed at the expected peak (a_P x SMASH K in the fine grid)
Aexp = C.AP * K1.sum()
R, hs, f = est(fixed=Aexp * KG); rec("gauss_fixA", R, hs, f, "the checker's Gaussian with the amplitude fixed at a_P x SMASH (%.0f rows)" % Aexp)
R, hs, f = est(fixed=Aexp * Kt0 / Kt0.sum()); rec("core_fixA", R, hs, f, "SMASH core shape at the same fixed amplitude")
# the same treatment of SMASH: the checker's estimator applied to SMASH bump + SMASH peak x a_P / f (the simulated peak
# with its tail, scaled as in the data)
R = C.fit_step(sp, T, g, free_cols=[KG]); hsf = C.h_lin(R, ysim + (C.AP / 0.42) * K1)
rec("chk_fairSMASH", R, hsf, R.h / hsf, "the checker's estimator, h_sim from SMASH bump + SMASH peak (a_P/f x C4+C5, with the SMASH tail)")
hsf2 = C.h_lin(R, ysim + (C.AP / 0.42) * Kt0)
rec("chk_fairCore", R, hsf2, R.h / hsf2, "the same, SMASH peak without its tail")
# the data template projected: L_h of the checker's estimator on the SMASH tail
RES["decomp"] = dec
prjk = dict(K=C.h_lin(Rk, C.AP * K1), core=C.h_lin(Rk, C.AP * Kc), tail=C.h_lin(Rk, C.AP * Ktl))
RES["checker_projections"] = prjk
L("   projections in the checker's estimator (rows, x a_P): SMASH peak %.1f = core %.1f + tail %.1f" % (prjk["K"], prjk["core"], prjk["tail"]))

# ------------------------------------------------------------------ 2. f(t)
L("\n2. f(t): the SMASH peak term with the tail scaled by t (rows kept), a_P fixed 0.418249 (nominal estimator)")
TS = [0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0]
ft = []
for t in TS:
    Kx = Kt(t)
    R, hs, f = est(fixed=C.AP * Kx)
    ft.append(dict(t=t, h=R.h, chi2=R.chi2, hs=hs, f=f, rel=f / f0 - 1))
    L("   t = %.2f: h %7.2f chi2 %5.2f f %.5f (%+6.2f %%)" % (t, R.h, R.chi2, f, 100 * (f / f0 - 1)))
RES["f_of_t"] = ft
# per class
L("   per class (the other class at t = 1):")
ftc = {}
for gg, nm in ((2, "C4"), (3, "C5")):
    ftc[nm] = []
    for t in (0.0, 0.5, 1.5, 2.0, 3.0):
        wx = w.copy(); wt = wt_t(t)
        wx[gS == gg] = wt[gS == gg]
        R, hs, f = est(fixed=C.AP * utmpl(wx))
        ftc[nm].append(dict(t=t, f=f, rel=f / f0 - 1))
        L("     %s t = %.1f: f %.5f (%+5.2f %%)" % (nm, t, f, 100 * (f / f0 - 1)))
RES["f_of_t_class"] = ftc
# the same with the free SMASH-shape amplitude (U8 aPfree analogue)
L("   with the peak amplitude free:")
ftf = []
for t in (0.0, 0.5, 1.0, 1.5, 2.0, 3.0):
    R, hs, f = est(free=[Kt(t)])
    ftf.append(dict(t=t, h=R.h, hs=hs, f=f, aP=float(R.p[3]), rel=f / f0 - 1, chi2=R.chi2))
    L("   t = %.1f: h %7.2f chi2 %5.2f a_P %.3f f %.5f (%+6.2f %%)" % (t, R.h, R.chi2, R.p[3], f, 100 * (f / f0 - 1)))
RES["f_of_t_free"] = ftf

# ------------------------------------------------------------------ 3. t free in the u fit
L("\n3. THE DATA CONSTRAINT FROM THE u FIT: t free (K(t) is linear in t: fixed a_P K_c N/N_c + t a_P (K_t - K_c N_t/N_c))")
def lin_parts(groups=(2, 3)):
    """K(t) = K0 + t K1 summed over the classes (normalisation per class)"""
    K0 = np.zeros(g.nb); KK = np.zeros(g.nb)
    for gg in groups:
        m = (gS == gg) & inR & (w > 0)
        N = w[m].sum(); Nt = (w * flag)[m].sum(); Nc = N - Nt
        kc = utmpl(w * (1 - flag), (gg,)); kt = utmpl(w * flag, (gg,))
        K0 += kc * N / Nc; KK += kt - kc * Nt / Nc
    return K0, KK
K0, KK = lin_parts()
tfree = {}
for lab, lo, hi, order in (("-45..+45 pol2", -0.045, 0.045, 2), ("-45..+55 pol2", -0.045, 0.055, 2), ("-45..+60 pol2", -0.045, 0.060, 2),
                           ("-45..+70 pol3", -0.045, 0.070, 3), ("-55..+70 pol3", -0.055, 0.070, 3), ("-45..+45 pol1", -0.045, 0.045, 1)):
    R = C.fit_step(sp, T, g, lo, hi, order, fixed=C.AP * K0, free_cols=[C.AP * KK])
    it = order + 1
    tt, et = float(R.p[it]), float(math.sqrt(R.C[it, it]))
    rho = float(R.C[it, -1] / math.sqrt(R.C[it, it] * R.C[-1, -1]))
    hs = C.h_lin(R, ysim)
    tfree[lab] = dict(t=tt, et=et, rho_th=rho, h=R.h, eh=R.eh, f=R.h / hs, chi2=R.chi2, ndf=R.ndf)
    L("   %-14s t = %6.2f +- %5.2f (corr. with h %+.2f), h %7.1f +- %5.1f, f %.4f, chi2 %.1f/%d" % (lab, tt, et, rho, R.h, R.eh, R.h / hs, R.chi2, R.ndf))
# with the peak amplitude free as well (a = a_P x amplitude scale, t free): columns K0 and KK free
for lab, lo, hi, order in (("-45..+60 pol2, amp. free", -0.045, 0.060, 2), ("-45..+70 pol3, amp. free", -0.045, 0.070, 3)):
    R = C.fit_step(sp, T, g, lo, hi, order, free_cols=[K0, KK])
    a0, a1 = float(R.p[order + 1]), float(R.p[order + 2])
    Cm = R.C[order + 1:order + 3, order + 1:order + 3]
    tt = a1 / a0
    et = abs(tt) * math.sqrt(Cm[1, 1] / a1 ** 2 + Cm[0, 0] / a0 ** 2 - 2 * Cm[0, 1] / (a0 * a1)) if a1 != 0 else float("nan")
    hs = C.h_lin(R, ysim)
    tfree[lab] = dict(t=tt, et=et, amp=a0, eamp=float(math.sqrt(Cm[0, 0])), h=R.h, eh=R.eh, f=R.h / hs, chi2=R.chi2, ndf=R.ndf)
    L("   %-26s amplitude %.3f +- %.3f, t = %6.2f +- %5.2f, h %7.1f +- %5.1f, f %.4f, chi2 %.1f/%d" % (lab, a0, math.sqrt(Cm[0, 0]), tt, et, R.h, R.eh, R.h / hs, R.chi2, R.ndf))
RES["t_free"] = tfree

# ------------------------------------------------------------------ 4. noise of the variant differences
L("\n4. NOISE of Delta f(t) (WP0 simulation replicas of the peak rows, 100; the data part is fixed for a fixed peak term)")
rows = np.where(pk & inR & (w > 0))[0]
RW = C.sim_replica_weights(S, rows, 100)
bb = g.bin(us[rows]); okb = bb >= 0
dfr = {t: [] for t in (0.0, 0.5, 2.0, 3.0)}
pr_rep = []
for r in range(100):
    wr = w.copy(); wr[rows] = w[rows] * RW[r]
    def Kr(t):
        out = wr.copy()
        for gg in (2, 3):
            m = (gS == gg) & inR & (wr > 0)
            N = wr[m].sum(); Nt = (wr * flag)[m].sum(); Nc = N - Nt
            c = (N - t * Nt) / Nc
            out[m] = wr[m] * ((1 - flag[m]) * c + t * flag[m])
        return utmpl(out)
    K1r = Kr(1.0)
    h1 = R0.h - C.h_lin(R0, C.AP * (K1r - K1))   # nominal with the replica template (linear in the fixed term)
    pr_rep.append(C.h_lin(R0, C.AP * K1r))
    for t in dfr:
        ht = R0.h - C.h_lin(R0, C.AP * (Kr(t) - K1))
        dfr[t].append((ht - h1) / hs0)
noise = {str(t): dict(sd=C.sd(v), mean=float(np.mean(v))) for t, v in dfr.items()}
RES["noise_rep"] = noise
RES["peak_projection_rep"] = dict(sd=C.sd(pr_rep), mean=float(np.mean(pr_rep)))
for t, v in dfr.items():
    L("   t = %.1f: Delta f replica mean %+.5f, sd %.5f (relative sd %.2f %%)" % (t, np.mean(v), C.sd(v), 100 * C.sd(v) / f0))
L("   projection of a_P K on h over the replicas: mean %.2f, sd %.2f rows" % (np.mean(pr_rep), C.sd(pr_rep)))
# bootstrap of the free-t fit (data noise) for the -45..+60 pol2 choice
rowsd = np.where(S0 & (g.bin(u) >= 0) & (C.ppbin_of(D["mpipi"]) >= 0))[0]
BW = C.data_boot_weights(D, rowsd, 200)
Dsub = {k: D[k][rowsd] for k in ("k", "mpipi", "mee")}
HR = C.hist_data(Dsub, u[rowsd], np.ones(len(rowsd), bool), g, weights=BW)
tb = []; fkb = []; f0b = []
for b in range(HR.shape[0]):
    spb = C.data_spec(HR[b], SL, "all", g)
    Rb = C.fit_step(spb, T, g, -0.045, 0.060, 2, fixed=C.AP * K0, free_cols=[C.AP * KK]); tb.append(float(Rb.p[3]))
    Rk_b = C.fit_step(spb, T, g, free_cols=[KG]); fkb.append(Rk_b.h / C.h_lin(Rk_b, ysim))
    Rn_b = C.fit_step(spb, T, g, fixed=C.AP * K1); f0b.append(Rn_b.h / C.h_lin(Rn_b, ysim))
dk = np.array(fkb) - np.array(f0b)
RES["boot"] = dict(t_sd=C.sd(tb), t_mean=float(np.mean(tb)), checker_minus_nom_sd=C.sd(dk), checker_minus_nom_mean=float(np.mean(dk)))
L("   bootstrap (200): t (-45..+60 pol2) sd %.2f (mean %.2f); f(checker) - f(nominal): mean %+.4f, sd %.4f (data noise of the checker's shift)" %
  (C.sd(tb), np.mean(tb), np.mean(dk), C.sd(dk)))
json.dump(RES, open(P.OUT + "peaktail_ufit.json", "w"), indent=1)
np.savez_compressed(P.OUT + "peaktail_ufitfig.npz", x=g.x, sig=sp.sig, N=sp.N, cb=sp.cb, cbv=np.diag(sp.cbcov), K1=K1, Kc=Kc, Ktl=Ktl, KG=KG,
                    Kt0=Kt0, Kt2=Kt(2.0), K4c=K4c, K4t=K4t, K5c=K5c, K5t=K5t, ysim=ysim, R0_idx=R0.idx, R0_model=R0.model, R0_smooth=R0.smooth,
                    Rk_model=Rk.model, Rk_smooth=Rk.smooth, Rk_amp=Rk.p[3], T=T, R0_h=R0.h, Rk_h=Rk.h)
L("\nwritten: peaktail_ufit.json, peaktail_ufit.log, peaktail_ufitfig.npz")
