#!/usr/bin/env python3
# ============================================================================
# edgeu_partII_datastep.py  --  Part II (data only): the bump step in the distance to its kinematic edge, u
# ============================================================================
# Written 2026-11-15 (prefix edgeu_partII_) for report_exp.tex (Part II, experimental data alone).
# u = M4 - E(M(pi+pi-), M(e+e-)), E = the frozen kinematic edge of defs.h (alldefs::Edge; edgeu_core.edge): pure kinematics,
# no simulation needed. This macro repeats the step fit of package edgeu (Part III) WITHOUT any simulated input:
#   - no SMASH peak term a_P K(u) (the pol2 absorbs the rise into the eta peak; in Part III the peak term moves h by only
#     3.86 rows per unit of a_P),
#   - the step width either left free (data decide) or fixed at the template width of Part III (6.70 MeV/c2, from SMASH;
#     quoted only for the comparison).
# The Part III nominal (peak term from SMASH, width 6.70) is reproduced first from the stored peak template
# (edgeu_figdata.npz), as a check of the machinery; it is not a Part II number.
# Data: S0 rows of the edgeu column dump (cor_ok, isBest >= 0), 0.30 <= M(pipi) < 0.40, M(ee) < m_pi0; N - CB per 2.5 MeV/c2
# of u, CB = 2 sqrt(N++ N--) in 10 MeV/c2 cells per M(pipi) slice (0.30-0.34-0.37-0.39-0.40), as in edgeu.
# Fit: N - CB = pol2(u) + h T(u) on -45 <= u < +45 MeV/c2 (GLS, V = cov(CB) + diag(model), iterated; edgeu_core.fit_step).
# Like-sign rows N++ and N-- alone: pol2 + step (null tests), as in edgeu.
# Outputs: /hdd3/przygoda/HADES/FAT_out/all_report/edgeu_partII/edgeu_partII.json, edgeu_partII.log;
#          docs/latex/figs_exp/edgeu_data_step.{pdf,png}.
# Run (from /home/przygoda/HADES/FAT/FAT): python3 -I plots/all_report/edgeu_partII_datastep.py
# ============================================================================
import json
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import edgeu_core as C

OUT = "/hdd3/przygoda/HADES/FAT_out/all_report/edgeu_partII/"
FIGX = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_exp/"
os.makedirs(OUT, exist_ok=True)
LOG = open(OUT + "edgeu_partII.log", "w")
def P(*a):
    s = " ".join(str(x) for x in a)
    print(s); LOG.write(s + "\n"); LOG.flush()

A = json.load(open(C.OUT + "edgeu_ana.json"))
FD = dict(np.load(C.OUT + "edgeu_figdata.npz"))
SIG_U = A["reso"]["sigma_core"]          # Part III template width (SMASH), 6.70 MeV/c2
g = C.GRID_NOM
D = C.load_data()
u = C.u_data(D)
S0 = (D["bits"] & 1) != 0
H = C.hist_data(D, u, S0, g)
SL = C.slices_of(0.30, 0.40)
sp = C.data_spec(H, SL, "all", g)
RES = {}

def fit(spec, sigma, lo=-0.045, hi=0.045, order=2, fixed=None, free=None):
    return C.fit_step(spec, C.step_gauss(g, sigma), g, lo, hi, order, fixed=fixed, free_cols=free)
def summ(R):
    return dict(h=R.h, eh=R.eh, chi2=R.chi2, chi20=R.chi20, ndf=R.ndf, dchi2=R.chi20 - R.chi2,
                z=math.sqrt(max(R.chi20 - R.chi2, 0.0)))
def free_width(spec, lo=-0.045, hi=0.045, order=2, fixed=None, free=None):
    """scan the step width; parabola at the chi2 minimum; h at the minimum with the width error propagated through dh/dsigma"""
    grid = np.arange(0.0030, 0.01201, 0.0002)
    rs = [fit(spec, s, lo, hi, order, fixed, free) for s in grid]
    ch = np.array([r.chi2 for r in rs]); hh = np.array([r.h for r in rs])
    i = int(np.argmin(ch)); i = min(max(i, 2), len(grid) - 3)
    c2, c1, c0 = np.polyfit(grid[i - 2:i + 3], ch[i - 2:i + 3], 2)
    if c2 <= 0:
        return dict(sigma=float("nan"), esigma=float("nan"), h=float("nan"), eh=float("nan"), chi2_grid=ch.tolist()), None
    smin = -c1 / (2 * c2); es = 1.0 / math.sqrt(c2)
    Rm = fit(spec, smin, lo, hi, order, fixed, free)
    dh = (fit(spec, smin + 0.0002, lo, hi, order, fixed, free).h - fit(spec, smin - 0.0002, lo, hi, order, fixed, free).h) / 0.0004
    eh_tot = math.sqrt(Rm.eh ** 2 + (dh * es) ** 2)
    out = summ(Rm); out.update(sigma=smin, esigma=es, eh_cond=Rm.eh, eh=eh_tot, dh_dsigma=dh)
    return out, Rm

# ---------------------------------------------------------------- 0. reproduce the Part III nominal (needs the SMASH peak term)
R3 = fit(sp, SIG_U, fixed=C.AP * FD["peakT"])
RES["partIII_nominal_repro"] = summ(R3)
P("0. Part III nominal reproduced (SMASH peak term a_P K, width %.2f MeV): h = %.1f +- %.1f, chi2 %.1f/%d (no step %.1f) [edgeu: %.1f +- %.1f]"
  % (SIG_U * 1e3, R3.h, R3.eh, R3.chi2, R3.ndf, R3.chi20, A["nom"]["hd"], A["nom"]["ehd"]))

# ---------------------------------------------------------------- 1. data only
RA = fit(sp, SIG_U)
RES["data_fixedwidth"] = summ(RA)
P("1a. data only, no peak term, width fixed %.2f MeV (Part III template width): h = %.1f +- %.1f, chi2 %.1f/%d, dchi2 %.1f (%.1f sigma)"
  % (SIG_U * 1e3, RA.h, RA.eh, RA.chi2, RA.ndf, RA.chi20 - RA.chi2, math.sqrt(RA.chi20 - RA.chi2)))
FB, RB = free_width(sp)
RES["data_freewidth"] = FB
P("1b. data only, no peak term, width free: sigma = %.2f +- %.2f MeV, h = %.1f +- %.1f (conditional %.1f; dh/dsigma %.0f rows/MeV), chi2 %.1f/%d, dchi2 %.1f (%.1f sigma)"
  % (FB["sigma"] * 1e3, FB["esigma"] * 1e3, FB["h"], FB["eh"], FB["eh_cond"], FB["dh_dsigma"] * 1e-3, FB["chi2"], FB["ndf"] - 1, FB["dchi2"], FB["z"]))
FB3, _ = free_width(sp, fixed=C.AP * FD["peakT"])
RES["partIII_freewidth_repro"] = FB3
P("    (with the SMASH peak term, width free: sigma = %.2f +- %.2f MeV [edgeu 7.02 +- 0.79], h = %.1f +- %.1f)" % (FB3["sigma"] * 1e3, FB3["esigma"] * 1e3, FB3["h"], FB3["eh"]))
# range / order variants without any simulated input
P("1c. data only, width %.2f MeV, range and order variants:" % (SIG_U * 1e3))
RES["data_variants"] = {}
for lab, lo, hi, order in [("-45..+45 pol2", -0.045, 0.045, 2), ("-45..+30 pol2", -0.045, 0.030, 2), ("-45..+20 pol2", -0.045, 0.020, 2),
                           ("-35..+35 pol2", -0.035, 0.035, 2), ("-40..+40 pol2", -0.040, 0.040, 2), ("-50..+50 pol2", -0.050, 0.050, 2),
                           ("-55..+55 pol2", -0.055, 0.055, 2), ("-45..+45 pol3", -0.045, 0.045, 3)]:
    R = fit(sp, SIG_U, lo, hi, order)
    RES["data_variants"][lab] = summ(R)
    P("    %-14s h = %7.1f +- %5.1f  chi2 %5.1f/%d  p %.2f" % (lab, R.h, R.eh, R.chi2, R.ndf, 1 - __import__("scipy.stats", fromlist=["chi2"]).chi2.cdf(R.chi2, R.ndf)))

# ---------------------------------------------------------------- 1d. a peak term built from the data themselves
# The eta peak enters the fit range at u > +15 MeV/c2. Its shape in u is the frozen data Gaussian in M4 (0.5439, 5.4 MeV/c2,
# Part II) shifted by -E(M(pipi), M(ee)) of each peak row. The (M(pipi), M(ee)) of the peak rows come from the data: unlike-sign
# rows of the peak window 0.53364-0.55416 (weight +1) minus the side bands 0.5225-0.53364 and 0.55416-0.56530 (weight -r,
# r = 20.52/22.28 = 0.921 for a linear background; the CB, 7-10 % and smooth, cancels in the same way). Amplitude free.
from scipy.special import erf
MPK, SPK = 0.5439, 0.0054
m4d = D["m4"].astype(np.float64)
pb_all = C.ppbin_of(D["mpipi"])
inpp = (pb_all >= C.PB(0.30)) & (pb_all < C.PB(0.40)) & (D["mee"] < C.MPI0)
base = S0 & (D["k"] == 0) & inpp
inP = base & (m4d >= 0.53364) & (m4d < 0.55416)
inSB = base & (((m4d >= 0.5225) & (m4d < 0.53364)) | ((m4d >= 0.55416) & (m4d < 0.56530)))
RSB = (0.55416 - 0.53364) / ((0.53364 - 0.5225) + (0.56530 - 0.55416))
Erow = C.edge(D["mpipi"], D["mee"])
def peak_template(wrows, Ec):
    """sum over rows of w x bin integral of Gauss(u; MPK - E, SPK); normalised to unit area"""
    cen = MPK - Ec
    K = np.zeros(g.nb)
    for j in range(g.nb):
        a, b = g.edges[j], g.edges[j + 1]
        K[j] = (wrows * 0.5 * (erf((b - cen) / (math.sqrt(2) * SPK)) - erf((a - cen) / (math.sqrt(2) * SPK)))).sum()
    return K / K.sum()
wr = np.concatenate([np.ones(inP.sum()), -RSB * np.ones(inSB.sum())])
Ec = np.concatenate([Erow[inP], Erow[inSB]])
KD = peak_template(wr, Ec)
RES["peak_data"] = dict(rows_P=int(inP.sum()), rows_SB=int(inSB.sum()), r=RSB, net=float(wr.sum()))
P("1d. data peak template: %d peak-window rows - %.3f x %d side-band rows = %.0f net rows (0.30-0.40, S0)" % (inP.sum(), RSB, inSB.sum(), wr.sum()))
RK = fit(sp, SIG_U, free=[KD])
RES["data_peakdata"] = summ(RK); RES["data_peakdata"]["A"] = float(RK.p[3]); RES["data_peakdata"]["eA"] = float(math.sqrt(RK.C[3, 3]))
P("    width %.2f, data peak term free: h = %.1f +- %.1f, chi2 %.1f/%d, dchi2 %.1f (%.1f sigma); peak amplitude %.0f +- %.0f rows (net side-band count %.0f / 0.943)"
  % (SIG_U * 1e3, RK.h, RK.eh, RK.chi2, RK.ndf, RK.chi20 - RK.chi2, math.sqrt(RK.chi20 - RK.chi2), RK.p[3], math.sqrt(RK.C[3, 3]), wr.sum()))
FK, RKf = free_width(sp, free=[KD])
RES["data_peakdata_freewidth"] = FK
P("    width free, data peak term free: sigma = %.2f +- %.2f MeV, h = %.1f +- %.1f, chi2 %.1f/%d" % (FK["sigma"] * 1e3, FK["esigma"] * 1e3, FK["h"], FK["eh"], FK["chi2"], FK["ndf"] - 1))
RES["data_peakdata_variants"] = {}
P("    data peak term free, width %.2f, range and order variants:" % (SIG_U * 1e3))
from scipy.stats import chi2 as chi2dist
for lab, lo, hi, order in [("-45..+45 pol2", -0.045, 0.045, 2), ("-40..+40 pol2", -0.040, 0.040, 2), ("-35..+35 pol2", -0.035, 0.035, 2),
                           ("-50..+50 pol2", -0.050, 0.050, 2), ("-55..+55 pol2", -0.055, 0.055, 2), ("-45..+30 pol2", -0.045, 0.030, 2),
                           ("-45..+45 pol3", -0.045, 0.045, 3)]:
    R = fit(sp, SIG_U, lo, hi, order, free=[KD])
    RES["data_peakdata_variants"][lab] = summ(R)
    P("      %-14s h = %7.1f +- %5.1f  chi2 %5.1f/%d  p %.2f" % (lab, R.h, R.eh, R.chi2, R.ndf, 1 - chi2dist.cdf(R.chi2, R.ndf)))
# the same range variants with the SMASH peak term (Part III), for comparison
RES["partIII_variants"] = {}
P("    for comparison, SMASH peak term (Part III), width %.2f:" % (SIG_U * 1e3))
for lab, lo, hi, order in [("-45..+45 pol2", -0.045, 0.045, 2), ("-40..+40 pol2", -0.040, 0.040, 2), ("-35..+35 pol2", -0.035, 0.035, 2),
                           ("-50..+50 pol2", -0.050, 0.050, 2), ("-55..+55 pol2", -0.055, 0.055, 2), ("-45..+30 pol2", -0.045, 0.030, 2)]:
    R = fit(sp, SIG_U, lo, hi, order, fixed=C.AP * FD["peakT"])
    RES["partIII_variants"][lab] = summ(R)
    P("      %-14s h = %7.1f +- %5.1f  chi2 %5.1f/%d  p %.2f" % (lab, R.h, R.eh, R.chi2, R.ndf, 1 - chi2dist.cdf(R.chi2, R.ndf)))

# ---------------------------------------------------------------- 2. bootstrap of the data-only fit (200 replicas, event groups)
rows = np.where(S0 & (g.bin(u) >= 0) & (C.ppbin_of(D["mpipi"]) >= 0))[0]
BW = C.data_boot_weights(D, rows, 200)
Dsub = {k: D[k][rows] for k in ("k", "mpipi", "mee")}
HR = C.hist_data(Dsub, u[rows], np.ones(len(rows), bool), g, weights=BW)
hb = []
for b in range(HR.shape[0]):
    spb = C.data_spec(HR[b], SL, "all", g)
    hb.append(fit(spb, SIG_U).h)
RES["data_fixedwidth"]["boot_sd"] = C.sd(hb)
P("2. bootstrap (200 replicas) of 1a: sd(h) = %.1f, ratio to the fit error %.2f" % (C.sd(hb), C.sd(hb) / RA.eh))

# ---------------------------------------------------------------- 3. like-sign rows alone
Z = np.zeros(g.nb); ZZ = np.zeros((g.nb, g.nb))
RES["ls"] = {}
for kk, kl in ((1, "pp"), (2, "mm")):
    for cl in ("all", "E1", "E2"):
        cs = C.CLS_SETS[cl]
        y = sum(H[kk, a:b][:, cs].sum((0, 1)) for a, b in SL)
        R = fit(C.Spec(y, Z, ZZ), SIG_U)
        RES["ls"]["%s_%s" % (kl, cl)] = summ(R)
        RES["ls"]["%s_%s" % (kl, cl)]["rows"] = float(y[R.idx].sum())
        P("3. like-sign %s %-3s: h = %6.1f +- %5.1f (chi2 %.1f/%d; %.1f sigma), rows in range %.0f" % (kl, cl, R.h, R.eh, R.chi2, R.ndf, R.h / R.eh, y[R.idx].sum()))
hcb = C.h_lin(RA, sp.cb)
RES["cb_step"] = dict(h=hcb, rel=hcb / RA.h)
P("4. step carried by the subtracted CB (data-only estimator): %.1f rows (%.3f of h)" % (hcb, hcb / RA.h))
# control region without bump: M(pipi) 0.42-0.46, pseudo edge E(M(pipi) - 0.10, M(ee)), no peak term
uc = C.u_data(D, mshift=-0.10)
Hc = C.hist_data(D, uc, S0, g)
SLc = C.slices_of(0.42, 0.46, [0.42, 0.46])
spc = C.data_spec(Hc, SLc, "all", g)
Rc = fit(spc, SIG_U)
RES["control"] = summ(Rc)
P("5. control M(pipi) 0.42-0.46, pseudo edge, data only: h = %.1f +- %.1f (chi2 %.1f/%d) [Part III with the SMASH peak term: %.1f +- %.1f]"
  % (Rc.h, Rc.eh, Rc.chi2, Rc.ndf, A["chk"]["control"]["hd"], A["chk"]["control"]["ehd"]))
# rows near the edge (bump share): for the plain-words comparison of the like-sign step
near = (g.x >= -0.010) & (g.x < 0.0)
RES["near_edge_per_bin"] = dict(sig=float(sp.sig[near].mean()), mm=float(sp.mm[near].mean()), pp=float(sp.pp[near].mean()))
P("6. rows per 2.5 MeV just below the edge (-10..0 MeV): N - CB %.0f, N-- %.0f, N++ %.0f" % (sp.sig[near].mean(), sp.mm[near].mean(), sp.pp[near].mean()))
json.dump(RES, open(OUT + "edgeu_partII.json", "w"), indent=1)

# ---------------------------------------------------------------- figure (Part II colours: N black, CB red, N - CB blue, e+e+ green, e-e- orange)
plt.rcParams.update({"font.size": 15, "axes.labelsize": 16, "axes.titlesize": 15, "legend.fontsize": 12.5, "xtick.labelsize": 14,
                     "ytick.labelsize": 12, "pdf.fonttype": 42})
BLUE, RED, GREEN, ORANGE, FIT = "#1f3fbf", "#e8202a", "#1a9641", "#ff7f00", "#c51b8a"
x = g.x * 1e3
fig = plt.figure(figsize=(12.5, 10.5))
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0])
axs = [fig.add_subplot(gs[0, :]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])]
R = RA; idx = R.idx
ax = axs[0]
m = (x > -100) & (x < 100)
ax.errorbar(x[m], sp.sig[m], yerr=np.sqrt(np.maximum(sp.N[m], 1) + np.diag(sp.cbcov)[m]), fmt="o", ms=3.5, color=BLUE, label="$N-\\mathrm{CB}$")
ax.plot(x[m], sp.cb[m], "-", color=RED, lw=1.5, label="CB $=2\\sqrt{N_{++}N_{--}}$")
ax.plot(x[idx], R.model, "-", color=FIT, lw=2.2, label="fit: pol2 + step at $u=0$")
ax.plot(x[idx], R.smooth, "--", color="#555555", lw=1.6, label="its smooth part (pol2)")
ax.axvline(0, color="k", lw=0.8, ls=":")
ax.axvspan(-45, 45, color="#eeeeee", zorder=0)
ax.set_xlim(-100, 100); ax.set_ylim(0, 6300)
ax.set_xlabel("$u=M_4-E(M_{\\pi\\pi},M_{ee})$ [MeV/$c^2$]"); ax.set_ylabel("rows per 2.5 MeV/$c^2$")
ax.set_title("(a) data, $0.30\\leq M(\\pi^+\\pi^-)<0.40$ GeV/$c^2$")
ax.text(-97, 5850, "bump ($u<0$)", fontsize=14); ax.text(50, 5850, "$\\eta$ peak", fontsize=14)
ax.legend(loc="upper left", bbox_to_anchor=(0.0, 0.92), fontsize=12.5)
ax = axs[1]
resid = sp.sig[idx] - R.smooth
err = np.sqrt(np.diag(R.V))
ax.errorbar(x[idx], resid, yerr=err, fmt="o", ms=4, color=BLUE, label="data minus its smooth part")
T = C.step_gauss(g, SIG_U)
ax.plot(x[idx], R.h * T[idx], "-", color=FIT, lw=2.2, label="step $h=%.0f\\pm%.0f$ (%.0f$\\sigma$)" % (R.h, R.eh, math.sqrt(R.chi20 - R.chi2)))
ax.axhline(0, color="k", lw=0.8); ax.axvline(0, color="k", lw=0.8, ls=":")
ax.set_xlabel("$u$ [MeV/$c^2$]"); ax.set_ylabel("rows per 2.5 MeV/$c^2$")
ax.set_title("(b) the step at $u=0$ (data only)")
ax.text(-43, -840, "no peak term, width 6.7 MeV/$c^2$: $\\chi^2/\\mathrm{ndf}=%.1f/%d$\nwidth free: $h=%.0f\\pm%.0f$, $\\sigma=%.1f\\pm%.1f$ MeV/$c^2$\ndata Gaussian peak term: $h=%.0f\\pm%.0f$ ($\\chi^2$ %.1f/%d)" %
        (R.chi2, R.ndf, FB["h"], FB["eh"], FB["sigma"] * 1e3, FB["esigma"] * 1e3, RK.h, RK.eh, RK.chi2, RK.ndf), fontsize=11.5, va="bottom")
ax.set_ylim(-880, 1850)
ax.legend(loc="upper right", fontsize=12)
ax = axs[2]
for kk, kl, col, mk, lab in ((1, "pp", GREEN, "s", "$e^+e^+$"), (2, "mm", ORANGE, "^", "$e^-e^-$")):
    y = sum(H[kk, a:b].sum((0, 1)) for a, b in SL)
    Rl = fit(C.Spec(y, Z, ZZ), SIG_U)
    rr = y[Rl.idx] - Rl.smooth
    off = 0.6 if kl == "pp" else -0.6
    ax.errorbar(x[Rl.idx] + off, rr, yerr=np.sqrt(np.maximum(Rl.model, 1)), fmt=mk, ms=4.5, color=col,
                label="%s minus its smooth part: $h=%.1f\\pm%.1f$" % (lab, Rl.h, Rl.eh))
    ax.plot(x[Rl.idx], Rl.h * T[Rl.idx], "-", color=col, lw=2)
ax.axhline(0, color="k", lw=0.8); ax.axvline(0, color="k", lw=0.8, ls=":")
ax.set_xlabel("$u$ [MeV/$c^2$]"); ax.set_ylabel("rows per 2.5 MeV/$c^2$")
ax.set_title("(c) like-sign rows alone (same fit)")
ax.set_ylim(-80, 185)
ax.legend(loc="upper right", fontsize=11.5)
fig.tight_layout()
fig.savefig(FIGX + "edgeu_data_step.pdf"); fig.savefig(FIGX + "edgeu_data_step.png", dpi=130)
P("figure written: %sedgeu_data_step.{pdf,png}" % FIGX)
