#!/usr/bin/env python3
# ============================================================================
# peaktail_figs.py  --  package peaktail: figures (Part III colours: data blue, C4 green, C5 cyan, fit green, bedge brown)
# ============================================================================
# Written 2026-11-15 (prefix peaktail). Inputs: the columns of peaktail_dump.C and the json / npz files of peaktail_sim,
# _ufit, _line, _rad, _prop, _summary. Outputs (new files): docs/latex/figs_all/peaktail_{sources,constraint,effect}.{pdf,png}
# Run: python3 -I plots/all_report/peaktail_figs.py
# ----------------------------------------------------------------------------
# _cnt copy (2026-11-15, FIGURE AGENT C of the integration run): the analyst calls the stored pi+pi-e+e- combinations
# COUNTS, not rows. Only label, legend, title and annotation strings of peaktail_sources ('rows' -> 'counts'; the template
# statistics 'C4: 180 simulated rows' -> 'C4: 180 simulated entries') and the output name differ from peaktail_figs.py;
# the inputs (the same cached peaktail columns and json / npz files) and the plotted content are the same. Only
# peaktail_sources is written (as peaktail_sources_cnt); peaktail_constraint and peaktail_effect are used in no report.
# Run: python3 -I plots/all_report/peaktail_figs_cnt.py
# ============================================================================
import json
import math
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import edgeu_core as C
import peaktail_core as P

plt.rcParams.update({"font.size": 13, "axes.labelsize": 14, "axes.titlesize": 13.5, "legend.fontsize": 10.5, "xtick.labelsize": 12,
                     "ytick.labelsize": 12, "pdf.fonttype": 42})
COL = dict(data="#1f3fbf", c4="#009900", c5="#00b5c8", fit="#009900", bedge="#8c564b", qed="#d62728", chk="#ff8c00", grey="#888888",
           tag="#9467bd", x17="#e377c2", ap="#1f3fbf", dr="#555555")
N = json.load(open(P.OUT + "peaktail_numbers.json"))
UF = json.load(open(P.OUT + "peaktail_ufit.json")); SJ = json.load(open(P.OUT + "peaktail_sim.json"))
S = P.load_sim(); w = P.sim_weight(S); gS = S["g"]
us, mpp, m4p = C.u_sim(S); d = P.resid(S)
M4r, J = P.m4_and_jacobian(S)
dl = C.S_NOM * (J[0] * (S["ep_p"] - S["ep_ptrue"]) + J[1] * (S["em_p"] - S["em_ptrue"])).astype(np.float64) * 1e-3
dpi_ = C.S_NOM * (J[2] * (S["pip_p"] - S["pip_ptrue"]) + J[3] * (S["pim_p"] - S["pim_ptrue"])).astype(np.float64) * 1e-3
flag = dl < -0.005
inR = (mpp >= 0.30) & (mpp < 0.40)
td, sd = N["t_data"]["t"], N["t_data"]["s"]

# ============================================================================== 1. sources (simulation)
fig, axs = plt.subplots(2, 2, figsize=(13.0, 10.0))
ax = axs[0, 0]
for gg, nm, col, bw in ((3, "C5 direct", COL["c5"], 2.0), (2, "C4 conversion", COL["c4"], 5.0)):
    bins = np.arange(-80, 40.01, bw)
    m = (gS == gg) & (w > 0) & inR
    tot = w[m].sum()
    ax.hist(1e3 * d[m], bins=bins, weights=w[m] / tot / bw, histtype="step", color=col, lw=2.2, label="%s: all counts" % nm)
    ax.hist(1e3 * d[m & flag], bins=bins, weights=w[m & flag] / tot / bw, histtype="stepfilled", color=col, alpha=0.30,
            label="%s: tail counts (%.0f %%)" % (nm, 100 * w[m & flag].sum() / tot))
    ax.hist(1e3 * d[m & ~flag], bins=bins, weights=w[m & ~flag] / tot / bw, histtype="step", color=col, lw=1.3, ls="--")
ax.set_yscale("log"); ax.set_ylim(2e-5, 8.0); ax.set_xlim(-80, 40)
ax.set_xlabel("$d=M_4'-0.99615\\,m_{4,\\mathrm{true}}$ [MeV/$c^2$]"); ax.set_ylabel("fraction per MeV/$c^2$")
ax.set_title("(a) SMASH peak counts, $0.30\\leq M(\\pi^+\\pi^-)<0.40$")
ax.text(39, 6.0, "dashed: counts without a lepton loss\nfilled: tail counts (the leptons' loss\nshifts $M_4$ by more than 5 MeV/$c^2$)\nC4: 180 simulated entries, 5 MeV bins", fontsize=8.8, va="top", ha="right")
ax.legend(loc="upper left", fontsize=8.8, handlelength=1.5)   # _cnt: 9.5 -> 8.8 and shorter handles, so that the longer word fits
ax = axs[0, 1]
labs = []; vals = []
for gg, nm in ((2, "C4"), (3, "C5")):
    for thr in (-0.015, -0.025):
        m = (gS == gg) & (w > 0) & inR & (d < thr)
        ww = w[m]
        vals.append([np.average(1e3 * dl[m], weights=ww), np.average(1e3 * dpi_[m], weights=ww), np.average(1e3 * (d - dl - dpi_)[m], weights=ww)])
        labs.append("%s, $d<%d$" % (nm, int(1e3 * thr)))
vals = np.array(vals); yy = np.arange(len(labs))
ax.barh(yy, vals[:, 0], color="#1f77b4", label="lepton momenta")
ax.barh(yy, vals[:, 1], left=vals[:, 0], color="#ff7f0e", label="pion momenta")
ax.barh(yy, vals[:, 2], left=vals[:, 0] + vals[:, 1], color="#7f7f7f", label="angles, rest")
ax.set_yticks(yy); ax.set_yticklabels(labs); ax.axvline(0, color="k", lw=0.8); ax.set_xlim(-52, 32)
ax.set_xlabel("mean shift of $M_4$ [MeV/$c^2$]")
ax.set_title("(b) what makes the tail: shift split by particle")
gen = SJ["generator"]
ax.text(0.63, 0.04, "generated $m_4$: a delta\nfunction (C5 548.000,\nC4 547.46 MeV/$c^2$):\nno radiative photon\nin the generator", transform=ax.transAxes, fontsize=10)
ax.legend(loc="upper right", fontsize=10)
ax = axs[1, 0]
ub = np.arange(-60, 80.01, 2.5)
for gg, nm, col in ((2, "C4", COL["c4"]), (3, "C5", COL["c5"])):
    m = (gS == gg) & (w > 0) & inR & np.isfinite(us)
    ax.hist(1e3 * us[m & ~flag], bins=ub, weights=C.AP * w[m & ~flag], histtype="step", color=col, lw=1.5, ls="--", label="%s core $\\times a_P$" % nm)
    ax.hist(1e3 * us[m & flag], bins=ub, weights=C.AP * w[m & flag], histtype="stepfilled", color=col, alpha=0.35, label="%s tail $\\times a_P$" % nm)
ax.axvspan(-45, 45, color="#eeeeee", zorder=0); ax.axvline(0, color="k", lw=1)
ax.set_xlabel("$u=M_4-E(M_{\\pi\\pi},M_{ee})$ [MeV/$c^2$]"); ax.set_ylabel("expected data counts per 2.5 MeV/$c^2$")
ax.set_title("(c) the peak in $u$: the tail reaches the edge")
ax.text(-58, ax.get_ylim()[1] * 0.85, "grey: fit range\nedgeu", fontsize=10.5)
ax.legend(loc="upper left", bbox_to_anchor=(0.0, 0.78), fontsize=10)
ax = axs[1, 1]
RR = np.load(P.OUT + "peaktail_rad_rows.npz")
ppb = np.array([0.26, 0.30, 0.34, 0.37, 0.40, 0.42, 0.46, 0.50])
f15, f15q, ee = [], [], []
for a, b in zip(ppb[:-1], ppb[1:]):
    m = ((gS == 2) | (gS == 3)) & (w > 0) & (mpp >= a) & (mpp < b)
    f15.append(P.frac(d[m] < -0.015, w[m]))
    m4_ = (gS == 2) & (w > 0) & (mpp >= a) & (mpp < b)
    mr = (RR["mpp"] >= a) & (RR["mpp"] < b)
    num = (w[m4_ & (d < -0.015)]).sum() + RR["w"][mr & (RR["d1"] < -0.015)].sum(); den = w[m4_].sum() + RR["w"][mr].sum()
    f15q.append(num / den)
    ee.append(1e3 * (C.ME ** 2 - (0.5 * (a + b) / C.EDGE_PIPI) ** 2) / (2 * C.ME))
xm = 0.5 * (ppb[1:] + ppb[:-1])
ax.step(ppb, np.r_[f15, f15[-1]], where="post", color="k", lw=2, label="SMASH (C4 + C5)")
ax.step(ppb, np.r_[f15q, f15q[-1]], where="post", color=COL["qed"], lw=2, ls="--", label="+ QED radiation of C5 (not simulated)")
ax.axvspan(0.30, 0.40, color="#e8f3e8", zorder=0); ax.axvspan(0.42, 0.46, color="#e8eef8", zorder=0)
ax.text(0.305, 0.165, "u fit", fontsize=11); ax.text(0.421, 0.165, "no bump", fontsize=11)
ax.set_xlabel("$M(\\pi^+\\pi^-)$ [GeV/$c^2$]"); ax.set_ylabel("fraction of peak counts with $d<-15$ MeV/$c^2$")
ax.set_ylim(0, 0.19); ax.set_xlim(0.26, 0.50)
ax.set_title("(d) the tail shrinks with $M(\\pi^+\\pi^-)$ (softer leptons)")
for x, e in zip(xm, ee):
    ax.text(x, 0.004, "%.0f" % e, ha="center", fontsize=9, color="#555555")
ax.text(0.262, 0.012, "$E^*_{ee}$ [MeV]:", fontsize=9, color="#555555")
ax.legend(loc="upper right", bbox_to_anchor=(1.0, 0.84), fontsize=10)
fig.tight_layout()
fig.savefig(P.FIG + "peaktail_sources_cnt.pdf"); fig.savefig(P.FIG + "peaktail_sources_cnt.png", dpi=110)
plt.close(fig)

# ============================================================================== 2. the data constraint
UFD = dict(np.load(P.OUT + "peaktail_ufitfig.npz")); LFD = dict(np.load(P.OUT + "peaktail_linefig.npz"))
fig = plt.figure(figsize=(13.0, 10.5))
gsp = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0])
ax = fig.add_subplot(gsp[0, 0])
x = UFD["x"] * 1e3; idx = UFD["R0_idx"]
resid = UFD["sig"][idx] - UFD["R0_smooth"] - UFD["R0_h"] * UFD["T"][idx]
err = np.sqrt(np.maximum(UFD["N"][idx], 1) + UFD["cbv"][idx])
ax.errorbar(x[idx], resid, yerr=err, fmt="o", ms=4, color=COL["data"], label="data minus the smooth part and the step")
ax.plot(x[idx], C.AP * UFD["K1"][idx], "-", color=COL["fit"], lw=2.2, label="SMASH peak, $t=1$ (edgeu nominal)")
ax.plot(x[idx], C.AP * UFD["Kt0"][idx], "--", color="k", lw=1.6, label="SMASH peak without tail, $t=0$")
ax.plot(x[idx], C.AP * UFD["Kt2"][idx], ":", color=COL["qed"], lw=2.0, label="doubled tail, $t=2$")
ax.plot(x[idx], UFD["Rk_amp"] * UFD["KG"][idx], "-.", color=COL["chk"], lw=1.8, label="checker's Gaussian (amplitude fitted)")
ax.axhline(0, color="k", lw=0.7); ax.axvline(0, color="k", lw=0.8, ls=":")
ax.set_xlim(-46, 46); ax.set_ylim(-190, 600)
ax.set_xlabel("$u$ [MeV/$c^2$]"); ax.set_ylabel("rows per 2.5 MeV/$c^2$")
ax.set_title("(a) u fit: what is left for the peak term")
ax.legend(loc="upper left", fontsize=9.5)
ax = fig.add_subplot(gsp[0, 1])
x4 = LFD["x"] * 1e3; i4 = LFD["idx"]
o = 2
res3 = LFD["sig"][i4] - LFD["s3"]
ax.errorbar(x4[i4], res3, yerr=LFD["err"][i4], fmt="o", ms=4, color=COL["data"], label="data minus the smooth part ($t=1$ fit)")
ax.plot(x4[i4], LFD["p3"][o + 1] * (LFD["K0"] + LFD["KK"])[i4], "-", color=COL["fit"], lw=2.2, label="peak, $t=1$: amp. %.2f" % LFD["p3"][o + 1])
ax.plot(x4[i4], LFD["s4"] + LFD["p4"][o + 1] * LFD["K0"][i4] - LFD["s3"], "--", color="k", lw=1.6, label="fit with $t=0$ (no tail)")
ax.axhline(0, color="k", lw=0.7); ax.axvline(543.9, color="k", lw=0.8, ls=":"); ax.set_ylim(-220, 345)
ax.set_xlabel("$M(\\pi^+\\pi^-e^+e^-)$ [MeV/$c^2$]"); ax.set_ylabel("rows per 2.5 MeV/$c^2$")
ax.set_title("(b) no bump: $0.42\\leq M(\\pi^+\\pi^-)<0.46$ GeV/$c^2$")
lc = N["t_ctrl"]["0.42-0.46|r1.000|0.49-0.60 pol2"]
ax.text(492, -150, "$t=%.2f\\pm%.2f$ (amplitude free)" % (lc["t"], lc["stat_scaled"]), fontsize=11, va="center")
ax.legend(loc="upper left", fontsize=9.5)
ax = fig.add_subplot(gsp[1, :])
items = []
for k, lab in (("-45..+45 pol2", "u: $-45..+45$, pol2 (edgeu range)"), ("-45..+55 pol2", "u: $-45..+55$, pol2"), ("-45..+60 pol2", "u: $-45..+60$, pol2"),
               ("-45..+70 pol3", "u: $-45..+70$, pol3"), ("-55..+70 pol3", "u: $-55..+70$, pol3")):
    v = N["t_u"]["SMASH|" + k]; items.append((lab, v["t"], math.hypot(v["stat"], v["mc"]), COL["data"]))
for k, lab in (("0.42-0.46|r1.000|0.49-0.60 pol2", "no bump: 490-600, pol2"), ("0.42-0.46|r1.000|0.48-0.62 pol3", "no bump: 480-620, pol3"),
               ("0.42-0.46|r1.000|0.50-0.59 pol2", "no bump: 500-590, pol2"), ("0.42-0.46|r1.000|0.47-0.62 pol3", "no bump: 470-620, pol3"),
               ("0.42-0.46|r0.893|0.49-0.60 pol2", "no bump: 490-600, data-like core width")):
    v = N["t_ctrl"][k]; items.append((lab, v["t"], v["stat_scaled"], "#2ca02c"))
for k, lab in (("E1|SMASH|-45..+60 pol2", "u, $M_{ee}<10$ (conversions + direct)"), ("E2|SMASH|-45..+60 pol2", "u, $M_{ee}\\geq10$ (direct only)")):
    v = N["t_u_class"][k]; items.append((lab, v["t"], v["stat"], "#7f7f7f"))
yy = np.arange(len(items))[::-1]
for (lab, t, e, col), y in zip(items, yy):
    ax.errorbar(t, y, xerr=e, fmt="o", ms=6, color=col, capsize=3)
ax.set_yticks(yy); ax.set_yticklabels([i[0] for i in items], fontsize=11)
ax.axvspan(td - sd, td + sd, color="#cccccc", alpha=0.6, zorder=0)
ax.axvline(td, color="k", lw=1.5); ax.axvline(1.0, color=COL["fit"], lw=2, ls="--"); ax.axvline(N["t_qed"]["t_combined"], color=COL["qed"], lw=2, ls=":")
tchk = 1 + N["checker"]["rel"] / N["f_u_slope"]
ax.axvline(0.0, color="k", lw=1.2, ls="-."); ax.axvline(tchk, color=COL["chk"], lw=1.8, ls="-.")
ax.set_xlim(-1.2, 3.3)
ax.set_xlabel("tail scale $t$ (1 = SMASH; tail = rows whose leptons' loss shifts $M_4$ by more than 5 MeV/$c^2$)")
ax.set_title("(c) the tail in the data: $t_{\\mathrm{data}}=%.2f\\pm%.2f$ (grey band)" % (td, sd))
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
hl = [Patch(color="#cccccc", alpha=0.6, label="data: $t_{\\mathrm{data}}\\pm1\\sigma$"), Line2D([], [], color=COL["fit"], lw=2, ls="--", label="SMASH, $t=1$"),
      Line2D([], [], color=COL["qed"], lw=2, ls=":", label="SMASH + QED radiation of C5, $t=%.2f$" % N["t_qed"]["t_combined"]),
      Line2D([], [], color="k", lw=1.2, ls="-.", label="no tail, $t=0$"),
      Line2D([], [], color=COL["chk"], lw=1.8, ls="-.", label="the checker's 0.384 ($t\\approx%.2f$)" % tchk),
      Line2D([], [], color=COL["data"], marker="o", ls="", label="u fit (stat + MC templates)"),
      Line2D([], [], color="#2ca02c", marker="o", ls="", label="no-bump line shape (stat, scaled)"),
      Line2D([], [], color="#7f7f7f", marker="o", ls="", label="u fit per $M_{ee}$ class (stat)")]
ax.legend(handles=hl, loc="upper right", fontsize=10)
ax.set_ylim(-1, len(items))
fig.tight_layout()
pass   # peaktail_constraint: not written by the _cnt copy
plt.close(fig)

# ============================================================================== 3. the effect
fig, axs = plt.subplots(1, 2, figsize=(13.0, 5.6))
ax = axs[0]
ft = UF["f_of_t"]
ax.plot([v["t"] for v in ft], [v["f"] for v in ft], "-o", color="k", ms=4, lw=2, label="$f_\\eta$(u-edge) (edgeu)")
tb = [("t0", 0.0), ("t_data-1s", td - sd), ("t_data", td), ("nominal", 1.0), ("t_data+1s", td + sd), ("t2", 2.0)]
ax.plot([t for _, t in tb], [N["f_bedge"][k]["f"] for k, _ in tb], "-s", color=COL["bedge"], ms=5, lw=1.8, label="$f_\\eta$(edge) per slice (bedge)")
ax.axvspan(td - sd, td + sd, color="#cccccc", alpha=0.6, zorder=0)
ax.axvline(1.0, color=COL["fit"], lw=1.5, ls="--")
ax.plot([N["t_qed"]["t_combined"]], [N["f_u"]["qed"]["f"]], "*", color=COL["qed"], ms=14, label="QED radiation of C5 added (%.3f)" % N["f_u"]["qed"]["f"])
ax.axhline(N["checker"]["f"], color=COL["chk"], lw=1.8, ls="-.", label="checker's Gaussian peak term: %.3f" % N["checker"]["f"])
ax.axhline(N["decomp"]["chk_fairSMASH"]["f"], color=COL["chk"], lw=1.2, ls=":", label="the same fit, SMASH peak in $h_{\\mathrm{sim}}$: %.3f" % N["decomp"]["chk_fairSMASH"]["f"])
ax.set_xlim(-0.1, 3.1); ax.set_ylim(0.37, 0.49)
ax.set_xlabel("tail scale $t$"); ax.set_ylabel("$\\eta$ scale of the bump $f_\\eta$")
ax.set_title("(a) the bump edge measurements vs the tail")
ax.text(td, 0.372, "data", ha="center", fontsize=10.5)
ax.legend(loc="upper left", fontsize=9.8)
ax = axs[1]
tq = [("t0", 0.0), ("t_data-1s", td - sd), ("t_data", td), ("nominal", 1.0), ("t_data+1s", td + sd), ("t2", 2.0)]
def g(dct, k, key):
    return 1.0 if k == "nominal" else dct[k][key]
ax.plot([t for _, t in tq], [100 * (N["f_u"][k]["f"] / N["f_u"]["nominal"]["f"] - 1) for k, _ in tq], "-o", color="k", lw=2, ms=4, label="$f_\\eta$(u-edge)")
ax.plot([t for _, t in tq], [100 * (g(N["sidebands"], k, "aP_rel") - 1) for k, _ in tq], "-s", color=COL["ap"], lw=2, ms=5, label="$a_P$ and $\\sigma_\\eta$ from the peak side bands")
ax.plot([t for _, t in tq], [0.0 if k == "nominal" else 100 * (N["DR_rel"][k] - 1) for k, _ in tq], "-^", color=COL["dr"], lw=1.6, ms=5,
        label="double ratio peak/bump")
ax.plot([t for _, t in tq], [0.0 if k == "nominal" else 100 * (N["tag_proxy"][k]["ftag_rel"] - 1) for k, _ in tq], "--d", color=COL["tag"], lw=1.6, ms=5,
        label="ECAL tag (proxy: Gaussian content + pol2)")
xk = [("t0", 0.0), ("t_data-1s", td - sd), ("nominal", 1.0), ("t_data+1s", td + sd), ("t2", 2.0)]
ax.plot([t for _, t in xk], [0.0 if k == "nominal" else 100 * (1 / N["x17_RX_rel"][k] - 1) for k, _ in xk], ":v", color=COL["x17"], lw=1.6, ms=5,
        label="X17 limit, normalisation (1)")
ax.axvspan(td - sd, td + sd, color="#cccccc", alpha=0.6, zorder=0); ax.axhline(0, color="k", lw=0.7); ax.axvline(1.0, color=COL["fit"], lw=1.5, ls="--")
ax.set_xlim(-0.1, 2.1); ax.set_ylim(-17, 23)
ax.set_xlabel("tail scale $t$"); ax.set_ylabel("change against $t=1$ [%]")
ax.set_title("(b) other Part III quantities")
ax.legend(loc="upper left", fontsize=9.8)
fig.tight_layout()
pass   # peaktail_effect: not written by the _cnt copy
plt.close(fig)
print("written: peaktail_sources_cnt (.pdf, .png)")
