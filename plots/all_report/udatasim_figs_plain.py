#!/usr/bin/env python3
# ============================================================================
# udatasim_figs_plain.py -- COPY of udatasim_figs.py (2026-11-16, short-report figure agent for report_all_short.tex): writes only
# udatasim_u_overlay_plain.{pdf,png}. The visible texts (suptitle, panel title (a), the "fit range" marks, the in-panel lambda line
# and the shared legend) name selections and SMASH classes in plain words (readability instruction I2: no S0, C2-C5, G_eta, G_D,
# G_C, V3, D38, edgeu, stage A inside a short-report figure). Same cached inputs, same drawing code, every number unchanged.
# (LAB_A / LAB_B and the per-panel labels are not displayed: the figure has only the shared legend; they are left as they were.)
# Run: python3 -I plots/all_report/udatasim_figs_plain.py
# ============================================================================
# ============================================================================
# udatasim_figs.py  --  package udatasim: figures "data against SMASH in u"
# ============================================================================
# Written 2026-11-15 (prefix udatasim). New file. Inputs: /hdd3/przygoda/HADES/FAT_out/all_report/udatasim/udatasim_figdata.npz
# and udatasim_overlay.json (udatasim_overlay.py). Outputs (new files only):
#   docs/latex/figs_all/udatasim_u_overlay.{pdf,png}         S0, 0.30 <= M(pipi) < 0.40, all M(ee): (a) SMASH as it is, (b) measured scales
#   docs/latex/figs_all/udatasim_u_overlay_allmpp.{pdf,png}  the same for all M(pipi) with a defined edge (0.27-0.4116)
#   docs/latex/figs_all/udatasim_u_overlay_mee.{pdf,png}     model (b) per M(ee) class (E1 < 10, E2 >= 10 MeV/c2), 0.30-0.40
# Part III colours: data blue, CB red, C2 orange, C3 magenta, C4 green, C5 cyan, G_eta violet, G_D light grey, G_C dark grey
# (as Fig. 12 / fix_C3). Unit: counts (stored pi+pi-e+e- combinations) per 2.5 MeV/c2; simulation in expected counts.
# Run: python3 -I plots/all_report/udatasim_figs.py
# ============================================================================
import json
import math

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

O = "/hdd3/przygoda/HADES/FAT_out/all_report/udatasim/"
FIG = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_all/"
FD = dict(np.load(O + "udatasim_figdata.npz"))
R = json.load(open(O + "udatasim_overlay.json"))
plt.rcParams.update({"font.size": 14, "axes.labelsize": 15, "axes.titlesize": 14.5, "legend.fontsize": 11.5, "xtick.labelsize": 13,
                     "ytick.labelsize": 13, "pdf.fonttype": 42})
COL = dict(data="#1f3fbf", cb="#e8202a", C2="#ff8c00", C3="#d400d4", C4="#009900", C5="#00b5c8", Geta="#7b3fbf", GD="#c8c8c8", GC="#7f7f7f")
x = FD["x"] * 1e3
ed = FD["edges"] * 1e3
F38, AP, AG = float(FD["F38"]), float(FD["AP"]), float(FD["AG"])
SA = R["stageA"]
STACK = [("GD", ("GDom", "GDnom")), ("GC", ("GCom", "GCnom")), ("Geta", ("Geta",)), ("C5", ("C5",)), ("C4", ("C4",)), ("C3", ("C3",)), ("C2", ("C2",))]
LAB_A = dict(GD="$G_D$ non-$\\eta$, Dalitz/vector", GC="$G_C$ non-$\\eta$, conversion", Geta="$G_\\eta$ other $\\eta$ decays",
             C5="C5 peak, direct $\\eta\\to\\pi^+\\pi^-e^+e^-$", C4="C4 peak, $\\pi^+\\pi^-\\gamma$, $\\gamma$ converted",
             C3="C3 bump, $\\pi^0\\to\\gamma\\gamma$, $\\gamma$ converted", C2="C2 bump, $\\pi^0$ Dalitz (D38 weights)")
LAB_B = dict(GD="$G_D\\times c_j(1+bx)$ (stage A)", GC="$G_C\\times c_j(1+bx)\\,\\kappa_B$, $\\kappa_B$ = %.2f" % SA["kB"],
             Geta="$G_\\eta\\times a_\\eta g$ = %.3f" % AG, C5="C5 $\\times\\,a_P$ = %.3f" % AP, C4="C4 $\\times\\,a_P$ = %.3f" % AP,
             C3="C3 $\\times\\,f_\\eta(u)$ = %.3f" % F38, C2="C2 (D38) $\\times\\,f_\\eta(u)$ = %.3f" % F38)
XLAB = "$u=M_4-E(M_{\\pi\\pi},M_{ee})$ [MeV/$c^2$]"
YLAB = "counts per 2.5 MeV/$c^2$"


def save(fig, name):
    if name != "udatasim_u_overlay":   # _plain: only the figure of the short report
        plt.close(fig); return
    fig.savefig(FIG + name + "_plain.pdf"); fig.savefig(FIG + name + "_plain.png", dpi=130); plt.close(fig)


def parts(tag, mk):
    return {k: sum(FD[tag + mk + "_" + g_] for g_ in gs) for k, gs in STACK}


def ratio_pts(tag, mk, M=None):
    d = FD[tag + "sig"]; M = sum(parts(tag, mk).values()) if M is None else M
    VD = np.diag(FD[tag + "VD"]); VM = np.diag(FD[tag + "VM" + mk])
    r = d / M
    e = np.sqrt(VD / M ** 2 + d ** 2 * VM / M ** 4)
    return r, e


def stairs_fill(ax, lo, hi, col, lab):
    ax.fill_between(ed, np.append(lo, lo[-1]), np.append(hi, hi[-1]), step="post", color=col, lw=0, label=lab)


def lam_model(tag, rk, cls):
    """model (b) with one free factor on G_D + G_C (fitted in -45...+45; display only)"""
    p = parts(tag, "b")
    lam = R[rk][cls]["b"]["chi2"]["fit"]["lam"]
    return sum(p.values()) + (lam - 1.0) * (p["GD"] + p["GC"]), lam


def mark(ax, top=False):
    ax.axvspan(-45, 45, color="#eeeeee", zorder=0)
    ax.axvline(0, color="k", lw=1.2, zorder=5)


def panel(axm, axr, tag, rk, cls, mk, ymax, rlo, rhi, title, legend=True, showlam=True):
    p = parts(tag, mk)
    cum = np.zeros_like(x)
    for k, _ in STACK:
        stairs_fill(axm, cum, cum + p[k], COL[k], (LAB_A if mk == "a" else LAB_B)[k])
        cum = cum + p[k]
    M = cum
    d = FD[tag + "sig"]
    axm.errorbar(x, d, yerr=np.sqrt(np.diag(FD[tag + "VD"])), fmt="o", ms=4.2, color=COL["data"], zorder=10, label="data S0: $N-$CB")
    axm.step(ed, np.append(FD[tag + "cb"], FD[tag + "cb"][-1]), where="post", color=COL["cb"], lw=1.6, ls="-", zorder=9,
             label="CB $=2\\sqrt{N_{++}N_{--}}$ (subtracted)")
    if mk == "b" and showlam:
        Ml, lam = lam_model(tag, rk, cls)
        axm.step(ed, np.append(Ml, Ml[-1]), where="post", color="k", lw=1.8, ls="--", zorder=9,
                 label="(b) with $G_D+G_C\\times\\lambda$ = %.2f (fit $-$45...+45, display)" % lam)
    mark(axm)
    axm.plot([-45, 45], [0.975 * ymax] * 2, color="#888888", lw=2.0, solid_capstyle="butt", zorder=6)
    axm.text(-30, 0.955 * ymax, "$u$-edge fit range", color="#555555", ha="center", va="top", fontsize=11.5)
    axm.plot([20, 60], [0.94 * ymax] * 2, color=COL["C4"], lw=5, solid_capstyle="butt", zorder=6)
    axm.text(62, 0.94 * ymax, "peak region", color=COL["C4"], ha="left", va="center", fontsize=11.5)
    axm.set_xlim(ed[0], ed[-1]); axm.set_ylim(0, ymax)
    axm.set_ylabel(YLAB); axm.set_title(title, loc="left")
    plt.setp(axm.get_xticklabels(), visible=False)
    o = R[rk][cls][mk]
    txt = ("data / model: %.2f (below the edge), %.2f (at the edge), %.2f (above)\n"
           "$\\chi^2$/ndf (stat): %.0f/%d in $-$45...+45, %.0f/%d in +20...+60" %
           (o["bands"]["below"]["R"], o["bands"]["edge"]["R"], o["bands"]["above"]["R"],
            o["chi2"]["fit"]["chi2_stat"], o["chi2"]["fit"]["n"], o["chi2"]["peak"]["chi2_stat"], o["chi2"]["peak"]["n"]))
    if mk == "b":
        txt += "\nwith the SMASH background $\\times\\lambda$: %.0f/%d and %.0f/%d" % (o["chi2"]["fit"]["chi2_lam"], o["chi2"]["fit"]["n"] - 1,
                                                                         o["chi2"]["peak"]["chi2_lam"], o["chi2"]["peak"]["n"] - 1)
    axm.text(0.015, 0.875, txt, transform=axm.transAxes, fontsize=11.2, va="top",
             bbox=dict(facecolor="white", edgecolor="#999999", alpha=0.92), zorder=12)
    # ratio pad
    r, e = ratio_pts(tag, mk)
    mark(axr)
    axr.axhline(1.0, color="k", lw=0.9, ls=":")
    if mk == "b":
        ea = np.sqrt(np.diag(FD[tag + "VA"])) / M
        axr.fill_between(ed, np.append(1 - ea, 1 - ea[-1]), np.append(1 + ea, 1 + ea[-1]), step="post", color="#b0b0b0", alpha=0.55, lw=0,
                         label="stage-A normalisation error of the model")
        if showlam:
            Ml, lam = lam_model(tag, rk, cls)
            rl, el = ratio_pts(tag, mk, Ml)
            axr.errorbar(x + 0.5, rl, yerr=el, fmt="s", ms=3.6, mfc="white", color="#444444", lw=0.9, label="data / [(b) with $\\lambda$]")
    axr.errorbar(x, r, yerr=e, fmt="o", ms=4.0, color=COL["data"], label="data / model, stat. (data $\\oplus$ simulation MC)")
    axr.set_xlim(ed[0], ed[-1]); axr.set_ylim(rlo, rhi)
    axr.set_xlabel(XLAB); axr.set_ylabel("data / model")


def ylim_of(tags):
    m = 0.0
    for tag, mk in tags:
        M = sum(parts(tag, mk).values())
        m = max(m, M.max(), (FD[tag + "sig"] + np.sqrt(np.diag(FD[tag + "VD"]))).max())
    return 1000.0 * math.ceil(1.30 * m / 1000.0)


def rlim_of(tags):
    lo, hi = 10.0, 0.0
    for tag, mk in tags:
        r, e = ratio_pts(tag, mk)
        lo = min(lo, (r - e).min()); hi = max(hi, (r + e).max())
    return math.floor(10 * (lo - 0.12)) / 10.0, math.ceil(10 * (hi + 0.05)) / 10.0


LAB_L = dict(GD="Dalitz and vector-meson pairs (no conversion); (b): $\\times\\,c_j(1+bx)$, normalised above the peak",
             GC="non-$\\eta$ conversion pairs; (b): $\\times\\,c_j(1+bx)\\,\\kappa_B$, $\\kappa_B$ = %.2f" % SA["kB"],
             Geta="other $\\eta$ decays; (b): $\\times\\,a_\\eta g$ = %.3f" % AG,
             C5="peak, direct $\\eta\\to\\pi^+\\pi^-e^+e^-$; (b): $\\times\\,a_P$ = %.3f" % AP,
             C4="peak, $\\eta\\to\\pi^+\\pi^-\\gamma$, $\\gamma$ converted; (b): $\\times\\,a_P$",
             C3="bump, $\\eta\\to\\pi^+\\pi^-\\pi^0$, $\\pi^0\\to\\gamma\\gamma$, one $\\gamma$ converted; (b): $\\times\\,f_\\eta(u)$",
             C2="bump, $\\eta\\to\\pi^+\\pi^-\\pi^0$, $\\pi^0$ Dalitz, QED weight per count; (b): $\\times\\,f_\\eta(u)$ = %.3f" % F38)


def shared_legend(fig, ax, with_lam=True):
    h = [Patch(color=COL[k], label=LAB_L[k]) for k, _ in reversed(STACK)]
    h += [Line2D([], [], marker="o", ls="none", color=COL["data"], ms=6, label="data, both pions in the vertex fit: $N-$CB (stat.)"),
          Line2D([], [], color=COL["cb"], lw=1.8, label="CB $=2\\sqrt{N_{++}N_{--}}$ of the data (already subtracted)")]
    if with_lam:
        h += [Line2D([], [], color="k", lw=1.8, ls="--", label="(b) with the SMASH background scaled by one free $\\lambda$ (fit in $-$45...+45; display only)"),
              Line2D([], [], marker="s", ls="none", mfc="white", color="#444444", ms=6, label="ratio pad: data / [(b) with $\\lambda$]")]
    h += [Patch(color="#b0b0b0", alpha=0.55, label="ratio pad (b): error of the background normalisation of the model"),
          Patch(color="#eeeeee", label="grey shading: $u$-edge fit range $-$45...+45 MeV/$c^2$")]
    ax.axis("off")
    ax.legend(handles=h, loc="center", ncol=2, fontsize=11.2, frameon=False, columnspacing=2.0, handlelength=2.2)


ORIGPOS = [(0.0693587302, 0.4932924065, 0.4197126984, 0.4117566131), (0.5584301587, 0.4932924065, 0.4197126984, 0.4117566131),
           (0.0693587302, 0.2396174243, 0.4197126984, 0.1830029392), (0.5584301587, 0.2396174243, 0.4197126984, 0.1830029392),
           (0.0693587302, 0.0262030888, 0.908784127, 0.1427422925)]


def two_panel(rk, cls, name, sup):
    tag = "%s_%s_" % (rk, cls)
    fig = plt.figure(figsize=(17.5, 13.6))
    gs = fig.add_gridspec(3, 2, height_ratios=[2.25, 1, 0.78])
    axs = np.array([[fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])], [None, None]], dtype=object)
    axs[1, 0] = fig.add_subplot(gs[1, 0], sharex=axs[0, 0]); axs[1, 1] = fig.add_subplot(gs[1, 1], sharex=axs[0, 1])
    axl = fig.add_subplot(gs[2, :])
    ymax = ylim_of([(tag, "a"), (tag, "b")])
    rlo, rhi = rlim_of([(tag, "a"), (tag, "b")])
    panel(axs[0, 0], axs[1, 0], tag, rk, cls, "a", ymax, rlo, rhi, "(a) SMASH as it is ($\\sigma_{\\mathrm{prod}}\\times$PDG BR$\\times\\eta$-decay weights, no scale)")
    panel(axs[0, 1], axs[1, 1], tag, rk, cls, "b", ymax, rlo, rhi, "(b) SMASH with the measured scales")
    shared_legend(fig, axl)
    fig.suptitle(sup, fontsize=15)
    fig.tight_layout(rect=(0, 0, 1, 0.975))
    if name == "udatasim_u_overlay":   # _plain: the axes positions of the original udatasim_u_overlay.pdf (its tight_layout result),
        for a_, pos in zip(fig.axes, ORIGPOS):   # so that the longer plain legend labels do not move the panels
            a_.set_position(pos)
    save(fig, name)
    return ymax, (rlo, rhi)


# the D40 band on this display (largest relative change of model (b) per bin in the plotted range)
for rk in ("r1", "all"):
    tag = "%s_all_" % rk
    Mb = sum(parts(tag, "b").values())
    dd = [float(np.max(np.abs(FD[tag + "Mb_t%.1f" % t] / Mb - 1.0))) for t in (0.5, 1.5)]
    i = int(np.argmax(np.abs(FD[tag + "Mb_t0.5"] / Mb - 1.0)))
    print("D40 on the display, %s: largest |model(t)/model(t = 1) - 1| = %.4f (t = 0.5), %.4f (t = 1.5), at u = %+.1f MeV/c2" % (rk, dd[0], dd[1], x[i]))

y1 = two_panel("r1", "all", "udatasim_u_overlay",
               "Data and SMASH in the edge variable $u$: both pions in the vertex fit, $0.30\\leq M(\\pi^+\\pi^-)<0.40$ GeV/$c^2$ (the $u$-edge fit region), all $M(e^+e^-)$")
y2 = two_panel("all", "all", "udatasim_u_overlay_allmpp",
               "Data and SMASH in $u$: S0, all $M(\\pi^+\\pi^-)$ with a kinematic edge ($0.28\\leq M(\\pi^+\\pi^-)<0.4116$ GeV/$c^2$), all $M(e^+e^-)$")
print("y ranges: r1 %s, all %s" % (y1, y2))

# per M(ee) class, model (b)
fig = plt.figure(figsize=(17.5, 13.6))
gs = fig.add_gridspec(3, 2, height_ratios=[2.25, 1, 0.78])
axs = np.array([[fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])], [None, None]], dtype=object)
axs[1, 0] = fig.add_subplot(gs[1, 0], sharex=axs[0, 0]); axs[1, 1] = fig.add_subplot(gs[1, 1], sharex=axs[0, 1])
axl = fig.add_subplot(gs[2, :])
tags = [("r1_E1_", "b"), ("r1_E2_", "b")]
ymax = ylim_of(tags); rlo, rhi = rlim_of(tags)
panel(axs[0, 0], axs[1, 0], "r1_E1_", "r1", "E1", "b", ymax, rlo, rhi, "(a) $M(e^+e^-)<10$ MeV/$c^2$ (E1): stage A extrapolated from E2")
panel(axs[0, 1], axs[1, 1], "r1_E2_", "r1", "E2", "b", ymax, rlo, rhi, "(b) $M(e^+e^-)\\geq10$ MeV/$c^2$ (E2): where stage A was fitted")
shared_legend(fig, axl)
fig.suptitle("SMASH with the measured scales, per $M(e^+e^-)$ class: S0, $0.30\\leq M(\\pi^+\\pi^-)<0.40$ GeV/$c^2$", fontsize=15)
fig.tight_layout(rect=(0, 0, 1, 0.97))
save(fig, "udatasim_u_overlay_mee")
print("mee y range %s, ratio %s" % (ymax, (rlo, rhi)))
