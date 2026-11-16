#!/usr/bin/env python3
# ============================================================================
# fix2_P1abc.py -- final review (prefix fix2, 2026-11-14): copy of vtx_P1abc.py with three changes in panel (c):
#   the total error of sigma is 38.3 % (the +-5 % period normalisation of Part II added to N; fix2_numbers.txt);
#   the arrow of the direct-decay route with the one-sided items starts at its weighted value 0.522 mb and ends at 0.588 mb
#   (0.58-0.59 mb; the per-class side-band value 0.621 mb would give 0.70 mb, written next to it);
#   the label no longer mixes the two definitions. Output: docs/latex/figs_all/fix2_P1abc.pdf/.png.
# Run: python3 -I plots/all_report/fix2_P1abc.py   (seconds). Original header follows.
# ----------------------------------------------------------------------------
# vtx_P1abc.py  --  Fig. P1 panels (a)-(c) of report_all after the measurement of the eta efficiency of the event-vertex cut
# (integration of 2026-11-14). Copy of the P1 part of fix_figs.py (review fix; same numbers, colours and layout); the only
# change is panel (c): the arrow of the one-sided items no longer contains the event-vertex factor x1.4-1.7 (measured:
# C_vtx = 1.006 +- 0.021, notes_all/ve_verify_note.md) and points to 0.48 mb (live time x sorter x cor_ok = 1.109-1.127,
# tables_all/vtx_numbers.txt); a second arrow shows the direct-decay route with the same items (0.58-0.70 mb).
# Output: docs/latex/figs_all/vtx_P1abc.pdf/.png (fix_P1abc.pdf is left unchanged).
# Run: python3 -I plots/all_report/vtx_P1abc.py   (seconds)
# ============================================================================
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FIG = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_all/"
TAB = "/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all/"
plt.rcParams.update({"font.size": 13, "axes.labelsize": 15, "axes.titlesize": 14, "legend.fontsize": 11, "xtick.labelsize": 12,
                     "ytick.labelsize": 12, "pdf.fonttype": 42})
PAL = dict(main="#1f3fbf", direct="#00a0a0", conv="#009900", fmain="#555555", b4="#ff6a00", edge="#8c564b", tag="#9467bd", pm="#d400d4",
           lit="#e8202a", smash="#ff8c00")

# ---------------------------------------------------------------- eta-scale estimates (PDG 2026; fix_numbers.txt)
# (label, value, stat, total, colour, marker, note)
SC = [("peak side bands $a_P$ (MAIN)", 0.4091, 0.0348, 0.4091 * math.sqrt(0.085 ** 2 + 0.355 ** 2 + 0.059 ** 2), PAL["main"], "o", "total incl. composition $\\pm32\\%$"),
      ("direct-decay route (C5), weighted", 0.497, 0.059, None, PAL["direct"], "D", "stat only"),
      ("$-$ side bands per $M_{ee}$ class, $\\rho=1$", 0.592, 0.139, None, PAL["direct"], "d", ""),
      ("$-$ X17 2D fit, $N_5/N_5^{\\mathrm{ref}}$", 0.476, 0.065, None, PAL["direct"], "d", ""),
      ("conversion route (C4), weighted", 0.237, 0.049, None, PAL["conv"], "s", "stat only"),
      ("F-main $a_\\eta$ (cross-check)", 0.3806, 0.0169, 0.3806 * math.sqrt(0.044 ** 2 + 0.555 ** 2 + 0.104 ** 2), PAL["fmain"], "o", ""),
      ("$f_\\eta$(edge) $=E_{\\mathrm{data}}/E_{\\mathrm{sim}}$", 0.3899, 0.0239, math.sqrt(0.0239 ** 2 + 0.039 ** 2 + 0.021 ** 2), PAL["edge"], "o", ""),
      ("$f_\\eta$(tag), emulated (UNVERIFIED)", 0.4520, 0.014, math.sqrt(0.014 ** 2 + 0.040 ** 2), PAL["tag"], "o", ""),
      ("$P_{\\mathrm{main}}$ $S_{\\mathrm{data}}/S_{\\mathrm{sim}}$, $r_{\\mathrm{res}}=1$", 0.3944, 0.0276, math.sqrt(0.0276 ** 2 + 0.0378 ** 2), PAL["pm"], "o", ""),
      ("$P_{\\mathrm{main}}$, $r_{\\mathrm{res}}=0.89$ (data width)", 0.3548, 0.0248, math.sqrt(0.0248 ** 2 + 0.0340 ** 2), PAL["pm"], "^", ""),
      ("bump B4 $a_B(E2)$: NOT DETERMINED", 0.3103, 0.018, None, PAL["b4"], "o", "syst 165 %")]

def draw_scales(ax, title=True):
    n = len(SC)
    for k, (lab, v, s, t, col, mk, note) in enumerate(SC):
        y = n - k
        if t is not None: ax.errorbar(v, y, xerr=t, fmt="none", ecolor=col, elinewidth=1.2, capsize=0)
        ax.errorbar(v, y, xerr=s, fmt=mk, color=col, ms=7, elinewidth=3.5, capsize=0, mfc=col if "NOT" not in lab else "white")
        if "NOT" in lab:
            ax.annotate("", xy=(0.95, y), xytext=(0.62, y), arrowprops=dict(arrowstyle="<->", color=col, lw=1.2))
            ax.text(0.785, y + 0.25, "syst 165 %", color=col, ha="center", fontsize=10)
    ax.axvline(0.4091, color=PAL["main"], ls=":", lw=1)
    ax.axvline(1.0, color="k", ls="--", lw=1); ax.text(0.99, n + 0.6, "SMASH", ha="right", fontsize=10)
    ax.set_yticks(range(n, 0, -1)); ax.set_yticklabels([s[0] for s in SC], fontsize=11)
    ax.set_xlim(0, 1.05); ax.set_ylim(0.3, n + 1.0)
    ax.set_xlabel("data / SMASH ($\\sigma_{\\mathrm{prod}}\\times$PDG 2026 BR), no $C_{\\mathrm{sec}}(1-f_\\eta)$")
    if title: ax.set_title("$\\eta$ scale estimates (thick: stat; thin: total)")

# ---------------------------------------------------------------- P1 (a)(b)(c)
fig = plt.figure(figsize=(14, 12.5))
gsP = fig.add_gridspec(2, 2, width_ratios=[1.35, 1.0], height_ratios=[1, 1.05], hspace=0.42, wspace=0.62)
axs = [fig.add_subplot(gsP[0, 0]), fig.add_subplot(gsP[0, 1]), fig.add_subplot(gsP[1, :])]
draw_scales(axs[0], title=False); axs[0].set_title("(a) $\\eta$ scale relative to SMASH")
ax = axs[1]
BR = [("BR$_G$ (F-gold, $M_{ee}\\geq20$)", 2.666, 0.371, math.sqrt(0.371 ** 2 + 3.383 ** 2 + 0.667 ** 2 + 0.09 ** 2)),
      ("BR$_A$ $=$ BR$_{\\mathrm{PDG}}\\,\\rho_d$ (F-main)", 3.540, 0.308, math.sqrt(0.308 ** 2 + 3.253 ** 2 + 0.609 ** 2 + 0.11 ** 2)),
      ("BR$_C$ (E1+E2 fit; withdrawn)", 3.6, 0.9, math.sqrt(0.9 ** 2 + 3.6 ** 2))]
ax.axvspan(2.76 - 0.16, 2.76 + 0.16, color="#bbbbbb", alpha=0.7, lw=0, label="PDG 2026: $2.76\\pm0.16$")
ax.axvline(2.68, color="k", ls="--", lw=1, label="PDG 2020 (simulation): 2.68")
for k, (lab, v, s, t) in enumerate(BR):
    y = 3 - k
    ax.errorbar(v, y, xerr=t, fmt="none", ecolor="k", elinewidth=1.2)
    ax.errorbar(v, y, xerr=s, fmt="o", color="k", ms=7, elinewidth=4, mfc="k" if k < 2 else "white")
ax.set_yticks([3, 2, 1]); ax.set_yticklabels([b[0] for b in BR], fontsize=11)
ax.set_ylim(0.4, 3.9); ax.set_xlim(-1.5, 8.5)
ax.set_xlabel("BR$(\\eta\\to\\pi^+\\pi^-e^+e^-)$ [$10^{-4}$]")
ax.set_title("(b) BR, cross-checks (closure failed)")
ax.legend(loc="upper right", fontsize=10, frameon=False)
ax.text(-1.3, 0.55, "total uncertainty 90-130 %: not a measurement", fontsize=10)
ax = axs[2]
pts = [("this work, peak side bands", 0.429, 0.036, 0.429 * 0.3826, PAL["main"], "o"),
       ("  direct-decay route (C5)", 0.522, 0.062, None, PAL["direct"], "D"),
       ("  conversion route (C4)", 0.249, 0.051, None, PAL["conv"], "s"),
       ("  F-main / edge / $P_{\\mathrm{main}}$ / tag", None, None, None, None, None),
       ("HADES 3.5 GeV, PYTHIA extrap.", 1.14, 0.20, None, PAL["lit"], "*"),
       ("HADES 3.5 GeV, UrQMD extrap.", 0.93, 0.14, None, PAL["lit"], "*"),
       ("same beam time, ECAL, prelim. (no errors)", None, None, None, None, None)]
n = len(pts)
for k, (lab, v, s, t, col, mk) in enumerate(pts):
    x = k + 1
    if v is None: continue
    if t is not None: ax.errorbar(x, v, yerr=t, fmt="none", ecolor=col, elinewidth=1.2)
    ax.errorbar(x, v, yerr=s, fmt=mk, color=col, ms=10 if mk == "*" else 7, elinewidth=3.5)
for v, mk, nm in ((0.399, "o", "F-main 0.399"), (0.409, "s", "edge 0.409"), (0.414, "^", "$P_{\\mathrm{main}}$ 0.414"), (0.474, "v", "tag 0.474")):
    ax.plot(4, v, mk, color="#555555", mfc="white", ms=7)
ax.text(4.2, 0.36, "F-main 0.399\nedge 0.409\n$P_{\\mathrm{main}}$ 0.414\ntag 0.474", fontsize=8.5, color="#555555", va="top")
ax.plot([7, 7], [1.15, 1.72], "*", color=PAL["lit"], mfc="white", ms=12)
ax.annotate("", xy=(1.25, 0.484), xytext=(1.25, 0.429), arrowprops=dict(arrowstyle="->", color=PAL["main"], lw=1.6))
ax.text(1.12, 1.62, "one-sided items (not applied):\nlive time, sorter, cor_ok: 0.48 mb;\nevent vertex (measured):\n$C_{\\mathrm{vtx}}=1.006\\pm0.021$, $\\leq$+4 % (95 % CL)",
        fontsize=9.5, color=PAL["main"], va="top")
ax.annotate("", xy=(2.25, 0.588), xytext=(2.25, 0.522), arrowprops=dict(arrowstyle="->", color=PAL["direct"], lw=1.6))
ax.text(2.35, 0.86, "with the one-sided items:\n0.58-0.59 mb (0.70 mb with\nthe per-class side-band value)", fontsize=9.5, color=PAL["direct"], va="top")
ax.annotate("", xy=(0.8, 0.383), xytext=(0.8, 0.429), arrowprops=dict(arrowstyle="->", color="#444444", lw=1.4))
ax.text(0.12, 0.27, "$r_{\\mathrm{res}}=0.89$\n(data width):\n0.383", fontsize=9, color="#444444")
ax.axhline(1.0694, color=PAL["smash"], ls="--", lw=1.4); ax.text(7.4, 1.085, "SMASH 1.069 mb", color=PAL["smash"], ha="right", fontsize=10)
ax.set_xticks(range(1, n + 1)); ax.set_xticklabels([p[0] for p in pts], rotation=35, ha="right", fontsize=10)
ax.set_xlim(0.0, n + 0.6); ax.set_ylim(0, 1.9)
ax.set_ylabel("$\\sigma(pp\\to\\eta X)$ [mb]")
ax.set_title("(c) $\\sigma(pp\\to\\eta X)$, $T=4.5$ GeV (3.5 GeV: published)")
ax.text(0.25, 1.82, "this work: PRELIMINARY, $\\pm$38.3 % (thin), stat (thick); routes: stat only", fontsize=9.5, va="top")
fig.savefig(FIG + "fix2_P1abc.pdf", bbox_inches="tight"); fig.savefig(FIG + "fix2_P1abc.png", dpi=70, bbox_inches="tight")

