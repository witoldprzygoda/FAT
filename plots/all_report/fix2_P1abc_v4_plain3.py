#!/usr/bin/env python3
# ============================================================================
# fix2_P1abc_v4_plain3.py (2026-11-16, short-report figure agent, second pass, for report_all_short.tex): copy of
# fix2_P1abc_v4_plain.py with the last fit-parameter tags of panels (a) and (b) also in words (instruction I2):
#   "peak side bands a_P (MAIN)" -> "peak side bands (MAIN)"; "per M_ee class, rho = 1" -> "per M_ee class, PDG rate";
#   "template fit a_eta" -> "template fit"; "BR_PDG rho_d, template fit" -> "BR, direct-decay scale"; "bump side band" (the
#   background-subtracted bump area, B4) -> "bump area". (fix2_P1abc_v4_plain2.py/.pdf: a first try whose longer labels were cut at the
#   left edge; not to be used.) Every number, colour, position and the page size are unchanged.
# Output: docs/latex/figs_all/fix2_P1abc_v4_plain3.pdf/.png.   Run: python3 -I plots/all_report/fix2_P1abc_v4_plain3.py   (seconds)
# ============================================================================
# fix2_P1abc_v4_plain.py (2026-11-16, short-report figure agent for report_all_short.tex): copy of fix2_P1abc_v4.py for the SHORT
# report (readability instruction I2: no internal labels inside the figure). ONLY label strings change, every number, colour, position
# and the page size stay as in fix2_P1abc_v4.pdf (the bounding box is fixed to the one of the original instead of "tight"):
#   (C5)/(C4) routes -> "direct decays alone" / "converted photons alone"; F-main -> template fit; P_main -> fixed-Gaussian fit;
#   r_res = 1 / 0.89 -> SMASH width / SMASH width x 0.89; "C2 angle corrected" -> "angle corrected"; B4 a_B(E2) -> bump side band;
#   BR_G (F-gold) / BR_A (F-main) / BR_C (E1+E2) -> BR, template fit, M_ee >= 20 / BR_PDG rho_d, template fit / BR, all-M_ee fit;
#   "routes" -> "parts of the peak"; cor_ok -> corr. momenta; N5/N5ref -> direct decays. Labels were shortened until they fit inside the original bounding box.
# Output: docs/latex/figs_all/fix2_P1abc_v4_plain.pdf/.png.   Run: python3 -I plots/all_report/fix2_P1abc_v4_plain.py   (seconds)
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
# ----------------------------------------------------------------------------
# _v3 copy (2026-11-15, FIGURE AGENT C of the fourth integration): ruling D40 (option (a): the peak-tail family t = 1 +- 0.5
# is counted; sigma(peak) and a_P +-8.75 %). Three changes, nothing else:
#   panel (a): the total of a_P gains the peak-tail term: syst sqrt(35.5^2 + 8.75^2) = 36.6 % (p3w_tableS1_v3.tex);
#   panel (c): the total of sigma(peak) is 39.2 % instead of 38.3 % (0.429 x 0.3925 = 0.168 mb; sigma_headline_v3.tex),
#              and the note on the panel says +-39.2 %.
# D38 does not enter this figure (the peak has no C2; 'edge 0.409' is the uncorrected per-slice edge, unchanged).
# Output: docs/latex/figs_all/fix2_P1abc_v3.pdf/.png.
# Run: python3 -I plots/all_report/fix2_P1abc_v3.py   (seconds)
# ----------------------------------------------------------------------------
# _v4 copy (2026-11-15, figure agent of the fifth integration): the bump-edge numbers brought to the final state (D35, D37,
# D38, D40); every other number of the three panels was checked against report_all.tex and is unchanged except one rounding.
#   panel (a): the row "f_eta(edge) = E_data/E_sim" (0.3899 +- 0.0239, total 0.050: the per-slice edge of WP8 without the C2
#              angle correction and with its old stat/MC errors) is replaced by two rows, both with the C2 angle corrected by the
#              per-count QED weights (D38), PDG 2026 (/1.00436, the convention of the panel and of bedge_scales_v3):
#                f_eta(u-edge), the bump measurement (D37): 0.42574 +- 0.01989 +- 0.03387 +- 0.03002 (total 11.6 %)
#                    -> 0.4239, stat 0.0198, total 0.0492   (edgeu_v3/edgeu_numbers_v3.json; report 0.426 +- 0.020 +- 0.034 +- 0.030)
#                f_eta(edge) per slice (D35, cross-check): 0.42345 +- 0.02275 +- 0.03806 +- 0.04803 (total 15.4 %)
#                    -> 0.4216, stat 0.0227, total 0.0651   (bedge_v3/bedge_numbers_v3.json fc_rec; report 0.423 +- 0.023 +- 0.038 +- 0.048)
#              the ECAL tag row reads "emulated (verification only)" instead of "(UNVERIFIED)" (report: acceptance verified,
#              photon-energy scale not; verification only); the F-main total uses the Table S1 syst 55.4 % (was 55.5 %).
#   panel (c): the cross-check group shows the bump through its edge with the C2 angle corrected: u-edge 0.444 mb and per-slice
#              edge 0.442 mb (Table bedge-sigma) instead of the uncorrected per-slice 0.409 mb; tick label names both edges.
# Output: docs/latex/figs_all/fix2_P1abc_v4.pdf/.png.
# Run: python3 -I plots/all_report/fix2_P1abc_v4.py   (seconds)
# ============================================================================
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
# _plain: the bounding boxes of the original fix2_P1abc_v4.pdf/.png (tight box + 0.1 in), so that the page size is unchanged
BB_PDF = Bbox.from_extents(-1.291667 - 0.1, -0.537430 - 0.1, 12.600000 + 0.1, 11.236111 + 0.1)
BB_PNG = Bbox.from_extents(-1.303376 - 0.1, -0.547682 - 0.1, 12.600000 + 0.1, 11.269048 + 0.1)

FIG = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_all/"
TAB = "/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all/"
plt.rcParams.update({"font.size": 13, "axes.labelsize": 15, "axes.titlesize": 14, "legend.fontsize": 11, "xtick.labelsize": 12,
                     "ytick.labelsize": 12, "pdf.fonttype": 42})
PAL = dict(main="#1f3fbf", direct="#00a0a0", conv="#009900", fmain="#555555", b4="#ff6a00", edge="#8c564b", tag="#9467bd", pm="#d400d4",
           lit="#e8202a", smash="#ff8c00")
PDG26 = 1.00436   # _v4: PDG 2020 -> PDG 2026 for the C2+C3-based scales (as bedge_figs_v3.py and the f_eta(tag) row)

# ---------------------------------------------------------------- eta-scale estimates (PDG 2026; fix_numbers.txt)
# (label, value, stat, total, colour, marker, note)
SC = [("peak side bands (MAIN)", 0.4091, 0.0348, 0.4091 * math.sqrt(0.085 ** 2 + 0.355 ** 2 + 0.0875 ** 2 + 0.059 ** 2), PAL["main"], "o", "total incl. composition $\\pm32\\%$"),
      ("direct decays alone, weighted", 0.497, 0.059, None, PAL["direct"], "D", "stat only"),
      ("$-$ side bands per $M_{ee}$ class, PDG rate", 0.592, 0.139, None, PAL["direct"], "d", ""),
      ("$-$ X17 2D fit, direct decays", 0.476, 0.065, None, PAL["direct"], "d", ""),
      ("converted photons alone, weighted", 0.237, 0.049, None, PAL["conv"], "s", "stat only"),
      ("template fit (cross-check)", 0.3806, 0.0169, 0.3806 * math.sqrt(0.044 ** 2 + 0.554 ** 2 + 0.104 ** 2), PAL["fmain"], "o", ""),
      # _v4: the bump through its edge, C2 angle corrected (D38), PDG 2026 (/1.00436); value, stat, total = stat (+) syst (+) MC
      ("$f_\\eta$($u$-edge), angle corrected (bump)", 0.42573810 / PDG26, 0.01988893 / PDG26,
       math.sqrt(0.01988893 ** 2 + 0.03386803 ** 2 + 0.03001558 ** 2) / PDG26, PAL["edge"], "o", ""),
      ("$f_\\eta$(edge) per slice, angle corrected", 0.423451 / PDG26, 0.022752 / PDG26,
       math.sqrt(0.022752 ** 2 + 0.03805759 ** 2 + 0.04802620 ** 2) / PDG26, PAL["edge"], "s", ""),
      ("$f_\\eta$(tag), emulated (verification only)", 0.4520, 0.014, math.sqrt(0.014 ** 2 + 0.040 ** 2), PAL["tag"], "o", ""),
      ("Gaussian fit $S_{\\mathrm{data}}/S_{\\mathrm{sim}}$, SMASH width", 0.3944, 0.0276, math.sqrt(0.0276 ** 2 + 0.0378 ** 2), PAL["pm"], "o", ""),
      ("Gaussian fit, SMASH width $\\times$ 0.89", 0.3548, 0.0248, math.sqrt(0.0248 ** 2 + 0.0340 ** 2), PAL["pm"], "^", ""),
      ("bump area: NOT DETERMINED", 0.3103, 0.018, None, PAL["b4"], "o", "syst 165 %")]

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
BR = [("BR, template fit, $M_{ee}\\geq20$", 2.666, 0.371, math.sqrt(0.371 ** 2 + 3.383 ** 2 + 0.667 ** 2 + 0.09 ** 2)),
      ("BR, direct-decay scale", 3.540, 0.308, math.sqrt(0.308 ** 2 + 3.253 ** 2 + 0.609 ** 2 + 0.11 ** 2)),
      ("BR, all-$M_{ee}$ fit (withdrawn)", 3.6, 0.9, math.sqrt(0.9 ** 2 + 3.6 ** 2))]
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
pts = [("this work, peak side bands", 0.429, 0.036, 0.429 * 0.3925, PAL["main"], "o"),
       ("  direct decays alone", 0.522, 0.062, None, PAL["direct"], "D"),
       ("  converted photons alone", 0.249, 0.051, None, PAL["conv"], "s"),
       ("  template fit / edges / Gaussian fit / tag", None, None, None, None, None),
       ("HADES 3.5 GeV, PYTHIA extrap.", 1.14, 0.20, None, PAL["lit"], "*"),
       ("HADES 3.5 GeV, UrQMD extrap.", 0.93, 0.14, None, PAL["lit"], "*"),
       ("same beam time, ECAL, prelim. (no errors)", None, None, None, None, None)]
n = len(pts)
for k, (lab, v, s, t, col, mk) in enumerate(pts):
    x = k + 1
    if v is None: continue
    if t is not None: ax.errorbar(x, v, yerr=t, fmt="none", ecolor=col, elinewidth=1.2)
    ax.errorbar(x, v, yerr=s, fmt=mk, color=col, ms=10 if mk == "*" else 7, elinewidth=3.5)
for v, mk, nm in ((0.399, "o", "F-main 0.399"), (0.444, "s", "$u$-edge 0.444"), (0.442, "D", "edge per slice 0.442"),
                  (0.414, "^", "$P_{\\mathrm{main}}$ 0.414"), (0.474, "v", "tag 0.474")):
    ax.plot(4, v, mk, color="#555555", mfc="white", ms=7)
ax.text(4.2, 0.36, "template fit 0.399\nGaussian fit 0.414\nedge per slice 0.442 (angle corr.)\n$u$-edge 0.444 (bump, angle corr.)\ntag 0.474",
        fontsize=8.5, color="#555555", va="top")
ax.plot([7, 7], [1.15, 1.72], "*", color=PAL["lit"], mfc="white", ms=12)
ax.annotate("", xy=(1.25, 0.484), xytext=(1.25, 0.429), arrowprops=dict(arrowstyle="->", color=PAL["main"], lw=1.6))
ax.text(1.12, 1.62, "one-sided items (not applied):\nlive time, sorter, corr. momenta: 0.48 mb;\nevent vertex (measured):\n$C_{\\mathrm{vtx}}=1.006\\pm0.021$, $\\leq$+4 % (95 % CL)",
        fontsize=9.5, color=PAL["main"], va="top")
ax.annotate("", xy=(2.25, 0.588), xytext=(2.25, 0.522), arrowprops=dict(arrowstyle="->", color=PAL["direct"], lw=1.6))
ax.text(2.35, 0.86, "with the one-sided items:\n0.58-0.59 mb (0.70 mb with\nthe per-class side-band value)", fontsize=9.5, color=PAL["direct"], va="top")
ax.annotate("", xy=(0.8, 0.383), xytext=(0.8, 0.429), arrowprops=dict(arrowstyle="->", color="#444444", lw=1.4))
ax.text(0.12, 0.27, "SMASH width\n$\\times$ 0.89:\n0.383", fontsize=9, color="#444444")
ax.axhline(1.0694, color=PAL["smash"], ls="--", lw=1.4); ax.text(7.4, 1.085, "SMASH 1.069 mb", color=PAL["smash"], ha="right", fontsize=10)
ax.set_xticks(range(1, n + 1)); ax.set_xticklabels([p[0] for p in pts], rotation=35, ha="right", fontsize=10)
ax.set_xlim(0.0, n + 0.6); ax.set_ylim(0, 1.9)
ax.set_ylabel("$\\sigma(pp\\to\\eta X)$ [mb]")
ax.set_title("(c) $\\sigma(pp\\to\\eta X)$, $T=4.5$ GeV (3.5 GeV: published)")
ax.text(0.25, 1.82, "this work: PRELIMINARY, $\\pm$39.2 % (thin), stat (thick); parts of the peak: stat only", fontsize=9.5, va="top")
fig.savefig(FIG + "fix2_P1abc_v4_plain3.pdf", bbox_inches=BB_PDF); fig.savefig(FIG + "fix2_P1abc_v4_plain3.png", dpi=70, bbox_inches=BB_PNG)
for lab, v, s, t, col, mk, note in SC:   # _v4: the drawn numbers, for the check against report_all.tex
    print("%-62s value %.4f  stat %.4f  total %s" % (lab, v, s, "%.4f (%.1f %%)" % (t, 100 * t / v) if t is not None else "--"))

