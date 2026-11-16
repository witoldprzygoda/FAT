#!/usr/bin/env python3
# ============================================================================
# fix_figs.py  --  review fix (prefix fix): figures that replace superseded panels of report_all, drawn from the published
# numbers (no refit):
#   fix_E4cd.pdf : Fig. E4 panels (c) displaced-window closure (from wp6/results.txt, labels clear of the markers, with the
#                  uncertainty of the rms) and (d) eta-scale estimates with the corrected errors, PDG 2026;
#   fix_P1abc.pdf: Fig. P1 panels (a) eta scale per estimator and route, (b) BR(eta -> pi+pi-e+e-) cross-checks against PDG 2026
#                  and 2020, (c) sigma(pp -> eta X) with the adopted total, the known-sign items and the published inclusive values;
#   fix_X2.pdf   : Fig. X2 with a narrowed y range, the variants in a ratio pad, both frozen medians named.
# Inputs: docs/latex/tables_all/fix_numbers.txt values (typed below with their source), /hdd3/.../wp6/results.txt,
#         tables_all/x17d_scan.txt, x17c_expected.txt, x17c_expected_asimov.txt.
# Run: python3 -I plots/all_report/fix_figs.py   (seconds)
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

# ---------------------------------------------------------------- E4 (c) + (d)
def disp_pulls():
    rows = [l.split() for l in open("/hdd3/przygoda/HADES/FAT_out/all_report/wp6/results.txt") if l.startswith("disp ")]
    c = np.array([float(r[2]) for r in rows])
    def col(name, which):   # pull after the k-th occurrence of "pull"
        out = []
        for r in rows:
            idx = [i for i, x in enumerate(r) if x == "pull"]
            out.append(float(r[idx[which] + 1]))
        return np.array(out)
    return c, col("pull", 0), col("pull", 1), col("pull", 2)

fig, axs = plt.subplots(1, 2, figsize=(14, 6.2), gridspec_kw=dict(width_ratios=[1, 1.25]))
ax = axs[0]
c, p_ref, p_nor, p_e2 = disp_pulls()
ax.axhspan(-1, 1, color="#dddddd", lw=0)
ax.axhline(0, color="k", ls=":", lw=1)
ax.plot(c - 0.0015, p_ref, "o", color="black", ms=8, label="all $M_{ee}$, window removed from CR-U, stage A refitted (rms %.2f)" % math.sqrt((p_ref ** 2).mean()))
ax.plot(c, p_nor, "o", color="#777777", mfc="white", ms=7, label="all $M_{ee}$, nominal stage A (rms %.2f)" % math.sqrt((p_nor ** 2).mean()))
ax.plot(c + 0.0015, p_e2, "s", color=PAL["main"], mfc="white", ms=7, label="$M_{ee}\\geq10$ MeV/$c^2$ only (rms %.2f)" % math.sqrt((p_e2 ** 2).mean()))
ax.set_xlim(0.568, 0.702); ax.set_ylim(-3.2, 7.4)
ax.set_xlabel("centre of the displaced window [GeV/$c^2$]"); ax.set_ylabel("pull $=(Y-Y_{\\mathrm{exp}})/\\sigma_{\\mathrm{stat}}$")
ax.set_title("(c) peak side-band estimator at 13 displaced windows")
ax.text(0.572, 7.2, "$\\sigma_{\\mathrm{tr}}=\\sigma_{\\mathrm{stat}}\\sqrt{\\mathrm{rms}^2-1}=239^{+254}_{-239}$ rows (4.7 % of $Y_{\\mathrm{peak}}$)\n"
        "windows overlap: $n_{\\mathrm{eff}}\\approx3.4$; non-overlapping subsets: rms 0.67-1.40\n"
        "$M_{ee}\\geq10$ only: rms 1.57, $\\sigma_{\\mathrm{tr}}(E2)=395$ rows", fontsize=10.5, va="top")
ax.legend(loc="upper left", fontsize=9.8, frameon=False, bbox_to_anchor=(0.0, 0.80))
ax2 = axs[1]
draw_scales(ax2, title=False)
ax2.set_title("(d) $\\eta$ scale estimates, PDG 2026 (thick: stat; thin: total)")
fig.tight_layout()
fig.savefig(FIG + "fix_E4cd.pdf"); fig.savefig(FIG + "fix_E4cd.png", dpi=80)

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
pts = [("this work, peak side bands", 0.429, 0.036, 0.163, PAL["main"], "o"),
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
ax.annotate("", xy=(1.25, 0.82), xytext=(1.25, 0.429), arrowprops=dict(arrowstyle="->", color=PAL["main"], lw=1.6))
ax.text(1.35, 0.80, "one-sided items (not applied):\nvertex efficiency $\\times$1.4-1.7,\nlive time, sorter, cor_ok: 0.67-0.82 mb",
        fontsize=9.5, color=PAL["main"], va="top")
ax.annotate("", xy=(0.8, 0.383), xytext=(0.8, 0.429), arrowprops=dict(arrowstyle="->", color="#444444", lw=1.4))
ax.text(0.12, 0.27, "$r_{\\mathrm{res}}=0.89$\n(data width):\n0.383", fontsize=9, color="#444444")
ax.axhline(1.0694, color=PAL["smash"], ls="--", lw=1.4); ax.text(7.4, 1.085, "SMASH 1.069 mb", color=PAL["smash"], ha="right", fontsize=10)
ax.set_xticks(range(1, n + 1)); ax.set_xticklabels([p[0] for p in pts], rotation=35, ha="right", fontsize=10)
ax.set_xlim(0.0, n + 0.6); ax.set_ylim(0, 1.9)
ax.set_ylabel("$\\sigma(pp\\to\\eta X)$ [mb]")
ax.set_title("(c) $\\sigma(pp\\to\\eta X)$, $T=4.5$ GeV (3.5 GeV: published)")
ax.text(0.25, 1.82, "this work: PRELIMINARY, $\\pm$37.9 % (thin), stat (thick); routes: stat only", fontsize=9.5, va="top")
fig.savefig(FIG + "fix_P1abc.pdf", bbox_inches="tight"); fig.savefig(FIG + "fix_P1abc.png", dpi=70, bbox_inches="tight")

# ---------------------------------------------------------------- X2
def scan():
    d = []
    for l in open(TAB + "x17d_scan.txt"):
        if l.startswith("#") or not l.strip(): continue
        p = l.replace("|", " ").split()
        d.append([float(x) for x in p])
    d = np.array(d)
    return dict(m=d[:, 0], ul1=d[:, 2], ul2=d[:, 3], ul3=d[:, 4], cos2=d[:, 11])
def asimov():
    d = []
    for l in open(TAB + "x17c_expected_asimov.txt"):
        if l.startswith("#") or not l.strip(): continue
        p = l.replace("|", " ").split()
        d.append([float(x) for x in p[:7]])
    d = np.array(d)
    return dict(m=d[:, 0], m2=d[:, 2], m1=d[:, 3], med=d[:, 4], p1=d[:, 5], p2=d[:, 6])
def toys():
    t = []
    for l in open(TAB + "x17c_expected.txt"):
        p = l.split()
        if len(p) >= 9 and p[1] == "1" and p[2] == "toys":
            t.append([float(p[0]), float(p[5]), float(p[6]), float(p[7])])
    return np.array(t)
S = scan(); A = asimov(); T = toys()
fig = plt.figure(figsize=(11, 9.5))
ax = fig.add_axes([0.12, 0.40, 0.84, 0.55]); ar = fig.add_axes([0.12, 0.08, 0.84, 0.30], sharex=ax)
ax.fill_between(A["m"], A["m2"], A["p2"], color="#ffd23f", lw=0, label="frozen expectation $\\pm2\\sigma$ (Asimov)")
ax.fill_between(A["m"], A["m1"], A["p1"], color="#3cc23c", lw=0, label="frozen expectation $\\pm1\\sigma$ (Asimov)")
ax.plot(A["m"], A["med"], "k--", lw=1.6, label="frozen Asimov median (2.76 at 17 MeV/$c^2$)")
ax.errorbar(T[:, 0], T[:, 2], yerr=[T[:, 2] - T[:, 1], T[:, 3] - T[:, 2]], fmt="o", color="#006400", mfc="white", ms=7,
            label="frozen toy median, 16-84 % (2.60 at 17 MeV/$c^2$)")
ax.plot(S["m"], S["ul1"], "k-", lw=2.4, marker=".", ms=7, label="observed, normalisation (1), PRIMARY")
tc = [(16.7, 3.06e-5 * 0.958), (16.85, 3.065e-5 * 1.002), (17.0, 3.06e-5 * 0.983)]
ax.plot([x for x, _ in tc], [y for _, y in tc], "*", color="#e8202a", ms=13, label="observed (1), toy CL$_s$")
ax.set_yscale("log"); ax.set_ylim(7e-6, 1.3e-4); ax.set_xlim(9.5, 40.5)
ax.set_ylabel("90 % CL$_s$ UL on BR$(\\eta\\to\\pi^+\\pi^-X)$BR$(X\\to e^+e^-)$")
ax.legend(loc="upper right", fontsize=10.5, frameon=False, ncol=1)
ax.text(10.0, 8.2e-6, "S0, isBest $\\geq0$, X2D; $L=6079.2$ nb$^{-1}$; BR$_{\\mathrm{dir}}=2.68\\times10^{-4}$ (PDG 2020, frozen); spin-0 hypothesis", fontsize=10)
plt.setp(ax.get_xticklabels(), visible=False)
med = np.interp(S["m"], A["m"], A["med"])
ar.fill_between(A["m"], A["m2"] / A["med"], A["p2"] / A["med"], color="#ffd23f", lw=0)
ar.fill_between(A["m"], A["m1"] / A["med"], A["p1"] / A["med"], color="#3cc23c", lw=0)
ar.axhline(1, color="k", ls="--", lw=1.2)
ar.plot(S["m"], S["ul1"] / med, "k-", lw=2.2, label="(1) $N_5$ fitted (primary)")
ar.plot(S["m"], S["ul2"] / med, "--", color="#1f3fbf", lw=1.8, label="(2) $a_\\eta=0.38$ (spec. fall-back primary)")
ar.plot(S["m"], S["ul3"] / med, ":", color="#d400d4", lw=2.2, label="(3) peak $a_P=0.41$")
ar.plot(S["m"], S["cos2"] / med, "-.", color="#e8202a", lw=1.8, label="(1), $e^+$ angle $1+\\cos^2\\theta$ (vector-like)")
ar.set_ylim(0.3, 2.6); ar.set_xlabel("$m_X$ [MeV/$c^2$]"); ar.set_ylabel("UL / frozen Asimov median")
ar.legend(loc="upper right", fontsize=9.5, frameon=False, ncol=2)
fig.savefig(FIG + "fix_X2.pdf"); fig.savefig(FIG + "fix_X2.png", dpi=80)
print("done")
