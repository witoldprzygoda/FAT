#!/usr/bin/env python3
# ============================================================================
# fix_after_plain.py (2026-11-16, short-report figure agent for report_all_short.tex): COPY of fix_after.py; writes
# docs/latex/figs_all/fix_C2_plain.{pdf,png} and fix_C3_plain.{pdf,png}; its table goes to the scratch directory (TAB below), the
# file tables_all/fix_tableC3_before.txt is not touched. Only label strings in plain words (readability instruction I2):
# stage A -> background normalisation; CR-U -> 0.575-0.70 (the normalisation fit range); B / P / CR-L windows -> bump window /
# peak window / 0.5225-0.53; S1..S5 -> the M(pipi) ranges; S0, isBest -> pions in vertex fit; G_D, G_C, G_eta, C2-C5 and T12 ->
# physics names; B_A -> normalised background. The axes of fix_C2 keep the positions of the original (fixed after
# tight_layout). Same cached inputs, same drawing code, every number unchanged.
# Run: python3 -I plots/all_report/fix_after_plain.py   (seconds)
# ============================================================================
# ============================================================================
# fix_after.py  --  review fix (prefix fix): redrawn Figs. C2 and C3 (stage-A renormalisation) in the common palette of the
# report, with a "before" series (SMASH unrenormalised, the same E2 class and the same unit), and the "before, E2" column of
# Table C3.
# Inputs (read only): /hdd3/.../all_report/v5/data_S0.txt, cbrep_S0.txt, sim_S0.txt (plain dumps of hists_v1.root by v5_dump.C)
#                     /hdd3/.../all_report/fix/stageA_dump.txt (frozen stage A of WP3, dumped by fix_dump_stageA.C)
# Definitions: after  = (N - k CB)/(B_A + 0.4 eta),  k = k_E2 of stage A = 0.95434, B_A = frozen stage-A background;
#              before = (N - CB)/T12, T12 = SMASH T1 + T2 (all groups, raw simulated rows, no renormalisation), k = 1.
#              eta = C2 + C3 + C4 + C5 + G_eta at SMASH x PDG BR (PDG 2020, as in the simulation), raw rows.
# Outputs: docs/latex/figs_all/fix_C2.{pdf,png}, fix_C3.{pdf,png}; docs/latex/tables_all/fix_tableC3_before.txt
# Run: python3 -I plots/all_report/fix_after.py   (seconds)
# ============================================================================
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

V5 = "/hdd3/przygoda/HADES/FAT_out/all_report/v5/"
FIX = "/hdd3/przygoda/HADES/FAT_out/all_report/fix/"
FIG = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_all/"
TAB = "/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ6/figplain_allshort/out/after/"   # _plain: the table goes to scratch
L = 6079.2; NB2RAW = L * 63. / 64.; RAW2NB = (64. / 63.) / L
KE2 = 0.95434
PAL = dict(N="black", CB="red", NmCB="#1f3fbf", GD="#d0d0d0", GC="#707070", Geta="#c9a3f5", C2="#ff6a00", C3="#d400d4",
           C4="#009900", C5="#00cccc", Gvo="#a6f23a", X17="#cc0077", B="#2b6fd6", P="#e8202a", U="#2ca02c", K="#d400d4", Lc="#ff8c00")
plt.rcParams.update({"font.size": 13, "axes.labelsize": 15, "axes.titlesize": 14, "legend.fontsize": 11, "xtick.labelsize": 12,
                     "ytick.labelsize": 12, "pdf.fonttype": 42})

edges = [0.400 + 0.005 * i for i in range(25)] + [0.520 + 0.0025 * i for i in range(1, 21)] + [0.570 + 0.005 * i for i in range(1, 27)]
edges = np.array([round(x, 6) for x in edges]); cen = 0.5 * (edges[1:] + edges[:-1]); wid = edges[1:] - edges[:-1]

def rows(fn):
    d = {}
    for line in open(fn):
        p = line.split()
        if p: d[p[0]] = np.array([float(x) for x in p[1:]])
    return d
D = rows(V5 + "data_S0.txt"); S = rows(V5 + "sim_S0.txt"); A = rows(FIX + "stageA_dump.txt")
R = {}
cur = None
for line in open(V5 + "cbrep_S0.txt"):
    p = line.split()
    if not p: continue
    R.setdefault(p[0], []).append([float(x) for x in p[2:]])   # key, replica index, 70 values
R = {k: np.array(v) for k, v in R.items()}
SIGC = ["C2", "C3", "C4", "C5", "Geta"]
H = ["F10_S%dE2" % j for j in range(1, 6)]
def g(h):
    N = D["N_" + h]; vN = D["Nerr_" + h] ** 2; CB = D["CB_" + h]
    vCB = R["CBrep_" + h].var(axis=0, ddof=1) if ("CBrep_" + h) in R else np.zeros(70)
    comp = {c: S[c + "_" + h] * NB2RAW for c in SIGC + ["GD", "Gconv"]}
    eta = sum(comp[c] for c in SIGC)
    T12 = eta + comp["GD"] + comp["Gconv"]
    vT12 = sum((S[c + "err_" + h] * NB2RAW) ** 2 for c in SIGC + ["GD", "Gconv"])
    BA = A["BA_" + h]; eBA = np.sqrt(A["BAstat_" + h] ** 2 + A["BAmc_" + h] ** 2)
    GDs, GCs = A["GDsmooth_" + h], A["GCsmooth_" + h]
    den_ = GDs + 3.0627053 * GCs
    fGD = np.divide(GDs, den_, out=np.zeros_like(GDs), where=den_ > 0)
    return dict(N=N, vN=vN, CB=CB, vCB=vCB, comp=comp, eta=eta, T12=T12, vT12=vT12, BA=BA, eBA=eBA, BAgd=BA * fGD, BAgc=BA * (1 - fGD))
G = {h: g(h) for h in H}
tot = {k: sum(G[h][k] for h in H) for k in ("N", "vN", "CB", "vCB", "eta", "T12", "vT12", "BA", "BAgd", "BAgc")}
tot["eBA"] = np.sqrt(sum(G[h]["eBA"] ** 2 for h in H))   # upper estimate (slices correlated through b, kB)
tot["comp"] = {c: sum(G[h]["comp"][c] for h in H) for c in SIGC + ["GD", "Gconv"]}

def shade(ax, ylo, yhi):
    for (a, b, col, ls) in ((0.455, 0.5225, PAL["B"], "-"), (0.53364, 0.55416, PAL["P"], "-"), (0.575, 0.700, PAL["U"], "-"), (0.5225, 0.5300, PAL["Lc"], "--")):
        ax.axvspan(a, b, color=col, alpha=0.10, lw=0)

# ---------------------------------------------------------------- Fig. C2 (per slice)
SLT = ["$M_{\\pi\\pi}<0.34$ GeV/$c^2$", "$M_{\\pi\\pi}$ 0.34$-$0.37 GeV/$c^2$", "$M_{\\pi\\pi}$ 0.37$-$0.39 GeV/$c^2$", "$M_{\\pi\\pi}$ 0.39$-$0.43 GeV/$c^2$", "$M_{\\pi\\pi}\\geq0.43$ GeV/$c^2$"]   # _plain: slice titles
fig, axs = plt.subplots(1, 5, figsize=(15, 4.9), sharey=False)
for j, h in enumerate(H):
    ax = axs[j]; q = G[h]
    den = q["BA"] + 0.4 * q["eta"]
    den = np.where(den > 1.0, den, np.nan)      # bins below the slice threshold (no model content) are not drawn
    rat = (q["N"] - KE2 * q["CB"]) / den
    err = np.sqrt(q["vN"] + KE2 ** 2 * q["vCB"]) / den
    band = q["eBA"] / den
    ax.fill_between(edges, np.append(1 - band, 1 - band[-1]), np.append(1 + band, 1 + band[-1]), step="post", color="#bbbbbb", alpha=0.8, lw=0, label="normalisation error (stat $\\oplus$ MC)")
    infit = (edges[:-1] >= 0.575 - 1e-9) & np.isfinite(rat)
    outf = (edges[:-1] < 0.575 - 1e-9) & np.isfinite(rat)
    ax.errorbar(cen[infit], rat[infit], yerr=err[infit], fmt="o", ms=4.5, color="black", mfc="black", label="in the normalisation fit (0.575$-$0.70)")
    ax.errorbar(cen[outf], rat[outf], yerr=err[outf], fmt="s", ms=4.5, color="#1f3fbf", mfc="white", label="outside the fit")
    exp1 = (q["BA"] + q["eta"]) / den
    ax.step(edges, np.append(exp1, exp1[-1]), where="post", color=PAL["CB"], lw=1.8, label="expected if $\\eta$ at SMASH ($a_\\eta=1$)")
    ax.axhline(1, color="k", ls=":", lw=1)
    ok = np.isfinite(rat) & np.isfinite(err) & (err < 0.5)
    lo = min(np.nanmin((rat - err)[ok]), np.nanmin(exp1[ok])); hi = max(np.nanmax((rat + err)[ok]), np.nanmax(exp1[ok]))
    ax.set_ylim(lo - 0.05 * (hi - lo), hi + 0.25 * (hi - lo))
    shade(ax, 0, 1)
    ax.set_xlim(0.40, 0.70); ax.set_xlabel("$M_4$ [GeV/$c^2$]")
    ax.set_title(SLT[j], fontsize=12)
    if j == 0: ax.set_ylabel("$(N-k\\,\\mathrm{CB})\\,/\\,(\\mathrm{bkg}+0.4\\,\\eta)$")
hh, ll = axs[0].get_legend_handles_labels()
fig.legend(hh, ll + [], loc="lower center", ncol=4, frameon=False, fontsize=12)
fig.text(0.5, 0.005, "all panels M(e\u207ae\u207b) \u2265 10 MeV/c\u00b2; shaded: bump window (blue), 0.5225\u20130.53 (orange), peak window (red), 0.575\u20130.70 (green, the normalisation fit range)", ha="center", fontsize=11)   # _plain: plain text (no mathtext: same descent as the original line)
fig.tight_layout(rect=[0, 0.10, 1, 1])
for a_, x0_ in zip(axs, (0.0568518519, 0.2517277778, 0.4466037037, 0.6414796296, 0.8363555556)):   # _plain: original positions
    a_.set_position((x0_, 0.2570975057, 0.1418111111, 0.6507936508))
fig.savefig(FIG + "fix_C2_plain.pdf"); fig.savefig(FIG + "fix_C2_plain.png", dpi=90)

# ---------------------------------------------------------------- Fig. C3 (all slices, E2)
fig = plt.figure(figsize=(10, 9.5))
ax = fig.add_axes([0.13, 0.36, 0.84, 0.60]); ar = fig.add_axes([0.13, 0.08, 0.84, 0.26], sharex=ax)
u = RAW2NB / wid      # raw rows per bin -> nb per GeV/c2
layers = [("BAgd", "Dalitz and vector-meson pairs, renormalised", PAL["GD"]), ("BAgc", "non-$\\eta$ conversion pairs, renormalised", PAL["GC"]),
          ("Geta", "other $\\eta$ decays (SMASH)", PAL["Geta"]), ("C2", "bump, $\\pi^0$ Dalitz (SMASH)", PAL["C2"]),
          ("C3", "bump, conversion (SMASH)", PAL["C3"]), ("C4", "peak, $\\eta\\to\\pi^+\\pi^-\\gamma$ conv. (SMASH)", PAL["C4"]),
          ("C5", "peak, $\\eta\\to\\pi^+\\pi^-e^+e^-$ (SMASH)", PAL["C5"])]
base = np.zeros(70)
for key, lab, col in layers:
    v = tot[key] if key.startswith("BA") else tot["comp"][key]
    top = base + v
    ax.fill_between(edges, np.append(base * u, base[-1] * u[-1]), np.append(top * u, top[-1] * u[-1]), step="post", color=col, lw=0, label=lab)
    base = top
tot04 = (tot["BA"] + 0.4 * tot["eta"]) * u
ax.step(edges, np.append(tot04, tot04[-1]), where="post", color="k", ls="--", lw=1.6, label="total with all $\\eta\\times0.4$")
y = (tot["N"] - KE2 * tot["CB"]) * u; ey = np.sqrt(tot["vN"] + KE2 ** 2 * tot["vCB"]) * u
ax.errorbar(cen, y, yerr=ey, fmt="o", ms=4, color=PAL["NmCB"], label="data $N-k\\,\\mathrm{CB}$, $k=0.954$")
ax.set_xlim(0.40, 0.70); ax.set_ylim(0, 1.25 * max(base * u))
ax.set_ylabel("$d\\sigma/dM_4$ [nb/(GeV/$c^2$)]")
ax.axvspan(0.575, 0.700, ymin=0, ymax=0.02, color=PAL["U"])
ax.text(0.637, 0.035 * 1.25 * max(base * u), "normalisation region", color=PAL["U"], ha="center", fontsize=11)
ax.legend(loc="upper left", frameon=False, fontsize=10.5, ncol=1)
ax.set_title("after normalising the background: pions in vertex fit, $M(e^+e^-)\\geq10$ MeV/$c^2$; $\\times64/63$, $L=6079.2$ nb$^{-1}$", fontsize=12)
plt.setp(ax.get_xticklabels(), visible=False)
den4 = tot["BA"] + 0.4 * tot["eta"]
r_after = (tot["N"] - KE2 * tot["CB"]) / den4; e_after = np.sqrt(tot["vN"] + KE2 ** 2 * tot["vCB"]) / den4
r_after1 = (tot["N"] - KE2 * tot["CB"]) / (tot["BA"] + tot["eta"])
r_before = (tot["N"] - tot["CB"]) / tot["T12"]
e_before = r_before * np.sqrt((tot["vN"] + tot["vCB"]) / (tot["N"] - tot["CB"]) ** 2 + tot["vT12"] / tot["T12"] ** 2)
ar.errorbar(cen, r_before, yerr=e_before, fmt="o", ms=3.5, color="#888888", mfc="white", label="before: $(N-\\mathrm{CB})$/SMASH total, not renormalised")
ar.errorbar(cen, r_after, yerr=e_after, fmt="o", ms=3.5, color="black", label="after: $(N-k\\,\\mathrm{CB})$/(normalised background $+\\,0.4\\,\\eta$)")
ar.step(edges, np.append(r_after1, r_after1[-1]), where="post", color=PAL["CB"], lw=1.4, label="after, $\\eta$ at SMASH: /(normalised background $+\\,\\eta$)")
ar.axhline(1, color="k", ls=":", lw=1)
ar.set_ylim(0.55, 1.65); ar.set_xlabel("$M(\\pi^+\\pi^-e^+e^-)$ [GeV/$c^2$]"); ar.set_ylabel("data / model")
ar.legend(loc="upper left", frameon=False, fontsize=9.5, ncol=1)
fig.savefig(FIG + "fix_C3_plain.pdf"); fig.savefig(FIG + "fix_C3_plain.png", dpi=90)

# ---------------------------------------------------------------- Table C3: "before, E2" column
regions = [("R0' 0.400-0.455", 0.400, 0.455, 13786), ("R1 = B 0.4550-0.5225", 0.455, 0.5225, 53241), ("R2 0.5225-0.5336", 0.5225, 0.53364, 11471),
           ("R3 = P 0.5336-0.5542", 0.53364, 0.55416, 27099), ("R4 0.5542-0.565", 0.55416, 0.565, 15519), ("R5 0.565-0.700", 0.565, 0.700, 357941),
           ("CR-U 0.575-0.700", 0.575, 0.700, 342098), ("CR-L 0.5225-0.5300", 0.5225, 0.5300, 7519)]
def frac(a, b):
    return np.clip((np.minimum(edges[1:], b) - np.maximum(edges[:-1], a)) / wid, 0, 1)
out = open(TAB + "fix_tableC3_before.txt", "w")
out.write("# fix_tableC3_before.txt (fix_after.py): 'before' ratio in the regions of Table C3, E2 (M(ee) >= 10 MeV/c2), all slices, S0, raw rows.\n")
out.write("# before = (N - CB)/T12 (SMASH T1+T2, no renormalisation, k = 1); after = (N - k CB)/(B_A + 0.4 eta) (WP3 Table C3). F10 bins split\n")
out.write("# uniformly at fractional window edges (N_frac vs the exact N of Table C3 given as a check).\n")
out.write("# region | N_frac (exact) | CB | T12 | before +- err | after (recomputed) | after (WP3)\n")
wp3 = [0.8276, 0.9518, 0.9325, 0.9687, 0.9604, 0.9981, 0.9993, 0.9227]
res = []
for (nm, a, b, Nex), w3 in zip(regions, wp3):
    f = frac(a, b)
    N = (f * tot["N"]).sum(); CB = (f * tot["CB"]).sum(); T = (f * tot["T12"]).sum()
    vN = (f * f * tot["vN"]).sum(); vCB = (f * f * tot["vCB"]).sum(); vT = (f * f * tot["vT12"]).sum()
    rb = (N - CB) / T; eb = rb * math.sqrt((vN + vCB) / (N - CB) ** 2 + vT / T ** 2)
    BA = (f * tot["BA"]).sum(); eta = (f * tot["eta"]).sum()
    ra = (N - KE2 * CB) / (BA + 0.4 * eta)
    out.write("%-22s | %9.0f (%6d) | %8.1f | %10.1f | %.3f +- %.3f | %.3f | %.3f\n" % (nm, N, Nex, CB, T, rb, eb, ra, w3))
    res.append((nm, rb, eb, w3))
out.close()
for r in res: print("%-22s before %.3f +- %.3f   after (WP3) %.3f" % r)
