#!/usr/bin/env python3
# ============================================================================
# fix_S1.py  --  review fix (prefix fix): Table S1 and Fig. S1 of WP7 brought up to date, from the published family values
# (tables_all/wp7_tableS1.txt) with three replacements:
#   (1) family S12 = the completed S12 of WP8 (tables_all/wp8_S12.txt), and for sigma_eta the composition of the peak taken from
#       the data (fix_numbers.txt: +-31.8 %, half the difference of the direct-decay and the conversion routes), family maximum;
#   (2) family S8 widened to the observed range of g (0.47-2.6) by linear extrapolation of the WP7 "g free" variant where that is
#       meaningful (Y_bump, Y_peak, sigma_eta, Y5, Y4, a_eta); other columns unchanged;
#   (3) values with PDG 2026 (BR_G, BR_A, DR, sigma_eta 0.429 mb); the X17 limit stays with the frozen PDG 2020 BR_dir.
# Totals: families in quadrature (rounded family values; differences to WP7 below 0.1 point). Rows "without R1" and "Barlow with
# stat only": WP7 values, recomputed only for sigma_eta (with the composition term).
# Output: docs/latex/tables_all/fix_tableS1.{txt,tex}, docs/latex/figs_all/fix_S1.{pdf,png}
# Run: python3 -I plots/all_report/fix_S1.py
# ============================================================================
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TAB = "/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all/"
FIG = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_all/"
Q = ["Y_bump", "Y_peak", "Y5", "Y4", "a_eta", "BR_G", "BR_A", "kappa", "DR", "sigma_eta"]
FAM = [("S1", "CB"), ("S2", "CR definition"), ("S3", "background model"), ("S4", "background template"), ("S5", "range, binning"),
       ("S6", "mass transfer"), ("S7", "physics weights"), ("S8", "G_eta"), ("S9", "signal templates"), ("S10", "closure"),
       ("S11", "selection"), ("S12", "efficiency model")]
rows = {}
for l in open(TAB + "wp7_tableS1.txt"):
    p = l.split()
    if not p: continue
    for f, nm in FAM:
        if l.startswith(f + " " + nm):
            vals = p[len(nm.split()) + 1:len(nm.split()) + 1 + 10]
            rows[f] = [float(v) if v != "n/a" else float("nan") for v in vals]
    if l.startswith("stat "): rows["stat"] = [float(v) for v in p[1:11]]
    if l.startswith("MC (template"): rows["MC"] = [float(v) for v in p[3:13]]
    if l.startswith("  without rule R1"): rows["noR1"] = [float(v) for v in p[3:13]]
    if l.startswith("  Barlow with stat only"): rows["barlow"] = [float(v) for v in p[4:14]]
assert all(len(rows[f]) == 10 for f, _ in FAM), rows
fam = {f: list(rows[f]) for f, _ in FAM}
# (1) S12 of WP8 and the composition for sigma_eta
S12 = dict(Y_bump=3.4, Y_peak=0.1, Y5=1.3, Y4=8.4, a_eta=12.2, BR_G=8.6, BR_A=12.1, kappa=3.4, DR=10.2, sigma_eta=31.8)
for k, q in enumerate(Q): fam["S12"][k] = S12[q]
# (2) S8 widened (fix_numbers.txt section 6; Y5 from the WP7 'g free' shift +0.68 % -> -1.0 % at g = 0.47)
S8 = dict(Y_bump=7.7, Y_peak=1.0, sigma_eta=1.0, a_eta=9.1, Y4=7.4, Y5=1.0)
for k, q in enumerate(Q):
    if q in S8: fam["S8"][k] = S8[q]
tot = []
for k in range(10):
    s = math.sqrt(sum(fam[f][k] ** 2 for f, _ in FAM if not math.isnan(fam[f][k])))
    tot.append(s)
noR1 = list(rows["noR1"]); bar = list(rows["barlow"])
noR1[9] = math.sqrt(18.6 ** 2 - 2.5 ** 2 + 31.8 ** 2); bar[9] = math.sqrt(17.7 ** 2 - 2.5 ** 2 + 31.8 ** 2)
N = ["--", "--", "--", "--", "+7.6/-7.6", "+3.2/-3.2", "+3.2/-3.2", "--", "+4.5/-4.5", "+8.2/-8.6"]
VAL = ["16389", "5088", "3198", "1847", "0.3823", "2.67e-4", "3.54e-4", "0.832", "1.373", "0.429 mb"]
VALT = ["16389", "5088", "3198", "1847", "0.3823", "$2.67\\cdot10^{-4}$", "$3.54\\cdot10^{-4}$", "0.832", "1.373", "0.429\\,mb"]
HEAD = ["$Y_\\mathrm{bump}$", "$Y_\\mathrm{peak}$", "$Y_5$", "$Y_4$", "$a_\\eta$", "BR$_G$", "BR$_A$", "$\\kappa$", "DR", "$\\sigma_\\eta$"]
def f1(x): return "--" if math.isnan(x) else "%.1f" % x
with open(TAB + "fix_tableS1.txt", "w") as o:
    o.write("# fix_tableS1.txt (fix_S1.py): Table S1 of WP7 updated: S12 of WP8 (+ composition of the peak from the data for sigma_eta),\n")
    o.write("# S8 widened to g = 0.47-2.6 (Y_bump, Y_peak, sigma_eta, Y5, Y4, a_eta), PDG 2026 values. Per cent of the value.\n")
    o.write("%-26s" % "" + "".join("%11s" % q for q in Q) + "\n")
    o.write("%-26s" % "value" + "".join("%11s" % v for v in VAL) + "\n")
    o.write("%-26s" % "stat" + "".join("%11.1f" % v for v in rows["stat"]) + "\n")
    for f, nm in FAM: o.write("%-26s" % (f + " " + nm) + "".join("%11s" % f1(v) for v in fam[f]) + "\n")
    o.write("%-26s" % "total syst" + "".join("%11.1f" % v for v in tot) + "\n")
    o.write("%-26s" % "  without rule R1" + "".join("%11.1f" % v for v in noR1) + "\n")
    o.write("%-26s" % "  Barlow, stat only" + "".join("%11.1f" % v for v in bar) + "\n")
    o.write("%-26s" % "MC" + "".join("%11.1f" % v for v in rows["MC"]) + "\n")
    o.write("%-26s" % "N" + "".join("%11s" % v for v in N) + "\n")
with open(TAB + "fix_tableS1.tex", "w") as o:
    o.write("% fix_tableS1.tex -- generated by plots/all_report/fix_S1.py (review fix of WP7 Table S1); relative uncertainties in per cent\n")
    o.write("\\begin{tabular}{l" + "r" * 10 + "}\n\\toprule\n & " + " & ".join(HEAD) + " \\\\\n\\midrule\n")
    o.write("value & " + " & ".join(VALT) + " \\\\\n")
    o.write("stat & " + " & ".join("%.1f" % v for v in rows["stat"]) + " \\\\\n\\midrule\n")
    for f, nm in FAM:
        nm2 = nm.replace("G_eta", "$G_\\eta$")
        cells = []
        for k, v in enumerate(fam[f]):
            t = f1(v)
            if (f == "S12" and Q[k] == "sigma_eta") or (f == "S8" and Q[k] in S8) or (f == "S12"): t = "\\textit{" + t + "}" if t != "--" else t
            cells.append(t)
        o.write("%s %s & " % (f, nm2) + " & ".join(cells) + " \\\\\n")
    o.write("\\midrule\n\\textbf{total syst} & " + " & ".join("\\textbf{%.1f}" % v if k in (1, 9) else "%.1f" % v for k, v in enumerate(tot)) + " \\\\\n")
    o.write("\\quad without rule R1 & " + " & ".join("%.1f" % v for v in noR1) + " \\\\\n")
    o.write("\\quad Barlow, stat only & " + " & ".join("%.1f" % v for v in bar) + " \\\\\n")
    o.write("MC & " + " & ".join("%.1f" % v for v in rows["MC"]) + " \\\\\n")
    o.write("N & " + " & ".join("--" if v == "--" else "$^{+%s}_{-%s}$" % tuple(v[1:].split("/-")) for v in N) + " \\\\\n")
    o.write("\\bottomrule\n\\end{tabular}\n")
print(open(TAB + "fix_tableS1.txt").read())

# ---------------------------------------------------------------- Fig. S1
plt.rcParams.update({"font.size": 13, "axes.labelsize": 15, "axes.titlesize": 15, "xtick.labelsize": 13, "ytick.labelsize": 13, "pdf.fonttype": 42})
TEXH = [r"$Y_\mathrm{bump}$ (B4, extrap.)", r"$Y_\mathrm{peak}$ (side bands)", r"$Y_5$ (F-main)", r"$Y_4$ (F-main)", r"$a_\eta$ (F-main)",
        r"BR$_G$ (F-gold)", r"BR$_A$ (F-main)", r"$\kappa$ (F-main)", r"DR (side bands/B4)", r"$\sigma_\eta$ (side bands)"]
fig = plt.figure(figsize=(16, 11))
gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 1.0], wspace=0.32)
ax = fig.add_subplot(gs[0, 0])
rowsH = [f for f, _ in FAM] + ["total syst", "stat", "MC"]
M = np.zeros((len(rowsH), 10))
for j in range(10):
    for i, (f, _) in enumerate(FAM): M[i, j] = fam[f][j]
    M[len(FAM), j] = tot[j]; M[len(FAM) + 1, j] = rows["stat"][j]; M[len(FAM) + 2, j] = rows["MC"][j]
Mc = np.where(np.isnan(M), 0, M)
im = ax.imshow(np.log10(np.clip(Mc, 0.1, 300)), cmap="YlOrRd", aspect="auto", vmin=-1, vmax=np.log10(300))
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        t = "n/a" if np.isnan(M[i, j]) else ("%.0f" % M[i, j] if M[i, j] >= 10 else "%.1f" % M[i, j])
        ax.text(j, i, t, ha="center", va="center", fontsize=10.5, color="black")
ax.set_xticks(range(10)); ax.set_xticklabels(TEXH, rotation=45, ha="right")
ax.set_yticks(range(len(rowsH)))
ax.set_yticklabels(["%s %s" % (f, nm.replace("G_eta", r"$G_\eta$")) for f, nm in FAM] + ["total syst", "stat", "MC"])
ax.axhline(len(FAM) - 0.5, color="k", lw=1.5)
ax.set_title("(a) relative uncertainty [%] per family")
cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02); cb.set_label(r"$\log_{10}$(uncertainty [%])")
ax2 = fig.add_subplot(gs[0, 1])
Nup = [0, 0, 0, 0, 7.6, 3.2, 3.2, 0, 4.5, 8.2]; Ndn = [0, 0, 0, 0, 7.6, 3.2, 3.2, 0, 4.5, 8.6]
for jj, j in enumerate(range(9, -1, -1)):
    st, sy, mc = rows["stat"][j] / 100, tot[j] / 100, rows["MC"][j] / 100
    nu, nd = Nup[j] / 100, Ndn[j] / 100
    lev = [("stat", st, st, "#08306b"), ("stat (+) syst", math.hypot(st, sy), math.hypot(st, sy), "#2171b5"),
           ("(+) MC", math.sqrt(st ** 2 + sy ** 2 + mc ** 2), math.sqrt(st ** 2 + sy ** 2 + mc ** 2), "#6baed6"),
           ("(+) N", math.sqrt(st ** 2 + sy ** 2 + mc ** 2 + nu ** 2), math.sqrt(st ** 2 + sy ** 2 + mc ** 2 + nd ** 2), "#c6dbef")]
    for name, up, dn, col in lev[::-1]:
        ax2.fill_betweenx([jj - 0.35, jj + 0.35], 1 - dn, 1 + up, color=col, label=name if jj == 0 else None, lw=0)
    ax2.plot([1, 1], [jj - 0.4, jj + 0.4], color="k", lw=1.5)
ax2.set_yticks(range(10)); ax2.set_yticklabels(TEXH[::-1])
ax2.set_xlim(-1.2, 3.2); ax2.set_xlabel("value / nominal value"); ax2.axvline(0, color="grey", lw=0.8, ls=":")
ax2.set_title("(b) nested uncertainty ribbons")
h, l = ax2.get_legend_handles_labels(); ax2.legend(h[::-1], l[::-1], loc="lower right", fontsize=12, frameon=True)
fig.text(0.01, 0.005, r"S0, isBest $\geq$ 0, L = 6079.2 nb$^{-1}$; PDG 2026. MAIN: $Y_\mathrm{peak}$, $\sigma_\eta$ (side bands); NOT DETERMINED: $Y_\mathrm{bump}$ (B4), DR; "
         r"CROSS-CHECK (closure failed): $Y_5$, $Y_4$, $a_\eta$, BR$_A$, $\kappa$, BR$_G$. S12: WP8 + composition of the peak from the data ($\sigma_\eta$); S8 widened.", fontsize=11)
for ext in ("pdf", "png"):
    fig.savefig(FIG + "fix_S1." + ext, dpi=80 if ext == "png" else None, bbox_inches="tight")
