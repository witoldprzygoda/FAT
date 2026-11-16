#!/usr/bin/env python3
# vflex_fig.py -- independent verification of package flex: one summary figure
#   left : the nominal F8 scales as a function of the weight lambda of the BB-lite template variance in the fit weights
#          (lambda = 0: Poisson-like weights, templates exact; 1: flex nominal); bands = flex total error at lambda = 1
#   right: hybrid GLS chi2 of 2 x 1000 independent RF pseudo-experiments (truth = F8 fit / SMASH as is) and the data value
# Run: python3 -B plots/all_report/vflex_fig.py
import sys, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import vflex_core as V

D = V.load(); L = V.Layout(D, (1,))
M = V.FlexModel(["C2", "C3", "C4", "C5", "Geta", "GD_om", "GD_nom", "GC"]); T, Vt = V.templates(D, L, M)
lams = [0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0]
TH = []
for lam in lams:
    R = V.fit(M, L, T, lam * Vt, egcov=L.eg, sandwich=False)
    TH.append(R['th'])
TH = np.array(TH)
fj = json.load(open(V.DIR + "/vflex_fit.json"))["F8"]
json.dump(dict(lams=lams, th=TH.tolist(), names=M.names), open(V.DIR + "/vflex_lamscan.json", "w"), indent=1)
for lam, th in zip(lams, TH):
    print(f"lambda {lam:4.2f}: " + " ".join(f"{n}={x:.3f}" for n, x in zip(M.names, th)))
# Part III palette (report_all common palette: C2 orange, C3+C4 magenta, C5 cyan, G_eta violet; non-eta greys)
sel = [("s_C2", "C2 (bump, $\\pi^0$ Dalitz)", "#ff6600"), ("s_C5", "C5 (direct $\\eta\\to\\pi^+\\pi^-e^+e^-$)", "#00b3b3"),
       ("s_Geta", "$G_\\eta$ (other $\\eta$)", "#8a5cd6"), ("s_GD_nom", "$G_D$, other pions", "#8c8c8c"),
       ("s_GD_om", "$G_D$, $\\omega/\\rho^0$ pions", "#4d4d4d"), ("s_GC", "$G_C$ (conversion pair)", "#000000")]
plt.rcParams.update({"font.size": 13})
fig, ax = plt.subplots(1, 2, figsize=(13, 5.2), gridspec_kw=dict(width_ratios=[1.25, 1]))
for n, lab, col in sel:
    k = M.names.index(n)
    ls = "--" if n == "s_GD_om" else "-"
    ax[0].plot(lams, TH[:, k], ls, color=col, marker="o", ms=4, lw=2, label=lab)
    ax[0].errorbar([1.0 + 0.03 * (sel.index((n, lab, col)) - 2.5)], [fj["th"][k]], yerr=[fj["tot"][k]], color=col, capsize=3, lw=1.5)
ax[0].axhline(1, color="k", ls=":", lw=1)
ax[0].axvline(1, color="0.6", lw=1)
ax[0].set_xlabel("weight $\\lambda$ of the template variance in the fit weights (1 = flex nominal)")
ax[0].set_ylabel("scale relative to SMASH")
ax[0].set_title("F8 scales depend on the fit weights (same data, same templates)", fontsize=12)
ax[0].legend(fontsize=10.5, loc="upper left", ncol=2, frameon=False)
ax[0].set_ylim(-0.6, 4.6)
g = []
for s in (0, 1):
    z = np.load(V.DIR + f"/vflex_toys_s{s}.npz")
    g.append(z["gl"])
bins = np.linspace(250, 650, 81)
ax[1].hist(g[0], bins=bins, histtype="step", lw=2, color="#1f4e9e", label="truth = F8 fit\n(1000 toys)")
ax[1].hist(g[1], bins=bins, histtype="step", lw=2, ls="--", color="#7f7f7f", label="truth = SMASH as is\n(1000 toys)")
ax[1].axvline(fj["gls"], color="#d62728", lw=2.5)
ax[1].text(fj["gls"] - 8, 35, f"data: {fj['gls']:.1f}", color="#d62728", rotation=90, ha="right", va="center", fontsize=12)
ax[1].set_xlabel("hybrid GLS $\\chi^2$ (342 degrees of freedom)")
ax[1].set_ylabel("pseudo-experiments per 5 units")
ax[1].set_title("goodness of fit of F8: no pseudo-experiment reaches the data", fontsize=12)
ax[1].legend(fontsize=11, frameon=False, loc="center", bbox_to_anchor=(0.6, 0.72))
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(f"/home/przygoda/HADES/FAT/FAT/docs/latex/figs_all/vflex_summary.{ext}", dpi=150)
print("written docs/latex/figs_all/vflex_summary.{pdf,png}")
