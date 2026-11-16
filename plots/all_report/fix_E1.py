#!/usr/bin/env python3
# ============================================================================
# fix_E1.py  --  review fix (prefix fix): Fig. E1 (F-main, CROSS-CHECK) redrawn in the common palette of report_all, with the
# legend outside the data. The fit is redone with the independent library of the verifier v5 (plots/all_report/v5_fit.py,
# read only), which reproduces WP5 to 1e-5 (a_eta 0.382284 vs 0.382283).
# Top row: F10 E2 (the fitted class), pulls (N - mu)/sigma with sigma^2 = event-group variance of N scaled to the model + k^2 var(CB)
#          + template variance. Bottom row: E1 (NOT fitted, fall-back "drop E1"): prediction with the E2 parameters and k_E1 = 1,
#          and data/model.
# Output: docs/latex/figs_all/fix_E1.{pdf,png}
# Run: python3 -I plots/all_report/fix_E1.py   (from /home/przygoda/HADES/FAT/FAT; under a minute)
# ============================================================================
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import v5_fit as F

FIG = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_all/"
plt.rcParams.update({"font.size": 12, "axes.labelsize": 13, "axes.titlesize": 12, "legend.fontsize": 12, "xtick.labelsize": 11, "ytick.labelsize": 11, "pdf.fonttype": 42})
PAL = dict(CB="#e8202a", bg="#b0b0b0", Geta="#c9a3f5", C34="#d400d4", C2="#ff6a00", C5="#00cccc")

T, _ = F.load_sim("S0"); X = F.load_data("S0")
S, _ = F.smooth_all(T, lambda c, h: T[(c, h)])
hE2 = ["F10_S%dE2" % j for j in range(1, 6)]
B = F.make_bins(X, T, S, hE2)
th0 = np.array([0.4, 1.0, 1.0, 1.0, 1, 1, 1, 1, 1, -0.07, 2.0, 1.0])
th, _ = F.fit_irls(B, F.model_main, F.CONS_MAIN, th0)
a, ka, rho, g = th[0:4]; c = th[4:9]; b = th[9]; kB = th[10]; kE2 = th[11]
print("F-main (v5 library): a_eta %.6f kappa %.4f rho %.4f g %.4f b %.4f kB %.3f kE2 %.4f" % (a, ka, rho, g, b, kB, kE2))
edges = F.edges; cen = 0.5 * (edges[1:] + edges[:-1]); w = (edges[1:] - edges[:-1]) * 1e3   # per MeV/c2
x = (cen - 0.55) / 0.1
fig = plt.figure(figsize=(17, 9.5))
gs = fig.add_gridspec(5, 5, height_ratios=[3, 1.1, 0.55, 3, 1.1], hspace=0.06, wspace=0.22, left=0.06, right=0.99, top=0.95, bottom=0.17)
for e, row in ((2, 0), (1, 3)):
    for j in range(1, 6):
        h = "F10_S%dE%d" % (j, e)
        N = X[("N", h)]; vN = X[("Nerr", h)] ** 2; CB = X[("CB", h)]
        vCB = np.var(X[("CBrep", h)], axis=0, ddof=1)
        k = kE2 if e == 2 else 1.0
        bg = c[j - 1] * (1 + b * x) * (S[("GD", h)] + kB * S[("Gconv", h)])
        comp = [("CB", k * CB), ("bg", bg), ("Geta", a * g * T[("Geta", h)]), ("C34", a * ka * (T[("C3", h)] + T[("C4", h)])),
                ("C2", a * T[("C2", h)]), ("C5", a * rho * T[("C5", h)])]
        mu = sum(v for _, v in comp)
        vs = (a * g) ** 2 * T[("Geta", h, "v")] + (a * ka) ** 2 * (T[("C3", h, "v")] + T[("C4", h, "v")]) + a ** 2 * T[("C2", h, "v")] + (a * rho) ** 2 * T[("C5", h, "v")]
        sN = np.where(N > 0, vN / np.maximum(N, 1) * np.maximum(mu, 1), np.maximum(mu, 1))
        sig = np.sqrt(sN + k * k * vCB + vs)
        ax = fig.add_subplot(gs[row, j - 1]); ap = fig.add_subplot(gs[row + 1, j - 1], sharex=ax)
        base = np.zeros(70)
        labels = dict(CB="CB, like-sign data (red)", bg="non-$\\eta$ background $c_j(1+bx)(G_D+\\kappa_BG_C)$, smoothed",
                      Geta="$G_\\eta$ ($a_\\eta g$)", C34="C3 + C4, conversions ($a_\\eta\\kappa$)", C2="C2 bump $\\pi^0$ Dalitz ($a_\\eta$)",
                      C5="C5 direct ($a_\\eta\\rho_d$)")
        for key, v in comp:
            top = base + v
            ax.fill_between(edges, np.append(base / w, base[-1] / w[-1]), np.append(top / w, top[-1] / w[-1]), step="post", color=PAL[key], lw=0,
                            label=labels[key] if (row == 0 and j == 1) else None)
            base = top
        ax.errorbar(cen, N / w, yerr=np.sqrt(vN) / w, fmt="o", ms=2.8, color="black", label="data $N$" if (row == 0 and j == 1) else None)
        lo = 0.0; hi = 1.12 * max(np.max((N + np.sqrt(vN)) / w), np.max(mu / w))
        ax.set_ylim(lo, hi); ax.set_xlim(0.40, 0.70)
        plt.setp(ax.get_xticklabels(), visible=False)
        ax.set_title("S%d, $M(e^+e^-)%s10$ MeV/$c^2$%s" % (j, "\\geq" if e == 2 else "<", "" if e == 2 else " (not fitted)"), fontsize=12)
        if j == 1: ax.set_ylabel("rows / (MeV/$c^2$)")
        if e == 2:
            pull = (N - mu) / sig
            ap.bar(cen, pull, width=w * 1e-3, color="#444444")
            ap.axhline(0, color="k", lw=0.8); ap.set_ylim(-4.5, 4.5)
            if j == 1: ap.set_ylabel("pull")
        else:
            r = N / np.where(mu > 0, mu, np.nan); er = np.sqrt(vN) / np.where(mu > 0, mu, np.nan)
            ap.errorbar(cen, r, yerr=er, fmt="o", ms=2.5, color="black")
            ap.axhline(1, color="k", lw=0.8, ls=":"); ap.set_ylim(0.5, 1.5)
            if j == 1: ap.set_ylabel("data/model")
        ap.set_xlabel("$M_4$ [GeV/$c^2$]" if row == 3 else "")
        if row == 0: plt.setp(ap.get_xticklabels(), visible=False)
handles, labs = fig.axes[0].get_legend_handles_labels()
fig.legend(handles, labs, loc="lower center", ncol=4, frameon=False, fontsize=12)
fig.text(0.5, 0.005, "F-main (CROSS-CHECK, failed the closure test): S0, isBest $\\geq0$, $L=6079.2$ nb$^{-1}$; $a_\\eta=%.3f$, $\\kappa=%.2f$, "
         "$\\rho_d=%.2f$, $g=%.2f$, $k_{E2}=%.3f$; lower half: E1 prediction with the E2 parameters and $k_{E1}=1$" % (a, ka, rho, g, kE2),
         ha="center", fontsize=12)
fig.savefig(FIG + "fix_E1.pdf"); fig.savefig(FIG + "fix_E1.png", dpi=60)
print("done")
