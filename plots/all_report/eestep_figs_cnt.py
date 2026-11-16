#!/usr/bin/env python3
# ============================================================================
# eestep_figs.py  --  package eestep: figures (from eestep_ana.json and eestep_figdata.npz)
# ============================================================================
# Written 2026-11-15 (prefix eestep). Output (pdf + png) in /hdd3/przygoda/HADES/FAT_out/all_report/eestep/:
#   eestep_spectra  u spectra of the like-sign rows per M(ee) band with the edgeu step fit
#   eestep_bands    step heights per M(ee) band: like-sign rows; N, CB, N - CB against SMASH x f
#   eestep_align    is the e-e- step tied to E(M(pipi), M(ee))? delta chi2 per alignment variable; free edge position
# Run: python3 -I plots/all_report/eestep_figs.py
# ============================================================================
# _cnt copy (2026-11-15, FIGURE AGENT D of the integration run): the analyst calls the stored pi+pi-e+e- combinations
# COUNTS, not rows. Only label, legend, title and annotation strings 'rows' -> 'counts' and the output names (suffix _cnt)
# differ from the original eestep_figs.py; the inputs and the plotted content are the same.
# Written to docs/latex/figs_all/<name>_cnt.{pdf,png}; eestep_align (no 'rows' in it) is not written.
import json, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
OUT = "/hdd3/przygoda/HADES/FAT_out/all_report/eestep/"
R = json.load(open(OUT + "eestep_ana.json"))
F = np.load(OUT + "eestep_figdata.npz")
plt.rcParams.update({"font.size": 11, "axes.labelsize": 12, "axes.titlesize": 11, "legend.fontsize": 9})
x = 1e3 * F["x"]

FIG = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_all/"

def save(fig, name):
    if name == "eestep_align":
        plt.close(fig); return
    fig.savefig(FIG + name + "_cnt.pdf"); fig.savefig(FIG + name + "_cnt.png", dpi=110); plt.close(fig)

# ------------------------------------------------------------------ spectra
fig, ax = plt.subplots(2, 2, figsize=(10, 7.2), sharex=True)
panels = [("mm_lo", "(a) e$^-$e$^-$, M(ee) < 50 MeV/c$^2$", "tab:red"), ("mm_mid", "(b) e$^-$e$^-$, 50-90 MeV/c$^2$", "tab:red"),
          ("mm_hi", "(c) e$^-$e$^-$, M(ee) $\\geq$ 90 MeV/c$^2$", "tab:red"), ("pp_hi", "(d) e$^+$e$^+$, M(ee) $\\geq$ 90 MeV/c$^2$", "tab:purple")]
band_of = {"mm_lo": (0.0, 0.05, "hmm"), "mm_mid": (0.05, 0.09, "hmm"), "mm_hi": (0.09, 0.1349768, "hmm"), "pp_hi": (0.09, 0.1349768, "hpp")}
for a, (key, title, col) in zip(ax.ravel(), panels):
    y = F["y_" + key]; idx = F["i_" + key]
    sel = (x > -60) & (x < 60)
    a.errorbar(x[sel], y[sel], yerr=np.sqrt(np.maximum(y[sel], 1)), fmt="o", ms=3.5, color=col, label="data counts")
    a.plot(x[idx], F["m_" + key], "-", color="k", lw=1.6, label="fit: smooth part + step")
    a.plot(x[idx], F["s_" + key], "--", color="0.45", lw=1.2, label="smooth part of the fit")
    a.axvline(0, color="0.6", lw=0.8, ls=":")
    lo, hi, hk = band_of[key]
    b = [z for z in R["bands"] if abs(z["lo"] - lo) < 1e-6 and abs(z["hi"] - hi) < 1e-6][0]
    eh = b["eh" + hk[1:]]
    a.set_title(title + ":  step h = %.1f $\\pm$ %.1f" % (b[hk], eh))
    a.set_ylabel("counts / 2.5 MeV/c$^2$")
    a.set_ylim(0, 1.35 * y[sel].max())
    a.legend(loc="upper left", frameon=False)
for a in ax[1]:
    a.set_xlabel("u = M$_4$ - E(M($\\pi\\pi$), M(ee))  [MeV/c$^2$]")
fig.suptitle("Like-sign counts (S0, 0.30 $\\leq$ M($\\pi^+\\pi^-$) < 0.40 GeV/c$^2$) in the distance to the bump edge", fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.96))
save(fig, "eestep_spectra")

# ------------------------------------------------------------------ bands
fine = [b for b in R["bands"] if (b["lo"], round(b["hi"], 4)) in [(0.0, 0.01), (0.01, 0.03), (0.03, 0.05), (0.05, 0.07), (0.07, 0.09), (0.09, 0.11), (0.11, 0.135)]]
xc = np.array([0.5e3 * (b["lo"] + b["hi"]) for b in fine]); xw = np.array([0.5e3 * (b["hi"] - b["lo"]) for b in fine])
fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
a = ax[0]
a.errorbar(xc - 1, [b["hmm"] for b in fine], xerr=xw, yerr=[b["ehmm"] for b in fine], fmt="o", color="tab:red", label="e$^-$e$^-$ counts")
a.errorbar(xc + 1, [b["hpp"] for b in fine], xerr=xw, yerr=[b["ehpp"] for b in fine], fmt="s", mfc="none", color="tab:purple", label="e$^+$e$^+$ counts")
a.axhline(0, color="0.5", lw=0.8)
a.set_xlabel("M(ee) [MeV/c$^2$]"); a.set_ylabel("step height h [counts / 2.5 MeV/c$^2$]")
a.set_title("(a) like-sign counts: the step sits at large M(ee)")
a.legend(frameon=False)
a = ax[1]
hi = [i for i, b in enumerate(fine) if b["lo"] >= 0.0499]
xh, wh = xc[hi], xw[hi]
fh = [fine[i] for i in hi]
a.errorbar(xh - 1.5, [b["hN"] for b in fh], xerr=wh, yerr=[b["ehN"] for b in fh], fmt="o", color="k", label="N (e$^+$e$^-$ counts)")
a.errorbar(xh + 1.5, [b["hd"] for b in fh], xerr=wh, yerr=[b["ehd"] for b in fh], fmt="o", color="tab:blue", label="N - CB")
a.step(np.r_[xh - wh, xh[-1] + wh[-1]], np.r_[[b["hcb"] for b in fh], fh[-1]["hcb"]], where="post", color="tab:red", lw=1.6, label="CB = 2$\\sqrt{N_{++}N_{--}}$")
a.step(np.r_[xh - wh, xh[-1] + wh[-1]], np.r_[[0.42 * b["hs"] for b in fh], 0.42 * fh[-1]["hs"]], where="post", color="tab:green", lw=1.6,
       label="SMASH bump x f (f = 0.420)")
a.axhline(0, color="0.5", lw=0.8)
a.set_xlim(48, 137)
a.set_xlabel("M(ee) [MeV/c$^2$]"); a.set_ylabel("step height h [counts / 2.5 MeV/c$^2$]")
r90 = R["readings"]["90"]
a.set_title("(b) e$^+$e$^-$ counts, M(ee) $\\geq$ 50 MeV/c$^2$: N - CB follows SMASH x f")
a.text(0.03, 0.04, "M(ee) $\\geq$ 90 MeV/c$^2$: N - CB = %.1f $\\pm$ %.1f, SMASH x f = %.1f;\nif the like-sign step were an artefact: %.1f (data %.1f$\\sigma$ away)" %
       (r90["hd"], r90["ehd"], r90["fhs"], r90["fhs"] - r90["hcb"], r90["z_art"]), transform=a.transAxes, fontsize=9)
a.legend(frameon=False, loc="upper right")
fig.tight_layout()
save(fig, "eestep_bands")

# ------------------------------------------------------------------ alignment
fig, ax = plt.subplots(1, 2, figsize=(11, 4.0))
a = ax[0]
names = list(R["align"].keys())
lab = {"u = M4 - E(m, q) (nominal)": "u = M$_4$ - E(m, q)", "M4 - E(m, 0) (no M(ee) term)": "M$_4$ - E(m, 0)", "M4 - 0.535 (fixed M4)": "M$_4$ - 0.535",
       "M4 - E(m, q) - 10 MeV": "u - 10 MeV", "M4 - E(m, q) + 10 MeV": "u + 10 MeV"}
vals = [R["align"][n]["dchi2"] for n in names]
a.barh(range(len(names)), vals, color=["tab:blue" if i == 0 else "0.6" for i in range(len(names))])
a.set_yticks(range(len(names))); a.set_yticklabels([lab.get(n, n) for n in names])
a.invert_yaxis()
a.set_xlabel("$\\Delta\\chi^2$ (no step - step)")
a.set_title("(a) e$^-$e$^-$, M(ee) $\\geq$ 70 MeV/c$^2$: which variable has the step?")
a = ax[1]
pr = np.array(R["u0_profile"])
a.plot(1e3 * pr[:, 0], pr[:, 1] - pr[:, 1].min(), "o-", color="tab:red")
a.axhline(1, color="0.6", ls=":"); a.axvline(0, color="0.6", ls=":")
a.set_xlabel("position of the step u$_0$ [MeV/c$^2$]"); a.set_ylabel("$\\chi^2$ - $\\chi^2_{min}$")
a.set_title("(b) free step position: the step is at the bump edge")
fig.tight_layout()
save(fig, "eestep_align")
print("figures written to", FIG, "(eestep_spectra_cnt, eestep_bands_cnt)")
