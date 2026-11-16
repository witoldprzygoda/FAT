#!/usr/bin/env python3
# ============================================================================
# eestep_partII_fig.py  --  package eestep: data-only version of the band figure for Part II
# ============================================================================
# Written 2026-11-15 (Amend step of the report integration). Reads only eestep_ana.json (no new analysis).
# Output (pdf + png) in /hdd3/przygoda/HADES/FAT_out/all_report/eestep/:
#   eestep_bands_data  (a) like-sign step heights per M(ee) band; (b) e+e- rows N, CB and N - CB per band (no simulation)
# Run: python3 -I plots/all_report/eestep_partII_fig.py
# ============================================================================
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
OUT = "/hdd3/przygoda/HADES/FAT_out/all_report/eestep/"
R = json.load(open(OUT + "eestep_ana.json"))
plt.rcParams.update({"font.size": 11, "axes.labelsize": 12, "axes.titlesize": 11, "legend.fontsize": 9})

fine = [b for b in R["bands"] if (b["lo"], round(b["hi"], 4)) in [(0.0, 0.01), (0.01, 0.03), (0.03, 0.05), (0.05, 0.07), (0.07, 0.09),
                                                                    (0.09, 0.11), (0.11, 0.135)]]
xc = np.array([0.5e3 * (b["lo"] + b["hi"]) for b in fine]); xw = np.array([0.5e3 * (b["hi"] - b["lo"]) for b in fine])
fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
a = ax[0]
a.errorbar(xc - 1, [b["hmm"] for b in fine], xerr=xw, yerr=[b["ehmm"] for b in fine], fmt="o", color="tab:red", label="e$^-$e$^-$ rows")
a.errorbar(xc + 1, [b["hpp"] for b in fine], xerr=xw, yerr=[b["ehpp"] for b in fine], fmt="s", mfc="none", color="tab:purple",
           label="e$^+$e$^+$ rows")
a.axhline(0, color="0.5", lw=0.8)
a.set_xlabel("M(ee) [MeV/c$^2$]"); a.set_ylabel("step height h [rows / 2.5 MeV/c$^2$]")
a.set_title("(a) like-sign rows: the step sits at large M(ee)")
a.legend(frameon=False)
a = ax[1]
hi = [i for i, b in enumerate(fine) if b["lo"] >= 0.0499]
xh, wh = xc[hi], xw[hi]
fh = [fine[i] for i in hi]
a.errorbar(xh - 1.5, [b["hN"] for b in fh], xerr=wh, yerr=[b["ehN"] for b in fh], fmt="o", color="k", label="N (e$^+$e$^-$ rows)")
a.errorbar(xh + 1.5, [b["hd"] for b in fh], xerr=wh, yerr=[b["ehd"] for b in fh], fmt="o", color="tab:blue", label="N - CB")
a.step(np.r_[xh - wh, xh[-1] + wh[-1]], np.r_[[b["hcb"] for b in fh], fh[-1]["hcb"]], where="post", color="tab:red", lw=1.6,
       label="CB = 2$\\sqrt{N_{++}N_{--}}$")
a.axhline(0, color="0.5", lw=0.8)
a.set_xlim(48, 137)
a.set_xlabel("M(ee) [MeV/c$^2$]"); a.set_ylabel("step height h [rows / 2.5 MeV/c$^2$]")
r90 = R["readings"]["90"]
a.set_title("(b) e$^+$e$^-$ rows, M(ee) $\\geq$ 50 MeV/c$^2$: N carries the CB step")
a.text(0.03, 0.04, "M(ee) $\\geq$ 90 MeV/c$^2$: N = %.1f $\\pm$ %.1f, CB = %.1f $\\pm$ %.1f,\nN - CB = %.1f $\\pm$ %.1f" %
       (r90["hN"], r90["ehN"], r90["hcb"], r90["ehCB"], r90["hd"], r90["ehd"]), transform=a.transAxes, fontsize=9)
a.legend(frameon=False, loc="upper right")
fig.tight_layout()
fig.savefig(OUT + "eestep_bands_data.pdf"); fig.savefig(OUT + "eestep_bands_data.png", dpi=110); plt.close(fig)
print("written", OUT + "eestep_bands_data.{pdf,png}")
