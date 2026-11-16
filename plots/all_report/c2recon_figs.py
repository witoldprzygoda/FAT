#!/usr/bin/env python3
# ============================================================================
# c2recon_figs.py  --  package c2recon: one summary figure (written to the c2recon output directory, not to docs/)
# ============================================================================
# (a) the C2 gamma* factor eps_QED/eps_iso from every source: the earlier numbers (toys, closure, applied values,
#     ecalacc window means) and the factor each estimator actually needs (this package);
# (b) the factor against the momentum of the softer lepton: full simulation (C2 rows, circle-averaged weights) and the
#     toy (p_e > 0.06 GeV/c, true weights), with the share of rows in each bin (full simulation vs toy).
# Output: /hdd3/przygoda/HADES/FAT_out/all_report/c2recon/c2recon_factors.{pdf,png}
# Run: python3 -I plots/all_report/c2recon_figs.py
# ============================================================================
import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import c2recon_weights as CW

O = "/hdd3/przygoda/HADES/FAT_out/all_report/c2recon/"
N = json.load(open(O + "c2recon_numbers.json"))
EU = N["edgeu"]; BE = N["bedge"]; TY = N["toy"]
plt.rcParams.update({"font.size": 12})
fig, ax = plt.subplots(1, 2, figsize=(15, 6.2), gridspec_kw=dict(width_ratios=[1.25, 1]))

# ---- (a)
D = EU["definitions"]
items = [
    ("WP8 toy, all M4 (p_e > 60 MeV/c)", 0.880, None, "grey"),
    ("WP8 closure lep/std (prompt pi0)", 0.891, None, "grey"),
    ("bedge edge toy (p_e > 60 MeV/c)", 0.8826, None, "grey"),
    ("verifier u-edge toy (p_e > 60 MeV/c)", 0.879, None, "grey"),
    ("applied in bedge (0.879/0.882), effective", BE["c_applied"]["all"], None, "orange"),
    ("applied in edgeu (0.879/0.8825), effective", (0.879 * EU["factor"]["E1"]["q_A"]["hiso"] + 0.8825 * EU["factor"]["E2"]["q_A"]["hiso"]) / EU["factor"]["all"]["q_A"]["hiso"], None, "orange"),
    ("full sim.: all C2 rows (all M4)", D["all C2 rows (S0, all M(pipi), all M4) | all"]["c"], 0.0027, "tab:blue"),
    ("full sim.: window B", D["window B (0.4550 <= M4 < 0.5225) | all"]["c"], 0.0031, "tab:blue"),
    ("full sim.: -15...+5 MeV/c2 from the edge", D["edge window -15...+5 MeV in u (ecalacc) | all"]["c"], 0.0053, "tab:blue"),
    ("full sim.: far below the edge (u < -60)", D["far below u < -60 MeV | all"]["c"], 0.0058, "tab:blue"),
    ("full sim.: ECAL-tagged rows", N["ecal"]["c_tag"], N["ecal"]["ec_tag"], "tab:green"),
    ("NEEDED by bedge (its own fit)", BE["c"]["all"], BE["stat"]["pois2000"]["all"], "red"),
    ("NEEDED by u-edge (its own fit)", EU["factor"]["all"]["q_A"]["c"], EU["stat"]["all"]["rms"], "red"),
]
y = np.arange(len(items))[::-1]
for yi, (lab, v, e, col) in zip(y, items):
    ax[0].errorbar([v], [yi], xerr=[e] if e else None, fmt="o", color=col, ms=8, capsize=4)
ax[0].set_yticks(y); ax[0].set_yticklabels([it[0] for it in items], fontsize=11)
ax[0].axvspan(BE["c_applied"]["all"] - 0.024, BE["c_applied"]["all"] + 0.024, color="orange", alpha=0.12, label="D35 allowance +-0.024")
ax[0].set_xlabel(r"C2 factor $\varepsilon_{\rm QED}/\varepsilon_{\rm iso}$")
ax[0].set_xlim(0.83, 0.91)
ax[0].legend(loc="upper left", fontsize=9)
ax[0].set_title("(a) the factor from every source", fontsize=13)
ax[0].grid(axis="x", alpha=0.3)

# ---- (b)
R = CW.load_rows(); Wt = np.load(O + "c2recon_weights.npz")
W = R["W108"]; q = Wt["q_A"]; s0 = R["s0"] > 0
ps = np.minimum(np.sqrt(R["ep_px"] ** 2 + R["ep_py"] ** 2 + R["ep_pz"] ** 2), np.sqrt(R["em_px"] ** 2 + R["em_py"] ** 2 + R["em_pz"] ** 2)) * 1e-3
pb = np.array([0.04, 0.06, 0.08, 0.10, 0.15, 0.20, 0.30, 0.50])
gk = R["tree"].astype(np.int64) * (1 << 32) + R["grp"].astype(np.int64)
uq, inv = np.unique(gk, return_inverse=True)
K = np.random.default_rng(7).poisson(1.0, size=(1000, uq.size))[:, inv]
xs, vs, es, sh = [], [], [], []
for a, b in zip(pb[:-1], pb[1:]):
    m = s0 & (ps >= a) & (ps < b)
    if m.sum() < 10:
        continue
    xs.append(0.5 * (a + b)); vs.append((W * q)[m].sum() / W[m].sum())
    es.append(((K * (W * q * m)).sum(1) / np.maximum((K * (W * m)).sum(1), 1e-9)).std()); sh.append(W[m].sum() / W[s0].sum())
ax[1].errorbar(np.array(xs) * 1e3, vs, yerr=es, fmt="o", color="red", ms=8, capsize=4, label="full simulation, C2 rows (S0)")
t = TY["0.06"]["vs_pmin"]
ax[1].plot([1e3 * 0.5 * (r[0] + r[1]) for r in t], [r[3] for r in t], "s--", color="grey", ms=7, label="toy, p_e > 60 MeV/c (true weight)")
ax[1].set_xscale("log")
ax[1].set_xlabel("momentum of the softer lepton [MeV/c]")
ax[1].set_ylabel(r"$\varepsilon_{\rm QED}/\varepsilon_{\rm iso}$")
ax[1].set_ylim(0.74, 0.96)
ax2 = ax[1].twinx()
tn = np.array([r[2] for r in t], float)
ax2.bar(np.array(xs) * 1e3, sh, width=np.diff(pb)[:len(xs)] * 1e3 * 0.8, alpha=0.18, color="red", label="share of rows, full simulation")
ax2.bar([1e3 * 0.5 * (r[0] + r[1]) for r in t], tn / tn.sum(), width=[1e3 * (r[1] - r[0]) * 0.4 for r in t], alpha=0.25, color="grey",
        label="share of rows, toy")
ax2.set_ylabel("share of rows")
ax2.set_ylim(0, 1.0)
h1, l1 = ax[1].get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
ax[1].legend(h1 + h2, l1 + l2, loc="upper right", fontsize=10)
ax[1].set_title("(b) the factor is set by the softer lepton", fontsize=13)
ax[1].grid(alpha=0.3)
fig.tight_layout()
fig.savefig(O + "c2recon_factors.pdf"); fig.savefig(O + "c2recon_factors.png", dpi=80)
print("written", O + "c2recon_factors.pdf")
