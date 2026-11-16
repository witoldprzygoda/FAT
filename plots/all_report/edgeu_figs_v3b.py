# edgeu_figs_v3b.py -- COPY of edgeu_figs_v3.py by the final checker (integ4, 2026-11-15): writes only edgeu_checks_v3b
# (same content as edgeu_checks_v3); the x label of panel (h) printed "5 %%", now "5 %".
#!/usr/bin/env python3
# ============================================================================
# edgeu_figs.py  --  package edgeu: figures (Part III colours: data black/blue, CB red, C2 orange, C3 magenta, C4 green, C5 cyan,
#                    other eta violet, non-eta greys; fit green as bedge_steps)
# ============================================================================
# Written 2026-11-15 (prefix edgeu). Inputs: /hdd3/.../edgeu/edgeu_figdata.npz, edgeu_ana.json (edgeu_ana.py), edgeu_toys.json
# (edgeu_toys.py), edgeu_numbers.json (edgeu_combine.py). Outputs docs/latex/figs_all/edgeu_{kinematics,spectra,fit,checks,syst,
# toys,compare}.{pdf,png}.
# Run: python3 -I plots/all_report/edgeu_figs.py
# ----------------------------------------------------------------------------
# _v3 copy (2026-11-15, FIGURE AGENT C of the fourth integration): rulings D38 (C2 gamma* angle corrected with the per-count
# QED weights of c2recon instead of the class factors 0.879 / 0.8825) and D40 (peak-tail family t = 1 +- 0.5 counted in U8),
# and the analyst's terminology (COUNTS, not rows). Differences to edgeu_figs.py:
#   inputs  edgeu_v3/edgeu_figdata_v3.npz, edgeu_v3/edgeu_ana_v3.json, edgeu_v3/edgeu_numbers_v3.json (edgeu_ana_v3.py,
#           edgeu_combine_v3.py); edgeu_toys.json unchanged. The uncorrected u readings of edgeu_checks (a) come from the
#           original edgeu_ana.json: in the v3 run the C2 weights carry the QED weight, so lin['C2'] = lin['C2c'] there,
#           and D38 does not change the uncorrected reading (edgeu_ana_v3.py header);
#   labels  'rows' -> 'counts'; edgeu_spectra (b): C2 'x0.88' -> 'QED weight per count, about 0.86';
#   values  edgeu_checks (b): per-slice edge E1 / E2 (C2 corr.) 0.3842 / 0.4752 -> 0.3859 / 0.4870 (bedge_v3/ana_results_v3.txt);
#           edgeu_syst: the D40 peak-tail term (+-tail_rel of edgeu_numbers_v3.json) drawn as two U8 points (t = 0.5, 1.5);
#           edgeu_compare: per-slice edge 0.437 +- 0.023 +- 0.036 +- 0.049 -> 0.442 +- 0.024 +- 0.040 +- 0.050 and the peak
#           syst 0.152 -> 0.157 mb (D40; tables_all/sigma_headline_v3.tex).
# Written: docs/latex/figs_all/edgeu_{fit,spectra,checks,syst,compare}_v3.{pdf,png}; edgeu_kinematics and edgeu_toys do not
# change (edgeu_kinematics_cnt.pdf from edgeu_figs_cnt.py; edgeu_toys has no counts label) and are not written.
# Run: python3 -I plots/all_report/edgeu_figs_v3.py
# ============================================================================
import json
import math
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import edgeu_core as C

plt.rcParams.update({"font.size": 13, "axes.labelsize": 14, "axes.titlesize": 13.5, "legend.fontsize": 11, "xtick.labelsize": 12,
                     "ytick.labelsize": 12, "pdf.fonttype": 42})
COL = dict(data="#1f3fbf", data0="black", cb="#e8202a", c2="#ff8c00", c3="#d400d4", c4="#009900", c5="#00b5c8", eta="#7b3fbf",
           gd="#7f7f7f", gc="#b5b5b5", fit="#009900", bedge="#8c564b", b4="#ff6a00", tag="#9467bd")
OUT3 = "/hdd3/przygoda/HADES/FAT_out/all_report/edgeu_v3/"
FD = dict(np.load(OUT3 + "edgeu_figdata_v3.npz"))
A = json.load(open(OUT3 + "edgeu_ana_v3.json"))
A0 = json.load(open(C.OUT + "edgeu_ana.json"))      # uncorrected readings (unchanged by D38)
TY = json.load(open(C.OUT + "edgeu_toys.json"))
NB = json.load(open(OUT3 + "edgeu_numbers_v3.json"))
x = FD["x"] * 1e3
f0 = A["nom"]["f"]
def save(fig, name):
    if name not in ("edgeu_checks",):
        plt.close(fig); return
    fig.savefig(C.FIG + name + "_v3b.pdf"); fig.savefig(C.FIG + name + "_v3b.png", dpi=130); plt.close(fig)

# ---------------------------------------------------------------- 1. kinematics
fig, axs = plt.subplots(2, 2, figsize=(14.5, 10.5))
ax = axs[0, 0]
m = np.linspace(0.279, 0.4116, 400)
for q, ls, lab in [(0.0, "-", "$M_{ee}$ = 0 (conversion C3)"), (0.020, "--", "20 MeV/$c^2$"), (0.050, "-.", "50 MeV/$c^2$"), (0.100, ":", "100 MeV/$c^2$")]:
    ax.plot(m, C.edge(m, np.full_like(m, q)), ls, color=COL["c2"] if q > 0 else COL["c3"], lw=2.2, label=lab)
ax.axhline(0.5439, color=COL["c4"], lw=1.5); ax.text(0.282, 0.5452, "$\\eta$ peak 0.5439 (data scale)", color=COL["c4"], fontsize=11)
for e in (0.30, 0.34, 0.37, 0.39, 0.40):
    ax.axvline(e, color="#999999", lw=0.8, ls=":")
ax.axvline(0.4116, color="k", lw=1.0); ax.text(0.4085, 0.475, "no bump\nabove 0.4116", fontsize=10.5, ha="right")
ax.set_xlabel("$M(\\pi^+\\pi^-)$ [GeV/$c^2$] (data scale)"); ax.set_ylabel("edge $E(M_{\\pi\\pi}, M_{ee})$ [GeV/$c^2$]")
ax.set_title("(a) kinematic edge of $\\eta\\to\\pi^+\\pi^-\\pi^0$, $\\pi^0\\to e^+e^-\\gamma$ ($\\gamma$ lost)")
ax.legend(loc="lower left", fontsize=10.5); ax.set_ylim(0.465, 0.56)
ax = axs[0, 1]
bins = np.arange(-160, 20.1, 2.5)
for gname, col in (("C2", COL["c2"]), ("C3", COL["c3"])):
    ut = FD["ut_" + gname] * 1e3; w = FD["utw_" + gname]
    ax.hist(ut, bins=bins, weights=w, histtype="step", color=col, lw=2, label="%s: $u_{\\mathrm{true}}$ with $M_{ee}$" % gname)
mt0 = C.EDGE_SCALE * (FD["m4t_C2"] - C.m4max_true(FD["mppt_C2"], np.zeros_like(FD["mppt_C2"])))
ax.hist(mt0 * 1e3, bins=bins, weights=FD["utw_C2"], histtype="step", color=COL["c2"], lw=1.3, ls="--", label="C2: $u_{\\mathrm{true}}$ with $E(m,0)$")
ax.axvline(0, color="k", lw=1)
ax.set_xlabel("$u_{\\mathrm{true}}=0.99274\\,[m_4^{\\mathrm{true}}-M_4^{\\max}(m_{\\pi\\pi}^{\\mathrm{true}},m_{ee}^{\\mathrm{true}})]$ [MeV/$c^2$]")
ax.set_ylabel("expected rows per 2.5 MeV/$c^2$"); ax.set_title("(b) simulated truth: no bump row above $u=0$")
ax.legend(loc="upper left", fontsize=10.5)
k = A["kin"]
ax.text(0.03, 0.55, "max excess: C2 %+.2f, C3 %+.2f MeV\nwith $E(m,0)$: %.1f %% of C2 above\nby $>1$ MeV" %
        (1e3 * k["truth_C2"]["max_excess"], 1e3 * k["truth_C3"]["max_excess"], 100 * k["truth_C2"]["frac_above_q0_1MeV"]), transform=ax.transAxes, fontsize=10.5)
ax = axs[1, 0]
bins = np.arange(-60, 40.1, 2.5)
sel = FD["umee_C2"] >= 0.020
ok = np.isfinite(FD["ur_C2"]) & sel; ok0 = np.isfinite(FD["ur0_C2"]) & sel
ax.hist(FD["ur_C2"][ok] * 1e3, bins=bins, weights=FD["utw_C2"][ok], histtype="step", color=COL["c2"], lw=2.2, label="$u$ with $E(m,M_{ee})$")
ax.hist(FD["ur0_C2"][ok0] * 1e3, bins=bins, weights=FD["utw_C2"][ok0], histtype="step", color="k", lw=1.5, ls="--", label="the same rows with $E(m,0)$")
ax.axvline(0, color="k", lw=1)
ax.set_xlabel("reconstructed $u$ [MeV/$c^2$]"); ax.set_ylabel("expected rows per 2.5 MeV/$c^2$")
ax.set_title("(c) simulated C2 with $M_{ee}\\geq 20$ MeV/$c^2$: the $M_{ee}$ term")
ax.legend(loc="upper right", fontsize=10.5)
ax.text(0.03, 0.68, "above +20 MeV/$c^2$:\n%.2f %% with $E(m,M_{ee})$\n%.2f %% with $E(m,0)$" % (100 * k["above_C2_E20_q"][1], 100 * k["above_C2_E20_q0"][1]), transform=ax.transAxes, fontsize=10.5)
ax = axs[1, 1]
d = FD["reso_d"] * 1e3; wd = FD["reso_w"]
bins = np.arange(-40, 40.1, 2.0)
h, e = np.histogram(d, bins=bins, weights=wd)
ax.step(0.5 * (e[1:] + e[:-1]), h, where="mid", color=COL["c2"], lw=2, label="C2+C3, $u_{\\mathrm{true}}>-15$ MeV/$c^2$")
xx = np.linspace(-40, 40, 400)
r = A["reso"]
for s_, ls, lab in [(r["sigma_core"] * 1e3, "-", "Gaussian core $\\sigma$ = %.2f MeV/$c^2$ (nominal)" % (r["sigma_core"] * 1e3)),
                    (r["sig_data"] * 1e3, "--", "data-like $\\sigma\\times0.893$ = %.2f MeV/$c^2$" % (r["sig_data"] * 1e3))]:
    ax.plot(xx, h.max() * np.exp(-0.5 * (xx - r["mean_core"] * 1e3) ** 2 / s_ ** 2), ls, color="k", lw=1.5, label=lab)
ax.set_xlabel("$u_{\\mathrm{reco}}-u_{\\mathrm{true}}$ [MeV/$c^2$]"); ax.set_ylabel("expected rows per 2 MeV/$c^2$")
ax.set_title("(d) $u$ resolution at the edge (simulation)")
ax.legend(loc="upper left", fontsize=10)
ax.text(0.62, 0.55, "rms %.1f MeV\nhalf 68 %%: %.1f MeV\nM4 part %.1f, edge part %.1f" % (r["rms"] * 1e3, r["half68"] * 1e3, r["rms_dm4"] * 1e3, r["rms_dE"] * 1e3), transform=ax.transAxes, fontsize=10.5)
fig.tight_layout(); save(fig, "edgeu_kinematics")

# ---------------------------------------------------------------- 2. spectra
fig, axs = plt.subplots(1, 2, figsize=(15.5, 6.4))
ax = axs[0]
m = (x > -150) & (x < 120)
ax.errorbar(x[m], FD["sig"][m], yerr=np.sqrt(FD["sigvar"][m]), fmt="o", ms=3, color=COL["data0"], label="data S0: $N-$CB")
ax.errorbar(x[m], FD["sig_pre"][m], fmt="s", ms=3, mfc="none", color=COL["data"], label="data pre: $N-$CB")
ax.plot(x[m], FD["cb"][m], "-", color=COL["cb"], lw=2, label="CB $=2\\sqrt{N_{++}N_{--}}$ (S0)")
ax.plot(x[m], FD["mm"][m], ":", color=COL["cb"], lw=1.6, label="$N_{--}$ (S0)")
ax.plot(x[m], FD["pp"][m], "-.", color=COL["cb"], lw=1.3, label="$N_{++}$ (S0)")
ax.set_yscale("log"); ax.set_ylim(5, 2e4)
ax.axvline(0, color="k", lw=1); ax.axvspan(-45, 45, color="#eeeeee", zorder=0)
ax.set_xlabel("$u=M_4-E(M_{\\pi\\pi},M_{ee})$ [MeV/$c^2$]"); ax.set_ylabel("counts per 2.5 MeV/$c^2$")
ax.set_title("(a) data, $0.30\\leq M(\\pi^+\\pi^-)<0.40$ GeV/$c^2$, $L$ = 6079.2 nb$^{-1}$")
ax.legend(loc="lower left", fontsize=10, ncol=2)
ax = axs[1]
for gname, col, ls, lab in [("C2", COL["c2"], "-", "C2 bump, $\\pi^0$ Dalitz (QED weight per count, $\\approx$0.86)"), ("C3", COL["c3"], "-", "C3 bump, conversion"),
                            ("C4", COL["c4"], "-", "C4 peak, $\\pi^+\\pi^-\\gamma$ conversion"), ("C5", COL["c5"], "-", "C5 peak, direct"),
                            ("Geta", COL["eta"], "-", "other $\\eta$ (G$_\\eta$)"), ("GDom", COL["gd"], "--", "G$_D$, $\\omega/\\rho^0$ pions"),
                            ("GDnom", COL["gd"], "-", "G$_D$, other pions"), ("GCom", COL["gc"], "--", "G$_C$, $\\omega/\\rho^0$ pions"), ("GCnom", COL["gc"], "-", "G$_C$, other pions")]:
    y = FD["sim_" + gname]
    ax.step(x[m], np.maximum(y[m], 0.3), where="mid", color=col, ls=ls, lw=1.8, label=lab)
ax.set_yscale("log"); ax.set_ylim(1, 1e4)
ax.axvline(0, color="k", lw=1); ax.axvspan(-45, 45, color="#eeeeee", zorder=0)
ax.set_xlabel("$u$ [MeV/$c^2$]"); ax.set_ylabel("expected counts per 2.5 MeV/$c^2$ (SMASH)")
ax.set_title("(b) SMASH, S0, mode A, $w\\times$V3, $M_4'=0.99615\\,M_4$")
ax.legend(loc="lower left", fontsize=9.5, ncol=2)
fig.tight_layout(); save(fig, "edgeu_spectra")

# ---------------------------------------------------------------- 3. the fit
fig, axs = plt.subplots(2, 1, figsize=(10.5, 10.5), gridspec_kw=dict(height_ratios=[1.25, 1]), sharex=True)
idx = FD["fit_idx"]; xi = x[idx]
ax = axs[0]
yd = FD["sig"][idx] - FD["fit_fx"]
ax.errorbar(xi, yd, yerr=np.sqrt(FD["sigvar"][idx]), fmt="o", ms=4, color=COL["data"], label="data $N-$CB$-a_P\\,K$ (S0)")
ax.plot(xi, FD["fit_model"], "-", color=COL["fit"], lw=2.2, label="fit: pol2 + $h\\times$step ($\\chi^2$/ndf %.1f/%d)" % (A["nom"]["chi2"], A["nom"]["ndf"]))
ax.plot(xi, FD["fit_smooth"], "--", color=COL["fit"], lw=1.6, label="the same without the step")
ax.fill_between(xi, 0, FD["fit_fx"], color=COL["c4"], alpha=0.35, step="mid", label="peak term $a_P\\,K$ (subtracted), $a_P$ = 0.418")
ax.axvline(0, color="k", lw=1)
ax.set_ylabel("counts per 2.5 MeV/$c^2$"); ax.set_ylim(0, 1.15 * yd.max())
ax.set_title("$h_{\\mathrm{data}}$ = %.0f $\\pm$ %.0f counts per 2.5 MeV/$c^2$, step significance $\\Delta\\chi^2$ = %.0f" % (A["nom"]["hd"], A["nom"]["ehd"], A["nom"]["chi20"] - A["nom"]["chi2"]))
ax.legend(loc="lower left", fontsize=10.5)
ax = axs[1]
rd = yd - FD["fit_smooth"]
rs = (FD["sim_y"] - FD["sim_fit_smooth"]) * f0
ax.axhline(0, color="k", lw=0.8, ls=":")
ax.errorbar(xi, rd, yerr=np.sqrt(FD["sigvar"][idx]), fmt="o", ms=4, color=COL["data"], label="data $-$ fitted smooth part")
ax.step(xi, rs, where="mid", color=COL["c2"], lw=2, label="SMASH bump C2+C3 $-$ its smooth part, $\\times f_\\eta$ = %.3f" % f0)
ax.plot(xi, A["nom"]["hd"] * FD["T"][idx], "-", color=COL["fit"], lw=2.2, label="$h_{\\mathrm{data}}\\times$step ($\\sigma$ = %.1f MeV/$c^2$)" % (A["reso"]["sigma_core"] * 1e3))
ax.axvline(0, color="k", lw=1)
ax.set_xlabel("$u=M_4-E(M_{\\pi\\pi},M_{ee})$ [MeV/$c^2$]"); ax.set_ylabel("counts per 2.5 MeV/$c^2$")
ax.legend(loc="lower left", fontsize=10.5)
fig.tight_layout(); save(fig, "edgeu_fit")

# ---------------------------------------------------------------- 4. checks
CH = A["chk"]
fig, axs = plt.subplots(2, 2, figsize=(15, 10.5))
ax = axs[0, 0]
bed = [(0.342, 0.044), (0.459, 0.039), (0.376, 0.027), (0.389, 0.062)]
for i, s in enumerate(CH["slices"]):
    s0 = A0["chk"]["slices"][i]     # uncorrected reading: original run (D38 does not change it)
    assert abs(s0["hd"] - s["hd"]) < 1e-9
    fu = s0["hd"] / (s0["lin"]["C2"] + s0["lin"]["C3"])
    ax.errorbar(i - 0.12, s["f"], yerr=s["statboot"], fmt="o", color=COL["data0"], capsize=3, label="$u$, C2 corrected (stat)" if i == 0 else None)
    ax.errorbar(i - 0.12, s["f"], yerr=math.hypot(s["statboot"], min(s["mc"], 0.3)), fmt="none", color=COL["data0"], lw=0.8)
    ax.errorbar(i + 0.0, fu, yerr=s0["statboot"] * fu / s0["f"], fmt="s", mfc="none", color=COL["data"], capsize=3, label="$u$, uncorrected (stat)" if i == 0 else None)
    ax.errorbar(i + 0.12, bed[i][0], yerr=bed[i][1], fmt="D", mfc="none", color=COL["bedge"], capsize=3, label="$M_4$ per slice (bedge, uncorr., stat)" if i == 0 else None)
ax.axhline(f0, color=COL["data0"], lw=1, ls="--"); ax.axhline(A0["nom"]["f_unc"], color=COL["data"], lw=1, ls=":")
ax.set_xticks(range(4)); ax.set_xticklabels(["0.30$-$0.34", "0.34$-$0.37", "0.37$-$0.39", "0.39$-$0.40"])
ax.set_xlabel("$M(\\pi^+\\pi^-)$ slice [GeV/$c^2$]"); ax.set_ylabel("$f_\\eta$"); ax.set_ylim(0.15, 0.85)
ax.set_title("(a) per slice (thin bars: stat $\\oplus$ MC; $\\chi^2$ %.1f/3)" % CH["slices_chi2"])
ax.legend(loc="upper left", fontsize=10)
ax = axs[0, 1]
labs = ["all", "E1 ($<$10)", "E2 ($\\geq$10)", "$\\geq$20 MeV/$c^2$"]
vals = [(f0, A["nom"]["statboot"], A["nom"]["mc"])] + [(CH["classes"][c]["f"], CH["classes"][c]["statboot"], CH["classes"][c]["mc"]) for c in ("E1", "E2", "E20")]
for i, (v, s_, mc) in enumerate(vals):
    ax.errorbar(i, v, yerr=s_, fmt="o", color=COL["data0"], capsize=4, lw=2)
    ax.errorbar(i, v, yerr=math.hypot(s_, mc), fmt="none", color=COL["data0"], lw=0.9)
ax.errorbar([1.15, 2.15], [0.3859, 0.4870], yerr=[0.0231, 0.0462], fmt="D", mfc="none", color=COL["bedge"], capsize=3, label="$M_4$ per slice (bedge, C2 corr., stat)")
ax.set_xticks(range(4)); ax.set_xticklabels(labs); ax.set_ylabel("$f_\\eta$ (C2 corrected)"); ax.set_ylim(0.2, 0.65)
ax.set_title("(b) per $M(e^+e^-)$ class; $\\kappa$ from E1/E2: %.2f $\\pm$ %.2f $\\pm$ %.2f" % (CH["kappa_free"]["kappa"], CH["kappa_free"]["kappa_stat"], CH["kappa_free"]["kappa_mc"]))
ax.legend(loc="upper right", fontsize=10)
ax = axs[1, 0]
labs = ["I", "IIA", "IIB", "III", "fold 0", "fold 1", "pre", "isBest=1"]
vals = [(p["fc"], p["statboot"]) for p in CH["periods"]] + [(fo["f"], fo["statboot"]) for fo in CH["folds"]] + \
       [(A["var"]["pre"]["f"], A["var"]["pre"]["statboot"]), (A["var"]["best"]["f"], A["var"]["best"]["statboot"])]
for i, (v, s_) in enumerate(vals):
    ax.errorbar(i, v, yerr=s_, fmt="o", color=COL["data0"] if i < 4 else COL["data"], capsize=4)
ax.axhline(f0, color="k", lw=1, ls="--")
ax.set_xticks(range(len(labs))); ax.set_xticklabels(labs); ax.set_ylabel("$f_\\eta$ (C2 corrected; periods $\\times C_{\\mathrm{sec}}$)")
ax.set_title("(c) periods ($\\chi^2$ %.1f/3, p = %.3f), folds (z = %.2f), selections (stat)" % (CH["periods_chi2"], CH["periods_p"], CH["folds_z"]))
ax = axs[1, 1]
nulls = [("like-sign $e^+e^+$", CH["null_ls_pp_all"]["h"], CH["null_ls_pp_all"]["eh"], COL["cb"]),
         ("like-sign $e^-e^-$", CH["null_ls_mm_all"]["h"], CH["null_ls_mm_all"]["eh"], COL["cb"]),
         ("$e^-e^-$, E1", CH["null_ls_mm_E1"]["h"], CH["null_ls_mm_E1"]["eh"], COL["cb"]),
         ("$e^-e^-$, E2", CH["null_ls_mm_E2"]["h"], CH["null_ls_mm_E2"]["eh"], COL["cb"]),
         ("CB (subtracted)", CH["null_cb"]["h_lin"], 0.0, COL["cb"]),
         ("SMASH G$_D$", CH["null_sim"]["GD"]["h"], CH["null_sim"]["GD"]["rep_rms"], COL["gd"]),
         ("SMASH G$_C$", CH["null_sim"]["GC"]["h"], CH["null_sim"]["GC"]["rep_rms"], COL["gc"]),
         ("SMASH G$_\\eta$", CH["null_sim"]["Geta"]["h"], CH["null_sim"]["Geta"]["rep_rms"], COL["eta"]),
         ("data, $M_{\\pi\\pi}$ 0.42$-$0.46\n(pseudo edge)", CH["control"]["hd"], CH["control"]["ehd"], COL["data0"])]
for i, (lab, v, e_, col) in enumerate(nulls):
    ax.errorbar(v, i, xerr=e_, fmt="o", color=col, capsize=4)
ax.axvline(0, color="k", lw=1)
ax.axvspan(-0.05 * A["nom"]["hd"], 0.05 * A["nom"]["hd"], color="#eeeeee", zorder=0)
ax.set_yticks(range(len(nulls))); ax.set_yticklabels([n[0] for n in nulls], fontsize=10.5); ax.invert_yaxis()
ax.set_xlabel("step height at $u=0$ [counts / 2.5 MeV/$c^2$]; grey: $\\pm$5 % of $h_{\\mathrm{data}}$" )
ax.set_title("(d) null tests (SMASH: replica rms)")
fig.tight_layout(); save(fig, "edgeu_checks")

# ---------------------------------------------------------------- 5. systematics
rows = [r for r in NB["rows"]]
# D40: the peak-tail term t = 1 +- 0.5 enters U8 by the family maximum (edgeu_combine_v3.py); f rises with t
rows += [dict(fam="U8", tag="peak tail $t=0.5$ (D40)", d=-NB["tail_rel"] * f0, noise_boot=0.0, counted=True),
         dict(fam="U8", tag="peak tail $t=1.5$ (D40)", d=+NB["tail_rel"] * f0, noise_boot=0.0, counted=True)]
fams = ["U1", "U2", "U3", "U4", "U5", "U6", "U7", "U8", "U9", "U10", "chk"]
cmap = plt.get_cmap("tab10")
fig, ax = plt.subplots(figsize=(10.5, 15))
yy = 0; ticks = []; labels = []
for fi, fam in enumerate(fams):
    for r in rows:
        if r["fam"] != fam:
            continue
        col = cmap(fi % 10)
        ax.errorbar(r["d"] / f0 * 100, yy, xerr=r["noise_boot"] / f0 * 100, fmt="o" if r["counted"] else "x", color=col, ms=6, capsize=2)
        ticks.append(yy); labels.append("%s %s" % (fam, r["tag"])); yy += 1
    yy += 0.6
ax.axvline(0, color="k", lw=1)
ax.set_yticks(ticks); ax.set_yticklabels(labels, fontsize=9.5); ax.invert_yaxis()
ax.set_xlabel("shift of $f_\\eta$(u edge) [%] (bars: bootstrap $\\oplus$ replica noise of the shift; x = not counted)")
ax.set_title("variants: total syst %.1f %% (U2 resolution %.1f %%, U8 peak and tail %.1f %%,\nU11 background step %.1f %%)" % (100 * NB["syst"] / f0, 100 * NB["fams"]["dl"]["U2"] / f0, 100 * NB["fams"]["dl"]["U8"] / f0, 100 * NB["fams"]["dl"]["U11"] / f0))
ax.set_xlim(-15, 12)
fig.tight_layout(); save(fig, "edgeu_syst")

# ---------------------------------------------------------------- 6. toys
order = ["a", "c1", "a1", "null", "q-", "q--", "q+", "om", "ty", "sA", "sm3", "sm5", "sm6", "datalike", "kap", "pk30", "pk55", "mass", "raw", "narrow", "wide"]
bedge_shift = {"q-": -2.1, "q--": -4.1, "q+": 2.2, "om": -1.7, "wide": 19.9, "narrow": 9.1}
fig, axs = plt.subplots(1, 2, figsize=(15.5, 8.5), gridspec_kw=dict(width_ratios=[2.2, 1]))
ax = axs[0]
for i, k in enumerate(order):
    o = TY[k]
    if k == "null":
        v1, v2 = 100 * o["mean_aa"] / 0.42, 100 * o["mean_rf"] / 0.42
        e1, e2 = 100 * o["rms_aa"] / 0.42, 100 * o["rms_rf"] / 0.42
    else:
        v1, v2 = o["bias_aa"], o["bias_rf"]; e1 = 100 * o["rms_aa"] / o["ftrue"]; e2 = 100 * o["rms_rf"] / o["ftrue"]
    ax.errorbar(v1, i - 0.15, xerr=e1, fmt="o", color=COL["data"], capsize=2, label="templates = truth sample (AA)" if i == 0 else None)
    ax.errorbar(v2, i + 0.15, xerr=e2, fmt="s", mfc="none", color=COL["c2"], capsize=2, label="templates = replica (RF)" if i == 0 else None)
    if k in bedge_shift:
        ax.plot(TY["a"]["bias_aa"] + bedge_shift[k], i, "D", color=COL["bedge"], ms=7, label="per-slice $M_4$ edge (bedge), shift vs (a)" if k == "q-" else None)
ax.axvspan(-5, 5, color="#eeeeee", zorder=0); ax.axvline(0, color="k", lw=1)
ax.set_yticks(range(len(order))); ax.set_yticklabels(order); ax.invert_yaxis()
ax.set_xlabel("bias of $f_\\eta$(u edge) [%] (null: $f$/0.42 in %; bars: toy rms)")
ax.set_title("(a) pseudo-data closure, 1000 toys per truth")
ax.text(0.98, 0.45, "B4 with q = $-$0.173: +150 %%\nu edge: %+.1f %% (vs (a))" % TY["q-"]["shift_vs_a_asimov"], transform=ax.transAxes, ha="right", fontsize=11,
        bbox=dict(fc="white", ec=COL["b4"]))
ax.legend(loc="upper right", fontsize=10)
ax = axs[1]
for i, k in enumerate(order):
    o = TY[k]
    ax.plot(o["pull_sd_aa"], i - 0.15, "o", color=COL["data"]); ax.plot(o["pull_sd_rf"], i + 0.15, "s", mfc="none", color=COL["c2"])
ax.axvline(1, color="k", lw=1)
ax.set_yticks(range(len(order))); ax.set_yticklabels(order); ax.invert_yaxis(); ax.set_xlim(0.4, 1.3)
ax.set_xlabel("pull SD (AA: stat; RF: stat $\\oplus$ MC)"); ax.set_title("(b) pull widths")
fig.tight_layout(); save(fig, "edgeu_toys")

# ---------------------------------------------------------------- 7. comparison of routes (sigma)
S_ = NB["sigma"]; SU = NB["sigma_unc"]; SK = NB["sigma_k48"]; cmpv = NB["cmp"]
routes = [("u edge, C2 corrected (this work)", S_["v"], S_["stat"], math.sqrt(S_["stat"] ** 2 + S_["syst"] ** 2 + S_["mc"] ** 2), COL["data0"]),
          ("u edge, uncorrected", SU["v"], SU["stat"], math.sqrt(SU["stat"] ** 2 + SU["syst"] ** 2 + SU["mc"] ** 2), COL["data"]),
          ("u edge, C2 corr., conversions x 0.48", SK["v"], SK["stat"], math.sqrt(SK["stat"] ** 2 + SK["syst"] ** 2 + SK["mc"] ** 2), COL["c3"]),
          ("per-slice M4 edge (bedge), C2 corr.", 0.442, 0.024, math.sqrt(0.024 ** 2 + 0.040 ** 2 + 0.050 ** 2), COL["bedge"]),
          ("peak side bands (Part III)", 0.429, 0.036, math.sqrt(0.036 ** 2 + 0.157 ** 2 + 0.025 ** 2), COL["c4"]),
          ("peak, direct-decay route", 0.522, 0.062, 0.062, COL["c5"]),
          ("peak, conversion route", 0.249, 0.051, 0.051, COL["c4"]),
          ("F-main (cross-check)", 0.399, None, None, COL["gd"]),
          ("B4 (not determined)", NB["b4_sigma"], None, None, COL["b4"]),
          ("ECAL tag (verification only)", 0.474, 0.045, 0.045, COL["tag"])]
fig, ax = plt.subplots(figsize=(11.5, 7.2))
for i, (lab, v, es, et, col) in enumerate(routes):
    if es is None:
        ax.plot(v, i, "o", mfc="none", color=col, ms=8)
    else:
        ax.errorbar(v, i, xerr=et, fmt="none", color=col, lw=1.2, capsize=3)
        ax.errorbar(v, i, xerr=es, fmt="o", color=col, lw=3.5, ms=7)
ax.set_yticks(range(len(routes))); ax.set_yticklabels([r[0] for r in routes]); ax.invert_yaxis()
ax.axvline(1.069, color=COL["c2"], lw=1.2, ls="--"); ax.text(1.06, len(routes) - 0.6, "SMASH 1.069", ha="right", color=COL["c2"], fontsize=11)
ax.set_xlim(0.1, 1.15)
ax.set_xlabel("$\\sigma(pp\\to\\eta X)$ [mb] (PDG 2026; thick: stat, thin: stat $\\oplus$ syst $\\oplus$ MC; no normalisation)")
ax.set_title("$\\sigma(pp\\to\\eta X)$ at T = 4.5 GeV: the bump in u against the other routes")
fig.tight_layout(); save(fig, "edgeu_compare")
print("figures written: edgeu_{fit,spectra,checks,syst,compare}_v3")
