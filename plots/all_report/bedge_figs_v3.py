#!/usr/bin/env python3
# ============================================================================
# bedge_figs.py  --  package bedge: figures (Part III colours: N - CB blue, C2 orange, C3 magenta, fit green, scale estimates
# as in fix_figs.py)
# ============================================================================
# Written 2026-11-14 (prefix bedge). Inputs: /hdd3/.../bedge/bedge_figdata.json (bedge_ana.C) and bedge_numbers.json
# (bedge_combine.py). Outputs docs/latex/figs_all/bedge_{steps,scales,sigma,syst,toys}.{pdf,png}.
# Run: python3 -I plots/all_report/bedge_figs.py
# ----------------------------------------------------------------------------
# _v3 copy (2026-11-15, FIGURE AGENT C of the fourth integration): rulings D38 (C2 gamma* angle corrected with the per-count
# QED weights of c2recon) and D40 (peak-tail family t = 1 +- 0.5, counted). Differences to bedge_figs.py:
#   inputs  bedge_v3/bedge_figdata_v3.json (identical to bedge_figdata.json) and bedge_v3/bedge_numbers_v3.json
#           (bedge_ana_v3.C, bedge_combine_v3.py): fc_rec, a_k48, the E1/E2 split and all sigma values with D38 + D40;
#   values  (D40, not in bedge_numbers_v3.json) the total of a_P in bedge_scales gains the peak-tail term 8.75 %
#           (syst sqrt(35.5^2 + 8.75^2) = 36.6 %), and the syst of sigma(peak) in bedge_sigma is 0.157 mb instead of 0.152
#           (tables_all/sigma_headline_v3.tex, p3w_tableS1_v3.tex).
# Written: docs/latex/figs_all/bedge_{scales,sigma}_v3.{pdf,png}. bedge_steps (uncorrected f; only its 'rows' label changes,
# see bedge_figs_cnt.py), bedge_syst (f_all variants) and bedge_toys do not change and are not written.
# Run: python3 -I plots/all_report/bedge_figs_v3.py
# ============================================================================
import json
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = "/hdd3/przygoda/HADES/FAT_out/all_report/bedge/"
FIG = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_all/"
plt.rcParams.update({"font.size": 13, "axes.labelsize": 14, "axes.titlesize": 13.5, "legend.fontsize": 11, "xtick.labelsize": 12,
                     "ytick.labelsize": 12, "pdf.fonttype": 42})
PAL = dict(main="#1f3fbf", direct="#00a0a0", conv="#009900", fmain="#555555", b4="#ff6a00", edge="#8c564b", tag="#9467bd", pm="#d400d4",
           lit="#e8202a", smash="#ff8c00", data="#1f3fbf", fit="#009900", c2="#ff8c00", c3="#d400d4")
D3 = "/hdd3/przygoda/HADES/FAT_out/all_report/bedge_v3/"
FD = json.load(open(D3 + "bedge_figdata_v3.json"))
NB = json.load(open(D3 + "bedge_numbers_v3.json"))
PDG26 = 1.00436

def save(fig, name):
    if name not in ("bedge_scales", "bedge_sigma"):
        plt.close(fig); return
    fig.savefig(FIG + name + "_v3.pdf"); fig.savefig(FIG + name + "_v3.png", dpi=130); plt.close(fig)

# ---------------------------------------------------------------- 1. edge step per slice, data vs SMASH (scaled)
f = FD["f"]
fig, axs = plt.subplots(2, 4, figsize=(17, 8.6), sharex=False)
for i, s in enumerate(FD["slices"]):
    tx, fy, fy0, ty = map(np.array, (s["tx"], s["fy"], s["fy0"], s["ty"]))
    x, y, ey, ys, y23 = map(np.array, (s["x"], s["y"], s["ey"], s["ys"], s["y23"]))
    sfy, sfy0 = np.array(s["sfy"]), np.array(s["sfy0"])
    m = (x >= tx.min() - 1e-6) & (x <= tx.max() + 1e-6)
    ax = axs[0, i]
    ax.errorbar(x[m] * 1e3, y[m], yerr=ey[m], fmt="o", ms=3.5, color=PAL["data"], label="data $N-$CB")
    ax.plot(tx * 1e3, fy, "-", color=PAL["fit"], lw=2, label="fit: pol2 [+$\\eta$] + $h\\times$step")
    ax.plot(tx * 1e3, fy0, "--", color=PAL["fit"], lw=1.5, label="the same without the step")
    ax.axvspan(s["edgeLo"] * 1e3, s["edgeHi"] * 1e3, color="#cccccc", alpha=0.5, lw=0)
    ax.set_title("$M(\\pi^+\\pi^-)$ %.2f$-$%.2f GeV/$c^2$" % (s["lo"], s["hi"]))
    ax.text(0.04, 0.05 if i < 3 else 0.80, "$h_{\\mathrm{data}}=%.0f\\pm%.0f$\n$\\chi^2$/ndf %.1f/%d" % (s["h"], s["eh"], s["chi2"], s["ndf"]), transform=ax.transAxes, fontsize=11.5)
    if i == 0:
        ax.set_ylabel("rows per 2.5 MeV/$c^2$")
    if i == 3:
        ax.legend(loc="lower right", fontsize=9.5)
    ax = axs[1, i]
    # residual step view: data minus the fitted smooth part; SMASH T1+T2 minus its own fitted smooth part, x f
    idx = [int(np.argmin(np.abs(x - t))) for t in tx]
    rd, ed = y[idx] - fy0, ey[idx]
    rs = (ys[idx] - sfy0) * f
    ax.axhline(0, color="k", lw=0.8, ls=":")
    ax.errorbar(tx * 1e3, rd, yerr=ed, fmt="o", ms=3.5, color=PAL["data"], label="data $-$ fitted smooth part")
    ax.step(tx * 1e3, rs, where="mid", color=PAL["smash"], lw=1.8, label="SMASH T1+T2 $-$ its smooth part, $\\times f_\\eta$=%.3f" % f)
    ax.plot(tx * 1e3, s["h"] * ty, "-", color=PAL["fit"], lw=2, label="$h_{\\mathrm{data}}\\times$step")
    ax.plot(tx * 1e3, f * s["hs"] * ty, "--", color=PAL["smash"], lw=1.5, label="$f_\\eta\\,h_{\\mathrm{sim}}\\times$step")
    ax.axvspan(s["edgeLo"] * 1e3, s["edgeHi"] * 1e3, color="#cccccc", alpha=0.5, lw=0)
    ax.set_xlabel("$M(\\pi^+\\pi^-e^+e^-)$ [MeV/$c^2$]")
    ax.text(0.04, 0.05, "$h_{\\mathrm{data}}/h_{\\mathrm{sim}}=%.3f$" % (s["h"] / s["hs"]), transform=ax.transAxes, fontsize=11.5)
    if i == 0:
        ax.set_ylabel("step view [rows per 2.5 MeV/$c^2$]")
    if i == 3:
        ax.legend(loc="upper right", fontsize=9.0)
fig.suptitle("Edge step of the $\\eta\\to\\pi^+\\pi^-\\pi^0$ bump: data S0 (isBest$\\geq$0, L = 6079.2 nb$^{-1}$) and SMASH scaled by "
             "$f_\\eta$(edge) = %.3f; grey: the edge range of the slice" % f, fontsize=13)
fig.tight_layout(rect=(0, 0, 1, 0.96))
save(fig, "bedge_steps")

# ---------------------------------------------------------------- 2. eta-scale comparison and kappa
R = NB["RES"]; S = NB["SIGMA"]; C = NB["CMP"]
def tot(x):
    return math.sqrt(x["stat"] ** 2 + x["syst"] ** 2 + x["mc"] ** 2)
rows = [
    ("bump edge, SMASH as is ($\\kappa$=1)", R["f_all"]["v"] / PDG26, R["f_all"]["stat"] / PDG26, tot(R["f_all"]) / PDG26, PAL["edge"], "o", True),
    ("bump edge, C2 angle corrected", R["fc_rec"]["v"] / PDG26, R["fc_rec"]["stat"] / PDG26, tot(R["fc_rec"]) / PDG26, PAL["edge"], "s", True),
    ("bump edge, C2 corr., conversions $\\times$0.48", R["a_k48"]["v"] / PDG26, R["a_k48"]["stat"] / PDG26, tot(R["a_k48"]) / PDG26, PAL["edge"], "D", True),
    ("bump edge, E1/E2 split ($\\kappa$ free)", NB["nomq"]["a_D"]["v"] / PDG26, NB["nomq"]["a_D"]["statboot"] / PDG26,
     math.hypot(NB["nomq"]["a_D"]["statboot"], NB["nomq"]["a_D"]["mc68"]) / PDG26, PAL["edge"], "^", False),
    ("ECAL tag (verification only)", C["tag_scale"][0] / PDG26, C["tag_scale"][1] / PDG26, math.hypot(C["tag_scale"][1], C["tag_scale"][2]) / PDG26, PAL["tag"], "o", False),
    ("peak side bands $a_P$ (MAIN, $\\kappa$=1)", C["aP"][0], C["aP"][1], C["aP"][0] * math.sqrt(0.085 ** 2 + 0.355 ** 2 + 0.0875 ** 2 + 0.059 ** 2), PAL["main"], "o", True),   # D40: + peak tail 8.75 %
    ("peak, direct-decay route (C5)", C["direct_scale"][0], C["direct_scale"][1], None, PAL["direct"], "D", True),
    ("peak, conversion route (C4)", C["conv_scale"][0], C["conv_scale"][1], None, PAL["conv"], "s", True),
    ("F-main $a_\\eta$ (cross-check)", 0.3806, 0.0169, 0.3806 * math.sqrt(0.044 ** 2 + 0.555 ** 2 + 0.104 ** 2), PAL["fmain"], "o", True),
    ("bump B4 $a_B$(E2): NOT DETERMINED", 0.3103, 0.018, None, PAL["b4"], "o", True)]
fig, axs = plt.subplots(1, 2, figsize=(16.5, 7.4), gridspec_kw=dict(width_ratios=[1.55, 1]))
ax = axs[0]
n = len(rows)
for k, (lab, v, s, t, col, mk, filled) in enumerate(rows):
    yv = n - k
    if t is not None:
        ax.errorbar(v, yv, xerr=t, fmt="none", ecolor=col, elinewidth=1.2)
    ax.errorbar(v, yv, xerr=s, fmt=mk, color=col, ms=8, elinewidth=3.5, mfc=col if filled else "white")
    if "NOT" in lab:
        ax.annotate("", xy=(0.97, yv), xytext=(0.64, yv), arrowprops=dict(arrowstyle="<->", color=col, lw=1.2))
        ax.text(0.80, yv + 0.25, "syst 165 %", color=col, ha="center", fontsize=10)
ax.axvline(1.0, color="k", ls="--", lw=1); ax.text(0.99, n + 0.55, "SMASH", ha="right", fontsize=10)
ax.axhline(n - 4.5, color="#999999", lw=0.8)
ax.set_yticks(range(n, 0, -1)); ax.set_yticklabels([r[0] for r in rows], fontsize=11.5)
ax.set_xlim(0, 1.05); ax.set_ylim(0.3, n + 1.0)
ax.set_xlabel("data / SMASH ($\\sigma\\times$PDG 2026 BR), no $C_{\\mathrm{sec}}(1-f_\\eta)$")
ax.set_title("(a) $\\eta$ scale relative to SMASH (thick: stat; thin: total)")
ax = axs[1]
krows = [("bump edge, E1/E2 split", NB["nomq"]["kappa"]["v"], NB["nomq"]["kappa"]["statboot"], math.hypot(NB["nomq"]["kappa"]["statboot"], NB["nomq"]["kappa"]["mc68"]), PAL["edge"], "^"),
         ("peak: conversion / direct route", 0.477, 0.11, None, PAL["main"], "o"),
         ("peak: side bands per $M_{ee}$ class", 0.37, 0.23, None, PAL["direct"], "d"),
         ("peak: X17 2D fit $N_4/N_5$", 0.51, None, None, PAL["direct"], "s"),
         ("F-main E1+E2", 0.504, 0.162, None, PAL["fmain"], "o")]
for k, (lab, v, s, t, col, mk) in enumerate(krows):
    yv = len(krows) - k
    if t is not None:
        ax.errorbar(v, yv, xerr=t, fmt="none", ecolor=col, elinewidth=1.2)
    ax.errorbar(v, yv, xerr=s if s else None, fmt=mk, color=col, ms=8, elinewidth=3.5)
ax.axvline(1.0, color="k", ls="--", lw=1); ax.text(0.98, len(krows) + 0.55, "SMASH/HGeant", ha="right", fontsize=10)
ax.axvline(0.48, color=PAL["edge"], ls=":", lw=1.2); ax.text(0.49, 0.45, "0.48 used for\nthe $\\kappa$-scaled edge", color=PAL["edge"], fontsize=10)
ax.set_yticks(range(len(krows), 0, -1)); ax.set_yticklabels([r[0] for r in krows], fontsize=11.5)
ax.set_xlim(-0.3, 1.6); ax.set_ylim(0.2, len(krows) + 1.0)
ax.set_xlabel("$\\kappa$ (conversion vs Dalitz/direct pairs)")
ax.set_title("(b) conversion scale $\\kappa$")
fig.tight_layout()
save(fig, "bedge_scales")

# ---------------------------------------------------------------- 3. sigma comparison
srows = [
    ("bump edge, SMASH as is", S["f_all"], PAL["edge"], "o", True),
    ("bump edge, C2 angle corrected", S["fc_rec"], PAL["edge"], "s", True),
    ("bump edge, C2 corr., conversions $\\times$0.48", S["a_k48"], PAL["edge"], "D", True)]
fig, ax = plt.subplots(figsize=(13.5, 7.2))
lab, vals = [], []
items = []
for name, s, col, mk, fl in srows:
    items.append((name, s["v"], s["stat"], s["tot"], col, mk, fl))
items.append(("bump edge, E1/E2 split ($\\kappa$ free; stat, MC68)", S["a_D"]["v"], S["a_D"]["stat"], math.hypot(S["a_D"]["stat"], S["a_D"]["mc68"]), PAL["edge"], "^", False))
items.append(("ECAL tag (verification only)", 0.474, 0.474 * 0.014 / 0.454, 0.474 * math.hypot(0.014, 0.040) / 0.454, PAL["tag"], "o", False))
pk = list(C["peak"])
assert abs(pk[2] - 0.152) < 1e-9
pk[2] = 0.157     # D40: syst 0.1523 (+) 0.0375 (peak tail) = 0.1569 mb (sigma_headline_v3.tex)
items.append(("peak side bands (adopted, PRELIMINARY)", pk[0], pk[1], math.sqrt(pk[1] ** 2 + pk[2] ** 2 + pk[3] ** 2 + pk[4] ** 2), PAL["main"], "o", True))
items.append(("peak, direct-decay route (stat)", C["direct"][0], C["direct"][1], None, PAL["direct"], "D", True))
items.append(("peak, conversion route (stat)", C["conv"][0], C["conv"][1], None, PAL["conv"], "s", True))
items.append(("F-main (cross-check)", C["fmain_sigma"], 0.018, 0.399 * 0.57, PAL["fmain"], "o", True))
items.append(("HADES 3.5 GeV inclusive, PYTHIA (lower energy)", 1.14, 0.20, None, "#888888", "v", False))
items.append(("HADES 3.5 GeV inclusive, UrQMD (lower energy)", 0.93, 0.14, None, "#888888", "v", False))
n = len(items)
for k, (name, v, s, t, col, mk, fl) in enumerate(items):
    yv = n - k
    if t is not None:
        ax.errorbar(v, yv, xerr=t, fmt="none", ecolor=col, elinewidth=1.2)
    ax.errorbar(v, yv, xerr=s, fmt=mk, color=col, ms=8, elinewidth=3.5, mfc=col if fl else "white")
    ax.text(1.62, yv, "%.3f" % v, va="center", fontsize=10.5, color=col)
ax.axvline(1.0694, color=PAL["smash"], ls="--", lw=1.2); ax.text(1.075, n + 0.55, "SMASH 1.069", color=PAL["smash"], fontsize=10.5)
ax.axhline(n - 8.5, color="#999999", lw=0.8)
ax.set_yticks(range(n, 0, -1)); ax.set_yticklabels([it[0] for it in items], fontsize=11.5)
ax.set_xlim(0, 1.72); ax.set_ylim(0.3, n + 1.0)
ax.set_xlabel("$\\sigma(pp\\to\\eta X)$ at $T$ = 4.5 GeV [mb]  (PDG 2026, $C_{\\mathrm{sec}}(1-f_\\eta)$ applied)")
ax.set_title("$\\sigma(pp\\to\\eta X)$ from the bump edge and the other routes (thick: stat; thin: stat $\\oplus$ syst $\\oplus$ MC $\\oplus$ norm.)")
fig.tight_layout()
save(fig, "bedge_sigma")

# ---------------------------------------------------------------- 4. systematics: variant shifts, per-slice ratios, periods
fig, axs = plt.subplots(1, 3, figsize=(18, 7.6), gridspec_kw=dict(width_ratios=[1.5, 1, 1]))
ax = axs[0]
rws = [r for r in NB["rows"]["f_all"] if r["fam"] != "CHK"]
rws.sort(key=lambda r: (r["fam"], r["tag"]))
famcol = {"E1": "#1f3fbf", "E2": "#009900", "E3": "#ff6a00", "E4": "#d62728", "E5": "#9467bd", "E6": "#555555"}
for k, r in enumerate(rws):
    yv = len(rws) - k
    ax.barh(yv, r["lit"], left=-r["lit"], height=0.7, color="#dddddd", lw=0)
    ax.plot(r["d"], yv, "o" if r["cnt"] else "x", color=famcol.get(r["fam"], "k"), ms=7, mfc=famcol.get(r["fam"], "k") if r["cnt"] else "none")
ax.axvline(0, color="k", lw=0.8)
ax.set_yticks(range(len(rws), 0, -1)); ax.set_yticklabels(["%s %s" % (r["fam"], r["tag"]) for r in rws], fontsize=9.5)
ax.set_xlabel("$\\Delta f_\\eta$(edge) = variant $-$ nominal")
ax.set_title("(a) variants: shift (points; x = not counted) and its noise (grey)")
ax.set_ylim(0.3, len(rws) + 0.7)
ax = axs[1]
mid = [0.32, 0.355, 0.38, 0.395]
for tag, col, mk, off in (("nominal", PAL["data"], "o", -0.0015), ("bump1aY", PAL["conv"], "s", 0.0), ("bumpPS", PAL["b4"], "^", 0.0015)):
    sh = NB["slice_shape"][tag]
    v = [a for a, _ in sh["r"]]; e = [b for _, b in sh["r"]]
    ax.errorbar(np.array(mid) + off, v, yerr=e, fmt=mk, color=col, ms=7,
                label={"nominal": "KLOE-2 Dalitz plot (nominal)", "bump1aY": "1 + aY", "bumpPS": "phase space"}[tag] + ": $\\chi^2$ %.1f/3" % sh["chi2"])
ax.axhline(NB["RES"]["f_all"]["v"], color=PAL["data"], ls=":", lw=1)
ax.set_xlabel("$M(\\pi^+\\pi^-)$ slice centre [GeV/$c^2$]"); ax.set_ylabel("$h_{\\mathrm{data}}/h_{\\mathrm{sim}}$")
ax.set_title("(b) per-slice ratio (data statistical errors)"); ax.legend(fontsize=10, loc="upper left"); ax.set_ylim(0.1, 0.85)
ax = axs[2]
pr = NB["periods"]["f_all_1"]
xs = np.arange(4)
ax.errorbar(xs, [a for a, _ in pr["vals"]], yerr=[b for _, b in pr["vals"]], fmt="o", color=PAL["data"], ms=8, label="period, $\\times C_{\\mathrm{sec}}(p)$")
fo = NB["folds"]
ax.errorbar([4.6, 5.0], [fo["f0"], fo["f1"]], yerr=[fo["e0"], fo["e1"]], fmt="s", color=PAL["fmain"], ms=7, label="folds (z = %.2f)" % fo["z"])
ax.axhline(pr["mean"], color=PAL["data"], ls=":", lw=1)
ax.set_xticks(list(xs) + [4.8]); ax.set_xticklabels(["I", "IIA", "IIB", "III", "folds"])
ax.set_ylabel("$f_\\eta$(edge)"); ax.set_title("(c) periods: $\\chi^2$ %.1f/3, p = %.2f" % (pr["chi2"], pr["p"]))
ax.legend(fontsize=10, loc="upper right"); ax.set_ylim(0.15, 0.7)
fig.tight_layout()
save(fig, "bedge_syst")

# ---------------------------------------------------------------- 5. pseudo-data closure
T = NB["TOY"]
sc = [("a1", "SMASH bkg x1, a = 1", "f_all"), ("c1", "SMASH bkg x1, a = 0.4", "f_all"), ("a", "WP4 (a): bkg x1.35, tilt", "f_all"),
      ("om", "WP4 (e): omega/rho0 x0.3, rest x1.6", "f_all"), ("sA", "stage-A bkg, F-main eta", "f_all"),
      ("q-", "(a) + q = -0.173 (B4: +150 %)", "f_all"), ("q--", "(a) + q = -0.346", "f_all"), ("q+", "(a) + q = +0.173", "f_all"),
      ("smooth", "(a), smoothed bkg in the truth", "f_all"), ("wide", "(a) x (1 + 0.10 G(490, 15 MeV))", "f_all"),
      ("narrow", "(a) x (1 + 0.05 G(505, 8 MeV))", "f_all"), ("mass", "truth s + 0.0009, r = 0.89", "f_all"),
      ("qed", "truth C2 x 0.8855 -> C2-corrected", "fc_rec"), ("kap", "truth conversions x 0.48 -> kappa 0.48", "a_k48_unc"),
      ("datalike", "a 0.55, C2 x 0.8855, conv. x 0.48 -> kappa 0.48", "a_k48")]
fig, axs = plt.subplots(1, 2, figsize=(16, 7.6), gridspec_kw=dict(width_ratios=[1.6, 1]))
ax = axs[0]
for k, (s, lab, q) in enumerate(sc):
    yv = len(sc) - k
    aa, rf = T["%s|AA|%s" % (s, q)], T["%s|RF|%s" % (s, q)]
    ax.errorbar(100 * aa["rel"], yv + 0.12, xerr=100 * aa["rms"] / aa["truth"], fmt="o", color=PAL["main"], ms=7, capsize=0, elinewidth=1,
                label="AA: template = truth sample" if k == 0 else None)
    ax.errorbar(100 * rf["rel"], yv - 0.12, xerr=100 * rf["rms"] / rf["truth"], fmt="s", color=PAL["b4"], ms=6, mfc="white", capsize=0, elinewidth=1,
                label="RF: template = replica" if k == 0 else None)
ax.axvline(0, color="k", lw=0.8); ax.axvspan(-5, 5, color="#eeeeee", lw=0)
ax.text(39, len(sc) - 5, "B4 moved by +150 %\nfor this term (WP6)", color=PAL["b4"], fontsize=10.5, ha="right", va="center")
ax.set_yticks(range(len(sc), 0, -1)); ax.set_yticklabels([s[1] for s in sc], fontsize=10.5)
ax.set_xlim(-40, 40); ax.set_xlabel("relative bias of the estimator [%] (bars: rms of the toys)")
ax.set_title("(a) closure on pseudo-data, 1000 toys per truth"); ax.legend(fontsize=10, loc="upper right")
ax = axs[1]
for k, (s, lab, q) in enumerate(sc):
    yv = len(sc) - k
    aa, rf = T["%s|AA|%s" % (s, q)], T["%s|RF|%s" % (s, q)]
    ax.plot(aa["psd"], yv + 0.12, "o", color=PAL["main"], ms=7)
    ax.plot(rf["psd"], yv - 0.12, "s", color=PAL["b4"], ms=6, mfc="white")
ax.axvline(1, color="k", ls="--", lw=0.8)
ax.set_yticks(range(len(sc), 0, -1)); ax.set_yticklabels([""] * len(sc))
ax.set_xlim(0.5, 1.6); ax.set_xlabel("pull width (AA: stat; RF: stat $\\oplus$ MC)"); ax.set_title("(b) pull widths")
fig.tight_layout()
save(fig, "bedge_toys")
print("figures written: bedge_scales_v3, bedge_sigma_v3")
