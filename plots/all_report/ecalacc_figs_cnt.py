#!/usr/bin/env python3
# ============================================================================
# ecalacc_figs.py -- package ecalacc (2026-11-15): figures of the ECAL-photon feasibility study
# ============================================================================
# Part III colours: data black / blue, CB red, bump C2 orange, C3 magenta, peak C4 green, C5 cyan, other eta violet,
# non-eta greys. Inputs: /hdd3/przygoda/HADES/FAT_out/all_report/ecalacc/{ecalacc_numbers.json, geo_cells.json,
# hist_data.txt, hist_sim.txt, tag_fits.txt, rows_c23.bin}. Outputs docs/latex/figs_all/ecalacc_{geometry,fractions,
# energy,tag,closure}.{pdf,png}.
# Run: python3 -I plots/all_report/ecalacc_figs.py   (after ecalacc_fill.C, ecalacc_geo.py, ecalacc_ana.py)
# ============================================================================
# _cnt copy (2026-11-15, FIGURE AGENT D of the integration run): the analyst calls the stored pi+pi-e+e- combinations
# COUNTS, not rows. Only label, legend, title and annotation strings 'rows' -> 'counts' and the output names (suffix _cnt)
# differ from the original ecalacc_figs.py; the inputs and the plotted content are the same.
# ecalacc_geometry is not written (superseded by ecalacc_geometry_v2, see ecalacc_geometry_v2_cnt.py).
import json
import math
from math import erf, sqrt
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

O = "/hdd3/przygoda/HADES/FAT_out/all_report/ecalacc/"
FIG = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_all/"
plt.rcParams.update({"font.size": 13, "axes.labelsize": 14, "axes.titlesize": 13.5, "legend.fontsize": 10.5, "xtick.labelsize": 12,
                     "ytick.labelsize": 12, "pdf.fonttype": 42})
PAL = dict(data="#000000", sig="#1f3fbf", cb="#e8202a", c2="#ff8c00", c3="#d400d4", c4="#009900", c5="#00a0a0", oeta="#7f3fbf",
           bkg="#888888", fit="#009900", dark="#444444")
NUM = json.load(open(O + "ecalacc_numbers.json"))


def save(fig, name):
    if name == "ecalacc_geometry":
        plt.close(fig); return
    fig.savefig(FIG + name + "_cnt.pdf"); fig.savefig(FIG + name + "_cnt.png", dpi=130); plt.close(fig)


def hread(fn):
    H = {}; L = open(fn).read().split("\n"); i = 0
    while i < len(L):
        if L[i].startswith("HIST"):
            _, nm, n, lo, hi = L[i].split(); n = int(n)
            v = np.array([list(map(float, L[i + 1 + k].split())) for k in range(n)])
            H[nm] = (float(lo), float(hi), v[:, 0], v[:, 1]); i += n + 1
        else:
            i += 1
    return H


HD = hread(O + "hist_data.txt"); HS = hread(O + "hist_sim.txt")


def centres(lo, hi, n):
    return lo + (np.arange(n) + 0.5) * (hi - lo) / n


# records
lay = open(O + "rows_c23_layout.txt").read().split("\n")
n, nd = int(lay[0].split()[1]), int(lay[0].split()[3])
names = lay[1].split()[1:]
A = np.fromfile(O + "rows_c23.bin", dtype=np.float64).reshape(n, nd)
R = {k: A[:, i] for i, k in enumerate(names)}

# =================================================================================================================
# 1. geometry: module map, photon phi and theta in data and simulation
# =================================================================================================================
cells = json.load(open(O + "geo_cells.json"))
fig = plt.figure(figsize=(17, 5.6))
ax = fig.add_axes([0.05, 0.13, 0.40, 0.78])
for c in cells:
    th, ph = np.array(c["th"]), np.array(c["ph"])
    # unwrap phi around the sector centre
    pc = (90.0 + 60.0 * c["sec"]) % 360.0
    ph = pc + ((ph - pc + 180.0) % 360.0 - 180.0)
    cen_t, cen_p = th.mean(), ph.mean()
    o = np.argsort(np.arctan2(th - cen_t, ph - cen_p))
    col = "#9ecae1" if c["status"] == 1 else ("#bbbbbb" if c["sec"] == 0 else "#e8202a")
    ax.add_patch(Polygon(np.c_[ph[o], th[o]], closed=True, fc=col, ec="#3a6fa0" if c["status"] == 1 else "#777777", lw=0.3))
ax.axhspan(18, 85, color="#cccccc", alpha=0.25, lw=0)
ax.axhline(18, color=PAL["dark"], ls=":", lw=1)
ax.text(358, 18.6, "18$^\\circ$: start of the tracking acceptance", fontsize=9.5, color=PAL["dark"], ha="right")
ax.set_xlim(0, 360); ax.set_ylim(3, 57)
ax.set_xlabel("$\\phi$ [deg] (seen from the target centre)"); ax.set_ylabel("$\\theta$ [deg]")
ax.set_title("(a) ECAL modules in the simulation and in the data parameters")
from matplotlib.patches import Patch
ax.legend(handles=[Patch(fc="#9ecae1", ec="#3a6fa0", label="active (807)"), Patch(fc="#bbbbbb", label="sector 0: off (163)"),
                   Patch(fc="#e8202a", label="dead single modules (8)")], loc="upper center", fontsize=10, ncol=3, frameon=False)
for s in range(6):
    ax.text((90 + 60 * s) % 360, 6.0, "sector %d" % s, ha="center", fontsize=10)

# (b) phi of the tagged photon, data vs simulation (shapes)
ax = fig.add_axes([0.51, 0.13, 0.22, 0.78])
lo, hi, d, dv = HD["data_us_gphi"]; _, _, s_, sv = HS["sim_T12_gphi"]
x = centres(lo, hi, len(d))
ax.step(x, d / d.sum(), where="mid", color=PAL["data"], lw=1.6, label="data, tagged counts")
ax.step(x, s_ / s_.sum(), where="mid", color=PAL["c2"], lw=1.4, label="simulation T1+T2, tagged")
ax.set_xlabel("photon $\\phi$ (cluster) [deg]"); ax.set_ylabel("share per 5$^\\circ$")
ax.set_title("(b) azimuth of the ECAL photon")
ax.legend(loc="upper center", fontsize=9.5)
ax.set_ylim(0, 1.45 * max(d.max() / d.sum(), s_.max() / s_.sum()))

# (c) theta: tagged photons (cluster) in data and simulation, and the TRUE photon of all accepted bump rows
ax = fig.add_axes([0.78, 0.13, 0.20, 0.78])
lo, hi, d, dv = HD["data_us_gtheta"]; _, _, s_, sv = HS["sim_T12_gtheta"]
x = centres(lo, hi, len(d))
th = NUM["thspec"]; tha = np.array(th["all"]); xe = np.array(th["edges"])
ax.step(0.5 * (xe[1:] + xe[:-1]), tha / tha.sum(), where="mid", color=PAL["c3"], lw=1.8, label="true photon of all\naccepted bump counts (C2+C3)")
ax.step(x, d / d.sum() * (tha[(xe[:-1] >= 11) & (xe[:-1] < 48)].sum() / tha.sum()), where="mid", color=PAL["data"], lw=1.4, label="data, tagged (scaled)")
ax.step(x, s_ / s_.sum() * (tha[(xe[:-1] >= 11) & (xe[:-1] < 48)].sum() / tha.sum()), where="mid", color=PAL["c2"], lw=1.2, ls="--", label="sim., tagged (scaled)")
ax.axvspan(11.5, 47.5, color="#9ecae1", alpha=0.3, lw=0)
ax.set_xlabel("photon $\\theta$ [deg]"); ax.set_ylabel("share per 2$^\\circ$")
sh = NUM["theta_shares"]
ax.set_title("(c) polar angle of the photon")
ax.set_ylim(0, 0.075)
ax.text(52, 0.043, "true photon:\n$<11.5^\\circ$: %.0f %%\n11.5-47.5$^\\circ$: %.0f %%\n$>47.5^\\circ$: %.0f %%" % (
    100 * sh["below11p5"], 100 * sh["band"], 100 * sh["above47p5"]), fontsize=10)
ax.legend(loc="upper right", fontsize=8.5)
ax.set_xlim(0, 90)
save(fig, "ecalacc_geometry")

# =================================================================================================================
# 2. fractions versus M(pi+pi-) and versus the distance to the edge
# =================================================================================================================
fig, axs = plt.subplots(1, 2, figsize=(15.5, 5.8))
curves = [("geo", "photon in an active ECAL module", "#3a6fa0", "o", "-"), ("geo_6sec", "the same, six full sectors", "#9ecae1", "s", ":"),
          ("geo_E100", "in a module and $E_\\gamma^{true}>100$ MeV", PAL["c3"], "^", "-"),
          ("found_quality", "found by the ECAL simulation (quality cuts)", PAL["dark"], "D", "-"),
          ("tag_true_m5c", "fully reconstructed (true photon, $m_{5c}$ window)", PAL["c2"], "o", "-"),
          ("tagTrueL3_m5c", "the same, cluster energies $\\times$1.285 (data-like)", PAL["c2"], "o", "--")]
for k, (key, xl, title) in enumerate([("vs_mpipi", "$M(\\pi^+\\pi^-)$ [GeV/$c^2$]", "(a) versus $M(\\pi^+\\pi^-)$"),
                                      ("vs_dist", "distance to the edge $M_4' - E(M(\\pi^+\\pi^-), M(e^+e^-))$ [MeV/$c^2$]",
                                       "(b) versus the distance to the bump edge")]):
    ax = axs[k]
    tab = NUM[key]
    lo = np.array([r["lo"] for r in tab]); hi = np.array([r["hi"] for r in tab])
    if key == "vs_mpipi":
        lo = np.where(lo == 0, 0.28, lo); xs = 0.5 * (lo + hi); xerr = 0.5 * (hi - lo)
    else:
        xs = 500 * (lo + hi); xerr = 500 * (hi - lo)
    for i, (q, lab, col, mk, ls) in enumerate(curves):
        y = np.array([r["C23"][q][0] for r in tab]); e = np.array([r["C23"][q][1] for r in tab])
        off = (i - 2.5) * (0.0012 if key == "vs_mpipi" else 0.8)
        ax.errorbar(xs + off, 100 * y, yerr=100 * e, xerr=xerr if i == 0 else None, fmt=mk, ls=ls, color=col, ms=5, lw=1.3, mfc="white" if ls == "--" else col,
                    label=lab, capsize=0)
    ax.set_xlabel(xl); ax.set_ylabel("share of the accepted bump counts C2+C3 [%]")
    ax.set_title(title); ax.set_ylim(0, 64)
    if key == "vs_dist":
        ax.axvline(0, color=PAL["dark"], ls=":", lw=1)
        ax.text(-205, 1.5, "far below the edge: hard photon", fontsize=10); ax.text(-8, 1.5, "at the edge: soft photon", fontsize=10, ha="right")
    else:
        ax.axvline(0.4116, color=PAL["dark"], ls=":", lw=1)
        ax.text(0.409, 1.5, "kinematic limit", fontsize=9, ha="right")
axs[0].legend(loc="upper center", fontsize=9.5, ncol=2)
axs[1].legend(loc="upper center", fontsize=9.5, ncol=2)
fig.tight_layout()
save(fig, "ecalacc_fractions")

# =================================================================================================================
# 3. photon energy: spectra and efficiency inside the modules
# =================================================================================================================
es = NUM["espec"]
ed = np.array(es["edges"]); ec = 0.5 * (ed[1:] + ed[:-1])
fig, axs = plt.subplots(1, 2, figsize=(15, 5.6))
ax = axs[0]
ax.step(ec, es["all"], where="mid", color=PAL["c3"], lw=1.8, label="all accepted bump counts (true $E_\\gamma$)")
ax.step(ec, es["acc"], where="mid", color="#3a6fa0", lw=1.8, label="photon in an active module")
ax.errorbar(ec, es["found"], yerr=np.sqrt(es["found_var"]), fmt="s", color=PAL["dark"], ms=4, label="found by the ECAL simulation (any cluster)")
ax.errorbar(ec, es["foundQ"], yerr=np.sqrt(es["foundQ_var"]), fmt="o", color=PAL["c2"], ms=4, label="found and quality cuts ($E_{cl}>100$ MeV, $\\beta$, PID)")
ax.axvline(100, color=PAL["dark"], ls=":", lw=1)
ax.set_xlabel("true photon energy $E_\\gamma$ [MeV]"); ax.set_ylabel("expected raw counts (SMASH scale) per 20 MeV")
ax.set_title("(a) photon energy of the bump counts C2+C3 (\"pre\")"); ax.legend(fontsize=9.5); ax.set_xlim(0, 600)
ax = axs[1]
# coarser bins (20 MeV units): 0-20, ..., 100-120, 120-160, 160-200, 200-260, 260-340, 340-600
grp = [[0], [1], [2], [3], [4], [5], [6, 7], [8, 9], [10, 11, 12], [13, 14, 15, 16], list(range(17, 30))]
acc = np.array(es["acc"]); fo = np.array(es["found"]); fov = np.array(es["found_var"]); fq = np.array(es["foundQ"]); fqv = np.array(es["foundQ_var"])
xc_, xe_, eff, ee, effQ, eeQ = [], [], [], [], [], []
for gg in grp:
    a_ = acc[gg].sum()
    if a_ <= 0:
        continue
    lo_, hi_ = ed[gg[0]], ed[gg[-1] + 1]
    xc_.append(0.5 * (lo_ + hi_)); xe_.append(0.5 * (hi_ - lo_))
    eff.append(fo[gg].sum() / a_); ee.append(math.sqrt(fov[gg].sum()) / a_)
    effQ.append(fq[gg].sum() / a_); eeQ.append(math.sqrt(fqv[gg].sum()) / a_)
ax.errorbar(xc_, eff, yerr=ee, xerr=xe_, fmt="s", color=PAL["dark"], ms=5, label="any cluster")
ax.errorbar(np.array(xc_) + 4, effQ, yerr=eeQ, xerr=xe_, fmt="o", color=PAL["c2"], ms=5, label="cluster passing the quality cuts")
ax.axvline(100, color=PAL["dark"], ls=":", lw=1); ax.axhline(1, color=PAL["dark"], lw=0.6)
cl = NUM["closure"]["ecl_over_etrue_q_quality"]
ax.text(230, 0.06, "simulated cluster energy / true energy:\nmedian %.3f (16-84 %%: %.2f-%.2f);\ndata: photon energy about 17 %% too high" % (cl[1], cl[0], cl[2]), fontsize=10)
ax.set_xlabel("true photon energy $E_\\gamma$ [MeV]"); ax.set_ylabel("efficiency inside the active modules")
ax.set_title("(b) ECAL simulation: efficiency inside the modules"); ax.set_ylim(0, 1.15); ax.set_xlim(0, 600); ax.legend(loc="upper right", fontsize=10)
fig.tight_layout()
save(fig, "ecalacc_energy")

# =================================================================================================================
# 4. the tag: data spectrum with the expectation, per-slice counts, and the data/expected ratios
# =================================================================================================================
spec = {}
for l in open(O + "tag_fits.txt"):
    p = l.split()
    if p and p[0] == "SPEC":
        spec["%s_%s" % (p[1], p[2])] = np.array(list(map(float, p[3:])))
    if p and p[0] == "FITPAR":
        spec["fit_%s_%s" % (p[1], p[2])] = p[3:]


def gfrac(a, b, m, s):
    return 0.5 * (erf((b - m) / (s * sqrt(2))) - erf((a - m) / (s * sqrt(2))))


fig, axs = plt.subplots(1, 3, figsize=(20, 5.8), gridspec_kw=dict(width_ratios=[1.15, 1.0, 1.0]))
ax = axs[0]
x = centres(0.40, 0.80, 160); bw = 0.0025
N = spec["0_data_N"]; CB = spec["0_data_CB"]; V = spec["0_data_sigVar"]
reb = 2
def rb(y):
    return y.reshape(-1, reb).sum(1)
xr = x.reshape(-1, reb).mean(1)
ax.errorbar(1e3 * xr, rb(N - CB), yerr=np.sqrt(rb(V)), fmt="o", color=PAL["sig"], ms=3.5, label="data $N-$CB (\"pre\")")
ax.step(1e3 * xr, rb(CB), where="mid", color=PAL["cb"], lw=1.2, label="CB")
fp = spec["fit_0_data"]; par = [float(v) for v in fp[:4]]; x0, dx = float(fp[5]), float(fp[7])
xf = np.linspace(0.49, 0.62, 300)
bk = sum(par[1 + j] * ((xf - x0) / dx) ** j for j in range(3))
tot = bk + par[0] * np.array([gfrac(v - bw / 2, v + bw / 2, 0.5463, 0.0062) for v in xf])
ax.plot(1e3 * xf, reb * tot, color=PAL["fit"], lw=2, label="frozen fit: $T$ = %.0f" % par[0])
ax.plot(1e3 * xf, reb * bk, color=PAL["fit"], lw=1.2, ls="--")
f_as = 0.3916
c2 = rb(spec["0_sim_C2"]); c3 = rb(spec["0_sim_C3"])
best = NUM["best_pre"]
lamF = NUM["lambda_fit"]["pre"]["L3_factor"]
ax.fill_between(1e3 * xr, 0, f_as * (c2 + c3), step="mid", color=PAL["c2"], alpha=0.35, label="SMASH C2+C3 tagged $\\times f_\\eta$=0.392")
ax.step(1e3 * xr, f_as * (c2 + c3) * lamF * NUM["resolution_factor"] * best["startf"], where="mid", color=PAL["c2"], lw=1.4, ls="--",
        label="the same, data-like ECAL calibration")
ax.set_xlim(450, 700); ax.set_xlabel("$m_{5c} = M(\\pi^+\\pi^-e^+e^-\\gamma)$, $\\gamma$ rescaled to $m_{\\pi^0}$ [MeV/$c^2$]")
ax.set_ylabel("counts per 5 MeV/$c^2$"); ax.set_title("(a) ECAL-tagged $\\eta$: data and expectation"); ax.legend(fontsize=9, loc="upper left"); ax.set_ylim(-100, 4200)
ax = axs[1]
tsn = ["<0.34", "0.34-0.37", "0.37-0.40", "0.40-0.43"]
xs = np.array([0.31, 0.355, 0.385, 0.415])
dT = np.array([NUM["tagfits"]["pre_" + t]["dT"] for t in tsn]); edT = np.array([NUM["tagfits"]["pre_" + t]["edT"] for t in tsn])
def exs(rd, L):
    e = [NUM["expected"]["pre_%s_%s_L%d" % (t, rd, L)] for t in tsn]
    v = np.array([x["value"] for x in e]); er = np.array([math.sqrt(x["stat"] ** 2 + x["syst"] ** 2 + x["mc"] ** 2) for x in e])
    return v, er
v0, e0 = exs("asis", 0); v3, e3 = exs("c2", 3)
v3 = v3 * NUM["resolution_factor"] * best["startf"]; e3 = e3 * NUM["resolution_factor"] * best["startf"]
ax.errorbar(xs, dT, yerr=edT, fmt="o", color=PAL["data"], ms=6, label="data $T$ (\"pre\")")
ax.errorbar(xs - 0.003, v0, yerr=e0, fmt="s", color=PAL["c2"], mfc="white", ms=6, label="expected: $f_\\eta$=0.392 $\\times$ simulated tag")
ax.errorbar(xs + 0.003, v3, yerr=e3, fmt="D", color=PAL["c2"], ms=6, label="expected: C2 corr., data-like calibration,\nresolution, START")
ax.set_xticks(xs); ax.set_xticklabels(["<0.34", "0.34-0.37", "0.37-0.40", "0.40-0.43"])
ax.set_xlabel("$M(\\pi^+\\pi^-)$ slice [GeV/$c^2$]"); ax.set_ylabel("tagged $\\eta$ (counts)"); ax.set_title("(b) per $M(\\pi^+\\pi^-)$ slice")
ax.legend(fontsize=9, loc="lower left"); ax.set_ylim(0, 2200)
ax = axs[2]
lab = []; val = []; err = []
for rd, L, t, extra in [("asis", 0, "as is, $f_\\eta$ = 0.392", 1.0), ("c2", 0, "C2 corrected, 0.418", 1.0), ("k48", 0, "conversions $\\times$0.48, 0.560", 1.0),
                        ("c2", 1, "C2 corr., $E\\times$1.10", 1.0), ("c2", 2, "C2 corr., $E\\times$1.20", 1.0), ("c2", 3, "C2 corr., $E\\times$1.285", 1.0),
                        ("c2", 3, "C2 corr., $E\\times$1.285, res., START", NUM["resolution_factor"] * best["startf"])]:
    e = NUM["expected"]["pre_all_%s_L%d" % (rd, L)]
    v = e["value"] * extra; et = math.sqrt(e["stat"] ** 2 + e["syst"] ** 2 + e["mc"] ** 2) * extra
    r = 3314.0 / v; er = r * math.sqrt((105 / 3314.) ** 2 + (et / v) ** 2)
    lab.append(t); val.append(r); err.append(er)
yy = np.arange(len(lab))[::-1]
ax.errorbar(val, yy, xerr=err, fmt="o", color=PAL["c2"], ms=6)
ax.axvline(1, color=PAL["dark"], lw=1)
ax.set_yticks(yy); ax.set_yticklabels(lab, fontsize=10)
ax.set_xlabel("data $T$ / expected $T$ (\"pre\", all $M(\\pi^+\\pi^-)$)"); ax.set_title("(c) data over expectation"); ax.set_xlim(0.4, 1.7)
fig.subplots_adjust(left=0.045, right=0.985, bottom=0.12, top=0.92, wspace=0.62)
save(fig, "ecalacc_tag")

# =================================================================================================================
# 5. closure and calibration
# =================================================================================================================
fig, axs = plt.subplots(1, 3, figsize=(17, 5.2))
ax = axs[0]
m = R["tMatch"] > 0
ax.hist(np.clip(1e4 * R["tCosRes"][m], -5, 5), bins=50, range=(-5, 5), color=PAL["sig"], histtype="step", lw=1.6)
ax.set_xlabel("$\\cos\\alpha_{true} - \\cos\\alpha_{circle}$ [$10^{-4}$]"); ax.set_ylabel("photons found by the ECAL simulation")
cl = NUM["closure"]
ax.set_title("(a) the found photons lie on the predicted circle")
ax.text(-4.8, 0.85 * ax.get_ylim()[1], "%d photons\ninside the predicted modules: %.1f %%\n$|\\Delta\\cos\\alpha| < 10^{-3}$: %.1f %%" % (
    cl["n_found"], 100 * cl["inside"], 100 * cl["frac_abs_cos_lt_1e3"]), fontsize=10)
ax = axs[1]
rat = (R["tEcl"] / np.maximum(R["tE"], 1e-9))[m]
ax.hist(rat, bins=60, range=(0, 1.5), color=PAL["c2"], histtype="step", lw=1.6, label="all found photons")
mq = m & (R["tPass"] > 0)
ax.hist((R["tEcl"] / np.maximum(R["tE"], 1e-9))[mq], bins=60, range=(0, 1.5), color=PAL["dark"], histtype="step", lw=1.4, label="passing the quality cuts")
ax.axvline(1, color=PAL["dark"], ls=":", lw=1); ax.axvline(1.17, color=PAL["data"], ls="--", lw=1.2)
ax.text(1.19, 0.8 * ax.get_ylim()[1], "data:\n+17 %", fontsize=10)
ax.set_xlabel("$E_{cluster} / E_{\\gamma}^{true}$ (simulation)"); ax.set_ylabel("photons"); ax.set_title("(b) energy scale of the simulated ECAL")
ax.legend(fontsize=9.5, loc="upper left")
ax = axs[2]
lo, hi, d, dv = HD["data_us_meeg"]; _, _, pp, _ = HD["data_pp_meeg"]; _, _, mm, _ = HD["data_mm_meeg"]
x = centres(lo, hi, len(d)); dd = d - 2 * np.sqrt(np.maximum(pp, 0) * np.maximum(mm, 0))
_, _, s0, _ = HS["sim_C23_meeg"]
g = (x > 0.06) & (x < 0.24)
ax.step(1e3 * x, dd / dd[g].max(), where="mid", color=PAL["sig"], lw=1.6, label="data: all e$^+$e$^-$ counts with a good photon ($N-$CB)")
x2 = x.reshape(-1, 2).mean(1); g2 = (x2 > 0.06) & (x2 < 0.24)
s2 = s0.reshape(-1, 2).sum(1)
ax.step(1e3 * x2, s2 / s2[g2].max(), where="mid", color=PAL["c2"], lw=1.6, label="simulation C2+C3, as simulated")
mm3 = (R["pre"] > 0) & (R["meegL3"] > 0) & (R["mee"] < 0.1)
h3, _ = np.histogram(R["meegL3"][mm3], bins=50, range=(0, 0.30), weights=R["w"][mm3])
ax.step(1e3 * x2, h3 / h3[g2].max(), where="mid", color=PAL["c2"], lw=1.4, ls="--", label="simulation C2+C3, energies $\\times$1.285")
ax.axvspan(119, 173, color="#cccccc", alpha=0.35, lw=0)
pz = NUM["pi0"]
ax.text(178, 0.80, "Gaussian fits:\ndata %.1f / %.1f MeV\nsim. %.1f / %.1f MeV\nsim. $\\times$1.285: %.1f / %.1f MeV" % (
    1e3 * pz["data"]["mean"], 1e3 * pz["data"]["sigma"], 1e3 * pz["sim_L0"]["mean"], 1e3 * pz["sim_L0"]["sigma"],
    1e3 * pz["simrec_L3"]["mean"], 1e3 * pz["simrec_L3"]["sigma"]), fontsize=9.5)
ax.set_xlim(60, 260); ax.set_ylim(0, 1.5)
ax.set_xlabel("$M(e^+e^-\\gamma)$ [MeV/$c^2$]"); ax.set_ylabel("normalised to the maximum"); ax.set_title("(c) $\\pi^0$ in $M(e^+e^-\\gamma)$ (grey: tag window)")
ax.legend(fontsize=8.5, loc="upper left")
fig.tight_layout()
save(fig, "ecalacc_closure")
print("figures written")
