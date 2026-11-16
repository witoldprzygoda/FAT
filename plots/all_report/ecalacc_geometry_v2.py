#!/usr/bin/env python3
# ============================================================================
# ecalacc_geometry_v2.py -- Part III integration (2026-11-15): corrected ECAL geometry figure and table of package ecalacc
# ============================================================================
# Replaces (as NEW files; the originals are kept) figs_all/ecalacc_geometry.{pdf,png} and tables_all/ecalacc_geometry.tex
# with the corrections of the verification vecalacc (ecalacc_results.md, Sect. 7):
#   - the module map and the polar-angle coverage are drawn as seen from the TARGET CENTRE z = -124 mm (the original
#     figure said "target centre" but used z = 0); module corners 10.9-45.6 deg;
#   - panel (c): the polar angle of the lost photon over the WHOLE sphere (0-180 deg); shares 9.0 / 54.5 / 36.5 %
#     (9.9 % backward) instead of 11 / 60 / 29 % (which were normalised to theta < 90 deg);
#   - colours: the true photon is drawn in brown and the simulated tag in grey, so that they are not read as the Part III
#     classes C3 (magenta) and C2 (orange).
# Inputs (read only): /hdd3/przygoda/HADES/FAT_out/all_report/vecalacc/{vemcgeo_stdout.txt, vfrac_rows.npz} (verifier's
# own geometry dump and lost-photon histogram), /hdd3/przygoda/HADES/FAT_out/all_report/ecalacc/{hist_data.txt,
# hist_sim.txt} (tagged-photon cluster angles of data and simulation).
# Outputs: docs/latex/figs_all/ecalacc_geometry_v2.{pdf,png}, docs/latex/tables_all/ecalacc_geometry_v2.tex
# Run: python3 -I plots/all_report/ecalacc_geometry_v2.py
# ============================================================================
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Patch
sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import vecalacc_geo as G

OE = "/hdd3/przygoda/HADES/FAT_out/all_report/ecalacc/"
OV = "/hdd3/przygoda/HADES/FAT_out/all_report/vecalacc/"
FIG = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_all/ecalacc_geometry_v2"
TAB = "/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all/ecalacc_geometry_v2.tex"
ZT = -124.0                      # target centre [mm] (simulated vertex uniform in -148...-99 mm; data median -125 mm)
plt.rcParams.update({"font.size": 13, "axes.labelsize": 14, "axes.titlesize": 13.5, "legend.fontsize": 10.5,
                     "xtick.labelsize": 12, "ytick.labelsize": 12, "pdf.fonttype": 42})
BROWN, GREY, BLUE = "#8c510a", "#666666", "#3a6fa0"


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


def angles(pts, z0):
    d = pts - np.array([0.0, 0.0, z0])
    th = np.degrees(np.arccos(d[:, 2] / np.linalg.norm(d, axis=1)))
    ph = np.degrees(np.arctan2(d[:, 1], d[:, 0])) % 360.0
    return th, ph


cells, _ = G.load_geo(OV + "vemcgeo_stdout.txt")
HD = hread(OE + "hist_data.txt"); HS = hread(OE + "hist_sim.txt")
d = np.load(OV + "vfrac_rows.npz")
w = d["w"]; share = lambda k: (w * d[k]).sum() / w.sum()
sh = {k: 100 * share(k) for k in ("fLow", "fBand", "fHigh", "fBack", "fGeo")}

# ---------------------------------------------------------------- per-sector table (corners from z = 0 and z = -124 mm)
rows = []
for s in range(6):
    cs = [c for c in cells if c[0] == s]
    P = np.vstack([c[5] for c in cs])
    t0, p0 = angles(P, 0.0); tt, _ = angles(P, ZT)
    pc = (90.0 + 60.0 * s) % 360.0
    pu = pc + ((p0 - pc + 180.0) % 360.0 - 180.0)
    rows.append((s, len(cs), sum(1 for c in cs if c[3] == 1), t0.min(), t0.max(), tt.min(), tt.max(), pu.min() % 360, pu.max() % 360))
allc = np.vstack([c[5] for c in cells]); tta, _ = angles(allc, ZT); t0a, _ = angles(allc, 0.0)
with open(TAB, "w") as f:
    f.write("% ecalacc_geometry_v2.tex -- plots/all_report/ecalacc_geometry_v2.py (corrected ecalacc_geometry.tex: polar angles from the\n")
    f.write("% target centre z = -124 mm instead of z = -60 mm; solid angle from the target centre; vecalacc geometry dump)\n")
    f.write("\\begin{tabular}{@{}c c c c c c@{}}\n\\toprule\n")
    f.write("sector & modules & active & $\\theta$ from $z=0$ [deg] & $\\theta$ from the target centre, $z=-124$\\,mm [deg] & $\\phi$ from $z=0$ [deg] \\\\\n\\midrule\n")
    for r in rows:
        f.write("%d & %d & %d & %.1f--%.1f & %.1f--%.1f & %.1f--%.1f \\\\\n" % r)
    f.write("\\midrule\n")
    f.write("all & %d & %d & %.1f--%.1f & %.1f--%.1f & \\\\\n" % (len(cells), sum(1 for c in cells if c[3] == 1), t0a.min(), t0a.max(), tta.min(), tta.max()))
    f.write("\\multicolumn{6}{@{}l}{solid angle of the active modules from the target centre: 8.45\\,\\% of $4\\pi$ (six full sectors: 10.3\\,\\%);"
            " from $z=0$: 9.2\\,\\% (11.2\\,\\%)} \\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n")
print("table:", TAB)
for r in rows:
    print("sector %d  %d/%d  theta z=0 %.2f-%.2f  z=-124 %.2f-%.2f  phi %.1f-%.1f" % r)
print("all corners z=-124: %.2f-%.2f ; z=0: %.2f-%.2f" % (tta.min(), tta.max(), t0a.min(), t0a.max()))
print("shares:", sh)

# ---------------------------------------------------------------- figure
fig = plt.figure(figsize=(17, 5.6))
ax = fig.add_axes([0.05, 0.13, 0.40, 0.78])
for c in cells:
    s, st = c[0], c[3]
    th, ph = angles(c[5], ZT)
    pc = (90.0 + 60.0 * s) % 360.0
    ph = pc + ((ph - pc + 180.0) % 360.0 - 180.0)
    o = np.argsort(np.arctan2(th - th.mean(), ph - ph.mean()))
    col = "#9ecae1" if st == 1 else ("#bbbbbb" if s == 0 else "#e8202a")
    ax.add_patch(Polygon(np.c_[ph[o], th[o]], closed=True, fc=col, ec=BLUE if st == 1 else "#777777", lw=0.3))
ax.axhspan(18, 85, color="#cccccc", alpha=0.25, lw=0)
ax.axhline(18, color="#444444", ls=":", lw=1)
ax.text(358, 18.6, "18$^\\circ$: start of the tracking acceptance", fontsize=9.5, color="#444444", ha="right")
ax.set_xlim(0, 360); ax.set_ylim(3, 57)
ax.set_xlabel("$\\phi$ [deg]"); ax.set_ylabel("$\\theta$ [deg] from the target centre")
ax.set_title("(a) ECAL modules, simulation and data parameters")
ax.legend(handles=[Patch(fc="#9ecae1", ec=BLUE, label="active (807)"), Patch(fc="#bbbbbb", label="ECAL sector 0: off (163)"),
                   Patch(fc="#e8202a", label="dead single modules (8)")], loc="upper center", fontsize=10, ncol=3, frameon=False)
for s in range(6):
    ax.text((90 + 60 * s) % 360, 6.0, "sector %d" % s, ha="center", fontsize=10)

ax = fig.add_axes([0.51, 0.13, 0.20, 0.78])
lo, hi, dd, _ = HD["data_us_gphi"]; _, _, ss, _ = HS["sim_T12_gphi"]
x = lo + (np.arange(len(dd)) + 0.5) * (hi - lo) / len(dd)
ax.step(x, dd / dd.sum(), where="mid", color="black", lw=1.6, label="data, tagged rows")
ax.step(x, ss / ss.sum(), where="mid", color=GREY, lw=1.4, ls="--", label="simulation T1+T2, tagged")
ax.axvspan(62.8, 117.1, color="#bbbbbb", alpha=0.35, lw=0)
ax.text(90, 0.002, "ECAL\nsector 0", ha="center", fontsize=9)
ax.set_xlabel("photon $\\phi$ (cluster) [deg]"); ax.set_ylabel("share per 5$^\\circ$")
ax.set_title("(b) azimuth of the ECAL photon")
ax.legend(loc="upper center", fontsize=9.5)
ax.set_ylim(0, 1.45 * max(dd.max() / dd.sum(), ss.max() / ss.sum()))

ax = fig.add_axes([0.765, 0.13, 0.225, 0.78])
edges = np.arange(0, 182, 2.0); thH = d["thH"]; W = thH.sum()
cen = 0.5 * (edges[1:] + edges[:-1])
ax.stairs(thH / W, edges, color=BROWN, lw=1.8, label="lost photon, all accepted\nbump rows C2+C3 (sim.)")
inb = (edges[:-1] >= 12) & (edges[:-1] < 48)
scale = thH[inb].sum() / W
lo, hi, dd, _ = HD["data_us_gtheta"]; _, _, ss, _ = HS["sim_T12_gtheta"]
x = lo + (np.arange(len(dd)) + 0.5) * (hi - lo) / len(dd)
ax.step(x, dd / dd.sum() * scale, where="mid", color="black", lw=1.4, label="data, tagged (scaled)")
ax.step(x, ss / ss.sum() * scale, where="mid", color=GREY, lw=1.2, ls="--", label="sim., tagged (scaled)")
ax.axvspan(10.9, 45.6, color="#9ecae1", alpha=0.35, lw=0, label="ECAL from the target\ncentre: 10.9-45.6$^\\circ$")
ax.set_xlabel("photon $\\theta$ [deg]"); ax.set_ylabel("share per 2$^\\circ$ (whole sphere = 1)")
ax.set_title("(c) polar angle of the photon")
ax.set_ylim(0, 0.075); ax.set_xlim(0, 180)
ax.text(70, 0.034, "lost photon, whole sphere:\n$<11.5^\\circ$: %.1f %%\n11.5-47.5$^\\circ$: %.1f %%\n$>47.5^\\circ$: %.1f %%\n"
        "  (backward, $>90^\\circ$: %.1f %%)\nin an active module: %.1f %%" % (sh["fLow"], sh["fBand"], sh["fHigh"], sh["fBack"], sh["fGeo"]),
        fontsize=9.5, va="top")
ax.legend(loc="upper right", fontsize=8.3)
fig.savefig(FIG + ".pdf"); fig.savefig(FIG + ".png", dpi=130)
print("figure:", FIG + ".pdf/.png")
