#!/usr/bin/env python3
# vecalacc_figs.py -- INDEPENDENT VERIFIER of package ecalacc (2026-11-15): verification figure
#   docs/latex/figs_all/vecalacc_theta.{pdf,png}
#   (a) polar angle of the lost photon of the accepted bump rows C2+C3 ("pre", own completion, own SMASH sample),
#       over the FULL range 0-180 deg, all photons and those in an active module; ECAL corner range from z = 0 and from
#       the target centre z = -124 mm; shares written for the whole sphere;
#   (b) covered azimuth fraction versus polar angle from z = 0 and z = -124 mm (own geometry);
#   (c) empirical check: true polar angle of target photons that made an ECAL cluster (all standard-set events).
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import vecalacc_geo as G

O = "/hdd3/przygoda/HADES/FAT_out/all_report/vecalacc"
F = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_all/vecalacc_theta"
d = np.load(O + "/vfrac_rows.npz")
thH, thHacc = d["thH"], d["thHacc"]
W = thH.sum()
edges = np.arange(0, 182, 2.0)
plt.rcParams.update({"font.size": 12})
fig, ax = plt.subplots(1, 3, figsize=(18, 5.2))
a = ax[0]
a.stairs(thH / W, edges, color="m", lw=2, label="lost photon, all accepted bump rows C2+C3")
a.stairs(thHacc / W, edges, color="tab:blue", lw=2, label="the same, photon in an active ECAL module")
a.axvspan(11.44, 47.56, color="tab:blue", alpha=0.12, label="ECAL corners seen from z = 0: 11.4-47.6 deg")
a.axvspan(10.90, 45.59, color="none", hatch="//", ec="grey", lw=0, label="seen from the target centre z = -124 mm: 10.9-45.6 deg")
c = 0.5 * (edges[1:] + edges[:-1])
wr = d["w"]; sh = lambda k: (wr * d[k]).sum() / wr.sum()          # exact per-row shares (not binned)
txt = (f"shares of ALL lost photons (whole sphere):\n< 11.5 deg: {100*sh('fLow'):.1f} %\n11.5-47.5 deg: {100*sh('fBand'):.1f} %\n"
       f"> 47.5 deg: {100*sh('fHigh'):.1f} %  (> 90 deg: {100*sh('fBack'):.1f} %)\nin an active module: {100*sh('fGeo'):.1f} %")
a.text(0.36, 0.50, txt, transform=a.transAxes, fontsize=10.5, va="top")
a.set_xlabel(r"polar angle of the lost photon $\theta_\gamma$ [deg]"); a.set_ylabel("share per 2 deg (all photons = 1)")
a.set_title("(a) lost photon of the bump rows (simulation)"); a.set_xlim(0, 180); a.set_ylim(0, 0.05); a.legend(fontsize=9, loc="upper right")
# (b)
g = G.Geo(G.load_geo()[0])
nphi = 3600; phis = np.radians((np.arange(nphi) + 0.5) * 360 / nphi)
ths = np.arange(8.0, 50.01, 0.25)
for z0, col, ls in ((0.0, "tab:blue", "-"), (-124.0, "black", "--")):
    prof = []
    for th in ths:
        t = np.radians(th)
        dd = np.vstack([np.sin(t) * np.cos(phis), np.sin(t) * np.sin(phis), np.full(nphi, np.cos(t))]).T
        V = np.tile([0.0, 0.0, z0], (nphi, 1))
        prof.append((g.hit(V, dd)[3].mean(), g.hit(V, dd, six=True)[3].mean()))
    prof = np.array(prof)
    ax[1].plot(ths, prof[:, 0], color=col, ls=ls, lw=2, label=f"active modules, from z = {z0:.0f} mm")
    ax[1].plot(ths, prof[:, 1], color=col, ls=ls, lw=1, alpha=0.5, label=f"six full sectors, from z = {z0:.0f} mm")
ax[1].set_xlabel(r"polar angle $\theta$ [deg]"); ax[1].set_ylabel("covered share of the azimuth")
ax[1].set_title("(b) ECAL azimuthal coverage (own geometry dump)"); ax[1].legend(fontsize=9, loc="lower center"); ax[1].set_ylim(0, 1)
# (c)
a3 = np.loadtxt(O + "/vclusters_std.txt", comments="#")
sid = a3[:, 9].astype(int); p = a3[:, 11:14]; v = a3[:, 14:17]
ph = (sid == 1) & (np.hypot(v[:, 0], v[:, 1]) < 10) & (np.abs(v[:, 2] + 124) < 40)
th = np.degrees(np.arccos(p[ph, 2] / np.linalg.norm(p[ph], axis=1)))
ax[2].hist(th, bins=np.arange(0, 60.5, 0.5), histtype="step", color="tab:gray", lw=2, label=f"true direction, target photons ({ph.sum()})")
ax[2].axvline(10.90, color="black", ls="--", lw=1); ax[2].axvline(45.59, color="black", ls="--", lw=1, label="corner range from z = -124 mm")
ax[2].axvline(11.44, color="tab:blue", ls=":", lw=1); ax[2].axvline(47.56, color="tab:blue", ls=":", lw=1, label="corner range from z = 0")
ax[2].set_xlabel(r"true polar angle $\theta_\gamma$ [deg]"); ax[2].set_ylabel("photons per 0.5 deg"); ax[2].set_ylim(0, 1700)
ax[2].set_title("(c) sim.: photons with an ECAL cluster"); ax[2].legend(fontsize=9, loc="upper right")
fig.tight_layout()
fig.savefig(F + ".pdf"); fig.savefig(F + ".png", dpi=110)
print("written", F + ".pdf/.png")
