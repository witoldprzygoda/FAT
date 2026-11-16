#!/usr/bin/env python3
"""vpt_qed.py -- verifier of peaktail: analytic leading-log estimate of the final-state radiation of eta -> pi+pi-e+e-
that the generator lacks (no toy): per C5 row and lepton, P(z > z0) = (c/2)[-2 ln z0 - 2(1 - z0) + (1 - z0^2)/2],
c = (alpha/pi)[(1 + b^2)/(2b) ln((1+b)/(1-b)) - 1] (b = lepton velocity in the e+e- frame), z0 = delta/E*_l (E*_l = lepton
energy in the eta frame; an M4 shift of delta needs a photon of about delta in the eta frame). New tail rows: rows that are
not tail rows (sub5) and radiate z > z0 with delta = 5 MeV. t_C5 = (N_tail + sum_core w P_row)/N_tail.
"""
import sys
import numpy as np
sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import vpt_core as V
import vpt_tmpl as TM
R = TM.rows(); S = R["S"]
ALPHA = 1 / 137.036
mass = {"ep": V.ME_E, "em": V.ME_E, "pip": V.MPI, "pim": V.MPI}
vt = {n: V.p4(S[n + "_pt"], S[n + "_tht"], S[n + "_ph"], mass[n]) for n in mass}   # true |p| and theta, reco phi
Et = sum(vt[n][0] for n in mass); Px = sum(vt[n][1] for n in mass); Py = sum(vt[n][2] for n in mass); Pz = sum(vt[n][3] for n in mass)
M = np.sqrt(np.maximum(Et ** 2 - Px ** 2 - Py ** 2 - Pz ** 2, 1e-12))
def estar(n):   # energy in the frame of the four-particle system
    return (Et * vt[n][0] - Px * vt[n][1] - Py * vt[n][2] - Pz * vt[n][3]) / M
q = np.maximum(S["eet"].astype(float), 2.0001 * V.ME_E)
b = np.sqrt(1 - 4 * V.ME_E ** 2 / q ** 2)
c = (ALPHA / np.pi) * ((1 + b * b) / (2 * b) * np.log((1 + b) / (1 - b)) - 1)
def P(z0, m=None):
    cc = c if m is None else c[m]
    z0 = np.clip(z0, 1e-6, 1.0)
    return np.where(z0 < 1, (cc / 2) * (-2 * np.log(z0) - 2 * (1 - z0) + 0.5 * (1 - z0 * z0)), 0.0)
for rn, a, bb in (("0.30-0.40", 0.30, 0.40), ("0.42-0.46", 0.42, 0.46), ("0.30-0.34", 0.30, 0.34), ("0.37-0.40", 0.37, 0.40)):
    m = R["m45"] & (S["g"] == 3) & TM.in_pp(R["mpp"], a, bb)
    w = R["w"][m]; fl = R["flag"]["sub5"][m]
    for delta in (0.005, 0.010):
        Prow = 1 - (1 - P(delta / estar("ep")[m], m)) * (1 - P(delta / estar("em")[m], m))
        Nt = w[fl].sum(); add = (w[~fl] * Prow[~fl]).sum()
        pl = (w * P(0.002 * np.ones(m.sum()), m)).sum() / w.sum()
        cls = {"E1": S["ee"][m] < 0.010, "E2": ~(S["ee"][m] < 0.010)}
        sub = " ".join("%s x%.3f" % (k, ((w[fl & cc]).sum() + (w[~fl & cc] * Prow[~fl & cc]).sum()) / w[fl & cc].sum()) for k, cc in cls.items())
        if delta == 0.005:
            print("C5 %s: median E*(e) %.0f MeV, median c %.4f; per-lepton P(z > 0.002) %.1f %%; delta 5 MeV: rows radiating %.2f %% -> t_C5 = x%.3f (%s)" % (
                rn, 1e3 * np.median(np.concatenate([estar("ep")[m], estar("em")[m]])), np.median(c[m]), 100 * pl, 100 * (w * Prow).sum() / w.sum(), (Nt + add) / Nt, sub))
        else:
            print("      delta 10 MeV (sub10 family): t_C5 = x%.3f" % (((R["w"][m][R["flag"]["sub10"][m]]).sum() + (w[~R["flag"]["sub10"][m]] * Prow[~R["flag"]["sub10"][m]]).sum()) / R["w"][m][R["flag"]["sub10"][m]].sum()))
# combined C4 + C5 in 0.30-0.40 (C4 unchanged)
m4 = R["m45"] & (S["g"] == 2) & TM.in_pp(R["mpp"], 0.30, 0.40); m5 = R["m45"] & (S["g"] == 3) & TM.in_pp(R["mpp"], 0.30, 0.40)
fl = R["flag"]["sub5"]
Prow = 1 - (1 - P(0.005 / estar("ep"))) * (1 - P(0.005 / estar("em")))
N4 = R["w"][m4 & fl].sum(); N5 = R["w"][m5 & fl].sum(); add = (R["w"][m5 & ~fl] * Prow[m5 & ~fl]).sum()
print("C4 + C5 0.30-0.40: tail rows C4 %.0f C5 %.0f + QED %.0f -> combined t = %.3f" % (N4, N5, add, (N4 + N5 + add) / (N4 + N5)))
