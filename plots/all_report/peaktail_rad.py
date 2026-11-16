#!/usr/bin/env python3
# ============================================================================
# peaktail_rad.py  --  package peaktail: the radiative tail that the simulation does not contain (photon emission in the
#                      decay eta -> pi+pi-e+e-, "internal bremsstrahlung"), a leading-logarithm toy on the simulated C5 rows;
#                      and the dependence of the peak side bands (Y_peak, a_P) and of the X17 counting window on the tail
# ============================================================================
# Written 2026-11-15 (prefix peaktail).
# Physics (standard QED, final-state radiation off a lepton pair of mass q = M(ee), collinear / leading-log approximation):
#   each lepton emits a photon carrying the fraction z of its energy with
#       dP/dz = (c/2) (1 + (1 - z)^2) / z,   c = (alpha/pi) [ (1 + b^2)/(2 b) ln((1 + b)/(1 - b)) - 1 ],  b = sqrt(1 - 4 m_e^2/q^2)
#   (b -> 1: c = (alpha/pi) [ln(q^2/m_e^2) - 1]; e.g. q = 25 MeV: c = 0.0158, q = 5 MeV: c = 0.0083). Leading-log accuracy,
#   about 20-30 % (non-logarithmic terms, radiation of the pions and its interference, multiple emission are left out).
#   The generator (PLUTO, "eta --> e+ + e- + pi+ pi-", no PHOTOS) and HGeant write no such photon: the generated four-body
#   mass of every simulated C5 row is 548.000 MeV/c2 (peaktail_sim.py). Class C4 (a real photon that converts) has no such
#   term; its radiation in the material is simulated by HGeant.
# The toy: every simulated C5 row (S0, mode A, w x V3) is copied NREP times; in each copy each lepton radiates with the
# probability P(z > z0) and z drawn from the spectrum; the reconstructed lepton momentum is scaled by (1 - z) (the photon
# is emitted at the vertex, before the detector); M4 and M(ee) are recomputed from the four-vectors; then the M4 transfer
# and u as in edgeu. Outputs: /hdd3/.../peaktail/peaktail_rad.json, peaktail_rad.log, peaktail_rad_rows.npz
# Run: python3 -I plots/all_report/peaktail_rad.py
# ============================================================================
import json
import math
import sys

import numpy as np

sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import edgeu_core as C
import peaktail_core as P

LOG = open(P.OUT + "peaktail_rad.log", "w")
def L(*a):
    s = " ".join(str(x) for x in a)
    print(s); LOG.write(s + "\n"); LOG.flush()
RES = {}
ALPHA = 1 / 137.035999
ME = 0.00051099895
Z0 = 0.002
NREP = 3
S = P.load_sim(); w = P.sim_weight(S); gS = S["g"]
us0, mpp0, m4p0 = C.u_sim(S)
M4r, J = P.m4_and_jacobian(S)
dl = C.S_NOM * (J[0] * (S["ep_p"] - S["ep_ptrue"]) + J[1] * (S["em_p"] - S["em_ptrue"])).astype(np.float64) * 1e-3
flag = (dl < -0.005)
c5 = np.where((gS == 3) & (w > 0))[0]
L("# peaktail_rad.py: %d simulated C5 rows (S0, mode A); leading-log final-state radiation, %d copies per row, z0 = %.3f" % (len(c5), NREP, Z0))
q = np.maximum(S["meet"][c5].astype(np.float64), 2.0001 * ME)
b = np.sqrt(1 - 4 * ME * ME / (q * q))
cpair = ALPHA / math.pi * ((1 + b * b) / (2 * b) * np.log((1 + b) / (1 - b)) - 1)
I0 = -1.5 - 2 * math.log(Z0) + 2 * Z0 - Z0 * Z0 / 2
Plep = 0.5 * cpair * I0
L("   pair coefficient c: weighted median %.4f; per-lepton emission probability (z > %.3f): weighted mean %.4f, E1 %.4f, E2 %.4f" %
  (P.wq(cpair, w[c5], [0.5])[0], Z0, np.average(Plep, weights=w[c5]), np.average(Plep[S["mee"][c5] < 0.010], weights=w[c5][S["mee"][c5] < 0.010]),
   np.average(Plep[S["mee"][c5] >= 0.010], weights=w[c5][S["mee"][c5] >= 0.010])))
# inverse CDF of z on [z0, 1)
zg = np.exp(np.linspace(math.log(Z0), 0, 4001))
dens = (1 + (1 - zg) ** 2) / zg
cdf = np.concatenate([[0], np.cumsum(0.5 * (dens[1:] + dens[:-1]) * np.diff(zg))]); cdf /= cdf[-1]
rng = np.random.default_rng(0xBADC0DE)
# reconstructed four-vectors of the C5 rows
parts = [("ep", ME), ("em", ME), ("pip", P.MPI), ("pim", P.MPI)]
V = {a: P.fourvec(S[a + "_p"][c5] * 1e-3, S[a + "_th"][c5], S[a + "_ph"][c5], m) for a, m in parts}
def mass(vs):
    E = sum(v[0] for v in vs); X = sum(v[1] for v in vs); Y = sum(v[2] for v in vs); Z = sum(v[3] for v in vs)
    return np.sqrt(np.maximum(E * E - X * X - Y * Y - Z * Z, 0))
def scaled(v, f, m):
    px, py, pz = v[1] * f, v[2] * f, v[3] * f
    return np.sqrt(px * px + py * py + pz * pz + m * m), px, py, pz
M4base = mass([V[a] for a, _ in parts]); Meebase = mass([V["ep"], V["em"]])
rows = []; wts = []; dM4 = []; m4new = []; meenew = []; zmax = []
for r in range(NREP):
    fz = []
    for a in ("ep", "em"):
        em = rng.random(len(c5)) < Plep
        z = np.where(em, np.interp(rng.random(len(c5)), cdf, zg), 0.0)
        fz.append(z)
    vp = scaled(V["ep"], 1 - fz[0], ME); vm = scaled(V["em"], 1 - fz[1], ME)
    M4n = mass([vp, vm, V["pip"], V["pim"]]); Meen = mass([vp, vm])
    rows.append(c5); wts.append(w[c5] / NREP); dM4.append(M4n - M4base)
    m4new.append(S["m4"][c5].astype(np.float64) + (M4n - M4base)); meenew.append(S["mee"][c5].astype(np.float64) * np.where(Meebase > 0, Meen / Meebase, 1.0))
    zmax.append(np.maximum(fz[0], fz[1]))
rows = np.concatenate(rows); wts = np.concatenate(wts); dM4 = np.concatenate(dM4); m4new = np.concatenate(m4new); meenew = np.concatenate(meenew)
zmax = np.concatenate(zmax)
# transferred M4', u
mt = S["m4t"][rows].astype(np.float64)
m4p_new = C.S_NOM * (mt + C.R_NOM * (m4new - mt))
mpp_r = mpp0[rows]
u_new = m4p_new - C.edge(mpp_r, meenew)
u_new[~(meenew < C.MPI0)] = np.nan
emitted = zmax > 0
L("   rows with an emission: %.4f (weighted); with M4 shift < -5 / -15 / -25 MeV from the radiation: %.4f / %.4f / %.4f" %
  (P.frac(emitted, wts), P.frac(dM4 < -0.005, wts), P.frac(dM4 < -0.015, wts), P.frac(dM4 < -0.025, wts)))
# the residual d and the tail share, with and without the radiation
d0 = P.resid(S)[rows]; d1 = d0 + C.S_NOM * dM4
fl0 = flag[rows]; fl1 = (dl[rows] + C.S_NOM * dM4) < -0.005
out = {}
for rn, (a, bb) in (("all", (0.26, 0.62)), ("fit 0.30-0.40", (0.30, 0.40)), ("0.415-0.46", (0.415, 0.46)), ("0.42-0.55", (0.42, 0.55))):
    m = (mpp_r >= a) & (mpp_r < bb)
    for cn, cm in (("all", np.ones(len(rows), bool)), ("E1", S["mee"][rows] < 0.010), ("E2", S["mee"][rows] >= 0.010)):
        mm = m & cm
        o = dict(tail0=P.frac(fl0[mm], wts[mm]), tail1=P.frac(fl1[mm], wts[mm]),
                 d15_0=P.frac(d0[mm] < -0.015, wts[mm]), d15_1=P.frac(d1[mm] < -0.015, wts[mm]),
                 d25_0=P.frac(d0[mm] < -0.025, wts[mm]), d25_1=P.frac(d1[mm] < -0.025, wts[mm]))
        o["t_eff"] = o["tail1"] / o["tail0"]
        out["%s|%s" % (rn, cn)] = o
        L("   C5 %-14s %-3s tail share (lepton shift < -5 MeV) SMASH %.4f -> with radiation %.4f (t_C5 = %.2f); d < -15 MeV %.4f -> %.4f; d < -25 MeV %.4f -> %.4f" %
          (rn, cn, o["tail0"], o["tail1"], o["t_eff"], o["d15_0"], o["d15_1"], o["d25_0"], o["d25_1"]))
RES["toy"] = out
RES["toy_setup"] = dict(z0=Z0, nrep=NREP, plep_mean=float(np.average(Plep, weights=w[c5])), cpair_median=P.wq(cpair, w[c5], [0.5])[0])
np.savez_compressed(P.OUT + "peaktail_rad_rows.npz", rows=rows, w=wts, m4p=m4p_new, u=u_new, mee=meenew, mpp=mpp_r, d1=d1, fl1=fl1)

# ------------------------------------------------------------------------------------------------ peak side bands, a_P, X17 window
L("\nPEAK SIDE BANDS (WP6 estimator Y = [(N - CB)_P - R' (N - CB)_SB] / (1 - R' f), f = (C4 + C5)(SB)/(C4 + C5)(P), R' = 1.1859)")
PW = (0.53364, 0.55416); SBW = ((0.5250, 0.53364), (0.55416, 0.5628))
Rp = 1.1859
pk = (gS == 2) | (gS == 3)
def windows(wt, x, sel):
    inP = (x >= PW[0]) & (x < PW[1]); inSB = ((x >= SBW[0][0]) & (x < SBW[0][1])) | ((x >= SBW[1][0]) & (x < SBW[1][1]))
    return float(wt[sel & inP].sum()), float(wt[sel & inSB].sum()), float(wt[sel].sum())
def wt_t(t):
    """tail scaled by t, rows kept, per class over all M(pipi) of the dump"""
    out = w.copy()
    for gg in (2, 3):
        m = (gS == gg) & (w > 0)
        N = w[m].sum(); Nt = (w * flag)[m].sum(); Nc = N - Nt
        cc = (N - t * Nt) / Nc
        out[m] = w[m] * ((~flag[m]) * cc + t * flag[m])
    return out
sb = {}
base = pk & (w > 0)
P1, SB1, N1 = windows(w, m4p0, base)
f1 = SB1 / P1
L("   SMASH (t = 1): P %.1f, SB %.1f, all M4 %.1f expected raw rows; f = %.4f [WP6 0.1544]; P / all = %.4f" % (P1, SB1, N1, f1, P1 / N1))
sb["t1"] = dict(P=P1, SB=SB1, N=N1, f=f1)
for t in (0.0, 0.5, 1.5, 2.0):
    wt = wt_t(t)
    Pt, SBt, Nt = windows(wt, m4p0, base)
    ft = SBt / Pt
    yr = (1 - Rp * f1) / (1 - Rp * ft)          # Y_peak(t) / Y_peak(1) for the same data
    aPr = yr * P1 / Pt                            # a_P(t) / a_P(1)
    sigr = yr * (P1 / N1) / (Pt / Nt)             # sigma_peak(t)/sigma_peak(1): Y / (P fraction) relative
    sb["t%.1f" % t] = dict(P=Pt, SB=SBt, N=Nt, f=ft, Y_rel=yr, aP_rel=aPr, sigma_rel=sigr)
    L("   t = %.1f: f = %.4f, P/all %.4f; Y_peak x %.4f, a_P x %.4f, sigma(peak) x %.4f" % (t, ft, Pt / Nt, yr, aPr, sigr))
# with the radiation toy (C5 replaced)
wr_all = w.copy(); wr_all[gS == 3] = 0.0
Pa, SBa, Na = windows(wr_all, m4p0, base & (gS == 2))
Pr, SBr, Nr = windows(wts, m4p_new, np.ones(len(rows), bool))
Pt, SBt, Nt = Pa + Pr, SBa + SBr, Na + Nr
ft = SBt / Pt; yr = (1 - Rp * f1) / (1 - Rp * ft); aPr = yr * P1 / Pt; sigr = yr * (P1 / N1) / (Pt / Nt)
sb["rad"] = dict(P=Pt, SB=SBt, N=Nt, f=ft, Y_rel=yr, aP_rel=aPr, sigma_rel=sigr)
L("   C5 with the radiation toy: f = %.4f, P/all %.4f; Y_peak x %.4f, a_P x %.4f, sigma(peak) x %.4f" % (ft, Pt / Nt, yr, aPr, sigr))
RES["sidebands"] = sb
# X17 counting window X17P: P x M(ee) 15.5-24.5 MeV (C5 rows at that M(ee) as the proxy of the X17 signal shape)
L("\nX17 COUNTING WINDOW (P x M(ee) 15.5-24.5 MeV): content of the simulated C5 rows of that M(ee) (proxy of the X17 signal, M(ee) = 17 MeV)")
mx = (gS == 3) & (w > 0) & (S["mee"] >= 0.0155) & (S["mee"] < 0.0245)
inPw = (m4p0 >= PW[0]) & (m4p0 < PW[1])
e1 = P.frac(inPw[mx], w[mx])
mxr = (meenew >= 0.0155) & (meenew < 0.0245)
inPr = (m4p_new >= PW[0]) & (m4p_new < PW[1])
e_r = float(wts[mxr & inPr].sum() / w[mx].sum())
x17 = dict(eff_P_smash=e1, eff_P_rad=e_r, rel=e_r / e1)
for t in (0.0, 2.0):
    wt = wt_t(t)
    x17["eff_P_t%.0f" % t] = float(wt[mx & inPw].sum() / wt[mx].sum())
L("   fraction in P: SMASH %.4f; with the radiation toy %.4f (x %.4f); tail scale t = 0: %.4f, t = 2: %.4f" % (e1, e_r, e_r / e1, x17["eff_P_t0"], x17["eff_P_t2"]))
RES["x17"] = x17
json.dump(RES, open(P.OUT + "peaktail_rad.json", "w"), indent=1)
L("\nwritten: peaktail_rad.json, peaktail_rad.log, peaktail_rad_rows.npz")
