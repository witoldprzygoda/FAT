#!/usr/bin/env python3
# ============================================================================
# vedgeu_ana.py  --  INDEPENDENT VERIFIER of the package edgeu: kinematics, u spectra, global step fit, f, stat, MC,
#                    variants, null tests (own code; vedgeu_core.py + vedgeu_fill.C columns)
# ============================================================================
# Run (from /home/przygoda/HADES/FAT/FAT):  python3 -I plots/all_report/vedgeu_ana.py > .../vedgeu/vedgeu_ana_stdout.txt
# Output: /hdd3/przygoda/HADES/FAT_out/all_report/vedgeu/vedgeu_ana.json, vedgeu_spec.npz
# ============================================================================
import json
import math
import sys
import time

import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import vedgeu_core as V

T0 = time.time()
OUT = {}
D, S = V.load()
print("rows: data %d, sim %d (%.0f s)" % (len(D["k"]), len(S["c"]), time.time() - T0))

# =============================================================================== 1. kinematics
print("\n1. KINEMATICS")
pts = [(0.30, 0.0), (0.39, 0.0), (0.39, 0.060), (0.4116, 0.0), (0.30, 0.100), (0.30, 0.050), (0.35, 0.020)]
OUT["edge_points"] = {"%.4f,%.3f" % p: float(V.E(p[0], p[1])) for p in pts}
print("   own E(m, q):", OUT["edge_points"])
# independent brute-force check of the limit: random decays eta -> pi+ pi- pi0 at fixed M(pipi) = m, pi0 -> gamma (e+e-)(q),
# maximum of M(pi+pi-e+e-) over 2e5 random photon directions (pi0 frame) -> must approach M4max(m, q) from below
rng = np.random.default_rng(12345)
def brute(m, q, n=200000):
    E0 = (V.MET ** 2 + V.MP0 ** 2 - m ** 2) / (2 * V.MET); p0 = math.sqrt(E0 ** 2 - V.MP0 ** 2)
    bet, gam = p0 / E0, E0 / V.MP0
    eg = (V.MP0 ** 2 - q ** 2) / (2 * V.MP0)            # photon energy in the pi0 frame
    c = rng.uniform(-1, 1, n)                             # cos of the photon angle to the pi0 flight direction
    egl = gam * eg * (1 + bet * c)                        # photon energy in the eta frame
    m4 = np.sqrt(V.MET ** 2 - 2 * V.MET * egl)
    return float(m4.max())
OUT["edge_brute"] = {}
for m, q in [(0.30, 0.0), (0.35, 0.03), (0.39, 0.06), (0.40, 0.10)]:
    b = brute(m, q); f = float(V.m4max(m, q))
    OUT["edge_brute"]["%.2f,%.2f" % (m, q)] = [b, f]
    print("   brute force max M4 (m %.2f, q %.3f): %.6f; formula %.6f; diff %.2e" % (m, q, b, f, b - f))
# simulated truth: excess over the true edge
gC = S["c"]; s0s = (S["bits"] & 1) != 0
mt4 = V.m4max(S["ppt"].astype(float), S["eet"].astype(float))
mt40 = V.m4max(S["ppt"].astype(float), np.zeros(len(S["c"])))
OUT["truth"] = {}
for g, nm in [(0, "C2"), (1, "C3")]:
    m = (gC == g) & s0s & (S["m4t"] > 0) & (S["ppt"] > 0)
    d = (S["m4t"] - mt4)[m]; d0 = (S["m4t"] - mt40)[m]; ww = S["w"][m]
    # edge reached? the largest m4_true near given (m_true, q_true) cells
    OUT["truth"][nm] = dict(rows=int(m.sum()), max_excess_MeV=1e3 * float(np.nanmax(d)), n_above=int(np.sum(d > 0)),
                            q0_frac_above_1MeV=float(np.sum(ww[d0 > 0.001]) / ww.sum()), q0_frac_above_5MeV=float(np.sum(ww[d0 > 0.005]) / ww.sum()),
                            frac_within_2MeV=float(np.sum(ww[d > -0.002]) / ww.sum()))
    print("   %s: %d rows, max(m4_t - M4max(mpp_t, mee_t)) = %+.3f MeV, rows above %d; q = 0: %.2f %% above by > 1 MeV, %.2f %% > 5 MeV; within 2 MeV of the edge %.2f %%" %
          (nm, m.sum(), OUT["truth"][nm]["max_excess_MeV"], OUT["truth"][nm]["n_above"], 100 * OUT["truth"][nm]["q0_frac_above_1MeV"],
           100 * OUT["truth"][nm]["q0_frac_above_5MeV"], 100 * OUT["truth"][nm]["frac_within_2MeV"]))
# tightness of the limit at a few (M(pipi)_true, M(ee)_true) points: largest m4_true among C2 + C3 rows in a small cell
OUT["truth_cells"] = []
for (m, q, dm, dq) in [(0.31, 0.0, 0.01, 0.004), (0.35, 0.0, 0.01, 0.004), (0.38, 0.0, 0.01, 0.004), (0.32, 0.03, 0.01, 0.01), (0.36, 0.05, 0.01, 0.01), (0.33, 0.08, 0.02, 0.02)]:
    m_ = (gC <= 1) & (S["m4t"] > 0) & (np.abs(S["ppt"] - m) < dm / 2) & (np.abs(S["eet"] - q) < dq / 2 + (dq / 2 if q == 0 else 0))
    if m_.sum() == 0:
        continue
    mx = float(S["m4t"][m_].max())
    lim_lo = float(V.m4max(m - dm / 2, q + dq / 2)); lim_hi = float(V.m4max(m - dm / 2, max(q - dq / 2, 0) if q > 0 else 0.0))
    # largest possible M4 in the cell: smallest M(pipi), largest M(ee)
    lim = float(np.nanmax(mt4[m_]))
    OUT["truth_cells"].append(dict(m=m, q=q, rows=int(m_.sum()), max_m4t=mx, max_limit_rows=lim, edge_centre=float(V.m4max(m, q))))
    print("   cell M(pipi)_t %.2f+-%.3f, M(ee)_t %.3f: %d rows, max m4_t %.5f, max M4max of these rows %.5f, M4max(centre) %.5f" %
          (m, dm / 2, q, m_.sum(), mx, lim, float(V.m4max(m, q))))

# =============================================================================== 2. u, histograms
print("\n2. SPECTRA (%.0f s)" % (time.time() - T0))
ud = V.u_data(D)
bD = V.ubin(ud); sD = V.slice_index(D["pp"])
s0d = (D["bits"] & 1) != 0
def dhists(mask, W=None):
    out = []
    for k in range(3):
        mk = mask & (D["k"] == k)
        if W is None:
            out.append(V.hist(bD[mk], sD[mk]))
        else:
            out.append(mk)
    return out
Hus, Hpp, Hmm = dhists(s0d)
N, CB, CBC = V.spectrum(Hus, Hpp, Hmm)
print("   data S0 rows in 0.30 <= M(pipi) < 0.40 and the u grid: us %d pp %d mm %d" % (Hus.sum(), Hpp.sum(), Hmm.sum()))
us0, mpp0 = V.u_sim(S)
bS = V.ubin(us0); sS = V.slice_index(mpp0)
wS = np.where(s0s, S["w"], 0.0)
cfac = np.where(S["ee"] < 0.010, V.CF1, V.CF2)
def stmpl(g, w=wS, b=bS, s=sS, c2=True):
    m = (gC == g)
    ww = w * (cfac if (g == 0 and c2) else 1.0)
    return V.hist(b[m], s[m], ww[m]).sum(0)
TM = {nm: stmpl(g) for g, nm in enumerate(["C2c", "C3", "C4", "C5", "Geta", "GD", "GC"])}
TM["C2"] = stmpl(0, c2=False)
for nm in ["C2", "C2c", "C3", "C4", "C5", "Geta", "GD", "GC"]:
    print("   sim %-4s in the grid: %.1f expected raw rows; in -45..+45: %.1f" % (nm, TM[nm].sum(), TM[nm][(V.X > -0.045) & (V.X < 0.045)].sum()))
K = TM["C4"] + TM["C5"]
# spectra around the edge
sel = (V.X > -0.050) & (V.X < 0.060)
print("   u bins (MeV): N, CB, N - CB, C2c + C3, peak")
for i in np.where(sel)[0][::2]:
    print("     %+6.2f  %7.0f %6.1f %7.1f  %7.1f  %6.1f" % (1e3 * V.X[i], N[i], CB[i], N[i] - CB[i], TM["C2c"][i] + TM["C3"][i], K[i]))
# reconstructed sim C2/C3 above the edge
for g, nm in [(0, "C2"), (1, "C3")]:
    m = (gC == g) & s0s & np.isfinite(us0) & (mpp0 >= 0.30) & (mpp0 < 0.40)
    print("   reco %s: fraction u > +10 MeV %.4f, > +20 MeV %.4f" % (nm, wS[m & (us0 > 0.010)].sum() / wS[m].sum(), wS[m & (us0 > 0.020)].sum() / wS[m].sum()))

# =============================================================================== 3. resolution
print("\n3. RESOLUTION (%.0f s)" % (time.time() - T0))
ut = V.ESC * (S["m4t"].astype(float) - mt4)
mr = (gC <= 1) & s0s & np.isfinite(us0) & np.isfinite(ut) & (mpp0 >= 0.30) & (mpp0 < 0.40) & (ut > -0.015) & (S["m4t"] > 0)
dd, wd = (us0 - ut)[mr], S["w"][mr]
gau = lambda x, a, mu, s: a * np.exp(-0.5 * (x - mu) ** 2 / s ** 2)
def core(lo, hi, bw):
    e = np.arange(lo, hi + 1e-9, bw)
    h, _ = np.histogram(dd, e, weights=wd); h2, _ = np.histogram(dd, e, weights=wd ** 2)
    x = 0.5 * (e[1:] + e[:-1])
    p, cv = curve_fit(gau, x, h, p0=[h.max(), 0, 0.007], sigma=np.sqrt(np.maximum(h2, h2[h2 > 0].min())), absolute_sigma=True)
    return abs(p[2]), math.sqrt(cv[2, 2]), p[1]
s_a = core(-0.040, 0.040, 0.002)
# own iterated +-2 sigma core
sg = 0.007; mu = 0.0
for _ in range(10):
    e = np.arange(mu - 2 * sg, mu + 2 * sg + 1e-9, 0.001)
    h, _ = np.histogram(dd, e, weights=wd); h2, _ = np.histogram(dd, e, weights=wd ** 2)
    x = 0.5 * (e[1:] + e[:-1])
    p, cv = curve_fit(gau, x, h, p0=[h.max(), mu, sg], sigma=np.sqrt(np.maximum(h2, h2[h2 > 0].min())), absolute_sigma=True)
    mu, sg = p[1], abs(p[2])
s_b = (sg, math.sqrt(cv[2, 2]), mu)
q16, q50, q84 = np.percentile(dd, [15.865, 50, 84.135])   # unweighted percentiles (rows)
OUT["reso"] = dict(rows=int(mr.sum()), core_40=s_a[0], ecore_40=s_a[1], mean_40=s_a[2], core_2s=s_b[0], ecore_2s=s_b[1],
                   half68_rows=(q84 - q16) / 2, rms=float(math.sqrt(np.average((dd - np.average(dd, weights=wd)) ** 2, weights=wd))))
print("   %d rows; Gaussian core (+-40 MeV, 2 MeV bins) %.2f +- %.2f MeV (mean %+.2f); iterated +-2 sigma core %.2f +- %.2f MeV; half 68 %% (rows) %.2f; rms %.2f MeV" %
      (mr.sum(), 1e3 * s_a[0], 1e3 * s_a[1], 1e3 * s_a[2], 1e3 * s_b[0], 1e3 * s_b[1], 1e3 * OUT["reso"]["half68_rows"], 1e3 * OUT["reso"]["rms"]))

# =============================================================================== 4. nominal fit
print("\n4. NOMINAL FIT (%.0f s)" % (time.time() - T0))
SIG = 0.00670
T = V.step(SIG)
def run(N_, CB_, CBC_, T_, ysim, Kp=K, aP=V.AP, **kw):
    R = V.gls_fit(N_, CB_, CBC_, T_, fixed=aP * Kp, **kw)
    return R, R.h / V.proj(R, ysim)
YS = TM["C2c"] + TM["C3"]
R0, f0 = run(N, CB, CBC, T, YS)
hs0 = V.proj(R0, YS)
lin = {k: V.proj(R0, TM[k]) for k in TM}
fu = R0.h / (lin["C2"] + lin["C3"])
OUT["nom"] = dict(hd=R0.h, ehd=R0.eh, chi2=R0.chi2, ndf=R0.ndf, chi20=R0.chi20, hs=hs0, f=f0, statfit=R0.eh / hs0, f_unc=fu, lin=lin,
                  LhK=V.proj(R0, K), a_k48=R0.h / (lin["C2c"] + 0.48 * lin["C3"]))
print("   h_data = %.2f +- %.2f; chi2 %.2f/%d; without step %.1f (delta %.1f); h_sim = %.2f; f = %.4f +- %.4f (fit); uncorrected %.4f; kappa 0.48 %.4f" %
      (R0.h, R0.eh, R0.chi2, R0.ndf, R0.chi20, R0.chi20 - R0.chi2, hs0, f0, R0.eh / hs0, fu, OUT["nom"]["a_k48"]))
print("   linear projections:", {k: round(v, 1) for k, v in lin.items()}, " Lh.K (peak, a_P = 1) = %.2f rows" % OUT["nom"]["LhK"])
# own resolution sigma
for sg_, lab in [(s_a[0], "core40"), (s_b[0], "core2s")]:
    Rx, fx = run(N, CB, CBC, V.step(sg_), YS)
    OUT["nom"]["f_sig_" + lab] = fx
    print("   with own sigma %s = %.2f MeV: h_data %.2f, f %.4f" % (lab, 1e3 * sg_, Rx.h, fx))
# a_P free
Rf = V.gls_fit(N, CB, CBC, T, extra=[K])
OUT["nom"]["aP_free"] = float(Rf.p[3]); OUT["nom"]["f_aP_free"] = Rf.h / V.proj(Rf, YS)
print("   a_P free: a_P = %.3f, f %.4f" % (Rf.p[3], OUT["nom"]["f_aP_free"]))
# free step width on the data
chs = []
sgrid = np.arange(0.0040, 0.0111, 0.0002)
for sg_ in sgrid:
    chs.append(V.gls_fit(N, CB, CBC, V.step(sg_), fixed=V.AP * K).chi2)
chs = np.array(chs); i = int(np.argmin(chs)); pf = np.polyfit(sgrid[max(0, i - 6):i + 7], chs[max(0, i - 6):i + 7], 2)
OUT["nom"]["free_sigma"] = [float(-pf[1] / (2 * pf[0])), float(1 / math.sqrt(pf[0]))]
print("   free step width on the data: %.2f +- %.2f MeV" % (1e3 * OUT["nom"]["free_sigma"][0], 1e3 * OUT["nom"]["free_sigma"][1]))
np.savez(V.DIR + "vedgeu_spec.npz", x=V.X, N=N, CB=CB, pp=Hpp.sum(0), mm=Hmm.sum(0), **{"t_" + k: v for k, v in TM.items()},
         fit_idx=R0.idx, fit_smooth=R0.smooth, fit_y=R0.y)

# =============================================================================== 5. bootstrap (stat) and replicas (MC)
print("\n5. BOOTSTRAP AND REPLICAS (%.0f s)" % (time.time() - T0))
NBT = 200
rows = np.where(s0d & (bD >= 0) & (sD >= 0))[0]
Wb = V.group_weights(D["k"][rows], D["grp"][rows], NBT, 0xBEE5ED6E)
kk = D["k"][rows]
HB = [V.hist_reps(bD[rows][kk == k], sD[rows][kk == k], Wb[:, kk == k]) for k in range(3)]
del Wb
fb, hb = [], []
for r in range(NBT):
    Nb, CBb, CBCb = V.spectrum(HB[0][r], HB[1][r], HB[2][r])
    Rb, fbb = run(Nb, CBb, CBCb, T, YS)
    fb.append(fbb); hb.append(Rb.h)
fb = np.array(fb)
OUT["nom"]["statboot"] = V.sd(fb); OUT["nom"]["boot_mean"] = float(fb.mean())
print("   stat (bootstrap, %d) = %.4f (fit %.4f, ratio %.2f); mean %.4f" % (NBT, V.sd(fb), R0.eh / hs0, V.sd(fb) / (R0.eh / hs0), fb.mean()))
NRS = 100
rs = np.where(s0s & (gC <= 4) & (bS >= 0) & (sS >= 0))[0]
Wr = V.group_weights(S["tree"][rs], S["grp"][rs], NRS, 0x5EED7E57)
def rep_t(g, c2=True):
    m = gC[rs] == g
    ww = wS[rs][m] * (cfac[rs][m] if (g == 0 and c2) else 1.0)
    return V.hist_reps(bS[rs][m], sS[rs][m], Wr[:, m] * ww[None, :]).sum(1)
RT = {nm: rep_t(g) for g, nm in enumerate(["C2c", "C3", "C4", "C5", "Geta"])}
RT["C2"] = rep_t(0, c2=False)
del Wr
fr, hsr, fur = [], [], []
for r in range(NRS):
    yd = N - CB - V.AP * (RT["C4"][r] + RT["C5"][r])
    hd_r = V.proj(R0, yd)
    hs_r = V.proj(R0, RT["C2c"][r] + RT["C3"][r])
    fr.append(hd_r / hs_r); hsr.append(hs_r); fur.append(hd_r / V.proj(R0, RT["C2"][r] + RT["C3"][r]))
fr = np.array(fr); hsr = np.array(hsr)
geta_r = np.array([V.proj(R0, RT["Geta"][r]) for r in range(NRS)])
OUT["nom"].update(mc=V.sd(fr), mc_mean=float(fr.mean()), hs_rep_rms=V.sd(hsr), hs_rep_mean=float(hsr.mean()), mc_unc=V.sd(fur),
                  geta=lin["Geta"], geta_rep_rms=V.sd(geta_r))
print("   MC (own replicas, %d) = %.4f (mean %.4f); h_sim replica rms %.1f (%.1f %%), mean %.1f; uncorrected MC %.4f" %
      (NRS, V.sd(fr), fr.mean(), V.sd(hsr), 100 * V.sd(hsr) / hs0, hsr.mean(), V.sd(fur)))
print("   G_eta step projection %.1f +- %.1f (own replica rms)" % (lin["Geta"], V.sd(geta_r)))

# =============================================================================== 6. variants
print("\n6. VARIANTS (%.0f s)" % (time.time() - T0))
def boot_delta(fun, nb=NBT):
    """bootstrap noise of a data-side variant: fun(Hus, Hpp, Hmm) -> f"""
    d = []
    for r in range(nb):
        d.append(fun(HB[0][r], HB[1][r], HB[2][r]) - fb[r])
    return V.sd(d)
VAR = {}
# (a) CB = N++ + N--
fcb = lambda a, b, c: run(*V.spectrum(a, b, c, mode="sum"), T, YS)[1]
VAR["cbsum"] = dict(f=fcb(Hus, Hpp, Hmm)); VAR["cbsum"]["noise_boot"] = boot_delta(fcb)
# (b) r = 0.89 / 1.10 (sim M4 resolution; bump and peak templates)
for rr in (0.89, 1.10):
    ur, mr_ = V.u_sim(S, r=rr)
    br, sr = V.ubin(ur), V.slice_index(mr_)
    TR = {nm: stmpl(g, b=br, s=sr) for g, nm in enumerate(["C2c", "C3", "C4", "C5"])}
    Rr, frr = run(N, CB, CBC, T, TR["C2c"] + TR["C3"], Kp=TR["C4"] + TR["C5"])
    VAR["r%.2f" % rr] = dict(f=frr, hs=V.proj(Rr, TR["C2c"] + TR["C3"]))
# (c) own smearing per sample: data template sigma x 5.44/6.09, sim template sigma nominal (sim fitted with V of the data fit)
Td = V.step(SIG * 5.44 / 6.09)
Rd = V.gls_fit(N, CB, CBC, Td, fixed=V.AP * K)
Rs = V.gls_fit(YS, np.zeros(V.NB), np.zeros((V.NB, V.NB)), T, Vfix=Rd.V)
VAR["sigData"] = dict(f=Rd.h / Rs.h, hd=Rd.h, hs=Rs.h)
# (d) q = 0 (E(m, 0) for data and simulation)
ud0 = V.u_data(D, useq=False); bD0 = V.ubin(ud0)
H0 = [V.hist(bD0[s0d & (D["k"] == k)], sD[s0d & (D["k"] == k)]) for k in range(3)]
us_q0, mpq0 = V.u_sim(S, useq=False); bq0 = V.ubin(us_q0)
TQ = {nm: stmpl(g, b=bq0) for g, nm in enumerate(["C2c", "C3", "C4", "C5"])}
VAR["q0"] = dict(f=run(*V.spectrum(*H0), T, TQ["C2c"] + TQ["C3"], Kp=TQ["C4"] + TQ["C5"])[1])
# (e) pol3, 5 MeV CB cells, 20 MeV CB cells
VAR["pol3"] = dict(f=run(N, CB, CBC, T, YS, order=3)[1])
VAR["cb5"] = dict(f=run(*V.spectrum(Hus, Hpp, Hmm, reb=2), T, YS)[1])
VAR["cb20"] = dict(f=run(*V.spectrum(Hus, Hpp, Hmm, reb=8), T, YS)[1])
# (f) T12 (h_sim from all simulated classes, the bedge convention)
VAR["T12"] = dict(f=R0.h / (hs0 + lin["Geta"] + lin["GD"] + lin["GC"]))
for k, v in VAR.items():
    v["delta_pct"] = 100 * (v["f"] / f0 - 1)
    print("   %-8s f %.4f (%+.2f %%)%s" % (k, v["f"], v["delta_pct"], ("  boot noise of delta %.4f" % v["noise_boot"]) if "noise_boot" in v else ""))
OUT["var"] = VAR

# =============================================================================== 7. null tests
print("\n7. NULL TESTS (%.0f s)" % (time.time() - T0))
NT = {}
idx = R0.idx
for lab, H in [("pp", Hpp.sum(0)), ("mm", Hmm.sum(0)), ("pp+mm", Hpp.sum(0) + Hmm.sum(0))]:
    h = V.proj(R0, H); e = math.sqrt(float(R0.Lh @ np.diag(np.maximum(H[idx], 1.0)) @ R0.Lh))
    Rn = V.gls_fit(H.astype(float), np.zeros(V.NB), np.zeros((V.NB, V.NB)), T)   # own fit, Poisson-like variance
    NT[lab] = dict(lin=h, lin_err=e, fit=Rn.h, fit_err=Rn.eh, chi2=Rn.chi2, chi20=Rn.chi20)
    print("   %-6s linear estimator %.1f +- %.1f; own fit %.1f +- %.1f (chi2 %.1f/%d, without step %.1f)" % (lab, h, e, Rn.h, Rn.eh, Rn.chi2, Rn.ndf, Rn.chi20))
NT["cb_sub"] = V.proj(R0, CB)
print("   subtracted CB: step projection %.1f rows (%.1f %% of h_data)" % (NT["cb_sub"], 100 * NT["cb_sub"] / R0.h))
# e-e- by M(ee) class and by period
eeD = D["ee"]
per = np.full(len(D["fid"]), -1)
per[(D["fid"] >= 34) & (D["fid"] <= 46)] = 0; per[(D["fid"] >= 47) & (D["fid"] <= 48)] = 1
per[(D["fid"] >= 49) & (D["fid"] <= 53)] = 2; per[(D["fid"] >= 59) & (D["fid"] <= 67)] = 3
for lab, mm_ in [("mm_E1", eeD < 0.010), ("mm_E2", eeD >= 0.010), ("pp_E1", eeD < 0.010), ("pp_E2", eeD >= 0.010),
                 ("mm_perI", per == 0), ("mm_perIIA", per == 1), ("mm_perIIB", per == 2), ("mm_perIII", per == 3),
                 ("mm_preNotS0", None)]:
    k = 2 if lab.startswith("mm") else 1
    if mm_ is None:
        msk = ((D["bits"] & 2) != 0) & ~s0d & (D["k"] == k)
    else:
        msk = s0d & mm_ & (D["k"] == k)
    H = V.hist(bD[msk], sD[msk]).sum(0)
    h = V.proj(R0, H); e = math.sqrt(float(R0.Lh @ np.diag(np.maximum(H[idx], 1.0)) @ R0.Lh))
    NT[lab] = dict(lin=h, lin_err=e)
    print("   %-12s %.1f +- %.1f" % (lab, h, e))
# sim non-eta backgrounds: step projection with own replicas (G_D, G_C)
rb = np.where(s0s & (gC >= 5) & (bS >= 0) & (sS >= 0))[0]
Wg = V.group_weights(S["tree"][rb], S["grp"][rb], 50, 0xB6B6B6)
for g, nm in [(5, "GD"), (6, "GC")]:
    m = gC[rb] == g
    HRg = V.hist_reps(bS[rb][m], sS[rb][m], Wg[:, m] * wS[rb][m][None, :]).sum(1)
    pr = np.array([V.proj(R0, HRg[r]) for r in range(HRg.shape[0])])
    NT["sim_" + nm] = dict(lin=lin[nm], rep_rms=V.sd(pr))
    print("   SMASH %s: step projection %.1f +- %.1f (own 50-replica rms)" % (nm, lin[nm], V.sd(pr)))
del Wg
OUT["null"] = NT

# =============================================================================== 8. yields and sigma
print("\n8. YIELDS (%.0f s)" % (time.time() - T0))
m4p = V.m4p_sim(S)
inB = (m4p >= 0.4550) & (m4p < 0.5225)
c2c = np.sum(wS[(gC == 0) & inB] * cfac[(gC == 0) & inB]); c2 = np.sum(wS[(gC == 0) & inB]); c3 = np.sum(wS[(gC == 1) & inB])
c2all = np.sum(wS[gC == 0] * cfac[gC == 0]); c3all = np.sum(wS[gC == 1])
YB = f0 * (c2c + c3) * 64 / 63; Yall = f0 * (c2all + c3all) * 64 / 63
sig = f0 / 1.00436 * 1.0694 * 1.0104 * 0.97
OUT["yield"] = dict(c2c_B=c2c, c2_B=c2, c3_B=c3, sum_B_6463=(c2 + c3) * 64 / 63, corr_B_6463=(c2c + c3) * 64 / 63, YB=YB, YB_nb=YB / 6079.2,
                    Yall=Yall, sigma=sig, sigma_unc=fu / 1.00436 * 1.0694 * 1.0104 * 0.97)
print("   (C2 + C3) in B x 64/63 = %.1f; (c C2 + C3) x 64/63 = %.1f; Y_B = %.0f rows = %.3f nb; Y all M4 = %.0f; sigma = %.4f mb (uncorr. %.4f)" %
      (OUT["yield"]["sum_B_6463"], OUT["yield"]["corr_B_6463"], YB, YB / 6079.2, Yall, sig, OUT["yield"]["sigma_unc"]))

json.dump(OUT, open(V.DIR + "vedgeu_ana.json", "w"), indent=1, default=float)
print("done in %.0f s" % (time.time() - T0))
