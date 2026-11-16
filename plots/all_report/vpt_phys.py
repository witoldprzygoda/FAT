#!/usr/bin/env python3
"""vpt_phys.py -- verifier of peaktail: physical origin of the low-mass tail of the simulated eta peak (own code).

Generator line shape, reconstruction residual, substitution of true momenta (instead of a linearised Jacobian),
momentum dependence of the lepton deficit (bremsstrahlung = fractional loss), conversion pairs of C4.
Output: stdout (tee to /hdd3/przygoda/HADES/FAT_out/all_report/vpt/vpt_phys.log) and vpt_phys.npz.
"""
import sys
import numpy as np
sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import vpt_core as V


def wq(x, w, qs):
    o = np.argsort(x); x = x[o]; c = np.cumsum(w[o]); c /= c[-1]
    return [float(np.interp(q, c, x)) for q in qs]


def frac(m, w, sel=None):
    s = np.ones(len(w), bool) if sel is None else sel
    return float(w[s & m].sum() / w[s].sum())


S = V.load_sim()
w = V.sim_weight(S, c2=False)
mpp = V.SPP * S["pp"].astype(float)
reg = (mpp.astype(np.float32) >= np.float32(0.30)) & (mpp.astype(np.float32) < np.float32(0.40))
cls = {"C4": S["g"] == 2, "C5": S["g"] == 3}
ok = w > 0
print("=" * 100)
print("A. generator: true four-body mass m4_true [MeV/c2], weighted quantiles 0.1 / 1 / 50 / 99 / 99.9 %")
for k, m in cls.items():
    for rn, rr in (("0.30-0.40", reg), ("all M(pipi)", np.ones(len(w), bool))):
        s = m & ok & rr & (S["m4t"] > 0)
        q = wq(S["m4t"][s].astype(float) * 1e3, w[s], [0.001, 0.01, 0.5, 0.99, 0.999])
        lo1 = frac(S["m4t"][s] * 1e3 < q[2] - 1.0, w[s]); lo5 = frac(S["m4t"][s] * 1e3 < q[2] - 5.0, w[s])
        print("  %s %-11s rows %7d  N_eff %8.0f  q: %s   below median -1 MeV: %.4f %%, -5 MeV: %.4f %%" %
              (k, rn, s.sum(), w[s].sum() ** 2 / (w[s] ** 2).sum(), " / ".join("%.3f" % x for x in q), 100 * lo1, 100 * lo5))
# ---------------------------------------------------------------- four-vectors
mass = {"ep": V.ME_E, "em": V.ME_E, "pip": V.MPI, "pim": V.MPI}
vr = {n: V.p4(S[n + "_p"], S[n + "_th"], S[n + "_ph"], mass[n]) for n in mass}
vt = {n: V.p4(S[n + "_pt"], S[n + "_th"], S[n + "_ph"], mass[n]) for n in mass}     # true |p|, reco angles
vu = {n: V.p4(S[n + "_pr"], S[n + "_th"], S[n + "_ph"], mass[n]) for n in mass}     # uncorrected |p|
m4r = V.minv(vr["ep"], vr["em"], vr["pip"], vr["pim"])
sel45 = (S["g"] >= 2) & (S["g"] <= 3) & ok & (S["m4t"] > 0)
print("\nB. M4 from (p, theta, phi): max |M4(rebuilt) - m4| over C4+C5 rows = %.3g MeV/c2 (corrected p); uncorrected p: mean shift %.2f MeV/c2"
      % (np.abs(m4r[sel45] - S["m4"][sel45]).max() * 1e3,
         1e3 * np.average((V.minv(vu["ep"], vu["em"], vu["pip"], vu["pim"]) - m4r)[sel45], weights=w[sel45])))
m4_lt = V.minv(vt["ep"], vt["em"], vr["pip"], vr["pim"])     # leptons with true |p|
m4_pt = V.minv(vr["ep"], vr["em"], vt["pip"], vt["pim"])     # pions with true |p|
m4_at = V.minv(vt["ep"], vt["em"], vt["pip"], vt["pim"])     # all four with true |p| (reco angles)
s = V.S_NOM
d = s * (S["m4"].astype(float) - S["m4t"].astype(float))     # M4' - s m4_true (r = 1)
dlep = s * (m4r - m4_lt); dpi = s * (m4r - m4_pt); dang = s * (m4_at - S["m4t"].astype(float))
print("\nC. residual d = M4' - 0.99615 m4_true and its sources by SUBSTITUTION (true |p| put in, reco angles kept), 0.30-0.40")
for k, m in cls.items():
    ss = m & sel45 & reg
    ww = w[ss]; dd = d[ss] * 1e3
    core = np.abs(dd) < 12
    mu = np.average(dd[core], weights=ww[core])
    # Gaussian core sigma from a fit to the central part (|d| < 12 MeV): truncated-Gaussian moment matching
    from scipy.optimize import brentq
    from scipy.stats import norm
    v = np.average((dd[core] - mu) ** 2, weights=ww[core])
    def tr(sg):
        a = 12 / sg
        return sg * sg * (1 - 2 * a * norm.pdf(a) / (2 * norm.cdf(a) - 1)) - v
    sg = brentq(tr, 1.0, 30.0)
    out = "  %s rows %6d N_eff %7.0f  core mu %.2f sigma %.2f MeV | below -10/-15/-25/-40: %s %% | above +10/+15: %s %%" % (
        k, ss.sum(), ww.sum() ** 2 / (ww ** 2).sum(), mu, sg,
        " / ".join("%.1f" % (100 * frac(dd < -x, ww)) for x in (10, 15, 25, 40)),
        " / ".join("%.1f" % (100 * frac(dd > x, ww)) for x in (10, 15)))
    print(out)
    for thr in (-15, -25):
        t = dd < thr
        print("     rows d < %d MeV (%.1f %%): mean d %.1f = leptons %.1f + pions %.1f + angles/rest %.1f (non-additive rest %.1f); "
              "a lepton p/p_true < 0.90: %.0f %%" % (
                  thr, 100 * frac(t, ww), np.average(dd[t], weights=ww[t]), 1e3 * np.average(dlep[ss][t], weights=ww[t]),
                  1e3 * np.average(dpi[ss][t], weights=ww[t]), 1e3 * np.average(dang[ss][t], weights=ww[t]),
                  np.average(dd[t], weights=ww[t]) - 1e3 * np.average((dlep + dpi + dang)[ss][t], weights=ww[t]),
                  100 * frac(np.minimum(S["ep_p"][ss] / S["ep_pt"][ss], S["em_p"][ss] / S["em_pt"][ss])[t] < 0.90, ww[t])))
    tl = dlep[ss] * 1e3 < -5
    print("     substitution tail (lepton |p| deficit shifts M4 by more than 5 MeV/c2): %.1f %% of the rows; rows without: "
          "below -10/-15/-20 %s %%, above +10/+15/+20 %s %%" % (
              100 * frac(tl, ww), " / ".join("%.2f" % (100 * frac(dd[~tl] < -x, ww[~tl])) for x in (10, 15, 20)),
              " / ".join("%.2f" % (100 * frac(dd[~tl] > x, ww[~tl])) for x in (10, 15, 20))))
# ---------------------------------------------------------------- lepton momentum deficit vs momentum
print("\nD. lepton momentum ratio p_reco/p_true (corrected p), C4 + C5 S0 mode A, all M(pipi); bremsstrahlung is a FRACTIONAL loss")
pr = np.concatenate([S["ep_p"][sel45] / S["ep_pt"][sel45], S["em_p"][sel45] / S["em_pt"][sel45]]).astype(float)
pt = np.concatenate([S["ep_pt"][sel45], S["em_pt"][sel45]]).astype(float)
pu = np.concatenate([S["ep_pr"][sel45] / S["ep_pt"][sel45], S["em_pr"][sel45] / S["em_pt"][sel45]]).astype(float)
ww2 = np.concatenate([w[sel45], w[sel45]])
for lo, hi in ((0, 100), (100, 150), (150, 200), (200, 300), (300, 500), (500, 2000)):
    m = (pt >= lo) & (pt < hi)
    if m.sum() < 20:
        continue
    print("  p_true %4d-%4d MeV/c: leptons %6d  median ratio %.4f (uncorr. %.4f)  P(ratio < 0.95 / 0.90 / 0.80) = %.2f / %.2f / %.2f %%  "
          "mean (1 - ratio | ratio < 0.9) %.3f" % (
              lo, hi, m.sum(), wq(pr[m], ww2[m], [0.5])[0], wq(pu[m], ww2[m], [0.5])[0],
              100 * frac(pr[m] < 0.95, ww2[m]), 100 * frac(pr[m] < 0.90, ww2[m]), 100 * frac(pr[m] < 0.80, ww2[m]),
              np.average(1 - pr[m][pr[m] < 0.9], weights=ww2[m][pr[m] < 0.9])))
print("  all: P(ratio < 0.90) = %.2f %%, P(ratio > 1.10) = %.2f %% (resolution is symmetric; the low side is the loss)" % (
    100 * frac(pr < 0.90, ww2), 100 * frac(pr > 1.10, ww2)))
# effective radiator thickness from the integrated Bethe-Heitler spectrum (complete screening)
def bh(y0):
    return (4. / 3.) * np.log(1 / y0) - (4. / 3.) * (1 - y0) + 0.5 * (1 - y0 * y0)
for y0 in (0.10, 0.20, 0.30):
    print("  P(1 - ratio > %.2f) = %.2f %%  -> effective x/X0 = %.2f %% (Bethe-Heitler, P = x/X0 [4/3 ln(1/y) - 4/3 (1-y) + (1-y^2)/2])"
          % (y0, 100 * frac(pr < 1 - y0, ww2), 100 * frac(pr < 1 - y0, ww2) / bh(y0)))
# ---------------------------------------------------------------- C4 conversion pairs
print("\nE. C4 (eta -> pi+pi-gamma, gamma converts): lepton creation, process and parent ids; pair mass at creation")
c4 = cls["C4"] & sel45
for n in ("ep", "em"):
    u, c = np.unique(S[n + "_proc"][c4], return_counts=True)
    up, cp = np.unique(S[n + "_par"][c4], return_counts=True)
    print("  %s: process ids %s; parent ids (top) %s" % (n, dict(zip(u.tolist(), c.tolist())),
                                                     dict(sorted(zip(up.tolist(), cp.tolist()), key=lambda x: -x[1])[:5])))
rr = np.maximum(S["ep_vr"][c4], 0).astype(float)
print("  creation radius r (e+) 5 / 50 / 95 %% = %s mm; z %s mm; |vz(e+) - vz(e-)| max %.3g mm" % (
    " / ".join("%.1f" % x for x in wq(rr, w[c4], [0.05, 0.5, 0.95])),
    " / ".join("%.0f" % x for x in wq(S["ep_vz"][c4].astype(float), w[c4], [0.05, 0.5, 0.95])),
    np.abs(S["ep_vz"][c4] - S["em_vz"][c4]).max()))
print("  true pair mass at creation M(ee)_true 50 / 95 / 99 %% = %s MeV/c2" % " / ".join(
    "%.2f" % x for x in wq(S["eet"][c4].astype(float) * 1e3, w[c4], [0.5, 0.95, 0.99])))
print("  m_eta(generator, C5 median) - m4_true(C4): 50 / 99.9 %% = %s MeV/c2 (the e+e- pair carries the photon four-momentum up to the nuclear recoil)" % (
    " / ".join("%.2f" % (548.0 - x) for x in wq(S["m4t"][c4].astype(float) * 1e3, w[c4], [0.5, 0.001]))))
c5 = cls["C5"] & sel45
for n in ("ep", "em"):
    u, c = np.unique(S[n + "_proc"][c5], return_counts=True)
    print("  C5 %s process ids %s" % (n, dict(zip(u.tolist(), c.tolist()))))
np.savez(V.OUT + "vpt_phys.npz", d=d, dlep=dlep, dpi=dpi, dang=dang)
