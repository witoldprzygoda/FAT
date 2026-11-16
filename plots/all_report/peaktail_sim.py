#!/usr/bin/env python3
# ============================================================================
# peaktail_sim.py  --  package peaktail, task 1: where the low-mass tail of the simulated eta peak comes from, and where
#                      the peak rows lie in u = M4 - E(M(pipi), M(ee))
# ============================================================================
# Written 2026-11-15 (prefix peaktail). Reads the columns of peaktail_dump.C. Writes
#   /hdd3/przygoda/HADES/FAT_out/all_report/peaktail/peaktail_sim.json, peaktail_sim.log, peaktail_simfig.npz
# Classes: C4 = eta -> pi+pi-gamma, photon converted (std set, 307 rows); C5 = eta -> pi+pi-e+e- (eta set, 571 238 rows).
# Selection: S0, mode A, w x V3 (108), M4' = 0.99615 M4 (r = 1), M(pipi)' = 0.9952 M(pipi).
# Run: python3 -I plots/all_report/peaktail_sim.py
# ============================================================================
import json
import math
import sys

import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import edgeu_core as EC
import peaktail_core as P

LOG = open(P.OUT + "peaktail_sim.log", "w")
def L(*a):
    s = " ".join(str(x) for x in a)
    print(s); LOG.write(s + "\n"); LOG.flush()
RES = {}
S = P.load_sim()
w = P.sim_weight(S)
g = S["g"]
u, mpp, m4p = EC.u_sim(S)
d = P.resid(S)
E1 = S["mee"] < 0.010
L("# peaktail_sim.py: %d simulated rows" % S["n"])

# ------------------------------------------------------------------------------------------------ 1. the generator
L("\n1. GENERATED four-body mass m4_true (weighted quantiles 0.1 / 1 / 50 / 99 / 99.9 %; MeV)")
gen = {}
for gg, nm in ((2, "C4"), (3, "C5")):
    m = (g == gg) & (w > 0)
    mt = S["m4t"][m].astype(np.float64)
    q = P.wq(mt, w[m], [0.001, 0.01, 0.5, 0.99, 0.999])
    low1 = P.frac(mt < q[2] - 0.001, w[m]); low5 = P.frac(mt < q[2] - 0.005, w[m])
    gen[nm] = dict(rows=int(m.sum()), q=[1e3 * x for x in q], below_1MeV=low1, below_5MeV=low5, min=1e3 * float(mt.min()))
    L("   %s: %d rows; %s; fraction more than 1 / 5 MeV/c2 below the median: %.2e / %.2e; minimum %.2f" %
      (nm, m.sum(), " ".join("%.3f" % (1e3 * x) for x in q), low1, low5, 1e3 * mt.min()))
L("   -> the generated four-body mass is a delta function in both classes: the generator writes no radiative photon (C5),")
L("      and the conversion pair carries the whole photon four-momentum at its creation point (C4).")
RES["generator"] = gen

# ------------------------------------------------------------------------------------------------ 2. the residual
def gcore(dd, ww):
    h, e = np.histogram(dd, bins=np.arange(-0.012, 0.0121, 0.001), weights=ww)
    h2, _ = np.histogram(dd, bins=e, weights=ww * ww)
    x = 0.5 * (e[1:] + e[:-1])
    f = lambda x, a, mu, s: a * np.exp(-0.5 * (x - mu) ** 2 / s ** 2)
    p, cv = curve_fit(f, x, h, p0=[h.max(), -0.002, 0.006], sigma=np.sqrt(np.maximum(h2, h2[h2 > 0].min())), absolute_sigma=True)
    return float(p[1]), float(abs(p[2])), float(math.sqrt(cv[2, 2]))
L("\n2. RECONSTRUCTION RESIDUAL d = M4' - 0.99615 m4_true (MeV/c2)")
REG = {"all": (0.26, 0.62), "fit 0.30-0.40": (0.30, 0.40), "0.40-0.46": (0.40, 0.46), "0.42-0.55": (0.42, 0.55), "0.46-0.55": (0.46, 0.55)}
res = {}
TM = {}
for gg, nm in ((2, "C4"), (3, "C5"), (-1, "C4+C5")):
    for rn, (a, b) in REG.items():
        m = ((g == gg) if gg >= 0 else ((g == 2) | (g == 3))) & (w > 0) & (mpp >= a) & (mpp < b)
        if w[m].sum() <= 0:
            continue
        dd, ww = d[m], w[m]
        q = P.wq(dd, ww, [0.01, 0.05, 0.16, 0.5, 0.84, 0.95, 0.99])
        fr = {("<%d" % k): P.frac(dd < -k * 1e-3, ww) for k in (10, 15, 20, 25, 30, 40, 50)}
        fr.update({(">+%d" % k): P.frac(dd > k * 1e-3, ww) for k in (10, 15, 20)})
        try:
            mu, sg, esg = gcore(dd, ww)
        except Exception:
            mu, sg, esg = float("nan"), float("nan"), float("nan")
        tm = P.TailModel(dd, ww)
        TM[(nm, rn)] = tm
        neff = ww.sum() ** 2 / (ww * ww).sum()
        res["%s|%s" % (nm, rn)] = dict(rows=int(m.sum()), sumw=float(ww.sum()), neff=float(neff), q=[1e3 * x for x in q], frac=fr, core_mu=1e3 * mu,
                                       core_sigma=1e3 * sg, e_core_sigma=1e3 * esg, mirror_mu=1e3 * tm.mu, ftail=tm.ftail)
        L("   %-5s %-14s rows %7d (N_eff %7.0f), sum %8.1f: quantiles 1/5/16/50/84/95/99 %% = %s; core mu %.2f sigma %.2f +- %.2f;" %
          (nm, rn, m.sum(), neff, ww.sum(), " ".join("%+.1f" % (1e3 * x) for x in q), 1e3 * mu, 1e3 * sg, 1e3 * esg))
        L("         below -10/-15/-20/-25/-30/-40 MeV: %s; above +10/+15/+20: %s; mirror tail fraction %.4f (mode %.2f)" %
          (" ".join("%.4f" % fr["<%d" % k] for k in (10, 15, 20, 25, 30, 40)), " ".join("%.4f" % fr[">+%d" % k] for k in (10, 15, 20)), tm.ftail, 1e3 * tm.mu))
RES["resid"] = res

# ------------------------------------------------------------------------------------------------ 3. sources
L("\n3. SOURCES of the residual (linearised: d = sum_i J_i (p_i - p_i,true) + rest; J_i = dM4/d|p_i| from the reconstructed")
L("   four-vectors; rest = angles + non-linear terms); x = p / p_true of each track")
sub = {k: v for k, v in S.items()}
M4r, J = P.m4_and_jacobian(sub)
dp = [(S[a + "_p"].astype(np.float64) - S[a + "_ptrue"].astype(np.float64)) * 1e-3 for a in ("ep", "em", "pip", "pim")]
dl = EC.S_NOM * (J[0] * dp[0] + J[1] * dp[1]); dpi = EC.S_NOM * (J[2] * dp[2] + J[3] * dp[3]); dang = d - dl - dpi
xl = np.minimum(S["ep_p"] / S["ep_ptrue"], S["em_p"] / S["em_ptrue"]).astype(np.float64)
xpi = np.minimum(S["pip_p"] / S["pip_ptrue"], S["pim_p"] / S["pim_ptrue"]).astype(np.float64)
src = {}
for gg, nm in ((2, "C4"), (3, "C5")):
    for rn in ("all", "fit 0.30-0.40"):
        a, b = REG[rn]
        base = (g == gg) & (w > 0) & (mpp >= a) & (mpp < b)
        out = {}
        for tn, tm in (("core |d| < 5", np.abs(d) < 0.005), ("tail d < -15", d < -0.015), ("tail d < -25", d < -0.025), ("high d > +15", d > 0.015)):
            m = base & tm
            ww = w[m]
            if ww.sum() <= 0:
                continue
            mean = lambda v: float(np.average(v[m], weights=ww))
            out[tn] = dict(share=P.frac(tm[base], w[base]), d=1e3 * mean(d), lep=1e3 * mean(dl), pi=1e3 * mean(dpi), ang=1e3 * mean(dang),
                           lepdom=P.frac((dl[m] < 0.5 * d[m]), ww), x_lep_lt095=P.frac(xl[m] < 0.95, ww), x_lep_lt090=P.frac(xl[m] < 0.90, ww),
                           x_pi_lt097=P.frac(xpi[m] < 0.97, ww), rows=int(m.sum()))
            L("   %-3s %-14s %-13s share %.4f (%5d rows): mean d %+6.1f = leptons %+6.1f + pions %+6.1f + angles/rest %+6.1f MeV; lepton part > half of d in %.2f;"
              " a lepton with p/p_true < 0.95 / 0.90: %.2f / %.2f; a pion < 0.97: %.2f" %
              (nm, rn, tn, out[tn]["share"], m.sum(), out[tn]["d"], out[tn]["lep"], out[tn]["pi"], out[tn]["ang"], out[tn]["lepdom"],
               out[tn]["x_lep_lt095"], out[tn]["x_lep_lt090"], out[tn]["x_pi_lt097"]))
        src["%s|%s" % (nm, rn)] = out
# lepton momentum ratio distributions (all rows of the class, both leptons)
L("   lepton p/p_true (both leptons, weighted): fraction below 0.98 / 0.95 / 0.90 / 0.80; pions below 0.98 / 0.95")
for gg, nm in ((2, "C4"), (3, "C5")):
    m = (g == gg) & (w > 0)
    xe = np.concatenate([S["ep_p"][m] / S["ep_ptrue"][m], S["em_p"][m] / S["em_ptrue"][m]]).astype(np.float64)
    xp = np.concatenate([S["pip_p"][m] / S["pip_ptrue"][m], S["pim_p"][m] / S["pim_ptrue"][m]]).astype(np.float64)
    ww = np.concatenate([w[m], w[m]])
    fe = [P.frac(xe < c, ww) for c in (0.98, 0.95, 0.90, 0.80)]; fp = [P.frac(xp < c, ww) for c in (0.98, 0.95)]
    src["x_%s" % nm] = dict(lep=fe, pi=fp, lep_median=P.wq(xe, ww, [0.5])[0], pi_median=P.wq(xp, ww, [0.5])[0])
    L("   %s: leptons %s (median %.4f); pions %s (median %.4f)" % (nm, " / ".join("%.4f" % x for x in fe), src["x_%s" % nm]["lep_median"],
                                                               " / ".join("%.4f" % x for x in fp), src["x_%s" % nm]["pi_median"]))
# creation points of the leptons
L("   lepton creation points (weighted): r [mm] quantiles 5/50/95 %, z [mm] 5/50/95 %; |ep_vz - em_vz| > 1 mm")
for gg, nm in ((2, "C4"), (3, "C5")):
    m = (g == gg) & (w > 0)
    qr = P.wq(S["ep_vr"][m].astype(np.float64), w[m], [0.05, 0.5, 0.95]); qz = P.wq(S["ep_vz"][m].astype(np.float64), w[m], [0.05, 0.5, 0.95])
    dz = P.frac(np.abs(S["ep_vz"][m] - S["em_vz"][m]) > 1.0, w[m])
    src["vtx_%s" % nm] = dict(r=qr, z=qz, dz_gt1=dz)
    L("   %s: r %s; z %s; e+ and e- created more than 1 mm apart: %.3f" % (nm, " ".join("%.2f" % x for x in qr), " ".join("%.1f" % x for x in qz), dz))
# tail fraction vs the lepton pair energy in the eta frame and vs M(pipi)
L("   tail (d < -15 MeV) and (d < -25 MeV) per M(pipi)' bin, C4+C5 and C5 alone; E*_ee = (m_eta^2 + M(ee)^2 - M(pipi)^2)/(2 m_eta)")
ppb = [0.26, 0.30, 0.34, 0.37, 0.40, 0.42, 0.46, 0.50, 0.55]
mdep = []
for a, b in zip(ppb[:-1], ppb[1:]):
    row = dict(lo=a, hi=b)
    for gs, nm in (((2, 3), "C4+C5"), ((3,), "C5"), ((2,), "C4")):
        m = np.isin(g, gs) & (w > 0) & (mpp >= a) & (mpp < b)
        if w[m].sum() <= 0:
            continue
        mid = 0.5 * (a + b) / EC.EDGE_PIPI
        row[nm] = dict(sumw=float(w[m].sum()), rows=int(m.sum()), f15=P.frac(d[m] < -0.015, w[m]), f25=P.frac(d[m] < -0.025, w[m]),
                       q16=1e3 * P.wq(d[m], w[m], [0.16])[0], ftail=P.TailModel(d[m], w[m]).ftail if m.sum() > 50 else float("nan"))
    row["Eee"] = 1e3 * (EC.ME ** 2 - (0.5 * (a + b) / EC.EDGE_PIPI) ** 2) / (2 * EC.ME)
    mdep.append(row)
    L("   M(pipi) %.2f-%.2f (E*_ee %.0f MeV): C4+C5 sum %.1f f15 %.4f f25 %.4f mirror %.4f | C5 f15 %.4f f25 %.4f | C4 rows %d f15 %.3f" %
      (a, b, row["Eee"], row.get("C4+C5", {}).get("sumw", 0), row.get("C4+C5", {}).get("f15", float("nan")), row.get("C4+C5", {}).get("f25", float("nan")),
       row.get("C4+C5", {}).get("ftail", float("nan")), row.get("C5", {}).get("f15", float("nan")), row.get("C5", {}).get("f25", float("nan")),
       row.get("C4", {}).get("rows", 0), row.get("C4", {}).get("f15", float("nan"))))
src["mpipi_dependence"] = mdep
RES["sources"] = src

# ------------------------------------------------------------------------------------------------ 4. mapping into u
L("\n4. WHERE THE PEAK ROWS LIE IN u (0.30 <= M(pipi)' < 0.40, M(ee) < m_pi0; expected raw rows; fractions of the class; tail = lepton-loss rows, dl < -5 MeV)")
UB = [(-1.0, -0.045), (-0.045, 0.0), (0.0, 0.010), (0.010, 0.020), (0.020, 0.045), (0.045, 1.0)]
tmP = P.TailModel(d[((g == 2) | (g == 3)) & (w > 0) & (mpp >= 0.30) & (mpp < 0.40)], w[((g == 2) | (g == 3)) & (w > 0) & (mpp >= 0.30) & (mpp < 0.40)])
ptail = (dl < -0.005).astype(np.float64)   # tail rows: the leptons' momentum deficit shifts M4 by more than 5 MeV (as in peaktail_ufit.py)
umap = {}
for gs, nm in (((2,), "C4"), ((3,), "C5"), ((2, 3), "C4+C5")):
    for cn, cm in (("all", np.ones(S["n"], bool)), ("E1", E1), ("E2", ~E1)):
        m = np.isin(g, gs) & cm & (w > 0) & (mpp >= 0.30) & (mpp < 0.40) & np.isfinite(u)
        if w[m].sum() <= 0:
            continue
        tot = w[m].sum(); tail = (w * ptail)[m].sum()
        fr = [float(w[m & (u >= a) & (u < b)].sum() / tot) for a, b in UB]
        frt = [float((w * ptail)[m & (u >= a) & (u < b)].sum() / max(tail, 1e-30)) for a, b in UB]
        frc = [float((w * (1 - ptail))[m & (u >= a) & (u < b)].sum() / max(tot - tail, 1e-30)) for a, b in UB]
        upk = P.wq(u[m], w[m], [0.16, 0.5, 0.84])
        umap["%s|%s" % (nm, cn)] = dict(sumw=float(tot), rows=int(m.sum()), tailfrac=float(tail / tot), frac=fr, frac_tail=frt, frac_core=frc, uq=[1e3 * x for x in upk])
        L("   %-5s %-3s sum %7.1f (%6d rows), u 16/50/84 %% %+.1f/%+.1f/%+.1f MeV; tail share %.3f; fractions in u <-45 | -45..0 | 0..10 | 10..20 | 20..45 | >45:"
          % (nm, cn, tot, m.sum(), 1e3 * upk[0], 1e3 * upk[1], 1e3 * upk[2], tail / tot))
        L("         all rows  %s" % " ".join("%.4f" % x for x in fr))
        L("         tail part %s" % " ".join("%.4f" % x for x in frt))
        L("         core part %s" % " ".join("%.4f" % x for x in frc))
RES["umap"] = umap
# the edge region |u| < 10 MeV: fraction of each class there and the share of tail rows among the peak rows there
L("   edge region -10 <= u < +10 MeV:")
for gs, nm in (((2,), "C4"), ((3,), "C5"), ((2, 3), "C4+C5")):
    for cn, cm in (("all", np.ones(S["n"], bool)), ("E1", E1), ("E2", ~E1)):
        m = np.isin(g, gs) & cm & (w > 0) & (mpp >= 0.30) & (mpp < 0.40) & np.isfinite(u)
        me = m & (u >= -0.010) & (u < 0.010)
        if w[m].sum() <= 0:
            continue
        umap["%s|%s" % (nm, cn)]["edge10"] = float(w[me].sum() / w[m].sum())
        umap["%s|%s" % (nm, cn)]["edge10_tailshare"] = float((w * ptail)[me].sum() / max(w[me].sum(), 1e-30))
        umap["%s|%s" % (nm, cn)]["edge10_rows_x_aP"] = float(EC.AP * w[me].sum())
        L("     %-5s %-3s fraction %.4f (%.0f expected data rows x a_P), of which tail rows %.2f" % (nm, cn, w[me].sum() / w[m].sum(), EC.AP * w[me].sum(),
          (w * ptail)[me].sum() / max(w[me].sum(), 1e-30)))
# rows without a lepton loss: symmetric residuals (resolution); with a lepton loss: the tail
L("   symmetry check: residual of the rows WITHOUT a lepton loss (dl > -5 MeV), 0.30-0.40: below -10/-15/-20 vs above +10/+15/+20 MeV")
for gg, nm in ((2, "C4"), (3, "C5")):
    m = (g == gg) & (w > 0) & (mpp >= 0.30) & (mpp < 0.40) & (dl >= -0.005)
    lo_ = [P.frac(d[m] < -k * 1e-3, w[m]) for k in (10, 15, 20)]; hi_ = [P.frac(d[m] > k * 1e-3, w[m]) for k in (10, 15, 20)]
    RES.setdefault("symmetry", {})[nm] = dict(low=lo_, high=hi_)
    L("     %s: %s  vs  %s" % (nm, " / ".join("%.4f" % x for x in lo_), " / ".join("%.4f" % x for x in hi_)))
RES["tail_model_fit"] = dict(mu=1e3 * tmP.mu, ftail=tmP.ftail)
L("   mirror tail model of C4+C5 in 0.30-0.40: mode %.2f MeV, tail fraction %.4f" % (1e3 * tmP.mu, tmP.ftail))
# the peak position in u, row by row: u_p = 0.5439 - E (data scale); the tail reaches u = 0 when the row lost u_p
m = ((g == 2) | (g == 3)) & (w > 0) & (mpp >= 0.30) & (mpp < 0.40) & np.isfinite(u)
up = P.MPEAK - EC.edge(mpp, S["mee"].astype(np.float64))
L("   peak position u_p = 0.5439 - E per row: 5/50/95 %% = %s MeV; a row reaches u < 0 when it lost more than u_p" %
  " ".join("%.1f" % (1e3 * x) for x in P.wq(up[m], w[m], [0.05, 0.5, 0.95])))
RES["u_peak_position"] = [1e3 * x for x in P.wq(up[m], w[m], [0.05, 0.5, 0.95])]

json.dump(RES, open(P.OUT + "peaktail_sim.json", "w"), indent=1)
# figure data
fd = {}
for gg, nm in ((2, "C4"), (3, "C5")):
    for rn in ("all", "fit 0.30-0.40", "0.42-0.55"):
        a, b = REG[rn]
        m = (g == gg) & (w > 0) & (mpp >= a) & (mpp < b)
        k = nm + "_" + rn.split()[-1].replace(".", "").replace("-", "_")
        fd["d_" + k] = d[m]; fd["w_" + k] = w[m]
        fd["dl_" + k] = dl[m]; fd["dpi_" + k] = dpi[m]; fd["dang_" + k] = dang[m]
    m = (g == gg) & (w > 0) & (mpp >= 0.30) & (mpp < 0.40) & np.isfinite(u)
    fd["u_" + nm] = u[m]; fd["uw_" + nm] = w[m]; fd["upt_" + nm] = ptail[m]; fd["ue1_" + nm] = E1[m]
fd["tm_x"] = tmP.x; fd["tm_h"] = tmP.hs; fd["tm_core"] = tmP.core; fd["tm_tail"] = tmP.tail; fd["tm_p"] = tmP.ptail
np.savez_compressed(P.OUT + "peaktail_simfig.npz", **fd)
L("\nwritten: peaktail_sim.json, peaktail_sim.log, peaktail_simfig.npz")
