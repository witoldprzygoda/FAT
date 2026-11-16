#!/usr/bin/env python3
# ============================================================================
# prevstyle_results_xs.py  --  package prevstyle of report_all, section 10: what the earlier analysis' approach gives for
#                              sigma(pp -> eta X) and for the X17 upper limit in our data
#                              DEMONSTRATION in the style of the earlier analysis (not a measurement; no result of the report uses it)
# ============================================================================
# Written 2026-11-16. New file. prevstyle_fill.C, prevstyle_ana.py, prevstyle_twoscale.py and their outputs are only read.
# Inputs : /hdd3/przygoda/HADES/FAT_out/all_report/prevstyle/prevstyle_data.json, prevstyle_sim.json (prevstyle_fill.C)
#          /hdd3/przygoda/HADES/FAT_out/all_report/prevstyle/prevstyle_numbers.json (section 6: k, global fit, X17 residuals)
#          /hdd3/przygoda/HADES/FAT_out/all_report/prevstyle/prevstyle_twoscale_numbers.json (section 9: two-scale demonstration)
# Outputs: /hdd3/przygoda/HADES/FAT_out/all_report/prevstyle/prevstyle_xs_numbers.json, prevstyle_xs.log,
#          prevstyle_xs_trolke_in.txt, prevstyle_xs_trolke_out.txt
# Macro  : plots/all_report/prevstyle_xs_trolke.C (ROOT TRolke, model 5), called once with all configurations
# Selection everywhere: pre + E1 (alpha_lab < 50 deg, alpha_eta > 140 deg, M(pipi) < 0.46 GeV/c2, OA(e+e-) > 2 deg);
#   "all" = isBest >= 0 (all stored combinations, NOMINAL), "best" = isBest == 1 (best combination only, earlier slides; comparison);
#   CB = 2 sqrt(N++ N--) from the like-sign data (k = 1), never from the simulation.
#
# TASK A (eta): the two scales of the two-scale demonstration (section 9: bump f_b, peak a_p) -> sigma(pp -> eta X) with the chain of
#   the report (no new correction):
#     bump: sigma = f / 1.00436 x 1.0694 mb x C_sec(B) 1.0104 x (1 - f_eta) 0.97          (report: 0.42574 -> 0.444 mb)
#     peak: sigma = a / 1.0224  x 1.0694 mb x C_sec(P) 1.0116 x (1 - f_eta) 0.97          (report: 0.418249 -> 0.429 mb)
#   errors: data statistics and SMASH statistics (split of the fit covariance, sandwich with the fit weights), c_b (section 9),
#   method (normalisation windows, fit ranges, background tilt of section 9; stated separately), normalisation as for the headline
#   (L 6.9 % + C_sec(1 - f_eta) +3.1/-4.1 % + PDG BRs 3.15 % (peak) / 3.17 % (bump) + period 5 %), production model +0.08/-0.03 mb
#   (relative +18.6/-7.0 %) separate.
# TASK B (X17): M(pi+pi-e+e-) eta window 530-557 MeV/c2 (Wpk), side bands 518-530 + 557-569 MeV/c2; X17 windows in M(e+e-):
#   15.5-24.5 (ours) and 12-22 MeV/c2 (earlier analysis); exact cells of prevstyle_fill.C.
#     expected counts in the window without X17:  B = CB_W + k (N - CB)_SB + a (C4 + C5)_W     (section-6 convention)
#     observed: x = N_W (integer);  residual = x - B
#     k: 1.064 (global fit), 1.075 (smoothed-template fit), 27/24 (width ratio); isBest == 1: 1.0595, 1.0865, 27/24
#     a: 0.418 (peak side bands), two-scale peak (0.652; best 0.689), global fit (0.697; best 0.743)
#     sigma_B^2 = sigma(CB_W)^2 + k^2 var(SB) + SB^2 sigma_k^2 (global-fit error of k) + (sigma_a (C4+C5)_W)^2 + a^2 var_MC(C4+C5)_W
#   90 % CL upper limits on the X17 counts: (i) TRolke model 5 (Poisson count, Gaussian background, known efficiency), bounded;
#   (ii) CLs for Poisson counting, background marginalised over its Gaussian (Gauss-Hermite, 64 nodes).
#   The X17 signal also leaks into the side bands: the counts that a BR produces in the residual are BR x eta x (C6_W - k C6_SB).
#   Conversion to BR(eta -> pi+pi- X17) x BR(X17 -> e+e-), efficiencies from SMASH (X17 at 17 MeV/c2, C6, per unit BR):
#     (a) to the eta Dalitz bump: BR_X BR_ee = BR(eta -> pi+pi-pi0, PDG 2026) x (eps_D / eps_X) x N_X / N_D,
#         N_D = f_b (C2 + C3) in the bump window B (two-scale demonstration), eps_D per eta -> pi+pi-pi0 decay
#     (b) to the direct decays: BR_X BR_ee = BR(eta -> pi+pi-e+e-, PDG 2026) x (eps_dir / eps_X) x N_X / N_dir,
#         N_dir = a (C5) in the eta window (all M(e+e-)), eps_dir per direct decay
#   (eps per produced eta decay: expected counts / (sigma_SMASH(eta) x L x 63/64 x BR used in SMASH (PDG 2020)).)
# Run: python3 -I plots/all_report/prevstyle_results_xs.py   (a few seconds; one short ROOT call)
# ============================================================================
import json
import math
import subprocess

import numpy as np
from scipy import optimize, special

OUT = "/hdd3/przygoda/HADES/FAT_out/all_report/prevstyle/"
MACRO = "/home/przygoda/HADES/FAT/FAT/plots/all_report/prevstyle_xs_trolke.C"
L = 6079.2
PT3 = 64. / 63.
NB2RAW = L / PT3
SIG_SMASH = 1.0694                     # mb, SMASH sigma(pp -> eta X) of the simulated sample (report, v8)
N_ETA = SIG_SMASH * 1e6 * L / PT3      # SMASH eta in the data units (expected, x 63/64), for the efficiencies
CSEC_B, CSEC_P, ONE_M_FETA = 1.0104, 1.0116, 0.97
PDG_B = 1.00436                        # BR(eta -> pi+pi-pi0) 2026 / 2020 (C2, C3)
PDG_P = 1.0224                         # peak mix C4 + C5 in P, S0 (v8; x17p26_limit.C)
PDG_C4, PDG_C5 = 1.01422, 2.76 / 2.68  # single classes (report sec. P2)
BR3PI_26 = 0.2302
BR3PI_20 = BR3PI_26 / PDG_B
BRDIR_26, BRDIR_20 = 2.76e-4, 2.68e-4
K_B = SIG_SMASH * CSEC_B * ONE_M_FETA / PDG_B     # mb per unit bump scale
K_P = SIG_SMASH * CSEC_P * ONE_M_FETA / PDG_P     # mb per unit peak scale
NORM_UP_P = math.sqrt(6.9 ** 2 + 3.1 ** 2 + 3.15 ** 2 + 5.0 ** 2) / 100.
NORM_DN_P = math.sqrt(6.9 ** 2 + 4.1 ** 2 + 3.15 ** 2 + 5.0 ** 2) / 100.
NORM_UP_B = math.sqrt(6.9 ** 2 + 3.1 ** 2 + 3.17 ** 2 + 5.0 ** 2) / 100.
NORM_DN_B = math.sqrt(6.9 ** 2 + 4.1 ** 2 + 3.17 ** 2 + 5.0 ** 2) / 100.
PROD_UP, PROD_DN = 0.08 / 0.429, 0.03 / 0.429
REF_AP, REF_FU = 0.418249, 0.42574    # report: a_P (side bands) and f_eta(u-edge), C2 corrected
REF_SIG_P, REF_SIG_B = 0.429, 0.444   # report: sigma from the peak and from the u-edge [mb]
F_BUMP0, A_P0, G_ETA = 0.426, 0.418, 0.656
A_P_TOT = math.sqrt(0.035 ** 2 + 0.025 ** 2 + (0.115 * 0.418) ** 2)   # a_P: stat, MC, syst 11.5 % (Y_peak)
FROZEN_UL, TROLKE_REF, BES = 3.16e-5, 2.6e-4, 1.5e-5
SELS = ["all", "best"]
SELNAME = {"all": "isBest >= 0 (all stored combinations; nominal)", "best": "isBest == 1 (best combination only; earlier slides; comparison)"}
XW = [(0, "15.5-24.5"), (1, "12-22")]
LOG = open(OUT + "prevstyle_xs.log", "w")


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    LOG.write(s + "\n")
    LOG.flush()


DJ = json.load(open(OUT + "prevstyle_data.json"))
SJ = json.load(open(OUT + "prevstyle_sim.json"))
NUM = json.load(open(OUT + "prevstyle_numbers.json"))
TS = json.load(open(OUT + "prevstyle_twoscale_numbers.json"))
PARTS = {k: tuple(v) for k, v in DJ["meta"]["parts"].items()}
WIN = [w[0] for w in DJ["meta"]["windows"]]
iW = {w: i for i, w in enumerate(WIN)}
MEEWIN = DJ["meta"]["meewin"]           # SBlo, Wpk, SBhi, PB
CHK = []


def check(what, val, ref, tol):
    ok = abs(val - ref) <= tol
    CHK.append(f"{what:<92s} {val:15.7g} ref {ref:15.7g}  {'OK' if ok else 'DIFFERS'}")
    return ok


def part(arr, tag):
    o, n = PARTS[tag]
    return np.array(arr[o:o + n], dtype=float)


def D(s, tag, what="sig"):
    d = DJ["sel"][s]
    N, NV = part(d["N"], tag), part(d["NV"], tag)
    CB, CBe = part(d["CB"], tag), part(d["CBrep"], tag)
    if what == "N":
        return N, NV
    if what == "CB":
        return CB, CBe ** 2
    return N - CB, NV + CBe ** 2


def S(s, comps, tag):
    v = 0.
    e = 0.
    for c in comps:
        x = SJ["sel"][s][c]
        v = v + part(x["s"], tag) * NB2RAW
        e = e + part(x["v"], tag) * NB2RAW ** 2
    return v, e


GROUPS = {"bump": ["C2", "C3"], "peak": ["C4", "C5"], "Geta": ["Geta"], "nonEta": ["GD", "GC"]}
XL = 0.30 + 0.005 * np.arange(160)


def b5(a):
    return a[1:161]


def bins_in(lo, hi):
    return (XL >= lo - 1e-9) & (XL + 0.005 <= hi + 1e-9)


# =========================================================================== TASK A: sigma from the two-scale demonstration
# re-implementation of the nominal two-scale fit of prevstyle_twoscale.py (window 0.575-0.70, fit 0.43-0.575), to split the fit
# covariance into data and SMASH statistics; the nominal numbers must reproduce section 9 exactly (checks)
def two_scale_nominal(s):
    dv, dvar = (b5(x) for x in D(s, "b5"))
    comp = {g: tuple(b5(x) for x in S(s, cs, "b5")) for g, cs in GROUPS.items()}
    compw = {g: S(s, cs, "win") for g, cs in GROUPS.items()}
    dw, dwv = D(s, "win")
    ic = iW["CR"]

    def cb_win(fB, aP):
        eta = fB * compw["bump"][0][ic] + aP * compw["peak"][0][ic] + G_ETA * compw["Geta"][0][ic]
        etav = fB ** 2 * compw["bump"][1][ic] + aP ** 2 * compw["peak"][1][ic] + G_ETA ** 2 * compw["Geta"][1][ic]
        ne, nev = compw["nonEta"][0][ic], compw["nonEta"][1][ic]
        cb = (dw[ic] - eta) / ne
        return cb, math.sqrt(dwv[ic] + etav + cb ** 2 * nev) / ne

    def fit(cb):
        ne, nev = comp["nonEta"]
        ge, gev = comp["Geta"]
        B = cb * ne + G_ETA * ge
        Bv = cb ** 2 * nev + G_ETA ** 2 * gev
        T = [comp["bump"][0], comp["peak"][0]]
        Tv = [comp["bump"][1], comp["peak"][1]]
        y = dv - B
        m = bins_in(0.43, 0.575) & (dvar > 0)
        a = np.array([F_BUMP0, A_P0])
        for _ in range(200):
            var = dvar + Bv + a[0] ** 2 * Tv[0] + a[1] ** 2 * Tv[1]
            w = 1. / var[m]
            A = np.array([[np.sum(T[i][m] * T[j][m] * w) for j in range(2)] for i in range(2)])
            b = np.array([np.sum(y[m] * T[i][m] * w) for i in range(2)])
            an = np.linalg.solve(A, b)
            if np.max(np.abs(an - a)) < 1e-11:
                a = an
                break
            a = an
        var = dvar + Bv + a[0] ** 2 * Tv[0] + a[1] ** 2 * Tv[1]
        w = 1. / var[m]
        A = np.array([[np.sum(T[i][m] * T[j][m] * w) for j in range(2)] for i in range(2)])
        cov = np.linalg.inv(A)
        # sandwich split: cov = A^-1 [T' W (V_data + V_SMASH) W T] A^-1
        vd = dvar[m]
        vmc = (Bv + a[0] ** 2 * Tv[0] + a[1] ** 2 * Tv[1])[m]
        Md = np.array([[np.sum(T[i][m] * T[j][m] * w * w * vd) for j in range(2)] for i in range(2)])
        Mm = np.array([[np.sum(T[i][m] * T[j][m] * w * w * vmc) for j in range(2)] for i in range(2)])
        cov_d = cov @ Md @ cov
        cov_m = cov @ Mm @ cov
        return dict(a=a, cov=cov, cov_d=cov_d, cov_m=cov_m)

    cb0, _ = cb_win(F_BUMP0, A_P0)
    f0 = fit(cb0)
    cb1, cbe1 = cb_win(*f0["a"])
    f1 = fit(cb1)
    fp, fm = fit(cb1 + cbe1), fit(cb1 - cbe1)
    ecb = 0.5 * np.abs(fp["a"] - fm["a"])
    return dict(cb=cb1, cbe=cbe1, a=f1["a"], e=np.sqrt(np.diag(f1["cov"])), e_data=np.sqrt(np.diag(f1["cov_d"])),
                e_mc=np.sqrt(np.diag(f1["cov_m"])), ecb=ecb, rho=float(f1["cov"][0, 1] / math.sqrt(f1["cov"][0, 0] * f1["cov"][1, 1])),
                split_sum=np.diag(f1["cov_d"] + f1["cov_m"]), cov=np.diag(f1["cov"]))


def q(*x):
    return math.sqrt(sum(v * v for v in x))


A = {}
for s in SELS:
    r = two_scale_nominal(s)
    tsn = TS[s]["nominal"]
    for i, nm in enumerate(["bump", "peak"]):
        check(f"{s}: two-scale {nm} scale reproduces section 9", float(r["a"][i]), tsn["f1"]["a"][i], 1e-9)
        check(f"{s}: two-scale {nm} fit error reproduces section 9", float(r["e"][i]), tsn["f1"]["e"][i], 1e-9)
        check(f"{s}: two-scale {nm} c_b error reproduces section 9", float(r["ecb"][i]), tsn["ecb"][i], 1e-9)
        check(f"{s}: {nm}: data part^2 + SMASH part^2 = fit variance (sandwich split)", float(r["split_sum"][i]), float(r["cov"][i]),
              1e-12)
    check(f"{s}: c_b after the iteration reproduces section 9", float(r["cb"]), tsn["cb1"], 1e-9)
    # method systematic: largest up and down deviation in each family (normalisation window, fit range, background tilt), quadrature
    fam = {"window": ["norm 0.575-0.65", "norm 0.60-0.70", "norm 0.575-0.75"],
           "fit range": ["fit 0.40-0.575", "fit 0.45-0.575", "fit 0.47-0.575", "fit 0.43-0.56"],
           "tilt": ["tilt +3%/0.1", "tilt -3%/0.1"]}
    meth = {}
    for i, nm in enumerate(["bump", "peak"]):
        a0 = float(r["a"][i])
        famd = {}
        for f, vs in fam.items():
            dev = [TS[s][v]["f1"]["a"][i] - a0 for v in vs]
            famd[f] = (max(0., max(dev)), max(0., -min(dev)), {v: TS[s][v]["f1"]["a"][i] for v in vs})
        up = q(*[famd[f][0] for f in fam])
        dn = q(*[famd[f][1] for f in fam])
        upwf = q(famd["window"][0], famd["fit range"][0])
        dnwf = q(famd["window"][1], famd["fit range"][1])
        meth[nm] = dict(up=up, dn=dn, fam=famd, up_no_tilt=upwf, dn_no_tilt=dnwf)
    r["meth"] = meth
    # sigma
    sig = {}
    for i, nm in enumerate(["bump", "peak"]):
        K = K_B if nm == "bump" else K_P
        nu, nd = (NORM_UP_B, NORM_DN_B) if nm == "bump" else (NORM_UP_P, NORM_DN_P)
        a0 = float(r["a"][i])
        x = a0 * K
        st, mc, ecb = float(r["e_data"][i]) * K, float(r["e_mc"][i]) * K, float(r["ecb"][i]) * K
        mu, md = meth[nm]["up"] * K, meth[nm]["dn"] * K
        tot_up = q(st, mc, ecb, mu, x * nu)
        tot_dn = q(st, mc, ecb, md, x * nd)
        tot_nometh = q(st, mc, ecb, x * nu), q(st, mc, ecb, x * nd)
        ref = REF_SIG_B if nm == "bump" else REF_SIG_P
        sig[nm] = dict(scale=a0, K=K, sigma=x, stat=st, mc=mc, cb=ecb, meth_up=mu, meth_dn=md, norm_up=x * nu, norm_dn=x * nd,
                       prod_up=x * PROD_UP, prod_dn=x * PROD_DN, tot_up=tot_up, tot_dn=tot_dn, tot_nometh=tot_nometh,
                       fit_tot=q(st, mc, ecb),
                       ref=ref, diff=x - ref, ratio=x / ref,
                       ref_other=(REF_SIG_P if nm == "bump" else REF_SIG_B))
    r["sig"] = sig
    A[s] = r

# chain checks against the report
check("peak chain: a_P 0.418249 -> 0.429 mb (report)", REF_AP * K_P, 0.42927, 5e-5)
check("bump chain: f(u-edge) 0.42574 -> 0.444 mb (report)", REF_FU * K_B, 0.44428, 5e-5)
check("peak chain: 0.418 -> 0.429 mb (rounded inputs)", round(0.418 * K_P, 3), 0.429, 1e-9)
check("normalisation error of the headline: +9.6 % / -10.0 % (0.429 +0.041/-0.043)", NORM_UP_P * 100, 9.60, 0.01)
check("normalisation error of the headline (down)", NORM_DN_P * 100, 9.97, 0.01)
# PDG mix of the peak in this selection (for information): C4, C5 in P
mixinfo = {}
for s in SELS:
    c4 = S(s, ["C4"], "win")[0][iW["P"]]
    c5 = S(s, ["C5"], "win")[0][iW["P"]]
    mixinfo[s] = dict(C4=float(c4), C5=float(c5), factor=float((c4 * PDG_C4 + c5 * PDG_C5) / (c4 + c5)))

P("==== TASK A: sigma(pp -> eta X) from the two-scale demonstration (DEMONSTRATION, not a measurement) ====")
P(f"chain: bump sigma = f x {K_B:.6f} mb (1.0694 x 1.0104 x 0.97 / 1.00436); peak sigma = a x {K_P:.6f} mb (1.0694 x 1.0116 x 0.97 / 1.0224)")
P(f"  report: 0.418249 x K_P = {REF_AP * K_P:.5f} mb (0.429); 0.42574 x K_B = {REF_FU * K_B:.5f} mb (0.444)")
P(f"  normalisation (headline): peak +{NORM_UP_P * 100:.2f} / -{NORM_DN_P * 100:.2f} %, bump +{NORM_UP_B * 100:.2f} / -{NORM_DN_B * 100:.2f} %; "
  f"production model +{PROD_UP * 100:.1f} / -{PROD_DN * 100:.1f} % (separate)")
for s in SELS:
    P(f"  PDG 2026 factor of the peak mix in P for this selection ({s}): C4 {mixinfo[s]['C4']:.1f}, C5 {mixinfo[s]['C5']:.1f} -> "
      f"{mixinfo[s]['factor']:.5f} (report, S0: 1.0224; difference {100 * (1.0224 / mixinfo[s]['factor'] - 1):+.2f} % on sigma; not applied)")
for s in SELS:
    r = A[s]
    P(f"\n[{SELNAME[s]}] c_b = {r['cb']:.4f} +- {r['cbe']:.4f}; correlation bump-peak {r['rho']:+.3f}")
    for i, nm in enumerate(["bump", "peak"]):
        x = r["sig"][nm]
        m = r["meth"][nm]
        P(f"  {nm}: scale {r['a'][i]:.4f} +- {r['e'][i]:.4f} (fit = data {r['e_data'][i]:.4f} (+) SMASH {r['e_mc'][i]:.4f}) +- {r['ecb'][i]:.4f} (c_b) "
          f"+{m['up']:.4f} -{m['dn']:.4f} (method)")
        for f, (u, dd, vals) in m["fam"].items():
            P(f"      method family {f:10s}: +{u:.4f} / -{dd:.4f}  " + ", ".join(f"{k} {v:.4f}" for k, v in vals.items()))
        P(f"      method without the tilt: +{m['up_no_tilt']:.4f} / -{m['dn_no_tilt']:.4f}")
        P(f"    sigma = {x['sigma']:.4f} +- {x['stat']:.4f} (data) +- {x['mc']:.4f} (SMASH) +- {x['cb']:.4f} (c_b) +{x['meth_up']:.4f} -{x['meth_dn']:.4f} (method) "
          f"+{x['norm_up']:.4f} -{x['norm_dn']:.4f} (norm.) mb; production model +{x['prod_up']:.3f} -{x['prod_dn']:.3f} mb separate")
        P(f"    total (measurement, without the model item): +{x['tot_up']:.4f} -{x['tot_dn']:.4f} mb; without the method: +{x['tot_nometh'][0]:.4f} "
          f"-{x['tot_nometh'][1]:.4f}; fit + c_b only: {x['fit_tot']:.4f}")
        P(f"    against the report's {x['ref']:.3f} mb: {x['diff']:+.4f} mb (ratio {x['ratio']:.3f}); against {x['ref_other']:.3f} mb: "
          f"{x['sigma'] - x['ref_other']:+.4f} mb")

# =========================================================================== TASK B: X17 in the earlier style
KVAR = {"all": [("global fit", NUM["all"]["sb"]["k"]), ("smoothed fit", NUM["all"]["sb"]["k_sm"]), ("width 27/24", 27. / 24.)],
        "best": [("global fit", NUM["best"]["sb"]["k"]), ("smoothed fit", NUM["best"]["sb"]["k_sm"]), ("width 27/24", 27. / 24.)]}
KERR = {s: NUM[s]["sb"]["ke"] for s in SELS}
AVAR = {}
for s in SELS:
    ts = A[s]
    peak_ts_err = q(float(ts["e"][1]), float(ts["ecb"][1]))
    AVAR[s] = [("side bands 0.418", A_P0, A_P_TOT),
               ("two-scale", float(ts["a"][1]), peak_ts_err),
               ("global fit", NUM[s]["gf_raw"]["eta"][0], NUM[s]["gf_raw"]["etae"][0])]
REFK, REFA = "global fit", "two-scale"      # reference configuration of this demonstration


def xr(s, qi, w, src="data"):
    i = qi * 2 + w
    if src == "data":
        N, NV = (x[i] for x in D(s, "xr", "N"))
        CB, CBv = (x[i] for x in D(s, "xr", "CB"))
        return float(N), float(NV), float(CB), float(CBv)
    v, e = S(s, src, "xr")
    return float(v[i]), float(e[i])


def gh_nodes(n=64):
    z, w = np.polynomial.hermite_e.hermegauss(n)
    return z, w / math.sqrt(2 * math.pi)


GHZ, GHW = gh_nodes()


def p_le(x, mu):
    """P(n <= x | mu), Poisson"""
    return special.gammaincc(x + 1, np.maximum(mu, 1e-300))


def cls_ul(x, B, sB, cl=0.90):
    """CLs upper limit for a Poisson count x with a Gaussian background B +- sB (marginalised; truncated at b > 0)"""
    if sB > 0:
        b = B + sB * GHZ
        w = GHW * (b > 0)
        w = w / w.sum()
        b = np.maximum(b, 0.)
    else:
        b, w = np.array([B]), np.array([1.])
    clb = float(np.sum(w * p_le(x, b)))

    def f(sv):
        return float(np.sum(w * p_le(x, b + sv))) / clb - (1. - cl)

    hi = max(10., 3. * (abs(x - B) + 5 * math.sqrt(x + sB ** 2 + 1)))
    return optimize.brentq(f, 0., hi, xtol=1e-6), clb


check("CLs code: x = 0, B = 0 -> -ln(0.1) = 2.302585", cls_ul(0, 0., 0.)[0], -math.log(0.1), 1e-5)

CONF = []          # one entry per configuration
BAL = {}
for s in SELS:
    fb = float(A[s]["a"][0])
    BW = S(s, ["C2", "C3"], "win")
    ND_win = BW[0][iW["B"]]                                    # SMASH C2 + C3 in the bump window B (expected counts)
    C5W = S(s, ["C5"], "win")[0][iW["Wpk"]]                    # SMASH C5 in the eta window, all M(ee)
    for w, wn in XW:
        NW, NVW, CBW, CBvW = xr(s, 1, w)
        Nl, NVl, CBl, CBvl = xr(s, 0, w)
        Nh, NVh, CBh, CBvh = xr(s, 2, w)
        SB = (Nl - CBl) + (Nh - CBh)
        SBv = NVl + NVh + CBvl + CBvh
        pk, pkv = xr(s, 1, w, ["C4", "C5"])
        c6w, c6wv = xr(s, 1, w, ["C6"])
        c6sb = xr(s, 0, w, ["C6"])[0] + xr(s, 2, w, ["C6"])[0]
        # extra (diagnostic, not in the requested ranges): for each k the "balanced" eta scale that makes the window minus k x side
        # bands equal to the SMASH eta peak over all M(ee): a_bal = (W - k SB)_all / (C4 + C5)_W,all
        Wall, Wall_e = NUM[s]["sb"]["W"], NUM[s]["sb"]["We"]
        SBall, SBall_e = NUM[s]["sb"]["SB"], NUM[s]["sb"]["SBe"]
        pka, pkav = S(s, ["C4", "C5"], "win")
        pka, pkav = float(pka[iW["Wpk"]]), float(pkav[iW["Wpk"]])
        for kn, k in KVAR[s]:
            abal = (Wall - k * SBall) / pka
            abal_e = q(math.sqrt(Wall_e ** 2 + k ** 2 * SBall_e ** 2) / pka, abal * math.sqrt(pkav) / pka)
            BAL[(s, kn)] = (abal, abal_e)
            for an, a, ae in AVAR[s] + [("balanced", abal, abal_e)]:
                B = CBW + k * SB + a * pk
                comps = dict(cb=math.sqrt(CBvW), sb=k * math.sqrt(SBv), k=SB * KERR[s], eta_scale=ae * pk, eta_mc=a * math.sqrt(pkv))
                sB = q(*comps.values())
                c6eff = c6w - k * c6sb                        # X17 counts in the residual per unit BR (side-band leakage subtracted)
                CONF.append(dict(sel=s, win=wn, k_name=kn, k=k, a_name=an, a=a, a_err=ae, x=int(round(NW)), NV=NVW, CB=CBW,
                                 SB=SB, SBv=SBv, peak=pk, B=B, sB=sB, sB_comps=comps, resid=NW - B, c6w=c6w, c6sb=c6sb, c6eff=c6eff,
                                 c6rel=math.sqrt(c6wv) / c6w, fb=fb, ND_win=ND_win, C5W=C5W, extra=(an == "balanced"),
                                 id=f"{s}|{wn}|{kn}|{an}".replace(" ", "_")))

# section-6 residuals reproduced (all M(ee) conventions of prevstyle_ana.py)
for s in SELS:
    for w, wn in XW:
        x6 = NUM[s]["x17"][wn]
        c = [c for c in CONF if c["sel"] == s and c["win"] == wn]
        g = next(c for c in c if c["k_name"] == "global fit" and c["a_name"] == "global fit")
        check(f"{s} {wn}: residual (k global fit, a global fit) = section 6 'excess'", g["resid"], x6["excess"], 1e-6)
        g = next(c for c in c if c["k_name"] == "width 27/24" and c["a_name"] == "global fit")
        check(f"{s} {wn}: residual (k = 27/24) = section 6", g["resid"], x6["excess_kw"], 1e-6)
        g = next(c for c in c if c["k_name"] == "smoothed fit" and c["a_name"] == "global fit")
        check(f"{s} {wn}: residual (smoothed k) = section 6", g["resid"], x6["excess_ksm"], 1e-6)
        g = next(c for c in c if c["k_name"] == "global fit" and c["a_name"] == "side bands 0.418")
        check(f"{s} {wn}: residual (a_P = 0.418) = section 6", g["resid"], x6["excess_aP"], 1e-6)

# ---- TRolke (one ROOT call)
with open(OUT + "prevstyle_xs_trolke_in.txt", "w") as f:
    f.write("# id x B sigma_B  (prevstyle_results_xs.py; TRolke model 5)\n")
    f.write("frozen_X17c_check 6424 6311.7 1195.5\n")
    for c in CONF:
        f.write(f"{c['id']} {c['x']} {c['B']:.6f} {c['sB']:.6f}\n")
res = subprocess.run(["root", "-l", "-b", "-q", f'{MACRO}("{OUT}prevstyle_xs_trolke_in.txt")'], capture_output=True, text=True,
                     cwd="/home/przygoda/HADES/FAT/FAT")
with open(OUT + "prevstyle_xs_trolke_out.txt", "w") as f:
    f.write(res.stdout)
TR = {}
for line in res.stdout.splitlines():
    p = line.split()
    if len(p) == 4 and p[1] in ("0", "1"):
        TR[(p[0], int(p[1]))] = (float(p[2]), float(p[3]))
check("TRolke call reproduces the frozen one-window count (N_X < 2083 counts, trolke_x17.txt)", TR[("frozen_X17c_check", 1)][1], 2083.1456, 0.01)

for c in CONF:
    lo, hi = TR[(c["id"], 1)]
    lo0, hi0 = TR[(c["id"], 0)]
    c["tr_lo"], c["tr_ul"], c["tr_ul_unb"], c["tr_lo_unb"] = lo, hi, hi0, lo0
    c["cls_ul"], c["clb"] = cls_ul(c["x"], c["B"], c["sB"])
    c["cls_ul_nob"] = cls_ul(c["x"], c["B"], 0.)[0]
    # conversions
    epsX = c["c6eff"] / N_ETA
    epsD = c["ND_win"] / (N_ETA * BR3PI_20)
    ND = c["fb"] * c["ND_win"]
    epsdir = c["C5W"] / (N_ETA * BRDIR_20)
    Ndir = c["a"] * c["C5W"]
    c.update(epsX=epsX, epsD=epsD, ND=ND, epsdir=epsdir, Ndir=Ndir)
    for m in ("tr", "cls"):
        ul = c[f"{m}_ul"]
        c[f"{m}_BRa"] = BR3PI_26 * (epsD / epsX) * ul / ND
        c[f"{m}_BRb"] = BRDIR_26 * (epsdir / epsX) * ul / Ndir
    c["tr_BRa_noleak"] = c["tr_BRa"] * c["c6eff"] / c["c6w"]
    # shortcut forms (check of the epsilon chain)
    c["chk_a"] = PDG_B * c["tr_ul"] / (c["fb"] * c["c6eff"])
    c["chk_b"] = PDG_C5 * c["tr_ul"] / (c["a"] * c["c6eff"])
cc = CONF[0]
check("epsilon chain (a) = 1.00436 N_X / (f_b C6eff)", cc["tr_BRa"], cc["chk_a"], 1e-12 * cc["chk_a"] + 1e-20)
check("epsilon chain (b) = (2.76/2.68) N_X / (a C6eff)", cc["tr_BRb"], cc["chk_b"], 1e-12 * cc["chk_b"] + 1e-20)

P("\n==== TASK B: X17 upper limits from a one-window count in the eta window (DEMONSTRATION, not a measurement) ====")
P("x = data N in the window cell; B = CB_W + k (N - CB)_SB + a (C4 + C5)_W; residual = x - B; limits on the X17 counts in the residual")
P(f"N_eta(SMASH, data units) = {N_ETA:.5e}; BR(eta -> pi+pi-pi0) 2026 {BR3PI_26} (SMASH 2020 {BR3PI_20:.5f}); BR(eta -> pi+pi-e+e-) 2026 {BRDIR_26} "
  f"(SMASH 2020 {BRDIR_20})")
for s in SELS:
    P(f"\n[{SELNAME[s]}]  f_b (two-scale bump) = {A[s]['a'][0]:.4f}; k variants " + ", ".join(f"{kn} {k:.4f}" for kn, k in KVAR[s])
      + f" (sigma_k {KERR[s]:.4f}); a variants " + ", ".join(f"{an} {a:.4f} +- {ae:.4f}" for an, a, ae in AVAR[s]))
    for w, wn in XW:
        cs = [c for c in CONF if c["sel"] == s and c["win"] == wn]
        c0 = cs[0]
        P(f"  M(ee) {wn} MeV/c2: x = N_W = {c0['x']}, NV/N = {c0['NV'] / c0['x']:.3f}, CB_W = {c0['CB']:.1f}, side bands N - CB = {c0['SB']:.1f} "
          f"+- {math.sqrt(c0['SBv']):.1f}, SMASH (C4 + C5)_W = {c0['peak']:.1f}; C6 per unit BR: W {c0['c6w']:.4e}, SB {c0['c6sb']:.4e} "
          f"(rel. MC {c0['c6rel']:.4f})")
        P(f"    N_D = f_b (C2 + C3)_B = {c0['ND']:.1f} (SMASH {c0['ND_win']:.1f}); eps_D = {c0['epsD']:.5e} per eta -> pi+pi-pi0; "
          f"SMASH C5 in W = {c0['C5W']:.1f}; eps_dir = {c0['epsdir']:.5e} per direct decay")
        P(f"    {'k':12s} {'a':17s} {'B':>9s} {'sB':>7s} {'resid':>8s} {'eps_X':>10s} | {'TRolke UL':>9s} {'(lo)':>6s} {'CLs UL':>8s} {'CLb':>5s} | "
          f"{'BR(a) TR':>9s} {'BR(a) CLs':>9s} {'BR(b) TR':>9s} {'BR(b) CLs':>9s} [1e-5]")
        for c in cs:
            P(f"    {c['k_name']:12s} {c['a_name']:17s} {c['B']:9.1f} {c['sB']:7.1f} {c['resid']:8.1f} {c['epsX']:10.4e} | {c['tr_ul']:9.1f} "
              f"{c['tr_lo']:6.0f} {c['cls_ul']:8.1f} {c['clb']:5.3f} | {1e5 * c['tr_BRa']:9.2f} {1e5 * c['cls_BRa']:9.2f} {1e5 * c['tr_BRb']:9.2f} "
              f"{1e5 * c['cls_BRb']:9.2f}")
        c = next(c for c in cs if c["k_name"] == REFK and c["a_name"] == REFA)
        P(f"    sigma_B parts (reference k {REFK}, a {REFA}): " + ", ".join(f"{k} {v:.1f}" for k, v in c["sB_comps"].items()))

# summaries: reference and spread over (k, a) per window, and over both windows
SUMM = {}
for s in SELS:
    SUMM[s] = {}
    for key in ["tr_ul", "cls_ul", "tr_BRa", "cls_BRa", "tr_BRb", "cls_BRb", "resid"]:
        SUMM[s][key] = {}
        for w, wn in XW:
            cs = [c for c in CONF if c["sel"] == s and c["win"] == wn and not c["extra"]]
            ref = next(c for c in cs if c["k_name"] == REFK and c["a_name"] == REFA)
            vals = [c[key] for c in cs]
            bal = [c[key] for c in CONF if c["sel"] == s and c["win"] == wn and c["extra"]]
            SUMM[s][key][wn] = dict(ref=ref[key], lo=min(vals), hi=max(vals), bal_lo=min(bal), bal_hi=max(bal))
        allv = [c[key] for c in CONF if c["sel"] == s and not c["extra"]]
        SUMM[s][key]["both"] = dict(lo=min(allv), hi=max(allv))
P("\n==== summary: reference (k global fit, a two-scale) and range over k x a (9 combinations) ====")
for s in SELS:
    P(f"[{SELNAME[s]}]")
    for key, lab, sc in [("resid", "residual [counts]", 1), ("tr_ul", "TRolke UL [counts]", 1), ("cls_ul", "CLs UL [counts]", 1),
                         ("tr_BRa", "TRolke BR (a) bump [1e-5]", 1e5), ("cls_BRa", "CLs BR (a) bump [1e-5]", 1e5),
                         ("tr_BRb", "TRolke BR (b) direct [1e-5]", 1e5), ("cls_BRb", "CLs BR (b) direct [1e-5]", 1e5)]:
        P(f"  {lab:30s} " + "; ".join(f"{wn}: {sc * SUMM[s][key][wn]['ref']:.2f} [{sc * SUMM[s][key][wn]['lo']:.2f}, {sc * SUMM[s][key][wn]['hi']:.2f}]"
                                       for _, wn in XW) + f"; both windows [{sc * SUMM[s][key]['both']['lo']:.2f}, {sc * SUMM[s][key]['both']['hi']:.2f}]")
P("balanced eta scales (window minus k x side bands = SMASH eta peak over all M(ee)): "
  + "; ".join(f"{s}: " + ", ".join(f"k {kn} -> a {BAL[(s, kn)][0]:.3f} +- {BAL[(s, kn)][1]:.3f}" for kn, _ in KVAR[s]) for s in SELS))
for s in SELS:
    for key, lab, sc in [("resid", "residual", 1), ("tr_BRa", "TRolke BR (a)", 1e5), ("cls_BRa", "CLs BR (a)", 1e5), ("tr_BRb", "TRolke BR (b)", 1e5),
                         ("cls_BRb", "CLs BR (b)", 1e5)]:
        P(f"  balanced pairs, {s}, {lab:14s}: " + "; ".join(f"{wn}: {sc * SUMM[s][key][wn]['bal_lo']:.2f}-{sc * SUMM[s][key][wn]['bal_hi']:.2f}" for _, wn in XW))
# ratios to the frozen limit, the one-window TRolke of the report and BESIII
P("\nreference / frozen 3.16e-5 (isBest >= 0, 15.5-24.5): "
  + ", ".join(f"{k} {SUMM['all'][k]['15.5-24.5']['ref'] / FROZEN_UL:.2f}" for k in ["tr_BRa", "cls_BRa", "tr_BRb", "cls_BRb"])
  + "; range over everything (all, both windows, both methods, both normalisations): "
  + f"{min(SUMM['all'][k]['both']['lo'] for k in ['tr_BRa', 'cls_BRa', 'tr_BRb', 'cls_BRb']) * 1e5:.2f}-"
  + f"{max(SUMM['all'][k]['both']['hi'] for k in ['tr_BRa', 'cls_BRa', 'tr_BRb', 'cls_BRb']) * 1e5:.2f} x 1e-5")
# side-band leakage of X17 and the plain k x sigma_B view
c = next(c for c in CONF if c["sel"] == "all" and c["win"] == "15.5-24.5" and c["k_name"] == REFK and c["a_name"] == REFA)
P(f"X17 leakage into the side bands (reference): C6eff / C6_W = {c['c6eff'] / c['c6w']:.4f} (without the leakage the BR limits would be "
  f"x {c['c6eff'] / c['c6w']:.3f})")
NX_FROZEN = {"a": FROZEN_UL * c["fb"] / PDG_B * c["c6eff"], "b": FROZEN_UL * c["a"] / PDG_C5 * c["c6eff"]}
P(f"an X17 at the frozen limit 3.16e-5 would give {NX_FROZEN['a']:.0f} counts (eta scale of normalisation (a)) or {NX_FROZEN['b']:.0f} counts "
  f"(normalisation (b), a two-scale) in the residual (reference, isBest >= 0, 15.5-24.5); B = {c['B']:.0f} +- {c['sB']:.0f}, x = {c['x']}; "
  f"background share of x {c['B'] / c['x']:.3f}")
kspread = [k for _, k in KVAR["all"]]
P(f"k spread 1.064 -> 1.125 moves B by {(max(kspread) - min(kspread)) * c['SB']:.0f} counts = {(max(kspread) - min(kspread)) * c['SB'] / c['sB']:.1f} sigma_B "
  f"(sigma_B {c['sB']:.0f}); a spread 0.418 -> 0.697 moves it by {(0.6973 - 0.418) * c['peak']:.0f} counts")

for x in CHK:
    P(x)
NBAD = sum("DIFFERS" in x for x in CHK)
P(f"checks: {len(CHK)} done, {NBAD} differ")


# ---------------------------------------------------------------- numbers
def clean(o):
    if isinstance(o, dict):
        return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, np.ndarray):
        return clean(o.tolist())
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    return o


OUTJ = dict(note="DEMONSTRATION in the style of the earlier analysis; not a measurement; no result of the report uses it",
            constants=dict(K_B=K_B, K_P=K_P, SIG_SMASH=SIG_SMASH, CSEC_B=CSEC_B, CSEC_P=CSEC_P, ONE_M_FETA=ONE_M_FETA, PDG_B=PDG_B, PDG_P=PDG_P,
                           NORM_UP_P=NORM_UP_P, NORM_DN_P=NORM_DN_P, NORM_UP_B=NORM_UP_B, NORM_DN_B=NORM_DN_B, PROD_UP=PROD_UP, PROD_DN=PROD_DN,
                           N_ETA=N_ETA, BR3PI_26=BR3PI_26, BRDIR_26=BRDIR_26, A_P_TOT=A_P_TOT, REFK=REFK, REFA=REFA),
            taskA={s: dict(cb=A[s]["cb"], cbe=A[s]["cbe"], a=A[s]["a"], e=A[s]["e"], e_data=A[s]["e_data"], e_mc=A[s]["e_mc"], ecb=A[s]["ecb"],
                           rho=A[s]["rho"], meth=A[s]["meth"], sig=A[s]["sig"], pdg_mix_P=mixinfo[s]) for s in SELS},
            taskB=dict(kvar=KVAR, kerr=KERR, avar=AVAR, balanced={f"{s}|{kn}": v for (s, kn), v in BAL.items()}, conf=CONF, summary=SUMM),
            nx_frozen_ref=NX_FROZEN, checks=CHK)
json.dump(clean(OUTJ), open(OUT + "prevstyle_xs_numbers.json", "w"), indent=1)
P("written prevstyle_xs_numbers.json")
