#!/usr/bin/env python3
# ============================================================================
# udatasim_overlay.py  --  package udatasim: data against SMASH in the edge variable u = M4 - E(M(pi+pi-), M(e+e-))
# ============================================================================
# Written 2026-11-15 (prefix udatasim). New file. edgeu_core.py and peaktail_core.py are imported and only read; no edgeu,
# c2recon, peaktail or wp3 output is written. Outputs go to /hdd3/przygoda/HADES/FAT_out/all_report/udatasim/.
#
# What it does:
#   1. rebuilds the nominal edgeu fit with the frozen edgeu code (edgeu_core: u, CB rule, step template, iterated GLS; the
#      configuration of edgeu_ana.py NOM, as c2recon_edgeu.py) and checks it against edgeu_ana.json;
#   2. joins the per-count QED weights of D38 (c2recon_weights.npz) to the simulated C2 counts and recomputes
#      f_eta(u-edge) = h_data / [h(C2 with D38 weights) + h(C3)] with the same linear estimator;
#   3. builds two SMASH models in u, per class C2 (D38 weights), C3, C4, C5, G_eta, G_D, G_C:
#        (a) SMASH as it is (sigma_prod x PDG BR x V3, no scale);
#        (b) C2, C3 x f_eta(u-edge, D38); C4, C5 x a_P = 0.418249; G_eta x a_eta g = 0.656 (Part III, peak side bands);
#            G_D, G_C x the frozen stage-A normalisation of Fig. 12 / Table C2 (stageA.txt): c_j (1 + b x) [G_D + kappa_B G_C],
#            x = (M4' - 0.55)/0.1, j = the M(pipi) slice of defs.h (S1 < 0.34, S2 0.34-0.37, S3 0.37-0.39, S4 0.39-0.43, S5);
#   4. data N - CB (CB = 2 sqrt(N++ N--) of the data, per M(pipi) slice in 10 MeV/c2 cells of u, the edgeu rule) for
#      S0, 0.30 <= M(pipi) < 0.40 (the edgeu fit region) and for all M(pipi) where the edge exists (0.27-0.415, u defined
#      below 0.4116);
#   5. covariances: data N from the event groups, CB from the like-sign counts (edgeu rule), simulation MC from the
#      simulated event groups (every scaled class together), the stage-A normalisation from its parameter covariance
#      (cov_total of stageA.txt); the D40 peak-tail band (t = 0.5 / 1.5) from the peaktail tail flag (peaktail dump);
#   6. chi2 of data against (a) and (b) in -45...+45 and +20...+60 MeV/c2, data/model in three u bands, per M(ee) class;
#      a check of the stage-A implementation in CR-U (M4 0.575-0.70, E2), where Fig. 12 gives 0.999.
# Units: counts = stored pi+pi-e+e- combinations; data raw counts per 2.5 MeV/c2; simulation = expected counts
# (sum w x 6079.2/221.05 x 63/64, already in the edgeu columns); both sides of every ratio in the same unit.
# Outputs: udatasim_overlay.json, udatasim_overlay.log, udatasim_figdata.npz
# Run: python3 -I plots/all_report/udatasim_overlay.py     (about 2 min, 6 GB)
# ============================================================================
import json
import math
import os
import sys
import time

import numpy as np
import scipy.sparse as sps
from scipy.optimize import curve_fit

sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import edgeu_core as C
import peaktail_core as PT

O = "/hdd3/przygoda/HADES/FAT_out/all_report/udatasim/"
os.makedirs(O, exist_ok=True)
T0 = time.time()
LOG = open(O + "udatasim_overlay.log", "w")


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s); LOG.write(s + "\n"); LOG.flush()


RES = {}
FD = {}
AG = 0.656            # a_eta g: the effective G_eta scale of Part III (peak side bands; F-main a 0.382 x g 1.72)
TBAND = (0.5, 1.5)    # D40: peak-tail family t = 1 +- 0.5
UPLOT = (-0.100, 0.100)
RANGES = {"fit": (-0.045, 0.045), "peak": (0.020, 0.060)}
BANDS = {"below": (-0.045, -0.005), "edge": (-0.005, 0.005), "above": (0.005, 0.045)}

# ---------------------------------------------------------------- stage A (frozen; Fig. 12 = fix_C3, Table C2)
def read_stageA(path="/hdd3/przygoda/HADES/FAT_out/all_report/wp3/stageA.txt"):
    L = open(path).read().splitlines()
    out = {}
    for i, l in enumerate(L):
        t = l.split()
        if not t:
            continue
        if t[0] == "names":
            out["names"] = t[1:]
        if t[0] == "values":
            out["values"] = np.array([float(v) for v in t[1:]])
        if t[0] == "err_total":
            out["err_total"] = np.array([float(v) for v in t[1:]])
        if t[0] == "cov_total":
            n = int(t[1])
            out["cov_total"] = np.array([[float(v) for v in L[i + 1 + r].split()] for r in range(n)])
    return out
SA = read_stageA()
assert SA["names"][:7] == ["c_1", "c_2", "c_3", "c_4", "c_5", "b", "kappa_B"], SA["names"]
CJ = SA["values"][:5]; BT = SA["values"][5]; KB = SA["values"][6]; KE2 = SA["values"][8]
COVA = SA["cov_total"][:7, :7]
P("# udatasim_overlay.py, 2026-11-15")
P("   stage A (stageA.txt values): c_j = %s, b = %.6f, kappa_B = %.5f, k_E2 = %.5f (Table C2: 0.9878 1.0012 0.9751 1.0552 0.9091, -0.0683, 3.0627, 0.9543)" %
  (" ".join("%.5f" % c for c in CJ), BT, KB, KE2))
RES["stageA"] = dict(c=CJ.tolist(), b=BT, kB=KB, kE2=KE2, err=SA["err_total"][:7].tolist())


def slice_defs(m):
    """defs.h Slice(): S1 < 0.34, S2 0.34-0.37, S3 0.37-0.39, S4 0.39-0.43, S5 >= 0.43 (data scale)"""
    j = np.zeros(m.shape, dtype=np.int64)
    for e in (0.34, 0.37, 0.39, 0.43):
        j += (m >= e)
    return j


# ---------------------------------------------------------------- columns
D = C.load_data()
S = C.load_sim()
gS = S["g"]
P("   data counts %d, simulated counts %d (%.0f s)" % (D["n"], S["n"], time.time() - T0))
w0 = C.sim_weight(S, c2corr=False)          # S0, mode A, w x V3 (code 108), expected counts, no C2 factor
us, mpp, m4p = C.u_sim(S)

# ---------------------------------------------------------------- 1. the nominal edgeu fit (edgeu_ana.py NOM, as c2recon_edgeu.py)
ut_s = C.EDGE_SCALE * (S["m4t"].astype(np.float64) - C.m4max_true(S["mpipit"].astype(np.float64), S["meet"].astype(np.float64)))
def reso_sample(near=-0.015):
    m = (gS <= 1) & np.isfinite(us) & np.isfinite(ut_s) & (w0 > 0) & (mpp >= 0.30) & (mpp < 0.40) & (ut_s > near)
    return (us - ut_s)[m], w0[m]
def gauss_core(d, ww):
    h, e = np.histogram(d, bins=np.arange(-0.040, 0.0401, 0.002), weights=ww)
    h2, _ = np.histogram(d, bins=e, weights=ww * ww)
    x = 0.5 * (e[1:] + e[:-1])
    f = lambda x, a, mu, s: a * np.exp(-0.5 * (x - mu) ** 2 / s ** 2)
    p, cv = curve_fit(f, x, h, p0=[h.max(), 0, 0.007], sigma=np.sqrt(np.maximum(h2, 1.0)), absolute_sigma=True)
    return float(abs(p[2]))
SIG_U = gauss_core(*reso_sample())
P("   step template width (edgeu rule): %.4f MeV/c2 (edgeu 6.70)" % (1e3 * SIG_U))

g = C.Grid(-0.200, 0.0025, 140, 4)
T = C.step_gauss(g, SIG_U)
ud = C.u_data(D)
S0d = (D["bits"] & 1) != 0
H = C.hist_data(D, ud, S0d, g)
SL1 = C.slices_of(0.30, 0.40, [0.30, 0.34, 0.37, 0.39, 0.40])
CF = np.array([C.CF_E1, C.CF_E2, C.CF_E2])
Hs0 = C.hist_sim(S, us, mpp, w0, g)
def simpart(Hs, groups, sl, cls=(0, 1, 2), c2fac=None):
    y = np.zeros(g.nb)
    for gname in groups:
        gi = C.G[gname]
        for (a, b) in sl:
            for k in cls:
                y = y + Hs[gi, a:b, k].sum(0) * (c2fac[k] if (gname == "C2" and c2fac is not None) else 1.0)
    return y
sp1 = C.data_spec(H, SL1, "all", g)
peakT = simpart(Hs0, ["C4", "C5"], SL1)
R0 = C.fit_step(sp1, T, g, -0.045, 0.045, 2, fixed=C.AP * peakT)
lin = dict(C2=C.h_lin(R0, simpart(Hs0, ["C2"], SL1)), C2c=C.h_lin(R0, simpart(Hs0, ["C2"], SL1, c2fac=CF)), C3=C.h_lin(R0, simpart(Hs0, ["C3"], SL1)))
ref = json.load(open(C.OUT + "edgeu_ana.json"))["nom"]
P("   nominal edgeu fit rebuilt: h_data %.3f +- %.3f, chi2 %.2f/%d; C2 raw %.3f, C2 x 0.879/0.8825 %.3f, C3 %.3f; f(D35) %.5f" %
  (R0.h, R0.eh, R0.chi2, R0.ndf, lin["C2"], lin["C2c"], lin["C3"], R0.h / (lin["C2c"] + lin["C3"])))
P("   edgeu_ana.json:            h_data %.3f, C2 raw %.3f, C2c %.3f, C3 %.3f; f %.5f" % (ref["hd"], ref["lin"]["C2"], ref["lin"]["C2c"], ref["lin"]["C3"], ref["f"]))
assert abs(R0.h - ref["hd"]) < 1e-6 * ref["hd"] and abs(lin["C3"] - ref["lin"]["C3"]) < 1e-6 * ref["lin"]["C3"]

# ---------------------------------------------------------------- 2. D38: per-count QED weights of C2, f_eta(u-edge)
Wt = np.load("/hdd3/przygoda/HADES/FAT_out/all_report/c2recon/c2recon_weights.npz")
key_w = {}
for i in range(len(Wt["tree"])):
    key_w[(int(Wt["tree"][i]), int(Wt["grp"][i]), float(np.float32(Wt["m4"][i])), float(np.float32(Wt["m4_true"][i])))] = i
Q = np.ones(S["n"])
c2idx = np.where(gS == 0)[0]
nj = 0
for i in c2idx:
    j = key_w.get((int(S["tree"][i]), int(S["grp"][i]), float(S["m4"][i]), float(S["m4t"][i])), -1)
    if j >= 0:
        Q[i] = Wt["q_A"][j]; nj += 1
nused = int(((w0 > 0) & (gS == 0)).sum())
nused_j = int(((w0 > 0) & (gS == 0) & (Q != 1.0)).sum())
P("   D38 weights: %d of %d simulated C2 counts joined (%d with S0, mode A, w > 0; joined %d); mean weight over S0 C2 %.4f" %
  (nj, len(c2idx), nused, nused_j, float((w0 * Q)[gS == 0].sum() / w0[gS == 0].sum())))
assert nj == len(c2idx)
wa = w0 * Q                                   # model (a): SMASH as it is, C2 with the D38 weights
Hsa = C.hist_sim(S, us, mpp, wa, g)
h2q = C.h_lin(R0, simpart(Hsa, ["C2"], SL1))
F38 = R0.h / (h2q + lin["C3"])
P("   f_eta(u-edge) with the D38 weights (frozen edgeu estimator): h(C2, D38) = %.3f, h(C3) = %.3f, h_sim = %.3f -> f = %.5f "
  "(c2recon_edgeu.json 0.42574; effective C2 factor %.4f)" % (h2q, lin["C3"], h2q + lin["C3"], F38, h2q / lin["C2"]))
RES["edgeu_repro"] = dict(hd=R0.h, ehd=R0.eh, chi2=R0.chi2, ndf=R0.ndf, C2=lin["C2"], C2c=lin["C2c"], C3=lin["C3"], f_D35=R0.h / (lin["C2c"] + lin["C3"]),
                          h2q=h2q, f_D38=F38, c2eff=h2q / lin["C2"], sigma_tmpl=SIG_U)

# ---------------------------------------------------------------- 3. model (b) weights
jsl = slice_defs(mpp)
xt = (m4p - 0.55) / 0.1
isGD = (gS == C.G["GDom"]) | (gS == C.G["GDnom"])
isGC = (gS == C.G["GCom"]) | (gS == C.G["GCnom"])
cA = CJ[np.clip(jsl, 0, 4)] * (1.0 + BT * xt)
sc = np.ones(S["n"])
sc[(gS == 0) | (gS == 1)] = F38
sc[(gS == 2) | (gS == 3)] = C.AP
sc[gS == 4] = AG
sc[isGD] = cA[isGD]
sc[isGC] = cA[isGC] * KB
wb = wa * sc
Hsb = C.hist_sim(S, us, mpp, wb, g)
# Jacobian of the non-eta part of (b) with respect to (c_1..c_5, b, kappa_B)
ne = isGD | isGC
base_ne = np.where(isGC, KB, 1.0)
jac_w = []
for jj in range(5):
    jac_w.append(np.where(ne & (jsl == jj), wa * (1.0 + BT * xt) * base_ne, 0.0))
jac_w.append(np.where(ne, wa * CJ[np.clip(jsl, 0, 4)] * xt * base_ne, 0.0))
jac_w.append(np.where(isGC, wa * cA, 0.0))
Hjac = [C.hist_sim(S, us, mpp, wj, g) for wj in jac_w]
P("   model (b) weights built (%.0f s)" % (time.time() - T0))

# ---------------------------------------------------------------- 4. D40 peak-tail variants (peaktail tail flag)
S2 = PT.load_sim()
pk2 = np.where((S2["g"] == 2) | (S2["g"] == 3))[0]
S2 = {k: (v[pk2] if isinstance(v, np.ndarray) else v) for k, v in S2.items()}
S2["n"] = len(pk2)
w2 = PT.sim_weight(S2, c2corr=False)
us2, mpp2, _ = C.u_sim(S2)
M4r, J = PT.m4_and_jacobian(S2)
dl = C.S_NOM * (J[0] * (S2["ep_p"] - S2["ep_ptrue"]) + J[1] * (S2["em_p"] - S2["em_ptrue"])).astype(np.float64) * 1e-3
flag = (dl < -0.005).astype(np.float64)
del M4r, J
def wt_t(t, lo, hi):
    """peaktail wt_t: tail of each peak class scaled by t, the number of counts of the class in lo <= M(pipi)' < hi kept"""
    out = w2.copy()
    inR = (mpp2 >= lo) & (mpp2 < hi)
    for gg in (2, 3):
        m = (S2["g"] == gg) & inR & (w2 > 0)
        N = w2[m].sum(); Nt = (w2 * flag)[m].sum(); Nc = N - Nt
        c = (N - t * Nt) / Nc
        out[m] = w2[m] * ((1 - flag[m]) * c + t * flag[m])
    return out
P("   peaktail counts (C4 + C5): %d; tail share (0.30-0.40) C4 %.4f, C5 %.4f (peaktail: 0.263, 0.240)" %
  (len(pk2), *[float((w2 * flag)[(S2["g"] == gg) & (mpp2 >= 0.30) & (mpp2 < 0.40)].sum() / w2[(S2["g"] == gg) & (mpp2 >= 0.30) & (mpp2 < 0.40)].sum()) for gg in (2, 3)]))


# ---------------------------------------------------------------- helpers: covariances
def region_bins(lo, hi):
    return C.PB(lo), C.PB(hi)
def cov_groups(keys, bins, wts, nb):
    """sum over groups of s_i s_j, s_i = sum of the weights of the group in bin i"""
    ok = (bins >= 0) & (wts != 0)
    k = keys[ok]; b = bins[ok]; ww = wts[ok]
    uq, inv = np.unique(k, return_inverse=True)
    M = sps.csr_matrix((ww, (inv, b)), shape=(uq.size, nb))
    return (M.T @ M).toarray()
def data_cov_N(region, cls):
    a, b = region_bins(*region)
    pb = C.ppbin_of(D["mpipi"]); bb = g.bin(ud)
    m = (D["k"] == 0) & S0d & (pb >= a) & (pb < b) & np.isin(C.cls_of(D["mee"]), C.CLS_SETS[cls])
    return cov_groups(D["grp"][m].astype(np.int64), bb[m], np.ones(int(m.sum())), g.nb)
def sim_cov(wts, region, cls):
    a, b = region_bins(*region)
    pb = C.ppbin_of(mpp.astype(np.float32)); bb = g.bin(us)
    m = (pb >= a) & (pb < b) & np.isin(C.cls_of(S["mee"]), C.CLS_SETS[cls]) & (wts != 0)
    key = S["tree"][m].astype(np.int64) * (1 << 32) + S["grp"][m].astype(np.int64)
    return cov_groups(key, bb[m], wts[m], g.nb)
def rsum(Hs, region, cls, gi=None):
    a, b = region_bins(*region)
    cs = C.CLS_SETS[cls]
    if gi is None:
        return Hs[:, a:b][:, :, cs].sum((1, 2))      # [9][nb]
    return Hs[gi, a:b][:, cs].sum((0, 1))
def kt_hist(t, region, cls):
    lo, hi = region
    wt = wt_t(t, lo, hi)
    a, b = region_bins(lo, hi)
    pb = C.ppbin_of(mpp2.astype(np.float32)); bb = g.bin(us2)
    m = (pb >= a) & (pb < b) & np.isin(C.cls_of(S2["mee"]), C.CLS_SETS[cls]) & (bb >= 0) & (wt != 0)
    return np.bincount(bb[m], weights=wt[m], minlength=g.nb).astype(np.float64)


# ---------------------------------------------------------------- 5. spectra, covariances and numbers per region
REGIONS = {"r1": dict(range=(0.30, 0.40), slices=SL1, text="0.30 <= M(pipi) < 0.40 (edgeu fit region)"),
           "all": dict(range=(0.27, 0.415), slices=[(C.PB(0.27), C.PB(0.30))] + SL1 + [(C.PB(0.40), C.PB(0.415))],
                       text="all M(pipi) with a defined edge (0.27-0.415; u exists below 0.4116)")}
ip = np.where((g.x >= UPLOT[0] - 1e-12) & (g.x < UPLOT[1] - 1e-12))[0]
FD["x"] = g.x[ip]; FD["edges"] = g.edges[ip[0]:ip[-1] + 2]
GN = C.GROUPS


def chi2_of(r, V, idx):
    Vi = np.linalg.inv(V[np.ix_(idx, idx)])
    rr = r[idx]
    return float(rr @ Vi @ rr), len(idx)
def ratio_band(Dv, Mv, VD, VM, VA, idx):
    o = np.zeros(len(Dv)); o[idx] = 1.0
    d = Dv @ o; mm = Mv @ o
    vd = o @ VD @ o; vm = o @ VM @ o; va = o @ VA @ o
    R = d / mm
    e_stat = R * math.sqrt(vd / d ** 2 + vm / mm ** 2)
    e_A = R * math.sqrt(va) / mm
    return dict(R=float(R), stat=float(e_stat), stageA=float(e_A), data=float(d), model=float(mm))
def idx_of(lo, hi):
    return np.where((g.x >= lo - 1e-12) & (g.x < hi - 1e-12))[0]


for rk, RG in REGIONS.items():
    RES[rk] = dict(text=RG["text"])
    for cls in ("all", "E1", "E2"):
        sp = C.data_spec(H, RG["slices"], cls, g)
        VN = data_cov_N(RG["range"], cls)
        VD = VN + sp.cbcov
        Ga = rsum(Hsa, RG["range"], cls); Gb = rsum(Hsb, RG["range"], cls)
        Ma = Ga.sum(0); Mb = Gb.sum(0)
        VMa = sim_cov(wa, RG["range"], cls); VMb = sim_cov(wb, RG["range"], cls)
        Jm = np.vstack([rsum(Hj, RG["range"], cls).sum(0) for Hj in Hjac]).T      # nb x 7
        VA = Jm @ COVA @ Jm.T
        # D40 peak-tail variants of (b)
        K1 = kt_hist(1.0, RG["range"], cls)
        Kref = Gb[C.G["C4"]] + Gb[C.G["C5"]]
        kmax = float(np.max(np.abs(C.AP * K1 - Kref)))
        Mbt = {t: Mb + C.AP * (kt_hist(t, RG["range"], cls) - K1) for t in TBAND}
        out = dict(K1_check_maxdiff=kmax, N_over_poisson=float(np.trace(VN[np.ix_(ip, ip)]) / sp.N[ip].sum()))
        for mk, Mv, VM, VAx in (("a", Ma, VMa, np.zeros_like(VA)), ("b", Mb, VMb, VA)):
            r = sp.sig - Mv
            res = {}
            for nm, (lo, hi) in RANGES.items():
                ii = idx_of(lo, hi)
                c2s, n = chi2_of(r, VD + VM, ii)
                c2A, _ = chi2_of(r, VD + VM + VAx, ii)
                # one free factor lam on the non-eta part (shape test of the smooth part): D - (M - B) = lam B
                Bv = (Gb if mk == "b" else Ga)[[C.G["GDom"], C.G["GDnom"], C.G["GCom"], C.G["GCnom"]]].sum(0)
                Vi = np.linalg.inv((VD + VM)[np.ix_(ii, ii)])
                yv = (sp.sig - (Mv - Bv))[ii]; bv = Bv[ii]
                lam = float(bv @ Vi @ yv / (bv @ Vi @ bv)); elam = float(1.0 / math.sqrt(bv @ Vi @ bv))
                rl = yv - lam * bv
                res[nm] = dict(chi2_stat=c2s, chi2_stageA=c2A, n=n, lam=lam, elam=elam, chi2_lam=float(rl @ Vi @ rl))
            bands = {nm: ratio_band(sp.sig, Mv, VD, VM, VAx, idx_of(lo, hi)) for nm, (lo, hi) in BANDS.items()}
            Bfull = (Gb if mk == "b" else Ga)[[C.G["GDom"], C.G["GDnom"], C.G["GCom"], C.G["GCnom"]]].sum(0)
            Mlam = Mv + (res["fit"]["lam"] - 1.0) * Bfull
            bands_lam = {nm: ratio_band(sp.sig, Mlam, VD, VM, VAx, idx_of(lo, hi)) for nm, (lo, hi) in BANDS.items()}
            out[mk] = dict(chi2=res, bands=bands, bands_lam=bands_lam)
            if mk == "b":
                out["b"]["bands_t"] = {str(t): {nm: ratio_band(sp.sig, Mbt[t], VD, VM, VAx, idx_of(lo, hi))["R"] for nm, (lo, hi) in BANDS.items()} for t in TBAND}
                out["b"]["chi2_t"] = {str(t): {nm: chi2_of(sp.sig - Mbt[t], VD + VM, idx_of(lo, hi))[0] for nm, (lo, hi) in RANGES.items()} for t in TBAND}
        # composition in the bands (model b)
        out["composition_b"] = {nm: {GN[i]: float(Gb[i][idx_of(lo, hi)].sum()) for i in range(9)} for nm, (lo, hi) in BANDS.items()}
        out["composition_a"] = {nm: {GN[i]: float(Ga[i][idx_of(lo, hi)].sum()) for i in range(9)} for nm, (lo, hi) in BANDS.items()}
        out["data_bands"] = {nm: dict(N=float(sp.N[idx_of(lo, hi)].sum()), CB=float(sp.cb[idx_of(lo, hi)].sum())) for nm, (lo, hi) in BANDS.items()}
        RES[rk][cls] = out
        tag = "%s_%s_" % (rk, cls)
        FD[tag + "N"] = sp.N[ip]; FD[tag + "cb"] = sp.cb[ip]; FD[tag + "sig"] = sp.sig[ip]
        FD[tag + "VD"] = VD[np.ix_(ip, ip)]
        for i in range(9):
            FD[tag + "a_" + GN[i]] = Ga[i][ip]; FD[tag + "b_" + GN[i]] = Gb[i][ip]
        FD[tag + "VMa"] = VMa[np.ix_(ip, ip)]; FD[tag + "VMb"] = VMb[np.ix_(ip, ip)]; FD[tag + "VA"] = VA[np.ix_(ip, ip)]
        for t in TBAND:
            FD[tag + "Mb_t%.1f" % t] = Mbt[t][ip]
        P("\n== region %s, M(ee) class %s: %s" % (rk, cls, RG["text"]))
        P("   peak template check: max |a_P K(t=1) [peaktail dump] - (C4 + C5) x a_P [edgeu columns]| = %.2e counts; data N event-group variance / N = %.4f" %
          (kmax, out["N_over_poisson"]))
        for mk in ("a", "b"):
            o = out[mk]
            P("   model (%s): chi2 (stat: data + MC) -45..+45: %.1f/%d, +20..+60: %.1f/%d; with the stage-A error: %.1f, %.1f" %
              (mk, o["chi2"]["fit"]["chi2_stat"], o["chi2"]["fit"]["n"], o["chi2"]["peak"]["chi2_stat"], o["chi2"]["peak"]["n"],
               o["chi2"]["fit"]["chi2_stageA"], o["chi2"]["peak"]["chi2_stageA"]))
            P("               free factor on G_D + G_C: -45..+45 lam %.4f +- %.4f (chi2 %.1f), +20..+60 lam %.4f +- %.4f (chi2 %.1f)" %
              (o["chi2"]["fit"]["lam"], o["chi2"]["fit"]["elam"], o["chi2"]["fit"]["chi2_lam"], o["chi2"]["peak"]["lam"], o["chi2"]["peak"]["elam"], o["chi2"]["peak"]["chi2_lam"]))
            P("               data/model: " + "; ".join("%s %.4f +- %.4f (stat) +- %.4f (stage A)" % (nm, o["bands"][nm]["R"], o["bands"][nm]["stat"], o["bands"][nm]["stageA"]) for nm in BANDS))
            P("               data/[model with G_D + G_C x lam(-45..+45)]: " + "; ".join("%s %.4f +- %.4f (stat)" % (nm, o["bands_lam"][nm]["R"], o["bands_lam"][nm]["stat"]) for nm in BANDS))
            if mk == "b":
                for t in TBAND:
                    P("               D40 t = %.1f: data/model %s; chi2 -45..+45 %.1f, +20..+60 %.1f" %
                      (t, ", ".join("%s %.4f" % (nm, o["bands_t"][str(t)][nm]) for nm in BANDS), o["chi2_t"][str(t)]["fit"], o["chi2_t"][str(t)]["peak"]))
        for nm in BANDS:
            cb = out["composition_b"][nm]; ca = out["composition_a"][nm]
            P("   band %-5s: data N %.0f, CB %.0f, N - CB %.0f | model (b): bump C2+C3 %.0f, peak C4+C5 %.0f, G_eta %.0f, G_D %.0f, G_C %.0f | model (a): bump %.0f, peak %.0f, G_eta %.0f, G_D %.0f, G_C %.0f" %
              (nm, out["data_bands"][nm]["N"], out["data_bands"][nm]["CB"], out["data_bands"][nm]["N"] - out["data_bands"][nm]["CB"],
               cb["C2"] + cb["C3"], cb["C4"] + cb["C5"], cb["Geta"], cb["GDom"] + cb["GDnom"], cb["GCom"] + cb["GCnom"],
               ca["C2"] + ca["C3"], ca["C4"] + ca["C5"], ca["Geta"], ca["GDom"] + ca["GDnom"], ca["GCom"] + ca["GCnom"]))
    P("   (%.0f s)" % (time.time() - T0))

# ---------------------------------------------------------------- 6. check of the stage-A implementation in CR-U (M4, E2), per slice
# S5 (M(pipi) >= 0.43) is open-ended and the edgeu columns stop at 0.46, so the all-slice numbers of Table C3 cannot be
# rebuilt; c_j is free per slice, so each of S1..S4 must come out close to 1 within the noise of the raw templates.
gM = C.Grid(0.400, 0.0025, 120, 4)
HM = C.hist_data(D, D["m4"].astype(np.float64), S0d, gM)
cs = C.CLS_SETS["E2"]
iu = np.where((gM.x >= 0.575) & (gM.x < 0.700))[0]
inU = np.zeros(gM.nb); inU[iu] = 1.0
HMa = C.hist_sim(S, m4p, mpp, wa, gM); HMb = C.hist_sim(S, m4p, mpp, wb, gM)
HJ = [C.hist_sim(S, m4p, mpp, wj, gM) for wj in jac_w]
pbS = C.ppbin_of(mpp.astype(np.float32)); bM = gM.bin(m4p)
keyS = S["tree"].astype(np.int64) * (1 << 32) + S["grp"].astype(np.int64)
RES["crU_check"] = {}
for lab, (lo, hi) in (("S1 0.27-0.34", (0.27, 0.34)), ("S2 0.34-0.37", (0.34, 0.37)), ("S3 0.37-0.39", (0.37, 0.39)), ("S4 0.39-0.43", (0.39, 0.43)),
                      ("S1-S4 0.27-0.43", (0.27, 0.43)), ("S5 part 0.43-0.46 (incomplete)", (0.43, 0.46))):
    a0, b0 = C.PB(lo), C.PB(hi)
    sls = [s_ for s_ in [(C.PB(0.27), C.PB(0.34)), (C.PB(0.34), C.PB(0.37)), (C.PB(0.37), C.PB(0.39)), (C.PB(0.39), C.PB(0.43)), (C.PB(0.43), C.PB(0.46))]
           if s_[0] >= a0 and s_[1] <= b0]
    spM = C.data_spec(HM, sls, "E2", gM, cbk=KE2, clsrule=1)
    msum = lambda Hs, gi: Hs[gi, a0:b0][:, cs].sum((0, 1))
    BA = sum(msum(HMb, C.G[n]) for n in ("GDom", "GDnom", "GCom", "GCnom"))
    ETA = sum(msum(HMa, C.G[n]) for n in ("C2", "C3", "C4", "C5", "Geta"))
    GETA = msum(HMa, C.G["Geta"])
    d_ = spM.sig[iu].sum()
    m1 = (BA + 0.4 * ETA)[iu].sum(); m2 = (BA + 0.4 * (ETA - GETA) + AG * GETA)[iu].sum()
    er = math.sqrt(spM.N[iu].sum() + np.diag(spM.cbcov)[iu].sum()) / d_
    mU = (pbS >= a0) & (pbS < b0) & np.isin(C.cls_of(S["mee"]), cs) & (bM >= 0)
    wU = np.where(ne, wb, 0.0) + np.where(gS <= 4, 0.4 * wa, 0.0)
    VU = cov_groups(keyS[mU], bM[mU], wU[mU], gM.nb)
    JU = np.vstack([Hj[:, a0:b0][:, :, cs].sum((1, 2)).sum(0) for Hj in HJ]).T
    e_mc = math.sqrt(inU @ VU @ inU) / m1
    e_A = math.sqrt((inU @ JU) @ COVA @ (JU.T @ inU)) / m1
    RES["crU_check"][lab] = dict(R_g04=d_ / m1, R_g0656=d_ / m2, rel_stat=er, rel_mc=e_mc, rel_stageA=e_A, N=float(spM.N[iu].sum()),
                                 CB=float(spM.cb[iu].sum() / KE2), Geta_share=float((AG * GETA[iu]).sum() / m2))
    P("   CR-U check, E2, %-31s: (N - k_E2 CB)/(B_A + 0.4 eta) = %.4f +- %.4f (data) +- %.4f (MC, raw templates) [stage-A error %.4f]; "
      "with G_eta x 0.656: %.4f (G_eta share %.3f); N %.0f, CB %.0f" %
      (lab, d_ / m1, d_ / m1 * er, d_ / m1 * e_mc, e_A, d_ / m2, RES["crU_check"][lab]["Geta_share"], spM.N[iu].sum(), spM.cb[iu].sum() / KE2))

json.dump(RES, open(O + "udatasim_overlay.json", "w"), indent=1, default=float)
FD["F38"] = F38; FD["AP"] = C.AP; FD["AG"] = AG; FD["SIG_U"] = SIG_U
np.savez_compressed(O + "udatasim_figdata.npz", **FD)
P("done (%.0f s)" % (time.time() - T0))
