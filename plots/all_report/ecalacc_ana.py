#!/usr/bin/env python3
# ============================================================================
# ecalacc_ana.py -- package ecalacc (2026-11-15): numbers and tables of the ECAL-photon feasibility study
# ============================================================================
# Inputs (/hdd3/przygoda/HADES/FAT_out/all_report/ecalacc/, written by ecalacc_fill.C and ecalacc_geo.py):
#   rows_c23.bin + rows_c23_layout.txt  (one record per simulated C2/C3 row of "pre": weights, conditional photon
#                                        acceptance, truth photon of the ECAL simulation, compact and recomputed tags)
#   tag_fits.txt                         (data and simulated tag fits, pre / S0, M(pi+pi-) slices; evtStart split)
#   hist_data.txt, hist_sim.txt          (photon maps, M(e+e-gamma), distance to the edge)
#   geo_numbers.json
# Outputs: ecalacc_numbers.json (same directory); docs/latex/tables_all/ecalacc_{results,tag}.{txt,tex}
# Typed inputs (sources): f_eta(edge) and its errors, C2 factor 0.8855, kappa 0.48: bedge_results.md (headline table);
#   T(pre) = 3314 +- 105 +- 124: report_exp.tex Sect. "ECAL-tagged eta"; MC error of the simulated tag 603/6804 (WP6 S0,
#   frozen fit, 100 replicas): wp6_results.md 2.2; visible bump yields: bedge_results.md 2.7.
# Run: python3 -I plots/all_report/ecalacc_ana.py
# ============================================================================
import json
import math
import numpy as np

O = "/hdd3/przygoda/HADES/FAT_out/all_report/ecalacc"
T = "/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all"

# ---- typed inputs ----------------------------------------------------------------------------------------------
F_EDGE = dict(asis=(0.3916, 0.0210, 0.0321, 0.0425), c2=(0.4184, 0.0225, 0.0348, 0.0471), k48=(0.5602, 0.0301, 0.0639, 0.0812))
C2FAC = 0.8855          # recommended C2 gamma* angle factor at the edge (bedge, mean of WP8 toy and closure)
KAPPA = 0.48
T_DATA_PRE = (3314.0, 105.0, 124.0)      # report_exp (stat, model)
T_SIM_MC_REL = 603.0 / 6804.0            # WP6: MC error of the simulated frozen-fit tag, S0 all M4
Y_BUMP_ALLM4 = 26959.0                   # bedge: visible bump yield (S0, all M4, f x (C2 + C3), rows x 64/63)
Y_BUMP_B = 21445.0

# ---- records ---------------------------------------------------------------------------------------------------
lay = open(O + "/rows_c23_layout.txt").read().split("\n")
n, nd = int(lay[0].split()[1]), int(lay[0].split()[3])
names = lay[1].split()[1:]
lam = [float(x) for x in lay[2].split()[1:]]
thr = [float(x) for x in lay[3].split()[1:]]
a = np.fromfile(O + "/rows_c23.bin", dtype=np.float64).reshape(n, nd)
R = {k: a[:, i] for i, k in enumerate(names)}
NEB, EBW, NTH = 30, 20.0, 45
EH = np.stack([R["eHist%d" % i] for i in range(NEB)], 1)
EHA = np.stack([R["eHistAcc%d" % i] for i in range(NEB)], 1)
TH = np.stack([R["thHist%d" % i] for i in range(NTH)], 1)
w = R["w"]
ok = R["ok"] > 0
fr5 = (R["m5cC"] >= 0.49) & (R["m5cC"] < 0.62)


def frac(x, m):
    """weighted mean of x over the rows m, with its MC error (independent rows)"""
    W = w[m].sum()
    if W <= 0:
        return (0.0, 0.0)
    F = (w[m] * x[m]).sum() / W
    v = (w[m] ** 2 * (x[m] - F) ** 2).sum() / W ** 2
    return (float(F), float(math.sqrt(v)))


def wsum(x, m):
    return float((w[m] * x[m]).sum()), float(math.sqrt((w[m] ** 2 * x[m] ** 2).sum()))


one = np.ones(n)
num = {"n_records": int(n), "lambda": lam, "thresholds": thr}
sels = {"pre": R["pre"] > 0, "S0": R["s0"] > 0}
classes = {"C2": R["cls"] == 2, "C3": R["cls"] == 3, "C23": (R["cls"] == 2) | (R["cls"] == 3)}
regions = {"all": np.ones(n, bool), "B": (R["m4p"] >= 0.455) & (R["m4p"] < 0.5225),
           "edge15": (R["dist"] >= -0.015) & (R["dist"] < 0.005), "far": (R["dist"] < -0.06) & (R["dist"] > -1)}

tagTrueW = (R["tagTrueC"] > 0) & fr5
tagW = (R["tagC"] > 0) & fr5
found = (R["tMatch"] > 0)
foundQ = found & (R["tPass"] > 0)

quant = {
    "geo": R["fGeo"], "geo_strict": R["fGeoStrict"], "geo_6sec": R["fGeo6"], "geo_flat": R["fGeoFlat"], "geo_kloe": R["fGeoKloe"],
    "geo_E60": R["fGeoThr60"], "geo_E80": R["fGeoThr80"], "geo_E100": R["fGeoThr100"], "geo_E110": R["fGeoThr110"], "geo_E130": R["fGeoThr130"],
    "found": found.astype(float), "found_quality": foundQ.astype(float), "tag_any": (R["tagC"] > 0).astype(float),
    "tag_any_m5c": tagW.astype(float), "tag_true_m5c": tagTrueW.astype(float),
}
for l in range(len(lam)):
    m5 = (R["m5cL%d" % l] >= 0.49) & (R["m5cL%d" % l] < 0.62)
    quant["tagL%d_m5c" % l] = ((R["tagL%d" % l] > 0) & m5).astype(float)
    quant["tagTrueL%d_m5c" % l] = ((R["tagTrueL%d" % l] > 0) & m5).astype(float)
fr = {}
for s, ms in sels.items():
    for c, mc in classes.items():
        for g, mg in regions.items():
            m = ms & mc & mg & ok
            key = "%s_%s_%s" % (s, c, g)
            fr[key] = {"W": float(w[ms & mc & mg].sum()), "Nrows": int((ms & mc & mg).sum()), "ok": frac(ok.astype(float), ms & mc & mg)[0]}
            for q, x in quant.items():
                fr[key][q] = frac(x, m if q.startswith("geo") else ms & mc & mg)
            if c == "C2":
                # QED gamma* weighting: acceptance fraction and the C2 angle factor of the selection
                xq = R["fGeoQed"] * R["meanQed"]
                Wq = (w[m] * R["meanQed"][m]).sum()
                fr[key]["geo_qed"] = (float((w[m] * xq[m]).sum() / Wq), 0.0)
                fr[key]["c_qed_circle"] = frac(R["meanQed"], m)
num["fractions"] = fr

# C2 angle factor of the tagged rows (exact for a found photon, circle average otherwise)
qtag = np.where(R["tMatch"] > 0, R["tQed"], R["meanQed"])
for s, ms in sels.items():
    m = ms & classes["C2"] & tagW
    num["c_tag_%s" % s] = frac(qtag, m)
    m2 = ms & classes["C2"] & tagTrueW
    num["c_tagtrue_%s" % s] = frac(R["tQed"], m2)

# per M(pipi) slice and per distance to the edge (pre, C2+C3 and C2)
slices = [(0.0, 0.34), (0.34, 0.37), (0.37, 0.40), (0.40, 0.43)]
dbins = [(-0.20, -0.10), (-0.10, -0.06), (-0.06, -0.04), (-0.04, -0.03), (-0.03, -0.02), (-0.02, -0.01), (-0.01, 0.0), (0.0, 0.02)]
def table_over(binsel):
    out = []
    for lo, hi, m0 in binsel:
        row = {"lo": lo, "hi": hi}
        for c in ("C23", "C2"):
            mm = sels["pre"] & classes[c] & m0
            row[c] = {"W": float(w[mm].sum()), "N": int(mm.sum())}
            for q in ("geo", "geo_6sec", "geo_E100", "found", "found_quality", "tag_true_m5c", "tag_any_m5c", "tagTrueL3_m5c"):
                x = quant[q]
                row[c][q] = frac(x, mm & ok if q.startswith("geo") else mm)
            row[c]["eMean"] = frac(R["eMean"], mm & ok)
            row[c]["eMeanAcc"] = frac(R["eMeanAcc"], mm & ok & (R["fGeo"] > 0))
        out.append(row)
    return out
num["vs_mpipi"] = table_over([(lo, hi, (R["mpp"] >= lo) & (R["mpp"] < hi)) for lo, hi in slices])
num["vs_dist"] = table_over([(lo, hi, (R["dist"] >= lo) & (R["dist"] < hi)) for lo, hi in dbins])

# photon energy spectra (pre, C2+C3): all photons, in the ECAL modules, found (any / quality), per 20 MeV of TRUE energy
mm = sels["pre"] & classes["C23"] & ok
eAll = (w[mm, None] * EH[mm]).sum(0)
eAcc = (w[mm, None] * EHA[mm]).sum(0)
def ehist(msk, field="tE"):
    h = np.zeros(NEB); v = np.zeros(NEB)
    idx = np.clip((R[field][msk] / EBW).astype(int), 0, NEB - 1)
    np.add.at(h, idx, w[msk]); np.add.at(v, idx, w[msk] ** 2)
    return h, v
eF, eFv = ehist(sels["pre"] & classes["C23"] & found & (R["tGeo"] == 1))
eFQ, eFQv = ehist(sels["pre"] & classes["C23"] & foundQ & (R["tGeo"] == 1))
num["espec"] = dict(edges=[i * EBW for i in range(NEB + 1)], all=eAll.tolist(), acc=eAcc.tolist(), found=eF.tolist(), found_var=eFv.tolist(),
                    foundQ=eFQ.tolist(), foundQ_var=eFQv.tolist())
eff = np.where(eAcc > 0, eF / np.maximum(eAcc, 1e-12), 0.0)
effQ = np.where(eAcc > 0, eFQ / np.maximum(eAcc, 1e-12), 0.0)
num["espec"]["eff"] = eff.tolist(); num["espec"]["effQ"] = effQ.tolist()
# theta spectrum of all photons of accepted rows
thAll = (w[mm, None] * TH[mm]).sum(0)
num["thspec"] = dict(edges=[2.0 * i for i in range(NTH + 1)], all=thAll.tolist())
thc = np.arange(NTH) * 2.0 + 1.0
Wt = thAll.sum() + 1e-12
num["theta_shares"] = dict(below11p5=float(thAll[thc < 11.5].sum() / Wt), band=float(thAll[(thc >= 11.5) & (thc < 47.5)].sum() / Wt),
                           above47p5=float(thAll[thc >= 47.5].sum() / Wt), hades18_85=float(thAll[(thc >= 18) & (thc < 85)].sum() / Wt))
# energy shares
ec = (np.arange(NEB) + 0.5) * EBW
num["energy_shares"] = dict(all_below100=float(eAll[ec < 100].sum() / eAll.sum()), acc_below100=float(eAcc[ec < 100].sum() / eAcc.sum()),
                            mean_all=float((eAll * ec).sum() / eAll.sum()), mean_acc=float((eAcc * ec).sum() / eAcc.sum()))
# efficiency inside the modules
for s in ("pre", "S0"):
    f_ = fr["%s_C23_all" % s]
    num["eff_in_acc_%s" % s] = dict(found=f_["found"][0] / f_["geo"][0], found_quality=f_["found_quality"][0] / f_["geo"][0],
                                     tag_true=f_["tag_true_m5c"][0] / f_["geo"][0])

# closure of the conditional completion with the photons found by the ECAL simulation
mF = found
num["closure"] = dict(n_found=int(mF.sum()), inside=float(np.mean(R["tGeo"][mF] == 1)), deadmodule=float(np.mean(R["tGeo"][mF] == 0.5)),
                      outside=float(np.mean(R["tGeo"][mF] == 0)),
                      res_pi0_q=[float(x) for x in np.percentile(R["tRes1"][mF], [5, 50, 95])],
                      res_eta_q=[float(x) for x in np.percentile(R["tRes2"][mF], [5, 50, 95])],
                      cos_q=[float(x) for x in np.percentile(R["tCosRes"][mF], [5, 50, 95])],
                      frac_abs_cos_lt_1e3=float(np.mean(np.abs(R["tCosRes"][mF]) < 1e-3)),
                      ecl_over_etrue_q=[float(x) for x in np.percentile((R["tEcl"] / np.maximum(R["tE"], 1e-9))[mF], [16, 50, 84])],
                      ecl_over_etrue_q_quality=[float(x) for x in np.percentile((R["tEcl"] / np.maximum(R["tE"], 1e-9))[foundQ], [16, 50, 84])],
                      conv_cluster_share=frac(R["tConv"], sels["pre"] & classes["C23"]))
# where is the true photon found by the ECAL simulation relative to the circle prediction: mean predicted acceptance
num["closure"]["mean_pred_geo_of_found"] = frac(R["fGeo"], found & ok)
# calibration of the circle weighting: rows binned in the predicted probability fGeo; the found share divided by the mean
# predicted probability must be the same in every bin (= the efficiency inside the modules) if the density along the
# circle (SMASH eta spectrum) is right; the flat-weight variant is shown next to it
calib = []
for lo_, hi_ in ((0.0, 0.05), (0.05, 0.2), (0.2, 0.4), (0.4, 0.6), (0.6, 0.8), (0.8, 1.01)):
    mb = ok & sels["pre"] & classes["C23"] & (R["fGeo"] >= lo_) & (R["fGeo"] < hi_)
    if not mb.any():
        continue
    fF = frac(found.astype(float), mb); pG = frac(R["fGeo"], mb); pF = frac(R["fGeoFlat"], mb)
    calib.append(dict(lo=lo_, hi=hi_, N=int(mb.sum()), found=fF, pred=pG[0], predFlat=pF[0], ratio=fF[0] / pG[0] if pG[0] > 0 else 0.0,
                      eratio=fF[1] / pG[0] if pG[0] > 0 else 0.0))
num["closure"]["calib"] = calib

# ---- tag fits --------------------------------------------------------------------------------------------------
tf = {}
spec = {}
for l in open(O + "/tag_fits.txt"):
    p = l.split()
    if not p:
        continue
    if p[0] == "TAG":
        s, j = int(p[1]), int(p[2]); v = list(map(float, p[3:]))
        tf[(s, j)] = dict(dT=v[0], edT=v[1], chi2=v[2], ndf=int(v[3]), dfree_mean=v[4], dfree_sig=v[5], dfree_T=v[6], sT=v[7], sTfree=v[8],
                          sfree_mean=v[9], sfree_sig=v[10], trW=v[11], tr2=v[12], tr3=v[13], trA=v[14], c23=v[15], c2=v[16], c3=v[17], c23B=v[18])
    elif p[0] == "TAGSTART":
        tf[("start", int(p[1]), int(p[2]))] = dict(dT=float(p[3]), edT=float(p[4]))
    elif p[0] == "SPEC":
        spec["%s_%s" % (p[1], p[2])] = [float(x) for x in p[3:]]
    elif p[0] == "FITPAR":
        spec["fit_%s_%s" % (p[1], p[2])] = p[3:]
tsn = ["all", "<0.34", "0.34-0.37", "0.37-0.40", "0.40-0.43", ">=0.43"]
num["tagfits"] = {"%s_%s" % (["pre", "S0"][s], tsn[j]): v for (s, j), v in [(k, x) for k, x in tf.items() if k[0] in (0, 1)]}
num["tag_start"] = {("start1" if k[1] == 2 else "start0") + "_" + tsn[k[2]]: v for k, x in tf.items() if k[0] == "start" for v in [x]}

# ---- the frozen fit in Python (copy of scan::FitPeakArrays for a spectrum without CB; reproduces the C++ T_sim exactly)
from math import erf, sqrt
def gfrac(a_, b_, m_, s_):
    return 0.5 * (erf((b_ - m_) / (s_ * sqrt(2))) - erf((a_ - m_) / (s_ * sqrt(2))))
def tagfit(y, mean=0.5463, sig=0.0062, lo=0.49, hi=0.62, order=2):
    nb = len(y); bw = 0.4 / nb; x1 = 0.40 + np.arange(nb) * bw; xc = x1 + bw / 2
    idx = np.where((xc >= lo) & (xc < hi))[0]
    x0, dx = (lo + hi) / 2, (hi - lo) / 2
    F = np.zeros((len(idx), order + 2))
    F[:, 0] = [gfrac(x1[i], x1[i] + bw, mean, sig) for i in idx]
    t = (xc[idx] - x0) / dx
    for j in range(order + 1):
        F[:, 1 + j] = t ** j
    yy = y[idx]; var = np.maximum(yy, 1.0)
    for it in range(7):
        Wt = 1 / var; A = F.T @ (F * Wt[:, None]); b = F.T @ (yy * Wt); par = np.linalg.solve(A, b)
        var = np.maximum(F @ par, 1.0)
    return float(par[0])
def binm5(mask, field):
    h = np.zeros(160)
    x = R[field][mask]
    i = np.floor((x - 0.40) / 0.0025).astype(int)
    g = (i >= 0) & (i < 160)
    np.add.at(h, i[g], w[mask][g])
    return h
lamfit = {}
for s, sname in ((0, "pre"), (1, "S0")):
    y12 = np.array(spec["%d_sim_T12" % s]); yc = np.array(spec["%d_sim_C2" % s]) + np.array(spec["%d_sim_C3" % s])
    msel = sels[sname] & classes["C23"]
    c0 = binm5(msel & (R["tagL0"] > 0), "m5cL0")
    res = dict(T12=tagfit(y12), check_C23_binned_vs_spec=float(np.abs(c0 - yc).max()))
    for l in range(len(lam)):
        cl_ = binm5(msel & (R["tagL%d" % l] > 0), "m5cL%d" % l)
        res["L%d" % l] = tagfit(y12 - c0 + cl_)
        res["L%d_factor" % l] = res["L%d" % l] / res["T12"]
    lamfit[sname] = res
num["lambda_fit"] = lamfit

# ---- pi0 position and width in M(e+e-gamma): data (e+e- rows with a good photon, pre) and simulation (C2 + C3)
from scipy.optimize import curve_fit
def hread(fn):
    H = {}; Ls = open(fn).read().split("\n"); i = 0
    while i < len(Ls):
        if Ls[i].startswith("HIST"):
            _, nm, nb_, lo_, hi_ = Ls[i].split(); nb_ = int(nb_)
            v = np.array([list(map(float, Ls[i + 1 + k].split())) for k in range(nb_)])
            H[nm] = (float(lo_), float(hi_), v[:, 0], v[:, 1]); i += nb_ + 1
        else:
            i += 1
    return H
HD = hread(O + "/hist_data.txt"); HS = hread(O + "/hist_sim.txt")
def gpol(x, A, m_, s_, c0_, c1_):
    return A * np.exp(-0.5 * ((x - m_) / s_) ** 2) + c0_ + c1_ * (x - 0.14)
def pi0fit(lo_, hi_, y, v, flo=0.095, fhi=0.20):
    nb_ = len(y); x = lo_ + (np.arange(nb_) + 0.5) * (hi_ - lo_) / nb_
    g = (x >= flo) & (x < fhi)
    p0 = [y[g].max(), 0.137, 0.010, np.median(y[g]) * 0.3, 0.0]
    try:
        par, cov = curve_fit(gpol, x[g], y[g], p0=p0, sigma=np.sqrt(np.maximum(v[g], 1e-9)) if v is not None else None, maxfev=20000)
        return dict(mean=float(par[1]), sigma=float(abs(par[2])), emean=float(math.sqrt(cov[1, 1])), esigma=float(math.sqrt(cov[2, 2])))
    except Exception as e:
        return dict(mean=0.0, sigma=0.0, err=str(e))
dlo, dhi, dy, dv = HD["data_us_meeg"]
_, _, py_, pv_ = HD["data_pp_meeg"]; _, _, my_, mv_ = HD["data_mm_meeg"]
cbm = 2 * np.sqrt(np.maximum(py_, 0) * np.maximum(my_, 0))
pi0 = {"data": pi0fit(dlo, dhi, dy - cbm, dv + py_ + my_)}
slo, shi, sy, sv = HS["sim_C23_meeg"]
pi0["sim_L0"] = pi0fit(slo, shi, sy, sv)
for l in range(len(lam)):
    mm_ = sels["pre"] & classes["C23"] & (R["meegL%d" % l] > 0) & (R["mee"] < 0.1)
    h = np.zeros(100); hv = np.zeros(100)
    i = np.floor(R["meegL%d" % l][mm_] / 0.003).astype(int); g = (i >= 0) & (i < 100)
    np.add.at(h, i[g], w[mm_][g]); np.add.at(hv, i[g], w[mm_][g] ** 2)
    pi0["simrec_L%d" % l] = pi0fit(0.0, 0.30, h, hv)
def winacc(mean_, sig_):
    return gfrac(0.119, 0.173, mean_, sig_) if sig_ > 0 else 0.0
for k in list(pi0.keys()):
    pi0[k]["window_acceptance"] = winacc(pi0[k]["mean"], pi0[k]["sigma"])
num["pi0"] = pi0
RESFAC = pi0["data"]["window_acceptance"] / pi0["simrec_L3"]["window_acceptance"] if pi0["simrec_L3"]["window_acceptance"] > 0 else 1.0
num["resolution_factor"] = RESFAC

# expected tag count in the data at the measured eta scale
def expect(s, j, reading, lamIdx=0):
    t = tf[(s, j)]
    msel = sels[["pre", "S0"][s]]
    jm = np.ones(n, bool) if j == 0 else (R["mpp"] >= [0, 0, 0.34, 0.37, 0.40, 0.43][j]) & (R["mpp"] < [9, 0.34, 0.37, 0.40, 0.43, 9][j])
    # truth tagged C2 / C3 in m5c 0.49-0.62 for the calibration factor lambda
    m5 = (R["m5cL%d" % lamIdx] >= 0.49) & (R["m5cL%d" % lamIdx] < 0.62)
    tg = (R["tagL%d" % lamIdx] > 0) & m5 & msel & jm
    m50 = (R["m5cL0"] >= 0.49) & (R["m5cL0"] < 0.62)
    tg0 = (R["tagL0"] > 0) & m50 & msel & jm
    T2 = (w * tg * classes["C2"]).sum(); T3 = (w * tg * classes["C3"]).sum()
    T20 = (w * tg0 * classes["C2"]).sum(); T30 = (w * tg0 * classes["C3"]).sum()
    cq = frac(qtag, tg & classes["C2"])[0] if (tg & classes["C2"]).any() else C2FAC
    lamFacTruth = (T2 + T3) / (T20 + T30) if (T20 + T30) > 0 else 1.0
    # calibration factor: the frozen fit of the modified simulated spectrum (all M(pi+pi-)); per slice the truth ratio
    # of the tagged C2 + C3 rows times (fit factor / truth factor of all M(pi+pi-))
    sname = ["pre", "S0"][s]
    if j == 0:
        lamFac = lamfit[sname]["L%d_factor" % lamIdx]
    else:
        mall = (R["tagL%d" % lamIdx] > 0) & m5 & msel & classes["C23"]
        mall0 = (R["tagL0"] > 0) & m50 & msel & classes["C23"]
        tAll = (w * mall).sum() / (w * mall0).sum()
        lamFac = lamFacTruth * lamfit[sname]["L%d_factor" % lamIdx] / tAll
    if reading == "asis":
        f, comp = F_EDGE["asis"], 1.0
    elif reading == "c2":
        f, comp = F_EDGE["c2"], (cq * T2 + T3) / (T2 + T3)
    else:
        f, comp = F_EDGE["k48"], (cq * T2 + KAPPA * T3) / (T2 + T3)
    val = f[0] * t["sT"] * comp * lamFac
    rel = [f[1] / f[0], f[2] / f[0], math.hypot(f[3] / f[0], T_SIM_MC_REL)]
    return dict(value=val, stat=val * rel[0], syst=val * rel[1], mc=val * rel[2], comp=comp, lamFac=lamFac, lamFacTruth=lamFacTruth, cq=cq, T2=T2, T3=T3,
                sT=t["sT"], dT=t["dT"], edT=t["edT"], ratio=t["dT"] / val if val > 0 else 0.0)
ex = {}
for s in (0, 1):
    for j in range(5):
        for rd in ("asis", "c2", "k48"):
            for li in range(len(lam)):
                ex["%s_%s_%s_L%d" % (["pre", "S0"][s], tsn[j], rd, li)] = expect(s, j, rd, li)
num["expected"] = ex
# f_tag (data / simulation, frozen fits) per slice, and the edge per slice for comparison (bedge 2.5, uncorrected)
num["ftag_slices"] = {k: dict(ftag=v["dT"] / v["sT"] if v["sT"] > 0 else 0.0, eftag=v["edT"] / v["sT"] if v["sT"] > 0 else 0.0) for k, v in num["tagfits"].items()}
num["edge_slices_bedge"] = {"0.30-0.34": (0.342, 0.044), "0.34-0.37": (0.459, 0.039), "0.37-0.39": (0.376, 0.027), "0.39-0.40": (0.389, 0.062)}

# evtStart: tag rate per e+e- row for events with and without START (pre)
logtxt = open(O + "/ecalacc_fill.txt").read()
import re
mS = re.search(r"evtStart == 1 (\d+), != 1 (\d+) .*?bump window B: (\d+) / (\d+)", logtxt)
if mS:
    n1, n0, b1, b0 = map(float, mS.groups())
    t1 = tf[("start", 2, 0)]; t0 = tf[("start", 3, 0)]
    r_all = (t0["dT"] / n0) / (t1["dT"] / n1)
    r_B = (t0["dT"] / b0) / (t1["dT"] / b1)
    er = t0["edT"] / t0["dT"]
    num["start"] = dict(rows_start1=n1, rows_start0=n0, rowsB_start1=b1, rowsB_start0=b0, T1=t1["dT"], eT1=t1["edT"], T0=t0["dT"], eT0=t0["edT"],
                        ratio_rate_all=r_all, ratio_rate_B=r_B, ratio_err_rel=er,
                        loss_if_start0_like_start1=t0["dT"] * (1.0 / r_B - 1.0))
mC = re.search(r"C2 angle factor of ALL C2 rows \(pre\) from the circle average of the QED gamma\* weight: all ([0-9.]+), E1 ([0-9.]+), E2 ([0-9.]+)", logtxt)
if mC:
    num["c_qed_circle_all"] = [float(x) for x in mC.groups()]
mM = re.search(r"m_eta ([0-9.]+), m_pi0 ([0-9.]+)", logtxt)
if mM:
    num["generator_masses"] = [float(x) for x in mM.groups()]
num["geo"] = json.load(open(O + "/geo_numbers.json"))
json.dump(num, open(O + "/ecalacc_numbers.json", "w"), indent=1)

# ---- tables ----------------------------------------------------------------------------------------------------
def pc(x):
    return "%.1f +- %.1f" % (100 * x[0], 100 * x[1])
L = []
L.append("ecalacc: the photon of the bump rows (simulation, SMASH x V3, owned rows; fractions of the weight, in per cent, +- MC)")
L.append("columns: C2 pre all M4 | C3 pre all M4 | C2+C3 pre all M4 | C2+C3 pre window B | C2+C3 pre near the edge (-15..+5 MeV) | C2+C3 S0 all M4")
keys = ["pre_C2_all", "pre_C3_all", "pre_C23_all", "pre_C23_B", "pre_C23_edge15", "S0_C23_all"]
L.append("rows (weight) | " + " | ".join("%d (%.0f)" % (fr[k]["Nrows"], fr[k]["W"]) for k in keys))
lab = [("geo", "photon hits an active ECAL module (conditional completion)"), ("geo_6sec", "  the same with all six sectors and no dead module"),
       ("geo_strict", "  the same with the 92 mm faces only (no 2 mm gaps)"), ("geo_flat", "  the same with a flat weight along the circle"),
       ("geo_kloe", "  the same with the KLOE-2 Dalitz plot along the circle"),
       ("geo_E80", "  and true E_gamma > 80 MeV"), ("geo_E100", "  and true E_gamma > 100 MeV"), ("geo_E130", "  and true E_gamma > 130 MeV"),
       ("found", "photon found by the ECAL simulation (cluster truth, any quality)"), ("found_quality", "  and the cluster passes ecal_quality (PID 1, 0.8 < beta < 1.2, E > 100 MeV)"),
       ("tag_any_m5c", "row tagged (any photon), m5c 0.49-0.62"), ("tag_true_m5c", "row tagged with the TRUE photon, m5c 0.49-0.62 (fully reconstructed)"),
       ("tagTrueL3_m5c", "  the same with the cluster energies x 1.285 (data-like calibration)")]
for q, t in lab:
    L.append("%-92s | " % t + " | ".join(pc(fr[k][q]) for k in keys))
L.append("efficiency inside the active modules (pre, C2+C3, all M4): found %.3f, found + quality %.3f, fully reconstructed %.3f" % (
    num["eff_in_acc_pre"]["found"], num["eff_in_acc_pre"]["found_quality"], num["eff_in_acc_pre"]["tag_true"]))
L.append("photon polar angle of all accepted C2+C3 rows (pre): theta < 11.5 deg %.3f, 11.5-47.5 deg %.3f, > 47.5 deg %.3f; in 18-85 deg %.3f" % (
    num["theta_shares"]["below11p5"], num["theta_shares"]["band"], num["theta_shares"]["above47p5"], num["theta_shares"]["hades18_85"]))
L.append("true photon energy: mean %.0f MeV (all), %.0f MeV (in the modules); share below 100 MeV: %.3f (all), %.3f (in the modules)" % (
    num["energy_shares"]["mean_all"], num["energy_shares"]["mean_acc"], num["energy_shares"]["all_below100"], num["energy_shares"]["acc_below100"]))
L.append("")
L.append("versus M(pi+pi-) (pre, C2+C3): slice | weight | in modules | and E > 100 MeV | found + quality | fully reconstructed | the same x 1.285 | <E_gamma> all / in modules [MeV]")
for r in num["vs_mpipi"]:
    c = r["C23"]
    L.append("%.2f-%.2f | %.0f | %s | %s | %s | %s | %s | %.0f / %.0f" % (r["lo"], r["hi"], c["W"], pc(c["geo"]), pc(c["geo_E100"]), pc(c["found_quality"]),
                                                                pc(c["tag_true_m5c"]), pc(c["tagTrueL3_m5c"]), c["eMean"][0], c["eMeanAcc"][0]))
L.append("versus the distance to the edge M4' - E(M(pipi), M(ee)) [GeV/c2] (pre, C2+C3): the same columns")
for r in num["vs_dist"]:
    c = r["C23"]
    L.append("%+.3f..%+.3f | %.0f | %s | %s | %s | %s | %s | %.0f / %.0f" % (r["lo"], r["hi"], c["W"], pc(c["geo"]), pc(c["geo_E100"]), pc(c["found_quality"]),
                                                                    pc(c["tag_true_m5c"]), pc(c["tagTrueL3_m5c"]), c["eMean"][0], c["eMeanAcc"][0]))
L.append("")
cl = num["closure"]
L.append("closure: %d photons found by the ECAL simulation; inside the predicted modules %.4f, in a dead module %.4f, outside %.4f; "
         "residuals (5/50/95 %%) M(e+e-gamma) - m_pi0 %s MeV, M(5) - m_eta %s MeV, cos(alpha) %s; E_cluster/E_true (16/50/84 %%) %s (with quality %s)" % (
             cl["n_found"], cl["inside"], cl["deadmodule"], cl["outside"], " ".join("%.4f" % x for x in cl["res_pi0_q"]),
             " ".join("%.4f" % x for x in cl["res_eta_q"]), " ".join("%.1e" % x for x in cl["cos_q"]), " ".join("%.3f" % x for x in cl["ecl_over_etrue_q"]),
             " ".join("%.3f" % x for x in cl["ecl_over_etrue_q_quality"])))
L.append("calibration of the circle prediction (pre, C2+C3; rows binned in the predicted probability): " + "; ".join(
    "%.2f-%.2f: N %d, found %.3f +- %.3f, predicted %.3f (flat %.3f), found/predicted %.2f +- %.2f" % (c["lo"], c["hi"], c["N"], c["found"][0], c["found"][1], c["pred"], c["predFlat"], c["ratio"], c["eratio"])
    for c in num["closure"]["calib"]))
L.append("C2 gamma* angle factor: all C2 rows (circle average) %.4f (E1 %.4f, E2 %.4f); near the edge %.4f; tagged C2 rows (pre) %.4f +- %.4f (bedge recommended 0.8855)" % (
    *num["c_qed_circle_all"], fr["pre_C2_edge15"]["c_qed_circle"][0], *num["c_tag_pre"]))
open(T + "/ecalacc_results.txt", "w").write("\n".join(L) + "\n")
print("\n".join(L))

# tag table
L = ["ecalacc: ECAL-tagged eta, data against the simulation at the measured eta scale (frozen fit 0.5463 / 6.2 MeV, pol2, 0.49-0.62)",
     "T_exp = f_eta(edge) x T_sim(frozen fit) x [composition] x [calibration]; errors: stat, syst of f, MC of f (+) MC of T_sim (%.1f %%, WP6)" % (100 * T_SIM_MC_REL)]
L.append("sel slice | data T | sim T (frozen) | f_tag = T/T_sim | expected: as is (f 0.392) | C2 corr. (0.418) | C2 corr. + calib. x 1.285 | kappa 0.48 (0.560) | data / expected (C2 corr.) | data / expected (C2 corr. + calib.)")
for s in ("pre", "S0"):
    for j in range(5):
        k = "%s_%s" % (s, tsn[j])
        t = num["tagfits"][k]
        e0 = ex["%s_%s_asis_L0" % (s, tsn[j])]; e1 = ex["%s_%s_c2_L0" % (s, tsn[j])]; e3 = ex["%s_%s_c2_L3" % (s, tsn[j])]; e4 = ex["%s_%s_k48_L0" % (s, tsn[j])]
        def fe(e):
            return "%.0f +- %.0f +- %.0f +- %.0f" % (e["value"], e["stat"], e["syst"], e["mc"])
        def rr(e):
            tot = math.sqrt((t["edT"] / t["dT"]) ** 2 + ((e["stat"] ** 2 + e["syst"] ** 2 + e["mc"] ** 2) / e["value"] ** 2)) if t["dT"] > 0 and e["value"] > 0 else 0
            return "%.3f +- %.3f" % (e["ratio"], e["ratio"] * tot)
        L.append("%s %s | %.0f +- %.0f | %.0f | %.3f +- %.3f | %s | %s | %s | %s | %s | %s" % (s, tsn[j], t["dT"], t["edT"], t["sT"], t["dT"] / t["sT"] if t["sT"] else 0,
                                                                         t["edT"] / t["sT"] if t["sT"] else 0, fe(e0), fe(e1), fe(e3), fe(e4), rr(e1), rr(e3)))
st = num.get("start")
if st:
    L.append("data, events without START (evtStart != 1): T = %.0f +- %.0f of %.0f; tag rate per e+e- row relative to events with START: %.2f +- %.2f (all M4), %.2f +- %.2f (window B); "
             "if they had the START rate: +%.0f rows (+%.1f %%)" % (st["T0"], st["eT0"], st["T0"] + st["T1"], st["ratio_rate_all"], st["ratio_rate_all"] * st["ratio_err_rel"],
                                                                st["ratio_rate_B"], st["ratio_rate_B"] * st["ratio_err_rel"], st["loss_if_start0_like_start1"],
                                                                100 * st["loss_if_start0_like_start1"] / (st["T0"] + st["T1"])))
L.append("calibration factor of the simulated tag (pre; frozen fit of T1+T2 with the C2+C3 rows re-tagged at cluster energies x lambda; truth ratio in brackets): " +
         ", ".join("lambda %.3f: x %.3f (%.3f)" % (lam[l], ex["pre_all_c2_L%d" % l]["lamFac"], ex["pre_all_c2_L%d" % l]["lamFacTruth"]) for l in range(len(lam))))
pd, ps0, ps3 = num["pi0"]["data"], num["pi0"]["sim_L0"], num["pi0"]["simrec_L3"]
L.append("pi0 in M(e+e-gamma): data %.4f / %.1f MeV (window 0.119-0.173 keeps %.3f); simulation as is %.4f / %.1f MeV (%.3f); simulation x 1.285 %.4f / %.1f MeV (%.3f); "
         "resolution factor data/sim(x 1.285) = %.3f" % (pd["mean"], 1e3 * pd["sigma"], pd["window_acceptance"], ps0["mean"], 1e3 * ps0["sigma"], ps0["window_acceptance"],
                                                          ps3["mean"], 1e3 * ps3["sigma"], ps3["window_acceptance"], RESFAC))
if st:
    startf = 1.0 - st["loss_if_start0_like_start1"] / (st["T0"] + st["T1"] + st["loss_if_start0_like_start1"])
    for sname in ("pre", "S0"):
        e3 = ex["%s_all_c2_L3" % sname]; t = num["tagfits"]["%s_all" % sname]
        best = e3["value"] * RESFAC * startf
        tot = math.sqrt(e3["stat"] ** 2 + e3["syst"] ** 2 + e3["mc"] ** 2) * RESFAC * startf
        num["best_%s" % sname] = dict(value=best, stat=e3["stat"] * RESFAC * startf, syst=e3["syst"] * RESFAC * startf, mc=e3["mc"] * RESFAC * startf,
                                       ratio=t["dT"] / best, eratio=t["dT"] / best * math.sqrt((t["edT"] / t["dT"]) ** 2 + (tot / best) ** 2), startf=startf)
        L.append("best estimate %s: f(C2 corr.) x T_sim(x 1.285) x resolution %.3f x START %.3f = %.0f +- %.0f (stat) +- %.0f (syst) +- %.0f (MC); data %.0f +- %.0f; data/expected %.3f +- %.3f" % (
            sname, RESFAC, startf, best, num["best_%s" % sname]["stat"], num["best_%s" % sname]["syst"], num["best_%s" % sname]["mc"], t["dT"], t["edT"],
            num["best_%s" % sname]["ratio"], num["best_%s" % sname]["eratio"]))
json.dump(num, open(O + "/ecalacc_numbers.json", "w"), indent=1)
open(T + "/ecalacc_tag.txt", "w").write("\n".join(L) + "\n")
print("\n".join(L))

# LaTeX tables
tex = [r"\begin{tabular}{@{}l r r r r r@{}}", r"\toprule",
       r" & C2 & C3 & C2+C3 & C2+C3, $B$ & C2+C3, edge \\", r"\midrule"]
rows_tex = [("geo", r"photon in an active ECAL module"), ("geo_6sec", r"\quad six full sectors"), ("geo_E100", r"\quad and $E_\gamma^{\mathrm{true}}>100$\,MeV"),
            ("found_quality", r"found by the ECAL simulation (quality cuts)"), ("tag_true_m5c", r"fully reconstructed (true photon, $m_{5c}$ window)"),
            ("tagTrueL3_m5c", r"\quad with data-like energy scale ($\times1.285$)")]
for q, t in rows_tex:
    tex.append(t + " & " + " & ".join("%.1f" % (100 * fr[k][q][0]) for k in ["pre_C2_all", "pre_C3_all", "pre_C23_all", "pre_C23_B", "pre_C23_edge15"]) + r" \\")
tex.append(r"\bottomrule"); tex.append(r"\end{tabular}")
open(T + "/ecalacc_results.tex", "w").write("\n".join(tex) + "\n")
tex = [r"\begin{tabular}{@{}l r r r r r@{}}", r"\toprule",
       r"$M(\pi^+\pi^-)$ [GeV/$c^2$] & $T$ data & $T_{\mathrm{sim}}$ & $f_{\mathrm{tag}}$ & $T_{\mathrm{exp}}$ (C2 corr.) & $T_{\mathrm{exp}}$ (C2 corr., $\times1.285$) \\", r"\midrule"]
for j in range(5):
    k = "pre_%s" % tsn[j]
    t = num["tagfits"][k]; e1 = ex["pre_%s_c2_L0" % tsn[j]]; e3 = ex["pre_%s_c2_L3" % tsn[j]]
    tex.append(r"%s & $%.0f\pm%.0f$ & %.0f & $%.3f\pm%.3f$ & $%.0f\pm%.0f$ & $%.0f\pm%.0f$ \\" % (
        tsn[j].replace(">=", r"$\ge$").replace("<", "$<$"), t["dT"], t["edT"], t["sT"], t["dT"] / t["sT"] if t["sT"] else 0, t["edT"] / t["sT"] if t["sT"] else 0,
        e1["value"], math.sqrt(e1["stat"] ** 2 + e1["syst"] ** 2 + e1["mc"] ** 2), e3["value"], math.sqrt(e3["stat"] ** 2 + e3["syst"] ** 2 + e3["mc"] ** 2)))
tex.append(r"\bottomrule"); tex.append(r"\end{tabular}")
open(T + "/ecalacc_tag.tex", "w").write("\n".join(tex) + "\n")
