#!/usr/bin/env python3
# ============================================================================
# edgeu_ana.py  --  package edgeu: the bump in u = M4 - E(M(pipi), M(ee)): kinematics checks, spectra, global fit,
#                   errors (bootstrap, replicas), systematic variants, checks, conversions
# ============================================================================
# Written 2026-11-15 (prefix edgeu). Reads the columns of edgeu_dump.C; writes
#   /hdd3/przygoda/HADES/FAT_out/all_report/edgeu/edgeu_ana.json   (every number; read by edgeu_tables.py, edgeu_figs.py)
#   /hdd3/przygoda/HADES/FAT_out/all_report/edgeu/edgeu_ana.log    (readable log)
#   /hdd3/przygoda/HADES/FAT_out/all_report/edgeu/edgeu_figdata.npz (spectra for the figures)
#
# Estimator (nominal):
#   data   S0, cor_ok == 1, isBest >= 0, 0.30 <= M(pipi) < 0.40, M(ee) < m_pi0 (u is defined only there; no bump row is lost),
#          u = M4 - E(M(pipi), M(ee)) (defs.h Edge, data scale 0.99274 / 0.9970), 2.5 MeV bins with an edge at u = 0;
#          N - CB, CB = 2 sqrt(N++ N--) in 10 MeV cells of u, interpolated (ScanEngine rule), per M(pipi) slice
#          (0.30-0.34-0.37-0.39-0.40, the bedge slices), k = 1.
#   model  N - CB - a_P K(u) = pol2(u) + h T(u), range -45 <= u < +45 MeV; K = simulated peak C4 + C5 in u, a_P = 0.418249
#          (fixed, the measured peak scale); T = bin average of 0.5 erfc(u/(sqrt2 sigma)), sigma = Gaussian core of the
#          simulated u resolution at the edge; V = cov(CB) + diag(max(model of N, 1)), iterated (scan::FitEdgeStep logic).
#   sim    h_sim = the same linear estimator (final data covariance) applied to the simulated bump C2 + C3 in u (mode A,
#          V3, C2 x 0.879 / 0.8825 in E1 / E2 for the gamma* angle, D35; M4' = 0.99615 M4, M(pipi)' = 0.9952 M(pipi)).
#   f_eta(u-edge) = h_data / h_sim.
#   stat   200 bootstrap replicas of the data (Poisson(1) per event group and tree, own seed), each refitted;
#   MC     100 simulation replicas (WP0 generator and seed), bump and peak templates, estimator of the nominal fit.
#
# Run (from /home/przygoda/HADES/FAT/FAT; about 6 GB, about 10 minutes):
#   python3 -I plots/all_report/edgeu_ana.py
# ============================================================================
import json
import math
import sys
import time

import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import chi2 as chi2dist

sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import edgeu_core as C

T0 = time.time()
LOG = open(C.OUT + "edgeu_ana.log", "w")
def P(*a):
    s = " ".join(str(x) for x in a)
    print(s); LOG.write(s + "\n"); LOG.flush()
RES = {}

D = C.load_data()
S = C.load_sim()
P("# edgeu_ana.py, 2026-11-15: data rows %d, simulated rows %d (%.0f s)" % (D["n"], S["n"], time.time() - T0))
PER = C.period_of_day(D["day"])

# =====================================================================================================================
# 0. constants measured here: the simulated u resolution at the edge, the simulated M(pipi) reco/true scale
# =====================================================================================================================
w_nom_raw = C.sim_weight(S, c2corr=False)
u_s0, mpp_s0, m4p_s0 = C.u_sim(S)
ut_s = C.EDGE_SCALE * (S["m4t"].astype(np.float64) - C.m4max_true(S["mpipit"].astype(np.float64), S["meet"].astype(np.float64)))
gS = S["g"]
def reso_sample(u, mpp, w, near=-0.015):
    m = (gS <= 1) & np.isfinite(u) & np.isfinite(ut_s) & (w > 0) & (mpp >= 0.30) & (mpp < 0.40) & (ut_s > near)
    return (u - ut_s)[m], w[m]
def gauss_core(d, ww):
    h, e = np.histogram(d, bins=np.arange(-0.040, 0.0401, 0.002), weights=ww)
    h2, _ = np.histogram(d, bins=e, weights=ww * ww)
    x = 0.5 * (e[1:] + e[:-1])
    f = lambda x, a, mu, s: a * np.exp(-0.5 * (x - mu) ** 2 / s ** 2)
    p, cv = curve_fit(f, x, h, p0=[h.max(), 0, 0.007], sigma=np.sqrt(np.maximum(h2, 1.0)), absolute_sigma=True)
    return float(p[1]), float(abs(p[2])), float(math.sqrt(cv[2, 2]))
dS, wdS = reso_sample(u_s0, mpp_s0, w_nom_raw)
mu_c, SIG_U, eSIG_U = gauss_core(dS, wdS)
q = np.percentile(dS, [15.865, 50, 84.135])
mean_d = float(np.average(dS, weights=wdS)); rms_d = float(np.sqrt(np.average((dS - mean_d) ** 2, weights=wdS)))
RPEAK = 5.44 / 6.09   # data / simulated eta-peak width (ScanEngine 5.44 +- 0.60 MeV; WP6/bedge 6.1 MeV)
SIG_DATA = RPEAK * SIG_U
P("0. simulated u resolution at the edge (C2 + C3, 0.30 <= M(pipi)' < 0.40, u_true > -15 MeV; %d rows, sum w %.0f):" % (len(dS), wdS.sum()))
P("   Gaussian core sigma %.2f +- %.2f MeV (mean %.2f); half of the central 68 %% %.2f MeV; rms %.2f MeV; median %.2f MeV" %
  (1e3 * SIG_U, 1e3 * eSIG_U, 1e3 * mu_c, 1e3 * (q[2] - q[0]) / 2, 1e3 * rms_d, 1e3 * q[1]))
P("   data-like smearing (x %.3f, the data/sim eta-peak width ratio): %.2f MeV" % (RPEAK, 1e3 * SIG_DATA))
m = (gS <= 1) & (S["mpipit"] > 0)
CPP = float(np.median(S["mpipi"][m] / S["mpipit"][m]))
C.CPP_SIM = CPP
P("   simulated bump M(pipi) reco/true median %.5f (K0S: 1.0018); M4 reco/true median %.5f (peak transfer: 0.99274/0.99615 = %.5f)" %
  (CPP, float(np.median(S["m4"][(gS <= 1) & (S["m4t"] > 0)] / S["m4t"][(gS <= 1) & (S["m4t"] > 0)])), C.EDGE_SCALE / C.S_NOM))
# decomposition of the u resolution
mm = (gS <= 1) & np.isfinite(u_s0) & np.isfinite(ut_s) & (w_nom_raw > 0) & (mpp_s0 >= 0.30) & (mpp_s0 < 0.40) & (ut_s > -0.015)
dm4 = (m4p_s0 - C.EDGE_SCALE * S["m4t"])[mm]; dE = (u_s0 - ut_s)[mm] - dm4
RESO = dict(sigma_core=SIG_U, e_sigma_core=eSIG_U, mean_core=mu_c, half68=float((q[2] - q[0]) / 2), rms=rms_d, median=float(q[1]),
            rows=int(len(dS)), sumw=float(wdS.sum()), sig_data=SIG_DATA, rpeak=RPEAK, cpp=CPP,
            rms_dm4=float(np.std(dm4)), rms_dE=float(np.std(dE)), corr=float(np.corrcoef(dm4, dE)[0, 1]))
P("   parts: M4' error rms %.2f MeV, edge-function error rms %.2f MeV, correlation %.2f" % (1e3 * RESO["rms_dm4"], 1e3 * RESO["rms_dE"], RESO["corr"]))
RES["reso"] = RESO

# =====================================================================================================================
# 1. kinematics with the simulated truth
# =====================================================================================================================
P("\n1. KINEMATICS (simulated truth)")
K = {}
K["edge_checks"] = {k: float(v) for k, v in zip(["E(0.30,0)", "E(0.39,0)", "E(0.39,60MeV)", "E(0.4116,0)", "E(0.30,100MeV)"],
                    C.edge(np.array([0.30, 0.39, 0.39, 0.4116, 0.30]), np.array([0, 0, 0.06, 0, 0.1])))}
P("   frozen edge function:", K["edge_checks"], "(report: 0.51859, 0.50077, 0.50958, 0.47390, 0.5326)")
mt4 = C.m4max_true(S["mpipit"].astype(np.float64), S["meet"].astype(np.float64))
mt40 = C.m4max_true(S["mpipit"].astype(np.float64), np.zeros(S["n"]))
for gg in (0, 1):
    m = (gS == gg) & (S["m4t"] > 0) & (S["mpipit"] > 0)
    d = (S["m4t"] - mt4)[m]; d0 = (S["m4t"] - mt40)[m]; ww = w_nom_raw[m]
    K["truth_%s" % C.GROUPS[gg]] = dict(rows=int(m.sum()), max_excess=float(np.nanmax(d)), frac_above_1e4=float(np.mean(d > 1e-4)),
                                        frac_above_q0_1MeV=float(np.average(d0 > 0.001, weights=ww + 1e-30)),
                                        frac_above_q0_5MeV=float(np.average(d0 > 0.005, weights=ww + 1e-30)),
                                        mee_true_q=[float(x) for x in np.percentile(S["meet"][m], [5, 25, 50, 75, 95])])
    P("   %s: %d rows; max(m4_true - M4max(M(pipi)_true, M(ee)_true)) = %+.3f MeV (never above); with q = 0: %.2f %% above by > 1 MeV, %.2f %% by > 5 MeV (weighted)" %
      (C.GROUPS[gg], m.sum(), 1e3 * np.nanmax(d), 100 * K["truth_%s" % C.GROUPS[gg]]["frac_above_q0_1MeV"], 100 * K["truth_%s" % C.GROUPS[gg]]["frac_above_q0_5MeV"]))
# simulation scale: u on the data scale (transfer + data edge) vs the simulation-scale edge on the raw simulated masses
u_rawscale = S["m4"].astype(np.float64) - (C.EDGE_SCALE / C.S_NOM) * C.m4max_true(S["mpipi"].astype(np.float64) * C.SPP_NOM / C.EDGE_PIPI, S["mee"].astype(np.float64))
m = np.isfinite(u_s0) & np.isfinite(u_rawscale) & (S["m4t"] > 0)
K["scale_equiv_maxdiff"] = float(np.max(np.abs(u_s0[m] - C.S_NOM * u_rawscale[m])))
P("   data-scale u (M4' = 0.99615 M4, M(pipi)' = 0.9952 M(pipi), data edge) = 0.99615 x u(raw simulated masses, simulation-scale edge %.5f M4max(M(pipi)/%.5f)): max |diff| %.2e GeV" %
  (C.EDGE_SCALE / C.S_NOM, C.EDGE_PIPI / C.SPP_NOM, K["scale_equiv_maxdiff"]))
# step sharpness: the spread of the edge inside the bedge slices vs the resolution
spread = {}
for a, b in [(0.30, 0.34), (0.34, 0.37), (0.37, 0.39), (0.39, 0.40)]:
    spread["%.2f-%.2f" % (a, b)] = float(1e3 * (C.edge(np.array([a]), np.array([0.]))[0] - C.edge(np.array([b]), np.array([0.]))[0]))
K["edge_spread_in_slices_MeV"] = spread
P("   edge spread inside the bedge slices (E(lo) - E(hi), MeV):", spread)
# reconstructed u of C2 and C3: fraction above the edge
for lab, uq in [("E(m,q)", True), ("E(m,0)", False)]:
    uu, mp_, _ = C.u_sim(S, useq=uq)
    for gg in (0, 1):
        for cl, cm in [("all", np.ones(S["n"], bool)), ("E1", S["mee"] < 0.010), ("E2", S["mee"] >= 0.010), ("E20", S["mee"] >= 0.020)]:
            m = (gS == gg) & cm & np.isfinite(uu) & (w_nom_raw > 0) & (mp_ >= 0.30) & (mp_ < 0.40)
            if w_nom_raw[m].sum() <= 0:
                continue
            fr10 = float(w_nom_raw[m & (uu > 0.010)].sum() / w_nom_raw[m].sum()); fr20 = float(w_nom_raw[m & (uu > 0.020)].sum() / w_nom_raw[m].sum())
            K["above_%s_%s_%s" % (C.GROUPS[gg], cl, "q" if uq else "q0")] = [fr10, fr20]
            P("   reco %s %s %s: fraction with u > +10 MeV %.4f, > +20 MeV %.4f" % (C.GROUPS[gg], cl, lab, fr10, fr20))
RES["kin"] = K

# =====================================================================================================================
# 2. machinery: configurations, caches, evaluation
# =====================================================================================================================
NOM = dict(sel="S0", best=False, period=None, fold=None, dscale=C.EDGE_SCALE, dspp=C.EDGE_PIPI, useq=True, qcut=True, mshift=0.0,
           grid=(-0.200, 0.0025, 140, 4), pplo=0.30, pphi=0.40, sledges=(0.30, 0.34, 0.37, 0.39, 0.40), cbslice=True,
           cbk=1.0, cbmode=0, clsrule=0, cbmap=False, cls="all",
           s=C.S_NOM, r=C.R_NOM, spp=C.SPP_NOM, rpp=None, code=108, mode="A", ssel=None, sbest=None, c2corr=True,
           sscale=C.EDGE_SCALE, sepp=C.EDGE_PIPI, suseq=None,
           lo=-0.045, hi=0.045, order=2, sigma=SIG_U, sigma_sim=None, tmpl="gauss", simfit="datacov",
           peak="fixed", aP=C.AP, peak_s=None, peak_r=None, simtarget="bump", simscale=1.0, cfE1=C.CF_E1, cfE2=C.CF_E2)
def cfg(**kw):
    c = dict(NOM); c.update(kw)
    return c
def grid_of(c):
    lo, bw, nb, reb = c["grid"]
    return C.Grid(lo, bw, nb, reb)

_UD, _US, _HD, _HS, _HDR, _HSR = {}, {}, {}, {}, {}, {}
def ud_key(c):
    return (c["dscale"], c["dspp"], c["useq"], c["qcut"], c["mshift"])
def U_D(c):
    k = ud_key(c)
    if k not in _UD:
        _UD[k] = C.u_data(D, c["dscale"], c["dspp"], c["useq"], c["qcut"], c["mshift"])
    return _UD[k]
def us_key(c, s=None, r=None):
    suseq = c["useq"] if c["suseq"] is None else c["suseq"]
    return (c["s"] if s is None else s, c["r"] if r is None else r, c["spp"], c["rpp"], c["sscale"], c["sepp"], suseq, c["qcut"], c["mshift"])
def U_S(key):
    if key not in _US:
        if len(_US) > 6:
            _US.pop(next(iter(_US)))
        s, r, spp, rpp, sc, epp, uq, qc, ms = key
        _US[key] = C.u_sim(S, s, r, spp, rpp, sc, epp, uq, qc, None, ms)
    return _US[key]
def dmask(c):
    b = D["bits"]
    m = (b & (1 if c["sel"] == "S0" else 2)) != 0
    if c["best"]:
        m &= (b & 4) != 0
    if c["period"] is not None:
        m &= PER == c["period"]
    if c["fold"] is not None:
        m &= D["fold"] == c["fold"]
    return m
def hd_key(c):
    return (ud_key(c), c["sel"], c["best"], c["period"], c["fold"], c["grid"])
def H_D(c):
    k = hd_key(c)
    if k not in _HD:
        _HD[k] = C.hist_data(D, U_D(c), dmask(c), grid_of(c))
    return _HD[k]
# bootstrap weights for all data rows (one matrix: the same replica b uses the same group weights in every variant)
NB_D = 200
BW_D = C.data_boot_weights(D, np.arange(D["n"]), NB_D)
P("   data bootstrap weights: %d replicas x %d rows (%.0f s)" % (NB_D, D["n"], time.time() - T0))
def H_DR(c):
    k = hd_key(c)
    if k not in _HDR:
        if len(_HDR) > 2:
            _HDR.pop(next(iter(_HDR)))
        u = U_D(c); g = grid_of(c)
        pb = C.ppbin_of(D["mpipi"]); bb = g.bin(u)
        rows = np.where(dmask(c) & (bb >= 0) & (pb >= 0))[0]
        Dsub = {kk: D[kk][rows] for kk in ("k", "mpipi", "mee")}
        _HDR[k] = C.hist_data(Dsub, u[rows], np.ones(len(rows), bool), g, weights=BW_D[:, rows])
    return _HDR[k]
NR_S = 100
ROWS_EP = np.where(gS <= 3)[0]                       # bump and peak classes
BW_S_EP = C.sim_replica_weights(S, ROWS_EP, NR_S)
P("   simulation replica weights (WP0 generator, seed 0x5157A11): %d replicas x %d bump/peak rows (%.0f s)" % (NR_S, len(ROWS_EP), time.time() - T0))
def sim_w(c, best=None):
    ssel = c["sel"] if c["ssel"] is None else c["ssel"]
    sb = c["best"] if c["sbest"] is None else c["sbest"]
    return C.sim_weight(S, c["code"], ssel, c["mode"], sb, c2corr=False) * c["simscale"]
def hs_key(c, s=None, r=None):
    ssel = c["sel"] if c["ssel"] is None else c["ssel"]
    sb = c["best"] if c["sbest"] is None else c["sbest"]
    return (us_key(c, s, r), c["code"], c["mode"], ssel, sb, c["simscale"], c["grid"])
def H_S(c, s=None, r=None):
    k = hs_key(c, s, r)
    if k not in _HS:
        if len(_HS) > 12:
            _HS.pop(next(iter(_HS)))
        u, mpp, _ = U_S(us_key(c, s, r))
        _HS[k] = C.hist_sim(S, u, mpp, sim_w(c), grid_of(c))
    return _HS[k]
def H_SR(c, s=None, r=None):
    k = hs_key(c, s, r)
    if k not in _HSR:
        if len(_HSR) > 3:
            _HSR.pop(next(iter(_HSR)))
        u, mpp, _ = U_S(us_key(c, s, r))
        w = sim_w(c)
        Ssub = {kk: S[kk][ROWS_EP] for kk in ("g", "mee")}
        _HSR[k] = C.hist_sim(Ssub, u[ROWS_EP], mpp[ROWS_EP], w[ROWS_EP], grid_of(c), reps=BW_S_EP)
    return _HSR[k]
def slices(c):
    e = list(c["sledges"])
    return C.slices_of(c["pplo"], c["pphi"], e)
def c2fac(c, cls):
    """C2 gamma* factor per M(ee) class index (0 E1, 1-2 E2)"""
    return np.array([c["cfE1"], c["cfE2"], c["cfE2"]]) if c["c2corr"] else np.ones(3)
def simpart(Hs, c, groups, cls=None, c2=True):
    """sum of simulated groups over the slices and classes of c; C2 x factor per class (c2 = True)"""
    cs = C.CLS_SETS[c["cls"] if cls is None else cls]
    y = np.zeros(Hs.shape[-1])
    fac = c2fac(c, cls)
    for gname in groups:
        gi = C.G[gname]
        for (a, b) in slices(c):
            for k in cs:
                y = y + Hs[gi, a:b, k].sum(0) * (fac[k] if (gname == "C2" and c2) else 1.0)
    return y
def tmpl(c, g, sim=False):
    sg = c["sigma"] if (not sim or c["sigma_sim"] is None) else c["sigma_sim"]
    if c["tmpl"] == "gauss":
        return C.step_gauss(g, sg)
    if c["tmpl"] == "centre":
        from scipy.special import erfc
        return 0.5 * erfc(g.x / (math.sqrt(2) * sg))
    if c["tmpl"] == "emp":
        u, mpp, _ = U_S(us_key(c))
        d, wd = reso_sample(u, mpp, sim_w(c))
        return C.step_emp(g, d, wd)
    raise ValueError(c["tmpl"])
_LSR = {}
def cb_m4map(c, wls=None):
    """variant: the CB of the frozen M4 rule (2 sqrt(N++ N--) in 10 MeV M4 cells, per slice) distributed over u with the
    like-sign u shape of each 50 MeV M4 band (bands from 0.30). Returns (cb, cov) on the u grid of c; wls: bootstrap
    weights of all data rows (None = 1)."""
    g = grid_of(c)
    key = (hd_key(c), tuple(slices(c)))
    if key not in _LSR:
        u = U_D(c); m = dmask(c) & (D["k"] >= 1)
        pb = C.ppbin_of(D["mpipi"]); bu = g.bin(u)
        b4 = np.floor((D["m4"].astype(np.float64) - 0.30) / 0.0025).astype(np.int64)
        b4[(b4 < 0) | (b4 >= 220)] = -1
        lst = []
        for (a, b) in slices(c):
            rows = np.where(m & (pb >= a) & (pb < b) & (b4 >= 0))[0]
            lst.append((rows, b4[rows], bu[rows], D["k"][rows]))
        _LSR[key] = lst
    gm4 = C.Grid(0.30, 0.0025, 220, 4)   # M4 0.30-0.85
    band = np.arange(gm4.nb) // 20
    cb = np.zeros(g.nb); cov = np.zeros((g.nb, g.nb))
    for rows, b4, bu, kk in _LSR[key]:
        wgt = np.ones(len(rows)) if wls is None else wls[rows].astype(np.float64)
        pp = np.bincount(b4[kk == 1], weights=wgt[kk == 1], minlength=gm4.nb)
        mm_ = np.bincount(b4[kk == 2], weights=wgt[kk == 2], minlength=gm4.nb)
        cbm, covm = C.cb_rule(pp, mm_, gm4, c["cbk"], c["cbmode"])
        okb = bu >= 0
        for bnd in range(11):
            inb = (b4 // 20 == bnd)
            n_b = wgt[inb].sum()
            if n_b <= 0:
                continue
            shape = np.bincount(bu[inb & okb], weights=wgt[inb & okb], minlength=g.nb) / n_b
            sel = band == bnd
            cbb = cbm[sel].sum(); vb = covm[np.ix_(sel, sel)].sum()
            cb += cbb * shape; cov += vb * np.outer(shape, shape)
    return cb, cov
def dspec(c, H=None, cbo=None):
    g = grid_of(c)
    H = H_D(c) if H is None else H
    sp = C.data_spec(H, slices(c), c["cls"], g, c["cbk"], c["cbmode"], c["clsrule"], cbslice=c["cbslice"])
    if cbo is not None:
        sp = C.Spec(sp.N, cbo[0], cbo[1], sp.pp, sp.mm)
    return sp

def evaluate(c, boot=True, rep=True, keep=False):
    """the estimator for configuration c: f = h_data / h_sim, its bootstrap (stat) and replica (MC) distributions,
    the linear decomposition of h_sim; returns a dict"""
    g = grid_of(c)
    T = tmpl(c, g); Ts = tmpl(c, g, sim=True)
    Hs = H_S(c)
    Hp = Hs if (c["peak_s"] is None and c["peak_r"] is None) else H_S(c, c["peak_s"], c["peak_r"])
    peakT = simpart(Hp, c, ["C4", "C5"])
    cbo = cb_m4map(c) if c["cbmap"] else None
    sp = dspec(c, cbo=cbo)
    fixed = c["aP"] * peakT if c["peak"] == "fixed" else None
    free = [peakT] if c["peak"] == "free" else None
    R = C.fit_step(sp, T, g, c["lo"], c["hi"], c["order"], fixed=fixed, free_cols=free)
    if c["simtarget"] == "bump":
        ysim = simpart(Hs, c, ["C2", "C3"])
    else:   # T12: everything simulated; the simulated peak is subtracted with its own scale 1
        ysim = simpart(Hs, c, ["C2", "C3", "Geta", "GDom", "GDnom", "GCom", "GCnom"])
    if c["simfit"] == "datacov" and np.allclose(Ts, T):
        hs = C.h_lin(R, ysim)
        Rs = None
    else:
        Rs = C.fit_step(None, Ts, g, c["lo"], c["hi"], c["order"], own_var=ysim) if c["simfit"] == "own" else \
            C.fit_step(None, Ts, g, c["lo"], c["hi"], c["order"], own_var=ysim, Vfix=R.V)
        hs = Rs.h
    out = dict(hd=R.h, ehd=R.eh, hs=hs, f=R.h / hs, statfit=R.eh / hs, chi2=R.chi2, chi20=R.chi20, ndf=R.ndf,
               aPfit=(float(R.p[c["order"] + 1]) if c["peak"] == "free" else None))
    # linear decomposition (data covariance)
    lin = {}
    for gname in ["C2", "C3", "C4", "C5", "Geta", "GDom", "GDnom", "GCom", "GCnom"]:
        lin[gname] = C.h_lin(R, simpart(Hs, c, [gname], c2=False))
    lin["C2c"] = C.h_lin(R, simpart(Hs, c, ["C2"], c2=True))
    out["lin"] = lin
    cfall = lin["C2c"] / lin["C2"] if lin["C2"] != 0 else 1.0
    out["f_unc"] = R.h / (lin["C2"] + lin["C3"]) if c["simtarget"] == "bump" else None
    out["a_k48"] = R.h / (lin["C2c"] + 0.48 * lin["C3"])
    out["a_k48_unc"] = R.h / (lin["C2"] + 0.48 * lin["C3"])
    out["c2share"] = lin["C2c"] / (lin["C2c"] + lin["C3"])
    out["cfall"] = cfall
    if keep:
        out["R"] = R; out["sp"] = sp; out["ysim"] = ysim; out["peakT"] = peakT; out["T"] = T; out["Hs"] = Hs
    # ---- data bootstrap: refit every replica, the same estimator on the simulation
    if boot:
        HR = H_DR(c)
        fb, fub, ak, hdb = [], [], [], []
        for b in range(NB_D):
            cbo_b = None
            if c["cbmap"]:
                cbo_b = cb_m4map(c, wls=BW_D[b])
            spb = dspec(c, H=HR[b], cbo=cbo_b)
            Rb = C.fit_step(spb, T, g, c["lo"], c["hi"], c["order"], fixed=fixed, free_cols=free)
            if Rs is None:
                hsb = C.h_lin(Rb, ysim); h2 = C.h_lin(Rb, simpart(Hs, c, ["C2"], c2=False)); h3 = C.h_lin(Rb, simpart(Hs, c, ["C3"], c2=False))
                h2c = C.h_lin(Rb, simpart(Hs, c, ["C2"]))
            else:
                hsb = hs; h2, h3, h2c = lin["C2"], lin["C3"], lin["C2c"]
            fb.append(Rb.h / hsb); fub.append(Rb.h / (h2 + h3)); ak.append(Rb.h / (h2c + 0.48 * h3)); hdb.append(Rb.h)
        out["boot"] = np.array(fb); out["boot_unc"] = np.array(fub); out["boot_k48"] = np.array(ak); out["boot_hd"] = np.array(hdb)
        out["statboot"] = C.sd(fb)
    # ---- simulation replicas: bump and peak templates, the estimator of the nominal fit
    if rep:
        HRs = H_SR(c)
        HRp = HRs if (c["peak_s"] is None and c["peak_r"] is None) else H_SR(c, c["peak_s"], c["peak_r"])
        fr, fur, akr, hsr, r2c, r3, rhd = [], [], [], [], [], [], []
        for r in range(NR_S):
            pk = simpart(HRp[r], c, ["C4", "C5"])
            yd = sp.sig - (c["aP"] * pk if c["peak"] == "fixed" else 0.)
            if c["peak"] == "free":
                Rr = C.fit_step(sp, T, g, c["lo"], c["hi"], c["order"], free_cols=[pk], Vfix=R.V)
                hdr = Rr.h; Lh = Rr.Lh; idx = Rr.idx
            else:
                hdr = float(R.Lh @ yd[R.idx]); Lh = R.Lh; idx = R.idx
            b2 = simpart(HRs[r], c, ["C2"], c2=False); b2c = simpart(HRs[r], c, ["C2"]); b3 = simpart(HRs[r], c, ["C3"])
            if c["simtarget"] == "bump":
                ys = b2c + b3
                hsr_ = float(Lh @ ys[idx]) if Rs is None else C.fit_step(None, Ts, g, c["lo"], c["hi"], c["order"], own_var=ys,
                                                                           Vfix=(R.V if c["simfit"] == "datacov" else None)).h
            else:
                hsr_ = float(Lh @ (b2c + b3)[idx]) + (hs - (lin["C2c"] + lin["C3"]))   # background part not replicated
            h2 = float(Lh @ b2[idx]); h3 = float(Lh @ b3[idx]); h2c = float(Lh @ b2c[idx])
            fr.append(hdr / hsr_); fur.append(hdr / (h2 + h3)); akr.append(hdr / (h2c + 0.48 * h3)); hsr.append(hsr_)
            r2c.append(h2c); r3.append(h3); rhd.append(hdr)
        out["rep"] = np.array(fr); out["rep_unc"] = np.array(fur); out["rep_k48"] = np.array(akr); out["rep_hs"] = np.array(hsr)
        out["rep_h2c"] = np.array(r2c); out["rep_h3"] = np.array(r3); out["rep_hd"] = np.array(rhd)
        out["mc"] = C.sd(fr); out["mc68"] = C.h68(fr)
    return out

# =====================================================================================================================
# 3. nominal
# =====================================================================================================================
P("\n3. NOMINAL (%.0f s)" % (time.time() - T0))
nom = evaluate(cfg(), keep=True)
R0 = nom["R"]
P("   h_data = %.2f +- %.2f rows / 2.5 MeV; chi2 %.1f / %d (without step %.1f); h_sim = %.2f (C2 + C3, C2 corrected)" %
  (nom["hd"], nom["ehd"], nom["chi2"], nom["ndf"], nom["chi20"], nom["hs"]))
P("   f_eta(u-edge) = %.4f +- %.4f (fit) +- %.4f (bootstrap; mean %.4f) +- %.4f (MC rms; robust %.4f; replica mean %.4f)" %
  (nom["f"], nom["statfit"], nom["statboot"], nom["boot"].mean(), nom["mc"], nom["mc68"], nom["rep"].mean()))
P("   uncorrected (no C2 factor): %.4f (boot %.4f, MC %.4f); kappa = 0.48 reading %.4f (unc. %.4f)" %
  (nom["f_unc"], C.sd(nom["boot_unc"]), C.sd(nom["rep_unc"]), nom["a_k48"], nom["a_k48_unc"]))
P("   linear decomposition of h_sim:", {k: round(v, 2) for k, v in nom["lin"].items()}, "C2 share %.3f" % nom["c2share"])
P("   step significance: delta chi2 = %.1f (%.1f sigma)" % (nom["chi20"] - nom["chi2"], math.sqrt(max(nom["chi20"] - nom["chi2"], 0))))
RES["nom"] = {k: (float(v) if isinstance(v, (float, np.floating)) else v) for k, v in nom.items() if k in
              ("hd", "ehd", "hs", "f", "statfit", "statboot", "mc", "mc68", "chi2", "chi20", "ndf", "f_unc", "a_k48", "a_k48_unc", "c2share", "cfall")}
RES["nom"]["lin"] = nom["lin"]
RES["nom"]["boot_mean"] = float(nom["boot"].mean()); RES["nom"]["rep_mean"] = float(nom["rep"].mean())
RES["nom"]["stat_unc"] = C.sd(nom["boot_unc"]); RES["nom"]["mc_unc"] = C.sd(nom["rep_unc"]); RES["nom"]["mc68_unc"] = C.h68(nom["rep_unc"])
RES["nom"]["stat_k48"] = C.sd(nom["boot_k48"]); RES["nom"]["mc_k48"] = C.sd(nom["rep_k48"])
RES["nom"]["hs_rep_rms"] = C.sd(nom["rep_hs"]); RES["nom"]["hs_rep_mean"] = float(nom["rep_hs"].mean())
RES["nom"]["fit_params"] = [float(x) for x in R0.p]
np.save(C.OUT + "nom_boot.npy", nom["boot"]); np.save(C.OUT + "nom_rep.npy", nom["rep"])

# figure data: spectra
g0 = grid_of(cfg())
Hs0 = nom["Hs"]
fd = dict(x=g0.x, N=nom["sp"].N, cb=nom["sp"].cb, sig=nom["sp"].sig, sigvar=nom["sp"].N + np.diag(nom["sp"].cbcov),
          peakT=nom["peakT"], T=nom["T"], fit_idx=R0.idx, fit_model=R0.model, fit_smooth=R0.smooth, fit_y=R0.y, fit_fx=R0.fx)
c0 = cfg()
for gname in C.GROUPS:
    fd["sim_" + gname] = simpart(Hs0, c0, [gname], c2=True)
fd["sim_C2raw"] = simpart(Hs0, c0, ["C2"], c2=False)
Hpre = H_D(cfg(sel="pre")); sppre = dspec(cfg(sel="pre"), H=Hpre)
fd["N_pre"] = sppre.N; fd["cb_pre"] = sppre.cb; fd["sig_pre"] = sppre.sig
fd["pp"] = nom["sp"].pp; fd["mm"] = nom["sp"].mm
# like-sign rows in u for the null test
# sim step view: the same fit on the simulated bump (own weights) for drawing
Rsim = C.fit_step(None, nom["T"], g0, -0.045, 0.045, 2, own_var=nom["ysim"], Vfix=R0.V)
fd["sim_fit_smooth"] = Rsim.smooth; fd["sim_fit_model"] = Rsim.model; fd["sim_y"] = Rsim.y
# resolution sample
fd["reso_d"] = dS; fd["reso_w"] = wdS
# truth u of C2/C3 and the reco u with and without q (for the kinematics figure)
for gg in (0, 1):
    m = (gS == gg) & np.isfinite(ut_s) & (w_nom_raw > 0) & (mpp_s0 >= 0.30) & (mpp_s0 < 0.40)
    fd["ut_%s" % C.GROUPS[gg]] = ut_s[m]; fd["utw_%s" % C.GROUPS[gg]] = w_nom_raw[m]
    fd["ur_%s" % C.GROUPS[gg]] = u_s0[m]; fd["umee_%s" % C.GROUPS[gg]] = S["mee"][m]
    u0q, _, _ = C.u_sim(S, useq=False)
    fd["ur0_%s" % C.GROUPS[gg]] = u0q[m]
    fd["m4p_%s" % C.GROUPS[gg]] = m4p_s0[m]; fd["mpp_%s" % C.GROUPS[gg]] = mpp_s0[m]
    fd["m4t_%s" % C.GROUPS[gg]] = S["m4t"][m]; fd["mppt_%s" % C.GROUPS[gg]] = S["mpipit"][m]; fd["meet_%s" % C.GROUPS[gg]] = S["meet"][m]
np.savez_compressed(C.OUT + "edgeu_figdata.npz", **fd)

# =====================================================================================================================
# 4. fit-range and order scan (choice of the nominal), free smearing
# =====================================================================================================================
P("\n4. RANGE / ORDER SCAN (no bootstrap; fit errors) (%.0f s)" % (time.time() - T0))
scan = []
for (lo, hi) in [(-0.030, 0.030), (-0.035, 0.035), (-0.040, 0.040), (-0.045, 0.045), (-0.050, 0.050), (-0.055, 0.055), (-0.060, 0.060),
                 (-0.070, 0.070), (-0.045, 0.030), (-0.045, 0.020), (-0.060, 0.045), (-0.030, 0.045), (-0.045, 0.060)]:
    for o in (1, 2, 3, 4):
        r = evaluate(cfg(lo=lo, hi=hi, order=o), boot=False, rep=False)
        p = 1 - chi2dist.cdf(r["chi2"], r["ndf"])
        scan.append(dict(lo=lo, hi=hi, order=o, f=r["f"], statfit=r["statfit"], chi2=r["chi2"], ndf=r["ndf"], p=p, hd=r["hd"], hs=r["hs"]))
        P("   %+.3f..%+.3f pol%d: f %.4f +- %.4f, chi2 %.1f/%d (p %.3f), h_d %.1f h_s %.1f" % (lo, hi, o, r["f"], r["statfit"], r["chi2"], r["ndf"], p, r["hd"], r["hs"]))
RES["scan"] = scan
# free smearing on the data and on the simulation (chi2 scan in sigma)
P("   free smearing (chi2 scan):")
sg = np.arange(0.0030, 0.0121, 0.0002)
Vh = C.hist_sim_var(S, u_s0, mpp_s0, sim_w(cfg()) * np.where((gS == 0), np.where(S["mee"] < 0.010, C.CF_E1, C.CF_E2), 1.0), grid_of(cfg()))
VSIM = simpart(Vh, cfg(), ["C2", "C3"], c2=False)   # event-group variance of the simulated bump (C2 factor in the weights)
fs = {}
for lab, own in [("data", None), ("sim", "sim")]:
    ch = []
    for s_ in sg:
        c_ = cfg(sigma=s_)
        g = grid_of(c_)
        T = C.step_gauss(g, s_)
        if own is None:
            sp = dspec(c_)
            R = C.fit_step(sp, T, g, -0.045, 0.045, 2, fixed=C.AP * nom["peakT"])
        else:
            R = C.fit_step(None, T, g, -0.045, 0.045, 2, own_var=nom["ysim"], Vfix=np.diag(np.maximum(VSIM[R0.idx], 1e-6)))
        ch.append(R.chi2)
    ch = np.array(ch); i = int(np.argmin(ch))
    a = np.polyfit(sg[max(0, i - 6):i + 7], ch[max(0, i - 6):i + 7], 2)
    smin = -a[1] / (2 * a[0]); es = 1 / math.sqrt(a[0]) if a[0] > 0 else float("nan")
    fs[lab] = dict(sigma=float(smin), err=float(es), chi2min=float(ch[i]), grid=sg.tolist(), chi2=ch.tolist())
    P("     %s: sigma = %.2f +- %.2f MeV (chi2 min %.1f)" % (lab, 1e3 * smin, 1e3 * es, ch[i]))
RES["free_sigma"] = fs

# =====================================================================================================================
# 5. variants (systematic families) with bootstrap and replica noise
# =====================================================================================================================
P("\n5. VARIANTS (%.0f s)" % (time.time() - T0))
# (tag, family, pair id, counted, r1 type, config overrides, text)
V = [
    # U1 fit model
    ("pol1", "U1", -1, 1, 1, dict(order=1), "polynomial order 1"),
    ("pol3", "U1", -1, 1, 1, dict(order=3), "polynomial order 3"),
    ("r35", "U1", 1, 1, 2, dict(lo=-0.035, hi=0.035), "range -35..+35 MeV"),
    ("r55", "U1", 1, 1, 2, dict(lo=-0.055, hi=0.055), "range -55..+55 MeV"),
    ("r40", "U1", 2, 1, 2, dict(lo=-0.040, hi=0.040), "range -40..+40 MeV"),
    ("r50", "U1", 2, 1, 2, dict(lo=-0.050, hi=0.050), "range -50..+50 MeV"),
    ("loM60", "U1", 3, 1, 2, dict(lo=-0.060), "lower end -60 MeV"),
    ("loM30", "U1", 3, 1, 2, dict(lo=-0.030), "lower end -30 MeV"),
    ("hiP60", "U1", 4, 1, 2, dict(hi=0.060), "upper end +60 MeV"),
    ("hiP30", "U1", 4, 1, 2, dict(hi=0.030), "upper end +30 MeV"),
    ("simown", "U1", -1, 1, 0, dict(simfit="own"), "simulation fitted with its own Poisson-like weights"),
    ("tcentre", "U1", -1, 1, 0, dict(tmpl="centre"), "step template at the bin centre (scan::StepTemplate) instead of the bin average"),
    ("T12", "U1", -1, 0, 0, dict(simtarget="T12"), "check: h_sim from all simulated classes (bedge/WP6 convention)"),
    # U2 step smearing
    ("sigData", "U2", -1, 1, 0, dict(sigma=SIG_DATA, sigma_sim=SIG_U), "own smearing per sample: data %.2f MeV (x 0.893), sim %.2f MeV" % (1e3 * SIG_DATA, 1e3 * SIG_U)),
    ("sigUp", "U2", 5, 1, 0, dict(sigma=SIG_U + 0.0010), "smearing +1 MeV (both)"),
    ("sigDn", "U2", 5, 1, 0, dict(sigma=SIG_U - 0.0010), "smearing -1 MeV (both)"),
    ("sig54", "U2", -1, 1, 0, dict(sigma=0.0054), "smearing 5.4 MeV (the frozen bedge/WP6 value)"),
    ("emp", "U2", -1, 1, 0, dict(tmpl="emp"), "empirical step: the simulated u resolution with its tails"),
    # U3 M(pipi) range and CB slices
    ("pp28", "U3", -1, 1, 0, dict(pplo=0.28), "M(pipi) from 0.28"),
    ("pp29", "U3", 6, 1, 0, dict(pplo=0.29), "M(pipi) from 0.29"),
    ("pp31", "U3", 6, 1, 0, dict(pplo=0.31), "M(pipi) from 0.31"),
    ("pp39", "U3", -1, 1, 0, dict(pphi=0.39), "M(pipi) up to 0.39"),
    ("pp41", "U3", -1, 1, 0, dict(pphi=0.41, sledges=(0.30, 0.34, 0.37, 0.39, 0.40, 0.41)), "M(pipi) up to 0.41 (boxes close inside the range)"),
    # U4 mass scales and edge function
    ("dsUp", "U4", 7, 1, 0, dict(dscale=C.EDGE_SCALE * 1.0009), "data edge scale x 1.0009"),
    ("dsDn", "U4", 7, 1, 0, dict(dscale=C.EDGE_SCALE / 1.0009), "data edge scale / 1.0009"),
    ("dpp52", "U4", -1, 1, 0, dict(dspp=0.9952), "data M(pipi) scale 0.9952 (eta-decay pions as in the simulation) instead of 0.9970"),
    ("dppUp", "U4", 8, 1, 0, dict(dspp=0.9974), "data M(pipi) scale 0.9974"),
    ("dppDn", "U4", 8, 1, 0, dict(dspp=0.9966), "data M(pipi) scale 0.9966"),
    ("q0", "U4", -1, 1, 0, dict(useq=False), "edge without the M(ee) term, E(m, 0) (data and simulation)"),
    ("q0all", "U4", -1, 0, 0, dict(useq=False, qcut=False), "check: E(m, 0) and all M(ee) rows (the bedge row set)"),
    ("sUp", "U4", 9, 1, 0, dict(s=C.S_NOM + C.S_ERR), "simulation s = 0.99705"),
    ("sDn", "U4", 9, 1, 0, dict(s=C.S_NOM - C.S_ERR), "simulation s = 0.99525"),
    ("sK0", "U4", -1, 1, 0, dict(s=0.99521), "simulation s = 0.99521 (K0S)"),
    ("r089", "U4", 10, 1, 0, dict(r=0.89), "simulation r = 0.89 (data-like M4 resolution)"),
    ("r110", "U4", 10, 1, 0, dict(r=1.10), "simulation r = 1.10"),
    ("sppUp", "U4", 11, 1, 0, dict(spp=C.SPP_NOM + 0.0007), "simulation M(pipi) x 0.9959"),
    ("sppDn", "U4", 11, 1, 0, dict(spp=C.SPP_NOM - 0.0007), "simulation M(pipi) x 0.9945"),
    ("rpp116", "U4", -1, 1, 0, dict(rpp=1.16), "simulation M(pipi) resolution x 1.16 (data/sim K0S width)"),
    # U5 CB
    ("cb095", "U5", 12, 1, 0, dict(cbk=0.95), "CB x 0.95"),
    ("cb105", "U5", 12, 1, 0, dict(cbk=1.05), "CB x 1.05"),
    ("cb1032", "U5", -1, 1, 0, dict(cbk=1.032), "CB x 1.032 (k from data mixing)"),
    ("cb5", "U5", 13, 1, 0, dict(grid=(-0.200, 0.0025, 140, 2)), "CB cells 5 MeV"),
    ("cb20", "U5", 13, 1, 0, dict(grid=(-0.200, 0.0025, 140, 8)), "CB cells 20 MeV"),
    ("cbsum", "U5", -1, 1, 0, dict(cbmode=1), "CB = N++ + N--"),
    ("cbglob", "U5", -1, 1, 0, dict(cbslice=False), "one CB for 0.30-0.40 (no slices)"),
    ("cb5mev", "U5", -1, 1, 0, dict(sledges=tuple(np.round(np.arange(0.30, 0.4001, 0.005), 3))), "CB per 5 MeV M(pipi) unit"),
    ("cbmap", "U5", -1, 1, 0, dict(cbmap=True), "CB of the frozen M4 rule, distributed over u by the like-sign shape per 50 MeV band"),
    ("cb25", "U5", -1, 1, 0, dict(grid=(-0.200, 0.0025, 140, 1)), "CB cells 2.5 MeV (no smoothing)"),
    ("cbPP", "U5", 16, 1, 0, dict(cbmode=2), "CB shape from e+e+ only (no edge step), normalised to the nominal CB"),
    ("cbMM", "U5", 16, 1, 0, dict(cbmode=3), "CB shape from e-e- only (with its edge step), normalised to the nominal CB"),
    # U6 binning
    ("bw20", "U6", -1, 1, 0, dict(grid=(-0.200, 0.0020, 175, 5)), "2.0 MeV bins"),
    ("bw50", "U6", -1, 1, 0, dict(grid=(-0.200, 0.0050, 70, 2)), "5.0 MeV bins"),
    ("off", "U6", -1, 1, 0, dict(grid=(-0.19875, 0.0025, 140, 4)), "grid shifted by 1.25 MeV"),
    # U7 physics weights, ownership
    ("wtsOff", "U7", -1, 0, 0, dict(code=0), "physics weights off (phase-space eta -> 3pi Dalitz plot)"),
    ("bump1aY", "U7", -1, 1, 0, dict(code=116), "eta -> 3pi Dalitz plot 1 + aY"),
    ("bumpPS", "U7", -1, 0, 0, dict(code=100), "eta -> 3pi phase space"),
    ("V2", "U7", -1, 1, 0, dict(code=107), "direct decay V2 (peak template)"),
    ("convOff", "U7", -1, 1, 0, dict(code=76), "conversion F_V off (peak template)"),
    ("modeB", "U7", -1, 1, 0, dict(mode="B"), "ownership mode B"),
    # U8 peak term
    ("aPup", "U8", 14, 1, 0, dict(aP=C.AP + math.hypot(math.hypot(C.AP_STAT, C.AP_SYST), C.AP_MC)), "a_P + 1 sigma (0.080)"),
    ("aPdn", "U8", 14, 1, 0, dict(aP=C.AP - math.hypot(math.hypot(C.AP_STAT, C.AP_SYST), C.AP_MC)), "a_P - 1 sigma"),
    ("aPfree", "U8", -1, 1, 0, dict(peak="free"), "peak amplitude free (simulated u shape)"),
    ("pkr089", "U8", -1, 1, 0, dict(peak_r=0.89), "peak template with r = 0.89 (data peak width)"),
    ("pksUp", "U8", -1, 1, 0, dict(peak_s=C.S_NOM + C.S_ERR), "peak template with s + 0.0009"),
    ("nopeak", "U8", -1, 0, 0, dict(peak="none"), "check: no peak term"),
    # U9 selection
    ("pre", "U9", -1, 1, 0, dict(sel="pre"), "selection pre (data and simulation)"),
    # checks
    ("best", "chk", -1, 0, 0, dict(best=True), "check: isBest == 1 in data and simulation"),
    ("bestD", "chk", -1, 0, 0, dict(best=True, sbest=False), "check: isBest == 1 in the data only (WP7 convention)"),
    ("c2off", "chk", -1, 0, 0, dict(c2corr=False), "check: no C2 angle correction"),
    ("cfUp", "U10", 15, 1, 0, dict(cfE1=C.CF_E1 + C.CF_ERR, cfE2=C.CF_E2 + C.CF_ERR), "C2 factor + 0.024"),
    ("cfDn", "U10", 15, 1, 0, dict(cfE1=C.CF_E1 - C.CF_ERR, cfE2=C.CF_E2 - C.CF_ERR), "C2 factor - 0.024"),
]
VR = {}
for (tag, fam, pair, counted, r1, ov, text) in V:
    t1 = time.time()
    r = evaluate(cfg(**ov))
    d = dict(fam=fam, pair=pair, counted=counted, r1=r1, text=text, f=r["f"], statfit=r["statfit"], statboot=r["statboot"], mc=r["mc"],
             mc68=r["mc68"], chi2=r["chi2"], ndf=r["ndf"], hd=r["hd"], hs=r["hs"], f_unc=r["f_unc"], a_k48=r["a_k48"],
             delta=r["f"] - nom["f"], nb=C.sd(r["boot"] - nom["boot"]), nr=C.sd(r["rep"] - nom["rep"]), nr68=C.h68(r["rep"] - nom["rep"]),
             d_unc=(r["f_unc"] - nom["f_unc"]) if r["f_unc"] is not None else None, aPfit=r["aPfit"],
             nb_unc=C.sd(r["boot_unc"] - nom["boot_unc"]), nr68_unc=C.h68(r["rep_unc"] - nom["rep_unc"]))
    VR[tag] = d
    P("   %-8s %-3s f %.4f (unc %s) delta %+.4f  noise boot %.4f rep %.4f (68: %.4f)  chi2 %.1f/%d  [%s] %.0f s" %
      (tag, fam, r["f"], ("%.4f" % r["f_unc"]) if r["f_unc"] is not None else "-", d["delta"], d["nb"], d["nr"], d["nr68"], r["chi2"], r["ndf"], text, time.time() - t1))
RES["var"] = VR

# =====================================================================================================================
# 6. checks: per slice, per M(ee) class, periods, folds, null tests, control region
# =====================================================================================================================
P("\n6. CHECKS (%.0f s)" % (time.time() - T0))
CHK = {}
# per M(pipi) slice (own polynomial each)
sl_res = []
for a, b in [(0.30, 0.34), (0.34, 0.37), (0.37, 0.39), (0.39, 0.40)]:
    r = evaluate(cfg(pplo=a, pphi=b), rep=True)
    sl_res.append(dict(lo=a, hi=b, hd=r["hd"], ehd=r["ehd"], hs=r["hs"], f=r["f"], statfit=r["statfit"], statboot=r["statboot"], mc=r["mc"],
                       chi2=r["chi2"], ndf=r["ndf"], boot_hd=r["boot_hd"].tolist(), rep_hd=r["rep_hd"].tolist(), rep_hs=r["rep_hs"].tolist(),
                       lin=r["lin"]))
    P("   slice %.2f-%.2f: h_d %.1f +- %.1f, h_s %.1f, f %.4f +- %.4f (boot) +- %.4f (MC), chi2 %.1f/%d" % (a, b, r["hd"], r["ehd"], r["hs"], r["f"], r["statboot"], r["mc"], r["chi2"], r["ndf"]))
fsl = np.array([x["f"] for x in sl_res]); esl = np.array([math.hypot(x["statboot"], x["mc"]) for x in sl_res]); estat = np.array([x["statboot"] for x in sl_res])
wm = np.sum(fsl / esl ** 2) / np.sum(1 / esl ** 2)
CHK["slices"] = sl_res
CHK["slices_chi2"] = float(np.sum((fsl - wm) ** 2 / esl ** 2)); CHK["slices_chi2_statonly"] = float(np.sum((fsl - np.sum(fsl / estat ** 2) / np.sum(1 / estat ** 2)) ** 2 / estat ** 2))
CHK["slices_sum_hd"] = float(sum(x["hd"] for x in sl_res)); CHK["slices_sum_hs"] = float(sum(x["hs"] for x in sl_res))
CHK["slices_ratio_of_sums"] = CHK["slices_sum_hd"] / CHK["slices_sum_hs"]
# the four-slice estimator in u (one pol2 per slice, f = sum h_data / sum h_sim): bootstrap (h_sim fixed) and replicas
fb4 = np.sum([x["boot_hd"] for x in sl_res], axis=0) / CHK["slices_sum_hs"]
fr4 = np.sum([x["rep_hd"] for x in sl_res], axis=0) / np.sum([x["rep_hs"] for x in sl_res], axis=0)
CHK["sl4"] = dict(f=CHK["slices_ratio_of_sums"], statboot=C.sd(fb4), mc=C.sd(fr4), mc68=C.h68(fr4),
                  f_unc=CHK["slices_sum_hd"] / sum(x["lin"]["C2"] + x["lin"]["C3"] for x in sl_res))
P("   four-slice estimator in u: f %.4f +- %.4f (boot) +- %.4f (MC rms; robust %.4f); uncorrected %.4f" %
  (CHK["sl4"]["f"], CHK["sl4"]["statboot"], CHK["sl4"]["mc"], CHK["sl4"]["mc68"], CHK["sl4"]["f_unc"]))
P("   slices: chi2 (stat+MC) %.2f/3, stat only %.2f/3; sum h_d %.1f, sum h_s %.1f, ratio %.4f" % (CHK["slices_chi2"], CHK["slices_chi2_statonly"], CHK["slices_sum_hd"], CHK["slices_sum_hs"], CHK["slices_ratio_of_sums"]))
# per M(ee) class
cl_res = {}
for cl in ("E1", "E2", "E20"):
    r = evaluate(cfg(cls=cl), rep=True)
    cl_res[cl] = dict(hd=r["hd"], ehd=r["ehd"], hs=r["hs"], f=r["f"], statfit=r["statfit"], statboot=r["statboot"], mc=r["mc"], mc68=r["mc68"],
                      chi2=r["chi2"], ndf=r["ndf"], f_unc=r["f_unc"], lin=r["lin"], boot=r["boot"].tolist(), rep=r["rep"].tolist(),
                      boot_hd=r["boot_hd"].tolist(), rep_hs=r["rep_hs"].tolist(), rep_h2c=r["rep_h2c"].tolist(), rep_h3=r["rep_h3"].tolist(),
                      rep_hd=r["rep_hd"].tolist())
    P("   class %s: h_d %.1f +- %.1f, h_s %.1f, f %.4f +- %.4f +- %.4f (unc %.4f), chi2 %.1f/%d; lin C2c %.1f C3 %.1f" %
      (cl, r["hd"], r["ehd"], r["hs"], r["f"], r["statboot"], r["mc"], r["f_unc"], r["chi2"], r["ndf"], r["lin"]["C2c"], r["lin"]["C3"]))
    # own-rule CB for the class
    r2 = evaluate(cfg(cls=cl, clsrule=1), boot=False, rep=False)
    cl_res[cl]["f_ownCB"] = r2["f"]
CHK["classes"] = cl_res
# kappa free from the E1/E2 split: h_d,i = a (c2 C2_i + kappa C3_i) (C2_i already corrected)
def kappa_solve(e1, e2, h1, h2):
    A1, B1 = e1["C2c"], e1["C3"]; A2, B2 = e2["C2c"], e2["C3"]
    rho = h1 / h2
    kap = (rho * A2 - A1) / (B1 - rho * B2)
    a = h2 / (A2 + kap * B2)
    return a, kap
a_, k_ = kappa_solve(cl_res["E1"]["lin"], cl_res["E2"]["lin"], cl_res["E1"]["hd"], cl_res["E2"]["hd"])
kb = [kappa_solve(cl_res["E1"]["lin"], cl_res["E2"]["lin"], x, y) for x, y in zip(cl_res["E1"]["boot_hd"], cl_res["E2"]["boot_hd"])]
kr_ = [kappa_solve(dict(C2c=cl_res["E1"]["rep_h2c"][i], C3=cl_res["E1"]["rep_h3"][i]), dict(C2c=cl_res["E2"]["rep_h2c"][i], C3=cl_res["E2"]["rep_h3"][i]),
                    cl_res["E1"]["rep_hd"][i], cl_res["E2"]["rep_hd"][i]) for i in range(NR_S)]
CHK["kappa_free"] = dict(a=a_, kappa=k_, a_stat=C.sd([x[0] for x in kb]), kappa_stat=C.sd([x[1] for x in kb]),
                         a_mc=C.sd([x[0] for x in kr_]), kappa_mc=C.sd([x[1] for x in kr_]), a_mc68=C.h68([x[0] for x in kr_]), kappa_mc68=C.h68([x[1] for x in kr_]))
P("   E1/E2 split: a_D %.3f +- %.3f (stat) +- %.3f (MC rms; robust %.3f), kappa %.3f +- %.3f (stat) +- %.3f (MC rms; robust %.3f)" %
  (a_, CHK["kappa_free"]["a_stat"], CHK["kappa_free"]["a_mc"], CHK["kappa_free"]["a_mc68"], k_, CHK["kappa_free"]["kappa_stat"], CHK["kappa_free"]["kappa_mc"], CHK["kappa_free"]["kappa_mc68"]))
# periods (simulation x L_p / L), folds (simulation x row fraction)
pr = []
for p in range(4):
    r = evaluate(cfg(period=p, simscale=C.LP[p] / C.L_DATA), rep=False)
    pr.append(dict(f=r["f"], statfit=r["statfit"], statboot=r["statboot"], fc=r["f"] * C.CSEC_P[p], Lp=C.LP[p]))
    P("   period %s: f %.4f +- %.4f (boot) x C_sec %.4f = %.4f" % (["I", "IIA", "IIB", "III"][p], r["f"], r["statboot"], C.CSEC_P[p], r["f"] * C.CSEC_P[p]))
fp = np.array([x["fc"] for x in pr]); ep = np.array([x["statboot"] * C.CSEC_P[i] for i, x in enumerate(pr)])
mp_ = np.sum(fp / ep ** 2) / np.sum(1 / ep ** 2); chp = float(np.sum((fp - mp_) ** 2 / ep ** 2))
CHK["periods"] = pr; CHK["periods_chi2"] = chp; CHK["periods_p"] = float(1 - chi2dist.cdf(chp, 3))
fpr = np.array([x["f"] for x in pr]); epr = np.array([x["statboot"] for x in pr])
mpr = np.sum(fpr / epr ** 2) / np.sum(1 / epr ** 2); CHK["periods_chi2_raw"] = float(np.sum((fpr - mpr) ** 2 / epr ** 2))
CHK["periods_p_raw"] = float(1 - chi2dist.cdf(CHK["periods_chi2_raw"], 3))
P("   periods: chi2 %.2f/3 (p %.3f) with C_sec; raw %.2f/3 (p %.3f)" % (chp, CHK["periods_p"], CHK["periods_chi2_raw"], CHK["periods_p_raw"]))
u_nom = U_D(cfg()); pbn = C.ppbin_of(D["mpipi"])
inr = (D["bits"] & 1 != 0) & (D["k"] == 0) & (pbn >= C.PB(0.30)) & (pbn < C.PB(0.40))
fo = []
for k_f in (0, 1):
    frac = float(np.sum(inr & (D["fold"] == k_f)) / np.sum(inr))
    r = evaluate(cfg(fold=k_f, simscale=frac), rep=False)
    fo.append(dict(f=r["f"], statboot=r["statboot"], statfit=r["statfit"], frac=frac))
    P("   fold %d (row fraction %.4f): f %.4f +- %.4f (boot)" % (k_f, frac, r["f"], r["statboot"]))
z = (fo[0]["f"] - fo[1]["f"]) / math.hypot(fo[0]["statboot"], fo[1]["statboot"])
CHK["folds"] = fo; CHK["folds_z"] = z; CHK["folds_p"] = float(1 - chi2dist.cdf(z * z, 1))
P("   folds: z = %.2f (p %.3f)" % (z, CHK["folds_p"]))
# null tests: step of the like-sign rows, of the CB, of the simulated non-eta and other-eta parts, control region
c0 = cfg(); g = grid_of(c0)
sp0 = nom["sp"]
ls_y = sp0.pp + sp0.mm
Rls = C.fit_step(C.Spec(ls_y, np.zeros(g.nb), np.zeros((g.nb, g.nb))), nom["T"], g, -0.045, 0.045, 2)
Rcb = C.fit_step(None, nom["T"], g, -0.045, 0.045, 2, own_var=sp0.cb, Vfix=R0.V)
CHK["null_ls"] = dict(h=Rls.h, eh=Rls.eh, chi2=Rls.chi2, ndf=Rls.ndf, rel=Rls.h / nom["hd"])
Hn = H_D(c0)
for kk, kl in ((1, "pp"), (2, "mm")):
    for cl in ("all", "E1", "E2"):
        cs = C.CLS_SETS[cl]
        y = sum(Hn[kk, a:b][:, cs].sum((0, 1)) for a, b in slices(c0))
        Rx = C.fit_step(C.Spec(y, np.zeros(g.nb), np.zeros((g.nb, g.nb))), nom["T"], g, -0.045, 0.045, 2)
        CHK["null_ls_%s_%s" % (kl, cl)] = dict(h=Rx.h, eh=Rx.eh, chi2=Rx.chi2, ndf=Rx.ndf, sum=float(y[Rx.idx].sum()))
        P("   null: like-sign %s %s: h = %.1f +- %.1f (chi2 %.1f/%d), rows in range %.0f" % (kl, cl, Rx.h, Rx.eh, Rx.chi2, Rx.ndf, y[Rx.idx].sum()))
for p_ in range(4):
    Hp_ = H_D(cfg(period=p_))
    y = sum(Hp_[2, a:b].sum((0, 1)) for a, b in slices(c0))
    Rx = C.fit_step(C.Spec(y, np.zeros(g.nb), np.zeros((g.nb, g.nb))), nom["T"], g, -0.045, 0.045, 2)
    CHK["null_ls_mm_per%d" % p_] = dict(h=Rx.h, eh=Rx.eh)
    P("   null: like-sign mm, period %d: h = %.1f +- %.1f" % (p_, Rx.h, Rx.eh))
CHK["null_cb"] = dict(h=Rcb.h, h_lin=C.h_lin(R0, sp0.cb), rel=C.h_lin(R0, sp0.cb) / nom["hd"])
P("   null: like-sign rows N++ + N--: h = %.1f +- %.1f (chi2 %.1f/%d), i.e. %.3f of h_data; CB itself (data estimator): h = %.1f (%.3f of h_data)" %
  (Rls.h, Rls.eh, Rls.chi2, Rls.ndf, CHK["null_ls"]["rel"], CHK["null_cb"]["h_lin"], CHK["null_cb"]["rel"]))
# simulated background parts: raw projection with the data estimator and its replica noise
ROWS_BG = np.where(gS >= 4)[0]
BW_S_BG = C.sim_replica_weights(S, ROWS_BG, NR_S)
u, mpp, _ = U_S(us_key(c0)); w = sim_w(c0)
Ssub = {kk: S[kk][ROWS_BG] for kk in ("g", "mee")}
HbgR = C.hist_sim(Ssub, u[ROWS_BG], mpp[ROWS_BG], w[ROWS_BG], g, reps=BW_S_BG)
del BW_S_BG
nb_ = {}
for lab, grp in [("Geta", ["Geta"]), ("GD", ["GDom", "GDnom"]), ("GC", ["GCom", "GCnom"]), ("nonEta", ["GDom", "GDnom", "GCom", "GCnom"])]:
    y = simpart(nom["Hs"], c0, grp); h = C.h_lin(R0, y)
    hr = [C.h_lin(R0, simpart(HbgR[r], c0, grp)) for r in range(NR_S)]
    Rown = C.fit_step(None, nom["T"], g, -0.045, 0.045, 2, own_var=y)
    nb_[lab] = dict(h=h, rep_rms=C.sd(hr), rep_mean=float(np.mean(hr)), rel_hs=h / nom["hs"], own=Rown.h, own_e=Rown.eh, sum=float(y[R0.idx].sum()))
    P("   null: simulated %s: h = %.1f +- %.1f (replica rms) rows/2.5 MeV = %.3f of h_sim; own-weight fit %.1f +- %.1f" % (lab, h, C.sd(hr), h / nom["hs"], Rown.h, Rown.eh))
CHK["null_sim"] = nb_
# control region without a bump: M(pipi) 0.42-0.46, pseudo edge E(M(pipi) - 0.10, q)
cc = cfg(pplo=0.42, pphi=0.46, sledges=(0.42, 0.46), mshift=-0.10)
rc = evaluate(cc, rep=False)
CHK["control"] = dict(hd=rc["hd"], ehd=rc["ehd"], hs=rc["hs"], f=rc["f"], chi2=rc["chi2"], ndf=rc["ndf"], statboot=rc["statboot"])
P("   control M(pipi) 0.42-0.46 with the pseudo edge E(M(pipi) - 0.10, M(ee)): h_data %.1f +- %.1f (chi2 %.1f/%d); simulated bump there h_sim %.1f" %
  (rc["hd"], rc["ehd"], rc["chi2"], rc["ndf"], rc["hs"]))
cc2 = cfg(pplo=0.30, pphi=0.40, mshift=-0.10, sledges=(0.30, 0.34, 0.37, 0.39, 0.40))
# a shifted pseudo edge in the bump region: the step must vanish (u = M4 - E(M(pipi) - 0.10)) is not a clean null (the
# bump edge moves to u = +E(m) - E(m - 0.1)); not used
RES["chk"] = CHK

# =====================================================================================================================
# 7. smoothed background projections (E7-like term) and step projection of the simulated background deformations
# =====================================================================================================================
P("\n7. BACKGROUND STEP PROJECTION (%.0f s)" % (time.time() - T0))
def smooth_exp(y, x, lo, hi, order):
    """Poisson-likelihood fit of exp(pol n) to a (weighted) histogram on [lo, hi); orders raised one by one, each started
    from the previous solution; returns the smooth function on the grid (0 outside [lo, hi))"""
    from scipy.optimize import minimize
    m = (x >= lo) & (x < hi)
    xx = (x[m] - 0.5 * (lo + hi)) / (0.5 * (hi - lo)); yy = np.maximum(y[m], 0)
    def nll(p):
        mu = np.exp(np.clip(np.polyval(p[::-1], xx), -50, 50))
        return float(np.sum(mu - yy * np.log(mu + 1e-300)))
    p = np.array([math.log(max(yy.mean(), 1e-6))])
    for o in range(1, order + 1):
        p = minimize(nll, np.r_[p, 0.0], method="BFGS").x
    out = np.zeros_like(x)
    out[m] = np.exp(np.polyval(p[::-1], xx))
    return out
bgp = {}
ybg = simpart(nom["Hs"], c0, ["GDom", "GDnom", "GCom", "GCnom"]); yeta = simpart(nom["Hs"], c0, ["Geta"])
for lo, hi in [(-0.150, 0.150), (-0.100, 0.100)]:
    for o in (2, 3, 4, 5, 6):
        sm = smooth_exp(ybg, g.x, lo, hi, o); h = C.h_lin(R0, sm)
        sme = smooth_exp(yeta, g.x, lo, hi, o); he = C.h_lin(R0, sme)
        bgp["%g_%g_pol%d" % (lo, hi, o)] = dict(h_bg=h, h_eta=he)
        P("   smoothed exp(pol%d) on %.3f..%.3f: non-eta step projection %.2f rows (%.4f of h_sim), G_eta %.2f" % (o, lo, hi, h, h / nom["hs"], he))
RES["bgproj"] = bgp
# deformations of the non-eta background in M4: (1 + b x + q x^2), x = (M4' - 0.55)/0.1; Gaussian structures in M4
x4 = (m4p_s0 - 0.55) / 0.1
def bg_weighted(fac):
    ww = sim_w(c0) * fac
    Hh = C.hist_sim(S, u_s0, mpp_s0, ww, g, groups=[4, 5, 6, 7, 8])
    return Hh
defs_ = {"x": x4, "x2": x4 * x4,
         "g490_15": np.exp(-0.5 * ((m4p_s0 - 0.490) / 0.015) ** 2), "g505_8": np.exp(-0.5 * ((m4p_s0 - 0.505) / 0.008) ** 2),
         "ypp": (mpp_s0 - 0.35) / 0.05}
dproj = {}
for k_, fac in defs_.items():
    Hh = bg_weighted(fac)
    ynon = simpart(Hh, c0, ["GDom", "GDnom", "GCom", "GCnom"])
    dproj[k_] = dict(raw=C.h_lin(R0, ynon))
    for o in (3, 4):
        dproj[k_]["sm%d" % o] = C.h_lin(R0, smooth_exp(ynon, g.x, -0.150, 0.150, o))
    P("   deformation %s: raw step projection %.2f, smoothed exp(pol3) %.2f, exp(pol4) %.2f rows (h_sim %.1f)" % (k_, dproj[k_]["raw"], dproj[k_]["sm3"], dproj[k_]["sm4"], nom["hs"]))
RES["defproj"] = dproj

# =====================================================================================================================
# 8. windows (visible yields) with the simulation replicas
# =====================================================================================================================
P("\n8. WINDOWS AND YIELDS (%.0f s)" % (time.time() - T0))
c0 = cfg()
ws = {}
wr = sim_w(c0)
c2f = np.where(S["mee"] < 0.010, C.CF_E1, C.CF_E2)
inB = (m4p_s0 >= 0.4550) & (m4p_s0 < 0.5225)
for lab, sel_ in [("B", inB), ("all", np.ones(S["n"], bool)), ("B_E2", inB & (S["mee"] >= 0.010)), ("B_E1", inB & (S["mee"] < 0.010))]:
    m2 = (gS == 0) & sel_; m3 = (gS == 1) & sel_
    W2 = float(wr[m2].sum()); W2c = float((wr * c2f)[m2].sum()); W3 = float(wr[m3].sum())
    # replicas
    sub = np.isin(ROWS_EP, np.where(m2 | m3)[0])
    rows = ROWS_EP[sub]; kk = BW_S_EP[:, sub]
    W2r = (kk * (wr * c2f)[rows] * (gS[rows] == 0)).sum(1); W3r = (kk * wr[rows] * (gS[rows] == 1)).sum(1)
    W2ur = (kk * wr[rows] * (gS[rows] == 0)).sum(1)
    ws[lab] = dict(C2=W2, C2c=W2c, C3=W3, rep_c=(W2r + W3r).tolist(), rep_u=(W2ur + W3r).tolist())
    P("   window %s: SMASH C2 %.1f (corrected %.1f), C3 %.1f raw rows; (C2 + C3) x 64/63 = %.1f" % (lab, W2, W2c, W3, (W2 + W3) * C.PT3))
RES["windows"] = ws

json.dump(RES, open(C.OUT + "edgeu_ana.json", "w"), indent=1, default=float)
P("\ndone in %.0f s" % (time.time() - T0))
