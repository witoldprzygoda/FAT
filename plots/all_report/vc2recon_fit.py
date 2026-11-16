#!/usr/bin/env python3
# ============================================================================
# vc2recon_fit.py -- independent verification of package c2recon, part 3: the effective C2 factor of the u-edge
#                    estimator (edgeu) and the row list / bootstrap of the per-slice edge estimator (bedge)
# ============================================================================
# Written 2026-11-15 (prefix vc2recon, verifier). Uses the FROZEN fit code of edgeu_core.py (data spectrum, CB rule,
# fit_step, h_lin) and the edgeu column files; the reweighted C2 templates are built here from the verifier's own
# per-row weights (vc2_w.npz of vc2recon_w.py). No c2recon file is read.
#   mode "edgeu": nominal reproduction, own C2 template (checked against hist_sim), the factor = h(C2 QED) / h(C2),
#                 per class, separate class fits, fit variants, group-Poisson replicas, form factor, f and sigma;
#                 writes vc2_bedge_rows.txt for vc2recon_bedge.C;
#   mode "bedge": reads vc2_bedge.txt (vc2recon_bedge.C) and does the row-level replicas of the per-slice factor.
# Output: /hdd3/.../vc2recon/vc2_fit_<mode>.json / .log
# Run: python3 -I plots/all_report/vc2recon_fit.py edgeu ; (root vc2recon_bedge.C) ; python3 -I ... bedge
# ============================================================================
import json
import math
import sys
import time

import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import edgeu_core as C

O = "/hdd3/przygoda/HADES/FAT_out/all_report/vc2recon/"
MODE = sys.argv[1] if len(sys.argv) > 1 else "edgeu"
T0 = time.time()
LOG = open(O + "vc2_fit_%s.log" % MODE, "w")
def P(*a):
    s = " ".join(str(x) for x in a)
    print(s); LOG.write(s + "\n"); LOG.flush()
RES = {}
V = np.load(O + "vc2_w.npz")
seq = V["seq"]
MEL, MPI0 = 0.51099895, 134.9768
NREP = 2000
rng = np.random.default_rng(0x7C2)

def variants():
    """weight multipliers per C2 row (order of seq)"""
    out = {"qed": V["w_circ"], "best": V["w_best"]}
    for k in ("A", "B", "flat", "C_coarse", "C_raw", "pdgmass", "kloeX", "noeps"):
        if "w_" + k in V.files:
            out[k] = V["w_" + k]
    # pi0 transition form factor on top of the QED weight: |F|^2 = (1 + a M^2/m^2)^2 / <..>_KW, a = 0.0332 (PDG)
    x = np.linspace(2 * MEL * 1.000001, MPI0 * 0.999999, 200001)
    kw = (2 / x) * (1 - x * x / MPI0 ** 2) ** 3 * (1 + 2 * MEL ** 2 / x ** 2) * np.sqrt(1 - 4 * MEL ** 2 / x ** 2)
    a = 0.0332
    F2 = lambda m: (1 + a * m * m / MPI0 ** 2) ** 2
    norm = np.trapz(kw * F2(x), x) / np.trapz(kw, x)
    mt = 1e3 * V["meet"]
    out["ff"] = F2(mt) / norm
    out["qed_ff"] = V["w_circ"] * F2(mt) / norm
    RES["ff_norm"] = float(norm)
    return out

def poisson_groups(tree, grp, n):
    key = tree.astype(np.int64) * (1 << 32) + grp.astype(np.int64)
    uq, inv = np.unique(key, return_inverse=True)
    G = rng.poisson(1.0, size=(n, len(uq))).astype(np.float64)
    return G[:, inv]

# ======================================================================================================================
if MODE == "edgeu":
    D = C.load_data(); S = C.load_sim()
    P("# vc2recon_fit.py edgeu: data %d rows, sim %d rows (%.0f s)" % (D["n"], S["n"], time.time() - T0))
    # ---- join
    okg = np.all(S["g"][seq] == 0) and (S["g"] == 0).sum() == len(seq)
    dm = max(np.abs(S["m4"][seq] - V["m4"]).max(), np.abs(S["mee"][seq] - V["mee"]).max(), np.abs(S["m4t"][seq] - V["m4t"]).max())
    dw = np.abs(S["w"][seq] * S["f108"][seq] - V["W"]).max() / V["W"].max()
    P("A. join: all %d edgeu C2 rows found by sequence number: %s; max |mass diff| %.2g GeV; max rel. weight diff %.2g; bits equal %s" %
      (len(seq), okg, dm, dw, np.all(S["bits"][seq] == V["bits"])))
    # ---- nominal (edgeu_ana.py NOM)
    w_nom = C.sim_weight(S, c2corr=False)
    u_s0, mpp_s0, m4p_s0 = C.u_sim(S)
    ut_s = C.EDGE_SCALE * (S["m4t"].astype(np.float64) - C.m4max_true(S["mpipit"].astype(np.float64), S["meet"].astype(np.float64)))
    gS = S["g"]
    m = (gS <= 1) & np.isfinite(u_s0) & np.isfinite(ut_s) & (w_nom > 0) & (mpp_s0 >= 0.30) & (mpp_s0 < 0.40) & (ut_s > -0.015)
    dS, wdS = (u_s0 - ut_s)[m], w_nom[m]
    h, e = np.histogram(dS, bins=np.arange(-0.040, 0.0401, 0.002), weights=wdS)
    h2, _ = np.histogram(dS, bins=e, weights=wdS * wdS)
    xc = 0.5 * (e[1:] + e[:-1])
    pf, _ = curve_fit(lambda x, a, mu, s: a * np.exp(-0.5 * (x - mu) ** 2 / s ** 2), xc, h, p0=[h.max(), 0, 0.007],
                      sigma=np.sqrt(np.maximum(h2, 1.0)), absolute_sigma=True)
    SIG = abs(pf[2])
    grid = C.Grid(-0.200, 0.0025, 140, 4)
    sl = C.slices_of(0.30, 0.40, [0.30, 0.34, 0.37, 0.39, 0.40])
    Hd = C.hist_data(D, C.u_data(D), (D["bits"] & 1) != 0, grid)
    Hs = C.hist_sim(S, u_s0, mpp_s0, w_nom, grid)
    def simpart(H, groups, cls="all"):
        return C.sim_sum(H, groups, sl, cls)
    def nomfit(cls="all", lo=-0.045, hi=0.045, order=2, sigma=SIG, Hd_=Hd, Hs_=Hs, grid_=grid):
        sp = C.data_spec(Hd_, sl, cls, grid_, 1.0, 0, 0)
        pk = C.sim_sum(Hs_, ["C4", "C5"], sl, cls)
        return C.fit_step(sp, C.step_gauss(grid_, sigma), grid_, lo, hi, order, fixed=C.AP * pk)
    R = nomfit()
    hC3 = C.h_lin(R, simpart(Hs, ["C3"]))
    hC2raw = C.h_lin(R, simpart(Hs, ["C2"]))
    P("B. nominal: sigma(T) %.4f MeV; h_data %.2f +- %.2f, chi2 %.2f/%d; h(C2 raw) %.3f, h(C3) %.3f   [c2recon/edgeu: 6.7019, 1305.22 +- 62.27, 27.34/32, 2020.11, 1327.02]" %
      (1e3 * SIG, R.h, R.eh, R.chi2, R.ndf, hC2raw, hC3))
    # ---- own C2 template and per-row coefficients
    ub = grid.bin(u_s0[seq]); pb = C.ppbin_of(mpp_s0[seq].astype(np.float32)); cl = C.cls_of(S["mee"][seq])
    inpp = (pb >= C.PB(0.30)) & (pb < C.PB(0.40))
    s0A = ((V["bits"] & 1) != 0) & ((V["bits"] & 8) != 0)
    W = np.where(s0A, S["w"][seq] * S["f108"][seq], 0.0)
    use = s0A & inpp & (ub >= 0)
    pos = -np.ones(grid.nb, int); pos[R.idx] = np.arange(len(R.idx))
    def coef(Rf, ubins=ub, use_=use):
        p = -np.ones(grid.nb, int); p[Rf.idx] = np.arange(len(Rf.idx))
        a = np.zeros(len(seq))
        q = use_ & (p[np.clip(ubins, 0, grid.nb - 1)] >= 0) & (ubins >= 0)
        a[q] = Rf.Lh[p[ubins[q]]]
        return a
    a = coef(R)
    y_iso = np.bincount(ub[use], weights=W[use], minlength=grid.nb)
    y_ref = simpart(Hs, ["C2"])
    P("C. own C2 template vs hist_sim C2: max |diff| %.3g (sum %.4f / %.4f); h from row coefficients %.4f (h_lin %.4f)" %
      (np.abs(y_iso - y_ref).max(), y_iso.sum(), y_ref.sum(), (a * W).sum(), hC2raw))
    E1 = cl == 0; E2 = cl >= 1
    app = np.where(E1, C.CF_E1, C.CF_E2)
    h_app = (a * W * app).sum()
    P("   applied class factors 0.879 / 0.8825: h(C2) %.3f [edgeu 1780.31]; effective factor %.5f [c2recon 0.8813]" % (h_app, h_app / hC2raw))
    var = variants()
    out = {}
    def fac(x, m=None):
        mm = np.ones(len(seq), bool) if m is None else m
        return float((a * W * x)[mm].sum() / (a * W)[mm].sum())
    # group-Poisson replicas of the simulated C2 rows
    G = poisson_groups(V["tree"], V["grp"], NREP)
    def fac_rep(x, aa=a, m=None):
        mm = np.ones(len(seq), bool) if m is None else m
        num = G[:, mm] @ (aa * W * x)[mm]; den = G[:, mm] @ (aa * W)[mm]
        return num / den
    for k, x in var.items():
        c_all = fac(x); c1 = fac(x, E1); c2 = fac(x, E2)
        hq = (a * W * x).sum()
        f_new = R.h / (hq + hC3) if k not in ("ff",) else R.h / (hq + hC3)
        out[k] = dict(c=c_all, cE1=c1, cE2=c2, h=float(hq), f=float(f_new))
    rep = fac_rep(var["qed"]); rep1 = fac_rep(var["qed"], m=E1); rep2 = fac_rep(var["qed"], m=E2)
    out["qed"].update(stat=float(rep.std()), statE1=float(rep1.std()), statE2=float(rep2.std()))
    f_app = R.h / (h_app + hC3)
    f_unc = R.h / (hC2raw + hC3)
    share2 = float((a * W * E2).sum() / (a * W).sum())
    P("D. u-edge factor (all-class fit, data covariance): QED %.5f +- %.5f (group-Poisson, %d replicas); E1 %.5f +- %.5f, E2 %.5f +- %.5f; E2 share of the C2 step %.3f" %
      (out["qed"]["c"], out["qed"]["stat"], NREP, out["qed"]["cE1"], out["qed"]["statE1"], out["qed"]["cE2"], out["qed"]["statE2"], share2))
    P("   [c2recon: 0.8607 +- 0.0071; E1 0.8696 +- 0.012, E2 0.8560 +- 0.009; E2 share 0.656]")
    P("   f_eta(u-edge): uncorrected %.5f, applied factors %.5f, per-row QED %.5f (%+.2f %%)   [c2recon: 0.38995, 0.42004, 0.42574 (+1.36 %%)]" %
      (f_unc, f_app, out["qed"]["f"], 100 * (out["qed"]["f"] / f_app - 1)))
    for k in var:
        P("   variant %-9s: factor %.5f (E1 %.4f, E2 %.4f), f %.5f (%+.2f %% vs applied)" %
          (k, out[k]["c"], out[k]["cE1"], out[k]["cE2"], out[k]["f"], 100 * (out[k]["f"] / f_app - 1)))
    # ---- separate class fits
    for cls_, msk in [("E1", E1), ("E2", E2)]:
        Rc = nomfit(cls_)
        ac = coef(Rc)
        cc = float((ac * W * var["qed"])[msk].sum() / (ac * W)[msk].sum())
        P("E. class fit %s: h_data %.2f; C2 factor %.5f" % (cls_, Rc.h, cc))
        out["classfit_" + cls_] = cc
    # ---- fit variants: the constant factor that each setting would need
    for lab, kw in [("range +-35", dict(lo=-0.035, hi=0.035)), ("range +-55", dict(lo=-0.055, hi=0.055)), ("upper +30", dict(hi=0.030)),
                    ("upper +60", dict(hi=0.060)), ("pol3", dict(order=3)), ("sigma 5.76", dict(sigma=0.00576)), ("sigma 7.10", dict(sigma=0.00710))]:
        Rv = nomfit(**kw)
        av = coef(Rv)
        cv = float((av * W * var["qed"]).sum() / (av * W).sum())
        repv = fac_rep(var["qed"], aa=av)
        dpair = float(np.std(repv - rep))
        P("F. variant %-11s: C2 factor %.5f (%+.4f vs nominal, paired replica rms %.4f)" % (lab, cv, cv - out["qed"]["c"], dpair))
        out["var_" + lab] = dict(c=cv, d=cv - out["qed"]["c"], pair=dpair)
    # E(m, 0): the edge without the M(ee) term (data and simulation)
    ud0 = C.u_data(D, useq=False); us0, mp0, _ = C.u_sim(S, useq=False)
    Hd0 = C.hist_data(D, ud0, (D["bits"] & 1) != 0, grid); Hs0 = C.hist_sim(S, us0, mp0, w_nom, grid)
    R0e = nomfit(Hd_=Hd0, Hs_=Hs0)
    ub0 = grid.bin(us0[seq]); use0 = s0A & inpp & (ub0 >= 0)
    a0 = coef(R0e, ubins=ub0, use_=use0)
    c0 = float((a0 * W * var["qed"]).sum() / (a0 * W).sum())
    P("   variant E(m,0)    : C2 factor %.5f (%+.4f)" % (c0, c0 - out["qed"]["c"]))
    out["var_Em0"] = c0
    # ---- where the estimator puts its weight, and the factor in M(ee)_true bins
    u_ = u_s0[seq]
    for lo_, hi_ in [(-0.045, -0.015), (-0.015, 0.0), (0.0, 0.015), (0.015, 0.045)]:
        mm = use & (u_ >= lo_) & (u_ < hi_)
        P("G. u %+.3f..%+.3f: sum a W %+8.1f (sum |a| W %7.1f), <w> (|a| W) %.4f" %
          (lo_, hi_, (a * W)[mm].sum(), (np.abs(a) * W)[mm].sum(), (np.abs(a) * W * var["qed"])[mm].sum() / max((np.abs(a) * W)[mm].sum(), 1e-9)))
    mt = 1e3 * V["meet"]
    for lo_, hi_ in [(0, 5), (5, 10), (10, 20), (20, 40), (40, 135)]:
        mm = use & (mt >= lo_) & (mt < hi_)
        P("   M(ee)_true %3d-%3d MeV: step share %.3f, factor %.4f; window mean (all S0 A rows) %.4f; noeps %.4f" %
          (lo_, hi_, (a * W)[mm].sum() / (a * W).sum(), fac(var["qed"], mm), float((W * var["qed"])[s0A & (mt >= lo_) & (mt < hi_)].sum() / W[s0A & (mt >= lo_) & (mt < hi_)].sum()),
           fac(var["noeps"], mm) if "noeps" in var else float("nan")))
    # window means for reference
    for lab, mm in [("all M4 (S0 A)", s0A), ("-15...+5 MeV/c2 from the edge", s0A & inpp & (u_ >= -0.015) & (u_ < 0.005)),
                    ("-45...+45, 0.30-0.40", s0A & inpp & (u_ >= -0.045) & (u_ < 0.045))]:
        P("   window mean %-32s: %.4f" % (lab, float((W * var["qed"])[mm].sum() / W[mm].sum())))
    # ---- sigma from the bump (u): scales with f
    sig_app = 0.4383
    P("H. sigma from the bump (u): applied %.4f mb -> per-row QED %.4f mb (x %.4f)" % (sig_app, sig_app * out["qed"]["f"] / f_app, out["qed"]["f"] / f_app))
    RES["edgeu"] = dict(sigma_T=SIG, hd=R.h, ehd=R.eh, chi2=R.chi2, ndf=R.ndf, hC2raw=hC2raw, hC3=hC3, h_app=h_app, c_app=h_app / hC2raw,
                        f_unc=f_unc, f_app=f_app, share_E2=share2, out=out, sigma_bump_new=sig_app * out["qed"]["f"] / f_app)
    # ---- row list for bedge (nominal configuration: S0, mode A, s = 0.99615, r = 1, M(pipi) x 0.9952)
    m4 = S["m4"][seq].astype(np.float64); m4t = S["m4t"][seq].astype(np.float64)
    m4p = np.where(m4t > 0, 0.99615 * (m4t + 1.0 * (m4 - m4t)), 0.99615 * m4)
    mpp = 0.9952 * S["mpipi"][seq].astype(np.float64)
    eD = np.round((0.27 + 0.005 * np.arange(41)) * 1000.) / 1000.
    p = np.searchsorted(eD, mpp, side="right") - 1
    p[(mpp < eD[0]) | (mpp >= eD[40]) | (p >= 40)] = -1
    ib = np.floor((m4p - 0.35) / 0.0025).astype(int)
    ib[(m4p < 0.35) | (m4p >= 0.95) | (ib >= 240)] = -1
    e = np.where(S["mee"][seq].astype(np.float64) < 0.010, 0, 1)
    Wb = np.where(s0A, S["w"][seq] * S["f108"][seq].astype(np.float64), 0.0)
    names = list(var.keys())
    with open(O + "vc2_bedge_rows.txt", "w") as fo:
        fo.write("# " + " ".join(names) + "\n")
        for j in range(len(seq)):
            if Wb[j] == 0:
                continue
            fo.write("%d %d %d %.12g %s\n" % (e[j], p[j], ib[j], Wb[j], " ".join("%.9g" % var[k][j] for k in names)))
    np.savez(O + "vc2_bedge_rowbins.npz", e=e, p=p, ib=ib, Wb=Wb)
    P("I. bedge row list written: %d rows, variants %s" % ((Wb != 0).sum(), names))

# ======================================================================================================================
if MODE == "bedge":
    B = np.load(O + "vc2_bedge_rowbins.npz")
    e, p, ib, Wb = B["e"], B["p"], B["ib"], B["Wb"]
    var = variants()
    coef = {}
    nom = {}
    lines = open(O + "vc2_bedge.txt").read().splitlines()
    for l in lines:
        t = l.split()
        if t[0] == "COEF":
            k, s, pa, pb, i, a = int(t[1]), int(t[2]), int(t[3]), int(t[4]), int(t[5]), float(t[6])
            coef[(k, s, i)] = (pa, pb, a)
        else:
            P(l)
    # slice ranges
    sl = {}
    for (k, s, i), (pa, pb, a) in coef.items():
        sl[(k, s)] = (pa, pb)
    def rowcoef(kfit):
        """coefficient of every row in the fit of class kfit (0 all, 1 E1, 2 E2)"""
        a = np.zeros(len(e))
        for j in range(len(e)):
            if Wb[j] == 0 or p[j] < 0 or ib[j] < 0:
                continue
            if kfit > 0 and e[j] != kfit - 1:
                continue
            for (k, s), (pa, pb) in sl.items():
                if k == kfit and pa <= p[j] < pb:
                    a[j] = coef.get((k, s, ib[j]), (0, 0, 0.0))[2]
        return a
    a_cls = rowcoef(1) + rowcoef(2)
    a_all = rowcoef(0)
    G = poisson_groups(V["tree"], V["grp"], NREP)
    out = {}
    for lab, aa in [("class fits", a_cls), ("all-class fit", a_all)]:
        x = var["qed"]
        c = (aa * Wb * x).sum() / (aa * Wb).sum()
        r = (G @ (aa * Wb * x)) / (G @ (aa * Wb))
        m1 = e == 0; m2 = e == 1
        c1 = (aa * Wb * x)[m1].sum() / (aa * Wb)[m1].sum(); c2 = (aa * Wb * x)[m2].sum() / (aa * Wb)[m2].sum()
        r1 = (G[:, m1] @ (aa * Wb * x)[m1]) / (G[:, m1] @ (aa * Wb)[m1]); r2 = (G[:, m2] @ (aa * Wb * x)[m2]) / (G[:, m2] @ (aa * Wb)[m2])
        sh2 = (aa * Wb)[m2].sum() / (aa * Wb).sum()
        P("B. per-slice factor (%s, row coefficients): %.5f +- %.5f; E1 %.5f +- %.5f, E2 %.5f +- %.5f; E2 share %.3f   [c2recon class fits: 0.8599 +- 0.0098; E1 0.8682 +- 0.012, E2 0.8539 +- 0.0155]" %
          (lab, c, r.std(), c1, r1.std(), c2, r2.std(), sh2))
        out[lab] = dict(c=float(c), stat=float(r.std()), cE1=float(c1), sE1=float(r1.std()), cE2=float(c2), sE2=float(r2.std()), shareE2=float(sh2))
    RES["bedge"] = out
json.dump(RES, open(O + "vc2_fit_%s.json" % MODE, "w"), indent=1, default=float)
P("done in %.0f s" % (time.time() - T0))
