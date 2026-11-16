#!/usr/bin/env python3
# ============================================================================
# fix_stageA_gls.py  --  review fix (prefix fix): goodness of fit of stage A with the HYBRID GLS of the verifier v5
#   Stage A = CR-U 0.575-0.700, F10 E2, model M1 c_j (1 + b x) [G_D + kappa_B G_C], k_E2 = 1 +- 0.05, eta terms FIXED.
#   chi2 definitions at the fitted point (as v5_gls.py for F-main):
#     data only : event-group covariance of N (WP3 file) scaled to the model + CB replica covariance;
#     GLS-rep   : data + sample covariance of the 100 replica predictions (the WP3 definition; rank <= 99 for 125 bins);
#     GLS-hyb   : data + EXACT diagonal variance of the fixed raw eta templates (BB-lite) + replica covariance of the
#                 smoothed non-eta background part only.
#   Nested D11 models (M0 scale only, M1 common tilt, M2 + quadratic, M3 per-slice tilts) are compared with the same
#   hybrid covariance (built at M1, kept fixed), so that the Delta chi2 of D11 / A4 can be read with it.
#   Two eta settings: WP3 frozen (a = 0.4, kappa = rho = g = 1) and the F-main values (WP6 rerun).
#   Smoothing: the v5 library (WP4/WP5 smoothing from 0.400, about 3.5 % different from WP3 on c_j; stated).
# Uses plots/all_report/v5_fit.py as a library (read only). Writes only /hdd3/.../all_report/fix/fix_stageA_gls.log.
# Run: python3 -I plots/all_report/fix_stageA_gls.py   (from /home/przygoda/HADES/FAT/FAT; a few minutes)
# ============================================================================
import sys, math
import numpy as np
from scipy.stats import chi2 as CHI2
sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import v5_fit as F

OUT = "/hdd3/przygoda/HADES/FAT_out/all_report/fix/"
LOG = open(OUT + "fix_stageA_gls.log", "w")
def P(*a):
    s = " ".join(str(x) for x in a); print(s); LOG.write(s + "\n"); LOG.flush()

def make_model(eta, kind):
    a, ka, rho, g = eta
    co = {"C2": a, "C3": a * ka, "C4": a * ka, "C5": a * rho, "Geta": a * g}
    def m(th, Bx):
        c = th[0:5][Bx.j]
        if kind == "M0":
            tl = np.ones(Bx.n); kB = th[5]; k = th[6]
        elif kind == "M1":
            tl = 1. + th[5] * Bx.x; kB = th[6]; k = th[7]
        elif kind == "M2":
            tl = 1. + th[5] * Bx.x + th[6] * Bx.x ** 2; kB = th[7]; k = th[8]
        else:  # M3 per-slice tilts
            tl = 1. + th[5:10][Bx.j] * Bx.x; kB = th[10]; k = th[11]
        mu = k * Bx.cb + c * tl * (Bx.gd + kB * Bx.gc) + sum(co[t] * Bx.s[t] for t in F.SIG)
        vs = sum(co[t] ** 2 * Bx.vs[t] for t in F.SIG)
        return mu, vs, k
    return m

def start(kind):
    if kind == "M0": return np.array([1, 1, 1, 1, 1, 3.0, 1.0]), [(6, 1.0, 0.05)]
    if kind == "M1": return np.array([1, 1, 1, 1, 1, -0.07, 3.0, 1.0]), [(7, 1.0, 0.05)]
    if kind == "M2": return np.array([1, 1, 1, 1, 1, -0.07, 0.0, 3.0, 1.0]), [(8, 1.0, 0.05)]
    return np.array([1, 1, 1, 1, 1, -0.07, -0.07, -0.07, -0.07, -0.07, 3.0, 1.0]), [(11, 1.0, 0.05)]

def bg_part(th, Bx, kind):
    c = th[0:5][Bx.j]
    if kind == "M1": return c * (1. + th[5] * Bx.x) * (Bx.gd + th[6] * Bx.gc)
    raise ValueError

def main():
    T, _ = F.load_sim("S0"); X = F.load_data("S0")
    S, _ = F.smooth_all(T, lambda c, h: T[(c, h)])
    hE2 = ["F10_S%dE2" % j for j in range(1, 6)]
    B = F.make_bins(X, T, S, hE2, 0.575, 0.700)
    P("stage-A bins: %d (expect 125)" % B.n)
    # data covariance: event groups (WP3 file, 700 x 700) + CB replicas
    eg3 = np.loadtxt(F.D + "egcov_wp3_S0_F10.txt")
    hidx = {h: q for q, h in enumerate(F.F10)}
    cells = np.array([hidx[B.h[q]] * 70 + B.i[q] for q in range(B.n)])
    eg = eg3[np.ix_(cells, cells)]
    idx = np.where((F.edges[:-1] >= 0.575 - 1e-9) & (F.edges[1:] <= 0.700 + 1e-9))[0]
    assert all((B.i[q] == idx[q % len(idx)]) for q in range(B.n))
    cbrep = np.concatenate([X[("CBrep", h)][:, idx] for h in hE2], axis=1)
    cbcov = np.cov(cbrep, rowvar=False, ddof=1)
    # replicas (smoothing redone per replica)
    R = F.load_simrep("S0")
    REPS = []
    for q in range(F.NREP):
        Tr = {}
        for c in F.COMPS:
            for h in F.HN: Tr[(c, h)] = R[(c, h)][q]; Tr[(c, h, "v")] = T[(c, h, "v")]
            Tr[(c, "win")] = R[(c, "win")][q]; Tr[(c, "tot")] = R[(c, "tot")][q]
        Sr, _ = F.smooth_all(Tr, lambda c, h: Tr[(c, h)])
        REPS.append(F.swap_templates(B, Tr, Sr, T))
    P("replicas smoothed: %d" % len(REPS))
    settings = (("WP3 frozen eta (a 0.4, kappa = rho = g = 1)", (0.4, 1.0, 1.0, 1.0)),
                ("F-main eta (a 0.3823, kappa 0.8323, rho 1.3152, g 1.7168)", (0.38228, 0.83234, 1.31520, 1.71677)))
    for tag, eta in settings:
        P("\n=== %s" % tag)
        m1 = make_model(eta, "M1"); th0, cons = start("M1")
        th, _ = F.fit_irls(B, m1, cons, th0)
        mu, vs, k = m1(th, B)
        pos = B.y > 0
        sf = np.where(pos, np.sqrt(np.maximum(mu, 1e-9) / np.where(pos, B.y, 1.)), 0.)
        Cd = eg * np.outer(sf, sf) + np.diag(np.where(pos, 0., np.maximum(mu, 1.))) + k * k * cbcov
        preds, bgs = [], []
        for Br in REPS:
            preds.append(m1(th, Br)[0]); bgs.append(bg_part(th, Br, "M1"))
        preds = np.array(preds); bgs = np.array(bgs)
        Crep = np.cov(preds, rowvar=False, ddof=1); Cbg = np.cov(bgs, rowvar=False, ddof=1)
        r = B.y - mu
        pen = sum(((th[kk] - v) / s) ** 2 for (kk, v, s) in cons)
        def chi(C, rr=r, pp=pen): return rr @ np.linalg.solve(C, rr) + pp
        ndf = B.n - len(th) + len(cons)
        P("M1 parameters: " + " ".join("%.4f" % v for v in th))
        out = [("data only", chi(Cd)), ("Pearson-type (diag data + BB-lite)", chi(np.diag(np.diag(Cd) + vs))),
               ("GLS-rep (WP3 definition: data + replica cov of the prediction)", chi(Cd + Crep)),
               ("GLS-rep, 50 replicas", chi(Cd + np.cov(preds[:50], rowvar=False, ddof=1))),
               ("GLS-hyb (data + diag BB-lite eta + replica cov of the background)", chi(Cd + np.diag(vs) + Cbg)),
               ("GLS-hyb, background replica cov from 50 replicas", chi(Cd + np.diag(vs) + np.cov(bgs[:50], rowvar=False, ddof=1)))]
        for nm, v in out:
            P("  %-72s chi2 = %7.1f / %d = %.3f  (p = %.2g)" % (nm, v, ndf, v / ndf, CHI2.sf(v, ndf)))
        evr = np.linalg.eigvalsh(Crep); evb = np.linalg.eigvalsh(Cbg)
        P("  rank: replica prediction cov %d, background replica cov %d (of %d bins); trace BB-lite eta %.1f, background cov %.1f, prediction cov %.1f"
          % ((evr > 1e-9 * evr.max()).sum(), (evb > 1e-9 * evb.max()).sum(), B.n, vs.sum(), np.trace(Cbg), np.trace(Crep)))
        # nested models with the hybrid covariance fixed at M1 (GLS minimisation would change the weights; the fixed
        # covariance is the quantity used for the D11 Delta chi2)
        Chyb = Cd + np.diag(vs) + Cbg
        Ldec = np.linalg.cholesky(np.linalg.inv(Chyb))
        res = {}
        for kind in ("M0", "M1", "M2", "M3"):
            mk = make_model(eta, kind); t0, cs = start(kind)
            # minimise the hybrid GLS chi2 with the covariance fixed
            from scipy.optimize import least_squares
            def resid(t):
                rr = B.y - mk(t, B)[0]
                return np.concatenate([Ldec.T @ rr, [(t[kk] - v) / s for (kk, v, s) in cs]])
            sol = least_squares(resid, t0, method="trf", x_scale="jac", xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=20000)
            c2 = float((sol.fun ** 2).sum()); nd = B.n - len(t0) + len(cs)
            res[kind] = (c2, nd)
            P("  %s minimised hybrid GLS (cov fixed at M1): chi2 = %7.1f / %d = %.3f" % (kind, c2, nd, c2 / nd))
        P("  Delta chi2 (hybrid): M0->M1 %.1f, M1->M2 %.1f, M1->M3 %.1f (4 dof, p = %.3g)"
          % (res["M0"][0] - res["M1"][0], res["M1"][0] - res["M2"][0], res["M1"][0] - res["M3"][0], CHI2.sf(res["M1"][0] - res["M3"][0], 4)))

if __name__ == "__main__":
    main()
