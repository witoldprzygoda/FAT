#!/usr/bin/env python3
# ============================================================================
# v5_quad.py -- does the data prefer the quadratic background term (WP5 variant FMAIN_quad, a_eta +39 %)?  prefix v5.
# Nominal (M1: c_j (1 + b x)) and quad (c_j (1 + b x + q x^2)) F-main E2 fits with the own IRLS; goodness of fit with
#   Pearson (V = BB-lite weights), GLS-rep (WP5 definition: replica covariance of the full prediction) and GLS-hyb
#   (exact diagonal BB-lite variance of the raw signal templates + replica covariance of the smoothed background part),
#   each at the fit point and minimised (Gauss-Newton with the covariance frozen, 3 refreshes).
# Run: python3 plots/all_report/v5_quad.py   (from /home/przygoda/HADES/FAT/FAT; about 5 min)
# ============================================================================
import sys, math
import numpy as np
sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import v5_fit as F

LOG = open(F.D + "v5_quad2.log", "w")
def P(*a):
    s = " ".join(str(x) for x in a); print(s); LOG.write(s + "\n"); LOG.flush()

NAMES = F.MAIN + ["q"]
def model_q(th, B):
    a, ka, rho, g = th[0], th[1], th[2], th[3]
    c = th[4:9][B.j]; b = th[9]; kB = th[10]; k = th[11]; q = th[12]
    tl = 1. + b * B.x + q * B.x * B.x
    co = {"C2": a, "C3": a * ka, "C4": a * ka, "C5": a * rho, "Geta": a * g}
    mu = k * B.cb + c * tl * (B.gd + kB * B.gc) + sum(co[t] * B.s[t] for t in F.SIG)
    vs = sum(co[t] ** 2 * B.vs[t] for t in F.SIG)
    return mu, vs, k
def bg_part(th, B, gd, gc):
    q = th[12] if len(th) > 12 else 0.
    c = th[4:9][B.j]
    return c * (1. + th[9] * B.x + q * B.x * B.x) * (gd + th[10] * gc)

def main():
    T, _ = F.load_sim("S0"); X = F.load_data("S0")
    S, _ = F.smooth_all(T, lambda c, h: T[(c, h)])
    hE2 = ["F10_S%dE2" % j for j in range(1, 6)]
    B = F.make_bins(X, T, S, hE2)
    th0 = np.array([0.4, 1.0, 1.0, 1.0, 1, 1, 1, 1, 1, -0.07, 2.0, 1.0])
    thN, _ = F.fit_irls(B, F.model_main, F.CONS_MAIN, th0)
    thQ, _ = F.fit_irls(B, model_q, F.CONS_MAIN, np.concatenate([thN, [0.0]]))
    EN = F.errors(B, F.model_main, F.CONS_MAIN, thN); EQ = F.errors(B, model_q, F.CONS_MAIN, thQ)
    P("nominal: a_eta %.4f Y_bump %.0f | Pearson %.1f" % (thN[0], F.derived_main(thN, T)["Y_bump"], EN["chi2P"]))
    P("quad   : a_eta %.4f Y_bump %.0f q = %.4f +- %.4f (Fisher) | Pearson %.1f" % (thQ[0], F.derived_main(thQ, T)["Y_bump"], thQ[12], math.sqrt(EQ["C"][12, 12]), EQ["chi2P"]))
    # replicas: raw signal and smoothed background templates in the E2 bins
    R = F.load_simrep("S0")
    reps = []
    for r in range(100):
        Tr = {}
        for c in F.COMPS:
            for h in F.HN: Tr[(c, h)] = R[(c, h)][r]; Tr[(c, h, "v")] = T[(c, h, "v")]
            Tr[(c, "win")] = R[(c, "win")][r]; Tr[(c, "tot")] = R[(c, "tot")][r]
        Sr, _ = F.smooth_all(Tr, lambda c, h: Tr[(c, h)])
        reps.append(F.swap_templates(B, Tr, Sr, T))
    eg3 = np.loadtxt(F.D + "egcov_wp3_S0_F10.txt")
    hidx = {h: q for q, h in enumerate(F.F10)}
    cells = np.array([hidx[B.h[q]] * 70 + B.i[q] for q in range(B.n)])
    eg = eg3[np.ix_(cells, cells)]
    cbcov = np.cov(np.concatenate([X[("CBrep", h)] for h in hE2], axis=1), rowvar=False, ddof=1)
    def covs(th, model):
        mu, vs, k = model(th, B)
        pos = B.y > 0
        sf = np.where(pos, np.sqrt(np.maximum(mu, 1e-9) / np.where(pos, B.y, 1.)), 0.)
        Cd = eg * np.outer(sf, sf) + np.diag(np.where(pos, 0., np.maximum(mu, 1.))) + k * k * cbcov
        preds = np.array([model(th, Br)[0] for Br in reps])
        bgs = np.array([bg_part(th, Br, Br.gd, Br.gc) for Br in reps])
        return Cd + np.cov(preds, rowvar=False, ddof=1), Cd + np.diag(vs) + np.cov(bgs, rowvar=False, ddof=1)
    def chi(th, model, W):
        r = B.y - model(th, B)[0]
        return r @ np.linalg.solve(W, r) + sum(((th[kk] - v) / s) ** 2 for (kk, v, s) in F.CONS_MAIN)
    def gls_min(th, model, which):
        th = th.copy()
        for outer in range(3):
            W = covs(th, model)[which]; Wi = np.linalg.inv(W)
            for it in range(40):
                J, mu = F.jac(B, model, th)
                r = B.y - mu
                A = J.T @ Wi @ J; bb = J.T @ Wi @ r
                for (kk, v, s) in F.CONS_MAIN: A[kk, kk] += 1 / s ** 2; bb[kk] -= (th[kk] - v) / s ** 2
                d = np.linalg.solve(A, bb); th = th + d
                if np.max(np.abs(d) / np.sqrt(np.diag(np.linalg.inv(A)))) < 1e-6: break
        return th, chi(th, model, W)
    for nm, th, model, ndf in (("nominal", thN, F.model_main, 341), ("quad", thQ, model_q, 340)):
        Wr, Wh = covs(th, model)
        cr, ch = chi(th, model, Wr), chi(th, model, Wh)
        tmr, cmr = gls_min(th, model, 0); tmh, cmh = gls_min(th, model, 1)
        P("%-8s GLS-rep at fit %.1f, minimised %.1f (a_eta %.4f) | GLS-hyb at fit %.1f (%.3f/ndf), minimised %.1f (a_eta %.4f%s) | ndf %d"
          % (nm, cr, cmr, tmr[0], ch, ch / ndf, cmh, tmh[0], (", q %.4f" % tmh[12]) if len(tmh) > 12 else "", ndf))
        for tg, t in (("IRLS", th), ("GLS-rep min", tmr), ("GLS-hyb min", tmh)):
            d = F.derived_main(t, T)
            P("   %-8s %-12s a_eta %.4f kappa %.3f rho_d %.3f g %.3f b %.4f kappa_B %.3f | Y_bump %.0f Y_peak %.0f DR %.3f"
              % (nm, tg, t[0], t[1], t[2], t[3], t[9], t[10], d["Y_bump"], d["Y_peak"], d["DR"]))

if __name__ == "__main__":
    main()
