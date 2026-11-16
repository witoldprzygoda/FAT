#!/usr/bin/env python3
# ============================================================================
# v5_gls.py  --  independent check of the WP5 goodness-of-fit claim (GLS chi2 883.7/341 = 2.59) and of the
#                rank of its MC covariance; prefix v5. Uses the own smoothing/fitter of v5_fit.py.
# chi2 variants at the nominal F-main (E2) fixed point theta_hat:
#   GLS-rep   : data covariance (event groups scaled to the model + CB replicas) + sample covariance of the 100 replica
#               predictions (WP5 definition; rank <= 99 for 350 bins);
#   GLS-hyb   : data covariance + EXACT diagonal variance of the raw signal templates (sum w^2, the BB-lite term) +
#               replica covariance of the smoothed non-eta background part only (low rank: 4 shape parameters x 10 fits);
#   GLS-half  : GLS-rep with the replica covariance from 50 replicas, to show the dependence on the replica number.
# Run: python3 plots/all_report/v5_gls.py   (from /home/przygoda/HADES/FAT/FAT; about 4 min)
# ============================================================================
import sys, math
import numpy as np
from scipy.stats import chi2 as CHI2
sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import v5_fit as F

LOG = open(F.D + "v5_gls.log", "w")
def P(*a):
    s = " ".join(str(x) for x in a); print(s); LOG.write(s + "\n"); LOG.flush()

def main():
    T, _ = F.load_sim("S0"); X = F.load_data("S0")
    S, _ = F.smooth_all(T, lambda c, h: T[(c, h)])
    hE2 = ["F10_S%dE2" % j for j in range(1, 6)]
    B = F.make_bins(X, T, S, hE2)
    th0 = np.array([0.4, 1.0, 1.0, 1.0, 1, 1, 1, 1, 1, -0.07, 2.0, 1.0])
    th, _ = F.fit_irls(B, F.model_main, F.CONS_MAIN, th0)
    mu, vs, k = F.model_main(th, B)
    eg3 = np.loadtxt(F.D + "egcov_wp3_S0_F10.txt")
    hidx = {h: q for q, h in enumerate(F.F10)}
    cells = np.array([hidx[B.h[q]] * 70 + B.i[q] for q in range(B.n)])
    eg = eg3[np.ix_(cells, cells)]
    cbcov = np.cov(np.concatenate([X[("CBrep", h)] for h in hE2], axis=1), rowvar=False, ddof=1)
    pos = B.y > 0
    sf = np.where(pos, np.sqrt(np.maximum(mu, 1e-9) / np.where(pos, B.y, 1.)), 0.)
    Cd = eg * np.outer(sf, sf) + np.diag(np.where(pos, 0., np.maximum(mu, 1.))) + k * k * cbcov
    r = B.y - mu
    pen = sum(((th[kk] - v) / s) ** 2 for (kk, v, s) in F.CONS_MAIN)
    R = F.load_simrep("S0")
    preds, bgs = [], []
    for q in range(100):
        Tr = {}
        for c in F.COMPS:
            for h in F.HN: Tr[(c, h)] = R[(c, h)][q]; Tr[(c, h, "v")] = T[(c, h, "v")]
            Tr[(c, "win")] = R[(c, "win")][q]; Tr[(c, "tot")] = R[(c, "tot")][q]
        Sr, _ = F.smooth_all(Tr, lambda c, h: Tr[(c, h)])
        Br = F.swap_templates(B, Tr, Sr, T)
        preds.append(F.model_main(th, Br)[0])
        c_ = th[4:9][Br.j]; bg = c_ * (1. + th[9] * Br.x) * (Br.gd + th[10] * Br.gc)
        bgs.append(bg)
    preds = np.array(preds); bgs = np.array(bgs)
    np.save(F.D + "v5_gls_preds.npy", preds); np.save(F.D + "v5_gls_bgs.npy", bgs)
    Crep = np.cov(preds, rowvar=False, ddof=1)
    Cbg = np.cov(bgs, rowvar=False, ddof=1)
    def chi(C): return r @ np.linalg.solve(C, r) + pen
    ndf = 341
    out = {}
    out["data only (full)"] = chi(Cd)
    out["Pearson-type: data diag + BB-lite diag"] = chi(np.diag(np.diag(Cd) + vs))
    out["GLS-rep (WP5 definition)"] = chi(Cd + Crep)
    out["GLS-rep, 50 replicas"] = chi(Cd + np.cov(preds[:50], rowvar=False, ddof=1))
    out["GLS-hyb: data + diag(BB-lite signal) + replica cov(background)"] = chi(Cd + np.diag(vs) + Cbg)
    out["GLS-hyb with diagonal data covariance"] = chi(np.diag(np.diag(Cd)) + np.diag(vs) + Cbg)
    out["GLS-hyb, background replica cov from 50 replicas"] = chi(Cd + np.diag(vs) + np.cov(bgs[:50], rowvar=False, ddof=1))
    out["GLS-hyb, background replica cov from replicas 50-99"] = chi(Cd + np.diag(vs) + np.cov(bgs[50:], rowvar=False, ddof=1))
    # signal noise with its event-group (off-diagonal) structure: replica covariance of the signal part, shrunk to the
    # exact diagonal (diag replaced by the BB-lite variances, off-diagonal kept)
    sig = preds - bgs
    Csig = np.cov(sig, rowvar=False, ddof=1)
    Csig_od = Csig - np.diag(np.diag(Csig)) + np.diag(vs)
    w_, v_ = np.linalg.eigh(Csig_od); Csig_od = (v_ * np.maximum(w_, 0.)) @ v_.T
    out["GLS-hyb2: data + signal (exact diag, replica off-diag) + background"] = chi(Cd + Csig_od + Cbg)
    tot = sig.sum(axis=1)
    P("signal part summed over the 350 bins: replica variance %.1f vs sum of diagonal BB-lite variances %.1f (ratio %.3f; > 1 = positive bin-bin correlations)" % (tot.var(ddof=1), vs.sum(), tot.var(ddof=1) / vs.sum()))
    for kk, v in out.items():
        P("%-66s chi2 = %8.1f / %d = %.3f  (p = %.2g)" % (kk, v, ndf, v / ndf, CHI2.sf(v, ndf)))
    ev = np.linalg.eigvalsh(Crep); evb = np.linalg.eigvalsh(Cbg)
    P("rank (eigenvalues > 1e-9 x max): replica prediction covariance %d, background-only %d (of %d bins)" % ((ev > 1e-9 * ev.max()).sum(), (evb > 1e-9 * evb.max()).sum(), B.n))
    P("trace: replica prediction cov %.1f; diag BB-lite signal %.1f + background replica cov %.1f" % (np.trace(Crep), vs.sum(), np.trace(Cbg)))
    # replica-sampled signal noise vs exact BB-lite: compare the trace of the signal part of the replica covariance
    sig = preds - bgs - k * B.cb[None, :]
    Csig = np.cov(sig, rowvar=False, ddof=1)
    P("signal-template part: trace of the replica covariance %.1f vs exact sum of BB-lite variances %.1f" % (np.trace(Csig), vs.sum()))

if __name__ == "__main__":
    main()
