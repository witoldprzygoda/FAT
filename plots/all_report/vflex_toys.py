#!/usr/bin/env python3
# ============================================================================
# vflex_toys.py -- independent verification of package flex: bias test of the nominal model F8 (one configuration)
# ============================================================================
# Own implementation of the WP4 "realistic fit" (RF) test with own random numbers (numpy PCG64, seed 20261115 + toy):
#   truth   : mu_T = k_T CB_data + sum_t s_T,t T_t(full sample), scales = the independent F8 fit of the data (scenario 0)
#             or all 1 (scenario 1);
#   toy     : N ~ Poisson(mu_T); constraint centre k_T + 0.05 G; measured CB = data CB replica (toy mod 200);
#   RF fit  : templates of bootstrap replica (toy mod 100), smoothing redone; method fit: full-sample templates;
#   errors  : sandwich (Poisson = diag mu, CB replica covariance, constraint) (+) sigma_MC (RMS of the 100 replica fits of the
#             Asimov data), as flex_toys.C;
#   plus    : the hybrid GLS chi2 of every RF toy fit (goodness-of-fit reference for the data value 603.6).
# Run: python3 -B plots/all_report/vflex_toys.py <scenario> <ntoy>     (6 worker processes)
# Output: /hdd3/przygoda/HADES/FAT_out/all_report/vflex/vflex_toys_s<scen>.npz, vflex_toys_s<scen>.log
# ============================================================================
import sys, time
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import vflex_core as V

scen = int(sys.argv[1]) if len(sys.argv) > 1 else 0
ntoy = int(sys.argv[2]) if len(sys.argv) > 2 else 500
D = V.load()
L = V.Layout(D, (1,))
M = V.FlexModel(["C2", "C3", "C4", "C5", "Geta", "GD_om", "GD_nom", "GC"])
T, Vt = V.templates(D, L, M)
TR = [V.templates(D, L, M, rep=r)[0] for r in range(V.NREPS)]
Rd = V.fit(M, L, T, Vt, egcov=L.eg)
truth = Rd['th'].copy() if scen == 0 else M.start.copy()
muT = M.eval(truth, L, T, Vt)[0]
consT = {M.np - 1: (truth[-1], 0.05)}
# Asimov replica fits -> sigma_MC and the replica shift
A = V.replica_fits(M, L, TR, Vt, truth, y=muT, cons=consT)
sigmc = np.nanstd(A, 0, ddof=1); ashift = np.nanmean(A, 0) - truth
RA = V.fit(M, L, T, Vt, y=muT, cons=consT, start=truth, egcov=None)


def one(t):
    rng = np.random.Generator(np.random.PCG64(20261115 + 7919 * scen + t))
    y = rng.poisson(muT).astype(float)
    cons = {M.np - 1: (truth[-1] + 0.05 * rng.standard_normal(), 0.05)}
    cb = L.cbrep[:, t % V.NREPC]
    out = []
    for mode in (0, 1):
        Tm = TR[t % V.NREPS] if mode == 0 else T
        R = V.fit(M, L, Tm, Vt, y=y, cb=cb, start=truth, cons=cons, egcov=None)
        g = np.nan
        if R['ok'] and mode == 0:
            import copy
            Lt = copy.copy(L); Lt.cb = cb
            Mt = copy.copy(M); Mt.cons = cons
            g = V.hybrid_gls(Mt, Lt, Tm, Vt, R, TR)[0]
        out.append((R['ok'], R['th'], np.sqrt(np.diag(R['Cstat'])), R['pearson'], g))
    return out


if __name__ == "__main__":
    t0 = time.time()
    with Pool(6) as p:
        res = p.map(one, range(ntoy), chunksize=10)
    ok = np.array([[r[m][0] for m in (0, 1)] for r in res])
    th = np.array([[r[m][1] for m in (0, 1)] for r in res])
    es = np.array([[r[m][2] for m in (0, 1)] for r in res])
    pe = np.array([[r[m][3] for m in (0, 1)] for r in res])
    gl = np.array([r[0][4] for r in res])
    np.savez(f"{V.DIR}/vflex_toys_s{scen}.npz", ok=ok, th=th, es=es, pe=pe, gl=gl, truth=truth, sigmc=sigmc, ashift=ashift, names=M.names)
    log = open(f"{V.DIR}/vflex_toys_s{scen}.log", "w")

    def P(s):
        print(s); log.write(s + "\n")
    P(f"scenario {scen}: {ntoy} toys, RF ok {ok[:, 0].sum()}, method ok {ok[:, 1].sum()} ({time.time() - t0:.0f} s); Asimov full templates returns truth to {np.max(np.abs(RA['th'] - truth)):.2e}")
    P(f"RF hybrid GLS: mean {np.nanmean(gl):.1f} rms {np.nanstd(gl, ddof=1):.1f} quantiles 5/50/95 % {np.nanpercentile(gl, 5):.1f} {np.nanpercentile(gl, 50):.1f} {np.nanpercentile(gl, 95):.1f} max {np.nanmax(gl):.1f}; n >= 603.55: {(gl >= 603.55).sum()}")
    P(f"RF Pearson mean {pe[:, 0].mean():.1f} rms {pe[:, 0].std(ddof=1):.1f}; method Pearson mean {pe[:, 1].mean():.1f} rms {pe[:, 1].std(ddof=1):.1f}")
    P(f"{'param':10s} {'truth':>8s} {'s_stat':>7s} {'s_MC':>7s} {'MCsh%':>7s} {'RFbias%':>8s} {'+-':>6s} {'b/stot':>7s} {'pullw':>6s} {'cov68':>6s} {'meth%':>7s} {'methw':>6s}")
    for k, n in enumerate(M.names):
        g = ok[:, 0]
        v = th[g, 0, k]; e = es[g, 0, k]
        et = np.sqrt(e ** 2 + sigmc[k] ** 2)
        pull = (v - truth[k]) / et
        bias = v.mean() - truth[k]
        berr = np.sqrt(v.var(ddof=1) / len(v) + sigmc[k] ** 2 / 100)
        sst = e.mean(); stot = np.sqrt(sst ** 2 + sigmc[k] ** 2)
        gm = ok[:, 1]
        vm = th[gm, 1, k]; em = es[gm, 1, k]
        rel = 100 * bias / truth[k] if truth[k] != 0 else np.nan
        P(f"{n:10s} {truth[k]:8.4f} {sst:7.4f} {sigmc[k]:7.4f} {100 * ashift[k] / truth[k]:7.1f} {rel:8.1f} {100 * berr / abs(truth[k]):6.1f} {bias / stot:7.2f} {pull.std(ddof=1):6.2f} {np.mean(np.abs(pull) < 1):6.3f} "
          f"{100 * (vm.mean() - truth[k]) / truth[k]:7.2f} {((vm - truth[k]) / em).std(ddof=1):6.2f}")
