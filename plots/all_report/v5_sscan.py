#!/usr/bin/env python3
# ============================================================================
# v5_sscan.py  --  independent check of the WP5 mass-transfer test ("free s within 2 sigma of 0.99615") and of the
#                  r (resolution) statement; prefix v5. Uses the own fitter/smoothing of v5_fit.py (imported).
# For every WP0 template variant (sK0 0.99521, sDn 0.99525, nom 0.99615, sUp 0.99705; r089, r110 at s = 0.99615):
#   (a) whole F-main E2 fit, fixed weights = V of the nominal fit;
#   (b) peak region 0.5225-0.565 (E2): a_eta, kappa, g, b, kappa_B fixed to the nominal F-main; rho_d, c_j, k_E2 free;
#       fixed weights = V of the nominal local fit (as WP5) and, as a check, V of the nominal F-main;
#   (c) peak region with ALL eta and background parameters free except kappa_B and b (to see whether fixing a_eta/g
#       drives the preference).
# Run: python3 plots/all_report/v5_sscan.py   (from /home/przygoda/HADES/FAT/FAT; about 1 min)
# ============================================================================
import sys, math
import numpy as np
from scipy.optimize import least_squares
sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import v5_fit as F

LOG = open(F.D + "v5_sscan.log", "w")
def P(*a):
    s = " ".join(str(x) for x in a); print(s); LOG.write(s + "\n"); LOG.flush()

def sub_model(model, th_full, free):
    free = list(free)
    def m(t, B):
        th = th_full.copy(); th[free] = t
        return model(th, B)
    return m

def fit_fixed(B, model, cons_full, th_full, free, Vfix):
    free = list(free)
    cons = [(free.index(k), v, s) for (k, v, s) in cons_full if k in free]
    mdl = sub_model(model, th_full, free)
    def res(t):
        return np.concatenate([(B.y - mdl(t, B)[0]) / np.sqrt(Vfix), [(t[k] - v) / s for (k, v, s) in cons]])
    r = least_squares(res, th_full[free], method="trf", x_scale="jac", xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=5000)
    th = th_full.copy(); th[free] = r.x
    return th, 2 * r.cost

def main():
    T, _ = F.load_sim("S0"); X = F.load_data("S0")
    S, _ = F.smooth_all(T, lambda c, h: T[(c, h)])
    hE2 = ["F10_S%dE2" % j for j in range(1, 6)]
    B = F.make_bins(X, T, S, hE2)
    th0 = np.array([0.4, 1.0, 1.0, 1.0, 1, 1, 1, 1, 1, -0.07, 2.0, 1.0])
    thI, _ = F.fit_irls(B, F.model_main, F.CONS_MAIN, th0)
    Vn, _ = F.V_of(B, thI, F.model_main)
    P("nominal F-main a_eta %.5f" % thI[0])
    sv = F.read_rows(F.D + "sim_S0_svar.txt")
    def templ(var):
        if var == "nom": return T
        Tv = {}
        for c in F.COMPS:
            for h in F.HN:
                Tv[(c, h)] = sv[var + ":" + c + "_" + h] * F.NB2RAW; Tv[(c, h, "v")] = (sv[var + ":" + c + "err_" + h] * F.NB2RAW) ** 2
        return Tv
    svals = {"sK0": 0.99521, "sDn": 0.99525, "nom": 0.99615, "sUp": 0.99705, "r089": 0.99615, "r110": 0.99615}
    # local nominal fit (peak region) and its fixed point V (IRLS restricted)
    Bp = F.make_bins(X, T, S, hE2, 0.5225, 0.565)
    freeL = [2, 4, 5, 6, 7, 8, 11]   # rho_d, c_1..c_5, k_E2
    mL = sub_model(F.model_main, thI, freeL)
    consL = [(freeL.index(k), v, s) for (k, v, s) in F.CONS_MAIN if k in freeL]
    tL, _ = F.fit_irls(Bp, mL, consL, thI[freeL])
    thL = thI.copy(); thL[freeL] = tL
    VL, _ = F.V_of(Bp, thL, F.model_main)
    # V of the nominal F-main in the same bins
    selp = np.array([(F.edges[i] >= 0.5225 - 1e-9) and (F.edges[i + 1] <= 0.565 + 1e-9) for i in B.i])
    VnP = Vn[selp]
    freeC = [0, 1, 2, 3, 4, 5, 6, 7, 8, 11]   # (c): b and kappa_B fixed, the rest free
    res = {}
    for var in ("sK0", "sDn", "nom", "sUp", "r089", "r110"):
        Tv = templ(var); Sv, _ = F.smooth_all(Tv, lambda c, h: Tv[(c, h)])
        Bw = F.make_bins(X, Tv, Sv, hE2)
        tw, cw = fit_fixed(Bw, F.model_main, F.CONS_MAIN, thI, range(12), Vn)
        Bv = F.make_bins(X, Tv, Sv, hE2, 0.5225, 0.565)
        tb, cb = fit_fixed(Bv, F.model_main, F.CONS_MAIN, thI, freeL, VL)
        tb2, cb2 = fit_fixed(Bv, F.model_main, F.CONS_MAIN, thI, freeL, VnP)
        tc, cc = fit_fixed(Bv, F.model_main, F.CONS_MAIN, thI, freeC, VL)
        res[var] = (cw, cb, cb2, cc)
        P("%-5s s %.5f | whole fit chi2 %.2f a_eta %.4f | peak (WP5 config, local V) chi2 %.2f rho_d %.4f | peak (F-main V) %.2f | peak (eta free) %.2f a_eta %.4f rho_d %.4f"
          % (var, svals[var], cw, tw[0], cb, tb[2], cb2, cc, tc[0], tc[2]))
    def parab(keys, idx):
        s = np.array([svals[k] for k in keys]); c = np.array([res[k][idx] for k in keys])
        A, Bq, C = np.polyfit(s - 0.99615, c, 2)
        sh = -Bq / (2 * A); sg = 1 / math.sqrt(A) if A > 0 else float("nan")
        return 0.99615 + sh, sg
    for idx, nm in ((0, "whole fit"), (1, "peak, WP5 config"), (2, "peak, F-main V"), (3, "peak, eta free")):
        s3, e3 = parab(["sDn", "nom", "sUp"], idx); s4, e4 = parab(["sK0", "sDn", "nom", "sUp"], idx)
        d = res["nom"][idx] - res["sDn"][idx]
        P("%-18s 3-point parabola s_hat %.5f +- %.5f (%.2f sigma) | 4-point %.5f +- %.5f (%.2f sigma) | chi2(nom) - chi2(sDn) = %+.2f, chi2(sK0) - chi2(sDn) = %+.2f"
          % (nm, s3, e3, (s3 - 0.99615) / e3, s4, e4, (s4 - 0.99615) / e4, d, res["sK0"][idx] - res["sDn"][idx]))
    for idx, nm in ((0, "whole fit"), (1, "peak, WP5 config"), (3, "peak, eta free")):
        P("r test %-18s chi2(r089) - chi2(nom) = %+.2f ; chi2(r110) - chi2(nom) = %+.2f" % (nm, res["r089"][idx] - res["nom"][idx], res["r110"][idx] - res["nom"][idx]))

if __name__ == "__main__":
    main()
