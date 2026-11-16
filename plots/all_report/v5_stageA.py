#!/usr/bin/env python3
# ============================================================================
# v5_stageA.py  --  independent stage-A refit (prefix v5): CR-U 0.575-0.700, F10 E2, model M1 c_j (1 + b x)
#                   [G_D + kappa_B G_C], k_E2 = 1 +- 0.05; eta terms FIXED. Purpose: (1) reproduce the frozen WP3
#                   values with the own code; (2) the "rerun at the WP5 value" that WP3/WP5 did not do (eta fixed to the
#                   F-main a_eta, kappa, rho_d, g). Does not write into the WP3 directory.
# Run: python3 plots/all_report/v5_stageA.py   (from /home/przygoda/HADES/FAT/FAT; seconds)
# ============================================================================
import sys
import numpy as np
sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import v5_fit as F

LOG = open(F.D + "v5_stageA.log", "w")
def P(*a):
    s = " ".join(str(x) for x in a); print(s); LOG.write(s + "\n"); LOG.flush()

def main():
    T, _ = F.load_sim("S0"); X = F.load_data("S0")
    S, _ = F.smooth_all(T, lambda c, h: T[(c, h)])
    hE2 = ["F10_S%dE2" % j for j in range(1, 6)]
    B = F.make_bins(X, T, S, hE2, 0.575, 0.700)
    def model_A(eta):
        a, ka, rho, g = eta
        def m(th, Bx):
            c = th[0:5][Bx.j]; b = th[5]; kB = th[6]; k = th[7]
            co = {"C2": a, "C3": a * ka, "C4": a * ka, "C5": a * rho, "Geta": a * g}
            mu = k * Bx.cb + c * (1. + b * Bx.x) * (Bx.gd + kB * Bx.gc) + sum(co[t] * Bx.s[t] for t in F.SIG)
            vs = sum(co[t] ** 2 * Bx.vs[t] for t in F.SIG)
            return mu, vs, k
        return m
    cons = [(7, 1.0, 0.05)]
    names = ["c_1", "c_2", "c_3", "c_4", "c_5", "b", "kappa_B", "k_E2"]
    wp3 = [0.98777581, 1.0012206, 0.97505082, 1.0551912, 0.90912575, -0.068343866, 3.0627053, 0.95434004]
    res = {}
    for tag, eta in (("a_eta(0) = 0.4, g = 1 (WP3 frozen)", (0.4, 1, 1, 1)),
                     ("WP5 F-main: a 0.3823, kappa 0.8323, rho 1.3152, g 1.7168", (0.38228, 0.83234, 1.31520, 1.71677)),
                     ("WP5 a_eta, g = 1", (0.38228, 0.83234, 1.31520, 1.0))):
        th, _ = F.fit_irls(B, model_A(eta), cons, np.array([1, 1, 1, 1, 1, -0.07, 3.0, 1.0]))
        E = F.errors(B, model_A(eta), cons, th)
        res[tag] = (th, np.sqrt(np.diag(E["Cstat"])))
        P("%-60s %s  Pearson %.1f/%d" % (tag, " ".join("%s %.4f" % (n, v) for n, v in zip(names, th)), E["chi2P"], B.n - 8 + 1))
    P("WP3 frozen values: " + " ".join("%s %.4f" % (n, v) for n, v in zip(names, wp3)))
    t0, e0 = res["a_eta(0) = 0.4, g = 1 (WP3 frozen)"]
    t1, _ = res["WP5 F-main: a 0.3823, kappa 0.8323, rho 1.3152, g 1.7168"]
    P("shift (WP5 values - frozen), in units of the stat error: " + " ".join("%s %+.4f (%+.2f)" % (n, d, d / e) for n, d, e in zip(names, t1 - t0, e0)))
    # background prediction in the bump window B and the peak window P (E2, all slices) for both
    Bb = F.make_bins(X, T, S, hE2)
    inB = (F.edges[Bb.i] >= 0.455 - 1e-9) & (F.edges[Bb.i + 1] <= 0.5225 + 1e-9)
    for tag in res:
        th = res[tag][0]
        c = th[0:5][Bb.j]; bg = c * (1 + th[5] * Bb.x) * (Bb.gd + th[6] * Bb.gc)
        P("stage-A background in B (E2, F10 cells) with %-55s %.1f rows" % (tag, bg[inB].sum()))

if __name__ == "__main__":
    main()
