#!/usr/bin/env python3
# ============================================================================
# x17p26_fold.py  --  X17 limit with PDG 2026 (analyst ruling D25, 15 November 2026), prefix x17p26, written 2026-11-15:
#   copy of fix2_x17_fold.py with K = BR(eta -> pi+pi-e+e-) = 2.76e-4 (PDG 2026) instead of the frozen 2.68e-4, and the factor k of the
#   PDG 2026 BR_dir uncertainty in the X17 yield nuisance (rerun of the frozen code, x17p26_limit.C: k = 1.00124 at 17.00, 1.00132 at
#   16.85 MeV/c2) applied to the result (the vx profile scans were made with the frozen s_X). The folding is linear in K.
#   Writes only /hdd3/.../all_report/x17p26/x17p26_fold.log and docs/latex/tables_all/x17_fold_pdg26.txt.
# Run: python3 -I plots/all_report/x17p26_fold.py   (seconds)
# --- header of fix2_x17_fold.py follows ---
# ============================================================================
# fix2_x17_fold.py  --  final review (prefix fix2): copy of fix_x17_fold.py with the total uncertainty of sigma
#   updated to 38.3 % (the +-5 % period normalisation of Part II added to N; fix2_numbers.txt). Original header follows.
# fix_x17_fold.py  --  review fix (prefix fix): limit on sigma(pp -> eta X) x BR(eta -> pi+pi- X17) x BR(X17 -> e+e-)
#   with the uncertainty of sigma(pp -> eta X) folded in as a log-normal nuisance (profiled), instead of the plain product
#   UL(BR) x sigma_eta.
#   Input: the profile-likelihood scan of the verifier vx (vx_prof_m1700_n1.txt, vx_prof_m1685_n1.txt): columns r, BR,
#          F_obs - F_min, F_Asimov - F_A0 (F = -2 ln L); q~ = (F - F_min)/qScale, qScale = 1.0299 (frozen), q~ = 0 for r < r_hat.
#   Check: the frozen asymptotic CLs limit (BR 3.0636e-5 at 17.0 MeV) is reproduced from the scan.
#   Folding: P = r x K x sigma, K = BR(eta -> pi+pi-e+e-) = 2.68e-4 (frozen, PDG 2020), sigma = sigma0 exp(d theta),
#            theta ~ N(0, 1);  q~_P(P) = min_theta [ q~(P / (K sigma0 e^{d theta})) + theta^2 ] (minus its global minimum),
#            q_A,P(P) likewise from the Asimov column; CLs = (1 - Phi(sqrt q~_P)) / Phi(sqrt q_A,P - sqrt q~_P) = 0.10.
#   d = sqrt(ln(1 + s^2)) for a relative standard deviation s of sigma0.
# Writes only /hdd3/.../all_report/fix/fix2_x17_fold.log and docs/latex/tables_all/fix2_x17_fold.txt.
# Run: python3 -I plots/all_report/fix2_x17_fold.py   (seconds)
# ============================================================================
import math
import numpy as np
from scipy.stats import norm
from scipy.optimize import minimize_scalar, brentq

VX = "/hdd3/przygoda/HADES/FAT_out/all_report/vx/"
OUT = "/hdd3/przygoda/HADES/FAT_out/all_report/x17p26/"
TAB = "/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all/x17_fold_pdg26.txt"
QS = 1.0299
K = 2.76e-4
KSX = {"17.00": 1.00124, "16.85": 1.00132}   # rerun factor of the PDG 2026 BR_dir uncertainty in s_X
LOG = open(OUT + "x17p26_fold.log", "w")
TX = open(TAB, "w")
def P(*a, tab=False):
    s = " ".join(str(x) for x in a); print(s); LOG.write(s + "\n"); LOG.flush()
    if tab: TX.write(s + "\n"); TX.flush()

def load(fn):
    d = np.array([[float(x) for x in l.split()] for l in open(fn) if l.strip() and not l.startswith(("vx_", "#"))])
    r, qo, qa = d[:, 0], d[:, 2] / QS, d[:, 3] / QS
    # unconditional minimum of the observed curve (r_hat >= 0 for q~): q~ = 0 below r_hat
    i0 = np.argmin(qo); rhat = r[i0]
    qt = np.where(r < rhat, 0.0, qo - qo[i0])
    qa = qa - qa[np.argmin(np.abs(r))]
    return r, qt, qa, rhat

def interp(r, y):
    def f(x):
        if x <= r[-1]: return float(np.interp(x, r, y))
        # quadratic extrapolation from the last 10 points
        c = np.polyfit(r[-10:], y[-10:], 2); return float(np.polyval(c, x))
    return f

def cls(qt, qa):
    # the formula used by the frozen limit and by vx's by-hand check (valid for q~ <= qA, which holds near the limit)
    qt = max(qt, 0.0); qa = max(qa, 1e-12)
    return (1 - norm.cdf(math.sqrt(qt))) / max(norm.cdf(math.sqrt(qa) - math.sqrt(qt)), 1e-300)

def ul_plain(r, qt, qa):
    fo, fa = interp(r, qt), interp(r, qa)
    g = lambda x: cls(fo(x), fa(x)) - 0.10
    return brentq(g, 0.03, 0.45)

def ul_fold(r, qt, qa, sig0, s_rel):
    d = math.sqrt(math.log(1 + s_rel ** 2))
    fo, fa = interp(r, qt), interp(r, qa)
    def prof(f, Pv):
        h = lambda th: f(Pv / (K * sig0 * math.exp(d * th))) + th * th
        res = minimize_scalar(h, bounds=(-4, 4), method="bounded", options={"xatol": 1e-7})
        return res.fun
    # global minimum of the observed folded curve is 0 (at r_hat, theta = 0); Asimov minimum 0 at P = 0
    g = lambda Pv: cls(prof(fo, Pv), prof(fa, Pv)) - 0.10
    lo, hi = 0.03 * K * sig0, 0.42 * K * sig0
    return brentq(g, lo, hi)

def main():
    P("x17p26_fold: 90 % CLs limits on sigma(pp -> eta X) x BR_X x BR_ee [nb], normalisation (1), K = 2.76e-4 (PDG 2026; ruling D25)", tab=True)
    P("columns: mass, UL(BR) from the scan, plain product UL x sigma0, folded with s_rel of sigma0; [x k]: with the PDG 2026 BR_dir uncertainty in s_X", tab=True)
    sig0 = 0.429e6   # nb
    srels = [("total before the final review 37.9 %", 0.379), ("adopted total with the period normalisation of Part II 38.3 %", 0.3826)]
    CV, CV95 = 0.99970 / 0.994, 0.99970 / (0.994 - 1.645 * 0.021)   # C_vtx and its 95 % CL bound (not applied)
    for m, fn in (("17.00", "vx_prof_m1700_n1.txt"), ("16.85", "vx_prof_m1685_n1.txt")):
        r, qt, qa, rhat = load(VX + fn)
        u = ul_plain(r, qt, qa)
        k = KSX[m]
        P("m_X %s MeV: r_hat %.4f; UL r %.5f -> BR %.4e [x k: %.4e] (PDG 2020 frozen: %.4e); plain product with sigma0 = 0.429 mb: %.2f nb [x k: %.2f nb]"
          % (m, rhat, u, u * K, u * K * k, u * 2.68e-4, u * K * sig0, u * K * sig0 * k), tab=True)
        for nm, s in srels:
            Pf = ul_fold(r, qt, qa, sig0, s)
            P("   folded (%s): UL = %.2f nb [x k: %.2f nb]  (x %.3f of the plain product)" % (nm, Pf, Pf * k, Pf / (u * K * sig0)), tab=True)
        for nm, fac in (("C_vtx applied", CV), ("C_vtx at its 95 % CL bound", CV95)):
            Pf = ul_fold(r, qt, qa, sig0 * fac, 0.3826)
            P("   with sigma_eta x %.4f (%s), folded 38.3 %%: UL = %.2f nb [x k: %.2f nb]" % (fac, nm, Pf, Pf * k), tab=True)

if __name__ == "__main__":
    main()
