#!/usr/bin/env python3
"""vpt_line.py -- verifier of peaktail: tail scale t from the eta peak line shape where the bump cannot contribute
(M(pi+pi-) >= 0.42 GeV/c2), own code. N - CB in M4 (CB per 20 MeV M(pipi) slice, 10 MeV cells from 0.30),
model pol(n) + A K(t), A free, chi2 profile in t. Template MC from simulated event-group bootstraps.
Output: stdout (tee vpt_line.log), vpt_line.json.
"""
import json
import math
import sys
import numpy as np
sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import vpt_core as V
import vpt_tmpl as TM
from vpt_nom import in_sl

gm = V.Grid(lo=0.30, bw=0.0025, nb=220, reb=4)
ZERO = np.zeros(gm.nb)
D = V.load_data(); S0 = (D["bits"] & 1) != 0
R = TM.rows(); Ssim = R["S"]
m4d = D["m4"].astype(float)


def dspec(slices, extra=None):
    N = np.zeros(gm.nb); CB = np.zeros(gm.nb); COV = np.zeros((gm.nb, gm.nb))
    for a, b in slices:
        m = S0 & in_sl(D["mpipi"], a, b)
        if extra is not None:
            m &= extra
        n = gm.h(m4d[m & (D["k"] == 0)]); p = gm.h(m4d[m & (D["k"] == 1)]); q = gm.h(m4d[m & (D["k"] == 2)])
        c, cv = V.cb_cells(p, q, gm); N += n; CB += c; COV += cv
    return N, CB, COV


def fitK(sig, CB, COV, K, lo, hi, od, fixedA=None):
    i = np.where((gm.x >= lo - 1e-12) & (gm.x < hi - 1e-12))[0]
    t = (gm.x[i] - 0.5 * (lo + hi)) / (0.5 * (hi - lo))
    cols = [t ** j for j in range(od + 1)]
    fx = np.zeros(len(i))
    if fixedA is None:
        cols.append(K[i])
    else:
        fx = fixedA * K[i]
    F = np.vstack(cols).T; y = sig[i] - fx; V0 = COV[np.ix_(i, i)]
    Vm = V0 + np.diag(np.maximum(sig[i] + CB[i], 1.0))
    for _ in range(7):
        p, C, L, chi2 = V.gls(y, F, Vm)
        Vm = V0 + np.diag(np.maximum(F @ p + fx + CB[i], 1.0))
    p, C, L, chi2 = V.gls(y, F, Vm)
    return p, C, chi2, len(i) - F.shape[1]


def profile(sig, CB, COV, T, lo, hi, od, tg=np.linspace(-1.5, 4.0, 111)):
    ch = np.array([fitK(sig, CB, COV, TM.K_of_t(T, t), lo, hi, od)[2] for t in tg])
    j = int(np.argmin(ch))
    # parabola through the minimum
    jj = np.clip(np.arange(j - 3, j + 4), 0, len(tg) - 1)
    a, b, c = np.polyfit(tg[jj], ch[jj], 2)
    tmin = -b / (2 * a) if a > 0 else tg[j]
    inside = tg[ch <= ch[j] + 1.0]
    return tmin, inside.min(), inside.max(), ch[j], ch, tg


out = {}
REG = {"0.42-0.46": [(0.42, 0.44), (0.44, 0.46)], "0.415-0.46": [(0.415, 0.42), (0.42, 0.44), (0.44, 0.46)],
       "0.42-0.50": [(0.42, 0.44), (0.44, 0.46), (0.46, 0.48), (0.48, 0.50)]}
FITS = [(0.49, 0.60, 2), (0.48, 0.62, 3), (0.50, 0.59, 2), (0.47, 0.62, 3), (0.49, 0.61, 2), (0.485, 0.605, 3)]
NB = 100
Wb = TM.group_boot(NB, seed=99)
for rn, sl in REG.items():
    a0, b0 = sl[0][0], sl[-1][1]
    N, CB, COV = dspec(sl)
    sig = N - CB
    regs = TM.in_pp(R["mpp"], a0, b0)
    for fl in ("sub5", "p95", "sub10"):
        T = TM.templates(R["m4p"], gm.e, regs, R["flag"][fl])
        print("%s  flag %s: C4 rows N_eff %.0f, C5 N_eff %.0f; tail share C4 %.3f C5 %.3f; peak in data (sum A K) ~ %.0f rows; N in 0.53-0.555: %.0f" % (
            rn, fl, 0, 0, T["C4"][5] / T["C4"][3], T["C5"][5] / T["C5"][3], 0, N[(gm.x > 0.53) & (gm.x < 0.555)].sum()))
        for lo, hi, od in FITS:
            tmin, tl, th, c2, ch, tg = profile(sig, CB, COV, T, lo, hi, od)
            p, C, c2b, ndf = fitK(sig, CB, COV, TM.K_of_t(T, tmin), lo, hi, od)
            A = p[-1]
            scale = math.sqrt(max(c2 / ndf, 1.0))
            dch = float(ch[np.argmin(np.abs(tg - 0.0))] - ch[np.argmin(np.abs(tg - 1.0))])
            # MC: refit the minimum with bootstrap templates (only for the first two fits, sub5)
            mc = float("nan")
            if fl == "sub5" and (lo, hi, od) in FITS[:2] and rn == "0.42-0.46":
                tm = []
                for r in range(NB):
                    Tr = TM.templates(R["m4p"], gm.e, regs, R["flag"][fl], w=R["w"] * Wb[r])
                    tm.append(profile(sig, CB, COV, Tr, lo, hi, od, tg=np.linspace(-1.5, 4.0, 56))[0])
                mc = float(np.std(tm, ddof=1))
            out["%s_%s_%g_%g_%d" % (rn, fl, lo, hi, od)] = dict(t=tmin, lo=tl, hi=th, chi2=c2, ndf=ndf, A=A, scale=scale, dchi2_0_1=dch, mc=mc)
            print("   %.3f-%.3f pol%d: t = %5.2f  [%5.2f, %5.2f] (Delta chi2 = 1)  chi2 %.1f/%d  -> sqrt(chi2/ndf) %.2f  A = %.3f  "
                  "chi2(t=0) - chi2(t=1) = %+.1f  MC %.2f" % (lo, hi, od, tmin, tl, th, c2, ndf, scale, A, dch, mc))
    # resolution variants (sub5): r = 0.893 on all rows, and on the core rows only
    T1 = TM.templates(R["m4p"], gm.e, regs, R["flag"]["sub5"])
    m4r = V.sim_m4p(Ssim, r=0.893)
    m4c = np.where(R["flag"]["sub5"], R["m4p"], m4r)
    for name, xx in (("r = 0.893 all rows", m4r), ("r = 0.893 core rows", m4c)):
        T = TM.templates(xx, gm.e, regs, R["flag"]["sub5"])
        for lo, hi, od in FITS[:2]:
            tmin, tl, th, c2, ch, tg = profile(sig, CB, COV, T, lo, hi, od)
            print("   %s %.3f-%.3f pol%d: t = %5.2f [%5.2f, %5.2f] chi2 %.1f" % (name, lo, hi, od, tmin, tl, th, c2))
            out["%s_%s_%g_%g_%d" % (rn, name, lo, hi, od)] = dict(t=tmin, lo=tl, hi=th, chi2=c2)
    # M(ee) < m_pi0 only
    if rn == "0.42-0.46":
        N2, CB2, COV2 = dspec(sl, extra=D["mee"] < V.MPI0)
        T2 = TM.templates(R["m4p"], gm.e, regs & (Ssim["ee"] < V.MPI0), R["flag"]["sub5"])
        for lo, hi, od in FITS[:2]:
            tmin, tl, th, c2, ch, tg = profile(N2 - CB2, CB2, COV2, T2, lo, hi, od)
            print("   M(ee) < m_pi0 %.3f-%.3f pol%d: t = %5.2f [%5.2f, %5.2f] chi2 %.1f" % (lo, hi, od, tmin, tl, th, c2))
json.dump(out, open(V.OUT + "vpt_line.json", "w"), indent=1)
