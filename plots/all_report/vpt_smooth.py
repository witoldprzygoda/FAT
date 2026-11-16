#!/usr/bin/env python3
"""vpt_smooth.py -- verifier of peaktail: tail scale t and its propagation with SMOOTHED C4 templates (each of the few
C4 rows replaced by 400 entries with the C5 residual of its M(pipi) bin), and with the peak scale a_P treated as a
measured quantity that itself depends on t (a_P(t) from the side bands, error 0.080).
Output: stdout (tee vpt_smooth.log), vpt_smooth.json.
"""
import json, math, sys
import numpy as np
sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import vpt_core as V
import vpt_tmpl as TM

g = V.Grid()
Z = np.load(V.OUT + "vpt_nom.npz")
N, CB, COV, T = Z["N"], Z["CB"], Z["COV"], Z["T"]
sig = N - CB
BUMP = Z["C2"] + Z["C3"]
R = TM.rows(); S = R["S"]
sm = TM.smooth_c4()
out = {}
AP_ERR = 0.080


# ------------------------------------------------------------------ a_P(t) / a_P(1) from the side bands (smoothed C4)
def ap_ratio_table(flag):
    x5 = R["m4p"]; m5 = R["m45"] & (S["g"] == 3)
    P5 = (x5 >= 0.53364) & (x5 < 0.55416); SB5 = ((x5 >= 0.5250) & (x5 < 0.5336)) | ((x5 >= 0.5542) & (x5 < 0.5628))
    x4 = sm["m4"]; P4 = (x4 >= 0.53364) & (x4 < 0.55416); SB4 = ((x4 >= 0.5250) & (x4 < 0.5336)) | ((x4 >= 0.5542) & (x4 < 0.5628))
    def q(t):
        Pc = SBc = 0.0
        for w, m, f, Pm, SBm in ((R["w"], m5, R["flag"][flag], P5, SB5), (sm["w"], np.ones(len(x4), bool), sm[flag], P4, SB4)):
            Nn = w[m].sum(); Nt = w[m & f].sum(); Nc = Nn - Nt
            wt = np.where(f, w * t, w * (Nn - t * Nt) / Nc)
            Pc += wt[m & Pm].sum(); SBc += wt[m & SBm].sum()
        return (1 - 1.1859 * SBc / Pc) * Pc, SBc / Pc
    q1, f1 = q(1.0)
    ts = np.linspace(-1.5, 4.0, 56)
    return ts, np.array([q1 / q(t)[0] for t in ts]), f1


def fit_c(sig_, K_t, A0, sA, lo=-0.045, hi=0.045, order=2, nit=7):
    """pol + A K + h T with a Gaussian constraint A = A0 +- sA (sA = 0: A fixed at A0; sA = inf: free)"""
    i = np.where((g.x >= lo - 1e-12) & (g.x < hi - 1e-12))[0]
    t = (g.x[i] - 0.5 * (lo + hi)) / (0.5 * (hi - lo))
    cb = CB[i]; V0 = COV[np.ix_(i, i)]
    if sA == 0:
        F = np.vstack([t ** j for j in range(order + 1)] + [T[i]]).T; y = sig_[i] - A0 * K_t[i]; fx = A0 * K_t[i]
        Vm = V0 + np.diag(np.maximum(sig_[i] + cb, 1.0))
        for _ in range(nit):
            p, C, L, chi2 = V.gls(y, F, Vm); Vm = V0 + np.diag(np.maximum(F @ p + fx + cb, 1.0))
        p, C, L, chi2 = V.gls(y, F, Vm)
        return chi2, p[-1], A0, L[-1]
    F = np.vstack([t ** j for j in range(order + 1)] + [K_t[i], T[i]]).T
    y = sig_[i].copy()
    Vm = V0 + np.diag(np.maximum(sig_[i] + cb, 1.0))
    for _ in range(nit):
        if np.isfinite(sA):
            Fa = np.vstack([F, np.eye(F.shape[1])[order + 1]]); ya = np.append(y, A0)
            Va = np.zeros((len(i) + 1, len(i) + 1)); Va[:-1, :-1] = Vm; Va[-1, -1] = sA * sA
        else:
            Fa, ya, Va = F, y, Vm
        p, C, L, chi2 = V.gls(ya, Fa, Va)
        Vm = V0 + np.diag(np.maximum(F @ p + cb, 1.0))
    return chi2, p[-1], p[order + 1], L[-1][:len(i)]


def prof(vals, ts):
    vals = np.asarray(vals); j = int(np.argmin(vals))
    jj = np.clip(np.arange(j - 3, j + 4), 0, len(ts) - 1)
    a, b, c = np.polyfit(ts[jj], vals[jj], 2)
    tmin = -b / (2 * a) if a > 0 else ts[j]
    ins = ts[vals <= vals[j] + 1]
    return tmin, ins.min(), ins.max(), vals[j]


ts = np.linspace(-1.5, 4.0, 56)
print("=" * 110)
print("A. u fit (0.30-0.40), smoothed C4. a_P treatments: fixed 0.418 | a_P(t) (side bands redone with K(t)) fixed | a_P(t) +- 0.080 | free")
Wb = TM.group_boot(100, seed=4242)
for flag in ("sub5", "p95", "sub10", "res15"):
    Tt = TM.templates_smooth("u", g.e, 0.30, 0.40, flag)
    tsr, apr, f1 = ap_ratio_table(flag)
    K1 = TM.K_of_t(Tt, 1.0)
    c0 = fit_c(sig, K1, V.AP, 0)
    f_nom_s = c0[1] / float(c0[3] @ BUMP[np.where((g.x >= -0.045 - 1e-12) & (g.x < 0.045 - 1e-12))[0]])
    if flag == "sub5":
        print("  smoothed-C4 nominal: h %.2f f %.5f (unsmoothed 0.42004); side-band leakage at t = 1: %.4f (unsmoothed 0.1544)" % (c0[1], f_nom_s, f1))
        out["f_nom_smooth"] = f_nom_s
    for lo, hi, od in ((-0.045, 0.045, 2), (-0.045, 0.060, 2), (-0.045, 0.070, 3)):
        res = {}
        for name, mk in (("fixed", lambda t, ap: (V.AP, 0)), ("ap(t)", lambda t, ap: (V.AP * ap, 0)),
                         ("ap(t)+-0.08", lambda t, ap: (V.AP * ap, AP_ERR * ap)), ("free", lambda t, ap: (V.AP, np.inf))):
            ch = []
            for t, ap in zip(ts, apr):
                A0, sA = mk(t, ap)
                ch.append(fit_c(sig, TM.K_of_t(Tt, t), A0, sA, lo, hi, od)[0])
            res[name] = prof(ch, ts)
        # template MC of the fixed-a_P t (bootstrap of C5 rows and C4 parent rows)
        mc = []
        if lo == -0.045 and hi == 0.045 or flag == "sub5":
            for r in range(60):
                Tr = TM.templates_smooth("u", g.e, 0.30, 0.40, flag, wC5=R["w"] * Wb[r], wC4row=Wb[r])
                mc.append(prof([fit_c(sig, TM.K_of_t(Tr, t), V.AP, 0, lo, hi, od)[0] for t in ts], ts)[0])
        mcs = float(np.std(mc, ddof=1)) if mc else float("nan")
        print("  %-5s %+3.0f..%+3.0f pol%d: " % (flag, lo * 1e3, hi * 1e3, od) + " | ".join(
            "%s t=%.2f [%.2f,%.2f]" % (k, v[0], v[1], v[2]) for k, v in res.items()) + " | MC(fixed) %.2f" % mcs)
        out["u_%s_%g_%g_%d" % (flag, lo, hi, od)] = dict({k: list(map(float, v)) for k, v in res.items()}, mc=mcs)
    # f(t) slope with smoothed templates, a_P fixed and a_P(t)
    fs = []
    i45 = np.where((g.x >= -0.045 - 1e-12) & (g.x < 0.045 - 1e-12))[0]
    for t in (0.0, 0.5, 1.0, 1.5, 2.0):
        ap = float(np.interp(t, tsr, apr))
        a = fit_c(sig, TM.K_of_t(Tt, t), V.AP, 0); b = fit_c(sig, TM.K_of_t(Tt, t), V.AP * ap, 0)
        fs.append((t, a[1] / float(a[3] @ BUMP[i45]), b[1] / float(b[3] @ BUMP[i45])))
    f1v = [x for x in fs if x[0] == 1.0][0][1]
    sl1 = np.polyfit([x[0] for x in fs], [x[1] / f1v for x in fs], 1)[0]; sl2 = np.polyfit([x[0] for x in fs], [x[2] / f1v for x in fs], 1)[0]
    print("        f(t)/f(1) a_P fixed: %s  slope %+.4f | a_P(t): %s slope %+.4f" % (
        " ".join("%.4f" % (x[1] / f1v) for x in fs), sl1, " ".join("%.4f" % (x[2] / f1v) for x in fs), sl2))
    out["fslope_%s" % flag] = [sl1, sl2]

print("\nB. line shape where no bump exists (M4, amplitude free), smoothed C4; position shift dm and width r of the non-loss part")
gm = V.Grid(lo=0.30, bw=0.0025, nb=220, reb=4)
D = V.load_data(); S0d = (D["bits"] & 1) != 0
m4d = D["m4"].astype(float)
from vpt_nom import in_sl


def dspec(slices):
    Nn = np.zeros(gm.nb); C_ = np.zeros(gm.nb); CV = np.zeros((gm.nb, gm.nb))
    for a, b in slices:
        m = S0d & in_sl(D["mpipi"], a, b)
        n = gm.h(m4d[m & (D["k"] == 0)]); p = gm.h(m4d[m & (D["k"] == 1)]); q = gm.h(m4d[m & (D["k"] == 2)])
        c, cv = V.cb_cells(p, q, gm); Nn += n; C_ += c; CV += cv
    return Nn, C_, CV


def fitK(sg, cb_, cov_, K, lo, hi, od):
    i = np.where((gm.x >= lo - 1e-12) & (gm.x < hi - 1e-12))[0]
    t = (gm.x[i] - 0.5 * (lo + hi)) / (0.5 * (hi - lo))
    F = np.vstack([t ** j for j in range(od + 1)] + [K[i]]).T; y = sg[i]; V0 = cov_[np.ix_(i, i)]
    Vm = V0 + np.diag(np.maximum(sg[i] + cb_[i], 1.0))
    for _ in range(7):
        p, C, L, chi2 = V.gls(y, F, Vm); Vm = V0 + np.diag(np.maximum(F @ p + cb_[i], 1.0))
    p, C, L, chi2 = V.gls(y, F, Vm)
    return p, C, chi2, len(i) - F.shape[1]


base = V.S_NOM * S["m4t"].astype(float)
Wb2 = TM.group_boot(60, seed=777)
for rn, sl in (("0.42-0.46", [(0.42, 0.44), (0.44, 0.46)]), ("0.415-0.46", [(0.415, 0.42), (0.42, 0.44), (0.44, 0.46)])):
    a0, b0 = sl[0][0], sl[-1][1]
    Nn, C_, CV = dspec(sl); sg = Nn - C_
    for flag in ("sub5", "p95", "sub10"):
        for lo, hi, od in ((0.49, 0.60, 2), (0.48, 0.62, 3), (0.50, 0.59, 2), (0.47, 0.62, 3)):
            tab = {}
            for r in (0.90, 1.00):
                for dm in (-2.0, -1.0, 0.0, 1.0, 2.0):
                    x5 = np.where(S["m4t"] > 0, base + r * (R["d"] - R["dlep"]) + R["dlep"], R["m4p"]) + dm * 1e-3
                    smx = dict(sm); smx["m4"] = sm["m4"] + dm * 1e-3      # smoothed C4 keeps the pool residual (r = 1)
                    Rsave = R["m4p"]; R["m4p"] = x5
                    Tt = TM.templates_smooth("m4", gm.e, a0, b0, flag, sm=smx)
                    R["m4p"] = Rsave
                    ch = [fitK(sg, C_, CV, TM.K_of_t(Tt, t), lo, hi, od)[2] for t in ts]
                    tab[(r, dm)] = prof(ch, ts)
            c0 = tab[(1.0, 0.0)]
            bestk = min(tab, key=lambda k: tab[k][3])
            mc = float("nan")
            if flag == "sub5" and od == 2 and lo == 0.49:
                Tt0 = TM.templates_smooth("m4", gm.e, a0, b0, flag)
                vv = []
                for rr_ in range(60):
                    Tr = TM.templates_smooth("m4", gm.e, a0, b0, flag, wC5=R["w"] * Wb2[rr_], wC4row=Wb2[rr_])
                    vv.append(prof([fitK(sg, C_, CV, TM.K_of_t(Tr, t), lo, hi, od)[2] for t in ts], ts)[0])
                mc = float(np.std(vv, ddof=1))
            print("  %-10s %-5s %.3f-%.3f pol%d: nominal t = %.2f [%.2f, %.2f] chi2 %.1f | t at r=1, dm=-2..+2: %s (chi2 %s) | best r %.2f dm %+.0f: t %.2f chi2 %.1f | MC %.2f" % (
                rn, flag, lo, hi, od, c0[0], c0[1], c0[2], c0[3],
                " ".join("%.2f" % tab[(1.0, d)][0] for d in (-2.0, -1.0, 0.0, 1.0, 2.0)),
                " ".join("%.1f" % tab[(1.0, d)][3] for d in (-2.0, -1.0, 0.0, 1.0, 2.0)),
                bestk[0], bestk[1], tab[bestk][0], tab[bestk][3], mc))
            out["line_%s_%s_%g_%g_%d" % (rn, flag, lo, hi, od)] = dict(nom=list(map(float, c0)), best=[bestk[0], bestk[1]] + list(map(float, tab[bestk])), mc=mc,
                                                                       dm=[float(tab[(1.0, d)][0]) for d in (-2.0, -1.0, 0.0, 1.0, 2.0)])
json.dump(out, open(V.OUT + "vpt_smooth.json", "w"), indent=1)
