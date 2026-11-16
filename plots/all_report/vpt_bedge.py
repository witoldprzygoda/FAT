#!/usr/bin/env python3
"""vpt_bedge.py -- verifier of peaktail: f_eta(edge) of the per-slice method (bedge) versus the tail scale t.
Own arrays from vpt_dump (data S0; SMASH C2, C3, T12, C4, C5 in the vbedge layout) and the independent bedge
re-implementation vbedge_core (fit, Cfg, data_spec; evaluate/quant re-coded here from vbedge_ana.py). The tail enters
through T12 = ... + C4(t) + C5(t). Output: stdout (tee vpt_bedge.log), vpt_bedge.json.
"""
import json, sys
import numpy as np
sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import vpt_core as V
import vpt_tmpl as TM
import vbedge_core as B

NM4, NPP, NCL = 240, 40, 3
D = V.load_data(); R = TM.rows(); S = R["S"]


def ppbin(m):
    m = m.astype(np.float64)
    b = np.floor((m - 0.27) / 0.005 + 1e-12).astype(int)
    b[(m < 0.27) | (b >= NPP)] = -1
    return b


def m4bin(m):
    m = m.astype(np.float64)
    b = np.floor((m - 0.35) / 0.0025).astype(int)
    b[(m < 0.35) | (b >= NM4)] = -1
    return b


def clsof(mee):
    c = np.where(mee < 0.010, 0, 1)
    c[~np.isfinite(mee)] = 2
    return c


def fill(m4, pp, mee, w, sel):
    a = np.zeros((NCL, NPP, NM4))
    ib, jb, cb = m4bin(m4), ppbin(pp), clsof(mee)
    ok = sel & (ib >= 0) & (jb >= 0)
    np.add.at(a, (cb[ok], jb[ok], ib[ok]), w[ok])
    return a


S0 = (D["bits"] & 1) != 0
Nd = np.stack([fill(D["m4"], D["mpipi"], D["mee"].astype(float), np.ones(len(D["k"])), S0 & (D["k"] == k)) for k in range(3)])
ws = V.sim_weight(S, c2=False)
mppS = (V.SPP * S["pp"].astype(np.float64))
eeS = S["ee"].astype(float)
SA = {}
for name, gs in (("C2", [0]), ("C3", [1]), ("T12", list(range(9))), ("C4", [2]), ("C5", [3])):
    SA[name] = fill(R["m4p"], mppS, eeS, ws, np.isin(S["g"], gs) & (ws > 0))
NOMSET = B.NOMSET
CF = {"toy": (0.880, 0.889, 0.868), "clos": (0.891, 0.869, 0.896)}
CF["rec"] = tuple(0.5 * (a + b) for a, b in zip(CF["toy"], CF["clos"]))


def simy(arr, cls, pa, pbb):
    return arr[:, pa:pbb].sum(axis=(0, 1)) if cls < 0 else arr[cls, pa:pbb].sum(0)


def evaluate(T12):
    res = []
    for k in range(3):
        cls = k - 1; Ed = Es = 0.0; lin = np.zeros(3)
        for lo, hi in NOMSET:
            pa, pbb = B.pb(lo), B.pb(hi)
            y, cb, cov, u = B.data_spec(Nd, cls, pa, pbb)
            cD = B.Cfg(lo, hi)
            h, eh, c2, nd, Vd, p = B.fit(y, cb, cov, cD)
            Es += B.fit(simy(T12, cls, pa, pbb), None, None, cD)[0]
            for g_, arr in enumerate((T12, SA["C2"], SA["C3"])):
                lin[g_] += B.fit(simy(arr, cls, pa, pbb), None, None, cD, Vfix=Vd)[0]
            Ed += h
        res.append((Ed, Es, lin))
    return res


def quant(res):
    Ed = [r[0] for r in res]; Es = [r[1] for r in res]; lin = [r[2] for r in res]
    t = lin[1][0] + lin[2][0]; c = CF["rec"]
    red = ((1 - c[1]) * lin[1][1] + (1 - c[2]) * lin[2][1]) / t
    return dict(f_all=Ed[0] / Es[0], fc_rec=Ed[0] / (Es[0] * (1 - red)), Es=Es[0], Ed=Ed[0], rest=lin[0][0] - lin[0][1] - lin[0][2])


q0 = quant(evaluate(SA["T12"]))
print("nominal (own arrays): E_data %.4f  E_sim %.4f  f_all %.5f  fc_rec %.5f  rest of E_sim (lin) %.1f   (bedge: 1141.8584, 2915.6983, 0.39162, 0.41841, 31.9)" % (
    q0["Ed"], q0["Es"], q0["f_all"], q0["fc_rec"], q0["rest"]))
out = {"nom": q0}


def dK_arrays(flag, local=True, smooth=False):
    """d T12 / dt in the vbedge layout: per class C4/C5 (and per 5 MeV M(pipi) bin when local) K_tail - (N_t/N_c) K_core"""
    out_ = np.zeros((NCL, NPP, NM4))
    if not smooth:
        for gi in (2, 3):
            m = R["m45"] & (S["g"] == gi)
            jb = ppbin(mppS)
            keys = np.unique(jb[m]) if local else [None]
            for b_ in keys:
                mm = m if b_ is None else m & (jb == b_)
                Nn = R["w"][mm].sum(); Nt = R["w"][mm & flag].sum(); Nc = Nn - Nt
                if Nc <= 0 or Nt <= 0:
                    continue
                wt = np.where(flag, R["w"], -R["w"] * Nt / Nc)
                out_ += fill(R["m4p"], mppS, eeS, wt, mm)
    else:
        sm = TM.smooth_c4()
        # C5 as above, C4 from the smoothed entries
        m = R["m45"] & (S["g"] == 3); jb = ppbin(mppS)
        for b_ in np.unique(jb[m]):
            mm = m & (jb == b_); Nn = R["w"][mm].sum(); Nt = R["w"][mm & flag].sum(); Nc = Nn - Nt
            if Nc <= 0 or Nt <= 0:
                continue
            out_ += fill(R["m4p"], mppS, eeS, np.where(flag, R["w"], -R["w"] * Nt / Nc), mm)
        f4 = sm[FLAGNAME]; jb4 = ppbin(sm["mpp"]); ee4 = eeS[sm["i"]]
        for b_ in np.unique(jb4):
            mm = jb4 == b_; Nn = sm["w"][mm].sum(); Nt = sm["w"][mm & f4].sum(); Nc = Nn - Nt
            if Nc <= 0 or Nt <= 0:
                continue
            out_ += fill(sm["m4"], sm["mpp"], ee4, np.where(f4, sm["w"], -sm["w"] * Nt / Nc), mm)
    return out_


for FLAGNAME in ("sub5", "p95", "sub10"):
    for local, smooth in ((True, False), (False, False), (True, True)):
        dk = dK_arrays(R["flag"][FLAGNAME], local, smooth)
        base = SA["T12"].copy()
        if smooth:   # replace the raw C4 rows by the smoothed C4 entries at t = 1
            sm = TM.smooth_c4()
            base = base - SA["C4"] + fill(sm["m4"], sm["mpp"], eeS[sm["i"]], sm["w"], np.ones(len(sm["w"]), bool))
        qs = {t: quant(evaluate(base + (t - 1) * dk)) for t in (0.0, 0.34, 0.72, 1.0, 1.10, 2.0)}
        f1 = qs[1.0]["fc_rec"]
        sl = np.polyfit(list(qs), [qs[t]["fc_rec"] / f1 for t in qs], 1)[0]
        print("%-5s %-6s %-8s: fc_rec(t=1) %.5f; " % (FLAGNAME, "local" if local else "global", "smoothC4" if smooth else "raw", f1) +
              "  ".join("t=%.2f %+.2f%% (E_sim %.1f)" % (t, 100 * (q["fc_rec"] / f1 - 1), q["Es"]) for t, q in qs.items()) + "  slope %+.4f" % sl)
        out["%s_%d_%d" % (FLAGNAME, local, smooth)] = dict(f1=f1, slope=sl, rel={str(t): q["fc_rec"] / f1 for t, q in qs.items()})
json.dump(out, open(V.OUT + "vpt_bedge.json", "w"), indent=1)
