#!/usr/bin/env python3
"""vpt_line2.py -- verifier of peaktail: no-bump line shape with a free peak position shift dm and a free width scale r of
the non-loss resolution (x = s m4_true + r (d - d_lep) + d_lep + dm), tail scale t profiled; also t versus the tail
definition. Output: stdout (tee vpt_line2.log), vpt_line2.json.
"""
import json, math, sys
import numpy as np
sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import vpt_core as V
import vpt_tmpl as TM
exec(open("/home/przygoda/HADES/FAT/FAT/plots/all_report/vpt_line.py").read().split("out = {}")[0])
Ss = R["S"]
base = V.S_NOM * Ss["m4t"].astype(float)
dd = R["d"]; dl = R["dlep"]
out = {}
for rn, sl in (("0.42-0.46", [(0.42, 0.44), (0.44, 0.46)]), ("0.415-0.46", [(0.415, 0.42), (0.42, 0.44), (0.44, 0.46)])):
    a0, b0 = sl[0][0], sl[-1][1]
    N, CB, COV = dspec(sl); sig = N - CB
    regs = TM.in_pp(R["mpp"], a0, b0)
    for lo, hi, od in ((0.49, 0.60, 2), (0.48, 0.62, 3)):
        best = None; tab = {}
        for r in (0.85, 0.90, 0.95, 1.00, 1.05):
            for dm in np.arange(-3.0, 4.01, 0.5):
                x = np.where(Ss["m4t"] > 0, base + r * (dd - dl) + dl, R["m4p"]) + dm * 1e-3
                T = TM.templates(x, gm.e, regs, R["flag"]["sub5"])
                tmin, tl, th, c2, ch, tg = profile(sig, CB, COV, T, lo, hi, od, tg=np.linspace(-1.5, 4.0, 56))
                tab[(r, round(dm, 1))] = (tmin, tl, th, c2)
                if best is None or c2 < best[3]:
                    best = (tmin, tl, th, c2, r, dm)
        c0 = tab[(1.00, 0.0)]
        print("%s %.2f-%.2f pol%d: nominal (r 1, dm 0): t %.2f [%.2f, %.2f] chi2 %.1f | best: r %.2f dm %+.1f MeV: t %.2f [%.2f, %.2f] chi2 %.1f (Delta chi2 %.1f)"
              % (rn, lo, hi, od, c0[0], c0[1], c0[2], c0[3], best[4], best[5], best[0], best[1], best[2], best[3], c0[3] - best[3]))
        # t profiled over dm at r = 1 and over r at the best dm; and t marginal over (r, dm) by chi2 weighting
        cs = np.array([v[3] for v in tab.values()]); ts = np.array([v[0] for v in tab.values()])
        wts = np.exp(-0.5 * (cs - cs.min()))
        print("    chi2-weighted t over the (r, dm) grid: %.2f, sd %.2f;  t at r = 1 for dm = -1, 0, +1, +2, +3: %s; t at dm = best for r = 0.85..1.05: %s" % (
            np.average(ts, weights=wts), math.sqrt(np.average((ts - np.average(ts, weights=wts)) ** 2, weights=wts)),
            " ".join("%.2f" % tab[(1.0, d)][0] for d in (-1.0, 0.0, 1.0, 2.0, 3.0)),
            " ".join("%.2f" % tab[(r, round(best[5], 1))][0] for r in (0.85, 0.90, 0.95, 1.00, 1.05))))
        print("    chi2 at r = 1 for dm = -1..+3: %s" % " ".join("%.1f" % tab[(1.0, d)][3] for d in (-1.0, 0.0, 1.0, 2.0, 3.0)))
        out["%s_%g_%g_%d" % (rn, lo, hi, od)] = dict(nom=c0, best=best, tw=float(np.average(ts, weights=wts)))
json.dump(out, open(V.OUT + "vpt_line2.json", "w"), indent=1, default=float)
