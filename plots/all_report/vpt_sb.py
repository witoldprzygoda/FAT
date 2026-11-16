#!/usr/bin/env python3
"""vpt_sb.py -- verifier of peaktail: peak side-band estimator (WP6: Y = [(N-CB)_P - R'(N-CB)_SB]/(1 - R' f)) versus the
tail scale t: leakage f(t) = (C4+C5)(SB)/(C4+C5)(P), P content, a_P(t)/a_P(1). Tail scaling per class and per 20 MeV
M(pipi) bin (counts conserved locally) and globally. Output: stdout (tee vpt_sb.log), vpt_sb.json.
"""
import json, sys
import numpy as np
sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import vpt_core as V
import vpt_tmpl as TM
R = TM.rows(); S = R["S"]; w = R["w"]; x = R["m4p"]
RP = 1.1859
P = (x >= 0.53364) & (x < 0.55416)
SB = ((x >= 0.5250) & (x < 0.5336)) | ((x >= 0.5542) & (x < 0.5628))   # defs.h kWin SBlo, SBhi
ppb = np.floor((R["mpp"] - 0.20) / 0.02).astype(int)


def weights(t, flag, local=True):
    wt = w.copy()
    for gi in (2, 3):
        m = R["m45"] & (S["g"] == gi)
        bins = np.unique(ppb[m]) if local else [None]
        for b in bins:
            mm = m if b is None else m & (ppb == b)
            Nn = w[mm].sum(); Nt = w[mm & flag].sum(); Nc = Nn - Nt
            if Nc <= 0:
                continue
            wt[mm & flag] = w[mm & flag] * t
            wt[mm & ~flag] = w[mm & ~flag] * (Nn - t * Nt) / Nc
    return wt


out = {}
m45 = R["m45"]
for fl in ("sub5", "p95", "sub10"):
    for local in (True, False):
        base = None
        rows = []
        for t in (0.0, 0.34, 0.5, 0.72, 1.0, 1.10, 1.5, 2.0):
            wt = weights(t, R["flag"][fl], local)
            Pc = wt[m45 & P].sum(); SBc = wt[m45 & SB].sum(); tot = wt[m45].sum()
            f = SBc / Pc
            q = (1 - RP * f) * Pc
            rows.append((t, f, Pc / tot, q))
        q1 = [r for r in rows if r[0] == 1.0][0][3]
        sl = np.polyfit([r[0] for r in rows], [q1 / r[3] for r in rows], 1)[0]
        print("%-6s %s: " % (fl, "per 20 MeV M(pipi) bin" if local else "global          ") +
              "  ".join("t=%.2f f=%.4f P/all=%.3f a_P x%.3f" % (t, f, pa, q1 / q) for t, f, pa, q in rows) + "  slope %+.3f / unit t" % sl)
        out["%s_%d" % (fl, local)] = dict(rows=rows, slope=sl)
json.dump(out, open(V.OUT + "vpt_sb.json", "w"), indent=1)
