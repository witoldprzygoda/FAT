#!/usr/bin/env python3
# vx_toys.py -- independent verification (prefix vx): READ-ONLY check of the frozen background-only toys of x17c
# (fraction of failed fits, quantiles of the UL, r_hat/sigma) and comparison with the frozen numbers and the observed UL.
import glob, numpy as np, ROOT
ROOT.gROOT.SetBatch(True)
obs = {16.7: 3.0631e-05, 16.85: 3.0646e-05, 17.0: 3.0636e-05, 10.0: 3.892e-05, 25.0: 2.286e-05, 40.0: 8.529e-06}
for m in (16.7, 16.85, 17.0, 10.0, 25.0, 40.0):
    files = sorted(glob.glob("/hdd3/przygoda/HADES/FAT_out/all_report/x17c/toys_m%04d_i0000_*.root" % round(100 * m)))
    files = [f for f in files if not f.endswith(("_9000.root", "_9100.root"))]
    ul = [[], [], []]; ok = [[], [], []]; pull = []; ids = []; nr = []
    for fn in files:
        f = ROOT.TFile.Open(fn); t = f.Get("toys")
        for e in t:
            ids.append(e.id); nr.append(e.nrows)
            for a in range(3):
                ok[a].append(e.ok[a]); ul[a].append(e.ul[a])
            if e.ok[0] and e.rerr[0] > 0: pull.append(e.rhat[0] / e.rerr[0])
        f.Close()
    ids = np.array(ids)
    s = "m %5.2f: %d toys (%d unique ids) from %d files;" % (m, len(ids), len(np.unique(ids)), len(files))
    for a in range(3):
        o = np.array(ok[a], bool); u = np.array(ul[a])[o & np.isfinite(np.array(ul[a]))]
        q = np.quantile(u, [0.02275, 0.15866, 0.5, 0.84134, 0.97725])
        s += " (%d) failed %d, q = %s;" % (a + 1, (~o).sum(), " ".join("%.3e" % v for v in q))
    u0 = np.array(ul[0])[np.array(ok[0], bool)]
    pl = np.array(pull)
    s += " obs UL(1) percentile %.3f; r_hat/sigma mean %+.3f rms %.3f; rows/toy mean %.0f" % ((u0 <= obs[m]).mean(), pl.mean(), pl.std(), np.mean(nr))
    print(s)
