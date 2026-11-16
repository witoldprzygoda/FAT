#!/usr/bin/env python3
"""vpt_ucons.py -- verifier of peaktail: u-fit tail scale with RAW templates and the peak scale constrained by the side
bands, A = a_P(t) +- 0.080 a_P(t)/a_P(1) (a_P(t) from vpt_sb.py logic, sub5, raw C4); f_eta(u-edge) at the profile
minimum (t and A free with the constraint). Output: stdout (tee vpt_ucons.log).
"""
import sys
import numpy as np
sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import vpt_core as V
import vpt_tmpl as TM
exec(open("/home/przygoda/HADES/FAT/FAT/plots/all_report/vpt_smooth.py").read().split('print("=" * 110)')[0])
ts = np.linspace(-1.5, 4.0, 111)
# raw a_P(t) table (sub5, per 20 MeV M(pipi) bin as vpt_sb)
x = R["m4p"]; w = R["w"]
Pm = (x >= 0.53364) & (x < 0.55416); SBm = ((x >= 0.5250) & (x < 0.5336)) | ((x >= 0.5542) & (x < 0.5628))
ppb = np.floor((R["mpp"] - 0.20) / 0.02).astype(int)
def apr_raw(t, fl):
    Pc = SBc = 0.0
    for gi in (2, 3):
        m = R["m45"] & (S["g"] == gi)
        for b in np.unique(ppb[m]):
            mm = m & (ppb == b); Nn = w[mm].sum(); Nt = w[mm & fl].sum(); Nc = Nn - Nt
            wt = np.where(fl, w * t, w * (Nn - t * Nt) / Nc)
            Pc += wt[mm & Pm].sum(); SBc += wt[mm & SBm].sum()
    return (1 - 1.1859 * SBc / Pc) * Pc
for flag in ("sub5", "p95"):
    fl = R["flag"][flag]
    q1 = apr_raw(1.0, fl)
    apr = np.array([q1 / apr_raw(t, fl) for t in ts])
    reg = TM.in_pp(R["mpp"], 0.30, 0.40) & np.isfinite(R["u"])
    Tt = TM.templates(R["u"], g.e, reg, fl)
    for lo, hi, od in ((-0.045, 0.045, 2), (-0.045, 0.055, 2), (-0.045, 0.060, 2), (-0.045, 0.070, 3), (-0.055, 0.070, 3)):
        i = np.where((g.x >= lo - 1e-12) & (g.x < hi - 1e-12))[0]
        rows_ = {}
        for name, mk in (("a_P fixed", lambda ap: (V.AP, 0)), ("a_P(t) fixed", lambda ap: (V.AP * ap, 0)),
                         ("a_P(t)+-0.08", lambda ap: (V.AP * ap, AP_ERR * ap)), ("A free", lambda ap: (V.AP, np.inf))):
            ch = []; fs = []
            for t, ap in zip(ts, apr):
                A0, sA = mk(ap)
                c2, h, A, Lh = fit_c(sig, TM.K_of_t(Tt, t), A0, sA, lo, hi, od)
                ch.append(c2); fs.append(h / float(Lh @ BUMP[i]))
            tm, tl, th, cmin = prof(ch, ts)
            fmin = float(np.interp(tm, ts, fs)); flo = float(np.interp(tl, ts, fs)); fhi = float(np.interp(th, ts, fs))
            rows_[name] = (tm, tl, th, fmin, flo, fhi)
        print("%-5s %+3.0f..%+3.0f pol%d: " % (flag, lo * 1e3, hi * 1e3, od) + " | ".join(
            "%s: t %.2f [%.2f, %.2f] f %.4f [%.4f, %.4f]" % (k, *v) for k, v in rows_.items()))
