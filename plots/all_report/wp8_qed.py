#!/usr/bin/env python3
# ============================================================================
# wp8_qed.py  --  WP8 (report_all): generator-level check of the pi0 Dalitz pair mass of class C2
# ============================================================================
# SPEC report_all, WP8 "eps2 QED check at M(ee) >= 20 MeV", part 1 (generator level). Written by WP8 (prefix wp8),
# 2026-11-14. Pure computation, no input files.
#
# Class C2 (eta -> pi+pi-pi0, pi0 -> e+e-gamma) comes from the STANDARD SMASH set, where the pi0 decays in HGeant.
# HGeant source (hgeant2-7.4, geant/gdecay.F, gdalet.F, uhbook.F; the local copy in /home/przygoda/hydra2/hydra2-git):
#   * gdecay.F: pi0 (IPART 7) -> GDALET(..., 9995) instead of the phase-space GDECA3;
#   * uhbook.F: histogram 9995 = 2000 bins of 0.1 MeV in q = M(ee) 0-200 MeV, content at the bin centre
#     1e5 x vdmform(q, 134.9766, 1) (vdm = .true.):  vdmform = 8.2e-4 (1 - q^2/M^2)^3 / (q/M) / (1 - q^2/m_rho^2)^2
#     x sqrt(1 - 4 m_e^2/q^2) (1 + 2 m_e^2/q^2), m_rho = 770 MeV, m_e = 0.511 MeV, zero below 2 m_e;
#   * gdalet.F: q = HRNDM1(9995) (sampled from the histogram, uniform inside a bin), y uniform in [-beta, beta]
#     => the gamma* -> e+e- decay angle is ISOTROPIC (QED: 1 + cos^2 + (1 - beta^2) sin^2).
# So the HGeant q spectrum is Kroll-Wada (per dq) x a VDM form factor (rho pole), binned at 0.1 MeV.
#
# Compared here (fractions of all produced Dalitz decays):
#   (a) QED Kroll-Wada (no form factor)
#   (b) Kroll-Wada x |F|^2, F = 1 + a q^2/m_pi0^2, PDG slope a = 0.0332 +- 0.0029   <- reference ("QED with the pi0 form factor")
#   (c) HGeant histogram 9995 as sampled by HRNDM1 (C2 generator)
# Output: docs/latex/tables_all/wp8_qed.txt (+ machine lines 'qed ...')
# Run: python3 -I plots/all_report/wp8_qed.py
# ============================================================================
import math
import numpy as np
from scipy import integrate

ME = 0.51099895          # MeV
MPI0 = 134.9768          # MeV
A_PDG, A_ERR = 0.0332, 0.0029
OUT = "/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all/wp8_qed.txt"


def kw(q, a=None, mrho=None, mass=MPI0):
    """Kroll-Wada dGamma/dq (arbitrary normalisation) with an optional form factor."""
    if q <= 2 * ME or q >= mass:
        return 0.0
    x = (q / mass) ** 2
    beta = math.sqrt(1.0 - 4 * ME * ME / (q * q))
    f = (1.0 - x) ** 3 / q * (1.0 + 2 * ME * ME / (q * q)) * beta
    if a is not None:
        f *= (1.0 + a * x) ** 2
    if mrho is not None:
        f *= 1.0 / (1.0 - (q / mrho) ** 2) ** 2
    return f


def frac_analytic(lo, hi, **kw_args):
    """fraction of decays with lo <= q < hi (MeV)"""
    pts = [2 * ME + 1e-9, 2 * ME * 1.0001, 1.1, 1.5, 2.0, 3.0, 5.0, 10.0, 20.0, 50.0, MPI0 - 1e-9]
    def integ(a, b):
        a = max(a, 2 * ME + 1e-12)
        b = min(b, MPI0)
        if b <= a:
            return 0.0
        brk = [p for p in pts if a < p < b]
        s = 0.0
        edges = [a] + brk + [b]
        for u, v in zip(edges[:-1], edges[1:]):
            s += integrate.quad(lambda q: kw(q, **kw_args), u, v, limit=400, epsabs=0, epsrel=1e-10)[0]
        return s
    tot = integ(0, MPI0)
    return integ(lo, hi) / tot


def hgeant_hist():
    """histogram 9995 of uhbook.F: 2000 bins of 0.1 MeV, value at the bin centre (low edge + 0.05)"""
    lo = np.arange(2000) * 0.1
    c = np.zeros(2000)
    for i in range(2000):
        m = lo[i] + 0.05
        if m < 2 * 0.511 or m >= 134.9766:
            continue
        f = math.sqrt(1. - 4. * (0.511 / m) ** 2) * (1. + 2. * (0.511 / m) ** 2)
        dalit = ((1. - (m / 134.9766) ** 2) ** 3 / (m / 134.9766)) * 1. / (1. - (m / 770.0) ** 2) ** 2
        c[i] = 1e5 * 8.2e-4 * dalit * f
    return lo, c


def frac_hist(lo_e, c, a, b):
    """fraction with a <= q < b for a histogram sampled uniformly inside each 0.1 MeV bin"""
    tot = c.sum()
    s = 0.0
    for i in range(len(c)):
        u, v = lo_e[i], lo_e[i] + 0.1
        ov = max(0.0, min(v, b) - max(u, a))
        s += c[i] * ov / 0.1
    return s / tot


def main():
    lo_e, c = hgeant_hist()
    cuts = [2, 8, 10, 17, 20, 30, 50, 100]
    classes = [("<8", 0, 8), ("8-20", 8, 20), ("20-30", 20, 30), (">=30", 30, 1e9), ("30-100", 30, 100),
               (">=10 (E2)", 10, 1e9), (">=20", 20, 1e9), ("15.5-24.5", 15.5, 24.5)]
    models = [("QED Kroll-Wada, no form factor", dict()),
              ("Kroll-Wada x (1 + a x)^2, a = 0.0332 (PDG)", dict(a=A_PDG)),
              ("Kroll-Wada x (1 + a x)^2, a = 0.0303", dict(a=A_PDG - A_ERR)),
              ("Kroll-Wada x (1 + a x)^2, a = 0.0361", dict(a=A_PDG + A_ERR)),
              ("Kroll-Wada x rho pole (770 MeV), exact", dict(mrho=770.0))]
    res = {}
    for name, kwa in models:
        res[name] = {cl[0]: frac_analytic(cl[1], cl[2], **kwa) for cl in classes}
        res[name].update({f">{q}": frac_analytic(q, 1e9, **kwa) for q in cuts})
    hname = "HGeant histogram 9995 (C2 generator)"
    res[hname] = {cl[0]: frac_hist(lo_e, c, cl[1], cl[2]) for cl in classes}
    res[hname].update({f">{q}": frac_hist(lo_e, c, q, 1e9) for q in cuts})
    # probability below 2 m_e in the HGeant histogram (bin 1.0-1.1 MeV straddles the threshold)
    below = frac_hist(lo_e, c, 0, 2 * ME)
    ref = "Kroll-Wada x (1 + a x)^2, a = 0.0332 (PDG)"
    with open(OUT, "w") as f:
        def P(s=""):
            print(s)
            f.write(s + "\n")
        P("# wp8_qed.txt -- plots/all_report/wp8_qed.py (WP8, 2026-11-14): generator-level pi0 Dalitz pair mass q = M(ee)")
        P("# fractions of ALL produced pi0 -> e+e-gamma decays; HGeant = histogram 9995 of hgeant2-7.4 uhbook.F (KW x rho-pole")
        P("# VDM form factor, 0.1 MeV bins, value at the bin centre, sampled uniformly inside a bin by HRNDM1)")
        P("# reference = Kroll-Wada with the PDG pi0 transition form factor slope a = 0.0332 +- 0.0029 (F = 1 + a q^2/m_pi0^2)")
        P("")
        hdr = "%-46s" % "model" + "".join("%11s" % k for k in [cl[0] for cl in classes])
        P(hdr)
        for name in [m[0] for m in models] + [hname]:
            P("%-46s" % name + "".join("%11.5f" % res[name][cl[0]] for cl in classes))
        P("")
        P("%-46s" % "fraction above q [MeV]" + "".join("%9s" % (">" + str(q)) for q in cuts))
        for name in [m[0] for m in models] + [hname]:
            P("%-46s" % name + "".join("%9.5f" % res[name][">" + str(q)] for q in cuts))
        P("")
        P("HGeant / reference (ratio of the fractions):")
        P("%-46s" % "" + "".join("%11s" % cl[0] for cl in classes))
        P("%-46s" % "HGeant / KW x PDG form factor" + "".join("%11.5f" % (res[hname][cl[0]] / res[ref][cl[0]]) for cl in classes))
        P("%-46s" % "KW no FF / KW x PDG form factor" + "".join("%11.5f" % (res[models[0][0]][cl[0]] / res[ref][cl[0]]) for cl in classes))
        P("")
        P("HGeant probability below 2 m_e (bin 1.0-1.1 MeV straddles the threshold): %.2e" % below)
        P("form-factor slope uncertainty on the >=20 MeV fraction: a -+ 0.0029 -> %.5f / %.5f (ref %.5f)" %
          (res[models[2][0]][">=20"], res[models[3][0]][">=20"], res[ref][">=20"]))
        P("")
        P("READING. The C2 generator (HGeant) samples the Kroll-Wada pair mass with a rho-pole form factor. Its fraction of")
        P("decays with M(ee) >= 20 MeV differs from Kroll-Wada x the PDG form factor by the ratio printed above; the q")
        P("spectrum of C2 is therefore QED at the per-mille level. What HGeant does NOT reproduce is the decay angle of")
        P("gamma* -> e+e- (isotropic instead of 1 + cos^2 theta + (1 - beta^2) sin^2 theta); see wp8_results.md (toy estimate).")
        P("")
        for name in [m[0] for m in models] + [hname]:
            P("qed %s | " % name.replace(" ", "_") + " ".join("%s=%.6f" % (cl[0].replace(" ", ""), res[name][cl[0]]) for cl in classes))
    return res


if __name__ == "__main__":
    main()
