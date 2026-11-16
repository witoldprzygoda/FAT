#!/usr/bin/env python3
# ============================================================================
# veestep_sim.py  --  INDEPENDENT verification of eestep (prefix veestep): simulated like-sign cross pairs, SMASH x f,
#                     material scaling, ring / energy-loss test of merged conversion pairs, peak-window content
# ============================================================================
# Written 2026-11-15 by the independent verifier. Inputs: edgeu simls_* columns (own loader), veestep_dump.C and
# veestep_pat2.C columns. The simulation is used ONLY to identify origins and to compare the size of the like-sign step;
# no CB or k number is taken from it (rule of the simulation study).
# Run: python3 -I plots/all_report/veestep_sim.py  (-> /hdd3/przygoda/HADES/FAT_out/all_report/veestep/veestep_sim.log)
# ============================================================================
import sys

import numpy as np

sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import veestep_core as V

LOG = open(V.OUT + "veestep_sim.log", "w")


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s); LOG.write(s + "\n"); LOG.flush()


np.seterr(all="ignore")
F = 0.4193   # own estimator (veestep_ana.py section 1); nominal 0.4200
L = V.load_sim_ls()
ok = ((L["bits"] & 1) != 0) & ((L["bits"] & 8) != 0)
w = np.where(ok, L["w"] * L["f108"], 0.0)
inpp = (L["mpp"] >= 0.30) & (L["mpp"] < 0.40) & np.isfinite(L["u"])


def org(proc, par, gpar):
    o = np.full(proc.shape, "X", dtype=object)
    o[(proc == 5) & (par == 7) & (gpar == 17)] = "D(eta)"
    o[(proc == 5) & (par == 7) & (gpar != 17)] = "D(other pi0)"
    o[(proc == 6) & (par == 1)] = "C"
    o[(proc == 7) & (par == 1)] = "Compton"
    o[proc == 10] = "delta"
    return o


o1 = org(L["proc1"], L["par1"], L["gpar1"]); o2 = org(L["proc2"], L["par2"], L["gpar2"])
pair = np.array(["+".join(sorted([a, b])) for a, b in zip(o1, o2)])
P("# veestep_sim.py, 2026-11-15 (std set: trees 1 e+e+, 2 e-e-; S0, ownA, w x f108; M(pipi)' 0.30-0.40; pions of one eta, pcls 2)")
P("  'edge-consistent' = true u = m4_true - M4max(mpipi_true, mee_true) <= +2 MeV (the leptons can come from the eta's pi0)")
P("\n== 1. origin of the simulated like-sign rows (pcls 2)")
for t, nm in [(2, "e-e-"), (1, "e+e+")]:
    m = (L["tree"] == t) & (L["pcls"] == 2) & ok & inpp
    P("  %s: rows %d, weighted %.1f expected raw rows" % (nm, m.sum(), w[m].sum()))
    for p in sorted(set(pair[m]), key=lambda x: -w[m & (pair == x)].sum()):
        mm = m & (pair == p)
        ec = mm & (L["utrue"] <= 0.002)
        near = ec & (np.abs(L["u"]) < 0.045)
        P("    %-26s rows %4d (w %7.1f); edge-consistent %3d (w %6.1f); of them |u| < 45 MeV: %3d (w %6.1f)" %
          (p, mm.sum(), w[mm].sum(), ec.sum(), w[ec].sum(), near.sum(), w[near].sum()))
# e- only processes
ism = (L["tree"] == 2) & (L["pcls"] == 2) & ok & inpp
cd = np.array([("Compton" in p) or ("delta" in p) for p in pair])
P("  share of the e-e- weight with a Compton or delta electron: all %.2f; edge-consistent |u| < 45 MeV %.2f" %
  (w[ism & cd].sum() / w[ism].sum(), w[ism & cd & (L["utrue"] <= 0.002) & (np.abs(L["u"]) < 0.045)].sum() /
   w[ism & (L["utrue"] <= 0.002) & (np.abs(L["u"]) < 0.045)].sum()))
pt = np.array([all(x in ("D(eta)", "D(other pi0)", "C") for x in p.split("+")) for p in pair])
for nm, mk in [("pair-type only (Dalitz / conversion legs)", pt), ("all pcls 2", np.ones(pt.size, bool))]:
    a = w[(L["tree"] == 1) & (L["pcls"] == 2) & ok & inpp & mk].sum(); b = w[ism & mk].sum()
    P("  e+e+/e-e- weighted, %s: %.3f" % (nm, a / b))

P("\n== 2. SMASH x f against the data step (data estimator: pol2 + step, weights from the data e-e- spectrum)")
D = V.load_data()
base = D["S0"] & np.isfinite(D["u"])
N, PP, MM, CB = V.spectra(D, base)
r = V.fit_pois(MM)
LM = V.Lin(np.maximum(np.interp(V.XC, V.XC[V.FIT], r["mu"]), 1))
rng = np.random.default_rng(11)
data = {0.0: (73.5, 15.6), 0.090: (57.0, 9.0)}
for t, nm in [(2, "e-e-"), (1, "e+e+")]:
    for sel_nm, sel in [("pcls 2", L["pcls"] == 2), ("pcls 2, edge-consistent", (L["pcls"] == 2) & (L["utrue"] <= 0.002))]:
        for q in [0.0, 0.090]:
            m = (L["tree"] == t) & sel & inpp & ok & (L["mee"] >= q)
            y = V.hist(L["u"][m], w[m])
            h = LM.h(y)
            uq, inv = np.unique(L["grp"][m].astype(np.int64), return_inverse=True)
            hb = [LM.h(V.hist(L["u"][m], w[m] * rng.poisson(1, uq.size)[inv])) for _ in range(500)]
            below = m & (L["u"] > -0.020) & (L["u"] <= 0)
            dens = w[below].sum() / 8
            nb = below.sum()
            line = "  %s %-25s M(ee) >= %3.0f: rows |u|<45 %3d; step x f %6.1f +- %5.1f; density below the edge x f %5.1f (from %d rows, +- %.0f %%)" % (
                nm, sel_nm, q * 1000, (m & (np.abs(L["u"]) < 0.045)).sum(), F * h, F * np.std(hb), F * dens, nb, 100 / np.sqrt(max(nb, 1)))
            if t == 2:
                dv, de = data[q]
                sh = F * np.std(hb)
                line += " | data %.1f +- %.1f: ratio (step) %.2f, z %.1f; ratio (density) %.2f, z %.1f" % (
                    dv, de, dv / max(F * h, 1e-9), (dv - F * h) / np.hypot(de, sh), dv / max(F * dens, 1e-9),
                    (dv - F * dens) / np.hypot(de, F * dens / np.sqrt(max(nb, 1))))
            P(line)

P("\n== 3. material scaling needed for the excess (e-e-, edge-consistent, |u| < 45 MeV; number of photon interactions 0/1/2)")
nmat = lambda pr, pa: np.where((pr == 5) & (pa == 7), 0, np.where(np.isin(pr, [6, 7]) & (pa == 1), 1, np.where(pr == 10, 1, 9)))
nm_ = nmat(L["proc1"], L["par1"]) + nmat(L["proc2"], L["par2"])
m = ism & (L["utrue"] <= 0.002) & (np.abs(L["u"]) < 0.045)
sh = [w[m & (nm_ == k)].sum() / w[m].sum() for k in range(3)]
P("  shares %.2f / %.2f / %.2f (rows %d)" % (sh[0], sh[1], sh[2], m.sum()))
for R in [2.3, 3.5]:
    a, b, c = sh[2], sh[1], sh[0] - R
    X = (-b + np.sqrt(b * b - 4 * a * c)) / (2 * a)
    P("  data/SMASH = %.1f needs photon-interaction probability x %.2f; the bump's own conversion scale kappa = 0.88 +- 0.21 +- 0.35 -> %.1f sigma" %
      (R, X, (X - 0.88) / np.hypot(0.21, 0.35)))

P("\n== 4. merged conversion pairs (one track with the whole photon): ring size and MDC energy loss of the e-e- step rows")
X2 = {}
for v in ["rich_padnum", "dedx_mdc"]:
    X2[v] = (np.fromfile(V.OUT + "d2_%s_1.f32" % v, np.float32), np.fromfile(V.OUT + "d2_%s_2.f32" % v, np.float32))
sel0 = D["S0"] & (D["mpipi"] >= 0.30) & (D["mpipi"] < 0.40) & (D["k"] == 2) & (D["mee"] >= 0.070) & np.isfinite(D["u"])
win = sel0 & (np.abs(D["u"]) < 0.045)
for v, (a, b) in X2.items():
    for lab, z in [("larger of the two", np.maximum(a, b)), ("ratio larger/smaller", np.maximum(a, b) / np.maximum(np.minimum(a, b), 1e-6))]:
        t = np.percentile(z[win], [33.33, 66.67])
        out = []
        for lo, hi in [(-1e30, t[0]), (t[0], t[1]), (t[1], 1e30)]:
            mm = sel0 & (z >= lo) & (z < hi)
            rr = V.fit_pois(V.hist(D["u"][mm]))
            out.append("%5.1f +- %4.1f (share %.2f)" % (rr["h"], rr["eh"], rr["h"] / max((mm & win).sum() / 36, 1)))
        P("  %-12s %-22s tertiles (cuts %.3g, %.3g): " % (v, lab, t[0], t[1]) + " | ".join(out))
P("  (a merged pair would carry about twice the ring and energy-loss signal: the step would sit in the top tertile)")

P("\n== 5. cross pairs in the peak windows (S0, 0.26 <= M(pipi) < 0.46, e-e- rows with u <= 0)")
b = D["S0"] & (D["mpipi"] >= 0.30) & (D["mpipi"] < 0.40) & (D["mee"] >= 0.090) & np.isfinite(D["u"]) & (D["k"] == 2)
dens = (b & (D["u"] > -0.015) & (D["u"] <= 0)).sum() / 6
share = 57.0 / dens
P("  cross-pair share of the e-e- rows just below the edge (M(ee) >= 90, -15..0 MeV): 57.0 / %.1f = %.2f" % (dens, share))
for nm, lo, hi in [("SBlo 0.5250-0.5336", 0.5250, 0.5336), ("R2 0.5225-0.5336", 0.5225, 0.53364), ("P 0.53364-0.55416", 0.53364, 0.55416),
                   ("SBhi 0.5542-0.5628", 0.5542, 0.5628)]:
    w_ = D["S0"] & (D["mpipi"] >= 0.26) & (D["mpipi"] < 0.46) & (D["m4"] >= lo) & (D["m4"] < hi) & (D["k"] == 2)
    nb = (w_ & np.isfinite(D["u"]) & (D["u"] <= 0)).sum()
    P("  %-20s e-e- rows with u <= 0: %4d -> about %.0f cross pairs; CB they carry (x 0.744) about %.0f rows" % (nm, nb, share * nb, 0.744 * share * nb))
P("  expected k_true of the CB rule for these pairs (factorised acceptance, observed rc = 0.127, rB = 0.229): 0.96 (no Compton) to"
  " 1.02 (half Compton): residual in P below 1.5 rows (0.03 % of Y_peak 5088)")
P("done")
