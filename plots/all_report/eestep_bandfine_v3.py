#!/usr/bin/env python3
# ============================================================================
# eestep_band_v3.py  --  ruling D38 (15 November 2026): the SMASH bump prediction f x h_sim in the large-M(ee) bands of the
#   eestep CB test, with the per-row C2 QED weights (c2recon) and f_eta(u-edge) = 0.42574 (edgeu_ana_v3.py), instead of the
#   class factors 0.879 / 0.8825 and f = 0.4200. The band estimator is eestep_ana.py's band_fit, copied unchanged; the data
#   side (h(N), h(CB), h(N - CB)) does not depend on the factor. Only the signal prediction and the derived z, k_cross change.
#   Reads the edgeu columns and c2recon_weights.npz; writes only /hdd3/.../all_report/edgeu_v3/eestep_band_v3.log.
# Run: python3 -I plots/all_report/eestep_band_v3.py   (about 1 min)
# ----------------------------------------------------------------------------
# eestep_bandfine_v3.py (2026-11-15, FIGURE AGENT C of the fourth integration): copy of eestep_band_v3.py for the figure
#   eestep_bands_v3: the same band estimator (band_fit, unchanged) evaluated in the fine M(ee) bands of panel (b) of
#   eestep_bands (50-70, 70-90, 90-110, 110-m_pi0 MeV/c^2) and in the band >= 90 MeV/c^2 of its annotation, with the
#   class factors of D35 (check: must reproduce h_sim of eestep_ana.json) and the per-count QED weights of D38.
#   Writes only the new files /hdd3/przygoda/HADES/FAT_out/all_report/eestep_v3/eestep_bandfine_v3.{log,json}.
# Run: python3 -I plots/all_report/eestep_bandfine_v3.py   (about 1 min)
# ============================================================================
import math, sys
import numpy as np
sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import edgeu_core as C
import eestep_core as E

import json, os
OUTF = "/hdd3/przygoda/HADES/FAT_out/all_report/eestep_v3/"
os.makedirs(OUTF, exist_ok=True)
LOG = open(OUTF + "eestep_bandfine_v3.log", "w")
def P(*a):
    s = " ".join(str(x) for x in a); print(s); LOG.write(s + "\n"); LOG.flush()
D = C.load_data(); uD = C.u_data(D)
S = C.load_sim(); uS, mppS, _ = C.u_sim(S)
wS = C.sim_weight(S, 108, "S0", "A", c2corr=False)
s0D = (D["bits"] & 1) != 0
SL = C.slices_of(0.30, 0.40)
pbD = C.ppbin_of(D["mpipi"]); pbS = C.ppbin_of(mppS.astype(np.float32))
G = E.G; T = E.T
A_P = C.AP
# per-row weights (as edgeu_ana_v3.py)
Wt = np.load("/hdd3/przygoda/HADES/FAT_out/all_report/c2recon/c2recon_weights.npz")
key_w = {(int(Wt["tree"][i]), int(Wt["grp"][i]), float(np.float32(Wt["m4"][i])), float(np.float32(Wt["m4_true"][i]))): i for i in range(len(Wt["tree"]))}
ISC2 = S["g"] == C.G["C2"]
QROW = np.ones(S["n"])
for i in np.where(ISC2)[0]:
    QROW[i] = Wt["q_A"][key_w[(int(S["tree"][i]), int(S["grp"][i]), float(S["m4"][i]), float(S["m4t"][i]))]]
def hk(k, mask):
    out = []
    b = G.bin(uD)
    for (a, c) in SL:
        m = mask & s0D & (D["k"] == k) & (pbD >= a) & (pbD < c) & (b >= 0)
        out.append(np.bincount(b[m], minlength=G.nb).astype(np.float64))
    return out
def hsim(groups, mask, mode):
    b = G.bin(uS)
    m = mask & np.isin(S["g"], [C.G[g] for g in groups]) & (pbS >= SL[0][0]) & (pbS < SL[-1][1]) & (b >= 0) & (wS != 0)
    w = wS.copy()
    if mode == "D35":
        w = np.where(ISC2 & (S["mee"] < 0.010), w * C.CF_E1, w)
        w = np.where(ISC2 & ~(S["mee"] < 0.010), w * C.CF_E2, w)
    elif mode == "D38":
        w = w * QROW
    return np.bincount(b[m], weights=w[m], minlength=G.nb)
def band_fit(dmask, smask):
    Nn = sum(hk(0, dmask)); pps = hk(1, dmask); mms = hk(2, dmask)
    CB = np.zeros(G.nb); COV = np.zeros((G.nb, G.nb))
    for p_, m_ in zip(pps, mms):
        cb, cov = C.cb_rule(p_, m_, G, 1.0, 0)
        CB += cb; COV += cov
    sp = C.Spec(Nn, CB, COV, sum(pps), sum(mms))
    K = hsim(["C4", "C5"], smask, "none")
    R = C.fit_step(sp, T, G, -0.045, 0.045, 2, fixed=A_P * K)
    hs = {md: C.h_lin(R, hsim(["C2", "C3"], smask, md)) for md in ("D35", "D38")}
    hcb = C.h_lin(R, CB); hN = C.h_lin(R, Nn - A_P * K)
    Lh = R.Lh; var_model = np.maximum(R.model + R.fx + sp.cb[R.idx], 1.0)
    ehN = float(np.sqrt((Lh ** 2 * var_model).sum())); ehCB = float(np.sqrt(Lh @ COV[np.ix_(R.idx, R.idx)] @ Lh))
    return dict(hd=R.h, ehd=R.eh, hs=hs, hcb=hcb, hN=hN, ehN=ehN, ehCB=ehCB)
F = dict(D35=0.4200, D38=0.42574)
allr = band_fit(np.ones(D["n"], bool), np.ones(S["n"], bool))
P("# eestep_band_v3.py: all M(ee): h(N - CB) %.2f; h_sim D35 %.2f (f %.5f), D38 %.2f (f %.5f)" % (allr["hd"], allr["hs"]["D35"], allr["hd"] / allr["hs"]["D35"],
  allr["hs"]["D38"], allr["hd"] / allr["hs"]["D38"]))
OUTJ = dict(F=F, all=dict(hd=allr["hd"], hs=allr["hs"]), bands=[])
for (a, b) in ((0.050, 0.070), (0.070, 0.090), (0.090, 0.110), (0.110, C.MPI0)):
    md = (D["mee"] >= a) & (D["mee"] < b); ms = (S["mee"] >= a) & (S["mee"] < b)
    x = band_fit(md, ms)
    OUTJ["bands"].append(dict(lo=a, hi=b, hd=x["hd"], ehd=x["ehd"], hcb=x["hcb"], hN=x["hN"], hs_D35=x["hs"]["D35"], hs_D38=x["hs"]["D38"]))
    P("fine band %.0f-%.0f MeV: h(N) %.2f, h(CB) %.2f, h(N - CB) %.2f +- %.2f; h_sim D35 %.3f, D38 %.3f (ratio %.4f)" %
      (1e3 * a, 1e3 * b, x["hN"], x["hcb"], x["hd"], x["ehd"], x["hs"]["D35"], x["hs"]["D38"], x["hs"]["D38"] / x["hs"]["D35"]))
for (a, b) in ((0.090, C.MPI0),):
    md = (D["mee"] >= a) & (D["mee"] < b); ms = (S["mee"] >= a) & (S["mee"] < b)
    x = band_fit(md, ms)
    OUTJ["r90"] = dict(hd=x["hd"], ehd=x["ehd"], hcb=x["hcb"], hN=x["hN"], ehN=x["ehN"], ehCB=x["ehCB"], hs_D35=x["hs"]["D35"], hs_D38=x["hs"]["D38"])
    for mode in ("D35", "D38"):
        fh = F[mode] * x["hs"][mode]
        z_phys = (x["hd"] - fh) / x["ehd"]; z_art = (x["hd"] - (fh - x["hcb"])) / x["ehd"]
        k = (x["hN"] - fh) / x["hcb"]; ek = math.sqrt(x["ehN"] ** 2 + (k * x["ehCB"]) ** 2) / x["hcb"]
        P("M(ee) >= %.0f MeV [%s]: h(N) %.1f +- %.1f, h(CB) %.1f, h(N - CB) %.1f +- %.1f, f x h_sim %.2f (h_sim %.2f); z (CB rule) %+.2f; "
          "artefact N - CB = %.1f, z %+.2f; k_cross %.2f +- %.2f; band window factor h_sim(D38)/h_sim(D35) %.4f" %
          (1e3 * a, mode, x["hN"], x["ehN"], x["hcb"], x["hd"], x["ehd"], fh, x["hs"][mode], z_phys, fh - x["hcb"], z_art, k, ek, x["hs"]["D38"] / x["hs"]["D35"]))
for mode in ("D35", "D38"):
    r = OUTJ["r90"]; fh = F[mode] * r["hs_" + mode]
    r["fhs_" + mode] = fh; r["z_art_" + mode] = (r["hd"] - (fh - r["hcb"])) / r["ehd"]; r["z_phys_" + mode] = (r["hd"] - fh) / r["ehd"]
json.dump(OUTJ, open(OUTF + "eestep_bandfine_v3.json", "w"), indent=1, default=float)
P("written", OUTF + "eestep_bandfine_v3.json")
