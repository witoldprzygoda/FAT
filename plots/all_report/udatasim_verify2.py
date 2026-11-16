#!/usr/bin/env python3
# ============================================================================
# udatasim_verify2.py  --  INDEPENDENT check of udatasim, part 2
# ============================================================================
# Written 2026-11-15 (prefix udatasim, verifier). New file.
#   1. bin-by-bin comparison of the verifier's own arrays (udatasim_verify_arrays.npz, from udatasim_verify.py) with the
#      arrays that udatasim_figs.py draws (udatasim_figdata.npz; read only, for the comparison);
#   2. own rebuild of the stage-A check in CR-U (M4 0.575-0.700, E2, per slice, (N - k_E2 CB)/(B_A + 0.4 eta));
#   3. size of the peak classes per bin (bound on the D40 tail variation shown on the display).
# Output: /hdd3/przygoda/HADES/FAT_out/all_report/udatasim/udatasim_verify2.log
# ============================================================================
import math

import numpy as np

IN = "/hdd3/przygoda/HADES/FAT_out/all_report/edgeu/"
O = "/hdd3/przygoda/HADES/FAT_out/all_report/udatasim/"
LOG = open(O + "udatasim_verify2.log", "w")


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True); LOG.write(s + "\n")


V = np.load(O + "udatasim_verify_arrays.npz")
F = np.load(O + "udatasim_figdata.npz")
x = V["x"]
ip = np.where((x > -0.100) & (x < 0.100))[0]
assert np.allclose(x[ip], F["x"])
GMAP = {"C2": ["C2"], "C3": ["C3"], "C4": ["C4"], "C5": ["C5"], "Geta": ["Geta"], "GD": ["GDom", "GDnom"], "GC": ["GCom", "GCnom"]}
for tv, tf in (("r1", "r1"), ("allmpp", "all")):
    for c in ("all", "E1", "E2"):
        a = "%s_%s_" % (tv, c); b = "%s_%s_" % (tf, c)
        dN = np.max(np.abs(V[a + "N"][ip] - F[b + "N"])); dC = np.max(np.abs(V[a + "CB"][ip] - F[b + "cb"])); dS = np.max(np.abs(V[a + "sig"][ip] - F[b + "sig"]))
        dm = {}
        for mk in ("a", "b"):
            for g, gs in GMAP.items():
                mine = V[a + mk + "_" + g][ip]; theirs = sum(F[b + mk + "_" + q] for q in gs)
                dm[mk + g] = float(np.max(np.abs(mine - theirs)) / max(np.max(np.abs(mine)), 1e-9))
        P("%-6s %-3s: max |dN| %.3g, max |dCB| %.3g, max |d(N-CB)| %.3g counts; model parts max rel. diff %.2e" % (tv, c, dN, dC, dS, max(dm.values())))

# ---------------------------------------------------------------- peak classes per bin, model (b)
for tv in ("r1", "allmpp"):
    a = "%s_all_" % tv
    M = sum(V[a + "b_" + g] for g in GMAP)
    pk = V[a + "b_C4"] + V[a + "b_C5"]
    m = (x > -0.1) & (x < 0.1)
    i = np.argmax(np.where(m, pk / M, 0))
    P("%s: largest (C4 + C5)/model (b) in -100..+100: %.4f at u = %+.1f MeV/c2; at u = +16.25: %.4f" % (tv, (pk / M)[i], 1e3 * x[i], (pk / M)[np.argmin(np.abs(x - 0.01625))]))

# ---------------------------------------------------------------- CR-U: own rebuild
MPI0 = 0.1349768
def col(n, dt):
    return np.fromfile(IN + n, dtype=dt)
dk = col("data_k.i8", np.int8); dm4 = col("data_m4.f32", np.float32); dpp = col("data_mpipi.f32", np.float32)
dee = col("data_mee.f32", np.float32); dbits = col("data_bits.u8", np.uint8)
sg = col("sim_g.i8", np.int8); sb = col("sim_bits.u8", np.uint8); sm4 = col("sim_m4.f32", np.float32); spp = col("sim_mpipi.f32", np.float32)
see = col("sim_mee.f32", np.float32); sw = col("sim_w.f64", np.float64); f108 = col("sim_f108.f32", np.float32)
sn = sg.size
W = np.load("/hdd3/przygoda/HADES/FAT_out/all_report/c2recon/c2recon_weights.npz")
st = col("sim_tree.i8", np.int8); sgrp = col("sim_grp.i32", np.int32)
kw = {(int(W["tree"][i]), int(W["grp"][i]), float(np.float32(W["mee"][i])), float(np.float32(W["m4"][i]))): W["q_A"][i] for i in range(len(W["tree"]))}
q = np.ones(sn)
for i in np.where(sg == 0)[0]:
    q[i] = kw[(int(st[i]), int(sgrp[i]), float(see[i]), float(sm4[i]))]
w0 = np.where(((sb & 1) != 0) & ((sb & 8) != 0), sw * f108, 0.0) * q
m4p = 0.99615 * sm4.astype(np.float64); mppp = (0.9952 * spp.astype(np.float64)).astype(np.float32)
SA = open("/hdd3/przygoda/HADES/FAT_out/all_report/wp3/stageA.txt").read().splitlines()
vals = [np.array([float(v) for v in l.split()[1:]]) for l in SA if l.startswith("values ")][0]
CJ, BT, KB, KE2 = vals[:5], vals[5], vals[6], vals[8]
jsl = np.searchsorted(np.array([0.34, 0.37, 0.39, 0.43]), mppp.astype(np.float64), side="right")
cA = CJ[jsl] * (1 + BT * (m4p - 0.55) / 0.1)
isGD = (sg == 5) | (sg == 6); isGC = (sg == 7) | (sg == 8)
BA_w = np.where(isGD, w0 * cA, 0.0) + np.where(isGC, w0 * cA * KB, 0.0)
ETA_w = np.where(sg <= 4, w0, 0.0); GETA_w = np.where(sg == 4, w0, 0.0)
dE2 = ~(dee < 0.010); sE2 = ~(see < 0.010)
dS0 = (dbits & 1) != 0
P("\nCR-U (M4 0.575-0.700), E2, S0; data N - k_E2 CB (k_E2 = %.4f; CB = 2 sqrt(N++ N--) per slice in 10 MeV/c2 M4 cells, E2 like-sign counts), model B_A + 0.4 eta (raw templates)" % KE2)
for lab, lo, hi in (("S1", 0.27, 0.34), ("S2", 0.34, 0.37), ("S3", 0.37, 0.39), ("S4", 0.39, 0.43), ("S1-S4", 0.27, 0.43)):
    edges = [e for e in (0.27, 0.34, 0.37, 0.39, 0.43) if lo - 1e-9 <= e <= hi + 1e-9]
    num = 0.0; Nsum = 0.0; CBs = 0.0
    for a, b in zip(edges[:-1], edges[1:]):
        base = dS0 & dE2 & (dpp >= np.float32(a)) & (dpp < np.float32(b))
        cells = np.arange(0.575, 0.7001, 0.010)   # 0.575-0.705 cells; take the M4 range exactly with fine bins below
        fine = lambda kk: np.histogram(dm4[base & (dk == kk)], bins=np.round(np.arange(0.400, 0.7001, 0.0025), 6))[0].astype(float)
        nU, pU, mU = fine(0), fine(1), fine(2)
        # 10 MeV cells from 0.400, CB per cell spread evenly to its four 2.5 MeV bins (cell constant)
        cbw = 2 * np.sqrt(pU.reshape(-1, 4).sum(1) * mU.reshape(-1, 4).sum(1))
        cbf = np.repeat(cbw / 4, 4)
        xf = 0.400 + 0.0025 * (np.arange(nU.size) + 0.5)
        iu = (xf > 0.575) & (xf < 0.700)
        Nsum += nU[iu].sum(); CBs += cbf[iu].sum()
    ms = sE2 & (mppp >= np.float32(lo)) & (mppp < np.float32(hi)) & (m4p >= 0.575) & (m4p < 0.700)
    BA = BA_w[ms].sum(); ET = ETA_w[ms].sum(); GE = GETA_w[ms].sum()
    d = Nsum - KE2 * CBs
    P("  %-5s: N %.0f, CB %.0f -> (N - k CB)/(B_A + 0.4 eta) = %.4f +- %.4f (data); with G_eta x 0.656: %.4f; G_eta share %.3f" %
      (lab, Nsum, CBs, d / (BA + 0.4 * ET), math.sqrt(Nsum + CBs) / (BA + 0.4 * ET), d / (BA + 0.4 * (ET - GE) + 0.656 * GE), 0.656 * GE / (BA + 0.4 * (ET - GE) + 0.656 * GE)))
LOG.close()
