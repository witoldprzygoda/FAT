#!/usr/bin/env python3
# ============================================================================
# vedgeu_angle.py  --  INDEPENDENT VERIFIER of the package edgeu: is the C2 gamma* angle factor the same at the u edge?
# ============================================================================
# Written 2026-11-15 (prefix vedgeu). edgeu applies the WP8 all-M4 factors (E1 0.879, E2 0.8825, means of the WP8 toy and
# closure) to C2 at the u edge without re-deriving them for u. This toy asks whether the ratio of the step heights
# with the QED gamma* -> e+e- angle (1 + cos^2 + (1 - beta^2) sin^2) and with the isotropic HGeant angle is the same
# at the u edge as over all M4.
# Generator: the WP8 toy event generator (eta (pT, y) from SMASH, KLOE-2 Dalitz plot, Kroll-Wada q, isotropic and QED
# weights on the same events), imported from bedge_angle_toy.py (functions load_gen, generate, accept; own random
# stream). Everything else is own code: u_true = 0.99274 [M4 - M4max(M(pipi), q)] with the TRUE q, smeared with
# (a) a Gaussian of 6.70 MeV, (b) the simulated u residuals of C2 + C3 near the edge (vedgeu columns); the step
# height is taken with the linear estimator of the nominal data fit (vedgeu_core, the way h_sim is taken).
# Run: python3 -I plots/all_report/vedgeu_angle.py > /hdd3/przygoda/HADES/FAT_out/all_report/vedgeu/vedgeu_angle_stdout.txt
# ============================================================================
import json
import math
import sys

import numpy as np

sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import vedgeu_core as V
import bedge_angle_toy as G

G.rng = np.random.default_rng(0xA9E1)       # own random stream for the generator
rng = np.random.default_rng(0xA9E2)

# nominal data fit (linear estimator Lh), as in vedgeu_ana.py
D, S = V.load()
ud = V.u_data(D); bD = V.ubin(ud); sD = V.slice_index(D["pp"]); s0d = (D["bits"] & 1) != 0
Hk = [V.hist(bD[s0d & (D["k"] == k)], sD[s0d & (D["k"] == k)]) for k in range(3)]
N, CB, CBC = V.spectrum(*Hk)
us, mpp = V.u_sim(S); bS = V.ubin(us); sS = V.slice_index(mpp); s0s = (S["bits"] & 1) != 0
w = np.where(s0s, S["w"], 0.0)
K = V.hist(bS[np.isin(S["c"], [2, 3])], sS[np.isin(S["c"], [2, 3])], w[np.isin(S["c"], [2, 3])]).sum(0)
R0 = V.gls_fit(N, CB, CBC, V.step(0.0067), fixed=V.AP * K)
# empirical residuals (sim C2 + C3 near the edge)
mt4 = V.m4max(S["ppt"].astype(float), S["eet"].astype(float)); ut = V.ESC * (S["m4t"].astype(float) - mt4)
mr = (S["c"] <= 1) & s0s & np.isfinite(us) & np.isfinite(ut) & (mpp >= 0.30) & (mpp < 0.40) & (ut > -0.015) & (S["m4t"] > 0)
RES, RESW = (us - ut)[mr], w[mr]
RESP = RESW / RESW.sum()

h = G.load_gen()
acc_h = {}
CL = [("all", 0.0, 1.0), ("E1", 0.0, 0.010), ("E2", 0.010, 1.0)]
SM = ["gauss", "emp"]
for ic in range(4):
    ev = G.generate(3_000_000, h)
    ok = G.accept(ev, 0.06)
    mppd = 0.9970 * ev["mpp"]
    q = ev["q"]
    u_t = V.ESC * (ev["m4"] - V.m4max(ev["mpp"], q))
    insl = (mppd >= 0.30) & (mppd < 0.40)
    for sm in SM:
        if sm == "gauss":
            uu = u_t + rng.normal(0, 0.0067, u_t.size)
        else:
            uu = u_t + RES[rng.choice(RES.size, size=u_t.size, p=RESP)]
        b = V.ubin(uu)
        for cn, lo, hi in CL:
            m = ok & insl & (q >= lo) & (q < hi) & (b >= 0)
            for lab, ww in (("iso", ev["wk"]), ("qed", ev["wk"] * ev["wq"]), ("v", ev["wk"] ** 2)):
                key = (sm, cn, lab)
                acc_h[key] = acc_h.get(key, 0) + np.bincount(b[m], weights=ww[m], minlength=V.NB)
            # all M4, all M(pipi) acceptance ratio
            m2 = ok & (q >= lo) & (q < hi)
            key = (sm, cn, "all")
            a = acc_h.get(key, np.zeros(2)); a += [ev["wk"][m2].sum(), (ev["wk"] * ev["wq"])[m2].sum()]; acc_h[key] = a
    print("chunk %d" % ic, flush=True)
OUT = {}
for sm in SM:
    for cn, lo, hi in CL:
        yi, yq, yv = acc_h[(sm, cn, "iso")], acc_h[(sm, cn, "qed")], acc_h[(sm, cn, "v")]
        Ru = V.proj(R0, yq) / V.proj(R0, yi)
        # own fit with Poisson-like weights of the toy, same template
        Rf = V.gls_fit(yi, np.zeros(V.NB), np.zeros((V.NB, V.NB)), V.step(0.0067), Vfix=np.diag(np.maximum(yv[R0.idx], 1e-9)))
        Rq = Rf.L[-1] @ yq[Rf.idx] / Rf.h
        ra = acc_h[(sm, cn, "all")]
        OUT["%s_%s" % (sm, cn)] = dict(R_u_datacov=Ru, R_u_own=float(Rq), R_all=float(ra[1] / ra[0]))
        print("%-6s %-4s R_u (data estimator) %.4f; R_u (own weights) %.4f; R_all (all M4, all M(pipi)) %.4f" % (sm, cn, Ru, Rq, ra[1] / ra[0]))
# effective factor applied by edgeu on the C2 step: weights of E1 / E2 in the C2 step projection
c2 = S["c"] == 0
fac = np.where(S["ee"] < 0.010, V.CF1, V.CF2)
h1 = V.proj(R0, V.hist(bS[c2 & (S["ee"] < 0.010)], sS[c2 & (S["ee"] < 0.010)], w[c2 & (S["ee"] < 0.010)]).sum(0))
h2 = V.proj(R0, V.hist(bS[c2 & (S["ee"] >= 0.010)], sS[c2 & (S["ee"] >= 0.010)], w[c2 & (S["ee"] >= 0.010)]).sum(0))
OUT["sim_C2_step_E1_E2"] = [h1, h2]
print("SMASH C2 step projection: E1 %.1f, E2 %.1f (E2 share %.3f); applied effective factor %.4f" % (h1, h2, h2 / (h1 + h2), (V.CF1 * h1 + V.CF2 * h2) / (h1 + h2)))
json.dump(OUT, open(V.DIR + "vedgeu_angle.json", "w"), indent=1)
