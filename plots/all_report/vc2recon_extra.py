#!/usr/bin/env python3
# ============================================================================
# vc2recon_extra.py -- independent verification of package c2recon, part 4: model checks of the effective factor
# ============================================================================
# Written 2026-11-15 (prefix vc2recon, verifier). Uses vc2_w.npz (own weights), the edgeu columns and the frozen
# edgeu fit (as vc2recon_fit.py), and the bedge unit responses of vc2recon_bedge.C (vc2_bedge.txt).
#   1. soft-lepton efficiency: rows with a softer true lepton below 100 MeV/c x 0.8 / x 1.2;
#   2. form factor: HGeant samples M(ee) from Kroll-Wada x the VDM rho pole (uhbook.F, m_rho = 770 MeV, i.e. slope
#      m_pi0^2/m_rho^2 = 0.0307); residual weight PDG slope (a = 0.0332) / VDM, normalised over Kroll-Wada;
#   3. "exact where found" versus the circle mean: the bias of mixing (Cov(found, w)).
# Output: /hdd3/.../vc2recon/vc2_extra.json / .log
# Run: python3 -I plots/all_report/vc2recon_extra.py
# ============================================================================
import json
import sys

import numpy as np

sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import edgeu_core as C

O = "/hdd3/przygoda/HADES/FAT_out/all_report/vc2recon/"
LOG = open(O + "vc2_extra.log", "w")
def P(*a):
    s = " ".join(str(x) for x in a)
    print(s); LOG.write(s + "\n"); LOG.flush()
V = np.load(O + "vc2_w.npz"); seq = V["seq"]
MEL, MPI0, MRHO = 0.51099895, 134.9766, 770.0
RES = {}
# ---- edgeu coefficients (same nominal as vc2recon_fit.py; sigma taken from its json)
J = json.load(open(O + "vc2_fit_edgeu.json"))["edgeu"]
D = C.load_data(); S = C.load_sim()
w_nom = C.sim_weight(S, c2corr=False)
u_s0, mpp_s0, _ = C.u_sim(S)
grid = C.Grid(-0.200, 0.0025, 140, 4)
sl = C.slices_of(0.30, 0.40, [0.30, 0.34, 0.37, 0.39, 0.40])
Hd = C.hist_data(D, C.u_data(D), (D["bits"] & 1) != 0, grid)
Hs = C.hist_sim(S, u_s0, mpp_s0, w_nom, grid)
sp = C.data_spec(Hd, sl, "all", grid, 1.0, 0, 0)
R = C.fit_step(sp, C.step_gauss(grid, J["sigma_T"]), grid, -0.045, 0.045, 2, fixed=C.AP * C.sim_sum(Hs, ["C4", "C5"], sl, "all"))
P("edgeu nominal refit: h_data %.3f (vc2recon_fit %.3f)" % (R.h, J["hd"]))
ub = grid.bin(u_s0[seq]); pb = C.ppbin_of(mpp_s0[seq].astype(np.float32))
s0A = ((V["bits"] & 1) != 0) & ((V["bits"] & 8) != 0)
W = np.where(s0A, S["w"][seq] * S["f108"][seq], 0.0)
use = s0A & (pb >= C.PB(0.30)) & (pb < C.PB(0.40)) & (ub >= 0)
p = -np.ones(grid.nb, int); p[R.idx] = np.arange(len(R.idx))
a_u = np.zeros(len(seq)); q = use & (p[np.clip(ub, 0, grid.nb - 1)] >= 0)
a_u[q] = R.Lh[p[ub[q]]]
# ---- bedge coefficients
B = np.load(O + "vc2_bedge_rowbins.npz")
coef = {}; slc = {}
for l in open(O + "vc2_bedge.txt"):
    t = l.split()
    if t and t[0] == "COEF":
        k, s, pa, pb_, i, a = int(t[1]), int(t[2]), int(t[3]), int(t[4]), int(t[5]), float(t[6])
        coef[(k, s, i)] = a; slc[(k, s)] = (pa, pb_)
a_b = np.zeros(len(seq))
for j in range(len(seq)):
    if B["Wb"][j] == 0 or B["p"][j] < 0 or B["ib"][j] < 0:
        continue
    kf = B["e"][j] + 1
    for (k, s), (pa, pb_) in slc.items():
        if k == kf and pa <= B["p"][j] < pb_:
            a_b[j] = coef.get((k, s, B["ib"][j]), 0.0)
Wb = B["Wb"]
x = V["w_circ"]
def fac(a, w, m):
    return float((a * w * x * m).sum() / (a * w * m).sum())
ones = np.ones(len(seq))
cu, cb = fac(a_u, W, ones), fac(a_b, Wb, ones)
P("factors: u-edge %.5f, per-slice %.5f" % (cu, cb))
# 1. soft-lepton efficiency
soft = V["psoft"] < 100.0
for s_ in (0.8, 1.2):
    m = np.where(soft, s_, 1.0)
    P("1. soft-lepton (< 100 MeV/c) efficiency x %.1f: u-edge %.5f (%+.4f), per-slice %.5f (%+.4f)" %
      (s_, fac(a_u, W, m), fac(a_u, W, m) - cu, fac(a_b, Wb, m), fac(a_b, Wb, m) - cb))
    RES["soft_%.1f" % s_] = [fac(a_u, W, m) - cu, fac(a_b, Wb, m) - cb]
# 2. form factor residual
xx = np.linspace(2 * MEL * 1.000001, MPI0 * 0.999999, 400001)
kw = (1 - xx * xx / MPI0 ** 2) ** 3 / xx * np.sqrt(1 - 4 * MEL ** 2 / xx ** 2) * (1 + 2 * MEL ** 2 / xx ** 2)
vdm = 1.0 / (1 - xx * xx / MRHO ** 2) ** 2
for aa in (0.0332, 0.0332 + 0.0029, 0.0332 - 0.0029):
    pdg = (1 + aa * xx * xx / MPI0 ** 2) ** 2
    norm = np.trapz(kw * vdm * pdg / vdm, xx) / np.trapz(kw * vdm, xx)
    mt = 1e3 * V["meet"]
    r = ((1 + aa * mt * mt / MPI0 ** 2) ** 2 * (1 - mt * mt / MRHO ** 2) ** 2) / norm
    P("2. form factor PDG slope a = %.4f instead of the HGeant VDM pole (a_eff = %.4f): effective C2 template x %.5f (u-edge), x %.5f (per-slice)" %
      (aa, MPI0 ** 2 / MRHO ** 2, float((a_u * W * x * r).sum() / (a_u * W * x).sum()), float((a_b * Wb * x * r).sum() / (a_b * Wb * x).sum())))
    RES["ff_%.4f" % aa] = [float((a_u * W * x * r).sum() / (a_u * W * x).sum()), float((a_b * Wb * x * r).sum() / (a_b * Wb * x).sum())]
# 3. exact-where-found mixing
f_ = V["found"]
P("3. found photons: %d rows; in the u step they carry %.3f of sum a W; exact - circle on them: %+.4f (W-weighted)" %
  (f_.sum(), float((a_u * W)[f_].sum() / (a_u * W).sum()), float(((W * (V["w_ex"] - x))[f_ & s0A]).sum() / W[f_ & s0A].sum())))
P("   mixing exact-where-found with the unconditional circle mean shifts the u factor by %+.4f and the per-slice factor by %+.4f;" %
  (fac(a_u, W, ones * 0 + 1) * 0 + float((a_u * W * V["w_best"]).sum() / (a_u * W).sum()) - cu,
   float((a_b * Wb * V["w_best"]).sum() / (a_b * Wb).sum()) - cb))
P("   this mix is biased by Cov(found, w | charged tracks): the circle mean alone (c2recon's choice) is the unbiased weight.")
json.dump(RES, open(O + "vc2_extra.json", "w"), indent=1)
