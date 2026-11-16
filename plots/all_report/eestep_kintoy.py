#!/usr/bin/env python3
# eestep_kintoy.py -- package eestep: toy check that like-sign "cross-photon" lepton pairs of eta -> pi+pi-pi0 obey the bump
# limit M4 <= E(M(pipi), M(ee)) (true masses, massless leptons). Two topologies: (A) pi0 -> gamma gamma, both photons convert,
# one e- from each photon; (B) pi0 -> gamma* gamma (Dalitz), the e- of gamma* and the e- of the converted real photon.
# Conversion: the e- takes a fraction x (uniform) of the photon 4-momentum (collinear). Prints the largest u_true and the
# M(ee) distribution of the pairs near the limit. Run: python3 -I plots/all_report/eestep_kintoy.py
import numpy as np
rng = np.random.default_rng(20261115)
ME, MPI0 = 0.547862, 0.1349768
N = 400000
def iso(n):
    c = rng.uniform(-1, 1, n); p = rng.uniform(0, 2 * np.pi, n); s = np.sqrt(1 - c * c)
    return np.stack([s * np.cos(p), s * np.sin(p), c], 1)
def boost(p4, b):   # boost 4-vectors p4 (n,4) by velocity b (n,3)
    b2 = (b * b).sum(1); g = 1 / np.sqrt(1 - b2); bp = (b * p4[:, 1:]).sum(1)
    g2 = np.where(b2 > 0, (g - 1) / b2, 0.0)
    E = g * (p4[:, 0] + bp)
    v = p4[:, 1:] + (g2 * bp + g * p4[:, 0])[:, None] * b
    return np.column_stack([E, v])
def mass(p4):
    return np.sqrt(np.maximum(p4[:, 0] ** 2 - (p4[:, 1:] ** 2).sum(1), 0))
m = rng.uniform(2 * 0.13957 + 1e-4, ME - MPI0 - 1e-4, N)              # M(pipi), flat (only the limit matters)
E0 = (ME ** 2 + MPI0 ** 2 - m ** 2) / (2 * ME); p0 = np.sqrt(E0 ** 2 - MPI0 ** 2)
Epi0 = np.column_stack([E0, p0[:, None] * iso(N)])                      # pi0 in the eta frame
Ppp = np.column_stack([ME - E0, -Epi0[:, 1:]])                           # pi+pi- system
bpi0 = Epi0[:, 1:] / E0[:, None]
def lim(m, q):
    E0 = (ME ** 2 + MPI0 ** 2 - m ** 2) / (2 * ME); p0 = np.sqrt(E0 ** 2 - MPI0 ** 2)
    return np.sqrt(ME ** 2 - ME * (E0 - p0) * (1 - q ** 2 / MPI0 ** 2))
out = {}
# (A) gamma gamma, both convert
d = iso(N); k1 = np.column_stack([np.full(N, MPI0 / 2), MPI0 / 2 * d]); k2 = np.column_stack([np.full(N, MPI0 / 2), -MPI0 / 2 * d])
k1 = boost(k1, bpi0); k2 = boost(k2, bpi0)
x1 = rng.uniform(0, 1, N); x2 = rng.uniform(0, 1, N)
l1 = x1[:, None] * k1; l2 = x2[:, None] * k2
out["A"] = (l1, l2)
# (B) Dalitz gamma* (mass m*, 1/m* spectrum between 2 m_e and m_pi0) + converted real photon
mst = np.exp(rng.uniform(np.log(0.0011), np.log(MPI0 * 0.999), N))
Est = (MPI0 ** 2 + mst ** 2) / (2 * MPI0); pst = (MPI0 ** 2 - mst ** 2) / (2 * MPI0)
d = iso(N)
gst = np.column_stack([Est, pst[:, None] * d]); gam = np.column_stack([pst, -pst[:, None] * d])
dd = iso(N); lst = np.column_stack([mst / 2, mst[:, None] / 2 * dd])     # e- of gamma* in the gamma* frame (massless)
lst = boost(lst, gst[:, 1:] / gst[:, 0:1])
gam = boost(gam, bpi0); lst = boost(lst, bpi0)
x = rng.uniform(0, 1, N)
out["B"] = (lst, x[:, None] * gam)
for k, (a, b) in out.items():
    q = mass(a + b); M4 = mass(Ppp + a + b)
    u = M4 - lim(m, q)
    near = u > -0.005
    print("topology %s: max u_true = %+.2e GeV (should be <= 0); fraction within 5 MeV of the limit %.4f; M(ee) of those: median %.0f MeV, 10-90 %% %.0f-%.0f MeV" %
          (k, u.max(), near.mean(), 1e3 * np.median(q[near]), 1e3 * np.percentile(q[near], 10), 1e3 * np.percentile(q[near], 90)))
