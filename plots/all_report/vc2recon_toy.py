#!/usr/bin/env python3
# ============================================================================
# vc2recon_toy.py -- independent verification of package c2recon, part 2: closure toy of the circle posterior
# ============================================================================
# Written 2026-11-15 (prefix vc2recon, verifier). Own generator (no c2recon code):
#   eta: SMASH momenta of list A (ecalacc extraction), mass = the HGeant mass measured in vc2recon_w.py;
#   eta -> pi+ pi- pi0 phase space (s12 from p0 q / sqrt(s12), isotropic pi0 in the eta frame, isotropic pi+ in the dipion
#   frame); pi0 -> gamma gamma*, M(ee) from Kroll-Wada, isotropic; gamma* -> e+e- isotropic (HGeant);
#   true QED weight from the two-step definition (e+ in the gamma* rest frame against the gamma* direction of the pi0 frame);
#   acceptance on the charged tracks only (18-85 deg, pions > 100 MeV/c, leptons > p_min).
# The posterior mean <w> of each accepted decay is computed with vc2recon_w.circle_mean and the density of the OTHER
# SMASH list C (own extraction). Closure: mean(true w) = mean(<w>) in every subset defined by the charged tracks.
# Output: /hdd3/.../vc2recon/vc2_toy.json, vc2_toy.log
# Run: python3 -I plots/all_report/vc2recon_toy.py
# ============================================================================
import json
import math
import sys
import time

import numpy as np

sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import vc2recon_w as VW

O = VW.O
T0 = time.time()
TAG = sys.argv[1] if len(sys.argv) > 1 else ""
LOG = open(O + "vc2_toy%s.log" % TAG, "w")
def P(*a):
    s = " ".join(str(x) for x in a)
    print(s); LOG.write(s + "\n"); LOG.flush()
MEL, MPI = VW.MEL, VW.MPI
J = json.load(open(O + "vc2_w.json"))
META, MPI0 = J["masses"]["meta"], J["masses"]["mpi0"]
rng = np.random.default_rng(20261115 if not TAG.startswith("_seed") else 1000 + int(TAG[5:]))

def iso(n):
    c = rng.uniform(-1, 1, n); ph = rng.uniform(0, 2 * np.pi, n); s = np.sqrt(1 - c * c)
    return np.stack([s * np.cos(ph), s * np.sin(ph), c], 1)
def kw_sample(n):
    x = np.linspace(2 * MEL * 1.0000001, MPI0 * 0.9999999, 200001)
    d = (2 / x) * (1 - x * x / MPI0 ** 2) ** 3 * (1 + 2 * MEL ** 2 / x ** 2) * np.sqrt(1 - 4 * MEL ** 2 / x ** 2)
    c = np.cumsum(d); c /= c[-1]
    return np.interp(rng.uniform(0, 1, n), c, x)
def two_body(M, m1, m2, n_hat, frameP):
    """daughters of a particle (4-vector frameP in the lab) decaying with daughter-1 direction n_hat in its rest frame"""
    p = np.sqrt(np.maximum((M * M - (m1 + m2) ** 2) * (M * M - (m1 - m2) ** 2), 0)) / (2 * M)
    d1 = np.concatenate([np.sqrt(p * p + m1 * m1)[:, None], p[:, None] * n_hat], 1)
    d2 = np.concatenate([np.sqrt(p * p + m2 * m2)[:, None], -p[:, None] * n_hat], 1)
    return VW.boost_from_rest(d1, frameP), VW.boost_from_rest(d2, frameP)

def generate(n, eta_list):
    idx = rng.integers(0, len(eta_list), n)
    p3 = eta_list[idx, 1:4] * 1e3
    eta = np.concatenate([np.sqrt((p3 * p3).sum(1) + META ** 2)[:, None], p3], 1)
    # s12 = M(pi+pi-)^2 by rejection with p0 q / m12
    lo, hi = 2 * MPI, META - MPI0
    m12 = np.empty(n); todo = np.arange(n)
    gmax = None
    while len(todo):
        m = rng.uniform(lo, hi, len(todo) * 2)
        p0 = np.sqrt(np.maximum((META ** 2 - (m + MPI0) ** 2) * (META ** 2 - (m - MPI0) ** 2), 0)) / (2 * META)
        q = np.sqrt(np.maximum(m * m / 4 - MPI ** 2, 0))
        g = p0 * q / m * (2 * m)     # d(s12) = 2 m dm
        if gmax is None:
            gmax = 1.2 * g.max()
        acc = rng.uniform(0, gmax, len(m)) < g
        m = m[acc][:len(todo)]
        m12[todo[:len(m)]] = m
        todo = todo[len(m):]
    pi0, dpi = two_body(np.full(n, META), MPI0, m12, iso(n), eta)        # pi0 and the dipion
    pip, pim = two_body(m12, MPI, MPI, iso(n), dpi)
    mee = kw_sample(n)
    gst, gam = two_body(np.full(n, MPI0), mee, 0.0, iso(n), pi0)          # gamma* and the photon
    ep, em = two_body(mee, MEL, MEL, iso(n), gst)
    # true weight: two-step definition
    Qp = VW.boost_to_rest(gst, pi0); epp = VW.boost_to_rest(ep, pi0)
    eq = VW.boost_to_rest(epp, Qp)[:, 1:]; ax = Qp[:, 1:]
    c = (eq * ax).sum(1) / np.linalg.norm(eq, axis=1) / np.linalg.norm(ax, axis=1)
    eps = 4 * MEL ** 2 / mee ** 2
    w = VW.qed_w(c, eps)
    return dict(ep=ep, em=em, pip=pip, pim=pim, gam=gam, w=w, mee=mee, m12=m12)

def theta(p):
    return np.degrees(np.arccos(p[:, 3] / np.linalg.norm(p[:, 1:], axis=1)))
def pabs(p):
    return np.linalg.norm(p[:, 1:], axis=1)

if __name__ == "__main__":
    A = np.loadtxt(VW.ETA_A, usecols=(0, 1, 2, 3))
    # variants (argv[1]): "" density C smoothed (nominal); "_Asame": density from the generation list A itself;
    # "_Craw": list C without smoothing; "_Cfine": list C, finer grid and less smoothing
    if TAG == "_Asame":
        densC = VW.Dens(VW.ETA_A)
    elif TAG == "_Craw":
        densC = VW.Dens(VW.ETA_C, smooth=0.0)
    elif TAG == "_Cfine":
        densC = VW.Dens(VW.ETA_C, ny=84, npt=84, smooth=0.8)
    else:
        densC = VW.Dens(VW.ETA_C)
    densF = VW.Dens(VW.ETA_C, flat=True)
    P("# vc2recon_toy.py: generation from SMASH list A (%d etas), posterior density from list C (%d etas); m_eta %.4f, m_pi0 %.4f MeV" %
      (len(A), densC.n, META, MPI0))
    out = {}
    G = generate(3_000_000 if not TAG.startswith("_seed") else 8_000_000, A)
    P("   generated %d decays: <w> (all, no acceptance) = %.5f (1 expected), <w> with M(ee) > 10 MeV %.5f" %
      (len(G["w"]), G["w"].mean(), G["w"][G["mee"] > 10].mean()))
    th = [theta(G[k]) for k in ("ep", "em", "pip", "pim")]
    pa = [pabs(G[k]) for k in ("ep", "em", "pip", "pim")]
    geo = np.all([(t >= 18) & (t <= 85) for t in th], 0) & (pa[2] > 100) & (pa[3] > 100)
    for pmin in ((60.0, 100.0) if TAG == "" else (100.0,)):
        acc = geo & (pa[0] > pmin) & (pa[1] > pmin)
        ii = np.where(acc)[0]
        Q = G["ep"][ii] + G["em"][ii]; Pv = Q + G["pip"][ii] + G["pim"][ii]
        t1 = time.time()
        wc, ok, fz = VW.circle_mean(Q, Pv, G["ep"][ii], densC, META, MPI0, "lab", nphi=360)
        wf, _, _ = VW.circle_mean(Q, Pv, G["ep"][ii], densF, META, MPI0, "lab", nphi=360)
        wtrue = G["w"][ii]
        # subsets by charged-track quantities only
        m4 = np.sqrt(VW.mdot(Pv, Pv)); mpp = np.sqrt(VW.mdot(G["pip"][ii] + G["pim"][ii], G["pip"][ii] + G["pim"][ii]))
        mee = G["mee"][ii]
        E0 = (META ** 2 + MPI0 ** 2 - mpp ** 2) / (2 * META); p0 = np.sqrt(np.maximum(E0 * E0 - MPI0 ** 2, 0))
        m4max = np.sqrt(META ** 2 - META * (E0 - p0) * (1 - mee ** 2 / MPI0 ** 2))
        u = m4 - m4max
        ps = np.minimum(pa[0][ii], pa[1][ii])
        subsets = [("all", np.ones(len(ii), bool)), ("E1 (< 10 MeV)", mee < 10), ("E2", mee >= 10), ("near edge u > -15 MeV", u > -15),
                   ("far u < -60 MeV", u < -60), ("0.30 <= M(pipi) < 0.40", (mpp >= 300) & (mpp < 400)),
                   ("soft lepton < 100", ps < 100), ("soft lepton 100-150", (ps >= 100) & (ps < 150)), ("soft lepton 150-250", (ps >= 150) & (ps < 250)),
                   ("soft lepton > 250", ps >= 250)]
        P("\n   p_min %.0f MeV/c: %d accepted (%.2f %%), circle solved %d, zero density %d (%.0f s)" %
          (pmin, len(ii), 100 * len(ii) / len(G["w"]), ok.sum(), fz.sum(), time.time() - t1))
        res = {}
        for lab, m in subsets:
            m = m & np.isfinite(wc)
            if m.sum() < 100:
                continue
            a, b, c_ = wtrue[m].mean(), wc[m].mean(), wf[m].mean()
            d = wtrue[m] - wc[m]
            e = d.std() / math.sqrt(m.sum())
            P("   %-24s n %7d: true <w> %.4f, circle (density C) %.4f, diff %+.4f +- %.4f (%.1f sigma); flat density %.4f" %
              (lab, m.sum(), a, b, a - b, e, (a - b) / e, c_))
            res[lab] = dict(n=int(m.sum()), true=float(a), circ=float(b), flat=float(c_), err=float(e))
        out["pmin%.0f" % pmin] = res
    json.dump(out, open(O + "vc2_toy%s.json" % TAG, "w"), indent=1)
    P("done in %.0f s" % (time.time() - T0))
