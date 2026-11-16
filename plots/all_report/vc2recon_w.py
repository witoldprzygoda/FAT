#!/usr/bin/env python3
# ============================================================================
# vc2recon_w.py -- independent verification of package c2recon, part 1: per-row QED weights of the simulated C2 rows
# ============================================================================
# Written 2026-11-15 (prefix vc2recon, verifier). Own code; no c2recon macro, weight or output is read.
# Input : /hdd3/.../vc2recon/vc2recon_rows.txt (vc2recon_dump.C: PAT truth of the 1649 C2 rows, edgeu order)
#         SMASH eta lists: C = own extraction (files 3001-3040, vc2_smash_eta_C.txt); A (ecalacc) and B (vecalacc) as
#         variants only.
# Method (own derivation, LAB frame):
#   photon k = |k| (1, n): 2 Q.k = A = m_pi0^2 - Q^2, 2 P.k = B = m_eta^2 - P^2 (Q = e+e-, P = pi+pi-e+e-)
#   => (A P_vec - B Q_vec).n = A P0 - B Q0 (n on a cone), |k| = A / (2 (Q0 - Q_vec.n)),
#   measure d^3k/|k| delta delta = |k(phi)| dphi / (2 |V|) along the cone azimuth phi (checked numerically against a flat
#   azimuth in the e+e- rest frame);
#   prior = SMASH invariant eta density E d^3N/dp^3 = dN/(dy d pT^2)/pi at P + k (HGeant: eta -> 3 pi phase space,
#   pi0 -> gamma gamma* isotropic, so no further factor);
#   weight w = 3 [1 + c^2 + e (1 - c^2)] / (4 + 2 e), e = 4 m_e^2 / Q^2, c = cos(e+, k) in the Q rest frame
#   (= e+ against the gamma* direction of the pi0 rest frame).
# Output: /hdd3/.../vc2recon/vc2_w.npz (per-row weights and kinematics), vc2_w.json, vc2_w.log
# Run: python3 -I plots/all_report/vc2recon_w.py
# ============================================================================
import json
import math
import sys
import time

import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.interpolate import RegularGridInterpolator

O = "/hdd3/przygoda/HADES/FAT_out/all_report/vc2recon/"
ETA_C = O + "vc2_smash_eta_C.txt"
ETA_A = "/hdd3/przygoda/HADES/FAT_out/all_report/ecalacc/smash_eta_lab.txt"
ETA_B = "/hdd3/przygoda/HADES/FAT_out/all_report/vecalacc/vsmash_eta_lab.txt"
MEL, MPI = 0.51099895, 139.57039      # MeV
NPHI = 720
T0 = time.time()
LOG = None
def P(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    if LOG is not None:
        LOG.write(s + "\n"); LOG.flush()
RES = {}

# ------------------------------------------------------------------ rows
def parse(fn):
    rows = []
    for l in open(fn):
        if l.startswith("#"):
            continue
        segs = [s.split() for s in l.split("|")]
        rows.append(segs)
    n = len(rows)
    a0 = np.array([[float(x) for x in r[0]] for r in rows])
    a1 = np.array([[float(x) for x in r[1]] for r in rows])
    a3 = np.array([[float(x) for x in r[3]] for r in rows])
    a4 = np.array([[float(x) for x in r[4]] for r in rows])
    a5 = np.array([[float(x) for x in r[5]] for r in rows])
    names0 = "seq tree idx set pat_entry grp bits isBest ownA ownB cls lcls pcls tmpl W0 f108 m4 m4t mpipi mpipit mee meet ep_sim_p pip_sim_p".split()
    R = {k: a0[:, i] for i, k in enumerate(names0)}
    for k in ("seq", "tree", "idx", "set", "pat_entry", "grp", "bits", "isBest", "ownA", "ownB", "cls", "lcls", "pcls", "tmpl"):
        R[k] = R[k].astype(np.int64)
    R["mom"] = a1[:, :12].reshape(n, 4, 3)       # ep em pip pim
    R["par"] = a1[:, 12:14]
    R["g_idx"], R["g_e"], R["m_eeg"], R["m5"], R["mm4"], R["vz"] = [a3[:, i] for i in range(6)]
    R["nmult"] = a4[:, 0]
    cl = a4[:, 1:].reshape(n, 5, 9)
    R["cl_id"], R["cl_ppid"] = cl[:, :, 0], cl[:, :, 1]
    R["cl_p"] = cl[:, :, 2:5]; R["cl_e"] = cl[:, :, 5]; R["cl_v"] = cl[:, :, 6:9]
    R["pv"] = a5
    R["n"] = n
    return R

def four(p3, m):
    return np.concatenate([np.sqrt((p3 * p3).sum(-1, keepdims=True) + m * m), p3], -1)
def mdot(a, b):
    return a[..., 0] * b[..., 0] - (a[..., 1:] * b[..., 1:]).sum(-1)
def boost_to_rest(p, R):
    """p (..., 4) boosted into the rest frame of R (..., 4) (broadcast)"""
    b = R[..., 1:] / R[..., :1]
    b2 = (b * b).sum(-1, keepdims=True)
    g = 1.0 / np.sqrt(1.0 - b2)
    bp = (b * p[..., 1:]).sum(-1, keepdims=True)
    gm1 = np.where(b2 > 0, (g - 1.0) / np.where(b2 > 0, b2, 1.0), 0.0)
    E = g * (p[..., :1] - bp)
    v = p[..., 1:] + gm1 * bp * b - g * b * p[..., :1]
    return np.concatenate([E, v], -1)
def boost_from_rest(p, R):
    """inverse: p given in the rest frame of R -> frame where R has its momentum"""
    b = -R[..., 1:] / R[..., :1]
    Rm = np.concatenate([R[..., :1], -R[..., 1:]], -1)
    return boost_to_rest(p, Rm)

def qed_w(c, eps):
    return 3.0 * (1.0 + c * c + eps * (1.0 - c * c)) / (4.0 + 2.0 * eps)

def costh(ep, k, Q):
    """cos(e+, k) in the Q rest frame"""
    a = boost_to_rest(ep, Q)[..., 1:]; b = boost_to_rest(k, Q)[..., 1:]
    return (a * b).sum(-1) / np.sqrt((a * a).sum(-1) * (b * b).sum(-1))

# ------------------------------------------------------------------ SMASH density
class Dens:
    """invariant density E d3N/dp3 ~ dN/(dy d pT^2) (GeV units), histogram + smoothing + bilinear interpolation"""
    def __init__(self, fn, ny=56, npt=56, ylo=-0.6, yhi=3.0, p2hi=1.5, smooth=1.2, flat=False):
        a = np.loadtxt(fn, usecols=(0, 1, 2, 3))
        E, px, py, pz = a.T
        m2 = E * E - px * px - py * py - pz * pz
        ok = (m2 > 0) & (E > np.abs(pz))
        E, px, py, pz = E[ok], px[ok], py[ok], pz[ok]
        self.n = int(ok.sum()); self.mass = float(np.sqrt(m2[ok]).mean())
        y = 0.5 * np.log((E + pz) / (E - pz)); p2 = px * px + py * py
        self.yb = np.linspace(ylo, yhi, ny + 1); self.pb = np.linspace(0.0, p2hi, npt + 1)
        H, _, _ = np.histogram2d(y, p2, bins=[self.yb, self.pb])
        self.outside = float(1 - H.sum() / len(y))
        if smooth > 0:
            H = gaussian_filter(H, sigma=smooth, mode="constant")
        yc = 0.5 * (self.yb[1:] + self.yb[:-1]); pc = 0.5 * (self.pb[1:] + self.pb[:-1])
        self.yc, self.pc = yc, pc
        self.f = RegularGridInterpolator((yc, pc), H, bounds_error=False, fill_value=0.0)
        self.flat = flat
        self.ylo, self.yhi, self.p2hi = yc[0], yc[-1], pc[-1]
    def __call__(self, p):   # p (..., 4) MeV
        if self.flat:
            return np.ones(p.shape[:-1])
        E, px, py, pz = (p[..., i] / 1e3 for i in range(4))
        with np.errstate(divide="ignore", invalid="ignore"):
            y = 0.5 * np.log((E + pz) / (E - pz))
        p2 = px * px + py * py
        p2c = np.clip(p2, self.pc[0], None)
        v = self.f(np.stack([np.clip(np.nan_to_num(y, nan=-9), -9, 9), p2c], -1))
        v[(y < self.ylo) | (y > self.yhi) | ~np.isfinite(y)] = 0.0
        return np.maximum(v, 0.0)

# ------------------------------------------------------------------ the cone in the lab
def cone_lab(Q, Pv, meta, mpi0, nphi=NPHI):
    """Q, Pv (n,4) -> k (n, nphi, 4), measure (n, nphi), ok (n)"""
    A = mpi0 ** 2 - mdot(Q, Q); B = meta ** 2 - mdot(Pv, Pv)
    V = A[:, None] * Pv[:, 1:] - B[:, None] * Q[:, 1:]
    Vn = np.linalg.norm(V, axis=1)
    ca = (A * Pv[:, 0] - B * Q[:, 0]) / Vn
    ok = (np.abs(ca) <= 1.0) & (A > 0) & (B > 0)
    ca = np.clip(ca, -1, 1); sa = np.sqrt(1 - ca * ca)
    ah = V / Vn[:, None]
    t = np.where((np.abs(ah[:, 2]) < 0.9)[:, None], np.array([0, 0, 1.0])[None, :], np.array([1.0, 0, 0])[None, :])
    e1 = np.cross(t, ah); e1 /= np.linalg.norm(e1, axis=1)[:, None]; e2 = np.cross(ah, e1)
    phi = 2 * np.pi * (np.arange(nphi) + 0.5) / nphi
    n = ca[:, None, None] * ah[:, None, :] + sa[:, None, None] * (np.cos(phi)[None, :, None] * e1[:, None, :] + np.sin(phi)[None, :, None] * e2[:, None, :])
    kk = A[:, None] / (2 * (Q[:, None, 0] - (n * Q[:, None, 1:]).sum(-1)))
    k = np.concatenate([kk[..., None], kk[..., None] * n], -1)
    return k, kk, ok

def cone_qframe(Q, Pv, meta, mpi0, nphi=NPHI):
    """alternative: flat azimuth around P in the Q rest frame; returns k in the lab (n, nphi, 4)"""
    A = mpi0 ** 2 - mdot(Q, Q); B = meta ** 2 - mdot(Pv, Pv)
    MQ = np.sqrt(mdot(Q, Q))
    Pq = boost_to_rest(Pv, Q)
    kq = A / (2 * MQ)
    pn = np.linalg.norm(Pq[:, 1:], axis=1)
    ct = (Pq[:, 0] - B / (2 * kq)) / pn
    ok = np.abs(ct) <= 1
    ct = np.clip(ct, -1, 1); st = np.sqrt(1 - ct * ct)
    ah = Pq[:, 1:] / pn[:, None]
    t = np.where((np.abs(ah[:, 2]) < 0.9)[:, None], np.array([0, 0, 1.0])[None, :], np.array([1.0, 0, 0])[None, :])
    e1 = np.cross(t, ah); e1 /= np.linalg.norm(e1, axis=1)[:, None]; e2 = np.cross(ah, e1)
    phi = 2 * np.pi * (np.arange(nphi) + 0.5) / nphi
    n = ct[:, None, None] * ah[:, None, :] + st[:, None, None] * (np.cos(phi)[None, :, None] * e1[:, None, :] + np.sin(phi)[None, :, None] * e2[:, None, :])
    k4 = np.concatenate([np.broadcast_to(kq[:, None, None], n.shape[:2] + (1,)), kq[:, None, None] * n], -1)
    klab = boost_from_rest(k4, np.broadcast_to(Q[:, None, :], k4.shape))
    return klab, ok

def circle_mean(Q, Pv, ep, dens, meta, mpi0, frame="lab", extra=None, nphi=NPHI, chunk=200, noeps=False):
    """posterior mean of w over the photon circle; extra(k) -> additional factor (e.g. an acceptance)"""
    n = len(Q)
    out = np.full(n, np.nan); okall = np.zeros(n, bool); fz = np.zeros(n, bool)
    eps = 4 * MEL ** 2 / mdot(Q, Q)
    if noeps:
        eps = np.zeros(n)
    for s in range(0, n, chunk):
        sl = slice(s, min(n, s + chunk))
        if frame == "lab":
            k, meas, ok = cone_lab(Q[sl], Pv[sl], meta, mpi0, nphi)
        else:
            k, ok = cone_qframe(Q[sl], Pv[sl], meta, mpi0, nphi)
            meas = np.ones(k.shape[:2])
        Qb = np.broadcast_to(Q[sl][:, None, :], k.shape); eb = np.broadcast_to(ep[sl][:, None, :], k.shape)
        c = costh(eb, k, Qb)
        w = qed_w(c, eps[sl][:, None])
        pe = np.broadcast_to(Pv[sl][:, None, :], k.shape) + k
        f = dens(pe) * meas
        if extra is not None:
            f = f * extra(k, sl)
        sf = f.sum(1)
        out[sl] = np.where(sf > 0, (f * w).sum(1) / np.where(sf > 0, sf, 1), np.nan)
        okall[sl] = ok; fz[sl] = sf <= 0
    return out, okall, fz

# ======================================================================================================================
if __name__ == "__main__":
    LOG = open(O + "vc2_w.log", "w")
    R = parse(O + "vc2recon_rows.txt")
    n = R["n"]
    P("# vc2recon_w.py, 2026-11-15: %d C2 rows (%.0f s)" % (n, time.time() - T0))
    m = R["mom"]
    ep, em, pip, pim = four(m[:, 0], MEL), four(m[:, 1], MEL), four(m[:, 2], MPI), four(m[:, 3], MPI)
    Q = ep + em; Pv = Q + pip + pim; PP = pip + pim
    MQ = np.sqrt(mdot(Q, Q)); MP = np.sqrt(mdot(Pv, Pv)); MPP = np.sqrt(mdot(PP, PP))
    d1 = np.abs(MQ - 1e3 * R["meet"]).max(); d2 = np.abs(MP - 1e3 * R["m4t"]).max(); d3 = np.abs(MPP - 1e3 * R["mpipit"]).max()
    P("A. truth check: max |M(ee) - mee_true| %.3g MeV, |M4 - m4_true| %.3g MeV, |M(pipi) - mpipi_true| %.3g MeV" % (d1, d2, d3))
    P("   lepton parents (GEANT id) e+: %s, e-: %s" % (dict(zip(*np.unique(R["par"][:, 0], return_counts=True))), dict(zip(*np.unique(R["par"][:, 1], return_counts=True)))))
    RES["truth_check"] = dict(dmee=float(d1), dm4=float(d2), dmpp=float(d3))

    # ---------------------------------------------------------------- B. photons in the ECAL clusters, HGeant masses
    cid = R["cl_id"]; cp = R["cl_p"]
    nm = np.clip(np.rint(R["nmult"]).astype(int), 0, 5)
    valid = np.arange(5)[None, :] < nm[:, None]
    kc = four(cp, 0.0)
    m0c = np.sqrt(np.maximum(mdot(Q[:, None, :] + kc, Q[:, None, :] + kc), 0))
    mec = np.sqrt(np.maximum(mdot(Pv[:, None, :] + kc, Pv[:, None, :] + kc), 0))
    isg = valid & (np.rint(cid) == 1) & (np.linalg.norm(cp, axis=-1) > 0)
    # the true masses of HGeant: the sharp peaks of M(e+e-gamma) and M(pi+pi-e+e-gamma) for photon clusters
    x0 = m0c[isg]; x1 = mec[isg]
    sel0 = np.abs(x0 - 135.0) < 2.0; sel1 = np.abs(x1 - 547.7) < 2.0
    mpi0_hg = float(np.median(x0[sel0 & sel1])); meta_hg = float(np.median(x1[sel0 & sel1]))
    P("B. photon clusters (GEANT id 1): %d in %d rows; M(e+e-gamma) and M5 within 2 MeV of 135.0 / 547.7: %d" % (isg.sum(), isg.any(1).sum(), (sel0 & sel1).sum()))
    P("   HGeant masses from these: m_pi0 = %.4f MeV (spread %.2g), m_eta = %.4f MeV (spread %.2g)  [PDG 134.9768, 547.862; GEANT3 0.13498, 0.54745]" %
      (mpi0_hg, np.std(x0[sel0 & sel1]), meta_hg, np.std(x1[sel0 & sel1])))
    tol = 0.05
    match = isg & (np.abs(m0c - mpi0_hg) < tol) & (np.abs(mec - meta_hg) < tol)
    found = match.any(1)
    kfound = np.zeros((n, 4))
    for i in np.where(found)[0]:
        j = np.where(match[i])[0][0]
        kfound[i] = kc[i, j]
    # how many photon clusters near but not inside the tolerance (other photons of the event)?
    near = isg & ~match & (np.abs(m0c - mpi0_hg) < 2.0)
    P("   own match (photon cluster with M(e+e-gamma) = m_pi0 and M5 = m_eta within %.2f MeV): %d rows; parent pid 7 among them: %.3f; near-misses (< 2 MeV in M(ee gamma) only): %d" %
      (tol, found.sum(), float(np.mean(np.rint(R["cl_ppid"][match]) == 7)), near.sum()))
    RES["masses"] = dict(mpi0=mpi0_hg, meta=meta_hg, nfound=int(found.sum()))
    # several clusters of the same photon?
    P("   rows with 2+ matching clusters (split photon): %d" % (match.sum(1) >= 2).sum())

    # ---------------------------------------------------------------- C. exact weights of the found photons; theta* definitions
    eps = 4 * MEL ** 2 / MQ ** 2
    cex = costh(ep, kfound, Q); w_ex = np.where(found, qed_w(cex, eps), np.nan)
    pi0 = Q + kfound
    # (ii) e+ in the gamma* rest frame against the gamma* direction in the pi0 rest frame (two-step boost)
    Qp = boost_to_rest(Q, pi0); epp_ = boost_to_rest(ep, pi0)
    eq = boost_to_rest(epp_, Qp)[:, 1:]; ax = Qp[:, 1:]
    c2 = (eq * ax).sum(1) / np.linalg.norm(eq, axis=1) / np.linalg.norm(ax, axis=1)
    # (iii) wrong axis: the lab direction of Q
    eql = boost_to_rest(ep, Q)[:, 1:]
    c3 = (eql * Q[:, 1:]).sum(1) / np.linalg.norm(eql, axis=1) / np.linalg.norm(Q[:, 1:], axis=1)
    # (iv) e- instead of e+
    c4 = costh(em, kfound, Q)
    f_ = found
    P("C. found photons: %d rows. theta* definitions: max ||c(e+,k)_Q| - |c(e+, gamma* axis of the pi0 frame)|| = %.2e; e- instead of e+: %.2e;" %
      (f_.sum(), np.nanmax(np.abs(np.abs(cex[f_]) - np.abs(c2[f_]))), np.nanmax(np.abs(np.abs(cex[f_]) - np.abs(c4[f_])))))
    W = R["W0"] * R["f108"]
    bits = R["bits"]
    s0A = ((bits & 1) != 0) & ((bits & 8) != 0)
    preA = ((bits & 2) != 0) & ((bits & 8) != 0)
    P("   <w> exact (unweighted) %.4f; the wrong axis (lab direction of Q) would give %.4f" % (np.nanmean(w_ex[f_]), np.mean(qed_w(c3[f_], eps[f_]))))

    # ---------------------------------------------------------------- D. density, circle, measure check
    dens = {}
    for tag, fn, kw in [("C", ETA_C, {}), ("A", ETA_A, {}), ("B", ETA_B, {}), ("C_coarse", ETA_C, dict(ny=28, npt=28, smooth=0.8)),
                        ("C_raw", ETA_C, dict(smooth=0.0)), ("flat", ETA_C, dict(flat=True))]:
        try:
            dens[tag] = Dens(fn, **kw)
            P("D. SMASH density %s: %d etas, mean mass %.5f GeV, outside the grid %.4f" % (tag, dens[tag].n, dens[tag].mass, dens[tag].outside))
        except Exception as ex:
            P("D. SMASH density %s not available: %s" % (tag, ex))
    dC = dens["C"]
    # measure check: lab cone with |k| dphi against a flat azimuth in the Q rest frame
    wl, okl, fzl = circle_mean(Q, Pv, ep, dC, meta_hg, mpi0_hg, "lab")
    wq, okq, fzq = circle_mean(Q, Pv, ep, dC, meta_hg, mpi0_hg, "q")
    P("   rows with a solution: lab %d, Q frame %d of %d; zero density on the whole circle: %d" % (okl.sum(), okq.sum(), n, fzl.sum()))
    P("   measure check: max |<w>_lab(|k| dphi) - <w>_Qframe(flat)| = %.2e (rms %.2e)" % (np.nanmax(np.abs(wl - wq)), np.sqrt(np.nanmean((wl - wq) ** 2))))
    # without the |k| factor (a wrong lab measure) for comparison
    wl_nok, _, _ = circle_mean(Q, Pv, ep, dC, meta_hg, mpi0_hg, "lab", extra=lambda k, sl: 1.0 / k[..., 0])
    P("   (a flat lab azimuth, i.e. a wrong measure, would shift the mean by %+.4f)" % np.nanmean(wl_nok - wl))
    # does the true photon lie on the circle? (found rows: distance of k_true to the nearest circle point)
    kl, _, _ = cone_lab(Q[f_], Pv[f_], meta_hg, mpi0_hg)
    dmin = np.min(np.linalg.norm(kl[..., 1:] - kfound[f_][:, None, 1:], axis=-1), axis=1)
    P("   found photons: distance to the nearest of the %d circle points: median %.3f MeV, max %.3f MeV (|k| median %.1f MeV)" %
      (NPHI, np.median(dmin), dmin.max(), np.median(kfound[f_, 0])))
    wv = {}
    for tag, d in dens.items():
        if tag == "C":
            wv[tag] = wl
            continue
        wv[tag], _, _ = circle_mean(Q, Pv, ep, d, meta_hg, mpi0_hg, "lab")
    # PDG masses instead of the HGeant ones
    wv["pdgmass"], okp, _ = circle_mean(Q, Pv, ep, dC, 547.862, 134.9768, "lab")
    # KLOE-2 X dependence along the circle (eta -> 3 pi matrix element; Y fixed by M(pipi))
    def kloe(k, sl):
        pe = Pv[sl][:, None, :] + k
        pi0_ = Q[sl][:, None, :] + k
        Tp = boost_to_rest(np.broadcast_to(pip[sl][:, None, :], k.shape), pe)[..., 0] - MPI
        Tm = boost_to_rest(np.broadcast_to(pim[sl][:, None, :], k.shape), pe)[..., 0] - MPI
        T0_ = boost_to_rest(pi0_, pe)[..., 0] - mpi0_hg
        Qe = meta_hg - 2 * MPI - mpi0_hg
        X = math.sqrt(3) * (Tp - Tm) / Qe; Y = 3 * T0_ / Qe - 1
        a, b, d_, f = -1.095, 0.145, 0.081, 0.141
        return 1 + a * Y + b * Y * Y + d_ * X * X + f * Y ** 3
    wv["kloeX"], _, _ = circle_mean(Q, Pv, ep, dC, meta_hg, mpi0_hg, "lab", extra=kloe)
    wv["noeps"], _, _ = circle_mean(Q, Pv, ep, dC, meta_hg, mpi0_hg, "lab", noeps=True)
    for tag in wv:
        P("   <w> per row, density %-9s: mean over rows (W, S0 mode A) %.5f; max |diff to C| %.4f" %
          (tag, np.nansum(W * s0A * wv[tag]) / np.nansum(W * s0A * np.isfinite(wv[tag])), np.nanmax(np.abs(wv[tag] - wl))))
    # spread of w along the circle for the median row
    k, meas, ok = cone_lab(Q, Pv, meta_hg, mpi0_hg)
    cc = costh(np.broadcast_to(ep[:, None, :], k.shape), k, np.broadcast_to(Q[:, None, :], k.shape))
    wc = qed_w(cc, eps[:, None])
    fpe = dC(np.broadcast_to(Pv[:, None, :], k.shape) + k) * meas
    rng = wc.max(1) - wc.min(1)
    # rms of w along the circle (posterior)
    sf = fpe.sum(1); mw = (fpe * wc).sum(1) / np.where(sf > 0, sf, 1)
    sw = np.sqrt(np.maximum((fpe * wc * wc).sum(1) / np.where(sf > 0, sf, 1) - mw ** 2, 0))
    P("   along the circle: range of w median %.3f (min..max of the median row %.3f..%.3f); posterior rms of w median %.4f" %
      (np.median(rng), np.median(wc.min(1)), np.median(wc.max(1)), np.median(sw)))
    RES["post_rms_median"] = float(np.median(sw))

    # ---------------------------------------------------------------- E. found photons: exact against the circle
    def wmean(x, ww, m):
        return float(np.nansum((ww * x)[m]) / np.nansum(ww[m]))
    for lab, msk in [("all rows", f_), ("pre, mode A", f_ & preA), ("S0, mode A", f_ & s0A)]:
        P("E. found photons, %-12s: %3d rows; exact <w> %.4f (W-weighted %.4f); circle (no ECAL condition) %.4f (W %.4f)" %
          (lab, msk.sum(), np.nanmean(w_ex[msk]), wmean(w_ex, W, msk), np.nanmean(wl[msk]), wmean(wl, W, msk)))
    # ECAL condition approximated by the polar angle of the photon (lab, from the target) and the six sectors:
    th = np.degrees(np.arccos(kfound[f_, 3] / kfound[f_, 0]))
    ph = np.degrees(np.arctan2(kfound[f_, 2], kfound[f_, 1])) % 60.0
    thlo, thhi = np.percentile(th, [0.5, 99.5])
    P("   found photons: lab polar angle 0.5-99.5 %% range %.1f-%.1f deg; azimuth within a sector (mod 60): %.1f-%.1f deg" %
      (thlo, thhi, *np.percentile(ph, [0.5, 99.5])))
    def ecal_acc(k, sl, lo=thlo, hi=thhi, pmin=0.0):
        t = np.degrees(np.arccos(np.clip(k[..., 3] / k[..., 0], -1, 1)))
        return ((t >= lo) & (t <= hi) & (k[..., 0] > pmin)).astype(float)
    emin = np.percentile(kfound[f_, 0], 0.5)
    wacc, _, _ = circle_mean(Q, Pv, ep, dC, meta_hg, mpi0_hg, "lab", extra=lambda k, sl: ecal_acc(k, sl, pmin=emin))
    dd = (w_ex - wacc)[f_]
    P("   circle restricted to the polar range %.1f-%.1f deg and E_gamma > %.0f MeV: %.4f against exact %.4f: difference exact - circle %+.4f +- %.4f (corr %.2f)" %
      (thlo, thhi, emin, np.nanmean(wacc[f_]), np.nanmean(w_ex[f_]), np.nanmean(dd), np.nanstd(dd) / math.sqrt(np.isfinite(dd).sum()),
       np.corrcoef(w_ex[f_][np.isfinite(wacc[f_])], wacc[f_][np.isfinite(wacc[f_])])[0, 1]))
    RES["found"] = dict(n=int(f_.sum()), exact=float(np.nanmean(w_ex[f_])), circ_all=float(np.nanmean(wl[f_])), circ_acc=float(np.nanmean(wacc[f_])),
                        diff=float(np.nanmean(dd)), ediff=float(np.nanstd(dd) / math.sqrt(np.isfinite(dd).sum())))
    # posterior calibration (PIT) for the found photons with the polar acceptance: is the true w a typical draw?
    kfa, measa, _ = cone_lab(Q[f_], Pv[f_], meta_hg, mpi0_hg)
    cca = costh(np.broadcast_to(ep[f_][:, None, :], kfa.shape), kfa, np.broadcast_to(Q[f_][:, None, :], kfa.shape))
    wca = qed_w(cca, eps[f_][:, None])
    fa = dC(np.broadcast_to(Pv[f_][:, None, :], kfa.shape) + kfa) * measa * ecal_acc(kfa, None, pmin=emin)
    sfa = fa.sum(1)
    pit = np.array([(fa[i] * (wca[i] < w_ex[f_][i])).sum() / sfa[i] if sfa[i] > 0 else np.nan for i in range(f_.sum())])
    hp, _ = np.histogram(pit[np.isfinite(pit)], bins=5, range=(0, 1))
    P("   PIT of the exact w in its restricted posterior (5 bins, uniform if calibrated): %s; mean %.3f (0.5 expected, +- %.3f)" %
      (hp.tolist(), np.nanmean(pit), 0.2887 / math.sqrt(np.isfinite(pit).sum())))
    RES["pit"] = dict(hist=hp.tolist(), mean=float(np.nanmean(pit)))

    # ---------------------------------------------------------------- F. the per-row weight: exact where known
    w_best = np.where(found, w_ex, wl)
    # soft-lepton dependence
    pe_ = np.linalg.norm(m[:, 0], axis=1); pm_ = np.linalg.norm(m[:, 1], axis=1)
    psoft = np.minimum(pe_, pm_)
    for a, b in [(0, 100), (100, 150), (150, 200), (200, 300), (300, 2000)]:
        mm = s0A & (psoft >= a) & (psoft < b)
        P("F. softer lepton %4d-%4d MeV/c: %4d rows (S0 A), <w>_circle (W) %.4f" % (a, b, mm.sum(), wmean(wl, W, mm)))
    for lab, mm in [("all S0 A", s0A), ("E1", s0A & (R["mee"] < 0.010)), ("E2", s0A & (R["mee"] >= 0.010))]:
        P("   %-8s: window mean (all M4) %.4f (exact where found: %.4f)" % (lab, wmean(wl, W, mm), wmean(w_best, W, mm)))
    # ---------------------------------------------------------------- G. ECAL-tagged C2 rows (frozen tag of ecalacc)
    kM = 0.1349768
    a_ = R["m_eeg"] ** 2 - R["mee"] ** 2
    fq = (kM * kM - R["mee"] ** 2) / np.where(a_ > 0, a_, 1)
    m5c = np.sqrt(np.maximum(R["m4"] ** 2 + fq * (R["m5"] ** 2 - R["m4"] ** 2), 0))
    tag = (R["g_idx"] > 0) & (R["mee"] < 0.100) & (R["m_eeg"] >= 0.1190) & (R["m_eeg"] < 0.1730) & (a_ > 0)
    tw = tag & preA & (m5c >= 0.49) & (m5c < 0.62)
    # is the tagged photon the true one? (a matched cluster with the tagged energy)
    tagtrue = np.array([bool(found[i]) and np.any(match[i] & (np.abs(R["cl_e"][i] - R["g_e"][i]) < 1e-3 * max(1.0, R["g_e"][i]))) for i in range(n)])
    P("G. ECAL-tagged C2 rows (pre, mode A, m5c 0.49-0.62): %d; with a found true photon %d; the tagged photon is the true one %d" %
      (tw.sum(), (tw & found).sum(), (tw & tagtrue).sum()))
    c_tag = wmean(w_best, W, tw); c_tag_circ = wmean(wl, W, tw)
    rs = []
    keyg = R["tree"] * (1 << 32) + R["grp"]
    rg = np.random.default_rng(7)
    uq, inv = np.unique(keyg, return_inverse=True)
    for _ in range(2000):
        g = rg.poisson(1.0, len(uq))[inv]
        rs.append(np.nansum((g * W * w_best)[tw]) / np.nansum((g * W)[tw]))
    P("   c_tag = %.4f +- %.4f (exact where found, circle otherwise); circle only %.4f; exact only on the tagged-true rows %.4f   [c2recon/ecalacc 0.876 +- 0.009; circle 0.869]" %
      (c_tag, np.std(rs), c_tag_circ, wmean(w_ex, W, tw & tagtrue)))
    RES["tag"] = dict(n=int(tw.sum()), nfound=int((tw & found).sum()), c=c_tag, stat=float(np.std(rs)), circ=c_tag_circ)
    np.savez(O + "vc2_w.npz", seq=R["seq"], tree=R["tree"], grp=R["grp"], bits=R["bits"], W=W, m4=R["m4"], m4t=R["m4t"], mpipi=R["mpipi"],
             mee=R["mee"], meet=R["meet"], w_circ=wl, w_q=wq, w_ex=w_ex, found=found, w_best=w_best,
             **{"w_" + k: v for k, v in wv.items() if k != "C"}, psoft=psoft, w_acc=wacc, kfound=kfound,
             g_idx=R["g_idx"], g_e=R["g_e"], m_eeg=R["m_eeg"], m5=R["m5"], cl_e=R["cl_e"], match=match)
    json.dump(RES, open(O + "vc2_w.json", "w"), indent=1)
    P("done in %.0f s" % (time.time() - T0))
