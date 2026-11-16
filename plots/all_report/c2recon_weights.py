#!/usr/bin/env python3
# ============================================================================
# c2recon_weights.py  --  package c2recon: the QED gamma* weight of every simulated C2 row
# ============================================================================
# Written 2026-11-15 (prefix c2recon). New file. Reads c2recon_rows.txt (c2recon_rows.C) and writes
#   /hdd3/przygoda/HADES/FAT_out/all_report/c2recon/c2recon_weights.npz   (per-row weights, all variants)
#   /hdd3/przygoda/HADES/FAT_out/all_report/c2recon/c2recon_weights.txt   (tree idx q_nom ... for c2recon_bedge.C)
#   /hdd3/przygoda/HADES/FAT_out/all_report/c2recon/c2recon_weights.log
#
# The weight. HGeant decays gamma* -> e+e- isotropically (cos theta* uniform); QED gives 1 + c^2 + e (1 - c^2),
# e = 4 m_e^2 / M(ee)^2, c = cos theta* (e+ in the gamma* rest frame against the gamma* direction, which in that frame is
# the line of the lost photon). Per row the reweighting factor is w = 3 [1 + c^2 + e (1 - c^2)] / (4 + 2 e) (mean 1 over
# the isotropic decay at fixed M(ee)).
#
# The lost photon is not stored. Given the four TRUE charged momenta (PAT), the photon k must satisfy
# (Q + k)^2 = m_pi0^2 and (P + k)^2 = m_eta^2 (Q = e+e-, P = pi+pi-e+e-). In the rest frame of Q this fixes |k| and the
# angle between k and P; one azimuth phi is free with a flat measure, and the generator density along it is the SMASH
# eta invariant density f(P + k) (HGeant: eta -> 3pi phase space, isotropic pi0 and gamma* decays). The event-by-event
# weight is therefore replaced by its conditional mean given everything the analysis sees,
#     <w>_row = sum_phi f(P + k(phi)) w(phi) / sum_phi f(P + k(phi)),
# which gives the same expected reweighted template as the true w for every quantity built from the charged tracks
# (c2recon_toy.py checks this closure on a toy). Where the ECAL simulation found the photon (cluster truth), the exact
# w is computed too.
#
# Variants of the density along the circle: A = SMASH etas of the ecalacc file (60 files, 140 572 etas; nominal),
# B = the verifier's SMASH etas (40 other files, 94 314 etas), A coarse (35 x 23 bins), flat (no density: extreme),
# A x KLOE-2 X dependence (physics of the data, not of the simulation; the analysis weights use Y only).
#
# Run: python3 -I plots/all_report/c2recon_weights.py
# ============================================================================
import math
import time

import numpy as np

O = "/hdd3/przygoda/HADES/FAT_out/all_report/c2recon/"
ETA_A = "/hdd3/przygoda/HADES/FAT_out/all_report/ecalacc/smash_eta_lab.txt"
ETA_B = "/hdd3/przygoda/HADES/FAT_out/all_report/vecalacc/vsmash_eta_lab.txt"
ECAL_ROWS = "/hdd3/przygoda/HADES/FAT_out/all_report/ecalacc/rows_c23.bin"
ECAL_LAY = "/hdd3/przygoda/HADES/FAT_out/all_report/ecalacc/rows_c23_layout.txt"
ME_E, MPI = 0.51099895, 139.57039            # MeV
META_G, MPI0_G = 547.45092, 134.97640        # HGeant masses measured by ecalacc from the found photons (MeV)
NPHI = 720
T0 = time.time()
LOG = None   # opened in main() only (this module is imported by the other c2recon scripts)


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    if LOG is not None:
        LOG.write(s + "\n"); LOG.flush()


COLS = ("tree idx set pat_entry grp s0 best cor_ok W108 W0 W116 W100 m4 m4_true mpipi mpipi_true mee mee_true e p ib tag m5c g_e "
        "ep_px ep_py ep_pz em_px em_py em_pz pip_px pip_py pip_pz pim_px pim_py pim_pz vx vy vz gfound g_px g_py g_pz g_ecl g_istag").split()


def load_rows():
    a = np.loadtxt(O + "c2recon_rows.txt", comments="#")
    assert a.shape[1] == len(COLS), a.shape
    return {k: a[:, i] for i, k in enumerate(COLS)}


# ---------------------------------------------------------------- four-vector helpers (E, px, py, pz), arrays (..., 4)
def p4(px, py, pz, m):
    return np.stack([np.sqrt(px * px + py * py + pz * pz + m * m), px, py, pz], -1)


def mass2(v):
    return v[..., 0] ** 2 - np.sum(v[..., 1:] ** 2, -1)


def boost(v, b):
    """boost the four-vectors v (..., 4) by the velocity b (..., 3) (v' = Lambda(b) v; use -beta to go to a rest frame)"""
    b2 = np.sum(b * b, -1)
    g = 1.0 / np.sqrt(1.0 - b2)
    bp = np.sum(b * v[..., 1:], -1)
    g2 = np.where(b2 > 0, (g - 1.0) / np.where(b2 > 0, b2, 1.0), 0.0)
    E = g * (v[..., 0] + bp)
    vec = v[..., 1:] + (g2 * bp)[..., None] * b + (g * v[..., 0])[..., None] * b
    return np.concatenate([E[..., None], vec], -1)


# ---------------------------------------------------------------- SMASH eta density E d3N/dp3 ~ dN/(dy dmT)/mT
class Dens:
    def __init__(self, fn, ny=70, nx=45, flat=False):
        self.flat = flat
        if flat:
            return
        a = np.loadtxt(fn)
        E, px, py, pz = a[:, 0], a[:, 1], a[:, 2], a[:, 3]
        m2 = E * E - px * px - py * py - pz * pz
        ok = (m2 > 0) & (E > np.abs(pz))
        E, px, py, pz, m2 = E[ok], px[ok], py[ok], pz[ok], m2[ok]
        m = np.sqrt(m2); mt = np.sqrt(px * px + py * py + m2)
        y = 0.5 * np.log((E + pz) / (E - pz))
        self.n = len(E); self.meanM = float(m.mean())
        self.ye = np.linspace(-0.5, 3.0, ny + 1); self.xe = np.linspace(0.0, 1.8, nx + 1)
        self.h, _, _ = np.histogram2d(y, mt - m, bins=[self.ye, self.xe], weights=1.0 / mt)
        self.yc = 0.5 * (self.ye[1:] + self.ye[:-1]); self.xc = 0.5 * (self.xe[1:] + self.xe[:-1])

    def __call__(self, pe):
        """pe: (..., 4) in GeV"""
        if self.flat:
            return np.ones(pe.shape[:-1])
        E, pz = pe[..., 0], pe[..., 3]
        okE = E > np.abs(pz)
        y = 0.5 * np.log(np.where(okE, (E + pz) / np.where(okE, E - pz, 1.0), 1.0))
        mt = np.sqrt(pe[..., 1] ** 2 + pe[..., 2] ** 2 + self.meanM ** 2)
        x = np.maximum(mt - self.meanM, self.xc[0])          # below the first centre: the first row (as ecalacc)
        ok = okE & (y > self.yc[0]) & (y < self.yc[-1]) & (x < self.xc[-1])
        dy = self.yc[1] - self.yc[0]; dx = self.xc[1] - self.xc[0]
        ty = np.clip((y - self.yc[0]) / dy, 0, len(self.yc) - 1.000001); tx = np.clip((x - self.xc[0]) / dx, 0, len(self.xc) - 1.000001)
        iy = ty.astype(int); ix = tx.astype(int); fy = ty - iy; fx = tx - ix
        iy1 = np.minimum(iy + 1, len(self.yc) - 1); ix1 = np.minimum(ix + 1, len(self.xc) - 1)
        h = self.h
        v = (1 - fy) * (1 - fx) * h[iy, ix] + fy * (1 - fx) * h[iy1, ix] + (1 - fy) * fx * h[iy, ix1] + fy * fx * h[iy1, ix1]
        return np.where(ok, np.maximum(v, 0.0), 0.0)


def kloe(X, Y):
    a, b, d, f, g = -1.095, 0.145, 0.081, 0.141, -0.044
    return 1 + a * Y + b * Y * Y + d * X * X + f * Y ** 3 + g * X * X * Y


def qed_w(c, aQ):
    return (1.0 + c * c + aQ * (1.0 - c * c)) * 3.0 / (4.0 + 2.0 * aQ)


def circle(ep, em, pip, pim, dens_list, meta=META_G, mpi0=MPI0_G, nphi=NPHI, seed=0xC2EC0, kloe_dens=None):
    """conditional mean of the QED weight along the photon circle for every row; ep.. (n, 4) in MeV.
    Returns dict: ok, and per density d: q[d] (mean w), plus the circle points of the nominal density for checks."""
    rng = np.random.default_rng(seed)
    n = ep.shape[0]
    Q = ep + em; Pv = Q + pip + pim
    Mq = np.sqrt(mass2(Q))
    bQ = Q[:, 1:] / Q[:, :1]
    Ps = boost(Pv, -bQ); eps = boost(ep, -bQ)
    w0 = (mpi0 * mpi0 - Mq * Mq) / (2.0 * Mq)
    c2 = 0.5 * (meta * meta - mass2(Pv))
    Pm = np.linalg.norm(Ps[:, 1:], axis=1)
    cosA = (Ps[:, 0] * w0 - c2) / (Pm * w0)
    ok = np.abs(cosA) <= 1.0
    cosA = np.clip(cosA, -1, 1); sinA = np.sqrt(1 - cosA * cosA)
    e1 = Ps[:, 1:] / Pm[:, None]
    tmp = np.where((np.abs(e1[:, 2]) < 0.9)[:, None], np.array([0, 0, 1.0]), np.array([1.0, 0, 0]))
    e2 = np.cross(tmp, e1); e2 /= np.linalg.norm(e2, axis=1)[:, None]
    e3 = np.cross(e1, e2)
    off = rng.random(n)
    ph = 2 * np.pi * (np.arange(nphi)[None, :] + off[:, None]) / nphi          # (n, nphi)
    ks = w0[:, None, None] * (cosA[:, None, None] * e1[:, None, :] + sinA[:, None, None] *
                              (np.cos(ph)[..., None] * e2[:, None, :] + np.sin(ph)[..., None] * e3[:, None, :]))
    aQ = 4.0 * ME_E * ME_E / (Mq * Mq)
    eu = eps[:, 1:] / np.linalg.norm(eps[:, 1:], axis=1)[:, None]
    ct = np.sum(eu[:, None, :] * ks, -1) / w0[:, None]
    wq = qed_w(ct, aQ[:, None])
    k4 = np.concatenate([w0[:, None, None] * np.ones((n, nphi, 1)), ks], -1)
    k4l = boost(k4, np.repeat(bQ[:, None, :], nphi, 1))
    pe = (Pv[:, None, :] + k4l) * 1e-3                                         # GeV
    out = {"ok": ok, "aQ": aQ, "wq_min": wq.min(1), "wq_max": wq.max(1)}
    for name, d in dens_list:
        f = d(pe)
        sf = f.sum(1)
        out["q_" + name] = np.where(sf > 0, (f * wq).sum(1) / np.where(sf > 0, sf, 1.0), np.nan)
        out["nz_" + name] = sf > 0
        if kloe_dens is not None and name == kloe_dens:
            # KLOE-2 X dependence along the circle (energies in the eta rest frame)
            be = pe[..., 1:] / pe[..., :1]
            a = boost(np.repeat(pip[:, None, :], nphi, 1) * 1e-3, -be); b = boost(np.repeat(pim[:, None, :], nphi, 1) * 1e-3, -be)
            cc = boost(np.repeat(Q[:, None, :], nphi, 1) * 1e-3 + k4l * 1e-3, -be)
            Qe = (meta - 2 * MPI - mpi0) * 1e-3
            X = math.sqrt(3) * ((a[..., 0] - MPI * 1e-3) - (b[..., 0] - MPI * 1e-3)) / Qe
            Y = 3 * (cc[..., 0] - mpi0 * 1e-3) / Qe - 1
            fk = f * kloe(X, Y)
            sk = fk.sum(1)
            out["q_" + name + "_kloe"] = np.where(sk > 0, (fk * wq).sum(1) / np.where(sk > 0, sk, 1.0), np.nan)
            out["Yspread"] = Y.max(1) - Y.min(1)
    return out


def exact(ep, em, gk):
    """exact QED weight from a known photon momentum gk (n, 3) MeV"""
    Q = ep + em
    bQ = Q[:, 1:] / Q[:, :1]
    g4 = np.concatenate([np.linalg.norm(gk, axis=1)[:, None], gk], 1)
    gs = boost(g4, -bQ); eps = boost(ep, -bQ)
    ct = np.sum(gs[:, 1:] * eps[:, 1:], 1) / (np.linalg.norm(gs[:, 1:], axis=1) * np.linalg.norm(eps[:, 1:], axis=1))
    Mq = np.sqrt(mass2(Q))
    return qed_w(ct, 4 * ME_E * ME_E / (Mq * Mq)), ct


def main():
    global LOG
    LOG = open(O + "c2recon_weights.log", "w")
    R = load_rows()
    n = len(R["tree"])
    P("# c2recon_weights.py, 2026-11-15: %d simulated C2 rows (pre, owned mode A); sum W108 = %.3f expected raw rows" % (n, R["W108"].sum()))
    ep = p4(R["ep_px"], R["ep_py"], R["ep_pz"], ME_E); em = p4(R["em_px"], R["em_py"], R["em_pz"], ME_E)
    pip = p4(R["pip_px"], R["pip_py"], R["pip_pz"], MPI); pim = p4(R["pim_px"], R["pim_py"], R["pim_pz"], MPI)
    dA = Dens(ETA_A); dB = Dens(ETA_B); dAc = Dens(ETA_A, ny=35, nx=23); dF = Dens(None, flat=True)
    P("   SMASH densities: A %d etas (mean mass %.5f), B %d etas (mean mass %.5f)" % (dA.n, dA.meanM, dB.n, dB.meanM))
    out = circle(ep, em, pip, pim, [("A", dA), ("B", dB), ("Ac", dAc), ("flat", dF)], kloe_dens="A")
    P("   circle: %d rows without a solution; zero density on the whole circle: A %d, B %d, Ac %d" %
      ((~out["ok"]).sum(), (~out["nz_A"]).sum(), (~out["nz_B"]).sum(), (~out["nz_Ac"]).sum()))
    # second phase offset (numerical check of the 720-point average)
    out2 = circle(ep, em, pip, pim, [("A", dA)], seed=0x5EC0, nphi=2880)
    dq = out2["q_A"] - out["q_A"]
    P("   720 vs 2880 phi points (other offsets): max |dq| %.2e, mean dq %.2e" % (np.nanmax(np.abs(dq)), np.nanmean(dq)))
    # exact weight where the photon is known
    gk = np.stack([R["g_px"], R["g_py"], R["g_pz"]], 1)
    fnd = R["gfound"] > 0
    qx = np.full(n, np.nan); ctx = np.full(n, np.nan)
    if fnd.any():
        qx[fnd], ctx[fnd] = exact(ep[fnd], em[fnd], gk[fnd])
        g4 = np.concatenate([np.linalg.norm(gk[fnd], axis=1)[:, None], gk[fnd]], 1)
        mp0 = np.sqrt(mass2(ep[fnd] + em[fnd] + g4)); meta = np.sqrt(mass2(ep[fnd] + em[fnd] + pip[fnd] + pim[fnd] + g4))
        P("   found photons: %d rows; M(e+e-gamma) median %.3f MeV (5-95 %%: %.3f-%.3f), M5 median %.3f MeV" %
          (fnd.sum(), np.median(mp0), np.percentile(mp0, 5), np.percentile(mp0, 95), np.median(meta)))
    W = R["W108"]
    s0 = R["s0"] > 0
    def wm(x, m):
        mm = m & np.isfinite(x)
        Wt = W[mm].sum()
        F = (W[mm] * x[mm]).sum() / Wt
        return F, math.sqrt((W[mm] ** 2 * (x[mm] - F) ** 2).sum()) / Wt
    for lab, m in (("pre", np.ones(n, bool)), ("S0", s0)):
        for key in ("q_A", "q_B", "q_Ac", "q_flat", "q_A_kloe"):
            P("   %-4s all C2 rows, %-9s: <w> = %.4f +- %.4f;  E1 %.4f, E2 %.4f" % (lab, key, *wm(out[key], m), wm(out[key], m & (R["mee"] < 0.010))[0],
                                                                         wm(out[key], m & (R["mee"] >= 0.010))[0]))
    # found-photon rows: exact vs circle (unrestricted circle; these rows have their photon in the ECAL)
    P("   found-photon rows (pre): exact <w> %.4f +- %.4f; circle A on the same rows %.4f" % (*wm(qx, fnd), wm(out["q_A"], fnd)[0]))
    # comparison with ecalacc (rows_c23.bin, meanQed) by (grp, m4_true)
    lay = open(ECAL_LAY).read().split("\n")
    nr, nd = int(lay[0].split()[1]), int(lay[0].split()[3])
    names = lay[1].split()[1:]
    a = np.fromfile(ECAL_ROWS, dtype=np.float64).reshape(nr, nd)
    E = {k: a[:, i] for i, k in enumerate(names)}
    c2e = E["cls"] == 2
    key_e = {(int(g), round(float(m), 6)): i for i, (g, m) in enumerate(zip(E["grp"][c2e], E["m4true"][c2e]))}
    idx_e = np.where(c2e)[0]
    jm = np.array([key_e.get((int(g), round(float(np.float32(m)), 6)), -1) for g, m in zip(R["grp"], R["m4_true"])])
    nmatch = (jm >= 0).sum()
    q_ecal = np.full(n, np.nan); tq_ecal = np.full(n, np.nan)
    q_ecal[jm >= 0] = E["meanQed"][idx_e[jm[jm >= 0]]]
    tq_ecal[jm >= 0] = np.where(E["tMatch"][idx_e[jm[jm >= 0]]] > 0, E["tQed"][idx_e[jm[jm >= 0]]], np.nan)
    d = out["q_A"] - q_ecal
    P("   ecalacc meanQed matched for %d of %d rows: own - ecalacc: mean %.2e, rms %.2e, max |d| %.2e; weighted <w>: own %.5f, ecalacc %.5f" %
      (nmatch, n, np.nanmean(d), np.nanstd(d), np.nanmax(np.abs(d)), wm(out["q_A"], jm >= 0)[0], wm(q_ecal, jm >= 0)[0]))
    dt = qx - tq_ecal
    P("   exact weight of the found photons: own - ecalacc tQed: max |d| %.2e (%d rows)" % (np.nanmax(np.abs(dt)), np.isfinite(dt).sum()))
    # per-row spread of the circle average: how much information the charged tracks carry about cos theta*
    P("   per-row <w> (A): 5/50/95 %% = %.3f / %.3f / %.3f; range on the circle (min..max of w) median %.3f..%.3f" %
      (*np.nanpercentile(out["q_A"], [5, 50, 95]), np.median(out["wq_min"]), np.median(out["wq_max"])))
    P("   M(pi+pi-) is fixed by the charged tracks: Y spread along the circle max %.2e (KLOE weight then varies only through X)" % np.nanmax(out["Yspread"]))
    np.savez_compressed(O + "c2recon_weights.npz", tree=R["tree"], idx=R["idx"], grp=R["grp"], m4=R["m4"], m4_true=R["m4_true"],
                        mee=R["mee"], W108=W, s0=R["s0"], q_A=out["q_A"], q_B=out["q_B"], q_Ac=out["q_Ac"], q_flat=out["q_flat"],
                        q_A_kloe=out["q_A_kloe"], q_A2880=out2["q_A"], q_exact=qx, ct_exact=ctx, gfound=R["gfound"], g_istag=R["g_istag"],
                        q_ecalacc=q_ecal, ok=out["ok"])
    with open(O + "c2recon_weights.txt", "w") as f:
        f.write("# tree idx q_A q_B q_Ac q_flat q_A_kloe (per C2 row; c2recon_weights.py)\n")
        for i in range(n):
            vals = [out[k][i] if np.isfinite(out[k][i]) else 1.0 for k in ("q_A", "q_B", "q_Ac", "q_flat", "q_A_kloe")]
            f.write("%d %d %.8f %.8f %.8f %.8f %.8f\n" % (int(R["tree"][i]), int(R["idx"][i]), *vals))
    P("   written c2recon_weights.npz and c2recon_weights.txt (%.0f s)" % (time.time() - T0))


if __name__ == "__main__":
    main()
