#!/usr/bin/env python3
# ============================================================================
# c2recon_toy.py  --  package c2recon: closure test of the circle-averaged QED weight on a toy
# ============================================================================
# Written 2026-11-15 (prefix c2recon). New file. A TOY, not a detector simulation.
#
# Question: c2recon_weights.py replaces the event-by-event QED weight w(cos theta*) of a simulated C2 row (the lost
# photon is not stored) by its conditional mean along the photon circle, given the four true charged momenta and the
# SMASH eta density. Does the reweighted template close, i.e. is sum <w>_circle = sum w_true for any selection made
# from the charged tracks?
#
# Toy (own code, modelled on the WP8 generator): eta lab momenta drawn from the SMASH etas of the verifier's file
# (B, 94 314 etas); eta -> pi+pi-pi0 phase space (HGeant); pi0 -> gamma gamma*, M(ee) from Kroll-Wada x (1 + a x)^2;
# gamma* -> e+e- isotropic (HGeant); w_true = 3 [1 + c^2 + e (1 - c^2)] / (4 + 2 e). Acceptance: four tracks 18-85 deg,
# pions 0.10 < p < 2.0 GeV/c, leptons p_min < p < 2.0 GeV/c (p_min 0.04 / 0.06 / 0.10 / 0.15).
# For the accepted decays the circle average is computed with the density of file A (independent SMASH files), of
# file B (the generator sample) and with a flat density. Subsets: all, E1/E2, near the edge (u_true > -15 MeV), far.
# Also: the factor against the momentum of the softer lepton (to explain toy vs full simulation).
#
# Output: /hdd3/przygoda/HADES/FAT_out/all_report/c2recon/c2recon_toy.json, c2recon_toy.log
# Run: python3 -I plots/all_report/c2recon_toy.py   (about 1-2 min)
# ============================================================================
import json
import math
import sys
import time

import numpy as np

sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import c2recon_weights as W

O = "/hdd3/przygoda/HADES/FAT_out/all_report/c2recon/"
ME, MPI, MPI0, META = 0.51099895e-3, 0.13957039, 0.1349768, 0.547862   # GeV (toy eta mass = PDG; circle uses the same)
A_FF = 0.0332
rng = np.random.default_rng(20261115 + 2)
LOG = open(O + "c2recon_toy.log", "w")
T0 = time.time()


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s); LOG.write(s + "\n"); LOG.flush()


def iso_dir(n):
    c = 2 * rng.random(n) - 1
    s = np.sqrt(1 - c * c)
    ph = 2 * np.pi * rng.random(n)
    return np.stack([s * np.cos(ph), s * np.sin(ph), c], axis=1)


def sample_q(n):
    out = np.empty(0)
    qmin, qmax = 2 * ME * 1.000001, MPI0 * 0.999999
    while out.size < n:
        m = 4 * (n - out.size) + 1000
        u = rng.random(m)
        q = qmin * (qmax / qmin) ** u
        x = (q / MPI0) ** 2
        beta = np.sqrt(1 - 4 * ME * ME / (q * q))
        f = (1 - x) ** 3 * (1 + 2 * ME * ME / (q * q)) * beta * (1 + A_FF * x) ** 2
        keep = rng.random(m) * 2.0 < f
        out = np.concatenate([out, q[keep]])
    return out[:n]


def two_body(parent_m, m1, m2, n):
    """momentum and isotropic directions of a two-body decay at rest"""
    pm = np.sqrt(np.maximum(0, (parent_m ** 2 - (m1 + m2) ** 2) * (parent_m ** 2 - (m1 - m2) ** 2))) / (2 * parent_m)
    d = iso_dir(n)
    return pm, d


def generate(n, etas):
    # eta lab four-momenta from the SMASH list (3-momentum kept, mass set to META)
    j = rng.integers(0, etas.shape[0], n)
    p3 = etas[j, 1:4]
    eta = np.column_stack([np.sqrt(np.sum(p3 * p3, 1) + META ** 2), p3])
    # eta -> pi+ pi- pi0 phase space: uniform in the Dalitz plot (s12 = M(pipi)^2, s13)
    s12lo, s12hi = (2 * MPI) ** 2, (META - MPI0) ** 2
    s12 = np.empty(0); s13 = np.empty(0)
    while s12.size < n:
        m = 3 * (n - s12.size) + 1000
        a = s12lo + (s12hi - s12lo) * rng.random(m)
        b = (MPI + MPI0) ** 2 + ((META - MPI) ** 2 - (MPI + MPI0) ** 2) * rng.random(m)
        E2s = np.sqrt(a) / 2
        E3s = (META ** 2 - a - MPI0 ** 2) / (2 * np.sqrt(a))
        p2s = np.sqrt(np.maximum(0, E2s ** 2 - MPI ** 2)); p3s = np.sqrt(np.maximum(0, E3s ** 2 - MPI0 ** 2))
        lo = (E2s + E3s) ** 2 - (p2s + p3s) ** 2
        hi = (E2s + E3s) ** 2 - (p2s - p3s) ** 2
        ok = (b >= lo) & (b <= hi)
        s12 = np.concatenate([s12, a[ok]]); s13 = np.concatenate([s13, b[ok]])
    s12 = s12[:n]
    mpp = np.sqrt(s12)
    E0 = (META ** 2 + MPI0 ** 2 - s12) / (2 * META)
    p0 = np.sqrt(np.maximum(0, E0 ** 2 - MPI0 ** 2))
    d0 = iso_dir(n)
    pi0 = np.column_stack([E0, p0[:, None] * d0])
    pp = np.column_stack([META - E0, -p0[:, None] * d0])
    kp = np.sqrt(np.maximum(0, mpp ** 2 / 4 - MPI ** 2))
    dpi = iso_dir(n)
    pip = np.column_stack([np.sqrt(kp ** 2 + MPI ** 2), kp[:, None] * dpi])
    pim = np.column_stack([np.sqrt(kp ** 2 + MPI ** 2), -kp[:, None] * dpi])
    bpp = pp[:, 1:] / pp[:, :1]
    pip = W.boost(pip, bpp); pim = W.boost(pim, bpp)
    # pi0 -> gamma gamma*: gamma* direction dg in the pi0 frame
    q = sample_q(n)
    ps = (MPI0 ** 2 - q * q) / (2 * MPI0)
    dg = iso_dir(n)
    gs = np.column_stack([np.sqrt(ps * ps + q * q), ps[:, None] * dg])
    gam = np.column_stack([ps, -ps[:, None] * dg])
    # gamma* -> e+e- isotropic in the gamma* frame; c = cos(e+, gamma* direction)
    k = np.sqrt(np.maximum(0, q * q / 4 - ME * ME))
    c = 2 * rng.random(n) - 1
    s = np.sqrt(1 - c * c)
    ph = 2 * np.pi * rng.random(n)
    tmp = np.where(np.abs(dg[:, :1]) < 0.9, np.array([[1.0, 0, 0]]), np.array([[0, 1.0, 0]]))
    ex = np.cross(tmp, dg); ex /= np.linalg.norm(ex, axis=1)[:, None]
    ey = np.cross(dg, ex)
    dl = c[:, None] * dg + (s * np.cos(ph))[:, None] * ex + (s * np.sin(ph))[:, None] * ey
    Ek = np.sqrt(k * k + ME * ME)
    ep = np.column_stack([Ek, k[:, None] * dl]); em = np.column_stack([Ek, -k[:, None] * dl])
    bg = gs[:, 1:] / gs[:, :1]
    ep = W.boost(ep, bg); em = W.boost(em, bg)
    aQ = 4 * ME * ME / (q * q)
    wq = W.qed_w(c, aQ)
    b0 = pi0[:, 1:] / pi0[:, :1]
    ep = W.boost(ep, b0); em = W.boost(em, b0); gam = W.boost(gam, b0)
    be = eta[:, 1:] / eta[:, :1]
    ep, em, pip, pim, gam = W.boost(ep, be), W.boost(em, be), W.boost(pip, be), W.boost(pim, be), W.boost(gam, be)
    P4 = ep + em + pip + pim
    m4 = np.sqrt(np.maximum(0, W.mass2(P4)))
    return dict(q=q, wq=wq, ep=ep, em=em, pip=pip, pim=pim, gam=gam, mpp=mpp, m4=m4)


def accept(ev, pmin_e):
    ok = np.ones(ev["q"].size, bool)
    for key, pmin in (("ep", pmin_e), ("em", pmin_e), ("pip", 0.10), ("pim", 0.10)):
        v = ev[key][:, 1:]
        p = np.linalg.norm(v, axis=1)
        th = np.degrees(np.arccos(np.clip(v[:, 2] / p, -1, 1)))
        ok &= (th > 18) & (th < 85) & (p > pmin) & (p < 2.0)
    return ok


def m4max(m, q):
    E0 = (META ** 2 + MPI0 ** 2 - m * m) / (2 * META)
    p2 = np.maximum(E0 * E0 - MPI0 ** 2, 0)
    return np.sqrt(META ** 2 - META * (E0 - np.sqrt(p2)) * (1 - q * q / MPI0 ** 2))


def main():
    etasA = np.loadtxt(W.ETA_A); etasB = np.loadtxt(W.ETA_B)
    dA = W.Dens(W.ETA_A); dB = W.Dens(W.ETA_B); dF = W.Dens(None, flat=True)
    nchunk, n = 40, 100000
    pmins = (0.04, 0.06, 0.10, 0.15)
    acc = {p: {k: [] for k in ("wq", "qA", "qB", "qF", "q", "u", "pmin_e", "mpp")} for p in pmins}
    allw = []
    for ic in range(nchunk):
        ev = generate(n, etasB)
        allw.append(ev["wq"].mean())
        ok_any = accept(ev, min(pmins))
        idx = np.where(ok_any)[0]
        if idx.size == 0:
            continue
        sub = {k: v[idx] for k, v in ev.items()}
        mev = lambda v: v * 1e3
        out = W.circle(mev(sub["ep"]), mev(sub["em"]), mev(sub["pip"]), mev(sub["pim"]), [("A", dA), ("B", dB), ("flat", dF)],
                       meta=META * 1e3, mpi0=MPI0 * 1e3, nphi=180, seed=ic + 1)
        pe = np.minimum(np.linalg.norm(sub["ep"][:, 1:], axis=1), np.linalg.norm(sub["em"][:, 1:], axis=1))
        u = sub["m4"] - m4max(sub["mpp"], sub["q"])
        for p in pmins:
            m = accept(sub, p)
            acc[p]["wq"].append(sub["wq"][m]); acc[p]["qA"].append(out["q_A"][m]); acc[p]["qB"].append(out["q_B"][m])
            acc[p]["qF"].append(out["q_flat"][m]); acc[p]["q"].append(sub["q"][m]); acc[p]["u"].append(u[m])
            acc[p]["pmin_e"].append(pe[m]); acc[p]["mpp"].append(sub["mpp"][m])
    P("# c2recon_toy.py: %d x %d toy decays (eta from SMASH file B); <w_true> over all decays %.5f (expected 1)" % (nchunk, n, np.mean(allw)))
    res = {}
    for p in pmins:
        A = {k: np.concatenate(v) for k, v in acc[p].items()}
        nacc = A["wq"].size
        res[p] = {}
        subsets = {"all": np.ones(nacc, bool), "E1": A["q"] < 0.010, "E2": A["q"] >= 0.010,
                   "edge (u_true > -15 MeV)": A["u"] > -0.015, "far (u_true < -60 MeV)": A["u"] < -0.060,
                   "M(pipi) 0.30-0.40": (A["mpp"] >= 0.30) & (A["mpp"] < 0.40)}
        P("\np_e > %.2f GeV/c: %d accepted decays" % (p, nacc))
        P("   subset                     | true <w>          | circle A (indep.)  | circle B (gen.)  | flat     | A - true (sigma)")
        for name, m in subsets.items():
            wt = A["wq"][m]; nn = m.sum()
            et = wt.std() / math.sqrt(nn)
            # paired error of the difference (same decays)
            dAe = (A["qA"][m] - wt).std() / math.sqrt(nn)
            r = dict(n=int(nn), true=float(wt.mean()), etrue=float(et), A=float(np.nanmean(A["qA"][m])), B=float(np.nanmean(A["qB"][m])),
                     flat=float(np.nanmean(A["qF"][m])), dA=float(np.nanmean(A["qA"][m]) - wt.mean()), edA=float(dAe))
            res[p][name] = r
            P("   %-26s | %.4f +- %.4f | %.4f            | %.4f          | %.4f   | %+.4f (%+.1f)" % (name, r["true"], et, r["A"], r["B"], r["flat"], r["dA"], r["dA"] / dAe))
        # factor against the softer lepton momentum
        pb = [0.04, 0.06, 0.08, 0.10, 0.15, 0.20, 0.30, 0.50, 2.0]
        rowsb = []
        for a, b in zip(pb[:-1], pb[1:]):
            m = (A["pmin_e"] >= a) & (A["pmin_e"] < b)
            if m.sum() < 100:
                continue
            rowsb.append((a, b, int(m.sum()), float(A["wq"][m].mean()), float(np.nanmean(A["qA"][m]))))
        res[p]["vs_pmin"] = rowsb
        P("   factor vs softer lepton momentum: " + "; ".join("%.2f-%.2f: %.3f (circle %.3f, n %d)" % (a, b, t, c, nn) for a, b, nn, t, c in rowsb))
    json.dump({str(k): v for k, v in res.items()}, open(O + "c2recon_toy.json", "w"), indent=1)
    P("\ndone in %.0f s" % (time.time() - T0))


if __name__ == "__main__":
    main()
