#!/usr/bin/env python3
# vecalacc_frac.py -- INDEPENDENT VERIFIER of package ecalacc (2026-11-15): where does the lost photon go?
# Own implementation (no code of ecalacc_fill.C / ecalacc_ana.py):
#   rows   : vrows.txt of vecalacc_rows.C (own selection with defs.h; true momenta read from the PAT files)
#   geometry: vecalacc_geo.py (own dump, own grid hit test, all six planes tried)
#   SMASH  : own extraction of the eta momenta from 40 OTHER converted files (1001-1040) of the standard set;
#            invariant density E d3N/dp3 = dN/(dy d pT^2) / pi estimated as a smoothed histogram in (y, pT^2)
#            (the package uses files 1..~60 and dN/(dy dmT)/mT); the package's file is a variant.
#   completion: solved in the LAB frame (the package works in the e+e- rest frame). With k = E n:
#            2 Q.k = m_pi0^2 - Q^2 = a, 2 P.k = m_eta^2 - P^2 = b  =>  (b p_Q - a p_P).n = b E_Q - a E_P,
#            i.e. n lies on a cone around c = b p_Q - a p_P, and E = a / (2 (E_Q - p_Q.n)).
#            Measure on the cone azimuth psi: d^3k/k delta delta -> E(psi) dpsi (derived in the results note and
#            checked numerically below against a flat azimuth in the e+e- rest frame).
#   weight : f_SMASH(P + k) x |M|^2, |M|^2 = 1 (HGeant: phase space, isotropic gamma*).
# Output: /hdd3/przygoda/HADES/FAT_out/all_report/vecalacc/vfrac_numbers.json, vfrac_rows.npz and printout.
import json, sys
import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.interpolate import RegularGridInterpolator

sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import vecalacc_geo as G

O = "/hdd3/przygoda/HADES/FAT_out/all_report/vecalacc"
ME, MPI = 0.51099895, 139.57039          # MeV
KMETA, KMPI0 = 0.547862, 0.1349768        # GeV, ScanEngine constants (edge function)
NPSI = 720


# ------------------------------------------------------------------ rows
def read_rows(fn=O + "/vrows.txt", c5every=10):
    rows = []
    n5 = 0
    with open(fn) as f:
        for line in f:
            if line[0] == "#":
                continue
            p = line.split("|")
            h = p[0].split()
            cls = int(h[2])
            if cls == 5:
                n5 += 1
                if n5 % c5every:
                    continue
            rows.append((list(map(float, h)), list(map(float, p[1].split())), list(map(float, p[2].split())),
                         np.array(list(map(float, p[3].split()))).reshape(18, 5)))
    return rows, c5every


def edge(m, q):
    mpr = m / 0.9970
    E0 = (KMETA ** 2 + KMPI0 ** 2 - mpr ** 2) / (2 * KMETA)
    p2 = E0 ** 2 - KMPI0 ** 2
    if p2 < 0:
        return -1.0
    a = KMETA ** 2 - KMETA * (E0 - np.sqrt(p2)) * (1 - q * q / KMPI0 ** 2)
    return 0.99274 * np.sqrt(a) if a > 0 else -1.0


def four(px, py, pz, m):
    return np.array([np.sqrt(px * px + py * py + pz * pz + m * m), px, py, pz])


def mdot(a, b):
    return a[..., 0] * b[..., 0] - (a[..., 1:] * b[..., 1:]).sum(-1)


def boost_to_rest(p, frame):
    """boost four-vectors p (...,4) into the rest frame of the four-vector frame (4,)"""
    m = np.sqrt(mdot(frame, frame))
    bvec = frame[1:] / frame[0]
    b2 = bvec @ bvec
    g = frame[0] / m
    bp = p[..., 1:] @ bvec
    g2 = (g - 1) / b2 if b2 > 0 else 0.0
    E = g * (p[..., 0] - bp)
    v = p[..., 1:] + (g2 * bp - g * p[..., 0])[..., None] * bvec
    return np.concatenate([E[..., None], v], -1)


# ------------------------------------------------------------------ SMASH density
class Dens:
    def __init__(self, fn, ny=50, npt=48):
        a = np.loadtxt(fn, usecols=(0, 1, 2, 3))
        E, px, py, pz = a.T
        y = 0.5 * np.log((E + pz) / (E - pz))
        pt2 = px * px + py * py
        self.n = len(a)
        self.yb = np.linspace(-0.5, 3.0, ny + 1); self.pb = np.linspace(0.0, 1.6, npt + 1)
        H, _, _ = np.histogram2d(y, pt2, bins=[self.yb, self.pb])
        H = gaussian_filter(H, sigma=0.8, mode="nearest")
        yc = 0.5 * (self.yb[1:] + self.yb[:-1]); pc = 0.5 * (self.pb[1:] + self.pb[:-1])
        self.f = RegularGridInterpolator((yc, pc), H, bounds_error=False, fill_value=0.0)
        self.ymin, self.ymax = yc[0], yc[-1]

    def __call__(self, p):  # p (...,4) in MeV
        E, px, py, pz = p[..., 0] / 1e3, p[..., 1] / 1e3, p[..., 2] / 1e3, p[..., 3] / 1e3
        with np.errstate(divide="ignore", invalid="ignore"):
            y = 0.5 * np.log((E + pz) / (E - pz))
        pt2 = px * px + py * py
        pt2c = np.clip(pt2, 0.5 * (self.pb[0] + self.pb[1]), None)   # below the first bin centre: first bin value
        v = self.f(np.stack([np.clip(y, -10, 10), pt2c], -1))
        v[(y < self.ymin) | (y > self.ymax) | ~np.isfinite(y)] = 0.0
        return np.maximum(v, 0.0)


# ------------------------------------------------------------------ completion of one row
def cone(Q, P, meta, mpi0):
    a = mpi0 ** 2 - mdot(Q, Q); b = meta ** 2 - mdot(P, P)
    c = b * Q[1:] - a * P[1:]
    cn = np.linalg.norm(c)
    cb = (b * Q[0] - a * P[0]) / cn
    return a, b, c / cn, cb


def complete(Q, P, meta, mpi0, psi):
    a, b, ch, cb = cone(Q, P, meta, mpi0)
    if not (abs(cb) <= 1.0) or a <= 0 or b <= 0:
        return None
    sb = np.sqrt(1 - cb * cb)
    t = np.array([0, 0, 1.0]) if abs(ch[2]) < 0.9 else np.array([1.0, 0, 0])
    e1 = np.cross(t, ch); e1 /= np.linalg.norm(e1); e2 = np.cross(ch, e1)
    n = cb * ch[None, :] + sb * (np.cos(psi)[:, None] * e1[None, :] + np.sin(psi)[:, None] * e2[None, :])
    E = a / (2 * (Q[0] - n @ Q[1:]))
    k = np.concatenate([E[:, None], E[:, None] * n], 1)
    return n, E, k


def main():
    geo = G.Geo(G.load_geo()[0])
    smash_mine = O + "/vsmash_eta_lab.txt"
    for a_ in sys.argv[1:]:                      # test runs: --mine=<file>
        if a_.startswith("--mine="):
            smash_mine = a_.split("=", 1)[1]
    smash_pkg = "/hdd3/przygoda/HADES/FAT_out/all_report/ecalacc/smash_eta_lab.txt"
    dens = {"mine": Dens(smash_mine), "pkg": Dens(smash_pkg)}
    print(f"SMASH etas: own files {dens['mine'].n}, package file {dens['pkg'].n}")
    rows, c5every = read_rows()
    print(f"rows read: {len(rows)} (C5 every {c5every}th row)")
    # generator masses from the found photons (own matching)
    def match(sc, nu):
        best, bE = -1, -1
        for k in range(min(5, int(round(sc[23])))):
            if int(round(nu[7, k])) != 1 or int(round(nu[8, k])) != 7:
                continue
            if nu[9, k] != sc[18] or nu[10, k] != sc[19] or nu[11, k] != sc[20]:
                continue
            if np.sqrt((nu[15, k] - sc[12]) ** 2 + (nu[16, k] - sc[13]) ** 2 + (nu[17, k] - sc[14]) ** 2) > 0.01:
                continue
            if nu[2, k] > bE:
                best, bE = k, nu[2, k]
        return best
    meta_v, mpi0_v = [], []
    for h, cp, sc, nu in rows:
        if int(h[2]) not in (2, 3):
            continue
        k = match(sc, nu)
        if k < 0:
            continue
        ep, em = four(*sc[0:3], ME), four(*sc[3:6], ME); pp, pm = four(*sc[6:9], MPI), four(*sc[9:12], MPI)
        g = four(nu[12, k], nu[13, k], nu[14, k], 0.0)
        tot = ep + em + pp + pm + g
        meta_v.append(np.sqrt(mdot(tot, tot)))
        if int(h[2]) == 2:
            t2 = ep + em + g; mpi0_v.append(np.sqrt(mdot(t2, t2)))
    meta, mpi0 = float(np.median(meta_v)), float(np.median(mpi0_v))
    print(f"generator masses from {len(meta_v)} / {len(mpi0_v)} found photons: m_eta {meta:.3f}, m_pi0 {mpi0:.3f} MeV")

    psi = (np.arange(NPSI) + 0.5) * 2 * np.pi / NPSI
    out = dict(meta=meta, mpi0=mpi0)
    R = []   # per-row results for C2/C3
    thH = np.zeros(90); thHacc = np.zeros(90)   # weighted polar-angle spectra of the lost photon, 0-180 deg, 2 deg
    nosol = 0
    # numerical check of the measure on the first 200 C2 rows: flat azimuth in the e+e- rest frame -> lab psi density ~ E
    chk = []
    for h, cp, sc, nu in rows:
        cls = int(h[2])
        if cls not in (2, 3):
            continue
        w = h[4]; s0 = int(h[3]); m4, mpipi, mee = h[5], h[7], h[8]
        ep, em = four(*sc[0:3], ME), four(*sc[3:6], ME); pp, pm = four(*sc[6:9], MPI), four(*sc[9:12], MPI)
        Q = ep + em; P = Q + pp + pm
        V = np.array(sc[12:15])
        res = complete(Q, P, meta, mpi0, psi)
        if res is None:
            nosol += 1
            continue
        n, E, k = res
        # measure check (rest frame flat phi -> lab), first 200 rows
        if len(chk) < 200:
            a_, b_, ch, cb = cone(Q, P, meta, mpi0)
            Qr = Q[0] * 0  # placeholder to keep names explicit
            Ps = boost_to_rest(P, Q); ws = (mpi0 ** 2 - mdot(Q, Q)) / (2 * np.sqrt(mdot(Q, Q)))
            ca = (Ps[0] * ws - 0.5 * (meta ** 2 - mdot(P, P))) / (np.linalg.norm(Ps[1:]) * ws)
            u = Ps[1:] / np.linalg.norm(Ps[1:]); t = np.array([0, 0, 1.0]) if abs(u[2]) < 0.9 else np.array([1.0, 0, 0])
            f1 = np.cross(t, u); f1 /= np.linalg.norm(f1); f2 = np.cross(u, f1)
            phis = np.random.default_rng(len(chk)).uniform(0, 2 * np.pi, 20000)
            sa = np.sqrt(max(0.0, 1 - ca * ca))
            ks = ws * np.concatenate([np.ones((len(phis), 1)), ca * u[None, :] + sa * (np.cos(phis)[:, None] * f1 + np.sin(phis)[:, None] * f2)], 1)
            # boost back to the lab: inverse boost = boost into the frame moving with -beta
            Qinv = np.array([Q[0], -Q[1], -Q[2], -Q[3]])
            kl = boost_to_rest(ks, Qinv)
            nl = kl[:, 1:] / kl[:, :1]
            t2 = np.array([0, 0, 1.0]) if abs(ch[2]) < 0.9 else np.array([1.0, 0, 0])
            e1 = np.cross(t2, ch); e1 /= np.linalg.norm(e1); e2 = np.cross(ch, e1)
            ps = np.arctan2(nl @ e2, nl @ e1) % (2 * np.pi)
            hh, _ = np.histogram(ps, bins=36, range=(0, 2 * np.pi))
            pred = np.interp((np.arange(36) + 0.5) * 2 * np.pi / 36, psi, E); pred = pred / pred.sum() * hh.sum()
            chk.append(((hh - pred) ** 2 / np.maximum(pred, 1)).sum() / 35)
        dist = 0.99615 * m4 - edge(0.9952 * mpipi, mee)
        # geometry
        Vr = np.tile(V, (NPSI, 1))
        _, _, _, act = geo.hit(Vr, n)
        _, _, _, act6 = geo.hit(Vr, n, six=True)
        _, _, _, actS = geo.hit(Vr, n, half=46.0)
        th = np.degrees(np.arccos(np.clip(n[:, 2], -1, 1)))
        meas = E                                   # lab measure E dpsi
        fr = {}
        for dn in ("mine", "pkg"):
            f = dens[dn](P[None, :] + k) * meas
            fr[dn] = f
        fflat = meas
        f = fr["mine"]; sf = f.sum()
        if not sf > 0:
            nosol += 1
            continue
        # QED gamma* weight (C2)
        if cls == 2:
            eps = boost_to_rest(ep, Q)[1:]; eps /= np.linalg.norm(eps)
            ks_ = boost_to_rest(k, Q)[:, 1:]; ks_ /= np.linalg.norm(ks_, axis=1)[:, None]
            ct = ks_ @ eps; aq = 4 * ME * ME / mdot(Q, Q)
            wq = (1 + ct * ct + aq * (1 - ct * ct)) * 3 / (4 + 2 * aq)
            qmean = (f * wq).sum() / sf; qacc = (f * wq * act).sum() / max((f * wq).sum(), 1e-300)
        else:
            qmean, qacc = np.nan, np.nan
        km = match(sc, nu)
        found = km >= 0
        rec = dict(cls=cls, s0=s0, w=w, m4p=0.99615 * m4, mpp=0.9952 * mpipi, mee=mee, dist=dist,
                   fGeo=(f * act).sum() / sf, fGeo6=(f * act6).sum() / sf, fGeoS=(f * actS).sum() / sf,
                   fGeoPkg=(fr["pkg"] * act).sum() / fr["pkg"].sum(), fGeoFlat=(fflat * act).sum() / fflat.sum(),
                   fThr80=(f * act * (E >= 80)).sum() / sf, fThr100=(f * act * (E >= 100)).sum() / sf,
                   fThr130=(f * act * (E >= 130)).sum() / sf,
                   fLow=(f * (th < 11.5)).sum() / sf, fBand=(f * ((th >= 11.5) & (th <= 47.5))).sum() / sf,
                   fHigh=(f * (th > 47.5)).sum() / sf,
                   fBandV=(f * ((th >= 10.9) & (th <= 45.6))).sum() / sf, fHighV=(f * (th > 45.6)).sum() / sf,
                   fTrk=(f * ((th >= 18) & (th <= 85))).sum() / sf, fBack=(f * (th > 90)).sum() / sf,
                   # the package's definition (2-deg bin centres of a 0-90 deg histogram, normalised to theta < 90 deg)
                   pLow=(f * (th < 12)).sum() / max((f * (th < 90)).sum(), 1e-300),
                   pBand=(f * ((th >= 12) & (th < 48))).sum() / max((f * (th < 90)).sum(), 1e-300),
                   pHigh=(f * ((th >= 48) & (th < 90))).sum() / max((f * (th < 90)).sum(), 1e-300),
                   pTrk=(f * ((th >= 18) & (th < 86))).sum() / max((f * (th < 90)).sum(), 1e-300),
                   eMean=(f * E).sum() / sf, eMeanAcc=(f * E * act).sum() / max((f * act).sum(), 1e-300),
                   eBelow100=(f * (E < 100)).sum() / sf, eBelow100Acc=(f * act * (E < 100)).sum() / max((f * act).sum(), 1e-300),
                   qmean=qmean, qacc=qacc, found=int(found))
        # closure for the found photon
        if found:
            g = np.array([nu[12, km], nu[13, km], nu[14, km]]); Eg = np.linalg.norm(g); ng = g / Eg
            a_, b_, ch, cb = cone(Q, P, meta, mpi0)
            rec.update(cRes=float(ng @ ch - cb), eRes=float(Eg - a_ / (2 * (Q[0] - ng @ Q[1:]))) , eTrue=Eg, eCl=nu[2, km],
                       tGeo=int(geo.hit(V[None, :], ng[None, :])[3][0]), tGeo6=int(geo.hit(V[None, :], ng[None, :], six=True)[3][0]),
                       tPass=int(int(round(nu[0, km])) == 1 and 0.8 < nu[1, km] < 1.2 and nu[2, km] > 100.))
        # tag from the compact columns (frozen definition of WP6), true-photon tag
        gidx, ge, gth, gph, meeg, m5 = cp
        tag = gidx > 0 and mee < 0.100 and 0.1190 <= meeg < 0.1730
        m5c = -1.0
        if tag:
            aa = meeg * meeg - mee * mee
            if aa > 0:
                ff = (KMPI0 ** 2 - mee * mee) / aa
                m2 = m4 * m4 + ff * (m5 * m5 - m4 * m4)
                m5c = np.sqrt(m2) if m2 > 0 else -1.0
        rec["tagW"] = int(tag and 0.49 <= m5c < 0.62)
        rec["tagTrue"] = int(rec["tagW"] and found and abs(ge - nu[2, km]) < 1e-3 * max(1.0, nu[2, km]))
        hb = np.minimum((th / 2).astype(int), 89)
        thH += w * np.bincount(hb, weights=f / sf, minlength=90); thHacc += w * np.bincount(hb, weights=f * act / sf, minlength=90)
        rec["foundPass"] = int(found and rec.get("tPass", 0) == 1)
        rec["nmult"] = sc[23]
        R.append(rec)
    print(f"C2+C3 rows completed: {len(R)}; no solution or zero density: {nosol}")
    chk = np.array(chk)
    print(f"measure check (rest-frame flat azimuth vs lab E dpsi, 36 bins, 20000 points per row): chi2/ndf median {np.median(chk):.2f}, "
          f"95 % {np.quantile(chk, 0.95):.2f} over {len(chk)} rows")
    out["measure_check_chi2ndf_median"] = float(np.median(chk))
    keys = list(R[0].keys())
    A = {k_: np.array([r.get(k_, np.nan) for r in R], float) for k_ in set().union(*[r.keys() for r in R])}
    np.savez(O + "/vfrac_rows.npz", thH=thH, thHacc=thHacc, **A)

    def wfrac(key, sel):
        w = A["w"][sel]; x = A[key][sel]
        m = np.isfinite(x); w = w[m]; x = x[m]
        f = (w * x).sum() / w.sum()
        e = np.sqrt((w * w * (x - f) ** 2).sum()) / w.sum()
        return f, e

    allm = np.ones(len(R), bool)
    sels = {"C2": A["cls"] == 2, "C3": A["cls"] == 3, "C2+C3": allm,
            "B": (A["m4p"] >= 0.4550) & (A["m4p"] < 0.5225),
            "edge": (A["dist"] >= -0.015) & (A["dist"] < 0.005), "S0": A["s0"] == 1}
    tab = {}
    for name, sel in sels.items():
        line = {}
        for key in ("fGeo", "fGeo6", "fGeoS", "fGeoPkg", "fGeoFlat", "fThr80", "fThr100", "fThr130", "fLow", "fBand", "fHigh",
                    "fBandV", "fHighV", "fTrk", "fBack", "pLow", "pBand", "pHigh", "pTrk", "found", "foundPass", "tagW", "tagTrue",
                    "eBelow100", "eBelow100Acc"):
            line[key] = wfrac(key, sel)
        line["eMean"] = wfrac("eMean", sel)[0]
        wa = A["w"][sel] * A["fGeo"][sel]
        line["eMeanAcc"] = (wa * A["eMeanAcc"][sel]).sum() / wa.sum()
        line["rows"] = int(sel.sum()); line["sumw"] = float(A["w"][sel].sum())
        tab[name] = line
        print(f"{name:6s} rows {line['rows']:5d} sum w {line['sumw']:9.1f} | in active module {100*line['fGeo'][0]:5.2f} +- {100*line['fGeo'][1]:.2f} "
              f"(six {100*line['fGeo6'][0]:.2f}, strict {100*line['fGeoS'][0]:.2f}, pkg SMASH {100*line['fGeoPkg'][0]:.2f}, flat {100*line['fGeoFlat'][0]:.2f}) | "
              f"module & E>80/100/130: {100*line['fThr80'][0]:.1f}/{100*line['fThr100'][0]:.1f}/{100*line['fThr130'][0]:.1f} | "
              f"found {100*line['found'][0]:.1f} +- {100*line['found'][1]:.1f} (+ quality {100*line['foundPass'][0]:.1f}) | tag(any, m5c win) {100*line['tagW'][0]:.1f} | true-photon tag "
              f"{100*line['tagTrue'][0]:.2f} +- {100*line['tagTrue'][1]:.2f}")
        print(f"        theta(direction) <11.5 {100*line['fLow'][0]:.1f}, 11.5-47.5 {100*line['fBand'][0]:.1f}, >47.5 {100*line['fHigh'][0]:.1f}; "
              f"coverage band from the vertex 10.9-45.6: {100*line['fBandV'][0]:.1f}, above 45.6: {100*line['fHighV'][0]:.1f}; 18-85: {100*line['fTrk'][0]:.1f}; "
              f"<E> {line['eMean']:.0f} MeV (in modules {line['eMeanAcc']:.0f}); E<100: {100*line['eBelow100'][0]:.0f} % (in modules {100*line['eBelow100Acc'][0]:.0f} %)")
        print(f"        backward (theta > 90): {100*line['fBack'][0]:.1f} %; package definition (bins of 2 deg, normalised to theta < 90): "
              f"<12 {100*line['pLow'][0]:.1f}, 12-48 {100*line['pBand'][0]:.1f}, 48-90 {100*line['pHigh'][0]:.1f}, 18-86 {100*line['pTrk'][0]:.1f}")
    out["table"] = tab
    # versus M(pipi) and distance
    print("versus M(pipi): slice | in module | module & E>100 | found | true-photon tag")
    for lo, hi in ((0, .34), (.34, .37), (.37, .40), (.40, .43)):
        s = (A["mpp"] >= lo) & (A["mpp"] < hi)
        print(f"   {lo:.2f}-{hi:.2f} w {A['w'][s].sum():8.0f} | {100*wfrac('fGeo', s)[0]:5.1f} | {100*wfrac('fThr100', s)[0]:5.1f} | "
              f"{100*wfrac('found', s)[0]:5.1f} | {100*wfrac('tagTrue', s)[0]:5.1f}")
    print("versus the distance to the edge [MeV]: bin | weight | in module | module & E>100 | true-photon tag | <E>")
    dtab = []
    for lo, hi in ((-200, -100), (-100, -60), (-60, -40), (-40, -30), (-30, -20), (-20, -10), (-10, 0), (0, 20)):
        s = (A["dist"] * 1e3 >= lo) & (A["dist"] * 1e3 < hi)
        r_ = (lo, hi, A["w"][s].sum(), wfrac("fGeo", s)[0], wfrac("fThr100", s)[0], wfrac("tagTrue", s)[0], wfrac("eMean", s)[0])
        dtab.append(r_)
        print(f"   {lo:5d}..{hi:4d} {r_[2]:8.0f} | {100*r_[3]:5.1f} | {100*r_[4]:5.1f} | {100*r_[5]:5.1f} | {r_[6]:4.0f}")
    out["dist_table"] = dtab
    # share of the true-photon tags by distance (where does the tagged sample sit?)
    wt = A["w"] * A["tagTrue"]
    for lo, hi in ((-1, -0.060), (-0.060, -0.020), (-0.020, 1)):
        s = (A["dist"] >= lo) & (A["dist"] < hi)
        print(f"   share of the true-photon tagged weight with distance in [{lo*1e3:.0f}, {hi*1e3:.0f}) MeV: {100*wt[s].sum()/wt.sum():.1f} %  "
              f"(all bump weight there: {100*A['w'][s].sum()/A['w'].sum():.1f} %)")
    # closure of the found photons
    fd = A["found"] == 1
    cr = np.abs(A["cRes"][fd]); er = A["eRes"][fd]
    print(f"closure: {fd.sum()} found photons; |cos residual| < 1e-3: {100*(cr < 1e-3).mean():.1f} %; energy residual 5/50/95 %: "
          f"{np.quantile(er, .05):.3f} / {np.quantile(er, .5):.3f} / {np.quantile(er, .95):.3f} MeV; in an active module of the own geometry: "
          f"{100*A['tGeo'][fd].mean():.2f} %, in any module (six sectors): {100*A['tGeo6'][fd].mean():.2f} %")
    # calibration of the predicted probability
    print("found share / predicted in-module probability (rows binned in the prediction):")
    for lo, hi in ((0, .05), (.05, .2), (.2, .4), (.4, .6), (.6, .8), (.8, 1.01)):
        s = (A["fGeo"] >= lo) & (A["fGeo"] < hi)
        if s.sum() == 0:
            continue
        nf = A["found"][s].sum(); pr = A["fGeo"][s].sum()
        print(f"   {lo:.2f}-{hi:.2f}: rows {s.sum():5d}, found {nf:4.0f}, predicted {pr:7.1f}, ratio {nf/pr:.3f} +- {np.sqrt(max(nf,1))/pr:.3f}")
    nm = A["nmult"]
    print(f"neutral multiplicity of the C2+C3 rows: >= 5 (PAT keeps at most 5): {100*(nm >= 5).mean():.1f} % of the rows; mean {nm.mean():.2f}")
    eff = wfrac("found", allm)[0] / wfrac("fGeo", allm)[0]
    print(f"efficiency inside the modules (weighted found / weighted in-module): {eff:.3f}")
    out["eff_in_module"] = eff
    # efficiency inside the modules versus the true energy, from the found photons and the predicted energy spectrum
    # C2 gamma* factor
    c2 = A["cls"] == 2
    for name, s in (("all C2", c2), ("C2 window B", c2 & sels["B"]), ("C2 edge -15..+5", c2 & sels["edge"]),
                    ("C2 far below (<-60)", c2 & (A["dist"] < -0.060)), ("C2 true-photon tagged", c2 & (A["tagTrue"] == 1))):
        f_, e_ = wfrac("qmean", s)
        print(f"C2 gamma* factor (circle average of the QED weight), {name}: {f_:.4f} +- {e_:.4f}")
        out[f"qed_{name}"] = (f_, e_)
    # the found C2 photons: exact per-row QED weight vs the circle average
    json.dump(out, open(O + "/vfrac_numbers.json", "w"), indent=1, default=float)

    # ---------------- peak classes C4 / C5: the e+e- pair (converted / virtual photon) is MEASURED
    print("\npeak classes: direction of the e+e- pair (C4: the converted real photon; C5: the virtual photon) from the true vertex")
    for cls in (4, 5):
        ws, act, act6, tag, low, band, high = [], [], [], [], [], [], []
        Vs, ns = [], []
        for h, cp, sc, nu in rows:
            if int(h[2]) != cls:
                continue
            q = np.array(sc[0:3]) + np.array(sc[3:6])
            Vs.append(sc[12:15]); ns.append(q / np.linalg.norm(q)); ws.append(h[4])
            gidx, ge, gth, gph, meeg, m5 = cp
            mee, m4 = h[8], h[5]
            t_ = gidx > 0 and mee < 0.1 and 0.119 <= meeg < 0.173
            m5c = -1
            if t_ and meeg * meeg - mee * mee > 0:
                ff = (KMPI0 ** 2 - mee * mee) / (meeg * meeg - mee * mee); m2 = m4 * m4 + ff * (m5 * m5 - m4 * m4)
                m5c = np.sqrt(m2) if m2 > 0 else -1
            tag.append(1.0 if (t_ and 0.49 <= m5c < 0.62) else 0.0)
        ws = np.array(ws); Vs = np.array(Vs); ns = np.array(ns); tag = np.array(tag)
        _, _, _, a = geo.hit(Vs, ns); _, _, _, a6 = geo.hit(Vs, ns, six=True)
        th = np.degrees(np.arccos(ns[:, 2]))
        W = ws.sum()
        print(f"  C{cls}: rows used {len(ws)} (sum w {W*(c5every if cls == 5 else 1):.0f}); pair direction in an active module {100*(ws*a).sum()/W:.1f} %, "
              f"six sectors {100*(ws*a6).sum()/W:.1f} %; theta <11.5 {100*(ws*(th<11.5)).sum()/W:.1f} %, 11.5-47.5 {100*(ws*((th>=11.5)&(th<=47.5))).sum()/W:.1f} %, "
              f">47.5 {100*(ws*(th>47.5)).sum()/W:.1f} %; ECAL tag (any photon, m5c window) {100*(ws*tag).sum()/W:.2f} %")
        out[f"peak_C{cls}"] = dict(active=float((ws * a).sum() / W), six=float((ws * a6).sum() / W), tag=float((ws * tag).sum() / W))
    json.dump(out, open(O + "/vfrac_numbers.json", "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
