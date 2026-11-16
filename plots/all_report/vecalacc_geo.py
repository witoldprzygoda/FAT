#!/usr/bin/env python3
# vecalacc_geo.py -- INDEPENDENT VERIFIER of package ecalacc (2026-11-15): ECAL coverage from the own geometry dump
# Input : /hdd3/przygoda/HADES/FAT_out/all_report/vecalacc/vemcgeo_stdout.txt  (vecalacc_emcgeo_vae25.C)
#         /hdd3/przygoda/HADES/FAT_out/all_report/vecalacc/vemcstat_stdout.txt (vecalacc_emcstat_vae25.C)
# Output: /hdd3/przygoda/HADES/FAT_out/all_report/vecalacc/vgeo.npz (grid description used by vecalacc_frac.py) and
#         vgeo_numbers.json; printout.
# Method (deliberately different from ecalacc_geo.py / ecalacc_fill.C):
#   - per sector a plane fitted (SVD) to all front-face corners; the modules are put on a regular (row, column) grid of
#     that plane (pitch measured from the centres); a ray hits module (row, col) if its plane point is within
#     +-pitch/2 (or +-46 mm strict) of that centre in both in-plane coordinates -- no polygon test;
#   - the sector of a ray is not taken from phi: all six planes are tried and the hit with the smallest positive
#     distance wins;
#   - solid angle by a Fibonacci-sphere quadrature (deterministic), not random rays.
import json, sys
import numpy as np

O = "/hdd3/przygoda/HADES/FAT_out/all_report/vecalacc"


def load_geo(fn=O + "/vemcgeo_stdout.txt"):
    cells = []
    raw_sim = None
    for line in open(fn):
        if line.startswith("VCELL"):
            a = line.split()
            s, c, pos, st = int(a[1]), int(a[2]), int(a[3]), int(a[4])
            ctr = np.array(list(map(float, a[5:8])))
            pts = np.array(list(map(float, a[8:8 + 24]))).reshape(8, 3)
            cells.append((s, c, pos, st, ctr, pts[:4], pts[4:]))
        elif line.startswith("RAWSTATUS_SIM"):
            raw_sim = np.array(list(map(int, line.split()[1:])))
    return cells, raw_sim


def load_real(fn=O + "/vemcstat_stdout.txt"):
    out = {}
    for line in open(fn):
        if line.startswith("RAWSTATUS_REAL"):
            a = line.split()
            out[int(a[1])] = np.array(list(map(int, a[2:])))
    return out


class Geo:
    def __init__(self, cells):
        self.sec = []
        for s in range(6):
            cs = [c for c in cells if c[0] == s]
            front = np.vstack([c[5] for c in cs])
            cen = front.mean(axis=0)
            u, sv, vt = np.linalg.svd(front - cen)
            n = vt[2]
            if n.dot(cen) < 0:
                n = -n
            # in-plane axes: a along the module rows (from the module grid), b = n x a
            ctrs = np.array([c[5].mean(axis=0) for c in cs])     # front-face centres
            # the columns of one row (same cell // 17) define the a axis
            rows = np.array([c[1] // 17 for c in cs]); cols = np.array([c[1] % 17 for c in cs])
            r0 = np.bincount(rows).argmax()
            sel = rows == r0
            k = np.argsort(cols[sel])
            d = ctrs[sel][k][-1] - ctrs[sel][k][0]
            a = d - n * n.dot(d); a /= np.linalg.norm(a)
            b = np.cross(n, a)
            uc = (ctrs - cen) @ a; vc = (ctrs - cen) @ b
            # regular grid: u = u0 + pitch * col, v = v0 + pitch * row (fit)
            A = np.vstack([np.ones_like(uc), cols]).T
            pu, *_ = np.linalg.lstsq(A, uc, rcond=None)
            B = np.vstack([np.ones_like(vc), rows]).T
            pv, *_ = np.linalg.lstsq(B, vc, rcond=None)
            res = max(np.abs(A @ pu - uc).max(), np.abs(B @ pv - vc).max())
            # face size from the corners
            sz = []
            for c in cs:
                f = c[5]
                sz.append(sorted([np.linalg.norm(f[i] - f[j]) for i in range(4) for j in range(i + 1, 4)])[:4])
            st = {(int(r), int(cc)): int(c[3]) for r, cc, c in zip(rows, cols, cs)}
            self.sec.append(dict(n=n, c=cen, a=a, b=b, u0=pu[0], du=pu[1], v0=pv[0], dv=pv[1], gridres=res,
                                 planerms=sv[2] / np.sqrt(len(front)), st=st, face=np.mean(sz), dist=n.dot(cen),
                                 ncell=len(cs), nact=sum(1 for c in cs if c[3] == 1)))

    def hit(self, V, d, half=None, six=False):
        """V (N,3) vertex mm, d (N,3) unit directions. Returns sector (or -1), row, col, active flag."""
        N = len(d)
        best_t = np.full(N, np.inf); sec = np.full(N, -1); row = np.zeros(N, int); col = np.zeros(N, int)
        act = np.zeros(N, bool)
        for s, S in enumerate(self.sec):
            h = (0.5 * S["du"]) if half is None else half
            den = d @ S["n"]
            with np.errstate(divide="ignore", invalid="ignore"):
                t = ((S["c"] - V) @ S["n"]) / den
            ok = (den > 1e-9) & (t > 0)
            X = V + t[:, None] * d - S["c"]
            u = X @ S["a"]; v = X @ S["b"]
            cc = np.rint((u - S["u0"]) / S["du"]).astype(int); rr = np.rint((v - S["v0"]) / S["dv"]).astype(int)
            uu = S["u0"] + S["du"] * cc; vv = S["v0"] + S["dv"] * rr
            inside = ok & (np.abs(u - uu) <= h) & (np.abs(v - vv) <= h)
            stv = np.array([S["st"].get((int(r_), int(c_)), -1) for r_, c_ in zip(rr[inside], cc[inside])])
            exists = np.zeros(N, bool); active = np.zeros(N, bool)
            idx = np.where(inside)[0]
            exists[idx] = stv >= 0
            active[idx] = stv == 1
            take = exists & (t < best_t)
            best_t[take] = t[take]; sec[take] = s; row[take] = rr[take]; col[take] = cc[take]
            act[take] = active[take] | (six & exists[take])
        return sec, row, col, act


def fib_sphere(n):
    i = np.arange(n) + 0.5
    z = 1 - 2 * i / n
    ph = np.pi * (1 + 5 ** 0.5) * i
    r = np.sqrt(1 - z * z)
    return np.vstack([r * np.cos(ph), r * np.sin(ph), z]).T


def main():
    cells, raw_sim = load_geo()
    real = load_real()
    g = Geo(cells)
    out = {}
    print("cells:", len(cells))
    for s, S in enumerate(g.sec):
        print(f"sector {s}: {S['ncell']} modules, {S['nact']} active; plane rms {S['planerms']:.4f} mm, grid residual "
              f"{S['gridres']:.4f} mm, pitch {S['du']:.3f} / {S['dv']:.3f} mm, face side {S['face']:.2f} mm, plane distance {S['dist']:.1f} mm")
    # status: sim vs real, and the cell -> position mapping
    for run, arr in real.items():
        print(f"real status run {run}: zeros {int((arr == 0).sum())}, identical to sim: {bool((arr == raw_sim).all())}")
    dead = [(c[0], c[1]) for c in cells if c[3] != 1 and c[0] != 0]
    print("dead single modules (sector, cell):", dead)
    out["dead_single"] = dead
    out["n_active"] = int(sum(1 for c in cells if c[3] == 1))
    # theta / phi ranges of the front-face corners from z = 0 and z = -60 mm
    for z0 in (0.0, -60.0):
        ths = []
        for s in range(6):
            fr = np.vstack([c[5] for c in cells if c[0] == s]) - np.array([0, 0, z0])
            th = np.degrees(np.arctan2(np.hypot(fr[:, 0], fr[:, 1]), fr[:, 2]))
            ph = np.degrees(np.arctan2(fr[:, 1], fr[:, 0])) % 360
            ths.append((th.min(), th.max()))
            if z0 == 0.0:
                lo = ph.min() if s != 5 else ph[ph < 180].min()
                hi = ph.max() if s != 5 else ph[ph < 180].max()
                print(f"  sector {s}: theta {th.min():.2f}-{th.max():.2f} deg, phi {lo:.1f}-{hi:.1f} deg (corners, z0 = 0)")
        lo = min(t[0] for t in ths); hi = max(t[1] for t in ths)
        print(f"corner theta range from z0 = {z0:.0f} mm: {lo:.2f} - {hi:.2f} deg")
        out[f"theta_range_z{int(z0)}"] = [lo, hi]
    # where along the polar angle is the coverage complete? theta scan of the covered azimuth fraction (z0 = 0)
    nphi = 7200
    phis = np.radians((np.arange(nphi) + 0.5) * 360 / nphi)
    prof = []
    for th in np.arange(10.0, 50.01, 0.25):
        t = np.radians(th)
        d = np.vstack([np.sin(t) * np.cos(phis), np.sin(t) * np.sin(phis), np.full(nphi, np.cos(t))]).T
        V = np.zeros_like(d)
        _, _, _, act = g.hit(V, d)
        _, _, _, act6 = g.hit(V, d, six=True)
        prof.append((th, act.mean(), act6.mean()))
    prof = np.array(prof)
    out["phi_profile"] = prof.tolist()
    print("azimuthal coverage versus theta (z0 = 0): theta, active, six full sectors")
    for th, a, a6 in prof[::8]:
        print(f"   {th:5.2f}  {a:.3f}  {a6:.3f}")
    nz = prof[prof[:, 2] > 0]
    print(f"   first / last theta with any module: {nz[0, 0]:.2f} / {nz[-1, 0]:.2f} deg; >= 50 % of the six-sector maximum: "
          f"{prof[prof[:, 2] >= 0.5 * prof[:, 2].max(), 0].min():.2f}-{prof[prof[:, 2] >= 0.5 * prof[:, 2].max(), 0].max():.2f} deg")
    # solid angle (Fibonacci sphere, 4e6 points) from z0 = 0 and z0 = -60 mm, pitch faces and strict 92 mm faces
    d = fib_sphere(4_000_000)
    for z0 in (0.0, -60.0):
        V = np.tile([0.0, 0.0, z0], (len(d), 1))
        s1, r1, c1, a1 = g.hit(V, d)
        _, _, _, a6 = g.hit(V, d, six=True)
        _, _, _, aS = g.hit(V, d, half=46.0)
        th = np.degrees(np.arccos(d[:, 2]))
        band = (th >= 11.5) & (th <= 47.5)
        print(f"solid angle from z0 = {z0:.0f} mm: active {a1.mean() * 100:.3f} % of 4 pi, six sectors {a6.mean() * 100:.3f} %, "
              f"strict 92 mm {aS.mean() * 100:.3f} %; active share of the 11.5-47.5 deg band {a1[band].mean() * 100:.2f} %")
        out[f"omega_z{int(z0)}"] = dict(active=a1.mean(), six=a6.mean(), strict=aS.mean(), band_share=a1[band].mean())
    np.savez(O + "/vgeo.npz", cells=np.array([(c[0], c[1], c[3]) for c in cells]))
    json.dump(out, open(O + "/vgeo_numbers.json", "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
