#!/usr/bin/env python3
# ============================================================================
# ecalacc_geo.py -- package ecalacc (2026-11-15): ECAL coverage from the geometry of the simulation
# ============================================================================
# Input: /hdd3/przygoda/HADES/FAT_out/all_report/ecalacc/emc_cells_sim_v2.txt (EmcGeomPar + EmcCellStatusPar of the gen4
#        simulation parameter file, written by ecalacc_emcgeo_vae25.C on vae25) and emcstat_real_vae25_stdout_v2.txt
#        (the same status in the data parameter file for eight runs).
# Output: docs/latex/tables_all/ecalacc_geometry.{txt,tex}; /hdd3/.../ecalacc/geo_numbers.json; geo_cells.json (for
#         the figures).
# Run: python3 -I plots/all_report/ecalacc_geo.py
# ============================================================================
import json
import math
import numpy as np

O = "/hdd3/przygoda/HADES/FAT_out/all_report/ecalacc"
T = "/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all"

rows = [l.split() for l in open(O + "/emc_cells_sim_v2.txt") if l.startswith("CELL")]
cells = []
for r in rows:
    s, c, st = int(r[1]), int(r[2]), int(r[3])
    pts = np.array([[float(r[9 + 3 * i + j]) for j in range(3)] for i in range(8)])
    d = np.linalg.norm(pts, axis=1)
    o = np.argsort(d)
    cells.append(dict(sec=s, cell=c, status=st, front=pts[o[:4]], back=pts[o[4:]], centre=np.array([float(r[6]), float(r[7]), float(r[8])])))


def thph(p, vz=0.0):
    x, y, z = p[..., 0], p[..., 1], p[..., 2] - vz
    th = np.degrees(np.arctan2(np.hypot(x, y), z))
    ph = np.degrees(np.arctan2(y, x)) % 360.0
    return th, ph


# sector planes
sec = {}
for s in range(6):
    f = np.vstack([c["front"] for c in cells if c["sec"] == s])
    cen = f.mean(0)
    u, sv, vt = np.linalg.svd(f - cen)
    n = vt[2] if np.dot(vt[2], cen) > 0 else -vt[2]
    dist = float(np.dot(cen, n))
    rms = float(np.sqrt(np.mean(((f - cen) @ n) ** 2)))
    sec[s] = dict(cen=cen, n=n, dist=dist, rms=rms)

# per-sector table
out = []
for s in range(6):
    cs = [c for c in cells if c["sec"] == s]
    act = [c for c in cs if c["status"] == 1]
    res = dict(sector=s, n=len(cs), active=len(act), dead=[c["cell"] for c in cs if c["status"] != 1] if len(act) > 0 else "all")
    for vz in (0.0, -60.0, -120.0):
        ths, phs = [], []
        for c in cs:
            th, ph = thph(c["front"], vz)
            ths += list(th); phs += list(ph)
        res["th_%d" % int(-vz)] = (min(ths), max(ths))
        phs = np.array(phs)
        # phi range around the sector centre
        pc = (90.0 + 60.0 * s) % 360.0
        dph = (phs - pc + 180.0) % 360.0 - 180.0
        res["ph_%d" % int(-vz)] = (pc + dph.min(), pc + dph.max(), dph.max() - dph.min())
    res["plane_dist_mm"] = sec[s]["dist"]
    res["plane_rms_mm"] = sec[s]["rms"]
    res["tilt_deg"] = math.degrees(math.acos(abs(np.dot(sec[s]["n"], sec[s]["cen"] / np.linalg.norm(sec[s]["cen"])))))
    out.append(res)

# module size and pitch
f0 = cells[0]["front"]
edges = sorted(np.linalg.norm(f0[i] - f0[j]) for i in range(4) for j in range(i + 1, 4))
cent = np.array([c["front"].mean(0) for c in cells if c["sec"] == 2])
dd = np.sort(np.linalg.norm(cent[:, None, :] - cent[None, :, :], axis=2), axis=1)[:, 1]
depth = float(np.linalg.norm(cells[0]["back"].mean(0) - cells[0]["front"].mean(0)))

# solid angle of the active modules seen from z = 0 and z = -60 mm (isotropic rays, Monte Carlo)
rng = np.random.default_rng(12345)
N = 4_000_000
cth = rng.uniform(-1, 1, N); phi = rng.uniform(0, 2 * np.pi, N)
sth = np.sqrt(1 - cth ** 2)
dirs = np.stack([sth * np.cos(phi), sth * np.sin(phi), cth], 1)


def covered(vz, active_only=True, enlarge=94.0 / 92.0):
    V = np.array([0, 0, vz])
    hit = np.zeros(N, bool)
    phd = np.degrees(phi) % 360
    for s in range(6):
        x = (phd - 60.0) % 360.0
        m = (np.floor(x / 60.0).astype(int) == s)
        S = sec[s]; n = S["n"]; c0 = S["cen"]
        den = dirs[m] @ n
        t = (np.dot(c0 - V, n)) / np.where(den > 1e-9, den, np.nan)
        X = V + t[:, None] * dirs[m] - c0
        a = np.array([0, 0, 1.0]) - n * n[2]; a /= np.linalg.norm(a); b = np.cross(n, a)
        u, v = X @ a, X @ b
        h = np.zeros(m.sum(), bool)
        for c in cells:
            if c["sec"] != s or (active_only and c["status"] != 1):
                continue
            fu = (c["front"] - c0) @ a; fv = (c["front"] - c0) @ b
            uc, vc = fu.mean(), fv.mean()
            o = np.argsort(np.arctan2(fv - vc, fu - uc))
            qu = uc + (fu[o] - uc) * enlarge; qv = vc + (fv[o] - vc) * enlarge
            inside = np.ones(len(u), bool); pos = np.zeros(len(u), bool); neg = np.zeros(len(u), bool)
            for i in range(4):
                j = (i + 1) % 4
                cr = (qu[j] - qu[i]) * (v - qv[i]) - (qv[j] - qv[i]) * (u - qu[i])
                pos |= cr > 0; neg |= cr < 0
            h |= ~(pos & neg)
        h &= np.isfinite(t) & (t > 0)
        hit[np.where(m)[0]] = h
    return hit


res_sa = {}
for vz in (0.0, -60.0):
    h = covered(vz)
    h6 = covered(vz, active_only=False)
    th = np.degrees(np.arccos(cth))
    band = (th > 11.5) & (th < 47.6)
    res_sa[str(int(vz))] = dict(frac4pi=float(h.mean()), frac4pi_6sec=float(h6.mean()), frac_band=float(h[band].mean()),
                                band4pi=float(band.mean()))

# data status (vae25, real parameters)
real = [l.split() for l in open(O + "/emcstat_real_vae25_stdout_v2.txt") if l.startswith("RUN")]
dead_real = sorted((int(l.split()[2]), int(l.split()[3])) for l in open(O + "/emcstat_real_vae25_stdout_v2.txt") if l.startswith("DEAD"))
dead_sim = sorted((c["sec"], c["cell"]) for c in cells if c["status"] != 1)
same = dead_real == dead_sim

num = dict(sectors=out, module_edges_mm=[float(e) for e in edges], pitch_mm=float(np.median(dd)), depth_mm=depth,
           solid_angle=res_sa, real_runs=[" ".join(r) for r in real], dead_sim=len(dead_sim), dead_real=len(dead_real),
           status_identical=same, n_cells=len(cells), n_active=sum(c["status"] == 1 for c in cells))
json.dump(num, open(O + "/geo_numbers.json", "w"), indent=1, default=lambda x: x.tolist() if hasattr(x, "tolist") else str(x))

# cells for the figure (theta/phi corners from z = 0 and z = -60)
cj = []
for c in cells:
    th0, ph0 = thph(c["front"], 0.0)
    cj.append(dict(sec=c["sec"], cell=c["cell"], status=c["status"], th=list(map(float, th0)), ph=list(map(float, ph0))))
json.dump(cj, open(O + "/geo_cells.json", "w"))

# ---- tables
L = []
L.append("ECAL coverage in the simulation (gen4 parameter file allParam_feb22_sim_run_17000_gen4_10032026.root, EmcGeomPar + EmcCellStatusPar)")
L.append("module front face %.1f x %.1f mm, pitch %.1f mm, depth %.1f mm; per sector one plane (rms %.4f mm); %d modules, %d active" % (
    edges[0], edges[1], float(np.median(dd)), depth, max(s["plane_rms_mm"] for s in out), len(cells), num["n_active"]))
L.append("sector | modules | active | theta range from z = 0 / -60 / -120 mm [deg] | phi range from z = 0 [deg] (width) | plane distance [mm] | tilt [deg]")
for s in out:
    L.append("%d | %d | %d | %.2f-%.2f / %.2f-%.2f / %.2f-%.2f | %.2f-%.2f (%.2f) | %.1f | %.2f" % (
        s["sector"], s["n"], s["active"], *s["th_0"], *s["th_60"], *s["th_120"], *s["ph_0"], s["plane_dist_mm"], s["tilt_deg"]))
L.append("status-0 modules (simulation): %s" % ", ".join("s%d/c%d" % k for k in dead_sim if k[0] != 0) + "; and all 163 of sector 0")
L.append("data (real) parameters, 8 runs of the four periods: identical status map: %s" % same)
for k, v in res_sa.items():
    L.append("solid angle of the active modules seen from z = %s mm: %.4f of 4 pi (all six sectors, any status: %.4f); band 11.5-47.6 deg = %.4f of 4 pi, of which %.3f covered" % (
        k, v["frac4pi"], v["frac4pi_6sec"], v["band4pi"], v["frac_band"]))
open(T + "/ecalacc_geometry.txt", "w").write("\n".join(L) + "\n")
print("\n".join(L))

tex = [r"\begin{tabular}{@{}c c c c c c@{}}", r"\toprule",
       r"sector & modules & active & $\theta$ from $z=0$ [deg] & $\theta$ from $z=-60$\,mm [deg] & $\phi$ from $z=0$ [deg] \\", r"\midrule"]
for s in out:
    tex.append(r"%d & %d & %d & %.1f--%.1f & %.1f--%.1f & %.1f--%.1f \\" % (s["sector"], s["n"], s["active"], *s["th_0"], *s["th_60"], s["ph_0"][0], s["ph_0"][1]))
tex.append(r"\midrule")
v0 = res_sa["0"]
tex.append(r"all & %d & %d & \multicolumn{3}{l}{solid angle: %.1f\,\%% of $4\pi$ (six full sectors: %.1f\,\%%)} \\" % (
    len(cells), num["n_active"], 100 * v0["frac4pi"], 100 * v0["frac4pi_6sec"]))
tex.append(r"\bottomrule")
tex.append(r"\end{tabular}")
open(T + "/ecalacc_geometry.tex", "w").write("\n".join(tex) + "\n")
