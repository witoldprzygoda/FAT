#!/usr/bin/env python3
# vecalacc_clusters.py -- INDEPENDENT VERIFIER of package ecalacc (2026-11-15): empirical ECAL acceptance from the
# neutral candidates of the standard-set PAT file (vclusters_std.txt of vecalacc_rows.C; one event written once).
#   1. photons (GEANT id 1) that made a neutral candidate: polar angle of the TRUE direction from the creation vertex;
#   2. module occupancy: which (sector, cell) carry clusters (sector 0 and the 8 dead modules must be empty);
#   3. position check: the module hit by the true photon direction in the own geometry (vecalacc_geo.py) against the
#      cluster's (sector, cell) -- same module or a direct neighbour.
import sys
import numpy as np
sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import vecalacc_geo as G

O = "/hdd3/przygoda/HADES/FAT_out/all_report/vecalacc"
a = np.loadtxt(O + "/vclusters_std.txt", comments="#")
# columns: entry k pid beta E theta phi sector cell sim_id parentpid px py pz vx vy vz
pid, beta, E, thc, phc, sec, cell, sid, ppid = a[:, 2], a[:, 3], a[:, 4], a[:, 5], a[:, 6], a[:, 7].astype(int), a[:, 8].astype(int), a[:, 9].astype(int), a[:, 10].astype(int)
p = a[:, 11:14]; v = a[:, 14:17]
print(f"neutral candidates: {len(a)}; with GEANT id 1 (photon): {(sid == 1).sum()}")
cells, raw = G.load_geo(); geo = G.Geo(cells)
status = {(c[0], c[1]): c[3] for c in cells}
# 2. occupancy
occ = {}
for s_, c_ in zip(sec, cell):
    occ[(s_, c_)] = occ.get((s_, c_), 0) + 1
n_act_occ = sum(1 for k, n in occ.items() if status.get(k, -1) == 1)
n_dead_occ = sum(n for k, n in occ.items() if status.get(k, -1) == 0)
n_bad = sum(n for k, n in occ.items() if k not in status)
print(f"modules with clusters: {len(occ)}; of them active {n_act_occ} (of 807); clusters in status-0 modules {n_dead_occ}; "
      f"clusters with an unknown (sector, cell) {n_bad}; sector-0 clusters {(sec == 0).sum()}")
low = sorted([(n, k) for k, n in occ.items() if status.get(k) == 1])[:5]
print("least occupied active modules (count, (sector, cell)):", low)
# 1. true photon directions (photons from near the target: |vz + 124| < 40 mm, r < 10 mm)
ph = (sid == 1) & (np.hypot(v[:, 0], v[:, 1]) < 10) & (np.abs(v[:, 2] + 124) < 40)
d = p[ph] / np.linalg.norm(p[ph], axis=1)[:, None]
th = np.degrees(np.arccos(d[:, 2]))
print(f"target photons with a cluster: {ph.sum()}; true polar angle 0.1/1/50/99/99.9 %: " +
      " / ".join(f"{np.quantile(th, q):.2f}" for q in (0.001, 0.01, 0.5, 0.99, 0.999)) + " deg")
# compare with the geometric limit from each photon's own vertex
s_h, r_h, c_h, act = geo.hit(v[ph], d)
_, _, _, any6 = geo.hit(v[ph], d, six=True)
print(f"their true direction hits an active module of the own geometry: {100*act.mean():.2f} %; any module: {100*any6.mean():.2f} %")
# 3. same module or neighbour
cs, cc = sec[ph], cell[ph]
same = (s_h == cs) & (r_h * 17 + c_h == cc)
nb = (s_h == cs) & (np.abs(r_h - cc // 17) <= 1) & (np.abs(c_h - cc % 17) <= 1)
print(f"cluster module = module hit by the true direction: {100*same.mean():.1f} %; same or direct neighbour: {100*nb.mean():.1f} %")
# theta of the cluster (PAT value) vs the module-centre theta from z = 0
print(f"cluster theta (PAT) range: {thc[thc > 0].min():.2f} - {thc.max():.2f} deg")
