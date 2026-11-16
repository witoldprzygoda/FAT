#!/usr/bin/env python3
# vecalacc_tag.py -- INDEPENDENT VERIFIER of package ecalacc (2026-11-15): the ECAL tag of the simulated bump rows
# recomputed from the PAT neutral candidates with the cluster energies x lambda (own code; FAT logic as documented in
# the package: quality cuts PID 1, 0.8 < beta < 1.2, lambda E > 100 MeV; the three most energetic; the one with
# M(e+e-gamma) closest to 0.140; tag if M(ee) < 0.10 and 0.119 <= M(e+e-gamma) < 0.173; m5c in 0.49-0.62).
# lambda = 1 is checked against the compact's own photon columns.
import sys
import numpy as np
sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
from vecalacc_frac import read_rows, four, mdot, ME, MPI, KMPI0

O = "/hdd3/przygoda/HADES/FAT_out/all_report/vecalacc"
LAM = (1.0, 1.10, 1.20, 1.285)


def p4(p, th, ph, m):
    t, f = np.radians(th), np.radians(ph)
    return four(p * np.sin(t) * np.cos(f), p * np.sin(t) * np.sin(f), p * np.cos(t), m)


def main():
    rows = []
    with open(O + "/vrows.txt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            p = line.split("|")
            h = list(map(float, p[0].split()))
            if int(h[2]) not in (2, 3):
                continue
            rows.append((h, list(map(float, p[1].split())), list(map(float, p[2].split())),
                         np.array(list(map(float, p[3].split()))).reshape(18, 5), list(map(float, p[4].split()))))
    print("C2+C3 rows:", len(rows))
    W = sum(r[0][4] for r in rows)
    sumT = {l: 0.0 for l in LAM}; sumTany = {l: 0.0 for l in LAM}
    nagree = 0; edge_like = {l: 0.0 for l in LAM}
    for h, cp, sc, nu, rk in rows:
        w, m4, mee = h[4], h[5], h[8]
        ep, em = p4(rk[0], rk[1], rk[2], ME), p4(rk[3], rk[4], rk[5], ME)
        pp, pm = p4(rk[6], rk[7], rk[8], MPI), p4(rk[9], rk[10], rk[11], MPI)
        ee = ep + em; v4 = ee + pp + pm
        n = min(5, int(round(sc[23]))) if sc[23] > 0 else 0
        # matched true photon (same definition as in vecalacc_frac.py)
        km, bE = -1, -1
        for k in range(n):
            if int(round(nu[7, k])) == 1 and int(round(nu[8, k])) == 7 and nu[9, k] == sc[18] and nu[10, k] == sc[19] and nu[11, k] == sc[20] \
               and np.sqrt((nu[15, k] - sc[12]) ** 2 + (nu[16, k] - sc[13]) ** 2 + (nu[17, k] - sc[14]) ** 2) <= 0.01 and nu[2, k] > bE:
                km, bE = k, nu[2, k]
        for l in LAM:
            ok = [k for k in range(n) if int(round(nu[0, k])) == 1 and 0.8 < nu[1, k] < 1.2 and nu[2, k] * l > 100.]
            ok = sorted(ok, key=lambda k: -nu[2, k])[:3]
            best, kb, meeg, m5 = 1e9, -1, -1, -1
            for k in ok:
                g = p4(nu[2, k] * l, nu[3, k], nu[4, k], 0.0)
                t = ee + g; m = np.sqrt(max(mdot(t, t), 0)) / 1e3
                if abs(m - 0.140) < best:
                    t5 = v4 + g
                    best, kb, meeg, m5 = abs(m - 0.140), k, m, np.sqrt(max(mdot(t5, t5), 0)) / 1e3
            tag = kb >= 0 and mee < 0.1 and 0.119 <= meeg < 0.173
            m5c = -1
            if tag and meeg * meeg - mee * mee > 0:
                f = (KMPI0 ** 2 - mee * mee) / (meeg * meeg - mee * mee); m2 = m4 * m4 + f * (m5 * m5 - m4 * m4)
                m5c = np.sqrt(m2) if m2 > 0 else -1
            inw = tag and 0.49 <= m5c < 0.62
            if inw:
                sumTany[l] += w
                if kb == km:
                    sumT[l] += w
            if l == 1.0:
                # compare with the compact's photon columns
                gidx, ge, gth, gph, meeg_c, m5_c = cp
                tag_c = gidx > 0 and mee < 0.1 and 0.119 <= meeg_c < 0.173
                same = (tag == tag_c) and (not tag or abs(meeg - meeg_c) < 2e-5)
                nagree += same
    print(f"lambda = 1: recomputed tag agrees with the compact for {nagree} of {len(rows)} rows")
    for l in LAM:
        print(f"lambda {l:.3f}: true-photon tag share {100*sumT[l]/W:.2f} %, any-photon tag share {100*sumTany[l]/W:.2f} %; ratio to lambda 1: "
              f"true {sumT[l]/sumT[1.0]:.3f}, any {sumTany[l]/sumTany[1.0]:.3f}")


if __name__ == "__main__":
    main()
