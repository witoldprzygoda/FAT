#!/usr/bin/env python3
"""vpt_nom.py -- verifier of peaktail: machinery check (own code reproduces the edgeu nominal) and the shared spectra.

Writes /hdd3/przygoda/HADES/FAT_out/all_report/vpt/vpt_spectra.npz (data N, N++, N-- per slice in u; sim C2c, C3, C4, C5
per row index lists are recomputed by the other scripts).
"""
import sys
import numpy as np
sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import vpt_core as V

SL = [(0.30, 0.34), (0.34, 0.37), (0.37, 0.39), (0.39, 0.40)]
g = V.Grid()


def in_sl(m, a, b):
    m = m.astype(np.float32)
    return (m >= np.float32(a)) & (m < np.float32(b))


def data_u(D):
    u = D["m4"].astype(float) - V.edge(D["mpipi"].astype(float), D["mee"].astype(float))
    u[~(D["mee"] < V.MPI0)] = np.nan
    return u


def data_spec(D, u, sel, slices=SL, extra=None, mode=0):
    """N - CB summed over slices, CB per slice"""
    N = np.zeros(g.nb); CB = np.zeros(g.nb); COV = np.zeros((g.nb, g.nb))
    for a, b in slices:
        m = sel & in_sl(D["mpipi"], a, b)
        if extra is not None:
            m &= extra
        n = g.h(u[m & (D["k"] == 0)]); p = g.h(u[m & (D["k"] == 1)]); q = g.h(u[m & (D["k"] == 2)])
        c, cv = V.cb_cells(p, q, g, mode)
        N += n; CB += c; COV += cv
    return N, CB, COV


def sim_u(S):
    m4p = V.sim_m4p(S); mpp = V.SPP * S["pp"].astype(float)
    u = m4p - V.edge(mpp, S["ee"].astype(float))
    u[~(S["ee"] < V.MPI0)] = np.nan
    return u, mpp, m4p


if __name__ == "__main__":
    D = V.load_data(); S = V.load_sim()
    u = data_u(D); S0 = (D["bits"] & 1) != 0
    N, CB, COV = data_spec(D, u, S0)
    us, mpp, m4p = sim_u(S)
    w = V.sim_weight(S)
    insl = np.zeros(len(us), bool)
    for a, b in SL:
        insl |= in_sl(mpp, a, b)
    y = {}
    for k, gi in (("C2", 0), ("C3", 1), ("C4", 2), ("C5", 3)):
        m = insl & (S["g"] == gi)
        y[k] = g.h(us[m], w[m])
    T = V.step_t(g, V.SIG_STEP)
    K = y["C4"] + y["C5"]
    R = V.fitstep(N - CB, CB, COV, T, g, fixed=V.AP * K)
    hs = V.hlin(R, y["C2"] + y["C3"])
    print("data rows S0 us/pp/mm in slices:", [int(((D["k"] == k) & S0 & np.isfinite(u)).sum()) for k in range(3)])
    print("nominal: h_data %.4f +- %.4f chi2 %.2f/%d  h_sim %.4f  f %.5f   (edgeu: 1305.2164 +- 62.2741, 27.34/32, 3107.3304, 0.42004)"
          % (R.h, R.eh, R.chi2, R.ndf, hs, R.h / hs))
    print("lin: C2c %.2f C3 %.2f C4 %.2f C5 %.2f  (edgeu C2c 1780.31 C3 1327.02 C4 40.95 C5 -37.11)" %
          tuple(V.hlin(R, y[k]) for k in ("C2", "C3", "C4", "C5")))
    print("sim sums in slices, all u: C2c %.1f C3 %.1f C4 %.1f C5 %.1f" % tuple(y[k].sum() for k in ("C2", "C3", "C4", "C5")))
    np.savez(V.OUT + "vpt_nom.npz", N=N, CB=CB, COV=COV, C2=y["C2"], C3=y["C3"], C4=y["C4"], C5=y["C5"], T=T)
