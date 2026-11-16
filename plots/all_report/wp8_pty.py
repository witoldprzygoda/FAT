#!/usr/bin/env python3
# ============================================================================
# wp8_pty.py  --  WP8 (report_all): (pT4, y4) reweighting variant of the eta production (3 x 3 cells)
# ============================================================================
# Written by WP8 (prefix wp8), 2026-11-14.
# Inputs:  /hdd3/przygoda/HADES/FAT_out/all_report/wp8/wp8_grid.txt (wp8_sim.C: data S0 rows in P and the side bands per
#          (pT4, y4) grid cell; simulated C4, C5 in P and C2, C3 in B, S0, V3, in nb), the SMASH standard-set generator
#          histogram of the produced eta (pT, y_cm) (/hdd3/.../wp8/remote/gen_eta_ptY_std40.txt, 40 converted files on vae25, read-only).
# Procedure (plain words):
#   1. Cells: 3 x 3 in (pT4, y4); inner edges = tertiles of the simulated peak signal C4 + C5 in P (rounded to the
#      20 MeV/c x 0.05 grid); the outer cells extend to 0 / infinity in pT and to -+infinity in y, so every produced eta
#      falls into a cell.
#   2. Data peak yield per cell = the peak side-band estimator of WP6 applied cell by cell with the global WP6 constants
#      (R' = 1.1859, leakage f = 0.1544): Y_c = (64/63) [(N - CB)_P,c - R' (N - CB)_SB,c] / (1 - R' f); CB_c = window CB
#      (WP6: 3279.0 in P, 2695.0 in the side bands) x the like-sign share of the cell.
#   3. r_c = Y_c / S_c, S_c = simulated C4 + C5 in P (expected raw rows x 64/63, SMASH sigma x PDG BR, V3).
#   4. Nominal scale (accepted-weighted): a_P = sum Y_c / sum S_c. Reweighted production: the data-like produced
#      distribution in cell c is r_c x G_c (G_c = produced SMASH etas in the cell), so sigma_data / sigma_SMASH =
#      a_P' = sum_c g_c r_c (g_c = produced fraction). The change a_P'/a_P is the (pT4, y4) variant of a_P and sigma_eta.
#   5. Bump (approximation, information): the same weights w_c = r_c / a_P' applied to the simulated C2 + C3 in B with
#      their RECONSTRUCTED (pT4, y4) as a proxy of the eta kinematics (the photon is missing): k_B = sum w_c S_B,c / S_B;
#      peak k_P = a_P/a_P'. DR' = DR x k_B / k_P, a_B' = a_B / k_B.
#   6. Statistical noise of the variant: 2000 Gaussian toys of Y_c (event-group factor 1.19 on the Poisson variance,
#      + CB), S_c (MC event-group variance) and g_c (multinomial, generator statistics).
# Output: docs/latex/tables_all/wp8_pty.txt (+ machine lines 'pty ...'), /hdd3/.../wp8/wp8_pty.json
# Run: python3 -I plots/all_report/wp8_pty.py
# ============================================================================
import json
import math
import numpy as np

GRID = "/hdd3/przygoda/HADES/FAT_out/all_report/wp8/wp8_grid.txt"
GEN = "/hdd3/przygoda/HADES/FAT_out/all_report/wp8/remote/gen_eta_ptY_std40.txt"   # copy of the output of wp8/remote/remote_eta_ptY.sh
OUT = "/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all/wp8_pty.txt"
OUTJ = "/hdd3/przygoda/HADES/FAT_out/all_report/wp8/wp8_pty.json"
NB2RAW = 6079.2 * 63.0 / 64.0
YF = 64.0 / 63.0
RP, FLEAK = 1.1859, 0.1544            # WP6 results.txt (side-band R' and leakage f)
CBP, CBSB = 3279.0, 2695.0            # WP6 window CB (spread rule)
EGF = 1.19                            # event-group variance factor (WP5: var(F10 total) = 1.19 N in S0)
NPT, NY = 75, 100                     # wp8_sim.C grid: pT 20 MeV/c from 0 (last bin = overflow), y 0.05 from -2.5 (clamped)
rng = np.random.default_rng(8)


def load():
    sim = {k: np.zeros((NPT, NY)) for k in ("C4_P", "C5_P", "C2_B", "C3_B", "Geta_P", "C45_SB")}
    simv = {k: np.zeros((NPT, NY)) for k in sim}
    dat = {(w, s): np.zeros((NPT, NY)) for w in ("P", "SB") for s in ("us", "pp", "mm")}
    with open(GRID) as f:
        for line in f:
            if line.startswith("#"):
                continue
            a = line.split()
            if a[0] == "sim":
                sim[a[1]][int(a[2]), int(a[3])] += float(a[4]); simv[a[1]][int(a[2]), int(a[3])] += float(a[5])
            else:
                dat[(a[1], a[2])][int(a[3]), int(a[4])] += float(a[5])
    gen = np.zeros((100, 100))
    hdr = ""
    with open(GEN) as f:
        for line in f:
            if line.startswith("#"):
                hdr = line.strip(); continue
            a = line.split()
            if len(a) == 3:
                gen[int(a[0]), int(a[1])] += float(a[2])
    return sim, simv, dat, gen, hdr


def tertiles(w1d, step):
    c = np.cumsum(w1d) / w1d.sum()
    e1 = (np.searchsorted(c, 1 / 3) + 1) * step
    e2 = (np.searchsorted(c, 2 / 3) + 1) * step
    return e1, e2


def cellsum(h, ptedges_bins, yedges_bins):
    """sum of a (pT, y) grid in 3 x 3 cells; edges given as bin indices [0, i1, i2, n]"""
    out = np.zeros((3, 3))
    for a in range(3):
        for b in range(3):
            out[a, b] = h[ptedges_bins[a]:ptedges_bins[a + 1], yedges_bins[b]:yedges_bins[b + 1]].sum()
    return out


def main():
    sim, simv, dat, gen, hdr = load()
    S45 = sim["C4_P"] + sim["C5_P"]
    V45 = simv["C4_P"] + simv["C5_P"]
    # 1. cells from the tertiles of the simulated peak signal
    p1, p2 = tertiles(S45.sum(axis=1), 1)     # in bins of 20 MeV/c (edge index)
    y1, y2 = tertiles(S45.sum(axis=0), 1)
    p1, p2, y1, y2 = int(p1), int(p2), int(y1), int(y2)
    pe = [0, p1, p2, NPT]; ye = [0, y1, y2, NY]
    pe_g = [0, p1, p2, 100]; ye_g = [0, y1, y2, 100]     # generator histogram: 100 pT bins (0-2 GeV/c), same y grid
    pt_edges = [0, 20 * p1, 20 * p2, float("inf")]
    y_edges = [-float("inf"), -2.5 + 0.05 * y1, -2.5 + 0.05 * y2, float("inf")]
    # 2. data per cell
    NP = cellsum(dat[("P", "us")], pe, ye); LSP = cellsum(dat[("P", "pp")] + dat[("P", "mm")], pe, ye)
    NS = cellsum(dat[("SB", "us")], pe, ye); LSS = cellsum(dat[("SB", "pp")] + dat[("SB", "mm")], pe, ye)
    cbP = CBP * LSP / LSP.sum(); cbS = CBSB * LSS / LSS.sum()
    den = 1.0 - RP * FLEAK
    Y = YF * ((NP - cbP) - RP * (NS - cbS)) / den
    # variance: event groups on N, Poisson-like CB share
    VY = (YF / den) ** 2 * (EGF * (NP + RP ** 2 * NS) + (cbP ** 2 / np.maximum(LSP, 1) + RP ** 2 * cbS ** 2 / np.maximum(LSS, 1)))
    # 3. simulation
    S = cellsum(S45, pe, ye) * NB2RAW * YF
    VS = cellsum(V45, pe, ye) * (NB2RAW * YF) ** 2
    SB2 = cellsum(sim["C2_B"] + sim["C3_B"], pe, ye) * NB2RAW * YF
    G = cellsum(gen, pe_g, ye_g)
    g = G / G.sum()
    r = Y / S
    aP_cells = Y.sum() / S.sum()
    aPp = float(np.sum(g * r))
    kP = aP_cells / aPp
    w = r / aPp
    kB = float(np.sum(w * SB2) / SB2.sum())
    s_acc = S / S.sum()
    # 6. toys
    nt = 2000
    dA, dKB = [], []
    for _ in range(nt):
        Yt = Y + rng.normal(size=Y.shape) * np.sqrt(VY)
        St = S + rng.normal(size=S.shape) * np.sqrt(VS)
        Gt = rng.multinomial(int(G.sum()), (G / G.sum()).ravel()).reshape(G.shape)
        gt = Gt / Gt.sum()
        rt = Yt / St
        a0 = Yt.sum() / St.sum(); a1 = float(np.sum(gt * rt))
        dA.append(a1 / a0 - 1)
        dKB.append(float(np.sum(rt / a1 * SB2) / SB2.sum()) / (a0 / a1) - 1)
    dA = np.array(dA); dKB = np.array(dKB)
    rel = aPp / aP_cells - 1
    sig = float(np.std(dA))
    counted = math.sqrt(max(0.0, rel * rel - sig * sig))
    rel_DR = kB / kP - 1
    sig_DR = float(np.std(dKB))
    lines = []
    lines_extra = []
    def P(s=""):
        print(s); lines.append(s)
    P("# wp8_pty.txt -- plots/all_report/wp8_pty.py (WP8, 2026-11-14): (pT4, y4) reweighting variant (3 x 3 cells) of the eta production")
    P("# data S0, isBest >= 0 (wp8_sim.C grid); simulation C4 + C5 in P, S0, V3, SMASH x PDG BR; generator: %s" % hdr)
    P("# cells: pT4 edges [MeV/c] %s ; y_cm edges %s (inner edges = tertiles of the simulated C4 + C5 in P)" %
      (", ".join("%.0f" % x if math.isfinite(x) else "inf" for x in pt_edges), ", ".join("%.2f" % x if math.isfinite(x) else "%sinf" % ("-" if x < 0 else "") for x in y_edges)))
    P("")
    P("%-24s %9s %9s %9s %9s %9s %9s %9s" % ("cell (pT, y)", "N_P", "N_SB", "Y_c", "+-", "S_c", "r_c", "g_c"))
    for a in range(3):
        for b in range(3):
            P("%-24s %9.0f %9.0f %9.1f %9.1f %9.1f %9.4f %9.4f" % ("pT%d y%d" % (a + 1, b + 1), NP[a, b], NS[a, b], Y[a, b], math.sqrt(VY[a, b]), S[a, b], r[a, b], g[a, b]))
    P("%-24s %9.0f %9.0f %9.1f %9.1f %9.1f" % ("sum", NP.sum(), NS.sum(), Y.sum(), math.sqrt(VY.sum()), S.sum()))
    chi2 = float(np.sum((Y - aP_cells * S) ** 2 / (VY + aP_cells ** 2 * VS)))
    P("")
    P("sum of the cell yields %.1f (WP6 side-band Y_peak 5088.0); accepted fraction of the simulated signal per cell = %s" %
      (Y.sum(), " ".join("%.3f" % x for x in s_acc.ravel())))
    P("flatness of r_c: chi2 = %.1f / 8 (p = %.3f) against a common a_P" % (chi2, 1 - __import__("scipy.stats", fromlist=["chi2"]).chi2.cdf(chi2, 8)))
    P("a_P (accepted-weighted, cells) = %.4f ; a_P' (produced-weighted, reweighted SMASH production) = %.4f ; ratio %.4f" % (aP_cells, aPp, aPp / aP_cells))
    P("  => sigma(pp -> eta X) and a_P change by %+.2f %% (toy stat noise of the variant +- %.2f %%; counted after the Barlow subtraction %.2f %%)" %
      (100 * rel, 100 * sig, 100 * counted))
    P("  peak efficiency factor k_P = %.4f ; bump (C2 + C3 in B, reconstructed (pT4, y4) as proxy) k_B = %.4f" % (kP, kB))
    P("  => DR x k_B/k_P = %+.2f %% (noise +- %.2f %%), a_B x 1/k_B = %+.2f %%, BR_G-type ratios (C5/C2) x k_B/k_P = %+.2f %%" %
      (100 * rel_DR, 100 * sig_DR, 100 * (1 / kB - 1), 100 * rel_DR))
    # information: one-dimensional variants (3 pT cells x all y; all pT x 3 y cells)
    for name, pe1, ye1, pg1, yg1 in (("pT only (3 x 1)", pe, [0, NY], pe_g, [0, 100]), ("y only (1 x 3)", [0, NPT], ye, [0, 100], ye_g)):
        def cs(h, a, b):
            out = np.zeros((len(a) - 1, len(b) - 1))
            for i in range(len(a) - 1):
                for j in range(len(b) - 1):
                    out[i, j] = h[a[i]:a[i + 1], b[j]:b[j + 1]].sum()
            return out
        NP1 = cs(dat[("P", "us")], pe1, ye1); LP1 = cs(dat[("P", "pp")] + dat[("P", "mm")], pe1, ye1)
        NS1 = cs(dat[("SB", "us")], pe1, ye1); LS1 = cs(dat[("SB", "pp")] + dat[("SB", "mm")], pe1, ye1)
        Y1 = YF * ((NP1 - CBP * LP1 / LP1.sum()) - RP * (NS1 - CBSB * LS1 / LS1.sum())) / den
        VY1 = (YF / den) ** 2 * EGF * (NP1 + RP ** 2 * NS1)
        S1 = cs(S45, pe1, ye1) * NB2RAW * YF; VS1 = cs(V45, pe1, ye1) * (NB2RAW * YF) ** 2
        G1 = cs(gen, pg1, yg1); g1 = G1 / G1.sum()
        a0 = Y1.sum() / S1.sum(); a1 = float(np.sum(g1 * Y1 / S1))
        d = []
        for _ in range(nt):
            Yt = Y1 + rng.normal(size=Y1.shape) * np.sqrt(VY1); St = S1 + rng.normal(size=S1.shape) * np.sqrt(VS1)
            d.append(float(np.sum(g1 * Yt / St)) / (Yt.sum() / St.sum()) - 1)
        SB1 = cs(sim["C2_B"] + sim["C3_B"], pe1, ye1) * NB2RAW * YF
        kP1 = a0 / a1; kB1 = float(np.sum((Y1 / S1) / a1 * SB1) / SB1.sum())
        r1 = a1 / a0 - 1; n1 = float(np.std(d)); c1 = math.sqrt(max(0.0, r1 * r1 - n1 * n1))
        P("INFO %-16s r = %s, g = %s: a_P'/a_P - 1 = %+.2f %% +- %.2f %% (noise), counted %.2f %%; k_B/k_P - 1 (DR) = %+.2f %%, 1/k_B - 1 (a_B) = %+.2f %%" %
          (name, " ".join("%.3f" % x for x in (Y1 / S1).ravel()), " ".join("%.3f" % x for x in g1.ravel()), 100 * r1, 100 * n1, 100 * c1,
           100 * (kB1 / kP1 - 1), 100 * (1 / kB1 - 1)))
        lines_extra.append("pty1d %s rel %.6f noise %.6f counted %.6f relDR %.6f aB_rel %.6f" % (name.split()[0], r1, n1, c1, kB1 / kP1 - 1, 1 / kB1 - 1))
    with open(OUT, "w") as f:
        f.write("\n".join(lines) + "\n")
        f.write("pty aP_cells %.6f aP_rw %.6f rel %.6f noise %.6f counted %.6f kP %.6f kB %.6f relDR %.6f noiseDR %.6f chi2 %.3f\n" %
                (aP_cells, aPp, rel, sig, counted, kP, kB, rel_DR, sig_DR, chi2))
        for l in lines_extra:
            f.write(l + "\n")
    json.dump(dict(pt_edges=[float(x) if math.isfinite(x) else None for x in pt_edges], y_edges=[float(x) if math.isfinite(x) else None for x in y_edges],
                   Y=Y.tolist(), VY=VY.tolist(), S=S.tolist(), VS=VS.tolist(), g=g.tolist(), r=r.tolist(), aP_cells=aP_cells, aP_rw=aPp,
                   rel=rel, noise=sig, counted=counted, kP=kP, kB=kB, relDR=rel_DR, noiseDR=sig_DR, chi2=chi2, aB_rel=1 / kB - 1), open(OUTJ, "w"), indent=1)


if __name__ == "__main__":
    main()
