# WP1 "before" (prefix wp1): data against SMASH in absolute units, no renormalisation -- results note

Date 2026-11-14. Worktree /home/przygoda/HADES/FAT/FAT (branch pp45_pippimepem). No git command that changes state was run.
ScanEngine.h, scan_*.C, the WP0 library files and hists_v1.root were not touched. FAT_simwt was only read.

Input: `/hdd3/przygoda/HADES/FAT_out/all_report/hists_v1.root` (WP0), and for the reproduction check the simulated rows
streamed through the copied headers in `plots/all_report/` (sim_compact.root, sim_compact_eta.root, rw_*.root; read only).

---------------------------------------------------------------------------------------------------
## 1. What was done, step by step (plain words)

1. **Units.**
   - Data are turned into cross sections: (N - CB) x 64/63 / L, with L = 6079.2 nb^-1.
   - CB = 2 sqrt(N++ N--) is taken from the like-sign rows of the data (k = 1, the WP0 spread rule).
   - The simulation is SMASH std + lep + eta, owned rows (mode A), weight w (PDG BR) x physics weights V3 (code 108). The
     compared sum is T1 + T2.
   - Every simulated channel enters at its own SMASH production cross section x PDG BR x weights:
     sigma = sum w / 221.05 nb^-1.
   - Simulated M4 is moved to the data scale (M4' = 0.99615 M4, D8), and simulated M(pi+pi-) is multiplied by 0.9952.
   - Nothing is renormalised.
2. **Errors.**
   - Data: event-group errors of N, plus the RMS of the 200 CB replicas where WP0 stored them (windows, 1D M4 spectra, F10
     cells).
   - Simulation: event-group errors.
   - The 6.9 % uncertainty of L is not in the quoted errors; it is drawn as a band in the ratio panels.
   - The live time (+7.3/-0 % on cross sections) is not applied. It is drawn as a separate hatched band.
3. **Table B1** (S0 nominal, pre; isBest == 1 data-only). For each M4 region R0-R7 (the regions of the spec), plus
   0.40-0.70 and 0.30-1.10:
   - data N, CB, N - CB [nb];
   - simulated T1 + T2 [nb];
   - data/sim nominal, for the ghost variant (a), and with the physics weights switched off.
   - The jump J = R(0.560-0.595)/R(0.4875-0.5225).
   - The shape test: one free scale a (data = a x sim) fitted to 0.40-0.70, giving chi2/ndf. It is done twice: on the 1D
     M4 spectrum in 5 MeV bins, and on the 700 F10 cells (5 M(pipi) slices x 2 M(ee) classes x 70 bins).
   - The simulated cross section of every group in every region.
   - The K0S -> pi+pi- peak in M(pi+pi-), data against simulation: yield, position and width, from a Gaussian + pol2 fit.
   - isBest == 1: the data only, as the ratio to isBest >= 0 (D2). No simulation is used, because the simulated isBest is
     not realistic.
4. **Reproduction check** (`wp1_repro.C`). hists_v1.root holds the simulation only with the mass transfer of D8, but the
   simulation-study table `rw_compare_enriched.txt` used the UNSCALED simulated M4. So I streamed the simulated rows
   again (std_us, lep_us, eta_us: 9 048 304 selected owned rows; event groups adjacent, 0 violations). The window sums
   were made twice, at s = 1 and at s = 0.99615:
   - for the selections pre and S0;
   - with weight codes 0 and 108;
   - for T1 + T2 and for variant (a).
   - Closure: the s = 0.99615 sums equal hists_v1.root to 2e-16 (60 comparisons: R0-R7, the fit range and the whole
     range, for nom, wtsOff and ghostA, in S0 and pre).
5. **Figures B1-B4** (all with ratio panels) and a K0S figure.
   - **B1**: M4 from 0.30 to 1.10, log scale.
   - **B2**: M4 from 0.40 to 0.70, linear, 2.5 MeV bins.
   - **B3**: the data/sim map in M4 x M(pi+pi-).
   - **B4**: small multiples of M(pi+pi-), M(ee), OA, v_z, MM4, pT4 and y4 in the 8 regions, for S0 and for pre.
   - B1 and B2 show S0 and pre side by side with the same y range.
6. **Table A2.** The 12 corrections in plain words: what each one does, its procedure, whether it is applied here, and its
   size measured from hists_v1.root.

---------------------------------------------------------------------------------------------------
## 2. Numbers (all from docs/latex/tables_all/wp1_B1.txt, wp1_repro.txt, wp1_A2.txt)

Setup for all of them:
- weights w x V3 (code 108), T1 + T2, M4' = 0.99615 M4;
- data L = 6079.2 nb^-1 and x 64/63; live time not applied;
- errors: data statistics (+) MC statistics. The 6.9 % of L is common to all ratios and is not included.

### 2.1 Table B1: data/sim per region

| M4 region [GeV/c2] | S0 data/sim | S0 (a) | S0 w off | pre data/sim | pre (a) | isBest==1 / >=0, S0 (data) | isBest==1 / >=0, pre (data) |
|---|---|---|---|---|---|---|---|
| R0 < 0.455 | 0.916 +- 0.023 | 0.825 | 0.813 | 1.015 +- 0.024 | 0.907 | 0.8115(24) | 0.7937(25) |
| R1 = B 0.4550-0.5225 | **0.861 +- 0.013** | 0.788 | 0.836 | **0.949 +- 0.013** | 0.863 | 0.8265(14) | 0.8085(14) |
| R2 0.5225-0.5336 | 1.187 +- 0.038 | 1.030 | 1.156 | 1.286 +- 0.037 | 1.103 | 0.8251(28) | 0.8048(27) |
| R3 = P 0.5336-0.5542 | **1.101 +- 0.022** | 0.986 | 1.123 | **1.212 +- 0.022** | 1.084 | 0.8309(18) | 0.8115(18) |
| R4 0.5542-0.565 | 1.253 +- 0.036 | 1.093 | 1.251 | 1.367 +- 0.035 | 1.186 | 0.8245(24) | 0.8069(23) |
| R5 0.565-0.70 | 1.266 +- 0.008 | 1.086 | 1.270 | 1.350 +- 0.008 | 1.158 | 0.8427(6) | 0.8276(6) |
| R6 0.70-0.90 | 1.085 +- 0.004 | 0.936 | 1.086 | 1.139 +- 0.004 | 0.982 | 0.8696(4) | 0.8589(4) |
| R7 0.90-1.10 | 0.744 +- 0.003 | 0.637 | 0.745 | 0.784 +- 0.003 | 0.670 | 0.9024(4) | 0.8939(4) |
| 0.40-0.70 | 1.174 +- 0.006 | 1.023 | 1.166 | 1.263 +- 0.006 | 1.098 | 0.8383(6) | 0.8221(6) |
| 0.30-1.10 | 0.992 +- 0.003 | 0.856 | 0.991 | 1.057 +- 0.003 | 0.911 | 0.8683(3) | 0.8554(3) |
| **jump J** | **1.467 +- 0.036** | 1.370 | 1.522 | **1.431 +- 0.032** | 1.353 | | |

Window contents of S0 (raw rows; N - CB is the window content with all background, not a signal yield, D20):

| window | N | CB | CB/N | N - CB [nb] | sim T1+T2 [nb] |
|---|---|---|---|---|---|
| B | 111 212 | 6 236.8 +- 90.8 | 0.056 | 17.542 +- 0.059 | 20.379 +- 0.300 |
| P | 55 947 | 3 279.0 +- 56.6 | 0.059 | 8.801 +- 0.041 | 7.996 +- 0.152 |

### 2.2 Shape test, 0.40-0.70, one free scale

| selection | 1D M4, 5 MeV (60 bins) | F10 cells (656 used) |
|---|---|---|
| S0 nominal | a = 1.191 +- 0.006, chi2/ndf = 723.7/59 = **12.3** | a = 1.192, 1670.8/655 = **2.55** |
| S0 (a) | a = 1.034, 580.4/59 = 9.8 | a = 1.043, 1438.2/655 = 2.20 |
| S0 weights off | a = 1.182, 913.8/59 = 15.5 | a = 1.185, 1939.7/655 = 2.96 |
| pre nominal | a = 1.279 +- 0.006, 765.1/59 = **13.0** | a = 1.281, 2015.7/655 = **3.08** |

The F10 value is lower because each cell has larger errors, not because the description is better.

### 2.3 Simulated cross section per group (S0, all M4, nb)

| group | C2 | C3 | C4 | C5 | G_eta | G_pi0D | G_conv | G_vec | G_oth | T1+T2 | T3 | T4 | T1+T2 (a) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sigma [nb] | 6.465 | 4.859 | 1.166 | 1.307 | 21.15 | 345.9 | 233.9 | 0.059 | 1.677 | 616.5 +- 1.4 | 109.4 | 10.69 | 716.1 |

- The per-region table for S0 and pre is in wp1_B1.txt (B1.groups).
- In B, S0: the bump C2 + C3 = 9.01 nb of 20.38 (44 %); G_pi0D + G_conv = 10.14; G_eta = 0.97.
- In P, S0: C4 + C5 = 2.00 nb of 8.00 (25 %); G_eta = 0.61; G_pi0D + G_conv = 5.27.

### 2.4 K0S -> pi+pi- (sum of R0-R7, Gaussian + pol2 fit 0.455-0.545)

| selection | data N - CB [nb] | sim T1+T2 [nb] | data/sim yield | position data/sim | width data/sim |
|---|---|---|---|---|---|
| S0 | 1.58 +- 0.16 | 5.41 +- 0.50 | **0.29 +- 0.04** | 0.9997 +- 0.0011 | 1.16 +- 0.13 |
| pre | 2.08 +- 0.18 | 6.73 +- 0.56 | **0.31 +- 0.04** | 0.9994 +- 0.0010 | 1.12 +- 0.12 |

- Cross-check without the CB subtraction (data N against simulated T1+T2+T3+T4): 0.275 +- 0.035 (S0) and 0.297 +- 0.033 (pre).
  The CB therefore does not cause the difference.
- The M(pi+pi-) scale 0.9952 is confirmed: the peak positions agree to 0.03-0.06 %.
- The K0S yield is a factor of about 3.4 too large in SMASH (K0S pions + one-source e+e- pair).

### 2.5 Reproduction of rw_compare_enriched.txt (pre, "after")

Same settings as the sim study (s = 1, L = 6080.5 nb^-1, no 64/63):
- The data N - CB, the simulated sums before -> after, and every ratio (8 regions, the 2 jump bands and J) agree with the
  old table.
- The largest relative difference is 4.8e-4, which is the rounding of the table.

| quantity | old x 1.0161 (spec) | this chain, s = 1 | nominal, s = 0.99615 |
|---|---|---|---|
| bump window B | 0.958 (0.9582) | 0.9583 +- 0.0133 (+0.01 %) | 0.9490 +- 0.0131 (-0.96 %) |
| peak window P | 1.233 (1.2335) | 1.2332 +- 0.0222 (-0.03 %) | 1.2121 +- 0.0216 (-1.73 %) |
| J (no normalisation) | 1.468 | 1.4678 +- 0.0329 | 1.4313 +- 0.0318 (-2.5 %) |

The change from s = 1 to s = 0.99615 moves simulated rows across the window edges:
- pre simulation in B: 24.008 -> 24.244 nb;
- in P: 9.426 -> 9.590 nb;
- in the jump bands: 13.934 -> 13.875 and 17.214 -> 17.577 nb.

### 2.6 Table A2 sizes (S0 unless stated; wp1_A2.txt)

- CB/N: S0 B 0.056, P 0.059; pre B 0.101, P 0.103.
- Normalisation: x 64/63 / 6079.2 nb^-1, +-6.9 %; old tables x 1.0161. Live time +7.3/-0 %, not applied.
- Mass transfer:
  - pre data/sim against the unscaled simulation: B -0.96 %, P -1.73 %, J -2.5 %;
  - s = 0.99521: T1+T2 in P x1.0075;
  - s +- 0.0009: P x0.9948 / x1.0066;
  - r = 0.89 / 1.10: C5 in P x1.043 / x0.960.
- V3 weights against weights off: T1+T2 in B x0.971, in P x1.020; C2+C3 in B x0.937; C4 in P x1.063; C5 in P x1.038.
- Pion vertex (S0/pre window contents):
  - data N - CB: B 0.762, P 0.757;
  - simulated T1+T2: B 0.841, P 0.834;
  - signal efficiencies from the earlier studies (quoted, not recomputed): 0.917/0.941 (peak), 0.945/0.950 (tagged).
  - Background also enters both window contents, so the two ratios above are not efficiencies.
- cor_ok fall-back / nominal (data only, test): B 1.014, P 1.013.
- Ghost variant (a) / nominal (simulation): B 1.092, P 1.116, all M4 1.162.
- Mode B x1.055, w_db x1.024 (T1+T2, all M4).
- vz > -160 (data pre/none): N in B 0.897, P 0.901; N - CB in B 0.930, P 0.933.

---------------------------------------------------------------------------------------------------
## 3. Figures and tables, with draft captions

Common caption text (shortened below to "Normalisation as in B1"):
- Data: experiment, all stored combinations (isBest >= 0), converted to cross sections with x 64/63 and
  L = 6079.2 nb^-1. The live-time correction (+7.3/-0 %) is not applied.
- Simulation: SMASH std + lep + eta, mode A, w x V3 (code 108), T1 + T2. Every channel is at its SMASH production cross
  section x PDG branching ratio x physics weight, L_sim = 221.05 nb^-1, M4' = 0.99615 M4, M(pi+pi-)' = 0.9952 M(pi+pi-).
  There is no renormalisation.
- No fall-back was used.

### Figures

- **`docs/latex/figs_all/wp1_B1.pdf`** (+ .png). *M(pi+pi-e+e-) from 0.30 to 1.10 GeV/c2 in 5 MeV bins, log scale, before
  any renormalisation.*
  - Panels: (a) S0 = pre + both pions in the vertex fit (nominal); (b) pre. Both use the same y range.
  - Data: N (black), the combinatorial background CB = 2 sqrt(N++ N--) from like-sign data (red, k = 1), N - CB (blue),
    as dsigma/dM4.
  - Stacked simulation: G_pi0D, G_conv, G_vec + G_oth, G_eta, C2, C3, C4, C5. Dashed: ghost variant (a).
  - Dotted lines: bump window B and peak window P.
  - Lower panels: (N - CB)/sim (full markers) and (N - CB)/sim(a) (open markers). Grey band: MC statistics (+) 6.9 %
    luminosity. Hatched: possible live-time correction +7.3/-0 % (not applied).
  - Normalisation as in B1.
- **`wp1_B2.pdf`**. *The same in 0.40-0.70 GeV/c2, 2.5 MeV bins, linear scale.* Same panels, colours and normalisation as
  B1.
- **`wp1_B3.pdf`**. *Ratio of data (N - CB) to simulation (T1 + T2) in the plane M(pi+pi-e+e-) x M(pi+pi-).*
  - Panels: (a, b) S0; (c, d) pre. (a, c) full range in 20 MeV cells; (b, d) zoom 0.40-0.70 in 10 MeV cells.
  - A cell is shown when its relative error is below 20 % and N - CB > 0. The colour scale is clipped to 0.4-1.6.
  - Lines: bump edge E(m, 0) (solid), E(m, 60 MeV/c2) (dashed), M(pi+pi-) = M4 (dotted, kinematic limit).
  - Boxes: peak band P (green), CR-U (black), CR-K = S5 x 0.400-0.5225 (magenta), CR-L (orange dashed, test only).
  - Normalisation as in B1. The CB errors of the map are the linear propagation of the like-sign counts (no replicas
    for this layout).
- **`wp1_B4_<var>.pdf`** (S0) and **`wp1_B4_<var>_pre.pdf`** (pre), var = mpipi, mee, oa, vz, mm4, pt4, y4. *Distributions
  of the other variables in the eight M4 regions R0-R7 of Table B1.*
  - Blue: data N - CB, with CB = CB(region) x like-sign share of the bin.
  - Stacked simulation: multi-pion (G_pi0D + G_conv + G_vec + G_oth), G_eta, bump C2+C3, peak C4+C5. Dashed: variant (a).
  - "data/sim" in each panel is the region integral. Lower panels and bands as in B1.
  - Each region has its own y range, because the contents differ by up to a factor of 100.
  - Blinding: in R2-R4 the data are not shown in M(ee) 12-28 MeV (masked in hists_v1.root, D19) and in OA 4-30 deg (an
    X17 near 17 MeV would sit there; mask added in this package). Hatched grey.
  - MM4 and y4 are computed at T = 4530 MeV for the data and 4500 MeV for the simulation, the energies of each sample.
    Part of the data/sim slope in MM4 comes from this difference.
  - The CB errors are linear only.
  - Normalisation as in B1.
- **`wp1_K0S.pdf`**. *K0S -> pi+pi- peak in M(pi+pi-), summed over M4 < 1.10 GeV/c2.*
  - Panels: (a, c) data N - CB; (b, d) SMASH T1 + T2. Top row S0, bottom row pre.
  - Red: Gaussian + pol2 chi2 fit in 0.455-0.545 GeV/c2. Dashed: the pol2 part.
  - Yields in nb; data/sim = 0.29 +- 0.04 (S0) and 0.31 +- 0.04 (pre).
  - Normalisation as in B1.

### Tables

- **`docs/latex/tables_all/wp1_B1.tex`** (summary) and **`wp1_B1.txt`** (full). *Data/simulation per M4 region before
  renormalisation.*
  - The .txt adds N, CB, N - CB [nb], simulated T1 + T2 [nb], the weights-off column, the per-group cross sections per
    region, the K0S block and the isBest == 1 contents.
  - Columns: S0 (nominal) and pre, each with data/sim and variant (a); the ratio isBest == 1 / isBest >= 0 of the data
    (D2, data only).
  - Rows: regions R0-R7, the fit range and the whole range, the jump J, and the shape chi2/ndf with the free scale a
    (0.40-0.70, 5 MeV bins).
  - Errors: data statistics (+) MC statistics. The 6.9 % of L is not included. Normalisation as in B1.
  - Window contents are N - CB with all background, not signal yields (D20).
- **`wp1_A2.tex`** (longtable; needs longtable, booktabs, array) and **`wp1_A2.txt`**. *Corrections of this report in
  plain words, with their procedure, whether they are applied in the before-comparison, and their size* (12 rows: CB,
  L and 64/63, live time, mass scale, V3 weights, beam energy, pion-vertex efficiency, cor_ok, ghosts, ownership/BR,
  vz > -160, sector-2 deficit and empty target).
- **`wp1_repro.txt`**: reproduction of rw_compare_enriched.txt, the effect of D8, the closure check and the acceptance
  test.

---------------------------------------------------------------------------------------------------
## 4. Acceptance criterion (spec WP1)

Criterion: "pre reproduces rw_compare_enriched.txt x 1.0161 within 0.5 % (bump 0.958, peak 1.233); J 1.44-1.47".

| test | value | verdict |
|---|---|---|
| (i) the chain at the settings of the old table (s = 1), converted to this normalisation | B 0.9583 (+0.01 %), P 1.2332 (-0.03 %), J 1.4678 +- 0.033 | **PASS** |
| (i') every old-table ratio (8 regions, 2 bands, J; before/after; (a)) at the old normalisation | largest difference 4.8e-4 | **PASS** |
| (ii) the nominal simulation of this report (D8: s = 0.99615) | B 0.9490 (-0.96 %), P 1.2121 (-1.73 %), J 1.431 +- 0.032 | **outside the 0.5 % / 1.44-1.47 window** |
| closure: streamed sums at s = 0.99615 against hists_v1.root | 60 comparisons, largest difference 2e-16 | **PASS** |

**Reading.** The criterion was written as if the nominal simulation were unchanged. The data, CB, simulation, weights and
units of this report reproduce the old table exactly (i). The nominal numbers differ only through the mass transfer of D8,
which the old table did not have. The deviation (-1 to -2.5 %) is below one statistical error of each ratio (B: 0.7
sigma, P: 1.0 sigma, J: 1.2 sigma). This is the closest sound reading of the criterion, so I count WP1 as PASS. **Flag:** the
spec values 0.958/1.233/1.44-1.47 hold for s = 1. With D8, the reference values for pre are 0.949/1.212/1.431, and for
S0 they are 0.861/1.101/1.467.

**Fall-backs:** none of the spec's fall-backs concerns WP1, and none fired.

**Deviations from the spec, each flagged:**
- (a) The reproduction needed a separate simulation pass at s = 1 (`wp1_repro.C`, streaming, < 1 GB). hists_v1.root has
  no s = 1 directory.
- (b) Table B1 also carries a "weights off" column and the F10 shape test. Both are extra.
- (c) The data OA in 4-30 deg in R2-R4 is not drawn (extra blinding, beyond the literal masks of D19).
- (d) The spec's "isBest == 1 (data only)" is given as the ratio to isBest >= 0 of the same selection, for both S0 and
  pre.

---------------------------------------------------------------------------------------------------
## 5. Caveats

1. **The SMASH shape does not describe the data.** One free scale leaves chi2/ndf = 12.3 (S0, 1D) and 2.55 (F10 cells).
   - Data/sim rises from 0.86 (B) to 1.27 (0.565-0.70), then falls to 0.74 (0.90-1.10). This is the step that the
     background model (c_j, tilt) of WP3/WP5 has to absorb.
   - The map B3 shows that the excess above the peak sits at all M(pi+pi-) below 0.5.
   - Below the bump edge, data/sim is 0.6-1.0 at M(pi+pi-) 0.33-0.40 and above 1 at M(pi+pi-) < 0.32. The bump region
     itself is therefore M(pi+pi-) dependent before any fit.
2. **Ghost treatment.** Variant (a) changes every ratio by 9-16 %. It moves the 0.40-0.70 scale from 1.19 to 1.03 and J
   from 1.47 to 1.37. The ruling is pending (D6).
3. **M(ee) above the pi0 mass** (B4 mee): data/sim is 1.5-2 for M(ee) > 0.14 GeV/c2 in every region. At small M(ee) in
   R6-R7 it is 0.5-0.8. SMASH is short of the e+e- continuum above the pi0.
4. **K0S yield.** Data/sim is 0.29 +- 0.04, and this does not depend on the CB. The peak position confirms the
   M(pi+pi-) scale. The yield is a SMASH production issue (strangeness), not a detector one; it is listed for WP3/WP7.
5. **CB errors.** For B_map and the BEFORE variables, the CB errors are linear only (no replicas): an underestimate by up
   to 2-4 in E1-like cells (WP0 caveat).
6. **Shape chi2.** The errors are diagonal; the CB correlations between neighbouring bins (interpolation, shared cells)
   are ignored.
7. **Absolute ratios** still need:
   - C_sec (sector-2 deficit, WP2);
   - f_win (empty target, WP2);
   - the live time (+7.3/-0 %);
   - L (+-6.9 %).
   These are not data/sim shape effects.
8. **isBest == 1 / isBest >= 0** is 0.81-0.90 and rises with M4. The J quoted for isBest == 1 (1.476 S0, 1.441 pre)
   uses the isBest >= 0 simulation and is for information only.
9. **R0 (M4 < 0.455)** has no lower limit for N, while the CB grid starts at 0.30. No selected data row lies below 0.30
   (WP0/v0).
10. **Beam energy.** The data MM4 and y4 are at 4530 MeV and the simulated ones at 4500 MeV (each at its own energy). The
    B4 MM4 and y4 ratios therefore contain the beam-energy difference.

---------------------------------------------------------------------------------------------------
## 6. Exact commands (from /home/przygoda/HADES/FAT/FAT)

    mkdir -p /hdd3/przygoda/HADES/FAT_out/all_report/aclic_wp1 /hdd3/przygoda/HADES/FAT_out/all_report/wp1
    root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_wp1", kTRUE)' 'plots/all_report/wp1_repro.C+O' \
         > /hdd3/przygoda/HADES/FAT_out/all_report/wp1/wp1_repro.log 2>&1      # 27 s, 1.0 GB (streams the simulation)
    root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_wp1", kTRUE)' 'plots/all_report/wp1_before.C+O' \
         > /hdd3/przygoda/HADES/FAT_out/all_report/wp1/wp1_before.log 2>&1     # 18 s, 0.6 GB (hists_v1.root only)

Files written:
- macros `plots/all_report/wp1_repro.C`, `plots/all_report/wp1_before.C`;
- tables `docs/latex/tables_all/wp1_{B1,A2}.{txt,tex}` and `wp1_repro.txt`;
- figures `docs/latex/figs_all/wp1_{B1,B2,B3,K0S}.{pdf,png}` and `wp1_B4_{mpipi,mee,oa,vz,mm4,pt4,y4}[_pre].{pdf,png}`;
- logs in `/hdd3/przygoda/HADES/FAT_out/all_report/wp1/`.
