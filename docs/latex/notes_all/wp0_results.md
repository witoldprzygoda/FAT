# WP0 histogram factory (prefix wp0) -- results note

Date 2026-11-14. Branch pp45_pippimepem (worktree /home/przygoda/HADES/FAT/FAT). Nothing committed, no git state changed;
ScanEngine.h and scan_*.C untouched; FAT_simwt only read.

**Main output:** `/hdd3/przygoda/HADES/FAT_out/all_report/hists_v1.root` (141.7 MB, write-protected).
**Sealed X17 cells:** `/hdd3/przygoda/HADES/FAT_out/all_report/sealed/hists_v1_x17_sealed.root` (write-protected; open only after x17_frozen.txt).
**Library:** `/home/przygoda/HADES/FAT/FAT/plots/all_report/` (README.md = naming scheme, origin and every change of the copies).

---------------------------------------------------------------------------------------------------
## 1. What was done, step by step

1. **Shared library.**
   - Copied from FAT_simwt (commit 7206630 of pp45_pippimepem_betanew_sim): `SimEngine.h`, `sim_eta_physics.h`,
     `s5_blind_ana_fit.h`, `s5_blind_ana_tmpl.h`, `src/sim_truth.h`, and in addition `sim_ghostvar.h` (it defines the ghost
     variant (a)). md5 of every source file is in README.md.
   - Changes, only in `SimEngine.h`:
     - `#include "../../src/sim_truth.h"` became `#include "sim_truth.h"` (the local copy).
     - `#include "../ScanEngine.h"` was kept: it now finds the DATA-branch `plots/ScanEngine.h`.
     - `kFigDir` / `kTabDir` now point to `/hdd3/przygoda/HADES/FAT_out/all_report/wp0/{figs,tables}`.
     - A comment says that `kLdata` (6080.5, unchanged) is not the normalisation of this report.
   - The data `scan::Row` extension (`pip_z`, `pim_z`, `dz_pi`) does not touch the simulation reader: `simeng::Row` has its
     own field list. All 10 510 151 owned simulated e+e- rows load, and every simulated total of the sim study is
     reproduced (Section 2.1).
2. **defs.h** (namespace `alldefs`) freezes the spec definitions:
   - L = 6079.2 nb^-1, 64/63, live time +7.3/-0 %, old-table factor 1.0161;
   - conversions: `kNbToRaw` = 5984.2 (nb -> expected raw rows), `kSimToRaw` = 27.0718 (sum w -> raw rows), `kRawToNb`;
   - periods by day; mass transfer s, r and the variants of s and r; M(pipi) x 0.9952;
   - the edge E(m, q), checked against the 5 values of Sect. 2;
   - the selections pre, S0 and their variants (templates usable for data and simulation rows);
   - the layouts F10, F5g, X2D, G3, BEFORE, fine;
   - 21 windows and control regions; the CB parameters;
   - the simulation components (22) and variants (15), physics-weight codes, file names.
   - `hists_io.h` is a small reader.
3. **One factory pass** (`wp0_fill.C`, 157 s, 7.4 GB). Data first, then the simulation. Every histogram is filled from the
   rows directly: integer cell index, flat arrays, double precision.
   - **Event-group variance** (data and simulation): for every event group and cell, the sum of its weights is squared
     and added.
   - **Halves:** half A = even event group, half B = odd.
   - **Simulation replicas:** 100 replicas. Each event group gets a Poisson(1) weight from a counter-based generator
     (reproducible).
   - **CB replicas:** 200 replicas. Each like-sign event group gets a Poisson(1) weight; the 10 MeV cells and the
     fractions f are then rebuilt.
4. **Data.**
   - Inputs: scan_v2 compact (us 4 783 467, pp 122 875, mm 371 524 rows).
   - 11 selections: S0, pre, dz40, PEAK, vz150, vz170, none, S0_best, pre_best, S0_corfb, pre_corfb.
   - 7 subsets: all, fold0, fold1, perI, perIIA, perIIB, perIII.
   - Stored quantities: N, Npp, Nmm (raw rows, event-group errors); CB; CB20; CBrep (subset all).
   - **CB:**
     - 2 sqrt(N++ N--) in 10 MeV cells on the grid from 0.30.
     - Linear interpolation to the 2.5 MeV fine bins, the ScanEngine rule.
     - Spread over the cells of a layout by f(cell | 50 MeV band) from N++ + N--.
     - Windows that are not on the fine grid take fractional fine bins.
     - BEFORE variables: CB(region) x f(bin | region).
5. **Simulation.**
   - Rows: std_us + lep_us + eta_us, owned rows.
   - Weights: w x physics code 108 (V3), mode A. C6 is X17 with weight 2 per row, i.e. per unit BR.
   - Mass transfer:
     - M4' = s [m4_true + r (m4 - m4_true)], with s = 0.99615 and r = 1. Rows without m4_true get s m4.
     - M(pipi) x 0.9952.
     - MM4 and y4 are recomputed at T = 4500 MeV. Validation: at 4530 MeV the same code reproduces the stored values to
       1.2e-7 GeV/c2 and 4.3e-8.
   - Units: nb per bin (sum w / 221.05) with event-group errors.
   - 34 directories:
     - S0 and pre, each with all 15 variants: nom, V2, bump1aY, bumpPS, convOff, GS, wtsOff, modeB, wdb, ghostA, sK0, sUp,
       sDn, r089, r110;
     - dz40, PEAK (MM4(4500) < 2.85), vz150 and vz170, each nominal only.
6. **Blinding (D19).**
   - The data N of X2D in M(ee) 12-28 MeV (all M4 bins) and of B_mee_R2..R4 in 12-28 MeV is set to -1 in hists_v1.root.
     Masks: `data/mask_X2D`, `data/mask_B_mee`.
   - The unmasked cells are written only to the sealed file.
   - N++, N-- and CB are not masked.
7. **Closure (A4)** is checked inside the factory before masking. Every block is tested:
   - F10 sum = fine 0.40-0.70; X2D (with M(ee) flows) = fine 0.5225-0.61; B_M4 = fine; aligned windows = fine;
   - CR-K = F10 S5 bins; G3 = fine 0.40-0.85; F5g = G3 class >= 20 MeV; F10 = G3;
   - B_map = fine; every B_<var>_R<k> (with flows) = window R_k.
   - The check covers N, Npp, Nmm and CB of every data block and all 22 components of every simulation directory.
8. **Read-back check** (`wp0_check.C`):
   - acceptance values re-read from the file;
   - F10 closure for 1056 (directory, quantity) pairs;
   - CB spread rule vs per-cell geometric mean;
   - replica validation;
   - 3 figures and 2 tables.

---------------------------------------------------------------------------------------------------
## 2. Numbers

### 2.1 Reproduction of the simulation study (pre, mode A, w x V3, all M4, s = 1 where windows are sim-scale)

| quantity | this | sim-study table | |
|---|---|---|---|
| T1+T2 = T12 | 701.0832 +- 1.4624 nb | 701.1 +- 1.5 (701.083) | PASS |
| C2+C3 in 0.4576-0.5255 (s = 1) | 9.6887 nb | 9.689 +- 0.23 | PASS |
| C4+C5 in 0.5333-0.5587 (s = 1) | 2.2617 nb | 2.262 | ok |
| C2 / C3 / C4 / C5 | 7.0456 / 5.1898 / 1.2420 / 1.4168 nb | 7.046 / 5.19 / 1.242 / 1.417 | ok |
| G_eta / G_pi0D / G_conv / G_vec / G_oth | 24.273 / 394.50 / 265.35 / 0.10657 / 1.9668 nb | 24.27 / 394.5 / 265.3 / 0.1066 / 1.967 | ok |
| T3 / T4 | 124.96 / 12.510 nb | 125.0 / 12.51 | ok |
| X17 x 1e-5 | 0.071265 nb | 0.07127 | ok |
| C5 in V2 / GS / weights off (all M4) | 1.6442 / 1.4243 / 1.3852 nb | 1.3852 x 1.1870 / 1.0282 / 1 (rw_weights_xsec) | ok |

Omega/rho0 split (pre, all M4, nb):
- G_pi0D_om 37.85, G_pi0D_nom 356.65;
- G_conv_om 23.86, G_conv_nom 241.49;
- G_vec_om 0;
- G_oth_om 0.0033.

The parts sum to the groups to 1e-9.

### 2.2 Simulation per directory

Full table in `docs/latex/tables_all/wp0_sim_totals.txt`. The windows are data-scale windows applied to M4'.

| directory | T12 all M4 [nb] | C2+C3 in B | C4 in P | C5 in P | G_eta in P |
|---|---|---|---|---|---|
| S0__nom | 616.53 +- 1.36 | 9.008 +- 0.225 | 0.955 +- 0.067 | 1.0462 +- 0.0029 | 0.610 +- 0.046 |
| pre__nom | 701.08 +- 1.46 | 9.754 +- 0.233 | 1.007 +- 0.069 | 1.1218 +- 0.0031 | 0.711 +- 0.049 |
| S0__V2 | 616.74 | 9.008 | 0.955 | 1.2147 | 0.610 |
| S0__bumpPS | 617.72 | 9.616 | 0.955 | 1.0462 | 0.576 |
| S0__convOff | 616.46 | 9.008 | 0.898 | 1.0462 | 0.601 |
| S0__modeB | 650.38 | 9.008 | 0.955 | 1.0462 | 0.653 |
| S0__wdb | 631.20 | 8.407 | 1.082 | 1.0462 | 0.617 |
| S0__ghostA | 716.13 | 9.008 | 0.955 | 1.0462 | 0.669 |
| S0__sK0 / sUp / sDn | 616.53 | 8.977 / 9.046 / 8.977 | 0.956 / 0.950 / 0.958 | 1.046 / 1.045 / 1.046 | 0.602 / 0.616 / 0.602 |
| S0__r089 / r110 | 616.53 | 9.067 / 8.952 | 0.982 / 0.929 | 1.092 / 1.005 | 0.613 / 0.607 |
| dz40 / PEAK / vz150 / vz170 | 631.60 / 628.23 / 605.91 / 616.83 | 9.210 / 8.799 / 8.765 / 9.033 | 0.969 / 0.941 / 0.926 / 0.955 | 1.082 / 1.047 / 1.023 / 1.048 | 0.662 / 0.613 / 0.588 / 0.611 |

- **MC statistics** (S0__nom, F10, T12): N_eff per histogram 1 782 (S3E1) to 7 141 (S1E2); C5 N_eff 8 938-60 689.
- **Replicas:** replica RMS / event-group error = 0.989-1.016 (median over bins, per histogram).
- **Halves:** half A / half B = 0.954-1.065 per histogram, consistent with N_eff.

### 2.3 Data (raw rows; CB = spread rule; CB errors = RMS of 200 replicas)

| selection | N in B | CB in B | N in P | CB in P (CB/N) | N in X17P | CB in X17P |
|---|---|---|---|---|---|---|
| S0 | 111 212 | 6 237 +- 91 | **55 947** | 3 279 +- 57 (**0.0586**) | 5 258 | 171.0 +- 6.2 |
| pre | 153 183 | 15 507 +- 164 | 77 513 | 7 949 +- 93 (0.103) | 7 550 | 459 +- 11 |
| dz40 | 124 327 | 9 184 +- 112 | 62 969 | 4 784 +- 68 (0.076) | 5 991 | 239 +- 7 |
| PEAK | 111 270 | 7 150 +- 92 | 56 532 | 4 066 +- 62 (0.072) | 5 483 | 191 +- 7 |
| vz150 | 103 430 | 4 824 +- 81 | 52 299 | 2 598 +- 49 (0.050) | 4 885 | 141 +- 6 |
| vz170 | 114 408 | 7 188 +- 96 | 57 289 | 3 727 +- 61 (0.065) | 5 439 | 194 +- 7 |
| none | 170 751 | 22 660 +- 179 | 86 011 | 11 488 +- 114 (0.134) | 8 616 | 641 +- 14 |
| S0_best | 91 899 | 5 133 +- 82 | 46 494 | 2 732 +- 48 (0.059) | 4 369 | 145 +- 6 |
| pre_best | 123 499 | 12 189 +- 131 | 62 759 | 6 308 +- 78 (0.101) | 6 125 | 361 +- 9 |
| S0_corfb | 112 791 | 6 336 +- 91 | 56 672 | 3 318 +- 57 (0.059) | 5 315 | 174 +- 6 |
| pre_corfb | 156 201 | 15 864 +- 165 | 78 860 | 8 099 +- 95 (0.103) | 7 677 | 467 +- 11 |

- **Cross-checks with the other teams' notes** (exact): pre B 153 183 and P 77 513; dz40 P 62 969; PEAK 56 532; pre_best P 62 759;
  none P 86 011 with CB/N 0.134; X17 S0_best 4 369; pre_corfb rows 4 449 982.
- **none, bump window:** W = N - CB = **148 091.1 +- 466.6**, with N = 170 751 and CB = 22 659.9.
- **S0, X17 counting window X17c** (spec: M4 0.5312-0.5566, M(ee) 15.5-24.5 MeV): N = 6 424, CB = 209.5 +- 7.1.
- **CB spread rule vs per-cell geometric mean**, data S0, values of 2 sqrt(N++ N--) / CB(spread):

  | cells | N--/N++ | 2 sqrt(N++ N--) / CB(spread) |
  |---|---|---|
  | F10 E1 (S1-S5) | 7.0-14.0 | 0.88 / 0.91 / 0.80 / 0.79 / 0.64 |
  | F10 E2 (S1-S5) | 3.6-6.3 | 1.12 / 1.09 / 1.06 / 1.01 / 0.87 |
  | windows B / P / CR-U | | 1.003 / 1.001 / 1.003 |
  | X17c / X17P | | 0.64 / 0.65 |
  | CR-K | 75 (N++ = 2) | 0.17 |

- **CB error:** replica RMS / linear error (median over F10 bins): E1 3.2-4.2, E2 1.2-1.6. The relative CB error reaches
  0.30-0.47 in S4/S5 at 0.40-0.50, where the bands hold few like-sign rows.
- No data row lacks an F10 cell, and no CB band is empty in any family.

---------------------------------------------------------------------------------------------------
## 3. Figures and tables (draft captions)

- `docs/latex/figs_all/wp0_definitions.pdf` (.png). *Definitions of the analysis plane M(pi+pi-e+e-) x M(pi+pi-).*
  - Grey scale: the simulated bump C2+C3 (pi0 Dalitz and conversion, photon lost).
    - Selection S0. SMASH, mode A, w x V3.
    - Mass transfer: M4' = 0.99615 M4, M(pipi)' = 0.9952 M(pipi).
    - Units: ub/(GeV/c2)^2, 10 MeV cells.
  - Slices S1-S5: dotted lines at 0.34, 0.37, 0.39, 0.43.
  - Windows: bump window B (blue) 0.4550-0.5225; peak window P (red) 0.53364-0.55416; CR-U (green) 0.565-0.700;
    CR-K (magenta) S5 x 0.400-0.5225; CR-L (orange, test only) 0.5225-0.5300.
  - Edges: E(m, 0) (solid) and E(m, 60 MeV) (dashed); no bump above m = 0.4116.
  - Simulation only, no data. Simulation units: L_sim = 221.05 nb^-1.
- `docs/latex/figs_all/wp0_cb_errors.pdf`. *Relative uncertainty of the combinatorial background per F10 cell.*
  - Data S0, raw rows, isBest >= 0.
  - CB = 2 sqrt(N++ N--) in 10 MeV cells, spread by the like-sign share of the 50 MeV band.
  - Red: RMS of 200 Poisson event-group replicas of the like-sign rows. Black dashed: linear propagation of the 10 MeV
    cell counts only.
  - The difference is the statistical error of the cell fractions. It dominates in E1 (M(ee) < 10 MeV) and in S4/S5
    below 0.50 GeV/c2.
- `docs/latex/figs_all/wp0_sim_mcerr.pdf`. *Relative MC error of the simulated sum T1+T2 per F10 cell.*
  - S0, SMASH std + lep + eta, mode A, w x V3, M4' = 0.99615 M4.
  - Red: RMS of 100 Poisson event-group replicas. Black dashed: event-group error.
  - The two agree (median ratio 0.99-1.02).
- `docs/latex/tables_all/wp0_acceptance.tex`: acceptance criteria A1-A4 with values.
- `docs/latex/tables_all/wp0_sim_totals.{txt,tex}`: simulated cross sections per directory.
  - Columns: T12 all M4, C2+C3 in B, C4/C5/G_eta in P, C6 per BR in X17c.
  - Units: nb, L_sim = 221.05 nb^-1; event-group errors.
- `docs/latex/tables_all/wp0_numbers.txt`: the factory's own record.
  - Data: window counts, CB (linear and replica errors) and N - CB for all 11 selections.
  - Simulation: every directory and component, with the total and 7 windows.
  - All checks.
- `docs/latex/tables_all/wp0_check.txt`: numbers read back from the file, the CB comparison and the replica validation.

---------------------------------------------------------------------------------------------------
## 4. Acceptance criteria (spec WP0)

| # | criterion | value | verdict |
|---|---|---|---|
| 1 | sim pre T1+T2 = 701.1 nb | 701.0832 +- 1.4624 nb | **PASS** |
| 1 | C2+C3 in the sim bump window = 9.689 nb (s = 1) | 9.6887 nb | **PASS** |
| 2 | data S0 peak window N = 55947 | 55 947 | **PASS** |
| 2 | data S0 peak window CB/N = 0.059 | 0.0586 (CB 3 279.0 +- 56.6) | **PASS** |
| 3 | data none bump-window W = 148091 | 148 091.1 +- 466.6 | **PASS** |
| 4 | sum of the cells = 1D to 1e-9 | data 3.5e-15, sim 2.9e-12 (worst of all blocks and layouts); read-back F10: 2.7e-15 over 1056 pairs | **PASS** |

**Edge function:**
- E(0.30, 0) = 0.51859, E(0.39, 0) = 0.50077, E(0.39, 60 MeV) = 0.50958, E(0.4116, 0) = 0.47390;
- E(0.30, 100 MeV) = 0.5326 (spec "0.533").

All 28 internal checks pass.

**Fall-backs:** none of the spec's fall-backs concerns WP0, and none fired.

**Deviations from the spec** (the closest sound choice in each case; all flagged here and in README.md):

| # | deviation | reason |
|---|---|---|
| a | CR-K upper edge = 0.5225 | The spec writes "< 0.522"; that is not an F10 bin edge, and the nearest edge is 0.5225. |
| b | Extra window X17P = P x M(ee) 15.5-24.5 MeV | The count 5258 quoted in D19 belongs to this window. The spec's counting window X17c (0.5312-0.5566) holds 6424. |
| c | Extra layout G3 (M4 0.40-0.85 in 2.5 MeV x 27 M(pipi) bins x 5 M(ee) classes) | Not in the spec, but F10 alone cannot serve S2 (CR-K at 0.425/0.44, CR-U edges 0.560/0.575, upper 0.75), S5 (2.5 MeV bins, range 0.40-0.75, slice edges +-5 MeV, M(ee) split 8/12 MeV), A4 (CR-U split at 0.42), A5 (K0S veto 0.485-0.507) or the S3 "range 0.85". It is filled for data (all) and S0__nom, pre__nom. |
| d | BEFORE variables and the map only in the sim directories nom, ghostA and wtsOff (S0 and pre) | Size and speed. Every directory has B_M4, fine, F10, F5g, X2D, win and tot. |
| e | Simulation templates stored in nb, not in expected raw rows | Conversion: x `alldefs::kNbToRaw` = 6079.2 x 63/64. |
| f | The 200 CB replicas are Poisson event-group replicas, not row-level Poisson | Rows of one event are correlated (1.2 rows per event). |
| g | pre treated as a "full parallel column" (D1) | It has every variant, halves, replicas, G3 and BEFORE, like S0. |

---------------------------------------------------------------------------------------------------
## 5. Caveats

1. **CB construction matters per cell.**
   - The spread rule assumes that the N--/N++ ratio of a 50 MeV band holds in every cell. In data it does not:
     N--/N++ is 7-18 at small M(ee) and in the X17 region, and 3.6-6.3 at M(ee) >= 10 MeV.
   - Over whole windows the two constructions agree (0.3 %). In single cells they do not: 0.64-1.12, X17 windows 0.65,
     CR-K 0.17.
   - Consequences:
     - k_E2 = 1 +- 0.05 may be too narrow for single slices (S1E2 1.12, S5E2 0.87).
     - The CR-K null test (A2 of WP3) and the X17 background normalisation (WP9: k CB_ij) depend strongly on the choice.
     - D19's "CB 112" is the geometric mean of the window; the spec's rule gives 171 +- 6.
   - The S1 family must include this variant. The per-cell geometric mean is computable from the stored Npp/Nmm.
2. **CB errors.** The stored CB bin errors exclude the error of f. Use `CBrep_*` (factor 1.2-4.2 larger).
3. **Blinding.**
   - Coarse M(ee) information is NOT masked: F10 E2, F5g, G3 classes 12-20 and >= 20 MeV, the X17c/X17P window counts.
     These are coarser than the 2 MeV shape, and their content is the already-known window count and peak share.
   - During development I opened the sealed TEST file once, built from the first 200 000 e+e- rows only, to verify that
     the copy works. It was one X2D cell (8 rows). The final sealed file was never opened.
4. **Simulation and the vertex cut.** The simulation has 0.15 % of its rows below vz = -160 mm (data: 7 %), so sim vz170
   differs from S0 by 0.05 % only. The vz variants therefore test the data, not the simulation.
5. **Sim MM4 and y4** are recomputed at 4500 MeV everywhere (BEFORE variables and PEAK). The data stay at 4530 MeV. The
   two MM4 spectra differ by the beam energy, which is physics; say so in the B4 captions.
6. **Data units.** Data histograms are raw rows: no 64/63, no L. Errors of N, N++, N-- are event-group errors; for Poisson
   fits use the contents.
7. **N - CB per window** is printed in wp0_numbers.txt (as for the A3 check). These are window contents including all
   background, not signal yields (D20).
8. **sim_eta_physics.h reads** rw_generator.root and rw_angles_eta.root from the simulation production directory. Angles
   were found for all 3 432 684 class-5 lookups.

---------------------------------------------------------------------------------------------------
## 6. Exact commands (from /home/przygoda/HADES/FAT/FAT)

    root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_wp0", kTRUE)' 'plots/all_report/wp0_fill.C+O(-1)' \
         > /hdd3/przygoda/HADES/FAT_out/all_report/wp0/wp0_fill.log 2>&1      # 157 s, 7.4 GB; argument N > 0 = test on N rows per tree
    root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_wp0", kTRUE)' 'plots/all_report/wp0_check.C+O'

Before a rerun: `chmod u+w` on hists_v1.root and the sealed file. Log of the final run:
`/hdd3/przygoda/HADES/FAT_out/all_report/wp0/wp0_fill_v1.log`.

## 7. Naming scheme (short; full text in plots/all_report/README.md)

    data/<sel>/<sub>/<what>_<layout>
       sel  S0 pre dz40 PEAK vz150 vz170 none S0_best pre_best S0_corfb pre_corfb
       sub  all fold0 fold1 perI perIIA perIIB perIII
       what N Npp Nmm CB | CB20 CBrep (all only)          [raw rows]
    sim/<sel>__<var>/<comp>_<layout>[_hA|_hB|_rep]       [nb per bin]
       sel__var  S0__{nom,V2,bump1aY,bumpPS,convOff,GS,wtsOff,modeB,wdb,ghostA,sK0,sUp,sDn,r089,r110}, pre__{same},
                 dz40__nom PEAK__nom vz150__nom vz170__nom
       comp C2 C3 C4 C5 C6 Geta Gpi0D(_om,_nom) Gconv(_om,_nom) Gvec(_om,_nom) Goth(_om,_nom) GD T12 T3 T4
    layout F10_S<1-5>E<1-2>  F5g_S<1-5>  X2D  G3  fine_M4  B_M4  B_map  B_<mpipi|mee|oa|vz|mm4|pt4|y4>_R<0-7>  win  tot
    win bins (alldefs::WinId): B P SBlo SBhi X17c CRU CRK CRL R0..R7 Jlo Jhi Fit Bef X17P
    blinded: data N of X2D M(ee) 12-28 MeV and B_mee_R2..R4 12-28 MeV = -1 (masks data/mask_X2D, data/mask_B_mee)

Reader: `#include "hists_io.h"`, e.g. `allio::Data(f, "S0", "all", "CB", allio::F10(2, 1))`,
`allio::Sim(f, "S0", "nom", alldefs::kC5, "X2D", "_rep")`.

---------------------------------------------------------------------------------------------------
## Verification (independent verifier, prefix v0, 2026-11-14)

**Verdict: WP0 PASSES.** I found no numerical error. 189 quantities were recomputed independently from the compact files,
and all 189 agree with hists_v1.root: integer counts exactly, CB to 1e-9 relative, simulation to 1e-9 relative. All four
acceptance criteria were re-derived and pass. Two comments were corrected (no value changed, hists_v1.root untouched).
Caveat 1 of WP0 (the CB spread rule) is confirmed and is the main open item for WP3/WP5/WP9.

### V1. What was checked, step by step
1. **Copies vs originals** (diff against FAT_simwt commit 7206630).
   - sim_eta_physics.h, s5_blind_ana_fit.h, s5_blind_ana_tmpl.h, sim_truth.h and sim_ghostvar.h are byte-identical (md5 as in
     README.md).
   - SimEngine.h has exactly the three announced changes:
     - include of the local sim_truth.h;
     - kFigDir/kTabDir moved to /hdd3/.../all_report/wp0/{figs,tables};
     - a comment at kLdata.
   - The copied headers write nothing else: the only other output, `SaveFig`, writes to a path the caller gives. The reads
     (rw_generator.root, rw_angles_eta.root, blind templates) are read-only on /hdd3.
   - The data-branch plots/ScanEngine.h that SimEngine.h now includes differs from the sim-branch version only by the
     scan_v2 paths, the three added columns (pip_z, pim_z, dz_pi) and a loader that tolerates their absence. It changes no
     constant and no function used by SimEngine.
   - ScanEngine.h and scan_*.C were not touched by WP0 (time stamps 04:13/04:14, before WP0 started at 04:56).
2. **defs.h against the spec, line by line.** Everything below agrees with the spec:
   - pre and S0, and every variant: dz40 and PEAK on pre; vz150/170 on S0; none = pre without vz > -160.
   - Data PT3: the production holds only trigbit 8192 (main.cc / setup_cuts.h). Sim PT3: (trigbit & 8192) != 0.
   - Windows: B, P, the side bands, X17c (= sim 0.5333-0.5587 x 0.99615 = 0.53125-0.55655), CR-U, CR-L, R0-R7, Jlo/Jhi.
   - CR-K upper edge 0.5225 instead of "0.522" (deviation a). I agree: 0.522 is not a bin edge, and 0.5225 is the edge of B.
   - Layouts:
     - F10: 24 + 20 + 26 = 70 bins.
     - X2D: 35 x 30.
     - G3: 27 M(pipi) bins, containing every S2/S5/A4/A5 edge; 5 M(ee) classes.
   - Slices with ">=" edges; E1/E2 split at 10 MeV; tilt x = (M4 - 0.55)/0.1.
   - Mass transfer: s = 0.99615, r = 1; variants 0.99521, 0.99705, 0.99525, r = 0.89/1.10; M(pipi) x 0.9952; M(ee) not scaled.
   - Physics-weight codes decoded bit by bit:
     - 108 = V3 + KLOE + conv + bkg;
     - 107 = V2; 116 = 1+aY; 100 = bump off; 76 = conversion F_V off; 236 = + Gounaris-Sakurai.
   - BkgGroup and the ghost rule (a) are identical to s5_blind_prep.C and sim_ghostvar.h (RelabelTmpl v = 1).
   - Units:
     - kNbToRaw = 6079.2 x 63/64 = 5984.21;
     - kSimToRaw = 6079.2/221.05 x 63/64 = 27.0718 (= the spec's "sum w x L/L_sim x 63/64");
     - kRawToNb = (64/63)/6079.2;
     - old-table factor (64/63)(6080.5/6079.2) = 1.01609.
   - Periods by day: the mapping matches the "period" column of lumi_startfree_perday.csv (034-046 I, 047-048 IIA,
     049-053 IIB, 059-067 III). kLper equals sf_luminosity_perperiod.csv.
   - **Edge E(m, q)**, recomputed in Python with m_eta = 0.547862 and m_pi0 = 0.1349768:
     - 0.51859 (0.30, 0), 0.50077 (0.39, 0), 0.50958 (0.39, 60 MeV), 0.47390 (0.4116, 0), 0.53262 (0.30, 100 MeV);
     - spec: 0.5186, 0.5008, 0.5096, 0.4739, 0.533. PASS.
     - The limit (m_eta - m_pi0) x 0.9970 = 0.41165.
3. **Own macro `plots/all_report/v0_verify.C`.** It uses neither defs.h nor hists_io.h nor the WP0 copies.
   - **Data:** read with raw TTree branches; my own selections, windows and slices.
   - **CB:** my own code. 2 sqrt(N++ N--) in 10 MeV cells from 0.30, linear interpolation between the cell centres, then
     fractional fine bins for the windows. For a cell, the CB is multiplied by f(cell | 50 MeV band) from N++ + N--.
   - **Simulation:** read through the ORIGINAL FAT_simwt headers by absolute include, streaming the rows (nothing loaded into
     memory).
     - My own ownership and weights.
     - My own M4' and r-variant.
     - My own MM4 at 4500 MeV, computed with TLorentzVector.
     - My own event-group variance, with a global map keyed by (tree, grp). This does not rely on adjacent rows.
   - **The file itself:** closure, component sums, halves and replicas.
   - Table: `docs/latex/tables_all/v0_verify.txt`.

### V2. Numbers (independent value = file value in every line; selection isBest >= 0 unless stated)
| quantity | value |
|---|---|
| data S0: N / N++ / N-- in P | 55 947 / 691 / 3 894 |
| data S0: CB in P (CB/N) | 3 279.024 (0.05861); replicas 56.6, linear 58.9 |
| data S0: N, CB in B | 111 212, 6 236.84 |
| data S0: N, CB in CR-U | 722 506, 37 078.72 |
| data S0: N, CB in CR-K (N++ = 2, N-- = 150) | 964, 199.52 (geometric mean of the window: 34.6) |
| data S0: N, CB in X17P (window of D19) | 5 258, 171.03 (geometric mean 111.7 = D19's "112") |
| data S0: N, CB in X17c | 6 424, 209.50 |
| data S0: F10_S2E2, M4 0.540-0.5425 | N 772, N++ 14, N-- 83, CB 79.030 |
| data S0: X2D, M4 0.5425-0.545 x M(ee) 30-32 MeV (unmasked) | N 84, N++ 2 |
| data S0: N in P, period IIA / fold0 | 3 505 / 27 892 |
| data S0: event-group error of N in P | 239.067 |
| data none: W = N - CB in B | 170 751 - 22 659.94 = **148 091.1** |
| data pre: N in B / P | 153 183 / 77 513 |
| data S0_best: N in P | 46 494 |
| data S0_corfb: N in P | 56 672 |
| sim pre: T1+T2, all M4 | **701.0832 nb** |
| sim pre: C2+C3 in the sim-scale window 0.4576-0.5255 | **9.6887 nb** |
| sim S0__nom: C5 in P | 1.046212 +- 0.002949 nb (event-group error identical) |
| sim S0__nom: C2+C3 in B | 9.00758 +- 0.22452 nb |
| sim S0__nom: G_eta in P | 0.610206 nb |
| sim S0__nom: T12 in CR-U | 90.44899 +- 0.54830 nb |
| sim S0__nom: G_conv in CR-U | 34.28537 nb |
| sim S0__nom: T12 in F10_S3E1, 0.5425-0.545 | 0.090849 nb |
| sim S0 variants: T12, all M4 | modeB 650.3789; ghostA 716.1308; w_db 631.2007 nb |
| sim PEAK: T12, all M4 | 628.2276 nb |
| sim S0: C5 in P, r = 0.89 / V2 | 1.091711 / 1.214721 nb |
| sim S0: C6 per unit BR in X17c | 3 501.597 nb |

Further checks:
- **Data structure.**
  - Event groups are adjacent in all three data trees: no group reappears after another group.
  - The dz_pi column equals max(|pip_z - vz|, |pim_z - vz|) in all 5.28 M rows.
  - No selected row has M4 < 0.30. So R0 and tot, whose N has no lower M4 limit while the CB grid starts at 0.30, are
    consistent.
- **Blinding.**
  - In all 11 selections x 7 subsets, the X2D bins with M(ee) 12-28 MeV are exactly -1, and no other bin is.
  - The same holds for B_mee_R2..R4: 85 470 cells checked, 0 wrong.
- **Closure.**
  - Sum of the 10 F10 cells = fine 0.40-0.70: data N 976 545, data CB 51 691.336, sim T12 131.5883 nb.
  - Simulation T12 = C2..C5 + G_eta + G_pi0D + G_conv + G_vec + G_oth; G_D = G_pi0D + G_vec + G_oth;
    G_conv = _om + _nom (to 1e-9).
- **Halves.** C5 in P: half A + half B = full; A/B = 0.994.
- **Replicas.**
  - **CB replicas:** the means are within 2.5 sigma_mean of the nominal (P +2.5, B +0.1, CR-K -1.1, X17P +1.4; the windows
    are correlated). The RMS / linear-error ratios are 0.96 (P), 0.98 (B), 2.10 (CR-K) and 2.00 (X17P).
  - **Generator:** I tested it standalone on 1.6e8 draws: mean 1.000016 +- 0.000079, variance 0.999985. The +2.5 is
    therefore a fluctuation, not a bias.
  - **Sim C5 replicas** in P: RMS / event-group error = 0.91 (100 replicas: +-7 %).

### V3. Acceptance criteria of WP0 (re-derived independently)
| # | criterion | independent value | verdict |
|---|---|---|---|
| 1 | sim pre T1+T2 = 701.1 nb; C2+C3 in the sim bump window = 9.689 nb | 701.0832; 9.6887 (window 0.4576-0.5255 of the sim study, s = 1) | PASS |
| 2 | data S0 P: N = 55947, CB/N 0.059 | 55 947; 0.05861 | PASS |
| 3 | data none B: W = 148091 | 148 091.1 | PASS |
| 4 | cells = 1D to 1e-9 | data N exact, CB 3e-16 relative, sim 2e-16 | PASS |

**Note on criterion 1.** The window 0.4576-0.5255 is the sim study's own bump window: it equals B/0.9943, not B/s =
0.4568-0.5245. For a reproduction check this is the right choice. With the spec's s = 0.99615 on the data-scale B, pre
C2+C3 = 9.754 nb.

**Fall-backs:** none fired; none concerns WP0.

### V4. Fixes made in place (comments only; no value, no histogram changed)
1. **defs.h, comment at kL.** The old wording "sector-2 deficit of days 034-044 already removed in L" was incomplete.
   - It is true that the elastic L contains the removal (factor S2fac of lumi_startfree_perday.csv), so L is unbiased.
   - But the same deficit still lowers the analysis' own PT3 yields. Those need C_sec (spec D16, WP2; A3: "NOT in L: handle
     it in the analysis' own trigger efficiency").
   - The new comment says both. Otherwise a WP8 agent could skip C_sec.
   - Note: the orchestrator's sentence "L already contains the sector-2 deficit removal" is consistent with this reading.
2. **SimEngine.h, comment at kLdata.** It named a non-existent `alldefs::kSimToRawRows`. It now names `kNbToRaw` and
   `kSimToRaw`.
   - The README's description of the change was already right.
   - Syntax check of defs.h + hists_io.h + sim_eta_physics.h after the edits: clean.

### V5. Findings and caveats for the consumers (ranked)
1. **The CB spread rule has a known bias direction** (confirms WP0 caveat 1).
   - The cause: f(cell | band) is built from N++ + N--. It therefore overestimates the CB where N--/N++ is above the band
     average, and underestimates it where N--/N++ is below.
     - Above the average: E1, the X17 region, CR-K (N--/N++ = 75).
     - Below the average: E2 in S1-S3.
   - Size, as 2 sqrt(N++ N--) / CB(spread):
     - CR-K 0.17: CB 199.5 = 20.7 % of N, vs 34.6 = 3.6 % of N. This decides the WP3 A2 null test.
     - X17P 0.65.
     - E1 0.64-0.91.
     - E2 0.87-1.12.
   - Whole windows B, P and CR-U agree to 0.3 %.
   - Recommendations:
     - WP3 should quote A2 with both constructions.
     - WP5/WP9 should carry the S1 variant "per-cell geometric mean" (the spec already lists it).
     - k_E2 = 1 +- 0.05 is narrower than the per-slice spread (+-12 %); the free c_j absorb part of it.
2. **Soft blinding (D19 is fulfilled literally).** Two unmasked data histograms still carry X17-sensitive information in
   the peak band:
   - the G3 class 12-20 MeV in 2.5 MeV M4 bins (S0: 8 888 rows in M4 0.5325-0.555);
   - B_oa_R2..R4 in 1-degree bins (an X17 at 17 MeV sits at OA of about 5-15 degrees).

   Agents should not form these in data before x17_frozen.txt. Masking them in a v2 is the analyst's call.
3. **CB errors outside "all".** The folds, the periods and G3 have only the linear CB error, without the error of f, and no
   replicas. In E1, X17 and CR-K cells this underestimates the CB error by a factor of 2-4. Fold/period consistency tests
   should inflate it, or use the "all" replica ratio.
4. **Simulation units.** The simulation is stored in nb (deviation e). Raw rows = x kNbToRaw = 5984.21. This is
   documented; a missed conversion is a factor of 6000, so it is obvious.
5. **Cosmetic.** The title of wp0_definitions.pdf renders "C _2+C_3" with a stray space. The figure content (windows,
   slices, E(m, 0) and E(m, 60 MeV), CR boxes) is correct.

### V6. Files and commands
- Macro: `/home/przygoda/HADES/FAT/FAT/plots/all_report/v0_verify.C` (argument: 1 data, 2 sim, 4 file; 7 = all; 20 s, 0.8 GB).
- Table: `/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all/v0_verify.txt` (189 OK, 0 MISMATCH).
- Logs and the generator test: `/hdd3/przygoda/HADES/FAT_out/all_report/v0/`.

      cd /home/przygoda/HADES/FAT/FAT
      root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_v0", kTRUE)' 'plots/all_report/v0_verify.C+O(7)'
