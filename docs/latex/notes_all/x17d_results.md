# WP9 X-d (X17 unblinding) -- results note (prefix x17d)

Date 2026-11-14. Branch pp45_pippimepem (worktree /home/przygoda/HADES/FAT/FAT). Nothing was committed and no git command
that changes state was run. I did not touch ScanEngine.h, scan_*.C, hists_v1.root, the sealed file, x17_frozen.txt,
x17c_config.txt or any file of WP0-WP6, x17ab or x17c. I only included those files and ran no macro of theirs.
After the run the md5 sums of x17c_model.h, x17c_lib.h, x17c_*.C, defs.h and hists_io.h are still the ones listed in
x17_frozen.txt.

**Precondition.** `/hdd3/przygoda/HADES/FAT_out/all_report/x17/x17_frozen.txt` (md5 60a5376a...) existed before I opened
the sealed file. Every md5 that it lists (code, hists_v1.root, x17_lines.root, wp5/results.txt, x17c_egcov.root) matched.

**Deviation from the frozen procedure: none.** The frozen text left four details open. I filled them as described in
Sect. 5 (items 1-4); none of them changes the primary limit.

----------------------------------------------------------------------------------------------------------------------
## Headline (S0, isBest >= 0, L = 6079.2 nb^-1, isotropic e+ angle)

90 % CLs upper limits on BR(eta -> pi+pi- X17) x BR(X17 -> e+e-):

| m_X [MeV/c2] | normalisation (1), PRIMARY | toy CLs on the data | frozen expectation: median [16 %, 84 %] | local Z |
|---|---|---|---|---|
| 16.70 | **3.06e-5** | 2.93e-5 | 2.64e-5 [1.52e-5, 4.84e-5] | 0.46 |
| 16.85 | **3.06e-5** | 3.07e-5 | 2.62e-5 [1.51e-5, 4.82e-5] | 0.48 |
| 17.00 | **3.06e-5** | 3.01e-5 | 2.60e-5 [1.50e-5, 4.77e-5] | 0.49 |

- **No signal.** The largest local significance on the 10-40 MeV/c2 grid is Z = 0.98 (p0 0.16), at 22 MeV/c2.
- At all 38 masses the observed limit lies inside the frozen +-1 sigma band. Relative to the frozen Asimov median it is
  0.74-1.38 times as large.
- At 17 MeV/c2, 61 % of the 2000 frozen background-only toys have a lower limit than the data.

----------------------------------------------------------------------------------------------------------------------
## 1. What was done, step by step (plain words)

1. **Check that nothing changed since the freeze.**
   - I compared the md5 sums of all frozen code and inputs with x17_frozen.txt Sect. 1: all equal.
   - The configuration is read with the frozen reader (x17cl::ReadConfig on x17c_config.txt).
2. **Full data set** (frozen Sect. 10, step 1).
   - Side-band cells (M(ee) 2-12 and 28-62 MeV/c2): taken from hists_v1.root with the frozen loader, which also checks that
     the mask is intact.
   - The 8 formerly masked rows (12-28 MeV/c2, all 35 M4 bins): taken from the sealed file
     `sealed/hists_v1_x17_sealed.root`, `data/S0/all/N_X2D`.
   - The sealed file holds the side-band cells as well. They are identical to hists_v1.root (largest difference 0), and
     the binning is identical.
   - Totals: 255 331 rows in the 1050 cells, of which 197 752 are side-band rows (as in X-c) and 57 579 are unblinded rows.
3. **Reproduction of the freeze** (my own check, done before any limit).
   - I reran the frozen side-band fit on the side-band cells only.
   - It reproduces B_pred(X17c) = 6311.7 rows (frozen 6311.7) and the Asimov median at 17 MeV/c2, 2.7589e-5 (frozen 2.7590e-5).
   - So inputs, code and configuration are those of the freeze.
4. **Full fit and goodness of fit** (frozen Sect. 10, step 4).
   - The frozen model is fitted to all 1050 cells with r = 0, and again with r free.
   - The Pearson chi2 is compared with the frozen fall-back rule (qScale is replaced by chi2/ndf if chi2/ndf > 1.5).
   - Per M(ee) row I record: the chi2, and the residual N - mu in the peak band (M4 0.5325-0.5575 GeV/c2).
   - The fitted background n_j of the unblinded rows is compared with the frozen prediction.
5. **Mass scan** (frozen Sect. 10, step 2). At each of the 38 masses of the frozen grid, and for each of normalisations
   (1), (2), (3):
   - the line is set with x17cl::SetLine;
   - the observed limit and r_hat come from x17cl::LimitBR (asymptotic CLs, q~_r divided by qScale = 1.0299);
   - the local p0 is computed as q0 = (F(r = 0) - F(r_hat)) / qScale for r_hat > 0, Z = sqrt(q0);
   - the e+ decay angle 1 + cos^2 limit is computed for normalisation (1), and quoted separately.
   - The Asimov bands built from the OBSERVED r = 0 fit come out of the same call. They are stored as information only. The
     expected bands of the report stay the FROZEN ones (x17c_state.root, and the x17c toys).
6. **Toy CLs on the data** (frozen Sect. 6), with the frozen x17c::ToyCLs at 16.7, 16.85 and 17.0 MeV/c2, all three
   normalisations.
   - 9 values of r per mass: 0.55-1.8 x the observed asymptotic limit.
   - 1000 toys per point and per hypothesis.
   - Toy UL from the frozen x17c::ToyUL (quadratic fit of ln CLs).
   - Frozen rule: the toy UL is primary if it differs from the asymptotic UL by more than 10 %.
7. **Feldman-Cousins cross-check** at 17 MeV/c2, normalisation (1).
   - Test statistic: t_r = (F(r) - F(max(r_hat, 0))) / qScale. This is q~ without the one-sided zeroing; the boundary
     r >= 0 is kept.
   - For each of 14 values of r (BR 0 to 3.75e-5): pseudo-data from the fit of the DATA at fixed r, over-dispersed by
     qScale, with the constraint centres drawn, exactly as in the frozen toys. 2000 toys per value.
   - t_c(r) = the 90 % quantile of the toy t_r. The interval is the set of r with t_obs(r) <= t_c(r). Its edges come from
     a linear interpolation of t_obs - t_c between grid points.
8. **Template statistics, 'MC'** (frozen Sect. 4 and 10, step 5).
   - The observed limit is recomputed with each of the 200 CB replicas (CBrep_X2D) and each of the 100 C4 replicas
     (C4_X2D_rep), at the three quoted masses and in all normalisations.
   - The RMS of the limit is quoted separately (D13).
9. **Counting cross-check in X17c** (frozen Sect. 6): N(X17c) - B_pred, with the frozen B_pred = 6311.7 and
   sigma_B = 1195.5. Then BR = excess / (a_eta s1), with a_eta = 0.3823 and s1 = 2.0954e7.
10. **Look-elsewhere.** A global p-value from toys is needed only if the largest local Z exceeds 3. It is 0.98, so this step
    was not run (frozen rule).
11. **Information only, not part of the limit.** The observed limit at 17 MeV/c2 under the 18 model variants that X-c used
    to size the frozen nuisances. These choices are already inside the frozen limit through the N5-bias and sX; the
    numbers are given for the UL_X17 column of Table S1 (WP7).
12. **Figures X1-X4 and Table D1.**

----------------------------------------------------------------------------------------------------------------------
## 2. Numbers (S0, isBest >= 0, raw rows; BR = BR(eta -> pi+pi-X17) BR(X17 -> e+e-); L = 6079.2 nb^-1)

### 2.1 Full fit (all 1050 cells, normalisation (1)); table `x17d_fit.txt`
**Goodness of fit.**

| fit | Pearson chi2 / ndf | value |
|---|---|---|
| r = 0 | 994.9 / 1002 | **0.993** |
| r free, 17 MeV/c2 | 994.7 / 1001 | 0.994 |
| unblinded cells alone (r = 0) | 245.0 / 280 cells | -- |
| side-band cells (r = 0) | 749.9 / 770 cells | -- |

The frozen fall-back (chi2/ndf > 1.5) did not fire, so qScale stays at 1.0299.

**Parameters of the r = 0 fit.**

| parameter | full fit | side-band fit of X-c |
|---|---|---|
| N5 | 3659 +- 503 (13.7 %) | 3380 +- 686 |
| N4 | 1648 +- 383 | 1803 |
| k_E1 | 1.04 +- 0.50 | |
| k_E2 | 0.998 +- 0.050 | |
| a_eta | 0.253 +- 0.138 (prior 0.382 +- 0.165) | |
| g | 1.48 +- 0.93 | |
| thC | -0.16 +- 0.97 | |

**Peak-band residual of the r = 0 fit** (data minus model, M4 0.5325-0.5575 GeV/c2):

| M(ee) range [MeV/c2] | N - mu [rows] | significance |
|---|---|---|
| 12-28 | +44 +- 117 | +0.38 sigma |
| 14-26 | +71 +- 95 | +0.74 sigma |

- Per row, |pull| <= 0.5 in 12-28 MeV/c2.
- The largest pull anywhere is -2.05 sigma, at 44-46 MeV/c2 (side band).

**Background n_j of the unblinded rows against the frozen prediction.**
- The sum over the 8 rows agrees: 53 126 fitted against 53 132 predicted (ratio 0.9999).
- Row by row the ratio is 0.81 (26-28 MeV/c2) to 1.10 (16-18 MeV/c2). These deviations follow the statistical
  fluctuations of the simulated reference rows: raw/smoothed simulation is 0.82 at 16-18 and 1.14 at 26-28 MeV/c2 (x17c
  freeze log).
- The observed 2D limit does not use the prediction, because n_j is free in every row (D18). Only the frozen expected
  limit and B_pred(X17c) depend on it (caveat 2).

### 2.2 Observed limits versus mass; table `x17d_scan.txt` (38 masses), figure X2

| m [MeV/c2] | UL (1) | UL (2) | UL (3) | r_hat (1) +- Hessian | local Z (1) | frozen Asimov median | UL (1), e+ angle 1 + cos^2 |
|---|---|---|---|---|---|---|---|
| 10 | 3.89e-5 | 4.77e-5 | 4.08e-5 | -0.83 +- 2.12 e-5 | 0 | 5.24e-5 | 4.66e-5 |
| 12 | 3.26e-5 | 3.86e-5 | 3.29e-5 | -0.26 +- 1.59 e-5 | 0 | 4.01e-5 | 3.90e-5 |
| 15 | 3.06e-5 | 3.41e-5 | 2.87e-5 | +0.31 +- 1.26 e-5 | 0.25 | 3.12e-5 | 3.67e-5 |
| 16.7 | 3.06e-5 | 3.29e-5 | 2.76e-5 | +0.52 +- 1.18 e-5 | 0.46 | 2.81e-5 | 3.67e-5 |
| 16.85 | 3.06e-5 | 3.29e-5 | 2.76e-5 | +0.54 +- 1.17 e-5 | 0.48 | 2.78e-5 | 3.67e-5 |
| 17.0 | 3.06e-5 | 3.28e-5 | 2.75e-5 | +0.55 +- 1.17 e-5 | 0.49 | 2.76e-5 | 3.67e-5 |
| 18 | 3.08e-5 | 3.23e-5 | 2.70e-5 | +0.66 +- 1.14 e-5 | 0.62 | 2.61e-5 | 3.69e-5 |
| 20 | 3.06e-5 | 3.13e-5 | 2.61e-5 | +0.82 +- 1.08 e-5 | 0.84 | 2.33e-5 | 3.67e-5 |
| 22 | 2.89e-5 | 2.94e-5 | 2.44e-5 | +0.86 +- 1.01 e-5 | **0.98** | 2.09e-5 | 3.47e-5 |
| 25 | 2.29e-5 | 2.41e-5 | 2.01e-5 | +0.61 +- 0.82 e-5 | 0.83 | 1.80e-5 | 2.74e-5 |
| 30 | 1.40e-5 | 1.58e-5 | 1.33e-5 | +0.12 +- 0.58 e-5 | 0.20 | 1.45e-5 | 1.67e-5 |
| 35 | 1.02e-5 | 1.20e-5 | 1.02e-5 | -0.07 +- 0.47 e-5 | 0 | 1.24e-5 | 1.22e-5 |
| 40 | 0.85e-5 | 1.04e-5 | 0.89e-5 | -0.05 +- 0.43 e-5 | 0 | 1.03e-5 | 1.02e-5 |

- All 38 observed limits of normalisation (1) are inside the frozen +-1 sigma Asimov band.
- The Asimov median of the OBSERVED r = 0 fit (information only) is 2.52e-5 at 17 MeV/c2. The frozen value is 2.76e-5; the
  difference comes from the masked-row prediction (caveat 2).

### 2.3 Toy CLs on the data (`toycls_m*_n*.txt`, `x17d_checks.txt`)

| m | toy/asymptotic UL, (1) | (2) | (3) |
|---|---|---|---|
| 16.70 | 0.958 | 1.008 | 1.011 |
| 16.85 | 1.002 | 1.034 | 1.018 |
| 17.00 | 0.983 | 1.041 | 1.013 |

- All ratios are within 10 %, so the **asymptotic limit stays primary** (frozen rule).
- CL_b > 0 at every point (CL_b 0.64-0.71).

### 2.4 Feldman-Cousins at 17 MeV/c2, normalisation (1)
- **90 % CL interval: BR in [0, 3.45e-5]** (lower edge at the boundary; t_obs(0) = 0.24 < t_c(0) = 1.56).
- t_c(r) is 1.56-2.85 (the asymptotic two-sided value is 2.71).
- The FC upper edge lies above the CLs limit (3.06e-5). This is expected: away from the boundary a two-sided 90 % interval
  has a 95 % upper edge, while the CLs limit is one-sided.

### 2.5 Local p0 and look-elsewhere
- **Largest local Z = 0.98 (p0 = 0.16) at 22 MeV/c2.** In the 16.7-17.0 MeV/c2 region Z = 0.46-0.49 (p0 0.31-0.32).
- **Empirical check at 17 MeV/c2** with the 2000 frozen background-only toys: the fraction of toys with
  r_hat/sigma >= 0.47 is 0.295, against an asymptotic p0 of 0.311. At 16.7 / 16.85: 0.307 / 0.302, against 0.323 / 0.317.
- Global p from toys: not needed (Z < 3).

### 2.6 Template statistics (MC, quoted separately; normalisation (1), 17 MeV/c2)

| replicas | RMS of the UL | RMS of r_hat |
|---|---|---|
| 200 CB replicas | 1.9e-6 | 1.2e-6 |
| 100 C4 replicas | 2.8e-6 | 1.7e-6 |
| combined | **3.4e-6 (11 %)** | |

- In normalisations (2) and (3) the combined RMS is 7.5 % and 7.2 %. 16.7 and 16.85 MeV/c2 give the same numbers.
- The mean limit over the C4 replicas is 2.88e-5, 6 % below the nominal value. The extra noise of a replica template lowers
  the limit (a template-attenuation effect).
- The MC error of r_hat (2.1e-6) is 0.18 of its statistical error (1.17e-5), so the result is not simulation-limited.

### 2.7 Counting cross-check in X17c (M4 0.5312-0.5566 GeV/c2, M(ee) 15.5-24.5 MeV/c2)
- **Data.**
  - N = 6465.6 rows, from the X2D cells with overlap fractions (the same rule as B_pred).
  - The exact row count is 6424, the value already known before the freeze.
  - CB 212.3 (overlap) / 209.5 (exact).
- **Prediction.** B_pred = 6311.7 +- 1195.5 (frozen).
- **Excess = +154 +- 1198 rows.**
  - BR = excess / (a_eta s1) = **(+1.9 +- 15.0) x 1e-5** (exact count: +112 rows, +1.4e-5).
  - The check is weak by construction (frozen statement), but it agrees with zero.

### 2.8 The three normalisations and the absolute limits at 17 MeV/c2 (Table D1)

| normalisation | 90 % CLs UL | N_X at the UL [raw rows] |
|---|---|---|
| (1) r x 2.68e-4, N5 fitted | **3.06e-5** | 563 |
| (2) N_X / (a_eta s_X2D), a_eta = 0.3823, s_X2D = 3.745e7 | 3.28e-5 | 469 |
| (3) relative to Y_peak, a_P = 0.4147 | 2.75e-5 | 427 |

- In normalisation (1), 563 raw rows correspond to 572 rows after the 64/63 PT3 correction.
- **sigma x BR x BR: PROVISIONAL.**
  - WP8 has not delivered sigma(pp -> eta X). I used the WP8 formula with the available inputs:
    sigma(pp -> eta X) = a_eta (WP5, 0.3823) x 1.0665 mb x C_sec (1 - f_win) (WP2, 0.842) = 0.343 mb.
  - This gives sigma(pp -> eta X) x BR_X x BR_ee < **10.5 nb** (90 % CL) at 16.7-17.0 MeV/c2.
  - Uncertainties: a_eta +-11 % (WP5 total; the frozen prior is +-43 %), C_sec (1 - f_win) +7/-6 %, L 6.9 %.
  - Information: the route through N_X and the simulated efficiency eps_X gives N_X^UL / (s_X2D / 1.0665 mb) x 0.842 =
    13.5 nb. It differs from 10.5 nb because it uses the fitted N5/N5ref = 0.46 as eta scale, instead of a_eta = 0.38.
  - When WP8 delivers sigma(pp -> eta X), the value scales linearly: sigma x BR x BR = 3.06e-5 x sigma(pp -> eta X).
- **Simulation expectation, for comparison only:** 1.9e-5 / 5.7e-6 (simulation study).

### 2.9 Information only: observed UL at 17 MeV/c2 under the X-c model variants (`x17d_variants.txt`)
Each entry is UL / nominal, with the nominal UL = 3.06e-5:

| family | variant: UL / nominal |
|---|---|
| background shape, simulated | exp(pol1) 1.006, exp(pol2) 1.004, exp(pol3) 0.863 (chi2/ndf 1.039), kappa_B 1.005 |
| background form | three classes 0.744; flat pol2 1.166; tilt only 1.154; no M(ee) dependence 1.156 |
| CB | factorised 0.983; per-row geometric mean 0.992 |
| CB factor | k_E1 free 1.018 |
| peak/direct templates | V2 1.184; GS 1.006; r = 0.89 / 1.10: 0.976 / 1.028; s up / down: 1.065 / 0.948 |

- The spread is 0.74-1.18. The frozen limit contains these choices through the N5-bias (15 %) and sX(1) = 0.247.
- This table is input for the UL_X17 column of Table S1. It is not a second limit.

----------------------------------------------------------------------------------------------------------------------
## 3. Figures and tables (draft captions)

**Fig. X1 `docs/latex/figs_all/x17d_X1_mee.pdf` (.png).**
- **Left.** M(e+e-) spectrum of the eta peak band (M4 0.5325-0.5575 GeV/c2), rows per MeV/c2, log scale.
  - Data N, S0, isBest >= 0: black.
  - Stacked: the frozen 2D model fitted with r = 0 to all 1050 cells:
    - CB from like-sign data (red);
    - eta-related terms C2+C3, G_eta (violet);
    - conversions eta -> pi+pi-gamma, C (blue);
    - direct eta -> pi+pi-e+e-, D (cyan);
    - other background n_j B(M4) (grey).
  - Dashed magenta: the same model plus the X17 line (17 MeV/c2) at the observed 90 % CLs limit (563 rows in the X2D range).
  - Lower pad: (N - mu_0)/sqrt(mu_0), with the X17 at the limit drawn as an expected pull. The yellow band is the M(ee)
    range that was blinded until the freeze.
- **Right.** The same peak band after the non-peak model terms (CB, other background, eta-related) are subtracted.
  - Compared with the fitted D and C (stacked), and with D + C + X17 at the limit (dashed).
  - Lower pad: the residual after D + C.
- L = 6079.2 nb^-1; raw rows (no 64/63).

**Fig. X2 `x17d_X2_limit.pdf`.**
- 90 % CLs upper limit on BR(eta -> pi+pi-X) BR(X -> e+e-) versus m_X, 10-40 MeV/c2 (38 masses).
- Observed:
  - normalisation (1), N5 fitted, PRIMARY: black line with points;
  - toy CLs at 16.7, 16.85 and 17.0 MeV/c2: red stars;
  - normalisation (2), a_eta = 0.38: blue dashed;
  - normalisation (3), peak a_P = 0.41: magenta;
  - (1) with e+ decay angle 1 + cos^2: red dash-dotted.
- Expected, FROZEN before unblinding: Asimov median (dashed) with +-1 sigma (green) and +-2 sigma (yellow); toy medians with
  16-84 % (open circles; 2000 toys at the quoted masses, 500 at 10/25/40).
- Grey dotted: Asimov median of the observed r = 0 fit (information).
- S0, isBest >= 0, X2D; L = 6079.2 nb^-1; 2D Poisson likelihood with the frozen nuisances.

**Fig. X3 `x17d_X3_p0.pdf`.**
- Top: local p0 versus m_X (asymptotic, q0 / qScale), normalisations (1)-(3). The three curves coincide.
- Bottom: fitted BR_X BR_ee (normalisation 1) with its Hessian error band.
- Largest local Z = 0.98 at 22 MeV/c2. No global p-value is needed (Z < 3).

**Fig. X4 `x17d_X4_valid.pdf`.** Validation at 17 MeV/c2:
- (a) CLs versus BR on the data: toys (1000 per point and hypothesis) against the asymptotic formula. Toy/asymptotic UL =
  0.983.
- (b) Feldman-Cousins construction: t_obs(r) against the 68 % and 90 % toy quantiles. 90 % CL interval [0, 3.45e-5].
- (c) Observed limit in the distribution of the 2000 frozen background-only toy limits (61 % of the toys lower).
- (d) Observed r_hat/sigma in the background-only toys: empirical p 0.295, asymptotic p0 0.311.
- (e) Observed limit recomputed with 200 CB and 100 C4 template replicas (RMS 0.19 and 0.28 x 1e-5).
- (f) Background n_j per M(ee) row: full fit divided by the frozen side-band fit. In 12-28 MeV/c2 this is measured
  divided by the frozen prediction.

**Tables** (`docs/latex/tables_all/`):

| file | content |
|---|---|
| `x17d_D1.txt`, `x17d_D1.tex` | **Table D1** at 16.7 / 16.85 / 17.0 MeV/c2 x normalisations (1)-(3): observed asymptotic and toy UL, frozen expected median [16 %, 84 %], r_hat +- error, local Z and p0, N_X^UL, MC RMS, 1 + cos^2 UL, provisional sigma x BR x BR; counting row; FC row |
| `x17d_scan.txt` | observed scan, 38 masses |
| `x17d_fit.txt` | full fit, GoF, n_j per row (fitted, predicted), counting |
| `x17d_checks.txt` | toy CLs, FC points, MC replicas, percentiles in the frozen toys |
| `x17d_variants.txt` | information only, Sect. 2.9 |

**Draft caption for Table D1.** "Observed 90 % CLs upper limits on BR(eta -> pi+pi-X17) BR(X17 -> e+e-) for three X17
masses and the three normalisations of the frozen procedure:
- (1) r x 2.68e-4 with the direct-decay yield N5 fitted (primary);
- (2) N_X/(a_eta s_X2D);
- (3) relative to the measured eta peak yield.

Expected = frozen background-only toys (median and 68 % band), computed before the M(ee) 12-28 MeV/c2 rows were opened. Data
S0, isBest >= 0, X2D, L = 6079.2 nb^-1, isotropic e+ decay angle. Template statistics (MC) and the 1 + cos^2 angular
hypothesis are quoted separately."

----------------------------------------------------------------------------------------------------------------------
## 4. Acceptance criteria and fall-backs

| item | value | verdict |
|---|---|---|
| precondition: x17_frozen.txt exists; frozen md5 sums unchanged | all equal | PASS |
| sealed side-band cells = hists_v1.root | max difference 0 | PASS |
| reproduction of the freeze (B_pred, Asimov median) | 6311.7 / 2.7589e-5 against 6311.7 / 2.7590e-5 | PASS |
| X-d: observed CLs (q~, r >= 0) | 38 masses x 3 normalisations, all fits converged | DONE |
| X-d: Feldman-Cousins cross-check | [0, 3.45e-5] at 17 MeV/c2 | DONE |
| X-d: scan 10-40 MeV/c2 (1 MeV; 0.5 MeV in 15-20) and 16.7, 16.85, 17.0 | 38 masses | DONE |
| X-d: local p0; global p if Z > 3 | Z_max 0.98, so no global p | DONE; rule not triggered |
| X-d: counting in the frozen window | +154 +- 1198 rows | DONE |
| frozen rule: toy vs asymptotic UL within 10 % | 0.958-1.041 | PASS (asymptotic primary) |
| frozen rule: CL_b = 0 | never | not triggered |
| fall-back: Pearson chi2/ndf of the full fit > 1.5 | 0.993 | not fired |
| fall-back: N5 error > 50 % | 13.7 % (full fit) | not fired; (1) stays primary |
| fall-back: slope off by > 5 % | 1.0385 (X-c injection; not re-tested on data) | not fired |
| simulation-limited (MC error > data stat. error) | r_hat MC/stat = 0.18 | not fired |
| observed limit inside the frozen +-1 sigma band | 38 of 38 masses | consistent with background only |

No fall-back fired, and no deviation from the frozen procedure.

----------------------------------------------------------------------------------------------------------------------
## 5. Caveats (flagged)

1. **Details that the frozen text left open.** Each is filled in the most direct way; none is a deviation.
   - **N(X17c) in the counting check.** I used the X2D overlap fractions, the same rule as the frozen B_pred, so that both
     terms use one estimator. The exact row count (6424) is given next to it: excess +112 instead of +154 rows.
   - **The r grid of the toy CLs.** It was not frozen. I used 9 points, 0.55-1.8 x the asymptotic UL.
   - **The FC construction.** Profile construction (pseudo-data from the data fit at fixed r) with t_r / qScale, the same
     toy generator as the frozen ToyCLs; 14 values of r, 2000 toys each.
   - **'MC'.** The RMS of the observed limit over 200 CB plus 100 C4 replicas, combined in quadrature.
2. **The masked-row prediction of X-c is off row by row** (-19 % to +10 %); its sum is right to 0.01 %.
   - The pattern follows the fluctuations of the simulated reference rows.
   - It does not enter the observed limit, where n_j is free per row.
   - It does enter the frozen expected limits and B_pred(X17c). The Asimov median built from the full data is 2.52e-5, not
     2.76e-5.
   - The frozen expectation is still what Fig. X2 and Table D1 show, as frozen.
3. **r_hat +- Hessian error is a fit result, not a measurement.**
   - The profile likelihood is about 40 % wider than the Hessian error at the limit (as X-c found with injection toys).
   - The interval to quote is the FC interval.
4. **Model dependence.** The limit assumes:
   - a prompt X17 -> e+e- decay;
   - an isotropic e+ angle;
   - three-body phase space for M(pi+pi-).

   With an e+ angle 1 + cos^2 the limit is 3.67e-5 (x 1.198), quoted separately.
5. **Template statistics (11 % of the UL, normalisation 1)** are quoted separately, as frozen (D13). The C4 replica mean
   is 6 % below nominal.
6. **The eta-scale inputs come from WP5, which is a cross-check (WP4 failed F-main).**
   - The a_eta prior and a_P are WP5 values; they are frozen and kept.
   - WP6 (side bands, the main result) gives Y_peak = 5088 instead of WP5's 5045. a_P would then be 0.418 instead of
     0.415, which changes UL (3) by about -1 %.
   - All WP6 eta-scale estimates agree with a_eta = 0.38.
7. **sigma x BR x BR is PROVISIONAL**: WP8 is missing; sigma(pp -> eta X) = 0.343 mb is built from WP5 and WP2 (Sect. 2.8).
   WP10 must replace it by WP8's value: multiply 3.06e-5 by sigma(pp -> eta X).
8. **Not used by the frozen procedure:** isBest == 1, ECAL, any OA cut. The limit is for isBest >= 0 only.
9. **Variants (Sect. 2.9)** are information. Three classes (0.74) and exp(pol3) (0.86) give lower limits; V2 and the
   background forms without the M(ee) dependence give limits up to 1.18 x higher. All are inside what the frozen N5-bias
   and sX cover.
10. **Machine.** Up to 6 parallel ROOT processes (26 jobs). The CPU package reached 79 C (limit 90 C). Nothing was deleted
    and no foreign file was overwritten. Small helper files of mine sit in /hdd3/.../x17d/: `ls_sealed.C`, `run_job.sh`,
    `jobs_x17d.txt`, `jobs_done.txt`.

----------------------------------------------------------------------------------------------------------------------
## 6. Exact commands (from /home/przygoda/HADES/FAT/FAT; ACLiC build dir /hdd3/przygoda/HADES/FAT_out/all_report/aclic_x17d)

    B='gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_x17d", kTRUE)'
    root -l -b -q -e "$B" 'plots/all_report/x17d_unblind.C+O'      # 25 s: reproduction, full fit, scan, p0, counting
                                                                     #   -> x17d/x17d_state.root, x17d_unblind.log,
                                                                     #      tables x17d_scan.txt, x17d_fit.txt
    # compile once, then 26 jobs, 6 in parallel (about 3 min in total):
    for m in x17d_toycls x17d_fc x17d_mc; do root -l -b -q -e "$B" -e ".L plots/all_report/$m.C+O"; done
    xargs -P 6 -L 1 sh /hdd3/przygoda/HADES/FAT_out/all_report/x17d/run_job.sh < /hdd3/przygoda/HADES/FAT_out/all_report/x17d/jobs_x17d.txt
    #   toycls m norm  -> x17d_toycls.C+O(m, norm, 1000)   (m = 16.7, 16.85, 17.0; norm 1-3)
    #   mc m           -> x17d_mc.C+O(m)                   (200 CB + 100 C4 replicas)
    #   fc r           -> x17d_fc.C+O(r, 2000)             (r = 0, 0.01, ..., 0.12, 0.14)
    root -l -b -q -e "$B" 'plots/all_report/x17d_variants.C+O'     # information table x17d_variants.txt
    root -l -b -q -e "$B" 'plots/all_report/x17d_figs.C+O'         # Figs X1-X4, x17d_D1.txt/.tex, x17d_checks.txt

**Code** (new, mine; `plots/all_report/`):
- `x17d_lib.h`: LoadFull (side bands + sealed rows), LocalP0, TStat/FCToys, Comps, GoodnessOfFit;
- `x17d_unblind.C`, `x17d_toycls.C`, `x17d_fc.C`, `x17d_mc.C`, `x17d_variants.C`, `x17d_figs.C`.

They include `x17c_lib.h` / `x17c_model.h` unchanged.

**Outputs.**
- `/hdd3/przygoda/HADES/FAT_out/all_report/x17d/`: x17d_state.root, x17d_unblind.log, toycls_m*_n*.txt,
  fc_m1700_n1_r*.txt, mc_m*.txt, log_*.txt.
- Figures: `docs/latex/figs_all/x17d_X{1_mee,2_limit,3_p0,4_valid}.{pdf,png}`.
- Tables: `docs/latex/tables_all/x17d_{D1.txt,D1.tex,scan.txt,fit.txt,checks.txt,variants.txt}`.

----------------------------------------------------------------------------------------------------------------------
## Verification (independent verifier, prefix vx, 2026-11-14 08:15-09:00)

**Verdict.** The unblinding followed x17_frozen.txt exactly, and the observed limit is confirmed:
3.06e-5 at 16.85 MeV/c2 (normalisation (1), 90 % CLs).
- An independent method that shares no code with the frozen model gives 2.7-3.4e-5.
- The conversion r -> BR and all normalisation constants were recomputed from the histograms; they are correct.
- The expected-band toys are sound.

I found three substantive errors in the TABLE D1 / Fig. X1 layer, plus two cosmetic ones. All are fixed in x17d_figs.C, and
its outputs were regenerated. None of them changes the BR limit itself. Some statements of Sect. 2.8, 3 and 5 above are
superseded (list in V.3).

### V.1 What I did, step by step (plain words)
1. **Frozen procedure against the X-d code.**
   - I re-read x17d_unblind.C, x17d_lib.h, x17d_toycls.C, x17d_fc.C, x17d_mc.C and x17d_figs.C against x17_frozen.txt
     Sect. 2-10.
   - After X-d, the md5 sums of the frozen code, the inputs and x17_frozen.txt itself (60a5376a) are unchanged.
2. **Data and blinding.** I made my own read-only dump of hists_v1.root and the sealed file (vx_dump.py).
   - The sealed side-band cells equal hists_v1.root (max |diff| 0).
   - The mask in hists_v1.root is intact (-1).
   - Totals: 255 331 / 197 752 / 57 579 rows, as in X-d.
3. **Normalisation from scratch** (vx_indep.py), from the S0__nom X2D templates with flows:
   - sigma_vis(C5) = 1.2464 nb, so N5ref = 7458.8 rows;
   - C6 per unit BR = 6257.5 nb, so s_X2D = 3.7446e7 rows and eps_X(S0, X2D) = 5.87e-3;
   - R_X = C6 x 2.68e-4 / C5 = 1.3455.

   All three equal the frozen numbers. The weight conventions agree:
   - The eta-set weight is 2 x BR per forced decay for C5, and the C6 weight is 2 (per unit BR).
   - C5 is at 1.0665 mb x 2.68e-4 = 285.8 nb.
   - WP0's check "X17 at 1e-5 = 0.07127 nb" passes.

   So the X17 rows a sigma_eta BR_X BR_ee eps_X L are equal to N5 r R_X. With r = BR_X BR_ee / BR_dir this gives BR = r x 2.68e-4.

   Other checks:
   - s1 = 3501.6 nb is WP0's exact C6_win(X17c).
   - a_P = 5045/12165 = 0.4147.
   - The line fractions X_j of x17_lines.root equal the C6 X2D fractions in every row (ratio 1.000).
   - The 1 + cos^2 variant is the x17_lines "C2" R_X.
4. **CLs arithmetic by hand** (vx_prof.C).
   - I scanned the profile likelihood F(r) on the data and on the Asimov set, with the frozen model.
   - From the scan I recomputed CLs with Cowan et al. eqs. 65-66.
   - UL: 3.0651e-5 against LimitAsym 3.0646e-5 (16.85 MeV/c2, (1)); 3.2861e-5 against 3.2861e-5 ((2)).
5. **Independent simple method** (vx_indep.py, numpy only).
   - (a) In every 2 MeV/c2 M(ee) row the eta-peak yield Y_j is fitted from the M4 spectrum: Poisson likelihood, CB fixed,
     the eta-related terms fixed at the frozen prior centre, exp(pol2) background free per row, C5 peak shape. The errors are
     multiplied by sqrt(qScale).
   - (b) The 30 yields are fitted with N5 [D_j + r R_X f_X X_j (1 + sX)^theta] + N4 C_j. Here f_X = 1.030 is the response of
     the extractor to the X17 M4 shape.
   - (c) Limits: a Bayesian one (flat prior in r >= 0; N5 and N4 marginalised analytically, theta marginalised), an
     asymptotic CLs one, and a pure counting variant (peak yields of the rows 14-24 MeV/c2 against N5 D + N4 C fitted
     outside 12-28 MeV/c2).
6. **Bayesian limit on the frozen likelihood** (vx_prof.C): flat prior in r >= 0 on exp(-(F(r) - F_min)/(2 qScale)).
7. **Direct limit on the number of X17 rows** (vx_nx.C). The frozen data, model and q~ are used, with POI = rows: the X17
   rows are r xA N5ref R_X, and the yield nuisance is off, since a row count needs no efficiency or BR factor.
8. **Frozen background-only toys** (vx_toys.py): all toy trees of x17c read again.
9. **Figures X1-X4 and Table D1.** I checked them number by number against the logs and the toy files.

### V.2 Numbers (S0, isBest >= 0, X2D, L = 6079.2 nb^-1, raw rows; 16.85 MeV/c2 unless stated)
**Observed limit and independent checks.**

| check | 90 % UL on BR_X BR_ee | comment |
|---|---|---|
| frozen, CLs asymptotic (1) | 3.065e-5 | toy CLs 3.07e-5 |
| CLs recomputed by hand from the F(r) scan (1) | 3.065e-5 | arithmetic confirmed |
| independent two-step method, exp(pol2) background, CLs (1) | 2.71e-5 | Asimov median 2.68e-5 (frozen Asimov 2.78e-5) |
| same, Bayesian with the X17 yield nuisance (1) | 3.10e-5 | 2.85e-5 without the yield nuisance |
| same, background pol1 / pol3 in the exponent, CLs | 2.85e-5 / 2.25e-5 | Bayesian 3.31e-5 / 2.72e-5 |
| same, normalisation (2), CLs / Bayesian | 3.07e-5 / 3.34e-5 | frozen (2): 3.29e-5 |
| Bayesian, flat prior, on the FROZEN profile likelihood (1) / (2) | 3.98e-5 / 3.89e-5 | +30 % / +18 % over CLs (caveat V.5.1) |
| Feldman-Cousins (X-d) | [0, 3.45e-5] at 17 MeV/c2 | |

- The fitted signal of the independent method is r_hat x 2.68e-4 = (+0.03 +- 1.26) x 1e-5. The frozen 2D fit gives
  (+0.54 +- 1.17) x 1e-5.
- Counting variant: rows 14-24 MeV/c2 hold an excess of +28 +- 167 peak rows, so N_X < 397 rows (statistics only).

**Number of X17 rows (Table D1, N_X).**
- Direct 90 % CLs limit: **N_X < 393 raw rows** (394 / 393 / 392 at 16.7 / 16.85 / 17.0 MeV/c2); 399 rows after x 64/63.
  - N_X_hat = +93 +- 192 rows.
  - The counting variant above gives the same, 397.
- The old Table D1 value, 563 rows, is r_UL x N5(r = 0 fit) x R_X: a conversion of the BR limit, not a limit on rows. The
  conditional fit at the UL has N5 = 2696 and theta_X = -0.46. So the model holds 375 X17 rows at the limit.

**Model nonlinearity (Asimov).** sigma_A(r) = r / sqrt(q_A) rises with r:

| r_hat x 2.68e-4 | 0.54e-5 | 2.7e-5 | 4.0e-5 |
|---|---|---|---|
| sigma_A (BR) | 1.33e-5 | 1.57e-5 | 1.80e-5 |

The cause is that normalisation (1) ties the X17 yield to N5, and N5 can trade against the line.

**Frozen expected toys.**
- At the quoted masses: 2000 toys, 2000 unique ids, 0 failed fits in all three normalisations.
- The quantiles reproduce x17_frozen.txt Sect. 9 to all printed digits.
- Mean number of rows per toy: 255 333 (data 255 331).
- The observed UL (1) sits at percentile 0.605 / 0.606 / 0.614.
- 38 of 38 observed limits are inside the frozen +-1 sigma Asimov band; ratio to the frozen median 0.742-1.384.

**Provisional sigma x BR x BR, corrected** (spec AMENDMENT 1, 08:20):
- sigma(pp -> eta X) = a_eta x 1.0665 mb x C_sec (1 - f_eta) = 0.3823 x 1.0665 x 1.011 x 0.97 = 0.3998 mb.
- This gives sigma x BR_X x BR_ee < **12.25 nb** (was 10.52 nb with (1 - f_win)).
- Information: the route through the rows limit, 392 / 3.745e7 x 1.0665 mb x 0.981, gives 10.9 nb.

### V.3 Errors found, changes made, superseded statements
**Changes** (x17d_figs.C only; md5 e7be8b73 -> 804259d1; every change is marked "[vx verification]" in the code).
- The file was rerun with ACLiC build dir aclic_vx.
- Regenerated: x17d_X1..X4 (.pdf, .png), x17d_D1.txt and x17d_D1.tex. x17d_checks.txt is regenerated with identical content.
- Originals are kept in /hdd3/przygoda/HADES/FAT_out/all_report/vx/backup_x17d_pre_vx/.
- No frozen file, no x17c file, x17d_unblind.C or another macro's output, and no BR number was changed.

1. **sigma x BR x BR used (1 - f_win) = 0.842.** Amendment 1 makes C_sec (1 - f_eta) = 0.981 binding for every eta yield.
   kCsecFwin = 1.0110 x 0.97, so 10.52 -> 12.25 nb. (The amendment came after X-d ran; this was not an error of X-d.)
2. **Table D1 column MC_RMS was the RMS of the 300 pooled replica limits.** CB and C4 are independent sources, so the right
   combination is the quadrature sum (as in x17d_checks.txt and Sect. 2.6 above). At 16.85 MeV/c2:

   | normalisation | old (pooled) | new (quadrature) |
   |---|---|---|
   | (1) | 2.39e-6 | 3.34e-6 (10.9 %) |
   | (2) | 1.69e-6 | 2.47e-6 |
   | (3) | 1.35e-6 | 1.98e-6 |
3. **Table D1 "N_X_UL" was not a limit on N_X** (spec WP9 asks for "limits on N_X").
   - N_X_UL is now the direct 90 % CLs limit: 394 / 393 / 392 rows, the same in all normalisations.
   - The old conversion stays as a new last column, rUL_N5_RX (562 / 470 / 427).
   - The x17d_D1.tex column N_X^UL now holds the direct limit.
   - Fig. X1 draws the X17 line with 392 rows (was 563); legends "at UL: N_X = 392" and "(392 rows)".
4. **Cosmetic.**
   - Fig. X1 right: the exponent of the UL label was clipped.
   - Fig. X4(f): the legend markers sat at y 1.22-1.27 inside the 12-28 MeV/c2 band and looked like data; the legend now sits
     in the empty lower right.

**Superseded statements in this note.**
- Headline/Sect. 2.8 "N_X at the UL 563 / 469 / 427 rows" and "563 raw rows correspond to 572": now **N_X < 392 raw rows (398
  after 64/63) at 17 MeV/c2**, for every normalisation.
- Sect. 2.8 "sigma x BR x BR < 10.5 nb (sigma(pp -> eta X) = 0.343 mb)": now **12.2 nb (0.400 mb)**. The uncertainty
  "C_sec (1 - f_win) +7/-6 %" is now +3/-4 %. The "information" route gives 10.9 nb, not 13.5 nb.
- Its explanation "the fitted N5/N5ref = 0.46": the factor actually used was N5(r = 0 fit)/N5ref = 3659/7459 = 0.49.
- Sect. 3, caption of Fig. X1: "(563 rows in the X2D range)" now reads "(392 rows: the 90 % CLs limit on the number of X17
  rows)".
- Sect. 3, Table D1 caption: the N_X column is the direct limit on rows. MC is CB (+) C4 in quadrature.

### V.4 Acceptance (verifier) and fall-backs
| item | value | verdict |
|---|---|---|
| frozen procedure applied unchanged | code, configuration, grid, nuisances, qScale rule, toy CLs, FC, p0 and counting all as frozen; md5 unchanged | PASS |
| blinding and sealed data | identical side bands; totals reproduced independently | PASS |
| observed UL, order of magnitude, independent method | 2.7-3.4e-5 against 3.06e-5 | PASS |
| CLs arithmetic | 3.0651e-5 by hand against 3.0646e-5 | PASS |
| r -> BR, R_X, N5ref, s_X2D, s1, a_P, BR_dir, units | recomputed from histograms, identical | PASS |
| expected toys and bands | exact reproduction, 0 failures, 38/38 masses inside +-1 sigma | PASS |
| Table D1 / Fig. X1 normalisation layer | 3 errors | FIXED |
| figures | 2 cosmetic defects | FIXED |

- The spec's Table D1 item (2) says N_X / (a_eta s1) with the window s1. The frozen procedure divides N_X over the whole X2D
  range by s_X2D, the matching efficiency. This is the correct adaptation, decided at the freeze; it is not an X-d
  deviation. s1 is used only in the counting check.
- No fall-back fired. I re-checked the frozen ones: chi2/ndf 0.994; N5 error 13.7 %; toy/asymptotic 0.96-1.04; Z_max 0.98.

### V.5 Caveats (flagged for WP10)
1. **Method dependence.** In normalisation (1) the likelihood in r is not Gaussian: the X17 yield is r N5, and N5 can drop
   (2696 at the UL against 3659 at r = 0).
   - A flat-prior Bayesian limit on the same likelihood is 3.98e-5 (+30 %).
   - The FC upper edge is 3.45e-5.
   - CLs is the frozen and spec method and is confirmed by toys, so 3.06e-5 stands. The report should say that other
     standard constructions give up to 4.0e-5.
2. **C4 template statistics.** The mean limit over the C4 replicas is 6 % below the nominal one. Replicas carry the nominal's
   noise twice, so the nominal limit itself is probably about 6 % optimistic, i.e. about 3.25e-5 with a template of
   infinite statistics. This is inside the separately quoted MC error of 11 %.
3. **The simulation expectation (1.9e-5 / 5.7e-6) is not like-for-like.** It assumed the SMASH eta yield (a_eta = 1) and the
   counting window. With a_eta = 0.38 the counting expectation would be about 2.6 times higher. Its caption must say this
   (comparison only).
4. **The independent method moves with the background form.** Its r_hat x 2.68e-4 is
   - -1.80e-5 +- 1.21e-5 with pol3 in the exponent;
   - about 0 with pol1 or pol2.

   The frozen 2D fit gives +0.55e-5. All of these are consistent with no signal; the spread is a background-model effect
   of about 1 sigma.

### V.6 Files and exact commands (from /home/przygoda/HADES/FAT/FAT)
**New files of mine.**
- Code: plots/all_report/vx_dump.py, vx_indep.py, vx_toys.py, vx_prof.C, vx_nx.C.
- Outputs in /hdd3/przygoda/HADES/FAT_out/all_report/vx/:
  - vx_inputs.npz, vx_indep.txt, vx_toys.txt, vx_prof_m1685_n1.txt, vx_prof_m1700_n1.txt, vx_prof_m1685_n2.txt, vx_nx.txt;
  - logs log_prof_*.txt, x17d_figs_rerun_stdout.txt;
  - backup_x17d_pre_vx/.

**Commands.**

    source /usr/local/root/bin/thisroot.sh; export PYTHONPATH=/usr/local/root/lib
    python3 plots/all_report/vx_dump.py; python3 plots/all_report/vx_indep.py; python3 plots/all_report/vx_toys.py
    B='gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_vx", kTRUE)'
    root -l -b -q -e "$B" 'plots/all_report/vx_prof.C+O(16.85, 1)'   # also (17.0, 1), (16.85, 2)
    root -l -b -q -e "$B" 'plots/all_report/vx_nx.C+O'
    root -l -b -q -e "$B" 'plots/all_report/x17d_figs.C+O'           # regenerates Figs X1-X4, Table D1 (fixed)
