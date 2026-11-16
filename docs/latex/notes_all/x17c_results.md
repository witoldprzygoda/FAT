# WP9 X-c (X17 freeze) -- results note (prefix x17c)

Date 2026-11-14. Branch pp45_pippimepem (worktree /home/przygoda/HADES/FAT/FAT). Nothing committed, no git state changed.
ScanEngine.h, scan_*.C, the files of WP0-WP5 and x17ab (x17_model.h included) untouched; FAT_simwt not used.

**Main deliverable:** `/hdd3/przygoda/HADES/FAT_out/all_report/x17/x17_frozen.txt` (md5 60a5376aeabda9aacabd6584f0f72572).
It holds the complete frozen procedure and the expected limits. The frozen numbers are also in machine-readable form in
`/hdd3/przygoda/HADES/FAT_out/all_report/x17c/x17c_config.txt`.

**Blinding (D19).**
- No data cell with M(ee) 12-28 MeV/c2 was read. In hists_v1.root these cells hold -1; the code checks only that value
  and never copies it.
- The sealed file `sealed/hists_v1_x17_sealed.root` was never opened.
- The data tree was read once (event-group covariance, `x17c_disp.C`). The very first operation on every row drops it if
  12 <= M(ee) < 28 MeV/c2.
- Already known before this package (D19): the X17c window count 6424 (CB 209.5) and the X17P count 5258. This note
  freezes a background prediction for X17c but does not compare it with the count.

----------------------------------------------------------------------------------------------------------------------
## 1. What was done, step by step (plain words)

1. **Event-group covariance** (`x17c_disp.C`).
   - One event gives several rows, so neighbouring cells are correlated.
   - For the 770 side-band cells of X2D (M(ee) 2-12 and 28-62 MeV/c2, S0) I summed, over event groups, the products of the
     row counts of every pair of cells.
   - Check: the counts and the diagonal equal hists_v1.root exactly (max difference 0 and 2e-13).
2. **Model.** The model is the X-a model (x17ab), copied into `x17c_model.h` with four additions:
   - **Normalisation switch.**
     - (1) X17 rows = r N5 R_X, with N5 fitted.
     - (2)/(3) X17 rows = r xA N5ref R_X, with a fixed eta scale xA: a_eta of WP5, or the peak-based a_P.
     - In (2)/(3) the eta scale of the X17 term is separate from the a_eta parameter of the eta-related background terms.
       Otherwise the side-band fit would move the X17 normalisation; in the 3-class variant it pulled a_eta to 0.09.
   - **X17 yield nuisance.** It is log-normal, (1 + sX)^thX.
   - **Line shape in the likelihood.** Two nuisances of the X17 line: width +-10 % and M(ee) scale +-2 % (spec WP9),
     each implemented as a half-difference template.
   - **Variance factor qScale.** The test statistic is divided by qScale. The toys are drawn over-dispersed, as
     qScale x Poisson(mu/qScale).

   Everything else is x17ab's validated model:
   - n_j free in every M(ee) row (no simulated background per row, D18);
   - background shape x (1 + (b + b'u) t + (c + c'u) t^2) per class (bkgMode 4);
   - k_E1 = 1 +- 0.5 and k_E2 = 1 +- 0.05;
   - the conversion-tail nuisance;
   - q~_r and asymptotic CLs.

   As the spec asks, the background M4 shapes B_E1, B_E2 are the WP5 shapes (`x17bkg` of wp5/results.txt).
3. **Inputs.**
   - X2D of hists_v1.root, selection S0 (isBest >= 0, vz > -160, cor_ok, both pions in the vertex fit), data CB
     (like-sign, spread rule), simulation S0__nom (V3) in raw rows (nb x 5984.2).
   - X17 line and R_X(m) from x17_lines.root (x17ab).
   - **a_eta prior: 0.382 +- 0.165.** This is the WP5 F-main total error x sqrt(chi2/ndf) = 1.61, plus (in quadrature)
     the shift 0.149 of the WP5 quadratic variant, which its chi2 does not reject.
   - **g prior: 1.72 +- 0.96**, built the same way.
   - **a_P = Y_peak(WP5)/Y_peak(sim) = 5045/12165 = 0.415 +- 0.109** (the same error rule).
4. **Side-band fit.** The model was fitted with r = 0 to the 770 side-band cells only; the masked rows are excluded from
   the likelihood. Goodness of fit: Pearson chi2 per class, per row, and the residual near m_eta.
5. **Variants of the side-band fit.** 22 variants:
   - background shape: simulated exp(pol1/2/3), kappa_B, flat pol2, three classes, tilt only, no M(ee) dependence;
   - CB: factorised, geometric mean;
   - k_E1 free;
   - masked-row prediction: pol1, pol3, like-sign, flat and smoothed references;
   - peak/direct templates: V2, GS, r = 0.89/1.10, s +- 0.0009.

   They give three things:
   - the **N5-bias** (half spread of N5 over the background-shape and CB variants). Rule fixed before the fits: drop a
     variant whose deviance is worse than nominal by more than 2 sqrt(2 ndf). None was dropped;
   - the direct-shape, peak-width and mass-scale terms of sX (change of N5_hat R_X, resp. N5ref R_X);
   - the effect of each choice on the expected limit.
6. **Variance factor.**
   - I computed phi = w'Sw / w'diag(N)w for X17-like weights, w = (X17 template)/(model), with lines at 8-44 MeV/c2 and
     only side-band cells.
   - Frozen: qScale = the largest phi = 1.0299. The Pearson fall-back of the spec would override it if chi2/ndf > 1.5.
7. **Masked-row prediction.**
   - The background of the 8 masked rows was predicted as n_j = sim_j x exp(pol2(y)), with sim_j the raw simulated
     G_D + G_C row sum. The pol2 in ln(n_j/sim_j) was fitted to the 22 side-band rows. Weights: the fit error of n_j
     combined with the simulation error of sim_j.
   - Two global smoothings of the simulated reference were tried first and rejected:
     - a single exp(pol4) over all rows (chi2/ndf 27);
     - G_D and G_C separately with exp(pol3); its G_C extrapolation diverged.
     Neither follows the M(ee) shape. A local smoothing is kept as a variant.
8. **Background-only Asimov set** = the model with the side-band parameters and the predicted rows, r = 0, on all 1050
   cells.
9. **Asimov expected limits** (asymptotic CLs, Asimov bands): 38 masses (10-40 MeV/c2 in 1 MeV/c2 steps, 0.5 MeV/c2 in
   15-20, plus 16.7 and 16.85), three normalisations.
10. **Background-only toys** (`x17c_toys.C`).
    - Pseudo-data from the Asimov expectation, over-dispersed by qScale; constraint centres drawn.
    - Each toy fitted on all 1050 cells exactly as X-d will fit the data; its asymptotic CLs limit recorded.
    - 2000 toys at each of 16.7, 16.85 and 17.0 MeV/c2; 500 at 10, 25 and 40.
11. **Injection toys.** X17 at 17 MeV/c2 with BR = 1, 3, 6 x 1e-5, 500 toys each: linearity, pulls, coverage.
12. **Toy CLs against asymptotic CLs** (`x17c_toycls.C`). On the Asimov set at 17 MeV/c2: 8 values of r, 3000 toys per
    value and hypothesis.
13. **Freeze.** `x17c_frozen.py` writes x17_frozen.txt with:
    - md5 sums of code and inputs;
    - selection, layout, model, priors, test statistic, CLs settings, mass grid, quoted masses, counting window;
    - the expected limits;
    - the step list for X-d.

----------------------------------------------------------------------------------------------------------------------
## 2. Numbers (S0, isBest >= 0, L = 6079.2 nb^-1, raw rows; BR = BR(eta -> pi+pi-X17) BR(X17 -> e+e-))

### 2.1 Side-band fit (r = 0), table `x17c_sidebands.txt`
- Data in the side bands: 197 752 rows in 770 cells.
- Pearson chi2/ndf = 749.0/730 = **1.026**.
  - By class: E1 99.1 for 140 cells, E2 650.0 for 630 cells.
  - Largest row: 49.3/35 (30-32 MeV/c2).
- Fitted parameters:
  - N5 = 3380 +- 686 (**20.3 %**), N4 = 1803 +- 429;
  - k_E1 = 1.04 +- 0.50, k_E2 = 1.000 +- 0.050;
  - a_eta = 0.265 +- 0.142, g = 1.40 +- 0.93.
- Residual N - model at M4 0.5425-0.5525: E1 -41 +- 119 rows, E2 -76 +- 76 rows (no peak-shape misfit near m_eta).
- Masked-row prediction:
  - ratio n_j/sim_j 1.41 (12-14 MeV/c2) to 1.27 (26-28 MeV/c2);
  - 57 572 predicted rows in the 8 masked rows; Asimov total 255 326 rows.

### 2.2 Variance factor
- phi for X17-like weights: 1.017-1.030 (lines at 8-44 MeV/c2). Maximum 1.0299, which is frozen as qScale.
- Within-cell dispersion 1.002; all side-band cells together (w = 1) 1.078.
- Effect on the expected median: x 1.021.

### 2.3 Nuisance sizes (frozen; `x17c_variants.txt`)
- **sX(1) = 0.247.** Components (in quadrature):

  | component | value |
  |---|---|
  | R_X with direct-like M(pi+pi-) | 0.021 |
  | R_X with the generator pion angle | 0.041 |
  | slice/C6 validation | 0.043 |
  | BR_dir | 0.041 |
  | **N5-bias** (N5 range 2993-4026 / 3380) | **0.153** |
  | **direct shape (V2)** | **0.158** |
  | peak width | 0.078 |
  | mass scale | 0.029 |

- **sX(2) = 0.436.** Dominated by the a_eta prior (0.431); the direct shape is 0, as expected (C6 per BR does not depend
  on the direct-decay weights).
- **sX(3) = 0.270**, with a_P at 0.262.
- **Conversion tail:** sC = 0.194 (RMS of the 100 C4 replicas).

### 2.4 Expected 90 % CLs limits (`x17c_expected.txt`, `x17c_expected_asimov.txt`, `x17c_expected.tex`)
Normalisation (1) is PRIMARY: BR = r x 2.68e-4 with N5 fitted.

| m [MeV/c2] | method | -2 sigma | -1 sigma | median | +1 sigma | +2 sigma | (2) median | (3) median |
|---|---|---|---|---|---|---|---|---|
| 16.70 | Asimov | 1.19e-5 | 1.74e-5 | 2.81e-5 | 5.09e-5 | 1.05e-4 | 2.68e-5 | 2.27e-5 |
| 16.70 | 2000 toys | 0.90e-5 | 1.52e-5 | **2.64e-5** | 4.84e-5 | 9.30e-5 | 2.64e-5 | 2.24e-5 |
| 16.85 | Asimov | 1.18e-5 | 1.72e-5 | 2.78e-5 | 5.06e-5 | 1.04e-4 | 2.66e-5 | 2.25e-5 |
| 16.85 | 2000 toys | 0.90e-5 | 1.51e-5 | **2.62e-5** | 4.82e-5 | 9.29e-5 | 2.61e-5 | 2.23e-5 |
| 17.00 | Asimov | 1.17e-5 | 1.70e-5 | 2.76e-5 | 5.01e-5 | 1.03e-4 | 2.63e-5 | 2.23e-5 |
| 17.00 | 2000 toys | 0.89e-5 | 1.50e-5 | **2.60e-5** | 4.77e-5 | 9.24e-5 | 2.59e-5 | 2.20e-5 |
| 10 | Asimov / 500 toys | | | 5.24e-5 / 4.91e-5 | | | | |
| 25 | Asimov / 500 toys | | | 1.80e-5 / 1.77e-5 | | | | |
| 40 | Asimov / 500 toys | | | 1.03e-5 / 1.01e-5 | | | | |

- **Frozen choice for the report:** at the quoted masses the toy median and toy bands; versus m_X the Asimov curve with
  the toy points overlaid.
- Toy/Asimov median: 0.94 at 16.7-17.0 MeV/c2, 0.94 at 10, 0.98 at 25 and 40. The toy -2 sigma band is lower
  (0.89e-5 against 1.17e-5).
- At the 17 MeV/c2 Asimov median, the limit corresponds to 457 X17 rows.
- X17 rows per unit BR at a_eta = 1: s_X2D = 3.745e7.
- **Comparison with the pre-data projection of x17ab:** 2.83e-5 at 17 MeV/c2 (simulation Asimov, a_eta 0.40, background
  1.35 x SMASH). The data-driven value 2.76e-5 is within 3 % of it.

### 2.5 Checks of the frozen fit
**Background-only toys.**
- r_hat/sigma: mean -0.14, RMS 0.94 (17 MeV/c2); -0.04 to -0.10 at the other masses.
- Converged: 100 % (2000/2000 and 500/500).

**Injection at 17 MeV/c2, normalisation (1).** Mean r_hat/injected:

| injected BR | mean r_hat / injected | coverage of the 90 % CLs UL |
|---|---|---|
| 1e-5 | 1.106 | 0.990 |
| 3e-5 | 1.059 | 0.888 |
| 6e-5 | 1.061 | 0.876 +- 0.015 |

- Slope after the background-only offset (+1.0e-6): **1.0385**.
- Pull RMS 0.97 / 1.19 / 1.22, mean -0.22 to -0.35: the Hessian error of r_hat is too small at large signal. The limit
  itself comes from q~, not from that error.

**Injection, normalisations (2) and (3).**
- Slopes 1.245 and 1.106. By construction they should be N5/(a N5ref) = 1.185 and 1.093, because the toys inject
  normalisation-(1) rows.
- The ratios are therefore 1.050 and 1.012.
- The remaining +5 % of (2) is the mean of a log-normal yield factor with randomised global observable:
  exp(ln^2(1.436)/2) = 1.067. The medians agree: 0.975.
- Coverage 0.94-0.99.

**Toy CLs against asymptotic CLs** (17 MeV/c2, Asimov):

| r | 0.06 | 0.08 | 0.09 | 0.10 | 0.11 | 0.12 | 0.14 | 0.17 |
|---|---|---|---|---|---|---|---|---|
| CLs toys | 0.290 | 0.189 | 0.154 | 0.108 | 0.077 | 0.075 | 0.043 | 0.024 |
| CLs asymptotic | 0.286 | 0.178 | 0.141 | 0.112 | 0.089 | 0.071 | 0.046 | 0.025 |

- The toy errors are 0.004-0.014.
- Toy UL 2.82e-5 against asymptotic 2.76e-5: ratio **1.023**.

**Template statistics (MC; quoted separately, D13).** On the Asimov set at 17 MeV/c2:
- sigma(r_hat) = 1.16e-5;
- r_hat RMS over 50 CB replicas 1.1e-6, over 50 C4 replicas 1.4e-6.

### 2.6 Effect of single choices on the expected median at 17 MeV/c2 (nominal 2.76e-5; ratio)

| choice | ratio |
|---|---|
| without qScale | 0.979 |
| without sX | 0.953 |
| without the line-shape nuisances | 0.991 |
| a_eta and g fixed | 0.966 |
| without the conversion tail | 0.997 |
| k_E2 fixed | 1.000 |
| line width x 1.10 / x 0.90 (fixed) | 1.097 / 0.907 |
| R_X with direct-like M(pi+pi-) | 1.020 |
| **e+ decay angle 1 + cos^2 (quoted separately)** | **1.198** |

**Variants with their own side-band fit:**

| variant | ratio (1) | ratio (2) |
|---|---|---|
| simulated background exp(pol2) / exp(pol1) | 1.003 / 1.001 | |
| exp(pol3) (deviance +41) | 0.784 | |
| three classes | 0.745 | |
| flat pol2 | 1.018 | 0.933 |
| tilt only | 1.011 | 0.933 |
| no M(ee) dependence | 1.038 | 0.930 |
| CB factorised / geometric mean | 0.982 / 0.975 | |
| k_E1 free | 1.012 | |
| fill pol1 / pol3 | 1.004 / 1.002 | |
| fill like-sign / flat / local smoothing | 1.122 / 1.104 / 1.046 | |
| V2 | 1.186 | |
| r = 0.89 / 1.10 | 1.015 / 0.989 | |
| s up / down | 1.044 / 0.961 | |

### 2.7 Counting window X17c (frozen prediction, not compared with the count)
- Window: M4 0.5312-0.5566 GeV/c2 x M(ee) 15.5-24.5 MeV/c2.
- **B_pred = 6311.7 rows** (CB part 212.3). It is the background-only expectation of the frozen model, integrated with
  the overlap fractions of the X2D cells.
- Variants: fill pol1 6403, pol3 6372, local smoothing 7507; like-sign 9622 and flat 9043 (information only).
- **Frozen sigma_B = 1195.5 rows** (largest deviation among pol1, pol3, local smoothing).
- s1 (C6 per unit BR in X17c) = 2.0954e7 raw rows at a_eta = 1.

----------------------------------------------------------------------------------------------------------------------
## 3. Figures and tables (draft captions)

**`docs/latex/figs_all/x17c_sb_mee.pdf`.**
- Data N per MeV/c2 (black) in the side bands M(ee) 2-12 and 28-62 MeV/c2, S0, isBest >= 0, M4 0.5225-0.610 GeV/c2,
  L = 6079.2 nb^-1.
- Stacked: the model fitted with r = 0 to the side bands only:
  - CB from like-sign data (red);
  - eta-related terms C2 + C3 + G_eta (violet);
  - eta peak D + C (cyan);
  - other background n_j B(M4) (grey).
- Hatched: 12-28 MeV/c2, blinded (D19). There the dashed line is the PREDICTED background, not data.
- Lower panel: data/model. It is 1 by construction, because n_j is free per row.

**`x17c_sb_m4.pdf`.**
- M4 projections of the two side bands (left 2-12, right 28-62 MeV/c2) with the same fitted components.
- Lower panels: (N - mu)/sqrt(mu).
- Pearson chi2/ndf of the side-band fit 749/730.

**`x17c_fill.pdf`.**
- Fitted n_j divided by the simulated G_D + G_C row sum (black; errors: fit (+) simulation).
- Prediction for the masked rows, a pol2 in ln of the ratio (open blue).
- Grey crosses: the same divided by the locally smoothed simulation (variant).

**`x17c_expected.pdf`.**
- EXPECTED 90 % CLs upper limit on BR(eta -> pi+pi-X) BR(X -> e+e-) versus m_X.
- Background-only Asimov set of the data side-band fit (no data of M(ee) 12-28 MeV/c2 used); S0, isBest >= 0, X2D,
  L = 6079.2 nb^-1, 2D likelihood, q~_r.
- Curves:
  - normalisation (1): Asimov median (dashed), +-1 sigma (green), +-2 sigma (yellow);
  - points: toy median and 16-84 % (2000 toys at 16.7/16.85/17.0, 500 at 10/25/40);
  - (2) a_eta normalisation (blue dashed);
  - (3) peak normalisation (magenta dotted);
  - pre-data projection of x17ab (grey dash-dotted).

**`x17c_toys_ul.pdf`.**
- Distribution of the 90 % CLs upper limit of 2000 background-only toys at m_X = 17 MeV/c2, normalisations (1) (black),
  (2) (blue) and (3) (magenta).
- Vertical lines: Asimov median (red) and +-1/2 sigma (orange) of (1).

**`x17c_inject.pdf`.**
- Mean fitted against injected BR at 17 MeV/c2. Truth = side-band fit (r = 0) + X17 injected in normalisation (1);
  500 toys per point.
- (2) and (3) read the same rows with their own eta scale. Their expected slopes are 1.185 and 1.093.

**Tables** (`docs/latex/tables_all/`):

| file | content |
|---|---|
| `x17c_sidebands.txt` | side-band fit, n_j per row (fitted / predicted), variance factors |
| `x17c_variants.txt` | 22 variants: chi2, N5, expected medians, B(X17c); nuisance composition; single choices |
| `x17c_expected_asimov.txt` | Asimov bands at 38 masses, 3 normalisations |
| `x17c_expected.txt` | toy quantiles against Asimov at 6 masses, pulls |
| `x17c_inject.txt` | injection tests |
| `x17c_expected.tex`, `x17c_nuisances.tex` | LaTeX snippets: Table "expected limits", Table "frozen priors" |

**Logs:** `/hdd3/przygoda/HADES/FAT_out/all_report/x17c/` (`x17c_freeze.log`, `x17c_toycls_m17.txt`, `x17c_disp.txt`,
`x17c_summary.log`, `toys_*.root`, `x17c_state.root`).

----------------------------------------------------------------------------------------------------------------------
## 4. Acceptance criteria and fall-backs

| criterion | value | verdict |
|---|---|---|
| X-c: fit r = 0 on the side bands | done; chi2/ndf 1.026; converged | PASS |
| X-c: Asimov limit | done, 38 masses x 3 normalisations | PASS |
| X-c: 2000 background-only toys (median, +-1/2 sigma) | 2000 at each quoted mass, 100 % converged | PASS |
| X-c: x17_frozen.txt | written, with the md5 of code and inputs | PASS |
| D19: no data of 12-28 MeV/c2 looked at | masked cells never read; sealed file not opened | PASS |
| X17 fall-back: Pearson chi2/ndf > 1.5 | 1.026 | not fired |
| X17 fall-back: slope off by > 5 % | 1.0385 (norm. 1) | not fired |
| X17 fall-back: N5 error > 50 % | 20.3 % | not fired; (1) stays primary |
| spec X-b rule: toy and asymptotic CLs within 10 % | 1.023 | PASS; asymptotic stays primary |
| coverage >= 0.88 (X-b criterion, repeated on the frozen model), norm. (1) | 0.990 / 0.888 / 0.876 +- 0.015 at 1 / 3 / 6 x 1e-5 | PASS / PASS / marginal (0.3 sigma below 0.88) |
| background-only pull | mean -0.14, RMS 0.94 | OK (slightly conservative) |

No fall-back fired.

----------------------------------------------------------------------------------------------------------------------
## 5. Caveats and deviations (flagged)

1. **Side bands.** They are 2-12 and 28-62 MeV/c2, not the spec's 2-13 and 27-62. With 2 MeV/c2 bins from 2 MeV/c2, the
   bins 12-14 and 26-28 overlap 13-27 and are masked in hists_v1.root.
2. **Additions to the x17ab model.** Each is flagged in x17c_model.h:
   - normalisations (2)/(3) decoupled from the fitted a_eta;
   - log-normal yield nuisance;
   - line width and M(ee) scale as likelihood nuisances;
   - the variance factor.
   Reason for the variance factor: the X2D rows are correlated through events (multiple combinations). For X17-like
   combinations phi is only 1.02-1.03, but it is frozen and stated.
3. **The spec's "B = WP5 shapes" is followed** (x17ab had used the simulated exp(pol2)). On the side bands both give the
   same chi2 (deviance difference 0.8) and expected limit (+0.3 %).
4. **V2 counts in the direct-shape family** (spec S7/WP9), so sX(1) = 0.247 rather than about 0.19. V2 has no Kroll-Wada
   angular terms, so this is conservative. Effect on the limit: about 1-2 %.
5. **N5-bias = 15 %.** The two ends are:
   - exp(pol3), deviance +41 at the same number of parameters; kept by the rule fixed in advance; conservative;
   - "no M(ee) dependence".
   The 3-class and exp(pol3) variants also move the expected limit by about -22 %. The masked rows will tell which
   background form holds.
6. **Masked-row prediction** (it affects only the expectation, not the observed 2D limit, where n_j is free).
   - Reference choices change the expected median by -0 to +12 %.
   - They change the counting prediction B_pred(X17c) much more: 6312 nominal against 7507 for the local smoothing
     (9622 with the like-sign reference).
   - The counting cross-check is therefore weak; frozen sigma_B = 1196 rows.
   - Two global smoothings of the simulated reference failed (Sect. 1, step 7) and are not used.
7. **Toy against Asimov band.** The toy expected median is 6 % below the Asimov one at 16.7-17 MeV/c2, and the toy -2 sigma
   band is about 25 % lower. Frozen: toys at the quoted masses, Asimov versus mass.
8. **Large injected signal.** At 6e-5 the coverage is 0.876 +- 0.015 (criterion >= 0.88) and the r_hat pull RMS is 1.2. The
   Hessian error of r_hat under-estimates the error when N5 and r are both free. The CLs limit uses q~ and is not affected
   at BR <= 3e-5, where the expected limit lies. X-d should not quote r_hat +- Hessian error as a measurement without a
   toy check.
9. **M(ee) scale and resolution priors (+-2 %, +-10 %)** are the spec values. The data check on the E1 conversion peak and
   the pi0 Dalitz shape was not possible: the WP0 layouts have 2 MeV/c2 M(ee) bins, so the conversion peak sits in one
   bin, and a fine-binned simulation pass was out of scope. It is an open item; the priors are frozen as they are. Fixed
   width x 1.10 would change the expected median by +9.7 %; as a profiled nuisance the effect is 1 %.
10. **Normalisation (2)** carries a 43 % a_eta prior. It includes the shift of the WP5 quadratic variant, because WP5 is a
    cross-check (WP4 failed). Its expected median (2.63e-5) is close to (1) only because the side-band N5 (3380) is 1.19
    x a_eta N5ref.
11. **Not from this package:**
    - the a_eta prior and a_P come from WP5 (cross-check values); WP6 (side bands) was not yet available;
    - the frozen file allows X-d to state which WP6 numbers would change (2)/(3) but keeps the frozen values.
12. **Machine.** The first toy batch failed because six ACLiC processes recompiled the same library concurrently. It was
    rerun after a single compile. Stale logs `toys_log_*.txt` and test files `toys_m1700_i0000_9000.root` and `_9100.root`
    remain in `/hdd3/.../x17c/` (no deletion allowed); the summary reads only the chunks 0/500/1000/1500. The CPU package
    reached 89 C (limit 90) during the toys.
13. **Header timestamp.** `x17c_lib.h` was back-dated with touch (07:41:14) so that running toy jobs would not recompile.
    The only change after the toy compile was the added WindowX17c function. x17c_config.txt was identical before and
    after (checked with diff).

----------------------------------------------------------------------------------------------------------------------
## 6. Exact commands (from /home/przygoda/HADES/FAT/FAT; ACLiC build dir /hdd3/przygoda/HADES/FAT_out/all_report/aclic_x17c)

    B='gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_x17c", kTRUE)'
    root -l -b -q -e "$B" 'plots/all_report/x17c_disp.C+O'       # 10 s, 3 GB: x17c/x17c_egcov.root (side-band cells only)
    root -l -b -q -e "$B" 'plots/all_report/x17c_freeze.C++O'    # 40 s: x17c_config.txt, x17c_state.root, tables, sb/fill figures
    root -l -b -q -e "$B" 'plots/all_report/x17c_toys.C+O(17.0, 0., 9100, 2)'   # compile ONCE before the parallel jobs
    # 18 jobs, <= 6 in parallel, about 100 s each (lines "mX inj first ntoys"):
    #   17.0/16.85/16.7 x first 0,500,1000,1500 (500 toys); 10.0, 25.0, 40.0 (500); 17.0 with inj 1e-5, 3e-5, 6e-5 (500)
    xargs -P 6 -L 1 sh -c 'root -l -b -q -e "$B" "plots/all_report/x17c_toys.C+O($0, $1, $2, $3)" > /hdd3/.../x17c/toys_log2_m$0_i$1_f$2.txt 2>&1' < jobs.txt
    for r in 0.06 0.08 0.09 0.10 0.11 0.12 0.14 0.17; do                      # <= 6 in parallel, 1 min each
      root -l -b -q -e "$B" "plots/all_report/x17c_toycls.C+O(17.0, $r, 3000)" | grep x17c_toycls: >> /hdd3/.../x17c/x17c_toycls_m17.txt
    done
    root -l -b -q -e "$B" 'plots/all_report/x17c_summary.C+O'    # x17c_expected.txt, x17c_inject.txt, expected/toys_ul/inject figures
    python3 plots/all_report/x17c_frozen.py                        # /hdd3/.../x17/x17_frozen.txt, x17c_expected.tex, x17c_nuisances.tex

**Code:** `plots/all_report/`: `x17c_model.h`, `x17c_lib.h`, `x17c_disp.C`, `x17c_freeze.C`, `x17c_toys.C`,
`x17c_toycls.C`, `x17c_summary.C`, `x17c_frozen.py`.

**For X-d:** Sect. 10 of x17_frozen.txt. The building blocks are:
- `x17cl::Build`, `SetLine`, `MakeConfig(ReadConfig)`, `LimitBR`, `WindowX17c`;
- `x17c::ToyCLs` / `ToyUL`.

Load the masked rows from the sealed file into the same 35 x 30 vector; nothing else changes.
