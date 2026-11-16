# WP5 main extraction, stage B (prefix wp5) -- results note

Date 2026-11-14. Branch pp45_pippimepem (worktree /home/przygoda/HADES/FAT/FAT). Nothing committed, no git state changed;
ScanEngine.h, scan_*.C, the WP0/WP3/WP4 files and hists_v1.root untouched; FAT_simwt not used directly (only through the WP0 copies).

**Status in one line.** WP4 FAILED F-main, F-free and F-gold. So every number in this note is a **CROSS-CHECK**, labelled as
such in every table and figure. The main results are B4 and the peak side bands of WP6. This package supplies:
- the "other eta" inputs of WP6;
- the a_eta constraint and the background shapes for X17;
- the fits and macros for the WP7 variants;
- the ratios for WP8.

**Machine-readable output:** `/hdd3/przygoda/HADES/FAT_out/all_report/wp5/results.txt`. It holds parameters, the stat/MC/total
covariances of every fit, derived quantities, the per-bin model, window sums, X17 inputs, profiles and acceptance lines.

---------------------------------------------------------------------------------------------------
## 1. What was done, step by step (plain words)

1. **Event-group covariances** (`wp5_rows.C`).
   - One event can give several rows, so the counts of neighbouring cells are correlated.
   - For every data sample used here I counted, cell by cell, how many rows of the same event group fall in each pair of
     cells. Samples: S0 all, folds and periods; S0 with isBest == 1; pre, for F10 and for F5g.
   - Checks: the cell counts and the diagonal equal hists_v1.root to 2e-12; the S0 and pre F10 matrices equal WP3's
     wp3_rows.root exactly.
   - The variance of the F10 total is 1.19 x N (S0) and 1.26 x N (pre). With isBest == 1 it is 1.00, because there is one row
     per event.
2. **Templates.**
   - Simulation S0__nom from hists_v1.root, converted to expected raw rows: nb x 5984.2 (L = 6079.2 nb^-1 x 63/64).
   - The non-eta templates G_D and G_C are smoothed exactly as in WP4 (`wp4::SmoothFit`: exp(pol3) in x = (M4 - 0.55)/0.1, one
     common scale per histogram; F5g G_C of S2-S4 takes the shape of its F10 parent).
   - The eta classes C2-C5 and G_eta stay raw. Their bin variance goes into the fit weights (BB-lite).
   - 100 WP0 bootstrap replicas of the simulation are loaded the same way, and their smoothing is redone.
3. **Fitter** (`wp5_lib.h`).
   - It is the WP4 estimator written in a general form: Poisson/IRLS fixed point by Fisher scoring, with variance
     V = max(mu, 1) + k^2 var(CB) + BB-lite.
   - Checked against `wp4::NLFit` on the configuration that WP4 validated (F-main, E1 + E2, the same 700 bins): a_eta, kappa,
     rho_d and k_E1 agree to all printed digits (a_eta 0.393166 in both). The only chi2 difference is the penalty term of the
     constraints, which wp4 does not add.
4. **Errors of every fit.**
   - **stat** is a sandwich error with three parts:
     - N: the full event-group covariance of the counts, scaled to the model as in WP3;
     - CB: the full 700 x 700 covariance of the CB from the 200 like-sign replicas;
     - one part for each Gaussian constraint (k_E2 = 1 +- 0.05, g = 1 +- 0.30, and kappa = 1 +- 0.2 in the nominal fit).
     For a data subset (folds, periods), the CB covariance of "all" is scaled with the CB ratio of the bin (flagged).
   - **MC**: the spread of 100 refits, each with one simulation replica as the template set.
   - **total** = stat (+) MC.
   - The replica mean minus the nominal value is the template "attenuation" shift found in WP4. The value "bias-corrected" =
     2 v - mean(replicas) is WP4's RB procedure, given as information only.
5. **Fall-back of WP3 applied first.**
   - F-main as specified (F10 E1 + E2, k_E1 free) gives **k_E1 = 7.63 +- 0.52 (stat) +- 1.64 (MC)**. This is outside 0.8-3,
     the same as WP3's 13.2.
   - So the fall-back "drop E1, kappa = 1 +- 0.2" FIRED. The **nominal F-main uses F10 E2 only (M(ee) >= 10 MeV, 350 bins)**:
     - free: a_eta, rho_d, c_1..c_5, b, kappa_B;
     - constrained: kappa = 1 +- 0.2, g = 1 +- 0.3, k_E2 = 1 +- 0.05;
     - the background model is the D11 choice frozen in stage A (M1: c_j (1 + b x), common tilt).
6. **Fits run** (all on S0 and isBest >= 0 unless stated):
   - F-main nominal (E2) and as specified (E1 + E2);
   - variants for the fall-back "stage B chi2 > 1.5": b constrained to stage A, + quadratic term, per-slice tilts, kappa free,
     g free, M4 range 0.42-0.70 and 0.45-0.70, raw non-eta templates, exp(pol2) smoothing;
   - per-slice a_eta;
   - F-free on E2 and on E1 + E2;
   - F-gold on F5g:
     - free: alpha_2, alpha_5, g, c_j, b, k_E2;
     - kappa_B and a kappa (C3 + C4) fixed to F-main, and fixed to its replica fit inside each replica;
     - the kappa_B input error from +-1 sigma refits;
   - folds 0/1, with the templates x the row fraction 0.50018 / 0.49982;
   - pre, with pre templates and pre replicas;
   - data isBest == 1, with the isBest >= 0 templates (D2);
   - the 4 periods for the WP2 stability test:
     - background shape, kappa and g fixed to F-main, templates x the L share;
     - free: a_eta, rho_d, one overall background scale.
7. **Mass scale** (free s).
   - Refits with the WP0 templates at s = 0.99521, 0.99525, 0.99615, 0.99705 (and r = 0.89, 1.10).
   - The comparison uses the fixed-weight chi2: the weights of the nominal fit are frozen, so the minimum of this chi2 is exactly
     the fit's fixed point. The plain Pearson chi2, whose weights change with the model, is not minimal at the fit and gave
     inconsistent profiles; it was replaced.
   - Two versions:
     - the whole fit;
     - the peak region only: M4 0.5225-0.565, E2, with a_eta, kappa, g, b and kappa_B fixed to F-main and rho_d, c_j, k_E2 free.
   - s_hat comes from a parabola through 0.99525 / 0.99615 / 0.99705.
8. **Profiles for Fig. E3:** a_eta, rho_d, kappa (without its constraint) and g (without its constraint), all with the
   fixed-weight chi2.
9. **GLS chi2 of the nominal fit, as in WP3.** The data covariance (event groups + CB, full or diagonal) plus the MC covariance
   of the predicted spectrum. The latter is the spread of the 100 replica templates at the fitted parameters.

---------------------------------------------------------------------------------------------------
## 2. Numbers (all CROSS-CHECKS; S0, isBest >= 0, L = 6079.2 nb^-1; yields = raw rows x 64/63, not efficiency corrected)

### 2.1 Nominal F-main (E2 only, kappa = 1 +- 0.2), Table E1

Each entry is value +- stat +- MC.

| param | value | stat | MC | total | x 1.61 (chi2-scaled) |
|---|---|---|---|---|---|
| a_eta | 0.3823 | 0.0169 | 0.0396 | 0.0431 | 0.0694 |
| kappa (constrained) | 0.832 | 0.145 | 0.267 | 0.304 | 0.489 |
| rho_d | 1.315 | 0.115 | 0.226 | 0.254 | 0.409 |
| g (constrained) | 1.717 | 0.163 | 0.366 | 0.401 | 0.645 |
| c_1..c_5 | 1.092 / 1.105 / 1.090 / 1.153 / 0.984 | 0.025-0.035 | 0.115-0.178 | | |
| b [per 0.1 GeV/c2] | -0.0111 | 0.0046 | 0.0206 | 0.0211 | |
| kappa_B | 1.947 | 0.173 | 0.912 | 0.928 | |
| k_E2 | 0.982 | 0.044 | 0.052 | 0.068 | |

- **Fit quality:**
  - Pearson chi2 1158/341 = 3.40;
  - GLS with the MC covariance 883.7/341 = 2.59 at the fit, 2.16 minimised;
  - WP4 predicted about 2.8 from template noise alone.
- **Where the misfit sits** (Pearson by region; full table in wp5_extract.log):
  - 0.400-0.455: 406 for 55 bins (S2 122, S3 216);
  - CR-U 0.575-0.700: 420 for 125 bins;
  - B: 181 for 70 bins (S5 alone 89);
  - 0.5225-0.565: 125 for 85 bins.
- **Replica shifts (template attenuation, cf. WP4 RF):**
  - a_eta -0.050 (-13 %), Y_bump -16 %, Y_peak -14 %;
  - bias-corrected values: a_eta 0.432, Y_bump 22 472, Y_peak 5 755, rho_d 1.286 (information only).
- **Simulation-limited:** MC/stat = 2.3 (a_eta), 2.3 (Y_bump), 2.0 (Y_peak).

### 2.2 Yields and ratios, Table E2

Error columns: stat (event groups) / CB (replicas + k_E2) / nuisance (g, kappa) / MC; S10 = WP4 RF bias in scenario (a).

| quantity | value | stat | CB | nuis. | MC | S10 |
|---|---|---|---|---|---|---|
| Y_bump = Y2 + Y3 in B | 19 426 | 687 | 196 | 677 | 2 274 | 11.6 % |
| Y2 in B (C2) | 11 943 | 383 | 115 | 346 | 1 246 | |
| Y3 in B (C3) | 7 483 | 516 | 92 | 983 | 2 021 | |
| Y_peak = Y4 + Y5 in P | 5 045 | 298 | 86 | 181 | 720 | 6.6 % |
| Y4 in P (C4) | 1 847 | 127 | 23 | 243 | 496 | |
| Y5 in P (C5) | 3 198 | 272 | 82 | 62 | 571 | |
| Y2 / Y3 / Y4 / Y5, all M4 | 15 025 / 9 398 / 2 255 / 3 995 | | | | | |
| DR | 1.169 | 0.061 | 0.020 | 0.006 | 0.131 | |
| rho_d | 1.315 | 0.109 | 0.035 | 0.015 | 0.226 | 9.7 % |
| BR_A = 2.68e-4 rho_d | 3.52e-4 | 0.29e-4 | 0.09e-4 | 0.04e-4 | 0.61e-4 | 9.7 % |
| BR_G (F-gold) | 2.65e-4 | 0.34e-4 | 0.13e-4 | 0.05e-4 | 0.70e-4 (mc68) | 14.8 % |
| R_34 (F-free E1 + E2) | 0.917 | 0.179 | 0.044 | 0.015 | 0.241 | |
| BR_C (F-free E1 + E2) | 3.6e-4 | 0.9e-4 | | | 3.6e-4 (mc68) | 66.9 % |

- Visible cross sections (Y/L): bump 3.20 +- 0.16 (stat) +- 0.37 (MC) nb; peak 0.83 +- 0.06 +- 0.12 nb.
- **BR_G against PDG 2.68 +- 0.11 e-4:** -0.0 sigma. Bias-corrected 2.53e-4.
- **BR_A:** +1.2 sigma (total), +0.8 sigma chi2-scaled.
- **BR_C is WITHDRAWN** in the nominal E2 configuration. Without E1, F-free cannot determine the conversion classes:
  - alpha_3 = -0.07 +- 0.24, alpha_4 = 0.29 +- 0.19;
  - so R_34 = -4.2 +- 4.8 and BR_C is negative.
  The E1 + E2 value above is a labelled cross-check: R_34 = 0.92 +- 0.30 is consistent with 1, but E1 is the class whose
  model failed (k_E1 = 7.6).
- **F-gold parameters:**
  - alpha_2 = 0.388 +- 0.019 +- 0.046, the same as a_eta;
  - alpha_5 = 0.385 +- 0.054 +- 0.072;
  - g = 2.19 +- 0.18 +- 0.33;
  - chi2 902/342.
  - The kappa_B input error on BR_G is 0.006e-4 (negligible).

### 2.3 kappa vs kappa_B

- **Nominal fit:**
  - kappa = 0.83 +- 0.30. The constraint 1 +- 0.2 dominates: its part of the stat error is 0.13 of 0.145.
  - kappa_B = 1.95 +- 0.93 (MC-dominated).
  - Kappa free: kappa = 0.33 +- 0.52, and a_eta rises to 0.426 (+11 %).
- **F-main E1 + E2:** kappa = 0.50 +- 0.16, kappa_B = 1.11 +- 0.05. They differ by 0.61 +- 0.17. But in E1 an unmodelled
  component is absorbed by k_E1 = 7.6.
- **E1 prediction from the E2 parameters** (lower row of Fig. E1): data/model = 0.75-0.90. kappa_B = 1.95 multiplies G_C,
  which dominates E1. WP3 found 0.65-0.80.
- **Conclusion:** no data-driven statement on the conversion scales is possible in this configuration.

### 2.4 Variants (full list: tables_all/wp5_acceptance.txt; change relative to nominal)

| variant | GLS chi2 (at fit) | a_eta | Y_bump | Y_peak |
|---|---|---|---|---|
| nominal | 883.7 | 0.382 | 19 426 | 5 045 |
| b constrained to stage A (-0.068 +- 0.025) | 888.3 | -0.8 % | -0.8 % | -0.8 % |
| + quadratic q x^2 (S3) | 796.7 (min 738.3 vs 737.5) | **+39 %** | +32 % | +5 % |
| per-slice tilts (S3) | 901.8 | -2.7 % | -2.9 % | +1.3 % |
| kappa free | 895.7 | +11.5 % | -14.2 % | -12.9 % |
| g free, g = 2.56 (S8) | 947.1 | -6.1 % | -1.8 % | +2.3 % |
| range 0.42-0.70 (S5) | 734.3/321 | -7.8 % | -5.4 % | -8.2 % |
| range 0.45-0.70 (S5) | 541.6/291 (1.86) | -5.0 % | -2.5 % | -11.1 % |
| raw non-eta templates (S4) | 3389 | +30 % | -16 % | +28 % |
| exp(pol2) smoothing (S4) | 2019 | +39 % | +50 % | +95 % |

- The raw and pol2 variants are rejected by their chi2.
- The quadratic term is NOT rejected: its minimised GLS chi2 equals the nominal one. Yet it moves a_eta from 0.38 to 0.53.
- The bump amplitude is therefore degenerate with the curvature of the background. This is the largest single model effect
  and belongs in S3 (WP7). It confirms WP3's "EXTRAPOLATED" label for the bump.

### 2.5 Consistency fits

- **Per-slice a_eta (S1..S5):**
  - values 0.356, 0.394, 0.406, 0.353, 0.506; stat 0.030 / 0.020 / 0.024 / 0.027 / 0.076;
  - p = 0.92 with stat + MC, 0.11 with stat only, 0.94 for S1-S4.
- **Folds:**
  - a_eta 0.356 +- 0.019 / 0.395 +- 0.020 (z = -1.69);
  - Y_bump z = -1.57, Y_peak z = -1.71, DR z = -0.92;
  - joint test of all free parameters with the data covariance: chi2 10.8/12, p = 0.55.
- **pre:**
  - a_eta 0.383 +- 0.018 +- 0.040, rho_d 1.36, DR 1.20 +- 0.07 +- 0.15;
  - k_E2 = 0.68 +- 0.03 (the CB in pre is over-estimated by the spread rule);
  - GLS chi2/ndf 3.86.
- **Data isBest == 1** (ratio to isBest >= 0, same templates):
  - a_eta 0.867 +- 0.055, Y_bump 0.868 +- 0.063, Y_peak 0.888 +- 0.090, DR 1.023 +- 0.081;
  - the stat errors are an independent-sample upper bound, because the samples are nested; MC errors 0.010-0.019.
- **Periods** (background shape fixed; for WP2):

  | period | a_eta | DR |
  |---|---|---|
  | I | 0.367 +- 0.012 | 1.29 +- 0.10 |
  | IIA | 0.423 +- 0.031 | 0.68 +- 0.21 |
  | IIB | 0.383 +- 0.019 | 1.02 +- 0.14 |
  | III | 0.368 +- 0.015 | 1.21 +- 0.11 |

  - a_eta chi2 3.2/3, p = 0.36; DR chi2 8.4/3, p = 0.038 (> 0.01).
  - Period I contains the sector-2 deficit days (C_sec not applied).

### 2.6 Mass scale

- **Whole fit:** s_hat = 0.99653 +- 0.00019 (+1.95 sigma; +1.06 sigma chi2-scaled).
- **Peak region:** s_hat = 0.99500 +- 0.00057 (-2.01 sigma; -1.79 sigma chi2-scaled). This is 0.6 MeV/c2 at the peak.
- The two pull in opposite directions. Both lie inside the S6 variation +-0.0009.
- r = 1.10 is preferred over r = 1 by a fixed-weight delta chi2 of 31; r = 0.89 is disfavoured by 9. This means a wider
  simulated resolution. [v5 correction: this is OPPOSITE to D8, where the data peak (5.44 +- 0.60 MeV) is narrower than the
  simulated one (6.09 MeV), which is why the variant r = 0.89 = 5.44/6.09 exists. The r = 1.10 preference comes from outside
  the peak: in the peak region alone (0.5225-0.565, configuration of the local s fit) r = 0.89 and r = 1.10 lower the
  fixed-weight chi2 by 2.9 and 1.8, i.e. no significant preference. See "Verification", item V7.]

---------------------------------------------------------------------------------------------------
## 3. Figures and tables (draft captions)

- **Fig. E1** `docs/latex/figs_all/wp5_E1.pdf` (.png). *Stage-B fit F-main (CROSS-CHECK: it failed the WP4 closure).*
  - Data S0, isBest >= 0, L = 6079.2 nb^-1, in rows per MeV/c2.
  - Each column is one M(pi+pi-) slice (S1 < 0.34, S2 0.34-0.37, S3 0.37-0.39, S4 0.39-0.43, S5 >= 0.43 GeV/c2).
  - Top row, M(ee) >= 10 MeV/c2, is the fitted class. Pull = (N - mu)/sigma, with sigma = model variance x event-group factor
    + CB + template variance.
  - Bottom row, M(ee) < 10 MeV/c2, is NOT fitted (fall-back "drop E1"). It shows the prediction with the E2 parameters and
    k_E1 = 1, and the ratio data/model.
  - Black: data N. Stacked: CB = 2 sqrt(N++ N--) (red), smoothed non-eta simulation c_j (1 + b x)(G_D + kappa_B G_C) (grey),
    G_eta, C3 + C4, C2, C5. Simulation at SMASH sigma_prod x PDG BR x V3 weights, scaled by the fitted amplitudes.
- **Fig. E1 (spec version)** `wp5_E1_spec.pdf`: the same for the fit as specified, E1 + E2 with k_E1 = 7.6 free; pulls in both
  rows.
- **Fig. E2** `wp5_E2.pdf`. *Data minus the measured CB minus the fitted non-eta background, per slice, M(ee) >= 10 MeV/c2
  (blue).*
  - Errors: event-group variance of N plus the CB variance.
  - Compared with the fitted eta classes (stacked): C2 x a_eta, C3 + C4 x a_eta kappa, C5 x a_eta rho_d, G_eta x a_eta g.
  - Dotted: all eta classes at SMASH x PDG BR (a = kappa = rho = g = 1). Its spikes are single heavy simulated rows.
  - Last panel: sum of the slices. CROSS-CHECK.
  - [v5 correction: in `wp5_E2.pdf` the white legend of the S5 panel hides the S5 data of M4 0.40-0.57 (negative residuals
    near 0.48-0.53). Use `docs/latex/figs_all/v5_wp5_E2_corrected.pdf` (same numbers, legend moved to the last panel); the
    legend placement is also fixed in `wp5_extract.C` for any future rerun.]
- **Fig. E3** `wp5_E3.pdf`. *Profile scans of the nominal F-main.*
  - Delta chi2 with the fit weights frozen at the nominal fit (black), for a_eta, rho_d, kappa and g (the last two without their
    constraint).
  - Dashed: parabolas with the stat (sandwich) and the stat (+) MC errors.
  - Bottom middle/right: mass-transfer scale s. Points are refits with the WP0 templates at s = 0.99521 / 0.99525 / 0.99615 /
    0.99705; the dashed line is the parabola; red marks the nominal 0.99615. Whole fit and peak region 0.5225-0.565.
- **Fig. E4 (extra)** `wp5_E4.pdf`. *a_eta of every fit:* nominal, as specified, the variants, per slice, folds, periods, pre,
  isBest == 1, and alpha_2 of F-gold and F-free. Thick: stat; thin: stat (+) MC; grey band: nominal stat. All are cross-checks.
- **Table E1** `docs/latex/tables_all/wp5_tableE1.{txt,tex}`: parameters of F-main nominal / as specified / pre / isBest == 1
  (value +- stat +- MC; c = constrained) and the stat + MC correlation matrix of the nominal fit.
  - Largest correlations: c_j-c_k about 0.97; c_j-kappa_B about -0.9; a_eta-b 0.05; rho_d-kappa 0.04.
    [v5 addition: the largest correlations of a_eta are with g (-0.53) and kappa (-0.44), stat + MC; with c_j -0.37 to -0.41.]
- **Table E2** `wp5_tableE2.{txt,tex}`: yields and ratios with stat / CB / nuisance / MC / S10 / bias-corrected. The syst column
  is for WP7. Also the isBest == 1 ratios.
- **Acceptance and variants:** `wp5_acceptance.txt`.
- **Logs:**
  - `/hdd3/przygoda/HADES/FAT_out/all_report/wp5/wp5_extract.log` (every fit, the chi2 by region, the scans);
  - `wp5_rows.txt` (the covariance checks).

---------------------------------------------------------------------------------------------------
## 4. Acceptance criteria (spec WP5) and fall-backs

| criterion | value | verdict |
|---|---|---|
| chi2/ndf <= 1.2 | GLS 2.59 (2.16 minimised); Pearson 3.40 | **FAIL** -> errors x sqrt(chi2_GLS/ndf) = x 1.61, given as a separate column (`error_scale_chi2` in results.txt) |
| \|rho(a_eta, b)\| < 0.7 | 0.05 (stat + MC), 0.28 (stat) | **PASS** |
| per-slice a_eta p > 0.01 | 0.92 (stat + MC), 0.11 (stat) | **PASS** |
| no pull run > 4 bins in 0.52-0.56 | longest same-sign run S1 8, S2 3, S3 5, S4 6, S5 4; slice sum 4 | **FAIL** (3 of 5 slices; random signs give P = 0.37 per histogram, 0.90 for at least one of 5). S1 is low by 1.6-2.5 sigma in 0.52-0.55: a level offset of the background, not a peak-shape effect. |
| free s within 2 sigma of 0.99615 | whole fit +1.95 sigma; peak region -2.01 sigma (-1.79 with the chi2 scaling the spec prescribes) | **PASS, marginal** |
| folds p > 0.05 | joint p = 0.55; a_eta z = -1.69 | **PASS** |

**Fall-backs:**

| fall-back | status | consequence |
|---|---|---|
| WP4 gate | **fired** | F-main, F-free and F-gold are cross-checks; B4 + the peak side bands (WP6) are the main result. |
| k_E1 outside 0.8-3 (k_E1 = 7.6) | **fired** | E1 dropped; kappa = 1 +- 0.2. |
| stage B chi2/ndf > 1.5, and pull runs | **fired** | S3/S4 variants run. Best GLS chi2/ndf 1.86 (range 0.45-0.70); none reaches 1.5. Unresolved, so the side bands are the main result (already the case). |
| MC error > stat error | **fired** | "Simulation-limited"; a larger std production is requested. |
| R_34 off by > 2 sigma | not fired | In E2, R_34 is undetermined (+-4.8), so BR_C is withdrawn in the nominal configuration. |
| \|rho(a_eta, b)\| > 0.7 | not fired | b fixed to stage A shown as a variant: -0.8 %. |
| per-slice p < 0.01 | not fired | |
| corr(rho_d, kappa) > 0.9 | not fired | max 0.07 |
| direct-decay bias > 20 % | not fired | WP4 RF +9.7 % in (a); only the stress scenario (e) reaches 22.6 %. |

---------------------------------------------------------------------------------------------------
## 5. Caveats (flagged)

1. **E1 is dropped, so every yield over all M(ee) is extrapolated** from M(ee) >= 10 MeV with the simulated M(ee) shape.
   This holds for Y_bump, Y_peak and the Y_c over all M4.
   - Y3 and Y4 are conversions (mostly E1). They therefore rest on the constraint kappa = 1 +- 0.2: its part of the error is
     983 of Y3's 1114 stat.
2. **g = 1.72 +- 0.40** is pulled 2.4 sigma(prior) above 1. g = 2.2 in F-main E1 + E2 and in F-gold; g free gives 2.56.
   - G_eta (an eta in the event, but not a C2-C5 signal pair) is the largest eta-related term in E2 above 0.56.
   - For WP6: the effective G_eta scale is a_eta g = 0.66, not a_eta.
   - Stage A was run with a_eta(0) = 0.4 and g = 1. With dtheta/da_eta of stageA.txt, the shift 0.4 -> 0.66 on G_eta moves
     c_1..c_5 by -0.016 / -0.025 / -0.024 / -0.014 / -0.001 (up to 0.6 sigma stat, 0.1 sigma total) and b by +0.0007.
   - **I did NOT rerun wp3_stageA.C.** Its C2-C5 dependence is negligible: a_eta 0.38 vs 0.40 changes c_j by < 0.002. The g
     effect is quantified above for WP6 to decide.
3. **Bump-background degeneracy:** adding a quadratic term (equally good GLS chi2) raises a_eta by 39 %. The tilt b of the
   full fit (-0.011 +- 0.021) differs from stage A (-0.068 +- 0.025) by 1.7 sigma (total). Stage A takes it from CR-U only.
4. **The template attenuation of WP4 is present in the data fit.** Replica shift -13 % on a_eta, against WP4 RF -10 %. The
   bias-corrected numbers are given, but WP4 found that their coverage fails: the replica RMS under-estimates their MC error by
   about 1.5.
5. **Ratios use the robust MC spread.**
   - BR_C, R_34 and BR_G take mc68, half of the central 68 % replica interval. A few replicas with alpha_3 or alpha_4 near 0
     make the RMS meaningless: BR_C RMS 1.9e-3.
   - The replica covariances in results.txt are RMS-based.
6. **Subsets.**
   - For folds and periods the CB covariance is that of "all" scaled per bin by CB_sub/CB_all (no replicas exist for subsets).
   - The isBest == 1 ratios use an independent-sample bound for stat.
7. **Periods** use the L shares of `alldefs::kLper` (all runs). C_sec and f_win are not applied (WP2).
8. **Not run:**
   - exp(pol4) smoothing (`wp4::Smooth` has 4 coefficients only); exp(pol2) is shown instead and is rejected by its chi2.
   - The DIRECT and PEAK selections and the vz -150/-170 variants: their templates exist in hists_v1.root, and `wp5_lib.h`
     runs them with `LoadTempl(fh, "<sel>", "nom", 1, T)` / `LoadData`. Left to WP7.
9. **Profiles:** the Pearson chi2 with model-dependent weights is not minimal at the fit, so the profiles and the s scan use
   the fixed-weight chi2. The parabola width then corresponds to the fitter covariance (BB-lite), not to the sandwich.
10. **Blinding:** only F10 and F5g (coarse) were used. No X2D or B_mee data cell of 12-28 MeV was read.

---------------------------------------------------------------------------------------------------
## 6. For the other packages (results.txt keys)

- **WP6 (other eta terms):**
  - `values FMAIN` with names a_eta kappa rho_d g c_1..c_5 b kappa_B k_E2;
  - `cov_stat|cov_mc|cov_total FMAIN`;
  - `simwin` (template window sums, raw rows at unit amplitudes);
  - `bkgwin` (F-main background in B, P, SBlo, SBhi, CR-U per class; E1 = prediction);
  - `bin` lines (per F10 bin: N, CB, k CB, background with stat and MC errors, C2..Geta, mu, in_fit).
- **X17:**
  - `x17_aeta 0.382284 stat 0.016925 total 0.043092 bias_corr 0.431915`;
  - `x17_g`, `x17_kappaB`, `x17_b`;
  - `x17bkg E1|E2` = background per X2D M4 bin (2.5 MeV, 0.5225-0.6100), summed over slices, with the MC spread (E1 =
    prediction);
  - `smooth` = exp(pol3) parameters per histogram.
- **WP7:** every fit block includes values and errors. The macro `wp5_extract.C` with `wp5_lib.h` refits any WP0 variant:
  `LoadTempl(fh, sel, var, 1, T)`, `SpecMain(Opt)`, `Run(...)`.
- **WP8:** `derived FMAIN` (Y*, DR, BR_A, sigvis), `derived FGOLD BR_G`, `derived FFREE_SPEC R_34/BR_C`, `best_ratio`.
- **WP2:** `derived PERIOD_<p>` and `periods_aeta`, `periods_DR`.

---------------------------------------------------------------------------------------------------
## 7. Exact commands (from /home/przygoda/HADES/FAT/FAT; ACLiC build dir aclic_wp5)

    B='gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_wp5", kTRUE)'
    root -l -b -q -e "$B" 'plots/all_report/wp5_rows.C+O'      # 12 s, 2.7 GB: wp5/wp5_egcov.root, wp5_rows.txt
    root -l -b -q -e "$B" 'plots/all_report/wp5_extract.C+O'   # 25 s, < 3 GB: results.txt, log, tables, figures

- Code: `plots/all_report/wp5_rows.C`, `wp5_lib.h` (includes wp4_lib.h), `wp5_extract.C`.
- Outputs: `/hdd3/przygoda/HADES/FAT_out/all_report/wp5/{results.txt, wp5_extract.log, wp5_extract_stdout.txt, wp5_egcov.root,
  wp5_rows.txt}`; `docs/latex/tables_all/wp5_{tableE1,tableE2}.{txt,tex}`, `wp5_acceptance.txt`;
  `docs/latex/figs_all/wp5_{E1,E1_spec,E2,E3,E4}.{pdf,png}`.

---------------------------------------------------------------------------------------------------
## Verification (independent verifier, prefix v5, 2026-11-14)

**Verdict.** The WP5 numbers are REPRODUCED by independent code. I found no error in the fits, the yields, the units, the
WP4 gate or the use of stage A. I made three corrections: two in the text and one in a figure. I add two qualifications to
the acceptance values and one important caveat on how much the estimator choice moves the results. All F-main / F-free /
F-gold numbers stay CROSS-CHECKS, as WP5 states.

Comparison table: `docs/latex/tables_all/v5_verification.txt`. Logs: `/hdd3/przygoda/HADES/FAT_out/all_report/v5/`.

### V1. Independent reproduction (own code, same histograms)
- **Input.** `v5_dump.C` writes hists_v1.root (data S0/pre/S0_best, sim S0/pre nominal + 100 replicas + s/r variants) and
  the two event-group covariance files as text. Nothing from wp4_lib.h, wp5_lib.h or s5_blind_ana_fit.h is used.
- **Smoothing.** My own exp(pol3) scaled-Poisson fit (one scale per histogram, 16-point Gauss-Legendre bin integrals, scipy
  BFGS, the same order and parent rules) gives the same smoothed bin integrals as WP5 within 2e-6 (all 30 fits).
- **Fit.** The fixed point of the spec estimator (V = max(mu,1) + k^2 var(CB) + BB-lite) is solved as iterated weighted
  least squares with scipy's trust-region solver, which is a different algorithm from WP5's Fisher scoring.

| quantity | WP5 | v5 |
|---|---|---|
| a_eta | 0.382284 | 0.382283 |
| kappa | 0.83234 | 0.83234 |
| rho_d | 1.31520 | 1.31521 |
| g | 1.71677 | 1.71676 |
| c_j, b, kappa_B, k_E2 | | equal to 1e-5 |
| Y_bump | 19 426.1 | 19 426.1 |
| Y_peak | 5 044.8 | 5 044.8 |
| Y2/Y3/Y4/Y5 in windows and all M4 | | equal to < 1e-5 |
| DR | 1.16896 | 1.16896 |
| BR_A | 3.5247e-4 | 3.5248e-4 |
| sigma_vis bump / peak | 3.1955 / 0.8298 nb | 3.1955 / 0.8298 nb |

- **Errors.**
  - The event-group (N), CB-replica and constraint parts of the sandwich agree to 3-4 digits. The CB part differs by
    0.3 %, from ddof = 1 against WP5's 1/n.
  - With my own refits of the 100 replicas (own smoothing), the MC errors and the replica shifts agree to 4-5 digits for
    every parameter and yield. Example: a_eta 0.03963 / -0.04963.
- **F-gold.**
  - alpha_2 = 0.38825, alpha_5 = 0.38457, g = 2.1850, BR_G = 2.6546e-4 +- 0.368e-4 (stat): identical.
  - MC rms 6.58e-5 against 6.64e-5, and mc68 6.95e-5 against 7.02e-5 (1 %, from the replica convergence tolerance).
  - Pearson chi2 902.9/342: identical.
- **Other numbers reproduced:**
  - Pearson 1158.96/341 and Baker-Cousins 2221.3;
  - the Pearson split by region (405.7 / 181.5 / 124.9 / 420.8);
  - the pull-run lengths 8/3/5/6/4;
  - corr(a_eta, b) = 0.276 (stat);
  - per-slice a_eta 0.356/0.394/0.406/0.353/0.506 (stat 0.030/0.020/0.024/0.027/0.076), p = 0.109 with the full stat
    covariance;
  - the quad variant (a_eta 0.5317, Y_bump 25 638).

### V2. Units and inputs
- **Units:**
  - templates nb x 5984.2125 (= 6079.2 x 63/64);
  - yields x 64/63, so that Y/L = a x sigma_vis,sim;
  - sigma_vis = Y/6079.2.
  All correct.
- **Window sums.** The sim "win" B window equals the sum of the F10 cells (C2: 30 754.03 raw rows both ways). The data
  N in B is 111 212 both ways.
- **WP0 numbers.** Data S0, window P: N = 55 947 and CB/N = 0.0586, as in the WP0 acceptance.
- **Event-group covariance.**
  - The WP3 file (egcov_S0) and the WP5 file (egcov_S0_all_F10) are identical.
  - Their diagonal equals the hists_v1 N errors^2 to 8e-12.
  - var(F10 total)/N = 1.194.

### V3. Stage A and the WP4 gate
- **WP4 gate.** wp4/verdict.txt says F-main, F-free and F-gold DECISION = FAIL. Every WP5 table, figure and results line
  carries "CROSS-CHECK". D20 and the fall-back "estimator fails WP4 -> cross-check only" are respected.
- **Stage A in WP5.** The D11 choice frozen in stage A (M1, common tilt) is the nominal stage-B model. The background
  stays free in stage B, as specified. The variant "b constrained to stage A" uses stageA.txt b = -0.0683 +- 0.0249
  (total), which is correct.
- **Stage A not rerun at the WP5 value** (spec WP3; WP5 caveat 2). My own stage-A refit (`v5_stageA.py`) uses CR-U
  0.575-0.700, E2, M1, with the eta terms fixed to the F-main values instead of a_eta = 0.4, g = 1:
  - c_j shift by -0.0095/-0.0176/-0.0157/-0.0072/+0.0040 (at most 0.4 sigma stat);
  - b shifts by +0.0004 and kappa_B by -0.085;
  - the stage-A background in the B window (E2) drops by 2.3 % (38 537 -> 37 648 rows), almost entirely through g. WP6 should
    take this into account.
- **Smoothing mismatch (flag for WP6).** With WP4/WP5 smoothing, my stage A does not exactly reproduce the frozen WP3
  values: c_1 0.952 against 0.988, kappa_B 3.29 against 3.06, b -0.066 against -0.068. WP3's smoothing starts each fit
  at a per-histogram threshold (second column of the `smooth` lines in stageA.txt); WP4/WP5 smooth from 0.400. As a
  result, WP6's B_A (WP3 smoothing) and WP5's background (WP4 smoothing) come from slightly different templates: about
  3.5 % on c_j and about 1 % on B_A(B).

### V4. Estimator dependence (caveat; not an error)
On the same histograms (Table v5_verification.txt):

| estimator | a_eta | Y_bump | Y_peak | BR_G |
|---|---|---|---|---|
| WP5 IRLS (reference) | 0.3823 | 19 426 | 5 045 | 2.65e-4 |
| Poisson likelihood with a Conway BB-lite nuisance per bin (variance fixed at the fit), the textbook form of D12 | 0.3813 (-0.3 %) | -1.1 % | -0.4 % | 2.69e-4 |
| plain Poisson likelihood without template variance | 0.3349 (-12 %) | -5 % | -4 % | 2.32e-4 (-13 %) |
| WP5's own GLS (data + replica MC covariance), minimised; WP5 computes it but reports only its chi2 | 0.485 (+27 %) | 24 228 (+25 %) | 5 541 (+10 %) | |
| hybrid GLS, minimised | 0.513 (+34 %) | 26 335 (+36 %) | 6 193 (+23 %) | |

- The BB-lite likelihood confirms that the WP5 IRLS is a faithful Poisson + BB-lite fit.
- In both GLS fits, rho_d falls to 1.07-1.12 and **g falls from 1.72 to 0.47-0.50**.
- **Conclusion.** Once the correlated template noise of the smoothed background may absorb shape (GLS), the decomposition
  bump / G_eta / background changes by about 2.5-3 sigma_MC. This is the same freedom that the quadratic variant (+39 %)
  exposes. With the hybrid GLS minimised, nominal and quad differ by only Delta chi2 = 1.2 (567.8/341 against 566.6/340),
  so the quad term is not required once the MC covariance is in the fit.
- **Consequences:**
  - The F-main a_eta, Y_bump and g are model- and estimator-limited at the 25-35 % level. The replica MC error of 10 %
    does not cover this; WP7 needs an S3/S4 entry for it.
  - WP6 takes the "other eta" term a g G_eta from F-main. In window B it is 3 853 rows with g = 1.72, but about 1 400 rows
    at the GLS solution. g is therefore not robust (0.47-2.56 over estimators and variants), and WP6/WP7 must vary it
    accordingly (S8).

### V5. Goodness of fit (qualification of the acceptance value)
- **Reproduced.** GLS-rep (WP5 definition) is 883.6/341 = 2.59.
- **Not converged in the number of replicas.**
  - With the first 50 replicas it is 1136/341 = 3.33.
  - The replica covariance has rank 99 for 350 bins and cannot represent the bin-by-bin noise of the raw signal templates
    (C2-C5, G_eta).
- **The exact diagonal (BB-lite) variance of these templates is confirmed by the replicas:**
  - trace 295.8k against 293.4k;
  - variance of the sum over bins / sum of the variances = 0.97, so there are no bin-to-bin correlations to speak of.
- **Hybrid GLS.** With that variance, plus the replica covariance of the smoothed background (rank 54, adequately
  sampled), the chi2 is 648/341 = **1.90** at the fit (2.00 with 50 replicas) and 1.67 minimised.
- **Consequences:**
  - The acceptance "chi2/ndf <= 1.2" still FAILS, and the fall-back "> 1.5" still FIRES.
  - The error scale factor sqrt(chi2/ndf) = 1.61 is conservative. With the hybrid GLS it is 1.38.
  - The same bias affects the GLS values of the variants (wp5_acceptance.txt).

### V6. Free s (qualification)
- **Reproduced.** Whole fit: s = 0.99653 +- 0.00019 (+1.95 sigma). Peak region: 0.99500 +- 0.00057 (-2.01 sigma). A
  4-point parabola gives -1.98 sigma.
- **The peak-region minimum is not bracketed.**
  - It lies below the lowest template (0.99521).
  - chi2(0.99521) is 0.09 ABOVE chi2(0.99525), so the curve does not keep falling.
  - What the data show is chi2(0.99525) - chi2(0.99615) = -3.8, i.e. a preference of about 2 sigma for a lower s.
- **Missing pieces.** The test has no MC error on s. The whole fit and the peak region disagree by 2.5 sigma.
- **Verdict.** "PASS only with the chi2-scaled error" is literally correct. I would label it MARGINAL / NOT DEMONSTRATED,
  because templates below 0.995 would be needed to bracket the minimum.

### V7. Resolution r (text error, corrected in Sect. 2.6)
- WP5 wrote that the whole-fit preference for r = 1.10 is "consistent with D8". It is the opposite: D8 has the data peak
  narrower (5.44 MeV) than the simulation (6.09 MeV), which is why r = 0.89 = 5.44/6.09 exists.
- Peak region only: r = 0.89 gives -2.9 and r = 1.10 gives -1.8 against r = 1. So there is no preference, and the
  Delta chi2 = -31 of r = 1.10 comes from outside the peak.

### V8. Figures
- **E1, E3, E4:** content checked against the numbers.
- **E2 (error, corrected).** In the S5 panel the white-filled legend hid the S5 data of 0.40-0.57, where the residuals are
  negative (about -5 to -8 rows/MeV in 0.48-0.53).
  - Corrected figure, same numbers (from results.txt `bin` lines and hists_v1), legend moved to the last panel:
    `docs/latex/figs_all/v5_wp5_E2_corrected.{pdf,png}` (`v5_figE2_prep.py`, `v5_figE2.C`).
  - wp5_E2.pdf was NOT overwritten (WP5 not rerun).
- **E1 (cosmetic, not fixed).** The legend text "non-eta background (G_D + kappa_B G_C)" runs into the S5 data, and the
  PNG renders "isBest >= 0" with a gap.

### V9. Changes made to WP5 files (all marked "[v5 ...]")
1. `wp5_results.md` Sect. 2.6: correction of the r statement (V7).
2. `wp5_results.md` Sect. 3, Table E1 item: the largest correlations of a_eta, which are with g (-0.53) and kappa (-0.44)
   (stat + MC).
3. `wp5_results.md` Sect. 3, Fig. E2 item: pointer to the corrected figure.
4. `plots/all_report/wp5_extract.C`, Fig. E2: the legend is moved to the "all slices" panel and made transparent.
   - The edit compiles: ACLiC into aclic_v5, not executed.
   - The outputs were not regenerated, so wp5_E2.pdf still has the old legend.
No number of WP5 was changed. results.txt and the tables are untouched.

### V10. Acceptance criteria (spec WP5), as verified

| criterion | WP5 | v5 |
|---|---|---|
| chi2/ndf <= 1.2 | FAIL (GLS 2.59) | FAIL confirmed; best estimate 1.90 at the fit, 1.67 minimised; scale x 1.38 instead of x 1.61 |
| \|rho(a_eta, b)\| < 0.7 | PASS | PASS (0.28 stat, 0.05 stat + MC) |
| per-slice a_eta p > 0.01 | PASS | PASS (p = 0.109 stat) |
| no pull run > 4 in 0.52-0.56 | FAIL | FAIL (8/3/5/6/4) |
| free s within 2 sigma | PASS (scaled) | MARGINAL: -2.0 sigma unscaled, minimum not bracketed |
| folds p > 0.05 | PASS | not refitted; the numbers are consistent |

Fall-backs that fired (WP4 gate; k_E1 = 7.63; stage-B chi2 > 1.5; MC > stat with MC/stat 2.34) are all confirmed. I
recommend no new fall-back. V4 belongs in WP6/WP7.

### V11. Commands (from /home/przygoda/HADES/FAT/FAT)

    root -l -b -q 'plots/all_report/v5_dump.C'               # 1 s: text dump -> /hdd3/.../all_report/v5/
    python3 plots/all_report/v5_fit.py 100                   # ~12 min: v5_fit.log, v5_fit.json (fits, errors, replicas, F-gold, s scan)
    python3 plots/all_report/v5_sscan.py                     # 10 s: v5_sscan.log (s and r, whole fit and peak region)
    python3 plots/all_report/v5_gls.py                       # 4 min: v5_gls.log (GLS variants, replica-rank test)
    python3 plots/all_report/v5_quad.py                      # 8 min: v5_quad2.log (nominal vs quad, GLS minimised fits)
    python3 plots/all_report/v5_stageA.py                    # 5 s: v5_stageA.log (stage A at the WP5 eta values)
    python3 plots/all_report/v5_figE2_prep.py && root -l -b -q 'plots/all_report/v5_figE2.C'   # corrected Fig. E2
