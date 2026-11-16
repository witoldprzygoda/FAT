# WP4 stress closure on simulation toys (prefix wp4) -- results note

Date 2026-11-14. Branch pp45_pippimepem (worktree /home/przygoda/HADES/FAT/FAT). No git command that changes state was run;
ScanEngine.h, scan_*.C, README.md and the WP0 files were not touched; FAT_simwt only read; no data yield was computed (D20):
the data enter only through the like-sign counts (CB-like component) and the CB replica variances.

**Verdict file:** `/hdd3/przygoda/HADES/FAT_out/all_report/wp4/verdict.txt` (machine-readable, lines ROW / EST / MC / S10 / GATE / X17).

---------------------------------------------------------------------------------------------------
## 0. Verdict in short

| estimator | SPEC test (BA) | method (AA, BB) | realistic MC (RF) | RF + bootstrap bias corr. (RB) | sigma_MC / sigma_stat | DECISION (as specified) |
|---|---|---|---|---|---|---|
| F-main | FAIL | FAIL (only (f) Y_bump -6.0 % in AA) | FAIL (a_eta -10 %, Y_bump -12 %) | FAIL (coverage 0.51-0.54) | 3.0 (a_eta), 3.3 (Y_bump) | **FAIL** |
| F-free | FAIL | FAIL (BR_C in (e), (f)) | FAIL | FAIL | 3.0 | **FAIL** |
| F-gold (BR_G) | FAIL | PASS | FAIL ((e) +23.9 %, (f) +22.4 %; coverage b+ 0.54, h 0.46) | FAIL marginal (coverage 0.77, 0.78) | 1.6 | **FAIL** |
| B4 (Y_bump) | FAIL | PASS | FAIL marginal ((e) -5.8 % vs 5 %) | PASS | 3.6 | **FAIL marginal**; PASS with bias correction |
| peak side bands (Y_peak) | FAIL | PASS | PASS | FAIL marginal (h: coverage 0.57) | 0.42 | **PASS** |
| X17 | not tested here | | | | | delegated to X-b (x17ab) |

Main findings:
1. **The SPEC test as written cannot be passed by any estimator.** With truth from half B and templates from half A, the
   two halves of the smoothed non-eta background differ by 5-10 % per 50 MeV band (MC statistics), and this alone moves
   a_eta by +32 % (BA) / -37 % (AB), B4 by +27 % / -31 %. The two halves differ by twice the MC error of the real analysis
   (variance 4 x var_full). The method itself closes exactly when the truth and the templates come from the same half
   (AA, BB: Asimov bias < 0.03 %, chi2 = 0.0).
2. **Realistic test (RF: truth = full sample, templates = bootstrap replicas, error = stat (+) sigma_MC):**
   - the raw signal templates (C2-C5, G_eta; N_eff 10-200 per histogram) bias every template fit DOWN through
     errors-in-variables attenuation: a_eta -10.1 %, Y_bump(F-main) -11.6 %, Y_peak(F-main) -6.6 %, BR_G +14.8 %.
     BB-lite does not remove it;
   - the smoothed background templates make the spread: sigma_MC = 8.8 % (a_eta), 5.6 % (Y_bump), 9.2 % (Y_peak),
     9.2 % (B4), but 3.7 % for the peak side bands;
   - the analysis is **simulation-limited** for F-main, F-free, F-gold and B4 (sigma_MC 1.6-3.6 x sigma_stat).
3. **Only the peak side bands pass everything as specified.** B4 fails by 0.8 point in scenario (e) only (within 2 sigma
   of the toy precision) and passes once the "other eta" terms it takes from F-main are bootstrap bias-corrected.
4. **Expected chi2/ndf from template MC noise alone** (perfect physics model, full-sample templates): F-main 2.83,
   stage A 2.10 (scenario a); 6.1 / 2.5 for a_eta = 1. A data chi2/ndf of 2-3 in WP3/WP5 is therefore NOT by itself
   evidence of background mismodelling. With half-sample templates: 8.8 / 4.7.

---------------------------------------------------------------------------------------------------
## 1. What was done, step by step (plain words)

1. **Simulation templates on toy cells** (`wp4_fill.C`, 30 s, 8 GB).
   - The simulated rows (std + lep + eta, owned rows, mode A, w x V3, selection S0) are put into "toy cells":
     - the 70 F10 mass bins, split at every window edge that is not an F10 edge (P, side bands, X17c): 77 cells;
     - x 5 M(pipi) slices x 3 M(ee) classes (< 10, 10-20, >= 20 MeV). F10 E2 is the sum of the last two; F5g is the last one.
   - With these cells every window of the estimators (B, P, side bands, CR-U, CR-K) is an exact sum.
   - Filled separately for half A (even event group) and half B (odd) and for four mass transfers: nominal; (f) s + 0.0009
     with r = 0.89; and, as checks, sUp and r089 alone.
   - Check: all 36 668 comparisons with hists_v1.root agree (F10 content and event-group errors, F5g, windows, totals; worst 1.7e-13).
2. **Smoothing of the non-eta templates** (D13; `wp4_lib.h` SmoothFit).
   - exp(pol3) in x = (M4 - 0.55)/0.1, one fit per F10 histogram (slice x E1/E2) and per F5g slice, integrated over each bin.
   - The likelihood counts each histogram as a scaled Poisson sample with ONE scale (sum of variances / sum of contents).
   - Correction made on the way: per-bin scales (variance/content of each bin) biased the normalisation low by 10-40 %,
     because a bin hit by a heavy std event group gets a small weight. The one-scale fit preserves the total exactly.
   - Sparse histograms: N_eff 20-50 get exp(pol2); below 20 the shape of the parent histogram is used, scaled to the
     histogram's own content (parent: for F5g G_C the F10 E2 G_C of the slice; for an omega split its total).
     This keeps the kinematic thresholds of S4/S5. Only F5g G_C of S2-S4 (N_eff 11-18) uses this rule.
3. **CB in the toys.**
   - The data like-sign counts N++ and N-- are taken per (2.5 MeV bin, slice, M(ee) class) from the S0 G3 histograms;
     below 0.40 and above 0.85 they come from fine_M4.
   - The WP0 rule is rebuilt from them: 2 sqrt(N++ N--) in 10 MeV cells, interpolated, spread by f(cell | 50 MeV band).
   - It reproduces CB_F10, CB_F5g and the window CBs of hists_v1 to 3.5e-16.
   - In every toy:
     - the like-sign counts are Poisson-fluctuated and the CB is rebuilt;
     - the true CB-like part of N is the data CB;
     - the fits use the replica variance of the CB (CBrep) per bin.
4. **Truth of a toy** (expected raw rows at L = 6079.2 nb^-1):
   - data CB;
   - plus the smoothed non-eta background x c (1 + b x), with the tilt integrated exactly over each cell;
   - plus a_eta [C2 + kappa (C3 + C4) + rho_d C5 + g G_eta] from the raw simulated rows.
   - Toy N = Poisson of the truth in every cell.
5. **Fitter.** The model is not linear in its parameters (c_j b, c_j kappa_B, a_eta g). It is fitted by Fisher scoring:
   - every step is one call of `bla::LinFit` (s5_blind_ana_fit.h) on the linearised model;
   - variance V = max(mu, 1) + k^2 var(CB) + sum coef^2 var(template) (BB-lite on the raw signal templates);
   - steps are limited to 3 sigma; convergence is reached when every step is below 1e-5 sigma (about 20 iterations);
   - about 160 000 runs of the estimator chain: 30 000 toys of BA, AA, AB, BB and RF, 1 200 replica fits, and
     6 000 toys x 21 template sets in RB. No F-main, F-free, stage-A or side-band fit failed. F-gold failed in 3 of 500
     toys of the swap configuration in the null scenario.
   - Errors: the covariance of the fitter (with BB-lite), the Fisher error with the data variance only, and the
     **sandwich** error C (J^T V^-1 Vd V^-1 J) C. The sandwich is the real statistical error of the BB-weighted fit.
     Without template mismatch it gives pull widths 0.83-1.12, mostly 0.9-1.07 (AA, BB, scenarios a-d, h). The BB-lite
     covariance over-covers (0.95), and the Fisher data-only error under-covers (0.56) for a_eta.
6. **Estimators** (implementation in `wp4_lib.h`, RunEstimators):
   - **F-main** on F10 (700 bins), 13 parameters.
     - Internally: beta2 = a_eta, beta34 = a_eta kappa, beta5 = a_eta rho_d. This is the same model, regular at a_eta = 0.
       kappa and rho_d come by the delta method.
     - Free: c_1..c_5, b, kappa_B, k_E1. Constrained: g = 1 +- 0.30, k_E2 = 1 +- 0.05.
     - Y_bump = (64/63)(a C2 + a kappa C3) in B; Y_peak = (64/63)(a kappa C4 + a rho C5) in P.
   - **F-free**: alpha_2..alpha_5 free; G_eta is multiplied by alpha_2 g.
     BR_C = 2.68e-4 alpha_5 (alpha_3/alpha_2)/alpha_4.
   - **Stage A**: F10 bins of CR-U (0.565-0.700, all slices) and CR-K (S5, 0.400-0.5225).
     - Free: c_j, b, kappa_B, k_E1; k_E2 constrained.
     - Eta terms fixed to F-main (the spec's "rerun at the WP5 value").
   - **F-gold** on F5g (350 bins): alpha_2, alpha_5, g, c_j, b and k_E2.
     - a kappa (C3 + C4) and kappa_B are fixed to F-main.
     - BR_G = 2.68e-4 alpha_5/alpha_2.
   - **B4**: Y = (64/63) {(N - CB)_B - [(N - CB)_U - eta_U] Bk_A(B)/Bk_A(U) - a (kappa C4 + rho C5 + g G_eta)_B}.
     - U = CR-U; Bk_A = stage-A background.
     - Deviation (flagged): the eta content of U (G_eta is 6 % of CR-U) is subtracted before the extrapolation. The
       formula as written does not do it, which would add about +25 % to Y_bump.
   - **Peak side bands**: Y = (64/63) [(N - CB)_P - R' (N - CB)_SB] / (1 - R' f).
     - R' = [Bk_A + a(C2 + kappa C3) + a g G_eta](P) / [same](SB).
     - f = (C4 + C5)(SB) / (C4 + C5)(P): the unweighted shape, as "iterated out with C4 + C5" of WP6. This is the closed
       form of the iteration.
     - With the fitted kappa and rho, f is undefined when the fitted peak is near zero (seen in scenario g).
   - **Errors of B4 and the side bands** are propagated linearly:
     - Poisson counts;
     - the window CB covariance from the 200 CB replicas;
     - the stage-A and F-main covariances, through numerical gradients.
     The correlations between these blocks are neglected; the toys show pull widths 0.94-1.06.
7. **Scenarios** (truth relative to the fit model):
   - (a) background x 1.35, tilt b = -0.15, a_eta 0.4, rho_d = kappa = g = 1;
   - (b-/b+) rho_d 0.5/1.5;
   - (c-/c+) kappa 0.8/1.2;
   - (d-/d+) g 0.5/1.5;
   - (e) background with omega/rho0 pions x 0.3 and the rest x 1.6 (tilt kept);
   - (f) truth with s = 0.99705 AND r = 0.89 (both at once), templates nominal;
   - (g) a_eta = 0;
   - (h) a_eta = 1;
   - (i, informational, not gated) the true CB distributed like the per-cell geometric mean of the band (S1 variant).
   500 toys per scenario and configuration, plus one Asimov data set.
8. **Configurations** (why several: see Sect. 5, caveat 1):
   - **BA** (SPEC): truth half B, templates half A;
   - **AB**: the swap;
   - **AA, BB** ("method"): truth and templates from the same half, so there is no MC mismatch;
   - **RF** ("realistic MC"): truth = the full sample. Toy t is fitted with templates from WP0 bootstrap replica t mod 100
     (hists_v1 `_rep`), so the template-to-truth difference has exactly the MC noise of the real analysis.
     sigma_MC = RMS of the 100 replica fits of the Asimov data; error = sandwich (+) sigma_MC;
   - **RB**: RF plus the bootstrap bias correction est = 2 est_p - mean_c est_pc.
     - Parent p = toy mod 100. Its 20 children are bootstrap replicas of the parent, k_pc ~ Poisson(k_p), made by
       `wp4_fill_boot.C`.
     - Error = sandwich (+) RMS of the children.
     - This is the procedure the data analysis can apply with its own replicas. It is not in the spec.
9. **Criteria** (SPEC WP4):
   - |bias| < max(0.5 sigma_stat, f |true|), with f = 5 % (a_eta, Y_bump), 10 % (Y_peak), 20 % (rho_d, BR);
   - 68 % coverage 0.60-0.76;
   - null scenario (g): |mean pull| < 0.5;
   - at most 5 failed fits in 500.
   Interpretations (flagged):
   - coverage is required only where the model contains the truth (a, b, c, d, g, h). In (e) and (f) a bias up to the
     allowed 5 % is up to 2 sigma for Y_bump, so coverage cannot hold there;
   - the ratios (rho_d, BR) are not gated in (g);
   - sigma_stat = toy RMS (BA, AA, BB) or the median sandwich error (RF, RB, whose RMS contains the MC part);
   - "FAIL marginal" = every failed criterion lies within 2 sigma of its limit: the bias error of the mean (with the MC
     part counted with 100 degrees of freedom), or the binomial error 0.021 of the coverage.
10. **Decision rule** (stated before the RB run; the SPEC result is reported next to it):
   - an estimator passes if it passes both method configurations (AA and BB) and the realistic test RF;
   - "with bias correction" uses RB instead of RF;
   - MC_LIMITED = sigma_MC > sigma_stat.

---------------------------------------------------------------------------------------------------
## 2. Numbers (selection S0, isBest >= 0, L = 6079.2 nb^-1, yields in raw rows x 64/63; 500 toys each)

### 2.1 Table V1 (docs/latex/tables_all/wp4_V1.txt; scenario (a) unless stated; per cent of the true value)

| estimator | quantity | sigma_stat | sigma_MC | MC shift | BA | AB | method | RF | RB | worst method | worst RF | worst RB |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| F-main | a_eta | 2.9 | 8.8 | -10.1 | +32.1 | -36.9 | +0.1 | -10.1 | -1.7 | -3.2 (e) | -11.9 (e) | -4.4 (f) |
| F-main | Y_bump | 1.7 | 5.6 | -11.5 | +8.0 | -33.4 | +0.1 | -11.6 | -1.8 | -6.0 (f) | -14.5 (f) | -5.4 (f) |
| F-main | Y_peak | 4.8 | 9.2 | -6.6 | +0.6 | -8.1 | +0.2 | -6.6 | -2.4 | +4.1 (e) | -9.2 (f) | -4.9 (f) |
| F-main | rho_d | 10.0 | 19.4 | +9.3 | -12.4 | +80.9 | +0.2 | +9.7 | -2.5 | +14.2 (e) | +22.6 (e) | +9.6 (e) |
| F-free | alpha_2 | 2.9 | 8.6 | -12.2 | +24.6 | -41.6 | +0.1 | -12.2 | -3.3 | -5.0 (f) | -14.3 (f) | -6.5 (f) |
| F-free | Y_bump | 1.7 | 5.9 | -10.4 | +10.8 | -31.6 | +0.0 | -10.5 | -1.3 | -5.0 (f) | -13.4 (f) | -4.7 (f) |
| F-free | Y_peak | 4.9 | 8.9 | -9.6 | -9.2 | -18.0 | +0.2 | -9.6 | -5.4 | -6.4 (f) | -13.2 (b-) | -8.2 (f) |
| F-free | BR_C | 27.7 | 51.6 | +63.8 | +100 | +553 | +1.0 | +66.9 | -5.0 | +49.4 (f) | +89.8 (f) | -21.9 (c-) |
| F-gold | BR_G | 14.9 | 23.7 | +15.0 | -25.3 | +95.5 | -0.2 | +14.8 | -0.3 | +12.3 (e) | +23.9 (e) | +7.5 (e) |
| B4 | Y_bump | 2.6 | 9.2 | -2.4 | +27.0 | -31.4 | -0.0 | -2.6 | -1.3 | -3.3 (f) | -5.8 (e) | -4.0 (e) |
| peak SB | Y_peak | 8.8 | 3.7 | -0.6 | +9.9 | -7.8 | -0.1 | -0.5 | -0.2 | +6.5 (f) | +5.5 (f) | +6.6 (f) |

- sigma_stat, sigma_MC and MC shift come from RF in scenario (a), as a fraction of the true value. The MC shift is the
  mean of the 100 replica fits of the Asimov data minus the true value.
- The true values in RF (a): a_eta 0.4, Y_bump 21 904, Y_peak 4 866 rows.
- Statistical precision of the biases: 0.1 % (BA, AB, AA, BB); 0.4-1 % (RF, RB: 100 independent template replicas).

### 2.2 Coverage (role error; full table in wp4_V1.txt, Fig. V1 coverage panel)
- **Method (AA/BB, sandwich error)**, scenarios a-d, g, h:
  - F-main 0.60-0.74;
  - F-gold 0.62-0.69;
  - B4 0.65-0.74;
  - peak SB 0.64-0.73;
  - F-free alpha_2 0.58 (BB d+) and 0.78-0.79 (h).
- **RF (sandwich (+) sigma_MC):**
  - F-main a_eta 0.33-0.49 and Y_bump 0.02-0.21. The error is right (pull width 1.01-1.03); the bias of -1.1 to
    -2.0 sigma is what fails;
  - B4 0.65-0.69;
  - peak SB 0.64-0.74;
  - F-gold 0.46 (h) to 0.65.
- **RB (bias-corrected):**
  - F-main a_eta 0.51-0.59 and Y_bump 0.51-0.56, with pull widths 1.50 and 1.46. The child RMS underestimates the MC
    error of the corrected estimator (bias correction adds variance: RMS of the corrected a_eta 0.053, uncorrected 0.045);
  - Y_peak 0.60-0.67; rho_d 0.63-0.68;
  - B4 0.62-0.67; peak SB 0.57-0.68;
  - F-gold 0.71-0.78.
- **Null scenario (g)**: mean pulls 0.01-0.10 in absolute value (AA, BB) and at most 0.23 (RF, RB) for every gated
  quantity. The null criterion passes everywhere except the SPEC configurations (BA/AB pulls up to 13).

### 2.3 Which template noise does what (docs/latex/tables_all/wp4_mcsplit.txt, Asimov of scenario (a), 100 replicas)

| part taken from a replica | F-main a_eta | F-main Y_bump | F-main Y_peak | rho_d | BR_G | B4 | peak SB |
|---|---|---|---|---|---|---|---|
| raw signal templates only | -10.2 (rms 3.4) | -11.5 (1.6) | -7.2 (2.0) | +8.0 (5.9) | +14.1 (6.2) | -1.4 (0.6) | -0.5 (3.3) |
| smoothed background only | +0.7 (9.9) | +0.2 (6.0) | +0.8 (9.3) | +1.4 (16.6) | +0.9 (18.7) | -0.9 (9.1) | -0.0 (1.0) |
| both | -10.1 (8.8) | -11.5 (5.6) | -6.6 (9.2) | +9.3 (19.4) | +15.0 (23.7) | -2.4 (9.2) | -0.6 (3.7) |

In short: the raw signal templates cause the bias, through attenuation (noisy regressors); the smoothed background
templates cause the spread. The side-band estimators do not scale a signal template, so they are not attenuated.

### 2.4 Other numbers
- **Fit quality in the toys** (mean chi2/ndf, F-main / stage A):
  - AA/BB 0.99-1.00 / 0.97-0.98;
  - RF 2.83 / 2.10 (a), 1.98 / 2.02 (g), 6.14 / 2.52 (h);
  - BA 8.8 / 4.7.
- **Correlations**: corr(rho_d, kappa) is 0.26 mean and 0.37 max in (a), 0.60 max in (h): the 0.9 fall-back is not
  reached. corr(a_eta, b) is 0.15-0.35, below the 0.7 limit of WP5.
- **Fisher scoring**: 17-26 iterations on average, at most 54; 0 failed fits.
- **Scenario (i)** (informational: true CB = per-cell geometric mean of the band):
  - a_eta +4.1/+4.8 % (method), -7.8 % (RF);
  - Y_bump +1.6/+1.8 %, B4 +0.5/+1.5 %, peak SB +0.1/-1.5 %.
  The fits absorb most of it through k_E1. It belongs in S1 (WP0 caveat 1).
- **Mismodelling scenarios in the method configurations** (the inputs for S3 and S6):
  - (e) omega: F-main Y_bump -3.7/-2.7 %, Y_peak +3.6/+4.1 %, rho_d +14 %; B4 -2.3/-2.0 %; peak SB -0.4/+0.1 %;
  - (f) mass transfer: F-main Y_bump -6.0/-4.7 %; B4 -1.6/-3.3 %; peak SB +5.6/+6.6 %; BR_G +2.7/+11.0 %.

---------------------------------------------------------------------------------------------------
## 3. Figures and tables (draft captions)

- **Fig. V1** `docs/latex/figs_all/wp4_V1.pdf` (.png). *Stress closure on simulation: relative bias of every gated quantity
  in the scenarios (a)-(h).*
  - Toys at L = 6079.2 nb^-1 (raw rows, isBest >= 0, selection S0); 500 per point.
  - Truth: simulated non-eta background (smoothed exp(pol3)) plus the raw eta classes plus the CB from the data like-sign
    pairs; toy N and the like-sign counts Poisson-fluctuated.
  - Open squares: SPEC test (truth half B, templates half A); open triangles: swap.
  - Full and open circles: truth and templates from the same half (method only).
  - Red squares: truth = full sample, templates = bootstrap replicas (the MC noise of the real analysis).
  - Magenta diamonds: the same with the bootstrap bias correction.
  - Green band: the SPEC tolerance (5 % a_eta, Y_bump; 10 % Y_peak; 20 % rho_d, BR).
  - Scenario g (a_eta = 0): bias in per cent of the true value of (a). Triangles at the frame edge are outside the range.
- **Fig. V1 (coverage)** `docs/latex/figs_all/wp4_V1_coverage.pdf`. *Fraction of toys whose error interval contains the truth.*
  - Errors: the fitter's error (SPEC), the sandwich error (method), and the sandwich error (+) sigma_MC from the
    template replicas (realistic MC).
  - Required band 0.60-0.76. In (e) and (f) the model does not contain the truth and coverage is not required.
- **Fig. V1 (pulls)** `docs/latex/figs_all/wp4_V1_pulls.pdf`. *Pull distributions of scenario (a).*
  - Method (blue), realistic MC (red), realistic MC with bias correction (magenta); dashed: N(0, 1) for 500 toys.
- **Table V1** `docs/latex/tables_all/wp4_V1.tex` (.txt with the per-scenario blocks of all configurations). Columns:
  - sigma_stat, sigma_MC and the replica mean shift (realistic MC, scenario a);
  - the bias in (a) for every configuration;
  - the largest bias over (a)-(h);
  - the verdicts. "fail (m)" = within 2 sigma of the limit.
- **Table** `docs/latex/tables_all/wp4_mcsplit.txt`: the split of the MC effects into signal and background templates
  (Sect. 2.3).
- `/hdd3/przygoda/HADES/FAT_out/all_report/wp4/wp4_summary.txt`: every number of every (configuration, scenario,
  quantity), including the Asimov value, the median errors of the three kinds, the pulls and the coverages.

---------------------------------------------------------------------------------------------------
## 4. Acceptance criteria (SPEC WP4) and fall-backs

| criterion | F-main | F-free | F-gold | B4 | peak SB |
|---|---|---|---|---|---|
| bias, SPEC (BA) | FAIL (+32 % a_eta) | FAIL | FAIL (-25 %) | FAIL (+27 %) | FAIL (+14 % in b-) |
| bias, method (AA, BB) | FAIL: (f) Y_bump -6.0 % (AA) | FAIL: BR_C +23/+49 % (e, f) | PASS (worst +12.3 %) | PASS (worst -3.3 %) | PASS (worst +6.5 %) |
| bias, realistic (RF) | FAIL (-10 to -15 %) | FAIL | FAIL (+23.9 % e, +22.4 % f) | FAIL marginal (-5.8 % e) | PASS (worst +5.5 %) |
| coverage 0.60-0.76 (role error) | FAIL in RF/RB | FAIL | FAIL in RF (0.46-0.54); RB 0.77/0.78 | PASS | PASS (RB h 0.57) |
| null (g): pull < 0.5 | PASS (method, RF, RB) | PASS | n/a | PASS | PASS |
| failed fits <= 1 % | PASS (0) | PASS (0) | PASS (3 of 500 in AB g) | PASS (0) | PASS (0) |
| **gate (method + RF)** | **FAIL** | **FAIL** | **FAIL** | **FAIL marginal** | **PASS** |
| with bias correction (method + RB) | FAIL (coverage) | FAIL | FAIL marginal | PASS | FAIL marginal |

**S10 entries** (|bias| of scenario (a)), recommended from RF (the estimator as specified):
- peak SB 0.5 %; B4 2.6 %;
- F-main a_eta 10.1 %, Y_bump 11.6 %, Y_peak 6.6 %, rho_d 9.7 %;
- BR_G 14.8 %.
With the bias correction (RB): 0.2 %, 1.3 %, 1.7 %, 1.8 %, 2.4 %, 2.5 %, 0.3 %. The SPEC (BA) values (32 % etc.) are
template-half MC noise and must not be used as S10.

**Fall-backs (SPEC Sect. 6) that fired:**
1. "An estimator fails WP4 -> cross-check only. If F-main fails, B4 + the peak side bands are the main result."
   - FIRED: F-main, F-free and F-gold are cross-checks.
   - Y_peak: the peak side bands (PASS).
   - Y_bump: B4, which fails marginally as specified (-5.8 % in the omega scenario against 5 %) and passes when its
     F-main inputs are bias-corrected with the replicas. Analyst's call.
2. "MC error > data statistical error -> Simulation-limited; request a larger std production." FIRED for F-main
   (3.0-3.3x), F-free, F-gold (1.6x) and B4 (3.6x); not for the peak side bands (0.42x).
3. "direct-decay bias > 20 % -> quote Y_peak only; BR from F-gold." Reached only in the stress scenario (e) of RF
   (rho_d +22.6 %); in (a) it is +9.7 %. Since F-gold also fails, no BR estimator is validated as specified.
4. Not fired:
   - corr(rho_d, kappa) > 0.9 (max 0.60);
   - |rho(a_eta, b)| > 0.7 (max 0.35).
5. Predicted, for WP3/WP5: the chi2/ndf fall-backs ("> 1.5") will fire from template MC noise alone (2.1 stage A,
   2.8 F-main, Sect. 2.4).

---------------------------------------------------------------------------------------------------
## 5. Caveats

1. **Why the SPEC configuration is not used as the gate (flagged deviation).**
   - Truth from half B and templates from half A make the template-to-truth MC difference twice the real one in sigma
     (variance 4 var_full against var_full).
   - It is a single realisation that does not average out over toys. In this simulation it is +-30 % on a_eta and B4.
   - So the SPEC bias criterion measures the MC statistics of one half-split, not the estimator.
   - The closest sound test with the same intent is RF: independent-looking templates through the bootstrap, at the real
     MC level, and the MC error that the analysis will quote. AA/BB isolate the method.
   - The SPEC result is reported in full (verdict.txt rows "BA").
2. **The bootstrap bias correction is not in the spec.** It is tested because it is free once WP5 runs its 100 replica
   fits.
   - It removes the attenuation: a_eta -10.1 % -> -1.7 %, Y_bump -11.6 % -> -1.8 %, BR_G +14.8 % -> -0.3 %.
   - But it inflates the MC error. The replica RMS then underestimates it by about 1.5 for a_eta and Y_bump, and by
     1.1-1.2 for Y_peak and rho_d.
   - If adopted for F-main, the MC error must be scaled by these factors, or the signal templates made less noisy
     (factorised or smoothed templates, S9), and WP4 rerun.
3. **Implementation choices that WP5/WP6 should copy** (or include `plots/all_report/wp4_lib.h`) so that this verdict
   applies to their numbers:
   - one-scale smoothing likelihood;
   - parent shapes for sparse histograms;
   - eta subtracted in U for B4;
   - unweighted C4 + C5 leakage for the side bands;
   - kappa_B and a kappa (C3 + C4) fixed from F-main in F-gold;
   - G_eta x alpha_2 g in F-free and F-gold;
   - stage A with eta fixed to F-main;
   - sandwich errors.
4. **No event-group clustering in the toys.** Toy rows are independent Poisson. The data event-group error is about
   1 % larger (P: 239.1 vs sqrt(55 947) = 236.5); the data fits should use the sandwich with event-group variances (D12).
5. **CB fluctuations** are Poisson per (2.5 MeV, slice, class) cell, not per event group. The fits use the data replica
   variance per bin, which ignores bin-to-bin CB correlations. The toys show this is harmless (B4 and side-band pulls
   0.94-1.06).
6. **Toy statistics.**
   - The coverage has a binomial error of 0.021.
   - RF and RB use 100 template replicas, so their bias error contains sigma_MC/10 (0.4-1 %).
   - The double-bootstrap parents give a slightly larger sigma_MC for a_eta than the WP0 replicas (11.9 % vs 9.7 % over
     100; 2 sigma).
7. **Scenario details (flagged).**
   - (e) keeps the tilt -0.15 of (a);
   - (f) applies s + 0.0009 and r = 0.89 together, and the yield truth is the content of the windows (resolution-dependent);
   - the toy totals in 0.40-0.70 are 0.86 (BA) to 0.93 (RF) of the data count. Scenario (a) fixes the background
     factor at 1.35.
8. **X17**: not run here. The X17 fit belongs to the x17ab package (x17_model.h, still changing while WP4 ran). Its X-b
   validation (1500 toys + CB, slope, coverage) is the WP4 equivalent for X17.

---------------------------------------------------------------------------------------------------
## 6. Files and exact commands (from /home/przygoda/HADES/FAT/FAT; ACLiC build dir aclic_wp4)

Code (plots/all_report/):
- `wp4_grid.h` (toy cells);
- `wp4_lib.h` (smoothing, CB, fitter, estimators, scenarios);
- `wp4_fill.C` (templates per half);
- `wp4_fill_boot.C` (double bootstrap);
- `wp4_toys.C` (toys);
- `wp4_eval.C` (verdict, tables, figures);
- `wp4_mcsplit.C` (Sect. 2.3).

Outputs (/hdd3/przygoda/HADES/FAT_out/all_report/wp4/):
- wp4_templates.root, wp4_fill_checks.txt;
- wp4_boot.root (100 parents x 20 children);
- toys/toys_{BA,AA,AB,BB,RF}.root, toys/toys_RB_m*.root (+ .log);
- wp4_summary.txt, verdict.txt.

    B='gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_wp4", kTRUE)'
    root -l -b -q -e "$B" 'plots/all_report/wp4_fill.C+O'                       # 30 s, 8 GB
    root -l -b -q -e "$B" 'plots/all_report/wp4_fill_boot.C+O(100, 20)'         # 4.5 min, 8 GB
    root -l -b -q -e "$B" -e '.L plots/all_report/wp4_toys.C+O'                 # compile ONCE before parallel runs
    for c in 0 1 2 3 4; do root -l -b -q -e "$B" "plots/all_report/wp4_toys.C+O($c, 500, 4095)" & done; wait   # 15-20 s each
    for m in 7 56 448 3584; do root -l -b -q -e "$B" "plots/all_report/wp4_toys.C+O(5, 500, $m, \"_m$m\")" & done; wait  # 80 s each
    root -l -b -q -e "$B" 'plots/all_report/wp4_eval.C+O'
    root -l -b -q -e "$B" 'plots/all_report/wp4_mcsplit.C+O'

Configurations: 0 BA, 1 AA, 2 AB, 3 BB, 4 RF, 5 RB.
Scenario bit mask: a = 1, b- = 2, ..., h = 1024, i = 2048.
