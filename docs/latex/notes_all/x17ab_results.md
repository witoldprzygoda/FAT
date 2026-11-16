# WP9 X-a (X17 model) and X-b (validation) -- results note (prefix x17ab)

Date 2026-11-14. Branch pp45_pippimepem (worktree /home/przygoda/HADES/FAT/FAT). Nothing committed, no git state changed.
ScanEngine.h, scan_*.C, the WP0 files and every other agent's file untouched; FAT_simwt only read.

**Blinding (D19).** No data row N was read anywhere in this package: not in X2D, not in the side bands. The only data input
is the like-sign CB of `hists_v1.root` (`data/S0/all/CB_X2D`, `CBrep_X2D`, `Npp_X2D`, `Nmm_X2D`), which is not blinded.
The sealed file `sealed/hists_v1_x17_sealed.root` was never opened. The sealed key of the SIMULATION blind test
(`.../blind/sealed/x17_key.root`, unblinded on 2026-11-13 by the simulation study) was read to get the injected BR of the 1500
toys. It was read only after all fits of a configuration were done, and once to select the BR = 0 toys for the parent check.

**Main deliverables.**
- `/home/przygoda/HADES/FAT/FAT/plots/all_report/x17_model.h`: model, likelihood, CLs, builders and the line of any mass. It is
  reusable by X-c and X-d.
- Macros `/home/przygoda/HADES/FAT/FAT/plots/all_report/x17ab_{line,model,toys,summary,parent,toycls,clsfig}.C` and `x17ab_tex.py`.
- Large outputs in `/hdd3/przygoda/HADES/FAT_out/all_report/x17ab/`.

----------------------------------------------------------------------------------------------------------------------
## 1. What was done, step by step (plain words)

### 1.1 The model (X-a), `x17_model.h`
The fit uses a 2D histogram:
- M4 = M(pi+pi-e+e-), 0.5225-0.610 GeV/c2, 35 bins;
- M(ee) 2-62 MeV/c2, 30 bins (layout X2D of WP0, selection S0).

In each cell the expected number of rows is the sum of five parts:

    mu_ij = k_e CB_ij                                                  (1) combinatorial background, from like-sign data
          + n_j B_c(i) [1 + (b_c + b'_c u_j) t_i + (c_c + c'_c u_j) t_i^2]   (2) other background
          + g_c(i) [N5 D_j + N4 C_j (1 + sC thC T_j)]                     (3) eta peak
          + r N5 R_X (1 + sX thX) gX(i) X_j                               (4) X17 line
          + a_eta [E23_ij + g Geta_ij]                                    (5) eta-related non-peak terms

Here t = (M4 - 0.55)/0.1 and u = (M(ee) - 30 MeV)/30 MeV. The terms in plain words:

1. **CB.** The measured CB (spread rule of WP0) times a factor k.
   - E2 (M(ee) >= 10 MeV): k = 1 +- 0.05.
   - E1: k = 1 +- 0.5 (deviation 1 below).
2. **Other background.**
   - n_j is free in every M(ee) bin, so no simulated background is used per bin (D18).
   - B_c is the M4 shape of the simulated G_D + G_C, summed over the rows of a class (E1, E2) and smoothed with exp(pol2).
   - A quadratic correction per class is free in the fit; its tilt and curvature may change linearly with M(ee).
3. **eta peak.**
   - M4 shape g_c: the direct decays C5 (V3) of the class.
   - D: M(ee) distribution of the direct decays (C5, V3).
   - C: M(ee) distribution of the conversions (C4, raw).
   - D and C are fractions of all rows of the M4 range, with the M(ee) flows in the total. N5 and N4 are therefore totals over
     the full M(ee) range, and both are free.
   - The conversion tail (M(ee) >= 12 MeV) has a nuisance sC = 0.194, from the 100 simulation replicas of C4.
4. **X17 line.**
   - X = the M(ee) line of the simulated X17 (C6, 17.0 MeV/c2); gX = its M4 shape.
   - R_X = eps_X / eps_5 = 1.3455: efficiency of eta -> pi+pi-X17(ee) over that of eta -> pi+pi-e+e-, for the same selection
     and range.
   - POI r = BR(eta -> pi+pi-X17) BR(X17 -> ee) / BR(eta -> pi+pi-e+e-), so BR_X BR_ee = r x 2.68e-4.
5. **eta-related terms.**
   - C2 + C3 and G_eta from the simulation, written as (row sum) x (smoothed M4 shape of the class).
   - a_eta has a Gaussian constraint: WP5 value; placeholder 0.40 +- 0.10 here.
   - g = 1 +- 0.30.

**Likelihood and limit.**
- Poisson likelihood of all 1050 cells (masks possible), with Gaussian penalty terms.
- Minimisation: Newton steps with the expected curvature (Fisher scoring), with step halving.
- Test statistic q~_r (Cowan et al. 2011), with r_hat >= 0 in the reference fit.
- 90 % CLs, computed two ways:
  - asymptotic: sigma from the background-only Asimov set;
  - toys: pseudo-data from the conditional fits at r and at 0, with the constraint centres drawn as well.
- Expected bands from the Asimov set.
- Helper for X-c, `FillMaskedRows`: it predicts n_j of masked M(ee) rows from the side-band rows before an Asimov set is built.

### 1.2 The X17 line for any mass (`x17ab_line.C`)
1. Only one X17 mass was simulated (17.0 MeV/c2). For a mass m the line is built from direct decays (C5) whose TRUE M(ee) is
   within +-0.5 MeV/c2 of m. Each such row is treated in two steps:
   - it is shifted in reconstructed M(ee) by (m - q_true);
   - it is reweighted from the generator kinematics to X17 kinematics at fixed q:
     - M(pi+pi-) from three-body phase space;
     - isotropic e+ angle;
     - isotropic pi+ angle (weight capped at |cos| = 0.97).
2. The weights average to 1 over the generated decays at fixed q. The slice efficiency is therefore (accepted weight) /
   (generated decays in the slice), and the generated decays come from the generator density g(q) H(q) of
   `sim_eta_physics.h`.
3. Validation on C6 at 17 MeV/c2 (Sect. 2.1): the slice line and efficiency agree with the real X17 sample.
4. The slices are statistically limited (about 1200 effective rows at 17 MeV). The model therefore uses:
   - the C6 line, morphed to mass m: its mean and width follow the slices' truncated mean and RMS, smoothed by pol2 in m;
   - R_X(m) = R_X(17, C6) x eps(m)/eps(17), with eps from a pol2 fit of the slice efficiencies;
   - at exactly 17.0 the C6 line itself.
5. Variants of the reweighting:
   - DL: M(pi pi) "direct-like";
   - C2: e+ decay angle 1 + cos^2;
   - NP: generator pion angle kept;
   - raw slice line instead of the morphed C6 line.

### 1.3 Projection and model choices (`x17ab_model.C`)
1. The model is built on the data-scale X2D, with sim/S0__nom and the data CB.
2. A simulation Asimov set ("projection", pre-data) is made:
   - n_j = 1.35 x simulated G_D + G_C per M(ee) bin;
   - peak and eta terms at a_eta = 0.40;
   - data CB.
3. The expected limit is computed at 17 MeV and for m = 10-40 MeV in 1 MeV steps.
4. Every model choice of the WP9 nuisance list is tested the same way: the nominal Asimov set is refitted with the changed
   model, and the bias r_hat and the limit are recorded. Sect. 2.3 lists the results.

### 1.4 Validation on toys (X-b, `x17ab_toys.C`, `x17ab_summary.C`, `x17ab_parent.C`)
1. **The 1500 toys of the simulation study, with a CB-like component added** (the spec's test).
   - The toys are Poisson fluctuations of half B. They hold T1 + T2 plus X17 at BR = 0, 1, 3, 10, 30 x 1e-6, 300 toys each,
     at L = 6080.5 nb^-1.
   - Toy grid: simulated M4 0.5230-0.6080 (34 bins) x M(ee) 2-62 MeV (30 bins).
   - Templates come from half A only.
   - CB-like component: Poisson(CB_data) rows are added to each toy. The fit uses as CB a Poisson event-group replica of the
     like-sign rows (`CBrep_X2D`, replica t mod 200), i.e. N++ and N-- fluctuated.
   - Variants: no CB; CB injected as a per-row geometric mean but fitted with the spread rule; background pol2 / three classes /
     fixed / tilt-only; no conversion-tail nuisance; line from the C5 slice; CB without fluctuation; X17-yield nuisance 10 %.
2. **"Smooth-parent" toys (own, added).**
   - Truth = the model itself, with the parameters of a fit to the half-A Asimov set. N5 and N4 are set to the half-A
     expectations; a_eta = g = 1.
   - X17 is added as r N5 R_X gX X; Poisson rows plus Poisson(CB) are generated, 300 per BR.
   - Reason: the simulation-study toys share one parent, whose heavy simulated rows (one std row = 55 expected rows in one cell)
     give a Pearson chi2/ndf of about 29 for any smooth model. The coverage of a smooth-background method can only be tested
     on toys without that structure.
   - Stress tests:
     - the truth uses THREE background classes (E1, 10-30, >= 30 MeV, each with its own simulated shape) and is fitted with
       two (cfg 13, 14);
     - cross-checks with three-class fits (cfg 16, 17) and without the M(ee) dependence (cfg 15).
3. **Parent check.** Two sets are fitted with the nominal model:
   - the mean of the 300 BR = 0 toys (the half-B parent);
   - the half-A Asimov set.
4. **Asymptotic against toy CLs** (`x17ab_toycls.C`, `x17ab_clsfig.C`). Eight data sets, 8 values of r, 1000 toys per value
   and hypothesis. The toy limit comes from a weighted quadratic fit of ln CLs(r).

### 1.5 How the nominal background form was chosen (all on simulation, before any data)
| background form | test | outcome |
|---|---|---|
| tilt only, B (1 + b t) (the literal spec form) | smooth-parent toys whose truth has a background shape that changes with M(ee) (cfg 13, with curvature but no M(ee) dependence) | offset -8.4e-6, coverage 0.39 / 0.33 at BR 1e-5 / 3e-5: FAIL |
| B (1 + (b + b'u) t + (c + c'u) t^2), two classes | the same truth (cfg 14) | offset -1.0e-6, coverage 0.94 / 0.91: PASS |
| the same form | 2-class truth (cfg 10) | PASS; cost in the projection limit: +13 % against tilt only |

The smoothing order of B matters as well:
- An exp(pol3) shape has 1 %-level wiggles under the peak. The peak itself is only 1-4 % of the background per M(ee) bin.
- With exp(pol3), the three-class variant of the projection moved r by -1.8e-5 and N5 by +51 %.
- With exp(pol2), the same variant moves r by -3.0e-6 and N5 by +22 %.
- exp(pol2) is therefore the nominal smoothing.

----------------------------------------------------------------------------------------------------------------------
## 2. Numbers

### 2.1 X17 line and efficiencies (S0, M4' 0.5225-0.610, all M(ee); simulation, mode A, L_sim = 221.05 nb^-1)
**C6 (X17 at 17.0 MeV/c2).**
- sigma_vis = 6257.5 nb per unit BR_X BR_ee; eps_X(17) = 5.867e-3.
- Line: mean 17.33 MeV/c2, RMS 4.54. Truncated (17 +- 12): 17.24 / 4.19.
- Double Gaussian: mu 17.26, sigma 3.69 / 5.53, core fraction 0.65. The simulation study quoted 17.25 / 4.04.

**Efficiency ratio.** C5 (V3): sigma_vis = 1.2464 nb, eps_5 = 4.361e-3, so R_X(17) = 1.3455.

**Validation of the slice method at 17 MeV/c2** (C5 rows with q in 16.5-17.5, reweighted and shifted):
| quantity | value |
|---|---|
| eps_slice / eps_C6 | 1.043 +- 0.030 (1.4 sigma) |
| line shape chi2 (1 MeV bins, 5-32 MeV) | 26.2 / 27 |
| line mean / RMS, slice vs C6 | 17.29 / 4.48 vs 17.33 / 4.54 |
| M4' shape chi2, slice vs C6 | 18.4 / 24 |
| M4' shape chi2, C6 vs all C5 | 238 / 35 |

The last line means that the X17 M4 shape differs from that of the direct decays. The model therefore takes gX from C6.

**R_X(m)** (model curve, pol2 fit chi2 57/72 over 8-45 MeV):

| m [MeV/c2] | 10 | 15 | 17 | 20 | 25 | 30 | 35 | 40 |
|---|---|---|---|---|---|---|---|---|
| R_X | 1.275 | 1.328 | 1.346 | 1.368 | 1.394 | 1.406 | 1.406 | 1.392 |

Variants at 17 MeV/c2:
- direct-like M(pi pi): -1.9 %;
- generator pion angle: -3.9 %;
- e+ decay angle 1 + cos^2: **-16.5 %**.

Line width: the truncated RMS grows from 3.8 (m = 8) to 4.3-4.4 MeV/c2 (m > 25). The truncated mean - m falls from +1.3 to
-0.4 MeV/c2.

### 2.2 Model inputs (data-scale X2D; raw rows at L = 6079.2 nb^-1 x 63/64)
- N5 (C5, a_eta = 1) = 7458.8 rows; N4 (C4) = 6656.3 rows.
- X17 rows per unit BR at a_eta = 1: s_X2D = N5 R_X / 2.68e-4 = **3.745e7**. At BR 1e-5 that is 374.5 rows.
  - This is the s of normalisation (2) for the 2D fit: N_X / (a_eta s_X2D).
  - The spec's s1 = 2.127e7 belongs to the X17 counting window.
- Fractions in M(ee) 2-62 MeV: D 0.855, C 0.970, X 0.9996. In 14-26 MeV: D 0.189, C 0.035, X 0.752.
- Conversion-tail fraction (>= 12 MeV) of C4: 0.088 +- 0.017 (100 replicas), i.e. sC = 0.194.
- Data CB in X2D: 6266 rows (spread rule).

### 2.3 Pre-data projection (simulation Asimov set; a_eta 0.40; background 1.35 x SMASH; data CB; NOT a result)
**Expected 90 % CLs at 17 MeV/c2, BR_X BR_ee:**

| -2 sigma | -1 sigma | median | +1 sigma | +2 sigma |
|---|---|---|---|---|
| 1.21e-5 | 1.76e-5 | **2.83e-5** | 5.07e-5 | 1.07e-4 |

- The limit corresponds to 424 X17 rows.
- The X-c recipe (side-band fit with rows 12-28 MeV masked, then `FillMaskedRows`, then Asimov) reproduces the same numbers.
  Without `FillMaskedRows` it gives a spurious 1.44e-5, because the masked rows have no information.
- Statistical error of N5 in this fit: 21.5 %.

**Median versus mass:**

| m [MeV/c2] | 10 | 12 | 15 | 17 | 20 | 25 | 30 | 35 | 40 |
|---|---|---|---|---|---|---|---|---|---|
| median [1e-5] | 5.29 | 4.09 | 3.19 | 2.83 | 2.42 | 1.88 | 1.49 | 1.25 | 1.07 |

**Model choices on the projection** (bias r_hat x 2.68e-4, and UL/UL0; table `x17ab_nuis.txt`):

| model choice | bias [1e-6] | UL/UL0 |
|---|---|---|
| X17 line width +10 % / -10 % | 0 | 1.10 / 0.90 |
| M(ee) scale +-2 % | 0 | < 1.003 |
| line from the C5 slice | 0 | 0.96 |
| R_X: direct-like M(pi pi) | 0 | 1.02 |
| R_X: generator pion angle | 0 | 1.04 |
| R_X: e+ decay angle 1 + cos^2 | 0 | **1.20** |
| direct-decay templates V2 | +0.02 | 1.19 (eps_5 of V2 is 16 % higher) |
| direct-decay templates GS | +0.005 | 1.005 |
| peak width r = 0.89 / 1.10 | -0.19 / +0.21 | 1.03 / 0.98 |
| mass scale s +- 0.0009 | < 0.05 | <= 1.01 |
| conversion replicas (20) | r_hat RMS **2.1** | -- |
| three background classes | **-3.0** | 0.68; N5 +22 % |
| background template exp(pol3) | +1.0 | 0.75; N5 +26 % |
| background template exp(pol1) | 0.0 | 1.00 |
| background pol2 only | -0.6 | 0.84 |
| tilt only | 0 | 0.88 |
| curvature only | 0 | 0.88 |
| background shape fixed | 0 | 0.88 |
| background template raw (not smoothed) | fit unstable (heavy simulated rows) | not a variant |
| kappa_B = 1.3 | 0 | 1.00 |
| CB factorised | +0.2 | 0.99 |
| CB as per-row geometric mean | +0.7 | 1.01 |
| k_E1 free; k_E2 fixed | 0 | 1.00 |
| no conversion-tail nuisance | 0 | 1.00 |
| a_eta and g fixed | 0 | 0.99 |
| X17-yield nuisance 10 % | 0 | 1.01 |

**N5-bias nuisance.** The half spread of N5 over the background and CB variants is **10.9 %**.

### 2.4 Toys (BR in 1e-6; tables `x17ab_toys.txt`, `.tex`)
The two normalisations:
- (1) BR1 = r_hat x 2.68e-4;
- (2) BR2 = N_X / s = r_hat N5_hat / N5_true x 2.68e-4.

`cov` = fraction of toys whose UL (asymptotic CLs 90 %) is >= the injected BR.

| cfg | toys and fit | slope (1) | slope (2) | offset (2) | cov at 10 / 30 | cov (2), common offset removed | chi2/ndf |
|---|---|---|---|---|---|---|---|
| 0 | sim-study toys + CB-like, NOMINAL | 0.673 +- 0.008 | **1.038 +- 0.012** | -19.5 | 0.000 / 0.000 | **1.000 / 0.957** | 28.8 |
| 1 | sim toys, no CB | 0.692 | 1.039 | -17.6 | 0 / 0 | 1.000 / 0.957 | 30.7 |
| 2 | sim toys, CB geometric mean injected, spread rule fitted | 0.676 | 1.042 | -19.8 | 0 / 0 | 1.000 / 0.960 | 29.0 |
| 7 | sim toys, line from the C5 slice | 0.638 | 0.993 | -20.5 | 0 / 0 | 1.000 / 0.950 | 28.8 |
| 18 | sim toys, background without M(ee) dependence | 0.703 | 1.037 | -16.3 | 0 / 0 | 1.000 / 0.940 | 28.9 |
| 10 | smooth parent + CB, NOMINAL | **1.024 +- 0.015** | **1.015 +- 0.012** | -0.3 | **0.937 / 0.917** | 0.953 / 0.950 | 1.01 |
| 11 | smooth parent, CB mismatch | 1.010 | 1.007 | -0.3 | 0.930 / 0.907 | -- | 1.00 |
| 14 | 3-class truth, NOMINAL fit | 1.056 +- 0.016 | 0.991 +- 0.012 | -1.0 | 0.937 / 0.913 | -- | 1.18 |
| 13 | 3-class truth, fit without M(ee) dependence | 0.910 | 0.994 | -8.4 | 0.393 / 0.333 | -- | 1.24 |
| 15 / 16 / 17 | other shape cross-tests | 0.97-1.00 | 0.97-1.00 | +-0.3 | 0.93 / 0.87-0.90 | -- | 1.0-1.1 |

More detail on these toys:
- **Smooth parent, nominal (cfg 10).**
  - Pulls: mean -0.09, RMS 0.98.
  - The +-1 sigma interval covers 64-70 % of the toys.
  - Median UL at BR = 0: 9.2e-6 (a_eta = 1, simulation background).
- **Simulation-study toys:**
  - N5_hat = 11168 at BR = 0, against a true 7200 (+55 %).
  - Median UL at BR = 0: 1.9e-6, because of the negative offset.
- **Parent check** (`x17ab_parent.txt`): both sets are fitted with the nominal model.

  | set | N5_hat | BR2 | Pearson chi2/ndf |
  |---|---|---|---|
  | mean of the 300 BR = 0 toys (half B) | 10876 (+51 %) | -17.8 | 29.7 (even for the mean of 300 toys) |
  | half-A Asimov set | 3005 (-58 %) | +18.4 | -- |

  The two parents give opposite offsets of about 18e-6. This is the size of the simulation-background noise projected on X17.
  The simulation study found 11.8e-6 to 23e-6 for counting.

### 2.5 Asymptotic against toy CLs (`x17ab_toycls.txt`)
| set | data set | UL asymptotic [1e-6] | UL toys [1e-6] | toys/asymptotic |
|---|---|---|---|---|
| 0 | smooth-parent Asimov | 9.22 | 9.60 | 1.042 |
| 1 | smooth toy, BR 0 | 9.58 | 9.70 | 1.013 |
| 2 | smooth toy, BR 0 | 6.54 | 6.51 | 0.995 |
| 3 | smooth toy, BR 1e-5 | 14.73 | 14.80 | 1.005 |
| 4 | smooth toy, BR 1e-5 | 9.49 | 9.39 | 0.990 |
| 6 | sim-study toy, BR 0 | 1.97 | 1.90 | 0.967 |
| 7 | data-scale projection Asimov | 28.29 | 29.58 | 1.046 |
| 5 | sim-study toy, BR 0 | 1.58 | undefined | -- |

Set 5: its r_hat lies 4.6 sigma below 0 (parent structure). The observed q~ is then above every toy under both hypotheses, so
CL_b = 0 and the toy CLs is undefined.

----------------------------------------------------------------------------------------------------------------------
## 3. Figures and tables (draft captions)

**`docs/latex/figs_all/x17ab_line_valid.pdf`**
- Left panel: reconstructed M(ee) of three samples, unit area per 1 MeV/c2:
  - the simulated X17 (C6, m_X = 17.0 MeV/c2, black);
  - the direct decays eta -> pi+pi-e+e- with true M(ee) 16.5-17.5 MeV/c2, shifted to 17.0 and reweighted to X17 kinematics
    (red: phase-space M(pi+pi-), isotropic e+ and pi+ angles);
  - the same with direct-like M(pi+pi-) (blue dashed).
- Right panel: M'(pi+pi-e+e-) of C6, of all direct decays (C5, V3) and of the slice.
- Selection S0; M4' 0.5225-0.610 GeV/c2; SMASH eta set, mode A.
- Simulation only, no data.

**`x17ab_line_eff.pdf`**
- Left panel: R_X = eps_X/eps_5 versus m_X.
  - Points: slices of +-0.5 MeV/c2 for four kinematic assumptions (phase space; direct-like M(pi pi); e+ angle 1 + cos^2;
    generator pion angle).
  - Star: the simulated X17.
  - Line: the model, a pol2 in m anchored on C6 at 17 MeV/c2.
- Right panel: truncated (m +- 12 MeV/c2) mean - m_X and RMS of the slice lines, with the pol2 curves that morph the C6 line.
- S0, X2D M4 range.

**`x17ab_model_mee.pdf`**
- M(ee) shapes of the 2D model: fraction of all rows of a component per 2 MeV/c2 bin (M(ee) flows in the total).
- Components:
  - direct decays D (C5, V3);
  - conversions C (C4; tail +-19 %);
  - X17 at 17 MeV/c2 (C6, R_X = 1.345);
  - X17 at 12 and 30 MeV/c2 (morphed).
- S0, M4' 0.5225-0.610 GeV/c2. Simulation templates; L = 6079.2 nb^-1.

**`x17ab_model_m4.pdf`**
- M4' shapes of the 2D model, unit area:
  - eta peak in E1 and E2 (C5);
  - X17 (C6);
  - background templates (simulated G_D + G_C, exp(pol2)-smoothed) in E1 and E2;
  - data CB per class (like-sign, spread rule; the steps at 0.55 and 0.60 GeV/c2 are the 50 MeV band edges).

**`x17ab_proj_mass.pdf`**
- PRE-DATA projection of the expected 90 % CLs upper limit on BR(eta -> pi+pi-X17) BR(X17 -> e+e-) versus m_X, with the
  +-1 sigma and +-2 sigma bands.
- Inputs: simulation Asimov set (a_eta = 0.40, background 1.35 x SMASH, data CB); L = 6079.2 nb^-1; S0, X2D.
- Superseded by X-c and X-d.

**`x17ab_toys_lin.pdf`**
- Linearity of the X17 fit:
  - left: the 1500 toys of the simulation study plus a CB-like component;
  - right: smooth-parent toys plus CB.
- Fitted against injected BR_X BR_ee in the two normalisations, (1) r x 2.68e-4 and (2) N_X/s; bars = RMS of 300 toys.
- Nominal model.
- Left panel: the common offset and the N5 bias come from the heavy-row structure of the simulated parent (Sect. 2.4).

**`x17ab_toys_ul.pdf`**
- Left: 90 % CLs limits at injected BR = 0. Right: pulls (r_fit - r_true)/sigma_r.
- Black: simulation-study toys. Blue: smooth-parent toys.

**`x17ab_toycls.pdf`**
- CLs versus r for four data sets: toys (1000 per point and hypothesis, binomial errors) against the asymptotic formula (red).
- Dashed line: CLs = 0.10.

**Tables** (all in `docs/latex/tables_all/`):

| file | content |
|---|---|
| `x17ab_lines.txt` | line, eps, R_X per mass and variant; validation numbers |
| `x17ab_model.txt` | model inputs and per-M(ee)-bin D, C, X, CB, simulated background |
| `x17ab_nuis.txt` / `.tex` | model choices on the projection |
| `x17ab_proj_mass.txt` | projection versus mass |
| `x17ab_toys.txt` / `.tex` | all toy configurations |
| `x17ab_parent.txt` | parent check |
| `x17ab_toycls.txt` / `.tex` | asymptotic against toy CLs |

----------------------------------------------------------------------------------------------------------------------
## 4. Acceptance criteria and fall-backs

| criterion (spec WP9) | value | verdict |
|---|---|---|
| X-a: other masses from C5 slices, validated on C6 | eps ratio 1.043 +- 0.030; line chi2 26/27; mean and RMS within 0.04 / 0.07 MeV/c2 | PASS |
| X-b: slope 1 +- 0.05, simulation-study toys + CB, normalisation (2) | 1.038 +- 0.012 | PASS |
| X-b: the same, normalisation (1) | 0.673 +- 0.008 | FAIL: N5 +55 % from the parent structure |
| X-b: slope, smooth-parent toys | (1) 1.024 +- 0.015, (2) 1.015 +- 0.012 | PASS |
| X-b: slope, three-class stress test | (1) 1.056 +- 0.016, (2) 0.991 +- 0.012 | PASS (2); (1) at the edge |
| X-b: 90 % coverage >= 0.88, simulation-study toys, raw | 1.00 / 1.00 / 0.11 / 0.00 / 0.00 | FAIL: common offset -19.5e-6 |
| X-b: the same, common offset removed, normalisation (2) | 1.00 / 1.00 / 1.00 / 1.00 / 0.957 | PASS (diagnostic) |
| X-b: coverage, smooth-parent toys | 1.00 / 1.00 / 1.00 / 0.937 / 0.917 | PASS |
| X-b: coverage, three-class stress test | 0.937 / 0.913 | PASS |
| asymptotic against toy CLs within 10 % | 0.967-1.046 on 7 sets (1 undefined) | PASS: asymptotic stays primary |

**Fall-backs** (spec Sect. 6, X17 row):
- **"Pearson chi2/ndf of n_j > 1.5: scale the errors" fires on the simulation-study toys** (28.8). It is a property of the
  toy parent: the mean of 300 toys still gives 29.7. It was not applied, because scaling by sqrt(28.8) = 5.4 would hide the
  problem. For the data, X-c must evaluate it on the side bands.
- **"Slope off by > 5 %" fires for normalisation (1) on the simulation-study toys**, so normalisation (2) is primary for them.
  - On the smooth toys slope (1) passes.
  - The decision for the data is X-c's.
  - My recommendation: carry (1) and (2) in parallel. (1) depends on N5 from the X2D fit; its statistical error is 21 % and
    its shape-variant spread 11 %.
- **"N5 error > 50 %": not fired** (21.5 % statistical, 10.9 % bias spread on the projection).
- No other fall-back concerns X-a/X-b.

**Deviations from the spec** (the closest sound choice; flagged):

| # | deviation | reason |
|---|---|---|
| 1 | k_E1 = 1 +- 0.5 instead of free | With n_j free per M(ee) bin, k_E1 is degenerate with the E1 rows; free, it ran to -19 on a toy. On the projection, free and constrained give the same limit (1.001). |
| 2 | Background: B (exp(pol2)-smoothed simulation, per class) x a quadratic correction whose coefficients change linearly with M(ee) | The literal "WP5 shape x n_j" (tilt only) fails coverage when the shape changes with M(ee) (cfg 13). X-c can insert the WP5 shapes into Inputs::B and keep bkgMode 4. Cost: +13 % on the expected limit. |
| 3 | C2 + C3 and G_eta factorised and smoothed (D13 says signal templates raw) | Raw heavy rows put single-cell spikes into the model. On a toy, raw G_eta pulled a_eta to 0.09 and N5 to +64 %. |
| 4 | C = raw C4 plus a tail nuisance (sC 0.194 from the 100 replicas) | Interpretation of "replica-averaged". The replica spread of r_hat (2.1e-6 in BR on the projection) should be added as a systematic. |
| 5 | Line of other masses = C6 morphed (mean, width from the slices); R_X(m) a pol2 anchored on C6 | The raw slices are statistically limited. The raw slice is a variant (UL -4 % at 17 MeV). |
| 6 | Toy grid in simulated M4, at the old L = 6080.5 and without 63/64 | This is how the simulation study produced the toys; irrelevant for the validation. The data CB is mapped bin by bin (offset 1.5 MeV in M4; the CB is smooth). |
| 7 | Own smooth-parent toys added | The simulation-study toys cannot test coverage (Sect. 1.4). |

----------------------------------------------------------------------------------------------------------------------
## 5. Caveats (for X-c / X-d and the report)

1. **The peak is 1-4 % of the background per M(ee) bin.**
   - A background-shape error of about 0.5 % near m_eta moves N5 by tens of per cent, and with it normalisation (1).
   - X-c must check the background shape on the data side-band rows: Pearson chi2 per class and per row, and the residuals
     near m_eta.
   - It must also evaluate the variants of Sect. 2.3. Avoid cubic-or-higher template shapes.
2. **Masked rows.**
   - Before the Asimov set of X-c is built from a side-band fit, call `x17m::FillMaskedRows`. The masked n_j stay at their
     start values otherwise.
   - Reference per row: e.g. the simulated G_D + G_C row sums, given in `x17ab_model.txt`.
3. **The e+ decay angle of the X17 (isotropic or 1 + cos^2) changes R_X by -16.5 % and the limit by +20 %.**
   - This is a model dependence, not a measurement error.
   - Recommendation: quote the isotropic (nominal) limit and the 1 + cos^2 limit separately, and keep it out of sX.
   - Suggested sX for normalisation (1), in quadrature:
     - DL 1.9 %;
     - pion angle 3.9 %;
     - slice/C6 validation 4.3 %;
     - N5 bias 10.9 %;
     - BR_dir 4.1 %;
     - total 13 %.
   - With the nuisance in the likelihood, the limit changes by about +1 % per 10 % of sX.
4. **The direct-decay model (V2 against V3) changes eps_5 by 16 % and the limit by 19 %** (normalisation (1)). This belongs to
   the "direct shape" nuisance of WP9.
5. **The two CB constructions agree in the X17 rows** (bias < 0.7e-6). The spread-rule CB has steps at the band edges 0.55 and
   0.60 GeV/c2 inside the X2D range (Fig. model_m4). No effect was seen.
6. **Toy CLs can be undefined** when the observed q~ is above every toy (seen on one simulation-study toy with r_hat at
   -4.6 sigma). If this happens in data, report it, check the goodness of fit, and use the asymptotic CLs.
7. **The projection is not a result.** a_eta = 0.40 and the 1.35 background scale are placeholders. Its limit (2.8e-5 at
   17 MeV) is weaker than the simulation study's counting value (5.7e-6 with B known, a_eta = 1), for three reasons:
   - a_eta = 0.4;
   - the background per M(ee) bin is profiled;
   - N5 is fitted.
8. **The generator q-density** (log-log interpolation of g(q)) integrates to 1.023, so the absolute slice efficiencies carry
   +-2 %. R_X is anchored on C6, so this matters only for the mass dependence (pol2 chi2 57/72).
9. **The sealed key of the simulation blind test was read** (that test was unblinded on 2026-11-13); the data sealed file
   never.

----------------------------------------------------------------------------------------------------------------------
## 6. Exact commands (from /home/przygoda/HADES/FAT/FAT; build dir /hdd3/przygoda/HADES/FAT_out/all_report/aclic_x17ab)

    B='gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_x17ab", kTRUE)'
    root -l -b -q -e "$B" 'plots/all_report/x17ab_line.C+O'                 # 2 min, 3 GB: x17_lines.root, x17ab_lines.txt, line figures
    root -l -b -q -e "$B" 'plots/all_report/x17ab_toys.C+O(99, 0, 0)'       # compile once (writes an empty toys_cfg99_0.root)
    for c in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18; do              # <= 6 in parallel, about 1 min each
      root -l -b -q -e "$B" "plots/all_report/x17ab_toys.C+O($c, 0, 1500)" > /hdd3/przygoda/HADES/FAT_out/all_report/x17ab/toys_cfg$c.log 2>&1
    done
    root -l -b -q -e "$B" 'plots/all_report/x17ab_summary.C+O'              # x17ab_toys.txt, toys_lin, toys_ul
    root -l -b -q -e "$B" 'plots/all_report/x17ab_parent.C+O'               # x17ab_parent.txt
    root -l -b -q -e "$B" 'plots/all_report/x17ab_model.C+O'                # 1.5 min: model, projection, nuisances, mass scan
    for w in 0 1 2 3 4 5 6 7; do                                             # 7 needs x17ab_model.root; 3-5 min each
      root -l -b -q -e "$B" "plots/all_report/x17ab_toycls.C+O($w, 1000)" > /hdd3/przygoda/HADES/FAT_out/all_report/x17ab/toycls_$w.log 2>&1
    done
    root -l -b -q -e "$B" 'plots/all_report/x17ab_clsfig.C+O'               # x17ab_toycls.txt and figure
    python3 plots/all_report/x17ab_tex.py                                    # .tex snippets

**Use of the model in X-c/X-d** (sketch):

    #include "x17_model.h"
    x17m::Inputs I;  x17m::HistOpt o;                  // sel S0, var nom, cbMode 0, nbcl 2, exp(pol2) B, smoothed eta terms
    x17m::BuildFromHists(f, I, o, alldefs::kNbToRaw);  // optional: replace I.B by the WP5 shapes (unit sum per class)
    x17m::LineFor(mX, I.G, I.X, I.RX);                 // any mass (17.0 = C6); variants: rxVar 1-3, scale, res, rawSlice
    x17m::Config C; C.kSig[0] = 0.5; C.kSig[1] = 0.05; C.a0 = <a_eta WP5>; C.aSig = <err>; C.gSig = 0.30;
    C.bkgMode = 4; C.cTailSig = 0.194; C.xSig = <sX>; C.use = <mask>;
    x17m::Model M(I, C);  x17m::Limit L = x17m::LimitAsym(M, n, x17m::GlobOf(C), 0.04);   // UL in r; BR = r x 2.68e-4
    // side bands (X-c): fit with C.use masking rows 12-28 MeV, x17m::FillMaskedRows(M, p, simRow), Asimov = M.Expected(p)
