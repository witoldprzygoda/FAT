# WP3 control regions, stage A (prefix wp3) -- results note

Date 2026-11-14. Branch pp45_pippimepem (worktree /home/przygoda/HADES/FAT/FAT). Nothing committed; no git state changed;
ScanEngine.h, scan_*.C and every WP0 file untouched; FAT_simwt only read (through the WP0 copies).

**Frozen output:** `/hdd3/przygoda/HADES/FAT_out/all_report/wp3/stageA.txt` (machine-readable; + `stageA.root`).
**Macros:** `plots/all_report/wp3_rows.C` (row-level parts), `plots/all_report/wp3_stageA.C` (fits, tests, tables, figures).

**Headline.** Four fall-backs fired. The frozen stage A is:
- **model** M1 = c_j (1 + b x), common tilt;
- **data** M(ee) >= 10 MeV only (E2);
- **control region** CR-U 0.575-0.700 only; CR-K is dropped.
- **Consequences:**
  - The bump background is **EXTRAPOLATED**: the tilt comes from CR-U only.
  - The templates are **simulation-limited**: the MC error is 4-5 x the data error on c_j.

Two deviations from the literal spec are flagged (Sect. 5):
- D11 and the chi2 <= 1.5 test use the chi2 **with** the MC covariance of the smoothed templates;
- E1 is dropped, but kappa_B stays free.

---------------------------------------------------------------------------------------------------
## 1. What was done, step by step (plain words)

1. **Inputs.**
   - The 10 F10 histograms (5 M(pi+pi-) slices x 2 M(ee) classes, M4 0.40-0.70), selection S0 (isBest >= 0), from hists_v1.root:
     - data N, and CB = 2 sqrt(N++ N--) by the spread rule;
     - the CB bin variance from the 200 CB replicas;
     - simulation (SMASH, mode A, w x V3), converted to expected raw data rows (nb x 5984.2).
   - From the data rows (wp3_rows.C), the event-group covariance of N over all 700 F10 cells:
     - S0 and pre;
     - grouping by adjacency;
     - the cell counts equal hists_v1 to 1e-12.
2. **Smoothing (D13).**
   - Per histogram, the non-eta templates G_D (= pi0 Dalitz + vector + other) and G_C (conversions) are replaced by exp(cubic
     polynomial in M4), integrated over each bin.
   - The fit is a Poisson-type likelihood. It keeps the integral exactly and starts at the first bin with simulated content (the
     kinematic threshold).
   - The 20 fits have chi2/ndf 0.6-1.4 (p 0.03-0.96).
   - The eta templates (C2-C5, G_eta) stay raw. Their bin variance enters the weights (BB-lite).
3. **Stage-A fit.**
   - Model, in every bin of the control region:

         N = k_e CB + c_j (1 + b x)[G_D + kappa_B G_C] + 0.4 x (all eta components at SMASH)

   - x = (M4 - 0.55)/0.1.
   - Poisson maximum likelihood, by iteratively re-weighted least squares.
   - Errors:
     - "stat": sandwich with the event-group covariance;
     - "MC": the RMS of the parameters over 100 simulation replicas, each smoothed again and refitted.
4. **Step 0, the spec as written** (E1 + E2, CR-U 0.565-0.700 + CR-K).
   - D11 with the data-only chi2 adds everything: quadratic, then per-slice tilts. chi2/ndf = 760/297.
   - **k_E1 = 13.19 +- 0.61.** The k_E1 term (CB shape in E1) absorbs an E1 excess:
     - E1 data/sim 1.30-1.43, against 1.12-1.29 in E2;
     - the fake-lepton template T3 (not in the model, D6) is 26.5 % of T12 in E1;
     - with the ghost-variant (a) templates k_E1 is still 8.1 +- 0.6.
   - -> Fall-back **"k_E1 outside 0.8-3: drop E1"** fired. All later fits use E2 only.
5. **The CR chain of A1/A2**, each configuration with E2 only:
   - spec CRs: CR-U 0.565 + CR-K S5;
   - fall-back 1: CR-U from 0.575 + CR-K with M(pipi) >= 0.44;
   - fall-back 2: CR-U from 0.575, no CR-K.

   The A1 and A2 values are in Sect. 4. The spec CRs fail, and so does fall-back 1. **Adopted: fall-back 2 (CR-K dropped).**
6. **Model choice (D11) with the template MC error included.**
   - The data-only chi2 has no template uncertainty, while the templates are 10-50 x less precise than the data in CR-U. The data-only
     dchi2 therefore rewards terms that follow MC fluctuations. Example: the omega split gives c_om = 3.4 +- 0.5 (stat) +- 2.2 (MC).
   - I therefore evaluate D11 and the chi2 criterion with a GLS chi2: data variance + MC covariance of the predicted spectrum (from the
     100 replicas, 3 iterations). This is also how the sim-study fits behind D11's "3.8 -> 1.3" were made (A2 note: "iterated GLS
     chi2 with template covariance").
   - In fall-back 2 (CR-U 0.575-0.700, 125 bins):

     | model | GLS chi2/ndf | GLS dchi2 | data-only chi2/ndf |
     |---|---|---|---|
     | M0 scale only | 138.8/119 | -- | 443.4/119 |
     | M1 common tilt (adopted) | 129.8/118 | 9.0 for the tilt | 326.4/118 |
     | M2 + quadratic | 127.9/117 | 1.9 | 252.3/117 |
     | M3 per-slice tilts | 127.5/114 | 2.4 | 271.8/114 |

   - -> **M1 adopted**, chi2/ndf = 1.10.
   - The data-only D11 would give M4 (per-slice tilts + quadratic), then the omega split: data-only chi2 139.0/112. It is kept as the
     variants D11data and D11dataOmega (an S3 family member).
   - The parameter VALUES come from the spec's estimator (Poisson ML, data variance). The GLS estimates agree within the MC error
     (stageA.txt `values_gls`).
7. **Tests A1-A5, variants, closure in CR-K and CR-L, Tables C1-C3, Figs C1-C3** (below).

---------------------------------------------------------------------------------------------------
## 2. Numbers

### 2.1 Frozen stage-A parameters

Model M1, E2 only, CR-U 0.575-0.700, S0, a_eta(0) = 0.4; 125 bins.

| parameter | value | stat (sandwich) | MC (100 replicas) | total | d/d a_eta |
|---|---|---|---|---|---|
| c_1 (S1) | 0.9878 | 0.0427 | 0.1929 | 0.1976 | -0.062 |
| c_2 (S2) | 1.0012 | 0.0440 | 0.1991 | 0.2039 | -0.097 |
| c_3 (S3) | 0.9751 | 0.0445 | 0.2066 | 0.2113 | -0.092 |
| c_4 (S4) | 1.0552 | 0.0429 | 0.1866 | 0.1915 | -0.054 |
| c_5 (S5) | 0.9091 | 0.0335 | 0.1425 | 0.1464 | -0.004 |
| b [per 0.1 GeV/c2] | -0.0683 | 0.0057 | 0.0242 | 0.0249 | +0.003 |
| kappa_B | 3.06 | 0.29 | 1.42 | 1.45 | -0.17 |
| k_E2 (constr. 1 +- 0.05) | 0.954 | 0.048 | 0.024 | 0.054 | +0.011 |

- **Goodness of fit:** chi2 (with the MC covariance) = 129.8/118 = 1.10. Data only: Pearson 326.4/118, Baker-Cousins 374.3.
- **Correlations (total):** c_j - c_k 0.97-0.98; c_j - kappa_B -0.91 to -0.93. In E2, G_C and G_D have similar shapes, so only their
  sum is constrained. The per-bin background B_A is much better known than c_j or kappa_B separately (Sect. 2.3).
- **Smaller error sources:**
  - Event groups: sandwich / model error = 1.002-1.007.
  - CB: 200 CB replicas give 0.007-0.009 on c_j. The diagonal BB-lite treatment gives half of that. Both are negligible against the
    MC error.
- **GLS estimates (cross-check):** c_j 0.887-1.018, b = -0.051 +- 0.023, kappa_B 3.19 +- 1.04.
- **Parameters at other a_eta(0)** (same model and CR; full rows in stageA.txt):

  | a_eta(0) | c_1 | c_5 | b | GLS chi2/118 |
  |---|---|---|---|---|
  | 0 | 1.007 | 0.906 | -0.0696 | -- |
  | 1 | 0.944 | 0.902 | -0.0669 | 153.4 |
  | 2 | 0.866 | 0.882 | -0.0606 | -- |

  The CR data prefer a small a_eta through G_eta (8.5 % of CR-U E2 at SMASH).
- **Tagged rerun.** The full rerun at a_eta(0) = 1 (`stageA_aeta1.txt`) reproduces `values_aeta1.0` exactly.

### 2.2 Variants

Same model and CR unless stated; stat errors; full list in stageA.txt and wp3_stageA_log.txt.

| variant | chi2/ndf (data only) | c_1 | c_5 | b | kappa_B |
|---|---|---|---|---|---|
| adopted | 326.4/118 | 0.988 | 0.909 | -0.068 | 3.06 |
| CR-U 0.565-0.700 | 349.7/133 | 0.958 | 0.885 | -0.064 | 3.24 |
| spec CRs (+ CR-K) | 403.5/148 | 0.920 | 0.855 | -0.065 | 3.56 |
| fall-back 1 (CR-K >= 0.44) | 375.1/133 | 0.930 | 0.865 | -0.071 | 3.53 |
| kappa_B = 1 +- 0.2 | 369.6/119 | 1.184 | 1.056 | -0.061 | 1.94 |
| CB per-cell geometric mean (S1) | 296.4/118 | 1.053 | 0.980 | -0.070 | 2.56 |
| ghost variant (a) | 295.1/118 | 0.887 | 0.857 | -0.084 | 2.78 |
| K0s veto (A5) | 350.8/118 | 0.971 | 0.940 | -0.073 | 3.19 |
| D11 data-only (M4: b_j 0.25-0.38, q = -0.21) | 186.0/113 | 1.023 | 0.831 | per slice | 2.29 |
| D11 data-only + omega (c_om = 3.36) | 139.0/112 | 0.831 | 0.671 | per slice | 2.68 |
| pre (parallel column) | 607.5/118 | 1.212 | 1.110 | -0.088 | 2.52 |
| step 0: E1+E2, k_E1 free (M4) | 760.1/297 | 1.243 | 0.979 | per slice | 0.95 (k_E1 13.2) |

### 2.3 Table C3: counted "after" display

- **What is shown:** data S0, E2, raw rows. B_A = stage-A background. eta(1) = all eta components at SMASH.
- **D20:** these are display counts, NOT yields; they are not to be quoted before WP4 passes.
- **Ratio** = (N - k CB)/(B_A + 0.4 eta).

| region | N | k CB | B_A +- stat +- MC | eta(1) | N - kCB - B_A | ratio |
|---|---|---|---|---|---|---|
| 0.400-0.455 | 13786 | 1379 | 12897 +- 188 +- 879 | 5237 | -490 +- 907 | 0.828 +- 0.050 |
| B 0.4550-0.5225 | 53241 | 5409 | 38937 +- 476 +- 1941 | 28290 | 8895 +- 2012 | 0.952 +- 0.038 |
| 0.5225-0.5336 | 11471 | 1281 | 10152 +- 107 +- 390 | 1943 | 39 +- 419 | 0.933 +- 0.036 |
| P 0.5336-0.5542 | 27099 | 2863 | 22173 +- 214 +- 669 | 7120 | 2064 +- 721 | 0.969 +- 0.028 |
| 0.5542-0.565 | 15519 | 1736 | 13698 +- 121 +- 300 | 1633 | 85 +- 347 | 0.960 +- 0.023 |
| CR-U 0.575-0.700 | 342098 | 31203 | 302314 +- 1713 +- 776 | 21961 | 8581 +- 1971 | 0.999 +- 0.006 |
| CR-L 0.5225-0.5300 | 7519 | 850 | 6690 +- 71 +- 264 | 1345 | -21 +- 287 | 0.923 +- 0.037 |
| CR-K S5 0.400-0.5225 | 331 | 146 | 292 +- 8 +- 52 | 7 | -107 +- 55 | 0.63 +- 0.13 |

### 2.4 Closure

Prediction from the adopted fit, E2.

| region | N | predicted | pull | S10 term max(\|pred - obs\|, sigma_pred) |
|---|---|---|---|---|
| CR-L | 7519 | 8078 +- 522 | -1.06 | 559 rows = 8.4 % of N - kCB |
| CR-K | 331 | 440 +- 32 | -3.0 | 109 rows |

- The CR-L term is the S10 entry required by fall-back 2.
- CR-K is rejected; there the S5 threshold region is mis-described (Sect. 5).

### 2.5 Row-level A2 inputs (wp3_rows.txt)

**Validation:**
- The old control slice (none, 0.42-0.46, edge 0.5150) is reproduced: 36.48 +- 19.29 (36.5 +- 19.3 quoted).
- ECAL tag, none: T = 3301.1 +- 107.8, S/B 0.667 (quoted 3301 +- 108, S/B 0.667).
- The spread-rule CB of CR-K equals WP0: 199.52.

**S0 results:**

| test | CB spread | CB own |
|---|---|---|
| edge step, 0.43-0.46 | h = 19.2 +- 1.9 (chi2/ndf 113/32) | 23.3 +- 2.0 |
| edge step, 0.44-0.46 | 10.0 +- 1.3 (67/32) | -- |
| edge step, 0.46-0.55 | -0.15 +- 0.97 (19/32) | 1.70 +- 0.98 |

| tagged sample | T |
|---|---|
| CR-K | -2.1 +- 4.0 |
| bump region (positive control) | 3252.6 +- 91.9 |
| all M4 | 3090 +- 94 (S/B 0.90) |

---------------------------------------------------------------------------------------------------
## 3. Figures and tables (draft captions)

- **Fig. C1** `docs/latex/figs_all/wp3_map.pdf` (.png).
  - Content: data N - CB in M(pi+pi-e+e-) x M(pi+pi-), S0, isBest >= 0, all M(ee).
  - Units: raw rows per 10 x 10 (MeV/c2)^2. CB = 2 sqrt(N++ N--), spread rule.
  - Boxes:
    - CR-U: adopted from 0.575, thick green; the spec's 0.565, thin;
    - CR-K: M(pipi) >= 0.43, M4 < 0.5225, magenta, "dropped";
    - CR-L, orange;
    - bump window B (blue, S1-S4) and peak window P (red).
  - Lines: the bump edge E(m, q) for q = 0 (solid) and q = 60 MeV/c2 (dashed); dotted lines = the slice boundaries.
  - No simulation.
- **Fig. C2** `wp3_ratio.pdf`.
  - Content: per slice S1-S5, M(ee) >= 10 MeV: (N - k CB)/(B_A + 0.4 eta), with B_A the stage-A background.
  - Filled points: bins of the fit (green band, CR-U 0.575-0.700). Open points: outside the fit.
  - Grey: +- stage-A error (stat (+) MC). Red: the ratio expected if the eta were at SMASH (a_eta = 1) instead of 0.4.
  - S0, isBest >= 0; simulation SMASH mode A, w x V3, L = 6079.2 nb^-1 x 63/64.
- **Fig. C3** `wp3_after_E2.pdf`.
  - Content: "after" spectrum, M(ee) >= 10 MeV, all slices. Data N - k CB (blue) in nb/(GeV/c2), x 64/63 / 6079.2 nb^-1.
  - Stack:
    - renormalised G_D and G_C (stage A);
    - G_eta, C2, C3, C4, C5 at SMASH x PDG BR (a_eta = 1);
    - dashed: total with all eta x 0.4.
  - Ratio panel: data / (B_A + eta at 0.4), black; data / (B_A + eta at SMASH), red.
  - Green bar: CR-U. In E2, G_D and G_C cannot be separated (kappa_B - c_j correlation -0.92); only their sum is constrained.
- **Smoothing (procedure figure)** `wp3_smoothing.pdf`: simulated G_D and G_C per F10 histogram (points, event-group errors) and the
  exp(pol3) smoothing (lines). Units: expected raw rows per MeV/c2. S0, mode A, w x V3.
- **A2 figures:**
  - `wp3_crk_edge.pdf`: N (black), CB (red), N - CB (blue) in the CR-K slices 0.43-0.46 and 0.46-0.55; S0, raw rows per 2.5 MeV/c2.
    - The pol2 + step fit (green) has the step at E(0.4116, 0) = 0.4739 (dotted). Left: spread-rule CB; right: CB of the slice itself.
    - The data are empty below the phase-space threshold of the slice (about 0.47 and 0.50). The step absorbs the threshold.
  - `wp3_crk_tag.pdf`: ECAL-tagged eta, N - CB in m5c. Left: CR-K, T = -2.1 +- 4.0. Right: M(pipi) < 0.43, M4 < 0.5225, T =
    3253 +- 92 (positive control).
- **Tables:**

  | table | file | content |
  |---|---|---|
  | C1 | `docs/latex/tables_all/wp3_tableC1.{txt,tex}` | candidate CRs (N, CB/N, C2-C5 and G_eta share, data/sim, verdict) |
  | C2 | `wp3_tableC2.{txt,tex}` | nested models in 4 configurations (data-only and GLS chi2), adopted parameters, MC and total errors, correlation, d/d a_eta |
  | C3 | `wp3_tableC3.{txt,tex}` | counted "after" display |
  | A1 detail | `wp3_A1_bins.txt` | every CR bin of the 3 configurations |
  | acceptance | `wp3_acceptance.txt` | the tests of Sect. 4 |

---------------------------------------------------------------------------------------------------
## 4. Acceptance criteria (spec WP3) and fall-backs

| test | value | verdict |
|---|---|---|
| chi2/ndf <= 1.5 | GLS 129.8/118 = 1.10 (data only 2.77) | **PASS** (data-only: FAIL, explained by the template MC error) |
| A1 sim C2-C5 < 1 % in every CR bin | spec CRs: max 7.7 % (8/151 bins); fall-back 1: 6.3 % (5/136); adopted: 2.77 % in 1 bin (S2E2 0.605-0.610 = ONE simulated C2 row, N_eff 1.0); relative to the fitted model (a_eta 0.4): max 0.88 % | literal **FAIL** (single MC row); PASS against the fitted model. Integrated: C2-C5 0.022 %, G_eta 8.5 % (max 26 % in a bin) |
| A2 CR-K edge step \|h\| < 2 sigma | 0.43-0.46: 10 sigma, chi2/ndf 3.5 (pol2 cannot describe the S5 threshold; only 15 rows lie below the edge + 2 sigma); 0.46-0.55: -0.15 +- 0.97 | **FAIL** (test ill-posed) / PASS |
| A2 tagged m5c in CR-K = 0 | T = -2.1 +- 4.0 (in CR-K m5c >= M(pipi) + m_pi0 = 0.565, so it is zero by kinematics) | **PASS** (trivially) |
| A3 eta injected x {0, 0.4, 1, 2}: c_j, b move < 0.3 sigma | C2-C5 injected: 0.045 sigma | **PASS** |
| A3 with G_eta also injected | 0.68 sigma (by design, D10) | information |
| A4 per-slice tilts consistent with common b | GLS dchi2 2.38/4, p = 0.67 (data only p = 4e-11) | **PASS** |
| A4 CR-U split at M(pipi) 0.42 | data/model 1.0010 +- 0.0138 vs 1.0013 +- 0.0152, -0.02 sigma | **PASS** |
| A5 K0s veto: shift < 1 sigma | max 0.68 sigma (k_E2) | **PASS** |

**Fall-backs that fired (Sect. 6):**
1. k_E1 = 13.2 is outside 0.8-3 -> **drop E1**. The "kappa = 1 +- 0.2" part applies to the stage-B signal factor; kappa_B stays free
   (Sect. 5b).
2. A1/A2 fail -> **CR-K at 0.44, CR-U from 0.575** -> still fail.
3. -> **Drop CR-K**:
   - b comes from CR-U only;
   - the bump is labelled **"extrapolated"**;
   - the CR-L term (559 rows, E2) is added to S10.
4. **MC error > data statistical error -> "simulation-limited"**, request a larger std production. The MC error is 0.19 against stat
   0.04 on c_j.

**Not fired:**
- stage-A chi2/ndf > 1.5: GLS 1.10, so no omega split and no error scaling. With the literal data-only chi2 it fires, and the chain
  ends at M4 + omega split, chi2/ndf 1.24; that is the variant D11dataOmega.

---------------------------------------------------------------------------------------------------
## 5. Caveats and deviations (flagged)

a. **D11 / chi2 with the MC covariance (deviation).**
   - The literal data-only chi2 makes D11 add a quadratic term, per-slice tilts and the omega split. All three follow template MC noise
     (GLS dchi2 1.9, 2.4 and 3.5, against data-only 74, 66 and 118).
   - c_om = 3.4 has an MC error 4 x its stat error.
   - The analyst may prefer the literal choice: it is fully stored (D11data, D11dataOmega) and belongs to the S3 family anyway.
b. **E1 dropped, kappa_B free.**
   - Constraining kappa_B = 1 +- 0.2 (variant kBcons) pulls it to 1.94 and worsens the fit (data-only chi2 369.6 against 326.4).
   - kappa_B = 3.1 +- 1.4 in E2 is nearly degenerate with c_j. Do not read it as a measured conversion scale.
   - E1 data cannot be described with the E2 parameters: the predicted/observed ratio is 0.65-0.80 (log, Table C3 "all M(ee)",
     information only).
c. **The stage-A background overshoots the data below the bump** (E2):
   - in 0.400-0.455 even the background alone exceeds N - k CB: 12897 vs 12407 (ratio 0.83 with 0.4 eta);
   - CR-K is low: pull -3.0;
   - CR-L is low by 1.1 sigma.
   The extrapolation into B (40-120 MeV below CR-U) is the weak point. Stage B must keep the background free (as specified), and WP4/WP6
   should test it. The B window display (N - kCB - B_A = 8895 +- 2012 against 0.4 eta(1) = 11316) is not a yield.
d. **CB in CR-K.**
   - The spread rule puts CB = 46 % of N in S5E2 (N++ = 2), with a band-edge jump at 0.500.
   - The per-cell geometric mean gives 6 x less.
   - The CR-K normalisation changes from c' = 0.67 to 1.11 between the two (A2 template test). Another reason CR-K is unusable.
e. **A2 as specified is ill-posed.**
   - In CR-K the M4 spectrum starts at the slice's phase-space threshold, right at E(0.4116, 0) = 0.474.
   - The sound statements are:
     - only 15 rows (S0) lie below the edge + 2 sigma;
     - the ECAL tag is empty (kinematically required).
   - Template version (stage-A shape + step): h = -0.35 +- 0.14 (spread CB) and -2.4 +- 0.5 (per-cell CB), chi2/ndf 2.8-5.3. A
     NEGATIVE step is not an eta leak. It shows that the threshold is mis-modelled.
f. **A1 literal** fails only on single heavy simulated rows (1 std row = 22-38 expected data rows here). The C5 radiative tail in CR-K
   (1-3 % per bin, N_eff 9-29) is real; CR-K is dropped anyway.
g. **Implementation details.**
   - CR-K at 0.44 (fall-back 1): an extra G3-built S5 histogram with M(pipi) >= 0.44, used only for the CR-K bins (same c_5).
   - A5: K0s veto and its reference are both rebuilt from G3. The rebuild alone shifts nothing (0.000 sigma).
   - The CR-U split at 0.42 uses raw G3 templates.
h. **pre (parallel column)** is fitted with the adopted model and CR only. It has no GLS and no MC replicas; data-only chi2/ndf 5.1.
i. **Blinding.** Only F10 E2 is used (coarse, unmasked). No X2D or B_mee cell of 12-28 MeV was read. D20: Table C3 / Fig C3 are
   displays, not yields.
j. **A one-sided variant.** The a_eta(0) = 1 GLS chi2 (153.4) is worse than at 0.4 (129.8). Rerun at the WP5 value with the tagged
   call below; stageA.txt also gives dtheta/da_eta.

---------------------------------------------------------------------------------------------------
## 6. Exact commands (from /home/przygoda/HADES/FAT/FAT)

Row-level part (15 s, 2.6 GB). It writes wp3_rows.{root,txt}, wp3_rows_summary.txt and figs wp3_crk_edge, wp3_crk_tag:

    root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_wp3", kTRUE)' 'plots/all_report/wp3_rows.C+O'

Stage A (3 s, 0.6 GB). It writes stageA.{txt,root}, wp3_stageA_log.txt, the tables and figs; the output is deterministic,
stageA.txt md5 cfc30e3a...:

    root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_wp3", kTRUE)' 'plots/all_report/wp3_stageA.C+O'

Rerun at another a_eta(0). It writes only stageA_<tag>.txt/.root and its log, for example:

    root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_wp3", kTRUE)' 'plots/all_report/wp3_stageA.C+O(1.0,"aeta1")'

**stageA.txt keys:**
- model and CR: `model`, `control_region`, `cr_u`, `cr_k`, `meec_classes_used`, `bump_background_label`;
- parameters: `names`, `fixed`, `constraints`, `values`, `values_gls`, `err_*`;
- covariances: `cov_stat`, `cov_mc`, `cov_total`;
- a_eta dependence: `dtheta_daeta`, `values_aeta{0,1,2}.0`;
- templates and background: `smooth <hist>` (exp(pol3) parameters of G_D and G_C), 700 `bin` lines (B_A, stat, MC, smoothed G_D,
  G_C, in-CR flag);
- `variant ...`;
- `closure CR-K/CR-L`.
