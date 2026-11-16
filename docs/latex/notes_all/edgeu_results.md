# edgeu: the eta bump in the distance to its kinematic edge, u = M4 - E(M(pi+pi-), M(e+e-)) -- results note

Date 2026-11-15. Branch pp45_pippimepem (worktree /home/przygoda/HADES/FAT/FAT). No git command was run. No existing file was
edited, moved or removed: report_all.tex, defs.h, SimEngine.h, wp4_lib.h, wp5_lib.h, wp6_ana.C, every bedge/vbedge file and
the data/simulation compact files were only read. All outputs are new files with the prefix edgeu. At most 2 processes ran at
a time (one ROOT dump of about 10 GB, otherwise numpy); CPU package 38-41 C.

**Machine-readable outputs** (`/hdd3/przygoda/HADES/FAT_out/all_report/edgeu/`): `edgeu_ana.json` (md5 81831c42...),
`edgeu_toys.json`, `edgeu_bgsmooth.json`, `edgeu_numbers.json`; readable logs `edgeu_ana.log`, `edgeu_toys.log`,
`edgeu_bgsmooth.txt`, `edgeu_combine_stdout.txt`, `edgeu_dump.log`.

---------------------------------------------------------------------------------------------------
## 0. In plain words

**What the variable u is, and why it helps.**
- In the bump, an eta decays into pi+ pi- pi0, the pi0 gives an e+e- pair, and one photon is not detected. The four detected
  particles then always have a mass M4 below a limit E. The limit is reached when the lost photon is as soft as possible.
- That limit depends on two things that are measured in every event:
  - M(pi+pi-): it fixes how fast the pi0 moves in the eta frame;
  - M(e+e-): it fixes how much energy the photon gets in the pi0 frame (a heavier e+e- pair leaves less for the photon).
- So every bump event carries its own limit. The new variable is the distance to that limit, u = M4 - E. For every bump event
  u <= 0, whatever its M(pi+pi-) and M(e+e-). In u all the edges sit on top of each other at u = 0.
- The per-slice method (bedge) had to cut M(pi+pi-) into four slices and fit each slice on its own, with an edge that moves
  inside the slice by 5-7 MeV/c2 and ignores M(e+e-). In u there is one sharp step at u = 0 and one fit.

**What was done, step by step.**
1. The frozen edge function of defs.h already contains the M(e+e-) term. It was checked against the report values and against
   the simulated truth: not one simulated bump row (C2 Dalitz, C3 conversion) lies above its true limit.
2. Data and simulated rows were put into u, with the data mass scale for the data and the transferred scale for SMASH.
3. The u spectra were drawn for the data (N - CB, S0 and pre), for the like-sign rows and for every simulated class.
4. One fit was made in u: a smooth background (pol2) plus a step of height h at u = 0, smeared with the simulated resolution.
   The eta peak, which sits at u = +20 to +55 MeV/c2, was taken from SMASH and scaled by the measured peak scale a_P.
5. The same fit, applied to the simulated bump alone, gave h_sim. The eta scale of the bump is f = h_data / h_sim.
6. Errors: 200 bootstrap copies of the data (stat), 100 simulation replicas (MC), 70 variants (systematics, Part III rules).
7. Checks: each M(pi+pi-) slice, each M(e+e-) class, periods, folds, isBest == 1, and null tests (like-sign rows, simulated
   backgrounds, a control region without bump).
8. 21 pseudo-data truths x 1000 toys, including the quadratic background term that destroyed B4.
9. Conversion to yields and sigma(pp -> eta X), and comparison with all other routes.

**What came out.**
- **The step is very clear.** h_data = 1305 +- 62 rows per 2.5 MeV/c2; without the step the fit is worse by delta chi2 = 439
  (21 sigma). The fit describes the data well (chi2 27.3/32).
- **f_eta(u-edge) = 0.420 +- 0.020 (stat) +- 0.031 (syst) +- 0.029 (MC)** (SMASH with the C2 gamma* correction, D35);
  total 0.047 = 11 %.
  - Without the C2 correction: 0.390 +- 0.018 +- 0.026 +- 0.027.
  - The per-slice edge (bedge) gave 0.418 +- 0.022 +- 0.035 +- 0.047 (total 15 %) and 0.392 without the correction.
    **The two methods agree to 0.5 %.**
- **sigma(pp -> eta X) from the bump = 0.438 +- 0.020 (stat) +- 0.032 (syst) +- 0.030 (MC) +0.042/-0.044 (norm.) mb**
  (bedge: 0.437 mb).
- **Visible bump yield** in window B, all M(ee): 21 443 +- 1 002 +- 1 559 +- 1 329 rows (3.53 nb); all M4: 26 955 rows.
- **The checks pass, except two marginal ones that were also marginal in bedge.**
  - The four slices agree (chi2 1.1/3) and so do the M(ee) classes (E1 0.409, E2 0.433).
  - Periods: p = 0.025 (above the WP7 limit 0.01, no term added). Folds: z = -2.1 (p = 0.03).
- **The null tests:**
  - **No step** in the simulated non-eta backgrounds, nor in a data control region without bump (3 +- 17 rows).
  - **A new finding: the e-e- like-sign rows do have a step at u = 0** (74 +- 16 rows, 4.7 sigma), only for M(ee) >= 10 MeV.
    The e+e+ rows have none (2 +- 7). Because of this, the CB carries 39 rows of step (3 % of h_data). It is covered by the CB
    systematic, but its origin is not known. The per-slice method has the same feature.
- **Pseudo-data:** the estimator closes within +-2 % for smooth backgrounds; pull widths are 0.95-1.07.
  - **The quadratic background term that moves B4 by +150 % moves the u edge by -2.2 %** (per-slice edge: -2.1 %).
  - Two things do bias it, as they bias the per-slice edge:
    - a mismatch of the mass resolution between data and SMASH (+6 % for r = 0.89);
    - a background structure as narrow as the edge itself (+8 % for 5 % at 8 MeV width, +17 % for 10 % at 15 MeV).

**What it means.**
- The bump measured in u confirms the per-slice edge: the same eta scale (0.420 against 0.418), the same sigma (0.438 against
  0.437 mb). It also agrees with the peak side bands (0.429 mb; the double ratio peak/bump is 0.98 +- 0.22).
- **Is u better? Yes, but by less than one might expect: the total error drops from 15 % to 11 %.**
  - stat: 4.7 % instead of 5.4 % (uncorrected readings), because every row is aligned on its own edge.
  - **MC: 7.0 % instead of 11 % (the largest gain).** It comes from taking h_sim from the bump classes alone. The raw simulated
    background made 57 % of the bedge MC variance. This choice could also be made in the per-slice method.
  - Fitting the four slices in u separately gives the same result and errors as the single global fit (0.424 +- 0.020 +-
    0.031). **The single fit in u by itself adds nothing; the gain comes from the alignment and from the bump-only templates.**
  - syst: 7.3 % against 8.3 %. The leading items are the same in both methods:
    - the data/SMASH resolution ratio (4.5 %);
    - the step that a smooth background projects onto the edge (4.0 %).
    The slice choices (bedge E2 family, 2.9 %) disappear, and the M(ee) dependence of the edge is treated exactly.
- **Model dependence is not smaller.** Both edge methods rest on the same two assumptions:
  - SMASH has the right resolution at the edge (the data step width, 7.0 +- 0.8 MeV/c2, agrees with SMASH, 6.4 +- 0.9);
  - the background has no structure on the 10 MeV scale near the edge.
  Both are tested, not proven.
- **Recommendation.** Quote the u edge as the bump measurement of Part III: 0.420 +- 0.047; sigma 0.438 +- 0.049 (stat (+)
  syst (+) MC) +0.042/-0.044 (norm.) mb. Keep the per-slice edge as the cross-check that agrees.

---------------------------------------------------------------------------------------------------
## Headline table (S0, isBest >= 0, L = 6079.2 nb^-1; yields = raw rows x 64/63; sigma with PDG 2026 BRs, C_sec (1 - f_eta) = 1.0104 x 0.97)

| quantity | uncorrected | **C2 corrected (D35, recommended)** | C2 corrected, conversions x 0.48 |
|---|---|---|---|
| f_eta(u-edge) | 0.390 +- 0.018 +- 0.026 +- 0.027 | **0.420 +- 0.020 +- 0.031 +- 0.029** | 0.540 +- 0.025 +- 0.051 +- 0.036 |
| Y_bump, window B, all M(ee) [rows] | -- | 21 443 +- 1 002 +- 1 559 +- 1 329 | -- |
| Y_bump / L, window B [nb] | -- | 3.53 +- 0.16 +- 0.26 +- 0.22 (+ L, period 8.5 %) | -- |
| Y_bump, all M4 [rows] | -- | 26 955 +- 1 260 +- 1 960 +- 1 690 (4.43 nb) | -- |
| sigma(pp -> eta X) [mb] | 0.407 +- 0.019 +- 0.028 +- 0.028 +0.039/-0.041 | **0.438 +- 0.020 +- 0.032 +- 0.030 +0.042/-0.044** | 0.563 +- 0.026 +- 0.053 +- 0.038 +0.054/-0.056 |
| total without / with norm. | 0.043 / 0.058 mb (11 % / 14 %) | 0.049 / 0.064 mb (11 % / 15 %) | 0.071 / 0.089 mb |

Errors: stat (bootstrap, 200 copies), syst (families, Part III rules), MC (rms of 100 replicas), normalisation (L 6.9 %,
period 5 %, BR 3.17 %, C_sec (1 - f_eta) +3.1/-4.1 %).

---------------------------------------------------------------------------------------------------
## 1. Kinematics

### 1.1 The formula, in words and in symbols
- The eta decays at rest (eta frame). The pi0 has energy E_pi0 = (m_eta^2 + m_pi0^2 - m^2)/(2 m_eta) and momentum p_pi0, with
  m = M(pi+pi-). A heavier pion pair leaves the pi0 slower.
- The pi0 decays into a photon and the e+e- pair (mass q = M(e+e-); conversion pairs have q close to 0). In the pi0 frame the
  photon has the energy (m_pi0^2 - q^2)/(2 m_pi0).
- The four detected particles carry everything except the photon: M4^2 = m_eta^2 - 2 m_eta E*_gamma (E*_gamma = photon energy
  in the eta frame).
- M4 is largest when E*_gamma is smallest, i.e. when the photon flies backwards against the pi0 motion:
  E*_gamma,min = (E_pi0 - p_pi0)(m_pi0^2 - q^2)/(2 m_pi0^2). Hence

      M4max(m, q)^2 = m_eta^2 - m_eta (E_pi0 - p_pi0)(1 - q^2/m_pi0^2).

- This is exactly the frozen edge of defs.h (alldefs::Edge, report eq. (edge)), on the data scale:
  E(m, q) = 0.99274 M4max(m/0.9970, q).
  - Checks: E(0.30, 0) = 0.51859, E(0.39, 0) = 0.50077, E(0.39, 60 MeV) = 0.50958, E(0.4116, 0) = 0.47390, E(0.30, 100 MeV) =
    0.53262. These equal the report values.
  - At q = 0 it equals scan::M4max (the bedge/WP6 edge). bedge used the scale 0.99277 = 0.5439/m_eta (scan::kScale) instead of
    0.99274; the difference is 0.015 MeV.
- The M(ee) term matters: at m = 0.30 the edge moves from 0.5186 (q = 0) to 0.5221 (q = 50 MeV) and 0.5326 (q = 100 MeV).
- There is no bump above m = 0.4116 (pi0 at rest). Near that limit the bump becomes a narrow box: at m = 0.40 the box spans
  57 MeV/c2 below the edge, at m = 0.41 only 21 MeV/c2. This is why the nominal M(pi+pi-) range ends at 0.40.

### 1.2 How the mass scales enter
- **Data:** u = M4 - E(M(pi+pi-), M(ee)). M(pi+pi-) is brought to the true scale with 0.9970 (K0S in the data), the limit back to
  the measured four-body scale with 0.99274 (eta peak 0.5439/m_eta). M(ee) is not scaled.
- **Simulation:** the rows are first transferred to the data scale (D8: M4' = 0.99615 [m4_t + r (m4 - m4_t)], r = 1, M(pipi)' =
  0.9952 M(pipi)), then the same E is applied.
  - This is identical to applying the simulation-scale edge (0.99658 M4max(M(pipi)/1.00181, q), i.e. the simulated eta peak
    and K0S scales) to the raw simulated masses. The only difference is a factor 0.99615 on u: max |difference| 1.7e-16 GeV.
- **Measured in the simulation:**
  - the simulated bump pions have M(pipi) reco/true = 0.99976, not the K0S value 1.0018;
  - the simulated bump M4 reco/true is 0.99731, against 0.99658 for the peak.
  If the data pions behave the same way, the data M(pipi) scale for the edge would be 0.9952 instead of 0.9970. This is the
  variant dpp52 (-1.0 %).

### 1.3 The simulated truth (Fig. edgeu_kinematics)
- C2 (1 649 rows) and C3 (1 187 rows): **max(m4_true - M4max(M(pipi)_true, M(ee)_true)) = -0.48 MeV and -0.51 MeV.** Not one
  row lies above its true edge.
- Without the M(ee) term (E(m, 0)), 3.8 % of C2 (weighted) lie above the edge by more than 1 MeV and 1.8 % by more than 5 MeV.
  C3 (median true M(ee) 2 MeV) is not affected.
- **Reconstructed u:** above +10 MeV/c2 lie 0.85 % of C2 and 1.04 % of C3; above +20 MeV/c2, 0.08 % and 0.16 %. This is the
  resolution tail, nothing more.
  - For C2 with M(ee) >= 20 MeV, above +20 MeV/c2: 0.21 % with E(m, q), but 1.05 % with E(m, 0).
- **Resolution of u at the edge** (C2 + C3, u_true > -15 MeV/c2, 645 rows):
  - Gaussian core sigma = 6.70 +- 0.21 MeV/c2; half of the central 68 %: 7.78; rms 10.2 (tails).
  - Almost all of it is the M4 error (rms 10.3); the edge-function error from M(pipi) and M(ee) is only 1.2 MeV/c2
    (correlation 0.17).
  - The M4 resolution grows towards the edge (half of the central 68 %: 4.9 MeV/c2 more than 100 MeV/c2 below it, 7.3 MeV/c2
    within 15 MeV/c2 of it). Near the edge the e+e- pair carries almost the whole pi0 momentum. There it is similar to the
    resolution at the eta peak (simulated peak: half 68 % 7.9 MeV/c2, Gaussian core of the peak fit 6.1 MeV/c2).
- **Why one sharp step.** In the per-slice method the edge moves inside a slice by 4.7 / 6.0 / 7.1 / 6.0 MeV/c2 (slices
  0.30-0.34 / 0.34-0.37 / 0.37-0.39 / 0.39-0.40). In u this spread is gone, and so is the upward smear of the Dalitz pairs with
  large M(ee). What is left is the resolution (6.7 MeV/c2).

---------------------------------------------------------------------------------------------------
## 2. Distributions and null tests (Fig. edgeu_spectra, edgeu_checks (d))

- **Data (0.30 <= M(pipi) < 0.40, S0).** N - CB rises up to u = -5 MeV/c2 (3 870 rows per 2.5 MeV/c2), drops to 3 420 at
  +6 MeV/c2 and rises again into the peak (4 200 at +30 to +40).
  - The CB is about 200-215 rows per bin in the fit range, i.e. 6 % of N.
  - pre has 38 % more rows in the fit range and the same shape.
- **SMASH:**
  - the bump classes end at u = 0 (C2 orange, C3 magenta);
  - the peak (C4 green, C5 cyan) lies at u = +15 to +60 MeV/c2. For each M(pipi) it sits at 0.5439 - E(m, q), which spans
    +24 MeV/c2 (m = 0.28) to +49 MeV/c2 (m = 0.40);
  - G_D and G_C (greys) and G_eta (violet) are smooth within their template noise.
  - N_eff per 5 MeV/c2 bin: 64-330 for G_D, 18-105 for G_C and 14-52 for G_eta. **The simulated backgrounds are very noisy.**
- **Peak treatment.** In the fit range the peak is a fixed term a_P x K(u) (K = simulated C4 + C5 in u, a_P = 0.418249 = the
  measured peak scale of WP6/WP7). In the step region it is small: 17-48 rows per bin at 0-+10 MeV/c2 against 3 500 rows.
  - With the amplitude free the fit gives a_P = 0.33 and the same f (0.4200).
  - a_P +- 0.080 (stat (+) syst (+) MC) moves f by 0.0001.
  - So the peak treatment does not matter for the step.
- **Null tests** (step height h of the nominal estimator; h_data = 1305):

  | sample | h [rows / 2.5 MeV/c2] | verdict |
  |---|---|---|
  | data like-sign e+e+ (N++) | 2.1 +- 6.5 | no step |
  | **data like-sign e-e- (N--)** | **73.5 +- 15.7** (E1: -5.4 +- 4.6; E2: 78.7 +- 15.0) | **step, 4.7 sigma** |
  | N++ + N-- | 75.7 +- 17.0 | step |
  | the CB as subtracted (2 sqrt(N++ N--), cells) | 39.2 (3.0 % of h_data); cells 2.5/5/20 MeV: 37/33/34 | |
  | SMASH G_D (non-eta, Dalitz/vector) | 12 +- 140 (replica rms) | no step |
  | SMASH G_C (non-eta, conversion) | 59 +- 148 | no step |
  | SMASH G_eta (other eta rows) | 123 +- 57 (smoothed shapes: -5 to -31) | 2.2 sigma, smoothed: none |
  | data control M(pipi) 0.42-0.46, pseudo edge E(M(pipi) - 0.10, M(ee)) | 3.1 +- 16.7 | no step |

- **The e-e- step (new).**
  - It sits only in M(ee) >= 10 MeV and mostly in S0 (pre-not-S0: 23 +- 18).
  - Per period: I 50.8 +- 10.8, IIA 1.0 +- 3.9, IIB 7.3 +- 6.9, III 14.7 +- 8.2.
  - The simulated like-sign rows show no such step (e-e-: -5 +- 14, Poisson-like error; edgeu_dumpls.C, scratch analysis).
  - Its origin is not identified. Two readings are possible:
    - genuine combinatorial pairs whose unlike-sign partners the CB removes correctly;
    - a like-sign artefact of bump events (for example a wrong charge sign of the Dalitz e+). Then the CB over-subtracts.
  - The CB variants bracket it:
    - CB shape from e+e+ only (no step): +1.1 %;
    - CB shape from e-e- only: -1.3 %;
    - CB = N++ + N-- (the full e-e- step): -2.6 %.
    It is counted in U5 (2.6 %).
  - It affects every CB-subtracted bump measurement in the same way, the per-slice edge included.

---------------------------------------------------------------------------------------------------
## 3. The fit and its errors (Fig. edgeu_fit)

**Model.** N - CB - a_P K(u) = pol2(u) + h T(u). The pieces:
- range: -45 <= u < +45 MeV/c2 (36 bins of 2.5 MeV/c2, with a bin edge at u = 0);
- T: the bin average of 0.5 erfc(u/(sqrt2 sigma)), sigma = 6.70 MeV/c2 (simulated core);
- V: cov(CB) + diag(max(model of N, 1)), iterated as scan::FitEdgeStep;
- CB: 2 sqrt(N++ N--) in 10 MeV/c2 cells of u, interpolated, per M(pipi) slice 0.30-0.34-0.37-0.39-0.40 (the ScanEngine rule
  applied to u, as scan_bump.C did for its aligned spectra);
- h_sim: the same linear estimator (final data covariance) applied to the simulated bump C2 + C3 (C2 x 0.879 / 0.8825 for
  E1 / E2).

**Choice of the range** (scan, Sect. 4 of edgeu_ana.log):
- pol1 is rejected at every range >= +-35 MeV (chi2 115.6/33 at +-45).
- pol2 describes +-35 to +-55 MeV with p = 0.47-0.70. +-60 gives p = 0.15, +-70 p = 0.03.
- pol3 gives the same f (0.418) with an error 1.8 times larger: the cubic term and the step are both odd around u = 0.
- **Nominal +-45 MeV/c2, pol2.** Reasons:
  - it is the half range of the per-slice method;
  - the lower end stays above -57 MeV/c2, where the first bump boxes (m = 0.40) end;
  - the upper end stays below the region where pol2 starts to fail.
- f is stable over +-40 to +-55 MeV (0.420-0.429).

**Nominal result.**
- h_data = 1305.2 +- 62.3 rows per 2.5 MeV/c2; chi2 27.3/32 (without step 466.6).
- h_sim = 3107.3 (C2 + C3, C2 corrected). Its linear decomposition: C2 1780.3 (raw 2020.1), C3 1327.0; the C2 share is 57 %.
- **f_eta(u-edge) = 0.4200.**
  - stat (fit) 0.0200; stat (bootstrap) 0.0196, bootstrap/fit 0.98;
  - MC (replica rms) 0.0292 (robust 0.0253); the replica mean is 0.4258;
  - the h_sim replica rms is 201 rows (6.5 %; bedge E_sim 288 rows, 9.9 %).
- Uncorrected: 0.3900 (stat 0.0182, MC 0.0267). kappa = 0.48 reading: 0.5400 (stat 0.0252, MC 0.0364).
- **Step width in the data:** free sigma = 7.02 +- 0.79 MeV/c2. The same fit on the Asimov spectrum of SMASH gives
  6.42 +- 0.86 MeV/c2.
  - The data/SMASH ratio from the step is 1.09 +- 0.12; from the eta peak width it is 0.893 +- 0.099. Combined: 0.97 +- 0.08.
  - So the resolution check on the data passes. The nominal (equal smearing) is kept; r = 0.89/1.10 and the own-smearing variant
    give the systematic.
- **Four-slice estimator in u** (one pol2 per slice; f = sum h_data / sum h_sim): 0.4235 +- 0.0199 (boot) +- 0.0307 (MC).
  This is the same as the global fit.

**Per slice and per class** (Fig. edgeu_checks (a), (b)):

| slice | h_data | h_sim | f (C2 corr.) | stat | MC | uncorr. | bedge uncorr. (data err.) |
|---|---|---|---|---|---|---|---|
| 0.30-0.34 | 325.5 +- 38.9 | 701.8 | 0.464 | 0.053 | 0.054 | 0.430 | 0.342 (0.044) |
| 0.34-0.37 | 436.7 +- 36.2 | 1077.5 | 0.405 | 0.038 | 0.050 | 0.375 | 0.459 (0.039) |
| 0.37-0.39 | 387.0 +- 27.6 | 1032.1 | 0.375 | 0.026 | 0.042 | 0.353 | 0.376 (0.027) |
| 0.39-0.40 | 145.4 +- 16.5 | 245.7 | 0.592 | 0.073 | 0.56 | 0.525 | 0.389 (0.062) |

- chi2 with stat (+) MC: 1.1/3; stat only: 9.1/3.
- The u step per slice is about 13 % higher than the M4 step of bedge (sum 1295 against 1142). The step in u is sharper, and the
  template width differs (6.7 against 5.4 MeV/c2). The ratio to SMASH is what counts.

| M(ee) class | f (C2 corr.) | stat | MC | uncorr. | own-class CB | bedge C2 corr. |
|---|---|---|---|---|---|---|
| E1 (< 10 MeV) | 0.409 | 0.026 | 0.047 | 0.389 | 0.409 | 0.384 |
| E2 (>= 10 MeV) | 0.433 | 0.031 | 0.041 | 0.391 | 0.432 | 0.475 |
| >= 20 MeV | 0.355 | 0.045 | 0.041 | 0.313 | 0.344 | -- |

- **E1/E2 split** (h_data,i = a (C2c_i + kappa C3_i)): a_D = 0.444 +- 0.046 (stat) +- 0.063 (MC), **kappa = 0.88 +- 0.21 (stat)
  +- 0.35 (MC)**. bedge had 0.51 +- 0.25 +- 0.48.
  - In u the conversion scale prefers about 0.9 but cannot exclude the peak value 0.48 (1.1 sigma).

**Periods, folds, selections** (Fig. edgeu_checks (c)):
- Periods (x C_sec(p) = 1.0241, 1, 1, 1; simulation x L_p/L): I 0.433 +- 0.029, IIA 0.391 +- 0.083, IIB 0.531 +- 0.045,
  III 0.356 +- 0.036.
  - chi2 9.4/3, p = 0.025 (raw p = 0.026). Above the WP7 limit p < 0.01, so no term is added.
  - It is the same pattern as bedge (IIB high, III low; p 0.07).
- Folds: 0.376 +- 0.028 / 0.463 +- 0.029, z = -2.14 (p = 0.032). **Marginal**, as in bedge (z 1.9-2.0, p 0.05-0.06). Folds
  are a random split, so a 2.1 sigma difference happens in 3 % of cases.
- pre (data and simulation): 0.408 (-2.8 %, noise 3 %). isBest == 1 in data and simulation: 0.449 (ratio 1.068; bedge 1.085).
  isBest == 1 in the data only (WP7 convention): 0.368 (ratio 0.877; WP7 0.885).

---------------------------------------------------------------------------------------------------
## 4. Systematics (Table tables_all/edgeu_syst.{txt,tex}; Fig. edgeu_syst)

**Part III rules (WP7).**
- delta^2 = max(0, Delta^2 - |stat_v^2 - stat_0^2| - MC(Delta)^2); MC(Delta) = the robust replica spread of the
  replica-by-replica difference.
- Next to it the bootstrap version: the data noise of Delta from the same 200 bootstrap copies.
- A +-1 sigma pair counts (|D+| + |D-|)/2. A family takes the maximum; families add in quadrature.
- R1: a variant the data reject is listed, not counted.

| family | content (largest item) | literal [%] | bootstrap [%] |
|---|---|---|---|
| U1 fit model | upper end +60/+30 MeV pair (hiP30 +7.0 %, noise 6.0 %, counted 3.5 %; pair 1.8 %) | 1.8 | 1.7 |
| **U2 resolution** | own smearing per sample (data 5.99, sim 6.70 MeV: -4.7 %); r = 0.89/1.10 pair 4.4 % | **4.5** | **4.7** |
| U3 M(pipi) range | 0.28 / 0.29 / 0.31 start, 0.39 / 0.41 end: all within noise | 0 | 0 |
| U4 scales, edge function | data M(pipi) scale 0.9952 (-1.0 %); data edge scale x 1.0009 (+0.8 %) | 1.0 | 0.8 |
| U5 CB | CB = N++ + N-- (-2.6 %, the e-e- step); M4-rule CB mapped to u (-1.4 %) | 2.6 | 2.4 |
| U6 binning | 5 MeV bins (-1.2 %) | 0.9 | 0.6 |
| U7 physics weights | conversion F_V off (peak template, +1.6 %); Dalitz plot 1 + aY (-1.4 %) | 1.6 | 1.6 |
| U8 peak term | peak template r = 0.89 (-0.8 %); a_P +-0.080: 0.02 % | 0.8 | 0.8 |
| U9 selection | pre: -2.8 % with noise 3.0 % | 0 | 0 |
| U10 C2 factor +-0.024 | +-1.6 % | 1.3 | 1.6 |
| **U11 background step** | stage-A-like smooth background in the pseudo-data: +4.0 % (envelope of 8 smooth truths) | **4.0** | **4.0** |
| **total** | | **7.3** | **7.4** |

**Listed, not counted:**
- pol1 (R1: chi2 rise 102 > 25).
- Weights off / phase-space Dalitz plot (+6.4 % / +7.9 %). They are excluded by the KLOE-2 Dalitz plot and by WP7 R1 on F-main
  (delta chi2 38 > 25), not by the edge itself. If counted, U7 would become 7.2 % and the total 10.1 %.
- q = 0 (E(m, 0) for data and simulation): +5.8 % with noise 3.9 %. The M(ee) term is exact kinematics, not an uncertainty; the
  row shows what the term does. With the bedge row set (all M(ee), E(m, 0)) the result is the same, +5.8 %.
- T12 (h_sim with the raw simulated background, the bedge convention): -5.9 %. This is the template noise of the raw
  background (G_eta 123 +- 57 and non-eta 71 +- 202 rows), not a systematic.
- isBest == 1, no C2 factor, no peak term: checks.

**Why U2 and U11 dominate.**
- U2: the step height depends on the assumed width. If the data and SMASH have different resolutions, the same template
  measures their steps differently. A change of the data/SMASH resolution ratio by 0.1 moves f by about 5 %. bedge had the same:
  E3, 4.3 %.
- U11: a smooth background has an odd (cubic) part across the edge that pol2 cannot absorb. It projects onto the step.
  - Its size depends on the background shape, and the simulated backgrounds are too noisy to fix it.
  - With smoothed SMASH shapes it is -0.5 % to +2.0 % (exp(pol3-6), c = 1 or 1.35); with the stage-A-like data background it
    is +4.0 %.
  - Linear projections (edgeu_bgsmooth.txt): non-eta exp(pol2-6) on +-100 MeV gives K = -35.6, -7.8, +18.6, +9.3, +17.3 rows
    (h_sim 3107); the smoother fits have p = 0.01-0.06.
  - bedge had E7 = 6.0 % for the same effect.

---------------------------------------------------------------------------------------------------
## 5. Pseudo-data (Fig. edgeu_toys; Table tables_all/edgeu_toys.{txt,tex}; edgeu_toys.py)

**Truth** (the WP4 recipe, in u):
- the data CB per slice;
- c x the non-eta background of SMASH (summed over the slices, smoothed by exp(pol4) on +-100 MeV), x (1 + b x + q x^2),
  x = (M4' - 0.55)/0.1, with optional local Gaussian structures in M4 and an M(pipi) tilt;
- a x (C2c + kappa C3) + a x G_eta (smoothed);
- a_Pt x (kappa C4 + C5).

**Toys:** Poisson N; Poisson like-sign; CB rebuilt per slice; the nominal estimator. Two template choices:
- AA: templates = the truth sample;
- RF: bump and peak templates from WP0 replica toy mod 100.

| truth | bias (Asimov) | shift vs (a) | AA pull SD | RF bias | RF pull SD |
|---|---|---|---|---|---|
| (a) c 1.35, b -0.15, a 0.42 | +2.0 % | -- | 0.99 | +2.9 % | 0.95 |
| c = 1, no tilt | +1.3 % | | 1.03 | +2.6 % | 1.01 |
| a = 1 | +0.9 % | | 1.02 | +1.9 % | 0.90 |
| a = 0 (no eta) | f = 0.009 | | 1.00 | | |
| **(a) + q = -0.173 (B4 +150 %)** | -0.2 % | **-2.2 %** | 1.03 | +1.0 % | 1.00 |
| (a) + q = -0.346 / +0.173 | -2.5 / +4.3 % | -4.5 / +2.2 % | 0.99 / 1.01 | | |
| WP4 (e) omega/rho0 x 0.3, rest x 1.6 | +2.5 % | +0.5 % | 1.02 | | |
| M(pipi) tilt 1 + 0.3 y | -1.3 % | -3.4 % | 1.05 | | |
| stage-A background (G_D + 2.93 G_C) | +4.0 % | +2.0 % | 1.04 | +5.5 % | 0.97 |
| smoother exp(pol3/5/6) | -0.5 / +0.8 / +1.0 % | | 1.02-1.03 | | |
| data-like (kappa 0.48, c 1.16) | +1.1 % | | 1.00 | +2.6 % | 0.98 |
| kappa 0.48 in the truth | +1.7 % (of what f measures) | | 0.98 | | |
| true peak scale 0.30 / 0.55 | +2.1 / +2.0 % | +0.1 / -0.1 % | | | |
| eta truth at s + 0.0009, r = 0.89 | +8.0 % | +6.0 % | 1.07 | | |
| raw SMASH background (template noise) | +12.1 % | +10.0 % | 1.02 | | |
| 5 % Gauss(0.505, 8 MeV) in M4 | +9.6 % | +7.6 % | 1.04 | | |
| 10 % Gauss(0.490, 15 MeV) in M4 | +19.3 % | +17.3 % | 1.05 | | |

- **Coverage of +-1 sigma:** 0.63-0.68 (AA) and 0.66-0.73 (RF) for the truths with |bias| < 2.5 %. **Pull widths:** AA
  0.98-1.07; RF 0.90-1.01 (null truth and the 15 MeV structure excepted).
- **Comparison with B4 and with the per-slice edge** (the same scenario definitions; bedge numbers from bedge_toys):

  | truth change | B4 | per-slice M4 edge | u edge |
  |---|---|---|---|
  | quadratic term q = -0.173 | +150 % | -2.1 % | -2.2 % |
  | q = -0.346 / +0.173 | -- | -4.1 / +2.2 % | -4.5 / +2.2 % |
  | omega/rho0 re-mix | -- | -1.7 % | +0.5 % |
  | 10 % at 15 MeV / 5 % at 8 MeV | -- | +20 / +9 % | +17 / +8 % |

- **Reading.**
  - Both edge methods are equally insensitive to the smooth deformation that destroyed B4: the shift is about 70 times smaller.
  - Both are equally sensitive to narrow structures and to a resolution mismatch.
  - The raw-background truth (+10 %) is the template noise of the simulated background placed into the truth (its step
    projection scatters by +-200 rows, 6 % of h_sim, from replica to replica). It does not enter the estimator, because h_sim
    uses the bump classes only. Nothing of this size is seen in the data control region (3 +- 17 rows).

---------------------------------------------------------------------------------------------------
## 6. Conversions and comparison (Fig. edgeu_compare; Table tables_all/edgeu_results.{txt,tex})

- **Visible bump yield** = f x (c C2 + C3)_SMASH x 64/63. The MC error comes from the same replicas (f_r x W_r):
  - window B, all M(ee): **21 443 +- 1 002 +- 1 559 +- 1 329 rows** = 3.53 nb (SMASH c C2 + C3: 51 048.8; C2 + C3 =
    54 758.9, the same as WP6/bedge);
  - all M4: 26 955 +- 1 260 +- 1 960 +- 1 690 rows = 4.43 nb;
  - B in E2: 9 552 (global f); 9 848 +- 699 (stat) +- 928 (MC) with the own E2 fit.
- **sigma(pp -> eta X) = f / 1.00436 x 1.0694 mb x 1.0104 x 0.97:**

  | reading | sigma [mb] |
  |---|---|
  | **u edge, C2 corrected** | **0.438 +- 0.020 +- 0.032 +- 0.030 +0.042/-0.044** |
  | u edge, uncorrected | 0.407 +- 0.019 +- 0.028 +- 0.028 +0.039/-0.041 |
  | u edge, C2 corrected, conversions x 0.48 (E9 +-0.033 on the scale included) | 0.563 +- 0.026 +- 0.053 +- 0.038 +0.054/-0.056 |
  | per-slice edge (bedge), C2 corrected / uncorrected / x 0.48 | 0.437 / 0.409 / 0.585 |
  | peak side bands (adopted, PRELIMINARY) | 0.429 +- 0.036 +- 0.152 +- 0.025 |
  | direct-decay route / conversion route | 0.522 +- 0.062 / 0.249 +- 0.051 |
  | B4 (not determined; scale 0.312) | 0.326 |
  | F-main (cross-check) | 0.399 |
  | ECAL tag (verification only; scale 0.454 +- 0.043) | 0.474 |
  | SMASH | 1.069 |

- **Double ratio peak/bump:** DR_u = a_P / f (PDG 2026) = 0.978 +- 0.094 (stat) +- 0.172 (syst) +- 0.090 (MC), with the full
  a_P errors of WP7 (bedge: 0.982).
- **Plain answer: is u better than the per-slice edge?**
  - **Smaller total error: yes, 11.1 % against 15.0 %** (stat 4.7 against 5.4 %, syst 7.3 against 8.3 %, MC 7.0 against 11.3 %).
  - **Why:**
    - (i) every row is measured against its own edge, so the step is sharper and the stat error 13 % smaller;
    - (ii) h_sim comes from the bump classes alone, which removes the template noise of the raw simulated background (the
      bedge caveat 1);
    - (iii) the slice choices disappear (bedge E2 2.9 % -> U3 0 %).
    - The single fit itself adds nothing: four separate fits in u give the same numbers.
  - **Less model dependence: no.** The two leading systematics are the same in both methods:
    - the data/SMASH resolution ratio at the edge;
    - the odd part of the smooth background across the edge.
    The M(ee) dependence of the edge is now treated exactly instead of being averaged.
  - Both edge methods are a factor of about 70 less sensitive than B4 to a smooth background deformation.

---------------------------------------------------------------------------------------------------
## 7. Figures and tables (draft captions)

**Fig. `figs_all/edgeu_kinematics.pdf`.**
- **What it shows: the edge of the bump depends on M(pi+pi-) and M(e+e-), and once both are taken into account no simulated
  bump event lies above it; the reconstructed u is smeared by 6.7 MeV/c2.**
- (a) The edge E(M(pipi), M(ee)) on the data scale for M(ee) = 0, 20, 50, 100 MeV/c2; green: the eta peak 0.5439; dotted: the
  bedge slice boundaries; no bump above M(pipi) = 0.4116 GeV/c2.
- (b) Simulated truth u_true = 0.99274 [m4_true - M4max(M(pipi)_true, M(ee)_true)] for C2 (orange) and C3 (magenta); dashed:
  C2 with the edge at M(ee) = 0 (rows above u = 0).
- (c) Reconstructed u of simulated C2 with M(ee) >= 20 MeV/c2, with the M(ee) term (orange) and without it (dashed).
- (d) u_reco - u_true near the edge (C2 + C3, u_true > -15 MeV/c2) with the Gaussian core (6.70 MeV/c2, nominal) and the
  data-like width x 0.893.
- SMASH mode A, w x V3, S0, M4' = 0.99615 M4, M(pipi)' = 0.9952 M(pipi), expected raw rows at L = 6079.2 nb^-1.

**Fig. `figs_all/edgeu_spectra.pdf`.**
- **What it shows: in u the data show one edge at u = 0 and the peak at u = +20 to +55 MeV/c2; in SMASH the bump classes end at
  u = 0 and no other class has an edge there.**
- (a) Data 0.30 <= M(pipi) < 0.40 GeV/c2: N - CB for S0 (black) and pre (blue squares), CB = 2 sqrt(N++ N--) (red), N-- (red
  dotted) and N++ (red dash-dotted), per 2.5 MeV/c2; grey band: fit range.
- (b) SMASH (expected raw rows): C2 (x 0.88 gamma* factor) orange, C3 magenta, C4 green, C5 cyan, G_eta violet, G_D dark grey,
  G_C light grey (dashed: omega/rho0 pions).

**Fig. `figs_all/edgeu_fit.pdf`.**
- **What it shows: one fit with a smooth background and a step at u = 0 describes the data (chi2 27/32); the step is 21 sigma
  and has the shape of the SMASH bump scaled by f = 0.420.**
- Top: data S0 N - CB - a_P K (blue; K = SMASH peak C4 + C5 in u, a_P = 0.418, green area); fit pol2 + h x step (green) and
  without the step (dashed).
- Bottom: data minus the fitted smooth part (blue); SMASH bump C2 + C3 minus its own smooth part (same estimator), x f (orange);
  h_data x step (green, smearing 6.7 MeV/c2).

**Fig. `figs_all/edgeu_checks.pdf`.**
- **What it shows: the slices, the M(e+e-) classes, the periods and the selections agree within their errors; there is no step
  in the backgrounds, except in the e-e- like-sign rows.**
- (a) f per M(pipi) slice: u, C2 corrected (black; thick: stat, thin: stat (+) MC), u uncorrected (blue), bedge per-slice M4
  ratios (brown, data errors only).
- (b) f per M(ee) class with bedge's C2-corrected E1/E2; kappa from the E1/E2 split.
- (c) Periods (x C_sec), folds, pre, isBest == 1 (stat).
- (d) Null tests: step height at u = 0 of the like-sign rows (e+e+, e-e-, e-e- by class), of the subtracted CB, of SMASH G_D,
  G_C, G_eta (bars: replica rms) and of the data control region M(pipi) 0.42-0.46 with a pseudo edge; grey: +-5 % of h_data.

**Fig. `figs_all/edgeu_syst.pdf`.**
- **What it shows: only the resolution variants move f by more than their noise; the excluded Dalitz-plot variants, q = 0 and
  isBest == 1 are listed as crosses.**
- Shift of every variant in per cent, with the bootstrap (+) replica noise of the shift; colours = families U1-U10.

**Fig. `figs_all/edgeu_toys.pdf`.**
- **What it shows: the u-edge estimator closes within 2 % for smooth backgrounds; the quadratic term that moves B4 by +150 %
  moves it by -2 %; only a resolution mismatch and background structures as narrow as the edge bias it.**
- (a) Bias of f for 21 truths, 1000 toys each: AA (blue), RF (orange squares); brown diamonds: the per-slice M4 edge (bedge)
  for the same truth change, placed relative to (a); grey band +-5 %.
- (b) Pull widths.

**Fig. `figs_all/edgeu_compare.pdf`.**
- **What it shows: sigma(pp -> eta X) from the bump in u (0.438 mb) agrees with the per-slice edge and with the peak side
  bands; all routes lie a factor 2-2.5 below SMASH.**
- PDG 2026, C_sec (1 - f_eta) = 1.0104 x 0.97 (bump); thick: stat; thin: stat (+) syst (+) MC; normalisation not included;
  open circles: values without errors (F-main cross-check; B4 not determined).

**Tables** (`docs/latex/tables_all/`; the three .tex files compile; tested in the scratchpad at \tiny, no error):
- `edgeu_results.{txt,tex}`: f, yields, sigma (three readings), DR, E1/E2 split, comparison;
- `edgeu_syst.{txt,tex}`: every variant with value, shift, noise (literal, bootstrap), counted delta, status; families;
- `edgeu_toys.{txt,tex}`: closure.

---------------------------------------------------------------------------------------------------
## 8. Caveats

1. **The e-e- step in the like-sign rows** (Sect. 2) is not understood.
   - If it is an artefact of bump events (for example the Dalitz e+ reconstructed with a negative charge), the CB
     over-subtracts. The migrated pairs would also be missing from the unlike-sign sample. f would then be higher by up to
     3 % (CB) + 5 % (migration, if not simulated).
   - If it is genuine combinatorics, the nominal is right.
   - U5 (2.6 %) covers the CB part only. **To be checked by the analyst** (lepton charge and track quality of the e-e- rows
     near u = 0). The same applies to bedge and to B4.
2. **Simulation statistics.** Only 1 649 (C2) and 1 187 (C3) simulated rows, 645 within 15 MeV/c2 of the edge. MC (7 %) is
   still larger than stat (4.7 %).
3. **U11** rests on the assumed background shapes (smoothed SMASH, stage A). The data cannot fix the cubic part across the edge
   better than the pol3 fit does (f 0.418 +- 0.036).
4. **U2:** the data/SMASH resolution ratio is known to +-8 % (peak width and step width combined). The variants use +-10 %.
5. **The C2 gamma* correction** is applied (decision D35). The uncorrected reading is given everywhere.
6. **Folds** z = -2.1 (p = 0.03) and **periods** p = 0.025 are both marginal. They point to the same data feature as in bedge,
   and the u fit does not remove it.
7. **The h_sim replica mean** is 1 % below the nominal (3078 vs 3107; the WP0 replica set; the verifier found the same set 1.3
   sigma low). Information only.
8. **Rows with M(ee) >= m_pi0** have no defined u and are not used. No bump row is lost; the bedge row set (all M(ee),
   E(m, 0)) gives the same as q = 0.
9. **Blinding:** only the M(ee) classes < 10, >= 10 and >= 20 MeV (F10/F5g edges) were formed in the data. M(ee) enters u
   only through E(m, q) row by row; no X2D cell or fine M(ee) shape was made.
10. **The absolute sigma** inherits the PRELIMINARY status of Part III. The one-sided items of the per-slice note are not
    included: live time +7.3 %, sorter +1.3 to +2.8 %, cor_ok +2.2 %, C_vtx, production model.

---------------------------------------------------------------------------------------------------
## 9. Exact commands (from /home/przygoda/HADES/FAT/FAT)

    B='gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_edgeu", kTRUE)'
    O=/hdd3/przygoda/HADES/FAT_out/all_report/edgeu
    root -l -b -q -e "$B" 'plots/all_report/edgeu_dump.C+O'    > $O/edgeu_dump_stdout.txt 2>&1     # 35 s, about 10 GB: column files
    root -l -b -q -e "$B" 'plots/all_report/edgeu_dumpls.C+O'  > $O/edgeu_dumpls_stdout.txt 2>&1   # 5 s: simulated like-sign rows (null test)
    python3 -I plots/all_report/edgeu_ana.py     > $O/edgeu_ana_stdout.txt 2>&1                    # 70 s: edgeu_ana.json, edgeu_ana.log, figdata
    python3 -I plots/all_report/edgeu_bgsmooth.py > $O/edgeu_bgsmooth.txt 2>&1                     # 20 s: smoother chi2 and step projections
    python3 -I plots/all_report/edgeu_toys.py 1000 > $O/edgeu_toys_stdout.txt 2>&1                 # 45 s: edgeu_toys.json, edgeu_toys.log
    python3 -I plots/all_report/edgeu_combine.py                                                   # tables_all/edgeu_*.{txt,tex}, edgeu_numbers.json
    python3 -I plots/all_report/edgeu_figs.py                                                      # figs_all/edgeu_*.{pdf,png}

**Code** (`plots/all_report/`):
- `edgeu_core.py`: columns, the edge function (= defs.h Edge, vectorised), u, histograms, the CB rule in u, step templates, the
  fit (the logic of scan::FitEdgeStep), and the WP0 replica generator (reproduces bedge::GroupPois exactly; checked on nine
  groups) with its data bootstrap.
- `edgeu_dump.C`: data and simulated e+e- rows. `edgeu_dumpls.C`: simulated like-sign rows.
- `edgeu_ana.py`, `edgeu_bgsmooth.py`, `edgeu_toys.py`, `edgeu_combine.py`, `edgeu_figs.py`.

**Typed inputs:**
- a_P = 0.418249 +- 0.035 +- 0.067 +- 0.025 (WP7/WP8);
- C2 factors 0.879 / 0.8825 (bedge recommended means of WP8) +- 0.024;
- 1.0694 mb, 1.0104, 0.97, 1.00436, the normalisation parts, kappa 0.48 +- 0.11, the comparison values (bedge_combine.py);
- L_p and C_sec(p) (bedge_ana.C);
- the data/SMASH peak width 5.44/6.09 MeV/c2.

---------------------------------------------------------------------------------------------------
## 10. Verification (independent verifier, prefix vedgeu, 2026-11-15)

### 10.0 In plain words

**What was done.**
- I rebuilt the u measurement with my own code. No edgeu macro, column file or JSON was used for any number below.
  - **Own row pass** (`vedgeu_fill.C`): own loops over the data compact trees (us, pp, mm) and over the simulated trees
    (std_us, lep_us, eta_us). The selection code (pre, S0, ownership A) is written from defs.h Sect. D.
  - **Own analysis** (`vedgeu_core.py`, numpy):
    - edge function and u;
    - CB rule in u;
    - step template (own numerical integration);
    - iterated GLS fit;
    - bootstrap and replica generators (own seeds; independent of the WP0 generator).
  - **Shared by design:** the frozen definitions (alldefs::CompsOf, simphys::RowFactor 108, 27.0718) and the typed Part III
    inputs (a_P = 0.418249, the C2 factors 0.879 / 0.8825, 1.0694 mb, 1.0104, 0.97, 1.00436). The step width 6.70 MeV/c2 was
    also taken from edgeu, so the numbers can be compared exactly; my own width estimates are run as variants.
- **What I recomputed:**
  - the edge function: at seven points; by a brute-force decay generator; against the simulated truth;
  - the u spectra of the data and of SMASH;
  - the global step fit (h_data, h_sim), f, its bootstrap and replica errors;
  - nine variants;
  - the null tests and the e-e- step;
  - the classes, the four slices and the periods;
  - the yields and sigma;
  - pseudo-data with the B4-killing quadratic term (1000 toys);
  - a new toy for the C2 gamma* angle factor at the u edge.
- Then I tried to refute the plain-words conclusions, and I rendered all seven figures.

**What came out.**
- **Every decisive number is reproduced, almost to the last digit:**
  - h_data = 1305.05 +- 62.27 (edgeu 1305.2 +- 62.3);
  - chi2 27.34/32, without the step 466.6;
  - h_sim = 3106.95 (3107.3);
  - **f = 0.4200 +- 0.0200 (fit)**, uncorrected 0.3899, kappa 0.48 reading 0.5400.
- **The errors agree:**
  - stat (my bootstrap) 0.0202, against 0.0196;
  - MC (my independent replicas) 0.0259, against 0.0292 (their robust value is 0.0253). The quoted MC is on the safe side.
- **The edge kinematics are right:**
  - a brute-force scan of the photon direction reproduces M4max(m, q) to 1e-6 GeV/c2;
  - no simulated C2 or C3 row lies above its true edge (largest excess -0.476 / -0.508 MeV/c2);
  - without the M(ee) term, 3.83 % of C2 lie above it by more than 1 MeV/c2.
- **The variants I chose agree exactly:**
  - CB = N++ + N--: -2.60 %;
  - r = 0.89 / 1.10: -5.55 / +3.28 % (pair 4.4 %);
  - own smearing per sample: -4.74 %;
  - q = 0: +5.85 %;
  - pol3: -0.55 %;
  - T12: -5.87 %.
- **The pseudo-data claim holds.** I built my own truth and estimator. The quadratic term that moves B4 by +150 % moves the u
  edge by -2.0 % (Asimov) or -1.9 % (1000 toys); edgeu has -2.2 %. The narrow structures give +7.0 % and +16.6 % (edgeu +7.6 %
  and +17.3 %). The pull widths are 0.99-1.02.
- **New check: the C2 gamma* factor transfers to u.** In a toy, the QED/isotropic ratio of the u-edge step is 0.902 (E1) and
  0.868 (E2). With the SMASH E1/E2 mix of the C2 step this gives 0.879, while edgeu applies 0.881. The effect on f is +0.14 %,
  far inside U10 (1.3-1.6 %).

**What it means.**
- **Verdict: CONFIRMED.** The u-edge result f = 0.420 and sigma = 0.438 mb are reproduced with independent code. The quoted
  errors are adequate.
- **Four statements need different wording.** None of them changes a number.
  1. **"11 % instead of 15 %" is not a like-for-like comparison of the two variables.**
     - The MC gain comes from computing h_sim from the bump classes only. The per-slice edge can use the same rule.
     - With that rule, its MC would be about 6.6 % (computed from the vbedge replica decomposition), not 11.3 %, and its total
       about 11.9 %.
     - So the gain that belongs to u itself is about 11.9 % -> 11.1 %. It comes from a stat error about 8 % smaller (not 13 %)
       and from the slice choices, which disappear.
  2. **The double ratio is quoted with too small an error.**
     - "The full a_P errors of WP7" is not what the code uses. It takes the a_P systematic without the composition term (16 %).
     - The adopted peak side-band sigma carries 35.5 % systematic. With it, DR_u = 0.98 +- 0.09 +- 0.35 +- 0.09 (total 0.37,
       not 0.22).
     - "Agrees with the peak side bands" is therefore a weak test.
  3. **The e-e- step is confirmed (73.5 +- 15.7 rows, 4.7 sigma), but two of the descriptions are not significant.**
     - "Only for M(ee) >= 10 MeV" mostly follows from where the e-e- rows are: 91.5 % of them are in E2. A step spread evenly
       over the rows would put 6 of the 74 rows in E1; the observed E1 value is -6 +- 5 (2.6 sigma below).
     - "Mostly in period I" is not significant: a step that follows L has chi2 4.9/3, p = 0.18.
     - It is large: about 30 % of the e-e- rows near u = 0, against about 36 % for the bump in the unlike-sign rows. That
       points to bump events in the like-sign sample, so caveat 1 (up to about +8 % on f) deserves the analyst's attention.
  4. **G_eta has "no step" only after smoothing.** A smoothed shape has no step by construction. The raw value, 123 +- 53 rows
     (my replicas, 2.3 sigma), has no visible structure; 644 heavily weighted std rows carry half of it. Treating it as noise is
     reasonable. If it were real, f would be 3.8 % lower, which is within U11.

### 10.1 Agreements (own code vs edgeu)

| quantity | own (vedgeu) | edgeu | verdict |
|---|---|---|---|
| E(0.30,0), E(0.39,0), E(0.39,60), E(0.4116,0), E(0.30,100) | 0.51859, 0.50077, 0.50958, 0.47390, 0.53262 | same | exact |
| brute-force max M4 (2e5 photon directions) vs M4max at (0.30,0) (0.35,30) (0.39,60) (0.40,100 MeV) | differences -8.6e-7, -2.1e-8, -7.8e-7, -3.2e-7 GeV | -- | formula exact |
| truth: max(m4_t - M4max(mpp_t, mee_t)) C2 / C3 | -0.476 / -0.508 MeV (1514 / 1108 S0 rows), 0 rows above | -0.48 / -0.51 (S0 + pre rows) | exact |
| truth with q = 0: C2 above by > 1 / > 5 MeV | 3.83 % / 1.80 % | 3.8 / 1.8 % | exact |
| truth, tightness: largest m4_t in (M(pipi)_t, M(ee)_t) cells | 0.8-1.2 MeV below the cell limit where the cell has > 100 rows | -- | the limit is reached |
| reco u above +10 / +20 MeV (C2; C3) | 0.85 / 0.08 %; 1.04 / 0.16 % | same | exact |
| u resolution: rows, half 68 %, rms | 645, 7.78, 10.22 MeV | 645, 7.78, 10.2 | exact |
| Gaussian core | 7.10 +- 0.26 (+-40 MeV, own error floor); 5.76 +- 0.31 (iterated +-2 sigma) | 6.70 +- 0.21 | method-dependent (10.2 item 4); f 0.4204 / 0.4188 |
| SMASH C2, C3 in the dump range (S0) | 38 688.6, 29 075.4 expected raw rows | same | exact |
| data S0 rows in 0.30-0.40 and the u grid (us, pp, mm) | 398 814, 4 187, 21 024 | -- | -- |
| h_data, chi2, chi2 without step | 1305.05 +- 62.27, 27.34/32, 466.6 | 1305.2 +- 62.3, 27.3/32, 466.6 | exact (0.01 %) |
| h_sim; linear split C2c / C3 / C2 raw | 3106.95; 1780.1 / 1326.9 / 2019.9 | 3107.3; 1780.3 / 1327.0 / 2020.1 | exact |
| f; uncorrected; kappa 0.48 | 0.42004; 0.3899; 0.5400 | 0.42004; 0.3900; 0.5400 | exact |
| peak term: projection of K on the step; a_P free | 3.86 rows per unit a_P (a_P +- 0.08 -> 0.31 rows, 0.02 %); a_P 0.329, f 0.4200 | 0.02 %; 0.33, 0.4200 | exact |
| free step width on the data | 7.02 +- 0.79 MeV | 7.02 +- 0.79 | exact |
| range / order scan | pol1 115.6/33; pol2 p 0.70 (+-35) ... 0.47 (+-55), 0.15 (+-60), 0.03 (+-70); f 0.423-0.429 for +-35..+-55; pol3 f 0.4177, error x 1.82 | same | exact |
| stat (bootstrap) | 0.0202 (200, own seed; bootstrap/fit 1.01) | 0.0196 | agree |
| MC (replicas) | 0.0259 (100, own generator); h_sim rms 192 (6.2 %), mean 3117.9 (+0.35 %) | 0.0292 (robust 0.0253); 201 (6.5 %), 3077.7 (-1 %) | agree; confirms caveat 7 (WP0 set low) |
| variants: CB = N++ + N--, r 0.89, r 1.10, own smearing, q = 0, pol3, CB 5 / 20 MeV cells, T12 | -2.60, -5.55, +3.28, -4.74, +5.85, -0.55, +0.43, +0.36, -5.87 % | -2.6, -5.5, +3.3, -4.7, +5.8, -0.55, +0.43, +0.36, -5.9 % | exact |
| bootstrap noise of the CB = N++ + N-- shift | 0.0040 | 0.0037 | agree |
| null tests e+e+, e-e-, sum | 2.1 +- 6.5, 73.5 +- 15.7, 75.7 +- 17.0 (own fit); linear 2.2, 74.3, 76.6 | 2.1, 73.5, 75.7 | exact |
| e-e- E1 / E2; periods I, IIA, IIB, III; pre-not-S0 | -5.6 +- 4.6 / 79.9 +- 15.1; 51.4, 1.2, 6.9, 14.7; 20.8 +- 18.4 | -5.4 / 78.7; 50.8, 1.0, 7.3, 14.7; 23 +- 18 | agree |
| step of the subtracted CB | 39.2 rows (3.0 %) | 39.2 | exact |
| SMASH G_D, G_C, G_eta step | 12.4 +- 102, 58.5 +- 180, 122.7 +- 53 (own replicas) | 12 +- 140, 59 +- 148, 123 +- 57 | agree |
| classes E1 / E2 (CB share rule); kappa split | 0.4089 / 0.4331; a_D 0.444, kappa 0.876 +- 0.19 (linear stat) | 0.409 / 0.433; 0.444, 0.88 +- 0.21 | exact |
| four slices: h_data; h_sim; f (sum / sum) | 325.4, 436.7, 387.0, 145.4; 701.8, 1077.4, 1031.9, 245.6; 0.4235 | same; 0.4235 | exact |
| periods (x C_sec) | 0.433, 0.391, 0.531, 0.356; chi2 8.75/3 (fit errors), p 0.033 | same; 9.4/3, p 0.025 (bootstrap errors) | agree |
| yields: (C2 + C3), (c C2 + C3) in B x 64/63; Y_B; Y all M4 | 54 758.9, 51 048.8; 21 443 rows = 3.527 nb; 26 954 | same; 21 443, 3.53 nb; 26 955 | exact |
| sigma C2 corrected / uncorrected / kappa 0.48 | 0.4383 / 0.4069 / 0.5635 mb | 0.438 / 0.407 / 0.563 | exact |
| pseudo-data (a): Asimov / toy bias, pull SD | +1.81 % / +1.73 % (1000 toys), 1.02 (kernel smoother: +1.70 / +1.75 %, 1.01) | +2.0 % / +1.5 %, 0.99 | agree |
| q = -0.173 shift vs (a) | -2.03 % (Asimov), -1.93 % (toys); kernel smoother -1.78 / -1.87 % | -2.2 % | agree |
| q = -0.346 / +0.173; 5 % G(505, 8 MeV); 10 % G(490, 15 MeV) | -4.1 / +2.0 %; +7.0 %; +16.6 % (Asimov) | -4.5 / +2.2 %; +7.6 %; +17.3 % | agree |
| step projections of the smoothed non-eta background, exp(pol2/3/5/6) | -35.6, -7.8, +9.3, +17.3 rows | -35.6, -7.8, +9.3, +17.3 | exact |
| C2 gamma* factor in u (new toy) | R_u E1 0.902, E2 0.868, all 0.887 (R_all 0.889 / 0.868 / 0.880); with the SMASH C2 step mix (E2 share 0.656): 0.879 against the applied 0.881 | not derived for u | f +0.14 %: inside U10 |
| tables | the three .tex compile (pdflatex, \tiny, landscape A4) | stated | confirmed |
| units and conventions | rows per 2.5 MeV/c2 bin for data and simulation; f a ratio of two such heights; sim in expected raw rows (27.0718); yields x 64/63; M4' = 0.99615 M4, M(pipi)' = 0.9952 M(pipi), edge 0.99274 M4max(m/0.9970, q) | same | correct |

### 10.2 Discrepancies and wording problems (with sizes)

1. **"Is u better? total error 11.1 % against 15.0 %" (sect. 0, 6; the package summary).**
   - The MC part of the gain (11.3 % -> 7.0 %) comes from computing h_sim from the bump classes only. The note says this
     choice "could also be made in the per-slice method", but the headline comparison still mixes the two conventions.
   - With the vbedge replica decomposition, the bump-only part of the per-slice E_sim has an rms of 189-194 rows (6.5-6.7 %
     of 2884). So the per-slice edge with bump-only templates would have MC about 6.6 % and a total of about 11.9 %.
   - **Like-for-like, u brings 11.9 % -> 11.1 % (about 7 % smaller):**
     - stat about 8 % smaller: fit errors 4.77 % against 5.19 %; with E(m, 0) in u it is 5.02 %, so the M(ee) alignment gives
       about half of it;
     - the slice-choice family (2.9 %) disappears.
   - **"stat error 13 % smaller" should read "about 8 % smaller".** The 13 % compares two bootstraps, whose ratios to their fits
     are 0.98 and 1.03 (noise of a 200-copy rms, about 5 %).
2. **DR_u (sect. 6, sect. 0 "What it means", Table edgeu_results).**
   - edgeu_combine.py line 171 uses AP_SYST/AP = 0.067/0.418 = 16 %. That is the a_P systematic without the peak-composition
     term, as in bedge, and not "the full a_P errors of WP7".
   - The adopted peak side-band sigma (0.429 +- 0.036 +- 0.152 +- 0.025 mb, writer log) carries 35.5 % systematic. With it,
     **DR_u = 0.98 +- 0.09 +- 0.35 +- 0.09, total +- 0.37 (38 %), not +- 0.22.**
   - With kappa = 0.48 in both readings the ratio is about 0.89 (vbedge item 4).
   - Say "consistent within the 38 % uncertainty of the double ratio" instead of "agrees".
3. **The e-e- step (sect. 0, 2, caveat 1, "Things to flag" 1).** The number is confirmed: 73.5 +- 15.7, 4.7 sigma, own fit
   chi2 37.1/32 with the step and 59.1 without. Two descriptions overstate what it shows:
   - **"Only for M(ee) >= 10 MeV":** 91.5 % of the e-e- rows in -45..+45 MeV/c2 are in E2. A step spread evenly over the rows
     would put 6.3 rows in E1; the observed value is -5.6 +- 4.6 (2.6 sigma below). The E2 location is mostly the row
     distribution.
   - **"Mostly in period I":** I 51.4 +- 10.8, IIA 1.2 +- 3.9, IIB 6.9 +- 7.0, III 14.7 +- 8.3. Against a share proportional to
     L (34.3, 4.6, 13.2, 22.2) the chi2 is 4.9/3 (p 0.18); proportional to the e-e- rows it is 4.6/3 (p 0.21). Not significant.
   - **What is significant is its size:** about 30 % of the e-e- rows per bin near u = 0 (step 74 on about 255 rows per bin),
     against about 36 % for the bump in the unlike-sign rows (1305 on about 3600). Combinatorial like-sign pairs have no reason
     to step at the eta edge. So the bump-related reading (a wrong charge sign, f up to about +8 %) is a real possibility. The
     recommendation to the analyst stands.
4. **Gaussian core width (sect. 1.3, Fig. edgeu_kinematics (d)).**
   - "sigma = 6.70 +- 0.21 MeV/c2" depends on the method: my +-40 MeV fit with another error floor gives 7.10 +- 0.26; an
     iterated +-2 sigma core gives 5.76 +- 0.31.
   - It does not matter for f: the same template enters data and SMASH, and f is 0.4204 / 0.4188 with these widths.
   - Quote it as "the template width (Gaussian fit on +-40 MeV/c2)" rather than as a measured resolution.
5. **G_eta "smoothed shapes: -5 to -31, no step" (sect. 2 table).**
   - A smooth exp(pol n) shape has no step by construction, so this does not test the raw value.
   - The raw value, 122.7 +- 53.3 (my independent replicas), is 2.3 sigma.
   - The G_eta spectrum near u = 0 shows no structure: it rises smoothly through u = 0. Its noise is dominated by 644 std rows
     with weights of about 24.
   - Treating it as template noise is reasonable. Say so explicitly, and say that if it were real, f would be 3.8 % lower,
     which U11 covers in size.
6. **Toys (Table edgeu_toys).**
   - For truth (a), the AA toy mean (+1.48 %) is 0.56 % below edgeu's own Asimov value (+2.04 %). That is 3.7 times the
     standard error of a 1000-toy mean.
   - My toys agree with my Asimov to 0.1 % (+1.73 / +1.81 %). The other truths differ by 0.1-0.3 %.
   - No conclusion changes. The Asimov column is the cleaner number to quote.
7. **C2 gamma* factor in u (sect. 3, caveat 5).**
   - edgeu applies the WP8 all-M4 factors without checking them at the u edge.
   - My toy (WP8 generator, own u-edge estimator) gives R_u = 0.902 (E1) and 0.868 (E2); with the SMASH C2 mix this is 0.879,
     while 0.881 is applied (f +0.14 %).
   - Add one sentence; no number changes.
8. **Minor.**
   - Fig. edgeu_syst: U11 (pseudo-data, 4.0 %) is in the title but not among the rows. Add "U11 not shown (from the
     pseudo-data)" to the caption.
   - Fig. edgeu_fit, top: the dashed line is the smooth part (pol2) of the same fit, not a refit without the step. The legend
     "the same without the step" is right; the caption should say "the smooth part".
   - The pol2 of the fit also absorbs the slope of the bump below the edge (SMASH C2c + C3 rises from 1370 to 2590 rows per bin
     between -44 and -14 MeV/c2). Only the step height is a physical quantity. One sentence would prevent a misreading of the
     dashed line as "the background".

### 10.3 Figures and tables

- All seven PNG files were rendered and read. Fonts are readable, the colours follow Part III, and the axes carry units
  (MeV/c2, GeV/c2, mb).
- **Numbers in the figures agree with mine:**
  - edgeu_fit: h_data, chi2, delta chi2;
  - edgeu_checks: slices, classes, periods, null tests;
  - edgeu_toys: the q-term shift;
  - edgeu_compare: the sigma values.
- **Remarks:** the edgeu_syst caption (U11) and the edgeu_fit dashed line, see 10.2 item 8. Panel (a) of edgeu_spectra is on
  a log scale, where the 30 % dip at u = 0 is hard to see; edgeu_fit shows it clearly.
- **Tables:** edgeu_results / edgeu_syst / edgeu_toys (.txt and .tex) agree with the note. The DR row inherits 10.2 item 2.
  The three .tex files compile (pdflatex, scratchpad).

### 10.4 Own code, outputs, exact commands (from /home/przygoda/HADES/FAT/FAT)

    O=/hdd3/przygoda/HADES/FAT_out/all_report/vedgeu
    root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/vedgeu/aclic", kTRUE)' 'plots/all_report/vedgeu_fill.C+O'   # 25 s, 1 GB: own row columns
    python3 -I plots/all_report/vedgeu_ana.py    > $O/vedgeu_ana_stdout.txt      # 6 s: kinematics, fit, bootstrap, replicas, variants, null tests, yields
    python3 -I plots/all_report/vedgeu_toys.py   > $O/vedgeu_toys_stdout.txt     # 9 s: pseudo-data (a), q-terms, structures, 1000 toys
    python3 -I plots/all_report/vedgeu_checks.py > $O/vedgeu_checks_stdout.txt   # stat of E(m,0) and pol3, classes, kappa, four slices, periods, e-e-
    python3 -I plots/all_report/vedgeu_angle.py  > $O/vedgeu_angle_stdout.txt    # 22 s, 4.4 GB: C2 gamma* factor at the u edge (toy, 1.2e7 decays)

- **Code** (`plots/all_report/`):
  - `vedgeu_fill.C`: own row pass. Its column layout is in the header.
  - `vedgeu_core.py`: own edge, u, CB, step template, GLS fit, Poisson replica weights.
  - `vedgeu_ana.py`, `vedgeu_toys.py`, `vedgeu_checks.py`.
  - `vedgeu_angle.py`: imports only the WP8 event generator (load_gen, generate, accept) from bedge_angle_toy.py.
- **Scratch:**
  - /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/vedgeu/geta_look.py (the
    G_eta spectrum near u = 0);
  - /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/vedgeu/scan.py (range and
    order scan);
  - /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/vedgeu/tabtest.tex (table
    compilation).
- **Outputs** in /hdd3/przygoda/HADES/FAT_out/all_report/vedgeu/: vedgeu_ana.json, vedgeu_toys.json, vedgeu_checks.json,
  vedgeu_angle.json, vedgeu_spec.npz, the *_stdout.txt files and vedgeu_fill.log.
- **Housekeeping:**
  - No existing file was edited; only new vedgeu_* files were written, plus this appended section. No git command was run.
  - One process ran at a time, and the CPU package stayed at 38-57 C.

### 10.5 Limits of this verification

- **Not redone:**
  - the 70-variant family combination: nine variants were reproduced exactly, and the U2 pair rule was checked by hand;
  - the RF toys and the 21-truth table: I reproduced (a), the three q-terms and the two narrow structures;
  - folds, isBest == 1 and pre;
  - the E2-only yield and the bgsmooth chi2 values.
- **Shared with edgeu by design:**
  - the frozen definitions (RowFactor 108, CompsOf);
  - the typed Part III inputs (a_P, the C2 factors and their +-0.024, kappa 0.48, C_sec, 1.0694 mb, the normalisation);
  - the template width 6.70 MeV/c2, used for the exact comparison (own widths as variants);
  - the WP8 event generator used in the angle toy.
