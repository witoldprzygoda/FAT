# peaktail: the low-mass tail of the eta peak and the bump edge -- results note

> **STATUS (final check, 15 November 2026, 17:50): independently verified; read the section "Verification (independent)" at the
> end first; the copy in docs/latex/notes_all/ is identical to this file.** The physics and the arithmetic of this note are confirmed. **Withdrawn** by the verification: "t_data = 0.72 +- 0.38",
> "the data have a tail of the simulated size", "no tail disfavoured at 1.9 sigma", "0.384 excluded by the data at 2.6 sigma",
> "the open item is closed", "the template noise is included in the t errors" (true for the u fit only) and the proposed terms
> +-2.6 / +-2.8 / +-6.0 %. The reports use the verified reading: the data allow t = 0.3 +- 1.0; proposed D40 with t = 1 +- 0.5
> (a physics band): +-3.4 % on f_eta(u-edge), +-3.6 % per slice, +-8.8 % on a_P and sigma(pp -> eta X) from the peak side bands.
> The figure peaktail_constraint.pdf and the grey t_data band of peaktail_effect.pdf are not used in any report.

Date 2026-11-15 (after the integration of 15:45). Branch pp45_pippimepem (worktree /home/przygoda/HADES/FAT/FAT). No git command
was run. No existing file was edited, moved or removed: the reports, defs.h, edgeu_core.py, edgeu_ana.py, bedge_ana.C,
bedge_common.h, the edgeu/bedge outputs and the compact files were only read (bedge_ana.C is #included unchanged by
peaktail_bedge.C). All new files carry the prefix peaktail. At most 2 processes ran at a time (one ROOT dump of about 12 GB
for 35 s, otherwise numpy); CPU package 38-47 C. **Not independently verified.**

**Machine-readable outputs** (`/hdd3/przygoda/HADES/FAT_out/all_report/peaktail/`): `peaktail_numbers.json` (all numbers below),
`peaktail_{sim,ufit,line,rad,prop}.json`, the readable logs `peaktail_{sim,ufit,line,rad,prop,summary}.log`,
`peaktail_bedge.txt`.

---------------------------------------------------------------------------------------------------
## 0. In plain words

**Key results.**
1. **The tail of the simulated eta peak is made by the leptons losing energy in the detector material (bremsstrahlung).**
   In the generator the four-body mass is a sharp line: there is no radiative photon in the simulated decays. About a quarter
   of the simulated peak rows (C4 26 %, C5 24 % for 0.30 <= M(pi+pi-) < 0.40) have leptons that lost enough momentum to shift
   M4 down by more than 5 MeV/c2. These rows form the tail.
2. **The tail reaches the bump edge.** The peak sits at u = +21 to +46 MeV/c2. Its tail rows fill the region between the edge and
   the peak: 71 % of the simulated peak rows within 10 MeV/c2 of the edge are tail rows.
3. **The data have a tail of the simulated size.** The tail scale t (t = 1: the SMASH tail; t = 0: no tail) is measured in two
   independent ways:
   - in the u fit itself, with the peak scale a_P fixed: t = 0.76 +- 0.47;
   - in the peak line shape where the bump cannot contribute (M(pi+pi-) 0.42-0.46 GeV/c2): t = 0.66 +- 0.64.
   - Together: **t_data = 0.72 +- 0.38.** "No tail" is disfavoured at 1.9 sigma.
   - The physics expectation is t = 1.10: the SMASH material, plus the QED radiation of the direct decay that the generator
     leaves out (t_C5 = 1.19).
4. **The checker's 0.384 is reproduced exactly (0.38388), but it is not a fair bound; it is an unphysical extreme.**
   - In the family of tail scales it corresponds to t = -0.26, i.e. a negative tail.
   - It mixes three effects:
     - no tail at all: -6.9 %;
     - peak positions taken from a side-band subtraction of the data instead of SMASH: about -3.7 % more;
     - the data tail is left to the step, while h_sim is taken from the bump alone.
   - The same fit treated in the same way for SMASH (the simulated peak included in h_sim) gives **0.414 (-1.4 %)**.
5. **Effect on the bump.**
   - f_eta(u-edge) changes by **+6.8 % per unit of t**, f_eta(edge) (bedge) by **+7.3 % per unit of t**.
   - With sigma(t) = 0.38 this gives **+-2.6 % (u-edge)** and **+-2.8 % (per-slice edge)**.
   - At t_data: -1.9 % and -2.0 %. With the QED tail: +0.6 % and +0.7 %.
   - The two edge methods move together (their difference changes by less than 0.2 %).
6. **The peak side bands are more sensitive than the edge.**
   - a_P and sigma(pp -> eta X) from the peak change by +16 % per unit of t, i.e. **+-6.0 %** for sigma(t).
   - Y_peak changes by +-2.8 %.
   - The double ratio, the ECAL-tag comparison and the X17 limit change by a few per cent at most.

**In short:** the SMASH tail is right within the 40 % precision that the data allow. Keep it as the nominal. A new
systematic term of +-2.6 % (u-edge) / +-2.8 % (per-slice edge) and +-6.0 % (sigma from the peak) is proposed (D40); no headline
number changes.

---------------------------------------------------------------------------------------------------
## Headline table

| quantity | nominal (SMASH tail, t = 1) | t_data = 0.72 | t = 1.10 (+1 sigma) | t = 0.34 (-1 sigma) | no tail t = 0 | t = 2 | QED tail of C5 | proposed term |
|---|---|---|---|---|---|---|---|---|
| f_eta(u-edge) | 0.4200 | 0.4121 (-1.9 %) | 0.4230 (+0.7 %) | 0.4011 (-4.5 %) | 0.3912 (-6.9 %) | 0.4487 (+6.8 %) | 0.4225 (+0.6 %) | **+-2.6 %** |
| f_eta(edge), bedge (fc_rec) | 0.4184 | 0.4102 (-2.0 %) | 0.4216 (+0.8 %) | 0.3995 (-4.5 %) | 0.3904 (-6.7 %) | 0.4513 (+7.9 %) | 0.4211 (+0.7 %) | **+-2.8 %** |
| a_P = sigma_eta(peak) (relative) | 1 | 0.954 | 1.018 | 0.898 | 0.853 | 1.209 | 1.019 | **+-6.0 %** |
| Y_peak (relative) | 1 | 0.980 | 1.008 | 0.955 | 0.934 | 1.093 | 1.009 | +-2.8 % |
| DR_u = a_P / f (relative) | 1 | 0.973 | 1.011 | 0.940 | 0.915 | 1.132 | 1.013 | +-3.6 % (negligible against 38 %) |
| f_eta(tag), proxy (relative) | 1 | 0.956 | 1.018 | 0.901 | 0.856 | 1.202 | 1.016 | listed only |
| X17 limit, normalisation (1) | 1 | -- | 1.005 | 0.969 | 0.955 | 1.059 | 1.010 (1.023 toy rows) | listed only |

- The checker's Gaussian peak term: f = 0.3839 (-8.6 %). The same estimator with the SMASH peak in h_sim: 0.4140 (-1.4 %).
- With a_P re-derived for each t (side bands), the u-edge shifts become -1.8 / +0.7 / -4.0 / -5.8 / +8.2 % for
  t_data / +1 sigma / -1 sigma / 0 / 2. The peak term's size hardly matters, so the change is small.

---------------------------------------------------------------------------------------------------
## 1. Where the tail comes from (simulation; peaktail_sim.py; Fig. peaktail_sources)

**In short:** the generated decays have no tail. The tail is the leptons' energy loss in the material: bremsstrahlung,
simulated by HGeant. Pion tracks and angles add little. Momentum-resolution tails are symmetric and small.

- **Generator (task 1, last item).** The generated four-body mass is a delta function.
  - C5 (eta -> pi+pi-e+e-): weighted quantiles 0.1 / 50 / 99.9 % = 547.998 / 548.000 / 548.002 MeV/c2. No row lies more than
    1 MeV below the median.
  - C4 (eta -> pi+pi-gamma, the photon converts): 546.5 / 547.46 / 547.85 MeV/c2.
  - So:
    - PLUTO's "eta --> e+ + e- + pi+ pi-" writes no radiative photon. The converter on vae25 (read only) has no PHOTOS. Its
      "brems" plugin is nucleon-nucleon bremsstrahlung, and "eta_physics" is disabled.
    - HGeant's three-body decays are phase space.
    - **The conversion pair carries the whole photon four-momentum at its creation point.** C4 rows do NOT have an e+e- that
      carries less than the photon energy; their tail arises after the conversion, like that of C5.
  - The conversion points of C4 are r = 2.4 / 15.4 / 165 mm (5 / 50 / 95 %), z = -134 / -106 / +83 mm.
- **Reconstruction residual** d = M4' - 0.99615 m4_true, 0.30 <= M(pi+pi-) < 0.40:
  - C5 (339 729 rows, N_eff 99 900):
    - Gaussian core sigma 6.07 MeV/c2;
    - below -10 / -15 / -20 / -25 / -30 / -40 MeV/c2: 18.8 / 10.9 / 7.1 / 5.0 / 3.6 / 2.0 %;
    - above +10 / +15 / +20: 3.9 / 1.4 / 0.6 %.
  - C4 (180 rows, N_eff 166):
    - core 4.84 +- 0.43 MeV/c2;
    - below -15 / -25 MeV/c2: 12.9 / 5.1 %;
    - above +15: 1.7 %.
  - The low side is several times heavier than the high side. The tail is a low-side phenomenon.
- **Sources** (linearised residual d = sum_i J_i (p_i - p_i,true) + rest, J_i = dM4/d|p_i| from the reconstructed four-vectors;
  checked: M4 from the stored momenta reproduces m4 to 1e-4 MeV/c2):

  | rows | share | mean d [MeV/c2] | leptons | pions | angles, rest | a lepton with p/p_true < 0.90 |
  |---|---|---|---|---|---|---|
  | C5, d < -15 | 10.9 % | -29.2 | -26.3 | -1.8 | -1.2 | 83 % |
  | C5, d < -25 | 5.0 % | -41.3 | -39.9 | -1.7 | +0.3 | 96 % |
  | C4, d < -15 | 12.9 % | -30.2 | -27.4 | -1.4 | -1.3 | 67 % |
  | C4, d < -25 | 5.1 % | -47.6 | -47.8 | -2.2 | +2.4 | 86 % |
  | C5, core abs(d) < 5 | 50.6 % | -0.3 | -1.2 | -0.0 | +0.9 | 3 % |

  - Leptons with p/p_true < 0.90: C4 6.1 %, C5 7.7 % of all leptons. Below 0.80: 2.6 % and 3.4 %. Pions below 0.95: 2 %.
- **Tail definition used everywhere below:** a row is a tail row when its leptons' momentum deficit shifts M4 by more than
  5 MeV/c2 (J_e+ dp_e+ + J_e- dp_e- < -5 MeV/c2).
  - Tail share in 0.30-0.40: C4 26.3 %, C5 24.0 %.
  - The rows without such a loss are symmetric. C5: below -10/-15/-20 MeV/c2 5.5 / 1.3 / 0.45 %, above +10/+15/+20
    5.0 / 1.7 / 0.76 %. C4: 6.4 / 2.0 / 1.0 % against 7.3 / 2.4 / 1.7 %.
  - This symmetric part is the momentum-resolution part. It is small and the same on both sides.
- **The tail shrinks with M(pi+pi-).** A heavier pion pair leaves the e+e- pair less energy in the eta frame (E*_ee).
  - Fraction of C4 + C5 with d < -15 MeV/c2: 15.5 % (0.30-0.34, E*_ee 180 MeV), 12.9 % (0.34-0.37), 8.2 % (0.37-0.40),
    8.7 % (0.40-0.42), 3.3 % (0.42-0.46, 96 MeV), 0.8 % (0.46-0.50).
  - So the region without bump has a much shorter tail in MeV/c2 than the fit region.
- **The QED radiation that SMASH lacks** (peaktail_rad.py). Final-state radiation off the e+e- pair of the direct decay,
  leading logarithm, collinear:
  - each lepton emits with dP/dz = (c/2)(1 + (1 - z)^2)/z;
  - c = (alpha/pi)[(1 + b^2)/(2b) ln((1 + b)/(1 - b)) - 1], b = the lepton velocity in the pair frame;
  - the photon is emitted at the vertex; the reconstructed lepton momentum is scaled by (1 - z); M4 and M(ee) are recomputed;
    3 copies per row.
  - Results:
    - per-lepton emission probability (z > 0.002) 6.8 % (E1 3.5 %, E2 8.5 %); 13.1 % of the rows emit;
    - 2.5 % of the rows are shifted by more than 15 MeV/c2;
    - **C5 tail share x 1.19 in 0.30-0.40 (E1 x 1.08, E2 x 1.24)**, x 1.23 in 0.415-0.46;
    - with C4 unchanged (a real photon has no such term): t = 1.10 for C4 + C5.
  - Accuracy of the leading-log estimate: about 20-30 %. Pion radiation, interference and multiple emission are left out.

---------------------------------------------------------------------------------------------------
## 2. Where the peak rows lie in u (0.30 <= M(pi+pi-) < 0.40, M(ee) < m_pi0; peaktail_sim.py Sect. 4)

**In short:** 87 % of the simulated peak rows lie inside the u-fit range (-45..+45 MeV/c2), mostly at u = +20..+45. The rows at
u < +20 (17 % of the peak; 670 data rows x a_P) are mostly tail rows (73 %). At the edge itself (abs(u) < 10 MeV/c2) there are
197 expected rows over 8 bins (5 % of the peak, 71 % of them tail rows). For comparison: the step is 1305 rows per bin.

Fractions of each class in u bins (all rows / tail rows / core rows):

| class | u < -45 | -45..0 | 0..10 | 10..20 | 20..45 | > 45 | in fit range | abs(u) < 10 (tail share) |
|---|---|---|---|---|---|---|---|---|
| C4 all | 1.1 / 4.2 / 0 % | 2.6 / 9.8 / 0 % | 1.4 / 2.6 / 1.0 % | 8.8 / 33.3 / 0 % | 70.3 / 44.2 / 79.7 % | 15.8 / 5.9 / 19.4 % | 83.1 % | 3.0 % (76 %) |
| C4 E1 | 1.4 % | 2.4 % | 1.5 % | 9.4 % | 72.6 % | 12.7 % | 85.9 % | 3.3 % (72 %) |
| C4 E2 | 0 | 3.5 % | 1.0 % | 6.3 % | 61.5 % | 27.8 % | 72.2 % | 2.1 % (100 %) |
| C5 all | 0.3 / 1.1 / 0 % | 3.6 / 13.6 / 0.5 % | 4.8 / 12.3 / 2.4 % | 11.2 / 23.4 / 7.3 % | 70.1 / 47.9 / 77.2 % | 10.1 / 1.7 / 12.7 % | 89.7 % | 6.7 % (68 %) |
| C5 E1 | 0.3 % | 2.6 % | 2.5 % | 6.8 % | 74.5 % | 13.3 % | 86.4 % | 3.7 % (96 %) |
| C5 E2 | 0.3 % | 4.1 % | 5.8 % | 13.2 % | 68.1 % | 8.5 % | 91.2 % | 8.1 % (62 %) |
| C4 + C5 | 0.6 / 2.6 / 0 % | 3.2 / 11.8 / 0.3 % | 3.3 / 7.8 / 1.7 % | 10.1 / 28.1 / 4.1 % | 70.2 / 46.2 / 78.3 % | 12.6 / 3.7 / 15.6 % | 86.7 % | 5.0 % (71 %) |

- The peak position u_p = 0.5439 - E per row: 21.1 / 32.6 / 46.0 MeV/c2 (5 / 50 / 95 %). A row reaches u = 0 when it lost u_p.
- **In the u fit, the core and the tail of the SMASH peak term act in opposite directions, and they cancel.**
  - Step rows per unit a_P: core +162.4, tail -158.5, total +3.8 (the edgeu value).
  - Split by class: C4 core +117.3, tail -76.4; C5 core +45.0, tail -82.1.
  - This is why the size of the peak term does not matter (3.8 rows per unit a_P), while its shape does.
  - Removing the tail removes the cancelling part. With the rows kept, the core grows by 1/(1 - 0.25), and the projection goes
    from +3.8 to +216.5 rows per unit a_P: 0.418 x 212.7 = 89 rows of step = -6.9 % of f.

---------------------------------------------------------------------------------------------------
## 3. The tail in the data (peaktail_ufit.py, peaktail_line.py, peaktail_prop.py; Fig. peaktail_constraint)

**In short:**
- The data measure t to about +-0.4: t_data = 0.72 +- 0.38. This is consistent with SMASH (t = 1, 0.7 sigma) and with the QED
  expectation (1.10, 1.0 sigma).
- In the u fit, the tail and the pol2 are degenerate unless the peak scale is fixed. With a free amplitude, t is not determined
  (+-0.9 to +-1.1).

**Model.** Per class, the peak template with the tail scaled by t keeps the number of peak rows (radiation moves rows from
the core into the tail): K(t) = K_core (N - t N_tail)/N_core + t K_tail. K(t) is linear in t, so t can be fitted as one extra
column.

**(i) The u fit itself** (0.30-0.40, the edgeu estimator, peak scale fixed at a_P = 0.418; MC = the WP0 replicas of the peak
templates):

| range, order | t | stat | MC templates | correlation with h | chi2/ndf |
|---|---|---|---|---|---|
| -45..+45, pol2 (edgeu) | 0.47 | 0.50 | 0.33 | +0.58 | 26.2/31 |
| -45..+55, pol2 | 0.68 | 0.40 | 0.25 | +0.37 | 27.9/35 |
| -45..+60, pol2 | 1.01 | 0.37 | 0.22 | +0.28 | 32.5/37 |
| -45..+70, pol3 | 0.69 | 0.44 | 0.26 | +0.61 | 31.6/40 |
| -55..+70, pol3 | 0.95 | 0.41 | 0.24 | +0.56 | 45.0/44 |
| -45..+60, pol2, amplitude free | -0.24 | 1.10 | -- | -- | 28.1/36 (amplitude 0.26 +- 0.07) |
| -45..+70, pol3, amplitude free | 0.00 | 0.91 | -- | -- | 29.4/39 (amplitude 0.30 +- 0.08) |

- The bootstrap (200 copies) confirms the fit error (-45..+60: sd 0.37).
- Mean of the five choices: **t_u = 0.76 +- 0.43 (stat (+) MC) +- 0.20 (choice) = 0.76 +- 0.47.**
- With the Part III composition of the peak (conversion route 0.237, direct route 0.497 instead of 0.418 for both), t_u = 0.57.
- **The degeneracy, in plain words.** Between u = 0 and +20 MeV/c2 the tail is a slow rise. A pol2 can mimic it.
  - With a_P fixed (the side bands fix the core), the data constrain t to +-0.4-0.5.
  - With a_P free, the fit trades core against tail and t is undetermined.
  - A fit without any peak term: h = 1308, f = 0.4210 (+0.2 %). The pol2 absorbs the whole SMASH peak, core and tail.
- **Per M(ee) class** (-45..+60, pol2):
  - E1 (conversions + low-mass direct decays): t = 0.39 +- 0.36.
  - E2 (direct decays only): t = 2.06 +- 0.75.
  - They differ by 2.0 sigma. With the Part III composition: -0.33 +- 0.53 and 2.41 +- 0.75.
  - M(ee) >= 20 MeV (direct decays only, no conversion part): t = 2.04 +- 1.20 (-45..+60), 0.08 +- 2.08 (-45..+45). This class
    has too few rows to constrain t; there f moves by only 3.5 % per unit t.
  - Read as "the conversions have less tail in the data than in SMASH, the direct decays more". The direct-decay excess is in
    the direction of the missing QED radiation (E2: t_C5 1.24 expected). This is not significant, and it cancels in f: with
    t_C4 = 0.4 and t_C5 = 2.0, f moves by 3.77 % x (-0.6) + 3.08 % x (+1.0) = +0.8 %.

**(ii) The line shape where the bump cannot contribute** (M(pi+pi-) 0.42-0.46 GeV/c2; the bump ends at 0.4116; SMASH has one C2
row at 0.415-0.42, so 0.42 is the lower end):
- Model: N - CB (CB per 20 MeV M(pi+pi-) slice) in M4, pol2/pol3 + amplitude x K(t) (SMASH C4 + C5 of the same region,
  M4' with r = 1).
- Signal: 1016 peak rows in P on N = 6967, S/B 0.18. The SMASH tail share in this region is 14.7 %. Only 3.3 % of the rows
  lie below -15 MeV/c2 (softer leptons).

| M4 range, order | t (amplitude free) | error x sqrt(chi2/ndf) | amplitude | t with a_P fixed |
|---|---|---|---|---|
| 490-600, pol2 | 0.86 | 0.44 -> 0.55 | 0.34 | 1.04 +- 0.35 |
| 480-620, pol3 | 0.61 | 0.63 -> 0.69 | 0.24 | 1.07 +- 0.36 |
| 500-590, pol2 | 0.73 | 0.54 -> 0.58 | 0.29 | 1.10 +- 0.36 |
| 470-620, pol3 | 0.43 | 0.60 -> 0.74 | 0.25 | 0.92 +- 0.35 |
| 490-600, pol2, data-like core width (r = 0.893 on the core rows) | 1.15 | 0.44 -> 0.56 | 0.33 | 1.39 +- 0.34 |

- chi2(t = 0) - chi2(t = 1), amplitude free, 490-600 pol2: +3.2 for 0.42-0.46 and +7.2 for 0.415-0.46. "No tail" is disfavoured
  here too.
- The bootstrap (200; 0.415-0.46, a_P fixed) agrees with the fit error (0.37 against 0.35).
- Mean of the four r = 1 choices: **t_c = 0.66 +- 0.55 (stat, scaled) +- 0.33 (choice and core width) = 0.66 +- 0.64.**
- Wider regions: 0.415-0.46 gives t = 1.28 +- 0.41 (amplitude free); 0.415-0.50 gives 1.16 +- 0.59.
- The fitted amplitude here (0.24-0.34) is below a_P = 0.418. That is a property of the peak at large M(pi+pi-), not of the tail.

**(iii) Combination.** The weighted mean of (i) and (ii) (independent data) is **t_data = 0.72 +- 0.38**. With the Part III
composition in (i): 0.60. QED expectation (SMASH material + the QED radiation of C5): 1.10.

---------------------------------------------------------------------------------------------------
## 4. The checker's peak-term check, reproduced and taken apart (peaktail_ufit.py Sect. 1)

**In short:** 0.384 is reproduced to all digits. It is the sum of three things: no tail, peak positions from a side-band
subtraction, and an h_sim that does not see the peak. Treated the same way on both sides, the fit gives 0.414.

- **Reproduction:** h_data = 1217.25 +- 66.29, chi2 24.57/31, Gaussian amplitude 2838 +- 717, h_sim = 3170.95, f = 0.38388
  (-8.61 %). These equal chk_hsim.py.
  - The data noise of the shift (bootstrap, 200): -0.0358 +- 0.0076. The shift is a model difference, not noise.

| variant (edgeu estimator, h_sim = the same linear estimator on C2c + C3) | h_data | f | change |
|---|---|---|---|
| SMASH peak, a_P fixed (edgeu nominal) | 1305.2 | 0.4200 | 0 |
| SMASH peak shape, amplitude free (a_P = 0.329) | 1305.6 | 0.4200 | 0.0 % |
| no peak term (Part II data-only fit) | 1307.9 | 0.4210 | +0.2 % |
| SMASH peak without its tail (t = 0), a_P fixed | 1215.8 | 0.3912 | -6.9 % |
| SMASH peak without its tail, amplitude free (0.289) | 1243.9 | 0.3955 | -5.8 % |
| Gaussian 5.4 MeV at the SMASH peak rows' positions, amplitude free | 1256.1 | 0.3996 | -4.9 % |
| Gaussian 6.09 MeV at the SMASH positions, amplitude free | 1261.8 | 0.4017 | -4.4 % |
| Gaussian 6.09 MeV at the data positions, amplitude free | 1221.9 | 0.3853 | -8.3 % |
| **the checker: Gaussian 5.4 MeV at the data positions (P window minus side bands), amplitude free** | **1217.3** | **0.3839** | **-8.6 %** |
| the checker's Gaussian with the amplitude fixed at a_P x SMASH (3926 rows) | 1183.2 | 0.3808 | -9.4 % |
| **the checker's estimator, h_sim from SMASH bump + SMASH peak x a_P/f (with its tail)** | 1217.3 | **0.4140** | **-1.4 %** |
| the same, SMASH peak without its tail | 1217.3 | 0.3939 | -6.2 % |

- In the checker's estimator, the SMASH peak x a_P projects -96.9 rows onto h (core -25.9, tail -71.0). A data peak with the
  SMASH tail therefore lowers h_data by about 70 rows. The checker's h_sim from the bump alone does not see this.
- About 3.7 % of the -8.6 % come from the peak positions. The data positions come from 37 090 P-window rows minus 0.921 x
  37 466 side-band rows (net 2584). In u they put the peak narrower and later than the SMASH rows. This is either a real
  difference of the peak's (M(pi+pi-), M(ee)) distribution or a side-band artefact: at small M(pi+pi-) and large M(ee) the
  bump edge reaches the low side band. Not resolved here.
- **Fair bound or unphysical extreme? An unphysical extreme.**
  - It equals t = -0.26 in the tail family (f(t) = 0.420 [1 + 0.068 (t - 1)]).
  - It is 2.6 sigma below t_data.
  - A zero tail is not possible physically: the material exists, and HGeant's bremsstrahlung tail is a property of it. The QED
    radiation (absent in SMASH) only adds to it.
  - A fair lower edge is "SMASH without tail" at the SMASH positions: -5.8 % (amplitude free) to -6.9 % (a_P fixed). That
    already lies 1.9 sigma below t_data.

---------------------------------------------------------------------------------------------------
## 5. Propagation into the two edge measurements (peaktail_prop.py, peaktail_bedge.C; Fig. peaktail_effect (a))

**In short:**
- Both edges depend on t with nearly the same slope (6.8 and 7.3 % per unit t). At t_data both move by -2 %.
- Each is computed with the estimator's own code:
  - edgeu_core.fit_step with the nominal settings of edgeu_ana.py;
  - bedge_ana.C (Evaluate, Quant) with the C4/C5 arrays of its nominal simulation configuration replaced in memory.

- **edgeu** (a_P fixed; h_sim = the same linear estimator on C2c + C3):
  - f(t) = 0.3912, 0.3985, 0.4057, 0.4129, 0.4200, 0.4272, 0.4344, 0.4487, 0.4629, 0.4771 for
    t = 0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3. This is linear: d ln f/dt = +0.0683.
  - Per class (the other class at t = 1): C4 3.77 % and C5 3.08 % per unit t.
  - MC noise of the shift (peak-template replicas): +-0.0040 at t = 0 and t = 2 (+-1 % of f). The data part cancels for a
    fixed peak term.
  - The chi2 of the fit rises with t: 27.1 (t = 0), 27.3 (1), 35.5 (2), 51.7 (3) for 32 ndf. The data disfavour t >= 2.
  - With the amplitude free as well: f = 0.3955 (t = 0), 0.4200 (1), 0.4449 (2), 0.4430 (3; a_P then falls to 0.10).
  - QED radiation of C5 (t_C5 = 1.19): f = 0.4225 (+0.6 %). With the toy rows themselves (shape, not only scale): 0.4214
    (+0.3 %).
  - Part III composition: 0.4234 (+0.8 %). With the composition and no tail: 0.3990 (-5.0 %). With the composition and QED:
    0.4249 (+1.2 %).
- **bedge** (nominal fc_rec 0.41841 reproduced; f_all 0.39162 = WP6):
  - t = 0: 0.3904 (-6.7 %); t = 0.34: 0.3995 (-4.5 %); t = 0.72: 0.4102 (-2.0 %); t = 1.10: 0.4216 (+0.8 %);
    t = 2: 0.4513 (+7.9 %). Slope 7.3 % per unit t.
  - QED: +0.7 % (scale) / +0.6 % (toy rows). Part III composition in E_sim: -1.0 %.
  - "No simulated peak in E_sim at all": -4.6 %.
  - The variant "simulation fitted with the data covariance" gives the same shifts (-1.9 / -6.6 / +7.7 %).
  - The tail enters bedge only through E_sim: the data fit has a tail-free Gaussian with a free amplitude. The rest of E_sim
    (peak + backgrounds) is 31.9 rows at t = 1, 223.9 at t = 0 and -160.1 at t = 2. So the SMASH tail projects -192 rows per
    unit t onto E_sim.
  - The bump-only reading of bedge (E_data/(c C2 + C3)) is 0.4276. It does not depend on t, because it ignores the peak in
    E_sim, as the checker's fit does.
- **edgeu against bedge:** the difference of the relative shifts is +0.06 % (t_data), -0.05 / +0.01 % (+-1 sigma), -0.17 %
  (t = 0) and -1.0 % (t = 2). **The comparison "the two agree to 0.5 %" does not change.**

---------------------------------------------------------------------------------------------------
## 6. Other Part III quantities (peaktail_summary.py; Fig. peaktail_effect (b))

**In short:**
- The peak side bands depend on the tail more than the edges do: sigma from the peak moves by +16 % per unit t. At sigma(t),
  that is +-6.0 %, or +-0.026 mb on 0.43 mb, against the measurement error of +-0.16 mb.
- The double ratio, the ECAL tag and X17 change by a few per cent at most.

- **Peak side bands (WP6 estimator Y = [(N - CB)_P - R'(N - CB)_SB]/(1 - R' f), R' = 1.1859).** The tail changes the leakage
  f = (C4+C5)(SB)/(C4+C5)(P) and the fraction of the peak inside P.
  - SMASH: f = 0.1552 (WP6 0.1544; all M(pi+pi-) of the dump), P / all = 0.809.

  | t | f | P / all | Y_peak | a_P = sigma(peak) |
  |---|---|---|---|---|
  | 0 | 0.1069 | 0.887 | x 0.934 | x 0.853 |
  | 0.34 | 0.1225 | 0.860 | x 0.955 | x 0.898 |
  | 0.72 | 0.1409 | 0.831 | x 0.980 | x 0.954 |
  | 1.10 | 0.1608 | 0.801 | x 1.008 | x 1.018 |
  | 2 | 0.2138 | 0.732 | x 1.093 | x 1.209 |
  | QED (C5 x 1.19 / toy rows) | 0.1611 / 0.1612 | 0.801 / 0.796 | x 1.009 | x 1.019 / x 1.026 |

  - Slope of a_P and sigma(peak) around t_data: +0.158 per unit t, x sigma(t) = +-6.0 %. The +-1 sigma(t) points give +6.7 / -5.9 %.
  - **sigma(pp -> eta X) = 0.43 mb (PRELIMINARY, peak side bands):** at t_data it would be 0.41 mb. The proposed term is
    +-0.026 mb. Added in quadrature to the systematic 0.152 mb, it gives 0.154 mb: the quoted 0.16 (meas.) does not change.
  - The WP7 term S6 (resolution r = 0.89/1.10: 9.0 % on sigma) also scales the tail, partly. The overlap is not resolved here.
- **Double ratio DR_u = a_P / f:** x 0.973 at t_data, x 1.011 / 0.940 at +-1 sigma (+-3.6 %), x 1.013 with QED. Negligible
  against its 38 % total.
- **ECAL tag (proxy).** The frozen tag fit is a Gaussian content + pol2 in m5c. Proxy: the same kind of fit (Gaussian 0.5439 /
  6.09 MeV + pol2, 0.49-0.62, equal weights because the data background dominates) on the simulated peak with the tail scaled.
  - T_sim changes by -17 % per unit t, so f(tag) = T_data/T_sim by +17 % per unit t: x 0.956 (t_data), +1.8 / -9.9 % at
    +-1 sigma, +1.6 % with QED.
  - It moves in the same direction as the edges, about 2.5 times more strongly. The tag/edge ratio changes by about
    -2.5 % at t_data.
  - Inside the tag's 9.5 % error ("verification only"). The m5c rows themselves were not available here: this is an estimate.
- **X17, normalisation (1)** (X17 rows = r N5 R_X, N5 = the direct decays fitted in the same cells).
  - The tail acts on the line and on N5 in the same way. Proxy: C5 rows with M(ee) 15.5-24.5 MeV/c2 stand for the line;
    R_X = (line in the peak band 0.5325-0.5575) / (C5 in the X2D cells).
  - The limit (proportional to 1/R_X) changes by -4.5 % (t = 0), -3.1 % (t = 0.34), +0.5 % (t = 1.10), +5.9 % (t = 2),
    and +1.0 % (QED, scale form) / +2.3 % (QED toy rows).
  - Negligible against sX(1) = 24.7 % and the factor-2 band. Without the normalisation to N5, QED alone would lower the line's
    P content by 3.3 %.

---------------------------------------------------------------------------------------------------
## 7. Recommendation (proposal D40, for the analyst; no number changed)

**In short:** keep the SMASH tail as the nominal. Add one systematic term per estimator from the measured tail scale. Replace
the checker's "-8.6 %, not in the systematic error" statement by the measured term and its explanation.

- **D40 (proposed): the low-mass tail of the eta peak.**
  1. **Nominal unchanged** (SMASH tail, t = 1). The data give t = 0.72 +- 0.38 and the physics expectation is 1.10. The
     nominal lies between them, within 1 sigma of both.
  2. **New systematic family "peak tail":**
     - **f_eta(u-edge): +-2.6 %** (slope 6.83 % per unit t x sigma(t) 0.38);
     - **f_eta(edge) per slice: +-2.8 %** (7.27 % x 0.38).
     - Asymmetric alternative around the data value: u-edge +0.7 / -4.5 %, bedge +0.8 / -4.5 %.
     - Effect on the totals: u-edge syst 7.3 -> 7.75 %, total 11.1 -> 11.5 % (0.420 +- 0.048); sigma(bump) 0.438, syst
       0.032 -> 0.034 mb. bedge syst 8.4 -> 8.8 %, total 15.0 -> 15.2 %.
     - Edgeu's U8 (peak template r = 0.89, 0.8 %) can be merged into the new family by taking the maximum.
  3. **sigma(pp -> eta X) and a_P from the peak side bands: +-6.0 %** (slope 15.8 % x 0.38; +-0.026 mb). Y_peak +-2.8 %. The quoted 0.43 +- 0.16
     (meas.) does not change at the quoted precision.
  4. **Listed, not counted:**
     - the double ratio (+-3.6 %);
     - the ECAL-tag comparison: same direction as the edges, ratio -2.5 % at t_data;
     - X17 normalisation (1): -3 / +1 % on the limit;
     - the QED radiation of C5 (+0.6 % u-edge, +0.7 % bedge, +1.9 % sigma(peak));
     - the Part III composition (+0.8 % u-edge, -1.0 % bedge).
  5. **The integration check (f = 0.384, -8.6 %) should be reworded** in Part III (Sect. 1.3, 12.10, 1.13) and in the short
     report as follows. It is a no-tail extreme (equivalent t = -0.26). The data exclude it at 2.6 sigma. The same fit with the
     simulated peak in h_sim gives 0.414 (-1.4 %), inside the new term. The open item "the tail in the data is not checked"
     (report_all.tex line 955) is closed by this measurement.
- **Edge methods compared:** unchanged (the shifts agree to 0.2 % for 0 <= t <= 1.1).
- **Y_peak, DR, ECAL tag, X17 fit:** as in Sect. 6. Only sigma(peak) / a_P gets a term worth quoting.
- **Open, for the lead analysis:**
  - the E1/E2 split of t (0.39 +- 0.36 against 2.06 +- 0.75, 2.0 sigma). If real, it points to fewer conversion tails and more
    direct-decay tails than in SMASH, in the direction of the conversion question (kappa) and of the missing QED radiation.
    It cancels in f (+0.8 %), but it would move the E1/E2 kappa split of the edge. Not evaluated.
  - a simulation with final-state radiation (PHOTOS or a dilepton channel with radiative corrections) would replace the
    leading-log toy.

---------------------------------------------------------------------------------------------------
## 8. Figures (new files, docs/latex/figs_all/; draft captions)

**Fig. `figs_all/peaktail_sources.pdf` (.png).**
- **What it shows:** the tail of the simulated eta peak comes from the leptons' energy loss in the material, not from the
  generator. It reaches the bump edge in u, and it shrinks at large M(pi+pi-).
- (a) Residual d = M4' - 0.99615 m4_true of the simulated peak rows (0.30 <= M(pi+pi-) < 0.40), C5 cyan (2 MeV/c2 bins),
  C4 green (180 rows, 5 MeV/c2 bins). Solid: all rows. Filled: tail rows (the leptons' momentum deficit shifts M4 by more than
  5 MeV/c2; 24 % and 26 %). Dashed: the other rows, symmetric.
- (b) Mean shift of M4 for rows with d < -15 and < -25 MeV/c2, split into lepton momenta (blue), pion momenta (orange) and
  angles/rest (grey). The generated m4 is a delta function: no radiative photon in the generator.
- (c) The peak in u (expected data rows x a_P = 0.418 per 2.5 MeV/c2): core rows dashed, tail rows filled. Grey: the u-fit
  range.
- (d) Fraction of peak rows with d < -15 MeV/c2 per M(pi+pi-) bin: SMASH (black) and SMASH + the QED radiation of C5 that the
  generator lacks (red dashed, leading-log toy). Numbers: the e+e- energy in the eta frame. Shaded: the u-fit region and the
  no-bump region.

**Fig. `figs_all/peaktail_constraint.pdf` (.png).**
- **What it shows:** the data contain a tail of about the simulated size: t_data = 0.72 +- 0.38 (t = 1: SMASH). "No tail" is
  disfavoured, and the checker's Gaussian peak term corresponds to a negative tail.
- (a) u fit (data S0, 0.30-0.40): data minus the smooth part and the step of the nominal fit, with the peak terms a_P K(t) for
  t = 1 (green, nominal), t = 0 (black dashed) and t = 2 (red dotted), and the checker's Gaussian (orange). The smooth part
  was fitted with t = 1, so the points lean towards t = 1; the fitted t is in (c).
- (b) The peak where the bump cannot contribute (0.42 <= M(pi+pi-) < 0.46): N - CB minus the fitted pol2, with the fits t = 1
  (green) and t = 0 (black dashed), amplitude free.
- (c) The tail scale from the u fit (five fit choices, blue; error stat (+) template MC), from the no-bump line shape (four
  choices and a data-like core width, green; error scaled by sqrt(chi2/ndf)) and per M(ee) class (grey). Grey band:
  t_data +- 1 sigma. Lines: SMASH (green), SMASH + QED radiation (red), no tail (black), the checker's 0.384 (orange).

**Fig. `figs_all/peaktail_effect.pdf` (.png).**
- **What it shows:** both bump-edge measurements rise by about 7 % per unit of tail scale, so the measured tail gives +-2.6 %
  (u-edge) and +-2.8 % (per-slice edge). The peak cross section reacts more strongly (+-6 %). The edge comparison, the double
  ratio, the ECAL tag and X17 move by a few per cent at most.
- (a) f_eta(u-edge) (black) and f_eta(edge) per slice (brown) against t; grey band t_data +- 1 sigma; red star: QED radiation
  of C5 added; orange: the checker's Gaussian peak term (dash-dot) and the same fit with the simulated peak in h_sim (dotted).
- (b) Change against t = 1 of f_eta(u-edge), of a_P and sigma_eta from the peak side bands, of the double ratio, of the ECAL
  tag (proxy: Gaussian content + pol2 on the simulated peak) and of the X17 limit (normalisation (1)).

---------------------------------------------------------------------------------------------------
## 9. Caveats

1. **Not independently verified.** All numbers are from this package. The edgeu and bedge nominal values, the checker's
   0.38388 and the WP6 leakage 0.1544 (here 0.1552) were reproduced as checks of the machinery.
2. **Tail definition.** t scales the rows whose leptons' momentum deficit shifts M4 by more than 5 MeV/c2. The shape of the tail
   (the momentum-loss spectrum) is taken from HGeant. A data tail of a different shape (not only size) is not tested.
3. **C4 template statistics.** C4 has only 180 simulated rows in 0.30-0.40 (49 in 0.42-0.46). The template MC noise is in the
   t errors (+-0.2-0.33) and in the shift noise (+-1 % of f at t = 0).
4. **The no-bump fit** has chi2/ndf 1.1-1.6 (pol2 fails at 590-600 MeV/c2); its errors are scaled. Its region has softer
   leptons, so its tail is short in MeV/c2. t is assumed to be the same there as in the fit region. For radiation it is
   (a fractional energy loss); the QED toy gives 1.23 there against 1.19.
5. **QED radiation:** leading-log, collinear toy, 20-30 % accuracy. It acts only on C5 (and on X17-like rows).
6. **ECAL-tag and X17 numbers are proxies** built from the simulated peak rows, not from the m5c rows or the X17 line template.
7. **a_P is held fixed in the u-fit variants.** With a_P re-derived per t, the shifts change by less than 1 % (Sect. Headline).

---------------------------------------------------------------------------------------------------
## 10. Exact commands and files (from /home/przygoda/HADES/FAT/FAT)

    O=/hdd3/przygoda/HADES/FAT_out/all_report/peaktail
    root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/peaktail/aclic", kTRUE)' 'plots/all_report/peaktail_dump.C+O' > $O/peaktail_dump_stdout.txt 2>&1   # 35 s, 12 GB: columns, M(pipi) 0.26-0.62
    python3 -I plots/all_report/peaktail_sim.py      # task 1: generator, residual, sources, u map            -> peaktail_sim.{json,log}, peaktail_simfig.npz
    python3 -I plots/all_report/peaktail_ufit.py     # u fit: checker reproduced and decomposed, f(t), t free   -> peaktail_ufit.{json,log}, peaktail_ufitfig.npz
    python3 -I plots/all_report/peaktail_line.py     # no-bump line shape, u fit per M(ee) class               -> peaktail_line.{json,log}, peaktail_linefig.npz
    python3 -I plots/all_report/peaktail_rad.py      # QED radiation toy of C5; side bands; X17 window         -> peaktail_rad.{json,log}, peaktail_rad_rows.npz
    python3 -I plots/all_report/peaktail_prop.py     # t combination, f(u-edge) cases, bedge input arrays      -> peaktail_prop.{json,log}, bedge_C45_{core,tail,rad}.bin
    root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/peaktail/aclic", kTRUE)' 'plots/all_report/peaktail_bedge.C+O'   # bedge's own code -> peaktail_bedge.txt
    python3 -I plots/all_report/peaktail_summary.py  # side bands, DR, tag proxy, X17; peaktail_numbers.json
    python3 -I plots/all_report/peaktail_figs.py     # docs/latex/figs_all/peaktail_{sources,constraint,effect}.{pdf,png}

- Code (`plots/all_report/`): `peaktail_dump.C`, `peaktail_core.py`, `peaktail_sim.py`, `peaktail_ufit.py`, `peaktail_line.py`,
  `peaktail_rad.py`, `peaktail_prop.py`, `peaktail_bedge.C` (#includes bedge_ana.C unchanged), `peaktail_summary.py`,
  `peaktail_figs.py`.
- Scratch: /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ3/peaktail/look{1..6}.py
  (first looks; every number used is reproduced by the package scripts).

---------------------------------------------------------------------------------------------------
## Verification (independent)

Date 2026-11-15 (verifier, prefix vpt). No peaktail_* code or output was used as input. My own row dump (vpt_dump.C, frozen
headers only), my own edge function, u, CB rule, step template and GLS fit (vpt_core.py), my own tail flags (by putting the
true momenta into the reconstructed four-vectors instead of a linearised Jacobian), my own simulated-event-group bootstrap, an
analytic QED estimate instead of the toy. For the per-slice edge I used my own arrays with the independent bedge
re-implementation vbedge_core.py (written by the bedge verifier). No report, no existing file was changed; this section is the
only edit (backup: integ3/backup/vpeaktail/). No git, rm or mv; one process at a time; CPU package 38-42 C.
**The copy docs/latex/notes_all/peaktail_results.md was not updated** (it now differs from this file by this section).

### V0. In plain words (key results first)

**Verdict: the physics and the arithmetic are CONFIRMED; the measured tail size t_data = 0.72 +- 0.38 is NOT CONFIRMED.**
The data constrain the tail about three times more weakly than stated, so the error assigned from them is too small, and the
statements "no tail disfavoured at 1.9 sigma", "the checker's value excluded at 2.6 sigma" and "the open item is closed" do not
hold.

1. **Where the tail comes from: confirmed.** The generator writes a sharp line (no radiative photon; the converter on vae25 has no
   PHOTOS, its "brems" plugins are nucleon-nucleon bremsstrahlung, eta_physics is disabled). The tail is the leptons' momentum
   loss after creation (HGeant material); the conversion pairs of C4 are created by pair production (GEANT process 6, parent a
   photon, in all 287 rows) and carry the photon's four-momentum to within 1.5 MeV/c2 in M4.
2. **The numbers of the study are reproduced** (simulated tail shares, the u map, the core/tail cancellation, the u fit with a_P
   fixed, the line-shape fit with the raw templates, f(t) for both edge methods, the side-band a_P(t), the checker's 0.38388 to
   all digits, the QED factor).
3. **But the data do not measure t to +-0.38.** Three things were left out:
   - **the peak scale is not a fixed number.** The u-fit value of t used a_P = 0.418 exactly. a_P has an error of 0.080 (19 %)
     and itself depends on t (the side-band leakage: a_P x 0.853 at t = 0, x 1.209 at t = 2). With a_P(t) +- 0.080 the u fit
     gives t = 0.0 +- 0.9 (raw templates) or 0.6 +- 1.0 (smoothed C4), not 0.76 +- 0.47;
   - **the C4 template is a handful of heavy rows**: 180 simulated rows in 0.30-0.40 and 49 rows (each worth about 32 data rows,
     58 % of the peak) in the no-bump region. Their bootstrap noise on the line-shape t is 0.56-0.81, and it is not in the quoted
     0.66 +- 0.64. With C4 smoothed (each row given the C5 residuals of its M(pipi) bin) the line-shape fit improves from
     chi2 60/40 to 44/40 and t moves to 0.61 / -0.68 / 0.30 / -1.61 for the four fit choices;
   - **position and tail cannot be told apart in the line shape**: shifting the SMASH peak by -2 ... +2 MeV/c2 moves the fitted
     t from -1.55 to +2.26 with the same chi2 (44.8 ... 43.8).
   - My reading of the data: **t = 0.3 +- 1.0.** Both the SMASH tail (t = 1) and no tail (t = 0) are allowed by the data.
4. **No tail is excluded by physics, not by the data.** The material is there and HGeant models it; the decay also radiates
   (my analytic leading-log estimate: t = 1.16 for C5, 1.09 for C4 + C5; the study's toy: 1.19 and 1.10).
5. **The checker's 0.384 is a combination, as the study says**: no tail (-5.8 to -6.9 %) plus peak positions from a side-band
   subtraction (-3.9 %). The "fair version 0.414" is not an independent check: it puts the SMASH tail (t = 1) back into the
   calibration of the checker's estimator, so it returns the nominal by construction (with t = 0 there: 0.395).
6. **Effect per unit of t: confirmed.** f_eta(u-edge) +6.8 % (6.4-6.5 % with smoothed C4 or a_P(t)), f_eta(edge) per slice
   +7.2 % (6.2-7.7 % depending on the template), a_P and sigma from the peak side bands +17.7 % (0 <= t <= 2).

**In short:** keep the SMASH tail as the nominal and add the proposed "peak tail" term, but take its size from physics
(sigma_t = 0.5, i.e. half the tail), not from the data: **+-3.4 % (u-edge), +-3.6 % (per-slice edge), +-8.8 % (sigma and a_P
from the peak side bands)**. No headline number changes; the totals grow slightly (Sect. V6).

### V1. Machinery checks (own code against the frozen results)

**In short:** my code reproduces every nominal number it touches to all printed digits.

| quantity | own code | reference | verdict |
|---|---|---|---|
| edgeu nominal: h_data, chi2, h_sim, f | 1305.2164 +- 62.2741, 27.34/32, 3107.3304, 0.42004 | edgeu_ana.json, same | exact |
| edgeu linear projections C2c / C3 / C4 / C5 | 1780.31 / 1327.02 / 40.95 / -37.11 | same | exact |
| bedge nominal: E_data, E_sim, f_all, fc_rec, rest of E_sim | 1141.8584, 2915.6983, 0.39162, 0.41841, 31.9 | bedge / vbedge, same | exact |
| WP6 leakage (C4+C5)(SB)/(C4+C5)(P), frozen windows SBlo 0.5250-0.5336, SBhi 0.5542-0.5628 | 0.1544 | WP6 0.1544 (peaktail 0.1552) | exact |
| the checker: h_data, chi2, amplitude, h_sim, f | 1217.25 +- 66.29, 24.57/31, 2838 +- 717, 3170.95, 0.38388 | chk_hsim.py / peaktail | exact |

### V2. Where the tail comes from (vpt_phys.py, own substitution method)

**In short:** confirmed in every point; the tail is the leptons' momentum loss after creation, modelled by HGeant.

- **Generator.** True four-body mass, weighted quantiles 0.1 / 50 / 99.9 %: C5 (339 729 rows) 547.998 / 548.000 / 548.002 MeV/c2;
  C4 (180 rows) 546.505 / 547.456 / 547.852. Not one C5 row lies 1 MeV below the median.
- **Converter (vae25, read only):** `/lustre/hades/user/hadesdst/sub/feb22/sim/p4500p/convert_smash/smashtrunk/`
  - database_decays_eta.dat: `channel=eta:e+,e-,pi+,pi-; active=1` (the pi0 and gamma channels are inactive);
  - pluto_decay.h: `Disable("eta_physics")`, `Exec("vdm")`, `Exec("elementary")`, `Exec("brems: kaptari; weighting")`,
    `Exec("brems: elastic")`, `Exec("brems: fsi")`, `Exec("dalitz_mod: krivoruchenko")`;
  - `grep -i photos` over the whole converter tree: no match. The "brems" plugins are the nucleon-nucleon bremsstrahlung models.
- **C4 pairs:** both leptons of all 287 C4 rows have GEANT process 6 (pair production) and parent id 1 (photon); creation radius
  2.4 / 15.4 / 164.7 mm (5 / 50 / 95 %); true M(ee) at creation 1.85 MeV/c2 (median; 6.6 at 95 %); m4_true lies 0.54 MeV/c2 below
  548.0 (median) and at most 1.5 MeV/c2 below it for 99.9 % of the rows. So the conversion pair takes the whole photon
  four-momentum; C4 gets its tail after the conversion, like C5. C5 leptons: process 0 (generator) in all rows.
- **Residual and sources** (0.30-0.40, true |p| substituted with the reconstructed angles kept; M4 rebuilt from (p, theta, phi)
  reproduces m4 to 2e-4 MeV/c2):

  | rows | share | mean d | leptons | pions | angles, rest |
  |---|---|---|---|---|---|
  | C5, d < -15 MeV/c2 | 10.9 % | -29.2 | -25.3 | -1.8 | -2.1 |
  | C5, d < -25 | 5.0 % | -41.3 | -38.1 | -1.7 | -1.5 |
  | C4, d < -15 | 12.9 % | -30.2 | -26.1 | -1.4 | -2.6 |
  | C4, d < -25 | 5.1 % | -47.6 | -44.8 | -2.2 | -0.5 |

  The lepton share is 2-3 MeV smaller than with the linearised Jacobian (-39.9 / -47.8); the conclusion is the same.
  - Tail shares (lepton loss shifts M4 by more than 5 MeV/c2): C4 26.3 %, C5 23.9 % (study: 26.3 / 24.0 %). Rows without such a
    loss are symmetric (C5: 5.49 / 1.33 / 0.45 % below -10 / -15 / -20, 4.97 / 1.73 / 0.76 % above).
- **Is it bremsstrahlung?** Lepton p/p_true below 0.90: 6.96 %; above 1.10: 0.71 % (one-sided loss). The integrated
  Bethe-Heitler spectrum gives an equivalent radiator of 2.9 / 2.0 / 1.2 % X0 for losses above 10 / 20 / 30 %: hard losses are
  rarer than in a thin radiator (such tracks are probably lost in the reconstruction). The loss probability grows with momentum
  (P(p/p_true < 0.90) = 1.4 % below 100 MeV/c, 6.2-6.7 % at 100-300 MeV/c, 12.5 % at 300-500 MeV/c), so it is not a purely
  fractional loss. Bremsstrahlung in the material is the natural reading and nothing contradicts it; that it is the only
  mechanism is not proven. For the analysis only the HGeant tail shape matters.
- **u map and cancellation:** reproduced. 86.7 % of the peak in -45..+45; 17.1 % at u < +20 (673 data rows x a_P, 73 % tail);
  5.0 % within 10 MeV/c2 of the edge (197 rows, 70 % tail). Step rows per unit a_P: core C4 +117.3, C5 +44.9; tail C4 -76.4,
  C5 -82.0; total +3.84; without tail +218.4 (study 216.5). Tail versus M(pipi) (d < -15): 15.5 / 12.9 / 8.2 / 8.7 / 3.3 / 0.8 %.
- **QED radiation that the generator lacks** (vpt_qed.py; analytic, no toy): per lepton and row
  P(z > z0) = (c/2)[-2 ln z0 - 2(1 - z0) + (1 - z0^2)/2], z0 = 5 MeV / E*_lepton (eta frame), c from the pair velocity.
  - Per-lepton emission z > 0.002: 7.1 % (study 6.8 %). t_C5 = 1.16 (E1 1.075, E2 1.21) in 0.30-0.40, 1.21 in 0.42-0.46;
    C4 + C5 combined 1.09. Study's toy: 1.19 (1.08 / 1.24), 1.23, 1.10. **Confirmed** within the 20-30 % leading-log accuracy.

### V3. The tail size in the data (vpt_ufit.py, vpt_ucons.py, vpt_line.py, vpt_line2.py, vpt_smooth.py)

**In short:** the study's fits are reproduced, but their errors leave out the peak-scale uncertainty, the C4 template noise
and the position degeneracy. With them, the data give t = 0.3 +- 1.0: they cannot tell the SMASH tail from no tail.

**(i) u fit with a_P fixed at 0.418, raw templates, tail flag as in the study (sub5): reproduced.**

| range, order | t (stat) own | MC, own group bootstrap (200) | study |
|---|---|---|---|
| -45..+45 pol2 | 0.47 +- 0.50 | 0.38 | 0.47 +- 0.50 +- 0.33 |
| -45..+55 pol2 | 0.68 +- 0.40 | 0.29 | 0.68 +- 0.40 +- 0.25 |
| -45..+60 pol2 | 1.01 +- 0.37 | 0.26 | 1.01 +- 0.37 +- 0.22 |
| -45..+70 pol3 | 0.69 +- 0.44 | 0.32 | 0.69 +- 0.44 +- 0.26 |
| -55..+70 pol3 | 0.95 +- 0.41 | 0.29 | 0.95 +- 0.41 +- 0.24 |

- Per M(ee) class (own-class CB): E1 0.25 +- 0.44 / 0.42 +- 0.36, E2 0.82 +- 1.14 / 2.00 +- 0.75 (-45..+45 / -45..+60).
- With the amplitude free, the chi2 profile is one-sided: t = -0.40 ... 0.00 with upper 1-sigma ends 0.50-0.75 (study: "-0.24
  +- 1.10"). The free amplitude (0.26-0.33) is 1-2 sigma below the side-band a_P.

**(ii) The peak scale is a measurement and depends on t.** a_P = 0.418 +- 0.080, and the side-band leakage changes with t
(reproduced: f = 0.1064 / 0.1217 / 0.1400 / 0.1544 / 0.1597 / 0.2126 and a_P x 0.853 / 0.898 / 0.954 / 1 / 1.018 / 1.209 for
t = 0 / 0.34 / 0.72 / 1 / 1.10 / 2; slope +17.7 % per unit t over 0-2). Profile of the u fit in t (Delta chi2 = 1 intervals):

| range, order | a_P fixed | a_P(t) fixed | a_P(t) +- 0.080 | amplitude free |
|---|---|---|---|---|
| raw: -45..+45 pol2 | 0.47 [0.00, 0.95] | -0.05 [-0.90, 0.60] | -0.19 [-1.25, 0.60] | -0.41 [-1.5, 0.50] |
| raw: -45..+55 pol2 | 0.68 [0.30, 1.05] | 0.05 [-0.70, 0.65] | -0.17 [-1.25, 0.55] | -0.32 [-1.5, 0.50] |
| raw: -45..+60 pol2 | 1.01 [0.65, 1.35] | 0.43 [-0.20, 0.95] | 0.03 [-1.05, 0.75] | -0.24 [-1.5, 0.65] |
| raw: -45..+70 pol3 | 0.69 [0.25, 1.10] | 0.24 [-0.45, 0.80] | 0.12 [-0.80, 0.80] | 0.00 [-1.15, 0.75] |
| raw: -55..+70 pol3 | 0.95 [0.55, 1.35] | 0.39 [-0.25, 0.95] | 0.11 [-0.90, 0.80] | -0.17 [-1.5, 0.70] |
| smoothed C4: -45..+45 pol2 | 0.40 [-0.40, 1.20] | 0.00 [-1.5, 1.00] | -0.01 [-1.5, 1.00] | -0.02 [-1.5, 1.10] |
| smoothed C4: -45..+60 pol2 | 1.46 [1.00, 1.90] | 0.93 [0.00, 1.60] | 0.90 [-0.20, 1.70] | 0.89 [-0.30, 1.70] |
| smoothed C4: -45..+70 pol3 | 1.13 [0.50, 1.70] | 0.89 [0.00, 1.60] | 0.93 [0.00, 1.60] | 0.96 [0.00, 1.70] |

- Smoothed C4: every C4 row is replaced by 400 entries with the residuals of C5 rows (M(ee) < 20 MeV/c2) of its 10 MeV M(pipi)
  bin; its template MC error is 0.13-0.17 instead of 0.26-0.38. It changes the u-fit t by up to +0.45 (a_P fixed) and +0.9
  (a_P constrained): the C4 rows matter.
- Lower bound 1.5 = the edge of my scan (open).

**(iii) The line shape where no bump exists (0.42-0.46 GeV/c2).**
- Raw templates: reproduced (pol2 490-600: 0.87 [0.45, 1.30]; pol3 480-620: 0.62; pol2 500-590: 0.73; pol3 470-620: 0.44).
- **The template noise is not in the study's error.** Bootstrap of the simulated event groups (60-100 copies): sd 0.56-0.81
  (individual copies from -0.05 to 2.74). The cause: C4 has 49 rows of about 32 data rows each, 58 % of the peak there.
- **With C4 smoothed:** chi2 60.3/40 -> 44.1/40 (the poor chi2 was the C4 granularity). t = 0.61 [-0.20, 1.30] (pol2 490-600),
  -0.68 [-1.5, 0.50] (pol3 480-620), 0.30 [-0.90, 1.30] (pol2 500-590), -1.61 [.., -0.40] (pol3 470-620). 0.415-0.46: 1.07,
  0.64, 1.13, -0.14.
- **Position and tail are degenerate.** Shifting the SMASH peak by dm = -2 / -1 / 0 / +1 / +2 MeV/c2 gives t = -1.55 / -0.46 /
  0.61 / 1.52 / 2.26 with chi2 44.8 / 44.0 / 44.1 / 44.0 / 43.8 (pol2 490-600). The global mass transfer is known to +-0.5
  MeV/c2 (s +- 0.0009), the M(pipi) dependence not at all: one MeV of position is one unit of t.
- Resolution variants (raw, r = 0.893): on all rows t = 2.15 / 2.00; on the non-loss part only 1.16 / 0.89.
- **The shape family matters.** Scaling only the rows with losses above 10 MeV/c2 (sub10) gives t = -1.2 ... -5 in the line
  shape (unphysical) and 0.3-1.3 in the u fit: the line-shape t follows the background shape and the position, not the tail.

**(iv) My combination.** The u fit with the side-band constraint is the only handle that does not hinge on a position assumption:
t = 0.0 +- 0.9 (raw) to 0.6 +- 1.0 (smoothed). **t_data = 0.3 +- 1.0.** The line shape adds no robust information.
- t = 1 (SMASH) is 0.7 sigma away; t = 0 is 0.3 sigma away; the physics expectation 1.09-1.10 (confirmed) is 0.8 sigma away.
- The study's statements to withdraw: "t_data = 0.72 +- 0.38", "no tail is disfavoured at 1.9 sigma", "the data have a tail of
  the simulated size", "the SMASH tail is right within the 40 % precision that the data allow", "the checker's value is
  excluded by the data at 2.6 sigma", "the open item is closed by this measurement", and in the caveats "the template noise is
  included in the t errors" (true for the u fit only, not for the line shape).

### V4. Propagation into the estimators (vpt_ufit.py, vpt_smooth.py, vpt_bedge.py, vpt_sb.py)

**In short:** confirmed. Per unit of t: u-edge +6.8 %, per-slice edge +7.2 %, peak a_P +17.7 %; the two edges move together.

- **f_eta(u-edge), a_P fixed, raw sub5:** f = 0.3913 / 0.4011 / 0.4057 / 0.4120 / 0.4200 / 0.4229 / 0.4344 / 0.4486 for
  t = 0 / 0.34 / 0.5 / 0.72 / 1 / 1.10 / 1.5 / 2; d ln f / dt = +0.0683 (study 0.0683); chi2 27.1 / 27.3 / 35.5 at t = 0 / 1 / 2.
  C4 alone at t = 0: -3.77 %, C5 alone -3.07 % (study 3.77 / 3.08). Template noise of f(0)/f(1) - 1: 0.94 % (study 1 %).
  - Slope with other choices: smoothed C4 6.41 %, a_P(t) instead of fixed 6.50 %; flag p95 6.76 % (smoothed 6.05 %); sub10
    3.6-4.3 %; res15 3.8-4.4 % (these scale a smaller part of the peak, so their t is less constrained still).
  - Smoothed C4 changes the nominal itself by -0.5 % (0.41782): C4 template noise, inside the MC error (7.0 %).
- **f_eta(edge) per slice (bedge):** global tail normalisation -6.62 / -4.48 / -1.95 / +0.72 / +7.73 % for t = 0 / 0.34 / 0.72 /
  1.10 / 2 (study -6.7 / -4.5 / -2.0 / +0.8 / +7.9); slope 7.18 % (study 7.27 %); normalised per 5 MeV M(pipi) bin 7.66 %;
  smoothed C4 6.15 %; p95 7.0-7.1 %. Smoothed C4 moves the bedge nominal by -0.9 % (0.41452), inside its MC error (11 %).
- **The two edges move together:** 6.8 against 7.2 % per unit t. "They agree to 0.5 %" stays.
- **Peak side bands:** a_P(t) above; slope 17.7 % over 0-2 (16.8 % locally at t = 0.72-1.10; the study's 15.8 % is the local
  slope at t_data). Y_peak x 0.935 / 0.965 / 1 / 1.042 / 1.092 for t = 0 / 0.5 / 1 / 1.5 / 2.

### V5. The checker's estimator (vpt_ufit.py Sect. 3)

**In short:** 0.384 = no tail (about -6 %) + data positions (about -4 %); the "fair 0.414" contains the SMASH tail by construction.

| variant | f | change |
|---|---|---|
| nominal (SMASH peak, a_P fixed) | 0.4200 | 0 |
| SMASH shape, amplitude free | 0.4200 | 0.0 % |
| no peak term | 0.4210 | +0.2 % |
| SMASH without tail, a_P fixed / free | 0.3913 / 0.3955 | -6.9 / -5.8 % |
| Gaussian 5.4 / 6.09 MeV at the SMASH positions, free | 0.3996 / 0.4017 | -4.9 / -4.4 % |
| Gaussian 6.09 MeV at the data positions, free | 0.3853 | -8.3 % |
| **the checker (Gaussian 5.4 MeV, data positions, free)** | **0.38388** | **-8.6 %** |
| the checker's Gaussian, amplitude fixed at a_P x SMASH | 0.3808 | -9.4 % |
| checker's estimator calibrated with a_P K(t): t = 1 / 0.72 / 0 / 2 | 0.4144 / 0.4089 / 0.3946 / 0.4342 | -1.3 / -2.7 / -6.0 / +3.4 % |

- The checker's estimator projects -96.9 rows of the SMASH peak x a_P onto h (core -26.0, tail -70.9), as the study says.
- The data positions (P window minus side bands) lie 1.0 MeV/c2 later in u than the SMASH rows (mean 34.34 against 33.30
  MeV/c2); this part is -3.9 % of f. Unexplained (a real difference of the peak's (M(pipi), M(ee)) or a side-band effect).
- Data noise of the checker's shift: sd 0.0106 (60 own bootstrap copies; study 0.0076 with 200).
- **Reading.** The "equivalent t = -0.26" is a formal mapping (it puts the position effect onto t). With t_data = 0.3 +- 1.0 the
  no-tail part of the check (t = 0) is not excluded by the data; it is excluded by physics (the material exists, the decay
  radiates). The calibrated value 0.414 is the nominal model expressed in the checker's estimator, not an independent bound.

### V6. Recommendation (for the analyst; no headline number changes)

**In short:** D40 as proposed, but sigma_t = 0.5 from physics instead of 0.38 from the data; reword the checker statements.

- **Nominal t = 1** (SMASH tail). Physics expectation 1.09-1.10 (confirmed). A lower value is possible if the conversion deficit
  kappa = 0.48 +- 0.11 is partly a material effect (pair production and bremsstrahlung both scale with x/X0); this link is a
  hypothesis, not tested here.
- **Size of the term: sigma_t = 0.5** (half the simulated tail; covers 0.5-1.5, i.e. a material deficit like kappa and the
  missing QED radiation with margin; the data, t = 0.3 +- 1.0, allow it and cannot narrow it). Terms:

  | estimator | slope per unit t | proposed term (sigma_t = 0.5) | effect on the quoted errors | option sigma_t = 1.0 (data alone) |
  |---|---|---|---|---|
  | f_eta(u-edge) | +6.8 % | **+-3.4 %** (U8 merged into it by the family maximum) | syst 7.3 -> 8.0 %, total 11.1 -> 11.6 % (0.420 +- 0.049); sigma(bump) syst 0.032 -> 0.035 mb | +-6.8 %; syst 9.9 %, total 13.0 % |
  | f_eta(edge) per slice | +7.2 % | **+-3.6 %** | syst 8.4 -> 9.1 %, total 15.0 -> 15.5 % | +-7.2 %; syst 11.1 %, total 16.7 % |
  | a_P, sigma(pp -> eta X) from the peak side bands | +17.7 % | **+-8.8 %** (+-0.038 mb) | syst 0.152 -> 0.157 mb; total (meas.) stays 0.16 mb | +-17.7 % (0.076 mb): 0.16 -> 0.18 mb |
  | Y_peak | +7.8 % | +-3.8 % | listed | |
  | DR_u = a_P / f | +11 % | +-5.5 % | negligible against 38 % | |

- **Listed, not counted:** no tail (t = 0): u-edge -6.9 %, per-slice edge -6.6 to -7.0 %, a_P -14.7 % (physically excluded,
  data-allowed); the checker's 0.384 (-8.6 %, = no tail + data positions); the QED radiation of C5 (+0.6 / +0.7 / +1.9 %); the
  data-position effect (-3.9 %, unexplained); ECAL-tag and X17 proxies (not re-checked here).
- **Open (for the lead analysis):** a larger simulated eta -> pi+pi-gamma sample (C4: 180 rows in the fit region, 49 in the
  no-bump region) and a simulation with final-state radiation would make the tail measurable; the peak position versus M(pipi)
  in the data would break the line-shape degeneracy; the E1/E2 difference of t (E1 0.42 +- 0.36, E2 2.00 +- 0.75 with a_P fixed)
  was not re-evaluated with the a_P constraint.

### V7. Exact commands and files (from /home/przygoda/HADES/FAT/FAT)

    root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/vpt/aclic", kTRUE)' 'plots/all_report/vpt_dump.C+O'   # 35 s
    python3 -I plots/all_report/vpt_nom.py      # machinery: edgeu nominal, spectra -> vpt_nom.npz
    python3 -I plots/all_report/vpt_phys.py     # generator, substitution sources, momentum loss, C4 pairs -> vpt_phys.log
    python3 -I plots/all_report/vpt_ufit.py     # u fit t (a_P fixed), f(t), the checker -> vpt_ufit.{log,json}
    python3 -I plots/all_report/vpt_line.py     # no-bump line shape, raw templates, template bootstrap -> vpt_line.{log,json}
    python3 -I plots/all_report/vpt_line2.py    # line shape with position shift and width (raw) -> vpt_line2.{log,json}
    python3 -I plots/all_report/vpt_sb.py       # side-band leakage and a_P(t) -> vpt_sb.{log,json}
    python3 -I plots/all_report/vpt_smooth.py   # smoothed C4: u fit with a_P(t) +- 0.08, line shape with dm, r -> vpt_smooth.{log,json}
    python3 -I plots/all_report/vpt_ucons.py    # raw templates, u fit with a_P(t) +- 0.08 -> vpt_ucons.log
    python3 -I plots/all_report/vpt_qed.py      # analytic leading-log FSR estimate -> vpt_qed.log
    python3 -I plots/all_report/vpt_bedge.py    # per-slice edge versus t (own arrays + vbedge_core) -> vpt_bedge.{log,json}

- Code: `plots/all_report/vpt_{dump.C,core.py,tmpl.py,nom.py,phys.py,ufit.py,line.py,line2.py,sb.py,smooth.py,ucons.py,qed.py,bedge.py}`.
- Outputs: `/hdd3/przygoda/HADES/FAT_out/all_report/vpt/` (column files of the dump, the logs and json files above).
- Scratch: `/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ3/vpt/`
  (look_line.py, look_umap.py, win.C: first looks; every number used above is printed by the package scripts, except the
  line-shape bootstrap spread 0.56 (look_line.py) and the u-map fractions (look_umap.py)).
