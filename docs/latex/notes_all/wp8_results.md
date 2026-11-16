# WP8 physics quantities (prefix wp8) -- results note

Date 2026-11-14. Branch pp45_pippimepem (worktree /home/przygoda/HADES/FAT/FAT). No git command that changes state was run.
ScanEngine.h, scan_*.C, hists_v1.root, the sealed file and every file of WP0-WP7/X17 were not touched; no package was re-run
(wp7_syst.C is #included by wp8_syst.C to reuse its fit chain, unchanged and not executed). FAT_simwt was only read (one table).
Nothing was deleted. Remote: two read-only awk passes over converted SMASH event files on vae25 (no remote writes).

**Machine-readable output:** `/hdd3/przygoda/HADES/FAT_out/all_report/wp8/results.txt` (P1, eps2, S12, P2 lines).

## Headline (S0, isBest >= 0, L = 6079.2 nb^-1; yields raw rows x 64/63)

| quantity | value | stat | syst | MC | norm. | status |
|---|---|---|---|---|---|---|
| sigma(pp -> eta X) [mb] (peak side bands) | **0.437** | 0.037 | 0.070 | 0.026 | +0.035/-0.037 | MAIN; live time +0.032/-0 apart [v8: 0.429 +- 0.036 +- 0.079 +- 0.025 +0.035/-0.037 after the verification (PDG 2026 BRs, peak-composition systematic); about 3 sigma below the published inclusive HADES value at 3.5 GeV; see "Verification"] |
| sigma / sigma_SMASH (1.0665 mb) | **0.410** | | total 0.085 | | | |
| BR_G (F-gold) [1e-4] | 2.65 | 0.37 | 3.37 | 0.66 | 0.09 | CROSS-CHECK; PDG 2.68 +- 0.11 (pull -0.0) |
| BR_A (F-main) [1e-4] | 3.52 | 0.31 | 3.24 | 0.61 | 0.11 | CROSS-CHECK; pull +0.26 |
| BR_C [1e-4] | 3.6 (E1+E2 fit only) | 0.9 | -- | 3.6 (mc68) | | WITHDRAWN (E2: R_34 undetermined) |
| R_34 | 0.92 (E1+E2); -4.2 +- 5.6 (E2) | 0.19 | | 0.29 | | consistent with 1 |
| DR (side bands / B4) | 1.40 | 0.15 | 0.94 | 0.53 | 0.06 | MAIN, bump EXTRAPOLATED |
| DR (F-main) | 1.17 | 0.06 | 0.66 | 0.13 | 0.05 | CROSS-CHECK |
| kappa / kappa_B (F-main) | 0.83 / 1.95 | 0.14 / 0.17 | 0.44 / -- | 0.27 / 0.91 | | constraint dominated |
| sigma_vis peak [nb] | 0.837 | 0.071 | 0.096 | 0.040 | 0.058 | MAIN |
| sigma_vis bump [nb] | 2.70 | 0.19 | 4.24 | 0.68 | 0.19 | MAIN, EXTRAPOLATED |

**Two findings of WP8 that the report must carry:**
1. **Class C2 (bump pi0 Dalitz, from HGeant) has the right pair-mass spectrum but the wrong gamma* decay angle.** HGeant
   (gdalet.F) decays gamma* -> e+e- isotropically; QED gives 1 + cos^2 theta*. PLUTO (lepton set) has the QED angle
   (measured here: <cos^2> = 0.403/0.397, QED 0.400). This explains the "unexplained" lep/std closure 0.891 of the simulation
   study. Consequence: the simulated C2 efficiency is 11-12 % too high (toy 0.880, closure 0.891). Refit effects:
   a_eta +9/+12 %, a_B +9/+12 %, BR_G -7/-9 %, BR_A -9/-12 %, DR -8/-10 %; sigma_eta (peak) unchanged (-0.1 %).
   Not applied (frozen procedure), entered as S12 (sign known) and as a labelled variant.
2. **The (pT4, y4) reweighting as specified (3 x 3) is noise-dominated** (+6.2 +- 16.1 %, counted 0), but the
   one-dimensional pT reweighting is significant: +21.0 +- 9.1 % (y only: -8.0 +- 4.4 %) on every absolute eta scale. Data
   show a softer pT4 spectrum than SMASH (2.3 sigma). Ratios (DR, BR) move by < 2 %. This is a production-model
   extrapolation uncertainty of sigma(pp -> eta X) of about +19/-7 % that is NOT in the total (analyst's decision, flagged).
   [v8: no analyst decision on this item is recorded in the spec or in the notes; leaving it out of the total is WP8's
   choice. It is listed next to the total as an asymmetric production-model item in tables_all/v8_tableP2_corr.tex.]

---------------------------------------------------------------------------------------------------
## 1. What was done, step by step (plain words)

1. **Simulated efficiencies (wp8_sim.C, one row pass, 30 s, ~9 GB).**
   - Simulation rows: std + lep + eta (e+e- trees), mode A, owned rows, S0 and pre.
   - For each class C2-C6 the visible cross section sum w / L_sim is filled with event-group errors. Variants:
     - physics weights V3 (code 108) and off;
     - windows all M4 / B / P / X17c-M4;
     - M(ee) classes all / E1 / E2 / >= 20 / 15.5-24.5 MeV (reconstructed) and true M(ee) >= 20 / >= 10.
   - Mass transfer as everywhere: M4' = 0.99615 M4, M(pipi)' = 0.9952 M(pipi).
   - eps = sigma_vis / sigma_prod with sigma_prod = 1.0665 mb x PDG BR. C3 is counted per eta -> pi+pi-pi0 (the conversion is
     inside eps); C6 is per unit BR.
   - Check: pre, all M4, with the sim study's 1.069 mb, reproduces Tab. s5eps for all five classes to 0.1 %
     (C5 0.4945 % vs 0.4944 %).
2. **Pion-vertex efficiency, data vs simulation**, recomputed with the report definitions.
   - Peak: Y_peak(S0)/Y_peak(pre) from the WP7 side-band values, with a nested-sample error, against (C4 + C5)(S0)/(pre).
   - Tag: T(S0)/T(pre) of the ECAL-tagged eta.
3. **eps2 QED check, part 1: generator.**
   - From the HGeant source (hgeant2-7.4 uhbook.F/gdalet.F/gdecay.F): the pi0 Dalitz pair mass is sampled from a 0.1 MeV
     histogram of Kroll-Wada x a rho-pole form factor.
   - Its fraction of decays with M(ee) >= 20 MeV is computed exactly (wp8_qed.py). It is compared with Kroll-Wada x the PDG
     pi0 form-factor slope (a = 0.0332).
   - PLUTO (lepton set): 124 378 pi0 Dalitz decays were read from 3 converted event files on vae25 (remote awk, read-only):
     the pair mass and the e+ decay angle theta* (to the gamma* direction, gamma* frame).
4. **eps2 part 2: decay angle.**
   - The source shows that HGeant gdalet.F samples cos theta* uniformly (isotropic).
   - Toy (wp8_toy_angle.py), using the same events with weights:
     - eta (pT, y) from the SMASH generator;
     - eta -> pi+pi-pi0 with the KLOE-2 Dalitz plot, pi0 -> gamma e+e- with Kroll-Wada;
     - acceptance: 18-85 deg, pions p > 0.1, leptons p > 0.04/0.06/0.10 GeV/c, optional OA > 2.
     It gives eps(QED angle)/eps(isotropic) per M(ee) class.
   - Second estimate: the full-simulation lep/std closure per M(ee) bin (simulation study, closure_ratio_bins.txt), divided by
     the PLUTO/HGeant pair-mass ratio.
   - Both factor sets were applied to the C2 templates (E1, E2, F5g; wp8_syst.C). The complete WP7 chain was refitted
     (F-main, F-gold, stage A, B4, side bands); the nominal reproduces WP7 to all digits.
5. **eps2 part 3: ECAL-tagged eta per M(ee) class** (< 8, 8-20, 20-30, 30-100 MeV; frozen m5c fit of WP6).
   - Data, the emulated simulated tag (same fit) and the truth C2 + C3 of the tagged simulated rows, in S0.
   - The data scan's numbers (0.077 / 0.144) are reproduced on its own selection (0.076 / 0.135).
6. **(pT4, y4) reweighting (wp8_pty.py).**
   - Produced eta (pT, y_cm) from 40 converted SMASH standard-set files (vae25, 94 517 etas; 0.023629 per collision =
     1.069 +- 0.004 mb, consistent with 1.0665).
   - 3 x 3 cells at the tertiles of the simulated peak signal.
   - Data peak yield per cell: WP6 side-band estimator with the global R' = 1.1859 and leakage f = 0.1544. The CB of a cell
     is the window CB x the cell's like-sign share. The sum is 5089 (WP6: 5088).
   - r_c = data/SMASH per cell. The reweighted production gives sigma'/sigma_SMASH = sum_c g_c r_c (g_c = produced
     fraction), against the accepted-weighted a_P.
   - Noise: 2000 toys.
   - The same weights are applied to the bump with its reconstructed (pT4, y4) as a proxy (photon missing), and refitted in
     wp8_syst.C. The pT-only and y-only binnings are information.
7. **Combination (wp8_physics.py).**
   - S12 = max(pion vertex, eps2 angle, pT/y 3x3 after Barlow subtraction): the spec's family rule. The quadrature version
     is given too.
   - New totals replace S12 in the WP7 totals. Table P2, Fig. P1, TeX snippets.
8. **Literature search** (web). Found: no published inclusive value near p_lab 5.4 GeV/c; exclusive values and models are
   listed in Sect. 2.6.

---------------------------------------------------------------------------------------------------
## 2. Numbers

### 2.1 Table P1 (S0, V3; eps in per cent; MC error)
| class | sigma_prod [nb] | all M4 | window | E1 | E2 | >= 20 | window, no weights |
|---|---|---|---|---|---|---|---|
| C2 (B) | 2870 | 0.2253(65) | 0.1791(59) | 0.0641 | 0.1150 | 0.0673 | 0.1934 |
| C3 (B) | 244 442 | 0.001988 | 0.001583 | 0.001244 | 0.000339 | 7e-6 | 0.001663 |
| C4 (P) | 45 006 | 0.00259 | 0.00212 | 0.00178 | 0.00034 | 0 | 0.00200 |
| C5 (P) | 285.8 | 0.4573(12) | 0.3660(10) | 0.1230 | 0.2430 | 0.1515 | 0.3525 |
| C6 (X17c, per unit BR) | 1.0665e6 | 0.6172 | 0.3283 | -- | -- | -- | 0.3283 |

- Simulated S0/pre per class: 0.914 (C2), 0.936, 0.948, 0.933, 0.930 (C6).
- Accepted C2 in B with M(ee) >= 20 MeV: 0.376 (reconstructed), 0.377 (true); generator fraction 0.269.
- **Pion-vertex efficiency S0/pre:**

  | sample | data | simulation | data/sim |
  |---|---|---|---|
  | eta peak | 0.941 +- 0.052 | 0.940 | 1.002 +- 0.056 |
  | ECAL-tagged eta | 0.933 +- 0.013 | 0.933 (emulated), 0.947 (truth) | 1.000 +- 0.014 (0.985) |

  The WP7 S12 value of +-2.5 % is therefore conservative. It is kept.

### 2.2 eps2 QED check (tables_all/wp8_eps2.txt)
- **Generator, fraction with M(ee) >= 20 MeV:**

  | Kroll-Wada | KW x PDG form factor | HGeant (C2) | PLUTO (lepton set) |
  |---|---|---|---|
  | 0.2679 | 0.2693 | 0.2688 | 0.2746 +- 0.0013 |

  - HGeant/QED = 0.998: **PASS**.
  - PLUTO is 2 % harder. PLUTO feeds G_pi0D, not C2.
- **Decay angle:**
  - HGeant is isotropic (source).
  - PLUTO <cos^2 theta*> = 0.4025 +- 0.0012 (< 10 MeV) and 0.3970 +- 0.0017 (>= 20): QED.
  - **FAIL for C2.**
- **eps(QED)/eps(iso):**

  | estimate | all | E1 | E2 | >= 20 |
  |---|---|---|---|---|
  | toy, p_e > 0.06 | 0.880 | 0.889 | 0.868 | 0.886 |
  | toy, p_e 0.04 / 0.10 | 0.900 / 0.858 | | | 0.911 / 0.856 |
  | closure (q-corrected) | 0.891 | 0.869 | 0.896 | 0.909 |
- **Tagged shares** (S0; data / sim fit / sim truth):

  | M(ee) [MeV] | data | sim fit | sim truth |
  |---|---|---|---|
  | < 8 | 0.393(15) | 0.423(29) | 0.459(30) |
  | 8-20 | 0.405(15) | 0.387(30) | 0.328(28) |
  | 20-30 | 0.071(8) | 0.074(17) | 0.049(12) |
  | 30-100 | 0.131(11) | 0.116(17) | 0.163(22) |
  | **>= 20** | **0.202(14)** | **0.191(24)** | **0.213(25)** |

  - Data/sim at >= 20 MeV = 1.06 +- 0.15 (fit), 0.95 +- 0.13 (truth): **PASS**, limited by the simulation statistics.
  - The test is blind to the angle effect, which scales C2 almost uniformly in M(ee): the toy share ratio is 1.007.
- **Refit shifts** (toy / closure factors): see the headline. In addition: Y_bump +3.4/+1.5 %, kappa -3.4/-2.6 %,
  Y4 +8.4/+6.5 %, Y5 -1.3/-1.0 %.
- **Supporting indications** (not proof): with the angle correction, DR_F moves 1.17 -> 1.09/1.10 (towards 1), and f_eta(edge)
  moves 0.39 -> about 0.42, close to a_P = 0.418.

### 2.3 S12 completed (tables_all/wp8_S12.txt, per cent)

| quantity | pivtx | eps2 (angle) | pT/y 3x3 | S12 new | syst old -> new |
|---|---|---|---|---|---|
| sigma_eta, a_P | 2.5 | 0.1 | 0 | 2.5 | 16.01 -> 16.01 |
| Y_peak | 0 | 0.1 | 0 | 0.1 | 11.43 -> 11.44 |
| a_eta (sigma_eta_F) | 2.5 | 12.2 | 0 | 12.2 | 53.7 -> 55.0 |
| BR_G | 2.0 | 8.6 | 0 | 8.6 | 126.6 -> 126.9 |
| BR_A, rho_d | 2.0 | 12.1 | 0 | 12.1 | 91.1 -> 91.9 |
| DR | 2.0 | 10.2 | 0 | 10.2 | 66.3 -> 67.1 |
| DR_F | 2.0 | 7.0 | 0 | 7.0 | 56.4 -> 56.8 |
| a_B | 2.5 | 12.4 | 0 | 12.4 | 167.8 -> 168.2 |
| Y_bump | 0 | 3.4 | 0 | 3.4 | 157.2 -> 157.3 |
| Y4 / Y5 / kappa | 0 / 0 / 2 | 8.4 / 1.3 / 3.4 | 0 | | 57.1 -> 57.7 / 88.7 / 52.8 -> 52.9 |

- The quadrature version changes the totals by < 0.1 point.
- **Outside the total:** the 1D pT/y production extrapolation, +18.9 / -6.7 % after the Barlow subtraction, on the absolute
  eta scales.

### 2.4 Table P2 (tables_all/wp8_tableP2.{txt,tex})

- **sigma(pp -> eta X) = 0.437 +- 0.037 (stat) +- 0.070 (syst) +- 0.026 (MC) +0.035/-0.037 (norm.) mb.**
  - Total 20.7 %; live time +0.032/-0. This is 0.410 +- 0.085 of SMASH.
- **Cross-checks of sigma:**

  | estimator | sigma [mb] | errors |
  |---|---|---|
  | F-main | 0.400 | +- 0.018 +- 0.220 +- 0.041; with the QED angle 0.437-0.449 |
  | f_eta(edge) | 0.410 | f = 0.392 +- 0.024 (stat, x 1.2 report_exp calibration) +- 0.039 +- 0.021 |
  | f_eta(tag), emulated, UNVERIFIED | 0.475 | |
  | P_main | 0.422 | S = 4388 +- 307 +- 420 (model = report_exp edge-safe rms); 0.403 +- 0.028 +- 0.039 of the simulated fit |
- **(pT4, y4) variant:**
  - 3 x 3: +6.2 +- 16.1 %; DR -1.8 +- 2.5 %; a_eta (refit) +8.1 %.
  - pT-only: +21.0 +- 9.1 %; y-only: -8.0 +- 4.4 %.
  - Per-cell data/SMASH is 0.23-0.77, chi2 15.0/8 (p = 0.06). It is highest at low pT and backward y.
- **sigma x BR** (from sigma x the PDG BR unless stated):

  | channel | sigma x BR |
  |---|---|
  | eta -> pi+pi-e+e- | 117 nb (+-21 %) |
  | eta -> pi+pi-e+e-, measured with F-main a_eta rho_d (cross-check) | 141 nb |
  | eta -> pi+pi-gamma | 18.5 ub |
  | eta -> pi+pi-pi0, from B4 (extrapolated, 168 %) | 75 ub |
- **Visible cross sections** (Y/L, C_sec(1 - f_eta) not applied):
  - peak: 0.837 nb in total; E1 0.391 +- 0.048 +- 0.027 nb, E2 0.457 +- 0.054 +- 0.029 nb;
  - bump (E2, not extrapolated): 1.287 +- 0.075 +- 0.316 nb;
  - F-main per class (cross-check): Y2_B 1.97, Y3_B 1.23, Y4_P 0.30, Y5_P 0.53 nb;
  - SMASH visible: C4 + C5 in P 2.001 nb, C2 + C3 in B 9.008 nb.
- **kappa twice:**

  | fit | kappa | kappa_B |
  |---|---|---|
  | F-main, constrained | 0.83 +- 0.14 (stat) +- 0.44 (syst) +- 0.27 (MC) | 1.95 +- 0.17 +- 0.91 |
  | kappa free | 0.33 +- 0.15 +- 0.49 | 1.90 |
  | E1 + E2 | 0.50 +- 0.05 +- 0.15 | 1.11 +- 0.02 +- 0.05 |

  - F-free E1 + E2: alpha_3/alpha_2 = 0.51, alpha_4/alpha_2 = 0.47.
  - No data-driven conversion scale is possible (WP5).

### 2.5 Fig. P1 (figs_all/wp8_P1.{pdf,png}), draft caption
*WP8 summary. Data S0, isBest >= 0, L = 6079.2 nb^-1; simulation SMASH x PDG BR, V3 weights, M4' = 0.99615 M4.*
- (a) eta scale sigma(pp -> eta X)/sigma_SMASH, including C_sec(1 - f_eta) = 0.981, per estimator. Error bars: thick = stat,
  medium = + syst (Table S1 with S12 of WP8), thin = + MC + normalisation. Triangles: with the QED gamma* decay angle of C2
  (toy, closure).
- (b) BR(eta -> pi+pi-e+e-) from F-gold, F-main and F-free (cross-checks, WP4 fail) against PDG (band).
- (c) sigma(pp -> eta X) from the side bands and from F-main, against SMASH. The exclusive pp -> pp eta points are shown for
  orientation only: HADES 3.5 GeV (preliminary 211 +- 47 ub; thesis 136.9 ub) and models at sqrt(s) = 3.46 GeV (81-137 ub).
  [v8: panel (c) lacks the published INCLUSIVE HADES point at 3.5 GeV (1.14 +- 0.20 / 0.93 +- 0.14 mb). WP10 should add it
  or quote it in the caption. The figure file was not changed.]
  Dotted: the 1D (pT, y) reweightings.
- (d) eps2 check: the share of the ECAL-tagged eta per M(ee) class, data (fit), the emulated simulated tag (same fit), the
  simulated truth C2 + C3, and the generator fraction before acceptance.
- (e) Decay angle of pi0 -> gamma e+e- in PLUTO (lepton set, 124k decays) against QED (3/8)(1 + cos^2) and the isotropic decay
  of HGeant (class C2).
- (f) Data/SMASH of the peak yield in 3 x 3 (pT4, y4) cells (side-band estimator per cell).

### 2.6 Literature
No published inclusive sigma(pp -> eta X) near p_lab = 5.4 GeV/c was found, so there is **no independent comparison**.
[v8 correction: the nearest published INCLUSIVE value exists, at a lower energy (T = 3.5 GeV, sqrt(s) = 3.18 GeV):
sigma(pp -> eta X) = 1.14 +- 0.20 mb (PYTHIA extrapolation), 0.93 +- 0.14 mb (UrQMD), G. Agakishiev et al. (HADES),
Eur. Phys. J. A 48 (2012) 64, arXiv:1112.3607, Table 1. The cross section rises with energy, so this is a soft lower
reference for 4.5 GeV. The result 0.43 mb is a factor 2.2-2.7 below it (about 3 sigma). See "Verification", V6.]
Context only:

| reference | content | value |
|---|---|---|
| K. Teilab (HADES), Int. J. Mod. Phys. A 26 (2011) 694, arXiv:1009.3442 | exclusive pp -> pp eta, T = 3.5 GeV, preliminary | 211 +- 47 ub |
| K. Teilab, PhD thesis, Goethe University Frankfurt (2011) | the same | 136.9 +- 0.9 +- 10.1 ub |
| P. Lebiedowicz, arXiv:2501.14952 [hep-ph] (2025), Table 1 | exclusive pp -> pp eta model predictions at sqrt(s) = 3.46 GeV (this energy), models 1 / 2 / 3 | 136.9 / 116.2 / 80.8 ub |
| J. Staudenmaier et al., Phys. Rev. C 98 (2018) 054908 | SMASH produces "too many eta mesons from pp -> pp eta for sqrt(s) > 3.25 GeV" | qualitative |

- Staudenmaier et al. agree qualitatively with 0.41 of SMASH.
  [v8 correction: the quoted sentence (verified, arXiv:1711.10297 p. 4) is about the EXCLUSIVE channel pp -> pp eta in the
  2018 SMASH version. It says nothing about the inclusive cross section and does not support the value 0.41.]
- G. Alexander et al., Phys. Rev. 154 (1967) 1284 (5.5 GeV/c bubble chamber): the abstract has no eta value, and no number was
  taken from it.
- An unpublished, preliminary HADES inclusive eta -> gamma gamma analysis of the same beam time exists. It is a talk by
  A. Strach, Kraków, 20 Sep 2024, and gives no final number. It is not used.
  [v8 correction: the talk exists ("Inclusive production of eta and omega in pp@4.5 GeV with HADES", Workshop at 1 GeV
  scale: From mesons to axions, Jagiellonian University, 20 Sep 2024, indico.meson.if.uj.edu.pl/event/5/contributions/430).
  Its slide 13 DOES show preliminary numbers: sigma_eta = 1.15 mb (+-1 sigma window, 176 000 eta) and 1.72 mb (+-2 sigma
  window, 226 000 eta), without uncertainties (eta -> pi+pi-pi0, pi0 -> gamma gamma in the ECAL). Unpublished, not a
  comparison value, but it must not be described as "no number": it is 2.7-4 times the result here.]

---------------------------------------------------------------------------------------------------
## 3. Figures and tables
- `docs/latex/figs_all/wp8_P1.{pdf,png}` (Fig. P1, caption 2.5).
- `docs/latex/tables_all/`:
  - `wp8_tableP1.{txt,tex}` (Table P1, with the pion vertex);
  - `wp8_eps2.{txt,tex}`;
  - `wp8_S12.txt`;
  - `wp8_tableP2.{txt,tex}`;
  - `wp8_qed.txt`, `wp8_toy_angle.txt`, `wp8_pty.txt` (inputs).
- The TeX compiles (tested with the sim-report macros \piz \pipi \ppee \Mee \MeVcc \nbinv; no overfull boxes).

## 4. Acceptance criteria and fall-backs
| item | value | verdict |
|---|---|---|
| sim efficiencies reproduce sim study (pre, all M4) | 5 classes within 0.1 % | PASS |
| WP7 chain reproduced in wp8_syst.C | a_eta, BR_G, Y_bump, Y_peak, DR, sigma exact | PASS |
| eps2: generator C2 pair-mass fraction vs KW x FF | 0.998 | PASS |
| eps2: decay angle | HGeant isotropic vs QED | FAIL -> S12 entry (sign known) and a variant, not applied |
| eps2: tagged share >= 20 MeV data/sim | 1.06 +- 0.15 (fit), 0.95 +- 0.13 (truth) | PASS (MC-limited) |
| pion vertex data/sim | 1.00 +- 0.06 (peak), 1.00 +- 0.014 (tag) | PASS (2.5 % kept) |
| BR vs PDG | BR_G -0.0 sigma, BR_A +0.26 sigma (total) | consistent (cross-checks) |
| R_34 | 0.92 +- 0.34 (E1+E2) | not fired |

Fall-backs:
- Inherited and respected: the WP4 gate (F-main/F-gold/F-free cross-checks); E1 dropped; bump EXTRAPOLATED; BR_C withdrawn
  in E2; simulation-limited (MC > stat for BR_G, BR_A, a_eta, DR).
- Not fired: "C_sec / f_win not determinable" (the WP2 + AMENDMENT 1 values are used).

Deviations (flagged):
- (a) The 3 x 3 (pT4, y4) variant has 16 % noise. The 1D reweightings are added as the closest sound test and are reported
  outside the total.
- (b) The eps2 angle effect is a systematic of known sign, not a correction.
- (c) S12 combines independent sub-sources by the spec's maximum rule; the quadrature version is given.

## 5. Caveats
1. The toy acceptance is crude (sharp cuts, no tracking/PID efficiencies). The closure-based factors come from prompt
   SMASH pi0s with different kinematics. Together they bracket 0.86-0.91. A proper correction needs the pi0 rest frame per
   C2 row: store the true photon or pi0 momentum in the sim compact, or re-simulate C2 with the QED angle.
2. The HGeant source is the local hgeant2-7.4 copy. That it is identical to the production build is ASSUMED. The measured
   closure 0.891 matches the toy 0.880, which supports the assumption.
3. **Per-cell side bands use the global R' and leakage f.** A pT dependence of the background shape could mimic the pT
   effect of the 1D variant. The event-group factor 1.19 is applied uniformly.
4. **The bump (pT4, y4) proxy** uses reconstructed 4-body kinematics (the photon is missing).
5. **Fit-based simulated tag shares** have no replicas: their MC error is taken from the relative truth error. The emulated
   tag is UNVERIFIED (WP6).
6. **The generator (pT, y)** comes from 40 of 100 000 standard-set files (94 517 etas).
7. **PDG BRs** are typed as in sim_eta_physics.h, not re-checked online.
   [v8: checked. They are the PDG 2020 values. The current edition (PDG 2026) has eta -> pi+pi-pi0 23.02 +- 0.25 %,
   pi+pi-gamma 4.28 +- 0.07 %, pi+pi-e+e- 2.76 +- 0.16 e-4 (S = 1.5); pi0 -> e+e-gamma is unchanged (1.174 +- 0.035 %).
   Effect: sigma -2.2 %, BR_G and BR_A +0.4 %, DR -1.8 %. See "Verification", V3.]
8. **sigma x BR(pi+pi-pi0)** from B4 is extrapolated (168 %): information only.

## 6. Exact commands (from /home/przygoda/HADES/FAT/FAT; ACLiC dir aclic_wp8)
    python3 -I plots/all_report/wp8_qed.py
    S="ssh -o RemoteCommand=none -o RequestTTY=no -o BatchMode=yes vae25 bash -s"; R=/hdd3/przygoda/HADES/FAT_out/all_report/wp8/remote
    $S < $R/remote_eta_ptY.sh > $R/gen_eta_ptY_std40.txt              # read-only, ~15 min
    $S < $R/remote_pi0dalitz_angle.sh > $R/pluto_pi0dalitz_lep3.txt   # read-only, ~10 min
    B='gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_wp8", kTRUE)'
    root -l -b -q -e "$B" 'plots/all_report/wp8_sim.C+O'     # 30 s, ~9 GB: wp8/{wp8_sim.root, wp8_sim.txt, wp8_grid.txt}
    python3 -I plots/all_report/wp8_toy_angle.py              # 7 s
    python3 -I plots/all_report/wp8_pty.py                    # 5 s: tables_all/wp8_pty.txt, wp8/wp8_pty.json
    root -l -b -q -e "$B" 'plots/all_report/wp8_syst.C+O'    # 16 s: wp8/wp8_syst_raw.txt
    python3 -I plots/all_report/wp8_physics.py                # tables, Fig. P1, wp8/results.txt

---------------------------------------------------------------------------------------------------
## Verification (independent verifier, prefix v8, 2026-11-14)

**Verdict.** The arithmetic of WP8 is CORRECT: every physics number is reproduced from the primary inputs (WP5 fit parameters,
WP6 yields and SMASH window sums, WP7 family deltas) with my own code, to 1e-5 or better. The WP8 finding on the gamma* decay
angle of class C2 is CONFIRMED (source and an independent toy). Several factors and parts of the error budget are NOT
correct, and the literature comparison was incomplete in a way that hides a large tension:
1. The branching ratios are PDG 2020, not the current edition (PDG 2026). Corrected: sigma -2.2 %, BR_G and BR_A +0.4 %, DR -1.8 %.
2. sigma from the peak side bands has no systematic for the C4 : C5 composition of the peak (conversion scale of C4).
   With the analysis' own prior kappa = 1 +- 0.2 this is +-9.6 %. Syst 16.0 % -> 18.5 %.
3. Small factor errors: SMASH sigma of the simulated sample is 1.069 mb, not 1.0665 (+0.28 %); C_sec of the peak is 1.0116,
   not the bump/peak mean 1.0110 (+0.06 %).
4. Three one-sided normalisation items are missing from Table P2 (START-based sorter loss, cor_ok rule, event-vertex
   efficiency). All would raise sigma.
5. Literature: all six WP8 citations exist and their numbers are right. But the published INCLUSIVE HADES value at 3.5 GeV
   (1.14 +- 0.20 / 0.93 +- 0.14 mb) was missed, the Staudenmaier sentence was over-interpreted, and the Krakow talk does give
   preliminary numbers (1.15 / 1.72 mb). The corrected result is about 3 sigma below the published 3.5 GeV inclusive value,
   although the cross section should rise with energy.

**Corrected headline (v8):** sigma(pp -> eta X) = **0.429 +- 0.036 (stat) +- 0.079 (syst) +- 0.025 (MC) +0.035/-0.037 (norm.) mb**
(total 22.9 %), data/SMASH = 0.401 +- 0.092. The one-sided items are listed apart. **Recommendation: label the absolute
sigma PRELIMINARY**, because the event-vertex efficiency of the mandatory vz > -500 mm cut is not measured for the eta signal
(V7). Ratios (BR_G, BR_A, DR) are not affected by this.

Files: `docs/latex/tables_all/v8_verification.txt` (all numbers), `v8_tableP2_corr.tex` (corrected sigma, one-sided items,
literature; compiles, no overfull box), `/hdd3/przygoda/HADES/FAT_out/all_report/v8/v8_verify.log`. No WP8 output was
overwritten.

### V1. Reproduction with my own arithmetic (plots/all_report/v8_verify.py; WP8 value in brackets)
| quantity | formula | v8 | WP8 |
|---|---|---|---|
| BR_G | 2.68e-4 x alpha_5/alpha_2 = 2.68e-4 x 0.38457/0.38825 | 2.6546e-4 | 2.6546e-4 |
| BR_A | 2.68e-4 x rho_d (1.31520) | 3.5247e-4 | 3.5247e-4 |
| a_P | Y_peak / SMASH (C4 + C5) in P = 5087.97 / 12164.94 | 0.418249 | 0.418249 |
| a_B (all M(ee)) | 16389.3 / 54758.9 | 0.299299 | (DR input) |
| a_B (E2) | 7822.17 / 25097.0 | 0.311677 | 0.311677 |
| DR | a_P / a_B(all) | 1.39743 | 1.39743 |
| DR_F | (5044.8/19426.1) / (12164.9/54758.9) | 1.16896 | 1.16896 |
| sigma | a_P x 1.0665 x 1.0110 x 0.97 | 0.43744 mb | 0.43744 |
| sigma_F | a_eta x 1.0665 x 1.0110 x 0.97 | 0.39982 mb | 0.39982 |
| sigma_vis peak / bump | Y / 6079.2 | 0.83695 / 2.6960 nb | 0.837 / 2.696 |
| sigma x BR | ppee / ppg / ppp0 (a_B E2) / ppee F-main | 117.2 / 18 460 / 74 714 / 140.9 nb | same |
| units | raw side-band limit 5008.5 x 64/63; raw / (SMASH x 63/64) gives the same a_P | 0.418249 | |

- The WP5 inputs (alpha_2, alpha_5, rho_d, a_eta, F-main yields) were themselves reproduced by v5 with independent code.
- Table P1 efficiencies agree with the WP6 window sums: C5 in P 6360.13 rows / 6079.2 = 1.04621 nb = 285.822 nb x 0.366036 %;
  C4 0.954863 nb; C2 + C3 in B 9.00759 nb.

### V2. Factors
| factor | check | verdict |
|---|---|---|
| 64/63 | data yields x 64/63, templates in raw rows = nb x 6079.2 x 63/64; a_P identical both ways | PASS |
| L = 6079.2 nb^-1, 6.9 % | WP2 6079.195 (analysed runs), A3 note | PASS |
| L_sim = 221.05 nb^-1 | 1e10 / 45.239 mb; 6079.2/221.05 x 63/64 = 27.0718 | PASS |
| SMASH sigma(eta) | 1.0665 +- 0.0016 mb is a 200-file sub-sample; the simulated sample has N_eta = 2.364e8 in 1e10 events, i.e. 1.0694 mb (the sim study's 1.069; WP8's own 40-file check 1.0690 +- 0.0035) | MINOR ERROR: sigma +0.28 %; Table P1 efficiencies 0.28 % too high |
| C_sec | sigma comes from the peak only: C_sec(peak) = 1.0116, not the mean 1.0110 | MINOR ERROR: +0.06 % |
| (1 - f_eta) = 0.97 +0.03/-0.04 | AMENDMENT 1 used, never (1 - f_win); C_sec (1 - f_eta) = 0.981 +3.1/-4.1 % | PASS |
| live time | one-sided, separate, correct sign (L contains the PT2 live fraction; PT3 is lower). +7.3 % = 1/0.932 (rounded "overall"); the L-weighted LF3/LF2 x 64/63 = 0.9334 (5675.7/6080.5) gives +7.1 % | PASS (quote +7.1-7.3 %) |
| BR values | see V3 | FAIL (outdated) |

### V3. PDG branching ratios (confirmed error, corrected)
- Verified on line (pdg.lbl.gov listings): the simulation and WP7/WP8 use **PDG 2020** (P.A. Zyla et al., PTEP 2020, 083C01):
  eta -> pi+pi-pi0 22.92 +- 0.28 %, pi+pi-gamma 4.22 +- 0.08 %, pi+pi-e+e- 2.68 +- 0.11 e-4, pi0 -> e+e-gamma 1.174 +- 0.035 %.
- Current edition **PDG 2026** (F. Takahashi et al. (PDG), Int. J. Mod. Phys. A 41, 2630011 (2026); listing of 2026-06-01):
  23.02 +- 0.25 %, 4.28 +- 0.07 %, **2.76 +- 0.16 e-4** (S = 1.5), 1.174 +- 0.035 %. (PDG 2022 and 2024 already had 23.02 and 4.28.)
- Rescaling of the simulated classes (no refit needed): C2, C3 x 1.00436; C4 x 1.01422; C5 x 1.02985. SMASH peak in P x 1.0224.
  - sigma: -2.19 %. a_P 0.41825 -> 0.40909.
  - BR_G and BR_A scale with BR(eta -> pi+pi-pi0) of the C2 normalisation (the 2.68e-4 factor cancels): 2.655 -> 2.666e-4 and
    3.525 -> 3.540e-4 (the F-main value is approximate, a_eta is C2/C3-dominated).
  - DR: 1.397 -> 1.373; DR_F: 1.169 -> 1.148.
  - The PDG reference for the BR comparison becomes 2.76 +- 0.16 e-4. Pulls stay below 0.1 sigma (total errors of 90-130 %).
  - N row: peak BR mix 2.33 -> 3.15 % (C4 share 0.473); N(sigma) +7.92/-8.38 -> +8.20/-8.64 %; bump/C2 BRs 3.22 -> 3.17 %;
    DR N 3.97 -> 4.47 %.

### V4. Error propagation and uncertainty classes
- **Reproduced:** N(sigma) = L 6.9 % (+) C_sec (1 - f_eta) +3.11/-4.14 % (+) BR mix 2.33 % = +7.92/-8.38 % (WP7 +7.915/-8.357);
  total 20.8 % (WP8 20.7 %); BR_G N 3.22 % (pi+pi-pi0 (+) pi0 Dalitz); DR N 3.98 %; visible cross sections carry L only.
  stat, syst, MC, N and live time are kept apart everywhere. L, sigma_SMASH and C_sec cancel in BR_G, BR_A, DR. PASS.
- **Missing systematic (error, added): composition of the peak.** a_P = Y_peak / (C4 + C5)_SMASH equals the eta scale only if
  the conversion part C4 (48 % of the simulated peak) is simulated correctly. Nothing in Table S1 varies it for the side-band
  estimator (the leakage-shape variant changes only the leakage, +0.2 %). With the analysis' own prior kappa = 1 +- 0.2
  (WP3/WP5 fall-back), a = a_P (C4 + C5)/(kappa C4 + C5) moves by +10.5/-8.7 %, i.e. **+-9.6 %**.
  - Added to S12 (efficiency model) with the spec's maximum rule: syst(sigma) 16.0 -> **18.5 %** (18.6 % as a separate family).
  - DR is almost immune (C3 in the bump moves the same way: about 1 %). BR_G is not affected (C4 = 0 at M(ee) >= 20 MeV).
  - **Data indication (not adopted):** the side-band yields per M(ee) class, E1 2380 +- 334 and E2 2780 +- 372 (stat (+) MC),
    against the SMASH shapes (C4 E1/E2 4946/942, C5 2201/4349 rows), give with rho = 1 (PDG BR): a = 0.59 +- 0.10, kappa = 0.37,
    i.e. sigma = 0.62 +- 0.10 mb. The data peak is more E2-like than SMASH (share 0.546 +- 0.053 vs 0.423). E1 is the class
    with the unmodelled component (k_E1 = 7.6), so this is only an indication, but it points in the same direction as the
    F-main kappa (0.83 +- 0.30 with the prior, 0.33 +- 0.52 free): the a_P route would then UNDER-estimate sigma.
- **One-sided items missing from Table P2 (all would raise sigma; listed apart, not added):**
  - START-based track-sorter loss: +1.3 to +2.8 % (report_exp Table syst; L is corrected for it, the yields and the
    simulation are not).
  - cor_ok rule: +2.2 +- 1.0 % (WP7 fallback test: Y_peak +111 rows, z = 2.3). The data drop cor_ok == 0 rows; the
    simulation has none (checked: 0 of 230 469 std and 7 441 490 lep rows of sim_compact.root have cor_ok == 0), so it does not
    lose them.
  - Event-vertex efficiency of vz > -500 mm: see V7.
- **pT-only reweighting (+19/-7 %) outside the total:** WP8 calls this the analyst's decision; no such decision is recorded in
  the spec or the notes. It is WP8's choice. I list it next to the total as an asymmetric production-model item.
- **Known-sign effects quoted as symmetric:** the QED-angle effect (BR_G -9 %, a_eta +12 %, DR -10 %) is a bias of known sign.
  WP8 states this; acceptable under the frozen procedure.

### V5. Gamma* decay angle of class C2 (WP8 finding 1): CONFIRMED
- Source: hgeant2-7.4/geant/gdalet.F lines 81-94. YINV = beta x (2 RNDM - 1) is uniform, and E2 = (E_gamma* + YINV p_gamma*)/2,
  so YINV = beta cos(theta*): cos(theta*) is uniform, i.e. the decay is isotropic. QED: 1 + cos^2 + (4 m_e^2/q^2) sin^2.
- Kroll-Wada fraction with M(ee) >= 20 MeV (own integration): 0.2679 (no form factor), 0.2693 (slope 0.0332). WP8: 0.2679 / 0.2693.
- **Own toy** (`v8_angle_toy.py`, written from scratch, 6e6 decays; eta m_T-exponential, flat Dalitz plot, 18-85 deg,
  p_pi > 0.10 GeV/c): eps(QED)/eps(iso)

  | p_e min [GeV/c] | all | < 10 MeV | >= 10 | >= 20 | 20-30 |
  |---|---|---|---|---|---|
  | 0.04 | 0.895 | 0.901 | 0.885 | 0.906 | 0.867 |
  | 0.06 | 0.875 | 0.886 | 0.860 | 0.878 | 0.840 |
  | 0.10 | 0.856 | 0.869 | 0.836 | 0.847 | 0.801 |

  WP8 toy at 0.06: 0.880 (all), 0.886 (>= 20). Agreement within 1 %, with different eta kinematics. The size of the effect
  (11-13 %) and the S12 entries of WP8 stand.

### V6. Literature (every citation checked on line, 2026-11-14)
| reference | exists | number in WP8 | verdict |
|---|---|---|---|
| K. Teilab (HADES), Int. J. Mod. Phys. A 26 (2011) 694, arXiv:1009.3442 | yes | 211 +- 47 ub (prelim.), pp -> pp eta, 3.5 GeV | correct (text: "(211+-47) ub (prelim.)") |
| K. Teilab, PhD thesis, Goethe-Universitaet Frankfurt (2011) | yes (hades.gsi.de thesis PDF) | 136.9 +- 0.9 +- 10.1 ub | correct (abstract and sect. on eta) |
| P. Lebiedowicz, arXiv:2501.14952 [hep-ph], Table 1 | yes (v2, 2025; no journal reference yet) | 136.9 / 116.2 / 80.8 ub at sqrt(s) = 3.46 GeV, models 1/2/3 | correct; 136.9 coincides numerically with the thesis value (different quantities) |
| J. Staudenmaier et al., Phys. Rev. C 98 (2018) 054908, arXiv:1711.10297 | yes | quote "too many eta mesons ... pp -> pp eta for sqrt(s) > 3.25 GeV" | quote correct; the conclusion "agrees with 0.41 of SMASH" is NOT supported (exclusive channel, 2018 SMASH) -> corrected in Sect. 2.6 |
| G. Alexander et al., Phys. Rev. 154 (1967) 1284 | yes (pp at 5.5 GeV/c) | no eta value | correct (INSPIRE abstract: omega 0.11 +- 0.02 mb, no eta) |
| A. Strach, talk, Krakow, 20 Sep 2024 | yes (indico contribution 430) | "gives no final number" | INCOMPLETE: slide 13 shows 1.15 mb (+-1 sigma window) and 1.72 mb (+-2 sigma), preliminary, no errors -> corrected in Sect. 2.6 |

- No citation had to be removed: all could be verified.
- **Missed reference (added):** G. Agakishiev et al. (HADES), "Inclusive dielectron spectra in p+p collisions at 3.5 GeV kinetic
  beam energy", Eur. Phys. J. A 48 (2012) 64, arXiv:1112.3607, Table 1: **sigma(pp -> eta X) = 1.14 +- 0.20 mb** (PYTHIA
  extrapolation) and **0.93 +- 0.14 mb** (UrQMD), T = 3.5 GeV (eta Dalitz decay, model-dependent extrapolation).
- **Tension.** The corrected 0.429 +- 0.098 mb is 0.38 / 0.46 of these values (3.2 / 2.9 sigma with the total errors), although
  the inclusive eta cross section rises from 3.5 to 4.5 GeV (Cassing-Sibirtsev-Mosel parametrisation, shown on Strach's slide 15).
  SMASH (1.07 mb) is close to the measured inclusive value; the statement "SMASH produces too many eta" is supported only for
  the exclusive channel. The preliminary same-beam-time ECAL analysis (1.15-1.72 mb) points the same way. The report must state
  this tension and must not present 0.41 x SMASH as confirmed by the literature.

### V7. Event-vertex efficiency of the mandatory vz > -500 mm cut (open; the prime suspect for the tension)
- The cut stays (analyst decision). What is missing is its EFFICIENCY for LH2 eta events in data compared with the simulation.
- Data (report_exp Sect. 3): the cut keeps 0.565 of the unlike-sign rows in the peak window (0.62-0.72 per day); with the
  like-sign survival 0.255/0.202 and CB/N 0.103 after the cut, about 0.68 of N - CB. The eta content of the removed rows is
  not measured (report_exp says so).
- Simulation: FAT rows / PAT rows (isBest >= 0, vz > -500, PT3) = 0.968 (std), 0.976 (lep), 0.959 (eta set) (sim report
  Tab. fat), so the simulated vertex survival is >= 0.96.
- If the eta events lost their vertex at the N - CB rate, sigma would rise by x 1.4 (x 1.7 at the raw row rate). This alone
  would bring 0.43 mb to 0.6-0.7 mb; with the other one-sided items (+11-12 %) to 0.7-0.8 mb.
- Weak counter-indication: the per-period eta scale does not follow the per-period vertex survival (a_eta x K: III/I =
  0.98 +- 0.05, expected about 1.10 if it tracked 0.655 -> 0.72). This tests only the period-dependent part, not a constant offset.
- Not determinable within WP8/v8 (needs a PAT-level pass of the rows without a vertex, as a diagnostic only, or the simulated
  vertex efficiency for eta events checked against data). Fall-back of the spec, applied by analogy: "absolute values
  PRELIMINARY with a stated extra uncertainty; ratios unaffected".

### V8. Corrected numbers (v8; v8_verification.txt, v8_tableP2_corr.tex)
| quantity | WP8 | v8 | change |
|---|---|---|---|
| sigma(pp -> eta X) [mb] | 0.437 +- 0.037 +- 0.070 +- 0.026 +0.035/-0.037 | **0.429 +- 0.036 +- 0.079 +- 0.025 +0.035/-0.037** | PDG 2026, sample sigma, C_sec(peak), peak composition |
| total | 0.091 (20.7 %) | 0.098 (22.9 %) | |
| data / SMASH | 0.410 +- 0.085 | 0.401 +- 0.092 | |
| sigma x BR(pi+pi-e+e-) | 117 nb (PDG 2.68e-4) | 118.5 nb (PDG 2026 2.76e-4) | |
| sigma x BR(pi+pi-gamma) | 18.5 ub | 18.4 ub | |
| BR_G [1e-4] | 2.65 (vs 2.68 +- 0.11) | 2.67 (vs PDG 2026 2.76 +- 0.16) | errors unchanged in per cent |
| BR_A [1e-4] | 3.52 | 3.54 | |
| DR (side bands) / DR_F | 1.397 / 1.169 | 1.373 / 1.148 | N 4.5 % |
| one-sided (not in total) | live time +7.3 % | live time +7.1-7.3 %, sorter +1.3-2.8 %, cor_ok +2.2 +- 1.0 %, event-vertex efficiency not measured (up to x 1.4-1.7), (pT, y) +19/-7 % | |

### V9. Checks of the verification (PASS/FAIL)
| check | value | verdict |
|---|---|---|
| BR_G, BR_A, DR, DR_F, sigma, sigma_F, sigma x BR, sigma_vis reproduced | to 1e-5 | PASS |
| 64/63, L, L_sim | consistent | PASS |
| SMASH sigma, C_sec | 1.0665 -> 1.0694; 1.0110 -> 1.0116 | MINOR ERRORS, corrected (+0.34 % together) |
| f_eta instead of f_win | AMENDMENT 1 respected | PASS |
| live-time statement | one-sided, separate, +7.1-7.3 % | PASS |
| PDG BRs | PDG 2020 used; PDG 2026 differs | FAIL -> corrected |
| N-row propagation | +7.92/-8.38 % | PASS |
| completeness of the sigma error budget | peak composition missing; 3 one-sided items missing | FAIL -> added |
| gamma* angle finding | source + own toy (0.875 vs 0.880) | PASS (confirmed) |
| literature: citations exist, numbers right | 6 of 6 | PASS |
| literature: completeness and interpretation | published inclusive 3.5 GeV value missed; Staudenmaier over-interpreted; Strach numbers omitted | FAIL -> corrected |

**Fall-backs.** Fired by analogy (recommended): "absolute values PRELIMINARY with a stated extra uncertainty" for sigma (V7).
WP8's fall-backs (WP4 gate, E1 dropped, bump EXTRAPOLATED, BR_C withdrawn, simulation-limited) are confirmed.

### V10. For WP10
- Use `tables_all/v8_tableP2_corr.tex` for the sigma line, the one-sided items and the literature block, or copy its numbers
  into Table P2.
- Quote BR_G/BR_A against PDG 2026 (2.76 +- 0.16 e-4) and state which PDG edition the simulation used (2020).
- State the tension with EPJA 48 (2012) 64 in the conclusions and the abstract of the absolute result.
- Fig. P1 panel (c): add the inclusive 3.5 GeV point (1.14 +- 0.20 / 0.93 +- 0.14 mb) or name it in the caption.

### V11. Changes made
1. `wp8_results.md`: [v8] markers in the headline (sigma row), finding 2 (pT decision), Sect. 2.5 (Fig. P1 (c)), Sect. 2.6
   (three corrections), caveat 7; this section appended. No WP8 number was changed in place.
2. New files (all v8): `plots/all_report/v8_verify.py`, `plots/all_report/v8_angle_toy.py`,
   `docs/latex/tables_all/v8_verification.txt`, `docs/latex/tables_all/v8_tableP2_corr.tex`,
   `/hdd3/przygoda/HADES/FAT_out/all_report/v8/v8_verify.log`; a read-only ROOT check of sim_compact.root (cor_ok, vz) from
   the scratchpad (`v8_corok.C`). Downloaded literature in the scratchpad folder `v8_web/` (Teilab, Lebiedowicz, Staudenmaier,
   EPJA 48 (2012) 64, PDG listings 2020-2026, Strach slides).
3. Not changed: WP8 tables, Fig. P1, results.txt; WP0-WP7 files; ScanEngine.h, scan_*.C; FAT_simwt. No git command.

### V12. Commands (from /home/przygoda/HADES/FAT/FAT)
    python3 -I plots/all_report/v8_verify.py > /hdd3/przygoda/HADES/FAT_out/all_report/v8/v8_verify.log   # < 1 s
    python3 -I plots/all_report/v8_angle_toy.py                                                            # about 1 min
    root -l -b -q <scratchpad>/v8_corok.C                                                                  # read-only, 1 min
