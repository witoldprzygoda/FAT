> **AMENDMENT 1 (lead analyst, 2026-11-14 08:20; binding, overrides every f_win formula below for eta quantities).** WP2 measured the material (empty-target) share of the N - CB rows, f_win = 16.8 % (B+P; 14.7 +- 1.1 +1.0/-2.0 % after the halo reweighting of the verification). Two independent verifiers showed that the eta SIGNAL has almost no material component: **f_eta = 3 +4/-3 %** (vz-slice and cell-centre tests, 16.8 % excluded at 3.4 sigma). Therefore: sigma(pp -> eta X) = a_eta x 1.0665 mb x C_sec x (1 - f_eta), with C_sec x (1 - f_eta) = 0.98 +0.03/-0.04; every eta yield, BR or cross section uses (1 - f_eta), never (1 - f_win) (that would remove the material twice, because the non-eta background model already absorbs it). f_win is quoted only as the composition of the non-eta background (data background vs SMASH: about 13-15 % of N - CB in 0.40-0.70 is material SMASH does not simulate). Per-period stability factors: K = C_sec(period) x 0.97. Details: work_all/wp2_results.md section 7 and verify_wp2/vw2a_note.md, vw2b_note.md.

# SPEC report_all: data described by simulation, eta bump and peak, X17 limit

**Basis.** Proposal 3 (ratios and 2D fits), cut down to Proposal 1's layout of 10 histograms (both judges).
- From P1: per-row mass transfer, stress closure, Barlow rule, systematic families, folds, CR-K null test.
- From P2: table of corrections in plain words, CR-L and displaced-window closure, fallbacks declared in advance, counted "after" display, smoothed background, three X17 normalisations, toy gate.
- Every "must fix" item of both judges is applied (map in Sect. 9).

**Names.**
- Signal classes: C2 bump pi0 Dalitz, C3 bump conversion, C4 peak conversion (eta -> pi+pi-gamma), C5 direct eta -> pi+pi-e+e-, C6 X17.
- T1-T4 and G_* as in the sim report. G_D = G_pi0D + G_vec + G_oth; G_C = G_conv.

**Locations.**
- Code: /home/przygoda/HADES/FAT/FAT/plots/all_report/. Untracked; no commits, no branches.
- Copies of SimEngine.h, sim_eta_physics.h, s5_blind_ana_fit.h, s5_blind_ana_tmpl.h and sim_truth.h from FAT_simwt, with the includes and kFigDir/kTabDir redirected. Never run the sim macros in place: they write into the data repository.
- Large outputs: /hdd3/przygoda/HADES/FAT_out/all_report/.
- LaTeX: docs/latex/report_all.tex, report_all_short.tex, figs_all/, tables_all/.

## 1. Analyst's decisions (listed in both reports, open for review)

| # | decision | reason |
|---|---|---|
| D1 | Nominal S0 = pre + both pions in the vertex fit, for bump, peak and X17. "pre" is a full parallel column. | CB/N 0.13 -> 0.06; S0 is the frozen X17 selection; sim F/F_best 0.95-0.98 for all targets; efficiencies cancel in ratios. |
| D2 | isBest >= 0 is nominal. isBest == 1 gives data-only yield ratios. | Analyst rule; sim isBest is unrealistic in the lepton set. |
| D3 | No upper cuts on MM4, M(pipi), alpha_eta-rest or OA. | Sculpting and band-pass effects were measured in data. |
| D4 | DIRECT (dz_pi < 40 mm) is a variant. PEAK (+ MM4 < 2.85) is a variant only after the fake-peak test, with sim MM4 recomputed at 4500 MeV. | MM4 sculpts; the simulation was generated at 4500 MeV. |
| D5 | CB = 2 sqrt(N++ N--), k = 1, never from simulation. Nuisances: k_E2 = 1 +- 0.05; k_E1 (M(ee) < 10 MeV) free. | Analyst rule; there is no 4-body k; the CB at small OA is only a lower limit. |
| D6 | Simulation = std + lep + eta, owned rows (mode A), w x V3 (kPhysMain 108), T1 + T2, every channel at SMASH sigma_prod x PDG BR. Ghost variant (a) is a systematic and is shown in "before". | Rules v2; the ghost ruling is still pending. |
| D7 | Data unit: (N - CB) x 64/63 / L, L = 6079.2 nb^-1. Sim unit: sum w / 221.05 nb^-1. Old tables x 1.0161. | Analysed runs; PT3 prescription; the old tables used 6080.5 and no 64/63. |
| D8 | Sim M4' = s [m4_true + r (m4 - m4_true)], s = 0.99615, r = 1; sim M(pipi) x 0.9952; M(ee) not scaled. | The peak shift (0.5439 vs 0.5460) is significant; the width difference (5.44 +- 0.60 vs 6.09 MeV) is 1.1 sigma. |
| D9 | Control region = the eta-free part of the M4 x M(pipi) plane (CR-U, CR-K). CR-L and the displaced windows are tests only. | M4 alone has no signal-free region below the bump; CR-K is kinematically bump-free. |
| D10 | Only the non-eta background is renormalised. Every eta component scales with a_eta; G_eta also gets g = 1 +- 0.30. | Indications: eta 0.3-0.5, background 1.3-1.4 of SMASH; G_eta is 6.4 % of CR-U and 37 % of the simulated peak. |
| D11 | Background: per-slice scale c_j x common tilt (1 + b x), x = (M4 - 0.55)/0.1. A quadratic term, then per-slice tilts, are added only if delta chi2 > 4 in stage A. | chi2/ndf 3.8 with a scale only, 1.3 with a tilt; c_j absorbs M(pipi) mismodelling. |
| D12 | Main extraction: Poisson likelihood of the unsubtracted N, 5 M(pipi) slices x 2 M(ee) classes, M4 0.400-0.700, fitter s5_blind_ana_fit.h (IRLS, BB-lite, sandwich errors). | The moving edge and the M(ee) split break the 1D tilt/bump degeneracy (a_bump 0.41 -> 0.03). |
| D13 | Non-eta templates smoothed (exp(pol3) per histogram). Signal templates raw with BB-lite. MC error quoted separately. | One std row = 27.5 data rows; G_conv N_eff is 18 in the X17 window. |
| D14 | Side-band estimators (B4, peak side bands) are reported next to the fit, not combined. P_main, E_step and the ECAL tag are cross-checks only. | Transparency; the sim ECAL tag is emulated (UNVERIFIED) and the photon scale is 8 % high. |
| D15 | BR(eta -> pi+pi-e+e-): primary = golden ratio (M(ee) >= 20 MeV, relative to the C2 bump); second = rho_d; third = conversion route. | The golden ratio is free of L, sigma_eta, PT3, live time and material. |
| D16 | sigma(pp -> eta X) is secondary and needs WP2. | L is LH2-only; the sector-2 deficit is not in L. |
| D17 | Live time is not applied; it is quoted as a one-sided +7.3/-0 % on cross sections. | Not known whether it applies to PT3 single-interaction events. |
| D18 | X17: 2D Poisson likelihood with the background free per M(ee) bin (from data). POI r = BR_X BR_ee / BR_dir in the same fit. 90 % CLs. | No simulated background per bin; nothing is used twice. |
| D19 | Blinding: the window count (5258, CB 112) and the M(ee) > 20 MeV peak share (0.306) are already known, and the reports say so. The 2 MeV M(ee) shape of the peak band in 13-27 MeV stays masked until x17_frozen.txt (procedure + Asimov limit) exists. | Honest blinding. |
| D20 | No data yield is shown before the estimator passes WP4. | The earlier blind test was not blind by construction. |

## 2. Common definitions (frozen in plots/all_report/defs.h)

- **pre** = isBest >= 0 && vz > -500 && PT3 (data trigbit == 8192, sim bit 8192) && vz > -160 && cor_ok == 1. **S0** = pre && pip_isUsedVertex && pim_isUsedVertex.
- **Inputs.**
  - Data: /hdd3/.../exp/pippimepem_isBestAll_T4530/scan_v2/scan_compact.root (us/pp/mm; it has dz_pi).
  - Sim: /hdd3/.../sim/pippimepem_betanew/scan/sim_compact.root + sim_compact_eta.root (+ rw_*.root).
  - Variances: event groups (grp).
- **Windows (data scale).**

  | region | M4 [GeV/c2] |
  |---|---|
  | bump B | 0.4550-0.5225 |
  | peak P | 0.53364-0.55416 |
  | peak side bands | 0.5250-0.5336 and 0.5542-0.5628 |
  | X17 counting cross-check | 0.5312-0.5566 (= sim 0.5333-0.5587 x 0.99615), M(ee) 15.5-24.5 MeV |
- **Bump edge, per row:**

      E(m, q) = 0.99274 sqrt(m_eta^2 - m_eta (E_pi0 - p_pi0)(1 - q^2/m_pi0^2))
      E_pi0 = (m_eta^2 + m_pi0^2 - m'^2)/(2 m_eta),   p_pi0 = sqrt(E_pi0^2 - m_pi0^2),   m' = m/0.9970

  - m = M(pipi) on the data scale, q = M(ee). There is no bump for m > 0.4116.
  - Checks: E = 0.5186 (m 0.30, q 0), 0.5008 (0.39, 0), 0.5096 (0.39, 60 MeV), 0.4739 (0.4116, 0).
  - At m = 0.30 and q = 100 MeV the edge is at 0.533, so C2 reaches under the peak; the templates carry it.
- **Layout F10.**
  - M(pipi) slices: S1 < 0.34, S2 0.34-0.37, S3 0.37-0.39, S4 0.39-0.43, S5 >= 0.43.
  - M(ee) classes: E1 < 10 MeV, E2 >= 10 MeV.
  - M4 bins: 5 MeV in 0.400-0.520, 2.5 MeV in 0.520-0.570, 5 MeV in 0.570-0.700 (70 bins).
  - **F5g** = the same slices, M(ee) >= 20 MeV only.
- **Control regions.**
  - CR-U: 0.565 <= M4 < 0.700, all slices.
  - CR-K: S5 with 0.400 <= M4 < 0.522. The limit is 0.4116; 0.43 adds 3 sigma of M(pipi) resolution.
  - CR-L (test only): 0.5225-0.5300.
- **CB.**
  - 2 sqrt(N++ N--) in 10 MeV M4 cells, interpolated, spread over the F10 cells by f(cell | 50 MeV band) taken from N++ + N--.
  - Error: 200 Poisson replicas of N++ and N--.
- **Model** (bin i of histogram h = (j, e)):

      mu_hi = k_e CB_hi + c_j (1 + b x_i)[G_D + kappa_B G_C]_hi + a_eta [C2 + kappa (C3 + C4) + rho_d C5 + g G_eta]_hi

  - Fits run on raw rows. Sim templates are in expected raw rows: sum w x (L/L_sim) x 63/64.
  - Free: a_eta, kappa, rho_d, c_1..c_5, b, kappa_B, k_E1. Constrained: g, k_E2.
  - **F-free**: alpha_2..alpha_5 replace a_eta, kappa and rho_d.
- **Yields** (rows x 64/63): Y_c = factor x C_c, in its window and over all M4; Y_bump = Y2 + Y3 in B; Y_peak = Y4 + Y5 in P; DR = (Y_peak/Y_bump)_data / (Y_peak/Y_bump)_sim.

## 3. Work packages

### WP0 Histogram factory (the only heavy pass; gate for everything else)
- **Output.** /hdd3/.../all_report/hists_v1.root, one directory per (selection x variant), sim memory 8-10 GB.
- **Layouts.**
  - F10 and F5g.
  - X2D: M4 0.5225-0.6100 in 2.5 MeV x M(ee) 2-62 MeV in 2 MeV.
  - BEFORE: M4 0.30-1.10 in 5 MeV; M(pipi), M(ee), OA, vz, MM4, pT4 and y4 in 8 M4 regions; an M4 x M(pipi) map.
- **Data:** N, N++, N--, CB; folds; periods I/IIA/IIB/III; isBest == 1.
- **Sim:** C2-C6 (C6 per unit BR), G_eta, G_pi0D, G_conv, G_vec, G_oth (the last three also split by omega/rho0 pions), T3, T4. Event-group variances, halves A/B, 100 Poisson event-group replicas (nominal only).
- **Variants.**
  - Weights: V2; bump 1+aY or phase space; conversion F_V off; Gounaris-Sakurai F_V; weights off.
  - Mode B; w_db; ghost (a).
  - Scales: s = 0.99521, s +- 0.0009; r = 0.89, 1.10.
  - Selections: pre; pre + dz_pi < 40; PEAK; vz > -150 and > -170.
  - Data only: cor_ok fallback.
- **Acceptance.**
  1. sim pre T1+T2 = 701.1 nb, and C2+C3 in the sim bump window = 9.689 nb (s = 1);
  2. data S0 peak window: N = 55947, CB/N 0.059;
  3. data "none" (pre without vz > -160) bump-window W = 148091;
  4. the sum of the cells equals the 1D spectrum to 1e-9.

### WP1 "Before" (absolute, per-channel cross sections, no renormalisation)
- **Figures.**
  - B1: M4 0.30-1.10, log. N black, CB red, N - CB blue; stacked sim groups; variant (a) dashed.
  - B2: 0.40-0.70, linear.
  - B3: data/sim map in M4 x M(pipi), with E(m, 0), the peak band and the CR boxes.
  - B4: small multiples of the other variables.
  - Ratio panels: band = MC (+) 6.9 %, and a separate one-sided +7.3 % live-time band.
- **Table B1** for S0, pre and pre isBest == 1 (data only).
  - data/sim in <0.455, B, 0.5225-0.5336, P, 0.5542-0.565, 0.565-0.70, 0.70-0.90, 0.90-1.10.
  - Jump J = R(0.560-0.595)/R(0.4875-0.5225).
  - Shape chi2/ndf (0.40-0.70) after one free scale.
  - sigma per group; the K0s peak, data/sim.
- **Table A2.** Corrections in plain words with their procedure: CB, 64/63 and L, scales, V3, beam energy, pion-vertex efficiency data/sim (0.917/0.941 peak, 0.945/0.950 tagged), cor_ok, ghosts.
- **Acceptance.** pre reproduces rw_compare_enriched.txt x 1.0161 within 0.5 % (bump 0.958, peak 1.233); J 1.44-1.47.

### WP2 Normalisation items (parallel)
- **Sector deficit.**
  - sector = (floor(phi/60) + 5) mod 6; verify on the phi gaps.
  - Rows per nb^-1, per day and per sector used, in S0.
  - C_sec = 1 + f_s (1/eps_s - 1) L_034-044/L: eps_s = the early/late rate ratio of rows using the deficient sector; f_s = the signal share with a track there (sim).
  - Variant: days 045-067 only, with their own L.
- **Empty target.**
  - Look for empty-target runs in the run lists and CSVs.
  - If there are none: fit the vz profile of the ECAL-tagged eta and of the bump-window N - CB with the sim LH2 shape plus two window components. f_win is subtracted with a 100 % error.
  - If the fit fails: f_win = 0.7-2.1 % (elastic w6 = 0.70 % x an A-dependence of 1-3).
- **Stability.** a_eta and DR per period, with the background fixed (after WP5); chi2 p > 0.01.
- **Output.** Table N1, Fig. N1.

### WP3 Control regions, stage A (after WP0)
- **Stage A.**
  - F10 restricted to the CR-U and CR-K bins, background free.
  - Eta components at a_eta(0) = 0.4; rerun at the WP5 value; a_eta = 1 shown.
  - The D11 rule is applied here and frozen.
- **Acceptance.**
  - chi2/ndf <= 1.5.
  - A1: sim C2-C5 < 1 % in every CR bin (G_eta share printed).
  - A2: CR-K null test, an edge step at E(0.4116, 0) in 0.43-0.46 and 0.46-0.55: |h| < 2 sigma (earlier 0.42-0.46: 36.5 +- 19.3). The tagged m5c peak in CR-K is consistent with 0.
  - A3: eta injected x {0, 0.4, 1, 2} moves c_j and b by < 0.3 sigma.
  - A4: per-slice tilts consistent with a common b (p > 0.01); the CR-U split at M(pipi) 0.42 is consistent.
  - A5: a K0s veto (0.485-0.507) changes the parameters by < 1 sigma.
- **Output.**
  - Table C1: candidate CRs with contamination and verdict, including the rejected 0.60-0.80, below the bump, vz < -160, like-sign and ECAL tag.
  - Table C2: parameters, covariance, nested delta chi2.
  - Figs. C1 (N - CB map with boxes and E for q = 0 and 60 MeV), C2 (data/renormalised sim per slice), C3 ("after" spectrum, eta at SMASH scale).
  - Table C3: counted "after" display in the Table B1 regions.

### WP4 Stress closure (after WP0; gate for showing WP5/WP6 numbers)
- **Toys.** Parametric, at L_data, 500 per scenario.
  - Truth: half B (smoothed background, raw signal) plus a CB-like component (data CB shape, with N++ and N-- Poisson-fluctuated).
  - Fit templates: half A.
- **Scenarios.**
  - (a) background x 1.35 with tilt -0.15/0.1, a_eta 0.4, rho_d 1, kappa 1, g 1;
  - (b) rho_d 0.5/1.5;
  - (c) kappa 0.8/1.2;
  - (d) g 0.5/1.5;
  - (e) omega x 0.3 and the rest x 1.6;
  - (f) s + 0.0009, r = 0.89;
  - (g) a_eta = 0;
  - (h) a_eta = 1.
- **Estimators.** F-main, F-free, F-gold, B4, the peak side bands, X17.
- **Pass.**
  - Bias below max(0.5 sigma_stat, 5 %) for a_eta and Y_bump; max(0.5 sigma_stat, 10 %) for Y_peak; 20 % for rho_d and BR_G.
  - 68 % coverage 0.60-0.76.
  - Null scenario: |mean pull| < 0.5.
- **Afterwards.** |bias| of (a) is the S10 entry; an estimator that fails is demoted. Output: Table V1, Fig. V1.

### WP5 Main extraction, stage B (after WP3; sealed until WP4 passes)
- **Fits.**
  - F-main.
  - F-free.
  - F-gold on F5g: alpha_2 and alpha_5 free, background free, C3 and C4 fixed from F-main (they are < 1 % there).
  - Per-slice a_eta,j; folds; pre; data isBest == 1.
- **Derived.**
  - BR_G = 2.68e-4 alpha_5/alpha_2 (F-gold).
  - BR_A = 2.68e-4 rho_d.
  - BR_C = 2.68e-4 alpha_5 (alpha_3/alpha_2)/alpha_4.
  - Closure R_34 = alpha_4/alpha_3 = 1 expected (the PDG ratio 0.184 is in the sim; the material cancels at first order).
  - kappa vs kappa_B.
- **Acceptance.**
  - chi2/ndf <= 1.2 (above that, errors x sqrt(chi2/ndf), stated).
  - |rho(a_eta, b)| < 0.7.
  - Per-slice a_eta p > 0.01.
  - No pull run longer than 4 bins in 0.52-0.56.
  - A free s within 2 sigma of 0.99615.
  - Folds p > 0.05.
- **Output.**
  - Tables E1 (parameters and correlations) and E2 (yields: stat, MC, CB, syst).
  - Figs. E1 (10 panels + pulls), E2 (N - CB - background with the eta templates per slice), E3 (profile scans).

### WP6 Side bands and cross-checks (after WP3; "other eta" terms from WP5)
- **Bump B4.**

      Y = (N - CB)_B - (N - CB)_U [B_A(B)/B_A(U)] - a_eta (C4 + C5 + g G_eta)_B

  - B_A = the stage-A background (tilt from CR-U + CR-K, never from CR-K alone); U = 0.565-0.700.
- **Peak.** The same with the peak side bands and B_A + bump tail; the peak leakage into the side bands is iterated out with C4 + C5.
- **Closure.**
  - 13 displaced windows 0.575-0.695, each removed from CR-U before the refit: sigma_tr = sigma_stat sqrt(max(0, rms^2 - 1)).
  - CR-K and CR-L prediction pulls.
- **Cross-checks.**
  - P_main on S0 (4388 +- 307 +- 1078).
  - f_eta(edge) = E_data/E_sim, with the same edge-step fit on sim T1+T2.
  - The ECAL tag against the emulated sim tag.
- **Output.** Table E3, Fig. E4.

### WP7 Systematics (after WP5 and WP6)
- Refit every Sect. 5 variant from the WP0 histograms.
- Table S1; Fig. S1 (ribbons).

### WP8 Physics quantities (after WP7; absolute values also after WP2)
- **Table P1.** Efficiencies per class, window and M(ee) class (V3 and unweighted); pion-vertex data/sim.
- **eps2 QED check at M(ee) >= 20 MeV.**
  - Generator C2 fraction against Kroll-Wada with the pi0 form factor (HGeant/PLUTO closure 0.891).
  - Tagged data shares 0.077 (20-30 MeV) and 0.144 (> 30 MeV) against sim.
- **Table P2.**
  - BR_G, BR_A, BR_C against PDG 2.68 +- 0.11 e-4.
  - DR; kappa twice; R_34.
  - sigma(pp -> eta X) = a_eta x 1.0665 mb x C_sec x (1 - f_win), against SMASH.
  - sigma x BR; visible cross sections.
  - (pT4, y4)-reweighting variant (3x3 bins, to the side-band-subtracted peak).
- **Literature.** Published inclusive pp -> eta X near p_lab 5.4 GeV/c, references supplied by the analyst. Otherwise "no independent comparison". Never invent numbers.
- **Fig. P1** summary.

### WP9 X17 (X-a, X-b after WP0; X-c after WP5; X-d after the freeze)
- **X-a Model** on X2D:

      mu_ij = k CB_ij + n_j B_e(j)(M4_i) + P_j g_e(j)(M4_i) + a_eta [C2 + C3 + g G_eta]_ij

  - B_E1 and B_E2 are the WP5 background shapes (pol2 as a variant); n_j is free per M(ee) bin.
  - P_j = N5 D_j + N4 C_j + r N5 (eps_X/eps_5) X_j(m_X). N5 and N4 are free over the full M(ee) range.
  - D = C5 (V3); C = C4, replica-averaged.
  - X = the C6 line (17.25 / 4.04 MeV, double Gaussian). Other masses come from C5 rows in true-M(ee) slices of +-0.5 MeV, reweighted to X17 kinematics, validated on C6.
  - N5-bias nuisance = the spread of N5 over the background-shape variants.
- **X-b.** The 1500 existing X17 toys plus a CB-like component: slope 1 +- 0.05, 90 % coverage >= 0.88. If asymptotic and toy CLs differ by more than 10 %, the toys are primary.
- **X-c.**
  - Fit with r = 0 on the M(ee) side bands 2-13 and 27-62 MeV.
  - Asimov limit and 2000 background-only toys (median, +-1/2 sigma).
  - Write x17_frozen.txt.
- **X-d.**
  - Observed CLs (q~_mu, r >= 0); Feldman-Cousins cross-check.
  - Scan 10-40 MeV (1 MeV; 0.5 MeV in 15-20); quote 16.7, 16.85 and 17.0 MeV.
  - Local p0; global p from toys if Z > 3.
  - Counting in the frozen window.
- **Nuisances.**
  - Direct shape (V3/V2/V4, Gounaris-Sakurai F_V); C replicas.
  - M(ee) scale +-2 % and resolution +-10 %, checked on the E1 conversion peak and the pi0 Dalitz shape.
  - Peak width.
  - eps_X/eps_5: M(pipi) phase space vs direct-like; decay angle isotropic vs 1 + cos^2.
  - k; N5 bias; BR_dir +-4 %; a_eta and g.
- **Table D1.**
  1. primary: r x 2.68e-4;
  2. N_X/(a_eta s1), s1 = 2.127e7 rows per unit BR at SMASH, rescaled to 6079.2 and x 63/64 for raw rows;
  3. relative to Y_peak.
  - Also limits on N_X and on sigma x BR x BR [nb]; the sim expectation (1.9e-5 / 5.7e-6) for comparison only.
- **Figs.** X1 (M(ee) spectrum with D, C and X17 at the limit), X2 (limit vs m_X with bands), X3 (p0), X4 (validation).

### WP10 Reports
- Skeleton first.
- Every number is traced to a tables_all/*.txt file and its macro.
- Captions give the selection, units, L and the fallbacks used.

## 4. Dependencies

- WP0 -> {WP1, WP2, WP3, WP4, X-a, X-b}.
- WP3 -> WP5. WP4 PASS unseals WP5 and WP6.
- WP5 -> {WP6, X-c, WP2 stability}. WP6 -> WP7 -> WP8. X-c -> freeze -> X-d. All -> WP10.
- **Day 1:** WP0, WP1, stage A, B4 + peak side bands, F-gold, the X17 Asimov limit.

## 5. Systematics (Table S1; columns: Y_bump, Y_peak, Y5, Y4, a_eta, BR_G, BR_A, kappa, DR, sigma_eta, UL_X17)

| # | family | variations (nominal first) |
|---|---|---|
| S1 | CB | k_E2 1 +- 0.05 / fixed 1 / 1.032; k_E1 free / 1; like-sign sum; 20 MeV cells; per-cell geometric mean; E1 dropped |
| S2 | CR definition | CR-U lower edge 0.560/0.565/0.575, upper 0.65/0.70/0.75; CR-K 0.425/0.43/0.44; K0s veto |
| S3 | background model | c_j(1+bx); + quadratic; per-slice tilt; exp(bx); common c + c1 y; omega split, range 0.85; kappa_B = kappa |
| S4 | background template | exp(pol3) / exp(pol4) / raw BB-lite; ghost (a); mode B; w_db |
| S5 | range, binning | 0.40-0.70 / 0.42-0.65 / 0.45-0.70 / 0.40-0.75; 2.5 MeV; slice edges +-5 MeV; M(ee) split 8/12 MeV |
| S6 | mass transfer | s +- 0.0009; 0.99521; r 0.89/1.10; M(pipi) scale +-0.0007; s free |
| S7 | physics weights | V3 / V2 / V4; KLOE / 1+aY / phase space; F_V on/off, GS (weights off shown only) |
| S8 | G_eta | width 0.30 -> 0.15/0.50; g free |
| S9 | signal templates | raw / factorised f(ee) g(slice) h(M4 given slice) |
| S10 | closure | WP4 (a) bias; sigma_tr (peak); max(pred - obs, sigma_pred) in CR-K/CR-L (bump) |
| S11 | selection | S0 / pre; pivtx / dz_pi < 40; vz -150/-170 (efficiency-corrected quantities only) |
| S12 | efficiency model | pion vertex +-2.5 %; eps2 QED at M(ee) >= 20; (pT4, y4) reweighting; X17 kinematics |
| MC | template statistics | 100 replicas, always added, quoted separately |
| N | normalisation (absolute only) | L 6.9 %; live time +7.3/-0 % (separate); C_sec; f_win; PDG BRs |
| - | folds, periods, isBest == 1, cor_ok fallback | consistency tests only (a fold counts if p < 0.05) |

**Combining.**
- delta^2 = max(0, Delta^2 - |sigma_v^2 - sigma_0^2|).
- Per family: the maximum over discrete alternatives; (|D+| + |D-|)/2 for a +-1 sigma pair.
- Families add in quadrature; MC and N are listed apart.
- In S11 a variation counts only if |Delta| > 2 sqrt(|sigma_v^2 - sigma_0^2|).
- Ratios are recomputed inside each variant, so correlations are automatic.

## 6. Fall-backs (the reports list which ones fired)

| trigger | response |
|---|---|
| A1 or A2 fails | CR-K at 0.44, CR-U from 0.575. If it still fails: drop CR-K, take b from CR-U only, label the bump "extrapolated" and add the CR-L term to S10. |
| Stage A chi2/ndf > 1.5 | Add terms (D11 order, then the omega split). If still > 1.5: errors x sqrt(chi2/ndf), "background-model limited". |
| Stage B chi2/ndf > 1.5, or pull runs in 0.52-0.56 | S3/S4 variants. If unresolved, the side-band estimators become the main result. |
| \|rho(a_eta, b)\| > 0.7, or a_eta ~ 0 with a significant E_step | Fix b to stage A (+-1 sigma propagated). If still degenerate: f_eta(edge) only, no absolute bump yield. |
| Per-slice a_eta p < 0.01 | Per-slice a_eta,j becomes nominal; the difference is a systematic. |
| corr(rho_d, kappa) > 0.9, or the direct-decay bias > 20 % | Quote Y_peak only; BR from F-gold. |
| R_34 off by > 2 sigma | Withdraw BR_C. |
| k_E1 outside 0.8-3, or like-sign null pull > 2 | Drop E1; kappa = 1 +- 0.2. |
| An estimator fails WP4 | Cross-check only. If F-main fails, B4 + the peak side bands are the main result. |
| X17: Pearson chi2/ndf of n_j > 1.5; slope off by > 5 %; N5 error > 50 % | Scale the errors; correct the slope and add a systematic; normalisation (2) becomes primary. |
| C_sec or f_win not determinable | Absolute values PRELIMINARY with a stated extra uncertainty; ratios unaffected. |
| MC error > data statistical error | "Simulation-limited"; request a larger std production. |

## 7. report_all.tex (style of the sim report: noteboxes, booktabs, colours)

0. Status notebox; Table of decisions D1-D20 (open for review); notation longtable.
1. Scope and inputs: production, L, sim sets, Table A2.
2. How SMASH describes the data before renormalisation: Figs B1-B4, Table B1.
3. Control regions and renormalisation: Tables C1-C3, Figs C1-C3, tests A1-A5.
4. Validation on simulation: Table V1, Fig V1.
5. Extraction of the bump and the peak: the model in plain words, Tables E1-E2, Figs E1-E3, F-free, F-gold, folds, isBest == 1.
6. Side bands and cross-checks: Table E3, Fig E4.
7. Systematic uncertainties: Table S1, Fig S1, the combining rule.
8. Physics results: Tables P1, P2, N1; Figs N1, P1.
9. X17 limit: model, blinding statement, Figs X1-X4, Table D1.
10. Conclusions and open points.

Appendices: fallback log; files and macros; normalisation details.

## 8. report_all_short.tex (8 pages; captions refer to Table nota)

1. p1: summary, condensed decisions, notation.
2. p2: Fig B1 with a zoom inset, short Table B1.
3. p3: Figs C1 + C2, short Table C2, the CR verdict.
4. p4: condensed Fig E1 (S2, S4, S5 x E2, peak zoom), Fig E2.
5. p5: Tables E2 + E3, short Table S1.
6. p6: Table P2.
7. p7: Figs X1 + X2, Table D1.
8. p8: open points and the fallbacks that fired.

## 9. Judge fixes and where each is handled

| fix | where |
|---|---|
| Empty target, sector 2 | WP2 |
| Live-time sign | D17 |
| Edge depends on M(ee) | E(m, q) |
| G_eta factor | D10, S8 |
| CB in E1 | D5, S1 |
| omega | S3, WP4 (e) |
| ECAL tag | D14 |
| Blinding, frozen Asimov | D19, X-c |
| X17 window transfer | Sect. 2 |
| No sim background per M(ee) bin, N5 bias | WP9 |
| Class names | header |
| PIVTX next to pre | D1 |
| cor_ok fallback not combined | Sect. 5 |
| Literature | WP8 |
| CR-K at 0.43 | Sect. 2 |
| No double use | D18 |
| MC statistics | D13 |
| Source of the tilt | WP6 |
| Distorted closure | WP4 |
| eps2 QED | WP8 |
| K0s | A5, S2 |
| MM4 | D4 |
| Ghosts | D6 |
| 6079.2 and 64/63 | D7 |
