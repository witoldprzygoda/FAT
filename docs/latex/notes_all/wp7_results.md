# WP7 systematics (prefix wp7) -- results note

Date 2026-11-14. Branch pp45_pippimepem (worktree /home/przygoda/HADES/FAT/FAT). No git command that changes state was run.
ScanEngine.h, scan_*.C, hists_v1.root, the sealed file and every file of WP0-WP6 and X17 were not touched; no other package was
re-run. FAT_simwt was not used (only the WP0 copies). Nothing was deleted.

**Machine-readable output:** `/hdd3/przygoda/HADES/FAT_out/all_report/wp7/syst.txt` (md5 4ec7f6d7...).
- `syst` lines: value, stat, the 12 family deltas, total syst, three alternative totals, MC, N up/down, fold term (absolute units).
- `rel` lines: the same in per cent.
- `item` lines: every variant of every quantity, with Delta, the counted delta and the reason when it is not counted.
- `period`, `fold`, `best`, `corfb` lines: the consistency tests.
- `x17` line: the UL_X17 column.

**Headline (S0, isBest >= 0, L = 6079.2 nb^-1; yields = raw rows x 64/63, not efficiency corrected).**

| quantity (estimator) | value | stat | syst | MC | N |
|---|---|---|---|---|---|
| Y_peak (peak side bands, MAIN) | 5088 | 432 (8.5 %) | 582 (11.4 %) | 245 (4.8 %) | -- |
| sigma(pp -> eta X) from the side bands (MAIN) | 0.437 mb | 0.037 | 0.070 (16.0 %) | 0.026 | +0.035/-0.037; live time +0.032/-0 |
| Y_bump (B4, all M(ee), MAIN, EXTRAPOLATED) | 16389 | 1129 | 25766 (157 %) | 4110 | -- |
| DR (side bands / B4, MAIN) | 1.40 | 0.15 | 0.93 (66 %) | 0.53 | 0.06 |
| a_eta (F-main, CROSS-CHECK) | 0.382 | 0.017 | 0.205 (54 %) | 0.040 | 0.029 |
| BR_G (F-gold, CROSS-CHECK) | 2.65e-4 | 0.37e-4 | 3.36e-4 (127 %) | 0.66e-4 | 0.09e-4 |
| BR_A (F-main, CROSS-CHECK) | 3.52e-4 | 0.31e-4 | 3.21e-4 (91 %) | 0.61e-4 | 0.11e-4 |
| UL_X17, 16.7-17.0 MeV/c2 (x17d) | 3.06e-5 | in the limit | in the frozen limit | 11 % | -- |

- **The peak is measured.** Y_peak and sigma_eta carry 11-16 % systematics. The largest sources are:
  - the simulated peak resolution (S6, r = 0.89/1.10: 5.7 % on Y_peak, 9.0 % on sigma_eta);
  - the side-band width (S2, 5.6 %);
  - the background model (S3, 5.7 %);
  - the direct-decay weights (S7 V2: 7.0 % on sigma_eta);
  - the displaced-window closure (S10, 4.7 %).
- **The absolute bump yield from B4 is not determined.** Its systematics exceed 100 %:
  - a quadratic background term or the WP3 data-only D11 model gives +139 % (S3);
  - the CR-U upper edge at 0.65 gives +41 % (S2);
  - raw or exp(pol4) templates give -44 % (S4);
  - the CR-L closure term adds 38 % (S10).
  This confirms the WP3/WP6 label EXTRAPOLATED. DR inherits it (66 %).
- **The F-main and F-gold decompositions (a_eta, Y4, Y5, kappa, BR_A, BR_G) are template- and model-limited at 50-130 %.**
  They stay cross-checks (WP4 FAIL).
- **The per-period stability test (WP2) PASSES** for every quantity and variant (p = 0.04-0.99).

---------------------------------------------------------------------------------------------------
## 1. What was done, step by step (plain words)

1. **Rows for every selection** (`wp7_rows.C`, 30 s).
   - Every data row (e+e-, e+e+, e-e-) with isBest >= 0, vz > -500, PT3 and vz > -170 mm is written to a small tree in file
     order, so that the event groups stay together.
   - Each row keeps the columns that the selection variants need: vz, dz_pi, MM4, the pion-vertex flag, cor_ok, isBest,
     period and fold.
   - For cor_ok == 0 rows the tree holds the reconstructed-pi- (_fb) kinematics. They are used only in the cor_ok fallback test.
   - Check: S0 has N = 55 947 in P and 111 212 in B, the same as hists_v1.root.
2. **A generic version of the validated code** (`wp7_lib.h`). Nothing in it is a new method; it generalises WP0/WP4/WP5/WP6
   code so that any variant can be refitted.
   - **Data.** For any layout, selection, subset and CB rule, the counts N, the CB (2 sqrt(N++ N--) in 10 MeV cells, spread rule)
     and the 200 CB replicas are rebuilt from the rows. The replicas use the WP0 generator and seed. The event-group covariance
     is rebuilt as well. This is wp6_ana.C LoadDataW/CBWin, generalised.
   - **Simulation.** Every hists_v1 directory is read in the F10/F5g layout. The G3 layout (S0__nom) is summed into any other
     layout: other range, binning, slice edges, M(ee) split, K0S veto, or simulated M(pipi) scale.
   - **Smoothing.** wp4::SmoothFit is used for F-main and F-gold. The WP3 smoothing of wp6_ana.C (copied, generalised to
     exp(pol n)) is used for stage A.
   - **Fitter.** A copy of wp5::Fit (Fisher scoring, BB-lite, sandwich errors) with optional model terms: exp(b x), an
     omega/rho0 scale, a common c (1 + c1 y), kappa_B = kappa, BB-lite on raw background templates.
   - **Estimators.** B4 and the peak side bands exactly as wp6_ana.C. Their stat error is N + CB + stage A + F-main, as in
     WP6.
3. **The nominal chain reproduces WP5 and WP6 to all printed digits** (wp7_syst.log, CHECK block):
   - a_eta 0.38228397, kappa, rho_d, g, Y4, Y5, DR_F, BR_A, BR_G 2.6546e-4;
   - B4 7822.17 (E2) / 16389.3 (all), side bands 5087.97, a_P 0.418249, DR 1.39743;
   - stage A as WP6;
   - stat errors to 1e-3 (the CB replica covariance uses 1/(n - 1));
   - MC errors to all digits: a_eta 0.039629, B4 1919.6 / 4109.8, side bands 245.12, BR_G 6.6354e-5.
   The G3-built templates give the same results to 1e-7.
4. **Every Sect. 5 variant is refitted** (`wp7_syst.C`, 65 configurations, 30 s). The full chain is rerun each time: F-main ->
   F-gold -> stage A at the F-main eta values -> B4 and side bands.
   - **S2 (control regions):** only stage A and the estimators change; F-main is held at its nominal fit.
   - **S3, S4, S8, S9 (background model and template, G_eta, signal templates):** the change is applied to F-main and to stage A
     together.
   - **Ratios** are recomputed inside each variant (DR, BR_G, BR_A).
5. **MC noise of every variant difference.**
   - For every variant that uses the S0__nom (or pre__nom) simulation sample in the standard layout, the whole chain is rerun
     with each of the 100 WP0 bootstrap replicas. This follows the WP5/WP6 procedure: the WP4 smoothing is redone, the WP3
     smoothing starts from the nominal parameters, and the fits start at the nominal values.
   - Replica r of the variant and replica r of the nominal use the same simulated events. So the per-replica difference
     measures the template noise of the variant shift. Its robust spread is MC(Delta), half of the central 68 %.
   - For variants with other simulation samples (V2, ghostA, s/r shifts, ...) or G3 layouts, no replicas exist and
     MC(Delta) = 0 (conservative).
6. **Combination** (`wp7_combine.py`). For each quantity and variant:

       Delta = variant - reference
       delta^2 = max(0, Delta^2 - |stat_v^2 - stat_0^2| - MC(Delta)^2)

   - The reference is the nominal, or the G3-built nominal for G3 variants.
   - A +-1 sigma pair gives (|D+| + |D-|)/2.
   - A family takes the maximum over its items. Families add in quadrature.
   - S11 counts only for efficiency-corrected quantities, and only if |Delta| > 2 sqrt(|stat_v^2 - stat_0^2| + MC(Delta)^2).
   - MC and N are listed apart.
   - Rule R1 (WP7 decision, Sect. 5): a variant rejected by the data is listed but not counted.
7. **External entries.**
   - S10 (closure): WP4 RF scenario (a) biases and the WP6 closure terms (sigma_tr = 239 rows; CR-L 37.5 %).
   - S2: the WP6 side-band width variants.
   - S12: pion-vertex efficiency.
   - N: L 6.9 %, PDG BRs, C_sec (1 - f_eta) of AMENDMENT 1.
   - UL_X17: the x17d variant table.
8. **Consistency tests:** folds, isBest == 1, cor_ok fallback, and the per-period stability test with the WP2 inputs.
   - "F-main type": the nominal background shape is kept, with one free scale and kappa and g fixed.
   - "Side-band type": the nominal stage-A shape and F-main eta terms are kept; only the period data change.
   - The data errors of each period come from its own rows (event groups and CB replicas).

---------------------------------------------------------------------------------------------------
## 2. Numbers

### 2.1 Table S1 (per cent of the value; `docs/latex/tables_all/wp7_tableS1.{txt,tex}`)

| | Y_bump | Y_peak | Y5 | Y4 | a_eta | BR_G | BR_A | kappa | DR | sigma_eta | UL_X17 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| value | 16389 | 5088 | 3198 | 1847 | 0.3823 | 2.65e-4 | 3.52e-4 | 0.832 | 1.397 | 0.437 mb | 3.06e-5 |
| stat | 6.9 | 8.5 | 9.1 | 14.9 | 4.4 | 13.9 | 8.7 | 17.4 | 10.9 | 8.5 | (in UL) |
| S1 CB | 5.8 | 1.1 | 2.6 | 11.0 | 2.7 | 0.0 | 1.8 | 7.2 | 3.1 | 1.1 | (1.8) |
| S2 CR definition | 41.3 | 5.6 | 0 | 0 | 0 | 0 | 0 | 0 | 7.7 | 5.6 | -- |
| S3 background model | 138.5 | 5.7 | 11.1 | 38.3 | 37.4 | 18.3 | 28.9 | 49.6 | 46.2 | 5.7 | (25.6) |
| S4 background template | 43.9 | 2.0 | 86.0 | 35.1 | 34.1 | 116.5 | 83.2 | 2.8 | 13.3 | 5.6 | (13.7) |
| S5 range, binning | 19.6 | 0.4 | 16.2 | 13.1 | 12.0 | 24.3 | 13.2 | 11.7 | 22.7 | 0.4 | -- |
| S6 mass transfer | 6.6 | 5.7 | 5.6 | 4.4 | 4.3 | 7.0 | 5.3 | 3.3 | 4.3 | 9.0 | (5.8) |
| S7 physics weights | 0.7 | 2.6 | 4.1 | 9.4 | 1.5 | 12.5 | 13.6 | 5.6 | 6.4 | 7.0 | (18.4) |
| S8 G_eta | 4.5 | 0.3 | 0 | 3.4 | 5.7 | 28.8 | 5.4 | 9.1 | 2.9 | 0.3 | -- |
| S9 signal templates | 0 | 0 | 5.9 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | -- |
| S10 closure | 37.6 | 4.7 | 1.6 | 12.1 | 10.1 | 14.8 | 9.7 | 1.8 | 37.9 | 4.7 | -- |
| S11 selection | n/a | n/a | n/a | n/a | 1.7 | 16.4 | 4.5 | 0.0 | 0.0 | 0.0 | -- |
| S12 efficiency model | 0 | 0 | 0 | 0 | 2.5 | 2.0 | 2.0 | 2.0 | 2.0 | 2.5 | (19.9) |
| **total syst** | **157.2** | **11.4** | **88.7** | **57.1** | **53.7** | **126.6** | **91.1** | **52.8** | **66.3** | **16.0** | in limit |
| MC | 25.1 | 4.8 | 17.8 | 26.9 | 10.4 | 25.0 | 17.2 | 32.1 | 38.2 | 5.9 | 11.1 |
| N | -- | -- | -- | -- | 7.6 | 3.2 | 3.2 | -- | 4.0 | +7.9/-8.4 | -- |

Notes on the table:
- **UL_X17:** values in brackets are information from x17d (observed UL under the X-c model variants). These choices are
  inside the frozen limit through the N5 bias (15 %) and sX (0.247), so they are not added. S12 = the e+ angle 1 + cos^2
  (3.67e-5), quoted separately as in x17d.
- **sigma_eta live time:** +7.3/-0 % (+0.032 mb), separate.
- **F-main quantities:** the stat errors are not scaled. WP5's chi2 rule would multiply them by 1.61 (GLS-rep) or 1.38 (v5
  hybrid GLS).
- **Alternative totals** (syst.txt):

  | quantity | Barlow with stat only | without rule R1 | sub-sources of S5/S6/S7 in quadrature |
  |---|---|---|---|
  | Y_bump | 171 | 157.5 | 157.4 |
  | Y_peak | 14.7 | 15.2 | 11.5 |
  | sigma_eta | 17.7 | 18.6 | 17.0 |
  | DR | 118 | 67.7 | 66.7 |
  | a_eta | 55.6 | 53.8 | 54.4 |
  | BR_G | 132 | 144 | 126.6 |

  Values in per cent.
- **Further quantities** (syst.txt), value +- stat +- syst +- MC:
  - Y_bump(E2) = 7822 +- 456 +- 12870 (165 %) +- 1920;
  - a_P = 0.418 +- 0.035 +- 0.067 (16.0 %) +- 0.025, N 7.3 %;
  - a_B = 0.312 +- 0.018 +- 0.523 (168 %) +- 0.076;
  - DR_F = 1.169 +- 0.065 +- 0.660 (56 %) +- 0.134;
  - sigma_eta from F-main = 0.400 +- 0.018 +- 0.215 (54 %) +- 0.041 mb;
  - g = 1.72 +- 0.16 +- 2.07 +- 0.37; rho_d = 1.32 +- 0.12 +- 1.20 +- 0.23.

### 2.2 The dominant variant per family (Table S1 quantities)

| quantity | dominant items |
|---|---|
| Y_peak | S3 D11data + omega (+6.6 % shift, 5.7 % counted); S6 r = 0.89/1.10 (-7.6/+4.7 %); S2 wide/narrow side bands (-6.7/-7.0 %, WP6); S7 conversion F_V off (+3.3 %); S4 GLS-hybrid F-main (+2.1 %); S10 4.7 % |
| sigma_eta | S6 r pair (-10.8/+8.4 %); S7 V2 (-7.7 %; the simulated C5 changes); S4 w_db (-6.3 %); S3 D11data + omega (+6.6 %); S2 side-band width; S10 4.7 %; S12 2.5 % |
| Y_bump | S3 D11data (+153 %), quad (+144 %); S4 raw (-49 %), exp(pol4) (-35 %); S2 CR-U to 0.65 (+47 %); S5 M(ee) split 8/12 MeV (-16/-24 %); S10 CR-L 37.5 % |
| DR | S3 D11data/quad (DR 0.58/0.60); S10 37.9 %; S5 M(ee) split; S4 w_db (13 %) |
| a_eta | S3 quad (+39 %); S4 GLS-hybrid (+34 %); S5 range 0.42-0.65 (+15 %); S10 10.1 % |
| BR_G | S4 raw BB-lite (+121 %); S8 g free (+32 %); S5 range 0.45-0.70 (-28 %); S3 quad (-25 %); S11 dz_pi < 40 (+17 %); S7 V2 (-14 %) |
| BR_A | S4 raw (+87 %); S3 quad (-31 %); S7 V2 (-14 %); S5 range 0.42-0.65 (-17 %) |
| Y5 / Y4 | S4 raw (+89 % / +6 %) and GLS-hybrid (+14 % / +38 %); S3 kappa_B = kappa (Y4 +50 %); S5 range; S1 per-cell geometric-mean CB (Y4 -12 %) |
| kappa | S3 kappa_B = kappa (1.42); S5 slice edges; S8; S1 geometric-mean CB (0.76) |

### 2.3 Variants rejected by the data (rule R1; listed, not counted)

- **Same bins:** a variant is rejected if (chi2_v - chi2_0)/(chi2_0/ndf_0) - (ndf_0 - ndf_v) > 25.

  | variant | fit | value | consequence |
  |---|---|---|---|
  | S3 common c (1 + c1 y) | F-main | 67 | all quantities |
  | S3 common c (1 + c1 y) | stage A | 144 | |
  | S7 bump phase space instead of KLOE | F-main | 38 | all quantities |
  | S7 bump phase space instead of KLOE | F-gold | 28 | |
  | S9 factorised signal templates | F-gold | 71 | BR_G only |
  | S8 g = 1 +- 0.15 | F-gold | 28 | BR_G only |

- **Other bins:** a variant is rejected if (chi2_v/ndf_v)/(chi2_0/ndf_0) > 1 + 5 sqrt(2/ndf_v).

  | variant | fit | ratio (limit) | Y_peak shift |
  |---|---|---|---|
  | S3 omega split with range 0.85 | F-main | 1.63 (1.32) | +10.9 % |
  | S3 omega split with range 0.85 | stage A | 1.55 (1.43) | |
  | S5 F-main range 0.40-0.75 | F-main | 1.37 (1.36), marginal | +2.0 % |

- Without R1 the totals are 15.2 % (Y_peak), 18.6 % (sigma_eta), 133 % (Y5) and 144 % (BR_G) (Sect. 2.1).

### 2.4 Per-period stability test (`tables_all/wp7_periods.{txt,tex}`; WP2 Table N1 inputs)

- Templates x L_p/L, with L_p = 2809.453 / 372.610 / 1078.022 / 1819.110 nb^-1 (sum 6079.2).
- K_p = C_sec(p) x 0.97, with C_sec(I) = 1.0241 +- 0.0077 and 1 for the other periods.

| quantity | factor | I | IIA | IIB | III | chi2/ndf | p | verdict |
|---|---|---|---|---|---|---|---|---|
| a_eta (F-main type) | raw | 0.367(12) | 0.423(32) | 0.383(19) | 0.368(15) | 3.17/3 | 0.37 | PASS |
| a_eta | K (f_eta) | 0.365 | 0.411 | 0.372 | 0.357 | 2.65/3 | 0.45 | PASS |
| a_eta | old K (f_win; IIA-IIB correlated) | 0.303 | 0.350 | 0.317 | 0.322 | 2.15/3 | 0.54 | PASS |
| a_P (side bands) | raw | 0.417(52) | 0.455(145) | 0.426(86) | 0.408(64) | 0.10/3 | 0.99 | PASS |
| a_P | K | | | | | 0.12/3 | 0.99 | PASS |
| a_B (B4) | raw | 0.318(16) | 0.410(44) | 0.302(26) | 0.289(19) | 6.76/3 | 0.080 | PASS |
| a_B | K | | | | | 7.28/3 | 0.064 | PASS |
| a_B | old K | | | | | 4.94/3 | 0.18 | PASS |
| DR_F (F-main type) | -- | 1.29(10) | 0.68(21) | 1.02(14) | 1.21(11) | 8.38/3 | 0.039 | PASS |
| DR (side bands) | -- | 1.37(18) | 1.16(39) | 1.47(32) | 1.47(25) | 0.55/3 | 0.91 | PASS |

- The F-main-type values agree with WP5 (0.36684 / 0.42339 / 0.38319 / 0.36782). The small differences come from the analysed
  L shares (WP2) instead of kLper, and from the period's own CB replicas instead of the scaled "all" covariance.
- Period IIA (6 % of L) is high in a_B (+2.2 sigma) and low in DR_F (-2.3 sigma). Neither reaches p < 0.01.

### 2.5 Other consistency tests (`tables_all/wp7_tests.txt`)

- **Folds** (event-group parity; yields x inverse row fraction): no quantity has p < 0.05.
  - Largest: Y_peak_F z = -1.55;
  - Y_peak / a_P / sigma_eta z = -1.43 (fold0 4479 +- 611, fold1 5713 +- 610);
  - a_eta z = -1.37.
  - No fold term is added.
- **isBest == 1 / isBest >= 0** (data only, the same isBest >= 0 templates, D2; the stat error is the independent-sample bound):

  | quantity | ratio |
  |---|---|
  | Y_peak | 0.915 +- 0.109 |
  | Y_bump | 0.885 +- 0.087 |
  | a_eta | 0.867 +- 0.055 |
  | DR | 1.03 +- 0.16 |
  | BR_G | 1.05 +- 0.21 |
  | BR_A | 1.04 +- 0.13 |

  This agrees with WP5 (0.888 / 0.868 / 0.867 / 1.023).
- **cor_ok fallback** (nested sample; z = Delta/sqrt|stat_v^2 - stat_0^2|):
  - the ratios are unchanged: DR z = -0.06, BR_G +0.18, BR_A +0.44;
  - the yields and absolute scales rise by about 2.2 % (Y_peak +111 rows, z = +2.3; a_eta z = +2.2; a_B z = +3.2).
  - This is expected: the fallback adds about 1.3 % e+e- rows (P: 56 672 against 55 947), and the simulation has no cor_ok == 0
    rows. Not combined (spec).

---------------------------------------------------------------------------------------------------
## 3. Figures and tables (draft captions)

**Fig. S1** `docs/latex/figs_all/wp7_S1.pdf` (.png). *Systematic uncertainties of the report_all quantities.*
- Data S0, isBest >= 0, L = 6079.2 nb^-1. Simulation SMASH mode A, w x V3, M4' = 0.99615 M4.
- (a) Relative uncertainty in per cent, per family S1-S12 (Sect. 5 of the spec), with the total systematic, statistical and
  template-statistics (MC) errors. Colour = log10 of the per cent value; "n/a": S11 is not applied to raw yields.
- (b) Nested ribbons around value/nominal = 1: stat (dark), stat (+) syst, (+) MC, (+) normalisation N (absolute quantities).
- MAIN: Y_bump (B4, extrapolated), Y_peak (peak side bands), DR, sigma_eta (from the side bands).
- CROSS-CHECK (WP4 FAIL): Y5, Y4, a_eta, BR_A, kappa (F-main); BR_G (F-gold).
- Family deltas after the Barlow subtraction of the statistical and template-statistical noise of each variant difference;
  family = maximum over its variants; variants rejected by the data (rule R1) are not counted.

**Table S1** `docs/latex/tables_all/wp7_tableS1.tex` (+ `.txt` with absolute values).
- Draft caption: "Systematic uncertainties in per cent of the value; S0, isBest >= 0, L = 6079.2 nb^-1; yields in raw rows x
  64/63.
  - Columns: Y_bump = B4 (M(ee) >= 10 MeV extrapolated to all M(ee), EXTRAPOLATED), Y_peak = peak side bands, Y5/Y4 = F-main C5/C4
    in P, a_eta, BR_A = 2.68e-4 rho_d and kappa from F-main, BR_G from F-gold, DR = (Y_peak/Y_bump)_data/(...)_sim from the side
    bands, sigma_eta = a_P x 1.0665 mb x C_sec (1 - f_eta), UL_X17 from the frozen X17 fit.
  - The rule: delta^2 = max(0, Delta^2 - |stat_v^2 - stat_0^2| - MC(Delta)^2); family = maximum over the variants (pairs: mean of
    |D+| and |D-|); families in quadrature; MC and N separate; S11 only for efficiency-corrected quantities.
  - UL_X17 entries in brackets are information: these choices are inside the frozen limit."

**Table** `tables_all/wp7_variants.txt`: every variant of every quantity (value, Delta, stat and MC of both, MC(Delta), counted delta,
R1 verdict).

**Table** `tables_all/wp7_periods.{txt,tex}`: per-period stability test.
- Draft caption: "a_eta (F-main, background shape fixed), a_P (peak side bands), a_B (B4), DR_F and DR per period with the
  background fixed to the all-data fit; data errors from the period's event groups and CB replicas; raw, x K_p = C_sec(p)
  (1 - f_eta), and x the old K_p = C_sec (1 - f_win,p) with IIA-IIB correlated; chi2 against the weighted mean."

**Table** `tables_all/wp7_tests.txt`: folds, isBest == 1, cor_ok fallback.

---------------------------------------------------------------------------------------------------
## 4. Acceptance criteria and fall-backs

| item | value | verdict |
|---|---|---|
| refit every Sect. 5 variant from the WP0 histograms | 59 variants + 6 tests, S1-S9 and S11 refitted; S10, S12, MC, N from WP4/WP6/WP1/WP2 | DONE (exceptions below) |
| nominal reproduces WP5/WP6 (values, errors, MC) | to all printed digits; stat to 1e-3 | PASS |
| G3-built nominal = F10 nominal | 1e-7 | PASS |
| Table S1 with the 11 columns | yes | DONE |
| Fig. S1 (ribbons) | yes | DONE |
| per-period stability (WP2): chi2 p > 0.01 | 0.039-0.99 | PASS |
| folds: a fold counts if p < 0.05 | min p 0.12 | none counts |

**Fall-backs:**
- **"MC error > data statistical error -> simulation-limited": fired.**
  - It fires for Y_bump (MC 25 % vs stat 6.9 %), DR (38 % vs 11 %), a_eta (10 % vs 4.4 %), Y4, Y5, BR_A, BR_G and kappa.
  - It does not fire for Y_peak and sigma_eta (MC 4.8-5.9 % vs stat 8.5 %).
- **Inherited and respected:**
  - E1 dropped (k_E1 outside 0.8-3);
  - CR-K dropped, bump EXTRAPOLATED, CR-L term in S10;
  - F-main/F-free/F-gold cross-checks, B4 and the side bands are the main result (WP4);
  - stage-B chi2/ndf > 1.5.
- **Not fired:** C_sec / f_win not determinable (WP2 determined both).

**Spec items not done, or done in a different way (flagged):**
1. **S7 "V4":** there is no V4 directory in hists_v1.root. Only V2 (and GS) exist, and S7 uses them.
2. **S12 eps2 QED at M(ee) >= 20 MeV and the (pT4, y4) reweighting:** these are WP8 items and need a new simulation pass. They are
   left to WP8 and marked "pending" in syst.txt and Table S1.
3. **S6 "M(pipi) scale +-0.0007":** no simulated directory exists. It is done with G3 by moving the simulated slice edges by
   -+0.07 % and sharing the 5 MeV M(pipi) bins linearly. Effect at most 0.4 % on the main quantities.
4. **S6 "s free":** linear interpolation of the sUp/sDn shifts to the WP5 free-s values (0.42 x sUp for 0.99653; 1.28 x sDn for
   0.99500); the larger one is taken.
5. **S1:**
   - "k_E1 free / 1" and "E1 dropped" do not apply after the E1 fall-back. The spec configuration E1 + E2 with k_E1 free is
     listed as information.
   - "1.032" is implemented as CB x 1.032 everywhere with k_E2 fixed at 1.
   - "like-sign sum" is read as N++ + N-- instead of 2 sqrt(N++ N--).
6. **S2:**
   - CR-K variants (0.43, 0.425, 0.44) are information only, because CR-K was rejected by A1/A2. 0.425/0.44 move the S4/S5 edge
     in the whole fit (G3), not only in CR-K.
   - The side-band width variants (narrow 5.6 MeV, wide 11.6 MeV) are taken from WP6, not refitted, with a stat-only Barlow
     subtraction.
7. **S3 "omega split, range 0.85"** is one variant: omega/rho0 scale free, F-main 0.40-0.85, CR-U 0.575-0.85. It is rejected by R1.
   The omega split on the standard range is counted.
8. **S4 "GLS hybrid":** the template-noise-aware estimator of the v5 verification (minimised hybrid GLS: a 0.5126, kappa 0.856,
   rho 1.118, g 0.465, kappa_B 2.871, from v5_quad2.log). It is added because v5 asked for an S3/S4 entry. F-main is not
   refitted by WP7; F-gold, stage A and the estimators are.
9. **S5:** the F-gold range follows the F-main range. 0.40-0.75 and 2.5 MeV bins are built from G3.
10. **sigma_eta** is computed from the peak side bands (a_P), not from F-main a_eta: F-main failed WP4 and the side bands are the
    main result. The F-main version (0.400 mb) is in syst.txt as `sigma_eta_F`. **For WP8 to decide.**
11. **Barlow subtraction with MC(Delta)** (WP7 decision, see Sect. 5). The literal stat-only version is given next to it.

---------------------------------------------------------------------------------------------------
## 5. Caveats

1. **Template noise is subtracted from every variant shift (decision).**
   - The analysis is simulation-limited: MC/stat 2-4 for B4 and F-main.
   - Without the subtraction, a variant that changes how the noisy templates enter counts the template noise twice: once in
     the MC row and once in the variant shift. Examples: raw templates instead of smoothed, exp(pol4) instead of exp(pol3).
   - Example: raw BB-lite moves Y_peak by -451 rows, but the replica spread of that difference is 1580 rows. So it is not
     counted.
   - The stat-only (literal) version raises Y_peak from 11.4 % to 14.7 %, sigma_eta to 17.7 % and DR to 118 %.
2. **Rule R1 (decision).** Variants that the data reject (Sect. 2.3) are not counted. The threshold is 25 in chi2 after scaling
   with the nominal chi2/ndf, i.e. more than 5 sigma. The S5 range 0.40-0.75 is rejected only marginally (1.37 vs 1.36).
   Totals without R1 are given.
3. **The family maximum is the spec's rule.** S5, S6 and S7 contain independent sub-sources. Their quadrature changes the
   totals by less than 1 point (`syst_quadSubsrc`).
4. **The bump systematics come from the extrapolation of the background from CR-U into B.**
   - S3 quadratic, S2 CR-U edges, S4 smoothing order, S10 CR-L, S5 M(ee) split and S11 pre (a_B -34 %) all move B4 by 20-150 %.
   - The absolute bump yield (and DR) should not be quoted as a measurement.
   - f_eta(edge) = 0.392 +- 0.049 (WP6) is the robust bump-scale check. It has no Table S1 column, and WP7 did not vary it.
5. **The F-main decomposition** is model-, template- and estimator-limited.
   - raw BB-lite gives rho_d 2.45 and BR_G 5.9e-4; GLS hybrid gives a_eta 0.51; kappa_B = kappa gives Y4 +50 %.
   - g ranges from 0.47 to 2.6 over the variants.
   - The WP5 values remain cross-checks; their errors would grow by x 1.38-1.61 with the chi2 rule of WP5.
6. **MC(Delta) is available only for variants on the S0__nom/pre__nom samples in the standard layout.** It is 0 elsewhere
   (conservative): simulation variants, s/r shifts, G3 layouts.
7. **S11 uses the pre__nom replicas for pre** (same event-group weights as S0). dz40/vz150/vz170 have no replicas.
   - S11 does not count for DR: the per-replica spread of the pre - S0 difference is 1.17, because the B4 denominator is
     close to 0 in some pre replicas.
   - It counts for a_B (pre: 0.205 vs 0.312, 32 %) and for BR_G (dz40: +17 %).
8. **G3-built templates** have variances summed per G3 cell (a slight under-estimate of the event-group variance). The G3
   nominal equals the F10 nominal to 1e-7, so the effect on the BB-lite weights is negligible.
9. **N row:**
   - PDG BRs: 3.2 % for the bump channels (eta -> pi+pi-pi0 22.92 +- 0.28 %, pi0 Dalitz 1.174 +- 0.035 %); 2.3 % for the peak mix
     (eta -> pi+pi-gamma 4.22 +- 0.08 %, eta -> pi+pi-e+e- 2.68 +- 0.11 e-4, C4 : C5 = 0.48 : 0.52 in P).
   - L 6.9 %; C_sec (1 - f_eta) +3.1/-4.1 % (AMENDMENT 1). f_win does not enter any eta quantity.
   - The PDG values are typed from memory of the PDG tables used in the simulation study (kBRdir, kBRppg in
     sim_eta_physics.h agree). WP8 should confirm them.
10. **S12 pion vertex** (WP1 Table A2, earlier studies, not recomputed): 2.5 % on absolute scales. For the ratios DR, BR_G, BR_A
    and kappa: 2.0 %, the peak/bump difference of the data/sim pion-vertex efficiency (0.917/0.941 vs 0.945/0.950).
11. **Period test:** the stage-A and F-main errors are common to all periods and are left out of the chi2; the template noise
    is common as well. IIA and IIB share the empty-target day 047 (old-K variant correlated).
12. **Blinding:** only F10, F5g, G3 coarse classes (8/10/12/20 MeV) and window counts were used. No X2D cell or 2 MeV M(ee) shape
    was formed.
13. **Machine:** one process at a time, at most 1.6 GB, about 30 s per full run; CPU package 42-43 C.

---------------------------------------------------------------------------------------------------
## 6. For the other packages

- **WP8:** syst.txt `syst <q>` lines (absolute) and `rel <q>` lines (per cent).
  - For sigma(pp -> eta X), take `sigma_eta` (side bands; 0.437 mb) or `sigma_eta_F` (F-main; 0.400 mb). The value already
    includes C_sec (1 - f_eta) = 0.981.
  - S12 eps2 QED and the (pT4, y4) reweighting are still to be added by WP8.
  - The live time (+7.3/-0 %) is quoted separately.
- **WP10:** Table S1 = `tables_all/wp7_tableS1.tex`, Fig. S1 = `figs_all/wp7_S1.pdf`, periods = `tables_all/wp7_periods.tex`. The
  combination rule, R1 and the MC(Delta) subtraction must be stated in the text (Sect. 1, items 5-6).

---------------------------------------------------------------------------------------------------
## 7. Exact commands (from /home/przygoda/HADES/FAT/FAT)

    B='gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_wp7", kTRUE)'
    root -l -b -q -e "$B" 'plots/all_report/wp7_rows.C+O'       # 30 s, 3 GB: wp7/wp7_rows.root, wp7_rows.txt
    root -l -b -q -e "$B" 'plots/all_report/wp7_syst.C+O(0)'    # 30 s, 1.6 GB: wp7/syst_raw.txt, syst_rep.txt, periods_raw.txt, wp7_syst.log
                                                                #   argument 1 = nominal check only (writes *_check files)
    python3 -I plots/all_report/wp7_combine.py                  # 5 s: wp7/syst.txt; tables_all/wp7_{tableS1,variants,periods,tests}.*; figs_all/wp7_S1.{pdf,png}

**Code** (`plots/all_report/`):
- `wp7_rows.C`: rows.
- `wp7_lib.h`: the generic layer (it includes wp5_lib.h and through it wp4_lib.h).
- `wp7_syst.C`: the configurations, the replica chains and the periods.
- `wp7_combine.py`: the combination, the tables and the figure.
  - The typed external inputs are listed with their source at the top: WP4 summary (RF scenario (a) biases), WP6 results.txt
    (side-band widths, sigma_tr, CR-L), WP1 (pion-vertex efficiencies), WP2 (C_sec, f_win per period, L per period), x17d
    (variant table, MC, angle), v5 (GLS-hybrid parameters).

**Outputs:**
- `/hdd3/przygoda/HADES/FAT_out/all_report/wp7/`: syst.txt, syst_raw.txt, syst_rep.txt, periods_raw.txt, wp7_rows.root,
  wp7_syst.log, wp7_syst_stdout.txt, wp7_combine_stdout.txt.
- `docs/latex/tables_all/wp7_tableS1.{txt,tex}`, `wp7_variants.txt`, `wp7_periods.{txt,tex}`, `wp7_tests.txt`.
- `docs/latex/figs_all/wp7_S1.{pdf,png}`.
