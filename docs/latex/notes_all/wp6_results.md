# WP6 side bands and cross-checks (prefix wp6) -- results note

Date 2026-11-14. Branch pp45_pippimepem (worktree /home/przygoda/HADES/FAT/FAT). No git command that changes state was run.
ScanEngine.h, scan_*.C, hists_v1.root and every WP0/WP3/WP4/WP5 file were not touched, and no other package was re-run.
WP3's smoothing code was COPIED into my macro (attributed). FAT_simwt was used only through the WP0 copies.

**Machine-readable output:** `/hdd3/przygoda/HADES/FAT_out/all_report/wp6/results.txt`
- `est` lines: value, stat (N, CB, stage A, F-main), MC, total, MC shift, bias-corrected, mc68.
- Also: `stageA`, `b4parts`, `psbparts`, `disp`, `crpull`, `xcheck`, `simsmash` lines.
- md5 d8dd8948...; the output is deterministic.

**Headline (S0, isBest >= 0, L = 6079.2 nb^-1; raw rows x 64/63; no efficiency correction).**

| quantity | value | stat | MC | status |
|---|---|---|---|---|
| Y_peak, side bands, all M(ee) | **5088** | 432 | 245 | WP4 PASS; + sigma_tr 239 (S10) |
| Y_bump, B4, M(ee) >= 10 MeV | **7822** | 456 | 1920 | WP4 FAIL marginal (bias-corrected 8638 passes); EXTRAPOLATED; S10 (CR-L) 37.5 % |
| Y_bump, B4, all M(ee) | 16389 | 1129 | 4110 | E2 / f_E2 (simulation), extrapolation +1196/-1356 |
| f_eta(edge) = E_data/E_sim | 0.392 | 0.020 | 0.039 | model 0.021; cross-check |
| f_eta(tag) = T_data/T_sim | 0.454 | 0.014 | 0.040 | emulated simulated tag, UNVERIFIED |
| P_main on S0 | 4388 | 307 | -- | spread 1078; reproduced exactly |

- Every eta-scale estimate agrees with F-main's a_eta = 0.382 +- 0.043 (Fig. E4 d):
  - peak / SMASH: 0.418 +- 0.043;
  - bump / SMASH (E2): 0.312 +- 0.079;
  - f_eta(edge): 0.392 +- 0.049;
  - f_eta(tag): 0.454 +- 0.043;
  - P_main: 0.403 +- 0.103.
- **Main caveat: the B4 bump background is EXTRAPOLATED from CR-U and not under control.** Two tests show it:
  - With the linear tilt the prediction lies above the data in every region below the peak: CR-L pull -1.88, CR-K -1.98,
    0.400-0.455 -2.78.
  - With an added quadratic term the prediction lies below the data (+1.62, -1.01, +1.74), and Y_bump(E2) becomes 19526
    (+150 %).

---------------------------------------------------------------------------------------------------
## 1. What was done, step by step (plain words)

1. **Row-level pass** (`wp6_fill.C`: data rows, then simulated rows; 40 s, about 10 GB).
   - **Data S0 rows.** All S0 rows (e+e-, e+e+, e-e-) are written to a small tree (`wp6_rows.root`), so that windows of any
     width can be counted with exact edges.
   - **Simulated rows.** Weighting:
     - SMASH, owned rows (mode A), w x V3 (code 108);
     - T1 + T2;
     - expected raw data rows = sum w x 27.0718.

     Mass transfer to the data scale: M4' = 0.99615 M4, M(pipi)' = 0.9952 M(pipi).
   - **Simulated window sums.** For each (slice S1-S5 x M(ee) class E1/E2) and each of C2, C3, C4, C5, G_eta, G_D, G_C and T12,
     the sums are filled in exact windows:
     - B, P, the side bands (nominal, narrow 5.6 MeV, wide 11.6 MeV);
     - U = 0.575-0.700 and 0.565-0.700, CR-L, the CR-K range, 0.400-0.455;
     - the 13 displaced windows (P_k, SBlo_k, SBhi_k).

     Each sum has its event-group variance and 100 replicas.
   - **Replicas.** I use WP0's generator and seed, so replica r here is the same as replica r of hists_v1. Checked: replica 0
     of C2 in B (E2) is 19706.7461 in both.
2. **Edge step, data** (procedure of scan_bump.C, frozen).
   - Slices M(pipi) 0.30-0.34, 0.34-0.37, 0.37-0.39 and 0.39-0.40.
   - Fit of N - CB: pol2 [+ eta Gaussian] + h x step, with the step at the fixed edge 0.99277 M4max(M(pipi)/0.9970) and a
     smearing of 5.4 MeV.
   - E = the sum of the four h.
3. **Edge step, simulation: the SAME fit** on the simulated T1 + T2 spectra.
   - f_eta(edge) = E_data / E_sim.
   - MC error: the same fit on the 100 replicas.
   - Model spread: the 18 scan_bump variations (order 2/3; half range 35/45/55 MeV; smearing 4.4/5.4/6.4 MeV) applied jointly
     to data and simulation. A varied slice fit with p <= 0.05 keeps its nominal h (scan_bump rule).
   - The same was done for E1 and E2 separately. There the class CB is the slice CB x the class share of the like-sign rows in
     the 50 MeV band.
4. **ECAL tag** (frozen definition).
   - Selection: good photon, M(e+e-gamma) 0.119-0.173, M(ee) < 0.10.
   - The data T is the m5c Gaussian content (0.5463 / 6.2 MeV, pol2, 0.49-0.62).
   - The emulated simulated tag uses the same row selection. The same fit is applied to the simulated T1 + T2 tagged rows,
     with free-shape and truth versions next to it.
   - Additional measurement: the simulated pi0 position in M(e+e-gamma). The tag window was then also re-centred on the
     simulation (variant).
5. **P_main on S0**: ScanEngine FitPeakMain + its 48-fit model spread. The same frozen fit is also applied to the simulated
   T1 + T2 (r = 1 and r = 0.89).
6. **F-main of WP5 refitted** with wp5_lib.h (F10 E2, kappa = 1 +- 0.2), together with its 100 replica fits.
   - It reproduces WP5 results.txt to 4e-8 (all 12 parameters).
   - Its (a_eta, kappa, rho_d, g) are the "other eta" inputs.
7. **Stage A re-run at the WP5 eta values** (spec WP3: "rerun at the WP5 value").
   - The frozen WP3 configuration: M1 = c_j (1 + b x), E2 only, CR-U 0.575-0.700, no CR-K, k_E2 = 1 +- 0.05.
   - The eta terms are fixed to F-main (a = 0.382, kappa = 0.832, rho = 1.315, g = 1.717).
   - The smoothing is WP3's own, copied; the fitter is wp5_lib.
   - B_A(window) = c_j x the integral of (1 + b x)[G_D + kappa_B G_C] over the exact window.
8. **B4** (WP4 implementation, on E2):

       Y = (64/63) {(N - CB)_B - [(N - CB)_U - eta_U] B_A(B)/B_A(U) - [a kappa C4 + a rho C5 + a g G_eta]_B}

   - eta_U = all eta classes in U at the F-main values.
   - All M(ee): Y_E2 / f_E2, with f_E2 = (C2 + kappa C3)_E2 / (C2 + kappa C3)_all in B = 0.4773 (simulation).
9. **Peak side bands** (all M(ee), one R'; WP4 implementation):

       Y = (64/63) [(N - CB)_P - R' (N - CB)_SB] / (1 - R' f)

   - R' = [B_A + a (C2 + kappa C3) + a g G_eta](P) / [same](SB).
   - f = (C4 + C5)(SB) / (C4 + C5)(P) = 0.1544 is the leakage of the peak into the side bands.
   - "Iterated out": Y_{n+1} = (N - CB)_P - R' [(N - CB)_SB - f Y_n] converges to the closed form; the iterations are listed.
10. **Errors of every estimator.**
    - stat = the sum in quadrature of:
      - data N (event-group covariance of all window cells);
      - CB (200 like-sign replicas, with the WP0 spread rule rebuilt for every window);
      - stage A (its sandwich covariance, numerical gradient);
      - F-main (its sandwich covariance; total derivative, i.e. stage A is refitted at a_eta, kappa, rho, g +- h).
    - MC = RMS over 100 replica chains: F-main replica r -> stage A with replica r (smoothing redone) -> the estimator with the
      replica-r window sums.
    - MC shift = replica mean - nominal; bias-corrected = 2 v - mean (WP4 "RB").
11. **Closure.**
    - **Displaced windows.**
      - 13 window centres, 0.575-0.695.
      - For each, the window and its side bands are removed from CR-U, stage A is refitted, and the side-band estimator is
        applied there (the leakage f is the one of the nominal peak).
      - pull = (Y - Y_exp)/sigma_stat; Y_exp = the simulated C4 + C5 content passed through the same estimator.
      - sigma_tr = sigma_stat sqrt(max(0, rms^2 - 1)).
    - **CR-L, CR-K and 0.400-0.455 (all E2).**
      - Prediction = k_E2 CB + B_A + eta(F-main).
      - Errors: CB, stage A, F-main and MC.
      - S10 term = max(|N - pred|, sigma_pred). It is transferred to the bump as the same relative error on the B4 background in
        B.
12. **Variants:** side-band definitions, B4 definitions, stage-A variants (quadratic term, WP4 smoothing, CR-U 0.565, WP3
    frozen stage A), and B4 per slice.

---------------------------------------------------------------------------------------------------
## 2. Numbers

All numbers: S0, isBest >= 0, raw data rows; yields x 64/63; L = 6079.2 nb^-1; simulation SMASH x PDG BR, w x V3.

### 2.1 Input checks (all exact)
- **Data window sums vs hists_v1 `win`** (B, P, SBlo, SBhi, CR-L, CR-U 0.565):
  - N is identical; CB agrees to 1e-10.
  - CB replica RMS: P 56.7 (hists_v1 56.6), B 91.1 (90.8).
- **E2 sums in B and U575 vs F10:** N and CB are identical.
- **Simulated window sums vs hists_v1:** C2 in B (E2) 19750.3420 in both; C5 in P 6260.7553 in both; the replicas are the same.
- **F-main refit:** a_eta 0.38228 +- 0.01692 (stat) +- 0.03963 (MC); max |difference to WP5| 3.9e-8.
- **Data E_step (S0):** 1141.9 +- 59.2 +- 236 (data team "pre + vertex": 1142 +- 59 +- 236).
- **Tag T (S0, all M4):** 3090.3 +- 94.3 (WP3: 3090 +- 94).
- **P_main:** 4387.8 +- 306.6 +- 1078.2 (spec: 4388 +- 307 +- 1078).

### 2.2 Stage A at the WP5 eta values (WP3 configuration and smoothing)

| | c_1 | c_2 | c_3 | c_4 | c_5 | b | kappa_B | k_E2 |
|---|---|---|---|---|---|---|---|---|
| this (eta = F-main) | 0.9816 | 0.9844 | 0.9505 | 1.0473 | 0.9145 | -0.0639 | 2.932 | 0.955 |
| stat | 0.041 | 0.042 | 0.042 | 0.041 | 0.033 | 0.0058 | 0.269 | 0.048 |
| this, eta of WP3 (a 0.4, g 1) | 0.9917 | 1.0050 | 0.9746 | 1.0582 | 0.9113 | -0.0661 | 3.024 | 0.954 |
| WP3 frozen | 0.9878 | 1.0012 | 0.9751 | 1.0552 | 0.9091 | -0.0683 | 3.063 | 0.954 |

- **Reproduction of WP3.** With WP3's eta, the frozen WP3 values are reproduced within 0.1 sigma(stat). The small rest comes
  from the weights: WP3 includes the eta-template variance in them.
- **WP3's closure predictions are reproduced:** CR-L 8060.2 (WP3 8078.0), CR-K 440.2 (440.4).
- **chi2:**
  - GLS with the template MC covariance: 172.8/118 = 1.46 (WP3 at a = 0.4: 1.10).
  - Data-only Pearson: 384/118.
- **WP4-smoothing variant** (as used in F-main): GLS 1.56, data-only Pearson 439/118.
- **Background in B, E2:**
  - B_A(B) = 37763, B_A(U) = 296527 (exact integral; the bin model gives 37763 / 296530).
  - B4 background in B = (N - CB - eta)_U x R = 37568.

### 2.3 Side-band estimators (Table E3)

| quantity | value | stat (N / CB / stage A / F-main) | MC | MC shift | bias-corr. |
|---|---|---|---|---|---|
| Y_bump B4, E2 | 7822 | 456 (250 / 87 / 368 / 47) | 1920 | -816 | 8638 |
| Y_bump B4, all M(ee) | 16389 | 1129 | 4110 | -1996 | 18385 |
| Y_peak side bands, all | 5088 | 432 (426 / 69 / 6 / 9) | 245 | -36 | 5124 |
| Y_peak, E1 (< 10 MeV) | 2380 | 291 | 164 | | |
| Y_peak, E2 (>= 10 MeV) | 2780 | 326 | 179 | | |
| peak share M(ee) >= 10 MeV | 0.546 | 0.046 | 0.027 | | |
| sigma_vis bump (all M(ee)) [nb] | 2.696 | 0.186 | 0.676 | | 3.024 |
| sigma_vis peak [nb] | 0.837 | 0.071 | 0.040 | | |
| bump / SMASH (E2, kappa = 1) | 0.312 | 0.018 | 0.076 | | 0.344 |
| peak / SMASH (C4 + C5) | 0.418 | 0.035 | 0.025 | | |
| DR (side bands; bump extrapolated) | 1.40 | 0.15 | 0.53 | +0.31 | 1.09 |
| DR (both in E2) | 1.73 | 0.23 | 0.61 | | 1.38 |

**B4 ingredients (E2):**
- N_B 53241, CB_B 5667.8; N_U 342098, CB_U 32696.3.
- eta_U 14399 (G_eta 14375 at a g = 0.656).
- R = B_A(B)/B_A(U) = 0.12735.
- Other eta in B: 2305 (C4 12, C5 115, G_eta 2178).
- Extrapolation to all M(ee): kappa = 0.832 +- 0.304 gives f_E2 0.4773 (0.4448 / 0.5203), i.e. Y_bump(all) +1196 / -1356.

**Side-band ingredients:**
- N_P 55947, CB_P 3279.0; N_SB 43656, CB_SB 2695.0.
- B_A(P) 59210, B_A(SB) 49859; bump tail + G_eta: 2435 in P, 2121 in SB.
- R' = 1.1859 (background alone 1.1875; flat background, the width ratio: 1.1930).
- Iteration (raw rows): 4091.5, 4840.6, 4977.7, 5002.8, 5007.4, 5008.3; the limit is 5008.5, x 64/63 = 5088.0.

### 2.4 Variants (value +- stat +- MC)

| variant | value | change |
|---|---|---|
| **side bands** | | |
| per class (E1 + E2 separately) | 5159 +- 438 +- 250 | +1.4 % |
| only the background normalised by the side bands | 5112 +- 432 +- 273 | +0.5 % |
| leakage shape kappa C4 + rho C5 | 5100 +- 433 +- 240 | +0.2 % |
| no leakage correction (f = 0) | 4156 +- 353 +- 192 | -18 % |
| narrow side bands (5.6 MeV) | 4733 +- 510 +- 301 | -7.0 % |
| wide side bands (11.6 MeV) | 4745 +- 388 +- 225 | -6.7 % |
| flat background (R' = width ratio) | 4733 +- 433 | -7.0 % |
| stage A + quadratic term | 5263 +- 431 +- 293 | +3.4 % |
| stage A with the WP4 smoothing | 5068 +- 432 +- 244 | -0.4 % |
| stage A on 0.565-0.700 | 5093 +- 432 +- 245 | +0.1 % |
| WP3 frozen stage A | 5073 +- 432 | -0.3 % |
| **B4, E2** | | |
| without the eta_U subtraction (spec literal) | 5959 +- 482 +- 2079 | -24 % |
| CB x k_E2 of stage A | 7891 +- 447 +- 1931 | +0.9 % |
| per-slice transfer | 7810 +- 461 +- 1919 | -0.2 % |
| normalised on U 0.565-0.700 (same stage A) | 7867 +- 449 +- 1895 | +0.6 % |
| stage A and U on 0.565-0.700 | 8213 +- 424 +- 1717 | +5.0 % |
| stage A with the WP4 smoothing | 8182 +- 460 +- 2056 | +4.6 % |
| WP3 frozen stage A (eta 0.4, g 1; MC from its parameter covariance) | 7389 +- 452 +- 1597 | -5.5 % |
| **stage A + quadratic term** (q = -0.173 +- 0.022 stat) | **19526 +- 1434 +- 3919** | **+150 %** |
| E1 + E2 direct (E1 shape predicted with the E2 parameters) | 12256 +- 884 +- 3738 | test only |
| per slice S1 / S2 / S3 / S4 | 1567 / 3426 / 2055 / 878 (MC 1147 / 673 / 424 / 387) | |
| per slice S5 (no bump, null test) | **-117 +- 21 +- 52** | -2.1 sigma |

### 2.5 Closure

**Displaced windows** (centres 0.575, 0.585, ..., 0.695; all M(ee); window removed from CR-U, stage A refitted):
- pulls: +0.97, +0.12, -0.49, -0.38, +0.84, -0.78, +1.94, +1.10, -1.17, +2.29, +0.64, +0.08, +1.53;
- rms 1.143 (mean +0.51, std. dev. 1.06). Without the refit: 1.181. E2-only estimator: 1.572.
- **sigma_tr = 431.6 x sqrt(1.143^2 - 1) = 239 rows = 4.7 % of Y_peak -> S10 (peak).**
- The 13 windows overlap (10 MeV steps, 37.8 MeV wide), so the pulls are correlated.
- At 0.575 the lower side band holds the tail of the real peak; Y_exp = -279 accounts for it.

**Prediction pulls (E2; prediction = k CB + B_A + eta(F-main)):**

| region | N | prediction (stat, MC) | pull | S10 term | as % of B_A | -> on Y_bump(E2) |
|---|---|---|---|---|---|---|
| CR-L 0.5225-0.5300 | 7519 | 8018 (59, 244) | -1.88 | 499 | 7.7 % | 2931 rows = **37.5 %** |
| CR-K (S5, 0.400-0.5225) | 331 | 441 (11, 51) | -1.98 | 110 | 37.8 % | 14436 = 185 % |
| 0.400-0.455 (below B) | 13786 | 16086 (173, 801) | -2.78 | 2300 | 18.4 % | 7025 = 90 % |

Pulls of the same regions with other stage-A variants:

| stage A | CR-L | CR-K | 0.400-0.455 |
|---|---|---|---|
| + quadratic | +1.62 | -1.01 | +1.74 |
| WP4 smoothing | -1.56 | -3.08 | -2.54 |

### 2.6 Cross-checks (D14: not combined)

**Edge step (S0), h in rows / 2.5 MeV:**

| M(pipi) | data h | sim h (T1 + T2) | data/sim |
|---|---|---|---|
| 0.30-0.34 | 287 +- 37 | 841 | 0.342 |
| 0.34-0.37 | 389 +- 33 | 847 | 0.459 |
| 0.37-0.39 | 331 +- 24 | 882 | 0.376 |
| 0.39-0.40 | 135 +- 22 | 346 | 0.389 |

- **f_eta(edge) = 1141.9 / 2915.7 = 0.392 +- 0.020 (stat) +- 0.039 (MC) +- 0.021 (model)**.
  - By class: E1 0.367 +- 0.024 +- 0.045 +- 0.076; E2 0.431 +- 0.036 +- 0.073 +- 0.073.
- **The simulated step comes from the bump classes:** E(T12 - C2 - C3) = 40.6 (1.4 % of E_sim), and C2 + C3 alone give 3010.
- **Control slices:**
  - 0.40-0.41: data 18.7 +- 14.7, sim 116 (0.16);
  - 0.42-0.46 (edge by hand): data 15.4 +- 15.0, sim -14.5.

**ECAL tag (S0):**
- **Data T:**
  - all M4: 3090.3 +- 94.3 (free shape: 0.5463, 5.8 MeV, T 2978);
  - bump region (M(pipi) < 0.43, M4 < 0.5225): 3252.6 +- 91.9;
  - CR-K: -2.1 +- 4.0.
- **Simulated T1 + T2, emulated tag:**
  - all M4: 6804 (frozen) +- 603 (MC), 6867 free (0.5479, 6.0 MeV);
  - truth C2 + C3: 9020 (8833 in 0.49-0.62);
  - bump region: 6995.
- **f_eta(tag) = 0.454 +- 0.014 (stat) +- 0.040 (MC).** Variants:
  - free shapes 0.434;
  - bump region 0.465;
  - tag window re-centred on the simulation (0.113-0.145): 0.413 (frozen fit), 0.372 (free).
- **The simulated pi0 in M(e+e-gamma)** (C2 + C3 with a good photon): 0.1287, sigma 8.0 MeV. Data: 0.1459, 13.25 MeV. The photon
  energy scale is -4.7 % in the simulation and +8 % in the data.
- **E2 share of the tagged eta:** data 0.458 +- 0.011; simulation 0.426 (fit), 0.439 (truth). This supports the f_E2
  extrapolation of the bump within 4-8 %.

**P_main:**
- Data: 4388 +- 307 +- 1078. In P x 64/63: 4202 (+- 294 stat, +- 1073 with the spread).
- The same frozen fit on the simulated T1 + T2: 10881 (r = 1), 12098 (r = 0.89). S_data/S_sim = 0.403 / 0.363.
- The truth C4 + C5 in P is 11856.

**F-main (WP5, cross-check):** Y_bump 19426 +- 2478, Y_peak 5045 +- 804, a_eta 0.382 +- 0.043.

---------------------------------------------------------------------------------------------------
## 3. Figures and tables (draft captions)

- **Fig. E4** `docs/latex/figs_all/wp6_E4.pdf` (.png). *Side-band estimators and cross-checks.* Data S0, isBest >= 0,
  L = 6079.2 nb^-1, raw rows.
  - (a) B4 in M(ee) >= 10 MeV/c^2.
    - Blue: data N - CB (CB = 2 sqrt(N++ N--), spread rule), per MeV/c^2, summed over the five M(pi+pi-) slices.
    - Grey: the stage-A background (smoothed simulation G_D + kappa_B G_C x c_j (1 + b x), fitted in CR-U 0.575-0.700),
      normalised to the data in U (green bar) and extrapolated into the bump window B (blue bar).
    - Violet: the other eta classes (C4, C5, G_eta at the F-main values).
    - Orange: the bump C2 + kappa C3 at the F-main values, for its shape.
    - Y_bump(E2) = data minus grey minus violet in B, x 64/63.
  - (b) Peak side bands, all M(ee): N - CB per MeV/c^2 in 2.5 MeV bins.
    - Green solid: background (B_A + bump tail + G_eta) scaled to the side bands (green bands) after the removal of the peak
      leakage.
    - Red dashed: plus the simulated C4 + C5 shape at the extracted yield.
    - Red band: window P.
  - (c) Closure: pulls of the side-band estimator at 13 displaced windows.
    - Full circles: the window removed from CR-U and stage A refitted. Open circles: nominal stage A. Squares:
      M(ee) >= 10 MeV only.
    - Grey band: +-1. sigma_tr = 239 rows. The CR-L and CR-K pulls are quoted in the panel.
  - (d) Eta scale estimates relative to SMASH sigma x PDG BR: F-main a_eta (WP5), bump B4 / simulation (E2), peak side bands /
    simulation, f_eta(edge), f_eta(tag) (emulated simulated tag, UNVERIFIED), P_main. Thick: stat; thin: total.
- **Fig.** `wp6_edge.pdf`. *Edge step of the bump in data and simulation.*
  - Top row: data S0 N - CB in M(pi+pi-) 0.30-0.34, 0.34-0.37, 0.37-0.39, 0.39-0.40, rows per 2.5 MeV/c^2.
  - Bottom row: simulated T1 + T2 at SMASH x PDG BR, expected raw rows (L = 6079.2 nb^-1 x 63/64), event-group errors; orange =
    the bump classes C2 + C3.
  - Green: the frozen fit pol2 + h x step at 0.99277 M4max(M(pipi)/0.9970), smearing 5.4 MeV; dashed: without the step.
  - h is printed per panel. f_eta(edge) = sum h_data / sum h_sim = 0.392.
- **Fig.** `wp6_tag.pdf`. *ECAL tag.*
  - (a) Data S0, ECAL-tagged rows: N - CB in m5c with the frozen fit (0.5463 / 6.2 MeV, pol2).
  - (b) The emulated simulated tag (UNVERIFIED): T1 + T2 with the same row definition; orange = C2 + C3.
  - (c) Simulated M(e+e-gamma) of the bump classes with a good photon.
    - Blue: the data pi0 window 0.119-0.173. Red dashed: the window re-centred on the simulation.
    - Dotted: the data pi0 position 0.1459.
- **Table E3** `docs/latex/tables_all/wp6_tableE3.{txt,tex}`: estimators (value, stat, MC, total, bias-corrected, status),
  closure (sigma_tr, CR-L/CR-K/below-B pulls and the S10 terms), cross-checks, and all variants.
- Logs:
  - `/hdd3/przygoda/HADES/FAT_out/all_report/wp6/wp6_fill.txt`: edge-step tables per slice and class, tag fits, P_main; KEY
    lines.
  - `/hdd3/przygoda/HADES/FAT_out/all_report/wp6/wp6_ana.log`.

---------------------------------------------------------------------------------------------------
## 4. Acceptance criteria and fall-backs

WP6 has no numeric criteria of its own. Checks and gates:

| item | value | verdict |
|---|---|---|
| inputs reproduced (windows, CB, replicas, F-main, E_step, T, P_main) | exact (Sect. 2.1) | PASS |
| P_main on S0 = 4388 +- 307 +- 1078 | 4387.8 +- 306.6 +- 1078.2 | PASS |
| WP4 gate, peak side bands | PASS (WP4) | main result (peak) |
| WP4 gate, B4 | FAIL marginal; PASS with bias correction (WP4) | reported both ways; analyst's call |
| stage A at the WP5 eta: chi2/ndf <= 1.5 | GLS 1.46 (WP3 smoothing); 1.56 with the WP4 smoothing | PASS (marginal) |
| displaced-window closure | rms 1.14 -> sigma_tr 239 rows (4.7 %) | S10 (peak) entry |
| CR-L / CR-K prediction pulls | -1.88 / -1.98 (0.400-0.455: -2.78) | within 2 sigma but one-sided; S10 (bump) 37.5 % |
| f_eta(edge) vs a_eta(F-main) | 0.392 +- 0.049 vs 0.382 +- 0.043 | consistent |
| S5 null test of B4 | -117 +- 56 | -2.1 sigma (same sign as the CR-K pull) |

**Fall-backs:**
1. **"F-main fails WP4 -> B4 + peak side bands are the main result"**: fired upstream (WP4/WP5). Applied: Table E3 labels B4 and
   the side bands MAIN and F-main CROSS-CHECK.
2. **"A1/A2 fail -> drop CR-K, b from CR-U only, bump EXTRAPOLATED, add the CR-L term to S10"**: fired in WP3. The CR-L term
   is computed here: 499 rows in CR-L = 7.7 % of B_A, i.e. **2931 rows = 37.5 % of Y_bump(E2)** for S10 (bump).
3. **"k_E1 outside 0.8-3 -> drop E1"**: fired in WP3/WP5. B4 is evaluated on E2 and extrapolated to all M(ee) with the
   simulated f_E2.
4. **"MC error > data stat -> simulation-limited"**:
   - fires for B4 (MC/stat = 4.2);
   - does not fire for the side bands (MC/stat = 0.57).
5. **Not fired:**
   - |rho(a_eta, b)| > 0.7 or a_eta ~ 0 with a significant E_step: a_eta = 0.38 and f_eta(edge) = 0.39 agree;
   - stage-A chi2/ndf > 1.5: GLS 1.46 with WP3's smoothing.

**Deviations from the spec (closest sound choice; flagged):**

| # | deviation | reason |
|---|---|---|
| a | B4 on E2, extrapolated to all M(ee) with f_E2 | Stage A exists only in E2 (E1 fall-back) |
| b | U = 0.575-0.700 instead of 0.565-0.700 | The adopted CR-U. Variant: +0.6 % |
| c | eta_U subtracted in U, as WP4 implemented it | The spec-literal formula gives -24 % because the eta content (mostly G_eta) is 4.7 % of (N - CB) in U |
| d | "other eta" in B with the F-main kappa and rho (WP4) | |
| e | Stage A re-run inside WP6 | Not with WP3's macro (forbidden to re-run). Same configuration and copied smoothing; reproduces WP3 within 0.1 sigma at a_eta(0) = 0.4 |
| f | Side-band leakage in the displaced windows fixed to the nominal peak's f | (C4 + C5) is about 0 there, so a local f is undefined |

---------------------------------------------------------------------------------------------------
## 5. Caveats

1. **The bump background is extrapolation-limited.**
   - The linear tilt fitted in CR-U predicts too much below the peak: CR-L -6 %, CR-K -25 %, 0.400-0.455 -14 % (data/prediction
     - 1).
   - A quadratic term over-corrects: Y_bump(E2) 7822 -> 19526.
   - The true background in B lies between the two. Linear interpolation to zero CR-L pull would give about 12900, but that is
     information only.
   - The S10 entry from CR-L (37.5 %) follows the spec's fall-back. The region below B suggests up to 90 %.
   - F-main (19426 all M(ee), about 9300 in E2) and B4 agree within these errors.
   - Recommendation: quote the absolute bump yield with this S10, and use f_eta(edge) (no background normalisation) as the
     robust bump-scale check.
2. **B4 is simulation-limited.**
   - MC 1920 = 25 % of Y_bump(E2) (WP4 found 9 % for all M(ee)). The cause is the MC error of the smoothed background shape
     ratio B_A(B)/B_A(U), about 5 %.
   - The replica mean is -10 % below the nominal (template attenuation, as in WP4/WP5); bias-corrected 8638.
3. **Stage A at the WP5 eta (a g = 0.66 on G_eta)** fits CR-U worse in the GLS chi2 (1.46) than WP3 at 0.4 (1.10).
   - My GLS covariance holds only the smoothed-background replicas, not the raw G_eta template noise, so 1.46 is an upper
     estimate.
   - The WP4 smoothing (start at 0.400 instead of the kinematic threshold) gives 1.56 and +4.6 % on B4.
4. **The side bands are robust to the background model** (the stage-A part of the error is 6 rows). They are sensitive to the
   side-band width (narrow/wide/flat: -7 %, within the correlated stat error) and to the leakage correction (+22 %; it rests
   on the simulated peak shape, whose width is 6.1 MeV in the simulation against 5.4 MeV in the data).
5. **Emulated simulated ECAL tag: UNVERIFIED.**
   - The simulated photon scale (-4.7 %) and the data scale (+8 %) differ by 13 %.
   - The 100 MeV cluster-energy cut and the ECAL efficiency are not validated in the simulation.
   - Re-centring the window changes f_tag by -9 % (frozen fit) to -14 % (free fit).
6. **Simulated edge-step fits** use the data-like Poisson weights of the expected rows ("the same fit"). Their chi2 (500-1500
   per 22-35 ndf) reflects the MC fluctuations; the MC error comes from the 100 replicas.
7. **P_main data/sim** uses the frozen data shape on the simulation, whose free fit gives 5.2 MeV at r = 1 (D8: 6.1 MeV). It has
   no MC error; information only.
8. **DR from the side bands** (1.40 +- 0.15 +- 0.53) is dominated by the B4 MC error. DR in E2 is 1.73 +- 0.65. F-main's DR is
   1.17 +- 0.15.
9. **Blinding (D19/D20):** only the coarse M(ee) classes (10 MeV split) and the tag cut M(ee) < 0.10 were used. No X2D or B_mee
   cell and no 2 MeV M(ee) shape in the peak band was formed.
10. **Correlations neglected,** as in WP4:
    - between the stage-A (or F-main) parameters and the data counts of the estimator windows (they share the U and E2 data);
    - between N and CB.
    WP4 toys showed pull widths of 0.94-1.06 for this treatment.

**For WP7 (S10 and others):**
- peak: sigma_tr 239 rows (4.7 %) + WP4 (a) 0.5 %;
- bump: CR-L 37.5 % + WP4 (a) 2.6 %;
- S3: quadratic +150 % (bump), +3.4 % (peak);
- S4: WP4 smoothing +4.6 % / -0.4 %;
- S2: CR-U 0.565 +5.0 % / +0.1 %.

The macro refits any WP0 variant: `BuildCtx(tag, lo, hi, drop, quad, withF, withRep, mode)`.

---------------------------------------------------------------------------------------------------
## 6. Exact commands (from /home/przygoda/HADES/FAT/FAT)

    B='gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_wp6", kTRUE)'
    root -l -b -q -e "$B" 'plots/all_report/wp6_fill.C+O' > /hdd3/przygoda/HADES/FAT_out/all_report/wp6/wp6_fill.log 2>&1      # 40 s, about 10 GB
    root -l -b -q -e "$B" 'plots/all_report/wp6_ana.C+O'  > /hdd3/przygoda/HADES/FAT_out/all_report/wp6/wp6_ana_stdout.txt 2>&1 # 20 s, about 3 GB

**Code** (`plots/all_report/`):
- `wp6_fill.C`: rows, edge step, tag, P_main, simulated window sums with replicas.
- `wp6_ana.C`: F-main refit, stage A, B4, side bands, closure, Table E3, figures. It includes `wp5_lib.h` (and through it
  `wp4_lib.h`).

**Outputs:**
- `/hdd3/przygoda/HADES/FAT_out/all_report/wp6/{wp6_rows.root, wp6_fill.root, wp6_fill.txt, wp6_fill.log, results.txt,
  wp6_ana.log, wp6_ana_stdout.txt}`;
- `docs/latex/tables_all/wp6_tableE3.{txt,tex}`;
- `docs/latex/figs_all/wp6_{E4,edge,tag}.{pdf,png}`.
