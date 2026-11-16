# vw2a: independent reproduction of the WP2 material share f_win (empty target)

Verdict: **CONFIRMED WITH CHANGES.**
- The WP2 number is reproduced exactly with independent code: f_win(N - CB, S0, B+P) = 16.76 +- 1.88 %.
- It is the material share of the selected N - CB rows. Quote it as 16.8 +2.8/-3.2 % (stat and syst).
- It is **not** the material share of the eta signal. The full-target data, split in vertex-z slices, give an
  eta-specific share f_eta = 5 +- 7 % (1-8 % depending on the model).
- So the 17 % reduction of sigma(pp -> eta X) is not supported. Use 1 - f_eta = 0.95 +0.05/-0.07, not 0.832.

Files are in `.../scratchpad/final/verify_wp2/vw2a/`:
- macros: vw2a_slim.C, vw2a_main.C, vw2a_extra.C, vw2a_eta.C, vw2a_pat.C, vw2a_patnorm.C, vw2a_fig.C;
- outputs in out/: vw2a_main.txt, vw2a_extra.txt, vw2a_eta.txt, vw2a_patnorm.txt, vw2a_main.root;
- figure: out/vw2a_fig.pdf (and .png).

All code is my own. No WP2 code and no ScanEngine are used. The definitions follow defs.h.
- Event groups: adjacent rows with equal (file_id, run, grp).
- CB: 2 sqrt(N++ N--) in 10 MeV cells, linearly interpolated (the frozen definition), plus two variants.
- Errors: 200 Poisson bootstrap replicas of the event groups, in both samples.

## 1. Flux ratio b (normaliser)
- **Reproduced.** b = U_full/U_empty, base rows us+pp+mm, vz in [-200,-170):
  - 10.350 +- 0.787 (034), 15.444 +- 0.921 (047), 0.932 +- 0.030 (065);
  - identical to WP2.
- **Neighbouring regions** (b / b_nominal for 034 / 047 / 065):

  | region | 034 | 047 | 065 |
  |---|---|---|---|
  | [-230,-170) | 1.007 | 0.998 | 1.005 |
  | [-250,-170) | 1.011 | 0.999 | 1.004 |
  | [-200,-180) | 1.17 | 1.10 | 0.96 |
  | [-190,-170) | 0.95 | 0.96 | 1.01 |
  | [-180,-170) | 0.93 | 0.96 | 1.02 |
  | isBest == 1 rows | 1.03 | 0.97 | 1.01 |
  | us only | 1.02 | 0.96 | 0.99 |
  | like-sign only | 0.96 | 1.11 | 1.02 |
  | both pions in the vertex fit | 1.52 +- 0.30 | 1.15 +- 0.18 | 1.11 +- 0.07 |

  f_win(B+P) with these regions: 16.0-18.5 %; 22.4 % with the pion-vertex rows.
- **Other hypothesis trees** (PAT rows, PT3, isBest == 1, same window), b / b_ch:
  - EpEm_ID 1.36 / 1.16 / 1.32. Its far-upstream F/(bE) is 2-2.8, so its rows contain LH2 tails. This agrees with WP2,
    which rejected it.
  - Three-track hadronic trees (PipPimEp, PipPimEm, PEpEm, EmEm): 1.03-1.21. On day 065: 1.08-1.11 +- 0.015, i.e. 3 sigma
    above the channel.
  - The channel's own PAT rows: 0.96-1.04.
  - A b about 10 % higher conflicts at about 2 sigma with positivity at [-170,-160) in the channel on day 065
    (F/E = 1.04 b).
- **Shape in front of the cell.** F/(bE) per 10 mm bin over [-260,-160):
  - flat within the errors: slope 1.8-4.9e-3 per mm, 1-2 sigma;
  - the excess starts only at [-170,-160): 1.15, 1.25, 1.04.
  - Radial share r < 4 mm in [-250,-170): full 0.005, empty 0.000. The normaliser rows are off axis.
- **LH2 leak into U.** The pion-vertex fraction of the window rows is 0.156 (full) against 0.142 (empty). That is at most
  about 2 % LH2 in U, so f_win changes by -2 %.
- **The normaliser follows the beam halo, not the luminosity.**
  - h = U_full/L is 31.4 (034), 35.6 (047) and 18.2 (065).
  - Per period: I 24.2, II 32.4, III 18.6. Within period I it ranges over 21.3-31.4.
  - The yield per luminosity Y/L stays flat within +-5 %.
  - Within one day, in-cell rows per U (C/U) change by 20-40 % (chi2 553/5, 2439/5 and 166/5 for blocks of the full runs).
  - The C/U of the full runs nearest to the empty runs gives timing factors on b: 0.95 +- 0.05, 1.03 +- 0.05,
    0.99 +- 0.06. Net effect on f_win: -2 %.
  - Day-034 empty runs: C/U is not stable (chi2 16.5/6, p = 0.011). Days 047 and 065 are stable (p = 0.29, 0.73).

## 2. f_win per day and combined (S0, N - CB, window CB as WP2)

| window | 034 | 047 | 065 | Lper (WP2) | Y-weighted | nearest day | pooled 3 days |
|---|---|---|---|---|---|---|---|
| all M4 | 12.31 | 12.84 | 10.39 | 11.86 +- 0.59 | 11.89 | 11.98 | 11.86 |
| B | 21.06 | 21.73 | 12.44 | 18.64 +- 2.03 | 18.73 | 18.79 | 18.37 |
| P | 20.75 | 10.12 | 12.61 | 15.78 +- 3.12 | 15.77 | 13.33 | 13.30 |
| B+P | 19.29 | 17.24 | 12.47 | 16.76 +- 1.88 | 16.81 | 16.29 | 16.05 |
| 0.40-0.70 | 13.66 | 14.20 | 11.42 | 13.12 +- 0.78 | 13.14 | 13.24 | 13.09 |

- Every WP2 variant is reproduced to 0.01: pre 22.44, S0 isBest==1 16.48, vz150 12.36, vz170 19.00, dz40 16.63, N 12.79,
  like-sign 29.67.
- **Weights.** WP2 uses the analysed luminosity per period: 2809.452 / 1450.632 / 1819.110 nb^-1. This is correct.
  Weighting by the full-target yield gives the same result.
- **Problem: the halo level of the days.**
  - f_d rises with h: slope 0.142 +- 0.051 %/unit for all M4 (2.8 sigma) and 0.36 +- 0.17 for B+P.
  - Day 034 has the highest h of period I (31.4 against 24.2), so it does not represent period I.
  - A linear model evaluated at the L-weighted mean h = 24.5 gives:
    - all M4 11.29 +- 0.37 (-5 %);
    - B+P 14.75 +- 1.15 (-12 %);
    - B 16.2, P 12.5, 0.40-0.70 12.45.
  - Scaling each day to its period's halo level gives B+P 14.6.

## 3. Where the material is (Fig. (a), (b))
Unlike-sign rows, both pions in the vertex fit, days 034 + 047 + 065:
- **Front material** at -200..-165 mm: off axis (r > 4 mm), hit by the beam halo. This is the normaliser region.
- **Like-sign-rich narrow peak** at -157.5 mm, also present in the empty runs. It is the source of the CB-rich rows seen
  in the data report.
- **Entrance window** at -149.5 mm: 30 % of its rows have r < 4 mm.
- **Exit window** at -88.5 mm: beam-core-like, r < 4 mm share 0.46 against 0.42 for LH2.
- **In-cell plateau** at 1.2 % of the LH2 rows ([-130,-115)).
- **Material share** (empty x b / full): 0.33-0.93 at the entrance (-170..-145), 0.9-3 % in the cell centre and 0.4-0.57
  at the exit (-95..-80).
- **f_win by region** (B+P): entrance 6.55, cell 5.88, exit 4.68. This is the same as WP2.
- **Gas.** If the plateau is H2 vapour, it is about 1 % of the LH2 density. Its bound is then at most about -1 % absolute.
  WP2's -3.6 overstates the downward band.

## 4. Eta in the material
- **Empty sample directly.**
  - Pmain-like fit on the scaled empty N - CB: S = 58 +- 98, against the full-target S = 221 +- 97 on the same three days.
    So f_eta = 0.26 +- 0.83: no sensitivity.
  - Bump edge: J_full = 1.81 +- 0.03, J_mat = 1.26 +- 0.20. f(M4) is higher below the edge by +6.1 +- 3.6 % absolute
    (the material spectrum is softer). Inconclusive.
  - An edge-step fit in M4 - E(M(pipi), M(ee)) failed: chi2/ndf far above 1.
- **Full data, all 29 days, in vz slices (Fig. (d)).**
  - Fit model: S_k = s_L (1 - m_k) Y_k + s_M m_k Y_k, with m_k from the empty days.
  - r = s_M/s_L = 0.07 +- 0.30 (chi2 8.5/10).
  - Split by region:
    - entrance r = -0.27 +- 0.34;
    - cell r = 0.2 +- 1.2;
    - exit r = 1.21 +- 0.64.
  - f_eta (vz > -160) = 1.1 +- 4.7 % (one material coefficient), 2.8 +- 4.9 % (entrance + cell / exit),
    7.9 +- 5.6 % (entrance / cell + exit) and 4.7 +- 7.0 % (three regions).
  - The N - CB material share of the same fit range is 13.3 %.
  - The rows from the entrance region carry no eta peak. The exit-foil rows are consistent with eta at the LH2 rate.

## 5. CB of the empty sample
- CB/N of the empty sample (window CB): 0.10-0.13 (all M4) and 0.13-0.16 (B+P). For the full target: 0.04 and 0.06.
- The material supplies 29.7 % of the like-sign rows (all M4) and 44 % in B+P.
- **Frozen cell CB.** Biased low at the empty counts: toys give mean/truth 0.44-0.61 on 034/047 and 0.92-0.96 on 065.
  It would give B+P 18.08 (+8 %). WP2 rightly used the window CB.
- **Window CB.** Residual toy bias 0.79 / 0.90 / 0.99 (B+P), i.e. +1-2 % on f. The extreme case is window B on day 034,
  where N++ = 0 gives CB = 0.
- **Kappa CB calibrated on the full data** (B+P factor 0.805): B+P 16.31 (-3 %).
- **Is N - CB the right quantity?** It is right for the background composition. It is not right for the eta yields: the
  fit absorbs the eta-free material background, so only f_eta applies (Sect. 4).

## Problems (size in % of f_win = 16.8, B+P)

| # | problem | effect on f_win |
|---|---|---|
| 1 | used for the eta yields instead of f_eta | about -70 % (16.8 -> 5 +- 7) |
| 2 | halo level of day 034 not representative of period I | -12 % (all M4: -5 %) |
| 3 | normaliser choice: hadronic trees / pion-vertex U / sub-windows | +3 to +13 % / -4 % |
| 4 | CB method | window-CB bias +1-2 %; kappa variant -3 % |
| 5 | gas bound too wide | -21 % (WP2) -> about -6 % |
| 6 | LH2 leak into U | about -2 % |
| 7 | timing of the empty runs | -2 % net; per day +-5-15 % |
| 8 | day-034 empty runs unstable (p = 0.011) | inside the day-034 statistical error (19 %) |
