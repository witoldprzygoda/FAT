# WP2 results: normalisation items (L content, C_sec, f_win)

Package WP2 of SPEC_report_all. First attempt 05:43-06:28 (it hung on an rm approval); its outputs were checked and reused. This
continuation finished the empty-target chain, the f_win measurement, Fig. N1 and Table N1.

**Headline**
- L already contains the sector-2 removal for PT2 and the elastic target-window removal. It is LH2 only.
- C_sec for the PT3 channel is small: 1.010-1.012 (all periods) and 1.022-1.027 (period I only).
- The empty-target (material) share f_win of the pi+pi-e+e- rows is LARGE and measured directly:
  - all M4: 11.9 +- 0.6 %;
  - bump window B: 18.6 +- 2.0 %;
  - peak window P: 15.8 +- 3.1 %;
  - recommended for the eta yields: 16.8 +5.0/-5.9 %.
- This is 8x the spec's fall-back range (0.7-2.1 %). sigma(pp -> eta X) gets C_sec x (1 - f_win) = 0.842 +0.060/-0.051.

---------------------------------------------------------------------------------------------------
## 1. What was done (plain words, step by step)

### 1.1 What L contains
- **Sources.** Header of the per-run luminosity CSV (`remote/lumi_startfree_perrun.csv`, chain "N", total 6080.5) and the
  Task E documents of the luminosity analysis on vae25 (`lgadfree_lumi/E/README_E.txt`, `N/EW/window_fraction_*.txt`).
  Copies are in `.../scratchpad/final/wp2work/taskE/taskE_dump.txt`.
- **L contains:**
  - the PT2 DAQ live fraction and the beam duty factor;
  - S_2: the sector-2 PT2 (elastic) trigger deficit of days 034-044, removed by replacing sector 2 with rho x the mean of the
    other sectors (S2fac 1.058-1.120 per day);
  - (1 - w6): elastic events from the target windows removed (w6 = 0.702 / 0.806 / 0.806 / 0.883 % in I / IIA / IIB / III).
  So L is LH2 only.
- **L does NOT contain:**
  - the PT3 word factor 64/63;
  - the PT3 live fraction (one-sided +7.3/-0 %, not applied, D17);
  - the sector-2 deficit of the PT3 channel rows (-> C_sec);
  - the material contribution to the PT3 channel rows (-> f_win).
  The CSV says so verbatim: "Sector 2 on days 034-044 has a hardware trigger deficit that is NOT in L" and "Yields must be
  empty-target subtracted (L is LH2 only)".
- **How w6 was measured (Task E).**
  - Same-day empty/full method on days 034, 047 and 065, applied to the periods.
  - Normaliser U1 = inclusive tracks with z in [-700,-400) mm.
  - Tight elastic cuts, so the elastic window share is only 0.40-0.46 % raw.
- **L of the analysed runs** (`wp2_lumi.py`, first attempt): 6079.195 nb^-1 (17 557 runs).
  - Per period: I 2809.453, IIA 372.610, IIB 1078.022, III 1819.110.
  - Deficit span 034-045: 2619.622 (034-044: 2449.224).
  - Clean days 046-067: 3459.573 (045-067: 3629.972).

### 1.2 Sector-2 deficit of the channel (C_sec), first attempt, checked
1. **Sector numbering.** sector = (floor(phi/60) + 5) mod 6.
   - Check 1: the acceptance gaps of the S0 e+e- tracks sit at k x 60 deg (content in +-1 deg / side bands = 0.001 at all six
     boundaries).
   - Check 2: the rate deficit of the early days lies exactly in 180-240 deg, which this formula calls sector 2 (Fig. N1a).
2. **Rates.** Rows per nb^-1 per day and per sector used (S0, pre, none; N and N - CB): `wp2_sector.txt`.
3. **Deficit efficiency eps_s2.**
   - Definition: rows that use a sector-2 track, compared between deficit days and clean days, as a DOUBLE ratio against rows
     without a sector-2 track. This removes L and global drifts.
   - Blocks of 25 nb^-1 (Fig. N1b) show the deficit is intermittent: 30 blocks on days 034-045 have eps < 1 - 4 sigma, with
     eps down to 0.58. Elsewhere there is only one such block (day 061).
   - Day 045 still has three deficit blocks (00:00-02:07, 02:08-04:05, 07:27-17:56; day eps 0.903 +- 0.007). Day 046 has none
     (1.008 +- 0.007).
   - So the nominal deficit span is 034-045 (the spec says 034-044). eps = 0.9266 +- 0.0024 (034-045 vs 046-067). Spec span:
     0.9290 +- 0.0025.
4. **Signal share with a sector-2 track, f_s.** From the simulation (S0, V3, mode A, M4' scale; `wp2_simfs.C`). Bump C2+C3 in B:
   0.326 +- 0.012. Peak C4+C5 in P: 0.362 +- 0.016.
5. **Correction.** C_sec = 1/(1 - f_s (1 - eps) L_e/L) (exact form; the spec's first-order form is +0.07 % higher).
   - Systematics: the loss computed per sector pattern (B0..B3) or per number of sector-2 tracks (signal mixes differ from
     inclusive rows); f_s from the data of the clean days; span 034-044. These are added in quadrature.
   - Period I alone: C_sec(I) with L_I.

### 1.3 Empty-target runs and the material share f_win
1. **Runs found.** Empty-target DSTs exist in `/lustre/hades/dst/feb22/gen4/4500MeV/3195A/emptytarget/` for days:
   - 034: 7 runs, 70 DSTs;
   - 047: 16 runs, 160 DSTs;
   - 065: 104 runs, 1040 DSTs.
   These are the same days the luminosity analysis used.
2. **PAT** (first attempt, on vae25 with SLURM arrays).
   - Binary built from the production tag `prod-feb22-betanew-2026-10-28`. main.cc, core/ and user/ are identical to the tag
     (checked here with diff). Same PID file (md5 b8f0cc2f...).
   - All 44 jobs ended with SUCCESS, 0 failed files.
   - Reproducibility (`wp2_patcheck.C`, run here): one full-target DST (be2203405213001) gives 38/38 PipPimEpEm_ID and
     764/764 EpEm_ID rows identical to the production file lep_034.root.
3. **Copy and merge.**
   - 171 day-065 outputs were missing or truncated locally (two files of 388 bytes). They were copied to `empty/pat_parts2/065`.
   - Sizes were checked against vae25 for all 1270 outputs.
   - Day 065 was merged into `empty/pat/lep_empty_065.root` (1040 inputs). Days 034 and 047 had been merged by the first
     attempt.
4. **FAT** (HEAD 63a9a76, build of the first attempt in `wp2work/fat_head`, main.cc and src identical to HEAD).
   - Same config as the data production (T 4530, PT3, corrected), on the three empty-target day files.
   - Outputs go to the NEW directory `empty/fat2` (`empty/fat` is the broken one). scan_nt rows: e+e- 31 479, e+e+ 2 250,
     e-e- 6 196.
5. **Compact file.** `wp2_scan_prepare_empty.C` (copy of the data scan_prepare.C with only the paths changed, now pointing to
   fat2/scan2) -> `empty/scan2/scan_compact_empty.root`. Internal read-back checks pass (0 mismatches).
6. **Flux normalisation b (full/empty, same day).**
   - The beam hits the material in front of the cell (vz -200..-170 mm) in both kinds of runs. So b = U_full/U_empty, where U
     counts the channel rows (base selection: isBest >= 0, vz > -500, PT3, cor_ok; us + pp + mm) with vz in [-200,-170) mm.
   - The simulated LH2 puts 1.7e-5 of its rows there, so the window is LH2-free.
   - b = 10.35 +- 0.79 (034), 15.44 +- 0.95 (047), 0.932 +- 0.030 (065).
   - The radial vertex profile in this window is the same in full and empty runs (mostly r > 8 mm). This supports "same
     material".
7. **f_win.** f_win = b x Y_empty / Y_full for every selection, M4 window and quantity (N, like-sign, N - CB).
   - Combined over the three days with the luminosity of the period each day represents: 034 -> I, 047 -> IIA + IIB,
     065 -> III. This is the Task E convention.
   - Errors: event-group variances of Y and U in both samples.
8. **Cross-check of b without a window: vz-profile fit** (the spec's fall-back method, used here as a test; `wp2_vzfit.C`).
   - F(z) = a x LH2-sim(z; shift, smearing) + b x empty(z), per day, in [-250,-60) mm.
   - b_fit = 8.76 +- 0.55, 13.87 +- 0.57, 0.816 +- 0.023, i.e. 0.85-0.90 of the window value. chi2/ndf = 1.6-2.7: the
     simulated vertex resolution is not perfect.
   - Narrow range [-220,-70): the same.
   - Without smearing: b closer to the window value (11.7, 14.9, 0.884).
9. **Upper bound on b (positivity).** In the entrance-window bin [-170,-160) the full-target rows are at least the material.
   So b <= F/E = 12.59, 17.99, 0.992 (1.22 / 1.17 / 1.06 x b).
10. **Rejected normalisers.**
    - EpEm_ID rows in front of the cell give b larger by 1.16-1.36 (near window) and 1.65-2.37 ([-700,-250) mm).
    - Both violate the positivity bound: subtracting would make the full-minus-material vz profile negative at the entrance
      window.
    - For day 065, the channel's own far-upstream rows agree with the near window (b 0.900 +- 0.106 vs 0.932). So the EpEm
      excess is LH2-related (2-track vertices with long z tails), not flux.
    - Task E's U1 normaliser, taken with an exact 1/20 sampling of its full-target S1 sample, corresponds to b = 11.8 / 16.2 /
      1.09. That is close to the channel normaliser (+14 / +5 / +17 %) and inside the band of item 9.
11. **Gas in the "empty" cell.** If the empty cell held H2 gas, its rows (absent in the full runs, where the liquid replaces
    the gas) would be wrongly subtracted.
    - The in-cell empty rows have an on-axis (r < 4 mm) share of 0.227-0.335, against 0.54 for LH2 and 0.00 for the front
      material. So at most a fraction g = 0.62 of them is beam-like (gas or on-axis foil tails).
    - g x (in-cell part of f_win) is taken as a one-sided bound (`wp2_wincomp.C`).
12. **Composition by vertex region (B+P, N - CB).** Entrance [-160,-145): 6.5 %; in-cell [-145,-105): 5.9 %; exit [-105,...):
    4.7 %. This is consistent with Fig. N1c: two window peaks plus a low in-cell plateau.
13. **The eta-specific share cannot be measured.** The empty-target peak signal is 12 +- 50 rows. f_win of the eta yields is
    therefore taken from the N - CB share in B + P, and the difference to the 0.40-0.70 range is added as a "composition"
    systematic.

---------------------------------------------------------------------------------------------------
## 2. Numbers (all S0 = pre + both pions in the vertex fit, isBest >= 0, unless stated)

### 2.1 C_sec (stat, syst; nominal span 034-045; Table N1)

| signal | f_s (sim) | C_sec all periods | C_sec period I only |
|---|---|---|---|
| bump C2+C3 in B | 0.326 +- 0.012 | 1.0104 +- 0.0005 +- 0.0035 | 1.0229 +- 0.0077 |
| peak C4+C5 in P | 0.362 +- 0.016 | 1.0116 +- 0.0006 +- 0.0029 | 1.0254 +- 0.0065 |
| C2 in B | 0.320 | 1.0102 +- 0.0035 | 1.0224 |
| C3 in B | 0.335 | 1.0107 +- 0.0035 | 1.0235 |
| C4 in P | 0.340 | 1.0109 +- 0.0036 | 1.0238 |
| C5 in P | 0.382 | 1.0122 +- 0.0026 | 1.0269 |
| T12 in CR-U | 0.444 | 1.0143 +- 0.0016 | 1.0314 |

- Spec variants:
  - first-order form: bump 1.0112, peak 1.0124;
  - span 034-044: bump 1.0094, peak 1.0105.
- Ratios: C_sec(peak)/C_sec(bump) = 1.0011 (DR); C5/C2 = 1.0020 (golden ratio). Both are 1 within 0.2 %, so no correction is
  needed for ratios.
- Variant "clean days only": L(046-067) = 3459.573, L(045-067) = 3629.972 (with their own yields, C_sec = 1).
  - Indicative check with the WP5 period fits: a_eta(IIA+IIB+III) = 0.380 +- 0.011 against a_eta(I) x C_sec(I) =
    0.376 +- 0.013, consistent.
  - The raw ratio I/(II+III) = 0.966 +- 0.042 against 1/C_sec(I) = 0.976.

### 2.2 f_win (N - CB, %, L-weighted over the 3 days; normaliser ch[-200,-170))

| window | f_win +- stat | fit-b (lower) | bound-b (upper) | day 034 | day 047 | day 065 |
|---|---|---|---|---|---|---|
| all M4 | 11.86 +- 0.59 | 10.29 | 13.80 | 12.31 +- 1.15 | 12.84 +- 0.98 | 10.39 +- 0.36 |
| bump B | 18.64 +- 2.04 | 16.15 | 21.84 | 21.06 +- 3.86 | 21.73 +- 3.98 | 12.45 +- 0.94 |
| peak P | 15.78 +- 3.07 | 13.58 | 18.49 | 20.75 +- 6.29 | 10.12 +- 3.79 | 12.61 +- 1.22 |
| B + P | 16.76 +- 1.88 | 14.50 | 19.61 | 19.29 +- 3.69 | 17.24 +- 3.16 | 12.47 +- 0.81 |
| 0.40-0.70 (fit range) | 13.12 +- 0.85 | 11.38 | 15.27 | 13.66 +- 1.67 | 14.20 +- 1.44 | 11.42 +- 0.46 |
| CR-U | 12.29 +- 0.87 | 10.66 | 14.30 | 12.69 +- 1.70 | 13.24 +- 1.49 | 10.92 +- 0.47 |

- **Excluded EpEm normalisers** (for the record), all M4: near 15.4 %, far 24.1 %.
- **N** (not N - CB), all M4: 12.79 +- 0.60. **Like-sign:** 29.7 +- 2.2 (the material has a higher CB fraction).
- **Selection variants, B+P (N - CB):**
  - S0 16.76 +- 1.88;
  - pre 22.44 +- 2.21;
  - S0 isBest == 1: 16.48 +- 1.95;
  - vz > -150: 12.36 +- 1.58;
  - vz > -170: 19.00 +- 2.02;
  - S0 + dz_pi < 40: 16.63 +- 1.90.
- **Selection variants, all M4:**
  - S0 11.86;
  - pre 16.06;
  - isBest == 1: 10.73;
  - vz150 9.24;
  - vz170 12.86;
  - dz40 11.87.
  The pion-vertex requirement and a tighter vz cut reduce the material. Every efficiency-corrected variant must use its own
  f_win.
- **M4 dependence** (Fig. N1e): f_win is about 0.11-0.12 above 0.56 GeV/c2 and rises to 0.15-0.25 below 0.52 (the material
  spectrum is softer).
- **Recommended for the eta yields** (a_eta, sigma(pp -> eta X)): f_win = 16.8 % +5.0/-5.9 (absolute).
  - stat: 1.9;
  - normaliser: +2.8 (positivity bound) / -2.3 (vz fit);
  - gas: one-sided -3.6;
  - composition: 3.6 (B+P vs 0.40-0.70).
  - So 1 - f_win = 0.832 +0.059/-0.050 (+7.1/-6.0 % relative).
- **Combined normalisation factor for sigma(pp -> eta X)** = a_eta x 1.0665 mb x C_sec x (1 - f_win):
  C_sec(eta) = 1.0110 +- 0.0035 (mean of bump and peak), so C_sec x (1 - f_win) = 0.842 +0.060/-0.051.
- **Ratios** (DR, golden ratio BR_G): (1 - f_P)/(1 - f_B) = 1.035 +- 0.046 (stat). This is consistent with 1, so it is NOT
  applied; quote 4.6 % as a systematic of the window-dependent material share if desired. The eta-specific shares cannot be
  separated.

### 2.3 Per-period inputs for the stability test (WP7; Table N1 "[per period]")

| period | L [nb^-1] | C_sec | f_win (B+P) [%] | K = C_sec (1 - f_win) |
|---|---|---|---|---|
| I (034-046) | 2809.453 | 1.0241 +- 0.0077 | 19.29 +- 3.69 (day 034) | 0.827 +- 0.038 |
| IIA | 372.610 | 1 | 17.24 +- 3.16 (day 047) | 0.828 +- 0.032 |
| IIB | 1078.022 | 1 | 17.24 +- 3.16 (day 047, correlated with IIA) | 0.828 +- 0.032 |
| III | 1819.110 | 1 | 12.47 +- 0.81 (day 065) | 0.875 +- 0.008 |

**Method.** a_eta,p (WP5, background shape fixed, templates x the L share) x K_p. Then chi2 against the weighted mean; pass if
p > 0.01. Do it with and without the per-period f_win (the f_win differences between days are 1.6-2 sigma).

**Indicative, with the WP5 period values** (WP7 must redo it with the correlations):
- raw: chi2 3.22/3 (p 0.36);
- x C_sec: 2.65/3 (p 0.45);
- x K_p: 2.04/3 (p 0.56).

DR is unchanged by these factors to first order (WP5: p = 0.038).

---------------------------------------------------------------------------------------------------
## 3. Figures and tables (draft captions)

**Fig. N1** = `docs/latex/figs_all/wp2_N1v2.pdf` (+ .png).
- The old `wp2_N1.pdf` (two panels, first attempt) was not overwritten (shell-safety rule). Use the v2 file.
- Draft caption: "Normalisation items.
  - (a) Track rate per nb^-1 of the S0 e+e- rows (isBest >= 0), days 034-045 divided by days 046-067, per 10 deg of azimuth,
    normalised to the bins outside 180-240 deg. The deficit is confined to sector 2.
  - (b) eps_s2: rows with a sector-2 track relative to rows without one, per block of 25 nb^-1 of the analysed runs,
    normalised to days 047-067. Shaded: the deficit span 034-045.
  - (c) Vertex z of the unlike-sign rows (isBest >= 0, PT3, cor_ok, both pions in the vertex fit, all M4), days 034 + 047 +
    065. Black: full target. Red: empty target x b (flux ratio of the same day from the rows in the shaded window
    -200..-170 mm). Blue: their difference. Green: LH2 simulation fitted to the full-target profile. Dashed line: the
    vz > -160 mm cut.
  - (d) f_win = material share of N - CB in S0 per M4 window, per empty-target day and L-weighted (stat errors). Grey boxes:
    the range of the flux ratio (vz-profile fit to positivity bound). Hatched: the elastic window fraction w6 already removed
    from L (0.70-0.88 %).
  - L = 6079.2 nb^-1."

**Fig. N1e** = `docs/latex/figs_all/wp2_N1e.pdf` (+ .png).
- Draft caption: "M(pi+pi-e+e-) of N - CB (S0, isBest >= 0), days 034 + 047 + 065.
  - Black: full target. Red: material = sum over days of b_d x empty target (20 MeV bins).
  - Lower pad: f_win(M4). The dashed line is the all-M4 value 11.9 %.
  - Shaded: bump window B and peak window P."

**Table N1**: `docs/latex/tables_all/wp2_N1.txt` (all numbers; machine lines `n1|...`) and `wp2_N1.tex` (booktabs, label
tab:N1, caption included).

**Supporting tables** (all in `docs/latex/tables_all/`):
- `wp2_lumi_perday.txt`, `wp2_lumi_perperiod.txt`: L per day and per group, S2fac, w6, live fractions;
- `wp2_runs.txt`: eps per day, per 25 nb^-1 block, step scan, deficit blocks;
- `wp2_sector.txt`, `wp2_sector_eps.txt`: phi gaps, per-sector rates, eps per category, trigger model, phi profile;
- `wp2_simfs.txt`: f_s per class, window and sector pattern;
- `wp2_csec.txt`, `wp2_csec_summary.txt`: C_sec with variants;
- `wp2_win.txt`, `wp2_win_summary.txt`: normalisers, f_win per selection (S0, pre, none, S0b1, vz150, vz170, dz40), window,
  quantity, day, normaliser;
- `wp2_vzfit.txt`: vz-profile fits;
- `wp2_wincomp.txt`: b bounds, on-axis shares, composition by vertex region, gas bound.

---------------------------------------------------------------------------------------------------
## 4. Acceptance criteria (spec WP2) and fall-backs

| item | result | verdict |
|---|---|---|
| sector formula verified on the phi gaps | gaps at k x 60 deg (ratio 0.001 at all 6 boundaries); deficit exactly in 180-240 deg = sector 2 | PASS |
| rows per nb^-1, per day and per sector, S0 | wp2_sector.txt, wp2_runs.txt | DONE |
| C_sec = 1 + f_s (1/eps - 1) L_e/L | exact form 1.0104 (bump), 1.0116 (peak); spec form 1.0112 / 1.0124 | DONE (deviation: span 034-045, see 5) |
| variant: days 045-067 only with their own L | L given (and 046-067); indicative check consistent (0.380 +- 0.011 vs 0.376 +- 0.013) | DONE (yields: WP7) |
| empty-target runs looked for | found: 127 runs, 1270 DSTs (034, 047, 065) | PASS -> direct method |
| f_win | 11.9 +- 0.6 % all M4; 16.8 +5.0/-5.9 % recommended for eta | DONE |
| vz-profile fit (spec method if no empty runs) | used as cross-check; b 0.85-0.90 of the window value; chi2/ndf 1.6-2.7 | lower bound of the normaliser band |
| stability chi2 p > 0.01 | WP7; indicative p = 0.56 (a_eta x K_p) | inputs and method given |
| Table N1, Fig. N1 | wp2_N1.txt/.tex, wp2_N1v2.pdf (+ wp2_N1e.pdf) | DONE |

**Fall-backs.**
- "C_sec or f_win not determinable -> PRELIMINARY": NOT fired.
- The f_win fall-back "0.7-2.1 % (elastic w6 x A-dependence)": NOT used. The measured value is about 8x larger.
- The rule "f_win subtracted with a 100 % error" belonged to the fit fall-back. Here f_win is subtracted with its measured
  uncertainty, +5.0/-5.9 % absolute.

---------------------------------------------------------------------------------------------------
## 5. Caveats and deviations from the spec (flagged)

1. **f_win is large.** About 12 % of N - CB overall, and 16-19 % in the bump and peak windows of S0. Consequences:
   - sigma(pp -> eta X) and every absolute yield drop by about 17 %.
   - The "background 1.3-1.4 of SMASH" indication (D10) and the WP1 "before" data/sim ratios (Table B1) include this
     material, which SMASH does not simulate. About 12-19 % of the data N - CB is not LH2.
   - The material M4 spectrum is softer: f_win(M4) is 0.15-0.25 below 0.52. It is absorbed by the c_j (1 + b x) background
     model of the fits, but it changes the physical meaning of the background renormalisation, also in CR-K.
   - Ratios (DR, BR_G) are unaffected within 4.6 % (stat).
2. **Deficit span 034-045 instead of 034-044.** The PT3 data show deficit blocks on day 045 until about 18:00 (Fig. N1b,
   wp2_runs.txt). The elastic analysis put S_2 only on 034-044. The spec span is given as a variant (C_sec 0.1 % lower) and
   enters the C_sec systematic.
3. **eps is a double ratio** (rows with a sector-2 track relative to rows without). The spec's single early/late ratio is
   0.9538 +- 0.0017 (S0, 034-044 vs 045-067) and is in wp2_sector.txt. The double ratio removes L and global efficiency
   drifts.
4. **The flux normaliser is the main systematic of f_win.**
   - The nominal window (channel rows at -200..-170 mm) is LH2-free in the simulation and has the same radial profile in full
     and empty runs.
   - The vz fit gives 13 % lower b; the positivity bound gives 16 % higher b. Both are in the error.
   - EpEm-based normalisers are rejected (they violate positivity).
   - Task E's U1 normaliser of L would give values inside the band. This assumes an exact 1/20 sampling of its S1 full
     sample, which was not verified.
5. **Gas in the empty cell.** Its state (vapour pressure, temperature) is unknown. The radial check limits the on-axis part of
   the in-cell empty rows to g <= 0.62. This gives a one-sided -3.6 % (B+P) bound. A run-log entry on the empty-cell state
   would remove this term.
6. **The eta-specific material share is not measurable.** The empty peak signal is 12 +- 50 rows. B+P N - CB is used as the
   proxy, plus a 3.6 % composition term.
7. **Day 034 empty runs** were taken 01:26-01:47, before the first full-target block of that day. The state of the
   intermittent sector-2 deficit during them is unknown. With 46 % of the rows using sector 2 and eps >= 0.6, this could
   change f_win(034) by up to about 15 % relative. That is below its 19 % statistical error.
8. **Per-period f_win** comes from one day per period, and IIA and IIB share day 047. Day-to-day differences (B+P: 19.3, 17.2,
   12.5) are 1.6-2 sigma, so the combined value is recommended. The per-period values are for the stability test only.
9. **PAT for the empty runs** was rebuilt from the production tag. Reproducibility was checked on one full-target DST only
   (100 % identical rows). FAT: the same HEAD build (main.cc and src identical to HEAD 63a9a76).
10. **Files.**
    - `/hdd3/.../wp2/empty/fat` (broken self-link) and `empty/scan/scan_compact_empty.root` (675 bytes, empty) are stale
      first-attempt leftovers. They were not touched; the valid outputs are fat2/ and scan2/.
    - `empty/pat_parts/065` contains two truncated files (be2206508405404, be2206509124003, 388 bytes). Their good copies are
      in pat_parts2/065.
11. **Remote clean-up NOT done.** `/lustre/hades/user/przygoda/claude_work/wp2_empty` (7.2 GB) is still on vae25. An rm
    command can trigger an approval prompt, which hung the first attempt, so it was skipped. Everything needed was copied:
    - outputs: `empty/pat_parts{,2}`, `empty/pat`;
    - lists, job arrays, logs, PAT sources and tests: `empty/remote_copy`.
    To clean up by hand: `ssh -o RemoteCommand=none -o RequestTTY=no vae25 'rm -rf /lustre/hades/user/przygoda/claude_work/wp2_empty'`.

---------------------------------------------------------------------------------------------------
## 6. Exact commands (from /home/przygoda/HADES/FAT/FAT unless stated)

    B='gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_wp2", kTRUE)'
    # --- first attempt (reused, not rerun) ---
    python3 -I plots/all_report/wp2_lumi.py                      # wp2_lumi_perday.txt, wp2_lumi_perperiod.txt
    root -l -b -q -e "$B" 'plots/all_report/wp2_runs.C+O'       # wp2/wp2_runs_S0.csv
    python3 -I plots/all_report/wp2_runs.py S0                   # wp2_runs.txt, wp2/wp2_blocks_S0.csv
    root -l -b -q -e "$B" 'plots/all_report/wp2_sector.C+O'     # wp2_sector.txt, wp2_sector_eps.txt, wp2/wp2_sector.root
    root -l -b -q -e "$B" 'plots/all_report/wp2_simfs.C+O'      # wp2_simfs.txt (21 s)
    python3 -I plots/all_report/wp2_csec.py                      # wp2_csec.txt, wp2_csec_summary.txt
    root -l -b -q -e "$B" 'plots/all_report/wp2_simvz.C+O'      # wp2/wp2_simvz.root
    #   PAT on vae25: /lustre/.../claude_work/wp2_empty/jobarray_wp2e_<day>.dat (SLURM arrays; copy in wp2/empty/remote_copy)
    # --- this run ---
    E=/hdd3/przygoda/HADES/FAT_out/all_report/wp2/empty
    rsync -a -e "ssh -o RemoteCommand=none -o RequestTTY=no -o BatchMode=yes" --files-from=<missing065.txt> \
          vae25:/lustre/hades/user/przygoda/claude_work/wp2_empty/out/065/ $E/pat_parts2/065/
    hadd -j 6 $E/pat/lep_empty_065.root @<list065.txt>           # 869 files of pat_parts/065 + 171 of pat_parts2/065, sorted
    (cd .../scratchpad/final/wp2work/fat_head && for c in config config_epep config_emem; do ./ana $E/cfg2/$c.json; done)
    root -l -b -q -e "$B" 'plots/all_report/wp2_scan_prepare_empty.C+O'   # empty/scan2/scan_compact_empty.root
    root -l -b -q -e "$B" 'plots/all_report/wp2_patcheck.C+O'   # wp2/wp2_patcheck.log
    root -l -b -q -e "$B" 'plots/all_report/wp2_upstream.C+O'   # wp2/wp2_upstream.root (15 s)
    root -l -b -q -e "$B" 'plots/all_report/wp2_win.C+O'        # wp2/wp2_win.root, wp2_win.txt, wp2_win_summary.txt (10 s)
    root -l -b -q -e "$B" 'plots/all_report/wp2_vzfit.C+O'      # wp2/wp2_vzfit.root, wp2_vzfit.txt
    root -l -b -q -e "$B" 'plots/all_report/wp2_wincomp.C+O'    # wp2_wincomp.txt
    python3 -I plots/all_report/wp2_tableN1.py                   # wp2_N1.txt, wp2_N1.tex
    root -l -b -q -e "$B" 'plots/all_report/wp2_fig.C+'         # figs_all/wp2_N1v2.pdf/.png, wp2_N1e.pdf/.png
    # remote copy of logs/lists/sources: rsync ... --exclude out/ vae25:.../wp2_empty/ $E/remote_copy/

**Macros** (all in `plots/all_report/`):
- from the first attempt: wp2_lumi.py, wp2_runs.C, wp2_runs.py, wp2_sector.C, wp2_simfs.C, wp2_csec.py, wp2_simvz.C,
  wp2_patcheck.C, wp2_upstream.C;
- edited here: wp2_scan_prepare_empty.C (paths to fat2/scan2), wp2_win.C (rewritten), wp2_fig.C (panels c, d, e);
- new: wp2_vzfit.C, wp2_wincomp.C, wp2_tableN1.py.

**Unused leftovers:** wp2_fig.C no longer calls the non-existent wp2_fig_win.C.

**Large outputs:** `/hdd3/przygoda/HADES/FAT_out/all_report/wp2/`:
- wp2_upstream.root, wp2_win.root, wp2_vzfit.root, wp2_patcheck.log, wp2_upstream.log;
- empty/{pat, pat_parts2, cfg2, fat2, scan2, remote_copy}.
