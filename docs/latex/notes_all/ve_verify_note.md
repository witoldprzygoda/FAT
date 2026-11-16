# ve_verify: independent check of the event-vertex efficiency for the eta (vertex z > -500 mm)

Date 2026-11-14. Role: adversarial verifier of the two reports ve_data and ve_sim_cause. All code is new and my own (it does
not use their macros, and it does not use plots/ScanEngine.h or src/sim_truth.h; it only follows the same written
conventions). No repository file was edited, nothing was deleted, no git command was run.

- Macros and text results: `/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/vtx_eff/ve_verify/`
  (vv_fill.C, vv_lib.h, vv_load.h, vv_ana.C, vv_ana2.C, vv_ana3.C, vv_sim.C, vv_k0s.C; results vv_ana.txt, vv_ana2.txt,
  vv_ana3.txt, vv_sim.txt, vv_k0s.txt)
- Rows, ACLiC build, logs: `/hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_verify/` (rows/d_00..28.root, rows/s_00.root, rows/s_01.root)

The cut itself is mandatory and is not questioned. Only its efficiency for the eta is measured.

---------------------------------------------------------------------------------------------------
## 0. Verdict in short

1. **The numbers of both agents are reproduced exactly with my own code.** S(V) = 4657 +- 390, S(Nv) = 902 +- 401,
   eps_S = 0.838 +- 0.063; T(V) = 3297 +- 110, T(Nv) = -4 +- 77; E(V) = 1241 +- 76, E(Nv) = 21 +- 80; simulation C5 0.99970,
   C4 307 of 307.
2. **The "eta in the no-vertex rows" is not an eta signal.** It sits in the like-sign CB, not in the unlike-sign rows:
   - a fit of the unlike-sign rows alone (no CB) gives S(Nv) = 165 +- 324;
   - the CB of Nv itself has a Gaussian dip of -736 +- 341 rows at m_eta, caused by the e+e+ rows (-256 +- 108, a 2.4 sigma
     downward fluctuation; e-e-: +73 +- 283);
   - with a smooth like-sign CB (the like-sign spectra fitted with a polynomial) S(Nv) = 263 +- 324 (0.8 sigma);
   - the free fit of the Nv "peak" is 3.75 +- 1.0 MeV wide, narrower than the detector resolution (5.4 MeV).
   Both readings of the two agents are therefore too strong: ve_sim_cause calls it real eta from steel (not supported),
   ve_data calls it a fluctuation of S (right, but it is a fluctuation of the CB).
3. **Target eta lose their vertex at the 1 % level at most.** eps_vtx(eta, data) = **0.994 +- 0.021** (smooth-CB S, T and E
   combined, chi2 0.7/2); with the nominal S instead: 0.990 +- 0.021 (chi2 3.3/2). One-sided 95 % CL: > 0.96.
4. **Simulation, same selection:** eps_vtx(C4+C5) = **0.99970 +- 0.00002** (0.99975 in the peak window); bump C2+C3 0.9986 +- 0.0007.
5. **Correction factor:** C_vtx = eps_sim/eps_data = **1.006 +- 0.021** for sigma(pp -> eta X), **1.00 +- 0.02** for the bump
   scale. Recommendation: do **not** apply it; quote it as a one-sided systematic of **+2 %** (+4 % at 95 % CL) on the absolute
   results. 0.43 mb stays 0.43 mb (0.433 mb if applied).
6. To reach 0.93 mb the data efficiency would have to be 0.46. Excluded by T (23 sigma), E (7.7 sigma), smooth-CB S
   (7.6 sigma) and even the nominal S (4.3 sigma). **The vertex cut is not the cause of the deficit.**

---------------------------------------------------------------------------------------------------
## 1. Step 1: rows straight from PAT (own reader, vv_fill.C)

- **Data input:** the 29 day files of `/hdd3/przygoda/HADES/files/exp/gen4/lep_betanew.list`, trees PipPimEpEm_ID (us),
  PipPimEpEp_ID (pp), PipPimEmEm_ID (mm). Rows with isBest >= 0 and trigbit == 8192, any vertex.
- **Simulation input:** `/hdd3/przygoda/HADES/files/sim/ppee_eta_all.root` and `ppee_std_all.root`, tree PipPimEpEm_ID,
  isBest >= 0 and (trigbit & 8192) != 0, with the truth columns of the four tracks.
- **What is computed per row:**
  - M4, M(ee), M(pipi) from the corrected momenta (lepton slot 1 `_p_corr_ep`, slot 2 `_p_corr_em`, pions `_p_corr_pip/pim`),
    masses 0.51099895 and 139.56995 MeV; cor_ok = all four finite;
  - photon: ECAL slots with neutr_pid == 1, 0.8 <= beta <= 1.2, cluster energy > 100 MeV (cluster theta/phi); the three most
    energetic; the one with M(ee gamma) closest to 0.140; m5c = M(pi+pi-e+e-gamma) with the photon energy scaled so that
    M(ee gamma) = m_pi0;
  - track z and r of the four tracks; median z (zp); number of own tracks in the finder window (z in (-170, -30), |r| < 20);
    crossing point of the two pion lines (z, r, doca);
  - event group: adjacent rows with equal (run, event) in one tree of one file.
- **Row counts (sum of 29 files):**

  | sample | rows (PT3) | with vertex | cor_ok | cor_ok with vertex |
  |---|---|---|---|---|
  | us | 7 063 032 | 4 783 467 (0.6773) | 7 003 545 | 4 753 629 (0.6787) |
  | pp | 416 407 | 122 875 (0.2951) | 411 512 | 121 789 (0.2960) |
  | mm | 1 543 953 | 371 524 (0.2406) | 1 523 048 | 367 954 (0.2416) |

  The "with vertex" and "cor_ok with vertex" numbers equal the FAT production exactly (spec_vz500.txt).

## 2. Step 2: S in V and Nv with my own P_main (vv_lib.h, vv_ana.C)

- **Method (written from scratch):** 2.5 MeV bins 0.35-0.95; CB = 2 sqrt(N++ N--) per 10 MeV, linear interpolation to the
  fine bins, CB covariance between fine bins; generalised least squares of N - CB with Gaussian (0.5439, 5.4 MeV) + pol2 in
  0.5225-0.61, variance = expected rows (iterated) + CB covariance. stat_g = stat x sqrt(event-group variance / N) in the
  peak window. Classes: V = vz > -500, Nv = vz <= -500; base PT3, isBest >= 0, cor_ok.
- **Results:**

  | class | N (peak window) | CB | CB/N | S | stat | stat_g | chi2/31 |
  |---|---|---|---|---|---|---|---|
  | all | 152 059 | 50 258 | 0.331 | 5533 | 560 | 572 | 38.5 |
  | V | 85 942 | 11 488 | 0.134 | **4657** | 390 | 397 | 33.5 |
  | Nv | 66 117 | 38 664 | 0.585 | **902** | 401 | 411 | 35.5 |
  | vz > -160 (report "pre") | 77 425 | 7 949 | 0.103 | 4671 | 367 | 372 | 37.2 |
  | -500 < vz <= -160 | 8 517 | 3 523 | 0.414 | -19 | 133 | 138 | 44.9 |

  - "pre" reproduces the report (4671 +- 367). S(V) + S(Nv) = 5560 against S(all) = 5533.
  - eps_S = 4657/(4657 + 902) = **0.838 +- 0.063** (stat_g). isBest == 1: 0.845 +- 0.063.
  - The vz > -160 cut keeps all the eta of V in the data (-19 +- 133 below it); the simulation keeps 0.9973 of C4+C5.
- **Model grid** (order 1-4, upper bound 0.58-0.65, lower bound 0.5225-0.5275, width +-0.5 MeV; 152 pairs with p > 0.05):
  S(V) mean 4572 (rms 478), S(Nv) mean 1249 (rms 498), eps mean 0.789, rms 0.060, range 0.60-0.87.

## 3. Step 3: is the eta in Nv real? (vv_ana.C, vv_ana2.C, vv_ana3.C)

1. **Unlike-sign rows alone.** The same frozen fit on N (no CB; pol2 / pol3) gives S(Nv) = 165 +- 324 / 393 +- 416. In V:
   4646 +- 370 / 4668 +- 476, the same as with the CB.
2. **The CB carries the "signal".** The frozen Gaussian fitted to the CB of Nv gives -736 +- 341 rows (a dip). The source is
   the e+e+ rows: N++ = -256 +- 108 (2.4 sigma), N-- = +73 +- 283, N++ + N-- = -183 +- 303. In V the CB has no dip (-18 +- 182).
   The like-sign rows are lost 3-4 times more often by the vertex cut than the e+e- rows, so Nv has a large like-sign sample
   and its CB is 58 % of N; a 2.4 sigma hole in N++ then moves S(Nv) by several hundred rows.
3. **Smooth CB.** N++ and N-- fitted with pol3-pol5 in 0.45/0.48-0.66/0.70 (fits good: 94.9/95 and 100.5/95), CB = 2 sqrt(f++ f--)
   per fine bin, P_main as before: S(V) = 4670 +- 370, **S(Nv) = 263 +- 324** (rms over 12 variants 30), eps = **0.947 +- 0.063**.
   CB in 20 MeV bins: S(Nv) = 557 +- 365.
4. **Position and width.** Free fit (profile chi2 grid): V 0.54325 +- 0.00037, 5.50 +- 0.50 MeV; Nv 0.54575 +- 0.0014,
   **3.75 +- 1.0 MeV**; Nv period III: width at the grid edge (3.0 MeV). A true eta, from the target or from upstream material,
   cannot be narrower than the resolution.
5. **Is a 2.2 sigma value unusual in Nv?** Frozen-width fits at 39 positions 0.575-0.865: pulls in Nv have mean 0.05, rms 0.93,
   2 of 39 beyond 2 sigma. The errors of Nv are right; at the eta position the pull is 2.25 (one-sided p about 1 %), and that
   pull comes from the CB (item 2).
6. **CB scale.** k = 0.85 ... 1.15 moves S(Nv) from 793 to 1011 (S(V) unchanged). N/CB in the Nv peak window is 0.96 for
   0.03 < M(ee) < 0.10 and 1.05 for M(ee) > 0.10, so the scale is right within about 5 %. A fit with a free CB scale is
   degenerate with pol2 (k = 2.6 +- 1.6 in V, -2.0 +- 1.3 in Nv) and not informative.
7. **Other physics content.** Nv has none of the target structures:
   - ECAL-tagged eta: T(V) = 3297 +- 110, **T(Nv) = -4 +- 77** -> eps_T = **1.001 +- 0.023** (ve_data: 27 +- 67, 0.992 +- 0.020).
   - Bump edge: E(V) = 1241 +- 76, **E(Nv) = 21 +- 80** -> eps_E = **0.983 +- 0.062**; over the 18 model variations applied to both
     classes eps_E = 0.965-1.043, rms 0.028.
   - The M4 structure at 0.73-0.75 (omega region; pulls up to 12.5 in V): S_V = 7483 +- 600, S_Nv = -914 +- 448 at 0.740;
     ratio Nv/V = -0.12 +- 0.06 (0.7325: -0.04 +- 0.11).
   - If the Nv eta were real (from the target or from material), T and E would be there in proportion; a material source,
     richer in conversions, would make E even larger relative to S.
8. **Does the photon tag work without a vertex?** Yes. pi0 in M(ee gamma) (M(ee) < 0.10, N - CB): V 539 887 +- 2342 at
   0.1416, sigma 13.5 MeV; Nv 61 648 +- 4364 at 0.1361, sigma 18.4 MeV. The Nv peak is 5.5 MeV lower and 36 % wider (the
   ECAL direction and beta use the cluster vertex when there is no reconstructed vertex, Hydra hemcneutralcandfinder.cc), but
   most of it is inside the tag window 0.119-0.173. For a target event without a vertex the cluster vertex is at the right
   place, so T is a valid test of the loss of target eta. Tag rate per (N - CB): V 0.165, Nv 0.080.

## 4. Step 4: where do the Nv rows come from? (checks of ve_sim_cause)

- **The zp argument of ve_sim_cause is circular.** A row has no vertex when at most one of its tracks is in the finder
  window, so the median track z of an Nv row lies outside the window almost by construction: only 108 124 of 2 249 916 Nv rows
  (4.8 %) have zp in (-170, -30), and Nv rows with >= 2 own tracks in the window are 1407. "The Nv eta sits at zp < -170"
  (875 +- 392) would hold for target events with lost vertices too. It is not evidence for an upstream source.
- **Non-circular evidence that the Nv rows are mostly not LH2 events** (I agree with this part): pion crossing points with
  doca < 10 mm lie at r >= 8 mm in 87 % of the Nv rows against 35 % of V rows (r is not used by the finder below 20 mm);
  ve_sim_cause: elastic pp per 2-track event upstream 0.028 +- 0.003 of the cell's, the Nv share follows the beam period,
  the median z does not change with polar angle; CB/N 0.585 against 0.134. This agrees with the material study vw2b
  (front rows off axis at r >= 8 mm, same composition in empty-target runs, no measurable eta in front of the cell).
- **Per period** (S(N - CB) per 1e3 us rows): V 1.07, 0.91, 0.89, 0.91 (I, IIA, IIB, III); Nv 0.02 +- 0.26, 0.89 +- 0.61,
  0.37 +- 0.37, 1.03 +- 0.38. The period-III hint (515 +- 191) is half CB: S(N) = 254 +- 154, CB dip -273 +- 161.
- **By pion crossing:** eps_S = 1.001 +- 0.037 (r < 5 mm), 1.044 +- 0.099 (5-8 mm), 0.738 +- 0.121 (r >= 8 mm), 0.690 +- 0.135
  (doca >= 10 mm); the nominal-S "loss" sits only where the Nv CB is.

## 5. Step 5: simulation, same selection (vv_sim.C, own truth tagging)

- Truth (GEANT ids e+ 2, e- 3, pi+ 8, pi- 9, pi0 7, eta 17; process 5 decay, 6 pair production):
  C5 = the four tracks are direct products of the forced eta decay (geninfo 1702030848.f, parent -1, process 0, same geninfo2);
  C4 = pions from eta (same vertex), leptons from one conversion with grandparent eta, equal generator codes;
  C2 = leptons from pi0 Dalitz with grandparent eta; C3 = leptons from one conversion with grandparent pi0.
- Base: isBest >= 0, PT3 bit, cor_ok (no ownership needed: C5 rows belong to the eta set, C2-C4 rows to the standard set).

  | class | rows | eps(vz > -500) | lost | eps(vz > -160) |
  |---|---|---|---|---|
  | C5, all M4 | 572 822 | 0.99970 +- 0.00002 | 171 | 0.99725 |
  | C5, peak window 0.5333-0.5587 | 484 140 | 0.99975 | 120 | 0.99774 |
  | C4, all M4 | 307 | 1 (> 0.990 at 95 % CL) | 0 | 1 |
  | C4+C5, all M4 | 573 129 | **0.99970 +- 0.00002** | 171 | 0.99725 |
  | C2 | 1 660 | 0.99940 +- 0.00060 | 1 | 0.99337 |
  | C3 | 1 192 | 0.99748 +- 0.00145 | 3 | 0.99581 |
  | C2+C3, bump window | 2 152 | **0.99814 +- 0.00093** | 4 | 0.99396 |

  This is the same as ve_sim_cause (572 822 / 171 lost; 307/307; 1660/1; 1192/3). The statement "the simulation keeps >= 0.96"
  is true but loose: it keeps 0.9997.

## 6. Step 6: combination and the correction factor

| measure of eps_vtx(eta, data) | value |
|---|---|
| S, P_main, nominal CB | 0.838 +- 0.063 (stat_g) +- 0.060 (model) |
| S, unlike-sign rows only (pol2) | 0.966 +- 0.065 |
| S, smooth like-sign CB | 0.947 +- 0.063 |
| T, ECAL-tagged eta | 1.001 +- 0.023 |
| E, bump edge step | 0.983 +- 0.062 +- 0.028 |
| **S (smooth CB) + T + E** | **0.994 +- 0.021**, chi2 0.67/2, > 0.96 at 95 % CL |
| S (nominal) + T + E | 0.990 +- 0.021, chi2 3.3/2 |
| S (smooth CB) + E, no ECAL | 0.964 +- 0.046 |

- **sigma(pp -> eta X):** C_vtx = 0.9997/0.994 = **1.006 +- 0.021** (upper limit 1.04 at 95 % CL). It is compatible with 1;
  physically it cannot be below 0.9997 (eps_data <= 1). **Do not apply it; quote +2 % (one-sided) as a
  systematic** of the absolute result. Applied, 0.43 mb would become 0.433 mb.
- **Bump scale:** T measures the same decay (eta -> pi+pi-pi0, pi0 Dalitz). T and E combined: 0.999 +- 0.022; simulation
  0.9986 (C2+C3, bump window). C_vtx(bump) = **1.00 +- 0.02**; same treatment.
- **What would be needed for the published value:** 0.43 mb x 0.9997/eps = 0.93 mb needs eps = 0.46 (1.14 mb: 0.38).
  Excluded by T at 23 sigma, E at 7.7 sigma, smooth-CB S at 7.6 sigma, nominal S at 4.3 sigma.
- **The worst single number,** nominal S alone (0.838), would give 0.51 mb, not 0.49 mb as written by ve_data (ve_data used
  eps_sim = 0.96 instead of 0.9997 in "sigma scales as eps_sim/eps_data"; its 0.42 mb should read 0.436 mb). The nominal-S
  value is a CB fluctuation and should not be used.

## 7. Findings on the two reports

- **ve_data:** numbers confirmed. Conclusion (eps about 0.99, cut not the cause) confirmed. Two corrections: (a) the S(Nv) hint
  is a hole in the e+e+ CB, not an excess of e+e- rows (its caveat 1 should say this); (b) the sigma arithmetic uses
  eps_sim = 0.96 instead of 0.9997 (0.42 -> 0.436 mb; 0.49 -> 0.51 mb). The combination of S with T and E is acceptable
  because T dominates.
- **ve_sim_cause:** simulation numbers and the vertex-finder description confirmed. Its claim that the Nv eta (906 +- 401) is
  real eta made in the steel is not supported: the zp split is circular (Section 4), the excess is in the CB, the width is
  below resolution, and T, E and the omega-region structure are all absent in Nv. Its final recommendation (keep the
  simulated eps; correction about +1 %, at most +3.5 %) agrees with mine.
- **K0S (side check, vv_k0s.C):** M(pi+pi-) of rows with M4 < 1.10, frozen Gaussian at 0.4950 + pol2/pol3: Y_V = 13 464 +- 715,
  Y_Nv = 2603 +- 551 (same width), eps 0.82-0.85 over ranges and orders; 0.77 if the Nv peak is 8.9 MeV wide (ve_sim_cause
  0.718 uses another method). Unlike the eta, the K0S is visible in Nv. Even 0.72 explains only a factor 1.4 of the K0S
  deficit (0.29 x SMASH). Not studied further.

## 8. Caveats

1. T assumes the photon tag works for a target event without a vertex. The pi0 peak exists in Nv (shifted, wider) and for a
   target event the fallback cluster vertex is at the right place, so this is well founded but not proven with target events.
2. E(V) is 16 sigma statistically but has a large model spread (+-319, linear sum over slices, ve_data); with the
   statistical error alone its limit is eps_E > 0.89 (95 % CL).
3. The smooth-CB variant is not the report convention; it is used here only to remove the statistical structure of the
   like-sign rows, which cannot contain an eta. With the nominal CB the combination changes from 0.994 to 0.990.
4. A strongly smeared eta component in Nv (target events with badly reconstructed tracks) would escape frozen-width fits.
   N-only fits with 8 and 12 MeV Gaussians give S(Nv) = 70 +- 471 and -150 +- 812: no sign of it.
5. Not checked: which vertex condition the elastic luminosity analysis used. If it differs from vz > -500, the luminosity and
   the yield could differ by about 1 % (the elastic loss estimate of ve_sim_cause), within the +2 % above.

## 9. Exact commands

```
D=/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/vtx_eff/ve_verify
W=/hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_verify
mkdir -p $W/rows $W/aclic $W/logs; cd $W
B='gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_verify/aclic", kTRUE)'
seq 0 28 | xargs -P 6 -I{} sh -c "root -l -b -q -e '$B' '$D/vv_fill.C+({},0)' > logs/fill_d{}.log 2>&1"    # data, 14 s
for i in 0 1; do root -l -b -q -e "$B" "$D/vv_fill.C+($i,1)" > logs/fill_s$i.log 2>&1 & done; wait       # simulation
cd $D
root -l -b -q -e "$B" 'vv_ana.C+'  > $W/logs/ana.log  2>&1    # -> vv_ana.txt  (S, eps, grid, free fits, null tests, T, E)
root -l -b -q -e "$B" 'vv_ana2.C+' > $W/logs/ana2.log 2>&1    # -> vv_ana2.txt (N-only fits, CB dip, pi0 tag, omega region)
root -l -b -q -e "$B" 'vv_ana3.C+' > $W/logs/ana3.log 2>&1    # -> vv_ana3.txt (smooth CB, 20 MeV CB, wide Gaussians)
root -l -b -q -e "$B" 'vv_sim.C+'                             # -> vv_sim.txt  (simulation eps per truth class)
root -l -b -q -e "$B" 'vv_k0s.C+'  > $W/logs/k0s.log  2>&1    # -> vv_k0s.txt  (K0S side check)
```
The CPU package stayed at 40-64 C; at most 6 processes ran at once.
