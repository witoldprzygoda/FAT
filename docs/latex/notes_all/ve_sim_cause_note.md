# ve_sim_cause: efficiency of the event-vertex cut (vertex z > -500 mm), simulation and cause

Date 2026-11-14. Package prefix `ve_sim_cause`. HADES p+p, T = 4.5 GeV (feb22), channel pi+pi-e+e-.
The vertex cut itself is not questioned. This note only measures what it does and explains why data and simulation differ.

Work directory (macros, text outputs): `/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/vtx_eff/ve_sim_cause/` (called `W/` below).
Big outputs (ROOT files, logs, ACLiC build): `/hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_sim_cause/` (called `O/` below).
No repository file was edited, nothing was deleted, no git command was run.

---------------------------------------------------------------------------------------------------
## 0. Short answer

1. **Simulation.** The simulation keeps essentially every signal row: eps_vtx = 0.9997 for the direct eta -> pi+pi-e+e-
   (C5), 307 of 307 rows for the conversion peak (C4), 0.9994 / 0.9975 for the bump classes C2 / C3, and 0.991-0.996 for
   K0S pion pairs.
2. **Data.**
   - The cut keeps 0.68 of the unlike-sign rows and 0.84 +- 0.06 of the eta signal (P_main).
   - The eta signal of the rows without a vertex, 906 +- 401, does NOT come from the LH2 cell. The median z of its four
     tracks lies upstream of the cell, at z < -170 mm (875 +- 392). The part with the median z inside the cell or next to
     it is 28 +- 83.
3. **Cause.**
   - What PAT stores is the DST vertex of `HParticleVertexFind`. It is set to -1000 when fewer than two tracks point into
     z (-170, -30) mm with r < 20 mm. The finder uses no START and no timing information.
   - In the data, 35 % of the rows have tracks that start upstream of the target cell, near z = -205 mm. This is where
     the steel and copper target holder and the steel beam tube sit (HGeant target geometry). Hardly any of these events
     are hydrogen: their pp-elastic content is 2.8 +- 0.3 % of that of the LH2 cell.
   - The simulation generates vertices only inside the LH2 cell and has no beam halo, so it has no such events.
   - The data/simulation difference is therefore non-target background, and the cut removes it correctly. It is not an
     inefficiency for target events.
4. **Real loss of target events in the data:** about 0.8 % (simulation 0.03 %), and at most 3.5 % (95 % CL) from the
   signal itself.
5. **Consequence for sigma(pp -> eta X):**
   - Using the simulated eps_vtx is correct to about 1 %.
   - A correction of at most +1 % (+3.5 % at 95 % CL) applies, i.e. 0.43 -> 0.435 mb.
   - The vertex cut does NOT explain the factor 0.41 x SMASH, nor the 3-sigma gap to 0.93-1.14 mb.
   - Dividing by the data value 0.84 would be wrong: it would count eta mesons produced in the steel.

---------------------------------------------------------------------------------------------------
## 1. How the event vertex is made, and when it is -1000

1. **PAT.**
   - In `core/hparticlepool.cc` (lines 264-270 of the pat3 production tag `prod-feb22-betanew-2026-10-28`; copy in
     `F/src/tag/core/`), PAT copies the vertex from the event header: `gHades->getCurrentEvent()->getHeader()->getVertexReco()`
     gives `eVertReco_x/y/z/chi2`.
   - PAT does not refit the vertex and never writes -1000 itself. The value comes from the DST.
2. **DST production, data and simulation alike.**
   - Data: `/lustre/hades/user/hadesdst/sub/feb22/dstreal/analysisDST.cc`. Simulation:
     `/lustre/hades/user/hadesdst/sub/feb22/sim/p4500p/gen4/smash/smash_standard/dst/analysisDST.cc`. Copies are in
     `W/src_dst/`.
   - Both run `HParticleTrackCleaner` (sets kIsUsed with the default `HParticleTrackSorter`), then
     `HParticleVertexFind("particlevertexfind","particlevertexfind","")` with `setCut(-170.,-30.,20.)`.
   - Both DST builds use Hydra 2.7.4 (`par1=.../hydra2-7.4/defall.sh`).
   - `particle/hparticlevertexfind.cc` is byte-identical in Hydra 2.7.3, 2.7.4 and 2.7.6 (md5 8c726694c549e09ddf099ffd32855bf6,
     checked on vae25 under /cvmfs/hadessoft.gsi.de/source). Copy in `W/src_hydra76/`.
3. **The algorithm (`readInput` and `doFit`).**
   - **Candidates:** `HParticleCand` with kIsUsed, not fake-rejected, and track z in (-170, -30) mm and r < 20 mm. Here z and
     r are the point of closest approach of the track to the beam axis, the same values PAT stores as `<x>_z`, `<x>_r`.
   - **More than 10 candidates (rare in p+p):** an extra filter on RK chi2 < 100, inner-segment chi2 < 6 and DST beta between
     0 and 1.2, followed by an iterative window around the mean z (z' relative to the x, y of the MDC cluster vertex).
     This is the only place where a timing quantity (the START-based DST beta) enters.
   - **Lines:** each surviving track is turned into a line in the lab. It is kept if its closest approach to the z axis
     has z in (-170, -30) and rho < 20 mm.
   - **Result:** with two or more lines, a least-squares point of closest approach is computed, followed by Tukey-weighted
     iterations (at most 100). With one line or none, the vertex is (-1000, -1000, -1000), chi2 = -1 and iterations = -1.
     Tracks with Tukey weight > 0 get kIsUsedVertex (PAT column `isUsedVertex`).
   - There is no requirement on the START hit, on the number of tracks beyond two, or on fit quality. Not converging after
     100 iterations only sets chi2 = -1; the vertex is still written.
4. **Check with the rows themselves.** nwin = number of the four row tracks with z in (-170, -30) and |r| < 20 mm.

   | nwin | data us: f(no vertex) (rows) | simulation std: f(no vertex) |
   |---|---|---|
   | 4 | 0.0000 (2 866 284) | 0.0000 |
   | 3 | 0.0000 (1 168 435) | 0.0000 |
   | 2 | 0.0028 (513 888) | 0.0010 |
   | 1 | 0.72 (724 402) | 0.38 |
   | 0 | 0.98 (1 790 023) | 0.89 |

   - The rule is the same in data and simulation: two tracks in the window always give a vertex.
   - The difference is how many rows have at most one track in the window: data 0.356 (us), 0.78 (pp), 0.81 (mm);
     simulation 0.003.

---------------------------------------------------------------------------------------------------
## 2. Simulation: eps_vtx per truth class

How it was measured:
- **Files.** The simulated PAT files `/hdd3/przygoda/HADES/files/sim/ppee_{std,lep,eta}_all.root`, tree `PipPimEpEm_ID`.
  These are the same files FAT reads (sim_*.list of FAT_simwt).
- **Base selection.** isBest >= 0 and the PT3 bit `(trigbit & 8192) != 0`. In the mass windows also cor_ok (all four
  corrected momenta finite). M4 is computed from the corrected momenta exactly as FAT does.
- **Truth.** Classes, weights W (rules v2, PDG BR) and ownership (mode A) come from `FAT_simwt/src/sim_truth.h`
  (md5 4e82009176ff14ab72fdcba92966796e), called exactly as in FAT `main.cc` (`EvaluateRow(set, trk, 2, 3, opt)`). Only
  owned rows are counted.
- **eps_vtx** = owned rows with vertex z > -500 divided by all owned rows.
- **Errors.** Event groups (adjacent rows with equal run and event), ratio estimator.
- **Windows.** Simulation-centred peak 0.5333-0.5587 and bump 0.4576-0.5255. The data windows applied to 0.99615 x M4 give the
  same numbers to 1e-4. K0S: M(pi+pi-) 0.485-0.512 GeV/c2.
- No physics re-weighting (V3) was applied. It changes weights inside a class only, and eps is 0.9997 either way
  (weighted and unweighted agree).

| class (set) | all M4 | peak window | bump window | rows (all M4) | vz > -160 given vertex |
|---|---|---|---|---|---|
| C5 eta -> pi+pi-e+e- (eta) | 0.99970 +- 0.00002 | 0.99975 +- 0.00002 | 0.9994 +- 0.0001 | 572 822 | 0.9975 |
| C4 eta -> pi+pi-gamma, conversion (std) | 307/307 = 1 (> 0.990 at 95 % CL) | 252/252 | 16/16 | 307 | 1 |
| C2 bump, pi0 Dalitz (std) | 0.9994 +- 0.0006 | - | 0.9992 +- 0.0008 | 1 660 | 0.994 |
| C3 bump, conversion (std) | 0.9975 +- 0.0015 | - | 0.9966 +- 0.0019 | 1 192 | 0.998 |
| K0S pion pair P6 (std + lep), unweighted | 0.9956 +- 0.0003 | - | - | 51 417 | 0.996 |
| K0S P6, weighted W | 0.9910 +- 0.0020 | K0S window 0.9915 +- 0.0021 | | | |
| T1+T2, weighted | 0.9989 +- 0.0001 | 0.9983 +- 0.0007 | 0.9964 +- 0.0008 | 9.25 M | 0.9987 |

- The claim "the simulation keeps >= 0.96" is true: it keeps 0.9997.
- Lost simulated C5 rows: 171 of 572 822.
- Table: `W/vtxana_tables.txt`, section A.

---------------------------------------------------------------------------------------------------
## 3. Data: rows and signals with and without the vertex

**Rows.** The 29 PAT day files of `lep_betanew.list`, trees us / pp / mm, isBest >= 0, trigbit == 8192.
- Rows that pass: us 4 783 467 of 7 063 032 (0.6773), pp 122 875 of 416 407 (0.2951), mm 371 524 of 1 543 953 (0.2406).
  The pass counts are identical to the FAT row counts in `spec_vz500.txt`.
- cor_ok rows: us 0.6787 +- 0.0002; bump window 0.5705 +- 0.0010; peak window 0.5651 +- 0.0013 (86 015 of 152 209 rows).

**Eta peak.**
- Method: P_main, the frozen Gaussian 0.5439 GeV/c2 / 5.4 MeV/c2 + pol2 in 0.5225-0.61, fitted to N - CB. CB = 2 sqrt(N++ N--)
  per 10 MeV of the same category of rows (k = 1). The code is the `scan::FitPeakMain` of `plots/ScanEngine.h`, fed with
  hand-made spectra in the ScanEngine binning (2.5 MeV, CB interpolated as in FillMany).

| category | S | N in peak window | CB/N |
|---|---|---|---|
| all PAT rows | 5537 +- 560 | 152 060 | 0.330 |
| vertex (vz > -500) | 4658 +- 390 | 85 943 | 0.134 |
| no vertex | 906 +- 401 | 66 117 | 0.585 |
| vertex, vz > -160 | 4672 +- 367 | 77 426 | 0.103 |
| vertex, vz <= -160 | -19 +- 133 | 8 518 | 0.414 |
| no vertex, evtStart 1 / 0 | 917 +- 366 / -11 +- 164 | | 0.585 / 0.582 |
| vertex, evtStart 1 / 0 | 4356 +- 373 / 298 +- 115 | | |

- eps_S(data) = S(vertex) / (S(vertex) + S(no vertex)) = 0.837 +- 0.061.
- Width check (`W/vtxwidth.txt`): the no-vertex S is 975 +- 475 (6.5 MeV) and 1024 +- 550 (7.5 MeV). The K0S peak is 1.4 times
  wider in the rows without a vertex, so a wider eta peak there is likely. With the wider peak, eps_S would be 0.82.

**K0S -> pi+pi-** (`W/kpass_fit.txt`).
- Method: free Gaussian + pol2 fit in 0.46-0.54 of M(pi+pi-), cor_ok us rows.
- Yields: with vertex 15 893 +- 881 (sigma 6.4 MeV); without vertex 6 235 +- 819 (sigma 8.9 MeV). eps_vtx(K0S, data) =
  0.718 +- 0.029; simulation 0.991-0.996.
- 5 343 +- 857 of the no-vertex K0S have zp < -170 mm (zp is defined in section 4). For K0S, zp is biased towards the cell
  because the pions are displaced downstream.

---------------------------------------------------------------------------------------------------
## 4. Where the rows without a vertex come from

Every number below is from the full data pass, unless a single day is named.

1. **The tracks start upstream.**
   - Day 049: 87.5 % of the tracks of the no-vertex rows have z <= -170 mm (simulation 9 %). The median is -209 mm,
     quantiles 2 % / 98 %: -545 / -112. The tracks of the rows with a vertex have median -130 mm.
     Files: `W/zd_data049.txt`, `W/ex_data049.txt`.
   - Pseudo-vertex zp = median of the four track z of the row. In the data, zp shows two populations
     (`W/zpprint.C`, output `W/zpprint.txt`):
     - the LH2 cell, zp -160 to -90 mm (peak 753 241 rows per 10 mm at -140);
     - a second peak at zp = -200 to -180 mm (384 000 rows per 10 mm), with a tail down to -700 mm.
   - Shares with zp < -170 (`W/zpfrac.txt`): data us rows 0.354, of which no-vertex rows 0.949 and vertex rows 0.070;
     pp 0.776; mm 0.821. Simulation 0.0001-0.0002.
2. **It is a fixed place, not a reconstruction distortion.** The median track z of the no-vertex rows is the same at all
   polar angles: -214, -210, -208, -206, -205, -204, -204, -203, -197 mm for theta 18-24 ... 66-72 deg. The vertex rows
   give -128 to -135 (`W/zth.C`, output `W/zth.txt`). A timing or drift-time error would change with angle.
3. **The place is the target holder and the beam tube.**
   - HGeant target geometry TARGET v87 (Oracle 04-FEB-2026): `/lustre/hades/user/hadesdst/sub/feb22/p4500p/hgeant/p4500p/smash_standard/hgeant/geom/target.geo`,
     copy `W/src_dst/hgeant_geom_target.txt`.
   - Lab z = RTAM z - 532.2 mm. RICH sits at -532.2 in CAVE and RTAM at 0 in RICH.
   - What is where along the beam:

     | lab z [mm] | element |
     |---|---|
     | -148.4 to -98.6 | LH2 cell (cone TARG); the simulated true vertices span -147.5 to -99.5 |
     | -174 to -148 | LH2 ring, r 8-12 mm (THB2) |
     | -229 to -174 | steel holder (TFM1), copper (TFO1) and an LH2 ring (THM1), r > 8 mm |
     | -243 to -148 | steel beam tube, inner radius 7.5 mm (TBS4) |
     | -252 to -229 | steel flanges (TVS1, TVS2) |
     | down to -521 | more steel tubes (TBS6, TBS7) |

   - The upstream population (z -260 to -175, peak near -200) lies on the holder and the flanges. Beam halo
     (r > 7.5 mm) hitting steel and copper is the natural source.
   - The simulation has none of this: its HGeant card uses `JVER 2 2 0` (random vertex in the target segments only), and
     start.geo is commented out.
4. **It is not hydrogen.**
   - Test: pp elastic in 2-track events (totalmult == 2, isBest == 1) of the hadron PAT tree PP_ID. Files
     `/hdd2/przygoda/HADES/pp45/exp/had_{036,049,062}_p1.root`, all rows. Signal: tan theta1 tan theta2 within
     0.2943 +- 0.04 and |dphi - 180| < 3 deg, side bands subtracted (`W/elastic.C`, `W/elsum.C`).
   - Results: elastic excess per 2-track event is 0.0198 in the LH2 cell (53 966 +- 324 of 2 723 368 events with a vertex)
     and 0.0006 upstream (593 +- 72 of 1 054 027 no-vertex events with zbar < -170). The ratio is 0.028 +- 0.003.
   - Per file: 0.042 +- 0.007 (036), 0.023 +- 0.005 (049), 0.018 +- 0.008 (062).
   - So 97 % of the upstream events are not on hydrogen. In 37-46 % of all 2-track hadron events there is no vertex.
5. **They are different events.**
   - CB/N in the peak window is 0.585 without vertex and 0.134 with vertex. The cut removes 70-76 % of the like-sign rows.
   - The no-vertex fraction changes with the beam period (us rows): days 034-045 0.31, 046-053 0.37, 059-067 0.275. This is a
     beam condition (halo), not a detector state. The sectors are flat: 0.317-0.327 for every sector and every track.
6. **Other variables** (`W/vtxana_tables.txt`, section C; data us; the simulation is 0.001-0.003 in every bin).
   - **Track multiplicity:** 0.28 / 0.30 / 0.41 / 0.41 / 0.32 / 0.22 / 0.14 for 4-10 tracks. With many tracks, one more
     track lands in the window.
   - **Polar angle:** the fraction rises from 0.25 (17 deg) to 0.39 (60 deg).
   - **Momentum:** lepton 0.36 (p < 200 MeV/c) -> 0.20 (1 GeV/c); pion 0.37 -> 0.28.
   - **Inner-segment chi2 of pi+:** 0.31 (< 3), 0.41 (3-6), 0.61 (6-9). 81 % of the no-vertex rows have chi2 < 3, so
     segment quality is not the cause.
   - **START:** no-vertex fraction 0.53 for evtStart = 0 and 0.30 for evtStart = 1. evtStart = 0 makes up 17 % of the
     no-vertex rows and 7 % of the vertex rows. The simulation always has evtStart = 1. Halo particles that miss the START
     fit this pattern. There is no eta signal in the no-vertex, evtStart = 0 rows (-11 +- 164).
   - **DST beta (START-based):** broader for the no-vertex rows; 49 % of e+ have 0.96 < beta < 1.04, against 69 % with a
     vertex. This follows from the origin: the path and time are computed for the wrong start point. It is not a cause,
     because the vertex finder does not use beta below 11 tracks.
7. **The eta signal without a vertex sits upstream.**
   - Signal of the no-vertex rows per zp bin (`W/vtxana_tables.txt` B2):

     | zp bin [mm] | S |
     |---|---|
     | < -400 | 209 +- 131 |
     | -400 to -300 | 224 +- 116 |
     | -300 to -250 | 64 +- 108 |
     | -250 to -200 | 122 +- 250 |
     | -200 to -180 | 181 +- 198 |
     | -180 to -170 | 103 +- 99 |
     | -170 to -160 | -1 +- 57 |
     | -160 to -150 | 33 +- 37 |
     | -150 to -100 | 3 +- 42 |
     | -100 to -50 | 11 +- 17 |
     | > -50 | -18 +- 18 |

   - One fit per side: zp < -170 gives 875 +- 392; zp >= -170 gives 28 +- 83.
   - For comparison, the vertex rows have S = 3972 +- 304 at zp -150 to -100 and only 77 +- 129 at zp -200 to -170.
   - Signal per non-CB row at zp -250 to -170 is 0.012-0.04, against 0.076 in the cell. These rows are not displaced target
     events. If they were, they would carry the target's signal fraction.

---------------------------------------------------------------------------------------------------
## 5. Real loss of target events in the data

1. **Signal-based limit.** No-vertex eta with zp >= -170: 28 +- 83 = 0.006 +- 0.018 of S(vertex). The one-sided 95 % limit is
   0.035.
2. **Resolution-based estimate** (`W/tgtres.C`, `W/tgtres.txt`). Rows with a vertex inside the cell; data days 036, 049, 062;
   simulation eta + std.
   - Per-track probability to fall outside the window:

     | vz slice [mm] | data | simulation |
     |---|---|---|
     | -150 to -140 | 0.160 | 0.080 |
     | -140 to -125 | 0.113 | 0.067 |
     | -125 to -110 | 0.087 | 0.057 |
     | -110 to -95 | 0.080 | 0.053 |

   - Track z - vz in the cell: rms 23.7 mm (data) against 17.9 mm (simulation). The data beam axis is also off by
     (2.7, -2.3) mm, while the finder measures r from (0, 0).
   - Share of target rows with at most one own track in the window: 0.0219 (data), 0.0009 (simulation).
   - In the simulation, a row with at most one track in the window loses the vertex with probability 0.33-0.38 (other
     event tracks rescue the rest). This reproduces the measured simulated loss: 0.0009 x 0.35 = 0.0003.
   - Applied to the data: 0.022 x 0.35 = about 0.008. The first slice still contains upstream events, so this is an upper
     estimate.
3. **Recommended value:** eps_vtx(target eta, data) = 0.992 +- 0.010, against the simulation's 0.9997.

---------------------------------------------------------------------------------------------------
## 6. What it means for sigma(pp -> eta X) = 0.43 mb

- The result divides by the simulated eps_vtx (about 1). For target events this is right to about 1 %. The correction is
  x 1.008 +- 0.010 (at most x 1.035 at 95 % CL), so 0.43 -> 0.435 mb.
- It must NOT be replaced by the data ratio 0.84 (or K0S 0.72). That 16 % of the PAT eta rows (and 28 % of K0S) are eta and
  K0S made in p + steel/copper upstream, and they are not part of L (LH2 only).
- The deficit (0.41 x SMASH; about 3 sigma below 0.93-1.14 mb at 3.5 GeV; K0S 0.29 x SMASH) has to be looked for elsewhere.
  The vertex cut accounts for at most 1-3.5 %.

---------------------------------------------------------------------------------------------------
## 7. Caveats

1. **zp is crude.** It is the median of four track z values. It is clean for the eta rows (prompt tracks), but biased for
   K0S rows (displaced pions). The upstream share of the no-vertex K0S is therefore a lower limit.
2. **CB of the no-vertex rows.** The no-vertex eta signal uses CB from the no-vertex like-sign rows with k = 1. The p + A
   events may need another k. S without vertex is only 2.3 sigma.
3. **Hydrogen test.** It uses the June hadron production (same DST vertex), 3 of 29 days and 2-track events only, with
   side-band subtraction and no fit. Upstream acceptance for elastic pairs was not modelled.
4. **Resolution-based estimate.** The 0.8 % assumes that the simulated conditional loss (0.33-0.38) holds in the data.
5. **Upstream events with a vertex.** Some upstream events keep a vertex: 7 % of the us rows and 33-35 % of the like-sign
   rows with a vertex have zp < -170. vz > -160 removes part of them, and their eta content is consistent with zero
   (77 +- 129 at zp -200 to -170). The extra CB they bring is in the like-sign rows and therefore in the CB. Not studied
   further.
6. **Luminosity.** L must refer to the same LH2 cell. Upstream elastic content is only 2.8 % of the cell's, so L is not
   strongly affected whatever vertex condition it used. Not checked in the L analysis itself.
7. **What was not checked.** Hydra 2.7.4 sources of `HParticleTrackCleaner` and `HParticleTrackSorter` were not compared with
   2.7.6; the vertex finder is identical. The DST-level kIsUsed was not checked directly. PAT redoes the same sorting, and
   rows with exactly 2 own tracks in the window have a vertex in 99.7 % (data) / 99.9 % (simulation), and rows with 3-4 in
   100 %.

---------------------------------------------------------------------------------------------------
## 8. Commands (run from W/ unless a path is given)

`A='gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_sim_cause/aclic",kTRUE)'`

1. Quick look, one data day and the std simulation:
   - `root -l -b -q -e "$A" 'explore.C+("/hdd3/przygoda/HADES/files/exp/gen4/lep_049.root","PipPimEpEm_ID","ep","em",0)' > ex_data049.txt`
   - The same with the std file and `1` gives `ex_simstd.txt`.
   - `zdist.C+` with the same arguments plus an output file gives `zd_*.txt`.
2. Main pass (simulation per set, data in 4 parts, 6 processes at most):
   - `for s in 0 1 2; do root -l -b -q -e "$A" "vtxpass.C+(\"sim\",$s)" > O/sim_part$s.log 2>&1; done`
   - `for k in 0 1 2 3; do root -l -b -q -e "$A" "vtxpass.C+(\"data\",$k,4)" > O/data_part$k.log 2>&1 & done; wait`
3. Tables: `root -l -b -q -e "$A" 'vtxana.C+'` writes `vtxana_tables.txt`. Then:
   - `root -l -b -q -e "$A" 'vtxwidth.C+'` (width check);
   - `root -l -b -q zpprint.C`, `zpfrac.C`, `zth.C`, `counts.C`.
4. K0S:
   - `for k in 0 1 2; do root -l -b -q -e "$A" "kpass.C+(\"fill\",$k,3)"; done`
   - `root -l -b -q -e "$A" 'kpass.C+("fit")'`
5. Hydrogen test:
   - `root -l -b -q -e "$A" "elastic.C+(\"/hdd2/przygoda/HADES/pp45/exp/had_049_p1.root\",-1,\"O/elastic_full_049_p1.root\")"`
     (also 036_p1 and 062_p1);
   - `root -l -b -q "elsum.C(\"O/elastic_full_036_p1.root,O/elastic_full_049_p1.root,O/elastic_full_062_p1.root\")"`
6. Target resolution:
   - `root -l -b -q -e "$A" 'tgtres.C+("<lep_036,lep_049,lep_062 paths>",0)'`
   - then with `"ppee_eta_all.root,ppee_std_all.root",1`.
7. Sources from vae25: `ssh -o RemoteCommand=none -o RequestTTY=no vae25 'cat <path>'`. The paths are in section 1 and 4.3;
   the copies are in `W/src_hydra76/` and `W/src_dst/`.

Run times: data pass 25 s per part, lep set 104 s, elastic 1-6 min per file. The CPU package stayed at 39-67 C.
