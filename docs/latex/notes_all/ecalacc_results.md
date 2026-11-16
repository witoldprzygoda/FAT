# ecalacc: can the eta bump be reconstructed in full with the ECAL photon? (prefix ecalacc) -- results note

Date 2026-11-15. Branch pp45_pippimepem (worktree /home/przygoda/HADES/FAT/FAT). No git command was run. No report .tex and no
existing file was edited; all outputs are new files with the prefix ecalacc. vae25 was only read (two ROOT macros piped through
ssh, output to stdout; nothing written there). Machine: at most 1 heavy process at a time, CPU package 39-44 C.

**Machine-readable outputs:** `/hdd3/przygoda/HADES/FAT_out/all_report/ecalacc/ecalacc_numbers.json` (md5 c39f6ffd...),
`ecalacc_fill.txt` (log), `tag_fits.txt`, `rows_c23.bin` + `rows_c23_layout.txt` (one record per simulated bump row),
`geo_numbers.json`, `emc_cells_sim_v2.txt` (the ECAL geometry and status of the simulation).

---------------------------------------------------------------------------------------------------
## 0. In plain words

**The question.** The bump is the decay eta -> pi+pi-pi0 with pi0 -> e+e-gamma (class C2) or pi0 -> gamma gamma with one photon
converted (C3), in which the remaining photon is NOT seen. If the ECAL sees that photon, the whole eta is reconstructed
("fully reconstructed eta", the ECAL tag of Part II, T = 3314 +- 105 +- 124 rows for "pre"). How much of the bump can this
full reconstruction see, given that the HADES ECAL covers only part of the polar angles, and what does one have to know
before it can be used as a measurement?

**What was done.**
1. **ECAL coverage from the geometry used in the simulation.** I read the ECAL geometry (EmcGeomPar) and the module status
   (EmcCellStatusPar) of the parameter file of the gen4 simulation production on vae25 with Hydra 2.7.6, and the same status
   in the data parameter file for eight runs of the four periods.
2. **Where does the lost photon go?** The simulation does not store the photon of a bump row. But it can be reconstructed
   from what is stored. The four charged particles are known (true momenta). The photon must give m_pi0 together with
   the e+e- pair and m_eta together with all four. These two conditions leave one free angle: the photon lies on a known
   circle. Its position on the circle is distributed according to the eta momentum spectrum of SMASH (read from the SMASH
   event files) and the decay matrix element (flat in HGeant). So every simulated bump row gets the full distribution of
   its photon direction and energy. I call this the "conditional completion".
3. **Check against the ECAL simulation.** For 598 rows the ECAL simulation did find the photon (cluster truth: a photon
   from the same pi0). These photons lie on the predicted circle (|Delta cos| < 1e-3 for 99.8 %). 99.5 % of them hit an
   active module of my geometry test. The predicted hit probability is calibrated: rows with a higher predicted
   probability have proportionally more found photons, with the same ratio 0.64-0.68 in every bin.
4. **Fractions.** For the accepted pi+pi-e+e- bump rows: photon in an ECAL module, photon energy above the 100 MeV cut,
   photon found by the ECAL simulation, row fully reconstructed. Each fraction is given against M(pi+pi-), against the
   photon energy and against the distance to the bump edge.
5. **Prediction for the data.** The simulated tag was scaled with the measured eta scale (f_eta(edge) = 0.392 as is,
   0.418 C2-corrected, 0.560 with conversions x 0.48) and compared with the data's T, also per M(pi+pi-) slice.
   - The ECAL energy is miscalibrated in a different way in data (+17 %) and in the simulation (-7 to -9 %). Therefore the
     tag was also recomputed with the simulated cluster energies x lambda (lambda = 1.10, 1.20, 1.285 = data-like).
   - Two effects that exist only in the data were added: the wider pi0 peak, and the START dependence of the photon
     beta cut.

**What came out.**
- **ECAL coverage.** 6 sectors x 163 lead-glass modules (92 x 92 mm, 420 mm deep, about 2.4 m from the target), polar angle
  **11.5-47.5 deg** (seen from z = 0; 11.2-46.6 deg from z = -60 mm), about 54 of every 60 deg in azimuth.
  - **Sector 0 (azimuth 63-117 deg) is switched off completely**, in the simulation and in the data parameters alike.
  - 8 single modules are dead; the status map is the same in the data (all four periods) and in the simulation.
  - Active modules: 807 of 978, 9.2 % of the full solid angle (11.2 % if all six sectors worked).
  - The data photons confirm this: their azimuth and polar-angle distributions agree bin by bin with the simulation, and
    sector 0 is empty in both (Fig. ecalacc_geometry).
- **Where the lost photon goes** (accepted bump rows C2+C3, "pre"):
  - 11 % fly below 11.5 deg (into the beam hole of the ECAL);
  - 29 % fly above 47.5 deg (where there is no ECAL; the tracking acceptance goes to 85 deg);
  - 60 % fly into the ECAL polar band, but sector 0, the sector gaps and the module edges cut this down;
  - **32.7 +- 0.4 % hit an active ECAL module** (40.3 % with six working sectors). This is almost the same for C2 and C3,
    the same in "pre" and S0, and flat in M(pi+pi-).
- **The photon is soft.** The mean true energy is 130 MeV (172 MeV inside the modules), and 51 % of the photons have less
  than 100 MeV (34 % inside the modules).
  - The analysis cut is E_cluster > 100 MeV. Only **21.6 % of the rows have their photon in a module and above 100 MeV
    (true energy)**.
  - Near the bump edge (M4 within -15...+5 MeV of the edge) this falls to 13 %. Far below the edge (more than 60 MeV
    below) it is 34 %.
  - This is kinematics: at the edge the photon takes the smallest possible energy.
- **What the ECAL simulation finds.**
  - Inside the active modules it finds a cluster of the true photon in 65 % of the cases, at all energies.
  - With the quality cuts (neutral candidate, 0.8 < beta < 1.2, E_cluster > 100 MeV) it keeps 38 %. The efficiency is 0
    below 80 MeV and about 60-75 % above 160 MeV.
- **Fully reconstructed** (tag with the true photon, m5c in 0.49-0.62):
  - **9.8 +- 0.6 % of the accepted bump rows** in the simulation as it is (14.1 % with a data-like energy scale);
  - in the data T / (visible bump yield) = 3314 / about 28 600 = **about 12 %**.
  - Per region:

    | region | as simulated | data-like energy scale |
    |---|---|---|
    | near the edge | 6.9 % | 11 % |
    | in the bump window B | 8.4 % | 13 % |
    | far below the edge | 16 % | 20 % |

  - 23 % of the simulated C2+C3 tags in the m5c window use a wrong photon. The frozen fit counts essentially only the
    true-photon tags: T_sim = 7295 against 7175 true-photon tags.
- **Data against the prediction** (pre, all M(pi+pi-); T_data = 3314 +- 105):

  | simulated photon-energy scale | expected T | data / expected |
  |---|---|---|
  | as is (f 0.392) | 2857 +- 153 (stat) +- 234 (syst) +- 400 (MC) | **1.16 +- 0.21** |
  | x 1.10 | | 0.94 |
  | x 1.20 | | 0.85 |
  | x 1.285 (data-like) | 4148 | 0.80 |
  | x 1.285 + data resolution (x 0.964) + START (x 0.978) | **3910 +- 210 +- 325 +- 560** | **0.85 +- 0.15** |

  - C2-corrected (0.418 x composition 0.928) and kappa = 0.48 (0.560 x 0.711) give the same as "as is" within 2 %.
  - The ratio data/simulation per M(pi+pi-) slice, f_tag = 0.392 / 0.505 / 0.452 (< 0.34 / 0.34-0.37 / 0.37-0.40), is 1.10-1.19
    times the edge's h_data/h_sim of the same slices. The tag/edge difference is a global factor, not a shape.
- **START.** The neutral-candidate beta is measured against the START time in the data. The 6.4 % of the e+e- rows without
  START have a tag rate of only 0.80 +- 0.12 (all M4) or 0.69 +- 0.11 (window B) of the others. The data lose about 2.3 %
  of T this way; the simulation does not have this loss.
- **By-product: the C2 gamma\* angle factor.** The circle average of the QED weight 1 + cos^2 + (4m^2/q^2) sin^2 gives the
  C2 factor of the whole simulated sample directly:
  - 0.863 +- 0.003 (all M4), 0.861 in B, **0.855 +- 0.005 at the edge**, 0.876 +- 0.009 for the tagged rows;
  - bedge uses 0.8855 (WP8 toy 0.880, closure 0.891). The new value is 2.5-3.5 % lower; at the edge it would raise
    f(C2 corr.) by +1.8 % (0.418 -> 0.425). This is within the E8 uncertainty (+-0.024 on the factor), but on its edge.

**What it means.**
- **Fraction of the bump the full reconstruction can see.**
  - Geometry: one third of the accepted bump rows has its photon in the ECAL. The limit is the polar-angle coverage
    (11.5-47.5 deg) and the dead sector 0.
  - After the 100 MeV threshold and the cluster efficiency, about 10-14 % are fully reconstructed (data: about 12 %).
  - The full reconstruction sees mainly the part of the bump far below the edge (hard photons). At the edge itself only
    7-11 % are seen.
  - So the tag and the edge measure different parts of the bump, which makes the tag a genuinely independent check.
- **What it costs in statistics: nothing, rather the opposite.**
  - The tagged sample keeps only about 12 % of the bump, but it is a clean peak (S/B = 0.73). Its statistical error
    (3.2 %) is smaller than that of the edge (5.4 %).
  - The simulation statistics limit both equally: MC error of T_sim 8.9 %, of the edge 11 %. The standard set has only
    2836 C2+C3 rows in "pre", of which about 390 are tagged.
- **The real cost is the systematics of the photon.** The data/expected ratio moves from 1.16 to 0.80 when the threshold
  is applied at the simulated or at a data-like energy scale (factor 1.46). Today the tag therefore cannot test the eta
  scale better than about +-20 %. It agrees with the edge within that.
- **What must be known from the simulation** (input for the work after the ECAL recalibration; Sect. 6):
  1. the photon energy scale in data and in simulation, because the 100 MeV threshold and the pi0 window act on it;
  2. the efficiency inside the modules versus the true energy (65 % any cluster, why not more?);
  3. the energy resolution (data pi0 13.7-14.3 MeV wide, simulation 9.6 MeV);
  4. the START dependence of the beta cut;
  5. a much larger simulation of eta -> pi+pi-pi0.
  The acceptance itself (geometry and dead modules) is known and agrees between data and simulation.

---------------------------------------------------------------------------------------------------
## 1. What was done, step by step

1. **Geometry** (`ecalacc_emcgeo_vae25.C`, `ecalacc_emcstat_vae25.C`, run on vae25; `ecalacc_geo.py` locally).
   - Source: `/lustre/hades/user/hadesdst/sub/feb22/sim/p4500p/gen4/smash/smash_lepton/dst/allParam_feb22_sim_run_17000_gen4_10032026.root`
     (the parameter file of the gen4 DST production of the SMASH simulation; analysisDST.cc there: paramSource root,
     release FEB22SIM_v4a, run 17000).
   - All 255 module slots of each sector were read. The lab position of every module was computed as in Hydra
     `emc/hemcclusterf.cc` (module transform, then the sector lab transform).
   - The status is indexed by the tile number (emc_cell2tile_map), as in `hemcdigitizer.cc`. The digitizer skips status-0
     modules, so they are dead in the simulation.
   - Data: `/cvmfs/hadessoft.gsi.de/param/real/feb22/allParam_Feb22_gen4_07032026.root` at eight runs (first runs of days
     034, 039, 045, 047, 049, 053, 059, 067).
   - Per sector the module fronts lie in one plane (rms 0.003 mm). Solid angle: 4e6 isotropic rays.
2. **SMASH eta spectrum.**
   - Read: the eta lines (GEANT id 17) of the converted SMASH files of the standard set, which are the HGeant input of
     the simulation (`/lustre/hades/dstsim/feb22/p4500p/smash/convert/smash_standard/1/p_p_4500MeV_100000evts_<i>.evt`).
   - 140 572 etas. The transfer stopped at my 15-min ssh timeout after about 60 of the 120 requested files; the partial
     last line was dropped.
   - Density: E d3N/dp3 ~ dN/(dy dmT)/mT in 70 x 45 bins (y_lab -0.5...3.0, mT - m 0-1.8 GeV), bilinear interpolation.
3. **Row pass** (`ecalacc_fill.C`, 31 s, about 10 GB).
   - **Data, frozen tag** (good photon, 0.119 <= M(e+e-gamma) < 0.173, M(ee) < 0.10; m5c; frozen fit 0.5463 / 6.2 MeV,
     pol2, 0.49-0.62) for pre and S0, all and per M(pi+pi-) slice, and split by evtStart.
     - Also: photon maps of the tagged rows and M(e+e-gamma).
     - Reproduces T(pre) = 3313.7 +- 104.9 and T(S0) = 3090.3 +- 94.3 exactly.
   - **Simulation, emulated tag** (SimEngine copy, mode A, w x V3 = code 108, owned rows, std + lep + eta, expected raw rows
     = sum w x 27.0718).
     - The same frozen fit; reproduces WP6's T_sim(S0) = 6803.9 exactly.
     - pre: T_sim = 7295.0.
   - **C2/C3 rows** (2836 in pre; all from the standard set, which owns them): the PAT entry was read from
     `/hdd3/przygoda/HADES/files/sim/ppee_std_all.root`.
     - Identity check against the compact: 0 failures.
     - Read: true momenta of e+, e-, pi+, pi-, the vertex, the generator codes, and up to five neutral candidates with
       their cluster and GEANT truth.
   - **Generator masses** from the found photons: m_eta = 547.451, m_pi0 = 134.976 MeV/c2 (the HGeant values). These are
     used in the constraints.
   - **Conditional completion.**
     - Work in the rest frame of Q = e+e- (C3: the converted photon's pair). Then |k| = (m_pi0^2 - M_ee^2)/(2 M_ee), and the
       angle between k and P = pi+pi-e+e- is fixed by (P + k)^2 = m_eta^2.
     - The remaining azimuth phi has a flat measure (d^3k/k with the two delta functions is const x dphi).
     - Sample 720 points in phi. Weight each point with f_SMASH(P + k) x |M|^2, where |M|^2 = 1 for the simulation (HGeant
       phase space, isotropic gamma*).
     - Variants:
       - KLOE-2 Dalitz plot f(X, Y) along the circle;
       - QED gamma* angle (normalised to mean 1);
       - flat weight (no SMASH density).
     - 1 row of 2836 has no solution and none has zero density.
   - **Geometric test.** A ray from the true vertex, intersected with the sector plane, hits a module if it crosses its
     front face enlarged to the 94 mm pitch (the strict 92 mm faces are a variant). Status-0 modules do not count.
   - **Truth photon of the ECAL simulation.**
     - Match: a cluster with GEANT id 1, parent pi0, the generator codes of the pions, and the vertex of the pions within
       0.01 mm.
     - Flags: quality cuts passed; whether the compact's tagged photon is this one.
     - Also recorded: clusters made by an e+- of the converted photon (3.8 % of the rows).
   - **Recomputed tag.**
     - From the PAT clusters (FAT logic: quality cuts, the three most energetic, closest to 0.140) with cluster energies
       x lambda, lambda = 1, 1.10, 1.20, 1.285.
     - For lambda = 1 it agrees with the compact for all 2836 rows.
4. **Numbers** (`ecalacc_ana.py`).
   - Fractions of the weight with independent-row MC errors.
   - The expected tag in the data: T_exp = f_eta(edge) x T_sim(frozen fit) x composition x calibration.
     - Composition: (c_tag C2 + C3)/(C2 + C3), with c_tag = 0.876 for the C2-corrected reading; (c_tag C2 + 0.48 C3)/(C2 + C3)
       for kappa = 0.48.
     - Calibration: the frozen fit (copied to Python; it reproduces 7295.0 / 6803.9) of the T1+T2 spectrum in which the
       C2+C3 rows are re-tagged at lambda.
     - Errors: stat, syst, MC of f from bedge; the MC of T_sim (8.9 %, WP6 replicas) is added in quadrature to the MC of f.
       This is conservative, because both come from the same standard-set rows.
   - Data-only factors:
     - resolution: the pi0 window acceptance of a Gaussian, data / simulation x 1.285 = 0.931/0.966 = 0.964;
     - START: 1 - 75/3389 = 0.978.
5. **Figures** (`ecalacc_figs.py`).

---------------------------------------------------------------------------------------------------
## 2. Numbers

### 2.1 ECAL coverage (`tables_all/ecalacc_geometry.{txt,tex}`)

| sector | modules | active | theta from z = 0 [deg] | theta from z = -60 mm | phi from z = 0 [deg] |
|---|---|---|---|---|---|
| 0 | 163 | **0** | 11.6-47.6 | 11.3-46.6 | 62.8-117.1 |
| 1 | 163 | 161 | 12.0-47.3 | 11.7-46.3 | 123.3-177.0 |
| 2 | 163 | 163 | 11.7-47.1 | 11.4-46.2 | 183.5-237.0 |
| 3 | 163 | 161 | 11.7-47.1 | 11.4-46.1 | 242.9-296.4 |
| 4 | 163 | 161 | 11.6-47.0 | 11.3-46.1 | 302.6-356.4 |
| 5 | 163 | 161 | 11.4-47.0 | 11.2-46.0 | 2.6-56.7 |

- Modules: 92.0 x 92.0 mm front, pitch 94.0 mm, depth 420 mm; front planes at 2357-2411 mm from the target centre.
- Dead single modules (sector/module index): 1/24, 1/195, 3/160, 3/237, 4/129, 4/157, 5/240, 5/244.
- The data status is identical at all eight runs.
- Solid angle of the active modules: 9.20 % of 4 pi from z = 0 (8.84 % from z = -60 mm); six full sectors 11.2 %; 60 % of the
  11.5-47.6 deg band.
- Data check (Fig. ecalacc_geometry b, c): the tagged-photon azimuth and polar angle of data and simulation agree bin by bin
  (both empty at 63-117 deg and above 47.5 deg).

### 2.2 The photon of the accepted bump rows (`tables_all/ecalacc_results.{txt,tex}`; pre unless stated; per cent; +- MC)

| quantity | C2 | C3 | C2+C3 | C2+C3, window B | C2+C3, edge (-15...+5 MeV) | C2+C3, S0 |
|---|---|---|---|---|---|---|
| rows (expected raw rows, SMASH scale) | 1649 (42 162) | 1187 (31 057) | 2836 (73 219) | 2152 (58 368) | 678 (19 633) | 2622 (67 764) |
| photon hits an active module | 32.9 +- 0.5 | 32.3 +- 0.5 | **32.7 +- 0.4** | 31.9 | 31.5 | 32.9 |
| the same, six full sectors | 40.5 | 40.1 | 40.3 | 39.5 | 38.5 | 40.6 |
| strict 92 mm faces | 31.5 | 31.0 | 31.3 | 30.5 | 30.2 | 31.5 |
| flat weight on the circle / KLOE / QED | 33.6 / 32.9 / 33.2 | 33.3 / 32.3 / - | 33.5 / 32.7 / - | | | |
| in a module and E_true > 80 / 100 / 130 MeV | 24.3 / 21.4 / 17.6 | 24.8 / 21.9 / 17.9 | 24.5 / **21.6** / 17.7 | 22.5 / 19.1 / 14.8 | 17.0 / **13.0** / 8.2 | - / 21.8 / - |
| found by the ECAL simulation (any cluster) | 21.5 +- 1.1 | 21.1 +- 1.3 | 21.3 +- 0.9 | 20.8 | 21.9 | 21.6 |
| found and quality cuts | 12.3 +- 0.9 | 12.6 +- 1.1 | 12.4 +- 0.7 | 10.5 | 7.6 | 12.7 |
| tagged, any photon, m5c 0.49-0.62 | 12.9 | 12.5 | 12.7 +- 0.7 | 11.4 | 10.2 | 13.0 |
| **fully reconstructed** (true photon, m5c window) | 9.9 +- 0.8 | 9.6 +- 0.9 | **9.8 +- 0.6** | 8.4 | **6.9** | 10.0 |
| the same, cluster energies x 1.285 | 14.0 | 14.3 | **14.1 +- 0.7** | 13.0 | 11.1 | 14.4 |

- Polar angle of the true photon (C2+C3): < 11.5 deg 10.8 %; 11.5-47.5 deg 60.4 %; > 47.5 deg 28.8 %. In 18-85 deg (the
  tracking acceptance): 76 %.
- True photon energy: mean 130 MeV (all), 172 MeV (in the modules). Below 100 MeV: 51 % (all), 34 % (in the modules).
- Efficiency inside the active modules (C2+C3): found 0.653, found + quality 0.380, fully reconstructed 0.299 (S0: 0.657,
  0.386, 0.304).
- Versus the true energy (Fig. ecalacc_energy b):
  - any cluster: about 0.55-0.75, flat;
  - with the quality cuts: 0 below 60 MeV, 0.02 at 60-80, 0.05 at 80-100, 0.16 at 100-120, 0.50 at 120-160, 0.64-0.75 above
    160 MeV.

**Versus M(pi+pi-)** (C2+C3, pre):

| M(pi+pi-) | weight | in modules | and E > 100 MeV | found + quality | fully reconstructed | the same x 1.285 | <E_gamma> all / in modules [MeV] |
|---|---|---|---|---|---|---|---|
| < 0.34 | 24 239 | 33.2 +- 0.5 | 23.3 | 13.0 +- 0.9 | 9.8 +- 0.8 | 14.2 | 152 / 194 |
| 0.34-0.37 | 25 419 | 33.3 +- 0.6 | 22.6 | 12.8 +- 1.2 | 9.4 +- 1.0 | 14.5 | 134 / 177 |
| 0.37-0.40 | 20 837 | 31.3 +- 0.8 | 19.3 | 11.9 +- 1.5 | 10.6 +- 1.4 | 14.5 | 115 / 147 |
| 0.40-0.43 | 2 683 | 32.5 +- 2.9 | 16.2 | 8.1 +- 3.9 | 6.1 +- 3.4 | 8.1 | 108 / 117 |

**Versus the distance to the edge** d = M4' - E(M(pi+pi-), M(ee)) (C2+C3, pre):

| d [MeV/c2] | weight | in modules | and E > 100 MeV | found + quality | fully reconstructed | the same x 1.285 | <E_gamma> all / in modules |
|---|---|---|---|---|---|---|---|
| -200...-100 | 1 913 | 37.5 | 37.1 | 23.5 | 18.6 | 19.8 | 319 / 404 |
| -100...-60 | 10 156 | 36.0 | 33.5 | 21.3 | 15.9 | 20.1 | 231 / 293 |
| -60...-40 | 12 943 | 32.7 | 26.8 | 14.6 | 11.1 | 16.1 | 165 / 216 |
| -40...-30 | 8 571 | 32.8 | 24.1 | 11.1 | 7.7 | 13.9 | 135 / 172 |
| -30...-20 | 11 366 | 31.3 | 19.5 | 12.6 | 10.9 | 14.3 | 108 / 136 |
| -20...-10 | 12 467 | 32.4 | 17.3 | 9.0 | 6.7 | 12.3 | 91 / 117 |
| -10...0 | 10 797 | 32.2 | 12.2 | 7.7 | 7.2 | 10.7 | 74 / 89 |
| 0...+20 | 4 508 | 28.2 | 9.4 | 5.4 | 5.1 | 7.7 | 71 / 85 |

### 2.3 Closure of the conditional completion

- 598 photons found by the ECAL simulation (cluster truth).
- **On the predicted circle:**
  - residuals (5/50/95 %): M(e+e-gamma) - m_pi0 = 0.0000 / 0.0000 / 0.15 MeV; M(5) - m_eta = -0.004 / 0.000 / 0.10 MeV;
  - |Delta cos alpha| < 1e-3 for 99.8 %.
- **In the geometry:** 99.50 % inside the predicted active modules, 0.33 % in a dead module (the cluster lands in a
  neighbour), 0.17 % (1 photon, at 36.4 deg) outside.
- **Calibration of the predicted probability** (rows binned in it; found share / predicted probability; this ratio is the
  efficiency inside the modules):

  | predicted probability | rows | found / predicted |
  |---|---|---|
  | 0-0.05 | 158 | 0 photons found (0.6 % predicted) |
  | 0.05-0.2 | 446 | 0.56 +- 0.10 |
  | 0.2-0.4 | 1338 | 0.68 +- 0.04 |
  | 0.4-0.6 | 726 | 0.64 +- 0.04 |
  | 0.6-0.8 | 132 | 0.68 +- 0.07 |
  | 0.8-1 | 35 | 0.63 +- 0.10 |

  Flat, so the photon-direction distribution along the circle is right.
- **Simulated energy scale:** E_cluster / E_true (16/50/84 %) = 0.752 / 0.914 / 1.043 for all found photons, and 0.808 /
  0.931 / 1.048 for those passing the quality cuts.

### 2.4 The tag: data against the prediction (`tables_all/ecalacc_tag.{txt,tex}`; Fig. ecalacc_tag)

| sel, M(pi+pi-) | data T | sim T (frozen) | f_tag | expected, as is (0.392) | expected, C2 corr. x 1.285 | data / exp. as is | data / exp. x 1.285 |
|---|---|---|---|---|---|---|---|
| pre, all | 3314 +- 105 | 7295 | 0.454 +- 0.014 | 2857 +- 153 +- 234 +- 400 | 4148 +- 223 +- 345 +- 594 | 1.16 +- 0.21 | 0.80 +- 0.14 |
| pre, < 0.34 | 1014 +- 71 | 2588 | 0.392 +- 0.027 | 1013 | 1472 | 1.00 | 0.69 |
| pre, 0.34-0.37 | 1204 +- 57 | 2385 | 0.505 +- 0.024 | 934 | 1410 | 1.29 | 0.85 |
| pre, 0.37-0.40 | 1116 +- 44 | 2468 | 0.452 +- 0.018 | 966 | 1358 | 1.15 | 0.82 |
| pre, 0.40-0.43 | 38 +- 15 | 123 | 0.31 +- 0.12 | 48 | 58 | 0.80 | 0.66 |
| S0, all | 3090 +- 94 | 6804 | 0.454 +- 0.014 | 2664 +- 143 +- 218 +- 373 | 3910 +- 210 +- 325 +- 560 | 1.16 +- 0.21 | 0.79 +- 0.14 |

- **Readings (pre, all):**

  | reading | expected T | data / expected |
  |---|---|---|
  | C2 corrected | 2832 | 1.17 |
  | kappa 0.48 | 2905 | 1.14 |
  | C2 corrected, x 1.10 | 3534 | 0.94 |
  | C2 corrected, x 1.20 | 3897 | 0.85 |
  | C2 corrected, x 1.285 | 4148 | 0.80 |
  | **best** (x 1.285, resolution 0.964, START 0.978) | **3910 +- 210 +- 325 +- 560** | **0.85 +- 0.15** |
  | S0 best | 3685 | 0.84 +- 0.15 |

- **Calibration factor of the simulated tag** (frozen fit; truth ratio of the tags with any photon in brackets):
  x 1.245 (1.192) at lambda 1.10, x 1.374 (1.279) at 1.20, x 1.464 (1.326) at 1.285. Per slice (true-photon tags): 1.45, 1.54,
  1.37, 1.33.
- **pi0 in M(e+e-gamma)** (Gaussian + linear fit, 0.095-0.20):

  | sample | position / width [MeV/c2] | window 0.119-0.173 keeps |
  |---|---|---|
  | data, all e+e- rows with a good photon | 141.9 / 14.3 | 0.931 |
  | data, tagged eta, from report_exp | 146.0 / 13.7 | 0.951 |
  | simulation C2+C3 as is | 129.4 / 9.6 | 0.860 |
  | simulation C2+C3, x 1.285 | 144.9 / 12.7 | 0.966 |

- **f_tag against the edge per slice** (bedge, uncorrected h_data/h_sim):
  - edge: 0.342 (0.30-0.34), 0.459 (0.34-0.37), 0.376 / 0.389 (0.37-0.39 / 0.39-0.40);
  - tag / edge: 1.15, 1.10, about 1.19. A global factor.
- **START** (data pre):
  - T with START = 3145 +- 102; T without START = 170 +- 26;
  - e+e- rows without START: 6.4 % (all) and 7.2 % (window B);
  - relative tag rate without START: 0.80 +- 0.12 (all M4), 0.69 +- 0.11 (window B);
  - if those events tagged like the others: +75 rows (+2.3 %).

### 2.5 C2 gamma\* angle factor from the circle average (by-product)

| sample | factor |
|---|---|
| all C2 rows (pre) | 0.8633 +- 0.0026 (E1 0.8647, E2 0.8624) |
| window B | 0.8608 +- 0.0029 |
| edge (-15...+5 MeV) | **0.8554 +- 0.0050** |
| far below (< -60 MeV) | 0.8724 +- 0.0057 |
| tagged C2 rows | 0.876 +- 0.009 (exact per row where the photon is known: 0.883 +- 0.011) |

- bedge: 0.8855 recommended (toy 0.880 / 0.883 at the edge, closure 0.891).
- Effect on f(C2 corrected) if the edge value 0.855 were used: 0.418 -> 0.425 (+1.8 %).
- The QED weight hardly changes the photon acceptance (32.9 -> 33.2 % for C2).

---------------------------------------------------------------------------------------------------
## 3. Figures and tables (draft captions)

**Fig. `figs_all/ecalacc_geometry.pdf`.**
- **What it shows: the ECAL covers 11.5-47.5 deg in polar angle in five of the six sectors, the data photons fill exactly
  this coverage, and two thirds of the lost bump photons fly outside it.**
- (a) The 978 ECAL modules seen from the target centre (theta, phi), from the EmcGeomPar and EmcCellStatusPar of the gen4
  simulation parameter file:
  - light blue: active (807);
  - grey: sector 0, switched off (163);
  - red: dead single modules (8).
  - The data parameters of all four periods have the same status. Dotted: 18 deg, where the tracking acceptance starts.
- (b) Azimuth of the ECAL photon of the tagged rows (good photon, 0.119 <= M(e+e-gamma) < 0.173, M(ee) < 0.10, "pre"): data
  (black) and simulation T1 + T2 (orange), each normalised to 1.
- (c) Polar angle:
  - magenta: the TRUE photon of all accepted bump rows C2+C3 (conditional completion, SMASH x V3);
  - black / orange dashed: the cluster polar angle of the tagged rows in data and simulation, scaled to the magenta share
    in 11-48 deg;
  - blue band: the ECAL coverage.

**Fig. `figs_all/ecalacc_fractions.pdf`.**
- **What it shows: one third of the bump photons hit the ECAL at every M(pi+pi-), but the photons near the bump edge are too
  soft for the 100 MeV threshold, so the full reconstruction sees mainly the part of the bump far below the edge.**
- Simulation, accepted pi+pi-e+e- rows of the bump classes C2+C3, "pre", SMASH x V3, MC errors.
- Curves:
  - photon in an active module (blue), the same with six full sectors (light blue);
  - photon in a module and true energy above 100 MeV (magenta);
  - found by the ECAL simulation with the quality cuts (dark grey);
  - fully reconstructed, i.e. tagged with the true photon and m5c in 0.49-0.62 (orange); the same with the simulated
    cluster energies x 1.285, the data-like energy scale (orange open, dashed).
- (a) Versus M(pi+pi-) (data scale; dotted: kinematic limit 0.4116).
- (b) Versus the distance to the bump edge, M4' - E(M(pi+pi-), M(ee)).

**Fig. `figs_all/ecalacc_energy.pdf`.**
- **What it shows: half of the lost photons have less than 100 MeV; inside the modules the simulated ECAL finds a cluster
  for about 65 % of the photons at every energy, but the quality cuts remove almost everything below 120 MeV.**
- (a) True photon energy of the accepted bump rows C2+C3 ("pre"):
  - all (magenta) and in an active module (blue), from the conditional completion;
  - photons found by the ECAL simulation (any cluster: grey; passing the quality cuts PID 1, 0.8 < beta < 1.2,
    E_cluster > 100 MeV: orange), from the cluster truth.
- (b) Efficiency inside the active modules = found / in a module, versus the true energy.
  - Text: the simulated cluster energy is 7-9 % low (median E_cluster/E_true = 0.931 for clusters passing the cuts); the
    data photon energy is about 17 % high (Part II).

**Fig. `figs_all/ecalacc_tag.pdf`.**
- **What it shows: the number of ECAL-tagged eta in the data agrees with the eta scale of the bump edge within the errors,
  but the prediction moves by almost 50 % with the assumed photon-energy scale of the 100 MeV threshold.**
- (a) Data, "pre", tagged rows: N - CB in m5c (blue; CB red) with the frozen fit (green; T = 3314).
  - Orange band: the simulated C2+C3 tagged rows x f_eta(edge) = 0.392.
  - Orange dashed: the same with the data-like energy scale, data resolution and START factors.
- (b) T per M(pi+pi-) slice:
  - data (black);
  - expected f_eta x T_sim (open squares);
  - best expectation (C2-corrected f = 0.418, cluster energies x 1.285, resolution 0.964, START 0.978; diamonds).
  - Errors of the expectation: stat (+) syst (+) MC of f and the MC of T_sim.
- (c) Data T / expected T (pre, all M(pi+pi-)) for the readings of the eta scale and the assumed photon-energy scales.

**Fig. `figs_all/ecalacc_closure.pdf`.**
- **What it shows: the photon predicted for each bump row is the photon that the ECAL simulation finds, and the ECAL energy
  scale differs between data and simulation.**
- (a) For the 598 photons found by the ECAL simulation: cos(alpha) of the true photon minus the predicted cos(alpha) of the
  circle (alpha = angle to P = pi+pi-e+e- in the e+e- rest frame).
- (b) Simulated E_cluster/E_true for all found photons (orange) and those passing the quality cuts (grey). Dashed: the data
  scale (+17 %).
- (c) pi0 in M(e+e-gamma):
  - data, all e+e- rows with a good photon, N - CB (blue);
  - simulation C2+C3 as simulated (orange) and with the cluster energies x 1.285 (orange dashed);
  - grey: the tag window 0.119-0.173.

**Tables** (`docs/latex/tables_all/`; the .tex versions compile with booktabs, tested in the scratchpad):
- `ecalacc_geometry.{txt,tex}`;
- `ecalacc_results.{txt,tex}` (fractions, including the versus-M(pi+pi-) and versus-distance lists in the .txt);
- `ecalacc_tag.{txt,tex}` (data vs expectation; the .txt also has all readings, calibration factors, the pi0 fits and START).

---------------------------------------------------------------------------------------------------
## 4. Caveats

1. **The expected T is MC-limited and correlated with f.** T_sim comes from the same 2836 standard-set rows as E_sim of the
   edge. Adding their MC errors in quadrature (as done) overestimates the error of the ratio.
2. **lambda is not a calibration.** It is one global factor on the cluster energy. The real miscalibration may depend on
   energy and module, and the simulated resolution was not adjusted (only the pi0-window factor 0.964).
   - The data-like value is 1.20-1.29, depending on which pi0 position is used: all e+e- rows with a photon give 141.9 MeV,
     the tagged eta's pi0 gives 146.0 MeV (report_exp), the simulation 129.4 MeV.
   - Over that range the prediction moves by 7 %.
3. **The 35 % inefficiency inside the modules is not explained.**
   - Known parts: about 4 % of the rows have a cluster made by an e+- of the converted photon (not counted as found).
   - Unknown parts: cluster-track matching (a neutral candidate needs no matched track), cluster splitting, module edges.
   - Only the cluster truth of the found photons is stored, so the lost ones cannot be followed.
4. **The SMASH eta density** comes from about 60 of the 100 000 converted files of the standard set (140 572 etas, about 6e6 SMASH events). The flat-weight
   variant changes the in-module fraction by only +0.8 points, so the result does not depend on the details.
5. **The geometric test uses the nominal module fronts.** Showers near a sector edge or a dead module can still give a
   cluster in a neighbour (0.3 % of the found photons). Front faces enlarged to the pitch are nominal; the strict 92 mm faces
   give -1.4 points.
6. **"Fully reconstructed" here = tagged with the true photon and m5c in 0.49-0.62.** The data T is the Gaussian content of
   the frozen fit. In the simulation the two agree (7295 vs 7175), but this was not checked per slice.
7. **The kinematic distance to the edge of the tagged eta could not be checked in the data.** The m5c side bands have a
   different M4 distribution than the peak window (the subtraction goes negative above the edge). This check needs a 2-D
   fit, not side bands.
8. **The C2 factor from the circle average** assumes that the SMASH eta spectrum and the HGeant decay model describe the
   simulated sample (they do by construction). It is not a measurement of the physics. It disagrees with 0.8855 by 2.5-3.5 %;
   reported for information, nothing was changed in bedge.
9. **Own mistakes, repaired** (the superseded files are kept, not deleted):
   - the first geometry dump read only module slots 0-162 (of 255) and indexed the status by module instead of tile. Files
     `emc_cells_sim.txt`, `emcgeo_vae25_stdout.txt`, `emcstat_real_vae25_stdout.txt` are superseded by the `_v2` files. It was
     caught by the closure (31 % of the found photons outside the coverage before, 0.2 % after);
   - test run outputs `ecalacc_fill_test*.txt` (100 000-eta density) are superseded.

---------------------------------------------------------------------------------------------------
## 5. Exact commands (from /home/przygoda/HADES/FAT/FAT)

    O=/hdd3/przygoda/HADES/FAT_out/all_report/ecalacc
    # 1. geometry and status (vae25, read only; output to stdout)
    (cat plots/all_report/ecalacc_emcgeo_vae25.C; echo ".q") | ssh -o RemoteCommand=none -o RequestTTY=no vae25 \
       'source /lustre/hades/user/przygoda/pat3/profile.sh >/dev/null 2>&1; cd /tmp; root -l -b 2>&1' > $O/emcgeo_vae25_stdout_v2.txt
    grep "^CELL" $O/emcgeo_vae25_stdout_v2.txt > $O/emc_cells_sim_v2.txt
    (cat plots/all_report/ecalacc_emcstat_vae25.C; echo ".q") | ssh -o RemoteCommand=none -o RequestTTY=no vae25 \
       'source /lustre/hades/user/przygoda/pat3/profile.sh >/dev/null 2>&1; cd /tmp; root -l -b 2>&1' > $O/emcstat_real_vae25_stdout_v2.txt
    # 2. SMASH eta momenta (converted standard-set files; stopped by a 900 s timeout after about 60 files)
    ssh -o RemoteCommand=none -o RequestTTY=no vae25 'd=/lustre/hades/dstsim/feb22/p4500p/smash/convert/smash_standard/1; \
       for i in $(seq 1 120); do awk "NF==9 && \$5==17 {print \$1, \$2, \$3, \$4}" $d/p_p_4500MeV_100000evts_$i.evt; done' > $O/smash_eta_lab_120files.txt
    awk 'NF==4' $O/smash_eta_lab_120files.txt > $O/smash_eta_lab.txt          # 140 572 etas, md5 028c03e8...
    # 3. rows, fits, truth (31 s, about 10 GB)
    root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_ecalacc", kTRUE)' \
       'plots/all_report/ecalacc_fill.C+O' > $O/ecalacc_fill_stdout.txt 2>&1
    # 4. numbers, tables, figures
    python3 -I plots/all_report/ecalacc_geo.py
    python3 -I plots/all_report/ecalacc_ana.py
    python3 -I plots/all_report/ecalacc_figs.py

**Code** (`plots/all_report/`): `ecalacc_emcgeo_vae25.C`, `ecalacc_emcstat_vae25.C`, `ecalacc_fill.C`, `ecalacc_geo.py`,
`ecalacc_ana.py`, `ecalacc_figs.py`.
- Typed inputs, with their sources at the top of `ecalacc_ana.py`: f_eta(edge) readings and their errors (bedge), C2 factor
  0.8855 and kappa 0.48 (bedge), T(pre) (report_exp), the WP6 MC error of T_sim, the bedge yields.

---------------------------------------------------------------------------------------------------
## 6. Input for the work after the ECAL recalibration (what the simulation and the data must provide)

1. **Energy scale, data and simulation.**
   - The cut E_cluster > 100 MeV and the M(e+e-gamma) window act on the cluster energy. Today the data are about +17 % and
     the simulation -7 to -9 % off the true energy, so the same cut is about 85 MeV true in the data and about 108 MeV in the
     simulation.
   - This alone changes the expected tag by up to x 1.46.
   - After the recalibration, the simulation must be produced (or reweighted) with the same energy scale and resolution as
     the data. Alternatively, define the threshold in TRUE energy units in both.
   - Check: the pi0 position and width in M(e+e-gamma) of the tagged eta must agree (today 146.0 / 13.7 MeV in the data
     against 129.4 / 9.6 MeV in the simulation).
2. **Efficiency inside the modules versus energy** (simulation 65 % for any cluster, 38 % with the cuts).
   - The losses should be understood: conversions in front of the ECAL, neutral-candidate definition (track matching),
     cluster splitting.
   - Validate in data, e.g. with the pi0 -> e+e-gamma photon of the inclusive e+e- sample.
3. **The beta cut 0.8-1.2 uses the START time in the data.** Events without START lose 20-30 % of their tags. Either make
   the cut START-free or simulate the START efficiency.
4. **The acceptance is known.** Same module status in data and simulation; sector 0 off; the data photon maps agree with the
   simulation. Only the vertex distribution (data tails to -260 mm) moves theta by about 1 deg per 6 cm.
5. **Simulation statistics.** The standard set provides 2836 C2+C3 rows in "pre" (about 390 tagged): 9 % MC on T_sim, the
   same order as for the edge (11 %). A forced eta -> pi+pi-pi0 (pi0 -> e+e-gamma and gamma gamma) production with the full
   ECAL simulation, at least 10 times the current statistics, would make both the tag and the edge data-limited.
6. **Physics weights of the tagged rows.**
   - C2 gamma\* factor of the tagged rows: 0.876 +- 0.009, not 0.8855.
   - KLOE-2 Dalitz-plot weights do not change the photon acceptance (32.9 % with and without).
7. **What the tag can give after that.** The tag sees mostly the low-M4 part of the bump (photon energy far above
   threshold), the edge sees the top. Together they test the bump shape along the "distance to the edge" coordinate:
   fully-reconstructed share 16-20 % far below the edge, 7-11 % at the edge.

---------------------------------------------------------------------------------------------------
## 7. Verification (independent verifier, prefix vecalacc, 2026-11-15)

### 7.0 In plain words

**What I did.** I rebuilt the decisive numbers with my own code. I did not use the package's macros, its geometry dump, its
row table or its SMASH file.
1. **Geometry.** I wrote my own dump of the ECAL geometry and status on vae25, with a different status lookup (the
   digitizer's `HEmcDetector::getPositionFromCell`). I also read the raw status array of the data at four runs. My own
   hit test uses a module grid per sector plane, not polygons.
2. **Empirical check.** As a check that does not use the parameter files, I took every neutral candidate of the
   standard-set PAT file and looked at which modules carry clusters and where the true photons that made them point.
3. **Rows.** I made my own selection of the C2, C3, C4 and C5 rows with the frozen definitions, and read their true
   momenta from the PAT files myself.
4. **Completion.** I solved the missing photon in the LAB frame instead of the e+e- rest frame. There the photon lies on a
   cone around b p_Q - a p_P, its energy follows from the cone position, and the measure along the cone is E_gamma d(psi).
   I checked this numerically against the package's flat azimuth in the e+e- rest frame (chi2/ndf 1.01).
5. **SMASH spectrum.** I extracted the SMASH eta momenta from 40 OTHER converted files (1001-1040; 94 314 etas) and
   estimated the density with a different estimator.
6. **Tag.** I recomputed the ECAL tag at lambda = 1 ... 1.285 from the PAT clusters.

**What came out.**
- **Every central number of the package is reproduced to within 0.1-0.2 points, and the result does not depend on the
  SMASH sample.** One third of the lost photons hit an active module: 32.6 % (package 32.7 %). After the threshold and the
  cluster efficiency, about 10-14 % of the rows are fully reconstructed. The prediction moves by about x 1.45 with the
  energy scale. The tag is specific to the bump: the peak classes are tagged in only 0.8-1.8 % of the rows.
- **Geometry and status are confirmed by the clusters themselves.** The simulated clusters sit in exactly the 807 active
  modules: none in sector 0 and none in the 8 dead modules. The status is identical in data and simulation.

**Two descriptive numbers are wrong. Neither changes any fraction, tag or ratio.**
1. **The polar-angle split of the lost photon.** The note says 11 % below 11.5 deg, 60 % in the band and 29 % above
   47.5 deg. It should say **9.0 %, 54.5 % and 36.5 %**.
   - The package divided by the photons with theta < 90 deg only. 9.9 % of the photons go backward and were left out.
   - The share in the tracking range 18-85 deg is **68.9 %**, not 76 %.
2. **"Seen from the target centre" actually means seen from z = 0.** The target centre is at z = -124 mm (simulated
   vertex uniform in -148...-99 mm, data median -125 mm).
   - From the target centre the module corners span **10.9-45.6 deg** (more than half of the azimuth is covered in
     11.5-43 deg).
   - The active modules cover **8.45 % of 4 pi** (10.3 % with six sectors), not 9.2 %.
   - The simulated photons that make a cluster point to 11.75-44.3 deg (0.1-99.9 % quantiles). This is what the
     detector really sees.

**One factual statement is wrong.** The by-product C2 factor at the edge is confirmed exactly (0.8553 +- 0.0050). But it
does NOT lie "within the ±0.024 allowance": 0.8855 - 0.855 = **0.030**.
- On f(C2 corrected) this is +1.8 %, against E8 = 1.4 %.
- It is also below the WP8 bracket 0.86-0.91.
- It remains a 1.3-sigma-of-the-allowance tension for the analyst to decide. It is not an error of bedge's arithmetic.

**Wording.**
- "Fully reconstructed mostly far below the edge" should be "preferentially far below the edge".
  - The fully reconstructed share is 2.5 times higher far below the edge than at it.
  - But the tagged sample itself lies 27.5 % below -60 MeV, 46.7 % between -60 and -20 MeV, and 25.8 % above -20 MeV.
  - So "the tag and the edge measure different parts of the bump" is too strong. They weight the same bump differently.
- Smaller points are listed in 7.4.

**Verdict: CONFIRMED with corrections.** The answer to the question is right and reproduced: one third of the bump has its
photon in the ECAL, 10-14 % is fully reconstructed, and the limit is the photon energy scale, not the statistics. Before
integration, correct the polar-angle split (Sect. 0, 2.2, Fig. ecalacc_geometry c), the "target centre" wording (Sect. 0,
2.1, Fig. ecalacc_geometry a) and the allowance sentence of the C2 factor (Sect. 0 and 2.5).

### 7.1 Geometry and status (own dump `vecalacc_emcgeo_vae25.C`, `vecalacc_emcstat_vae25.C`; own analysis `vecalacc_geo.py`)

| quantity | package | verifier | agreement |
|---|---|---|---|
| modules / active | 978 / 807 | 978 / 807 | exact |
| sector 0 | all 163 status 0 | all 163 status 0 | exact |
| dead single modules | 1/24, 1/195, 3/160, 3/237, 4/129, 4/157, 5/240, 5/244 | same list | exact |
| status sim vs data | identical (8 runs) | raw 6 x 163 arrays identical (sim + 4 data runs 443850846, 444793253, 445484442, 446694022; 171 zeros each); status via `getPositionFromCell` | exact |
| plane flatness, pitch, face | 0.003 mm, 94 mm, 92 mm | 0.0003 mm, 94.000 mm (grid residual 0.0005 mm), 92.00 mm | exact |
| front planes | 2357-2411 mm | 2357.1-2410.9 mm | exact |
| corner theta from z = 0 | 11.5-47.5 deg (per sector table) | 11.44-47.56 deg; per sector equal to 0.1 deg | exact |
| corner phi per sector | 62.8-117.1, ... | same to 0.1 deg | exact |
| active solid angle from z = 0 (z = -60 mm) | 9.20 % (8.84 %) | 9.16 % (8.81 %) | -0.04 points |
| six sectors / strict 92 mm faces | 11.2 % / - | 11.17 % / 8.78 % | exact |
| active share of the 11.5-47.5 deg band | 60 % | 60.2 % | exact |
| **from the target centre z = -124 mm** | not given | **corners 10.90-45.59 deg; >= 50 % azimuth 11.5-43.0 deg; active 8.45 % of 4 pi (six 10.30 %)** | see 7.0 |

- Vertex: the simulated pion vertex of the C2+C3 rows is uniform in -148.4...-98.6 mm (mean -124.3, rms 14.3 mm). The data
  "pre" vz has quantiles 5/50/95 % = -153 / -125 / -92 mm. The z = -60 mm column of `ecalacc_geometry` is therefore not
  the target either.
- Azimuthal coverage at fixed polar angle (own geometry): active modules 0.55-0.74 of 2 pi, six sectors 0.73-0.89 (Fig.
  vecalacc_theta b). The phrase "about 54 deg of each 60 deg" is the corner extent; at a given theta a working sector
  covers 45-53 deg.
- **Empirical acceptance** (`vecalacc_clusters.py`; standard-set PAT file, every event once):
  - 147 972 neutral candidates; 97 561 with GEANT id 1.
  - Clusters occupy exactly 807 modules, all active. No cluster is in sector 0, in a status-0 module or at an unknown
    (sector, cell).
  - 66 169 target photons with a cluster: true polar angle 0.1/1/50/99/99.9 % = 11.75 / 13.14 / 27.73 / 43.03 / 44.31 deg.
  - 99.16 % of them point to an active module of my geometry (99.58 % to any module).
  - The cluster module equals the predicted module in 81.6 % of the cases, and the same or a direct neighbour in 98.9 %.
  - So the HGeant geometry, the EmcGeomPar geometry and the cell-index mapping agree.

### 7.2 Photon fractions of the accepted bump rows (own rows `vecalacc_rows.C`, own completion `vecalacc_frac.py`)

- Row selection (ownA, "pre", template 1, w x RowFactor(108) x kSimToRaw):
  - C2 1649 rows (42 162.3) and C3 1187 rows (31 056.6), identical to the package; identity check 0 failures;
  - C4 307 rows (7 432.5) and C5 571 246 rows (8 478.3; every 10th used).
- Generator masses from the 598 found photons: 547.452 / 134.976 MeV (package 547.451 / 134.976). 1 row has no solution
  (as in the package; weight about 58, in the 0.40-0.43 slice).
- Lab-frame measure E_gamma d(psi) against the package's flat azimuth in the e+e- frame: chi2/ndf median 1.01 (95 %: 1.40)
  over 200 rows x 20 000 points.

| quantity (pre, per cent of the weight) | package | verifier (own SMASH 1001-1040) | verifier with the package's SMASH file |
|---|---|---|---|
| C2+C3 photon in an active module | 32.7 +- 0.4 | **32.63 +- 0.35** | 32.63 |
| C2 / C3 | 32.9 / 32.3 | 32.88 / 32.29 | 32.87 / 32.30 |
| window B / edge (-15...+5) / S0 | 31.9 / 31.5 / 32.9 | 31.83 / 31.48 / 32.83 | |
| six sectors / strict faces / flat weight | 40.3 / 31.3 / 33.5 | 40.27 / 31.24 / 33.45 | |
| in a module and E_true > 80 / 100 / 130 MeV | 24.5 / 21.6 / 17.7 | 24.5 / 21.6 / 17.7 | |
| the same at the edge | 17.0 / 13.0 / 8.2 | 17.0 / 13.0 / 8.2 | |
| found (any cluster) / + quality cuts | 21.3 / 12.4 | 21.4 / 12.4 | |
| tagged, any photon, m5c window | 12.7 | 12.7 | |
| fully reconstructed (true photon) all / B / edge / S0 | 9.8 / 8.4 / 6.9 / 10.0 | 9.77 +- 0.62 / 8.42 / 6.92 / 9.99 | |
| the same with E_cl x 1.285 | 14.1 | 14.15 (`vecalacc_tag.py`) | |
| <E_gamma> all / in modules [MeV]; E < 100 MeV all / in modules | 130 / 172; 51 / 34 % | 133 / 178; 51 / 35 % | |
| efficiency inside the modules | 0.653 | 0.654 | |
| **polar angle < 11.5 / 11.5-47.5 / > 47.5 deg** | **10.8 / 60.4 / 28.8** | **9.0 / 54.5 / 36.5 (> 90 deg: 9.9)** | DISCREPANCY |
| **18-85 deg** | **76** | **68.9** | DISCREPANCY |
| the same with the package's definition (2-deg bin centres, divided by theta < 90 deg) | | 10.4 / 59.3 / 29.9; 76.9 | explanation confirmed |
| band of the coverage seen from the target (10.9-45.6 deg) / above | | 53.2 / 38.6 | |

- The cause of the discrepancy is in `ecalacc_ana.py`, lines 148-153. The theta histogram of `ecalacc_fill.C` stops at
  90 deg, and the shares are divided by its own sum. So the 9.9 % backward photons drop out of the denominator. The band
  limits are also 2-deg bin centres (12 and 48 deg, not 11.5 and 47.5 deg).
- The same normalisation is in `ecalacc_numbers.json` (theta_shares) and in the annotation and magenta curve of Fig.
  ecalacc_geometry (c). The in-module fractions are not affected.
- Versus M(pi+pi-) and versus the distance to the edge, all entries agree to within 0.2 points. For example, the
  fully reconstructed share by distance bin is 18.6 / 15.9 / 11.1 / 7.7 / 10.9 / 6.7 / 7.2 / 5.1 %.
- **Where the tagged sample sits:** of the true-photon tagged weight,

  | distance d | share of the tagged weight | share of the bump weight |
  |---|---|---|
  | d < -60 MeV | 27.5 % | 16.5 % |
  | -60 <= d < -20 MeV | 46.7 % | 44.9 % |
  | d >= -20 MeV | 25.8 % | 38.2 % |
- Closure of the 598 found photons:
  - all lie on the lab cone (energy residual 5/50/95 % = 0.000 / 0.000 / 0.53 MeV);
  - 99.50 % are in an active module of my geometry (package 99.50 %).
- Found / predicted per bin of the predicted probability: 0.57 / 0.66 / 0.63 / 0.66 / 0.60 (package 0.56 / 0.68 / 0.64 /
  0.68 / 0.63).
- Truncation of the PAT neutral list (at most 5) matters for only 0.1 % of the rows, so it does not explain the 35 % loss.
- Tag at lambda = 1 recomputed from the PAT clusters: agrees with the compact for 2836 of 2836 rows.
  - Any-photon tag ratios at lambda 1.10 / 1.20 / 1.285 = 1.192 / 1.279 / 1.326 (package truth ratios identical).
  - True-photon ratio at 1.285 = 1.449 (package frozen-fit factor 1.464).
- **Peak classes** (the e+e- pair is the measured converted photon in C4 and the virtual photon in C5):

  | class | pair direction in an active module | six sectors | in 11.5-47.5 deg | rows with an ECAL tag (m5c window) |
  |---|---|---|---|---|
  | C4 | 56.8 % | 69.2 % | 74.7 % | 0.78 % |
  | C5 | 54.0 % | 64.8 % | 77.2 % | 1.78 % |

  - For comparison, C2+C3 are tagged in 12.7 % of the rows.
  - So the peak rows make about 2.2 % of the simulated tagged weight in the m5c window, and the tag is specific to the
    bump.
  - The lost bump photon reaches the ECAL less often (32.6 %) than a measured pair direction (54-57 %). The pair is
    accepted in 18-85 deg itself, while the lost photon is spread over the whole sphere.

### 7.3 C2 gamma\* factor, data/expected arithmetic

- Circle average of the QED weight, my lab-frame code: all C2 0.8633 +- 0.0026, B 0.8608, **edge 0.8553 +- 0.0050**, far
  below 0.8724, tagged 0.8734 +- 0.0080. The package has 0.8633 / 0.8608 / 0.8554 / 0.8724 / 0.876. Agreement is exact
  (the tagged value differs by 0.3 sigma).
- **Statement corrected:**
  - 0.8855 - 0.8554 = 0.030. That is 1.25 times the E8 allowance of ±0.024, and outside the WP8 bracket 0.86-0.91.
  - The f(C2 corr.) shift (+1.8 %, 0.418 -> 0.425; I reproduce the composition 0.9386 -> 0.9223 at a C2 step share of
    53.6 %) exceeds E8 on f (1.4 %).
  - The note's "within the E8 uncertainty, but on its edge" should read "just outside the E8 allowance (0.030 vs 0.024)".
- Expected T and ratios, recomputed by hand from the typed inputs:
  - 0.3916 x 7295 = 2857;
  - errors 5.36 % (stat), 8.2 % (syst) and 14.0 % (MC: 10.85 % for f and 8.86 % for T_sim, in quadrature) give
    153 / 234 / 400;
  - 3314 / 2857 = 1.160 +- 0.206, including the data's 105 (stat) and 124 (model);
  - best 4148 x 0.964 x 0.978 = 3910, so 0.848 +- 0.151;
  - C2 corr. 0.4184 x 0.928 x 7295 = 2832 and kappa 0.5602 x 0.711 x 7295 = 2905.
  - All agree.
- "Visible bump yield about 28 600" is not computed in the code.
  - It equals f (0.3916) x the C2+C3 "pre" weight (73 219) = 28 670.
  - The typed constant Y_BUMP_ALLM4 = 26 959 (S0) is unused.
  - T / yield = 11.5 %, so "about 12 %" is fine.
- START: "+75 rows (+2.3 %)" uses the window-B share of START-less rows (7.2 %). With the all-M4 share (6.4 %) it is +46
  rows (+1.4 %). Quote 1.4-2.3 %.
- Not re-derived by me: the data tag fits (T = 3313.7, per slice, START split), the pi0 fits of M(e+e-gamma), and the
  data-side figure histograms.

### 7.4 Figures, tables, units, conventions

All five figures were rendered (PNG) and read. The numbers in them match the note.
- `ecalacc_geometry`:
  - (a) the axis says "seen from the target centre", but it is z = 0;
  - (c) the magenta curve is normalised to theta < 90 deg, and the annotation 11 / 60 / 29 % should be 9.0 / 54.5 /
    36.5 %. The blue band is the coverage from z = 0; from the target it is 10.9-45.6 deg.
- `ecalacc_fractions`, `ecalacc_energy`: correct.
  - Colour convention: magenta is used for "module and E > 100 MeV" and for the photon of all C2+C3, and orange for the
    simulated T1+T2 (also in `ecalacc_geometry`, `ecalacc_closure`).
  - In the Part III palette orange is C2 and magenta is C3. This is minor, but a reader of Part III may read the
    magenta curve as C3.
- `ecalacc_tag`: correct. There is a wide empty gap between panels (a) and (b) (layout).
- `ecalacc_closure`: (a) shows only ±5e-4, while the text quotes |Delta cos| < 1e-3.
- Tables: the three .tex files compile with booktabs (tested in the scratchpad). Their numbers equal the note.
  - `ecalacc_tag.tex` quotes the total error (stat + syst + MC) of T_exp; the caption should say so.
  - `ecalacc_geometry.tex`: the z = -60 mm column should be replaced by z = -124 mm (target centre).
- Units: MeV/c2 and GeV/c2 are consistent. The distance in MeV/c2 uses M4' = 0.99615 M4 and M(pi+pi-) x 0.9952 as in
  defs.h, and I reproduce the same bins.
- Also checked: "pre"/S0, the weight w x RowFactor(108) x 27.0718, and C2+C3 owned by the standard set. All correct.

### 7.5 Verification figure (draft caption)

**Fig. `figs_all/vecalacc_theta.pdf`.**
- **What it shows: one third of the lost bump photons hit an active ECAL module; over the whole sphere 9 % fly below the
  ECAL, 55 % into its polar band and 37 % above it (10 % backward), and the ECAL seen from the real target (z = -124 mm)
  covers 10.9-45.6 deg.**
- (a) Polar angle of the lost photon of the accepted bump rows C2+C3 ("pre", own conditional completion, own SMASH
  sample), all (magenta) and in an active module (blue). Shaded: module corners seen from z = 0; hatched: from the target
  centre.
- (b) Covered share of the azimuth versus polar angle, active modules and six full sectors, from z = 0 and z = -124 mm.
- (c) True polar angle of target photons that made an ECAL cluster (all standard-set events); dashed / dotted: the corner
  range from z = -124 mm and from z = 0.

### 7.6 Files and commands (verifier)

Code in `/home/przygoda/HADES/FAT/FAT/plots/all_report/`:
- `vecalacc_emcgeo_vae25.C`, `vecalacc_emcstat_vae25.C` (vae25, read only);
- `vecalacc_rows.C`, `vecalacc_geo.py`, `vecalacc_frac.py`, `vecalacc_clusters.py`, `vecalacc_tag.py`, `vecalacc_figs.py`.

Figure: `/home/przygoda/HADES/FAT/FAT/docs/latex/figs_all/vecalacc_theta.{pdf,png}`.

Outputs in `/hdd3/przygoda/HADES/FAT_out/all_report/vecalacc/`:
- `vemcgeo_stdout.txt`, `vemcstat_stdout.txt`, `vgeo_stdout.txt`, `vgeo_numbers.json`;
- `vrows.txt`, `vclusters_std.txt`;
- `vsmash_eta_lab.txt` (94 314 etas, md5 370d56ed...);
- `vfrac_stdout.txt`, `vfrac_numbers.json`, `vfrac_rows.npz`;
- `vclusters_stdout.txt`, `vtag_stdout.txt`.

    cd /home/przygoda/HADES/FAT/FAT; V=/hdd3/przygoda/HADES/FAT_out/all_report/vecalacc
    (cat plots/all_report/vecalacc_emcgeo_vae25.C; echo ".q") | ssh -o RemoteCommand=none -o RequestTTY=no vae25 \
      'source /lustre/hades/user/przygoda/pat3/profile.sh >/dev/null 2>&1; cd /tmp; root -l -b 2>&1' > $V/vemcgeo_stdout.txt
    (cat plots/all_report/vecalacc_emcstat_vae25.C; echo ".q") | ssh -o RemoteCommand=none -o RequestTTY=no vae25 \
      'source /lustre/hades/user/przygoda/pat3/profile.sh >/dev/null 2>&1; cd /tmp; root -l -b 2>&1' > $V/vemcstat_stdout.txt
    ssh -o RemoteCommand=none -o RequestTTY=no vae25 'd=/lustre/hades/dstsim/feb22/p4500p/smash/convert/smash_standard/1; \
      for i in $(seq 1001 1040); do awk "NF==9 && \$5==17 {print \$1, \$2, \$3, \$4}" $d/p_p_4500MeV_100000evts_$i.evt; echo "ENDFILE $i"; done' > $V/vsmash_eta_lab_raw.txt
    grep -v ENDFILE $V/vsmash_eta_lab_raw.txt | awk 'NF==4' > $V/vsmash_eta_lab.txt
    root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_vecalacc", kTRUE)' \
      'plots/all_report/vecalacc_rows.C+O(1)'          # (1): also writes vclusters_std.txt; about 40 s, 10 GB
    python3 -I plots/all_report/vecalacc_geo.py > $V/vgeo_stdout.txt
    python3 -I plots/all_report/vecalacc_clusters.py > $V/vclusters_stdout.txt
    python3 -I plots/all_report/vecalacc_frac.py > $V/vfrac_stdout.txt
    python3 -I plots/all_report/vecalacc_tag.py > $V/vtag_stdout.txt
    python3 -I plots/all_report/vecalacc_figs.py

No report .tex, no package file and no existing file was edited (only this section was appended to the note). No git
command was run. On vae25 I only read. At most 2 processes ran at a time, and the CPU package stayed at 39-44 C.
