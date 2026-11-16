# c2recon: which C2 gamma* factor does each estimator need? -- results note

Date 2026-11-15. Branch pp45_pippimepem (worktree /home/przygoda/HADES/FAT/FAT). No report was edited and no git command was
run. No existing file was edited, moved or removed: bedge_ana.C, bedge_common.h, edgeu_core.py, defs.h, SimEngine.h, the
bedge / edgeu / ecalacc outputs and the PAT files were only read. All new files have the prefix c2recon. At most one heavy
process ran at a time; CPU package 39 C.

**Machine-readable outputs** (`/hdd3/przygoda/HADES/FAT_out/all_report/c2recon/`): `c2recon_numbers.json` (everything),
`c2recon_table.txt` (the factor table), `c2recon_edgeu.json`, `c2recon_bedge.txt` (KEY lines), `c2recon_toy.json`,
`c2recon_weights.npz` / `.txt` (one QED weight per simulated C2 row), `c2recon_rows.txt` (the C2 rows with their true
momenta); logs `*.log`; figure `c2recon_factors.{pdf,png}` (not in docs/).

---------------------------------------------------------------------------------------------------
## 0. In plain words

**The question.**
- In the bump (class C2: eta -> pi+pi-pi0, pi0 -> e+e-gamma, photon lost), HGeant lets the virtual photon decay into e+e-
  with a flat angular distribution. QED gives 1 + cos^2 theta* + e sin^2 theta* (e = 4 m_e^2 / M(ee)^2). The simulated C2
  efficiency is therefore too high, and the simulated C2 yield must be multiplied by a factor eps_QED / eps_iso.
- Several values of this factor were in use or quoted: 0.8855 (bedge), 0.881 (edgeu), 0.879 (the verifier's toy),
  0.855 / 0.863 / 0.876 (the ecalacc circle method), 0.875-0.880 and 0.891 (HGeant note, report_sim.tex).
- Which one is right for which measurement?

**What was done.**
1. **Every simulated C2 row got its own QED weight.** w = 3 [1 + c^2 + e (1 - c^2)] / (4 + 2 e), with c = cos theta*.
   - The lost photon is not stored. But the true momenta of the four charged particles are (PAT file).
   - They fix the photon up to one angle: the photon lies on a small circle (the pi0 mass and the eta mass must come out).
   - The weight was averaged along that circle, with the eta momentum spectrum of SMASH as the density. This is the
     conditional mean of w given everything the analysis sees. For any spectrum built from the charged tracks it gives
     the same template as the true w.
   - Three checks:
     - a toy with the true photon: the circle average reproduces the true weight in every subset, to 0.0004 or better;
     - 359 simulated rows in which the ECAL found the photon: exact weight 0.866 against 0.867 from the circle (difference
       -0.001 +- 0.003);
     - the per-row weights agree with the independent ecalacc code (rms difference 2e-6).
2. **The reweighted templates were passed through each estimator's own fit code**, unchanged:
   - the per-slice edge (bedge: bedge_ana.C included as it is, the full evaluation repeated);
   - the u edge (edgeu: edgeu_core.py, the nominal fit reproduced to the last digit);
   - the ECAL tag expectation (ecalacc).
   The factor an estimator needs is the ratio of the simulated step it fits with and without the QED weights.
3. **The factor's statistical error** came from 2000 Poisson replicas of the simulated event groups (and the 100 WP0
   replicas). Its dependence on the density, on the soft-lepton efficiency and on the estimator settings was also tested.

**What came out.**
- **Both edge estimators need about 0.86, not 0.88:**
  - u edge: **0.861 +- 0.007 (sim. stat.) +- 0.004 (model)**; per class E1 0.870, E2 0.856;
  - per-slice edge: **0.860 +- 0.010 +- 0.004**; E1 0.868, E2 0.854.
- **What is applied now:**
  - bedge applies 0.879 / 0.882 per class. On its C2 step this is effectively **0.881** (the "0.8855" quoted for bedge
    is the all-class mean, not what acts on the step).
  - edgeu applies 0.879 / 0.8825, effectively **0.881**.
- **The gap is 0.021 for both estimators.** That is inside the ±0.024 allowance of D35, but its sign is known: the applied
  factor is 2.4 % too high.
- **The ECAL-tagged rows need 0.876 +- 0.009.** The ecalacc value is right as it stands, because it is already computed row
  by row (exact where the photon is known).

**Why the earlier numbers differ.**
- **The toys (0.875-0.883, the verifier's 0.879) use a sharp 60 MeV/c lepton threshold.** The factor depends mainly on the
  momentum of the softer lepton:
  - it is 0.91 at 60-100 MeV/c, 0.87 at 100-150 MeV/c and 0.80-0.82 above 200 MeV/c;
  - in the toy 34 % of the accepted decays have a softer lepton below 80 MeV/c; in the full simulation only 5 %.
  - The toy with a 100 MeV/c threshold gives 0.862, the full-simulation value.
- **The closure 0.891** compares the lepton and standard sets for prompt pi0 Dalitz pairs:
  - it is another pi0 population with another selection;
  - PLUTO and HGeant differ in more than the angle (M(ee) spectrum 2 % harder);
  - its E1/E2 pattern (0.869 / 0.896) is the reverse of every C2 estimate.
  The mean of toy and closure (0.879 / 0.882) inherits both offsets. Most of the gap (0.017 of 0.021) comes from the E2
  class, where the closure (0.896) pulls the applied factor up.
- **The ecalacc numbers (0.863 all M4, 0.855 "at the edge") are window means** of the same per-row weights.
  - The -15...+5 MeV/c2 window catches the rows with the lowest factor.
  - The step estimators also weight the rows between -45 and -15 MeV/c2, and give a negative weight to rows that leak
    above the edge. Their effective factor (0.860-0.861) lies between the edge window (0.855 +- 0.005) and the all-M4
    mean (0.864 +- 0.003).
  - So "at the edge" is not the right definition; the estimator's own weighting is.
- **The +1.8 % on f quoted by the ecalacc verifier** compared 0.855 with 0.8855. The right comparison (0.860 against the
  effective 0.881) gives **+1.2 % (per-slice edge) and +1.4 % (u edge)**.

**What it means.**
- **The bump results change by +1.2 to +1.4 %, far inside their errors (11-15 %):**
  - f_eta(u-edge) 0.420 -> **0.426**; sigma from the bump 0.438 -> **0.444 mb** (with conversions x 0.48: 0.563 ->
    0.573 mb);
  - f_eta(per-slice edge) 0.418 -> **0.423**; sigma 0.437 -> 0.442 mb;
  - visible yields almost unchanged (window B +0.1 %);
  - ECAL data/expected (best reading) 0.85 -> 0.84;
  - the headline sigma (0.429 mb, peak side bands) does not change.
- **Proposal for the analyst (D38, the lead analysis):** replace the constant WP8 class factors by the per-row QED weight
  (effective 0.86) and shrink the allowance from ±0.024 to ±0.004. The reports should not be changed before the ruling.
  Until then they can say in one sentence that the applied factor is 2.4 % high, so f is 1.2-1.4 % low.

---------------------------------------------------------------------------------------------------
## Headline table: the factor eps_QED / eps_iso per estimator and definition

Simulated C2 rows of the standard set, S0, owned in mode A, W = w x RowFactor(108) x 27.0718; per-row weight = circle mean
with SMASH density A. stat = rms of 2000 group-Poisson replicas.

| estimator / definition | all | E1 (< 10 MeV) | E2 (>= 10 MeV) | stat (all) | applied now |
|---|---|---|---|---|---|
| **u-edge estimator (edgeu global fit, linear h, data covariance)** | **0.8607** | 0.8696 (+-0.012) | 0.8560 (+-0.009) | 0.0071 | 0.879 / 0.8825 -> 0.8813 |
| u-edge class fits (E1 and E2 fitted separately) | -- | 0.8694 | 0.8562 | -- | same |
| **per-slice edge estimator (bedge, 4 slices, linear h, data covariance)** | **0.8599** | 0.8682 (+-0.012) | 0.8539 (+-0.016) | 0.0098 | 0.879 / 0.882 -> 0.8807 ("0.8855") |
| **ECAL-tagged C2 rows (pre, m5c 0.49-0.62; exact photon in 175 of 217 rows)** | **0.876** | -- | -- | 0.009 | 0.876 (ecalacc) |
| window mean: all C2 rows, all M4 (S0) | 0.8637 | 0.8651 | 0.8629 | 0.0027 | |
| window mean: 0.30 <= M(pipi) < 0.40, all M4 | 0.8647 | 0.8657 | 0.8641 | -- | |
| window mean: window B (0.4550-0.5225) | 0.8613 | 0.8642 | 0.8597 | 0.0031 | |
| window mean: -15...+5 MeV/c2 from the edge (ecalacc "edge") | 0.8547 | 0.8580 | 0.8530 | 0.0053 | (pre: 0.8554, = ecalacc) |
| window mean: -45...+45 MeV/c2, 0.30-0.40 (the edgeu fit range) | 0.8605 | 0.8653 | 0.8576 | 0.0034 | |
| window mean: far below the edge (u < -60 MeV/c2) | 0.8717 | 0.8554 | 0.8822 | 0.0058 | |
| toy, p_e > 60 MeV/c (WP8-like acceptance), all / edge | 0.8840 / 0.8879 | 0.8922 | 0.8720 | 0.0002 | WP8 toy 0.880 |
| toy, p_e > 100 MeV/c, all / edge | 0.8619 / 0.8617 | 0.8735 | 0.8449 | 0.0003 | bedge edge toy 0.859 |

Fig. `c2recon_factors.pdf`: (a) every source of the factor, (b) the factor against the softer-lepton momentum.

---------------------------------------------------------------------------------------------------
## 1. Method

### 1.1 The weight of one row
- QED: dN/dc ~ 1 + c^2 + e (1 - c^2), e = 4 m_e^2 / M(ee)^2, c = cos theta* (e+ in the gamma* rest frame against the gamma*
  direction). HGeant: dN/dc = const (gdalet.F, WP8 V5).
- Per row the reweighting factor is w = 3 [1 + c^2 + e (1 - c^2)] / (4 + 2 e). Its mean over the isotropic decay at fixed
  M(ee) is 1; the M(ee) spectrum of HGeant (Kroll-Wada) is already right (WP8: 0.998).
- In the gamma* rest frame the gamma* direction is the line of the lost photon, so c = cos(e+, photon) in the e+e- rest frame
  (only c^2 enters).

### 1.2 The lost photon: conditional mean along its circle
- The PAT file has no photon and no pi0 for a bump row. The four charged particles are known (true momenta, read by
  `c2recon_rows.C`, same rows and selection as ecalacc / bedge; identity check against the compact: 0 failures).
- The photon k must satisfy (Q + k)^2 = m_pi0^2 and (P + k)^2 = m_eta^2 (Q = e+e-, P = pi+pi-e+e-; HGeant masses 134.976 and
  547.451 MeV/c2, measured from the found photons: their M(e+e-gamma) and M5 are exactly these values).
  - In the rest frame of Q this fixes |k| and the angle between k and P.
  - One azimuth phi is free. Its measure is flat, and the generator density along it is the SMASH eta invariant density
    f(P + k) (HGeant: eta -> 3pi phase space, isotropic decays).
- **<w>_row = sum_phi f(P + k) w(phi) / sum_phi f(P + k)** (720 points; 2880 points change it by < 3e-6).
  - This is the conditional mean of w given the four charged momenta. So for every quantity built from the charged tracks
    (M4, M(pipi), M(ee), u, acceptance, selection), the template weighted with <w> has the same expectation as the
    template weighted with the true w, and a smaller variance.
- **How much does the circle leave open?** Little. Along the circle w varies only between 0.78 and 0.89 for the median row.
  The four charged tracks almost fix cos theta*.
- Every row has a solution, and no row has zero density on its whole circle.

### 1.3 Checks of the method
1. **Toy closure** (`c2recon_toy.py`, 4e6 decays): etas from SMASH file B; phase-space eta -> 3pi; Kroll-Wada; isotropic
   gamma*; true w known. Acceptance: four tracks 18-85 deg, pions > 0.10 GeV/c, leptons > p_min.
   - The circle average with the density of the INDEPENDENT file A reproduces the true mean in every subset (all, E1, E2,
     near the edge, far below, 0.30-0.40) and for every p_min (0.04-0.15). Largest difference 0.0004 (1.8 sigma).
   - Example, p_min 0.06: true 0.8840, circle 0.8840; near the edge 0.8879 / 0.8879.
2. **Full simulation, the 359 C2 rows with a found photon** (cluster truth): the exact w averages 0.8659. The circle
   restricted to the active ECAL modules gives 0.8671 (difference -0.0012 +- 0.0034); the per-row correlation is 0.90.
   This tests the posterior with the real detector and selection.
3. **Independent code:** the per-row means agree with ecalacc's meanQed for all 1649 rows (rms difference 1.8e-6), and the
   exact weights of the found photons to 6e-8.
4. **Density variants** (all C2 rows, pre):
   - file A 0.8633; file B (other SMASH files) 0.8632; coarse binning 0.8633;
   - KLOE-2 X dependence along the circle 0.8633 (Y is fixed by M(pipi), so the analysis's Y-only weights need no
     change);
   - no density at all (flat, an extreme bound) 0.8663.

### 1.4 The estimators, run with their own code
- **u edge (`c2recon_edgeu.py`).**
  - edgeu_core.py is imported and the nominal configuration of edgeu_ana.py rebuilt (template width 6.7019 MeV/c2).
  - Reproduced exactly: h_data 1305.22 +- 62.27, chi2 27.34/32, h_sim 3107.33 (C2 raw 2020.11, C2 x 0.879/0.8825 1780.31,
    C3 1327.02), f 0.42004.
  - The edgeu C2 rows (1514 in S0, mode A) were joined one to one with the per-row weights.
  - h_sim is linear in the template (fixed data covariance), so h(C2 QED) = sum_rows a_row W_row <w>_row. The sum of
    a_row W_row reproduces 2020.107.
  - Class fits (E1 and E2 fitted separately), the kappa split, windows and estimator variants were redone.
- **Per-slice edge (`c2recon_bedge.C`).**
  - bedge_ana.C is #included unchanged; bedge_ana() is never called and nothing is written into the bedge directory.
  - Reproduced: E_data 1141.8584 +- 59.2243, E_sim 2915.6983, f 0.39162, f_rec 0.41841.
  - The C2 template [nom] was rebuilt from the rows: equal to sim_S.bin bin by bin (max difference 1.1e-7).
  - The QED-weighted template was substituted and the full evaluation repeated:
    - the linear C2 step per class (data covariance) -> the factor;
    - E_sim with the WP6 own weights -> f directly;
    - the bedge formula with these class factors -> f as bedge quotes it.
  - The linear functional of every slice x class fit was exported for the row-level statistics (`c2recon_combine.py`).
- **ECAL tag (`c2recon_combine.py`).** c_tag recomputed with the own weights: exact w where the photon is known, the
  circle mean otherwise (= ecalacc's qtag).

---------------------------------------------------------------------------------------------------
## 2. Results per estimator

### 2.1 u edge (edgeu)
- **c_u = 0.8607 +- 0.0071 (stat)**; E1 0.8696 +- 0.0120, E2 0.8560 +- 0.0092.
  - The class fits give the same: E1 0.8694, E2 0.8562.
  - E2 makes 65.6 % of the isotropic C2 step (the verifier's 0.656).
- Applied now: 0.879 (E1) / 0.8825 (E2), effectively 0.8813. Gap -0.0206: E1 part -0.0032, **E2 part -0.0174**.
- **f_eta(u-edge): 0.42004 -> 0.42574 (+1.36 %)**; uncorrected 0.38995 (unchanged); the C2 share of the step 0.567.
- **Estimator settings.** Event-by-event weights adapt to each setting automatically. A constant factor would need:
  - range ±35 MeV/c2: 0.8555 (-0.0053 +- 0.0053 paired);
  - range ±55: 0.8571;
  - upper end +30 / +60: 0.8590 / 0.8607;
  - template width 5.76 / 7.10 MeV/c2: 0.8603 / 0.8609;
  - r = 0.89 / 1.10: 0.8596 / 0.8603;
  - E(m, 0): 0.8571 (-0.0036 +- 0.0050);
  - pol3: 0.8395 (-0.021 +- 0.012, 1.8 sigma; its coefficients swing more strongly across the edge).
- **Where the estimator puts its weight** (sum of a W over the C2 rows):
  - u -45...-15 MeV/c2: +686 (sum |a| W 1886; <w> 0.864);
  - -15...0: +1620 (<w> 0.851);
  - 0...+15: -277 (resolution leak above the edge, <w> 0.856);
  - +15...+45: -8.
- **Sensitivity to the soft-lepton efficiency** (rows with a softer lepton below 100 MeV/c x 0.8 / x 1.2): -0.0023 /
  +0.0021.
- Other densities: B 0.8609, coarse 0.8607, KLOE-X 0.8607, flat 0.8634.

### 2.2 Per-slice edge (bedge)
- **c_edge = 0.8599 +- 0.0098 (stat)**; E1 0.8682 +- 0.0120, E2 0.8539 +- 0.0155 (WP0 replicas: the same).
- Applied now: 0.879 / 0.882 (the WP8 means; "0.8855" is the all-class mean), effectively **0.8807** on the C2 step.
  - Gap -0.0207: E1 -0.0046, **E2 -0.0161**.
- **f (the bedge formula, as quoted): 0.41841 -> 0.42345 (+1.20 %).**
  - Full re-evaluation with the reweighted template (E_sim with the WP6 own weights):
    - constant class factors 0.41933;
    - QED weights 0.42458 (+1.25 %).
  - The 0.2 % between the formula and the direct fit is a method difference that does not depend on the factor.
- Other constants through the same formula:
  - 0.8855 for both classes: 0.41727;
  - the edge-toy class factors 0.8983 / 0.8628: 0.41908.
- Other densities: B 0.8597, coarse 0.8597, KLOE-X 0.8599, flat 0.8633. Soft-lepton efficiency ±20 %: ∓0.002.

### 2.3 ECAL tag (ecalacc)
- **c_tag = 0.8762 +- 0.0093** (217 tagged C2 rows; 175 with the exact photon). This equals ecalacc's 0.876 +- 0.009:
  that value is already row by row and needs no change.
  - The circle alone gives 0.869 for these rows. A found photon sits in the ECAL, so its position on the circle is not
    typical; the exact weights are the right ones.
- **The expected tag changes only through f:**
  - best reading (energies x 1.285, resolution x 0.964, START x 0.978): expected 3910 -> **3957** with the per-slice f
    (the ecalacc input) or 3978 with the u-edge f;
  - data/expected 0.848 +- 0.150 -> **0.837 +- 0.148** (0.833 +- 0.147 with the u-edge f);
  - the "as is" reading (uncorrected f x uncorrected simulation; 2857, ratio 1.16) does not depend on the factor.

---------------------------------------------------------------------------------------------------
## 3. Why the earlier numbers differ

| source | value | what it measures | why it is not the factor of the edge estimators |
|---|---|---|---|
| WP8 toy | 0.880 (E1 0.889, E2 0.868) | all M4, toy acceptance, leptons > 60 MeV/c (sharp) | too many soft leptons: the factor is 0.91 for a softer lepton of 60-100 MeV/c and 0.80-0.84 above 200 MeV/c; 34 % of the toy rows have it below 80 MeV/c, 5 % of the simulated rows |
| WP8 closure (report_sim 0.891 +- 0.003) | 0.891 (E1 0.869, E2 0.896) | lepton/standard sets, prompt SMASH pi0 Dalitz pairs | another pi0 population and selection; PLUTO and HGeant differ also in M(ee) (2 % harder); E1/E2 pattern reversed |
| bedge "recommended" | 0.8855 (classes 0.879 / 0.882) | mean of toy and closure | inherits both; effective on the C2 step 0.8807 |
| bedge edge toy | 0.8826 (E1 0.898, E2 0.863) | WP8 toy through the bedge fit | same toy acceptance; with leptons > 100 MeV/c it gave 0.859, the full-simulation value |
| edgeu applied | 0.879 / 0.8825 -> 0.8813 | the WP8 means | as bedge |
| vedgeu toy | 0.879 (E1 0.902, E2 0.868) | WP8 generator through a u-edge estimator | toy acceptance |
| HGeant note / Part III | 0.875-0.880 | toys with HADES-like acceptance | toy acceptance |
| ecalacc circle | 0.863 all M4, 0.861 B, 0.855 "edge", 0.872 far, 0.876 tagged | full simulation, the same per-row weights as here | window means; the estimators weight the rows differently (Sect. 2.1) |
| **this note** | **0.861 (u edge), 0.860 (per-slice), 0.876 (tag)** | each estimator's own fit on the reweighted full simulation | -- |

In detail:
1. **Toy acceptance.** In the full simulation the factor also falls faster with the softer-lepton momentum than in the toy
   (0.842 against 0.853 at 150-200 MeV/c, 0.815 against 0.839 at 200-300 MeV/c). So the toy misses more than the momentum
   threshold; the cause was not studied (candidates: the acceptance of the harder lepton and the pair cuts of the
   selection). No sharp-cut toy reproduces this. Only the full simulation does.
2. **E1/E2 mix.** C2 at the edge is mostly E2 (65.6 % of the u step, 57.3 % of the per-slice step), and E2 has the lower
   factor (0.854-0.856). The applied E2 factor (0.882 / 0.8825) is the mean of the toy (0.868) and the closure (0.896), and
   it is the largest single cause of the gap (0.016-0.017 of 0.021).
3. **Definition of "at the edge".** A window mean (-15...+5 MeV/c2: 0.855) is not what a step fit uses.
   - The fit weights the rows just below the edge positively (80 % of h), the rows down to -45 MeV/c2 with mixed signs,
     and the rows that leak above the edge negatively.
   - Its factor (0.860-0.861) therefore lies between the edge window and the all-M4 mean (0.864).
4. **Statistics.**
   - The edge window rests on 364 rows (±0.005).
   - The estimator factors carry ±0.007 (u) and ±0.010 (per-slice), because their coefficients change sign.
   - The 0.855 / 0.861 difference is about one sigma of correlated rows.
   - The toy factors have negligible statistical errors; their problem is the acceptance, not the statistics.

---------------------------------------------------------------------------------------------------
## 4. Recommended values

| use | factor | uncertainty | applied now | comment |
|---|---|---|---|---|
| **u edge (edgeu)** | **0.861** (E1 0.870, E2 0.856) | ±0.007 (sim. stat.) ±0.004 (model) = ±0.008 | 0.881 | best: the per-row weight in the template |
| **per-slice edge (bedge)** | **0.860** (E1 0.868, E2 0.854) | ±0.010 ±0.004 = ±0.011 | 0.881 | same |
| **ECAL tag** | **0.876** | ±0.009 | 0.876 | unchanged |
| windows (any other use) | all M4 0.864; window B 0.861; -15...+5 MeV/c2 0.855 | ±0.003 / ±0.003 / ±0.005 | -- | definition matters at the 1 % level |

- **Model uncertainty ±0.004:**
  - density along the circle: ±0.002 (two thirds of the extreme "no density" shift +0.003; the two SMASH samples agree to
    0.0002);
  - soft-lepton efficiency ±20 % below 100 MeV/c: ±0.002;
  - the posterior check with found photons: ±0.003 (its own precision; centre -0.001).
- **Event-by-event or constant?**
  - **Event-by-event (recommended):** every simulated C2 row is multiplied by its <w> (`c2recon_weights.txt`).
    - The MC replicas then carry the factor's statistical error automatically.
    - Every systematic variant gets its own factor (range ±35: 0.856; pol3: 0.840).
    - Only the model part (±0.004) stays as U10/E8.
  - **Constant:** 0.861 / 0.860 with ±0.008 / ±0.011, and the spread of the variants is not followed.
- On f (C2 share 0.57 / 0.54) the model allowance ±0.004 is ±0.27 % / ±0.25 %. That replaces U10 1.3 % (u) and E8 1.4 %
  (per-slice). The syst totals become 7.2 % (u; was 7.3 %) and 8.2 % (per-slice; was 8.3 %).

---------------------------------------------------------------------------------------------------
## 5. Effects (event-by-event weights instead of the applied class factors)

| quantity | applied (D35) | event-by-event QED weights | change |
|---|---|---|---|
| **f_eta(u-edge)** | 0.420 +- 0.020 +- 0.031 +- 0.029 | **0.426 +- 0.020 +- 0.031 +- 0.030** (syst 7.2 %) | +1.36 % |
| f_eta(u-edge), conversions x 0.48 | 0.540 | 0.549 | +1.75 % |
| u edge, E1 / E2 class fits | 0.409 / 0.433 | 0.411 / 0.444 | +0.4 / +2.5 % |
| u edge, E1/E2 split: a_D, kappa | 0.444, 0.876 | 0.459, 0.835 | within ±0.21 (stat) ±0.35 (MC) |
| **sigma(pp -> eta X) from the bump (u)** | 0.438 +- 0.020 +- 0.032 +- 0.030 +0.042/-0.044 mb | **0.444 +- 0.021 +- 0.032 +- 0.031 +0.043/-0.044 mb** | +1.36 % |
| sigma from the bump, conversions x 0.48 (u) | 0.563 mb | 0.573 mb | +1.75 % |
| double ratio peak/bump DR_u (0.429 / sigma_u) | 0.98 | 0.97 | -1.3 % (error ±0.37) |
| **f_eta(per-slice edge)** (bedge formula) | 0.418 +- 0.022 +- 0.035 +- 0.047 | **0.423 +- 0.023 +- 0.035 +- 0.048** (syst 8.2 %) | +1.20 % |
| per-slice, conversions x 0.48 | 0.560 | 0.569 | +1.62 % |
| per-slice, E1/E2 split: a_D, kappa | 0.551, 0.509 | 0.574, 0.470 | within ±0.25 / ±0.48 |
| sigma from the bump (per-slice) | 0.437 mb | 0.442 mb | +1.20 % |
| visible bump yield, window B (u) | 21 443 rows (3.53 nb) | 21 468 rows (3.53 nb) | +0.1 % |
| visible bump yield, all M4 (u) | 26 955 rows | 27 027 rows | +0.3 % |
| ECAL tag, best reading: expected / data over expected | 3910 / 0.848 +- 0.150 | 3957 / 0.837 +- 0.148 (u-edge f: 3978 / 0.833) | -1.2 % (-1.7 %) |
| ECAL tag, "as is" reading | 2857 / 1.16 | unchanged | 0 |
| sigma(pp -> eta X) headline (peak side bands) 0.429 mb | -- | unchanged (C2 does not enter the peak) | 0 |

- The uncorrected readings (0.390 u edge, 0.392 per-slice) do not change.
- The C2 correction as a whole is now +9.2 % (u: 0.390 -> 0.426) and +8.1 % (per-slice: 0.392 -> 0.423), instead of
  +7.7 % / +6.8 %.
- **Simulated C2 efficiency.** With the factor 0.86 it is about 16 % too high (1/0.861 = 1.16), not "12-14 %"; that range
  came from the toy and the closure.
- **D31 statement.** "The bump would give 0.56 mb with conversions x 0.48" becomes 0.57 mb (u edge).

---------------------------------------------------------------------------------------------------
## 6. Proposal D38 for the analyst (the lead analysis; not applied in the reports)

**D38 (proposed): the HGeant gamma* correction of class C2 is applied row by row.**
- Every simulated C2 row is weighted with the conditional mean of the QED weight along its photon circle, computed from its
  true charged momenta (`c2recon_weights.txt`). This replaces the constant WP8 class factors 0.879 / 0.882 (bedge) and
  0.879 / 0.8825 (edgeu).
- The effective factors are 0.861 (u edge) and 0.860 (per-slice edge). The allowance U10/E8 shrinks from ±0.024 to ±0.004.
- It changes f_eta and sigma from the bump by +1.4 % (u edge) and +1.2 % (per-slice edge):
  - 0.420 -> 0.426 and 0.438 -> 0.444 mb (u);
  - 0.418 -> 0.423 and 0.437 -> 0.442 mb (per-slice).
  These shifts are well inside the quoted errors, and the visible yields hardly change. ECAL data/expected goes
  0.85 -> 0.84.
- **Reasons:**
  1. the applied factors come from a toy with a sharp lepton threshold and from a closure of another pi0 sample;
  2. the full simulation itself gives the factor once the lost photon is integrated out;
  3. the method closes on a toy to < 0.05 % and is confirmed by the rows with a known photon (-0.1 +- 0.3 %).

**If ruled yes:**
- In Part III, rerun edgeu_ana.py / bedge_ana.C with the per-row weights, or apply the numbers of Sect. 5 directly.
- Then change in the text:
  - "0.875-0.891" -> 0.86;
  - "12-14 % too high" -> about 16 %;
  - "0.8855, checked at the edge within 1 %" -> per-row weights, effective 0.860-0.861;
  - U10 / E8 -> ±0.004.
- Part I keeps 0.891 as the closure of prompt pi0 pairs. Its later-finding box can add "for the eta bump rows the full
  simulation gives 0.86".

**If ruled no:**
- Keep D35. Add one sentence: "the full simulation, reweighted row by row, needs 0.861 (u edge) / 0.860 (per-slice edge)
  instead of the applied 0.881; f would be 1.4 % / 1.2 % higher, inside the U10/E8 allowance of ±0.024".
- Correct the ecalacc verifier's "+1.8 %" to +1.2-1.4 %.

---------------------------------------------------------------------------------------------------
## 7. Caveats

1. **The density along the circle is the SMASH eta spectrum (marginal).** If the selection or the trigger depended on the
   other particles of the event in a way that correlates with the eta momentum, the posterior would be slightly
   different. The extreme bound (no density) moves the factor by +0.003; the found-photon check (-0.001 +- 0.003) shows
   no bias.
2. **Leading-order QED angle only**, as in WP8 (no radiative corrections to the gamma* angle). The pi0 Dalitz M(ee)
   spectrum of HGeant is right (WP8).
3. **The simulated lepton efficiency at low momentum** sets the factor (Fig. (b)). A ±20 % error below 100 MeV/c moves it
   by ±0.002. A larger mismodelling of the soft-lepton efficiency would also bias the efficiency itself, which is not a
   question of this factor.
4. **Only C2 is corrected here.** The same HGeant defect affects any other pi0 Dalitz pair of the standard set:
   - it does not enter h_sim of the u edge (C2 + C3 only);
   - in the per-slice E_sim (T12) the rest (G_*) is 1.1 % of the step, so the effect is < 0.2 %.
   The eta Dalitz closure of Part I (0.879) could be checked in the same way (photon on a sphere instead of a circle); this
   was not done.
5. **Simulation statistics:** 1514 C2 rows in S0. The factor's error (±0.007 / ±0.010) is part of the MC error of f when
   the weights are applied row by row.

---------------------------------------------------------------------------------------------------
## 8. Exact commands (from /home/przygoda/HADES/FAT/FAT)

    O=/hdd3/przygoda/HADES/FAT_out/all_report/c2recon
    root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/c2recon/aclic", kTRUE)' \
         'plots/all_report/c2recon_rows.C+O'   > $O/c2recon_rows_stdout.txt 2>&1     # 27 s, 10 GB: C2 rows + PAT truth
    python3 -I plots/all_report/c2recon_weights.py                                     # 2 s: per-row QED weights (5 densities)
    python3 -I plots/all_report/c2recon_toy.py                                         # 27 s: closure toy (4e6 decays)
    python3 -I plots/all_report/c2recon_edgeu.py                                       # 10 s: u-edge factor, f, classes, yields
    root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/c2recon/aclic", kTRUE)' \
         'plots/all_report/c2recon_bedge.C+O'  > $O/c2recon_bedge_stdout.txt 2>&1     # 10 s: bedge factor, f (bedge code)
    python3 -I plots/all_report/c2recon_combine.py                                     # bedge stat, sigma, ECAL, table
    python3 -I plots/all_report/c2recon_figs.py                                        # c2recon_factors.{pdf,png}

**Code** (`plots/all_report/`, all new):
- `c2recon_rows.C`: row pass (SimEngine, the ecalacc PAT reader, bedge bins);
- `c2recon_weights.py`: circle mean and exact weights, density variants, comparison with ecalacc;
- `c2recon_toy.py`: closure toy;
- `c2recon_edgeu.py`: the u-edge estimator via edgeu_core;
- `c2recon_bedge.C`: the per-slice estimator via bedge_ana.C (#included);
- `c2recon_combine.py`, `c2recon_figs.py`.

**Read, not changed:**
- bedge_ana.C, bedge_common.h, edgeu_core.py;
- the bedge arrays, the edgeu columns and JSON, ecalacc rows_c23.bin / ecalacc_numbers.json;
- the SMASH eta lists of ecalacc (A) and vecalacc (B);
- ppee_std_all.root.

---------------------------------------------------------------------------------------------------
## Verification (independent)

Verifier, 2026-11-15 (evening). Own code with the prefix vc2recon; no c2recon macro, weight file or output was used as an
input. Frozen fit code was called (edgeu_core.py; bedge_ana.C #included), but every reweighted C2 template was built by
the verifier from its own rows and weights. No report, no existing macro and no existing output was changed. Outputs:
`/hdd3/przygoda/HADES/FAT_out/all_report/vc2recon/`. At most two heavy processes ran at a time; CPU package 39-58 C.

**Verdict: CONFIRMED.** The verifier's factors agree with c2recon to 0.0004 or better, well within the ±0.004 test band:

| quantity | c2recon | verifier (own code) | difference |
|---|---|---|---|
| u-edge factor (all-class fit, data covariance) | 0.8607 ± 0.0071 | **0.8606 ± 0.0070** | -0.0001 |
| u edge, E1 / E2 | 0.8696 / 0.8560 | 0.8696 / 0.8559 | 0.0000 / -0.0001 |
| u edge, separate class fits E1 / E2 | 0.8694 / 0.8562 | 0.8694 / 0.8561 | 0.0000 / -0.0001 |
| per-slice edge factor (bedge class fits) | 0.8599 ± 0.0098 | **0.8596 ± 0.0097** | -0.0003 |
| per-slice, E1 / E2 | 0.8682 / 0.8539 | 0.8681 / 0.8533 | -0.0001 / -0.0006 |
| ECAL-tagged C2 rows (pre, m5c 0.49-0.62; 217 rows) | 0.876 ± 0.009 | **0.8764 ± 0.0091** | +0.0004 |
| factor applied now: u edge / per-slice | 0.8813 / 0.8807 | 0.8813 / 0.8807 | 0 |
| f_eta(u-edge): applied -> per-row | 0.42004 -> 0.42574 (+1.36 %) | 0.42004 -> **0.42577 (+1.36 %)** | +0.00003 |
| sigma from the bump (u) | 0.438 -> 0.444 mb | 0.4383 -> **0.4443 mb** | 0 |
| f_eta(per-slice), bedge formula | 0.41841 -> 0.42345 (+1.20 %) | 0.41841 -> **0.42354 (+1.23 %)** | +0.00009 |
| f_eta(per-slice), own-weight E_sim refit | 0.41933 -> 0.42458 | 0.41933 -> 0.42467 | +0.00009 |
| sigma from the bump (per-slice) | 0.437 -> 0.442 mb | 0.4370 -> **0.4423 mb** | 0 |
| window mean, all M4 (S0, mode A) | 0.8637 | 0.8636 | -0.0001 |

**In short:** the per-row QED weight, the photon circle, the HGeant masses and both propagations are right. D38's
numbers (0.861 / 0.860, ±0.004; f and sigma +1.4 % / +1.2 %) stand as they are.

### V1. What was done (own chain)
1. **Rows** (`vc2recon_dump.C`): the SimEngine compact read with the filter and the order of edgeu_dump.C (2 995 797 rows,
   identical count), so that row j is row j of the edgeu columns. The 1649 C2 rows (all from the standard set) got their
   PAT truth from an own reader of ppee_std_all.root (0 identity failures). Join with the edgeu columns: masses equal to
   5e-9 GeV, weights to 5e-8, bits equal.
2. **Photon circle in the LAB frame** (`vc2recon_w.py`), not in the e+e- rest frame as c2recon.
   - The photon k = |k|(1, n) obeys 2Q.k = m_pi0^2 - Q^2 and 2P.k = m_eta^2 - P^2. So n lies on a cone around
     A P_vec - B Q_vec, and |k| = A / (2 (Q0 - Q_vec.n)).
   - Along the cone azimuth the measure is |k| dphi (own derivation from d^3k/|k| and the two delta functions).
     Numerical check against a flat azimuth in the e+e- rest frame: max difference of the row means 6e-6. A flat lab
     azimuth (wrong measure) would shift the mean by -0.005.
   - Prior: the SMASH invariant eta density at P + k. It comes from an **own extraction of 40 further converted SMASH
     files** (3001-3040, 94 863 etas; c2recon used files 1-60 and 1001-1040), as a smoothed (y, pT^2) histogram.
   - Every row has a solution; no row has zero density on its circle.
3. **HGeant masses measured from the found photons**: m_pi0 = 134.9764 MeV (spread 5e-5), m_eta = 547.4507 MeV (spread
   0.03), as c2recon. With the PDG masses the factors move by +0.0001.
4. **Own photon match**, kinematic instead of the geninfo/vertex rule: a photon cluster whose truth gives M(e+e-gamma) =
   m_pi0 and M5 = m_eta within 0.05 MeV. It finds 354 rows (c2recon 359), all with parent pi0.
5. **Propagation** (`vc2recon_fit.py`, `vc2recon_bedge.C`).
   - u edge: the nominal fit is reproduced to every digit (sigma_T 6.7019 MeV; h_data 1305.22 ± 62.27; chi2 27.34/32;
     h(C2) 2020.107, h(C3) 1327.021). The own C2 template equals hist_sim to 5e-13. The applied factors give 1780.309.
   - per-slice edge: the own iso C2 template equals sim_S.bin [nom][C2] to 5e-6 in every bin. The template is
     substituted (SS - T_iso + T_QED) and the full bedge evaluation repeated. Row-level replicas use the bedge unit
     responses (417 non-zero coefficients). They reproduce the substituted refit (0.85963) exactly.
   - stat: 2000 group-Poisson replicas (own seed).

### V2. theta* definition, angle model, M(ee) dependence
- **Definition.** cos(e+, photon) in the e+e- rest frame equals the angle of the e+ in the gamma* rest frame against the
  gamma* direction of the pi0 rest frame, reached by the two-step boost lab -> pi0 -> gamma*. Max difference 1.4e-10 for
  the 354 found photons; the e- gives the same |c| (5e-10).
  - A wrong axis (the lab direction of the e+e- pair) would give 0.840 instead of 0.869 for the found photons.
  - So the definition matters at the 3 % level, and c2recon uses the right one.
- **HGeant angle, read from the source** (hgeant2-7.4/geant/gdalet.F, lines 81-84).
  - YINV = beta (2 RNDM - 1) and E2 = (m0/4) [x + 1 + YINV (1 - x)], with x = M(ee)^2 / m_pi0^2 and beta the lepton
    velocity in the gamma* frame.
  - With the boost of the gamma*, E2 = (m0/4) [1 + x + beta cos theta* (1 - x)]. So YINV = beta cos theta*: cos theta* is
    uniform (isotropic), as WP8 and c2recon state.
- **QED weight.** 1 + c^2 + e (1 - c^2) is normalised to mean 1 at fixed M(ee) ((4 + 2e)/3); it is isotropic at threshold
  (e -> 1). The e term matters: without it the factors would be 0.8505 (u) and 0.8477 (per-slice).
- **Form factor.** HGeant samples M(ee) from Kroll-Wada x the VDM rho pole (uhbook.F, histogram 9995, vdm = .true.,
  m_rho = 770 MeV, i.e. an effective slope 0.0307). Replacing it by the PDG slope a = 0.0332 ± 0.0029 changes the effective
  C2 template by 1.0002 (u) and 1.0000 (per-slice), at most 1.0004. The form factor does not enter the angle.
- **M(ee) dependence of the u factor** (true M(ee), share of the step in brackets): 0-5 MeV 0.911 (0.20), 5-10 MeV 0.820
  (0.16), 10-20 MeV 0.854 (0.19), 20-40 MeV 0.850 (0.20), above 40 MeV 0.860 (0.24).

### V3. Is the circle mean the right weight? Three tests
1. **Closure toy** (`vc2recon_toy.py`, own generator).
   - Generation: etas from SMASH list A; eta -> 3 pi phase space; Kroll-Wada; isotropic decays. The true w comes from the
     two-step definition.
   - Posterior: density from the own list C. Cuts on the charged tracks only (18-85 deg, pions > 100 MeV/c, leptons >
     p_min).
   - Two fresh samples of 8e6 decays (267 643 and 267 813 accepted, p_min 100 MeV/c): true minus circle mean +0.0001 ±
     0.0001 for all decays. In every subset (E1, E2, near and far from the edge, soft-lepton bins) it is 0.0003 or less,
     except one subset at +0.0006 ± 0.0002.
   - A first sample of 3e6 decays gave -0.0005 ± 0.0002 (E2 -0.0013 ± 0.0004). The same events with the density of the
     generating list itself, or on a finer grid, gave the same values. So that offset is a fluctuation of that sample,
     not a density effect.
   - With p_min 60 MeV/c: true 0.8838, circle 0.8840 (c2recon's toy: 0.8840 / 0.8840).
2. **Found photons with the exact weight.**
   - 354 rows: exact <w> 0.8693. The circle without an ECAL condition gives 0.8645 for the same rows: a found photon
     sits at a non-typical place on its circle.
   - The circle restricted to the photon's polar range (12.6-43.1 deg, E_gamma > 24 MeV; no module map) gives 0.8691:
     difference exact - circle **+0.0005 ± 0.0035**, correlation 0.86 (c2recon, with the module map: -0.0012 ± 0.0034).
   - PIT of the exact w in its restricted posterior: 79/74/74/64/62 in five bins, chi2 3.0 for 4 degrees of freedom.
   - Caveat: along the circle w varies little (median posterior rms 0.014). So this test checks the kinematics more than
     the isotropic prior. The prior itself is fixed by the gdalet.F source (V2).
3. **Density and model variants** (u / per-slice factor):
   - lists A, B, C: 0.8606 / 0.8596, all within 0.0001;
   - coarse grid, no smoothing: within 0.0001 (per-slice +0.0002);
   - KLOE-2 X dependence along the circle: +0.0000;
   - flat density (extreme bound): +0.0027 / +0.0038;
   - soft-lepton efficiency (softer lepton < 100 MeV/c) x 0.8 / x 1.2: -0.0023 / +0.0021 (u), -0.0021 / +0.0020
     (per-slice).
   These support the model allowance of ±0.004.

### V4. Points c2recon handled correctly that could have gone wrong
- **Exact weights only where the photon conditions the selection.**
  - For the edge templates c2recon uses the circle mean for every row. Taking the exact w wherever a photon was found,
    and the unconditional circle elsewhere, would bias the u factor by +0.0012 (per-slice +0.0001). The reason is that
    finding the photon is correlated with w (V3.2).
  - For the ECAL tag the rows are selected by the found photon, so the exact w is the right weight there: 0.8764
    (circle alone 0.8691).
- **The factor depends on the estimator settings** (the constant factor that each variant would need; c2recon in
  brackets):
  - range ±35 MeV/c2: 0.8555 (0.8555); ±55: 0.8569 (0.8571); upper end +30 / +60: 0.8589 / 0.8606 (0.8590 / 0.8607);
  - template width 5.76 / 7.10 MeV/c2: 0.8602 / 0.8608 (0.8603 / 0.8609);
  - E(m, 0): 0.8568 (0.8571); pol3: 0.8395 (0.8395).
  - So event-by-event weights, as D38 rules, are the right implementation.
- **Where the u estimator puts its weight** (sum a W over the C2 rows): u -45...-15: +686; -15...0: +1620; 0...+15: -277;
  +15...+45: -8 (c2recon: +686, +1620, -277, -8).
- **Window means** (for reference): -45...+45 MeV/c2 in u, 0.30-0.40: 0.8603 (c2recon 0.8605, same definition).
  - The verifier's -15...+5 MeV/c2 window, taken in u (with the M(ee) term, 0.30-0.40), gives 0.8566.
  - c2recon's 0.8547 is the ecalacc "edge" window, which is defined differently, so the two numbers cannot be
    compared one to one.
  - Neither window mean is used by any estimator.

### V5. Small differences, none of them relevant
- Found photons: 354 (own kinematic match) against 359 (geninfo/vertex rule); 10 photon clusters close to within 2 MeV
  in M(e+e-gamma) only. Tagged rows with a found photon: 173 against 175.
- The per-slice E2 factor is 0.0006 lower (0.8533 against 0.8539); all classes 0.8596 against 0.8599.
  - It is not the density: with c2recon's list A the verifier also gets 0.8532 (E2) and 0.8596 (all).
  - The cause was not found. It is 0.04 of the statistical error (±0.015) and 0.15 of the ±0.004 allowance.
  - The u factors agree to 1e-4.

### V6. Files (all new)
- Code (`plots/all_report/`): `vc2recon_dump.C`, `vc2recon_w.py`, `vc2recon_toy.py`, `vc2recon_fit.py`,
  `vc2recon_bedge.C`, `vc2recon_extra.py`.
- Outputs (`/hdd3/przygoda/HADES/FAT_out/all_report/vc2recon/`):
  - `vc2recon_rows.txt` (rows with PAT truth), `vc2_smash_eta_C.txt` (own SMASH list, md5 1a4a84fa...);
  - `vc2_w.{npz,json,log}`, `vc2_toy*.{json,log}`, `vc2_fit_edgeu.{json,log}`, `vc2_bedge.txt`, `vc2_fit_bedge.{json,log}`,
    `vc2_extra.{json,log}`.
- Commands (from /home/przygoda/HADES/FAT/FAT):

      root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/vc2recon/aclic", kTRUE)' 'plots/all_report/vc2recon_dump.C+O'
      python3 -I plots/all_report/vc2recon_w.py
      python3 -I plots/all_report/vc2recon_toy.py ; python3 -I plots/all_report/vc2recon_toy.py _seed2 ; python3 -I plots/all_report/vc2recon_toy.py _seed3
      python3 -I plots/all_report/vc2recon_fit.py edgeu
      root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/vc2recon/aclic", kTRUE)' 'plots/all_report/vc2recon_bedge.C+O'
      python3 -I plots/all_report/vc2recon_fit.py bedge
      python3 -I plots/all_report/vc2recon_extra.py
- The copy `docs/latex/notes_all/c2recon_results.md` was not changed by the verifier (it still equals the note before this
  section); backup of this note: `integ4/backup/vc2recon/c2recon_results.md.pre_vc2recon`.
