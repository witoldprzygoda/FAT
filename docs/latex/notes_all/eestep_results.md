# eestep: where the step at u = 0 in the like-sign e-e- rows comes from -- results note

> **STATUS (final check, 15 November 2026, 17:50): independently verified, CONFIRMED WITH CORRECTIONS; read the section
> "Verification (independent)" at the end first; the copy in docs/latex/notes_all/ is identical to this file.** Corrections used in the reports: e-e- dominance also from Compton/knock-on
> electrons (negative only; SMASH e+e+/e-e- 0.08-0.09, 0.11 without them); the pair bookkeeping (expected k 0.96-1.02);
> 2.7 sigma is statistical (2.4-2.9 sigma with the MC error of the 17 simulated rows); k_cross 1.05 +- 0.39 (about +-0.5;
> independent 0.96 +- 0.48); flip limit < 1.6 % (1.9 % for either charge); excess over SMASH 2.3-3.5 times (2.3-3.7 sigma,
> 11-14 simulated rows), not "at least 3.5 times, about 3 sigma"; about 180 cross pairs in the lower side band (255 is the range
> 0.5225-0.5336).

Date 2026-11-15. Branch pp45_pippimepem (worktree /home/przygoda/HADES/FAT/FAT). No git command was run; no report and no
existing file was edited. All new files carry the prefix eestep_ (macros and scripts in plots/all_report/, outputs in
/hdd3/przygoda/HADES/FAT_out/all_report/eestep/). Inputs were only read: the data compact rows (scan_v2), the PAT files
/hdd3/przygoda/HADES/files/exp/gen4/lep_<day>.root (local copies of the production inputs; vae25 was not needed), the
simulation compact rows and the edgeu column files. At most one heavy process ran at a time; CPU package 38-42 C.

---------------------------------------------------------------------------------------------------
## 0. In plain words

**Verdict: the mechanism is found, and it is physics, not a detector artefact.**

- **What the step rows are.** The e-e- rows of the step are pairs of electrons from the *two different photons* of the pi0 in
  eta -> pi+ pi- pi0:
  - one e- from the Dalitz photon (gamma*) and one e- from the conversion of the other photon,
  - or one e- from each of two converted photons.
  The analysis calls such pairs "cross-photon pairs". In SMASH they are the class T4 ("different sources"). The analysis leaves
  them to the CB (they are never in the simulated templates).
- **Why they sit exactly at u = 0.** The lost particles of such a pair are again pi0 decay products. A short calculation shows
  that, for a given M(e e), the largest possible M4 is *the same* E(M(pipi), M(ee)) as for the bump. A toy check of both
  topologies confirms it: the true u never exceeds 0 (largest values -0.0055 and -0.016 MeV/c2).
- **Why they appear at large M(ee).** The two electrons come from two photons with a wide angle between them. Near the limit
  their mass is large: in the toy 113-121 MeV/c2 (median). In the data the e-e- step sits at the same place:
  - 57.0 +- 9.1 rows (6.3 sigma) for M(ee) >= 90 MeV/c2;
  - 13.0 +- 9.1 rows for 50-90 MeV/c2;
  - 3.6 +- 8.9 rows (none) below 50 MeV/c2, where 90 % of the bump lives.
  Opening angle: the step is mostly in pairs of 20-60 deg (none below 10 deg).
- **Why e-e- and hardly e+e+.** These pairs contain soft leptons. HADES accepts soft e- far better than soft e+, and the
  simulation does the same: for these pairs SMASH gives e+e+/e-e- = 0.08. Data, M(ee) >= 90 MeV/c2: e+e+ 7.3 +- 4.3 rows, a
  ratio 0.13 +- 0.08. So the e+e+ rows have the step too, at the size expected. "e+e+ has no step" (2.1 +- 6.5 over all M(ee))
  was a statement about small numbers.
- **The CB does with them what it is designed to do.** Before acceptance, each such event has one e-e-, one e+e+ and two e+e-
  cross pairs.
  The e+e- ones are in N. CB = 2 sqrt(N++ N--) subtracts them. **This was tested in the data**, in M(ee) >= 90 MeV/c2, where
  the bump is small:
  - the step of N is 52.6 +- 10.9 rows. That is the SMASH bump x f (15.0) plus the CB step (35.8);
  - N - CB = 16.8 +- 13.9 against 15.0 expected (z = +0.13);
  - if the like-sign step had no e+e- partners (an artefact), N - CB would be -20.8. The data are 2.7 sigma away from that;
  - the scale of the cross pairs in N relative to the CB rule is k_cross = 1.05 +- 0.39 (1 = the CB rule).
- **The charge-flip reading is ruled out.** A Dalitz e+ reconstructed as e- would keep its momentum. It would appear at the
  M(ee) of the bump: below 50 MeV/c2 there is no e-e- step (3.6 +- 8.9 against 1210 rows of bump step in N).
  - Flip rate r < 1.2 % (90 % CL);
  - simulation: 0.01-0.04 % per true pair;
  - no wrong-sign copies of the e+ (e+e- rows with an e-e- twin: step -4.1 +- 7.4);
  - the step rows have normal track quality (all quality variables flat, p = 0.14-0.99);
  - no sector, no period and no day dependence (sector shares 0.30-0.47, day p = 0.26).
- **Effects.**
  - f_eta(u-edge) and f_eta(bedge): the nominal values are right. The earlier caveat "f could be higher by up to about 8 %" is
    withdrawn:
    - a flip could raise f by less than +1.6 % (90 % CL);
    - k_cross instead of 1 moves f(u-edge) by -0.1 +- 1.2 % (M(ee) >= 90 MeV/c2 value) or +0.6 +- 1.3 % (>= 70);
    - both are inside U5 (2.6 %). For bedge: inside its E4 variant (+1.3 %).
  - Y_peak: no effect, because the cross pairs reach the peak window only for M(ee) >= 102-120 MeV/c2 and are subtracted
    with their partners. Only in the artefact reading, which the data disfavour at 2.7 sigma, would Y_peak move: -0.7 % (P
    alone) or +3.0 % (with the side bands at R' = 1).
  - X17: no effect. For M(ee) <= 24.5 MeV/c2 the limit is E <= 0.5217 GeV/c2, 12 MeV/c2 (1.8 sigma) below the peak window.
    Below 50 MeV/c2 the e-e- step is zero anyway.
- **One open point (it changes no number).** SMASH x f predicts fewer like-sign cross pairs at the edge than the data show.
  - SMASH x f: at most 21 +- 7 rows per 2.5 MeV/c2. Data: 73.5 +- 15.7. Ratio >= 3.5, about 3 sigma.
  - The SMASH sample is very small: 10 simulated e-e- rows near the edge.
  - Possible reasons (not tested): conversion material that gives single leptons but no reconstructed pairs; the close-pair
    question of D32.
  - It does not affect the CB argument. The CB comes from the data, and the data show partners at the CB-rule level
    (k_cross).
- **Proposed handling (D39, for the analyst to rule on).**
  - Record the origin.
  - Keep all headline numbers and U5 unchanged.
  - Replace the "up to +8 %" caveat by the flip limit (< +1.6 %, 90 % CL) and the k_cross variant (-0.1 +- 1.2 %).
  - Details in Sect. 5.

---------------------------------------------------------------------------------------------------
## 1. What was done

1. **Data dump with track details and event links** (`eestep_dump.C`).
   - The rows: the compact rows of edgeu (cor_ok == 1, pre, 0.26 <= M(pipi) < 0.46, 0.30 <= M4 < 0.85).
     - every like-sign row: 37 835 e+e+ and 127 656 e-e-;
     - the e+e- rows of the edge region: 0.30-0.40, |u| < 100 MeV/c2, 385 800 rows.
   - Each row was joined to its PAT entry: the file lep_<day>.root, the entry pat_entry of its own tree.
     - Check: for all 551 291 rows the compact momenta equal the PAT momenta exactly.
     - Added per lepton: sector, segment chi2, RK chi2, META and RICH match, ring index, beta, ...
     - Added per event: lpos_mult, lneg_mult, hpos_mult, hneg_mult, START.
   - Links between the three PAT trees of one event.
     - Event key: (runnumber, event, start_time, eVertClust_z, totalmult). The event counter alone is not unique (it restarts
       per sub-file); the two reconstructed floats make the key unique.
     - Check: every key is exactly one block of adjacent entries in every tree (0 exceptions in 29 days x 3 trees).
     - Tracks across trees: equal (theta, phi, p).
     - A "twin" of an e-e- row is an e+e- entry of the same event with the same pi+, pi- and one of the two e-.
2. **Simulation dump with truth and twins** (`eestep_simdump.C`).
   - All three sets, the same selection as edgeu, the same key logic (generator event + vertex floats; 0 exceptions).
   - The pion-pair class of the like-sign rows was taken from the edgeu file simls_pcls.i8; the row order was checked identical.
   - Used only to identify the mechanism, never for a CB or k number (rule of the simulation study).
3. **Analysis** (`eestep_ana.py`, `eestep_core.py`). Estimators:
   - the edgeu ones: pol2 + h T(u) on -45..+45 MeV/c2, sigma 6.70 MeV/c2, S0, 0.30 <= M(pipi) < 0.40;
   - like-sign rows: the edgeu null-test fit;
   - N - CB: the nominal edgeu fit (CB per slice in 10 MeV/c2 cells of u, peak term a_P K(u) fixed, h_sim = the same linear
     estimator applied to SMASH C2 + C3 with the C2 gamma* factors).
   - Reproduced exactly: e-e- 73.5 +- 15.7, e+e+ 2.1 +- 6.5, h_data 1305.2, h(CB) 39.2, f 0.4200.
   - "Step share" = h / (rows per 2.5 MeV/c2 bin in -15..0 MeV/c2).
4. **Kinematic toy** (`eestep_kintoy.py`): the limit of cross-photon pairs (Sect. 2.4).
5. **Figures** (`eestep_figs.py`): eestep_spectra, eestep_bands, eestep_align (pdf + png in the output directory).

---------------------------------------------------------------------------------------------------
## 2. Results

### 2.1 Event content (Sect. 2 of eestep_ana.log)
In short: the step rows are events with two e- and no e+ at all. They are not copies of an e+e- bump row.

| e-e- rows, S0, 0.30-0.40 | step h [rows / 2.5 MeV/c2] | rows (-45..45) |
|---|---|---|
| all | 73.5 +- 15.7 | 8964 |
| no e+e- entry in the event | 66.4 +- 14.5 | 7645 |
| e+e- twin in the event (same pi+, pi-, one e-) | 6.8 +- 6.0 | 1319 |
| no positive lepton in the event (lpos_mult 0) | 65.4 +- 14.4 | 7622 |
| lneg_mult >= 3 | 0.6 +- 2.8 | 299 |
| the two e- share a RICH ring | -2.7 +- 4.7 | 812 |
| isBest == 1 / == 0 | 52.6 +- 14.2 / 20.9 +- 6.5 | 7402 / 1562 |

- 85 % of all e-e- rows near the edge have no positive lepton in the event. The step rows are not different in this.
- e+e- rows of the edge region that have an e-e- twin (a negative copy of the e+): step -4.1 +- 7.4, i.e. a clone rate of
  -0.003 +- 0.006 of the bump step. With an e+e+ twin: -2.6 +- 5.6.

### 2.2 The step per M(ee) band (Fig. eestep_bands; Sect. 3 of the log)
In short: the like-sign step lives at M(ee) >= 90 MeV/c2. There the e+e- rows carry the matching step, as the CB rule assumes.

| M(ee) [MeV/c2] | e+e+ | e-e- | N | CB | N - CB | SMASH x f |
|---|---|---|---|---|---|---|
| 0-10 | 1.3 +- 1.3 | -5.4 +- 4.6 | 683.6 | 2.6 | 681.0 +- 43.5 | 698.7 |
| 10-30 | -1.7 +- 1.4 | 4.2 +- 5.8 | 422.7 | -1.1 | 423.8 +- 33.9 | 363.6 |
| 30-50 | -0.4 +- 1.9 | 4.7 +- 5.0 | 103.3 | -0.5 | 103.8 +- 17.1 | 113.2 |
| 50-70 | -2.6 +- 2.6 | 10.1 +- 6.1 | 52.3 | -0.1 | 52.4 +- 14.1 | 78.1 |
| 70-90 | -2.0 +- 3.0 | 3.0 +- 6.7 | 30.0 | 2.8 | 27.2 +- 12.4 | 37.8 |
| 90-110 | 7.6 +- 3.2 | 28.8 +- 6.7 | 31.1 | 22.9 | 8.2 +- 10.3 | 10.9 |
| 110-135 | -0.2 +- 2.9 | 27.3 +- 6.1 | 19.7 | 12.8 | 6.9 +- 9.2 | 3.5 |
| **0-50** | -0.7 +- 2.7 | **3.6 +- 8.9** | 1209.7 | 0.2 | 1209.4 +- 57.7 | 1176.0 |
| **90-135** | **7.3 +- 4.3** | **57.0 +- 9.1** | 52.6 +- 10.9 | 35.8 +- 8.6 | **16.8 +- 13.9** | **15.0** |
| 70-135 | 5.1 +- 5.3 | 60.1 +- 11.3 | 83.8 +- 15.2 | 40.0 +- 10.7 | 43.8 +- 18.6 | 52.4 |

(steps in rows per 2.5 MeV/c2; f = 0.420; the CB of each band comes from that band's own like-sign rows.)

- **The test of the two readings (M(ee) >= 90 MeV/c2).** The model is: step of N = f h_sim + k_cross h(CB).
  - CB rule right (the e+e- sample contains the cross pairs): N - CB = 16.8 +- 13.9 against 15.0, z = +0.13;
    k_cross = 1.05 +- 0.39.
  - Artefact (the like-sign step has no e+e- partners): expected N - CB = 15.0 - 35.8 = -20.8; the data are z = 2.7 away.
  - With M(ee) >= 70 MeV/c2: z = -0.46 against z = 1.69; k_cross = 0.79 +- 0.44.
- Opening angle of the e-e- step rows: 0-10 deg: -8.1 +- 7.1; 10-20: 9.1 +- 6.0; 20-35: 25.1 +- 8.7; 35-60: 40.9 +- 8.5;
  60-180: 5.9 +- 3.3. The bump step of N is at 0-10 deg (1170 of 1349).

### 2.3 The step follows E(M(pipi), M(ee)) (Fig. eestep_align; Sect. 4)
In short: the step is tied to the edge of eta -> pi+pi-pi0, including its M(ee) term. It is not a fixed structure in M4.

- e-e- rows, M(ee) >= 70 MeV/c2. Step size and delta chi2 (step against none) for five variables:

  | variable | step | delta chi2 |
  |---|---|---|
  | u = M4 - E(m, q) | 60.1 +- 11.3 | 28.3 |
  | M4 - E(m, 0) (no M(ee) term) | -7.7 +- 10.8 | 0.5 |
  | M4 - 0.535 (fixed M4) | -- | 10.3 (chi2 60/32) |
  | u - 10 MeV | -- | 3.5 |
  | u + 10 MeV | -- | 19.8 |

- Free position of the step: best at u0 = -2.5 MeV/c2. In chi2, u0 = 0 is 2.0 higher, u0 = -5 is 1.6 higher, u0 = +2.5 is
  6.9 higher and u0 = +5 is 13.4 higher.
- Width: not constrained (4-13 MeV/c2 all give chi2 35-42).
- Per M(pipi) slice the step grows toward the high slice: 8.7 +- 6.7, 17.8 +- 6.4, 30.9 +- 6.4.
- Control region M(pipi) 0.42-0.46, where eta -> pi+pi-pi0 is not possible, with the edgeu pseudo edge:
  - M(ee) >= 70 MeV/c2: -2.0 +- 3.1 (none);
  - all M(ee): -12.4 +- 5.7 (no positive step).

### 2.4 Kinematic limit of cross-photon pairs (eestep_kintoy.py, eestep_kintoy.log)
In short: both topologies have exactly the bump edge, and the pairs near the edge are heavy.

- Let photon 1 (energy E1) give a lepton with fraction x1, and photon 2 (E2) give a lepton with fraction x2. Then:
  - M(ee)^2 = x1 x2 m_pi0^2;
  - M4^2 = m_eta^2 - 2 m_eta [(1 - x1) E1 + (1 - x2) E2] + (1 - x1)(1 - x2) m_pi0^2 (eta frame).
- At fixed M(ee) the maximum is at x2 = 1 and E1 = (E0 - p0)/2. It gives m_eta^2 - m_eta (E0 - p0)(1 - M(ee)^2/m_pi0^2). This
  is the bump limit of defs.h (before the scale 0.99274).
- Toy, 400 000 decays per topology, with the true masses:
  - (A) two conversions: largest u_true -0.0055 MeV/c2. Pairs within 5 MeV/c2 of the limit: M(ee) median 121 MeV/c2 (10-90 %:
    84-131 MeV/c2).
  - (B) Dalitz e- + conversion e-: largest u_true -0.016 MeV/c2. M(ee) median 113 MeV/c2 (49-130 MeV/c2).

### 2.5 Track and row properties (Sect. 5 of the log; e-e-, M(ee) >= 70 MeV/c2)
In short: the step rows are ordinary e-e- rows. No sign of bad or flipped tracks.

- Step share in tertiles. A flat share means that the step rows look like the other e-e- rows.
  - Flat: RK chi2 (p 0.94), inner-segment chi2 (0.47), outer-segment chi2 (0.74), META match (0.92), RICH match (0.14), ring
    pads (0.91), beta (0.99), p_min (0.99), p_max (0.40), systems (0.55), pi+/pi- momenta, vertex z, isBest (0.54), day (0.26).
  - Same sector or not: 0.39 / 0.43.
  - Rows with a sector-2 lepton: 0.39 (others 0.40).
  - Sector by sector: 0.30-0.47.
- Not flat: M(ee) (p 0.01: 0.04 below 89 MeV/c2, 0.52-0.57 above), as expected from Sect. 2.2.
- Polar angle (p 0.06): larger at large theta, as expected for wide pairs.

### 2.6 Charge flip (Sect. 6)
In short: a flip would put the step at small M(ee); there is none.

- M(ee) < 50 MeV/c2: e-e- step 3.6 +- 8.9 against 1209.7 rows of bump step in N.
  - Flip rate r = 0.003 +- 0.007, r < 1.24 % (90 % CL).
- At the limit, f would be higher by less than +1.6 %:
  - 1.2 % of lost bump rows;
  - 0.35 % of CB over-subtraction.
- Simulation (true e+ and e- of one source in a like-sign row, per true e+e- pair, weighted): std 0.04 %, lep 0.01 %, eta
  0.008 %.
  - Per true e+ (lep set, any partner, by reconstructed momentum): 1.7 % below 100 MeV/c, 0.46 % at 100-200, 0.10 % at
    200-300, <= 0.06 % above.

### 2.7 Peak window and X17 window (Sect. 7)
In short: the cross pairs reach the peak window only at very large M(ee), and the X17 band not at all.

- E(m, q) reaches the peak window (0.53364) only for M(ee) >= 101.5 / 104.0 / 110.5 / 119.5 MeV/c2 at M(pipi) = 0.28 / 0.30 /
  0.35 / 0.40.
- Cross-pair share of the e-e- rows just below the edge (M(ee) >= 90 MeV/c2): 0.56 +- 0.09. Applied to the e-e- rows with
  u <= 0 in each window:

  | window | e-e- rows with u <= 0 | cross pairs | CB rows carried |
  |---|---|---|---|
  | SBlo | 454 | about 255 | 184 |
  | P | 89 | about 50 | 37 |
  | SBhi | 0 | 0 | 0 |

- With partners (the CB rule, supported at 2.7 sigma) Y_peak is unchanged.
- Without partners (the artefact reading), Y_peak would move:
  - -0.7 % from P alone;
  - +3.0 % with the side bands at R' = 1.
- X17: for M(ee) <= 24.5 MeV/c2 the limit is E <= 0.5217, 12.0 MeV/c2 (1.8 sigma) below P. The X17P window holds 207 e-e- and
  13 e+e+ rows (S0, M(pipi) 0.26-0.46; D19: 223 / 14); none can be a cross pair below its edge. The e-e- step at 10-30 MeV/c2
  is 4.2 +- 5.8.

### 2.8 Simulation of the cross pairs (Sect. 8)
In short: SMASH has the same kind of like-sign pairs, but too few rows to measure their step, and fewer than the data.

- Selection: std set (the only set with eta -> pi+pi-pi0), pions of one eta (pcls 2), both leptons from pi0 photons.
- e-e-: 110 of 733 rows (weighted 3184 of 19124 expected raw rows); median M(ee) 98 MeV/c2, median OA 37 deg.
- e+e+: 9 of 115 rows (weighted 262); weighted ratio e+e+/e-e- 0.08.
- Near the edge: 10 e-e- rows below u = 0 and 5 above. The step fit is not possible (-9 +- 35 weighted).
- SMASH x f gives at most 21.1 +- 7.1 rows per 2.5 MeV/c2 below the edge. Data step: 73.5 +- 15.7. Ratio >= 3.5, z = 3.0.

---------------------------------------------------------------------------------------------------
## 3. Answers to the six questions of the brief

1. **Per-row properties.** The step rows are normal e-e- rows: track quality, momenta and sectors are all flat (Sect. 2.5).
   They differ only in M(ee) (>= 90 MeV/c2) and opening angle (20-60 deg). The unlike-sign bump rows are at small M(ee) and
   at 0-10 deg.
   - "Which lepton would be the flipped e+": neither. Both e- of a step row are ordinary electrons from the two photons.
2. **Is there an e+ in the event?** In 85 % of the e-e- rows (step rows alike) there is no positive lepton at all and no
   e+e- entry. The step is in these events (66.4 +- 14.5). Events with an e+e- twin show none (6.8 +- 6.0). The links use
   event keys with reconstructed floats; each key is one adjacent block.
3. **Charge-flip rate.**
   - Simulation: 0.01-0.04 % per true pair; 1.7 % per e+ below 100 MeV/c, <= 0.1 % above 200 MeV/c.
   - Data: r < 1.2 % (90 % CL) for the bump e+.
   - No sector, period or day effect: sector shares 0.30-0.47; sector 2 0.39; days p 0.26; edgeu periods p 0.18.
4. **Other mechanisms.**
   - Ghosts or duplicates: no step with a shared RICH ring (-2.7 +- 4.7), in twins (6.8 +- 6.0) or with three or more e-
     (0.6 +- 2.8).
   - Pi- taken as e-: the step rows have electron-like beta and RICH quality (flat), and a pion could not follow E(m, q).
   - **Conversion of the second photon: yes. This is the mechanism** (Dalitz + conversion, or two conversions).
5. **Why e-e- and not e+e+.** The acceptance for the soft leptons of these pairs is much smaller for e+:
   - SMASH e+e+/e-e- = 0.08;
   - data step ratio 0.13 +- 0.08 at M(ee) >= 90 MeV/c2;
   - e+e+ step there 7.3 +- 4.3.
6. **Peak.** The cross pairs reach P only at M(ee) >= 102-120 MeV/c2; most of them sit in SBlo (Sect. 2.7).
   - In the CB-rule reading nothing changes, for Y_peak or for the peak CB.
   - The X17 CB is not touched.

---------------------------------------------------------------------------------------------------
## 4. Effects on the measurements

| quantity | nominal (CB rule, supported) | variant k_cross (data) | flip (90 % CL limit) | artefact reading (disfavoured, 2.7 sigma) |
|---|---|---|---|---|
| f_eta(u-edge) 0.420 | unchanged | -0.1 +- 1.2 % (M(ee) >= 90); +0.6 +- 1.3 % (>= 70) | < +1.6 % | +3.0 % (the CB step of 39.2 rows removed) |
| f_eta(bedge) 0.418 | unchanged | inside its E4 variant (+1.3 %) | < +1.6 % | not computed (same order) |
| Y_peak 5088 | unchanged | -- | -- | -0.7 % (P only) / +3.0 % (with side bands) |
| X17 (r, BR x BR) | unchanged | -- | -- | none (out of reach) |

---------------------------------------------------------------------------------------------------
## 5. Proposed handling (D39, for the analyst to rule on; no headline number changes)

1. **Record the origin.**
   - The step of the like-sign e-e- rows at u = 0 is made of cross-photon lepton pairs of eta -> pi+pi-pi0 (Dalitz e- +
     conversion e-, or two conversion e-).
   - They have the same kinematic limit as the bump and lie at M(ee) >= 90 MeV/c2.
   - They are seen mostly as e-e- because of the acceptance for soft e+ (SMASH 0.08, data 0.13 +- 0.08).
   - It is not a charge flip: r < 1.2 % (90 % CL).
2. **Keep the nominal CB** (2 sqrt(N++ N--)) and **U5 (2.6 %) as they are.** The data test in M(ee) >= 90 MeV/c2 supports the
   CB rule for these pairs (k_cross = 1.05 +- 0.39; the artefact reading is disfavoured at 2.7 sigma).
3. **Withdraw the caveat** "f could be higher by up to about 8 % (wrong charge sign)" in all three parts. Replace it by:
   - the flip limit (< +1.6 % on f, 90 % CL);
   - the k_cross variant (-0.1 +- 1.2 %), quoted as a check covered by U5, not as a new term.
4. **Correct the wording** "the e+e+ rows have no step". Proposed text: "the e+e+ rows have the same step at the size expected
   from the e+/e- acceptance (7.3 +- 4.3 rows at M(ee) >= 90 MeV/c2)". Also the "30 % against 36 % share" argument: the e-e-
   step lives at M(ee) >= 90 MeV/c2, where the bump is small, so the coincidence of shares means nothing.
5. **Optional, for the open point:** SMASH under-predicts the like-sign cross pairs at the edge (ratio >= 3.5, about 3 sigma,
   10 simulated rows). A larger cross-pair sample (std set or a dedicated eta -> pi+pi-pi0 sample with GEANT conversions)
   would tell whether this is statistics, material, or the close-pair effect of D32. It affects no number of the reports.

Where it belongs:
- Part II (data only): the like-sign null test, Sect. 2.1-2.3, 2.5-2.7.
- Part III: the CB test with SMASH x f (Sect. 2.2), the simulation comparison (2.8) and the decisions.
- Part I: the simulated flip rates (2.6) and the cross-pair population (2.8), if the analyst wants them there.

---------------------------------------------------------------------------------------------------
## 6. Limits of this study

- The k_cross test uses SMASH x f for the bump in M(ee) >= 90 MeV/c2 (15.0 rows). A 20 % change of that prediction moves z by
  about 0.2.
- The share 0.56 +- 0.09 used for the peak windows is measured at the edge in 0.30-0.40. Deeper below the edge the share is not
  known; the window numbers are estimates.
- The simulated cross-pair sample is small (110 e-e- rows); no SMASH number is used for any correction.
- Large M(ee) bands (>= 50 MeV/c2) are outside the blinded X17 region; the X17P like-sign counts were only counted.

---------------------------------------------------------------------------------------------------
## 7. Files and how to run (from /home/przygoda/HADES/FAT/FAT)

    B='gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_eestep", kTRUE)'
    O=/hdd3/przygoda/HADES/FAT_out/all_report/eestep
    root -l -b -q -e "$B" 'plots/all_report/eestep_dump.C+O'    > $O/eestep_dump_stdout.txt 2>&1    # 75 s, data + PAT
    root -l -b -q -e "$B" 'plots/all_report/eestep_simdump.C+O' > $O/eestep_simdump_stdout.txt 2>&1 # 33 s, simulation
    python3 -I plots/all_report/eestep_ana.py                      # 4 s -> eestep_ana.log, eestep_ana.json, eestep_figdata.npz
    python3 -I plots/all_report/eestep_kintoy.py > $O/eestep_kintoy.log
    python3 -I plots/all_report/eestep_figs.py                     # eestep_spectra, eestep_bands, eestep_align (.pdf, .png)

- Code: plots/all_report/eestep_dump.C, eestep_simdump.C, eestep_core.py, eestep_ana.py, eestep_kintoy.py, eestep_figs.py
  (eestep_simls_look.py: first look at the edgeu simulated like-sign columns).
- Outputs in /hdd3/przygoda/HADES/FAT_out/all_report/eestep/:
  - tables ls_*, usr_*, sus_*, sls_* (.f32; column lists in eestep_layout.txt, eestep_simlayout.txt);
  - logs eestep_dump.log, eestep_simdump.log, eestep_ana.log, eestep_kintoy.log;
  - eestep_ana.json (md5 5d7cf70f13f554a54ccad100bb02a37b);
  - figures eestep_spectra, eestep_bands, eestep_align (.pdf, .png).
- Figures for the reports must be copied by the writers. Each caption should start with "What it shows":
  - eestep_bands: the like-sign step sits at M(ee) >= 90 MeV/c2, and there N - CB follows SMASH x f;
  - eestep_spectra: e-e- u spectra per M(ee) band, with no step below 50 MeV/c2 and the step at >= 90 MeV/c2;
  - eestep_align: the step is tied to E(m, q) and sits at u = 0.

## 8. Self-checks
- edgeu numbers reproduced exactly: e-e- 73.5 +- 15.7, e+e+ 2.1 +- 6.5, h_data 1305.2 (1305.3), f 0.4200, CB step 39.2.
- PAT join: 551 291 rows, momenta identical (max |diff| 0). Event keys: one adjacent block each (0 exceptions, data and
  simulation).
- Simulated like-sign rows: 238 040 in both dumps, identical M4 column (row order equal).

---------------------------------------------------------------------------------------------------
## Verification (independent)

Date 2026-11-15, 16:10. Independent verifier (prefix veestep). No eestep_* macro, output or number was used as an input; the
note above was read only to know the claims. Own code: a new data dump with PAT joins, an own step estimator (Poisson
likelihood), an own CB construction, an own bootstrap, an own kinematic toy and an own reading of the simulated like-sign
rows. Inputs: the data compact rows (scan_v2) with the frozen selection of defs.h, the PAT files lep_<day>.root (read only),
the edgeu simulation column files (sim_*, simls_*). No report, no existing file and no git state was touched; the backup of
this note is in final/integ3/backup/veestep/. CPU package 39-43 C, never more than one heavy process.

### V0. In plain words

**Verdict: CONFIRMED WITH CORRECTIONS.** The mechanism is right, the numbers of the data side are reproduced, and nothing in the
headline results changes. Four statements need softer or more complete wording.

- **The mechanism holds.** The e-e- step at u = 0 is made of lepton pairs from the two different photons of the pi0 in
  eta -> pi+pi-pi0. I proved the edge in general (see V1): any pair of pi0 decay products, with the rest lost, can never pass
  E(M(pipi), M(ee)). My own toy confirms it for every pair type (true u never above 0 in 600 000 decays per topology).
- **The toy also predicts where the step sits.** For these pairs the density falls in the last few MeV below the edge (reaching the
  limit needs a second condition: the partner of the other lepton must carry no energy). After the resolution this looks like a
  step at u = -2.5 MeV/c2. The data put the step exactly there; for the bump classes the toy puts it at 0.
- **The data side is reproduced row by row.** Same rows as edgeu (1 910 080, identical columns), the same steps (e-e- 73.5 +- 15.6,
  e+e+ 2.1 +- 6.5), the same M(ee) and opening-angle pattern, the same event content (85 % of the e-e- rows have no positive
  lepton), the same flip limit.
- **The CB rule handles these pairs.** At M(ee) >= 90 MeV/c2 the unlike-sign rows have the matching step: N steps by 52.7 +- 10.4
  rows, where the bump alone would give 14.7. k_cross = 0.96 +- 0.48 (eestep: 1.05 +- 0.39). A simple acceptance model gives an
  expected k between 0.96 and 1.02. So f_eta is right to about 0.1 %.
- **Corrections (wording only, no number of the reports changes):**
  1. **Why mostly e-e-.** It is not only the soft-e+ acceptance. In SMASH about one third (34-36 %) of the simulated e-e- pairs of
     this kind contain a Compton electron or a knock-on (delta) electron of the photon. Such electrons are always negative and have
     no positive partner. Without them SMASH gives e+e+/e-e- = 0.11 (0.08-0.09 with them, as quoted).
  2. **"Each decay gives one e-e-, one e+e+ and two e+e-."** This holds for Dalitz and conversion legs only. A Compton or delta
     electron gives one e-e- and one e+e- pair, no e+e+. The CB rule still works (expected k 0.96-1.02).
  3. **Size of the CB test.** k_cross has an error of about +-0.5 when the fluctuation of the CB step and the simulation error are
     included (only 17 simulated bump rows at M(ee) >= 90 MeV/c2: the prediction 14.7 has a 46 % MC error). The artefact reading is
     disfavoured at 2.4-2.9 sigma with these errors (2.6-3.7 sigma statistical only). The note's 2.7 sigma is a fair middle value. Its
     sentence "a 20 % change of that prediction moves z by about 0.2" understates the MC error.
  4. **The open point (excess over SMASH).** The data have about 2.3 to 3.5 times more such like-sign pairs at the edge than SMASH x
     f (2.3-3.7 sigma, depending on the estimator; 11-14 simulated rows). "At least 3.5 times" is the upper end. Neither suggested
     explanation fits on its own:
     - conversion material would need 1.8-2.4 times more photon interactions than HGeant; the bump's own conversion scale
       (kappa = 0.88 +- 0.41) disfavours that at 2.3-3.7 sigma;
     - the close-pair question of D32 can raise the expectation only through the normalisation, by at most x 1.2;
     - pairs merged into one track (a whole photon in one lepton) are not seen: the step rows have no excess of large rings or
       of double energy loss.
     The excess stays open (statistics of the simulation, or more Compton/delta electrons than in HGeant).
- **Small corrections.**
  - A flip of either charge (also e- -> e+, which the CB amplifies 3.25 times) raises f by less than 1.9 % (90 % CL); the e+ -> e-
    flip alone by less than 1.6 % (as quoted).
  - The lower peak side band SBlo (0.5250-0.5336) holds about 180 cross pairs, not 255; 255 is the wider range 0.5225-0.5336.
- **Y_peak and X17 are not touched** (confirmed). Cross pairs reach the peak window only at M(ee) >= 101.5-119.5 MeV/c2, about 50
  of them, carrying 37 CB rows; with the expected k the residual is below 1.5 rows (0.03 % of Y_peak). Below M(ee) = 24.5 MeV/c2
  their edge is at most 0.5211 (0.5223 at M(pipi) = 0.26), 11-12.5 MeV/c2 below the peak window. None of the 207 e-e- rows of the
  X17 window lies within 12 MeV/c2 above its edge.

### V1. The edge of cross-photon pairs (a)
In short: the edge is exactly E(M(pipi), M(ee)) for any pair of pi0 decay products; the toy confirms it and predicts the step
position seen in the data.

- **General proof (pi0 rest frame).** Let Q be the observed pair (Q^2 = q^2) and L = P_pi0 - Q the rest (lost particles, L^2 >= 0,
  L^0 >= 0). M4^2 = (P_pipi + Q)^2 = m_pipi^2 + q^2 + 2 P_pipi.Q. At fixed q this is largest when Q points against P_pipi and
  |Q*| is as large as allowed. L^2 >= 0 means |Q*| <= m_pi0 - Q0*. The maximum is at |Q*| = m_pi0 - Q0*, i.e. Q0* =
  (m_pi0^2 + q^2)/(2 m_pi0): the two-body decay pi0 -> (q) + massless. That is exactly the bump limit
  m_eta^2 - m_eta (E0 - p0)(1 - q^2/m_pi0^2). So every pair whose leptons come from the eta's pi0 (Dalitz + conversion, two
  conversions, Compton electrons up to the electron mass) has the same edge; pairs with a lepton from elsewhere do not.
- **Toy** (veestep_toy.py; own generator: 3-body phase space, Kroll-Wada, 1 + cos^2, Bethe-Heitler splitting; lab frame
  pp -> pp eta at 4.5 GeV; crude acceptance 18-85 deg, p > 80 MeV/c; charge-symmetric).

  | pair type (600 000 decays each) | largest true u [MeV/c2] |
  |---|---|
  | C2 bump (Dalitz pair, photon lost) | -0.0003 |
  | C3 bump (conversion pair, other photon lost) | -0.0000 |
  | Dalitz + conversion: e-e- / e+e+ / e+e- / e-e+ | -0.108 / -0.075 / -0.027 / -0.082 |
  | two conversions: e-e- / e+e+ / e+e- / e-e+ | -0.030 / -0.009 / -0.030 / -0.043 |

- **Shape near the edge** (accepted, density per 5 MeV/c2 from -40 to 0, relative to -40..-35):
  - bump C2: 1.00 ... 1.37 1.49 (rises to the edge: a sharp step);
  - cross pairs (Dalitz + conversion, e-e-): 1.00 ... 1.04 0.87 0.38 (falls in the last 5-10 MeV/c2);
  - with sigma 6.7 MeV/c2 and the data estimator: best step position +0.0 MeV/c2 for the bump classes, **-2.5 MeV/c2 for the cross
    pairs**; data: best -2.5 MeV/c2 (veestep_ana, Sect. 6).
- **Pairs near the edge** (u > -5 MeV/c2, accepted): M(ee) 10/50/90 % = 107/125/131 (Dalitz + conversion), 113/128/132 MeV/c2 (two
  conversions); opening angle 21/37/61 deg; 90 % of these rows have no detected e+ in the event.
- **M(ee) split of the step.** The toy puts 76-88 % of the cross-pair step at M(ee) >= 110 MeV/c2. The data have 28.8 (90-110) and
  27.3 (110-135). The full simulation (12 SMASH rows within -15..+5 MeV/c2 of the edge) has 0.22 / 0.38 / 0.40 below 90 / 90-110 /
  >= 110, the same split as the data (0.23 / 0.40 / 0.38). The difference is the crude toy acceptance, not a problem of the
  mechanism.

### V2. Steps per M(ee) band and opening angle (b)
In short: every number of Sect. 2.2 and 2.3 is reproduced; my CB construction differs, so N - CB differs within its error.

Estimator: Poisson likelihood, pol2 + step (sigma 6.7 MeV/c2), -45..+45 MeV/c2, S0, 0.30-0.40; N - CB with a linear estimator,
errors from 300 event-group bootstrap replicas. CB per slice and band in 10 MeV/c2 cells, shared over the fine bins by the like-sign
counts (edgeu interpolates linearly).

| M(ee) [MeV/c2] | e+e+ | e-e- | N (peak term removed) | CB | N - CB | SMASH x f |
|---|---|---|---|---|---|---|
| 0-50 | -0.7 +- 2.8 | 3.6 +- 8.9 | 1209.7 | -1.2 | 1210.8 +- 56.3 | 1174.6 |
| 50-70 | -2.6 +- 2.6 | 10.1 +- 6.1 | 52.3 | 1.3 | 51.0 +- 13.8 | 77.2 |
| 70-90 | -2.0 +- 2.9 | 3.0 +- 6.7 | 30.1 | 0.6 | 29.5 +- 13.0 | 37.3 |
| 90-110 | 7.6 +- 3.3 | 28.8 +- 6.8 | 31.0 | 25.6 | 5.4 +- 11.4 | 10.7 |
| 110-135 | -0.2 +- 2.8 | 27.3 +- 6.0 | 21.7 | 14.4 | 7.3 +- 10.0 | 4.0 |
| **90-135** | **7.3 +- 4.3** | **57.0 +- 9.0** | **52.7 +- 10.4** | **39.6 +- 9.8** | **13.1 +- 14.6** | **14.7 +- 6.7 (MC)** |

- e+e+/e-e- step ratio at M(ee) >= 90: 0.127 +- 0.079 (eestep 0.13 +- 0.08).
- Opening angle (e-e-, all M(ee)): 0-10 deg -8.1 +- 7.0; 10-20 9.1 +- 5.8; 20-35 25.1 +- 8.7; 35-60 40.8 +- 8.5; 60-180 5.9 +- 3.4.
  At M(ee) >= 90: 0-20 deg 0.0 +- 1.4; 20-60 deg 50.7; 60-180 deg 5.9 +- 3.3. The e+e+ rows show no significant step in any class.
- Alignment (M(ee) >= 70): step in u 60.1 +- 11.2 (delta(-2lnL) 28.5); in M4 - E(m, 0) -7.7 +- 10.7 (0.5). Free position: best
  -2.5 MeV/c2 (0: +2.1; +2.5: +7.0; +5: +13.7). Width 4 / 6.7 / 10 / 13 MeV/c2: -2lnL 35.2 / 36.6 / 39.6 / 41.8.
- M(pipi) slices (M(ee) >= 70): 8.7 +- 6.6, 17.8 +- 6.5, 30.9 +- 6.2. isBest == 1 / != 1: 52.6 +- 14.2 / 20.9 +- 6.5.
- Control region 0.42-0.46 with the pseudo edge E(m - 0.10, q): -12.4 +- 5.8 (all M(ee)); 1.1 +- 2.4 (50-70).

### V3. Event content (c)
In short: confirmed with my own links; the events were found by key blocks of adjacent entries, and every key was one block.

- Join: all 165 491 like-sign rows found at their PAT entry (lepton and pion angles equal; phi taken modulo 360 deg). The e+e-
  and e-e- trees of each day were cut into blocks of adjacent entries with the same (runnumber, event, start_time, eVertReco_z,
  totalmult); a key seen again after other keys would have been excluded: 0 cases in 29 days.
- e-e- rows (S0, 0.30-0.40), step all M(ee) | M(ee) >= 90:
  - lpos_mult == 0: 65.4 +- 14.3 | 54.4 +- 8.5; lpos_mult >= 1: 7.8 +- 6.1 | 2.6 +- 3.1;
  - no e+e- entry in the event: 66.4 +- 14.4; with an e+e- twin (same pi+, pi-, one e-): 6.8 +- 6.0 (every event with an e+e-
    entry has a twin);
  - lneg_mult >= 3: 0.6 +- 2.8; two e- with a shared ring: -2.7 +- 4.8; group with three or more leptons: -0.3 +- 2.8.
- Share of the e-e- rows near the edge with no positive lepton: 0.850 (toy: 0.90 for cross pairs near the edge).
- **Merged conversion pairs** (a new test, connected with D32): a conversion pair reconstructed as ONE lepton would carry the whole
  photon and sit at the edge with a sharp step; it would have about twice the ring and the MDC energy-loss signal. Step share by
  tertile of the larger value of the two e-: ring pads 0.59 / 0.53 / 0.30; MDC dE/dx 0.24 / 0.79 / 0.36. The top tertile is not
  enriched: no sign of merged pairs. (rich_amp and rich_avg_ringcharge are 0 in the PAT files and could not be used.)

### V4. Charge flip (d)
In short: confirmed; the limit for both charges is slightly larger.

- M(ee) < 50: e-e- step 3.6 +- 8.9, e+e+ -0.7 +- 2.8, bump step of N - CB - peak 1210.8. e+ -> e-: r = 0.0030 +- 0.0073, r < 1.23 %
  (90 % CL); effect on f < +1.6 % (lost row + CB rise sqrt(N++/N--) = 0.308).
- e- -> e+: r = -0.0006 +- 0.0023. An e+e+ excess raises the CB by sqrt(N--/N++) = 3.25 times its size, so it matters more for f.
  Both together: f higher by +0.1 +- 1.4 %, < +1.9 % (90 % CL).
- Negative copies of the e+ (e+e- rows with an e-e- entry of the same pi+, pi-, e- and a second e- within 3 deg of the e+):
  2 102 of 466 088 rows; step -4.8 +- 3.1 rows (-0.0035 of the N step). With any e-e- entry in the event: -4.1 +- 7.4.

### V5. The CB test at large M(ee) (e)
In short: the unlike-sign partners are there; k_cross is consistent with 1 and with the expectation; the stated significance is
statistical only.

- **Direct test.** At M(ee) >= 90 the step of N (peak term removed) is 52.7 +- 10.4 rows; the SMASH bump x f gives 14.7 +- 6.7 (MC:
  17 simulated rows) +- 1.6 (f). The excess, 38 rows, is the e+e- cross pairs. Artefact reading (no partners): 3.7 sigma statistical,
  2.9 sigma with the MC error and a 20 % peak-term scale error.
- **The note's form** (N - CB against pred - h_CB): 2.6 sigma statistical (eestep 2.7); 2.4 sigma with the MC error. It is
  conservative, because the error of N - CB contains the CB fluctuation, which cancels in this test.
- **k_cross** = (h_N - f h_sim)/h_CB = 0.96 +- 0.48 (bootstrap of the ratio); +- 0.52 with the MC and peak-term errors (eestep 1.05 +-
  0.39). Other bands: 90-110 0.79 +- 0.52; >= 100 0.74 +- 0.46.
- **The peak term matters at large M(ee).** At M(ee) >= 90 the fixed peak term (direct decays just above the edge) has a step of
  -20.1 rows, as large as the bump prediction. With the direct-decay scale 0.50 instead of a_P = 0.418 it is -24.1 (k 1.06). The note
  does not mention it; it changes no conclusion.
- **Expected k of the CB rule** (acceptance factorised per lepton; observed rc = N++/N-- step ratio 0.127, background rB = 0.229 in
  the cells): CB response 0.744 x (e-e- step). Partners in N: 2 sqrt(rc) = 0.713 for Dalitz/conversion legs only, so k = 0.96; with a
  third of Compton/delta legs (no e+e+), k = 0.98; with half, 1.02. Effect on f: at most +0.13 % (the CB step is 40.9 rows of 1303.8).
- Own f: 1303.8 / 3109.2 = 0.4193 (nominal 0.4200). Using k = 0.96 +- 0.48 instead of 1: f +0.1 +- 1.5 % (eestep -0.1 +- 1.2 %).

### V6. The excess over SMASH (f)
In short: real at the 2.3-3.7 sigma level, but resting on 11-14 simulated rows; material and D32 do not explain it.

- **Origins in SMASH** (standard set, pcls 2, S0, ownA): e-e- 133 rows (weighted 3748); leading origins C+D(eta) 915, C+C 839,
  C+Compton 630, C+D(other pi0) 459, D(eta)+D(other) 202, C+delta 181, Compton+D(eta) 163, Compton+Compton 67. **34 % of the weight
  (36 % near the edge) has a Compton or delta electron.** e+e+ 11 rows (weighted 327), none with a Compton or delta electron.
  e+e+/e-e- = 0.087 (all), 0.106 (Dalitz/conversion legs only).
- **SMASH x f against the data** (e-e-, rows whose true u <= +2 MeV, i.e. leptons that can come from the eta's pi0):

  | M(ee) | SMASH x f: step (data estimator) | SMASH x f: density below the edge | data step | ratio | z |
  |---|---|---|---|---|---|
  | all | 20.9 +- 13.1 | 21.3 (14 rows, +- 27 %) | 73.5 +- 15.6 | 3.5 / 3.5 | 2.6 / 3.1 |
  | >= 90 | 25.2 +- 10.3 | 18.2 (11 rows, +- 30 %) | 57.0 +- 9.0 | 2.3 / 3.1 | 2.3 / 3.7 |

- **Conversion material.** Near the edge the SMASH cross pairs have 0 / 1 / 2 photon interactions in 7 / 58 / 35 % of the weight. A
  ratio of 2.3 (3.5) needs 1.83 (2.40) times more photon interactions than HGeant. The bump's own conversion scale, kappa = 0.88 +-
  0.21 +- 0.35, disfavours that at 2.3 (3.7) sigma; the converted-photon scale of the peak (0.24 against 0.50 for the direct
  decays) even points the other way.
- **D32 (close pairs).** Rows are all combinations, so a lost close partner does not add like-sign rows. D32 can act only through the
  normalisation: if close pairs are lost in the data, f understates the eta yield for the wide pairs; the direct-decay scale (0.50
  against 0.42) bounds this at x 1.2. A merged pair (one track) would add sharp-step pairs, but V3 shows no such leptons.
- **Conclusion.** The excess is not explained. It may be the small simulated sample, or more Compton/delta electrons than in HGeant.
  It does not affect any number: the CB comes from the data and the data show the partners.

### V7. Y_peak and X17 (g)
In short: not touched; confirmed.

- E(m, q) reaches the peak window (0.53364) for M(ee) >= 101.5 / 104.0 / 110.5 / 119.5 MeV/c2 at M(pipi) 0.28 / 0.30 / 0.35 / 0.40;
  at most 0.5439 (q = m_pi0), so the upper side band SBhi gets none.
- Cross-pair share just below the edge (M(ee) >= 90, -15..0 MeV/c2): 57.0 / 101.5 = 0.56. e-e- rows with u <= 0 (S0, 0.26-0.46):
  SBlo (0.5250-0.5336) 320 (about 180 cross pairs), P 89 (about 50; 37 CB rows), SBhi 0. **The note's "SBlo 454 / 255" is the
  range 0.5225-0.5336 (R2), not SBlo.** This enters only the artefact scenario.
- With the expected k (0.96-1.02) the residual in P is below 1.5 rows (0.03 % of Y_peak 5088); with k = 0.96 +- 0.48 it is within
  +-18 rows (+-0.4 %), before the side bands, which carry the same pairs and reduce it further.
- X17: E(m, 24.5 MeV/c2) <= 0.5211 at M(pipi) = 2 m_pi (0.5223 at 0.26, the lowest row): 11.4-12.5 MeV/c2 below P. X17P holds 207
  e-e- and 13 e+e+ rows (S0); none has u <= +12 MeV/c2. The e-e- step at 10-30 MeV/c2 is 4.2 +- 5.8. The X17 fit is not touched.

### V8. Agreement table

| quantity | eestep | independent | verdict |
|---|---|---|---|
| rows (edgeu range) | 551 291 joined | 1 910 080 identical to edgeu (k, m4, mee, bits) | agree |
| e-e- / e+e+ step, all M(ee) | 73.5 +- 15.7 / 2.1 +- 6.5 | 73.5 +- 15.6 / 2.1 +- 6.5 | agree |
| f (own estimator) | 0.4200 | 0.4193 | agree |
| edge of cross pairs | toy: u <= -0.0055, -0.016 MeV | general proof; toy: u <= -0.0003 ... -0.11 MeV | agree |
| step position | best -2.5 MeV | data -2.5; toy predicts -2.5 for cross pairs, 0 for the bump | agree (new support) |
| step M(ee) < 50 / >= 90 | 3.6 +- 8.9 / 57.0 +- 9.1 | 3.6 +- 8.9 / 57.0 +- 9.0 | agree |
| e+e+ at >= 90; ratio | 7.3 +- 4.3; 0.13 +- 0.08 | 7.3 +- 4.3; 0.127 +- 0.079 | agree |
| opening angle | 20-60 deg | 20-35: 25.1, 35-60: 40.8; toy 21-61 deg | agree |
| no e+ in the event | 85 % | 85.0 % (toy 90 %) | agree |
| twin / shared ring / >= 3 e- | 6.8 / -2.7 / 0.6 | 6.8 +- 6.0 / -2.7 +- 4.8 / 0.6 +- 2.8 | agree |
| flip r, effect on f | < 1.24 %; < +1.6 % | < 1.23 %; < +1.6 % (e+ -> e-), < +1.9 % (both charges) | agree, wording |
| N, CB, N - CB at >= 90 | 52.6, 35.8, 16.8 +- 13.9 vs 15.0 | 52.7, 39.6, 13.1 +- 14.6 vs 14.7 +- 6.7 (MC) | agree |
| artefact reading | 2.7 sigma | 2.6 (stat, same form); 3.7 (N alone, stat); 2.4-2.9 with MC error | agree, wording |
| k_cross | 1.05 +- 0.39 | 0.96 +- 0.48 (+- 0.52 with MC); expected 0.96-1.02 | agree, error larger |
| why e-e- | soft-e+ acceptance (SMASH 0.08) | acceptance AND Compton/delta electrons (34-36 % of SMASH e-e-) | **correction** |
| pairs per decay | 1 e-e-, 1 e+e+, 2 e+e- | true for Dalitz/conversion legs; Compton/delta: 1 e-e-, 1 e+e-, no e+e+ | **correction** |
| data / SMASH x f | >= 3.5, about 3 sigma, 10 rows | 2.3-3.5, 2.3-3.7 sigma, 11-14 rows | **correction** |
| material, D32 | "could explain it" | material x 1.8-2.4 needed, disfavoured 2.3-3.7 sigma; D32 <= x 1.2; no merged pairs | **correction** |
| SBlo content | 454 rows, about 255 cross pairs | SBlo 320 / about 180; 454 is 0.5225-0.5336 | small correction |
| X17 edge | <= 0.5217, 12.0 MeV below P | <= 0.5211-0.5223, 11.4-12.5 MeV below P | agree |
| Y_peak | unchanged | unchanged (residual < 0.03 % with expected k) | agree |

### V9. Files and commands (from /home/przygoda/HADES/FAT/FAT)

    B='gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_veestep", kTRUE)'
    root -l -b -q -e "$B" 'plots/all_report/veestep_dump.C+O'   # 30 s: rows + PAT joins -> veestep/d_*, veestep_dump.log
    root -l -b -q -e "$B" 'plots/all_report/veestep_pat2.C+O'   # 2 s: ring pads, dE/dx of the like-sign rows -> veestep/d2_*
    python3 -I plots/all_report/veestep_ana.py                    # 5 s: Sect. V2-V5, V7 -> veestep_ana.log, veestep_ana.json
    python3 -I plots/all_report/veestep_toy.py                    # 3 s: Sect. V1 -> veestep_toy.log
    python3 -I plots/all_report/veestep_sim.py                    # 10 s: Sect. V3 (rings), V6, V7 -> veestep_sim.log

- Code: plots/all_report/veestep_dump.C, veestep_pat2.C, veestep_core.py, veestep_ana.py, veestep_toy.py, veestep_sim.py.
- Outputs: /hdd3/przygoda/HADES/FAT_out/all_report/veestep/ (column files d_*, d2_*, veestep_layout.txt, logs, veestep_ana.json).
- The copy of this note in docs/latex/notes_all/ was not updated (for the integrator).
