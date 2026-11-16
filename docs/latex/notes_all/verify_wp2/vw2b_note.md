# vw2b: adversarial method audit of WP2 f_win (material share of the pi+pi-e+e- channel)

Auditor prefix vw2b, 2026-11-14. Read-only on WP2 and the repository. Code and outputs:
`/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/verify_wp2/vw2b/`
(vw2b_dump*.C, vw2b_[a-n].py, out_*.txt).

## Verdict

- **As the material share of N - CB rows in B+P: CONFIRMED WITH CHANGES.**
  - Value: **14.7 +- 1.1 (stat) +1.0/-2.0 (syst) %**, band 12.5-16.5 %. WP2 gives 16.8 +5.0/-5.9 %.
  - All M4: 11.3 +- 0.4 (stat) +0.8/-1.5 %.
- **As the correction of sigma(pp -> eta X): REFUTED.**
  - The eta signal itself has almost no material component.
  - Its vertex-z distribution in the full-target data matches the LH2 simulation.
  - f_eta = 1.5 +- 4.5 % with a template-independent method; 0.0-2.7 % (+-1.3..2.5) from shape fits; 95 % CL upper limit 3-5 %.
  - The 16.8 % hypothesis is excluded at 3.4 sigma (cell-centre method) and at more than 5 sigma (shape fits).
  - Recommended: **f_eta = 2 +3/-2 %**, so 1 - f = 0.98 +0.02/-0.03 instead of 0.832.
  - With this, sigma(pp -> eta X) is about 17 % higher than WP2's value.
- **Why WP2 is wrong here.** The material rows are mostly non-eta pi+pi-e+e- combinations. The background model of the eta fits already absorbs them. Multiplying a_eta by (1 - f_win(N-CB)) therefore subtracts that material a second time.

## Checks

**(1) Is the "empty" target empty? First principles.**

Target, from the feb22 gen4 simulation geometry (`pat3/docs/geometry` on vae25, target.geo):
- LH2 cone from z = -148.4 to -98.6 mm (49.8 mm long), r 12 -> 10 mm.
- Mylar windows, each 0.1 mm: entrance (TWS4, r < 7.5 mm) at z = -148.5 mm; cell exit at -98.6 mm; outer exit window (TKW1) at -96.7 mm.
- LH2 annulus at r 8-12 mm from -174.2 to -148.4 mm.
- Carbon housing, 0.5 mm wall, at r = 20 mm, starting at z = -176.7 mm.

Areal densities and expected shares:
- LH2: 0.353 g/cm2 = 2.11e23 p/cm2. Mylar (0.3 mm): 0.0417 g/cm2 = 2.51e22 nucleons/cm2. Nucleon ratio 11.9 % (exit windows only: 7.9 %).
- With A^alpha scaling, the per-nucleon factor of Mylar is 0.45 / 0.61 / 1.0 for alpha = 2/3 / 0.8 / 1. Expected share: 5.3 / 7.3 / 11.9 % (0.2 mm of Mylar: 3.6 / 4.9 / 7.9 %).
- The measured row share (11.3-11.9 % in all M4) is therefore plausible.

Gas in the "empty" cell:
- Saturated H2 vapour (20 K, 1 atm) would be 1.9 % of the LH2 density.
- The elastic data exclude this. The uniform in-cell part of the empty-target elastic yield (Task E) is 0.08-0.12 % of the LH2 yield. That is also the limit on any H2 gas, and it is consistent with warm gas at 1 atm (0.12 %).
- The channel's own on-axis cell-centre ratio gives a looser limit of 0.8-1.0 %.
- **WP2's one-sided gas term (-3.6 % absolute in B+P) is too large by a factor of 3-30.**

**(2) Is the flux normaliser sound?**
- The front rows sit off axis (89-93 % at r >= 8 mm, 66-69 % at r >= 13 mm), so they come from beam halo on structures. Their composition is the same in full and empty runs: like-sign/unlike-sign 0.455 vs 0.459.
- F/(b E) from -300 to -160 mm on day 065 is flat at 0.92-1.10. There is no LH2 leak (simulation: 1.7e-5) and no annulus excess.
- **Problem: the halo-to-core ratio is not stable in time.**
  - h = 100 x U/N_S0 (front rows per channel row) ranges from 3.15 to 5.69 over the beam time (rms 19 %).
  - 2-hour blocks within a day scatter by 8.3 % beyond Poisson.
  - The material rows have a halo part: the per-day f_win follows f = 7.2 + 0.98 h (all M4; chi2 0.01/1) and f = 4.6 + 2.43 h (B+P; chi2 0.37/1). The r >= 8 mm part scales with h; the r < 4 mm part does not.
  - The empty-target days are not typical of their periods:
    - day 034: h = 5.31, period I mean 4.15;
    - day 047: h = 5.69, period II mean 5.31;
    - day 065: h = 3.26, period III mean 3.30.
  - WP2's day-to-period weighting implies h = 4.79, while the true L-weighted value is 4.17.
  - Reweighted with the halo model: B+P **14.70 +- 1.08**, all M4 11.29, fit range 12.45. Inverse-variance mean: B+P 13.04 +- 0.76.
- **Within-day mismatch.**
  - Day 034: the empty runs (01:26-01:47) came before any full-target run of that day. The first block has h = 4.77 against a day mean of 5.31, so the core part may be 11 % too high.
  - Day 047: the empty runs (14:10-14:50) fall in a transition of h from 6.9 to 4.8, giving -3..-17 %.
  - Day 065 is stable (h = 3.1-3.6 over the day, and the empty runs span 8.5 h).
- **Normaliser band, now dominated by day 065:** b_fit/b = 0.876, b_max/b = 1.064, i.e. -12 / +6 % of f_win. Task E's U1 normaliser (b x 1.17) violates the positivity bound of day 065 by 2.4 sigma.

**(3) Trigger, live time and run conditions.**
- Normaliser and yield come from the same PT3 events, so live time cancels. Only the trigger word 8192 is kept in both samples.
- The day-to-day differences in f_win are explained by the halo level, not by changing material. The elastic per-day trend goes the other way (0.70 -> 0.88 %).

**(4) The elastic w6 does not contradict 12-19 %.**
- Tight elastic cuts see only free hydrogen. In 0.3 mm of Mylar that is 0.50 % of the LH2 protons (measured raw share 0.40-0.46 %, 96 % of it at the exit).
- Nucleons per free proton in Mylar = 24. So 0.5 % x 24 x (0.45..1) = 5-12 %.
- The design's 0.7-2.1 % expectation left out this factor of 24.

**(5) "No eta below -160 mm" is consistent with material rows.**
- The expected material eta at -200..-170 mm is about 60 events, against a measured -19 +- 76: not sensitive.
- The decisive slice is -160..-150 mm, where 71-83 % of the rows are material.
  - Tagged eta there: T = 93 +- 29 (my refit; T for all z = 3279 +- 107, report 3312 +- 105). This is 2.9 % of T.
  - The LH2 simulation tail alone gives 2.4-2.8 %.
  - So the material rows carry r = 0.0 +- 0.1 times the eta per row of LH2.
- Over all slices (peak S, pre and S0, plus tagged T; 29 days), the fit of (1 - f) x LH2-sim + f x empty-row shape gives the f_eta values quoted in the verdict.
- Template-independent version: the cell-centre share [-140,-110) mm is 0.571 +- 0.024 (T) and 0.566 +- 0.045 (S), against 0.577 / 0.575 in the simulation. A 16.8 % material eta share would give 0.497.

**(6) Double counting.**
- L(1 - w6) and f_win act on different quantities, and C_sec cancels in the same-day ratio. There is no double counting there.
- **There is double counting in applying f_win(N-CB) to a_eta** (see the verdict).

## Problems (size as % of f_win = 16.8)

| # | problem | size |
|---|---|---|
| P1 | f_win(N-CB, B+P) used as the eta-specific share | -88 % (16.8 -> 2) |
| P2 | days 034/047 are not representative of their periods' halo | -12 % (-> 14.7); -22 % with inverse-variance weighting |
| P3 | gas term overstated | lower band -3.6 -> <= -1.0 (16-21 % of f_win) |
| P4 | within-day halo mismatch on days 034 and 047 | +-5 % after reweighting |
| P5 | normaliser band | -12 / +6 % (day 065) |
| P6 | design expectation 0.7-2.1 % was flawed | the measured row share is physical |
| P7 | small-count CB bias of the per-window CB (days 034/047) | about +1 % |
