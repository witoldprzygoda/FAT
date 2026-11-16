# Glossary of internal labels: p+p → π+π−e+e− reports, Parts I–III (HADES, feb22, T = 4.5 GeV)

Written 15–16 November 2026 for the readability pass (analyst's instructions I1–I4). **Version 2 (16 November 2026):** every
definition checked against the code (defs.h, the Part II and Part I selection code, the FAT configuration and cuts, the truth and
weight code) and against the six reports; the corrections are listed in §0b, the disagreements are re-examined in §16. It is a lookup tool for the editors: every
internal label that appears in the six reports, with (a) its exact definition, cuts and values, taken from the code and the notation
sections, with file:line; (b) a plain-words phrase that can be pasted into a caption; (c) for classes, a short legend name and a
longer name; (d) for windows and regions, the mass ranges. Section 15 has ready-made "Selection:" sentences, section 16 lists every
place where the reports (or a report and the code) disagree on a definition, and section 17 counts where each label is printed.

Nothing in this file changes a number. All numbers are copied from the code or the reports. "Count" always means one stored
π+π−e+e− combination (one ntuple entry); with all combinations kept, one event can give more than one count. The simulation never
gives the combinatorial background (CB).

Coverage: 351 numbered entries in sections 1–13 (614 individual labels and symbols when grouped rows such as L0–L11 or D1–D40 are counted per member) plus 41 entries added in version 2 (§13b, rows 352–392: the rows of the three notation tables that version 1 did not list), a disambiguation table of 34 multi-meaning tokens (§14), 36 ready-made "Selection:" sentences (§15), 23 disagreements with a status each (§16) and a per-report inventory of 288 label patterns (§17).

## 0. Key to the sources

| key | file |
|---|---|
| defs | /home/przygoda/HADES/FAT/FAT/plots/all_report/defs.h (frozen Part III definitions) |
| SimE | /home/przygoda/HADES/FAT/FAT/plots/all_report/SimEngine.h (copy of the simulation-branch engine; "pre" and the named presets for data and simulation) |
| Scan | /home/przygoda/HADES/FAT/FAT/plots/ScanEngine.h (Part II engine: windows, P_main, edge) |
| spec_c, res_c, opt_c, cb_c | /home/przygoda/HADES/FAT/FAT/plots/exp_report/{spec,res,opt,cb}_common.h (Part II selections) |
| e1f | /home/przygoda/HADES/FAT/FAT/plots/all_report/e1_fill.C (Part III cut sets E1–E4, X1–X4) |
| s5prep | /home/przygoda/HADES/FAT/FAT_simwt/plots/eta_sim/s5_blind_prep.C (Part I windows, groups, blind test) |
| frz_eta, frz_x17, frz_x17v1 | /home/przygoda/HADES/FAT/FAT_simwt/plots/eta_sim/tables/s5_opt_eta_frozen.txt, s5_opt_x17_frozen.txt, s5_opt_x17_frozen_v1.txt |
| truth | /home/przygoda/HADES/FAT/FAT_simwt/src/sim_truth.h (copy in plots/all_report/sim_truth.h) |
| phys | /home/przygoda/HADES/FAT/FAT_simwt/plots/eta_sim/sim_eta_physics.h (physics weights V0–V4, codes) |
| I, I5, Is | /home/przygoda/HADES/FAT/FAT/docs/latex/report_sim.tex, sec5_sim.tex, report_sim_short.tex (identical copies in FAT_simwt/docs/latex) |
| II, IIs | /home/przygoda/HADES/FAT/FAT/docs/latex/report_exp.tex, report_exp_short.tex |
| III, IIIs | /home/przygoda/HADES/FAT/FAT/docs/latex/report_all.tex, report_all_short.tex |
| tS1, tU, tE | docs/latex/tables_all/p3w_tableS1_v3.tex (families S1–S13), edgeu_syst_v3_cnt.tex (U1–U11), III:3940–3957 (E1–E10) |
| SPEC | /home/przygoda/HADES/FAT/FAT/docs/latex/notes_all/SPEC_report_all.md (Part III specification) |
| FATcfg | /home/przygoda/HADES/FAT/FAT/config.json (and config_epep.json, config_emem.json), src/setup_cuts.h, main.cc (FAT production cuts and flags; added in v2) |
| b4I | /home/przygoda/HADES/FAT/FAT_simwt/plots/eta_sim/s5_blind_ana_eta.C (Part I extraction methods M1–D5; added in v2) |
| spec5 | /home/przygoda/HADES/FAT/FAT_simwt/plots/eta_sim/s5_spec.C (Part I spectra and groups; added in v2) |
| wp3, wp6, flex, edgeu, e1a | /home/przygoda/HADES/FAT/FAT/plots/all_report/wp3_stageA.C, wp6_ana.C, flex_fit.C / flex_info.C / flex_lib.h, edgeu_core.py / edgeu_ana_v3.py, e1_ana.py (added in v2) |

Line numbers are those of the files on 16 November 2026 (after the last builds of 15 November, 23:46).

How to use the phrases: a caption of a short report takes the plain phrase only; a caption of a long report takes the phrase and
may add the label in parentheses, e.g. "Selection: ... both pions used in the event-vertex fit (S0)".

## 0b. Version 2: what was checked and what changed

**Checked.** The code behind every definition was opened and compared (defs.h, SimEngine.h, ScanEngine.h, spec/res/opt/cb_common.h,
e1_fill.C, s5_blind_prep.C, the frozen tables, sim_truth.h, sim_eta_physics.h), and in addition the FAT production itself
(config.json, src/setup_cuts.h, main.cc), the Part I extraction code (s5_blind_ana_eta.C, s5_spec.C) and the Part III package code
named in the reports (wp3_stageA.C, wp6_ana.C, flex_fit.C, flex_lib.h, edgeu_core.py, e1_ana.py). All 696 file:line citations of
version 1 (410 to the reports and the specification, 286 to the code) were printed next to the cited lines and read: they point to
the right content (one, spec_c:93, was off by one line and is fixed; SPEC:245 is now SPEC:246). The values of the pion-vertex flags were read from the data compact file (one ROOT pass over
four integer branches). Seven labels of §17 were recounted in all six PDFs (S0, PEAK, CR-U, F-main, B4, D40, isBest: 42 counts, all
equal) and the text extract of report_all_short.pdf was confirmed to be current. No number of the reports was changed.

**Corrected or completed in this version** (row numbers of §§1–13):
- Row 7 and §16 item 15: the "!= 0" and "== 1" pion-vertex tests are equivalent; the flags are 0 or 1 in every data entry.
- Rows 1, 2, 6, 9, 37, 46, 100: the FAT source of the trigger word, the vertex cut, isBest ≥ 0, the beam energy 4530 MeV, the
  old-code flag sel_pass_cor, the strict OA flags of X2old and the ECAL photon quality cuts.
- Row 11, 99, 100: the other CB grids (ECAL tag in 20 MeV/c² bins; u-edge in 10 MeV/c² cells of u per M(π+π−) slice).
- Row 56: the limit 0.413 GeV/c² is on the true scale (0.4116 on the data scale).
- Rows 57, 99, Q21, Q26, Q33: u exists only for M(e+e−) below m_π0 = 135 MeV/c² (the code drops the rest).
- Row 69: the Part I regions "below" and "above 1.00" are 0.30–0.4576 and 1.00–1.20 GeV/c² in the code.
- Rows 70, 71: widths of the two side-band sets; the 0.04 MeV/c² slivers between P and the rounded Part III side bands.
- Row 72: 0.5312–0.5566 is not exactly "× 0.99615, rounded" (that gives 0.5312–0.5565); it is slightly wider than both exact windows.
- Rows 87, 88: CR-U has three uses in the code (0.565–0.700 definition and flexible fits; 0.575–0.700 in E2 for stage A and B4;
  0.575–0.70 over all M(e+e−) for the typical-cuts display). Version 1 wrongly put the §e1 display under "M(e+e−) ≥ 10 MeV/c²".
- Row 98: origin of the two edge scales 0.99277 and 0.99274; Part II does state both (version 1 said no report does).
- Rows 115–117: T3 (same particle in both lepton slots) takes priority over the signal; T2 includes misidentified pions.
- Rows 124, 125: G_oth checked in the code: leptons of a forced η → π+π−e+e− with pions not from that η (X17 entries weigh 0);
  the plain phrase of G_D no longer says "non-η".
- Row 129: multi-pion confirmed in the code; G_vec is 0.02 % of it.
- Row 135: w and W mean different things in Part I and in Part III / the code (new §16 item 21).
- Row 144: the constant factors of D35 named per M(e+e−) class.
- Rows 151, 156: the Part I side-band ranges of M3 and B4 from the code.
- Row 167: the 48 variations use four fit ranges (both ends move), not only four lower bounds.
- Rows 181, 188, 196: free parameters of stage A; BR_C withdrawn; the definition and value of f_E2.
- Rows 230–242: each family S1–S13 now lists the variations of the specification (SPEC:280–293); S13 is not in the specification.
- Rows 275, 285: outcomes of D5 and D15 as written in III.
- §13b (new): 41 notation rows of the three long reports that version 1 did not list.
- §14: CR-U and side-band rows updated; two tokens added (w; E2a/E2b).
- §15: a rule for short reports (drop every parenthesised label), plain SMASH1 and Q36, completed Q21, Q26, Q28, Q33, Q35.
- §16: every disagreement re-examined against the code, with a status (settled by the code / documented in the reports / open);
  three new items (21–23).

---------------------------------------------------------------------------------------------------------------------------------

## 1. Building blocks (plain phrases used in every selection)

| # | item | exact definition (source) | plain phrase for captions |
|---|---|---|---|
| 1 | PT3 trigger, trigger bit | Data: the production keeps only trigbit == 8192 (SimE:404–407, 428; res_c:19). Simulation: (trigbit & 8192) != 0, simulated PT3 entries carry 12288 (SimE:407, 427; I:1789). Data yields are multiplied by 64/63 because 1/64 of the PT3 events carry the word 4096 (defs:52; III:6667). FAT: config.json trigger.selection = "PT3" → cut trigbit == 8192 (src/setup_cuts.h:32–34; main.cc:252–259). | "PT3 trigger (at least two sectors with an inner-TOF hit in coincidence with TOF or RPC, and at least three inner-TOF paddles hit; trigger bit 8192)"; for normalised data add "counts × 64/63, because only the first-arriving trigger bit is written and 1/64 of the PT3 events carry the PT2 bit 4096" (definition: P. Tlusty, LVL1 trigger in Feb22 p+p 4.5 GeV, 13 Dec 2024) |
| 2 | event vertex, vz > −500 mm | A reconstructed event vertex (eVertReco_z) exists; PAT writes vz = −1000 mm when the DST vertex finder (window −170 < z < −30 mm, r < 20 mm from the beam axis) finds fewer than two tracks (defs:115; SimE:430; II:3183–3185, 3459). Mandatory (applied in the FAT production: cut vertex_z > −500, src/setup_cuts.h:30; main.cc:250). | "a reconstructed event vertex" |
| 3 | target region, vz > −160 mm | Event vertex z > −160 mm (defs:148; spec_c:72; res_c:94). Target spans −152 to −98 mm at half maximum (I:1643–1644). | "event vertex inside the target region (z > −160 mm)" |
| 4 | PID | PAT, START-free velocity β_new; every track inside the polygon of its slot in (momentum, β_new): p and π± at 3σ, e± at 5σ; data cut file pp45_pid_cuts_exp_betanew_v1.root, simulation pp45_pid_cuts_sim_betanew_v6.root; lepton slots take only tracks with a correlated RICH ring, pion slots only tracks without one (a ring track outside the 3σ RICH–MDC matching ellipse goes to the hadrons) (II:1279–1281, 1318; I:871–872, 1372–1373). | "particles identified by the START-free velocity β_new (pions within 3σ, leptons within 5σ of the expected band; leptons with a RICH ring)" |
| 5 | cor_ok = 1 | All four energy-loss-corrected momenta are finite; the corrected π− momentum is missing for 0.7 % of the counts (π− at θ < 16°), such counts are dropped (defs:115; spec_c:78; I:1645–1647; II:3458). | "energy-loss-corrected momenta for all four tracks" |
| 6 | isBest ≥ 0 / isBest = 1 | PAT ranks the combinations of an event by the reduced χ² of the times of flight: best = 1, others 0; isBest ≥ 0 keeps all stored combinations (nominal everywhere), isBest = 1 only the best (defs:141; II:1282–1286; III:6658). FAT keeps isBest ≥ 0 (cut isBestAll, isBest > −0.5; src/setup_cuts.h:25–29; main.cc:249); every PAT entry has isBest 0 or 1 (II:814). | "all stored π+π−e+e− combinations of an event (an event can give more than one count)" / "only the best-ranked combination of each event" |
| 7 | both pions in the vertex fit ("vertex", "pions in vertex fit") | pip_isUsedVertex == 1 && pim_isUsedVertex == 1 (defs:143; SimE:437; res_c:95; cb_c:52; s5prep:157). Part II spec code writes != 0 (spec_c:13, 79). Checked in v2 on scan_v2/scan_compact.root: pip/pim/ep/em_isUsedVertex take only the values 0 and 1 in all 4 783 467 (us), 122 875 (pp) and 371 524 (mm) entries, and "both != 0" selects the same 3 640 366 / 51 544 / 156 154 entries as "both == 1": the two codings are equivalent (§16 item 15, settled). | "both pion tracks used in the fit of the event vertex" |
| 8 | Δz_π < 40 mm (dz_pi) | Δz_π = max(abs(z_π+ − vz), abs(z_π− − vz)), z of a pion track at its closest approach to the beam axis (defs:121, 138; s5prep:148; II:3502, Eq. at II:3794). | "both pion tracks within 40 mm of the event vertex along the beam (Δz_π < 40 mm)" |
| 9 | MM4 < 2.85 GeV/c² | Missing mass of the four tracks, mass of (beam + target − π+π−e+e−), T = 4530 MeV in the data (Scan:64; FAT config.json beam.kinetic_energy 4530.0); simulated MM4 at 4500 MeV in Part III (defs:118; III:1209–1210), at 4530 MeV in Part I (I:2055–2057). | "missing mass of π+π−e+e− below 2.85 GeV/c²" |
| 10 | OA, α_lab, α_η | OA: e+e− opening angle (lab). α_lab (oa_lab): angle between the π+π− and the e+e− momenta in the lab. α_η (oa_eta, "η-rest angle"): the same angle in the rest frame of a hypothetical η (measured π+π−e+e− momentum, mass set to m_η) (II:3499–3501; I:1823–1825). | "e+e− opening angle"; "angle between the π+π− and e+e− systems in the laboratory"; "the same angle in the η rest frame" |
| 11 | CB | Part II: 2√(N++ N−−) of the like-sign hypotheses π+π−e+e+ and π+π−e−e− with the same cuts, in 10 MeV/c² bins of M(π+π−e+e−), linearly interpolated to 2.5 MeV/c², k = 1 (Scan:13, 77–82, 357–400; II:3514). Part III: the same in 10 MeV/c² cells from 0.30 GeV/c², and spread over the cells of a layout (M(π+π−) slice × M(e+e−) class) by the like-sign share f(cell given 50 MeV band); error from 200 Poisson event-group replicas (defs:283–297; III:1307–1320). The two grids have the same 10 MeV/c² cell edges (Part II starts at 0.35, Part III at 0.30 GeV/c²). Other CB grids: the ECAL-tag spectrum in m5c uses 20 MeV/c² bins (res_c:383–384); the u-edge uses 10 MeV/c² cells of u in each M(π+π−) slice (row 99); the vertex study a smooth like-sign CB (row 178). Never from simulation. | "combinatorial background CB = 2√(N++ N−−) from the like-sign data (π+π−e+e+, π+π−e−e−) with the same cuts, in 10 MeV/c² cells, k = 1" (Part III layouts add: "shared out over the M(π+π−) and M(e+e−) cells by the like-sign share of each 50 MeV/c² band") |
| 12 | data units | Raw counts; visible yield or cross section (N − CB) × 64/63 / L, L = 6079.2 nb⁻¹ (analysed runs, ±6.9 %) (defs:42–58; III:1189–1190). Part I used L_data = 6080.5 nb⁻¹ (all runs, then unconfirmed) as a scale (I:1785). | "raw counts" / "(N − CB) × 64/63 / L, L = 6079.2 nb⁻¹" |
| 13 | SMASH as used in Part III | SMASH p+p 4.5 GeV, 10^10 events in three sets of the same events (standard; lepton: dilepton decays forced; η: η decays forced), each L_sim = 221.05 nb⁻¹; ownership mode A; weight w (PDG BR) × physics weights V3 (code 108); M4' = 0.99615 [m4_true + 1.0 (m4 − m4_true)], M(π+π−)' = 0.9952 M(π+π−), M(e+e−) not scaled; MM4, y4 at 4500 MeV; compared sum T1 + T2 (defs:51, 85–95, 365; III:1193–1217). | "SMASH simulation (p+p, 4.5 GeV), every channel at its SMASH cross section × PDG branching ratio × η-decay physics weights, masses moved to the data scale" |
| 14 | SMASH as used in Part I | The same sets, mode A, w (rules v2) × V3, no mass transfer, windows centred on the simulated peak, counts at L_data = 6080.5 nb⁻¹ (I:1005–1137; I5:28–42). | "SMASH simulation only (no combinatorial background), η-decay physics weights V3" |
| 15 | errors | Event groups: adjacent entries of one event (equal file_id, run, event); variance of a sum = sum over groups of the squared group content (II:3452; III:6657). MC: rms of 100 bootstrap replicas of the simulation (Poisson(1) weight per simulated event group) (III:1357–1359). | "errors from whole events (event groups)"; "simulation-statistics error (MC)" |

---------------------------------------------------------------------------------------------------------------------------------

## 2. Selections

All selections contain building blocks 1–6 unless the row says otherwise ("pre" = blocks 1–6).

| # | label (where used) | exact definition (source) | plain phrase |
|---|---|---|---|
| 16 | "pre" (I, II, III, all shorts) | isBest ≥ 0 && vz > −500 mm && PT3 && vz > −160 mm && cor_ok == 1, after the PAT PID (defs:115, 130, 149; SimE:404–408, 429–431; spec_c:11, 78; res_c:94, 106; opt_c:151–152; cb_c:51; frz_eta header). Notation: I:1790, 1636–1650; II:3462, 1451–1457; III:6662; Is:664; IIs:417, 802; IIIs:1378. Identical to BUMP. | "PT3 trigger; event vertex in the target region (z > −160 mm); β_new particle identification; corrected momenta for all four tracks; all stored combinations per event (isBest ≥ 0)" |
| 17 | S0 (III, IIIs; D1) | "pre" && pip_isUsedVertex == 1 && pim_isUsedVertex == 1 (defs:117, 129, 148; SimE kPrePiVtx 409, 437). Nominal selection of Part III for bump, peak and X17. Notation III:6663; tab:sel III:1229–1252 (row III:1236). | "… as "pre", plus both pions used in the event-vertex fit" |
| 18 | "pre + vertex" (II, IIs), "pre + pions in vertex fit" (I, SimE:412), "prevtx" (spec_c:96) | Same as S0 (res_c:107–108; opt_c:158; cb_c:67). | same as S0 |
| 19 | dz40 (III) = DIRECT (I, II) | "pre" && Δz_π < 40 mm (defs:118, 131, 150; res_c:122–123; opt_c:161; spec_c:108–109; cb_c:68; frz_eta). Part I: frozen DIRECT selection, optimised for class 5 in the simulated peak window with class 4 as background (I5:608–609, 745–747). Part III: variant only (D4; family S11). | "… "pre", plus both pion tracks within 40 mm of the event vertex along the beam" |
| 20 | PEAK (I, II, III) | "pre" && Δz_π < 40 mm && MM4 < 2.85 GeV/c² (defs:118, 132, 151; res_c:120–121; opt_c:162; spec_c:106–107; cb_c:69; frz_eta). Part I target: classes 4+5 in the simulated peak window. Simulated MM4: 4530 MeV in Part I, 4500 MeV in Part III (defs:118). Part II's recommended η-peak selection (II:2652–2666). | "… "pre", plus both pion tracks within 40 mm of the event vertex and missing mass of π+π−e+e− below 2.85 GeV/c²" |
| 21 | BUMP (I, II) | = "pre", no further cut (opt_c:163; spec_c:19, 94–95; frz_eta). Part I target: classes 2+3 in the simulated bump window. | "the base selection only ("pre"; no further cut)" |
| 22 | X17 (frozen X17 selection, version 2; I, II) | "pre" + both pions in the vertex fit + 15.5 ≤ M(e+e−) < 24.5 MeV/c² + 0.5333 ≤ M4 < 0.5587 GeV/c² (simulation peak window) (frz_x17:7–21; opt_c:164; II:1518–1528; Is:500). | "… "pre", plus both pions in the event-vertex fit, 15.5 ≤ M(e+e−) < 24.5 MeV/c² and M(π+π−e+e−) 0.5333–0.5587 GeV/c² (simulation scale)" |
| 23 | X17 selection v1 (I5 only) | "pre" + both pions in the vertex fit + 13.0 ≤ M(e+e−) < 26.5 MeV/c² + 0.5330 ≤ M4 < 0.5555 (frz_x17v1; s5prep:103; I5:803, 907). Used in the Part I blind test, replaced by version 2. | "… "pre", plus both pions in the event-vertex fit, 13.0–26.5 MeV/c² in M(e+e−), M(π+π−e+e−) 0.5330–0.5555 GeV/c²" |
| 24 | X17 band, X17 region (II, IIs, III) | "pre" + both pions in the vertex fit + 15.5 ≤ M(e+e−) < 24.5 MeV/c², no M4 window (res_c:98, 125–126; cb_c:58, 76; II:3507, 2619–2621; IIs:418). The M4 window is applied separately (§4 rows 72–75). | "… "pre", plus both pions in the event-vertex fit and 15.5 ≤ M(e+e−) < 24.5 MeV/c²" |
| 25 | vz150, vz170 (III) | S0 with vz > −150 / vz > −170 mm instead of −160 (defs:119, 133–134, 152–153). | "as S0 but with the event vertex at z > −150 mm (or z > −170 mm)" |
| 26 | none (III) | "pre" without vz > −160 mm: isBest ≥ 0, vz > −500, PT3, cor_ok (defs:119, 135, 154). Part II: "pre without vz > −160" (opt_c:153). | "as "pre" but without the target-region cut (any reconstructed event vertex)" |
| 27 | fat (II code, cb_c:64: "FAT only (vz > −500)") | isBest ≥ 0 after the FAT production cuts only (vz > −500, PT3); cor_ok applied by the engine (cb_c:7–8). | "trigger, event vertex and PID only" |
| 28 | S0_best, pre_best (III tab:sel), isBest = 1 variants (II) | S0 or "pre" with isBest == 1 (III:1243–1244; SimE:451; res_c Best1). Data only (simulated isBest is unrealistic, D2). | "as S0 (or "pre") but only the best-ranked combination of each event" |
| 29 | S0_corfb, pre_corfb, corfb (III) | cor_ok fall-back: counts with cor_ok == 0 kept with the reconstructed (not energy-loss-corrected) π− momentum (*_fb columns) (defs:120; Scan:20–24, 273 Data::SetFallback; III:1245–1246, 6665). Test only. | "as S0 (or "pre") and, where the corrected π− momentum is missing, the uncorrected one" |
| 30 | pre+OA>2, "+ OA > 2°" variants (I, II) | "pre" && OA > 2° (SimE:436; cb_c:66; spec_c:114; res_c OA2). Never a fixed cut (II:3499). | "… plus e+e− opening angle above 2°" |
| 31 | peak nominal (I, II) | "pre" && α_lab < 55° && M(π+π−) > 0.30 GeV/c² (SimE:440; res_c:96, 109; opt_c:154; spec_c:80; II:1494–1501). Chosen on these data in the first cut-scan study (2026-11-11). | "… "pre", plus π+π−/e+e− lab angle below 55° and M(π+π−) above 0.30 GeV/c²" |
| 32 | peak nominal + vertex (II) | peak nominal && both pions in the vertex fit (res_c:111–113; opt_c:155). | "… plus both pions in the event-vertex fit" |
| 33 | bump nominal (I, II) | "pre" && α_lab < 40° && 0.30 < M(π+π−) < 0.41 GeV/c² (SimE:441; res_c:97, 114; opt_c:156; spec_c:81). | "… "pre", plus π+π−/e+e− lab angle below 40° and 0.30 < M(π+π−) < 0.41 GeV/c²" |
| 34 | bump nominal + vertex (II) | bump nominal && both pions in the vertex fit (res_c:116–118; opt_c:157). | "… plus both pions in the event-vertex fit" |
| 35 | peak / bump nominal + Δz_π < 40 (II) | peak or bump nominal && Δz_π < 40 mm (opt_c:166–167). | "… plus both pion tracks within 40 mm of the vertex" |
| 36 | slides (I, II, III) | "pre" && α_lab < 50° && M(π+π−) < 0.46 GeV/c² && α_η > 140° && OA > 2° (SimE:438; opt_c:168; I:1683; II:1505). Same cuts as cut set E1. | "… "pre", plus lab angle below 50°, M(π+π−) below 0.46 GeV/c², η-rest angle above 140° and opening angle above 2° (the cuts of earlier slides)" |
| 37 | old code (I, II, III) | "pre" && α_lab < 50° && M(π+π−) < 0.42 GeV/c² && α_η > 140° (= FAT flag sel_pass_cor: oa_pippim_epem_lab < 50°, m_pippim < 0.420 GeV/c², oa_pippim_epem_eta_rest > 140°, on the corrected values; src/setup_cuts.h:59–64, main.cc:437–443) (SimE:439; opt_c:169; II:1493, 1506). Same cuts as cut set E3. | "… "pre", plus lab angle below 50°, M(π+π−) below 0.42 GeV/c² and η-rest angle above 140° (cuts of the earlier code)" |
| 38 | E1 (cut set; III §e1, sans-serif) | on "pre" (plain name) or on S0 ("S0 + E1"): α_lab < 50°, α_η > 140°, M(π+π−) < 0.46, OA > 2° (e1f:26, 83, 96–98; III:1687, 6707). Simulation cut on M(π+π−)' = 0.9952 M(π+π−) (III:1627–1628). | as row 36 ("the cuts of earlier slides"), and say "on top of "pre"" or "on top of S0" |
| 39 | E2 (cut set) | E1 without OA > 2°: α_lab < 50°, α_η > 140°, M(π+π−) < 0.46 (e1f:27, 84, 99; III:1688). | "… lab angle below 50°, η-rest angle above 140°, M(π+π−) below 0.46 GeV/c²" |
| 40 | E3 (cut set) | α_lab < 50°, α_η > 140°, M(π+π−) < 0.42 (old code) (e1f:27, 85, 100; III:1689). | as row 37 |
| 41 | E3oa | E3 + OA > 2° (e1f:27, 86, 101; III:1690). | "… plus opening angle above 2°" |
| 42 | E4 (cut set, "η band") | α_lab < 50°, α_η > 140° only (e1f:28, 87, 102; III:1691). | "… lab angle below 50° and η-rest angle above 140°" |
| 43 | E4oa | E4 + OA > 2° (e1f:28, 103; III:1692). | "… plus opening angle above 2°" |
| 44 | X1 (cut set) | E2 + OA < 6° (e1f:29, 104; III:1693, 6708). | "… E2 cuts plus opening angle below 6°" (write the cuts out) |
| 45 | X2 (cut set) | E2 + 2° < OA < 6° (e1f:29, 105; III:1694). | "… plus opening angle between 2° and 6°" |
| 46 | X2old (cut set) | E3 + 2° < OA < 6° (old code with its OA flags: FAT pass_oa_min = OA > 2° and pass_oa_max = OA < 6°, both strict, not the inclusive FAT window oa_pass 2–6°; src/setup_cuts.h:38–48, main.cc:446–447; Scan:349–350 SelOldCodeOA) (e1f:30, 106; III:1695). | "… old-code cuts plus opening angle between 2° and 6°" |
| 47 | X3 (cut set) | S0 + OA < 6° (e1f:30 "S0_O6"; III:1696). | "S0 plus opening angle below 6°" |
| 48 | X4 (cut set) | "pre" + OA < 6° (e1f:30 "pre_O6"; III:1696). | ""pre" plus opening angle below 6°" |
| 49 | V, Nv (vertex study; II, III, shorts) | PT3, isBest ≥ 0, cor_ok, read from the 29 PAT day files without a vertex requirement; V: vz > −500 mm (event vertex; the FAT entries), Nv: vz = −1000 mm (whole events are split; each class gets its own CB from its own like-sign counts, k = 1); no vz > −160 cut (II:3192–3196; III:6673; IIIs:1382). | "counts with (V) and without (Nv) a reconstructed event vertex, PT3, all combinations, corrected momenta; no target-region cut" |

---------------------------------------------------------------------------------------------------------------------------------

## 3. Ranges of M(π+π−) and classes of M(e+e−)

| # | label | definition (source) | plain phrase / range |
|---|---|---|---|
| 50 | S1 (slice) | M(π+π−) < 0.34 GeV/c² on the data scale (simulation: M(π+π−)' = 0.9952 M(π+π−)) (defs:162–170; III:1288, 6700; IIIs:1396). Not family S1 (§9) and not sector 0–5 (§14). | "M(π+π−) below 0.34 GeV/c²" |
| 51 | S2 (slice) | 0.34 ≤ M(π+π−) < 0.37 | "M(π+π−) 0.34–0.37 GeV/c²" |
| 52 | S3 (slice) | 0.37 ≤ M(π+π−) < 0.39 | "M(π+π−) 0.37–0.39 GeV/c²" |
| 53 | S4 (slice) | 0.39 ≤ M(π+π−) < 0.43 | "M(π+π−) 0.39–0.43 GeV/c²" |
| 54 | S5 (slice) | M(π+π−) ≥ 0.43; above the kinematic limit of the bump (edge limit 0.4116, defs:99, 110) | "M(π+π−) above 0.43 GeV/c² (no bump possible there)" |
| 55 | S1E1 … S5E2 (F10 cells, III) | slice j × class E1/E2, e.g. S2E2 = 0.34–0.37 × M(e+e−) ≥ 10 MeV/c² (defs:397 F10Name; III:1393, 2093, 2100). | "M(π+π−) 0.34–0.37 GeV/c² and M(e+e−) ≥ 10 MeV/c²" (spell out) |
| 56 | four M(π+π−) slices of the per-slice edge (II, III §bedge) | 0.30–0.34, 0.34–0.37, 0.37–0.39, 0.39–0.40 GeV/c²; slice 0.40–0.41 fitted but not used; control slice 0.42–0.46 (above the kinematic limit m_η − m_π0 = 0.413 GeV/c² on the true scale, 0.4116 on the data scale; II:1610–1611; defs:110); in each slice the step is averaged over 21 values of M(π+π−) (Scan:1035–1046; II:1604) (II:1596–1611; III:3893). Note: not the S1–S4 slices of F10. | "in the four M(π+π−) slices 0.30–0.34, 0.34–0.37, 0.37–0.39 and 0.39–0.40 GeV/c²" |
| 57 | u-edge range | 0.30 ≤ M(π+π−) < 0.40 GeV/c² (III:4321–4322; II:3530 for h_u; edgeu_ana_v3.py:165). u exists only for M(e+e−) < m_π0 = 0.135 GeV/c²: the code sets u to "not a number" above (edgeu_core.py:146, 166); Part II says so (II:2865), Part III does not (§16 item 23). "All M(e+e−)" in a u caption therefore means M(e+e−) below 135 MeV/c². | "0.30 ≤ M(π+π−) < 0.40 GeV/c², M(e+e−) below 135 MeV/c²" |
| 58 | E1 (class) | M(e+e−) < 10 MeV/c² (defs:171–173 kEsplit = 0.010; III:1289, 6701). Not cut set E1, not family E1, not spec. Fig. E1. | "M(e+e−) below 10 MeV/c²" |
| 59 | E2 (class) | M(e+e−) ≥ 10 MeV/c². The class of every stage-A, F-main and B4 fit after the E1 fall-back (III:2073–2082). | "M(e+e−) of 10 MeV/c² or more" |
| 60 | M(e+e−) ≥ 20 MeV/c² (F5g, "golden") | defs:174 kF5gMee; III:1291. | "M(e+e−) of 20 MeV/c² or more" |
| 61 | G3 classes (III, variants only) | M(e+e−) edges 0, 8, 10, 12, 20 MeV/c², 2 GeV/c² (5 classes); M(π+π−) bins: 0, then 5 MeV/c² from 0.330 to 0.445, then 0.485, 0.507, 2.0 (27 bins) (defs:189–201). | spell out the class used |
| 62 | X17 M(e+e−) band | 15.5 ≤ M(e+e−) < 24.5 MeV/c² (res_c:98; opt_c:78; defs:244). | "15.5 ≤ M(e+e−) < 24.5 MeV/c²" |
| 63 | X17 blinded region and fit side bands (III) | masked 12 ≤ M(e+e−) < 28 MeV/c² until the freeze; side bands of the r = 0 fit 2–12 and 28–62 MeV/c² (specification: 2–13 and 27–62) (III:5881–5897; SPEC:246). | "M(e+e−) 12–28 MeV/c² (kept hidden until the procedure was frozen)" |
| 64 | K0S veto (III, variant A5) | 0.485 ≤ M(π+π−) < 0.507 GeV/c² (defs:191; III:2214). | "M(π+π−) 0.485–0.507 GeV/c² removed (K0S)" |
| 65 | Part II M(e+e−) and OA classes | M(e+e−) 0–5, 5–10, 10–20, 20–50, 50–100, 100–270 MeV/c²; OA 0–1, 1–2, 2–4, 4–6, 6–10, 10–20, > 20° (tables_exp/res_mee.tex). | spell out |

---------------------------------------------------------------------------------------------------------------------------------

## 4. Windows, regions and layouts in M(π+π−e+e−) (GeV/c², data scale unless stated; simulation enters through M4')

| # | label | range and definition (source) | plain phrase |
|---|---|---|---|
| 66 | B (bump window; III), "data bump window" (II) | 0.4550 ≤ M4 < 0.5225 (defs:226, 233, 240; Scan:112–113; II:1544; III:1261). Upper end = lower bound of the P_main fit range. Also region R1. | "bump window, M(π+π−e+e−) 0.4550–0.5225 GeV/c²" |
| 67 | P (peak window; III), "data peak window" (II), R_pk (II) | 0.53364 ≤ M4 < 0.55416 = 0.5439 ± 1.9 × 5.4 MeV/c² (0.9426 of the frozen Gaussian) (defs:226, 234, 241; Scan:107–109; II:1543). Also region R3. Not P (purity, Part I), not pion class P0–P9. | "peak window, M(π+π−e+e−) 0.53364–0.55416 GeV/c²" |
| 68 | simulation peak / bump window (I; II "sim. windows") | peak 0.5333 ≤ M4 < 0.5587 (0.5460 ± 1.9 × 6.7 MeV/c², centred on the simulated peak); bump 0.4576 ≤ M4 < 0.5255 (data bump window × 1.0057) (s5prep:102; spec_c:68–69; I5:46–60; II:1547). | "M(π+π−e+e−) 0.5333–0.5587 GeV/c² (simulation scale)" / "0.4576–0.5255 GeV/c² (simulation scale)" |
| 69 | other regions of Part I | below the bump window, 0.30–0.4576; gap 0.5255–0.5333; 0.5587–0.70; 0.70–1.00; above 1.00, i.e. 1.00–1.20 (I5:61 prints only "below 0.4576" and "above 1.00"; the bounds are in the code, s5prep:363–364, 659, 973–974). | spell out with the bounds |
| 70 | SB, peak side bands (III, IIIs) | 0.5250 ≤ M4 < 0.5336 (SBlo) and 0.5542 ≤ M4 < 0.5628 (SBhi) (defs:226, 242–243; III:1263, 6694; IIIs:1395). Each 8.6 MeV/c² wide. The rounded inner edges leave two slivers of 0.04 MeV/c², 0.5336–0.53364 and 0.55416–0.5542, in neither P nor SB (defs:241–243; no effect at the printed precision). | "side bands 0.5250–0.5336 and 0.5542–0.5628 GeV/c²" |
| 71 | side bands, R_SB (II) | 0.5225 ≤ M4 < 0.53364 and 0.55416 ≤ M4 < 0.56530 (res_c:76–78; II:1545, 1627–1637, 3510; tables_exp/res_mee.tex). Each 11.14 MeV/c² wide, touching P; the lower one starts at the P_main fit bound 0.5225. Differs from row 70 (§16 item 1). | "side bands 0.5225–0.53364 and 0.55416–0.5653 GeV/c²" |
| 72 | X17c (III) | 0.5312 ≤ M4 < 0.5566 × 15.5 ≤ M(e+e−) < 24.5 MeV/c² (defs:227, 244; III:1264, 6695). 6424 counts (S0). defs.h and III:1264 call it "simulation 0.5333–0.5587 × 0.99615"; the exact products are 0.531247–0.556549, so the code window is slightly wider at both ends (the upper edge 0.5566 is the 4-digit rounding of Part II's 0.556551, not of 0.556549); this is why it holds 6424 counts against Part II's 6395 (III:5892–5893; §16 item 9). | "M(π+π−e+e−) 0.5312–0.5566 GeV/c² and M(e+e−) 15.5–24.5 MeV/c²" |
| 73 | X17 M4 window "moved to the data scale" (II, recommended) | 0.531249 ≤ M4 < 0.556551 (simulation window × 0.5439/0.5460), printed 0.5312–0.5566 (opt_c:174–176; II:2622–2623, 2632–2644). 6395 counts. | as row 72 |
| 74 | X17P (III) | P × 15.5 ≤ M(e+e−) < 24.5 MeV/c² (defs:231–232, 260; III:1265). The window of the count 5258 (D19). | "peak window and M(e+e−) 15.5–24.5 MeV/c²" |
| 75 | η peak band PB; "X17 window" of the typical-cuts section (III §e1) | PB = 0.5325 ≤ M4 < 0.5575 (fine bins 93–102) (e1f:146, 152–153; III:679, 1854, 6188); the §e1 "X17 window" = PB × 15.5 ≤ M(e+e−) < 24.5 MeV/c² (e1f:155; III:1647, 1671). | "η peak band, M(π+π−e+e−) 0.5325–0.5575 GeV/c²" |
| 76 | X2D (III) | M4 0.5225–0.6100 in 35 bins of 2.5 MeV/c² × M(e+e−) 2–62 MeV/c² in 30 bins of 2 MeV/c² (defs:186–188; III:5901). The X17 fit plane. | "M(π+π−e+e−) 0.5225–0.6100 GeV/c² × M(e+e−) 2–62 MeV/c²" |
| 77 | P_main fit range (II) | 0.5225–0.6100 (Scan:102–103; II:1546). | "fit range 0.5225–0.6100 GeV/c²" |
| 78 | F10 range, "fit range" (III) | 0.400–0.700 (defs:185, 258; III:1268). | "M(π+π−e+e−) 0.40–0.70 GeV/c²" |
| 79 | R0 | M4 < 0.455 (defs:230, 248) | "below 0.455 GeV/c²" |
| 80 | R1 (region) | = B, 0.4550–0.5225 (defs:249). Not rule R1 (§9). | "bump window" |
| 81 | R2 | 0.5225–0.53364 (defs:250) | "0.5225–0.5336 GeV/c²" |
| 82 | R3 | = P, 0.53364–0.55416 (defs:251) | "peak window" |
| 83 | R4 | 0.55416–0.565 (defs:252) | "0.5542–0.565 GeV/c²" |
| 84 | R5 | 0.565–0.700 (defs:253) | "0.565–0.70 GeV/c²" |
| 85 | R6 | 0.700–0.900 (defs:254) | "0.70–0.90 GeV/c²" |
| 86 | R7 | 0.900–1.100 (defs:255) | "0.90–1.10 GeV/c²" |
| 87 | U (upper region of B4, III) | 0.575–0.700, M(e+e−) ≥ 10 MeV/c², all five M(π+π−) slices (wp6_ana.C:17, 897; III:2765–2772, 6697). The specification had U = 0.565–0.700 (SPEC:205). Not U1–U11, not Part II's halo counts U. | "upper region 0.575–0.70 GeV/c²" |
| 88 | CR-U (III) | defined as 0.565 ≤ M4 < 0.700, all slices (defs:228, 245, 274; SPEC:74). Adopted after the fall-backs for stage A and B4: 0.575–0.700 with M(e+e−) ≥ 10 MeV/c² (wp3_stageA.C:901–902, 1401; wp6_ana.C:14–17, 941; III:2088–2092, 1269, 6703; IIIs:1397). The flexible fits keep 0.565–0.700 (flex_fit.C:144; flex_info.C:28; III:757, 5460; IIIs:981). The typical-cuts display (b) of §e1 normalises in 0.575–0.70 over ALL M(e+e−) and all M(π+π−) (e1f:140, window "CR"; e1_ana.py:218), although III:1639 calls it "the control region of stage A" (§16 items 2, 22). | stage A, B4: "control region above the peak, M(π+π−e+e−) 0.575–0.70 GeV/c² and M(e+e−) ≥ 10 MeV/c²"; flexible fits: "M(π+π−e+e−) 0.565–0.70 GeV/c²"; §e1 display: "0.575–0.70 GeV/c², all M(e+e−)" |
| 89 | CR-K (III) | slice S5 (M(π+π−) ≥ 0.43) × 0.400 ≤ M4 < 0.5225 (specification "0.522"); fall-back 1: M(π+π−) ≥ 0.44; then dropped (defs:228, 246, 275; III:1270, 2087–2101). | "kinematic control region M(π+π−) ≥ 0.43 GeV/c², M(π+π−e+e−) 0.40–0.5225 GeV/c² (not used: it failed its tests)" |
| 90 | CR-L (III) | 0.5225 ≤ M4 < 0.5300, test only (defs:229, 247, 276; III:1271). | "test region 0.5225–0.53 GeV/c²" |
| 91 | J (jump; III) | J = R(0.560–0.595) / R(0.4875–0.5225) of the "before" ratio data/SMASH (defs:231, 256–257; III:6698). | "ratio of data/SMASH in 0.560–0.595 to that in 0.4875–0.5225 GeV/c²" |
| 92 | BEFORE range and layout (III) | M4 0.30–1.10 in 5 MeV/c² (160 bins); the other variables per region R0–R7; map M4 0.30–1.10 × M(π+π−) 0.25–1.05 in 10 MeV/c² (defs:202–221, 259). | "M(π+π−e+e−) 0.30–1.10 GeV/c², 5 MeV/c² bins" |
| 93 | fine grid (III) | M4 0.30–1.10 in 2.5 MeV/c² (320 bins) (defs:202–205). | "2.5 MeV/c² bins" |
| 94 | 50 MeV/c² bands of the CB spread (III) | 16 bands 0.30–0.35, 0.35–0.40, … 1.05–1.10 (defs:293–296). | "each 50 MeV/c² band of M(π+π−e+e−)" |
| 95 | F10 (III layout) | ten histograms S1–S5 × E1/E2, M4 0.400–0.700 in 70 bins: 5 MeV/c² in 0.400–0.520 and 0.570–0.700, 2.5 MeV/c² in 0.520–0.570 (defs:175–185; III:1288–1290). | "the five M(π+π−) slices × the two M(e+e−) classes, M(π+π−e+e−) 0.40–0.70 GeV/c²" |
| 96 | F5g (III layout) | the five slices with M(e+e−) ≥ 20 MeV/c² only, F10 binning (defs:174; III:1291). | "the five M(π+π−) slices, M(e+e−) ≥ 20 MeV/c²" |
| 97 | G3 (III layout, variants) | M4 0.400–0.850 in 2.5 MeV/c² (180) × 27 M(π+π−) bins × 5 M(e+e−) classes (defs:189–201; III:1292–1294). | "fine M(π+π−e+e−) × M(π+π−) × M(e+e−) grid" |
| 98 | E(m, q) (III; I "later"), M_edge (II) | bump edge on the data scale: E(m,q) = 0.99274 √(m_η² − m_η (E_π0 − p_π0)(1 − q²/m_π0²)), E_π0 = (m_η² + m_π0² − m'²)/(2 m_η), m' = m/0.9970, m = M(π+π−), q = M(e+e−); no bump for m > 0.4116 (defs:96–110; III:1272–1286; II:3529). Part II per-slice edge (also Part III §bedge): m_edge = 0.99277 M4max(M(π+π−)/0.9970) with M(e+e−) → 0, printed 0.9927 in Part II (Scan:120, 141–162; wp6_fill.C:13; II:1586–1604; III:531, 2906, 2997, 3891). The two scales: 0.99277 = 0.5439/m_η, the rounded frozen peak position (scan::kScale); 0.99274 ≈ 0.54388/m_η, the free-fit position (Part II quotes 0.9927 ± 0.0009 and the difference of 4 × 10⁻⁵, II:2800–2801, and uses 0.99274 for its M_edge of u, II:2859–2864; res_peak.C:401; defs:100). Effect on the edge about 0.015 MeV/c² (notes_all/edgeu_results.md:119–120; §16 item 14). | "the kinematic upper edge of the bump (largest M(π+π−e+e−) of η → π+π−π0, π0 → e+e−γ with the photon lost) for the measured M(π+π−) and M(e+e−)" |
| 99 | u | u = M(π+π−e+e−) − E(M(π+π−), M(e+e−)), defined for M(e+e−) < m_π0 only (edgeu_core.py:146, 166); fit range −45 ≤ u < +45 MeV/c² in 36 bins of 2.5 MeV/c² with a bin edge at u = 0; CB = 2√(N++N−−) in 10 MeV/c² cells of u, separately in the M(π+π−) slices 0.30–0.34, 0.34–0.37, 0.37–0.39, 0.39–0.40 GeV/c² (edgeu_core.py:65, 264, 343; edgeu_ana_v3.py:165; III:4321–4322, 4345–4347, 6724; II:2864–2868). Not u of the X17 model (§14). | "distance to the kinematic edge of the bump" |
| 100 | ECAL-tag windows (II, III) | good ECAL photon (ecal_pid == 1, 0.8 < β < 1.2, cluster energy > 100 MeV: FAT cut set ecal_quality, src/setup_cuts.h:50–54), of up to three such photons the one whose M(e+e−γ) is closest to 0.140 GeV/c² (Scan:219–220), with 0.119 ≤ M(e+e−γ) < 0.173 GeV/c² and M(e+e−) < 0.10 GeV/c²; m5c = M(π+π−e+e−γ) with the photon energy rescaled to M(e+e−γ) = m_π0; T fitted in m5c 0.49–0.62 (2.5 MeV/c² bins, CB in 20 MeV/c² bins) with a frozen Gaussian 0.5463 / 6.2 MeV/c² + pol2 (res_c:375–398; II:1613–1618). ECAL sector 0 is switched off (row 347). | "counts with an ECAL photon that makes a π0 with the e+e− pair (M(e+e−γ) 0.119–0.173 GeV/c²)" |
| 101 | vertex-finder window | −170 < z < −30 mm, r < 20 mm (HParticleVertexFind setCut(−170, −30, 20)) (II:3183–3185; III:6672). | "the vertex finder looks for at least two tracks in −170 < z < −30 mm" |
| 102 | front (halo) region, U of Part II | −200 ≤ vz < −170 mm (II:3476). | "counts with the vertex 170–200 mm upstream" |

---------------------------------------------------------------------------------------------------------------------------------

## 5. Simulation: truth classes, templates, groups (legend names)

Classification code: truth:216–263 (lepton pair, pion pair, four-particle class), truth:439–455 (template), defs:303–363 and
s5prep:322–333 (groups). Part I tables: I:1213–1246; I5:43–76.

| # | label | definition (source) | short legend name | longer name |
|---|---|---|---|---|
| 103 | class 0 / 1 (Part I) | four-particle class 0: a misidentified particle; 1: a pion decayed into a muon (truth:260; I:1232) | "misidentified" / "π → μ" | "a misidentified particle in a slot" / "a pion that decayed into a muon" |
| 104 | class 2 = C2 = Gsig2 | η → π+π−π0, π0 → e+e−γ (Dalitz), photon lost: bump (truth:260; defs:303, 320; s5prep:325). Physics weight: KLOE-2 Dalitz plot; per-count QED γ* angle weight at the bump edges (D38). | "η→π+π−π0, π0 Dalitz" | "η → π+π−π0 with the π0 Dalitz decay π0 → e+e−γ, the photon not detected (bump)" |
| 105 | class 3 = C3 = Gsig3 | η → π+π−π0, a π0 photon converted into e+e− in the material, the other photon lost: bump (truth:260; defs:303). | "η→π+π−π0, converted photon" | "η → π+π−π0, π0 → γγ with one photon converted into e+e− in the detector material (bump)" |
| 106 | class 4 = C4 = Gsig4 | η → π+π−γ with the photon converted: peak (truth:261; defs:303). | "η→π+π−γ, converted photon" | "η → π+π−γ, the photon converted into e+e− in the material (peak)" |
| 107 | class 5 = C5 = Gsig5 | direct η → π+π−e+e− (η set, physics weights V3): peak (truth:261; defs:304). | "η→π+π−e+e− (direct)" | "the direct decay η → π+π−e+e− through a virtual photon (peak)" |
| 108 | class 6 = C6 | η → π+π−X17, X17 → e+e− (η set), weight 2 per simulated entry = per unit BR, never in sums (truth:261; defs:304, 321). | "X17 (17 MeV/c²)" | "hypothetical X17 boson in η → π+π−X17, X17 → e+e−, per unit branching ratio" |
| 109 | class 7 / 8 (Part I) | 7: pions with a correlated e+e− pair from one source ("multi-pion"); 8: leptons from different sources ("CB-like") (truth:262). | "pions + true pair" / "two sources" | as named |
| 110 | C1 | not a class in the reports: only "spec. Table C1 / Fig. C1" (§12). | – | – |
| 111 | C2c (III:4451, 4831) | C2 weighted with its per-count QED γ* angle weight (D38). | "η→π+π−π0, π0 Dalitz (angle corrected)" | "C2 with the QED decay angle of the virtual photon restored per count" |
| 112 | E23 (III X17 model) | smoothed C2 + C3 term of the X17 fit (III:5927–5929, 6755). | "η bump (π0 Dalitz + conversion)" | "the two bump classes together, smoothed" |
| 113 | L0–L11 (lepton-pair class, I) | L0 misidentified; L1 / L2 same π0 / η Dalitz decay; L3 same conversion; L4 ρ0, ω, φ → e+e−; L5 ω → π0e+e−; L6 Dalitz decay of the π0 of ω → π0e+e−; L7 two different Dalitz decays; L8 Dalitz + conversion; L9 other different sources; L10 leptons of a forced η → π+π−e+e−; L11 X17 leptons (truth:217–219; I:1218–1223). | spell out | spell out |
| 114 | P0–P9 (pion-pair class, I) | P0 not a pion; P1 pion decayed into a muon; P2 / P3 / P4 / P5 / P6 both pions from the same η / ω / ρ0 / φ / K0S; P7 both primary (SMASH final state); P8 mixed or secondary; P9 pions of a forced η → π+π−e+e− (truth:242–243; I:1226–1229). Not "P3" of the flexible fits, not spec. P1/P2. | spell out | spell out |
| 115 | T1 | signal: four-particle classes 2–6 (also with a ghost lepton track) (truth:440–441, 452; I:1238). Compared sum uses classes 2–5 only. A count with the same particle in both lepton slots is T3 even when its class is 2–6 (TemplateClass, truth:450–451). | "η signal" | "the η signal channels (and X17) with the four true particles" |
| 116 | T2 | physical background: a true e+e− pair from one source (L1–L6; η set also L10, L11 outside the signal) with any pions (misidentified pions and π → μ included), no fake lepton (truth:442–443, 454; I:1239–1240). | "physical background" | "true e+e− pair from one decay or conversion, with any two pions" |
| 117 | T3 | fake lepton: L0, a ghost track in a lepton slot, or the same particle in both slots (this last case takes priority over the signal); not in the compared sum (truth:444–445, 451–453; I:1241). | "fake lepton" | "a lepton slot holds a misidentified or ghost track" |
| 118 | T4 | the two leptons from different sources (L7–L9); shown apart (truth:446–447; I:1242–1243). | "two sources" | "the two leptons come from different decays" |
| 119 | T12 (III) | T1 (C2–C5) + T2, never X17 (defs:307, 324). | "SMASH total" | "SMASH signal + physical background" |
| 120 | G_η (G_eta; Geta in code) | template-2 counts with both pions from one η (P2) or an η Dalitz pair (L2), first match (defs:305, 330; s5prep:323; I5:64–66). | "other η decays" | "η produced but not a signal combination (e.g. η Dalitz pair or η pions with another pair)" |
| 121 | G_π0D | remaining template-2 counts with a π0 Dalitz pair (L1) (defs:305, 331; I5:67). | "non-η π0 Dalitz" | "a π0 Dalitz e+e− pair with any two pions, not from an η" |
| 122 | G_conv | remaining template-2 counts with a conversion pair (L3) (defs:305, 332; I5:68). | "non-η conversions" | "an e+e− pair from a converted photon, with any two pions, not from an η" |
| 123 | G_vec | remaining template-2 counts with L4–L6 (ρ0, ω, φ → e+e−; ω → π0e+e−; second-decay π0 Dalitz) (defs:306, 333; I5:69–70). | "vector mesons" | "ρ0, ω, φ → e+e− and ω Dalitz pairs, with any two pions" |
| 124 | G_oth | all other template-2 counts (in practice the e+e− pair of a direct η → π+π−e+e− decay with pions from elsewhere) (defs:306, 334; I5:71–72). Code check (v2): template 2 holds lepton classes L1–L6, L10, L11 (truth:453–454); after G_η (P2 or L2), G_π0D (L1), G_conv (L3) and G_vec (L4–L6) only L10 and L11 remain, and L11 (X17) entries weigh 0 (sim_truth.h Opt brX17 = 0; SimE:469–470). So G_oth = leptons of a forced η → π+π−e+e− whose pions are not both from that η (foreign, misidentified or decayed pions): it is of η origin. | "other true pairs" | "other single-source pairs (mostly a direct-decay pair with foreign pions)" |
| 125 | G_D (III; GD in figures) | G_π0D + G_vec + G_oth (defs:307, 324, 352–355; III:1214, 6687). | "Dalitz and vector-meson pairs (no conversion)" (do not write "non-η": its G_oth part is of η origin, row 124) | "SMASH background without a conversion pair: π0 Dalitz, vector-meson and other single-source pairs (the last mostly direct-decay η pairs with foreign pions)" |
| 126 | G_C (III; GC in figures) | = G_conv (defs:307; III:6687). | "non-η conversion pairs" | "SMASH background with a conversion pair" |
| 127 | G_D(ω/ρ), G_D(other); suffixes _om / _nom (III flex) | G_D split into counts whose two pions come from one ω or ρ0 (pion class P3 or P4) and the rest (defs:305–306, 336; III:6734). | "π0 Dalitz etc., pions from ω/ρ0" / "… other pions" | as named |
| 128 | G_fake, G_diff (I) | template 3 and template 4 (variants only) (I5:73). | "fake lepton" / "two sources" | as T3 / T4 |
| 129 | multi-pion (I, Is, III) | G_π0D + G_conv + G_vec (I:1801; I5:74; code spec5:27); Part I short: "π0 Dalitz and conversion pairs with pions (no η)" (Is:424, 675). G_vec is 0.107 of the 659.9 nb after "pre" (I5:88–90), so the short wording is a harmless simplification (§16 item 11, settled). | "multi-pion background" | "π0 Dalitz, conversion and vector-meson pairs with pions, no η" |
| 130 | Gsig2–Gsig5 (I figure legends, s5_spec) | groups of the Part I blind test: classes 2, 3, 4, 5 (s5prep:325–326). | as classes 2–5 | – |
| 131 | ghost | reconstructed track flagged by PAT as not matching one GEANT particle (sim_isghost) (truth:48; I:1802). | – | "ghost track" |
| 132 | variant (a), ghostA (III) | ghost-flagged true pairs moved from T3 into T2 (and its groups) (defs:308, 341; III:6688). Not display (a), not scenario (a), not option (a) (§14). | – | "ghost-flagged true pairs counted as physical background" |

---------------------------------------------------------------------------------------------------------------------------------

## 6. Simulation weights, ownership, mass transfer and their variants

| # | label | definition (source) | plain phrase |
|---|---|---|---|
| 133 | std, lep, η set | three simulation sets of the same 10^10 SMASH events: natural decays; every π0, η, ρ0, ω, φ forced into dileptons; every η forced into π+π−e+e− or π+π−X17 (I:1783; III:6680). | "standard, lepton and η simulation sets" |
| 134 | mode A / mode B | ownership: mode A (default) the lepton set provides every count with a forced-decay product of a π0, η, ρ0, ω or φ (the standard set drops them), the η set every count with a forced-η product; mode B leaves the π0 and η Dalitz counts in the standard set (I:1024–1028; III:1197–1199, 6682). | "each simulated count taken from one set only" |
| 135 | w, W, w_db | Part I notation: w = generator weight stored on a track by the converter (n_ch × database BR on forced-decay products, 2 on X17 products; I:829–832, 1794–1795); W = physical weight of a count from the rules v2 with PDG branching ratios (I:1007–1028, 1795). Code and Part III: w = the count's physical weight with PDG BR, i.e. Part I's W (SimE row fields w = sim_w, w_db = sim_w_db; III:1199, 6682; D7 "Σw/221.05"); w_db = the same with the converter-database BR (variant, family S4). See §16 item 21. | "weighted to the PDG branching ratios" |
| 136 | rules v1 / v2 (I) | weight rules; v2 uses the exact HGeant branching ratios of the card ETA (I:1796). | – |
| 137 | pre-rw (I) | numbers before the physics reweighting of the η decays (I:1803, 1760). | "before the η-decay weights" |
| 138 | V0–V4 (I; phys) | direct-decay weight variants: V0 = Kroll–Wada × p*³ λ^{3/2} (no form factors, no angles); V1 = V0 × abs(F_V)²; V2 = V1 × abs(F_VMD)²; V3 = V2 × decay-angle factors (main); V4 = V3 × plane term R_φ (check) (phys:64–67). Not spec. Fig./Table V1 (§12), not PID file v1. | "physics weights of the η decays" |
| 139 | V3, code 108 (I, III) | kPhysMain = 108: direct decay V3 + KLOE-2 bump Dalitz plot + conversion abs(F_V)² photon spectrum + background rules (phys:156; defs:365; III:1199–1202, 6683). | "η-decay physics weights (QED virtual-photon spectrum with form factors and angles for η → π+π−e+e−, measured Dalitz plot for η → π+π−π0, pion form factor for η → π+π−γ)" |
| 140 | V2, code 107; GS, code 236; 116; 100; 76; 0 | variants: V2 direct decay without decay angles (107); F_V Gounaris–Sakurai (236); bump Dalitz plot 1 + aY (116); bump phase space (100); conversion F_V off (76); physics weights off (0) (defs:365, 367–387; phys:157, 164). | spell out the variant |
| 141 | variant tags nom, V2, bump1aY, bumpPS, convOff, GS, wtsOff, modeB, wdb, ghostA, sK0, sUp, sDn, r089, r110 (code only) | defs:367–389 (SimVarText 374–380): nominal; …; s = 0.99521 (K0S-based); s = 0.99705; s = 0.99525; r = 0.89; r = 1.10. | spell out |
| 142 | s, r_res, M4' (III; "mass transfer", D8) | M4' = s [m4_true + r_res (m4 − m4_true)], s = 0.99615, r_res = 1 (entries without a true mass: s m4); M(π+π−)' = 0.9952 M(π+π−); M(e+e−) not scaled (defs:85–95; III:1206–1210, 6689). r_res is written r in package tables. | "simulated masses moved to the data mass scale (M(π+π−e+e−) × 0.99615, M(π+π−) × 0.9952)" |
| 143 | L_sim, L_data, L, L_SF, L_an, L_LGAD | L_sim = 221.05 nb⁻¹ per set; L_data = 6080.5 nb⁻¹ (Part I, all runs, then unconfirmed); L = L_an = 6079.2 nb⁻¹ (17 557 analysed runs, START-free, ±6.9 %); L_SF = 6080.5 (all 17 615 runs); L_LGAD = 5587.5 (LGAD START chain) (defs:42–51; II:3467–3469; I:1785). | "luminosity of the analysed runs, 6079.2 nb⁻¹" |
| 144 | per-count C2 weight (III; D38) | QED/isotropic weight of the γ* decay angle of a simulated C2 entry, averaged over the circle of its lost photon; effective 0.861 (u) / 0.860 (per slice); replaces the constant factors of D35 (0.879 / 0.882 for M(e+e−) < 10 / ≥ 10 MeV/c² per slice, 0.879 / 0.8825 in u; 0.881 effective; III:3911, 6610) (III:6728; IIIs:1366). | "simulated π0 Dalitz decays weighted to the QED decay angle of the virtual photon" |
| 145 | peak tail, t; family "peak tail" (III; D40) | low-mass tail of the simulated η peak (combinations whose leptons' momentum loss shifts M4 by more than 5 MeV/c²), scale t (t = 1 SMASH, t = 0 no tail), systematic t = 1 ± 0.5 (III:6730). Not t_i of the X17 model. | "the simulated low-mass tail of the η peak (lepton radiation in the material)" |
| 146 | k_cross, cross-photon pairs (III; D39) | lepton pairs from the two different photons of the π0 in η → π+π−π0; their scale in N relative to the CB rule (III:1103–1104, 6729). | "lepton pairs from the two photons of one π0" |

---------------------------------------------------------------------------------------------------------------------------------

## 7. Signal measures, estimators and fits

### 7.1 Part I (simulation blind test; I5:930–945, I:1840–1841)

| # | label | definition | plain phrase |
|---|---|---|---|
| 147 | M1 (I) | Gaussian fixed at 0.5460 GeV/c², σ = 6.7 MeV/c², + second-order polynomial, calibrated with the response to the half-A peak template. Not stage-A model M1 of Part III (§14). | "Gaussian + polynomial fit of the peak (calibrated)" |
| 148 | M2 (I) | template fit: free normalisations of the peak, bump and background templates, binned Poisson likelihood with the template variance added. | "template fit" |
| 149 | M2b (I) | M2 with the background split into G_η and the rest. | "template fit with the other-η background free" |
| 150 | M2d (I) | M2 with the peak split into class 5 (free) and class 4 (constrained to its half-A value ± MC error). | "template fit of the direct decay, conversion constrained" |
| 151 | M3 (I) | linear side bands 0.5255–0.5333 and 0.5587–0.5665 GeV/c² (simulation scale; 7.8 MeV/c² on each side of the simulation peak window; b4I:61–63; s5_blind_ana_dev.C:11), signal leakage λ corrected: S = [N_win − (N_L + N_U)/2 · 25.4/7.8]/(1 − λ). | "linear side-band subtraction" |
| 152 | M3c (I) | M3 with the curvature of the background taken from half A. | "side bands with the simulated curvature" |
| 153 | B1 (I) | bifurcated Gaussian + η Gaussian + polynomial (parametric bump fit). Not spec. Table B1 (§12). | "parametric fit of the bump" |
| 154 | B2 (I) | the bump of M2. | "bump of the template fit" |
| 155 | B3 (I) | edge step in M(π+π−e+e−). | "step at the kinematic edge" |
| 156 | B4 (I) | side band 0.5665–0.70 GeV/c² (simulation scale; it starts where the upper M3 side band ends; about 0.564–0.697 on the data scale) with the background shape from half A (b4I:64; I5:940). The same idea as B4 of Part III with another side band and another background shape (§16 item 3). | "bump from an upper side band with the simulated background shape" |
| 157 | D5 (I, method) | direct-decay yield from a split at M(e+e−) = 10 MeV/c² with a 2 × 2 unfolding. Not decision D5 (§14). | "direct decays from the M(e+e−) split" |
| 158 | X1 (I, estimator) | fit of M(e+e−) in the window with free normalisations of the direct decay and the rest. Not cut set X1, not spec. Fig. X1. | "fit of the M(e+e−) spectrum" |
| 159 | X2 (I, estimator) | X1 with the raw half-A template. | – |
| 160 | X3 (I, estimator) | X1 on 8–40 MeV/c². | – |
| 161 | X4 (I, estimator) | X1 with the direct-decay normalisation constrained to 1 ± 0.10. | – |
| 162 | X5 (I, estimator) | counting with B from half A. | "simple count" |
| 163 | half A / half B (I) | simulated events split by event-group number (even/odd): A gives the templates, B the pseudo-data (I:1837). Not mode A/B, not class A/B (§14). | "templates from one half of the simulation, pseudo-data from the other" |
| 164 | S, B, P (purity), F, F_best, Q, UL, UL_stat, s_1 (I) | S signal and B other T1+T2 counts in a window at L_data; purity P = S/(S+B); F = S/√(S+B); F_best its maximum; Q = ε_S/√ε_B; UL expected 90 % CL limit with the MC error of B, UL_stat without; s_1 X17 counts per unit BR (I:1833–1842). | spell out |
| 165 | plateau, freeze rule (I) | cuts with F ≥ 0.95 F_best on the full sample and both halves; among them the fewest active cuts, then the best worst-neighbour F (I5:622–626). | – |

### 7.2 Part II (data only; II:3521–3543, 1553–1640)

| # | label | definition | plain phrase |
|---|---|---|---|
| 166 | P_main (II, III, shorts) | χ² fit of N − CB with a frozen Gaussian (mean 0.5439 GeV/c², σ 5.4 MeV/c²) + pol2 in 0.5225–0.6100 GeV/c² (Scan:84–103; II:3522, Eq. at II:1564). | "fit of the η peak with a fixed Gaussian (0.5439 GeV/c², width 5.4 MeV/c²) plus a second-order polynomial" |
| 167 | S, stat, model | S: content of the whole Gaussian of P_main; stat: fit error; model: rms of S over the edge-safe grid of 144 fit variations (polynomial order, fit range with lower bound ≥ 0.5225 GeV/c², width) with χ² probability p > 0.05 (II:3523; IIs:818); old envelope: half range of S over the 48 variations of the chapter tables with p > 0.05: polynomial order 1–4 × the four fit ranges 0.520–0.590, 0.515–0.610, 0.510–0.620 and 0.505–0.650 GeV/c² × width 4.9 / 5.4 / 5.9 MeV/c² (Scan:125–132, 944–990; II:3524 names only the lower bounds). Not S (signal) of Part I. | "η-peak content of the fit" |
| 168 | P_pois (II) | Poisson-likelihood fit of the unsubtracted N (control) (II:3527). | "Poisson fit without subtraction" |
| 169 | E, h (II; III "per-slice edge") | edge-step visibility: sum of the step heights h in the four M(π+π−) slices (row 56), counts per 2.5 MeV/c², not a yield; fit pol2 [+ η Gaussian] + h × erfc step at the fixed edge, smeared 5.4 MeV/c² (II:1596–1607, 3528). Not E1–E10, not E(m,q). | "summed height of the step at the kinematic edge of the bump in four M(π+π−) slices" |
| 170 | u, h_u, σ_u (II) | u = M4 − M_edge(M(π+π−), M(e+e−)); h_u step height at u = 0 ("pre" + vertex, 0.30 ≤ M(π+π−) < 0.40), counts per 2.5 MeV/c²; σ_u = 6.7 MeV/c² or free (II:3529–3530). | "height of the step at the edge in u" |
| 171 | T, m5c (II; III) | ECAL-tagged η → π+π−π0 content (row 100). Not the beam energy, not T1–T4. | "η → π+π−π0 counts with the lost photon found in the ECAL" |
| 172 | D, S_SB, r_SB, f_pk, f_SB (II) | D = (N−CB)_Rpk − r_SB (N−CB)_RSB; S_SB = D/(f_pk − r_SB f_SB); r_SB ratio of the P_main pol2 integrals (0.9163 for "pre"); f_pk = 0.9426, f_SB = 0.0574 (II:1627–1637, 3542). Not decision Dn, not D of the X17 model. | "η-peak counts after a side-band subtraction" |
| 173 | W, N_b, Q, gain_E, gain_T (II) | W = N − CB in the bump window; N_b counts there; Q = E/√N_b; gain = kept share of E or T / √(kept share of N_b), relative to "pre" (II:3531, 3539). | spell out |
| 174 | S_fake, flags G, M (II) | S_fake: Gaussian content a cut makes from a smooth spectrum; flag G: like-sign turn-on in the fit range; flag M: abs(S_fake) > 92 counts (II:3540). | "sculpting test" |
| 175 | rules best-F, loose, pure, nominal; F_A, F_B, b(F) (II) | fixed rules that choose a working point on one half of the data; F on half A / half B; optimisation bias b(F) = ½[(F_A(c_A) − F_B(c_A)) + (F_B(c_B) − F_A(c_B))] (II:3541, 3543). | – |
| 176 | k, k_mix, k_ee; GM, AM (II) | k CB factor (convention 1); k_mix event-mixing estimate in M4; k_ee the same in M(e+e−); GM = 2√(N++N−−), AM = N++ + N−− (II:3515–3516). | "CB factor from event mixing" |
| 177 | Y, Y_N, Y_S, Y_T (II) | visible yield X × 64/63 / L in visible counts per nb⁻¹ (counts, not events; not a cross section) (II:3472). | "visible yield per luminosity" |
| 178 | ε_vtx, C_vtx, smooth like-sign CB (II, III) | ε_vtx = X_V/(X_V + X_Nv) for a measure X (S, T, E); C_vtx = ε_sim/ε_data = 1.006 ± 0.021, not applied; smooth like-sign CB: N++ and N−− fitted with pol3–pol5 in 12 range variants, CB = 2√(f++ f−−) (II:3460–3461; III:6674–6675). | "η efficiency of the event-vertex cut" |

### 7.3 Part III: stage A, stage B and the template fits

| # | label | definition | plain phrase |
|---|---|---|---|
| 179 | model (Eq. model) | μ = k_e CB + c_j (1 + b x)[G_D + κ_B G_C] + a_η [C2 + κ (C3 + C4) + ρ_d C5 + g G_η], x = (M4 − 0.55)/0.1 (III:1338–1345; defs:278). | "SMASH background with a free scale per M(π+π−) slice and a common tilt, plus the SMASH η classes with a common scale" |
| 180 | c_j, b, x, κ_B, a_η, κ, ρ_d, g, k_e (k_E1, k_E2) | background scale per slice; common tilt; x; conversion scale of the background; η scale; conversion scale of the signal (C3, C4); direct-decay scale; G_η scale; CB factor per M(e+e−) class (E1 / E2) (III:6714–6716, 6661). | spell out |
| 181 | stage A, B_A (III) | background fit in the control region with the η terms fixed (a_η = 0.4, κ = ρ_d = g = 1); free: c_1 … c_5, b, κ_B and the CB factor, k_E2 = 1 ± 0.05 (III:2056–2058); frozen: model c_j(1 + bx), class E2 only, CR-U 0.575–0.700 only; B_A its background prediction (III:2043–2072, 6717). | "renormalisation of the SMASH non-η background in the control region above the peak" |
| 182 | M0–M3 (III, stage-A models) | M0 scale c_j only; M1 + common tilt (adopted); M2 + quadratic; M3 per-slice tilts (III:2139–2153, 6589). Not the Part I methods M1–M3 (§14, §16). | "background scale only / with a common linear tilt / with a quadratic term / with a tilt per slice" |
| 183 | Step 0 (III) | the specification as written (E1 + E2, CR-U 0.565–0.700 + CR-K, data-only χ²) before the E1 fall-back (III:2073–2082). | – |
| 184 | stage B (III) | the template fits of the bump and the peak (§sec:fit, III:2476–2499). | "template fits" |
| 185 | F-main (III) | nominal template fit of F10 E2 (350 bins): a_η, ρ_d, c_1…c_5, b, κ_B free; κ = 1 ± 0.2, g = 1 ± 0.3, k_E2 = 1 ± 0.05; model M1 (III:2500–2504, 6718). Cross-check (failed the closure test). | "template fit of the ten M(π+π−) × M(e+e−) spectra (only M(e+e−) ≥ 10 MeV/c² fitted)" |
| 186 | F-free (III) | F-main with free α_2 … α_5 (one scale per η class) (III:2576, 6716). | "template fit with a free scale per η class" |
| 187 | F-gold (III) | fit of F5g (M(e+e−) ≥ 20 MeV/c²): α_2, α_5, g, c_j, b, k_E2 free, κ_B and aκ(C3 + C4) from F-main; gives BR_G (III:2568–2575). | "template fit of M(e+e−) ≥ 20 MeV/c²" |
| 188 | α_2–α_5, BR_G, BR_A, BR_C, R_34 (III) | free class scales; BR(η → π+π−e+e−) from F-gold (α_5/α_2), from ρ_d, from the conversion route (BR_C withdrawn; outcome of D15, III:6575); R_34 = α_4/α_3 (III:6716, 6742). | spell out |
| 189 | BB-lite, sandwich, GLS, hybrid GLS (III) | template variance in the fit weights (D13); error with the true data covariance; χ² with the template covariance (replicas); hybrid GLS: data covariance + exact diagonal of the raw templates + replica covariance of the smoothed ones (III:6736, 6747). | "χ² including the simulation-statistics covariance" |
| 190 | Y_c, Y_peak, Y_bump, Y_bump(edge), Y_2…Y_5 (III) | yields in raw counts × 64/63: class c; C2 + C3 in B; C4 + C5 in P; Y_bump(edge) = f_η(edge) × (C2 + C3)_B (III:1349–1352, 6720, 3917–3919). | "η-peak yield" / "bump yield" |
| 191 | a_P, a_B(E2), a_B(all) (III) | a_P = Y_peak/(C4 + C5)_SMASH in P; a_B(E2) = Y_bump(E2)/(C2 + C3)_E2; a_B(all) = Y_bump(all)/(C2 + C3)_all (III:6721). | "peak yield relative to SMASH" / "bump yield relative to SMASH" |
| 192 | DR, DR_F, DR_edge, DR_u (III) | DR = a_P/a_B(all) (side bands); DR_F from F-main; DR_edge = a_P/f_η(edge) (without the peak-composition term); DR_u = a_P/f_η(u-edge) (with it) (III:6722–6726). | "double ratio peak/bump, data over SMASH" |
| 193 | routes (direct-decay route, conversion route) (III) | η scale from C5 alone (direct decays, PDG BR) or from C4 alone (converted photons) (III:6741). | "η scale from the direct decays alone / from the converted photons alone" |
| 194 | K_p (III) | per-period normalisation factor C_sec(p)(1 − f_η) (III:6740). | – |

### 7.4 Part III side-band estimators (main results) and edges

| # | label | definition | plain phrase |
|---|---|---|---|
| 195 | B4 (III) | Y_bump = 64/63 {(N−CB)_B − [(N−CB)_U − η_U] B_A(B)/B_A(U) − [aκ C4 + aρ C5 + ag G_η]_B}, in E2, extrapolated to all M(e+e−) with f_E2 = 0.4773 (III:2765–2780, Eq. b4). Not determined. | "bump area: data minus CB in the bump window minus the background carried down from 0.575–0.70 GeV/c²" |
| 196 | η_U, R, f_E2 (III) | η content of U subtracted (D27); R = B_A(B)/B_A(U) (0.12735); f_E2 = (C2 + κC3)_E2/(C2 + κC3)_all in B = 0.4773 at κ = 0.832 (III:2775–2776, 2780–2781, 6737). | – |
| 197 | peak side bands (III) | Y_peak = 64/63 [(N−CB)_P − R′ (N−CB)_SB]/(1 − R′ f); R′ = ratio of [B_A + a(C2 + κC3) + ag G_η] in P and SB; f = (C4 + C5)_SB/(C4 + C5)_P = 0.1544 (III:2830–2846). Main peak result. | "η peak from the peak window minus the scaled side bands, corrected for the peak's leakage into them" |
| 198 | σ_tr (III) | closure term of the side bands from 13 displaced windows (III:6738). | – |
| 199 | f_η(edge), per-slice edge (III) | E_data/E_sim with the same per-slice fit (row 169) on the data and on SMASH T1 + T2; with or without the C2 angle correction (III:3893–3899, 6744). Cross-check since D37. | "η scale of the bump from the step at its kinematic edge, per M(π+π−) slice, data/SMASH" |
| 200 | f_η(tag) (III) | T_data/T_sim with the same fit (III:6744). Verification only. | "η scale from the ECAL-tagged bump, data/SMASH" |
| 201 | f_η(u-edge), h_data, h_sim (III) | step heights at u = 0 in the data and in SMASH C2 + C3 (one linear estimator); their ratio is the bump measurement (D37) (III:4345–4351, 6725). | "η scale of the bump from the step at its kinematic edge in u, data/SMASH" |
| 202 | κ (bump edge) (III) | scale of the simulated conversion pairs (C3) relative to the Dalitz pairs; 0.48 ± 0.11 from the peak routes, or free (III:3915–3917, 6732). | – |
| 203 | full reconstruction, T_exp, λ (ECAL; III) | bump combination with its ECAL photon found (m5c at the η mass, 0.49–0.62); T_exp = f_η × T_sim; λ = scale of the simulated cluster energies (1.285 data-like) (III:6731). | "bump counts fully reconstructed with the ECAL photon" |

### 7.5 Part III flexible SMASH fits (III:5438–5503)

| # | label | definition | plain phrase |
|---|---|---|---|
| 204 | ASIS | SMASH as is, no fit (III:5491). | "SMASH as it is" |
| 205 | P3 (flex grid) | Part III F-main reproduced (III:5492). Not pion class P3. | "the template fit of Part III" |
| 206 | F7 | 7 free scales: C2, C3, C4, C5, G_η, G_D, G_C (III:5493). | "one free scale per SMASH channel (7)" |
| 207 | F7b / F7bj / F7cj | F7 + common tilt / + per-slice tilts (alternative fit F7bj) / Part III background c_j(1+bx)[G_D + s_GC G_C] with free η classes (III:5494, 5498). | – |
| 208 | F8 (nominal) | 8 scales: F7 with G_D split into ω/ρ0 pions and the rest; k_E2 = 1 ± 0.05; S0, isBest ≥ 0, M(e+e−) ≥ 10 MeV/c², M4 0.40–0.70, five slices (350 bins) (III:5495, 5440–5449, 751). | "free scale for each SMASH channel (8 channels)" |
| 209 | F8b / F8bj; F9, F9b, F9bj | F8 + common / per-slice tilts; 9 scales (G_C split too), no / common / per-slice tilt (III:5496–5497). | – |
| 210 | s_X, λ (flex) | scale of component X; weight of the template noise in the fit weights (λ = 1 nominal, 0 Poisson) (III:6733–6735). In Fig. res-udatasim λ is a free factor on G_D + G_C (III:1107–1108). | – |

### 7.6 Part III X17 (III:5855–6300)

| # | label | definition | plain phrase |
|---|---|---|---|
| 211 | X-a, X-b, X-c, X-d | steps: model; validation on toys; freeze (r = 0 fit on the M(e+e−) side bands, x17_frozen.txt); unblinding (III:5881–5890; SPEC:234–250). | "the X17 procedure was fixed before the data were looked at" |
| 212 | r, BR_X BR_ee, r_UL, BR_dir | POI r = BR(η → π+π−X17) BR(X17 → e+e−)/BR(η → π+π−e+e−); BR = r × BR_dir, BR_dir = BR(η → π+π−e+e−) = 2.76e-4 (PDG 2026; frozen file 2.68e-4) (III:6752, 5921–5925). | "X17 strength relative to the direct decay" |
| 213 | N_5, N_4, n_j, D, C, X, R_X, s_X, s_X2D, s_1, q-scale, N5/N5ref | fitted direct and conversion peak yields; free background per M(e+e−) row; M(e+e−) shapes of C5, C4 and the X17 line; R_X = ε_X/ε_5 = 1.3455; log-normal yield nuisance; X17 counts per unit BR in X2D and in the counting window; variance factor 1.0299 (III:6753–6755); N5/N5ref (legend of the η-scale summary figure, III Fig. P1 and IIIs Fig. P1s): the η scale from the X17 fit, fitted N_5 over its SMASH value. | spell out |
| 214 | normalisation (1), (2), (3) | X17 yield tied to N_5 (primary, D30); to a_η (specification's fall-back); to a_P (III:6756). Not D31/D40 options. | "X17 yield tied to the fitted direct-decay yield" |
| 215 | CL_s, q̃_r, Asimov, FC | limit method, test statistic, median data set, Feldman–Cousins (III:6757). | "90 % CL upper limit (CL_s method)" |

### 7.7 Part III comparison after the typical cuts (III:1621–1700; not used for results)

| # | label | definition | plain phrase |
|---|---|---|---|
| 216 | displays (a), (b), (c); c_b | (a) SMASH as it is; (b) η channels × 0.40 and non-η background × c_b, c_b = [data − 0.40 η]/non-η in 0.575–0.70; (c) display fit a_η η + c non-η in 0.40–0.70 (60 bins), not a yield (III:1636–1641, 6709). | "SMASH as it is / with the η scaled to 0.40 and the background normalised above the peak / fitted" |
| 217 | δ (mismatch) | fraction of the model added in quadrature to every 5 MeV/c² bin error to bring χ²/ndf to 1 (III:1643–1646). Not the variant shift δ of the systematics. | "average mismatch of SMASH per bin" |
| 218 | proxy, Z, S/B gain, N_X | in the §e1 X17 window: proxy √N/(ε_X/ε_X(S0))/√N_S0; Z = N_X/√(N + var CB); S/B gain (ε_X/ε_X(S0))/(N/N_S0); N_X = 392 × C6(sel, window)/C6(S0, X2D) (III:1647–1655). | spell out |

---------------------------------------------------------------------------------------------------------------------------------

## 8. Closure (pseudo-data) configurations, scenarios and tests (Part III)

| # | label | definition (source) | plain phrase |
|---|---|---|---|
| 219 | BA | truth from half B of the simulation, templates from half A (the specification's gate) (III:2360–2361, 6749). | "pseudo-data from one half of the simulation, templates from the other" |
| 220 | AB | the swap of BA. | same |
| 221 | AA, BB ("method") | truth and templates from the same half. | "pseudo-data and templates from the same simulated events" |
| 222 | RF ("realistic MC") | truth from the full sample, toy t fitted with the templates of bootstrap replica t mod 100. | "pseudo-data with the template noise of the real analysis" |
| 223 | RB | RF with a bootstrap bias correction (2v − mean of 20 child replicas). Adopted gate: AA, BB + RF (D22). | "… with a bias correction" |
| 224 | scenarios (a)–(i) | truth relative to the fit model: (a) background × 1.35, tilt −0.15, a_η = 0.4; (b±) ρ_d 0.5/1.5; (c±) κ 0.8/1.2; (d±) g 0.5/1.5; (e) ω/ρ0-pion background × 0.3, rest × 1.6; (f) s + 0.0009, r_res 0.89; (g) a_η = 0; (h) a_η = 1; (i) CB like the per-cell geometric mean (III:2355–2358). | spell out |
| 225 | attenuation | bias of a fit towards zero when its signal templates are noisy (III:6748). | "a fit with noisy simulated signal shapes comes out low" |
| 226 | A1–A5 (stage-A tests) | A1 simulated C2–C5 < 1 % in every control-region bin; A2 null edge step and zero tagged η in CR-K; A3 injected η × {0, 0.4, 1, 2} moves c_j, b by < 0.3σ; A4 per-slice tilts consistent with a common b, and CR-U split at M(π+π−) = 0.42; A5 K0S veto 0.485–0.507 moves the parameters by < 1σ (III:2197–2214). Not the defs.h comments "A2 note", "A3 note", not spec. Table A2. | spell out |
| 227 | fall-back 1, fall-back 2; configurations 1–3 (III) | configuration 1 = specification (CR-U 0.565–0.700 + CR-K); 2 = fall-back 1 (CR-U from 0.575 + CR-K with M(π+π−) ≥ 0.44); 3 = fall-back 2 (CR-U 0.575–0.700, no CR-K) (III:2086–2090). | – |
| 228 | null tests S(++), S(−−), S(Σ), E_ls (II) | frozen peak fit and edge-step sum applied to N++, N−−, N++ + N−− (II:3534). | "the same fit applied to the like-sign counts (should be zero)" |
| 229 | toy, bias, pull, coverage, stress test | one pseudo-data set; mean (extracted − true); (extracted − true)/error; fraction of toys within the error (I:1838–1839; I:1846). | – |

---------------------------------------------------------------------------------------------------------------------------------

## 9. Systematic families and error labels

### 9.1 Families S1–S13 of the Part III quantities (tS1; III:3061–3100, 6638)

| # | label | content | plain phrase |
|---|---|---|---|
| 230 | S1 (family) | CB: k_E2 1 ± 0.05 / fixed 1 / 1.032; k_E1 free / 1; like-sign sum; 20 MeV/c² cells; per-cell geometric mean (SPEC:282) | "combinatorial background" |
| 231 | S2 (family) | control-region definition: CR-U lower edge 0.560 / 0.565 / 0.575, upper edge 0.65 / 0.70 / 0.75; CR-K 0.425 / 0.43 / 0.44; K0S veto (SPEC:283); for the peak also the side-band width (III:3110) | "control-region definition" |
| 232 | S3 (family) | background model: c_j(1 + bx); + quadratic; per-slice tilt; exp(bx); common c + c1 y; ω split with the range to 0.85; κ_B = κ (SPEC:284); D11data, D11data + ω (III:2155) | "background model" |
| 233 | S4 (family) | background template: exp(pol3) / exp(pol4) / raw BB-lite; ghost variant (a); mode B; w_db (SPEC:285); hybrid-GLS F-main (III:3108) | "background shape from SMASH" |
| 234 | S5 (family) | range and binning: 0.40–0.70 / 0.42–0.65 / 0.45–0.70 / 0.40–0.75 GeV/c²; 2.5 MeV/c² bins; slice edges ± 5 MeV/c²; M(e+e−) split at 8 / 12 MeV/c² (SPEC:286) | "fit range and binning" |
| 235 | S6 (family) | mass transfer: s ± 0.0009; s = 0.99521; r_res 0.89 / 1.10; M(π+π−) scale ± 0.0007; s free (SPEC:287) | "mass scale and resolution of the simulation" |
| 236 | S7 (family) | physics weights: direct decay V3 / V2 / V4; bump KLOE / 1 + aY / phase space; conversion F_V on / off, GS (weights off shown only) (SPEC:288) | "η-decay physics weights" |
| 237 | S8 (family) | G_η: prior width of g 0.30 → 0.15 / 0.50; g free (SPEC:289); widened in the review to the observed range of g, 0.47–2.6 (III:3067) | "scale of the other η decays" |
| 238 | S9 (family) | signal templates: raw / factorised f(ee) g(slice) h(M4 given slice) (SPEC:290) | "simulated signal shapes" |
| 239 | S10 (family) | closure: abs(bias) of scenario (a) in RF (WP4); for the peak σ_tr of the 13 displaced windows (4.7 % of Y_peak); for the bump max(abs(N − pred), σ_pred) in CR-L (37.5 % of Y_bump(E2)) (SPEC:291; III:2464, 2866, 2888) | "pseudo-data closure" |
| 240 | S11 (family) | selection: S0 / "pre"; pions in the vertex fit / Δz_π < 40 mm; vz > −150 / −170 mm (SPEC:292); counted only for efficiency-corrected quantities and only if its shift exceeds twice its noise (III:3043–3044) | "selection" |
| 241 | S12 (family) | efficiency model: pion-vertex efficiency ± 2.5 %; ε2 QED at M(e+e−) ≥ 20 MeV/c²; (p_T4, y4) reweighting; X17 kinematics (SPEC:293); for σ_η also the composition of the peak from the data, ± 31.8 % (III:3066–3067) | "efficiency model" |
| 242 | S13 (family) | peak tail t = 1 ± 0.5 (D40). Not in the specification (S1–S12, SPEC:280–293): added in the fourth integration as row S13 of tables_all/p3w_tableS1_v3.tex (III:3075–3078; IIIs:893, 912) | "low-mass tail of the η peak" |

### 9.2 Families E1–E10 of the per-slice bump edge (III:3934–3957; "family E1" etc.)

| # | label | content | plain phrase |
|---|---|---|---|
| 243 | E1 (family) | edge model | "edge-fit model" |
| 244 | E2 (family) | M(π+π−) slices (fit range from 0.29 etc.) | "M(π+π−) slices" |
| 245 | E3 (family) | mass transfer, resolution | "mass scale and resolution" |
| 246 | E4 (family) | CB | "combinatorial background" |
| 247 | E5 (family) | physics weights (η → 3π Dalitz plot) | "physics weights" |
| 248 | E6 (family) | selection | "selection" |
| 249 | E7 (family) | step of a smooth background (symmetric envelope 6.0 %) | "step faked by a smooth background" |
| 250 | E8 (family) | C2 factor, model allowance ± 0.004 (D38) | "QED angle correction of the π0 Dalitz decays" |
| 251 | E9 (family) | κ of the peak (± 0.11) | "conversion scale" |
| 252 | E10 (family) | peak tail t = 1 ± 0.5 (D40) | "low-mass tail of the η peak" |

### 9.3 Families U1–U11 of the u-edge (tU; III:4765–4800, 6727)

| # | label | content (tU) | plain phrase |
|---|---|---|---|
| 253 | U1 | fit model (polynomial order, range, simulation weights, template evaluation) | "fit model" |
| 254 | U2 | resolution: step smearing and the data/simulation resolution ratio (r = 0.89 / 1.10, own smearing per sample) | "resolution" |
| 255 | U3 | M(π+π−) range | "M(π+π−) range" |
| 256 | U4 | mass scales and edge function (data edge and M(π+π−) scales, simulated s, M(π+π−) transfer) | "mass scales" |
| 257 | U5 | CB (normalisation k, cells, rule, slices, like-sign shape) | "combinatorial background" |
| 258 | U6 | u binning | "binning" |
| 259 | U7 | physics weights (η → 3π Dalitz plot), peak form factors, ownership | "physics weights" |
| 260 | U8 (and U8′) | peak term (a_P, its template) and its low-mass tail t = 1 ± 0.5 (D40; family maximum); U8′ = the tail family in an earlier draft | "η-peak term and its tail" |
| 261 | U9 | selection (pre vs S0) | "selection" |
| 262 | U10 | C2 γ* angle: model allowance of the per-count QED weights (± 0.004, D38) | "QED angle correction" |
| 263 | U11 | step projection of the smooth background (pseudo-data envelope) | "step faked by a smooth background" |

### 9.4 Error labels and rules

| # | label | definition | plain phrase |
|---|---|---|---|
| 264 | stat, syst, MC, N | statistical (sandwich: event groups + CB replicas); systematic (families, quadrature); template statistics (100 bootstrap replicas of the simulation); normalisation (L 6.9 %, period 5 %, BR, C_sec(1 − f_η)) (III:6745). | "statistical / systematic / simulation-statistics / normalisation" |
| 265 | MC(Δ) | template noise of a variant shift (replica spread of the difference), subtracted (D24) (III:6746). | – |
| 266 | R1 (rule) | WP7 rule: a variant rejected by the data is listed, not counted (III:6746). Not region R1. | "variants that the data reject are not counted" |
| 267 | WP7 rules, "literal", Barlow | δ² = max(0, Δ² − abs(σ²_stat,v − σ²_stat,0) − MC(Δ)²); a ± pair counts the mean of its two shifts; family = largest member; families in quadrature; "Barlow, stat only" = without MC(Δ) and R1 (III:3029–3060, 3903–3909). | – |
| 268 | D11data, D11data + ω (S3 variants) | the data-only χ² chain of D11 (per-slice tilts + quadratic (+ ω split)) (III:2153–2155). | – |
| 269 | period normalisation | ± 5 % on integrated yields per L (rate steps between beam periods) (III:6669). | "± 5 % normalisation between beam periods" |
| 270 | one-sided items (live time, sorter, cor_ok, C_vtx) | not applied; live time +7.3/−0 % (D17); START-sorter loss; cor_ok loss; C_vtx +2 % (III:1119, 6668). | "corrections not applied, quoted one-sided" |

---------------------------------------------------------------------------------------------------------------------------------

## 10. Decisions, amendments, fall-backs, status words

Full table: III:6513–6616 (tab:decisions); condensed: IIIs:1286–1330. Short reports may carry at most ONE sentence pointing to
the decisions appendix (I1); the plain phrases below are for the long reports, or to replace a D-number in a short report by its
content.

| # | label | content | status (15 Nov 2026) |
|---|---|---|---|
| 271 | D1 | nominal selection S0 = "pre" + both pions in the vertex fit; "pre" in parallel | applied |
| 272 | D2 | isBest ≥ 0 nominal; isBest = 1 only for data-only ratios | applied |
| 273 | D3 | no upper cuts on MM4, M(π+π−), the η-rest angle or OA | applied |
| 274 | D4 | DIRECT a variant; PEAK only after a fake-peak test, simulated MM4 at 4500 MeV | applied |
| 275 | D5 | CB = 2√(N++N−−), k = 1, never from simulation; k_E2 = 1 ± 0.05, k_E1 free | applied; fall-back fired (k_E1 = 13.2 in stage A, 7.6 in stage B): E1 dropped; the X17 fit uses k_E1 = 1 ± 0.5 (III:6555–6556) |
| 276 | D6 | simulation std + lep + η, owned entries (mode A), w × V3 (108), T1 + T2; ghost variant (a) as systematic | applied; ghost ruling open |
| 277 | D7 | data unit (N − CB) × 64/63 / L, L = 6079.2; simulation Σw/221.05; old tables × 1.0161 | applied |
| 278 | D8 | mass transfer s = 0.99615, r_res = 1; M(π+π−) × 0.9952 | applied |
| 279 | D9 | control regions in the M4 × M(π+π−) plane (CR-U, CR-K); CR-L a test | fall-back: CR-K dropped, CR-U from 0.575 |
| 280 | D10 | only the non-η background renormalised; g = 1 ± 0.30 | applied |
| 281 | D11 | background c_j(1 + bx); quadratic, then per-slice tilts only if Δχ² > 4 | applied with the template covariance (deviation) |
| 282 | D12 | main extraction = 2D Poisson template fit (F-main) | failed closure: cross-check |
| 283 | D13 | non-η templates smoothed exp(pol3); signal templates raw with BB-lite | applied |
| 284 | D14 | side bands next to the fit; P_main, edge, ECAL tag as cross-checks | side bands became the main result |
| 285 | D15 | BR(η → π+π−e+e−): golden ratio, then ρ_d, then conversion route | F-gold failed closure: BR_G a cross-check; BR_C withdrawn; no BR estimator validated (III:6575–6576) |
| 286 | D16 | σ(pp → ηX) secondary, needs the normalisation items | done; PRELIMINARY |
| 287 | D17 | live time not applied, +7.3/−0 % one-sided | applied |
| 288 | D18 | X17: 2D likelihood, background free per M(e+e−) row, 90 % CL_s | applied |
| 289 | D19 | blinding of the X17 M(e+e−) shape 13–27 MeV/c² until the freeze | blind in the code, not in the information |
| 290 | D20 | no yield before its estimator passes the closure test | applied |
| 291 | D21 | Amendment 1: η yields × (1 − f_η), f_η = 3 (+4 −3) %, never × (1 − f_win) | approved |
| 292 | D22 | closure gate: method configurations + RF instead of BA | approved |
| 293 | D23 | D11 and the stage-A χ² criterion with the template covariance | approved |
| 294 | D24 | systematics rules MC(Δ) and R1 | approved |
| 295 | D25 | PDG 2026 everywhere, also for the X17 limit | ruled |
| 296 | D26 | sector-2 deficit span 034–045 for C_sec | approved |
| 297 | D27 | B4: η_U subtracted before the extrapolation | approved |
| 298 | D28 | σ(pp → ηX) labelled PRELIMINARY | ruled: label stays |
| 299 | D29 | production-model item (+19/−7 %) kept separate from the total | ruled |
| 300 | D30 | X17: k_E1 = 1 ± 0.5, M(e+e−)-dependent background, side bands; normalisation (1) primary | approved |
| 301 | D31 | make-up of the η peak from the data (routes) instead of κ = 1 ± 0.2 | ruled, option (a), provisional |
| 302 | D32 | data lose close e+e− pairs more than the simulation: to be tested in a separate analysis | ruled |
| 303 | D33 | vertex study: smooth like-sign CB; S, T, E independent; C_vtx not applied | approved |
| 304 | D34 | typical cuts not applied; shown as comparison only | ruled |
| 305 | D35 | per-slice edge scale with the C2 γ*-angle correction recommended | approved |
| 306 | D36 | flexible SMASH: F8 with three added conditions; scales are a description | approved |
| 307 | D37 | the u-edge is the bump measurement, the per-slice edge the cross-check | approved |
| 308 | D38 | C2 factor from per-count QED weights (0.861 / 0.860, ± 0.004) | approved |
| 309 | D39 | e−e− step at u = 0 = pairs from the two π0 photons; CB kept | approved |
| 310 | D40 | low-mass tail of the η peak kept (t = 1), family t = 1 ± 0.5 | approved, option (a) |
| 311 | Dn, "Later" (Part I) | "Later" marks a finding of Part II or III made after Part I; Dn a numbered decision of Part III (I:1845; Is:705) | – |
| 312 | Amendment 1 | = D21 (III:6518). | "η yields corrected for the target-wall share of the η signal" |
| 313 | options (a)/(b) of D31, D40 | the analyst's chosen alternative (III:6597, 6612). Not variant (a), display (a), scenario (a). | – |
| 314 | status words | PRELIMINARY (σ(pp → ηX)); MAIN; CROSS-CHECK; NOT DETERMINED (B4 area, DR); PASS / FAIL / FAIL marginal (closure); simulation-limited (MC error > stat); extrapolated (bump background); UNCONFIRMED (Part I L_data); final / check / frozen / model / later / open (Part I status column, I:512–514) | in short reports write in words: "preliminary", "cross-check only", "not determined" |
| 315 | fall-backs (1)–(6) (IIIs:1251–1256; III:6823–6850) | (1) k_E1 outside 0.8–3 → E1 dropped; (2) A1/A2 fail → CR-K dropped, bump extrapolated, CR-L term in S10; (3) estimators fail closure → template fits cross-checks; (4) stage-B χ²/ndf > 1.5; (5) MC error > stat → simulation-limited; (6) σ PRELIMINARY | – |

---------------------------------------------------------------------------------------------------------------------------------

## 11. Work packages, packages and verifiers (long reports only; I1: never in short reports)

| # | label | meaning (SPEC:91–267; III:1170–1173, 6889–6960) |
|---|---|---|
| 316 | WP0 | histogram factory and frozen definitions (defs.h, hists_v1.root) |
| 317 | WP1 | "before": absolute per-channel comparison, no renormalisation; Table A2 of corrections |
| 318 | WP2 | normalisation items: sector-2 deficit (C_sec), empty target (f_win, f_η) |
| 319 | WP3 | control regions, stage A |
| 320 | WP4 | stress closure on pseudo-data (configurations BA … RB) |
| 321 | WP5 | stage B template fits (F-main, F-free, F-gold) |
| 322 | WP6 | side-band estimators (B4, peak side bands) and cross-checks |
| 323 | WP7 | systematics (rules MC(Δ), R1) |
| 324 | WP8 | physics quantities (efficiencies, QED check, (p_T, y) reweighting, family S12) |
| 325 | WP9 | X17 (steps X-a … X-d); WP10 = reports (SPEC:267) |
| 326 | packages | e1 (typical cuts), bedge (per-slice bump edge), flex (flexible SMASH), edgeu (u-edge), eestep (like-sign step), ecalacc (ECAL acceptance), peaktail (peak tail), c2recon (per-count C2 weights), x17ab / x17c / x17d (X17 steps), vtx (vertex study), udatasim (data vs SMASH in u), fix / fix2 / fix3 / p3w (review and integration redraws) |
| 327 | verifiers | v0 (WP0), v5 (WP5), v8 (WP8), vx (X17), e1v, vbedge, vflex, vedgeu, vc2recon, veestep, vpt, vecalacc, ve_verify (vertex study), vw2a / vw2b (empty target, halo) |
| 328 | first study, first-study table | the first data cut-scan study of 2026-11-11 and its table rw_compare_enriched.txt (6080.5 nb⁻¹, no 64/63) (II:1491–1493; IIIs:1392) |
| 329 | old (LGAD-based) production | earlier production with the LGAD START condition (II:3738–3788) |

---------------------------------------------------------------------------------------------------------------------------------

## 12. Specification figure and table labels (Part III; written "spec. Fig. X" / "spec. Table X")

| # | label | content (III line of the caption) |
|---|---|---|
| 330 | spec. Table A2 | corrections: what, procedure, applied or not, size (III:1362–1370) |
| 331 | spec. Table B1, Fig. B1–B4 | data/SMASH per M4 region before renormalisation (1440); M4 0.30–1.10 before (1453); same 0.40–0.70 (1467); ratio map M4 × M(π+π−) (1502); M(π+π−), M(e+e−) per region R0–R7 (1526) |
| 332 | spec. Tables C1–C3, Figs. C1–C3 | candidate control regions (2106); nested stage-A models, four configurations (2182); counted "after" display per region (2258); N − CB map (2223); per-slice ratio to B_A + 0.4η (2234); "after" spectrum (2246) |
| 333 | spec. Table D1 | observed 90 % CL_s limits at three masses, three normalisations (6136). Not decision D1. |
| 334 | spec. Tables E1–E3, Figs. E1–E4 | F-main parameters (2508); template-fit yields (2556); side-band estimators (2961); F-main fit (2529); data minus CB minus background (2543); profile scans (2674); side-band estimators (2977). Not classes E1, E2, cut sets, families. |
| 335 | spec. Table N1, Fig. N1 | normalisation items (3335, 3372) |
| 336 | spec. Fig. P1, Tables P1, P2 | summary of the physics quantities (3832); efficiencies (wp8_tableP1.tex); physics quantities (3728). Not pion classes P1, P2. |
| 337 | spec. Table S1, Fig. S1 | systematic uncertainties (3061, 3087). Not slice S1 or family S1. |
| 338 | spec. Table V1, Fig. V1 | stress closure (2382, 2396). Not weight variant V1, not PID file v1. |
| 339 | spec. Figs. X1–X4 | M(e+e−) of the η peak band (6188); limit against m_X (6202); local p_0 (6217); validation at 17 MeV/c² (6218). Not cut sets X1–X4, not Part I estimators X1–X5. |

---------------------------------------------------------------------------------------------------------------------------------

## 13. Other bookkeeping labels

| # | label | definition (source) | plain phrase |
|---|---|---|---|
| 340 | period I, IIA, IIB, III | beam-time days 034–046, 047–048, 049–053, 059–067 (days 054–058 absent) (defs:69–77; II:3450). Not Parts I–III. | "days 034–046" etc. |
| 341 | fold | parity of the event-group index; two independent halves of the data (II:3464). | "two independent halves of the data" |
| 342 | us, pp, mm | trees of the unlike-sign (π+π−e+e−) and like-sign (π+π−e+e+, π+π−e−e−) hypotheses (II:3456). | – |
| 343 | day file, file_id | one merged PAT output per day (lep_034 … lep_067) and its day number (II:3449). | – |
| 344 | no_hf run | run with L = 0 in the luminosity table (II:3463). | – |
| 345 | class A, class B (II) | counts with both leptons β_DST ≥ 0.9 / at least one below; r_B share of the η signal in class B (II:3491). | – |
| 346 | lepton class n (II) | number of leptons of an event, n(e+) + n(e−) = 2, 3, ≥ 4 (II:3517). | – |
| 347 | ECAL sector 0 | ECAL sector at azimuth 63–117°, off in data and simulation (II:3537). Not tracking sector 2 of the trigger deficit. | – |
| 348 | sector-2 deficit, PT2 check, C_sec, ε_s2 | hardware trigger deficit in sector 2 (φ = 180–240°) on days 034–044 (elastic/PT2), PT3 correction span 034–045 (D26); C_sec its correction of the η yields (II:3480–3481; III:6670). | "correction for a trigger loss in one sector on the first beam days" |
| 349 | target walls, f_win, f_η, empty target | windows and wall of the LH2 cell; f_win share of the walls in N − CB of a window (13–15 %); f_η share in the η signal, 3 (+4 −3) % (II:3475–3479; III:6671). | "target-wall share" |
| 350 | v1, v6 | PID cut files: data betanew_v1, simulation betanew_v6 (II:3488; I:1816). | – |
| 351 | figure-legend codes LS, COR, N++NRICH, FD | LS = like-sign (II figure "LS shape"); COR = corrected momenta (II cut-flow figure); FD = forced decay (I technical note) | write in words |

---------------------------------------------------------------------------------------------------------------------------------

## 13b. Further notation of the three long reports (added in version 2)

Rows of the notation tables of Part I (I:1779–1850), Part II (II:3448–3544) and Part III (III:6653–6757) that version 1 did not
list. Most are physics or bookkeeping symbols rather than catalogue labels; they are given so that a caption can spell them out.

| # | label | definition (source) | plain phrase |
|---|---|---|---|
| 352 | DST, PAT, FAT | reconstructed data; PostDST Analysis Tool (hypotheses, PID, relative times, ranking); Final Analysis Tool (selection, kinematics, ntuples) (II:3448; I:1781–1782) | "the reconstruction and the two analysis programs" |
| 353 | count | one stored combination π+π−ℓℓ of tracks of a hypothesis (one ntuple entry); with isBest ≥ 0 an event can give more than one count (II:3451; III:6654; I:1786) | "count (one stored π+π−e+e− combination)" |
| 354 | n_g, var/N | counts of one event group in a window; Σ_g n_g² over the counts: event-group factor of the variance (II:3453) | – |
| 355 | T_beam, √s; T (Part I) | beam kinetic energy 4530 MeV in all data kinematics (FATcfg beam.kinetic_energy; Scan:64); √s = 3.467 GeV (II:3454); Part I: T = beam kinetic energy, SMASH generated at 4500 MeV (I:1784; defs:66–67) | "beam kinetic energy 4.53 GeV (data analysis) / 4.50 GeV (SMASH generation)" |
| 356 | f_LF, LF2, LF3 | PT3/PT2 live-fraction ratio LF3/LF2 × 64/63 as defined in the luminosity files (0.932 overall; 0.9334 with the analysis' L weighting); not applied (D17) (II:859–865, 3470) | "trigger live-time fraction (not applied)" |
| 357 | R (Part II) | counts of "pre" per day × 64/63 / L, the rate per day (II:3473) | "count rate per day" |
| 358 | σ_int | elastic cross section integrated over the acceptance, used for L (2.90699 mb) (II:3474) | – |
| 359 | U/L, h (halo level) | front counts per nb⁻¹; h = 100 U/N, front counts per 100 counts of "pre" + vertex (II:3477) | "upstream (beam-halo) counts per 100 counts" |
| 360 | R_sort, s_lost, T_DS | sorter correction of the luminosity analysis, R_sort = 1 + s_lost loss/((1 − loss) T_DS), s_lost = 0.544, T_DS = 0.960–0.970; R_sort − 1 = 1.3 %, one-sided, not applied (II:3482, 3714–3717) | "possible event loss of the START-based track sorter (not applied)" |
| 361 | β_new, β_DST | START-free velocity from the time differences of the tracks of a hypothesis (II:3485, Eq. tnew; I:1814); DST velocity measured against the START time (II:3486) | "velocity without the START detector" / "velocity against the START time" |
| 362 | reference track | the track whose time offset is shared with the others (e− in π+π−e+e−) (II:3487) | – |
| 363 | ring, kIsLepton; lepton / hadron pass | RICH ring correlated with the track (lepton slot); track-sorter lepton flag; the two passes of the DST track sorter (DST β ≥ 0.9 / β ≥ 0 and tof ≤ 60 ns) (II:3489–3490; I:1817) | – |
| 364 | mass scale | fitted position / nominal mass: of the η peak in M4 (0.9927 ± 0.0009), of K0S in M(π+π−) (0.99691 ± 0.00036, used as 0.9970) (II:174, 3495; Scan:115–121) | "measured mass scale" |
| 365 | Δz_ππ, DOCA | z of the crossing point of the two pion straight lines minus vz; their distance of closest approach (K0S selection: Δz_ππ > 20 mm, DOCA < 10 mm) (II:3496; Scan:117–118; SimE:380–398) | – |
| 366 | p_T,4, y4 | transverse momentum (lab) and c.m. rapidity of π+π−e+e− (II:3498; I:1822); simulated y4 recomputed at 4500 MeV in Part III (III:1209) | – |
| 367 | bump, peak | η → π+π−π0 seen below m_η (photon lost); η at m_η (η → π+π−e+e− or η → π+π−γ with conversion) (II:3511; I:1829) | "the bump (η → π+π−π0 with a lost photon)" / "the η peak" |
| 368 | X (Part II) | expected excess of the plain CB formula applied to isBest = 1 counts (II:3518) | – |
| 369 | B, S/B (Part II) | pol2 + CB in the data peak window; S in the window over B (II:3525) | "background under the peak (polynomial + CB)" |
| 370 | F, F_g, F_tot (Part II) | S/stat (compares selections, not a significance); the same with the event-group term; S/√(stat² + model²) (II:3526) | "peak content over its error" |
| 371 | M4max, m_edge (Part II) | largest M4 of η → π+π−π0, π0 → e+e−γ (photon lost) at given M(π+π−) for M(e+e−) → 0, Eq. m4max; its measured position (II:1589, 3532; Scan:141–147) | "kinematic upper edge of the bump" |
| 372 | Δ_3%, ΔS_CB | change of N − CB, of S, for a 3 % change of the CB (II:3533; spec_c:65 kCBsyst = 0.03) | "change for a 3 % change of the CB" |
| 373 | f, f++, f−− (Part II) | pass fraction of the counts for a cut; "CB kept" = √(f++ f−−) (II:3535) | – |
| 374 | ε_S, ε_B, ε_T | kept share of the peak signal (pass/fail fits), of pol2 + CB in the window, of T (II:3538) | – |
| 375 | X17 | hypothetical 17 MeV/c² boson in η → π+π−X17, X17 → e+e− (II:3544) | "hypothetical X17 boson (17 MeV/c²)" |
| 376 | Part I, II, III | simulation report (report_sim.pdf of the simulation branch); data-only report (report_exp.pdf); data + SMASH report (report_all.pdf) (III:6653) | "the simulation report / the data report / the data–simulation report" |
| 377 | close pairs | e+e− pairs with a small opening angle and small M(e+e−) (conversions, π0 Dalitz pairs); D32 (III:6676) | "close e+e− pairs" |
| 378 | M4, M4', MM4, m5c | M(π+π−e+e−); its simulated, mass-transferred value; missing mass of π+π−e+e−; M(π+π−e+e−γ) of the ECAL tag with the photon rescaled (III:6696) | – |
| 379 | SMASH, PLUTO; converter | event generator of p+p collisions; decay generator used by the converter; convert-smash: SMASH output → HGeant input, forces the decays of the lepton and η sets (I:1779–1780) | – |
| 380 | N, b_max (Part I) | events per set (10^10); largest impact parameter (1.2 fm), L_sim = N/(π b_max²) = 221.05 nb⁻¹ (I:844–845, 1784) | – |
| 381 | geninfo, geninfo1, geninfo2 | generator codes of a GEANT track: decay code, code of the direct parent, counter of the decay in the event; PAT stores them as float (I:1793, 1868–1872) | – |
| 382 | forced decay | a decay chosen by the converter instead of the natural one, with a stored weight w (I:1794) | – |
| 383 | σ_prod, BR, ε, N_gen | SMASH production cross section (mean number per collision × π b_max²); branching ratio (PDG unless stated); counts after "pre" per generated decay; number of generated decays (I:1806–1808) | – |
| 384 | N_eff | effective number of simulated entries, (ΣW)²/Var (I:1810; I5:32–33) | "effective number of simulated entries" |
| 385 | charge × momentum | horizontal axis of the PID plots (negative for negative tracks) (I:1815) | – |
| 386 | named selections (Part I) | pre, pre + OA > 2, slides, old code, peak nominal, bump nominal; the engine also has "pre + pions in vertex fit" (I:1828; SimE:409–425) | – |
| 387 | template noise; stress test | statistical fluctuation of the simulated templates; Part III's test of its estimators on pseudo-data with realistic template noise (I:1846) | – |
| 388 | target walls, beam halo (Part I) | windows and wall of the hydrogen cell; beam particles off the axis that hit material around the target (neither is simulated) (I:1847) | – |
| 389 | ECAL | electromagnetic calorimeter (6 sectors × 163 modules); sector 0 at azimuth 63–117° switched off (I:1849) | – |
| 390 | lost photon on its circle | the photon of a simulated bump combination reconstructed from the four true charged momenta, m_π0 and m_η, up to one angle (I:1850); basis of the per-count C2 weight (row 144) | – |
| 391 | E2a, E2b; s_T3 (III flexible-fit variants) | the class E2 split into 10–20 and ≥ 20 MeV/c² of M(e+e−) (flex_lib.h:21, 57, 108–109); s_T3 = free scale of the fake-lepton template T3 (III:5738–5740) | "M(e+e−) 10–20 MeV/c² / 20 MeV/c² or more" |
| 392 | BR_1, BR_2 (III, X17 validation) | toy estimates of the X17 BR under normalisation (1), BR_1 = r̂ × 2.68 × 10⁻⁴, and (2), BR_2 = N_X/s = r̂ N̂_5/N_5^true × 2.68 × 10⁻⁴ (III:5990–5991) | – |

---------------------------------------------------------------------------------------------------------------------------------

## 14. Labels with several meanings (disambiguation)

| token | meanings (where) |
|---|---|
| S0 | selection S0 (III); sector 0 (II figures of the sector study, II PDF text "Days Tracks S0 S1 … S5", III Fig. N1(a)) |
| S1–S5 | M(π+π−) slices (III); systematic families S1–S13 (III); sectors S0–S5 (II, III Fig. N1); signal classes in a Part I figure legend ("S = S4 + S5", "S = S2 + S3": classes 4+5, 2+3; I PDF text) |
| E1, E2 | M(e+e−) classes (III); cut sets E1–E4 (III §e1, sans-serif); bump-edge families E1–E10 (III §bedge); spec. Figs./Tables E1–E4 (III); Part II's edge-step sum E, Part III's edge function E(m, q); E_π0 in Eq. (edge) |
| X1–X5 | Part I X17 estimators X1–X5 (I5); Part III cut sets X1–X4, X2old (III §e1); spec. Figs. X1–X4 (III §x17) |
| M0–M3 | Part III stage-A background models M0–M3; Part I extraction methods M1, M2, M2b, M2d, M3, M3c |
| B | bump window (III); Part II background pol2 + CB in the peak window (II:3525); Part I background counts B (I:1833); B_A stage-A background; B4 (estimator, two definitions); B1–B3 (Part I methods); spec. Table/Figs. B1–B4; half B; class B |
| B4 | Part I side band 0.5665–0.70 with the half-A shape; Part III bump estimator from U = 0.575–0.70 with the stage-A shape |
| P | peak window (III); purity (Part I); pion-pair classes P0–P9 (I); "P3" model of the flexible grid; spec. Fig. P1, Tables P1, P2 |
| T | ECAL-tagged content (II, III); beam energy T (I) / T_beam (II); templates T1–T4, T12; T4530 (production name); T_DS, T_exp |
| R | rate per day (II); B4 transfer ratio B_A(B)/B_A(U) (III); R′ side-band ratio; R_X X17 efficiency ratio; R_34; R0–R7 regions; rule R1; R_pk, R_SB (II); R_sort (II) |
| U | upper region 0.575–0.70 of B4 (III); Part II halo counts at −200 ≤ vz < −170 mm; families U1–U11 |
| u | distance to the bump edge (II, III); u_j = (M(e+e−) − 30 MeV)/30 MeV in the X17 model (III:5909) |
| t | peak-tail scale (III); t_i = (M4 − 0.55)/0.1 in the X17 model (III:5909) |
| r | X17 POI (III); r_res mass-transfer resolution factor (written r in package tables); Part II r = η per count in the walls over hydrogen (II:3479), r_B; r_SB |
| D | decisions D1–D40; Part I method D5; Part II side-band-subtracted counts D; X17 M(e+e−) shape D of C5 |
| λ | template-noise weight (flex); free factor on G_D + G_C (Fig. res-udatasim); ECAL cluster-energy scale (1.285); signal leakage of Part I M3 |
| κ | signal conversion scale (model); κ (bump edge) of C3; κ_B background conversion scale |
| (a) | ghost variant (a); SMASH display (a); closure scenario (a); option (a) of D31/D40; reading (a) of κ in §bedge |
| f | Part II pass fraction of a cut; Part III side-band leakage f; f_E2; f_η; f_win; f_s; f(cell given band) |
| h | step height (II, III); Part II halo level h = 100 U/N |
| b | common tilt (III); Part II flux ratio full/empty; b_max (I) |
| F | Part I F = S/√(S + B) (F_I); Part II F = S/stat; flex F7–F9; layouts F10, F5g; F-main/F-free/F-gold; F_V, F_VMD form factors |
| V | vertex class V (with vertex); weight variants V0–V4; spec. V1; PID files v1, v6; X17 selection v1 |
| A | half A (I); mode A (I, III); class A (II); tests A1–A5 (III); spec. Table A2; "A2 note", "A3 note" in defs.h comments |
| N | counts N; normalisation error N (III); N = 10^10 events per set (I) |
| σ_vis | Part I visible cross section (nb); Part III visible yield per L, explicitly "not a cross section" |
| CR-U | 0.565–0.700 (definition: defs, SPEC; flexible fits: flex code, III:757, 5460, IIIs:981) and 0.575–0.700 with M(e+e−) ≥ 10 MeV/c² (adopted for stage A and B4: wp3/wp6 code; IIIs:1397); the §e1 display normalises in 0.575–0.70 over all M(e+e−) |
| side bands | Part II 0.5225–0.53364 and 0.55416–0.56530; Part III 0.5250–0.5336 and 0.5542–0.5628; X17 M(e+e−) side bands 2–12 and 28–62 MeV/c²; Part I M3 side bands 0.5255–0.5333 and 0.5587–0.5665 (simulation scale, 7.8 MeV/c² wide); Part I B4 side band 0.5665–0.70 |
| X17 window | 0.5333–0.5587 (Part I, simulation scale); 0.531249–0.556551 (Part II moved); 0.5312–0.5566 (III X17c); P (X17P); 0.5325–0.5575 (III §e1, peak band); X2D (fit) |
| peak window | data 0.53364–0.55416; simulation 0.5333–0.5587 |
| bump window | data 0.4550–0.5225; simulation 0.4576–0.5255 |
| "vertex" | event vertex (vz > −500, target region vz > −160); "+ vertex" = both pions in the event-vertex fit |
| w | Part I: generator weight stored on a track (W = weight of a count); code and Part III: the weight of a count (= Part I's W) |
| E2a, E2b | flexible-fit variant classes M(e+e−) 10–20 and ≥ 20 MeV/c² (III:5740); not the cut sets E1–E4 and not the class E2 |

---------------------------------------------------------------------------------------------------------------------------------

## 15. Ready-made "Selection:" sentences

Building blocks used below (copy as written). **Rule for the SHORT reports (v2):** leave out every label in parentheses — (isBest ≥ 0),
(S0), (V), (Nv), (V3), (pre), and so on; the words in front of them are complete. The long reports may keep them.

- BASE = "PT3 trigger, an event vertex in the target region (z > −160 mm), particles identified by their START-free velocity
  β_new (pions within 3σ, leptons within 5σ; leptons with a RICH ring), energy-loss-corrected momenta for all four tracks, all stored
  π+π−e+e− combinations of each event (isBest ≥ 0)"
- PIVTX = "both pions used in the event-vertex fit"
- CB = "combinatorial background CB = 2√(N++ N−−) from the like-sign data (π+π−e+e+, π+π−e−e−) with the same cuts, in 10 MeV/c²
  cells of M(π+π−e+e−), k = 1" (Part III layouts add: "shared out over the M(π+π−) and M(e+e−) cells by the like-sign share of each
  50 MeV/c² band")
- SMASH3 = "SMASH simulation of p+p at 4.5 GeV (standard, lepton and η sets, 221.05 nb⁻¹ each), every channel at its SMASH cross
  section × PDG branching ratio × η-decay physics weights, simulated M(π+π−e+e−) × 0.99615 and M(π+π−) × 0.9952 to the data scale,
  with the same cuts"
- SMASH1 = "SMASH simulation only (no combinatorial background; η-decay physics weights of the main model: virtual-photon
  spectrum with form factors and decay angles for η → π+π−e+e−, measured Dalitz plot for η → π+π−π0, pion form factor for
  η → π+π−γ), counts expected at 6080.5 nb⁻¹" (long report: add "(V3)")

### 15.1 Base selections (all M(π+π−), all M(e+e−))

| # | label | sentence |
|---|---|---|
| Q1 | S0 (III) | "Selection: data, BASE, PIVTX; all M(π+π−) and M(e+e−); L = 6079.2 nb⁻¹." |
| Q2 | "pre" (I, II, III) | "Selection: data, BASE; no further cut." |
| Q3 | "pre + vertex" (II) | same as Q1 (write "pre + vertex" in parentheses only in the long Part II) |
| Q4 | PEAK (II) | "Selection: data, BASE, both pion tracks within 40 mm of the event vertex along the beam, and missing mass of π+π−e+e− below 2.85 GeV/c²." |
| Q5 | DIRECT (II) | "Selection: data, BASE, both pion tracks within 40 mm of the event vertex along the beam." |
| Q6 | BUMP (II) | "Selection: data, BASE (no further cut)." |
| Q7 | X17 band (II) | "Selection: data, BASE, PIVTX, 15.5 ≤ M(e+e−) < 24.5 MeV/c²." Add the window: "… in M(π+π−e+e−) 0.5312–0.5566 GeV/c² (the simulation window moved to the data mass scale)" |
| Q8 | peak nominal (+ vertex) (II) | "Selection: data, BASE, angle between the π+π− and e+e− momenta in the lab below 55°, M(π+π−) > 0.30 GeV/c² [, PIVTX]." |
| Q9 | bump nominal (+ vertex) (II) | "Selection: data, BASE, angle between the π+π− and e+e− momenta in the lab below 40°, 0.30 < M(π+π−) < 0.41 GeV/c² [, PIVTX]." |
| Q10 | isBest = 1 variant | replace "all stored π+π−e+e− combinations of each event (isBest ≥ 0)" by "only the best-ranked combination of each event (isBest = 1)" |
| Q11 | OA > 2° variant | append "and an e+e− opening angle above 2°" |
| Q12 | vz150 / vz170 / none | replace "(z > −160 mm)" by "(z > −150 mm)" / "(z > −170 mm)" / drop the target-region cut ("any reconstructed event vertex") |
| Q13 | vertex study V / Nv | "Selection: data, PT3 trigger, β_new PID, corrected momenta, all combinations; counts with an event vertex (V) and without one (Nv); no target-region cut." |
| Q14 | Part III simulation under S0 | "Simulation: SMASH3, Q1 cuts (the simulated PT3 trigger is the trigger bit 8192)." |
| Q15 | Part I simulation, "pre" | "Simulation: SMASH1; BASE (simulated trigger bit 8192); windows on the simulation mass scale (peak 0.5333–0.5587, bump 0.4576–0.5255 GeV/c²)." |
| Q16 | Part I frozen PEAK / DIRECT / BUMP / X17 | append to Q15: PEAK "both pion tracks within 40 mm of the vertex and missing mass below 2.85 GeV/c²"; DIRECT "both pion tracks within 40 mm of the vertex"; BUMP "no further cut"; X17 "PIVTX, 15.5 ≤ M(e+e−) < 24.5 MeV/c², M(π+π−e+e−) 0.5333–0.5587 GeV/c²" |

### 15.2 Ranges of M(π+π−) and classes of M(e+e−) (append to any base sentence)

| # | condition | words to append |
|---|---|---|
| Q17 | all | "all M(π+π−) and M(e+e−)" |
| Q18 | slice S1 / S2 / S3 / S4 / S5 | "M(π+π−) below 0.34" / "0.34–0.37" / "0.37–0.39" / "0.39–0.43" / "above 0.43 GeV/c²" (simulation: on 0.9952 × M(π+π−)) |
| Q19 | the five slices | "in the five M(π+π−) slices below 0.34, 0.34–0.37, 0.37–0.39, 0.39–0.43 and above 0.43 GeV/c²" |
| Q20 | per-slice edge | "in the four M(π+π−) slices 0.30–0.34, 0.34–0.37, 0.37–0.39 and 0.39–0.40 GeV/c²" |
| Q21 | u-edge | "0.30 ≤ M(π+π−) < 0.40 GeV/c², M(e+e−) below 135 MeV/c² (the π0 mass; above it the distance to the edge is not defined)" |
| Q22 | E1 / E2 | "M(e+e−) below 10 MeV/c²" / "M(e+e−) of 10 MeV/c² or more" |
| Q23 | F5g | "M(e+e−) of 20 MeV/c² or more" |
| Q24 | X17 band / X2D | "15.5 ≤ M(e+e−) < 24.5 MeV/c²" / "M(π+π−e+e−) 0.5225–0.6100 GeV/c² and M(e+e−) 2–62 MeV/c²" |

### 15.3 The most frequent full sentences (written out)

| # | where | sentence |
|---|---|---|
| Q25 | III, F10 E2 per slice (stage A, F-main, flex) | "Selection: data, PT3 trigger, an event vertex in the target region (z > −160 mm), β_new particle identification (pions 3σ, leptons 5σ, leptons with a RICH ring), corrected momenta for all four tracks, both pions used in the event-vertex fit, all stored combinations per event (isBest ≥ 0); M(e+e−) of 10 MeV/c² or more; in the five M(π+π−) slices below 0.34, 0.34–0.37, 0.37–0.39, 0.39–0.43 and above 0.43 GeV/c²; M(π+π−e+e−) 0.40–0.70 GeV/c²." |
| Q26 | III, u-edge | "Selection: data, PT3 trigger, an event vertex in the target region (z > −160 mm), β_new particle identification, corrected momenta for all four tracks, both pions used in the event-vertex fit, all stored combinations per event; 0.30 ≤ M(π+π−) < 0.40 GeV/c², M(e+e−) below 135 MeV/c²; distance to the bump edge −45 to +45 MeV/c² in 2.5 MeV/c² bins; CB = 2√(N++N−−) from the like-sign data in 10 MeV/c² cells of the distance, separately in 0.30–0.34, 0.34–0.37, 0.37–0.39 and 0.39–0.40 GeV/c² of M(π+π−)." |
| Q27 | III, peak side bands | "Selection: as Q1; peak window M(π+π−e+e−) 0.53364–0.55416 GeV/c², side bands 0.5250–0.5336 and 0.5542–0.5628 GeV/c², all M(e+e−)." |
| Q28 | III, B4 | "Selection: as Q1 with M(e+e−) of 10 MeV/c² or more; bump window 0.4550–0.5225 GeV/c², background from 0.575–0.70 GeV/c² carried down with the simulated background shape; extrapolated to all M(e+e−) with the simulated share of the bump at M(e+e−) ≥ 10 MeV/c² (0.4773)." |
| Q29 | III, X17 fit | "Selection: as Q1; M(π+π−e+e−) 0.5225–0.6100 GeV/c² and M(e+e−) 2–62 MeV/c² (M(e+e−) 12–28 MeV/c² hidden until the procedure was frozen)." |
| Q30 | III, η peak band | "Selection: as Q1; η peak band M(π+π−e+e−) 0.5325–0.5575 GeV/c²." |
| Q31 | III, typical cuts E1 on "pre" | "Selection: data, BASE; plus the cuts of earlier slides: angle between π+π− and e+e− in the lab below 50°, the same angle in the η rest frame above 140°, M(π+π−) below 0.46 GeV/c², e+e− opening angle above 2°." (S0 + E1: add PIVTX) |
| Q32 | II, per-slice edge | "Selection: data, BASE (= frozen bump selection), all M(e+e−); in the four M(π+π−) slices 0.30–0.34, 0.34–0.37, 0.37–0.39 and 0.39–0.40 GeV/c²; CB." |
| Q33 | II, u | "Selection: data, BASE, PIVTX; 0.30 ≤ M(π+π−) < 0.40 GeV/c², M(e+e−) below 135 MeV/c²; CB in 10 MeV/c² cells of the distance to the edge, per M(π+π−) slice." |
| Q34 | II, η peak, PEAK | "Selection: data, BASE, both pion tracks within 40 mm of the event vertex along the beam, missing mass of π+π−e+e− below 2.85 GeV/c²; peak window 0.53364–0.55416 GeV/c²; fit with a fixed Gaussian (0.5439 GeV/c², 5.4 MeV/c²) plus a second-order polynomial in 0.5225–0.6100 GeV/c²." |
| Q35 | II, ECAL tag | "Selection: data, BASE; an ECAL photon (identified as a photon, 0.8 < β < 1.2, above 100 MeV; of several, the one with M(e+e−γ) closest to 0.140 GeV/c²) with 0.119 ≤ M(e+e−γ) < 0.173 GeV/c² and M(e+e−) < 0.10 GeV/c²; M(π+π−e+e−γ) with the photon rescaled to the π0 mass, fitted in 0.49–0.62 GeV/c² (2.5 MeV/c² bins, CB in 20 MeV/c² bins)." |
| Q36 | I, frozen peak selection | "Simulation only (SMASH, main η-decay physics weights; long report: V3); BASE; both pion tracks within 40 mm of the event vertex, missing mass below 2.85 GeV/c²; peak window 0.5333–0.5587 GeV/c² (simulation scale)." |

---------------------------------------------------------------------------------------------------------------------------------

## 16. Disagreements: status after the code check (version 2)

Each item has a status:
- **SETTLED (code)**: the code decides what was done; the glossary rows say it, nothing has to change in a report.
- **DOCUMENTED**: the reports already say it correctly; captions must still give the range or value.
- **OPEN**: a report says something the code does not do, or uses one name for two things. The fix is a wording fix in the report
  named (no number changes); short reports avoid the labels anyway (I2).

Ordered as in version 1; items 21–23 are new.

1. **Peak side bands: two estimators, two ranges. SETTLED (code), by design.** Part II's side-band method (res_c:76–78) uses
   0.5225–0.53364 and 0.55416–0.56530 GeV/c² (11.14 MeV/c² each, touching P, scaled with the P_main pol2 ratio r_SB; II:1545,
   1627–1637, 3510). Part III's peak estimator (defs:242–243; wp6_ana.C) uses 0.5250–0.5336 and 0.5542–0.5628 GeV/c² (8.6 MeV/c²
   each, scaled with the stage-A ratio R′ and corrected for the leakage f; III:1263, 2830–2846, 6694; IIIs:1395). Both are right for
   their method. Captions must give the ranges; "side bands" alone is ambiguous. Detail: the rounded Part III edges leave
   0.04 MeV/c² slivers next to P (row 70).
2. **CR-U: one name, three regions in the code. SETTLED (code); OPEN (wording).** 0.565–0.700 all M(e+e−) is the definition
   (defs:228, 245, 274; SPEC:74) and is what the flexible fits use (flex_fit.C:144; flex_info.C:28; III:757, 5460; IIIs:981).
   0.575–0.700 with M(e+e−) ≥ 10 MeV/c² is what stage A and B4 use after the fall-backs (wp3_stageA.C:901–902, 1401; wp6_ana.C:14–17,
   941; III:2046–2051, 2088–2092, 1269, 6703). The §e1 display (b) uses 0.575–0.70 with all M(e+e−) (item 22). Open: the IIIs notation
   (IIIs:1397) gives only 0.575–0.700 while IIIs:981 calls 0.565–0.70 "CR-U". Recommendation: write the range in every caption.
3. **B4: one idea, two implementations. SETTLED (code); OPEN (wording in Part I).** Part I (b4I:64; I5:940): side band 0.5665–0.70
   GeV/c² on the simulation scale (≈ 0.564–0.697 on the data scale), background shape from simulation half A, all M(e+e−), on
   simulated pseudo-data. Part III (wp6_ana.C; III:2765–2781): U = 0.575–0.700 GeV/c² on the data scale, M(e+e−) ≥ 10 MeV/c², the
   stage-A shape ratio B_A(B)/B_A(U) = 0.12735, η_U subtracted (D27), extrapolated with f_E2 = 0.4773. Part I and Is say "later: B4 is
   not determined in the data" (I:232, 571; Is:162, 588); this is true of the Part III implementation and should say "the data
   version of B4 (Part III)".
4. **M1–M3: Part I methods and Part III stage-A models. OPEN (label clash).** Part I: M1 Gaussian + pol2, M2 template fit, M3 linear
   side bands (I5:930–939; b4I enum mM1…mD5). Part III: M0 scale only, M1 + common tilt (adopted), M2 + quadratic, M3 per-slice
   tilts (III:2139–2153, 2744). The "symbols with two meanings" paragraph of III (III:6636–6648) does not list it.
5. **X1–X5: three meanings. OPEN (label clash).** Part I X17 estimators X1–X5 (I5:941–944); Part III cut sets X1–X4, X2old
   (e1f:29–30; III:1693–1696, sans-serif); spec. Figs. X1–X4 (III:6188–6218). III separates the last two, not Part I's estimators.
6. **D5: a Part I method and decision D5. OPEN (label clash).** Part I: direct-decay yield from the M(e+e−) split at 10 MeV/c² with a
   2 × 2 unfolding (I5:941–942; b4I:65 kMeeSplit); Part III: D5 = the CB rule (III:6555).
7. **P3: pion class and flexible-fit model. OPEN (label clash).** P3 = both pions from one ω (truth:242, 252; I:1227); "P3" in the
   flexible ladder = Part III F-main reproduced (III:5492).
8. **V1: four meanings. OPEN (label clash).** Weight variant V1 = V0 × abs(F_V)² (phys:64); spec. Table/Fig. V1 = stress closure
   (III:2382, 2396; Part I: "its Fig. V1", I:266, Is:220); PID cut file v1 of the data (II:3488); X17 selection v1 (frz_x17v1; I5:907).
9. **X17 window in M(π+π−e+e−): five windows. SETTLED (code); OPEN (one wording).** All five are in the code: Part I frozen
   0.5333–0.5587 on the simulation scale (frz_x17:8–9; the blind test used v1, 0.5330–0.5555 with M(e+e−) 13.0–26.5, frz_x17v1:6–9,
   s5prep:103); Part II moved window 0.531249–0.556551 (opt_c:174–176; 6395 counts); Part III X17c 0.5312–0.5566 (defs:244; 6424
   counts), wider than both exact transfers (0.5333–0.5587 × 0.99615 = 0.531247–0.556549), so the text "× 0.99615" of defs:227 and
   III:1264 is exact only to 0.05 MeV/c² at the upper edge; X17P = P (defs:260; 5258 counts); the §e1 "X17 window" = the η peak band
   PB 0.5325–0.5575 × M(e+e−) 15.5–24.5 (e1f:146, 153, 155, 182–184; III:1647, 1671). Open: in §e1 call it "peak band × X17 M(e+e−)
   band", and the windows table (III:1255–1272) could list PB.
10. **σ_vis: "visible cross section" (Part I) and "not a cross section" (Part III). OPEN (naming).** Part I: ΣW/L_sim of counts =
    σ_prod × BR × ε, in nb (I:1809; Is:670). Part III: Y/L in visible counts per nb⁻¹, "not a cross section, although written with σ"
    (III:6739; IIIs:1408). Both are counts (combinations, more than one per event possible) per luminosity without efficiency
    correction, so the quantity is of the same kind; only the name differs. Part II writes Y. Recommendation: "visible yield" in all
    parts.
11. **Multi-pion. SETTLED (code).** spec5:27 = G_π0D + G_conv + G_vec, as I:1801 and I5:74. Is (Is:424, 675) omits the vector
    mesons, which are 0.107 of 659.9 nb (I5:88–90): a harmless simplification.
12. **G_D. SETTLED (code); OPEN (wording III:1109).** defs:306–307, 352–355: G_D = G_π0D + G_vec + G_oth. The code (BkgGroup,
    TemplateClass, brX17 = 0) shows G_oth = the leptons of a forced η → π+π−e+e− with pions not both from that η, i.e. of η origin
    (row 124). III:1109 "non-η Dalitz and vector-meson pairs" omits G_oth and is not strictly true for it; IIIs:1389 ("non-η background
    without a conversion pair") includes it but keeps "non-η". The phrase of row 125 avoids "non-η".
13. **Number of families. SETTLED (specification); OPEN (one word in IIIs).** The specification has S1–S12 (SPEC:280–293); S13 (peak
    tail) was added by D40 (tables_all/p3w_tableS1_v3.tex; III:3075–3078, 6638; IIIs:893, 912). The IIIs notation (IIIs:1417) still
    says "S1–S12".
14. **Edge scale factor. SETTLED (code); DOCUMENTED in Part II; OPEN (minor) in Part III.** 0.99277 = 0.5439/m_η (rounded frozen
    peak position; Scan:120 kScale) for the per-slice edge (wp6_fill.C:13; III:531, 2906, 2997, 3891; II prints 0.9927); 0.99274 ≈
    0.54388/m_η (free-fit position; defs:100 kEdgeScale; res_peak.C:401) for E(m, q) and u (III:1278, 4343; II:2859–2864). Part II
    states both values and their difference of 4 × 10⁻⁵ (II:2800–2801), so version 1's "no report states it" was wrong; Part III uses
    both without a word. Effect about 0.015 MeV/c² at the edge (notes_all/edgeu_results.md:119–120).
15. **Pion-vertex flag: "!= 0" against "== 1". SETTLED (data).** In the data compact file every pip/pim/ep/em_isUsedVertex value is
    0 or 1 (4 783 467 us, 122 875 pp, 371 524 mm entries; ROOT check verify_v2/vtxflag_check.C), and both tests select 3 640 366 /
    51 544 / 156 154 entries. spec_c:13, 79 and the other codes give the same selection. No action.
16. **Simulated missing mass of PEAK: 4530 MeV (Part I) and 4500 MeV (Part III). DOCUMENTED.** I:2055–2057; defs:66–67, 118;
    III:1209–1210; D4. A Part III caption of PEAK should say "simulated missing mass at 4500 MeV".
17. **"CB" names several constructions. DOCUMENTED.** Part II: 10 MeV/c² bins from 0.35 GeV/c², interpolated (Scan:77–82); Part III:
    the same cell edges from 0.30 plus the 50 MeV/c²-band spread (defs:283–297); ECAL tag: 20 MeV/c² bins in m5c (res_c:383–384);
    u-edge: 10 MeV/c² cells of u per M(π+π−) slice (edgeu_core.py:264; III:4346–4347); vertex study: smooth like-sign fits (D33).
    Window CB of X17P: 112 (geometric mean of the window), 171 ± 6 (spread rule), 102 (Part II rule) (III:6581). Captions must say which.
18. **S0–S5 as sector labels. DOCUMENTED.** The figure shared by Part II (II:1129) and Part III (Fig. N1(a); III:6637) labels the
    sectors S0–S5 and says so. Short reports do not print the selection label S0 anyway (I2).
19. **S2–S5 in Part I figure legends. OPEN (figure).** Fig. 25 of report_sim.pdf (Q = ε_S/ε_B of eight cuts) has the legends
    "S = S4 + S5, B = T2", "S = S2 + S3" and "B = T2 + S4 + S5" (classes 4 + 5, 2 + 3), which clash with slices and families; only a
    redraw with "classes 4 + 5" etc. removes it.
20. **X17 side bands in M(e+e−). DOCUMENTED deviation.** Specification 2–13 and 27–62 MeV/c² (SPEC:246); used 2–12 and 28–62 because
    the 2 MeV/c² bins of X2D that overlap 13–27 are masked (III:5895–5897; D30). No action.
21. **(new) w and W. OPEN (notation).** Part I: w = generator weight stored on a track by the converter, W = physical weight of a
    count (I:829–832, 1007–1028, 1794–1795). The code (SimE row fields w = sim_w, w_db) and Part III (III:1199, 6682, D7 "Σw/221.05",
    "w × V3" in many captions) use w for the count weight, i.e. Part I's W. Readers of Part III who look up w in Part I find another
    quantity. A sentence in the Part III notation ("w = W of Part I") would close it.
22. **(new) Normalisation region of the §e1 display. OPEN (wording III:1639).** III:1639 says c_b comes from 0.575–0.70 GeV/c² "(the
    control region of stage A)"; the code normalises over all M(e+e−) and all M(π+π−) in that M4 range (e1f:140 window "CR";
    e1_ana.py:218), while stage A uses M(e+e−) ≥ 10 MeV/c² only. Only the M4 range is shared.
23. **(new) M(e+e−) range of u. OPEN (wording Part III and IIIs).** The u-edge code drops M(e+e−) ≥ m_π0 (edgeu_core.py:146, 166);
    Part II says u is defined for M(e+e−) < m_π0 (II:2865); Part III and IIIs do not say it, and the u figure captions write "all
    M(e+e−)" (e.g. III:4321–4322, where it describes the window yields). Captions of u spectra and fits should say "M(e+e−) below
    135 MeV/c²".

Numbers were not changed. One inconsistency of numbers, for the numbers editor: the D40 term of the per-slice edge is ±3.6 % in
IIIs:1322 and III:1077 ("±3.6 % (3.7 % with the per-count C2 weights)"), ±3.7 % in IIIs:1367, and E10 = 3.4 / 3.7 / 4.9 % (without
the C2 correction / with it / with κ = 0.48) in tab:bedge-syst (III:3955; tables_all/bedge_syst_v3.tex). The 3.6 % is the value of
the peaktail package (notes_all/peaktail_results.md:8, 505, 677: slope +7.2 % per unit t × 0.5), the final combination
(bedge_syst_v3) gives 3.4 / 3.7 / 4.9 % for its three readings. Which reading "per slice ±3.6 %" refers to, and whether it should
read 3.4 or 3.7, has to be settled by whoever owns the numbers (not checked further here).

---------------------------------------------------------------------------------------------------------------------------------

## 17. Inventory: where each label is printed (rendered text of the six built PDFs)

Counts are regex matches in pdftotext -layout of report_sim.pdf, report_sim_short.pdf, report_exp.pdf, report_exp_short.pdf,
report_all.pdf, report_all_short.pdf (body, tables, captions and the text layer of the vector figures). Single-letter labels (B, P, U,
J, T, E) are not counted. A label counted under one token may have several meanings (§14); "." = not printed. Script:
/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ6/glossary_work/inventory.py.

Labels still printed in the SHORT reports (non-zero in the "short" columns) are the ones the editors must replace (I2).

Version 2 check: the counts of S0, PEAK, CR-U, F-main, B4, D40 and isBest were recounted independently with grep on the six text
extracts (all 42 numbers equal), and the extract of report_all_short.pdf is byte-identical to a fresh pdftotext of the PDF. The
reports were not changed after the extraction (newest build 15 November 23:46).

| label | I | I short | II | II short | III | III short |
|---|---:|---:|---:|---:|---:|---:|
| "pre" | 116 | 32 | 551 | 79 | 162 | 28 |
| S0 | 5 | 1 | 15 | . | 488 | 71 |
| S0_best / pre_best | . | . | . | . | 3 | . |
| S0_corfb / pre_corfb / corfb | . | . | . | . | 3 | . |
| dz40 | . | . | . | . | 4 | . |
| vz150 / vz170 | . | . | . | . | 6 | . |
| PEAK | 39 | 3 | 115 | 24 | 9 | 1 |
| DIRECT | 25 | 2 | 79 | 11 | 3 | . |
| BUMP | 24 | 2 | 36 | 6 | . | . |
| X17 band / X17 region | 1 | . | 62 | 15 | 4 | 1 |
| "pre + vertex" | . | . | 63 | 11 | 3 | . |
| peak nominal | 4 | . | 135 | 13 | 1 | . |
| bump nominal | 4 | . | 120 | 11 | . | . |
| slides (named selection) | 52 | 20 | 9 | . | 12 | 7 |
| old code | 4 | . | 10 | . | 8 | 4 |
| E3oa / E4oa | . | . | . | . | 27 | 2 |
| X2old | . | . | . | . | 24 | 3 |
| Nv (no event vertex) | . | . | 28 | 2 | 48 | 26 |
| half A / half B | 39 | 4 | 10 | 3 | 5 | . |
| period IIA / IIB | . | . | 77 | 2 | 48 | 2 |
| fold(s) | . | . | 24 | 1 | 25 | 3 |
| isBest | 23 | 6 | 360 | 59 | 147 | 29 |
| C1 | . | . | . | . | 2 | . |
| C2 | 25 | 14 | 8 | 2 | 466 | 95 |
| C3 | 14 | 8 | . | . | 214 | 31 |
| C4 | 10 | 9 | 1 | . | 207 | 34 |
| C5 | 14 | 9 | 1 | . | 251 | 39 |
| C6 | 3 | . | . | . | 44 | 1 |
| C2c | . | . | . | . | 1 | . |
| E23 | . | . | . | . | 2 | . |
| Gsig2-Gsig5 (figure legend) | 8 | 8 | . | . | . | . |
| class(es) 0-8 (Part I numbering) | 132 | 12 | 11 | 1 | 3 | 1 |
| L0-L11 | 30 | . | . | . | . | . |
| P0-P9 (also flex P3, spec. P1/P2) | 22 | 2 | . | . | 8 | 2 |
| T1 | 51 | 4 | . | . | 56 | 12 |
| T2 | 77 | 6 | . | . | 56 | 12 |
| T3 | 32 | 6 | . | . | 11 | . |
| T4 | 28 | 5 | . | . | 5 | . |
| T12 | . | . | . | . | 37 | 2 |
| G_eta | 13 | 3 | . | . | 194 | 29 |
| G_pi0D | 12 | 2 | . | . | 14 | 2 |
| G_conv | 15 | 4 | . | . | 22 | 2 |
| G_vec | 14 | 4 | . | . | 17 | 2 |
| G_oth | 8 | 3 | . | . | 17 | 2 |
| G_D | . | . | . | . | 171 | 21 |
| G_C | . | . | . | . | 164 | 23 |
| G_fake / G_diff | 9 | 2 | . | . | . | . |
| multi-pion | 14 | 6 | . | . | 8 | . |
| ghost | 16 | 3 | . | . | 33 | 9 |
| variant (a) | 1 | . | . | 1 | 19 | 3 |
| mode A / mode B | 12 | 1 | . | . | 22 | . |
| V0-V4 | 65 | 6 | . | . | 76 | 13 |
| rules v1 / v2 | 12 | 2 | . | . | . | . |
| pre-rw | 5 | . | . | . | . | . |
| GS (Gounaris-Sakurai) | . | . | . | . | 4 | . |
| w_db | . | . | . | . | 7 | . |
| v1 / v6 (PID cut files; also X17 selection v1) | 39 | 4 | 6 | . | . | . |
| s, r_res (mass transfer) | 3 | . | 1 | . | 48 | 9 |
| SB (side bands) | . | . | 3 | . | 11 | 1 |
| SBlo / SBhi | . | . | . | . | . | . |
| X17c | . | . | . | . | 10 | . |
| X17P | . | . | . | . | 12 | . |
| R0 | . | . | . | . | 20 | . |
| R1 | . | . | . | . | 44 | 4 |
| R2 | . | . | . | . | 14 | . |
| R3 | . | . | . | . | 10 | . |
| R4 | . | . | . | . | 14 | . |
| R5 | . | . | . | . | 11 | . |
| R6 | . | . | . | . | 8 | . |
| R7 | . | . | . | . | 18 | . |
| CR-U | . | . | 1 | . | 104 | 30 |
| CR-K | . | . | . | . | 74 | 16 |
| CR-L | . | . | . | . | 33 | 9 |
| F10 | 1 | . | . | . | 20 | . |
| F5g | . | . | . | . | 7 | . |
| X2D | . | . | . | . | 25 | 2 |
| G3 | . | . | . | . | 9 | . |
| BEFORE | . | . | . | . | 1 | . |
| peak band | 4 | . | . | . | 44 | 7 |
| R_pk / R_SB | . | . | 21 | . | . | . |
| S1E1 ... S5E2 (cell) | . | . | . | . | 4 | . |
| S1 | 1 | . | 4 | . | 46 | 7 |
| S2 | 2 | . | 4 | . | 36 | 4 |
| S3 | 2 | . | 4 | . | 51 | 4 |
| S4 | 4 | . | 4 | . | 48 | 4 |
| S5 | 4 | . | 5 | . | 64 | 9 |
| S6 | . | . | . | . | 10 | 1 |
| S7 | . | . | . | . | 9 | 1 |
| S8 | . | . | . | . | 20 | 2 |
| S9 | . | . | . | . | 2 | 1 |
| S10 | . | . | . | . | 25 | 2 |
| S11 | . | . | . | . | 8 | 1 |
| S12 | 1 | . | . | . | 20 | 3 |
| S13 | . | . | . | . | 9 | 2 |
| E1 | 1 | . | . | . | 297 | 26 |
| E2 | 1 | . | . | . | 344 | 39 |
| E3 | . | . | . | . | 52 | 5 |
| E4 | . | . | . | . | 47 | 4 |
| E5 | . | . | . | . | 10 | . |
| E6 | . | . | . | . | 4 | . |
| E7 | . | . | . | . | 10 | . |
| E8 | . | . | . | . | 9 | . |
| E9 | . | . | . | . | 4 | . |
| E10 | . | . | . | . | 10 | . |
| U1 | . | . | . | . | 19 | 2 |
| U2 | . | . | . | . | 12 | . |
| U3 | . | . | . | . | 6 | . |
| U4 | . | . | . | . | 13 | . |
| U5 | . | . | . | . | 23 | 1 |
| U6 | . | . | . | . | 4 | . |
| U7 | . | . | . | . | 8 | . |
| U8 | . | . | . | . | 23 | 1 |
| U9 | . | . | . | . | 2 | . |
| U10 | . | . | . | . | 6 | 1 |
| U11 | . | . | . | . | 15 | 3 |
| U8' (U8 prime) | . | . | . | . | 1 | . |
| M0 | . | . | . | . | 3 | . |
| M1 | 17 | 5 | . | . | 14 | . |
| M2 | 26 | 9 | . | . | 5 | . |
| M2b | 12 | 5 | . | . | . | . |
| M2d | 16 | 5 | . | . | . | . |
| M3 | 14 | 3 | . | . | 3 | . |
| M3c | 10 | 3 | . | . | . | . |
| B1 | 13 | 5 | . | . | 3 | . |
| B2 | 1 | . | . | . | 1 | . |
| B3 | 15 | 3 | . | . | 1 | . |
| B4 | 26 | 8 | 1 | . | 138 | 37 |
| D5 | 9 | 3 | . | . | 8 | 3 |
| X1 | 14 | 2 | . | . | 28 | 4 |
| X2 | 6 | 2 | . | . | 35 | 3 |
| X3 | 7 | 2 | . | . | 19 | 2 |
| X4 | 7 | 2 | . | . | 21 | 3 |
| X5 | 11 | 3 | . | . | . | . |
| X-a | . | . | . | . | 1 | . |
| X-b | . | . | . | . | 3 | . |
| X-c | . | . | . | . | 6 | . |
| X-d | . | . | . | . | 6 | . |
| P_main | 1 | 1 | 76 | 12 | 35 | 13 |
| P_pois | . | . | 2 | . | . | . |
| F-main | 7 | 5 | . | . | 166 | 28 |
| F-free | 5 | 5 | . | . | 33 | 2 |
| F-gold | 3 | 2 | . | . | 44 | 5 |
| F7 | . | . | . | . | 9 | . |
| F7b | . | . | . | . | 1 | . |
| F7bj | . | . | . | . | 21 | 2 |
| F7cj | . | . | . | . | 1 | . |
| F8 | . | . | . | . | 66 | 12 |
| F8b | . | . | . | . | 2 | . |
| F8bj | . | . | . | . | 3 | . |
| F9 | . | . | . | . | 6 | . |
| F9b | . | . | . | . | 1 | . |
| F9bj | . | . | . | . | 4 | . |
| ASIS | . | . | . | . | 1 | . |
| stage A | . | . | . | . | 85 | 16 |
| stage B | . | . | . | . | 8 | . |
| Step 0 | . | . | . | . | 1 | . |
| BA | 2 | 1 | . | . | 67 | 14 |
| AB | . | . | . | . | 6 | . |
| AA | . | . | . | . | 23 | 2 |
| BB | 2 | 1 | . | . | 24 | 2 |
| RF | 2 | 2 | . | . | 65 | 8 |
| RB | . | . | . | . | 22 | 4 |
| BB-lite | . | . | . | . | 14 | . |
| GLS / hybrid GLS | . | . | . | . | 74 | 4 |
| normalisation (1)/(2)/(3) | . | . | . | . | 35 | 8 |
| direct-decay / conversion route | 1 | . | 1 | . | 53 | 12 |
| A1 | . | . | 1 | . | 16 | 1 |
| A2 | . | . | 1 | . | 16 | 2 |
| A3 | . | . | . | . | 4 | 1 |
| A4 | . | . | . | . | 4 | 1 |
| A5 | . | . | . | . | 3 | 1 |
| fall-back 1 / 2 | . | . | . | . | 6 | . |
| DR | . | . | . | . | 43 | 7 |
| DR_F | . | . | . | . | 10 | . |
| DR_edge | . | . | . | . | 15 | 4 |
| DR_u | . | . | . | . | 17 | 5 |
| a_P | 4 | 4 | . | . | 121 | 28 |
| a_B | . | . | . | . | 40 | 7 |
| a_eta | 2 | 1 | . | . | 194 | 24 |
| Y_peak | 2 | 2 | . | . | 93 | 16 |
| Y_bump | 3 | 3 | . | . | 85 | 13 |
| f_eta(edge) / f_eta(tag) / f_eta(u-edge) | . | . | 3 | . | 108 | 29 |
| kappa_B | . | . | . | . | 55 | 7 |
| rho_d | 1 | . | . | . | 42 | 4 |
| k_E1 / k_E2 | . | . | . | . | 85 | 11 |
| k_cross | . | . | 1 | . | 25 | 3 |
| k_mix | . | . | 17 | 6 | . | . |
| BR_G / BR_A / BR_C | 4 | 1 | . | . | 85 | 15 |
| R_34 | . | . | . | . | 10 | 1 |
| MC(Delta) | . | . | . | . | 15 | 1 |
| D11data | . | . | . | . | 6 | . |
| Amendment 1 | . | . | 1 | . | 20 | 3 |
| K_p | . | . | . | . | 11 | . |
| flags G, M (sculpting) | . | . | 4 | 1 | . | . |
| best-F / loose / pure rules | . | . | 10 | 2 | . | . |
| class A / class B (DST beta) | . | . | 10 | . | . | . |
| lepton class n | . | . | 2 | . | . | . |
| GM / AM | . | . | 2 | . | . | . |
| LS (like-sign, figure legend) | . | . | 6 | 1 | . | . |
| COR (figure label) | . | . | 1 | . | . | . |
| FD (forced decay) | 1 | . | . | . | . | . |
| first study | . | . | 16 | 2 | 7 | 1 |
| C_sec | 1 | 1 | 5 | . | 68 | 8 |
| f_win | . | . | 23 | 6 | 48 | 10 |
| f_eta | 1 | . | 24 | 9 | 236 | 71 |
| C_vtx / eps_vtx | . | . | 16 | 7 | 58 | 18 |
| T_exp (ECAL tag expectation) | . | . | . | . | 7 | . |
| D1-D40 (all) | 65 | 24 | 70 | 35 | 744 | 243 |
| D1 | . | . | . | . | 19 | 3 |
| D2 | . | . | . | . | 3 | 1 |
| D3 | . | . | . | . | 7 | 4 |
| D4 | . | . | . | . | 1 | . |
| D5 | 9 | 3 | . | . | 8 | 3 |
| D6 | 2 | 1 | . | . | 11 | 2 |
| D7 | . | . | . | . | 2 | 1 |
| D8 | . | . | . | . | 18 | 2 |
| D9 | . | . | . | . | 2 | 2 |
| D10 | . | . | . | . | 1 | 1 |
| D11 | . | . | . | . | 23 | 5 |
| D12 | . | . | . | . | 3 | 2 |
| D13 | . | . | . | . | 3 | . |
| D14 | . | . | . | . | 2 | 1 |
| D15 | . | . | . | . | 2 | 1 |
| D16 | . | . | . | . | 1 | 1 |
| D17 | . | . | . | . | 4 | 1 |
| D18 | . | . | . | . | 2 | 1 |
| D19 | . | . | 1 | . | 10 | 2 |
| D20 | . | . | . | . | 17 | 3 |
| D21 | 1 | . | 10 | 5 | 33 | 11 |
| D22 | . | . | . | . | 17 | 8 |
| D23 | . | . | . | . | 6 | 3 |
| D24 | . | . | . | . | 10 | 2 |
| D25 | 5 | 4 | 5 | 3 | 35 | 14 |
| D26 | . | . | 4 | 2 | 11 | 4 |
| D27 | . | . | . | . | 6 | 1 |
| D28 | 1 | . | 3 | 2 | 14 | 7 |
| D29 | 2 | . | 1 | . | 19 | 8 |
| D30 | . | . | . | . | 24 | 6 |
| D31 | 1 | . | . | . | 19 | 10 |
| D32 | 2 | 1 | 4 | 2 | 35 | 9 |
| D33 | . | . | 5 | 3 | 21 | 4 |
| D34 | 2 | 1 | . | . | 11 | 3 |
| D35 | 2 | 1 | . | . | 61 | 22 |
| D36 | . | . | . | . | 22 | 11 |
| D37 | 12 | 4 | 7 | 3 | 27 | 15 |
| D38 | 10 | 4 | 7 | 3 | 92 | 31 |
| D39 | 6 | 1 | 9 | 6 | 22 | 7 |
| D40 | 10 | 4 | 14 | 6 | 120 | 31 |
| WP0 | . | . | . | . | 16 | . |
| WP1 | . | . | . | . | 12 | . |
| WP2 | . | . | 13 | 2 | 23 | 3 |
| WP3 | . | . | . | . | 29 | 1 |
| WP4 | 1 | . | . | . | 39 | . |
| WP5 | . | . | . | . | 32 | 2 |
| WP6 | . | . | . | . | 26 | 1 |
| WP7 | . | . | . | . | 54 | 8 |
| WP8 | 1 | . | . | . | 50 | 2 |
| WP9 | . | . | . | . | 3 | . |
| spec. Fig./Table <label> | . | . | . | . | 44 | 2 |
| package edgeu | 3 | 1 | 4 | . | 23 | 4 |
| package eestep | 2 | . | 6 | 1 | 14 | 2 |
| package ecalacc | 3 | . | . | . | 8 | . |
| package peaktail | 3 | 1 | 1 | . | 10 | . |
| package c2recon | . | . | . | . | 13 | . |
| package bedge | . | . | 2 | . | 14 | . |
| package flex | . | . | . | . | 12 | . |
| package e1 | 1 | . | . | . | 12 | . |
| package e1v | . | . | . | . | 13 | 1 |
| package x17ab | . | . | . | . | 13 | . |
| package x17c | . | . | . | . | 11 | . |
| package x17d | . | . | . | . | 10 | . |
| package udatasim | . | . | . | . | 5 | . |
| package vx | . | . | . | . | 16 | . |
| package v0 | . | . | . | . | 7 | . |
| package v5 | . | . | . | . | 23 | . |
| package v8 | . | . | . | . | 36 | . |
| package vbedge | . | . | 1 | . | 13 | . |
| package vflex | . | . | . | . | 16 | . |
| package vedgeu | . | . | . | . | 12 | . |
| package vc2recon | . | . | . | . | 10 | . |
| package veestep | . | . | 2 | . | 13 | . |
| package vpt | . | . | 1 | . | 12 | . |
| package ve_verify | . | . | 3 | . | 9 | . |
| package vw2a | . | . | 4 | 1 | 2 | 1 |
| package vw2b | . | . | 8 | . | . | . |
| package vecalacc | . | . | . | . | 12 | . |
| package fix | 3 | 1 | 4 | 3 | 9 | 1 |
| package fix2 | . | . | . | . | 1 | . |

288 catalogue entries (regexes in the script).

### 17.1 Labels printed in the short reports (counts per short report)

(Not found in any of the six reports: "B0"; not printed as labels: SBlo, SBhi, which exist only in defs.h.)

The list includes standard notation that a short report may keep when it is named in words at first use (isBest, "pre" as a quoted
word in running text is a label and must go; C_sec, f_win, f_eta are physics symbols). Everything else here is an internal label
under I2.

- "pre": I short 32, II short 79, III short 28
- S0: I short 1, III short 71
- PEAK: I short 3, II short 24, III short 1
- DIRECT: I short 2, II short 11
- BUMP: I short 2, II short 6
- X17 band / X17 region: II short 15, III short 1
- "pre + vertex": II short 11
- peak nominal: II short 13
- bump nominal: II short 11
- slides (named selection): I short 20, III short 7
- old code: III short 4
- E3oa / E4oa: III short 2
- X2old: III short 3
- Nv (no event vertex): II short 2, III short 26
- half A / half B: I short 4, II short 3
- period IIA / IIB: II short 2, III short 2
- fold(s): II short 1, III short 3
- isBest: I short 6, II short 59, III short 29
- C2: I short 14, II short 2, III short 95
- C3: I short 8, III short 31
- C4: I short 9, III short 34
- C5: I short 9, III short 39
- C6: III short 1
- Gsig2-Gsig5 (figure legend): I short 8
- class(es) 0-8 (Part I numbering): I short 12, II short 1, III short 1
- P0-P9 (also flex P3, spec. P1/P2): I short 2, III short 2
- T1: I short 4, III short 12
- T2: I short 6, III short 12
- T3: I short 6
- T4: I short 5
- T12: III short 2
- G_eta: I short 3, III short 29
- G_pi0D: I short 2, III short 2
- G_conv: I short 4, III short 2
- G_vec: I short 4, III short 2
- G_oth: I short 3, III short 2
- G_D: III short 21
- G_C: III short 23
- G_fake / G_diff: I short 2
- multi-pion: I short 6
- ghost: I short 3, III short 9
- variant (a): II short 1, III short 3
- mode A / mode B: I short 1
- V0-V4: I short 6, III short 13
- rules v1 / v2: I short 2
- v1 / v6 (PID cut files; also X17 selection v1): I short 4
- s, r_res (mass transfer): III short 9
- SB (side bands): III short 1
- R1: III short 4
- CR-U: III short 30
- CR-K: III short 16
- CR-L: III short 9
- X2D: III short 2
- peak band: III short 7
- S1: III short 7
- S2: III short 4
- S3: III short 4
- S4: III short 4
- S5: III short 9
- S6: III short 1
- S7: III short 1
- S8: III short 2
- S9: III short 1
- S10: III short 2
- S11: III short 1
- S12: III short 3
- S13: III short 2
- E1: III short 26
- E2: III short 39
- E3: III short 5
- E4: III short 4
- U1: III short 2
- U5: III short 1
- U8: III short 1
- U10: III short 1
- U11: III short 3
- M1: I short 5
- M2: I short 9
- M2b: I short 5
- M2d: I short 5
- M3: I short 3
- M3c: I short 3
- B1: I short 5
- B3: I short 3
- B4: I short 8, III short 37
- D5: I short 3, III short 3
- X1: I short 2, III short 4
- X2: I short 2, III short 3
- X3: I short 2, III short 2
- X4: I short 2, III short 3
- X5: I short 3
- P_main: I short 1, II short 12, III short 13
- F-main: I short 5, III short 28
- F-free: I short 5, III short 2
- F-gold: I short 2, III short 5
- F7bj: III short 2
- F8: III short 12
- stage A: III short 16
- BA: I short 1, III short 14
- AA: III short 2
- BB: I short 1, III short 2
- RF: I short 2, III short 8
- RB: III short 4
- GLS / hybrid GLS: III short 4
- normalisation (1)/(2)/(3): III short 8
- direct-decay / conversion route: III short 12
- A1: III short 1
- A2: III short 2
- A3: III short 1
- A4: III short 1
- A5: III short 1
- DR: III short 7
- DR_edge: III short 4
- DR_u: III short 5
- a_P: I short 4, III short 28
- a_B: III short 7
- a_eta: I short 1, III short 24
- Y_peak: I short 2, III short 16
- Y_bump: I short 3, III short 13
- f_eta(edge) / f_eta(tag) / f_eta(u-edge): III short 29
- kappa_B: III short 7
- rho_d: III short 4
- k_E1 / k_E2: III short 11
- k_cross: III short 3
- k_mix: II short 6
- BR_G / BR_A / BR_C: I short 1, III short 15
- R_34: III short 1
- MC(Delta): III short 1
- Amendment 1: III short 3
- flags G, M (sculpting): II short 1
- best-F / loose / pure rules: II short 2
- LS (like-sign, figure legend): II short 1
- first study: II short 2, III short 1
- C_sec: I short 1, III short 8
- f_win: II short 6, III short 10
- f_eta: II short 9, III short 71
- C_vtx / eps_vtx: II short 7, III short 18
- D1-D40 (all): I short 24, II short 35, III short 243
- D1: III short 3
- D2: III short 1
- D3: III short 4
- D5: I short 3, III short 3
- D6: I short 1, III short 2
- D7: III short 1
- D8: III short 2
- D9: III short 2
- D10: III short 1
- D11: III short 5
- D12: III short 2
- D14: III short 1
- D15: III short 1
- D16: III short 1
- D17: III short 1
- D18: III short 1
- D19: III short 2
- D20: III short 3
- D21: II short 5, III short 11
- D22: III short 8
- D23: III short 3
- D24: III short 2
- D25: I short 4, II short 3, III short 14
- D26: II short 2, III short 4
- D27: III short 1
- D28: II short 2, III short 7
- D29: III short 8
- D30: III short 6
- D31: III short 10
- D32: I short 1, II short 2, III short 9
- D33: II short 3, III short 4
- D34: I short 1, III short 3
- D35: I short 1, III short 22
- D36: III short 11
- D37: I short 4, II short 3, III short 15
- D38: I short 4, II short 3, III short 31
- D39: I short 1, II short 6, III short 7
- D40: I short 4, II short 6, III short 31
- WP2: II short 2, III short 3
- WP3: III short 1
- WP5: III short 2
- WP6: III short 1
- WP7: III short 8
- WP8: III short 2
- spec. Fig./Table <label>: III short 2
- package edgeu: I short 1, III short 4
- package eestep: II short 1, III short 2
- package peaktail: I short 1
- package e1v: III short 1
- package vw2a: II short 1, III short 1
- package fix: I short 1, II short 3, III short 1

