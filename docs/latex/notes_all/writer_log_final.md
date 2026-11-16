# WRITER log -- report_all.tex and report_all_short.tex (Part III)

Date 2026-11-14. Worktree /home/przygoda/HADES/FAT/FAT (branch pp45_pippimepem). No git command that changes state was run; no
macro of any package was run; no file of another agent was edited (tables_all/*.tex were only \input, figures only included).
Files owned and written: docs/latex/report_all.tex, docs/latex/report_all_short.tex (and their build products .aux/.log/.out/.toc/.pdf).
Working copies of the parts of the long report: scratchpad/final/writer_all/p1_front.tex ... p11_end.tex (report_all.tex = their
concatenation, in this order: p1_front p2_scope p3_before p4_cr p5_val p6_fit p7_sb p8_syst p9_phys p10_x17 p11_end).
Shell safety: one `rm -f` of my own (not yet existing) page-preview PNGs in scratchpad/final/writer_all/pg was issued once by mistake;
it deleted nothing (the files were created by the next command). No other destructive command.

## 1. What was done (step by step)
1. Read the spec (SPEC_report_all.md incl. Amendment 1), all results notes of work_all (wp0 with v0 verification, wp1, wp2 with its
   verification section 7, wp3, wp4, wp5 with v5, wp6, wp7, wp8 with v8, x17ab, x17c, x17d with vx), the review fixes of the data report
   (work/writer_log.md "Review fixes"), report_exp.tex (notation, section numbers), and the style models (FAT_simwt report.tex,
   report_short.tex; report_exp_short.tex).
2. Viewed the key figures (PNG previews) to write the captions; checked page sizes/rotation (ROOT PDFs with /Rotate 90 render correctly).
3. Wrote the long report in 11 parts following spec Sect. 7: status notebox; Table D1-D20 with a "status in this report" column;
   notation longtable (3 columns, with a row on symbols that have two meanings); sections 1-10; appendices A fall-back log, B files and
   commands, C normalisation details (+ Table of simulated cross sections per variant), D additional figures, E reconciliation table.
4. Tables: package snippets are \input (wp0_acceptance, wp0_sim_totals, wp1_A2, wp1_B1, wp3_tableC1/C2/C3, wp4_V1, wp5_tableE1/E2,
   wp6_tableE3, wp7_tableS1, wp7_periods, wp8_tableP1, wp8_eps2, v8_tableP2_corr, x17ab_toycls, x17c_nuisances, x17c_expected, x17d_D1),
   wrapped in my own captions (and \resizebox where they were too wide). Own tables (numbers copied from the notes): selections,
   windows, K0S, D11 nested models, A1-A5, stage-A variants, verdict, F-main variants, estimator dependence (v5), WP5 acceptance,
   stage A at WP5 eta, B4/side-band variants, bump closure, edge steps, dominant variants, alternative totals, S12 (WP8), Table N1
   (own, with Amendment 1), Table P2 (own, consolidated with the v8 corrections), R_X, X17 toys, X17 scan, vx checks, fall-backs,
   reconciliation.
5. Short report (8 pages, spec Sect. 8): p1 summary notebox + condensed decisions + two-column notation; p2-3 Fig B1 + clipped
   panel (a) of Fig B2 as the zoom + short Table B1; p3-4 short Table C2 + CR verdict + Figs C1/C2; p4-5 condensed Fig E1 (top row,
   the fitted E2 class, clipped from wp5_E1.pdf) + Fig E4 (side bands, incl. the peak zoom); p5-6 yields table (E2+E3 merged) and short
   Table S1; p6-7 Table P2 short + literature + C2 angle finding; p7-8 Figs X1, X2 + Table D1 short; p8 fall-backs and open points.
6. Built both twice with pdflatex, fixed all errors, all overfull boxes > 10 pt (equations split, tables resized/p-columns) and all
   undefined references; rendered the pages to PNG and checked the layout.

## 2. Build result
- report_all.pdf: 80 pages, 0 errors, 0 undefined references; one overfull box of 8.3 pt left inside tables_all/wp8_tableP1.tex
  (WP8's own table, below the 10 pt limit; not my file).
- report_all_short.pdf: 8 pages, 0 errors, 0 overfull boxes, 0 undefined references.
- Command: `cd /home/przygoda/HADES/FAT/FAT/docs/latex && pdflatex report_all && pdflatex report_all && pdflatex report_all_short && pdflatex report_all_short`
- Sections of report_all: Decisions (p3), Notation (p5), 1 Scope and inputs (p7), 2 Before renormalisation (p14), 3 Control regions /
  stage A (p20), 4 Stress closure (p27), 5 Template fits, cross-checks (p31), 6 Side bands, main results (p39), 7 Systematics (p45),
  8 Physics results (p49), 9 X17 (p58), 10 Conclusions and open points (p70), App. A fall-back log (p71), B files (p72),
  C normalisation (p73), D additional figures (p75), E reconciliation (p79).

## 3. Flagged numbers (disagreements between sources and which value is used)
| # | item | original (source) | used in the reports (source) |
|---|---|---|---|
| 1 | sigma(pp -> eta X) | 0.437 +- 0.037 +- 0.070 +- 0.026 mb (WP7/WP8) | 0.429 +- 0.036 +- 0.079 +- 0.025 +0.035/-0.037 mb, PRELIMINARY (v8: PDG 2026, sigma_SMASH 1.0694, C_sec(peak), peak composition) |
| 2 | syst of sigma | 16.0 % (WP7 Table S1, Fig S1) | 18.5 % (v8); Table S1/Fig S1 shown as made, with a caption note |
| 3 | sigma/sigma_SMASH | 0.410 +- 0.085 | 0.401 +- 0.092 |
| 4 | BR_G, BR_A, DR, DR_F | 2.65, 3.52, 1.397, 1.169 (PDG 2020) | 2.67, 3.54, 1.373, 1.148 (v8, PDG 2026); original errors kept |
| 5 | PDG reference BR(eta -> pi+pi-e+e-) | 2.68 +- 0.11 e-4 (PDG 2020, WP5/WP8) | 2.76 +- 0.16 e-4 (PDG 2026), both stated |
| 6 | eta normalisation | C_sec (1 - f_win) = 0.842 (WP2 headline, wp2_N1.tex) | C_sec (1 - f_eta) = 0.98 +0.03/-0.04 (Amendment 1); own Table N1 |
| 7 | f_win(B+P) | 16.8 +5.0/-5.9 % (WP2) | 14.7 +- 1.1 +1.0/-2.0 % (verification), composition of the non-eta background only |
| 8 | Table P1 efficiencies | with 1.0665 mb | 0.28 % too high (v8: simulated sample 1.0694 mb); flagged, not changed |
| 9 | live time | +7.3 % | +7.1 to +7.3 % (v8) |
| 10 | literature | "no published inclusive value" (WP8 main text) | HADES 3.5 GeV 1.14 +- 0.20 / 0.93 +- 0.14 mb; Krakow talk 1.15/1.72 mb; Staudenmaier statement only for exclusive (v8) |
| 11 | f_eta(edge) stat | 0.020 (WP6, Table E3, Fig E4) | 0.024 (x 1.2 Part II calibration, as WP8 did); total 0.050 (writer arithmetic: sqrt(0.024^2+0.039^2+0.021^2)) |
| 12 | P_main model uncertainty | 1078 (old envelope, WP6/spec) | 420 (Part II edge-safe rms of "pre + vertex"); 1078 kept as "old envelope" |
| 13 | data-only Part II numbers | old envelope; E fit error | PEAK S = 4364 +- 310 +- 477, pre 4671 +- 367 +- 453, E = 1199 +- 85 +- 287 (report_exp, cited) |
| 14 | X17 sigma x BR x BR | < 10.5 nb (x17d, 1 - f_win, sigma 0.343 mb) | < 12.2 nb (vx, F-main sigma 0.400 mb) and < 13.1 nb with the adopted 0.429 mb (writer arithmetic 3.06e-5 x 0.429 mb; uncertainty of sigma not folded in) |
| 15 | X17 N_X^UL | 563 rows (x17d, a conversion) | 392-394 rows direct limit (vx); 398 after 64/63 |
| 16 | X17 MC column | 2.39e-6 (pooled) | 3.34e-6 (CB (+) C4, vx), 11 % |
| 17 | X17 a_eta prior / a_P | from WP5 (cross-check); a_P = 0.4147 | kept (frozen); WP6 a_P 0.418 would change UL(3) by about -1 % (stated) |
| 18 | WP5 resolution text | "r = 1.10 consistent with D8" | opposite to D8; no preference in the peak region (v5) |
| 19 | free mass scale s | PASS (WP5) | MARGINAL, minimum not bracketed (v5) |
| 20 | WP5 chi2 error scale | 1.61 (GLS, replicas) | 1.61 and 1.38 (v5 hybrid) both quoted, not applied |
| 21 | Fig E2 | wp5_E2.pdf (legend hides S5 data) | v5_wp5_E2_corrected.pdf |
| 22 | WP1 acceptance references | 0.958 / 1.233 / J 1.44-1.47 (spec) | valid for s = 1 (PASS there); with D8 0.949 / 1.212 / 1.431 (pre), 0.861 / 1.101 / 1.467 (S0) |
| 23 | D19 "CB 112" | 112 (spec) | geometric mean of X17P; spread rule 171 +- 6 (WP0/v0) |
| 24 | production-model item (+19/-7 %) | "analyst's decision" to leave it out (WP8) | no such decision recorded (v8); listed next to the total, open for review |
| 25 | Fig P1 | shows WP8 values (0.437 mb, PDG 2020); lacks the inclusive 3.5 GeV point | shown as made, caption flag |
| 26 | Table E3 (wp6_tableE3.tex) | f_eta(edge) stat 0.020, P_main spread 1078 | caption flag; corrected values in the text |
| 27 | Table S1 sigma_eta column / Fig S1 | before S12 completion and v8 | caption notes; S12 table of WP8 added |

## 4. Deviations of the reports from the spec outline (closest sound choice)
1. Short report p2 "Fig B1 with a zoom inset": no inset figure exists; the zoom is panel (a) of Fig B2, clipped with \includegraphics trim
   and placed next to the short Table B1.
2. Short report p4 "condensed Fig E1 (S2, S4, S5 x E2, peak zoom), Fig E2": the condensed E1 is the top (fitted E2) row of wp5_E1.pdf with
   all five slices (clipping single non-adjacent panels is not robust); Fig E2 (template-fit residuals, a cross-check) is replaced by Fig
   E4 (side-band estimators incl. the peak zoom), because by fall-back the side bands are the main result. E2 is in the long report.
3. Short report p5: Tables E2 and E3 are merged into one table of main results and cross-checks.
4. Table P2 and Table N1 are the writer's consolidated tables (all numbers from WP8/v8 and WP2/verification), not the package snippets,
   because the snippets carry superseded values (wp8_tableP2.tex: 0.437 mb, PDG 2020; wp2_N1.tex: 0.842). v8_tableP2_corr.tex is input
   as the literature table.
5. Table A2 (wp1_A2.tex) is a caption-less longtable; it is introduced by a bold paragraph heading.

## 5. Acceptance (writer task)
| item | value | verdict |
|---|---|---|
| report_all follows spec Sect. 7 (status box, D1-D20 table, notation, sections 1-10, appendices incl. fall-back log) | yes | PASS |
| report_all_short about 8 pages as spec Sect. 8, two-column notation | 8 pages | PASS (deviations 1-3 above) |
| every number from the notes/tables; disagreements flagged | Sect. 3 above; App. E of the report; red [..] flags in the text | PASS |
| corrections with procedures in plain words | Table A2, Sects. 1.5-1.6 (CB, model), 3.1, 6.2-6.3 (B4, side bands), 7.1 (combination), 8.1 (C_sec, f_win, f_eta) | PASS |
| captions with selection, variant, normalisation, L | all data figures/tables | PASS |
| build: no errors, no undefined refs, no overfull > 10 pt | 0 / 0 / largest 8.3 pt (in an agent's table) | PASS |

## 6. Caveats
- The long report is 80 pages, mostly because every main figure is included at full width with [H] placement (some figure pages have
  white space below the float).
- Writer arithmetic (flagged in the text): the total error 0.050 of f_eta(edge) with the calibrated stat error (0.020 x 1.2 = 0.024),
  and 3.06e-5 x 0.429 mb = 13.1 nb.
- The data-only numbers of Part II were not recomputed; they are cited from report_exp.tex.

---------------------------------------------------------------------------------------------------
## Review fixes (FIXER, 2026-11-14)

Scope: the findings of the three reviewers (physics/statistics, numbers, clarity) applied to docs/latex/report_all.tex and
report_all_short.tex; numbers or figures that were wrong were recomputed or redrawn into NEW files with prefix `fix` (no package file was
overwritten, no package macro re-run in place; ScanEngine.h, scan_*.C, report_exp* untouched; no git command; no rm/mv of foreign files).
The writer's part files scratchpad/final/writer_all/p*.tex are now STALE: report_all.tex was edited directly.
Backups of the pre-review versions: scratchpad/final/fixer/report_all_before_fix.{tex,pdf}, report_all_short_before_fix.{tex,pdf}.

### New files (all created in this run)
- Macros: plots/all_report/fix_numbers.py (all recomputed numbers -> tables_all/fix_numbers.txt), fix_x17_fold.py (-> tables_all/fix_x17_fold.txt),
  fix_stageA_gls.py (hybrid GLS of stage A, uses v5_fit.py read-only), fix_dump_stageA.C, fix_after.py (Figs C2/C3, before/after column),
  fix_figs.py (Fig E4 (c,d), Fig P1 (a,b,c), Fig X2), fix_S1.py (Table/Fig S1), fix_E1.py (Fig E1, v5 library), fix_figE2.C (copy of v5_figE2.C,
  colours only), fix_x17X1.C (X1 block of x17d_figs.C, colours only), fix_wp1.C (B1 legend, B3 box colours; copy of wp1_before.C), fix_wp0defs.C (Fig 1 title).
- Figures docs/latex/figs_all/: fix_definitions, fix_B1, fix_B3, fix_C2, fix_C3, fix_E1, fix_E2, fix_E4cd, fix_wp6_E4_unrot (pdftocairo copy of wp6_E4,
  used trimmed), fix_P1abc, fix_S1, fix_X1, fix_X2 (.pdf + .png); fix_x17d_X1_unrot.pdf (unused).
- Tables docs/latex/tables_all/: fix_A2.tex, fix_tableE3.tex, fix_tableS1.{tex,txt}, fix_tC1.tex, fix_tC3.tex, fix_tE2.tex, fix_periods.tex (LaTeX notation,
  no number changed), fix_tableC3_before.txt, fix_numbers.txt, fix_x17_fold.txt.
- docs/latex/notes_all/: copies of all wp*/x17* results notes, SPEC_report_all.md, verify_wp2/*.md (R3.21).
- Outputs /hdd3/przygoda/HADES/FAT_out/all_report/fix/ (logs, stageA_dump.txt); ACLiC dir aclic_fix.

### Key numbers changed
- Composition of the eta peak taken from the data (R1.2): direct-decay route a = 0.497 +- 0.059 (side bands per M(ee) class 0.592 +- 0.139 with the
  E2 closure term; X17 fit N5/N5ref 0.476 +- 0.065), conversion route 0.237 +- 0.049 (0.218 +- 0.097; N4/N4ref 0.244 +- 0.057), kappa 0.37-0.51;
  composition systematic = half the difference / a_P = +-31.8 % (X17 alone 28 %, side bands alone 46 %).
- sigma(pp -> eta X) = 0.429 +- 0.036 +- 0.152 +- 0.025 +0.035/-0.037 mb (syst 35.5 %, total 37.9 %), 0.401 +- 0.152 of SMASH; PRELIMINARY, lower value;
  alternatives 0.52 +- 0.06 (direct route), 0.25 +- 0.05 (conversion route); "3 sigma below" removed everywhere.
- Known-sign items (R1.14): r_res = 0.89 -> 0.383 mb; f_eta = 12 % -> 0.389; vertex x one-sided 0.67-0.82 mb, with the direct route 0.81-1.19 mb.
- Y_peak syst 11.4 -> 11.5 % (S8 widened); next to it 15.2 % (no R1) / 14.7 % (literal Barlow). Y_peak = 5088 +- 432 +- 584 +- 245.
- S8 widened to g = 0.47-2.6 by linear extrapolation: Y_bump 7.7, Y_peak 1.0, a_eta 9.1, Y4 7.4, Y5 1.0, sigma 1.0 %.
- Cross-check sigma values with PDG 2026 / 1.0694 mb: F-main 0.399, edge 0.409, P_main 0.414 (0.372 at r_res 0.89), tag 0.474; SMASH visible 2.046 / 9.047 nb.
- a_B(E2) = 0.312, a_B(all) = 0.299, DR = a_P/a_B(all) = 1.397 (1.373 PDG 2026); kappa = 1 everywhere 1.34 (1.32), kappa = 0.832 everywhere 1.41 (1.38);
  sigma x BR(pi+pi-pi0) 72 / 75 ub (a_B(all) / a_B(E2)), not determined.
- X17: r_UL = 0.114; 3.16e-5 with PDG 2026; 3.25e-5 (3.35e-5) with the template offset; range 3.06-3.28e-5 over normalisations (1)-(2); 3.67e-5 vector-like;
  68 % band at 17.0 MeV [1.50, 4.77]; sigma x BR x BR < 15.0 nb with the sigma uncertainty folded (log-normal, profiled in vx's profile likelihood;
  frozen limit reproduced exactly), 15.5 nb PDG 2026, 21-26 nb if the vertex efficiency is at the N - CB / raw-row rate (plain product 13.1 nb).
- Stage A with the hybrid GLS (R1.16, v5 smoothing): chi2/ndf 1.21 (WP3 eta) / 1.24 (F-main eta); D11 Delta chi2 M0->M1 3.5, M1->M2 4.2-4.5 (> 4:
  borderline), M1->M3 2.1 (A4 p = 0.72).
- Displaced windows (R1.17): n_eff 3.4, rms 1.14 +0.37/-0.59, sigma_tr 239 +254/-239, non-overlapping subsets 0-422 rows; sigma_tr(E2) = 395 rows.
- Before/after in E2 (R3.11): B 0.833 -> 0.952, P 1.066 -> 0.969, CR-U 1.196 -> 0.999, R0' 0.899 -> 0.828.
- Small: folds z -1.37/-1.22/-1.55; R_34 0.92 +- 0.18 +- 0.24 (E2 +-4.8); edge per class stat 0.029/0.043; peak SB (e) -0.4/-0.1, (f) +5.5/+6.5 %;
  ghost -8.5 to -14.5 % (J -6.6 %); material 13-15 %; tag 3090 vs Part II 3088; positive-control explanation; vertex survival 0.565 (P) / 0.677 (all,
  0.62-0.72 per day); F-main stat of the yields in the short report 985/359.

### Findings: status
Reviewer 1: 1 fixed (no 3 sigma, lower value, one-sided items, tension may disappear); 2 fixed (data composition, routes, N4/N5 quoted); 3 fixed
(Part II X17-band numbers listed, "blind in the code, not in the information"; X17 notes do not cite Part II, so whether the agents read them cannot be
established); 4 fixed (caveats in 2, 8.4, 9.5, 10; fold of sigma into the nb limit); 5 fixed; 6 fixed (status box, verdict notebox, D22); 7 fixed
(next to the headline + justification, D24); 8 fixed (A3 FAIL for G_eta, S8 widened, stage-A borderline; S8 widening is a linear extrapolation, not
a refit: stated); 9 fixed (r_UL, editions); 10 fixed (3.25e-5 quoted next to the headline; not folded into the likelihood); 11 fixed; 12 fixed (range
3.06-3.28, D30); 13 fixed (softened; no Delta z test run: diagnostic left as an open item); 14 fixed (known-sign table both directions); 15 fixed
(reasoning corrected; the binding f_eta kept, 12 % quoted as a known-sign item instead of widening the error: the amendment is binding and the
measurement puts 12 % at 2.3 sigma); 16 fixed (hybrid GLS computed); 17 fixed (rms uncertainty, subsets, E2 sigma_tr applied to the composition);
18 fixed (labels MAIN / NOT DETERMINED / DERIVED); 19 fixed (P1 (a)-(c), E4 (c)-(d), S1 and E3 regenerated as fix_* files).
Reviewer 2: 1-17 all fixed (1 a_B(E2)/a_B(all), kappa conventions, 75 ub relabelled and recomputed; 2 PDG note + App E row; 3 one main yield named,
165 % added; 4 panel (d) redrawn; 5 985/359; 6 one edition per table/column, labelled; 7 Part II CB 101.7 (about +-14 from Poisson; the "+-13.6" of the
review was not found in Part II) cited; 8 both MC values explained, WP7 value used in S1/P2/short; 9 WP5 values everywhere; 10 WP7 z values;
11 x 1.2; 12; 13; 14 versions labelled; 15 both; 16; 17 both).
Reviewer 3: 1, 2, 4-13, 15-21 fixed; 3 partly: one palette defined in the notation and the figures that broke it were redrawn (Fig 1 title, B3 boxes,
C2, C3, E1, E2, X1, E4 (c,d), P1 (a-c)); older package figures not named by the reviewer keep their own legends (stated in the notation); 14 fixed
for 60+ captions of the long report and all captions of the short report (a few procedure/appendix figures carry "No fall-back involved");
17: Fig 1, E1 legend, E4 (c), Table E3 labels, B1 sqrt fixed; the B4 (other variables) figures were not redrawn.
Skipped: none entirely. Partial: R3.3 (only the listed figures redrawn), R1.8 S8 (extrapolated, not refitted), R1.10 (quoted, not folded into the
likelihood), R1.13 (no diagnostic test run).

### Build
- report_all.pdf: 92 pages (was 80), 0 errors, 0 undefined references, one overfull box 8.3 pt inside WP8's own wp8_tableP1.tex (unchanged).
- report_all_short.pdf: 10 pages (was 8; the corrections box and the completed notation were requested), 0 errors, 0 overfull, 0 undefined.
- Commands: cd /home/przygoda/HADES/FAT/FAT/docs/latex && pdflatex report_all && pdflatex report_all && pdflatex report_all_short && pdflatex report_all_short
- Fix macros: see report_all App. B ("review fixes") and the headers of plots/all_report/fix_*.

---------------------------------------------------------------------------------------------------
## Final integration and review (FINAL FIXER, 2026-11-14)

Scope: the findings of the two final reviewers (R1: numbers across the four documents; R2: physics, logic and over-claiming) on the
four documents after the integration of the event-vertex study (notes_all/ve_verify_note.md, ve_data_note.md, ve_sim_cause_note.md).
All findings were applied; none was skipped. No git command was run. Backups made before any edit:
scratchpad/final/integration/backup/{report_exp,report_exp_short,report_all,report_all_short}_before_final.{tex,pdf}.
Edits were made with exact-string replacement scripts (every old string had to occur exactly once):
scratchpad/final/integration/scripts/fix2_edit_lib.py, fix2_edit_report_all.py (+ _b ... _f, fix2_edit_g/h/i.py),
fix2_edit_report_all_short.py, fix2_edit_report_exp.py, fix2_edit_report_exp_short.py.
Shell-safety note: one `rm -f` was issued on my own preview PNGs (scratchpad/final/integration/fix2/ra-*.png, made minutes before by pdftoppm);
nothing else was deleted, no foreign file was touched, no mv/chmod/kill.

### New files (prefix fix2; no existing file overwritten)
- plots/all_report/fix2_numbers.py -> docs/latex/tables_all/fix2_numbers.txt: E errors x1.2 and the combinations, distances for 0.93/1.14 mb,
  N of sigma with the +-5 % period normalisation, the direct-decay route split by its two values, loss ratios, crossing-point distances,
  target box, resolution estimate, D32 inputs with their status.
- plots/all_report/fix2_x17_fold.py -> tables_all/fix2_x17_fold.txt (+ /hdd3/.../all_report/fix/fix2_x17_fold.log): X17 nb limit refolded with 38.3 %.
- plots/all_report/fix2_vtx_fig.C -> figs_all/ and figs_exp/fix2_vtx_fig.{pdf,png}, tables_all/fix2_vtx_fig.txt (ACLiC in
  /hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_verify/aclic_fix2, log logs/fix2_vtx_fig.log): copy of vtx_fig.C, 2 x 3 layout, text 7.5-8.6 pt at
  full width (was 4-5.5 pt), numbers of the tables (stat_g, 12-variant mean 263, E stat x1.2 and model shown apart, combination with the nominal CB added).
- plots/all_report/fix2_P1abc.py -> figs_all/fix2_P1abc.{pdf,png}: Fig. P1(c) with total 38.3 % and the weighted direct-route arrow (0.58-0.59 mb).
- plots/all_report/fix2_x17count.C (copy of scratchpad/final/integration/fix2/fix2_x17count.C): S0 rows of the X17 window, rounded edges 6424,
  exact edges (0.531249-0.556551) 6395 -- confirms R1-17.
- tables_all/fix2_A2.tex (event-vertex row), fix2_tableE3.tex (a_P syst 0.067), fix2_tableP2_corr.tex (caption: not all items raise sigma;
  N with the period term), fix2_tableS1.tex (N of a_eta 9.1 %, sigma +9.6/-10.0 %).

### Findings: status (R1 = reviewer 1, R2 = reviewer 2)
- R1-1 Table A2 vertex row: FIXED (fix2_A2.tex: measured 0.994 +- 0.021, 0.990 with nominal CB, simulation 0.9997, C_vtx 1.006 +- 0.021 one-sided).
- R1-2 a_P syst 0.048: FIXED -> 0.067 (16.0 %, sigma systematic without the composition term; 0.148 with it) in Table E3 (long) and Table 4 (short);
  section 7 adds N 8.8 % with the period term.
- R1-3 period normalisation +-5 % missing: FIXED (applied, not dropped): N of sigma +0.041/-0.043 mb (+9.6/-10.0 %), total 38.3 %,
  sigma/sigma_SMASH 0.401 +- 0.153, sigma_vis N +-0.071 (L (+) period), Table S1 N row, Table N1 row, App. C, notation, X17 fold (15.04 nb).
- R1-4 per-period eps_T from the ve_data code: FIXED, labelled "code of ve_data (overall 0.992 +- 0.020)" in report_exp and report_all.
- R1-5 sigma_vis "visible cross section [nb]": FIXED -> visible yield, visible rows per nb^-1 (Part II's Y; rows, not events) in all places
  (status box, notation, sections 5/6, Table P2 incl. SMASH/F-main visible, conclusions, short notation and Table 7).
- R1-6 three C_sec values: FIXED (Amendment 1, Table N1, App. C say sigma uses 1.0116 x 0.97 = 0.981, 1.0110 is the bump/peak mean;
  the short report keeps the single rounding 1.012).
- R1-7 range 0.58-0.70 mixing two definitions: FIXED everywhere -> 0.58-0.59 mb (weighted route), 0.70 mb (per-class side-band value);
  with pT 0.70-0.71 / 0.83-0.85; with C_vtx 95 % 0.60-0.61 / 0.72-0.73.
- R1-8 / R2-11 "3-4 times": FIXED -> 2.2-2.4 times (70.5 % and 75.9 % against 32.3 %).
- R1-9 0.9934 vs 0.994: FIXED; with E stat x1.2 (R2-12) the combination is 0.9937 +- 0.0211 -> 0.994 (App. C explains the earlier 0.9934).
- R1-10 figure/table numbers: FIXED by the new figure (event-group errors, 12-variant mean, E stat x1.2 and model apart) and corrected captions
  (S errors: stat_g only for the nominal CB; plain fit errors otherwise).
- R1-11 / R2-19 "7 cm upstream": FIXED -> 5-6 cm in front of the cell's upstream end (7-8 cm from its centre); mode -197.5, mean -207 mm.
- R1-12 period shares: FIXED with the day spans (034-045 0.31, 046-053 0.37, 059-067 0.275).
- R1-13 / R2-17 "152-fit grid of Part II": FIXED -> the verifier's grid (order 1-4, upper 0.58-0.65, lower 0.5225-0.5275, width +-0.5; 152 pairs p > 0.05).
- R1-14 V < pre: FIXED (sentence in Table 34/49 captions and in the short Part II bullet; V adds -500 < vz <= -160 with S = -19 +- 133).
- R1-15 / R2-6 smooth like-sign CB procedure: FIXED (procedure in report_exp 7.9 and report_all 8.2; notation rows in all four documents).
- R1-16 E2 share 0.423/0.425: FIXED, both named (0.425 C4/C5 shapes of the side-band fit; 0.423 WP8's sums).
- R1-17 6424 vs 6395: FIXED; recounted (fix2_x17count.C): rounded edges 6424, exact edges 6395; both documents state the edges.
- R1-18 short Table 7 missing N parts: FIXED (DR +-0.06; sigma_vis +-0.071 L (+) period).
- R1-19 / R2-22 pion-vertex efficiency base: FIXED (report_exp 6.5, Table effcmp, section 8 item 1; report_all 8.3).
- R1-20 D21 +16 %: FIXED -> +17 %.
- R1-21 figure text size: FIXED (new figure at full width in the short report; 7.5-8.6 pt measured with pdftotext -bbox).
- R2-1 "lower value": FIXED, label dropped everywhere (status boxes, 8.4, Table P2, Fig. P1 caption, 8.5, conclusions, short p.1 and Table 7);
  replaced by "one-sided +11-13 % (C_vtx <= +4 % more) and known-sign -8 to -11 %"; Table 38 caption corrected (fix2_tableP2_corr.tex).
- R2-2 "do not remove the difference": FIXED -> "the one-sided items alone do not remove it; with the per-class side-band route value and the pT
  item the value reaches 0.85 mb, compatible with 0.93 +- 0.14 mb" (8.5, short literature, App. E).
- R2-3 D32 evidence: FIXED (long notebox and short paragraph): B4 marked NOT DETERMINED, f_eta(edge) 1.4 sigma below the direct route,
  f_eta(tag) 0.454 +- 0.043 added (sits with the direct route), F-main labelled cross-check with PDG 2026 0.488, direct-route systematic
  (82 % N5; 10.9 % spread, +22 % 3-class, +55 % on the spec toys), "suggestive, not conclusive"; vtx_numbers.txt section 6 (generated file) not
  edited -- the corrected inputs are in fix2_numbers.txt section 6.
- R2-4 author/status of D32: FIXED (D32 = the analyst's interpretation in every place; D21 listed as his binding decision; the others await his ruling).
- R2-5 vertex-study choices: FIXED -> new D33 in Part III and decision (ix) in Part II; 0.990 +- 0.022 (nominal CB) quoted next to 0.994 in all four.
- R2-7 circular evidence: FIXED (independent evidence first: r >= 8 mm, elastic, CB/N, START, target box 0.15 %; position marked as partly forced
  by the selection; the angle argument rephrased).
- R2-8 eta efficiency used for non-eta ratios: FIXED (target box 0.15 % and the resolution-estimate procedure, 0.0219 x 0.35 = 0.8 %, in report_all
  2 and conclusions and in the short report).
- R2-9 "at most about 1 %": FIXED -> "about 1 % if L used no vertex; not estimated if it required the DST vertex; to be checked".
- R2-10 halo stated as fact: FIXED -> hypothesis in all four; absence of eta unexplained; elastic-test caveats added.
- R2-12 E errors not calibrated: FIXED (x1.2: E(V) 1241 +- 91, E(Nv) 21 +- 96, eps_E 0.983 +- 0.075 +- 0.028; 0.46 excluded at 6.6 sigma
  (1.14 mb: 7.6 sigma); S + E 0.961 +- 0.049; nominal combination 0.990 +- 0.022; bump T + E 1.000 +- 0.022).
- R2-13 "every period": FIXED -> consistent with 1 within 1.7 sigma (IIB lowest), source named, period-III S hint (0.700 +- 0.088) and its
  CB-dip explanation added (report_exp 2.5, 7 stability, 7.9; short Part II; report_all 8.2).
- R2-14 width argument: FIXED -> 1.65 sigma below the resolution, weak support only.
- R2-15 X17 headline vs frozen fall-back: FIXED (short p.1, conclusions item 6).
- R2-16 open point 1: FIXED ("can confirm the direct-decay route; if flat, the question stays open"; also 8.5, D32 box, short).
- R2-18 K0S reason: FIXED (ve_data box split 2.2 vs 18.1 per 1e4 rows cited as support).
- R2-20 "whole" -> "nearly all": FIXED. R2-21 short Part II eps_vtx definition: FIXED.
Skipped: none.

### Layout
- The longer status box of report_all moved to page 2 after the first pass; it was shortened and fits on page 1 again.
- The vertex figure is now [tp] in report_all and report_exp (with [H] it left half a page empty); width 0.97 textwidth (long), textwidth (short).
- Known-sign table and Table N1 rows split/shortened to remove new overfull boxes.

### Build (pdflatex -interaction=nonstopmode, three passes each)
- report_exp.pdf 95 pages (94 before this step), report_exp_short.pdf 9 (9), report_all.pdf 103 (98), report_all_short.pdf 12 (11);
  0 errors, 0 undefined references; remaining overfull boxes are the old ones (report_exp 6.5 / 9.9 pt pid tables; report_all 8.3 pt in
  wp8_tableP1.tex).

### Final headline numbers
- sigma(pp -> eta X) = 0.429 +- 0.036 (stat) +- 0.152 (syst) +- 0.025 (MC) +0.041/-0.043 (norm., incl. the +-5 % period normalisation) mb,
  total 38.3 %; 0.401 +- 0.153 of SMASH; PRELIMINARY; one-sided +11-13 % (C_vtx <= +4 % more), known-sign -8 to -11 %; routes 0.52 +- 0.06 /
  0.25 +- 0.05 mb (stat); with one-sided items 0.48 mb, direct route 0.58-0.59 mb (0.70 per-class side bands), up to 0.85 mb with pT.
- Y_peak = 5088 +- 432 +- 584 +- 245 rows; sigma_vis = Y_peak/L = 0.837 +- 0.071 +- 0.096 +- 0.040 +- 0.071 visible rows per nb^-1; a_P = 0.418 +- 0.035 +- 0.067 +- 0.025.
- Part II: S(PEAK) = 4364 +- 310 +- 477; E("pre") = 1199 +- 71 (calibrated +-85) +- 287; T("pre") = 3314 +- 105 +- 124.
- X17: BR(eta -> pi+pi- X17) BR(X17 -> e+e-) < 3.06e-5 (90 % CLs, normalisation (1), kept against the fall-back; 3.06-3.28; 3.16 PDG 2026;
  3.67 vector-like); sigma x BR x BR < 15.0 nb (15.04; 15.1 with C_vtx; 15.7 at its 95 % bound).
- f_eta = 3 +4/-3 %; C_sec(1 - f_eta) = 0.98 +0.03/-0.04 (sigma: 1.0116 x 0.97 = 0.981).
- eps_vtx(eta) = 0.994 +- 0.021 (0.990 +- 0.022 with the nominal CB; T alone 1.001 +- 0.023); simulation 0.99970; C_vtx = 1.006 +- 0.021, not
  applied, +2 % (+4 % at 95 % CL); 0.46 (for 0.93 mb) excluded by T at 23 sigma, by E at 6.6 sigma.
