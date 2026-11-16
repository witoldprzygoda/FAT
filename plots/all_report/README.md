# plots/all_report -- shared library and histogram factory of the report "data described by simulation" (SPEC report_all)

Owner of this directory: WP0 (2026-11-14). Untracked; no commits, no branches.
Every agent of the report includes `defs.h` (frozen definitions) and reads its histograms from
`/hdd3/przygoda/HADES/FAT_out/all_report/hists_v1.root` (naming scheme below; helper `hists_io.h`).
Own macros of the other packages: `plots/all_report/<prefix>_*.C`, ACLiC build dir
`/hdd3/przygoda/HADES/FAT_out/all_report/aclic_<prefix>`.

## 1. Origin of the copied headers and every change

Source: worktree `/home/przygoda/HADES/FAT/FAT_simwt`, branch `pp45_pippimepem_betanew_sim`,
commit `7206630c2ccfda3b78017e0f5de7062c1268ebb9` ("Add the pi+pi-e+e- reports and move the guides to docs/").

| file here | copied from | md5 of the source | change |
|---|---|---|---|
| SimEngine.h | plots/eta_sim/SimEngine.h | f3603216f1a29cf4df9ea3de5ad43b18 | (1) `#include "../../src/sim_truth.h"` -> `#include "sim_truth.h"` (the local copy; src/sim_truth.h exists only on the sim branch). `#include "../ScanEngine.h"` kept: it now resolves to the DATA-branch `plots/ScanEngine.h` (scan_v2, with pip_z, pim_z, dz_pi). (2) `kFigDir` -> `/hdd3/przygoda/HADES/FAT_out/all_report/wp0/figs`, `kTabDir` -> `/hdd3/przygoda/HADES/FAT_out/all_report/wp0/tables` (were `/home/przygoda/HADES/FAT/FAT/plots/eta_sim/{figs,tables}`, i.e. inside this repository). (3) a comment next to `kLdata` (value 6080.5 unchanged) saying that it is NOT the normalisation of this report; use `alldefs::kL`, `alldefs::kNbToRaw`, `alldefs::kSimToRaw`. Nothing else. |
| sim_eta_physics.h | plots/eta_sim/sim_eta_physics.h | 4d65ce54412c3cf305f951eaed8d6c69 | none (includes `"SimEngine.h"`, i.e. the local copy; reads `rw_generator.root` and `rw_angles_eta.root` from /hdd3, read-only) |
| s5_blind_ana_fit.h | plots/eta_sim/s5_blind_ana_fit.h | d7b9a7b6be8ff449a325ec0d17133215 | none (no hard-coded output path; `SaveFig(c, base)` takes the full path) |
| s5_blind_ana_tmpl.h | plots/eta_sim/s5_blind_ana_tmpl.h | 16f1fe879e9ccafb65b8896c394ba27f | none. NOTE: it reads the half-A templates of the simulation blind test (`.../blind/templates/templates.root`, sim-centred windows). Templates of this report come from hists_v1.root (`_hA`, `_hB`). |
| sim_truth.h | src/sim_truth.h | 4e82009176ff14ab72fdcba92966796e | none |
| sim_ghostvar.h | plots/eta_sim/sim_ghostvar.h | d177b2536cff87be90f4407910384872 | none (extra copy: `GhostPair` defines the ghost variant (a)) |

Check of the data-branch ScanEngine.h extension: `scan::Row` gained `pip_z, pim_z, dz_pi` (scan_v2). The simulation rows are
read by `simeng::BindRow` into `simeng::Row` (own field list, unchanged), so reading `sim_compact.root` / `sim_compact_eta.root`
is not affected; the templated selections (`simeng::PassPre`, `alldefs::PassSel`) compile for both row types. Verified by
wp0_fill.C: all 10 510 151 owned e+e- rows (std_us 230 469, lep_us 7 091 146, eta_us 3 188 536) load, and the simulation
totals of the sim study are reproduced exactly (pre: T1+T2 701.083 nb, see wp0_numbers.txt).

## 2. Files written by WP0

| file | content |
|---|---|
| defs.h | frozen definitions (namespace `alldefs`): L, 64/63, conversions; periods; mass transfer s, r, M(pipi) scale; E(m, q); selections; layouts F10, F5g, X2D, G3, BEFORE; windows and control regions; CB parameters; simulation components and variants; names |
| hists_io.h | small reader (namespace `allio`) |
| wp0_fill.C | the histogram factory (one pass over data and simulation) |
| wp0_check.C | re-reads hists_v1.root: acceptance criteria, closure, validation figures |

Run (from /home/przygoda/HADES/FAT/FAT; about 10 min, 8 GB):

    root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_wp0", kTRUE)' 'plots/all_report/wp0_fill.C+O(-1)'
    root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_wp0", kTRUE)' 'plots/all_report/wp0_check.C+O'

## 3. hists_v1.root: directories and names

    hists_v1.root
      README                         TNamed, short form of this section
      data/mask_X2D, data/mask_B_mee 1 = data N blinded (see 3.4)
      data/<sel>/info                TNamed: selection
      data/<sel>/<sub>/<what>_<layout>
      sim/<sel>__<var>/info          TNamed: selection, variant, units
      sim/<sel>__<var>/<comp>_<layout>[_hA | _hB | _rep]
    sealed/hists_v1_x17_sealed.root  (separate file) data/<sel>/<sub>/N_X2D, data/<sel>/all/N_B_mee_R2..R4 unmasked

### 3.1 Data
- `<sel>`: `S0` (NOMINAL: pre + both pions in the vertex fit), `pre`, `dz40` (pre + dz_pi < 40 mm), `PEAK` (pre + dz_pi < 40 +
  MM4 < 2.85), `vz150`, `vz170` (S0 with vz > -150 / -170), `none` (pre without vz > -160 and without the pion-vertex cut),
  `S0_best`, `pre_best` (isBest == 1), `S0_corfb`, `pre_corfb` (cor_ok fallback: cor_ok == 0 rows with the reconstructed pi-;
  data-only test, never combined).
- `<sub>`: `all`, `fold0`, `fold1` (parity of the event group), `perI` (days 034-046), `perIIA` (047-048), `perIIB` (049-053),
  `perIII` (059-067).
- `<what>`: `N` (e+e- rows), `Npp`, `Nmm` (like-sign rows), `CB` (k = 1), `CB20` (20 MeV cells; sub all only),
  `CBrep` (sub all only: 200 replicas; for a 1D layout a TH2D with y = replica index 0..199, for X2D a TH3D with z = replica).
  Units: RAW ROWS (no 64/63, no L). N, Npp, Nmm errors = event-group errors (sqrt of the sum over event groups of k^2).
  CB errors: sub all = RMS of the 200 replicas is in CBrep; the bin error stored in CB is the LINEAR propagation of the
  like-sign counts of the 10 MeV cells (without the error of f) for every sub.
- CB construction: 2 sqrt(N++ N--) in 10 MeV M4 cells (grid from 0.30), linear interpolation to the 2.5 MeV fine bins
  (ScanEngine rule), spread over the cells of a layout by f(cell | 50 MeV band) = share of N++ + N-- of the band in the cell;
  windows not on the fine grid with fractional fine bins; BEFORE variables: CB(region) x f(bin | region). The cells of every
  layout sum to the 1D CB (closure checked to 1e-9). Replicas: every like-sign event group gets a Poisson(1) weight per replica
  (counter-based generator, reproducible), and CB is rebuilt (cells and f).

### 3.2 Simulation
- `<sel>__<var>`: `S0__<v>` and `pre__<v>` for all 15 variants; `dz40__nom`, `PEAK__nom`, `vz150__nom`, `vz170__nom`.
- `<var>`: `nom` (mode A, w x V3 code 108, s = 0.99615, r = 1), `V2` (107), `bump1aY` (116), `bumpPS` (100), `convOff` (76),
  `GS` (236), `wtsOff` (0), `modeB`, `wdb`, `ghostA` (ghost-flagged true pairs moved from T3 into T2 and its groups),
  `sK0` (s = 0.99521), `sUp` (0.99705), `sDn` (0.99525), `r089`, `r110`.
- `<comp>`: `C2 C3 C4 C5 C6 Geta Gpi0D Gpi0D_om Gpi0D_nom Gconv Gconv_om Gconv_nom Gvec Gvec_om Gvec_nom Goth Goth_om Goth_nom
  GD T12 T3 T4` (alldefs::Comp). `_om` = both pions from one omega or rho0 (pcls 3, 4), `_nom` = the rest; GD = G_pi0D + G_vec
  + G_oth; T12 = C2-C5 + T2 (= the simulated sum; X17 never inside); C6 = X17 per UNIT BR (weight 2 per row).
- Units: nb per bin = sum w / L_sim (L_sim = 221.05 nb^-1), errors = event-group errors (N_eff = content^2 / error^2).
  Expected RAW data rows = nb x `alldefs::kNbToRaw` (= 6079.2 x 63/64 = 5984.2).
- Mass transfer applied to every simulated row: M4' = s [m4_true + r (m4 - m4_true)] (rows without m4_true: s m4);
  M(pipi) x 0.9952 (slices, G3, BEFORE, map); M(ee), OA, vz, pT4 unchanged; MM4 and y4 recomputed at T = 4500 MeV (the
  simulated beam; validated: the same code at 4530 MeV reproduces the stored values to 1e-7). PEAK uses MM4(4500) < 2.85.
- `_hA`, `_hB` (S0__nom, pre__nom; layouts F10, F5g, X2D, fine_M4, B_M4, win, tot): half A = even event groups, half B = odd,
  NOT rescaled (each is about half of the full; multiply by 2 for a full-luminosity template).
- `_rep` (S0__nom, pre__nom; layouts F10, F5g, X2D, win, tot): 100 Poisson event-group replicas (TH2D y = replica index; X2D
  TH3D z = replica).

### 3.3 Layouts (`<layout>`)
| name | definition |
|---|---|
| `F10_S<j>E<e>` | M4 70 bins: 5 MeV 0.400-0.520, 2.5 MeV 0.520-0.570, 5 MeV 0.570-0.700; slice S1 < 0.34, S2 0.34-0.37, S3 0.37-0.39, S4 0.39-0.43, S5 >= 0.43 (data-scale M(pipi)); E1 M(ee) < 10 MeV, E2 >= 10 MeV |
| `F5g_S<j>` | as F10, M(ee) >= 20 MeV only |
| `X2D` | TH2D M4 0.5225-0.6100 in 2.5 MeV (35) x M(ee) 0.002-0.062 GeV/c2 in 2 MeV (30); M(ee) under/overflow filled |
| `G3` | TH3D M4 0.400-0.850 in 2.5 MeV (180) x M(pipi) edges 0, 0.330, 0.335, ..., 0.445, 0.485, 0.507, 2.0 (27) x M(ee) edges 0, 0.008, 0.010, 0.012, 0.020, 2.0 GeV/c2; superset of F10/F5g for the range, binning, slice-edge (+-5 MeV), M(ee)-split (8/12 MeV), CR-K (0.425/0.44), CR-U split at 0.42 and K0S-veto (0.485-0.507) variants. Data: sub all; simulation: S0__nom, pre__nom |
| `fine_M4` | 0.30-1.10 in 2.5 MeV (320): the 1D master spectrum |
| `B_M4` | 0.30-1.10 in 5 MeV (160) |
| `B_map` | TH2D M4 0.30-1.10 (80) x M(pipi) 0.25-1.05 (80), 10 MeV |
| `B_<var>_R<k>` | var = mpipi (0.25-1.05, 2.5 MeV), mee (0-0.5, 2 MeV), oa (0-180 deg, 1 deg), vz (-320-80 mm, 2 mm), mm4 (1.6-3.4, 10 MeV), pt4 (0-1200 MeV/c, 20), y4 (-2-2, 0.05); k = region R0 < 0.455, R1 = B, R2 0.5225-0.53364, R3 = P, R4 0.55416-0.565, R5 0.565-0.70, R6 0.70-0.90, R7 0.90-1.10 (by the row's M4) |
| `win` | one bin per window (alldefs::WinId; bin labels): B, P, SBlo, SBhi, X17c (M4 0.5312-0.5566 and M(ee) 15.5-24.5 MeV), CRU, CRK, CRL, R0-R7, Jlo, Jhi, Fit (0.40-0.70), Bef (0.30-1.10), X17P (P and M(ee) 15.5-24.5 MeV: the window of the count 5258 of D19). N and sim by the row's exact M4; CB with fractional fine bins |
| `tot` | 1 bin: all selected rows at any M4 (CB: 0.30-1.10) |

BEFORE layouts (`B_map`, `B_<var>_R<k>`) exist for data sub all and for the simulation directories `*__nom`, `*__ghostA`,
`*__wtsOff` of S0 and pre. `G3` only for data sub all and S0__nom, pre__nom.

### 3.4 Blinding (D19)
The data N of the X2D bins with M(ee) 12-28 MeV (the 2 MeV bins overlapping 13-27 MeV; all M4 bins) and of `B_mee_R2`,
`B_mee_R3`, `B_mee_R4` in 12-28 MeV is set to -1 (error 0) in hists_v1.root, in every data selection and subset. N++, N--, CB
are not masked (no signal). The unmasked cells are only in `sealed/hists_v1_x17_sealed.root`: open it only after
`x17_frozen.txt` (procedure + Asimov limit, WP9 X-c) exists. Coarse classes (F10 E2, F5g, G3 M(ee) 12-20 and >= 20 MeV, the
window X17c) are not masked: their content (window count, M(ee) > 20 MeV share) is already known (D19).

### 3.5 Model and yield conventions (spec Sect. 2, for the fitting packages)
- mu_hi = k_e CB_hi + c_j (1 + b x_i)[G_D + kappa_B G_C]_hi + a_eta [C2 + kappa (C3 + C4) + rho_d C5 + g G_eta]_hi with
  x = (M4 - 0.55)/0.1 (`alldefs::TiltX`), G_C = Gconv; fits on raw rows: sim templates x kNbToRaw.
- Data yields: (N - CB) x 64/63; cross sections (N - CB) x 64/63 / L (`alldefs::kRawToNb`).

## 4. Things every consumer must know (WP0 findings, 2026-11-14)
- hists_v1.root and the sealed file are write-protected (chmod a-w). A rerun of wp0_fill.C needs `chmod u+w` first and
  creates a new version (never overwrite v1 silently: copy it to hists_v1_old.root and say so).
- CB errors: use the replicas (`CBrep_*`). The stored bin error of `CB_*` is the linear propagation of the 10 MeV cell counts
  WITHOUT the error of f(cell | band); in F10 the replica RMS is 3.2-4.2 x the linear error in E1 and 1.2-1.6 x in E2
  (median over bins, data S0), up to 0.47 relative in S4/S5 at 0.40-0.50 (few like-sign rows in those bands).
- The spread rule assumes that N--/N++ of a 50 MeV band holds in every cell. It does not: N--/N++ is 7-14 in E1
  (M(ee) < 10 MeV) and 3.6-6.3 in E2 (S0). Summed over all cells (windows B, P, CR-U) the spread CB and 2 sqrt(N++ N--) agree
  to 0.3 %, but per cell 2 sqrt(N++ N--)/CB(spread) = 0.64-0.91 (E1), 0.87-1.12 (E2), 0.65 in the X17 windows, 0.17 in CR-K
  (N++ = 2). This is the S1 variant "per-cell geometric mean" and it is large in E1, CR-K and the X17 region
  (table in docs/latex/tables_all/wp0_check.txt).
- The count 5258 of D19 is the window X17P (P x M(ee) 15.5-24.5 MeV, S0), reproduced exactly; D19's CB 112 is 2 sqrt(14 x 223)
  of that window; the spread rule gives 171.0 +- 6.2. The spec's X17 counting window X17c (0.5312-0.5566) has N = 6424.
- CR-K upper edge: 0.5225 (F10 grid edge; the spec writes 0.522).
- The simulation has almost no rows below vz = -160 mm (0.15 % of the compact rows; data 7 %): sim vz170 = S0 within 0.05 %.
- Independent verification (v0, 2026-11-14): `v0_verify.C` recomputed 189 contents of hists_v1.root from the compact files
  (own data selections and CB code; original FAT_simwt simulation headers) -- all agree; table
  `docs/latex/tables_all/v0_verify.txt`, notes in `work_all/wp0_results.md` (section "Verification"). Comment fixes only:
  defs.h (kL: L contains the sector-2 removal S2fac, the analysis yields still need C_sec of WP2) and SimEngine.h (kLdata
  comment names kNbToRaw / kSimToRaw).
