# ve_data: how much eta is in the rows WITHOUT an event vertex? (experimental data, 29 PAT day files)

Date 2026-11-14. Package ve_data of the vertex-efficiency study. Data only. No repository file was edited, no git command
was run, nothing was deleted. Everything is new, in:
- macros, text results, figure: `/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/vtx_eff/ve_data/`
- big files (rows, ROOT spectra, logs, ACLiC build): `/hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_data/`

The cut vz > -500 mm stays as it is (mandatory). This note only MEASURES what it does to the eta.

---------------------------------------------------------------------------------------------------
## 0. Answer in five lines

1. The rows without an event vertex (Nv, 32 % of the e+e- rows) are NOT target events that lost their vertex. Their pion
   pairs cross at z = -200 mm and 10 mm from the beam axis: about 7 cm upstream of the LH2 target and off the axis.
2. These rows carry almost no eta. The cleanest measure is the ECAL-tagged eta: T(Nv) = 27 +- 67 rows, against
   T(V) = 3301 +- 108. The bump edge step gives the same answer: E(Nv) = 21 +- 80, against E(V) = 1241 +- 76.
3. The eta-peak content gives S(Nv) = 902 +- 411 (stat) +- 416 (model). This is a 1.6 sigma hint. T and E do not confirm it.
4. eps_vtx(eta, data) = S_V/(S_V+S_Nv):
   - from S: 0.838 +- 0.063 (stat) +- 0.056 (model);
   - from T: 0.992 +- 0.020;
   - from E: 0.983 +- 0.062 (stat) +- 0.028 (model);
   - the three combined: **0.986 +- 0.019** (chi2 2.3/2), above 0.955 at 95 % CL;
   - in the target box: 0.994 +- 0.005.
   This agrees with the simulation (>= 0.96).
5. The vertex cut does NOT explain the low sigma(pp -> eta X). With the data efficiency, 0.43 mb becomes 0.42 mb.
   Even the lowest single estimate (from S) gives only 0.49 mb. To reach 0.93 mb, the efficiency would have to be 0.44.
   That is excluded by every measure (T: by 23 sigma).

---------------------------------------------------------------------------------------------------
## 1. What was done, step by step

### Step 1: rows straight from PAT (`ve_fill.C`)
- **Input.** Each of the 29 PAT day files of `/hdd3/przygoda/HADES/files/exp/gen4/lep_betanew.list`, read directly
  (before FAT). Trees: PipPimEpEm_ID (us), PipPimEpEp_ID (pp), PipPimEmEm_ID (mm).
- **Rows kept.** Every row with isBest >= 0 and PT3 (trigbit == 8192), whatever its vertex. They are written in the
  column layout of scan::Row (plots/ScanEngine.h), so that the analysis engine can be used unchanged.
- **The quantities are computed as FAT and scan_prepare.C compute them:**
  - corrected momenta: pions `<pi>_p_corr_<pi>`, lepton slot 1 `<p1>_p_corr_ep`, slot 2 `<p2>_p_corr_em`;
  - masses: pi 139.56995, e 0.51099895 MeV/c2;
  - cor_ok: all four corrected momenta are finite;
  - M4, M(ee), M(pipi), OA, oa_lab, MM4 (T = 4530 MeV);
  - the photon: ECAL neutral slots with neutr_pid == 1, 0.8 <= beta <= 1.2 and cluster energy > 100 MeV. Of the three
    most energetic, the one with M(ee gamma) closest to 0.140 is taken. From it: m_eeg, m5;
  - pion track z, and the crossing point of the two pion straight lines (dz, r, doca; scan_prepare.C Crossing).
- **Event groups.** Adjacent rows with equal (run, event) in one tree of one file, counted on all rows.
- **Extra columns.** eVertClust_z, evtVertCand, evtVertClust, eVert_z, totalmult.

**Closure with FAT.** All 29 files took 47 s with 6 parallel jobs. The sums of the files give:

| sample | rows | PT3 | PT3 + vz > -500 | + cor_ok | cor_ok without vertex |
|---|---|---|---|---|---|
| us | 7 186 155 | 7 063 032 | **4 783 467** | **4 753 629** | 2 249 916 |
| pp | 423 230 | 416 407 | **122 875** | **121 789** | 289 723 |
| mm | 1 569 416 | 1 543 953 | **371 524** | **367 954** | 1 155 094 |

The bold numbers are exactly the FAT production row counts (scan_prepare expRows) and the cut flow of
`tables_exp/spec_vz500.txt`.

### Step 2: where do the no-vertex rows come from? (`ve_diag.C`, `ve_diag2.C`)
All numbers are for e+e- rows with cor_ok.

- **evtVertCand.** It is 0 in all 2 249 916 Nv rows. It is 1 in 99.5 % of the V rows.
- **Cluster vertex.** eVertClust_z sits in the target segments (-160 to -110 mm) for both classes. It says nothing
  about the event origin.
- **Event vertex in V.** vz has its bulk from -170 to -90 mm (the LH2 target).
- **Pion-pair crossing point** (doca < 10 mm; 3.44 M V rows, 1.18 M Nv rows):
  - V: z from -165 to -95 mm, as vz. Distance from the axis 2-6 mm.
  - Nv: z peaks at **-210 to -190 mm**, width about 25 mm. Distance from the axis peaks at **8-14 mm**.
  - The single pion track z and the lepton crossing z show the same shift.

  The no-vertex rows therefore come from about 7 cm upstream of the target centre and about 1 cm off the beam axis.
  These are not LH2 p+p events. A likely source is the beam halo on some structure at z about -200 mm; the target
  geometry has to be confirmed by the producers.
- **Boxes defined from this:**
  - target box X: -170 <= z_x < -85 mm and r_x < 8 mm;
  - upstream box U: -260 <= z_x < -170 mm.
- **Rows in the boxes:**
  - In X: V 2 182 631 rows, Nv **3 339**. When the pions cross inside the target, the vertex is missing in only 0.15 %
    of the rows.
  - In U: V 265 847, Nv 896 758.
- **Nv share of the rows per period:** I 0.31-0.37, IIA/IIB 0.37-0.38, III 0.27-0.28. The per-day values are in
  logs/diag2_us.log.

### Step 3: eta peak S, windows, null tests, K0S (`ve_ana.C`)
**Base of every class:** PT3, isBest >= 0, cor_ok. This means no vz > -160 cut, because the Nv rows have vz = -1000.

**Methods** (ScanEngine.h, res_common.h, res_bump.C, fix_checks.C):
- CB = 2 sqrt(N++ N--) per 10 MeV from the like-sign rows of the same class (k = 1).
- P_main: Gaussian frozen at 0.5439 GeV/c2 / 5.4 MeV + pol2, fitted in 0.5225-0.61.
- stat_g = the fit error x sqrt(var/N) of the peak window.
- model = the rms over the 144-fit edge-safe grid.

**Closure with the report.** The class "Vpre" (vz > -160) reproduces the report numbers exactly:
- S = 4671 +- 367 (res_peak.txt: 4671 +- 367);
- N = 77 513 and CB = 7949 in the peak window;
- T = 3314 +- 105 +- 124;
- E = 1199 +- 71 +- 287.

**Windows** (N, CB, N - CB):

| class | N peak | CB peak | N - CB peak | CB/N | N bump | CB bump | N - CB bump | CB/N |
|---|---|---|---|---|---|---|---|---|
| V | 86 011 | 11 488 +- 112 | 74 523 +- 319 | 0.134 | 170 751 | 22 660 +- 185 | 148 091 +- 469 | 0.133 |
| Nv | 66 194 | 38 664 +- 214 | 27 530 +- 339 | 0.584 | 128 570 | 77 194 +- 360 | 51 376 +- 524 | 0.600 |
| all | 152 205 | 50 258 +- 241 | 101 947 +- 466 | 0.330 | 299 321 | 100 170 +- 405 | 199 151 +- 704 | 0.335 |

**P_main fits:**

| class | S | stat (stat_g) | model | chi2/ndf | P_pois S | free fit: mean, sigma |
|---|---|---|---|---|---|---|
| V | **4657** | 390 (397) | 449 | 33.5/31 | 4646 | 0.54316 +- 0.00051, 5.42 +- 0.61 MeV |
| Nv | **902** | 401 (411) | 416 | 35.5/31 | 165 | no minimum (no peak to fit) |
| all | 5533 | 560 (572) | 635 | 38.5/31 | 4805 | 0.54380, 4.94 MeV |
| V_X (target box) | 2370 | 204 | 242 | 28.5/31 | 2433 | 0.54392, 5.25 MeV |
| Nv_X | 14 | 12 | 8 | 38.1/31 | -5 | - |
| V_U | -168 | 118 | 69 | 36.5/31 | -78 | - |
| Nv_U | 157 | 224 | 202 | 24.5/31 | -109 | - |
| Nv_rest | 718 | 342 | 237 | 28.8/31 | 280 | - |
| V, isBest == 1 | 4167 | 349 | 448 | | | |
| Nv, isBest == 1 | 767 | 361 | 413 | | | |

S(V) + S(Nv) = 5560 against S(all) = 5533.

**Null tests and robustness of S in Nv:**
- The frozen fit on the like-sign rows of Nv gives:
  - e+e+: -256 +- 109;
  - e-e-: 73 +- 283;
  - **N++ + N--: -183 +- 303 (-0.6 sigma)**.

  For V the like-sign sum gives -89 +- 156. The like-sign spectra have no structure at m_eta.
- CB x k on Nv, for k = 0.85 to 1.15: S(Nv) moves from 793 to 1011.
- CB bins of 5, 10 and 20 MeV: S(Nv) = 924, 902 and 577.
- S(Nv) by M(ee): 499 +- 260 (< 0.03), 289 +- 244 (0.03-0.10), 91 +- 182 (> 0.10).
- At displaced positions S(Nv) is largest, 983 +- 403, at 0.545.

**Caveat: the CB in Nv.** N - CB in M(ee) becomes negative above 50 MeV (N/CB = 0.90-0.92 at 60-100 MeV). In V and in
V_U, N >= CB there. The like-sign CB of Nv is therefore about 10 % too high at large M(ee). The M4 fit is insensitive to
this (the k test above), but a CB-based N - CB of Nv (window contents) is not a signal.

### Step 4: ECAL-tagged eta T and bump edge step E (`ve_ana.C`, `ve_checks.C`)
**The ECAL columns exist in PAT** (neutr_*_1..5).

**T** (same definition as the report):
- V: 3301 +- 108 +- 125 (free fit: 0.5462, 6.15 MeV).
- **Nv: 27 +- 67 +- 19**.
- V_X: 1776; Nv_X: 0 +- 4; V_noX: 1494; Nv_noX: 26 +- 67.

**Is the tag usable in Nv?** The pi0 peak in M(ee gamma) is present in Nv at the same place. But the pi0 tags per
(N - CB) of M(ee) < 0.10 rows are:
- V: 0.119;
- Nv: 0.035;
- V_U (upstream rows WITH a vertex): 0.065;
- V_X: 0.128.

So the low tag rate goes with the upstream source, whose pairs are mostly conversions: N - CB in M(ee) peaks at 9-11 MeV
in Nv, against 6-9 MeV in V. A vertex-dependent ECAL loss cannot be excluded completely, which is why E is the
ECAL-free cross-check.

**E** (four M(pipi) slices, the edge fixed):

| class | E | slice 0.30-0.34 | 0.34-0.37 | 0.37-0.39 | 0.39-0.40 | control 0.42-0.46 |
|---|---|---|---|---|---|---|
| V | 1241 +- 76 +- 319 | 355 +- 48 | 402 +- 43 | 360 +- 30 | 124 +- 27 | 36 +- 19 |
| **Nv** | **21 +- 80 +- 186** | -29 +- 52 | 4 +- 45 | 41 +- 31 | 5 +- 27 | 5 +- 24 |

With the same model variation applied to both classes (18 variations), eps_E = 1.003, rms 0.028 (0.965-1.044).

### Step 5: efficiencies

eps = X_V / (X_V + X_Nv). The classes are disjoint whole events, so the errors are independent.

| measure | eps_vtx(eta) | Nv/V ratio |
|---|---|---|
| S (P_main) | **0.838 +- 0.063 (stat_g) +- 0.056 (model, 116 grid pairs)**; CB x0.97/x1.03: 0.841/0.834; CB 20 MeV: 0.889 | 0.194 +- 0.120 |
| T (ECAL-tagged) | **0.992 +- 0.020** (isBest==1: 1.006 +- 0.022) | 0.008 +- 0.021 |
| E (edge step) | **0.983 +- 0.062 +- 0.028** (isBest==1: 0.971 +- 0.062) | 0.017 +- 0.071 |
| common ratio from S, T, E | **0.986 +- 0.019**, chi2 2.3/2, > 0.955 at 95 % CL | 0.014 +- 0.020 |
| S and E only (no ECAL) | 0.941 +- 0.054, > 0.86 at 95 % CL | 0.062 +- 0.061 |
| target box X, from S / T / E | 0.994 +- 0.005 / 1.000 +- 0.002 / 0.998 +- 0.005 | |
| outside the box, from S / T | 0.720 +- 0.098 +- 0.091 / 0.983 +- 0.044 | |
| isBest == 1, from S | 0.845 +- 0.063 +- 0.070 | |

**Expected Nv content, if the cut kept only the share of the rows** (eps 0.677, Nv/V ratio 0.477):
- S(Nv) about 2222;
- T(Nv) about 1575 (observed 27 +- 67, i.e. 23 sigma away);
- E(Nv) about 592.

**Expected Nv content at the simulated eps 0.96** (Nv/V ratio 0.042):
- S about 194 (observed 902 +- 580, 1.2 sigma);
- T about 138;
- E about 52.

**The no-vertex rows carry eta at the low rate that the simulation expects, not at the rate of the rows.**

### Step 6: per period and per day
Columns:
- L per period, in nb^-1, from res_common.h LumiPerDay;
- stat = stat_g;
- Y in rows per nb^-1, x 64/63.

| period | L | S_V | S_Nv | eps_S (stat, model) | T_V | T_Nv | eps_T | E_V | E_Nv | eps_E | Y_V | Y_Nv |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| I (34-46) | 2809.5 | 2358 +- 272 | 25 +- 276 | 0.990 +- 0.115 +- 0.094 | 1393 | -18 | 1.013 +- 0.033 | 560 +- 52 | 2 +- 54 | 0.996 +- 0.096 | 0.853 | 0.009 |
| IIA (47-48) | 372.6 | 291 +- 103 | 169 +- 120 | 0.632 +- 0.185 +- 0.148 | 193 | 1 | 0.993 +- 0.097 | 71 +- 20 | -11 +- 23 | 1.18 +- 0.45 | 0.793 | 0.462 |
| IIB (49-53) | 1078.0 | 807 +- 175 | 195 +- 200 | 0.805 +- 0.164 +- 0.072 | 661 | 55 | 0.924 +- 0.045 | 263 +- 33 | 43 +- 39 | 0.859 +- 0.110 | 0.761 | 0.184 |
| III (59-67) | 1819.1 | 1202 +- 207 | 515 +- 196 | 0.700 +- 0.088 +- 0.041 | 1058 | -6 | 1.006 +- 0.030 | 346 +- 40 | -15 +- 38 | 1.04 +- 0.12 | 0.671 | 0.288 |

- The S hint of Nv sits mostly in period III (2.6 sigma stat). There, T and E show nothing.
- With the target box, eps_S per period is 0.990 +- 0.007, 0.999 +- 0.022, 1.019 +- 0.018 and 0.993 +- 0.007.
- **Per day:** S_V is 25-290 per day with errors of 46-90, so eps per day is not meaningful. The table of all 29 days is
  in ve_ana.txt (section E). The per-day yields are constant:
  - V: chi2/ndf 24.8/28 (p 0.64);
  - Nv: chi2/ndf 40.2/28 (p 0.06).

### Step 7: K0S -> pi+pi-
Method: M(pi+pi-) of the rows with M4 < 1.10, N - CB (CB per 10 MeV of M(pipi)), Gaussian + pol2 in 0.455-0.545.

**Peak in V and Nv:**
- V: Y = 13 603 +- 723 (free fit: 0.49484 +- 0.00033, 6.07 MeV).
- Nv, with the V shape: Y = 2625 +- 557 (free fit: 0.49599 +- 0.00146, 7.6 +- 1.5 MeV).

**Efficiency of the vertex cut:**
- **eps_K0S = 0.838 +- 0.030**.
- Poisson fit on N without CB: 0.789 +- 0.021.
- free fits: 0.796 +- 0.050.
- isBest==1: 0.854 +- 0.034.

**K0S per 10^4 rows:**
- V: 28.6;
- Nv: 11.7;
- Nv_U: 2.2;
- Nv_rest: 18.1.

**Reading.** The K0S, unlike the eta, does lose about 16-21 % to the vertex cut. Its pions are displaced, so the
vertex then rests on the e+e- pair alone. Even if the simulation kept every K0S, the cut would explain only a factor of
0.84, not 0.29: the data/SMASH ratio before the cut would be 0.29/0.84, about 0.35 +- 0.05.

---------------------------------------------------------------------------------------------------
## 2. What it means for sigma(pp -> eta X)
- sigma scales as eps_sim/eps_data. With eps_sim >= 0.96 and the data values:

  | eps_data | from | sigma |
  |---|---|---|
  | 0.986 | combined | 0.42 mb |
  | 0.992 | T | 0.42 mb |
  | 0.983 | E | 0.42 mb |
  | 0.838 | S alone, the most pessimistic | 0.49 mb |

- Reaching 0.93 mb would need eps_data = 0.44. Every measure excludes that: S by 4.7 sigma (stat + model), and T and E
  by far more.
- **The event-vertex cut is not the cause of the deficit.** The removed 32 % of the e+e- rows are background from an
  upstream, off-axis source, with CB/N = 0.58 in the peak window. They contain no measurable eta (T, E), so the cut is
  also right for the analysis.

## 3. Caveats
1. **S(Nv) is fragile.**
   - The CB is 58 % of N in the peak window.
   - P_pois gives 165 against 902 for P_main.
   - The free fit finds no minimum.
   - With 20 MeV CB bins it is 577.

   The 1.6 sigma hint is not confirmed by T or E. An eta from an upstream nuclear source would also show in T and E, so
   it is most likely a fluctuation or a background-shape effect.
2. **T in Nv assumes the ECAL photon selection does not depend on the event vertex.** The pi0 tag rate per (N - CB) is
   lower in Nv (0.035 against 0.119). The rows of the upstream source that do have a vertex are also low (0.065), so this
   follows the source, but it is not a proof. E is ECAL-independent and agrees.
3. **The target box tests only the events whose pions cross well** (about half of the eta in V). For the rest, S gives
   0.72 +- 0.13 and T gives 0.98 +- 0.04.
4. **The like-sign CB of Nv overestimates N by about 10 % at M(ee) > 50 MeV.** Window values of N - CB in Nv are
   therefore not signal counts. The fitted S is stable against k = 0.85-1.15.
5. **The upstream origin is inferred from track geometry only.** The material at z = -200 mm, r = 10 mm should be
   confirmed with the target drawing.
6. **Errors of eps** use stat_g (event groups). The model error of eps_S uses the same grid variation in both classes.
   The combination treats S, T and E as independent; S and T share about 30 % of rows, which is a small effect.

---------------------------------------------------------------------------------------------------
## 4. Files and exact commands
Abbreviations used below:
- S = `/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/vtx_eff/ve_data`
- W = `/hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_data`
- B = `-e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_data/aclic", kTRUE)'`

Commands:

    mkdir -p $W/aclic $W/logs; cd $W
    root -l -b -q $B "$S/ve_fill.C+(20)"                         # compile + test one file
    seq 0 28 | xargs -P 6 -I{} sh -c "root -l -b -q -e 'gSystem->SetBuildDir(\"$W/aclic\", kTRUE)' '$S/ve_fill.C+({})' > logs/fill_{}.log 2>&1"   # 47 s
    root -l -b -q "$S/ve_diag.C(\"us\")"  > logs/diag_us.log      # vertex columns, pion z (TChain, 2 min)
    root -l -b -q "$S/ve_diag2.C(\"us\")" > logs/diag2_us.log     # pion crossing z, r; Nv share per day
    root -l -b -q $B "$S/ve_ana.C+"     > logs/ana.log            # -> $S/ve_ana.txt, $W/ve_ana.root (30 s)
    root -l -b -q $B "$S/ve_checks.C+"  > logs/checks.log         # -> $S/ve_checks.txt
    root -l -b -q $B "$S/ve_checks2.C+" > logs/checks2.log        # -> $S/ve_checks2.txt (ECAL tag, M(ee), M(ee gamma))
    root -l -b -q $B "$S/ve_checks3.C+" > logs/checks3.log        # -> $S/ve_checks3.txt (periods with model, E, k scan)
    root -l -b -q $B "$S/ve_fig.C+"     > logs/fig.log            # -> $S/ve_fig.png, ve_fig.pdf

Outputs:
- rows: `$W/rows/rows_00..28.root` (2.1 GB);
- spectra: `$W/ve_ana.root`;
- results: `$S/ve_ana.txt`, `ve_checks.txt`, `ve_checks2.txt`, `ve_checks3.txt`;
- figure: `$S/ve_fig.png`, `ve_fig.pdf`. Its panels: (a, b) pion crossing point z and r; (c) eta peak, (d) tagged eta,
  (e) K0S, each as N - CB - pol2 for V and Nv; (f) eps per period from S, T and E.
