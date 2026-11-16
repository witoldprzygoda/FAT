#!/usr/bin/env python3
# ============================================================================
# x17c_frozen.py  --  WP9 X-c (prefix x17c): writes /hdd3/przygoda/HADES/FAT_out/all_report/x17/x17_frozen.txt (the
#                     complete frozen X17 procedure and the expected limits) and the LaTeX snippets of the x17c tables.
# Inputs: x17c_config.txt, x17c_freeze.log, tables_all/x17c_{sidebands,variants,expected_asimov,expected,inject}.txt.
# Run: python3 plots/all_report/x17c_frozen.py   (after x17c_freeze.C, x17c_toys.C, x17c_summary.C)
# ============================================================================
import datetime
import hashlib
import os
import re

OUT = "/hdd3/przygoda/HADES/FAT_out/all_report/x17c"
FROZEN = "/hdd3/przygoda/HADES/FAT_out/all_report/x17/x17_frozen.txt"
TAB = "/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all"
CODE = "/home/przygoda/HADES/FAT/FAT/plots/all_report"


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read(path):
    with open(path) as f:
        return f.read()


cfg = {}
for line in read(f"{OUT}/x17c_config.txt").split("\n"):
    p = line.split()
    if len(p) == 2:
        cfg[p[0]] = float(p[1])
log = read(f"{OUT}/x17c_freeze.log")
sb = read(f"{TAB}/x17c_sidebands.txt")
var = read(f"{TAB}/x17c_variants.txt")
asi = read(f"{TAB}/x17c_expected_asimov.txt")
toy = read(f"{TAB}/x17c_expected.txt")
inj = read(f"{TAB}/x17c_inject.txt")


def grab(pattern, text, cast=float, group=1):
    m = re.search(pattern, text, re.M)
    if not m:
        raise RuntimeError(f"pattern not found: {pattern}")
    return cast(m.group(group))


# side-band fit
chi2 = grab(r"Pearson chi2 ([\d.]+) / ndf", sb)
ndf = grab(r"/ ndf (\d+)", sb, int)
n5 = grab(r"^N5 ([\d.]+) \+-", sb)
n5e = grab(r"^N5 [\d.]+ \+- ([\d.]+)", sb)
aeta_sb = grab(r"a_eta ([\d.]+) \+- [\d.]+; g", sb)
g_sb = grab(r"; g ([\d.]+) \+-", sb)
smooth_line = re.search(r"local smoothing.*", log).group(0)
masked = grab(r"background-only Asimov set: (\d+) rows, of which (\d+)", log, float, 2)
asim_tot = grab(r"background-only Asimov set: (\d+) rows", log)
bx17c = grab(r"counting window X17c: ([\d.]+) rows", log)
bx17c_cb = grab(r"counting window X17c: [\d.]+ rows \(CB part ([\d.]+)\)", log)
fill_var = {m.group(1): float(m.group(2)) for m in re.finditer(r"variant (fill\w+)\s+expected.*rows in X17c ([\d.]+)", log)}
phi_line = re.search(r"variance factor: phi_max.*", log).group(0)
sx_line1 = re.search(r"X17 yield nuisance \(1\).*", log).group(0)
sx_line2 = re.search(r"X17 yield nuisance \(2\).*", log).group(0)
rep_line = re.search(r"Asimov at 17 MeV: sigma.*", log).group(0)
prior_line = re.search(r"priors: .*", log).group(0)
ap_line = re.search(r"Y_peak sim .*", log).group(0)
ctail_line = re.search(r"conversion tail .*", log).group(0)

# Asimov table rows
arows = []
for line in asi.split("\n"):
    if not line.strip() or line.startswith("#"):
        continue
    t = line.replace("|", " ").split()
    arows.append(dict(m=float(t[0]), RX=float(t[1]), e1=[float(x) for x in t[2:7]], e2=float(t[7]), e3=float(t[8]), nX=float(t[9]), s=float(t[10])))
# toy table
trows = {}
for line in toy.split("\n"):
    t = line.split()
    if len(t) >= 9 and t[2] == "toys":
        trows[(float(t[0]), int(t[1]))] = dict(n=int(t[3]), q=[float(x) for x in t[4:9]], ratio=float(t[9]) if len(t) > 9 else float("nan"))
pull = {float(m.group(1)): (float(m.group(2)), float(m.group(3))) for m in re.finditer(r"^\s*([\d.]+)\s+background-only r_hat/sigma \(norm 1\): mean ([+-][\d.]+), RMS ([\d.]+)", toy, re.M)}
injl = [line for line in inj.split("\n") if line and not line.startswith("#")]
# toy CLs against asymptotic CLs on the Asimov set at 17 MeV/c2 (x17c_toycls.C)
tcl = []
for line in read(f"{OUT}/x17c_toycls_m17.txt").split("\n"):
    m = re.search(r"r ([\d.]+) BR [\de.+-]+ CLs_toy ([\d.]+) \+- ([\d.]+).*?(\d+) toys\) CLs_asym ([\d.]+)", line)
    if m:
        tcl.append((float(m.group(1)), float(m.group(2)), float(m.group(3)), int(m.group(4)), float(m.group(5))))
tcl.sort()
import numpy as np
_r = np.array([t[0] for t in tcl]); _c = np.array([t[1] for t in tcl]); _e = np.array([t[2] for t in tcl])
_p = np.polyfit(_r, np.log(_c), 2, w=_c / _e)
_roots = [x.real for x in np.roots([_p[0], _p[1], _p[2] - np.log(0.1)]) if abs(x.imag) < 1e-12 and _r.min() <= x.real <= _r.max()]
toy_ul_r = _roots[0] if _roots else float("nan")

date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
files = ["x17c_model.h", "x17c_lib.h", "x17c_disp.C", "x17c_freeze.C", "x17c_toys.C", "x17c_summary.C", "x17c_frozen.py", "defs.h", "hists_io.h"]
md5s = {f: md5(f"{CODE}/{f}") for f in files}
md5_lines = md5("/hdd3/przygoda/HADES/FAT_out/all_report/x17ab/x17_lines.root")
md5_hists = md5("/hdd3/przygoda/HADES/FAT_out/all_report/hists_v1.root")
md5_wp5 = md5("/hdd3/przygoda/HADES/FAT_out/all_report/wp5/results.txt")
md5_egcov = md5(f"{OUT}/x17c_egcov.root")


def a17(norm=1, q=2):
    for r in arows:
        if abs(r["m"] - 17.0) < 1e-6:
            return r["e1"][q] if norm == 1 else (r["e2"] if norm == 2 else r["e3"])


L = []
w = L.append
w("=" * 118)
w("x17_frozen.txt -- FROZEN X17 PROCEDURE of the report_all analysis (SPEC report_all, WP9 X-c; D18, D19)")
w("=" * 118)
w(f"Written {date} by package x17c (x17c_frozen.py). HADES p+p, T = 4.5 GeV (feb22). eta -> pi+pi- X17, X17 -> e+e-.")
w("This file freezes the procedure BEFORE the data of the blinded M(ee) rows are looked at (D19). It contains the complete")
w("procedure (selection, windows, model, nuisances and priors, test statistic, CLs settings, mass grid, quoted masses) and the")
w("expected limits. X-d (unblinding) must apply it unchanged; any deviation must be reported as such next to the frozen result.")
w("")
w("-" * 118)
w("0. BLINDING STATEMENT")
w("-" * 118)
w("- Masked: the data rows N of the X2D layout with M(ee) 12-28 MeV/c2 (the 2 MeV/c2 bins that overlap the 13-27 MeV/c2 of D19),")
w("  all M4 bins; in hists_v1.root they hold -1 (data/mask_X2D); the true counts are only in")
w("  /hdd3/przygoda/HADES/FAT_out/all_report/sealed/hists_v1_x17_sealed.root, which was NEVER opened by X-a, X-b or X-c.")
w("- Used by X-c: data N of the side-band rows M(ee) 2-12 and 28-62 MeV/c2 (770 cells, 197 752 rows), like-sign rows and CB")
w("  (not blinded), simulation, WP5 results. The event-group covariance (x17c_disp.C) read the data tree with a filter that")
w("  drops every row with 12 <= M(ee) < 28 MeV/c2 as its first operation.")
w("- Already known before the freeze (D19, stated in the reports): the X17P window count 5258 (CB 112 geometric / 171 spread),")
w("  the X17c counting-window count 6424 (CB 209.5), and the peak share at M(ee) > 20 MeV of 0.306. The 2 MeV/c2 M(ee) shape")
w("  of the peak band in 12-28 MeV/c2 is unknown to the analysis.")
w("- Deviation from the spec text: the side bands are 2-12 and 28-62 MeV/c2, not 2-13 and 27-62: with 2 MeV/c2 bins from")
w("  2 MeV/c2 the bins 12-14 and 26-28 MeV/c2 overlap 13-27 and are masked in hists_v1.root.")
w("")
w("-" * 118)
w("1. CODE AND INPUTS (md5)")
w("-" * 118)
for f in files:
    w(f"  {md5s[f]}  {CODE}/{f}")
w(f"  {md5_hists}  /hdd3/przygoda/HADES/FAT_out/all_report/hists_v1.root (WP0)")
w(f"  {md5_lines}  /hdd3/przygoda/HADES/FAT_out/all_report/x17ab/x17_lines.root (X17 lines and R_X(m), x17ab_line.C)")
w(f"  {md5_wp5}  /hdd3/przygoda/HADES/FAT_out/all_report/wp5/results.txt (a_eta, g, background shapes, Y_peak)")
w(f"  {md5_egcov}  {OUT}/x17c_egcov.root (event-group covariance of the side-band cells)")
w("  Frozen numbers: /hdd3/przygoda/HADES/FAT_out/all_report/x17c/x17c_config.txt (read by x17cl::ReadConfig):")
for k in ["a0", "aSig", "g0", "gSig", "kE1", "kE2", "cTail", "xSig1", "xSig2", "xSig3", "qScale", "aP", "aPSig", "bkgMode", "bshape", "fillOrder"]:
    w(f"    {k:9s} {cfg[k]:.6g}")
w("  The model is x17c_model.h = x17_model.h of x17ab (namespace x17m, validated in X-b) with four additions: normalisation")
w("  switch xNorm/xA, log-normal X17 yield nuisance, line-shape nuisances (width, M(ee) scale), variance factor qScale.")
w("")
w("-" * 118)
w("2. SELECTION AND DATA LAYOUT")
w("-" * 118)
w("- Rows: data scan_compact.root (scan_v2, all stored combinations), tree us (unlike-sign); CB from pp and mm.")
w("- Selection S0 (D1, frozen in defs.h): isBest >= 0 && vz > -500 mm && PT3 (trigbit 8192) && vz > -160 mm && cor_ok == 1")
w("  && pip_isUsedVertex == 1 && pim_isUsedVertex == 1. isBest >= 0 is nominal (D2); isBest == 1 is not used for the limit.")
w("  OA(e+e-) is NOT cut. ECAL is not used.")
w("- Layout X2D (hists_v1.root data/S0/all/N_X2D): M4 = M(pi+pi-e+e-) on the data scale 0.5225-0.6100 GeV/c2 in 35 bins of")
w("  2.5 MeV/c2 (a row goes to the X2D bin of the centre of its 2.5 MeV fine bin) x M(ee) 2-62 MeV/c2 in 30 bins of 2 MeV/c2.")
w("  Fit region = all 1050 cells (X-d). Rows outside 2-62 MeV/c2 are not in the likelihood (their share enters the totals N5, N4).")
w("- Classes: E1 = M(ee) < 10 MeV/c2 (CB factor k_E1, background shape class 0), E2 = M(ee) >= 10 MeV/c2 (k_E2, class 1).")
w("- CB = 2 sqrt(N++ N--) in 10 MeV M4 cells, spread over the X2D cells with the like-sign fractions of its 50 MeV band")
w("  (WP0 spread rule, data/S0/all/CB_X2D). Never from simulation.")
w("- Simulation templates sim/S0__nom (mode A, w x V3 = kPhysMain 108, T1 + T2, SMASH sigma_prod x PDG BR), mass transfer")
w("  s = 0.99615, r = 1 (D8), in expected RAW rows: nb x 5984.2 (L = 6079.2 nb^-1 x 63/64).")
w("- Counting cross-check window X17c (spec Sect. 2): M4 0.5312-0.5566 GeV/c2 x M(ee) 15.5-24.5 MeV/c2.")
w("")
w("-" * 118)
w("3. MODEL (expected raw rows in cell i (M4), j (M(ee)); x17c_model.h, bkgMode 4)")
w("-" * 118)
w("  mu_ij = k_e(j) CB_ij                                                       (CB, measured)")
w("        + n_j B_c(j)(i) [1 + (b_c + b'_c u_j) t_i + (c_c + c'_c u_j) t_i^2]          (other background)")
w("        + gP_c(j)(i) [N5 D_j + N4 C_j (1 + sC thC T_j)]                        (eta peak: direct + conversion)")
w("        + r Y (1 + sX)^thX gX(i) X_j(m_X; thW, thS)                            (X17 line)")
w("        + a_eta [E23_ij + g Geta_ij]                                           (eta-related non-peak terms)")
w("  t = (M4 - 0.55)/0.1, u = (M(ee) - 30 MeV)/30 MeV. n_j FREE in every M(ee) row (no simulated background per row, D18).")
w("  B_c = M4 shape of the WP5 F-main background (x17bkg E1, E2 of wp5/results.txt: c_j (1 + b x)(G_D + kappa_B G_C), exp(pol3)")
w("    smoothed per slice, summed over the slices), unit sum per class; b, c, b', c' free per class (4 x 2 parameters).")
w("  gP_c = M4 shape of the direct decays C5 (V3) of the class; D_j = M(ee) fraction of C5 (V3), C_j = of C4 (conversion),")
w("    fractions of all rows of the M4 range (M(ee) flows in the totals): N5, N4 = totals over the full M(ee) range, both free.")
w("  T_j = 1 for M(ee) >= 12 MeV/c2 (conversion tail).")
w("  E23 = C2 + C3, Geta = G_eta: factorised (row sum x exp(pol3)-smoothed M4 shape per class).")
w("  X17 line: X_j(m) = the simulated X17 line (C6, 17.0 MeV/c2) morphed to mass m (mean and width from the reweighted C5")
w("    slices, pol2 in m; at 17.0 the C6 line itself); gX = M4 shape of C6; R_X(m) = eps_X/eps_5 = 1.3455 x eps(m)/eps(17)")
w("    (x17_lines.root). Isotropic e+ decay angle, 3-body phase space for M(pi+pi-) (nominal X17 hypothesis).")
w("  POI r >= 0. Three normalisations (the same fit, only Y and sX change):")
w("    (1) PRIMARY: Y = N5 R_X (N5 fitted):            BR(eta -> pi+pi-X17) BR(X17 -> e+e-) = r x 2.68e-4")
w(f"    (2) Y = a_eta N5ref R_X, a_eta = {cfg['a0']:.4f} FIXED in Y (WP5), N5ref = 7458.8 direct rows at SMASH: BR = r x 2.68e-4")
w(f"        (= N_X / (a_eta s_X2D), s_X2D = N5ref R_X / 2.68e-4 = C6 per unit BR = 3.745e7 rows at 17 MeV/c2)")
w(f"    (3) as (2) with the peak-based eta scale a_P = Y_peak(data)/Y_peak(sim) = {cfg['aP']:.4f} (WP5 F-main Y_peak 5045, sim 12165)")
w("    In (2)/(3) the eta-yield uncertainty is inside sX; the a_eta parameter of the E23/Geta terms is NOT tied to Y.")
w("  Fall-back of the spec (N5 error > 50 % -> (2) primary): NOT fired (N5 error 20.3 % in the side-band fit). (1) stays primary.")
w("")
w("-" * 118)
w("4. NUISANCE PARAMETERS AND PRIORS (Gaussian penalty terms in -2 ln L; global observables randomised in toys)")
w("-" * 118)
w(f"  k_E1      CB factor, M(ee) < 10 MeV        1 +- {cfg['kE1']:.2f}  (x17ab: free k_E1 is degenerate with the free n_j)")
w(f"  k_E2      CB factor, M(ee) >= 10 MeV       1 +- {cfg['kE2']:.2f}  (D5)")
w(f"  a_eta     eta-related terms                {cfg['a0']:.4f} +- {cfg['aSig']:.4f}  (WP5 F-main 0.3823 +- 0.0431 x sqrt(chi2/ndf) 1.61 (+) S3 quadratic shift 0.149)")
w(f"  g         G_eta factor                     {cfg['g0']:.3f} +- {cfg['gSig']:.3f}  (WP5 1.717 +- 0.401 x 1.61 (+) S3 shift 0.711)")
w(f"  thC       conversion tail (>= 12 MeV)      N(0,1) x sC = {cfg['cTail']:.3f} (RMS of the 100 C4 replicas)")
w("  thW       X17 line width                   N(0,1) x 10 % (vertical morphing, half difference of width x 1.1 / 0.9)")
w("  thS       M(ee) scale of the X17 line      N(0,1) x 2 %")
w(f"  thX       X17 yield, LOG-NORMAL (1+sX)^thX sX(1) = {cfg['xSig1']:.4f}, sX(2) = {cfg['xSig2']:.4f}, sX(3) = {cfg['xSig3']:.4f}")
w("    composition (quadrature):")
w("      " + sx_line1.replace("x17c: ", ""))
w("      " + sx_line2.replace("x17c: ", ""))
w("      DL = R_X with direct-like M(pi+pi-); NP = generator pion angle (max over 10-40 MeV/c2); valid = slice/C6 efficiency")
w("      check of x17ab (1.043 +- 0.030); BR_dir = 0.11/2.68; N5bias = half spread of N5 over the background-shape and CB")
w("      variants of the side-band fit (Bsim3 and noMeeDep at the ends); dir = V2/GS direct-decay templates (V2 gives eps_5 +16 %);")
w("      width/scale = peak templates r = 0.89/1.10 and s +- 0.0009 (change of N5_hat R_X, resp. N5ref R_X).")
w("  Free: n_j (30), N5, N4, b_c, c_c, b'_c, c'_c (8), r.")
w("  NOT in the likelihood, quoted separately:")
w("    - e+ decay angle 1 + cos^2 instead of isotropic: R_X -16.5 %, limit x 1.198 (model dependence of the X17 hypothesis);")
w("    - MC statistics of the templates: X-d recomputes r_hat and the limit with the 100 C4 replicas (C4_X2D_rep) and with")
w("      the 200 CB replicas (CBrep_X2D) and quotes the RMS of the UL as 'MC' (D13). On the Asimov set at 17 MeV/c2:")
w("      " + rep_line.replace("x17c: ", ""))
w("")
w("-" * 118)
w("5. VARIANCE FACTOR (event groups) AND GOODNESS OF FIT")
w("-" * 118)
w("- One event gives several rows; the cell counts are correlated. x17c_disp.C built the event-group covariance S of the side-band")
w("  cells; phi = w'Sw / w'diag(N)w for X17-like weights w = s/mu (lines at 8-44 MeV/c2 restricted to the side-band cells).")
w("  " + phi_line.replace("x17c: ", ""))
w(f"- FROZEN: q~ is divided by qScale = {cfg['qScale']:.4f} in every CLs computation (observed, Asimov, toys); toys are drawn as")
w("  qScale x Poisson(mu/qScale). Within-cell dispersion is 1.002; var(total)/N of all side-band cells 1.078.")
w(f"- Side-band fit (r = 0): Pearson chi2/ndf = {chi2:.1f}/{ndf} = {chi2/ndf:.3f} -> fall-back 'Pearson chi2/ndf of n_j > 1.5' NOT fired.")
w("  X-d: if the full fit (all rows) has Pearson chi2/ndf > 1.5, qScale := max(qScale, chi2/ndf) (spec fall-back, frozen rule).")
slope1 = grab(r"norm 1: offset .* = ([\d.]+)", inj)
w("- Fall-backs of the spec (Sect. 6, X17 row), status at the freeze:")
w(f"    Pearson chi2/ndf of the n_j fit > 1.5: NOT fired ({chi2/ndf:.3f}, side bands).")
w(f"    slope off by > 5 %: NOT fired (injection toys, normalisation (1): slope {slope1:.4f}, statistical error about 0.02).")
w(f"    N5 error > 50 %: NOT fired ({100*n5e/n5:.1f} %).")
w("")
w("-" * 118)
w("6. TEST STATISTIC, CLs, LIMITS")
w("-" * 118)
w("- Likelihood: Poisson of all 1050 cells + Gaussian penalties; minimiser Fisher scoring with step halving (x17c::Model::Fit).")
w("- Test statistic q~_r (Cowan, Cranmer, Gross, Vitells, EPJC 71 (2011) 1554): r_hat may be negative in the free fit; reference")
w("  likelihood at max(r_hat, 0); q~_r = 0 for r_hat > r; q~ / qScale.")
w("- 90 % CL upper limit: CLs = CL_{s+b}/CL_b = 0.10, asymptotic formulae (eqs. 65-66 of Cowan et al.), sigma_A from the")
w("  background-only Asimov set built from the r = 0 fit of the OBSERVED data (x17c::LimitAsym).")
w("- Toy CLs (x17c::ToyCLs, 1000 toys per point and hypothesis, conditional fits at r and 0, global observables drawn) at the")
w("  three quoted masses; if toy and asymptotic UL differ by more than 10 %, the TOY limit is primary (spec X-b; x17ab found")
w("  0.97-1.05). If CL_b = 0 (observed q~ above every toy), report it and use the asymptotic CLs.")
w("- Feldman-Cousins cross-check (spec X-d): unified interval in r from toys of the profile-likelihood ratio, at 17 MeV/c2.")
w("- Local p0 at every mass (q0 of the r_hat > 0 fit, asymptotic, / qScale); global p from background-only toys over the mass")
w("  grid (look-elsewhere) if the largest local Z > 3.")
w("- Counting cross-check in X17c: B_pred = background-only expectation of the frozen model, i.e. the side-band fit (r = 0)")
w("  with the masked rows predicted (below), integrated over X17c with the overlap fractions of the X2D cells:")
w(f"    B_pred(X17c) = {bx17c:.1f} rows (CB part {bx17c_cb:.1f}); masked-row prediction variants: " + ", ".join(f"{k} {v:.1f}" for k, v in fill_var.items()))
sigB = max(abs(fill_var[k] - bx17c) for k in ("fill1", "fill3", "fillSmooth") if k in fill_var)
w(f"    FROZEN sigma_B(X17c) = max |B_variant - B_pred| over fill1, fill3, fillSmooth = {sigB:.1f} rows (the like-sign and flat")
w("    references do not follow the background M(ee) shape; they are listed for information only).")
w("    X-d: excess = N(X17c) - B_pred; error = sqrt(N) (+) sigma_B; BR = excess / (a_eta s1), s1 = C6 per unit BR in X17c")
w("    = 3501.6 nb x 5984.2 = 2.0954e7 raw rows at a_eta = 1. The counting check is weak by construction: its background comes")
w("    from the interpolation across the masked rows, which the 2D fit (n_j free per row) does not need.")
w("")
w("-" * 118)
w("7. MASS GRID AND QUOTED MASSES")
w("-" * 118)
w("- Scan: m_X = 10-40 MeV/c2 in 1 MeV/c2 steps, 0.5 MeV/c2 steps in 15-20 MeV/c2, plus 16.7 and 16.85 (38 masses).")
w("- Quoted: 16.7, 16.85 and 17.0 MeV/c2 (17.0 = the simulated X17; the others use the morphed line and R_X(m)).")
w("- Primary number of the report: normalisation (1), isotropic X17 decay, at the three quoted masses, with the toy CLs check.")
w("")
w("-" * 118)
w("8. HOW THE EXPECTED LIMITS WERE MADE (X-c)")
w("-" * 118)
w("1. The model above was fitted with r = 0 to the data side-band cells only (rows 12-28 MeV/c2 masked in the likelihood).")
w(f"   Result: N5 = {n5:.0f} +- {n5e:.0f} rows, a_eta = {aeta_sb:.3f}, g = {g_sb:.2f}; chi2/ndf {chi2/ndf:.3f}.")
w("2. The background n_j of the 8 masked rows was predicted: n_j = sim_j x exp(pol2(y)), sim_j = simulated G_D + G_C rows,")
w("   the pol2 in ln(n_j/sim_j) fitted to the side-band rows (weights: fit error of n_j (+) simulation error of sim_j).")
w(f"   Predicted {masked:.0f} rows in the masked rows (all components {asim_tot:.0f} rows in X2D).")
w("   " + smooth_line.replace("x17c: ", "(variant reference) "))
w("3. Background-only Asimov set = the model with these parameters, r = 0, on all 1050 cells.")
w("4. Asimov expected limits: x17c::LimitAsym on the Asimov set (bands from sigma_A), every mass, the three normalisations.")
w("5. 2000 background-only toys at each quoted mass (500 at 10, 25, 40 MeV/c2): pseudo-data from the Asimov expectation,")
w("   over-dispersed by qScale, constraint centres drawn; each toy fitted on all cells like the data; its asymptotic CLs UL.")
w("   The quantiles 2.3/15.9/50/84.1/97.7 % are the expected band from toys.")
w("6. Injection toys at 17 MeV/c2 (BR 1, 3, 6 x 1e-5, 500 each) test linearity and coverage of the frozen fit.")
w("")
w("-" * 118)
w("9. EXPECTED 90 % CLs UPPER LIMITS on BR(eta -> pi+pi-X17) BR(X17 -> e+e-)  [S0, isBest >= 0, L = 6079.2 nb^-1]")
w("-" * 118)
w("  FROZEN choice for the report: at the quoted masses the expected limit is the TOY median with the toy +-1/2 sigma bands")
w("  (they contain the non-asymptotic behaviour; toy/Asimov median 0.94 at 16.7-17.0 MeV/c2, 0.98-0.99 at 25-40); versus m_X the")
w("  Asimov median and bands are drawn (Fig. X2), with the toy points overlaid.")
w("  Quoted masses (normalisation (1) primary):")
w(f"  {'m [MeV/c2]':>10s} {'norm':>4s} {'method':>7s} {'-2s':>10s} {'-1s':>10s} {'median':>10s} {'+1s':>10s} {'+2s':>10s}")
for m in (16.7, 16.85, 17.0, 10.0, 25.0, 40.0):
    for norm in (1, 2, 3):
        r = next((x for x in arows if abs(x["m"] - m) < 1e-6), None)
        if r and norm == 1:
            w(f"  {m:10.2f} {norm:4d} {'Asimov':>7s} " + " ".join(f"{v:10.3e}" for v in r["e1"]))
        elif r:
            w(f"  {m:10.2f} {norm:4d} {'Asimov':>7s} {'':10s} {'':10s} {(r['e2'] if norm == 2 else r['e3']):10.3e}")
        t = trows.get((m, norm))
        if t:
            w(f"  {m:10.2f} {norm:4d} {'toys':>7s} " + " ".join(f"{v:10.3e}" for v in t["q"]) + f"   ({t['n']} toys; toy/Asimov median {t['ratio']:.3f})")
for m, (pm, pr) in sorted(pull.items()):
    w(f"  background-only toys at {m:.2f} MeV/c2: r_hat/sigma mean {pm:+.3f}, RMS {pr:.3f}")
w("")
w("  Toy CLs against the asymptotic formula on the background-only Asimov set, 17 MeV/c2, normalisation (1) (x17c_toycls.C):")
w(f"  {'r':>7s} {'BR':>10s} {'CLs toys':>9s} {'err':>7s} {'toys':>5s} {'CLs asym':>9s}")
for t in tcl:
    w(f"  {t[0]:7.4f} {t[0]*2.68e-4:10.3e} {t[1]:9.4f} {t[2]:7.4f} {t[3]:5d} {t[4]:9.4f}")
a17v = next(x for x in arows if abs(x["m"] - 17.0) < 1e-6)["e1"][2]
w(f"  toy UL (weighted quadratic fit of ln CLs) = {toy_ul_r*2.68e-4:.3e}; asymptotic = {a17v:.3e}; toy/asymptotic = {toy_ul_r*2.68e-4/a17v:.3f}")
w("  -> within 10 %: the asymptotic CLs stays primary (to be re-checked on the data by X-d, Sect. 6).")
w("")
w("  Asimov expectation on the whole mass grid (normalisation (1) bands; (2), (3) medians; X17 rows at the (1) median):")
w(f"  {'m':>6s} {'R_X':>7s} {'-2s':>10s} {'-1s':>10s} {'median':>10s} {'+1s':>10s} {'+2s':>10s} {'(2) med':>10s} {'(3) med':>10s} {'rows':>7s}")
for r in arows:
    w(f"  {r['m']:6.2f} {r['RX']:7.4f} " + " ".join(f"{v:10.3e}" for v in r["e1"]) + f" {r['e2']:10.3e} {r['e3']:10.3e} {r['nX']:7.1f}")
w("")
w("  Injection test at 17 MeV/c2 (x17c_inject.txt; truth = side-band fit + X17 in normalisation (1)):")
for line in injl:
    w("    " + line)
w("  Note: injected rows = BR/2.68e-4 x N5(side bands) R_X; normalisations (2)/(3) read the same rows with a_eta N5ref resp.")
w(f"  a_P N5ref, so their slope is N5/(a N5ref) = {n5/(cfg['a0']*7458.8):.3f} resp. {n5/(cfg['aP']*7458.8):.3f} by construction (not a bias of the method).")
w("")
w("  Model choices (x17c_variants.txt; expected median at 17 MeV/c2 from each variant's own side-band fit):")
for line in var.split("\n"):
    if line.startswith("# single choices"):
        w("    single choices, refitting the nominal Asimov set (expected median at 17 MeV/c2, ratio to nominal):")
    elif line and not line.startswith("#") and not line.startswith("variant"):
        w("    " + line)
w("")
w("-" * 118)
w("10. WHAT X-d DOES (unchanged procedure)")
w("-" * 118)
w("1. Load the full data X2D: side-band cells from hists_v1.root data/S0/all/N_X2D and the masked rows 12-28 MeV/c2 from the")
w("   sealed file (same layout); CB, simulation, WP5 inputs and x17c_config.txt exactly as here (x17cl::Build, MakeConfig).")
w("2. For every mass of Sect. 7: x17cl::SetLine(I, m); for normalisations 1, 2, 3: x17cl::LimitBR(I, MakeConfig(z, norm), n, true)")
w("   -> observed UL, r_hat; the expected bands of Sect. 9 stay the FROZEN expectation (the Asimov of the observed fit is shown")
w("   in addition, labelled).")
w("3. Toy CLs at 16.7, 16.85, 17.0 (Sect. 6); FC cross-check at 17.0; local p0 on the grid, global p if Z > 3.")
w("4. Goodness of fit of the full fit (Pearson chi2/ndf; fall-back of Sect. 5); pulls per M(ee) row.")
w("5. MC: UL with the C4 and CB replicas (Sect. 4), quoted separately; 1 + cos^2 limit quoted separately.")
w("6. Counting cross-check in X17c with B_pred of Sect. 6.")
w("7. Table D1: (1) r x 2.68e-4; (2) N_X / (a_eta s_X2D); (3) relative to Y_peak (a_P); limits on N_X (rows) and on")
w("   sigma x BR x BR [nb] = BR_X BR_ee x sigma(pp -> eta X) (only with the WP8 cross section). The simulation expectation")
w("   (1.9e-5 / 5.7e-6) for comparison only.")
w("=" * 118)

os.makedirs(os.path.dirname(FROZEN), exist_ok=True)
with open(FROZEN, "w") as f:
    f.write("\n".join(L) + "\n")
print("written", FROZEN)

# ---------------------------------------------------------------- LaTeX snippets
def fmt(v):
    e = int(f"{v:.2e}".split("e")[1])
    return f"{v / 10**e:.2f}"


def tex_expected():
    rows = []
    for m in (16.7, 16.85, 17.0, 10.0, 25.0, 40.0):
        r = next((x for x in arows if abs(x["m"] - m) < 1e-6), None)
        t = trows.get((m, 1))
        t2 = trows.get((m, 2))
        t3 = trows.get((m, 3))
        if not r:
            continue
        rows.append(f"{m:.2f} & {r['e1'][2]*1e5:.2f} & [{r['e1'][1]*1e5:.2f}, {r['e1'][3]*1e5:.2f}] & [{r['e1'][0]*1e5:.2f}, {r['e1'][4]*1e5:.2f}] & "
                    + (f"{t['q'][2]*1e5:.2f} & [{t['q'][1]*1e5:.2f}, {t['q'][3]*1e5:.2f}] & [{t['q'][0]*1e5:.2f}, {t['q'][4]*1e5:.2f}] & {t['n']}" if t else "-- & -- & -- & --")
                    + f" & {r['e2']*1e5:.2f} & {r['e3']*1e5:.2f} \\\\")
    s = ["% x17c_expected.tex -- written by plots/all_report/x17c_frozen.py (X-c). Expected 90 % CLs UL x 1e5.",
         "\\begin{tabular}{r ccc cccr cc}", "\\toprule",
         " & \\multicolumn{3}{c}{(1) Asimov} & \\multicolumn{4}{c}{(1) background-only toys} & (2) & (3) \\\\",
         "$m_X$ [MeV/$c^2$] & median & $\\pm1\\sigma$ & $\\pm2\\sigma$ & median & $\\pm1\\sigma$ & $\\pm2\\sigma$ & toys & median & median \\\\", "\\midrule"]
    s += rows
    s += ["\\bottomrule", "\\end{tabular}"]
    return "\n".join(s) + "\n"


with open(f"{TAB}/x17c_expected.tex", "w") as f:
    f.write(tex_expected())


def tex_nuis():
    s = ["% x17c_nuisances.tex -- written by plots/all_report/x17c_frozen.py (X-c). Frozen priors of the X17 likelihood.",
         "\\begin{tabular}{l l l}", "\\toprule", "parameter & prior & origin \\\\", "\\midrule",
         f"$k_{{E1}}$ (CB, $M_{{ee}}<10$ MeV/$c^2$) & $1\\pm{cfg['kE1']:.2f}$ & x17ab (free $k_{{E1}}$ degenerate with $n_j$) \\\\",
         f"$k_{{E2}}$ (CB, $M_{{ee}}\\ge10$ MeV/$c^2$) & $1\\pm{cfg['kE2']:.2f}$ & D5 \\\\",
         f"$a_\\eta$ ($\\eta$-related terms) & ${cfg['a0']:.3f}\\pm{cfg['aSig']:.3f}$ & WP5 $\\times1.61$ $\\oplus$ S3 shift \\\\",
         f"$g$ ($G_\\eta$) & ${cfg['g0']:.2f}\\pm{cfg['gSig']:.2f}$ & WP5 $\\times1.61$ $\\oplus$ S3 shift \\\\",
         f"conversion tail ($\\ge12$ MeV/$c^2$) & $\\pm{100*cfg['cTail']:.1f}$\\,\\% & 100 C4 replicas \\\\",
         "X17 line width & $\\pm10$\\,\\% & spec WP9 \\\\",
         "X17 $M_{ee}$ scale & $\\pm2$\\,\\% & spec WP9 \\\\",
         f"X17 yield (1), log-normal & $\\pm{100*cfg['xSig1']:.1f}$\\,\\% & kinematics, $N_5$ bias, direct shape, peak width, BR$_\\mathrm{{dir}}$ \\\\",
         f"X17 yield (2), log-normal & $\\pm{100*cfg['xSig2']:.1f}$\\,\\% & kinematics, $a_\\eta$ \\\\",
         f"X17 yield (3), log-normal & $\\pm{100*cfg['xSig3']:.1f}$\\,\\% & kinematics, $a_P$ \\\\",
         f"variance factor & $q/{cfg['qScale']:.3f}$ & event groups (side bands) \\\\",
         "\\bottomrule", "\\end{tabular}"]
    return "\n".join(s) + "\n"


with open(f"{TAB}/x17c_nuisances.tex", "w") as f:
    f.write(tex_nuis())
print("written", f"{TAB}/x17c_expected.tex", f"{TAB}/x17c_nuisances.tex")
