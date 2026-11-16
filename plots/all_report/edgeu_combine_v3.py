#!/usr/bin/env python3
# ============================================================================
# edgeu_combine_v3.py  --  copy of edgeu_combine.py for the rulings D38 and D40 (15 November 2026):
#   D38: reads edgeu_v3/edgeu_ana_v3.json (edgeu_ana_v3.py: per-row C2 QED weights; U10 = model allowance +-0.004);
#        the uncorrected reading (no C2 factor) is taken from the original edgeu_ana.json, which D38 does not change;
#   D40: the peak-term family U8 takes the maximum of its items and the peak-tail term t = 1 +- 0.5
#        (0.5 x d ln f / dt, d ln f / dt = 0.06834 from peaktail_numbers.json; h_sim does not depend on t, so the same
#        relative term holds for every reading of the u edge).
#   Tables: tables_all/edgeu_{syst,results}_v3.{tex,txt}; json edgeu_v3/edgeu_numbers_v3.json.
# edgeu_combine.py  --  package edgeu: systematic families, conversions (yields, sigma), comparisons, tables
# ============================================================================
# Written 2026-11-15 (prefix edgeu). Reads /hdd3/.../edgeu/edgeu_ana.json, edgeu_toys.json, edgeu_bgsmooth.json; writes
#   docs/latex/tables_all/edgeu_syst.{txt,tex}     every variant: value, shift, noise, counted delta, family; family totals
#   docs/latex/tables_all/edgeu_results.{txt,tex}  f_eta(u-edge), yields, sigma (three readings) and the comparison of routes
#   docs/latex/tables_all/edgeu_toys.{txt,tex}     pseudo-data closure
#   /hdd3/.../edgeu/edgeu_numbers.json             numbers for edgeu_figs.py and the note
# Part III rules (WP7): delta^2 = max(0, Delta^2 - |stat_v^2 - stat_0^2| - MC(Delta)^2) ("literal"; stat = fit error of f),
# MC(Delta) = robust spread (half of the central 68 %) of the replica-by-replica difference; next to it the bootstrap version
# (data noise of Delta = spread of the bootstrap-by-bootstrap difference). A +-1 sigma pair counts (|D+| + |D-|)/2; family =
# maximum over its items; families in quadrature; MC and normalisation apart; rule R1: a variant the data reject (chi2 rise
# > 25 after scaling with chi2_0/ndf_0, or for other bins chi2/ndf above 1 + 5 sqrt(2/ndf) times the nominal) is listed,
# not counted. Checks (isBest == 1, no C2 factor, q = 0, T12, no peak term) are listed, not counted.
# Run: python3 -I plots/all_report/edgeu_combine.py
# ============================================================================
import json
import math
import sys

import numpy as np

sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import edgeu_core as C

OUT3 = "/hdd3/przygoda/HADES/FAT_out/all_report/edgeu_v3/"
A = json.load(open(OUT3 + "edgeu_ana_v3.json"))
A_OLD = json.load(open(C.OUT + "edgeu_ana.json"))          # uncorrected reading (unchanged by D38)
PT = json.load(open("/hdd3/przygoda/HADES/FAT_out/all_report/peaktail/peaktail_numbers.json"))
TAIL_REL = 0.5 * PT["f_u_slope"]                         # D40: t = 1 +- 0.5
TY = json.load(open(C.OUT + "edgeu_toys.json"))
BS = json.load(open(C.OUT + "edgeu_bgsmooth.json"))
NOM = A["nom"]; VR = A["var"]; CHK = A["chk"]
NOM_OLD = A_OLD["nom"]; VR_OLD = A_OLD["var"]
for k_ in ("f_unc", "stat_unc", "mc_unc", "mc68_unc"):
    NOM[k_] = NOM_OLD[k_]
D38_SANITY = (VR["c2off"]["f"], NOM_OLD["f_unc"])
f0 = NOM["f"]
out = []
def P(s=""):
    print(s); out.append(s)

# ---------------------------------------------------------------- inputs (typed, with their source)
SIG_SMASH = 1.0694      # mb, sigma(pp -> eta X) of the simulated sample (v8)
CSEC_B, FETA = 1.0104, 0.97    # C_sec of the bump window and (1 - f_eta) (WP2 / Amendment 1)
PDG26 = 1.00436         # C2, C3 rescaling to PDG 2026 (v8 V3)
N_PARTS = dict(L=6.9, period=5.0, BR=3.17, csec_up=3.1, csec_dn=4.1)   # per cent (WP7 N row, fix2, v8)
KAP_PEAK, KAP_PEAK_ERR = 0.48, 0.11
NORM_UP = math.sqrt(N_PARTS["L"] ** 2 + N_PARTS["period"] ** 2 + N_PARTS["BR"] ** 2 + N_PARTS["csec_up"] ** 2)
NORM_DN = math.sqrt(N_PARTS["L"] ** 2 + N_PARTS["period"] ** 2 + N_PARTS["BR"] ** 2 + N_PARTS["csec_dn"] ** 2)
NORM_YL = math.hypot(N_PARTS["L"], N_PARTS["period"])
# comparison values (bedge_results.md, bedge_combine.py CMP, WP6/WP7/WP8)
CMP = dict(bedge_fc=(0.42345, 0.02275, 0.0381, 0.0480), bedge_f=(0.392, 0.021, 0.035, 0.042), bedge_sigma=(0.442, 0.024, 0.040, 0.050, 0.042, 0.044),
           bedge_sigma_unc=0.409, bedge_k48=(0.569, 0.031, 0.071, 0.084), bedge_sigma_k48=0.594,   # v3: bedge_combine_v3.py (D38 + D40)
           peak=(0.429, 0.036, 0.152, 0.025, 0.041), direct=(0.522, 0.062), conv=(0.249, 0.051), tag_scale=(0.454, 0.043), tag_sigma=0.474,
           fmain_scale=0.382, fmain_sigma=0.399, b4_scale=0.312, aP26=(0.4091, 0.0348), aP=(C.AP, C.AP_STAT, C.AP_SYST, C.AP_MC),
           bedge_YB=21353, bedge_kappa=(0.51, 0.25, 0.48), bedge_aD=(0.551, 0.085, 0.231))
def to_sigma(scale):
    return scale / PDG26 * SIG_SMASH * CSEC_B * FETA

# ---------------------------------------------------------------- families
FAM = {k: v["fam"] for k, v in VR.items()}
for k in ("r089", "r110", "rpp116"):
    FAM[k] = "U2"
FAMTEXT = {"U1": "fit model (polynomial order, range, simulation weights, template evaluation)",
           "U2": "resolution: step smearing and the data/simulation resolution ratio (r = 0.89 / 1.10, own smearing per sample)",
           "U3": "M(pipi) range", "U4": "mass scales and edge function (data edge and M(pipi) scales, simulated s, M(pipi) transfer)",
           "U5": "CB (normalisation k, cells, rule, slices, like-sign shape)", "U6": "u binning", "U7": "physics weights (eta -> 3pi Dalitz plot), peak form factors, ownership",
           "U8": "peak term (a_P, its template) and its low-mass tail t = 1 +- 0.5 (D40; family maximum)", "U9": "selection (pre vs S0)",
           "U10": "C2 gamma* angle: model allowance of the per-row QED weights (+-0.004, D38)",
           "U11": "step projection of the smooth background (pseudo-data envelope)"}
NOTCOUNT = {"pol1": "R1: chi2 115.6/33 against 27.3/32 (rise 102 after scaling > 25)",
            "wtsOff": "R1 (WP7): phase-space eta -> 3pi Dalitz plot, excluded by KLOE-2 and by the F-main fit (delta chi2 38 > 25)",
            "bumpPS": "R1 (WP7): phase-space eta -> 3pi Dalitz plot, excluded by KLOE-2 and by the F-main fit (delta chi2 38 > 25)",
            "q0": "check: E(m, 0) is the kinematically wrong edge for Dalitz pairs (the M(ee) term is exact); shows what the q term does",
            "q0all": "check: E(m, 0) and every M(ee) row (the bedge row set)", "T12": "check: h_sim with the raw simulated background (bedge convention)",
            "nopeak": "check: no peak term (chi2 40.2/32)", "best": "check", "bestD": "check", "c2off": "check (the uncorrected reading)"}
def r1_reject(tag, v):
    n0, c0 = NOM["ndf"], NOM["chi2"]
    if v["r1"] == 1:
        return (v["chi2"] - c0) / (c0 / n0) - (n0 - v["ndf"]) > 25
    if v["r1"] == 2:
        return (v["chi2"] / v["ndf"]) / (c0 / n0) > 1 + 5 * math.sqrt(2.0 / v["ndf"])
    return False
def combine(which="f"):
    """rows and family values for f (C2 corrected) or f_unc"""
    rows = []
    VRx = VR if which == "f" else VR_OLD
    NOMx = NOM if which == "f" else NOM_OLD
    s0 = NOMx["statfit"] if which == "f" else NOMx["statfit"] * NOMx["f_unc"] / NOMx["f"]
    for tag, v in VRx.items():
        if which == "unc" and (v["d_unc"] is None or FAM[tag] == "U10"):
            continue
        d = v["delta"] if which == "f" else v["d_unc"]
        sv = v["statfit"] if which == "f" else v["statfit"] * (v["f_unc"] / v["f"] if v["f_unc"] else 1)
        nr = v["nr68"] if which == "f" else v["nr68_unc"]
        nb = v["nb"] if which == "f" else v["nb_unc"]
        lit2 = abs(sv ** 2 - s0 ** 2) + nr ** 2
        boot2 = nb ** 2 + nr ** 2
        dl = math.sqrt(max(0., d * d - lit2)); db = math.sqrt(max(0., d * d - boot2))
        rej = r1_reject(tag, v)
        counted = v["counted"] == 1 and not rej and tag not in NOTCOUNT
        why = NOTCOUNT.get(tag, "R1" if rej else ("check" if v["counted"] == 0 else ""))
        rows.append(dict(tag=tag, fam=FAM[tag], pair=v["pair"], f=v["f"] if which == "f" else v["f_unc"], d=d, noise_lit=math.sqrt(lit2),
                         noise_boot=math.sqrt(boot2), dl=dl, db=db, counted=counted, why=why, text=v["text"], chi2=v["chi2"], ndf=v["ndf"]))
    fams = {}
    for version in ("dl", "db"):
        fv = {}
        for fam in FAMTEXT:
            items = [r for r in rows if r["fam"] == fam and r["counted"]]
            vals = []
            pairs = {}
            for r in items:
                if r["pair"] >= 0:
                    pairs.setdefault(r["pair"], []).append(r[version])
                else:
                    vals.append(r[version])
            for p_, l in pairs.items():
                vals.append(sum(l) / len(l) if len(l) == 2 else l[0])
            fv[fam] = max(vals) if vals else 0.
        fams[version] = fv
    return rows, fams
rows_f, fams_f = combine("f")
rows_u, fams_u = combine("unc")
FAMS_U8_OLD = {k: v["U8"] for k, v in fams_f.items()}
# U11: step projection of the smooth background, envelope of the pseudo-data truths with smooth backgrounds (Asimov)
smooth_tags = ["a", "c1", "sA", "sm3", "sm5", "sm6", "om", "datalike"]
u11_rel = max(abs(TY[t]["bias_asimov"]) for t in smooth_tags) / 100.
for fv in fams_f.values():
    fv["U11"] = u11_rel * f0
for fv in fams_u.values():
    fv["U11"] = u11_rel * NOM["f_unc"]
for fv in fams_u.values():
    fv["U10"] = 0.
# D38: U10 is a model allowance (+-0.004: density along the circle, soft-lepton efficiency, posterior check), not a variant with
# data noise; like U11 it is taken as the direct shift (mean of the pair), because the literal rule would subtract the trivial
# 1/h_sim scaling of the fit error from a pure simulation scaling (literal result 0.0 %, bootstrap 0.26 %).
U10_DIRECT = 0.5 * (abs(VR["cfUp"]["delta"]) + abs(VR["cfDn"]["delta"]))
for fv in fams_f.values():
    fv["U10"] = U10_DIRECT
# D38 only (before the D40 merge), for the double ratio (its tail term is a correlated ratio term, listed)
syst_f_d38 = {k: math.sqrt(sum(x * x for x in fv.values())) for k, fv in fams_f.items()}
syst_u_d38 = {k: math.sqrt(sum(x * x for x in fv.values())) for k, fv in fams_u.items()}
# D40: U8 = max(U8 items, peak-tail term)
for fv in fams_f.values():
    fv["U8"] = max(fv["U8"], TAIL_REL * f0)
for fv in fams_u.values():
    fv["U8"] = max(fv["U8"], TAIL_REL * NOM["f_unc"])
syst_f = {k: math.sqrt(sum(x * x for x in fv.values())) for k, fv in fams_f.items()}
syst_u = {k: math.sqrt(sum(x * x for x in fv.values())) for k, fv in fams_u.items()}
SYST = syst_f["dl"]; SYST_REL = SYST / f0
SYST_U = syst_u["dl"]
# kappa = 0.48 reading: relative syst of f (+) the kappa term
lin = NOM["lin"]
def a_kap(k):
    return NOM["hd"] / (lin["C2c"] + k * lin["C3"])
E9 = 0.5 * abs(a_kap(KAP_PEAK - KAP_PEAK_ERR) - a_kap(KAP_PEAK + KAP_PEAK_ERR))
SYST_K48 = math.hypot(SYST_REL * NOM["a_k48"], E9)

P("edgeu_combine_v3.py (D38 + D40): f_eta(u-edge) = %.4f +- %.4f (stat, bootstrap) +- %.4f (syst) +- %.4f (MC rms)" % (f0, NOM["statboot"], SYST, NOM["mc"]))
P("  check: c2off variant of the v3 run %.5f = original uncorrected %.5f; peak-tail term %.4f (%.2f %%); U8 items before the merge %.4f" %
  (D38_SANITY[0], D38_SANITY[1], TAIL_REL * f0, 100 * TAIL_REL, FAMS_U8_OLD["dl"]))
P("  D38 only (U8 without the tail): syst %.4f (%.2f %%) | bootstrap %.4f (%.2f %%); uncorrected D38-only syst %.4f" %
  (syst_f_d38["dl"], 100 * syst_f_d38["dl"] / f0, syst_f_d38["db"], 100 * syst_f_d38["db"] / f0, syst_u_d38["dl"]))
P("  families (literal | bootstrap), absolute and %:")
for fam in FAMTEXT:
    P("   %-4s %.4f %.4f  (%.1f %% | %.1f %%)  %s" % (fam, fams_f["dl"][fam], fams_f["db"][fam], 100 * fams_f["dl"][fam] / f0, 100 * fams_f["db"][fam] / f0, FAMTEXT[fam]))
P("   total %.4f (%.1f %%) | bootstrap version %.4f (%.1f %%)" % (SYST, 100 * SYST_REL, syst_f["db"], 100 * syst_f["db"] / f0))
P("  uncorrected: f %.4f +- %.4f +- %.4f +- %.4f" % (NOM["f_unc"], NOM["stat_unc"], SYST_U, NOM["mc_unc"]))
P("  kappa = 0.48 reading: %.4f +- %.4f +- %.4f (E9 %.4f) +- %.4f" % (NOM["a_k48"], NOM["stat_k48"], SYST_K48, E9, NOM["mc_k48"]))

# ---------------------------------------------------------------- yields and sigma
W = A["windows"]
rep_f = np.load(OUT3 + "nom_rep.npy")
def yld(f, st, sy_rel, key, corr=True):
    w = W[key]
    base = (w["C2c"] + w["C3"]) if corr else (w["C2"] + w["C3"])
    y = f * base * C.PT3
    reps = np.array(w["rep_c"] if corr else w["rep_u"])
    mc = C.sd(rep_f * reps * C.PT3) if corr else None
    return dict(v=y, stat=y * st / f, syst=y * sy_rel, mc=mc, nb=y / C.L_DATA, base=base * C.PT3)
Y = {k: yld(f0, NOM["statboot"], SYST_REL, k) for k in ("B", "all", "B_E2", "B_E1")}
fE2 = CHK["classes"]["E2"]
Y["B_E2_own"] = dict(v=fE2["f"] * (W["B_E2"]["C2c"] + W["B_E2"]["C3"]) * C.PT3)
Y["B_E2_own"]["stat"] = Y["B_E2_own"]["v"] * fE2["statboot"] / fE2["f"]; Y["B_E2_own"]["mc"] = Y["B_E2_own"]["v"] * fE2["mc"] / fE2["f"]
SIG = dict(v=to_sigma(f0), stat=to_sigma(NOM["statboot"]), syst=to_sigma(SYST), mc=to_sigma(NOM["mc"]))
SIG["up"] = SIG["v"] * NORM_UP / 100; SIG["dn"] = SIG["v"] * NORM_DN / 100
SIG_U = dict(v=to_sigma(NOM["f_unc"]), stat=to_sigma(NOM["stat_unc"]), syst=to_sigma(SYST_U), mc=to_sigma(NOM["mc_unc"]))
SIG_K = dict(v=to_sigma(NOM["a_k48"]), stat=to_sigma(NOM["stat_k48"]), syst=to_sigma(SYST_K48), mc=to_sigma(NOM["mc_k48"]))
for s_ in (SIG_U, SIG_K):
    s_["up"] = s_["v"] * NORM_UP / 100; s_["dn"] = s_["v"] * NORM_DN / 100
P("  Y_bump(B, all M(ee)) = %.0f +- %.0f +- %.0f +- %.0f rows (%.3f nb); all M4: %.0f (%.3f nb); B E2: %.0f (own E2 fit %.0f)" %
  (Y["B"]["v"], Y["B"]["stat"], Y["B"]["syst"], Y["B"]["mc"], Y["B"]["nb"], Y["all"]["v"], Y["all"]["nb"], Y["B_E2"]["v"], Y["B_E2_own"]["v"]))
P("  sigma(pp -> eta X) = %.3f +- %.3f +- %.3f +- %.3f +%.3f/-%.3f mb (C2 corrected); uncorrected %.3f; kappa 0.48 %.3f" %
  (SIG["v"], SIG["stat"], SIG["syst"], SIG["mc"], SIG["up"], SIG["dn"], SIG_U["v"], SIG_K["v"]))
tot = math.sqrt(NOM["statboot"] ** 2 + SYST ** 2 + NOM["mc"] ** 2)
tot_b = math.sqrt(CMP["bedge_fc"][1] ** 2 + CMP["bedge_fc"][2] ** 2 + CMP["bedge_fc"][3] ** 2)
sl4 = CHK["sl4"]
DR = CMP["aP26"][0] / (f0 / PDG26)
DR_stat = DR * math.hypot(C.AP_STAT / C.AP, NOM["statboot"] / f0); DR_syst = DR * math.hypot(C.AP_SYST / C.AP, syst_f_d38["dl"] / f0); DR_mc = DR * math.hypot(C.AP_MC / C.AP, NOM["mc"] / f0)
SIG_PEAK_SYST_REL = 0.152 / 0.429      # sigma(peak) syst with the composition term (D31), before D40
DR_syst_comp = DR * math.hypot(SIG_PEAK_SYST_REL, syst_f_d38["dl"] / f0)
DR_tail = DR * 0.5 * (0.177 - PT["f_u_slope"])      # D40, listed: d ln DR / dt = 0.177 (a_P, vpt) - 0.068 (f)
P("  total (stat+syst+MC) %.4f = %.1f %% ; per-slice edge (bedge, C2 corrected) %.4f = %.1f %%" % (tot, 100 * tot / f0, tot_b, 100 * tot_b / CMP["bedge_fc"][0]))
P("  DR_u = a_P / f (PDG 2026) = %.3f +- %.3f +- %.3f +- %.3f; with the peak-composition term syst %.3f, total %.3f; peak tail (listed) +-%.3f" %
  (DR, DR_stat, DR_syst, DR_mc, DR_syst_comp, math.sqrt(DR_stat ** 2 + DR_syst_comp ** 2 + DR_mc ** 2), DR_tail))

# ---------------------------------------------------------------- tables
def w(path, lines):
    open(path, "w").write("\n".join(lines) + "\n")
L = ["# edgeu_syst_v3.txt -- edgeu_combine_v3.py (D38 per-row C2 weights, D40 peak tail in U8): variants of f_eta(u-edge) (C2 corrected); nominal %.5f, stat (fit) %.5f, stat (boot) %.5f" % (f0, NOM["statfit"], NOM["statboot"]),
     "# columns: tag family pair f delta noise_literal noise_bootstrap counted_literal counted_bootstrap status | text",
     "# literal: delta^2 - |stat_v^2 - stat_0^2| - MC(Delta)^2; bootstrap: delta^2 - boot(Delta)^2 - MC(Delta)^2; MC(Delta) = robust replica spread"]
for r in rows_f:
    L.append("%-8s %-4s %3d %.5f %+.5f %.5f %.5f %.5f %.5f %-9s | %s%s" % (r["tag"], r["fam"], r["pair"], r["f"], r["d"], r["noise_lit"], r["noise_boot"], r["dl"], r["db"],
             "counted" if r["counted"] else "listed", r["text"], ("  [" + r["why"] + "]") if r["why"] else ""))
L.append("# families (literal | bootstrap), absolute and per cent of f:")
for fam in FAMTEXT:
    L.append("FAM %-4s %.5f %.5f %5.2f %5.2f  %s" % (fam, fams_f["dl"][fam], fams_f["db"][fam], 100 * fams_f["dl"][fam] / f0, 100 * fams_f["db"][fam] / f0, FAMTEXT[fam]))
L.append("TOTAL %.5f %.5f %5.2f %5.2f" % (SYST, syst_f["db"], 100 * SYST_REL, 100 * syst_f["db"] / f0))
L.append("# U11 = max |Asimov bias| of the smooth-background truths %s = %.2f %%" % (",".join(smooth_tags), 100 * u11_rel))
L.append("# uncorrected reading: families %s ; total %.5f" % (" ".join("%s %.4f" % (k, v) for k, v in fams_u["dl"].items()), SYST_U))
w(C.TAB + "edgeu_syst_v3.txt", L)
T = ["% edgeu_syst_v3.tex -- edgeu_combine_v3.py (rulings D38, D40): systematic families of f_eta(u-edge), per cent",
     "\\begin{tabular}{llrr}", "\\hline", " family & content & literal [\\%] & bootstrap [\\%] \\\\", "\\hline"]
for fam in FAMTEXT:
    T.append(" %s & %s & %.1f & %.1f \\\\" % (fam, FAMTEXT[fam].replace("_", "\\_").replace("->", "$\\to$").replace("+-", "$\\pm$"), 100 * fams_f["dl"][fam] / f0, 100 * fams_f["db"][fam] / f0))
T += ["\\hline", " total & & %.1f & %.1f \\\\" % (100 * SYST_REL, 100 * syst_f["db"] / f0), "\\hline", "\\end{tabular}"]
w(C.TAB + "edgeu_syst_v3.tex", T)
def pm(v, *e, d=3):
    return ("%%.%df" % d) % v + "".join(("$\\pm$%%.%df" % d) % x for x in e)
R = ["# edgeu_results_v3.txt -- edgeu_combine_v3.py (D38 + D40); S0, isBest >= 0, L = 6079.2 nb^-1; yields in counts (stored combinations) x 64/63; errors stat (bootstrap), syst, MC (replica rms)",
     "f_eta(u-edge) C2 per-row (D38)   %.4f %.4f %.4f %.4f   total %.4f (%.1f %%)" % (f0, NOM["statboot"], SYST, NOM["mc"], tot, 100 * tot / f0),
     "f_eta(u-edge) uncorrected          %.4f %.4f %.4f %.4f" % (NOM["f_unc"], NOM["stat_unc"], SYST_U, NOM["mc_unc"]),
     "kappa = 0.48 reading (C2 corr.)    %.4f %.4f %.4f %.4f   (E9 %.4f)" % (NOM["a_k48"], NOM["stat_k48"], SYST_K48, NOM["mc_k48"], E9),
     "four-slice estimator in u          %.4f %.4f  -  %.4f" % (sl4["f"], sl4["statboot"], sl4["mc"]),
     "Y_bump B all M(ee) [counts]          %.0f %.0f %.0f %.0f   %.3f nb   (SMASH c C2 + C3: %.1f)" % (Y["B"]["v"], Y["B"]["stat"], Y["B"]["syst"], Y["B"]["mc"], Y["B"]["nb"], Y["B"]["base"]),
     "Y_bump all M4 [counts]              %.0f %.0f %.0f %.0f   %.3f nb   (SMASH %.1f)" % (Y["all"]["v"], Y["all"]["stat"], Y["all"]["syst"], Y["all"]["mc"], Y["all"]["nb"], Y["all"]["base"]),
     "Y_bump B E2 [counts] (global f)      %.0f %.0f %.0f %.0f ; own E2 fit %.0f %.0f - %.0f" % (Y["B_E2"]["v"], Y["B_E2"]["stat"], Y["B_E2"]["syst"], Y["B_E2"]["mc"], Y["B_E2_own"]["v"], Y["B_E2_own"]["stat"], Y["B_E2_own"]["mc"]),
     "sigma C2 corrected [mb]            %.4f %.4f %.4f %.4f +%.4f -%.4f" % (SIG["v"], SIG["stat"], SIG["syst"], SIG["mc"], SIG["up"], SIG["dn"]),
     "sigma uncorrected [mb]             %.4f %.4f %.4f %.4f +%.4f -%.4f" % (SIG_U["v"], SIG_U["stat"], SIG_U["syst"], SIG_U["mc"], SIG_U["up"], SIG_U["dn"]),
     "sigma kappa 0.48 [mb]              %.4f %.4f %.4f %.4f +%.4f -%.4f" % (SIG_K["v"], SIG_K["stat"], SIG_K["syst"], SIG_K["mc"], SIG_K["up"], SIG_K["dn"]),
     "DR_u = a_P/f (PDG 2026)            %.3f %.3f %.3f %.3f" % (DR, DR_stat, DR_syst, DR_mc),
     "E1/E2 split: a_D %.3f stat %.3f MC %.3f ; kappa %.3f stat %.3f MC %.3f (robust %.3f)" % (CHK["kappa_free"]["a"], CHK["kappa_free"]["a_stat"], CHK["kappa_free"]["a_mc"],
                                                                                            CHK["kappa_free"]["kappa"], CHK["kappa_free"]["kappa_stat"], CHK["kappa_free"]["kappa_mc"], CHK["kappa_free"]["kappa_mc68"]),
     "# comparison (sigma in mb): per-slice edge (bedge_v3, C2 per-row) 0.442; uncorr. 0.409; kappa 0.48 0.594; peak side bands 0.429; direct 0.522; conversion 0.249;",
     "#   B4 (not determined) %.3f; F-main 0.399; ECAL tag (verification) 0.474; SMASH 1.069" % to_sigma(CMP["b4_scale"])]
w(C.TAB + "edgeu_results_v3.txt", R)
RT = ["% edgeu_results_v3.tex -- edgeu_combine_v3.py (rulings D38: per-row C2 weights; D40: peak tail in U8). Counts = stored pi+pi-e+e- combinations.", "\\begin{tabular}{llll}", "\\hline",
      " quantity & uncorrected & C2 corrected (per-row, D38) & C2 corr., $\\kappa=0.48$ \\\\", "\\hline",
      " $f_\\eta$(u edge) & $%s$ & $%s$ & $%s$ \\\\" % (pm(NOM["f_unc"], NOM["stat_unc"], SYST_U, NOM["mc_unc"]).replace("$", ""), pm(f0, NOM["statboot"], SYST, NOM["mc"]).replace("$", ""),
                                                   pm(NOM["a_k48"], NOM["stat_k48"], SYST_K48, NOM["mc_k48"]).replace("$", "")),
      " $Y_{\\mathrm{bump}}$(B, all $M_{ee}$) [counts] & -- & $%.0f\\pm%.0f\\pm%.0f\\pm%.0f$ & -- \\\\" % (Y["B"]["v"], Y["B"]["stat"], Y["B"]["syst"], Y["B"]["mc"]),
      " $Y/L$ (B) [visible counts per nb$^{-1}$] & -- & $%.2f\\pm%.2f\\pm%.2f\\pm%.2f$ & -- \\\\" % (Y["B"]["nb"], Y["B"]["stat"] / C.L_DATA, Y["B"]["syst"] / C.L_DATA, Y["B"]["mc"] / C.L_DATA),
      " $\\sigma(pp\\to\\eta X)$ [mb] & $%.3f\\pm%.3f\\pm%.3f\\pm%.3f^{+%.3f}_{-%.3f}$ & $%.3f\\pm%.3f\\pm%.3f\\pm%.3f^{+%.3f}_{-%.3f}$ & $%.3f\\pm%.3f\\pm%.3f\\pm%.3f^{+%.3f}_{-%.3f}$ \\\\" %
      (SIG_U["v"], SIG_U["stat"], SIG_U["syst"], SIG_U["mc"], SIG_U["up"], SIG_U["dn"], SIG["v"], SIG["stat"], SIG["syst"], SIG["mc"], SIG["up"], SIG["dn"],
       SIG_K["v"], SIG_K["stat"], SIG_K["syst"], SIG_K["mc"], SIG_K["up"], SIG_K["dn"]),
      " DR$_u=a_P/f_\\eta$(u edge), with the peak-composition term & -- & $%.2f\\pm%.2f\\pm%.2f\\pm%.2f$ (total $\\pm%.2f$; peak tail $\\pm%.2f$ listed) & -- \\\\" %
      (DR, DR_stat, DR_syst_comp, DR_mc, math.sqrt(DR_stat ** 2 + DR_syst_comp ** 2 + DR_mc ** 2), DR_tail),
      "\\hline", "\\end{tabular}"]
w(C.TAB + "edgeu_results_v3.tex", RT)
TT = ["# edgeu_toys.txt -- edgeu_combine.py (prefix edgeu): closure of f_eta(u-edge) on pseudo-data (edgeu_toys.py)",
      "# tag truth Asimov_bias[%] AA_bias[%] AA_pullSD AA_cov RF_bias[%] RF_pullSD RF_cov shift_vs_a[%] | truth"]
for k, o in TY.items():
    if k == "asimov_free_sigma":
        continue
    TT.append("%-8s %.4f %+7.2f %+7.2f %.2f %.2f %+7.2f %.2f %.2f %s | %s" % (k, o["ftrue"], o["bias_asimov"], o["bias_aa"], o["pull_sd_aa"], o["cov_aa"], o["bias_rf"], o["pull_sd_rf"], o["cov_rf"],
                                                                         ("%+6.2f" % o["shift_vs_a_asimov"]) if "shift_vs_a_asimov" in o else "   -  ", o["text"]))
TT.append("# free smearing on the Asimov of (a): %.2f +- %.2f MeV (data %.2f +- %.2f MeV)" % (1e3 * TY["asimov_free_sigma"]["sigma"], 1e3 * TY["asimov_free_sigma"]["err"],
                                                                                         1e3 * A["free_sigma"]["data"]["sigma"], 1e3 * A["free_sigma"]["data"]["err"]))
pass  # toys table unchanged (not rewritten)
TX = ["% edgeu_toys.tex -- edgeu_combine.py (prefix edgeu): pseudo-data closure of f_eta(u-edge); 1000 toys per truth", "\\begin{tabular}{lrrrrr}", "\\hline",
      " truth & bias (Asimov) [\\%] & shift vs (a) [\\%] & pull SD (AA) & bias RF [\\%] & pull SD (RF) \\\\", "\\hline"]
for k in ["a", "c1", "a1", "q-", "q--", "q+", "om", "ty", "sA", "sm3", "sm5", "sm6", "raw", "mass", "wide", "narrow", "kap", "pk30", "pk55", "datalike"]:
    o = TY[k]
    TX.append(" %s & %+.1f & %s & %.2f & %+.1f & %.2f \\\\" % (k.replace("_", "\\_"), o["bias_asimov"], ("%+.1f" % o["shift_vs_a_asimov"]) if "shift_vs_a_asimov" in o else "--", o["pull_sd_aa"], o["bias_rf"], o["pull_sd_rf"]))
TX += ["\\hline", "\\end{tabular}"]
pass  # toys table unchanged (not rewritten)
NUM = dict(f=f0, stat=NOM["statboot"], statfit=NOM["statfit"], syst=SYST, syst_boot=syst_f["db"], mc=NOM["mc"], mc68=NOM["mc68"], total=tot, total_rel=tot / f0,
           fams=fams_f, fams_unc=fams_u, f_unc=NOM["f_unc"], stat_unc=NOM["stat_unc"], syst_unc=SYST_U, mc_unc=NOM["mc_unc"],
           a_k48=NOM["a_k48"], stat_k48=NOM["stat_k48"], syst_k48=SYST_K48, mc_k48=NOM["mc_k48"], E9=E9, u11_rel=u11_rel,
           Y=Y, sigma=SIG, sigma_unc=SIG_U, sigma_k48=SIG_K, DR=[DR, DR_stat, DR_syst, DR_mc], sl4=sl4, cmp=CMP, b4_sigma=to_sigma(CMP["b4_scale"]),
           norm_up=NORM_UP, norm_dn=NORM_DN, rows=rows_f, tail_rel=TAIL_REL, syst_d38=syst_f_d38, syst_unc_d38=syst_u_d38,
           DR_syst_comp=DR_syst_comp, DR_tail=DR_tail, sig_d38_syst=to_sigma(syst_f_d38["dl"]))
json.dump(NUM, open(OUT3 + "edgeu_numbers_v3.json", "w"), indent=1, default=float)
open(OUT3 + "edgeu_combine_v3_stdout.txt", "w").write("\n".join(out) + "\n")
