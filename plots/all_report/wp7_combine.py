#!/usr/bin/env python3
# ============================================================================
# wp7_combine.py -- WP7 systematics (SPEC report_all Sect. 5), prefix wp7: combination, Table S1, Fig. S1, tests
# ============================================================================
# Written by WP7, 2026-11-14. Inputs (read-only):
#   /hdd3/przygoda/HADES/FAT_out/all_report/wp7/syst_raw.txt, periods_raw.txt   (wp7_syst.C)
#   numbers of other packages typed below with their source (WP4 verdict/summary, WP5/WP6 results.txt, WP2, x17d).
# Rules (spec Sect. 5):
#   delta^2 = max(0, Delta^2 - |sigma_v^2 - sigma_0^2|)  (sigma = statistical error of the variant and of its reference);
#   per family: maximum over the items; an item is a discrete alternative or a +-1 sigma pair, (|D+| + |D-|)/2;
#   families in quadrature; MC and N listed apart; S11: a variation counts only if |Delta| > 2 sqrt(|sigma_v^2 - sigma_0^2|),
#   and only for efficiency-corrected quantities.
# WP7 decisions (flagged in the results note):
#   (R1) a variant with the same fit bins is "rejected by the data" (listed, not counted) when its data-only Pearson chi2
#        exceeds the nominal one by more than 25 after scaling with the nominal chi2/ndf: (chi2_v - chi2_0)/(chi2_0/ndf_0)
#        - (ndf_0 - ndf_v) > 25. The fit used: F-main for the F-main quantities, F-gold for BR_G, stage A for the
#        side-band quantities (a variant is rejected for a quantity if any fit that feeds it fails: F-main for all, F-gold for
#        BR_G, stage A for the side-band quantities). Variants with other bins (CR, range, binning): rejected if
#        (chi2_v/ndf_v)/(chi2_0/ndf_0) > 1 + 5 sqrt(2/ndf_v). Not applied to raw BB-lite (variance differs) and S11 (other data).
#   (R2) within a family the maximum is the nominal combination (spec); the quadrature of independent sub-sources is
#        given as an alternative in syst.txt.
# Run: python3 -I plots/all_report/wp7_combine.py
# ============================================================================
import math, os, sys
import numpy as np

OUT = "/hdd3/przygoda/HADES/FAT_out/all_report/wp7"
TAB = "/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all"
FIG = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_all"

# ---------------------------------------------------------------------------------------------- read the raw results
cfg, Q, MCQ = {}, {}, {}
for line in open(os.path.join(OUT, "syst_raw.txt")):
    if line.startswith("#"):
        continue
    w = line.split()
    if not w:
        continue
    if w[0] == "cfg":
        if len(w) > 3 and w[3] == "FAILED":
            cfg[w[1]] = dict(fam=w[2], failed=True)
            continue
        txt = line.split("|", 1)[1].strip() if "|" in line else ""
        cfg[w[1]] = dict(fam=w[2], failed=False, chiF=float(w[4]), ndfF=int(w[5]), chiA=float(w[7]), ndfA=int(w[8]),
                         chiG=float(w[10]), ndfG=int(w[11]), text=txt)
    elif w[0] == "q":
        Q.setdefault(w[1], {})[w[2]] = (float(w[3]), float(w[4]))
        if len(w) > 14 and w[13] == "mc":
            MCQ.setdefault(w[1], {})[w[2]] = float(w[14])

REP = {}
for line in open(os.path.join(OUT, "syst_rep.txt")):
    if line.startswith("rep "):
        w = line.split()
        REP.setdefault(w[1], {})[w[2]] = np.array([float(x) for x in w[3:]])

def mcdiff(tag, q):
    """robust MC noise of the difference variant - reference: half of the central 68 % interval of the per-replica
    differences (same replica index = same simulation replica); None if not available"""
    r = ref(tag)
    a, b = REP.get(tag, {}).get(q), REP.get(r, {}).get(q)
    if a is None or b is None:
        return None
    d = a - b
    d = d[np.isfinite(d)]
    if len(d) < 20:
        return None
    lo, hi = np.percentile(d, [15.87, 84.13])
    return 0.5 * (hi - lo)

def val(tag, q):
    return Q.get(tag, {}).get(q, (float("nan"), float("nan")))

# ---------------------------------------------------------------------------------------------- quantities
# key, label, unit, efficiency-corrected (S11 applies), fit (for R1), estimator status
QUANT = [
    ("Y_bump",    "Y_bump",    "rows", False, "A", "B4, all M(ee) (E2 extrapolated); MAIN, EXTRAPOLATED"),
    ("Y_peak",    "Y_peak",    "rows", False, "A", "peak side bands; MAIN"),
    ("Y5",        "Y5",        "rows", False, "F", "F-main C5 in P; CROSS-CHECK"),
    ("Y4",        "Y4",        "rows", False, "F", "F-main C4 in P; CROSS-CHECK"),
    ("a_eta",     "a_eta",     "",     True,  "F", "F-main; CROSS-CHECK"),
    ("BR_G",      "BR_G",      "",     True,  "G", "F-gold; CROSS-CHECK"),
    ("BR_A",      "BR_A",      "",     True,  "F", "2.68e-4 rho_d (F-main); CROSS-CHECK"),
    ("kappa",     "kappa",     "",     True,  "F", "F-main (constrained 1 +- 0.2); CROSS-CHECK"),
    ("DR",        "DR",        "",     True,  "A", "side bands / B4 (bump extrapolated); MAIN"),
    ("sigma_eta", "sigma_eta", "mb",   True,  "A", "a_P x 1.0665 mb x C_sec (1 - f_eta), a_P = peak side bands / SMASH; MAIN"),
]
EXTRA = [
    ("Y_bump_E2",   "rows", False, "A", "B4, M(ee) >= 10 MeV; MAIN"),
    ("a_P",         "",     True,  "A", "peak side bands / SMASH (C4 + C5 in P)"),
    ("a_B",         "",     True,  "A", "B4 (E2) / SMASH (C2 + C3, E2, B)"),
    ("DR_F",        "",     True,  "F", "F-main double ratio; CROSS-CHECK"),
    ("sigma_eta_F", "mb",   True,  "F", "a_eta (F-main) x 1.0665 mb x C_sec (1 - f_eta); CROSS-CHECK"),
    ("Y_bump_F",    "rows", False, "F", "F-main Y2 + Y3 in B; CROSS-CHECK"),
    ("Y_peak_F",    "rows", False, "F", "F-main Y4 + Y5 in P; CROSS-CHECK"),
    ("g",           "",     True,  "F", "G_eta factor (F-main)"),
    ("rho_d",       "",     True,  "F", "direct-decay scale (F-main)"),
]
ALLQ = [(q[0], q[2], q[3], q[4], q[5]) for q in QUANT] + EXTRA
QINFO = {q[0]: q for q in ALLQ}

# G3-built variants are compared with the G3-built nominal (identical to the nominal to 1e-7, see syst_raw.txt)
G3 = {"S2_Uhi750", "S2_CRK044", "S2_CRK0425", "S2_K0veto", "S3_omega085", "S5_r4075", "S5_bin25", "S5_sliceUp", "S5_sliceDn",
      "S5_mee08", "S5_mee12", "S6_pipiUp", "S6_pipiDn"}
FIXF = {"S2_U560", "S2_U565", "S2_Uhi650", "S2_Uhi750", "S2_CRK043", "S2_CRK044", "S2_CRK0425", "S2_K0veto"}   # F-main kept nominal
SAMEBINS = {"S1_kE2fix1", "S1_k1032", "S1_LSsum", "S1_CB20", "S1_geo", "S3_quad", "S3_tiltj", "S3_expb", "S3_commonC", "S3_omega",
            "S3_kBtie", "S3_D11data", "S3_D11dataOm", "S4_pol4", "S4_ghostA", "S4_modeB", "S4_wdb", "S6_sUp", "S6_sDn", "S6_sK0",
            "S6_r089", "S6_r110", "S6_pipiUp", "S6_pipiDn", "S7_V2", "S7_bump1aY", "S7_bumpPS", "S7_convOff", "S7_GS",
            "S8_g015", "S8_g050", "S8_gfree", "S9_fact"}

def ref(tag):
    return "nomG3" if tag in G3 else "nom"

DIFFBINS = {"S2_U560", "S2_U565", "S2_Uhi650", "S2_Uhi750", "S2_K0veto", "S3_omega085", "S5_r4265", "S5_r4570", "S5_r4075", "S5_bin25",
            "S5_sliceUp", "S5_sliceDn", "S5_mee08", "S5_mee12"}
FEEDS = {"F": ("F",), "G": ("F", "G"), "A": ("F", "A")}   # fits that feed a quantity

def fitcheck(tag, fit):
    """(R1) returns (rejected, z or ratio) for one fit of a variant; None if not applicable."""
    if tag not in cfg or (tag not in SAMEBINS and tag not in DIFFBINS):
        return None
    if fit == "F" and tag in FIXF:
        return None
    key = {"F": ("chiF", "ndfF"), "A": ("chiA", "ndfA"), "G": ("chiG", "ndfG")}[fit]
    c, c0 = cfg[tag], cfg[ref(tag)]
    if c0[key[1]] <= 0 or c[key[1]] <= 0:
        return None
    s2 = c0[key[0]] / c0[key[1]]
    if tag in SAMEBINS:
        z = (c[key[0]] - c0[key[0]]) / s2 - (c0[key[1]] - c[key[1]])
        return (z > 25., "%s dchi2/s2 %.0f" % (fit, z))
    r = (c[key[0]] / c[key[1]]) / s2
    lim = 1. + 5. * math.sqrt(2. / c[key[1]])
    return (r > lim, "%s chi2/ndf ratio %.2f (limit %.2f)" % (fit, r, lim))

def rejected(tag, q):
    """(R1) data rejection: any fit that feeds the quantity fails; returns (flag, text)."""
    out, txt = False, []
    for fit in FEEDS[QINFO[q][3]]:
        r = fitcheck(tag, fit)
        if r is None:
            continue
        out = out or r[0]
        if r[0]:
            txt.append(r[1])
    return out, ", ".join(txt)

def mcof(tag, q):
    """replica MC error of a configuration; the nominal one where the chain was not run (other simulation samples, G3)"""
    m = MCQ.get(tag, {}).get(q, float("nan"))
    if m != m:
        m = MCQ.get("nom", {}).get(q, float("nan"))
    return m if m == m else 0.

def item(tag, q):
    """Delta, Barlow delta (sigma = stat (+) MC; and stat only), sigma_v, sigma_0, flags for one variant."""
    v, sv = val(tag, q)
    r = ref(tag)
    v0, s0 = val(r, q)
    mv, m0 = mcof(tag, q), mcof(r, q)
    if (v != v) and tag in FIXF:       # F-main / F-gold quantities are unchanged by construction
        v, sv, mv = v0, s0, m0
    D = v - v0
    dstat2 = abs(sv * sv - s0 * s0) if (sv == sv and s0 == s0) else 0.
    md = mcdiff(tag, q)
    dmc2 = md * md if md is not None else 0.
    dtot2 = dstat2 + dmc2
    d = math.sqrt(max(0., D * D - dtot2)) if D == D else float("nan")
    ds = math.sqrt(max(0., D * D - dstat2)) if D == D else float("nan")
    rej, z = rejected(tag, q)
    return dict(tag=tag, v=v, sv=sv, v0=v0, s0=s0, mv=mv, m0=m0, md=md, D=D, d=d, d_stat=ds, dstat=math.sqrt(dstat2), dtot=math.sqrt(dtot2), rej=rej, z=z)

# ---------------------------------------------------------------------------------------------- external inputs
nomv = {q: val("nom", q)[0] for q in QINFO}
noms = {q: val("nom", q)[1] for q in QINFO}
# WP6 side-band width variants (wp6/results.txt; same estimator, nominal templates): value, stat
SBW = {"narrow 5.6 MeV": (4732.71, 510.391, 300.845), "wide 11.6 MeV": (4745.16, 387.984, 224.546)}
PSB0 = (5087.97, 431.626, 245.122)
PSB_DERIVED = {"Y_peak": 1., "a_P": 1., "sigma_eta": 1., "DR": 1.}   # quantities proportional to the side-band yield
# S10 closure (relative): WP4 RF scenario (a) biases (wp4_summary.txt) and the WP6 closure terms
b4w, psbw = 0.02552, 0.00524
sig_tr = 239.0 / 5087.97          # WP6 displaced windows (4.7 %)
crl = 0.375                       # WP6 CR-L S10 term, relative to Y_bump(E2)
dr_bias = (1 - psbw) / (1 - b4w) - 1
S10 = {"Y_bump": math.hypot(b4w, crl), "Y_bump_E2": math.hypot(b4w, crl), "a_B": math.hypot(b4w, crl),
       "Y_peak": math.hypot(psbw, sig_tr), "a_P": math.hypot(psbw, sig_tr), "sigma_eta": math.hypot(psbw, sig_tr),
       "DR": math.sqrt(dr_bias ** 2 + sig_tr ** 2 + crl ** 2),
       "Y5": 0.0157, "Y4": 0.1208, "a_eta": 0.1006, "BR_G": 0.1482, "BR_A": 0.0969, "kappa": 0.0178,
       "DR_F": (1 - 0.0658) / (1 - 0.1158) - 1, "sigma_eta_F": 0.1006, "Y_bump_F": 0.1158, "Y_peak_F": 0.0658, "g": 0.3275,
       "rho_d": 0.0969}
S10TXT = {"Y_bump": "CR-L 37.5 % (+) WP4 B4 2.6 %", "Y_bump_E2": "CR-L 37.5 % (+) WP4 B4 2.6 %", "a_B": "as Y_bump(E2)",
          "Y_peak": "sigma_tr 239 rows 4.7 % (+) WP4 0.5 %", "a_P": "as Y_peak", "sigma_eta": "as Y_peak",
          "DR": "sigma_tr (+) CR-L (+) WP4 2.1 %"}
# S12 efficiency model (relative): pion vertex 2.5 % (absolute scales), 2.0 % peak/bump difference of the data/sim
# pion-vertex efficiency (0.917/0.941 peak vs 0.945/0.950 tagged; WP1 Table A2) for the ratios
S12 = {"a_eta": 0.025, "a_P": 0.025, "a_B": 0.025, "sigma_eta": 0.025, "sigma_eta_F": 0.025,
       "DR": 0.020, "DR_F": 0.020, "BR_G": 0.020, "BR_A": 0.020, "kappa": 0.020, "rho_d": 0.020}
# MC (nominal replica RMS; WP5 results.txt derived FMAIN / FGOLD mc68 for BR_G; WP6 results.txt est)
MC = {"Y_bump": 4109.78, "Y_bump_E2": 1919.6, "Y_peak": 245.122, "a_P": 0.0248145, "a_B": 0.0763396, "DR": 0.533303,
      "sigma_eta": 0.0248145 / 0.418249 * nomv["sigma_eta"], "Y5": 570.66, "Y4": 495.93, "a_eta": 0.039629,
      "BR_G": 6.6354e-05, "BR_A": 6.0679e-05, "kappa": 0.26701, "DR_F": 0.13367,
      "sigma_eta_F": 0.039629 / 0.38228397 * nomv["sigma_eta_F"], "Y_bump_F": 2274.1, "Y_peak_F": 719.8, "g": 0.36576,
      "rho_d": 0.22642}
# N normalisation (relative, absolute quantities only): L 6.9 %; PDG BRs (eta -> pi+pi-pi0 22.92 +- 0.28 %, pi0 Dalitz
# 1.174 +- 0.035 %: 3.2 % for the bump channels; eta -> pi+pi-gamma 4.22 +- 0.08 %, eta -> pi+pi-e+e- 2.68 +- 0.11 e-4 with the
# simulated peak mix C4 : C5 = 0.48 : 0.52 in P: 2.3 %); C_sec (1 - f_eta) = 0.98 +0.03/-0.04 (+3.1/-4.1 %); live time +7.3/-0 % separate
L_ = 0.069; pdgB = math.hypot(0.0122, 0.0298); pdgP = math.hypot(0.477 * 0.019, 0.523 * 0.041)
Nn = {"a_eta": (math.hypot(L_, pdgB), math.hypot(L_, pdgB)), "a_P": (math.hypot(L_, pdgP),) * 2, "a_B": (math.hypot(L_, pdgB),) * 2,
      "sigma_eta": (math.sqrt(L_ ** 2 + pdgP ** 2 + 0.031 ** 2), math.sqrt(L_ ** 2 + pdgP ** 2 + 0.041 ** 2)),
      "sigma_eta_F": (math.sqrt(L_ ** 2 + pdgB ** 2 + 0.031 ** 2), math.sqrt(L_ ** 2 + pdgB ** 2 + 0.041 ** 2)),
      "BR_G": (pdgB, pdgB), "BR_A": (pdgB, pdgB), "DR": (math.hypot(pdgB, pdgP),) * 2, "DR_F": (math.hypot(pdgB, pdgP),) * 2}
NTXT = {"a_eta": "L 6.9 % (+) PDG 3.2 %", "a_P": "L (+) PDG 2.3 %", "a_B": "L (+) PDG 3.2 %",
        "sigma_eta": "L (+) PDG 2.3 % (+) C_sec(1-f_eta) +3.1/-4.1 %; live time +7.3/-0 % separate",
        "sigma_eta_F": "L (+) PDG 3.2 % (+) C_sec(1-f_eta) +3.1/-4.1 %; live time separate",
        "BR_G": "PDG BR(3pi) BR(pi0 Dalitz) 3.2 %", "BR_A": "PDG BR(3pi) BR(pi0 Dalitz) 3.2 %",
        "DR": "PDG peak/bump BR ratio 3.9 % (comparison with 1)", "DR_F": "PDG peak/bump BR ratio 3.9 %"}

# ---------------------------------------------------------------------------------------------- families
FAM = [
    ("S1", "CB", [["S1_kE2fix1"], ["S1_k1032"], ["S1_LSsum"], ["S1_CB20"], ["S1_geo"]], ["S1_E1kept"]),
    ("S2", "CR definition", [["S2_U560"], ["S2_U565"], ["S2_Uhi650"], ["S2_Uhi750"], ["S2_K0veto"], ["SBW:narrow 5.6 MeV"], ["SBW:wide 11.6 MeV"]],
     ["S2_CRK043", "S2_CRK044", "S2_CRK0425"]),
    ("S3", "background model", [["S3_quad"], ["S3_tiltj"], ["S3_expb"], ["S3_commonC"], ["S3_omega"], ["S3_omega085"], ["S3_kBtie"],
                                ["S3_D11data"], ["S3_D11dataOm"]], []),
    ("S4", "background template", [["S4_pol4"], ["S4_raw"], ["S4_ghostA"], ["S4_modeB"], ["S4_wdb"], ["S4_GLShyb"]], []),
    ("S5", "range, binning", [["S5_r4265"], ["S5_r4570"], ["S5_r4075"], ["S5_bin25"], ["S5_sliceUp", "S5_sliceDn"], ["S5_mee08", "S5_mee12"]], []),
    ("S6", "mass transfer", [["S6_sUp", "S6_sDn"], ["S6_sK0"], ["SFREE"], ["S6_r089", "S6_r110"], ["S6_pipiUp", "S6_pipiDn"]], []),
    ("S7", "physics weights", [["S7_V2"], ["S7_bump1aY"], ["S7_bumpPS"], ["S7_convOff"], ["S7_GS"]], ["S7_wtsOff"]),
    ("S8", "G_eta", [["S8_g015"], ["S8_g050"], ["S8_gfree"]], []),
    ("S9", "signal templates", [["S9_fact"]], []),
    ("S10", "closure", [], []),
    ("S11", "selection", [["S11_pre"], ["S11_dz40"], ["S11_vz150", "S11_vz170"]], ["S11_PEAK"]),
    ("S12", "efficiency model", [], []),
]
SUBSRC = {"S6": [["S6_sUp", "S6_sK0", "SFREE"], ["S6_r089"], ["S6_pipiUp"]],          # first element of each item group
          "S5": [["S5_r4265", "S5_r4570", "S5_r4075"], ["S5_bin25"], ["S5_sliceUp"], ["S5_mee08"]],
          "S7": [["S7_V2"], ["S7_bump1aY", "S7_bumpPS"], ["S7_convOff", "S7_GS"]]}

def eval_item(group, q, fam):
    """returns dict with delta (counted), details; group = list of 1 (discrete) or 2 (pair) tags"""
    if group[0].startswith("SBW:"):
        if q not in PSB_DERIVED:
            return dict(name=group[0], d=0., D=0., na=True, rej=False, det="n/a (side-band width: peak quantities only)")
        v, sv, mv = SBW[group[0][4:]]
        D = (v - PSB0[0]) / PSB0[0]
        dds = abs(sv ** 2 - PSB0[1] ** 2) / PSB0[0] ** 2
        dd = dds   # no per-replica difference available (WP6 numbers): stat-only Barlow (conservative)
        d = math.sqrt(max(0., D * D - dd)); ds = math.sqrt(max(0., D * D - dds))
        return dict(name="side bands " + group[0][4:], d=d * abs(nomv[q]), d_stat=ds * abs(nomv[q]), D=D * nomv[q], na=False, rej=False,
                    det="WP6 value %.1f +- %.1f (stat) +- %.1f (MC) (rel. %+.2f %%)" % (v, sv, mv, 100 * D))
    if group[0] == "SFREE":
        # free s (WP5): whole fit s_hat = 0.99653 (0.42 of +0.0009), peak region 0.99500 (1.28 of -0.0009): linear interpolation
        a, b = item("S6_sUp", q), item("S6_sDn", q)
        d = max(0.42 * a["d"], 1.28 * b["d"]); ds = max(0.42 * a["d_stat"], 1.28 * b["d_stat"])
        return dict(name="s free (interpolated)", d=d, d_stat=ds, D=d, na=False, rej=False, det="max(0.42 d(sUp), 1.28 d(sDn))")
    its = [item(t, q) for t in group]
    rej = any(i["rej"] for i in its)
    if fam == "S11":
        if not QINFO[q][2]:
            return dict(name="/".join(group), d=0., D=its[0]["D"], na=True, rej=False, det="n/a (raw yield, not efficiency corrected)")
        ds, dss = [], []
        for i in its:
            cnt = abs(i["D"]) > 2 * i["dtot"]
            i["count"] = cnt
            ds.append(i["d"] if cnt else 0.)
            dss.append(i["d_stat"] if abs(i["D"]) > 2 * i["dstat"] else 0.)
        d = sum(ds) / len(ds) if len(its) == 2 else ds[0]
        d_s = sum(dss) / len(dss) if len(its) == 2 else dss[0]
    else:
        d = sum(i["d"] for i in its) / len(its)
        d_s = sum(i["d_stat"] for i in its) / len(its)
    det = "; ".join("%s %.5g->%.5g (D %+.4g, stat %.3g/%.3g, MC %.3g/%.3g, MC of D %s, d %.4g, d(stat-only Barlow) %.4g%s%s)" % (i["tag"], i["v0"], i["v"], i["D"], i["sv"], i["s0"], i["mv"], i["m0"], ("%.3g" % i["md"]) if i["md"] is not None else "n/a", i["d"], i["d_stat"],
                    (", REJECTED (R1: %s)" % i["z"]) if i["rej"] else "", (", not significant" if fam == "S11" and not i.get("count", True) else ""))
                    for i in its)
    return dict(name="/".join(group), d=0. if rej else d, d_norej=d, d_stat=0. if rej else d_s, D=its[0]["D"], na=False, rej=rej, det=det)

RES = {}
for q, unit, effc, fit, status in ALLQ:
    v0, s0 = nomv[q], noms[q]
    fams = {}
    for fam, title, items, info in FAM:
        rows = []
        if fam == "S10":
            fams[fam] = dict(d=S10.get(q, 0.) * abs(v0), rows=[dict(name="closure", d=S10.get(q, 0.) * abs(v0), det=S10TXT.get(q, "WP4 RF scenario (a) bias"))], quad=S10.get(q, 0.) * abs(v0))
            continue
        if fam == "S12":
            e = S12.get(q, 0.) * abs(v0)
            fams[fam] = dict(d=e, rows=[dict(name="pion vertex", d=e, det="2.5 % (absolute) / 2.0 % (peak/bump ratio)" if e else "n/a")], quad=e,
                             pending="eps2 QED at M(ee) >= 20 MeV and (pT4, y4) reweighting: WP8" if q in ("BR_G", "a_eta", "sigma_eta", "a_P", "DR", "BR_A") else "")
            continue
        for g in items:
            rows.append(eval_item(g, q, fam))
        for t in info:
            it = item(t, q)
            rows.append(dict(name=t + " (INFO)", d=0., D=it["D"], na=False, rej=False, info=True,
                             det="%s %.5g->%.5g (D %+.4g, d %.4g) not counted" % (t, it["v0"], it["v"], it["D"], it["d"])))
        dmax = max([r["d"] for r in rows if not r.get("info")] + [0.])
        dmax_norej = max([r.get("d_norej", r["d"]) for r in rows if not r.get("info")] + [0.])
        dmax_stat = max([r.get("d_stat", r["d"]) for r in rows if not r.get("info")] + [0.])
        quad = dmax
        if fam in SUBSRC:
            sub = []
            for grp in SUBSRC[fam]:
                ds = [r["d"] for r in rows if not r.get("info") and r["name"].split("/")[0] in grp or r["name"].startswith("s free") and "SFREE" in grp]
                sub.append(max(ds + [0.]))
            quad = math.sqrt(sum(x * x for x in sub))
        fams[fam] = dict(d=dmax, d_norej=dmax_norej, d_stat=dmax_stat, rows=rows, quad=quad)
    syst = math.sqrt(sum(f["d"] ** 2 for f in fams.values()))
    syst_norej = math.sqrt(sum(f.get("d_norej", f["d"]) ** 2 for f in fams.values()))
    syst_quad = math.sqrt(sum(f["quad"] ** 2 for f in fams.values()))
    syst_stat = math.sqrt(sum(f.get("d_stat", f["d"]) ** 2 for f in fams.values()))
    mc = MCQ["nom"].get(q, MC.get(q, float("nan")))
    if abs(mc - MC.get(q, mc)) > 1e-3 * abs(mc):
        print("WARNING: MC of %s from the chain %.6g differs from WP5/WP6 %.6g" % (q, mc, MC[q]))
    nup, ndn = Nn.get(q, (0., 0.))
    RES[q] = dict(v=v0, stat=s0, fams=fams, syst=syst, syst_norej=syst_norej, syst_quad=syst_quad, syst_stat=syst_stat, mc=mc, Nup=nup * abs(v0), Ndn=ndn * abs(v0),
                  unit=unit, status=status)

# ---------------------------------------------------------------------------------------------- X17 (x17d, information)
X17 = dict(v=3.06e-5, mc=3.4e-6, rows={
    "S1": [("factorised CB", 0.983), ("per-row geometric mean", 0.992), ("k_E1 free", 1.018)],
    "S3": [("three classes", 0.744), ("flat pol2", 1.166), ("tilt only", 1.154), ("no M(ee) dependence", 1.156)],
    "S4": [("exp(pol1)", 1.006), ("exp(pol2)", 1.004), ("exp(pol3)", 0.863), ("kappa_B", 1.005)],
    "S6": [("r 0.89/1.10 (pair)", (abs(0.976 - 1) + abs(1.028 - 1)) / 2 + 1), ("s up/down (pair)", (abs(1.065 - 1) + abs(0.948 - 1)) / 2 + 1)],
    "S7": [("V2", 1.184), ("GS", 1.006)],
    "S12": [("e+ angle 1 + cos^2 (quoted separately)", 3.67 / 3.06)]})
X17F = {f: max(abs(r - 1) for _, r in rows) for f, rows in X17["rows"].items()}

# ---------------------------------------------------------------------------------------------- consistency tests
tests = []
def z2p(z):
    return math.erfc(abs(z) / math.sqrt(2))
# folds: yields scaled by the inverse row fraction (templates were scaled with it); scales and ratios as they are
f0 = Q.get("T_fold0", {}); f1 = Q.get("T_fold1", {})
fr = {}
for t in ("T_fold0", "T_fold1"):
    txt = cfg[t]["text"]
    fr[t] = float(txt.split("x ")[1].split(" ")[0])
FOLDSYS = {}
for q, unit, effc, fit, status in ALLQ:
    if q not in f0 or q not in f1:
        continue
    s0_, s1_ = (1. / fr["T_fold0"], 1. / fr["T_fold1"]) if unit == "rows" else (1., 1.)
    a, sa = f0[q][0] * s0_, f0[q][1] * s0_
    b, sb = f1[q][0] * s1_, f1[q][1] * s1_
    z = (a - b) / math.sqrt(sa * sa + sb * sb) if sa > 0 else float("nan")
    p = z2p(z)
    FOLDSYS[q] = abs(a - b) / 2 if p < 0.05 else 0.
    tests.append(("fold", q, a, sa, b, sb, z, p))
for q in RES:
    RES[q]["fold"] = FOLDSYS.get(q, 0.)
# isBest == 1 / isBest >= 0 (data-only ratios, D2; independent-sample bound for stat)
BEST = []
for q in ("Y_bump", "Y_bump_E2", "Y_peak", "a_P", "a_B", "DR", "a_eta", "Y5", "BR_G", "BR_A", "DR_F"):
    b, sb = val("T_best", q); a, sa = val("nom", q)
    r = b / a
    er = abs(r) * math.sqrt((sb / b) ** 2 + (sa / a) ** 2)
    BEST.append((q, r, er))
# cor_ok fallback (nested sample): z = Delta / sqrt(|sigma_v^2 - sigma_0^2|)
CORFB = []
for q in ("Y_bump", "Y_peak", "a_P", "a_B", "DR", "a_eta", "BR_G", "BR_A", "Y5"):
    b, sb = val("T_corfb", q); a, sa = val("nom", q)
    dd = math.sqrt(abs(sb * sb - sa * sa))
    CORFB.append((q, a, b, b - a, dd, (b - a) / dd if dd > 0 else float("nan")))

# periods (WP2 stability test): K_p = C_sec(p) x (1 - f_eta); old variant K_p = C_sec(p) (1 - f_win,p) with IIA/IIB correlated
PER = []
for line in open(os.path.join(OUT, "periods_raw.txt")):
    if line.startswith("#"):
        continue
    w = line.split()
    PER.append(dict(p=w[0], L=float(w[1]), f=float(w[2]), a_eta=(float(w[3]), float(w[4])), DR_F=(float(w[5]), float(w[6])),
                    a_P=(float(w[7]), float(w[8])), a_B=(float(w[9]), float(w[10])), Y_peak=(float(w[11]), float(w[12])),
                    Y_bump=(float(w[13]), float(w[14])), DR=(float(w[15]), float(w[16]))))
Csec = {"I": (1.0241, 0.0077), "IIA": (1.0, 0.0), "IIB": (1.0, 0.0), "III": (1.0, 0.0)}
fwin = {"I": (0.1929, 0.0369), "IIA": (0.1724, 0.0316), "IIB": (0.1724, 0.0316), "III": (0.1247, 0.0081)}
def chi2test(v, C):
    v = np.array(v); C = np.array(C)
    Ci = np.linalg.inv(C)
    one = np.ones(len(v))
    w = Ci @ one / (one @ Ci @ one)
    m = w @ v
    em = math.sqrt(1. / (one @ Ci @ one))
    r = v - m
    c2 = float(r @ Ci @ r)
    from math import exp
    # p-value of chi2 with n-1 dof
    k = len(v) - 1
    try:
        from scipy.stats import chi2 as _c
        p = float(_c.sf(c2, k))
    except Exception:
        p = float("nan")
    return m, em, c2, k, p
PERTEST = []
for q, kind in (("a_eta", "scale"), ("a_P", "scale"), ("a_B", "scale"), ("DR_F", "ratio"), ("DR", "ratio")):
    for variant in ("raw", "K (f_eta)", "K old (f_win, IIA-IIB correlated)"):
        if kind == "ratio" and variant != "raw":
            continue
        v, C = [], np.zeros((4, 4))
        Ks, Kc = [], np.zeros((4, 4))
        for i, P_ in enumerate(PER):
            x, sx = P_[q]
            if variant == "raw":
                K, sK = 1., 0.
            elif variant == "K (f_eta)":
                K, sK = Csec[P_["p"]][0] * 0.97, Csec[P_["p"]][1] * 0.97
            else:
                cs, scs = Csec[P_["p"]]; fw, sfw = fwin[P_["p"]]
                K = cs * (1 - fw); sK = math.hypot(scs * (1 - fw), cs * sfw)
            Ks.append((K, sK))
            v.append(x * K)
            C[i][i] = (sx * K) ** 2
        if variant.startswith("K old"):
            for i in range(4):
                for j in range(4):
                    pi, pj = PER[i]["p"], PER[j]["p"]
                    same = (i == j) or ({pi, pj} == {"IIA", "IIB"})
                    if same:
                        C[i][j] += PER[i][q][0] * PER[j][q][0] * Ks[i][1] * Ks[j][1]
        elif variant == "K (f_eta)":
            for i in range(4):
                C[i][i] += (PER[i][q][0] * Ks[i][1]) ** 2
        m, em, c2, k, p = chi2test(v, C)
        PERTEST.append((q, variant, v, [math.sqrt(C[i][i]) for i in range(4)], m, em, c2, k, p))

# ---------------------------------------------------------------------------------------------- outputs
def fmt(x, q):
    if x != x:
        return "n/a"
    if q in ("BR_G", "BR_A"):
        return "%.3e" % x
    if q in ("Y_bump", "Y_peak", "Y5", "Y4", "Y_bump_E2", "Y_bump_F", "Y_peak_F"):
        return "%.0f" % x
    return "%.4g" % x

FAMS = [f[0] for f in FAM]
lines = []
lines.append("# syst.txt -- WP7 systematics of SPEC report_all (plots/all_report/wp7_syst.C + wp7_combine.py), 2026-11-14")
lines.append("# S0, isBest >= 0, L = 6079.2 nb^-1; yields = raw rows x 64/63 (not efficiency corrected).")
lines.append("# rules: delta^2 = max(0, Delta^2 - |sig_v^2 - sig_0^2|); family = max over items (pair: (|D+|+|D-|)/2); families in quadrature;")
lines.append("#        MC and N apart; S11 only efficiency-corrected quantities and only if |Delta| > 2 sqrt(|sig_v^2 - sig_0^2|);")
lines.append("#        WP7 rule R1: same-bin variants with (chi2_v - chi2_0)/(chi2_0/ndf_0) - (ndf_0 - ndf_v) > 25 are listed, not counted.")
lines.append("# Barlow: delta^2 = max(0, Delta^2 - |stat_v^2 - stat_0^2| - MC(Delta)^2); MC(Delta) = robust spread (half of the central 68 %) of the")
lines.append("#        per-replica differences variant - nominal (100-replica chains; only variants on the S0__nom sample in the standard layout; else 0).")
lines.append("#        syst_statBarlow = the same with sigma = stat only (spec literal if sigma is read as statistical only).")
lines.append("# columns of 'syst' lines: quantity value stat " + " ".join(FAMS) + " syst(total) syst_noR1 syst_quadSubsrc syst_statBarlow MC N_up N_dn fold")
for q, unit, effc, fit, status in ALLQ:
    r = RES[q]
    lines.append("syst %s %.8g %.6g %s %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g" % (
        q, r["v"], r["stat"], " ".join("%.6g" % r["fams"][f]["d"] for f in FAMS), r["syst"], r["syst_norej"], r["syst_quad"], r["syst_stat"], r["mc"],
        r["Nup"], r["Ndn"], r["fold"]))
lines.append("# relative (%) of the same")
for q, unit, effc, fit, status in ALLQ:
    r = RES[q]
    a = abs(r["v"])
    lines.append("rel %s stat %.3f %s syst %.3f MC %.3f N +%.3f/-%.3f" % (q, 100 * r["stat"] / a, " ".join("%s %.3f" % (f, 100 * r["fams"][f]["d"] / a) for f in FAMS),
                                                                     100 * r["syst"] / a, 100 * r["mc"] / a, 100 * r["Nup"] / a, 100 * r["Ndn"] / a))
lines.append("# UL_X17 (x17d, information: these choices are inside the frozen limit through the N5 bias 15 % and sX 0.247)")
lines.append("x17 UL %.3e MC %.2e %s" % (X17["v"], X17["mc"], " ".join("%s %.4f" % (f, X17F[f]) for f in X17F)))
lines.append("# items: quantity family item delta(counted) details")
for q, unit, effc, fit, status in ALLQ:
    for f in FAMS:
        for row in RES[q]["fams"][f]["rows"]:
            lines.append("item %s %s %s | %.6g | %s" % (q, f, row["name"].replace(" ", "_"), row["d"], row.get("det", "")))
lines.append("# periods: quantity variant values... errors... mean err chi2 ndf p")
for q, variant, v, e, m, em, c2, k, p in PERTEST:
    lines.append("period %s %s | %s | %s | mean %.5g +- %.3g chi2 %.2f ndf %d p %.3g" % (q, variant.replace(" ", "_"), " ".join("%.5g" % x for x in v),
                                                                                    " ".join("%.3g" % x for x in e), m, em, c2, k, p))
lines.append("# folds: quantity fold0 err fold1 err z p (yields scaled by the inverse row fraction)")
for t in tests:
    lines.append("fold %s %.5g %.3g %.5g %.3g z %+.2f p %.3g" % t[1:])
lines.append("# isBest==1 / isBest>=0 (data only, D2): quantity ratio err(independent-sample bound)")
for q, r, er in BEST:
    lines.append("best %s %.4f %.4f" % (q, r, er))
lines.append("# cor_ok fallback: quantity nominal fallback Delta sqrt|sv^2-s0^2| z")
for x in CORFB:
    lines.append("corfb %s %.5g %.5g %+.4g %.4g %+.2f" % x)
open(os.path.join(OUT, "syst.txt"), "w").write("\n".join(lines) + "\n")

# ---------------------------------------------------------------------------------------------- Table S1 (txt, tex)
TS1 = ["Y_bump", "Y_peak", "Y5", "Y4", "a_eta", "BR_G", "BR_A", "kappa", "DR", "sigma_eta"]
HEAD = {"Y_bump": "Y_bump", "Y_peak": "Y_peak", "Y5": "Y5", "Y4": "Y4", "a_eta": "a_eta", "BR_G": "BR_G", "BR_A": "BR_A",
        "kappa": "kappa", "DR": "DR", "sigma_eta": "sigma_eta"}
def pct(x, q):
    a = abs(RES[q]["v"])
    return 100. * x / a
T = []
T.append("Table S1 (WP7): systematic uncertainties in per cent of the value. S0, isBest >= 0, L = 6079.2 nb^-1; yields raw rows x 64/63.")
T.append("Y_bump = B4 (all M(ee), E2 extrapolated; MAIN, EXTRAPOLATED); Y_peak = peak side bands (MAIN); Y5, Y4, a_eta, BR_A, kappa = F-main,")
T.append("BR_G = F-gold (CROSS-CHECKS: WP4 FAIL); DR = side bands / B4 (MAIN); sigma_eta = a_P x 1.0665 mb x C_sec(1 - f_eta) (MAIN: side bands).")
T.append("")
hdr = "%-26s" % "" + "".join("%11s" % HEAD[q] for q in TS1) + "%11s" % "UL_X17"
T.append(hdr)
T.append("%-26s" % "value" + "".join("%11s" % fmt(RES[q]["v"], q) for q in TS1) + "%11s" % ("%.2e" % X17["v"]))
T.append("%-26s" % "stat" + "".join("%11.1f" % pct(RES[q]["stat"], q) for q in TS1) + "%11s" % "(in UL)")
for fam, title, items, info in FAM:
    xs = ""
    for q in TS1:
        d = RES[q]["fams"][fam]["d"]
        na = all(r.get("na", False) for r in RES[q]["fams"][fam]["rows"] if not r.get("info")) if RES[q]["fams"][fam]["rows"] else False
        xs += "%11s" % ("n/a" if (fam == "S11" and not QINFO[q][2]) else "%.1f" % pct(d, q))
    xv = ("%.1f" % (100 * X17F[fam]) + "*") if fam in X17F else "-"
    T.append("%-26s" % ("%s %s" % (fam, title)) + xs + "%11s" % xv)
T.append("%-26s" % "total syst (quadrature)" + "".join("%11.1f" % pct(RES[q]["syst"], q) for q in TS1) + "%11s" % "in limit")
T.append("%-26s" % "  without rule R1" + "".join("%11.1f" % pct(RES[q]["syst_norej"], q) for q in TS1) + "%11s" % "")
T.append("%-26s" % "  Barlow with stat only" + "".join("%11.1f" % pct(RES[q]["syst_stat"], q) for q in TS1) + "%11s" % "")
T.append("%-26s" % "MC (template statistics)" + "".join("%11.1f" % pct(RES[q]["mc"], q) for q in TS1) + "%11.1f" % (100 * X17["mc"] / X17["v"]))
T.append("%-26s" % "N (absolute only)" + "".join("%11s" % ("-" if RES[q]["Nup"] == 0 else "+%.1f/-%.1f" % (pct(RES[q]["Nup"], q), pct(RES[q]["Ndn"], q))) for q in TS1) + "%11s" % "-")
T.append("%-26s" % "folds (p < 0.05 only)" + "".join("%11s" % ("%.1f" % pct(RES[q]["fold"], q) if RES[q]["fold"] > 0 else "-") for q in TS1) + "%11s" % "-")
T.append("")
T.append("Absolute: value +- stat +- syst +- MC [+- N]:")
for q in TS1:
    r = RES[q]
    s = "  %-10s %s +- %s (stat) +- %s (syst) +- %s (MC)" % (q, fmt(r["v"], q), fmt(r["stat"], q), fmt(r["syst"], q), fmt(r["mc"], q))
    if r["Nup"] > 0:
        s += " +%s/-%s (N)" % (fmt(r["Nup"], q), fmt(r["Ndn"], q))
    if q == "sigma_eta":
        s += " mb; live time +%.3f/-0 mb" % (0.073 * r["v"])
    T.append(s + "   [" + r["status"] + "]")
T.append("  UL_X17     3.06e-5 (90 % CLs, 16.7-17.0 MeV/c2, normalisation (1)); MC 3.4e-6 (11 %); e+ angle 1 + cos^2: 3.67e-5 (+20 %, S12);")
T.append("             * = x17d information only: these variants are inside the frozen limit (N5 bias 15 %, sX 0.247).")
T.append("S12 pending (WP8): eps2 QED at M(ee) >= 20 MeV (BR_G), (pT4, y4) reweighting (absolute scales).")
open(os.path.join(TAB, "wp7_tableS1.txt"), "w").write("\n".join(T) + "\n")

# LaTeX snippet (booktabs)
TEXH = {"Y_bump": r"$Y_\mathrm{bump}$", "Y_peak": r"$Y_\mathrm{peak}$", "Y5": r"$Y_5$", "Y4": r"$Y_4$", "a_eta": r"$a_\eta$",
        "BR_G": r"$\mathrm{BR}_G$", "BR_A": r"$\mathrm{BR}_A$", "kappa": r"$\kappa$", "DR": r"DR", "sigma_eta": r"$\sigma_\eta$"}
tex = []
tex.append(r"% wp7_tableS1.tex -- generated by plots/all_report/wp7_combine.py (WP7); relative uncertainties in per cent")
tex.append(r"\begin{tabular}{l" + "r" * (len(TS1) + 1) + "}")
tex.append(r"\toprule")
tex.append(" & " + " & ".join(TEXH[q] for q in TS1) + r" & UL$_{X17}$ \\")
tex.append(r"\midrule")
def texv(q):
    v = RES[q]["v"]
    if q in ("BR_G", "BR_A"):
        return r"$%.2f\cdot10^{-4}$" % (v * 1e4)
    if q == "sigma_eta":
        return r"%.3f\,mb" % v
    return fmt(v, q)
tex.append("value & " + " & ".join(texv(q) for q in TS1) + r" & $3.06\cdot10^{-5}$ \\")
tex.append("stat & " + " & ".join("%.1f" % pct(RES[q]["stat"], q) for q in TS1) + r" & -- \\")
tex.append(r"\midrule")
for fam, title, items, info in FAM:
    cells = []
    for q in TS1:
        if fam == "S11" and not QINFO[q][2]:
            cells.append("--")
        else:
            cells.append("%.1f" % pct(RES[q]["fams"][fam]["d"], q))
    xv = ("(%.1f)" % (100 * X17F[fam])) if fam in X17F else "--"
    tex.append("%s %s & " % (fam, title.replace("G_eta", r"$G_\eta$")) + " & ".join(cells) + " & " + xv + r" \\")
tex.append(r"\midrule")
tex.append(r"total syst & " + " & ".join("%.1f" % pct(RES[q]["syst"], q) for q in TS1) + r" & in limit \\")
tex.append(r"MC & " + " & ".join("%.1f" % pct(RES[q]["mc"], q) for q in TS1) + r" & 11.0 \\")
tex.append(r"N & " + " & ".join(("--" if RES[q]["Nup"] == 0 else r"$^{+%.1f}_{-%.1f}$" % (pct(RES[q]["Nup"], q), pct(RES[q]["Ndn"], q))) for q in TS1) + r" & -- \\")
tex.append(r"\bottomrule")
tex.append(r"\end{tabular}")
open(os.path.join(TAB, "wp7_tableS1.tex"), "w").write("\n".join(tex) + "\n")

# all variants (txt)
V = ["WP7 variants: every item, every quantity (Delta = variant - reference; d = counted delta after the Barlow subtraction).", ""]
for q, unit, effc, fit, status in ALLQ:
    r = RES[q]
    V.append("=== %s = %s +- %s (stat)  [%s]" % (q, fmt(r["v"], q), fmt(r["stat"], q), status))
    for f in FAMS:
        V.append("  %s: counted %.4g (%.2f %%)%s" % (f, r["fams"][f]["d"], pct(r["fams"][f]["d"], q),
                 ("; " + r["fams"][f].get("pending", "")) if r["fams"][f].get("pending") else ""))
        for row in r["fams"][f]["rows"]:
            V.append("      %-28s d = %-10.4g %s%s" % (row["name"], row["d"], "REJECTED (R1) " if row.get("rej") else "", row.get("det", "")))
    V.append("")
open(os.path.join(TAB, "wp7_variants.txt"), "w").write("\n".join(V) + "\n")

# periods table
PT = ["Per-period stability test (WP2 inputs; WP7). Background fixed: F-main type = nominal background shape x one free scale, kappa and g",
      "fixed (a_eta, rho_d free); side-band type = nominal stage-A shape and F-main eta terms, only the period data change.",
      "Errors: data (event groups + CB replicas of the period). K = C_sec(period) x (1 - f_eta); old K = C_sec (1 - f_win) per period.",
      "Templates x L_p / L with the analysed L per period (WP2 Table N1: 2809.453, 372.610, 1078.022, 1819.110 nb^-1).", ""]
PT.append("%-8s %-36s %s   %s" % ("qty", "variant", "  ".join("%14s" % P_["p"] for P_ in PER), "mean, chi2/ndf, p"))
for q, variant, v, e, m, em, c2, k, p in PERTEST:
    PT.append("%-8s %-36s %s   %.4g +- %.3g, %.2f/%d, p = %.3g %s" % (q, variant, "  ".join("%7.4f+-%.4f" % (x, y) for x, y in zip(v, e)), m, em, c2, k, p,
                                                                    "PASS" if p > 0.01 else "FAIL"))
open(os.path.join(TAB, "wp7_periods.txt"), "w").write("\n".join(PT) + "\n")
ptex = [r"% wp7_periods.tex -- generated by plots/all_report/wp7_combine.py", r"\begin{tabular}{llrrrrl}", r"\toprule",
        r"quantity & factor & I & IIA & IIB & III & $\chi^2$/ndf ($p$) \\", r"\midrule"]
for q, variant, v, e, m, em, c2, k, p in PERTEST:
    ptex.append(r"%s & %s & %s & %.2f/%d (%.2g) \\" % (q.replace("_", r"\_"), variant.replace("_", r"\_"), " & ".join(r"$%.3f\pm%.3f$" % (x, y) for x, y in zip(v, e)), c2, k, p))
ptex += [r"\bottomrule", r"\end{tabular}"]
open(os.path.join(TAB, "wp7_periods.tex"), "w").write("\n".join(ptex) + "\n")

# consistency tests table
CT = ["Consistency tests (not combined unless a fold has p < 0.05). S0, isBest >= 0.", "", "Folds (event-group parity; yields x inverse row fraction):"]
for t in tests:
    CT.append("  %-12s fold0 %10.4g +- %-8.3g fold1 %10.4g +- %-8.3g z %+.2f p %.3g%s" % (t[1], t[2], t[3], t[4], t[5], t[6], t[7], "  -> counts" if t[7] < 0.05 else ""))
CT.append("")
CT.append("isBest == 1 / isBest >= 0 (data only, same templates; error = independent-sample bound):")
for q, r, er in BEST:
    CT.append("  %-12s %.4f +- %.4f" % (q, r, er))
CT.append("")
CT.append("cor_ok fallback (nested sample; z = Delta / sqrt|sig_v^2 - sig_0^2|):")
for x in CORFB:
    CT.append("  %-12s nominal %.5g  fallback %.5g  Delta %+.4g  sqrt|dsig2| %.4g  z %+.2f" % x)
open(os.path.join(TAB, "wp7_tests.txt"), "w").write("\n".join(CT) + "\n")

# ---------------------------------------------------------------------------------------------- Fig. S1
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 13, "axes.labelsize": 15, "axes.titlesize": 15, "xtick.labelsize": 13, "ytick.labelsize": 13,
                     "pdf.fonttype": 42})
fig = plt.figure(figsize=(16, 11))
gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 1.0], wspace=0.32)
ax = fig.add_subplot(gs[0, 0])
rowsH = FAMS + ["total syst", "stat", "MC"]
M = np.zeros((len(rowsH), len(TS1)))
for j, q in enumerate(TS1):
    r = RES[q]
    for i, f in enumerate(FAMS):
        M[i, j] = pct(r["fams"][f]["d"], q) if not (f == "S11" and not QINFO[q][2]) else np.nan
    M[len(FAMS), j] = pct(r["syst"], q)
    M[len(FAMS) + 1, j] = pct(r["stat"], q)
    M[len(FAMS) + 2, j] = pct(r["mc"], q)
Mc = np.where(np.isnan(M), 0, M)
im = ax.imshow(np.log10(np.clip(Mc, 0.1, 300)), cmap="YlOrRd", aspect="auto", vmin=-1, vmax=np.log10(300))
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        t = "n/a" if np.isnan(M[i, j]) else ("%.0f" % M[i, j] if M[i, j] >= 10 else "%.1f" % M[i, j])
        ax.text(j, i, t, ha="center", va="center", fontsize=10.5, color="black")
ax.set_xticks(range(len(TS1)))
ax.set_xticklabels([TEXH[q] for q in TS1], rotation=45, ha="right")
ax.set_yticks(range(len(rowsH)))
ax.set_yticklabels(["%s %s" % (f[0], f[1].replace("G_eta", r"$G_\eta$")) for f in [(x[0], x[1]) for x in FAM]] + ["total syst", "stat", "MC"])
ax.axhline(len(FAMS) - 0.5, color="k", lw=1.5)
ax.set_title("(a) relative uncertainty [%] per family")
cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
cb.set_label(r"$\log_{10}$(uncertainty [%])")
# (b) ribbons: nested bands around 1 (value / value)
ax2 = fig.add_subplot(gs[0, 1])
labs = []
for j, q in enumerate(TS1[::-1]):
    r = RES[q]
    a = abs(r["v"])
    st = r["stat"] / a; sy = r["syst"] / a; mc = r["mc"] / a
    nu = r["Nup"] / a; nd = r["Ndn"] / a
    lev = [("stat", st, st, "#08306b"), ("stat (+) syst", math.hypot(st, sy), math.hypot(st, sy), "#2171b5"),
           ("(+) MC", math.sqrt(st ** 2 + sy ** 2 + mc ** 2), math.sqrt(st ** 2 + sy ** 2 + mc ** 2), "#6baed6"),
           ("(+) N", math.sqrt(st ** 2 + sy ** 2 + mc ** 2 + nu ** 2), math.sqrt(st ** 2 + sy ** 2 + mc ** 2 + nd ** 2), "#c6dbef")]
    for name, up, dn, col in lev[::-1]:
        ax2.fill_betweenx([j - 0.35, j + 0.35], 1 - dn, 1 + up, color=col, label=name if j == 0 else None, lw=0)
    ax2.plot([1, 1], [j - 0.4, j + 0.4], color="k", lw=1.5)
    labs.append(TEXH[q])
ax2.set_yticks(range(len(TS1)))
ax2.set_yticklabels(labs)
ax2.set_xscale("linear")
ax2.set_xlim(-1.2, 3.2)
ax2.set_xlabel("value / nominal value")
ax2.axvline(0, color="grey", lw=0.8, ls=":")
ax2.set_title("(b) nested uncertainty ribbons")
h, l = ax2.get_legend_handles_labels()
ax2.legend(h[::-1], l[::-1], loc="lower right", fontsize=12, frameon=True)
fig.text(0.01, 0.005, r"S0, isBest $\geq$ 0, L = 6079.2 nb$^{-1}$. MAIN: $Y_\mathrm{bump}$ (B4, extrapolated), $Y_\mathrm{peak}$ (side bands), DR, $\sigma_\eta$ (side bands); "
         r"CROSS-CHECK (WP4 FAIL): $Y_5$, $Y_4$, $a_\eta$, BR$_A$, $\kappa$ (F-main), BR$_G$ (F-gold). UL$_{X17}$: see Table S1.", fontsize=11.5)
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(FIG, "wp7_S1." + ext), dpi=110 if ext == "png" else None, bbox_inches="tight")
print(open(os.path.join(TAB, "wp7_tableS1.txt")).read())
print(open(os.path.join(TAB, "wp7_periods.txt")).read())
print(open(os.path.join(TAB, "wp7_tests.txt")).read())
