#!/usr/bin/env python3
# ============================================================================
# bedge_combine_v3.py  --  copy of bedge_combine.py for the rulings D38 and D40 (15 November 2026):
#   D38: reads bedge_v3/ana_results_v3.txt (bedge_ana_v3.C: C2 class factors of the per-row QED weights, E1 0.868178,
#        E2 0.853877); E8 = the model allowance +-0.004 on the effective 0.860; the factor's simulation-statistics error
#        (+-0.0098, c2recon) is added to the MC error of the C2-corrected quantities;
#   D40: new family E10 'peak tail', t = 1 +- 0.5 (physics band): the symmetric half difference f(1.5) - f(0.5) of
#        peaktail_bedge_v3.C (D38 factors), per quantity.
#   Tables: tables_all/bedge_{syst,results}_v3.{tex,txt}; json bedge_v3/bedge_numbers_v3.json. Toys and the angle toy
#   are read from the bedge directory unchanged.
# bedge_combine.py  --  package bedge: systematic families, conversions (yield, sigma), tests and tables
# ============================================================================
# Written 2026-11-14 (prefix bedge). Reads /hdd3/przygoda/HADES/FAT_out/all_report/bedge/ana_results.txt (bedge_ana.C) and
# toys_results.txt (bedge_toys.C); writes
#   docs/latex/tables_all/bedge_syst.{txt,tex}     every variant of f_eta(edge): value, shift, noise, counted delta, family
#   docs/latex/tables_all/bedge_results.{txt,tex}  the results (scale, yields, sigma) with stat / syst / MC / normalisation
#   docs/latex/tables_all/bedge_toys.{txt,tex}     pseudo-data closure
#   /hdd3/.../bedge/bedge_numbers.json            numbers for bedge_figs.py
# Part III rules (WP7): delta^2 = max(0, Delta^2 - |stat_v^2 - stat_0^2| - MC(Delta)^2), MC(Delta) = robust spread of the
# replica-by-replica difference; a +-1 sigma pair counts (|D+| + |D-|)/2; family = maximum over its items; families in
# quadrature; MC and normalisation separate; variants rejected by the data (rule R1) listed, not counted. Next to it the
# "bootstrap" version, where the data part of the noise is the bootstrap spread of the difference (same data replicas).
# Run: python3 -I plots/all_report/bedge_combine.py
# ============================================================================
import json
import math
from scipy.stats import chi2 as chi2dist

D = "/hdd3/przygoda/HADES/FAT_out/all_report/bedge/"
D3 = "/hdd3/przygoda/HADES/FAT_out/all_report/bedge_v3/"
TAB = "/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all/"
R = [l.split() for l in open(D3 + "ana_results_v3.txt") if not l.startswith("#")]
# D40: peak-tail terms (t = 1 +- 0.5) from peaktail_bedge_v3.C, D38 factors
PT3 = {}
for l in open(D3 + "peaktail_bedge_v3.txt"):
    a = l.split()
    if a and a[0] == "PT3" and a[1] == "D38":
        PT3[float(a[3])] = dict(f_all=float(a[7]), fc_rec=float(a[10]), a_k48=float(a[13]), a_k48_unc=float(a[16]))
TAIL_REL = {q: 0.5 * (PT3[1.5][q] - PT3[0.5][q]) / PT3[1.0][q] for q in ("f_all", "fc_rec", "a_k48", "a_k48_unc")}
C2STAT = 0.0098          # c2recon: simulation-statistics error of the per-slice factor (2000 group-Poisson replicas)
C2MODEL = 0.004          # D38 model allowance
C2EFF = 0.8599           # c2recon: effective factor on the per-slice C2 step

# ---------------------------------------------------------------- inputs (typed, with their source)
SIG_SMASH = 1.0694      # mb, sigma(pp -> eta X) of the simulated sample (v8)
CSEC_B = 1.0104         # C_sec of the bump window (WP2 / Amendment 1, task statement)
FETA = 0.97             # 1 - f_eta, +0.03/-0.04 (Amendment 1)
PDG26_C23 = 1.00436     # C2, C3 rescaling to PDG 2026 (BR(eta -> pi+pi-pi0) 22.92 -> 23.02 %), v8 V3
L = 6079.2              # nb^-1
PT3 = 64.0 / 63.0
N_PARTS = dict(L=6.9, period=5.0, BR=3.17, csec_up=3.1, csec_dn=4.1)   # per cent (WP7 N row, fix2 period term, v8 BR of the bump)
KAP_PEAK, KAP_PEAK_ERR = 0.48, 0.11   # conversion scale of the peak: 0.237 +- 0.049 / 0.497 +- 0.059 (fix_numbers.txt)
STAGEA = dict(c=[0.9816, 0.9844, 0.9505, 1.0473], b=-0.0639, kB=2.932)  # stage A at the WP5 eta values (wp6_results.md 2.2)
# comparison values
CMP = dict(peak=(0.429, 0.036, 0.152, 0.025, 0.041), direct=(0.522, 0.062), conv=(0.249, 0.051),
           tag_scale=(0.454, 0.014, 0.040), fmain_scale=(0.382, 0.043), fmain_sigma=0.399, b4_scale=0.312,
           aP=(0.4091, 0.0348), direct_scale=(0.497, 0.059), conv_scale=(0.237, 0.049))

# ---------------------------------------------------------------- parse
nomq, var, vq, vsl, per, fold, win, yld, bkg, bkgalt, bkgsl, wincfg = {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}
for r in R:
    t = r[0]
    if t == "NOMQ":
        nomq[r[1]] = dict(v=float(r[2]), statfit=float(r[4]), statboot=float(r[6]), boot68=float(r[8]), bootmean=float(r[10]),
                          mc=float(r[12]), mc68=float(r[14]), repmean=float(r[16]))
    elif t == "VAR":
        bar = r.index("|")
        var[r[1]] = dict(fam=r[3], pair=int(r[5]), counted=int(r[7]), r1=int(r[9]), chi2=float(r[11]), ndf=int(r[13]),
                         chi2_0=float(r[15]), ndf_0=int(r[17]), ok=int(r[19]), boot=int(r[21]), text=" ".join(r[bar + 1:]))
    elif t == "VQ":
        vq[(r[1], r[2])] = dict(v=float(r[3]), statfit=float(r[5]), statboot=float(r[7]), mc=float(r[9]), d=float(r[11]),
                                nb=float(r[13]), nr=float(r[15]), nb68=float(r[17]), nr68=float(r[19]))
    elif t == "VSL":
        vals = list(map(float, r[3:]))
        vsl[(r[1], int(r[2]))] = [vals[i:i + 3] for i in range(0, len(vals), 3)]
    elif t == "PER":
        d = {r[i]: (float(r[i + 1]), float(r[i + 2])) for i in range(6, len(r), 3)}
        d["Lp"] = float(r[3]); d["Csec"] = float(r[5]); per[r[1]] = d
    elif t == "FOLD":
        d = {r[i]: (float(r[i + 1]), float(r[i + 2])) for i in range(4, len(r), 3)}
        fold[int(r[1])] = d
    elif t == "WIN":
        win[(r[1], int(r[2]))] = dict(C2=float(r[4]), C3=float(r[6]), mc=float(r[14]))
    elif t == "WINCFG":
        wincfg[(r[1], r[2])] = list(map(float, r[3:]))
    elif t == "YIELD":
        d = dict(W2=float(r[4]), W3=float(r[6]), c2=float(r[8]))
        for j in range(3):
            o = 9 + 8 * j
            d["Y%d" % (j + 1)] = dict(v=float(r[o + 1]), stat=float(r[o + 3]), mc=float(r[o + 5]), mc68=float(r[o + 7]))
        yld[(r[1], int(r[2]))] = d
    elif t == "BKG":
        bkg[int(r[1])] = {r[i]: float(r[i + 1]) for i in range(2, len(r), 2)}
    elif t == "BKGALT":
        bkgalt[int(r[1])] = {r[i]: float(r[i + 1]) for i in range(2, len(r), 2)}
    elif t == "BKGSL":
        v = list(map(float, r[2:])); bkgsl[int(r[1])] = [v[i:i + 9] for i in range(0, len(v), 9)]
    elif t == "NOM":
        pass
nom_lin = {}
for r in R:
    if r[0] == "NOM":
        nom_lin[int(r[1])] = dict(Ed=float(r[3]), eEd=float(r[5]), Es=float(r[7]), T12=float(r[13]), C2=float(r[15]), C3=float(r[17]))

out = []
def P(s=""):
    print(s); out.append(s)

# ---------------------------------------------------------------- family assignment (smearMatch is a resolution item: E3)
FAM = {k: v["fam"] for k, v in var.items()}
FAM["smearMatch"] = "E3"
FAMTEXT = {"E1": "edge model (smearing, step shape, polynomial, ranges, edge position, fit weights)", "E2": "M(pipi) slices",
           "E3": "mass transfer s, M(pipi) scale, resolution r_res", "E4": "CB normalisation k, cells, rule",
           "E5": "physics weights (eta -> 3pi Dalitz plot), ownership", "E6": "selection (pre vs S0)",
           "E7": "background step projection (odd part of the smooth background)", "E8": "C2 angle: model allowance of the per-row QED weights (+-0.004, D38)",
           "E9": "conversion scale kappa imported from the peak", "E10": "peak tail, t = 1 +- 0.5 (D40)"}
# weights-off / phase space: excluded by the measured Dalitz plot (WP7 rejected phase space by rule R1)
NOTCOUNT = {"wtsOff": "R1: includes the phase-space Dalitz plot; per-slice data/simulation pattern chi2 65/3 (nominal 4.5/3), WP7 R1",
            "bumpPS": "R1: phase-space Dalitz plot; per-slice data/simulation pattern chi2 67/3 (nominal 4.5/3), WP7 R1"}

def r1_reject(v):
    if v["r1"] == 0:
        return False
    c0, n0, cv, nv = v["chi2_0"], v["ndf_0"], v["chi2"], v["ndf"]
    if v["r1"] == 1:
        return (cv - c0) / (c0 / n0) - (n0 - nv) > 25
    return (cv / nv) / (c0 / n0) > 1 + 5 * math.sqrt(2.0 / nv)

def combine(q):
    """systematic items and families for quantity q; returns (families literal, families bootstrap, rows)"""
    nq = nomq[q]
    rows = []
    for tag, v in var.items():
        x = vq.get((tag, q))
        if x is None:
            continue
        d = x["d"]
        lit_noise2 = abs(x["statfit"] ** 2 - nq["statfit"] ** 2) + x["nr68"] ** 2
        if v["boot"]:
            boot_noise2 = x["nb68"] ** 2 + x["nr68"] ** 2
        else:
            boot_noise2 = lit_noise2
        dl = math.sqrt(max(0.0, d * d - lit_noise2))
        db = math.sqrt(max(0.0, d * d - boot_noise2))
        rej = r1_reject(v)
        cnt = v["counted"] == 1 and not rej and tag not in NOTCOUNT
        why = "" if cnt else ("R1: data reject (chi2 %.1f/%d vs %.1f/%d)" % (v["chi2"], v["ndf"], v["chi2_0"], v["ndf_0"]) if rej else
                              NOTCOUNT.get(tag, "check only"))
        rows.append(dict(tag=tag, fam=FAM[tag], pair=v["pair"], v=x["v"], d=d, lit=math.sqrt(lit_noise2), boot=math.sqrt(boot_noise2),
                         dl=dl, db=db, cnt=cnt, why=why, text=v["text"]))
    fams_l, fams_b = {}, {}
    for f in sorted(set(r["fam"] for r in rows)):
        if f == "CHK":
            continue
        items = [r for r in rows if r["fam"] == f and r["cnt"]]
        vals_l, vals_b, seen = [], [], set()
        for r in items:
            if r["pair"] >= 0:
                if r["pair"] in seen:
                    continue
                pr = [s for s in items if s["pair"] == r["pair"]]
                seen.add(r["pair"])
                vals_l.append(sum(s["dl"] for s in pr) / len(pr)); vals_b.append(sum(s["db"] for s in pr) / len(pr))
            else:
                vals_l.append(r["dl"]); vals_b.append(r["db"])
        fams_l[f] = max(vals_l) if vals_l else 0.0
        fams_b[f] = max(vals_b) if vals_b else 0.0
    return fams_l, fams_b, rows

# ---------------------------------------------------------------- E7: background step projection
def kproj(k):
    """expected bias of f for class k from the step projection of a smooth background (rows per 2.5 MeV)."""
    b = bkg[k]
    B = nom_lin[k]["T12"] - (b["GDr"] + b["GCr"])                    # bump + other eta (linear decomposition)
    Ks = b["GDs"] + b["GCs"]                                            # SMASH-like, exp(pol3) 0.40-0.70
    Kp2 = bkgalt[k]["GDp2"] + bkgalt[k]["GCp2"]                         # exp(pol2)
    Kr43 = bkgalt[k]["GDr43"] + bkgalt[k]["GCr43"]                      # exp(pol3) on 0.43-0.63
    KA = 0.0                                                            # stage A: c_j (1 + b x)(G_D + kappa_B G_C)
    for j, sl in enumerate(bkgsl[k]):
        GDs, GDsx, GCs, GCsx = sl[0], sl[1], sl[2], sl[3]
        KA += STAGEA["c"][j] * (GDs + STAGEA["b"] * GDsx + STAGEA["kB"] * (GCs + STAGEA["b"] * GCsx))
    f = nomq["f_all" if k == 0 else ("f_E1" if k == 1 else "f_E2")]["v"]
    items = {"SMASH-like, exp(pol3)": (1 - f) * Ks / B, "SMASH-like, exp(pol2)": (1 - f) * Kp2 / B,
             "SMASH-like, 0.43-0.63": (1 - f) * Kr43 / B, "stage A (c_j, b, kappa_B 2.93)": (KA - f * Ks) / B}
    return dict(B=B, Ks=Ks, Kp2=Kp2, Kr43=Kr43, KA=KA, f=f, items=items)

P("# bedge_combine.py -- f_eta(edge) as a measurement (prefix bedge, 2026-11-14)")
P("")
P("## 1. background step projection (E7): rows per 2.5 MeV summed over the 4 slices; bias of f = (K_data - f K_sim)/B")
E7 = {}
for k, name in ((0, "all"), (1, "E1"), (2, "E2")):
    kp = kproj(k)
    E7[k] = kp
    P("class %-3s B %.1f  K smooth %.1f (pol2 %.1f, range %.1f), stage A %.1f;  bias of f: %s" % (
        name, kp["B"], kp["Ks"], kp["Kp2"], kp["Kr43"], kp["KA"], ", ".join("%s %+.4f" % (a, b) for a, b in kp["items"].items())))

# ---------------------------------------------------------------- C2 correction factor uncertainty (E8)
# spread of the factor: WP8 toy vs closure (half difference) and the toy's lepton-momentum threshold (0.04 / 0.10 GeV/c:
# edge R 0.904 / 0.859 against 0.883, bedge_angle.txt): +-0.023 on c; on fc the relative effect is s_C2 x dc / c_eff
ang = {}
for l in open(D + "bedge_angle.txt"):
    if l.startswith("R "):
        a = l.split(); ang[(a[1], a[2])] = dict(all=float(a[4]), edge=float(a[6]), err=float(a[8]))
dc_old = max(abs(ang[("p_e_>_0.04", "all")]["edge"] - ang[("nominal", "all")]["edge"]),
         abs(ang[("p_e_>_0.10", "all")]["edge"] - ang[("nominal", "all")]["edge"]), 0.5 * abs(0.880 - 0.891))
dc = C2MODEL      # D38
sC2 = nom_lin[0]["C2"] / nom_lin[0]["T12"]
cav = C2EFF      # D38 (was the WP8 mean 0.8855)
def e8_rel(q):
    if q in ("fc_rec", "a_k48", "fc_E1", "fc_E2", "a_D"):
        den = cav * sC2 + (1 - sC2) if q != "a_k48" else cav * sC2 + KAP_PEAK * nom_lin[0]["C3"] / nom_lin[0]["T12"] + (1 - sC2 - nom_lin[0]["C3"] / nom_lin[0]["T12"])
        return sC2 * dc / den
    return 0.0
P("")
P("## 2. C2 angle correction (D38): allowance dc = %.4f on c = %.4f (D35 spread was %.4f; toy edge R: nominal %.4f, p_e 0.04 %.4f, p_e 0.10 %.4f)" % (
    dc, cav, dc_old, ang[("nominal", "all")]["edge"], ang[("p_e_>_0.04", "all")]["edge"], ang[("p_e_>_0.10", "all")]["edge"]))
P("   edge-specific toy factors (all / E1 / E2): %.4f / %.4f / %.4f; all-M4 toy %.4f / %.4f / %.4f" % (
    ang[("nominal", "all")]["edge"], ang[("nominal", "E1")]["edge"], ang[("nominal", "E2")]["edge"],
    ang[("nominal", "all")]["all"], ang[("nominal", "E1")]["all"], ang[("nominal", "E2")]["all"]))
P("   C2 share of the simulated step (all M(ee), data covariance): %.4f; C3 share %.4f" % (sC2, nom_lin[0]["C3"] / nom_lin[0]["T12"]))

# ---------------------------------------------------------------- families per quantity
QUANT = ["f_all", "fc_rec", "a_k48", "a_k48_unc", "f_E1", "f_E2", "fc_E1", "fc_E2"]
RES = {}
P("")
P("## 3. systematic families (absolute; literal WP7 rule | bootstrap-noise version)")
for q in QUANT:
    fl, fb, rows = combine(q)
    v = nomq[q]["v"]
    k = 1 if q in ("f_E1", "fc_E1") else (2 if q in ("f_E2", "fc_E2") else 0)
    e7 = max(abs(x) for x in E7[k]["items"].values()) * (v / E7[k]["f"])
    fl["E7"] = fb["E7"] = e7
    e8 = e8_rel(q) * v
    if e8 > 0:
        fl["E8"] = fb["E8"] = e8
    if q in ("a_k48",):
        s3 = nom_lin[0]["C3"] / nom_lin[0]["T12"]
        den = cav * sC2 + KAP_PEAK * s3 + (1 - sC2 - s3)
        e9 = v * s3 * KAP_PEAK_ERR / den
        fl["E9"] = fb["E9"] = e9
    if q == "a_k48_unc":
        s3 = nom_lin[0]["C3"] / nom_lin[0]["T12"]
        den = sC2 + KAP_PEAK * s3 + (1 - sC2 - s3)
        fl["E9"] = fb["E9"] = v * s3 * KAP_PEAK_ERR / den
    if q in TAIL_REL:
        fl["E10"] = fb["E10"] = TAIL_REL[q] * v
    mc_add = (e8_rel(q) * C2STAT / C2MODEL * v) if e8_rel(q) > 0 else 0.0     # factor statistics into MC (D38)
    tl = math.sqrt(sum(x * x for x in fl.values())); tb = math.sqrt(sum(x * x for x in fb.values()))
    RES[q] = dict(v=v, stat=nomq[q]["statboot"], statfit=nomq[q]["statfit"], mc=math.hypot(nomq[q]["mc"], mc_add), mc68=math.hypot(nomq[q]["mc68"], mc_add), mc_rep=nomq[q]["mc"], mc_add=mc_add,
                  fams=fl, famsb=fb, syst=tl, systb=tb, rows=rows, repmean=nomq[q]["repmean"], bootmean=nomq[q]["bootmean"])
    tl_nt = math.sqrt(sum(x * x for k_, x in fl.items() if k_ != "E10"))
    RES[q]["syst_notail"] = tl_nt
    P("%-10s %.4f  stat %.4f (fit %.4f)  MC %.4f (replicas %.4f + factor %.4f; 68 %.4f)  syst %.4f (%.1f %%; without E10 %.4f = %.1f %%) | boot-noise %.4f" % (
        q, v, nomq[q]["statboot"], nomq[q]["statfit"], RES[q]["mc"], nomq[q]["mc"], mc_add, nomq[q]["mc68"], tl, 100 * tl / v, tl_nt, 100 * tl_nt / v, tb))
    P("           families: " + ", ".join("%s %.4f" % (f, fl[f]) for f in sorted(fl)))

# ---------------------------------------------------------------- variant table (f_all, fc_rec, a_k48)
P("")
P("## 4. variants of f_eta(edge) (all M(ee)); Delta, its noise (literal | bootstrap), counted delta (literal)")
for r in sorted(RES["f_all"]["rows"], key=lambda x: (x["fam"], x["tag"])):
    P("%-4s %-11s f %.4f  D %+.4f  noise %.4f | %.4f  delta %.4f%s   %s" % (
        r["fam"], r["tag"], r["v"], r["d"], r["lit"], r["boot"], r["dl"], "" if r["cnt"] else "  [not counted: %s]" % r["why"], r["text"]))

# ---------------------------------------------------------------- tests
P("")
P("## 5. tests")
calib = nomq["f_all"]["statboot"] / nomq["f_all"]["statfit"]
for q in ("f_all", "fc_rec", "f_E1", "f_E2"):
    vals = [(p, per[p][q][0], per[p][q][1] * calib, per[p]["Csec"]) for p in ("perI", "perIIA", "perIIB", "perIII")]
    for use_c in (0, 1):
        xs = [(v * (c if use_c else 1.0), e * (c if use_c else 1.0)) for (_, v, e, c) in vals]
        w = [1 / e ** 2 for _, e in xs]
        m = sum(wi * x for wi, (x, _) in zip(w, xs)) / sum(w)
        ch = sum(((x - m) / e) ** 2 for x, e in xs)
        p = chi2dist.sf(ch, 3)
        RES.setdefault("periods", {})[(q, use_c)] = dict(vals=xs, mean=m, chi2=ch, p=p)
        P("periods %-7s %s: %s; mean %.4f, chi2 %.2f/3, p = %.3f %s" % (q, "x C_sec" if use_c else "raw    ",
          ", ".join("%.3f(%.0f)" % (x, 1000 * e) for x, e in xs), m, ch, p, "PASS" if p > 0.05 else ("MARGINAL" if p > 0.01 else "FAIL")))
f0, e0 = fold[0]["f_all"][0], fold[0]["f_all"][1] * calib
f1, e1 = fold[1]["f_all"][0], fold[1]["f_all"][1] * calib
zf = (f1 - f0) / math.sqrt(e0 ** 2 + e1 ** 2)
P("folds: %.4f +- %.4f / %.4f +- %.4f, z = %.2f, p = %.3f" % (f0, e0, f1, e1, zf, 2 * (1 - 0.5 * (1 + math.erf(abs(zf) / math.sqrt(2))))))
RES["folds"] = dict(f0=f0, e0=e0, f1=f1, e1=e1, z=zf)
for tag in ("best", "bestData", "pre"):
    x = vq[(tag, "f_all")]
    z = x["d"] / math.sqrt(x["nb68"] ** 2 + x["nr68"] ** 2)
    P("%-9s f %.4f, ratio to S0 %.3f, Delta %+.4f, noise %.4f, z %.2f   (%s)" % (tag, x["v"], x["v"] / nomq["f_all"]["v"], x["d"],
      math.sqrt(x["nb68"] ** 2 + x["nr68"] ** 2), z, var[tag]["text"]))
    RES[tag] = dict(v=x["v"], ratio=x["v"] / nomq["f_all"]["v"], z=z)
P("bootstrap / fit statistical error: S0 %.3f, pre %.3f, S0best %.3f (Part II halving of 'pre': 1.20)" % (
    calib, vq[("pre", "f_all")]["statboot"] / vq[("pre", "f_all")]["statfit"], vq[("best", "f_all")]["statboot"] / vq[("best", "f_all")]["statfit"]))
# per-slice ratio shape test (data stat only): nominal and the Dalitz-plot variants
P("per-slice data/simulation ratios (h_data/h_sim) and chi2 against a constant (data statistical errors):")
RES["slice_shape"] = {}
for tag in ("nominal", "bump1aY", "bumpPS", "wtsOff"):
    if tag == "nominal":
        sl = [(float(r[6]), float(r[8]), float(r[14])) for r in R if r[0] == "SLICE" and r[1] == "0"]
    else:
        sl = vsl[(tag, 0)]
    rat = [(h / hs, eh / hs) for h, eh, hs in sl]
    w = [1 / e ** 2 for _, e in rat]
    m = sum(wi * x for wi, (x, _) in zip(w, rat)) / sum(w)
    ch = sum(((x - m) / e) ** 2 for x, e in rat)
    RES["slice_shape"][tag] = dict(r=rat, chi2=ch, p=chi2dist.sf(ch, len(rat) - 1))
    P("   %-8s %s  chi2 %.2f/%d (p %.2f)" % (tag, " ".join("%.3f(%.0f)" % (x, 1000 * e) for x, e in rat), ch, len(rat) - 1, chi2dist.sf(ch, len(rat) - 1)))

# ---------------------------------------------------------------- toys
TOY = {}
for l in open(D + "toys_results.txt"):
    if l.startswith("#"):
        continue
    a = l.split()
    TOY[(a[0], a[1], a[2])] = dict(truth=float(a[3]), mean=float(a[4]), bias=float(a[5]), rel=float(a[6]), rms=float(a[7]),
                                   pm=float(a[8]), psd=float(a[9]), cov=float(a[10]))
P("")
P("## 6. pseudo-data closure (1000 toys per scenario; AA = template = truth sample, RF = template from replica)")
for sc in ("a1", "c1", "a", "null", "sA", "q-", "q--", "q+", "om", "wide", "narrow", "smooth", "mass", "qed", "kap", "datalike"):
    s = []
    for q in ("f_all", "fc_rec", "a_k48"):
        for m in ("AA", "RF"):
            t = TOY[(sc, m, q)]
            s.append("%s %s %+.3f (sd %.2f, cov %.2f)" % (q, m, t["rel"] if t["truth"] else t["bias"], t["psd"], t["cov"]))
    P("%-8s " % sc + "; ".join(s))

# ---------------------------------------------------------------- conversions
P("")
P("## 7. conversions")
def nrel():
    up = math.sqrt(N_PARTS["L"] ** 2 + N_PARTS["period"] ** 2 + N_PARTS["BR"] ** 2 + N_PARTS["csec_up"] ** 2)
    dn = math.sqrt(N_PARTS["L"] ** 2 + N_PARTS["period"] ** 2 + N_PARTS["BR"] ** 2 + N_PARTS["csec_dn"] ** 2)
    return up / 100, dn / 100
Nup, Ndn = nrel()
SIGMA = {}
for q, lab in (("f_all", "f_eta(edge), SMASH as is (kappa = 1), uncorrected"), ("fc_rec", "C2 angle corrected (per-row QED weights, D38), kappa = 1"),
               ("a_k48", "C2 corrected, conversions x 0.48 (peak kappa)"), ("a_k48_unc", "uncorrected, conversions x 0.48")):
    x = RES[q]
    fac = SIG_SMASH * CSEC_B * FETA / PDG26_C23
    s = x["v"] * fac
    SIGMA[q] = dict(v=s, stat=x["stat"] * fac, syst=x["syst"] * fac, mc=x["mc"] * fac, mc68=x["mc68"] * fac, nup=s * Nup, ndn=s * Ndn,
                    scale26=x["v"] / PDG26_C23, label=lab)
    tot = math.sqrt(x["stat"] ** 2 + x["syst"] ** 2 + x["mc"] ** 2) * fac
    SIGMA[q]["tot_noN"] = tot
    SIGMA[q]["tot"] = math.sqrt(tot ** 2 + (s * 0.5 * (Nup + Ndn)) ** 2)
    P("%-10s scale %.4f (PDG 2026 %.4f) -> sigma = %.3f +- %.3f (stat) +- %.3f (syst) +- %.3f (MC) +%.3f/-%.3f (norm.) mb; total %.3f (%.0f %%)  [%s]" % (
        q, x["v"], x["v"] / PDG26_C23, s, SIGMA[q]["stat"], SIGMA[q]["syst"], SIGMA[q]["mc"], SIGMA[q]["nup"], SIGMA[q]["ndn"], SIGMA[q]["tot"],
        100 * SIGMA[q]["tot"] / s, lab))
# a_D, kappa from the E1/E2 split (no systematics: MC-limited)
P("E1/E2 split (C2 corrected): a_D = %.3f +- %.3f (stat) +- %.3f (MC68; rms %.2f), kappa = %.3f +- %.3f +- %.3f (MC68); uncorrected a_D %.3f, kappa %.3f" % (
    nomq["a_D"]["v"], nomq["a_D"]["statboot"], nomq["a_D"]["mc68"], nomq["a_D"]["mc"], nomq["kappa"]["v"], nomq["kappa"]["statboot"], nomq["kappa"]["mc68"],
    nomq["a_D_unc"]["v"], nomq["kappa_unc"]["v"]))
SIGMA["a_D"] = dict(v=nomq["a_D"]["v"] * SIG_SMASH * CSEC_B * FETA / PDG26_C23, stat=nomq["a_D"]["statboot"] * SIG_SMASH * CSEC_B * FETA / PDG26_C23,
                    mc68=nomq["a_D"]["mc68"] * SIG_SMASH * CSEC_B * FETA / PDG26_C23)
P("   -> sigma(a_D) = %.3f +- %.3f (stat) +- %.3f (MC68) mb" % (SIGMA["a_D"]["v"], SIGMA["a_D"]["stat"], SIGMA["a_D"]["mc68"]))
# yields
YL = {}
for (w, k), d in sorted(yld.items()):
    for j, lab in ((1, "f x (C2 + C3)"), (2, "fc x (c C2 + C3)"), (3, "a_k48 x (c C2 + 0.48 C3)")):
        y = d["Y%d" % j]
        rel_syst = RES[["f_all", "fc_rec", "a_k48"][j - 1]]["syst"] / RES[["f_all", "fc_rec", "a_k48"][j - 1]]["v"]
        Y = y["v"] * PT3
        YL[(w, k, j)] = dict(Y=Y, stat=y["stat"] * PT3, mc=y["mc"] * PT3, mc68=y["mc68"] * PT3, syst=Y * rel_syst, nb=Y / L,
                             nb_stat=y["stat"] * PT3 / L, nb_mc=y["mc"] * PT3 / L, nb_syst=Y * rel_syst / L)
        if j in (1, 2) or w == "B":
            P("yield %-3s class %d  %-26s Y = %.0f +- %.0f (stat) +- %.0f (syst) +- %.0f (MC; 68 %.0f) rows = %.3f +- %.3f +- %.3f +- %.3f nb" % (
                w, k, lab, Y, y["stat"] * PT3, Y * rel_syst, y["mc"] * PT3, y["mc68"] * PT3, Y / L, y["stat"] * PT3 / L, Y * rel_syst / L, y["mc"] * PT3 / L))
P("   (B4 for comparison: Y_bump(all M(ee)) 16389 +- 1129 +- 25766 +- 4110, NOT DETERMINED; F-main 19426 +- 2478 (cross-check))")

# double ratio peak / bump with the edge (PDG 2026 both): DR_edge = a_P / f(edge)
aP, aP_stat, aP_syst_rel, aP_mc = 0.4091, 0.0348, 0.160, 0.025 * 0.4091 / 0.418   # fix_numbers / WP7 (syst without the composition term)
DRE = {}
for q in ("f_all", "fc_rec"):
    x = RES[q]; f26 = x["v"] / PDG26_C23
    dr = aP / f26
    st = dr * math.sqrt((aP_stat / aP) ** 2 + (x["stat"] / x["v"]) ** 2)
    sy = dr * math.sqrt(aP_syst_rel ** 2 + (x["syst_notail"] / x["v"]) ** 2)   # D40: the tail term of DR is listed (correlated ratio)
    mc = dr * math.sqrt((aP_mc / aP) ** 2 + (x["mc"] / x["v"]) ** 2)
    DRE[q] = dict(v=dr, stat=st, syst=sy, mc=mc)
    P("DR_edge (%s) = a_P / f = %.3f +- %.3f (stat) +- %.3f (syst) +- %.3f (MC)   [B4: 1.373 NOT DETERMINED; F-main 1.148]" % (q, dr, st, sy, mc))

# ---------------------------------------------------------------- tables
def fmt(x, n=3):
    return ("%." + str(n) + "f") % x
with open(TAB + "bedge_syst_v3.txt", "w") as f:
    f.write("# bedge_syst_v3.txt -- bedge_combine_v3.py (D38 + D40); variants of f_eta(edge) = E_data/E_sim (all M(ee)); S0, isBest >= 0\n")
    f.write("# Delta = variant - nominal (0.3916); noise = sqrt(|stat_v^2 - stat_0^2| + MC(Delta)^2) (literal WP7) | bootstrap version; delta = counted\n")
    for r in sorted(RES["f_all"]["rows"], key=lambda x: (x["fam"], x["tag"])):
        f.write("%-4s %-11s %.4f %+.4f %.4f %.4f %.4f %s | %s\n" % (r["fam"], r["tag"], r["v"], r["d"], r["lit"], r["boot"], r["dl"],
                                                               "counted" if r["cnt"] else "NOT COUNTED (%s)" % r["why"], r["text"]))
    for q in ("f_all", "fc_rec", "a_k48"):
        x = RES[q]
        f.write("FAMILIES %s: " % q + ", ".join("%s %.4f" % (k, v) for k, v in sorted(x["fams"].items())) + "; total %.4f (%.1f %%); bootstrap-noise total %.4f\n" % (
            x["syst"], 100 * x["syst"] / x["v"], x["systb"]))
TEXN = {"E1": "E1 edge model", "E2": "E2 $M(\\pi\\pi)$ slices", "E3": "E3 mass transfer, resolution", "E4": "E4 CB", "E5": "E5 physics weights",
        "E6": "E6 selection", "E7": "E7 background step", "E8": "E8 C2 factor, model ($\\pm0.004$, D38)", "E9": "E9 $\\kappa$ of the peak",
        "E10": "E10 peak tail, $t=1\\pm0.5$ (D40)"}
with open(TAB + "bedge_syst_v3.tex", "w") as f:
    f.write("% bedge_syst_v3.tex -- bedge_combine_v3.py (rulings D38, D40). Systematic families of the per-slice edge estimator, per cent of the value.\n")
    f.write("\\begin{tabular}{lrrr}\n\\hline\n family & $f_\\eta$(edge) & C2 corrected & $\\kappa=0.48$ \\\\\n\\hline\n")
    for fam in ("E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10"):
        cells = []
        for q in ("f_all", "fc_rec", "a_k48"):
            x = RES[q]
            cells.append("%.1f" % (100 * x["fams"][fam] / x["v"]) if fam in x["fams"] else "--")
        f.write(" %s & %s \\\\\n" % (TEXN[fam], " & ".join(cells)))
    f.write("\\hline\n total syst & %s \\\\\n" % " & ".join("%.1f" % (100 * RES[q]["syst"] / RES[q]["v"]) for q in ("f_all", "fc_rec", "a_k48")))
    f.write(" stat (bootstrap) & %s \\\\\n" % " & ".join("%.1f" % (100 * RES[q]["stat"] / RES[q]["v"]) for q in ("f_all", "fc_rec", "a_k48")))
    f.write(" MC (replicas; C2 factor statistics included) & %s \\\\\n" % " & ".join("%.1f" % (100 * RES[q]["mc"] / RES[q]["v"]) for q in ("f_all", "fc_rec", "a_k48")))
    f.write("\\hline\n\\end{tabular}\n")
with open(TAB + "bedge_results_v3.txt", "w") as f:
    f.write("# bedge_results_v3.txt -- bedge_combine_v3.py (D38 + D40); S0, isBest >= 0, L = 6079.2 nb^-1, PDG 2020 in the simulation (sigma rescaled to PDG 2026)\n")
    for q in ("f_all", "fc_rec", "a_k48", "a_k48_unc", "f_E1", "f_E2", "fc_E1", "fc_E2"):
        x = RES[q]
        f.write("scale %-10s %.4f stat %.4f syst %.4f MC %.4f (68 %.4f)\n" % (q, x["v"], x["stat"], x["syst"], x["mc"], x["mc68"]))
    for q, s in SIGMA.items():
        if "syst" in s:
            f.write("sigma %-10s %.4f mb stat %.4f syst %.4f MC %.4f norm +%.4f/-%.4f total %.4f | %s\n" % (q, s["v"], s["stat"], s["syst"], s["mc"], s["nup"], s["ndn"], s["tot"], s["label"]))
        else:
            f.write("sigma %-10s %.4f mb stat %.4f MC68 %.4f (E1/E2 split, no systematics)\n" % (q, s["v"], s["stat"], s["mc68"]))
    for (w, k, j), y in sorted(YL.items()):
        f.write("yield %s cls %d Y%d %.1f stat %.1f syst %.1f MC %.1f (68 %.1f) counts; %.4f visible counts per nb^-1\n" % (w, k, j, y["Y"], y["stat"], y["syst"], y["mc"], y["mc68"], y["nb"]))
with open(TAB + "bedge_results_v3.tex", "w") as f:
    f.write("% bedge_results_v3.tex -- bedge_combine_v3.py (rulings D38: per-row C2 factors; D40: peak-tail family). Counts = stored pi+pi-e+e- combinations.\n")
    f.write("\\begin{tabular}{lllll}\n\\hline\n quantity & SMASH as is ($\\kappa=1$) & C2 corrected & C2 corr., $\\kappa=0.48$ & \\\\\n\\hline\n")
    def cell(x, n=3):
        return "$%s\\pm%s\\pm%s\\pm%s$" % (fmt(x["v"], n), fmt(x["stat"], n), fmt(x["syst"], n), fmt(x["mc"], n))
    f.write(" $\\eta$ scale (bump edge) & %s & %s & %s & stat, syst, MC \\\\\n" % (cell(RES["f_all"]), cell(RES["fc_rec"]), cell(RES["a_k48"])))
    def ycell(y):
        return "$%.0f\\pm%.0f\\pm%.0f\\pm%.0f$" % (y["Y"], y["stat"], y["syst"], y["mc"])
    f.write(" $Y_{\\mathrm{bump}}$(B, all $M_{ee}$) [counts] & %s & %s & %s & stat, syst, MC \\\\\n" % (ycell(YL[("B", 0, 1)]), ycell(YL[("B", 0, 2)]), ycell(YL[("B", 0, 3)])))
    def ncell(y):
        return "$%.2f\\pm%.2f\\pm%.2f\\pm%.2f$" % (y["nb"], y["nb_stat"], y["nb_syst"], y["nb_mc"])
    f.write(" $Y/L$ (B, all $M_{ee}$) [visible counts per nb$^{-1}$] & %s & %s & %s & + L 8.5\\,\\%% \\\\\n" % (ncell(YL[("B", 0, 1)]), ncell(YL[("B", 0, 2)]), ncell(YL[("B", 0, 3)])))
    def scell(s):
        return "$%.3f\\pm%.3f\\pm%.3f\\pm%.3f^{+%.3f}_{-%.3f}$" % (s["v"], s["stat"], s["syst"], s["mc"], s["nup"], s["ndn"])
    f.write(" $\\sigma(pp\\to\\eta X)$ [mb] & %s & %s & %s & + norm. \\\\\n" % (scell(SIGMA["f_all"]), scell(SIGMA["fc_rec"]), scell(SIGMA["a_k48"])))
    f.write("\\hline\n\\end{tabular}\n")
with open("/dev/null", "w") as f:
    f.write("# bedge_toys.txt -- bedge_toys.C / bedge_combine.py; 1000 toys per scenario; relative bias (absolute for a = 0), pull SD, coverage of +-1 sigma\n")
    for (sc, m, q), t in sorted(TOY.items()):
        f.write("%-8s %s %-8s truth %.3f mean %.4f relbias %+.4f rms %.4f pullmean %+.3f pullsd %.3f cover %.3f\n" % (sc, m, q, t["truth"], t["mean"], t["rel"], t["rms"], t["pm"], t["psd"], t["cov"]))
TOYTEX = [("a1", "SMASH background $\\times1$, $a=1$ (exact closure expected)"), ("c1", "SMASH background $\\times1$, $a=0.4$"),
          ("a", "WP4 (a): background $\\times1.35$, tilt $-0.15$, $a=0.4$"), ("null", "(a) with $a=0$ (absolute bias)"),
          ("q-", "(a) + quadratic term $q=-0.173$ (B4: $+150$\\,\\%)"), ("q--", "(a) + $q=-0.346$"), ("q+", "(a) + $q=+0.173$"),
          ("om", "WP4 (e): $\\omega/\\rho^0$-pion background $\\times0.3$, rest $\\times1.6$"), ("sA", "stage-A background, F-main $\\eta$ values"),
          ("smooth", "(a) with the smoothed background (template noise)"), ("wide", "(a) $\\times(1+0.10\\,G(0.490, 15\\,\\mathrm{MeV}))$"),
          ("narrow", "(a) $\\times(1+0.05\\,G(0.505, 8\\,\\mathrm{MeV}))$"), ("mass", "truth at $s+0.0009$, $r=0.89$"),
          ("qed", "truth C2 $\\times0.8855$ (QED angle)"), ("kap", "truth conversions $\\times0.48$"), ("datalike", "$a=0.55$, C2 $\\times0.8855$, conversions $\\times0.48$")]
with open("/dev/null", "w") as f:
    f.write("% bedge_toys.tex -- bedge_combine.py (prefix bedge): relative bias in per cent (absolute for a = 0); in brackets the pull SD\n")
    f.write("\\begin{tabular}{lrrl}\n\\hline\n truth & $f_\\eta$(edge), AA & $f_\\eta$(edge), RF & matching estimator (AA / RF) \\\\\n\\hline\n")
    MATCH = {"qed": ("fc_rec", "C2 corrected"), "kap": ("a_k48_unc", "$\\kappa=0.48$"), "datalike": ("a_k48", "C2 corr., $\\kappa=0.48$")}
    for sc, lab in TOYTEX:
        a, rf = TOY[(sc, "AA", "f_all")], TOY[(sc, "RF", "f_all")]
        if sc in MATCH:
            q, nm = MATCH[sc]
            ma, mr = TOY[(sc, "AA", q)], TOY[(sc, "RF", q)]
            mcell = "%s: $%+.1f$ (%.2f) / $%+.1f$ (%.2f)" % (nm, 100 * ma["rel"], ma["psd"], 100 * mr["rel"], mr["psd"])
        else:
            mcell = "--"
        if a["truth"] == 0:
            f.write(" %s & $%+.3f$ (%.2f) & $%+.3f$ (%.2f) & %s \\\\\n" % (lab, a["bias"], a["psd"], rf["bias"], rf["psd"], mcell))
        else:
            f.write(" %s & $%+.1f$ (%.2f) & $%+.1f$ (%.2f) & %s \\\\\n" % (lab, 100 * a["rel"], a["psd"], 100 * rf["rel"], rf["psd"], mcell))
    f.write("\\hline\n\\end{tabular}\n")

# ---------------------------------------------------------------- json for the figures
J = dict(RES={q: {k: v for k, v in x.items() if k != "rows"} for q, x in RES.items() if q in QUANT},
         rows={q: RES[q]["rows"] for q in ("f_all", "fc_rec", "a_k48")}, SIGMA=SIGMA, CMP=CMP,
         YL={"%s_%d_%d" % k: v for k, v in YL.items()}, E7={k: dict(items=v["items"], Ks=v["Ks"], KA=v["KA"], B=v["B"]) for k, v in E7.items()},
         periods={"%s_%d" % k: v for k, v in RES["periods"].items()}, folds=RES["folds"], slice_shape=RES["slice_shape"],
         TOY={"%s|%s|%s" % k: v for k, v in TOY.items()}, nomq=nomq, nom_lin=nom_lin, ang={"%s|%s" % k: v for k, v in ang.items()},
         best=RES["best"], bestData=RES["bestData"], pre=RES["pre"], DRE=DRE)
J["TAIL_REL"] = TAIL_REL
json.dump(J, open(D3 + "bedge_numbers_v3.json", "w"), indent=1, default=float)
open(D3 + "bedge_combine_v3_stdout.txt", "w").write("\n".join(out) + "\n")
