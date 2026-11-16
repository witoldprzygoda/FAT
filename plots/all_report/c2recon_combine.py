#!/usr/bin/env python3
# ============================================================================
# c2recon_combine.py  --  package c2recon: statistics of the bedge factor, effects on f, sigma and the ECAL tag, tables
# ============================================================================
# Written 2026-11-15 (prefix c2recon). New file. Reads the outputs of c2recon_rows.C, c2recon_weights.py, c2recon_toy.py,
# c2recon_edgeu.py and c2recon_bedge.C, plus the published numbers of bedge / edgeu / ecalacc (read only), and writes
#   /hdd3/przygoda/HADES/FAT_out/all_report/c2recon/c2recon_numbers.json and c2recon_combine.log
#   /hdd3/przygoda/HADES/FAT_out/all_report/c2recon/c2recon_table.txt   (the factor table, readable)
# Run: python3 -I plots/all_report/c2recon_combine.py
# ============================================================================
import json
import math
import sys

import numpy as np

sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import c2recon_weights as CW
import edgeu_core as C

O = "/hdd3/przygoda/HADES/FAT_out/all_report/c2recon/"
LOG = open(O + "c2recon_combine.log", "w")


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s); LOG.write(s + "\n"); LOG.flush()


def sd(a):
    return C.sd(a)


NUM = {}
# ---------------------------------------------------------------- bedge: KEY lines and linear functionals
key = {}
for l in open(O + "c2recon_bedge.txt"):
    if not l.startswith("KEY"):
        continue
    t = l.split()
    if t[1] == "nom":
        v = list(map(float, [x for x in t[2:] if not x[0].isalpha()]))
        key["nom"] = dict(Ed=v[0:3], Es=v[3:6], linC2=v[6:9], linC3=v[9:12], linT12=v[12:15], f=v[15], frec=v[16])
    else:
        name = t[2]
        d = {}
        i = 3
        while i < len(t):
            k = t[i]; i += 1
            vals = []
            while i < len(t) and not t[i][0].isalpha():
                vals.append(float(t[i])); i += 1
            d[k] = vals if len(vals) > 1 else vals[0]
        key[name] = d
nomb = key["nom"]
P("# c2recon_combine.py, 2026-11-15")
P("bedge nominal: E_data %.4f, E_sim %.4f, f %.5f, f_rec %.5f" % (nomb["Ed"][0], nomb["Es"][0], nomb["f"], nomb["frec"]))
lin = {}
for l in open(O + "c2recon_bedge_lin.txt"):
    if l.startswith("#"):
        continue
    t = l.split()
    k, i, pa, pb = int(t[0]), int(t[1]), int(t[2]), int(t[3])
    lin[(k, i)] = (pa, pb, np.array(list(map(float, t[4:]))))
R = CW.load_rows()
Wt = np.load(O + "c2recon_weights.npz")
q = Wt["q_A"]
s0 = R["s0"] > 0
e = R["e"].astype(int); p = R["p"].astype(int); ib = R["ib"].astype(int)
W = R["W108"]
n = len(W)
a = np.zeros((3, n))
for (k, i), (pa, pb, arr) in lin.items():
    m = s0 & (p >= pa) & (p < pb) & (ib >= 0)
    if k == 1:
        m &= e == 0
    elif k == 2:
        m &= e == 1
    a[k, m] = arr[ib[m]]
for k in range(3):
    P("   bedge class %d: sum a W = %.4f (lin C2 %.4f); sum a W q = %.4f (QED run %.4f)" % (k, (a[k] * W).sum(), nomb["linC2"][k], (a[k] * W * q).sum(), key["A"]["linC2"][k]))
# replicas: 2000 group-Poisson (own seed) and the 100 WP0 replicas (seed 0x5157A11)
gk = R["tree"].astype(np.int64) * (1 << 32) + R["grp"].astype(np.int64)
uq, inv = np.unique(gk, return_inverse=True)
rng = np.random.default_rng(0xBED6EC2)
K = rng.poisson(1.0, size=(2000, uq.size))[:, inv]
K0 = C.group_pois(0x5157A11, (uq >> 32), (uq & 0xFFFFFFFF), 100)[:, inv].astype(float)


def bedge_quant(c1, c2):
    red = ((1 - c1) * nomb["linC2"][1] + (1 - c2) * nomb["linC2"][2]) / (nomb["linT12"][1] + nomb["linT12"][2])
    return nomb["Ed"][0] / (nomb["Es"][0] * (1 - red))


bst = {}
for lab, KK in (("pois2000", K), ("wp0", K0)):
    cA = (KK * (a[0] * W * q)).sum(1) / (KK * (a[0] * W)).sum(1)
    c1 = (KK * (a[1] * W * q)).sum(1) / (KK * (a[1] * W)).sum(1)
    c2 = (KK * (a[2] * W * q)).sum(1) / (KK * (a[2] * W)).sum(1)
    ce = (KK * ((a[1] + a[2]) * W * q)).sum(1) / (KK * ((a[1] + a[2]) * W)).sum(1)
    bst[lab] = dict(all=sd(cA), E1=sd(c1), E2=sd(c2), E12=sd(ce), mean_all=float(cA.mean()))
    P("   bedge factor stat (%s): all %.4f, E1 %.4f, E2 %.4f, E1+E2 %.4f (mean all %.4f)" % (lab, sd(cA), sd(c1), sd(c2), sd(ce), cA.mean()))
c_b = dict(all=key["A"]["cAll"], E1=key["A"]["cE1"], E2=key["A"]["cE2"], E12=key["A"]["cE12"])
# lepton-efficiency sensitivity (bedge)
pe = np.sqrt(R["ep_px"] ** 2 + R["ep_py"] ** 2 + R["ep_pz"] ** 2) * 1e-3
pm = np.sqrt(R["em_px"] ** 2 + R["em_py"] ** 2 + R["em_pz"] ** 2) * 1e-3
ps = np.minimum(pe, pm)
sens_b = {}
for lab, pc, fc in (("p_soft < 0.10 x 0.8", 0.10, 0.8), ("p_soft < 0.10 x 1.2", 0.10, 1.2)):
    sw = np.where(ps < pc, fc, 1.0)
    sens_b[lab] = float((a[0] * W * sw * q).sum() / (a[0] * W * sw).sum())
    P("   bedge sensitivity %s: c %.4f (%+.4f)" % (lab, sens_b[lab], sens_b[lab] - c_b["all"]))
# kappa split in bedge with QED class factors (bedge Quant logic)
def bedge_kappa(c1, c2):
    Es1, Es2 = nomb["Es"][1], nomb["Es"][2]
    sC2 = [None, nomb["linC2"][1] / nomb["linT12"][1], nomb["linC2"][2] / nomb["linT12"][2]]
    sC3 = [None, nomb["linC3"][1] / nomb["linT12"][1], nomb["linC3"][2] / nomb["linT12"][2]]
    A1 = Es1 * (c1 * sC2[1] + 1 - sC2[1] - sC3[1]); B1 = Es1 * sC3[1]
    A2 = Es2 * (c2 * sC2[2] + 1 - sC2[2] - sC3[2]); B2 = Es2 * sC3[2]
    rho = nomb["Ed"][1] / nomb["Ed"][2]
    kap = (rho * A2 - A1) / (B1 - rho * B2)
    return nomb["Ed"][2] / (A2 + kap * B2), kap
def bedge_k48(c1, c2):
    t = nomb["linT12"][1] + nomb["linT12"][2]
    cc = (c1 * nomb["linC2"][1] + c2 * nomb["linC2"][2]) / (nomb["linC2"][1] + nomb["linC2"][2])
    s2 = (nomb["linC2"][1] + nomb["linC2"][2]) / t; s3 = (nomb["linC3"][1] + nomb["linC3"][2]) / t
    return nomb["Ed"][0] / (nomb["Es"][0] * (cc * s2 + 0.48 * s3 + 1 - s2 - s3))
kb_app = bedge_kappa(0.879, 0.882); kb_q = bedge_kappa(c_b["E1"], c_b["E2"])
k48_app = bedge_k48(0.879, 0.882); k48_q = bedge_k48(c_b["E1"], c_b["E2"])
P("   bedge E1/E2 split: applied a_D %.4f kappa %.4f; QED a_D %.4f kappa %.4f;  kappa = 0.48 reading: applied %.4f, QED %.4f (%+.2f %%)" %
  (*kb_app, *kb_q, k48_app, k48_q, 100 * (k48_q / k48_app - 1)))
fb_app = bedge_quant(0.879, 0.882); fb_q = bedge_quant(c_b["E1"], c_b["E2"])
P("   bedge f (bedge formula): applied %.5f, QED %.5f (%+.2f %%); direct re-evaluation: applied-constant %.5f, QED %.5f (%+.2f %%)" %
  (fb_app, fb_q, 100 * (fb_q / fb_app - 1), key["const_rec"]["fDir"], key["A"]["fDir"], 100 * (key["A"]["fDir"] / key["const_rec"]["fDir"] - 1)))
NUM["bedge"] = dict(c=c_b, c_applied=dict(all=key["const_rec"]["cAll"], E1=0.879, E2=0.882), stat=bst, sens=sens_b,
                    f_app=fb_app, f_qed=fb_q, fdir_app=key["const_rec"]["fDir"], fdir_qed=key["A"]["fDir"], kappa_app=kb_app, kappa_qed=kb_q,
                    k48_app=k48_app, k48_qed=k48_q,
                    variants={k: key[k]["cAll"] for k in ("A", "B", "Ac", "flat", "A_kloe")}, const_8855=key["const_8855"], const_edgetoy=key["const_edgetoy"])

# ---------------------------------------------------------------- edgeu
EU = json.load(open(O + "c2recon_edgeu.json"))
NUM["edgeu"] = EU
c_u = EU["factor"]["all"]["q_A"]["c"]
P("edgeu: c_u %.4f (E1 %.4f, E2 %.4f), stat %.4f; f applied %.5f -> QED %.5f (%+.2f %%)" %
  (c_u, EU["factor"]["E1"]["q_A"]["c"], EU["factor"]["E2"]["q_A"]["c"], EU["stat"]["all"]["rms"], EU["f"]["applied"], EU["f"]["qed"], 100 * EU["f"]["rel"]))

# ---------------------------------------------------------------- sigma, yields, ECAL
EN = json.load(open("/hdd3/przygoda/HADES/FAT_out/all_report/edgeu/edgeu_numbers.json"))
sig_u = EN["sigma"]; sig_u48 = EN["sigma_k48"]
rel_u = EU["f"]["qed"] / EU["f"]["applied"]; rel_u48 = EU["k48"]["qed"] / EU["k48"]["applied"]
BN = json.load(open("/hdd3/przygoda/HADES/FAT_out/all_report/bedge/bedge_numbers.json"))
sig_b = BN["SIGMA"]["fc_rec"]; sig_b48 = BN["SIGMA"]["a_k48"]
rel_b = fb_q / fb_app; rel_b48 = k48_q / k48_app
sig = dict(u=dict(app=sig_u["v"], qed=sig_u["v"] * rel_u, rel=rel_u - 1), u_k48=dict(app=sig_u48["v"], qed=sig_u48["v"] * rel_u48, rel=rel_u48 - 1),
           b=dict(app=sig_b["v"], qed=sig_b["v"] * rel_b, rel=rel_b - 1), b_k48=dict(app=sig_b48["v"], qed=sig_b48["v"] * rel_b48, rel=rel_b48 - 1))
for k, v in sig.items():
    P("   sigma(pp -> eta X) from the bump, %-6s: %.4f -> %.4f mb (%+.2f %%)" % (k, v["app"], v["qed"], 100 * v["rel"]))
NUM["sigma"] = sig
# double ratio peak / bump (u): DR_u = 0.429 / sigma_u
NUM["DR_u"] = dict(app=0.4291 / sig_u["v"], qed=0.4291 / (sig_u["v"] * rel_u))
P("   DR_u = 0.429 / sigma_u: %.3f -> %.3f" % (NUM["DR_u"]["app"], NUM["DR_u"]["qed"]))
EC = json.load(open("/hdd3/przygoda/HADES/FAT_out/all_report/ecalacc/ecalacc_numbers.json"))
bp = EC["best_pre"]
ec = {}
for lab, rel in (("bedge f (ecalacc input)", rel_b), ("u-edge f", (EU["f"]["qed"] / 0.4184))):
    val = bp["value"] * rel
    ratio = 3313.7094 / val
    er = ratio * bp["eratio"] / bp["ratio"]
    ec[lab] = dict(expected=val, ratio=ratio, eratio=er)
    P("   ECAL tag, best reading, with %s: expected %.0f (was %.0f), data/expected %.3f +- %.3f (was %.3f +- %.3f)" %
      (lab, val, bp["value"], ratio, er, bp["ratio"], bp["eratio"]))
# c_tag with the own weights (exact where the photon is known, circle otherwise): tagged rows in the m5c window (pre)
qt = np.where(Wt["gfound"] > 0, Wt["q_exact"], q)
mt = (R["tag"] > 0) & (R["m5c"] >= 0.49) & (R["m5c"] < 0.62)
ctag = float((W * qt)[mt].sum() / W[mt].sum())
ectag = float(math.sqrt((W[mt] ** 2 * (qt[mt] - ctag) ** 2).sum()) / W[mt].sum())
ctag_circ = float((W * q)[mt].sum() / W[mt].sum())
P("   ECAL tagged C2 rows (pre, compact tag, m5c 0.49-0.62): c_tag %.4f +- %.4f (exact where found, circle otherwise; ecalacc %.4f +- %.4f); circle only %.4f; %d rows, %d with the exact photon" %
  (ctag, ectag, EC["c_tag_pre"][0], EC["c_tag_pre"][1], ctag_circ, mt.sum(), (mt & (Wt["gfound"] > 0)).sum()))
NUM["ecal"] = dict(c_tag=ctag, ec_tag=ectag, c_tag_circle=ctag_circ, c_tag_ecalacc=EC["c_tag_pre"], best=ec, best_ref=bp)
ex = EC["expected"]["pre_all_asis_L0"]
P("   ECAL 'as is' reading (f uncorrected x T_sim uncorrected): expected %.0f, data/expected %.3f -- independent of the C2 factor" % (ex["value"], ex["ratio"]))

# ---------------------------------------------------------------- toy and the soft lepton
TY = json.load(open(O + "c2recon_toy.json"))
NUM["toy"] = TY
for pm_ in ("0.04", "0.06", "0.1", "0.15"):
    t = TY[pm_]
    P("   toy p_e > %s: true %.4f, circle %.4f (all); edge %.4f / %.4f; E1 %.4f, E2 %.4f" % (pm_, t["all"]["true"], t["all"]["A"], t["edge (u_true > -15 MeV)"]["true"],
                                                                                     t["edge (u_true > -15 MeV)"]["A"], t["E1"]["true"], t["E2"]["true"]))
# ---------------------------------------------------------------- posterior check in the full simulation (found photons)
lay = open(CW.ECAL_LAY).read().split("\n")
nr, nd = int(lay[0].split()[1]), int(lay[0].split()[3])
names = lay[1].split()[1:]
ea = np.fromfile(CW.ECAL_ROWS, dtype=np.float64).reshape(nr, nd)
E = {k: ea[:, i] for i, k in enumerate(names)}
mf = (E["cls"] == 2) & (E["tMatch"] > 0) & (E["fGeo"] > 0) & (E["ok"] > 0)
wf = E["w"][mf]
ein = (E["fGeoQed"] * E["meanQed"] / np.where(E["fGeo"] > 0, E["fGeo"], 1.0))[mf]
exq = E["tQed"][mf]
dd = exq - ein
md = float((wf * dd).sum() / wf.sum()); emd = float(math.sqrt((wf ** 2 * (dd - md) ** 2).sum()) / wf.sum())
NUM["posterior_check"] = dict(n=int(mf.sum()), exact=float((wf * exq).sum() / wf.sum()), restricted=float((wf * ein).sum() / wf.sum()),
                              unrestricted=float((wf * E["meanQed"][mf]).sum() / wf.sum()), diff=md, ediff=emd, corr=float(np.corrcoef(exq, ein)[0, 1]))
P("   posterior check, %d C2 rows whose photon the ECAL simulation found: exact <w> %.4f; circle restricted to the active modules %.4f "
  "(difference %+.4f +- %.4f); unrestricted circle %.4f; per-row correlation %.2f" % (NUM["posterior_check"]["n"], NUM["posterior_check"]["exact"],
                                                                                    NUM["posterior_check"]["restricted"], md, emd,
                                                                                    NUM["posterior_check"]["unrestricted"], NUM["posterior_check"]["corr"]))
json.dump(NUM, open(O + "c2recon_numbers.json", "w"), indent=1, default=float)

# ---------------------------------------------------------------- the factor table
D = EU["definitions"]
Lt = []
Lt.append("c2recon: the C2 gamma* factor eps_QED / eps_iso, by estimator and definition (simulated C2 rows of the standard set, S0, mode A, W = w x RowFactor(108))")
Lt.append("event-by-event QED weight = conditional mean along the photon circle (SMASH density A); stat = rms of 2000 group-Poisson replicas")
Lt.append("estimator / definition | all | E1 | E2 | stat (all) | note")
Lt.append("u-edge estimator (edgeu, global fit, linear h with the data covariance) | %.4f | %.4f | %.4f | %.4f | applied 0.879 / 0.8825 -> effective %.4f" %
          (c_u, EU["factor"]["E1"]["q_A"]["c"], EU["factor"]["E2"]["q_A"]["c"], EU["stat"]["all"]["rms"],
           (0.879 * EU["factor"]["E1"]["q_A"]["hiso"] + 0.8825 * EU["factor"]["E2"]["q_A"]["hiso"]) / EU["factor"]["all"]["q_A"]["hiso"]))
Lt.append("u-edge class fits (E1, E2 fitted separately) | -- | %.4f | %.4f | -- | " % (EU["classes"]["E1"]["cq"], EU["classes"]["E2"]["cq"]))
Lt.append("per-slice edge estimator (bedge, 4 slices, linear h with the data covariance) | %.4f | %.4f | %.4f | %.4f | applied 0.879 / 0.882 -> effective %.4f (quoted as 0.8855)" %
          (c_b["all"], c_b["E1"], c_b["E2"], bst["pois2000"]["all"], key["const_rec"]["cAll"]))
for nm in ("all C2 rows (S0, all M(pipi), all M4)", "0.30 <= M(pipi) < 0.40, all M4", "window B (0.4550 <= M4 < 0.5225)",
           "edge window -15...+5 MeV in u (ecalacc)", "fit range -45...+45 MeV, 0.30-0.40", "far below u < -60 MeV"):
    Lt.append("window mean: %s | %.4f | %.4f | %.4f | -- | " % (nm, D[nm + " | all"]["c"], D[nm + " | E1"]["c"], D[nm + " | E2"]["c"]))
Lt.append("ECAL-tagged C2 rows (pre, m5c 0.49-0.62; exact photon where found) | %.4f | -- | -- | %.4f | ecalacc 0.876 +- 0.009" % (ctag, ectag))
Lt.append("toy (WP8-like acceptance, p_e > 0.06), all / edge | %.4f / %.4f | %.4f | %.4f | -- | WP8 toy 0.880 (0.889 / 0.868)" %
          (TY["0.06"]["all"]["true"], TY["0.06"]["edge (u_true > -15 MeV)"]["true"], TY["0.06"]["E1"]["true"], TY["0.06"]["E2"]["true"]))
Lt.append("toy, p_e > 0.10, all / edge | %.4f / %.4f | %.4f | %.4f | -- | bedge edge toy at p_e > 0.10: 0.859" %
          (TY["0.1"]["all"]["true"], TY["0.1"]["edge (u_true > -15 MeV)"]["true"], TY["0.1"]["E1"]["true"], TY["0.1"]["E2"]["true"]))
open(O + "c2recon_table.txt", "w").write("\n".join(Lt) + "\n")
P("\n" + "\n".join(Lt))
