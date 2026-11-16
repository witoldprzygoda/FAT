"""vbedge_ana.py -- INDEPENDENT VERIFIER of the package bedge: decisive numbers with own code.

Written 2026-11-15 (prefix vbedge). Uses only vbedge_core.py and the arrays of vbedge_fill.C.
  1. nominal f_eta(edge) = E_data / E_sim (all M(ee), E1, E2)
  2. statistical error: own data bootstrap (200 replicas, own generator)
  3. MC error: own simulation replicas (100, own generator)
  4. variants: r = 0.89 (simulation), lower start 0.29 (slices), smearing matched (sim 6.1 MeV), CB = N++ + N--
  5. linear decomposition of the simulated step (data covariance), C2 angle correction, kappa = 0.48 reading,
     kappa free (E1/E2)
  6. visible yield and sigma(pp -> eta X)
Output: /hdd3/przygoda/HADES/FAT_out/all_report/vbedge/vbedge_ana.json and stdout.
Run: python3 -I plots/all_report/vbedge_ana.py   (from /home/przygoda/HADES/FAT/FAT)
"""
import json
import sys

import numpy as np

sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
from vbedge_core import *  # noqa

d = load()
N, S, R, B = d["N"], d["S"], d["R"], d["B"]
CF = {"toy": (0.880, 0.889, 0.868), "clos": (0.891, 0.869, 0.896)}
CF["rec"] = tuple(0.5 * (a + b) for a, b in zip(CF["toy"], CF["clos"]))
KAP = 0.48
PT3 = 64.0 / 63.0
L = 6079.2


def simy(arr, cls, pa, pbb):
    return arr[:, pa:pbb].sum(axis=(0, 1)) if cls < 0 else arr[cls, pa:pbb].sum(0)


def evaluate(Ndata, sset, scfg=0, s_sim=PEAKS, cbopt=None, Sarr=None, need_lin=True):
    """returns dict per class k (0 all, 1 E1, 2 E2): Ed, vEd, Es, lin (T12, C2, C3), Vd list, cfgs"""
    Sarr = S if Sarr is None else Sarr
    out = []
    for k in range(3):
        cls = k - 1
        Ed = vEd = Es = 0.0; lin = np.zeros(3); Vd = []; cfgs = []; chi2 = 0; ndf = 0
        for lo, hi in sset:
            pa, pbb = pb(lo), pb(hi)
            y, cb, cov, u = data_spec(Ndata, cls, pa, pbb, cbopt)
            cD = Cfg(lo, hi)
            cS = Cfg(lo, hi, s=s_sim) if s_sim != PEAKS else cD
            h, eh, c2, nd, V, p = fit(y, cb, cov, cD)
            ys = simy(Sarr[scfg, 2], cls, pa, pbb)
            hs = fit(ys, None, None, cS)[0]
            if need_lin:
                for g, gi in enumerate((2, 0, 1)):
                    lin[g] += fit(simy(Sarr[scfg, gi], cls, pa, pbb), None, None, cS, Vfix=V)[0]
            Ed += h; vEd += eh * eh; Es += hs; Vd.append(V); cfgs.append((cD, cS)); chi2 += c2; ndf += nd
        out.append(dict(Ed=Ed, vEd=vEd, Es=Es, lin=lin, Vd=Vd, cfgs=cfgs, chi2=chi2, ndf=ndf))
    return out


def quant(Ed, Es, lin):
    """Ed, Es: [3]; lin: [3][3] (class; T12, C2, C3). Returns dict of the quantities."""
    q = {}
    q["f_all"], q["f_E1"], q["f_E2"] = Ed[0] / Es[0], Ed[1] / Es[1], Ed[2] / Es[2]
    t = lin[1][0] + lin[2][0]
    for tag in ("toy", "clos", "rec"):
        c = CF[tag]
        red = ((1 - c[1]) * lin[1][1] + (1 - c[2]) * lin[2][1]) / t
        q["fc_" + tag] = Ed[0] / (Es[0] * (1 - red))
    c = CF["rec"]
    q["fc_E1"] = Ed[1] / (Es[1] * (1 - (1 - c[1]) * lin[1][1] / lin[1][0]))
    q["fc_E2"] = Ed[2] / (Es[2] * (1 - (1 - c[2]) * lin[2][1] / lin[2][0]))
    s2 = (lin[1][1] + lin[2][1]) / t; s3 = (lin[1][2] + lin[2][2]) / t
    c2w = (c[1] * lin[1][1] + c[2] * lin[2][1]) / (lin[1][1] + lin[2][1])
    q["a_k48"] = Ed[0] / (Es[0] * (c2w * s2 + KAP * s3 + 1 - s2 - s3))
    q["a_k48_unc"] = Ed[0] / (Es[0] * (s2 + KAP * s3 + 1 - s2 - s3))
    # kappa free: Ed_i = a Es_i (c_i s2_i + kappa s3_i + 1 - s2_i - s3_i), i = E1, E2
    for corr in (0, 1):
        c1, cc2 = (c[1], c[2]) if corr == 0 else (1.0, 1.0)
        A, Bv = [], []
        for i, ci in ((1, c1), (2, cc2)):
            s2i, s3i = lin[i][1] / lin[i][0], lin[i][2] / lin[i][0]
            A.append(Es[i] * (ci * s2i + 1 - s2i - s3i)); Bv.append(Es[i] * s3i)
        rho = Ed[1] / Ed[2]
        kap = (rho * A[1] - A[0]) / (Bv[0] - rho * Bv[1])
        a = Ed[2] / (A[1] + kap * Bv[1])
        q["a_D" + ("" if corr == 0 else "_unc")] = a; q["kappa" + ("" if corr == 0 else "_unc")] = kap
    return q


def arr(res, key):
    return np.array([r[key] for r in res])


def h68(a):
    a = np.sort(np.asarray(a))
    return 0.5 * (np.quantile(a, 0.84134) - np.quantile(a, 0.15866))


print("=" * 100)
nom = evaluate(N, NOMSET)
Ed0, Es0 = arr(nom, "Ed"), arr(nom, "Es")
lin0 = [r["lin"] for r in nom]
q0 = quant(Ed0, Es0, lin0)
for k, r in enumerate(nom):
    print(f"class {k}: E_data {r['Ed']:.4f} +- {r['vEd']**.5:.4f}  E_sim {r['Es']:.4f}  f {r['Ed']/r['Es']:.5f}  chi2 {r['chi2']:.1f}/{r['ndf']}"
          f"  lin T12 {r['lin'][0]:.1f} C2 {r['lin'][1]:.1f} C3 {r['lin'][2]:.1f} rest {r['lin'][0]-r['lin'][1]-r['lin'][2]:.1f}")
print("C2 share (all) %.4f, C3 share %.4f" % (nom[0]["lin"][1] / nom[0]["lin"][0], nom[0]["lin"][2] / nom[0]["lin"][0]))
print("quantities:", {k: round(v, 5) for k, v in q0.items()})

# ------------------------------------------------------------------------------------------- stat: bootstrap
qb = []
Edb = np.zeros((NRD, 3))
for b in range(NRD):
    Nb = B[b]
    for k in range(3):
        cls = k - 1; e = 0
        for (lo, hi), (cD, cS) in zip(NOMSET, nom[k]["cfgs"]):
            y, cb, cov, u = data_spec(Nb, cls, pb(lo), pb(hi))
            e += fit(y, cb, cov, cD)[0]
        Edb[b, k] = e
    qb.append(quant(Edb[b], Es0, lin0))
# ------------------------------------------------------------------------------------------- MC: simulation replicas
qr = []
Esr = np.zeros((NRS, 3)); linr = np.zeros((NRS, 3, 3))
for r in range(NRS):
    Sr = R[r]   # [cfg][grp C2 C3 T12][cls][pp][m4]
    for k in range(3):
        cls = k - 1
        for i, ((lo, hi), (cD, cS)) in enumerate(zip(NOMSET, nom[k]["cfgs"])):
            pa, pbb = pb(lo), pb(hi)
            Esr[r, k] += fit(simy(Sr[0, 2], cls, pa, pbb), None, None, cS)[0]
            for g, gi in enumerate((2, 0, 1)):
                linr[r, k, g] += fit(simy(Sr[0, gi], cls, pa, pbb), None, None, cS, Vfix=nom[k]["Vd"][i])[0]
    qr.append(quant(Ed0, Esr[r], linr[r]))
print("\nquantity       value    stat_fit  stat_boot (68)   boot_mean   MC_rms  (68)   rep_mean")
statfit = {}
for key in q0:
    # fit-error propagation of the independent class E_data errors
    v = 0
    for k in range(3):
        Ed = Ed0.copy(); Ed[k] += nom[k]["vEd"] ** .5
        v += (quant(Ed, Es0, lin0)[key] - q0[key]) ** 2
    statfit[key] = v ** .5
    bb = np.array([x[key] for x in qb]); rr = np.array([x[key] for x in qr])
    print(f"{key:12s} {q0[key]:9.5f} {statfit[key]:9.5f} {bb.std(ddof=1):9.5f} ({h68(bb):7.5f}) {bb.mean():9.5f} {rr.std(ddof=1):9.5f} ({h68(rr):7.5f}) {rr.mean():9.5f}")
print("E_sim replicas (all): rms %.2f mean %.2f" % (Esr[:, 0].std(ddof=1), Esr[:, 0].mean()))
print("bootstrap/fit ratio (f_all): %.3f" % (np.std([x['f_all'] for x in qb], ddof=1) / statfit['f_all']))

# ------------------------------------------------------------------------------------------- variants
res = {"nominal": q0}
def shift_report(tag, qv, qbv=None, qrv=None):
    D = qv["f_all"] - q0["f_all"]
    nb = np.std([a["f_all"] - b["f_all"] for a, b in zip(qbv, qb)], ddof=1) if qbv else 0.0
    nr = np.std([a["f_all"] - b["f_all"] for a, b in zip(qrv, qr)], ddof=1) if qrv else 0.0
    nr68 = h68([a["f_all"] - b["f_all"] for a, b in zip(qrv, qr)]) if qrv else 0.0
    print(f"variant {tag:14s}: f {qv['f_all']:.5f}  fc_rec {qv['fc_rec']:.5f}  Delta {D:+.5f}  noise boot {nb:.5f}  noise rep {nr:.5f} (68 {nr68:.5f})")
    res[tag] = dict(f=qv["f_all"], fc=qv["fc_rec"], D=D, nb=nb, nr=nr, nr68=nr68)

# (a) r = 0.89: simulation configuration 1
v_r = evaluate(N, NOMSET, scfg=1)
q_r = quant(arr(v_r, "Ed"), arr(v_r, "Es"), [r["lin"] for r in v_r])
qr_r = []
for r in range(NRS):
    Sr = R[r]; Es = np.zeros(3); lin = np.zeros((3, 3))
    for k in range(3):
        cls = k - 1
        for i, ((lo, hi), (cD, cS)) in enumerate(zip(NOMSET, v_r[k]["cfgs"])):
            pa, pbb = pb(lo), pb(hi)
            Es[k] += fit(simy(Sr[1, 2], cls, pa, pbb), None, None, cS)[0]
            for g, gi in enumerate((2, 0, 1)):
                lin[k, g] += fit(simy(Sr[1, gi], cls, pa, pbb), None, None, cS, Vfix=v_r[k]["Vd"][i])[0]
    qr_r.append(quant(arr(v_r, "Ed"), Es, lin))
shift_report("r089", q_r, None, qr_r)

# (b) lower start 0.29 (slices; data and simulation change)
SET29 = [(0.29, 0.34), (0.34, 0.37), (0.37, 0.39), (0.39, 0.40)]
v_l = evaluate(N, SET29)
q_l = quant(arr(v_l, "Ed"), arr(v_l, "Es"), [r["lin"] for r in v_l])
qb_l, qr_l = [], []
Es_l, lin_l = arr(v_l, "Es"), [r["lin"] for r in v_l]
for b in range(NRD):
    Ed = np.zeros(3)
    for k in range(3):
        for (lo, hi), (cD, cS) in zip(SET29, v_l[k]["cfgs"]):
            y, cb, cov, u = data_spec(B[b], k - 1, pb(lo), pb(hi))
            Ed[k] += fit(y, cb, cov, cD)[0]
    qb_l.append(quant(Ed, Es_l, lin_l))
for r in range(NRS):
    Sr = R[r]; Es = np.zeros(3); lin = np.zeros((3, 3))
    for k in range(3):
        for i, ((lo, hi), (cD, cS)) in enumerate(zip(SET29, v_l[k]["cfgs"])):
            pa, pbb = pb(lo), pb(hi)
            Es[k] += fit(simy(Sr[0, 2], k - 1, pa, pbb), None, None, cS)[0]
            for g, gi in enumerate((2, 0, 1)):
                lin[k, g] += fit(simy(Sr[0, gi], k - 1, pa, pbb), None, None, cS, Vfix=v_l[k]["Vd"][i])[0]
    qr_l.append(quant(arr(v_l, "Ed"), Es, lin))
shift_report("lowstart029", q_l, qb_l, qr_l)

# (c) smearing matched: simulation fitted with 6.1 MeV (its own peak width), data 5.4 MeV
v_s = evaluate(N, NOMSET, s_sim=0.0061)
q_s = quant(arr(v_s, "Ed"), arr(v_s, "Es"), [r["lin"] for r in v_s])
shift_report("smearMatch", q_s)
# (d) CB = N++ + N-- (data only)
v_c = evaluate(N, NOMSET, cbopt={"mode": 1}, need_lin=True)
q_c = quant(arr(v_c, "Ed"), arr(v_c, "Es"), [r["lin"] for r in v_c])
shift_report("cbSum", q_c)

# ------------------------------------------------------------------------------------------- yields and sigma
W, WR = d["W"], d["WR"]       # [grp C2 C3][cls][win B, any]
c = CF["rec"]
out = {}
for wi, wn in ((0, "B"), (1, "allM4")):
    C2 = W[0, :, wi]; C3 = W[1, :, wi]
    C2r = WR[:, 0, :, wi]; C3r = WR[:, 1, :, wi]
    T = (C2 + C3).sum(); Tc = c[1] * C2[0] + c[2] * C2[1] + C2[2] + C3.sum()
    Tk = c[1] * C2[0] + c[2] * C2[1] + C2[2] + KAP * C3.sum()
    Y1 = q0["f_all"] * T * PT3; Y2 = q0["fc_rec"] * Tc * PT3; Y3 = q0["a_k48"] * Tk * PT3
    # MC of Y1 from the same replicas (E_sim and the window sum together)
    y1r = np.array([qr[r]["f_all"] * (C2r[r].sum() + C3r[r].sum()) * PT3 for r in range(NRS)])
    y1b = np.array([qb[b]["f_all"] * T * PT3 for b in range(NRD)])
    E2 = q0["f_all"] * (C2[1] + C3[1]) * PT3
    print(f"window {wn}: SMASH C2+C3 = {T:.1f} raw (x 64/63 = {T*PT3:.1f}); Y(f) = {Y1:.0f} +- {y1b.std(ddof=1):.0f} (stat) +- {y1r.std(ddof=1):.0f} (MC) rows = {Y1/L:.3f} nb;"
          f" Y(fc) = {Y2:.0f}; Y(kappa 0.48) = {Y3:.0f}; E2 part f x (C2+C3)_E2 = {E2:.0f}")
    out[wn] = dict(T=T, Y1=Y1, Y2=Y2, Y3=Y3, stat=y1b.std(ddof=1), mc=y1r.std(ddof=1))
# sigma
K = 1.0694 * 1.0104 * 0.97 / 1.00436
st_b = {k: np.std([x[k] for x in qb], ddof=1) for k in q0}
mc_r = {k: np.std([x[k] for x in qr], ddof=1) for k in q0}
for key in ("f_all", "fc_rec", "a_k48", "a_k48_unc"):
    print(f"sigma({key}) = {q0[key]*K:.4f} +- {st_b[key]*K:.4f} (stat) +- {mc_r[key]*K:.4f} (MC) mb  [scale {q0[key]:.5f}; PDG 2026 scale {q0[key]/1.00436:.4f}]")
normup = (6.9**2 + 5.0**2 + 3.17**2 + 3.1**2) ** .5; normdn = (6.9**2 + 5.0**2 + 3.17**2 + 4.1**2) ** .5
print(f"norm: +{normup:.2f} % / -{normdn:.2f} %  -> fc_rec: +{q0['fc_rec']*K*normup/100:.4f} / -{q0['fc_rec']*K*normdn/100:.4f} mb")
json.dump(dict(q0=q0, statfit=statfit, stat_boot=st_b, mc=mc_r,
               mc68={k: h68([x[k] for x in qr]) for k in q0}, repmean={k: float(np.mean([x[k] for x in qr])) for k in q0},
               variants=res, yields=out, Esr=Esr[:, 0].tolist(), fb=[x["f_all"] for x in qb], fr=[x["f_all"] for x in qr],
               nom=[dict(Ed=r["Ed"], eEd=r["vEd"] ** .5, Es=r["Es"], lin=list(r["lin"])) for r in nom]),
          open(f"{OUT}/vbedge_ana.json", "w"), indent=1, default=float)
print("written", f"{OUT}/vbedge_ana.json")
