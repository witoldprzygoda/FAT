#!/usr/bin/env python3
# vflex_ties.py -- independent verification of package flex: attempts to refute the plain-words conclusions
#   (1) "all eta channels are below SMASH by similar factors, one common eta scale": F8 with tied eta scales
#   (2) "C5/C2 shows no deviation from the PDG BR": F8 with s_C5 = s_C2
#   (3) "G_eta is free / degenerate": F8 with s_Geta = s_C5, s_Geta = s_C2, s_Geta = 1
#   (4) background: s_GC = s_GD_nom, s_GD_om = s_GD_nom
#   Every tied fit: the same estimator, hybrid GLS (data + BB + replica covariance of the smoothed part at the fitted theta),
#   delta chi2 against F8; Wald tests from the F8 total covariance (stat + MC); window balance; information shares.
# Run: python3 -B plots/all_report/vflex_ties.py
import sys, json
import numpy as np
sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import vflex_core as V


class Tied:
    def __init__(self, base, groups, fixed=None, tag=""):
        """groups: list of lists of base parameter indices that share one parameter; fixed: {base index: value}"""
        self.base = base; self.fixedv = fixed or {}
        n = base.np
        cols = []; used = set()
        for g in groups:
            c = np.zeros(n); c[g] = 1; cols.append(c); used |= set(g)
        for i in range(n):
            if i not in used and i not in self.fixedv:
                c = np.zeros(n); c[i] = 1; cols.append(c)
        self.A = np.array(cols).T
        self.off = np.zeros(n)
        for i, v in self.fixedv.items():
            self.off[i] = v
        self.np = self.A.shape[1]
        self.names = ["+".join(base.names[j] for j in np.where(self.A[:, q] > 0)[0]) for q in range(self.np)]
        self.start = np.linalg.lstsq(self.A, base.start - self.off, rcond=None)[0]
        kq = int(np.where(self.A[base.np - 1] > 0)[0][0])
        self.cons = {kq: base.cons[base.np - 1]}
        self.nscale = self.np - 1
        self.comps, self.bg = base.comps, base.bg

    def full(self, th):
        return self.A @ th + self.off

    def eval(self, th, L, T, Vt):
        mu, J, vbb, k, bg = self.base.eval(self.full(th), L, T, Vt)
        return mu, J @ self.A, vbb, k, bg


D = V.load(); L = V.Layout(D, (1,))
M = V.FlexModel(["C2", "C3", "C4", "C5", "Geta", "GD_om", "GD_nom", "GC"])
T, Vt = V.templates(D, L, M)
TR = [V.templates(D, L, M, rep=r)[0] for r in range(V.NREPS)]
N = {n: i for i, n in enumerate(M.names)}
R8 = V.fit(M, L, T, Vt, egcov=L.eg)
g8 = V.hybrid_gls(M, L, T, Vt, R8, TR)
out = {"F8": dict(th=R8['th'].tolist(), gls=g8[0], glsd=g8[1], pearson=R8['pearson'])}
# MC covariance of F8 for the Wald tests
THr = V.replica_fits(M, L, TR, Vt, R8['th'])
Cmc = np.cov(THr.T, ddof=1); Ct = R8['Cstat'] + Cmc
print(f"F8: hybrid GLS {g8[0]:.2f}, data-only {g8[1]:.1f}")
tests = {
    "eta C2=C3=C4=C5": ([[N['s_C2'], N['s_C3'], N['s_C4'], N['s_C5']]], None),
    "eta C2=C4=C5 (C3 free)": ([[N['s_C2'], N['s_C4'], N['s_C5']]], None),
    "eta all five (C2..C5, Geta)": ([[N['s_C2'], N['s_C3'], N['s_C4'], N['s_C5'], N['s_Geta']]], None),
    "C5=C2": ([[N['s_C2'], N['s_C5']]], None),
    "Geta=C5": ([[N['s_C5'], N['s_Geta']]], None),
    "Geta=C2": ([[N['s_C2'], N['s_Geta']]], None),
    "Geta=1": ([], {N['s_Geta']: 1.0}),
    "Geta=0.5": ([], {N['s_Geta']: 0.5}),
    "GC=GD_nom": ([[N['s_GD_nom'], N['s_GC']]], None),
    "GD_om=GD_nom": ([[N['s_GD_om'], N['s_GD_nom']]], None),
    "GC=1": ([], {N['s_GC']: 1.0}),
}
xc = 0.5 * (L.lo + L.hi)
for tag, (groups, fixed) in tests.items():
    Mt = Tied(M, groups, fixed, tag)
    R = V.fit(Mt, L, T, Vt, egcov=L.eg)
    g = V.hybrid_gls(Mt, L, T, Vt, R, TR)
    thf = Mt.full(R['th'])
    ddof = M.np - Mt.np
    s = (xc > 0.565) & (xc < 0.70)
    nonEta = sum(thf[N[n]] * T[s, N[n]].sum() for n in ("s_GD_om", "s_GD_nom", "s_GC")) / sum(T[s, N[n]].sum() for n in ("s_GD_om", "s_GD_nom", "s_GC"))
    print(f"[{tag}] ok {R['ok']} dof {ddof}: hybrid GLS {g[0]:.2f} (d {g[0] - g8[0]:+.2f}), data-only {g[1]:.1f} (d {g[1] - g8[1]:+.1f}), Pearson {R['pearson']:.1f}; "
          f"model/data fit range {R['mu'].sum() / L.y.sum() - 1:+.4f}; non-eta CR-U scale {nonEta:.3f}")
    print("     " + "  ".join(f"{n}={thf[i]:.3f}" for i, n in enumerate(M.names)))
    out[tag] = dict(ok=R['ok'], dof=ddof, gls=g[0], dgls=g[0] - g8[0], glsd=g[1], th_full=thf.tolist(), nonEta_CRU=nonEta)


# Wald tests from the F8 total covariance
def wald_common(idx, C=Ct, th=R8['th']):
    s = th[idx]; Ci = np.linalg.inv(C[np.ix_(idx, idx)]); one = np.ones(len(idx))
    a = (one @ Ci @ s) / (one @ Ci @ one)
    r = s - a
    return a, 1 / np.sqrt(one @ Ci @ one), r @ Ci @ r, len(idx) - 1


from scipy.stats import chi2 as CH
for tag, names in {"C2..C5": ["s_C2", "s_C3", "s_C4", "s_C5"], "C2,C4,C5": ["s_C2", "s_C4", "s_C5"], "C2..C5,Geta": ["s_C2", "s_C3", "s_C4", "s_C5", "s_Geta"],
                   "C2,C5": ["s_C2", "s_C5"], "C5,Geta": ["s_C5", "s_Geta"], "C2,Geta": ["s_C2", "s_Geta"]}.items():
    idx = [N[n] for n in names]
    a, ea, c2, nd = wald_common(idx)
    print(f"Wald (stat+MC) common scale of {tag}: {a:.3f} +- {ea:.3f}, chi2 {c2:.2f} / {nd} (p {CH.sf(c2, nd):.3g})")
    out["wald " + tag] = dict(a=a, ea=ea, chi2=c2, ndf=nd, p=CH.sf(c2, nd))
# Wald with the attenuation-corrected values (divide by 1 + RF bias of the independent toys, scenario 0)
z = np.load(V.DIR + "/vflex_toys_s0.npz")
att = 1 + (z['th'][z['ok'][:, 0], 0].mean(0) - z['truth']) / z['truth']
thc = R8['th'] / att
Cc = Ct / np.outer(att, att)
for tag, names in {"C2,C4,C5 (att.-corr.)": ["s_C2", "s_C4", "s_C5"], "C2,C5 (att.-corr.)": ["s_C2", "s_C5"], "C5,Geta (att.-corr.)": ["s_C5", "s_Geta"],
                   "C2,C4,C5,Geta (att.-corr.)": ["s_C2", "s_C4", "s_C5", "s_Geta"]}.items():
    idx = [N[n] for n in names]
    a, ea, c2, nd = wald_common(idx, Cc, thc)
    print(f"Wald {tag}: {a:.3f} +- {ea:.3f}, chi2 {c2:.2f} / {nd} (p {CH.sf(c2, nd):.3g})")
    out["wald " + tag] = dict(a=a, ea=ea, chi2=c2, ndf=nd, p=CH.sf(c2, nd))
print("attenuation factors (1 + RF bias):", dict(zip(M.names, np.round(att, 3))))
# ratio C5/C2 with delta method (stat+MC covariance)
i2, i5 = N['s_C2'], N['s_C5']
rr = R8['th'][i5] / R8['th'][i2]
gr = np.zeros(M.np); gr[i5] = 1 / R8['th'][i2]; gr[i2] = -R8['th'][i5] / R8['th'][i2] ** 2
print(f"s_C5/s_C2 = {rr:.3f} +- {np.sqrt(gr @ Ct @ gr):.3f} (stat {np.sqrt(gr @ R8['Cstat'] @ gr):.3f}, MC {np.sqrt(gr @ Cmc @ gr):.3f}); "
      f"replica spread of the ratio {np.std(THr[:, i5] / THr[:, i2], ddof=1):.3f}")
# information shares (Fisher T^2 / V at the F8 fit) by M4 region
Vv = R8['V']
regs = dict(low=(0.40, 0.455), B=(0.455, 0.5225), gap=(0.5225, 0.5336), P=(0.5336, 0.5542), gap2=(0.5542, 0.565), CRU=(0.565, 0.70))
print("information shares T^2/V by region:")
for k, n in enumerate(M.comps):
    I = T[:, k] ** 2 / Vv
    sh = {r: I[(xc > a) & (xc < b)].sum() / I.sum() for r, (a, b) in regs.items()}
    sl = [I[L.j == j].sum() / I.sum() for j in range(5)]
    print(f"   {n:7s} " + " ".join(f"{r} {100 * v:5.1f}%" for r, v in sh.items()) + " | slices " + " ".join(f"{100 * v:4.1f}" for v in sl))
json.dump(out, open(V.DIR + "/vflex_ties.json", "w"), indent=1, default=float)
