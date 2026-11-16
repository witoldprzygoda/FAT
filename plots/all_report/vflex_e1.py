#!/usr/bin/env python3
# vflex_e1.py -- independent verification of package flex: the variant V_E1k (E1 + E2, k_E1 free, k_E2 = 1 +- 0.05, model F8)
#   and the same fit with k_E1 fixed to 1 (to show how much the conversion scales depend on the E1 CB model)
# Run: python3 -B plots/all_report/vflex_e1.py
import sys, json
import numpy as np
sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import vflex_core as V


class FlexE(V.FlexModel):
    """F-type model on E1 + E2 with one CB scale per class: k_E1 free (or fixed), k_E2 constrained"""
    def __init__(self, comps, kE1fixed=None):
        super().__init__(comps)
        m = len(comps)
        self.names = ["s_" + c for c in comps] + ["k_E1", "k_E2"]
        self.np = m + 2
        self.start = np.array([1.0] * m + [1.0, 1.0])
        self.cons = {m + 1: (1.0, 0.05)}
        if kE1fixed is not None:
            self.cons[m] = (kE1fixed, 1e-6)
        self.nscale = m

    def eval(self, th, L, T, Vt):
        m = len(self.comps)
        s = th[:m]
        k = np.where(L.e == 0, th[m], th[m + 1])
        mu = k * L.cb + T @ s
        J = np.zeros((L.n, self.np)); J[:, :m] = T
        J[:, m] = np.where(L.e == 0, L.cb, 0); J[:, m + 1] = np.where(L.e == 1, L.cb, 0)
        vbb = (s ** 2 * Vt).sum(1)
        bgpart = (s * T * self.bg).sum(1)
        return mu, J, vbb, k, bgpart


D = V.load()
L = V.Layout(D, (0, 1))
print(f"E1 + E2 bins {L.n}, N {L.y.sum():.0f} (E1 {L.y[L.e == 0].sum():.0f}), CB E1 {L.cb[L.e == 0].sum():.0f}")
comps = ["C2", "C3", "C4", "C5", "Geta", "GD_om", "GD_nom", "GC"]
out = {}
for tag, kf in (("V_E1k", None), ("V_E1k1 (k_E1 = 1)", 1.0)):
    M = FlexE(comps, kf)
    T, Vt = V.templates(D, L, M)
    R = V.fit(M, L, T, Vt, egcov=L.eg)
    TR = [V.templates(D, L, M, rep=r)[0] for r in range(V.NREPS)]
    TH = V.replica_fits(M, L, TR, Vt, R['th'])
    Cmc = np.cov(TH.T, ddof=1)
    g = V.hybrid_gls(M, L, T, Vt, R, TR)
    st = np.sqrt(np.diag(R['Cstat'])); mc = np.sqrt(np.diag(Cmc)); tot = np.sqrt(st ** 2 + mc ** 2)
    print(f"[{tag}] ok {R['ok']} Pearson {R['pearson']:.1f} hybrid GLS {g[0]:.1f} data-only {g[1]:.1f}")
    for i, n in enumerate(M.names):
        print(f"   {n:9s} {R['th'][i]:8.4f} stat {st[i]:.4f} MC {mc[i]:.4f} tot {tot[i]:.4f}")
    out[tag] = dict(names=M.names, th=R['th'].tolist(), stat=st.tolist(), mc=mc.tolist(), gls=g[0], pearson=R['pearson'])
    th = R['th']
    for a, b in ((3, 0), (2, 0), (1, 0), (2, 1)):
        r = th[a] / th[b]
        gr = np.zeros(M.np); gr[a] = 1 / th[b]; gr[b] = -th[a] / th[b] ** 2
        es = np.sqrt(gr @ R['Cstat'] @ gr); em = np.std(TH[:, a] / TH[:, b], ddof=1)
        print(f"   ratio {M.names[a]}/{M.names[b]} = {r:.3f} +- {es:.3f} (stat) +- {em:.3f} (MC replica spread) -> total {np.hypot(es, em):.3f}")
        out[tag][f"ratio {M.names[a]}/{M.names[b]}"] = (r, es, em)
json.dump(out, open(V.DIR + "/vflex_e1.json", "w"), indent=1, default=float)
