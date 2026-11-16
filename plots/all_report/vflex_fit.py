#!/usr/bin/env python3
# ============================================================================
# vflex_fit.py -- independent verification of package flex: Part III anchor, SMASH as is, nominal F8, alternative F7bj
# ============================================================================
# Run (from /home/przygoda/HADES/FAT/FAT):  python3 -B plots/all_report/vflex_fit.py
# Output: /hdd3/przygoda/HADES/FAT_out/all_report/vflex/vflex_fit.json, vflex_fit.log
# ============================================================================
import json, sys, time
import numpy as np
sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import vflex_core as V

OUT = V.DIR
log = open(OUT + "/vflex_fit.log", "w")


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s); log.write(s + "\n"); log.flush()


t0 = time.time()
D = V.load()
L = V.Layout(D, (1,))
P(f"bins {L.n}; data N {L.y.sum():.0f}, CB {L.cb.sum():.1f}; event-group var/N {np.trace(L.eg) / L.y.sum():.4f}")
res = {}


def run(tag, model, rawbg=False, reps=True):
    T, Vt = V.templates(D, L, model, rawbg=rawbg)
    R = V.fit(model, L, T, Vt, egcov=L.eg)
    out = dict(names=model.names, th=R['th'].tolist(), ok=R['ok'], pearson=R['pearson'], ls_dmax=float(np.max(np.abs(R['ls_th'] - R['th'][R['free']]))))
    P(f"[{tag}] least-squares root vs polished fixed point: max |d theta| = {out['ls_dmax']:.2e} (nfev {R['ls_nfev']}, polish iterations {R['fp_it']})")
    stat = np.sqrt(np.diag(R['Cstat']))
    if reps:
        Treps = [V.templates(D, L, model, rep=r, rawbg=rawbg)[0] for r in range(V.NREPS)]
        TH = V.replica_fits(model, L, Treps, Vt, R['th'])
        good = np.all(np.isfinite(TH), 1)
        Cmc = np.cov(TH[good].T, ddof=1)
        shift = TH[good].mean(0) - R['th']
        gls, glsd = V.hybrid_gls(model, L, T, Vt, R, Treps if not rawbg else [])
        Ctot = R['Cstat'] + Cmc
        et = np.sqrt(np.diag(Ctot))
        cor = Ctot / np.outer(et, et)
        out.update(stat=stat.tolist(), mc=np.sqrt(np.diag(Cmc)).tolist(), tot=et.tolist(), shift=shift.tolist(),
                   cor=cor.tolist(), gls=gls, glsd=glsd, nrep_ok=int(good.sum()), Cstat=R['Cstat'].tolist(), Cmc=Cmc.tolist(),
                   eN=np.sqrt(np.diag(R['Cn'])).tolist(), eCB=np.sqrt(np.diag(R['Ccb'])).tolist())
        P(f"[{tag}] ok {R['ok']} Pearson {R['pearson']:.1f} hybrid GLS {gls:.1f} data-only GLS {glsd:.1f} ndf {L.n - model.np + len(model.cons)} reps ok {good.sum()}")
        for i, nme in enumerate(model.names):
            P(f"   {nme:10s} {R['th'][i]:9.4f} stat {stat[i]:.4f} (N {np.sqrt(R['Cn'][i, i]):.4f} CB {np.sqrt(R['Ccb'][i, i]):.4f}) MC {np.sqrt(Cmc[i, i]):.4f} tot {et[i]:.4f} shift {shift[i]:+.4f}")
    else:
        P(f"[{tag}] ok {R['ok']} Pearson {R['pearson']:.1f}")
        for i, nme in enumerate(model.names):
            P(f"   {nme:10s} {R['th'][i]:9.4f} stat {stat[i]:.4f}")
    res[tag] = out
    return R, T, Vt


# 1. Part III F-main anchor
P3 = V.P3Model()
run("P3", P3)

# 2. SMASH as is (no fit): Pearson and hybrid GLS
M8 = V.FlexModel(["C2", "C3", "C4", "C5", "Geta", "GD_om", "GD_nom", "GC"])
M7 = V.FlexModel(["C2", "C3", "C4", "C5", "Geta", "GD", "GC"])
for tag, M in (("ASIS", M7), ("ASIS9", M8)):
    T, Vt = V.templates(D, L, M)
    th = M.start.copy()
    mu, J, vbb, k, bgp = M.eval(th, L, T, Vt)
    Vd = np.maximum(mu, 1.0) + L.vcb; Vv = Vd + vbb
    pear = np.sum((L.y - mu) ** 2 / Vv)
    R = dict(th=th, mu=mu, V=Vv, Vd=Vd, y=L.y, SC=L.ccb, eg=L.eg)
    Treps = [V.templates(D, L, M, rep=r)[0] for r in range(V.NREPS)]
    gls, glsd = V.hybrid_gls(M, L, T, Vt, R, Treps)
    P(f"[{tag}] as is: Pearson {pear:.1f} hybrid GLS {gls:.1f} data-only {glsd:.1f}; SMASH sum {(mu - L.cb).sum():.0f}")
    res[tag] = dict(pearson=pear, gls=gls, glsd=glsd)
    # windows (bins whose centre is inside)
    if tag == "ASIS9":
        xc = 0.5 * (L.lo + L.hi)
        for wn, (a, b) in dict(fit=(0.40, 0.70), B=(0.455, 0.5225), P=(0.5336, 0.5542), CRU=(0.565, 0.70)).items():
            s = (xc > a) & (xc < b)
            P(f"   window {wn}: N {L.y[s].sum():.0f} CB {L.cb[s].sum():.0f} SMASH {(mu - L.cb)[s].sum():.0f} (N-CB)/SMASH {(L.y - L.cb)[s].sum() / (mu - L.cb)[s].sum():.4f}")
        for a, c in enumerate(M.comps):
            P(f"   SMASH {c}: {T[:, a].sum():.0f} rows")

# 3. ladder members needed: F7, F8 (nominal), F7bj (alternative), F8bj, F9
run("F7", M7)
R8, T8, Vt8 = run("F8", M8)
run("F7bj", V.FlexModel(["C2", "C3", "C4", "C5", "Geta", "GD", "GC"], tilt=2))
run("F8b", V.FlexModel(["C2", "C3", "C4", "C5", "Geta", "GD_om", "GD_nom", "GC"], tilt=1))
run("F9", V.FlexModel(["C2", "C3", "C4", "C5", "Geta", "GD_om", "GD_nom", "GC_om", "GC_nom"]))
info = {}
V.templates(D, L, M8, info=info)
for kx, vx in sorted(info.items()):
    P("   smoothing", kx, "order %d neff %.1f scaled %s" % vx)
res["smooth_info"] = {str(k): v for k, v in info.items()}
json.dump(res, open(OUT + "/vflex_fit.json", "w"), indent=1)
P(f"done {time.time() - t0:.0f} s")
