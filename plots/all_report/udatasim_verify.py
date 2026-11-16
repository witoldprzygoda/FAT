#!/usr/bin/env python3
# ============================================================================
# udatasim_verify.py  --  INDEPENDENT check of the udatasim figures (data against SMASH in u)
# ============================================================================
# Written 2026-11-15 (prefix udatasim, verifier). New file. It does NOT import edgeu_core, peaktail_core or any udatasim
# script and does not read any udatasim output: everything is rebuilt from the raw edgeu column files, the D38 weight file
# (c2recon_weights.npz) and stageA.txt, with its own edge function, binning, CB rule, covariances and GLS.
# Output (new files): /hdd3/przygoda/HADES/FAT_out/all_report/udatasim/udatasim_verify.json and udatasim_verify.log
# Run: python3 -I plots/all_report/udatasim_verify.py
# ============================================================================
import json
import math
import time

import numpy as np
import scipy.sparse as sps
from scipy.special import erfc

T0 = time.time()
IN = "/hdd3/przygoda/HADES/FAT_out/all_report/edgeu/"
OUTD = "/hdd3/przygoda/HADES/FAT_out/all_report/udatasim/"
LOG = open(OUTD + "udatasim_verify.log", "w")


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True); LOG.write(s + "\n"); LOG.flush()


RES = {}
# ---------------------------------------------------------------- constants (defs.h)
MET, MPI0 = 0.547862, 0.1349768
KEDGE, KPP = 0.99274, 0.9970
L, LSIM, PT3 = 6079.2, 221.05, 64.0 / 63.0
S4, SPP = 0.99615, 0.9952
P("# udatasim_verify.py (independent check), sim -> raw counts factor L/Lsim/PT3 = %.4f (dump uses 27.0718)" % (L / LSIM / PT3))


def edge(m, q):
    mp = np.asarray(m, dtype=np.float64) / KPP
    q = np.asarray(q, dtype=np.float64)
    E0 = (MET ** 2 + MPI0 ** 2 - mp ** 2) / (2 * MET)
    with np.errstate(invalid="ignore"):
        p0 = np.sqrt(E0 ** 2 - MPI0 ** 2)
        a = MET ** 2 - MET * (E0 - p0) * (1 - q ** 2 / MPI0 ** 2)
        return np.where(a > 0, KEDGE * np.sqrt(a), np.nan)


chk = [(0.30, 0, 0.51859), (0.39, 0, 0.50077), (0.39, 0.060, 0.50958), (0.4116, 0, 0.47390), (0.30, 0.100, 0.53262)]
P("edge check: " + ", ".join("E(%.4f,%.3f) = %.5f (defs.h %.5f)" % (m, q, float(edge(m, q)), r) for m, q, r in chk))
P("            E(0.4120,0) = %s (defs.h: -1, i.e. undefined)" % float(edge(0.412, 0)))


def col(n, dt):
    return np.fromfile(IN + n, dtype=dt)


lay = {l.split()[0]: int(l.split()[1]) for l in open(IN + "dump_layout.txt") if not l.startswith("#")}
# ---------------------------------------------------------------- data
dk = col("data_k.i8", np.int8); dgrp = col("data_grp.i32", np.int32); dm4 = col("data_m4.f32", np.float32)
dpp = col("data_mpipi.f32", np.float32); dee = col("data_mee.f32", np.float32); dbits = col("data_bits.u8", np.uint8)
assert dk.size == lay["data"]
du = dm4.astype(np.float64) - edge(dpp, dee)
du[~(dee < MPI0)] = np.nan
dS0 = (dbits & 1) != 0
# ---------------------------------------------------------------- sim
sn = lay["sim"]
st = col("sim_tree.i8", np.int8); sgrp = col("sim_grp.i32", np.int32); sg = col("sim_g.i8", np.int8); sbits = col("sim_bits.u8", np.uint8)
sm4 = col("sim_m4.f32", np.float32); sm4t = col("sim_m4true.f32", np.float32); spp = col("sim_mpipi.f32", np.float32)
see = col("sim_mee.f32", np.float32); sw = col("sim_w.f64", np.float64); f108 = col("sim_f108.f32", np.float32)
assert st.size == sn and sw.size == sn
selS = ((sbits & 1) != 0) & ((sbits & 8) != 0)
w0 = np.where(selS, sw * f108, 0.0)
GT = ["C2", "C3", "C4", "C5", "Geta", "GDom", "GDnom", "GCom", "GCnom"]
ref_dump = [38688.6, 29075.4, 6896.2, 7671.0, 74670.7, 90662.6, 463029.9, 53572.5, 319518.4]
P("sim S0 mode A V3 sum per group (dump log in brackets): " + ", ".join("%s %.1f (%.1f)" % (GT[i], w0[sg == i].sum(), ref_dump[i]) for i in range(9)))
m4p = S4 * sm4.astype(np.float64)            # r = 1: M4' = 0.99615 [m4_t + (m4 - m4_t)] = 0.99615 m4
mppp = SPP * spp.astype(np.float64)
su = m4p - edge(mppp, see)
su[~(see < MPI0)] = np.nan
mppp32 = mppp.astype(np.float32)

# ---------------------------------------------------------------- D38 weights: own join by (tree, grp, mee, m4)
W = np.load("/hdd3/przygoda/HADES/FAT_out/all_report/c2recon/c2recon_weights.npz")
c2 = np.where(sg == 0)[0]
kw = {}
for i in range(len(W["tree"])):
    kk = (int(W["tree"][i]), int(W["grp"][i]), float(np.float32(W["mee"][i])), float(np.float32(W["m4"][i])))
    assert kk not in kw
    kw[kk] = i
q = np.ones(sn)
nmiss = 0
for i in c2:
    j = kw.get((int(st[i]), int(sgrp[i]), float(see[i]), float(sm4[i])), -1)
    if j < 0:
        nmiss += 1
    else:
        q[i] = W["q_A"][j]
c2s0 = c2[w0[c2] > 0]
P("D38 join: %d C2 counts, %d not found; npz W108 vs dump w x f108 max |diff| %.2e; S0-mode-A weighted mean q_A %.4f (plain mean %.4f, n %d)" %
  (len(c2), nmiss, float(np.max(np.abs(W["W108"] - (sw * f108)[c2]))) if len(W["W108"]) == len(c2) else -1,
   float((w0[c2] * q[c2]).sum() / w0[c2].sum()), float(q[c2s0].mean()), len(c2s0)))
RES["d38"] = dict(n=len(c2), miss=nmiss, wmean=float((w0[c2] * q[c2]).sum() / w0[c2].sum()))
wa = w0 * q

# ---------------------------------------------------------------- grid, M(pipi) regions, M(ee) classes
LO, BW, NB = -0.200, 0.0025, 140
XC = LO + BW * (np.arange(NB) + 0.5)
def ubin(u):
    b = np.full(u.shape, -1, dtype=np.int64)
    ok = np.isfinite(u)
    bb = np.floor((u[ok] - LO) / BW).astype(np.int64)
    bb[(bb < 0) | (bb >= NB)] = -1
    b[ok] = bb
    return b
dub = ubin(du); sub = ubin(su)
def inpp(m32, lo, hi):
    return (m32 >= np.float32(lo)) & (m32 < np.float32(hi))
dE1 = dee < 0.010; sE1 = see < 0.010
CLS = {"all": (np.ones(dk.size, bool), np.ones(sn, bool)), "E1": (dE1, sE1), "E2": (~dE1, ~sE1)}
NC, RB = NB // 4, 4
CC = LO + 0.010 * (np.arange(NC) + 0.5)
# interpolation matrix (cell -> fine bin), own construction with np.interp on unit vectors
AI = np.zeros((NB, NC))
for j in range(NC):
    e = np.zeros(NC); e[j] = 1.0
    AI[:, j] = np.interp(XC, CC, e) / RB


def hist(b, m, wts=None):
    ok = m & (b >= 0)
    return np.bincount(b[ok], weights=None if wts is None else wts[ok], minlength=NB).astype(np.float64)


def data_spec(slices, cls):
    cm = CLS[cls][0]
    N = np.zeros(NB); CB = np.zeros(NB); VCB = np.zeros((NB, NB))
    shm = (XC >= -0.100) & (XC < 0.100)
    for lo, hi in slices:
        base = dS0 & inpp(dpp, lo, hi)
        N += hist(dub, base & (dk == 0) & cm)
        pa = hist(dub, base & (dk == 1)); ma = hist(dub, base & (dk == 2))
        pc = hist(dub, base & (dk == 1) & cm); mc = hist(dub, base & (dk == 2) & cm)
        pw = pa.reshape(NC, RB).sum(1); mw = ma.reshape(NC, RB).sum(1)
        cbw = 2 * np.sqrt(pw * mw); vw = pw + mw
        sh = (pc + mc)[shm].sum() / (pa + ma)[shm].sum()
        CB += sh * (AI @ cbw); VCB += sh * sh * (AI * vw) @ AI.T
    # event-group covariance of N
    m = dS0 & (dk == 0) & cm & (dub >= 0) & np.zeros(dk.size, bool)
    for lo, hi in slices:
        m |= dS0 & (dk == 0) & cm & (dub >= 0) & inpp(dpp, lo, hi)
    uq, inv = np.unique(dgrp[m], return_inverse=True)
    Mg = sps.csr_matrix((np.ones(int(m.sum())), (inv, dub[m])), shape=(uq.size, NB))
    VN = (Mg.T @ Mg).toarray()
    return N, CB, VCB, VN


def sim_mask(slices, cls):
    m = np.zeros(sn, bool)
    for lo, hi in slices:
        m |= inpp(mppp32, lo, hi)
    return m & CLS[cls][1] & (sub >= 0)


def sim_cov(wts, m):
    ok = m & (wts != 0)
    key = st[ok].astype(np.int64) * (1 << 32) + sgrp[ok].astype(np.int64)
    uq, inv = np.unique(key, return_inverse=True)
    Mg = sps.csr_matrix((wts[ok], (inv, sub[ok])), shape=(uq.size, NB))
    return (Mg.T @ Mg).toarray()


def gls(y, F, V):
    Vi = np.linalg.inv(V)
    C = np.linalg.inv(F.T @ Vi @ F)
    Lm = C @ F.T @ Vi
    p = Lm @ y
    r = y - F @ p
    return p, C, Lm, float(r @ Vi @ r)


def sel(lo, hi):
    return np.where((XC > lo) & (XC < hi))[0]


SL1 = [(0.30, 0.34), (0.34, 0.37), (0.37, 0.39), (0.39, 0.40)]
SLA = [(0.27, 0.30)] + SL1 + [(0.40, 0.415)]

# ---------------------------------------------------------------- 1. f_eta(u) with an own step fit (edgeu definition)
N1, CB1, VCB1, VN1 = data_spec(SL1, "all")
sig1 = N1 - CB1
mk1 = sim_mask(SL1, "all")
def simh(wts, groups, m):
    return hist(sub, m & np.isin(sg, groups), wts)
peak = simh(w0, [2, 3], mk1)
# template width: Gaussian core of u_reco - u_true of C2 + C3 near the edge (edgeu rule)
from scipy.optimize import curve_fit
def m4max(m, qq):
    E0 = (MET ** 2 + MPI0 ** 2 - m ** 2) / (2 * MET)
    with np.errstate(invalid="ignore"):
        p0 = np.sqrt(E0 ** 2 - MPI0 ** 2)
        a = MET ** 2 - MET * (E0 - p0) * (1 - qq ** 2 / MPI0 ** 2)
        return np.where(a > 0, np.sqrt(a), np.nan)
ut = KEDGE * (sm4t.astype(np.float64) - m4max(col("sim_mpipitrue.f32", np.float32).astype(np.float64), col("sim_meetrue.f32", np.float32).astype(np.float64)))
mr = (sg <= 1) & np.isfinite(su) & np.isfinite(ut) & (w0 > 0) & (mppp >= 0.30) & (mppp < 0.40) & (ut > -0.015)
hh, ee = np.histogram((su - ut)[mr], bins=np.arange(-0.040, 0.0401, 0.002), weights=w0[mr])
hh2, _ = np.histogram((su - ut)[mr], bins=ee, weights=w0[mr] ** 2)
xx = 0.5 * (ee[1:] + ee[:-1])
pg, _ = curve_fit(lambda x, a, mu, s: a * np.exp(-0.5 * (x - mu) ** 2 / s ** 2), xx, hh, p0=[hh.max(), 0, 0.007], sigma=np.sqrt(np.maximum(hh2, 1)))
SIG = abs(pg[2])
P("template width (own Gaussian core fit): %.4f MeV/c2 (note: 6.7019)" % (SIG * 1e3))
def step_T(s):
    def Gi(x):
        return 0.5 * (x * erfc(x / (math.sqrt(2) * s)) - math.sqrt(2 / math.pi) * s * np.exp(-x * x / (2 * s * s)))
    ed = LO + BW * np.arange(NB + 1)
    return (Gi(ed[1:]) - Gi(ed[:-1])) / BW
AP = 0.418249
def fstep(sigma):
    T = step_T(sigma)
    ii = sel(-0.045, 0.045)
    t = XC[ii] / 0.045
    F = np.vstack([np.ones_like(t), t, t * t, T[ii]]).T
    fx = AP * peak[ii]; y = sig1[ii] - fx
    V0 = VCB1[np.ix_(ii, ii)]
    V = V0 + np.diag(np.maximum(N1[ii], 1.0))
    for it in range(8):
        p, C, Lm, c2v = gls(y, F, V)
        V = V0 + np.diag(np.maximum(F @ p + fx + CB1[ii], 1.0))
    p, C, Lm, c2v = gls(y, F, V)
    hC2 = Lm[-1] @ simh(wa, [0], mk1)[ii]; hC2raw = Lm[-1] @ simh(w0, [0], mk1)[ii]; hC3 = Lm[-1] @ simh(w0, [1], mk1)[ii]
    c2fac = np.where(see < 0.010, 0.879, 0.8825)
    hC2c = Lm[-1] @ simh(w0 * c2fac, [0], mk1)[ii]
    return dict(h=float(p[-1]), eh=float(math.sqrt(C[-1, -1])), chi2=c2v, ndf=len(ii) - 4, hC2=float(hC2), hC2raw=float(hC2raw), hC3=float(hC3),
                hC2c=float(hC2c), f38=float(p[-1] / (hC2 + hC3)), f35=float(p[-1] / (hC2c + hC3)))
FS = fstep(0.0067019); FS2 = fstep(SIG)
P("own step fit (sigma 6.7019): h_data %.3f +- %.3f, chi2 %.2f/%d; C2 raw %.3f, C2 x 0.879/0.8825 %.3f, C2 (D38) %.3f, C3 %.3f; f(D35) %.5f, f(D38) %.5f; C2 eff. factor %.4f" %
  (FS["h"], FS["eh"], FS["chi2"], FS["ndf"], FS["hC2raw"], FS["hC2c"], FS["hC2"], FS["hC3"], FS["f35"], FS["f38"], FS["hC2"] / FS["hC2raw"]))
P("own step fit (own sigma %.4f): h %.3f, f(D38) %.5f" % (SIG * 1e3, FS2["h"], FS2["f38"]))
RES["fstep"] = FS; RES["fstep_ownsigma"] = dict(sigma=SIG, **FS2)
F38 = FS["f38"]

# ---------------------------------------------------------------- 2. model weights
SA = open("/hdd3/przygoda/HADES/FAT_out/all_report/wp3/stageA.txt").read().splitlines()
vals = None; cov = None
for i, l in enumerate(SA):
    t = l.split()
    if t and t[0] == "values":
        vals = np.array([float(v) for v in t[1:]])
    if t and t[0] == "cov_total":
        n = int(t[1]); cov = np.array([[float(v) for v in SA[i + 1 + r].split()] for r in range(n)])
CJ, BT, KB = vals[:5], vals[5], vals[6]
COV = cov[:7, :7]
P("stage A: c %s, b %.6f, kappa_B %.5f, err %s" % (np.round(CJ, 5), BT, KB, np.round(np.sqrt(np.diag(COV)), 4)))
AG = 0.656
jsl = np.searchsorted(np.array([0.34, 0.37, 0.39, 0.43]), mppp, side="right")
xt = (m4p - 0.55) / 0.1
isGD = (sg == 5) | (sg == 6); isGC = (sg == 7) | (sg == 8)
cA = CJ[jsl] * (1 + BT * xt)
scb = np.ones(sn)
scb[sg <= 1] = F38; scb[(sg == 2) | (sg == 3)] = AP; scb[sg == 4] = AG
scb[isGD] = cA[isGD]; scb[isGC] = cA[isGC] * KB
wb = wa * scb
base = np.where(isGC, KB, 1.0)
JW = [np.where((isGD | isGC) & (jsl == j), wa * (1 + BT * xt) * base, 0.0) for j in range(5)]
JW.append(np.where(isGD | isGC, wa * CJ[jsl] * xt * base, 0.0))
JW.append(np.where(isGC, wa * cA, 0.0))

GRPS = {"C2": [0], "C3": [1], "C4": [2], "C5": [3], "Geta": [4], "GD": [5, 6], "GC": [7, 8]}
BANDS = {"below": (-0.045, -0.005), "edge": (-0.005, 0.005), "above": (0.005, 0.045)}
RANGES = {"fit": (-0.045, 0.045), "peak": (0.020, 0.060)}


def analyse(tag, slices, cls):
    N, CB, VCB, VN = data_spec(slices, cls)
    sig = N - CB
    VD = VN + VCB
    m = sim_mask(slices, cls)
    out = dict(N_over_var=float(np.trace(VN) / N.sum()))
    comp = {}
    for mk, wts in (("a", wa), ("b", wb)):
        parts = {g: simh(wts, gs, m) for g, gs in GRPS.items()}
        M = sum(parts.values()); B = parts["GD"] + parts["GC"]
        VM = sim_cov(wts, m)
        if mk == "b":
            Jm = np.vstack([hist(sub, m, jw) for jw in JW]).T
            VA = Jm @ COV @ Jm.T
        else:
            VA = np.zeros((NB, NB))
        o = {}
        for rn, (lo, hi) in RANGES.items():
            ii = sel(lo, hi)
            Vt = (VD + VM)[np.ix_(ii, ii)]; Vi = np.linalg.inv(Vt)
            r = (sig - M)[ii]
            c2s = float(r @ Vi @ r)
            ViA = np.linalg.inv(Vt + VA[np.ix_(ii, ii)])
            c2A = float(r @ ViA @ r)
            y = (sig - M + B)[ii]; bb = B[ii]
            lam = float(bb @ Vi @ y / (bb @ Vi @ bb)); elam = float(1 / math.sqrt(bb @ Vi @ bb))
            rl = y - lam * bb
            o[rn] = dict(n=len(ii), chi2=c2s, chi2A=c2A, lam=lam, elam=elam, chi2lam=float(rl @ Vi @ rl))
        lamf = o["fit"]["lam"]
        Ml = M + (lamf - 1) * B
        ii = sel(0.020, 0.060)
        rr = (sig - Ml)[ii]
        o["peak"]["chi2_lamfit_fixed"] = float(rr @ np.linalg.inv((VD + VM)[np.ix_(ii, ii)]) @ rr)
        # level significance in counts: (d - M) / sqrt(var_d + var_MC + var_A) per band
        o["z_counts"] = {}
        for bn, (lo, hi) in BANDS.items():
            e = np.zeros(NB); e[sel(lo, hi)] = 1
            o["z_counts"][bn] = float((sig - M) @ e / math.sqrt(e @ (VD + VM + VA) @ e))
        ipl = sel(-0.100, 0.100); ifr = sel(-0.045, 0.045)
        o["geta_share"] = dict(fit=float(parts["Geta"][ifr].sum() / M[ifr].sum()), plot=float(parts["Geta"][ipl].sum() / M[ipl].sum()))
        bd = {}
        for bn, (lo, hi) in BANDS.items():
            ii = sel(lo, hi); e = np.zeros(NB); e[ii] = 1
            d = sig @ e; mm = M @ e; ml = Ml @ e
            R = d / mm
            bd[bn] = dict(R=float(R), stat=float(R * math.sqrt((e @ VD @ e) / d ** 2 + (e @ VM @ e) / mm ** 2)), A=float(R * math.sqrt(e @ VA @ e) / mm),
                          Rlam=float(d / ml), data=float(d), CB=float(CB @ e),
                          comp={g: float(parts[g] @ e) for g in parts})
        o["bands"] = bd
        out[mk] = o
        comp[mk] = parts
    RES[tag] = out
    a, b = out["a"], out["b"]
    P("\n== %s (%s)" % (tag, cls))
    for mk in ("a", "b"):
        o = out[mk]
        P("  (%s) chi2 -45..+45 %.1f/%d (+A %.1f), lam %.3f +- %.3f -> %.1f | +20..+60 %.1f/%d (+A %.1f), lam %.3f -> %.1f" %
          (mk, o["fit"]["chi2"], o["fit"]["n"], o["fit"]["chi2A"], o["fit"]["lam"], o["fit"]["elam"], o["fit"]["chi2lam"],
           o["peak"]["chi2"], o["peak"]["n"], o["peak"]["chi2A"], o["peak"]["lam"], o["peak"]["chi2lam"]))
        P("      peak chi2 with the -45..+45 lambda kept fixed: %.2f; z (counts, stat+MC+A): %s; G_eta share fit range %.3f, plotted range %.3f" %
          (o["peak"]["chi2_lamfit_fixed"], ", ".join("%s %.2f" % (bn, o["z_counts"][bn]) for bn in BANDS), o["geta_share"]["fit"], o["geta_share"]["plot"]))
        P("      bands: " + "; ".join("%s %.3f +- %.3f +- %.3f (A) [lam %.3f]" % (bn, o["bands"][bn]["R"], o["bands"][bn]["stat"], o["bands"][bn]["A"], o["bands"][bn]["Rlam"]) for bn in BANDS))
        for bn in BANDS:
            c = o["bands"][bn]["comp"]
            P("      %-5s data N-CB %.0f (CB %.0f): C2+C3 %.0f, C4+C5 %.0f, Geta %.0f, GD %.0f, GC %.0f" %
              (bn, o["bands"][bn]["data"], o["bands"][bn]["CB"], c["C2"] + c["C3"], c["C4"] + c["C5"], c["Geta"], c["GD"], c["GC"]))
    P("  N group var / N = %.4f" % out["N_over_var"])
    return N, CB, sig, comp


mr1 = sim_mask(SL1, "all") & (np.abs(su) < 0.1)
for nm, gs in (("G_D", [5, 6]), ("G_C", [7, 8]), ("G_eta", [4]), ("C2", [0])):
    m_ = mr1 & np.isin(sg, gs) & (w0 > 0)
    P("weight per simulated count in r1, |u| < 100 MeV/c2, %s: n %d, mean %.2f, median %.2f, sum w^2/sum w %.2f expected counts" %
      (nm, int(m_.sum()), float(w0[m_].mean()), float(np.median(w0[m_])), float((w0[m_] ** 2).sum() / w0[m_].sum())))
SPEC = {}
for tag, sl in (("r1", SL1), ("allmpp", SLA)):
    for cls in ("all", "E1", "E2"):
        SPEC[(tag, cls)] = analyse(tag + "_" + cls, sl, cls)
    P("   (%.0f s)" % (time.time() - T0))

# ---------------------------------------------------------------- 3. bin-by-bin arrays for the comparison with the figure data
np.savez_compressed(OUTD + "udatasim_verify_arrays.npz", x=XC,
                    **{"%s_%s_%s" % (t, c, n): v for (t, c), (N, CB, sig, comp) in SPEC.items() for n, v in (("N", N), ("CB", CB), ("sig", sig))},
                    **{"%s_%s_%s_%s" % (t, c, mk, g): comp[mk][g] for (t, c), (N, CB, sig, comp) in SPEC.items() for mk in ("a", "b") for g in GRPS})
json.dump(RES, open(OUTD + "udatasim_verify.json", "w"), indent=1)
P("done (%.0f s)" % (time.time() - T0))
