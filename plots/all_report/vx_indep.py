#!/usr/bin/env python3
# ============================================================================
# vx_indep.py -- INDEPENDENT VERIFICATION of the X17 limit (prefix vx, 2026-11-14).
# A deliberately simple two-step method that shares NO code with x17c_model.h / x17c_lib.h / x17d_*.C:
#   step 1: in every M(ee) row j (2 MeV/c2, 2-62 MeV/c2) the eta-peak yield Y_j is extracted from the M4 spectrum
#           (X2D range 0.5225-0.6100 GeV/c2, 35 bins) by a Poisson fit
#              mu_i = CB_i (like-sign data, fixed) + E_i (C2+C3+g G_eta at the frozen prior centre a 0.3823, g 1.717, fixed)
#                     + exp(c0 + c1 t + c2 t^2) (smooth other background, free per row) + Y_j gP_i (peak shape = C5 M4
#                     shape of the M(ee) class, from simulation);
#           errors from the inverse Fisher matrix, inflated by sqrt(qScale = 1.0299) (event-group variance factor).
#   step 2: the 30 yields are fitted with Y_j = N5 [D_j + r R_X f_X X_j (1 + sX)^theta] + N4 C_j (D, C, X = simulated
#           M(ee) fractions; f_X = response of the step-1 extractor to the X17 M4 shape; N5, N4 free; theta ~ N(0,1)).
#   limits (90 %): (a) Bayesian, flat prior in r >= 0, N5 and N4 marginalised analytically, theta marginalised;
#                  (b) asymptotic CLs with q~ (profile over N5, N4, theta), sigma_A from the Asimov set of the r = 0 fit;
#                  (c) Gaussian CLs from r_hat +- sigma (Hessian at r_hat).
#   BR(eta -> pi+pi-X17) BR(X17 -> e+e-) = r x 2.68e-4 (normalisation (1)); normalisation (2): X17 rows = r a_eta N5ref R_X.
# Inputs: vx_inputs.npz (vx_dump.py, read-only dump of hists_v1.root, the sealed file and x17_lines.root).
# Output: /hdd3/przygoda/HADES/FAT_out/all_report/vx/vx_indep.txt (+ stdout).
# ============================================================================
import sys
import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm as Nrm

IN = "/hdd3/przygoda/HADES/FAT_out/all_report/vx/vx_inputs.npz"
OUT = "/hdd3/przygoda/HADES/FAT_out/all_report/vx/vx_indep.txt"
d = np.load(IN, allow_pickle=True)
K = 6079.2 * 63.0 / 64.0          # nb -> raw rows
BRDIR = 2.68e-4
A0, G0 = 0.382284, 1.71677         # frozen prior centres (x17c_config.txt)
QS = 1.02988                       # frozen variance factor
SX1, SX2 = 0.246705, 0.435652      # frozen X17 yield nuisances (1), (2)
xe, ye = d["xe"], d["ye"]
nx, ny = len(xe) - 1, len(ye) - 1
xc = 0.5 * (xe[1:] + xe[:-1]); yc = 0.5 * (ye[1:] + ye[:-1])
t = (xc - 0.55) / 0.1
lines = []
def P(s=""):
    print(s); lines.append(s)

# ---------------------------------------------------------------- data
N = d["N_sb"].copy()
masked = np.zeros(ny, bool); masked[5:13] = True
assert np.all(N[:, masked] == -1), "mask not intact"
Ns = d["N_sealed"]
P("vx: sealed side-band cells equal hists_v1.root: max |diff| = %g" % np.max(np.abs(Ns[:, ~masked] - N[:, ~masked])))
N[:, masked] = Ns[:, masked]
P("vx: data rows: all %.0f, side bands %.0f, unblinded rows 12-28 MeV %.0f" % (N.sum(), N[:, ~masked].sum(), N[:, masked].sum()))
CB = d["CB_X2D"]

# ---------------------------------------------------------------- simulation (with flows: index [i+1, j+1])
def inner(a): return a[1:nx + 1, 1:ny + 1]
C2, C3, C4, C5, C6, Ge = (d["sim_" + c] for c in ["C2", "C3", "C4", "C5", "C6", "Geta"])
def tot_mee(a): return a[1:nx + 1, :].sum()     # X2D M4 range, all M(ee) incl. flows
t5, t4, t6 = tot_mee(C5), tot_mee(C4), tot_mee(C6)
RX = t6 * BRDIR / t5
N5ref = t5 * K
P("vx: normalisation from the histograms: sigma_vis(C5) %.4f nb -> N5ref %.1f rows (frozen 7458.8); C6 per unit BR %.1f nb"
  " -> s_X2D %.4e rows (frozen 3.745e7); R_X = eps_X/eps_5 = %.4f (frozen 1.3455)" % (t5, N5ref, t6, t6 * K, RX))
D = inner(C5).sum(0) / t5; Cc = inner(C4).sum(0) / t4; X = inner(C6).sum(0) / t6
P("vx: fractions inside 2-62 MeV: D %.4f, C %.4f, X %.4f" % (D.sum(), Cc.sum(), X.sum()))
E = (inner(C2) + inner(C3) + G0 * inner(Ge)) * K * A0
cls = (yc >= 0.010).astype(int)
gP = np.zeros((2, nx))
for c in (0, 1):
    gP[c] = inner(C5)[:, cls == c].sum(1); gP[c] /= gP[c].sum()
gX = C6[1:nx + 1, :].sum(1); gX /= gX.sum()
gC4 = inner(C4).sum(1); gC4 /= gC4.sum()
def mom(g): m = (g * xc).sum(); return m, np.sqrt((g * (xc - m) ** 2).sum())
P("vx: M4 shapes (mean, rms in MeV/c2): direct E1 %.2f/%.2f, direct E2 %.2f/%.2f, conversion %.2f/%.2f, X17 %.2f/%.2f" %
  tuple(1e3 * v for v in (*mom(gP[0]), *mom(gP[1]), *mom(gC4), *mom(gX))))
# X17c counting window by overlap fractions
def win(a2d):
    fx = np.clip((np.minimum(xe[1:], 0.5566) - np.maximum(xe[:-1], 0.5312)) / np.diff(xe), 0, 1)
    fy = np.clip((np.minimum(ye[1:], 0.0245) - np.maximum(ye[:-1], 0.0155)) / np.diff(ye), 0, 1)
    return (a2d * np.outer(fx, fy)).sum()
P("vx: X17c window (overlap fractions): C6 per unit BR %.1f nb -> s1 %.4e rows (frozen 2.0954e7 = 3501.6 nb x 5984.2); data N %.1f"
  " (exact N_win %.0f), CB %.1f" % (win(inner(C6)), win(inner(C6)) * K, win(N), d["Nwin"][4], win(CB)))

# ---------------------------------------------------------------- step 1: per-row peak yields
def rowfit(n, cb, e, g, order=2, start=None):
    npar = order + 2
    T = np.vstack([t ** k for k in range(order + 1)])      # (order+1, nx)
    def mu(p):
        return cb + e + np.exp(p[1:] @ T) + p[0] * g
    def nll(p):
        m = mu(p)
        if np.any(m <= 0): return 1e30
        return np.sum(m - n * np.log(m))
    def grad(p):
        m = mu(p); w = 1 - n / m; b = np.exp(p[1:] @ T)
        return np.concatenate([[np.sum(w * g)], T @ (w * b)])
    if start is None:
        sb = np.maximum(n - cb - e, 1.0)
        c = np.polyfit(t, np.log(sb), order)[::-1]
        start = np.concatenate([[0.05 * n.sum()], c])
    best = None
    for meth in ("BFGS", "Nelder-Mead"):
        r = minimize(nll, start, jac=grad if meth == "BFGS" else None, method=meth,
                     options={"maxiter": 20000, "gtol": 1e-8} if meth == "BFGS" else {"maxiter": 40000, "xatol": 1e-8, "fatol": 1e-10})
        if best is None or r.fun < best.fun: best = r
        start = best.x
    p = best.x
    m = mu(p); b = np.exp(p[1:] @ T)
    J = np.vstack([g, T * b])                               # (npar, nx)
    F = (J / m) @ J.T
    cov = np.linalg.inv(F)
    return p[0], np.sqrt(cov[0, 0]), p, m

def yields(order=2, addX=None, Ninput=None):
    Y = np.zeros(ny); sY = np.zeros(ny); chi = np.zeros(ny)
    NN = N if Ninput is None else Ninput
    for j in range(ny):
        n = NN[:, j].astype(float)
        y, s, p, m = rowfit(n, CB[:, j], E[:, j], gP[cls[j]], order)
        Y[j], sY[j] = y, s * np.sqrt(QS)
        chi[j] = np.sum((n - m) ** 2 / m)
    return Y, sY, chi

# response of the extractor to the X17 M4 shape: add 1000 rows of gX to an Asimov row and refit
def fX_response(order=2):
    out = []
    for j in (7, 8, 9):
        y0, s0, p, m = rowfit(N[:, j].astype(float), CB[:, j], E[:, j], gP[cls[j]], order)
        y1, _, _, _ = rowfit(m + 1000 * gX, CB[:, j], E[:, j], gP[cls[j]], order, start=p)
        y2, _, _, _ = rowfit(m, CB[:, j], E[:, j], gP[cls[j]], order, start=p)
        out.append((y1 - y2) / 1000.)
    return float(np.mean(out))

# ---------------------------------------------------------------- step 2: M(ee) fit and limits
class MeeFit:
    def __init__(self, Y, sY, Xline, fX, RXv, norm=1, sx=SX1, rows=None):
        sel = np.ones(ny, bool) if rows is None else rows
        self.Y, self.w = Y[sel], 1.0 / sY[sel] ** 2
        self.D, self.C, self.X = D[sel], Cc[sel], Xline[sel] * fX * RXv
        self.norm, self.sx = norm, sx
    def design(self, r, th):
        xf = (1 + self.sx) ** th
        if self.norm == 1:
            A = np.vstack([self.D + r * self.X * xf, self.C]).T; off = 0.0
        else:   # X17 rows independent of N5
            A = np.vstack([self.D, self.C]).T; off = r * A0 * N5ref * self.X * xf
        return A, off
    def chi2_lin(self, r, th):
        A, off = self.design(r, th)
        W = self.w
        M = (A.T * W) @ A
        b = (A.T * W) @ (self.Y - off)
        p = np.linalg.solve(M, b)
        res = self.Y - off - A @ p
        return np.sum(W * res ** 2), np.linalg.slogdet(M)[1], p
    def chi2_prof(self, r, nuis=True):
        if not nuis: return self.chi2_lin(r, 0.0)[0]
        f = lambda th: self.chi2_lin(r, th[0])[0] / QS + th[0] ** 2
        res = minimize(f, [0.0], method="Nelder-Mead", options={"xatol": 1e-6, "fatol": 1e-9})
        return res.fun
    def posterior(self, rgrid, nuis=True):
        ths = np.linspace(-5, 5, 81) if nuis else np.array([0.0])
        wt = np.exp(-0.5 * ths ** 2) if nuis else np.array([1.0])
        post = np.zeros_like(rgrid)
        for k, r in enumerate(rgrid):
            s = 0.0
            for th, w in zip(ths, wt):
                c2, ld, _ = self.chi2_lin(r, th)
                s += w * np.exp(-0.5 * c2 / QS - 0.5 * ld)
            post[k] = s
        return post

def bayes_ul(F, nuis=True, rmax=0.6, n=601):
    rg = np.linspace(0, rmax, n)
    # subtract a reference to avoid underflow
    post = F.posterior(rg, nuis)
    cdf = np.cumsum(0.5 * (post[1:] + post[:-1]) * np.diff(rg)); cdf = np.concatenate([[0], cdf]) / cdf[-1]
    return np.interp(0.90, cdf, rg)

def rhat_and_sigma(F):
    rg = np.linspace(-0.4, 0.6, 2001)
    c = np.array([F.chi2_lin(r, 0.0)[0] / QS for r in rg])
    k = np.argmin(c); rh = rg[k]
    # sigma from Delta chi2 = 1 (stat only, theta = 0)
    lo = np.interp(c[k] + 1, c[:k + 1][::-1], rg[:k + 1][::-1]); hi = np.interp(c[k] + 1, c[k:], rg[k:])
    return rh, 0.5 * (hi - lo), c[k]

def cls_asym_ul(F, nuis=True):
    rg = np.linspace(-0.4, 0.8, 1201)
    c = np.array([F.chi2_prof(r, nuis) for r in rg])
    k = np.argmin(c); rh = rg[k]
    cref = c[k] if rh >= 0 else np.interp(0.0, rg, c)
    # Asimov: Y_A = model at the r = 0 fit (theta = 0)
    _, _, p0 = F.chi2_lin(0.0, 0.0)
    A0m, off0 = F.design(0.0, 0.0)
    FA = MeeFit.__new__(MeeFit); FA.__dict__ = dict(F.__dict__); FA.Y = A0m @ p0 + off0
    cA0 = FA.chi2_prof(0.0, nuis)
    def cls(r):
        q = max(0.0, (np.interp(r, rg, c) - cref)) if r > rh else 0.0
        qA = max(1e-12, FA.chi2_prof(r, nuis) - cA0)
        sq, sA = np.sqrt(q), np.sqrt(qA)
        if q <= qA: ps, pb = 1 - Nrm.cdf(sq), Nrm.cdf(sA - sq)
        else: ps, pb = 1 - Nrm.cdf((q + qA) / (2 * sA)), 1 - Nrm.cdf((q - qA) / (2 * sA))
        return ps / pb
    lo, hi = 1e-4, 0.8
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if cls(mid) > 0.10: lo = mid
        else: hi = mid
    ul = 0.5 * (lo + hi)
    # expected median (Asimov): sigma_A(mu) z with z = Phi^-1(0.95)
    def qAf(r): return max(1e-12, FA.chi2_prof(r, nuis) - cA0)
    lo, hi = 1e-4, 0.8
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if np.sqrt(qAf(mid)) < Nrm.ppf(0.95): lo = mid
        else: hi = mid
    return ul, 0.5 * (lo + hi), rh

def gauss_cls(rh, s):
    # 90 % CLs for a Gaussian estimate with the boundary r >= 0 (= Bayesian flat prior)
    f = lambda mu: (1 - Nrm.cdf((mu - rh) / s)) / (1 - Nrm.cdf(-rh / s)) - 0.10
    lo, hi = 0.0, rh + 10 * s
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0: lo = mid
        else: hi = mid
    return 0.5 * (lo + hi)

def line_for(m):
    # C6 line (17.0 MeV) shifted to mass m (mean shift only, width scaled like the frozen pol2 rms(m)/rms(17)); fractions of
    # the C6 total in the X2D M(ee) rows. At 17.0 the histogram fraction X is used directly.
    if abs(m - 17.0) < 1e-9: return X.copy()
    ed, cont = d["c6_edges"], d["c6_cont"]
    T = cont.sum()
    fg = d["fgrid"]
    mean6 = d["c6mom"][0]
    dmean = (m + np.interp(m, fg, d["fMean"])) - (17.0 + np.interp(17.0, fg, d["fMean"]))
    wr = np.interp(m, fg, d["fRms"]) / np.interp(17.0, fg, d["fRms"])
    lo = 1e-3 * (mean6 + dmean + wr * (ed[:-1] - mean6)); hi = 1e-3 * (mean6 + dmean + wr * (ed[1:] - mean6))
    out = np.zeros(ny)
    for c, a, b in zip(cont[1:-1], lo, hi):
        if c == 0: continue
        ov = np.clip(np.minimum(b, ye[1:]) - np.maximum(a, ye[:-1]), 0, None) / (b - a)
        out += c * ov
    return out / T

def rx_for(m):
    fg = d["fgrid"]
    return RX * np.interp(m, fg, d["fEps_nom"]) / np.interp(17.0, fg, d["fEps_nom"])

# ================================================================ run
P("")
P("vx: STEP 1 (per-row peak yields, background exp(pol2) per row)")
Y, sY, chi = yields(2)
fX = fX_response(2)
P("vx: extractor response to the X17 M4 shape f_X = %.4f (1 = same as the direct-decay shape)" % fX)
P("#  Mee[MeV]   Y_j     sigma   row chi2/31")
for j in range(ny):
    P("  %4.0f-%4.0f %8.1f %7.1f   %5.1f %s" % (1e3 * ye[j], 1e3 * ye[j + 1], Y[j], sY[j], chi[j], "UNBL" if masked[j] else ""))
P("vx: sum of Y_j over 2-62 MeV = %.0f +- %.0f rows" % (Y.sum(), np.sqrt((sY ** 2).sum())))
res = {}
for order in (1, 2, 3):
    if order != 2:
        Yo, sYo, _ = yields(order); fXo = fX_response(order)
    else:
        Yo, sYo, fXo = Y, sY, fX
    for m in (16.7, 16.85, 17.0):
        Xl = line_for(m); rx = rx_for(m)
        F1 = MeeFit(Yo, sYo, Xl, fXo, rx, norm=1)
        F2 = MeeFit(Yo, sYo, Xl, fXo, rx, norm=2, sx=SX2)
        rh, s, cmin = rhat_and_sigma(F1)
        _, _, p = F1.chi2_lin(max(rh, 0), 0)
        ulB = bayes_ul(F1, True); ulB0 = bayes_ul(F1, False)
        ulC, expC, _ = cls_asym_ul(F1, True)
        ulG = gauss_cls(rh, s)
        ulB2 = bayes_ul(F2, True); ulC2, expC2, _ = cls_asym_ul(F2, True)
        res[(order, m)] = (ulB, ulC, expC, ulG, ulB0, rh, s, ulB2, ulC2, p)
        P("vx: order %d m %.2f: N5 %.0f N4 %.0f chi2_min/ndf %.1f/%d | r_hat %+.4f +- %.4f -> BR %+.2e +- %.2e | 90%% UL (1): Bayes %.3e (no yield nuis. %.3e),"
          " CLs-asym %.3e (Asimov median %.3e), Gauss-CLs stat only %.3e | (2): Bayes %.3e CLs %.3e (med %.3e)" %
          (order, m, p[0], p[1], cmin, ny - 3, rh, s, rh * BRDIR, s * BRDIR, ulB * BRDIR, ulB0 * BRDIR, ulC * BRDIR, expC * BRDIR,
           ulG * BRDIR, ulB2 * BRDIR, ulC2 * BRDIR, expC2 * BRDIR))
# a pure counting variant: sum of the peak yields in the rows 14-24 MeV/c2 against N5 D + N4 C fitted on the other rows
P("")
P("vx: COUNTING variant: peak yields of the rows 14-24 MeV/c2 against the prediction N5 D + N4 C fitted on the rows outside 12-28 MeV/c2")
sig = (ye[:-1] >= 0.014 - 1e-9) & (ye[1:] <= 0.024 + 1e-9)
out = ~masked
Fo = MeeFit(Y, sY, X, fX, RX, norm=1, rows=out)
_, _, po = Fo.chi2_lin(0.0, 0.0)
A = np.vstack([D, Cc]).T
Mo = (np.vstack([D[out], Cc[out]]) * Fo.w) @ np.vstack([D[out], Cc[out]]).T
covp = np.linalg.inv(Mo)
pred = float((A[sig] @ po).sum())
vpred = A[sig].sum(0) @ covp @ A[sig].sum(0)
obs = Y[sig].sum(); vobs = (sY[sig] ** 2).sum()
xr = (X[sig] * fX).sum()     # fraction of the X17 rows in the counting rows
ex = obs - pred; eex = np.sqrt(vobs + vpred)
P("vx: rows 14-24: observed peak yield %.0f +- %.0f, predicted %.0f +- %.0f (N5 %.0f, N4 %.0f from the outer rows), excess %+.0f +- %.0f rows;"
  " X17 fraction in these rows x f_X = %.3f" % (obs, np.sqrt(vobs), pred, np.sqrt(vpred), po[0], po[1], ex, eex, xr))
nxul = gauss_cls(ex, eex) / xr
P("vx: COUNTING 90%% CLs (Gaussian, boundary, statistics only) on the X17 rows in X2D: %.0f rows -> BR (1) = N_X/(N5 R_X) x 2.68e-4 = %.2e"
  " (N5 %.0f); BR (2) = N_X/(a_eta s_X2D) = %.2e" % (nxul, nxul / (po[0] * RX) * BRDIR, po[0], nxul / (A0 * t6 * K)))
open(OUT, "w").write("\n".join(lines) + "\n")
print("written", OUT)
