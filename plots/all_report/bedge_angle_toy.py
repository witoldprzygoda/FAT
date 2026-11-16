#!/usr/bin/env python3
# ============================================================================
# bedge_angle_toy.py  --  package bedge: the HGeant gamma* decay-angle effect of class C2 AT THE KINEMATIC EDGE
# ============================================================================
# Written 2026-11-14 (prefix bedge). A TOY, not a detector simulation.
#
# Question: WP8 found that the simulated C2 efficiency (pi0 Dalitz in the bump) is 11-14 % too high because HGeant
# decays gamma* -> e+e- isotropically (QED: 1 + cos^2 + (1 - beta^2) sin^2). Its factors (toy 0.880, closure 0.891)
# are averages over all M4. The edge estimator uses only the rows just below the kinematic edge
# M4max(M(pipi)), where the missing photon is softest and the gamma* is fastest. Is the factor the same there?
#
# Generator: a COPY of the event generator of wp8_toy_angle.py (WP8; functions load_gen, sample_eta, iso_dir, boost,
# kloe_weight, sample_q, generate; unchanged except that the true M(pi+pi-) and the true M4 are returned):
#   eta (pT, y_cm) from the SMASH std generator histogram (40 files); eta -> pi+pi-pi0 phase space x KLOE-2 weight;
#   pi0 -> gamma gamma*, q from Kroll-Wada x (1 + a q^2/m^2)^2; gamma* -> e+e- isotropic (HGeant) and QED (weight,
#   same events). Acceptance (wp8_toy_angle.py accept): 4 tracks 18-85 deg, pions p > 0.10, leptons p > p_min.
# New here: M4 and M(pipi) are put on the data scale (x 0.99277, x 0.9970), M4 is smeared with a Gaussian
# (5.4 MeV nominal, 6.1 MeV variant), and the frozen edge-step fit (pol2 + h x step at the fixed edge, smearing
# 5.4 MeV, range edge -+ 45 MeV (-+ 35 MeV from 0.39 on), no eta term) is applied to the weighted toy spectra of
# the four M(pipi) slices, once with the isotropic and once with the QED weights. Ratio R = sum h(QED)/sum h(iso).
#
# Output: /hdd3/przygoda/HADES/FAT_out/all_report/bedge/bedge_angle.txt (lines "R <variant> <class> <value> <err>")
#         and docs/latex/tables_all/bedge_angle.txt (readable table).
# Run: python3 -I plots/all_report/bedge_angle_toy.py      (about 2 min, 4 x 3e6 decays)
# ============================================================================
import math
import numpy as np
from scipy.special import erfc

GEN = "/hdd3/przygoda/HADES/FAT_out/all_report/wp8/remote/gen_eta_ptY_std40.txt"
OUT1 = "/hdd3/przygoda/HADES/FAT_out/all_report/bedge/bedge_angle.txt"
OUT2 = "/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all/bedge_angle.txt"
ME, MPI, MPI0, META, MP = 0.000510999, 0.13957039, 0.1349768, 0.547862, 0.938272
A_FF = 0.0332
SCALE, SCALE_PP = 0.5439 / 0.547862, 0.9970          # ScanEngine kScale, kScalePiPi (data scale)
rng = np.random.default_rng(20261115)


# ------------------------------------------------------------------ copy of wp8_toy_angle.py (WP8)
def load_gen():
    h = np.zeros((100, 100))
    with open(GEN) as f:
        for line in f:
            if line.startswith("#"):
                continue
            a = line.split()
            if len(a) == 3:
                h[int(a[0]), int(a[1])] += float(a[2])
    return h


def sample_eta(n, h):
    p = h.ravel() / h.sum()
    idx = rng.choice(len(p), size=n, p=p)
    ip, iy = idx // 100, idx % 100
    pt = (ip + rng.random(n)) * 0.020
    y = -2.5 + (iy + rng.random(n)) * 0.05
    return pt, y


def iso_dir(n):
    c = 2 * rng.random(n) - 1
    s = np.sqrt(1 - c * c)
    ph = 2 * np.pi * rng.random(n)
    return np.stack([s * np.cos(ph), s * np.sin(ph), c], axis=1)


def boost(p4, b):
    b2 = np.sum(b * b, axis=1)
    g = 1.0 / np.sqrt(1.0 - b2)
    bp = np.sum(b * p4[:, 1:], axis=1)
    g2 = np.where(b2 > 0, (g - 1.0) / np.where(b2 > 0, b2, 1.0), 0.0)
    E = g * (p4[:, 0] + bp)
    vec = p4[:, 1:] + (g2 * bp)[:, None] * b + (g * p4[:, 0])[:, None] * b
    return np.column_stack([E, vec])


def kloe_weight(mpp2, s13, s23):
    Q = META - 2 * MPI - MPI0
    E0 = (META ** 2 + MPI0 ** 2 - mpp2) / (2 * META)
    T0 = E0 - MPI0
    Ep = (META ** 2 + MPI ** 2 - s23) / (2 * META)
    Em = (META ** 2 + MPI ** 2 - s13) / (2 * META)
    Tp, Tm = Ep - MPI, Em - MPI
    X = math.sqrt(3) * (Tp - Tm) / Q
    Y = 3 * T0 / Q - 1
    a, b, d, f, g = -1.095, 0.145, 0.081, 0.141, -0.044
    return 1 + a * Y + b * Y * Y + d * X * X + f * Y ** 3 + g * X * X * Y


def sample_q(n):
    out = np.empty(0)
    qmin, qmax = 2 * ME * 1.000001, MPI0 * 0.999999
    while out.size < n:
        m = 4 * (n - out.size) + 1000
        u = rng.random(m)
        q = qmin * (qmax / qmin) ** u
        x = (q / MPI0) ** 2
        beta = np.sqrt(1 - 4 * ME * ME / (q * q))
        f = (1 - x) ** 3 * (1 + 2 * ME * ME / (q * q)) * beta * (1 + A_FF * x) ** 2
        keep = rng.random(m) * 2.0 < f
        out = np.concatenate([out, q[keep]])
    return out[:n]


def generate(n, h):
    pt, ycm = sample_eta(n, h)
    Eb = 4.5 + MP
    pzb = math.sqrt(Eb * Eb - MP * MP)
    yc = 0.5 * math.log((Eb + MP + pzb) / (Eb + MP - pzb))
    y = ycm + yc
    phi = 2 * np.pi * rng.random(n)
    mt = np.sqrt(META ** 2 + pt * pt)
    eta = np.column_stack([mt * np.cosh(y), pt * np.cos(phi), pt * np.sin(phi), mt * np.sinh(y)])
    s12lo, s12hi = (2 * MPI) ** 2, (META - MPI0) ** 2
    s12 = np.empty(0); s13 = np.empty(0)
    while s12.size < n:
        m = 3 * (n - s12.size) + 1000
        a = s12lo + (s12hi - s12lo) * rng.random(m)
        b = (MPI + MPI0) ** 2 + ((META - MPI) ** 2 - (MPI + MPI0) ** 2) * rng.random(m)
        E2s = np.sqrt(a) / 2
        E3s = (META ** 2 - a - MPI0 ** 2) / (2 * np.sqrt(a))
        p2s = np.sqrt(np.maximum(0, E2s ** 2 - MPI ** 2)); p3s = np.sqrt(np.maximum(0, E3s ** 2 - MPI0 ** 2))
        lo = (E2s + E3s) ** 2 - (p2s + p3s) ** 2
        hi = (E2s + E3s) ** 2 - (p2s - p3s) ** 2
        ok = (b >= lo) & (b <= hi)
        s12 = np.concatenate([s12, a[ok]]); s13 = np.concatenate([s13, b[ok]])
    s12, s13 = s12[:n], s13[:n]
    s23 = META ** 2 + 2 * MPI ** 2 + MPI0 ** 2 - s12 - s13
    wk = kloe_weight(s12, s13, s23)
    E0 = (META ** 2 + MPI0 ** 2 - s12) / (2 * META)
    p0 = np.sqrt(np.maximum(0, E0 ** 2 - MPI0 ** 2))
    d0 = iso_dir(n)
    pi0 = np.column_stack([E0, p0[:, None] * d0])
    pp = np.column_stack([META - E0, -p0[:, None] * d0])
    mpp = np.sqrt(s12)
    kp = np.sqrt(np.maximum(0, mpp ** 2 / 4 - MPI ** 2))
    dpi = iso_dir(n)
    pip = np.column_stack([np.sqrt(kp ** 2 + MPI ** 2), kp[:, None] * dpi])
    pim = np.column_stack([np.sqrt(kp ** 2 + MPI ** 2), -kp[:, None] * dpi])
    bpp = pp[:, 1:] / pp[:, :1]
    pip = boost(pip, bpp); pim = boost(pim, bpp)
    q = sample_q(n)
    ps = (MPI0 ** 2 - q * q) / (2 * MPI0)
    dg = iso_dir(n)
    gs = np.column_stack([np.sqrt(ps * ps + q * q), ps[:, None] * dg])
    k = np.sqrt(np.maximum(0, q * q / 4 - ME * ME))
    beta = 2 * k / q
    c = 2 * rng.random(n) - 1
    s = np.sqrt(1 - c * c)
    ph = 2 * np.pi * rng.random(n)
    ez = dg
    tmp = np.where(np.abs(ez[:, :1]) < 0.9, np.array([[1.0, 0, 0]]), np.array([[0, 1.0, 0]]))
    ex = np.cross(tmp, ez); ex /= np.linalg.norm(ex, axis=1)[:, None]
    ey = np.cross(ez, ex)
    dl = (c[:, None] * ez + (s * np.cos(ph))[:, None] * ex + (s * np.sin(ph))[:, None] * ey)
    Ek = np.sqrt(k * k + ME * ME)
    ep = np.column_stack([Ek, k[:, None] * dl]); em = np.column_stack([Ek, -k[:, None] * dl])
    bg = gs[:, 1:] / gs[:, :1]
    ep = boost(ep, bg); em = boost(em, bg)
    wq = 2 * (1 + c * c + (1 - beta ** 2) * (1 - c * c)) / (8.0 / 3 + 4.0 / 3 * (1 - beta ** 2))
    b0 = pi0[:, 1:] / pi0[:, :1]
    ep = boost(ep, b0); em = boost(em, b0)
    be = eta[:, 1:] / eta[:, :1]
    ep, em, pip, pim = boost(ep, be), boost(em, be), boost(pip, be), boost(pim, be)
    P4 = ep + em + pip + pim                      # new here: true M4 = M(pi+ pi- e+ e-)
    m4 = np.sqrt(np.maximum(0, P4[:, 0] ** 2 - np.sum(P4[:, 1:] ** 2, axis=1)))
    return dict(q=q, wk=wk, wq=wq, ep=ep, em=em, pip=pip, pim=pim, mpp=mpp, m4=m4)


def accept(ev, pmin_e, oa_min=0.0):
    ok = np.ones(ev["q"].size, bool)
    for key, pmin in (("ep", pmin_e), ("em", pmin_e), ("pip", 0.10), ("pim", 0.10)):
        v = ev[key][:, 1:]
        p = np.linalg.norm(v, axis=1)
        th = np.degrees(np.arccos(np.clip(v[:, 2] / p, -1, 1)))
        ok &= (th > 18) & (th < 85) & (p > pmin) & (p < 2.0)
    if oa_min > 0:
        a, b = ev["ep"][:, 1:], ev["em"][:, 1:]
        cosoa = np.sum(a * b, axis=1) / (np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1))
        ok &= np.degrees(np.arccos(np.clip(cosoa, -1, 1))) > oa_min
    return ok


# ------------------------------------------------------------------ the frozen edge-step fit (ScanEngine FitEdgeStep, no CB, no eta term)
def edge_meas(mp):
    mpr = mp / SCALE_PP
    E0 = (META ** 2 + MPI0 ** 2 - mpr ** 2) / (2 * META)
    p2 = E0 * E0 - MPI0 ** 2
    if p2 < 0:
        return -1
    a = META ** 2 - META * (E0 - math.sqrt(p2))
    return SCALE * math.sqrt(a)


def step_tmpl(x, lo, hi, s):
    v = np.zeros_like(x)
    n = 21
    for i in range(n):
        mp = lo + (i + 0.5) * (hi - lo) / n
        e = edge_meas(mp)
        v += 0.5 * erfc((x - e) / (math.sqrt(2) * s))
    return v / n


SLICES = [(0.30, 0.34, 0.045), (0.34, 0.37, 0.045), (0.37, 0.39, 0.045), (0.39, 0.40, 0.035)]
BW, MLO, NB = 0.0025, 0.35, 240


def fit_step(y, v, sl, s=0.0054):
    lo_s, hi_s, hr = sl
    edge = edge_meas(0.5 * (lo_s + hi_s))
    lo, hi = edge - hr, edge + hr
    # no eta term in the toy; keep the range (the upper extension of the frozen fit only adds the eta Gaussian)
    lo = MLO + BW * math.floor((lo - MLO) / BW + 0.5); hi = MLO + BW * math.floor((hi - MLO) / BW + 0.5)
    x = MLO + (np.arange(NB) + 0.5) * BW
    m = (x >= lo) & (x < hi)
    xs, ys, vs = x[m], y[m], np.maximum(v[m], 1e-12)
    x0, dx = 0.5 * (lo + hi), 0.5 * (hi - lo)
    t = (xs - x0) / dx
    F = np.column_stack([np.ones_like(t), t, t * t, step_tmpl(xs, lo_s, hi_s, s)])
    W = 1.0 / vs
    A = F.T @ (F * W[:, None])
    cov = np.linalg.inv(A)
    L = np.zeros(NB); L[m] = (cov @ (F * W[:, None]).T)[3]      # h = L . y (fixed weights)
    return float(L @ y), L


def main():
    h = load_gen()
    nchunk, n = 4, 3_000_000
    variants = [("nominal (p_e > 0.06, smear 5.4 MeV)", 0.06, 0.0054), ("p_e > 0.04", 0.04, 0.0054), ("p_e > 0.10", 0.10, 0.0054),
                ("smear 6.1 MeV", 0.06, 0.0061)]
    classes = [("all", 0.0, 1.0), ("E1", 0.0, 0.010), ("E2", 0.010, 1.0)]
    # accumulate weighted histograms per variant, class, slice: iso (w), qed (w wq), and the variances (w^2) for the fit
    H = {}
    tot = {}
    for ic in range(nchunk):
        ev = generate(n, h)
        mpp_d = SCALE_PP * ev["mpp"]
        for vname, pmin, sm in variants:
            acc = accept(ev, pmin)
            m4s = SCALE * ev["m4"] + rng.normal(0, sm, ev["m4"].size)
            for cname, qlo, qhi in classes:
                mc = acc & (ev["q"] >= qlo) & (ev["q"] < qhi)
                wi = ev["wk"]; wqd = ev["wk"] * ev["wq"]
                key = (vname, cname)
                t = tot.setdefault(key, np.zeros(4))
                t += [wi[mc].sum(), wqd[mc].sum(), 0, 0]
                for isl, (lo_s, hi_s, hr) in enumerate(SLICES):
                    ms = mc & (mpp_d >= lo_s) & (mpp_d < hi_s)
                    b = np.floor((m4s[ms] - MLO) / BW).astype(int)
                    ok = (b >= 0) & (b < NB)
                    arr = H.setdefault(key + (isl,), np.zeros((4, NB)))
                    for j, w in enumerate((wi[ms][ok], wqd[ms][ok], wi[ms][ok] ** 2, (wqd[ms][ok] - wi[ms][ok]) ** 2)):
                        arr[j] += np.bincount(b[ok], weights=w, minlength=NB)
        print("chunk %d done" % ic, flush=True)
    lines = []
    def P(s=""):
        print(s); lines.append(s)
    P("# bedge_angle.txt -- plots/all_report/bedge_angle_toy.py (prefix bedge, 2026-11-14). TOY (generator copied from wp8_toy_angle.py).")
    P("# R_edge = sum over the 4 slices of h(QED weights) / h(isotropic), h from the frozen edge-step fit (pol2 + step, 5.4 MeV);")
    P("# R_all = acceptance ratio over all M4 and all M(pipi) (= the WP8 toy quantity); classes in TRUE M(ee); %d decays" % (nchunk * n))
    P("")
    P("%-38s %-4s %8s %8s %8s %8s %8s %8s   %s" % ("variant", "cls", "R_all", "R_edge", "+-", "R_sl1", "R_sl2", "R_sl3", "R_sl4"))
    mach = []
    for vname, pmin, sm in variants:
        for cname, qlo, qhi in classes:
            key = (vname, cname)
            hi_sum = hq_sum = 0.0
            rs = []
            # the same fixed weights (isotropic variances) for both: h = L . y, so h(QED) - h(iso) = L . (y_QED - y_iso)
            # with the per-bin variance sum (w (wq - 1))^2 of the same events -> error of the ratio
            vd = 0.0
            for isl, sl in enumerate(SLICES):
                arr = H[key + (isl,)]
                h_i, L = fit_step(arr[0], arr[2], sl)
                h_q, _ = fit_step(arr[1], arr[2], sl)
                vd += float((L * L) @ arr[3])
                hi_sum += h_i; hq_sum += h_q
                rs.append(h_q / h_i if h_i != 0 else float("nan"))
            R_edge = hq_sum / hi_sum
            eR = math.sqrt(vd) / abs(hi_sum)
            t = tot[key]
            R_all = t[1] / t[0]
            P("%-38s %-4s %8.4f %8.4f %8.4f " % (vname, cname, R_all, R_edge, eR) + " ".join("%8.4f" % r for r in rs))
            mach.append("R %s %s all %.5f edge %.5f err %.5f slices %s" % (vname.split(" (")[0].replace(" ", "_"), cname, R_all, R_edge, eR,
                                                                         " ".join("%.5f" % r for r in rs)))
    P("")
    P("WP8 factors for comparison (all M4): toy all 0.880, E1 0.889, E2 0.868; closure all 0.891, E1 0.869, E2 0.896")
    with open(OUT1, "w") as f:
        f.write("\n".join(lines) + "\n")
        f.write("\n".join(mach) + "\n")
    with open(OUT2, "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
