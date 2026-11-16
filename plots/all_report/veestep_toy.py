#!/usr/bin/env python3
# ============================================================================
# veestep_toy.py  --  INDEPENDENT verification of eestep (prefix veestep): kinematic toy of cross-photon lepton pairs
# ============================================================================
# Written 2026-11-15 by the independent verifier (own generator; eestep_kintoy.py not used or read).
# eta -> pi+ pi- pi0 (3-body phase space), pi0 -> gamma gamma or gamma e+e- (Kroll-Wada q^2, 1 + cos^2 decay angle),
# photon conversion gamma -> e+ e- with the Bethe-Heitler fraction (complete screening, dsigma/dx ~ 1 - 4/3 x(1 - x)),
# collinear and massless conversion leptons. For every lepton pair (true pairs = the bump classes C2, C3 for reference;
# cross pairs = one lepton of each photon) the true u = M4 - M4max(M(pipi), M(ee)) is computed.
# Lab frame: p p -> p p eta at T = 4.5 GeV (3-body phase space), crude HADES acceptance (18-85 deg, p > 80 MeV/c for
# leptons, p > 100 MeV/c for pions), charge-symmetric (the e+/e- difference of the real detector is not modelled).
# Answers: (1) the largest u of every topology (the edge), (2) the shape of dN/du near the edge (step or ramp),
# (3) M(ee) and opening angle of the pairs near the edge, (4) what the data estimator (pol2 + step, sigma 6.7 MeV)
# returns for such a shape.
# Run: python3 -I plots/all_report/veestep_toy.py  (-> /hdd3/przygoda/HADES/FAT_out/all_report/veestep/veestep_toy.log)
# ============================================================================
import math
import sys

import numpy as np

sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import veestep_core as V

ME_ = 0.000511
MPIC = 0.13957
MP = 0.938272
MEta, MPi0 = V.MEta, V.MPi0
rng = np.random.default_rng(20261115)
LOG = open(V.OUT + "veestep_toy.log", "w")


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s); LOG.write(s + "\n"); LOG.flush()


def iso(n):
    c = rng.uniform(-1, 1, n); ph = rng.uniform(0, 2 * np.pi, n); s = np.sqrt(1 - c * c)
    return np.stack([s * np.cos(ph), s * np.sin(ph), c], 1)


def boost(p4, b):
    """boost 4-vectors p4 (n x 4: E, px, py, pz) by velocity b (n x 3)"""
    b2 = (b * b).sum(1)
    g = 1 / np.sqrt(1 - b2)
    bp = (b * p4[:, 1:]).sum(1)
    g2 = np.where(b2 > 0, (g - 1) / np.where(b2 > 0, b2, 1), 0)
    E = g * (p4[:, 0] + bp)
    p = p4[:, 1:] + (g2 * bp + g * p4[:, 0])[:, None] * b
    return np.concatenate([E[:, None], p], 1)


def two_body(M, m1, m2, n, dirs=None):
    M = np.broadcast_to(np.asarray(M, float), (n,)); m1 = np.broadcast_to(np.asarray(m1, float), (n,))
    m2 = np.broadcast_to(np.asarray(m2, float), (n,))
    p = np.sqrt(np.maximum((M * M - (m1 + m2) ** 2) * (M * M - (m1 - m2) ** 2), 0)) / (2 * M)
    d = iso(n) if dirs is None else dirs
    E1 = np.sqrt(p * p + m1 * m1); E2 = np.sqrt(p * p + m2 * m2)
    a = np.concatenate([E1[:, None], (p[:, None] * d)], 1)
    b = np.concatenate([E2[:, None], (-p[:, None] * d)], 1)
    return a, b


def mass(p4):
    return np.sqrt(np.maximum(p4[:, 0] ** 2 - (p4[:, 1:] ** 2).sum(1), 0))


def three_body_dalitz(M, m1, m2, m3, n):
    """uniform Dalitz-plot sampling: returns s12 (= m(12)^2) and s23"""
    out12 = []; out23 = []
    need = n
    while need > 0:
        k = int(need * 2.5) + 100
        s12 = rng.uniform((m1 + m2) ** 2, (M - m3) ** 2, k)
        s23 = rng.uniform((m2 + m3) ** 2, (M - m1) ** 2, k)
        m12 = np.sqrt(s12)
        E2 = (s12 - m1 * m1 + m2 * m2) / (2 * m12); E3 = (M * M - s12 - m3 * m3) / (2 * m12)
        p2 = np.sqrt(np.maximum(E2 * E2 - m2 * m2, 0)); p3 = np.sqrt(np.maximum(E3 * E3 - m3 * m3, 0))
        lo = (E2 + E3) ** 2 - (p2 + p3) ** 2; hi = (E2 + E3) ** 2 - (p2 - p3) ** 2
        ok = (s23 >= lo) & (s23 <= hi)
        out12.append(s12[ok]); out23.append(s23[ok]); need -= ok.sum()
    return np.concatenate(out12)[:n], np.concatenate(out23)[:n]


def kroll_wada_q(n, M=MPi0):
    """q from dGamma/dq^2 ~ (1/q^2)(1 - q^2/M^2)^3 (1 + 2 m^2/q^2) sqrt(1 - 4 m^2/q^2); sampled in ln q^2"""
    lo, hi = math.log(4 * ME_ * ME_ * 1.0001), math.log(M * M)
    out = []; need = n
    while need > 0:
        k = int(need * 4) + 100
        l = rng.uniform(lo, hi, k); q2 = np.exp(l)
        f = (1 - q2 / M ** 2) ** 3 * (1 + 2 * ME_ ** 2 / q2) * np.sqrt(np.maximum(1 - 4 * ME_ ** 2 / q2, 0))
        ok = rng.uniform(0, 1.0001, k) < f
        out.append(np.sqrt(q2[ok])); need -= ok.sum()
    return np.concatenate(out)[:n]


def bh_x(n):
    out = []; need = n
    while need > 0:
        k = int(need * 1.6) + 100
        x = rng.uniform(0, 1, k)
        ok = rng.uniform(0, 1, k) < (1 - 4.0 / 3.0 * x * (1 - x))
        out.append(x[ok]); need -= ok.sum()
    return np.concatenate(out)[:n]


def gen(n):
    """one eta per entry, eta rest frame: P_pipi (n x 4), pions, photon 1 = gamma* (Dalitz) or gamma, photon 2 = gamma"""
    s12, s23 = three_body_dalitz(MEta, MPIC, MPIC, MPi0, n)    # 1 = pi+, 2 = pi-, 3 = pi0
    mpp = np.sqrt(s12)
    pi0, pipi = two_body(MEta, MPi0, mpp, n)
    # pions inside the pipi system (isotropic, enough for the lab acceptance)
    pip, pim = two_body(mpp, MPIC, MPIC, n)
    bpp = pipi[:, 1:] / pipi[:, :1]
    pip = boost(pip, bpp); pim = boost(pim, bpp)
    return mpp, pipi, pip, pim, pi0


def pi0_gg(pi0):
    n = pi0.shape[0]
    g1, g2 = two_body(MPi0, 0.0, 0.0, n)
    b = pi0[:, 1:] / pi0[:, :1]
    return boost(g1, b), boost(g2, b)


def pi0_dalitz(pi0):
    """returns e+ (Dalitz), e- (Dalitz), real photon; gamma* angle 1 + cos^2"""
    n = pi0.shape[0]
    q = kroll_wada_q(n)
    gs, g = two_body(MPi0, q, np.zeros(n), n)
    # lepton angle in the gamma* frame relative to the gamma* direction: 1 + cos^2 (+ mass term neglected)
    c = []; need = n
    while need > 0:
        k = int(need * 1.6) + 100
        x = rng.uniform(-1, 1, k); ok = rng.uniform(0, 2, k) < 1 + x * x
        c.append(x[ok]); need -= ok.sum()
    c = np.concatenate(c)[:n]
    ax = gs[:, 1:] / np.linalg.norm(gs[:, 1:], axis=1)[:, None]
    tmp = np.where(np.abs(ax[:, 2:3]) < 0.9, np.array([[0, 0, 1.0]]), np.array([[1.0, 0, 0]]))
    e1 = np.cross(ax, tmp); e1 /= np.linalg.norm(e1, axis=1)[:, None]; e2 = np.cross(ax, e1)
    ph = rng.uniform(0, 2 * np.pi, n); s = np.sqrt(1 - c * c)
    d = c[:, None] * ax + (s * np.cos(ph))[:, None] * e1 + (s * np.sin(ph))[:, None] * e2
    ep, em = two_body(q, ME_, ME_, n, dirs=d)
    bgs = gs[:, 1:] / gs[:, :1]
    ep = boost(ep, bgs); em = boost(em, bgs)
    b = pi0[:, 1:] / pi0[:, :1]
    return boost(ep, b), boost(em, b), boost(g, b), q


def convert(g):
    x = bh_x(g.shape[0])
    em = g * x[:, None]; ep = g * (1 - x)[:, None]
    return ep, em, x


def u_true(pipi, mpp, a, b):
    Q = a + b
    q = mass(Q)
    m4 = mass(pipi + Q)
    return m4 - V.m4max(mpp, q), q


def lab_boost_eta(n):
    """eta 4-momentum in the lab for p p -> p p eta at T = 4.5 GeV (3-body phase space); returns velocity of the eta (n x 3)"""
    T = 4.5
    pb = math.sqrt(T * T + 2 * T * MP); Eb = T + MP
    sq = math.sqrt((Eb + MP) ** 2 - pb * pb)
    s12, s23 = three_body_dalitz(sq, MP, MP, MEta, n)    # 3 = eta
    m12 = np.sqrt(s12)
    eta, pp = two_body(sq, MEta, m12, n)
    bcm = np.tile(np.array([[0, 0, pb / (Eb + MP)]]), (n, 1))
    eta = boost(eta, bcm)
    return eta[:, 1:] / eta[:, :1]


def theta_p(p4):
    p = np.linalg.norm(p4[:, 1:], axis=1)
    th = np.degrees(np.arccos(np.clip(p4[:, 3] / np.maximum(p, 1e-12), -1, 1)))
    return th, p


def acc(p4, pmin):
    th, p = theta_p(p4)
    return (th > 18) & (th < 85) & (p > pmin)


def oa(a, b):
    pa = a[:, 1:]; pb = b[:, 1:]
    c = (pa * pb).sum(1) / np.maximum(np.linalg.norm(pa, axis=1) * np.linalg.norm(pb, axis=1), 1e-12)
    return np.degrees(np.arccos(np.clip(c, -1, 1)))


def shape(u, w, lo=-0.060):
    e = np.arange(lo, 0.0001, 0.0025)
    h, _ = np.histogram(u, bins=e, weights=w)
    return e, h


def main():
    n = 600000
    mpp, pipi, pip, pim, pi0 = gen(n)
    # --- topology B: pi0 Dalitz (gamma* -> e+e-) + conversion of the real photon
    Dp, Dm, g, q_d = pi0_dalitz(pi0)
    Cp, Cm, x = convert(g)
    # --- topology A: two conversions
    g1, g2 = pi0_gg(pi0)
    Ap, Am, x1 = convert(g1)
    Bp, Bm, x2 = convert(g2)
    # lab
    beta = lab_boost_eta(n)
    L = {nm: boost(v, beta) for nm, v in dict(Dp=Dp, Dm=Dm, Cp=Cp, Cm=Cm, Ap=Ap, Am=Am, Bp=Bp, Bm=Bm, pip=pip, pim=pim, g=g, g1=g1, g2=g2, pipi=pipi).items()}
    pions = acc(L["pip"], 0.100) & acc(L["pim"], 0.100)
    det = {k: acc(L[k], 0.080) for k in ["Dp", "Dm", "Cp", "Cm", "Ap", "Am", "Bp", "Bm"]}
    pairs = {
        "C2 bump (D+ D-, photon lost)": ("Dp", "Dm", None),
        "C3 bump (A+ A-, other photon lost)": ("Ap", "Am", None),
        "B like-sign (D- C-)": ("Dm", "Cm", None), "B like-sign (D+ C+)": ("Dp", "Cp", None),
        "B unlike (D+ C-)": ("Dp", "Cm", None), "B unlike (D- C+)": ("Dm", "Cp", None),
        "A like-sign (A- B-)": ("Am", "Bm", None), "A like-sign (A+ B+)": ("Ap", "Bp", None),
        "A unlike (A+ B-)": ("Ap", "Bm", None), "A unlike (A- B+)": ("Am", "Bp", None)}
    P("# veestep_toy.py: %d eta -> pi+pi-pi0 per topology (independent generator)" % n)
    P("\n== 1. largest true u per pair type (no acceptance) [MeV/c2]; defs.h limit without scale factors")
    store = {}
    for nm, (a, b, _) in pairs.items():
        u, q = u_true(L["pipi"], mpp, L[a], L[b])     # invariant: lab = eta frame
        u0, _ = u_true(pipi, mpp, {"Dp": Dp, "Dm": Dm, "Cp": Cp, "Cm": Cm, "Ap": Ap, "Am": Am, "Bp": Bp, "Bm": Bm}[a],
                       {"Dp": Dp, "Dm": Dm, "Cp": Cp, "Cm": Cm, "Ap": Ap, "Am": Am, "Bp": Bp, "Bm": Bm}[b])
        store[nm] = (u, q, a, b)
        ok = np.isfinite(u)
        P("  %-44s max u %.4f MeV (eta frame %.4f), fraction with u > 0: %.1e, max |u_lab - u_eta| %.1e MeV" %
          (nm, 1000 * np.nanmax(u), 1000 * np.nanmax(u0), np.mean(u[ok] > 1e-7), 1000 * np.nanmax(np.abs(u - u0))))
    P("\n== 2. shape near the edge, with the crude acceptance (pions and the two leptons of the pair detected)")
    P("  density in 5 MeV bins from u = -40 to 0 MeV, relative to the bin -40..-35 (a step: flat to the edge; a ramp: falls to 0)")
    for nm, (u, q, a, b) in store.items():
        sel = pions & det[a] & det[b] & np.isfinite(u)
        e = np.arange(-0.040, 0.0001, 0.005)
        h, _ = np.histogram(u[sel], bins=e)
        rel = h / max(h[0], 1)
        P("  %-44s N(acc, |u|<40) %6d: " % (nm, h.sum()) + " ".join("%.2f" % r for r in rel) +
          "  | M(ee) median at u > -5 MeV: %.0f MeV, OA median %.0f deg" %
          (1000 * np.median(q[sel & (u > -0.005)]) if np.any(sel & (u > -0.005)) else np.nan,
           np.median(oa(L[a], L[b])[sel & (u > -0.005)]) if np.any(sel & (u > -0.005)) else np.nan))
    P("\n== 3. cross pairs near the edge (u > -5 MeV/c2, accepted): M(ee) and opening-angle quantiles")
    for nm in ["B like-sign (D- C-)", "A like-sign (A- B-)", "C2 bump (D+ D-, photon lost)", "C3 bump (A+ A-, other photon lost)"]:
        u, q, a, b = store[nm]
        sel = pions & det[a] & det[b] & np.isfinite(u) & (u > -0.005)
        o = oa(L[a], L[b])[sel]
        P("  %-44s rows %5d: M(ee) 10/50/90 %% = %.0f / %.0f / %.0f MeV; share M(ee) >= 90: %.2f; OA 10/50/90 %% = %.0f / %.0f / %.0f deg" %
          ((nm, sel.sum()) + tuple(1000 * np.percentile(q[sel], [10, 50, 90])) + (np.mean(q[sel] >= 0.090),) + tuple(np.percentile(o, [10, 50, 90]))))
    # partner content of the like-sign rows: are the e+ partners detected?
    for nm, pa in [("B like-sign (D- C-)", ["Dp", "Cp"]), ("A like-sign (A- B-)", ["Ap", "Bp"])]:
        u, q, a, b = store[nm]
        sel = pions & det[a] & det[b] & np.isfinite(u)
        near = sel & (u > -0.010)
        nopos = ~(det[pa[0]] | det[pa[1]])
        P("  %s: share with no detected e+ in the event: all u %.2f, u > -10 MeV %.2f" % (nm, nopos[sel].mean(), nopos[near].mean()))
    P("\n== 4. what the data estimator sees: shape (with the acceptance) smeared with sigma 6.7 MeV on a flat background")
    for nm in ["C2 bump (D+ D-, photon lost)", "C3 bump (A+ A-, other photon lost)", "B like-sign (D- C-)", "A like-sign (A- B-)"]:
        u, q, a, b = store[nm]
        sel = pions & det[a] & det[b] & np.isfinite(u) & (u > -0.12)
        if nm.startswith("B") or nm.startswith("A"):
            sel &= q >= 0.090
        us = u[sel] + rng.normal(0, 0.0067, sel.sum())
        h = V.hist(us)
        # scale so that the density at -45..-30 MeV equals 60 rows per bin, on a flat background of 300 rows per bin
        dens = h[(V.XC > -0.045) & (V.XC < -0.030)].mean()
        y = 300 + 60 * h / max(dens, 1e-9)
        r = V.fit_pois(y)
        prof = [(u0 * 1000, V.fit_pois(y, u0=u0)["chi2"]) for u0 in np.arange(-0.0125, 0.0051, 0.0025)]
        best = min(prof, key=lambda t: t[1])
        # true density just below the edge (-5..0 MeV) relative to -45..-30
        d0 = np.histogram(u[sel], bins=[-0.005, 0])[0][0] / 2 / max(np.histogram(u[sel], bins=[-0.045, -0.030])[0][0] / 6, 1e-9)
        P("  %-44s%s step h = %.1f (density below the edge was set to 60): h/60 = %.2f; best step position %+.1f MeV; true density -5..0 / -45..-30 = %.2f" %
          (nm, " (M(ee) >= 90)" if nm[0] in "AB" else "", r["h"], r["h"] / 60, best[0], d0))
    P("\n== 5. share of the cross-pair step per M(ee) band (estimator applied to the smeared accepted shape; linear, flat weights)")
    Lf = V.Lin(np.full(V.NB, 300.0))
    bands = [(0, 50), (50, 70), (70, 90), (90, 110), (110, 135)]
    for nm in ["B like-sign (D- C-)", "A like-sign (A- B-)"]:
        u, q, a, b = store[nm]
        sel = pions & det[a] & det[b] & np.isfinite(u)
        us = u + rng.normal(0, 0.0067, u.size)
        hs = []
        for lo, hi in bands:
            m = sel & (q * 1000 >= lo) & (q * 1000 < hi)
            hs.append(Lf.h(V.hist(us[m])))
        tot = sum(hs)
        P("  %-24s " % nm + "  ".join("%d-%d: %.2f" % (lo, hi, h / tot) for (lo, hi), h in zip(bands, hs)))
    P("  data e-e- (veestep_ana, rows per bin): 0-50: 3.6, 50-70: 10.1, 70-90: 3.0, 90-110: 28.8, 110-135: 27.3 (sum 72.8)")
    P("  variants (B like-sign): M(pipi) 0.30-0.40 (data range / 0.9970) and lepton p_min 50 / 80 / 150 MeV/c")
    inr = (mpp >= 0.30 / 0.997) & (mpp < 0.40 / 0.997)
    u, q, a, b = store["B like-sign (D- C-)"]
    us = u + rng.normal(0, 0.0067, u.size)
    for pmin in [0.050, 0.080, 0.150]:
        sel = pions & acc(L[a], pmin) & acc(L[b], pmin) & np.isfinite(u) & inr
        hs = [Lf.h(V.hist(us[sel & (q * 1000 >= lo) & (q * 1000 < hi)])) for lo, hi in bands]
        tot = sum(hs)
        P("    p_min %3.0f MeV/c: " % (pmin * 1000) + "  ".join("%d-%d: %.2f" % (lo, hi, h / tot) for (lo, hi), h in zip(bands, hs)))
    P("done")


if __name__ == "__main__":
    main()
