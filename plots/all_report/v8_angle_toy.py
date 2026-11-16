#!/usr/bin/env python3
# v8_angle_toy.py -- independent verifier (v8), 2026-11-14. Own toy, written from scratch (no WP8 code or inputs used):
# effect of the gamma* -> e+e- decay angle (isotropic as in HGeant gdalet.F vs QED 1 + cos^2 + (4 m_e^2/q^2) sin^2) on the
# acceptance of eta -> pi+pi-pi0, pi0 -> gamma e+e- in a HADES-like acceptance (18-85 deg, p_pi > 0.10, p_e > p_min GeV/c).
# eta kinematics: m_T-exponential (T = 0.10 GeV) x Gaussian rapidity around y_cm = 1.2214 (sigma_y 0.55) -- a simple stand-in.
# eta -> 3 pi: flat Dalitz plot; pi0 Dalitz pair mass: Kroll-Wada x (1 + a q^2/m^2)^2, a = 0.0332.
# Run: python3 -I plots/all_report/v8_angle_toy.py      (about 30 s)
import numpy as np
rng = np.random.default_rng(20261114)
ME, MPI, MPI0, META = 0.000511, 0.13957, 0.134977, 0.547862
YCM = 1.2214

def boost(p, b):          # p: (N,4) E,px,py,pz ; b: (N,3) velocity
    b2 = (b * b).sum(1); g = 1.0 / np.sqrt(1 - b2)
    bp = (b * p[:, 1:]).sum(1)
    g2 = np.where(b2 > 0, (g - 1) / np.where(b2 > 0, b2, 1), 0.0)
    E = g * (p[:, 0] + bp)
    v = p[:, 1:] + (g2 * bp + g * p[:, 0])[:, None] * b
    return np.column_stack([E, v])

def iso_dirs(n):
    c = rng.uniform(-1, 1, n); ph = rng.uniform(0, 2 * np.pi, n); s = np.sqrt(1 - c * c)
    return np.column_stack([s * np.cos(ph), s * np.sin(ph), c])

def two_body(M, m1, m2, dirs):
    p = np.sqrt(np.maximum((M * M - (m1 + m2) ** 2) * (M * M - (m1 - m2) ** 2), 0)) / (2 * M)
    E1 = np.sqrt(p * p + m1 * m1); E2 = np.sqrt(p * p + m2 * m2)
    v = dirs * p[:, None]
    return np.column_stack([E1, v]), np.column_stack([E2, -v])

def sample_kw(n):
    # Kroll-Wada in q: (2/q)(1-x)^3 (1+2me^2/q^2) sqrt(1-4me^2/q^2) FF^2, x = q^2/m^2; sample ln q uniformly, accept-reject
    lo, hi = np.log(2 * ME * 1.0000001), np.log(MPI0)
    out = np.empty(0)
    while out.size < n:
        lq = rng.uniform(lo, hi, 4 * n); q = np.exp(lq); x = (q / MPI0) ** 2
        w = 2 * (1 - x) ** 3 * (1 + 2 * ME * ME / q ** 2) * np.sqrt(np.clip(1 - 4 * ME * ME / q ** 2, 0, None)) * (1 + 0.0332 * x) ** 2
        out = np.concatenate([out, q[rng.uniform(0, 2.0, q.size) < w]])
    return out[:n]

def three_body_flat(n):
    # eta -> pi+ pi- pi0 uniform in the Dalitz plot (m12^2, m23^2)
    res = []
    tot = 0
    while tot < n:
        k = 4 * n
        s12 = rng.uniform((MPI + MPI) ** 2, (META - MPI0) ** 2, k)
        s23 = rng.uniform((MPI + MPI0) ** 2, (META - MPI) ** 2, k)
        m12 = np.sqrt(s12)
        E2s = (s12 - MPI * MPI + MPI * MPI) / (2 * m12); E3s = (META ** 2 - s12 - MPI0 ** 2) / (2 * m12)
        p2 = np.sqrt(np.clip(E2s ** 2 - MPI ** 2, 0, None)); p3 = np.sqrt(np.clip(E3s ** 2 - MPI0 ** 2, 0, None))
        lo = (E2s + E3s) ** 2 - (p2 + p3) ** 2; hi = (E2s + E3s) ** 2 - (p2 - p3) ** 2
        ok = (s23 > lo) & (s23 < hi) & (E3s > MPI0)
        res.append(s12[ok]); tot += ok.sum()
    s12 = np.concatenate(res)[:n]
    m12 = np.sqrt(s12)
    # eta rest frame: (pi+pi-) system vs pi0, then pi+ pi- in their frame, isotropic
    pp, p0 = two_body(np.full(n, META), m12, np.full(n, MPI0), iso_dirs(n))
    a, b = two_body(m12, np.full(n, MPI), np.full(n, MPI), iso_dirs(n))
    bet = pp[:, 1:] / pp[:, :1]
    return boost(a, bet), boost(b, bet), p0

def run(n, pmins=(0.04, 0.06, 0.10)):
    # eta in the lab
    mT = META + rng.exponential(0.10, n); pT = np.sqrt(mT ** 2 - META ** 2)
    y = rng.normal(YCM, 0.55, n); ph = rng.uniform(0, 2 * np.pi, n)
    eta = np.column_stack([mT * np.cosh(y), pT * np.cos(ph), pT * np.sin(ph), mT * np.sinh(y)])
    beta_eta = eta[:, 1:] / eta[:, :1]
    pip, pim, pi0 = three_body_flat(n)
    # pi0 -> gamma gamma* in the pi0 frame
    q = sample_kw(n)
    gam, gst = two_body(np.full(n, MPI0), np.zeros(n), q, iso_dirs(n))
    # gamma* -> e+ e- in the gamma* frame; angle relative to the gamma* direction (in the pi0 frame)
    u = gst[:, 1:] / np.linalg.norm(gst[:, 1:], axis=1)[:, None]
    # orthonormal basis around u
    t = np.where(np.abs(u[:, 2:3]) < 0.9, np.array([[0, 0, 1.0]]), np.array([[1.0, 0, 0]]))
    e1 = np.cross(u, t); e1 /= np.linalg.norm(e1, axis=1)[:, None]; e2 = np.cross(u, e1)
    res = {}
    eps_iso = 4 * ME * ME / q ** 2
    for mode in ("iso", "qed"):
        if mode == "iso":
            c = rng.uniform(-1, 1, n)
        else:
            c = np.empty(n); todo = np.ones(n, bool)
            while todo.any():
                k = todo.sum(); cc = rng.uniform(-1, 1, k)
                w = (1 + cc * cc + eps_iso[todo] * (1 - cc * cc)) / 2.0
                acc = rng.uniform(0, 1, k) < w
                idx = np.where(todo)[0][acc]; c[idx] = cc[acc]; todo[idx] = False
        phi = rng.uniform(0, 2 * np.pi, n); s = np.sqrt(1 - c * c)
        d = c[:, None] * u + (s * np.cos(phi))[:, None] * e1 + (s * np.sin(phi))[:, None] * e2
        ep, em = two_body(q, np.full(n, ME), np.full(n, ME), d)
        bg = gst[:, 1:] / gst[:, :1]
        ep, em = boost(ep, bg), boost(em, bg)                    # -> pi0 frame
        b0 = pi0[:, 1:] / pi0[:, :1]
        ep, em = boost(ep, b0), boost(em, b0)                    # -> eta frame
        lab = [boost(v, beta_eta) for v in (pip, pim, ep, em)]   # -> lab
        def th_p(v):
            p = np.linalg.norm(v[:, 1:], axis=1); return np.degrees(np.arccos(v[:, 3] / p)), p
        acc_pi = np.ones(n, bool)
        for v in lab[:2]:
            th, p = th_p(v); acc_pi &= (th > 18) & (th < 85) & (p > 0.10)
        (the, pe), (thm, pm) = th_p(lab[2]), th_p(lab[3])
        for pmin in pmins:
            acc = acc_pi & (the > 18) & (the < 85) & (thm > 18) & (thm < 85) & (pe > pmin) & (pm > pmin)
            res[(mode, pmin)] = acc
    return q, res

if __name__ == "__main__":
    N = 2_000_000
    tot = {}
    qs = []
    for it in range(3):
        q, r = run(N)
        qs.append(q)
        for k, v in r.items():
            tot.setdefault(k, []).append(v)
    q = np.concatenate(qs)
    print("# v8_angle_toy: eps(QED angle)/eps(isotropic), %d decays; classes in true M(ee)" % q.size)
    print("%-10s %8s %8s %8s %8s %8s   %s" % ("p_e min", "all", "<10", ">=10", ">=20", "20-30", "eps_iso(all)"))
    for pmin in (0.04, 0.06, 0.10):
        a_iso = np.concatenate(tot[("iso", pmin)]); a_qed = np.concatenate(tot[("qed", pmin)])
        row = []
        for sel in (np.ones_like(q, bool), q < 0.010, q >= 0.010, q >= 0.020, (q >= 0.020) & (q < 0.030)):
            ni, nq = a_iso[sel].sum(), a_qed[sel].sum()
            row.append(nq / ni)
        print("%-10.2f %8.4f %8.4f %8.4f %8.4f %8.4f   %.4f" % (pmin, *row, a_iso.mean()))
    print("statistical precision of a ratio: about %.4f (all), %.4f (>=20)" % (
        1 / np.sqrt(np.concatenate(tot[("iso", 0.06)]).sum()), 1 / np.sqrt(np.concatenate(tot[("iso", 0.06)])[q >= 0.02].sum())))
