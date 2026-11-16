#!/usr/bin/env python3
# ============================================================================
# wp8_toy_angle.py  --  WP8 (report_all): toy estimate of the gamma* -> e+e- decay-angle effect on class C2
# ============================================================================
# Written by WP8 (prefix wp8), 2026-11-14. A TOY, not a detector simulation.
#
# Why: the C2 generator (HGeant gdalet.F) samples the pair mass from Kroll-Wada correctly (wp8_qed.py) but decays
# the virtual photon ISOTROPICALLY in its rest frame; QED gives dN/dcos(theta) ~ 1 + cos^2 + (1 - beta^2) sin^2,
# theta = lepton angle to the gamma* direction (pi0 rest frame), beta = lepton velocity in the gamma* frame.
# More asymmetric pairs (|cos| -> 1) mean a softer lepton, hence a slightly lower acceptance.
#
# What the toy does (vectorised numpy, same events for both angle models = weights):
#   eta (pT, y_cm) from the SMASH standard-set generator histogram (/hdd3/.../wp8/remote/gen_eta_ptY_std40.txt, 40 files) at
#   T = 4500 MeV; eta -> pi+ pi- pi0 phase space (KLOE-2 Dalitz weight applied); pi0 -> gamma gamma*,
#   q from Kroll-Wada x (1 + a q^2/m^2)^2 (a = 0.0332); gamma* -> e+e- isotropic (HGeant) or QED (weight).
#   Acceptance: all four tracks 18 < theta_lab < 85 deg; pions p > 0.10 GeV/c; leptons p > p_min (0.06 nominal,
#   variants 0.04 / 0.10 GeV/c); variant: + OA(e+e-) > 2 deg.
# Output: docs/latex/tables_all/wp8_toy_angle.txt: eps(QED)/eps(isotropic) per true M(ee) class and the shares.
# Run: python3 -I plots/all_report/wp8_toy_angle.py
# ============================================================================
import math
import numpy as np

GEN = "/hdd3/przygoda/HADES/FAT_out/all_report/wp8/remote/gen_eta_ptY_std40.txt"   # copy of the output of wp8/remote/remote_eta_ptY.sh
OUT = "/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all/wp8_toy_angle.txt"
ME, MPI, MPI0, META, MP = 0.000510999, 0.13957039, 0.1349768, 0.547862, 0.938272
A_FF = 0.0332
rng = np.random.default_rng(20261114)


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
    """boost four-vectors p4 (n,4: E,px,py,pz) by velocity b (n,3)"""
    b2 = np.sum(b * b, axis=1)
    g = 1.0 / np.sqrt(1.0 - b2)
    bp = np.sum(b * p4[:, 1:], axis=1)
    g2 = np.where(b2 > 0, (g - 1.0) / np.where(b2 > 0, b2, 1.0), 0.0)
    E = g * (p4[:, 0] + bp)
    vec = p4[:, 1:] + (g2 * bp)[:, None] * b + (g * p4[:, 0])[:, None] * b
    return np.column_stack([E, vec])


def two_body(M, m1, m2, n, axis=None):
    """momentum of particle 1 in the rest frame of M (isotropic unless axis/cos given)"""
    p = np.sqrt(np.maximum(0, (M * M - (m1 + m2) ** 2) * (M * M - (m1 - m2) ** 2))) / (2 * M)
    return p


def kloe_weight(mpp2, s13, s23):
    # KLOE-2 2016 eta -> pi+pi-pi0 Dalitz parameters; X = sqrt3 (T+ - T-)/Q, Y = 3 T0/Q - 1
    Q = META - 2 * MPI - MPI0
    E0 = (META ** 2 + MPI0 ** 2 - mpp2) / (2 * META)
    T0 = E0 - MPI0
    Ep = (META ** 2 + MPI ** 2 - s23) / (2 * META)   # s23 = M(pi- pi0)^2 -> pi+ energy
    Em = (META ** 2 + MPI ** 2 - s13) / (2 * META)
    Tp, Tm = Ep - MPI, Em - MPI
    X = math.sqrt(3) * (Tp - Tm) / Q
    Y = 3 * T0 / Q - 1
    a, b, d, f, g = -1.095, 0.145, 0.081, 0.141, -0.044
    return 1 + a * Y + b * Y * Y + d * X * X + f * Y ** 3 + g * X * X * Y


def sample_q(n):
    """Kroll-Wada x (1 + a x)^2 in q [GeV], by rejection on a log-uniform proposal"""
    out = np.empty(0)
    qmin, qmax = 2 * ME * 1.000001, MPI0 * 0.999999
    while out.size < n:
        m = 4 * (n - out.size) + 1000
        u = rng.random(m)
        q = qmin * (qmax / qmin) ** u          # density ~ 1/q
        x = (q / MPI0) ** 2
        beta = np.sqrt(1 - 4 * ME * ME / (q * q))
        f = (1 - x) ** 3 * (1 + 2 * ME * ME / (q * q)) * beta * (1 + A_FF * x) ** 2   # KW x q (divided by the 1/q proposal)
        keep = rng.random(m) * 2.0 < f
        out = np.concatenate([out, q[keep]])
    return out[:n]


def generate(n, h):
    pt, ycm = sample_eta(n, h)
    # beam 4500 MeV: y_cm of the initial state
    Eb = 4.5 + MP
    pzb = math.sqrt(Eb * Eb - MP * MP)
    yc = 0.5 * math.log((Eb + MP + pzb) / (Eb + MP - pzb))
    y = ycm + yc
    phi = 2 * np.pi * rng.random(n)
    mt = np.sqrt(META ** 2 + pt * pt)
    eta = np.column_stack([mt * np.cosh(y), pt * np.cos(phi), pt * np.sin(phi), mt * np.sinh(y)])
    # eta -> pi+ pi- pi0: Dalitz plot uniform (rejection), KLOE weight
    s12lo, s12hi = (2 * MPI) ** 2, (META - MPI0) ** 2
    s12 = np.empty(0); s13 = np.empty(0)
    while s12.size < n:
        m = 3 * (n - s12.size) + 1000
        a = s12lo + (s12hi - s12lo) * rng.random(m)
        b = (MPI + MPI0) ** 2 + ((META - MPI) ** 2 - (MPI + MPI0) ** 2) * rng.random(m)
        # Dalitz boundary for s13 at given s12
        E2s = np.sqrt(a) / 2                                   # pi- energy in the (pi+pi-) frame
        E3s = (META ** 2 - a - MPI0 ** 2) / (2 * np.sqrt(a))    # pi0 energy in the (pi+pi-) frame
        p2s = np.sqrt(np.maximum(0, E2s ** 2 - MPI ** 2)); p3s = np.sqrt(np.maximum(0, E3s ** 2 - MPI0 ** 2))
        lo = (E2s + E3s) ** 2 - (p2s + p3s) ** 2
        hi = (E2s + E3s) ** 2 - (p2s - p3s) ** 2
        ok = (b >= lo) & (b <= hi)
        s12 = np.concatenate([s12, a[ok]]); s13 = np.concatenate([s13, b[ok]])
    s12, s13 = s12[:n], s13[:n]
    s23 = META ** 2 + 2 * MPI ** 2 + MPI0 ** 2 - s12 - s13
    wk = kloe_weight(s12, s13, s23)
    # pi0 in the eta frame
    E0 = (META ** 2 + MPI0 ** 2 - s12) / (2 * META)
    p0 = np.sqrt(np.maximum(0, E0 ** 2 - MPI0 ** 2))
    d0 = iso_dir(n)
    pi0 = np.column_stack([E0, p0[:, None] * d0])
    pp = np.column_stack([META - E0, -p0[:, None] * d0])          # pi+pi- system
    mpp = np.sqrt(s12)
    # pions in the pi+pi- frame (direction isotropic; the Dalitz-plot angle is not needed for this toy)
    kp = np.sqrt(np.maximum(0, mpp ** 2 / 4 - MPI ** 2))
    dpi = iso_dir(n)
    pip = np.column_stack([np.sqrt(kp ** 2 + MPI ** 2), kp[:, None] * dpi])
    pim = np.column_stack([np.sqrt(kp ** 2 + MPI ** 2), -kp[:, None] * dpi])
    bpp = pp[:, 1:] / pp[:, :1]
    pip = boost(pip, bpp); pim = boost(pim, bpp)
    # pi0 -> gamma gamma*: gamma* direction isotropic in the pi0 frame
    q = sample_q(n)
    ps = (MPI0 ** 2 - q * q) / (2 * MPI0)
    dg = iso_dir(n)
    gs = np.column_stack([np.sqrt(ps * ps + q * q), ps[:, None] * dg])
    # gamma* -> e+ e-: lepton angle to the gamma* direction, isotropic; QED weight
    k = np.sqrt(np.maximum(0, q * q / 4 - ME * ME))
    beta = 2 * k / q
    c = 2 * rng.random(n) - 1
    s = np.sqrt(1 - c * c)
    ph = 2 * np.pi * rng.random(n)
    # local frame around dg
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
    # pi0 frame -> eta frame -> lab
    b0 = pi0[:, 1:] / pi0[:, :1]
    ep = boost(ep, b0); em = boost(em, b0)
    be = eta[:, 1:] / eta[:, :1]
    ep, em, pip, pim = boost(ep, be), boost(em, be), boost(pip, be), boost(pim, be)
    return dict(q=q, wk=wk, wq=wq, ep=ep, em=em, pip=pip, pim=pim)


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


def main():
    h = load_gen()
    n = 3_000_000
    ev = generate(n, h)
    q = ev["q"]
    classes = [("all", 0, 1), ("<10", 0, 0.010), (">=10", 0.010, 1), ("<20", 0, 0.020), (">=20", 0.020, 1),
               ("20-30", 0.020, 0.030), (">=30", 0.030, 1)]
    lines = []
    def P(s=""):
        print(s); lines.append(s)
    P("# wp8_toy_angle.txt -- plots/all_report/wp8_toy_angle.py (WP8, 2026-11-14). TOY, not a detector simulation.")
    P("# eta (pT, y_cm) from the SMASH std generator (40 files), eta -> pi+pi-pi0 phase space x KLOE-2, pi0 -> gamma e+e- Kroll-Wada x PDG FF;")
    P("# acceptance: 4 tracks 18-85 deg, pions p > 0.10 GeV/c, leptons p > p_min; classes in TRUE M(ee).")
    P("# eps ratio = acceptance with the QED decay angle / acceptance with the isotropic angle of HGeant (same q spectrum)")
    P("# share(>=20) = accepted fraction with M(ee) >= 20 MeV among all accepted; N = %d toy decays" % n)
    P("")
    P("%-26s" % "variant" + "".join("%10s" % c[0] for c in classes) + "%14s%14s" % ("sh>=20 QED/iso", "sh>=20|>=10"))
    res = {}
    for name, pmin, oa in (("p_e > 0.06 (nominal)", 0.06, 0), ("p_e > 0.04", 0.04, 0), ("p_e > 0.10", 0.10, 0),
                           ("p_e > 0.06, OA > 2 deg", 0.06, 2.0)):
        acc = accept(ev, pmin, oa)
        row = []
        for cname, lo, hi in classes:
            m = acc & (q >= lo) & (q < hi)
            e_iso = np.sum(ev["wk"][m]) ; e_qed = np.sum(ev["wk"][m] * ev["wq"][m])
            row.append(e_qed / e_iso if e_iso > 0 else float("nan"))
        a_all = acc
        s_iso = np.sum(ev["wk"][a_all & (q >= 0.02)]) / np.sum(ev["wk"][a_all])
        s_qed = np.sum((ev["wk"] * ev["wq"])[a_all & (q >= 0.02)]) / np.sum((ev["wk"] * ev["wq"])[a_all])
        s10_iso = np.sum(ev["wk"][a_all & (q >= 0.02)]) / np.sum(ev["wk"][a_all & (q >= 0.01)])
        s10_qed = np.sum((ev["wk"] * ev["wq"])[a_all & (q >= 0.02)]) / np.sum((ev["wk"] * ev["wq"])[a_all & (q >= 0.01)])
        eff = np.sum(ev["wk"][a_all]) / np.sum(ev["wk"])
        res[name] = (row, s_qed / s_iso, s10_qed / s10_iso, eff)
        P("%-26s" % name + "".join("%10.4f" % x for x in row) + "%14.4f%14.4f" % (s_qed / s_iso, s10_qed / s10_iso))
    P("")
    P("toy acceptance (isotropic): " + ", ".join("%s %.4f" % (k, v[3]) for k, v in res.items()))
    # statistical precision of the ratio (weights are on the same events): bootstrap-free estimate via the weight spread
    acc = accept(ev, 0.06, 0)
    m = acc & (q >= 0.02)
    w = ev["wq"][m] * ev["wk"][m] / np.mean(ev["wk"][m])
    P("statistical precision of the >=20 MeV ratio (nominal): +- %.4f" % (np.std(w) / math.sqrt(m.sum()) / np.mean(w)))
    with open(OUT, "w") as f:
        f.write("\n".join(lines) + "\n")
        for k, v in res.items():
            f.write("toy %s | " % k.replace(" ", "_") + " ".join("%s=%.5f" % (c[0], x) for c, x in zip(classes, v[0])) +
                    " share20=%.5f share20_10=%.5f eff=%.5f\n" % (v[1], v[2], v[3]))


if __name__ == "__main__":
    main()
