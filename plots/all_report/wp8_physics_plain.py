#!/usr/bin/env python3
# wp8_physics_plain.py -- COPY of plots/all_report/wp8_physics.py made by the short-report figure agent, redo shortfig_r (integ6, 2026-11-16):
# internal catalogue labels in the figure text replaced by plain names (glossary_v2.md); outputs redirected:
# figures of the scope -> docs/latex/<dir>/<base>_plain*.{pdf,png}; everything else -> /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ6/shortfig_r/run/wp8
# Built by /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ6/shortfig_r/mk_plain_r.py; the physics code is unchanged.
# ============================================================================
# wp8_physics.py  --  WP8 (report_all): Table P1, the eps2 QED check, the S12 update, Table P2 and Fig. P1
# ============================================================================
# Written by WP8 (prefix wp8), 2026-11-14. Combination only: every number is read from a file written by a macro.
# Inputs (all read-only):
#   /hdd3/.../all_report/wp8/wp8_sim.txt          wp8_sim.C: simulated visible cross sections (KEY eff), ECAL tag per M(ee) class (KEY tag)
#   /hdd3/.../all_report/wp8/wp8_syst_raw.txt     wp8_syst.C: WP7 chain with the S12 variants (decay angle, (pT4, y4))
#   /hdd3/.../all_report/wp8/remote/pluto_pi0dalitz_lep3.txt   PLUTO pi0 Dalitz decays of the lepton set (3 converted files, vae25)
#   docs/latex/tables_all/wp8_qed.txt, wp8_toy_angle.txt, wp8_pty.txt   (wp8_qed.py, wp8_toy_angle.py, wp8_pty.py)
#   /hdd3/.../all_report/wp7/syst.txt              WP7 systematics (Table S1)
#   /hdd3/.../all_report/wp5/results.txt, wp6/results.txt   WP5 fits (cross-checks), WP6 side bands (MAIN)
#   simulation-study table FAT_simwt/plots/eta_sim/tables/closure_ratio_bins.txt (lep/std closure of pi0 Dalitz pairs)
# Typed external numbers (with their source): PDG branching ratios (as in sim_eta_physics.h / A2 note), SMASH sigma, C_sec and
#   f_eta (WP2 + AMENDMENT 1), the data-report revisions (report_exp: P_main model rms 420 for pre + vertex = S0, E_step stat x 1.2),
#   literature values (references in LIT below).
# Outputs: docs/latex/tables_all/wp8_{tableP1,eps2,S12,tableP2}.{txt,tex}, figs_all/wp8_P1.{pdf,png}, wp8_eps2_fig.{pdf,png},
#          /hdd3/.../all_report/wp8/results.txt (machine-readable)
# Run: python3 -I plots/all_report/wp8_physics.py
# ============================================================================
import math
import os
import re

import numpy as np

W8 = "/hdd3/przygoda/HADES/FAT_out/all_report/wp8"
AR = "/hdd3/przygoda/HADES/FAT_out/all_report"
TAB = "/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all"
FIG = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_all"
OUTW = "/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ6/shortfig_r/run/wp8"   # _plain copy: every table / results file goes here
os.makedirs(OUTW, exist_ok=True)
CLOS = "/home/przygoda/HADES/FAT/FAT_simwt/plots/eta_sim/tables/closure_ratio_bins.txt"

# ---- typed constants (source in the comment) ---------------------------------------------------------------------
SIG_SMASH = 1.0665e6          # nb, SMASH sigma(pp -> eta X) (200 converted files; A2 note); 40 std files here: 0.023629 x 45.239 mb
L = 6079.2                    # nb^-1 (D7)
YF = 64.0 / 63.0
BR = dict(ppp0=0.2292, p0D=0.01174, ppg=0.0422, ppee=2.68e-4)   # PDG values of the simulation (sim_eta_physics.h, A2 note)
BR_PPEE_ERR = 0.11e-4
CNORM, CNORM_UP, CNORM_DN = 1.0110 * 0.97, 0.03, 0.04            # C_sec x (1 - f_eta) (AMENDMENT 1): 0.98 +0.03/-0.04
PMAIN_S, PMAIN_STAT, PMAIN_MODEL = 4387.8, 306.6, 420.0          # report_exp review: model = rms of the 144 edge-safe fits (pre + vertex)
EDGE_STAT_SCALE = 1.2                                            # report_exp review: E_step fit errors x 1.2
LIT = [
    ("excl", "pp -> pp eta, HADES, T = 3.5 GeV (sqrt(s) = 3.16 GeV), preliminary", 0.211, 0.047,
     "K. Teilab (HADES), Int. J. Mod. Phys. A 26 (2011) 694, arXiv:1009.3442"),
    ("excl", "pp -> pp eta, HADES, T = 3.5 GeV, PhD thesis", 0.1369, math.hypot(0.0009, 0.0101),
     "K. Teilab, PhD thesis, Goethe University Frankfurt (2011)"),
    ("model", "pp -> pp eta at sqrt(s) = 3.46 GeV (this energy), models 3 / 2 / 1", (0.0808, 0.1162, 0.1369), None,
     "P. Lebiedowicz, arXiv:2501.14952 [hep-ph] (2025), Table 1"),
]


def rd(path):
    with open(path) as f:
        return f.read().splitlines()


# ---- parsers --------------------------------------------------------------------------------------------------------
def parse_eff():
    eff = {}
    for l in rd(os.path.join(W8, "wp8_sim.txt")):
        if not l.startswith("KEY eff"):
            continue
        a = l.split()
        sel, w, comp, win = a[2], a[3], a[4], a[5]
        t = a[6:]
        for i in range(0, len(t), 3):
            eff[(sel, w, comp, win, t[i])] = (float(t[i + 1]), float(t[i + 2]))
    return eff


def parse_tag():
    tag = {}
    for l in rd(os.path.join(W8, "wp8_sim.txt")):
        if not l.startswith("KEY tag"):
            continue
        a = l.split()
        g, k, name = int(a[2]), int(a[3]), a[4]
        d = dict(data=float(a[6]), edata=float(a[7]), sim=float(a[13]), C2=float(a[17]), eC2=float(a[18]), C3=float(a[20]), eC3=float(a[21]))
        tag[(g, name)] = d
    return tag


def parse_syst_raw(path):
    q = {}
    for l in rd(path):
        if l.startswith("q "):
            a = l.split()
            q[(a[1], a[2])] = dict(v=float(a[3]), stat=float(a[4]))
    return q


def parse_wp7():
    s = {}
    for l in rd(os.path.join(AR, "wp7/syst.txt")):
        if l.startswith("syst "):
            a = l.split()
            v = [float(x) for x in a[2:]]
            s[a[1]] = dict(v=v[0], stat=v[1], fam=v[2:14], syst=v[14], noR1=v[15], quad=v[16], statB=v[17], mc=v[18], nup=v[19], ndn=v[20])
    return s


def parse_kv(path, prefix):
    out = {}
    for l in rd(path):
        if l.startswith(prefix):
            a = l.split()
            out[a[1] if len(a) > 1 else ""] = l
    return out


def wp5_derived(block, name):
    for l in rd(os.path.join(AR, "wp5/results.txt")):
        if l.startswith("derived %s %s " % (block, name)) or re.match(r"derived %s %s\s" % (block, re.escape(name)), l):
            a = l.split()
            d = dict(v=float(a[3]))
            for i in range(4, len(a) - 1, 2):
                try:
                    d[a[i]] = float(a[i + 1])
                except ValueError:
                    pass
            return d
    return None


def wp6_est(prefix):
    for l in rd(os.path.join(AR, "wp6/results.txt")):
        if l.startswith("est " + prefix):
            m = re.search(r"\s([-\d.e+]+) stat ([-\d.e+]+) .* mc ([-\d.e+]+) tot", l)
            return dict(v=float(m.group(1)), stat=float(m.group(2)), mc=float(m.group(3)))
    return None


def wp6_xcheck(name):
    for l in rd(os.path.join(AR, "wp6/results.txt")):
        if l.startswith("xcheck " + name):
            return l.split()
    return None


# ---- small helpers --------------------------------------------------------------------------------------------------
def fmt(x, e, nd=3):
    return ("%%.%df +- %%.%df" % (nd, nd)) % (x, e)


def ratio_nested(S, eS, A, eA):
    """S = subset of A (nested): r = S/A with the independent 'fail' part F = A - S, var F = eA^2 - eS^2"""
    F = A - S
    vF = max(eA * eA - eS * eS, 0.0)
    r = S / A
    v = (F * F * eS * eS + S * S * vF) / A ** 4
    return r, math.sqrt(v)


def main():
    out = []          # machine-readable lines
    eff = parse_eff()
    tag = parse_tag()
    w7 = parse_wp7()
    v8 = parse_syst_raw(os.path.join(W8, "wp8_syst_raw.txt"))

    # =================================================================================================================
    # 1. TABLE P1: efficiencies
    # =================================================================================================================
    sprod = {"C2": SIG_SMASH * BR["ppp0"] * BR["p0D"], "C3": SIG_SMASH * BR["ppp0"], "C4": SIG_SMASH * BR["ppg"],
             "C5": SIG_SMASH * BR["ppee"], "C6": SIG_SMASH * 1.0}
    cls_txt = {"C2": "C2 bump, pi0 Dalitz", "C3": "C3 bump, conversion", "C4": "C4 peak, eta->pi+pi-gamma conv.",
               "C5": "C5 direct eta->pi+pi-e+e-", "C6": "C6 X17 (per unit BR)"}
    cls_win = {"C2": "B", "C3": "B", "C4": "P", "C5": "P", "C6": "X17cM4"}
    L1 = []
    def P1(s=""):
        print(s); L1.append(s)
    P1("# wp8_tableP1.txt -- plots/all_report/wp8_physics.py (WP8, 2026-11-14): Table P1, acceptance x efficiency per class")
    P1("# eps = sigma_vis(sim) / sigma_prod; sigma_vis = sum w / L_sim (wp8_sim.C, mode A, owned rows, S0 unless stated, M4' = 0.99615 M4,")
    P1("# M(pipi)' = 0.9952 M(pipi)); sigma_prod = SMASH sigma(pp -> eta X) = 1.0665 mb x PDG BR (C2: 0.2292 x 0.01174; C3: 0.2292, i.e. per")
    P1("# eta -> pi+pi-pi0 with the conversion probability inside eps; C4: 0.0422; C5: 2.68e-4; C6: 1 = per unit BR, at m_X = 17 MeV).")
    P1("# The simulation study (sim report Tab. s5eps) used N_eta/L_sim = 1.069 mb, so its eps are 0.24 % lower. Errors: event-group MC.")
    P1("# windows: C2, C3 in B (0.4550-0.5225), C4, C5 in P (0.53364-0.55416), C6 in X17c (M4 0.5312-0.5566 AND M(ee) 15.5-24.5 MeV);")
    P1("# M(ee) classes reconstructed: E1 < 10, E2 >= 10, >= 20 MeV; 'true >= 20' = true pair mass of the accepted rows.")
    P1("")
    P1("%-30s %11s | %-20s | %-20s %-20s %-20s %-20s | %-20s %8s" % ("class", "sigma_prod", "eps all M4 (V3)", "eps window (V3)", "E1 < 10 MeV",
                                                                  "E2 >= 10 MeV", ">= 20 MeV", "window, no weights", "V3/off"))
    def e_(sel, w, c, win, m):
        v, e = eff[(sel, w, c, win, m)]
        return v / sprod[c], e / sprod[c]
    for c in ("C2", "C3", "C4", "C5", "C6"):
        win = cls_win[c]
        m_all = "x17" if c == "C6" else "all"
        a = e_("S0", "V3", c, "all", "all")
        w = e_("S0", "V3", c, win, m_all)
        e1 = e_("S0", "V3", c, win, "E1") if c != "C6" else (float("nan"), float("nan"))
        e2 = e_("S0", "V3", c, win, "E2") if c != "C6" else (float("nan"), float("nan"))
        g20 = e_("S0", "V3", c, win, "ge20") if c != "C6" else (float("nan"), float("nan"))
        off = e_("S0", "off", c, win, m_all)
        sc = 1e2  # per cent
        P1("%-30s %11.4g | %-20s | %-20s %-20s %-20s %-20s | %-20s %8.4f" % (
            cls_txt[c], sprod[c], "%.4g +- %.2g %%" % (sc * a[0], sc * a[1]), "%.4g +- %.2g %%" % (sc * w[0], sc * w[1]),
            "%.4g +- %.2g %%" % (sc * e1[0], sc * e1[1]), "%.4g +- %.2g %%" % (sc * e2[0], sc * e2[1]), "%.4g +- %.2g %%" % (sc * g20[0], sc * g20[1]),
            "%.4g +- %.2g %%" % (sc * off[0], sc * off[1]), w[0] / off[0]))
        out.append("P1 %s sprod %.6g eff_all %.6g %.3g eff_win %.6g %.3g eff_E1 %.6g eff_E2 %.6g eff_ge20 %.6g eff_off %.6g" % (
            c, sprod[c], a[0], a[1], w[0], w[1], e1[0], e2[0], g20[0], off[0]))
    P1("")
    P1("selection pre (for comparison) and the simulated pion-vertex efficiency S0/pre per class (window, all M(ee), V3):")
    for c in ("C2", "C3", "C4", "C5", "C6"):
        win = cls_win[c]; m_all = "x17" if c == "C6" else "all"
        s0 = eff[("S0", "V3", c, win, m_all)][0]; pr = eff[("pre", "V3", c, win, m_all)][0]
        P1("  %-28s eps(pre, window) = %.4g %%   S0/pre = %.4f" % (cls_txt[c], 100 * pr / sprod[c], s0 / pr))
    # check vs the simulation study (pre, all M4, V3; Tab. s5eps uses 1.069 mb): C5 0.4944 %, C2 0.245 %, C3 0.00212 %, C4 0.00275 %, C6 0.6665 %
    chk = {c: eff[("pre", "V3", c, "all", "all")][0] / (sprod[c] * 1.069 / 1.0665) for c in ("C2", "C3", "C4", "C5", "C6")}
    P1("  check, pre all M4 with the sim-study normalisation 1.069 mb: C2 %.4f %% (0.245), C3 %.5f %% (0.00212), C4 %.5f %% (0.00275), C5 %.4f %% (0.4944), C6 %.4f %% (0.6665)" %
       tuple(100 * chk[c] for c in ("C2", "C3", "C4", "C5", "C6")))
    # C2 migration: accepted C2 in B, reconstructed vs true >= 20 MeV
    rg, tg, al = eff[("S0", "V3", "C2", "B", "ge20")][0], eff[("S0", "V3", "C2", "B", "tge20")][0], eff[("S0", "V3", "C2", "B", "all")][0]
    P1("  C2 in B: share of accepted rows with M(ee) >= 20 MeV: reconstructed %.4f, true %.4f (generator 0.2688)" % (rg / al, tg / al))
    out.append("P1 C2_share_ge20 reco %.5f true %.5f" % (rg / al, tg / al))
    P1("")
    # ---- pion-vertex efficiency, data vs simulation (recomputed with the report_all definitions)
    Yp0, eYp0 = w7["Y_peak"]["v"], w7["Y_peak"]["stat"]
    q7 = {}
    for l in rd(os.path.join(AR, "wp7/syst_raw.txt")):
        if l.startswith("q S11_pre "):
            a = l.split(); q7[a[2]] = (float(a[3]), float(a[4]))
    Ypp, eYpp = q7["Y_peak"]
    rP, erP = ratio_nested(Yp0, eYp0, Ypp, eYpp)
    sP0 = eff[("S0", "V3", "C4", "P", "all")][0] + eff[("S0", "V3", "C5", "P", "all")][0]
    sPp = eff[("pre", "V3", "C4", "P", "all")][0] + eff[("pre", "V3", "C5", "P", "all")][0]
    T0, eT0, Tp, eTp = tag[(0, "all")]["data"], tag[(0, "all")]["edata"], tag[(1, "all")]["data"], tag[(1, "all")]["edata"]
    rT, erT = ratio_nested(T0, eT0, Tp, eTp)
    sT0, sTp = tag[(0, "all")]["sim"], tag[(1, "all")]["sim"]
    tr0, trp = tag[(0, "all")]["C2"] + tag[(0, "all")]["C3"], tag[(1, "all")]["C2"] + tag[(1, "all")]["C3"]
    P1("PION-VERTEX EFFICIENCY S0/pre (both pions in the vertex fit), data vs simulation (recomputed here; WP1 Table A2 quoted 0.917/0.941")
    P1("(peak) and 0.945/0.950 (tagged) from the earlier studies):")
    P1("  eta peak: data Y_peak(S0)/Y_peak(pre) = %.4f +- %.4f (side bands, nested-sample error); sim (C4 + C5 in P) %.4f; data/sim %.4f +- %.4f" %
       (rP, erP, sP0 / sPp, rP / (sP0 / sPp), erP / (sP0 / sPp)))
    P1("  ECAL-tagged eta: data T(S0)/T(pre) = %.4f +- %.4f; sim emulated tag %.4f, truth C2 + C3 %.4f; data/sim %.4f +- %.4f (emulated), %.4f (truth)" %
       (rT, erT, sT0 / sTp, tr0 / trp, rT / (sT0 / sTp), erT / (sT0 / sTp), rT / (tr0 / trp)))
    out.append("P1 pivtx peak data %.5f %.5f sim %.5f tagged data %.5f %.5f sim %.5f truth %.5f" % (rP, erP, sP0 / sPp, rT, erT, sT0 / sTp, tr0 / trp))
    with open(os.path.join(OUTW, "wp8_tableP1.txt"), "w") as f:
        f.write("\n".join(L1) + "\n")
    # tex
    with open(os.path.join(OUTW, "wp8_tableP1.tex"), "w") as f:
        f.write("% wp8_tableP1.tex -- plots/all_report/wp8_physics.py (WP8)\n")
        f.write("\\begin{table}[H]\n\\centering\n\\caption{Acceptance $\\times$ efficiency $\\varepsilon=\\sigma_{\\mathrm{vis}}/\\sigma_{\\mathrm{prod}}$ of the "
                "signal classes in selection S0 ($\\mathrm{isBest}\\geq0$), in per cent. $\\sigma_{\\mathrm{prod}}$ = SMASH $\\sigma(pp\\to\\eta X)=1.0665$\\,mb "
                "$\\times$ PDG BR (C3: per $\\eta\\to\\pipi\\piz$, conversion probability inside $\\varepsilon$; C6: per unit BR, $m_X=17$\\,\\MeVcc). Windows: "
                "C2, C3 in B; C4, C5 in P; C6 in X17c ($M_4$ 0.5312--0.5566, \\Mee 15.5--24.5\\,\\MeVcc). Physics weights V3 (code 108) and none; "
                "\\Mee classes reconstructed. Simulation errors: event groups. Mass transfer $M_4'=0.99615\\,M_4$.}\n\\label{tab:P1}\n\\scriptsize\n"
                "\\setlength{\\tabcolsep}{2.5pt}\n\\begin{tabular}{@{}l r r r r r r r@{}}\n\\toprule\n"
                "class & $\\sigma_{\\mathrm{prod}}$ [nb] & all $M_4$ & window & $E1$ ($<10$) & $E2$ ($\\geq10$) & $\\geq20$\\,\\MeVcc & window, no wts \\\\\n\\midrule\n")
        for c in ("C2", "C3", "C4", "C5", "C6"):
            win = cls_win[c]; m_all = "x17" if c == "C6" else "all"
            a = e_("S0", "V3", c, "all", "all"); w = e_("S0", "V3", c, win, m_all); off = e_("S0", "off", c, win, m_all)
            if c != "C6":
                e1 = e_("S0", "V3", c, win, "E1"); e2 = e_("S0", "V3", c, win, "E2"); g20 = e_("S0", "V3", c, win, "ge20")
                mid = "%.4g & %.4g & %.4g" % (100 * e1[0], 100 * e2[0], 100 * g20[0])
            else:
                mid = "-- & -- & --"
            texlab = {"C2": "C2 bump, \\piz Dalitz", "C3": "C3 bump, conversion", "C4": "C4 peak, $\\eta\\to\\pipi\\gamma$ conv.",
                      "C5": "C5 direct $\\eta\\to\\ppee$", "C6": "C6 X17 (per unit BR)"}[c]
            f.write("%s & %.0f & %.4g (%.1f) & %.4g (%.1f) & %s & %.4g \\\\\n" % (texlab,
                    sprod[c], 100 * a[0], 100 * a[1] / a[0] if a[0] > 0 else 0, 100 * w[0], 100 * w[1] / w[0] if w[0] > 0 else 0, mid, 100 * off[0]))
        f.write("\\bottomrule\n\\end{tabular}\n\\par\\smallskip\\footnotesize Numbers in brackets: relative MC error in per cent (window and all-$M_4$ columns); "
                "pion-vertex errors in units of the last digit (nested-sample error).\\\\\n"
                "Pion-vertex efficiency S0/pre: $\\eta$ peak data %.3f(%.0f), simulation %.3f; ECAL-tagged $\\eta$ data %.3f(%.0f), simulation %.3f (emulated tag), %.3f (truth).\n\\end{table}\n" %
                (rP, 1000 * erP, sP0 / sPp, rT, 1000 * erT, sT0 / sTp, tr0 / trp))

    # =================================================================================================================
    # 2. eps2 QED CHECK
    # =================================================================================================================
    L2 = []
    def P2(s=""):
        print(s); L2.append(s)
    qed = {}
    for l in rd(os.path.join(TAB, "wp8_qed.txt")):
        if l.startswith("qed ") and " | " in l:
            name, rest = l[4:].split(" | ")
            qed[name] = {kv.rsplit("=", 1)[0]: float(kv.rsplit("=", 1)[1]) for kv in rest.split()}
    toy = {}
    for l in rd(os.path.join(TAB, "wp8_toy_angle.txt")):
        if l.startswith("toy ") and " | " in l:
            name, rest = l[4:].split(" | ")
            toy[name] = {kv.rsplit("=", 1)[0]: float(kv.rsplit("=", 1)[1]) for kv in rest.split()}
    # PLUTO check
    pl = rd(os.path.join(W8, "remote/pluto_pi0dalitz_lep3.txt"))
    plc = {0: [0] * 20, 1: [0] * 20, 2: [0] * 20}; plq = [0] * 70; npl = 0
    for l in pl:
        a = l.split()
        if a[0] == "cos":
            plc[int(a[1])][int(a[2])] = int(a[3])
        elif a[0] == "q":
            plq[int(a[1])] = int(a[2])
        elif l.startswith("# events"):
            npl = int(a[4])
    cbin = np.array([-0.95 + 0.1 * i for i in range(20)])
    c2pl = {}
    for k in range(3):
        h = np.array(plc[k], float)
        c2 = float(np.sum(h * (cbin ** 2 + 1.0 / 1200.0)) / h.sum())   # bin-centre <c^2> + bin-width correction (0.1^2/12)
        c2pl[k] = (c2, math.sqrt(max((np.sum(h * cbin ** 4) / h.sum() - c2 ** 2), 0) / h.sum()))
    qall = sum(plq); q20 = sum(plq[10:]); q10 = sum(plq[5:])
    fpl20 = q20 / qall; efpl20 = math.sqrt(fpl20 * (1 - fpl20) / qall)
    # closure of the simulation study
    clos = []
    for l in rd(CLOS):
        a = l.split()
        if len(a) == 8 and a[0][0].isdigit():
            clos.append((float(a[0]), float(a[1]), float(a[4]), float(a[5])))
        if l.startswith("## pi0 Dalitz pair, pre+OA>2"):
            break
    clos = clos[:15]
    def clos_cls(lo, hi):
        rows = [(n, r) for (a, b, n, r) in clos if a >= lo - 1e-9 and b <= hi + 1e-9]
        return sum(n * r for n, r in rows) / sum(n for n, r in rows)
    cE1, cE2, c20, call = clos_cls(0, 0.010), clos_cls(0.010, 1), clos_cls(0.020, 1), clos_cls(0, 1)
    hg = qed["HGeant_histogram_9995_(C2_generator)"]; kwff = qed["Kroll-Wada_x_(1_+_a_x)^2,_a_=_0.0332_(PDG)"]; kw0 = qed["QED_Kroll-Wada,_no_form_factor"]
    rq20 = fpl20 / hg[">=20"]; rq10 = (q10 / qall) / hg[">=10(E2)"]; rqE1 = (1 - q10 / qall) / (1 - hg[">=10(E2)"])
    P2("# wp8_eps2.txt -- plots/all_report/wp8_physics.py (WP8, 2026-11-14): eps2 QED check of class C2 (bump pi0 Dalitz) at M(ee) >= 20 MeV")
    P2("")
    P2("A. GENERATOR LEVEL (fractions of all produced pi0 -> e+e-gamma decays; wp8_qed.py, wp8/remote/pluto_pi0dalitz_lep3.txt)")
    P2("   %-52s %9s %9s %9s %9s %9s" % ("", "<8", "8-20", "20-30", ">=30", ">=20"))
    for nm, d in (("QED Kroll-Wada", kw0), ("Kroll-Wada x PDG form factor (a = 0.0332)", kwff), ("HGeant (C2 generator, source uhbook.F/gdalet.F)", hg)):
        P2("   %-52s %9.4f %9.4f %9.4f %9.4f %9.4f" % (nm, d["<8"], d["8-20"], d["20-30"], d[">=30"], d[">=20"]))
    P2("   %-52s %9.4f %9.4f %9.4f %9.4f %9.4f   (+- %.4f)" % ("PLUTO, lepton set (%d decays, 3 converted files)" % qall, sum(plq[:4]) / qall, sum(plq[4:10]) / qall,
       sum(plq[10:15]) / qall, sum(plq[15:]) / qall, fpl20, efpl20))
    P2("   HGeant / (KW x PDG FF) at >= 20 MeV: %.4f -> the C2 pair-mass spectrum is QED (difference %.2f %%)." % (hg[">=20"] / kwff[">=20"], 100 * (hg[">=20"] / kwff[">=20"] - 1)))
    P2("   PLUTO / HGeant at >= 20 MeV: %.4f (PLUTO's form factor is harder; PLUTO feeds G_pi0D, not C2)" % rq20)
    P2("")
    P2("B. DECAY ANGLE of gamma* -> e+e- (lepton angle theta* to the gamma* direction, gamma* rest frame)")
    P2("   QED: dN/dcos ~ 1 + cos^2 + (1 - beta^2) sin^2  -> <cos^2> = 0.400 for beta -> 1; isotropic: <cos^2> = 0.333")
    P2("   HGeant gdalet.F: YINV uniform in [-beta, beta] = cos(theta*) uniform -> ISOTROPIC (source)")
    P2("   PLUTO (lepton set, measured on the converted events): <cos^2> = %.4f +- %.4f (M(ee) < 10), %.4f +- %.4f (10-20), %.4f +- %.4f (>= 20 MeV) -> QED" %
       (c2pl[0][0], c2pl[0][1], c2pl[1][0], c2pl[1][1], c2pl[2][0], c2pl[2][1]))
    P2("   => the lep/std closure of pi0 Dalitz pairs (0.891 +- 0.003, called 'unexplained' in the simulation study) is the decay-angle difference")
    P2("")
    P2("C. EFFECT ON THE C2 EFFICIENCY: eps(QED angle) / eps(HGeant isotropic)")
    P2("   %-56s %8s %8s %8s %8s" % ("", "all", "E1 <10", "E2 >=10", ">=20"))
    for nm in toy:
        t = toy[nm]
        P2("   %-56s %8.4f %8.4f %8.4f %8.4f" % ("toy, " + nm.replace("_", " "), t["all"], t["<10"], t[">=10"], t[">=20"]))
    P2("   %-56s %8.4f %8.4f %8.4f %8.4f" % ("lep/std closure, full simulation (prompt SMASH pi0, pre)", call, cE1, cE2, c20))
    P2("   %-56s %8.4f %8.4f %8.4f %8.4f" % ("closure / (PLUTO/HGeant q spectrum)  [used: ang_clos]", call, cE1 / rqE1, cE2 / rq10, c20 / rq20))
    P2("   toy nominal [used: ang_toy]: all 0.8800, E1 0.8886, E2 0.8675, >= 20 0.8864")
    out.append("eps2 closure all %.5f E1 %.5f E2 %.5f ge20 %.5f qcorr E1 %.5f E2 %.5f ge20 %.5f" % (call, cE1, cE2, c20, cE1 / rqE1, cE2 / rq10, c20 / rq20))
    out.append("eps2 pluto cos2 %.5f %.5f %.5f frac20 %.5f %.5f" % (c2pl[0][0], c2pl[1][0], c2pl[2][0], fpl20, efpl20))
    P2("")
    # tagged shares
    P2("D. ECAL-TAGGED eta, SHARES PER M(ee) CLASS (S0, all M4; frozen m5c fit; share = T(class) / sum of the four disjoint classes)")
    classes = ["<8", "8-20", "20-30", "30-100"]
    def shares(get, geterr):
        T = np.array([get(c) for c in classes]); E = np.array([geterr(c) for c in classes])
        S = T.sum()
        sh = T / S
        es = np.sqrt(((S - T) ** 2 * E ** 2 + T ** 2 * (np.sum(E ** 2) - E ** 2)) / S ** 4)
        return sh, es
    def agg(sh, es, idx):
        s = sum(sh[i] for i in idx)
        return s, math.sqrt(sum(es[i] ** 2 for i in idx))   # approximation (classes nearly independent)
    shD, esD = shares(lambda c: tag[(0, c)]["data"], lambda c: tag[(0, c)]["edata"])
    # simulated: MC error of the fit proxied by the relative MC error of the truth C2 + C3 in the class
    def truth_rel(c):
        d = tag[(0, c)]; v = d["C2"] + d["C3"]
        return math.hypot(d["eC2"], d["eC3"]) / v if v > 0 else 0
    shS, esS = shares(lambda c: tag[(0, c)]["sim"], lambda c: tag[(0, c)]["sim"] * truth_rel(c))
    shT, esT = shares(lambda c: tag[(0, c)]["C2"] + tag[(0, c)]["C3"], lambda c: math.hypot(tag[(0, c)]["eC2"], tag[(0, c)]["eC3"]))
    shC2 = np.array([tag[(0, c)]["C2"] for c in classes]) / np.array([tag[(0, c)]["C2"] + tag[(0, c)]["C3"] for c in classes])
    P2("   %-10s %18s %18s %18s %10s %16s %16s" % ("M(ee)", "data", "sim (fit)", "sim (truth C2+C3)", "C2 frac", "data/sim(fit)", "data/sim(truth)"))
    rows = [(classes[i], shD[i], esD[i], shS[i], esS[i], shT[i], esT[i], shC2[i]) for i in range(4)]
    for idx, nm in (([2, 3], ">=20"), ([1, 2, 3], ">=8")):
        a = agg(shD, esD, idx); b = agg(shS, esS, idx); c = agg(shT, esT, idx)
        rows.append((nm, a[0], a[1], b[0], b[1], c[0], c[1], float("nan")))
    for r in rows:
        rf = r[1] / r[3]; erf = rf * math.hypot(r[2] / r[1], r[4] / r[3])
        rt = r[1] / r[5]; ert = rt * math.hypot(r[2] / r[1], r[6] / r[5])
        P2("   %-10s %8.4f +- %6.4f %8.4f +- %6.4f %8.4f +- %6.4f %10.3f %8.3f +- %5.3f %8.3f +- %5.3f" % (r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], rf, erf, rt, ert))
        out.append("eps2 tagshare %s data %.5f %.5f simfit %.5f %.5f truth %.5f %.5f ratio_fit %.4f %.4f ratio_truth %.4f %.4f" % (r[0], r[1], r[2], r[3], r[4], r[5], r[6], rf, erf, rt, ert))
    # >= 20 among >= 10
    d20 = tag[(0, ">=20")]["data"] / tag[(0, ">=10")]["data"]
    s20 = tag[(0, ">=20")]["sim"] / tag[(0, ">=10")]["sim"]
    t20 = (tag[(0, ">=20")]["C2"] + tag[(0, ">=20")]["C3"]) / (tag[(0, ">=10")]["C2"] + tag[(0, ">=10")]["C3"])
    ed20 = d20 * math.sqrt((tag[(0, ">=20")]["edata"] / tag[(0, ">=20")]["data"]) ** 2 * (1 - d20) ** 2 + (math.hypot(tag[(0, ">=10")]["edata"], 0) / tag[(0, ">=10")]["data"]) ** 2 * d20 ** 2 * 0)
    ed20 = math.sqrt(max(d20 * (1 - d20), 0) * 0 + (tag[(0, ">=20")]["edata"] / tag[(0, ">=10")]["data"]) ** 2 * (1 - d20) ** 2 +
                     (math.sqrt(max(tag[(0, ">=10")]["edata"] ** 2 - tag[(0, ">=20")]["edata"] ** 2, 0)) / tag[(0, ">=10")]["data"]) ** 2 * d20 ** 2)
    P2("   share of >= 20 MeV among the tagged eta with M(ee) >= 10 MeV: data %.4f +- %.4f, sim fit %.4f, truth %.4f (C2 = %.3f of truth)" %
       (d20, ed20, s20, t20, tag[(0, ">=20")]["C2"] / (tag[(0, ">=20")]["C2"] + tag[(0, ">=20")]["C3"])))
    P2("   reproduction of the data scan (A1 note, bump start selection): share 20-30 %.4f (0.077), > 30 %.4f (0.144)" %
       tuple(tag[(2, c)]["data"] / sum(tag[(2, x)]["data"] for x in classes) for c in ("20-30", "30-100")))
    out.append("eps2 share20of10 data %.5f %.5f simfit %.5f truth %.5f" % (d20, ed20, s20, t20))
    P2("")
    # refit effects
    P2("E. EFFECT ON THE REPORT QUANTITIES (wp8_syst.C: the WP7 chain with C2 x the angle factors; relative change to the nominal)")
    qn = ["a_eta", "Y_bump", "Y_bump_E2", "a_B", "DR", "DR_F", "BR_G", "BR_A", "rho_d", "kappa", "Y5", "Y4", "Y_peak", "a_P", "sigma_eta", "sigma_eta_F"]
    P2("   %-12s %12s %10s %10s" % ("quantity", "nominal", "ang_toy", "ang_clos"))
    angsh = {}
    for q in qn:
        n = v8[("nom", q)]["v"]; a = v8[("ang_toy", q)]["v"]; b = v8[("ang_clos", q)]["v"]
        angsh[q] = (a / n - 1, b / n - 1)
        P2("   %-12s %12.5g %+9.2f%% %+9.2f%%" % (q, n, 100 * (a / n - 1), 100 * (b / n - 1)))
        out.append("eps2 shift %s nom %.8g toy %.6f clos %.6f" % (q, n, a / n - 1, b / n - 1))
    with open(os.path.join(OUTW, "wp8_eps2.txt"), "w") as f:
        f.write("\n".join(L2) + "\n")
    with open(os.path.join(OUTW, "wp8_eps2.tex"), "w") as f:
        f.write("% wp8_eps2.tex -- plots/all_report/wp8_physics.py (WP8)\n")
        f.write("\\begin{table}[H]\n\\centering\n\\caption{The $\\varepsilon_2$ QED check of class C2 ($\\eta\\to\\pipi\\piz$, $\\piz\\to\\gamma e^+e^-$). "
                "(a) Generator: fraction of all produced Dalitz decays with $\\Mee\\geq20$\\,\\MeVcc. (b) Decay angle of $\\gamma^*\\to e^+e^-$: "
                "$\\langle\\cos^2\\theta^*\\rangle$ (QED 0.400, isotropic 0.333). (c) C2 efficiency with the QED angle relative to HGeant (isotropic): toy and "
                "full-simulation lep/std closure. (d) ECAL-tagged $\\eta$ in S0: share with $\\Mee\\geq20$\\,\\MeVcc, data and simulation (same frozen fit; "
                "simulated tag emulated, UNVERIFIED).}\n\\label{tab:eps2}\n\\footnotesize\n\\begin{tabular}{@{}l l@{}}\n\\toprule\n")
        f.write("(a) Kroll--Wada / KW $\\times$ PDG form factor / HGeant (C2) / PLUTO (lep set) & %.4f / %.4f / %.4f / %.4f(13) \\\\\n" %
                (kw0[">=20"], kwff[">=20"], hg[">=20"], fpl20))
        f.write("(b) HGeant \\texttt{gdalet.F} (C2) / PLUTO measured ($\\Mee<10$, $\\geq20$) & 0.333 (source) / %.3f, %.3f \\\\\n" % (c2pl[0][0], c2pl[2][0]))
        f.write("(c) $\\varepsilon_{\\mathrm{QED}}/\\varepsilon_{\\mathrm{iso}}$: all / E1 / E2 / $\\geq20$, toy & 0.880 / 0.889 / 0.868 / 0.886 \\\\\n")
        f.write("(c) the same, lep/std closure ($q$-spectrum corrected) & %.3f / %.3f / %.3f / %.3f \\\\\n" % (call, cE1 / rqE1, cE2 / rq10, c20 / rq20))
        a_ = agg(shD, esD, [2, 3]); b_ = agg(shS, esS, [2, 3]); c_ = agg(shT, esT, [2, 3])
        f.write("(d) tagged share $\\Mee\\geq20$: data / sim (fit) / sim (truth C2+C3) & %.3f(%.0f) / %.3f(%.0f) / %.3f(%.0f) \\\\\n" %
                (a_[0], 1000 * a_[1], b_[0], 1000 * b_[1], c_[0], 1000 * c_[1]))
        f.write("(e) effect of the QED angle: $a_\\eta$, BR$_G$, BR$_A$, DR, $\\sigma_\\eta$ (side bands) & %+.0f, %+.0f, %+.0f, %+.0f, %+.1f\\,\\%% (toy) \\\\\n" %
                tuple(100 * angsh[q][0] for q in ("a_eta", "BR_G", "BR_A", "DR", "sigma_eta")))
        f.write("\\bottomrule\n\\end{tabular}\n\\end{table}\n")

    # =================================================================================================================
    # 3. S12 UPDATE (eps2 decay angle, (pT4, y4) reweighting) and new totals
    # =================================================================================================================
    pty = {}
    for l in rd(os.path.join(TAB, "wp8_pty.txt")):
        if l.startswith("pty "):
            a = l.split(); pty["3x3"] = {a[i]: float(a[i + 1]) for i in range(1, len(a) - 1, 2)}
        if l.startswith("pty1d "):
            a = l.split(); pty[a[1]] = {a[i]: float(a[i + 1]) for i in range(2, len(a) - 1, 2)}
    noise_abs = pty["3x3"]["noise"]          # noise of the produced-weighted absolute scale
    noise_ratio = pty["3x3"]["noiseDR"]      # noise of the peak/bump ratios
    L3 = []
    def P3(s=""):
        print(s); L3.append(s)
    P3("# wp8_S12.txt -- plots/all_report/wp8_physics.py (WP8, 2026-11-14): S12 (efficiency model) completed by WP8; per cent of the value")
    P3("# items: pion vertex (WP7: 2.5 % absolute / 2.0 % ratios); eps2 = C2 decay angle (max of the toy and the closure refit, wp8_syst.C;")
    P3("# sign known: see wp8_eps2.txt E); pT4,y4 = 3 x 3 reweighting (wp8_pty.py + wp8_syst.C refit), Barlow-subtracted with its toy noise")
    P3("# (%.1f %% absolute scales, %.1f %% ratios). S12 = max over the items (spec rule); 'quad' = items in quadrature. New total = old total with S12 replaced." %
       (100 * noise_abs, 100 * noise_ratio))
    P3("# NOT in the total (one decision for the analyst): the one-dimensional (pT4 only / y4 only) reweightings, which are significant:")
    P3("#   pT only %+.1f %% +- %.1f %% (Barlow %.1f %%), y only %+.1f %% +- %.1f %% (Barlow %.1f %%) on the absolute eta scales (sigma_eta, a_P, a_eta, a_B)." %
       (100 * pty["pT"]["rel"], 100 * pty["pT"]["noise"], 100 * pty["pT"]["counted"], 100 * pty["y"]["rel"], 100 * pty["y"]["noise"], 100 * pty["y"]["counted"]))
    P3("")
    P3("%-12s %9s %9s %9s %9s %9s %9s %9s %10s %10s %9s" % ("quantity", "S12 old", "pivtx", "eps2", "pT4,y4", "S12 new", "S12 quad", "syst old", "syst new",
                                                           "syst quad", "eps2 sign"))
    quant = ["Y_bump", "Y_peak", "Y5", "Y4", "a_eta", "BR_G", "BR_A", "kappa", "DR", "sigma_eta", "Y_bump_E2", "a_P", "a_B", "DR_F", "sigma_eta_F", "rho_d", "g"]
    absq = {"a_eta", "sigma_eta", "a_P", "a_B", "sigma_eta_F"}
    ratioq = {"BR_G", "BR_A", "kappa", "DR", "DR_F", "rho_d"}
    newsyst = {}
    for q in quant:
        s = w7[q]; val = s["v"]
        old12 = 100 * s["fam"][11] / abs(val)
        piv = old12
        if ("nom", q) in v8:
            n = v8[("nom", q)]["v"]
            a1 = v8[("ang_toy", q)]["v"] / n - 1; a2 = v8[("ang_clos", q)]["v"] / n - 1
            e2 = 100 * max(abs(a1), abs(a2))
            sign = "+" if (a1 + a2) > 0 else "-"
            p = v8[("pty3x3", q)]["v"] / n - 1
            nz = noise_abs if q in absq else noise_ratio
            pt = 100 * math.sqrt(max(0.0, p * p - nz * nz))
        else:
            e2, pt, sign = 0.0, 0.0, ""
        if q == "Y_peak":   # raw yield: the side-band estimator does not use the C2 efficiency (shift -0.06 %); keep the computed numbers
            pass
        new12 = max(piv, e2, pt)
        quad12 = math.sqrt(piv ** 2 + e2 ** 2 + pt ** 2)
        oldtot = 100 * s["syst"] / abs(val)
        newtot = math.sqrt(max(oldtot ** 2 - old12 ** 2 + new12 ** 2, 0))
        quadtot = math.sqrt(max(oldtot ** 2 - old12 ** 2 + quad12 ** 2, 0))
        newsyst[q] = dict(val=val, old12=old12, piv=piv, eps2=e2, pty=pt, new12=new12, quad12=quad12, oldtot=oldtot, newtot=newtot, quadtot=quadtot, sign=sign)
        P3("%-12s %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %10.2f %10.2f %9s" % (q, old12, piv, e2, pt, new12, quad12, oldtot, newtot, quadtot, sign))
        out.append("S12 %s old %.4f pivtx %.4f eps2 %.4f pty %.4f new %.4f quad %.4f syst_old %.4f syst_new %.4f syst_quad %.4f sign %s" %
                   (q, old12, piv, e2, pt, new12, quad12, oldtot, newtot, quadtot, sign or "0"))
    with open(os.path.join(OUTW, "wp8_S12.txt"), "w") as f:
        f.write("\n".join(L3) + "\n")

    # =================================================================================================================
    # 4. TABLE P2
    # =================================================================================================================
    L4 = []
    def P4(s=""):
        print(s); L4.append(s)
    def row(name, val, stat, syst, mc, nup=None, ndn=None, unit="", status="", extra=""):
        s = "%-46s %12.5g  stat %10.4g  syst %10.4g  MC %10.4g" % (name, val, stat, syst, mc)
        if nup is not None:
            s += "  N +%.4g/-%.4g" % (nup, ndn)
        P4(s + ("  [%s]" % unit if unit else "") + ("  %s" % status if status else "") + ("  %s" % extra if extra else ""))
    sy = lambda q: w7[q]["v"] * newsyst[q]["newtot"] / 100.0
    P4("# wp8_tableP2.txt -- plots/all_report/wp8_physics.py (WP8, 2026-11-14): Table P2, physics quantities")
    P4("# S0, isBest >= 0, L = 6079.2 nb^-1; yields raw rows x 64/63. syst = WP7 Table S1 total with S12 completed by WP8 (wp8_S12.txt, spec rule);")
    P4("# MC = template statistics (100 replicas); N = normalisation (absolute quantities only: L 6.9 %, C_sec (1 - f_eta) +3.1/-4.1 %, PDG BRs);")
    P4("# live time +7.3/-0 % on cross sections, separate (not applied, D17). MAIN = estimators that passed WP4 (peak side bands; B4 marginal,")
    P4("# EXTRAPOLATED); CROSS-CHECK = F-main / F-gold / F-free (WP4 FAIL).")
    P4("")
    P4("1. BR(eta -> pi+pi-e+e-)   [PDG 2.68 +- 0.11 e-4]")
    pdg, epdg = 2.68e-4, 0.11e-4
    br = {}
    for q, nm, st in (("BR_G", "BR_G = 2.68e-4 alpha_5/alpha_2 (F-gold, M(ee) >= 20)", "CROSS-CHECK (primary BR method of D15)"),
                      ("BR_A", "BR_A = 2.68e-4 rho_d (F-main)", "CROSS-CHECK")):
        s = w7[q]
        tot = math.sqrt(s["stat"] ** 2 + sy(q) ** 2 + s["mc"] ** 2)
        pull = (s["v"] - pdg) / math.hypot(tot, epdg)
        br[q] = (s["v"], s["stat"], sy(q), s["mc"], s["nup"])
        row(nm, s["v"], s["stat"], sy(q), s["mc"], s["nup"], s["ndn"], status=st, extra="pull vs PDG %+.2f (total)" % pull)
        out.append("P2 %s %.6g stat %.4g syst %.4g mc %.4g N %.4g pull %.3f" % (q, s["v"], s["stat"], sy(q), s["mc"], s["nup"], pull))
    bc = wp5_derived("FFREE_SPEC", "BR_C"); r34s = wp5_derived("FFREE_SPEC", "R_34"); r34 = wp5_derived("FFREE", "R_34")
    P4("%-46s %12.5g  stat %10.4g  MC(mc68) %10.4g  [E1 + E2 fit, k_E1 = 7.6; WITHDRAWN in the nominal E2 configuration]" %
       ("BR_C = 2.68e-4 alpha_5 (alpha_3/alpha_2)/alpha_4 (F-free)", bc["v"], bc["stat"], bc["mc68"]))
    for q, lab in (("BR_G", "BR_G"), ("BR_A", "BR_A")):
        n = v8[("nom", q)]["v"]
        P4("   %s with the QED decay angle of C2 (variant, not applied): %.3g (toy), %.3g (closure)  [%+.1f %%, %+.1f %%]" %
           (lab, v8[("ang_toy", q)]["v"], v8[("ang_clos", q)]["v"], 100 * (v8[("ang_toy", q)]["v"] / n - 1), 100 * (v8[("ang_clos", q)]["v"] / n - 1)))
    P4("")
    P4("2. R_34 = alpha_4/alpha_3 (conversion closure, 1 expected)")
    P4("   F-free E2 (nominal configuration): %.2f +- %.2f (stat) +- %.2f (MC) -> undetermined" % (r34["v"], r34["stat"], r34["mc"]))
    P4("   F-free E1 + E2 (as specified, k_E1 = 7.6): %.3f +- %.3f (stat) +- %.3f (MC); consistent with 1 (fall-back 'R_34 off by > 2 sigma' not fired)" %
       (r34s["v"], r34s["stat"], r34s["mc"]))
    P4("")
    P4("3. kappa (twice): conversion scale of the signal (kappa) and of the non-eta background (kappa_B)")
    s = w7["kappa"]
    row("kappa (F-main, constraint 1 +- 0.2)", s["v"], s["stat"], sy("kappa"), s["mc"], status="CROSS-CHECK, constraint dominated")
    kf = wp5_derived("FMAIN_kfree", "kappa"); kbf = wp5_derived("FMAIN_kfree", "kappa_B"); kb = wp5_derived("FMAIN", "kappa_B")
    ks = wp5_derived("FMAIN_SPEC", "kappa"); kbs = wp5_derived("FMAIN_SPEC", "kappa_B")
    P4("   kappa free (no constraint): %.3f +- %.3f (stat) +- %.3f (MC); kappa_B = %.2f +- %.2f +- %.2f" % (kf["v"], kf["stat"], kf["mc"], kbf["v"], kbf["stat"], kbf["mc"]))
    P4("   kappa_B (F-main nominal): %.2f +- %.2f (stat) +- %.2f (MC)" % (kb["v"], kb["stat"], kb["mc"]))
    P4("   E1 + E2 fit (k_E1 = 7.6): kappa = %.3f +- %.3f +- %.3f, kappa_B = %.3f +- %.3f +- %.3f (difference %.2f; E1 carries an unmodelled component)" %
       (ks["v"], ks["stat"], ks["mc"], kbs["v"], kbs["stat"], kbs["mc"], kbs["v"] - ks["v"]))
    a2s, a3s, a4s = wp5_derived("FFREE_SPEC", "alpha_2"), wp5_derived("FFREE_SPEC", "alpha_3"), wp5_derived("FFREE_SPEC", "alpha_4")
    P4("   F-free E1 + E2: alpha_3/alpha_2 = %.3f (bump conversion), alpha_4/alpha_2 = %.3f (peak conversion)" % (a3s["v"] / a2s["v"], a4s["v"] / a2s["v"]))
    P4("")
    P4("4. DR = (Y_peak/Y_bump)_data / (Y_peak/Y_bump)_SMASH")
    s = w7["DR"]
    row("DR, side bands / B4 (bump EXTRAPOLATED)", s["v"], s["stat"], sy("DR"), s["mc"], s["nup"], s["ndn"], status="MAIN (B4 marginal)")
    s = w7["DR_F"]
    row("DR, F-main", s["v"], s["stat"], sy("DR_F"), s["mc"], s["nup"], s["ndn"], status="CROSS-CHECK")
    P4("   with the QED decay angle of C2 (variant): DR %.3f / %.3f, DR_F %.3f / %.3f (toy / closure)" %
       (v8[("ang_toy", "DR")]["v"], v8[("ang_clos", "DR")]["v"], v8[("ang_toy", "DR_F")]["v"], v8[("ang_clos", "DR_F")]["v"]))
    P4("")
    P4("5. sigma(pp -> eta X) = a x 1.0665 mb x C_sec (1 - f_eta), C_sec (1 - f_eta) = 0.981 +0.03/-0.04  [SMASH: 1.0665 mb]")
    s = w7["sigma_eta"]
    sig, sst, ssy, smc, snu, snd = s["v"], s["stat"], sy("sigma_eta"), s["mc"], s["nup"], s["ndn"]
    row("sigma(pp -> eta X), peak side bands (a_P)", sig, sst, ssy, smc, snu, snd, unit="mb", status="MAIN", extra="live time +%.3f/-0" % (0.073 * sig))
    totsig = math.sqrt(sst ** 2 + ssy ** 2 + smc ** 2 + ((snu + snd) / 2) ** 2)
    P4("   = %.3f +- %.3f (stat) +- %.3f (syst) +- %.3f (MC) +%.3f/-%.3f (norm.) mb; total %.3f mb (%.1f %%); data/SMASH = %.3f +- %.3f" %
       (sig, sst, ssy, smc, snu, snd, totsig, 100 * totsig / sig, sig / 1.0665, totsig / 1.0665))
    P4("   production-model extrapolation (not in the total, analyst's decision): pT-only reweighting %+.1f %% (Barlow %.1f %%), y-only %+.1f %% (%.1f %%);" %
       (100 * pty["pT"]["rel"], 100 * pty["pT"]["counted"], 100 * pty["y"]["rel"], 100 * pty["y"]["counted"]))
    P4("   the specified 3 x 3 variant: %+.1f %% +- %.1f %% (noise) -> 0 after the Barlow subtraction" % (100 * pty["3x3"]["rel"], 100 * pty["3x3"]["noise"]))
    out.append("P2 sigma_eta %.6f stat %.5f syst %.5f mc %.5f Nup %.5f Ndn %.5f live %.5f total %.5f ptyPT %.4f ptyY %.4f pty3x3 %.4f %.4f" %
               (sig, sst, ssy, smc, snu, snd, 0.073 * sig, totsig, pty["pT"]["rel"], pty["y"]["rel"], pty["3x3"]["rel"], pty["3x3"]["noise"]))
    s = w7["sigma_eta_F"]
    row("sigma(pp -> eta X), F-main a_eta", s["v"], s["stat"], sy("sigma_eta_F"), s["mc"], s["nup"], s["ndn"], unit="mb", status="CROSS-CHECK")
    P4("   F-main with the QED decay angle of C2 (variant): %.3f / %.3f mb (toy / closure)" % (v8[("ang_toy", "sigma_eta_F")]["v"], v8[("ang_clos", "sigma_eta_F")]["v"]))
    fe = wp6_xcheck("f_edge")
    fedge, fst, fmc, fmod = float(fe[2]), float(fe[3]) * EDGE_STAT_SCALE, float(fe[4]), float(fe[5])
    P4("   f_eta(edge) = %.3f +- %.3f (stat, x 1.2 calibrated) +- %.3f (MC) +- %.3f (model) -> sigma = %.3f mb (bump classes; with the QED angle of C2 about +7 %%)" %
       (fedge, fst, fmc, fmod, fedge * 1.0665 * CNORM))
    ft = wp6_xcheck("f_tag")
    P4("   f_eta(tag) = %.3f +- %.3f +- %.3f (emulated simulated tag, UNVERIFIED) -> %.3f mb" % (float(ft[2]), float(ft[3]), float(ft[4]), float(ft[2]) * 1.0665 * CNORM))
    pm = wp6_xcheck("pmain")
    simS = float(pm[9])
    P4("   P_main (S0): S = %.0f +- %.0f (stat) +- %.0f (model, rms of the edge-safe grid, report_exp) ; / sim frozen fit %.0f = %.3f +- %.3f +- %.3f -> %.3f mb" %
       (PMAIN_S, PMAIN_STAT, PMAIN_MODEL, simS, PMAIN_S / simS, PMAIN_STAT / simS, PMAIN_MODEL / simS, PMAIN_S / simS * 1.0665 * CNORM))
    P4("   SMASH sigma check: 40 converted standard-set files give 0.023629 eta per collision = %.4f +- %.4f mb (200 files: 1.0665 mb)" %
       (0.023629 * 45.239, 0.023629 * 45.239 / math.sqrt(94517)))
    P4("")
    P4("6. sigma x BR [nb] (sigma from this analysis; BR as stated)")
    sxb = sig * 1e6 * 2.68e-4
    P4("   sigma(pp -> eta X) x BR_PDG(eta -> pi+pi-e+e-) = %.1f nb (+- %.1f %% from sigma, +- 4.1 %% from the BR)" % (sxb, 100 * totsig / sig))
    s5 = wp5_derived("FMAIN", "Y5_tot"); s4 = wp5_derived("FMAIN", "Y4_tot"); s2 = wp5_derived("FMAIN", "Y2_tot")
    a = w7["a_eta"]["v"]; rho = w7["rho_d"]["v"]; kap = w7["kappa"]["v"]
    sxb5 = a * rho * SIG_SMASH * BR["ppee"] * CNORM
    P4("   sigma x BR(eta -> pi+pi-e+e-) measured (F-main a_eta rho_d, CROSS-CHECK): %.1f nb (stat %.1f %%, syst %.0f %%, MC %.0f %%)" %
       (sxb5, 100 * math.hypot(w7["a_eta"]["stat"] / a, w7["rho_d"]["stat"] / rho), newsyst["BR_A"]["newtot"], 100 * w7["rho_d"]["mc"] / rho))
    sxb3 = w7["a_B"]["v"] * SIG_SMASH * BR["ppp0"] * CNORM
    P4("   sigma x BR(eta -> pi+pi-pi0) from B4 (a_B, EXTRAPOLATED): %.0f nb (syst %.0f %%)" % (sxb3, newsyst["a_B"]["newtot"]))
    sxbg = sig * 1e6 * BR["ppg"]
    P4("   sigma x BR_PDG(eta -> pi+pi-gamma) = %.0f nb; F-main a kappa x SMASH x BR: %.0f nb (kappa constraint dominated)" % (sxbg, a * kap * SIG_SMASH * BR["ppg"] * CNORM))
    out.append("P2 sxBR ppee %.3f ppee_F %.3f ppp0_B4 %.2f ppg %.2f" % (sxb, sxb5, sxb3, sxbg))
    P4("")
    P4("7. VISIBLE CROSS SECTIONS sigma_vis = Y x 64/63 / L (selection S0, no efficiency correction, C_sec (1 - f_eta) NOT applied) [nb]")
    for q, nm, st in (("Y_peak", "peak window P, side bands, all M(ee)", "MAIN"), ("Y_bump", "bump window B, B4, all M(ee)", "MAIN (EXTRAPOLATED)")):
        s = w7[q]
        row("sigma_vis " + nm, s["v"] / L, s["stat"] / L, sy(q) / L, s["mc"] / L, L * 0 + s["v"] / L * 0.069, s["v"] / L * 0.069, unit="nb", status=st)
    e1 = wp6_est("PSB_Ypeak_E1"); e2 = wp6_est("PSB_Ypeak_E2"); be2 = wp6_est("B4_Ybump_E2")
    P4("   peak E1 (M(ee) < 10): %.3f +- %.3f (stat) +- %.3f (MC); E2: %.3f +- %.3f +- %.3f; bump E2 (B4, not extrapolated): %.3f +- %.3f +- %.3f" %
       (e1["v"] / L, e1["stat"] / L, e1["mc"] / L, e2["v"] / L, e2["stat"] / L, e2["mc"] / L, be2["v"] / L, be2["stat"] / L, be2["mc"] / L))
    for nm in ("Y2_B", "Y3_B", "Y4_P", "Y5_P"):
        d = wp5_derived("FMAIN", nm)
        P4("   F-main %s / L = %.3f +- %.3f (stat) +- %.3f (MC) nb (CROSS-CHECK)" % (nm, d["v"] / L, d["stat"] / L, d["mc"] / L))
    P4("   SMASH visible (S0, V3): C4 + C5 in P = %.4f nb, C2 + C3 in B = %.4f nb" %
       (eff[("S0", "V3", "C4", "P", "all")][0] + eff[("S0", "V3", "C5", "P", "all")][0], eff[("S0", "V3", "C2", "B", "all")][0] + eff[("S0", "V3", "C3", "B", "all")][0]))
    P4("")
    P4("8. (pT4, y4) REWEIGHTING VARIANT (3 x 3 cells; wp8_pty.txt)")
    P4("   a_P, sigma_eta: %+.1f %% +- %.1f %% (noise); DR, BR ratios: %+.1f %% +- %.1f %%; a_eta (refit) %+.1f %%; a_B %+.1f %%" %
       (100 * pty["3x3"]["rel"], 100 * pty["3x3"]["noise"], 100 * pty["3x3"]["relDR"], 100 * pty["3x3"]["noiseDR"],
        100 * (v8[("pty3x3", "a_eta")]["v"] / v8[("nom", "a_eta")]["v"] - 1), 100 * (v8[("pty3x3", "a_B")]["v"] / v8[("nom", "a_B")]["v"] - 1)))
    P4("   information: pT only %+.1f %% +- %.1f %%, y only %+.1f %% +- %.1f %% (absolute scales); ratios change by < 2 %%" %
       (100 * pty["pT"]["rel"], 100 * pty["pT"]["noise"], 100 * pty["y"]["rel"], 100 * pty["y"]["noise"]))
    P4("")
    P4("9. LITERATURE (published inclusive sigma(pp -> eta X) near p_lab = 5.4 GeV/c: none found -> no independent comparison)")
    for kind, txt, v, e, ref in LIT:
        if isinstance(v, tuple):
            P4("   %s: %s mb  [%s]" % (txt, " / ".join("%.4f" % x for x in v), ref))
        else:
            P4("   %s: %.4f +- %.4f mb  [%s]" % (txt, v, e, ref))
    P4("   SMASH: 'too many eta mesons are produced from pp -> pp eta for sqrt(s) > 3.25 GeV' [J. Staudenmaier et al., Phys. Rev. C 98 (2018) 054908]")
    P4("   These are EXCLUSIVE pp -> pp eta values (a lower bound of the inclusive cross section), at T = 3.5 GeV or model predictions; not comparable directly.")
    with open(os.path.join(OUTW, "wp8_tableP2.txt"), "w") as f:
        f.write("\n".join(L4) + "\n")
    # ---- tex of Table P2 (short)
    def pm3(v, a, b, c, d=None, nd=3, sc=1.0):
        s = ("%%.%df \\pm %%.%df \\pm %%.%df \\pm %%.%df" % (nd, nd, nd, nd)) % (v * sc, a * sc, b * sc, c * sc)
        if d is not None:
            s += (" ^{+%%.%df}_{-%%.%df}" % (nd, nd)) % (d[0] * sc, d[1] * sc)
        return "$" + s + "$"
    with open(os.path.join(OUTW, "wp8_tableP2.tex"), "w") as f:
        f.write("% wp8_tableP2.tex -- plots/all_report/wp8_physics.py (WP8)\n")
        f.write("\\begin{table}[H]\n\\centering\n\\caption{Physics quantities (Table P2). S0, $\\mathrm{isBest}\\geq0$, $L=6079.2$\\,\\nbinv. Errors: "
                "statistical $\\pm$ systematic (Table S1, S12 completed by WP8) $\\pm$ template statistics (MC) [$\\pm$ normalisation: $L$ 6.9\\,\\%, "
                "$C_{\\mathrm{sec}}(1-f_\\eta)$ $^{+3.1}_{-4.1}$\\,\\%, PDG BRs]. Live time $+7.3/-0$\\,\\% on cross sections not included. MAIN: estimators that passed "
                "the closure test (WP4); CROSS-CHECK: F-main/F-gold/F-free (WP4 failed). Last column: change if class C2 had the QED $\\gamma^*$ decay angle "
                "(HGeant: isotropic; toy / closure).}\n\\label{tab:P2}\n\\scriptsize\n\\setlength{\\tabcolsep}{2.5pt}\n"
                "\\begin{tabular}{@{}l l l r@{}}\n\\toprule\nquantity & value & status & QED angle of C2 \\\\\n\\midrule\n")
        def angtxt(q):
            n = v8[("nom", q)]["v"]
            return "%+.0f / %+.0f\\,\\%%" % (100 * (v8[("ang_toy", q)]["v"] / n - 1), 100 * (v8[("ang_clos", q)]["v"] / n - 1))
        s = w7["BR_G"]; f.write("$\\mathrm{BR}_G\\ [10^{-4}]$ & %s & cross-check & %s \\\\\n" % (pm3(s["v"], s["stat"], sy("BR_G"), s["mc"], (s["nup"], s["ndn"]), 2, 1e4), angtxt("BR_G")))
        s = w7["BR_A"]; f.write("$\\mathrm{BR}_A\\ [10^{-4}]$ & %s & cross-check & %s \\\\\n" % (pm3(s["v"], s["stat"], sy("BR_A"), s["mc"], (s["nup"], s["ndn"]), 2, 1e4), angtxt("BR_A")))
        f.write("$\\mathrm{BR}_C\\ [10^{-4}]$ & $%.1f \\pm %.1f \\pm %.1f$ (E1+E2) & withdrawn & -- \\\\\n" % (1e4 * bc["v"], 1e4 * bc["stat"], 1e4 * bc["mc68"]))
        f.write("PDG $\\mathrm{BR}(\\eta\\to\\ppee)\\ [10^{-4}]$ & $2.68\\pm0.11$ & reference & \\\\\n")
        f.write("$R_{34}$ (E1+E2; E2: undetermined) & $%.2f\\pm%.2f\\pm%.2f$ & cross-check & \\\\\n" % (r34s["v"], r34s["stat"], r34s["mc"]))
        s = w7["kappa"]; f.write("$\\kappa$ (constrained $1\\pm0.2$) & %s & cross-check & %s \\\\\n" % (pm3(s["v"], s["stat"], sy("kappa"), s["mc"]), angtxt("kappa")))
        f.write("$\\kappa$ free / $\\kappa_B$ & $%.2f\\pm%.2f\\pm%.2f$ / $%.2f\\pm%.2f\\pm%.2f$ & cross-check & \\\\\n" % (kf["v"], kf["stat"], kf["mc"], kb["v"], kb["stat"], kb["mc"]))
        s = w7["DR"]; f.write("DR (side bands / B4) & %s & main (extrap.) & %s \\\\\n" % (pm3(s["v"], s["stat"], sy("DR"), s["mc"], (s["nup"], s["ndn"]), 2), angtxt("DR")))
        s = w7["DR_F"]; f.write("DR (F-main) & %s & cross-check & %s \\\\\n" % (pm3(s["v"], s["stat"], sy("DR_F"), s["mc"], (s["nup"], s["ndn"]), 2), angtxt("DR_F")))
        s = w7["sigma_eta"]; f.write("$\\sigma(pp\\to\\eta X)$ [mb] & %s & main & %s \\\\\n" % (pm3(s["v"], s["stat"], sy("sigma_eta"), s["mc"], (s["nup"], s["ndn"])), angtxt("sigma_eta")))
        s = w7["sigma_eta_F"]; f.write("$\\sigma(pp\\to\\eta X)$, F-main [mb] & %s & cross-check & %s \\\\\n" % (pm3(s["v"], s["stat"], sy("sigma_eta_F"), s["mc"], (s["nup"], s["ndn"])), angtxt("sigma_eta_F")))
        f.write("SMASH $\\sigma(pp\\to\\eta X)$ [mb] & $1.0665$ & model & \\\\\n")
        f.write("$\\sigma\\times\\mathrm{BR}_{\\mathrm{PDG}}(\\eta\\to\\ppee)$ [nb] & $%.0f$ ($\\pm%.0f$\\,\\%%) & main & \\\\\n" % (sxb, 100 * totsig / sig))
        s = w7["Y_peak"]; f.write("$\\sigma_{\\mathrm{vis}}$ peak [nb] & %s & main & \\\\\n" % pm3(s["v"] / L, s["stat"] / L, sy("Y_peak") / L, s["mc"] / L, None, 3))
        s = w7["Y_bump"]; f.write("$\\sigma_{\\mathrm{vis}}$ bump [nb] & %s & main (extrap.) & %s \\\\\n" % (pm3(s["v"] / L, s["stat"] / L, sy("Y_bump") / L, s["mc"] / L, None, 2), angtxt("Y_bump")))
        f.write("$(p_T,y)$ reweighting, $3\\times3$ & $%+.0f\\pm%.0f$\\,\\%% on $\\sigma$; $%+.1f\\pm%.1f$\\,\\%% on DR & variant & \\\\\n" %
                (100 * pty["3x3"]["rel"], 100 * pty["3x3"]["noise"], 100 * pty["3x3"]["relDR"], 100 * pty["3x3"]["noiseDR"]))
        f.write("\\bottomrule\n\\end{tabular}\n\\end{table}\n")

    # =================================================================================================================
    # 5. FIG. P1 (summary)
    # =================================================================================================================
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 13, "axes.labelsize": 15, "axes.titlesize": 14, "xtick.labelsize": 12, "ytick.labelsize": 12, "legend.fontsize": 10.5})
    fig, ax = plt.subplots(2, 3, figsize=(19, 11.5))
    # (a) eta scale relative to SMASH (x 0.981)
    A = ax[0, 0]
    ests = []
    def add(lbl, v, st, syt, mc, nrel, ang=None, col="k"):
        ests.append((lbl, v, st, syt, mc, nrel, ang, col))
    s = w7["a_P"]; add("peak side bands (MAIN)", s["v"] * CNORM, s["stat"] * CNORM, sy("a_P") * CNORM, s["mc"] * CNORM, 0.079, None, "C0")
    s = w7["a_B"]; add("bump B4 (extrapolated)", s["v"] * CNORM, s["stat"] * CNORM, sy("a_B") * CNORM, s["mc"] * CNORM, 0.079,
                       [v8[("ang_toy", "a_B")]["v"] * CNORM, v8[("ang_clos", "a_B")]["v"] * CNORM], "C1")
    s = w7["a_eta"]; add("F-main $a_\\eta$ (cross-check)", s["v"] * CNORM, s["stat"] * CNORM, sy("a_eta") * CNORM, s["mc"] * CNORM, 0.079,
                         [v8[("ang_toy", "a_eta")]["v"] * CNORM, v8[("ang_clos", "a_eta")]["v"] * CNORM], "C2")
    add("$f_\\eta$(edge) (cross-check)", fedge * CNORM, fst * CNORM, fmod * CNORM, fmc * CNORM, 0.079, None, "C3")
    add("$f_\\eta$(tag), emulated (cross-check)", float(ft[2]) * CNORM, float(ft[3]) * CNORM, 0, float(ft[4]) * CNORM, 0.079, None, "C4")
    add("$P_{\\mathrm{main}}$ (cross-check)", PMAIN_S / simS * CNORM, PMAIN_STAT / simS * CNORM, PMAIN_MODEL / simS * CNORM, 0, 0.079, None, "C5")
    for i, (lbl, v, st, syt, mc, nr, ang, col) in enumerate(ests):
        y = len(ests) - i
        t1 = st; t2 = math.hypot(st, syt); t3 = math.sqrt(st ** 2 + syt ** 2 + mc ** 2 + (nr * v) ** 2)
        A.plot([v - t3, v + t3], [y, y], color=col, lw=1.2)
        A.plot([v - t2, v + t2], [y, y], color=col, lw=4, alpha=0.6, solid_capstyle="butt")
        A.plot([v - t1, v + t1], [y, y], color=col, lw=8, solid_capstyle="butt")
        A.plot(v, y, "o", color="w", mec=col, ms=7, zorder=5)
        if ang:
            A.plot(ang, [y + 0.28, y + 0.28], "v", color=col, mfc="none", ms=7)
    A.set_yticks(range(len(ests), 0, -1)); A.set_yticklabels([e[0] for e in ests])
    A.axvline(1.0, color="gray", ls="--"); A.text(1.02, len(ests) + 0.35, "SMASH", ha="left", color="gray")
    A.set_xlim(0, 1.6); A.set_ylim(0.4, len(ests) + 0.7)
    A.set_xlabel("$\\sigma(pp\\to\\eta X)\\,/\\,\\sigma_{\\mathrm{SMASH}}$  [incl. $C_{\\mathrm{sec}}(1-f_\\eta)=0.981$]")
    A.set_title("(a) $\\eta$ scale relative to SMASH")
    A.text(0.98, 0.03, "thick: stat\nmedium: + syst\nthin: + MC + norm.\n$\\nabla$: QED angle of C2\n(toy, closure)", transform=A.transAxes, ha="right", va="bottom", fontsize=10)
    # (b) BR
    B = ax[0, 1]
    B.axvspan(2.68 - 0.11, 2.68 + 0.11, color="0.8"); B.axvline(2.68, color="0.4"); B.text(2.72, 0.55, "PDG", color="0.3")
    items = []
    s = w7["BR_G"]; items.append(("BR$_G$ (F-gold, $M_{ee}\\geq20$)", 1e4 * s["v"], 1e4 * s["stat"], 1e4 * sy("BR_G"), 1e4 * s["mc"], [1e4 * v8[("ang_toy", "BR_G")]["v"], 1e4 * v8[("ang_clos", "BR_G")]["v"]]))
    s = w7["BR_A"]; items.append(("BR$_A$ = 2.68 $\\rho_d$ (F-main)", 1e4 * s["v"], 1e4 * s["stat"], 1e4 * sy("BR_A"), 1e4 * s["mc"], [1e4 * v8[("ang_toy", "BR_A")]["v"], 1e4 * v8[("ang_clos", "BR_A")]["v"]]))
    items.append(("BR$_C$ (E1+E2 fit; withdrawn)", 1e4 * bc["v"], 1e4 * bc["stat"], 0, 1e4 * bc["mc68"], None))
    for i, (lbl, v, st, syt, mc, ang) in enumerate(items):
        y = len(items) - i
        B.plot([v - math.sqrt(st ** 2 + syt ** 2 + mc ** 2), v + math.sqrt(st ** 2 + syt ** 2 + mc ** 2)], [y, y], "k", lw=1.2)
        B.plot([v - math.hypot(st, syt), v + math.hypot(st, syt)], [y, y], "k", lw=4, alpha=0.5, solid_capstyle="butt")
        B.plot([v - st, v + st], [y, y], "k", lw=8, solid_capstyle="butt")
        B.plot(v, y, "o", color="w", mec="k", ms=7, zorder=5)
        if ang:
            B.plot(ang, [y + 0.25, y + 0.25], "v", color="C3", mfc="none", ms=7)
    B.set_yticks(range(len(items), 0, -1)); B.set_yticklabels([x[0] for x in items])
    B.set_xlim(-1, 9); B.set_ylim(0.4, len(items) + 0.7)
    B.set_xlabel("BR($\\eta\\to\\pi^+\\pi^-e^+e^-$)  [$10^{-4}$]")
    B.set_title("(b) BR($\\eta\\to\\pi^+\\pi^-e^+e^-$), cross-checks (WP4 fail)")
    # (c) sigma(pp -> eta X)
    C = ax[0, 2]
    C.axhline(1.0665, color="C1", ls="--"); C.text(0.02, 1.0665 + 0.02, "SMASH 1.0665 mb", color="C1")
    pts = [("this work\n(side bands)", sig, sst, ssy, smc, (snu + snd) / 2, "C0"),
           ("this work\n(F-main)", w7["sigma_eta_F"]["v"], w7["sigma_eta_F"]["stat"], sy("sigma_eta_F"), w7["sigma_eta_F"]["mc"], w7["sigma_eta_F"]["nup"], "C2")]
    for i, (lbl, v, st, syt, mc, nn, col) in enumerate(pts):
        x = i
        C.errorbar([x], [v], yerr=[[math.sqrt(st ** 2 + syt ** 2 + mc ** 2 + nn ** 2)], [math.sqrt(st ** 2 + syt ** 2 + mc ** 2 + nn ** 2 + (0.073 * v) ** 2)]], fmt="none", ecolor=col, lw=1.2, capsize=4)
        C.errorbar([x], [v], yerr=math.hypot(st, syt), fmt="none", ecolor=col, lw=4, alpha=0.6)
        C.errorbar([x], [v], yerr=st, fmt="o", color=col, mfc="w", lw=8, ms=8)
    C.plot([0, 0], [sig * (1 + pty["y"]["rel"]), sig * (1 + pty["pT"]["rel"])], color="C0", ls=":", lw=2)
    C.text(0.12, sig * (1 + pty["pT"]["rel"]), "1D $(p_T, y)$\nreweightings", color="C0", fontsize=10, va="center")
    # literature (exclusive)
    C.errorbar([1.95], [0.211], yerr=0.047, fmt="s", color="k", ms=7, capsize=3)
    C.errorbar([2.25], [0.1369], yerr=math.hypot(0.0009, 0.0101), fmt="D", color="k", mfc="w", ms=7, capsize=3)
    C.fill_between([2.95, 3.25], [0.0808, 0.0808], [0.1369, 0.1369], color="0.7")
    C.set_xticks([0, 1, 2.1, 3.1]); C.set_xticklabels(["incl.\nside bands", "incl.\nF-main", "$pp\\to pp\\eta$\n3.5 GeV", "$pp\\to pp\\eta$\nmodels"], fontsize=10)
    C.set_xlim(-0.5, 3.5); C.set_ylim(0, 1.2)
    C.set_ylabel("$\\sigma$ [mb]")
    C.set_title("(c) $\\sigma(pp\\to\\eta X)$, $T = 4.5$ GeV")
    C.text(0.98, 0.62, "exclusive values: Teilab 2011 (prelim., thesis);\nLebiedowicz 2025 (models)\nnot directly comparable", transform=C.transAxes, ha="right", fontsize=9.5, color="0.3")
    # (d) tagged shares
    D = ax[1, 0]
    xs = np.arange(4)
    D.errorbar(xs - 0.12, shD, yerr=esD, fmt="o", color="k", ms=8, label="data (ECAL tag, fit)")
    D.errorbar(xs, shS, yerr=esS, fmt="s", color="C3", mfc="none", ms=8, label="sim T1+T2, same fit (emulated)")
    D.errorbar(xs + 0.12, shT, yerr=esT, fmt="^", color="C0", ms=8, label="sim truth C2 + C3")
    gen4 = [hg["<8"], hg["8-20"], hg["20-30"], hg[">=30"]]
    D.plot(xs, gen4, "_", color="0.4", ms=30, mew=2, label="generator (KW), before acceptance")
    D.set_xticks(xs); D.set_xticklabels(["< 8", "8-20", "20-30", "30-100"])
    D.set_xlabel("$M(e^+e^-)$ class [MeV/$c^2$]"); D.set_ylabel("share of the tagged $\\eta$")
    D.set_ylim(0, 0.8); D.legend(loc="upper right", fontsize=9.5)
    D.set_title("(d) $\\varepsilon_2$ check: tagged $\\eta$ per $M(e^+e^-)$ (S0)")
    # (e) decay angle
    E = ax[1, 1]
    for k, col, lab in ((0, "k", "PLUTO, $M_{ee}<10$"), (2, "C0", "PLUTO, $M_{ee}\\geq20$")):
        h = np.array(plc[k], float); h = h / h.sum() * 20 / 2.0
        E.errorbar(cbin, h, yerr=np.sqrt(np.array(plc[k], float)) / sum(plc[k]) * 10, fmt="o", color=col, ms=5, label=lab)
    cc = np.linspace(-1, 1, 101)
    E.plot(cc, 3.0 / 8.0 * (1 + cc ** 2), "C3-", label="QED $\\frac{3}{8}(1+\\cos^2\\theta^*)$")
    E.plot(cc, 0.5 + 0 * cc, "C2--", label="isotropic (HGeant, $\\pi^0$ Dalitz bump)")
    E.set_xlabel("$\\cos\\theta^*$ ($e^+$ vs $\\gamma^*$ direction, $\\gamma^*$ frame)"); E.set_ylabel("$dN/d\\cos\\theta^*$ (normalised)")
    E.set_ylim(0.2, 0.9); E.legend(loc="upper center", fontsize=10)
    E.set_title("(e) $\\pi^0\\to\\gamma e^+e^-$ decay angle")
    E.text(0.5, 0.05, "efficiency of the $\\pi^0$ Dalitz bump\nx %.3f (toy) / %.3f (closure)\nat $M_{ee}\\geq20$: x %.3f / %.3f" % (0.8800, call, 0.8864, c20 / rq20),
           transform=E.transAxes, ha="center", fontsize=10.5)
    # (f) (pT4, y4) data/sim per cell
    F = ax[1, 2]
    import json
    pj = json.load(open(os.path.join(W8, "wp8_pty.json")))
    Y = np.array(pj["Y"]); VY = np.array(pj["VY"]); S = np.array(pj["S"]); g = np.array(pj["g"])
    r = Y / S; er = np.sqrt(VY) / S
    for a_ in range(3):
        F.errorbar(np.arange(3) + 0.1 * (a_ - 1), r[a_], yerr=er[a_], fmt="o-", ms=7, label="$p_T$ cell %d (%s)" % (a_ + 1, ["< 500", "500-680", "> 680"][a_]))
    F.axhline(pj["aP_cells"], color="k", ls="--"); F.text(2.25, pj["aP_cells"] + 0.03, "$a_P$", fontsize=11)
    F.set_xticks(range(3)); F.set_xticklabels(["$y_{cm}<-0.25$", "$-0.25$ to $0$", "$y_{cm}>0$"])
    F.set_ylabel("data / SMASH (peak, per cell)"); F.set_ylim(0, 1.2); F.legend(loc="upper right", fontsize=10)
    F.set_title("(f) $(p_T, y)$ of the peak: data/SMASH")
    F.text(0.03, 0.04, "reweighted production: $\\sigma$ x %.3f ($\\pm$%.0f%%)" % (1 + pty["3x3"]["rel"], 100 * pty["3x3"]["noise"]), transform=F.transAxes, fontsize=10)
    fig.suptitle("WP8 physics summary: data S0, isBest $\\geq$ 0, $L$ = 6079.2 nb$^{-1}$; SMASH $\\times$ PDG BR, V3 weights", fontsize=15)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG, "wp8_P1_plain." + ext), dpi=100 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)

    # =================================================================================================================
    # 6. machine-readable results
    # =================================================================================================================
    with open(os.path.join(OUTW, "results.txt"), "w") as f:
        f.write("# wp8 results.txt -- plots/all_report/wp8_physics.py (WP8, 2026-11-14)\n")
        f.write("\n".join(out) + "\n")
    print("written: tables_all/wp8_{tableP1,eps2,S12,tableP2}.*, figs_all/wp8_P1.*, wp8/results.txt")


if __name__ == "__main__":
    main()
