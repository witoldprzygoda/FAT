# x17p26_X2_plain3.py -- COPY of x17p26_X2_plain2.py (2026-11-16, writer of report_all_short.tex): legend words only:
# "frozen ..." -> "expected ..." (the caption says the expectation was fixed before the signal region was opened), "primary" -> "nominal",
# "fall-back" -> "alternative", "toy" -> "pseudo-experiments". Output: docs/latex/figs_all/fix_X2_pdg26_plain3.pdf/.png. Every number unchanged.
#!/usr/bin/env python3
# x17p26_X2_plain2.py (2026-11-16, short-report figure agent, second pass, for report_all_short.tex): COPY of x17p26_X2_plain.py;
# the two fit-parameter symbols of the lower legend in words (instruction I2): "eta scale a_eta = 0.38" -> "the eta scale 0.38",
# "peak a_P = 0.41" -> "the peak scale 0.41". Every number unchanged. Output: docs/latex/figs_all/fix_X2_pdg26_plain3.pdf/.png.
# ============================================================================
# x17p26_X2_plain.py (2026-11-16, short-report figure agent for report_all_short.tex): COPY of x17p26_X2.py; only label strings
# in plain words (readability instruction I2): normalisation "(1)" / "(2)" / "(3)" -> "tied to fitted direct decays (primary)" /
# "tied to eta scale (alternative)" / "tied to peak"; "S0, isBest >= 0, X2D" -> "both pions in the vertex fit" (the fit plane is
# given in the caption); "spec. fall-back primary" -> "fall-back". Same tables, same drawing code and axes, every number unchanged.
# Output: docs/latex/figs_all/fix_X2_pdg26_plain.pdf/.png.   Run: python3 -I plots/all_report/x17p26_X2_plain.py   (seconds)
# ============================================================================
# ============================================================================
# x17p26_X2.py  --  X17 limit with PDG 2026 (analyst ruling D25, 15 November 2026), prefix x17p26, written 2026-11-15.
# Copy of the X2 block of fix_figs.py (Fig. X2: limit against m_X, narrowed y range, variants in a ratio pad) drawn from the
# PDG 2026 tables of x17p26_tables.py: x17d_scan_pdg26.txt, x17c_expected_asimov_pdg26.txt, x17c_expected_pdg26.txt, x17d_D1_pdg26.txt.
# Output: docs/latex/figs_all/fix_X2_pdg26.pdf/.png (fix_X2.pdf is not touched).
# Run: python3 -I plots/all_report/x17p26_X2.py   (seconds)
# ============================================================================
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FIG = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_all/"
TAB = "/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all/"
plt.rcParams.update({"font.size": 13, "axes.labelsize": 15, "axes.titlesize": 14, "legend.fontsize": 11, "xtick.labelsize": 12,
                     "ytick.labelsize": 12, "pdf.fonttype": 42})

def rows(fn):
    for l in open(TAB + fn):
        if l.startswith("#") or not l.strip(): continue
        yield l.replace("|", " ").split()
def scan():
    d = np.array([[float(x) for x in p[:13]] for p in rows("x17d_scan_pdg26.txt")])
    return dict(m=d[:, 0], ul1=d[:, 2], ul2=d[:, 3], ul3=d[:, 4], cos2=d[:, 11])
def asimov():
    d = np.array([[float(x) for x in p[:7]] for p in rows("x17c_expected_asimov_pdg26.txt")])
    return dict(m=d[:, 0], m2=d[:, 2], m1=d[:, 3], med=d[:, 4], p1=d[:, 5], p2=d[:, 6])
def toys():
    t = []
    for l in open(TAB + "x17c_expected_pdg26.txt"):
        p = l.split()
        if len(p) >= 9 and p[1] == "1" and p[2] == "toys":
            t.append([float(p[0]), float(p[5]), float(p[6]), float(p[7])])
    return np.array(t)
def toycls():
    t = []
    for l in open(TAB + "x17d_D1_pdg26.txt"):
        p = l.split()
        if len(p) > 4 and not l.startswith("#") and p[1] == "1":
            try: t.append((float(p[0]), float(p[3])))
            except ValueError: pass
    return t
S = scan(); A = asimov(); T = toys(); TC = toycls()
i17 = int(np.argmin(np.abs(A["m"] - 17.0))); t17 = T[np.argmin(np.abs(T[:, 0] - 17.0))]
fig = plt.figure(figsize=(11, 9.5))
ax = fig.add_axes([0.12, 0.40, 0.84, 0.55]); ar = fig.add_axes([0.12, 0.08, 0.84, 0.30], sharex=ax)
ax.fill_between(A["m"], A["m2"], A["p2"], color="#ffd23f", lw=0, label="expected, $\\pm2\\sigma$ (Asimov)")
ax.fill_between(A["m"], A["m1"], A["p1"], color="#3cc23c", lw=0, label="expected, $\\pm1\\sigma$ (Asimov)")
ax.plot(A["m"], A["med"], "k--", lw=1.6, label="expected: Asimov median (%.2f at 17 MeV/$c^2$)" % (1e5 * A["med"][i17]))
ax.errorbar(T[:, 0], T[:, 2], yerr=[T[:, 2] - T[:, 1], T[:, 3] - T[:, 2]], fmt="o", color="#006400", mfc="white", ms=7,
            label="expected: median of pseudo-experiments, 16-84 %% (%.2f at 17 MeV/$c^2$)" % (1e5 * t17[2]))
ax.plot(S["m"], S["ul1"], "k-", lw=2.4, marker=".", ms=7, label="observed (nominal: tied to the fitted direct decays)")
ax.plot([x for x, _ in TC], [y for _, y in TC], "*", color="#e8202a", ms=13, label="observed (nominal), CL$_s$ from pseudo-experiments")
ax.set_yscale("log"); ax.set_ylim(7e-6, 1.3e-4); ax.set_xlim(9.5, 40.5)
ax.set_ylabel("90 % CL$_s$ UL on BR$(\\eta\\to\\pi^+\\pi^-X)$BR$(X\\to e^+e^-)$")
ax.legend(loc="upper right", fontsize=10.5, frameon=False, ncol=1)
ax.text(10.0, 8.2e-6, "both pions in the vertex fit; $L=6079.2$ nb$^{-1}$; BR$_{\\mathrm{dir}}=2.76\\times10^{-4}$ (PDG 2026); spin-0 hypothesis", fontsize=10)
plt.setp(ax.get_xticklabels(), visible=False)
med = np.interp(S["m"], A["m"], A["med"])
ar.fill_between(A["m"], A["m2"] / A["med"], A["p2"] / A["med"], color="#ffd23f", lw=0)
ar.fill_between(A["m"], A["m1"] / A["med"], A["p1"] / A["med"], color="#3cc23c", lw=0)
ar.axhline(1, color="k", ls="--", lw=1.2)
ar.plot(S["m"], S["ul1"] / med, "k-", lw=2.2, label="tied to the fitted direct decays (nominal)")
ar.plot(S["m"], S["ul2"] / med, "--", color="#1f3fbf", lw=1.8, label="tied to the $\\eta$ scale 0.38 (alternative)")
ar.plot(S["m"], S["ul3"] / med, ":", color="#d400d4", lw=2.2, label="tied to the peak scale 0.41")
ar.plot(S["m"], S["cos2"] / med, "-.", color="#e8202a", lw=1.8, label="nominal, $e^+$ angle $1+\\cos^2\\theta$ (vector-like)")
ar.set_ylim(0.3, 2.6); ar.set_xlabel("$m_X$ [MeV/$c^2$]"); ar.set_ylabel("UL / expected Asimov median")
ar.legend(loc="upper right", fontsize=9.5, frameon=False, ncol=2)
fig.savefig(FIG + "fix_X2_pdg26_plain3.pdf"); fig.savefig(FIG + "fix_X2_pdg26_plain3.png", dpi=80)
print("x17p26_X2: done; Asimov median 17: %.4e, toy median 17: %.4e, toy CLs: %s" % (A["med"][i17], t17[2], TC))
