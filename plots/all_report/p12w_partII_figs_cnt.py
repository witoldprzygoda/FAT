#!/usr/bin/env python3
# ============================================================================
# p12w_partII_figs_cnt.py  --  COPY of p12w_partII_figs.py made by figure agent A (integ4, 2026-11-15):
# the analyst's terminology "counts" instead of "rows" in the figure labels. Everything else as in the original:
# the two macros are read, the listed strings are replaced (each must match exactly the stated number of times),
# and the result is executed. Differences to the original wrapper:
#   - extra replacements rows -> counts in the axis titles, panel titles, legends of the two Part II figures;
#   - logs, json and all eestep figures go to the scratch folder OUT below (nothing of /hdd3 is touched);
#   - written into docs/latex/figs_exp/ ONLY: p12w_edgeu_data_step_cnt.{pdf,png} and p12w_eestep_spectra_cnt.{pdf,png}
#     (p12w_eestep_align has no "rows" and is not copied).
# The numbers printed into the new log must be identical to /hdd3/przygoda/HADES/FAT_out/all_report/p12w/p12w_edgeu_partII.log.
# Run: python3 -I plots/all_report/p12w_partII_figs_cnt.py
# ---------------------------------------------------------------------------- original header follows
# p12w_partII_figs.py  --  Part II figures with the edge function written M_edge(m, q)
# Written 2026-11-15 (prefix p12w, Parts I/II writer). Part II uses the letter E for its per-slice bump measure
# (E = 1199 rows); the kinematic edge function E(m, q) of Part III is therefore written M_edge(m, q) in Part II.
# This script re-runs, unchanged except for the axis labels and the output names, the two macros that made the
# Part II figures with an "E(...)" label:
#   plots/all_report/edgeu_partII_datastep.py -> docs/latex/figs_exp/p12w_edgeu_data_step.{pdf,png}
#   plots/all_report/eestep_figs.py           -> docs/latex/figs_exp/p12w_eestep_{spectra,align}.{pdf,png}
# ============================================================================
import os
import shutil
import sys

HERE = "/home/przygoda/HADES/FAT/FAT/plots/all_report/"
OUT = "/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figA/work/p12w/"
FIGX = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_exp/"
os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, HERE)


def patched(name, reps):
    src = open(HERE + name).read()
    for old, new, n in reps:
        c = src.count(old)
        if c != n:
            raise SystemExit("%s: '%s' found %d times, expected %d" % (name, old, c, n))
        src = src.replace(old, new)
    return src


def run(name, reps):
    src = patched(name, reps)
    ns = {"__name__": "__main__", "__file__": HERE + name}
    exec(compile(src, HERE + name + " (p12w cnt patch)", "exec"), ns)


# ---------------------------------------------------------------- 1. edgeu_partII_datastep.py
for ext in ("pdf", "png"):
    if os.path.exists(FIGX + "p12w_edgeu_data_step_cnt." + ext):
        print("note: overwriting own earlier output p12w_edgeu_data_step_cnt." + ext)
run("edgeu_partII_datastep.py", [
    ('OUT = "/hdd3/przygoda/HADES/FAT_out/all_report/edgeu_partII/"', 'OUT = "%s"' % OUT, 1),
    ('"edgeu_partII.log"', '"p12w_edgeu_partII.log"', 1),
    ('"edgeu_partII.json"', '"p12w_edgeu_partII.json"', 1),
    ('FIGX + "edgeu_data_step.pdf"', 'FIGX + "p12w_edgeu_data_step_cnt.pdf"', 1),
    ('FIGX + "edgeu_data_step.png"', 'FIGX + "p12w_edgeu_data_step_cnt.png"', 1),
    (r'"$u=M_4-E(M_{\\pi\\pi},M_{ee})$ [MeV/$c^2$]"', r'"$u=M_4-M_{\\mathrm{edge}}(M_{\\pi\\pi},M_{ee})$ [MeV/$c^2$]"', 1),
    # cnt: figure labels only
    ('ax.set_ylabel("rows per 2.5 MeV/$c^2$")', 'ax.set_ylabel("counts per 2.5 MeV/$c^2$")', 3),
    ('ax.set_title("(c) like-sign rows alone (same fit)")', 'ax.set_title("(c) like-sign counts alone (same fit)")', 1),
])

# ---------------------------------------------------------------- 2. eestep_figs.py (reads eestep_ana.json, eestep_figdata.npz)
run("eestep_figs.py", [
    ('fig.savefig(OUT + name + ".pdf"); fig.savefig(OUT + name + ".png", dpi=110)',
     'fig.savefig("%sp12w_" + name + ".pdf"); fig.savefig("%sp12w_" + name + ".png", dpi=110)' % (OUT, OUT), 1),
    (r'"u = M$_4$ - E(M($\\pi\\pi$), M(ee))  [MeV/c$^2$]"', r'"u = M$_4$ - M$_{\\rm edge}$(M($\\pi\\pi$), M(ee))  [MeV/c$^2$]"', 1),
    (r'"u = M$_4$ - E(m, q)"', r'"u = M$_4$ - M$_{\\rm edge}$(m, q)"', 1),
    (r'"M$_4$ - E(m, 0)"', r'"M$_4$ - M$_{\\rm edge}$(m, 0)"', 1),
    # cnt: labels of the spectra figure only (the other eestep figures stay in the scratch folder)
    ('label="data rows")', 'label="data counts")', 1),
    ('a.set_ylabel("rows / 2.5 MeV/c$^2$")', 'a.set_ylabel("counts / 2.5 MeV/c$^2$")', 1),
    ('fig.suptitle("Like-sign rows (S0,', 'fig.suptitle("Like-sign counts (S0,', 1),
])
for nm in ("eestep_spectra",):
    for ext in ("pdf", "png"):
        dst = FIGX + "p12w_" + nm + "_cnt." + ext
        shutil.copyfile(OUT + "p12w_" + nm + "." + ext, dst)
print("p12w cnt: figures written to", FIGX)
