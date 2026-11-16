#!/usr/bin/env python3
# ============================================================================
# p12w_partII_figs.py  --  Part II figures with the edge function written M_edge(m, q)
# ============================================================================
# Written 2026-11-15 (prefix p12w, Parts I/II writer). Part II uses the letter E for its per-slice bump measure
# (E = 1199 rows); the kinematic edge function E(m, q) of Part III is therefore written M_edge(m, q) in Part II.
# This script re-runs, unchanged except for the axis labels and the output names, the two macros that made the
# Part II figures with an "E(...)" label:
#   plots/all_report/edgeu_partII_datastep.py -> docs/latex/figs_exp/p12w_edgeu_data_step.{pdf,png}
#   plots/all_report/eestep_figs.py           -> docs/latex/figs_exp/p12w_eestep_{spectra,align}.{pdf,png}
# The source of each macro is read, the listed strings are replaced (each replacement must match exactly the
# stated number of times), and the result is executed. Logs and json go to /hdd3/.../all_report/p12w/ with the
# prefix p12w_; nothing of the original packages is overwritten. The numbers printed into the new log must be
# identical to /hdd3/przygoda/HADES/FAT_out/all_report/edgeu_partII/edgeu_partII.log (checked with diff).
# Run: python3 -I plots/all_report/p12w_partII_figs.py
# ============================================================================
import os
import shutil
import sys

HERE = "/home/przygoda/HADES/FAT/FAT/plots/all_report/"
OUT = "/hdd3/przygoda/HADES/FAT_out/all_report/p12w/"
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
    exec(compile(src, HERE + name + " (p12w patch)", "exec"), ns)


# ---------------------------------------------------------------- 1. edgeu_partII_datastep.py
run("edgeu_partII_datastep.py", [
    ('OUT = "/hdd3/przygoda/HADES/FAT_out/all_report/edgeu_partII/"', 'OUT = "%s"' % OUT, 1),
    ('"edgeu_partII.log"', '"p12w_edgeu_partII.log"', 1),
    ('"edgeu_partII.json"', '"p12w_edgeu_partII.json"', 1),
    ('FIGX + "edgeu_data_step.pdf"', 'FIGX + "p12w_edgeu_data_step.pdf"', 1),
    ('FIGX + "edgeu_data_step.png"', 'FIGX + "p12w_edgeu_data_step.png"', 1),
    (r'"$u=M_4-E(M_{\\pi\\pi},M_{ee})$ [MeV/$c^2$]"', r'"$u=M_4-M_{\\mathrm{edge}}(M_{\\pi\\pi},M_{ee})$ [MeV/$c^2$]"', 1),
])

# ---------------------------------------------------------------- 2. eestep_figs.py (reads eestep_ana.json, eestep_figdata.npz)
run("eestep_figs.py", [
    ('fig.savefig(OUT + name + ".pdf"); fig.savefig(OUT + name + ".png", dpi=110)',
     'fig.savefig("%sp12w_" + name + ".pdf"); fig.savefig("%sp12w_" + name + ".png", dpi=110)' % (OUT, OUT), 1),
    (r'"u = M$_4$ - E(M($\\pi\\pi$), M(ee))  [MeV/c$^2$]"', r'"u = M$_4$ - M$_{\\rm edge}$(M($\\pi\\pi$), M(ee))  [MeV/c$^2$]"', 1),
    (r'"u = M$_4$ - E(m, q)"', r'"u = M$_4$ - M$_{\\rm edge}$(m, q)"', 1),
    (r'"M$_4$ - E(m, 0)"', r'"M$_4$ - M$_{\\rm edge}$(m, 0)"', 1),
])
for nm in ("eestep_spectra", "eestep_align"):
    for ext in ("pdf", "png"):
        dst = FIGX + "p12w_" + nm + "." + ext
        if os.path.exists(dst):
            raise SystemExit("refusing to overwrite " + dst)
        shutil.copyfile(OUT + "p12w_" + nm + "." + ext, dst)
print("p12w: figures written to", FIGX)
