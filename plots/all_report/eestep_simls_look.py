#!/usr/bin/env python3
# eestep_simls_look.py -- package eestep: first look at the SIMULATED like-sign rows (edgeu_dumpls.C columns):
# which true particles sit in the two lepton slots of the e-e- and e+e+ rows (charge-flipped leptons, ghosts, pions)
import sys, numpy as np
sys.path.insert(0, "/home/przygoda/HADES/FAT/FAT/plots/all_report")
import edgeu_core as C
O = C.OUT
n = int(open(O + "dumpls_layout.txt").read().split()[1])
def col(nm, dt): return np.fromfile(O + "simls_" + nm, dtype=dt)
S = {k: col(f, d) for k, f, d in [("tree", "tree.i8", np.int8), ("grp", "grp.i32", np.int32), ("tmpl", "tmpl.i8", np.int8),
     ("cls", "cls.i8", np.int8), ("lcls", "lcls.i8", np.int8), ("pcls", "pcls.i8", np.int8), ("flags", "flags.i32", np.int32),
     ("ghost", "ghost.i32", np.int32), ("bits", "bits.u8", np.uint8), ("id1", "id1.i32", np.int32), ("id2", "id2.i32", np.int32),
     ("par1", "par1.i32", np.int32), ("par2", "par2.i32", np.int32), ("gpar1", "gpar1.i32", np.int32), ("gpar2", "gpar2.i32", np.int32),
     ("proc1", "proc1.i32", np.int32), ("proc2", "proc2.i32", np.int32), ("m4", "m4.f32", np.float32), ("m4t", "m4true.f32", np.float32),
     ("mpipi", "mpipi.f32", np.float32), ("mpipit", "mpipitrue.f32", np.float32), ("mee", "mee.f32", np.float32),
     ("meet", "meetrue.f32", np.float32), ("w", "w.f64", np.float64), ("f108", "f108.f32", np.float32)]}
S["n"] = n
u, mpp, m4p = C.u_sim(S)
w = S["w"] * S["f108"]
s0 = (S["bits"] & 1) != 0; ownA = (S["bits"] & 8) != 0
base = s0 & ownA & (mpp >= 0.30) & (mpp < 0.40) & np.isfinite(u)
tn = {1: "std_pp", 2: "std_mm", 4: "lep_pp", 5: "lep_mm", 7: "eta_pp", 8: "eta_mm"}
for t in (2, 5, 8, 1, 4, 7):
    m = base & (S["tree"] == t)
    print("tree %s: rows %d (S0, ownA, M(pipi)' 0.30-0.40), sum w x f108 = %.1f raw rows" % (tn[t], m.sum(), w[m].sum()))
    pairs = {}
    for a, b in zip(S["id1"][m], S["id2"][m]):
        pairs[(a, b)] = pairs.get((a, b), 0) + 1
    top = sorted(pairs.items(), key=lambda x: -x[1])[:8]
    print("   (id1,id2) counts:", top)
    near = m & (u > -0.045) & (u < 0.045)
    print("   in -45..45 MeV: rows %d, w %.1f; lcls counts:" % (near.sum(), w[near].sum()), np.bincount(S["lcls"][near], minlength=12))
    wrong = near & (((S["id1"] == 2) | (S["id2"] == 2)) if t in (2, 5, 8) else ((S["id1"] == 3) | (S["id2"] == 3)))
    print("   rows with a wrong-charge true lepton (e+ in an e- slot or vice versa): %d rows, w %.2f; their parents:" % (wrong.sum(), w[wrong].sum()),
          list(zip(S["id1"][wrong][:20], S["par1"][wrong][:20], S["id2"][wrong][:20], S["par2"][wrong][:20])))
    gh = near & ((S["flags"] & (1 << 10)) != 0)
    same = near & ((S["flags"] & (1 << 12)) != 0)
    print("   ghost-lepton rows %d (w %.1f), same-particle rows %d (w %.1f)" % (gh.sum(), w[gh].sum(), same.sum(), w[same].sum()))
