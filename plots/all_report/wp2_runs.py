#!/usr/bin/env python3
# ============================================================================
# wp2_runs.py  --  WP2: time span of the sector-2 deficit from the per-run counts of wp2_runs.C
# ============================================================================
# eps(block) = (n_s2/n_no2)(block) / (n_s2/n_no2)(reference), reference = days 047-067 (no deficit in L, and
# no deficit seen here); n = S0 unlike-sign rows (and like-sign rows as a second sample). No luminosity needed.
# Blocks: runs in time order (from the DST stem be22DDDhhmmss), one block = consecutive runs with ~ BLOCK_L nb^-1.
# Output: docs/latex/tables_all/wp2_runs.txt and /hdd3/.../wp2/wp2_blocks.csv (for Fig. N1)
# Usage: python3 -I plots/all_report/wp2_runs.py [S0|pre|none]
# ============================================================================
import csv, math, sys

SEL = sys.argv[1] if len(sys.argv) > 1 else "S0"
REMOTE = "/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/remote"
OUT = "/hdd3/przygoda/HADES/FAT_out/all_report/wp2"
TAB = "/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all"
BLOCK_L = 25.0

lum = {}
with open(REMOTE + "/lumi_startfree_perrun.csv") as f:
    for r in csv.DictReader(l for l in f if not l.startswith("#")):
        lum[int(r["runnumber"])] = (r["dst_stem"], float(r["L6_SF"]))

runs = []
with open(OUT + "/wp2_runs_%s.csv" % SEL) as f:
    for r in csv.DictReader(f):
        run = int(r["run"])
        if int(r["day"]) == 0: continue   # a run whose rows all failed the selection (entry made by the variance flush)
        stem, L = lum.get(run, ("?", 0.0))
        runs.append(dict(run=run, day=int(r["day"]), stem=stem, L=L,
                         s2=float(r["us_s2"]), v2=float(r["us_vs2"]), n0=float(r["us_no2"]), v0=float(r["us_vno2"]),
                         a12=float(r["us_A12"]), va12=float(r["us_vA12"]), a34=float(r["us_A34"]), va34=float(r["us_vA34"]),
                         b012=float(r["us_B012"]), vb012=float(r["us_vB012"]), b3=float(r["us_B3"]), vb3=float(r["us_vB3"]),
                         ls2=float(r["pp_s2"]) + float(r["mm_s2"]), lv2=float(r["pp_vs2"]) + float(r["mm_vs2"]),
                         ln0=float(r["pp_no2"]) + float(r["mm_no2"]), lv0=float(r["pp_vno2"]) + float(r["mm_vno2"])))
missing = [r for r in runs if r["stem"] == "?"]
runs.sort(key=lambda r: r["stem"])

def odds(rs, k2="s2", k0="n0", v2="v2", v0="v0"):
    a = sum(r[k2] for r in rs); b = sum(r[k0] for r in rs)
    va = sum(r[v2] for r in rs); vb = sum(r[v0] for r in rs)
    o = a / b
    return o, o * math.sqrt(va / a ** 2 + vb / b ** 2), a, b

ref = [r for r in runs if r["day"] >= 47]
oref, eref, _, _ = odds(ref)
olref, elref, _, _ = odds(ref, "ls2", "ln0", "lv2", "lv0")

out = open(TAB + "/wp2_runs.txt", "w")
out.write("# wp2_runs.py (%s): eps = (n_s2/n_no2)(group) / (n_s2/n_no2)(days 047-067); S0 rows; no luminosity used\n" % SEL)
out.write("# reference odds (us) %.5f +- %.5f; like-sign %.5f +- %.5f; runs without a luminosity entry: %d\n" % (oref, eref, olref, elref, len(missing)))
out.write("\n[per day] day  L  n_s2  n_no2  eps(us) +- stat   eps(like-sign) +- stat\n")
days = sorted(set(r["day"] for r in runs))
perday = {}
for d in days:
    rs = [r for r in runs if r["day"] == d]
    o, e, a, b = odds(rs)
    ol, el, _, _ = odds(rs, "ls2", "ln0", "lv2", "lv0")
    eps = o / oref; ee = eps * math.sqrt((e / o) ** 2 + (eref / oref) ** 2)
    epl = ol / olref; eel = epl * math.sqrt((el / ol) ** 2 + (elref / olref) ** 2)
    perday[d] = (eps, ee)
    out.write("  %03d %8.3f %7.0f %7.0f   %.4f +- %.4f   %.4f +- %.4f\n" % (d, sum(r["L"] for r in rs), a, b, eps, ee, epl, eel))

# groups
def grp(name, sel):
    rs = [r for r in runs if sel(r)]
    o, e, a, b = odds(rs)
    eps = o / oref; ee = eps * math.sqrt((e / o) ** 2 + (eref / oref) ** 2)
    out.write("  %-34s L %8.3f  eps %.4f +- %.4f\n" % (name, sum(r["L"] for r in rs), eps, ee))
    return eps, ee
out.write("\n[groups] (reference 047-067)\n")
G = {}
G["034-044"] = grp("days 034-044", lambda r: 34 <= r["day"] <= 44)
G["034-045"] = grp("days 034-045", lambda r: 34 <= r["day"] <= 45)
G["045"] = grp("day 045", lambda r: r["day"] == 45)
G["046"] = grp("day 046", lambda r: r["day"] == 46)
G["045-067"] = grp("days 045-067", lambda r: r["day"] >= 45)

# time blocks
blocks = []
cur = []; Lc = 0.0
for r in runs:
    if cur and (Lc >= BLOCK_L or r["day"] != cur[-1]["day"]):
        blocks.append(cur); cur = []; Lc = 0.0
    cur.append(r); Lc += r["L"]
if cur: blocks.append(cur)
out.write("\n[time blocks] ~%.0f nb^-1 each, within a day; first/last DST stem, L, eps(us)\n" % BLOCK_L)
bo = open(OUT + "/wp2_blocks_%s.csv" % SEL, "w")
bo.write("idx,day,first,last,L,eps,err\n")
for i, b in enumerate(blocks):
    o, e, a, bb = odds(b)
    eps = o / oref; ee = eps * math.sqrt((e / o) ** 2 + (eref / oref) ** 2)
    bo.write("%d,%03d,%s,%s,%.4f,%.5f,%.5f\n" % (i, b[0]["day"], b[0]["stem"], b[-1]["stem"], sum(r["L"] for r in b), eps, ee))
    if 43 <= b[0]["day"] <= 47:
        out.write("  %03d %s..%s  L %6.2f  eps %.4f +- %.4f\n" % (b[0]["day"], b[0]["stem"][7:], b[-1]["stem"][7:], sum(r["L"] for r in b), eps, ee))
bo.close()
# transition: last deficit run in day 045/046 by a scan of the split point (max likelihood of a step)
out.write("\n[step scan] days 044-046: split point t (DST stem); eps before / after; chi2 of a two-level model (before = e1, after = 1)\n")
seg = [r for r in runs if 44 <= r["day"] <= 46]
best = None
for k in range(5, len(seg) - 5, 5):
    o1, e1, _, _ = odds(seg[:k]); o2, e2, _, _ = odds(seg[k:])
    eps1 = o1 / oref; eps2 = o2 / oref
    # chi2 of 'after = 1' plus 'before = eps1' (free) -> only the after part tests
    c2 = ((eps2 - 1.0) / (eps2 * math.sqrt((e2 / o2) ** 2 + (eref / oref) ** 2))) ** 2
    if best is None or c2 < best[0]:
        best = (c2, k, seg[k - 1]["stem"], seg[k]["stem"], eps1, eps2)
    if k % 25 == 0:
        out.write("  split after %s: eps(044..split) %.4f, eps(split..046) %.4f\n" % (seg[k - 1]["stem"], eps1, eps2))
# blocks with a significant deficit anywhere (eps < 1 - 4 sigma)
out.write("\n[deficit blocks] all days: blocks with eps < 1 - 4 sigma (reference 047-067)\n")
nlow = {}
for b in blocks:
    o, e, a, bb = odds(b)
    eps = o / oref; ee = eps * math.sqrt((e / o) ** 2 + (eref / oref) ** 2)
    d = b[0]["day"]
    nlow.setdefault(d, [0, 0])
    nlow[d][1] += 1
    if eps < 1 - 4 * ee:
        nlow[d][0] += 1
        out.write("  %03d %s..%s  L %6.2f  eps %.4f +- %.4f\n" % (d, b[0]["stem"][7:], b[-1]["stem"][7:], sum(r["L"] for r in b), eps, ee))
out.write("# per day: deficit blocks / blocks: " + " ".join("%03d:%d/%d" % (d, nlow[d][0], nlow[d][1]) for d in sorted(nlow)) + "\n")
# trigger-model categories in the deficit blocks of days 034-045 vs the reference: A12 and B3 relative to A34
out.write("\n[categories in deficit blocks] days 034-045, blocks with eps < 0.85: odds relative to A34, normalised to days 047-067\n")
low = [b for b in blocks if b[0]["day"] <= 45 and odds(b)[0] / oref < 0.85]
rl = [r for b in low for r in b]
for k, vk in (("a12", "va12"), ("b012", "vb012"), ("b3", "vb3"), ("s2", "v2")):
    o1, e1, _, _ = odds(rl, k, "a34", vk, "va34"); o0, e0, _, _ = odds(ref, k, "a34", vk, "va34")
    out.write("  %-5s eps %.4f +- %.4f  (%d blocks, L %.1f)\n" % (k, o1 / o0, (o1 / o0) * math.sqrt((e1 / o1) ** 2 + (e0 / o0) ** 2), len(low), sum(r["L"] for r in rl)))
out.write("# split with the 'after' part closest to 1: after %s (next %s): eps before %.4f, after %.4f, chi2(after = 1) %.2f\n" % (best[2], best[3], best[4], best[5], best[0]))
out.close()
print(open(TAB + "/wp2_runs.txt").read())
