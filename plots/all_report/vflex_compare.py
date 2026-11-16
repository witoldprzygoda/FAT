#!/usr/bin/env python3
# vflex_compare.py -- compare the independent fits (vflex_fit.json) with the flex machine-readable results (flex_results.txt)
# Run: python3 -B plots/all_report/vflex_compare.py
import json
import numpy as np

FR = "/hdd3/przygoda/HADES/FAT_out/all_report/flex/flex_results.txt"
VJ = "/hdd3/przygoda/HADES/FAT_out/all_report/vflex/vflex_fit.json"


def parse_flex():
    out = {}
    lines = open(FR).read().split("\n")
    i = 0
    while i < len(lines):
        t = lines[i].split()
        if len(t) >= 2 and t[0] in ("values", "err_stat", "err_mc", "err_total", "mc_shift", "err_N", "err_CB"):
            out.setdefault(t[1], {})[t[0]] = np.array([float(x) for x in t[2:]])
        elif len(t) >= 2 and t[0] == "names":
            out.setdefault(t[1], {})["names"] = t[2:]
        elif len(t) >= 2 and t[0] == "gof":
            d = dict(zip(t[2::2], t[3::2]))
            out.setdefault(t[1], {})["gof"] = {k: float(v) for k, v in d.items()}
        elif len(t) == 3 and t[0] == "cor_total":
            n = int(t[2])
            out.setdefault(t[1], {})["cor"] = np.array([[float(x) for x in lines[i + 1 + q].split()] for q in range(n)])
            i += n
        i += 1
    return out


if __name__ == "__main__":
    F = parse_flex()
    Vj = json.load(open(VJ))
    for tag in ("P3", "F7", "F8", "F7bj", "F8b", "F9"):
        if tag not in F or tag not in Vj:
            continue
        f, v = F[tag], Vj[tag]
        assert f["names"] == v["names"], (f["names"], v["names"])
        th = np.array(v["th"]); st = np.array(v["stat"]); mc = np.array(v["mc"]); to = np.array(v["tot"]); sh = np.array(v["shift"])
        cor = np.array(v["cor"])
        print(f"{tag}: max |d value| {np.max(np.abs(th - f['values'])):.2e}, max rel d stat {np.max(np.abs(st / f['err_stat'] - 1)):.2e}, "
              f"max rel d MC {np.max(np.abs(mc / f['err_mc'] - 1)):.2e}, max rel d total {np.max(np.abs(to / f['err_total'] - 1)):.2e}, "
              f"max |d shift| {np.max(np.abs(sh - f['mc_shift'])):.2e}, max |d corr| {np.max(np.abs(cor - f['cor'])):.2e}; "
              f"hybrid GLS {v['gls']:.3f} vs {f['gof']['chi2GLShyb']:.3f}, data-only {v['glsd']:.3f} vs {f['gof']['chi2GLSdata']:.3f}, "
              f"Pearson {v['pearson']:.3f} vs {f['gof']['chi2P']:.3f}")
