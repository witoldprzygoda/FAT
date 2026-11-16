#!/usr/bin/env python3
# vx_dump.py -- independent X17 verification (prefix vx): READ-ONLY dump of the X2D inputs into a numpy file.
# Reads hists_v1.root (side bands, CB, like-sign, simulation templates S0__nom), the sealed file (rows 12-28 MeV/c2) and
# x17_lines.root (C6 line). Writes /hdd3/przygoda/HADES/FAT_out/all_report/vx/vx_inputs.npz. Nothing else is written.
import numpy as np
import ROOT
ROOT.gROOT.SetBatch(True)
H = "/hdd3/przygoda/HADES/FAT_out/all_report/hists_v1.root"
S = "/hdd3/przygoda/HADES/FAT_out/all_report/sealed/hists_v1_x17_sealed.root"
OUT = "/hdd3/przygoda/HADES/FAT_out/all_report/vx/vx_inputs.npz"

def arr(h, flows=False):
    nx, ny = h.GetNbinsX(), h.GetNbinsY()
    if flows:
        a = np.zeros((nx + 2, ny + 2)); e = np.zeros((nx + 2, ny + 2))
        for i in range(nx + 2):
            for j in range(ny + 2):
                a[i, j] = h.GetBinContent(i, j); e[i, j] = h.GetBinError(i, j)
    else:
        a = np.zeros((nx, ny)); e = np.zeros((nx, ny))
        for i in range(nx):
            for j in range(ny):
                a[i, j] = h.GetBinContent(i + 1, j + 1); e[i, j] = h.GetBinError(i + 1, j + 1)
    return a, e

f = ROOT.TFile.Open(H, "READ")
fs = ROOT.TFile.Open(S, "READ")
d = {}
hN = f.Get("data/S0/all/N_X2D")
d["xe"] = np.array([hN.GetXaxis().GetBinLowEdge(i) for i in range(1, hN.GetNbinsX() + 2)])
d["ye"] = np.array([hN.GetYaxis().GetBinLowEdge(j) for j in range(1, hN.GetNbinsY() + 2)])
d["N_sb"], _ = arr(hN)
d["N_sealed"], _ = arr(fs.Get("data/S0/all/N_X2D"))
for k in ["CB_X2D", "Npp_X2D", "Nmm_X2D"]:
    d[k], d[k + "_err"] = arr(f.Get("data/S0/all/" + k))
for c in ["C2", "C3", "C4", "C5", "C6", "Geta", "GD", "Gconv"]:
    a, e = arr(f.Get("sim/S0__nom/%s_X2D" % c), flows=True)
    d["sim_" + c] = a; d["sim_" + c + "_err"] = e
# exact window counts known before the freeze (N_win, CB_win)
hw = f.Get("data/S0/all/N_win"); hc = f.Get("data/S0/all/CB_win")
d["Nwin"] = np.array([hw.GetBinContent(i) for i in range(1, hw.GetNbinsX() + 1)])
d["CBwin"] = np.array([hc.GetBinContent(i) for i in range(1, hc.GetNbinsX() + 1)])
d["win_labels"] = np.array([hw.GetXaxis().GetBinLabel(i) for i in range(1, hw.GetNbinsX() + 1)])
fl = ROOT.TFile.Open("/hdd3/przygoda/HADES/FAT_out/all_report/x17ab/x17_lines.root", "READ")
c6 = fl.Get("C6_ee")
d["c6_edges"] = np.array([c6.GetXaxis().GetBinLowEdge(i) for i in range(1, c6.GetNbinsX() + 2)])
d["c6_cont"] = np.array([c6.GetBinContent(i) for i in range(0, c6.GetNbinsX() + 2)])
cs = fl.Get("consts"); d["lines_consts"] = np.array([cs.GetBinContent(i) for i in range(1, cs.GetNbinsX() + 1)])
d["lines_consts_labels"] = np.array([cs.GetXaxis().GetBinLabel(i) for i in range(1, cs.GetNbinsX() + 1)])
hm = fl.Get("c6moments"); d["c6mom"] = np.array([hm.GetBinContent(i) for i in range(1, hm.GetNbinsX() + 1)])
for nm in ["fMean", "fRms", "fEps_nom", "fEps_C2"]:
    fn = fl.Get(nm)
    d[nm] = np.array([fn.Eval(m) for m in np.arange(8.0, 45.01, 0.05)])
d["fgrid"] = np.arange(8.0, 45.01, 0.05)
np.savez(OUT, **d)
print("written", OUT)
print("N_sb masked rows value:", np.unique(d["N_sb"][:, 5:13]))
print("window labels:", list(d["win_labels"]))
print("consts:", list(zip(d["lines_consts_labels"], d["lines_consts"])), "c6mom", d["c6mom"])
