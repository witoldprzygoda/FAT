// ============================================================================
// x17c_toys.C  --  WP9 X-c (prefix x17c): background-only (and injection) toys of the FROZEN X17 procedure
// ============================================================================
// SPEC report_all WP9 X-c. Written 2026-11-14.
// Truth: the background-only expectation of the data side-band fit (r = 0) with the masked rows 12-28 MeV/c2 predicted
// (x17cl::FitSideBands, nominal frozen configuration of x17c_config.txt), optionally plus an X17 line at mass mX with
// BR_X BR_ee = inj (normalisation (1): rows = inj/2.68e-4 x N5 R_X).
// Pseudo-data: x17c::ToyData = qScale x Poisson(mu / qScale) per cell (over-dispersed by the variance factor), and the
// centres of all Gaussian constraints drawn around the truth (global observables).
// Fit: every toy is fitted with the full model on ALL cells (as X-d will fit the data) in the three normalisations; the
// 90 % CLs upper limit is the asymptotic one (q~_r / qScale, Asimov from the toy's own r = 0 fit).
// Output: /hdd3/przygoda/HADES/FAT_out/all_report/x17c/toys_m<100 mX>_i<inj in 1e-7>_<first>.root (tree "toys")
// Run: root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_x17c", kTRUE)' \
//        'plots/all_report/x17c_toys.C+O(17.0, 0., 0, 500)'
// ============================================================================
#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

#include "TFile.h"
#include "TStopwatch.h"
#include "TTree.h"

#include "x17c_lib.h"

void x17c_toys(double mX = 17.0, double inj = 0., int first = 0, int ntoys = 500, unsigned long long seed = 20261114ULL) {
    using namespace x17cl;
    TStopwatch sw;
    TFile* f = TFile::Open(alldefs::kHistFile, "READ");
    Wp5 w = ReadWp5();
    Frozen z;
    if (!w.ok || !ReadConfig(z)) return;
    x17c::Inputs I;
    if (!Build(f, I, w, BuildOpt(), z.a0)) return;
    if (!SetLine(I, mX)) return;
    std::vector<double> n;
    std::vector<char> use;
    if (!LoadSideBands(f, n, use)) return;
    const std::vector<double> simRow = SimRow(f, I);          // nominal reference of the masked-row prediction
    const std::vector<double> simRel = SimRowRelErr(f, I);
    x17c::Config C1 = MakeConfig(z, 1);
    x17c::Config Csb = C1; Csb.use = use;
    SideBand SB = FitSideBands(I, Csb, n, simRow, z.fillOrder, &simRel);
    x17c::Model M1(I, C1);
    std::vector<double> ptrue = SB.pFill;
    ptrue[M1.L.iR] = inj / kBR;
    const x17c::Glob gobs = x17c::GlobOf(C1);
    x17c::Config Cn[3] = {MakeConfig(z, 1), MakeConfig(z, 2), MakeConfig(z, 3)};
    printf("x17c_toys: m %.2f, inj %.2e, toys %d-%d, qScale %.4f, sX %.3f/%.3f/%.3f; side-band fit ok %d, N5 %.1f; truth X17 rows %.1f\n", mX, inj, first,
           first + ntoys - 1, z.qScale, z.xSig1, z.xSig2, z.xSig3, SB.fit.ok, SB.fit.p[M1.L.iN5], inj / kBR * SB.fit.p[M1.L.iN5] * I.RX);
    const std::string out = Form("%s/toys_m%04d_i%04d_%04d.root", kOutDir, (int)std::lround(100 * mX), (int)std::lround(inj * 1e7), first);
    TFile* fo = TFile::Open(out.c_str(), "RECREATE");
    TTree t("toys", "X17 toys of the frozen procedure (BR units)");
    int id, ok[3];
    double ul[3], rhat[3], rerr[3], N5h, mXo = mX, injo = inj, nrows;
    t.Branch("id", &id); t.Branch("mX", &mXo); t.Branch("inj", &injo);
    t.Branch("ok", ok, "ok[3]/I"); t.Branch("ul", ul, "ul[3]/D"); t.Branch("rhat", rhat, "rhat[3]/D"); t.Branch("rerr", rerr, "rerr[3]/D");
    t.Branch("N5", &N5h); t.Branch("nrows", &nrows);
    std::vector<double> nt;
    x17c::Glob gt;
    for (int k = first; k < first + ntoys; ++k) {
        std::mt19937_64 rng(seed ^ (0x9e3779b97f4a7c15ULL * (unsigned long long)(k + 1)));   // toy k is the same data set for every mX
        x17c::ToyData(M1, ptrue, rng, nt, gt, gobs);
        id = k;
        nrows = 0; for (double v : nt) nrows += v;
        for (int a = 0; a < 3; ++a) {
            Lim l = LimitBR(I, Cn[a], nt, false, 0.1, &gt);
            // global observables of the toy: the same draws for all normalisations (constraints are identical except sX)
            ok[a] = l.ok; ul[a] = l.ul; rhat[a] = l.rhat; rerr[a] = l.rerr;
            if (a == 0) N5h = l.N5;
        }
        t.Fill();
        if ((k - first) % 100 == 0)
            printf("x17c_toys: toy %d: UL %.3e / %.3e / %.3e, rhat %.3e +- %.3e, N5 %.0f, %.0f s\n", k, ul[0], ul[1], ul[2], rhat[0], rerr[0], N5h, sw.RealTime()),
                sw.Continue();
    }
    t.Write();
    fo->Close();
    f->Close();
    printf("x17c_toys: done, %s, %.0f s\n", out.c_str(), sw.RealTime());
}
