// ============================================================================
// x17ab_parent.C  --  WP9 X-b: the common offset of the simulation-study X17 toys (parent-MC structure)
// ============================================================================
// The 1500 toys of the simulation study are Poisson fluctuations of ONE parent (half B of the simulation, odd event
// groups, x 2 x 27.5 at L = 6080.5 nb^-1). One std-set row of the parent is about 55 expected rows in a single toy cell, so
// the parent carries spikes and 1 %-level shape structure that every toy shares. The model templates come from the other
// half (A). This macro fits, with the nominal toy-grid model (x17ab_toys.C cfg 0):
//   (a) the mean of the 300 toys at BR = 0 (a pseudo-Asimov set of the half-B parent; no CB),
//   (b) the half-A Asimov set (sumT12 of the templates x 27.5; its own parent),
// and reports N5, r_hat (BR1 = r x 2.68e-4, BR2 = r N5_hat/N5_true x 2.68e-4) and the Pearson chi2/ndf. The difference of (a)
// and (b) measures the size of the parent-MC structure projected on the X17 estimate.
// Output: docs/latex/tables_all/x17ab_parent.txt
// Run: root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_x17ab", kTRUE)' 'plots/all_report/x17ab_parent.C+O'
// ============================================================================
#include <cmath>
#include <cstdio>
#include <vector>

#include "TFile.h"
#include "TH2D.h"
#include "TTree.h"

#include "defs.h"
#include "x17_model.h"

void x17ab_parent() {
    x17m::ToyGridOpt to;
    x17m::Inputs I;
    if (!x17m::BuildToyGrid(I, to)) return;
    x17m::Config C;
    C.kSig[0] = 0.; C.kSig[1] = 0.;   // no CB in these sets
    C.a0 = 1.; C.aSig = 0.10; C.gSig = 0.30; C.bkgMode = 4;
    {
        TFile* f = TFile::Open(x17m::kToyTemplates, "READ");
        std::vector<double> e2;
        std::vector<double> v = x17m::ToyRebin((TH2*)f->Get("m4mee_X17_conv4"), to, &e2);
        double s = 0, ee = 0;
        for (int j = 0; j < I.G.ny; ++j) if (I.Ctail[j] > 0) for (int i = 0; i < I.G.nx; ++i) { s += v[I.G.idx(i, j)]; ee += e2[I.G.idx(i, j)]; }
        C.cTailSig = s > 0 ? std::sqrt(ee) / s : 0.;
        f->Close(); delete f;
    }
    // (a) mean of the BR = 0 toys
    std::vector<double> avg(I.G.n(), 0.);
    int nav = 0;
    {
        TFile* fk = TFile::Open(x17m::kToyKey, "READ");
        TTree* tk = (TTree*)fk->Get("key");
        int id; double br;
        tk->SetBranchAddress("id", &id); tk->SetBranchAddress("br", &br);
        std::vector<int> ids;
        for (Long64_t i = 0; i < tk->GetEntries(); ++i) { tk->GetEntry(i); if (br == 0) ids.push_back(id); }
        fk->Close(); delete fk;
        TFile* fx = TFile::Open(x17m::kToyFile, "READ");
        for (int id0 : ids) {
            TH2* h = (TH2*)fx->Get(Form("toy_%d_m4mee", id0));
            std::vector<double> n = x17m::ToyRebin(h, to);
            for (size_t b = 0; b < n.size(); ++b) avg[b] += n[b];
            ++nav; delete h;
        }
        fx->Close(); delete fx;
        for (auto& v : avg) v /= nav;
    }
    // (b) half-A Asimov
    std::vector<double> asA;
    {
        TFile* f = TFile::Open(x17m::kToyTemplates, "READ");
        asA = x17m::ToyRebin((TH2*)f->Get("m4mee_X17_sumT12"), to);
        for (auto& v : asA) v *= to.scale;
        f->Close(); delete f;
    }
    FILE* t = fopen("/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all/x17ab_parent.txt", "w");
    fprintf(t, "# x17ab_parent.txt -- parent-MC structure of the simulation-study X17 toys (x17ab_parent.C), nominal toy-grid model (bkgMode 4,\n");
    fprintf(t, "# B = exp(pol2)-smoothed half A), no CB. N5_true = %.1f rows (half A). BR in 1e-6. Pearson chi2/ndf with the per-toy Poisson variance\n", I.N5start);
    fprintf(t, "# (for the mean of %d toys the Poisson part is 1/%d of it: chi2/ndf x %d is the structure-to-Poisson ratio).\n", nav, nav, nav);
    fprintf(t, "%-34s %8s %9s %8s %8s %10s\n", "data set", "N5_hat", "BR1", "BR2", "err", "chi2/ndf");
    for (int k = 0; k < 2; ++k) {
        const std::vector<double>& n = k == 0 ? avg : asA;
        x17m::Model M(I, C);
        x17m::Profile P = x17m::ProfileFits(M, n, x17m::GlobOf(C));
        double chi = 0; int nd = 0;
        for (int j = 0; j < I.G.ny; ++j) for (int i = 0; i < I.G.nx; ++i) {
            const double mu = M.Mu(P.hat.p, i, j);
            if (mu > 0) { chi += std::pow(n[I.G.idx(i, j)] - mu, 2) / mu; ++nd; }
        }
        int nf = 0; for (char c : M.FreeMask(true)) nf += c;
        nd -= nf;
        const double N5 = P.hat.p[M.L.iN5];
        const char* nm = k == 0 ? "mean of the BR = 0 toys (half B)" : "half-A Asimov set";
        fprintf(t, "%-34s %8.0f %9.2f %8.2f %8.2f %10.2f\n", nm, N5, P.rhat * x17m::kBRdir * 1e6, P.rhat * N5 / I.N5start * x17m::kBRdir * 1e6,
                P.rhatErr * x17m::kBRdir * 1e6, chi / nd);
        printf("x17parent: %s: N5 %.0f (true %.0f), BR1 %.2fe-6, BR2 %.2fe-6, err %.2fe-6, chi2/ndf %.2f (%d)\n", nm, N5, I.N5start, P.rhat * x17m::kBRdir * 1e6,
               P.rhat * N5 / I.N5start * x17m::kBRdir * 1e6, P.rhatErr * x17m::kBRdir * 1e6, chi / nd, nd);
    }
    fclose(t);
}
