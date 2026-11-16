// vflex_g3.C -- independent verification of package flex: SMASH as is vs data above M4 = 0.70 from the G3 layout of
// hists_v1.root (M4 0.40-0.85 in 2.5 MeV x M(pi+pi-) x M(ee) classes), M(ee) >= 10 MeV, all slices; (N - CB) / SMASH per
// 25 MeV band and the share of each SMASH component (independent of the flex X-layout files).
// Run: root -l -b -q plots/all_report/vflex_g3.C
#include <cstdio>
#include "TFile.h"
#include "TH3.h"
void vflex_g3() {
    TFile* f = TFile::Open("/hdd3/przygoda/HADES/FAT_out/all_report/hists_v1.root");
    TH3* N = (TH3*)f->Get("data/S0/all/N_G3");
    TH3* C = (TH3*)f->Get("data/S0/all/CB_G3");
    const char* comps[] = {"C2", "C3", "C4", "C5", "Geta", "GD", "Gconv"};
    TH3* S[7];
    for (int q = 0; q < 7; ++q) S[q] = (TH3*)f->Get(Form("sim/S0__nom/%s_G3", comps[q]));
    const double nb2raw = 6079.2 * 63.0 / 64.0;
    const int ny = N->GetNbinsY(), nz = N->GetNbinsZ();
    printf("band        N        CB      SMASH  (N-CB)/SMASH | shares C2 C3 C4 C5 Geta GD GC (%%)\n");
    for (int b = 0; b < 18; ++b) {
        const int x0 = 1 + 10 * b, x1 = x0 + 9;
        double n = 0, c = 0, s[7] = {0}, st = 0;
        for (int x = x0; x <= x1; ++x)
            for (int y = 1; y <= ny; ++y)
                for (int z = 3; z <= nz; ++z) {
                    n += N->GetBinContent(x, y, z); c += C->GetBinContent(x, y, z);
                    for (int q = 0; q < 7; ++q) s[q] += S[q]->GetBinContent(x, y, z) * nb2raw;
                }
        for (int q = 0; q < 7; ++q) st += s[q];
        printf("%.3f-%.3f %8.0f %8.0f %9.0f %7.3f |", N->GetXaxis()->GetBinLowEdge(x0), N->GetXaxis()->GetBinUpEdge(x1), n, c, st, (n - c) / st);
        for (int q = 0; q < 7; ++q) printf(" %4.1f", 100 * s[q] / st);
        printf("\n");
    }
    printf("M(ee) edges of G3 z bins:"); for (int z = 1; z <= nz + 1; ++z) printf(" %.3f", N->GetZaxis()->GetBinLowEdge(z)); printf("\n");
}
