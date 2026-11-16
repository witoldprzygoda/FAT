// fix_dump_stageA.C -- review fix (prefix fix): dumps the frozen stage-A histograms of wp3/stageA.root (read only) to text
// for the redrawn Figs. C2/C3 (fix_after.py). Output: /hdd3/przygoda/HADES/FAT_out/all_report/fix/stageA_dump.txt
// Run: root -l -b -q plots/all_report/fix_dump_stageA.C
#include "TFile.h"
#include "TH1D.h"
#include <cstdio>
void fix_dump_stageA() {
    TFile* f = TFile::Open("/hdd3/przygoda/HADES/FAT_out/all_report/wp3/stageA.root", "READ");
    FILE* o = fopen("/hdd3/przygoda/HADES/FAT_out/all_report/fix/stageA_dump.txt", "w");
    const char* kinds[5] = {"GDsmooth", "GCsmooth", "BA", "BAstat", "BAmc"};
    for (int j = 1; j <= 5; ++j) for (int e = 1; e <= 2; ++e) for (int k = 0; k < 5; ++k) {
        TH1D* h = (TH1D*)f->Get(Form("%s_F10_S%dE%d", kinds[k], j, e));
        if (!h) continue;
        fprintf(o, "%s_F10_S%dE%d", kinds[k], j, e);
        for (int b = 1; b <= h->GetNbinsX(); ++b) fprintf(o, " %.8g", (k >= 3) ? h->GetBinError(b) : h->GetBinContent(b));
        fprintf(o, "\n");
        if (j == 1 && e == 1 && k == 0) { printf("bins %d: ", h->GetNbinsX()); for (int b = 1; b <= h->GetNbinsX() + 1; ++b) printf("%.4f ", h->GetBinLowEdge(b)); printf("\n"); }
    }
    fclose(o); f->Close();
}
