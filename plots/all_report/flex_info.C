// ============================================================================
// flex_info.C  --  package flex of report_all: what constrains each scale of the flexible fit, and what each scale changes
// ============================================================================
// For every free scale s_t of a fit (flex_fit.root directory), the Fisher information of the bins, I_i = T_i^2 / V_i (T = template
// = fitted part / s, V = the fit variance incl. BB-lite), is split by M4 region (< 0.455, bump B 0.455-0.5225, 0.5225-0.5336,
// peak P 0.5336-0.5542, 0.5542-0.565, CR-U 0.565-0.70) and by M(pi+pi-) slice: the shares say which part of the data fixes the
// scale (before correlations with the other scales). Second block: rows of every component in the windows at SMASH and after
// the fit (from flex_results.txt "win" lines, reformatted).
// Run: root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_flex", kTRUE)' 'plots/all_report/flex_info.C+O("F8")'
// Output: docs/latex/tables_all/flex_info.txt
// ============================================================================
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1D.h"

void flex_info(const char* tag = "F8") {
    TFile* f = TFile::Open("/hdd3/przygoda/HADES/FAT_out/all_report/flex/flex_fit.root");
    TDirectory* d = (TDirectory*)f->Get(tag);
    TH1D* th = (TH1D*)d->Get("theta");
    FILE* o = fopen("/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all/flex_info.txt", "w");
    fprintf(o, "# flex_info.txt -- plots/all_report/flex_info.C: what constrains each scale of the flexible fit %s (F10 E2, M4 0.40-0.70)\n", tag);
    fprintf(o, "# share (%%) of the Fisher information sum_i T_i^2 / V_i of each scale by M4 region and by M(pi+pi-) slice (T = template, V = fit variance)\n");
    const int nR = 6;
    const double lo[nR] = {0.40, 0.4550, 0.5225, 0.53364, 0.55416, 0.565}, hi[nR] = {0.4550, 0.5225, 0.53364, 0.55416, 0.565, 0.70};
    const char* rn[nR] = {"<0.455", "B", "0.5225-0.534", "P", "0.554-0.565", "CR-U"};
    const char* comps[] = {"C2", "C3", "C4", "C5", "Geta", "GD", "GD_om", "GD_nom", "GC", "GC_om", "GC_nom"};
    fprintf(o, "%-8s %8s |", "scale", "value");
    for (int r = 0; r < nR; ++r) fprintf(o, " %12s", rn[r]);
    fprintf(o, " |");
    for (int j = 1; j <= 5; ++j) fprintf(o, "     S%d", j);
    fprintf(o, "\n");
    for (const char* c : comps) {
        int k = -1;
        for (int b = 1; b <= th->GetNbinsX(); ++b) if (std::string(th->GetXaxis()->GetBinLabel(b)) == std::string("s_") + c) k = b;
        if (k < 0) continue;
        const double s = th->GetBinContent(k);
        double reg[nR] = {0}, sl[5] = {0}, tot = 0;
        for (int j = 1; j <= 5; ++j) {
            TH1D* p = (TH1D*)d->Get(Form("%s_S%dE2", c, j));
            TH1D* mu = (TH1D*)d->Get(Form("mu_S%dE2", j));
            if (!p || !mu) continue;
            for (int b = 1; b <= p->GetNbinsX(); ++b) {
                const double V = std::pow(mu->GetBinError(b), 2);
                if (!(V > 0) || s == 0) continue;
                const double T = p->GetBinContent(b) / s;
                const double I = T * T / V;
                const double x = p->GetBinCenter(b);
                for (int r = 0; r < nR; ++r) if (x >= lo[r] && x < hi[r]) reg[r] += I;
                sl[j - 1] += I; tot += I;
            }
        }
        fprintf(o, "s_%-6s %8.3f |", c, s);
        for (int r = 0; r < nR; ++r) fprintf(o, " %11.1f%%", tot > 0 ? 100 * reg[r] / tot : 0.);
        fprintf(o, " |");
        for (int j = 0; j < 5; ++j) fprintf(o, " %5.1f%%", tot > 0 ? 100 * sl[j] / tot : 0.);
        fprintf(o, "\n");
    }
    fclose(o);
    printf("flex_info: written\n");
}
