// ============================================================================
// vflex_dump.C  --  independent verification of package flex (prefix vflex): plain-text dump of the Part III inputs
// ============================================================================
// What it does: reads hists_v1.root (WP0) and wp5_egcov.root DIRECTLY (not the flex X-layout files, not the flex library)
// and writes plain text for the independent Python fitter vflex_fit.py:
//   data_F10.txt   : j e bin N Nerr CB            (S0, all, F10 slices 1-5, E1/E2, 70 bins)
//   cbrep_F10.txt  : j e bin r_0 .. r_199         (like-sign CB replicas)
//   egcov_F10.txt  : 700 x 700 event-group covariance (cell = (2 j + e) * 70 + bin - 1)
//   sim_F10.txt    : comp j e bin content error   (sim/S0__nom, nb per bin)
//   simrep_F10.txt : comp j e bin r_0 .. r_99     (WP0 bootstrap replicas, nb per bin)
// Run (from /home/przygoda/HADES/FAT/FAT):  root -l -b -q plots/all_report/vflex_dump.C
// ============================================================================
#include <cstdio>
#include <string>
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMatrixDSym.h"
#include "TSystem.h"

void vflex_dump() {
    const char* out = "/hdd3/przygoda/HADES/FAT_out/all_report/vflex";
    gSystem->mkdir(out, kTRUE);
    TFile* f = TFile::Open("/hdd3/przygoda/HADES/FAT_out/all_report/hists_v1.root");
    TFile* g = TFile::Open("/hdd3/przygoda/HADES/FAT_out/all_report/wp5/wp5_egcov.root");
    if (!f || !g) { printf("vflex_dump: input missing\n"); return; }
    FILE* fd = fopen(Form("%s/data_F10.txt", out), "w");
    FILE* fr = fopen(Form("%s/cbrep_F10.txt", out), "w");
    for (int j = 0; j < 5; ++j)
        for (int e = 0; e < 2; ++e) {
            const std::string h = Form("F10_S%dE%d", j + 1, e + 1);
            TH1* n = (TH1*)f->Get(("data/S0/all/N_" + h).c_str());
            TH1* c = (TH1*)f->Get(("data/S0/all/CB_" + h).c_str());
            TH2* r = (TH2*)f->Get(("data/S0/all/CBrep_" + h).c_str());
            if (!n || !c || !r) { printf("vflex_dump: data %s missing\n", h.c_str()); return; }
            if (j == 0 && e == 0) printf("F10 bins %d, edges %.4f-%.4f, CB replicas %d\n", n->GetNbinsX(), n->GetXaxis()->GetXmin(), n->GetXaxis()->GetXmax(), r->GetNbinsY());
            for (int x = 1; x <= n->GetNbinsX(); ++x) {
                fprintf(fd, "%d %d %d %.17g %.17g %.17g %.9g %.9g\n", j, e, x, n->GetBinContent(x), n->GetBinError(x), c->GetBinContent(x),
                        n->GetXaxis()->GetBinLowEdge(x), n->GetXaxis()->GetBinUpEdge(x));
                fprintf(fr, "%d %d %d", j, e, x);
                for (int q = 1; q <= r->GetNbinsY(); ++q) fprintf(fr, " %.17g", r->GetBinContent(x, q));
                fprintf(fr, "\n");
            }
        }
    fclose(fd); fclose(fr);
    TMatrixDSym* M = (TMatrixDSym*)g->Get("egcov_S0_all_F10");
    if (!M) { printf("vflex_dump: egcov missing\n"); return; }
    FILE* fe = fopen(Form("%s/egcov_F10.txt", out), "w");
    for (int a = 0; a < M->GetNrows(); ++a) { for (int b = 0; b < M->GetNcols(); ++b) fprintf(fe, "%.17g ", (*M)(a, b)); fprintf(fe, "\n"); }
    fclose(fe);
    const char* comps[] = {"C2", "C3", "C4", "C5", "Geta", "GD", "Gconv", "Gpi0D", "Gvec", "Goth", "Gpi0D_om", "Gvec_om", "Goth_om",
                           "Gpi0D_nom", "Gvec_nom", "Goth_nom", "Gconv_om", "Gconv_nom", "T3"};
    FILE* fs = fopen(Form("%s/sim_F10.txt", out), "w");
    FILE* fq = fopen(Form("%s/simrep_F10.txt", out), "w");
    for (const char* cp : comps)
        for (int j = 0; j < 5; ++j)
            for (int e = 0; e < 2; ++e) {
                const std::string h = Form("sim/S0__nom/%s_F10_S%dE%d", cp, j + 1, e + 1);
                TH1* s = (TH1*)f->Get(h.c_str());
                TH2* r = (TH2*)f->Get((h + "_rep").c_str());
                if (!s || !r) { printf("vflex_dump: %s missing\n", h.c_str()); return; }
                for (int x = 1; x <= s->GetNbinsX(); ++x) {
                    fprintf(fs, "%s %d %d %d %.17g %.17g\n", cp, j, e, x, s->GetBinContent(x), s->GetBinError(x));
                    fprintf(fq, "%s %d %d %d", cp, j, e, x);
                    for (int q = 1; q <= r->GetNbinsY(); ++q) fprintf(fq, " %.17g", r->GetBinContent(x, q));
                    fprintf(fq, "\n");
                }
            }
    fclose(fs); fclose(fq);
    printf("vflex_dump: written to %s\n", out);
}
