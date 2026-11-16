// ============================================================================
// v5_dump.C  --  independent verification of WP5 (prefix v5): dump the inputs of the F-main / F-gold fits
// ============================================================================
// Reads ONLY hists_v1.root (WP0) and the two event-group covariance files (WP3 wp3_rows.root, WP5 wp5_egcov.root)
// and writes plain-text tables to /hdd3/przygoda/HADES/FAT_out/all_report/v5/ for the independent Python fitter
// v5_fit.py. No code of wp4_lib.h / wp5_lib.h is used. Nothing is written anywhere else.
// Run (from /home/przygoda/HADES/FAT/FAT):
//   root -l -b -q 'plots/all_report/v5_dump.C'
// ============================================================================
#include <cstdio>
#include <string>
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMatrixDSym.h"

static void Row(FILE* f, const std::string& tag, TH1* h, bool err) {
    fprintf(f, "%s", tag.c_str());
    for (int b = 1; b <= h->GetNbinsX(); ++b) fprintf(f, " %.12g", err ? h->GetBinError(b) : h->GetBinContent(b));
    fprintf(f, "\n");
}
static void Rep(FILE* f, const std::string& tag, TH2* h) {
    for (int r = 1; r <= h->GetNbinsY(); ++r) {
        fprintf(f, "%s %d", tag.c_str(), r - 1);
        for (int b = 1; b <= h->GetNbinsX(); ++b) fprintf(f, " %.12g", h->GetBinContent(b, r));
        fprintf(f, "\n");
    }
}
static void Mat(const char* fn, TMatrixDSym* m) {
    FILE* f = fopen(fn, "w");
    for (int a = 0; a < m->GetNrows(); ++a) { for (int b = 0; b < m->GetNcols(); ++b) fprintf(f, "%.12g ", (*m)(a, b)); fprintf(f, "\n"); }
    fclose(f);
}

void v5_dump() {
    const std::string out = "/hdd3/przygoda/HADES/FAT_out/all_report/v5/";
    TFile* fh = TFile::Open("/hdd3/przygoda/HADES/FAT_out/all_report/hists_v1.root");
    std::vector<std::string> hn;
    for (int j = 1; j <= 5; ++j) for (int e = 1; e <= 2; ++e) hn.push_back(Form("F10_S%dE%d", j, e));
    for (int j = 1; j <= 5; ++j) hn.push_back(Form("F5g_S%d", j));
    // ---- data
    for (const char* sel : {"S0", "pre", "S0_best"}) {
        FILE* f = fopen((out + "data_" + sel + ".txt").c_str(), "w");
        FILE* fr = fopen((out + "cbrep_" + sel + ".txt").c_str(), "w");
        for (const auto& h : hn) {
            for (const char* w : {"N", "CB", "Npp", "Nmm"}) {
                TH1* x = (TH1*)fh->Get(Form("data/%s/all/%s_%s", sel, w, h.c_str()));
                if (!x) { printf("missing data/%s/all/%s_%s\n", sel, w, h.c_str()); continue; }
                Row(f, std::string(w) + "_" + h, x, false);
                if (std::string(w) == "N") Row(f, "Nerr_" + h, x, true);
            }
            TH2* r = (TH2*)fh->Get(Form("data/%s/all/CBrep_%s", sel, h.c_str()));
            if (r) Rep(fr, "CBrep_" + h, r); else printf("missing CBrep %s %s\n", sel, h.c_str());
        }
        for (const char* w : {"N", "CB"}) for (const char* l : {"win", "tot"}) {
            TH1* x = (TH1*)fh->Get(Form("data/%s/all/%s_%s", sel, w, l));
            if (x) Row(f, std::string(w) + "_" + l, x, false);
        }
        fclose(f); fclose(fr);
    }
    // ---- simulation
    const char* comps[] = {"C2", "C3", "C4", "C5", "Geta", "GD", "Gconv", "Gpi0D", "Gvec", "Goth", "T12"};
    for (const char* sel : {"S0", "pre"}) {
        FILE* f = fopen((out + Form("sim_%s.txt", sel)).c_str(), "w");
        FILE* fr = fopen((out + Form("simrep_%s.txt", sel)).c_str(), "w");
        for (const char* c : comps) {
            for (const auto& h : hn) {
                TH1* x = (TH1*)fh->Get(Form("sim/%s__nom/%s_%s", sel, c, h.c_str()));
                if (!x) { printf("missing sim %s %s %s\n", sel, c, h.c_str()); continue; }
                Row(f, std::string(c) + "_" + h, x, false);
                Row(f, std::string(c) + "err_" + h, x, true);
                TH2* r = (TH2*)fh->Get(Form("sim/%s__nom/%s_%s_rep", sel, c, h.c_str()));
                if (r) Rep(fr, std::string(c) + "_" + h, r);
            }
            for (const char* l : {"win", "tot"}) {
                TH1* x = (TH1*)fh->Get(Form("sim/%s__nom/%s_%s", sel, c, l));
                if (x) { Row(f, std::string(c) + "_" + l, x, false); Row(f, std::string(c) + "err_" + l, x, true); }
                TH2* r = (TH2*)fh->Get(Form("sim/%s__nom/%s_%s_rep", sel, c, l));
                if (r) Rep(fr, std::string(c) + "_" + l, r);
            }
        }
        fclose(f); fclose(fr);
    }
    // ---- s / r variants of the simulation (F10 only, signal + non-eta components)
    {
        FILE* f = fopen((out + "sim_S0_svar.txt").c_str(), "w");
        for (const char* v : {"sK0", "sUp", "sDn", "r089", "r110"})
            for (const char* c : {"C2", "C3", "C4", "C5", "Geta", "GD", "Gconv"})
                for (const auto& h : hn) {
                    TH1* x = (TH1*)fh->Get(Form("sim/S0__%s/%s_%s", v, c, h.c_str()));
                    if (!x) { printf("missing sim S0__%s %s %s\n", v, c, h.c_str()); continue; }
                    Row(f, std::string(v) + ":" + c + "_" + h, x, false);
                    Row(f, std::string(v) + ":" + c + "err_" + h, x, true);
                }
        fclose(f);
    }
    // ---- event-group covariances
    TFile* f3 = TFile::Open("/hdd3/przygoda/HADES/FAT_out/all_report/wp3/wp3_rows.root");
    Mat((out + "egcov_wp3_S0_F10.txt").c_str(), (TMatrixDSym*)f3->Get("egcov_S0"));
    TFile* f5 = TFile::Open("/hdd3/przygoda/HADES/FAT_out/all_report/wp5/wp5_egcov.root");
    Mat((out + "egcov_wp5_S0_F10.txt").c_str(), (TMatrixDSym*)f5->Get("egcov_S0_all_F10"));
    Mat((out + "egcov_wp5_S0_F5g.txt").c_str(), (TMatrixDSym*)f5->Get("egcov_S0_all_F5g"));
    printf("v5_dump done\n");
}
