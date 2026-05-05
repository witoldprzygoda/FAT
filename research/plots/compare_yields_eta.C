// research/plots/compare_yields_eta.C — exp vs sim shape comparison for η.
//
// Reads:
//   research/fit_results_eta_<exp_fl>.root        (exp, this branch — default 'rec')
//   research/fit_results_eta_<exp_fl>_sim.root    (sim — copied over from
//                                                  pp45_epem_sim with _sim suffix
//                                                  to avoid filename collision)
//
// Both files come from plots/fit_eta.C — same TTree schema with
// yield_data_{1,2,3}sig{,_err} per OA slice. Both pipelines now use 0.5°
// slicing for η, so OA bin centers match between exp and sim.
//
// Plot 1 — yield_eta_compare_norm_1sig_<fl>.{pdf,png}
//   Normalized η yield (μ ± 1σ) vs opening angle, exp & sim overlaid.
//   Each spectrum scaled so Σ_i y_i = 1 — pure shape comparison.
//
// Plot 2 — yield_eta_compare_ratio_norm_<fl>.{pdf,png}
//   Three series: ratio (norm_exp / norm_sim) for ±1σ, ±2σ, ±3σ windows.
//   Per-bin ratio formula: relative error = √((σ_y_e/y_e)² + (σ_y_s/y_s)²).
//
// Usage (from research/):
//   root -l -b -q plots/compare_yields_eta.C            # exp REC vs sim REC
//   root -l -b -q 'plots/compare_yields_eta.C("cor")'   # exp COR vs sim COR
//   root -l -b -q 'plots/compare_yields_eta.C("tru")'   # exp TRU(?)
//                                                       # — only sensible if
//                                                       #   tru exists on both
//                                                       #   sides; usually
//                                                       #   sim-only

#include <TFile.h>
#include <TTree.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TString.h>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>

namespace {
    struct SliceRow {
        double oa_c, oa_w;
        double y1, e1, y2, e2, y3, e3;
    };

    bool readTree(const std::string& path, std::vector<SliceRow>& rows) {
        TFile* f = TFile::Open(path.c_str(), "READ");
        if (!f || f->IsZombie()) {
            std::cerr << "Cannot open " << path << "\n"; return false;
        }
        auto* t = (TTree*)f->Get("fit_results");
        if (!t) { std::cerr << "No TTree 'fit_results' in " << path << "\n"; return false; }

        int   panel_idx;
        float oa_lo, oa_hi;
        float y1, e1, y2, e2, y3, e3;
        t->SetBranchAddress("panel_idx",           &panel_idx);
        t->SetBranchAddress("oa_lo",               &oa_lo);
        t->SetBranchAddress("oa_hi",               &oa_hi);
        t->SetBranchAddress("yield_data_1sig",     &y1);
        t->SetBranchAddress("yield_data_1sig_err", &e1);
        t->SetBranchAddress("yield_data_2sig",     &y2);
        t->SetBranchAddress("yield_data_2sig_err", &e2);
        t->SetBranchAddress("yield_data_3sig",     &y3);
        t->SetBranchAddress("yield_data_3sig_err", &e3);

        rows.clear();
        for (Long64_t ev = 0; ev < t->GetEntries(); ++ev) {
            t->GetEntry(ev);
            if (panel_idx == 0) continue;
            SliceRow r;
            r.oa_c = 0.5 * (oa_lo + oa_hi);
            r.oa_w = 0.5 * (oa_hi - oa_lo);
            r.y1 = y1;  r.e1 = e1;
            r.y2 = y2;  r.e2 = e2;
            r.y3 = y3;  r.e3 = e3;
            rows.push_back(r);
        }
        f->Close();
        return true;
    }

    double sumCol(const std::vector<SliceRow>& v, int which) {
        double s = 0.0;
        for (const auto& r : v) {
            switch (which) {
                case 1: s += r.y1; break;
                case 2: s += r.y2; break;
                case 3: s += r.y3; break;
            }
        }
        return s;
    }
}

void compare_yields_eta(const char* exp_flavour = "rec") {

    std::string fl = exp_flavour;
    for (auto& c : fl) c = std::tolower(c);
    if (fl != "rec" && fl != "cor" && fl != "tru") {
        std::cerr << "Unknown flavour '" << exp_flavour << "'\n"; return;
    }

    const std::string path_exp = "fit_results_eta_" + fl + ".root";
    const std::string path_sim = "fit_results_eta_" + fl + "_sim.root";

    std::vector<SliceRow> exp_rows, sim_rows;
    if (!readTree(path_exp, exp_rows)) return;
    if (!readTree(path_sim, sim_rows)) return;

    if (exp_rows.size() != sim_rows.size()) {
        std::cerr << "Slice count mismatch: exp=" << exp_rows.size()
                  << " sim=" << sim_rows.size() << "\n";
        return;
    }
    const int N = (int)exp_rows.size();

    for (int i = 0; i < N; ++i) {
        if (std::fabs(exp_rows[i].oa_c - sim_rows[i].oa_c) > 1e-3) {
            std::cerr << "OA bin mismatch at i=" << i
                      << ": exp=" << exp_rows[i].oa_c
                      << " sim=" << sim_rows[i].oa_c << "\n";
            return;
        }
    }

    const double S_e1 = sumCol(exp_rows, 1);
    const double S_e2 = sumCol(exp_rows, 2);
    const double S_e3 = sumCol(exp_rows, 3);
    const double S_s1 = sumCol(sim_rows, 1);
    const double S_s2 = sumCol(sim_rows, 2);
    const double S_s3 = sumCol(sim_rows, 3);

    if (S_e1 <= 0 || S_s1 <= 0) {
        std::cerr << "Empty 1σ totals (exp=" << S_e1 << " sim=" << S_s1 << ")\n";
        return;
    }

    // --- Plot 1: normalized 1σ shape, exp vs sim ----------------------------
    std::vector<double> X(N), EX(N), Ye_n1(N), Ee_n1(N), Ys_n1(N), Es_n1(N);
    for (int i = 0; i < N; ++i) {
        X [i] = exp_rows[i].oa_c;
        EX[i] = exp_rows[i].oa_w;
        Ye_n1[i] = exp_rows[i].y1 / S_e1;
        Ee_n1[i] = exp_rows[i].e1 / S_e1;
        Ys_n1[i] = sim_rows[i].y1 / S_s1;
        Es_n1[i] = sim_rows[i].e1 / S_s1;
    }

    auto* g_exp = new TGraphErrors(N, X.data(), Ye_n1.data(), EX.data(), Ee_n1.data());
    auto* g_sim = new TGraphErrors(N, X.data(), Ys_n1.data(), EX.data(), Es_n1.data());

    g_exp->SetMarkerColor(kBlack); g_exp->SetLineColor(kBlack);
    g_exp->SetMarkerStyle(20);     g_exp->SetMarkerSize(0.9);
    g_sim->SetMarkerColor(kRed);   g_sim->SetLineColor(kRed);
    g_sim->SetMarkerStyle(24);     g_sim->SetMarkerSize(0.9);

    auto* mg1 = new TMultiGraph();
    mg1->Add(g_exp, "P");
    mg1->Add(g_sim, "P");
    mg1->SetTitle(TString::Format(
        "Normalized #eta signal yield (#pm 1#sigma): exp(%s) vs sim;"
        "OA(e^{+}e^{-}) [deg];yield / #Sigma yield",
        fl.c_str()));

    gStyle->SetOptStat(0);
    gSystem->mkdir("plots/output", kTRUE);

    auto* c1 = new TCanvas("c_eta_norm_1sig", "norm η yield 1sig", 1000, 700);
    c1->SetMargin(0.13, 0.05, 0.12, 0.08);
    c1->SetGrid();
    mg1->Draw("A");
    mg1->GetXaxis()->SetLimits(0.0, 15.0);

    auto* leg1 = new TLegend(0.62, 0.78, 0.94, 0.90);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    leg1->SetTextSize(0.035);
    leg1->AddEntry(g_exp, TString::Format("exp (%s)", fl.c_str()).Data(), "lpe");
    leg1->AddEntry(g_sim, "sim",                                          "lpe");
    leg1->Draw();

    c1->Update();
    const std::string out1 = "plots/output/yield_eta_compare_norm_1sig_" + fl;
    c1->SaveAs((out1 + ".pdf").c_str());
    c1->SaveAs((out1 + ".png").c_str());

    // --- Plot 2: normalized exp/sim ratio for 1σ, 2σ, 3σ -------------------
    auto buildRatio = [&](int Nsig, double S_e, double S_s,
                          std::vector<double>& R, std::vector<double>& ER) {
        R.resize(N);  ER.resize(N);
        for (int i = 0; i < N; ++i) {
            const double y_e = (Nsig == 1) ? exp_rows[i].y1
                              : (Nsig == 2) ? exp_rows[i].y2 : exp_rows[i].y3;
            const double e_e = (Nsig == 1) ? exp_rows[i].e1
                              : (Nsig == 2) ? exp_rows[i].e2 : exp_rows[i].e3;
            const double y_s = (Nsig == 1) ? sim_rows[i].y1
                              : (Nsig == 2) ? sim_rows[i].y2 : sim_rows[i].y3;
            const double e_s = (Nsig == 1) ? sim_rows[i].e1
                              : (Nsig == 2) ? sim_rows[i].e2 : sim_rows[i].e3;

            if (y_e <= 0 || y_s <= 0) { R[i] = 0; ER[i] = 0; continue; }
            const double r = (y_e / S_e) / (y_s / S_s);
            const double rel = std::sqrt((e_e/y_e)*(e_e/y_e) + (e_s/y_s)*(e_s/y_s));
            R [i] = r;
            ER[i] = r * rel;
        }
    };

    std::vector<double> R1, ER1, R2, ER2, R3, ER3;
    buildRatio(1, S_e1, S_s1, R1, ER1);
    buildRatio(2, S_e2, S_s2, R2, ER2);
    buildRatio(3, S_e3, S_s3, R3, ER3);

    auto* gr1 = new TGraphErrors(N, X.data(), R1.data(), EX.data(), ER1.data());
    auto* gr2 = new TGraphErrors(N, X.data(), R2.data(), EX.data(), ER2.data());
    auto* gr3 = new TGraphErrors(N, X.data(), R3.data(), EX.data(), ER3.data());

    gr1->SetMarkerColor(kBlue);      gr1->SetLineColor(kBlue);
    gr1->SetMarkerStyle(20);         gr1->SetMarkerSize(0.9);
    gr2->SetMarkerColor(kGreen + 2); gr2->SetLineColor(kGreen + 2);
    gr2->SetMarkerStyle(21);         gr2->SetMarkerSize(0.9);
    gr3->SetMarkerColor(kRed);       gr3->SetLineColor(kRed);
    gr3->SetMarkerStyle(22);         gr3->SetMarkerSize(1.0);

    auto* mg2 = new TMultiGraph();
    mg2->Add(gr1, "P");
    mg2->Add(gr2, "P");
    mg2->Add(gr3, "P");
    mg2->SetTitle(TString::Format(
        "Normalized #eta yield ratio: exp(%s) / sim;"
        "OA(e^{+}e^{-}) [deg];"
        "(y_{exp}/#Sigma_{exp}) / (y_{sim}/#Sigma_{sim})",
        fl.c_str()));

    auto* c2 = new TCanvas("c_eta_ratio", "exp/sim η ratio", 1000, 700);
    c2->SetMargin(0.13, 0.05, 0.12, 0.08);
    c2->SetGrid();
    mg2->Draw("A");
    mg2->GetXaxis()->SetLimits(0.0, 15.0);

    auto* lref = new TLine(0.0, 1.0, 15.0, 1.0);
    lref->SetLineColor(kGray + 2);
    lref->SetLineStyle(2);
    lref->SetLineWidth(2);
    lref->Draw();

    auto* leg2 = new TLegend(0.62, 0.72, 0.94, 0.90);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetTextSize(0.035);
    leg2->AddEntry(gr1,  "ratio (#pm 1#sigma)", "lpe");
    leg2->AddEntry(gr2,  "ratio (#pm 2#sigma)", "lpe");
    leg2->AddEntry(gr3,  "ratio (#pm 3#sigma)", "lpe");
    leg2->AddEntry(lref, "shape match (= 1)",   "l");
    leg2->Draw();

    c2->Update();
    const std::string out2 = "plots/output/yield_eta_compare_ratio_norm_" + fl;
    c2->SaveAs((out2 + ".pdf").c_str());
    c2->SaveAs((out2 + ".png").c_str());

    std::cout << "Wrote:\n"
              << "  " << out1 << ".{pdf,png}\n"
              << "  " << out2 << ".{pdf,png}\n"
              << "  (" << N << " OA slices)\n"
              << "  exp totals: S_1σ=" << S_e1 << " S_2σ=" << S_e2 << " S_3σ=" << S_e3 << "\n"
              << "  sim totals: S_1σ=" << S_s1 << " S_2σ=" << S_s2 << " S_3σ=" << S_s3 << "\n";
}
