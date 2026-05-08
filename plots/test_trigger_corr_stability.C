// Verify PT3 trigger-bias correction stability over time/event count.
//
// Reads the input chain (single .root or .list of files), counts trigbits per
// chunk (chunk_size events) and cumulatively, and reports
//
//   w = (63 · N_PT2) / N_PT3
//
// per chunk vs. the running cumulative average. A stable run gives a flat
// curve; visible drift between chunks would mean the correction needs to be
// applied per-time-window, not per-run.
//
// Usage (from FAT/):
//   root -l -b -q 'plots/test_trigger_corr_stability.C(\
//       "/hdd2/przygoda/hades/pp45/exp/GEN4/LEPTONS/pp45_leptons_tot.list",\
//       "EpEm_ID", 1000000)'

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TString.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>

void test_trigger_corr_stability(
    const char* input_source = "/hdd2/przygoda/hades/pp45/exp/GEN4/LEPTONS/pp45_leptons_tot.list",
    const char* tree_name    = "EpEm_ID",
    Long64_t    chunk_size   = 1000000)
{
    auto* chain = new TChain(tree_name);
    std::string src = input_source;
    if (src.size() >= 5 && src.substr(src.size() - 5) == ".list") {
        std::ifstream fin(src);
        std::string line;
        while (std::getline(fin, line)) {
            const size_t q1 = line.find('"');
            const size_t q2 = line.rfind('"');
            if (q1 != std::string::npos && q2 > q1) {
                chain->Add(line.substr(q1 + 1, q2 - q1 - 1).c_str());
            } else if (!line.empty()) {
                chain->Add(line.c_str());
            }
        }
    } else {
        chain->Add(src.c_str());
    }

    Float_t trigbit_f = 0.0f;
    chain->SetBranchStatus("*", 0);
    chain->SetBranchStatus("trigbit", 1);
    chain->SetBranchAddress("trigbit", &trigbit_f);

    const Long64_t N = chain->GetEntries();
    std::cout << "Source: " << src << "\n"
              << "Tree:   " << tree_name << "\n"
              << "Total entries: " << N << "\n"
              << "Chunk size:    " << chunk_size << "\n\n";

    std::vector<double> X_mid, W_chunk, W_cumul;

    Long64_t chunk_n_PT3 = 0, chunk_n_PT2 = 0;
    Long64_t cumul_n_PT3 = 0, cumul_n_PT2 = 0;
    Long64_t both = 0;

    auto wOf = [](Long64_t n3, Long64_t n2) -> double {
        return (n3 > 0 && n2 > 0) ? (63.0 * static_cast<double>(n2))
                                          / static_cast<double>(n3) : 0.0;
    };

    Long64_t chunk_start = 0;
    for (Long64_t i = 0; i < N; ++i) {
        chain->GetEntry(i);
        const int trigbit = static_cast<int>(trigbit_f);
        if      (trigbit == 8192)  { ++chunk_n_PT3; ++cumul_n_PT3; }
        else if (trigbit == 4096)  { ++chunk_n_PT2; ++cumul_n_PT2; }
        else if (trigbit == 12288) { ++both; }

        if ((i + 1) % chunk_size == 0 || i == N - 1) {
            const double w_c = wOf(chunk_n_PT3, chunk_n_PT2);
            const double w_T = wOf(cumul_n_PT3, cumul_n_PT2);
            X_mid .push_back(0.5 * (chunk_start + i));
            W_chunk.push_back(w_c);
            W_cumul.push_back(w_T);
            printf("  evt [%10lld .. %10lld]  chunk w=%.5f  cumul w=%.5f"
                   "  (n_PT3=%lld n_PT2=%lld)\n",
                   (Long64_t)chunk_start, (Long64_t)i, w_c, w_T,
                   (Long64_t)chunk_n_PT3, (Long64_t)chunk_n_PT2);
            chunk_n_PT3 = chunk_n_PT2 = 0;
            chunk_start = i + 1;
        }
    }

    const double w_final = wOf(cumul_n_PT3, cumul_n_PT2);
    std::cout << "\nFinal cumul:  n_PT3=" << cumul_n_PT3
              << "  n_PT2=" << cumul_n_PT2
              << "  both(12288)=" << both
              << "  w=" << w_final << "\n";

    gStyle->SetOptStat(0);
    auto* g_chunk = new TGraph((int)X_mid.size(), X_mid.data(), W_chunk.data());
    auto* g_cumul = new TGraph((int)X_mid.size(), X_mid.data(), W_cumul.data());
    g_chunk->SetMarkerStyle(20); g_chunk->SetMarkerColor(kBlue+1);
    g_chunk->SetLineColor(kBlue+1); g_chunk->SetMarkerSize(0.9);
    g_cumul->SetMarkerStyle(21); g_cumul->SetMarkerColor(kRed+1);
    g_cumul->SetLineColor(kRed+1); g_cumul->SetLineWidth(2);

    auto* c = new TCanvas("c_w_stab","trigger_corr stability",1100,700);
    c->SetMargin(0.13, 0.05, 0.12, 0.08); c->SetGrid();
    g_chunk->SetTitle("PT3 trigger-bias correction stability;event index;w = (63#upoint N_{PT2})/N_{PT3}");
    g_chunk->Draw("APL");
    g_cumul->Draw("PL SAME");

    auto* lref = new TLine((double)0, w_final, (double)N, w_final);
    lref->SetLineColor(kGray + 2); lref->SetLineStyle(2); lref->SetLineWidth(1);
    lref->Draw();

    auto* leg = new TLegend(0.55, 0.74, 0.94, 0.92);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.034);
    leg->AddEntry(g_chunk, TString::Format("per chunk (%lld events)", chunk_size).Data(), "lp");
    leg->AddEntry(g_cumul, "cumulative",                                                  "lp");
    leg->AddEntry(lref,    TString::Format("final w = %.4f", w_final).Data(),             "l");
    leg->Draw();

    gSystem->mkdir("plots/output", kTRUE);
    c->SaveAs("plots/output/trigger_corr_stability.pdf");
    c->SaveAs("plots/output/trigger_corr_stability.png");
    std::cout << "\n  wrote plots/output/trigger_corr_stability.{pdf,png}\n";
}
