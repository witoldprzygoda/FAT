// mass_ee_adaptive_exp.C — variable-width binning of M_ee from ntuples,
// driven by PT2 statistics (the limiting sample) so each bin holds enough
// events for a meaningful PT2/PT3 ratio. Bin width caps at max_bin_width
// (default 0.2 GeV/c²) — beyond that the algorithm tiles uniformly even
// where statistics fall short.
//
// Three panels on one canvas:
//   1) PT3 spectrum   (trigbit == 8192)
//   2) PT2 spectrum   (trigbit == 4096)
//   3) 63 · N_PT2 / N_PT3 ratio with shared edges
//
// Adaptive-binning algorithm (single forward pass over a 5 MeV/bin seed
// histogram of PT2 events):
//   * accumulate PT2 counts;
//   * when the running sum reaches `N_min_pt2`, OR the running width
//     reaches `max_bin_width`, close the bin at the current 5 MeV edge;
//   * continue until the end of the range; flush any unfinished bin to 1.4.
//
// The edges are then used to fill BOTH PT3 and PT2 from the same ntuple,
// so dividing bin-by-bin is straightforward. Bin contents are RAW COUNTS
// (not counts/GeV) — variable-width inflation of the visual is acceptable
// here because the goal is to see per-bin statistics, not spectral density.
//
// Usage:
//   root -l -b -q plots/mass_ee_adaptive_exp.C
//   root -l -b -q 'plots/mass_ee_adaptive_exp.C(100, 0.15)'   // tighter
//
// Args:
//   N_min_pt2     — minimum PT2 events per bin (default 50)
//   max_bin_width — bin width cap in GeV/c²    (default 0.2)
//   in_file       — source ROOT file           (default output_epem_exp.root)
//   tree_name     — ntuple name                (default dilepton_nt)

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

void mass_ee_adaptive_exp(double      N_min_pt2     = 50.0,
                          double      max_bin_width = 0.2,
                          const char* in_file       = "output_epem_exp.root",
                          const char* tree_name     = "dilepton_nt") {
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    TFile* f = TFile::Open(in_file, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Cannot open " << in_file << "\n";
        return;
    }
    TTree* nt = dynamic_cast<TTree*>(f->Get(tree_name));
    if (!nt) {
        std::cerr << "Tree '" << tree_name << "' not found in " << in_file << "\n";
        f->Close(); return;
    }

    const double mass_max  = 1.4;
    const int    fine_nb   = 280;          // 5 MeV/bin seed
    const double fine_w    = mass_max / fine_nb;

    // -------- Step 1: seed PT2 histogram (fine fixed binning) ----------
    // Leave SetDirectory enabled while Draw fills (Draw's >>name lookup
    // needs the histogram registered in gDirectory). Detach afterwards.
    TH1D* h_seed = new TH1D("h_seed_pt2_fine", "", fine_nb, 0.0, mass_max);
    nt->Draw("m_ee>>h_seed_pt2_fine", "trigbit==4096", "goff");
    h_seed->SetDirectory(nullptr);
    std::cout << "PT2 seed entries: " << h_seed->GetEntries() << "\n";

    // -------- Step 2: walk fine bins, derive adaptive edges -----------
    std::vector<double> edges;
    edges.reserve(fine_nb + 1);
    edges.push_back(0.0);
    double accum_pt2 = 0.0;
    for (int b = 1; b <= fine_nb; ++b) {
        accum_pt2 += h_seed->GetBinContent(b);
        const double upper = h_seed->GetBinLowEdge(b + 1);
        const double width = upper - edges.back();
        if (accum_pt2 >= N_min_pt2 || width >= max_bin_width - 0.5 * fine_w) {
            edges.push_back(upper);
            accum_pt2 = 0.0;
        }
    }
    if (edges.back() < mass_max - 0.5 * fine_w) {
        edges.push_back(mass_max);  // flush trailing partial bin
    }
    const int nb = static_cast<int>(edges.size()) - 1;

    std::cout << "Adaptive binning summary:\n"
              << "  N_min_pt2     = " << N_min_pt2 << " events/bin\n"
              << "  max_bin_width = " << max_bin_width << " GeV/c²\n"
              << "  result        = " << nb << " bins on [0, "
              << mass_max << "]\n";

    delete h_seed;

    // -------- Step 3: PT3 and PT2 with adaptive edges -----------------
    TH1D* h_pt3 = new TH1D(
        "h_m_ee_pt3_adaptive",
        "M_{e^{+}e^{-}} (PT3, adaptive bins);M_{e^{+}e^{-}} [GeV/c^{2}];Counts / bin",
        nb, edges.data());
    h_pt3->Sumw2();
    nt->Draw("m_ee>>h_m_ee_pt3_adaptive", "trigbit==8192", "goff");
    h_pt3->SetDirectory(nullptr);

    TH1D* h_pt2 = new TH1D(
        "h_m_ee_pt2_adaptive",
        "M_{e^{+}e^{-}} (PT2, adaptive bins);M_{e^{+}e^{-}} [GeV/c^{2}];Counts / bin",
        nb, edges.data());
    h_pt2->Sumw2();
    nt->Draw("m_ee>>h_m_ee_pt2_adaptive", "trigbit==4096", "goff");
    h_pt2->SetDirectory(nullptr);

    // -------- Step 4: bin-by-bin 63 · N_PT2/N_PT3 ---------------------
    TH1D* h_ratio = new TH1D(
        "h_m_ee_ratio_pt2_pt3",
        "63 #upoint N_{PT2}/N_{PT3};M_{e^{+}e^{-}} [GeV/c^{2}];63 #upoint N_{PT2}/N_{PT3}",
        nb, edges.data());
    h_ratio->SetDirectory(nullptr);
    double ratio_global_num = 0.0, ratio_global_den = 0.0;
    for (int b = 1; b <= nb; ++b) {
        const double n2 = h_pt2->GetBinContent(b);
        const double e2 = h_pt2->GetBinError(b);
        const double n3 = h_pt3->GetBinContent(b);
        const double e3 = h_pt3->GetBinError(b);
        ratio_global_num += n2;
        ratio_global_den += n3;
        if (n3 > 0) {
            const double r   = 63.0 * n2 / n3;
            const double rel = std::sqrt(
                (n2 > 0 ? std::pow(e2 / n2, 2) : 1.0 / std::max(1.0, n2 + 1)) +
                std::pow(e3 / n3, 2));
            h_ratio->SetBinContent(b, r);
            h_ratio->SetBinError(b, std::abs(r) * rel);
        }
    }
    const double r_global = (ratio_global_den > 0)
        ? 63.0 * ratio_global_num / ratio_global_den : 0.0;

    // -------- Step 5: render 3 panels --------------------------------
    auto styleSeries = [](TH1D* h) {
        h->SetMarkerStyle(20);
        h->SetMarkerSize(0.65);
        h->SetMarkerColor(kBlack);
        h->SetLineColor(kBlack);
        h->GetXaxis()->SetTitleSize(0.050);
        h->GetYaxis()->SetTitleSize(0.050);
        h->GetXaxis()->SetLabelSize(0.045);
        h->GetYaxis()->SetLabelSize(0.045);
        h->GetXaxis()->SetTitleOffset(1.05);
        h->GetYaxis()->SetTitleOffset(1.15);
    };
    styleSeries(h_pt3);
    styleSeries(h_pt2);
    styleSeries(h_ratio);

    TCanvas* c = new TCanvas("c_mass_ee_adaptive",
                             "M_{ee} adaptive binning (PT3 / PT2 / ratio)",
                             1800, 600);
    c->Divide(3, 1, 0.001, 0.001);

    // Pad 1: PT3
    c->cd(1);
    gPad->SetLogy();
    gPad->SetMargin(0.13, 0.04, 0.13, 0.09);
    h_pt3->GetYaxis()->SetRangeUser(0.5, h_pt3->GetMaximum() * 4.0);
    h_pt3->Draw("E");

    // Pad 2: PT2
    c->cd(2);
    gPad->SetLogy();
    gPad->SetMargin(0.13, 0.04, 0.13, 0.09);
    h_pt2->GetYaxis()->SetRangeUser(0.5, h_pt2->GetMaximum() * 4.0);
    h_pt2->Draw("E");

    // Pad 3: 63 · PT2/PT3 ratio with auto Y range (with safety margin).
    c->cd(3);
    gPad->SetMargin(0.13, 0.04, 0.13, 0.09);
    double y_max_seen = 0.0;
    double y_min_seen = std::numeric_limits<double>::infinity();
    for (int b = 1; b <= nb; ++b) {
        const double v = h_ratio->GetBinContent(b);
        const double e = h_ratio->GetBinError(b);
        if (h_ratio->GetBinContent(b) != 0.0 || h_ratio->GetBinError(b) != 0.0) {
            y_max_seen = std::max(y_max_seen, v + e);
            y_min_seen = std::min(y_min_seen, std::max(0.0, v - e));
        }
    }
    if (!std::isfinite(y_min_seen)) y_min_seen = 0.0;
    if (y_max_seen <= 0.0)          y_max_seen = 5.0;
    const double y_lo = std::max(0.0, y_min_seen - 0.5);
    const double y_hi = y_max_seen * 1.20;
    h_ratio->GetYaxis()->SetRangeUser(y_lo, y_hi);
    h_ratio->Draw("E");

    // Reference lines: global ratio (= calibration <w> for this channel)
    // and the no-correction reference at 1.0.
    TLine* l_ref = new TLine(0.0, r_global, mass_max, r_global);
    l_ref->SetLineStyle(2);
    l_ref->SetLineColor(kBlue + 1);
    l_ref->SetLineWidth(2);
    l_ref->Draw();
    TLine* l_one = new TLine(0.0, 1.0, mass_max, 1.0);
    l_one->SetLineStyle(3);
    l_one->SetLineColor(kGray + 2);
    l_one->Draw();

    TLegend* leg = new TLegend(0.42, 0.74, 0.95, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.038);
    leg->AddEntry(h_ratio, "63 #upoint N_{PT2}/N_{PT3} per bin", "lpe");
    leg->AddEntry(l_ref,
        Form("global = %.3f (calibration #LTw#GT)", r_global), "l");
    leg->AddEntry(l_one,  "1.0 (no correction)", "l");
    leg->Draw();

    // ------------------------------------------------------------------
    // Diagnostics: print bin edges + counts so user can verify the
    // binning is reasonable (and tune N_min_pt2 / max_bin_width).
    // ------------------------------------------------------------------
    std::cout << "\nBin   M_low   M_high  width    PT3      PT2     63·PT2/PT3\n"
              << "-----------------------------------------------------------\n";
    for (int b = 1; b <= nb; ++b) {
        const double lo = h_ratio->GetBinLowEdge(b);
        const double hi = h_ratio->GetBinLowEdge(b + 1);
        const double n3 = h_pt3->GetBinContent(b);
        const double n2 = h_pt2->GetBinContent(b);
        const double r  = h_ratio->GetBinContent(b);
        std::cout << std::setw(3) << b << "  "
                  << std::fixed << std::setprecision(3)
                  << std::setw(6) << lo  << "  "
                  << std::setw(6) << hi  << "  "
                  << std::setw(6) << (hi - lo) << "  "
                  << std::setw(8) << static_cast<long>(n3) << "  "
                  << std::setw(7) << static_cast<long>(n2) << "  "
                  << std::setw(7) << std::setprecision(3) << r << "\n";
    }
    std::cout << "\nGlobal 63 · PT2/PT3 = " << r_global
              << "   (this channel's calibrated <w>)\n";

    gSystem->mkdir("plots/output", true);
    c->SaveAs("plots/output/mass_ee_adaptive.pdf");
    c->SaveAs("plots/output/mass_ee_adaptive.png");
    std::cout << "\nSaved: plots/output/mass_ee_adaptive.{pdf,png}\n";

    f->Close();
}
