// research/plots/slices.C — draw m_epemg with CB extraction.
//
// Reads three files produced by `meson_research` (one per CB channel):
//   outputs/research_epem.root   (signal — opposite-sign)
//   outputs/research_epep.root   (++ pairs — combinatorial)
//   outputs/research_emem.root   (-- pairs — combinatorial)
//
// Plots:
//   1) the integrated full spectrum (m_epemg_full) — first, then
//   2) one panel per OA slice (m_epemg_oa_X_Y).
//
// For each panel: black (all) / red (CB = 2√(N++ N--)) / blue (signal).
//
// Output: per-panel PDF/PNG in plots/output/, plus a single multipage
// PDF plots/output/slices_m_epemg.pdf for browsing.
// All paths are relative to the research/ directory.
//
// Usage (from research/):
//   root -l -b -q plots/slices.C

#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>
#include <TSystem.h>
#include <TStyle.h>
#include <iostream>
#include <cmath>
#include <string>
#include <vector>

namespace {
    // Match the slicing of meson_research.cc — keep these in sync if the
    // research binary's grid changes.
    constexpr double kSliceMin  = 0.0;
    constexpr double kSliceMax  = 10.0;
    constexpr double kSliceStep = 0.2;

    constexpr const char* kVarStem = "m_epemg";   // histogram base name in files

    std::string fmtEdge(double x) {
        char buf[32];
        snprintf(buf, sizeof(buf), "%.1f", x);
        std::string s(buf);
        for (auto& c : s) if (c == '.') c = 'p';
        return s;
    }
}

void slices() {

    TFile* f_all = TFile::Open("outputs/research_epem.root", "READ");
    TFile* f_pp  = TFile::Open("outputs/research_epep.root", "READ");
    TFile* f_mm  = TFile::Open("outputs/research_emem.root", "READ");

    if (!f_all || f_all->IsZombie()) { std::cerr << "Cannot open research_epem.root\n"; return; }
    if (!f_pp  || f_pp ->IsZombie()) { std::cerr << "Cannot open research_epep.root\n"; return; }
    if (!f_mm  || f_mm ->IsZombie()) { std::cerr << "Cannot open research_emem.root\n"; return; }

    gStyle->SetOptStat(0);

    const int n_slices = static_cast<int>(std::round((kSliceMax - kSliceMin) / kSliceStep));
    const int n_panels = 1 + n_slices;   // 1 full spectrum + n slice spectra

    gSystem->mkdir("plots/output", kTRUE);
    const std::string pdf_multi = "plots/output/slices_m_epemg.pdf";

    // We collect canvases over the full panel set and only append to the
    // multipage PDF in a separate second pass. Mixing per-panel SaveAs
    // with multipage Print confuses ROOT's TPDF state machine and yields
    // "An Object is already open" / "No Object currently opened" warnings.
    std::vector<TCanvas*> canvases;
    canvases.reserve(n_panels);

    // -- One-panel plotter (CB extraction + canvas + per-panel save).
    //    Multipage Print is deferred to the second pass below.
    auto plot_one = [&](const std::string& hname,
                        const std::string& ctitle,
                        const std::string& save_base)
    {
        auto* h_all_in = (TH1D*)f_all->Get(hname.c_str());
        auto* h_pp_in  = (TH1D*)f_pp ->Get(hname.c_str());
        auto* h_mm_in  = (TH1D*)f_mm ->Get(hname.c_str());

        if (!h_all_in || !h_pp_in || !h_mm_in) {
            std::cerr << "missing " << hname << " in one of the inputs — skip\n";
            return;
        }

        // Detach clones so closing input files later won't kill them.
        auto* h_all = (TH1D*)h_all_in->Clone((hname + "_all").c_str());
        auto* h_pp  = (TH1D*)h_pp_in ->Clone((hname + "_pp" ).c_str());
        auto* h_mm  = (TH1D*)h_mm_in ->Clone((hname + "_mm" ).c_str());
        h_all->SetDirectory(nullptr);
        h_pp ->SetDirectory(nullptr);
        h_mm ->SetDirectory(nullptr);

        // CB = 2 * sqrt(N_++ * N_--), bin-by-bin, with proper error propagation.
        TH1D* h_cb = (TH1D*)h_all->Clone((hname + "_cb").c_str());
        h_cb->Reset();
        h_cb->SetDirectory(nullptr);
        for (int b = 1; b <= h_cb->GetNbinsX(); ++b) {
            const double n_pp = h_pp->GetBinContent(b);
            const double n_mm = h_mm->GetBinContent(b);
            const double e_pp = h_pp->GetBinError(b);
            const double e_mm = h_mm->GetBinError(b);

            double cb = 0.0, err = 0.0;
            if (n_pp > 0 && n_mm > 0) {
                cb = 2.0 * std::sqrt(n_pp * n_mm);
                const double t1 = e_pp * std::sqrt(n_mm / n_pp);
                const double t2 = e_mm * std::sqrt(n_pp / n_mm);
                err = std::sqrt(t1 * t1 + t2 * t2);
            }
            h_cb->SetBinContent(b, cb);
            h_cb->SetBinError(b, err);
        }

        // Signal = all - CB
        TH1D* h_sig = (TH1D*)h_all->Clone((hname + "_sig").c_str());
        h_sig->SetDirectory(nullptr);
        h_sig->Add(h_cb, -1.0);

        // Style: black/red/blue markers with error bars.
        for (auto* h : {h_all, h_cb, h_sig}) {
            h->SetMarkerStyle(20);
            h->SetMarkerSize(0.6);
        }
        h_all->SetMarkerColor(kBlack); h_all->SetLineColor(kBlack);
        h_cb ->SetMarkerColor(kRed);   h_cb ->SetLineColor(kRed);
        h_sig->SetMarkerColor(kBlue);  h_sig->SetLineColor(kBlue);

        const TString cname = TString::Format("c_%s", hname.c_str());
        auto* c = new TCanvas(cname, hname.c_str(), 800, 600);
        c->SetMargin(0.12, 0.05, 0.12, 0.08);

        h_all->SetTitle(ctitle.c_str());
        const double y_max = std::max({h_all->GetMaximum(),
                                       h_cb->GetMaximum(),
                                       h_sig->GetMaximum()});
        h_all->SetMaximum(y_max * 1.2);
        h_all->SetMinimum(0.0);

        h_all->Draw("E");
        h_cb ->Draw("E SAME");
        h_sig->Draw("E SAME");

        auto* leg = new TLegend(0.62, 0.72, 0.94, 0.90);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.035);
        leg->AddEntry(h_all, "all (e^{+}e^{-})",            "lpe");
        leg->AddEntry(h_cb,  "CB (2#sqrt{N_{++}N_{--}})",   "lpe");
        leg->AddEntry(h_sig, "signal (all - CB)",           "lpe");
        leg->Draw();

        c->Update();

        // Per-panel PDF + PNG (multipage append happens in second pass).
        const std::string per = "plots/output/" + save_base;
        c->SaveAs((per + ".pdf").c_str());
        c->SaveAs((per + ".png").c_str());

        std::cout << "  " << hname
                  << ":  all=" << h_all->Integral()
                  << "   CB="  << h_cb ->Integral()
                  << "   sig=" << h_sig->Integral() << "\n";

        canvases.push_back(c);  // kept alive for multipage pass below
    };

    // === PASS 1: build canvases + per-panel PDF/PNG ========================

    // Panel 0: full spectrum (integrated over OA ∈ slice range)
    {
        const std::string hname = std::string(kVarStem) + "_full";
        const TString ctitle = TString::Format(
            "M(e^{+}e^{-}#gamma), OA #in [%.1f, %.1f] deg (full);"
            "M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Counts",
            kSliceMin, kSliceMax);
        plot_one(hname, ctitle.Data(), "full_" + hname);
    }

    // Panels 1..n_slices: per-slice
    for (int i = 0; i < n_slices; ++i) {
        const double oa_lo = kSliceMin + i * kSliceStep;
        const double oa_hi = kSliceMin + (i + 1) * kSliceStep;
        const std::string hname = std::string(kVarStem) + "_oa_"
                                + fmtEdge(oa_lo) + "_" + fmtEdge(oa_hi);
        const TString ctitle = TString::Format(
            "M(e^{+}e^{-}#gamma), OA #in [%.1f, %.1f] deg;"
            "M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Counts",
            oa_lo, oa_hi);
        plot_one(hname, ctitle.Data(), "slice_" + hname);
    }

    // === PASS 2: assemble multipage PDF ====================================
    // Done with all per-panel TPDF activity; the multipage TPDF stays open
    // across this loop without interference.
    for (size_t i = 0; i < canvases.size(); ++i) {
        auto* c = canvases[i];
        if (i == 0)                          c->Print((pdf_multi + "(").c_str(), "pdf");
        else if (i == canvases.size() - 1)   c->Print((pdf_multi + ")").c_str(), "pdf");
        else                                 c->Print(pdf_multi.c_str(),         "pdf");
    }

    for (auto* c : canvases) delete c;

    f_all->Close();
    f_pp ->Close();
    f_mm ->Close();

    std::cout << "\nMultipage PDF: " << pdf_multi
              << "  (" << canvases.size() << " pages)\n";
    std::cout << "Per-panel PDFs/PNGs in plots/output/\n";
}
