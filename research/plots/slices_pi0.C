// research/plots/slices_pi0.C — draw m_epemg per OA slice (EXP DATA).
//
// All/CB/signal triple per slice with the π⁰ Dalitz region zoomed in.
// Reads three files produced by `meson_research` (one per CB channel):
//   research_epem.root   (signal — opposite-sign)
//   research_epep.root   (++ pairs — combinatorial)
//   research_emem.root   (-- pairs — combinatorial)
//
// CB definition: bin-by-bin
//     CB(b)        = 2 √(N_++(b) · N_--(b))
//     σ_CB(b)²     = (e_++ √(N_--/N_++))² + (e_-- √(N_++/N_--))²
// Signal = all − CB (errors via TH1::Add with c2=-1).
//
// Plots, in order:
//   1) integrated full spectrum (m_epemg_full), then
//   2) one panel per OA slice (m_epemg_oa_X_Y).
//
// Each panel overlays:
//   black markers  — all (e+e-)
//   red markers    — CB (combinatorial)
//   blue markers   — signal (all − CB)
//
// Output: per-panel PDF/PNG in plots/output/ + multipage PDF
// plots/output/slices_pi0_<flavour>.pdf for browsing.
//
// Usage (from research/):
//   root -l -b -q plots/slices_pi0.C            # REC (default)
//   root -l -b -q 'plots/slices_pi0.C("cor")'   # COR

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
    // Match the slicing of meson_research.cc — keep these in sync.
    constexpr double kSliceMin  = 0.0;
    constexpr double kSliceMax  = 15.0;
    constexpr double kSliceStep = 0.2;

    // π⁰ Dalitz zoom on the X axis (display only — histograms keep [0, 0.8]).
    constexpr double kZoomMin = 0.0;
    constexpr double kZoomMax = 0.46;

    std::string fmtEdge(double x) {
        char buf[32];
        snprintf(buf, sizeof(buf), "%.1f", x);
        std::string s(buf);
        for (auto& c : s) if (c == '.') c = 'p';
        return s;
    }
}

void slices_pi0(const char* flavour = "rec") {

    // Pick mass flavour ('rec' → 'm_epemg', 'cor' → 'm_epemg_cor').
    std::string fl = flavour;
    for (auto& c : fl) c = std::tolower(c);
    std::string var_stem;
    if      (fl == "rec") var_stem = "m_epemg";
    else if (fl == "cor") var_stem = "m_epemg_cor";
    else { std::cerr << "Unknown flavour '" << flavour
                     << "' (expected 'rec' or 'cor')\n"; return; }

    TFile* f_all = TFile::Open("research_epem.root", "READ");
    TFile* f_pp  = TFile::Open("research_epep.root", "READ");
    TFile* f_mm  = TFile::Open("research_emem.root", "READ");
    if (!f_all || f_all->IsZombie()) { std::cerr << "Cannot open research_epem.root\n"; return; }
    if (!f_pp  || f_pp ->IsZombie()) { std::cerr << "Cannot open research_epep.root\n"; return; }
    if (!f_mm  || f_mm ->IsZombie()) { std::cerr << "Cannot open research_emem.root\n"; return; }

    gStyle->SetOptStat(0);

    const int n_slices = static_cast<int>(std::round((kSliceMax - kSliceMin) / kSliceStep));

    gSystem->mkdir("plots/output", kTRUE);
    const std::string pdf_multi =
        "plots/output/slices_pi0_" + fl + ".pdf";

    std::vector<TCanvas*> canvases;
    canvases.reserve(1 + n_slices);

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

        auto* h_all = (TH1D*)h_all_in->Clone((hname + "_all").c_str());
        auto* h_pp  = (TH1D*)h_pp_in ->Clone((hname + "_pp" ).c_str());
        auto* h_mm  = (TH1D*)h_mm_in ->Clone((hname + "_mm" ).c_str());
        h_all->SetDirectory(nullptr);
        h_pp ->SetDirectory(nullptr);
        h_mm ->SetDirectory(nullptr);

        // CB = 2 √(N_++ N_--) per bin, with proper error propagation.
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

        // Signal = all - CB (TH1::Add propagates errors in quadrature).
        TH1D* h_sig = (TH1D*)h_all->Clone((hname + "_sig").c_str());
        h_sig->SetDirectory(nullptr);
        h_sig->Add(h_cb, -1.0);

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
        h_all->GetXaxis()->SetRangeUser(kZoomMin, kZoomMax);
        const double y_max = std::max({h_all->GetMaximum(),
                                       h_cb ->GetMaximum(),
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
        leg->AddEntry(h_sig, "signal (all #minus CB)",      "lpe");
        leg->Draw();

        c->Update();
        const std::string per = "plots/output/" + save_base;
        c->SaveAs((per + ".pdf").c_str());
        c->SaveAs((per + ".png").c_str());

        std::cout << "  " << hname
                  << ":  all=" << h_all->Integral()
                  << "  CB="   << h_cb ->Integral()
                  << "  sig="  << h_sig->Integral() << "\n";

        canvases.push_back(c);
    };

    // === PASS 1: per-panel canvases ========================================
    {
        const std::string hname = var_stem + "_full";
        const TString ctitle = TString::Format(
            "M(e^{+}e^{-}#gamma) %s, OA #in [%.1f, %.1f] deg (full);"
            "M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Counts",
            fl.c_str(), kSliceMin, kSliceMax);
        plot_one(hname, ctitle.Data(), "full_pi0_" + fl + "_" + hname);
    }

    for (int i = 0; i < n_slices; ++i) {
        const double oa_lo = kSliceMin + i * kSliceStep;
        const double oa_hi = kSliceMin + (i + 1) * kSliceStep;
        const std::string hname = var_stem + "_oa_" + fmtEdge(oa_lo) + "_" + fmtEdge(oa_hi);
        const TString ctitle = TString::Format(
            "M(e^{+}e^{-}#gamma) %s, OA #in [%.1f, %.1f] deg;"
            "M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Counts",
            fl.c_str(), oa_lo, oa_hi);
        plot_one(hname, ctitle.Data(), "slice_pi0_" + fl + "_" + hname);
    }

    // === PASS 2: multipage PDF ==============================================
    for (size_t i = 0; i < canvases.size(); ++i) {
        auto* c = canvases[i];
        if (i == 0)                          c->Print((pdf_multi + "(").c_str(), "pdf");
        else if (i == canvases.size() - 1)   c->Print((pdf_multi + ")").c_str(), "pdf");
        else                                 c->Print(pdf_multi.c_str(),         "pdf");
    }

    for (auto* c : canvases) delete c;
    f_all->Close();  f_pp->Close();  f_mm->Close();

    std::cout << "\nMultipage PDF: " << pdf_multi
              << "  (" << canvases.size() << " pages)\n";
    std::cout << "Per-panel PDFs/PNGs in plots/output/\n";
}
