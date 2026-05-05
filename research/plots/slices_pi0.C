// research/plots/slices_pi0.C — draw m_epemg per OA slice (SIMULATION).
//
// Histogram range matches meson_research.cc: m_epemg in [0, 0.46] GeV/c²
// at 230 bins (2 MeV/bin), the resolution and zoom intended for the
// pi0-Dalitz study (peak ~0.135 GeV).
//
// Reads a single file produced by `meson_research`:
//   research_sim.root
//
// Plots, in order:
//   1) the integrated full spectrum (m_epemg_full + m_epemg_sim_full), then
//   2) one panel per OA slice (m_epemg_oa_X_Y + m_epemg_sim_oa_X_Y).
//
// Each panel overlays:
//   black markers — REC mass (m_epemg)
//   brown step    — SIM-truth mass (m_epemg_sim)
//
// All histograms are already weighted by sim_genweight at fill time inside
// the binary, so no further per-bin reweighting is applied here.
//
// Output: per-panel PDF/PNG in plots/output/, plus a single multipage PDF
// plots/output/slices_m_epemg.pdf for browsing.
// All paths are relative to the research/ directory.
//
// Usage (from research/):
//   root -l -b -q plots/slices_pi0.C

#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TColor.h>
#include <iostream>
#include <cmath>
#include <string>
#include <vector>

namespace {
    // Match the slicing of meson_research.cc — keep these in sync if the
    // research binary's grid changes.
    constexpr double kSliceMin  = 0.0;
    constexpr double kSliceMax  = 15.0;
    constexpr double kSliceStep = 0.2;

    constexpr const char* kVarStem    = "m_epemg";       // REC base name
    constexpr const char* kVarStemSim = "m_epemg_tru";   // TRU (truth) base name

    std::string fmtEdge(double x) {
        char buf[32];
        snprintf(buf, sizeof(buf), "%.1f", x);
        std::string s(buf);
        for (auto& c : s) if (c == '.') c = 'p';
        return s;
    }

    // Brown step line for the SIM-truth overlay — matches the convention
    // used by the JointPlotter on pp45_pippimepem_compare / pp45_epem_compare.
    int simBrown() {
        static int c = TColor::GetColor("#8B4513");   // SaddleBrown
        return c;
    }
}

void slices_pi0() {

    TFile* f = TFile::Open("research_sim.root", "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Cannot open research_sim.root\n";
        return;
    }

    gStyle->SetOptStat(0);

    const int n_slices = static_cast<int>(std::round((kSliceMax - kSliceMin) / kSliceStep));
    const int n_panels = 1 + n_slices;

    gSystem->mkdir("plots/output", kTRUE);
    const std::string pdf_multi = "plots/output/slices_m_epemg.pdf";

    // We collect canvases over the full panel set and only append to the
    // multipage PDF in a separate second pass. Mixing per-panel SaveAs with
    // multipage Print confuses ROOT's TPDF state machine.
    std::vector<TCanvas*> canvases;
    canvases.reserve(n_panels);

    auto plot_one = [&](const std::string& hname_rec,
                        const std::string& hname_sim,
                        const std::string& ctitle,
                        const std::string& save_base)
    {
        auto* h_rec_in = (TH1D*)f->Get(hname_rec.c_str());
        auto* h_sim_in = (TH1D*)f->Get(hname_sim.c_str());

        if (!h_rec_in) {
            std::cerr << "missing " << hname_rec << " — skip\n";
            return;
        }

        // Detach clones so closing the input file later won't kill them.
        auto* h_rec = (TH1D*)h_rec_in->Clone((hname_rec + "_c").c_str());
        h_rec->SetDirectory(nullptr);
        TH1D* h_sim = nullptr;
        if (h_sim_in) {
            h_sim = (TH1D*)h_sim_in->Clone((hname_sim + "_c").c_str());
            h_sim->SetDirectory(nullptr);
        }

        // Style: REC = black markers, SIM truth = brown step line.
        h_rec->SetMarkerStyle(20);
        h_rec->SetMarkerSize(0.6);
        h_rec->SetMarkerColor(kBlack);
        h_rec->SetLineColor(kBlack);
        if (h_sim) {
            h_sim->SetMarkerStyle(0);
            h_sim->SetMarkerSize(0);
            h_sim->SetLineColor(simBrown());
            h_sim->SetLineWidth(2);
            h_sim->SetFillStyle(0);
        }

        const TString cname = TString::Format("c_%s", hname_rec.c_str());
        auto* c = new TCanvas(cname, hname_rec.c_str(), 800, 600);
        c->SetMargin(0.12, 0.05, 0.12, 0.08);

        h_rec->SetTitle(ctitle.c_str());
        // Linear Y: tight headroom (1.2×) and floor at 0.
        const double y_max = std::max(h_rec->GetMaximum(),
                                      h_sim ? h_sim->GetMaximum() : 0.0);
        h_rec->SetMaximum(y_max * 1.2);
        h_rec->SetMinimum(0.0);

        h_rec->Draw("E");
        if (h_sim) h_sim->Draw("HIST SAME");
        // Re-draw markers on top of the step so they aren't visually hidden.
        h_rec->Draw("E SAME");

        auto* leg = new TLegend(0.62, 0.80, 0.94, 0.90);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.035);
        leg->AddEntry(h_rec, "REC (m_{e^{+}e^{-}#gamma})", "lpe");
        if (h_sim) leg->AddEntry(h_sim, "SIM truth", "l");
        leg->Draw();

        c->Update();

        // Per-panel PDF + PNG (multipage append happens in second pass).
        const std::string per = "plots/output/" + save_base;
        c->SaveAs((per + ".pdf").c_str());
        c->SaveAs((per + ".png").c_str());

        std::cout << "  " << hname_rec
                  << ":  REC=" << h_rec->Integral()
                  << "   SIM=" << (h_sim ? h_sim->Integral() : 0.0) << "\n";

        canvases.push_back(c);   // keep alive for multipage pass below
    };

    // === PASS 1: build canvases + per-panel PDF/PNG ========================

    // Panel 0: full spectrum (integrated over OA ∈ slice range)
    {
        const std::string hname_rec = std::string(kVarStem)    + "_full";
        const std::string hname_sim = std::string(kVarStemSim) + "_full";
        const TString ctitle = TString::Format(
            "M(e^{+}e^{-}#gamma), OA #in [%.1f, %.1f] deg (full);"
            "M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Counts (sim_genweight)",
            kSliceMin, kSliceMax);
        plot_one(hname_rec, hname_sim, ctitle.Data(), "full_" + hname_rec);
    }

    // Panels 1..n_slices: per-slice
    for (int i = 0; i < n_slices; ++i) {
        const double oa_lo = kSliceMin + i * kSliceStep;
        const double oa_hi = kSliceMin + (i + 1) * kSliceStep;
        const std::string suff = "oa_" + fmtEdge(oa_lo) + "_" + fmtEdge(oa_hi);
        const std::string hname_rec = std::string(kVarStem)    + "_" + suff;
        const std::string hname_sim = std::string(kVarStemSim) + "_" + suff;
        const TString ctitle = TString::Format(
            "M(e^{+}e^{-}#gamma), OA #in [%.1f, %.1f] deg;"
            "M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Counts (sim_genweight)",
            oa_lo, oa_hi);
        plot_one(hname_rec, hname_sim, ctitle.Data(), "slice_" + hname_rec);
    }

    // === PASS 2: assemble multipage PDF ====================================
    for (size_t i = 0; i < canvases.size(); ++i) {
        auto* c = canvases[i];
        if (i == 0)                          c->Print((pdf_multi + "(").c_str(), "pdf");
        else if (i == canvases.size() - 1)   c->Print((pdf_multi + ")").c_str(), "pdf");
        else                                 c->Print(pdf_multi.c_str(),         "pdf");
    }

    for (auto* c : canvases) delete c;

    f->Close();

    std::cout << "\nMultipage PDF: " << pdf_multi
              << "  (" << canvases.size() << " pages)\n";
    std::cout << "Per-panel PDFs/PNGs in plots/output/\n";
}
