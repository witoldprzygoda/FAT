// pippimepemg_spectra_tru_sim.C — eta-candidate quality study using SIMULATED-truth kinematics.
// Mirror of pippimepemg_spectra_cor_sim.C; reads m_pippimepem_sim / m_pippimepemg_sim
// from pippimepem_nt_cor and applies sim_genweight.
//
// Caveat: the photon used to build epemg has no MC-truth field in the SMASH
// ntuple — we mirror its RECONSTRUCTED 4-vector into SIMULATED in main.cc, so
// "sim" here means "leptons + pions from Geant truth, photon from ECAL". The
// resolution improvement vs _cor reflects the lepton/pion-side smearing only.
//
// Two sections (sim mode = no CB):
//   (A) OVERLAY — for each cut configuration, draw M(pippimepem)_sim for the
//       full sample and overlay the M(pippimepem)_sim of the pippimepemg-tagged
//       subset rescaled to match the full-sample integral in M ~ [0.45, 0.55] GeV.
//   (B) STANDALONE — M(pippimepemg)_sim under the same cut configurations.
//
// All output PDFs/PNGs land in plots/output/ with prefix "pippimepemg_*_tru_sim".
//
// Usage: root -l -b -q plots/pippimepemg_spectra_tru_sim.C

#include "PlotUtils.h"
#include <sstream>

namespace {
    constexpr double kProbeLo = 0.45;
    constexpr double kProbeHi = 0.55;

    // Helper: compose a weighted TTree::Draw cut. Per-event sim_genweight from
    // the ntuple is multiplied by an optional boolean filter so all spectra
    // reflect the SMASH luminosity normalisation.
    std::string wcut(const std::string& filter = "") {
        return filter.empty() ? std::string("sim_genweight")
                              : "(" + filter + ")*sim_genweight";
    }
}

void pippimepemg_spectra_tru_sim() {

    PlotUtils pu("output_pippimepem_sim.root");   // single-file (sim, no CB)

    const char* NT = "pippimepem_nt_cor";

    // -----------------------------------------------------------------
    // SECTION A: OVERLAY — M(pippimepem)_sim, full sample vs pippimepemg-tagged
    // -----------------------------------------------------------------
    auto drawOverlay = [&](const std::string& cut_base, const std::string& title,
                           const std::string& tag, int nbins) {
        std::string cut_tagged = cut_base.empty()
                                 ? std::string("pippimepemg_pass_narrow==1")
                                 : cut_base + " && pippimepemg_pass_narrow==1";

        auto* h_full = pu.drawNtupleSingle(
            NT, "m_pippimepem_sim", nbins, 0.2, 1.4, wcut(cut_base),
            ";M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}] (truth);a.u.");
        auto* h_tagged = pu.drawNtupleSingle(
            NT, "m_pippimepem_sim", nbins, 0.2, 1.4, wcut(cut_tagged),
            ";M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}] (truth);a.u.");

        if (!h_full || !h_tagged) return;

        // Rescale tagged to full in the probe window
        int blo = h_full->FindBin(kProbeLo);
        int bhi = h_full->FindBin(kProbeHi);
        double y_full   = h_full->Integral(blo, bhi);
        double y_tagged = h_tagged->Integral(blo, bhi);
        double scale    = (y_tagged > 0) ? y_full / y_tagged : 1.0;
        h_tagged->Scale(scale);

        TCanvas* cv = new TCanvas(("c_overlay_" + tag + "_tru_sim").c_str(), title.c_str(), 800, 600);
        cv->SetMargin(0.12, 0.05, 0.12, 0.08);

        h_full  ->SetMarkerStyle(20); h_full  ->SetMarkerSize(0.7);
        h_full  ->SetMarkerColor(kBlack);    h_full  ->SetLineColor(kBlack);
        h_tagged->SetMarkerStyle(24); h_tagged->SetMarkerSize(0.8);
        h_tagged->SetMarkerColor(kGreen+2);  h_tagged->SetLineColor(kGreen+2);

        h_full->SetTitle(title.c_str());
        h_full->SetMaximum(h_full->GetMaximum() * 1.25);
        h_full->SetMinimum(0);

        h_full  ->Draw("E");
        h_tagged->Draw("E SAME");

        TLegend* leg = new TLegend(0.45, 0.74, 0.93, 0.90);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.030);
        leg->AddEntry(h_full,   "M(#pi^{+}#pi^{-}e^{+}e^{-}) full sample (truth)", "lpe");
        leg->AddEntry(h_tagged, Form("#pi^{+}#pi^{-}#pi^{0}-tagged #times %.3g", scale), "lpe");
        leg->Draw();

        cv->Update();
        pu.save(cv, "pippimepemg_overlay_" + tag + "_tru_sim");

        std::cout << "[overlay " << tag << " truth] probe ["
                  << kProbeLo << ", " << kProbeHi << "] GeV/c^2:"
                  << "  full=" << y_full
                  << "  tagged=" << y_tagged
                  << "  scale=" << scale << "\n";
    };

    // -----------------------------------------------------------------
    // SECTION B: STANDALONE — M(pippimepemg)_sim
    // -----------------------------------------------------------------
    auto drawStandalone = [&](const std::string& cut_base, const std::string& title,
                              const std::string& tag, int nbins) {
        std::string cut_g = cut_base.empty() ? std::string("pippimepemg_pass_narrow==1")
                                             : cut_base + " && pippimepemg_pass_narrow==1";
        auto* h = pu.drawNtupleSingle(
            NT, "m_pippimepemg_sim", nbins, 0.2, 1.4, wcut(cut_g),
            ";M_{#pi^{+}#pi^{-}#pi^{0}} [GeV/c^{2}] (truth);a.u.");
        auto* cv = pu.drawSingle(h, title, "c_pippimepemg_" + tag + "_tru_sim");
        pu.save(cv, "pippimepemg_" + tag + "_tru_sim");
        std::cout << "[standalone " << tag << " truth] entries=" << h->GetEntries()
                  << "  integral=" << h->Integral() << "\n";
    };

    // -----------------------------------------------------------------
    // SECTION C: MM(pippimepemg)_sim — beam+target − pi+pi-e+e-γ (truth)
    // -----------------------------------------------------------------
    auto drawStandaloneMM = [&](const std::string& cut_base, const std::string& title,
                                 const std::string& tag, int nbins) {
        std::string cut_g = cut_base.empty() ? std::string("pippimepemg_pass_narrow==1")
                                             : cut_base + " && pippimepemg_pass_narrow==1";
        auto* h = pu.drawNtupleSingle(
            NT, "mm_pippimepemg_sim", nbins, 0.0, 4.0, wcut(cut_g),
            ";MM(#pi^{+}#pi^{-}e^{+}e^{-}#gamma) [GeV/c^{2}] (truth);a.u.");
        auto* cv = pu.drawSingle(h, title, "c_pippimepemg_mm_" + tag + "_tru_sim");
        pu.save(cv, "pippimepemg_mm_" + tag + "_tru_sim");
        std::cout << "[MM " << tag << " truth] entries=" << h->GetEntries()
                  << "  integral=" << h->Integral() << "\n";
    };

    // -----------------------------------------------------------------
    // Cut configurations
    // -----------------------------------------------------------------
    auto run = [&](const std::string& cut_base, const std::string& title,
                   const std::string& tag, int nbins = 200) {
        drawOverlay(cut_base, title + " (overlay, truth)", tag, nbins);
        drawStandalone(cut_base, title + " (#pi^{+}#pi^{-}#pi^{0}, truth)", tag, nbins);
        drawStandaloneMM(cut_base, "MM(#pi^{+}#pi^{-}e^{+}e^{-}#gamma) " + title + " (truth)",
                         tag, nbins);
    };

    run("",                                 "M (no cut)",                                   "base",                  200);
    run("sel_pass==1",                      "M (after selection)",                          "selected",              200);
    run("sel_pass==1 && cut2d_pass==1",     "M (after selection + cut_2d)",                 "selected_cut2d",        200);

    auto runSlice = [&](double lo, double hi, const std::string& tag_suffix, bool with_cut2d) {
        std::ostringstream cut;
        cut << "sel_pass==1 && mm_pippimepem>=" << lo
            << " && mm_pippimepem<=" << hi;
        if (with_cut2d) cut << " && cut2d_pass==1";

        std::ostringstream title;
        title << "M, sel" << (with_cut2d ? "+cut_2d" : "")
              << ", MM #in [" << lo << ", " << hi << "] GeV/c^{2}";

        std::string tag = "slice_" + tag_suffix + (with_cut2d ? "_cut2d" : "");
        run(cut.str(), title.str(), tag, 100);
    };

    // Slices without cut_2d
    runSlice(2.0, 2.2, "20_22", false);
    runSlice(2.2, 2.4, "22_24", false);
    runSlice(2.4, 2.6, "24_26", false);
    runSlice(2.6, 2.8, "26_28", false);
    runSlice(2.8, 3.0, "28_30", false);

    // Slices with cut_2d
    runSlice(2.0, 2.2, "20_22", true);
    runSlice(2.2, 2.4, "22_24", true);
    runSlice(2.4, 2.6, "24_26", true);
    runSlice(2.6, 2.8, "26_28", true);
    runSlice(2.8, 3.0, "28_30", true);

    std::cout << "\nDone. Plots saved to plots/output/pippimepemg_*_tru_sim\n";
}
