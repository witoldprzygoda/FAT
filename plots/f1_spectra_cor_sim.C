// f1_spectra_cor_sim.C — f1(1285) -> pi+pi- eta search on CORRECTED kinematics.
// Simulation mode (single file, no CB).
//
// Two reconstruction modes for the eta candidate:
//
//   (A) Dalitz eta: ECAL N_gamma==1, M(e+e-gamma) in [0.5, 0.6] GeV
//       eta -> e+e-gamma  (BR ~7e-3)
//       f1 candidate mass = M(pi+pi- e+e-gamma) = m_pippimepemg
//       Cut flag: eta_dalitz_pass==1
//
//   (B) gamma gamma eta: ECAL N_gamma==2, M(gamma gamma) in [0.5, 0.6] GeV
//       eta -> gamma gamma  (BR ~39%)
//       f1 candidate mass = M(pi+pi- gamma gamma) = m_pippim_gg
//       Cut flag: eta_gg_pass==1
//
// Plots produced (each as .pdf + .png, both log and linear Y):
//   m_epemg_cor_sim                — control: M(e+e-gamma) for ECAL N_gamma==1
//   m_gg_cor_sim                   — control: M(gamma gamma) for ECAL N_gamma==2
//   f1_pippim_eta_dalitz_cor_sim   — M(pi+pi-e+e-gamma) under eta_dalitz_pass==1
//   f1_pippim_eta_gg_cor_sim       — M(pi+pi-gg) under eta_gg_pass==1
//
// Reads single sim output from CWD:
//   output_pippimepem_sim.root  (no like-sign samples, no CB extraction)
//
// Usage:  root -l -b -q plots/f1_spectra_cor_sim.C

#include "PlotUtils.h"
#include <TH1.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <iostream>
#include <string>

// Helper: compose a weighted TTree::Draw cut. Per-event sim_genweight from
// the ntuple is multiplied by an optional boolean filter so all spectra
// reflect the SMASH luminosity normalisation.
namespace { std::string wcut(const std::string& filter = "") {
    return filter.empty() ? std::string("sim_genweight")
                          : "(" + filter + ")*sim_genweight";
}}


namespace {

void printIntegral(const char* label, TH1D* h) {
    std::cout << "  " << label
              << "  entries=" << (h ? h->GetEntries() : 0)
              << "  integral=" << (h ? h->Integral() : 0) << "\n";
}

// Force Y-axis: log -> 3 * data_max,  linear -> 1.2 * data_max.
// drawSingle's auto-cap is buggy when h has been SetMaximum'ed by a previous call.
void forceYCap(TH1D* h, double data_max, TCanvas* cv, bool logy) {
    h->SetMaximum(data_max * (logy ? 3.0 : 1.2));
    h->SetMinimum(logy ? 0.5 : 0.0);
    cv->Update();
}

// drawNtupleSingle + drawSingle (log + optional linear), with proper Y-cap.
void plotSingle(PlotUtils& pu,
                const char* nt_name, const char* varexpr,
                int nbins, double xmin, double xmax,
                const std::string& cut_filter,
                const char* title_with_axes,
                const char* basename,
                bool also_linear = false,
                bool linear_only = false)
{
    TH1D* h = pu.drawNtupleSingle(nt_name, varexpr, nbins, xmin, xmax,
                                   wcut(cut_filter), "");
    if (!h) {
        std::cerr << "plotSingle: drawNtupleSingle failed for " << basename << "\n";
        return;
    }

    double data_max = h->GetMaximum();

    if (linear_only) {
        auto* cv = pu.drawSingle(h, title_with_axes,
                                 (std::string("c_") + basename).c_str(),
                                 /*logy=*/false);
        forceYCap(h, data_max, cv, /*logy=*/false);
        pu.save(cv, basename);
    } else {
        auto* cv_log = pu.drawSingle(h, title_with_axes,
                                     (std::string("c_") + basename + "_log").c_str(),
                                     /*logy=*/true);
        forceYCap(h, data_max, cv_log, /*logy=*/true);
        pu.save(cv_log, also_linear ? (std::string(basename) + "_log").c_str() : basename);

        if (also_linear) {
            auto* cv_lin = pu.drawSingle(h, title_with_axes,
                                         (std::string("c_") + basename + "_lin").c_str(),
                                         /*logy=*/false);
            forceYCap(h, data_max, cv_lin, /*logy=*/false);
            pu.save(cv_lin, (std::string(basename) + "_lin").c_str());
        }
    }
    printIntegral(basename, h);
}

}  // namespace

void f1_spectra_cor_sim() {

    PlotUtils pu("output_pippimepem_sim.root");   // single-file (sim, no CB)

    std::cout << "\n=== Control plots (CORRECTED, sim) ===\n";

    // -- Control: M(e+e-gamma) for ECAL N_gamma==1 (eta candidate from Dalitz)
    plotSingle(pu, "pippimepem_nt_cor", "m_epemg",
               200, 0.0, 1.0,
               "m_epemg>0",   // -1 default for events without mult==1 epemg
               "M(e^{+}e^{-}#gamma)  (ECAL N_{#gamma}=1, control);"
               "M_{e^{+}e^{-}#gamma} [GeV/c^{2}];a.u.",
               "m_epemg_cor_sim",
               /*also_linear=*/true);

    // -- Control: M(gamma gamma) for ECAL N_gamma==2 (eta candidate from gg)
    plotSingle(pu, "pippimepem_nt_cor", "m_gg",
               200, 0.0, 1.0,
               "m_gg>0",   // -1 default for events without mult==2 gg
               "M(#gamma#gamma)  (ECAL N_{#gamma}=2, control);"
               "M_{#gamma#gamma} [GeV/c^{2}];a.u.",
               "m_gg_cor_sim",
               /*also_linear=*/true);

    std::cout << "\n=== f1(1285) candidates (CORRECTED, sim) ===\n";

    // -- f1 from Dalitz eta: M(pi+pi-e+e-gamma) under eta cut on M(epemg)
    plotSingle(pu, "pippimepem_nt_cor", "m_pippimepemg",
               360, 0.0, 1.8, "eta_dalitz_pass==1",
               "f_{1}(1285) #rightarrow #pi^{+}#pi^{-} #eta(e^{+}e^{-}#gamma)  "
               "[ECAL N_{#gamma}=1, M(e^{+}e^{-}#gamma) #in #eta];"
               "M_{#pi^{+}#pi^{-}e^{+}e^{-}#gamma} [GeV/c^{2}];a.u.",
               "f1_pippim_eta_dalitz_cor_sim",
               /*also_linear=*/true);

    // -- f1 from gamma gamma eta: M(pi+pi-gg) under eta cut on M(gg)
    plotSingle(pu, "pippimepem_nt_cor", "m_pippim_gg",
               360, 0.0, 1.8, "eta_gg_pass==1",
               "f_{1}(1285) #rightarrow #pi^{+}#pi^{-} #eta(#gamma#gamma)  "
               "[ECAL N_{#gamma}=2, M(#gamma#gamma) #in #eta];"
               "M_{#pi^{+}#pi^{-}#gamma#gamma} [GeV/c^{2}];a.u.",
               "f1_pippim_eta_gg_cor_sim",
               /*also_linear=*/true);

    std::cout << "\nDone. Plots saved under plots/output/\n";
}
