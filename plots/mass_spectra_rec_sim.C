// mass_spectra_rec_sim.C — Dilepton invariant mass: no OA cut vs opening_angle_4
//                          Reads RECONSTRUCTED kinematics from dilepton_nt.
// Simulation mode (single file, no CB extraction).
// Usage: root -l -b -q plots/mass_spectra_rec_sim.C

#include "PlotUtils.h"

// Helper: compose a weighted TTree::Draw cut. Per-event sim_genweight from
// the ntuple is multiplied by an optional boolean filter so all spectra
// reflect the SMASH luminosity normalisation.
namespace { std::string wcut(const std::string& filter = "") {
    return filter.empty() ? std::string("sim_genweight")
                          : "(" + filter + ")*sim_genweight";
}}

void printIntegrals(const char* label, TH1D* h) {
    double i_full = h->Integral();
    int bin_pi0 = h->FindBin(0.1401);
    double i_above = h->Integral(bin_pi0, h->GetNbinsX());

    std::cout << "\n=== " << label << " ===\n";
    std::cout << "  Full range: " << i_full
              << "    M > 0.14 GeV/c^2: " << i_above << "\n";
}

void mass_spectra_rec_sim() {

    PlotUtils pu("output_epem_sim.root");   // single-file (sim, no CB)

    // --- 1. Mass spectrum without OA cut ---
    TH1D* h1 = pu.drawNtupleSingle("dilepton_nt", "m_ee",
                                   160, 0, 0.8, wcut(""),
                                   ";M_{e^{+}e^{-}} [GeV/c^{2}];a.u.");
    auto* c1 = pu.drawSingle(h1,
                             "M_{e^{+}e^{-}} (no OA cut)", "c_mass_no_oa_rec_sim",
                             /*logy=*/true);
    pu.save(c1, "mass_ee_no_oa_rec_sim");
    printIntegrals("No OA cut", h1);

    double ymax = h1->GetMaximum();
    double ymin = h1->GetMinimum();

    // --- 2. Mass spectrum with opening_angle_4 (oa>4) cut applied (same Y range) ---
    TH1D* h2 = pu.drawNtupleSingle("dilepton_nt", "m_ee",
                                   160, 0, 0.8, wcut("oa_pass==1"),
                                   ";M_{e^{+}e^{-}} [GeV/c^{2}];a.u.");
    auto* c2 = pu.drawSingle(h2,
                             "M_{e^{+}e^{-}} (OA > 4#circ)", "c_mass_oa4_rec_sim",
                             /*logy=*/true);
    h2->SetMaximum(ymax);
    h2->SetMinimum(ymin);
    c2->Update();
    pu.save(c2, "mass_ee_oa4_rec_sim");
    printIntegrals("OA > 4 deg", h2);

    std::cout << "\nDone. Check plots/output/\n";
}
