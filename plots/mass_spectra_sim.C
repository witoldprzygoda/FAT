// mass_spectra_sim.C — Dilepton invariant mass: no OA cut vs opening_angle_4 (>4 deg)
//                       Reads SIMULATED-truth values from pippimepem_nt_cor
//                       (fields m_ee_sim) and applies the per-event sim_genweight.
// Single-file mode (no CB extraction — simulation has no like-sign sample).
//
// Usage: root -l -b -q plots/mass_spectra_sim.C

#include "PlotUtils.h"

void printIntegrals_sim(const char* label, TH1D* h) {
    double i_full = h->Integral();
    int bin_pi0 = h->FindBin(0.1401);
    double i_above = h->Integral(bin_pi0, h->GetNbinsX());

    std::cout << "\n=== " << label << " ===\n";
    std::cout << "  Full range: " << i_full
              << "    M > 0.14 GeV/c^2: " << i_above << "\n";
}

void mass_spectra_sim() {

    PlotUtils pu("output_pippimepem_sim.root");

    // --- 1. Mass spectrum (sim) without OA cut ---
    TH1D* h1 = pu.drawNtupleSingle("pippimepem_nt_cor", "m_ee_sim",
                                   160, 0, 0.8, "sim_genweight",
                                   ";M_{e^{+}e^{-}} [GeV/c^{2}] (sim);a.u.");
    auto* c1 = pu.drawSingle(h1,
                             "M_{e^{+}e^{-}} (sim, no OA cut)", "c_mass_no_oa_sim",
                             /*logy=*/true);
    pu.save(c1, "mass_ee_no_oa_sim");
    printIntegrals_sim("No OA cut (sim)", h1);

    double ymax = h1->GetMaximum();
    double ymin = h1->GetMinimum();

    // --- 2. Mass spectrum (sim) with opening_angle_4 (oa>4) cut applied ---
    // oa_pass==1 is the boolean from passMinCut("opening_angle_4", oa) in main.cc.
    TH1D* h2 = pu.drawNtupleSingle("pippimepem_nt_cor", "m_ee_sim",
                                   160, 0, 0.8, "(oa_pass==1)*sim_genweight",
                                   ";M_{e^{+}e^{-}} [GeV/c^{2}] (sim);a.u.");
    auto* c2 = pu.drawSingle(h2,
                             "M_{e^{+}e^{-}} (sim, OA > 4#circ)", "c_mass_oa4_sim",
                             /*logy=*/true);
    h2->SetMaximum(ymax);
    h2->SetMinimum(ymin);
    c2->Update();
    pu.save(c2, "mass_ee_oa4_sim");
    printIntegrals_sim("OA > 4 deg (sim)", h2);

    std::cout << "\nDone. Plots in plots/output/ (suffix _sim)\n";
}
