// mass_spectra_cor_exp.C — Dilepton invariant mass: no OA cut vs OA > 9 deg
//                          (CORRECTED, exp data).
//
// Reads from `dilepton_nt_cor` — same field layout as `dilepton_nt` but `m_ee`
// is the energy-loss-corrected invariant mass.
//
// Usage: root -l -b -q plots/mass_spectra_cor_exp.C

#include "PlotUtils.h"

void printIntegrals_cor(const char* label, TH1D* all, TH1D* cb, TH1D* sig) {
    double i_all = all->Integral();
    double i_cb  = cb->Integral();
    double i_sig = sig->Integral();

    int bin_pi0 = all->FindBin(0.1401);
    int bin_max = all->GetNbinsX();
    double i_all_above = all->Integral(bin_pi0, bin_max);
    double i_cb_above  = cb->Integral(bin_pi0, bin_max);
    double i_sig_above = sig->Integral(bin_pi0, bin_max);

    std::cout << "\n=== " << label << " ===\n";
    std::cout << "  Full range:    all = " << i_all
              << "  CB = " << i_cb
              << "  sig = " << i_sig << "\n";
    std::cout << "  M > 0.14:      all = " << i_all_above
              << "  CB = " << i_cb_above
              << "  sig = " << i_sig_above << "\n";
}

void mass_spectra_cor_exp() {

    PlotUtils pu("output_epem_exp.root",
                 "output_epep_exp.root",
                 "output_emem_exp.root");

    // --- 1. Mass spectrum without OA cut ---
    TH1D *all1, *cb1, *sig1;
    std::tie(all1, cb1, sig1) = pu.drawSignal("dilepton_nt_cor", "m_ee",
                                              160, 0, 0.8, "",
                                              ";M_{e^{+}e^{-}} [GeV/c^{2}];Counts");
    auto* c1 = pu.drawTriple(all1, cb1, sig1,
                             "M_{e^{+}e^{-}} (cor, no OA cut)", "c_mass_no_oa_cor_exp",
                             /*logy=*/true);
    pu.save(c1, "mass_ee_no_oa_cor_exp");
    printIntegrals_cor("No OA cut (cor)", all1, cb1, sig1);

    double ymax = all1->GetMaximum();
    double ymin = all1->GetMinimum();

    // --- 2. Mass spectrum with OA > 9 deg cut (same Y range) ---
    TH1D *all2, *cb2, *sig2;
    std::tie(all2, cb2, sig2) = pu.drawSignal("dilepton_nt_cor", "m_ee",
                                              160, 0, 0.8, "oa>9",
                                              ";M_{e^{+}e^{-}} [GeV/c^{2}];Counts");
    auto* c2 = pu.drawTriple(all2, cb2, sig2,
                             "M_{e^{+}e^{-}} (cor, OA > 9#circ)", "c_mass_oa9_cor_exp",
                             /*logy=*/true);
    all2->SetMaximum(ymax);
    all2->SetMinimum(ymin);
    c2->Update();
    pu.save(c2, "mass_ee_oa9_cor_exp");
    printIntegrals_cor("OA > 9 deg (cor)", all2, cb2, sig2);

    std::cout << "\nDone. Check plots/output/\n";
}
