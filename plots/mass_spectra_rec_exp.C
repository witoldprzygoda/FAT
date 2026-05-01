// mass_spectra_rec_exp.C — Dilepton invariant mass: no OA cut vs OA > 9 deg (RECONSTRUCTED, exp data)
// Usage: root -l -b -q plots/mass_spectra_rec_exp.C

#include "PlotUtils.h"

void printIntegrals(const char* label, TH1D* all, TH1D* cb, TH1D* sig) {
    // Full range
    double i_all = all->Integral();
    double i_cb  = cb->Integral();
    double i_sig = sig->Integral();

    // Above pi0: M > 0.14 GeV/c^2
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

void mass_spectra_rec_exp() {

    PlotUtils pu("output_pippimepem_exp.root",
                 "output_pippimepep_exp.root",
                 "output_pippimemem_exp.root");

    // --- 1. Mass spectrum without OA cut ---
    TH1D *all1, *cb1, *sig1;
    std::tie(all1, cb1, sig1) = pu.drawSignal("pippimepem_nt", "m_ee",
                                              160, 0, 0.8, "",
                                              ";M_{e^{+}e^{-}} [GeV/c^{2}];Counts");
    auto* c1 = pu.drawTriple(all1, cb1, sig1,
                             "M_{e^{+}e^{-}} (no OA cut)", "c_mass_no_oa",
                             /*logy=*/true);
    pu.save(c1, "mass_ee_no_oa_rec_exp");
    printIntegrals("No OA cut", all1, cb1, sig1);

    // Capture Y-axis range from first plot
    double ymax = all1->GetMaximum();
    double ymin = all1->GetMinimum();

    // --- 2. Mass spectrum with OA > 9 deg cut (same Y range) ---
    TH1D *all2, *cb2, *sig2;
    std::tie(all2, cb2, sig2) = pu.drawSignal("pippimepem_nt", "m_ee",
                                              160, 0, 0.8, "oa>9",
                                              ";M_{e^{+}e^{-}} [GeV/c^{2}];Counts");
    auto* c2 = pu.drawTriple(all2, cb2, sig2,
                             "M_{e^{+}e^{-}} (OA > 9#circ)", "c_mass_oa9",
                             /*logy=*/true);
    all2->SetMaximum(ymax);
    all2->SetMinimum(ymin);
    c2->Update();
    pu.save(c2, "mass_ee_oa9_rec_exp");
    printIntegrals("OA > 9 deg", all2, cb2, sig2);

    std::cout << "\nDone. Check plots/output/\n";
}
