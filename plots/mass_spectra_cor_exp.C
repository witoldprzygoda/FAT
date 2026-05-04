// mass_spectra_cor_exp.C — Dilepton invariant mass: no OA cut vs OA > 4 deg
//                          (CORRECTED, exp data). OA > 4 is the active
//                          analysis cut (oa_pass==1 flag from main.cc).
//
// Reads from `dilepton_nt_cor` — same field layout as `dilepton_nt` but `m_ee`
// is the energy-loss-corrected invariant mass.
//
// Each spectrum is saved in both log and linear Y variants. For LOG plots
// the OA-cut panel reuses the no-OA Y range so the conversion-peak
// suppression is directly visible. For LINEAR plots the OA-cut panel
// gets its own auto-range — otherwise the conversion-peak in the no-OA
// reference would set the cap too high and the cut spectrum would be
// flattened to invisibility.
//
// Usage: root -l -b -q plots/mass_spectra_cor_exp.C

#include "PlotUtils.h"
#include <string>

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

    // Binning: 280 bins over [0, 1.4] → 5 MeV/bin (same width as before).
    const int    nbins = 280;
    const double xmin = 0.0;
    const double xmax = 1.4;

    // --- Build histograms for both cuts ---
    TH1D *a1, *c1, *s1, *a2, *c2, *s2;
    std::tie(a1, c1, s1) = pu.drawSignal("dilepton_nt_cor", "m_ee",
                                         nbins, xmin, xmax, "",
                                         ";M_{e^{+}e^{-}} [GeV/c^{2}];Counts");
    std::tie(a2, c2, s2) = pu.drawSignal("dilepton_nt_cor", "m_ee",
                                         nbins, xmin, xmax, "oa_pass==1",
                                         ";M_{e^{+}e^{-}} [GeV/c^{2}];Counts");

    const char* t_no = "M_{e^{+}e^{-}} (cor, no OA cut)";
    const char* t_oa = "M_{e^{+}e^{-}} (cor, OA > 4#circ, active)";

    // === No-OA: log (sets reference Y), then linear ===
    auto* cv_no_log = pu.drawTriple(a1, c1, s1, t_no, "c_mass_no_oa_cor_exp_log", /*logy=*/true);
    pu.save(cv_no_log, "mass_ee_no_oa_cor_exp_log");
    const double ymax_log = a1->GetMaximum();
    const double ymin_log = a1->GetMinimum();
    printIntegrals_cor("No OA cut (cor, log)", a1, c1, s1);

    // Linear Y: PlotUtils::drawTriple multiplies the histogram's CURRENT
    // GetMaximum() by 1.2 — but after the previous (log) call the max is
    // already inflated 3×, so its second multiply gives an absurd Y cap.
    // Re-cap manually from raw bin contents.
    const double a1_data_max = a1->GetBinContent(a1->GetMaximumBin());
    auto* cv_no_lin = pu.drawTriple(a1, c1, s1, t_no, "c_mass_no_oa_cor_exp_lin", /*logy=*/false);
    a1->SetMaximum(a1_data_max * 1.2); a1->SetMinimum(0.0);
    cv_no_lin->Update();
    pu.save(cv_no_lin, "mass_ee_no_oa_cor_exp_lin");

    // === OA cut: log (Y matched to no-OA log), then linear (own auto-range) ===
    auto* cv_oa_log = pu.drawTriple(a2, c2, s2, t_oa, "c_mass_oa4_cor_exp_log", /*logy=*/true);
    a2->SetMaximum(ymax_log); a2->SetMinimum(ymin_log);
    cv_oa_log->Update();
    pu.save(cv_oa_log, "mass_ee_oa4_cor_exp_log");
    printIntegrals_cor("OA > 4 deg (cor, active, log)", a2, c2, s2);

    const double a2_data_max = a2->GetBinContent(a2->GetMaximumBin());
    auto* cv_oa_lin = pu.drawTriple(a2, c2, s2, t_oa, "c_mass_oa4_cor_exp_lin", /*logy=*/false);
    a2->SetMaximum(a2_data_max * 1.2); a2->SetMinimum(0.0);
    cv_oa_lin->Update();
    pu.save(cv_oa_lin, "mass_ee_oa4_cor_exp_lin");

    std::cout << "\nDone. Check plots/output/\n";
}
