// mm_pepem.C — Missing mass of pe+e- and M_ee with MM proton window
// Usage: root -l -b -q plots/mm_pepem.C

#include "PlotUtils.h"

void printIntegrals(const char* label, TH1D* all, TH1D* cb, TH1D* sig) {
    double i_all = all->Integral();
    double i_cb  = cb->Integral();
    double i_sig = sig->Integral();

    std::cout << "\n=== " << label << " ===\n";
    std::cout << "  Integral:  all = " << i_all
              << "  CB = " << i_cb
              << "  sig = " << i_sig << "\n";
}

void printIntegralsRange(const char* label, TH1D* all, TH1D* cb, TH1D* sig,
                         double xmin, double xmax) {
    int bin_lo = all->FindBin(xmin + 0.0001);
    int bin_hi = all->FindBin(xmax - 0.0001);
    double i_all = all->Integral(bin_lo, bin_hi);
    double i_cb  = cb->Integral(bin_lo, bin_hi);
    double i_sig = sig->Integral(bin_lo, bin_hi);

    std::cout << "  " << label << ":  all = " << i_all
              << "  CB = " << i_cb
              << "  sig = " << i_sig << "\n";
}

void mm_pepem() {

    PlotUtils pu("output_pepem.root", "output_pepep.root", "output_pemem.root");

    // --- 1. Missing mass of pe+e- (OA > 9, M_ee > 0.14) ---
    TH1D *all1, *cb1, *sig1;
    std::tie(all1, cb1, sig1) = pu.drawSignal("dilepton_nt", "mm_pepem_mass",
                                                200, 0, 2.0, "oa>9 && m_ee>0.14",
                                                ";MM(pe^{+}e^{-}) [GeV/c^{2}];Counts");

    auto* c1 = pu.drawTriple(all1, cb1, sig1,
                              "MM(pe^{+}e^{-}) (OA > 9#circ, M_{ee} > 0.14)",
                              "c_mm_pepem");

    // Linear Y-axis, fixed ranges
    c1->SetLogy(0);
    all1->GetXaxis()->SetRangeUser(0.6, 1.5);
    all1->SetMaximum(200);
    all1->SetMinimum(0);
    c1->Update();

    pu.save(c1, "mm_pepem");
    printIntegrals("MM(pe+e-) full range", all1, cb1, sig1);
    printIntegralsRange("MM 0.88-1.02", all1, cb1, sig1, 0.88, 1.02);

    // --- 2. M_ee with OA > 9 and MM(pe+e-) proton window ---
    TH1D *all2, *cb2, *sig2;
    std::tie(all2, cb2, sig2) = pu.drawSignal("dilepton_nt", "m_ee",
                                                160, 0, 0.8,
                                                "oa>9 && mm_pepem_mass>=0.88 && mm_pepem_mass<=1.02",
                                                ";M_{e^{+}e^{-}} [GeV/c^{2}];Counts");

    auto* c2 = pu.drawTriple(all2, cb2, sig2,
                              "M_{e^{+}e^{-}} (OA > 9#circ, MM proton)",
                              "c_mee_mm_pepem");
    pu.save(c2, "m_ee_mm_pepem");

    printIntegrals("M_ee (OA>9, MM proton) full range", all2, cb2, sig2);

    int bin_pi0 = all2->FindBin(0.1401);
    int bin_max = all2->GetNbinsX();
    double i_all_above = all2->Integral(bin_pi0, bin_max);
    double i_cb_above  = cb2->Integral(bin_pi0, bin_max);
    double i_sig_above = sig2->Integral(bin_pi0, bin_max);
    std::cout << "  M > 0.14:      all = " << i_all_above
              << "  CB = " << i_cb_above
              << "  sig = " << i_sig_above << "\n";

    // --- 3. Invariant mass of pe+e- (OA > 9, MM proton window) ---
    TH1D *all3, *cb3, *sig3;
    std::tie(all3, cb3, sig3) = pu.drawSignal("dilepton_nt", "pepem_mass",
                                                100, 0.8, 1.8,
                                                "oa>9 && m_ee>0.14 && mm_pepem_mass>=0.88 && mm_pepem_mass<=1.02",
                                                ";M(pe^{+}e^{-}) [GeV/c^{2}];Counts");

    auto* c3 = pu.drawTriple(all3, cb3, sig3,
                              "M(pe^{+}e^{-}) (OA > 9#circ, M_{ee} > 0.14, MM proton)",
                              "c_pepem_inv_mass");
    c3->SetLogy(0);
    all3->SetMaximum(40);
    all3->SetMinimum(0);
    c3->Update();
    pu.save(c3, "pepem_inv_mass");

    printIntegrals("M(pe+e-) (OA>9, MM proton) full range", all3, cb3, sig3);

    std::cout << "\nDone. Check plots/output/\n";
}
