// pippimepem_spectra.C — all / CB / signal for the pi+pi-e+e- observables:
//   M(e+e-), M(pi+pi-), M(pi+pi-e+e-),
//   MM(pi+pi-), MM(pi+pi-e+e-),
//   M(pi+pi-e+e-) after the pippimepem_selection cut chain
//
// CB = 2 * sqrt(N++ * N--) reconstructed from the like-sign samples
// (output_pippimepep.root, output_pippimemem.root). Signal = all - CB.
//
// Usage: root -l -b -q plots/pippimepem_spectra.C

#include "PlotUtils.h"
#include <sstream>

void printIntegrals(const char* label, TH1D* all, TH1D* cb, TH1D* sig) {
    double i_all = all->Integral();
    double i_cb  = cb->Integral();
    double i_sig = sig->Integral();

    std::cout << "\n=== " << label << " ===\n";
    std::cout << "  all = " << i_all
              << "   CB = " << i_cb
              << "   sig = " << i_sig << "\n";
}

void pippimepem_spectra() {

    PlotUtils pu("output_pippimepem.root",
                 "output_pippimepep.root",
                 "output_pippimemem.root");

    // --- 1. M(e+e-) -------------------------------------------------------
    TH1D *a1, *c1, *s1;
    std::tie(a1, c1, s1) = pu.drawSignal(
        "pippimepem_nt", "m_ee",
        160, 0.0, 0.8, "",
        ";M_{e^{+}e^{-}} [GeV/c^{2}];Counts");
    auto* cv1 = pu.drawTriple(a1, c1, s1, "M_{e^{+}e^{-}}", "c_m_ee", /*logy=*/true);
    pu.save(cv1, "m_ee");
    printIntegrals("M(e+e-)", a1, c1, s1);

    // --- 2. M(pi+pi-) -----------------------------------------------------
    TH1D *a2, *c2, *s2;
    std::tie(a2, c2, s2) = pu.drawSignal(
        "pippimepem_nt", "m_pippim",
        200, 0.0, 2.0, "",
        ";M_{#pi^{+}#pi^{-}} [GeV/c^{2}];Counts");
    auto* cv2 = pu.drawTriple(a2, c2, s2, "M_{#pi^{+}#pi^{-}}", "c_m_pippim");
    pu.save(cv2, "m_pippim");
    printIntegrals("M(pi+pi-)", a2, c2, s2);

    // --- 3. M(pi+pi-e+e-) -------------------------------------------------
    TH1D *a3, *c3, *s3;
    std::tie(a3, c3, s3) = pu.drawSignal(
        "pippimepem_nt", "m_pippimepem",
        200, 0.0, 2.0, "",
        ";M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}];Counts");
    auto* cv3 = pu.drawTriple(a3, c3, s3,
                              "M_{#pi^{+}#pi^{-}e^{+}e^{-}}", "c_m_pippimepem");
    pu.save(cv3, "m_pippimepem");
    printIntegrals("M(pi+pi-e+e-)", a3, c3, s3);

    // --- 4. MM(pi+pi-) ----------------------------------------------------
    TH1D *a4, *c4, *s4;
    std::tie(a4, c4, s4) = pu.drawSignal(
        "pippimepem_nt", "mm_pippim",
        200, 0.0, 4.0, "",
        ";MM(#pi^{+}#pi^{-}) [GeV/c^{2}];Counts");
    auto* cv4 = pu.drawTriple(a4, c4, s4,
                              "MM(#pi^{+}#pi^{-})", "c_mm_pippim");
    pu.save(cv4, "mm_pippim");
    printIntegrals("MM(pi+pi-)", a4, c4, s4);

    // --- 5. MM(pi+pi-e+e-) ------------------------------------------------
    TH1D *a5, *c5, *s5;
    std::tie(a5, c5, s5) = pu.drawSignal(
        "pippimepem_nt", "mm_pippimepem",
        200, 0.0, 4.0, "",
        ";MM(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}];Counts");
    auto* cv5 = pu.drawTriple(a5, c5, s5,
                              "MM(#pi^{+}#pi^{-}e^{+}e^{-})", "c_mm_pippimepem");
    pu.save(cv5, "mm_pippimepem");
    printIntegrals("MM(pi+pi-e+e-)", a5, c5, s5);

    // --- 6. M(pi+pi-e+e-) after pippimepem_selection cut chain ------------
    // Selection: OA_LAB((pippim),(epem)) < 50 deg
    //            M(pi+pi-) < 0.420 GeV/c^2
    //            OA_rest((pippim),(epem)) > 140 deg
    TH1D *a6, *c6, *s6;
    std::tie(a6, c6, s6) = pu.drawSignal(
        "pippimepem_nt", "m_pippimepem",
        200, 0.0, 2.0, "sel_pass==1",
        ";M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}];Counts");
    auto* cv6 = pu.drawTriple(a6, c6, s6,
                              "M_{#pi^{+}#pi^{-}e^{+}e^{-}} (after selection)",
                              "c_m_pippimepem_selected");
    pu.save(cv6, "m_pippimepem_selected");
    printIntegrals("M(pi+pi-e+e-) after selection", a6, c6, s6);

    // --- 7-12. M(pi+pi-e+e-) in MM(pi+pi-e+e-) slice windows --------------
    // Each slice plot applies: sel_pass==1 && lo <= mm_pippimepem <= hi
    // (matches the in-code RangeCut semantics, inclusive on both endpoints).
    auto drawSlice = [&](double lo, double hi, const std::string& tag) {
        std::ostringstream cut;
        cut << "sel_pass==1 && mm_pippimepem>=" << lo
            << " && mm_pippimepem<=" << hi;
        std::ostringstream title;
        title << "M_{#pi^{+}#pi^{-}e^{+}e^{-}} sel, MM #in [" << lo << ", " << hi << "] GeV/c^{2}";

        TH1D *a, *c, *s;
        std::tie(a, c, s) = pu.drawSignal(
            "pippimepem_nt", "m_pippimepem",
            100, 0.0, 2.0, cut.str(),
            ";M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}];Counts");
        auto* cv = pu.drawTriple(a, c, s,
                                  title.str(),
                                  "c_m_pippimepem_slice_" + tag);
        pu.save(cv, "m_pippimepem_slice_" + tag);

        std::ostringstream lbl;
        lbl << "M(pi+pi-e+e-) sel, MM in [" << lo << ", " << hi << "]";
        printIntegrals(lbl.str().c_str(), a, c, s);
    };

    drawSlice(2.0, 2.2, "20_22");
    drawSlice(2.2, 2.4, "22_24");
    drawSlice(2.4, 2.6, "24_26");
    drawSlice(2.6, 2.8, "26_28");
    drawSlice(2.8, 3.0, "28_30");

    std::cout << "\nDone. Check plots/output/\n";
}
