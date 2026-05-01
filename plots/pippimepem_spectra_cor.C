// pippimepem_spectra_cor.C — CORRECTED-kinematics counterpart of pippimepem_spectra.C.
//
// Reads from `pippimepem_nt_cor` (mirror ntuple where compound observables are
// computed from KinematicType::CORRECTED). Field names match the RECONSTRUCTED
// ntuple, so cut expressions stay identical; only the underlying values differ.
//
// All output filenames are suffixed with "_cor" to keep the figures from the
// two macros side by side in plots/output/ without overwrite collisions.
//
// Usage: root -l -b -q plots/pippimepem_spectra_cor.C

#include "PlotUtils.h"
#include <sstream>

void printIntegrals_cor(const char* label, TH1D* all, TH1D* cb, TH1D* sig) {
    double i_all = all->Integral();
    double i_cb  = cb->Integral();
    double i_sig = sig->Integral();

    std::cout << "\n=== " << label << " ===\n";
    std::cout << "  all = " << i_all
              << "   CB = " << i_cb
              << "   sig = " << i_sig << "\n";
}

void pippimepem_spectra_cor() {

    PlotUtils pu("output_pippimepem_exp.root",
                 "output_pippimepep_exp.root",
                 "output_pippimemem_exp.root");

    const char* NT = "pippimepem_nt_cor";

    // --- 1. M(e+e-) -------------------------------------------------------
    TH1D *a1, *c1, *s1;
    std::tie(a1, c1, s1) = pu.drawSignal(
        NT, "m_ee",
        160, 0.0, 0.8, "",
        ";M_{e^{+}e^{-}} [GeV/c^{2}];Counts");
    auto* cv1 = pu.drawTriple(a1, c1, s1, "M_{e^{+}e^{-}} (cor)", "c_m_ee_cor", /*logy=*/true);
    pu.save(cv1, "m_ee_cor");
    printIntegrals_cor("M(e+e-) cor", a1, c1, s1);

    // --- 2. M(pi+pi-) -----------------------------------------------------
    TH1D *a2, *c2, *s2;
    std::tie(a2, c2, s2) = pu.drawSignal(
        NT, "m_pippim",
        200, 0.0, 2.0, "",
        ";M_{#pi^{+}#pi^{-}} [GeV/c^{2}];Counts");
    auto* cv2 = pu.drawTriple(a2, c2, s2, "M_{#pi^{+}#pi^{-}} (cor)", "c_m_pippim_cor");
    pu.save(cv2, "m_pippim_cor");
    printIntegrals_cor("M(pi+pi-) cor", a2, c2, s2);

    // --- 3. M(pi+pi-e+e-) -------------------------------------------------
    TH1D *a3, *c3, *s3;
    std::tie(a3, c3, s3) = pu.drawSignal(
        NT, "m_pippimepem",
        200, 0.2, 1.4, "",
        ";M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}];Counts");
    auto* cv3 = pu.drawTriple(a3, c3, s3,
                              "M_{#pi^{+}#pi^{-}e^{+}e^{-}} (cor)", "c_m_pippimepem_cor");
    pu.save(cv3, "m_pippimepem_cor");
    printIntegrals_cor("M(pi+pi-e+e-) cor", a3, c3, s3);

    // --- 4. MM(pi+pi-) ----------------------------------------------------
    TH1D *a4, *c4, *s4;
    std::tie(a4, c4, s4) = pu.drawSignal(
        NT, "mm_pippim",
        200, 0.0, 4.0, "",
        ";MM(#pi^{+}#pi^{-}) [GeV/c^{2}];Counts");
    auto* cv4 = pu.drawTriple(a4, c4, s4,
                              "MM(#pi^{+}#pi^{-}) (cor)", "c_mm_pippim_cor");
    pu.save(cv4, "mm_pippim_cor");
    printIntegrals_cor("MM(pi+pi-) cor", a4, c4, s4);

    // --- 5. MM(pi+pi-e+e-) ------------------------------------------------
    TH1D *a5, *c5, *s5;
    std::tie(a5, c5, s5) = pu.drawSignal(
        NT, "mm_pippimepem",
        200, 0.0, 4.0, "",
        ";MM(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}];Counts");
    auto* cv5 = pu.drawTriple(a5, c5, s5,
                              "MM(#pi^{+}#pi^{-}e^{+}e^{-}) (cor)", "c_mm_pippimepem_cor");
    pu.save(cv5, "mm_pippimepem_cor");
    printIntegrals_cor("MM(pi+pi-e+e-) cor", a5, c5, s5);

    // --- 5a-5c. MM observables with the 2D graphical cut (cut_2d) ---------
    {
        TH1D *a, *c, *s;
        std::tie(a, c, s) = pu.drawSignal(
            NT, "mm_epem",
            200, 0.0, 4.0, "cut2d_pass==1",
            ";MM(e^{+}e^{-}) [GeV/c^{2}];Counts");
        auto* cv = pu.drawTriple(a, c, s,
                                 "MM(e^{+}e^{-}) after cut_2d (cor)", "c_mm_epem_cut2d_cor");
        pu.save(cv, "mm_epem_cut2d_cor");
        printIntegrals_cor("MM(e+e-) after cut_2d cor", a, c, s);
    }
    {
        TH1D *a, *c, *s;
        std::tie(a, c, s) = pu.drawSignal(
            NT, "mm_pippim",
            200, 0.0, 4.0, "cut2d_pass==1",
            ";MM(#pi^{+}#pi^{-}) [GeV/c^{2}];Counts");
        auto* cv = pu.drawTriple(a, c, s,
                                 "MM(#pi^{+}#pi^{-}) after cut_2d (cor)", "c_mm_pippim_cut2d_cor");
        pu.save(cv, "mm_pippim_cut2d_cor");
        printIntegrals_cor("MM(pi+pi-) after cut_2d cor", a, c, s);
    }
    {
        TH1D *a, *c, *s;
        std::tie(a, c, s) = pu.drawSignal(
            NT, "mm_pippimepem",
            200, 0.0, 4.0, "cut2d_pass==1",
            ";MM(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}];Counts");
        auto* cv = pu.drawTriple(a, c, s,
                                 "MM(#pi^{+}#pi^{-}e^{+}e^{-}) after cut_2d (cor)", "c_mm_pippimepem_cut2d_cor");
        pu.save(cv, "mm_pippimepem_cut2d_cor");
        printIntegrals_cor("MM(pi+pi-e+e-) after cut_2d cor", a, c, s);
    }

    // --- 5d. OA observables driving the pippimepem_selection cut chain ----
    {
        TH1D *a, *c, *s;
        std::tie(a, c, s) = pu.drawSignal(
            NT, "oa_pippim_epem_lab",
            180, 0.0, 180.0, "",
            ";OA_{LAB}((#pi^{+}#pi^{-}),(e^{+}e^{-})) [deg];Counts");
        auto* cv = pu.drawTriple(a, c, s,
                                 "OA_{LAB}((#pi^{+}#pi^{-}),(e^{+}e^{-})) (cor)",
                                 "c_oa_pippim_epem_lab_cor");
        pu.save(cv, "oa_pippim_epem_lab_cor");
        printIntegrals_cor("OA_LAB((pi+pi-),(e+e-)) cor", a, c, s);
    }
    {
        TH1D *a, *c, *s;
        std::tie(a, c, s) = pu.drawSignal(
            NT, "oa_pippim_epem_eta_rest",
            180, 0.0, 180.0, "",
            ";OA_{m_{#eta}-rest}((#pi^{+}#pi^{-}),(e^{+}e^{-})) [deg];Counts");
        auto* cv = pu.drawTriple(a, c, s,
                                 "OA((#pi^{+}#pi^{-}),(e^{+}e^{-})) in m_{#eta}-rest frame (cor)",
                                 "c_oa_pippim_epem_eta_rest_cor");
        pu.save(cv, "oa_pippim_epem_eta_rest_cor");
        printIntegrals_cor("OA_eta_rest((pi+pi-),(e+e-)) cor", a, c, s);
    }

    // --- 6. M(pi+pi-e+e-) after pippimepem_selection cut chain ------------
    TH1D *a6, *c6, *s6;
    std::tie(a6, c6, s6) = pu.drawSignal(
        NT, "m_pippimepem",
        200, 0.2, 1.4, "sel_pass==1",
        ";M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}];Counts");
    auto* cv6 = pu.drawTriple(a6, c6, s6,
                              "M_{#pi^{+}#pi^{-}e^{+}e^{-}} (after selection, cor)",
                              "c_m_pippimepem_selected_cor");
    pu.save(cv6, "m_pippimepem_selected_cor");
    printIntegrals_cor("M(pi+pi-e+e-) after selection cor", a6, c6, s6);

    // --- 6a. Same as plot 6 but with cut_2d additionally applied ----------
    {
        TH1D *a, *c, *s;
        std::tie(a, c, s) = pu.drawSignal(
            NT, "m_pippimepem",
            200, 0.2, 1.4, "sel_pass==1 && cut2d_pass==1",
            ";M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}];Counts");
        auto* cv = pu.drawTriple(a, c, s,
                                 "M_{#pi^{+}#pi^{-}e^{+}e^{-}} (after selection + cut_2d, cor)",
                                 "c_m_pippimepem_selected_cut2d_cor");
        pu.save(cv, "m_pippimepem_selected_cut2d_cor");
        printIntegrals_cor("M(pi+pi-e+e-) after selection + cut_2d cor", a, c, s);
    }

    // --- 7-12. M(pi+pi-e+e-) in MM(pi+pi-e+e-) slice windows --------------
    auto drawSlice = [&](double lo, double hi, const std::string& tag, bool cut2d = false) {
        std::ostringstream cut;
        cut << "sel_pass==1 && mm_pippimepem>=" << lo
            << " && mm_pippimepem<=" << hi;
        if (cut2d) cut << " && cut2d_pass==1";

        std::ostringstream title;
        title << "M_{#pi^{+}#pi^{-}e^{+}e^{-}} sel" << (cut2d ? "+cut_2d" : "")
              << ", MM #in [" << lo << ", " << hi << "] GeV/c^{2} (cor)";

        std::string suffix = tag + (cut2d ? "_cut2d" : "") + "_cor";

        TH1D *a, *c, *s;
        std::tie(a, c, s) = pu.drawSignal(
            NT, "m_pippimepem",
            100, 0.2, 1.4, cut.str(),
            ";M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}];Counts");
        auto* cv = pu.drawTriple(a, c, s,
                                  title.str(),
                                  "c_m_pippimepem_slice_" + suffix);
        pu.save(cv, "m_pippimepem_slice_" + suffix);

        std::ostringstream lbl;
        lbl << "M(pi+pi-e+e-) sel" << (cut2d ? "+cut_2d" : "")
            << ", MM in [" << lo << ", " << hi << "] cor";
        printIntegrals_cor(lbl.str().c_str(), a, c, s);
    };

    // Slice plots without 2D cut
    drawSlice(2.0, 2.2, "20_22");
    drawSlice(2.2, 2.4, "22_24");
    drawSlice(2.4, 2.6, "24_26");
    drawSlice(2.6, 2.8, "26_28");
    drawSlice(2.8, 3.0, "28_30");

    // Same slices, additionally with cut_2d
    drawSlice(2.0, 2.2, "20_22", true);
    drawSlice(2.2, 2.4, "22_24", true);
    drawSlice(2.4, 2.6, "24_26", true);
    drawSlice(2.6, 2.8, "26_28", true);
    drawSlice(2.8, 3.0, "28_30", true);

    std::cout << "\nDone. Check plots/output/ for *_cor.{pdf,png}\n";
}
