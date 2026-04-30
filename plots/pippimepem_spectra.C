// pippimepem_spectra.C — pi+pi-e+e- observables for SIMULATION (no CB).
// Single-file mode: reads only output_pippimepem.root and draws each curve
// with PlotUtils::drawSingle. Unlike the experimental analysis there is no
// CB extraction (no like-sign samples in simulation).
//
// Plots produced:
//   M(e+e-), M(pi+pi-), M(pi+pi-e+e-),
//   MM(pi+pi-), MM(pi+pi-e+e-),
//   MM observables after cut_2d,
//   M(pi+pi-e+e-) after pippimepem_selection chain,
//   M(pi+pi-e+e-) in MM slice windows (with and without cut_2d).
//
// Usage: root -l -b -q plots/pippimepem_spectra.C

#include "PlotUtils.h"
#include <sstream>

// Helper: compose a weighted TTree::Draw cut. Per-event sim_genweight from
// the ntuple is multiplied by an optional boolean filter so all spectra
// reflect the SMASH luminosity normalisation.
namespace { std::string wcut(const std::string& filter = "") {
    return filter.empty() ? std::string("sim_genweight")
                          : "(" + filter + ")*sim_genweight";
}}


void printIntegral(const char* label, TH1D* h) {
    std::cout << "\n=== " << label << " ===\n";
    std::cout << "  entries = " << h->GetEntries()
              << "   integral = " << h->Integral() << "\n";
}

void pippimepem_spectra() {

    PlotUtils pu("output_pippimepem.root");   // single-file (sim, no CB)

    // --- 1. M(e+e-) -------------------------------------------------------
    auto* h1 = pu.drawNtupleSingle(
        "pippimepem_nt", "m_ee",
        160, 0.0, 0.8, wcut(""),
        ";M_{e^{+}e^{-}} [GeV/c^{2}];a.u.");
    auto* cv1 = pu.drawSingle(h1, "M_{e^{+}e^{-}}", "c_m_ee", /*logy=*/true);
    pu.save(cv1, "m_ee");
    printIntegral("M(e+e-)", h1);

    // --- 2. M(pi+pi-) -----------------------------------------------------
    auto* h2 = pu.drawNtupleSingle(
        "pippimepem_nt", "m_pippim",
        200, 0.0, 2.0, wcut(""),
        ";M_{#pi^{+}#pi^{-}} [GeV/c^{2}];a.u.");
    auto* cv2 = pu.drawSingle(h2, "M_{#pi^{+}#pi^{-}}", "c_m_pippim");
    pu.save(cv2, "m_pippim");
    printIntegral("M(pi+pi-)", h2);

    // --- 3. M(pi+pi-e+e-) -------------------------------------------------
    auto* h3 = pu.drawNtupleSingle(
        "pippimepem_nt", "m_pippimepem",
        200, 0.2, 1.4, wcut(""),
        ";M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}];a.u.");
    auto* cv3 = pu.drawSingle(h3,
                              "M_{#pi^{+}#pi^{-}e^{+}e^{-}}", "c_m_pippimepem");
    pu.save(cv3, "m_pippimepem");
    printIntegral("M(pi+pi-e+e-)", h3);

    // --- 4. MM(pi+pi-) ----------------------------------------------------
    auto* h4 = pu.drawNtupleSingle(
        "pippimepem_nt", "mm_pippim",
        200, 0.0, 4.0, wcut(""),
        ";MM(#pi^{+}#pi^{-}) [GeV/c^{2}];a.u.");
    auto* cv4 = pu.drawSingle(h4, "MM(#pi^{+}#pi^{-})", "c_mm_pippim");
    pu.save(cv4, "mm_pippim");
    printIntegral("MM(pi+pi-)", h4);

    // --- 5. MM(pi+pi-e+e-) ------------------------------------------------
    auto* h5 = pu.drawNtupleSingle(
        "pippimepem_nt", "mm_pippimepem",
        200, 0.0, 4.0, wcut(""),
        ";MM(#pi^{+}#pi^{-}e^{+}e^{+}) [GeV/c^{2}];a.u.");
    auto* cv5 = pu.drawSingle(h5, "MM(#pi^{+}#pi^{-}e^{+}e^{-})",
                              "c_mm_pippimepem");
    pu.save(cv5, "mm_pippimepem");
    printIntegral("MM(pi+pi-e+e-)", h5);

    // --- 5a-5c. MM observables after cut_2d -------------------------------
    {
        auto* h = pu.drawNtupleSingle(
            "pippimepem_nt", "mm_epem",
            200, 0.0, 4.0, wcut("cut2d_pass==1"),
            ";MM(e^{+}e^{-}) [GeV/c^{2}];a.u.");
        auto* cv = pu.drawSingle(h, "MM(e^{+}e^{-}) after cut_2d", "c_mm_epem_cut2d");
        pu.save(cv, "mm_epem_cut2d");
        printIntegral("MM(e+e-) after cut_2d", h);
    }
    {
        auto* h = pu.drawNtupleSingle(
            "pippimepem_nt", "mm_pippim",
            200, 0.0, 4.0, wcut("cut2d_pass==1"),
            ";MM(#pi^{+}#pi^{-}) [GeV/c^{2}];a.u.");
        auto* cv = pu.drawSingle(h, "MM(#pi^{+}#pi^{-}) after cut_2d", "c_mm_pippim_cut2d");
        pu.save(cv, "mm_pippim_cut2d");
        printIntegral("MM(pi+pi-) after cut_2d", h);
    }
    {
        auto* h = pu.drawNtupleSingle(
            "pippimepem_nt", "mm_pippimepem",
            200, 0.0, 4.0, wcut("cut2d_pass==1"),
            ";MM(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}];a.u.");
        auto* cv = pu.drawSingle(h, "MM(#pi^{+}#pi^{-}e^{+}e^{-}) after cut_2d",
                                 "c_mm_pippimepem_cut2d");
        pu.save(cv, "mm_pippimepem_cut2d");
        printIntegral("MM(pi+pi-e+e-) after cut_2d", h);
    }

    // --- 6. M(pi+pi-e+e-) after pippimepem_selection cut chain ------------
    auto* h6 = pu.drawNtupleSingle(
        "pippimepem_nt", "m_pippimepem",
        200, 0.2, 1.4, wcut("sel_pass==1"),
        ";M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}];a.u.");
    auto* cv6 = pu.drawSingle(h6,
                              "M_{#pi^{+}#pi^{-}e^{+}e^{-}} (after selection)",
                              "c_m_pippimepem_selected");
    pu.save(cv6, "m_pippimepem_selected");
    printIntegral("M(pi+pi-e+e-) after selection", h6);

    // --- 6a. Same as plot 6 but with cut_2d additionally applied ----------
    {
        auto* h = pu.drawNtupleSingle(
            "pippimepem_nt", "m_pippimepem",
            200, 0.2, 1.4, wcut("sel_pass==1 && cut2d_pass==1"),
            ";M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}];a.u.");
        auto* cv = pu.drawSingle(h,
                                 "M_{#pi^{+}#pi^{-}e^{+}e^{-}} (after selection + cut_2d)",
                                 "c_m_pippimepem_selected_cut2d");
        pu.save(cv, "m_pippimepem_selected_cut2d");
        printIntegral("M(pi+pi-e+e-) after selection + cut_2d", h);
    }

    // --- 7-12. M(pi+pi-e+e-) in MM(pi+pi-e+e-) slice windows --------------
    auto drawSlice = [&](double lo, double hi, const std::string& tag, bool cut2d = false) {
        std::ostringstream cut;
        cut << "sel_pass==1 && mm_pippimepem>=" << lo
            << " && mm_pippimepem<=" << hi;
        if (cut2d) cut << " && cut2d_pass==1";

        std::ostringstream title;
        title << "M_{#pi^{+}#pi^{-}e^{+}e^{-}} sel" << (cut2d ? "+cut_2d" : "")
              << ", MM #in [" << lo << ", " << hi << "] GeV/c^{2}";

        std::string suffix = tag + (cut2d ? "_cut2d" : "");

        auto* h = pu.drawNtupleSingle(
            "pippimepem_nt", "m_pippimepem",
            100, 0.2, 1.4, wcut(cut.str()),
            ";M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}];a.u.");
        auto* cv = pu.drawSingle(h, title.str(),
                                 "c_m_pippimepem_slice_" + suffix);
        pu.save(cv, "m_pippimepem_slice_" + suffix);

        std::ostringstream lbl;
        lbl << "M(pi+pi-e+e-) sel" << (cut2d ? "+cut_2d" : "")
            << ", MM in [" << lo << ", " << hi << "]";
        printIntegral(lbl.str().c_str(), h);
    };

    drawSlice(2.0, 2.2, "20_22");
    drawSlice(2.2, 2.4, "22_24");
    drawSlice(2.4, 2.6, "24_26");
    drawSlice(2.6, 2.8, "26_28");
    drawSlice(2.8, 3.0, "28_30");

    drawSlice(2.0, 2.2, "20_22", true);
    drawSlice(2.2, 2.4, "22_24", true);
    drawSlice(2.4, 2.6, "24_26", true);
    drawSlice(2.6, 2.8, "26_28", true);
    drawSlice(2.8, 3.0, "28_30", true);

    std::cout << "\nDone. Check plots/output/\n";
}
