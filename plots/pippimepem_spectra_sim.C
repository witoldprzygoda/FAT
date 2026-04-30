// pippimepem_spectra_sim.C — SIMULATED-truth counterpart of pippimepem_spectra(_cor).C
//
// Reads SIMULATED fields from pippimepem_nt_cor (m_*_sim, mm_*_sim, etc.) and
// weights each event by sim_genweight via the TTree::Draw cut expression.
// Selection / cut_2d decisions are taken from CORRECTED kinematics (sel_pass,
// cut2d_pass) — same event population as the _cor macro, but plotted with
// generator-truth values.
//
// Output filenames are suffixed with "_sim" so REC / COR / SIM figures coexist
// in plots/output/.
//
// Usage: root -l -b -q plots/pippimepem_spectra_sim.C

#include "PlotUtils.h"
#include <sstream>

void printIntegral_sim(const char* label, TH1D* h) {
    std::cout << "\n=== " << label << " ===\n";
    std::cout << "  entries = " << h->GetEntries()
              << "   integral = " << h->Integral() << "\n";
}

// Compose a TTree::Draw cut string that filters by `filter` AND weights by
// sim_genweight. ROOT interprets the cut as `(filter)*weight` per entry.
namespace {
    std::string wcut(const std::string& filter = "") {
        return filter.empty() ? std::string("sim_genweight")
                              : "(" + filter + ")*sim_genweight";
    }
}

void pippimepem_spectra_sim() {

    PlotUtils pu("output_pippimepem_sim.root");

    const char* NT = "pippimepem_nt_cor";

    // --- 1. M(e+e-) (sim) ------------------------------------------------
    auto* h1 = pu.drawNtupleSingle(
        NT, "m_ee_sim",
        160, 0.0, 0.8, wcut(),
        ";M_{e^{+}e^{-}} [GeV/c^{2}] (sim);a.u.");
    auto* cv1 = pu.drawSingle(h1, "M_{e^{+}e^{-}} (sim)", "c_m_ee_sim", /*logy=*/true);
    pu.save(cv1, "m_ee_sim");
    printIntegral_sim("M(e+e-) sim", h1);

    // --- 2. M(pi+pi-) (sim) ----------------------------------------------
    auto* h2 = pu.drawNtupleSingle(
        NT, "m_pippim_sim",
        200, 0.0, 2.0, wcut(),
        ";M_{#pi^{+}#pi^{-}} [GeV/c^{2}] (sim);a.u.");
    auto* cv2 = pu.drawSingle(h2, "M_{#pi^{+}#pi^{-}} (sim)", "c_m_pippim_sim");
    pu.save(cv2, "m_pippim_sim");
    printIntegral_sim("M(pi+pi-) sim", h2);

    // --- 3. M(pi+pi-e+e-) (sim) ------------------------------------------
    auto* h3 = pu.drawNtupleSingle(
        NT, "m_pippimepem_sim",
        200, 0.2, 1.4, wcut(),
        ";M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}] (sim);a.u.");
    auto* cv3 = pu.drawSingle(h3,
                              "M_{#pi^{+}#pi^{-}e^{+}e^{-}} (sim)",
                              "c_m_pippimepem_sim");
    pu.save(cv3, "m_pippimepem_sim");
    printIntegral_sim("M(pi+pi-e+e-) sim", h3);

    // --- 4. MM(pi+pi-) (sim) ---------------------------------------------
    auto* h4 = pu.drawNtupleSingle(
        NT, "mm_pippim_sim",
        200, 0.0, 4.0, wcut(),
        ";MM(#pi^{+}#pi^{-}) [GeV/c^{2}] (sim);a.u.");
    auto* cv4 = pu.drawSingle(h4, "MM(#pi^{+}#pi^{-}) (sim)", "c_mm_pippim_sim");
    pu.save(cv4, "mm_pippim_sim");
    printIntegral_sim("MM(pi+pi-) sim", h4);

    // --- 5. MM(pi+pi-e+e-) (sim) -----------------------------------------
    auto* h5 = pu.drawNtupleSingle(
        NT, "mm_pippimepem_sim",
        200, 0.0, 4.0, wcut(),
        ";MM(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}] (sim);a.u.");
    auto* cv5 = pu.drawSingle(h5,
                              "MM(#pi^{+}#pi^{-}e^{+}e^{-}) (sim)",
                              "c_mm_pippimepem_sim");
    pu.save(cv5, "mm_pippimepem_sim");
    printIntegral_sim("MM(pi+pi-e+e-) sim", h5);

    // --- 5a-5c. MM observables after cut_2d (decision from CORRECTED) ----
    {
        auto* h = pu.drawNtupleSingle(
            NT, "mm_epem_sim",
            200, 0.0, 4.0, wcut("cut2d_pass==1"),
            ";MM(e^{+}e^{-}) [GeV/c^{2}] (sim);a.u.");
        auto* cv = pu.drawSingle(h, "MM(e^{+}e^{-}) after cut_2d (sim)",
                                 "c_mm_epem_cut2d_sim");
        pu.save(cv, "mm_epem_cut2d_sim");
        printIntegral_sim("MM(e+e-) after cut_2d sim", h);
    }
    {
        auto* h = pu.drawNtupleSingle(
            NT, "mm_pippim_sim",
            200, 0.0, 4.0, wcut("cut2d_pass==1"),
            ";MM(#pi^{+}#pi^{-}) [GeV/c^{2}] (sim);a.u.");
        auto* cv = pu.drawSingle(h, "MM(#pi^{+}#pi^{-}) after cut_2d (sim)",
                                 "c_mm_pippim_cut2d_sim");
        pu.save(cv, "mm_pippim_cut2d_sim");
        printIntegral_sim("MM(pi+pi-) after cut_2d sim", h);
    }
    {
        auto* h = pu.drawNtupleSingle(
            NT, "mm_pippimepem_sim",
            200, 0.0, 4.0, wcut("cut2d_pass==1"),
            ";MM(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}] (sim);a.u.");
        auto* cv = pu.drawSingle(h, "MM(#pi^{+}#pi^{-}e^{+}e^{-}) after cut_2d (sim)",
                                 "c_mm_pippimepem_cut2d_sim");
        pu.save(cv, "mm_pippimepem_cut2d_sim");
        printIntegral_sim("MM(pi+pi-e+e-) after cut_2d sim", h);
    }

    // --- 6. M(pi+pi-e+e-) after pippimepem_selection (sim) ---------------
    auto* h6 = pu.drawNtupleSingle(
        NT, "m_pippimepem_sim",
        200, 0.2, 1.4, wcut("sel_pass==1"),
        ";M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}] (sim);a.u.");
    auto* cv6 = pu.drawSingle(h6,
                              "M_{#pi^{+}#pi^{-}e^{+}e^{-}} (sim, after selection)",
                              "c_m_pippimepem_selected_sim");
    pu.save(cv6, "m_pippimepem_selected_sim");
    printIntegral_sim("M(pi+pi-e+e-) after selection sim", h6);

    // --- 6a. Same as 6 with cut_2d additionally applied ------------------
    {
        auto* h = pu.drawNtupleSingle(
            NT, "m_pippimepem_sim",
            200, 0.2, 1.4, wcut("sel_pass==1 && cut2d_pass==1"),
            ";M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}] (sim);a.u.");
        auto* cv = pu.drawSingle(h,
                                 "M_{#pi^{+}#pi^{-}e^{+}e^{-}} (sim, after sel + cut_2d)",
                                 "c_m_pippimepem_selected_cut2d_sim");
        pu.save(cv, "m_pippimepem_selected_cut2d_sim");
        printIntegral_sim("M(pi+pi-e+e-) after sel + cut_2d sim", h);
    }

    // --- 7-12. M(pi+pi-e+e-) in MM(pi+pi-e+e-) slice windows -------------
    auto drawSlice = [&](double lo, double hi, const std::string& tag, bool cut2d = false) {
        std::ostringstream filter;
        filter << "sel_pass==1 && mm_pippimepem>=" << lo
               << " && mm_pippimepem<=" << hi;
        if (cut2d) filter << " && cut2d_pass==1";

        std::ostringstream title;
        title << "M_{#pi^{+}#pi^{-}e^{+}e^{-}} (sim) sel" << (cut2d ? "+cut_2d" : "")
              << ", MM #in [" << lo << ", " << hi << "] GeV/c^{2}";

        std::string suffix = tag + (cut2d ? "_cut2d" : "") + "_sim";

        auto* h = pu.drawNtupleSingle(
            NT, "m_pippimepem_sim",
            100, 0.2, 1.4, wcut(filter.str()),
            ";M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}] (sim);a.u.");
        auto* cv = pu.drawSingle(h, title.str(),
                                 "c_m_pippimepem_slice_" + suffix);
        pu.save(cv, "m_pippimepem_slice_" + suffix);

        std::ostringstream lbl;
        lbl << "M(pi+pi-e+e-) sel" << (cut2d ? "+cut_2d" : "")
            << ", MM in [" << lo << ", " << hi << "] sim";
        printIntegral_sim(lbl.str().c_str(), h);
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

    std::cout << "\nDone. Plots in plots/output/ (suffix _sim)\n";
}
