// pippimepem_spectra_cor_sim.C — CORRECTED-kinematics counterpart of pippimepem_spectra_rec_sim.C.
// Simulation mode (single file, no CB).
//
// Reads from `pippimepem_nt_cor` (mirror ntuple where compound observables are
// computed from KinematicType::CORRECTED). Field names match the RECONSTRUCTED
// ntuple, so cut expressions stay identical; only the underlying values differ.
//
// All output filenames are suffixed with "_cor_sim" to keep figures from the two
// macros side by side in plots/output/ without overwrite collisions.
//
// Usage: root -l -b -q plots/pippimepem_spectra_cor_sim.C

#include "PlotUtils.h"
#include <sstream>

// Helper: compose a weighted TTree::Draw cut. Per-event sim_genweight from
// the ntuple is multiplied by an optional boolean filter so all spectra
// reflect the SMASH luminosity normalisation.
namespace { std::string wcut(const std::string& filter = "") {
    return filter.empty() ? std::string("sim_genweight")
                          : "(" + filter + ")*sim_genweight";
}}


void printIntegral_cor_sim(const char* label, TH1D* h) {
    std::cout << "\n=== " << label << " ===\n";
    std::cout << "  entries = " << h->GetEntries()
              << "   integral = " << h->Integral() << "\n";
}

void pippimepem_spectra_cor_sim() {

    PlotUtils pu("output_pippimepem_sim.root");   // single-file (sim, no CB)

    const char* NT = "pippimepem_nt_cor";

    // --- 1. M(e+e-) -------------------------------------------------------
    auto* h1 = pu.drawNtupleSingle(
        NT, "m_ee",
        160, 0.0, 0.8, wcut(""),
        ";M_{e^{+}e^{-}} [GeV/c^{2}];a.u.");
    auto* cv1 = pu.drawSingle(h1, "M_{e^{+}e^{-}} (cor)", "c_m_ee_cor_sim", /*logy=*/true);
    pu.save(cv1, "m_ee_cor_sim");
    printIntegral_cor_sim("M(e+e-) cor", h1);

    // --- 2. M(pi+pi-) -----------------------------------------------------
    auto* h2 = pu.drawNtupleSingle(
        NT, "m_pippim",
        200, 0.0, 2.0, wcut(""),
        ";M_{#pi^{+}#pi^{-}} [GeV/c^{2}];a.u.");
    auto* cv2 = pu.drawSingle(h2, "M_{#pi^{+}#pi^{-}} (cor)", "c_m_pippim_cor_sim");
    pu.save(cv2, "m_pippim_cor_sim");
    printIntegral_cor_sim("M(pi+pi-) cor", h2);

    // --- 3. M(pi+pi-e+e-) -------------------------------------------------
    auto* h3 = pu.drawNtupleSingle(
        NT, "m_pippimepem",
        200, 0.2, 1.4, wcut(""),
        ";M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}];a.u.");
    auto* cv3 = pu.drawSingle(h3,
                              "M_{#pi^{+}#pi^{-}e^{+}e^{-}} (cor)",
                              "c_m_pippimepem_cor_sim");
    pu.save(cv3, "m_pippimepem_cor_sim");
    printIntegral_cor_sim("M(pi+pi-e+e-) cor", h3);

    // --- 4. MM(pi+pi-) ----------------------------------------------------
    auto* h4 = pu.drawNtupleSingle(
        NT, "mm_pippim",
        200, 0.0, 4.0, wcut(""),
        ";MM(#pi^{+}#pi^{-}) [GeV/c^{2}];a.u.");
    auto* cv4 = pu.drawSingle(h4, "MM(#pi^{+}#pi^{-}) (cor)", "c_mm_pippim_cor_sim");
    pu.save(cv4, "mm_pippim_cor_sim");
    printIntegral_cor_sim("MM(pi+pi-) cor", h4);

    // --- 5. MM(pi+pi-e+e-) ------------------------------------------------
    auto* h5 = pu.drawNtupleSingle(
        NT, "mm_pippimepem",
        200, 0.0, 4.0, wcut(""),
        ";MM(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}];a.u.");
    auto* cv5 = pu.drawSingle(h5,
                              "MM(#pi^{+}#pi^{-}e^{+}e^{-}) (cor)",
                              "c_mm_pippimepem_cor_sim");
    pu.save(cv5, "mm_pippimepem_cor_sim");
    printIntegral_cor_sim("MM(pi+pi-e+e-) cor", h5);

    // --- 5a-5c. MM observables after cut_2d -------------------------------
    {
        auto* h = pu.drawNtupleSingle(
            NT, "mm_epem",
            200, 0.0, 4.0, wcut("cut2d_pass==1"),
            ";MM(e^{+}e^{-}) [GeV/c^{2}];a.u.");
        auto* cv = pu.drawSingle(h, "MM(e^{+}e^{-}) after cut_2d (cor)",
                                 "c_mm_epem_cut2d_cor_sim");
        pu.save(cv, "mm_epem_cut2d_cor_sim");
        printIntegral_cor_sim("MM(e+e-) after cut_2d cor", h);
    }
    {
        auto* h = pu.drawNtupleSingle(
            NT, "mm_pippim",
            200, 0.0, 4.0, wcut("cut2d_pass==1"),
            ";MM(#pi^{+}#pi^{-}) [GeV/c^{2}];a.u.");
        auto* cv = pu.drawSingle(h, "MM(#pi^{+}#pi^{-}) after cut_2d (cor)",
                                 "c_mm_pippim_cut2d_cor_sim");
        pu.save(cv, "mm_pippim_cut2d_cor_sim");
        printIntegral_cor_sim("MM(pi+pi-) after cut_2d cor", h);
    }
    {
        auto* h = pu.drawNtupleSingle(
            NT, "mm_pippimepem",
            200, 0.0, 4.0, wcut("cut2d_pass==1"),
            ";MM(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}];a.u.");
        auto* cv = pu.drawSingle(h, "MM(#pi^{+}#pi^{-}e^{+}e^{-}) after cut_2d (cor)",
                                 "c_mm_pippimepem_cut2d_cor_sim");
        pu.save(cv, "mm_pippimepem_cut2d_cor_sim");
        printIntegral_cor_sim("MM(pi+pi-e+e-) after cut_2d cor", h);
    }

    // --- 5d. OA observables driving the pippimepem_selection cut chain ----
    {
        auto* h = pu.drawNtupleSingle(
            NT, "oa_pippim_epem_lab",
            180, 0.0, 180.0, wcut(""),
            ";OA_{LAB}((#pi^{+}#pi^{-}),(e^{+}e^{-})) [deg];a.u.");
        auto* cv = pu.drawSingle(h,
                                 "OA_{LAB}((#pi^{+}#pi^{-}),(e^{+}e^{-})) (cor)",
                                 "c_oa_pippim_epem_lab_cor_sim");
        pu.save(cv, "oa_pippim_epem_lab_cor_sim");
        printIntegral_cor_sim("OA_LAB((pi+pi-),(e+e-)) cor", h);
    }
    {
        auto* h = pu.drawNtupleSingle(
            NT, "oa_pippim_epem_eta_rest",
            180, 0.0, 180.0, wcut(""),
            ";OA_{m_{#eta}-rest}((#pi^{+}#pi^{-}),(e^{+}e^{-})) [deg];a.u.");
        auto* cv = pu.drawSingle(h,
                                 "OA((#pi^{+}#pi^{-}),(e^{+}e^{-})) in m_{#eta}-rest frame (cor)",
                                 "c_oa_pippim_epem_eta_rest_cor_sim");
        pu.save(cv, "oa_pippim_epem_eta_rest_cor_sim");
        printIntegral_cor_sim("OA_eta_rest((pi+pi-),(e+e-)) cor", h);
    }

    // --- 6. M(pi+pi-e+e-) after pippimepem_selection cut chain ------------
    auto* h6 = pu.drawNtupleSingle(
        NT, "m_pippimepem",
        200, 0.2, 1.4, wcut("sel_pass==1"),
        ";M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}];a.u.");
    auto* cv6 = pu.drawSingle(h6,
                              "M_{#pi^{+}#pi^{-}e^{+}e^{-}} (after selection, cor)",
                              "c_m_pippimepem_selected_cor_sim");
    pu.save(cv6, "m_pippimepem_selected_cor_sim");
    printIntegral_cor_sim("M(pi+pi-e+e-) after selection cor", h6);

    // --- 6a. Same as plot 6 but with cut_2d additionally applied ----------
    {
        auto* h = pu.drawNtupleSingle(
            NT, "m_pippimepem",
            200, 0.2, 1.4, wcut("sel_pass==1 && cut2d_pass==1"),
            ";M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}];a.u.");
        auto* cv = pu.drawSingle(h,
                                 "M_{#pi^{+}#pi^{-}e^{+}e^{-}} (after selection + cut_2d, cor)",
                                 "c_m_pippimepem_selected_cut2d_cor_sim");
        pu.save(cv, "m_pippimepem_selected_cut2d_cor_sim");
        printIntegral_cor_sim("M(pi+pi-e+e-) after selection + cut_2d cor", h);
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

        std::string suffix = tag + (cut2d ? "_cut2d" : "") + "_cor_sim";

        auto* h = pu.drawNtupleSingle(
            NT, "m_pippimepem",
            100, 0.2, 1.4, wcut(cut.str()),
            ";M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}];a.u.");
        auto* cv = pu.drawSingle(h, title.str(),
                                 "c_m_pippimepem_slice_" + suffix);
        pu.save(cv, "m_pippimepem_slice_" + suffix);

        std::ostringstream lbl;
        lbl << "M(pi+pi-e+e-) sel" << (cut2d ? "+cut_2d" : "")
            << ", MM in [" << lo << ", " << hi << "] cor";
        printIntegral_cor_sim(lbl.str().c_str(), h);
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

    std::cout << "\nDone. Check plots/output/ for *_cor_sim.{pdf,png}\n";
}
