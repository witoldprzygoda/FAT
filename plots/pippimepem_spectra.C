// pippimepem_spectra.C — Joint EXP + SIM for the pi+pi-e+e- observables (RECONSTRUCTED).
// Each canvas: data 3-curve (all/CB/signal) + sim line rescaled to data signal
// in a per-plot control region (right tail / sideband).
//
// Plots: M(e+e-), M(pi+pi-), M(pi+pi-e+e-), MM(pi+pi-), MM(pi+pi-e+e-),
//        the same MM observables after cut_2d, M(pi+pi-e+e-) after selection,
//        and M(pi+pi-e+e-) in MM slice windows (with and without cut_2d).
// Output: plots/output/joint_*.{pdf,png}
//
// Usage: root -l -b -q plots/pippimepem_spectra.C

#include "PlotUtils.h"
#include "JointPlotter.h"
#include <sstream>

namespace {
    void printIntegrals(const char* label, TH1D* all, TH1D* cb, TH1D* sig, TH1D* sim,
                        double lo, double hi) {
        std::cout << "\n=== " << label << " ===\n";
        std::cout << "  exp: all=" << all->Integral()
                  << "   CB=" << cb->Integral()
                  << "   sig=" << sig->Integral() << "\n";
        std::cout << "  sim: " << sim->Integral()
                  << "   (rescaled to exp signal in [" << lo << ", " << hi << "])\n";
    }
}

void pippimepem_spectra() {

    PlotUtils exp("output_pippimepem_exp.root",
                  "output_pippimepep_exp.root",
                  "output_pippimemem_exp.root");
    JointPlotter::SimSource sim("output_pippimepem_sim.root");

    // Helper for a single joint plot. Picks per-plot control region for
    // sim->data signal rescaling; passes through cut + binning + axis labels.
    auto plot = [&](const std::string& var, int nbins, double xmin, double xmax,
                    const std::string& cut_filter, const std::string& title,
                    const std::string& basename, double ctrl_lo, double ctrl_hi,
                    bool logy = false) {
        std::string axis = ";" + title + " [GeV/c^{2}];a.u.";

        TH1D *a, *c, *s;
        std::tie(a, c, s) = exp.drawSignal(
            "pippimepem_nt", var, nbins, xmin, xmax, cut_filter, axis);

        TH1D* h_sim = sim.draw("pippimepem_nt", var, nbins, xmin, xmax, cut_filter);
        JointPlotter::styleSimLine(h_sim);
        double scale = JointPlotter::rescaleSimToData(h_sim, s);

        auto* cv = JointPlotter::drawJoint(a, c, s, h_sim, title,
                                           "c_" + basename, logy, scale);
        JointPlotter::save(cv, basename);
        printIntegrals(basename.c_str(), a, c, s, h_sim, ctrl_lo, ctrl_hi);
    };

    // --- M(e+e-) deliberately NOT here — duplicate of mass_spectra.C plots
    //     (joint_mass_ee_no_oa / joint_mass_ee_oa4). This macro focuses on
    //     pippim/compound observables to avoid generating identical figures.

    // --- 1. M(pi+pi-)      ctrl = [0.7, 2.0]  (above rho)              ------
    plot("m_pippim", 200, 0.0, 2.0, "",
         "M_{#pi^{+}#pi^{-}}", "m_pippim", 0.70, 2.00);

    // --- 3. M(pi+pi-e+e-)  ctrl = [1.0, 1.4]  (above eta region)       ------
    plot("m_pippimepem", 200, 0.2, 1.4, "",
         "M_{#pi^{+}#pi^{-}e^{+}e^{-}}", "m_pippimepem", 1.00, 1.40);

    // --- 4. MM(pi+pi-)     ctrl = [3.0, 4.0]  (above mass thresholds)  ------
    plot("mm_pippim", 200, 0.0, 4.0, "",
         "MM(#pi^{+}#pi^{-})", "mm_pippim", 3.00, 4.00);

    // --- 5. MM(pi+pi-e+e-) ctrl = [3.0, 4.0]                            ------
    plot("mm_pippimepem", 200, 0.0, 4.0, "",
         "MM(#pi^{+}#pi^{-}e^{+}e^{-})", "mm_pippimepem", 3.00, 4.00);

    // --- 5a-5c. MM observables after cut_2d ----------------------------------
    plot("mm_epem",       200, 0.0, 4.0, "cut2d_pass==1",
         "MM(e^{+}e^{-}) after cut_2d", "mm_epem_cut2d", 3.00, 4.00);
    plot("mm_pippim",     200, 0.0, 4.0, "cut2d_pass==1",
         "MM(#pi^{+}#pi^{-}) after cut_2d", "mm_pippim_cut2d", 3.00, 4.00);
    plot("mm_pippimepem", 200, 0.0, 4.0, "cut2d_pass==1",
         "MM(#pi^{+}#pi^{-}e^{+}e^{-}) after cut_2d", "mm_pippimepem_cut2d", 3.00, 4.00);

    // --- 6. M(pi+pi-e+e-) after pippimepem_selection -------------------------
    plot("m_pippimepem", 200, 0.2, 1.4, "sel_pass==1",
         "M_{#pi^{+}#pi^{-}e^{+}e^{-}} (after selection)",
         "m_pippimepem_selected", 1.00, 1.40);

    // --- 6a. Same with cut_2d additionally applied ---------------------------
    plot("m_pippimepem", 200, 0.2, 1.4, "sel_pass==1 && cut2d_pass==1",
         "M_{#pi^{+}#pi^{-}e^{+}e^{-}} (after selection + cut_2d)",
         "m_pippimepem_selected_cut2d", 1.00, 1.40);

    // --- 7-12. M(pi+pi-e+e-) in MM slice windows -----------------------------
    auto plotSlice = [&](double lo, double hi, const std::string& tag, bool cut2d = false) {
        std::ostringstream cut;
        cut << "sel_pass==1 && mm_pippimepem>=" << lo
            << " && mm_pippimepem<=" << hi;
        if (cut2d) cut << " && cut2d_pass==1";

        std::ostringstream title;
        title << "M_{#pi^{+}#pi^{-}e^{+}e^{-}} sel" << (cut2d ? "+cut_2d" : "")
              << ", MM #in [" << lo << ", " << hi << "] GeV/c^{2}";

        std::string suffix = tag + (cut2d ? "_cut2d" : "");
        plot("m_pippimepem", 100, 0.2, 1.4, cut.str(), title.str(),
             "m_pippimepem_slice_" + suffix, 1.00, 1.40);
    };

    plotSlice(2.0, 2.2, "20_22");
    plotSlice(2.2, 2.4, "22_24");
    plotSlice(2.4, 2.6, "24_26");
    plotSlice(2.6, 2.8, "26_28");
    plotSlice(2.8, 3.0, "28_30");

    plotSlice(2.0, 2.2, "20_22", true);
    plotSlice(2.2, 2.4, "22_24", true);
    plotSlice(2.4, 2.6, "24_26", true);
    plotSlice(2.6, 2.8, "26_28", true);
    plotSlice(2.8, 3.0, "28_30", true);

    std::cout << "\nDone. Plots in plots/output/joint_*\n";
}
