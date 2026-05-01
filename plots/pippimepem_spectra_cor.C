// pippimepem_spectra_cor.C — Joint EXP + SIM (CORRECTED kinematics).
// Mirror of pippimepem_spectra.C but reading the *_nt_cor ntuple and
// suffixing output filenames with "_cor".
//
// Per-plot config: norm_lo/hi for sim->exp-signal integral matching, and
// optional display_lo/hi to zoom the X axis (binning unchanged).
// norm_lo == norm_hi means "fall back to findBestScale (auto)".
//
// Usage: root -l -b -q plots/pippimepem_spectra_cor.C

#include "PlotUtils.h"
#include "JointPlotter.h"
#include <sstream>

namespace {
    void printIntegrals_cor(const char* label, TH1D* all, TH1D* cb, TH1D* sig, TH1D* sim) {
        std::cout << "\n=== " << label << " ===\n";
        std::cout << "  exp: all=" << all->Integral()
                  << "   CB="  << cb ->Integral()
                  << "   sig=" << sig->Integral() << "\n";
        std::cout << "  sim: " << sim->Integral() << "\n";
    }
}

void pippimepem_spectra_cor() {

    PlotUtils exp("output_pippimepem_exp.root",
                  "output_pippimepep_exp.root",
                  "output_pippimemem_exp.root");
    JointPlotter::SimSource sim("output_pippimepem_sim.root");

    const char* NT = "pippimepem_nt_cor";

    // Per-plot helper. Norm window from analysis spec; display zoom optional.
    auto plot = [&](const std::string& var, int nbins, double xmin, double xmax,
                    const std::string& cut_filter, const std::string& title,
                    const std::string& basename,
                    double norm_lo, double norm_hi,
                    double display_lo = 0.0, double display_hi = 0.0,
                    bool logy = false) {
        std::string axis = ";" + title + " [GeV/c^{2}];a.u.";

        TH1D *a, *c, *s;
        std::tie(a, c, s) = exp.drawSignal(NT, var, nbins, xmin, xmax, cut_filter, axis);

        TH1D* h_sim = sim.draw(NT, var, nbins, xmin, xmax, cut_filter);
        JointPlotter::styleSimLine(h_sim);
        double scale = (norm_lo < norm_hi)
            ? JointPlotter::rescaleSimInWindow(h_sim, s, norm_lo, norm_hi)
            : JointPlotter::rescaleSimToData(h_sim, s);

        auto* cv = JointPlotter::drawJoint(a, c, s, h_sim, title,
                                           "c_" + basename, logy, scale,
                                           display_lo, display_hi);
        JointPlotter::save(cv, basename);
        printIntegrals_cor(basename.c_str(), a, c, s, h_sim);
    };

    // ---- M(e+e-)  norm [0.25, 0.40]  (above pi0 Dalitz, below resonances) ----
    plot("m_ee",          160, 0.0, 0.8, "",
         "M_{e^{+}e^{-}} (cor)", "m_ee_cor",
         /*norm*/ 0.25, 0.40, /*disp*/ 0, 0, /*logy*/ true);

    // ---- M(pi+pi-)  display [0, 1.4],  norm [0.6, 1.2] ----
    plot("m_pippim",      200, 0.0, 2.0, "",
         "M_{#pi^{+}#pi^{-}} (cor)", "m_pippim_cor",
         /*norm*/ 0.60, 1.20, /*disp*/ 0.0, 1.4);

    // ---- M(pi+pi-e+e-)  display [0.2, 1.9],  norm [1.0, 1.2] ----
    plot("m_pippimepem",  200, 0.2, 1.9, "",
         "M_{#pi^{+}#pi^{-}e^{+}e^{-}} (cor)", "m_pippimepem_cor",
         /*norm*/ 1.00, 1.20, /*disp*/ 0.2, 1.9);

    // ---- MM(pi+pi-)  display [1.5, 4],  norm [2.8, 3.2] ----
    plot("mm_pippim",     200, 0.0, 4.0, "",
         "MM(#pi^{+}#pi^{-}) (cor)", "mm_pippim_cor",
         /*norm*/ 2.80, 3.20, /*disp*/ 1.5, 4.0);

    // ---- MM(pi+pi-e+e-)  display [1, 3.5],  norm [2.6, 3.0] ----
    plot("mm_pippimepem", 200, 0.0, 4.0, "",
         "MM(#pi^{+}#pi^{-}e^{+}e^{-}) (cor)", "mm_pippimepem_cor",
         /*norm*/ 2.60, 3.00, /*disp*/ 1.0, 3.5);

    // ---- After cut_2d ----
    // mm_epem_cut2d:    display [1.5, 4],  norm [3.3, 3.5]
    plot("mm_epem",       200, 0.0, 4.0, "cut2d_pass==1",
         "MM(e^{+}e^{-}) after cut_2d (cor)", "mm_epem_cut2d_cor",
         /*norm*/ 3.30, 3.50, /*disp*/ 1.5, 4.0);
    // mm_pippim_cut2d:  display [1.5, 4],  norm [2.8, 3.2]
    plot("mm_pippim",     200, 0.0, 4.0, "cut2d_pass==1",
         "MM(#pi^{+}#pi^{-}) after cut_2d (cor)", "mm_pippim_cut2d_cor",
         /*norm*/ 2.80, 3.20, /*disp*/ 1.5, 4.0);
    // mm_pippimepem_cut2d:  display [1, 3.5],  norm [2.6, 3.0]
    plot("mm_pippimepem", 200, 0.0, 4.0, "cut2d_pass==1",
         "MM(#pi^{+}#pi^{-}e^{+}e^{-}) after cut_2d (cor)", "mm_pippimepem_cut2d_cor",
         /*norm*/ 2.60, 3.00, /*disp*/ 1.0, 3.5);

    // ---- M(pippimepem) after selection:  norm [0.8, 1.0]  (no display zoom) ----
    plot("m_pippimepem", 200, 0.2, 1.4, "sel_pass==1",
         "M_{#pi^{+}#pi^{-}e^{+}e^{-}} (sel, cor)",
         "m_pippimepem_selected_cor",
         /*norm*/ 0.80, 1.00);

    // ---- selected + cut_2d (no analysis spec yet) — fallback findBestScale ----
    plot("m_pippimepem", 200, 0.2, 1.4, "sel_pass==1 && cut2d_pass==1",
         "M_{#pi^{+}#pi^{-}e^{+}e^{-}} (sel + cut_2d, cor)",
         "m_pippimepem_selected_cut2d_cor",
         /*norm*/ 0, 0);

    // ---- Slices (no per-plot spec yet — fallback findBestScale) ----
    auto plotSlice = [&](double lo, double hi, const std::string& tag, bool cut2d = false) {
        std::ostringstream cut;
        cut << "sel_pass==1 && mm_pippimepem>=" << lo
            << " && mm_pippimepem<=" << hi;
        if (cut2d) cut << " && cut2d_pass==1";

        std::ostringstream title;
        title << "M_{#pi^{+}#pi^{-}e^{+}e^{-}} sel" << (cut2d ? "+cut_2d" : "")
              << ", MM #in [" << lo << ", " << hi << "] GeV/c^{2} (cor)";

        std::string suffix = tag + (cut2d ? "_cut2d" : "") + "_cor";
        plot("m_pippimepem", 100, 0.2, 1.4, cut.str(), title.str(),
             "m_pippimepem_slice_" + suffix, /*norm*/ 0, 0);
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
