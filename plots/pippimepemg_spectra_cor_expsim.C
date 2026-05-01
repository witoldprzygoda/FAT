// pippimepemg_spectra_cor_expsim.C — Joint EXP + SIM for the pi+pi-pi0 candidate
// (CORRECTED kinematics).
//
// pippimepemg = pippim + (epem + gamma_ECAL) under:
//   - ecal_mult == 1
//   - ecal_quality CutSet
//   - M(e+e-gamma) in [0.125, 0.145] GeV/c^2 (pi0 Dalitz, narrow — ACTIVE).
//     Selection driven by pippimepemg_pass_narrow==1; the wide [0.10, 0.18]
//     cut (pippimepemg_pass) is also stored for fallback / systematics.
//
// For each cut configuration (selection, cut_2d, MM slices) we produce two
// joint exp+sim plots:
//   M(pippimepemg)   — output basename "pippimepemg_<tag>_cor_expsim"
//   MM(pippimepemg)  — output basename "pippimepemg_mm_<tag>_cor_expsim"
//                       (= beam+target − pi+pi-e+e-γ, range [0, 4] GeV/c²)
// The earlier "overlay rescaled M(pippimepem) full vs pi+pi-pi0-tagged" plots
// are intentionally dropped — they were comparing two subsets of the same
// dataset, redundant in the joint exp+sim view.
//
// Output: plots/output/joint_pippimepemg_*_cor_expsim.{pdf,png}
//
// Usage: root -l -b -q plots/pippimepemg_spectra_cor_expsim.C

#include "PlotUtils.h"
#include "JointPlotter.h"
#include <sstream>

namespace {
    void printIntegrals_cor_expsim(const char* label, TH1D* all, TH1D* cb, TH1D* sig, TH1D* sim,
                            double lo, double hi) {
        std::cout << "\n=== " << label << " ===\n";
        std::cout << "  exp: all=" << all->Integral()
                  << "   CB=" << cb->Integral()
                  << "   sig=" << sig->Integral() << "\n";
        std::cout << "  sim: " << sim->Integral()
                  << "   (rescaled to exp signal in [" << lo << ", " << hi << "])\n";
    }
}

void pippimepemg_spectra_cor_expsim() {

    PlotUtils exp("output_pippimepem_exp.root",
                  "output_pippimepep_exp.root",
                  "output_pippimemem_exp.root");
    JointPlotter::SimSource sim("output_pippimepem_sim.root");

    const char* NT = "pippimepem_nt_cor";

    // Always require pippimepemg_pass_narrow==1 (mult==1 ECAL gamma, narrow pi0 window).
    auto plot = [&](const std::string& cut_base, const std::string& title,
                    const std::string& basename, int nbins,
                    double xmin = 0.2, double xmax = 1.4,
                    double ctrl_lo = 1.00, double ctrl_hi = 1.40) {
        std::string cut_filter = cut_base.empty()
                                 ? std::string("pippimepemg_pass_narrow==1")
                                 : cut_base + " && pippimepemg_pass_narrow==1";
        std::string axis = ";M_{#pi^{+}#pi^{-}#pi^{0}} [GeV/c^{2}];a.u.";

        TH1D *a, *c, *s;
        std::tie(a, c, s) = exp.drawSignal(
            NT, "m_pippimepemg", nbins, xmin, xmax, cut_filter, axis);

        TH1D* h_sim = sim.draw(NT, "m_pippimepemg", nbins, xmin, xmax, cut_filter);
        JointPlotter::styleSimLine(h_sim);
        double scale = JointPlotter::rescaleSimToData(h_sim, s);

        auto* cv = JointPlotter::drawJoint(a, c, s, h_sim, title,
                                           "c_pippimepemg_" + basename + "_cor_expsim", false, scale);
        JointPlotter::save(cv, "pippimepemg_" + basename + "_cor_expsim");
        printIntegrals_cor_expsim(("pippimepemg " + basename + " cor_expsim").c_str(), a, c, s, h_sim, ctrl_lo, ctrl_hi);
    };

    // MM(pippimepemg) — beam+target − pi+pi-e+e-γ, full kinematic range [0, 4] GeV/c²
    auto plotMM = [&](const std::string& cut_base, const std::string& title,
                      const std::string& basename, int nbins) {
        std::string cut_filter = cut_base.empty()
                                 ? std::string("pippimepemg_pass_narrow==1")
                                 : cut_base + " && pippimepemg_pass_narrow==1";
        std::string axis = ";MM(#pi^{+}#pi^{-}e^{+}e^{-}#gamma) [GeV/c^{2}];a.u.";

        TH1D *a, *c, *s;
        std::tie(a, c, s) = exp.drawSignal(
            NT, "mm_pippimepemg", nbins, 0.0, 4.0, cut_filter, axis);

        TH1D* h_sim = sim.draw(NT, "mm_pippimepemg", nbins, 0.0, 4.0, cut_filter);
        JointPlotter::styleSimLine(h_sim);
        double scale = JointPlotter::rescaleSimToData(h_sim, s);

        auto* cv = JointPlotter::drawJoint(a, c, s, h_sim, title,
                                           "c_pippimepemg_mm_" + basename + "_cor_expsim",
                                           false, scale);
        JointPlotter::save(cv, "pippimepemg_mm_" + basename + "_cor_expsim");
        printIntegrals_cor_expsim(("pippimepemg MM " + basename + " cor_expsim").c_str(),
                                   a, c, s, h_sim, 0.0, 4.0);
    };

    // Cut configurations (mirror those in pippimepem_spectra_cor.C).
    // No-cut plot uses a wider X range [0.4, 1.8] to capture the f1(1285) tail;
    // post-selection plots stay on [0.2, 1.4] (default) where statistics live.
    plot("",                                "M (no cut, #pi^{+}#pi^{-}#pi^{0})",
         "base", 200, /*xmin*/ 0.4, /*xmax*/ 1.8);
    plot("sel_pass==1",                     "M (after selection, #pi^{+}#pi^{-}#pi^{0})",
         "selected", 200);
    plot("sel_pass==1 && cut2d_pass==1",    "M (after selection + cut_2d, #pi^{+}#pi^{-}#pi^{0})",
         "selected_cut2d", 200);

    plotMM("",                              "MM (no cut, #pi^{+}#pi^{-}#pi^{0})",
           "base", 200);
    plotMM("sel_pass==1",                   "MM (after selection, #pi^{+}#pi^{-}#pi^{0})",
           "selected", 200);
    plotMM("sel_pass==1 && cut2d_pass==1",  "MM (after selection + cut_2d, #pi^{+}#pi^{-}#pi^{0})",
           "selected_cut2d", 200);

    // Slice configurations — same MM(pippimepem) windows as the main macro
    auto plotSlice = [&](double lo, double hi, const std::string& tag, bool cut2d = false) {
        std::ostringstream cut;
        cut << "sel_pass==1 && mm_pippimepem>=" << lo
            << " && mm_pippimepem<=" << hi;
        if (cut2d) cut << " && cut2d_pass==1";
        std::string cut_str = cut.str();

        std::ostringstream title_m;
        title_m << "M (#pi^{+}#pi^{-}#pi^{0}) sel" << (cut2d ? "+cut_2d" : "")
                << ", MM #in [" << lo << ", " << hi << "] GeV/c^{2}";
        std::ostringstream title_mm;
        title_mm << "MM (#pi^{+}#pi^{-}#pi^{0}) sel" << (cut2d ? "+cut_2d" : "")
                 << ", MM #in [" << lo << ", " << hi << "] GeV/c^{2}";

        std::string suffix = "slice_" + tag + (cut2d ? "_cut2d" : "");
        plot  (cut_str, title_m.str(),  suffix, 100);
        plotMM(cut_str, title_mm.str(), suffix, 100);
    };

    plotSlice(2.0, 2.2, "20_22", false);
    plotSlice(2.2, 2.4, "22_24", false);
    plotSlice(2.4, 2.6, "24_26", false);
    plotSlice(2.6, 2.8, "26_28", false);
    plotSlice(2.8, 3.0, "28_30", false);

    plotSlice(2.0, 2.2, "20_22", true);
    plotSlice(2.2, 2.4, "22_24", true);
    plotSlice(2.4, 2.6, "24_26", true);
    plotSlice(2.6, 2.8, "26_28", true);
    plotSlice(2.8, 3.0, "28_30", true);

    std::cout << "\nDone. Plots in plots/output/joint_pippimepemg_*\n";
}
