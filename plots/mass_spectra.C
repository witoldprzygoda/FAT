// mass_spectra.C — Joint EXP + SIM dilepton invariant mass (REC and COR).
// Each canvas: data 3-curve (all/CB/signal) + sim line rescaled to data
// signal in a control region (above pi0 Dalitz, below resonance pile-up).
//
// Plots: m_ee with no OA cut, then with opening_angle_4 (>4 deg) applied.
// Two ntuple sources:
//   pippimepem_nt      — RECONSTRUCTED kinematics → joint_mass_ee_*.{pdf,png}
//   pippimepem_nt_cor  — CORRECTED kinematics    → joint_mass_ee_*_cor.{pdf,png}
// Field names match across the two ntuples; only the underlying values differ.
//
// Usage: root -l -b -q plots/mass_spectra.C

#include "PlotUtils.h"
#include "JointPlotter.h"
#include <string>

namespace {
    // Norm window for sim->data signal rescaling: M(ee) ∈ [0.25, 0.40] —
    // above pi0 Dalitz peak, below rho/omega/phi resonance region; smooth
    // continuum dominates and meson contamination is small.
    constexpr double kNormLo = 0.25;
    constexpr double kNormHi = 0.40;

    void printIntegrals(const char* label, TH1D* all, TH1D* cb, TH1D* sig, TH1D* sim) {
        std::cout << "\n=== " << label << " ===\n";
        std::cout << "  exp: all=" << all->Integral()
                  << "   CB=" << cb->Integral()
                  << "   sig=" << sig->Integral() << "\n";
        std::cout << "  sim: " << sim->Integral()
                  << "   (rescaled to exp signal in [" << kNormLo << ", " << kNormHi << "])\n";
    }
}

void mass_spectra() {

    PlotUtils exp("output_pippimepem_exp.root",
                  "output_pippimepep_exp.root",
                  "output_pippimemem_exp.root");
    JointPlotter::SimSource sim("output_pippimepem_sim.root");

    // Helper: one joint plot for a chosen (ntuple, output suffix, kinematics tag).
    // suffix=""  → REC outputs (joint_mass_ee_*),
    // suffix="_cor" → COR outputs (joint_mass_ee_*_cor).
    auto makeMee = [&](const std::string& nt, const std::string& suffix,
                       const std::string& title_tag) {
        // --- 1. Mass spectrum without OA cut ---
        TH1D *a1, *c1, *s1;
        std::tie(a1, c1, s1) = exp.drawSignal(
            nt, "m_ee", 160, 0, 0.8, "",
            ";M_{e^{+}e^{-}} [GeV/c^{2}];a.u.");

        TH1D* sim1 = sim.draw(nt, "m_ee", 160, 0, 0.8, "");
        JointPlotter::styleSimLine(sim1);
        double scale1 = JointPlotter::rescaleSimInWindow(sim1, s1, kNormLo, kNormHi);

        auto* cv1 = JointPlotter::drawJoint(
            a1, c1, s1, sim1,
            "M_{e^{+}e^{-}} (no OA cut" + title_tag + ")",
            "c_mass_no_oa" + suffix, /*logy=*/true, scale1);
        JointPlotter::save(cv1, "mass_ee_no_oa" + suffix);
        printIntegrals(("No OA cut" + title_tag).c_str(), a1, c1, s1, sim1);

        // --- 2. Mass spectrum with opening_angle_4 (oa>4) cut applied ---
        TH1D *a2, *c2, *s2;
        std::tie(a2, c2, s2) = exp.drawSignal(
            nt, "m_ee", 160, 0, 0.8, "oa_pass==1",
            ";M_{e^{+}e^{-}} [GeV/c^{2}];a.u.");

        TH1D* sim2 = sim.draw(nt, "m_ee", 160, 0, 0.8, "oa_pass==1");
        JointPlotter::styleSimLine(sim2);
        double scale2 = JointPlotter::rescaleSimInWindow(sim2, s2, kNormLo, kNormHi);

        auto* cv2 = JointPlotter::drawJoint(
            a2, c2, s2, sim2,
            "M_{e^{+}e^{-}} (OA > 4#circ" + title_tag + ")",
            "c_mass_oa4" + suffix, /*logy=*/true, scale2);
        JointPlotter::save(cv2, "mass_ee_oa4" + suffix);
        printIntegrals(("OA > 4 deg" + title_tag).c_str(), a2, c2, s2, sim2);
    };

    makeMee("pippimepem_nt",     "",     "");
    makeMee("pippimepem_nt_cor", "_cor", ", cor");

    std::cout << "\nDone. Plots in plots/output/joint_*\n";
}
