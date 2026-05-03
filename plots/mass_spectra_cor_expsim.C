// mass_spectra_cor_expsim.C — Joint EXP + SIM dilepton invariant mass.
// CORRECTED kinematics — reads m_ee from dilepton_nt_cor on both sides.
//
// Plots: m_ee with no OA cut, then with the active OA>4 cut (oa_pass==1 flag
// from passMinCut("opening_angle_4") in main.cc).
// Output: plots/output/joint_mass_ee_*_cor_expsim.{pdf,png}
//
// Usage: root -l -b -q plots/mass_spectra_cor_expsim.C

#include "PlotUtils.h"
#include "JointPlotter.h"
#include <string>

namespace {
    constexpr double kNormLo = 0.25;
    constexpr double kNormHi = 0.40;

    void printIntegrals(const char* label, TH1D* all, TH1D* cb, TH1D* sig, TH1D* sim) {
        std::cout << "\n=== " << label << " ===\n";
        std::cout << "  exp: all=" << all->Integral()
                  << "   CB=" << cb->Integral()
                  << "   sig=" << sig->Integral() << "\n";
        if (sim) {
            std::cout << "  sim: " << sim->Integral()
                      << "   (rescaled to exp signal in [" << kNormLo << ", " << kNormHi << "])\n";
        } else {
            std::cout << "  sim: (no contribution)\n";
        }
    }
}

void mass_spectra_cor_expsim() {

    PlotUtils exp("output_epem_exp.root",
                  "output_epep_exp.root",
                  "output_emem_exp.root");
    JointPlotter::SimSource sim("output_epem_sim.root");

    const char* NT = "dilepton_nt_cor";

    // --- 1. Mass spectrum (cor) without OA cut ---
    TH1D *a1, *c1, *s1;
    std::tie(a1, c1, s1) = exp.drawSignal(
        NT, "m_ee", 280, 0, 1.4, "",
        ";M_{e^{+}e^{-}} [GeV/c^{2}];Counts");

    TH1D* sim1 = sim.draw(NT, "m_ee", 280, 0, 1.4, "");
    JointPlotter::styleSimLine(sim1);
    double scale1 = JointPlotter::rescaleSimInWindow(sim1, s1, kNormLo, kNormHi);

    auto* cv1 = JointPlotter::drawJoint(
        a1, c1, s1, sim1,
        "M_{e^{+}e^{-}} (cor, no OA cut)",
        "c_mass_no_oa_cor_expsim", /*logy=*/true, scale1);
    JointPlotter::save(cv1, "mass_ee_no_oa_cor_expsim");
    printIntegrals("No OA cut (cor)", a1, c1, s1, sim1);

    // --- 2. Mass spectrum (cor) with active OA cut (OA > 4 deg, oa_pass flag) ---
    TH1D *a2, *c2, *s2;
    std::tie(a2, c2, s2) = exp.drawSignal(
        NT, "m_ee", 280, 0, 1.4, "oa_pass==1",
        ";M_{e^{+}e^{-}} [GeV/c^{2}];Counts");

    TH1D* sim2 = sim.draw(NT, "m_ee", 280, 0, 1.4, "oa_pass==1");
    JointPlotter::styleSimLine(sim2);
    double scale2 = JointPlotter::rescaleSimInWindow(sim2, s2, kNormLo, kNormHi);

    auto* cv2 = JointPlotter::drawJoint(
        a2, c2, s2, sim2,
        "M_{e^{+}e^{-}} (cor, OA > 4#circ, active)",
        "c_mass_oa4_cor_expsim", /*logy=*/true, scale2);
    JointPlotter::save(cv2, "mass_ee_oa4_cor_expsim");
    printIntegrals("OA > 4 deg (cor, active)", a2, c2, s2, sim2);

    std::cout << "\nDone. Plots in plots/output/joint_mass_ee_*_cor_expsim.{pdf,png}\n";
}
