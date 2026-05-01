// mass_spectra_rec_expsim.C — Joint EXP + SIM dilepton invariant mass.
// RECONSTRUCTED kinematics — reads m_ee from pippimepem_nt.
// Each canvas: data 3-curve (all/CB/signal) + sim line rescaled to data
// signal in a control region (above pi0 Dalitz, below resonance pile-up).
//
// Plots: m_ee with no OA cut, then with opening_angle_4 (>4 deg) applied.
// Outputs: joint_mass_ee_*_rec_expsim.{pdf,png}
//
// Usage: root -l -b -q plots/mass_spectra_rec_expsim.C

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

void mass_spectra_rec_expsim() {

    PlotUtils exp("output_pippimepem_exp.root",
                  "output_pippimepep_exp.root",
                  "output_pippimemem_exp.root");
    JointPlotter::SimSource sim("output_pippimepem_sim.root");

    const char* NT = "pippimepem_nt";

    // --- 1. Mass spectrum without OA cut ---
    TH1D *a1, *c1, *s1;
    std::tie(a1, c1, s1) = exp.drawSignal(
        NT, "m_ee", 160, 0, 0.8, "",
        ";M_{e^{+}e^{-}} [GeV/c^{2}];a.u.");

    TH1D* sim1 = sim.draw(NT, "m_ee", 160, 0, 0.8, "");
    JointPlotter::styleSimLine(sim1);
    double scale1 = JointPlotter::rescaleSimInWindow(sim1, s1, kNormLo, kNormHi);

    auto* cv1 = JointPlotter::drawJoint(
        a1, c1, s1, sim1,
        "M_{e^{+}e^{-}} (no OA cut)",
        "c_mass_no_oa_rec_expsim", /*logy=*/true, scale1);
    JointPlotter::save(cv1, "mass_ee_no_oa_rec_expsim");
    printIntegrals("No OA cut", a1, c1, s1, sim1);

    // --- 2. Mass spectrum with opening_angle_4 (oa>4) cut applied ---
    TH1D *a2, *c2, *s2;
    std::tie(a2, c2, s2) = exp.drawSignal(
        NT, "m_ee", 160, 0, 0.8, "oa_pass==1",
        ";M_{e^{+}e^{-}} [GeV/c^{2}];a.u.");

    TH1D* sim2 = sim.draw(NT, "m_ee", 160, 0, 0.8, "oa_pass==1");
    JointPlotter::styleSimLine(sim2);
    double scale2 = JointPlotter::rescaleSimInWindow(sim2, s2, kNormLo, kNormHi);

    auto* cv2 = JointPlotter::drawJoint(
        a2, c2, s2, sim2,
        "M_{e^{+}e^{-}} (OA > 4#circ)",
        "c_mass_oa4_rec_expsim", /*logy=*/true, scale2);
    JointPlotter::save(cv2, "mass_ee_oa4_rec_expsim");
    printIntegrals("OA > 4 deg", a2, c2, s2, sim2);

    std::cout << "\nDone. Plots in plots/output/joint_mass_ee_*_rec_expsim.{pdf,png}\n";
}
