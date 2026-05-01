// mass_spectra.C — Joint EXP + SIM dilepton invariant mass.
// Each canvas: data 3-curve (all/CB/signal) + sim line rescaled to data
// signal in the right tail (control region, no signal contamination).
//
// Plots: m_ee with no OA cut, then with opening_angle_4 (>4 deg) applied.
// Output: plots/output/joint_mass_ee_*.{pdf,png}
//
// Usage: root -l -b -q plots/mass_spectra.C

#include "PlotUtils.h"
#include "JointPlotter.h"

namespace {
    // Control region for sim->data signal rescaling — right-tail of M(ee)
    // where any meson contribution has died off and integrated counts are
    // dominated by combinatorial-corrected continuum.
    constexpr double kCtrlLo = 0.40;
    constexpr double kCtrlHi = 0.80;

    void printIntegrals(const char* label, TH1D* all, TH1D* cb, TH1D* sig, TH1D* sim) {
        std::cout << "\n=== " << label << " ===\n";
        std::cout << "  exp: all=" << all->Integral()
                  << "   CB=" << cb->Integral()
                  << "   sig=" << sig->Integral() << "\n";
        std::cout << "  sim: " << sim->Integral()
                  << "   (rescaled to exp signal in [" << kCtrlLo << ", " << kCtrlHi << "])\n";
    }
}

void mass_spectra() {

    PlotUtils exp("output_pippimepem_exp.root",
                  "output_pippimepep_exp.root",
                  "output_pippimemem_exp.root");
    JointPlotter::SimSource sim("output_pippimepem_sim.root");

    // --- 1. Mass spectrum without OA cut ---
    TH1D *a1, *c1, *s1;
    std::tie(a1, c1, s1) = exp.drawSignal(
        "pippimepem_nt", "m_ee", 160, 0, 0.8, "",
        ";M_{e^{+}e^{-}} [GeV/c^{2}];a.u.");

    TH1D* sim1 = sim.draw("pippimepem_nt", "m_ee", 160, 0, 0.8, "");
    JointPlotter::styleSimLine(sim1);
    // Norm window per analysis spec: M(ee) ∈ [0.25, 0.40] — above pi0 Dalitz
    // peak, below rho/omega/phi resonance region; smooth continuum dominates.
    double scale1 = JointPlotter::rescaleSimInWindow(sim1, s1, 0.25, 0.40);

    auto* cv1 = JointPlotter::drawJoint(a1, c1, s1, sim1,
                                        "M_{e^{+}e^{-}} (no OA cut)",
                                        "c_mass_no_oa", /*logy=*/true, scale1);
    JointPlotter::save(cv1, "mass_ee_no_oa");
    printIntegrals("No OA cut", a1, c1, s1, sim1);

    // --- 2. Mass spectrum with opening_angle_4 (oa>4) cut applied ---
    TH1D *a2, *c2, *s2;
    std::tie(a2, c2, s2) = exp.drawSignal(
        "pippimepem_nt", "m_ee", 160, 0, 0.8, "oa_pass==1",
        ";M_{e^{+}e^{-}} [GeV/c^{2}];a.u.");

    TH1D* sim2 = sim.draw("pippimepem_nt", "m_ee", 160, 0, 0.8, "oa_pass==1");
    JointPlotter::styleSimLine(sim2);
    double scale2 = JointPlotter::rescaleSimInWindow(sim2, s2, 0.25, 0.40);

    auto* cv2 = JointPlotter::drawJoint(a2, c2, s2, sim2,
                                        "M_{e^{+}e^{-}} (OA > 4#circ)",
                                        "c_mass_oa4", /*logy=*/true, scale2);
    JointPlotter::save(cv2, "mass_ee_oa4");
    printIntegrals("OA > 4 deg", a2, c2, s2, sim2);

    std::cout << "\nDone. Plots in plots/output/joint_*\n";
}
