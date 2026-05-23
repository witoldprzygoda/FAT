// mass_spectra_pt2_exp.C — dilepton mass spectrum (all / CB / signal) built
// from the PT2-selected histograms (`dilepton/mass_ee_pt2` in each file).
//
// "Set 2" of three parallel views (see mass_spectra_pt3_exp.C for the
// PT3 set and mass_spectra_ratio_exp.C for the ratio).
//
// The PT2 trigger is downscaled 64× in HADES, so the per-bin statistics
// here are roughly ÷27 vs. PT3 for the epem channel (and ÷50 for like-sign);
// expect noisier tails. The bin layout matches PT3 exactly because the
// _pt2 twins are clones produced by Manager::createPT2Clones().
//
// Usage:  root -l -b -q plots/mass_spectra_pt2_exp.C

#include "PlotUtils.h"
#include <string>
#include <tuple>

void mass_spectra_pt2_exp() {
    PlotUtils pu("output_epem_exp.root",
                 "output_epep_exp.root",
                 "output_emem_exp.root");

    // ------------------------------------------------------------------
    // No OA cut: mass_ee_pt2 is filled for every PT2 event that passes
    // the same event-level cuts as PT3 (isBest, vertex_z, start_detector).
    // ------------------------------------------------------------------
    TH1D *a1, *c1, *s1;
    std::tie(a1, c1, s1) = pu.getSignal("dilepton/mass_ee_pt2");
    if (!a1 || !c1 || !s1) {
        std::cerr << "Missing dilepton/mass_ee_pt2 in one of the inputs.\n";
        return;
    }

    const char* t_no = "M_{e^{+}e^{-}} (PT2, no OA cut)";
    auto* cv_log = pu.drawTriple(a1, c1, s1, t_no, "c_mass_pt2_log", /*logy=*/true);
    pu.save(cv_log, "mass_ee_pt2_log");

    const double a1_lin = a1->GetBinContent(a1->GetMaximumBin());
    auto* cv_lin = pu.drawTriple(a1, c1, s1, t_no, "c_mass_pt2_lin", /*logy=*/false);
    a1->SetMaximum(a1_lin * 1.2);
    a1->SetMinimum(0.0);
    cv_lin->Update();
    pu.save(cv_lin, "mass_ee_pt2_lin");

    // ------------------------------------------------------------------
    // OA > 4 deg variant.
    // ------------------------------------------------------------------
    TH1D *a2, *c2, *s2;
    std::tie(a2, c2, s2) = pu.getSignal("dilepton/mass_ee_after_oa_pt2");
    if (a2 && c2 && s2) {
        const char* t_oa = "M_{e^{+}e^{-}} (PT2, OA > 4#circ)";
        auto* cv2_log = pu.drawTriple(a2, c2, s2, t_oa, "c_mass_pt2_oa_log", /*logy=*/true);
        pu.save(cv2_log, "mass_ee_pt2_oa_log");

        const double a2_lin = a2->GetBinContent(a2->GetMaximumBin());
        auto* cv2_lin = pu.drawTriple(a2, c2, s2, t_oa, "c_mass_pt2_oa_lin", /*logy=*/false);
        a2->SetMaximum(a2_lin * 1.2);
        a2->SetMinimum(0.0);
        cv2_lin->Update();
        pu.save(cv2_lin, "mass_ee_pt2_oa_lin");
    }

    // Integrals.
    std::cout << "\n=== PT2 integrals (full mass_ee_pt2) ===\n"
              << "  all = " << a1->Integral()
              << "  CB = "  << c1->Integral()
              << "  sig = " << s1->Integral() << "\n";
    if (a2 && c2 && s2) {
        std::cout << "=== PT2 integrals (mass_ee_after_oa_pt2) ===\n"
                  << "  all = " << a2->Integral()
                  << "  CB = "  << c2->Integral()
                  << "  sig = " << s2->Integral() << "\n";
    }
}
