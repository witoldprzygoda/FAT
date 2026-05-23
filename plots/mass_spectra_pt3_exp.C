// mass_spectra_pt3_exp.C — dilepton mass spectrum (all / CB / signal) built
// from the PT3-selected histograms (`dilepton/mass_ee` in each output file).
//
// "Set 1" of three parallel views of the trigger-handled data:
//     PT3   — biased dilepton trigger (this macro)             histograms H
//     PT2   — downscaled unbiased reference                    histograms H_pt2
//     ratio — 63 · PT2/PT3 trigger-correction shape per bin    derived
//
// Reads pre-built histograms (not ntuples), so no rebinning is done — the
// bin layout is whatever setup_histograms.h fixed for "mass_ee" (200 bins
// over 0..1 GeV/c²). The matching PT2 partner is "mass_ee_pt2", and the
// _ratio macro divides them.
//
// Usage:  root -l -b -q plots/mass_spectra_pt3_exp.C

#include "PlotUtils.h"
#include <string>
#include <tuple>

void mass_spectra_pt3_exp() {
    PlotUtils pu("output_epem_exp.root",
                 "output_epep_exp.root",
                 "output_emem_exp.root");

    // ------------------------------------------------------------------
    // No OA cut: mass_ee is filled for every event passing event-level
    // cuts (isBest, vertex_z, start_detector) — full spectrum.
    // ------------------------------------------------------------------
    TH1D *a1, *c1, *s1;
    std::tie(a1, c1, s1) = pu.getSignal("dilepton/mass_ee");
    if (!a1 || !c1 || !s1) {
        std::cerr << "Missing dilepton/mass_ee in one of the inputs.\n";
        return;
    }

    const char* t_no = "M_{e^{+}e^{-}} (PT3, no OA cut)";
    auto* cv_log = pu.drawTriple(a1, c1, s1, t_no, "c_mass_pt3_log", /*logy=*/true);
    pu.save(cv_log, "mass_ee_pt3_log");

    // Linear variant — drawTriple already applied logY scaling once, so
    // re-cap the maximum from raw bin contents.
    const double a1_lin = a1->GetBinContent(a1->GetMaximumBin());
    auto* cv_lin = pu.drawTriple(a1, c1, s1, t_no, "c_mass_pt3_lin", /*logy=*/false);
    a1->SetMaximum(a1_lin * 1.2);
    a1->SetMinimum(0.0);
    cv_lin->Update();
    pu.save(cv_lin, "mass_ee_pt3_lin");

    // ------------------------------------------------------------------
    // OA > 4 deg active analysis selection (mass_ee_after_oa).
    // ------------------------------------------------------------------
    TH1D *a2, *c2, *s2;
    std::tie(a2, c2, s2) = pu.getSignal("dilepton/mass_ee_after_oa");
    if (a2 && c2 && s2) {
        const char* t_oa = "M_{e^{+}e^{-}} (PT3, OA > 4#circ)";
        auto* cv2_log = pu.drawTriple(a2, c2, s2, t_oa, "c_mass_pt3_oa_log", /*logy=*/true);
        pu.save(cv2_log, "mass_ee_pt3_oa_log");

        const double a2_lin = a2->GetBinContent(a2->GetMaximumBin());
        auto* cv2_lin = pu.drawTriple(a2, c2, s2, t_oa, "c_mass_pt3_oa_lin", /*logy=*/false);
        a2->SetMaximum(a2_lin * 1.2);
        a2->SetMinimum(0.0);
        cv2_lin->Update();
        pu.save(cv2_lin, "mass_ee_pt3_oa_lin");
    }

    // Integrals — useful sanity numbers for cross-check with PT2 / ratio.
    std::cout << "\n=== PT3 integrals (full mass_ee) ===\n"
              << "  all = " << a1->Integral()
              << "  CB = "  << c1->Integral()
              << "  sig = " << s1->Integral() << "\n";
    if (a2 && c2 && s2) {
        std::cout << "=== PT3 integrals (mass_ee_after_oa) ===\n"
                  << "  all = " << a2->Integral()
                  << "  CB = "  << c2->Integral()
                  << "  sig = " << s2->Integral() << "\n";
    }
}
