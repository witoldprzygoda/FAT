// pippimepemg_spectra_cor.C — eta-candidate quality study (CORRECTED only).
//
// pippimepemg = pippim + (epem + gamma_ECAL) under:
//   - ecal_mult == 1
//   - ecal_quality CutSet
//   - M(e+e-gamma) in [0.10, 0.18] GeV/c^2 (pi0 Dalitz window)
// Fields read from pippimepem_nt_cor: m_pippimepem, m_pippimepemg, pippimepemg_pass,
// sel_pass, cut2d_pass, mm_pippimepem.
//
// Two sections:
//   (A) OVERLAY — for each cut configuration, draw M(pippimepem) all/CB/signal and
//       overlay the signal of M(pippimepemg) rescaled to match the pippimepem signal
//       integral in the M ~ [0.45, 0.55] GeV/c^2 window.
//   (B) STANDALONE — M(pippimepemg) all/CB/signal under the same cut configurations
//       (no overlay), produced exactly like pippimepem_spectra_cor.C does for M(pippimepem).
//
// All output PDFs/PNGs land in plots/output/ with prefix "pippimepemg_*".
//
// Usage: root -l -b -q plots/pippimepemg_spectra_cor.C

#include "PlotUtils.h"
#include <sstream>

namespace {
    // probe window for the rescaling (centered at M ~ 0.5 GeV)
    constexpr double kProbeLo = 0.45;
    constexpr double kProbeHi = 0.55;
}

void pippimepemg_print(const char* label, TH1D* a, TH1D* c, TH1D* s) {
    std::cout << "\n=== " << label << " ===\n";
    std::cout << "  all = " << a->Integral()
              << "   CB = " << c->Integral()
              << "   sig = " << s->Integral() << "\n";
}

void pippimepemg_spectra_cor() {

    PlotUtils pu("output_pippimepem.root",
                 "output_pippimepep.root",
                 "output_pippimemem.root");

    const char* NT = "pippimepem_nt_cor";

    // -----------------------------------------------------------------
    // SECTION A: OVERLAY — same observable M(pippimepem), two samples
    // -----------------------------------------------------------------
    // First triple: M(pippimepem) all/CB/signal on the FULL sample (cut_base only).
    // Second curve: same M(pippimepem) signal but ONLY on events tagged as having
    // a good pippimepemg candidate (pippimepemg_pass==1). Both spectra plot the
    // identical observable — they differ only in which events feed them, so the
    // background shape should match and the rescale at M~0.5 GeV/c^2 is meaningful.
    auto drawOverlay = [&](const std::string& cut_base, const std::string& title,
                           const std::string& tag, int nbins) {
        std::string cut_tagged = cut_base.empty()
                                 ? std::string("pippimepemg_pass==1")
                                 : cut_base + " && pippimepemg_pass==1";

        TH1D *a1, *c1, *s1, *a2, *c2, *s2;
        std::tie(a1, c1, s1) = pu.drawSignal(
            NT, "m_pippimepem", nbins, 0.2, 1.4, cut_base,
            ";M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}];Counts");
        // SAME variable, SAME binning — just the "tagged" subset
        std::tie(a2, c2, s2) = pu.drawSignal(
            NT, "m_pippimepem", nbins, 0.2, 1.4, cut_tagged,
            ";M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}];Counts");

        // Rescale the tagged signal to match full-sample signal in the probe window
        int blo = s1->FindBin(kProbeLo);
        int bhi = s1->FindBin(kProbeHi);
        double y_full   = s1->Integral(blo, bhi);
        double y_tagged = s2->Integral(blo, bhi);
        double scale    = (y_tagged > 0) ? y_full / y_tagged : 1.0;
        s2->Scale(scale);

        // Custom canvas with overlay
        TCanvas* cv = new TCanvas(("c_overlay_" + tag).c_str(), title.c_str(), 800, 600);
        cv->SetMargin(0.12, 0.05, 0.12, 0.08);

        a1->SetMarkerStyle(20); a1->SetMarkerSize(0.7);
        a1->SetMarkerColor(kBlack); a1->SetLineColor(kBlack);
        c1->SetMarkerStyle(20); c1->SetMarkerSize(0.7);
        c1->SetMarkerColor(kRed);   c1->SetLineColor(kRed);
        s1->SetMarkerStyle(20); s1->SetMarkerSize(0.7);
        s1->SetMarkerColor(kBlue);  s1->SetLineColor(kBlue);
        s2->SetMarkerStyle(24); s2->SetMarkerSize(0.8);
        s2->SetMarkerColor(kGreen+2); s2->SetLineColor(kGreen+2);

        a1->SetTitle(title.c_str());
        a1->SetMaximum(a1->GetMaximum() * 1.25);
        a1->SetMinimum(0);

        a1->Draw("E");
        c1->Draw("E SAME");
        s1->Draw("E SAME");
        s2->Draw("E SAME");

        TLegend* leg = new TLegend(0.45, 0.62, 0.93, 0.90);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.030);
        leg->AddEntry(a1, "M(#pi^{+}#pi^{-}e^{+}e^{-}) all", "lpe");
        leg->AddEntry(c1, "CB (2#sqrt{N_{++}N_{--}})", "lpe");
        leg->AddEntry(s1, "Signal (full sample)", "lpe");
        leg->AddEntry(s2, Form("Signal (#pi^{+}#pi^{-}#pi^{0}-tagged) #times %.3g", scale), "lpe");
        leg->Draw();

        cv->Update();
        pu.save(cv, "pippimepemg_overlay_" + tag);

        std::cout << "[overlay " << tag << "] probe ["
                  << kProbeLo << ", " << kProbeHi << "] GeV/c^2:"
                  << "  signal(full)=" << y_full
                  << "  signal(tagged)=" << y_tagged
                  << "  scale=" << scale << "\n";
    };

    // -----------------------------------------------------------------
    // SECTION B: STANDALONE — M(pippimepemg) all/CB/signal, no overlay
    // -----------------------------------------------------------------
    auto drawStandalone = [&](const std::string& cut_base, const std::string& title,
                              const std::string& tag, int nbins) {
        std::string cut_g = cut_base.empty() ? std::string("pippimepemg_pass==1")
                                             : cut_base + " && pippimepemg_pass==1";
        TH1D *a, *c, *s;
        std::tie(a, c, s) = pu.drawSignal(
            NT, "m_pippimepemg", nbins, 0.2, 1.4, cut_g,
            ";M_{#pi^{+}#pi^{-}#pi^{0}} [GeV/c^{2}];Counts");
        auto* cv = pu.drawTriple(a, c, s, title, "c_pippimepemg_" + tag);
        pu.save(cv, "pippimepemg_" + tag);
        pippimepemg_print(("pippimepemg standalone " + tag).c_str(), a, c, s);
    };

    // -----------------------------------------------------------------
    // Cut configurations (mirror those in pippimepem_spectra_cor.C)
    // -----------------------------------------------------------------
    auto run = [&](const std::string& cut_base, const std::string& title,
                   const std::string& tag, int nbins = 200) {
        drawOverlay(cut_base, title + " (overlay)", tag, nbins);
        drawStandalone(cut_base, title + " (#pi^{+}#pi^{-}#pi^{0})", tag, nbins);
    };

    run("",                                 "M (no cut)",                                   "base",                  200);
    run("sel_pass==1",                      "M (after selection)",                          "selected",              200);
    run("sel_pass==1 && cut2d_pass==1",     "M (after selection + cut_2d)",                 "selected_cut2d",        200);

    // Slice configurations — same MM(pippimepem) windows as the main macro
    auto runSlice = [&](double lo, double hi, const std::string& tag_suffix, bool with_cut2d) {
        std::ostringstream cut;
        cut << "sel_pass==1 && mm_pippimepem>=" << lo
            << " && mm_pippimepem<=" << hi;
        if (with_cut2d) cut << " && cut2d_pass==1";

        std::ostringstream title;
        title << "M, sel" << (with_cut2d ? "+cut_2d" : "")
              << ", MM #in [" << lo << ", " << hi << "] GeV/c^{2}";

        std::string tag = "slice_" + tag_suffix + (with_cut2d ? "_cut2d" : "");
        run(cut.str(), title.str(), tag, 100);
    };

    // Slices without cut_2d
    runSlice(2.0, 2.2, "20_22", false);
    runSlice(2.2, 2.4, "22_24", false);
    runSlice(2.4, 2.6, "24_26", false);
    runSlice(2.6, 2.8, "26_28", false);
    runSlice(2.8, 3.0, "28_30", false);

    // Slices with cut_2d
    runSlice(2.0, 2.2, "20_22", true);
    runSlice(2.2, 2.4, "22_24", true);
    runSlice(2.4, 2.6, "24_26", true);
    runSlice(2.6, 2.8, "26_28", true);
    runSlice(2.8, 3.0, "28_30", true);

    std::cout << "\nDone. Plots saved to plots/output/pippimepemg_*\n";
}
