// mass_spectra_ratio_exp.C — bin-by-bin 63 · N_PT2 / N_PT3 ratio of the
// dilepton mass spectrum for each of {all, CB, signal}.
//
// "Set 3" of three parallel views (see mass_spectra_pt3_exp.C and
// mass_spectra_pt2_exp.C for the PT3 and PT2 sets).
//
// Per-bin ratio interpretation:
//     63·N_PT2/N_PT3 ≡ "trigger correction" — the multiplicative factor
//     that converts PT3-observed counts back to the PT2-unbiased reference.
//     For "all" (epem-dominated) it sits near the per-channel <w_em> ≈ 2.33;
//     for "CB" near sqrt(<w_++>·<w_-->) ≈ 1.28; for "sig" it can swing
//     because both numerator and denominator vanish where CB ≈ all
//     (high-mass tail and π⁰/ω peak vicinities).
//
// Y range is auto-determined from the data with safety margin so the
// curves are not cropped at the top (a complaint about previous plots).
//
// Usage:  root -l -b -q plots/mass_spectra_ratio_exp.C

#include "PlotUtils.h"
#include <TCanvas.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <tuple>

namespace {

// Bin-by-bin r = 63 · num / den with conservative independent-variable
// error propagation: σ_r/r = sqrt((σ_n/n)² + (σ_d/d)²). Bins with d==0
// (or non-finite values) are left at 0 ± 0 — they are not meaningful
// trigger-correction estimates.
TH1D* makeRatio(TH1D* num, TH1D* den, const std::string& name) {
    if (!num || !den) return nullptr;
    TH1D* r = static_cast<TH1D*>(num->Clone(name.c_str()));
    r->SetDirectory(nullptr);
    r->Reset();
    const int nb = num->GetNbinsX();
    for (int b = 1; b <= nb; ++b) {
        const double n  = num->GetBinContent(b);
        const double en = num->GetBinError(b);
        const double d  = den->GetBinContent(b);
        const double ed = den->GetBinError(b);
        if (d == 0.0 || !std::isfinite(d) || !std::isfinite(n)) continue;
        const double val = 63.0 * n / d;
        const double rel_n = (n != 0.0) ? en / n : 0.0;
        const double rel_d = ed / d;
        const double er    = std::abs(val) * std::sqrt(rel_n * rel_n + rel_d * rel_d);
        r->SetBinContent(b, val);
        r->SetBinError(b, er);
    }
    return r;
}

// Robust max-search ignoring bins with huge error bars / non-finite values.
// Returns the largest |bin content| within the visible-range bin loop.
double safeMaxAbs(const TH1D* h) {
    double mx = 0.0;
    const int nb = h->GetNbinsX();
    for (int b = 1; b <= nb; ++b) {
        const double v = std::abs(h->GetBinContent(b));
        if (std::isfinite(v) && v > mx) mx = v;
    }
    return mx;
}

} // anonymous namespace

void mass_spectra_ratio_exp() {
    gStyle->SetOptStat(0);

    PlotUtils pu("output_epem_exp.root",
                 "output_epep_exp.root",
                 "output_emem_exp.root");

    TH1D *a_p3, *c_p3, *s_p3;
    std::tie(a_p3, c_p3, s_p3) = pu.getSignal("dilepton/mass_ee");
    TH1D *a_p2, *c_p2, *s_p2;
    std::tie(a_p2, c_p2, s_p2) = pu.getSignal("dilepton/mass_ee_pt2");

    if (!a_p3 || !c_p3 || !s_p3 || !a_p2 || !c_p2 || !s_p2) {
        std::cerr << "Cannot build ratio — missing histograms.\n";
        return;
    }

    TH1D* r_all = makeRatio(a_p2, a_p3, "r_all");
    TH1D* r_cb  = makeRatio(c_p2, c_p3, "r_cb");
    TH1D* r_sig = makeRatio(s_p2, s_p3, "r_sig");

    pu.styleAll(r_all);
    pu.styleCB(r_cb);
    pu.styleSignal(r_sig);

    // Auto Y range with 20% headroom (user requested generous defaults).
    const double y_top = 1.2 * std::max({
        safeMaxAbs(r_all), safeMaxAbs(r_cb), safeMaxAbs(r_sig)});
    const double y_max = (y_top > 0.0 && std::isfinite(y_top)) ? y_top : 10.0;

    const double xmin = r_all->GetXaxis()->GetXmin();
    const double xmax = r_all->GetXaxis()->GetXmax();

    TCanvas* c = new TCanvas("c_mass_ratio_pt2_pt3",
                             "63 #upoint N_{PT2}/N_{PT3} for mass_{ee}",
                             1000, 700);
    c->SetMargin(0.12, 0.05, 0.12, 0.08);

    r_all->SetTitle(
        "63 #upoint N_{PT2} / N_{PT3} (M_{e^{+}e^{-}});"
        "M_{e^{+}e^{-}} [GeV/c^{2}];"
        "63 #upoint N_{PT2} / N_{PT3}");
    r_all->GetYaxis()->SetRangeUser(0.0, y_max);
    r_all->Draw("E");
    r_cb ->Draw("E SAME");
    r_sig->Draw("E SAME");

    // Reference line at 1 (= no correction).
    TLine* l_one = new TLine(xmin, 1.0, xmax, 1.0);
    l_one->SetLineStyle(3);
    l_one->SetLineColor(kGray + 2);
    l_one->Draw();

    TLegend* leg = new TLegend(0.45, 0.66, 0.92, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.034);
    leg->AddEntry(r_all, "63 #upoint N_{PT2}/N_{PT3}  (all = e^{+}e^{-})", "lpe");
    leg->AddEntry(r_cb,  "63 #upoint N_{PT2}/N_{PT3}  (CB)",                  "lpe");
    leg->AddEntry(r_sig, "63 #upoint N_{PT2}/N_{PT3}  (sig = all - CB)",      "lpe");
    leg->AddEntry(l_one, "1.0 (no correction reference)", "l");
    leg->Draw();

    pu.save(c, "mass_ee_ratio_pt2_over_pt3");

    // ------------------------------------------------------------------
    // Same with OA > 4 deg cut, if those histograms exist.
    // ------------------------------------------------------------------
    TH1D *a2_p3, *c2_p3, *s2_p3;
    std::tie(a2_p3, c2_p3, s2_p3) = pu.getSignal("dilepton/mass_ee_after_oa");
    TH1D *a2_p2, *c2_p2, *s2_p2;
    std::tie(a2_p2, c2_p2, s2_p2) = pu.getSignal("dilepton/mass_ee_after_oa_pt2");
    if (a2_p3 && c2_p3 && s2_p3 && a2_p2 && c2_p2 && s2_p2) {
        TH1D* r2_all = makeRatio(a2_p2, a2_p3, "r2_all");
        TH1D* r2_cb  = makeRatio(c2_p2, c2_p3, "r2_cb");
        TH1D* r2_sig = makeRatio(s2_p2, s2_p3, "r2_sig");
        pu.styleAll(r2_all);
        pu.styleCB(r2_cb);
        pu.styleSignal(r2_sig);

        const double y2_top = 1.2 * std::max({
            safeMaxAbs(r2_all), safeMaxAbs(r2_cb), safeMaxAbs(r2_sig)});
        const double y2_max = (y2_top > 0.0 && std::isfinite(y2_top)) ? y2_top : 10.0;

        TCanvas* c2 = new TCanvas("c_mass_ratio_pt2_pt3_oa",
                                  "63 #upoint N_{PT2}/N_{PT3} (OA > 4 deg)",
                                  1000, 700);
        c2->SetMargin(0.12, 0.05, 0.12, 0.08);
        r2_all->SetTitle(
            "63 #upoint N_{PT2} / N_{PT3} (M_{e^{+}e^{-}}, OA > 4#circ);"
            "M_{e^{+}e^{-}} [GeV/c^{2}];"
            "63 #upoint N_{PT2} / N_{PT3}");
        r2_all->GetYaxis()->SetRangeUser(0.0, y2_max);
        r2_all->Draw("E");
        r2_cb ->Draw("E SAME");
        r2_sig->Draw("E SAME");

        TLine* l_one2 = new TLine(xmin, 1.0, xmax, 1.0);
        l_one2->SetLineStyle(3);
        l_one2->SetLineColor(kGray + 2);
        l_one2->Draw();

        TLegend* leg2 = new TLegend(0.45, 0.66, 0.92, 0.88);
        leg2->SetBorderSize(0);
        leg2->SetFillStyle(0);
        leg2->SetTextSize(0.034);
        leg2->AddEntry(r2_all, "63 #upoint N_{PT2}/N_{PT3}  (all, OA > 4#circ)", "lpe");
        leg2->AddEntry(r2_cb,  "63 #upoint N_{PT2}/N_{PT3}  (CB)",                "lpe");
        leg2->AddEntry(r2_sig, "63 #upoint N_{PT2}/N_{PT3}  (sig)",               "lpe");
        leg2->AddEntry(l_one2, "1.0",                                              "l");
        leg2->Draw();

        pu.save(c2, "mass_ee_ratio_pt2_over_pt3_oa");
    }

    std::cout << "\nDone. Plots in plots/output/.\n";
}
