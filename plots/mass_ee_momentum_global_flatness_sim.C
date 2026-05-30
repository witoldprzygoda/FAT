// mass_ee_momentum_global_flatness_sim.C — 2D scan of (ep_p_rec, em_p_rec)
// momentum cuts on the SIM sample using a SINGLE global figure of merit:
//
//   FOM = chi2/ndf of a constant fit to N_PT2/N_PT3 vs m_ee_sim
//         inside the range [0.1, 1.1] GeV/c²
//
// Rationale for the fit range:
//   • lower bound 0.1: excludes the steeply-falling π⁰ region (not the
//     trigger story).
//   • upper bound 1.1: excludes the high-mass tail with huge errors
//     (also not the trigger story).
//
// What is scanned (ASYMMETRIC ep / em, STRICT >):
//   ep_p_rec > ep_cut    AND    em_p_rec > em_cut
//   Default grid: cut ∈ [0, 300] step 10 MeV/c
//                 (31 values per leg → 961 combinations).
//
// NO pre-filter is applied — only the momentum cuts themselves are used.
//
// Outputs (4 files total):
//   plots/output/scan_momentum_global_flatness_sim_heatmap.{pdf,png}
//   plots/output/scan_momentum_global_flatness_sim_best.{pdf,png}
//
// Best-case canvas is 1×2:
//   Left:  m_ee_sim spectrum at best cuts — PT3 (black filled square) overlaid
//          with PT2 (red filled circle), log Y.
//   Right: ratio N_PT2/N_PT3 with constant fit drawn (line in [0.1, 1.1]),
//          linear Y [0, 3], y=1 ref line, annotation listing a ± err,
//          chi²/ndf.
//
// Usage:
//   root -l -b -q plots/mass_ee_momentum_global_flatness_sim.C
//   root -l -b -q 'plots/mass_ee_momentum_global_flatness_sim.C(0, 300, 10)'
//
// @author Witold Przygoda (witold.przygoda@uj.edu.pl)
// @date 2026

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TMarker.h>
#include <TBox.h>
#include <TStyle.h>
#include <TSystem.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {

// -------- m_ee binning (= 20 MeV) and fit limits --------
constexpr int    kNb    = 70;
constexpr double kXmin  = 0.0;
constexpr double kXmax  = 1.4;
constexpr double kFitLo = 0.10;
constexpr double kFitHi = 1.10;

// -------- Trigger encoding --------
enum Trigger : int { TRG_PT2 = 0, TRG_PT3 = 1, NUM_TRG = 2 };
const std::array<std::string, NUM_TRG> kTrgName  = {"PT2", "PT3"};
const std::array<std::string, NUM_TRG> kTrgLabel = {"PT2 trigger",
                                                    "PT3 trigger"};

// -------- Best-point bookkeeping (single global) --------
struct BestPoint {
    int    ep_cut    = 0;     // MeV/c
    int    em_cut    = 0;     // MeV/c
    double a         = 0.0;
    double a_err     = 0.0;
    double chi2_ndf  = std::numeric_limits<double>::infinity();
    int    ndf       = 0;
    double sum_w     = 0.0;   // sum of weights retained (PT2+PT3, in fit range)
};

// -------- Thousands-separator helper --------
inline std::string fmtApos(double v) {
    const long long iabs = static_cast<long long>(std::llround(std::abs(v)));
    std::string s = std::to_string(iabs);
    for (int pos = static_cast<int>(s.size()) - 3; pos > 0; pos -= 3)
        s.insert(pos, "'");
    return (v < 0 ? "-" + s : s);
}

// -------- Constant-fit helper --------
struct FitResult { double a; double e; double chi2_ndf; int ndf; };

FitResult fitConstRange(TH1D* h, double xlo, double xhi) {
    const std::string fname = std::string(h->GetName()) + "_fc";
    TF1* f = new TF1(fname.c_str(), "[0]", xlo, xhi);
    f->SetLineColor(kBlue + 2);
    f->SetLineWidth(2);
    h->Fit(f, "RQ0");                  // "0" = don't draw — we attach later
    FitResult r;
    r.a        = f->GetParameter(0);
    r.e        = f->GetParError(0);
    r.ndf      = f->GetNDF();
    r.chi2_ndf = (r.ndf > 0) ? f->GetChisquare() / r.ndf : 0.0;
    return r;
}

TH1D* makeRatio(TH1D* num, TH1D* den, const std::string& name) {
    TH1D* r = static_cast<TH1D*>(num->Clone(name.c_str()));
    r->SetDirectory(nullptr);
    r->Reset();
    for (int b = 1; b <= num->GetNbinsX(); ++b) {
        const double n  = num->GetBinContent(b);
        const double en = num->GetBinError(b);
        const double d  = den->GetBinContent(b);
        const double ed = den->GetBinError(b);
        if (d != 0.0 && std::isfinite(d) && std::isfinite(n)) {
            const double val = n / d;
            const double rel = std::sqrt(
                (n != 0.0 ? std::pow(en / n, 2) : 0.0) +
                std::pow(ed / d, 2));
            r->SetBinContent(b, val);
            r->SetBinError  (b, std::abs(val) * rel);
        }
    }
    return r;
}

void styleSpec(TH1D* h, Color_t color, Style_t marker) {
    h->SetLineColor(color);
    h->SetMarkerColor(color);
    h->SetMarkerStyle(marker);
    h->SetMarkerSize(0.8);
    h->SetLineWidth(2);
    h->GetXaxis()->SetTitleSize(0.050);
    h->GetYaxis()->SetTitleSize(0.050);
    h->GetXaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetLabelSize(0.045);
    h->GetXaxis()->SetTitleOffset(1.05);
    h->GetYaxis()->SetTitleOffset(1.20);
}

}  // anonymous namespace

void mass_ee_momentum_global_flatness_sim(int cut_min = 0,
                                          int cut_max = 300,
                                          int cut_step = 10) {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetErrorX(0.0);
    gStyle->SetPalette(kBird);

    const int n_cuts = (cut_max - cut_min) / cut_step + 1;
    auto cutVal = [cut_min, cut_step](int i) {
        return cut_min + i * cut_step;
    };

    std::cout << "Scanning ep_p_rec, em_p_rec on ["
              << cut_min << ", " << cut_max << "] step " << cut_step
              << " MeV/c  (n_cuts=" << n_cuts << ", "
              << n_cuts * n_cuts << " combinations)\n"
              << "Convention: STRICT inequality, p_rec > cut_value passes.\n"
              << "FOM: single chi2/ndf of constant fit to N_PT2/N_PT3 "
              << "vs m_ee_sim in [" << kFitLo << ", " << kFitHi
              << "] GeV/c² (LOWER = FLATTER).\n";

    // ============================================================================
    // Allocate histograms: h_epem[ei][mi][trg] over 20 MeV mass grid.
    // ============================================================================
    using HistMat3 = std::vector<std::vector<std::array<TH1D*, NUM_TRG>>>;
    HistMat3 h_epem(n_cuts,
        std::vector<std::array<TH1D*, NUM_TRG>>(n_cuts));
    for (int i = 0; i < n_cuts; ++i) {
        for (int j = 0; j < n_cuts; ++j) {
            for (int t = 0; t < NUM_TRG; ++t) {
                h_epem[i][j][t] = new TH1D(
                    Form("h_epem_e%d_m%d_t%d",
                         cutVal(i), cutVal(j), t),
                    "", kNb, kXmin, kXmax);
                h_epem[i][j][t]->Sumw2();
                h_epem[i][j][t]->SetDirectory(nullptr);
            }
        }
    }

    // ============================================================================
    // Event loop — single pass through dilepton_nt_cor (truth m_ee_sim)
    // NO pre-filter; only the momentum cuts are imposed inside the inner loop.
    // ============================================================================
    std::cout << "\nReading output_epem_sim.root → dilepton_nt_cor...\n";
    Long64_t n_total = 0, n_use = 0;
    double sum_w_baseline = 0.0;  // sum of weights with no momentum cut
    {
        TFile f("output_epem_sim.root", "READ");
        if (f.IsZombie()) {
            std::cerr << "Cannot open output_epem_sim.root\n"; return;
        }
        TTreeReader r("dilepton_nt_cor", &f);
        TTreeReaderValue<float> mee   (r, "m_ee_sim");
        TTreeReaderValue<float> v_pt3 (r, "pt3");
        TTreeReaderValue<float> v_pt2 (r, "pt2");
        TTreeReaderValue<float> v_w   (r, "sim_genweight");
        TTreeReaderValue<float> v_ep_p(r, "ep_p_rec");
        TTreeReaderValue<float> v_em_p(r, "em_p_rec");

        while (r.Next()) {
            ++n_total;
            const bool is_pt3 = (static_cast<int>(*v_pt3) == 1);
            const bool is_pt2 = (static_cast<int>(*v_pt2) == 1);
            if (!is_pt3 && !is_pt2) continue;

            const float P_ep = *v_ep_p;
            const float P_em = *v_em_p;

            // STRICT inequality: cut_value passes if P > cut_value.
            // Max passing cut index k: cutVal(k) < P  ⇔ k < (P - cut_min)/step.
            // Hence max_e = ceil((P - cut_min)/step) - 1.
            int max_e = static_cast<int>(
                std::ceil((P_ep - cut_min) / (double)cut_step)) - 1;
            int max_m = static_cast<int>(
                std::ceil((P_em - cut_min) / (double)cut_step)) - 1;
            max_e = std::min(max_e, n_cuts - 1);
            max_m = std::min(max_m, n_cuts - 1);
            if (max_e < 0 || max_m < 0) continue;

            ++n_use;
            const double m = *mee;
            const double w = *v_w;

            // Track baseline sum of weights (cut = cut_min on both legs).
            // Baseline includes any event with max_e>=0 && max_m>=0.
            if (m >= kFitLo && m < kFitHi) {
                sum_w_baseline += w;
            }

            for (int e = 0; e <= max_e; ++e) {
                auto& row = h_epem[e];
                for (int mc = 0; mc <= max_m; ++mc) {
                    if (is_pt3) row[mc][TRG_PT3]->Fill(m, w);
                    if (is_pt2) row[mc][TRG_PT2]->Fill(m, w);
                }
            }
        }
    }
    std::cout << "  scanned " << n_total << ", kept " << n_use
              << " (" << (100.0 * n_use / std::max<Long64_t>(1, n_total))
              << " %)\n";

    // ============================================================================
    // χ²/ndf landscape — single global FOM
    // ============================================================================
    std::cout << "\nComputing global flatness landscape (chi2/ndf in ["
              << kFitLo << ", " << kFitHi << "])...\n";

    TH2D* h_chi2 = new TH2D(
        "h_chi2_global",
        Form("Momentum flatness landscape - single fit in [%.1f, %.1f] GeV/c^{2};"
             "ep_p_rec > [MeV/c];em_p_rec > [MeV/c];#chi^{2}/ndf",
             kFitLo, kFitHi),
        n_cuts, cut_min - 0.5 * cut_step, cut_max + 0.5 * cut_step,
        n_cuts, cut_min - 0.5 * cut_step, cut_max + 0.5 * cut_step);
    h_chi2->SetDirectory(nullptr);

    BestPoint best;
    BestPoint baseline;          // value at (cut_min, cut_min) for context

    for (int ei = 0; ei < n_cuts; ++ei) {
        for (int mi = 0; mi < n_cuts; ++mi) {
            TH1D* h_PT3 = h_epem[ei][mi][TRG_PT3];
            TH1D* h_PT2 = h_epem[ei][mi][TRG_PT2];
            TH1D* r_ratio = makeRatio(h_PT2, h_PT3,
                Form("r_tmp_e%d_m%d", ei, mi));

            const FitResult fr = fitConstRange(r_ratio, kFitLo, kFitHi);
            h_chi2->SetBinContent(ei + 1, mi + 1, fr.chi2_ndf);

            // Sum of weights in fit range at this point (PT3 used as
            // representative — same sample basis as PT2 just different
            // trigger acceptance).
            double sw = 0.0;
            for (int b = 1; b <= h_PT3->GetNbinsX(); ++b) {
                const double xc = h_PT3->GetBinCenter(b);
                if (xc >= kFitLo && xc < kFitHi) {
                    sw += h_PT3->GetBinContent(b) +
                          h_PT2->GetBinContent(b);
                }
            }

            if (ei == 0 && mi == 0) {
                baseline.ep_cut   = cutVal(ei);
                baseline.em_cut   = cutVal(mi);
                baseline.a        = fr.a;
                baseline.a_err    = fr.e;
                baseline.chi2_ndf = fr.chi2_ndf;
                baseline.ndf      = fr.ndf;
                baseline.sum_w    = sw;
            }

            if (fr.chi2_ndf < best.chi2_ndf && fr.ndf > 0) {
                best.ep_cut   = cutVal(ei);
                best.em_cut   = cutVal(mi);
                best.a        = fr.a;
                best.a_err    = fr.e;
                best.chi2_ndf = fr.chi2_ndf;
                best.ndf      = fr.ndf;
                best.sum_w    = sw;
            }
            delete r_ratio;
        }
    }

    // ============================================================================
    // Summary
    // ============================================================================
    const bool at_boundary =
        (best.ep_cut == cut_min || best.ep_cut == cut_max ||
         best.em_cut == cut_min || best.em_cut == cut_max);

    const double retention = (baseline.sum_w > 0.0)
        ? 100.0 * best.sum_w / baseline.sum_w
        : 0.0;

    std::cout << "\n=== Global FOM: best (ep_cut, em_cut) — "
              << "minimize chi2/ndf in [" << kFitLo << ", " << kFitHi
              << "] GeV/c² ===\n\n";
    std::cout << "  baseline (ep_cut="  << baseline.ep_cut
              << ", em_cut="  << baseline.em_cut << "):\n"
              << "      a = " << baseline.a << " +/- " << baseline.a_err
              << "   chi2/ndf = " << baseline.chi2_ndf
              << " / " << baseline.ndf << "\n\n";
    std::cout << "  BEST     (ep_cut="  << best.ep_cut
              << ", em_cut="  << best.em_cut << "):\n"
              << "      a = " << best.a << " +/- " << best.a_err
              << "   chi2/ndf = " << best.chi2_ndf
              << " / " << best.ndf << "\n"
              << "      sample retention vs baseline (in fit range): "
              << retention << " %\n";
    if (at_boundary) {
        std::cout << "\n  WARNING: best lands at scan-grid boundary "
                  << "(ep=" << best.ep_cut << ", em=" << best.em_cut
                  << "; range [" << cut_min << ", " << cut_max
                  << "]). Consider widening the scan range.\n";
    }
    std::cout << "\n";

    gSystem->mkdir("plots/output", true);

    // ============================================================================
    // OUTPUT 1 — Heatmap canvas (single panel)
    // ============================================================================
    TCanvas* c_heat = new TCanvas("c_momentum_global_flatness_heatmap",
        "Momentum global flatness landscape (sim)", 900, 800);
    c_heat->SetMargin(0.13, 0.16, 0.13, 0.10);
    c_heat->SetLogz(true);
    h_chi2->Draw("COLZ");

    TMarker* mk = new TMarker(best.ep_cut, best.em_cut, 29);
    mk->SetMarkerSize(2.5);
    mk->SetMarkerColor(kRed);
    mk->Draw();

    {
        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.030);
        tx.DrawLatex(0.16, 0.93,
            Form("best: (ep>%d, em>%d) MeV/c   #chi^{2}/ndf = %.2f / %d",
                 best.ep_cut, best.em_cut, best.chi2_ndf, best.ndf));
        tx.SetTextSize(0.026);
        tx.SetTextColor(kGray + 3);
        tx.DrawLatex(0.16, 0.89,
            Form("baseline (ep>%d, em>%d): #chi^{2}/ndf = %.2f / %d",
                 baseline.ep_cut, baseline.em_cut,
                 baseline.chi2_ndf, baseline.ndf));
        if (at_boundary) {
            tx.SetTextColor(kRed + 1);
            tx.SetTextSize(0.030);
            tx.DrawLatex(0.16, 0.85, "WARNING: best at grid boundary");
        }
    }
    c_heat->SaveAs("plots/output/scan_momentum_global_flatness_sim_heatmap.pdf");
    c_heat->SaveAs("plots/output/scan_momentum_global_flatness_sim_heatmap.png");
    std::cout << "Saved: plots/output/scan_momentum_global_flatness_sim_heatmap.{pdf,png}\n";

    // ============================================================================
    // OUTPUT 2 — Best-case canvas (1×2)
    // ============================================================================
    const int ei = (best.ep_cut - cut_min) / cut_step;
    const int mi = (best.em_cut - cut_min) / cut_step;

    TH1D* h_PT3 = static_cast<TH1D*>(
        h_epem[ei][mi][TRG_PT3]->Clone("c_PT3_global"));
    TH1D* h_PT2 = static_cast<TH1D*>(
        h_epem[ei][mi][TRG_PT2]->Clone("c_PT2_global"));
    h_PT3->SetDirectory(nullptr);
    h_PT2->SetDirectory(nullptr);
    TH1D* r_ratio = makeRatio(h_PT2, h_PT3, "c_ratio_global");

    TF1* f_fit = new TF1("f_fit_global", "[0]", kFitLo, kFitHi);
    f_fit->SetLineColor(kBlue + 2); f_fit->SetLineWidth(3);
    r_ratio->Fit(f_fit, "RQ0");

    styleSpec(h_PT3, kBlack,    21);     // filled square
    styleSpec(h_PT2, kRed + 1,  20);     // filled circle
    styleSpec(r_ratio, kBlue + 1, 20);

    TCanvas* c = new TCanvas(
        "c_momentum_global_flatness_best",
        Form("Best — ep_p_rec>%d  em_p_rec>%d MeV/c",
             best.ep_cut, best.em_cut),
        1700, 600);
    c->Divide(2, 1, 0.001, 0.001);

    // ── Left: m_ee_sim spectrum (log Y) at best cuts ──
    c->cd(1);
    gPad->SetLogy(true);
    gPad->SetMargin(0.13, 0.04, 0.13, 0.10);
    const double y_max = std::max(h_PT3->GetMaximum(),
                                  h_PT2->GetMaximum());
    h_PT3->SetTitle(Form(
        "m_{ee}^{sim} at best cuts ep>%d em>%d MeV/c;"
        "M_{e^{+}e^{-}}^{sim} [GeV/c^{2}];"
        "weighted entries / 20 MeV",
        best.ep_cut, best.em_cut));
    h_PT3->GetYaxis()->SetRangeUser(std::max(1e-2, y_max * 1e-6),
                                    y_max * 5.0);
    h_PT3->Draw("E1");
    h_PT2->Draw("E1 SAME");

    TBox* box = new TBox(kFitLo, std::max(1e-2, y_max * 1e-6),
                         kFitHi, y_max * 5.0);
    box->SetFillColorAlpha(kGreen - 9, 0.20);
    box->SetLineColor(kGreen + 2); box->SetLineStyle(2);
    box->Draw("SAME");
    h_PT3->Draw("E1 SAME");
    h_PT2->Draw("E1 SAME");

    TLegend* leg = new TLegend(0.55, 0.72, 0.95, 0.88);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.034);
    leg->AddEntry(h_PT3, Form("PT3  (#sum w = %.0f)",
                              h_PT3->Integral()), "lpe");
    leg->AddEntry(h_PT2, Form("PT2  (#sum w = %.0f)",
                              h_PT2->Integral()), "lpe");
    leg->AddEntry(box,   Form("fit range [%.2f, %.2f]",
                              kFitLo, kFitHi), "f");
    leg->Draw();

    // ── Right: ratio N_PT2/N_PT3 with constant fit drawn ──
    c->cd(2);
    gPad->SetMargin(0.13, 0.04, 0.13, 0.10);
    r_ratio->SetTitle(
        ";M_{e^{+}e^{-}}^{sim} [GeV/c^{2}];N_{PT2} / N_{PT3}");
    r_ratio->GetYaxis()->SetRangeUser(0.0, 3.0);
    r_ratio->Draw("E1");

    TBox* box2 = new TBox(kFitLo, 0.0, kFitHi, 3.0);
    box2->SetFillColorAlpha(kGreen - 9, 0.15);
    box2->SetLineColor(kGreen + 2); box2->SetLineStyle(2);
    box2->Draw("SAME");
    r_ratio->Draw("E1 SAME");
    f_fit->Draw("SAME");

    TLine* lref = new TLine(kXmin, 1.0, kXmax, 1.0);
    lref->SetLineStyle(3); lref->SetLineColor(kGray + 2); lref->Draw();

    TLatex tx; tx.SetNDC(); tx.SetTextAlign(11); tx.SetTextSize(0.034);
    tx.DrawLatex(0.16, 0.90,
        Form("best: ep>%d  em>%d MeV/c",
             best.ep_cut, best.em_cut));
    tx.SetTextColor(kBlue + 2);
    tx.SetTextSize(0.030);
    tx.DrawLatex(0.16, 0.85,
        Form("[%.2f-%.2f]: a = %.4f #pm %.4f   #chi^{2}/ndf = %.2f / %d",
             kFitLo, kFitHi, best.a, best.a_err,
             best.chi2_ndf, best.ndf));
    tx.SetTextColor(kGray + 3);
    tx.DrawLatex(0.16, 0.81,
        Form("baseline (ep>%d, em>%d): #chi^{2}/ndf = %.2f / %d",
             baseline.ep_cut, baseline.em_cut,
             baseline.chi2_ndf, baseline.ndf));
    tx.SetTextColor(kBlack);
    tx.DrawLatex(0.16, 0.77,
        Form("retention vs baseline: %.1f %%", retention));
    if (at_boundary) {
        tx.SetTextColor(kRed + 1);
        tx.SetTextSize(0.034);
        tx.DrawLatex(0.16, 0.73, "WARNING: best at grid boundary");
    }

    const std::string base =
        "plots/output/scan_momentum_global_flatness_sim_best";
    c->SaveAs((base + ".pdf").c_str());
    c->SaveAs((base + ".png").c_str());
    std::cout << "  saved: " << base << ".{pdf,png}  "
              << "(ep_cut=" << best.ep_cut
              << ", em_cut=" << best.em_cut
              << ", chi2/ndf=" << best.chi2_ndf << ")\n";

    std::cout << "\nDone. 2 output files saved under plots/output/.\n";
}
