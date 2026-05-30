// mass_ee_padnum_global_flatness_sim.C — 2D scan of (ep_rich_padnum,
// em_rich_padnum) cuts on the SIM sample, scored by a SINGLE global FOM
// computed over a fixed fit range.
//
// Motivation
// ----------
// Earlier per-window scans were judged unphysical: leg-quality cuts must
// not depend on m_ee. This macro replaces the 4-window FOM with one
// global figure of merit:
//
//   FOM = chi2/ndf of a constant fit to N_PT2 / N_PT3 vs m_ee_sim
//         inside [0.1, 1.1] GeV/c².
//
// Fit range rationale:
//   • lower bound 0.1 GeV/c²  — excludes the steeply falling π⁰ region.
//   • upper bound 1.1 GeV/c²  — excludes the high-mass tail with huge
//                                statistical errors.
//
// Convention
// ----------
//   ep_rich_padnum >= ep_cut   AND   em_rich_padnum >= em_cut
//   Default grid: cut ∈ [5, 25] (21 values per leg → 441 combinations).
//
// Outputs (2 canvases × {pdf,png}):
//   plots/output/scan_padnum_global_flatness_sim_heatmap.{pdf,png}
//   plots/output/scan_padnum_global_flatness_sim_best.{pdf,png}
//
// Heatmap: single panel χ²/ndf landscape with red marker at minimum.
// Best canvas: 1×2 — left m_ee_sim spectra (PT3 black, PT2 red) at best
// cuts, log Y; right N_PT2/N_PT3 ratio with constant fit drawn in
// [0.1, 1.1].
//
// Usage:
//   root -l -b -q plots/mass_ee_padnum_global_flatness_sim.C
//   root -l -b -q 'plots/mass_ee_padnum_global_flatness_sim.C(5, 25)'
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

// -------- Single global best point --------
struct Best {
    int    ep_cut    = 0;
    int    em_cut    = 0;
    double a         = 0.0;
    double a_err     = 0.0;
    double chi2_ndf  = std::numeric_limits<double>::infinity();
    int    ndf       = 0;
    double sum_w_pt3 = 0.0;
    double sum_w_pt2 = 0.0;
};

// -------- Constant-fit helper --------
struct FitResult { double a; double e; double chi2_ndf; int ndf; };

FitResult fitConstRange(TH1D* h, double xlo, double xhi) {
    const std::string fname = std::string(h->GetName()) + "_fc";
    TF1* f = new TF1(fname.c_str(), "[0]", xlo, xhi);
    f->SetLineColor(kBlue + 2);
    f->SetLineWidth(2);
    h->Fit(f, "RQ0");
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

void mass_ee_padnum_global_flatness_sim(int cut_min = 5, int cut_max = 25) {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetErrorX(0.0);
    gStyle->SetPalette(kBird);

    const int n_cuts = cut_max - cut_min + 1;
    std::cout << "Scanning ep_rich_padnum, em_rich_padnum on ["
              << cut_min << ", " << cut_max << "] step 1  (n_cuts="
              << n_cuts << ", " << n_cuts * n_cuts << " combinations)\n"
              << "Convention: padnum >= cut_value passes.\n"
              << "FOM: chi2/ndf of constant fit to N_PT2/N_PT3 vs m_ee_sim "
                 "in [" << kFitLo << ", " << kFitHi
              << "] GeV/c² (LOWER = FLATTER).\n";

    // ============================================================================
    // Allocate histograms — one per (ep_cut, em_cut, trigger).
    // ============================================================================
    using HistMat3 = std::vector<std::vector<std::array<TH1D*, NUM_TRG>>>;
    HistMat3 h_epem(n_cuts,
        std::vector<std::array<TH1D*, NUM_TRG>>(n_cuts));
    for (int i = 0; i < n_cuts; ++i) {
        for (int j = 0; j < n_cuts; ++j) {
            for (int t = 0; t < NUM_TRG; ++t) {
                h_epem[i][j][t] = new TH1D(
                    Form("h_epem_e%d_m%d_t%d", i + cut_min, j + cut_min, t),
                    "", kNb, kXmin, kXmax);
                h_epem[i][j][t]->Sumw2();
                h_epem[i][j][t]->SetDirectory(nullptr);
            }
        }
    }

    // Baseline (no padnum requirement) histograms for context comparison.
    std::array<TH1D*, NUM_TRG> h_base;
    for (int t = 0; t < NUM_TRG; ++t) {
        h_base[t] = new TH1D(Form("h_base_t%d", t), "",
                             kNb, kXmin, kXmax);
        h_base[t]->Sumw2();
        h_base[t]->SetDirectory(nullptr);
    }

    // ============================================================================
    // Event loop — single pass through dilepton_nt_cor (truth m_ee_sim).
    // ============================================================================
    std::cout << "\nReading output_epem_sim.root → dilepton_nt_cor...\n";
    Long64_t n_total = 0, n_use = 0;
    double   sum_w_baseline_pt3 = 0.0;
    double   sum_w_baseline_pt2 = 0.0;
    {
        TFile f("output_epem_sim.root", "READ");
        if (f.IsZombie()) {
            std::cerr << "Cannot open output_epem_sim.root\n"; return;
        }
        TTreeReader r("dilepton_nt_cor", &f);
        TTreeReaderValue<float> mee  (r, "m_ee_sim");
        TTreeReaderValue<float> v_pt3(r, "pt3");
        TTreeReaderValue<float> v_pt2(r, "pt2");
        TTreeReaderValue<float> v_w  (r, "sim_genweight");
        TTreeReaderValue<float> v_ep (r, "ep_rich_padnum");
        TTreeReaderValue<float> v_em (r, "em_rich_padnum");

        while (r.Next()) {
            ++n_total;
            const bool is_pt3 = (static_cast<int>(*v_pt3) == 1);
            const bool is_pt2 = (static_cast<int>(*v_pt2) == 1);
            if (!is_pt3 && !is_pt2) continue;

            const double m = *mee;
            const double w = *v_w;

            // Baseline (no padnum cut, just trigger).
            if (is_pt3) { h_base[TRG_PT3]->Fill(m, w); sum_w_baseline_pt3 += w; }
            if (is_pt2) { h_base[TRG_PT2]->Fill(m, w); sum_w_baseline_pt2 += w; }

            const int P_ep = static_cast<int>(*v_ep);
            const int P_em = static_cast<int>(*v_em);
            const int max_e = std::min(cut_max, P_ep);
            const int max_m = std::min(cut_max, P_em);
            if (max_e < cut_min || max_m < cut_min) continue;
            ++n_use;
            for (int e = cut_min; e <= max_e; ++e) {
                auto& row = h_epem[e - cut_min];
                for (int mc = cut_min; mc <= max_m; ++mc) {
                    if (is_pt3) row[mc - cut_min][TRG_PT3]->Fill(m, w);
                    if (is_pt2) row[mc - cut_min][TRG_PT2]->Fill(m, w);
                }
            }
        }
    }
    std::cout << "  scanned " << n_total << ", kept (passed lowest cut) "
              << n_use
              << " (" << (100.0 * n_use / std::max<Long64_t>(1, n_total))
              << " %)\n";

    // ============================================================================
    // χ²/ndf landscape — single global FOM per (ep_cut, em_cut).
    // ============================================================================
    std::cout << "\nComputing global flatness landscape (constant fit chi2/ndf in ["
              << kFitLo << ", " << kFitHi << "])...\n";

    TH2D* h_chi2 = new TH2D(
        "h_chi2_global",
        "Padnum flatness landscape — single fit in [0.1, 1.1] GeV/c^{2};"
        "ep_rich_padnum #geq;em_rich_padnum #geq;#chi^{2}/ndf",
        n_cuts, cut_min - 0.5, cut_max + 0.5,
        n_cuts, cut_min - 0.5, cut_max + 0.5);
    h_chi2->SetDirectory(nullptr);

    Best best;

    for (int ei = 0; ei < n_cuts; ++ei) {
        for (int mi = 0; mi < n_cuts; ++mi) {
            TH1D* h_PT3 = h_epem[ei][mi][TRG_PT3];
            TH1D* h_PT2 = h_epem[ei][mi][TRG_PT2];
            TH1D* r_ratio = makeRatio(h_PT2, h_PT3,
                Form("r_tmp_e%d_m%d", ei, mi));
            const FitResult fr = fitConstRange(r_ratio, kFitLo, kFitHi);
            h_chi2->SetBinContent(ei + 1, mi + 1, fr.chi2_ndf);

            if (fr.chi2_ndf < best.chi2_ndf && fr.ndf > 0) {
                best.ep_cut    = cut_min + ei;
                best.em_cut    = cut_min + mi;
                best.a         = fr.a;
                best.a_err     = fr.e;
                best.chi2_ndf  = fr.chi2_ndf;
                best.ndf       = fr.ndf;
                best.sum_w_pt3 = h_PT3->Integral();
                best.sum_w_pt2 = h_PT2->Integral();
            }
            delete r_ratio;
        }
    }

    // Baseline FOM (no-cut) for context.
    TH1D* r_base = makeRatio(h_base[TRG_PT2], h_base[TRG_PT3], "r_baseline");
    const FitResult f_base = fitConstRange(r_base, kFitLo, kFitHi);

    // Boundary check.
    const bool at_boundary =
        (best.ep_cut == cut_min || best.ep_cut == cut_max ||
         best.em_cut == cut_min || best.em_cut == cut_max);

    // ============================================================================
    // Summary
    // ============================================================================
    const double retention_pt3 = (sum_w_baseline_pt3 > 0)
        ? 100.0 * best.sum_w_pt3 / sum_w_baseline_pt3 : 0.0;
    const double retention_pt2 = (sum_w_baseline_pt2 > 0)
        ? 100.0 * best.sum_w_pt2 / sum_w_baseline_pt2 : 0.0;

    std::cout << "\n=== Global flatness FOM summary ===\n"
              << "  fit range: [" << kFitLo << ", " << kFitHi << "] GeV/c^2\n"
              << "  baseline (no padnum cut):\n"
              << "    a = " << f_base.a << " +/- " << f_base.e
              << "   chi2/ndf = " << f_base.chi2_ndf
              << " / " << f_base.ndf << "\n"
              << "  best:\n"
              << "    ep_cut = " << best.ep_cut
              << "   em_cut = " << best.em_cut << "\n"
              << "    a = " << best.a << " +/- " << best.a_err
              << "   chi2/ndf = " << best.chi2_ndf
              << " / " << best.ndf << "\n"
              << "  retention (sum_weight at best / baseline):\n"
              << "    PT3: " << retention_pt3 << " %\n"
              << "    PT2: " << retention_pt2 << " %\n";
    if (at_boundary) {
        std::cout << "  WARNING: best cut lands at grid boundary "
                  << "(ep=" << best.ep_cut << ", em=" << best.em_cut
                  << ", grid=[" << cut_min << ", " << cut_max
                  << "]). Scan range is too narrow — extend it.\n";
    } else {
        std::cout << "  best is interior — scan range is OK.\n";
    }
    std::cout << "\n";

    gSystem->mkdir("plots/output", true);
    delete r_base;

    // ============================================================================
    // OUTPUT 1 — Heatmap canvas (single panel)
    // ============================================================================
    TCanvas* c_heat = new TCanvas("c_padnum_global_flatness_sim_heatmap",
        "Padnum global flatness landscape (sim)", 900, 760);
    c_heat->cd();
    gPad->SetMargin(0.13, 0.16, 0.13, 0.10);
    gPad->SetLogz(true);
    h_chi2->Draw("COLZ");
    TMarker* mk = new TMarker(best.ep_cut, best.em_cut, 29);
    mk->SetMarkerSize(2.4);
    mk->SetMarkerColor(kRed);
    mk->Draw();

    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.030);
    tx.DrawLatex(0.16, 0.93,
        Form("best: (%d, %d)  #chi^{2}/ndf = %.2f / %d   a = %.4f #pm %.4f",
             best.ep_cut, best.em_cut,
             best.chi2_ndf, best.ndf, best.a, best.a_err));
    tx.SetTextSize(0.026);
    tx.SetTextColor(kGray + 3);
    tx.DrawLatex(0.16, 0.89,
        Form("baseline (no cut): #chi^{2}/ndf = %.2f / %d   a = %.4f #pm %.4f",
             f_base.chi2_ndf, f_base.ndf, f_base.a, f_base.e));
    if (at_boundary) {
        tx.SetTextColor(kRed + 1);
        tx.SetTextSize(0.030);
        tx.DrawLatex(0.16, 0.85, "WARNING: minimum at grid boundary");
    }

    c_heat->SaveAs("plots/output/scan_padnum_global_flatness_sim_heatmap.pdf");
    c_heat->SaveAs("plots/output/scan_padnum_global_flatness_sim_heatmap.png");
    std::cout << "Saved: plots/output/scan_padnum_global_flatness_sim_heatmap.{pdf,png}\n";

    // ============================================================================
    // OUTPUT 2 — Best-case canvas (1×2)
    // ============================================================================
    {
        const int ei = best.ep_cut - cut_min;
        const int mi = best.em_cut - cut_min;

        TH1D* h_PT3 = static_cast<TH1D*>(h_epem[ei][mi][TRG_PT3]->Clone("best_PT3"));
        TH1D* h_PT2 = static_cast<TH1D*>(h_epem[ei][mi][TRG_PT2]->Clone("best_PT2"));
        h_PT3->SetDirectory(nullptr);
        h_PT2->SetDirectory(nullptr);
        TH1D* r_ratio = makeRatio(h_PT2, h_PT3, "best_ratio");

        TF1* f_fit = new TF1("f_best", "[0]", kFitLo, kFitHi);
        f_fit->SetLineColor(kBlue + 2);
        f_fit->SetLineWidth(3);
        r_ratio->Fit(f_fit, "RQ0");

        styleSpec(h_PT3, kBlack,    21);
        styleSpec(h_PT2, kRed + 1,  20);
        styleSpec(r_ratio, kBlue + 1, 20);

        TCanvas* c = new TCanvas("c_padnum_global_flatness_best",
            Form("Best — ep_padnum>=%d  em_padnum>=%d",
                 best.ep_cut, best.em_cut),
            1700, 600);
        c->Divide(2, 1, 0.001, 0.001);

        // -- Left: m_ee_sim spectrum (log Y) at best cuts --
        c->cd(1);
        gPad->SetLogy(true);
        gPad->SetMargin(0.13, 0.04, 0.13, 0.10);
        const double y_max = std::max(h_PT3->GetMaximum(), h_PT2->GetMaximum());
        h_PT3->SetTitle(Form(
            "m_{ee}^{sim} at best cuts ep#geq%d em#geq%d;"
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
        leg->AddEntry(h_PT3, Form("PT3  (#sum w = %.0f)", h_PT3->Integral()), "lpe");
        leg->AddEntry(h_PT2, Form("PT2  (#sum w = %.0f)", h_PT2->Integral()), "lpe");
        leg->AddEntry(box,   Form("fit range [%.2f, %.2f]", kFitLo, kFitHi), "f");
        leg->Draw();

        // -- Right: ratio N_PT2/N_PT3 with constant fit drawn --
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

        TLatex tx2; tx2.SetNDC(); tx2.SetTextAlign(11); tx2.SetTextSize(0.034);
        tx2.DrawLatex(0.16, 0.90,
            Form("best ep#geq%d  em#geq%d", best.ep_cut, best.em_cut));
        tx2.SetTextColor(kBlue + 2);
        tx2.SetTextSize(0.030);
        tx2.DrawLatex(0.16, 0.85,
            Form("[%.2f, %.2f]: a = %.4f #pm %.4f   #chi^{2}/ndf = %.2f / %d",
                 kFitLo, kFitHi, best.a, best.a_err,
                 best.chi2_ndf, best.ndf));
        tx2.SetTextColor(kGray + 3);
        tx2.DrawLatex(0.16, 0.81,
            Form("baseline (no cut): a = %.4f #pm %.4f   #chi^{2}/ndf = %.2f / %d",
                 f_base.a, f_base.e, f_base.chi2_ndf, f_base.ndf));

        c->SaveAs("plots/output/scan_padnum_global_flatness_sim_best.pdf");
        c->SaveAs("plots/output/scan_padnum_global_flatness_sim_best.png");
        std::cout << "Saved: plots/output/scan_padnum_global_flatness_sim_best.{pdf,png}\n";
    }

    std::cout << "\nDone. 2 output canvases saved under plots/output/.\n";
}
