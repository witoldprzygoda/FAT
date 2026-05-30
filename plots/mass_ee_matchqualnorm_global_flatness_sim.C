// mass_ee_matchqualnorm_global_flatness_sim.C — 2D scan of
// (ep_richmatchqualitynorm, em_richmatchqualitynorm) cuts on the SIM sample.
//
// SINGLE GLOBAL FOM (no per-window logic, no pre-filters, no TOP marker):
//   FOM = chi2/ndf of constant fit to N_PT2/N_PT3 vs m_ee_sim,
//         fit range [0.1, 1.1] GeV/c²  (LOWER = FLATTER).
//
// Why this range:
//   - lower bound 0.1: excludes the steeply-falling π⁰ region
//   - upper bound 1.1: excludes high-mass region with huge errors
//
// Cut convention: ep_richmatchqualitynorm <= ep_cut
//                 em_richmatchqualitynorm <= em_cut       (lower = tighter)
//   For float Q, smallest passing int cut = ceil(Q).
//
// Default grid: cut ∈ [0, 22] step 1 → 23 values per leg → 529 combinations.
//
// Outputs (2 files):
//   plots/output/scan_matchqualnorm_global_flatness_sim_heatmap.{pdf,png}
//   plots/output/scan_matchqualnorm_global_flatness_sim_best.{pdf,png}
//
// Best canvas is 1×2:
//   Left:  m_ee_sim spectrum at best cuts — PT3 (black) + PT2 (red), log Y.
//   Right: ratio N_PT2/N_PT3 with constant fit drawn in [0.1, 1.1].
//
// Usage:
//   root -l -b -q plots/mass_ee_matchqualnorm_global_flatness_sim.C
//   root -l -b -q 'plots/mass_ee_matchqualnorm_global_flatness_sim.C(0, 22)'
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

// -------- m_ee binning (= 20 MeV) --------
constexpr int    kNb    = 70;
constexpr double kXmin  = 0.0;
constexpr double kXmax  = 1.4;

// -------- Global FOM fit range --------
constexpr double kFitLo = 0.1;
constexpr double kFitHi = 1.1;

// -------- Trigger encoding --------
enum Trigger : int { TRG_PT2 = 0, TRG_PT3 = 1, NUM_TRG = 2 };

// -------- Single global best-point --------
struct BestPoint {
    int    ep_cut    = 0;
    int    em_cut    = 0;
    double a         = 0.0;
    double a_err     = 0.0;
    double chi2_ndf  = std::numeric_limits<double>::infinity();
    int    ndf       = 0;
    double sum_w_PT3 = 0.0;
    double sum_w_PT2 = 0.0;
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
    delete f;
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

void mass_ee_matchqualnorm_global_flatness_sim(int cut_min = 0, int cut_max = 22) {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetErrorX(0.0);
    gStyle->SetPalette(kBird);

    const int n_cuts = cut_max - cut_min + 1;
    std::cout << "Scanning ep_richmatchqualitynorm, em_richmatchqualitynorm on ["
              << cut_min << ", " << cut_max << "] step 1  (n_cuts="
              << n_cuts << ", " << n_cuts * n_cuts << " combinations)\n"
              << "Convention: matchqualitynorm <= cut_value passes (lower=tighter).\n"
              << "FOM: chi2/ndf of constant fit to N_PT2/N_PT3 vs m_ee_sim "
                 "in [" << kFitLo << ", " << kFitHi << "] GeV/c² "
                 "(LOWER = FLATTER).\n";

    // ============================================================================
    // Allocate histograms: h_epem[ei][mi][trg]
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

    // Baseline (no-cut) histograms — fills for every event that passes trigger,
    // independent of the matchqualitynorm cut. Used for retention denominator and
    // baseline FOM context.
    TH1D* h_base_PT3 = new TH1D("h_base_PT3", "", kNb, kXmin, kXmax);
    TH1D* h_base_PT2 = new TH1D("h_base_PT2", "", kNb, kXmin, kXmax);
    h_base_PT3->Sumw2(); h_base_PT2->Sumw2();
    h_base_PT3->SetDirectory(nullptr); h_base_PT2->SetDirectory(nullptr);

    // ============================================================================
    // Event loop — single pass through dilepton_nt_cor
    // ============================================================================
    std::cout << "\nReading output_epem_sim.root → dilepton_nt_cor...\n";
    Long64_t n_total = 0, n_use = 0;
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
        TTreeReaderValue<float> v_ep (r, "ep_richmatchqualitynorm");
        TTreeReaderValue<float> v_em (r, "em_richmatchqualitynorm");

        while (r.Next()) {
            ++n_total;
            const bool is_pt3 = (static_cast<int>(*v_pt3) == 1);
            const bool is_pt2 = (static_cast<int>(*v_pt2) == 1);
            if (!is_pt3 && !is_pt2) continue;

            const float Q_ep = *v_ep;
            const float Q_em = *v_em;
            const double m = *mee;
            const double w = *v_w;

            // Baseline (no-cut): fill once per event for retention reference.
            if (is_pt3) h_base_PT3->Fill(m, w);
            if (is_pt2) h_base_PT2->Fill(m, w);

            // Smallest passing int cut (since cut is "<="): ceil(Q), clamped at 0.
            const int min_e = std::max(cut_min, (int)std::ceil(std::max(0.0f, Q_ep)));
            const int min_m = std::max(cut_min, (int)std::ceil(std::max(0.0f, Q_em)));
            if (min_e > cut_max || min_m > cut_max) continue;
            ++n_use;

            for (int e = min_e; e <= cut_max; ++e) {
                auto& row = h_epem[e - cut_min];
                for (int mc = min_m; mc <= cut_max; ++mc) {
                    if (is_pt3) row[mc - cut_min][TRG_PT3]->Fill(m, w);
                    if (is_pt2) row[mc - cut_min][TRG_PT2]->Fill(m, w);
                }
            }
        }
    }
    std::cout << "  scanned " << n_total << ", kept " << n_use
              << " (" << (100.0 * n_use / std::max<Long64_t>(1, n_total))
              << " %)\n";

    const double sum_w_base_PT3 = h_base_PT3->Integral();
    const double sum_w_base_PT2 = h_base_PT2->Integral();
    const double sum_w_base_tot = sum_w_base_PT3 + sum_w_base_PT2;

    // ============================================================================
    // χ²/ndf landscape (single panel)
    // ============================================================================
    std::cout << "\nComputing flatness landscape (constant fit chi2/ndf)...\n";

    TH2D* h_chi2 = new TH2D(
        "h_chi2_global",
        "Matchqualnorm flatness landscape — single fit in [0.1, 1.1] GeV/c^{2};"
        "ep_richmatchqualitynorm #leq;em_richmatchqualitynorm #leq;#chi^{2}/ndf",
        n_cuts, cut_min - 0.5, cut_max + 0.5,
        n_cuts, cut_min - 0.5, cut_max + 0.5);
    h_chi2->SetDirectory(nullptr);

    BestPoint best;

    for (int ei = 0; ei < n_cuts; ++ei) {
        for (int mi = 0; mi < n_cuts; ++mi) {
            TH1D* h_PT3 = h_epem[ei][mi][TRG_PT3];
            TH1D* h_PT2 = h_epem[ei][mi][TRG_PT2];
            TH1D* r_ratio = makeRatio(h_PT2, h_PT3,
                Form("r_tmp_e%d_m%d", ei, mi));

            const FitResult fw = fitConstRange(r_ratio, kFitLo, kFitHi);
            h_chi2->SetBinContent(ei + 1, mi + 1, fw.chi2_ndf);

            if (fw.chi2_ndf < best.chi2_ndf && fw.ndf > 0) {
                best.ep_cut    = cut_min + ei;
                best.em_cut    = cut_min + mi;
                best.a         = fw.a;
                best.a_err     = fw.e;
                best.chi2_ndf  = fw.chi2_ndf;
                best.ndf       = fw.ndf;
                best.sum_w_PT3 = h_PT3->Integral();
                best.sum_w_PT2 = h_PT2->Integral();
            }
            delete r_ratio;
        }
    }

    // Baseline (no-cut) FOM for context.
    TH1D* r_base = makeRatio(h_base_PT2, h_base_PT3, "r_base_global");
    const FitResult f_base = fitConstRange(r_base, kFitLo, kFitHi);

    // ============================================================================
    // Summary printout
    // ============================================================================
    const double sum_w_best_tot = best.sum_w_PT3 + best.sum_w_PT2;
    const double retention_pct  = (sum_w_base_tot > 0.0)
        ? 100.0 * sum_w_best_tot / sum_w_base_tot : 0.0;
    const bool   at_boundary    =
        (best.ep_cut == cut_min || best.ep_cut == cut_max ||
         best.em_cut == cut_min || best.em_cut == cut_max);

    std::cout << "\n=== Single global FOM scan summary ===\n"
              << "  Fit range:           [" << kFitLo << ", " << kFitHi
              << "] GeV/c²\n"
              << "  Baseline (no cut):   a = " << f_base.a
              << " ± " << f_base.e
              << "   chi2/ndf = " << f_base.chi2_ndf
              << " / " << f_base.ndf << "\n"
              << "  Best (ep_cut, em_cut): (" << best.ep_cut
              << ", " << best.em_cut << ")\n"
              << "  Best a (level):      " << best.a << " ± " << best.a_err << "\n"
              << "  Best chi2/ndf:       " << best.chi2_ndf
              << " / " << best.ndf << "\n"
              << "  Sum weights @ best:  PT3 = " << best.sum_w_PT3
              << "   PT2 = " << best.sum_w_PT2 << "\n"
              << "  Sum weights @ base:  PT3 = " << sum_w_base_PT3
              << "   PT2 = " << sum_w_base_PT2 << "\n"
              << "  Sample retention:    " << retention_pct
              << " %  (best_total / baseline_total)\n";

    if (at_boundary) {
        std::cout << "  *** WARNING *** best lands at grid boundary "
                     "— scan range may be too narrow.\n";
    }
    if (retention_pct < 1.0) {
        std::cout << "  *** RED FLAG *** retention < 1 % at best — "
                     "statistical noise from tiny sample is likely "
                     "deflating chi2/ndf artificially.\n";
    }

    gSystem->mkdir("plots/output", true);

    // ============================================================================
    // OUTPUT 1 — Heatmap canvas (single panel)
    // ============================================================================
    TCanvas* c_heat = new TCanvas("c_matchqualnorm_global_flatness_sim_heatmap",
        "Matchqualnorm flatness landscape (sim) — global FOM", 800, 700);
    c_heat->cd();
    gPad->SetMargin(0.13, 0.16, 0.13, 0.11);
    gPad->SetLogz(true);
    h_chi2->Draw("COLZ");

    TMarker* mk = new TMarker(best.ep_cut, best.em_cut, 29);
    mk->SetMarkerSize(2.5);
    mk->SetMarkerColor(kRed);
    mk->Draw();

    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.034);
    tx.DrawLatex(0.16, 0.93,
        Form("best: (%d, %d)  #chi^{2}/ndf = %.2f / %d",
             best.ep_cut, best.em_cut, best.chi2_ndf, best.ndf));
    tx.DrawLatex(0.16, 0.89,
        Form("a = %.4f #pm %.4f   retention = %.2f %%",
             best.a, best.a_err, retention_pct));
    if (at_boundary) {
        tx.SetTextColor(kRed + 1);
        tx.DrawLatex(0.16, 0.85, "WARNING: best at grid boundary");
    }
    if (retention_pct < 1.0) {
        tx.SetTextColor(kRed + 1);
        tx.DrawLatex(0.16, 0.81, "RED FLAG: retention < 1 % (noise artifact)");
    }

    c_heat->SaveAs("plots/output/scan_matchqualnorm_global_flatness_sim_heatmap.pdf");
    c_heat->SaveAs("plots/output/scan_matchqualnorm_global_flatness_sim_heatmap.png");
    std::cout << "\nSaved: plots/output/scan_matchqualnorm_global_flatness_sim_heatmap.{pdf,png}\n";

    // ============================================================================
    // OUTPUT 2 — Best-case canvas (1×2)
    // ============================================================================
    {
        const int ei = best.ep_cut - cut_min;
        const int mi = best.em_cut - cut_min;

        TH1D* h_PT3 = static_cast<TH1D*>(h_epem[ei][mi][TRG_PT3]->Clone(
            "c_PT3_global"));
        TH1D* h_PT2 = static_cast<TH1D*>(h_epem[ei][mi][TRG_PT2]->Clone(
            "c_PT2_global"));
        h_PT3->SetDirectory(nullptr);
        h_PT2->SetDirectory(nullptr);
        TH1D* r_ratio = makeRatio(h_PT2, h_PT3, "c_ratio_global");

        // Fit to attach for visualization
        TF1* f_fit = new TF1("f_fit_global", "[0]", kFitLo, kFitHi);
        f_fit->SetLineColor(kBlue + 2);
        f_fit->SetLineWidth(3);
        r_ratio->Fit(f_fit, "RQ0");

        styleSpec(h_PT3,   kBlack,    21);     // filled square
        styleSpec(h_PT2,   kRed + 1,  20);     // filled circle
        styleSpec(r_ratio, kBlue + 1, 20);

        TCanvas* c = new TCanvas(
            "c_matchqualnorm_global_flatness_best",
            Form("Best — ep#leq%d em#leq%d", best.ep_cut, best.em_cut),
            1700, 600);
        c->Divide(2, 1, 0.001, 0.001);

        // ── Left: m_ee_sim spectrum (log Y) at best cuts ──
        c->cd(1);
        gPad->SetLogy(true);
        gPad->SetMargin(0.13, 0.04, 0.13, 0.10);
        const double y_max = std::max(h_PT3->GetMaximum(), h_PT2->GetMaximum());
        const double y_lo  = std::max(1e-2, y_max * 1e-6);
        h_PT3->SetTitle(Form(
            "m_{ee}^{sim} at best cuts ep#leq%d em#leq%d;"
            "M_{e^{+}e^{-}}^{sim} [GeV/c^{2}];weighted entries / 20 MeV",
            best.ep_cut, best.em_cut));
        h_PT3->GetYaxis()->SetRangeUser(y_lo, y_max * 5.0);
        h_PT3->Draw("E1");
        h_PT2->Draw("E1 SAME");

        TLegend* leg = new TLegend(0.55, 0.74, 0.95, 0.88);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.034);
        leg->AddEntry(h_PT3, Form("PT3  (#sum w = %.0f)", h_PT3->Integral()), "lpe");
        leg->AddEntry(h_PT2, Form("PT2  (#sum w = %.0f)", h_PT2->Integral()), "lpe");
        leg->Draw();

        // ── Right: ratio N_PT2/N_PT3 with constant fit drawn ──
        c->cd(2);
        gPad->SetMargin(0.13, 0.04, 0.13, 0.10);
        r_ratio->SetTitle(
            ";M_{e^{+}e^{-}}^{sim} [GeV/c^{2}];N_{PT2} / N_{PT3}");
        r_ratio->GetYaxis()->SetRangeUser(0.0, 3.0);
        r_ratio->Draw("E1");

        // Shade fit region.
        TBox* box = new TBox(kFitLo, 0.0, kFitHi, 3.0);
        box->SetFillColorAlpha(kGreen - 9, 0.15);
        box->SetLineColor(kGreen + 2); box->SetLineStyle(2);
        box->Draw("SAME");
        r_ratio->Draw("E1 SAME");
        f_fit->Draw("SAME");

        TLine* lref = new TLine(kXmin, 1.0, kXmax, 1.0);
        lref->SetLineStyle(3); lref->SetLineColor(kGray + 2); lref->Draw();

        TLatex tx2; tx2.SetNDC(); tx2.SetTextAlign(11); tx2.SetTextSize(0.034);
        tx2.DrawLatex(0.16, 0.90,
            Form("best ep#leq%d  em#leq%d", best.ep_cut, best.em_cut));
        tx2.SetTextColor(kBlue + 2);
        tx2.SetTextSize(0.030);
        tx2.DrawLatex(0.16, 0.86,
            Form("[%.2f-%.2f]: a = %.4f #pm %.4f   #chi^{2}/ndf = %.2f / %d",
                 kFitLo, kFitHi, best.a, best.a_err, best.chi2_ndf, best.ndf));
        tx2.SetTextColor(kGray + 3);
        tx2.DrawLatex(0.16, 0.82,
            Form("baseline: a = %.4f #pm %.4f   #chi^{2}/ndf = %.2f / %d",
                 f_base.a, f_base.e, f_base.chi2_ndf, f_base.ndf));
        tx2.SetTextColor(kBlack);
        tx2.DrawLatex(0.16, 0.78,
            Form("retention = %.2f %%", retention_pct));
        if (at_boundary) {
            tx2.SetTextColor(kRed + 1);
            tx2.DrawLatex(0.16, 0.74, "WARNING: best at grid boundary");
        }
        if (retention_pct < 1.0) {
            tx2.SetTextColor(kRed + 1);
            tx2.DrawLatex(0.16, 0.70, "RED FLAG: retention < 1 %");
        }

        const std::string base =
            "plots/output/scan_matchqualnorm_global_flatness_sim_best";
        c->SaveAs((base + ".pdf").c_str());
        c->SaveAs((base + ".png").c_str());
        std::cout << "Saved: " << base << ".{pdf,png}\n";
    }

    delete r_base;

    std::cout << "\nDone. 2 output files saved under plots/output/.\n";
}
