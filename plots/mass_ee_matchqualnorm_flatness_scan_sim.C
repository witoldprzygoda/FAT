// mass_ee_matchqualnorm_flatness_scan_sim.C — 2D scan of
// (ep_richmatchqualitynorm, em_richmatchqualitynorm) cuts on the SIM sample
// to find combinations that FLATTEN the PT2/PT3 trigger-ratio modulation in
// m_ee_sim space.
//
// Cut convention (DIFFERS from padnum scan):
//   ep_richmatchqualitynorm <= ep_cut    AND
//   em_richmatchqualitynorm <= em_cut
//   Lower cut value = TIGHTER; cut = 22 = LOOSEST (all events pass).
//   This means heatmap reads "left-right = tight → loose" instead of loose → tight.
//
// Differences vs the exp matchqualnorm scan:
//   • NO combinatorial background (sim has the MC purity gate baked in
//     during analysis ⇒ every event in `output_epem_sim.root` is already a
//     clean truth-matched e⁺e⁻ pair; no like-sign sibling files needed).
//   • The FOM is the χ²/ndf of a constant fit to N_PT2/N_PT3 vs m_ee_sim
//     inside the window — LOWER is BETTER (= less sinus modulation).
//   • Per-window best cut: minimizes window χ²/ndf.
//   • TOP marker: the window whose best cut gives the smallest χ²/ndf in
//     the above-π⁰ range [0.14, 1.4] GeV/c² (mirrors the exp convention).
//
// What is scanned (ASYMMETRIC ep / em):
//   Default grid: cut ∈ [0, 22] (23 values per leg → 529 combinations).
//   The branches ep_richmatchqualitynorm / em_richmatchqualitynorm carry
//   fractional float values, so a per-event Q maps to the smallest passing
//   integer cut via ceil(max(0, Q)).
//
// Mass windows:
//   A) π⁰         0.000 - 0.135 GeV/c²
//   B) η Dalitz   0.135 - 0.600 GeV/c²
//   C) high mass  0.600 - 1.400 GeV/c²
//   D) full       0.000 - 1.400 GeV/c²
//   TOP-range used for ★ selection: 0.140 - 1.400 GeV/c² (above π⁰)
//
// Outputs (5 files total — names distinct from exp scans by `_flatness_sim`):
//   plots/output/scan_matchqualnorm_flatness_sim_heatmap.{pdf,png}
//   plots/output/scan_matchqualnorm_flatness_sim_best_<window>.{pdf,png}  (×4)
//
// Each best-case canvas is 1×2:
//   Left:  m_ee_sim spectrum at best cuts — PT3 (black filled square) overlaid
//          with PT2 (red filled circle), log Y, target window highlighted.
//   Right: ratio N_PT2/N_PT3 with constant fit drawn (line in window range),
//          linear Y [0, 3], y=1 ref line, annotation listing
//          a ± err, χ²/ndf for the WINDOW and the TOP range.
//
// Usage:
//   root -l -b -q plots/mass_ee_matchqualnorm_flatness_scan_sim.C
//   root -l -b -q 'plots/mass_ee_matchqualnorm_flatness_scan_sim.C(0, 22)'
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
constexpr double kTopLo = 0.140;
constexpr double kTopHi = 1.400;

// -------- Trigger encoding --------
enum Trigger : int { TRG_PT2 = 0, TRG_PT3 = 1, NUM_TRG = 2 };
const std::array<std::string, NUM_TRG> kTrgName  = {"PT2", "PT3"};
const std::array<std::string, NUM_TRG> kTrgLabel = {"PT2 trigger",
                                                    "PT3 trigger"};

// -------- Mass windows --------
struct Window {
    std::string short_id;
    std::string label;
    double      lo;
    double      hi;
};
const std::vector<Window> kWindows = {
    {"A_pi0",      "A: #pi^{0} (0-0.135)",         0.000, 0.135},
    {"B_etaDal",   "B: #eta Dalitz (0.135-0.6)",   0.135, 0.600},
    {"C_highmass", "C: high mass (0.6-1.4)",       0.600, 1.400},
    {"D_full",     "D: full (0-1.4)",              0.000, 1.400},
};

// -------- Best-point bookkeeping --------
struct BestPoint {
    int    ep_cut       = 0;
    int    em_cut       = 0;
    double a            = 0.0;
    double a_err        = 0.0;
    double chi2_ndf     = std::numeric_limits<double>::infinity();
    int    ndf          = 0;
    double top_chi2_ndf = std::numeric_limits<double>::infinity();
    int    top_ndf      = 0;
    double top_a        = 0.0;
    double top_a_err    = 0.0;
    bool   is_top       = false;
};

// -------- Thousands-separator helper (kept consistent with exp scans) --------
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

void mass_ee_matchqualnorm_flatness_scan_sim(int cut_min = 0, int cut_max = 22) {
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
              << "Convention: qual <= cut_value passes  "
                 "(cut=0 tightest, cut=22 loosest).\n"
              << "FOM: chi2/ndf of constant fit to N_PT2/N_PT3 vs m_ee_sim "
                 "in each window (LOWER = FLATTER).\n";

    // ============================================================================
    // Allocate histograms: h_epem[ei][mi][trg] over 5 MeV grid for fine bins.
    // Use 70 × 20 MeV bins to match the per-channel macro and ease cross-check.
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

    // ============================================================================
    // Event loop — single pass through dilepton_nt_cor (truth m_ee_sim)
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
            // qual <= cut passes iff cut >= ceil(qual). Use ceil on the
            // non-negative quality value (clip negatives to 0).
            const int min_e = std::max(cut_min,
                static_cast<int>(std::ceil(std::max(0.0f, *v_ep))));
            const int min_m = std::max(cut_min,
                static_cast<int>(std::ceil(std::max(0.0f, *v_em))));
            if (min_e > cut_max || min_m > cut_max) continue;
            ++n_use;
            const double m = *mee;
            const double w = *v_w;
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

    // ============================================================================
    // χ²/ndf landscape per window
    // ============================================================================
    std::cout << "\nComputing flatness landscape (constant fit chi2/ndf)...\n";

    std::vector<std::vector<TH2D*>> h_chi2(
        kWindows.size(), std::vector<TH2D*>(1, nullptr));   // single-trg compatibility
    for (size_t w = 0; w < kWindows.size(); ++w) {
        h_chi2[w][0] = new TH2D(
            Form("h_chi2_%s", kWindows[w].short_id.c_str()),
            Form("Flatness #chi^{2}/ndf — %s;ep_richmatchqualitynorm #leq;em_richmatchqualitynorm #leq;#chi^{2}/ndf",
                 kWindows[w].label.c_str()),
            n_cuts, cut_min - 0.5, cut_max + 0.5,
            n_cuts, cut_min - 0.5, cut_max + 0.5);
        h_chi2[w][0]->SetDirectory(nullptr);
    }

    std::vector<BestPoint> best(kWindows.size());

    for (int ei = 0; ei < n_cuts; ++ei) {
        for (int mi = 0; mi < n_cuts; ++mi) {
            // Build ratio once per (ep, em) — used for window and TOP fits.
            TH1D* h_PT3 = h_epem[ei][mi][TRG_PT3];
            TH1D* h_PT2 = h_epem[ei][mi][TRG_PT2];
            TH1D* r_ratio = makeRatio(h_PT2, h_PT3,
                Form("r_tmp_e%d_m%d", ei, mi));

            // TOP-range fit (used to mark ★ TOP later).
            const FitResult f_top = fitConstRange(r_ratio, kTopLo, kTopHi);

            for (size_t w = 0; w < kWindows.size(); ++w) {
                const FitResult fw = fitConstRange(
                    r_ratio, kWindows[w].lo, kWindows[w].hi);
                h_chi2[w][0]->SetBinContent(ei + 1, mi + 1, fw.chi2_ndf);
                if (fw.chi2_ndf < best[w].chi2_ndf && fw.ndf > 0) {
                    best[w] = BestPoint{
                        cut_min + ei, cut_min + mi,
                        fw.a, fw.e, fw.chi2_ndf, fw.ndf,
                        f_top.chi2_ndf, f_top.ndf,
                        f_top.a, f_top.e,
                        false};
                }
            }
            delete r_ratio;
        }
    }

    // TOP = window whose best cut yields lowest top_chi2_ndf
    {
        double min_top = std::numeric_limits<double>::infinity();
        size_t tw = 0;
        for (size_t w = 0; w < kWindows.size(); ++w)
            if (best[w].top_chi2_ndf < min_top) {
                min_top = best[w].top_chi2_ndf;
                tw = w;
            }
        best[tw].is_top = true;
        std::cout << "TOP criterion (min χ²/ndf in [" << kTopLo << ", "
                  << kTopHi << "] GeV): window " << kWindows[tw].short_id
                  << " at (ep_cut=" << best[tw].ep_cut
                  << ", em_cut="     << best[tw].em_cut
                  << ")  top_chi2/ndf = " << min_top << "\n";
    }

    // ============================================================================
    // Summary table
    // ============================================================================
    auto fmtFloat = [](double v, int prec = 1) {
        std::ostringstream os;
        os << std::fixed << std::setprecision(prec) << v;
        return os.str();
    };
    auto pad = [](const std::string& s, int w) {
        if ((int)s.size() >= w) return s;
        return std::string(w - s.size(), ' ') + s;
    };

    std::cout << "\n=== Best (ep_cut, em_cut) per window — minimize window χ²/ndf ===\n"
              << "  (★ = TOP = lowest χ²/ndf in [" << kTopLo << ", "
              << kTopHi << "] GeV across all 4 cases)\n\n";
    std::cout << "    " << std::left << std::setw(30) << "window"
              << "ep  em  |  "
              << pad("a (level)", 10) << "  "
              << pad("a_err",      8) << "  "
              << pad("win_chi2/ndf",  14) << "  |  "
              << pad("top_a",      10) << "  "
              << pad("top_chi2/ndf",  14) << "\n";
    std::cout << "    " << std::string(110, '-') << "\n";
    for (size_t w = 0; w < kWindows.size(); ++w) {
        const auto& b = best[w];
        std::cout
            << (b.is_top ? "  ★ " : "    ")
            << std::left << std::setw(30) << kWindows[w].label
            << pad(std::to_string(b.ep_cut), 2) << "  "
            << pad(std::to_string(b.em_cut), 2) << "  |  "
            << pad(fmtFloat(b.a,     4), 10) << "  "
            << pad(fmtFloat(b.a_err, 4),  8) << "  "
            << pad(fmtFloat(b.chi2_ndf, 2) + " / "
                       + std::to_string(b.ndf), 14) << "  |  "
            << pad(fmtFloat(b.top_a, 4), 10) << "  "
            << pad(fmtFloat(b.top_chi2_ndf, 2) + " / "
                       + std::to_string(b.top_ndf), 14) << "\n";
    }
    std::cout << "\n";

    gSystem->mkdir("plots/output", true);

    // ============================================================================
    // OUTPUT 1 — Heatmap canvas (1 row × 4 cols)
    // ============================================================================
    TCanvas* c_heat = new TCanvas("c_matchqualnorm_flatness_sim_heatmap",
        "Matchqualnorm flatness landscape (sim) — 4 windows", 1900, 500);
    c_heat->Divide(kWindows.size(), 1, 0.001, 0.001);
    for (size_t w = 0; w < kWindows.size(); ++w) {
        c_heat->cd(w + 1);
        gPad->SetMargin(0.13, 0.16, 0.13, 0.11);
        // Log color scale to handle wide dynamic range of chi2/ndf.
        gPad->SetLogz(true);
        h_chi2[w][0]->Draw("COLZ");
        const auto& b = best[w];
        TMarker* mk = new TMarker(b.ep_cut, b.em_cut, 29);
        mk->SetMarkerSize(b.is_top ? 3.2 : 2.0);
        mk->SetMarkerColor(b.is_top ? kBlack : kRed);
        mk->Draw();
        if (b.is_top) {
            TMarker* halo = new TMarker(b.ep_cut, b.em_cut, 24);
            halo->SetMarkerSize(5.0);
            halo->SetMarkerColor(kYellow + 2); halo->Draw();
            TMarker* halo2 = new TMarker(b.ep_cut, b.em_cut, 24);
            halo2->SetMarkerSize(4.0);
            halo2->SetMarkerColor(kOrange + 7); halo2->Draw();
            mk->Draw();
        }
        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.040);
        tx.DrawLatex(0.16, 0.93,
            Form("best: (%d, %d)  win #chi^{2}/ndf=%.2f",
                 b.ep_cut, b.em_cut, b.chi2_ndf));
        tx.SetTextSize(0.034);
        tx.SetTextColor(kGray + 3);
        tx.DrawLatex(0.16, 0.88,
            Form("top #chi^{2}/ndf=%.2f", b.top_chi2_ndf));
        if (b.is_top) {
            tx.SetTextColor(kRed + 1);
            tx.SetTextSize(0.044);
            tx.DrawLatex(0.55, 0.93, "#bigstar TOP");
        }
    }
    c_heat->SaveAs("plots/output/scan_matchqualnorm_flatness_sim_heatmap.pdf");
    c_heat->SaveAs("plots/output/scan_matchqualnorm_flatness_sim_heatmap.png");
    std::cout << "Saved: plots/output/scan_matchqualnorm_flatness_sim_heatmap.{pdf,png}\n";

    // ============================================================================
    // OUTPUT 2-5 — Per-window best-case canvases (1 row × 2 cols)
    // ============================================================================
    for (size_t w = 0; w < kWindows.size(); ++w) {
        const auto& b = best[w];
        const int ei = b.ep_cut - cut_min;
        const int mi = b.em_cut - cut_min;

        TH1D* h_PT3 = static_cast<TH1D*>(h_epem[ei][mi][TRG_PT3]->Clone(
            Form("c_PT3_%s", kWindows[w].short_id.c_str())));
        TH1D* h_PT2 = static_cast<TH1D*>(h_epem[ei][mi][TRG_PT2]->Clone(
            Form("c_PT2_%s", kWindows[w].short_id.c_str())));
        h_PT3->SetDirectory(nullptr);
        h_PT2->SetDirectory(nullptr);
        TH1D* r_ratio = makeRatio(h_PT2, h_PT3,
            Form("c_ratio_%s", kWindows[w].short_id.c_str()));

        // Fits to attach for visualization
        TF1* f_win = new TF1(
            Form("f_win_%s", kWindows[w].short_id.c_str()),
            "[0]", kWindows[w].lo, kWindows[w].hi);
        f_win->SetLineColor(kBlue + 2); f_win->SetLineWidth(3);
        r_ratio->Fit(f_win, "RQ0");
        TF1* f_top = new TF1(
            Form("f_top_%s", kWindows[w].short_id.c_str()),
            "[0]", kTopLo, kTopHi);
        f_top->SetLineColor(kGray + 2); f_top->SetLineWidth(2); f_top->SetLineStyle(2);
        r_ratio->Fit(f_top, "RQ0+");

        styleSpec(h_PT3, kBlack,    21);     // filled square
        styleSpec(h_PT2, kRed + 1,  20);     // filled circle
        styleSpec(r_ratio, kBlue + 1, 20);

        TCanvas* c = new TCanvas(
            Form("c_matchqualnorm_flatness_best_%s", kWindows[w].short_id.c_str()),
            Form("Best %s — ep_qual<=%d  em_qual<=%d",
                 kWindows[w].label.c_str(), b.ep_cut, b.em_cut),
            1700, 600);
        c->Divide(2, 1, 0.001, 0.001);

        // ── Left: m_ee_sim spectrum (log Y) at best cuts ──
        c->cd(1);
        gPad->SetLogy(true);
        gPad->SetMargin(0.13, 0.04, 0.13, 0.10);
        const double y_max = std::max(h_PT3->GetMaximum(), h_PT2->GetMaximum());
        h_PT3->SetTitle(Form(
            "m_{ee}^{sim} at best cuts ep#leq%d em#leq%d;M_{e^{+}e^{-}}^{sim} [GeV/c^{2}];weighted entries / 20 MeV",
            b.ep_cut, b.em_cut));
        h_PT3->GetYaxis()->SetRangeUser(std::max(1e-2, y_max * 1e-6),
                                        y_max * 5.0);
        h_PT3->Draw("E1");
        h_PT2->Draw("E1 SAME");
        TBox* box = new TBox(kWindows[w].lo, std::max(1e-2, y_max * 1e-6),
                             kWindows[w].hi, y_max * 5.0);
        box->SetFillColorAlpha(kGreen - 9, 0.20);
        box->SetLineColor(kGreen + 2); box->SetLineStyle(2);
        box->Draw("SAME");
        h_PT3->Draw("E1 SAME");
        h_PT2->Draw("E1 SAME");

        TLegend* leg = new TLegend(0.55, 0.72, 0.95, 0.88);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.034);
        leg->AddEntry(h_PT3, Form("PT3  (#sum w = %.0f)", h_PT3->Integral()), "lpe");
        leg->AddEntry(h_PT2, Form("PT2  (#sum w = %.0f)", h_PT2->Integral()), "lpe");
        leg->AddEntry(box,   Form("window [%.2f, %.2f]",
                                  kWindows[w].lo, kWindows[w].hi), "f");
        leg->Draw();

        // ── Right: ratio N_PT2/N_PT3 with both fits drawn ──
        c->cd(2);
        gPad->SetMargin(0.13, 0.04, 0.13, 0.10);
        r_ratio->SetTitle(
            ";M_{e^{+}e^{-}}^{sim} [GeV/c^{2}];N_{PT2} / N_{PT3}");
        r_ratio->GetYaxis()->SetRangeUser(0.0, 3.0);
        r_ratio->Draw("E1");

        // Shade window region.
        TBox* box2 = new TBox(kWindows[w].lo, 0.0,
                              kWindows[w].hi, 3.0);
        box2->SetFillColorAlpha(kGreen - 9, 0.15);
        box2->SetLineColor(kGreen + 2); box2->SetLineStyle(2);
        box2->Draw("SAME");
        r_ratio->Draw("E1 SAME");

        f_win->Draw("SAME");
        f_top->Draw("SAME");

        TLine* lref = new TLine(kXmin, 1.0, kXmax, 1.0);
        lref->SetLineStyle(3); lref->SetLineColor(kGray + 2); lref->Draw();

        TLatex tx; tx.SetNDC(); tx.SetTextAlign(11); tx.SetTextSize(0.034);
        tx.DrawLatex(0.16, 0.90,
            Form("best ep#leq%d  em#leq%d", b.ep_cut, b.em_cut));
        tx.SetTextColor(kBlue + 2);
        tx.SetTextSize(0.030);
        tx.DrawLatex(0.16, 0.85,
            Form("window: a = %.4f #pm %.4f   #chi^{2}/ndf = %.2f / %d",
                 b.a, b.a_err, b.chi2_ndf, b.ndf));
        tx.SetTextColor(kGray + 3);
        tx.DrawLatex(0.16, 0.81,
            Form("[%.2f-%.2f]: a = %.4f #pm %.4f   #chi^{2}/ndf = %.2f / %d",
                 kTopLo, kTopHi, b.top_a, b.top_a_err,
                 b.top_chi2_ndf, b.top_ndf));
        if (b.is_top) {
            tx.SetTextColor(kRed + 1);
            tx.SetTextSize(0.044);
            tx.DrawLatex(0.18, 0.22, "#bigstar TOP");
        }

        const std::string base = Form(
            "plots/output/scan_matchqualnorm_flatness_sim_best_%s",
            kWindows[w].short_id.c_str());
        c->SaveAs((base + ".pdf").c_str());
        c->SaveAs((base + ".png").c_str());
        std::cout << "  saved: " << base << ".{pdf,png}  "
                  << "(ep_cut=" << b.ep_cut
                  << ", em_cut=" << b.em_cut
                  << ", win_chi2/ndf=" << b.chi2_ndf
                  << ", top_chi2/ndf=" << b.top_chi2_ndf << ")\n";
    }

    std::cout << "\nDone. 5 output files saved under plots/output/.\n";
}
