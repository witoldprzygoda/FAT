// mass_ee_trigger_ratio_truth_vs_rec_sim.C — H3 diagnostic for the SMASH
// simulation trigger correction factor N_PT2 / N_PT3.
//
// HYPOTHESIS H3 (reconstruction/digitisation smearing causes the sinus shape):
//   In experimental data the PT2/PT3 ratio (with the factor 63 absorbed) sits
//   around ~2.5 and is reasonably flat in m_ee. In the SMASH lepton sim the
//   same ratio (without the factor 63 — PT2 is not downscaled in sim) is much
//   lower (~1.6) AND shows a sinus-like modulation with a minimum near
//   m_ee ≈ 0.5 GeV/c² and a maximum near 1.2 GeV/c². One candidate cause is
//   that detector/digitisation smearing turns sharp threshold structures (in
//   p, OA, etc.) into smooth m_ee-dependent dips and bumps after migration.
//
//   This macro tests H3 directly by computing the SAME ratio in THREE
//   alternative m_ee representations and comparing the χ²/ndf of a constant
//   fit on [0.1, 1.4] GeV/c²:
//
//     REC (m_ee from dilepton_nt = raw reconstructed kinematics)
//     COR (m_ee from dilepton_nt_cor = with electron energy-loss correction)
//     SIM (m_ee_sim from dilepton_nt_cor = MC truth, no detector smearing)
//
//   Interpretation:
//     * If the sinus mostly disappears in SIM (lowest χ²/ndf, flattest) then
//       H3 is supported — the modulation is reco/smearing-induced.
//     * If all three look identically wavy then H3 is disfavoured — the
//       shape is intrinsic to the SMASH process mix (H1) or to a true
//       kinematic feature visible already at truth level (H2/H4).
//     * COR sitting between REC and SIM (or close to REC) tells whether
//       the energy-loss correction alone undoes part of the smearing.
//
// CANVAS LAYOUT (1700 × 900, 2 rows × 3 columns):
//
//   Row 1 — m_ee spectra (log Y), PT3 = black filled circle,
//                                 PT2 = red filled square:
//     [1] m_ee (REC)        [2] m_ee (COR)        [3] m_ee (SIM truth)
//
//   Row 2 — Trigger correction factor N_PT2 / N_PT3 (linear Y, [0, 3]),
//           with a dashed gray reference line at y = 1, kBlue+1 markers,
//           and a constant fit (kBlue+2) on [0.1, 1.4] GeV/c²:
//     [4] ratio (REC)       [5] ratio (COR)       [6] ratio (SIM truth)
//
// All draws use the per-event SMASH luminosity weight `sim_genweight`. PT2
// and PT3 selection use the auto-stamped event flags pt2, pt3.
//
// Arguments (all optional, sentinel −1 ⇒ use macro-level defaults [0.1, 1.4]):
//   1) fit_xmin   — fit range lower bound [GeV/c²], default 0.1
//   2) fit_xmax   — fit range upper bound [GeV/c²], default 1.4
//
// Usage:
//   root -l -b -q plots/mass_ee_trigger_ratio_truth_vs_rec_sim.C
//   root -l -b -q 'plots/mass_ee_trigger_ratio_truth_vs_rec_sim.C(0.1, 1.0)'
//
// @author Witold Przygoda (witold.przygoda@uj.edu.pl)
// @date 2026

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

namespace {

// -------- Input --------
constexpr const char* kInputFile  = "output_epem_sim.root";
constexpr const char* kNtRec      = "dilepton_nt";       // m_ee = m_ee_rec
constexpr const char* kNtCor      = "dilepton_nt_cor";   // m_ee = m_ee_cor, has m_ee_sim

// -------- Binning (matches the rest of the sim macros) --------
constexpr int    kNb    = 70;
constexpr double kXmin  = 0.0;
constexpr double kXmax  = 1.4;     // 70 bins × 20 MeV

// -------- Constant-fit default range --------
constexpr double kFitLo = 0.1;
constexpr double kFitHi = 1.4;

// -------- Colors / markers --------
constexpr Color_t kColPT3      = kBlack;
constexpr Color_t kColPT2      = kRed + 1;
constexpr Style_t kMkPT3       = 20;        // filled circle
constexpr Style_t kMkPT2       = 21;        // filled square
constexpr Color_t kColRatio    = kBlue + 1;
constexpr Color_t kColRatioFit = kBlue + 2;
constexpr Style_t kMkRatio     = 20;        // filled circle

// -------- Hardcoded sim weight expression (sim_genweight is always > 0 here) --------
constexpr const char* kWeight  = "sim_genweight";

// Draw a histogram from a tree given the m_ee expression to project and the
// (cut * weight) string. The histogram is detached from the file directory.
TH1D* drawFromNt(TTree* t, const std::string& mee_expr,
                 const std::string& cut_w, const std::string& name) {
    TH1D* h = new TH1D(name.c_str(), "", kNb, kXmin, kXmax);
    h->Sumw2();
    t->Draw((mee_expr + ">>" + name).c_str(), cut_w.c_str(), "goff");
    h->SetDirectory(nullptr);
    return h;
}

// Per-bin ratio with relative-error propagation. Bins with den=0 stay empty.
TH1D* makeRatio(TH1D* num, TH1D* den, const std::string& name) {
    TH1D* r = static_cast<TH1D*>(num->Clone(name.c_str()));
    r->SetDirectory(nullptr);
    r->Reset();
    const int nb = num->GetNbinsX();
    for (int b = 1; b <= nb; ++b) {
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

void styleSpectrum(TH1D* h, Color_t color, Style_t marker) {
    h->SetMarkerStyle(marker);
    h->SetMarkerSize(0.9);
    h->SetMarkerColor(color);
    h->SetLineColor(color);
    h->GetXaxis()->SetTitleSize(0.050);
    h->GetYaxis()->SetTitleSize(0.050);
    h->GetXaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetLabelSize(0.045);
    h->GetXaxis()->SetTitleOffset(1.05);
    h->GetYaxis()->SetTitleOffset(1.25);
}

struct FitResult { double a; double e; double chi2_ndf; int ndf; };

// Weighted constant fit (option "RQ"). Color controls the drawn line color.
FitResult fitConst(TH1D* h, double fit_xmin, double fit_xmax, Color_t color) {
    const std::string fname = std::string(h->GetName()) + "_fconst";
    TF1* f = new TF1(fname.c_str(), "[0]", fit_xmin, fit_xmax);
    f->SetParName(0, "a");
    f->SetLineColor(color);
    f->SetLineWidth(2);
    f->SetLineStyle(1);
    h->Fit(f, "RQ");
    FitResult r;
    r.a        = f->GetParameter(0);
    r.e        = f->GetParError(0);
    r.ndf      = f->GetNDF();
    r.chi2_ndf = (r.ndf > 0) ? f->GetChisquare() / r.ndf : 0.0;
    return r;
}

// For log-Y range; ignore non-finite or non-positive bins.
double safePosMin(TH1D* h) {
    double mn = 0.0;
    const int nb = h->GetNbinsX();
    for (int b = 1; b <= nb; ++b) {
        const double v = h->GetBinContent(b);
        if (std::isfinite(v) && v > 0.0 && (mn == 0.0 || v < mn)) mn = v;
    }
    return mn;
}

double safeMax(TH1D* h) {
    double mx = 0.0;
    const int nb = h->GetNbinsX();
    for (int b = 1; b <= nb; ++b) {
        const double v = h->GetBinContent(b) + h->GetBinError(b);
        if (std::isfinite(v) && v > mx) mx = v;
    }
    return mx;
}

}  // anonymous namespace

void mass_ee_trigger_ratio_truth_vs_rec_sim(double fit_xmin = -1.0,
                                            double fit_xmax = -1.0) {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetErrorX(0.0);    // no horizontal error bars

    // Resolve fit range — sentinel −1 ⇒ default [kFitLo, kFitHi].
    const double fxlo = (fit_xmin >= 0.0 && fit_xmin < kXmax) ? fit_xmin : kFitLo;
    const double fxhi = (fit_xmax >  0.0 && fit_xmax <= kXmax) ? fit_xmax : kFitHi;
    std::cout << "Fit range: [" << fxlo << ", " << fxhi << "] GeV/c^2\n";

    // -------- Open input file and grab both trees --------
    TFile* f = TFile::Open(kInputFile, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Cannot open " << kInputFile << "\n";
        return;
    }
    TTree* nt_rec = dynamic_cast<TTree*>(f->Get(kNtRec));
    TTree* nt_cor = dynamic_cast<TTree*>(f->Get(kNtCor));
    if (!nt_rec) { std::cerr << kNtRec << " missing\n"; return; }
    if (!nt_cor) { std::cerr << kNtCor << " missing\n"; return; }

    // -------- Compose weight expressions --------
    // PT2 and PT3 are auto-stamped 0/1 flags; sim_genweight carries the
    // SMASH luminosity weight. No factor 63 in sim (PT2 not downscaled).
    const std::string w_pt3 = std::string("(pt3==1) * ") + kWeight;
    const std::string w_pt2 = std::string("(pt2==1) * ") + kWeight;
    std::cout << "Weight expressions:\n"
              << "  PT3 = " << w_pt3 << "\n"
              << "  PT2 = " << w_pt2 << "\n";

    // -------- Build all six base spectra --------
    // REC: m_ee from dilepton_nt (raw reconstructed)
    TH1D* h_pt3_rec = drawFromNt(nt_rec, "m_ee",     w_pt3, "h_PT3_REC");
    TH1D* h_pt2_rec = drawFromNt(nt_rec, "m_ee",     w_pt2, "h_PT2_REC");
    // COR: m_ee from dilepton_nt_cor (energy-loss corrected)
    TH1D* h_pt3_cor = drawFromNt(nt_cor, "m_ee",     w_pt3, "h_PT3_COR");
    TH1D* h_pt2_cor = drawFromNt(nt_cor, "m_ee",     w_pt2, "h_PT2_COR");
    // SIM truth: m_ee_sim from dilepton_nt_cor (no detector smearing)
    TH1D* h_pt3_sim = drawFromNt(nt_cor, "m_ee_sim", w_pt3, "h_PT3_SIM");
    TH1D* h_pt2_sim = drawFromNt(nt_cor, "m_ee_sim", w_pt2, "h_PT2_SIM");

    // -------- Build the three ratios N_PT2 / N_PT3 --------
    TH1D* r_rec = makeRatio(h_pt2_rec, h_pt3_rec, "r_PT2_PT3_REC");
    TH1D* r_cor = makeRatio(h_pt2_cor, h_pt3_cor, "r_PT2_PT3_COR");
    TH1D* r_sim = makeRatio(h_pt2_sim, h_pt3_sim, "r_PT2_PT3_SIM");

    // -------- Styling --------
    styleSpectrum(h_pt3_rec, kColPT3, kMkPT3);
    styleSpectrum(h_pt2_rec, kColPT2, kMkPT2);
    styleSpectrum(h_pt3_cor, kColPT3, kMkPT3);
    styleSpectrum(h_pt2_cor, kColPT2, kMkPT2);
    styleSpectrum(h_pt3_sim, kColPT3, kMkPT3);
    styleSpectrum(h_pt2_sim, kColPT2, kMkPT2);
    styleSpectrum(r_rec, kColRatio, kMkRatio);
    styleSpectrum(r_cor, kColRatio, kMkRatio);
    styleSpectrum(r_sim, kColRatio, kMkRatio);

    // -------- Canvas: 2 rows × 3 cols --------
    TCanvas* c = new TCanvas("c_mass_ee_trigger_ratio_truth_vs_rec_sim",
                             "H3 diagnostic: PT2/PT3 vs m_ee — REC/COR/SIM",
                             1700, 900);
    c->Divide(3, 2, 0.001, 0.001);

    // Common pad margins.
    const double pad_lm = 0.14, pad_rm = 0.04, pad_bm = 0.13, pad_tm = 0.10;

    // ---- Row 1 — Spectra (log Y), one helper to keep it DRY -----------------
    auto drawSpectrumPad =
        [&](int idx, TH1D* h_pt3, TH1D* h_pt2, const char* tag) {
            c->cd(idx);
            gPad->SetMargin(pad_lm, pad_rm, pad_bm, pad_tm);
            gPad->SetLogy(1);

            // Set Y range from both histograms; pad ×3 above max, /3 below min.
            const double y_max = std::max(safeMax(h_pt3), safeMax(h_pt2));
            const double y_pos = std::min(
                std::max(safePosMin(h_pt3), 0.0),
                std::max(safePosMin(h_pt2), 0.0));
            const double y_min_eff = (y_pos > 0.0) ? y_pos / 3.0 : 1e-3;
            const double y_max_eff = (y_max > 0.0) ? y_max * 3.0 : 1.0;

            h_pt3->SetTitle(
                (std::string("m_{ee} (") + tag + ");"
                 "M_{e^{+}e^{-}} [GeV/c^{2}];"
                 "weighted counts (sim_genweight)").c_str());
            h_pt3->GetYaxis()->SetRangeUser(y_min_eff, y_max_eff);

            h_pt3->Draw("E1");
            h_pt2->Draw("E1 SAME");

            // Yields for legend annotation.
            const double n_pt3 = h_pt3->Integral();
            const double n_pt2 = h_pt2->Integral();

            TLegend* leg = new TLegend(0.50, 0.74, 0.95, 0.88);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.035);
            leg->AddEntry(h_pt3,
                Form("PT3   N = %.3g", n_pt3), "lpe");
            leg->AddEntry(h_pt2,
                Form("PT2   N = %.3g", n_pt2), "lpe");
            leg->Draw();
        };

    drawSpectrumPad(1, h_pt3_rec, h_pt2_rec, "REC");
    drawSpectrumPad(2, h_pt3_cor, h_pt2_cor, "COR");
    drawSpectrumPad(3, h_pt3_sim, h_pt2_sim, "SIM truth");

    // ---- Row 2 — Ratio panels with constant fit -----------------------------
    auto drawRatioPad =
        [&](int idx, TH1D* r, const char* tag) -> FitResult {
            c->cd(idx);
            gPad->SetMargin(pad_lm, pad_rm, pad_bm, pad_tm);
            gPad->SetLogy(0);

            r->SetTitle(
                (std::string("Trigger correction factor — N_{PT2}/N_{PT3} (") +
                 tag + ");"
                 "M_{e^{+}e^{-}} [GeV/c^{2}];"
                 "N_{PT2} / N_{PT3}").c_str());
            r->GetYaxis()->SetRangeUser(0.0, 3.0);

            r->Draw("E1");

            // Reference line at y = 1 (dashed gray).
            TLine* lref = new TLine(kXmin, 1.0, kXmax, 1.0);
            lref->SetLineStyle(2);
            lref->SetLineColor(kGray + 2);
            lref->SetLineWidth(2);
            lref->Draw();

            const FitResult fit = fitConst(r, fxlo, fxhi, kColRatioFit);

            // Annotate fit values in the upper-left of the pad.
            TLatex tex;
            tex.SetNDC();
            tex.SetTextAlign(11);
            tex.SetTextSize(0.050);
            tex.SetTextColor(kColRatioFit);
            tex.DrawLatex(0.18, 0.85,
                Form("a = %.4f #pm %.4f", fit.a, fit.e));
            tex.SetTextSize(0.040);
            tex.SetTextColor(kGray + 3);
            tex.DrawLatex(0.18, 0.79,
                Form("#chi^{2}/ndf = %.1f / %d",
                     fit.chi2_ndf * fit.ndf, fit.ndf));
            return fit;
        };

    std::cout << "\nConstant fits y = a on N_PT2/N_PT3:\n";
    const FitResult fit_rec = drawRatioPad(4, r_rec, "REC");
    const FitResult fit_cor = drawRatioPad(5, r_cor, "COR");
    const FitResult fit_sim = drawRatioPad(6, r_sim, "SIM truth");

    std::cout << "  REC: a = " << fit_rec.a << " ± " << fit_rec.e
              << "   chi2/ndf = " << (fit_rec.chi2_ndf * fit_rec.ndf)
              << " / " << fit_rec.ndf
              << "   (chi2/ndf = " << fit_rec.chi2_ndf << ")\n"
              << "  COR: a = " << fit_cor.a << " ± " << fit_cor.e
              << "   chi2/ndf = " << (fit_cor.chi2_ndf * fit_cor.ndf)
              << " / " << fit_cor.ndf
              << "   (chi2/ndf = " << fit_cor.chi2_ndf << ")\n"
              << "  SIM: a = " << fit_sim.a << " ± " << fit_sim.e
              << "   chi2/ndf = " << (fit_sim.chi2_ndf * fit_sim.ndf)
              << " / " << fit_sim.ndf
              << "   (chi2/ndf = " << fit_sim.chi2_ndf << ")\n";

    // Tag the flattest variant — primary H3 indicator.
    const char* flattest = "REC";
    double best = fit_rec.chi2_ndf;
    if (fit_cor.chi2_ndf < best) { best = fit_cor.chi2_ndf; flattest = "COR"; }
    if (fit_sim.chi2_ndf < best) { best = fit_sim.chi2_ndf; flattest = "SIM"; }
    std::cout << "\nFlattest (lowest χ²/ndf) variant: " << flattest
              << "   (χ²/ndf = " << best << ")\n"
              << "  → if SIM ⇒ H3 (reco smearing) supported;\n"
              << "    if REC≈COR≈SIM all bad ⇒ H3 disfavoured.\n";

    // -------- Save --------
    gSystem->mkdir("plots/output", true);
    const std::string base =
        "plots/output/mass_ee_trigger_ratio_truth_vs_rec_sim";
    c->SaveAs((base + ".pdf").c_str());
    c->SaveAs((base + ".png").c_str());
    std::cout << "\nSaved: " << base << ".{pdf,png}\n";

    f->Close();
}
