// mass_ee_pt_ratio_purity_cascade_sim_3piece_v5.C
// =========================================================================
// V5 of the piecewise fit purity cascade study (sim) — UPDATED.
//
// Changes vs previous v5:
//
//   1) LOW piece [0.00, 0.15]: NO FIT. Data points are still drawn, but
//      no TF1 is attached and no annotation line is shown for LOW.
//      Binning in LOW is doubled in density:
//        [0.00, 0.15] : 60 bins of 2.5 MeV  (was 30 bins of 5 MeV)
//
//   2) Variable-width binning (220 bins total):
//        [0.00, 0.15] : 60 bins of 2.5 MeV         (2× finer than before)
//        [0.15, 0.85] : 140 bins of 5 MeV          (unchanged density)
//        [0.85, 1.10] :  12 bins of (0.25/12) GeV  (unchanged)
//        [1.10, 1.40] :   8 bins of (0.30/8) GeV   (unchanged)
//
//   3) HIGH quad fit on [0.70, 1.40] is FREE: y = a + b*x + c*x^2.
//      Exclusion windows applied (zero out bins in masked clone):
//        [0.73, 0.83]  ω band
//        [0.96, 1.08]  φ band
//        [1.30, 1.40]  high-mass tail (NEW)
//
//   4) MID quad fit on [0.15, 0.70] is now STITCHED to HIGH at x=0.70:
//        y(x) = HIGH(0.70) + b_mid*(x - 0.70) + c_mid*(x - 0.70)^2
//      HIGH(0.70) is a numeric constant; b_mid and c_mid are free.
//      Two free parameters now (slope at 0.70 and curvature).
//
// All other behaviour (5 cumulative cut stages, 2x2 layout, color coding,
// overlay canvases, width-normalized spectra) is unchanged.
//
// Usage:
//   root -l -b -q plots/mass_ee_pt_ratio_purity_cascade_sim_3piece_v5.C
//
// @author Witold Przygoda (witold.przygoda@uj.edu.pl)
// @date 2026

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace {

constexpr const char* kInputFile = "output_epem_sim.root";
constexpr const char* kNtRec     = "dilepton_nt";
constexpr const char* kNtTruth   = "dilepton_nt_cor";

constexpr int    kNb    = 220;   // 60 + 140 + 12 + 8 variable-width bins
constexpr double kXmin  = 0.0;
constexpr double kXmax  = 1.4;

// Piecewise fit ranges (GeV/c²)
// LOW piece: NO fit, [0.00, 0.15] just shows data.
constexpr double kMdLo  = 0.15, kMdHi = 0.70;   // MID   piece — quad stitched at 0.70
constexpr double kHiLo  = 0.70, kHiHi = 1.40;   // HIGH  piece — quad FREE

// Exclusion windows applied to HIGH fit (ω, φ, and high-mass regions)
constexpr double kExcl1Lo = 0.73, kExcl1Hi = 0.83;   // ω peak
constexpr double kExcl2Lo = 0.96, kExcl2Hi = 1.08;   // φ peak
constexpr double kExcl3Lo = 1.30, kExcl3Hi = 1.40;   // high-mass tail

constexpr Color_t kColPT3   = kBlack;
constexpr Color_t kColPT2   = kRed + 1;
constexpr Color_t kColRatio = kBlue + 1;

// Distinct color per fitted piece (LOW no longer fitted, color removed)
const Color_t kColMid  = kAzure  + 1;
const Color_t kColHigh = kOrange + 8;

struct Scenario {
    std::string id;
    std::string label;
    std::string cut;
};

// 5 CUMULATIVE cut stages — each shown in both REC and TRUTH variants
const std::vector<Scenario> kScenarios = {
    {"step1_raw",     "1) RAW — no cuts",
                      "1"},
    {"step2_vertex",  "2) +eVertReco_z > -500",
                      "(eVertReco_z>-500)"},
    {"step3_isBest",  "3) +isBest == 1",
                      "(eVertReco_z>-500)*(isBest==1)"},
    {"step4_simID",   "4) +ep_sim_id==2 && em_sim_id==3",
                      "(eVertReco_z>-500)*(isBest==1)*(ep_sim_id==2)*(em_sim_id==3)"},
    {"step5_samevtx", "5) +epem_same_vertex==1 (full old purity gate)",
                      "(eVertReco_z>-500)*(isBest==1)*(ep_sim_id==2)*(em_sim_id==3)*(epem_same_vertex==1)"},
};

const std::array<Color_t, 5> kOverlayColors = {
    kRed + 1, kOrange + 7, kGreen + 2, kAzure + 1, kViolet + 2
};

// Build the variable-width bin edge array for updated v5 binning.
//   [0.00, 0.15] :  60 bins of 2.5 MeV
//   [0.15, 0.85] : 140 bins of 5 MeV
//   [0.85, 1.10] :  12 bins of (0.25/12) GeV
//   [1.10, 1.40] :   8 bins of (0.30/8) GeV
const double* getEdges() {
    static double edges[kNb + 1];
    static bool   built = false;
    if (!built) {
        for (int i = 0; i <= 60;  ++i) edges[i]       = i * 0.0025;
        for (int i = 1; i <= 140; ++i) edges[60 + i]  = 0.15 + i * 0.005;
        for (int i = 1; i <= 12;  ++i) edges[200 + i] = 0.85 + i * (0.25 / 12.0);
        for (int i = 1; i <= 8;   ++i) edges[212 + i] = 1.10 + i * (0.30 / 8.0);
        built = true;
    }
    return edges;
}

TH1D* drawWithWeight(TTree* t, const std::string& mass_var,
                     const std::string& weight, const std::string& name) {
    TH1D* h = new TH1D(name.c_str(), "", kNb, getEdges());
    h->Sumw2();
    t->Draw((mass_var + ">>" + name).c_str(),
            weight.c_str(), "goff");
    h->SetDirectory(nullptr);
    return h;
}

// Divide each bin content and error by its bin width (for spectrum display).
void normalizeByBinWidth(TH1D* h) {
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
        const double w = h->GetBinWidth(b);
        if (w > 0.0) {
            h->SetBinContent(b, h->GetBinContent(b) / w);
            h->SetBinError  (b, h->GetBinError(b)   / w);
        }
    }
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

void styleHist(TH1D* h, Color_t c, Style_t m) {
    h->SetLineColor(c);
    h->SetMarkerColor(c);
    h->SetMarkerStyle(m);
    h->SetMarkerSize(0.8);
    h->SetLineWidth(2);
    h->GetXaxis()->SetTitleSize(0.050);
    h->GetYaxis()->SetTitleSize(0.050);
    h->GetXaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetLabelSize(0.045);
    h->GetXaxis()->SetTitleOffset(1.05);
    h->GetYaxis()->SetTitleOffset(1.20);
}

// Quadratic result for the MID piece, stitched at x=0.70 to HIGH's value.
//   y = high_join_y + b*(x - 0.70) + c*(x - 0.70)^2
// reported: free slope b, free curvature c, and the join-value constant.
struct QuadStitchedResult {
    double b, b_err;       // free slope at the stitch point
    double c, c_err;       // free curvature
    double high_join_y;    // numeric constant baked into the formula
    double chi2_ndf;
    int    ndf;
    TF1*   f;   // keep around so we can evaluate at the stitch point
};

// Free quadratic result for the HIGH piece: y = a + b*x + c*x*x
struct QuadResult {
    double a, a_err;
    double b, b_err;
    double c, c_err;
    double chi2_ndf;
    int    ndf;
    TF1*   f;
};

struct PieceFits {
    QuadStitchedResult mid;
    QuadResult         hi;
};

// FREE quadratic fit on the HIGH piece, with ω, φ and high-mass tail
// exclusion windows applied by zeroing the bins on a clone.
// The fit TF1 is NOT attached to the ORIGINAL histogram — it will be drawn
// explicitly on top of markers at draw time via TF1::Draw("SAME").
QuadResult fitQuadFreeWithExcl(TH1D* h_orig,
                               double xlo, double xhi,
                               Color_t color, const char* tag) {
    const std::string clone_name = std::string(h_orig->GetName()) + "_fitclone_" + tag;
    TH1D* h_fit = static_cast<TH1D*>(h_orig->Clone(clone_name.c_str()));
    h_fit->SetDirectory(nullptr);
    for (int b = 1; b <= h_fit->GetNbinsX(); ++b) {
        const double xc = h_fit->GetXaxis()->GetBinCenter(b);
        if ((xc >= kExcl1Lo && xc <= kExcl1Hi) ||
            (xc >= kExcl2Lo && xc <= kExcl2Hi) ||
            (xc >= kExcl3Lo && xc <= kExcl3Hi)) {
            h_fit->SetBinContent(b, 0.0);
            h_fit->SetBinError  (b, 0.0);
        }
    }

    const std::string fname = std::string(h_orig->GetName()) + "_fq_" + tag;
    TF1* f = new TF1(fname.c_str(), "[0]+[1]*x+[2]*x*x", xlo, xhi);
    f->SetLineColor(color);
    f->SetLineWidth(3);
    // "RQ0" — Range, Quiet, "0" prevents auto-draw and does NOT attach the
    // function to h_fit. We will Draw("SAME") the function ourselves at the
    // very end of the draw sequence so the line is on TOP of markers.
    h_fit->Fit(f, "RQ0");

    QuadResult r;
    r.a        = f->GetParameter(0);
    r.a_err    = f->GetParError(0);
    r.b        = f->GetParameter(1);
    r.b_err    = f->GetParError(1);
    r.c        = f->GetParameter(2);
    r.c_err    = f->GetParError(2);
    r.ndf      = f->GetNDF();
    r.chi2_ndf = (r.ndf > 0) ? f->GetChisquare() / r.ndf : 0.0;
    r.f        = f;
    return r;
}

// MID quadratic fit constrained to pass through (0.70, high_join_y):
//   y = high_join_y + b*(x - 0.70) + c*(x - 0.70)^2,   2 free params (b, c).
// Uses "RQ0" so the function is NOT auto-attached/drawn — we draw the TF1
// explicitly at the end of the draw sequence (after markers) so the line
// sits on TOP of the data markers.
QuadStitchedResult fitQuadStitchedRange(TH1D* h, double xlo, double xhi,
                                        double high_join_x, double high_join_y,
                                        Color_t color, const char* tag,
                                        const char* fit_opt) {
    // The stitch x is fixed at 0.70 per spec — build the formula with
    // high_join_y as a baked-in numeric constant.
    char mid_formula[256];
    std::snprintf(mid_formula, sizeof(mid_formula),
        "%.10g + [0]*(x - 0.70) + [1]*(x - 0.70)*(x - 0.70)",
        high_join_y);
    (void) high_join_x;  // stitch x is fixed at 0.70 by the formula above

    const std::string fname = std::string(h->GetName()) + "_fq_" + tag;
    TF1* f = new TF1(fname.c_str(), mid_formula, xlo, xhi);
    f->SetLineColor(color);
    f->SetLineWidth(3);
    h->Fit(f, fit_opt);
    QuadStitchedResult r;
    r.b           = f->GetParameter(0);
    r.b_err       = f->GetParError(0);
    r.c           = f->GetParameter(1);
    r.c_err       = f->GetParError(1);
    r.high_join_y = high_join_y;
    r.ndf         = f->GetNDF();
    r.chi2_ndf    = (r.ndf > 0) ? f->GetChisquare() / r.ndf : 0.0;
    r.f           = f;
    return r;
}

// Build the piecewise fits on a ratio histogram. LOW is no longer fitted.
// HIGH is fitted FIRST as a free quadratic; MID is then stitched to HIGH's
// value at x = 0.70.
PieceFits fitTwoPiece(TH1D* r, const std::string& tag_prefix) {
    PieceFits pf;
    pf.hi  = fitQuadFreeWithExcl(r, kHiLo, kHiHi, kColHigh,
                                 (tag_prefix + "hi").c_str());
    const double high_join_x = kHiLo;                 // 0.70
    const double high_join_y = pf.hi.f->Eval(high_join_x);
    // "RQ0" — Range, Quiet, "0" prevents auto-attach/draw. We Draw("SAME")
    // the TF1 explicitly at the end of the draw sequence so the line sits on
    // TOP of the data markers.
    pf.mid = fitQuadStitchedRange(r, kMdLo, kMdHi,
                                  high_join_x, high_join_y,
                                  kColMid,
                                  (tag_prefix + "mid").c_str(), "RQ0");
    // Defensive: if anything attached a function to r during fitting, strip
    // it — we want a clean histogram so the markers draw without TF1 below
    // them, and we will Draw("SAME") the TF1s ourselves at the very end.
    if (r->GetListOfFunctions()) {
        r->GetListOfFunctions()->Clear();
    }
    return pf;
}

// Draw the m_ee spectrum panel (top-left or bottom-left).
void drawSpectrumPanel(TH1D* h_PT3, TH1D* h_PT2,
                       const std::string& title,
                       double Y3, double Y2) {
    gPad->SetLogy(true);
    gPad->SetMargin(0.13, 0.04, 0.13, 0.10);
    const double y_max = std::max(h_PT3->GetMaximum(), h_PT2->GetMaximum());
    h_PT3->SetTitle(title.c_str());
    h_PT3->GetYaxis()->SetRangeUser(std::max(1e-2, y_max * 1e-6), y_max * 5.0);
    h_PT3->Draw("E1");
    h_PT2->Draw("E1 SAME");
    TLegend* leg = new TLegend(0.55, 0.74, 0.95, 0.88);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.034);
    leg->AddEntry(h_PT3, Form("PT3  (#sum w = %.0f)", Y3), "lpe");
    leg->AddEntry(h_PT2, Form("PT2  (#sum w = %.0f)", Y2), "lpe");
    leg->Draw();
}

// Draw the ratio panel with the MID + stitched-HIGH fits + ω/φ exclusion
// shading + annotations.
//
// CRITICAL draw order (so fit lines are NOT hidden under data markers):
//   1) r->Draw("E1")       — establish the frame/axes
//   2) TBoxes Draw("SAME") — gray exclusion windows under everything
//   3) r->Draw("E1 SAME")  — re-draw markers ON TOP of the boxes
//   4) f_mid->Draw("SAME") — fit lines drawn LAST so they sit ON TOP of
//   5) f_hi ->Draw("SAME")   the data markers
void drawRatioPanel(TH1D* r,
                    const std::string& title,
                    const PieceFits& pf) {
    gPad->SetMargin(0.13, 0.04, 0.13, 0.10);
    r->SetTitle(title.c_str());
    r->GetYaxis()->SetRangeUser(0.5, 2.0);

    // 1) Frame/axes first.
    r->Draw("E1");

    // 2) Gray exclusion windows.
    TBox* bx1 = new TBox(kExcl1Lo, 0.5, kExcl1Hi, 2.0);
    bx1->SetFillColorAlpha(kGray, 0.18);
    bx1->SetLineColor(0);
    bx1->Draw("SAME");
    TBox* bx2 = new TBox(kExcl2Lo, 0.5, kExcl2Hi, 2.0);
    bx2->SetFillColorAlpha(kGray, 0.18);
    bx2->SetLineColor(0);
    bx2->Draw("SAME");
    TBox* bx3 = new TBox(kExcl3Lo, 0.5, kExcl3Hi, 2.0);
    bx3->SetFillColorAlpha(kGray, 0.18);
    bx3->SetLineColor(0);
    bx3->Draw("SAME");

    // 3) Re-draw markers on TOP of the gray boxes.
    r->Draw("E1 SAME");

    // y=1 reference under fits but over markers — fine either way; place it
    // here so it doesn't visually compete with the fit curves below.
    TLine* lref = new TLine(kXmin, 1.0, kXmax, 1.0);
    lref->SetLineStyle(2); lref->SetLineColor(kGray + 2); lref->Draw();

    // 4) & 5) Fit lines LAST — they appear ON TOP of markers.
    if (pf.mid.f) pf.mid.f->Draw("SAME");
    if (pf.hi.f)  pf.hi.f->Draw("SAME");

    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.026);

    tx.SetTextColor(kColHigh);
    tx.DrawLatex(0.16, 0.88,
        Form("[%.2f-%.2f] quad FREE: a=%.3f b=%.3f c=%.3f  #chi^{2}/ndf = %.1f / %d  (excl #omega, #phi, hi)",
             kHiLo, kHiHi,
             pf.hi.a, pf.hi.b, pf.hi.c,
             pf.hi.chi2_ndf * pf.hi.ndf, pf.hi.ndf));

    tx.SetTextColor(kColMid);
    tx.DrawLatex(0.16, 0.84,
        Form("[%.2f-%.2f] quad STITCHED at 0.70: b=%.3f c=%.3f  #chi^{2}/ndf = %.1f / %d",
             kMdLo, kMdHi,
             pf.mid.b, pf.mid.c,
             pf.mid.chi2_ndf * pf.mid.ndf, pf.mid.ndf));
}

void printFits(const PieceFits& pf, const char* kind) {
    printf("  [%s]\n", kind);
    printf("  HIGH quad FREE [%.2f, %.2f] (excl ω,φ,hi):  a=%.4f±%.4f  b=%.4f±%.4f  c=%.4f±%.4f  chi2/ndf = %.1f / %d = %.2f\n",
           kHiLo, kHiHi,
           pf.hi.a, pf.hi.a_err, pf.hi.b, pf.hi.b_err, pf.hi.c, pf.hi.c_err,
           pf.hi.chi2_ndf * pf.hi.ndf, pf.hi.ndf, pf.hi.chi2_ndf);
    printf("  MID  quad stch [%.2f, %.2f] (HIGH(0.70)=%.4f):  b=%.4f±%.4f  c=%.4f±%.4f  chi2/ndf = %.1f / %d = %.2f\n",
           kMdLo, kMdHi, pf.mid.high_join_y,
           pf.mid.b, pf.mid.b_err, pf.mid.c, pf.mid.c_err,
           pf.mid.chi2_ndf * pf.mid.ndf, pf.mid.ndf, pf.mid.chi2_ndf);
}

void drawOverlay(const std::vector<TH1D*>& ratios,
                 const std::vector<PieceFits>& fits,
                 const std::string& title,
                 const std::string& canvas_name,
                 const std::string& out_base) {
    TCanvas* c_ovl = new TCanvas(canvas_name.c_str(),
        title.c_str(), 1500, 700);
    gPad->SetMargin(0.13, 0.04, 0.13, 0.10);
    ratios[0]->SetTitle(title.c_str());
    ratios[0]->GetYaxis()->SetTitleOffset(1.05);
    ratios[0]->GetYaxis()->SetRangeUser(0.5, 2.5);
    ratios[0]->Draw("E1");
    for (size_t i = 1; i < ratios.size(); ++i)
        ratios[i]->Draw("E1 SAME");

    TLine* lref = new TLine(kXmin, 1.0, kXmax, 1.0);
    lref->SetLineStyle(2); lref->SetLineColor(kGray + 2); lref->Draw();

    TLegend* leg = new TLegend(0.55, 0.62, 0.96, 0.88);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.022);
    for (size_t i = 0; i < kScenarios.size(); ++i) {
        leg->AddEntry(ratios[i],
            Form("%s   MID: b=%.3f c=%.3f #chi^{2}/ndf=%.1f",
                 kScenarios[i].label.c_str(),
                 fits[i].mid.b, fits[i].mid.c, fits[i].mid.chi2_ndf),
            "lpe");
    }
    leg->Draw();

    c_ovl->SaveAs((out_base + ".pdf").c_str());
    c_ovl->SaveAs((out_base + ".png").c_str());
}

}  // anonymous namespace

void mass_ee_pt_ratio_purity_cascade_sim_3piece_v5() {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetErrorX(0.0);

    TFile* fs = TFile::Open(kInputFile, "READ");
    if (!fs || fs->IsZombie()) {
        std::cerr << "Cannot open " << kInputFile << "\n"; return;
    }
    TTree* nt_rec   = dynamic_cast<TTree*>(fs->Get(kNtRec));
    TTree* nt_truth = dynamic_cast<TTree*>(fs->Get(kNtTruth));
    if (!nt_rec) {
        std::cerr << kNtRec << " missing in " << kInputFile << "\n"; return;
    }
    if (!nt_truth) {
        std::cerr << kNtTruth << " missing in " << kInputFile << "\n"; return;
    }
    std::cout << "Input: " << kInputFile << "\n";
    std::cout << "  REC tree   : " << kNtRec << " (" << nt_rec->GetEntries() << " entries)\n";
    std::cout << "  TRUTH tree : " << kNtTruth << " (" << nt_truth->GetEntries() << " entries)\n";
    std::cout << "Binning (variable, " << kNb << " bins):\n";
    std::cout << "  [0.00, 0.15] :  60 x 2.5 MeV\n";
    std::cout << "  [0.15, 0.85] : 140 x 5 MeV\n";
    std::cout << "  [0.85, 1.10] :  12 x ~20.833 MeV\n";
    std::cout << "  [1.10, 1.40] :   8 x 37.5 MeV\n";
    std::cout << "Fit: HIGH quad FREE ["
              << kHiLo << "," << kHiHi << "]"
              << " (excl [" << kExcl1Lo << "," << kExcl1Hi << "], ["
              << kExcl2Lo << "," << kExcl2Hi << "] and ["
              << kExcl3Lo << "," << kExcl3Hi << "])"
              << "  MID quad STITCHED at " << kHiLo << " ["
              << kMdLo << "," << kMdHi << "]\n";
    std::cout << "LOW [0.00, 0.15]: no fit (data points only)\n";
    std::cout << "Spectrum panels normalized by bin width (entries / GeV/c²)\n";
    std::cout << "Weight: sim_genweight × cumulative cut\n\n";

    gSystem->mkdir("plots/output", true);

    std::vector<TH1D*> ratios_rec_overlay, ratios_truth_overlay;
    std::vector<PieceFits> fits_rec, fits_truth;
    std::vector<double> yields_PT3_rec, yields_PT2_rec;
    std::vector<double> yields_PT3_truth, yields_PT2_truth;

    for (size_t i = 0; i < kScenarios.size(); ++i) {
        const auto& S = kScenarios[i];
        std::cout << "\n=== Scenario " << S.label << " ===\n";

        const std::string w_PT3 = "(pt3==1)*" + S.cut + "*sim_genweight";
        const std::string w_PT2 = "(pt2==1)*" + S.cut + "*sim_genweight";
        std::cout << "  PT3 weight: " << w_PT3 << "\n";

        // ---- REC: m_ee from dilepton_nt ----
        TH1D* h_PT3_rec = drawWithWeight(nt_rec, "m_ee", w_PT3,
            Form("h_PT3_rec_%s", S.id.c_str()));
        TH1D* h_PT2_rec = drawWithWeight(nt_rec, "m_ee", w_PT2,
            Form("h_PT2_rec_%s", S.id.c_str()));

        // ---- TRUTH: m_ee_sim from dilepton_nt_cor ----
        TH1D* h_PT3_truth = drawWithWeight(nt_truth, "m_ee_sim", w_PT3,
            Form("h_PT3_truth_%s", S.id.c_str()));
        TH1D* h_PT2_truth = drawWithWeight(nt_truth, "m_ee_sim", w_PT2,
            Form("h_PT2_truth_%s", S.id.c_str()));

        // Yields computed BEFORE width normalization (integral of raw entries).
        const double Y3_rec = h_PT3_rec->Integral();
        const double Y2_rec = h_PT2_rec->Integral();
        const double Y3_tr  = h_PT3_truth->Integral();
        const double Y2_tr  = h_PT2_truth->Integral();
        std::cout << "  REC   yields: PT3 = " << Y3_rec << "  PT2 = " << Y2_rec << "\n";
        std::cout << "  TRUTH yields: PT3 = " << Y3_tr  << "  PT2 = " << Y2_tr  << "\n";
        yields_PT3_rec.push_back(Y3_rec);
        yields_PT2_rec.push_back(Y2_rec);
        yields_PT3_truth.push_back(Y3_tr);
        yields_PT2_truth.push_back(Y2_tr);

        // Ratio built from RAW (non-width-normalized) histograms — the ratio
        // is unit-less and unaffected by width division anyway.
        TH1D* r_rec   = makeRatio(h_PT2_rec,   h_PT3_rec,
            Form("r_rec_%s", S.id.c_str()));
        TH1D* r_truth = makeRatio(h_PT2_truth, h_PT3_truth,
            Form("r_truth_%s", S.id.c_str()));

        // Now apply width normalization to the spectrum histograms only.
        normalizeByBinWidth(h_PT3_rec);
        normalizeByBinWidth(h_PT2_rec);
        normalizeByBinWidth(h_PT3_truth);
        normalizeByBinWidth(h_PT2_truth);

        styleHist(h_PT3_rec,   kColPT3, 20);
        styleHist(h_PT2_rec,   kColPT2, 21);
        styleHist(h_PT3_truth, kColPT3, 20);
        styleHist(h_PT2_truth, kColPT2, 21);
        styleHist(r_rec,   kColRatio, 20);
        styleHist(r_truth, kColRatio, 20);

        PieceFits pf_rec   = fitTwoPiece(r_rec,   "rec_");
        PieceFits pf_truth = fitTwoPiece(r_truth, "truth_");

        // Stitch-continuity sanity check for the very first scenario / REC.
        if (i == 0) {
            const double high70 = pf_rec.hi.f->Eval(kHiLo);
            const double mid70  = pf_rec.mid.f->Eval(kHiLo);
            std::cout << "  [stitch sanity step1 REC] HIGH(0.70) = " << high70
                      << "   MID(0.70) = " << mid70
                      << "   diff = " << (mid70 - high70) << "\n";
        }

        printFits(pf_rec,   "REC");
        printFits(pf_truth, "TRUTH");
        fits_rec.push_back(pf_rec);
        fits_truth.push_back(pf_truth);

        // -------- Canvas: 2×2 --------
        TCanvas* c = new TCanvas(
            Form("c_purity_cascade_3piece_v5_%s", S.id.c_str()),
            S.label.c_str(),
            1500, 900);
        c->Divide(2, 2, 0.001, 0.001);

        // Top-left: REC spectrum
        c->cd(1);
        drawSpectrumPanel(h_PT3_rec, h_PT2_rec,
            Form("m_{ee} REC — %s;M_{e^{+}e^{-}} [GeV/c^{2}];weighted entries / GeV/c^{2}",
                 S.label.c_str()),
            Y3_rec, Y2_rec);

        // Top-right: REC ratio
        c->cd(2);
        drawRatioPanel(r_rec,
            Form("N_{PT2}/N_{PT3} REC — %s;M_{e^{+}e^{-}} [GeV/c^{2}];N_{PT2}/N_{PT3}",
                 S.label.c_str()),
            pf_rec);

        // Bottom-left: TRUTH spectrum
        c->cd(3);
        drawSpectrumPanel(h_PT3_truth, h_PT2_truth,
            Form("m_{ee} TRUTH (m_{ee}^{sim}) — %s;M_{e^{+}e^{-}}^{sim} [GeV/c^{2}];weighted entries / GeV/c^{2}",
                 S.label.c_str()),
            Y3_tr, Y2_tr);

        // Bottom-right: TRUTH ratio
        c->cd(4);
        drawRatioPanel(r_truth,
            Form("N_{PT2}/N_{PT3} TRUTH — %s;M_{e^{+}e^{-}}^{sim} [GeV/c^{2}];N_{PT2}/N_{PT3}",
                 S.label.c_str()),
            pf_truth);

        // ----- Save -----
        const std::string base =
            std::string("plots/output/mass_ee_pt_ratio_purity_cascade_sim_3piece_v5_") + S.id;
        c->SaveAs((base + ".pdf").c_str());
        c->SaveAs((base + ".png").c_str());
        std::cout << "  saved: " << base << ".{pdf,png}\n";

        // Store ratios for overlay canvases
        TH1D* r_rec_ovl = static_cast<TH1D*>(r_rec->Clone(
            Form("r_rec_ovl_%s", S.id.c_str())));
        r_rec_ovl->SetDirectory(nullptr);
        r_rec_ovl->GetListOfFunctions()->Clear();
        styleHist(r_rec_ovl, kOverlayColors[i], 20 + (int)i);
        ratios_rec_overlay.push_back(r_rec_ovl);

        TH1D* r_truth_ovl = static_cast<TH1D*>(r_truth->Clone(
            Form("r_truth_ovl_%s", S.id.c_str())));
        r_truth_ovl->SetDirectory(nullptr);
        r_truth_ovl->GetListOfFunctions()->Clear();
        styleHist(r_truth_ovl, kOverlayColors[i], 20 + (int)i);
        ratios_truth_overlay.push_back(r_truth_ovl);
    }

    // ============================================================================
    // Overlay canvases — REC and TRUTH separately
    // ============================================================================
    drawOverlay(ratios_rec_overlay, fits_rec,
        "PT2/PT3 cumulative purity cascade — sim REC (v5 HIGH-free, MID-stitched);"
        "M_{e^{+}e^{-}} [GeV/c^{2}];N_{PT2}/N_{PT3}",
        "c_purity_cascade_3piece_v5_overlay_rec",
        "plots/output/mass_ee_pt_ratio_purity_cascade_sim_3piece_v5_overlay_rec");

    drawOverlay(ratios_truth_overlay, fits_truth,
        "PT2/PT3 cumulative purity cascade — sim TRUTH (v5 HIGH-free, MID-stitched);"
        "M_{e^{+}e^{-}}^{sim} [GeV/c^{2}];N_{PT2}/N_{PT3}",
        "c_purity_cascade_3piece_v5_overlay_truth",
        "plots/output/mass_ee_pt_ratio_purity_cascade_sim_3piece_v5_overlay_truth");

    // ============================================================================
    // Summary
    // ============================================================================
    std::cout << "\n\n=== SUMMARY (piecewise fits, v5 flipped stitch — REC vs TRUTH) ===\n";
    std::cout << "  LOW  [0.00, 0.15] : no fit (data points only)\n";
    std::cout << "  HIGH quad FREE  y=a+bx+cx^2 in [" << kHiLo << ", " << kHiHi
              << "]  excluding [" << kExcl1Lo << "," << kExcl1Hi << "], ["
              << kExcl2Lo << "," << kExcl2Hi << "] and ["
              << kExcl3Lo << "," << kExcl3Hi << "]\n";
    std::cout << "  MID  quad STITCHED y=HIGH(0.70)+b·(x-0.70)+c·(x-0.70)^2 in ["
              << kMdLo << ", " << kMdHi << "]\n\n";

    for (size_t i = 0; i < kScenarios.size(); ++i) {
        std::cout << "  [" << kScenarios[i].label << "]\n";
        std::cout << "    REC   yields: PT3=" << yields_PT3_rec[i]
                  << "  PT2=" << yields_PT2_rec[i] << "\n";
        std::cout << "    TRUTH yields: PT3=" << yields_PT3_truth[i]
                  << "  PT2=" << yields_PT2_truth[i] << "\n";

        printf("    REC   HIGH(FREE): a=%+.4f b=%+.4f c=%+.4f chi2/ndf=%.2f   MID(stitched quad): b=%+.4f c=%+.4f chi2/ndf=%.2f  HIGH(0.70)=%.4f\n",
               fits_rec[i].hi.a, fits_rec[i].hi.b, fits_rec[i].hi.c, fits_rec[i].hi.chi2_ndf,
               fits_rec[i].mid.b, fits_rec[i].mid.c, fits_rec[i].mid.chi2_ndf,
               fits_rec[i].mid.high_join_y);
        printf("    TRUTH HIGH(FREE): a=%+.4f b=%+.4f c=%+.4f chi2/ndf=%.2f   MID(stitched quad): b=%+.4f c=%+.4f chi2/ndf=%.2f  HIGH(0.70)=%.4f\n",
               fits_truth[i].hi.a, fits_truth[i].hi.b, fits_truth[i].hi.c, fits_truth[i].hi.chi2_ndf,
               fits_truth[i].mid.b, fits_truth[i].mid.c, fits_truth[i].mid.chi2_ndf,
               fits_truth[i].mid.high_join_y);
    }
    std::cout << "\n";

    fs->Close();
}
