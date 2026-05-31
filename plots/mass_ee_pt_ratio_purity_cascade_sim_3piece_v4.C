// mass_ee_pt_ratio_purity_cascade_sim_3piece_v4.C
// =========================================================================
// V4 of the three-piece piecewise fit purity cascade study (sim).
//
// IDENTICAL fit setup to v3:
//
//   piece   | range [GeV/c²]  | shape       | comment
//   --------+-----------------+-------------+------------------------------
//   LOW     | [0.00, 0.20]    | pol5        | smooth description only
//   MID     | [0.20, 0.70]    | linear      | y = a + b*x
//   HIGH    | [0.65, 1.40]    | quadratic   | y = a + b*x + c*x^2 with
//                                           | EXCLUSIONS of ω [0.75,0.83]
//                                           | and φ-region [1.00,1.15]
//
// CHANGES vs v3:
//   - 5 CUMULATIVE cut stages (no separate step0):
//       step1: RAW (no cuts)
//       step2: +eVertReco_z > -500
//       step3: +isBest == 1
//       step4: +ep_sim_id==2 && em_sim_id==3
//       step5: +epem_same_vertex==1 (full purity gate)
//   - Each per-step canvas is now 2×2 (1500×900): REC (top) + TRUTH (bottom).
//     Left column = m_ee spectrum, right column = ratio with 3-piece fits.
//     REC uses m_ee from dilepton_nt; TRUTH uses m_ee_sim from dilepton_nt_cor.
//     Cuts are applied identically to both variants.
//   - Two overlay canvases instead of one:
//       overlay_rec   : all 5 REC ratios overlaid
//       overlay_truth : all 5 TRUTH ratios overlaid
//
// Usage:
//   root -l -b -q plots/mass_ee_pt_ratio_purity_cascade_sim_3piece_v4.C
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
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace {

constexpr const char* kInputFile = "output_epem_sim.root";
constexpr const char* kNtRec     = "dilepton_nt";
constexpr const char* kNtTruth   = "dilepton_nt_cor";

constexpr int    kNb    = 70;
constexpr double kXmin  = 0.0;
constexpr double kXmax  = 1.4;

// Three-piece fit ranges (GeV/c²) — same as v3
constexpr double kLoLo  = 0.00, kLoHi = 0.20;   // LOW   piece — pol5
constexpr double kMdLo  = 0.20, kMdHi = 0.70;   // MID   piece — linear
constexpr double kHiLo  = 0.65, kHiHi = 1.40;   // HIGH  piece — quadratic

// Exclusion windows applied to HIGH fit (ω and φ regions)
constexpr double kExcl1Lo = 0.75, kExcl1Hi = 0.83;   // ω peak
constexpr double kExcl2Lo = 1.00, kExcl2Hi = 1.15;   // φ peak + wing

constexpr Color_t kColPT3   = kBlack;
constexpr Color_t kColPT2   = kRed + 1;
constexpr Color_t kColRatio = kBlue + 1;

// Distinct color per piece
const Color_t kColLow  = kViolet + 2;
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

// Same palette as v3 overlay (excluding the gray that was step0).
const std::array<Color_t, 5> kOverlayColors = {
    kRed + 1, kOrange + 7, kGreen + 2, kAzure + 1, kViolet + 2
};

TH1D* drawWithWeight(TTree* t, const std::string& mass_var,
                     const std::string& weight, const std::string& name) {
    TH1D* h = new TH1D(name.c_str(), "", kNb, kXmin, kXmax);
    h->Sumw2();
    t->Draw((mass_var + ">>" + name).c_str(),
            weight.c_str(), "goff");
    h->SetDirectory(nullptr);
    return h;
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

// pol5 result: y = sum_{k=0..5} p_k * x^k
struct Pol5Result {
    double p[6];
    double p_err[6];
    double chi2_ndf;
    int    ndf;
};

// Quadratic result: y = a + b*x + c*x^2
struct QuadResult {
    double a, a_err;
    double b, b_err;
    double c, c_err;
    double chi2_ndf;
    int    ndf;
};

// Linear result: y = a + b*x
struct LinResult {
    double a, a_err;
    double b, b_err;
    double chi2_ndf;
    int    ndf;
};

struct PieceFits {
    Pol5Result  lo;
    LinResult   mid;
    QuadResult  hi;
};

Pol5Result fitPol5Range(TH1D* h, double xlo, double xhi,
                        Color_t color, const char* tag,
                        const char* fit_opt) {
    const std::string fname = std::string(h->GetName()) + "_fp5_" + tag;
    TF1* f = new TF1(fname.c_str(), "pol5", xlo, xhi);
    f->SetLineColor(color);
    f->SetLineWidth(3);
    h->Fit(f, fit_opt);
    Pol5Result r;
    for (int k = 0; k < 6; ++k) {
        r.p[k]     = f->GetParameter(k);
        r.p_err[k] = f->GetParError(k);
    }
    r.ndf      = f->GetNDF();
    r.chi2_ndf = (r.ndf > 0) ? f->GetChisquare() / r.ndf : 0.0;
    return r;
}

LinResult fitLinRange(TH1D* h, double xlo, double xhi,
                      Color_t color, const char* tag,
                      const char* fit_opt) {
    const std::string fname = std::string(h->GetName()) + "_fl_" + tag;
    TF1* f = new TF1(fname.c_str(), "[0]+[1]*x", xlo, xhi);
    f->SetLineColor(color);
    f->SetLineWidth(3);
    h->Fit(f, fit_opt);
    LinResult r;
    r.a        = f->GetParameter(0);
    r.a_err    = f->GetParError(0);
    r.b        = f->GetParameter(1);
    r.b_err    = f->GetParError(1);
    r.ndf      = f->GetNDF();
    r.chi2_ndf = (r.ndf > 0) ? f->GetChisquare() / r.ndf : 0.0;
    return r;
}

// Quadratic fit performed on a masked clone of the histogram, then
// re-attached to the ORIGINAL histogram for drawing.
QuadResult fitQuadRangeWithExcl(TH1D* h_orig, double xlo, double xhi,
                                Color_t color, const char* tag) {
    const std::string clone_name = std::string(h_orig->GetName()) + "_fitclone_" + tag;
    TH1D* h_fit = static_cast<TH1D*>(h_orig->Clone(clone_name.c_str()));
    h_fit->SetDirectory(nullptr);
    for (int b = 1; b <= h_fit->GetNbinsX(); ++b) {
        const double xc = h_fit->GetXaxis()->GetBinCenter(b);
        if ((xc >= kExcl1Lo && xc <= kExcl1Hi) ||
            (xc >= kExcl2Lo && xc <= kExcl2Hi)) {
            h_fit->SetBinContent(b, 0.0);
            h_fit->SetBinError  (b, 0.0);
        }
    }

    const std::string fname = std::string(h_orig->GetName()) + "_fq_" + tag;
    TF1* f = new TF1(fname.c_str(), "[0]+[1]*x+[2]*x*x", xlo, xhi);
    f->SetLineColor(color);
    f->SetLineWidth(3);
    h_fit->Fit(f, "RQ");

    TF1* f_draw = static_cast<TF1*>(f->Clone((fname + "_draw").c_str()));
    f_draw->SetLineColor(color);
    f_draw->SetLineWidth(3);
    h_orig->GetListOfFunctions()->Add(f_draw);

    QuadResult r;
    r.a        = f->GetParameter(0);
    r.a_err    = f->GetParError(0);
    r.b        = f->GetParameter(1);
    r.b_err    = f->GetParError(1);
    r.c        = f->GetParameter(2);
    r.c_err    = f->GetParError(2);
    r.ndf      = f->GetNDF();
    r.chi2_ndf = (r.ndf > 0) ? f->GetChisquare() / r.ndf : 0.0;
    return r;
}

// Build the three-piece fits on a ratio histogram, with a unique tag so
// REC and TRUTH function names do not collide.
PieceFits fitThreePiece(TH1D* r, const std::string& tag_prefix) {
    PieceFits pf;
    pf.lo  = fitPol5Range(r, kLoLo, kLoHi, kColLow,
                          (tag_prefix + "lo").c_str(),  "RQ");
    pf.mid = fitLinRange (r, kMdLo, kMdHi, kColMid,
                          (tag_prefix + "mid").c_str(), "RQ+");
    pf.hi  = fitQuadRangeWithExcl(r, kHiLo, kHiHi, kColHigh,
                                  (tag_prefix + "hi").c_str());
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

// Draw the ratio panel with the three-piece fits + ω/φ exclusion shading
// + annotations. The TF1s are already attached to r.
void drawRatioPanel(TH1D* r,
                    const std::string& title,
                    const PieceFits& pf) {
    gPad->SetMargin(0.13, 0.04, 0.13, 0.10);
    r->SetTitle(title.c_str());
    r->GetYaxis()->SetRangeUser(0.5, 2.0);
    r->Draw("E1");

    TBox* bx1 = new TBox(kExcl1Lo, 0.5, kExcl1Hi, 2.0);
    bx1->SetFillColorAlpha(kGray, 0.18);
    bx1->SetLineColor(0);
    bx1->Draw("SAME");
    TBox* bx2 = new TBox(kExcl2Lo, 0.5, kExcl2Hi, 2.0);
    bx2->SetFillColorAlpha(kGray, 0.18);
    bx2->SetLineColor(0);
    bx2->Draw("SAME");

    r->Draw("E1 SAME");

    TLine* lref = new TLine(kXmin, 1.0, kXmax, 1.0);
    lref->SetLineStyle(2); lref->SetLineColor(kGray + 2); lref->Draw();

    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.026);

    tx.SetTextColor(kColLow);
    tx.DrawLatex(0.16, 0.88,
        Form("[%.2f-%.2f] pol5:   #chi^{2}/ndf = %.1f / %d",
             kLoLo, kLoHi,
             pf.lo.chi2_ndf * pf.lo.ndf, pf.lo.ndf));

    tx.SetTextColor(kColMid);
    tx.DrawLatex(0.16, 0.84,
        Form("[%.2f-%.2f] lin:    a=%.3f b=%.3f    #chi^{2}/ndf = %.1f / %d",
             kMdLo, kMdHi,
             pf.mid.a, pf.mid.b,
             pf.mid.chi2_ndf * pf.mid.ndf, pf.mid.ndf));

    tx.SetTextColor(kColHigh);
    tx.DrawLatex(0.16, 0.80,
        Form("[%.2f-%.2f] quad:   a=%.3f b=%.3f c=%.3f    #chi^{2}/ndf = %.1f / %d  (excl #omega, #phi)",
             kHiLo, kHiHi,
             pf.hi.a, pf.hi.b, pf.hi.c,
             pf.hi.chi2_ndf * pf.hi.ndf, pf.hi.ndf));
}

void printFits(const PieceFits& pf, const char* kind) {
    printf("  [%s]\n", kind);
    printf("  LOW  pol5 [%.2f, %.2f]:  chi2/ndf = %.1f / %d = %.2f\n",
           kLoLo, kLoHi,
           pf.lo.chi2_ndf * pf.lo.ndf, pf.lo.ndf, pf.lo.chi2_ndf);
    printf("           p0=%+.4f±%.4f p1=%+.4f±%.4f p2=%+.4f±%.4f\n",
           pf.lo.p[0], pf.lo.p_err[0], pf.lo.p[1], pf.lo.p_err[1],
           pf.lo.p[2], pf.lo.p_err[2]);
    printf("           p3=%+.4f±%.4f p4=%+.4f±%.4f p5=%+.4f±%.4f\n",
           pf.lo.p[3], pf.lo.p_err[3], pf.lo.p[4], pf.lo.p_err[4],
           pf.lo.p[5], pf.lo.p_err[5]);
    printf("  MID  lin  [%.2f, %.2f]:  a=%.4f±%.4f  b=%.4f±%.4f                       chi2/ndf = %.1f / %d = %.2f\n",
           kMdLo, kMdHi,
           pf.mid.a, pf.mid.a_err, pf.mid.b, pf.mid.b_err,
           pf.mid.chi2_ndf * pf.mid.ndf, pf.mid.ndf, pf.mid.chi2_ndf);
    printf("  HIGH quad [%.2f, %.2f] (excl ω,φ):  a=%.4f±%.4f  b=%.4f±%.4f  c=%.4f±%.4f  chi2/ndf = %.1f / %d = %.2f\n",
           kHiLo, kHiHi,
           pf.hi.a, pf.hi.a_err, pf.hi.b, pf.hi.b_err, pf.hi.c, pf.hi.c_err,
           pf.hi.chi2_ndf * pf.hi.ndf, pf.hi.ndf, pf.hi.chi2_ndf);
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
            Form("%s   MID: a=%.3f b=%.3f #chi^{2}/ndf=%.1f",
                 kScenarios[i].label.c_str(),
                 fits[i].mid.a, fits[i].mid.b, fits[i].mid.chi2_ndf),
            "lpe");
    }
    leg->Draw();

    c_ovl->SaveAs((out_base + ".pdf").c_str());
    c_ovl->SaveAs((out_base + ".png").c_str());
}

}  // anonymous namespace

void mass_ee_pt_ratio_purity_cascade_sim_3piece_v4() {
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
    std::cout << "Fit: three pieces — LOW pol5 ["
              << kLoLo << "," << kLoHi << "]  MID lin ["
              << kMdLo << "," << kMdHi << "]  HIGH quad ["
              << kHiLo << "," << kHiHi << "]"
              << " (excl [" << kExcl1Lo << "," << kExcl1Hi << "] and ["
              << kExcl2Lo << "," << kExcl2Hi << "])\n";
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

        TH1D* r_rec   = makeRatio(h_PT2_rec,   h_PT3_rec,
            Form("r_rec_%s", S.id.c_str()));
        TH1D* r_truth = makeRatio(h_PT2_truth, h_PT3_truth,
            Form("r_truth_%s", S.id.c_str()));

        styleHist(h_PT3_rec,   kColPT3, 20);
        styleHist(h_PT2_rec,   kColPT2, 21);
        styleHist(h_PT3_truth, kColPT3, 20);
        styleHist(h_PT2_truth, kColPT2, 21);
        styleHist(r_rec,   kColRatio, 20);
        styleHist(r_truth, kColRatio, 20);

        PieceFits pf_rec   = fitThreePiece(r_rec,   "rec_");
        PieceFits pf_truth = fitThreePiece(r_truth, "truth_");

        printFits(pf_rec,   "REC");
        printFits(pf_truth, "TRUTH");
        fits_rec.push_back(pf_rec);
        fits_truth.push_back(pf_truth);

        // -------- Canvas: 2×2 --------
        TCanvas* c = new TCanvas(
            Form("c_purity_cascade_3piece_v4_%s", S.id.c_str()),
            S.label.c_str(),
            1500, 900);
        c->Divide(2, 2, 0.001, 0.001);

        // Top-left: REC spectrum
        c->cd(1);
        drawSpectrumPanel(h_PT3_rec, h_PT2_rec,
            Form("m_{ee} REC — %s;M_{e^{+}e^{-}} [GeV/c^{2}];weighted entries / 20 MeV",
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
            Form("m_{ee} TRUTH (m_{ee}^{sim}) — %s;M_{e^{+}e^{-}}^{sim} [GeV/c^{2}];weighted entries / 20 MeV",
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
            std::string("plots/output/mass_ee_pt_ratio_purity_cascade_sim_3piece_v4_") + S.id;
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
        "PT2/PT3 cumulative purity cascade — sim REC (3-piece v4 fit);"
        "M_{e^{+}e^{-}} [GeV/c^{2}];N_{PT2}/N_{PT3}",
        "c_purity_cascade_3piece_v4_overlay_rec",
        "plots/output/mass_ee_pt_ratio_purity_cascade_sim_3piece_v4_overlay_rec");

    drawOverlay(ratios_truth_overlay, fits_truth,
        "PT2/PT3 cumulative purity cascade — sim TRUTH (3-piece v4 fit);"
        "M_{e^{+}e^{-}}^{sim} [GeV/c^{2}];N_{PT2}/N_{PT3}",
        "c_purity_cascade_3piece_v4_overlay_truth",
        "plots/output/mass_ee_pt_ratio_purity_cascade_sim_3piece_v4_overlay_truth");

    // ============================================================================
    // Summary
    // ============================================================================
    std::cout << "\n\n=== SUMMARY (three piecewise fits, v4 — REC vs TRUTH) ===\n";
    std::cout << "  LOW  pol5            in [" << kLoLo << ", " << kLoHi << "]\n";
    std::cout << "  MID  lin  y=a+bx     in [" << kMdLo << ", " << kMdHi << "]\n";
    std::cout << "  HIGH quad y=a+bx+cx² in [" << kHiLo << ", " << kHiHi
              << "]  excluding [" << kExcl1Lo << "," << kExcl1Hi << "] and ["
              << kExcl2Lo << "," << kExcl2Hi << "]\n\n";

    for (size_t i = 0; i < kScenarios.size(); ++i) {
        std::cout << "  [" << kScenarios[i].label << "]\n";
        std::cout << "    REC   yields: PT3=" << yields_PT3_rec[i]
                  << "  PT2=" << yields_PT2_rec[i] << "\n";
        std::cout << "    TRUTH yields: PT3=" << yields_PT3_truth[i]
                  << "  PT2=" << yields_PT2_truth[i] << "\n";

        printf("    REC   LOW : chi2/ndf = %.2f   MID : a=%+.4f b=%+.4f chi2/ndf=%.2f   HIGH: a=%+.4f b=%+.4f c=%+.4f chi2/ndf=%.2f\n",
               fits_rec[i].lo.chi2_ndf,
               fits_rec[i].mid.a, fits_rec[i].mid.b, fits_rec[i].mid.chi2_ndf,
               fits_rec[i].hi.a, fits_rec[i].hi.b, fits_rec[i].hi.c, fits_rec[i].hi.chi2_ndf);
        printf("    TRUTH LOW : chi2/ndf = %.2f   MID : a=%+.4f b=%+.4f chi2/ndf=%.2f   HIGH: a=%+.4f b=%+.4f c=%+.4f chi2/ndf=%.2f\n",
               fits_truth[i].lo.chi2_ndf,
               fits_truth[i].mid.a, fits_truth[i].mid.b, fits_truth[i].mid.chi2_ndf,
               fits_truth[i].hi.a, fits_truth[i].hi.b, fits_truth[i].hi.c, fits_truth[i].hi.chi2_ndf);
    }
    std::cout << "\n";

    fs->Close();
}
