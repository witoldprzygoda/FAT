// mass_ee_pt_ratio_signal_varbins_exp_step5pts.C
// =========================================================================
// EXP analysis: 63 · N_PT2 / N_PT3 ratio with variable binning (signal).
// Like mass_ee_pt_ratio_signal_varbins_exp_step5.C, BUT the sim shape
// is taken from MEASURED step5 sim PT2/PT3 RATIO POINTS (cascade-v5 220-bin
// variable layout) — NOT from any analytic fit.
//
// The sim ratio shape is then globally scaled by ONE factor k_simpts:
//
//   k_simpts = I_expRatio[0.15, 0.70] / I_simpoints[0.15, 0.70]
//
// where I_simpoints is computed via Σ (sim_ratio(b) · binwidth(b)) over
// sim bins fully contained in [0.15, 0.70].
//
// Left panel:
//   PT3 black markers (raw)
//   63 · PT2 red filled circles (measured)
//   63 · PT2 hypothetical BLUE OPEN SQUARES = PT3(x) · k_simpts ·
//     sim_ratio_at(x) on the FULL range [0, 1.4] (NOT restricted),
//     where sim_ratio_at(x) is nearest-neighbour lookup of the
//     scaled sim ratio TH1.
//
// Right panel:
//   exp signal ratio markers (measured, 58-bin layout)
//   native MID free quad DASHED on [0.10, 0.80] (reference only)
//   scaled sim POINTS (220-bin layout) overlaid as gray triangles
//   red open sim circles in [0.00, 0.15] (same scaling) as separate
//     annotation of the LOW π⁰ region
//   y range [0.5, 4.0]
//
// Output:
//   plots/output/mass_ee_pt_ratio_varbins_exp_signal_step5pts.{pdf,png}
//
// Step5 (same_vertex) cumulative purity gate:
//   "(eVertReco_z>-500)*(isBest==1)*(ep_sim_id==2)*(em_sim_id==3)*"
//   "(epem_same_vertex==1)"
//
// Usage:
//   root -l -b -q plots/mass_ee_pt_ratio_signal_varbins_exp_step5pts.C
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
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

namespace {

constexpr const char* kFileEm = "output_epem_exp.root";
constexpr const char* kFilePp = "output_epep_exp.root";
constexpr const char* kFileMm = "output_emem_exp.root";
constexpr const char* kTree   = "dilepton_nt";

// SIM inputs — provides the sim ratio POINTS used to drive the correction.
constexpr const char* kFileSim    = "output_epem_sim.root";
constexpr const char* kTreeSimRec = "dilepton_nt";

// Step5 (same_vertex) cumulative purity gate.
constexpr const char* kSimStep5Cut =
    "(eVertReco_z>-500)*(isBest==1)*(ep_sim_id==2)*(em_sim_id==3)*"
    "(epem_same_vertex==1)";

// SIM variable binning (cascade v5: 220 bins on [0, 1.4]).
//   [0.00, 0.15] :  60 bins of 2.5 MeV
//   [0.15, 0.85] : 140 bins of 5 MeV
//   [0.85, 1.10] :  12 bins of (0.25/12) GeV
//   [1.10, 1.40] :   8 bins of (0.30/8) GeV
constexpr int kSimNb = 220;
const double* getSimEdges() {
    static double edges[kSimNb + 1];
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
// Minimum sim PT3 counts per bin to count the bin as defined.
constexpr double kSimMinDen = 5.0;

// EXP 58-bin variable layout (matches the existing ratio macros).
constexpr int    kNb   = 58;
constexpr double kXmin = 0.0;
constexpr double kXmax = 1.4;

// Trigger downscale: PT2 is downscaled by 64, so we multiply by 63.
constexpr double kTrigCorr = 63.0;

// Trigger bit cuts (EXP).
constexpr const char* kCutPT3 = "trigbit==8192";
constexpr const char* kCutPT2 = "trigbit==4096";

// Native MID free quad reference fit range.
constexpr double kMdLo = 0.10, kMdHi = 0.80;

// Integral-matching range for k_simpts.
constexpr double kIntLo = 0.15, kIntHi = 0.70;

constexpr Color_t kColPT3   = kBlack;
constexpr Color_t kColPT2   = kRed   + 1;
constexpr Color_t kColRatio = kBlue  + 1;

// Native MID free quad colour & style (reference).
const Color_t kColMid       = kAzure  + 1;   // DASHED

// Build the EXP 58-bin variable-width edges.
//   [0.00, 0.05] : 20 bins of 2.5 MeV       (DENSE LOW)
//   [0.05, 1.40] : 38 bins of ~35.526 MeV   (UNIFORM)
const double* getEdges() {
    static double edges[kNb + 1];
    static bool   built = false;
    if (!built) {
        int k = 0;
        for (int i = 0; i <= 20; ++i) edges[k++] = i * 0.0025;
        const double w = (1.4 - 0.05) / 38.0;
        for (int i = 1; i <= 38; ++i) edges[k++] = 0.05 + i * w;
        built = true;
    }
    return edges;
}

// Fill a histogram from a TTree using TTree::Draw with cut. Sumw2 first.
TH1D* drawWithCut(TTree* t, const char* expr,
                  const char* cut, const std::string& name) {
    TH1D* h = new TH1D(name.c_str(), "", kNb, getEdges());
    h->Sumw2();
    t->Draw((std::string(expr) + ">>" + name).c_str(), cut, "goff");
    h->SetDirectory(nullptr);
    return h;
}

// Compute CB = 2*sqrt(N_++ · N_--) per bin with proper error propagation.
TH1D* makeCB(TH1D* h_pp, TH1D* h_mm, const std::string& name) {
    TH1D* cb = static_cast<TH1D*>(h_pp->Clone(name.c_str()));
    cb->SetDirectory(nullptr);
    cb->Reset();
    for (int b = 1; b <= h_pp->GetNbinsX(); ++b) {
        const double n_pp = h_pp->GetBinContent(b);
        const double e_pp = h_pp->GetBinError(b);
        const double n_mm = h_mm->GetBinContent(b);
        const double e_mm = h_mm->GetBinError(b);
        if (n_pp <= 0.0 || n_mm <= 0.0 ||
            !std::isfinite(n_pp) || !std::isfinite(n_mm)) {
            cb->SetBinContent(b, 0.0);
            cb->SetBinError  (b, 0.0);
            continue;
        }
        const double val = 2.0 * std::sqrt(n_pp * n_mm);
        const double rel_pp = e_pp / n_pp;
        const double rel_mm = e_mm / n_mm;
        const double err    = val * 0.5 *
            std::sqrt(rel_pp * rel_pp + rel_mm * rel_mm);
        cb->SetBinContent(b, val);
        cb->SetBinError  (b, err);
    }
    return cb;
}

// Signal = epem − CB with quadrature errors.
TH1D* makeSignal(TH1D* h_em, TH1D* h_cb, const std::string& name) {
    TH1D* sig = static_cast<TH1D*>(h_em->Clone(name.c_str()));
    sig->SetDirectory(nullptr);
    sig->Reset();
    for (int b = 1; b <= h_em->GetNbinsX(); ++b) {
        const double n_em = h_em->GetBinContent(b);
        const double e_em = h_em->GetBinError(b);
        const double n_cb = h_cb->GetBinContent(b);
        const double e_cb = h_cb->GetBinError(b);
        const double val  = n_em - n_cb;
        const double err  = std::sqrt(e_em * e_em + e_cb * e_cb);
        sig->SetBinContent(b, val);
        sig->SetBinError  (b, err);
    }
    return sig;
}

// Bin-by-bin ratio = factor · num / den with relative-error propagation.
TH1D* makeRatio(TH1D* num, TH1D* den, double factor,
                const std::string& name) {
    TH1D* r = static_cast<TH1D*>(num->Clone(name.c_str()));
    r->SetDirectory(nullptr);
    r->Reset();
    for (int b = 1; b <= num->GetNbinsX(); ++b) {
        const double n  = num->GetBinContent(b);
        const double en = num->GetBinError(b);
        const double d  = den->GetBinContent(b);
        const double ed = den->GetBinError(b);
        if (d == 0.0 || !std::isfinite(d) || !std::isfinite(n)) continue;
        const double val = factor * n / d;
        const double rel = std::sqrt(
            (n != 0.0 ? std::pow(en / n, 2) : 0.0) +
            std::pow(ed / d, 2));
        r->SetBinContent(b, val);
        r->SetBinError  (b, std::abs(val) * rel);
    }
    return r;
}

// Divide each bin (content & error) by its bin width — spectra display.
void normalizeByBinWidth(TH1D* h) {
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
        const double w = h->GetBinWidth(b);
        if (w > 0.0) {
            h->SetBinContent(b, h->GetBinContent(b) / w);
            h->SetBinError  (b, h->GetBinError(b)   / w);
        }
    }
}

void styleHist(TH1D* h, Color_t c, Style_t m) {
    h->SetLineColor(c);
    h->SetMarkerColor(c);
    h->SetMarkerStyle(m);
    h->SetMarkerSize(0.9);
    h->SetLineWidth(2);
    h->GetXaxis()->SetTitleSize(0.050);
    h->GetYaxis()->SetTitleSize(0.050);
    h->GetXaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetLabelSize(0.045);
    h->GetXaxis()->SetTitleOffset(1.05);
    h->GetYaxis()->SetTitleOffset(1.20);
}

// Free quadratic [0]+[1]*x+[2]*x*x on [xlo, xhi] (reference fit only).
struct QuadResult {
    double a, a_err;
    double b, b_err;
    double c, c_err;
    double chi2_ndf;
    int    ndf;
    TF1*   f;
};

QuadResult fitQuadFree(TH1D* h, double xlo, double xhi,
                       Color_t color, Style_t style, const char* tag) {
    const std::string fname = std::string(h->GetName()) + "_fq_" + tag;
    TF1* f = new TF1(fname.c_str(), "[0]+[1]*x+[2]*x*x", xlo, xhi);
    f->SetLineColor(color);
    f->SetLineWidth(3);
    f->SetLineStyle(style);
    h->Fit(f, "RQ0");
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

// Integral helper: sum over hist bins fully contained in [xlo, xhi] of
// (bin content * bin width).
double histIntegralXBinWidth(TH1D* h, double xlo, double xhi) {
    double s = 0.0;
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
        const double lo = h->GetXaxis()->GetBinLowEdge(b);
        const double hi = h->GetXaxis()->GetBinUpEdge (b);
        if (lo < xlo - 1e-9) continue;
        if (hi > xhi + 1e-9) continue;
        s += h->GetBinContent(b) * h->GetBinWidth(b);
    }
    return s;
}

// ---------------------------------------------------------------------------
// Build the FULL sim ratio histogram on the cascade-v5 220-bin layout.
// Returns pointer to TH1D (caller takes ownership). Uses step5 (same_vertex)
// sim cut and sim_genweight. Bins with insufficient PT3 are left at 0.
// ---------------------------------------------------------------------------
TH1D* buildSimFullRatio(const std::string& tag) {
    TFile* f_sim = TFile::Open(kFileSim, "READ");
    if (!f_sim || f_sim->IsZombie()) {
        std::cerr << "[sim] Cannot open sim file: " << kFileSim << "\n";
        if (f_sim) f_sim->Close();
        return nullptr;
    }
    TTree* t_sim = dynamic_cast<TTree*>(f_sim->Get(kTreeSimRec));
    if (!t_sim) {
        std::cerr << "[sim] Tree '" << kTreeSimRec
                  << "' missing in sim file\n";
        f_sim->Close();
        return nullptr;
    }

    const std::string nm3 = "h_sim_full_PT3_" + tag;
    const std::string nm2 = "h_sim_full_PT2_" + tag;
    TH1D* h_sim_PT3 = new TH1D(nm3.c_str(), "", kSimNb, getSimEdges());
    TH1D* h_sim_PT2 = new TH1D(nm2.c_str(), "", kSimNb, getSimEdges());
    h_sim_PT3->Sumw2();
    h_sim_PT2->Sumw2();

    const std::string w_PT3 =
        std::string("(pt3==1)*") + kSimStep5Cut + "*sim_genweight";
    const std::string w_PT2 =
        std::string("(pt2==1)*") + kSimStep5Cut + "*sim_genweight";
    t_sim->Draw(("m_ee>>" + nm3).c_str(), w_PT3.c_str(), "goff");
    t_sim->Draw(("m_ee>>" + nm2).c_str(), w_PT2.c_str(), "goff");
    h_sim_PT3->SetDirectory(nullptr);
    h_sim_PT2->SetDirectory(nullptr);

    TH1D* r = static_cast<TH1D*>(
        h_sim_PT3->Clone(("h_sim_full_ratio_" + tag).c_str()));
    r->SetDirectory(nullptr);
    r->Reset();
    int n_drawn = 0;
    for (int b = 1; b <= kSimNb; ++b) {
        const double n2 = h_sim_PT2->GetBinContent(b);
        const double e2 = h_sim_PT2->GetBinError  (b);
        const double n3 = h_sim_PT3->GetBinContent(b);
        const double e3 = h_sim_PT3->GetBinError  (b);
        if (n3 < kSimMinDen || !std::isfinite(n3) || !std::isfinite(n2)) {
            r->SetBinContent(b, 0.0);
            r->SetBinError  (b, 0.0);
            continue;
        }
        const double val = n2 / n3;
        const double rel = std::sqrt(
            (n2 != 0.0 ? std::pow(e2 / n2, 2) : 0.0) +
            std::pow(e3 / n3, 2));
        r->SetBinContent(b, val);
        r->SetBinError  (b, std::abs(val) * rel);
        ++n_drawn;
    }
    printf("[sim] FULL 220 bins drawn (PT3>=%.0f): %d / %d\n",
           kSimMinDen, n_drawn, kSimNb);

    delete h_sim_PT3;
    delete h_sim_PT2;
    f_sim->Close();
    return r;
}

// Nearest-neighbour lookup: get sim ratio at mass x from the sim ratio TH1.
// If the corresponding sim bin has 0 content (undefined), walk outward to
// find the nearest defined bin.
double simRatioAt(TH1D* h_sim_ratio, double x) {
    if (!h_sim_ratio) return 0.0;
    int b = h_sim_ratio->FindBin(x);
    if (b < 1) b = 1;
    if (b > h_sim_ratio->GetNbinsX()) b = h_sim_ratio->GetNbinsX();
    if (h_sim_ratio->GetBinContent(b) > 0.0) {
        return h_sim_ratio->GetBinContent(b);
    }
    // Walk outward to find the nearest defined bin.
    const int nb = h_sim_ratio->GetNbinsX();
    for (int d = 1; d <= nb; ++d) {
        const int bl = b - d;
        const int bu = b + d;
        if (bl >= 1 && h_sim_ratio->GetBinContent(bl) > 0.0) {
            return h_sim_ratio->GetBinContent(bl);
        }
        if (bu <= nb && h_sim_ratio->GetBinContent(bu) > 0.0) {
            return h_sim_ratio->GetBinContent(bu);
        }
    }
    return 0.0;
}

// ---------------------------------------------------------------------------
// Produce one full output file (PDF + PNG) for the given numerator/denominator
// histograms and the corresponding display spectra.
// ---------------------------------------------------------------------------
void makeOnePlot(TH1D* h_num, TH1D* h_den,
                 TH1D* h_disp_PT3, TH1D* h_disp_PT2,
                 const char* left_title,
                 const char* right_title,
                 const char* right_ylabel,
                 const std::string& base_path,
                 const std::string& tag,
                 const std::string& variant_label) {

    std::cout << "\n========== " << variant_label << " ==========\n";

    // ---------- Ratio = 63 · num / den ----------
    TH1D* r_ratio = makeRatio(h_num, h_den, kTrigCorr,
                              "r_ratio_" + tag);
    styleHist(r_ratio, kColRatio, 20);

    // ---------- Native MID free quad (reference only, dashed) ----------
    QuadResult mid = fitQuadFree(r_ratio, kMdLo, kMdHi,
                                 kColMid, 2 /* dashed */,
                                 (tag + "_mid").c_str());

    // ---------- Build sim ratio (FULL 220-bin layout) ----------
    TH1D* h_sim_ratio = buildSimFullRatio(tag);
    if (!h_sim_ratio) {
        std::cerr << "Sim ratio could not be built; aborting.\n";
        return;
    }

    // ---------- INTEGRAL-MATCHED k_simpts over [kIntLo, kIntHi] ----------
    const double I_expRatio_mid =
        histIntegralXBinWidth(r_ratio, kIntLo, kIntHi);
    const double I_simpoints_mid =
        histIntegralXBinWidth(h_sim_ratio, kIntLo, kIntHi);
    const double k_simpts = (I_simpoints_mid != 0.0)
        ? I_expRatio_mid / I_simpoints_mid
        : 0.0;

    printf("\nINTEGRAL-MATCHED k_simpts (range [%.2f, %.2f]):\n",
           kIntLo, kIntHi);
    printf("  I_expRatio  (sum ratio(b) x binwidth(b))   = %.6f\n",
           I_expRatio_mid);
    printf("  I_simpoints (sum sim_ratio(b) x binwidth)  = %.6f\n",
           I_simpoints_mid);
    printf("  k_simpts    = I_expRatio / I_simpoints     = %.4f\n",
           k_simpts);

    // ---------- Build SCALED sim ratio (clone, * k_simpts) ----------
    TH1D* h_sim_ratio_scaled = static_cast<TH1D*>(
        h_sim_ratio->Clone(("h_sim_ratio_scaled_" + tag).c_str()));
    h_sim_ratio_scaled->SetDirectory(nullptr);
    h_sim_ratio_scaled->Scale(k_simpts);

    // ---------- Diagnostic yields ----------
    const double y_num_full = h_num->Integral();
    const double y_den_full = h_den->Integral();
    printf("Full-range yields:\n");
    printf("  PT3 (entering ratio)         : %.3f\n", y_den_full);
    printf("  PT2 (entering ratio, raw)    : %.3f\n", y_num_full);

    printf("\nNative MID FREE quad [%.2f, %.2f] (reference):\n",
           kMdLo, kMdHi);
    printf("  a = %+.4f +- %.4f\n", mid.a, mid.a_err);
    printf("  b = %+.4f +- %.4f\n", mid.b, mid.b_err);
    printf("  c = %+.4f +- %.4f\n", mid.c, mid.c_err);
    printf("  chi2/ndf = %.2f / %d = %.3f\n",
           mid.chi2_ndf * mid.ndf, mid.ndf, mid.chi2_ndf);

    // ---------- Build LEFT panel display spectra ----------
    TH1D* h_PT3_w = static_cast<TH1D*>(
        h_disp_PT3->Clone(("h_PT3_w_" + tag).c_str()));
    h_PT3_w->SetDirectory(nullptr);
    TH1D* h_PT2_w = static_cast<TH1D*>(
        h_disp_PT2->Clone(("h_PT2_w_x63_" + tag).c_str()));
    h_PT2_w->SetDirectory(nullptr);
    h_PT2_w->Scale(kTrigCorr);
    normalizeByBinWidth(h_PT3_w);
    normalizeByBinWidth(h_PT2_w);
    styleHist(h_PT3_w, kColPT3, 20);
    styleHist(h_PT2_w, kColPT2, 21);

    // ---------- Hypothetical 63 * PT2 = PT3(x) * k_simpts * sim_ratio_at(x) ----------
    // FULL mass range — no LOW skip.
    TH1D* h_PT2_hypo_w = static_cast<TH1D*>(
        h_disp_PT3->Clone(("h_PT2_hypo_w_" + tag).c_str()));
    h_PT2_hypo_w->SetDirectory(nullptr);
    h_PT2_hypo_w->Reset();
    for (int b = 1; b <= h_disp_PT3->GetNbinsX(); ++b) {
        const double n3  = h_disp_PT3->GetBinContent(b);
        const double e3  = h_disp_PT3->GetBinError  (b);
        if (n3 <= 0.0 || !std::isfinite(n3)) continue;
        const double xc  = h_disp_PT3->GetXaxis()->GetBinCenter(b);
        const double sr  = simRatioAt(h_sim_ratio, xc);
        const double sf  = k_simpts * sr;
        const double val = sf * n3;
        const double err = std::abs(sf) * e3;
        h_PT2_hypo_w->SetBinContent(b, val);
        h_PT2_hypo_w->SetBinError  (b, err);
    }
    normalizeByBinWidth(h_PT2_hypo_w);
    h_PT2_hypo_w->SetMarkerStyle(25);          // open square
    h_PT2_hypo_w->SetMarkerSize(1.2);
    h_PT2_hypo_w->SetMarkerColor(kBlue + 1);
    h_PT2_hypo_w->SetLineColor  (kBlue + 1);
    h_PT2_hypo_w->SetLineWidth  (1);

    // ---------- Canvas: 1×2 ----------
    TCanvas* c = new TCanvas(("c_" + tag).c_str(), variant_label.c_str(),
                             1700, 600);
    c->Divide(2, 1, 0.001, 0.001);

    // ----- Left: spectra -----
    c->cd(1);
    gPad->SetLogy(true);
    gPad->SetMargin(0.13, 0.04, 0.13, 0.10);

    double y_max_pt3  = h_PT3_w->GetMaximum();
    double y_max_pt2  = h_PT2_w->GetMaximum();
    double y_max_glob = std::max(y_max_pt3, y_max_pt2);
    double y_min_pos  = 1e300;
    for (int b = 1; b <= h_PT3_w->GetNbinsX(); ++b) {
        const double v3 = h_PT3_w->GetBinContent(b);
        if (v3 > 0.0 && v3 < y_min_pos) y_min_pos = v3;
        const double v2 = h_PT2_w->GetBinContent(b);
        if (v2 > 0.0 && v2 < y_min_pos) y_min_pos = v2;
    }
    if (!std::isfinite(y_min_pos) || y_min_pos > 1e299) y_min_pos = 0.1;

    const double y_max = 1.5  * y_max_glob;
    const double y_min = std::max(0.1, 0.3 * y_min_pos);

    h_PT3_w->SetTitle(
        Form("%s;M_{e^{+}e^{-}} [GeV/c^{2}];weighted entries / GeV/c^{2}",
             left_title));
    h_PT3_w->GetXaxis()->SetRangeUser(kXmin, kXmax);
    h_PT3_w->GetYaxis()->SetRangeUser(y_min, y_max);
    h_PT3_w->Draw("E1");
    h_PT2_w->Draw("E1 SAME");
    h_PT2_hypo_w->Draw("E1 SAME");

    TLegend* leg = new TLegend(0.40, 0.68, 0.95, 0.88);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.028);
    leg->AddEntry(h_PT3_w,
        Form("PT3 (#sum = %.0f)", h_disp_PT3->Integral()), "lpe");
    leg->AddEntry(h_PT2_w,
        Form("63 #upoint PT2 (red dots) - measured (#sum_{PT2} = %.0f)",
             h_disp_PT2->Integral()), "lpe");
    leg->AddEntry(h_PT2_hypo_w,
        "63 #upoint PT2 hypothetical (blue open sq., full range) "
        "= PT3 #times k_{simpts} #upoint sim_ratio(x)", "lpe");
    leg->Draw();

    // ----- Right: ratio + sim points overlay -----
    c->cd(2);
    gPad->SetMargin(0.13, 0.04, 0.13, 0.10);
    r_ratio->SetTitle(
        Form("%s;M_{e^{+}e^{-}} [GeV/c^{2}];%s", right_title, right_ylabel));
    r_ratio->GetXaxis()->SetRangeUser(kXmin, kXmax);
    r_ratio->GetYaxis()->SetRangeUser(0.5, 4.0);

    r_ratio->Draw("E1");

    TLine* lref = new TLine(kXmin, 1.0, kXmax, 1.0);
    lref->SetLineStyle(2); lref->SetLineColor(kGray + 2);
    lref->Draw();

    // ----- Style and draw the SCALED SIM POINTS overlay (full 220-bin) -----
    // Gray-filled triangles to be visible as a dense cloud but not drown out
    // the exp markers.
    h_sim_ratio_scaled->SetMarkerStyle(22);   // filled triangle up
    h_sim_ratio_scaled->SetMarkerSize(0.5);
    h_sim_ratio_scaled->SetMarkerColor(kGray + 2);
    h_sim_ratio_scaled->SetLineColor  (kGray + 1);
    h_sim_ratio_scaled->SetLineWidth  (1);
    h_sim_ratio_scaled->Draw("P SAME");

    // ----- LOW [0, 0.15] separate red open-circle annotation -----
    // Build by cloning the sim ratio (already scaled) and zeroing bins
    // outside the LOW range.
    TH1D* h_sim_low_scaled = static_cast<TH1D*>(
        h_sim_ratio_scaled->Clone(("h_sim_low_scaled_" + tag).c_str()));
    h_sim_low_scaled->SetDirectory(nullptr);
    for (int b = 1; b <= h_sim_low_scaled->GetNbinsX(); ++b) {
        const double xc = h_sim_low_scaled->GetXaxis()->GetBinCenter(b);
        if (xc > 0.15 + 1e-9) {
            h_sim_low_scaled->SetBinContent(b, 0.0);
            h_sim_low_scaled->SetBinError  (b, 0.0);
        }
    }
    h_sim_low_scaled->SetMarkerStyle(24);   // open circle
    h_sim_low_scaled->SetMarkerSize(1.0);
    h_sim_low_scaled->SetMarkerColor(kRed);
    h_sim_low_scaled->SetLineColor  (kRed);
    h_sim_low_scaled->SetLineWidth  (1);
    h_sim_low_scaled->Draw("P SAME");

    // Redraw exp markers on top so they're not buried.
    r_ratio->Draw("E1 SAME");

    // Native MID quad reference (dashed, kept for reference).
    if (mid.f) mid.f->Draw("SAME");

    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.026);

    tx.SetTextColor(kColMid);
    tx.DrawLatex(0.16, 0.88,
        Form("[%.2f-%.2f] native quad FREE (reference): a=%.3f b=%.3f c=%.3f  "
             "#chi^{2}/ndf = %.1f / %d",
             kMdLo, kMdHi,
             mid.a, mid.b, mid.c,
             mid.chi2_ndf * mid.ndf, mid.ndf));

    tx.SetTextColor(kGray + 3);
    tx.DrawLatex(0.16, 0.84,
        Form("step5 POINTS (no fit), k_{simpts} = %.3f  "
             "(norm: integral [%.2f, %.2f])",
             k_simpts, kIntLo, kIntHi));

    tx.SetTextColor(kGray + 3);
    tx.DrawLatex(0.16, 0.80,
        Form("I_{exp}^{ratio} = %.4f,  I_{sim}^{points} = %.4f",
             I_expRatio_mid, I_simpoints_mid));

    tx.SetTextColor(kRed);
    tx.DrawLatex(0.16, 0.76,
        Form("LOW [0, 0.15] sim points (open red circles) = "
             "k_{simpts} #upoint sim PT2/PT3 (step5)"));

    TLegend* leg_r = new TLegend(0.55, 0.62, 0.95, 0.78);
    leg_r->SetBorderSize(0); leg_r->SetFillStyle(0); leg_r->SetTextSize(0.026);
    leg_r->AddEntry(mid.f,
        "native MID free quad [0.10-0.80] (ref.)", "l");
    leg_r->AddEntry(h_sim_ratio_scaled,
        Form("scaled sim POINTS (step5, k_{simpts}=%.3f)", k_simpts), "pe");
    leg_r->AddEntry(h_sim_low_scaled,
        "scaled sim POINTS in LOW [0, 0.15]", "pe");
    leg_r->Draw();

    // ---------- High-mass bin comparison ----------
    printf("\nHigh-mass comparison (raw bin contents, NOT width-normalized):\n");
    printf("    x_target   bin[lo, hi]            PT3        63*PT2_actual   "
           "63*PT2_hypo  hypo/actual\n");
    const double sample_xs[] = {1.03, 1.06, 1.10, 1.20};
    for (double xs : sample_xs) {
        const int b = h_disp_PT3->FindBin(xs);
        const double lo = h_disp_PT3->GetXaxis()->GetBinLowEdge(b);
        const double hi = h_disp_PT3->GetXaxis()->GetBinUpEdge (b);
        const double xc = h_disp_PT3->GetXaxis()->GetBinCenter (b);
        const double n3 = h_disp_PT3->GetBinContent(b);
        const double n2 = h_disp_PT2->GetBinContent(b) * kTrigCorr;
        const double sr = simRatioAt(h_sim_ratio, xc);
        const double sf = k_simpts * sr;
        const double n2_hypo = sf * n3;
        const double rr =
            (std::abs(n2) > 1e-12) ? n2_hypo / n2 : 0.0;
        printf("    x=%.2f   [%.4f, %.4f]   %8.3f    %10.3f     %10.3f    %.3f\n",
               xs, lo, hi, n3, n2, n2_hypo, rr);
    }

    // ---------- Save ----------
    c->SaveAs((base_path + ".pdf").c_str());
    c->SaveAs((base_path + ".png").c_str());
    std::cout << "Saved: " << base_path << ".{pdf,png}\n";
}

}  // anonymous namespace

void mass_ee_pt_ratio_signal_varbins_exp_step5pts() {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetErrorX(0.0);

    // ---------- Open input files ----------
    TFile* f_em = TFile::Open(kFileEm, "READ");
    TFile* f_pp = TFile::Open(kFilePp, "READ");
    TFile* f_mm = TFile::Open(kFileMm, "READ");
    if (!f_em || f_em->IsZombie() ||
        !f_pp || f_pp->IsZombie() ||
        !f_mm || f_mm->IsZombie()) {
        std::cerr << "Cannot open one of the EXP input files.\n";
        return;
    }
    TTree* t_em = dynamic_cast<TTree*>(f_em->Get(kTree));
    TTree* t_pp = dynamic_cast<TTree*>(f_pp->Get(kTree));
    TTree* t_mm = dynamic_cast<TTree*>(f_mm->Get(kTree));
    if (!t_em || !t_pp || !t_mm) {
        std::cerr << "Tree " << kTree << " missing in one of the EXP files.\n";
        return;
    }
    std::cout << "Input:\n";
    std::cout << "  epem : " << kFileEm << " (" << t_em->GetEntries() << " entries)\n";
    std::cout << "  epep : " << kFilePp << " (" << t_pp->GetEntries() << " entries)\n";
    std::cout << "  emem : " << kFileMm << " (" << t_mm->GetEntries() << " entries)\n";
    std::cout << "Binning (EXP variable, " << kNb << " bins):\n";
    std::cout << "  [0.00, 0.05] : 20 x 2.5 MeV (DENSE)\n";
    std::cout << "  [0.05, 1.40] : 38 x ~35.526 MeV (UNIFORM)\n";
    std::cout << "Binning (SIM variable, " << kSimNb
              << " bins, cascade v5):\n";
    std::cout << "  [0.00, 0.15] :  60 x 2.5 MeV\n";
    std::cout << "  [0.15, 0.85] : 140 x 5 MeV\n";
    std::cout << "  [0.85, 1.10] :  12 x (0.25/12) GeV\n";
    std::cout << "  [1.10, 1.40] :   8 x (0.30/8) GeV\n";
    std::cout << "Trigger cuts: PT3='" << kCutPT3
              << "'  PT2='" << kCutPT2 << "'\n";
    std::cout << "Correction factor (PT2 downscale 64 -> factor): "
              << kTrigCorr << "\n\n";

    gSystem->mkdir("plots/output", true);

    // ---------- Build raw mass histograms per trigger ----------
    TH1D* h_em_PT3 = drawWithCut(t_em, "m_ee", kCutPT3, "h_em_PT3");
    TH1D* h_em_PT2 = drawWithCut(t_em, "m_ee", kCutPT2, "h_em_PT2");
    TH1D* h_pp_PT3 = drawWithCut(t_pp, "m_ee", kCutPT3, "h_pp_PT3");
    TH1D* h_pp_PT2 = drawWithCut(t_pp, "m_ee", kCutPT2, "h_pp_PT2");
    TH1D* h_mm_PT3 = drawWithCut(t_mm, "m_ee", kCutPT3, "h_mm_PT3");
    TH1D* h_mm_PT2 = drawWithCut(t_mm, "m_ee", kCutPT2, "h_mm_PT2");

    std::cout << "Raw yields:\n";
    std::cout << "  PT3: epem=" << h_em_PT3->Integral()
              << "  ++="       << h_pp_PT3->Integral()
              << "  --="       << h_mm_PT3->Integral() << "\n";
    std::cout << "  PT2: epem=" << h_em_PT2->Integral()
              << "  ++="       << h_pp_PT2->Integral()
              << "  --="       << h_mm_PT2->Integral() << "\n";

    // ---------- CB per trigger, signal per trigger ----------
    TH1D* h_cb_PT3  = makeCB(h_pp_PT3, h_mm_PT3, "h_cb_PT3");
    TH1D* h_cb_PT2  = makeCB(h_pp_PT2, h_mm_PT2, "h_cb_PT2");
    TH1D* h_sig_PT3 = makeSignal(h_em_PT3, h_cb_PT3, "h_sig_PT3");
    TH1D* h_sig_PT2 = makeSignal(h_em_PT2, h_cb_PT2, "h_sig_PT2");

    std::cout << "Signal yields (epem - CB):\n";
    std::cout << "  PT3 signal = " << h_sig_PT3->Integral() << "\n";
    std::cout << "  PT2 signal = " << h_sig_PT2->Integral() << "\n";

    // =====================================================================
    // OUTPUT: SIGNAL (CB-subtracted) — single output, single PDF/PNG
    // =====================================================================
    makeOnePlot(
        h_sig_PT2, h_sig_PT3,
        h_sig_PT3, h_sig_PT2,
        "m_{ee} CB-subtracted signal (PT3 black; 63 #upoint PT2 red)",
        "63 #upoint N_{PT2}^{sig} / N_{PT3}^{sig} (exp, step5 POINTS)",
        "63 #upoint N_{PT2}^{sig} / N_{PT3}^{sig}",
        "plots/output/mass_ee_pt_ratio_varbins_exp_signal_step5pts",
        "signal_step5pts",
        "SIGNAL (CB-subtracted, step5 same_vertex POINTS)");

    f_em->Close();
    f_pp->Close();
    f_mm->Close();
}
