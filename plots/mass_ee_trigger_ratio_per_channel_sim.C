// mass_ee_trigger_ratio_per_channel_sim.C — H1 diagnostic for the sinusoidal
// shape observed in the SIM trigger-correction ratio N_PT2 / N_PT3.
//
// ───────────────────────────────────────────────────────────────────────────
//  HYPOTHESIS UNDER TEST  (H1: SMASH process-mix mismatch)
// ───────────────────────────────────────────────────────────────────────────
// The total SIM ratio N_PT2 / N_PT3 vs m_ee (truth) shows a clear sinus-like
// modulation: a minimum near m_ee ≈ 0.5 GeV/c² and a maximum near
// m_ee ≈ 1.2 GeV/c². If this modulation is intrinsic to ANY single physics
// channel (kinematic threshold, OA cut, conversion contamination, smearing —
// the H2/H3/H4 hypotheses) then each per-channel ratio must also wobble.
//
// Conversely, if every individual channel is FLAT (low χ²/ndf for a constant
// fit) but the total still wobbles, the modulation must come from the
// CHANNEL-WEIGHT MIX changing across m_ee (different channels dominate
// different m_ee regions, and they sit at different absolute ratio levels).
// That is the signature of H1 — a SMASH process-mix mismatch.
//
// We use m_ee_sim (TRUTH mass) on dilepton_nt_cor for this test because:
//   • The hard MC purity gate in main.cc guarantees all 286M events are
//     clean truth e⁺e⁻ pairs from the same decay vertex.
//   • Using truth m_ee removes reconstruction/digitization smearing as a
//     confound (so it isolates the H1 question from H3).
//
// ───────────────────────────────────────────────────────────────────────────
//  LAYOUT  (1 row × 2 columns, 1700 × 600 px)
// ───────────────────────────────────────────────────────────────────────────
//   Col 1 (log Y): m_ee_sim spectra @ PT3==1 — Total (black HIST) overlaid
//                  with π⁰ Dal., η Dal., ω all per the channel colour code.
//   Col 2 (linear Y, range [0, 3]): N_PT2 / N_PT3 ratios overlaid with
//                  matching-colour constant fits, dashed grey ref line at 1.0.
//                  Multi-line TLatex prints  a ± err   χ²/ndf  per curve;
//                  comparing per-channel χ²/ndf vs total χ²/ndf IS the H1 test.
//
// Channels (ep_sim_geninfo1) — IDs verified by m_ee peak / width matching
// the PLUTO / HGeant particle numbering used by the SMASH HADES output:
//    7051   π⁰ Dalitz       (~92.8 % weighted)   — peak m_ee≈0, endpoint m_π⁰=0.135
//   17051   η  Dalitz       (~ 6.0 %)            — peak m_ee≈0, endpoint m_η=0.548
//      41   ρ⁰ → e⁺e⁻        (~0.03 %)            — broad Breit-Wigner, peak 0.75 GeV, Γ~240 MeV
//      52   ω → e⁺e⁻ direct  (~0.05 %)            — narrow peak at m_ω=0.7826 (smeared)
//   52051   ω → π⁰ e⁺e⁻ Dal. (~0.22 %)            — broad, endpoint = m_ω-m_π⁰≈0.65 (smeared)
//      55   φ → e⁺e⁻ direct  (~0.006 %)           — narrow peak at m_φ=1.0195
//
// Earlier (incorrect) version of this macro labelled channel 41 as "ω all".
// That was wrong — 41 is the ρ⁰ direct decay; ω is split between 52 (direct)
// and 52051 (Dalitz via π⁰).  The empirical m_ee shape (broad Breit-Wigner for
// 41, sharp peak at m_ω for 52, Dalitz continuum to ~0.65 for 52051) confirms
// the PLUTO/HGeant convention.
//
// Binning: 70 bins on [0, 1.4] GeV/c² (= 20 MeV bins).
// Default constant-fit range: [kFitLo=0.1, kFitHi=1.4] GeV/c².
//
// Usage:
//   root -l -b -q plots/mass_ee_trigger_ratio_per_channel_sim.C
//   root -l -b -q 'plots/mass_ee_trigger_ratio_per_channel_sim.C(0.1, 0.8)'
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
#include <limits>
#include <string>

namespace {

constexpr const char* kInputFile = "output_epem_sim.root";
constexpr const char* kNtName    = "dilepton_nt_cor";  // truth m_ee_sim lives here

// 70 bins on [0, 1.4] GeV/c² = 20 MeV per bin (per task spec).
constexpr int    kNb    = 70;
constexpr double kXmin  = 0.0;
constexpr double kXmax  = 1.4;

// Default constant-fit range [kFitLo, kFitHi] (overridable by macro args).
constexpr double kFitLo = 0.1;
constexpr double kFitHi = 1.4;

// ── Colour code ───────────────────────────────────────────────────────────
// Related ω-channels share the kAzure family (direct = brighter, Dalitz = darker dashed)
// so the reader visually groups them. ρ⁰ gets distinct orange, φ distinct red.
constexpr Color_t kColTot         = kBlack;
constexpr Color_t kColTotFit      = kGray + 2;
constexpr Color_t kColPi0         = kGreen   + 2;
constexpr Color_t kColPi0Fit      = kGreen   + 3;
constexpr Color_t kColEta         = kMagenta + 1;
constexpr Color_t kColEtaFit      = kMagenta + 2;
constexpr Color_t kColRho         = kOrange  + 7;
constexpr Color_t kColRhoFit      = kOrange  + 8;
constexpr Color_t kColOmegaDir    = kAzure   + 1;
constexpr Color_t kColOmegaDirFit = kAzure   + 2;
constexpr Color_t kColOmegaDal    = kAzure   - 3;
constexpr Color_t kColOmegaDalFit = kAzure   - 4;
constexpr Color_t kColPhi         = kRed     + 2;
constexpr Color_t kColPhiFit      = kRed     + 3;

// ─────────────────────────────────────────────────────────────────────────────
// Helpers (mirrors patterns in mass_ee_ratio_const_fit_*.C)
// ─────────────────────────────────────────────────────────────────────────────

// Decide whether sim_genweight is non-trivially filled. Returns the
// TTree::Draw weight expression to multiply into per-event cuts.
std::string detectWeightExpr(TTree* nt) {
    TH1D* h = new TH1D("h_wsum_detect_perch", "", 1, 0, 2);
    nt->Draw("1>>h_wsum_detect_perch", "sim_genweight", "goff");
    const double sum = h->Integral();
    delete h;
    if (sum > 0.0) {
        std::cout << "[sim] sim_genweight integral = " << sum
                  << "  → using weighted fills (sim_genweight)\n";
        return "sim_genweight";
    }
    std::cout << "[sim] sim_genweight integral = 0  → unweighted (1.0)\n";
    return "1.0";
}

// Draw m_ee_sim from dilepton_nt_cor into a (kNb, kXmin, kXmax) TH1D.
TH1D* drawFromNt(TTree* t, const std::string& weight_expr,
                 const std::string& name) {
    TH1D* h = new TH1D(name.c_str(), "", kNb, kXmin, kXmax);
    h->Sumw2();
    t->Draw(("m_ee_sim>>" + name).c_str(), weight_expr.c_str(), "goff");
    h->SetDirectory(nullptr);
    return h;
}

// Per-bin error propagation:  r = num / den, σ_r = |r| · √((σn/n)² + (σd/d)²).
TH1D* makeRatio(TH1D* num, TH1D* den, const std::string& name,
                double scale = 1.0) {
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
            const double val = scale * n / d;
            const double rel = std::sqrt(
                (n != 0.0 ? std::pow(en / n, 2) : 0.0) +
                std::pow(ed / d, 2));
            r->SetBinContent(b, val);
            r->SetBinError  (b, std::abs(val) * rel);
        }
    }
    return r;
}

void styleDot(TH1D* h, Color_t color, Style_t marker, double size = 0.9) {
    h->SetMarkerStyle(marker);
    h->SetMarkerSize(size);
    h->SetMarkerColor(color);
    h->SetLineColor(color);
    h->GetXaxis()->SetTitleSize(0.050);
    h->GetYaxis()->SetTitleSize(0.050);
    h->GetXaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetLabelSize(0.045);
    h->GetXaxis()->SetTitleOffset(1.05);
    h->GetYaxis()->SetTitleOffset(1.25);
}

void styleLine(TH1D* h, Color_t color, int width = 2, int style = 1) {
    h->SetLineColor(color);
    h->SetLineWidth(width);
    h->SetLineStyle(style);
    h->SetFillStyle(0);
    h->GetXaxis()->SetTitleSize(0.050);
    h->GetYaxis()->SetTitleSize(0.050);
    h->GetXaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetLabelSize(0.045);
    h->GetXaxis()->SetTitleOffset(1.05);
    h->GetYaxis()->SetTitleOffset(1.25);
}

struct FitResult { double a; double e; double chi2_ndf; int ndf; };

// Constant fit y = a, drawn in `color`. Uses weighted χ² (option "RQ" — no W).
FitResult fitConst(TH1D* h, double xlo, double xhi, Color_t color) {
    const std::string fname = std::string(h->GetName()) + "_fconst_perch";
    TF1* f = new TF1(fname.c_str(), "[0]", xlo, xhi);
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

// Max bin value + error across a histogram (for auto-Y scaling).
double safeMax(TH1D* h) {
    double mx = 0.0;
    const int nb = h->GetNbinsX();
    for (int b = 1; b <= nb; ++b) {
        const double v = h->GetBinContent(b) + h->GetBinError(b);
        if (std::isfinite(v) && v > mx) mx = v;
    }
    return mx;
}

// Smallest POSITIVE bin content across a histogram (for log-Y floor).
double smallestPositive(TH1D* h) {
    double mn = std::numeric_limits<double>::infinity();
    const int nb = h->GetNbinsX();
    for (int b = 1; b <= nb; ++b) {
        const double v = h->GetBinContent(b);
        if (v > 0.0 && v < mn) mn = v;
    }
    return std::isfinite(mn) ? mn : 1e-3;
}

}  // anonymous namespace

// ─────────────────────────────────────────────────────────────────────────────
void mass_ee_trigger_ratio_per_channel_sim(double fit_xmin = -1.0,
                                           double fit_xmax = -1.0) {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetErrorX(0.0);   // no horizontal error bars on the ratio markers

    // Resolve fit range; -1 sentinel ⇒ use the default [kFitLo, kFitHi].
    const bool   custom_range =
        (fit_xmin >= 0.0 && fit_xmin < kXmax) ||
        (fit_xmax >  0.0 && fit_xmax <= kXmax);
    const double fxlo = (fit_xmin >= 0.0 && fit_xmin < kXmax) ? fit_xmin : kFitLo;
    const double fxhi = (fit_xmax >  0.0 && fit_xmax <= kXmax) ? fit_xmax : kFitHi;

    // -------- Open input --------
    TFile* fs = TFile::Open(kInputFile, "READ");
    if (!fs || fs->IsZombie()) {
        std::cerr << "Cannot open " << kInputFile << "\n";
        return;
    }
    TTree* nt = dynamic_cast<TTree*>(fs->Get(kNtName));
    if (!nt) {
        std::cerr << kNtName << " missing in " << kInputFile << "\n";
        return;
    }

    // -------- Compose weight expressions --------
    const std::string w_expr = detectWeightExpr(nt);

    // Total (no channel filter) and per-channel weight strings.
    const std::string w_pt3_tot  = "(pt3==1) * "                              + w_expr;
    const std::string w_pt2_tot  = "(pt2==1) * "                              + w_expr;
    const std::string w_pt3_pi0  = "(pt3==1) * (ep_sim_geninfo1==7051) * "    + w_expr;
    const std::string w_pt2_pi0  = "(pt2==1) * (ep_sim_geninfo1==7051) * "    + w_expr;
    const std::string w_pt3_eta  = "(pt3==1) * (ep_sim_geninfo1==17051) * "   + w_expr;
    const std::string w_pt2_eta  = "(pt2==1) * (ep_sim_geninfo1==17051) * "   + w_expr;
    const std::string w_pt3_rho  = "(pt3==1) * (ep_sim_geninfo1==41) * "      + w_expr;
    const std::string w_pt2_rho  = "(pt2==1) * (ep_sim_geninfo1==41) * "      + w_expr;
    const std::string w_pt3_omgD = "(pt3==1) * (ep_sim_geninfo1==52) * "      + w_expr;
    const std::string w_pt2_omgD = "(pt2==1) * (ep_sim_geninfo1==52) * "      + w_expr;
    const std::string w_pt3_omgL = "(pt3==1) * (ep_sim_geninfo1==52051) * "   + w_expr;
    const std::string w_pt2_omgL = "(pt2==1) * (ep_sim_geninfo1==52051) * "   + w_expr;
    const std::string w_pt3_phi  = "(pt3==1) * (ep_sim_geninfo1==55) * "      + w_expr;
    const std::string w_pt2_phi  = "(pt2==1) * (ep_sim_geninfo1==55) * "      + w_expr;

    std::cout << "Tree:   " << kNtName << "   (truth m_ee_sim)\n"
              << "Bins:   " << kNb << " on [" << kXmin << ", " << kXmax << "]"
              << "  (= " << (1000.0 * (kXmax - kXmin) / kNb) << " MeV / bin)\n"
              << "Fit range: [" << fxlo << ", " << fxhi << "] GeV/c^2"
              << (custom_range ? "  (custom)" : "  (default)") << "\n\n"
              << "PT3 weight strings:\n"
              << "  total       = " << w_pt3_tot  << "\n"
              << "  π⁰  Dalitz  = " << w_pt3_pi0  << "\n"
              << "  η   Dalitz  = " << w_pt3_eta  << "\n"
              << "  ρ⁰  (41)    = " << w_pt3_rho  << "\n"
              << "  ω   (52)    = " << w_pt3_omgD << "\n"
              << "  ω-D (52051) = " << w_pt3_omgL << "\n"
              << "  φ   (55)    = " << w_pt3_phi  << "\n";

    // -------- Draw m_ee_sim spectra (this is the slow step) --------
    std::cout << "\nDrawing PT3 / PT2 m_ee_sim spectra (Total + 6 channels)... "
              << std::flush;
    TH1D* h_PT3_tot  = drawFromNt(nt, w_pt3_tot,  "h_PT3_tot");
    TH1D* h_PT2_tot  = drawFromNt(nt, w_pt2_tot,  "h_PT2_tot");
    TH1D* h_PT3_pi0  = drawFromNt(nt, w_pt3_pi0,  "h_PT3_pi0");
    TH1D* h_PT2_pi0  = drawFromNt(nt, w_pt2_pi0,  "h_PT2_pi0");
    TH1D* h_PT3_eta  = drawFromNt(nt, w_pt3_eta,  "h_PT3_eta");
    TH1D* h_PT2_eta  = drawFromNt(nt, w_pt2_eta,  "h_PT2_eta");
    TH1D* h_PT3_rho  = drawFromNt(nt, w_pt3_rho,  "h_PT3_rho");
    TH1D* h_PT2_rho  = drawFromNt(nt, w_pt2_rho,  "h_PT2_rho");
    TH1D* h_PT3_omgD = drawFromNt(nt, w_pt3_omgD, "h_PT3_omgD");
    TH1D* h_PT2_omgD = drawFromNt(nt, w_pt2_omgD, "h_PT2_omgD");
    TH1D* h_PT3_omgL = drawFromNt(nt, w_pt3_omgL, "h_PT3_omgL");
    TH1D* h_PT2_omgL = drawFromNt(nt, w_pt2_omgL, "h_PT2_omgL");
    TH1D* h_PT3_phi  = drawFromNt(nt, w_pt3_phi,  "h_PT3_phi");
    TH1D* h_PT2_phi  = drawFromNt(nt, w_pt2_phi,  "h_PT2_phi");
    std::cout << "done.\n";

    // Channel fractions (weighted PT3 integrals over the histogram window).
    const double I_tot  = h_PT3_tot ->Integral();
    const double I_pi0  = h_PT3_pi0 ->Integral();
    const double I_eta  = h_PT3_eta ->Integral();
    const double I_rho  = h_PT3_rho ->Integral();
    const double I_omgD = h_PT3_omgD->Integral();
    const double I_omgL = h_PT3_omgL->Integral();
    const double I_phi  = h_PT3_phi ->Integral();
    auto pct = [&](double x) { return (I_tot > 0.0) ? 100.0 * x / I_tot : 0.0; };
    std::cout << "PT3-weighted integrals (m_ee_sim ∈ [0, 1.4]):\n"
              << "  total          = " << I_tot  << "\n"
              << "  π⁰  Dalitz     = " << I_pi0  << "   (" << pct(I_pi0)  << " %)\n"
              << "  η   Dalitz     = " << I_eta  << "   (" << pct(I_eta)  << " %)\n"
              << "  ρ⁰  direct (41)= " << I_rho  << "   (" << pct(I_rho)  << " %)\n"
              << "  ω   direct (52)= " << I_omgD << "   (" << pct(I_omgD) << " %)\n"
              << "  ω   Dalitz (52051)= " << I_omgL << "   (" << pct(I_omgL) << " %)\n"
              << "  φ   direct (55)= " << I_phi  << "   (" << pct(I_phi)  << " %)\n";

    // -------- Ratios N_PT2 / N_PT3 --------
    TH1D* r_tot  = makeRatio(h_PT2_tot,  h_PT3_tot,  "r_tot");
    TH1D* r_pi0  = makeRatio(h_PT2_pi0,  h_PT3_pi0,  "r_pi0");
    TH1D* r_eta  = makeRatio(h_PT2_eta,  h_PT3_eta,  "r_eta");
    TH1D* r_rho  = makeRatio(h_PT2_rho,  h_PT3_rho,  "r_rho");
    TH1D* r_omgD = makeRatio(h_PT2_omgD, h_PT3_omgD, "r_omgD");
    TH1D* r_omgL = makeRatio(h_PT2_omgL, h_PT3_omgL, "r_omgL");
    TH1D* r_phi  = makeRatio(h_PT2_phi,  h_PT3_phi,  "r_phi");

    // -------- Style --------
    // Spectra: total as thick black; channels coloured. ω-direct and ω-Dalitz
    // share the kAzure family with different line styles (solid vs dashed)
    // to highlight that they're parts of the same parent ω.
    styleLine(h_PT3_tot,  kColTot,      3, 1);   // bold black solid (total)
    styleLine(h_PT3_pi0,  kColPi0,      2, 1);
    styleLine(h_PT3_eta,  kColEta,      2, 1);
    styleLine(h_PT3_rho,  kColRho,      2, 1);
    styleLine(h_PT3_omgD, kColOmegaDir, 2, 1);   // solid for direct ω
    styleLine(h_PT3_omgL, kColOmegaDal, 2, 2);   // dashed for ω-Dalitz
    styleLine(h_PT3_phi,  kColPhi,      2, 1);

    // Ratios: per-channel markers with distinct shapes.
    styleDot(r_tot,  kColTot,      20, 0.9);   // filled circle
    styleDot(r_pi0,  kColPi0,      25, 1.1);   // open square
    styleDot(r_eta,  kColEta,      24, 1.1);   // open circle
    styleDot(r_rho,  kColRho,      26, 1.3);   // open triangle-up
    styleDot(r_omgD, kColOmegaDir, 22, 1.1);   // filled triangle-up
    styleDot(r_omgL, kColOmegaDal, 32, 1.3);   // open triangle-down
    styleDot(r_phi,  kColPhi,      29, 1.3);   // full star

    // Titles
    h_PT3_tot->SetTitle(
        ";M_{e^{+}e^{-}} [GeV/c^{2}];weighted entries / 20 MeV (PT3==1)");
    r_tot->SetTitle(
        ";M_{e^{+}e^{-}} [GeV/c^{2}];N_{PT2} / N_{PT3}");

    // -------- Fits (drawn-with-data; computed BEFORE the canvas Draw) --------
    std::cout << "\nConstant fit y = a on N_PT2 / N_PT3:\n";
    const FitResult f_tot  = fitConst(r_tot,  fxlo, fxhi, kColTotFit);
    const FitResult f_pi0  = fitConst(r_pi0,  fxlo, fxhi, kColPi0Fit);
    const FitResult f_eta  = fitConst(r_eta,  fxlo, fxhi, kColEtaFit);
    const FitResult f_rho  = fitConst(r_rho,  fxlo, fxhi, kColRhoFit);
    const FitResult f_omgD = fitConst(r_omgD, fxlo, fxhi, kColOmegaDirFit);
    const FitResult f_omgL = fitConst(r_omgL, fxlo, fxhi, kColOmegaDalFit);
    const FitResult f_phi  = fitConst(r_phi,  fxlo, fxhi, kColPhiFit);
    auto pf = [](const char* tag, const FitResult& f) {
        printf("  %-18s a = %.4f +/- %.4f   chi2/ndf = %10.1f / %d  =  %.2f\n",
               tag, f.a, f.e, f.chi2_ndf * f.ndf, f.ndf, f.chi2_ndf);
    };
    pf("total",            f_tot);
    pf("π⁰  Dalitz (7051)", f_pi0);
    pf("η   Dalitz (17051)",f_eta);
    pf("ρ⁰  direct (41)",   f_rho);
    pf("ω   direct (52)",   f_omgD);
    pf("ω   Dalitz (52051)",f_omgL);
    pf("φ   direct (55)",   f_phi);

    // ── H1 verdict (printed to stdout for the orchestrator to capture) ──
    const double c_tot = f_tot.chi2_ndf;
    const double c_max_ch = std::max({f_pi0.chi2_ndf,  f_eta.chi2_ndf,
                                      f_rho.chi2_ndf,  f_omgD.chi2_ndf,
                                      f_omgL.chi2_ndf, f_phi.chi2_ndf});
    std::cout << "\nH1 verdict (process-mix mismatch):\n"
              << "  total chi2/ndf      = " << c_tot << "\n"
              << "  max channel chi2/ndf = " << c_max_ch << "\n"
              << "  ratio (total / max_ch) = "
              << (c_max_ch > 0.0 ? c_tot / c_max_ch : 0.0) << "\n";

    // -------- Canvas --------
    TCanvas* c = new TCanvas("c_mass_ee_trigger_ratio_per_channel_sim",
                             "Trigger ratio per channel (H1 diagnostic)",
                             1700, 600);
    c->Divide(2, 1, 0.001, 0.001);

    // ── Col 1: m_ee_sim spectra @ PT3==1, log Y ──
    c->cd(1);
    gPad->SetMargin(0.13, 0.04, 0.13, 0.10);
    gPad->SetLogy();

    // Auto-Y for log scale.
    //   y_top = max bin content (with error) × 5  ← headroom for legend
    //   y_bot = smallest positive bin content across ALL curves × 0.5
    // The smallest-positive logic ensures rare channels (φ, ω all are
    // ~10⁴× smaller than total at most m_ee) stay visible.
    const double y_top = std::max({safeMax(h_PT3_tot),  safeMax(h_PT3_pi0),
                                   safeMax(h_PT3_eta),  safeMax(h_PT3_rho),
                                   safeMax(h_PT3_omgD), safeMax(h_PT3_omgL),
                                   safeMax(h_PT3_phi)});
    const double y_bot_raw = std::min({smallestPositive(h_PT3_tot),
                                       smallestPositive(h_PT3_pi0),
                                       smallestPositive(h_PT3_eta),
                                       smallestPositive(h_PT3_rho),
                                       smallestPositive(h_PT3_omgD),
                                       smallestPositive(h_PT3_omgL),
                                       smallestPositive(h_PT3_phi)});
    const double y_bot = std::max(1e-3, y_bot_raw * 0.5);
    h_PT3_tot ->GetYaxis()->SetRangeUser(y_bot, y_top * 5.0);
    h_PT3_tot ->Draw("HIST");
    h_PT3_pi0 ->Draw("HIST SAME");
    h_PT3_eta ->Draw("HIST SAME");
    h_PT3_rho ->Draw("HIST SAME");
    h_PT3_omgD->Draw("HIST SAME");
    h_PT3_omgL->Draw("HIST SAME");
    h_PT3_phi ->Draw("HIST SAME");

    TLegend* leg = new TLegend(0.50, 0.45, 0.95, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.029);
    leg->AddEntry(h_PT3_tot,  "total",                                                "l");
    leg->AddEntry(h_PT3_pi0,  Form("#pi^{0} Dalitz (7051):  %.2f %%",  pct(I_pi0)),   "l");
    leg->AddEntry(h_PT3_eta,  Form("#eta Dalitz (17051):  %.2f %%",    pct(I_eta)),   "l");
    leg->AddEntry(h_PT3_rho,  Form("#rho^{0} (41) direct:  %.3f %%",   pct(I_rho)),   "l");
    leg->AddEntry(h_PT3_omgD, Form("#omega (52) direct:  %.3f %%",     pct(I_omgD)),  "l");
    leg->AddEntry(h_PT3_omgL, Form("#omega-Dalitz (52051):  %.3f %%",  pct(I_omgL)),  "l");
    leg->AddEntry(h_PT3_phi,  Form("#phi (55) direct:  %.4f %%",       pct(I_phi)),   "l");
    leg->Draw();

    // ── Col 2: ratios N_PT2/N_PT3, linear Y in [0, 3] ──
    c->cd(2);
    gPad->SetMargin(0.13, 0.04, 0.13, 0.10);

    r_tot ->GetYaxis()->SetRangeUser(0.0, 3.0);
    r_tot ->Draw("E1");
    r_pi0 ->Draw("E1 SAME");
    r_eta ->Draw("E1 SAME");
    r_rho ->Draw("E1 SAME");
    r_omgD->Draw("E1 SAME");
    r_omgL->Draw("E1 SAME");
    r_phi ->Draw("E1 SAME");

    // Reference y = 1 dashed grey line.
    TLine* lref = new TLine(kXmin, 1.0, kXmax, 1.0);
    lref->SetLineStyle(2);
    lref->SetLineColor(kGray + 2);
    lref->SetLineWidth(1);
    lref->Draw();

    // Multi-line TLatex listing each fit. χ²/ndf comparison IS the H1 test.
    TLatex tex;
    tex.SetNDC();
    tex.SetTextAlign(11);
    tex.SetTextSize(0.026);
    double yT = 0.88;
    auto annot = [&](Color_t col, const char* name, const FitResult& f) {
        tex.SetTextColor(col);
        tex.DrawLatex(0.16, yT,
            Form("%s: a = %.3f #pm %.3f   #chi^{2}/ndf = %.0f / %d",
                 name, f.a, f.e, f.chi2_ndf * f.ndf, f.ndf));
        yT -= 0.035;
    };
    annot(kColTotFit,      "total",     f_tot);
    annot(kColPi0Fit,      "#pi^{0} Dal", f_pi0);
    annot(kColEtaFit,      "#eta Dal",    f_eta);
    annot(kColRhoFit,      "#rho^{0} (41)", f_rho);
    annot(kColOmegaDirFit, "#omega (52)", f_omgD);
    annot(kColOmegaDalFit, "#omega-D (52051)", f_omgL);
    annot(kColPhiFit,      "#phi (55)", f_phi);

    // Legend (right side of pad 2).
    TLegend* leg2 = new TLegend(0.55, 0.50, 0.95, 0.88);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetTextSize(0.026);
    leg2->AddEntry(r_tot,  "total",                       "lpe");
    leg2->AddEntry(r_pi0,  "#pi^{0} Dalitz (7051)",       "lpe");
    leg2->AddEntry(r_eta,  "#eta Dalitz (17051)",         "lpe");
    leg2->AddEntry(r_rho,  "#rho^{0} direct (41)",        "lpe");
    leg2->AddEntry(r_omgD, "#omega direct (52)",          "lpe");
    leg2->AddEntry(r_omgL, "#omega Dalitz (52051)",       "lpe");
    leg2->AddEntry(r_phi,  "#phi direct (55)",            "lpe");
    leg2->Draw();

    // -------- Save --------
    gSystem->mkdir("plots/output", true);
    const std::string fr_suffix =
        custom_range ? Form("_fit%g-%g", fxlo, fxhi) : std::string{};
    const std::string base =
        "plots/output/mass_ee_trigger_ratio_per_channel_sim" + fr_suffix;
    c->SaveAs((base + ".pdf").c_str());
    c->SaveAs((base + ".png").c_str());
    std::cout << "\nSaved: " << base << ".{pdf,png}\n";

    fs->Close();
}
