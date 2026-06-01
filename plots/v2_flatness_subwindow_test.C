// v2_flatness_subwindow_test.C
// =========================================================================
// Three-subwindow flatness test for the trigger-correction constant V2.
//
// DESIGN PRINCIPLE
// ----------------
//   Fits are performed on the EXTENDED window [0.15, 0.80] using a doubled
//   resolution uniform binning (36 bins, ~18 MeV/bin). This gives more bins
//   and a longer lever arm for any slope.
//
//   The resulting CORRECTION is then applied to the PT3 signal histogram
//   on the SAME reference range [0.15, 0.70] that was used by the previous
//   V2 analysis, with the SAME 58-bin variable layout:
//       20 x 2.5 MeV bins on [0.00, 0.05]
//       38 x 35.5 MeV bins on [0.05, 1.40]
//
//   That way N_corr (this macro) is directly comparable to the previous
//   V2 result N_corr_V2_previous = 577,617.805 +/- 13,911.529.
//
// FITS  (all on uniform-36 layout, 63 * sig_PT2 / sig_PT3)
// --------
//   - C_global   from pol0 on [0.15, 0.80]
//   - C_A        from pol0 on [0.15, 0.40]
//   - C_B        from pol0 on [0.40, 0.60]
//   - C_C        from pol0 on [0.60, 0.80]
//   - pol1: a + b*x  on [0.15, 0.80]
//
// CONSISTENCY
// -----------
//   Pairwise n_sigma between (A,B), (B,C), (A,C); each vs. global;
//   pol1 slope significance |b|/sigma_b; chi^2/ndf for each.
//
// CORRECTED YIELD on reference range [0.15, 0.70]  (58-bin layout)
// ----------------------------------------------------------------
//   (1) N_corr_pol0_extended : N_raw * C_global
//   (2) N_corr_piecewise     : per-bin c(b) from subwindow C that
//                              contains the bin center
//                              (C_C is extrapolated from [0.60,0.80] to
//                              cover bins in [0.60,0.70]).
//   (3) N_corr_pol1          : per-bin c(b) = a + b * x_center
//                              with full (a,b) covariance.
//   (4) Reference V2 previous: 577,617.805 +/- 13,911.529.
//
// PLOT
// ----
//   Data markers, three blue subwindow lines+bands, red global pol0
//   line+band, dashed black pol1, vertical line at x=0.70 marking the
//   reference range upper edge, yields summarised in the legend.
//
// Usage:
//   root -l -b -q plots/v2_flatness_subwindow_test.C
//
// @author Witold Przygoda (witold.przygoda@uj.edu.pl)
// @date 2026

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TMatrixDSym.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TBox.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TMath.h>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

namespace {

// -------------- Inputs --------------
constexpr const char* kFileEm = "output_epem_exp.root";
constexpr const char* kFilePp = "output_epep_exp.root";
constexpr const char* kFileMm = "output_emem_exp.root";
constexpr const char* kTree   = "dilepton_nt";

// -------------- Trigger cuts --------------
constexpr const char* kCutPT3 = "trigbit==8192";
constexpr const char* kCutPT2 = "trigbit==4096";
constexpr double kTrigCorr = 63.0;

// -------------- Uniform binning for the FITS: 36 bins on [0.15, 0.80] -----
constexpr int    kNb   = 36;
constexpr double kXmin = 0.15;
constexpr double kXmax = 0.80;

// Sub-fit windows
constexpr double kAlo = 0.15, kAhi = 0.40;
constexpr double kBlo = 0.40, kBhi = 0.60;
constexpr double kClo = 0.60, kChi = 0.80;

// -------------- Variable 58-bin layout for the YIELD computation ---------
//   Matches mass_ee_pt_ratio_signal_varbins_exp.C and other V2 macros.
constexpr int    kNvar = 58;
constexpr double kRefLo = 0.15;   // reference range for N_corr
constexpr double kRefHi = 0.70;

// -------------- Previous V2 reference --------------
constexpr double kN_prev   = 577617.805;
constexpr double kSN_prev  = 13911.529;

// -------------- Build the 58-bin variable edges --------------
const double* getVarEdges() {
    static double edges[kNvar + 1];
    static bool built = false;
    if (!built) {
        int k = 0;
        for (int i = 0; i <= 20; ++i) edges[k++] = i * 0.0025;
        const double w = (1.4 - 0.05) / 38.0;
        for (int i = 1; i <= 38; ++i) edges[k++] = 0.05 + i * w;
        built = true;
    }
    return edges;
}

// -------------- Histogram helpers --------------
TH1D* drawWithCutUniform(TTree* t, const char* expr, const char* cut,
                         const std::string& name,
                         int nb, double xlo, double xhi) {
    TH1D* h = new TH1D(name.c_str(), "", nb, xlo, xhi);
    h->Sumw2();
    t->Draw((std::string(expr) + ">>" + name).c_str(), cut, "goff");
    h->SetDirectory(nullptr);
    return h;
}

TH1D* drawWithCutVar(TTree* t, const char* expr, const char* cut,
                     const std::string& name,
                     int nb, const double* edges) {
    TH1D* h = new TH1D(name.c_str(), "", nb, edges);
    h->Sumw2();
    t->Draw((std::string(expr) + ">>" + name).c_str(), cut, "goff");
    h->SetDirectory(nullptr);
    return h;
}

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
        const double err = val * 0.5 *
            std::sqrt(rel_pp * rel_pp + rel_mm * rel_mm);
        cb->SetBinContent(b, val);
        cb->SetBinError  (b, err);
    }
    return cb;
}

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

const char* verdict(double nsig) {
    if (nsig < 1.0) return "CONSISTENT";
    if (nsig < 2.0) return "MARGINAL";
    if (nsig < 3.0) return "TENSION";
    return "DISCREPANT";
}

}  // namespace

void v2_flatness_subwindow_test() {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    printf("=========================================================\n");
    printf("V2 flatness sub-window test  (fits on [%.2f, %.2f] GeV/c^2)\n",
           kXmin, kXmax);
    printf("Fit binning: %d uniform bins, width = %.6f GeV (%.2f MeV)\n",
           kNb, (kXmax - kXmin) / kNb, 1000.0 * (kXmax - kXmin) / kNb);
    printf("Yield binning: 58 variable bins, reference [%.2f, %.2f]\n",
           kRefLo, kRefHi);
    printf("=========================================================\n\n");

    // ----- Open EXP inputs -----
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
        std::cerr << "Tree '" << kTree << "' missing in EXP files.\n";
        return;
    }

    // ===================================================================
    // Build the 36-bin uniform PT2/PT3 ratio for the FITS.
    // ===================================================================
    TH1D* h_em_PT3 = drawWithCutUniform(t_em, "m_ee", kCutPT3,
                                        "h_em_PT3_36", kNb, kXmin, kXmax);
    TH1D* h_em_PT2 = drawWithCutUniform(t_em, "m_ee", kCutPT2,
                                        "h_em_PT2_36", kNb, kXmin, kXmax);
    TH1D* h_pp_PT3 = drawWithCutUniform(t_pp, "m_ee", kCutPT3,
                                        "h_pp_PT3_36", kNb, kXmin, kXmax);
    TH1D* h_pp_PT2 = drawWithCutUniform(t_pp, "m_ee", kCutPT2,
                                        "h_pp_PT2_36", kNb, kXmin, kXmax);
    TH1D* h_mm_PT3 = drawWithCutUniform(t_mm, "m_ee", kCutPT3,
                                        "h_mm_PT3_36", kNb, kXmin, kXmax);
    TH1D* h_mm_PT2 = drawWithCutUniform(t_mm, "m_ee", kCutPT2,
                                        "h_mm_PT2_36", kNb, kXmin, kXmax);

    TH1D* h_cb_PT3 = makeCB(h_pp_PT3, h_mm_PT3, "h_cb_PT3_36");
    TH1D* h_cb_PT2 = makeCB(h_pp_PT2, h_mm_PT2, "h_cb_PT2_36");
    TH1D* h_pt3_signal_36 =
        makeSignal(h_em_PT3, h_cb_PT3, "h_pt3_signal_36");
    TH1D* h_pt2_signal_36 =
        makeSignal(h_em_PT2, h_cb_PT2, "h_pt2_signal_36");

    TH1D* h_exp_ratio_36 = makeRatio(h_pt2_signal_36, h_pt3_signal_36,
                                     kTrigCorr, "h_exp_ratio_36");

    // ===================================================================
    // Build the 58-bin variable PT3 signal for the YIELD on [0.15, 0.70].
    // ===================================================================
    const double* edges = getVarEdges();
    TH1D* h_em_PT3_58 = drawWithCutVar(t_em, "m_ee", kCutPT3,
                                       "h_em_PT3_58", kNvar, edges);
    TH1D* h_pp_PT3_58 = drawWithCutVar(t_pp, "m_ee", kCutPT3,
                                       "h_pp_PT3_58", kNvar, edges);
    TH1D* h_mm_PT3_58 = drawWithCutVar(t_mm, "m_ee", kCutPT3,
                                       "h_mm_PT3_58", kNvar, edges);
    TH1D* h_cb_PT3_58 = makeCB(h_pp_PT3_58, h_mm_PT3_58, "h_cb_PT3_58");
    TH1D* h_pt3_signal_58 =
        makeSignal(h_em_PT3_58, h_cb_PT3_58, "h_pt3_signal_58");

    // Identify the bins in the reference range and compute N_raw.
    std::vector<int> ref_bins;
    ref_bins.reserve(kNvar);
    double N_raw = 0.0, var_raw = 0.0;
    for (int b = 1; b <= h_pt3_signal_58->GetNbinsX(); ++b) {
        const double xc = h_pt3_signal_58->GetXaxis()->GetBinCenter(b);
        if (xc < kRefLo - 1e-9 || xc > kRefHi + 1e-9) continue;
        ref_bins.push_back(b);
        const double n = h_pt3_signal_58->GetBinContent(b);
        const double e = h_pt3_signal_58->GetBinError(b);
        N_raw   += n;
        var_raw += e * e;
    }
    const double sN_raw = std::sqrt(var_raw);
    printf("Reference-range yield (raw PT3 signal, 58-bin layout):\n");
    printf("  number of bins in [%.2f, %.2f] = %zu\n",
           kRefLo, kRefHi, ref_bins.size());
    if (!ref_bins.empty()) {
        printf("  first bin center = %.5f, last bin center = %.5f\n",
               h_pt3_signal_58->GetXaxis()->GetBinCenter(ref_bins.front()),
               h_pt3_signal_58->GetXaxis()->GetBinCenter(ref_bins.back()));
    }
    printf("  N_raw [0.15, 0.70] = %.3f  +/-  %.3f\n\n",
           N_raw, sN_raw);

    // ========================================================================
    // ANALYSIS A - Global V2 over [0.15, 0.80]
    // ========================================================================
    printf("=========================================================\n");
    printf("ANALYSIS A - Global V2 over [%.2f, %.2f]\n", kXmin, kXmax);
    printf("=========================================================\n");

    TF1* f_glob = new TF1("f_glob", "[0]", kXmin, kXmax);
    f_glob->SetParameter(0, 1.5);
    h_exp_ratio_36->Fit(f_glob, "RQ0", "", kXmin, kXmax);
    const double C_glob    = f_glob->GetParameter(0);
    const double sC_glob   = f_glob->GetParError(0);
    const int    ndf_glob  = f_glob->GetNDF();
    const double chi2_glob = f_glob->GetChisquare();
    const double cn_glob   = (ndf_glob > 0) ? chi2_glob / ndf_glob : 0.0;

    printf("  C_global       = %.6f\n", C_glob);
    printf("  sigma_C_global = %.6f\n", sC_glob);
    printf("  chi2 / ndf     = %.4f / %d = %.4f\n",
           chi2_glob, ndf_glob, cn_glob);

    // ========================================================================
    // ANALYSIS B - Three subwindow pol0 fits
    // ========================================================================
    printf("\n=========================================================\n");
    printf("ANALYSIS B - Three-subwindow pol0 fits\n");
    printf("=========================================================\n");

    auto fitSub = [&](const char* name, double xlo, double xhi,
                      double& C, double& sC, double& chi2, int& ndf,
                      int& nbins_used) {
        TF1* f = new TF1(name, "[0]", xlo, xhi);
        f->SetParameter(0, 1.5);
        h_exp_ratio_36->Fit(f, "RQ0", "", xlo, xhi);
        C    = f->GetParameter(0);
        sC   = f->GetParError(0);
        ndf  = f->GetNDF();
        chi2 = f->GetChisquare();
        nbins_used = 0;
        for (int b = 1; b <= h_exp_ratio_36->GetNbinsX(); ++b) {
            const double xc = h_exp_ratio_36->GetXaxis()->GetBinCenter(b);
            if (xc >= xlo - 1e-9 && xc <= xhi + 1e-9) ++nbins_used;
        }
    };

    double C_A, sC_A, chi2_A; int ndf_A, n_A;
    double C_B, sC_B, chi2_B; int ndf_B, n_B;
    double C_C, sC_C, chi2_C; int ndf_C, n_C;
    fitSub("f_A", kAlo, kAhi, C_A, sC_A, chi2_A, ndf_A, n_A);
    fitSub("f_B", kBlo, kBhi, C_B, sC_B, chi2_B, ndf_B, n_B);
    fitSub("f_C", kClo, kChi, C_C, sC_C, chi2_C, ndf_C, n_C);

    auto report = [](const char* lab, double xlo, double xhi,
                     double C, double sC, double chi2, int ndf, int nb) {
        const double cn = (ndf > 0) ? chi2 / ndf : 0.0;
        printf("  Subwindow %s [%.2f, %.2f]: bins=%d\n", lab, xlo, xhi, nb);
        printf("    C   = %.6f +/- %.6f\n", C, sC);
        printf("    chi2/ndf = %.4f / %d = %.4f\n", chi2, ndf, cn);
    };
    report("A", kAlo, kAhi, C_A, sC_A, chi2_A, ndf_A, n_A);
    report("B", kBlo, kBhi, C_B, sC_B, chi2_B, ndf_B, n_B);
    report("C", kClo, kChi, C_C, sC_C, chi2_C, ndf_C, n_C);

    // ========================================================================
    // ANALYSIS C - Pairwise n_sigma consistency tests
    // ========================================================================
    printf("\n=========================================================\n");
    printf("ANALYSIS C - Pairwise consistency tests\n");
    printf("=========================================================\n");

    auto nsig = [](double C1, double s1, double C2, double s2) {
        const double d   = std::fabs(C1 - C2);
        const double s   = std::sqrt(s1 * s1 + s2 * s2);
        return std::make_pair(d, (s > 0.0) ? d / s : 0.0);
    };

    auto compare = [&](const char* tag, double C1, double s1,
                       double C2, double s2) {
        double d, n;
        std::tie(d, n) = nsig(C1, s1, C2, s2);
        const double sc = std::sqrt(s1 * s1 + s2 * s2);
        printf("  %s: |%.4f - %.4f| = %.4f, sigma_comb = %.4f, n_sigma = %.3f [%s]\n",
               tag, C1, C2, d, sc, n, verdict(n));
        return n;
    };

    printf("  Subwindow vs. subwindow:\n");
    const double n_AB = compare("A vs B", C_A, sC_A, C_B, sC_B);
    const double n_BC = compare("B vs C", C_B, sC_B, C_C, sC_C);
    const double n_AC = compare("A vs C", C_A, sC_A, C_C, sC_C);

    printf("  Subwindow vs. global:\n");
    const double n_Ag = compare("A vs Global", C_A, sC_A, C_glob, sC_glob);
    const double n_Bg = compare("B vs Global", C_B, sC_B, C_glob, sC_glob);
    const double n_Cg = compare("C vs Global", C_C, sC_C, C_glob, sC_glob);

    // ========================================================================
    // ANALYSIS D - pol1 fit as a shape test
    // ========================================================================
    printf("\n=========================================================\n");
    printf("ANALYSIS D - pol1 fit (slope significance)\n");
    printf("=========================================================\n");

    TF1* f_lin = new TF1("f_lin", "[0]+[1]*x", kXmin, kXmax);
    f_lin->SetParameter(0, 1.5);
    f_lin->SetParameter(1, 0.0);
    // Use "S" to retrieve the full TFitResult (covariance matrix).
    TFitResultPtr fitres_lin =
        h_exp_ratio_36->Fit(f_lin, "RQ0S", "", kXmin, kXmax);
    const double a_lin    = f_lin->GetParameter(0);
    const double sa_lin   = f_lin->GetParError(0);
    const double b_lin    = f_lin->GetParameter(1);
    const double sb_lin   = f_lin->GetParError(1);
    const int    ndf_lin  = f_lin->GetNDF();
    const double chi2_lin = f_lin->GetChisquare();
    const double cn_lin   = (ndf_lin > 0) ? chi2_lin / ndf_lin : 0.0;
    const double slope_sig = (sb_lin > 0.0) ? std::fabs(b_lin) / sb_lin : 0.0;

    double cov_ab = 0.0;
    if (fitres_lin.Get()) {
        TMatrixDSym cov = fitres_lin->GetCovarianceMatrix();
        if (cov.GetNrows() >= 2 && cov.GetNcols() >= 2) {
            cov_ab = cov(0, 1);
        }
    }

    printf("  a (intercept) = %.6f +/- %.6f\n", a_lin, sa_lin);
    printf("  b (slope)     = %.6f +/- %.6f\n", b_lin, sb_lin);
    printf("  cov(a,b)      = %.6e\n", cov_ab);
    printf("  chi2/ndf      = %.4f / %d = %.4f\n",
           chi2_lin, ndf_lin, cn_lin);
    printf("  |b|/sigma_b   = %.3f sigma\n", slope_sig);
    if (slope_sig > 2.0) {
        printf("  -> slope > 2 sigma  : shape is mass-dependent.\n");
    } else if (slope_sig < 1.0) {
        printf("  -> slope < 1 sigma  : shape is consistent with flat.\n");
    } else {
        printf("  -> slope in [1, 2] sigma: marginal hint of structure.\n");
    }

    // ========================================================================
    // YIELDS over reference [0.15, 0.70] on the 58-bin layout.
    // ========================================================================
    printf("\n=========================================================\n");
    printf("YIELDS over reference [%.2f, %.2f]  (58-bin layout)\n",
           kRefLo, kRefHi);
    printf("=========================================================\n");

    // ---- (1) N_corr_pol0_extended = N_raw * C_global -----------------------
    const double N_corr_glob = N_raw * C_glob;
    // Conservative: ignore N_raw error => sigma = N_raw * sigma_C_global.
    const double sN_corr_glob = N_raw * sC_glob;

    // ---- (2) N_corr_piecewise: per-bin subwindow C -------------------------
    //   C_A on [0.00, 0.40), C_B on [0.40, 0.60), C_C on [0.60, ...]
    //   (C_C is fit on [0.60, 0.80] but used here only on [0.60, 0.70]).
    //   Error propagation: subwindow C errors are uncorrelated with each
    //   other AND with h_pt3_signal_58 (the ratio is built from PT2 + the
    //   same PT3, but the subwindow-C errors are taken as exogenous here,
    //   consistent with how V2 was used previously).
    double N_corr_pw = 0.0;
    double var_pw    = 0.0;
    int    n_pw_A = 0, n_pw_B = 0, n_pw_C = 0;
    for (int b : ref_bins) {
        const double xc = h_pt3_signal_58->GetXaxis()->GetBinCenter(b);
        const double sig = h_pt3_signal_58->GetBinContent(b);
        double cval = 0.0, csig = 0.0;
        if (xc < kAhi) {
            cval = C_A; csig = sC_A; ++n_pw_A;
        } else if (xc < kBhi) {
            cval = C_B; csig = sC_B; ++n_pw_B;
        } else {
            cval = C_C; csig = sC_C; ++n_pw_C;
        }
        N_corr_pw += sig * cval;
        // Per-bin variance contribution: ( sig * sigma_C )^2
        // (treating sig as exact, like for the global case).
        var_pw += (sig * csig) * (sig * csig);
    }
    const double sN_corr_pw = std::sqrt(var_pw);

    // ---- (3) N_corr_pol1: per-bin c(b) = a + b*x_center --------------------
    //   For full error: per-bin c-variance is
    //       sigma_c^2(x) = sa^2 + 2 x cov(a,b) + x^2 sb^2.
    //   For the SUM = Sum_b sig_b * (a + b * x_b) we have
    //       Sum = a * S0 + b * S1 ,  where S0 = sum sig_b, S1 = sum sig_b * x_b.
    //   Then sigma_Sum^2 = S0^2 sa^2 + 2 S0 S1 cov(a,b) + S1^2 sb^2.
    double S0 = 0.0, S1 = 0.0;
    for (int b : ref_bins) {
        const double xc = h_pt3_signal_58->GetXaxis()->GetBinCenter(b);
        const double sig = h_pt3_signal_58->GetBinContent(b);
        S0 += sig;
        S1 += sig * xc;
    }
    const double N_corr_lin = a_lin * S0 + b_lin * S1;
    const double var_lin =
        S0 * S0 * sa_lin * sa_lin
        + 2.0 * S0 * S1 * cov_ab
        + S1 * S1 * sb_lin * sb_lin;
    const double sN_corr_lin = (var_lin > 0.0) ? std::sqrt(var_lin) : 0.0;

    printf("  (1) N_corr_pol0_extended  (C_global on [0.15,0.80]):\n");
    printf("        N = %.3f  +/-  %.3f       sigma/N = %.3f %%\n",
           N_corr_glob, sN_corr_glob,
           100.0 * sN_corr_glob / std::fabs(N_corr_glob));

    printf("  (2) N_corr_piecewise  (bins: A=%d, B=%d, C=%d):\n",
           n_pw_A, n_pw_B, n_pw_C);
    printf("        N = %.3f  +/-  %.3f       sigma/N = %.3f %%\n",
           N_corr_pw, sN_corr_pw,
           100.0 * sN_corr_pw / std::fabs(N_corr_pw));

    printf("  (3) N_corr_pol1  (S0=sum sig=%.3f, S1=sum sig*x=%.3f):\n",
           S0, S1);
    printf("        N = %.3f  +/-  %.3f       sigma/N = %.3f %%\n",
           N_corr_lin, sN_corr_lin,
           100.0 * sN_corr_lin / std::fabs(N_corr_lin));

    printf("  (4) Reference (previous V2 [0.15, 0.70]):\n");
    printf("        N = %.3f  +/-  %.3f       sigma/N = %.3f %%\n",
           kN_prev, kSN_prev, 100.0 * kSN_prev / kN_prev);

    // ----- Deltas vs the previous V2 -----------------------------------------
    auto delta_report = [](const char* tag, double Nnew, double Snew) {
        const double d = Nnew - kN_prev;
        const double sc = std::sqrt(Snew * Snew + kSN_prev * kSN_prev);
        const double nsig = (sc > 0.0) ? std::fabs(d) / sc : 0.0;
        printf("  Delta %s vs V2_prev = %+.3f   sigma_comb = %.3f"
               "   |Delta|/sigma = %.3f  [%s]\n",
               tag, d, sc, nsig, verdict(nsig));
        return std::make_pair(d, nsig);
    };

    printf("\n  Yield deltas vs previous V2:\n");
    double dGlob, nGlob;
    std::tie(dGlob, nGlob) = delta_report("(pol0_ext)", N_corr_glob, sN_corr_glob);
    double dPw,   nPw;
    std::tie(dPw,   nPw)   = delta_report("(piecewise)", N_corr_pw,  sN_corr_pw);
    double dLin,  nLin;
    std::tie(dLin,  nLin)  = delta_report("(pol1     )", N_corr_lin, sN_corr_lin);

    // ========================================================================
    // PLOT
    // ========================================================================
    gSystem->mkdir("plots/output", true);

    TCanvas* c1 = new TCanvas("c1", "v2_flatness_subwindow_test",
                              1100, 750);
    c1->SetLeftMargin(0.12);
    c1->SetBottomMargin(0.12);
    c1->SetRightMargin(0.05);
    c1->SetTopMargin(0.08);
    c1->SetGridx();
    c1->SetGridy();

    double ymin = 1e30, ymax = -1e30;
    for (int b = 1; b <= h_exp_ratio_36->GetNbinsX(); ++b) {
        const double y = h_exp_ratio_36->GetBinContent(b);
        const double e = h_exp_ratio_36->GetBinError(b);
        if (h_exp_ratio_36->GetBinContent(b) == 0.0 && e == 0.0) continue;
        if (y - e < ymin) ymin = y - e;
        if (y + e > ymax) ymax = y + e;
    }
    if (!(ymin < ymax)) { ymin = 0.0; ymax = 3.0; }
    const double yspan = ymax - ymin;
    ymin -= 0.15 * yspan;
    ymax += 0.45 * yspan;   // leave more headroom for legend + text box

    h_exp_ratio_36->SetTitle(
        "V2 flatness test, m_{ee} #in [0.15, 0.80] GeV/c^{2}");
    h_exp_ratio_36->GetXaxis()->SetTitle("m_{ee} [GeV/c^{2}]");
    h_exp_ratio_36->GetYaxis()->SetTitle("63 #times N_{PT2} / N_{PT3}");
    h_exp_ratio_36->GetXaxis()->SetRangeUser(kXmin, kXmax);
    h_exp_ratio_36->GetYaxis()->SetRangeUser(ymin, ymax);
    h_exp_ratio_36->SetMarkerStyle(20);
    h_exp_ratio_36->SetMarkerSize(0.9);
    h_exp_ratio_36->SetMarkerColor(kBlack);
    h_exp_ratio_36->SetLineColor(kBlack);
    h_exp_ratio_36->Draw("E1");

    // Global pol0: red error band + solid red line over full range.
    TBox* box_glob = new TBox(kXmin, C_glob - sC_glob,
                              kXmax, C_glob + sC_glob);
    box_glob->SetFillColorAlpha(kRed, 0.18);
    box_glob->SetLineColor(0);
    box_glob->Draw("same");

    TLine* line_glob = new TLine(kXmin, C_glob, kXmax, C_glob);
    line_glob->SetLineColor(kRed);
    line_glob->SetLineWidth(3);
    line_glob->Draw("same");

    // Subwindow pol0s: blue lines + error bands over their fit ranges.
    auto drawSub = [&](double xlo, double xhi, double C, double sC) {
        TBox* bx = new TBox(xlo, C - sC, xhi, C + sC);
        bx->SetFillColorAlpha(kBlue, 0.18);
        bx->SetLineColor(0);
        bx->Draw("same");
        TLine* ln = new TLine(xlo, C, xhi, C);
        ln->SetLineColor(kBlue + 1);
        ln->SetLineWidth(3);
        ln->Draw("same");
    };
    drawSub(kAlo, kAhi, C_A, sC_A);
    drawSub(kBlo, kBhi, C_B, sC_B);
    drawSub(kClo, kChi, C_C, sC_C);

    // pol1: dashed black line over [0.15, 0.80].
    TF1* f_lin_draw = new TF1("f_lin_draw", "[0]+[1]*x", kXmin, kXmax);
    f_lin_draw->SetParameter(0, a_lin);
    f_lin_draw->SetParameter(1, b_lin);
    f_lin_draw->SetLineColor(kBlack);
    f_lin_draw->SetLineStyle(2);
    f_lin_draw->SetLineWidth(2);
    f_lin_draw->Draw("same");

    // Vertical line at x=0.70 marking the reference range upper edge.
    TLine* line_ref = new TLine(kRefHi, ymin, kRefHi, ymax);
    line_ref->SetLineColor(kGreen + 2);
    line_ref->SetLineStyle(7);
    line_ref->SetLineWidth(2);
    line_ref->Draw("same");

    // Re-draw data on top so markers are visible.
    h_exp_ratio_36->Draw("E1 same");

    // Legend (top right) with fit constants + slope.
    TLegend* leg = new TLegend(0.55, 0.62, 0.945, 0.91);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.028);
    leg->AddEntry(h_exp_ratio_36, "63 #times sig_{PT2} / sig_{PT3} (data)",
                  "lep");
    leg->AddEntry(line_glob,
                  Form("C_{global}[0.15,0.80] = %.3f #pm %.3f (#chi^{2}/ndf=%.2f)",
                       C_glob, sC_glob, cn_glob), "l");
    {
        TLine* lA = new TLine(); lA->SetLineColor(kBlue + 1);
        lA->SetLineWidth(3);
        leg->AddEntry(lA,
            Form("C_{A}[0.15,0.40] = %.3f #pm %.3f", C_A, sC_A), "l");
        TLine* lB = new TLine(); lB->SetLineColor(kBlue + 1);
        lB->SetLineWidth(3);
        leg->AddEntry(lB,
            Form("C_{B}[0.40,0.60] = %.3f #pm %.3f", C_B, sC_B), "l");
        TLine* lC = new TLine(); lC->SetLineColor(kBlue + 1);
        lC->SetLineWidth(3);
        leg->AddEntry(lC,
            Form("C_{C}[0.60,0.80] = %.3f #pm %.3f", C_C, sC_C), "l");
    }
    leg->AddEntry(f_lin_draw,
                  Form("pol1: a=%.3f, b=%.3f #pm %.3f (%.2f#sigma)",
                       a_lin, b_lin, sb_lin, slope_sig), "l");
    leg->AddEntry(line_ref,
                  "reference upper edge x=0.70", "l");
    leg->Draw();

    // Yields text box (top left).
    TLatex tx;
    tx.SetNDC(true);
    tx.SetTextFont(42);
    tx.SetTextSize(0.024);
    double y_tx = 0.89;
    const double dy = 0.032;
    tx.DrawLatex(0.14, y_tx,
                 Form("N_{raw}[0.15,0.70] = %.0f", N_raw));               y_tx -= dy;
    tx.DrawLatex(0.14, y_tx,
                 Form("N_{V2 prev} = %.0f #pm %.0f",
                      kN_prev, kSN_prev));                                y_tx -= dy;
    tx.DrawLatex(0.14, y_tx,
                 Form("N_{pol0 ext} = %.0f #pm %.0f  (#Delta=%+.0f, %.2f#sigma)",
                      N_corr_glob, sN_corr_glob, dGlob, nGlob));          y_tx -= dy;
    tx.DrawLatex(0.14, y_tx,
                 Form("N_{piecewise} = %.0f #pm %.0f  (#Delta=%+.0f, %.2f#sigma)",
                      N_corr_pw, sN_corr_pw, dPw, nPw));                  y_tx -= dy;
    tx.DrawLatex(0.14, y_tx,
                 Form("N_{pol1} = %.0f #pm %.0f  (#Delta=%+.0f, %.2f#sigma)",
                      N_corr_lin, sN_corr_lin, dLin, nLin));

    const std::string out_pdf = "plots/output/v2_flatness_subwindow_test.pdf";
    const std::string out_png = "plots/output/v2_flatness_subwindow_test.png";
    c1->SaveAs(out_pdf.c_str());
    c1->SaveAs(out_png.c_str());

    // ========================================================================
    // Summary block
    // ========================================================================
    printf("\n=========================================================\n");
    printf("SUMMARY - Fits (extended range, 36-bin uniform)\n");
    printf("=========================================================\n");
    printf("  C_global [0.15, 0.80] = %.4f +/- %.4f  (chi2/ndf = %.3f)\n",
           C_glob, sC_glob, cn_glob);
    printf("  C_A      [0.15, 0.40] = %.4f +/- %.4f  (chi2/ndf = %.3f)\n",
           C_A, sC_A, (ndf_A > 0) ? chi2_A / ndf_A : 0.0);
    printf("  C_B      [0.40, 0.60] = %.4f +/- %.4f  (chi2/ndf = %.3f)\n",
           C_B, sC_B, (ndf_B > 0) ? chi2_B / ndf_B : 0.0);
    printf("  C_C      [0.60, 0.80] = %.4f +/- %.4f  (chi2/ndf = %.3f)\n",
           C_C, sC_C, (ndf_C > 0) ? chi2_C / ndf_C : 0.0);
    printf("  pol1 a = %.4f +/- %.4f   b = %.4f +/- %.4f   |b|/sigma = %.2f\n",
           a_lin, sa_lin, b_lin, sb_lin, slope_sig);
    printf("  Pairwise n_sigma: A-B=%.2f, B-C=%.2f, A-C=%.2f\n",
           n_AB, n_BC, n_AC);
    printf("  Vs global       : A-G=%.2f, B-G=%.2f, C-G=%.2f\n",
           n_Ag, n_Bg, n_Cg);

    printf("\n=========================================================\n");
    printf("FINAL TABLE - N_corr over reference [0.15, 0.70]\n");
    printf("=========================================================\n");
    auto pct = [](double s, double n) {
        return (n != 0.0) ? 100.0 * s / std::fabs(n) : 0.0;
    };
    printf("  %-22s | %-22s | %12s | %12s | %7s | %12s\n",
           "Method", "C used", "N_corr", "sigma", "sigma/N%", "Delta vs prev");
    printf("  %-22s | %-22s | %12.0f | %12.0f | %7.3f | %12s\n",
           "V2 previous", "C_prev~1.858", kN_prev, kSN_prev, pct(kSN_prev, kN_prev), "0");
    printf("  %-22s | %-22s | %12.0f | %12.0f | %7.3f | %+12.0f\n",
           "V2 extended (pol0)",
           Form("C_global=%.3f", C_glob),
           N_corr_glob, sN_corr_glob, pct(sN_corr_glob, N_corr_glob), dGlob);
    printf("  %-22s | %-22s | %12.0f | %12.0f | %7.3f | %+12.0f\n",
           "Piecewise (3 sub)",
           Form("C_A=%.3f,B=%.3f,C=%.3f", C_A, C_B, C_C),
           N_corr_pw, sN_corr_pw, pct(sN_corr_pw, N_corr_pw), dPw);
    printf("  %-22s | %-22s | %12.0f | %12.0f | %7.3f | %+12.0f\n",
           "pol1 (linear)",
           Form("a=%.3f,b=%.3f", a_lin, b_lin),
           N_corr_lin, sN_corr_lin, pct(sN_corr_lin, N_corr_lin), dLin);

    printf("\n  n_sigma shifts vs previous V2:\n");
    printf("    pol0_ext  : %.3f sigma\n", nGlob);
    printf("    piecewise : %.3f sigma\n", nPw);
    printf("    pol1      : %.3f sigma\n", nLin);

    printf("\n  Plots: %s , %s\n", out_pdf.c_str(), out_png.c_str());
    printf("=========================================================\n");
    printf("DONE.\n");
    printf("=========================================================\n");

    f_em->Close();
    f_pp->Close();
    f_mm->Close();
}
