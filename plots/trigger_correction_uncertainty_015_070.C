// trigger_correction_uncertainty_015_070.C
// =========================================================================
// Clean three-variant uncertainty analysis for the trigger correction in
// the [0.15, 0.70] m_ee window using the 58-bin variable layout.
//
//   VARIANT 1 — Bin-by-bin EXP correction
//       corrected(b) = h_pt3_signal(b) * exp_ratio(b)
//       per-bin error from full propagation, summed in quadrature.
//
//   VARIANT 2 — Constant C from pol0 fit to EXP over [0.15, 0.70]
//       N_corr = N_raw_range * C
//
//   VARIANT 3 — Three sim shape variants (integral-matched to EXP)
//         3a) step4 FIT   : k_step4_fit  * simShape_step4(x)
//         3b) step5 FIT   : k_step5_fit  * simShape_step5(x)
//         3c) step5 PTS   : k_step5_pts  * sim_ratio_step5_at(x)
//       Spread of the three -> systematic from sim shape choice.
//
// Usage:
//   root -l -b -q plots/trigger_correction_uncertainty_015_070.C
//
// @author Witold Przygoda (witold.przygoda@uj.edu.pl)
// @date 2026

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TMath.h>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

namespace {

// -------------- Inputs --------------
constexpr const char* kFileEm = "output_epem_exp.root";
constexpr const char* kFilePp = "output_epep_exp.root";
constexpr const char* kFileMm = "output_emem_exp.root";
constexpr const char* kFileSim = "output_epem_sim.root";
constexpr const char* kTree   = "dilepton_nt";

// -------------- Trigger cuts --------------
constexpr const char* kCutPT3 = "trigbit==8192";
constexpr const char* kCutPT2 = "trigbit==4096";
constexpr double kTrigCorr = 63.0;

// Sim cuts.
constexpr const char* kSimStep5Cut =
    "(eVertReco_z>-500)*(isBest==1)*(ep_sim_id==2)*(em_sim_id==3)*"
    "(epem_same_vertex==1)";

// -------------- EXP binning (58 var bins, same as ratio macro) --------------
constexpr int    kNb   = 58;
constexpr double kXmin = 0.0;
constexpr double kXmax = 1.4;

// -------------- Sim cascade-v5 220-bin layout --------------
constexpr int kSimNb = 220;
constexpr double kSimMinDen = 5.0;

// -------------- Integration window --------------
constexpr double kIntLo = 0.15;
constexpr double kIntHi = 0.70;

// -------------- SIM FIT shapes (step4 and step5) --------------
constexpr double SIM_JOIN_X = 0.70;
// step4 REC piecewise:
constexpr double S4_HI_A  = 1.080;
constexpr double S4_HI_B  = 0.039;
constexpr double S4_HI_C  = 0.117;
constexpr double S4_MID_B = 0.264;
constexpr double S4_MID_C = 0.250;
const double S4_JOIN_Y =
    S4_HI_A + S4_HI_B * SIM_JOIN_X + S4_HI_C * SIM_JOIN_X * SIM_JOIN_X;

// step5 REC piecewise:
constexpr double S5_HI_A  =  1.467;
constexpr double S5_HI_B  = -0.803;
constexpr double S5_HI_C  =  0.738;
constexpr double S5_MID_B =  0.806;
constexpr double S5_MID_C =  1.444;
const double S5_JOIN_Y =
    S5_HI_A + S5_HI_B * SIM_JOIN_X + S5_HI_C * SIM_JOIN_X * SIM_JOIN_X;

double simShapeStep4(double x) {
    if (x < SIM_JOIN_X) {
        const double d = x - SIM_JOIN_X;
        return S4_JOIN_Y + S4_MID_B * d + S4_MID_C * d * d;
    }
    return S4_HI_A + S4_HI_B * x + S4_HI_C * x * x;
}

double simShapeStep5(double x) {
    if (x < SIM_JOIN_X) {
        const double d = x - SIM_JOIN_X;
        return S5_JOIN_Y + S5_MID_B * d + S5_MID_C * d * d;
    }
    return S5_HI_A + S5_HI_B * x + S5_HI_C * x * x;
}

// -------------- Build the 58-bin variable edges --------------
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

// -------------- Build cascade-v5 220-bin sim edges --------------
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

// -------------- Histogram helpers --------------
TH1D* drawWithCut(TTree* t, const char* expr, const char* cut,
                  const std::string& name, int nb, const double* edges) {
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

// Closed-form antiderivative for the piecewise sim shape.
double antiMid(double x, double joinY, double midB, double midC) {
    const double d = x - SIM_JOIN_X;
    return joinY * x + 0.5 * midB * d * d + (1.0 / 3.0) * midC * d * d * d;
}
double antiHigh(double x, double hiA, double hiB, double hiC) {
    return hiA * x + 0.5 * hiB * x * x + (1.0 / 3.0) * hiC * x * x * x;
}
double simIntegral(double a, double b,
                   double joinY,
                   double midB, double midC,
                   double hiA,  double hiB, double hiC) {
    if (b <= a) return 0.0;
    const double J = SIM_JOIN_X;
    if (b <= J) {
        return antiMid(b, joinY, midB, midC) - antiMid(a, joinY, midB, midC);
    } else if (a >= J) {
        return antiHigh(b, hiA, hiB, hiC) - antiHigh(a, hiA, hiB, hiC);
    } else {
        return (antiMid(J, joinY, midB, midC) - antiMid(a, joinY, midB, midC))
             + (antiHigh(b, hiA, hiB, hiC)    - antiHigh(J, hiA, hiB, hiC));
    }
}

double simIntegralStep4(double a, double b) {
    return simIntegral(a, b, S4_JOIN_Y, S4_MID_B, S4_MID_C,
                       S4_HI_A, S4_HI_B, S4_HI_C);
}
double simIntegralStep5(double a, double b) {
    return simIntegral(a, b, S5_JOIN_Y, S5_MID_B, S5_MID_C,
                       S5_HI_A, S5_HI_B, S5_HI_C);
}

// Build sim ratio histogram (cascade-v5 220-bin layout).
TH1D* buildSimFullRatio(TTree* t_sim, const char* sim_cut,
                        const std::string& tag) {
    const std::string nm3 = "h_sim_full_PT3_" + tag;
    const std::string nm2 = "h_sim_full_PT2_" + tag;
    TH1D* h_sim_PT3 = new TH1D(nm3.c_str(), "", kSimNb, getSimEdges());
    TH1D* h_sim_PT2 = new TH1D(nm2.c_str(), "", kSimNb, getSimEdges());
    h_sim_PT3->Sumw2();
    h_sim_PT2->Sumw2();
    const std::string w_PT3 =
        std::string("(pt3==1)*") + sim_cut + "*sim_genweight";
    const std::string w_PT2 =
        std::string("(pt2==1)*") + sim_cut + "*sim_genweight";
    t_sim->Draw(("m_ee>>" + nm3).c_str(), w_PT3.c_str(), "goff");
    t_sim->Draw(("m_ee>>" + nm2).c_str(), w_PT2.c_str(), "goff");
    h_sim_PT3->SetDirectory(nullptr);
    h_sim_PT2->SetDirectory(nullptr);

    TH1D* r = static_cast<TH1D*>(
        h_sim_PT3->Clone(("h_sim_full_ratio_" + tag).c_str()));
    r->SetDirectory(nullptr);
    r->Reset();
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
    }
    delete h_sim_PT3;
    delete h_sim_PT2;
    return r;
}

double simRatioAt(TH1D* h_sim, double x) {
    if (!h_sim) return 0.0;
    const int b = h_sim->FindBin(x);
    if (b < 1 || b > h_sim->GetNbinsX()) return 0.0;
    return h_sim->GetBinContent(b);
}

}  // anonymous namespace

void trigger_correction_uncertainty_015_070() {
    printf("=========================================================\n");
    printf("Three-variant trigger correction uncertainty in [%.2f, %.2f]\n",
           kIntLo, kIntHi);
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

    // ----- Open SIM input -----
    TFile* f_sim = TFile::Open(kFileSim, "READ");
    if (!f_sim || f_sim->IsZombie()) {
        std::cerr << "Cannot open SIM input " << kFileSim << "\n";
        return;
    }
    TTree* t_sim = dynamic_cast<TTree*>(f_sim->Get(kTree));
    if (!t_sim) {
        std::cerr << "Tree '" << kTree << "' missing in SIM file.\n";
        return;
    }

    const double* edges = getEdges();

    // ----- Build PT3/PT2 epem, ++ ,  -- on 58-bin layout -----
    TH1D* h_em_PT3 = drawWithCut(t_em, "m_ee", kCutPT3, "h_em_PT3", kNb, edges);
    TH1D* h_em_PT2 = drawWithCut(t_em, "m_ee", kCutPT2, "h_em_PT2", kNb, edges);
    TH1D* h_pp_PT3 = drawWithCut(t_pp, "m_ee", kCutPT3, "h_pp_PT3", kNb, edges);
    TH1D* h_pp_PT2 = drawWithCut(t_pp, "m_ee", kCutPT2, "h_pp_PT2", kNb, edges);
    TH1D* h_mm_PT3 = drawWithCut(t_mm, "m_ee", kCutPT3, "h_mm_PT3", kNb, edges);
    TH1D* h_mm_PT2 = drawWithCut(t_mm, "m_ee", kCutPT2, "h_mm_PT2", kNb, edges);

    TH1D* h_cb_PT3  = makeCB(h_pp_PT3, h_mm_PT3, "h_cb_PT3");
    TH1D* h_cb_PT2  = makeCB(h_pp_PT2, h_mm_PT2, "h_cb_PT2");
    TH1D* h_pt3_signal_58 = makeSignal(h_em_PT3, h_cb_PT3, "h_pt3_signal_58");
    TH1D* h_pt2_signal_58 = makeSignal(h_em_PT2, h_cb_PT2, "h_pt2_signal_58");

    // 63 * sig_PT2 / sig_PT3 ratio
    TH1D* h_exp_ratio_58 =
        makeRatio(h_pt2_signal_58, h_pt3_signal_58, kTrigCorr, "h_exp_ratio_58");

    // ----- Identify range bins (bin center in [0.15, 0.70]) -----
    std::vector<int> range_bins;
    for (int b = 1; b <= h_exp_ratio_58->GetNbinsX(); ++b) {
        const double xc = h_exp_ratio_58->GetXaxis()->GetBinCenter(b);
        if (xc < kIntLo - 1e-9 || xc > kIntHi + 1e-9) continue;
        range_bins.push_back(b);
    }
    printf("[setup] range bins (center in [%.2f, %.2f]): n = %zu\n",
           kIntLo, kIntHi, range_bins.size());
    if (!range_bins.empty()) {
        printf("[setup] first bin = %d (x_cent = %.4f), last bin = %d (x_cent = %.4f)\n",
               range_bins.front(),
               h_exp_ratio_58->GetXaxis()->GetBinCenter(range_bins.front()),
               range_bins.back(),
               h_exp_ratio_58->GetXaxis()->GetBinCenter(range_bins.back()));
    }

    // =====================================================================
    // VARIANT 1 — Bin-by-bin EXP correction
    // =====================================================================
    printf("\n=========================================================\n");
    printf("VARIANT 1 — Bin-by-bin EXP correction\n");
    printf("=========================================================\n");

    double N_corr_V1   = 0.0;
    double var_V1_sum  = 0.0;
    printf("  %4s  %8s  %12s  %10s  %12s  %10s  %14s  %14s\n",
           "bin", "x_cent",
           "sig_PT3", "s_sig",
           "exp_ratio", "s_r",
           "corrected", "err_bin");
    for (int b : range_bins) {
        const double xc      = h_exp_ratio_58->GetXaxis()->GetBinCenter(b);
        const double sig     = h_pt3_signal_58->GetBinContent(b);
        const double sig_err = h_pt3_signal_58->GetBinError  (b);
        const double r       = h_exp_ratio_58->GetBinContent(b);
        const double r_err   = h_exp_ratio_58->GetBinError  (b);
        const double corr    = sig * r;
        const double err2    = sig * sig * r_err * r_err
                             + r   * r   * sig_err * sig_err;
        const double err_bin = std::sqrt(err2);
        N_corr_V1  += corr;
        var_V1_sum += err2;
        printf("  %4d  %8.4f  %12.3f  %10.3f  %12.4f  %10.4f  %14.3f  %14.3f\n",
               b, xc, sig, sig_err, r, r_err, corr, err_bin);
    }
    const double sigma_V1 = std::sqrt(var_V1_sum);
    const double rel_V1   = (N_corr_V1 != 0.0) ? sigma_V1 / N_corr_V1 : 0.0;
    printf("\n  N_corr_V1 = %.3f\n", N_corr_V1);
    printf("  sigma_V1  = sqrt(sum err^2) = %.3f\n", sigma_V1);
    printf("  Relative  = sigma_V1 / N_corr_V1 = %.6f  (= %.3f %%)\n",
           rel_V1, 100.0 * rel_V1);

    // =====================================================================
    // VARIANT 2 — pol0 fit constant C over [0.15, 0.70]
    // =====================================================================
    printf("\n=========================================================\n");
    printf("VARIANT 2 — Constant C from pol0 fit\n");
    printf("=========================================================\n");

    TF1* f_const = new TF1("f_const", "[0]", kIntLo, kIntHi);
    f_const->SetParameter(0, 1.5);
    h_exp_ratio_58->Fit(f_const, "RQ");
    const double C        = f_const->GetParameter(0);
    const double sigma_C  = f_const->GetParError(0);
    const int    ndf      = f_const->GetNDF();
    const double chi2     = f_const->GetChisquare();
    const double chi2_ndf = (ndf > 0) ? chi2 / ndf : 0.0;

    printf("  C            = %.6f\n", C);
    printf("  sigma_C      = %.6f   (from pol0 fit)\n", sigma_C);
    printf("  chi2 / ndf   = %.3f / %d = %.4f\n", chi2, ndf, chi2_ndf);

    // N_raw_range = sum of h_pt3_signal over range bins.
    double N_raw_range  = 0.0;
    double var_N_raw    = 0.0;
    for (int b : range_bins) {
        N_raw_range += h_pt3_signal_58->GetBinContent(b);
        const double e = h_pt3_signal_58->GetBinError(b);
        var_N_raw += e * e;
    }
    const double sigma_N_raw = std::sqrt(var_N_raw);
    const double N_corr_V2   = N_raw_range * C;
    const double rel_C       = (C != 0.0) ? sigma_C / C : 0.0;
    const double rel_N_raw   = (N_raw_range != 0.0)
                               ? sigma_N_raw / N_raw_range : 0.0;
    const double sigma_V2    = N_corr_V2 *
                               std::sqrt(rel_C * rel_C + rel_N_raw * rel_N_raw);
    const double rel_V2      = (N_corr_V2 != 0.0) ? sigma_V2 / N_corr_V2 : 0.0;

    printf("  N_raw_range  = sum h_pt3_signal over range bins = %.3f\n",
           N_raw_range);
    printf("  sigma_N_raw  = sqrt(sum err^2) = %.3f\n", sigma_N_raw);
    printf("  N_corr_V2    = N_raw_range * C                 = %.3f\n",
           N_corr_V2);
    printf("  sigma_V2     = N * sqrt((sigma_C/C)^2 + (sigma_Nraw/Nraw)^2)\n");
    printf("               = %.3f\n", sigma_V2);
    printf("  Relative     = sigma_V2 / N_corr_V2            = %.6f  (= %.3f %%)\n",
           rel_V2, 100.0 * rel_V2);

    // =====================================================================
    // VARIANT 3 — Three sim shape variants, integral-matched to EXP
    // =====================================================================
    printf("\n=========================================================\n");
    printf("VARIANT 3 — Three sim correction variants\n");
    printf("=========================================================\n");

    const double I_expRatio = histIntegralXBinWidth(h_exp_ratio_58,
                                                    kIntLo, kIntHi);
    // Statistical error of I_exp from per-bin errors of exp_ratio (58-bin).
    // I_exp = Sum ratio(b) * binwidth(b)  over bins fully in [kIntLo, kIntHi]
    // sigma(I_exp)^2 = Sum sigma(ratio(b))^2 * binwidth(b)^2
    double var_I_expRatio = 0.0;
    for (int b = 1; b <= h_exp_ratio_58->GetNbinsX(); ++b) {
        const double lo = h_exp_ratio_58->GetXaxis()->GetBinLowEdge(b);
        const double hi = h_exp_ratio_58->GetXaxis()->GetBinUpEdge (b);
        if (lo < kIntLo - 1e-9) continue;
        if (hi > kIntHi + 1e-9) continue;
        const double er  = h_exp_ratio_58->GetBinError(b);
        const double bw  = h_exp_ratio_58->GetBinWidth(b);
        var_I_expRatio += er * er * bw * bw;
    }
    const double sigma_I_expRatio = std::sqrt(var_I_expRatio);
    const double rel_I_expRatio   = (I_expRatio != 0.0)
                                    ? sigma_I_expRatio / I_expRatio : 0.0;
    printf("  I_expRatio [%.2f, %.2f] = %.6f +/- %.6f  (rel = %.3f %%)\n",
           kIntLo, kIntHi, I_expRatio, sigma_I_expRatio,
           100.0 * rel_I_expRatio);

    // (3a) step4 FIT
    const double I_simFit_step4 = simIntegralStep4(kIntLo, kIntHi);
    const double k_step4_fit    = I_expRatio / I_simFit_step4;
    // (3b) step5 FIT
    const double I_simFit_step5 = simIntegralStep5(kIntLo, kIntHi);
    const double k_step5_fit    = I_expRatio / I_simFit_step5;
    // (3c) step5 POINTS
    TH1D* h_sim_step5 = buildSimFullRatio(t_sim, kSimStep5Cut, "step5_pts");
    if (!h_sim_step5) {
        std::cerr << "Failed to build sim step5 POINTS ratio -- abort.\n";
        return;
    }
    const double I_simPts_step5 =
        histIntegralXBinWidth(h_sim_step5, kIntLo, kIntHi);
    const double k_step5_pts =
        (I_simPts_step5 != 0.0) ? I_expRatio / I_simPts_step5 : 0.0;

    printf("  I_simFit_step4 = %.6f   ->  k_step4_fit = %.6f\n",
           I_simFit_step4, k_step4_fit);
    printf("  I_simFit_step5 = %.6f   ->  k_step5_fit = %.6f\n",
           I_simFit_step5, k_step5_fit);
    printf("  I_simPts_step5 = %.6f   ->  k_step5_pts = %.6f\n",
           I_simPts_step5, k_step5_pts);

    // ---- Per-bin tabulation over range bins ----
    printf("\n  Per-bin correction values c_v(b) and corrected_v(b):\n");
    printf("  %4s  %8s | %12s %12s %12s | %12s %12s %12s\n",
           "bin", "x_cent",
           "c_3a_s4fit", "c_3b_s5fit", "c_3c_s5pts",
           "N_3a", "N_3b", "N_3c");

    double N_V3a = 0.0, N_V3b = 0.0, N_V3c = 0.0;
    for (int b : range_bins) {
        const double xc  = h_exp_ratio_58->GetXaxis()->GetBinCenter(b);
        const double sig = h_pt3_signal_58->GetBinContent(b);
        const double c_a = k_step4_fit * simShapeStep4(xc);
        const double c_b = k_step5_fit * simShapeStep5(xc);
        const double c_c = k_step5_pts * simRatioAt(h_sim_step5, xc);
        const double n_a = sig * c_a;
        const double n_b = sig * c_b;
        const double n_c = sig * c_c;
        N_V3a += n_a;
        N_V3b += n_b;
        N_V3c += n_c;
        printf("  %4d  %8.4f | %12.4f %12.4f %12.4f | %12.3f %12.3f %12.3f\n",
               b, xc, c_a, c_b, c_c, n_a, n_b, n_c);
    }

    const double N_arr[3] = { N_V3a, N_V3b, N_V3c };
    double N_V3_mean = 0.0;
    for (double v : N_arr) N_V3_mean += v;
    N_V3_mean /= 3.0;
    double var_s = 0.0;
    for (double v : N_arr) var_s += (v - N_V3_mean) * (v - N_V3_mean);
    var_s /= 2.0;  // sample stddev: n-1 = 2
    const double sigma_V3_sample = std::sqrt(var_s);
    const double N_V3_min = *std::min_element(N_arr, N_arr + 3);
    const double N_V3_max = *std::max_element(N_arr, N_arr + 3);
    const double spread_V3 = N_V3_max - N_V3_min;
    const double rel_V3    = (N_V3_mean != 0.0)
                             ? sigma_V3_sample / N_V3_mean : 0.0;

    printf("\n  N_V3a (step4 FIT)  = %.3f\n", N_V3a);
    printf("  N_V3b (step5 FIT)  = %.3f\n", N_V3b);
    printf("  N_V3c (step5 PTS)  = %.3f\n", N_V3c);
    printf("  N_V3_mean          = %.3f\n", N_V3_mean);
    printf("  sigma_V3 (sample, n-1) = %.3f   <-- shape\n", sigma_V3_sample);
    printf("  spread_V3 (max-min)    = %.3f\n", spread_V3);
    printf("  Relative (shape) = sigma_V3_shape / N_V3_mean = %.6f  (= %.3f %%)\n",
           rel_V3, 100.0 * rel_V3);

    // ---- Normalization error from sigma(I_exp) ----
    // k = I_exp / I_sim  -> sigma(k)/k = sigma(I_exp)/I_exp (I_sim stats ~ 0)
    // Thus N_corr_v scales by k, so:
    //   sigma_norm(v) = (sigma(I_exp)/I_exp) * N_corr_v
    // Common to all three V3 variants.
    const double sigma_norm_3a   = rel_I_expRatio * N_V3a;
    const double sigma_norm_3b   = rel_I_expRatio * N_V3b;
    const double sigma_norm_3c   = rel_I_expRatio * N_V3c;
    const double sigma_norm_mean = rel_I_expRatio * N_V3_mean;
    const double sigma_V3_shape  = sigma_V3_sample;
    const double sigma_V3_total  = std::sqrt(sigma_V3_shape * sigma_V3_shape +
                                             sigma_norm_mean * sigma_norm_mean);
    const double rel_V3_total    = (N_V3_mean != 0.0)
                                   ? sigma_V3_total / N_V3_mean : 0.0;

    printf("\n  --- Normalization error from sigma(I_exp) ---\n");
    printf("  rel(I_exp)         = %.6f  (= %.3f %%)\n",
           rel_I_expRatio, 100.0 * rel_I_expRatio);
    printf("  sigma_norm(3a)     = %.3f\n", sigma_norm_3a);
    printf("  sigma_norm(3b)     = %.3f\n", sigma_norm_3b);
    printf("  sigma_norm(3c)     = %.3f\n", sigma_norm_3c);
    printf("  sigma_norm(mean)   = %.3f\n", sigma_norm_mean);
    printf("  sigma_V3_shape     = %.3f\n", sigma_V3_shape);
    printf("  sigma_V3_total     = sqrt(shape^2 + norm^2) = %.3f\n",
           sigma_V3_total);
    printf("  Relative (total)   = sigma_V3_total / N_V3_mean = %.6f  (= %.3f %%)\n",
           rel_V3_total, 100.0 * rel_V3_total);

    // =====================================================================
    // COMPARISON TABLE
    // =====================================================================
    printf("\n=========================================================\n");
    printf("COMPARISON TABLE\n");
    printf("=========================================================\n\n");

    printf("| %-22s | %-30s | %12s | %12s | %10s |\n",
           "Variant", "Type", "N_corr", "sigma (abs)", "sigma/N %");
    printf("|%s|%s|%s|%s|%s|\n",
           "------------------------",
           "--------------------------------",
           "--------------",
           "--------------",
           "------------");
    printf("| %-22s | %-30s | %12.3f | %12.3f | %10.3f |\n",
           "V1 EXP bin-by-bin", "data stat propagated",
           N_corr_V1, sigma_V1, 100.0 * rel_V1);
    printf("| %-22s | %-30s | %12.3f | %12.3f | %10.3f |\n",
           "V2 EXP constant fit", "fit uncertainty",
           N_corr_V2, sigma_V2, 100.0 * rel_V2);
    printf("| %-22s | %-30s | %12.3f | %12s | %10s |\n",
           "V3a step4 sim FIT", "(value only)",
           N_V3a, "--", "--");
    printf("| %-22s | %-30s | %12.3f | %12s | %10s |\n",
           "V3b step5 sim FIT", "(value only)",
           N_V3b, "--", "--");
    printf("| %-22s | %-30s | %12.3f | %12s | %10s |\n",
           "V3c step5 sim PTS", "(value only)",
           N_V3c, "--", "--");
    printf("| %-22s | %-30s | %12.3f | %12.3f | %10.3f |\n",
           "V3 shape only", "syst from shape choice",
           N_V3_mean, sigma_V3_shape, 100.0 * rel_V3);
    printf("| %-22s | %-30s | %12.3f | %12.3f | %10.3f |\n",
           "V3 norm only", "from sigma(I_exp)",
           N_V3_mean, sigma_norm_mean, 100.0 * (sigma_norm_mean / N_V3_mean));
    printf("| %-22s | %-30s | %12.3f | %12.3f | %10.3f |\n",
           "V3 TOTAL", "sqrt(shape^2 + norm^2)",
           N_V3_mean, sigma_V3_total, 100.0 * rel_V3_total);

    // ---- Narrative ----
    printf("\nNARRATIVE\n---------\n");
    // Which uncertainty is biggest? (Use V3 TOTAL now.)
    const char* biggest = "V1";
    double biggest_val = sigma_V1;
    if (sigma_V2 > biggest_val) { biggest = "V2"; biggest_val = sigma_V2; }
    if (sigma_V3_total > biggest_val) {
        biggest = "V3"; biggest_val = sigma_V3_total;
    }
    printf("Biggest absolute uncertainty: %s (sigma = %.3f).\n",
           biggest, biggest_val);
    printf("  sigma_V1 = %.3f, sigma_V2 = %.3f, sigma_V3_total = %.3f"
           "  (shape = %.3f, norm = %.3f)\n",
           sigma_V1, sigma_V2, sigma_V3_total,
           sigma_V3_shape, sigma_norm_mean);

    // Consistency: do their central values agree within their errors?
    // Quantify in number of sigma:
    //   nsig = |Na - Nb| / sqrt(sigma_a^2 + sigma_b^2)
    //   <1 sigma -> CONSISTENT
    //   1-2 sig  -> MARGINAL
    //   2-3 sig  -> TENSION
    //   >3 sig   -> DISCREPANT
    auto verdict = [](double nsig) {
        if (nsig < 1.0) return "CONSISTENT";
        if (nsig < 2.0) return "MARGINAL";
        if (nsig < 3.0) return "TENSION";
        return "DISCREPANT";
    };

    const double diff_12 = std::fabs(N_corr_V1 - N_corr_V2);
    const double comb_12 = std::sqrt(sigma_V1 * sigma_V1 +
                                     sigma_V2 * sigma_V2);
    const double nsig_12 = (comb_12 > 0.0) ? diff_12 / comb_12 : 0.0;

    const double diff_13 = std::fabs(N_corr_V1 - N_V3_mean);
    const double comb_13 = std::sqrt(sigma_V1 * sigma_V1 +
                                     sigma_V3_total * sigma_V3_total);
    const double nsig_13 = (comb_13 > 0.0) ? diff_13 / comb_13 : 0.0;

    const double diff_23 = std::fabs(N_corr_V2 - N_V3_mean);
    const double comb_23 = std::sqrt(sigma_V2 * sigma_V2 +
                                     sigma_V3_total * sigma_V3_total);
    const double nsig_23 = (comb_23 > 0.0) ? diff_23 / comb_23 : 0.0;

    printf("Consistency (using sigma_V3_TOTAL = shape (+) norm):\n");
    printf("  V1 vs V2:   |%.3f - %.3f| = %.3f   vs  %.3f"
           "   -> %.3f sigma  [%s]\n",
           N_corr_V1, N_corr_V2, diff_12, comb_12, nsig_12, verdict(nsig_12));
    printf("  V1 vs V3:   |%.3f - %.3f| = %.3f   vs  %.3f"
           "   -> %.3f sigma  [%s]\n",
           N_corr_V1, N_V3_mean, diff_13, comb_13, nsig_13, verdict(nsig_13));
    printf("  V2 vs V3:   |%.3f - %.3f| = %.3f   vs  %.3f"
           "   -> %.3f sigma  [%s]\n",
           N_corr_V2, N_V3_mean, diff_23, comb_23, nsig_23, verdict(nsig_23));

    const bool all_consistent = (nsig_12 < 1.0) && (nsig_13 < 1.0) &&
                                (nsig_23 < 1.0);
    printf("Overall: V1, V2 and V3 are %s within their quoted errors; "
           "the dominant uncertainty is %s.\n",
           all_consistent ? "CONSISTENT" : "NOT all <1 sigma",
           biggest);

    printf("\n=========================================================\n");
    printf("DONE.\n");
    printf("=========================================================\n");

    f_em->Close();
    f_pp->Close();
    f_mm->Close();
    f_sim->Close();
}
