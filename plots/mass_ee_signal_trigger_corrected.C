// mass_ee_signal_trigger_corrected.C
// =========================================================================
// EXP analysis: Fine-binned (280 x 5 MeV) PT3 CB-subtracted m_ee signal,
// drawn alongside the SAME spectrum with a trigger-bias correction applied
// bin-by-bin.
//
//   LOW region [0, 0.15]    : per-bin correction = 63 * sig_PT2 / sig_PT3
//                              measured on the variable 58-bin layout from
//                              mass_ee_pt_ratio_signal_varbins_exp.C
//                              (trigger cuts ONLY).
//   HIGH region (0.15, 1.4] : correction = k_orig * simShape(x), where
//                              k_orig is refitted on-the-fly to the 58-bin
//                              signal ratio over [0.10, 1.40] using the
//                              same scaled-SIM piecewise recipe.
//
// Input ntuples (dilepton_nt) already carry the upstream event-level cuts
// isBest==1, eVertReco_z>-500, start_iteration==3 baked in. Therefore the
// only remaining cut applied at TTree::Draw level is the trigger bit. The
// user's nominal cut spec was:
//   cut_pt3 && cut_vertex && cut_isBest && cut_start
// but for dilepton_nt this collapses to cut_pt3 only.
//
// Output: plots/output/mass_ee_signal_trigger_corrected.{pdf,png}
//
// Usage: root -l -b -q plots/mass_ee_signal_trigger_corrected.C
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
#include <TStyle.h>
#include <TSystem.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>

namespace {

constexpr const char* kFileEm = "output_epem_exp.root";
constexpr const char* kFilePp = "output_epep_exp.root";
constexpr const char* kFileMm = "output_emem_exp.root";
constexpr const char* kTree   = "dilepton_nt";

constexpr const char* kCutPT3 = "trigbit==8192";
constexpr const char* kCutPT2 = "trigbit==4096";

constexpr double kTrigCorr = 63.0;

// -------- Fine binning: 280 uniform bins x 5 MeV over [0, 1.4] --------
constexpr int    kFineNb   = 280;
constexpr double kFineXmin = 0.0;
constexpr double kFineXmax = 1.4;

// -------- Sparse binning: 58 variable bins (same as ratio macro) --------
constexpr int    kSparseNb   = 58;
constexpr double kSparseXmin = 0.0;
constexpr double kSparseXmax = 1.4;

// LOW/HIGH split.
constexpr double kSplitX = 0.15;

// SCALED SIM RECIPE constants (sim cascade v5 step4 REC).
constexpr double SIM_JOIN_X  = 0.70;
constexpr double SIM_HI_A    = 1.080;
constexpr double SIM_HI_B    = 0.039;
constexpr double SIM_HI_C    = 0.117;
constexpr double SIM_MID_B   = 0.264;
constexpr double SIM_MID_C   = 0.250;
const double SIM_JOIN_Y =
    SIM_HI_A + SIM_HI_B * SIM_JOIN_X + SIM_HI_C * SIM_JOIN_X * SIM_JOIN_X;

// Piecewise sim shape (no scale factor).
double simShape(double x) {
    if (x < SIM_JOIN_X) {
        const double d = x - SIM_JOIN_X;
        return SIM_JOIN_Y + SIM_MID_B * d + SIM_MID_C * d * d;
    }
    return SIM_HI_A + SIM_HI_B * x + SIM_HI_C * x * x;
}

// TF1-compatible: y = k * simShape(x), where k = par[0].
double simShapeScaled(double* x, double* par) {
    return par[0] * simShape(x[0]);
}

// Build the 58-bin variable-width bin edge array (matches ratio macro).
const double* getSparseEdges() {
    static double edges[kSparseNb + 1];
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

// Draw m_ee with the given cut into a histogram of the specified layout.
TH1D* drawFine(TTree* t, const char* cut, const std::string& name) {
    TH1D* h = new TH1D(name.c_str(), "", kFineNb, kFineXmin, kFineXmax);
    h->Sumw2();
    t->Draw(("m_ee>>" + name).c_str(), cut, "goff");
    h->SetDirectory(nullptr);
    return h;
}
TH1D* drawSparse(TTree* t, const char* cut, const std::string& name) {
    TH1D* h = new TH1D(name.c_str(), "", kSparseNb, getSparseEdges());
    h->Sumw2();
    t->Draw(("m_ee>>" + name).c_str(), cut, "goff");
    h->SetDirectory(nullptr);
    return h;
}

// CB = 2*sqrt(N++ * N--) per bin with proper error propagation.
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

// Signal = epem - CB with quadrature error propagation.
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

// Bin-by-bin ratio = factor * num / den with relative-error propagation.
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

// Divide each bin (content & error) by its bin width.
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
    h->SetMarkerSize(0.7);
    h->SetLineWidth(2);
    h->GetXaxis()->SetTitleSize(0.045);
    h->GetYaxis()->SetTitleSize(0.045);
    h->GetXaxis()->SetLabelSize(0.040);
    h->GetYaxis()->SetLabelSize(0.040);
    h->GetXaxis()->SetTitleOffset(1.05);
    h->GetYaxis()->SetTitleOffset(1.20);
}

// Sum bin contents over [xlo, xhi] (using bin centers to decide inclusion).
double sumByCenter(TH1D* h, double xlo, double xhi) {
    double s = 0.0;
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
        const double xc = h->GetXaxis()->GetBinCenter(b);
        if (xc >= xlo && xc <= xhi) s += h->GetBinContent(b);
    }
    return s;
}

}  // anonymous namespace

void mass_ee_signal_trigger_corrected() {
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
    std::cout << "  epem : " << kFileEm
              << " (" << t_em->GetEntries() << " entries)\n";
    std::cout << "  epep : " << kFilePp
              << " (" << t_pp->GetEntries() << " entries)\n";
    std::cout << "  emem : " << kFileMm
              << " (" << t_mm->GetEntries() << " entries)\n";

    // ---------- The user's nominal cut spec ----------
    // cut_pt3 && cut_vertex && cut_isBest && cut_start
    // For dilepton_nt, vertex/isBest/start are upstream-applied; only trigbit
    // is selectable here.
    std::cout << "Cuts applied at TTree::Draw level:\n";
    std::cout << "  PT3: '" << kCutPT3 << "'  (vertex/isBest/start are "
                 "already pre-filtered upstream)\n";
    std::cout << "  PT2: '" << kCutPT2 << "'\n\n";

    gSystem->mkdir("plots/output", true);

    // =====================================================================
    // SECTION A: Fine-binned PT3 signal mass spectrum (280 x 5 MeV)
    // =====================================================================
    std::cout << "============ SECTION A: fine-binned PT3 signal ============\n";

    TH1D* h_pt3_epem = drawFine(t_em, kCutPT3, "h_pt3_epem");
    TH1D* h_pt3_epep = drawFine(t_pp, kCutPT3, "h_pt3_epep");
    TH1D* h_pt3_emem = drawFine(t_mm, kCutPT3, "h_pt3_emem");

    TH1D* h_pt3_cb     = makeCB(h_pt3_epep, h_pt3_emem, "h_pt3_cb");
    TH1D* h_pt3_signal = makeSignal(h_pt3_epem, h_pt3_cb, "h_pt3_signal");

    std::cout << "Fine PT3 raw yields:\n";
    std::cout << "  epem = " << h_pt3_epem->Integral() << "\n";
    std::cout << "  epep = " << h_pt3_epep->Integral() << "\n";
    std::cout << "  emem = " << h_pt3_emem->Integral() << "\n";
    std::cout << "  CB     = " << h_pt3_cb->Integral() << "\n";
    std::cout << "  signal = " << h_pt3_signal->Integral() << "\n";

    // =====================================================================
    // SECTION B: Trigger-bias correction lookup (sparse 58-bin ratio)
    // =====================================================================
    std::cout << "\n========== SECTION B: trigger correction lookup ==========\n";

    TH1D* h_sp_em_PT3 = drawSparse(t_em, kCutPT3, "h_sp_em_PT3");
    TH1D* h_sp_em_PT2 = drawSparse(t_em, kCutPT2, "h_sp_em_PT2");
    TH1D* h_sp_pp_PT3 = drawSparse(t_pp, kCutPT3, "h_sp_pp_PT3");
    TH1D* h_sp_pp_PT2 = drawSparse(t_pp, kCutPT2, "h_sp_pp_PT2");
    TH1D* h_sp_mm_PT3 = drawSparse(t_mm, kCutPT3, "h_sp_mm_PT3");
    TH1D* h_sp_mm_PT2 = drawSparse(t_mm, kCutPT2, "h_sp_mm_PT2");

    TH1D* h_sp_cb_PT3  = makeCB(h_sp_pp_PT3, h_sp_mm_PT3, "h_sp_cb_PT3");
    TH1D* h_sp_cb_PT2  = makeCB(h_sp_pp_PT2, h_sp_mm_PT2, "h_sp_cb_PT2");
    TH1D* h_sp_sig_PT3 = makeSignal(h_sp_em_PT3, h_sp_cb_PT3, "h_sp_sig_PT3");
    TH1D* h_sp_sig_PT2 = makeSignal(h_sp_em_PT2, h_sp_cb_PT2, "h_sp_sig_PT2");

    // ratio(b) = 63 * sig_PT2(b) / sig_PT3(b)
    TH1D* h_sp_ratio = makeRatio(h_sp_sig_PT2, h_sp_sig_PT3, kTrigCorr,
                                 "h_sp_ratio");

    // ---------- Refit k_orig with the scaled SIM recipe over [0.10, 1.40] ----------
    TF1* f_simref = new TF1("f_simref_local",
                            simShapeScaled, 0.10, 1.40, 1);
    f_simref->SetParameter(0, 2.0);
    f_simref->SetParName(0, "k");
    h_sp_ratio->Fit(f_simref, "RQ0");
    const double k_orig    = f_simref->GetParameter(0);
    const double k_orig_err= f_simref->GetParError(0);
    const int    sim_ndf   = f_simref->GetNDF();
    const double sim_chi2  = f_simref->GetChisquare();
    const double sim_cndf  = (sim_ndf > 0) ? sim_chi2 / sim_ndf : 0.0;

    printf("Refitted k_orig from 58-bin signal ratio over [0.10, 1.40]:\n");
    printf("  k_orig    = %.4f +- %.4f\n", k_orig, k_orig_err);
    printf("  chi2/ndf  = %.2f / %d = %.3f\n", sim_chi2, sim_ndf, sim_cndf);

    // =====================================================================
    // SECTION C: Apply correction to fine-binned PT3 signal
    // =====================================================================
    std::cout << "\n========== SECTION C: apply correction ==========\n";

    TH1D* h_pt3_signal_corrected = static_cast<TH1D*>(
        h_pt3_signal->Clone("h_pt3_signal_corrected"));
    h_pt3_signal_corrected->SetDirectory(nullptr);
    h_pt3_signal_corrected->Reset();

    int n_low_filled  = 0;
    int n_low_missing = 0;
    int n_high_filled = 0;
    printf("\nFINE-BIN LOW CORRECTION DUMP (step4):\n");
    printf("   bin   x_center     correction_factor   raw_signal\n");
    printf("   ---   ---------    -----------------   -----------\n");
    for (int b = 1; b <= h_pt3_signal->GetNbinsX(); ++b) {
        const double xc  = h_pt3_signal->GetXaxis()->GetBinCenter(b);
        const double sv  = h_pt3_signal->GetBinContent(b);
        const double sev = h_pt3_signal->GetBinError  (b);
        double corr = 0.0;
        double corr_err = 0.0;
        if (xc <= kSplitX) {
            // LOW: look up the sparse measurement bin containing xc.
            const int sb = h_sp_ratio->FindBin(xc);
            if (sb < 1 || sb > h_sp_ratio->GetNbinsX()) {
                h_pt3_signal_corrected->SetBinContent(b, 0.0);
                h_pt3_signal_corrected->SetBinError  (b, 0.0);
                ++n_low_missing;
                continue;
            }
            const double rc = h_sp_ratio->GetBinContent(sb);
            const double re = h_sp_ratio->GetBinError  (sb);
            if (rc <= 0.0 || !std::isfinite(rc)) {
                h_pt3_signal_corrected->SetBinContent(b, 0.0);
                h_pt3_signal_corrected->SetBinError  (b, 0.0);
                ++n_low_missing;
                continue;
            }
            corr     = rc;
            corr_err = re;
            ++n_low_filled;
        } else {
            // HIGH: scaledSimFit(x) = k_orig * simShape(x).
            corr     = k_orig * simShape(xc);
            corr_err = 0.0;   // fit-side error negligible vs sig stat error
            ++n_high_filled;
        }
        if (b <= 30) {
            printf("  %4d   %.6f   %.15f   %.6f\n", b, xc, corr, sv);
        }
        const double val = sv * corr;
        // err = sqrt((corr*sev)^2 + (sv*corr_err)^2)
        const double t1  = corr * sev;
        const double t2  = sv   * corr_err;
        const double err = std::sqrt(t1 * t1 + t2 * t2);
        h_pt3_signal_corrected->SetBinContent(b, val);
        h_pt3_signal_corrected->SetBinError  (b, err);
    }
    printf("Correction applied:\n");
    printf("  LOW  bins filled  : %d\n", n_low_filled);
    printf("  LOW  bins missing : %d (sparse ratio undefined - zeroed)\n",
           n_low_missing);
    printf("  HIGH bins filled  : %d\n", n_high_filled);

    // =====================================================================
    // SECTION E: Integrals
    // =====================================================================
    std::cout << "\n========== SECTION E: integrals ==========\n";
    const double I_raw       = h_pt3_signal->Integral();
    const double I_corrected = h_pt3_signal_corrected->Integral();
    const double delta       = I_corrected - I_raw;
    const double ratio       = (I_raw != 0.0) ? I_corrected / I_raw : 0.0;
    printf("I_raw       = %.3f\n", I_raw);
    printf("I_corrected = %.3f\n", I_corrected);
    printf("delta       = I_corrected - I_raw = %.3f\n", delta);
    if (ratio > 1.0) {
        printf("ratio       = I_corrected / I_raw = %.3f  -> corrected is "
               "%.1f%% higher\n",
               ratio, 100.0 * (ratio - 1.0));
    } else if (ratio < 1.0 && ratio > 0.0) {
        printf("ratio       = I_corrected / I_raw = %.3f  -> corrected is "
               "%.1f%% lower\n",
               ratio, 100.0 * (1.0 - ratio));
    } else {
        printf("ratio       = I_corrected / I_raw = %.3f\n", ratio);
    }

    // LOW/HIGH split (use a small epsilon to avoid double-counting the
    // boundary bin; the fine binning has a bin edge at exactly 0.15).
    const double xlo = kFineXmin + 1e-9;
    const double xhi = kFineXmax - 1e-9;
    const double I_raw_LOW  = sumByCenter(h_pt3_signal, xlo, kSplitX);
    const double I_raw_HIGH = sumByCenter(h_pt3_signal, kSplitX, xhi);
    const double I_cor_LOW  = sumByCenter(h_pt3_signal_corrected, xlo, kSplitX);
    const double I_cor_HIGH = sumByCenter(h_pt3_signal_corrected, kSplitX, xhi);
    printf("\nSplit at %.2f:\n", kSplitX);
    printf("  I_raw_LOW  = %.3f   I_cor_LOW  = %.3f   (delta %+.3f, ratio %.3f)\n",
           I_raw_LOW, I_cor_LOW, I_cor_LOW - I_raw_LOW,
           (I_raw_LOW != 0.0 ? I_cor_LOW / I_raw_LOW : 0.0));
    printf("  I_raw_HIGH = %.3f   I_cor_HIGH = %.3f   (delta %+.3f, ratio %.3f)\n",
           I_raw_HIGH, I_cor_HIGH, I_cor_HIGH - I_raw_HIGH,
           (I_raw_HIGH != 0.0 ? I_cor_HIGH / I_raw_HIGH : 0.0));

    // =====================================================================
    // SECTION D: Drawing
    // =====================================================================
    std::cout << "\n========== SECTION D: drawing ==========\n";

    // Width-normalize clones for display.
    TH1D* h_pt3_signal_w_norm = static_cast<TH1D*>(
        h_pt3_signal->Clone("h_pt3_signal_w_norm"));
    h_pt3_signal_w_norm->SetDirectory(nullptr);
    TH1D* h_pt3_signal_corrected_w_norm = static_cast<TH1D*>(
        h_pt3_signal_corrected->Clone("h_pt3_signal_corrected_w_norm"));
    h_pt3_signal_corrected_w_norm->SetDirectory(nullptr);
    normalizeByBinWidth(h_pt3_signal_w_norm);
    normalizeByBinWidth(h_pt3_signal_corrected_w_norm);

    styleHist(h_pt3_signal_w_norm,           kBlack, 20);
    styleHist(h_pt3_signal_corrected_w_norm, kRed,   20);

    TCanvas* c = new TCanvas("c_trigcorr",
                             "PT3 signal raw vs trigger-bias corrected",
                             1100, 700);
    c->SetMargin(0.13, 0.04, 0.13, 0.10);
    gPad->SetLogy(true);

    // Y range: use max of corrected as upper, smallest positive as lower.
    double y_max_raw = h_pt3_signal_w_norm->GetMaximum();
    double y_max_cor = h_pt3_signal_corrected_w_norm->GetMaximum();
    double y_max     = 1.5 * std::max(y_max_raw, y_max_cor);
    double y_min_pos = 1e300;
    for (int b = 1; b <= h_pt3_signal_w_norm->GetNbinsX(); ++b) {
        const double v1 = h_pt3_signal_w_norm->GetBinContent(b);
        const double v2 = h_pt3_signal_corrected_w_norm->GetBinContent(b);
        if (v1 > 0.0 && v1 < y_min_pos) y_min_pos = v1;
        if (v2 > 0.0 && v2 < y_min_pos) y_min_pos = v2;
    }
    if (!std::isfinite(y_min_pos) || y_min_pos > 1e299) y_min_pos = 0.1;
    double y_min = std::max(0.1, 0.3 * y_min_pos);

    h_pt3_signal_w_norm->SetTitle(
        "m_{ee} PT3 signal - raw vs trigger-bias corrected;"
        "m_{ee} [GeV/c^{2}];weighted entries / GeV/c^{2}");
    h_pt3_signal_w_norm->GetXaxis()->SetRangeUser(kFineXmin, kFineXmax);
    h_pt3_signal_w_norm->GetYaxis()->SetRangeUser(y_min, y_max);
    h_pt3_signal_w_norm->Draw("E1");
    h_pt3_signal_corrected_w_norm->Draw("E1 SAME");

    TLegend* leg = new TLegend(0.55, 0.72, 0.95, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.030);
    leg->AddEntry(h_pt3_signal_w_norm,
                  "PT3 signal (raw)", "lpe");
    leg->AddEntry(h_pt3_signal_corrected_w_norm,
                  "PT3 trigger bias corrected", "lpe");
    leg->Draw();

    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.028);
    tx.DrawLatex(0.16, 0.86,
        Form("k_{orig} = %.3f #pm %.3f  (#chi^{2}/ndf = %.2f / %d)",
             k_orig, k_orig_err, sim_chi2, sim_ndf));
    tx.DrawLatex(0.16, 0.82,
        Form("I_{raw} = %.0f,  I_{cor} = %.0f,  ratio = %.3f",
             I_raw, I_corrected, ratio));

    // ---------- Save ----------
    const std::string out_base = "plots/output/mass_ee_signal_trigger_corrected";
    c->SaveAs((out_base + ".pdf").c_str());
    c->SaveAs((out_base + ".png").c_str());
    std::cout << "Saved: " << out_base << ".{pdf,png}\n";

    f_em->Close();
    f_pp->Close();
    f_mm->Close();
}
