// mass_ee_signal_trigger_corrected_persegment.C
// =========================================================================
// EXP analysis: Fine-binned (280 x 5 MeV) PT3 CB-subtracted m_ee signal,
// compared with TWO trigger-bias correction variants:
//
//   VARIANT A (V2 GLOBAL, reference):
//       LOW  [0, 0.15] : data-driven 63 * sig_PT2 / sig_PT3 sparse-bin ratio
//       HIGH (0.15, 1.4]: single constant C_full = 1.835
//
//   VARIANT B (PER-SEGMENT, new):
//       LOW  [0, 0.15] : per-segment data-driven sparse ratio
//                        corrected(b) = Sum_s [signal_s(b) * ratio_s(b)]
//                        Fallback to full-dataset sparse ratio when
//                        per-segment ratio is undefined.
//       HIGH (0.15, 1.4]: per-segment C_s applied to per-segment signals
//                         then summed: corrected(b) = Sum_s [signal_s(b) * C_s]
//
// C_s values are hard-coded from a previous v2_per_segment run on [0.15, 0.80].
//
// epem trees carry the seg_idx branch. epep/emem have seg_idx==0 only, so the
// CB contribution for HIGH is split across segments by the PT3-fraction
// f_s = n_PT3_s / Sum n_PT3, using the n_PT3 values from pt3_calibration_epem.
//
// Output: plots/output/mass_ee_signal_trigger_corrected_persegment.{pdf,png}
//
// Usage: root -l -b -q plots/mass_ee_signal_trigger_corrected_persegment.C
//
// @author Witold Przygoda (witold.przygoda@uj.edu.pl)
// @date 2026

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
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

constexpr double kTrigCorrLow = 63.0;

// -------- Fine binning: 280 uniform bins x 5 MeV over [0, 1.4] --------
constexpr int    kFineNb   = 280;
constexpr double kFineXmin = 0.0;
constexpr double kFineXmax = 1.4;

// -------- Sparse binning: 58 variable bins --------
constexpr int    kSparseNb   = 58;
constexpr double kSparseXmin = 0.0;
constexpr double kSparseXmax = 1.4;

constexpr double kSplitX = 0.15;

// Per-segment C_s constants from the v2_per_segment run on [0.15, 0.80].
constexpr int    kNSeg = 5;
const double kC_s    [kNSeg] = {2.0488, 1.7511, 1.4749, 1.7185, 1.9976};
const double kC_s_err[kNSeg] = {0.0799, 0.0645, 0.1003, 0.0615, 0.0774};

// V2 global value on [0.15, 0.80].
constexpr double kC_full     = 1.835;
constexpr double kC_full_err = 0.038;

// Per-segment n_PT3 counts (from pt3_calibration_epem.root)
const double kN_PT3[kNSeg] = {5610337.0, 9171792.0, 2363366.0,
                              8105986.0, 7310671.0};

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

// CB = 2 sqrt(N++ * N--) per bin with proper error prop.
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

// CB scaled by a constant factor f (per-segment fraction). Error scales the
// same way (errors on f neglected — it is a fixed n_PT3 ratio).
TH1D* makeCBScaled(TH1D* h_cb_total, double f, const std::string& name) {
    TH1D* cb = static_cast<TH1D*>(h_cb_total->Clone(name.c_str()));
    cb->SetDirectory(nullptr);
    cb->Reset();
    for (int b = 1; b <= h_cb_total->GetNbinsX(); ++b) {
        cb->SetBinContent(b, f * h_cb_total->GetBinContent(b));
        cb->SetBinError  (b, f * h_cb_total->GetBinError  (b));
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
    h->SetMarkerSize(0.85);
    h->SetLineWidth(2);
    h->GetXaxis()->SetTitleSize(0.045);
    h->GetYaxis()->SetTitleSize(0.045);
    h->GetXaxis()->SetLabelSize(0.040);
    h->GetYaxis()->SetLabelSize(0.040);
    h->GetXaxis()->SetTitleOffset(1.05);
    h->GetYaxis()->SetTitleOffset(1.20);
}

double sumByCenter(TH1D* h, double xlo, double xhi) {
    double s = 0.0;
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
        const double xc = h->GetXaxis()->GetBinCenter(b);
        if (xc >= xlo && xc <= xhi) s += h->GetBinContent(b);
    }
    return s;
}

}  // anonymous namespace

void mass_ee_signal_trigger_corrected_persegment() {
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

    std::cout << "Per-segment C_s values (from previous v2_per_segment run):\n";
    for (int s = 0; s < kNSeg; ++s) {
        printf("  C_s[%d] = %.4f +- %.4f   n_PT3 = %.0f\n",
               s, kC_s[s], kC_s_err[s], kN_PT3[s]);
    }
    printf("  C_full (V2 global on [0.15, 0.80]) = %.3f +- %.3f\n",
           kC_full, kC_full_err);

    // Compute per-segment fractions from n_PT3.
    double n_pt3_sum = 0.0;
    for (int s = 0; s < kNSeg; ++s) n_pt3_sum += kN_PT3[s];
    double f_seg[kNSeg];
    printf("\nPer-segment PT3 fractions f_s = n_PT3_s / Sum n_PT3:\n");
    double f_check = 0.0;
    for (int s = 0; s < kNSeg; ++s) {
        f_seg[s] = kN_PT3[s] / n_pt3_sum;
        f_check += f_seg[s];
        printf("  f_seg[%d] = %.6f\n", s, f_seg[s]);
    }
    printf("  Sum f_seg = %.6f (should be 1)\n", f_check);

    gSystem->mkdir("plots/output", true);

    // =====================================================================
    // SECTION A: Per-segment fine-binned PT3 signal spectra
    // =====================================================================
    std::cout << "\n========== SECTION A: per-segment PT3 spectra ==========\n";

    // Full-stat epep/emem PT3 spectra (no segment cut; epep/emem trees
    // have only seg_idx==0 entries but we just want the totals).
    TH1D* h_pt3_epep = drawFine(t_pp, kCutPT3, "h_pt3_epep_full");
    TH1D* h_pt3_emem = drawFine(t_mm, kCutPT3, "h_pt3_emem_full");
    TH1D* h_pt3_cb_full =
        makeCB(h_pt3_epep, h_pt3_emem, "h_pt3_cb_full");

    // Full-stat epem PT3 spectrum (no segment cut).
    TH1D* h_pt3_epem_full = drawFine(t_em, kCutPT3, "h_pt3_epem_full");
    TH1D* h_signal_full   =
        makeSignal(h_pt3_epem_full, h_pt3_cb_full, "h_signal_full");

    printf("Full PT3 raw yields:\n");
    printf("  epem = %.0f\n", h_pt3_epem_full->Integral());
    printf("  epep = %.0f\n", h_pt3_epep->Integral());
    printf("  emem = %.0f\n", h_pt3_emem->Integral());
    printf("  CB   = %.1f\n", h_pt3_cb_full->Integral());
    printf("  sig  = %.1f\n", h_signal_full->Integral());

    // Per-segment epem PT3 spectra and per-segment signals.
    TH1D* h_pt3_epem_s[kNSeg] = {nullptr};
    TH1D* h_pt3_cb_s  [kNSeg] = {nullptr};
    TH1D* h_signal_s  [kNSeg] = {nullptr};
    double sum_signals = 0.0;
    for (int s = 0; s < kNSeg; ++s) {
        const std::string cut =
            std::string(kCutPT3) + " && seg_idx==" + std::to_string(s);
        const std::string name_em = "h_pt3_epem_s" + std::to_string(s);
        const std::string name_cb = "h_pt3_cb_s"   + std::to_string(s);
        const std::string name_sg = "h_signal_s"   + std::to_string(s);
        h_pt3_epem_s[s] = drawFine(t_em, cut.c_str(), name_em);
        h_pt3_cb_s  [s] = makeCBScaled(h_pt3_cb_full, f_seg[s], name_cb);
        h_signal_s  [s] = makeSignal(h_pt3_epem_s[s], h_pt3_cb_s[s], name_sg);
        const double i_em  = h_pt3_epem_s[s]->Integral();
        const double i_cb  = h_pt3_cb_s  [s]->Integral();
        const double i_sig = h_signal_s  [s]->Integral();
        sum_signals += i_sig;
        printf("  seg %d : epem=%.0f  CB(scaled)=%.1f  signal=%.1f\n",
               s, i_em, i_cb, i_sig);
    }
    printf("Sum_s signal_s = %.1f   full signal = %.1f  (should match)\n",
           sum_signals, h_signal_full->Integral());

    // =====================================================================
    // SECTION B: LOW-region trigger correction (sparse 58-bin ratio)
    // =====================================================================
    std::cout << "\n========== SECTION B: LOW correction (sparse ratio) ==========\n";

    // Full-dataset sparse spectra
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

    TH1D* h_sp_ratio = makeRatio(h_sp_sig_PT2, h_sp_sig_PT3, kTrigCorrLow,
                                 "h_sp_ratio");
    printf("Sparse 58-bin ratio (63 * sig_PT2/sig_PT3) built (full dataset).\n");

    // -----------------------------------------------------------------
    // Per-segment sparse ratios.
    // For each segment s, build sparse PT3 and PT2 signals and the
    // 63*sig_PT2/sig_PT3 ratio.  We need per-segment CB fractions for
    // both PT3 and PT2: f_s_PT3 = n_PT3_s / n_PT3_tot (sparse,epem),
    // f_s_PT2 = n_PT2_s / n_PT2_tot (sparse,epem).
    // -----------------------------------------------------------------

    // Sparse epem totals (used as denominators for f_s).
    const double n_sp_em_PT3_tot = h_sp_em_PT3->Integral();
    const double n_sp_em_PT2_tot = h_sp_em_PT2->Integral();

    TH1D* h_sp_em_PT3_s [kNSeg] = {nullptr};
    TH1D* h_sp_em_PT2_s [kNSeg] = {nullptr};
    TH1D* h_sp_cb_PT3_s [kNSeg] = {nullptr};
    TH1D* h_sp_cb_PT2_s [kNSeg] = {nullptr};
    TH1D* h_sp_sig_PT3_s[kNSeg] = {nullptr};
    TH1D* h_sp_sig_PT2_s[kNSeg] = {nullptr};
    TH1D* h_sp_ratio_s  [kNSeg] = {nullptr};

    printf("\nPer-segment sparse ratio build:\n");
    for (int s = 0; s < kNSeg; ++s) {
        const std::string cutPT3 =
            std::string(kCutPT3) + " && seg_idx==" + std::to_string(s);
        const std::string cutPT2 =
            std::string(kCutPT2) + " && seg_idx==" + std::to_string(s);

        h_sp_em_PT3_s[s] = drawSparse(t_em, cutPT3.c_str(),
                                      "h_sp_em_PT3_s" + std::to_string(s));
        h_sp_em_PT2_s[s] = drawSparse(t_em, cutPT2.c_str(),
                                      "h_sp_em_PT2_s" + std::to_string(s));

        const double n_sp_em_PT3_s = h_sp_em_PT3_s[s]->Integral();
        const double n_sp_em_PT2_s = h_sp_em_PT2_s[s]->Integral();
        const double f_s_PT3 =
            (n_sp_em_PT3_tot > 0.0) ? n_sp_em_PT3_s / n_sp_em_PT3_tot : 0.0;
        const double f_s_PT2 =
            (n_sp_em_PT2_tot > 0.0) ? n_sp_em_PT2_s / n_sp_em_PT2_tot : 0.0;

        h_sp_cb_PT3_s[s] = makeCBScaled(h_sp_cb_PT3, f_s_PT3,
                                        "h_sp_cb_PT3_s" + std::to_string(s));
        h_sp_cb_PT2_s[s] = makeCBScaled(h_sp_cb_PT2, f_s_PT2,
                                        "h_sp_cb_PT2_s" + std::to_string(s));

        h_sp_sig_PT3_s[s] = makeSignal(h_sp_em_PT3_s[s], h_sp_cb_PT3_s[s],
                                       "h_sp_sig_PT3_s" + std::to_string(s));
        h_sp_sig_PT2_s[s] = makeSignal(h_sp_em_PT2_s[s], h_sp_cb_PT2_s[s],
                                       "h_sp_sig_PT2_s" + std::to_string(s));

        h_sp_ratio_s[s] = makeRatio(h_sp_sig_PT2_s[s], h_sp_sig_PT3_s[s],
                                    kTrigCorrLow,
                                    "h_sp_ratio_s" + std::to_string(s));

        printf("  seg %d: n_sp_em_PT3=%.0f (f_PT3=%.4f)  "
               "n_sp_em_PT2=%.0f (f_PT2=%.4f)\n",
               s, n_sp_em_PT3_s, f_s_PT3, n_sp_em_PT2_s, f_s_PT2);
    }

    // -----------------------------------------------------------------
    // Per-segment LOW diagnostics: count defined vs undefined sparse
    // bins in [0, 0.15] and compute width-weighted average r_avg_s.
    // -----------------------------------------------------------------
    printf("\nPer-segment LOW [0, 0.15] sparse-bin status & width-weighted avg:\n");
    for (int s = 0; s < kNSeg; ++s) {
        int n_def   = 0;
        int n_undef = 0;
        double w_sum  = 0.0;
        double rw_sum = 0.0;
        for (int sb = 1; sb <= h_sp_ratio_s[s]->GetNbinsX(); ++sb) {
            const double xlo = h_sp_ratio_s[s]->GetXaxis()->GetBinLowEdge(sb);
            const double xup = h_sp_ratio_s[s]->GetXaxis()->GetBinUpEdge(sb);
            // sparse bin entirely or mostly within [0, kSplitX]?
            const double xc = h_sp_ratio_s[s]->GetXaxis()->GetBinCenter(sb);
            if (xc > kSplitX) continue;
            const double n_pt3 = h_sp_sig_PT3_s[s]->GetBinContent(sb);
            const double rc    = h_sp_ratio_s[s]->GetBinContent(sb);
            const bool defined =
                (n_pt3 > 0.0) && std::isfinite(rc) && rc != 0.0;
            if (defined) {
                ++n_def;
                const double bw = xup - xlo;
                w_sum  += bw;
                rw_sum += rc * bw;
            } else {
                ++n_undef;
            }
        }
        const double r_avg = (w_sum > 0.0) ? rw_sum / w_sum : 0.0;
        printf("  seg %d: defined=%d  undefined=%d  total=%d  "
               "r_avg(LOW,width-weighted)=%.4f\n",
               s, n_def, n_undef, n_def + n_undef, r_avg);
    }

    // =====================================================================
    // SECTION C: Build VARIANT A (V2 global) and VARIANT B (per-segment).
    // =====================================================================
    std::cout << "\n========== SECTION C: apply corrections ==========\n";

    TH1D* h_corrected_A = static_cast<TH1D*>(
        h_signal_full->Clone("h_signal_full_V2"));
    h_corrected_A->SetDirectory(nullptr);
    h_corrected_A->Reset();

    TH1D* h_corrected_B = static_cast<TH1D*>(
        h_signal_full->Clone("h_signal_persegment"));
    h_corrected_B->SetDirectory(nullptr);
    h_corrected_B->Reset();

    int n_low_filled_A  = 0;
    int n_low_missing_A = 0;
    int n_low_filled_B  = 0;
    int n_low_missing_B = 0;
    int n_high_filled = 0;
    // Per-segment fallback counters in LOW (B variant)
    int n_fallback_seg[kNSeg] = {0};
    int n_defined_seg [kNSeg] = {0};

    for (int b = 1; b <= h_signal_full->GetNbinsX(); ++b) {
        const double xc = h_signal_full->GetXaxis()->GetBinCenter(b);

        if (xc <= kSplitX) {
            // -------- VARIANT A LOW: full-dataset sparse-bin lookup. --------
            const int sb = h_sp_ratio->FindBin(xc);
            const double rc_full = (sb >= 1 && sb <= h_sp_ratio->GetNbinsX())
                ? h_sp_ratio->GetBinContent(sb) : 0.0;
            const double re_full = (sb >= 1 && sb <= h_sp_ratio->GetNbinsX())
                ? h_sp_ratio->GetBinError(sb) : 0.0;
            const bool full_defined =
                (sb >= 1 && sb <= h_sp_ratio->GetNbinsX()) &&
                std::isfinite(rc_full) && rc_full > 0.0;

            if (full_defined) {
                const double sv  = h_signal_full->GetBinContent(b);
                const double sev = h_signal_full->GetBinError  (b);
                const double val = sv * rc_full;
                const double t1  = rc_full * sev;
                const double t2  = sv      * re_full;
                const double err = std::sqrt(t1 * t1 + t2 * t2);
                h_corrected_A->SetBinContent(b, val);
                h_corrected_A->SetBinError  (b, err);
                ++n_low_filled_A;
            } else {
                ++n_low_missing_A;
            }

            // -------- VARIANT B LOW: per-segment sparse-bin lookup. --------
            double val_B = 0.0;
            double var_B = 0.0;
            bool any_filled = false;
            for (int s = 0; s < kNSeg; ++s) {
                const int sb_s = h_sp_ratio_s[s]->FindBin(xc);
                bool seg_defined = false;
                double rc_s = 0.0, re_s = 0.0;
                if (sb_s >= 1 && sb_s <= h_sp_ratio_s[s]->GetNbinsX()) {
                    const double n_pt3_s =
                        h_sp_sig_PT3_s[s]->GetBinContent(sb_s);
                    rc_s = h_sp_ratio_s[s]->GetBinContent(sb_s);
                    re_s = h_sp_ratio_s[s]->GetBinError  (sb_s);
                    if (n_pt3_s > 0.0 && std::isfinite(rc_s) && rc_s != 0.0) {
                        seg_defined = true;
                    }
                }
                double c_s_use, c_s_err_use;
                if (seg_defined) {
                    c_s_use     = rc_s;
                    c_s_err_use = re_s;
                    ++n_defined_seg[s];
                } else {
                    // fallback to full-dataset sparse ratio
                    c_s_use     = rc_full;
                    c_s_err_use = re_full;
                    ++n_fallback_seg[s];
                }
                if (!std::isfinite(c_s_use) || c_s_use == 0.0) continue;
                const double svs  = h_signal_s[s]->GetBinContent(b);
                const double sevs = h_signal_s[s]->GetBinError  (b);
                val_B += svs * c_s_use;
                const double t1 = c_s_use * sevs;
                const double t2 = svs     * c_s_err_use;
                var_B += t1 * t1 + t2 * t2;
                any_filled = true;
            }
            if (any_filled) {
                h_corrected_B->SetBinContent(b, val_B);
                h_corrected_B->SetBinError  (b, std::sqrt(var_B));
                ++n_low_filled_B;
            } else {
                ++n_low_missing_B;
            }
        } else {
            // HIGH:
            // Variant A: C_full applied to full signal.
            const double sv  = h_signal_full->GetBinContent(b);
            const double sev = h_signal_full->GetBinError  (b);
            const double val_A = sv * kC_full;
            const double t1A = kC_full     * sev;
            const double t2A = sv          * kC_full_err;
            const double err_A = std::sqrt(t1A * t1A + t2A * t2A);
            h_corrected_A->SetBinContent(b, val_A);
            h_corrected_A->SetBinError  (b, err_A);

            // Variant B: sum_s [signal_s(b) * C_s] with proper error prop.
            double val_B = 0.0;
            double var_B = 0.0;
            for (int s = 0; s < kNSeg; ++s) {
                const double svs  = h_signal_s[s]->GetBinContent(b);
                const double sevs = h_signal_s[s]->GetBinError  (b);
                val_B += svs * kC_s[s];
                const double t1 = kC_s[s]     * sevs;
                const double t2 = svs         * kC_s_err[s];
                var_B += t1 * t1 + t2 * t2;
            }
            h_corrected_B->SetBinContent(b, val_B);
            h_corrected_B->SetBinError  (b, std::sqrt(var_B));
            ++n_high_filled;
        }
    }
    printf("Correction applied:\n");
    printf("  V_A LOW  bins filled  : %d   missing : %d\n",
           n_low_filled_A, n_low_missing_A);
    printf("  V_B LOW  bins filled  : %d   missing : %d\n",
           n_low_filled_B, n_low_missing_B);
    printf("  HIGH bins filled      : %d\n", n_high_filled);
    printf("Per-segment LOW fine-bin lookup outcomes "
           "(defined vs fallback to full-dataset ratio):\n");
    for (int s = 0; s < kNSeg; ++s) {
        printf("  seg %d: defined=%d  fallback=%d  total=%d\n",
               s, n_defined_seg[s], n_fallback_seg[s],
               n_defined_seg[s] + n_fallback_seg[s]);
    }

    // =====================================================================
    // SECTION D: Integrals
    // =====================================================================
    std::cout << "\n========== SECTION D: integrals ==========\n";

    const double xlo = kFineXmin + 1e-9;
    const double xhi = kFineXmax - 1e-9;

    const double I_raw_LOW   = sumByCenter(h_signal_full, xlo, kSplitX);
    const double I_raw_HIGH  = sumByCenter(h_signal_full, kSplitX, xhi);
    const double I_raw_tot   = I_raw_LOW + I_raw_HIGH;

    const double I_V2_LOW    = sumByCenter(h_corrected_A, xlo, kSplitX);
    const double I_V2_HIGH   = sumByCenter(h_corrected_A, kSplitX, xhi);
    const double I_V2_tot    = I_V2_LOW + I_V2_HIGH;

    const double I_PS_LOW    = sumByCenter(h_corrected_B, xlo, kSplitX);
    const double I_PS_HIGH   = sumByCenter(h_corrected_B, kSplitX, xhi);
    const double I_PS_tot    = I_PS_LOW + I_PS_HIGH;

    auto dpct = [](double a, double b) {
        return (b != 0.0) ? 100.0 * (a - b) / b : 0.0;
    };

    printf("\n");
    printf("| Region            | I_raw       | I_V2 (global)  | "
           "I_PS (per-seg)  | dPS-V2 abs   | dPS-V2 %%    |\n");
    printf("|-------------------|-------------|----------------|"
           "-----------------|--------------|------------|\n");
    printf("| LOW [0, 0.15]     | %11.2f | %14.2f | %15.2f | %12.4f | %9.4f%% |\n",
           I_raw_LOW,  I_V2_LOW,  I_PS_LOW,
           I_PS_LOW - I_V2_LOW,   dpct(I_PS_LOW,  I_V2_LOW));
    printf("| HIGH (0.15, 1.4]  | %11.2f | %14.2f | %15.2f | %12.4f | %9.4f%% |\n",
           I_raw_HIGH, I_V2_HIGH, I_PS_HIGH,
           I_PS_HIGH - I_V2_HIGH, dpct(I_PS_HIGH, I_V2_HIGH));
    printf("| TOTAL             | %11.2f | %14.2f | %15.2f | %12.4f | %9.4f%% |\n",
           I_raw_tot,  I_V2_tot,  I_PS_tot,
           I_PS_tot - I_V2_tot,   dpct(I_PS_tot,  I_V2_tot));

    printf("\nLOW delta PS-V2 = %.4f  (%.4f%%)\n",
           I_PS_LOW - I_V2_LOW, dpct(I_PS_LOW, I_V2_LOW));

    // =====================================================================
    // SECTION E: Drawing
    // =====================================================================
    std::cout << "\n========== SECTION E: drawing ==========\n";

    TH1D* h_signal_full_w = static_cast<TH1D*>(
        h_signal_full->Clone("h_signal_full_w"));
    h_signal_full_w->SetDirectory(nullptr);
    TH1D* h_corrected_A_w = static_cast<TH1D*>(
        h_corrected_A->Clone("h_signal_full_V2_w"));
    h_corrected_A_w->SetDirectory(nullptr);
    TH1D* h_corrected_B_w = static_cast<TH1D*>(
        h_corrected_B->Clone("h_signal_persegment_w"));
    h_corrected_B_w->SetDirectory(nullptr);
    normalizeByBinWidth(h_signal_full_w);
    normalizeByBinWidth(h_corrected_A_w);
    normalizeByBinWidth(h_corrected_B_w);

    styleHist(h_signal_full_w, kBlack, 20);
    styleHist(h_corrected_A_w, kRed,   22);
    styleHist(h_corrected_B_w, kBlue,  23);

    TCanvas* c = new TCanvas("c_trigcorr_ps",
                             "PT3 signal: raw vs V2 vs per-segment",
                             1100, 700);
    c->SetMargin(0.13, 0.04, 0.13, 0.10);
    gPad->SetLogy(true);

    double y_max = 1.5 * std::max(
        h_signal_full_w->GetMaximum(),
        std::max(h_corrected_A_w->GetMaximum(),
                 h_corrected_B_w->GetMaximum()));
    double y_min_pos = 1e300;
    for (int b = 1; b <= h_signal_full_w->GetNbinsX(); ++b) {
        const double v1 = h_signal_full_w ->GetBinContent(b);
        const double v2 = h_corrected_A_w ->GetBinContent(b);
        const double v3 = h_corrected_B_w ->GetBinContent(b);
        if (v1 > 0.0 && v1 < y_min_pos) y_min_pos = v1;
        if (v2 > 0.0 && v2 < y_min_pos) y_min_pos = v2;
        if (v3 > 0.0 && v3 < y_min_pos) y_min_pos = v3;
    }
    if (!std::isfinite(y_min_pos) || y_min_pos > 1e299) y_min_pos = 0.1;
    double y_min = std::max(0.1, 0.3 * y_min_pos);

    h_signal_full_w->SetTitle(
        "m_{ee} PT3 signal: raw vs V2 global vs per-segment trigger correction;"
        "m_{ee} [GeV/c^{2}];weighted entries / GeV/c^{2}");
    h_signal_full_w->GetXaxis()->SetRangeUser(kFineXmin, kFineXmax);
    h_signal_full_w->GetYaxis()->SetRangeUser(y_min, y_max);
    h_signal_full_w->Draw("E1");
    h_corrected_A_w->Draw("E1 SAME");
    h_corrected_B_w->Draw("E1 SAME");

    TLegend* leg = new TLegend(0.55, 0.70, 0.95, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.030);
    leg->AddEntry(h_signal_full_w, "PT3 signal (raw)",                  "lpe");
    leg->AddEntry(h_corrected_A_w, "V2 global (C_{full} = 1.835)",      "lpe");
    leg->AddEntry(h_corrected_B_w, "per-segment (#Sigma_{s} C_{s}#upoint sig_{s})",
                  "lpe");
    leg->Draw();

    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.028);
    tx.DrawLatex(0.16, 0.86,
        Form("I_{raw,tot} = %.0f,  I_{V2,tot} = %.0f,  "
             "I_{PS,tot} = %.0f", I_raw_tot, I_V2_tot, I_PS_tot));
    tx.DrawLatex(0.16, 0.82,
        Form("LOW: #DeltaPS-V2 = %.2f (%.3f%%)   "
             "HIGH: #DeltaPS-V2 = %.2f (%.3f%%)",
             I_PS_LOW  - I_V2_LOW,  dpct(I_PS_LOW,  I_V2_LOW),
             I_PS_HIGH - I_V2_HIGH, dpct(I_PS_HIGH, I_V2_HIGH)));

    const std::string out_base =
        "plots/output/mass_ee_signal_trigger_corrected_persegment";
    c->SaveAs((out_base + ".pdf").c_str());
    c->SaveAs((out_base + ".png").c_str());
    std::cout << "Saved: " << out_base << ".{pdf,png}\n";

    f_em->Close();
    f_pp->Close();
    f_mm->Close();
}
