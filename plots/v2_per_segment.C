// plots/v2_per_segment.C
// =========================================================================
// Per-segment V2 (PT2/PT3 trigger correction) analysis.
//
// Background
// ----------
// The previous V2 analysis fitted a single global pol0 constant C to the
// ratio   63 · N_PT2^sig / N_PT3^sig   over [0.15, 0.80] on the full epem
// dataset (one number). The reprocessed output_*_exp.root files now carry
// a `seg_idx` field (and `w_seg`) per row, partitioning the epem statistics
// into 5 segments (0..4). Same-sign channels (epep, emem) are not segmented
// here: each carries 1 segment (seg_idx == 0).
//
// What this macro does
// --------------------
// 1. For each EPEM segment s in 0..4:
//      a. Build h_epem_s_PT3 = m_ee histogram with trigbit==8192 && seg_idx==s
//      b. Build h_epem_s_PT2 = m_ee histogram with trigbit==4096 && seg_idx==s
//      c. Same 58-bin variable binning over [0.0, 1.4] used by other
//         V2 macros (20×2.5 MeV [0, 0.05] + 38×~35.526 MeV [0.05, 1.4]).
//
// 2. Combinatorial background uses the TOTAL same-sign yields (epep/emem
//    have only 1 segment). To attribute CB to each EPEM segment we use a
//    pragmatic per-bin fraction f_s = n_PT3_seg_s / sum_s' n_PT3_seg_s',
//    which scales the magnitude only (the shape is assumed roughly stable
//    across segments). Then:
//        CB_s_PT3(b) = f_s * 2 * sqrt(N_epep_PT3(b) * N_emem_PT3(b))
//        CB_s_PT2(b) = f_s * 2 * sqrt(N_epep_PT2(b) * N_emem_PT2(b))
//        h_sig_s_PT3 = h_epem_s_PT3 - CB_s_PT3
//        h_sig_s_PT2 = h_epem_s_PT2 - CB_s_PT2
//    (errors propagated bin-by-bin)
//
// 3. Per-segment ratio:
//        h_ratio_s = 63 * h_sig_s_PT2 / h_sig_s_PT3
//    Fit pol0 over [0.15, 0.80]  ->  C_s ± σ_C_s, χ²/ndf.
//
// 4. PT3-count-weighted combination across segments:
//        C_combined = Σ_s (n_PT3_s * C_s) / Σ_s n_PT3_s
//        σ²_C_comb  = Σ_s (n_PT3_s / Σ n_PT3)² * σ_C_s²
//    n_PT3 values are read from pt3_calibration_epem.root (one entry per
//    seg_idx; the calibration TTree has redundant rows so we keep the
//    unique seg_idx -> n_pt3 map).
//
// 5. Compare to the V2 extended pol0 result on the full dataset:
//        C_full = 1.835 ± 0.038  on [0.15, 0.80]
//        (χ²/ndf = 0.77 on 36 bins × 18.06 MeV, from
//         v2_flatness_subwindow_test).
//
// 6. Apply to integrated yield:
//        N_corr_segmented = Σ_s [ Σ_bin h_sig_s_PT3(b) * C_s ]  on [0.15, 0.80]
//    Compared to N_corr_V2_full = 570,323 ± 11,748 (V2 extended).
//
// Output
// ------
//   plots/output/v2_per_segment.{pdf,png}
//     - top: 5 ratio histograms overlaid (one per segment)
//     - per-segment horizontal pol0 lines with error bands
//     - annotations with C_combined and comparison to C_full.
//
// Pre-requisites
// --------------
// REQUIRES reprocessed output_*_exp.root with seg_idx field. To regenerate:
//     cd /home/przygoda/HADES/FAT/FAT
//     make
//     ./run_parallel.sh config_epem.json
//     ./run_parallel.sh config_epep.json
//     ./run_parallel.sh config_emem.json
//
// Usage
// -----
//     root -l -b -q plots/v2_per_segment.C
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
#include <TBox.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TColor.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace {

constexpr const char* kFileEm = "output_epem_exp.root";
constexpr const char* kFilePp = "output_epep_exp.root";
constexpr const char* kFileMm = "output_emem_exp.root";
constexpr const char* kTree   = "dilepton_nt";

constexpr const char* kFilePt3Cal = "pt3_calibration_epem.root";
constexpr const char* kTreePt3Cal = "pt3_calibration";

// 58-bin variable binning over [0.0, 1.4].
constexpr int    kNb   = 58;
constexpr double kXmin = 0.0;
constexpr double kXmax = 1.4;

// Trigger downscale: PT2 is downscaled by 64, so we multiply by 63.
constexpr double kTrigCorr = 63.0;

// Trigger bit cuts.
constexpr const char* kCutPT3 = "trigbit==8192";
constexpr const char* kCutPT2 = "trigbit==4096";

// Fit range for pol0 V2 fit.
constexpr double kFitLo = 0.15;
constexpr double kFitHi = 0.80;

// Number of EPEM segments expected.
constexpr int kNSeg = 5;

// Reference global V2 result on [0.15, 0.80] (from
// v2_flatness_subwindow_test extended pol0 fit).
constexpr double kC_full      = 1.835;
constexpr double kSC_full     = 0.038;
// Note: the previous hard-coded reference yield (570,323) was incorrect; it
// was N_raw_signal in the OLD range [0.15, 0.70] (= 310,886) times C_full.
// The proper reference is computed at runtime as
//   N_raw_signal_FULL_015_080 * C_full
// (see kFitLo, kFitHi and the FULL signal histogram built below).

// Build the variable-width bin edge array (58 bins total).
//   [0.00, 0.05] : 20 bins of 2.5 MeV
//   [0.05, 1.40] : 38 bins of ~35.526 MeV
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

// Compute CB = 2*sqrt(N_++ · N_--) per bin with relative error propagation.
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

// Scale a CB histogram by a scalar fraction (per-segment fraction f_s),
// preserving relative errors. Returns a new histogram.
TH1D* scaleCB(TH1D* h_cb, double frac, const std::string& name) {
    TH1D* h = static_cast<TH1D*>(h_cb->Clone(name.c_str()));
    h->SetDirectory(nullptr);
    h->Reset();
    for (int b = 1; b <= h_cb->GetNbinsX(); ++b) {
        const double v = h_cb->GetBinContent(b);
        const double e = h_cb->GetBinError(b);
        h->SetBinContent(b, frac * v);
        h->SetBinError  (b, std::abs(frac) * e);
    }
    return h;
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

// Bin-by-bin ratio = factor · num / den with full relative-error propagation.
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

// Pol0 fit over [xlo, xhi]; returns C, σ_C, χ²/ndf, ndf.
struct Pol0Result {
    double C;
    double sC;
    double chi2;
    int    ndf;
    double chi2_ndf;
};
Pol0Result fitPol0(TH1D* h, double xlo, double xhi, const std::string& tag) {
    const std::string fname = std::string(h->GetName()) + "_pol0_" + tag;
    TF1* f = new TF1(fname.c_str(), "[0]", xlo, xhi);
    f->SetParameter(0, 2.0);
    h->Fit(f, "RQ0");
    Pol0Result r;
    r.C       = f->GetParameter(0);
    r.sC      = f->GetParError(0);
    r.chi2    = f->GetChisquare();
    r.ndf     = f->GetNDF();
    r.chi2_ndf = (r.ndf > 0) ? r.chi2 / r.ndf : 0.0;
    delete f;
    return r;
}

// Sum h_sig over bins fully inside [xlo, xhi] -> (sum, σ_sum).
void sumWithError(TH1D* h, double xlo, double xhi,
                  double& sum, double& sum_err) {
    sum     = 0.0;
    double var = 0.0;
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
        const double lo = h->GetXaxis()->GetBinLowEdge(b);
        const double hi = h->GetXaxis()->GetBinUpEdge(b);
        if (lo >= xlo - 1e-9 && hi <= xhi + 1e-9) {
            sum += h->GetBinContent(b);
            const double e = h->GetBinError(b);
            var += e * e;
        }
    }
    sum_err = std::sqrt(var);
}

// Read pt3_calibration_epem.root: seg_idx -> n_pt3 (unique map).
// The calibration tree has one row per file with seg_idx repeated; we
// fold to the unique (seg_idx -> n_pt3) value (n_pt3 is constant within a
// segment per construction).
std::map<int, long long> readPT3Calibration() {
    std::map<int, long long> out;
    TFile* f = TFile::Open(kFilePt3Cal, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "[pt3 cal] WARNING: cannot open " << kFilePt3Cal << "\n";
        if (f) delete f;
        return out;
    }
    TTree* t = dynamic_cast<TTree*>(f->Get(kTreePt3Cal));
    if (!t) {
        std::cerr << "[pt3 cal] WARNING: tree '" << kTreePt3Cal
                  << "' not found in " << kFilePt3Cal << "\n";
        f->Close(); delete f;
        return out;
    }
    Int_t     seg_idx = -1;
    Long64_t  n_pt3   = 0;
    t->SetBranchAddress("seg_idx", &seg_idx);
    t->SetBranchAddress("n_pt3",   &n_pt3);
    const Long64_t N = t->GetEntries();
    for (Long64_t i = 0; i < N; ++i) {
        t->GetEntry(i);
        if (seg_idx < 0) continue;
        // Keep first value seen per seg_idx (they are constant per segment).
        if (out.find(seg_idx) == out.end()) {
            out[seg_idx] = n_pt3;
        }
    }
    f->Close(); delete f;
    return out;
}

// Style helper for hist markers / lines.
void styleHist(TH1D* h, Color_t c, Style_t m) {
    h->SetLineColor(c);
    h->SetMarkerColor(c);
    h->SetMarkerStyle(m);
    h->SetMarkerSize(0.9);
    h->SetLineWidth(2);
    h->GetXaxis()->SetTitleSize(0.048);
    h->GetYaxis()->SetTitleSize(0.048);
    h->GetXaxis()->SetLabelSize(0.042);
    h->GetYaxis()->SetLabelSize(0.042);
    h->GetXaxis()->SetTitleOffset(1.05);
    h->GetYaxis()->SetTitleOffset(1.15);
}

}  // namespace

void v2_per_segment() {
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
        std::cerr << "  Expected: " << kFileEm << ", " << kFilePp
                  << ", " << kFileMm << "\n";
        std::cerr << "  Re-run the analysis with the new seg_idx field first.\n";
        return;
    }
    TTree* t_em = dynamic_cast<TTree*>(f_em->Get(kTree));
    TTree* t_pp = dynamic_cast<TTree*>(f_pp->Get(kTree));
    TTree* t_mm = dynamic_cast<TTree*>(f_mm->Get(kTree));
    if (!t_em || !t_pp || !t_mm) {
        std::cerr << "Tree " << kTree
                  << " missing in one of the EXP files.\n";
        return;
    }

    // Verify seg_idx branch exists on the EPEM tree.
    if (!t_em->GetBranch("seg_idx")) {
        std::cerr << "ERROR: 'seg_idx' branch missing on " << kFileEm
                  << ":" << kTree << ".\n"
                  << "  Reprocess the EXP files with the new pipeline:\n"
                  << "      make\n"
                  << "      ./run_parallel.sh config_epem.json\n"
                  << "      ./run_parallel.sh config_epep.json\n"
                  << "      ./run_parallel.sh config_emem.json\n";
        f_em->Close(); f_pp->Close(); f_mm->Close();
        return;
    }

    std::cout << "Input:\n";
    std::cout << "  epem : " << kFileEm
              << " (" << t_em->GetEntries() << " entries)\n";
    std::cout << "  epep : " << kFilePp
              << " (" << t_pp->GetEntries() << " entries)\n";
    std::cout << "  emem : " << kFileMm
              << " (" << t_mm->GetEntries() << " entries)\n";
    std::cout << "Binning: 58 variable bins on [0.0, 1.4]\n";
    std::cout << "Trigger cuts: PT3='" << kCutPT3
              << "'  PT2='" << kCutPT2 << "'\n";
    std::cout << "Trigger correction factor: " << kTrigCorr << "\n";
    std::cout << "V2 fit range: [" << kFitLo << ", " << kFitHi << "]\n\n";

    gSystem->mkdir("plots/output", true);

    // ---------- Read n_PT3 per segment from calibration file ----------
    std::map<int, long long> seg_n_pt3 = readPT3Calibration();
    if (seg_n_pt3.empty()) {
        std::cerr << "ERROR: no segments read from " << kFilePt3Cal << "\n";
        f_em->Close(); f_pp->Close(); f_mm->Close();
        return;
    }
    long long sum_n_pt3 = 0;
    for (const auto& kv : seg_n_pt3) sum_n_pt3 += kv.second;
    std::cout << "PT3 calibration (n_pt3 per segment):\n";
    for (const auto& kv : seg_n_pt3) {
        std::cout << "  seg " << kv.first << " : n_pt3 = " << kv.second << "\n";
    }
    std::cout << "  total n_pt3 across segments = " << sum_n_pt3 << "\n\n";

    // ---------- Build same-sign TOTAL CB (no segment filter) ----------
    TH1D* h_pp_PT3 = drawWithCut(t_pp, "m_ee", kCutPT3, "h_pp_PT3");
    TH1D* h_pp_PT2 = drawWithCut(t_pp, "m_ee", kCutPT2, "h_pp_PT2");
    TH1D* h_mm_PT3 = drawWithCut(t_mm, "m_ee", kCutPT3, "h_mm_PT3");
    TH1D* h_mm_PT2 = drawWithCut(t_mm, "m_ee", kCutPT2, "h_mm_PT2");

    TH1D* h_cb_PT3_total = makeCB(h_pp_PT3, h_mm_PT3, "h_cb_PT3_total");
    TH1D* h_cb_PT2_total = makeCB(h_pp_PT2, h_mm_PT2, "h_cb_PT2_total");

    std::cout << "Same-sign totals (all segments):\n";
    std::cout << "  PT3:  ++=" << h_pp_PT3->Integral()
              << "  --=" << h_mm_PT3->Integral()
              << "  CB="  << h_cb_PT3_total->Integral() << "\n";
    std::cout << "  PT2:  ++=" << h_pp_PT2->Integral()
              << "  --=" << h_mm_PT2->Integral()
              << "  CB="  << h_cb_PT2_total->Integral() << "\n\n";

    // ---------- Build per-segment EPEM histograms + signals ----------
    // First pass: gather raw EPEM PT3 yields per segment, to compute the
    // bin-magnitude fraction f_s = n_em_PT3_s / Σ n_em_PT3_s (use the
    // histogram integrals for self-consistency with the binned analysis).
    struct SegHists {
        TH1D* h_em_PT3 = nullptr;
        TH1D* h_em_PT2 = nullptr;
        TH1D* h_cb_PT3 = nullptr;
        TH1D* h_cb_PT2 = nullptr;
        TH1D* h_sig_PT3 = nullptr;
        TH1D* h_sig_PT2 = nullptr;
        TH1D* h_ratio   = nullptr;
        double frac     = 0.0;
        double n_em_PT3 = 0.0;
        double n_em_PT2 = 0.0;
        Pol0Result fit{};
    };
    std::vector<SegHists> segs(kNSeg);

    double sum_em_PT3 = 0.0;
    for (int s = 0; s < kNSeg; ++s) {
        const std::string cutPT3 =
            std::string(kCutPT3) + " && seg_idx==" + std::to_string(s);
        const std::string cutPT2 =
            std::string(kCutPT2) + " && seg_idx==" + std::to_string(s);
        segs[s].h_em_PT3 = drawWithCut(
            t_em, "m_ee", cutPT3.c_str(),
            "h_em_s" + std::to_string(s) + "_PT3");
        segs[s].h_em_PT2 = drawWithCut(
            t_em, "m_ee", cutPT2.c_str(),
            "h_em_s" + std::to_string(s) + "_PT2");
        segs[s].n_em_PT3 = segs[s].h_em_PT3->Integral();
        segs[s].n_em_PT2 = segs[s].h_em_PT2->Integral();
        sum_em_PT3 += segs[s].n_em_PT3;
    }

    std::cout << "Per-segment EPEM raw yields (histogram integrals):\n";
    for (int s = 0; s < kNSeg; ++s) {
        segs[s].frac = (sum_em_PT3 > 0.0)
            ? segs[s].n_em_PT3 / sum_em_PT3 : 0.0;
        printf("  seg %d : N_em_PT3 = %10.0f  N_em_PT2 = %10.0f  frac = %.4f\n",
               s, segs[s].n_em_PT3, segs[s].n_em_PT2, segs[s].frac);
    }
    std::cout << "  sum  = " << sum_em_PT3 << "\n\n";

    // Build per-segment CB (scaled total) + signal + ratio + pol0 fit.
    for (int s = 0; s < kNSeg; ++s) {
        const std::string sg = "s" + std::to_string(s);
        segs[s].h_cb_PT3 = scaleCB(h_cb_PT3_total, segs[s].frac,
                                   "h_cb_PT3_" + sg);
        segs[s].h_cb_PT2 = scaleCB(h_cb_PT2_total, segs[s].frac,
                                   "h_cb_PT2_" + sg);
        segs[s].h_sig_PT3 = makeSignal(segs[s].h_em_PT3, segs[s].h_cb_PT3,
                                       "h_sig_PT3_" + sg);
        segs[s].h_sig_PT2 = makeSignal(segs[s].h_em_PT2, segs[s].h_cb_PT2,
                                       "h_sig_PT2_" + sg);
        segs[s].h_ratio = makeRatio(segs[s].h_sig_PT2, segs[s].h_sig_PT3,
                                    kTrigCorr, "h_ratio_" + sg);
        segs[s].fit = fitPol0(segs[s].h_ratio, kFitLo, kFitHi, sg);
    }

    // ---------- Build FULL (all segments combined) signal histogram ----------
    // Reference: N_raw_signal_FULL in [kFitLo, kFitHi] is computed on the
    // entire EPEM dataset (no segment filter) using the same 58-bin layout
    // and the same CB subtraction as the per-segment analysis. The proper
    // V2-corrected reference yield is then
    //     N_corr_V2_full_proper = N_raw_signal_FULL * C_full
    // (C_full = 1.835 +- 0.038 from v2_flatness_subwindow_test extended fit).
    TH1D* h_em_PT3_full = drawWithCut(t_em, "m_ee", kCutPT3, "h_em_PT3_full");
    TH1D* h_sig_PT3_full =
        makeSignal(h_em_PT3_full, h_cb_PT3_total, "h_sig_PT3_full");

    // N_raw_signal_FULL integrated over the OLD range [0.15, 0.70] (for
    // documentation: should be ~310,886).
    double N_raw_full_070     = 0.0;
    double sN_raw_full_070    = 0.0;
    sumWithError(h_sig_PT3_full, 0.15, 0.70,
                 N_raw_full_070, sN_raw_full_070);

    // N_raw_signal_FULL integrated over the EXTENDED range [0.15, 0.80]
    // (this is the proper N to scale by C_full).
    double N_raw_full_080     = 0.0;
    double sN_raw_full_080    = 0.0;
    sumWithError(h_sig_PT3_full, kFitLo, kFitHi,
                 N_raw_full_080, sN_raw_full_080);

    // Proper full-dataset V2-corrected reference yield on [0.15, 0.80].
    // Error propagation: relative error of C_full dominates; we add the
    // raw-yield Poisson error in quadrature.
    const double N_corr_full_proper = N_raw_full_080 * kC_full;
    const double sN_corr_full_proper = std::sqrt(
        (sN_raw_full_080 * kC_full) * (sN_raw_full_080 * kC_full)
        + (N_raw_full_080 * kSC_full) * (N_raw_full_080 * kSC_full));

    std::cout << "FULL (all segments) signal yield (PT3, no seg filter):\n";
    printf("  N_raw_signal_FULL in [0.15, 0.70] = %10.3f +- %8.3f\n",
           N_raw_full_070, sN_raw_full_070);
    printf("  N_raw_signal_FULL in [0.15, 0.80] = %10.3f +- %8.3f\n",
           N_raw_full_080, sN_raw_full_080);
    printf("  additional signal in [0.70, 0.80] = %10.3f\n",
           N_raw_full_080 - N_raw_full_070);
    printf("  N_corr_V2_full_proper [0.15,0.80] = %10.3f +- %8.3f\n\n",
           N_corr_full_proper, sN_corr_full_proper);

    // ---------- PT3-count-weighted combination ----------
    // Use n_PT3 from calibration (long long) as weights.
    double num = 0.0;
    double w_sum = 0.0;
    for (int s = 0; s < kNSeg; ++s) {
        auto it = seg_n_pt3.find(s);
        if (it == seg_n_pt3.end()) {
            std::cerr << "[combine] WARNING: seg " << s
                      << " missing from calibration map\n";
            continue;
        }
        const double w = static_cast<double>(it->second);
        num    += w * segs[s].fit.C;
        w_sum  += w;
    }
    const double C_comb  = (w_sum > 0.0) ? num / w_sum : 0.0;
    double var_comb = 0.0;
    for (int s = 0; s < kNSeg; ++s) {
        auto it = seg_n_pt3.find(s);
        if (it == seg_n_pt3.end()) continue;
        const double w   = static_cast<double>(it->second);
        const double rho = (w_sum > 0.0) ? (w / w_sum) : 0.0;
        var_comb += rho * rho * segs[s].fit.sC * segs[s].fit.sC;
    }
    const double sC_comb = std::sqrt(var_comb);

    // ---------- Apply to integrated yield (segmented) ----------
    // N_corr_segmented = Σ_s [ Σ_bin h_sig_s_PT3(b) * C_s ]   on [0.15, 0.80]
    // Variance treats C_s as independent gaussians times the integrated bin
    // sum n_s; per-bin Poisson errors contribute via σ²_n_s * C_s².
    double N_corr_seg  = 0.0;
    double var_N_seg   = 0.0;
    for (int s = 0; s < kNSeg; ++s) {
        double n_s = 0.0;
        double sn_s = 0.0;
        sumWithError(segs[s].h_sig_PT3, kFitLo, kFitHi, n_s, sn_s);
        const double C_s   = segs[s].fit.C;
        const double sC_s  = segs[s].fit.sC;
        N_corr_seg  += n_s * C_s;
        // var = (sn * C)^2 + (n * sC)^2  (treating sn and sC independent)
        var_N_seg   += (sn_s * C_s) * (sn_s * C_s)
                     + (n_s  * sC_s) * (n_s  * sC_s);
    }
    const double sN_corr_seg = std::sqrt(var_N_seg);

    // ---------- Print summary table ----------
    std::cout << "\n=================================================="
              << "==================================================\n";
    std::cout << "PER-SEGMENT V2 FIT SUMMARY (pol0 on [" << kFitLo
              << ", " << kFitHi << "])\n";
    std::cout << "=================================================="
              << "==================================================\n";
    printf("| %-22s | %12s | %8s | %8s | %10s |\n",
           "Segment", "n_PT3", "C_s", "sigma_C", "chi2/ndf");
    printf("|------------------------|--------------|----------|"
           "----------|------------|\n");
    for (int s = 0; s < kNSeg; ++s) {
        auto it = seg_n_pt3.find(s);
        const long long n = (it != seg_n_pt3.end()) ? it->second : 0;
        printf("| seg %-18d | %12lld | %8.4f | %8.4f | %5.2f/%-4d |\n",
               s, n, segs[s].fit.C, segs[s].fit.sC,
               segs[s].fit.chi2, segs[s].fit.ndf);
    }
    printf("|------------------------|--------------|----------|"
           "----------|------------|\n");
    // Per-segment n_PT3 counted in the fit range (EPEM PT3 raw histogram).
    {
        auto countInRangeLocal = [](TH1D* h, double xlo, double xhi) -> double {
            double sum = 0.0;
            for (int b = 1; b <= h->GetNbinsX(); ++b) {
                const double lo = h->GetXaxis()->GetBinLowEdge(b);
                const double hi = h->GetXaxis()->GetBinUpEdge(b);
                if (lo >= xlo - 1e-9 && hi <= xhi + 1e-9) {
                    sum += h->GetBinContent(b);
                }
            }
            return sum;
        };
        printf("n_PT3 in [%.2f, %.2f] per segment (EPEM raw):\n",
               kFitLo, kFitHi);
        for (int s = 0; s < kNSeg; ++s) {
            const double n_in =
                countInRangeLocal(segs[s].h_em_PT3, kFitLo, kFitHi);
            printf("  seg %d : n_PT3_in_range = %.0f\n", s, n_in);
        }
    }
    printf("| %-22s | %12lld | %8.4f | %8.4f | %10s |\n",
           "combined (PT3-weighted)", sum_n_pt3, C_comb, sC_comb, "-");
    printf("| %-22s | %12s | %8.4f | %8.4f | %10s |\n",
           "full V2 ext [0.15,0.80]", "-", kC_full, kSC_full, "-");
    std::cout << "\n";

    // Comparison vs full
    {
        const double d  = C_comb - kC_full;
        const double sc = std::sqrt(sC_comb * sC_comb + kSC_full * kSC_full);
        const double nsig = (sc > 0.0) ? std::abs(d) / sc : 0.0;
        printf("C_combined - C_full = %+.4f +- %.4f  (n_sigma = %.2f)\n",
               d, sc, nsig);
    }

    // ---------- Consistency tests across segments ----------
    // (A) chi² vs constant hypothesis: use the PT3-count-weighted combined
    //     value as the constant; chi² = Σ_s (C_s - C_comb)² / σ_C_s².
    //     n_segments = kNSeg, ndf = kNSeg - 1.
    {
        double chi2_const = 0.0;
        for (int s = 0; s < kNSeg; ++s) {
            const double d  = segs[s].fit.C - C_comb;
            const double sc = segs[s].fit.sC;
            if (sc > 0.0) chi2_const += (d * d) / (sc * sc);
        }
        const int ndf_const = kNSeg - 1;
        printf("\n[consistency] chi2 vs constant hypothesis "
               "(C = C_combined = %.4f):\n", C_comb);
        printf("  chi2/ndf = %.3f/%d = %.3f\n",
               chi2_const, ndf_const,
               (ndf_const > 0) ? chi2_const / ndf_const : 0.0);
    }

    // (B) Largest pairwise n_sigma between any two segments.
    {
        int    best_i = -1, best_j = -1;
        double best_n = 0.0;
        for (int i = 0; i < kNSeg; ++i) {
            for (int j = i + 1; j < kNSeg; ++j) {
                const double d  = segs[i].fit.C - segs[j].fit.C;
                const double sc = std::sqrt(
                    segs[i].fit.sC * segs[i].fit.sC +
                    segs[j].fit.sC * segs[j].fit.sC);
                const double nsig = (sc > 0.0) ? std::abs(d) / sc : 0.0;
                if (nsig > best_n) {
                    best_n = nsig;
                    best_i = i;
                    best_j = j;
                }
            }
        }
        if (best_i >= 0) {
            printf("[consistency] largest pairwise n_sigma: "
                   "seg %d vs seg %d -> n_sigma = %.2f\n",
                   best_i, best_j, best_n);
        }
    }

    // (C) Per-segment n_sigma vs C_full reference.
    {
        printf("[consistency] per-segment n_sigma vs C_full = %.4f +- %.4f:\n",
               kC_full, kSC_full);
        for (int s = 0; s < kNSeg; ++s) {
            const double d  = segs[s].fit.C - kC_full;
            const double sc = std::sqrt(
                segs[s].fit.sC * segs[s].fit.sC +
                kSC_full * kSC_full);
            const double nsig = (sc > 0.0) ? std::abs(d) / sc : 0.0;
            printf("  seg %d : C_s - C_full = %+.4f +- %.4f  (n_sigma = %.2f)\n",
                   s, d, sc, nsig);
        }
    }

    // Yield comparison vs full (PROPER reference)
    // The proper reference is N_raw_signal_FULL on [0.15, 0.80] times C_full,
    // NOT the old hard-coded 570,323 (which was N_raw [0.15,0.70] * C_full).
    {
        const double d  = N_corr_seg - N_corr_full_proper;
        const double sc = std::sqrt(sN_corr_seg * sN_corr_seg
                                    + sN_corr_full_proper
                                      * sN_corr_full_proper);
        const double nsig = (sc > 0.0) ? std::abs(d) / sc : 0.0;
        printf("\n--- YIELD COMPARISON [0.15, 0.80]  (PROPER reference) ---\n");
        printf("  N_raw_signal_FULL  [0.15, 0.70] = %10.3f  (doc only)\n",
               N_raw_full_070);
        printf("  N_raw_signal_FULL  [0.15, 0.80] = %10.3f  (proper)\n",
               N_raw_full_080);
        printf("  C_full (extended)               = %.4f +- %.4f\n",
               kC_full, kSC_full);
        printf("  N_corr_segmented  [0.15, 0.80]  = %10.3f +- %8.3f\n",
               N_corr_seg, sN_corr_seg);
        printf("  N_corr_V2_full_proper  [0.15,0.80] = %10.3f +- %8.3f\n",
               N_corr_full_proper, sN_corr_full_proper);
        printf("  Delta (seg - full_proper)       = %+10.3f +- %8.3f"
               "  (n_sigma = %.2f)\n",
               d, sc, nsig);
    }

    // ---------- Plot styling ----------
    // Color cycle for the 5 segments.
    const Color_t kSegCols[kNSeg] = {
        kBlack, kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1
    };
    const Style_t kSegMrks[kNSeg] = {20, 21, 22, 23, 33};

    const double ymin = 0.5, ymax = 4.0;

    // Style ratio histograms.
    for (int s = 0; s < kNSeg; ++s) {
        styleHist(segs[s].h_ratio, kSegCols[s], kSegMrks[s]);
        segs[s].h_ratio->SetMarkerSize(1.0);
    }

    // Count n_PT3 entries in the fit range per segment (for legend on
    // individual plots). Use the EPEM PT3 raw histogram.
    auto countInRange = [](TH1D* h, double xlo, double xhi) -> double {
        double sum = 0.0;
        for (int b = 1; b <= h->GetNbinsX(); ++b) {
            const double lo = h->GetXaxis()->GetBinLowEdge(b);
            const double hi = h->GetXaxis()->GetBinUpEdge(b);
            if (lo >= xlo - 1e-9 && hi <= xhi + 1e-9) {
                sum += h->GetBinContent(b);
            }
        }
        return sum;
    };

    // ---------- Plots 1-5: individual per-segment plots ----------
    for (int s = 0; s < kNSeg; ++s) {
        const std::string cname = "c_v2_seg_" + std::to_string(s);
        TCanvas* cs = new TCanvas(cname.c_str(),
                                  ("Per-segment V2 seg " +
                                   std::to_string(s)).c_str(),
                                  1100, 800);
        cs->SetMargin(0.12, 0.04, 0.12, 0.10);

        const std::string title =
            "63#upointPT2/PT3 signal ratio #minus segment "
            + std::to_string(s)
            + ";M_{e^{+}e^{-}} [GeV/c^{2}];"
            + "63 #upoint N_{PT2}^{sig} / N_{PT3}^{sig}";

        TH1D* h_axis_s = static_cast<TH1D*>(
            segs[s].h_ratio->Clone(("h_axis_s" + std::to_string(s)).c_str()));
        h_axis_s->SetDirectory(nullptr);
        h_axis_s->Reset();
        h_axis_s->SetTitle(title.c_str());
        h_axis_s->GetXaxis()->SetRangeUser(kXmin, kXmax);
        h_axis_s->GetYaxis()->SetRangeUser(ymin, ymax);
        h_axis_s->Draw("AXIS");

        // Reference y=1 dashed gray line.
        TLine* l_one = new TLine(kXmin, 1.0, kXmax, 1.0);
        l_one->SetLineColor(kGray + 1);
        l_one->SetLineStyle(2);
        l_one->SetLineWidth(1);
        l_one->Draw();

        // Vertical dotted lines at fit range edges.
        TLine* lv_lo = new TLine(kFitLo, ymin, kFitLo, ymax);
        TLine* lv_hi = new TLine(kFitHi, ymin, kFitHi, ymax);
        for (TLine* lv : {lv_lo, lv_hi}) {
            lv->SetLineColor(kGray + 1);
            lv->SetLineStyle(3);
            lv->SetLineWidth(1);
            lv->Draw();
        }

        // Pol0 fit horizontal line over [kFitLo, kFitHi].
        TLine* l_fit = new TLine(kFitLo, segs[s].fit.C,
                                 kFitHi, segs[s].fit.C);
        l_fit->SetLineColor(kSegCols[s]);
        l_fit->SetLineStyle(1);
        l_fit->SetLineWidth(2);
        l_fit->Draw();

        // Ratio markers + error bars on top.
        segs[s].h_ratio->Draw("E1 SAME");

        // Legend top-right.
        TLegend* leg_s = new TLegend(0.55, 0.70, 0.96, 0.90);
        leg_s->SetBorderSize(0);
        leg_s->SetFillStyle(0);
        leg_s->SetTextSize(0.030);
        leg_s->AddEntry(segs[s].h_ratio,
            Form("segment %d", s), "pl");
        leg_s->AddEntry(l_fit,
            Form("C_{%d} = %.3f #pm %.3f",
                 s, segs[s].fit.C, segs[s].fit.sC), "l");
        leg_s->AddEntry((TObject*)nullptr,
            Form("#chi^{2}/ndf = %.1f/%d",
                 segs[s].fit.chi2, segs[s].fit.ndf), "");
        const double n_in_range = countInRange(
            segs[s].h_em_PT3, kFitLo, kFitHi);
        // Comma-format the integer count.
        char nbuf[64];
        {
            long long nll = static_cast<long long>(n_in_range);
            // Reverse comma insertion.
            char tmp[64];
            int  ti = 0;
            if (nll == 0) { tmp[ti++] = '0'; }
            while (nll > 0) {
                if (ti > 0 && (ti % 4) == 3) { tmp[ti++] = ','; }
                tmp[ti++] = static_cast<char>('0' + (nll % 10));
                nll /= 10;
            }
            int ni = 0;
            for (int j = ti - 1; j >= 0; --j) nbuf[ni++] = tmp[j];
            nbuf[ni] = '\0';
        }
        leg_s->AddEntry((TObject*)nullptr,
            Form("n_{PT3} (in [%.2f, %.2f]) = %s",
                 kFitLo, kFitHi, nbuf), "");
        leg_s->Draw();

        const std::string base_s = "plots/output/v2_per_segment_s"
                                 + std::to_string(s);
        cs->SaveAs((base_s + ".pdf").c_str());
        cs->SaveAs((base_s + ".png").c_str());
        std::cout << "Saved: " << base_s << ".{pdf,png}\n";
    }

    // ---------- Plot 6: overlay of all 5 segments ----------
    TCanvas* c_ov = new TCanvas("c_v2_per_segment_overlay",
                                "Per-segment V2 overlay",
                                1400, 900);
    c_ov->SetMargin(0.10, 0.04, 0.12, 0.10);

    TH1D* h_axis_ov = static_cast<TH1D*>(segs[0].h_ratio->Clone("h_axis_ov"));
    h_axis_ov->SetDirectory(nullptr);
    h_axis_ov->Reset();
    h_axis_ov->SetTitle(
        "63#upointPT2/PT3 signal ratio #minus all 5 segments overlay"
        ";M_{e^{+}e^{-}} [GeV/c^{2}];"
        "63 #upoint N_{PT2}^{sig} / N_{PT3}^{sig}");
    h_axis_ov->GetXaxis()->SetRangeUser(kXmin, kXmax);
    h_axis_ov->GetYaxis()->SetRangeUser(ymin, ymax);
    h_axis_ov->Draw("AXIS");

    // Reference y=1 dashed gray line.
    TLine* l_one_ov = new TLine(kXmin, 1.0, kXmax, 1.0);
    l_one_ov->SetLineColor(kGray + 1);
    l_one_ov->SetLineStyle(2);
    l_one_ov->SetLineWidth(1);
    l_one_ov->Draw();

    // For each segment: draw a connecting LINE (HIST L) - no errors/markers.
    std::vector<TH1D*> h_lines_ov(kNSeg, nullptr);
    for (int s = 0; s < kNSeg; ++s) {
        TH1D* hl = static_cast<TH1D*>(segs[s].h_ratio->Clone(
            ("h_ov_line_s" + std::to_string(s)).c_str()));
        hl->SetDirectory(nullptr);
        hl->SetLineColor(kSegCols[s]);
        hl->SetLineWidth(2);
        hl->SetLineStyle(1);
        hl->SetMarkerStyle(1);   // no visible marker
        hl->SetMarkerSize(0);
        // Clear errors so HIST L doesn't include them visually.
        for (int b = 1; b <= hl->GetNbinsX(); ++b) hl->SetBinError(b, 0.0);
        hl->Draw("HIST L SAME");
        h_lines_ov[s] = hl;
    }

    // Per-segment pol0 horizontal lines over [kFitLo, kFitHi].
    std::vector<TLine*> hlines_ov;
    for (int s = 0; s < kNSeg; ++s) {
        TLine* ll = new TLine(kFitLo, segs[s].fit.C,
                              kFitHi, segs[s].fit.C);
        ll->SetLineColor(kSegCols[s]);
        ll->SetLineWidth(2);
        ll->SetLineStyle(1);
        ll->Draw();
        hlines_ov.push_back(ll);
    }

    // Legend top-right.
    TLegend* leg_ov = new TLegend(0.55, 0.62, 0.96, 0.90);
    leg_ov->SetBorderSize(0);
    leg_ov->SetFillStyle(0);
    leg_ov->SetTextSize(0.028);
    leg_ov->SetHeader(
        Form("V2 pol0 fit on [%.2f, %.2f]", kFitLo, kFitHi), "C");
    for (int s = 0; s < kNSeg; ++s) {
        leg_ov->AddEntry(h_lines_ov[s],
            Form("seg %d: C_{%d} = %.3f #pm %.3f",
                 s, s, segs[s].fit.C, segs[s].fit.sC), "l");
    }
    leg_ov->Draw();

    const std::string base_ov = "plots/output/v2_per_segment_overlay";
    c_ov->SaveAs((base_ov + ".pdf").c_str());
    c_ov->SaveAs((base_ov + ".png").c_str());
    std::cout << "Saved: " << base_ov << ".{pdf,png}\n";

    f_em->Close();
    f_pp->Close();
    f_mm->Close();
}
