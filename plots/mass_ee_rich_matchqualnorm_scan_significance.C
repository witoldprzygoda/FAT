// mass_ee_rich_matchqualnorm_scan_significance.C — sister of
// mass_ee_rich_padnum_scan_significance.C, but scanning thresholds on
// `richmatchqualitynorm` (track ↔ ring match-quality, NORMALIZED) instead
// of `rich_padnum`. Smaller `richmatchqualitynorm` = BETTER match, so the
// cut convention here is `qual <= cut_value` (lower cut = tighter; higher
// cut = looser). With cut = cut_max all events pass.
//
// What is scanned:
//   • EpEm pair (signal candidates): asymmetric cut grid
//       ep_richmatchqualitynorm <= ep_cut   AND
//       em_richmatchqualitynorm <= em_cut
//   • EpEp pair (like-sign + + for CB): SYNCHRONIZED cut on ep_cut
//       (both legs must satisfy qual <= ep_cut — both legs are e+)
//   • EmEm pair (like-sign − − for CB): SYNCHRONIZED cut on em_cut
//       (both legs must satisfy qual <= em_cut — both legs are e−)
//   The synchronization rationale is identical to the padnum macro.
//
// Default cut grid: cut ∈ [0, 22] (23 values per leg → 529 combinations).
// Convention: qual <= cut_value passes; cut=0 is tightest, cut=22 is loosest.
//
// What is computed at each grid point (ep_cut, em_cut), separately for
// PT2 and PT3:
//   N_epem(m_ee)  = entries from EpEm passing both cuts
//   N_epep(m_ee)  = entries from EpEp passing ep_cut on both legs
//   N_emem(m_ee)  = entries from EmEm passing em_cut on both legs
//   CB(m_ee)      = 2 · √(N_epep · N_emem)
//   sig(m_ee)     = N_epem − CB
//   significance per window per trigger = ΣS / √(ΣS + 2·ΣCB)
//
// Per-window selection: a SINGLE (ep_cut, em_cut) is chosen for both PT2 and
// PT3 (shared cuts). Objective = sig_PT2 + sig_PT3 (sum of per-trigger window
// significances at the same cuts). With PT3 ~5× more significant than PT2 in
// most windows, the chosen cuts are typically very close to the PT3-only
// optimum, but PT2 still contributes to the figure of merit.
//
// Four m_ee windows, one best (ep_cut, em_cut) per window (4 best cases):
//   A) π⁰         0.000 - 0.135 GeV/c²
//   B) η Dalitz   0.135 - 0.600 GeV/c²
//   C) high mass  0.600 - 1.400 GeV/c²
//   D) full       0.000 - 1.400 GeV/c²
//
// TOP determination: a separate "above-π⁰" significance is computed on
// [0.14, 1.4] GeV at the chosen cuts of each window (per trigger), then
// summed: top_sig_sum = top_sig_PT2 + top_sig_PT3. The window with the
// largest top_sig_sum gets the ★ TOP mark.
//
// Outputs (5 files total):
//   1.   plots/output/scan_rich_matchqualnorm_heatmap_significance.{pdf,png}
//        — 2×4 grid (rows: PT3/PT2; cols: A/B/C/D) of per-trigger
//          significance vs (ep_cut, em_cut). Best grid-point (shared
//          between PT2 and PT3 rows for a given column) marked with a ★.
//
//   2-5. plots/output/scan_rich_matchqualnorm_best_<window>.{pdf,png}
//        — One 2×2 four-panel canvas per window:
//          Top-left:    m_ee at PT3 (all=black, CB=red, sig=blue), window
//                       highlighted in green.
//          Top-right:   m_ee at PT2 (same colour scheme), window highlighted.
//          Bottom-left: 63·N_PT2 / N_PT3 (trigger correction factor), const
//                       fit on signal in [0.1, 1.4], Y axis fixed [0, 5].
//          Bottom-right: N_PT3 / (63·N_PT2) (efficiency), const fit in
//                        [0.1, 1.4], Y axis fixed [0, 1].
//
// Caveats:
//   • PT2 statistics are ~27× smaller than PT3 (downscale 64 in raw data).
//     PT2-based significances will be lower in absolute value; their best
//     cuts will tend to be LOOSER.
//   • Significance formula assumes Poisson on all three samples and uses
//     2·CB as the effective background variance (the 2× factor is the
//     standard like-sign uncertainty inflation).
//   • Signal can go negative in low-stats bins (CB > N_epem); the window
//     integral may still come out positive even if individual bins don't.
//
// Usage:
//   root -l -b -q plots/mass_ee_rich_matchqualnorm_scan_significance.C
//   root -l -b -q 'plots/mass_ee_rich_matchqualnorm_scan_significance.C(0, 22)'
//
// @author Witold Przygoda (witold.przygoda@uj.edu.pl)
// @date 2026

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TMarker.h>
#include <TBox.h>
#include <TStyle.h>
#include <TSystem.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

// -------- m_ee binning (5 MeV bins — windows align exactly with bin edges) -
constexpr int    kNb   = 280;
constexpr double kXmin = 0.0;
constexpr double kXmax = 1.4;
constexpr double kBinW = (kXmax - kXmin) / kNb;        // = 0.005 GeV/c²

// -------- TOP-criterion range (above π⁰): used to pick the "TOP" case ------
// (in earlier iteration this was the full 0 – 1.4 range, dominated by π⁰).
constexpr double kTopLo = 0.140;
constexpr double kTopHi = 1.400;

// -------- Fit range for constant fits on the trigger ratios ----------------
constexpr double kFitLo = 0.100;
constexpr double kFitHi = 1.400;

// -------- Thousands separator (apostrophe) — for stdout summary table ------
inline std::string fmtApos(double v) {
    const long long iabs = static_cast<long long>(std::llround(std::abs(v)));
    std::string s = std::to_string(iabs);
    for (int pos = static_cast<int>(s.size()) - 3; pos > 0; pos -= 3)
        s.insert(pos, "'");
    return (v < 0 ? "-" + s : s);
}

// -------- Trigger encoding ---------------------------------------------
enum Trigger : int { TRG_PT2 = 0, TRG_PT3 = 1, NUM_TRG = 2 };
const std::array<std::string, NUM_TRG> kTrgName       = {"PT2", "PT3"};
const std::array<std::string, NUM_TRG> kTrgLabel      = {"PT2 trigger",
                                                         "PT3 trigger"};
const std::array<int,         NUM_TRG> kTrgBit        = {4096, 8192};

// -------- Mass windows -------------------------------------------------
struct Window {
    std::string short_id;
    std::string label;
    double      lo;
    double      hi;
};
const std::vector<Window> kWindows = {
    {"A_pi0",      "A: #pi^{0} (0 - 0.135)",         0.000, 0.135},
    {"B_etaDal",   "B: #eta Dal. (0.135 - 0.6)",     0.135, 0.600},
    {"C_highmass", "C: high mass (0.6 - 1.4)",       0.600, 1.400},
    {"D_full",     "D: full (0 - 1.4)",              0.000, 1.400},
};

// -------- 1D scan storage shapes ---------------------------------------
// Indexed by (cut - cut_min) on [0, n_cuts). For EpEm there's a second
// dimension (em side). For EpEp only ep_cut matters; for EmEm only em_cut.
//
// h_epem[ep_idx][em_idx][trg]  — full 2D × trigger
// h_epep[ep_idx]      [trg]    — like-sign + + filtered with ep_cut on BOTH legs
// h_emem[em_idx]      [trg]    — like-sign − − filtered with em_cut on BOTH legs

using HistMat3 = std::vector<std::vector<std::array<TH1D*, NUM_TRG>>>;
using HistMat2 = std::vector<std::array<TH1D*, NUM_TRG>>;

// -------- Significance result per window -----------------------------
// One (ep_cut, em_cut) is shared by both PT2 and PT3 (chosen to maximize the
// SUM of the per-trigger window significances). Stores per-trigger window
// stats AND per-trigger top-range stats (range 0.14-1.4) AT those common cuts.
struct BestPoint {
    int    ep_cut       = 0;
    int    em_cut       = 0;
    std::array<double, NUM_TRG> S          = {0.0, 0.0};   // window S per trigger
    std::array<double, NUM_TRG> B          = {0.0, 0.0};   // window CB per trigger
    std::array<double, NUM_TRG> signif     = {0.0, 0.0};   // window S/sqrt(S+2B)
    std::array<double, NUM_TRG> top_S      = {0.0, 0.0};   // top-range S per trigger
    std::array<double, NUM_TRG> top_B      = {0.0, 0.0};   // top-range CB per trigger
    std::array<double, NUM_TRG> top_signif = {0.0, 0.0};   // top-range signif
    double signif_sum     = 0.0;   // sig_PT2 + sig_PT3  — objective being maximized
    double top_signif_sum = 0.0;   // top_sig_PT2 + top_sig_PT3 — used to pick TOP
    bool   is_top         = false; // window with max top_signif_sum
};

// -------- Trigger flag extraction --------------------------------------
inline int trgFromBit(int trigbit) {
    if (trigbit == kTrgBit[TRG_PT3]) return TRG_PT3;
    if (trigbit == kTrgBit[TRG_PT2]) return TRG_PT2;
    return -1;
}

// -------- Window integral helper -------------------------------------
// Sum bins fully inside [lo, hi] of m_ee. Uses bin-edge alignment because
// our 5 MeV grid contains 0.135, 0.6, 1.4 as exact edges.
inline std::pair<double,double> windowIntegralAndErr(TH1D* h,
                                                     double lo, double hi) {
    if (!h) return {0.0, 0.0};
    const int b_lo = h->GetXaxis()->FindFixBin(lo + 0.5 * kBinW);
    const int b_hi = h->GetXaxis()->FindFixBin(hi - 0.5 * kBinW);
    if (b_hi < b_lo) return {0.0, 0.0};
    double err = 0.0;
    const double val = h->IntegralAndError(b_lo, b_hi, err);
    return {val, err};
}

// -------- CB and signal builders (in fine bins) ----------------------
TH1D* makeCB(TH1D* pp, TH1D* mm, const std::string& name) {
    TH1D* cb = static_cast<TH1D*>(pp->Clone(name.c_str()));
    cb->SetDirectory(nullptr);
    cb->Reset();
    for (int b = 1; b <= cb->GetNbinsX(); ++b) {
        const double np = pp->GetBinContent(b);
        const double ep = pp->GetBinError(b);
        const double nm = mm->GetBinContent(b);
        const double em = mm->GetBinError(b);
        if (np > 0 && nm > 0) {
            const double cval = 2.0 * std::sqrt(np * nm);
            const double rel  = std::sqrt(std::pow(ep / np, 2) +
                                           std::pow(em / nm, 2));
            cb->SetBinContent(b, cval);
            cb->SetBinError  (b, cval * 0.5 * rel);
        }
    }
    return cb;
}

TH1D* makeSignal(TH1D* all, TH1D* cb, const std::string& name) {
    TH1D* sig = static_cast<TH1D*>(all->Clone(name.c_str()));
    sig->SetDirectory(nullptr);
    sig->Add(cb, -1.0);
    return sig;
}

// Per-bin ratio with full error propagation (used for trigger ratios).
TH1D* makeRatio(TH1D* num, TH1D* den, const std::string& name, double scale) {
    TH1D* r = static_cast<TH1D*>(num->Clone(name.c_str()));
    r->SetDirectory(nullptr);
    r->Reset();
    for (int b = 1; b <= r->GetNbinsX(); ++b) {
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

// -------- Rebin a fine 5 MeV hist down to 20 MeV for plotting -------
TH1D* rebinFor20(TH1D* h, const std::string& name) {
    TH1D* r = static_cast<TH1D*>(h->Rebin(4, name.c_str()));
    r->SetDirectory(nullptr);
    return r;
}

void styleAll(TH1D* h, Color_t c, Style_t m) {
    h->SetLineColor(c); h->SetMarkerColor(c); h->SetMarkerStyle(m);
    h->SetMarkerSize(0.8); h->SetLineWidth(1);
    h->GetXaxis()->SetTitleSize(0.050);
    h->GetYaxis()->SetTitleSize(0.050);
    h->GetXaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetLabelSize(0.045);
    h->GetXaxis()->SetTitleOffset(1.05);
    h->GetYaxis()->SetTitleOffset(1.20);
}

struct FitResult { double a; double e; double chi2_ndf; int ndf; };

FitResult fitConst(TH1D* h, double xlo, double xhi, Color_t color) {
    const std::string fname = std::string(h->GetName()) + "_fc";
    TF1* f = new TF1(fname.c_str(), "[0]", xlo, xhi);
    f->SetLineColor(color); f->SetLineWidth(2);
    h->Fit(f, "RQ");
    FitResult r;
    r.a        = f->GetParameter(0);
    r.e        = f->GetParError(0);
    r.ndf      = f->GetNDF();
    r.chi2_ndf = (r.ndf > 0) ? f->GetChisquare() / r.ndf : 0.0;
    return r;
}

}  // anonymous namespace

void mass_ee_rich_matchqualnorm_scan_significance(int cut_min = 0,
                                                  int cut_max = 22) {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetErrorX(0.0);
    gStyle->SetPalette(kBird);

    const int n_cuts = cut_max - cut_min + 1;
    std::cout << "Scanning ep_richmatchqualitynorm, em_richmatchqualitynorm on ["
              << cut_min << ", " << cut_max << "]  (n_cuts=" << n_cuts
              << ", " << n_cuts * n_cuts << " combinations)\n"
              << "Convention: qual <= cut_value passes  "
              << "(cut=0 tightest, cut=" << cut_max << " loosest).\n";

    // ============================================================================
    // Allocate storage
    // ============================================================================
    HistMat3 h_epem(n_cuts,
        std::vector<std::array<TH1D*, NUM_TRG>>(n_cuts));
    HistMat2 h_epep(n_cuts);
    HistMat2 h_emem(n_cuts);

    for (int i = 0; i < n_cuts; ++i) {
        for (int t = 0; t < NUM_TRG; ++t) {
            h_epep[i][t] = new TH1D(Form("h_epep_e%d_t%d", i + cut_min, t),
                                    "", kNb, kXmin, kXmax);
            h_epep[i][t]->Sumw2();
            h_epep[i][t]->SetDirectory(nullptr);
            h_emem[i][t] = new TH1D(Form("h_emem_m%d_t%d", i + cut_min, t),
                                    "", kNb, kXmin, kXmax);
            h_emem[i][t]->Sumw2();
            h_emem[i][t]->SetDirectory(nullptr);
        }
        for (int j = 0; j < n_cuts; ++j) {
            for (int t = 0; t < NUM_TRG; ++t) {
                h_epem[i][j][t] = new TH1D(
                    Form("h_epem_e%d_m%d_t%d", i + cut_min, j + cut_min, t),
                    "", kNb, kXmin, kXmax);
                h_epem[i][j][t]->Sumw2();
                h_epem[i][j][t]->SetDirectory(nullptr);
            }
        }
    }

    // ============================================================================
    // Pass 1 — EpEm: fill all 16×16 (ep_cut, em_cut) combinations the event
    //                passes (both legs ≥ cut).
    // ============================================================================
    std::cout << "\nPass 1/3 — EpEm (output_epem_exp.root)...\n";
    Long64_t n_em_total = 0, n_em_use = 0;
    {
        TFile f("output_epem_exp.root", "READ");
        TTreeReader r("dilepton_nt", &f);
        TTreeReaderValue<float> mee(r, "m_ee");
        TTreeReaderValue<float> trigb(r, "trigbit");
        TTreeReaderValue<float> ep_q(r, "ep_richmatchqualitynorm");
        TTreeReaderValue<float> em_q(r, "em_richmatchqualitynorm");

        while (r.Next()) {
            ++n_em_total;
            const int trg = trgFromBit(static_cast<int>(*trigb));
            if (trg < 0) continue;
            // Convention: qual <= cut passes. Lowest passing cut for an event
            // with quality Q is ceil(Q); highest passing cut is cut_max.
            const int min_e = std::max(cut_min,
                static_cast<int>(std::ceil(std::max(0.0f, *ep_q))));
            const int min_m = std::max(cut_min,
                static_cast<int>(std::ceil(std::max(0.0f, *em_q))));
            if (min_e > cut_max || min_m > cut_max) continue;
            ++n_em_use;
            const double m = *mee;
            for (int e = min_e; e <= cut_max; ++e) {
                auto& row = h_epem[e - cut_min];
                for (int mc = min_m; mc <= cut_max; ++mc) {
                    row[mc - cut_min][trg]->Fill(m);
                }
            }
        }
    }
    std::cout << "  scanned " << n_em_total << ", kept " << n_em_use
              << " (" << (100.0 * n_em_use / std::max<Long64_t>(1, n_em_total))
              << " %)\n";

    // ============================================================================
    // Pass 2 — EpEp: both legs share ep_cut (synchronized).
    // ============================================================================
    std::cout << "Pass 2/3 — EpEp (output_epep_exp.root)...\n";
    Long64_t n_pp_total = 0, n_pp_use = 0;
    {
        TFile f("output_epep_exp.root", "READ");
        TTreeReader r("dilepton_nt", &f);
        TTreeReaderValue<float> mee(r, "m_ee");
        TTreeReaderValue<float> trigb(r, "trigbit");
        TTreeReaderValue<float> ep_q(r, "ep_richmatchqualitynorm");
        TTreeReaderValue<float> em_q(r, "em_richmatchqualitynorm");

        while (r.Next()) {
            ++n_pp_total;
            const int trg = trgFromBit(static_cast<int>(*trigb));
            if (trg < 0) continue;
            // Both legs are e+; both must satisfy qual <= ep_cut.
            // Lowest passing cut = ceil(max(Q1, Q2)).
            const int min_c = std::max(cut_min,
                static_cast<int>(std::ceil(
                    std::max(0.0f, std::max(*ep_q, *em_q)))));
            if (min_c > cut_max) continue;
            ++n_pp_use;
            const double m = *mee;
            for (int e = min_c; e <= cut_max; ++e) {
                h_epep[e - cut_min][trg]->Fill(m);
            }
        }
    }
    std::cout << "  scanned " << n_pp_total << ", kept " << n_pp_use
              << " (" << (100.0 * n_pp_use / std::max<Long64_t>(1, n_pp_total))
              << " %)\n";

    // ============================================================================
    // Pass 3 — EmEm: both legs share em_cut (synchronized).
    // ============================================================================
    std::cout << "Pass 3/3 — EmEm (output_emem_exp.root)...\n";
    Long64_t n_mm_total = 0, n_mm_use = 0;
    {
        TFile f("output_emem_exp.root", "READ");
        TTreeReader r("dilepton_nt", &f);
        TTreeReaderValue<float> mee(r, "m_ee");
        TTreeReaderValue<float> trigb(r, "trigbit");
        TTreeReaderValue<float> ep_q(r, "ep_richmatchqualitynorm");
        TTreeReaderValue<float> em_q(r, "em_richmatchqualitynorm");

        while (r.Next()) {
            ++n_mm_total;
            const int trg = trgFromBit(static_cast<int>(*trigb));
            if (trg < 0) continue;
            // Both legs are e−; both must satisfy qual <= em_cut.
            const int min_c = std::max(cut_min,
                static_cast<int>(std::ceil(
                    std::max(0.0f, std::max(*ep_q, *em_q)))));
            if (min_c > cut_max) continue;
            ++n_mm_use;
            const double m = *mee;
            for (int e = min_c; e <= cut_max; ++e) {
                h_emem[e - cut_min][trg]->Fill(m);
            }
        }
    }
    std::cout << "  scanned " << n_mm_total << ", kept " << n_mm_use
              << " (" << (100.0 * n_mm_use / std::max<Long64_t>(1, n_mm_total))
              << " %)\n";

    // ============================================================================
    // Scan: at each (ep_cut, em_cut, trigger) compute (S, B, signif) per window
    // ============================================================================
    std::cout << "\nComputing significance landscape...\n";

    // 2D maps of significance per window × trigger (rows: em_cut, cols: ep_cut).
    std::vector<std::vector<TH2D*>> h_signif(
        kWindows.size(), std::vector<TH2D*>(NUM_TRG, nullptr));
    for (size_t w = 0; w < kWindows.size(); ++w) {
        for (int t = 0; t < NUM_TRG; ++t) {
            auto* h = new TH2D(
                Form("h_signif_%s_%s",
                     kWindows[w].short_id.c_str(), kTrgName[t].c_str()),
                Form("Significance — %s / %s;ep_richmatchqualitynorm #leq;em_richmatchqualitynorm #leq;S / #sqrt{S + 2B}",
                     kWindows[w].label.c_str(), kTrgLabel[t].c_str()),
                n_cuts, cut_min - 0.5, cut_max + 0.5,
                n_cuts, cut_min - 0.5, cut_max + 0.5);
            h->SetDirectory(nullptr);
            h_signif[w][t] = h;
        }
    }

    // Single best per window. Objective: sig_PT2 + sig_PT3.
    std::vector<BestPoint> best(kWindows.size());

    for (int ei = 0; ei < n_cuts; ++ei) {
        for (int mi = 0; mi < n_cuts; ++mi) {
            // Pre-compute CB per trigger (reused across all 4 windows).
            std::array<TH1D*, NUM_TRG> cb_t;
            for (int t = 0; t < NUM_TRG; ++t)
                cb_t[t] = makeCB(h_epep[ei][t], h_emem[mi][t],
                                 Form("tmp_cb_e%d_m%d_t%d", ei, mi, t));

            for (size_t w = 0; w < kWindows.size(); ++w) {
                std::array<double, NUM_TRG> Sw{}, Bw{}, sigw{};
                for (int t = 0; t < NUM_TRG; ++t) {
                    const auto [N_all, _eA] = windowIntegralAndErr(
                        h_epem[ei][mi][t], kWindows[w].lo, kWindows[w].hi);
                    const auto [N_cb, _eC]  = windowIntegralAndErr(
                        cb_t[t], kWindows[w].lo, kWindows[w].hi);
                    Sw[t] = N_all - N_cb;
                    Bw[t] = N_cb;
                    const double d = Sw[t] + 2.0 * Bw[t];
                    sigw[t] = (d > 0.0) ? Sw[t] / std::sqrt(d) : 0.0;
                    h_signif[w][t]->SetBinContent(ei + 1, mi + 1, sigw[t]);
                }
                const double sum_sig = sigw[TRG_PT2] + sigw[TRG_PT3];

                if (sum_sig > best[w].signif_sum) {
                    best[w].ep_cut     = cut_min + ei;
                    best[w].em_cut     = cut_min + mi;
                    best[w].S          = Sw;
                    best[w].B          = Bw;
                    best[w].signif     = sigw;
                    best[w].signif_sum = sum_sig;
                }
            }
            for (int t = 0; t < NUM_TRG; ++t) delete cb_t[t];
        }
    }

    // ============================================================================
    // Compute TOP-criterion significance (range 0.14 – 1.4 GeV) at each best
    // ----------------------------------------------------------------------------
    // The TOP-marker is awarded to the (window, trigger) case whose best
    // (ep_cut, em_cut) gives the highest CB-subtracted signal significance
    // in the m_ee range [0.14, 1.4] (above the π⁰ pump). Excluding π⁰ makes
    // the criterion sensitive to the genuinely interesting η-Dalitz and
    // vector-meson regions, instead of being trivially dominated by the
    // huge π⁰ peak.
    // ============================================================================
    {
        for (size_t w = 0; w < kWindows.size(); ++w) {
            auto& b = best[w];
            const int ei = b.ep_cut - cut_min;
            const int mi = b.em_cut - cut_min;
            for (int t = 0; t < NUM_TRG; ++t) {
                TH1D* cb_t = makeCB(h_epep[ei][t], h_emem[mi][t],
                                    Form("tmp_topCB_w%zu_t%d", w, t));
                const auto [N_all, _eA] = windowIntegralAndErr(
                    h_epem[ei][mi][t], kTopLo, kTopHi);
                const auto [N_cb, _eC]  = windowIntegralAndErr(cb_t, kTopLo, kTopHi);
                b.top_S[t] = N_all - N_cb;
                b.top_B[t] = N_cb;
                const double d = b.top_S[t] + 2.0 * b.top_B[t];
                b.top_signif[t] = (d > 0.0) ? b.top_S[t] / std::sqrt(d) : 0.0;
                delete cb_t;
            }
            b.top_signif_sum = b.top_signif[TRG_PT2] + b.top_signif[TRG_PT3];
        }
        double max_top = -1.0;
        size_t tw = 0;
        for (size_t w = 0; w < kWindows.size(); ++w)
            if (best[w].top_signif_sum > max_top) {
                max_top = best[w].top_signif_sum;
                tw = w;
            }
        best[tw].is_top = true;
        std::cout << "TOP criterion (range [" << kTopLo << ", " << kTopHi
                  << "] GeV, sum PT2+PT3): window "
                  << kWindows[tw].short_id
                  << " — top_signif_sum = " << max_top
                  << "  at (ep_cut=" << best[tw].ep_cut
                  << ", em_cut="     << best[tw].em_cut << ")\n";
    }

    // ============================================================================
    // Summary table to stdout — thousands separated with apostrophes
    // ============================================================================
    std::cout << "\n=== Best (ep_cut, em_cut) per window × trigger ===\n";
    std::cout << "  (★ = top sig. in [" << kTopLo << ", " << kTopHi
              << "] GeV across all cases)\n\n";
    auto pad = [](const std::string& s, int w) {
        if ((int)s.size() >= w) return s;
        return std::string(w - s.size(), ' ') + s;
    };
    auto fmtSig = [](double v) {
        std::ostringstream os;
        os << std::fixed << std::setprecision(1) << v;
        return os.str();
    };
    std::cout
        << "    " << std::left << std::setw(28) << "window"
        << std::left << std::setw(5)  << "trg"
        << "ep  em  |  "
        << pad("S",          12) << "  "
        << pad("CB",         12) << "  "
        << pad("win_sig",     9) << "  |  "
        << pad("top_S",      12) << "  "
        << pad("top_CB",     12) << "  "
        << pad("top_sig",     9) << "\n";
    std::cout << "    " << std::string(116, '-') << "\n";
    for (size_t w = 0; w < kWindows.size(); ++w) {
        const auto& b = best[w];
        for (int t = 0; t < NUM_TRG; ++t) {
            const bool first_row = (t == 0);
            std::cout
                << (b.is_top && first_row ? "  ★ " : "    ")
                << std::left << std::setw(28)
                    << (first_row ? kWindows[w].label : std::string{})
                << std::left << std::setw(5)  << kTrgName[t]
                << pad(first_row ? std::to_string(b.ep_cut) : std::string{}, 2)
                << "  "
                << pad(first_row ? std::to_string(b.em_cut) : std::string{}, 2)
                << "  |  "
                << pad(fmtApos(b.S[t]),       12) << "  "
                << pad(fmtApos(b.B[t]),       12) << "  "
                << pad(fmtSig (b.signif[t]),    9) << "  |  "
                << pad(fmtApos(b.top_S[t]),   12) << "  "
                << pad(fmtApos(b.top_B[t]),   12) << "  "
                << pad(fmtSig (b.top_signif[t]),9) << "\n";
        }
        // Per-window summary line
        std::cout
            << "                          sum |  "
            << pad("",                12) << "  "
            << pad("",                12) << "  "
            << pad(fmtSig(b.signif_sum), 9) << "  |  "
            << pad("",                12) << "  "
            << pad("",                12) << "  "
            << pad(fmtSig(b.top_signif_sum), 9) << "\n";
    }
    std::cout << "\n";

    gSystem->mkdir("plots/output", true);

    // ============================================================================
    // OUTPUT 1 — heatmap canvas (2 rows × 4 cols)
    // ============================================================================
    std::cout << "Drawing heatmap canvas...\n";
    TCanvas* c_heat = new TCanvas("c_scan_matchqualnorm_heatmap",
        "Significance landscape (matchqualnorm scan) — 4 windows × 2 triggers",
        1900, 900);
    c_heat->Divide(kWindows.size(), NUM_TRG, 0.001, 0.001);

    for (int t = 0; t < NUM_TRG; ++t) {
        for (size_t w = 0; w < kWindows.size(); ++w) {
            // Row = trigger (top=PT3 for readability), col = window.
            const int row = (NUM_TRG - 1) - t;       // PT3 on top
            const int idx = row * kWindows.size() + w + 1;
            c_heat->cd(idx);
            gPad->SetMargin(0.13, 0.16, 0.13, 0.11);
            TH2D* h = h_signif[w][t];
            h->Draw("COLZ");
            // Mark the COMMON best point (same for PT2 and PT3 rows).
            const auto& b = best[w];
            TMarker* mk = new TMarker(b.ep_cut, b.em_cut, 29);
            mk->SetMarkerSize(b.is_top ? 3.2 : 2.0);
            mk->SetMarkerColor(b.is_top ? kBlack : kRed);
            mk->Draw();
            if (b.is_top) {
                TMarker* halo = new TMarker(b.ep_cut, b.em_cut, 24);
                halo->SetMarkerSize(5.0);
                halo->SetMarkerColor(kYellow + 2);
                halo->Draw();
                TMarker* halo2 = new TMarker(b.ep_cut, b.em_cut, 24);
                halo2->SetMarkerSize(4.0);
                halo2->SetMarkerColor(kOrange + 7);
                halo2->Draw();
                mk->Draw();
            }
            TLatex tx; tx.SetNDC(); tx.SetTextSize(0.038);
            tx.DrawLatex(0.16, 0.92,
                Form("best: (%d, %d)  sig_%s=%.1f",
                     b.ep_cut, b.em_cut,
                     kTrgName[t].c_str(), b.signif[t]));
            tx.SetTextSize(0.032);
            tx.SetTextColor(kGray + 3);
            tx.DrawLatex(0.16, 0.87,
                Form("top_sig_%s [%.2f-%.2f] = %.1f",
                     kTrgName[t].c_str(), kTopLo, kTopHi, b.top_signif[t]));
            if (b.is_top && t == TRG_PT3) {
                // Only annotate TOP in the PT3 row (avoid duplication).
                tx.SetTextColor(kRed + 1);
                tx.SetTextSize(0.040);
                tx.DrawLatex(0.50, 0.92,
                    Form("#bigstar TOP (sum=%.1f)", b.top_signif_sum));
            }
        }
    }
    c_heat->SaveAs("plots/output/scan_rich_matchqualnorm_heatmap_significance.pdf");
    c_heat->SaveAs("plots/output/scan_rich_matchqualnorm_heatmap_significance.png");
    std::cout << "  saved: plots/output/scan_rich_matchqualnorm_heatmap_significance.{pdf,png}\n";

    // ============================================================================
    // OUTPUT 2-5 — best-case canvases (4 panels each: 2x2 layout)
    //
    // Per window: ONE canvas with the SAME (ep_cut, em_cut) for PT2 and PT3.
    //   TL: m_ee at PT3 (all/cb/sig), TR: m_ee at PT2 (all/cb/sig),
    //   BL: 63·PT2/PT3, BR: PT3/(63·PT2) — both ratios computed on SIGNAL
    //   from the per-trigger CB-subtracted spectra at the chosen cuts.
    // ============================================================================
    std::cout << "\nDrawing 4 best-case 2x2 canvases...\n";
    for (size_t w = 0; w < kWindows.size(); ++w) {
        const auto& b = best[w];
        const int ei = b.ep_cut - cut_min;
        const int mi = b.em_cut - cut_min;

        // ---- Build per-trigger histograms (fine 5 MeV → display 20 MeV) ----
        std::array<TH1D*, NUM_TRG> h_all{}, h_cb{}, h_sig{};
        std::array<TH1D*, NUM_TRG> sig_fine{};      // for ratios (5 MeV)
        for (int t = 0; t < NUM_TRG; ++t) {
            TH1D* cb_f = makeCB(h_epep[ei][t], h_emem[mi][t],
                Form("cb_w%zu_t%d_fine", w, t));
            TH1D* sg_f = makeSignal(h_epem[ei][mi][t], cb_f,
                Form("sg_w%zu_t%d_fine", w, t));
            sig_fine[t] = sg_f;
            h_all[t] = rebinFor20(h_epem[ei][mi][t],
                Form("h_all_w%zu_t%d_20", w, t));
            h_cb [t] = rebinFor20(cb_f,
                Form("h_cb_w%zu_t%d_20",  w, t));
            h_sig[t] = rebinFor20(sg_f,
                Form("h_sig_w%zu_t%d_20", w, t));
            delete cb_f;
            styleAll(h_all[t], kBlack,    20);
            styleAll(h_cb [t], kRed + 1,  21);
            styleAll(h_sig[t], kBlue + 1, 22);
        }

        // ---- Ratios at the same best cuts ----
        TH1D* sig_PT2 = rebinFor20(sig_fine[TRG_PT2],
            Form("sig_PT2_best_%s", kWindows[w].short_id.c_str()));
        TH1D* sig_PT3 = rebinFor20(sig_fine[TRG_PT3],
            Form("sig_PT3_best_%s", kWindows[w].short_id.c_str()));
        TH1D* r_corr = makeRatio(sig_PT2, sig_PT3,
            Form("r_corr_best_%s", kWindows[w].short_id.c_str()), 63.0);
        TH1D* r_eff  = makeRatio(sig_PT3, sig_PT2,
            Form("r_eff_best_%s",  kWindows[w].short_id.c_str()), 1.0 / 63.0);
        styleAll(r_corr, kBlue + 1, 20);
        styleAll(r_eff,  kBlue + 1, 20);

        const FitResult f_corr = fitConst(r_corr, kFitLo, kFitHi, kBlue + 2);
        const FitResult f_eff  = fitConst(r_eff,  kFitLo, kFitHi, kBlue + 2);

        // ---- Canvas ----
        TCanvas* c = new TCanvas(
            Form("c_matchqualnorm_best_%s", kWindows[w].short_id.c_str()),
            Form("Best %s (matchqualnorm) — ep_cut=%d, em_cut=%d",
                 kWindows[w].label.c_str(), b.ep_cut, b.em_cut),
            1500, 1000);
        c->Divide(2, 2, 0.001, 0.001);

        // Helper: draw m_ee spectrum panel for a given trigger.
        auto drawMeePanel = [&](int pad_idx, int t) {
            c->cd(pad_idx);
            gPad->SetLogy(true);
            gPad->SetMargin(0.13, 0.04, 0.13, 0.10);
            const double y_max = std::max({h_all[t]->GetMaximum(),
                                           h_cb [t]->GetMaximum(), 1.0});
            h_all[t]->SetTitle(Form(
                "m_{ee} at %s (best cuts ep#leq%d em#leq%d);"
                "M_{e^{+}e^{-}} [GeV/c^{2}];counts / 20 MeV",
                kTrgLabel[t].c_str(), b.ep_cut, b.em_cut));
            h_all[t]->GetYaxis()->SetRangeUser(0.5, y_max * 5.0);
            h_all[t]->Draw("E1");
            h_cb [t]->Draw("E1 SAME");
            h_sig[t]->Draw("E1 SAME");

            // Shade the chosen window.
            TBox* box = new TBox(kWindows[w].lo, 0.5,
                                 kWindows[w].hi, y_max * 5.0);
            box->SetFillColorAlpha(kGreen - 9, 0.20);
            box->SetLineColor(kGreen + 2); box->SetLineStyle(2);
            box->Draw("SAME");
            h_all[t]->Draw("E1 SAME");
            h_cb [t]->Draw("E1 SAME");
            h_sig[t]->Draw("E1 SAME");

            TLegend* leg = new TLegend(0.62, 0.74, 0.95, 0.89);
            leg->SetBorderSize(0); leg->SetFillStyle(0);
            leg->SetTextSize(0.034);
            leg->AddEntry(h_all[t], Form("all (N=%.0f)",
                h_all[t]->Integral()), "lpe");
            leg->AddEntry(h_cb [t], Form("CB (N=%.0f)",
                h_cb [t]->Integral()), "lpe");
            leg->AddEntry(h_sig[t], Form("signal (N=%.0f)",
                h_sig[t]->Integral()), "lpe");
            leg->AddEntry(box,      Form("window: [%.2f, %.2f]",
                kWindows[w].lo, kWindows[w].hi), "f");
            leg->Draw();

            TLatex tex; tex.SetNDC(); tex.SetTextSize(0.034);
            tex.DrawLatex(0.16, 0.92,
                Form("win:  S=%.0f  B=%.0f  S/#sqrt{S+2B}=%.1f",
                     b.S[t], b.B[t], b.signif[t]));
            tex.SetTextSize(0.028);
            tex.SetTextColor(kGray + 3);
            tex.DrawLatex(0.16, 0.88,
                Form("top [%.2f-%.2f]: S=%.0f  B=%.0f  sig=%.1f",
                     kTopLo, kTopHi, b.top_S[t], b.top_B[t], b.top_signif[t]));
            // Bottom-left red TOP banner — only on PT3 panel of the TOP window.
            if (b.is_top && t == TRG_PT3) {
                tex.SetTextColor(kRed + 1);
                tex.SetTextSize(0.046);
                tex.DrawLatex(0.18, 0.22, "#bigstar TOP");
                tex.SetTextSize(0.028);
                tex.DrawLatex(0.18, 0.18,
                    Form("(top_sig_sum in [%.2f,%.2f] = %.1f)",
                         kTopLo, kTopHi, b.top_signif_sum));
            }
        };

        // ----- TL: m_ee PT3 -----
        drawMeePanel(1, TRG_PT3);
        // ----- TR: m_ee PT2 -----
        drawMeePanel(2, TRG_PT2);

        // ----- BL: 63·N_PT2 / N_PT3 correction factor -----
        c->cd(3);
        gPad->SetMargin(0.13, 0.04, 0.13, 0.10);
        r_corr->SetTitle(Form(
            "Trigger correction factor (signal, best cuts);"
            "M_{e^{+}e^{-}} [GeV/c^{2}];63 #upoint N_{PT2}/N_{PT3}"));
        r_corr->GetYaxis()->SetRangeUser(0.0, 5.0);
        r_corr->Draw("E1");
        TLine* ln1 = new TLine(kXmin, 1.0, kXmax, 1.0);
        ln1->SetLineStyle(3); ln1->SetLineColor(kGray + 2); ln1->Draw();
        TLatex tx2; tx2.SetNDC(); tx2.SetTextSize(0.034);
        tx2.SetTextColor(kBlue + 2);
        tx2.DrawLatex(0.16, 0.86,
            Form("const fit [%.2f-%.2f]:  a = %.3f #pm %.3f",
                 kFitLo, kFitHi, f_corr.a, f_corr.e));
        tx2.SetTextSize(0.028);
        tx2.SetTextColor(kGray + 3);
        tx2.DrawLatex(0.16, 0.82,
            Form("#chi^{2}/ndf = %.1f / %d",
                 f_corr.chi2_ndf * f_corr.ndf, f_corr.ndf));

        // ----- BR: N_PT3 / (63·N_PT2) efficiency -----
        c->cd(4);
        gPad->SetMargin(0.13, 0.04, 0.13, 0.10);
        r_eff->SetTitle(Form(
            "Trigger efficiency (signal, best cuts);"
            "M_{e^{+}e^{-}} [GeV/c^{2}];N_{PT3}/(63 #upoint N_{PT2})"));
        r_eff->GetYaxis()->SetRangeUser(0.0, 1.0);
        r_eff->Draw("E1");
        TLine* ln2 = new TLine(kXmin, 1.0, kXmax, 1.0);
        ln2->SetLineStyle(3); ln2->SetLineColor(kGray + 2); ln2->Draw();
        TLatex tx3; tx3.SetNDC(); tx3.SetTextSize(0.034);
        tx3.SetTextColor(kBlue + 2);
        tx3.DrawLatex(0.16, 0.86,
            Form("const fit [%.2f-%.2f]:  a = %.3f #pm %.3f",
                 kFitLo, kFitHi, f_eff.a, f_eff.e));
        tx3.SetTextSize(0.028);
        tx3.SetTextColor(kGray + 3);
        tx3.DrawLatex(0.16, 0.82,
            Form("#chi^{2}/ndf = %.1f / %d",
                 f_eff.chi2_ndf * f_eff.ndf, f_eff.ndf));

        const std::string base = Form(
            "plots/output/scan_rich_matchqualnorm_best_%s",
            kWindows[w].short_id.c_str());
        c->SaveAs((base + ".pdf").c_str());
        c->SaveAs((base + ".png").c_str());
        std::cout << "  saved: " << base << ".{pdf,png}  "
                  << "(ep_cut=" << b.ep_cut
                  << ", em_cut=" << b.em_cut
                  << ", sig_sum=" << b.signif_sum << ")\n";

        for (int t = 0; t < NUM_TRG; ++t) delete sig_fine[t];
    }

    std::cout << "\nDone. 5 output files saved under plots/output/.\n";
}
