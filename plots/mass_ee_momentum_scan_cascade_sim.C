// mass_ee_momentum_scan_cascade_sim.C
// =========================================================================
// 2D momentum-cut scan (ASYMMETRIC ep, em) over the 6-scenario PURITY
// CASCADE on the SIM sample.  For each of the 6 cascade scenarios:
//
//   1. Loop over the (ep_p_cut, em_p_cut) 2D grid (STRICT >, p in MeV/c).
//   2. At each grid point, fill h_PT3 and h_PT2 m_ee histograms weighted
//      by sim_genweight × scenario_cut × (cut_var > grid_val).
//   3. Compute ratio = h_PT2 / h_PT3 with full error propagation.
//   4. Fit a LINEAR function a + b*x on the MID range [0.20, 0.70] GeV.
//      The FOM is chi2/ndf of this MID linear fit.
//   5. Track per-grid (a, b, chi2/ndf, retention); find min chi2/ndf.
//
// Per scenario two canvases are produced:
//   - heatmap of chi2/ndf vs (ep_cut, em_cut)
//   - best-case 1×2 panel (m_ee spectrum left, ratio + linear fit right,
//     ratio Y range [0.5, 2.0]).
//
// IMPORTANT — momentum branches DIFFER per scenario:
//   - scenario 0 (RAW truth) uses dilepton_nt_cor with m_ee_sim AND
//     TRUTH momenta ep_p_sim / em_p_sim.
//   - scenarios 1..5 use dilepton_nt with m_ee AND reco momenta
//     ep_p_rec / em_p_rec.
//
// Scan grid: p_cut ∈ [0, 300] step 15 MeV/c  → 21 values per leg
//             ⇒ 441 (ep, em) combinations per scenario.
//
// Performance: use the "effective-index forward-fill" pattern.  Per event
// compute eff_ep_idx / eff_em_idx (the LARGEST grid index whose value is
// still strictly < the event momentum), fill h_max[eff_ep][eff_em] ONCE,
// then build cumulative h_PT3[ep][em] = sum over (ep' >= ep, em' >= em).
//
// Outputs:
//   plots/output/scan_momentum_cascade_sim_step{0..5}_heatmap.{pdf,png}
//   plots/output/scan_momentum_cascade_sim_step{0..5}_best.{pdf,png}
//
// Usage:
//   root -l -b -q plots/mass_ee_momentum_scan_cascade_sim.C
//   root -l -b -q 'plots/mass_ee_momentum_scan_cascade_sim.C(0, 300, 15)'
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
#include <TBox.h>
#include <TMarker.h>
#include <TStyle.h>
#include <TSystem.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {

constexpr const char* kInputFile = "output_epem_sim.root";

// m_ee binning (20 MeV bins) and MID-range fit window.
constexpr int    kNb    = 70;
constexpr double kXmin  = 0.0;
constexpr double kXmax  = 1.4;
constexpr double kMdLo  = 0.20;
constexpr double kMdHi  = 0.70;

constexpr Color_t kColPT3   = kBlack;
constexpr Color_t kColPT2   = kRed + 1;
constexpr Color_t kColRatio = kBlue + 1;
constexpr Color_t kColLin   = kAzure + 1;

// -------- Trigger encoding --------
enum Trigger : int { TRG_PT2 = 0, TRG_PT3 = 1, NUM_TRG = 2 };

// -------- Scenario record (mirrors cascade reference) --------
struct Scenario {
    std::string id;
    std::string label;
    std::string tree;
    std::string mass_var;
    std::string ep_p_var;
    std::string em_p_var;
    // Functional cut applied to every event in addition to triggers and
    // momentum.  Each lambda receives the live TTreeReaderValue<float>
    // pointers (see set up below) so it can read the relevant branches.
};

// We hard-code the scenario "extra" cuts inside the event loop because the
// scenarios are evaluated against TTreeReader branches; doing it via lambdas
// would require a uniform reader set per tree.  See bookkeeping inside main.

const std::vector<Scenario> kScenarios = {
    {"step0_raw_truth", "0) RAW truth — no cuts",                           "dilepton_nt_cor", "m_ee_sim", "ep_p_sim", "em_p_sim"},
    {"step1_raw",       "1) RAW — no cuts",                                 "dilepton_nt",     "m_ee",     "ep_p_rec", "em_p_rec"},
    {"step2_vertex",    "2) +eVertReco_z > -500",                           "dilepton_nt",     "m_ee",     "ep_p_rec", "em_p_rec"},
    {"step3_isBest",    "3) +isBest == 1",                                  "dilepton_nt",     "m_ee",     "ep_p_rec", "em_p_rec"},
    {"step4_simID",     "4) +ep_sim_id==2 && em_sim_id==3",                 "dilepton_nt",     "m_ee",     "ep_p_rec", "em_p_rec"},
    {"step5_samevtx",   "5) +epem_same_vertex==1 (full old purity gate)",   "dilepton_nt",     "m_ee",     "ep_p_rec", "em_p_rec"},
};

// -------- helpers --------
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

void styleSpec(TH1D* h, Color_t color, Style_t marker) {
    h->SetLineColor(color);
    h->SetMarkerColor(color);
    h->SetMarkerStyle(marker);
    h->SetMarkerSize(0.8);
    h->SetLineWidth(2);
    h->GetXaxis()->SetTitleSize(0.050);
    h->GetYaxis()->SetTitleSize(0.050);
    h->GetXaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetLabelSize(0.045);
    h->GetXaxis()->SetTitleOffset(1.05);
    h->GetYaxis()->SetTitleOffset(1.20);
}

struct LinResult {
    double a = 0.0, a_err = 0.0;
    double b = 0.0, b_err = 0.0;
    double chi2_ndf = std::numeric_limits<double>::infinity();
    int    ndf      = 0;
};

LinResult fitLinRange(TH1D* h, double xlo, double xhi,
                      Color_t color, const char* tag,
                      const char* opt) {
    const std::string fname = std::string(h->GetName()) + "_fl_" + tag;
    TF1* f = new TF1(fname.c_str(), "[0]+[1]*x", xlo, xhi);
    f->SetLineColor(color);
    f->SetLineWidth(3);
    h->Fit(f, opt);
    LinResult r;
    r.a        = f->GetParameter(0);
    r.a_err    = f->GetParError(0);
    r.b        = f->GetParameter(1);
    r.b_err    = f->GetParError(1);
    r.ndf      = f->GetNDF();
    r.chi2_ndf = (r.ndf > 0) ? f->GetChisquare() / r.ndf : 0.0;
    return r;
}

struct BestPoint {
    int    ep_cut       = 0;
    int    em_cut       = 0;
    double a            = 0.0;
    double a_err        = 0.0;
    double b            = 0.0;
    double b_err        = 0.0;
    double chi2_ndf     = std::numeric_limits<double>::infinity();
    int    ndf          = 0;
    double yield_pt3    = 0.0;
    double yield_pt2    = 0.0;
    double retention    = 0.0;   // PT3 yield at cut / PT3 yield at no-cut (0,0)
};

// Build cumulative histogram h_PT[ep][em] from per-(eff_ep,eff_em) h_max.
// Convention (> cut):  event at eff_idx (e_eff, m_eff) passes cuts at
// every (ep_idx, em_idx) with ep_idx <= e_eff AND em_idx <= m_eff.
// So h_PT[ep][em] = sum over (e_eff >= ep, m_eff >= em) of h_max[e_eff][m_eff].
// We use 2D backwards-cumsum (start at top-right, accumulate).
void buildCumulativeBackward(
    const std::vector<std::vector<TH1D*>>& h_max,
    std::vector<std::vector<TH1D*>>& h_pt,
    int n_cuts, int n_bins,
    const std::string& prefix) {
    // Allocate destination hists.
    for (int e = 0; e < n_cuts; ++e) {
        for (int m = 0; m < n_cuts; ++m) {
            h_pt[e][m] = new TH1D(
                Form("%s_e%d_m%d", prefix.c_str(), e, m),
                "", kNb, kXmin, kXmax);
            h_pt[e][m]->Sumw2();
            h_pt[e][m]->SetDirectory(nullptr);
        }
    }
    // Accumulate per mass bin (vectorized over (e,m) using 2D backward sums).
    // For each bin b, value[e][m] = h_max[e][m]'s bin b content;
    // err2 similarly.  Then cumsum backwards.
    std::vector<std::vector<double>> v(n_cuts, std::vector<double>(n_cuts, 0.0));
    std::vector<std::vector<double>> v2(n_cuts, std::vector<double>(n_cuts, 0.0));
    for (int b = 1; b <= n_bins; ++b) {
        for (int e = 0; e < n_cuts; ++e) {
            for (int m = 0; m < n_cuts; ++m) {
                v[e][m]  = h_max[e][m]->GetBinContent(b);
                const double er = h_max[e][m]->GetBinError(b);
                v2[e][m] = er * er;
            }
        }
        // backward cumulative: rows from bottom-up, cols from right-to-left
        // h_pt[e][m] = sum_{e'>=e, m'>=m} h_max[e'][m']
        // first sweep along m: c[e][m] = sum_{m'>=m} v[e][m']
        std::vector<std::vector<double>> c(n_cuts, std::vector<double>(n_cuts, 0.0));
        std::vector<std::vector<double>> c2(n_cuts, std::vector<double>(n_cuts, 0.0));
        for (int e = 0; e < n_cuts; ++e) {
            c[e][n_cuts - 1]  = v[e][n_cuts - 1];
            c2[e][n_cuts - 1] = v2[e][n_cuts - 1];
            for (int m = n_cuts - 2; m >= 0; --m) {
                c[e][m]  = c[e][m + 1]  + v[e][m];
                c2[e][m] = c2[e][m + 1] + v2[e][m];
            }
        }
        // then sweep along e: out[e][m] = sum_{e'>=e} c[e'][m]
        for (int m = 0; m < n_cuts; ++m) {
            double acc = 0.0, acc2 = 0.0;
            for (int e = n_cuts - 1; e >= 0; --e) {
                acc  += c[e][m];
                acc2 += c2[e][m];
                h_pt[e][m]->SetBinContent(b, acc);
                h_pt[e][m]->SetBinError  (b, std::sqrt(acc2));
            }
        }
    }
}

inline int effectiveIdxStrict(float P, int cut_min, int cut_step, int n_cuts) {
    // For convention p > cut: event with momentum P passes all cuts with
    // value strictly < P, i.e. for grid values g_i = cut_min + i*cut_step,
    // passes if g_i < P.  The LARGEST passing i is:
    //   i_max = ceil((P - cut_min)/step) - 1
    int idx = static_cast<int>(
        std::ceil((P - cut_min) / static_cast<double>(cut_step))) - 1;
    if (idx < 0)        return -1;
    if (idx >= n_cuts)  return n_cuts - 1;
    return idx;
}

}  // anonymous namespace

void mass_ee_momentum_scan_cascade_sim(int cut_min = 0,
                                       int cut_max = 300,
                                       int cut_step = 15) {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetErrorX(0.0);
    gStyle->SetPalette(kBird);

    if (cut_step <= 0 || cut_max < cut_min) {
        std::cerr << "Bad scan grid args.\n"; return;
    }
    const int n_cuts = (cut_max - cut_min) / cut_step + 1;
    auto cutVal = [cut_min, cut_step](int i) {
        return cut_min + i * cut_step;
    };

    std::cout << "Momentum scan cascade — sim\n"
              << "  grid: p > cut ∈ [" << cut_min << ", " << cut_max
              << "] step " << cut_step << " MeV/c   ("
              << n_cuts << " values per leg, "
              << n_cuts * n_cuts << " combinations)\n"
              << "  m_ee binning: " << kNb << " × ["
              << kXmin << ", " << kXmax << "]   (= 20 MeV/bin)\n"
              << "  FOM: chi2/ndf of TF1(a + b*x) on ratio in ["
              << kMdLo << ", " << kMdHi << "]\n\n";

    gSystem->mkdir("plots/output", true);

    // ============================================================================
    // Open input file and cache trees
    // ============================================================================
    TFile* fs = TFile::Open(kInputFile, "READ");
    if (!fs || fs->IsZombie()) {
        std::cerr << "Cannot open " << kInputFile << "\n"; return;
    }
    TTree* t_nt     = dynamic_cast<TTree*>(fs->Get("dilepton_nt"));
    TTree* t_nt_cor = dynamic_cast<TTree*>(fs->Get("dilepton_nt_cor"));
    if (!t_nt || !t_nt_cor) {
        std::cerr << "Missing dilepton_nt / dilepton_nt_cor in input.\n"; return;
    }
    std::cout << "Input: " << kInputFile << "\n"
              << "  dilepton_nt     entries = " << t_nt->GetEntries() << "\n"
              << "  dilepton_nt_cor entries = " << t_nt_cor->GetEntries() << "\n\n";

    // ============================================================================
    // PASS A — single sweep through dilepton_nt_cor for scenario 0 (RAW truth).
    // ============================================================================
    std::cout << "PASS A — sweeping dilepton_nt_cor for scenario 0 (truth p)...\n";
    using HistMat = std::vector<std::vector<TH1D*>>;
    HistMat h_max_pt3_s0(n_cuts, std::vector<TH1D*>(n_cuts, nullptr));
    HistMat h_max_pt2_s0(n_cuts, std::vector<TH1D*>(n_cuts, nullptr));
    for (int e = 0; e < n_cuts; ++e) {
        for (int m = 0; m < n_cuts; ++m) {
            h_max_pt3_s0[e][m] = new TH1D(
                Form("h_max_pt3_s0_e%d_m%d", e, m), "", kNb, kXmin, kXmax);
            h_max_pt3_s0[e][m]->Sumw2();
            h_max_pt3_s0[e][m]->SetDirectory(nullptr);
            h_max_pt2_s0[e][m] = new TH1D(
                Form("h_max_pt2_s0_e%d_m%d", e, m), "", kNb, kXmin, kXmax);
            h_max_pt2_s0[e][m]->Sumw2();
            h_max_pt2_s0[e][m]->SetDirectory(nullptr);
        }
    }
    {
        TTreeReader r("dilepton_nt_cor", fs);
        TTreeReaderValue<float> v_m    (r, "m_ee_sim");
        TTreeReaderValue<float> v_pt3  (r, "pt3");
        TTreeReaderValue<float> v_pt2  (r, "pt2");
        TTreeReaderValue<float> v_w    (r, "sim_genweight");
        TTreeReaderValue<float> v_ep_p (r, "ep_p_sim");
        TTreeReaderValue<float> v_em_p (r, "em_p_sim");

        Long64_t n_tot = 0, n_use = 0;
        while (r.Next()) {
            ++n_tot;
            const bool is_pt3 = (static_cast<int>(*v_pt3) == 1);
            const bool is_pt2 = (static_cast<int>(*v_pt2) == 1);
            if (!is_pt3 && !is_pt2) continue;
            const int ei = effectiveIdxStrict(*v_ep_p, cut_min, cut_step, n_cuts);
            if (ei < 0) continue;
            const int mi = effectiveIdxStrict(*v_em_p, cut_min, cut_step, n_cuts);
            if (mi < 0) continue;
            ++n_use;
            const double m = *v_m;
            const double w = *v_w;
            if (is_pt3) h_max_pt3_s0[ei][mi]->Fill(m, w);
            if (is_pt2) h_max_pt2_s0[ei][mi]->Fill(m, w);
        }
        std::cout << "  pass A: scanned " << n_tot << " used " << n_use
                  << " (" << (100.0 * n_use / std::max<Long64_t>(1, n_tot))
                  << " % survival, requires p_sim within bounds)\n";
    }

    // ============================================================================
    // PASS B — single sweep through dilepton_nt for scenarios 1..5 (reco p).
    // We fill 5 parallel h_max tables (one per scenario) since the scenario
    // cuts are cumulative gates and the momentum effective idx is shared.
    // ============================================================================
    std::cout << "PASS B — sweeping dilepton_nt for scenarios 1..5 (reco p)...\n";
    // Index by scenario s = 1..5  ⇒ use 0-based local index [0..4].
    constexpr int kNSReco = 5;
    std::vector<HistMat> h_max_pt3_r(kNSReco,
        HistMat(n_cuts, std::vector<TH1D*>(n_cuts, nullptr)));
    std::vector<HistMat> h_max_pt2_r(kNSReco,
        HistMat(n_cuts, std::vector<TH1D*>(n_cuts, nullptr)));
    for (int s = 0; s < kNSReco; ++s) {
        for (int e = 0; e < n_cuts; ++e) {
            for (int m = 0; m < n_cuts; ++m) {
                h_max_pt3_r[s][e][m] = new TH1D(
                    Form("h_max_pt3_s%d_e%d_m%d", s + 1, e, m), "",
                    kNb, kXmin, kXmax);
                h_max_pt3_r[s][e][m]->Sumw2();
                h_max_pt3_r[s][e][m]->SetDirectory(nullptr);
                h_max_pt2_r[s][e][m] = new TH1D(
                    Form("h_max_pt2_s%d_e%d_m%d", s + 1, e, m), "",
                    kNb, kXmin, kXmax);
                h_max_pt2_r[s][e][m]->Sumw2();
                h_max_pt2_r[s][e][m]->SetDirectory(nullptr);
            }
        }
    }
    {
        TTreeReader r("dilepton_nt", fs);
        TTreeReaderValue<float> v_m     (r, "m_ee");
        TTreeReaderValue<float> v_pt3   (r, "pt3");
        TTreeReaderValue<float> v_pt2   (r, "pt2");
        TTreeReaderValue<float> v_w     (r, "sim_genweight");
        TTreeReaderValue<float> v_ep_p  (r, "ep_p_rec");
        TTreeReaderValue<float> v_em_p  (r, "em_p_rec");
        TTreeReaderValue<float> v_vtxz  (r, "eVertReco_z");
        TTreeReaderValue<float> v_best  (r, "isBest");
        TTreeReaderValue<float> v_epid  (r, "ep_sim_id");
        TTreeReaderValue<float> v_emid  (r, "em_sim_id");
        TTreeReaderValue<float> v_samev (r, "epem_same_vertex");

        Long64_t n_tot = 0, n_use = 0;
        Long64_t pass_s[5] = {0, 0, 0, 0, 0};
        while (r.Next()) {
            ++n_tot;
            const bool is_pt3 = (static_cast<int>(*v_pt3) == 1);
            const bool is_pt2 = (static_cast<int>(*v_pt2) == 1);
            if (!is_pt3 && !is_pt2) continue;
            const int ei = effectiveIdxStrict(*v_ep_p, cut_min, cut_step, n_cuts);
            if (ei < 0) continue;
            const int mi = effectiveIdxStrict(*v_em_p, cut_min, cut_step, n_cuts);
            if (mi < 0) continue;
            ++n_use;
            const double m = *v_m;
            const double w = *v_w;

            // Evaluate cumulative scenario cuts.
            // s1: no extra cut.
            const bool pass1 = true;
            const bool pass2 = pass1 && (*v_vtxz > -500.0f);
            const bool pass3 = pass2 && (static_cast<int>(*v_best) == 1);
            const bool pass4 = pass3 &&
                (static_cast<int>(*v_epid) == 2) &&
                (static_cast<int>(*v_emid) == 3);
            const bool pass5 = pass4 && (static_cast<int>(*v_samev) == 1);
            const bool passes[5] = {pass1, pass2, pass3, pass4, pass5};

            for (int s = 0; s < kNSReco; ++s) {
                if (!passes[s]) continue;
                ++pass_s[s];
                if (is_pt3) h_max_pt3_r[s][ei][mi]->Fill(m, w);
                if (is_pt2) h_max_pt2_r[s][ei][mi]->Fill(m, w);
            }
        }
        std::cout << "  pass B: scanned " << n_tot << " used " << n_use
                  << " (" << (100.0 * n_use / std::max<Long64_t>(1, n_tot))
                  << " % survival)\n";
        for (int s = 0; s < kNSReco; ++s) {
            std::cout << "    scenario " << (s + 1)
                      << " survivors (after cumulative cut, before p-cut): "
                      << pass_s[s] << "\n";
        }
    }

    // ============================================================================
    // Build cumulative histograms per scenario (forward in cuts == backward
    // in eff-idx).
    // ============================================================================
    std::cout << "\nBuilding cumulative PT3/PT2 histograms per scenario...\n";
    std::vector<HistMat> h_PT3(kScenarios.size(),
        HistMat(n_cuts, std::vector<TH1D*>(n_cuts, nullptr)));
    std::vector<HistMat> h_PT2(kScenarios.size(),
        HistMat(n_cuts, std::vector<TH1D*>(n_cuts, nullptr)));

    buildCumulativeBackward(h_max_pt3_s0, h_PT3[0], n_cuts, kNb, "h_PT3_s0");
    buildCumulativeBackward(h_max_pt2_s0, h_PT2[0], n_cuts, kNb, "h_PT2_s0");
    for (int s = 0; s < kNSReco; ++s) {
        buildCumulativeBackward(h_max_pt3_r[s], h_PT3[s + 1], n_cuts, kNb,
            Form("h_PT3_s%d", s + 1));
        buildCumulativeBackward(h_max_pt2_r[s], h_PT2[s + 1], n_cuts, kNb,
            Form("h_PT2_s%d", s + 1));
    }

    // ============================================================================
    // Scan FOM landscape and find best per scenario.
    // ============================================================================
    std::cout << "\nScanning FOM landscape per scenario...\n";

    std::vector<BestPoint> best_per_scenario(kScenarios.size());
    std::vector<LinResult> baseline_per_scenario(kScenarios.size());  // (ep_cut=0, em_cut=0)
    std::vector<double> baseline_yield_pt3(kScenarios.size(), 0.0);
    std::vector<TH2D*> h_chi2(kScenarios.size(), nullptr);

    for (size_t s = 0; s < kScenarios.size(); ++s) {
        const auto& S = kScenarios[s];
        h_chi2[s] = new TH2D(
            Form("h_chi2_%s", S.id.c_str()),
            Form("MID-fit #chi^{2}/ndf — %s;%s > [MeV/c];%s > [MeV/c];#chi^{2}/ndf",
                 S.label.c_str(), S.ep_p_var.c_str(), S.em_p_var.c_str()),
            n_cuts, cut_min - 0.5 * cut_step, cut_max + 0.5 * cut_step,
            n_cuts, cut_min - 0.5 * cut_step, cut_max + 0.5 * cut_step);
        h_chi2[s]->SetDirectory(nullptr);

        // Baseline yield (ep_cut=0, em_cut=0) is the "all pass" PT3 yield used
        // for retention calculations.
        baseline_yield_pt3[s] = h_PT3[s][0][0]->Integral();

        BestPoint best;
        for (int ei = 0; ei < n_cuts; ++ei) {
            for (int mi = 0; mi < n_cuts; ++mi) {
                TH1D* h3 = h_PT3[s][ei][mi];
                TH1D* h2 = h_PT2[s][ei][mi];
                TH1D* r = makeRatio(h2, h3,
                    Form("r_scan_%s_e%d_m%d", S.id.c_str(), ei, mi));
                const LinResult lf = fitLinRange(r, kMdLo, kMdHi, kColLin,
                    "scan", "RQ0");
                h_chi2[s]->SetBinContent(ei + 1, mi + 1, lf.chi2_ndf);

                if (ei == 0 && mi == 0) {
                    baseline_per_scenario[s] = lf;
                }
                if (lf.ndf > 0 && lf.chi2_ndf < best.chi2_ndf) {
                    best.ep_cut    = cutVal(ei);
                    best.em_cut    = cutVal(mi);
                    best.a         = lf.a;
                    best.a_err     = lf.a_err;
                    best.b         = lf.b;
                    best.b_err     = lf.b_err;
                    best.chi2_ndf  = lf.chi2_ndf;
                    best.ndf       = lf.ndf;
                    best.yield_pt3 = h3->Integral();
                    best.yield_pt2 = h2->Integral();
                    best.retention = (baseline_yield_pt3[s] > 0.0)
                        ? best.yield_pt3 / baseline_yield_pt3[s]
                        : 0.0;
                }
                delete r;
            }
        }
        best_per_scenario[s] = best;
    }

    // ============================================================================
    // Save heatmaps and best-case canvases per scenario.
    // ============================================================================
    auto fmt = [](double v, int prec = 2) {
        std::ostringstream os;
        os << std::fixed << std::setprecision(prec) << v;
        return os.str();
    };

    for (size_t s = 0; s < kScenarios.size(); ++s) {
        const auto& S = kScenarios[s];
        const BestPoint& b = best_per_scenario[s];
        const LinResult& base = baseline_per_scenario[s];

        // ---- Heatmap canvas ----
        TCanvas* c_heat = new TCanvas(
            Form("c_mom_scan_cascade_sim_heatmap_%s", S.id.c_str()),
            Form("Momentum scan heatmap — %s", S.label.c_str()),
            900, 800);
        gPad->SetMargin(0.13, 0.16, 0.13, 0.11);
        gPad->SetLogz(true);
        h_chi2[s]->Draw("COLZ");

        TMarker* mk = new TMarker(b.ep_cut, b.em_cut, 29);
        mk->SetMarkerSize(2.4); mk->SetMarkerColor(kRed + 1);
        mk->Draw();

        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.030);
        tx.DrawLatex(0.14, 0.94,
            Form("best: %s>%d, %s>%d  #chi^{2}/ndf=%.2f / %d  (baseline=%.2f)",
                 S.ep_p_var.c_str(), b.ep_cut,
                 S.em_p_var.c_str(), b.em_cut,
                 b.chi2_ndf, b.ndf, base.chi2_ndf));

        const std::string base_h =
            "plots/output/scan_momentum_cascade_sim_" + S.id + "_heatmap";
        c_heat->SaveAs((base_h + ".pdf").c_str());
        c_heat->SaveAs((base_h + ".png").c_str());

        // ---- Best-case 1×2 canvas ----
        const int ei = (b.ep_cut - cut_min) / cut_step;
        const int mi = (b.em_cut - cut_min) / cut_step;
        TH1D* h3 = static_cast<TH1D*>(h_PT3[s][ei][mi]->Clone(
            Form("c_best_PT3_%s", S.id.c_str())));
        TH1D* h2 = static_cast<TH1D*>(h_PT2[s][ei][mi]->Clone(
            Form("c_best_PT2_%s", S.id.c_str())));
        h3->SetDirectory(nullptr); h2->SetDirectory(nullptr);
        TH1D* r = makeRatio(h2, h3,
            Form("c_best_ratio_%s", S.id.c_str()));

        styleSpec(h3, kColPT3, 21);
        styleSpec(h2, kColPT2, 20);
        styleSpec(r , kColRatio, 20);

        // Linear fit for drawing.
        LinResult lf_draw = fitLinRange(r, kMdLo, kMdHi, kColLin, "best", "RQ+");

        TCanvas* c = new TCanvas(
            Form("c_mom_scan_cascade_sim_best_%s", S.id.c_str()),
            Form("Best %s — %s>%d  %s>%d MeV/c",
                 S.label.c_str(),
                 S.ep_p_var.c_str(), b.ep_cut,
                 S.em_p_var.c_str(), b.em_cut),
            1500, 600);
        c->Divide(2, 1, 0.001, 0.001);

        // Left: m_ee spectrum (log Y) at best cuts.
        c->cd(1);
        gPad->SetLogy(true);
        gPad->SetMargin(0.13, 0.04, 0.13, 0.10);
        const double y_max = std::max(h3->GetMaximum(), h2->GetMaximum());
        h3->SetTitle(Form(
            "m_{ee} spectrum at best cuts — %s;%s [GeV/c^{2}];weighted entries / 20 MeV",
            S.label.c_str(), S.mass_var == std::string("m_ee_sim")
                ? "M_{e^{+}e^{-}}^{sim}" : "M_{e^{+}e^{-}}"));
        h3->GetYaxis()->SetRangeUser(std::max(1e-2, y_max * 1e-6), y_max * 5.0);
        h3->Draw("E1");
        h2->Draw("E1 SAME");

        TLegend* leg = new TLegend(0.55, 0.72, 0.95, 0.88);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.030);
        leg->AddEntry(h3, Form("PT3  (#sum w = %.0f)", h3->Integral()), "lpe");
        leg->AddEntry(h2, Form("PT2  (#sum w = %.0f)", h2->Integral()), "lpe");
        leg->Draw();

        TLatex tx2; tx2.SetNDC(); tx2.SetTextSize(0.030);
        tx2.DrawLatex(0.16, 0.90,
            Form("%s > %d MeV/c, %s > %d MeV/c",
                 S.ep_p_var.c_str(), b.ep_cut,
                 S.em_p_var.c_str(), b.em_cut));
        tx2.DrawLatex(0.16, 0.86,
            Form("retention = %.1f %%   (PT3 yield %.0f / baseline %.0f)",
                 100.0 * b.retention, b.yield_pt3, baseline_yield_pt3[s]));

        // Right: ratio with linear fit attached.
        c->cd(2);
        gPad->SetMargin(0.13, 0.04, 0.13, 0.10);
        r->SetTitle(Form(
            "N_{PT2}/N_{PT3} — %s;%s [GeV/c^{2}];N_{PT2}/N_{PT3}",
            S.label.c_str(),
            S.mass_var == std::string("m_ee_sim")
                ? "M_{e^{+}e^{-}}^{sim}" : "M_{e^{+}e^{-}}"));
        r->GetYaxis()->SetRangeUser(0.5, 2.0);
        r->Draw("E1");
        TLine* lref = new TLine(kXmin, 1.0, kXmax, 1.0);
        lref->SetLineStyle(2); lref->SetLineColor(kGray + 2); lref->Draw();

        TLatex tx3; tx3.SetNDC(); tx3.SetTextSize(0.030);
        tx3.SetTextColor(kColLin);
        tx3.DrawLatex(0.16, 0.88,
            Form("[%.2f-%.2f] lin: a=%+.4f #pm %.4f  b=%+.4f #pm %.4f",
                 kMdLo, kMdHi,
                 lf_draw.a, lf_draw.a_err,
                 lf_draw.b, lf_draw.b_err));
        tx3.DrawLatex(0.16, 0.84,
            Form("#chi^{2}/ndf = %.2f / %d = %.2f   (baseline=%.2f)",
                 lf_draw.chi2_ndf * lf_draw.ndf, lf_draw.ndf,
                 lf_draw.chi2_ndf, base.chi2_ndf));

        const std::string base_b =
            "plots/output/scan_momentum_cascade_sim_" + S.id + "_best";
        c->SaveAs((base_b + ".pdf").c_str());
        c->SaveAs((base_b + ".png").c_str());

        std::cout << "  saved: " << base_h << ".{pdf,png}\n"
                  << "  saved: " << base_b << ".{pdf,png}\n";
    }

    // ============================================================================
    // Summary table.
    // ============================================================================
    std::cout << "\n\n=== SUMMARY — momentum scan cascade (sim) ===\n";
    std::cout << "  grid: p > cut ∈ [" << cut_min << ", " << cut_max
              << "] step " << cut_step << " MeV/c\n";
    std::cout << "  FOM:  chi2/ndf of linear fit (a + b*x) on ratio in ["
              << kMdLo << ", " << kMdHi << "]\n\n";

    auto padR = [](const std::string& s, int w) {
        if ((int)s.size() >= w) return s;
        return std::string(w - s.size(), ' ') + s;
    };

    std::cout << "  scenario                                          "
              << "ep_cut   em_cut   |   chi2/ndf   |     a          b        "
              << "| retention | baseline chi2/ndf | flag\n";
    std::cout << "  " << std::string(170, '-') << "\n";
    for (size_t s = 0; s < kScenarios.size(); ++s) {
        const auto& S = kScenarios[s];
        const BestPoint& b = best_per_scenario[s];
        const LinResult& base = baseline_per_scenario[s];
        const bool hit_max_ep = (b.ep_cut == cut_max);
        const bool hit_max_em = (b.em_cut == cut_max);
        const bool hit_min_ep = (b.ep_cut == cut_min);
        const bool hit_min_em = (b.em_cut == cut_min);
        std::string flag;
        if (hit_max_ep || hit_max_em || hit_min_ep || hit_min_em) {
            flag = "[RED FLAG: boundary";
            if (hit_max_ep) flag += " ep=max";
            if (hit_max_em) flag += " em=max";
            if (hit_min_ep) flag += " ep=min";
            if (hit_min_em) flag += " em=min";
            flag += "]";
        }
        std::cout << "  " << std::left << std::setw(50) << S.label
                  << padR(std::to_string(b.ep_cut), 6) << "    "
                  << padR(std::to_string(b.em_cut), 6) << "   |   "
                  << padR(fmt(b.chi2_ndf), 8) << "  |  "
                  << padR(fmt(b.a, 4), 9) << "  "
                  << padR(fmt(b.b, 4), 9) << "  |  "
                  << padR(fmt(100.0 * b.retention, 2) + " %", 9) << "  |  "
                  << padR(fmt(base.chi2_ndf), 8) << "         | "
                  << flag << "\n";
    }
    std::cout << "\n  Note: 'baseline chi2/ndf' = MID linear fit at (ep_cut=0, em_cut=0)\n";
    std::cout << "        (i.e. NO momentum cut applied; only scenario cumulative cut).\n";

    std::cout << "\n=== Cross-check: truth (scenario 0) vs reco (scenario 1) ===\n";
    {
        const BestPoint& b0 = best_per_scenario[0];
        const BestPoint& b1 = best_per_scenario[1];
        std::cout << "  scenario 0 (truth p): best ep=" << b0.ep_cut
                  << " em=" << b0.em_cut << "  chi2/ndf=" << b0.chi2_ndf << "\n";
        std::cout << "  scenario 1 (reco  p): best ep=" << b1.ep_cut
                  << " em=" << b1.em_cut << "  chi2/ndf=" << b1.chi2_ndf << "\n";
        std::cout << "  Δ(ep)=" << (b0.ep_cut - b1.ep_cut)
                  << "  Δ(em)=" << (b0.em_cut - b1.em_cut)
                  << "  Δ(chi2/ndf)=" << (b0.chi2_ndf - b1.chi2_ndf)
                  << "\n";
        const bool boundary300 =
            (b0.ep_cut == 300) || (b0.em_cut == 300) ||
            (b1.ep_cut == 300) || (b1.em_cut == 300);
        std::cout << "  300 MeV/c boundary hit: "
                  << (boundary300 ? "YES — consider extending scan upper bound."
                                  : "no") << "\n";
    }

    fs->Close();
    std::cout << "\nDone. " << (kScenarios.size() * 2)
              << " output canvases saved under plots/output/.\n";
}
