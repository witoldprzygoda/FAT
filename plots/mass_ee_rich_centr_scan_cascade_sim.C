// mass_ee_rich_centr_scan_cascade_sim.C
// =========================================================================
// 2D scan of (ep_rich_centr, em_rich_centr) cuts on the SIM sample across
// the 6-step purity-cascade scenarios.  At each grid point a LINEAR fit
// y = a + b*x is performed on the PT2/PT3 ratio in the MID range
// [0.20, 0.70] GeV/c² — chi²/ndf of that fit is the FOM (lower = flatter).
//
// rich_centr = radial offset of fitted ring center from ideal Cherenkov
//   center.  Smaller value = cleaner, more circular ring; larger value =
//   deformed.  Cut convention (lower = tighter):
//
//     ep_rich_centr <= ep_cut AND em_rich_centr <= em_cut
//
//   Cut range [5, 30] step 1 → 26 values per leg, 676 combinations.
//   Justification (from data, after vertex + isBest gates):
//     mean = 23.3, RMS = 1.66, p50 = 23.3, p95 = 26.0, p99 = 27.0,
//     max ≈ 600.  The scan covers the actionable bulk; cut = 30 keeps
//     essentially all events (~baseline).
//
// 6 scenarios (matched verbatim to the 3piece cascade reference macro):
//   step0_raw_truth  — no cuts, dilepton_nt_cor, m_ee_sim
//   step1_raw        — no cuts, dilepton_nt, m_ee
//   step2_vertex     — eVertReco_z > -500
//   step3_isBest     — + isBest == 1
//   step4_simID      — + ep_sim_id==2 && em_sim_id==3
//   step5_samevtx    — + epem_same_vertex == 1 (full old purity gate)
//
// Per scenario: max-index fill + forward 2D cumulative sum.
//   eff_ep_idx = clamp(ceil(max(0, ep_rich_centr)) - cut_min, [0, n_cuts-1])
//   eff_em_idx = clamp(ceil(max(0, em_rich_centr)) - cut_min, [0, n_cuts-1])
//   skip event if eff > cut_max.
//   h_max[scenario][eff_ep_idx][eff_em_idx][trg] += w
//   After loop:
//     h_PT3[ep_cut_idx][em_cut_idx][trg] = sum over (i<=ep_cut_idx,
//                                                    j<=em_cut_idx) h_max[i][j][trg]
//
// Outputs (24 files + this macro):
//   plots/output/scan_rich_centr_cascade_sim_step{0..5}_heatmap.{pdf,png}
//   plots/output/scan_rich_centr_cascade_sim_step{0..5}_best.{pdf,png}
//
// Usage:
//   root -l -b -q plots/mass_ee_rich_centr_scan_cascade_sim.C
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
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {

constexpr const char* kInputFile = "output_epem_sim.root";

// m_ee binning (20 MeV) — matches cascade reference
constexpr int    kNb    = 70;
constexpr double kXmin  = 0.0;
constexpr double kXmax  = 1.4;

// MID linear-fit range (the FOM range)
constexpr double kFitLo = 0.20;
constexpr double kFitHi = 0.70;

// Trigger encoding
enum Trigger : int { TRG_PT2 = 0, TRG_PT3 = 1, NUM_TRG = 2 };

struct Scenario {
    std::string id;
    std::string label;
    std::string cut;       // baked into the *cut* TString below (Bool kept as float * 1)
    std::string tree;
    std::string mass_var;
    // Branches actually needed (parsed by the event loop). The driver here
    // treats the cut as a runtime predicate over a fixed set of branches —
    // see makeScenarioPasses() below.
};

const std::vector<Scenario> kScenarios = {
    {"step0_raw_truth", "0) RAW truth — no cuts",                          "1", "dilepton_nt_cor", "m_ee_sim"},
    {"step1_raw",       "1) RAW — no cuts",                                "1", "dilepton_nt",     "m_ee"},
    {"step2_vertex",    "2) +eVertReco_z > -500",                          "(eVertReco_z>-500)", "dilepton_nt", "m_ee"},
    {"step3_isBest",    "3) +isBest == 1",                                 "(eVertReco_z>-500)*(isBest==1)", "dilepton_nt", "m_ee"},
    {"step4_simID",     "4) +ep_sim_id==2 && em_sim_id==3",                "(eVertReco_z>-500)*(isBest==1)*(ep_sim_id==2)*(em_sim_id==3)", "dilepton_nt", "m_ee"},
    {"step5_samevtx",   "5) +epem_same_vertex==1 (full old purity gate)",  "(eVertReco_z>-500)*(isBest==1)*(ep_sim_id==2)*(em_sim_id==3)*(epem_same_vertex==1)", "dilepton_nt", "m_ee"},
};

// -------- Best-point bookkeeping per scenario --------
struct BestPoint {
    int    ep_cut    = 0;
    int    em_cut    = 0;
    double a         = 0.0;
    double a_err     = 0.0;
    double b         = 0.0;
    double b_err     = 0.0;
    double chi2_ndf  = std::numeric_limits<double>::infinity();
    int    ndf       = 0;
    double yield_PT3 = 0.0;
    double yield_PT2 = 0.0;
};

// -------- Linear fit helper --------
struct LinFitResult {
    double a, a_err;
    double b, b_err;
    double chi2_ndf;
    int    ndf;
};

LinFitResult fitLinRange(TH1D* h, double xlo, double xhi,
                         Color_t color, const char* tag) {
    const std::string fname = std::string(h->GetName()) + "_fl_" + tag;
    TF1* f = new TF1(fname.c_str(), "[0]+[1]*x", xlo, xhi);
    f->SetLineColor(color);
    f->SetLineWidth(3);
    h->Fit(f, "RQ0");                  // "0" = don't auto-draw
    LinFitResult r;
    r.a        = f->GetParameter(0);
    r.a_err    = f->GetParError(0);
    r.b        = f->GetParameter(1);
    r.b_err    = f->GetParError(1);
    r.ndf      = f->GetNDF();
    r.chi2_ndf = (r.ndf > 0) ? f->GetChisquare() / r.ndf : 0.0;
    return r;
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

}  // anonymous namespace

void mass_ee_rich_centr_scan_cascade_sim(int cut_min = 5, int cut_max = 30) {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetErrorX(0.0);
    gStyle->SetPalette(kBird);

    const int n_cuts = cut_max - cut_min + 1;
    std::cout << "rich_centr 2D scan — cascade (sim, 6 scenarios)\n"
              << "Cut convention: ep_rich_centr <= ep_cut  AND  em_rich_centr <= em_cut\n"
              << "Grid: ep_cut, em_cut ∈ [" << cut_min << ", " << cut_max
              << "] step 1  (n_cuts=" << n_cuts << ", "
              << n_cuts * n_cuts << " combinations per scenario)\n"
              << "FOM: chi2/ndf of LINEAR fit y=a+b*x to PT2/PT3 ratio in ["
              << kFitLo << ", " << kFitHi << "] (LOWER = FLATTER).\n\n";

    TFile* fs = TFile::Open(kInputFile, "READ");
    if (!fs || fs->IsZombie()) {
        std::cerr << "Cannot open " << kInputFile << "\n"; return;
    }
    TTree* nt     = dynamic_cast<TTree*>(fs->Get("dilepton_nt"));
    TTree* nt_cor = dynamic_cast<TTree*>(fs->Get("dilepton_nt_cor"));
    if (!nt || !nt_cor) {
        std::cerr << "Tree(s) missing in " << kInputFile << "\n"; return;
    }
    std::cout << "Input: " << kInputFile << "\n"
              << "  dilepton_nt:     " << nt->GetEntries()     << " entries\n"
              << "  dilepton_nt_cor: " << nt_cor->GetEntries() << " entries\n\n";

    gSystem->mkdir("plots/output", true);

    // ------------------------------------------------------------------------
    // Allocate per-scenario h_max[ei][mi][trg] mass histograms (single fill
    // per event at the *effective* cut index).  After the event loop they
    // will be forward-summed into cumulative h_cum[ep_idx][em_idx][trg]
    // representing the histogram for cut at (ep_idx, em_idx).
    // ------------------------------------------------------------------------
    using HistMat3 = std::vector<std::vector<std::array<TH1D*, NUM_TRG>>>;
    std::vector<HistMat3> h_max(kScenarios.size());
    for (size_t s = 0; s < kScenarios.size(); ++s) {
        h_max[s] = HistMat3(n_cuts,
            std::vector<std::array<TH1D*, NUM_TRG>>(n_cuts));
        for (int i = 0; i < n_cuts; ++i) {
            for (int j = 0; j < n_cuts; ++j) {
                for (int t = 0; t < NUM_TRG; ++t) {
                    h_max[s][i][j][t] = new TH1D(
                        Form("h_max_s%zu_e%d_m%d_t%d", s, i + cut_min, j + cut_min, t),
                        "", kNb, kXmin, kXmax);
                    h_max[s][i][j][t]->Sumw2();
                    h_max[s][i][j][t]->SetDirectory(nullptr);
                }
            }
        }
    }

    // ------------------------------------------------------------------------
    // Event loop pass 1: dilepton_nt_cor for step0 (truth, m_ee_sim).
    //   No analysis cuts apply.
    // ------------------------------------------------------------------------
    Long64_t n_tot_cor = 0, n_use_cor = 0;
    {
        std::cout << "Loop 1/2: dilepton_nt_cor (step0 truth)...\n";
        TTreeReader r("dilepton_nt_cor", fs);
        TTreeReaderValue<float> v_mee    (r, "m_ee_sim");
        TTreeReaderValue<float> v_pt3    (r, "pt3");
        TTreeReaderValue<float> v_pt2    (r, "pt2");
        TTreeReaderValue<float> v_w      (r, "sim_genweight");
        TTreeReaderValue<float> v_ep_rc  (r, "ep_rich_centr");
        TTreeReaderValue<float> v_em_rc  (r, "em_rich_centr");

        while (r.Next()) {
            ++n_tot_cor;
            const bool is_pt3 = (static_cast<int>(*v_pt3) == 1);
            const bool is_pt2 = (static_cast<int>(*v_pt2) == 1);
            if (!is_pt3 && !is_pt2) continue;

            // Effective cut value: smallest integer cut the leg passes
            // under "rich_centr <= cut" (after clipping negatives to 0).
            const int eff_ep = std::max(cut_min,
                static_cast<int>(std::ceil(std::max(0.0f, *v_ep_rc))));
            const int eff_em = std::max(cut_min,
                static_cast<int>(std::ceil(std::max(0.0f, *v_em_rc))));
            if (eff_ep > cut_max || eff_em > cut_max) continue;
            ++n_use_cor;

            const int eff_ep_idx = eff_ep - cut_min;
            const int eff_em_idx = eff_em - cut_min;
            const double m = *v_mee;
            const double w = *v_w;
            // Scenario 0 only — no other cuts apply.
            if (is_pt3) h_max[0][eff_ep_idx][eff_em_idx][TRG_PT3]->Fill(m, w);
            if (is_pt2) h_max[0][eff_ep_idx][eff_em_idx][TRG_PT2]->Fill(m, w);
        }
    }
    std::cout << "  cor: scanned " << n_tot_cor << ", kept " << n_use_cor
              << " (" << (100.0 * n_use_cor /
                          std::max<Long64_t>(1, n_tot_cor)) << " %)\n";

    // ------------------------------------------------------------------------
    // Event loop pass 2: dilepton_nt for steps 1..5 (reco m_ee, cumulative
    // analysis cuts).  Decide per event which scenarios it satisfies.
    // ------------------------------------------------------------------------
    Long64_t n_tot = 0, n_use = 0;
    {
        std::cout << "Loop 2/2: dilepton_nt (steps 1..5 reco)...\n";
        TTreeReader r("dilepton_nt", fs);
        TTreeReaderValue<float> v_mee    (r, "m_ee");
        TTreeReaderValue<float> v_pt3    (r, "pt3");
        TTreeReaderValue<float> v_pt2    (r, "pt2");
        TTreeReaderValue<float> v_w      (r, "sim_genweight");
        TTreeReaderValue<float> v_ep_rc  (r, "ep_rich_centr");
        TTreeReaderValue<float> v_em_rc  (r, "em_rich_centr");
        TTreeReaderValue<float> v_vz     (r, "eVertReco_z");
        TTreeReaderValue<float> v_isBest (r, "isBest");
        TTreeReaderValue<float> v_ep_sid (r, "ep_sim_id");
        TTreeReaderValue<float> v_em_sid (r, "em_sim_id");
        TTreeReaderValue<float> v_samevtx(r, "epem_same_vertex");

        while (r.Next()) {
            ++n_tot;
            const bool is_pt3 = (static_cast<int>(*v_pt3) == 1);
            const bool is_pt2 = (static_cast<int>(*v_pt2) == 1);
            if (!is_pt3 && !is_pt2) continue;

            const int eff_ep = std::max(cut_min,
                static_cast<int>(std::ceil(std::max(0.0f, *v_ep_rc))));
            const int eff_em = std::max(cut_min,
                static_cast<int>(std::ceil(std::max(0.0f, *v_em_rc))));
            if (eff_ep > cut_max || eff_em > cut_max) continue;
            ++n_use;

            const int eff_ep_idx = eff_ep - cut_min;
            const int eff_em_idx = eff_em - cut_min;
            const double m = *v_mee;
            const double w = *v_w;

            // Build scenario pass flags (steps 1..5).  Cuts are cumulative.
            const bool p1 = true;
            const bool p2 = p1 && (*v_vz > -500.0f);
            const bool p3 = p2 && (static_cast<int>(*v_isBest) == 1);
            const bool p4 = p3 && (static_cast<int>(*v_ep_sid) == 2)
                               && (static_cast<int>(*v_em_sid) == 3);
            const bool p5 = p4 && (static_cast<int>(*v_samevtx) == 1);
            const bool pass[6] = {false, p1, p2, p3, p4, p5};

            for (int s = 1; s <= 5; ++s) {
                if (!pass[s]) continue;
                if (is_pt3) h_max[s][eff_ep_idx][eff_em_idx][TRG_PT3]->Fill(m, w);
                if (is_pt2) h_max[s][eff_ep_idx][eff_em_idx][TRG_PT2]->Fill(m, w);
            }
        }
    }
    std::cout << "  reco: scanned " << n_tot << ", kept " << n_use
              << " (" << (100.0 * n_use / std::max<Long64_t>(1, n_tot))
              << " %)\n\n";

    // ------------------------------------------------------------------------
    // Forward 2D cumulative sum per scenario: h_cum[i][j] = sum over
    //   (ip <= i, jm <= j) of h_max[ip][jm].
    //
    // IMPORTANT: do NOT use 2D inclusion-exclusion with subtraction — that
    // inflates Sumw2 errors at every step (each subtraction ADDS errors² in
    // quadrature instead of cancelling them).  Instead build the cumulative
    // via two ordered 1D passes (row-cumulative, then column-cumulative).
    // Each step is a pure ADD which keeps Sumw2 correct since the underlying
    // h_max cells are disjoint (each event fills exactly one (eff_ep, eff_em)).
    // ------------------------------------------------------------------------
    std::cout << "Building cumulative histograms per scenario...\n";
    std::vector<HistMat3> h_cum(kScenarios.size());
    for (size_t s = 0; s < kScenarios.size(); ++s) {
        // Allocate clones of h_max (these will become h_cum).
        h_cum[s] = HistMat3(n_cuts,
            std::vector<std::array<TH1D*, NUM_TRG>>(n_cuts));
        for (int i = 0; i < n_cuts; ++i) {
            for (int j = 0; j < n_cuts; ++j) {
                for (int t = 0; t < NUM_TRG; ++t) {
                    h_cum[s][i][j][t] = static_cast<TH1D*>(
                        h_max[s][i][j][t]->Clone(
                            Form("h_cum_s%zu_e%d_m%d_t%d", s,
                                 i + cut_min, j + cut_min, t)));
                    h_cum[s][i][j][t]->SetDirectory(nullptr);
                }
            }
        }
        // Row-cumulative along j (for each i, accumulate j: cum[i][j] += cum[i][j-1])
        for (int i = 0; i < n_cuts; ++i) {
            for (int j = 1; j < n_cuts; ++j) {
                for (int t = 0; t < NUM_TRG; ++t)
                    h_cum[s][i][j][t]->Add(h_cum[s][i][j - 1][t]);
            }
        }
        // Column-cumulative along i (for each j, accumulate i: cum[i][j] += cum[i-1][j])
        for (int j = 0; j < n_cuts; ++j) {
            for (int i = 1; i < n_cuts; ++i) {
                for (int t = 0; t < NUM_TRG; ++t)
                    h_cum[s][i][j][t]->Add(h_cum[s][i - 1][j][t]);
            }
        }
    }
    std::cout << "  done.\n\n";

    // ------------------------------------------------------------------------
    // For each scenario: compute chi²/ndf landscape, find best, also baseline
    // (loosest = (cut_max, cut_max)).
    // ------------------------------------------------------------------------
    std::vector<TH2D*>      h_chi2(kScenarios.size(), nullptr);
    std::vector<BestPoint>  best(kScenarios.size());
    std::vector<LinFitResult> baseline(kScenarios.size());
    std::vector<double>       baseline_PT3(kScenarios.size(), 0.0);
    std::vector<double>       baseline_PT2(kScenarios.size(), 0.0);

    // Per-scenario baseline yield (= yield at the loosest grid corner).
    // Used both as denominator for retention and as gate for "best" selection.
    for (size_t s = 0; s < kScenarios.size(); ++s) {
        const int top_idx_pre = n_cuts - 1;
        baseline_PT3[s] = h_cum[s][top_idx_pre][top_idx_pre][TRG_PT3]->Integral();
        baseline_PT2[s] = h_cum[s][top_idx_pre][top_idx_pre][TRG_PT2]->Integral();
    }
    // Minimum-retention gate to consider a cell "best".  Very tight cuts
    // produce tiny yields where huge bin errors make chi² nearly zero —
    // mathematically a "perfect" linear fit but physically meaningless.
    constexpr double kMinRetention = 0.10;   // require ≥ 10% of baseline (PT3+PT2)

    for (size_t s = 0; s < kScenarios.size(); ++s) {
        const auto& S = kScenarios[s];
        std::cout << "Scenario " << S.label << " — scanning landscape...\n";

        h_chi2[s] = new TH2D(
            Form("h_chi2_%s", S.id.c_str()),
            Form("Linearity #chi^{2}/ndf — %s;ep_rich_centr #leq;em_rich_centr #leq;#chi^{2}/ndf",
                 S.label.c_str()),
            n_cuts, cut_min - 0.5, cut_max + 0.5,
            n_cuts, cut_min - 0.5, cut_max + 0.5);
        h_chi2[s]->SetDirectory(nullptr);

        const int top_idx = n_cuts - 1;   // (cut_max, cut_max) = loosest = baseline
        const double Y_base_tot = baseline_PT3[s] + baseline_PT2[s];
        const double Y_min      = kMinRetention * Y_base_tot;

        // Number of bins inside the fit range (used to require sufficient
        // data points before accepting a fit as "best").
        const int bin_lo = h_chi2[s]->GetNbinsX();   // placeholder; reset below
        (void)bin_lo;
        const int fit_bin_lo = 1 + static_cast<int>(std::floor(
            (kFitLo - kXmin) / (kXmax - kXmin) * kNb));
        const int fit_bin_hi = static_cast<int>(std::ceil(
            (kFitHi - kXmin) / (kXmax - kXmin) * kNb));
        const int n_fit_bins = fit_bin_hi - fit_bin_lo + 1;
        // Require at least 80% of bins to have a nonzero ratio entry — a
        // very tight cut otherwise leaves only a few high-mass bins and
        // produces a degenerate fit with chi²≈0.
        const int min_nonzero = std::max(5,
            static_cast<int>(0.8 * n_fit_bins));

        for (int ei = 0; ei < n_cuts; ++ei) {
            for (int mi = 0; mi < n_cuts; ++mi) {
                TH1D* h_PT3 = h_cum[s][ei][mi][TRG_PT3];
                TH1D* h_PT2 = h_cum[s][ei][mi][TRG_PT2];
                TH1D* r_ratio = makeRatio(h_PT2, h_PT3,
                    Form("r_tmp_s%zu_e%d_m%d", s, ei, mi));

                // Count non-empty bins in the fit range.
                int n_nonzero = 0;
                for (int b = fit_bin_lo; b <= fit_bin_hi; ++b)
                    if (r_ratio->GetBinError(b) > 0.0) ++n_nonzero;

                const LinFitResult fr = fitLinRange(
                    r_ratio, kFitLo, kFitHi, kBlue + 2, "scan");
                h_chi2[s]->SetBinContent(ei + 1, mi + 1, fr.chi2_ndf);

                const double Y_cell = h_PT3->Integral() + h_PT2->Integral();
                const bool fit_ok = (fr.ndf > 0)
                                  && (n_nonzero >= min_nonzero)
                                  && (Y_cell >= Y_min);
                if (fit_ok && fr.chi2_ndf < best[s].chi2_ndf) {
                    best[s] = BestPoint{
                        cut_min + ei, cut_min + mi,
                        fr.a, fr.a_err, fr.b, fr.b_err,
                        fr.chi2_ndf, fr.ndf,
                        h_PT3->Integral(), h_PT2->Integral()};
                }
                if (ei == top_idx && mi == top_idx) {
                    baseline[s] = fr;   // yields already populated above
                }
                delete r_ratio;
            }
        }
    }

    // ------------------------------------------------------------------------
    // Summary table
    // ------------------------------------------------------------------------
    auto fmtFloat = [](double v, int prec = 3) {
        std::ostringstream os;
        os << std::fixed << std::setprecision(prec) << v;
        return os.str();
    };
    auto pad = [](const std::string& s, int w) {
        if (static_cast<int>(s.size()) >= w) return s;
        return std::string(w - s.size(), ' ') + s;
    };

    std::cout << "\n=== rich_centr scan — best per scenario (min chi²/ndf, lin fit ["
              << kFitLo << ", " << kFitHi << "]) ===\n\n";
    std::cout << "    " << std::left << std::setw(48) << "scenario"
              << "ep  em  |  "
              << pad("a",          10) << "  "
              << pad("b",          10) << "  "
              << pad("chi2/ndf",   12) << "  |  "
              << pad("base chi2/ndf", 14) << "  "
              << pad("retention[%]", 12) << "  flag\n";
    std::cout << "    " << std::string(140, '-') << "\n";
    for (size_t s = 0; s < kScenarios.size(); ++s) {
        const auto& b = best[s];
        const double Y_base_PT3 = baseline_PT3[s];
        const double Y_base_PT2 = baseline_PT2[s];
        const double Y_base_tot = Y_base_PT3 + Y_base_PT2;
        const double Y_best_tot = b.yield_PT3 + b.yield_PT2;
        const double ret = (Y_base_tot > 0.0)
                            ? 100.0 * Y_best_tot / Y_base_tot : 0.0;
        const bool flag_lo = (b.ep_cut == cut_min) || (b.em_cut == cut_min);
        const bool flag_hi = (b.ep_cut == cut_max) || (b.em_cut == cut_max);
        std::string flag;
        if (flag_lo) flag += " RED_FLAG(low_bound)";
        if (flag_hi) flag += " RED_FLAG(high_bound)";
        std::cout
            << "    " << std::left << std::setw(48) << kScenarios[s].label
            << pad(std::to_string(b.ep_cut), 2) << "  "
            << pad(std::to_string(b.em_cut), 2) << "  |  "
            << pad(fmtFloat(b.a, 4),  10) << "  "
            << pad(fmtFloat(b.b, 4),  10) << "  "
            << pad(fmtFloat(b.chi2_ndf, 2) + " / "
                       + std::to_string(b.ndf), 12) << "  |  "
            << pad(fmtFloat(baseline[s].chi2_ndf, 2) + " / "
                       + std::to_string(baseline[s].ndf), 14) << "  "
            << pad(fmtFloat(ret, 2), 12) << flag << "\n";
    }
    std::cout << "\n";

    // ------------------------------------------------------------------------
    // OUTPUT — heatmaps and best-case canvases per scenario
    // ------------------------------------------------------------------------
    for (size_t s = 0; s < kScenarios.size(); ++s) {
        const auto& S = kScenarios[s];
        const auto& b = best[s];

        // --- heatmap ---
        TCanvas* c_h = new TCanvas(
            Form("c_rich_centr_heat_%s", S.id.c_str()),
            Form("rich_centr scan heatmap — %s", S.label.c_str()),
            900, 750);
        gPad->SetMargin(0.13, 0.17, 0.13, 0.11);
        gPad->SetLogz(true);
        h_chi2[s]->Draw("COLZ");
        TMarker* mk = new TMarker(b.ep_cut, b.em_cut, 29);
        mk->SetMarkerSize(2.4); mk->SetMarkerColor(kRed); mk->Draw();
        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.030);
        tx.DrawLatex(0.15, 0.94,
            Form("%s   best: (%d, %d)  #chi^{2}/ndf = %.2f",
                 S.label.c_str(), b.ep_cut, b.em_cut, b.chi2_ndf));
        tx.DrawLatex(0.15, 0.90,
            Form("baseline (no cut): #chi^{2}/ndf = %.2f",
                 baseline[s].chi2_ndf));
        const std::string h_base =
            std::string("plots/output/scan_rich_centr_cascade_sim_") + S.id + "_heatmap";
        c_h->SaveAs((h_base + ".pdf").c_str());
        c_h->SaveAs((h_base + ".png").c_str());
        std::cout << "  saved heatmap: " << h_base << ".{pdf,png}\n";

        // --- best-case 1×2 canvas ---
        const int ei = b.ep_cut - cut_min;
        const int mi = b.em_cut - cut_min;
        TH1D* h_PT3 = static_cast<TH1D*>(h_cum[s][ei][mi][TRG_PT3]->Clone(
            Form("best_PT3_%s", S.id.c_str())));
        TH1D* h_PT2 = static_cast<TH1D*>(h_cum[s][ei][mi][TRG_PT2]->Clone(
            Form("best_PT2_%s", S.id.c_str())));
        h_PT3->SetDirectory(nullptr);
        h_PT2->SetDirectory(nullptr);
        TH1D* r_ratio = makeRatio(h_PT2, h_PT3,
            Form("best_ratio_%s", S.id.c_str()));

        // Attach the linear fit (drawn) plus a draw-only TF1 for the line.
        TF1* f_draw = new TF1(
            Form("f_lin_%s", S.id.c_str()),
            "[0]+[1]*x", kFitLo, kFitHi);
        f_draw->SetLineColor(kBlue + 2);
        f_draw->SetLineWidth(3);
        r_ratio->Fit(f_draw, "RQ0+");

        styleSpec(h_PT3, kBlack,   21);   // filled square
        styleSpec(h_PT2, kRed + 1, 20);   // filled circle
        styleSpec(r_ratio, kBlue + 1, 20);

        TCanvas* c = new TCanvas(
            Form("c_rich_centr_best_%s", S.id.c_str()),
            Form("Best %s — ep_rc<=%d em_rc<=%d",
                 S.label.c_str(), b.ep_cut, b.em_cut),
            1700, 600);
        c->Divide(2, 1, 0.001, 0.001);

        // Left: m_ee spectrum (log Y)
        c->cd(1);
        gPad->SetLogy(true);
        gPad->SetMargin(0.13, 0.04, 0.13, 0.10);
        const double y_max = std::max(h_PT3->GetMaximum(), h_PT2->GetMaximum());
        h_PT3->SetTitle(Form(
            "m_{ee} spectrum — %s (ep#leq%d em#leq%d);M_{e^{+}e^{-}} [GeV/c^{2}];weighted entries / 20 MeV",
            S.label.c_str(), b.ep_cut, b.em_cut));
        h_PT3->GetYaxis()->SetRangeUser(std::max(1e-2, y_max * 1e-6),
                                        y_max * 5.0);
        h_PT3->Draw("E1");
        h_PT2->Draw("E1 SAME");

        TLegend* leg = new TLegend(0.55, 0.74, 0.95, 0.88);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.034);
        leg->AddEntry(h_PT3, Form("PT3  (#sum w = %.0f)", h_PT3->Integral()), "lpe");
        leg->AddEntry(h_PT2, Form("PT2  (#sum w = %.0f)", h_PT2->Integral()), "lpe");
        leg->Draw();

        // Right: ratio + linear fit
        c->cd(2);
        gPad->SetMargin(0.13, 0.04, 0.13, 0.10);
        r_ratio->SetTitle(Form(
            "N_{PT2}/N_{PT3} — %s (lin fit in [%.2f, %.2f]);M_{e^{+}e^{-}} [GeV/c^{2}];N_{PT2}/N_{PT3}",
            S.label.c_str(), kFitLo, kFitHi));
        r_ratio->GetYaxis()->SetRangeUser(0.5, 2.0);
        r_ratio->Draw("E1");
        f_draw->Draw("SAME");
        TLine* lref = new TLine(kXmin, 1.0, kXmax, 1.0);
        lref->SetLineStyle(2); lref->SetLineColor(kGray + 2); lref->Draw();

        TLatex tx2; tx2.SetNDC(); tx2.SetTextSize(0.032);
        tx2.DrawLatex(0.16, 0.88,
            Form("best ep_rc#leq%d  em_rc#leq%d", b.ep_cut, b.em_cut));
        tx2.SetTextColor(kBlue + 2);
        tx2.DrawLatex(0.16, 0.83,
            Form("a = %.4f #pm %.4f   b = %.4f #pm %.4f",
                 b.a, b.a_err, b.b, b.b_err));
        tx2.DrawLatex(0.16, 0.78,
            Form("#chi^{2}/ndf = %.2f / %d = %.3f",
                 b.chi2_ndf * b.ndf, b.ndf, b.chi2_ndf));
        tx2.SetTextColor(kGray + 3);
        tx2.DrawLatex(0.16, 0.73,
            Form("baseline (cut=%d): #chi^{2}/ndf = %.3f",
                 cut_max, baseline[s].chi2_ndf));

        const std::string base =
            std::string("plots/output/scan_rich_centr_cascade_sim_") + S.id + "_best";
        c->SaveAs((base + ".pdf").c_str());
        c->SaveAs((base + ".png").c_str());
        std::cout << "  saved best : " << base << ".{pdf,png}"
                  << "  (ep=" << b.ep_cut << ", em=" << b.em_cut
                  << ", chi2/ndf=" << b.chi2_ndf
                  << ", baseline=" << baseline[s].chi2_ndf << ")\n";
    }

    std::cout << "\nDone. 6 heatmaps + 6 best-case canvases (24 files) "
                 "saved under plots/output/.\n";
    fs->Close();
}
