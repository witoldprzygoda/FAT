// mass_ee_padnum_scan_cascade_sim.C
// =========================================================================
// Cascade purity-step 2-D scan of (ep_rich_padnum, em_rich_padnum) cuts on
// the SIM sample.  Same six cascade scenarios as
//   mass_ee_pt_ratio_purity_cascade_sim_3piece_v3.C
// but, instead of fitting the PT2/PT3 ratio with a three-piece function,
// each scenario does an asymmetric 2-D scan
//
//     ep_rich_padnum >= ep_cut   AND   em_rich_padnum >= em_cut
//     ep_cut, em_cut in {5, 6, ..., 25}   (21 × 21 = 441 combinations)
//
// At every grid point we build PT3 and PT2 mass spectra, form the ratio
// (with proper error propagation), and fit a LINEAR function
//     y = a + b * x
// on the MID range [0.20, 0.70] GeV/c².  The figure-of-merit is the
// chi²/ndf of that linear fit (LOWER = more constant / less modulated).
//
// Performance trick (438 M entries):
//   - Per event we know the LARGEST (ep_cut, em_cut) pair the event still
//     passes:  eff_ep = min(P_ep, cut_max);  eff_em = min(P_em, cut_max).
//     Skip if eff_ep < cut_min OR eff_em < cut_min.
//   - Fill ONE bin of an auxiliary "max" histogram
//       h_max[scen][eff_ep_idx][eff_em_idx][trg]
//   - AFTER the event loop, the cumulative histogram
//       h_cum[ep_cut][em_cut] = sum_{ep' >= ep_cut, em' >= em_cut} h_max[ep'][em']
//     is built by a cheap 2-D forward sum on grid indices.
//   - One pass per tree (dilepton_nt and dilepton_nt_cor) covers all six
//     scenarios because the scenario cuts are evaluated inline per event.
//
// Outputs per scenario (24 image files in total):
//   plots/output/scan_padnum_cascade_sim_step{0..5}_heatmap.{pdf,png}
//   plots/output/scan_padnum_cascade_sim_step{0..5}_best.{pdf,png}
//
// Usage:
//   root -l -b -q plots/mass_ee_padnum_scan_cascade_sim.C
//   root -l -b -q 'plots/mass_ee_padnum_scan_cascade_sim.C(5, 25)'
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

// ------- m_ee binning (same as cascade reference) -------
constexpr int    kNb   = 70;
constexpr double kXmin = 0.0;
constexpr double kXmax = 1.4;

// ------- MID linear-fit range -------
constexpr double kMdLo = 0.20;
constexpr double kMdHi = 0.70;

// ------- Trigger encoding -------
enum Trigger : int { TRG_PT2 = 0, TRG_PT3 = 1, NUM_TRG = 2 };

// ------- Scenario layout — VERBATIM from cascade reference -------
struct Scenario {
    std::string id;
    std::string label;
    std::string cut;
    std::string tree;
    std::string mass_var;
};
const std::vector<Scenario> kScenarios = {
    {"step0_raw_truth", "0) RAW truth -- no cuts",                          "1",                                                                                "dilepton_nt_cor", "m_ee_sim"},
    {"step1_raw",       "1) RAW -- no cuts",                                "1",                                                                                "dilepton_nt",     "m_ee"},
    {"step2_vertex",    "2) +eVertReco_z > -500",                           "(eVertReco_z>-500)",                                                               "dilepton_nt",     "m_ee"},
    {"step3_isBest",    "3) +isBest == 1",                                  "(eVertReco_z>-500)*(isBest==1)",                                                   "dilepton_nt",     "m_ee"},
    {"step4_simID",     "4) +ep_sim_id==2 && em_sim_id==3",                 "(eVertReco_z>-500)*(isBest==1)*(ep_sim_id==2)*(em_sim_id==3)",                     "dilepton_nt",     "m_ee"},
    {"step5_samevtx",   "5) +epem_same_vertex==1 (full old purity gate)",   "(eVertReco_z>-500)*(isBest==1)*(ep_sim_id==2)*(em_sim_id==3)*(epem_same_vertex==1)","dilepton_nt",     "m_ee"},
};

// ------- Helpers -------

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
    double a = 0.0,  a_err = 0.0;
    double b = 0.0,  b_err = 0.0;
    double chi2_ndf = std::numeric_limits<double>::infinity();
    int    ndf      = 0;
};

LinResult fitLinearMid(TH1D* h) {
    const std::string fname = std::string(h->GetName()) + "_flin";
    TF1* f = new TF1(fname.c_str(), "[0]+[1]*x", kMdLo, kMdHi);
    h->Fit(f, "RQ0");
    LinResult r;
    r.a        = f->GetParameter(0);
    r.a_err    = f->GetParError(0);
    r.b        = f->GetParameter(1);
    r.b_err    = f->GetParError(1);
    r.ndf      = f->GetNDF();
    r.chi2_ndf = (r.ndf > 0) ? f->GetChisquare() / r.ndf : 0.0;
    delete f;
    return r;
}

// Best-point bookkeeping per scenario.
struct BestPoint {
    int    ep_cut      = -1;
    int    em_cut      = -1;
    double a           = 0.0,  a_err = 0.0;
    double b           = 0.0,  b_err = 0.0;
    double chi2_ndf    = std::numeric_limits<double>::infinity();
    int    ndf         = 0;
    double retention   = 0.0;        // PT3 yield at best / PT3 yield at (cut_min,cut_min)
    double sum_w_PT3   = 0.0;
    double sum_w_PT2   = 0.0;
};

}  // anonymous namespace

void mass_ee_padnum_scan_cascade_sim(int cut_min = 5, int cut_max = 25) {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetErrorX(0.0);
    gStyle->SetPalette(kBird);

    const int n_cuts = cut_max - cut_min + 1;
    const size_t n_scen = kScenarios.size();
    std::cout << "Cascade padnum scan (sim)\n"
              << "  scenarios   : " << n_scen << "\n"
              << "  cut range   : ep, em in [" << cut_min << ", "
              << cut_max << "] step 1 (" << n_cuts << " values; "
              << n_cuts * n_cuts << " combos)\n"
              << "  FOM         : chi2/ndf of linear fit y=a+bx on PT2/PT3 "
                 "in [" << kMdLo << ", " << kMdHi << "] GeV/c^2\n"
              << "  convention  : padnum >= cut passes\n\n";

    // ============================================================================
    // Allocate per-event "max" histograms
    //   h_max[scen][ei][mi][trg]    : 6 * 21 * 21 * 2 = 5292 hists * 70 bins
    // ============================================================================
    using HMat = std::vector<std::vector<std::vector<std::array<TH1D*, NUM_TRG>>>>;
    HMat h_max(n_scen,
        std::vector<std::vector<std::array<TH1D*, NUM_TRG>>>(n_cuts,
            std::vector<std::array<TH1D*, NUM_TRG>>(n_cuts)));
    for (size_t s = 0; s < n_scen; ++s) {
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

    // ============================================================================
    // Two single-pass event loops (one per tree).  For each event we evaluate
    // each scenario's cut (true/false) and add it to the corresponding
    // (eff_ep_idx, eff_em_idx) bin.
    // ============================================================================
    TFile* fs = TFile::Open("output_epem_sim.root", "READ");
    if (!fs || fs->IsZombie()) {
        std::cerr << "Cannot open output_epem_sim.root\n"; return;
    }

    // Map: which scenarios use which tree.
    std::vector<size_t> idx_cor, idx_nt;
    for (size_t s = 0; s < n_scen; ++s) {
        if (kScenarios[s].tree == "dilepton_nt_cor") idx_cor.push_back(s);
        else                                          idx_nt.push_back(s);
    }

    // ---- Pass 1: dilepton_nt_cor (m_ee_sim) ----
    if (!idx_cor.empty()) {
        std::cout << "Reading dilepton_nt_cor (m_ee_sim) ...\n";
        TTreeReader r("dilepton_nt_cor", fs);
        TTreeReaderValue<float> v_m  (r, "m_ee_sim");
        TTreeReaderValue<float> v_pt3(r, "pt3");
        TTreeReaderValue<float> v_pt2(r, "pt2");
        TTreeReaderValue<float> v_w  (r, "sim_genweight");
        TTreeReaderValue<float> v_ep (r, "ep_rich_padnum");
        TTreeReaderValue<float> v_em (r, "em_rich_padnum");

        Long64_t n_tot = 0, n_kept = 0;
        while (r.Next()) {
            ++n_tot;
            const bool is_pt3 = (static_cast<int>(*v_pt3) == 1);
            const bool is_pt2 = (static_cast<int>(*v_pt2) == 1);
            if (!is_pt3 && !is_pt2) continue;

            // Effective grid indices: for >= cuts, event with padnum P passes
            // all cuts c <= P.  Cap at cut_max; demote to -1 if below cut_min.
            const int P_ep = static_cast<int>(*v_ep);
            const int P_em = static_cast<int>(*v_em);
            const int eff_ep_idx = std::min(cut_max, P_ep) - cut_min;
            const int eff_em_idx = std::min(cut_max, P_em) - cut_min;
            if (eff_ep_idx < 0 || eff_em_idx < 0) continue;

            const double m = *v_m;
            const double w = *v_w;
            ++n_kept;

            // step0_raw_truth: cut == "1" → always passes.  Other cor-tree
            // scenarios (none in this default list) would be added here.
            for (size_t s : idx_cor) {
                // All current cor-tree scenarios have cut == "1".  Keep the
                // structure flexible — explicit `true` here.
                bool pass = true;
                if (kScenarios[s].cut != "1") {
                    // No other cor-tree scenario in the default list; bail.
                    pass = false;
                }
                if (!pass) continue;
                if (is_pt3) h_max[s][eff_ep_idx][eff_em_idx][TRG_PT3]->Fill(m, w);
                if (is_pt2) h_max[s][eff_ep_idx][eff_em_idx][TRG_PT2]->Fill(m, w);
            }
        }
        std::cout << "  scanned " << n_tot << ", kept " << n_kept << "\n";
    }

    // ---- Pass 2: dilepton_nt (m_ee, all cascade quality cuts) ----
    if (!idx_nt.empty()) {
        std::cout << "Reading dilepton_nt (m_ee) ...\n";
        TTreeReader r("dilepton_nt", fs);
        TTreeReaderValue<float> v_m      (r, "m_ee");
        TTreeReaderValue<float> v_pt3    (r, "pt3");
        TTreeReaderValue<float> v_pt2    (r, "pt2");
        TTreeReaderValue<float> v_w      (r, "sim_genweight");
        TTreeReaderValue<float> v_ep     (r, "ep_rich_padnum");
        TTreeReaderValue<float> v_em     (r, "em_rich_padnum");
        TTreeReaderValue<float> v_vz     (r, "eVertReco_z");
        TTreeReaderValue<float> v_isbest (r, "isBest");
        TTreeReaderValue<float> v_ep_sid (r, "ep_sim_id");
        TTreeReaderValue<float> v_em_sid (r, "em_sim_id");
        TTreeReaderValue<float> v_samevx (r, "epem_same_vertex");

        Long64_t n_tot = 0, n_kept = 0;
        while (r.Next()) {
            ++n_tot;
            const bool is_pt3 = (static_cast<int>(*v_pt3) == 1);
            const bool is_pt2 = (static_cast<int>(*v_pt2) == 1);
            if (!is_pt3 && !is_pt2) continue;

            const int P_ep = static_cast<int>(*v_ep);
            const int P_em = static_cast<int>(*v_em);
            const int eff_ep_idx = std::min(cut_max, P_ep) - cut_min;
            const int eff_em_idx = std::min(cut_max, P_em) - cut_min;
            if (eff_ep_idx < 0 || eff_em_idx < 0) continue;

            // Pre-compute cascade booleans once.
            const bool c_vertex  = (*v_vz > -500.0f);
            const bool c_isbest  = (static_cast<int>(*v_isbest) == 1);
            const bool c_simid   = (static_cast<int>(*v_ep_sid) == 2) &&
                                   (static_cast<int>(*v_em_sid) == 3);
            const bool c_samevtx = (static_cast<int>(*v_samevx) == 1);

            const double m = *v_m;
            const double w = *v_w;
            ++n_kept;

            for (size_t s : idx_nt) {
                bool pass = true;
                const std::string& id = kScenarios[s].id;
                if      (id == "step1_raw")     pass = true;
                else if (id == "step2_vertex")  pass = c_vertex;
                else if (id == "step3_isBest")  pass = c_vertex && c_isbest;
                else if (id == "step4_simID")   pass = c_vertex && c_isbest && c_simid;
                else if (id == "step5_samevtx") pass = c_vertex && c_isbest && c_simid && c_samevtx;
                if (!pass) continue;
                if (is_pt3) h_max[s][eff_ep_idx][eff_em_idx][TRG_PT3]->Fill(m, w);
                if (is_pt2) h_max[s][eff_ep_idx][eff_em_idx][TRG_PT2]->Fill(m, w);
            }
        }
        std::cout << "  scanned " << n_tot << ", kept " << n_kept << "\n";
    }

    fs->Close();

    // ============================================================================
    // 2-D forward-sum cumulative histograms:
    //   h_cum[scen][ep_idx][em_idx][trg](mass_bin) =
    //       sum over (ep' >= ep_idx, em' >= em_idx) of h_max[ep'][em'][trg](mass_bin)
    //
    // IMPORTANT: do NOT use 2D inclusion-exclusion with subtraction — that
    // inflates Sumw2 errors at every step (each subtraction ADDS errors² in
    // quadrature instead of cancelling them).  Instead build the cumulative
    // via two ordered 1D passes (row-cumulative, then column-cumulative).
    // Each step is a pure ADD which keeps Sumw2 correct since the underlying
    // h_max cells are disjoint (each event fills exactly one (eff_ep, eff_em)).
    // For >= convention we accumulate in DECREASING index order.
    // ============================================================================
    std::cout << "Building cumulative histograms (2-D forward sum)...\n";
    HMat h_cum(n_scen,
        std::vector<std::vector<std::array<TH1D*, NUM_TRG>>>(n_cuts,
            std::vector<std::array<TH1D*, NUM_TRG>>(n_cuts)));
    for (size_t s = 0; s < n_scen; ++s) {
        // Allocate clones of h_max (these will become h_cum).
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
        // Row-cumulative along j (decreasing): for each i, cum[i][j] += cum[i][j+1]
        for (int i = 0; i < n_cuts; ++i) {
            for (int j = n_cuts - 2; j >= 0; --j) {
                for (int t = 0; t < NUM_TRG; ++t)
                    h_cum[s][i][j][t]->Add(h_cum[s][i][j + 1][t]);
            }
        }
        // Column-cumulative along i (decreasing): for each j, cum[i][j] += cum[i+1][j]
        for (int j = 0; j < n_cuts; ++j) {
            for (int i = n_cuts - 2; i >= 0; --i) {
                for (int t = 0; t < NUM_TRG; ++t)
                    h_cum[s][i][j][t]->Add(h_cum[s][i + 1][j][t]);
            }
        }
    }

    // ============================================================================
    // Per-scenario landscape: chi2/ndf of linear fit, plus retention.
    // ============================================================================
    std::cout << "\nComputing chi2/ndf landscape per scenario...\n";
    gSystem->mkdir("plots/output", true);

    std::vector<TH2D*> h_chi2(n_scen, nullptr);
    std::vector<TH2D*> h_ret (n_scen, nullptr);
    std::vector<BestPoint> best(n_scen);
    std::vector<LinResult> base_fit(n_scen);  // baseline @ (cut_min, cut_min)
    std::vector<double>    base_sumPT3(n_scen, 0.0);

    for (size_t s = 0; s < n_scen; ++s) {
        h_chi2[s] = new TH2D(
            Form("h_chi2_%s", kScenarios[s].id.c_str()),
            Form("%s;ep_rich_padnum #geq;em_rich_padnum #geq;#chi^{2}/ndf (linear mid fit)",
                 kScenarios[s].label.c_str()),
            n_cuts, cut_min - 0.5, cut_max + 0.5,
            n_cuts, cut_min - 0.5, cut_max + 0.5);
        h_chi2[s]->SetDirectory(nullptr);

        h_ret[s] = new TH2D(
            Form("h_ret_%s", kScenarios[s].id.c_str()),
            Form("%s;ep_rich_padnum #geq;em_rich_padnum #geq;PT3 retention",
                 kScenarios[s].label.c_str()),
            n_cuts, cut_min - 0.5, cut_max + 0.5,
            n_cuts, cut_min - 0.5, cut_max + 0.5);
        h_ret[s]->SetDirectory(nullptr);

        // Baseline reference = (cut_min, cut_min)
        base_sumPT3[s] = h_cum[s][0][0][TRG_PT3]->Integral();

        for (int ei = 0; ei < n_cuts; ++ei) {
            for (int mi = 0; mi < n_cuts; ++mi) {
                TH1D* hPT3 = h_cum[s][ei][mi][TRG_PT3];
                TH1D* hPT2 = h_cum[s][ei][mi][TRG_PT2];
                TH1D* r = makeRatio(hPT2, hPT3,
                    Form("r_tmp_s%zu_e%d_m%d", s, ei, mi));
                const LinResult lf = fitLinearMid(r);

                const double y3 = hPT3->Integral();
                const double ret = (base_sumPT3[s] > 0.0)
                                       ? (y3 / base_sumPT3[s]) : 0.0;
                h_chi2[s]->SetBinContent(ei + 1, mi + 1, lf.chi2_ndf);
                h_ret [s]->SetBinContent(ei + 1, mi + 1, ret);

                // Track baseline (cut_min, cut_min) for the summary.
                if (ei == 0 && mi == 0) base_fit[s] = lf;

                if (lf.ndf > 0 && lf.chi2_ndf < best[s].chi2_ndf) {
                    best[s] = BestPoint{
                        cut_min + ei, cut_min + mi,
                        lf.a, lf.a_err, lf.b, lf.b_err,
                        lf.chi2_ndf, lf.ndf,
                        ret,
                        y3, hPT2->Integral()};
                }
                delete r;
            }
        }
    }

    // ============================================================================
    // Summary table (per scenario: best + baseline + red-flag check)
    // ============================================================================
    auto fmtF = [](double v, int p = 3) {
        std::ostringstream os; os << std::fixed << std::setprecision(p) << v;
        return os.str();
    };

    std::cout << "\n=== Cascade padnum scan — best per scenario ===\n"
              << "  FOM = chi2/ndf of linear y=a+bx on PT2/PT3 in ["
              << kMdLo << ", " << kMdHi << "] GeV/c^2\n"
              << "  Baseline    = (ep_cut=" << cut_min
              << ", em_cut="     << cut_min << ")  [no cut tightening]\n"
              << "  red flag    = best (ep,em) hits scan boundary ["
              << cut_min << " or " << cut_max << "]\n\n";

    std::cout << std::left << std::setw(48) << "scenario"
              << "ep  em  | "
              << std::setw(11) << "chi2/ndf"
              << std::setw(11) << "(base)"
              << std::setw(14) << "a +- err"
              << std::setw(14) << "b +- err"
              << std::setw(10) << "retent"
              << "flag\n";
    std::cout << std::string(140, '-') << "\n";
    for (size_t s = 0; s < n_scen; ++s) {
        const auto& b = best[s];
        const bool red = (b.ep_cut == cut_min || b.ep_cut == cut_max ||
                          b.em_cut == cut_min || b.em_cut == cut_max);
        std::cout
            << std::left << std::setw(48) << kScenarios[s].label
            << std::right << std::setw(2) << b.ep_cut << "  "
            << std::setw(2) << b.em_cut << "  | "
            << std::left << std::setw(11) << fmtF(b.chi2_ndf, 3)
            << std::setw(11) << fmtF(base_fit[s].chi2_ndf, 3)
            << std::setw(14) << (fmtF(b.a, 3) + "+-" + fmtF(b.a_err, 3))
            << std::setw(14) << (fmtF(b.b, 3) + "+-" + fmtF(b.b_err, 3))
            << std::setw(10) << (fmtF(100.0 * b.retention, 1) + "%")
            << (red ? "RED FLAG (boundary)" : "")
            << "\n";
    }
    std::cout << "\n";

    // ============================================================================
    // Per-scenario canvases: (1) heatmap of chi2/ndf, (2) best-case 1x2
    // ============================================================================
    for (size_t s = 0; s < n_scen; ++s) {
        const auto& S  = kScenarios[s];
        const auto& bp = best[s];

        // ---------- (1) Heatmap ----------
        TCanvas* c_h = new TCanvas(
            Form("c_heatmap_%s", S.id.c_str()),
            Form("Heatmap chi2/ndf — %s", S.label.c_str()),
            900, 800);
        gPad->SetMargin(0.13, 0.16, 0.13, 0.10);
        gPad->SetLogz(true);
        h_chi2[s]->Draw("COLZ");

        TMarker* mk = new TMarker(bp.ep_cut, bp.em_cut, 29);
        mk->SetMarkerSize(2.5);
        mk->SetMarkerColor(kRed);
        mk->Draw();

        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.030);
        tx.DrawLatex(0.16, 0.92,
            Form("best: (%d, %d)  #chi^{2}/ndf = %.3f  (baseline %.3f)",
                 bp.ep_cut, bp.em_cut, bp.chi2_ndf, base_fit[s].chi2_ndf));
        tx.SetTextSize(0.026);
        tx.DrawLatex(0.16, 0.88,
            Form("retention vs (%d,%d): %.1f%%   #sum w_{PT3} = %.0f",
                 cut_min, cut_min, 100.0 * bp.retention, bp.sum_w_PT3));
        const bool red = (bp.ep_cut == cut_min || bp.ep_cut == cut_max ||
                          bp.em_cut == cut_min || bp.em_cut == cut_max);
        if (red) {
            tx.SetTextColor(kRed + 1); tx.SetTextSize(0.034);
            tx.DrawLatex(0.55, 0.92, "RED FLAG: boundary");
        }

        const std::string base_h = "plots/output/scan_padnum_cascade_sim_"
                                   + S.id + "_heatmap";
        c_h->SaveAs((base_h + ".pdf").c_str());
        c_h->SaveAs((base_h + ".png").c_str());
        std::cout << "  saved: " << base_h << ".{pdf,png}\n";

        // ---------- (2) Best-case 1x2 ----------
        const int ei = bp.ep_cut - cut_min;
        const int mi = bp.em_cut - cut_min;
        TH1D* hPT3 = static_cast<TH1D*>(h_cum[s][ei][mi][TRG_PT3]->Clone(
            Form("c_PT3_%s", S.id.c_str())));
        TH1D* hPT2 = static_cast<TH1D*>(h_cum[s][ei][mi][TRG_PT2]->Clone(
            Form("c_PT2_%s", S.id.c_str())));
        hPT3->SetDirectory(nullptr);
        hPT2->SetDirectory(nullptr);
        TH1D* r = makeRatio(hPT2, hPT3,
            Form("c_r_%s", S.id.c_str()));

        // Attach the linear fit to the ratio histogram so it is drawn.
        TF1* f = new TF1(Form("f_lin_%s", S.id.c_str()),
                         "[0]+[1]*x", kMdLo, kMdHi);
        f->SetLineColor(kBlue + 2);
        f->SetLineWidth(3);
        r->Fit(f, "RQ+");

        styleSpec(hPT3, kBlack,   21);
        styleSpec(hPT2, kRed + 1, 20);
        styleSpec(r,    kBlue + 1, 20);

        TCanvas* c_b = new TCanvas(
            Form("c_best_%s", S.id.c_str()),
            Form("Best %s — ep>=%d em>=%d", S.label.c_str(),
                 bp.ep_cut, bp.em_cut),
            1500, 600);
        c_b->Divide(2, 1, 0.001, 0.001);

        // Left: m_ee spectrum (log Y)
        c_b->cd(1);
        gPad->SetLogy(true);
        gPad->SetMargin(0.13, 0.04, 0.13, 0.10);
        const double y_max = std::max(hPT3->GetMaximum(), hPT2->GetMaximum());
        hPT3->SetTitle(Form(
            "m_{ee} spectrum — %s  (ep#geq%d em#geq%d);"
            "M_{e^{+}e^{-}} [GeV/c^{2}];weighted entries / 20 MeV",
            S.label.c_str(), bp.ep_cut, bp.em_cut));
        hPT3->GetYaxis()->SetRangeUser(std::max(1e-2, y_max * 1e-6),
                                       y_max * 5.0);
        hPT3->Draw("E1");
        hPT2->Draw("E1 SAME");
        TLegend* leg = new TLegend(0.55, 0.74, 0.95, 0.88);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.034);
        leg->AddEntry(hPT3, Form("PT3  (#sum w = %.0f)", hPT3->Integral()), "lpe");
        leg->AddEntry(hPT2, Form("PT2  (#sum w = %.0f)", hPT2->Integral()), "lpe");
        leg->Draw();

        // Right: ratio + linear fit
        c_b->cd(2);
        gPad->SetMargin(0.13, 0.04, 0.13, 0.10);
        r->SetTitle(Form(
            "N_{PT2}/N_{PT3} — %s;M_{e^{+}e^{-}} [GeV/c^{2}];N_{PT2}/N_{PT3}",
            S.label.c_str()));
        r->GetYaxis()->SetRangeUser(0.5, 2.0);
        r->Draw("E1");

        TLine* lref = new TLine(kXmin, 1.0, kXmax, 1.0);
        lref->SetLineStyle(2); lref->SetLineColor(kGray + 2); lref->Draw();

        TLatex tx2; tx2.SetNDC(); tx2.SetTextSize(0.030);
        tx2.DrawLatex(0.16, 0.90,
            Form("best ep#geq%d em#geq%d   retention = %.1f%%",
                 bp.ep_cut, bp.em_cut, 100.0 * bp.retention));
        tx2.SetTextColor(kBlue + 2);
        tx2.DrawLatex(0.16, 0.85,
            Form("[%.2f, %.2f] linear:  a = %.4f #pm %.4f   "
                 "b = %.4f #pm %.4f", kMdLo, kMdHi, bp.a, bp.a_err,
                 bp.b, bp.b_err));
        tx2.DrawLatex(0.16, 0.81,
            Form("#chi^{2}/ndf = %.3f / %d = %.3f  (baseline %.3f)",
                 bp.chi2_ndf * bp.ndf, bp.ndf, bp.chi2_ndf,
                 base_fit[s].chi2_ndf));

        const std::string base_b = "plots/output/scan_padnum_cascade_sim_"
                                   + S.id + "_best";
        c_b->SaveAs((base_b + ".pdf").c_str());
        c_b->SaveAs((base_b + ".png").c_str());
        std::cout << "  saved: " << base_b << ".{pdf,png}\n";
    }

    std::cout << "\nDone. " << (n_scen * 2 * 2)
              << " output files saved under plots/output/.\n";
}
