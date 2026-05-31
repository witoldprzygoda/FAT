// mass_ee_matchqualnorm_scan_cascade_sim.C
// =========================================================================
// 2D scan of (ep_richmatchqualitynorm, em_richmatchqualitynorm) cuts on the
// SIM sample, evaluated across the 6 purity-cascade scenarios used in
// mass_ee_pt_ratio_purity_cascade_sim_3piece_v3.C.
//
// For EACH scenario:
//   • Scan a 2D grid (asymmetric) of ep_cut, em_cut on [0, 100] step 5
//     (21 values per leg, 441 combinations).
//   • Cut convention: ep_richmatchqualitynorm <= ep_cut AND
//                     em_richmatchqualitynorm <= em_cut   (lower = tighter)
//   • At each grid point compute h_PT3 and h_PT2 weighted by
//        sim_genweight × scenario_cut × cut_var_cut
//     and form the ratio r = h_PT2 / h_PT3 with full error propagation.
//   • FOM = chi²/ndf of LINEAR fit y = a + b·x on the MID range
//        [0.20, 0.70] GeV/c²  (option "RQ0").
//   • Track (a, b, chi²/ndf, retention) per grid point and report the
//     best (minimum chi²/ndf) cut combination.
//
// Outputs per scenario (12 files; 24 total incl. png twins):
//   plots/output/scan_matchqualnorm_cascade_sim_step{0..5}_heatmap.{pdf,png}
//   plots/output/scan_matchqualnorm_cascade_sim_step{0..5}_best.{pdf,png}
//
// Performance pattern:
//   - One pass over each input tree, computing the "effective cut index"
//     for ep / em. For <= convention, an event with quality Q passes all
//     cuts >= ceil(Q). We bin Q into the integer-step grid and fill
//     h_max[scenario][eff_ep_idx][eff_em_idx]. After the loop we build
//     cumulative h_PT3 / h_PT2 by summing over (ep'>=ep, em'>=em). This
//     avoids 441-way fills per event.
//   - For step0_raw_truth the macro reads dilepton_nt_cor / m_ee_sim;
//     all other scenarios use dilepton_nt / m_ee with the scenario gate
//     applied at fill time.
//
// Usage:
//   root -l -b -q plots/mass_ee_matchqualnorm_scan_cascade_sim.C
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

// -------- m_ee binning and MID fit range --------
constexpr int    kNb    = 70;
constexpr double kXmin  = 0.0;
constexpr double kXmax  = 1.4;
constexpr double kFitLo = 0.20;
constexpr double kFitHi = 0.70;

// -------- 2D scan grid (asymmetric ep / em) --------
constexpr int    kCutMin  = 0;
constexpr int    kCutMax  = 100;
constexpr int    kCutStep = 5;
constexpr int    kNCuts   = (kCutMax - kCutMin) / kCutStep + 1;  // 21

// -------- Scenarios (verbatim from the cascade reference macro) --------
struct Scenario {
    std::string id;
    std::string label;
    enum Sel : int {
        SEL_RAW_TRUTH = 0,   // step0  no cuts, dilepton_nt_cor, m_ee_sim
        SEL_RAW       = 1,   // step1  no cuts, dilepton_nt,     m_ee
        SEL_VTX       = 2,   // step2  +eVertReco_z>-500
        SEL_BEST      = 3,   // step3  +isBest==1
        SEL_SIMID     = 4,   // step4  +ep_sim_id==2 && em_sim_id==3
        SEL_SAMEVTX   = 5    // step5  +epem_same_vertex==1
    };
    Sel sel;
};

const std::vector<Scenario> kScenarios = {
    {"step0_raw_truth", "0) RAW truth - no cuts",                         Scenario::SEL_RAW_TRUTH},
    {"step1_raw",       "1) RAW - no cuts",                               Scenario::SEL_RAW},
    {"step2_vertex",    "2) +eVertReco_z > -500",                         Scenario::SEL_VTX},
    {"step3_isBest",    "3) +isBest == 1",                                Scenario::SEL_BEST},
    {"step4_simID",     "4) +ep_sim_id==2 && em_sim_id==3",               Scenario::SEL_SIMID},
    {"step5_samevtx",   "5) +epem_same_vertex==1 (full old purity gate)", Scenario::SEL_SAMEVTX},
};

constexpr int NUM_TRG = 2;
enum Trigger : int { TRG_PT2 = 0, TRG_PT3 = 1 };

// -------- Best-point bookkeeping --------
struct BestPoint {
    int    ep_cut    = 0;
    int    em_cut    = 0;
    double a         = 0.0;
    double a_err     = 0.0;
    double b         = 0.0;
    double b_err     = 0.0;
    double chi2_ndf  = std::numeric_limits<double>::infinity();
    int    ndf       = 0;
    double retention = 0.0;
};

// -------- Linear-fit helper --------
struct FitResult {
    double a;
    double a_err;
    double b;
    double b_err;
    double chi2_ndf;
    int    ndf;
    bool   ok;
};

FitResult fitLinRange(TH1D* h, double xlo, double xhi) {
    FitResult r{0.0, 0.0, 0.0, 0.0,
                std::numeric_limits<double>::infinity(), 0, false};
    // Count usable bins in range (non-zero error) to decide if the fit
    // can succeed and ndf > 0.
    int n_used = 0;
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
        const double xc = h->GetXaxis()->GetBinCenter(b);
        if (xc < xlo || xc > xhi) continue;
        if (h->GetBinError(b) > 0.0) ++n_used;
    }
    if (n_used < 3) return r;   // need at least 3 points for ndf >= 1

    const std::string fname = std::string(h->GetName()) + "_flin";
    TF1* f = new TF1(fname.c_str(), "[0]+[1]*x", xlo, xhi);
    f->SetLineColor(kBlue + 2);
    f->SetLineWidth(2);
    const int status = h->Fit(f, "RQ0");
    if (status != 0) { delete f; return r; }
    r.a        = f->GetParameter(0);
    r.a_err    = f->GetParError(0);
    r.b        = f->GetParameter(1);
    r.b_err    = f->GetParError(1);
    r.ndf      = f->GetNDF();
    r.chi2_ndf = (r.ndf > 0) ? f->GetChisquare() / r.ndf : 0.0;
    r.ok       = (r.ndf > 0);
    delete f;
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

// Map a non-negative float quality value Q to the smallest grid index whose
// cut value passes (cut value = kCutMin + idx*kCutStep, and Q passes the cut
// iff Q <= cut). The smallest passing cut is ceil(max(0, Q)) rounded up to
// the next multiple of kCutStep. Returns -1 if Q exceeds the largest cut.
inline int effIdx(float Q) {
    const float Qpos = std::max(0.0f, Q);
    const int q_int  = static_cast<int>(std::ceil(Qpos));
    if (q_int <= kCutMin) return 0;
    if (q_int > kCutMax)  return -1;
    // ceil((q_int - kCutMin) / kCutStep)
    const int idx = (q_int - kCutMin + kCutStep - 1) / kCutStep;
    return (idx > kNCuts - 1) ? -1 : idx;
}

inline int cutValue(int idx) {
    return kCutMin + idx * kCutStep;
}

}  // anonymous namespace

void mass_ee_matchqualnorm_scan_cascade_sim() {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetErrorX(0.0);
    gStyle->SetPalette(kBird);

    std::cout << "Scanning ep_richmatchqualitynorm, em_richmatchqualitynorm on ["
              << kCutMin << ", " << kCutMax << "] step " << kCutStep
              << "  (n_cuts=" << kNCuts << ", " << kNCuts * kNCuts
              << " combinations / scenario)\n"
              << "Convention: qual <= cut passes (lower = tighter; cut="
              << kCutMax << " keeps essentially everything in the low region).\n"
              << "FOM: chi2/ndf of LINEAR fit y=a+b*x to N_PT2/N_PT3 in ["
              << kFitLo << ", " << kFitHi << "] GeV/c^2 (LOWER = FLATTER).\n";

    gSystem->mkdir("plots/output", true);

    const int n_scen = static_cast<int>(kScenarios.size());

    // ========================================================================
    // Allocate per-scenario per-(eff_ep_idx, eff_em_idx) per-trigger
    // "max-grid" histograms. They will hold the m_ee distribution of events
    // whose SMALLEST passing cut indices are (ei, mi). The cumulative
    // h_PT3 / h_PT2 at any cut point (ei0, mi0) is the sum over ei >= ei0
    // and mi >= mi0 (because for <= convention, an event with eff idx X
    // passes all cuts X..max).
    //
    // NB: total number of TH1D allocated = 6 * 21 * 21 * 2 = 5292. Each is
    // 70 bins; memory is well under 1 GB.
    // ========================================================================
    using HistGrid =
        std::vector<std::vector<std::vector<std::array<TH1D*, NUM_TRG>>>>;
    HistGrid h_max(n_scen,
        std::vector<std::vector<std::array<TH1D*, NUM_TRG>>>(kNCuts,
            std::vector<std::array<TH1D*, NUM_TRG>>(kNCuts)));
    for (int s = 0; s < n_scen; ++s) {
        for (int ei = 0; ei < kNCuts; ++ei) {
            for (int mi = 0; mi < kNCuts; ++mi) {
                for (int t = 0; t < NUM_TRG; ++t) {
                    TH1D* h = new TH1D(
                        Form("hmax_s%d_e%d_m%d_t%d", s, ei, mi, t),
                        "", kNb, kXmin, kXmax);
                    h->Sumw2();
                    h->SetDirectory(nullptr);
                    h_max[s][ei][mi][t] = h;
                }
            }
        }
    }

    // ========================================================================
    // Single pass over each input tree.
    //
    // step0 uses dilepton_nt_cor with m_ee_sim.
    // step1..step5 use dilepton_nt with m_ee and accumulate at successively
    // tighter scenario cuts. We do step1..step5 in ONE pass through
    // dilepton_nt, branching on the scenario gate per event.
    //
    // Pre-cut totals (for retention computation): total summed weight of
    // events that pass the scenario gate AND have either PT2 or PT3.
    // ========================================================================
    std::vector<double> total_w(n_scen, 0.0);

    // ---- step0: dilepton_nt_cor / m_ee_sim ----
    {
        TFile f("output_epem_sim.root", "READ");
        if (f.IsZombie()) {
            std::cerr << "Cannot open output_epem_sim.root\n"; return;
        }
        TTreeReader r("dilepton_nt_cor", &f);
        TTreeReaderValue<float> v_mee (r, "m_ee_sim");
        TTreeReaderValue<float> v_pt3 (r, "pt3");
        TTreeReaderValue<float> v_pt2 (r, "pt2");
        TTreeReaderValue<float> v_w   (r, "sim_genweight");
        TTreeReaderValue<float> v_ep  (r, "ep_richmatchqualitynorm");
        TTreeReaderValue<float> v_em  (r, "em_richmatchqualitynorm");

        Long64_t n_total = 0, n_kept = 0;
        std::cout << "\nReading dilepton_nt_cor (step0 raw truth)...\n";
        while (r.Next()) {
            ++n_total;
            const bool is_pt3 = (static_cast<int>(*v_pt3) == 1);
            const bool is_pt2 = (static_cast<int>(*v_pt2) == 1);
            if (!is_pt3 && !is_pt2) continue;
            const double w = *v_w;
            total_w[0] += w;
            const int ei = effIdx(*v_ep);
            const int mi = effIdx(*v_em);
            if (ei < 0 || mi < 0) continue;
            ++n_kept;
            const double m = *v_mee;
            if (is_pt3) h_max[0][ei][mi][TRG_PT3]->Fill(m, w);
            if (is_pt2) h_max[0][ei][mi][TRG_PT2]->Fill(m, w);
        }
        std::cout << "  scanned " << n_total << ", kept " << n_kept
                  << " in-grid (" << (100.0 * n_kept / std::max<Long64_t>(1, n_total))
                  << " %); total weight (PT2||PT3) = " << total_w[0] << "\n";
    }

    // ---- steps 1..5: dilepton_nt / m_ee, ONE pass with per-step gating ----
    {
        TFile f("output_epem_sim.root", "READ");
        if (f.IsZombie()) {
            std::cerr << "Cannot open output_epem_sim.root\n"; return;
        }
        TTreeReader r("dilepton_nt", &f);
        TTreeReaderValue<float> v_mee  (r, "m_ee");
        TTreeReaderValue<float> v_pt3  (r, "pt3");
        TTreeReaderValue<float> v_pt2  (r, "pt2");
        TTreeReaderValue<float> v_w    (r, "sim_genweight");
        TTreeReaderValue<float> v_ep   (r, "ep_richmatchqualitynorm");
        TTreeReaderValue<float> v_em   (r, "em_richmatchqualitynorm");
        TTreeReaderValue<float> v_vz   (r, "eVertReco_z");
        TTreeReaderValue<float> v_best (r, "isBest");
        TTreeReaderValue<float> v_epid (r, "ep_sim_id");
        TTreeReaderValue<float> v_emid (r, "em_sim_id");
        TTreeReaderValue<float> v_sv   (r, "epem_same_vertex");

        Long64_t n_total = 0;
        std::cout << "Reading dilepton_nt (steps 1..5)...\n";
        while (r.Next()) {
            ++n_total;
            const bool is_pt3 = (static_cast<int>(*v_pt3) == 1);
            const bool is_pt2 = (static_cast<int>(*v_pt2) == 1);
            if (!is_pt3 && !is_pt2) continue;

            // Cumulative scenario gates
            const bool g_raw     = true;
            const bool g_vertex  = g_raw    && (*v_vz   > -500.0f);
            const bool g_best    = g_vertex && (static_cast<int>(*v_best) == 1);
            const bool g_simid   = g_best   && (static_cast<int>(*v_epid) == 2)
                                            && (static_cast<int>(*v_emid) == 3);
            const bool g_samevtx = g_simid  && (static_cast<int>(*v_sv) == 1);

            const double w = *v_w;
            const double m = *v_mee;

            // Track scenario pre-cut weight totals (denominator for retention)
            if (g_raw)     total_w[1] += w;
            if (g_vertex)  total_w[2] += w;
            if (g_best)    total_w[3] += w;
            if (g_simid)   total_w[4] += w;
            if (g_samevtx) total_w[5] += w;

            const int ei = effIdx(*v_ep);
            const int mi = effIdx(*v_em);
            if (ei < 0 || mi < 0) continue;

            auto fillIfGate = [&](int sidx, bool gate) {
                if (!gate) return;
                if (is_pt3) h_max[sidx][ei][mi][TRG_PT3]->Fill(m, w);
                if (is_pt2) h_max[sidx][ei][mi][TRG_PT2]->Fill(m, w);
            };
            fillIfGate(1, g_raw);
            fillIfGate(2, g_vertex);
            fillIfGate(3, g_best);
            fillIfGate(4, g_simid);
            fillIfGate(5, g_samevtx);
        }
        std::cout << "  scanned " << n_total << " entries; scenario weights:";
        for (int s = 1; s < n_scen; ++s)
            std::cout << "  step" << s << "=" << total_w[s];
        std::cout << "\n";
    }

    // ========================================================================
    // Build cumulative h_PT3 / h_PT2 from h_max via inclusive forward-sum
    // on (ei, mi): cum[ei][mi] = sum_{ei' >= ei, mi' >= mi} h_max[ei'][mi'].
    //
    // We build the cumulative directly into the storage by doing a reverse-
    // direction sweep: for each fixed mi, sweep ei from high to low; then
    // for each fixed ei, sweep mi from high to low. (Two passes of a 2D
    // suffix-sum on histograms.)
    // ========================================================================
    // PREFIX sum: cumulative at cut idx (ec, mc) = sum over eff_ep_idx <= ec
    // AND eff_em_idx <= mc, because for <= convention an event whose smallest
    // passing cut index is X passes all cut indices in [X, kNCuts-1] (i.e.,
    // appears in every cumulative bin with index >= X).
    auto buildCumulative = [&](int sidx) {
        // pass 1: along ei (forward)
        for (int mi = 0; mi < kNCuts; ++mi) {
            for (int ei = 1; ei < kNCuts; ++ei) {
                for (int t = 0; t < NUM_TRG; ++t) {
                    h_max[sidx][ei][mi][t]->Add(h_max[sidx][ei - 1][mi][t]);
                }
            }
        }
        // pass 2: along mi (forward)
        for (int ei = 0; ei < kNCuts; ++ei) {
            for (int mi = 1; mi < kNCuts; ++mi) {
                for (int t = 0; t < NUM_TRG; ++t) {
                    h_max[sidx][ei][mi][t]->Add(h_max[sidx][ei][mi - 1][t]);
                }
            }
        }
    };
    for (int s = 0; s < n_scen; ++s) buildCumulative(s);

    // ========================================================================
    // Per-scenario heatmap of chi2/ndf, with best point identified.
    // ========================================================================
    std::vector<TH2D*> h_chi2(n_scen, nullptr);
    std::vector<BestPoint> best(n_scen);
    std::vector<FitResult> baseline(n_scen);   // (ei=kNCuts-1, mi=kNCuts-1) = no cut
    std::vector<double>    baseline_chi2(n_scen, std::numeric_limits<double>::quiet_NaN());

    for (int s = 0; s < n_scen; ++s) {
        h_chi2[s] = new TH2D(
            Form("h_chi2_%s", kScenarios[s].id.c_str()),
            Form("MID linear fit #chi^{2}/ndf - %s;"
                 "ep_richmatchqualitynorm #leq;"
                 "em_richmatchqualitynorm #leq;"
                 "#chi^{2}/ndf",
                 kScenarios[s].label.c_str()),
            kNCuts, cutValue(0) - 0.5 * kCutStep,
                    cutValue(kNCuts - 1) + 0.5 * kCutStep,
            kNCuts, cutValue(0) - 0.5 * kCutStep,
                    cutValue(kNCuts - 1) + 0.5 * kCutStep);
        h_chi2[s]->SetDirectory(nullptr);

        for (int ei = 0; ei < kNCuts; ++ei) {
            for (int mi = 0; mi < kNCuts; ++mi) {
                TH1D* h_PT3 = h_max[s][ei][mi][TRG_PT3];
                TH1D* h_PT2 = h_max[s][ei][mi][TRG_PT2];

                TH1D* r = makeRatio(h_PT2, h_PT3,
                    Form("r_tmp_s%d_e%d_m%d", s, ei, mi));
                const FitResult fr = fitLinRange(r, kFitLo, kFitHi);
                delete r;

                if (fr.ok)
                    h_chi2[s]->SetBinContent(ei + 1, mi + 1, fr.chi2_ndf);

                // Retention = (PT2+PT3 weighted yield kept) / total scenario weight
                const double kept = h_PT3->Integral() + h_PT2->Integral();
                const double ret  = (total_w[s] > 0.0)
                                  ? kept / total_w[s] : 0.0;

                if (fr.ok && fr.chi2_ndf < best[s].chi2_ndf) {
                    best[s].ep_cut    = cutValue(ei);
                    best[s].em_cut    = cutValue(mi);
                    best[s].a         = fr.a;
                    best[s].a_err     = fr.a_err;
                    best[s].b         = fr.b;
                    best[s].b_err     = fr.b_err;
                    best[s].chi2_ndf  = fr.chi2_ndf;
                    best[s].ndf       = fr.ndf;
                    best[s].retention = ret;
                }
            }
        }

        // Baseline = loosest cut (everything passes within scan)
        {
            const int ei = kNCuts - 1, mi = kNCuts - 1;
            TH1D* h_PT3 = h_max[s][ei][mi][TRG_PT3];
            TH1D* h_PT2 = h_max[s][ei][mi][TRG_PT2];
            TH1D* r = makeRatio(h_PT2, h_PT3,
                Form("r_base_s%d", s));
            baseline[s] = fitLinRange(r, kFitLo, kFitHi);
            delete r;
            baseline_chi2[s] = baseline[s].ok ? baseline[s].chi2_ndf
                                              : std::numeric_limits<double>::quiet_NaN();
        }
    }

    // ========================================================================
    // Summary table
    // ========================================================================
    auto fmtFloat = [](double v, int prec = 3) {
        std::ostringstream os;
        os << std::fixed << std::setprecision(prec) << v;
        return os.str();
    };

    std::cout << "\n=== SUMMARY: best (ep_cut, em_cut) per cascade scenario "
                 "(min chi2/ndf of MID linear fit) ===\n"
              << "  Fit range: [" << kFitLo << ", " << kFitHi << "] GeV/c^2\n"
              << "  Grid: cut in [" << kCutMin << ", " << kCutMax
              << "] step " << kCutStep << " (" << kNCuts << " points/leg)\n\n";

    std::cout << "  " << std::left << std::setw(50) << "scenario"
              << std::right << std::setw(7)  << "ep_cut"
              << std::setw(7)  << "em_cut"
              << std::setw(12) << "chi2/ndf"
              << std::setw(12) << "baseline"
              << std::setw(12) << "a"
              << std::setw(12) << "b"
              << std::setw(11) << "retention"
              << "   flags\n";
    std::cout << "  " << std::string(132, '-') << "\n";
    for (int s = 0; s < n_scen; ++s) {
        const auto& b = best[s];
        const bool hits_edge =
               (b.ep_cut == kCutMin) || (b.ep_cut == kCutMax)
            || (b.em_cut == kCutMin) || (b.em_cut == kCutMax);
        const bool low_retention = (b.retention < 0.01);
        std::string flags;
        if (hits_edge)     flags += " [EDGE-RED-FLAG]";
        if (low_retention) flags += " [LOW-RET-RED-FLAG]";
        std::cout << "  " << std::left << std::setw(50) << kScenarios[s].label
                  << std::right << std::setw(7)  << b.ep_cut
                  << std::setw(7)  << b.em_cut
                  << std::setw(12) << fmtFloat(b.chi2_ndf, 3)
                  << std::setw(12) << (std::isnan(baseline_chi2[s])
                                         ? std::string("NA")
                                         : fmtFloat(baseline_chi2[s], 3))
                  << std::setw(12) << fmtFloat(b.a, 4)
                  << std::setw(12) << fmtFloat(b.b, 4)
                  << std::setw(11) << (fmtFloat(100.0 * b.retention, 2) + "%")
                  << "   " << flags << "\n";
    }
    std::cout << "\n";

    // ========================================================================
    // OUTPUTS — per scenario heatmap + best-case 1x2 canvas
    // ========================================================================
    for (int s = 0; s < n_scen; ++s) {
        const auto& Sc = kScenarios[s];
        const auto& b  = best[s];

        // ---- Heatmap ----
        {
            TCanvas* c_heat = new TCanvas(
                Form("c_heat_%s", Sc.id.c_str()),
                Form("Matchqualnorm scan chi2/ndf landscape - %s",
                     Sc.label.c_str()),
                900, 750);
            gPad->SetMargin(0.13, 0.16, 0.13, 0.11);
            gPad->SetLogz(true);
            h_chi2[s]->Draw("COLZ");

            TMarker* mk = new TMarker(b.ep_cut, b.em_cut, 29);
            mk->SetMarkerSize(2.4);
            mk->SetMarkerColor(kRed + 1);
            mk->Draw();

            TLatex tx; tx.SetNDC(); tx.SetTextSize(0.030);
            tx.DrawLatex(0.14, 0.94,
                Form("best: (ep#leq%d, em#leq%d)  #chi^{2}/ndf=%.3f  ret=%.2f%%",
                     b.ep_cut, b.em_cut, b.chi2_ndf, 100.0 * b.retention));
            tx.SetTextSize(0.026);
            tx.SetTextColor(kGray + 3);
            tx.DrawLatex(0.14, 0.905,
                Form("baseline (ep#leq%d, em#leq%d) #chi^{2}/ndf=%.3f",
                     cutValue(kNCuts - 1), cutValue(kNCuts - 1),
                     std::isnan(baseline_chi2[s]) ? 0.0 : baseline_chi2[s]));

            const std::string base = Form(
                "plots/output/scan_matchqualnorm_cascade_sim_%s_heatmap",
                Sc.id.c_str());
            c_heat->SaveAs((base + ".pdf").c_str());
            c_heat->SaveAs((base + ".png").c_str());
            std::cout << "  saved: " << base << ".{pdf,png}\n";
        }

        // ---- Best-case 1x2 canvas ----
        {
            const int ei = (b.ep_cut - kCutMin) / kCutStep;
            const int mi = (b.em_cut - kCutMin) / kCutStep;
            TH1D* h_PT3 = static_cast<TH1D*>(h_max[s][ei][mi][TRG_PT3]->Clone(
                Form("c_PT3_%s", Sc.id.c_str())));
            TH1D* h_PT2 = static_cast<TH1D*>(h_max[s][ei][mi][TRG_PT2]->Clone(
                Form("c_PT2_%s", Sc.id.c_str())));
            h_PT3->SetDirectory(nullptr);
            h_PT2->SetDirectory(nullptr);
            TH1D* r_ratio = makeRatio(h_PT2, h_PT3,
                Form("c_ratio_%s", Sc.id.c_str()));

            // Re-fit so we can attach the function for drawing
            TF1* f_lin = new TF1(
                Form("f_lin_%s", Sc.id.c_str()),
                "[0]+[1]*x", kFitLo, kFitHi);
            f_lin->SetLineColor(kBlue + 2);
            f_lin->SetLineWidth(3);
            r_ratio->Fit(f_lin, "RQ0+");

            styleSpec(h_PT3, kBlack,   21);
            styleSpec(h_PT2, kRed + 1, 20);
            styleSpec(r_ratio, kBlue + 1, 20);

            TCanvas* c = new TCanvas(
                Form("c_best_%s", Sc.id.c_str()),
                Form("Best matchqualnorm - %s  ep<=%d  em<=%d",
                     Sc.label.c_str(), b.ep_cut, b.em_cut),
                1700, 600);
            c->Divide(2, 1, 0.001, 0.001);

            // -- Left: m_ee spectrum log Y --
            c->cd(1);
            gPad->SetLogy(true);
            gPad->SetMargin(0.13, 0.04, 0.13, 0.10);
            const double y_max = std::max(h_PT3->GetMaximum(),
                                          h_PT2->GetMaximum());
            const std::string mass_var_label =
                (Sc.sel == Scenario::SEL_RAW_TRUTH)
                    ? "M_{e^{+}e^{-}}^{sim}"
                    : "M_{e^{+}e^{-}}";
            h_PT3->SetTitle(Form(
                "%s at best cuts ep#leq%d em#leq%d;%s [GeV/c^{2}];weighted entries / 20 MeV",
                Sc.label.c_str(), b.ep_cut, b.em_cut,
                mass_var_label.c_str()));
            h_PT3->GetYaxis()->SetRangeUser(std::max(1e-2, y_max * 1e-6),
                                            y_max * 5.0);
            h_PT3->Draw("E1");
            h_PT2->Draw("E1 SAME");

            TLegend* leg = new TLegend(0.55, 0.72, 0.95, 0.88);
            leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.034);
            leg->AddEntry(h_PT3,
                Form("PT3  (#sum w = %.0f)", h_PT3->Integral()), "lpe");
            leg->AddEntry(h_PT2,
                Form("PT2  (#sum w = %.0f)", h_PT2->Integral()), "lpe");
            leg->Draw();

            // -- Right: ratio + linear fit, Y in [0.5, 2.0] --
            c->cd(2);
            gPad->SetMargin(0.13, 0.04, 0.13, 0.10);
            r_ratio->SetTitle(Form(
                "N_{PT2}/N_{PT3} - %s;%s [GeV/c^{2}];N_{PT2}/N_{PT3}",
                Sc.label.c_str(), mass_var_label.c_str()));
            r_ratio->GetYaxis()->SetRangeUser(0.5, 2.0);
            r_ratio->Draw("E1");

            // Shade MID fit window
            TBox* box = new TBox(kFitLo, 0.5, kFitHi, 2.0);
            box->SetFillColorAlpha(kAzure + 6, 0.15);
            box->SetLineColor(kAzure + 1); box->SetLineStyle(2);
            box->Draw("SAME");
            r_ratio->Draw("E1 SAME");
            f_lin->Draw("SAME");

            TLine* lref = new TLine(kXmin, 1.0, kXmax, 1.0);
            lref->SetLineStyle(3); lref->SetLineColor(kGray + 2); lref->Draw();

            TLatex tx; tx.SetNDC(); tx.SetTextAlign(11); tx.SetTextSize(0.032);
            tx.DrawLatex(0.16, 0.90,
                Form("best ep#leq%d  em#leq%d   retention = %.2f%%",
                     b.ep_cut, b.em_cut, 100.0 * b.retention));
            tx.SetTextColor(kBlue + 2);
            tx.SetTextSize(0.030);
            tx.DrawLatex(0.16, 0.85,
                Form("MID lin [%.2f, %.2f]: a=%.4f #pm %.4f  b=%.4f #pm %.4f",
                     kFitLo, kFitHi, b.a, b.a_err, b.b, b.b_err));
            tx.DrawLatex(0.16, 0.81,
                Form("#chi^{2}/ndf = %.3f  (ndf = %d)",
                     b.chi2_ndf, b.ndf));
            tx.SetTextColor(kGray + 3);
            tx.DrawLatex(0.16, 0.77,
                Form("baseline (loosest) #chi^{2}/ndf = %.3f",
                     std::isnan(baseline_chi2[s]) ? 0.0 : baseline_chi2[s]));

            // Red flag overlay
            const bool hits_edge =
                   (b.ep_cut == kCutMin) || (b.ep_cut == kCutMax)
                || (b.em_cut == kCutMin) || (b.em_cut == kCutMax);
            const bool low_retention = (b.retention < 0.01);
            if (hits_edge || low_retention) {
                tx.SetTextColor(kRed + 1);
                tx.SetTextSize(0.034);
                std::string msg = "RED FLAG:";
                if (hits_edge)     msg += " boundary";
                if (low_retention) msg += " low-ret";
                tx.DrawLatex(0.16, 0.72, msg.c_str());
            }

            const std::string base = Form(
                "plots/output/scan_matchqualnorm_cascade_sim_%s_best",
                Sc.id.c_str());
            c->SaveAs((base + ".pdf").c_str());
            c->SaveAs((base + ".png").c_str());
            std::cout << "  saved: " << base << ".{pdf,png}\n";
        }
    }

    std::cout << "\nDone. "
              << (n_scen * 2) << " canvases (" << (n_scen * 4)
              << " files) saved under plots/output/.\n";
}
