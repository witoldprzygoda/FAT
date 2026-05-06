// research/plots/compare_yields_pi0.C — exp vs sim comparison for π⁰ Dalitz,
// with sim normalized to exp in the high-OA tail.
//
// Reads:
//   research/fit_results_<fl>.root        (exp, this branch — default 'rec')
//   research/fit_results_<fl>_sim.root    (sim — copied from
//                                          pp45_epem_sim with _sim suffix)
//
// Both files come from plots/fit_pi0.C and contain five yield variants per
// OA slice: ±1σ, ±2σ, ±3σ, [μ−5σ, μ+3σ], and full-fit-range.
//
// Normalization strategy:
//   Single anchor scale, computed from the ±1σ yield in the high-OA tail:
//       S_exp_1σ = Σ_{oa_lo ≥ 9}  yield_exp_1σ(slice)
//       S_sim_1σ = Σ_{oa_lo ≥ 9}  yield_sim_1σ(slice)
//       scale    = S_exp_1σ / S_sim_1σ
//   This single `scale` is then applied uniformly to all sim variants —
//   by construction the ±1σ ratio in the high-OA region equals 1; other
//   variants (±2σ, ±3σ, asym, full) deviate from 1 by however much their
//   integrated yield differs in shape between exp and sim.
//
// Plot 1 — yield_pi0_compare_norm_1sig_<fl>.{pdf,png}
//   Raw exp(1σ) overlaid with sim(1σ) × scale.
//   Per-slice (X, Y_exp, Y_sim×scale, errors) printed to terminal.
//
// Plot 2 — yield_pi0_compare_ratio_norm_<fl>.{pdf,png}
//   Per-slice ratio exp_t / (scale · sim_t) for all five variants.
//   Legend shows the mean ratio in the high-OA window for each variant.
//
// Usage (from research/):
//   root -l -b -q plots/compare_yields_pi0.C            # exp REC vs sim REC
//   root -l -b -q 'plots/compare_yields_pi0.C("cor")'   # exp COR vs sim COR

#include <TFile.h>
#include <TTree.h>
#include <TAxis.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TString.h>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>

namespace {
    constexpr int kNTypes = 5;
    constexpr const char* kBrY[kNTypes] = {
        "yield_data_1sig",     "yield_data_2sig",     "yield_data_3sig",
        "yield_data_m5p3sig",  "yield_data_full"
    };
    constexpr const char* kBrE[kNTypes] = {
        "yield_data_1sig_err",     "yield_data_2sig_err",     "yield_data_3sig_err",
        "yield_data_m5p3sig_err",  "yield_data_full_err"
    };
    constexpr const char* kLabel[kNTypes] = {
        "#mu #pm 1#sigma",   "#mu #pm 2#sigma",   "#mu #pm 3#sigma",
        "[#mu#minus5#sigma, #mu+3#sigma]",   "full fit range"
    };

    // High-OA window [oa_lo ≥ kScaleOAmin] used for the sim → exp normalization.
    constexpr double kScaleOAmin = 9.0;

    struct SliceRow {
        double oa_c, oa_w;
        double oa_lo, oa_hi;
        double y[kNTypes], e[kNTypes];
    };

    bool readTree(const std::string& path, std::vector<SliceRow>& rows) {
        TFile* f = TFile::Open(path.c_str(), "READ");
        if (!f || f->IsZombie()) {
            std::cerr << "Cannot open " << path << "\n"; return false;
        }
        auto* t = (TTree*)f->Get("fit_results");
        if (!t) { std::cerr << "No TTree 'fit_results' in " << path << "\n"; return false; }

        int   panel_idx;
        float oa_lo, oa_hi;
        float y[kNTypes], e[kNTypes];
        t->SetBranchAddress("panel_idx", &panel_idx);
        t->SetBranchAddress("oa_lo",     &oa_lo);
        t->SetBranchAddress("oa_hi",     &oa_hi);
        for (int i = 0; i < kNTypes; ++i) {
            t->SetBranchAddress(kBrY[i], &y[i]);
            t->SetBranchAddress(kBrE[i], &e[i]);
        }

        rows.clear();
        for (Long64_t ev = 0; ev < t->GetEntries(); ++ev) {
            t->GetEntry(ev);
            if (panel_idx == 0) continue;
            SliceRow r;
            r.oa_lo = oa_lo;  r.oa_hi = oa_hi;
            r.oa_c  = 0.5 * (oa_lo + oa_hi);
            r.oa_w  = 0.5 * (oa_hi - oa_lo);
            for (int i = 0; i < kNTypes; ++i) { r.y[i] = y[i];  r.e[i] = e[i]; }
            rows.push_back(r);
        }
        f->Close();
        // Sort by oa_lo so exp/sim arrays line up even if the producing
        // macro wrote rows out of order (e.g. exp has a reverse pass that
        // appends slices 1.6–2.4 at the end of the TTree).
        std::sort(rows.begin(), rows.end(),
                  [](const SliceRow& a, const SliceRow& b) {
                      return a.oa_lo < b.oa_lo;
                  });
        return true;
    }

    // Σ y_t over slices with oa_lo ≥ oa_min, for one yield type.
    double sumHighOA(const std::vector<SliceRow>& v, int type, double oa_min) {
        double s = 0.0;
        for (const auto& r : v) {
            if (r.oa_lo >= oa_min - 1e-6) s += r.y[type];
        }
        return s;
    }
}

void compare_yields_pi0(const char* exp_flavour = "rec") {

    std::string fl = exp_flavour;
    for (auto& c : fl) c = std::tolower(c);
    if (fl != "rec" && fl != "cor" && fl != "tru") {
        std::cerr << "Unknown flavour '" << exp_flavour << "'\n"; return;
    }

    const std::string path_exp = "fit_results_" + fl + ".root";
    const std::string path_sim = "fit_results_" + fl + "_sim.root";

    std::vector<SliceRow> exp_rows, sim_rows;
    if (!readTree(path_exp, exp_rows)) return;
    if (!readTree(path_sim, sim_rows)) return;

    if (exp_rows.size() != sim_rows.size()) {
        std::cerr << "Slice count mismatch: exp=" << exp_rows.size()
                  << " sim=" << sim_rows.size() << "\n";
        return;
    }
    const int N = (int)exp_rows.size();

    for (int i = 0; i < N; ++i) {
        if (std::fabs(exp_rows[i].oa_c - sim_rows[i].oa_c) > 1e-3) {
            std::cerr << "OA bin mismatch at i=" << i
                      << ": exp=" << exp_rows[i].oa_c
                      << " sim=" << sim_rows[i].oa_c << "\n";
            return;
        }
    }

    // Sums over high-OA window per yield type (used for legend "ratio").
    double S_exp[kNTypes], S_sim[kNTypes];
    for (int t = 0; t < kNTypes; ++t) {
        S_exp[t] = sumHighOA(exp_rows, t, kScaleOAmin);
        S_sim[t] = sumHighOA(sim_rows, t, kScaleOAmin);
    }

    // Single anchor scale: from ±1σ yields in OA ≥ kScaleOAmin window.
    if (S_sim[0] <= 0 || S_exp[0] <= 0) {
        std::cerr << "Empty 1σ high-OA totals (S_exp=" << S_exp[0]
                  << " S_sim=" << S_sim[0] << ") — cannot normalize.\n";
        return;
    }
    const double scale = S_exp[0] / S_sim[0];

    std::cout << "Anchor scale (1σ, OA #geq " << kScaleOAmin
              << " deg): S_exp = " << S_exp[0]
              << "   S_sim = " << S_sim[0]
              << "   scale = " << scale << "\n";
    std::cout << "High-OA mean ratio per variant (with single scale = "
              << scale << "):\n";
    for (int t = 0; t < kNTypes; ++t) {
        const double r_wide = S_exp[t] / (scale * S_sim[t]);
        std::cout << "  " << kBrY[t] << ":   S_exp = " << S_exp[t]
                  << "   S_sim = " << S_sim[t]
                  << "   ratio_high = " << r_wide << "\n";
    }

    gStyle->SetOptStat(0);
    gSystem->mkdir("plots/output", kTRUE);

    // X coordinates shared across plots.
    std::vector<double> X(N), EX(N);
    for (int i = 0; i < N; ++i) { X[i] = exp_rows[i].oa_c; EX[i] = exp_rows[i].oa_w; }

    // ------------------------------------------------------------------------
    // Plot 1 — raw exp(1σ) + sim(1σ) scaled to exp at oa_lo ≥ 9°
    // ------------------------------------------------------------------------
    {
        std::vector<double> Ye(N), EYe(N), Ys(N), EYs(N);
        const int t = 0;   // 1σ
        for (int i = 0; i < N; ++i) {
            Ye [i] = exp_rows[i].y[t];
            EYe[i] = exp_rows[i].e[t];
            Ys [i] = sim_rows[i].y[t] * scale;
            EYs[i] = sim_rows[i].e[t] * scale;
        }

        // --- Per-slice point dump for the 1σ overlay ---------------------
        std::cout << "\nPlot 1 (#mu#pm1#sigma) point dump — sim already scaled by "
                  << scale << ":\n";
        std::cout << "  "
                  << std::setw(11) << "OA_center"
                  << std::setw(15) << "Y_exp"
                  << std::setw(15) << "EY_exp"
                  << std::setw(15) << "Y_sim*scale"
                  << std::setw(15) << "EY_sim*scale" << "\n";
        std::cout << "  " << std::string(11+15*4, '-') << "\n";
        for (int i = 0; i < N; ++i) {
            std::cout << "  "
                      << std::setw(11) << std::fixed << std::setprecision(2) << X[i]
                      << std::setw(15) << std::fixed << std::setprecision(3) << Ye[i]
                      << std::setw(15) << std::fixed << std::setprecision(3) << EYe[i]
                      << std::setw(15) << std::fixed << std::setprecision(3) << Ys[i]
                      << std::setw(15) << std::fixed << std::setprecision(3) << EYs[i]
                      << "\n";
        }
        std::cout.unsetf(std::ios::fixed);
        std::cout << "\n";

        auto* g_exp = new TGraphErrors(N, X.data(), Ye.data(), EX.data(), EYe.data());
        auto* g_sim = new TGraphErrors(N, X.data(), Ys.data(), EX.data(), EYs.data());

        g_exp->SetMarkerColor(kBlack); g_exp->SetLineColor(kBlack);
        g_exp->SetMarkerStyle(20);     g_exp->SetMarkerSize(0.9);
        g_sim->SetMarkerColor(kRed);   g_sim->SetLineColor(kRed);
        g_sim->SetMarkerStyle(24);     g_sim->SetMarkerSize(0.9);

        auto* mg = new TMultiGraph();
        mg->Add(g_exp, "P");
        mg->Add(g_sim, "P");
        mg->SetTitle(TString::Format(
            "#pi^{0} yield (%s, #pm 1#sigma): exp vs sim "
            "[sim normalized to exp at OA #geq %.0f#circ];"
            "OA(e^{+}e^{-}) [deg];yield(#mu #pm 1#sigma)",
            fl.c_str(), kScaleOAmin));

        auto* c1 = new TCanvas("c_pi0_norm_1sig", "exp vs scaled sim", 1000, 700);
        c1->SetMargin(0.13, 0.05, 0.12, 0.08);
        c1->SetGrid();
        mg->Draw("A");
        mg->GetXaxis()->SetLimits(0.0, 15.0);

        auto* leg = new TLegend(0.55, 0.78, 0.94, 0.90);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.035);
        leg->AddEntry(g_exp, TString::Format("exp (%s)", fl.c_str()).Data(), "lpe");
        leg->AddEntry(g_sim,
            TString::Format("sim #times %.4g  (norm. at OA #geq %.0f#circ)",
                            scale, kScaleOAmin).Data(), "lpe");
        leg->Draw();

        c1->Update();
        const std::string out1 = "plots/output/yield_pi0_compare_norm_1sig_" + fl;
        c1->SaveAs((out1 + ".pdf").c_str());
        c1->SaveAs((out1 + ".png").c_str());
    }

    // ------------------------------------------------------------------------
    // Plot 2 — per-slice ratio across the full OA range, all five variants,
    //   single anchor scale (1σ-based at OA ≥ kScaleOAmin) used for every
    //   variant. By construction the average exp/sim ratio over OA ≥ 9°
    //   equals 1 for the 1σ series; deviations elsewhere show shape change.
    // ------------------------------------------------------------------------
    {
        const int colors[kNTypes] = {kBlue, kGreen + 2, kRed, kMagenta + 1, kBlack};
        const int styles[kNTypes] = {20, 21, 22, 29, 33};
        const double sizes[kNTypes] = {0.9, 0.9, 1.0, 1.1, 1.2};

        std::vector<TGraphErrors*> grs;
        grs.reserve(kNTypes);

        for (int t = 0; t < kNTypes; ++t) {
            std::vector<double> R(N), ER(N);
            for (int i = 0; i < N; ++i) {
                const double y_e = exp_rows[i].y[t];
                const double e_e = exp_rows[i].e[t];
                const double y_s = sim_rows[i].y[t];
                const double e_s = sim_rows[i].e[t];
                if (y_e <= 0 || y_s <= 0) {
                    R[i] = 0.0;  ER[i] = 0.0;  continue;
                }
                const double r   = y_e / (scale * y_s);
                const double rel = std::sqrt((e_e/y_e)*(e_e/y_e) + (e_s/y_s)*(e_s/y_s));
                R [i] = r;
                ER[i] = r * rel;
            }

            auto* g = new TGraphErrors(N, X.data(), R.data(), EX.data(), ER.data());
            g->SetMarkerColor(colors[t]); g->SetLineColor(colors[t]);
            g->SetMarkerStyle(styles[t]); g->SetMarkerSize(sizes[t]);
            grs.push_back(g);
        }

        auto* mg = new TMultiGraph();
        for (auto* g : grs) mg->Add(g, "P");
        mg->SetTitle(TString::Format(
            "#pi^{0} yield ratio: exp(%s) / (scale #times sim),   "
            "scale = %.4g  (anchored on 1#sigma at OA #geq %.0f#circ);"
            "OA(e^{+}e^{-}) [deg];"
            "yield_{exp} / (scale #times yield_{sim})",
            fl.c_str(), scale, kScaleOAmin));

        auto* c2 = new TCanvas("c_pi0_ratio", "exp/sim π⁰ ratio", 1100, 700);
        c2->SetMargin(0.13, 0.05, 0.12, 0.08);
        c2->SetGrid();
        mg->Draw("A");
        mg->GetXaxis()->SetLimits(0.0, 15.0);

        auto* lref = new TLine(0.0, 1.0, 15.0, 1.0);
        lref->SetLineColor(kGray + 2);
        lref->SetLineStyle(2);
        lref->SetLineWidth(2);
        lref->Draw();

        auto* leg = new TLegend(0.45, 0.62, 0.94, 0.90);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.030);
        for (int t = 0; t < kNTypes; ++t) {
            const double r_wide = (S_sim[t] > 0)
                                   ? S_exp[t] / (scale * S_sim[t]) : 0.0;
            leg->AddEntry(grs[t],
                TString::Format("%s  (#LTratio#GT_{OA #geq 9}=%.3f)", kLabel[t], r_wide).Data(),
                "lpe");
        }
        leg->AddEntry(lref, "match (= 1)", "l");
        leg->Draw();

        c2->Update();
        const std::string out2 = "plots/output/yield_pi0_compare_ratio_norm_" + fl;
        c2->SaveAs((out2 + ".pdf").c_str());
        c2->SaveAs((out2 + ".png").c_str());
    }

    std::cout << "Wrote:\n"
              << "  plots/output/yield_pi0_compare_norm_1sig_" << fl << ".{pdf,png}\n"
              << "  plots/output/yield_pi0_compare_ratio_norm_" << fl << ".{pdf,png}\n"
              << "  (" << N << " OA slices)\n";
}
