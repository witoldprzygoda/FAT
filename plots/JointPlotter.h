/**
 * @file JointPlotter.h
 * @brief Helpers for joint EXP + SIM plotting on the comparison branch.
 *
 * Each canvas shows:
 *   exp (3-curve all/CB/sig from PlotUtils CB extraction)
 *   sim (single blue connecting line, rescaled to match exp signal in a
 *        user-chosen control region — typically right tail).
 *
 * Required input files in CWD (renamed per-branch to avoid collision):
 *   output_pippimepem_exp.root      (exp signal, like-sign-pair "all")
 *   output_pippimepep_exp.root      (exp CB++)
 *   output_pippimemem_exp.root      (exp CB--)
 *   output_pippimepem_sim.root      (sim, weighted by sim_genweight)
 *
 * Usage in joint_*.C macros:
 *   PlotUtils  exp("output_pippimepem_exp.root",
 *                  "output_pippimepep_exp.root",
 *                  "output_pippimemem_exp.root");
 *   JointPlotter::SimSource sim("output_pippimepem_sim.root");
 *
 *   auto [a, c, s] = exp.drawSignal("pippimepem_nt", "m_ee", 200, 0, 1, "");
 *   auto* h_sim    = sim.draw("pippimepem_nt", "m_ee", 200, 0, 1, "");
 *
 *   double scale = JointPlotter::rescaleSimToData(h_sim, s, 0.5, 1.0);
 *   auto* cv = JointPlotter::drawJoint(a, c, s, h_sim,
 *                                      "M_{e^{+}e^{-}}", "c_m_ee_joint",
 *                                      true, scale);
 *   JointPlotter::save(cv, "m_ee");
 *
 * @author Witold Przygoda (witold.przygoda@uj.edu.pl)
 * @date 2025
 */

#ifndef JOINTPLOTTER_H
#define JOINTPLOTTER_H

#include "PlotUtils.h"
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TColor.h>
#include <TSystem.h>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

namespace JointPlotter {

// ============================================================================
// SimSource — wraps the sim ROOT file and produces TH1D from a sim ntuple
// with proper sim_genweight cut. Mirrors PlotUtils::drawNtupleSingle but
// without the styling (we apply our own line style via styleSimLine).
// ============================================================================
class SimSource {
public:
    explicit SimSource(const std::string& fname) {
        f_ = TFile::Open(fname.c_str(), "READ");
        if (!f_ || f_->IsZombie()) {
            std::cerr << "JointPlotter::SimSource: cannot open " << fname << "\n";
        }
        counter_ = 0;
    }
    ~SimSource() { if (f_) f_->Close(); }

    /// Draw `var` from `ntuple_name` in the sim file. `filter` is an optional
    /// boolean expression; sim_genweight is automatically applied.
    TH1D* draw(const std::string& ntuple_name,
               const std::string& var,
               int nbins, double xmin, double xmax,
               const std::string& filter = "")
    {
        if (!f_) return nullptr;
        auto* t = dynamic_cast<TTree*>(f_->Get(ntuple_name.c_str()));
        if (!t) {
            std::cerr << "JointPlotter::SimSource: tree '" << ntuple_name
                      << "' not in " << f_->GetName() << "\n";
            return nullptr;
        }
        std::string hname = "h_sim_" + var + "_" + std::to_string(counter_++);
        // strip non-name characters
        for (auto& ch : hname) if (ch == '/' || ch == ' ' || ch == '(' || ch == ')') ch = '_';

        auto* h = new TH1D(hname.c_str(), "", nbins, xmin, xmax);
        h->Sumw2();
        std::string drawcmd = var + ">>" + hname;
        std::string cut     = filter.empty() ? std::string("sim_genweight")
                                             : "(" + filter + ")*sim_genweight";
        t->Draw(drawcmd.c_str(), cut.c_str(), "goff");
        return h;
    }

private:
    TFile* f_ = nullptr;
    int counter_;
};

// ============================================================================
// Sim styling — drawn as a brown histogram step line on top of a light-gray
// stat-error band ((sim+err) and (sim-err) envelope, also stair-stepped).
//
//   styleSimLine(h)      — brown step line, no fill, no markers (drawn HIST)
//   styleSimBand(h)      — light-gray solid fill, no line, no markers
//                          (drawn E2 — per-bin error rectangles)
//
// The two stylers are applied to the same data (the band is a clone of the
// scaled sim with Sumw2 errors). Order in drawJoint is: band → line, so the
// brown step appears as the central value of the gray band.
// ============================================================================
namespace {
    inline int simBrown() {
        static int c = TColor::GetColor("#8B4513");   // SaddleBrown
        return c;
    }
}

inline void styleSimLine(TH1* h) {
    if (!h) return;
    h->SetMarkerStyle(0);
    h->SetMarkerSize(0);
    h->SetLineColor(simBrown());
    h->SetLineWidth(2);
    h->SetFillStyle(0);
}

inline void styleSimBand(TH1* h) {
    if (!h) return;
    h->SetMarkerStyle(0);
    h->SetMarkerSize(0);
    h->SetFillColor(kGray);
    h->SetFillStyle(1001);   // solid fill
    h->SetLineColor(kGray);
    h->SetLineWidth(0);
}

// ============================================================================
// findBestScale — pick the largest sim scale α such that α·sim stays under
// the EXP SIGNAL in MOST of the histogram, with two physically motivated
// guards.
//
//   1) "Yield window" — only consider bins inside the central 80% of the
//      data integral (default yield_threshold = 0.10 on each side). Tails
//      with vanishing statistics are excluded — they would otherwise drive
//      the constraint into noise.
//
//   2) "Overshoot budget" — within the yield window, allow α·sim to exceed
//      data in at most `overshoot_frac` (default 0.20) of the bins. The
//      largest α satisfying this is found by sorting per-bin breaking
//      ratios r[i] = data[i] / sim[i] and picking the m-th smallest where
//      m = floor(overshoot_frac * N_active). At α = r[m] exactly m bins
//      are at the breaking threshold and the rest of α·sim is ≤ data.
//
// Returns the chosen α (1.0 on degenerate input). NOT applied — caller
// must Scale() the histogram.
// ============================================================================
inline double findBestScale(TH1* h_sim, TH1* h_data,
                            double overshoot_frac  = 0.20,
                            double yield_threshold = 0.10)
{
    if (!h_sim || !h_data) return 1.0;
    const int nbins = h_data->GetNbinsX();
    const double total = h_data->Integral();
    if (total <= 0) return 1.0;

    // 1) Active region: bins where the running cumulative yield is in
    //    [yield_threshold, 1 - yield_threshold] of the total integral.
    int b_lo = -1, b_hi = -1;
    double cum = 0.0;
    for (int i = 1; i <= nbins; ++i) {
        cum += h_data->GetBinContent(i);
        const double frac = cum / total;
        if (b_lo < 0 && frac >= yield_threshold)        b_lo = i;
        if (b_hi < 0 && frac >= 1.0 - yield_threshold)  b_hi = i;
    }
    if (b_lo < 0) b_lo = 1;
    if (b_hi < 0) b_hi = nbins;

    // 2) Per-bin breaking ratios within the active region.
    std::vector<double> ratios;
    ratios.reserve(b_hi - b_lo + 1);
    for (int i = b_lo; i <= b_hi; ++i) {
        const double d = h_data->GetBinContent(i);
        const double s = h_sim ->GetBinContent(i);
        if (s > 0 && d > 0) ratios.push_back(d / s);
    }
    if (ratios.empty()) return 1.0;

    std::sort(ratios.begin(), ratios.end());

    // 3) Pick the m-th smallest breaking ratio where m = overshoot budget.
    //    α = ratios[m] yields exactly m overshoots (bins with α > r[i] for
    //    i = 0..m-1) and the remaining bins satisfy α·sim <= data.
    const int max_overshoots = static_cast<int>(overshoot_frac * ratios.size());
    int idx = max_overshoots;
    if (idx >= (int)ratios.size()) idx = ratios.size() - 1;
    if (idx < 0) idx = 0;
    const double scale = ratios[idx];

    std::cout << "  findBestScale: active bins=" << ratios.size()
              << "  max_overshoots=" << max_overshoots
              << "  α = " << scale
              << "  (yield window [" << b_lo << ", " << b_hi << "] of " << nbins << ")\n";
    return scale;
}

// Convenience: compute scale via findBestScale and apply it to h_sim.
inline double rescaleSimToData(TH1* h_sim, TH1* h_data,
                               double overshoot_frac  = 0.20,
                               double yield_threshold = 0.10)
{
    const double s = findBestScale(h_sim, h_data, overshoot_frac, yield_threshold);
    if (h_sim) h_sim->Scale(s);
    return s;
}

// ============================================================================
// rescaleSimInWindow — explicit integral matching:
//   α = ∫[lo,hi] data / ∫[lo,hi] sim
// Use when the analysis specifies a normalisation control region (e.g. a
// recoil-mass tail outside the resonance peak). Returns and applies α.
// ============================================================================
inline double rescaleSimInWindow(TH1* h_sim, TH1* h_data, double lo, double hi) {
    if (!h_sim || !h_data) return 1.0;
    const int blo = h_data->FindBin(lo);
    const int bhi = h_data->FindBin(hi);
    const double y_data = h_data->Integral(blo, bhi);
    const double y_sim  = h_sim ->Integral(blo, bhi);
    const double scale  = (y_sim > 0) ? y_data / y_sim : 1.0;
    h_sim->Scale(scale);
    std::cout << "  rescaleSimInWindow: norm [" << lo << ", " << hi
              << "]  ∫data=" << y_data << "  ∫sim=" << y_sim
              << "  α = " << scale << "\n";
    return scale;
}

// ============================================================================
// Draw exp triple + sim (band + step line) on one canvas.
//
// Drawing order (back → front), so h_sig stays on top and is not covered:
//   1) h_all  — black markers + error bars (sets axes)
//   2) h_cb   — red markers + error bars
//   3) h_sim_band — light-gray solid stat-error band (clone of h_sim, "E2")
//   4) h_sim  — brown step histogram line ("HIST")
//   5) h_sig  — blue markers + error bars (drawn LAST, on top of sim)
//
// Y-axis cap: 3x for log, 1.2x for linear (max of h_all and h_sim).
// ============================================================================
inline TCanvas* drawJoint(TH1* h_all, TH1* h_cb, TH1* h_sig, TH1* h_sim,
                          const std::string& canvasTitle,
                          const std::string& canvasName,
                          bool logy = false,
                          double sim_scale = 1.0,
                          double display_lo = 0.0, double display_hi = 0.0)
{
    TCanvas* c = new TCanvas(canvasName.c_str(), canvasTitle.c_str(), 800, 600);
    c->SetMargin(0.12, 0.05, 0.12, 0.08);
    if (logy) c->SetLogy();

    if (!canvasTitle.empty()) {
        h_all->SetTitle(canvasTitle.c_str());
    }

    // Optional X-axis zoom (display only — binning unchanged).
    const bool zoom = (display_lo < display_hi);
    int b_lo = 1, b_hi = h_all->GetNbinsX();
    if (zoom) {
        h_all->GetXaxis()->SetRangeUser(display_lo, display_hi);
        b_lo = h_all->FindBin(display_lo);
        b_hi = h_all->FindBin(display_hi);
    }

    // Y-axis range: scan only bins inside the visible window. Include sim
    // upper band edge (sim + sim_err) so it doesn't get clipped.
    double y_max = 0.0;
    for (int i = b_lo; i <= b_hi; ++i) {
        y_max = std::max(y_max, h_all->GetBinContent(i));
        if (h_sim) y_max = std::max(y_max, h_sim->GetBinContent(i) + h_sim->GetBinError(i));
    }
    h_all->SetMaximum(y_max * (logy ? 3.0 : 1.2));
    h_all->SetMinimum(logy ? 0.5 : 0.0);

    // 1, 2: data all + CB
    h_all->Draw("E");
    h_cb ->Draw("E SAME");

    // 3, 4: sim stat-error band + brown step line (under data signal)
    TH1D* h_sim_band = nullptr;
    if (h_sim) {
        h_sim_band = (TH1D*)h_sim->Clone((std::string(h_sim->GetName()) + "_band").c_str());
        styleSimBand(h_sim_band);
        h_sim_band->Draw("E2 SAME");

        styleSimLine(h_sim);   // brown step
        h_sim->Draw("HIST SAME");
    }

    // 5: signal LAST so it isn't covered by sim band/step
    h_sig->Draw("E SAME");

    TLegend* leg = new TLegend(0.55, 0.66, 0.94, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.032);
    leg->AddEntry(h_all, "exp: all (#pi^{+}#pi^{-}e^{+}e^{-})", "lpe");
    leg->AddEntry(h_cb,  "exp: CB (2#sqrt{N_{++}N_{--}})",      "lpe");
    leg->AddEntry(h_sig, "exp: signal (all - CB)",              "lpe");
    if (h_sim) {
        char sim_lbl[128];
        snprintf(sim_lbl, sizeof(sim_lbl), "sim #times %.3g", sim_scale);
        leg->AddEntry(h_sim,      sim_lbl,             "l");
        leg->AddEntry(h_sim_band, "sim stat. error",   "f");
    }
    leg->Draw();

    c->Update();
    return c;
}

// ============================================================================
// Save canvas as plots/output/joint_{basename}.{pdf,png}
// ============================================================================
inline void save(TCanvas* c, const std::string& basename) {
    gSystem->mkdir("plots/output", kTRUE);
    std::string pdf = "plots/output/joint_" + basename + ".pdf";
    std::string png = "plots/output/joint_" + basename + ".png";
    c->SaveAs(pdf.c_str());
    c->SaveAs(png.c_str());
    std::cout << "Saved: " << pdf << "\n";
}

} // namespace JointPlotter

#endif // JOINTPLOTTER_H
