// research/plots/fit_pi0.C — Crystal Ball + bg fit on the SIGNAL
// histogram (all − CB) per OA slice (EXP DATA).
//
// Reads three files produced by `meson_research`:
//   research_epem.root   (opposite-sign all)
//   research_epep.root   (++ pairs)
//   research_emem.root   (-- pairs)
// Builds bin-by-bin
//   CB    = 2 √(N_++ N_--)         (with quadratic error propagation)
//   sig   = all − CB                (TH1::Add → errors in quadrature)
// and fits sig with the same model as on the SMASH simulation:
//
//   F(m) = N · CrystalBall_left(m; μ, σ, α, n)
//        + a₀ + a₁·m + a₂·m²
//        + b · exp(-c·m)               (+ a₃·m³ at very high OA)
//
// Procedure:
//   (1) Pre-fit Gaussian on the apex window [0.125, 0.145] — seeds μ.
//   (2) Single full fit on [kFitMin, kFitMax] (OA-dependent kFitMin).
// CB tail constraints + α tightening at high OA mirror the SIM-side macro.
//
// Per-slice canvases + a separate residual canvas (sig − bg) with the CB
// overlay and the true ±{1,2,3}σ data yields. fit_results.root TTree.
//
// Usage (from research/):
//   root -l -b -q plots/fit_pi0.C            # REC (default)
//   root -l -b -q 'plots/fit_pi0.C("cor")'   # COR

#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <utility>
#include <vector>

namespace {
    constexpr double kSliceMin  = 0.0;
    constexpr double kSliceMax  = 15.0;
    constexpr double kSliceStep = 0.2;

    // Fit range: [0.04, 0.45] for low/mid OA. At very high OA the data
    // is cut off below ~0.05–0.06 and forcing the fit to start at 0.04
    // distorts the bg shape, so we push the lower edge upward there.
    constexpr double kFitMin = 0.04;
    constexpr double kFitMax = 0.45;

    constexpr double kOA_FitMin005 = 11.0;   // OA ≥ this: fitMin = 0.05
    constexpr double kOA_FitMin006 = 13.8;   // OA ≥ this: fitMin = 0.06
    // α tightened across all OA — without this, low-OA fits found a
    // degenerate minimum (small α + wider σ + μ pinned at 0.140 upper
    // limit) where the long CB tail compensates by pulling μ up.
    constexpr double kOA_AlphaTighten  = 0.0;    // always tighten
    constexpr double kAlphaTightMin    = 1.5;
    constexpr double kAlphaTightMaxMid = 1.65;   // OA < 9.0
    constexpr double kAlphaTightMaxHi  = 2.0;    // OA ≥ 9.0
    // P3 enabled globally — a₃·m³ gives the bg polynomial enough freedom to
    // form a "broad hump" through [0.03, 0.20] (rising slope on the left,
    // peak ~0.05–0.10, falling into the right sideband). Without it the
    // P2 + exp bg behaves nearly linearly through [0.02, 0.08] and cuts
    // straight across the data points there.
    constexpr double kOA_EnableP3      = 0.0;    // always on

    // Gauss-bg term enabled only at low OA where the bg has a visible
    // convex hump in the low-m sideband. From OA ≥ this threshold the
    // simpler CB + P3 + exp shape (concave) describes the bg better,
    // and adding a Gauss-bg degrades the fit visually.
    constexpr double kOA_DisableGaussBg = 6.0;

    // Pre-fit Gaussian apex window (seeds μ).
    constexpr double kPreMin = 0.125;
    constexpr double kPreMax = 0.145;

    constexpr double kPi0Mass = 0.135;

    std::string fmtEdge(double x) {
        char buf[32];
        snprintf(buf, sizeof(buf), "%.1f", x);
        std::string s(buf);
        for (auto& c : s) if (c == '.') c = 'p';
        return s;
    }
}

// ============================================================================
// Crystal Ball — single-sided, tail on the LEFT (m < μ when α > 0).
// ============================================================================
double CrystalBallLeft(double *x, double *par) {
    const double m     = x[0];
    const double N     = par[0];
    const double mu    = par[1];
    const double sigma = par[2];
    const double alpha = par[3];
    const double n     = par[4];

    if (sigma <= 0.0 || alpha <= 0.0 || n <= 0.0) return 0.0;

    const double t = (m - mu) / sigma;
    const double absAlpha = std::fabs(alpha);
    if (t > -absAlpha) return N * std::exp(-0.5 * t * t);
    const double A = std::pow(n / absAlpha, n) * std::exp(-0.5 * absAlpha * absAlpha);
    const double B = n / absAlpha - absAlpha;
    return N * A * std::pow(B - t, -n);
}

// ============================================================================
// Full model: CB + P3 + exp + broad Gauss-bg.
//   par[0..4]  = CB (N, μ, σ, α, n)
//   par[5..7]  = P2 (a0, a1, a2)
//   par[8..9]  = exp (b, c)  →  b·exp(-c·m)
//   par[10]    = a3 cubic coefficient
//   par[11..13]= broad Gauss bg (N_bg_g, m_bg_g, σ_bg_g) — captures the
//                "broad hump" in the [0.03, 0.20] sideband (conversions
//                + bremsstrahlung pairs); too broad to mimic the π⁰ peak.
// ============================================================================
double FitFunction(double *x, double *par) {
    const double m  = x[0];
    const double cb = CrystalBallLeft(x, par);
    const double p3 = par[5] + par[6]*m + par[7]*m*m + par[10]*m*m*m;
    const double ex = par[8] * std::exp(-par[9]*m);
    const double sigG = par[13];
    const double bg_g = (sigG > 0)
        ? par[11] * std::exp(-0.5 * ((m - par[12]) / sigG) * ((m - par[12]) / sigG))
        : 0.0;
    return cb + p3 + ex + bg_g;
}

// ============================================================================
// Pre-fit Gaussian on the apex; seeds μ for the full fit.
// ============================================================================
struct PrefitRes { double mu = 0.135, sigma = 0.005, N = 0; bool ok = false; };

PrefitRes prefitGauss(TH1D* h) {
    PrefitRes r;
    const int b_lo = h->FindBin(kPreMin);
    const int b_hi = h->FindBin(kPreMax);

    int peak_bin = b_lo;
    double peak_max = h->GetBinContent(b_lo);
    for (int b = b_lo; b <= b_hi; ++b) {
        if (h->GetBinContent(b) > peak_max) {
            peak_max = h->GetBinContent(b);
            peak_bin = b;
        }
    }
    if (peak_max <= 0) {
        r.mu = h->GetBinCenter(peak_bin);
        return r;
    }

    auto* g = new TF1("gpre", "gaus", kPreMin, kPreMax);
    g->SetParameter(0, peak_max);
    g->SetParameter(1, h->GetBinCenter(peak_bin));
    g->SetParameter(2, 0.005);
    g->SetParLimits(0, 0.0, 1e12);
    g->SetParLimits(1, kPreMin, kPreMax);
    g->SetParLimits(2, 0.001, 0.020);

    int status = h->Fit(g, "RQN");
    if (status == 0) {
        r.N     = g->GetParameter(0);
        r.mu    = g->GetParameter(1);
        r.sigma = g->GetParameter(2);
        r.ok    = true;
    } else {
        r.mu = h->GetBinCenter(peak_bin);
    }
    delete g;
    return r;
}

// ============================================================================
// Build signal hist = all − CB, with CB = 2 √(N_++ N_--) per bin.
// ============================================================================
TH1D* buildSignal(TFile* f_all, TFile* f_pp, TFile* f_mm,
                  const std::string& hname, const std::string& tag)
{
    auto* h_all_in = (TH1D*)f_all->Get(hname.c_str());
    auto* h_pp_in  = (TH1D*)f_pp ->Get(hname.c_str());
    auto* h_mm_in  = (TH1D*)f_mm ->Get(hname.c_str());
    if (!h_all_in || !h_pp_in || !h_mm_in) return nullptr;

    auto* h_all = (TH1D*)h_all_in->Clone((hname + "_all_" + tag).c_str());
    auto* h_pp  = (TH1D*)h_pp_in ->Clone((hname + "_pp_"  + tag).c_str());
    auto* h_mm  = (TH1D*)h_mm_in ->Clone((hname + "_mm_"  + tag).c_str());
    h_all->SetDirectory(nullptr);
    h_pp ->SetDirectory(nullptr);
    h_mm ->SetDirectory(nullptr);

    auto* h_cb = (TH1D*)h_all->Clone((hname + "_cb_" + tag).c_str());
    h_cb->Reset();
    h_cb->SetDirectory(nullptr);
    for (int b = 1; b <= h_cb->GetNbinsX(); ++b) {
        const double n_pp = h_pp->GetBinContent(b);
        const double n_mm = h_mm->GetBinContent(b);
        const double e_pp = h_pp->GetBinError(b);
        const double e_mm = h_mm->GetBinError(b);
        double cb = 0.0, err = 0.0;
        if (n_pp > 0 && n_mm > 0) {
            cb = 2.0 * std::sqrt(n_pp * n_mm);
            const double t1 = e_pp * std::sqrt(n_mm / n_pp);
            const double t2 = e_mm * std::sqrt(n_pp / n_mm);
            err = std::sqrt(t1 * t1 + t2 * t2);
        }
        h_cb->SetBinContent(b, cb);
        h_cb->SetBinError(b, err);
    }

    auto* h_sig = (TH1D*)h_all->Clone((hname + "_sig_" + tag).c_str());
    h_sig->SetDirectory(nullptr);
    h_sig->Add(h_cb, -1.0);

    delete h_all;  delete h_pp;  delete h_mm;  delete h_cb;
    return h_sig;
}

// ============================================================================
// One slice fit + residual canvas + true integrals.
// ============================================================================
struct FitRes {
    double mu = 0,    mu_err = 0;
    double sigma = 0, sigma_err = 0;
    double alpha = 0, alpha_err = 0;
    double n = 0,     n_err = 0;
    double N_cb = 0,  N_cb_err = 0;
    double yield_full = 0, yield_full_err = 0;
    double yield_3sig = 0, yield_3sig_err = 0;
    double yield_data_1sig = 0, yield_data_1sig_err = 0;
    double yield_data_2sig = 0, yield_data_2sig_err = 0;
    double yield_data_3sig = 0, yield_data_3sig_err = 0;
    // Full-fit-range integral of residual — captures the entire CB left
    // tail that symmetric ±Nσ windows miss for asymmetric peaks.
    double yield_data_full = 0, yield_data_full_err = 0;
    // Asymmetric window [μ − 5σ, μ + 3σ] — wider on the left to catch the
    // CB tail, tighter on the right where the peak is Gaussian-like.
    double yield_data_m5p3sig = 0, yield_data_m5p3sig_err = 0;
    double chi2 = 0;  int ndf = 0;
    double a0 = 0, a1 = 0, a2 = 0, a3 = 0, b_exp = 0, c_exp = 0;
    double n_bg_g = 0, m_bg_g = 0, sig_bg_g = 0;
    bool   ok = false;
};

FitRes fitOneHist(TH1D* h, const std::string& base_name,
                  const std::string& canvas_title,
                  const std::string& out_basename,
                  double oa_lo,
                  std::vector<TCanvas*>& canvases,
                  std::vector<TCanvas*>& canvases_res,
                  double prev_alpha = -1.0,
                  double prev_sigma = -1.0,
                  double prev_c_exp = -1.0,
                  double prev_n_bg_g = -1.0,
                  double prev_m_bg_g = -1.0,
                  double prev_sig_bg_g = -1.0,
                  const FitRes* prev_full = nullptr)
{
    FitRes r;
    if (!h || h->GetEntries() == 0) return r;

    const double bw = h->GetBinWidth(1);

    double fitMin = kFitMin;
    if      (oa_lo >= kOA_FitMin006) fitMin = 0.06;
    else if (oa_lo >= kOA_FitMin005) fitMin = 0.05;

    double alpha_min = 0.3, alpha_max = 5.0;
    if (oa_lo >= kOA_AlphaTighten) {
        alpha_min = kAlphaTightMin;
        alpha_max = (oa_lo >= kOA_EnableP3) ? kAlphaTightMaxHi
                                            : kAlphaTightMaxMid;
    }

    PrefitRes pre = prefitGauss(h);

    auto* fit = new TF1(("fit_" + base_name).c_str(),
                        FitFunction, fitMin, kFitMax, 14);

    const bool useP3 = (oa_lo >= kOA_EnableP3);

    fit->SetParName(0,  "N_{CB}");
    fit->SetParName(1,  "#mu");
    fit->SetParName(2,  "#sigma");
    fit->SetParName(3,  "#alpha");
    fit->SetParName(4,  "n");
    fit->SetParName(5,  "a_{0}");
    fit->SetParName(6,  "a_{1}");
    fit->SetParName(7,  "a_{2}");
    fit->SetParName(8,  "b_{exp}");
    fit->SetParName(9,  "c_{exp}");
    fit->SetParName(10, "a_{3}");
    fit->SetParName(11, "N_{bg,G}");
    fit->SetParName(12, "m_{bg,G}");
    fit->SetParName(13, "#sigma_{bg,G}");

    // Scale bg starting values with the data magnitude — for slices ~10k
    // counts the legacy values were OK, but the integrated full-OA panel
    // is ~3M counts and MIGRAD converged to a wrong minimum (μ pinned at
    // 0.140). Use bin content near m=0.40 as the bg-level reference.
    const double bg_ref = std::max(1.0, h->GetBinContent(h->FindBin(0.40)));

    fit->SetParameter(0, h->GetMaximum());
    fit->SetParameter(1, pre.mu);
    fit->SetParameter(2, 0.010);
    fit->SetParameter(3, 1.5);
    fit->SetParameter(4, 3.0);
    fit->SetParameter(5, bg_ref);              // a0 ~ data at right sideband
    fit->SetParameter(6, -bg_ref * 5.0);       // a1
    fit->SetParameter(7,  bg_ref * 5.0);       // a2
    fit->SetParameter(8,  bg_ref * 10.0);      // b_exp
    fit->SetParameter(9, 8.0);
    fit->SetParameter(10, 0.0);
    const bool useGaussBg = (oa_lo < kOA_DisableGaussBg);
    fit->SetParameter(11, useGaussBg ? bg_ref * 0.5 : 0.0);
    fit->SetParameter(12, 0.09);   // m_bg_g — fixed center
    fit->SetParameter(13, 0.07);   // σ_bg_g — fixed width

    fit->SetParLimits(0, 0.0, 1e9);
    fit->SetParLimits(1, 0.128, 0.150);
    fit->SetParLimits(2, 0.002, 0.030);
    fit->SetParLimits(3, alpha_min, alpha_max);
    fit->FixParameter(4, 3.0);

    // Sequential handoff from previous slice — when prev_alpha/prev_sigma
    // are provided, seed the fit at those values and tighten the window so
    // a single MIGRAD jump can't move the peak parameters too far. Prevents
    // isolated outlier slices from changing α/σ much vs. their neighbours.
    // Start values are nudged off the boundary; ROOT can release a param
    // sitting exactly at a limit.
    if (prev_alpha > 0) {
        const double a_lo = std::max(alpha_min, prev_alpha - 0.07);
        const double a_hi = std::min(alpha_max, prev_alpha + 0.07);
        if (a_hi - a_lo > 1e-4) {
            fit->SetParLimits(3, a_lo, a_hi);
            fit->SetParameter(3, std::min(std::max(prev_alpha,
                                                   a_lo + 1e-3),
                                          a_hi - 1e-3));
        }
    }
    if (prev_sigma > 0) {
        const double s_lo = std::max(0.005, prev_sigma * 0.95);
        const double s_hi = std::min(0.030, prev_sigma * 1.05);
        if (s_hi - s_lo > 1e-5) {
            fit->SetParLimits(2, s_lo, s_hi);
            fit->SetParameter(2, std::min(std::max(prev_sigma,
                                                   s_lo + 1e-4),
                                          s_hi - 1e-4));
        }
    }
    if (prev_c_exp > 0) {
        const double c_lo = std::max(0.0,   prev_c_exp * 0.50);
        const double c_hi = std::min(100.0, prev_c_exp * 1.50);
        if (c_hi - c_lo > 1e-3) {
            fit->SetParLimits(9, c_lo, c_hi);
            fit->SetParameter(9, std::min(std::max(prev_c_exp,
                                                   c_lo + 1e-3),
                                          c_hi - 1e-3));
        }
    }
    (void)prev_m_bg_g; (void)prev_sig_bg_g;
    fit->SetParLimits(8, 0.0, 1e9);
    fit->SetParLimits(9, 0.0, 100.0);
    if (useP3) fit->SetParLimits(10, -1e6, 1e6);
    else       fit->FixParameter(10, 0.0);

    // Gauss-bg (par 11..13) — only N_bg_g is fitted; the Gauss shape
    // (center, width) is held FIXED globally so the bg profile stays
    // visually consistent slice-to-slice. With center/width free, MIGRAD
    // hops between local minima with similar χ² but very different
    // appearance (e.g. m_bg_g jumping 0.04↔0.11 between neighbours).
    // From OA ≥ kOA_DisableGaussBg the whole term is zeroed.
    if (useGaussBg) {
        fit->SetParLimits(11, bg_ref * 0.1, 1e9);
        fit->FixParameter(12, 0.09);
        fit->FixParameter(13, 0.07);
    } else {
        fit->FixParameter(11, 0.0);
        fit->FixParameter(12, 0.09);
        fit->FixParameter(13, 0.07);
    }

    // True double-fit when the caller supplies the previous slice's full
    // FitRes:
    //   Pass 1 — bg params FROZEN at prev values; only signal (CB) is free.
    //            This anchors the bg shape exactly to the previous slice
    //            and lets MIGRAD relax just the peak.
    //   Pass 2 — bg released with ±20% bounds around prev; final refine.
    // This stops MIGRAD from hopping to a different bg shape minimum even
    // when bounds alone fail to restrain it.
    int status = -1;
    if (prev_full) {
        // Pass 1: freeze bg
        fit->SetParameter(5,  prev_full->a0);    fit->FixParameter(5,  prev_full->a0);
        fit->SetParameter(6,  prev_full->a1);    fit->FixParameter(6,  prev_full->a1);
        fit->SetParameter(7,  prev_full->a2);    fit->FixParameter(7,  prev_full->a2);
        fit->SetParameter(8,  prev_full->b_exp); fit->FixParameter(8,  prev_full->b_exp);
        fit->SetParameter(9,  prev_full->c_exp); fit->FixParameter(9,  prev_full->c_exp);
        if (useP3) {
            fit->SetParameter(10, prev_full->a3); fit->FixParameter(10, prev_full->a3);
        }
        if (useGaussBg) {
            fit->SetParameter(11, prev_full->n_bg_g);
            fit->FixParameter(11, prev_full->n_bg_g);
        }
        h->Fit(fit, "RQ");

        // Pass 2: release bg with ±20% bounds around prev (still tight).
        auto releaseRel = [&](int p, double v, double rel) {
            if (v == 0.0) { fit->ReleaseParameter(p); return; }
            const double lo = (v > 0) ? v * (1 - rel) : v * (1 + rel);
            const double hi = (v > 0) ? v * (1 + rel) : v * (1 - rel);
            const double a = std::min(lo, hi);
            const double b = std::max(lo, hi);
            fit->ReleaseParameter(p);
            fit->SetParLimits(p, a, b);
        };
        releaseRel(5, prev_full->a0,    0.20);
        releaseRel(6, prev_full->a1,    0.20);
        releaseRel(7, prev_full->a2,    0.20);
        releaseRel(8, prev_full->b_exp, 0.20);
        releaseRel(9, prev_full->c_exp, 0.20);
        if (useP3)      releaseRel(10, prev_full->a3,    0.20);
        if (useGaussBg) releaseRel(11, prev_full->n_bg_g, 0.30);
        status = h->Fit(fit, "RQ");
    } else {
        // Standard single-pass fit + one refine for slices without prev.
        status = h->Fit(fit, "RQ");
        h->Fit(fit, "RQ");
    }
    r.ok = (status == 0);

    r.N_cb    = fit->GetParameter(0); r.N_cb_err  = fit->GetParError(0);
    r.mu      = fit->GetParameter(1); r.mu_err    = fit->GetParError(1);
    r.sigma   = fit->GetParameter(2); r.sigma_err = fit->GetParError(2);
    r.alpha   = fit->GetParameter(3); r.alpha_err = fit->GetParError(3);
    r.n       = fit->GetParameter(4); r.n_err     = fit->GetParError(4);
    r.a0      = fit->GetParameter(5);
    r.a1      = fit->GetParameter(6);
    r.a2      = fit->GetParameter(7);
    r.b_exp   = fit->GetParameter(8);
    r.c_exp   = fit->GetParameter(9);
    r.a3      = fit->GetParameter(10);
    r.n_bg_g  = fit->GetParameter(11);
    r.m_bg_g  = fit->GetParameter(12);
    r.sig_bg_g= fit->GetParameter(13);
    r.chi2    = fit->GetChisquare();
    r.ndf     = fit->GetNDF();

    auto* cbOnly = new TF1(("cb_" + base_name).c_str(),
                           CrystalBallLeft, fitMin, kFitMax, 5);
    for (int i = 0; i < 5; ++i) cbOnly->SetParameter(i, fit->GetParameter(i));
    cbOnly->SetLineColor(kBlue);
    cbOnly->SetLineWidth(2);
    cbOnly->SetLineStyle(2);

    r.yield_full     = cbOnly->Integral(fitMin, kFitMax) / bw;
    r.yield_full_err = (r.N_cb > 0) ? r.yield_full * (r.N_cb_err / r.N_cb) : 0.0;
    const double w_lo = std::max(fitMin, r.mu - 3.0 * r.sigma);
    const double w_hi = std::min(kFitMax, r.mu + 3.0 * r.sigma);
    r.yield_3sig     = cbOnly->Integral(w_lo, w_hi) / bw;
    r.yield_3sig_err = (r.N_cb > 0) ? r.yield_3sig * (r.N_cb_err / r.N_cb) : 0.0;

    auto* bgOnly = new TF1(("bg_" + base_name).c_str(),
        "[0] + [1]*x + [2]*x*x + [5]*x*x*x + [3]*exp(-[4]*x)"
        " + [6]*exp(-0.5*((x-[7])/[8])*((x-[7])/[8]))",
        fitMin, kFitMax);
    bgOnly->SetParameter(0, r.a0);
    bgOnly->SetParameter(1, r.a1);
    bgOnly->SetParameter(2, r.a2);
    bgOnly->SetParameter(3, r.b_exp);
    bgOnly->SetParameter(4, r.c_exp);
    bgOnly->SetParameter(5, r.a3);
    bgOnly->SetParameter(6, r.n_bg_g);
    bgOnly->SetParameter(7, r.m_bg_g);
    bgOnly->SetParameter(8, r.sig_bg_g);
    bgOnly->SetLineColor(kGreen + 2);
    bgOnly->SetLineWidth(2);
    bgOnly->SetLineStyle(3);

    auto* c = new TCanvas(("c_" + base_name).c_str(), canvas_title.c_str(), 900, 700);
    c->SetMargin(0.12, 0.05, 0.12, 0.08);
    h->SetMarkerStyle(20);
    h->SetMarkerSize(0.7);
    h->SetMarkerColor(kBlue);   // signal markers in blue (consistent with slices_pi0.C)
    h->SetLineColor(kBlue);
    h->SetTitle(canvas_title.c_str());
    h->GetXaxis()->SetRangeUser(0.0, 0.46);
    // Y axis from 0 — without this, ROOT auto-scales below 0 when the
    // signal residual fluctuates negative (e.g. at low-stat slices like
    // 12.0–12.2 the auto-Y dropped the canvas baseline well below 0).
    h->SetMinimum(0.0);

    h->Draw("E");
    fit->SetLineColor(kRed);
    fit->SetLineWidth(2);
    fit->Draw("same");
    cbOnly->Draw("same");
    bgOnly->Draw("same");

    auto* leg = new TLegend(0.55, 0.62, 0.94, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.030);
    const char* bg_label_full =
        useP3 ? (useGaussBg ? "CB + P_{3} + exp + G_{bg}" : "CB + P_{3} + exp")
              : (useGaussBg ? "CB + P_{2} + exp + G_{bg}" : "CB + P_{2} + exp");
    const char* bg_label_only =
        useP3 ? (useGaussBg ? "bg: P_{3} + exp + G_{bg}" : "bg: P_{3} + exp")
              : (useGaussBg ? "bg: P_{2} + exp + G_{bg}" : "bg: P_{2} + exp");
    leg->AddEntry(h,      "signal (all #minus CB)",  "lpe");
    leg->AddEntry(fit,    bg_label_full,             "l");
    leg->AddEntry(cbOnly, "Crystal Ball: #pi^{0}",   "l");
    leg->AddEntry(bgOnly, bg_label_only,             "l");
    leg->Draw();

    auto* tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.030);
    tex->DrawLatex(0.15, 0.86, TString::Format("#mu = %.4f #pm %.4f GeV/c^{2}", r.mu, r.mu_err));
    tex->DrawLatex(0.15, 0.82, TString::Format("#sigma = %.4f #pm %.4f GeV/c^{2}", r.sigma, r.sigma_err));
    tex->DrawLatex(0.15, 0.78, TString::Format("#alpha = %.2f #pm %.2f", r.alpha, r.alpha_err));
    tex->DrawLatex(0.15, 0.74, TString::Format("n = %.2f (fixed)", r.n));
    tex->DrawLatex(0.15, 0.68, TString::Format("yield(fit range) = %.3g #pm %.3g",
                                               r.yield_full, r.yield_full_err));
    tex->DrawLatex(0.15, 0.64, TString::Format("yield(#mu#pm3#sigma) = %.3g #pm %.3g",
                                               r.yield_3sig, r.yield_3sig_err));
    tex->DrawLatex(0.15, 0.58, TString::Format("#chi^{2}/ndf = %.2f / %d = %.2f",
                                               r.chi2, r.ndf,
                                               r.ndf > 0 ? r.chi2 / r.ndf : 0.0));

    c->Update();
    const std::string per = "plots/output/fit_pi0_" + out_basename;
    c->SaveAs((per + ".pdf").c_str());
    c->SaveAs((per + ".png").c_str());
    canvases.push_back(c);

    // ------------------------------------------------------------------------
    // Residual: signal − bg (true CB shape) + ±Nσ true integrals.
    // ------------------------------------------------------------------------
    TH1D* h_res = (TH1D*)h->Clone(("res_" + base_name).c_str());
    h_res->SetDirectory(nullptr);

    for (int b = 1; b <= h_res->GetNbinsX(); ++b) {
        const double m = h_res->GetBinCenter(b);
        const double bg_val = bgOnly->Eval(m);
        h_res->SetBinContent(b, h->GetBinContent(b) - bg_val);
    }

    // Clip integration window to the fit range — outside [fitMin, kFitMax]
    // the bg model is extrapolated and the residual contains structure
    // (e.g. m → 0 dropoff) that would falsely inflate yields if included.
    auto trueIntegral = [&](double m_lo, double m_hi) {
        const double m_lo_c = std::max(m_lo, fitMin);
        const double m_hi_c = std::min(m_hi, kFitMax);
        if (m_hi_c <= m_lo_c) return std::make_pair(0.0, 0.0);
        const int b_lo = std::max(1, h_res->FindBin(m_lo_c));
        const int b_hi = std::min(h_res->GetNbinsX(), h_res->FindBin(m_hi_c));
        double err = 0.0;
        const double v = h_res->IntegralAndError(b_lo, b_hi, err);
        return std::make_pair(v, err);
    };

    {
        auto p = trueIntegral(r.mu - 1.0 * r.sigma, r.mu + 1.0 * r.sigma);
        r.yield_data_1sig = p.first;  r.yield_data_1sig_err = p.second;
    }
    {
        auto p = trueIntegral(r.mu - 2.0 * r.sigma, r.mu + 2.0 * r.sigma);
        r.yield_data_2sig = p.first;  r.yield_data_2sig_err = p.second;
    }
    {
        auto p = trueIntegral(r.mu - 3.0 * r.sigma, r.mu + 3.0 * r.sigma);
        r.yield_data_3sig = p.first;  r.yield_data_3sig_err = p.second;
    }
    {
        auto p = trueIntegral(fitMin, kFitMax);
        r.yield_data_full = p.first;  r.yield_data_full_err = p.second;
    }
    {
        auto p = trueIntegral(r.mu - 5.0 * r.sigma, r.mu + 3.0 * r.sigma);
        r.yield_data_m5p3sig = p.first;  r.yield_data_m5p3sig_err = p.second;
    }

    if (auto* funcs = h_res->GetListOfFunctions()) funcs->Clear();
    for (int b = 1; b <= h_res->GetNbinsX(); ++b) {
        if (h_res->GetBinContent(b) < 0.0) h_res->SetBinContent(b, 0.0);
    }

    auto* c2 = new TCanvas(("c_res_" + base_name).c_str(),
                           ("residual " + canvas_title).c_str(), 900, 700);
    c2->SetMargin(0.12, 0.05, 0.12, 0.08);
    h_res->SetMarkerStyle(20);
    h_res->SetMarkerSize(0.7);
    h_res->SetMarkerColor(kBlue);
    h_res->SetLineColor(kBlue);
    h_res->SetTitle(("residual: " + canvas_title).c_str());
    h_res->GetXaxis()->SetRangeUser(0.0, 0.30);
    h_res->SetMinimum(0.0);
    h_res->SetMaximum(h_res->GetMaximum() * 1.15);

    h_res->Draw("E");
    cbOnly->Draw("same");

    {
        const double y_lo = 0.0;
        const double y_hi = h_res->GetMaximum();
        const int colors[3] = {kGray + 1, kGray + 2, kGray + 3};
        for (int n = 1; n <= 3; ++n) {
            for (int sign = -1; sign <= 1; sign += 2) {
                auto* line = new TLine(r.mu + sign * n * r.sigma, y_lo,
                                       r.mu + sign * n * r.sigma, y_hi);
                line->SetLineColor(colors[n - 1]);
                line->SetLineStyle(2);
                line->SetLineWidth(1);
                line->Draw();
            }
        }
    }

    auto* leg2 = new TLegend(0.55, 0.78, 0.94, 0.90);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetTextSize(0.030);
    leg2->AddEntry(h_res,  "signal #minus bg",      "lpe");
    leg2->AddEntry(cbOnly, "Crystal Ball: #pi^{0}", "l");
    leg2->Draw();

    auto* tex2 = new TLatex();
    tex2->SetNDC();
    tex2->SetTextSize(0.030);
    tex2->DrawLatex(0.15, 0.86,
        TString::Format("yield(#mu#pm1#sigma) = %.3g #pm %.3g (sig#minusbg)",
                        r.yield_data_1sig, r.yield_data_1sig_err));
    tex2->DrawLatex(0.15, 0.82,
        TString::Format("yield(#mu#pm2#sigma) = %.3g #pm %.3g (sig#minusbg)",
                        r.yield_data_2sig, r.yield_data_2sig_err));
    tex2->DrawLatex(0.15, 0.78,
        TString::Format("yield(#mu#pm3#sigma) = %.3g #pm %.3g (sig#minusbg)",
                        r.yield_data_3sig, r.yield_data_3sig_err));
    tex2->DrawLatex(0.15, 0.74,
        TString::Format("yield([#mu#minus5#sigma, #mu+3#sigma]) = %.3g #pm %.3g (sig#minusbg)",
                        r.yield_data_m5p3sig, r.yield_data_m5p3sig_err));
    tex2->DrawLatex(0.15, 0.70,
        TString::Format("yield(full fit range) = %.3g #pm %.3g (sig#minusbg)",
                        r.yield_data_full, r.yield_data_full_err));
    tex2->DrawLatex(0.15, 0.64,
        TString::Format("#mu = %.4f, #sigma = %.4f GeV/c^{2}", r.mu, r.sigma));

    c2->Update();
    const std::string per_res = "plots/output/fit_pi0_residual_" + out_basename;
    c2->SaveAs((per_res + ".pdf").c_str());
    c2->SaveAs((per_res + ".png").c_str());
    canvases_res.push_back(c2);

    return r;
}

void fit_pi0(const char* flavour = "rec") {

    std::string fl = flavour;
    for (auto& c : fl) c = std::tolower(c);
    std::string var_stem;
    if      (fl == "rec") var_stem = "m_epemg";
    else if (fl == "cor") var_stem = "m_epemg_cor";
    else { std::cerr << "Unknown flavour '" << flavour
                     << "' (expected 'rec' or 'cor')\n"; return; }

    TFile* f_all = TFile::Open("research_epem.root", "READ");
    TFile* f_pp  = TFile::Open("research_epep.root", "READ");
    TFile* f_mm  = TFile::Open("research_emem.root", "READ");
    if (!f_all || f_all->IsZombie()) { std::cerr << "Cannot open research_epem.root\n"; return; }
    if (!f_pp  || f_pp ->IsZombie()) { std::cerr << "Cannot open research_epep.root\n"; return; }
    if (!f_mm  || f_mm ->IsZombie()) { std::cerr << "Cannot open research_emem.root\n"; return; }

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    gSystem->mkdir("plots/output", kTRUE);
    const std::string pdf_multi     = "plots/output/fit_pi0_" + fl + "_all.pdf";
    const std::string pdf_multi_res = "plots/output/fit_pi0_" + fl + "_residual_all.pdf";

    const std::string fres_path = "fit_results_" + fl + ".root";
    TFile* fres = TFile::Open(fres_path.c_str(), "RECREATE");
    auto* tres = new TTree("fit_results",
        "Crystal Ball + bg fits per OA slice (signal = all − CB)");

    char   name[64]   = {0};
    int    panel_idx  = 0;
    float  oa_lo = 0, oa_hi = 0;
    float  mu = 0, mu_err = 0, sigma = 0, sigma_err = 0;
    float  alpha = 0, alpha_err = 0, npar = 0, npar_err = 0;
    float  N_cb = 0, N_cb_err = 0;
    float  yield_full = 0, yield_full_err = 0;
    float  yield_3sig = 0, yield_3sig_err = 0;
    float  yield_data_1sig = 0, yield_data_1sig_err = 0;
    float  yield_data_2sig = 0, yield_data_2sig_err = 0;
    float  yield_data_3sig = 0, yield_data_3sig_err = 0;
    float  yield_data_full = 0, yield_data_full_err = 0;
    float  yield_data_m5p3sig = 0, yield_data_m5p3sig_err = 0;
    float  chi2 = 0; int ndf = 0; int ok = 0;
    float  a0 = 0, a1 = 0, a2 = 0, a3 = 0, b_exp = 0, c_exp = 0;
    float  n_bg_g = 0, m_bg_g = 0, sig_bg_g = 0;

    tres->Branch("name",            name,             "name/C");
    tres->Branch("panel_idx",       &panel_idx,       "panel_idx/I");
    tres->Branch("oa_lo",           &oa_lo,           "oa_lo/F");
    tres->Branch("oa_hi",           &oa_hi,           "oa_hi/F");
    tres->Branch("mu",              &mu,              "mu/F");
    tres->Branch("mu_err",          &mu_err,          "mu_err/F");
    tres->Branch("sigma",           &sigma,           "sigma/F");
    tres->Branch("sigma_err",       &sigma_err,       "sigma_err/F");
    tres->Branch("alpha",           &alpha,           "alpha/F");
    tres->Branch("alpha_err",       &alpha_err,       "alpha_err/F");
    tres->Branch("n",               &npar,            "n/F");
    tres->Branch("n_err",           &npar_err,        "n_err/F");
    tres->Branch("N_cb",            &N_cb,            "N_cb/F");
    tres->Branch("N_cb_err",        &N_cb_err,        "N_cb_err/F");
    tres->Branch("yield_full",      &yield_full,      "yield_full/F");
    tres->Branch("yield_full_err",  &yield_full_err,  "yield_full_err/F");
    tres->Branch("yield_3sig",      &yield_3sig,      "yield_3sig/F");
    tres->Branch("yield_3sig_err",  &yield_3sig_err,  "yield_3sig_err/F");
    tres->Branch("yield_data_1sig",     &yield_data_1sig,     "yield_data_1sig/F");
    tres->Branch("yield_data_1sig_err", &yield_data_1sig_err, "yield_data_1sig_err/F");
    tres->Branch("yield_data_2sig",     &yield_data_2sig,     "yield_data_2sig/F");
    tres->Branch("yield_data_2sig_err", &yield_data_2sig_err, "yield_data_2sig_err/F");
    tres->Branch("yield_data_3sig",     &yield_data_3sig,     "yield_data_3sig/F");
    tres->Branch("yield_data_3sig_err", &yield_data_3sig_err, "yield_data_3sig_err/F");
    tres->Branch("yield_data_full",     &yield_data_full,     "yield_data_full/F");
    tres->Branch("yield_data_full_err", &yield_data_full_err, "yield_data_full_err/F");
    tres->Branch("yield_data_m5p3sig",     &yield_data_m5p3sig,     "yield_data_m5p3sig/F");
    tres->Branch("yield_data_m5p3sig_err", &yield_data_m5p3sig_err, "yield_data_m5p3sig_err/F");
    tres->Branch("chi2",            &chi2,            "chi2/F");
    tres->Branch("ndf",             &ndf,             "ndf/I");
    tres->Branch("ok",              &ok,              "ok/I");
    tres->Branch("a0",              &a0,              "a0/F");
    tres->Branch("a1",              &a1,              "a1/F");
    tres->Branch("a2",              &a2,              "a2/F");
    tres->Branch("a3",              &a3,              "a3/F");
    tres->Branch("b_exp",           &b_exp,           "b_exp/F");
    tres->Branch("c_exp",           &c_exp,           "c_exp/F");
    tres->Branch("n_bg_g",          &n_bg_g,          "n_bg_g/F");
    tres->Branch("m_bg_g",          &m_bg_g,          "m_bg_g/F");
    tres->Branch("sig_bg_g",        &sig_bg_g,        "sig_bg_g/F");

    auto fillRow = [&](const FitRes& r, const std::string& nm,
                       int idx, double lo, double hi) {
        std::snprintf(name, sizeof(name), "%s", nm.c_str());
        panel_idx = idx;
        oa_lo = lo;  oa_hi = hi;
        mu = r.mu;             mu_err = r.mu_err;
        sigma = r.sigma;       sigma_err = r.sigma_err;
        alpha = r.alpha;       alpha_err = r.alpha_err;
        npar = r.n;            npar_err = r.n_err;
        N_cb = r.N_cb;         N_cb_err = r.N_cb_err;
        yield_full = r.yield_full;     yield_full_err = r.yield_full_err;
        yield_3sig = r.yield_3sig;     yield_3sig_err = r.yield_3sig_err;
        yield_data_1sig = r.yield_data_1sig;  yield_data_1sig_err = r.yield_data_1sig_err;
        yield_data_2sig = r.yield_data_2sig;  yield_data_2sig_err = r.yield_data_2sig_err;
        yield_data_3sig = r.yield_data_3sig;  yield_data_3sig_err = r.yield_data_3sig_err;
        yield_data_full = r.yield_data_full;  yield_data_full_err = r.yield_data_full_err;
        yield_data_m5p3sig = r.yield_data_m5p3sig;  yield_data_m5p3sig_err = r.yield_data_m5p3sig_err;
        chi2 = r.chi2;         ndf = r.ndf;
        ok = r.ok ? 1 : 0;
        a0 = r.a0;  a1 = r.a1;  a2 = r.a2;  a3 = r.a3;
        b_exp = r.b_exp;  c_exp = r.c_exp;
        n_bg_g = r.n_bg_g;  m_bg_g = r.m_bg_g;  sig_bg_g = r.sig_bg_g;
        tres->Fill();
    };

    std::vector<TCanvas*> canvases;
    std::vector<TCanvas*> canvases_res;
    canvases.reserve(60);
    canvases_res.reserve(60);

    const int n_slices = static_cast<int>(std::round((kSliceMax - kSliceMin) / kSliceStep));

    // Panel 0: full-range integrated (QA only).
    {
        const std::string hname = var_stem + "_full";
        TH1D* h_sig = buildSignal(f_all, f_pp, f_mm, hname, fl);
        if (h_sig) {
            const TString ctitle = TString::Format(
                "M(e^{+}e^{-}#gamma) %s — full OA range [%.1f, %.1f] deg "
                "(signal = all #minus CB);"
                "M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Counts",
                fl.c_str(), kSliceMin, kSliceMax);
            FitRes r = fitOneHist(h_sig, hname, ctitle.Data(),
                                  fl + "_full", kSliceMin,
                                  canvases, canvases_res);
            fillRow(r, hname, 0, kSliceMin, kSliceMax);
            std::cout << "  full: yield(fit) = " << r.yield_full
                      << "   chi2/ndf = "
                      << (r.ndf > 0 ? r.chi2/r.ndf : 0.0) << "\n";
        }
    }

    // Smoothed sequential α/σ handoff — each slice seeds the fit with the
    // moving average of the last few converged α and σ values. Smoothing
    // suppresses isolated MIGRAD jumps without freezing the natural drift
    // of σ with OA (resolution improves slightly toward higher OA).
    constexpr int kHandoffWindow = 3;
    std::vector<double> recent_alpha, recent_sigma, recent_c_exp;
    std::vector<double> recent_m_bg_g, recent_sig_bg_g;

    // Local reverse-engineering: slices in [kReverseRangeMin, kReverseRangeMax)
    // are SKIPPED in the forward pass. After the forward pass completes, we
    // fit them in REVERSE order seeded by the slice immediately after the
    // range (i.e. the 2.4-2.6 anchor → 2.2-2.4, then 2.2-2.4 → 2.0-2.2).
    // Used only for this isolated 2-slice anomaly; rest of the OA range
    // fits independently without forward propagation (which broke 3.6-3.8).
    constexpr double kReverseRangeMin = 1.6;
    constexpr double kReverseRangeMax = 2.4;

    // Pinpoint forward propagation: at this slice, copy bg shape from the
    // immediately preceding fitted slice (with amplitudes scaled by data
    // magnitude). Used to anchor an isolated slice whose free fit
    // converges to a degenerate minimum (here: bg cut too aggressively).
    constexpr double kForwardFixLo = 5.8;
    struct PendingSlice {
        int    i;
        TH1D*  h;
        std::string hname;
        std::string ctitle;
        double lo, hi;
    };
    std::vector<PendingSlice> pending;
    FitRes anchor_res;     // result for the slice just after the reverse range
    bool   anchor_ok = false;
    FitRes last_fwd_res;   // most recent successful forward-fit (for kForwardFixLo)
    bool   last_fwd_ok = false;
    auto smoothedAvg = [](const std::vector<double>& v) {
        if (v.empty()) return -1.0;
        double s = 0.0;
        for (double x : v) s += x;
        return s / v.size();
    };

    for (int i = 0; i < n_slices; ++i) {
        const double lo = kSliceMin + i * kSliceStep;
        const double hi = kSliceMin + (i + 1) * kSliceStep;
        const std::string suff = "oa_" + fmtEdge(lo) + "_" + fmtEdge(hi);
        const std::string hname = var_stem + "_" + suff;

        TH1D* h_sig = buildSignal(f_all, f_pp, f_mm, hname, fl);
        if (!h_sig) { std::cerr << "  missing " << hname << "\n"; continue; }

        const TString ctitle = TString::Format(
            "M(e^{+}e^{-}#gamma) %s, OA #in [%.1f, %.1f] deg (signal = all #minus CB);"
            "M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Counts",
            fl.c_str(), lo, hi);

        // Skip slices in the reverse range — fitted in a backward pass below.
        if (lo >= kReverseRangeMin - 1e-6 && lo < kReverseRangeMax - 1e-6) {
            pending.push_back({i, h_sig, hname, ctitle.Data(), lo, hi});
            continue;
        }

        // Pinpoint forward propagation at kForwardFixLo: scale bg amplitudes
        // by data ratio (m=0.40 probe) and pass last_fwd_res to fitOneHist.
        const bool use_fwd_fix =
            (std::fabs(lo - kForwardFixLo) < 1e-4) && last_fwd_ok;
        FitRes scaled_prev;
        const FitRes* fwd_prev_arg = nullptr;
        if (use_fwd_fix) {
            const double pend_bgref = std::max(1.0,
                h_sig->GetBinContent(h_sig->FindBin(0.40)));
            const double m0 = 0.40;
            const double anchor_bg_at = last_fwd_res.a0
                + last_fwd_res.a1 * m0
                + last_fwd_res.a2 * m0 * m0
                + last_fwd_res.a3 * m0 * m0 * m0
                + last_fwd_res.b_exp * std::exp(-last_fwd_res.c_exp * m0);
            const double scale_amp = (anchor_bg_at > 1.0)
                ? pend_bgref / anchor_bg_at : 1.0;
            scaled_prev = last_fwd_res;
            scaled_prev.a0     *= scale_amp;
            scaled_prev.a1     *= scale_amp;
            scaled_prev.a2     *= scale_amp;
            scaled_prev.a3     *= scale_amp;
            scaled_prev.b_exp  *= scale_amp;
            scaled_prev.n_bg_g *= scale_amp;
            fwd_prev_arg = &scaled_prev;
            std::cout << "  [fwd-fix at oa=" << lo
                      << "] using prev with bg×" << scale_amp << "\n";
        }

        FitRes r = fitOneHist(h_sig, hname, ctitle.Data(),
                              fl + "_slice_" + suff, lo,
                              canvases, canvases_res,
                              smoothedAvg(recent_alpha),
                              smoothedAvg(recent_sigma),
                              smoothedAvg(recent_c_exp),
                              -1.0,
                              smoothedAvg(recent_m_bg_g),
                              smoothedAvg(recent_sig_bg_g),
                              fwd_prev_arg);
        fillRow(r, hname, 1 + i, lo, hi);

        if (r.ok) {
            last_fwd_res = r;
            last_fwd_ok  = true;
        }

        if (r.ok && r.alpha > 0 && r.sigma > 0) {
            recent_alpha.push_back(r.alpha);
            recent_sigma.push_back(r.sigma);
            if (r.c_exp > 0) recent_c_exp.push_back(r.c_exp);
            if (r.n_bg_g > 0 && r.m_bg_g > 0 && r.sig_bg_g > 0) {
                recent_m_bg_g.push_back(r.m_bg_g);
                recent_sig_bg_g.push_back(r.sig_bg_g);
            }
            auto trim = [&](std::vector<double>& v) {
                if ((int)v.size() > kHandoffWindow) v.erase(v.begin());
            };
            trim(recent_alpha);  trim(recent_sigma);  trim(recent_c_exp);
            trim(recent_m_bg_g); trim(recent_sig_bg_g);
        }

        // Capture the anchor slice — the first one fitted just past the
        // reverse range — to use as bg-shape template for backward fitting.
        if (!anchor_ok && lo >= kReverseRangeMax - 1e-6 && r.ok) {
            anchor_res = r;
            anchor_ok  = true;
        }

        std::cout << "  " << hname
                  << ": yield(fit) = " << std::setw(10) << r.yield_full
                  << "   yield(#pm3#sigma) = " << std::setw(10) << r.yield_3sig
                  << "   μ = " << r.mu << "   σ = " << r.sigma
                  << "   χ²/ndf = "
                  << (r.ndf > 0 ? r.chi2/r.ndf : 0.0) << "\n";
    }

    // ===========================================================================
    // Reverse pass: fit pending slices in DESCENDING OA order, each seeded by
    // the previously fitted slice's full FitRes (anchor for the last one in the
    // pending list). Locks bg shape to the post-range anchor's bg.
    // ===========================================================================
    if (anchor_ok && !pending.empty()) {
        std::cout << "\n=== Reverse pass for OA #in [" << kReverseRangeMin
                  << ", " << kReverseRangeMax << ") — anchor: oa_lo = "
                  << kReverseRangeMax << " ===\n";
        // Scale anchor's bg amplitudes by ratio of data magnitudes between
        // pending slice and anchor (bg shape preserved, amplitude tracks
        // local data statistics — at the peak the bg is larger than in
        // the surrounding slices used as anchor).
        const double anchor_bgref = std::max(1.0,
            anchor_res.a0 + anchor_res.b_exp);  // approximate bg level
        FitRes prev_for_reverse = anchor_res;
        for (auto it = pending.rbegin(); it != pending.rend(); ++it) {
            const auto& p = *it;
            const std::string suff = "oa_" + fmtEdge(p.lo) + "_" + fmtEdge(p.hi);

            // Estimate scale = bg-level ratio at m=0.40 between pending
            // slice and anchor's frozen reference. Anchor used bg_ref from
            // its own data, so we use the same probe for the pending one.
            const double pend_bgref = std::max(1.0,
                p.h->GetBinContent(p.h->FindBin(0.40)));
            const double anchor_probe = std::max(1.0,
                p.h->GetBinContent(p.h->FindBin(0.40)));  // same hist in anchor probe
            (void)anchor_probe;
            const double scale = pend_bgref /
                std::max(1.0, anchor_res.a0 + anchor_res.b_exp);
            // simpler: scale = pend.bg_at_0.40 / anchor's same.
            // anchor's bg at m=0.40 ≈ a0 + a1·0.4 + a2·0.16 + a3·0.064
            //                         + b_exp·exp(-c·0.4) + n_bg_g·tinyGauss
            const double m0 = 0.40;
            const double anchor_bg_at = anchor_res.a0
                + anchor_res.a1 * m0
                + anchor_res.a2 * m0 * m0
                + anchor_res.a3 * m0 * m0 * m0
                + anchor_res.b_exp * std::exp(-anchor_res.c_exp * m0);
            const double scale_amp = (anchor_bg_at > 1.0)
                ? pend_bgref / anchor_bg_at : 1.0;

            FitRes scaled = anchor_res;
            scaled.a0    *= scale_amp;
            scaled.a1    *= scale_amp;
            scaled.a2    *= scale_amp;
            scaled.a3    *= scale_amp;
            scaled.b_exp *= scale_amp;
            scaled.n_bg_g*= scale_amp;
            // c_exp, m_bg_g, sig_bg_g (shape) unchanged.
            (void)anchor_bgref;

            FitRes r = fitOneHist(p.h, p.hname, p.ctitle,
                                  fl + "_slice_" + suff, p.lo,
                                  canvases, canvases_res,
                                  -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
                                  &scaled);
            fillRow(r, p.hname, 1 + p.i, p.lo, p.hi);
            if (r.ok) prev_for_reverse = r;
            std::cout << "  " << p.hname
                      << " (reverse, bg×" << scale_amp << "): μ = " << r.mu
                      << "   σ = " << r.sigma
                      << "   χ²/ndf = "
                      << (r.ndf > 0 ? r.chi2/r.ndf : 0.0) << "\n";
        }
        (void)prev_for_reverse;
    }

    for (size_t i = 0; i < canvases.size(); ++i) {
        auto* c = canvases[i];
        if (i == 0)                        c->Print((pdf_multi + "(").c_str(), "pdf");
        else if (i == canvases.size() - 1) c->Print((pdf_multi + ")").c_str(), "pdf");
        else                               c->Print(pdf_multi.c_str(),         "pdf");
    }
    for (size_t i = 0; i < canvases_res.size(); ++i) {
        auto* c = canvases_res[i];
        if (i == 0)                            c->Print((pdf_multi_res + "(").c_str(), "pdf");
        else if (i == canvases_res.size() - 1) c->Print((pdf_multi_res + ")").c_str(), "pdf");
        else                                   c->Print(pdf_multi_res.c_str(),         "pdf");
    }

    fres->cd();
    tres->Write();
    fres->Close();

    for (auto* c : canvases)     delete c;
    for (auto* c : canvases_res) delete c;
    f_all->Close();  f_pp->Close();  f_mm->Close();

    std::cout << "\nResults TTree: research/" << fres_path
              << "  (TTree 'fit_results')\n";
    std::cout << "Multipage PDFs:\n"
              << "  fit:      " << pdf_multi
              << "  (" << canvases.size() << " pages)\n"
              << "  residual: " << pdf_multi_res
              << "  (" << canvases_res.size() << " pages)\n";
}
