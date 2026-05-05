// research/plots/fit_eta.C — Crystal Ball + (P2 + exp) fit on the SIGNAL
// histogram (all − CB) per OA slice, in the η-Dalitz mass region.
//
// Same family as fit_pi0.C — single-sided Crystal Ball with the tail on the
// LEFT to capture bremsstrahlung losses on REC/COR; asymmetric peak shape
// describes the low-mass side better than a symmetric Gaussian.
//
//   F(m) = N · CrystalBall_left(m; μ, σ, α, n)
//        + a₀ + a₁·m + a₂·m²
//        + b · exp(-c·m)
//
// Inputs (reads from research/):
//   research_epem.root, research_epep.root, research_emem.root
// Builds CB = 2 √(N_++ N_--), signal = all − CB per bin.
//
// Usage (from research/):
//   root -l -b -q plots/fit_eta.C            # REC (default)
//   root -l -b -q 'plots/fit_eta.C("cor")'   # COR

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
    constexpr double kSliceStep = 0.5;   // coarser binning for η (low-stats)

    // Fit window targets the η peak.
    constexpr double kFitMin = 0.32;
    constexpr double kFitMax = 0.70;

    // Pre-fit Gaussian apex window — anchors μ.
    constexpr double kPreMin = 0.535;
    constexpr double kPreMax = 0.560;

    // CB peak constraints (η-tuned). η PDG = 0.5478; resolution wider
    // than π⁰ at this mass scale, so σ window is generous.
    constexpr double kMuMin     = 0.535;
    constexpr double kMuMax     = 0.565;
    constexpr double kSigInit   = 0.020;
    constexpr double kSigMin    = 0.005;
    constexpr double kSigMax    = 0.050;

    // CB tail params. n fixed (typical recipe); α generous since η peak shape
    // depends on flavour (REC/COR has brems tail, less so the further into m).
    constexpr double kAlphaMin  = 0.5;
    constexpr double kAlphaMax  = 5.0;
    constexpr double kAlphaInit = 1.5;
    constexpr double kFixedN    = 3.0;

    constexpr double kEtaMass   = 0.547;   // PDG η: 0.5478

    // Display zoom for canvases.
    constexpr double kDispMin = 0.30;
    constexpr double kDispMax = 0.80;

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
//   par[0] = N, par[1] = μ, par[2] = σ, par[3] = α, par[4] = n
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
// Full model: CB + P2 + exp (10 free params, n fixed → 9 effective).
//   par[0..4] = CB (N, μ, σ, α, n)
//   par[5..7] = P2 (a0, a1, a2)
//   par[8..9] = exp (b, c)  →  b·exp(-c·m)
// ============================================================================
double FitFunction(double *x, double *par) {
    const double m  = x[0];
    const double cb = CrystalBallLeft(x, par);
    const double p2 = par[5] + par[6]*m + par[7]*m*m;
    const double ex = par[8] * std::exp(-par[9]*m);
    return cb + p2 + ex;
}

// ============================================================================
// Pre-fit Gaussian on the apex; seeds μ for the full fit.
// ============================================================================
struct PrefitRes { double mu = kEtaMass, sigma = 0.015, N = 0; bool ok = false; };

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

    auto* g = new TF1("gpre_eta", "gaus", kPreMin, kPreMax);
    g->SetParameter(0, peak_max);
    g->SetParameter(1, h->GetBinCenter(peak_bin));
    g->SetParameter(2, 0.010);
    g->SetParLimits(0, 0.0, 1e12);
    g->SetParLimits(1, kPreMin, kPreMax);
    g->SetParLimits(2, 0.003, 0.025);

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
    // Full-fit-range integral of residual (data − bg) — captures the entire
    // CB left tail that symmetric ±Nσ windows miss for asymmetric peaks.
    double yield_data_full = 0, yield_data_full_err = 0;
    // Asymmetric window [μ − 5σ, μ + 3σ] — wider on the left to catch the
    // CB tail, tighter on the right where the peak is Gaussian-like.
    double yield_data_m5p3sig = 0, yield_data_m5p3sig_err = 0;
    double chi2 = 0;  int ndf = 0;
    double a0 = 0, a1 = 0, a2 = 0, b_exp = 0, c_exp = 0;
    bool   ok = false;
};

FitRes fitOneHist(TH1D* h, const std::string& base_name,
                  const std::string& canvas_title,
                  const std::string& out_basename,
                  std::vector<TCanvas*>& canvases,
                  std::vector<TCanvas*>& canvases_res,
                  double alpha_ref = -1.0)
{
    FitRes r;
    if (!h || h->GetEntries() == 0) return r;

    const double bw = h->GetBinWidth(1);

    PrefitRes pre = prefitGauss(h);

    auto* fit = new TF1(("fit_" + base_name).c_str(),
                        FitFunction, kFitMin, kFitMax, 10);

    fit->SetParName(0, "N_{CB}");
    fit->SetParName(1, "#mu");
    fit->SetParName(2, "#sigma");
    fit->SetParName(3, "#alpha");
    fit->SetParName(4, "n");
    fit->SetParName(5, "a_{0}");
    fit->SetParName(6, "a_{1}");
    fit->SetParName(7, "a_{2}");
    fit->SetParName(8, "b_{exp}");
    fit->SetParName(9, "c_{exp}");

    // Find peak height inside the apex window for amplitude seed.
    int b_apex = h->FindBin(pre.mu);
    double peak_val = h->GetBinContent(b_apex);
    if (peak_val <= 0) peak_val = h->GetMaximum();

    // α: free [kAlphaMin, kAlphaMax] for the full-panel fit (alpha_ref < 0).
    // For per-slice fits we pass the full-panel α as a tight prior to
    // prevent the tail from running away to compensate bg fluctuations
    // in low-stats slices. Asymmetric range — leans toward larger α
    // (= shorter / less prominent tail) since user prefers smaller tail
    // contribution in slices.
    const double alpha_init = (alpha_ref > 0) ? alpha_ref : kAlphaInit;
    const double alpha_lo   = (alpha_ref > 0)
                              ? std::max(kAlphaMin, alpha_ref - 0.3)
                              : kAlphaMin;
    const double alpha_hi   = (alpha_ref > 0)
                              ? std::min(kAlphaMax, alpha_ref + 1.0)
                              : kAlphaMax;

    fit->SetParameter(0, peak_val);
    fit->SetParameter(1, pre.mu);
    fit->SetParameter(2, kSigInit);
    fit->SetParameter(3, alpha_init);
    fit->SetParameter(4, kFixedN);
    fit->SetParameter(5, 50.0);
    fit->SetParameter(6, 0.0);
    fit->SetParameter(7, 0.0);
    fit->SetParameter(8, 100.0);
    fit->SetParameter(9, 5.0);

    fit->SetParLimits(0, 0.0, 1e9);
    fit->SetParLimits(1, kMuMin, kMuMax);
    fit->SetParLimits(2, kSigMin, kSigMax);
    fit->SetParLimits(3, alpha_lo, alpha_hi);
    fit->FixParameter(4, kFixedN);
    fit->SetParLimits(8, 0.0, 1e9);
    fit->SetParLimits(9, 0.0, 100.0);

    int status = h->Fit(fit, "RQ");
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
    r.chi2    = fit->GetChisquare();
    r.ndf     = fit->GetNDF();

    auto* cbOnly = new TF1(("cb_" + base_name).c_str(),
                           CrystalBallLeft, kFitMin, kFitMax, 5);
    for (int i = 0; i < 5; ++i) cbOnly->SetParameter(i, fit->GetParameter(i));
    cbOnly->SetLineColor(kBlue);
    cbOnly->SetLineWidth(2);
    cbOnly->SetLineStyle(2);

    r.yield_full     = cbOnly->Integral(kFitMin, kFitMax) / bw;
    r.yield_full_err = (r.N_cb > 0) ? r.yield_full * (r.N_cb_err / r.N_cb) : 0.0;
    const double w_lo = std::max(kFitMin, r.mu - 3.0 * r.sigma);
    const double w_hi = std::min(kFitMax, r.mu + 3.0 * r.sigma);
    r.yield_3sig     = cbOnly->Integral(w_lo, w_hi) / bw;
    r.yield_3sig_err = (r.N_cb > 0) ? r.yield_3sig * (r.N_cb_err / r.N_cb) : 0.0;

    auto* bgOnly = new TF1(("bg_" + base_name).c_str(),
        "[0] + [1]*x + [2]*x*x + [3]*exp(-[4]*x)", kFitMin, kFitMax);
    bgOnly->SetParameters(r.a0, r.a1, r.a2, r.b_exp, r.c_exp);
    bgOnly->SetLineColor(kGreen + 2);
    bgOnly->SetLineWidth(2);
    bgOnly->SetLineStyle(3);

    auto* c = new TCanvas(("c_" + base_name).c_str(), canvas_title.c_str(), 900, 700);
    c->SetMargin(0.12, 0.05, 0.12, 0.08);
    h->SetMarkerStyle(20);
    h->SetMarkerSize(0.7);
    h->SetMarkerColor(kBlue);
    h->SetLineColor(kBlue);
    h->SetTitle(canvas_title.c_str());
    h->GetXaxis()->SetRangeUser(kDispMin, kDispMax);

    // Y range: 1.3 × max bin in display window only. Without this auto-Y
    // would be inflated by the π⁰ peak that sits inside the signal hist
    // (CB subtraction doesn't remove real π⁰ — only like-sign combinatorial).
    {
        const int b_lo = h->FindBin(kDispMin);
        const int b_hi = h->FindBin(kDispMax);
        double y_max = 0.0;
        for (int b = std::max(1, b_lo); b <= std::min(h->GetNbinsX(), b_hi); ++b) {
            const double v = h->GetBinContent(b);
            if (v > y_max) y_max = v;
        }
        if (y_max > 0.0) h->SetMaximum(1.3 * y_max);
        h->SetMinimum(0.0);
    }

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
    leg->AddEntry(h,      "signal (all #minus CB)",  "lpe");
    leg->AddEntry(fit,    "CB + P_{2} + exp",        "l");
    leg->AddEntry(cbOnly, "Crystal Ball: #eta",      "l");
    leg->AddEntry(bgOnly, "background: P_{2} + exp", "l");
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
    const std::string per = "plots/output/fit_eta_" + out_basename;
    c->SaveAs((per + ".pdf").c_str());
    c->SaveAs((per + ".png").c_str());
    canvases.push_back(c);

    // Residual plot — signal − bg + ±Nσ true integrals.
    TH1D* h_res = (TH1D*)h->Clone(("res_" + base_name).c_str());
    h_res->SetDirectory(nullptr);

    for (int b = 1; b <= h_res->GetNbinsX(); ++b) {
        const double m = h_res->GetBinCenter(b);
        const double bg_val = bgOnly->Eval(m);
        h_res->SetBinContent(b, h->GetBinContent(b) - bg_val);
    }

    // Clip integration window to the fit range — outside [kFitMin, kFitMax]
    // the bg model is extrapolated and the residual contains real signal
    // from other peaks (π⁰), which would falsely inflate the η yield.
    auto trueIntegral = [&](double m_lo, double m_hi) {
        const double m_lo_c = std::max(m_lo, kFitMin);
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
        auto p = trueIntegral(kFitMin, kFitMax);
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
    h_res->GetXaxis()->SetRangeUser(kDispMin, kDispMax);
    h_res->SetMinimum(0.0);

    // Y range: 1.3 × max bin in display window only.
    {
        const int b_lo = h_res->FindBin(kDispMin);
        const int b_hi = h_res->FindBin(kDispMax);
        double y_max = 0.0;
        for (int b = std::max(1, b_lo); b <= std::min(h_res->GetNbinsX(), b_hi); ++b) {
            const double v = h_res->GetBinContent(b);
            if (v > y_max) y_max = v;
        }
        if (y_max > 0.0) h_res->SetMaximum(1.3 * y_max);
    }

    h_res->Draw("E");
    cbOnly->Draw("same");

    {
        const double y_lo = 0.0;
        const double y_hi = h_res->GetMaximum();
        const int colors[3] = {kGray + 1, kGray + 2, kGray + 3};
        for (int nn = 1; nn <= 3; ++nn) {
            for (int sign = -1; sign <= 1; sign += 2) {
                auto* line = new TLine(r.mu + sign * nn * r.sigma, y_lo,
                                       r.mu + sign * nn * r.sigma, y_hi);
                line->SetLineColor(colors[nn - 1]);
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
    leg2->AddEntry(h_res,     "signal #minus bg", "lpe");
    leg2->AddEntry(cbOnly,    "Crystal Ball: #eta",   "l");
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
    const std::string per_res = "plots/output/fit_eta_residual_" + out_basename;
    c2->SaveAs((per_res + ".pdf").c_str());
    c2->SaveAs((per_res + ".png").c_str());
    canvases_res.push_back(c2);

    return r;
}

void fit_eta(const char* flavour = "rec") {

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
    const std::string pdf_multi     = "plots/output/fit_eta_" + fl + "_all.pdf";
    const std::string pdf_multi_res = "plots/output/fit_eta_" + fl + "_residual_all.pdf";

    const std::string fres_path = "fit_results_eta_" + fl + ".root";
    TFile* fres = TFile::Open(fres_path.c_str(), "RECREATE");
    auto* tres = new TTree("fit_results",
        "Crystal Ball + bg fits per OA slice (signal = all − CB), η region");

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
    float  a0 = 0, a1 = 0, a2 = 0, b_exp = 0, c_exp = 0;

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
    tres->Branch("b_exp",           &b_exp,           "b_exp/F");
    tres->Branch("c_exp",           &c_exp,           "c_exp/F");

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
        a0 = r.a0;  a1 = r.a1;  a2 = r.a2;
        b_exp = r.b_exp;  c_exp = r.c_exp;
        tres->Fill();
    };

    std::vector<TCanvas*> canvases;
    std::vector<TCanvas*> canvases_res;
    canvases.reserve(60);
    canvases_res.reserve(60);

    const int n_slices = static_cast<int>(std::round((kSliceMax - kSliceMin) / kSliceStep));

    // Panel 0: full integrated. We capture α from this fit and pass it as
    // a tight prior to all per-slice fits below — at low stats (high OA
    // slices) bg fluctuations let the CB tail run away unphysically.
    double alpha_ref = -1.0;
    {
        const std::string hname = var_stem + "_full";
        TH1D* h_sig = buildSignal(f_all, f_pp, f_mm, hname, fl);
        if (h_sig) {
            const TString ctitle = TString::Format(
                "M(e^{+}e^{-}#gamma) %s — full OA range [%.1f, %.1f] deg "
                "(signal = all #minus CB, #eta region);"
                "M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Counts",
                fl.c_str(), kSliceMin, kSliceMax);
            FitRes r = fitOneHist(h_sig, hname, ctitle.Data(),
                                  fl + "_full", canvases, canvases_res);
            fillRow(r, hname, 0, kSliceMin, kSliceMax);
            alpha_ref = r.alpha;   // hand off to per-slice fits
            std::cout << "  full: yield(fit) = " << r.yield_full
                      << "   chi2/ndf = "
                      << (r.ndf > 0 ? r.chi2/r.ndf : 0.0)
                      << "   α_ref = " << alpha_ref << "\n";
        }
    }

    for (int i = 0; i < n_slices; ++i) {
        const double lo = kSliceMin + i * kSliceStep;
        const double hi = kSliceMin + (i + 1) * kSliceStep;
        const std::string suff = "oa_" + fmtEdge(lo) + "_" + fmtEdge(hi);
        const std::string hname = var_stem + "_" + suff;

        TH1D* h_sig = buildSignal(f_all, f_pp, f_mm, hname, fl);
        if (!h_sig) { std::cerr << "  missing " << hname << "\n"; continue; }

        const TString ctitle = TString::Format(
            "M(e^{+}e^{-}#gamma) %s, OA #in [%.1f, %.1f] deg (signal, #eta region);"
            "M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Counts",
            fl.c_str(), lo, hi);

        FitRes r = fitOneHist(h_sig, hname, ctitle.Data(),
                              fl + "_slice_" + suff, canvases, canvases_res,
                              alpha_ref);
        fillRow(r, hname, 1 + i, lo, hi);

        std::cout << "  " << hname
                  << ": yield(fit) = " << std::setw(10) << r.yield_full
                  << "   yield(#pm3#sigma) = " << std::setw(10) << r.yield_3sig
                  << "   μ = " << r.mu << "   σ = " << r.sigma
                  << "   χ²/ndf = "
                  << (r.ndf > 0 ? r.chi2/r.ndf : 0.0) << "\n";
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
