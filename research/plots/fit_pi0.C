// research/plots/fit_pi0.C — Crystal Ball + P2 + exp fit per OA slice.
//
// Model:
//   F(m) = N · CrystalBall_left(m; μ, σ, α, n)
//        + a₀ + a₁·m + a₂·m²
//        + b · exp(-c·m)
//
// Single full fit per slice. The only assist is a tight Gaussian
// pre-fit on the apex window [0.125, 0.145] used to seed μ.
// Background starts from hard-coded values (the simple original recipe);
// we don't pre-fit sidebands.
//
// Per-slice canvases + fit_results.root TTree as before.
//
// Usage (from research/):
//   root -l -b -q plots/fit_pi0.C

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

    // Default REC; overridden per-flavour in fit_pi0(flavour) below.
    constexpr const char* kVarStemDefault = "m_epemg";

    // Fit range — od 0.02 do końca histogramu (low-OA default).
    // Higher OA slices have rising/cut data at low m, see overrides below.
    constexpr double kFitMin = 0.02;
    constexpr double kFitMax = 0.45;

    // OA-dependent overrides: low-m region behaves differently with OA.
    //  • At OA ≥ 5.0° data curls upward in [0.03, 0.06] → start at 0.03.
    //  • At OA ≥ 7.0° data is cut off below ~0.04         → start at 0.04.
    //  • At OA ≥ 3.8° tighten α so the CB tail doesn't run away
    //    and steal background normalization in the low-m region.
    constexpr double kOA_FitMin003     = 5.0;
    constexpr double kOA_FitMin004     = 7.0;
    constexpr double kOA_AlphaTighten  = 3.8;
    constexpr double kAlphaTightMin    = 1.5;
    constexpr double kAlphaTightMaxMid = 1.65;  // 3.8 ≤ OA < 9.0 — neighbours
                                                // converge to α∈[1.45, 1.71];
                                                // 1.65 kills the local α≈1.7
                                                // dip at OA 6.0–6.2°.
    constexpr double kAlphaTightMaxHi  = 2.0;   // OA ≥ 9.0 — kills the
                                                // spurious α≈2.6 minimum
                                                // seen at OA 9.2–9.4°.

    // Add cubic term a3·m^3 to background only at very high OA, where the
    // 0.05–0.17 region needs an extra inflection that P2 + exp can't make.
    constexpr double kOA_EnableP3      = 9.0;

    // Pre-fit Gaussian range (apex only).
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
// Full model: CB + (P2 or P3) + exp. Always 11 params; par[10] = a3 is
// fixed to 0 below kOA_EnableP3, free above.
//   par[0..4]  = CB (N, μ, σ, α, n)
//   par[5..7]  = P2 (a0, a1, a2)
//   par[8..9]  = exp (b, c)  →  b·exp(-c·m)
//   par[10]    = a3 cubic coefficient (active for OA ≥ kOA_EnableP3)
// ============================================================================
double FitFunction(double *x, double *par) {
    const double m  = x[0];
    const double cb = CrystalBallLeft(x, par);
    const double p3 = par[5] + par[6]*m + par[7]*m*m + par[10]*m*m*m;
    const double ex = par[8] * std::exp(-par[9]*m);
    return cb + p3 + ex;
}

// ============================================================================
// Pre-fit Gauss on the apex window — seeds μ for the full fit.
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
// Single-pass full fit.
// ============================================================================
struct FitRes {
    double mu = 0,    mu_err = 0;
    double sigma = 0, sigma_err = 0;
    double alpha = 0, alpha_err = 0;
    double n = 0,     n_err = 0;
    double N_cb = 0,  N_cb_err = 0;
    double yield_full = 0, yield_full_err = 0;
    double yield_3sig = 0, yield_3sig_err = 0;
    // True data integrals on bg-subtracted histogram (data−bg per bin),
    // summed within μ ± Nσ. Errors via TH1::IntegralAndError (Sumw2-aware).
    double yield_data_1sig = 0, yield_data_1sig_err = 0;
    double yield_data_2sig = 0, yield_data_2sig_err = 0;
    double yield_data_3sig = 0, yield_data_3sig_err = 0;
    double chi2 = 0;  int ndf = 0;
    double a0 = 0, a1 = 0, a2 = 0, a3 = 0, b_exp = 0, c_exp = 0;
    bool   ok = false;
};

FitRes fitOneHist(TH1D* h, const std::string& base_name,
                  const std::string& canvas_title,
                  const std::string& out_basename,
                  double oa_lo,
                  std::vector<TCanvas*>& canvases,
                  std::vector<TCanvas*>& canvases_res)
{
    FitRes r;
    if (!h || h->GetEntries() == 0) return r;

    const double bw = h->GetBinWidth(1);

    // OA-dependent fit-window left edge.
    double fitMin = kFitMin;
    if      (oa_lo >= kOA_FitMin004) fitMin = 0.04;
    else if (oa_lo >= kOA_FitMin003) fitMin = 0.03;

    // OA-dependent α constraint — tighten at higher OA so the CB tail
    // stays well-behaved (no runaway pulling up the low-m region).
    // Mid-OA range gets generous max=3.0; very high OA narrows to 2.0
    // because of an extra spurious minimum near α≈2.6 there.
    double alpha_min = 0.3, alpha_max = 5.0;
    if (oa_lo >= kOA_AlphaTighten) {
        alpha_min = kAlphaTightMin;
        alpha_max = (oa_lo >= kOA_EnableP3) ? kAlphaTightMaxHi
                                            : kAlphaTightMaxMid;
    }

    // Stage 1 — pre-fit Gaussian on apex; seeds μ.
    PrefitRes pre = prefitGauss(h);

    // Stage 2 — single full fit.
    auto* fit = new TF1(("fit_" + base_name).c_str(),
                        FitFunction, fitMin, kFitMax, 11);

    const bool useP3 = (oa_lo >= kOA_EnableP3);

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
    fit->SetParName(10, "a_{3}");

    // Starting values — bg hard-coded (jak w przykładzie),
    // signal: μ z prefit, reszta sensowne defaults.
    fit->SetParameter(0, h->GetMaximum());   // N_CB
    fit->SetParameter(1, pre.mu);            // mu — z prefit
    fit->SetParameter(2, 0.010);             // sigma
    fit->SetParameter(3, 1.5);               // alpha
    fit->SetParameter(4, 3.0);               // n
    fit->SetParameter(5, 100.0);             // a0
    fit->SetParameter(6, -500.0);            // a1
    fit->SetParameter(7, 500.0);             // a2
    fit->SetParameter(8, 1000.0);            // b_exp
    fit->SetParameter(9, 8.0);               // c_exp
    fit->SetParameter(10, 0.0);              // a3

    // Więzy na kształt piku sygnału.
    fit->SetParLimits(0, 0.0, 1e9);
    fit->SetParLimits(1, 0.130, 0.140);
    fit->SetParLimits(2, 0.002, 0.030);
    fit->SetParLimits(3, alpha_min, alpha_max);
    fit->FixParameter(4, 3.0);

    // Eksponenta dodatnia i malejąca.
    fit->SetParLimits(8, 0.0, 1e9);
    fit->SetParLimits(9, 0.0, 100.0);

    // a3 — tylko dla wysokiego OA, w przeciwnym razie zerowane.
    if (useP3) fit->SetParLimits(10, -1e6, 1e6);
    else       fit->FixParameter(10, 0.0);

    // Chi2 fit (Sumw2 weighted hist; "L" Poisson byłaby błędna).
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
    r.a3      = fit->GetParameter(10);
    r.chi2    = fit->GetChisquare();
    r.ndf     = fit->GetNDF();

    // CB-only TF1 dla yieldów + overlayu.
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

    // Bg-only TF1.
    auto* bgOnly = new TF1(("bg_" + base_name).c_str(),
        "[0] + [1]*x + [2]*x*x + [5]*x*x*x + [3]*exp(-[4]*x)",
        fitMin, kFitMax);
    bgOnly->SetParameter(0, r.a0);
    bgOnly->SetParameter(1, r.a1);
    bgOnly->SetParameter(2, r.a2);
    bgOnly->SetParameter(3, r.b_exp);
    bgOnly->SetParameter(4, r.c_exp);
    bgOnly->SetParameter(5, r.a3);
    bgOnly->SetLineColor(kGreen + 2);
    bgOnly->SetLineWidth(2);
    bgOnly->SetLineStyle(3);

    auto* c = new TCanvas(("c_" + base_name).c_str(), canvas_title.c_str(), 900, 700);
    c->SetMargin(0.12, 0.05, 0.12, 0.08);
    h->SetMarkerStyle(20);
    h->SetMarkerSize(0.7);
    h->SetMarkerColor(kBlack);
    h->SetLineColor(kBlack);
    h->SetTitle(canvas_title.c_str());
    h->GetXaxis()->SetRangeUser(0.0, 0.46);

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
    const char* bg_label_full = useP3 ? "CB + P_{3} + exp"        : "CB + P_{2} + exp";
    const char* bg_label_only = useP3 ? "background: P_{3} + exp" : "background: P_{2} + exp";
    leg->AddEntry(h,      "data (sim_genweight)",    "lpe");
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

    // ========================================================================
    // Residual: data minus background model. Plot with CB overlay only,
    // and compute TRUE data integrals (sum of bin contents) inside ±Nσ.
    // ========================================================================
    TH1D* h_res = (TH1D*)h->Clone(("res_" + base_name).c_str());
    h_res->SetDirectory(nullptr);

    // Subtract bg, KEEPING signed values for unbiased integral.
    for (int b = 1; b <= h_res->GetNbinsX(); ++b) {
        const double m = h_res->GetBinCenter(b);
        const double bg_val = bgOnly->Eval(m);
        h_res->SetBinContent(b, h->GetBinContent(b) - bg_val);
        // Bin error stays = data error (bg model has no stat uncertainty).
    }

    // True integrals from the (still signed) residual — unbiased estimator.
    auto trueIntegral = [&](double m_lo, double m_hi) {
        const int b_lo = std::max(1, h_res->FindBin(m_lo));
        const int b_hi = std::min(h_res->GetNbinsX(), h_res->FindBin(m_hi));
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

    // For the plot only: zero negative bins and clear the inherited fit
    // function list (the cloned hist still carries the red full-fit TF1 on
    // its function list — clearing prevents it from being auto-drawn).
    if (auto* funcs = h_res->GetListOfFunctions()) funcs->Clear();
    for (int b = 1; b <= h_res->GetNbinsX(); ++b) {
        if (h_res->GetBinContent(b) < 0.0) h_res->SetBinContent(b, 0.0);
    }

    auto* c2 = new TCanvas(("c_res_" + base_name).c_str(),
                           ("residual " + canvas_title).c_str(), 900, 700);
    c2->SetMargin(0.12, 0.05, 0.12, 0.08);
    h_res->SetMarkerStyle(20);
    h_res->SetMarkerSize(0.7);
    h_res->SetMarkerColor(kBlack);
    h_res->SetLineColor(kBlack);
    h_res->SetTitle(("residual: " + canvas_title).c_str());
    h_res->GetXaxis()->SetRangeUser(0.0, 0.46);
    h_res->SetMinimum(0.0);
    h_res->SetMaximum(h_res->GetMaximum() * 1.15);

    h_res->Draw("E");
    cbOnly->Draw("same");

    // Vertical guide-lines at μ ± {1,2,3}σ — show integration windows.
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
    leg2->AddEntry(h_res,  "data #minus bg",        "lpe");
    leg2->AddEntry(cbOnly, "Crystal Ball: #pi^{0}", "l");
    leg2->Draw();

    auto* tex2 = new TLatex();
    tex2->SetNDC();
    tex2->SetTextSize(0.030);
    tex2->DrawLatex(0.15, 0.86,
        TString::Format("yield(#mu#pm1#sigma) = %.3g #pm %.3g (data#minusbg)",
                        r.yield_data_1sig, r.yield_data_1sig_err));
    tex2->DrawLatex(0.15, 0.82,
        TString::Format("yield(#mu#pm2#sigma) = %.3g #pm %.3g (data#minusbg)",
                        r.yield_data_2sig, r.yield_data_2sig_err));
    tex2->DrawLatex(0.15, 0.78,
        TString::Format("yield(#mu#pm3#sigma) = %.3g #pm %.3g (data#minusbg)",
                        r.yield_data_3sig, r.yield_data_3sig_err));
    tex2->DrawLatex(0.15, 0.72,
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
    else if (fl == "tru") var_stem = "m_epemg_tru";
    else { std::cerr << "Unknown flavour '" << flavour
                     << "' (expected 'rec', 'cor' or 'tru')\n"; return; }
    (void)kVarStemDefault;

    TFile* f = TFile::Open("research_sim.root", "READ");
    if (!f || f->IsZombie()) { std::cerr << "Cannot open research_sim.root\n"; return; }

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    gSystem->mkdir("plots/output", kTRUE);
    const std::string pdf_multi     = "plots/output/fit_pi0_" + fl + "_sim_all.pdf";
    const std::string pdf_multi_res = "plots/output/fit_pi0_" + fl + "_sim_residual_all.pdf";

    const std::string fres_path = "fit_results_" + fl + "_sim.root";
    TFile* fres = TFile::Open(fres_path.c_str(), "RECREATE");
    auto* tres = new TTree("fit_results",
        "Crystal Ball + P2 + exp fits per OA slice");
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
    float  chi2 = 0; int ndf = 0; int ok = 0;
    float  a0 = 0, a1 = 0, a2 = 0, a3 = 0, b_exp = 0, c_exp = 0;

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
    tres->Branch("chi2",            &chi2,            "chi2/F");
    tres->Branch("ndf",             &ndf,             "ndf/I");
    tres->Branch("ok",              &ok,              "ok/I");
    tres->Branch("a0",              &a0,              "a0/F");
    tres->Branch("a1",              &a1,              "a1/F");
    tres->Branch("a2",              &a2,              "a2/F");
    tres->Branch("a3",              &a3,              "a3/F");
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
        chi2 = r.chi2;         ndf = r.ndf;
        ok = r.ok ? 1 : 0;
        a0 = r.a0;  a1 = r.a1;  a2 = r.a2;  a3 = r.a3;
        b_exp = r.b_exp;  c_exp = r.c_exp;
        tres->Fill();
    };

    std::vector<TCanvas*> canvases;
    std::vector<TCanvas*> canvases_res;
    canvases.reserve(60);
    canvases_res.reserve(60);

    const int n_slices = static_cast<int>(std::round((kSliceMax - kSliceMin) / kSliceStep));

    // Panel 0: full integrated spectrum (QA only).
    {
        const std::string hname = var_stem + "_full";
        auto* h_in = (TH1D*)f->Get(hname.c_str());
        if (h_in) {
            auto* h = (TH1D*)h_in->Clone((hname + "_c").c_str());
            h->SetDirectory(nullptr);
            const TString ctitle = TString::Format(
                "M(e^{+}e^{-}#gamma) %s — full OA range [%.1f, %.1f] deg "
                "(QA fit only);"
                "M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Counts (sim_genweight)",
                fl.c_str(), kSliceMin, kSliceMax);
            FitRes r = fitOneHist(h, hname, ctitle.Data(),
                                  fl + "_sim_full",
                                  kSliceMin, canvases, canvases_res);
            fillRow(r, hname, 0, kSliceMin, kSliceMax);
            std::cout << "  full: yield(fit) = " << r.yield_full
                      << " ± " << r.yield_full_err
                      << "   yield(±3σ) = " << r.yield_3sig
                      << "   chi2/ndf = "
                      << (r.ndf > 0 ? r.chi2/r.ndf : 0.0) << "\n";
        }
    }

    for (int i = 0; i < n_slices; ++i) {
        const double lo = kSliceMin + i * kSliceStep;
        const double hi = kSliceMin + (i + 1) * kSliceStep;
        const std::string suff = "oa_" + fmtEdge(lo) + "_" + fmtEdge(hi);
        const std::string hname = var_stem + "_" + suff;

        auto* h_in = (TH1D*)f->Get(hname.c_str());
        if (!h_in) { std::cerr << "  missing " << hname << "\n"; continue; }
        auto* h = (TH1D*)h_in->Clone((hname + "_c").c_str());
        h->SetDirectory(nullptr);

        const TString ctitle = TString::Format(
            "M(e^{+}e^{-}#gamma) %s, OA #in [%.1f, %.1f] deg;"
            "M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Counts (sim_genweight)",
            fl.c_str(), lo, hi);

        FitRes r = fitOneHist(h, hname, ctitle.Data(),
                              fl + "_sim_slice_" + suff,
                              lo, canvases, canvases_res);
        fillRow(r, hname, 1 + i, lo, hi);

        std::cout << "  " << hname
                  << ": yield(fit) = " << std::setw(10) << r.yield_full
                  << "   yield(±3σ) = " << std::setw(10) << r.yield_3sig
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
    f->Close();

    std::cout << "\nResults TTree: research/" << fres_path
              << "  (TTree 'fit_results')\n";
    std::cout << "Multipage PDFs:\n"
              << "  fit:      " << pdf_multi
              << "  (" << canvases.size() << " pages)\n"
              << "  residual: " << pdf_multi_res
              << "  (" << canvases_res.size() << " pages)\n";
    std::cout << "Per-panel PDFs/PNGs:\n"
              << "  plots/output/fit_pi0_" << fl << "_sim_slice_*\n"
              << "  plots/output/fit_pi0_" << fl << "_sim_residual_slice_*\n";
}
