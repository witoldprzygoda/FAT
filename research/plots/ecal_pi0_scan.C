// research/plots/ecal_pi0_scan.C — π⁰ peak shift scan vs ECAL variables.
//
// Goal: identify whether the π⁰ peak position μ depends on photon-cluster
// kinematic/quality variables (E_γ, θ_γ, φ_γ mod 60°, ncells, χ²). A μ that
// varies with E_γ is a signature of imperfect ECAL energy-scale calibration
// — the input we need to build a correction map.
//
// Pipeline per scan variable bin:
//   1) read meson_dalitz_nt from output_{epem,epep,emem}_exp.root via
//      TTree::Draw with cut: ecal_quality_pass==1 && var_lo<=v<var_hi
//   2) build all / ++ / -- mass histograms of m_epemg (or m_epemg_cor)
//   3) CB(b) = 2 √(N₊₊(b) · N₋₋(b)),   signal = all − CB
//   4) fit signal with CB + P3 + exp on [0.04, 0.45]
//   5) record μ, σ, α, yield(asym), and scale factor s = (m_PDG/μ)²
//
// Output: ecal_pi0_scan_<var>_<fl>.root with:
//   - TTree "scan_results" (one row per bin)
//   - per-bin signal histograms + canvas snapshots
//   - multipage PDF of fits (plots/output/ecal_pi0_scan_<var>_<fl>.pdf)
//
// Usage (from research/):
//   root -l -b -q plots/ecal_pi0_scan.C                       # E_γ, REC
//   root -l -b -q 'plots/ecal_pi0_scan.C("E","cor")'          # E_γ, COR
//   root -l -b -q 'plots/ecal_pi0_scan.C("theta","rec")'      # θ_γ, REC
//   root -l -b -q 'plots/ecal_pi0_scan.C("phi_local","rec")'  # φ mod 60°
//   root -l -b -q 'plots/ecal_pi0_scan.C("ncells","rec")'     # cluster size
//   root -l -b -q 'plots/ecal_pi0_scan.C("chi2","rec")'       # cluster χ²

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TString.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>

namespace {
    // π⁰ PDG mass — used to derive correction scale s = (m_PDG/μ)².
    constexpr double kPi0PDG = 0.13498;

    // Per-bin m_epemg histogram binning (matches meson_research convention:
    // 5 MeV/bin over [0, 0.8] GeV/c²).
    constexpr int    kHistNBins = 160;
    constexpr double kHistMin   = 0.0;
    constexpr double kHistMax   = 0.8;

    // π⁰ fit window — single, simple range. No OA/E-dependent tweaks here:
    // each bin is fit independently with the same model.
    constexpr double kFitMin = 0.04;
    constexpr double kFitMax = 0.45;

    // Pre-fit Gaussian apex window.
    constexpr double kPreMin = 0.115;
    constexpr double kPreMax = 0.155;

    // Scan-variable definitions: name, TTree::Draw expression,
    // bin edges, axis title.
    struct ScanDef {
        std::string         name;
        std::string         expr;
        std::vector<double> edges;
        std::string         x_title;
    };

    ScanDef makeScan(const std::string& which) {
        // neutr_cluster_energy is in MeV → convert to GeV for plots.
        if (which == "E") {
            return { "E",
                     "neutr_cluster_energy/1000.0",
                     // GeV
                     {0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.70, 1.00, 1.50, 2.50},
                     "E_{#gamma} [GeV]" };
        }
        if (which == "theta") {
            return { "theta",
                     "neutr_cluster_theta",
                     // deg
                     {15.0, 25.0, 35.0, 45.0, 55.0, 65.0, 75.0, 85.0},
                     "#theta_{#gamma} [deg]" };
        }
        if (which == "phi_local") {
            // Fold to local sector [0, 60°). HADES has 6 sectors of 60° each.
            return { "phi_local",
                     "fmod(neutr_cluster_phi+360.0,60.0)",
                     // deg, finer at sector edges where shifts often appear
                     {0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0},
                     "(#phi_{#gamma} mod 60°) [deg]" };
        }
        if (which == "ncells") {
            return { "ncells",
                     "neutr_cluster_ncells",
                     // cluster size
                     {0.5, 1.5, 2.5, 3.5, 4.5, 6.5, 9.5, 30.5},
                     "N_{cells}" };
        }
        if (which == "chi2") {
            return { "chi2",
                     "neutr_chi2",
                     {0.0, 1.0, 2.0, 4.0, 8.0, 16.0, 100.0},
                     "#chi^{2}_{cluster}" };
        }
        return { "", "", {}, "" };
    }
}

// ============================================================================
// Single-sided Crystal Ball (left tail) — NOT normalised. par[0] is plain N.
// ============================================================================
double CBLeft(double *x, double *par) {
    const double m     = x[0];
    const double N     = par[0];
    const double mu    = par[1];
    const double sigma = par[2];
    const double alpha = par[3];
    const double n     = par[4];
    if (sigma <= 0.0 || alpha <= 0.0 || n <= 0.0) return 0.0;
    const double t = (m - mu) / sigma;
    const double aA = std::fabs(alpha);
    if (t > -aA) return N * std::exp(-0.5 * t * t);
    const double A = std::pow(n / aA, n) * std::exp(-0.5 * aA * aA);
    const double B = n / aA - aA;
    return N * A * std::pow(B - t, -n);
}

// CB + P3 + exp.  par[0..4]=CB, par[5..8]=P3, par[9..10]=b·exp(-c·m).
double FitFn(double *x, double *par) {
    const double m  = x[0];
    const double cb = CBLeft(x, par);
    const double p3 = par[5] + par[6]*m + par[7]*m*m + par[8]*m*m*m;
    const double ex = par[9] * std::exp(-par[10]*m);
    return cb + p3 + ex;
}

// ============================================================================
// Pre-fit Gaussian on apex window — seeds μ, σ for the full fit.
// ============================================================================
struct PreRes { double mu = 0.135, sigma = 0.005, N = 0; bool ok = false; };

PreRes prefitApex(TH1D* h) {
    PreRes r;
    const int b_lo = h->FindBin(kPreMin);
    const int b_hi = h->FindBin(kPreMax);
    int    peak_b = b_lo;
    double peak_v = h->GetBinContent(b_lo);
    for (int b = b_lo; b <= b_hi; ++b) {
        if (h->GetBinContent(b) > peak_v) {
            peak_v = h->GetBinContent(b);
            peak_b = b;
        }
    }
    if (peak_v <= 0) {
        r.mu = h->GetBinCenter(peak_b);
        return r;
    }
    auto* g = new TF1("gpre_scan", "gaus", kPreMin, kPreMax);
    g->SetParameter(0, peak_v);
    g->SetParameter(1, h->GetBinCenter(peak_b));
    g->SetParameter(2, 0.005);
    g->SetParLimits(0, 0.0, 1e12);
    g->SetParLimits(1, kPreMin, kPreMax);
    g->SetParLimits(2, 0.001, 0.020);
    if (h->Fit(g, "RQN") == 0) {
        r.N     = g->GetParameter(0);
        r.mu    = g->GetParameter(1);
        r.sigma = g->GetParameter(2);
        r.ok    = true;
    } else {
        r.mu = h->GetBinCenter(peak_b);
    }
    delete g;
    return r;
}

// ============================================================================
// Build N_all, N_++, N_-- mass histograms for one bin via TTree::Draw,
// CB-subtract, return signal (all − CB).
// ============================================================================
TH1D* buildSignal(TTree* t_all, TTree* t_pp, TTree* t_mm,
                  const std::string& mass_var,
                  const std::string& scan_expr,
                  double v_lo, double v_hi,
                  const std::string& tag)
{
    auto draw = [&](TTree* t, const std::string& hname) -> TH1D* {
        // NOTE: histogram MUST live in gDirectory before Draw, otherwise
        // TTree::Draw with `expr>>name` cannot find it and silently creates
        // a new empty default-binned hist, so the named one stays empty.
        // Detach AFTER Draw fills it.
        TH1D* h = new TH1D(hname.c_str(),
                           ";M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Counts",
                           kHistNBins, kHistMin, kHistMax);
        h->Sumw2();
        TString cut = TString::Format("ecal_quality_pass==1 && (%s)>=%g && (%s)<%g",
                                      scan_expr.c_str(), v_lo,
                                      scan_expr.c_str(), v_hi);
        TString draw_cmd = TString::Format("%s>>%s",
                                           mass_var.c_str(), hname.c_str());
        t->Draw(draw_cmd, cut, "goff");
        h->SetDirectory(nullptr);
        return h;
    };

    TH1D* h_all = draw(t_all, "h_all_" + tag);
    TH1D* h_pp  = draw(t_pp,  "h_pp_"  + tag);
    TH1D* h_mm  = draw(t_mm,  "h_mm_"  + tag);

    TH1D* h_cb = (TH1D*)h_all->Clone(("h_cb_" + tag).c_str());
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

    TH1D* h_sig = (TH1D*)h_all->Clone(("h_sig_" + tag).c_str());
    h_sig->SetDirectory(nullptr);
    h_sig->Add(h_cb, -1.0);

    delete h_all; delete h_pp; delete h_mm; delete h_cb;
    return h_sig;
}

// ============================================================================
// Fit one bin — single MIGRAD pass with sensible bounds.
// ============================================================================
struct BinFit {
    double mu = 0,    mu_err = 0;
    double sigma = 0, sigma_err = 0;
    double alpha = 0, n = 0;
    double yield_asym = 0, yield_asym_err = 0;
    double yield_full = 0, yield_full_err = 0;
    double chi2 = 0;  int ndf = 0;
    bool   ok = false;
    // Full parameter vector (length 11) so the canvas can render the
    // EXACT curve from the fit without doing a second blind refit.
    double pars[11] = {0};
};

BinFit fitBin(TH1D* h_sig, const std::string& base) {
    BinFit r;
    if (!h_sig || h_sig->GetEntries() == 0) return r;

    const double bw = h_sig->GetBinWidth(1);

    PreRes pre = prefitApex(h_sig);
    const int b_pre = h_sig->FindBin(pre.mu);
    const double bg_ref = std::max(0.0, h_sig->GetBinContent(b_pre) * 0.05);

    auto* fit = new TF1(("fit_" + base).c_str(),
                        FitFn, kFitMin, kFitMax, 11);

    fit->SetParName(0, "N_{CB}");
    fit->SetParName(1, "#mu");
    fit->SetParName(2, "#sigma");
    fit->SetParName(3, "#alpha");
    fit->SetParName(4, "n");
    fit->SetParName(5, "a_{0}");
    fit->SetParName(6, "a_{1}");
    fit->SetParName(7, "a_{2}");
    fit->SetParName(8, "a_{3}");
    fit->SetParName(9, "b_{exp}");
    fit->SetParName(10,"c_{exp}");

    fit->SetParameter(0, std::max(1.0, pre.N));
    fit->SetParameter(1, pre.mu);
    fit->SetParameter(2, std::max(0.004, std::min(0.020, pre.sigma)));
    fit->SetParameter(3, 1.2);
    fit->SetParameter(4, 5.0);
    fit->SetParameter(5,  bg_ref);
    fit->SetParameter(6,  0.0);
    fit->SetParameter(7,  0.0);
    fit->SetParameter(8,  0.0);
    fit->SetParameter(9,  bg_ref);
    fit->SetParameter(10, 5.0);

    fit->SetParLimits(0, 0.0, 1e12);
    fit->SetParLimits(1, 0.110, 0.155);
    fit->SetParLimits(2, 0.003, 0.025);
    fit->SetParLimits(3, 0.30,  5.0);
    fit->SetParLimits(4, 1.0,   50.0);
    fit->SetParLimits(9, 0.0,   1e12);
    fit->SetParLimits(10, 0.0,  100.0);

    int status = h_sig->Fit(fit, "RQ");
    if (status != 0) {
        // one retry with looser σ start
        fit->SetParameter(2, 0.008);
        status = h_sig->Fit(fit, "RQ");
    }

    r.mu        = fit->GetParameter(1);
    r.mu_err    = fit->GetParError(1);
    r.sigma     = fit->GetParameter(2);
    r.sigma_err = fit->GetParError(2);
    r.alpha     = fit->GetParameter(3);
    r.n         = fit->GetParameter(4);
    r.chi2      = fit->GetChisquare();
    r.ndf       = fit->GetNDF();
    r.ok        = (status == 0);
    for (int i = 0; i < 11; ++i) r.pars[i] = fit->GetParameter(i);

    // CB-only integral for yields.
    auto* cbOnly = new TF1(("cb_" + base).c_str(),
                           CBLeft, kFitMin, kFitMax, 5);
    for (int i = 0; i < 5; ++i) cbOnly->SetParameter(i, fit->GetParameter(i));
    r.yield_full = cbOnly->Integral(kFitMin, kFitMax) / bw;

    // Asymmetric data-driven yield in [μ−5σ, μ+3σ] from signal histogram.
    const double w_lo = std::max(kFitMin, r.mu - 5.0 * r.sigma);
    const double w_hi = std::min(kFitMax, r.mu + 3.0 * r.sigma);
    const int b_lo = h_sig->FindBin(w_lo);
    const int b_hi = h_sig->FindBin(w_hi);
    double err = 0.0;
    r.yield_asym     = h_sig->IntegralAndError(b_lo, b_hi, err);
    r.yield_asym_err = err;

    delete cbOnly;
    delete fit;
    return r;
}

// ============================================================================
// Main entry — orchestrates the scan.
// ============================================================================
void ecal_pi0_scan(const char* scan_var = "E", const char* flavour = "rec") {

    std::string fl = flavour;
    for (auto& c : fl) c = std::tolower(c);
    std::string mass_var;
    if      (fl == "rec") mass_var = "m_epemg";
    else if (fl == "cor") mass_var = "m_epemg_cor";
    else { std::cerr << "Unknown flavour '" << flavour
                     << "' (expected 'rec' or 'cor')\n"; return; }

    ScanDef sd = makeScan(scan_var);
    if (sd.edges.empty()) {
        std::cerr << "Unknown scan variable '" << scan_var << "'\n"
                  << "Choices: E, theta, phi_local, ncells, chi2\n";
        return;
    }

    TFile* f_all = TFile::Open("../output_epem_exp.root", "READ");
    TFile* f_pp  = TFile::Open("../output_epep_exp.root", "READ");
    TFile* f_mm  = TFile::Open("../output_emem_exp.root", "READ");
    if (!f_all || f_all->IsZombie()) { std::cerr << "Cannot open output_epem_exp.root\n"; return; }
    if (!f_pp  || f_pp ->IsZombie()) { std::cerr << "Cannot open output_epep_exp.root\n"; return; }
    if (!f_mm  || f_mm ->IsZombie()) { std::cerr << "Cannot open output_emem_exp.root\n"; return; }

    auto* t_all = (TTree*)f_all->Get("meson_dalitz_nt");
    auto* t_pp  = (TTree*)f_pp ->Get("meson_dalitz_nt");
    auto* t_mm  = (TTree*)f_mm ->Get("meson_dalitz_nt");
    if (!t_all || !t_pp || !t_mm) {
        std::cerr << "meson_dalitz_nt missing in one of the inputs\n"; return;
    }

    gStyle->SetOptStat(0);
    gSystem->mkdir("plots/output", kTRUE);

    const std::string out_root = "ecal_pi0_scan_" + sd.name + "_" + fl + ".root";
    const std::string out_pdf  = "plots/output/ecal_pi0_scan_" + sd.name + "_" + fl + ".pdf";

    TFile* fout = TFile::Open(out_root.c_str(), "RECREATE");
    if (!fout || fout->IsZombie()) { std::cerr << "Cannot create " << out_root << "\n"; return; }

    TTree* tout = new TTree("scan_results",
                            ("π⁰ scan vs " + sd.name + " (" + fl + ")").c_str());
    int    bin_idx = 0;
    float  v_lo = 0, v_hi = 0, v_ctr = 0;
    float  mu = 0, mu_err = 0;
    float  sigma = 0, sigma_err = 0;
    float  alpha = 0, n_par = 0;
    float  yield_asym = 0, yield_asym_err = 0;
    float  yield_full = 0;
    float  chi2 = 0;
    int    ndf = 0;
    int    fit_ok = 0;
    float  s_corr = 0;     // (m_PDG / mu)²
    tout->Branch("bin_idx",        &bin_idx);
    tout->Branch("v_lo",           &v_lo);
    tout->Branch("v_hi",           &v_hi);
    tout->Branch("v_ctr",          &v_ctr);
    tout->Branch("mu",             &mu);
    tout->Branch("mu_err",         &mu_err);
    tout->Branch("sigma",          &sigma);
    tout->Branch("sigma_err",      &sigma_err);
    tout->Branch("alpha",          &alpha);
    tout->Branch("n_par",          &n_par);
    tout->Branch("yield_asym",     &yield_asym);
    tout->Branch("yield_asym_err", &yield_asym_err);
    tout->Branch("yield_full",     &yield_full);
    tout->Branch("chi2",           &chi2);
    tout->Branch("ndf",            &ndf);
    tout->Branch("fit_ok",         &fit_ok);
    tout->Branch("s_corr",         &s_corr);

    std::vector<TCanvas*> canvases;
    canvases.reserve(sd.edges.size() - 1);

    std::cout << "ECAL π⁰ scan: var=" << sd.name
              << "  expr='" << sd.expr << "'"
              << "  mass=" << mass_var
              << "  bins=" << (sd.edges.size() - 1) << "\n";

    for (size_t i = 0; i + 1 < sd.edges.size(); ++i) {
        const double lo = sd.edges[i];
        const double hi = sd.edges[i + 1];
        const std::string tag = TString::Format("%s_b%02zu", sd.name.c_str(), i).Data();

        TH1D* h_sig = buildSignal(t_all, t_pp, t_mm, mass_var, sd.expr,
                                  lo, hi, tag);

        const double n_entries = h_sig->Integral(h_sig->FindBin(kFitMin),
                                                 h_sig->FindBin(kFitMax));
        std::cout << "  bin " << i << " [" << lo << ", " << hi << "): "
                  << "signal_integral≈" << n_entries << "\n";

        BinFit r = fitBin(h_sig, tag);

        bin_idx        = static_cast<int>(i);
        v_lo           = lo;
        v_hi           = hi;
        v_ctr          = 0.5 * (lo + hi);
        mu             = r.mu;
        mu_err         = r.mu_err;
        sigma          = r.sigma;
        sigma_err      = r.sigma_err;
        alpha          = r.alpha;
        n_par          = r.n;
        yield_asym     = r.yield_asym;
        yield_asym_err = r.yield_asym_err;
        yield_full     = r.yield_full;
        chi2           = r.chi2;
        ndf            = r.ndf;
        fit_ok         = r.ok ? 1 : 0;
        s_corr         = (r.mu > 0) ? std::pow(kPi0PDG / r.mu, 2.0) : 0.0;
        tout->Fill();

        // Per-bin canvas — show signal + fit + μ_PDG marker.
        const TString ctitle = TString::Format(
            "%s #in [%g, %g);  #mu=%.4f #pm %.4f;"
            "M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Signal counts",
            sd.x_title.c_str(), lo, hi, r.mu, r.mu_err);
        h_sig->SetTitle(ctitle);
        h_sig->SetMarkerStyle(20);
        h_sig->SetMarkerSize(0.7);
        h_sig->SetMarkerColor(kBlue);
        h_sig->SetLineColor(kBlue);
        h_sig->GetXaxis()->SetRangeUser(0.0, 0.30);

        const TString cname = TString::Format("c_scan_%s_%s_b%02zu",
                                              sd.name.c_str(), fl.c_str(), i);
        auto* c = new TCanvas(cname, cname, 800, 600);
        c->SetMargin(0.13, 0.05, 0.12, 0.08);
        h_sig->Draw("E");

        // Overlay fit using EXACT params from fitBin — do NOT refit here, that
        // would relaunch MIGRAD with default starting values and produce a
        // visually wrong curve while the recorded numbers stay correct.
        auto* fdraw = new TF1(("fdraw_" + tag).c_str(),
                              FitFn, kFitMin, kFitMax, 11);
        for (int p = 0; p < 11; ++p) fdraw->SetParameter(p, r.pars[p]);
        fdraw->SetLineColor(kRed);
        fdraw->SetLineWidth(2);
        fdraw->Draw("SAME");

        auto* lpdg = new TLine(kPi0PDG, 0, kPi0PDG, h_sig->GetMaximum() * 1.05);
        lpdg->SetLineColor(kBlack);
        lpdg->SetLineStyle(2);
        lpdg->Draw();

        auto* leg = new TLegend(0.55, 0.72, 0.94, 0.90);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.030);
        leg->AddEntry(h_sig, "signal (all - CB)",   "lpe");
        leg->AddEntry(fdraw, "CB + P3 + exp",       "l");
        leg->AddEntry(lpdg,  TString::Format("PDG #pi^{0} = %.5f", kPi0PDG).Data(), "l");
        leg->Draw();

        c->Update();
        canvases.push_back(c);

        // Save signal hist into output ROOT
        fout->cd();
        h_sig->Write(("h_sig_" + tag).c_str());
    }

    // Multipage PDF.
    for (size_t i = 0; i < canvases.size(); ++i) {
        auto* c = canvases[i];
        if (canvases.size() == 1)            c->SaveAs(out_pdf.c_str());
        else if (i == 0)                     c->Print((out_pdf + "(").c_str(), "pdf");
        else if (i == canvases.size() - 1)   c->Print((out_pdf + ")").c_str(), "pdf");
        else                                 c->Print(out_pdf.c_str(),         "pdf");
    }

    // -------- Summary plot: μ vs scan variable --------
    const int N = (int)canvases.size();
    std::vector<double> X(N), EX(N), MU(N), EMU(N), SIG(N), SCALE(N);
    for (int i = 0; i < N; ++i) {
        const double lo = sd.edges[i], hi = sd.edges[i + 1];
        X[i]    = 0.5 * (lo + hi);
        EX[i]   = 0.5 * (hi - lo);
        // pull from TTree
    }
    // simpler: re-loop over TTree
    {
        int    bi; float v_lo_, v_hi_, mu_, mu_err_, sig_, scale_;
        tout->SetBranchAddress("bin_idx", &bi);
        tout->SetBranchAddress("v_lo",    &v_lo_);
        tout->SetBranchAddress("v_hi",    &v_hi_);
        tout->SetBranchAddress("mu",      &mu_);
        tout->SetBranchAddress("mu_err",  &mu_err_);
        tout->SetBranchAddress("sigma",   &sig_);
        tout->SetBranchAddress("s_corr",  &scale_);
        for (int i = 0; i < N; ++i) {
            tout->GetEntry(i);
            X[i]     = 0.5 * (v_lo_ + v_hi_);
            EX[i]    = 0.5 * (v_hi_ - v_lo_);
            MU[i]    = mu_;
            EMU[i]   = mu_err_;
            SIG[i]   = sig_;
            SCALE[i] = scale_;
        }
        // restore default branches just in case
        tout->ResetBranchAddresses();
    }

    auto* g_mu = new TGraphErrors(N, X.data(), MU.data(), EX.data(), EMU.data());
    g_mu->SetMarkerStyle(20);
    g_mu->SetMarkerSize(1.0);
    g_mu->SetMarkerColor(kBlue);
    g_mu->SetLineColor(kBlue);
    g_mu->SetTitle(TString::Format(
        "#pi^{0} #mu vs %s (%s);%s;#mu_{CB} [GeV/c^{2}]",
        sd.name.c_str(), fl.c_str(), sd.x_title.c_str()));

    auto* c_mu = new TCanvas(("c_mu_" + sd.name + "_" + fl).c_str(),
                             ("mu vs " + sd.name).c_str(), 1000, 700);
    c_mu->SetMargin(0.13, 0.05, 0.12, 0.08);
    c_mu->SetGrid();
    g_mu->Draw("AP");

    auto* lref = new TLine(sd.edges.front(), kPi0PDG, sd.edges.back(), kPi0PDG);
    lref->SetLineColor(kRed);
    lref->SetLineStyle(2);
    lref->SetLineWidth(2);
    lref->Draw();

    auto* leg_mu = new TLegend(0.55, 0.78, 0.94, 0.90);
    leg_mu->SetBorderSize(0);
    leg_mu->SetFillStyle(0);
    leg_mu->SetTextSize(0.035);
    leg_mu->AddEntry(g_mu, "fit #mu (bars = stat)", "lpe");
    leg_mu->AddEntry(lref, TString::Format("PDG #pi^{0} = %.5f", kPi0PDG).Data(), "l");
    leg_mu->Draw();

    c_mu->Update();
    const std::string mu_base = "plots/output/ecal_pi0_scan_mu_" + sd.name + "_" + fl;
    c_mu->SaveAs((mu_base + ".pdf").c_str());
    c_mu->SaveAs((mu_base + ".png").c_str());

    // -------- Summary plot: scale s vs scan variable --------
    auto* g_s = new TGraphErrors(N, X.data(), SCALE.data(), EX.data(), nullptr);
    g_s->SetMarkerStyle(20);
    g_s->SetMarkerSize(1.0);
    g_s->SetMarkerColor(kMagenta + 1);
    g_s->SetLineColor(kMagenta + 1);
    g_s->SetTitle(TString::Format(
        "ECAL energy-scale s = (m_{PDG}/#mu)^{2} vs %s (%s);%s;s",
        sd.name.c_str(), fl.c_str(), sd.x_title.c_str()));

    auto* c_s = new TCanvas(("c_s_" + sd.name + "_" + fl).c_str(),
                            ("s vs " + sd.name).c_str(), 1000, 700);
    c_s->SetMargin(0.13, 0.05, 0.12, 0.08);
    c_s->SetGrid();
    g_s->Draw("AP");

    auto* lone = new TLine(sd.edges.front(), 1.0, sd.edges.back(), 1.0);
    lone->SetLineColor(kBlack);
    lone->SetLineStyle(2);
    lone->Draw();

    c_s->Update();
    const std::string s_base = "plots/output/ecal_pi0_scan_scale_" + sd.name + "_" + fl;
    c_s->SaveAs((s_base + ".pdf").c_str());
    c_s->SaveAs((s_base + ".png").c_str());

    // -------- Terminal summary table --------
    std::cout << "\n=== ECAL π⁰ scan summary (" << sd.name
              << ", " << fl << ") ===\n";
    std::cout.precision(5);
    std::cout << " bin   v_lo      v_hi      v_ctr      mu        mu_err    sigma     s_corr\n";
    for (int i = 0; i < N; ++i) {
        const double lo = sd.edges[i], hi = sd.edges[i + 1];
        std::cout << "  " << i << "   "
                  << std::fixed
                  << std::setw(8) << lo  << "  "
                  << std::setw(8) << hi  << "  "
                  << std::setw(8) << X[i] << "  "
                  << std::setw(8) << MU[i] << "  "
                  << std::setw(8) << EMU[i] << "  "
                  << std::setw(8) << SIG[i] << "  "
                  << std::setw(8) << SCALE[i] << "\n";
    }

    fout->cd();
    tout->Write();
    fout->Close();

    f_all->Close();
    f_pp ->Close();
    f_mm ->Close();

    std::cout << "\nWrote:\n"
              << "  " << out_root << "  (TTree scan_results, " << N << " bins)\n"
              << "  " << out_pdf  << "  (multipage)\n"
              << "  " << mu_base  << ".{pdf,png}\n"
              << "  " << s_base   << ".{pdf,png}\n";
}
