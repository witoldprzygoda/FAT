// research/plots/ecal_pi0_2dscan.C — 2D diagnostic scan μ_π⁰(E_γ, θ_γ).
//
// Builds a (E_γ × θ_γ) grid of m_epemg histograms with full CB subtraction,
// fits each cell with the same CB+P3+exp model, and produces:
//
//   1) heatmap μ_fit(E_γ, θ_γ)        — primary diagnostic
//   2) heatmap s = (m_PDG/μ)²(E_γ,θ_γ) — proposed correction lookup
//   3) heatmap σ_fit(E_γ, θ_γ)        — peak width
//   4) heatmap signal yield(E_γ, θ_γ) — statistics per cell
//   5) multipage PDF of per-cell fits — visual validation
//
// Cell-by-cell statistics are lower than 1D bins, so the cell grid is
// chosen conservatively. Cells with too few signal counts or unstable
// fits are dropped from the heatmaps via fit_ok flag.
//
// Usage (from research/):
//   root -l -b -q plots/ecal_pi0_2dscan.C            # REC, default grid
//   root -l -b -q 'plots/ecal_pi0_2dscan.C("cor")'   # COR

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TString.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>

namespace {
    constexpr double kPi0PDG = 0.13498;

    // m_epemg per-cell histogram binning.
    constexpr int    kHistNBins = 160;
    constexpr double kHistMin   = 0.0;
    constexpr double kHistMax   = 0.8;

    // π⁰ fit window — same as 1D scan.
    constexpr double kFitMin = 0.04;
    constexpr double kFitMax = 0.45;
    constexpr double kPreMin = 0.115;
    constexpr double kPreMax = 0.155;

    // 2D grid — finer than the first pass. Pass-2 fits can lean on
    // neighbour-cell parameters to stabilise low-stat cells, so we don't
    // need to be conservative with the binning.
    // E_γ bins (GeV).
    const std::vector<double> kE_edges = {
        0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00,
        1.25, 1.50, 2.00, 2.50
    };
    // θ_γ bins (deg) — HADES ECAL acceptance ≈ 15–55°.
    const std::vector<double> kTheta_edges = {
        15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0
    };

    // Minimum signal counts in [kFitMin, kFitMax] to attempt a fit.
    constexpr double kMinSignal = 1500.0;
}

// --- Single-sided CB (left tail) + P3 + exp -----------------------------
double CBLeft(double *x, double *par) {
    const double m = x[0], N = par[0], mu = par[1], s = par[2],
                 a = par[3], n = par[4];
    if (s <= 0 || a <= 0 || n <= 0) return 0.0;
    const double t = (m - mu) / s, aA = std::fabs(a);
    if (t > -aA) return N * std::exp(-0.5 * t * t);
    const double A = std::pow(n / aA, n) * std::exp(-0.5 * aA * aA);
    return N * A * std::pow(n / aA - aA - t, -n);
}
double FitFn(double *x, double *par) {
    const double m  = x[0];
    const double cb = CBLeft(x, par);
    const double p3 = par[5] + par[6]*m + par[7]*m*m + par[8]*m*m*m;
    const double ex = par[9] * std::exp(-par[10]*m);
    return cb + p3 + ex;
}

// --- Pre-fit Gaussian seed ----------------------------------------------
struct PreRes { double mu = 0.135, sigma = 0.005, N = 0; };
PreRes prefit(TH1D* h) {
    PreRes r;
    const int b_lo = h->FindBin(kPreMin), b_hi = h->FindBin(kPreMax);
    int    pb = b_lo;
    double pv = h->GetBinContent(b_lo);
    for (int b = b_lo; b <= b_hi; ++b)
        if (h->GetBinContent(b) > pv) { pv = h->GetBinContent(b); pb = b; }
    if (pv <= 0) { r.mu = h->GetBinCenter(pb); return r; }
    auto* g = new TF1("gp_2d", "gaus", kPreMin, kPreMax);
    g->SetParameter(0, pv);
    g->SetParameter(1, h->GetBinCenter(pb));
    g->SetParameter(2, 0.005);
    g->SetParLimits(0, 0.0, 1e12);
    g->SetParLimits(1, kPreMin, kPreMax);
    g->SetParLimits(2, 0.001, 0.020);
    if (h->Fit(g, "RQN") == 0) {
        r.N = g->GetParameter(0);
        r.mu = g->GetParameter(1);
        r.sigma = g->GetParameter(2);
    } else {
        r.mu = h->GetBinCenter(pb);
    }
    delete g;
    return r;
}

// --- Single-channel weighted slice (SIM: no like-sign CB to subtract) -----
// SMASH simulation has no random-charge dilepton background, so the CB
// machinery used on experiment data is unnecessary here. We just project
// the weighted m_epemg distribution into the (E_γ, θ_γ) cell.
TH1D* buildSignal(TTree* t,
                  const std::string& mass_var,
                  double e_lo, double e_hi,
                  double th_lo, double th_hi,
                  const std::string& tag)
{
    TH1D* h = new TH1D(("h_sig_" + tag).c_str(),
                       ";M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Counts (sim_genweight)",
                       kHistNBins, kHistMin, kHistMax);
    h->Sumw2();
    // Weighted draw: TTree::Draw treats the second arg as
    // (cut_expr) × (weight_expr). On sim every fill carries sim_genweight.
    TString cut = TString::Format(
        "(ecal_quality_pass==1"
        " && (neutr_cluster_energy/1000.0)>=%g"
        " && (neutr_cluster_energy/1000.0)<%g"
        " && neutr_cluster_theta>=%g"
        " && neutr_cluster_theta<%g)*sim_genweight",
        e_lo, e_hi, th_lo, th_hi);
    TString cmd = TString::Format("%s>>h_sig_%s",
                                  mass_var.c_str(), tag.c_str());
    t->Draw(cmd, cut, "goff");
    h->SetDirectory(nullptr);
    return h;
}

// --- Single-cell fit -----------------------------------------------------
struct CellFit {
    double mu = 0,    mu_err = 0;
    double sigma = 0, sigma_err = 0;
    double alpha = 0, n_par = 0;
    double yield = 0, yield_err = 0;
    int    n_signal_bins = 0;
    bool   ok = false;
    double pars[11] = {0};
};

// Optional seed parameters for pass-2 fits — supplies σ, α, n from
// already-fitted neighbour cells when MIGRAD struggles unaided.
struct SeedHint {
    double sigma = -1, alpha = -1, n_par = -1;
    bool   freeze_shape = false;   // if true, fix σ, α, n to seed values
};

CellFit fitCell(TH1D* h_sig, const std::string& base,
                const SeedHint& hint = SeedHint{})
{
    CellFit r;
    if (!h_sig || h_sig->GetEntries() == 0) return r;

    const double bw = h_sig->GetBinWidth(1);
    const int b_fit_lo = h_sig->FindBin(kFitMin);
    const int b_fit_hi = h_sig->FindBin(kFitMax);
    const double signal_in_range =
        h_sig->Integral(b_fit_lo, b_fit_hi);
    if (signal_in_range < kMinSignal) return r;

    PreRes pre = prefit(h_sig);
    const int b_pre = h_sig->FindBin(pre.mu);
    const double bg_ref = std::max(0.0, h_sig->GetBinContent(b_pre) * 0.05);

    // Shape parameters: prefer neighbour-derived hints when supplied.
    const double seed_sigma = (hint.sigma > 0)
        ? hint.sigma
        : std::max(0.004, std::min(0.020, pre.sigma));
    const double seed_alpha = (hint.alpha > 0) ? hint.alpha : 1.2;
    const double seed_n     = (hint.n_par > 0) ? hint.n_par : 5.0;

    auto* fit = new TF1(("fit_" + base).c_str(),
                        FitFn, kFitMin, kFitMax, 11);
    fit->SetParameter(0,  std::max(1.0, pre.N));
    fit->SetParameter(1,  pre.mu);
    fit->SetParameter(2,  seed_sigma);
    fit->SetParameter(3,  seed_alpha);
    fit->SetParameter(4,  seed_n);
    fit->SetParameter(5,  bg_ref);
    fit->SetParameter(6,  0.0);
    fit->SetParameter(7,  0.0);
    fit->SetParameter(8,  0.0);
    fit->SetParameter(9,  bg_ref);
    fit->SetParameter(10, 5.0);
    fit->SetParLimits(0, 0.0,    1e12);
    fit->SetParLimits(1, 0.115,  0.150);
    fit->SetParLimits(2, 0.003,  0.025);
    fit->SetParLimits(3, 0.30,   5.0);
    fit->SetParLimits(4, 1.0,    50.0);
    fit->SetParLimits(9, 0.0,    1e12);
    fit->SetParLimits(10, 0.0,   100.0);

    if (hint.freeze_shape) {
        fit->FixParameter(2, seed_sigma);
        fit->FixParameter(3, seed_alpha);
        fit->FixParameter(4, seed_n);
    }

    int status = h_sig->Fit(fit, "RQ");
    if (status != 0) {
        if (!hint.freeze_shape) fit->SetParameter(2, 0.008);
        status = h_sig->Fit(fit, "RQ");
    }
    r.mu        = fit->GetParameter(1);
    r.mu_err    = fit->GetParError(1);
    r.sigma     = fit->GetParameter(2);
    r.sigma_err = fit->GetParError(2);
    r.alpha     = fit->GetParameter(3);
    r.n_par     = fit->GetParameter(4);
    for (int i = 0; i < 11; ++i) r.pars[i] = fit->GetParameter(i);

    // Physically-motivated success criterion. MIGRAD often returns non-zero
    // status when a parameter (frequently α or n) is pinned at a bound,
    // even though the peak-position fit is fine. Accept the fit if:
    //   - μ is well inside its allowed interval (not pinned)
    //   - σ is not at the lower bound (would mean degenerate peak)
    //   - mu_err is small enough that the Hessian is invertible
    const bool mu_in_range    = (r.mu > 0.117 && r.mu < 0.148);
    const bool sigma_sensible = (r.sigma > 0.0035 && r.sigma < 0.0245);
    const bool err_reasonable = (r.mu_err > 0.0 && r.mu_err < 0.005);
    r.ok = mu_in_range && sigma_sensible && err_reasonable;
    (void)bw;

    // Asymmetric data-driven yield in [μ−5σ, μ+3σ] from signal histogram.
    const double w_lo = std::max(kFitMin, r.mu - 5.0 * r.sigma);
    const double w_hi = std::min(kFitMax, r.mu + 3.0 * r.sigma);
    const int b_lo = h_sig->FindBin(w_lo);
    const int b_hi = h_sig->FindBin(w_hi);
    double err = 0.0;
    r.yield     = h_sig->IntegralAndError(b_lo, b_hi, err);
    r.yield_err = err;
    r.n_signal_bins = b_hi - b_lo + 1;

    delete fit;
    return r;
}

// =========================================================================
void ecal_pi0_2dscan(const char* flavour = "rec") {

    std::string fl = flavour;
    for (auto& c : fl) c = std::tolower(c);
    std::string mass_var;
    // SIM has 3 mass flavours: REC (raw HADES tracking), COR (energy-loss
    // corrected), TRU (MC truth, branch m_epemg_sim).
    if      (fl == "rec") mass_var = "m_epemg";
    else if (fl == "cor") mass_var = "m_epemg_cor";
    else if (fl == "tru") mass_var = "m_epemg_sim";
    else { std::cerr << "Unknown flavour '" << flavour << "'\n"; return; }

    // Single sim channel — no like-sign combinatorial.
    TFile* f_all = TFile::Open("../output_epem_sim.root", "READ");
    if (!f_all || f_all->IsZombie()) { std::cerr << "Cannot open output_epem_sim.root\n"; return; }
    auto* t_all = (TTree*)f_all->Get("meson_dalitz_nt");
    if (!t_all) {
        std::cerr << "meson_dalitz_nt missing in input\n"; return;
    }

    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);
    gSystem->mkdir("plots/output", kTRUE);

    const int nE  = (int)kE_edges.size() - 1;
    const int nTh = (int)kTheta_edges.size() - 1;

    // Heatmap histograms — uniform-bin for cosmetics, but we override
    // the bin edges so they align with the variable-width grid.
    auto make2D = [&](const char* name, const char* title) {
        TH2D* h = new TH2D(name, title,
                           nE,  kE_edges.data(),
                           nTh, kTheta_edges.data());
        h->SetDirectory(nullptr);
        h->GetXaxis()->SetTitle("E_{#gamma} [GeV]");
        h->GetYaxis()->SetTitle("#theta_{#gamma} [deg]");
        return h;
    };
    TH2D* h_mu    = make2D(("h_mu_"    + fl).c_str(),
        TString::Format("#mu_{#pi^{0}} (%s);E_{#gamma} [GeV];#theta_{#gamma} [deg];#mu [GeV/c^{2}]", fl.c_str()));
    TH2D* h_mu_err= make2D(("h_muerr_" + fl).c_str(),
        TString::Format("#mu_{err} (%s);E_{#gamma} [GeV];#theta_{#gamma} [deg];#mu err", fl.c_str()));
    TH2D* h_sig   = make2D(("h_sigma_" + fl).c_str(),
        TString::Format("#sigma_{#pi^{0}} (%s);E_{#gamma} [GeV];#theta_{#gamma} [deg];#sigma [GeV/c^{2}]", fl.c_str()));
    TH2D* h_s     = make2D(("h_s_"     + fl).c_str(),
        TString::Format("s = (m_{PDG}/#mu)^{2} (%s);E_{#gamma} [GeV];#theta_{#gamma} [deg];s", fl.c_str()));
    TH2D* h_yield = make2D(("h_yield_" + fl).c_str(),
        TString::Format("Signal yield [#mu#minus5#sigma, #mu+3#sigma] (%s);E_{#gamma} [GeV];#theta_{#gamma} [deg];yield", fl.c_str()));

    const std::string out_root = "ecal_pi0_2dscan_" + fl + "_sim.root";
    const std::string out_pdf  = "plots/output/ecal_pi0_2dscan_" + fl + "_sim.pdf";
    TFile* fout = TFile::Open(out_root.c_str(), "RECREATE");

    // TTree of cell-level results.
    TTree* tout = new TTree("scan2d_results",
                            ("π⁰ 2D scan E×θ (" + fl + ")").c_str());
    int   eb = 0, tb = 0;
    float e_lo = 0, e_hi = 0, th_lo = 0, th_hi = 0;
    float mu = 0, mu_err = 0, sigma = 0, yield = 0, yield_err = 0, s_corr = 0;
    int   fit_ok = 0;
    tout->Branch("eb",          &eb);
    tout->Branch("tb",          &tb);
    tout->Branch("e_lo",        &e_lo);
    tout->Branch("e_hi",        &e_hi);
    tout->Branch("th_lo",       &th_lo);
    tout->Branch("th_hi",       &th_hi);
    tout->Branch("mu",          &mu);
    tout->Branch("mu_err",      &mu_err);
    tout->Branch("sigma",       &sigma);
    tout->Branch("yield",       &yield);
    tout->Branch("yield_err",   &yield_err);
    tout->Branch("s_corr",      &s_corr);
    tout->Branch("fit_ok",      &fit_ok);

    std::cout << "ECAL π⁰ 2D scan: E_γ × θ_γ  flavour=" << fl
              << "  cells=" << (nE * nTh) << "\n";

    std::vector<TCanvas*> canvases;
    canvases.reserve(nE * nTh);

    // -------- Pass 1: build per-cell signal hists and run unaided fits ---
    std::vector<std::vector<TH1D*>>     cell_hist(nE, std::vector<TH1D*>(nTh, nullptr));
    std::vector<std::vector<CellFit>>   cell_fit (nE, std::vector<CellFit>(nTh));
    std::vector<std::vector<double>>    cell_int (nE, std::vector<double>(nTh, 0.0));

    std::cout << "\n--- Pass 1: independent fits ---\n";
    for (int ei = 0; ei < nE; ++ei) {
        for (int ti = 0; ti < nTh; ++ti) {
            const double el = kE_edges[ei],     eh = kE_edges[ei + 1];
            const double tl = kTheta_edges[ti], th = kTheta_edges[ti + 1];
            const std::string tag =
                TString::Format("e%02d_t%02d", ei, ti).Data();

            TH1D* h_signal = buildSignal(t_all, mass_var,
                                         el, eh, tl, th, tag);
            const double sig_int = h_signal->Integral(
                h_signal->FindBin(kFitMin), h_signal->FindBin(kFitMax));
            cell_hist[ei][ti] = h_signal;
            cell_int[ei][ti]  = sig_int;
            cell_fit[ei][ti]  = fitCell(h_signal, tag);
        }
    }

    // -------- Pass 2: refit failed cells using neighbour-averaged shape --
    auto neighbourSeed = [&](int ei, int ti) {
        SeedHint h;
        double sum_s = 0, sum_a = 0, sum_n = 0;
        int cnt = 0;
        for (int de = -1; de <= 1; ++de) {
            for (int dt = -1; dt <= 1; ++dt) {
                if (de == 0 && dt == 0) continue;
                const int ne = ei + de, nt = ti + dt;
                if (ne < 0 || ne >= nE || nt < 0 || nt >= nTh) continue;
                const CellFit& cf = cell_fit[ne][nt];
                if (!cf.ok) continue;
                sum_s += cf.sigma;
                sum_a += cf.alpha;
                sum_n += cf.n_par;
                ++cnt;
            }
        }
        if (cnt > 0) {
            h.sigma = sum_s / cnt;
            h.alpha = sum_a / cnt;
            h.n_par = sum_n / cnt;
        }
        return h;
    };

    int n_recovered = 0;
    std::cout << "\n--- Pass 2: neighbour-seeded refits for failed cells ---\n";
    for (int ei = 0; ei < nE; ++ei) {
        for (int ti = 0; ti < nTh; ++ti) {
            if (cell_fit[ei][ti].ok) continue;
            if (cell_int[ei][ti] < kMinSignal) continue;

            SeedHint hint = neighbourSeed(ei, ti);
            if (hint.sigma <= 0) continue;   // no neighbours available

            const std::string tag =
                TString::Format("e%02d_t%02d_p2", ei, ti).Data();

            // First try with neighbour-seeded but free shape parameters.
            CellFit r2 = fitCell(cell_hist[ei][ti], tag, hint);
            // If still failing, freeze σ, α, n to neighbour values.
            if (!r2.ok) {
                hint.freeze_shape = true;
                const std::string tag_f =
                    TString::Format("e%02d_t%02d_p2f", ei, ti).Data();
                r2 = fitCell(cell_hist[ei][ti], tag_f, hint);
            }
            if (r2.ok) {
                cell_fit[ei][ti] = r2;
                ++n_recovered;
            }
        }
    }
    std::cout << "Pass-2 recovered " << n_recovered << " cells.\n";

    // -------- Final reporting + heatmaps + canvases ----------------------
    std::cout << "\n--- Final per-cell results ---\n";
    for (int ei = 0; ei < nE; ++ei) {
        for (int ti = 0; ti < nTh; ++ti) {
            const double el = kE_edges[ei],     eh = kE_edges[ei + 1];
            const double tl = kTheta_edges[ti], th = kTheta_edges[ti + 1];
            const CellFit& r = cell_fit[ei][ti];
            const double sig_int = cell_int[ei][ti];
            TH1D* h_signal = cell_hist[ei][ti];

            std::cout << "  E[" << el << "," << eh << ") × θ[" << tl << "," << th << ")"
                      << "  sig=" << std::fixed << std::setprecision(0) << sig_int
                      << "  μ=" << std::setprecision(4) << r.mu
                      << " ± " << r.mu_err
                      << "  σ=" << r.sigma
                      << "  ok=" << (r.ok ? "Y" : "N") << "\n";

            eb = ei; tb = ti;
            e_lo = el; e_hi = eh; th_lo = tl; th_hi = th;
            mu = r.mu; mu_err = r.mu_err; sigma = r.sigma;
            yield = r.yield; yield_err = r.yield_err;
            s_corr = (r.mu > 0) ? std::pow(kPi0PDG / r.mu, 2.0) : 0.0;
            fit_ok = r.ok ? 1 : 0;
            tout->Fill();

            if (r.ok && r.mu > 0 && sig_int > kMinSignal) {
                const int gx = ei + 1, gy = ti + 1;
                h_mu   ->SetBinContent(gx, gy, r.mu);
                h_mu_err->SetBinContent(gx, gy, r.mu_err);
                h_sig  ->SetBinContent(gx, gy, r.sigma);
                h_s    ->SetBinContent(gx, gy, s_corr);
                h_yield->SetBinContent(gx, gy, r.yield);
            }

            const TString ctitle = TString::Format(
                "E_{#gamma}#in[%g,%g) GeV, #theta_{#gamma}#in[%g,%g)°  #mu=%.4f#pm%.4f;"
                "M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Signal counts",
                el, eh, tl, th, r.mu, r.mu_err);
            h_signal->SetTitle(ctitle);
            h_signal->SetMarkerStyle(20);
            h_signal->SetMarkerSize(0.7);
            h_signal->SetMarkerColor(kBlue);
            h_signal->SetLineColor(kBlue);
            h_signal->GetXaxis()->SetRangeUser(0.0, 0.30);

            const std::string tag_disp =
                TString::Format("e%02d_t%02d", ei, ti).Data();
            const TString cname = TString::Format(
                "c_2dscan_%s_%s", fl.c_str(), tag_disp.c_str());
            auto* c = new TCanvas(cname, cname, 800, 600);
            c->SetMargin(0.13, 0.05, 0.12, 0.08);
            h_signal->Draw("E");

            if (r.ok) {
                auto* fdraw = new TF1(("fdraw_" + tag_disp).c_str(),
                                      FitFn, kFitMin, kFitMax, 11);
                for (int p = 0; p < 11; ++p) fdraw->SetParameter(p, r.pars[p]);
                fdraw->SetLineColor(kRed);
                fdraw->SetLineWidth(2);
                fdraw->Draw("SAME");
            }
            auto* lpdg = new TLine(kPi0PDG, 0, kPi0PDG,
                                   h_signal->GetMaximum() * 1.05);
            lpdg->SetLineColor(kBlack);
            lpdg->SetLineStyle(2);
            lpdg->Draw();
            c->Update();
            canvases.push_back(c);

            fout->cd();
            h_signal->Write(("h_sig_" + tag_disp).c_str());
        }
    }

    // Multipage PDF.
    for (size_t i = 0; i < canvases.size(); ++i) {
        auto* c = canvases[i];
        if (canvases.size() == 1)            c->SaveAs(out_pdf.c_str());
        else if (i == 0)                     c->Print((out_pdf + "(").c_str(), "pdf");
        else if (i == canvases.size() - 1)   c->Print((out_pdf + ")").c_str(), "pdf");
        else                                 c->Print(out_pdf.c_str(),         "pdf");
    }

    // Heatmap canvases.
    auto drawHeat = [&](TH2D* h, const std::string& base, const char* opt) {
        TString cname = TString::Format("c_%s", h->GetName());
        auto* c = new TCanvas(cname, cname, 1100, 800);
        c->SetMargin(0.13, 0.15, 0.12, 0.10);
        c->SetGrid();
        h->SetMarkerSize(1.4);
        h->Draw(opt);
        c->Update();
        c->SaveAs(("plots/output/" + base + "_" + fl + ".pdf").c_str());
        c->SaveAs(("plots/output/" + base + "_" + fl + ".png").c_str());
    };

    h_mu   ->SetMinimum(0.130);
    h_mu   ->SetMaximum(0.150);
    h_s    ->SetMinimum(0.85);
    h_s    ->SetMaximum(1.05);
    h_sig  ->SetMinimum(0.005);
    h_sig  ->SetMaximum(0.020);

    drawHeat(h_mu,    "ecal_pi0_2dscan_sim_mu",    "COLZ TEXT45");
    drawHeat(h_s,     "ecal_pi0_2dscan_sim_scale", "COLZ TEXT45");
    drawHeat(h_sig,   "ecal_pi0_2dscan_sim_sigma", "COLZ TEXT45");
    drawHeat(h_yield, "ecal_pi0_2dscan_sim_yield", "COLZ TEXT");

    // Terminal summary table.
    std::cout << "\n=== ECAL π⁰ 2D scan summary (" << fl << ") ===\n";
    std::cout << "  rows = θ bins,  cols = E_γ bins\n      ";
    for (int ei = 0; ei < nE; ++ei)
        std::cout << "  E[" << kE_edges[ei] << "," << kE_edges[ei + 1] << ")";
    std::cout << "\n";
    for (int ti = 0; ti < nTh; ++ti) {
        std::cout << "θ[" << kTheta_edges[ti]
                  << "," << kTheta_edges[ti + 1] << "):  ";
        for (int ei = 0; ei < nE; ++ei) {
            const double m = h_mu->GetBinContent(ei + 1, ti + 1);
            if (m > 0)
                std::cout << "  " << std::fixed << std::setprecision(4) << m;
            else
                std::cout << "    -   ";
        }
        std::cout << "\n";
    }

    fout->cd();
    h_mu->Write();  h_mu_err->Write();  h_sig->Write();
    h_s->Write();   h_yield->Write();
    tout->Write();
    fout->Close();

    f_all->Close();

    std::cout << "\nWrote:\n"
              << "  " << out_root << "  (TTree + heatmaps)\n"
              << "  " << out_pdf  << "  (multipage cells)\n"
              << "  plots/output/ecal_pi0_2dscan_sim_{mu,scale,sigma,yield}_" << fl << ".{pdf,png}\n";
}
