// research/plots/ecal_pi0_2dscan_corr.C — π⁰ peak scan in (E_γ, θ_γ) cells
// AFTER applying the ECAL energy-scale correction.
//
// Mirrors ecal_pi0_2dscan.C cell by cell, but instead of filling
// histograms with `m_epemg` from the ntuple, it loops events, looks up
// the correction scale s(E_γ, θ_γ) bilinearly from the calibration map
// produced by ecal_pi0_2dscan.C, and fills with
//
//     m_corr = √(m_ee² + 2·s·E_γ·D)        D = gamma_D from ntuple
//
// The output heatmap μ_corr(E_γ, θ_γ) is the validation: if the
// calibration is self-consistent, μ_corr should sit at PDG 0.13498
// across the entire grid (modulo statistical jitter). Any residual
// pattern indicates either:
//   - the 1D s(E,θ) lookup is too coarse (need finer grid), or
//   - the correction model itself misses a hidden axis (φ, ncells, …).
//
// Usage (from research/):
//   root -l -b -q plots/ecal_pi0_2dscan_corr.C            # REC
//   root -l -b -q 'plots/ecal_pi0_2dscan_corr.C("cor")'   # COR

#include "ecal_apply_lookup.h"

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

    constexpr int    kHistNBins = 160;
    constexpr double kHistMin   = 0.0;
    constexpr double kHistMax   = 0.8;

    constexpr double kFitMin = 0.04;
    constexpr double kFitMax = 0.45;
    constexpr double kPreMin = 0.115;
    constexpr double kPreMax = 0.155;

    const std::vector<double> kE_edges = {
        0.20, 0.30, 0.40, 0.55, 0.75, 1.00, 1.50, 2.50
    };
    const std::vector<double> kTheta_edges = {
        15.0, 25.0, 35.0, 45.0, 55.0
    };

    constexpr double kMinSignal = 5000.0;

    int findBin(const std::vector<double>& edges, double v) {
        if (v < edges.front() || v >= edges.back()) return -1;
        for (size_t i = 0; i + 1 < edges.size(); ++i)
            if (v >= edges[i] && v < edges[i + 1]) return (int)i;
        return -1;
    }
}

// --- Single-sided CB + P3 + exp -----------------------------------------
double CBLeftCorr(double *x, double *par) {
    const double m = x[0], N = par[0], mu = par[1], s = par[2],
                 a = par[3], n = par[4];
    if (s <= 0 || a <= 0 || n <= 0) return 0.0;
    const double t = (m - mu) / s, aA = std::fabs(a);
    if (t > -aA) return N * std::exp(-0.5 * t * t);
    const double A = std::pow(n / aA, n) * std::exp(-0.5 * aA * aA);
    return N * A * std::pow(n / aA - aA - t, -n);
}
double FitFnCorr(double *x, double *par) {
    const double m  = x[0];
    const double cb = CBLeftCorr(x, par);
    const double p3 = par[5] + par[6]*m + par[7]*m*m + par[8]*m*m*m;
    const double ex = par[9] * std::exp(-par[10]*m);
    return cb + p3 + ex;
}

struct PreResC { double mu = 0.135, sigma = 0.005, N = 0; };
PreResC prefitCorr(TH1D* h) {
    PreResC r;
    const int b_lo = h->FindBin(kPreMin), b_hi = h->FindBin(kPreMax);
    int    pb = b_lo;
    double pv = h->GetBinContent(b_lo);
    for (int b = b_lo; b <= b_hi; ++b)
        if (h->GetBinContent(b) > pv) { pv = h->GetBinContent(b); pb = b; }
    if (pv <= 0) { r.mu = h->GetBinCenter(pb); return r; }
    auto* g = new TF1("gpc", "gaus", kPreMin, kPreMax);
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

struct CellFitC {
    double mu = 0, mu_err = 0, sigma = 0;
    bool   ok = false;
    double pars[11] = {0};
};

CellFitC fitCellCorr(TH1D* h, const std::string& base) {
    CellFitC r;
    if (!h || h->GetEntries() == 0) return r;
    const int b_fit_lo = h->FindBin(kFitMin);
    const int b_fit_hi = h->FindBin(kFitMax);
    if (h->Integral(b_fit_lo, b_fit_hi) < kMinSignal) return r;

    PreResC pre = prefitCorr(h);
    // For the CORRECTED mass we know μ should be near PDG. Seed there
    // and use tighter bounds — prevents MIGRAD from settling on a
    // bg+sideband local minimum (observed in some cells where the prefit
    // happened to pick a high-stat bin off-peak).
    const double mu_seed = kPi0PDG;
    const int b_pre = h->FindBin(mu_seed);
    const double bg_ref = std::max(0.0, h->GetBinContent(b_pre) * 0.05);

    auto* fit = new TF1(("fitc_" + base).c_str(),
                        FitFnCorr, kFitMin, kFitMax, 11);
    fit->SetParameter(0,  std::max(1.0, pre.N));
    fit->SetParameter(1,  mu_seed);
    fit->SetParameter(2,  std::max(0.005, std::min(0.018, pre.sigma)));
    fit->SetParameter(3,  1.2);
    fit->SetParameter(4,  5.0);
    fit->SetParameter(5,  bg_ref);
    fit->SetParameter(6,  0.0);
    fit->SetParameter(7,  0.0);
    fit->SetParameter(8,  0.0);
    fit->SetParameter(9,  bg_ref);
    fit->SetParameter(10, 5.0);
    fit->SetParLimits(0, 0.0,    1e12);
    fit->SetParLimits(1, 0.125,  0.148);
    fit->SetParLimits(2, 0.003,  0.025);
    fit->SetParLimits(3, 0.30,   5.0);
    fit->SetParLimits(4, 1.0,    50.0);
    fit->SetParLimits(9, 0.0,    1e12);
    fit->SetParLimits(10, 0.0,   100.0);

    int status = h->Fit(fit, "RQ");
    if (status != 0) {
        fit->SetParameter(2, 0.008);
        status = h->Fit(fit, "RQ");
    }
    r.mu     = fit->GetParameter(1);
    r.mu_err = fit->GetParError(1);
    r.sigma  = fit->GetParameter(2);
    for (int i = 0; i < 11; ++i) r.pars[i] = fit->GetParameter(i);

    const bool mu_in_range    = (r.mu > 0.112 && r.mu < 0.153);
    const bool sigma_sensible = (r.sigma > 0.0035 && r.sigma < 0.0245);
    const bool err_reasonable = (r.mu_err > 0.0 && r.mu_err < 0.003);
    r.ok = mu_in_range && sigma_sensible && err_reasonable;

    delete fit;
    return r;
}

// =========================================================================
// Loop the single sim channel's meson_dalitz_nt and fill per-cell mass
// histos using the ECAL-corrected mass. SMASH per-event weight
// (sim_genweight) is applied as the fill weight.
// =========================================================================
void fillChannel(const std::string& fpath,
                 const std::string& /*mass_var*/,
                 const EcalLookup& look,
                 std::vector<std::vector<TH1D*>>& cells,
                 const std::string& tag_prefix)
{
    TFile* f = TFile::Open(fpath.c_str(), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Cannot open " << fpath << "\n"; return;
    }
    auto* t = (TTree*)f->Get("meson_dalitz_nt");
    if (!t) { std::cerr << "meson_dalitz_nt missing in " << fpath << "\n"; return; }

    // Activate only the branches we read — speeds up the loop dramatically.
    // Closed-form correction is used here, so gamma_D is NOT required.
    t->SetBranchStatus("*", 0);
    for (const char* b : {
            "m_ee",
            "m_epemg",
            "ecal_quality_pass",
            "neutr_cluster_energy",
            "neutr_cluster_theta",
            "sim_genweight"
         }) t->SetBranchStatus(b, 1);

    float m_ee=0, m_eg=0, ecal_q=0, ne_E=0, ne_th=0, w=1.0f;
    t->SetBranchAddress("m_ee",                 &m_ee);
    t->SetBranchAddress("m_epemg",              &m_eg);
    t->SetBranchAddress("ecal_quality_pass",    &ecal_q);
    t->SetBranchAddress("neutr_cluster_energy", &ne_E);
    t->SetBranchAddress("neutr_cluster_theta",  &ne_th);
    t->SetBranchAddress("sim_genweight",        &w);

    const Long64_t N = t->GetEntries();
    std::cout << "  " << fpath << " (" << tag_prefix << "): " << N << " entries\n";

    for (Long64_t ev = 0; ev < N; ++ev) {
        t->GetEntry(ev);
        if (ecal_q != 1.0f) continue;

        const double E_GeV   = ne_E / 1000.0;
        const double theta   = ne_th;
        const int    ei      = findBin(kE_edges,     E_GeV);
        const int    ti      = findBin(kTheta_edges, theta);
        if (ei < 0 || ti < 0) continue;

        const double m_corr = look.m_corrected_closed(m_ee, m_eg, E_GeV, theta);
        if (m_corr <= 0) continue;

        cells[ei][ti]->Fill(m_corr, w);
    }

    f->Close();
}

// =========================================================================
void ecal_pi0_2dscan_corr(const char* flavour = "rec") {

    std::string fl = flavour;
    for (auto& c : fl) c = std::tolower(c);
    std::string mass_var;
    // SIM has REC/COR/TRU; the closed-form correction here uses m_ee (REC),
    // which is fine for the REC flavour. Extending to COR/TRU would need
    // m_ee_cor / m_ee_sim in the formula and a per-flavour lookup map.
    if      (fl == "rec") mass_var = "m_epemg";
    else if (fl == "cor") mass_var = "m_epemg_cor";
    else { std::cerr << "Unknown flavour '" << flavour << "'\n"; return; }

    // Load calibration map from the previous (uncorrected) sim scan.
    EcalLookup look;
    const std::string map_path = "ecal_pi0_2dscan_" + fl + "_sim.root";
    const std::string map_hist = "h_s_" + fl;
    if (!look.load(map_path, map_hist)) {
        std::cerr << "Run ecal_pi0_2dscan.C first to produce " << map_path << "\n";
        return;
    }

    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);
    gSystem->mkdir("plots/output", kTRUE);

    const int nE  = (int)kE_edges.size() - 1;
    const int nTh = (int)kTheta_edges.size() - 1;

    // Allocate per-channel cell hists.
    auto makeCells = [&](const std::string& chan) {
        std::vector<std::vector<TH1D*>> cells(nE,
            std::vector<TH1D*>(nTh, nullptr));
        for (int ei = 0; ei < nE; ++ei)
            for (int ti = 0; ti < nTh; ++ti) {
                const TString hn = TString::Format(
                    "h_%s_corr_%s_e%02d_t%02d", chan.c_str(), fl.c_str(), ei, ti);
                cells[ei][ti] = new TH1D(hn,
                    ";M_{e^{+}e^{-}#gamma}^{corr} [GeV/c^{2}];Counts",
                    kHistNBins, kHistMin, kHistMax);
                cells[ei][ti]->Sumw2();
                cells[ei][ti]->SetDirectory(nullptr);
            }
        return cells;
    };

    auto cells_all = makeCells("all");

    std::cout << "ECAL π⁰ corrected 2D scan (SIM): flavour=" << fl
              << "  cells=" << (nE * nTh) << "\n";
    std::cout << "Map: " << map_path << " (" << map_hist << ")\n";

    fillChannel("../output_epem_sim.root", mass_var, look, cells_all, "sim");

    // SIM has no like-sign CB subtraction — the corrected histogram IS
    // the signal. Helper kept as a thin clone wrapper to mirror the
    // shape of the experiment-side macro and keep the rest of the
    // pipeline identical.
    auto buildSignalCorr = [&](TH1D* h_all,
                               const std::string& tag) -> TH1D* {
        TH1D* h_sig = (TH1D*)h_all->Clone(("h_sig_corr_" + tag).c_str());
        h_sig->SetDirectory(nullptr);
        return h_sig;
    };

    auto make2D = [&](const char* name, const char* title) {
        TH2D* h = new TH2D(name, title,
                           nE,  kE_edges.data(),
                           nTh, kTheta_edges.data());
        h->SetDirectory(nullptr);
        h->GetXaxis()->SetTitle("E_{#gamma} [GeV]");
        h->GetYaxis()->SetTitle("#theta_{#gamma} [deg]");
        return h;
    };
    TH2D* h_mu_corr = make2D(("h_mu_corr_" + fl).c_str(),
        TString::Format("#mu_{#pi^{0}} after correction (%s);"
                        "E_{#gamma} [GeV];#theta_{#gamma} [deg];#mu_{corr} [GeV/c^{2}]",
                        fl.c_str()));
    TH2D* h_dmu = make2D(("h_dmu_" + fl).c_str(),
        TString::Format("#mu_{corr} #minus m_{PDG} (%s);"
                        "E_{#gamma} [GeV];#theta_{#gamma} [deg];#Delta#mu [GeV/c^{2}]",
                        fl.c_str()));

    const std::string out_root = "ecal_pi0_2dscan_corr_" + fl + "_sim.root";
    const std::string out_pdf  = "plots/output/ecal_pi0_2dscan_corr_" + fl + "_sim.pdf";
    TFile* fout = TFile::Open(out_root.c_str(), "RECREATE");

    std::vector<TCanvas*> canvases;
    canvases.reserve(nE * nTh);

    std::cout << "\n=== Cell-by-cell μ_corr ===\n";
    for (int ei = 0; ei < nE; ++ei) {
        for (int ti = 0; ti < nTh; ++ti) {
            const std::string tag =
                TString::Format("%s_e%02d_t%02d", fl.c_str(), ei, ti).Data();
            TH1D* h_sig = buildSignalCorr(cells_all[ei][ti], tag);
            const double sig_int = h_sig->Integral(
                h_sig->FindBin(kFitMin), h_sig->FindBin(kFitMax));
            CellFitC r = fitCellCorr(h_sig, tag);

            std::cout << "  E[" << kE_edges[ei] << "," << kE_edges[ei + 1] << ")"
                      << " × θ[" << kTheta_edges[ti] << "," << kTheta_edges[ti + 1] << ")"
                      << "  sig=" << std::fixed << std::setprecision(0) << sig_int
                      << "  μ_corr=" << std::setprecision(4) << r.mu
                      << " ± " << r.mu_err
                      << "  Δ=" << std::setprecision(4) << (r.mu - kPi0PDG)
                      << "  ok=" << (r.ok ? "Y" : "N") << "\n";

            if (r.ok && r.mu > 0 && sig_int > kMinSignal) {
                h_mu_corr->SetBinContent(ei + 1, ti + 1, r.mu);
                h_dmu    ->SetBinContent(ei + 1, ti + 1, r.mu - kPi0PDG);
            }

            const TString ctitle = TString::Format(
                "E_{#gamma}#in[%g,%g) GeV, #theta_{#gamma}#in[%g,%g)°  "
                "#mu_{corr}=%.4f#pm%.4f;"
                "M_{e^{+}e^{-}#gamma}^{corr} [GeV/c^{2}];Signal counts",
                kE_edges[ei], kE_edges[ei + 1],
                kTheta_edges[ti], kTheta_edges[ti + 1],
                r.mu, r.mu_err);
            h_sig->SetTitle(ctitle);
            h_sig->SetMarkerStyle(20);
            h_sig->SetMarkerSize(0.7);
            h_sig->SetMarkerColor(kBlue);
            h_sig->SetLineColor(kBlue);
            h_sig->GetXaxis()->SetRangeUser(0.0, 0.30);

            const TString cname = TString::Format(
                "c_2dscan_corr_%s_e%02d_t%02d", fl.c_str(), ei, ti);
            auto* c = new TCanvas(cname, cname, 800, 600);
            c->SetMargin(0.13, 0.05, 0.12, 0.08);
            h_sig->Draw("E");

            if (r.ok) {
                auto* fdraw = new TF1(("fdrawc_" + tag).c_str(),
                                      FitFnCorr, kFitMin, kFitMax, 11);
                for (int p = 0; p < 11; ++p) fdraw->SetParameter(p, r.pars[p]);
                fdraw->SetLineColor(kRed);
                fdraw->SetLineWidth(2);
                fdraw->Draw("SAME");
            }
            auto* lpdg = new TLine(kPi0PDG, 0, kPi0PDG,
                                   h_sig->GetMaximum() * 1.05);
            lpdg->SetLineColor(kBlack);
            lpdg->SetLineStyle(2);
            lpdg->Draw();
            c->Update();
            canvases.push_back(c);
            fout->cd();
            h_sig->Write(("h_sig_corr_" + tag).c_str());
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

    auto drawHeat = [&](TH2D* h, const std::string& base, const char* opt,
                        double zmin, double zmax) {
        TString cname = TString::Format("c_%s", h->GetName());
        auto* c = new TCanvas(cname, cname, 1100, 800);
        c->SetMargin(0.13, 0.15, 0.12, 0.10);
        c->SetGrid();
        h->SetMinimum(zmin);
        h->SetMaximum(zmax);
        h->Draw(opt);
        c->Update();
        c->SaveAs(("plots/output/" + base + "_" + fl + ".pdf").c_str());
        c->SaveAs(("plots/output/" + base + "_" + fl + ".png").c_str());
    };

    drawHeat(h_mu_corr, "ecal_pi0_2dscan_corr_sim_mu",
             "COLZ TEXT45", 0.130, 0.140);
    drawHeat(h_dmu,     "ecal_pi0_2dscan_corr_sim_dmu",
             "COLZ TEXT45", -0.005, 0.005);

    // Terminal table.
    std::cout << "\n=== ECAL π⁰ 2D scan (CORRECTED) summary (" << fl << ") ===\n";
    std::cout << "  rows = θ bins,  cols = E_γ bins\n      ";
    for (int ei = 0; ei < nE; ++ei)
        std::cout << "  E[" << kE_edges[ei] << "," << kE_edges[ei + 1] << ")";
    std::cout << "\n";
    for (int ti = 0; ti < nTh; ++ti) {
        std::cout << "θ[" << kTheta_edges[ti]
                  << "," << kTheta_edges[ti + 1] << "):  ";
        for (int ei = 0; ei < nE; ++ei) {
            const double m = h_mu_corr->GetBinContent(ei + 1, ti + 1);
            if (m > 0)
                std::cout << "  " << std::fixed << std::setprecision(4) << m;
            else
                std::cout << "    -   ";
        }
        std::cout << "\n";
    }

    fout->cd();
    h_mu_corr->Write();  h_dmu->Write();
    fout->Close();

    std::cout << "\nWrote:\n"
              << "  " << out_root << "  (heatmaps + sig hists)\n"
              << "  " << out_pdf  << "  (multipage cells)\n"
              << "  plots/output/ecal_pi0_2dscan_corr_sim_{mu,dmu}_" << fl << ".{pdf,png}\n";
}
