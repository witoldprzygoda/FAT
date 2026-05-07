// research/plots/ecal_eta_validate_corr.C — INDEPENDENT validation of the
// ECAL energy-scale correction on η Dalitz.
//
// The 2D lookup map s(E_γ, θ_γ) was trained on π⁰. If the underlying
// physics (ECAL energy mis-calibration) is real and not a coincidence
// of π⁰ kinematics, applying the same correction to η decays should
// move the η peak toward its PDG value 0.5478 GeV/c² as well.
//
// What this macro does:
//   1) loads the calibration map ecal_pi0_2dscan_<fl>.root
//   2) loops meson_dalitz_nt across the 3 channels (epem, epep, emem)
//      and fills BOTH:
//        - h_eta_uncorr_*  (m_epemg as written by main.cc)
//        - h_eta_corr_*    (m_corr from the lookup formula)
//   3) builds CB-subtracted signal for each
//   4) fits both with the same η model (CB + P2 + exp) and reports
//      μ_uncorr vs μ_corr — winner: closer to 0.5478.
//
// Usage (from research/):
//   root -l -b -q plots/ecal_eta_validate_corr.C            # REC
//   root -l -b -q 'plots/ecal_eta_validate_corr.C("cor")'   # COR

#include "ecal_apply_lookup.h"

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
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
    constexpr double kEtaPDG = 0.5478;

    // η mass-window histogram binning.
    constexpr int    kHistNBins = 200;
    constexpr double kHistMin   = 0.30;
    constexpr double kHistMax   = 0.80;

    // η fit window — wider on the left for the CB tail.
    constexpr double kFitMin = 0.40;
    constexpr double kFitMax = 0.70;

    // Pre-fit Gaussian apex window — looks for peak near PDG.
    constexpr double kPreMin = 0.520;
    constexpr double kPreMax = 0.580;
}

// CB(left) + P2 + exp — simpler model since the η bg is smoother.
double CBLeftEta(double *x, double *par) {
    const double m = x[0], N = par[0], mu = par[1], s = par[2],
                 a = par[3], n = par[4];
    if (s <= 0 || a <= 0 || n <= 0) return 0.0;
    const double t = (m - mu) / s, aA = std::fabs(a);
    if (t > -aA) return N * std::exp(-0.5 * t * t);
    const double A = std::pow(n / aA, n) * std::exp(-0.5 * aA * aA);
    return N * A * std::pow(n / aA - aA - t, -n);
}
double FitFnEta(double *x, double *par) {
    const double m  = x[0];
    const double cb = CBLeftEta(x, par);
    const double p2 = par[5] + par[6]*m + par[7]*m*m;
    const double ex = par[8] * std::exp(-par[9]*m);
    return cb + p2 + ex;
}

struct EtaFit {
    double mu = 0,    mu_err = 0;
    double sigma = 0, sigma_err = 0;
    double pars[10] = {0};
    bool   ok = false;
};

EtaFit fitEta(TH1D* h, const std::string& base) {
    EtaFit r;
    if (!h || h->GetEntries() == 0) return r;

    // Pre-find apex on prefit window.
    const int b_lo = h->FindBin(kPreMin), b_hi = h->FindBin(kPreMax);
    int    pb = b_lo;
    double pv = h->GetBinContent(b_lo);
    for (int b = b_lo; b <= b_hi; ++b)
        if (h->GetBinContent(b) > pv) { pv = h->GetBinContent(b); pb = b; }

    auto* fit = new TF1(("fit_eta_" + base).c_str(),
                        FitFnEta, kFitMin, kFitMax, 10);
    fit->SetParameter(0,  std::max(1.0, pv));
    fit->SetParameter(1,  h->GetBinCenter(pb));
    fit->SetParameter(2,  0.020);
    fit->SetParameter(3,  1.5);
    fit->SetParameter(4,  5.0);
    fit->SetParameter(5,  0.0);
    fit->SetParameter(6,  0.0);
    fit->SetParameter(7,  0.0);
    fit->SetParameter(8,  pv * 0.1);
    fit->SetParameter(9,  3.0);
    fit->SetParLimits(0, 0.0,    1e12);
    fit->SetParLimits(1, 0.500,  0.590);
    fit->SetParLimits(2, 0.005,  0.060);
    fit->SetParLimits(3, 0.30,   5.0);
    fit->SetParLimits(4, 1.0,    50.0);
    fit->SetParLimits(8, 0.0,    1e12);
    fit->SetParLimits(9, 0.0,    50.0);

    int status = h->Fit(fit, "RQ");
    if (status != 0) {
        fit->SetParameter(2, 0.025);
        h->Fit(fit, "RQ");
    }
    r.mu        = fit->GetParameter(1);
    r.mu_err    = fit->GetParError(1);
    r.sigma     = fit->GetParameter(2);
    r.sigma_err = fit->GetParError(2);
    for (int i = 0; i < 10; ++i) r.pars[i] = fit->GetParameter(i);
    r.ok = (r.mu > 0.510 && r.mu < 0.585 && r.mu_err > 0 && r.mu_err < 0.005);
    delete fit;
    return r;
}

// Loop one channel and fill 2 histograms (uncorrected + corrected).
void fillEtaChannel(const std::string& fpath,
                    const EcalLookup& look,
                    TH1D* h_uncorr, TH1D* h_corr,
                    const std::string& tag)
{
    TFile* f = TFile::Open(fpath.c_str(), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Cannot open " << fpath << "\n"; return;
    }
    auto* t = (TTree*)f->Get("meson_dalitz_nt");
    if (!t) { std::cerr << "meson_dalitz_nt missing in " << fpath << "\n"; return; }

    t->SetBranchStatus("*", 0);
    for (const char* b : {
            "m_ee",
            "m_epemg",
            "ecal_quality_pass",
            "neutr_cluster_energy",
            "neutr_cluster_theta",
            "neutr_cluster_phi",
            "gamma_D",
            "oa_epem"
         }) t->SetBranchStatus(b, 1);

    float m_ee=0, m_eg=0, ecal_q=0, ne_E=0, ne_th=0, ne_ph=0, gD=0, oa=0;
    t->SetBranchAddress("m_ee",                 &m_ee);
    t->SetBranchAddress("m_epemg",              &m_eg);
    t->SetBranchAddress("ecal_quality_pass",    &ecal_q);
    t->SetBranchAddress("neutr_cluster_energy", &ne_E);
    t->SetBranchAddress("neutr_cluster_theta",  &ne_th);
    t->SetBranchAddress("neutr_cluster_phi",    &ne_ph);
    t->SetBranchAddress("gamma_D",              &gD);
    t->SetBranchAddress("oa_epem",              &oa);

    const Long64_t N = t->GetEntries();
    std::cout << "  " << fpath << " (" << tag << "): " << N << " entries\n";

    // η Dalitz analysis enhances signal/bg by requiring OA(e+e-) > 4°.
    // This same window is used by fit_eta.C / yields_eta.C and is what
    // separates η from the π⁰ Dalitz (which is at smaller OA) and from
    // the broad combinatorial bg at very small OA.
    constexpr double kEtaOAMin = 4.0;
    constexpr double kEtaOAMax = 15.0;

    for (Long64_t ev = 0; ev < N; ++ev) {
        t->GetEntry(ev);
        if (ecal_q != 1.0f) continue;
        if (oa < kEtaOAMin || oa >= kEtaOAMax) continue;
        const double E_GeV = ne_E / 1000.0;
        const double m_corr = look.m_corrected(m_ee, E_GeV, gD, ne_th, ne_ph);
        h_uncorr->Fill(m_eg);
        h_corr  ->Fill(m_corr);
    }
    f->Close();
}

// Build CB-subtracted signal (all - 2√(N₊₊N₋₋)) per uncorr/corr stream.
TH1D* makeSignal(TH1D* h_all, TH1D* h_pp, TH1D* h_mm, const std::string& tag) {
    TH1D* h_cb = (TH1D*)h_all->Clone(("h_cb_" + tag).c_str());
    h_cb->Reset();
    h_cb->SetDirectory(nullptr);
    for (int b = 1; b <= h_cb->GetNbinsX(); ++b) {
        const double n_pp = h_pp->GetBinContent(b);
        const double n_mm = h_mm->GetBinContent(b);
        const double e_pp = h_pp->GetBinError(b);
        const double e_mm = h_mm->GetBinError(b);
        if (n_pp > 0 && n_mm > 0) {
            const double cb = 2.0 * std::sqrt(n_pp * n_mm);
            const double t1 = e_pp * std::sqrt(n_mm / n_pp);
            const double t2 = e_mm * std::sqrt(n_pp / n_mm);
            h_cb->SetBinContent(b, cb);
            h_cb->SetBinError(b, std::sqrt(t1 * t1 + t2 * t2));
        }
    }
    TH1D* h_sig = (TH1D*)h_all->Clone(("h_sig_" + tag).c_str());
    h_sig->SetDirectory(nullptr);
    h_sig->Add(h_cb, -1.0);
    delete h_cb;
    return h_sig;
}

void ecal_eta_validate_corr(const char* flavour = "rec",
                            const char* map_dim = "2d") {

    std::string fl = flavour;
    for (auto& c : fl) c = std::tolower(c);
    std::string mass_var;
    if      (fl == "rec") mass_var = "m_epemg";
    else if (fl == "cor") mass_var = "m_epemg_cor";
    else { std::cerr << "Unknown flavour '" << flavour << "'\n"; return; }
    (void)mass_var;

    std::string md = map_dim ? map_dim : "2d";
    for (auto& c : md) c = std::tolower(c);

    EcalLookup look;
    std::string map_path, map_hist;
    bool loaded = false;
    if (md == "3d") {
        map_path = "ecal_pi0_3dscan_" + fl + ".root";
        map_hist = "h_s_" + fl + "_3d";
        loaded = look.load3D(map_path, map_hist);
    } else {
        map_path = "ecal_pi0_2dscan_" + fl + ".root";
        map_hist = "h_s_" + fl;
        loaded = look.load(map_path, map_hist);
    }
    if (!loaded) {
        std::cerr << "Run ecal_pi0_" << md << "scan.C first to produce "
                  << map_path << "\n";
        return;
    }
    const std::string map_tag = (md == "3d") ? "_3dmap" : "";

    gStyle->SetOptStat(0);
    gSystem->mkdir("plots/output", kTRUE);

    auto mk = [&](const std::string& nm) {
        TH1D* h = new TH1D(nm.c_str(),
                           ";M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Counts",
                           kHistNBins, kHistMin, kHistMax);
        h->Sumw2();
        h->SetDirectory(nullptr);
        return h;
    };
    TH1D* h_all_u = mk("h_eta_all_uncorr");
    TH1D* h_all_c = mk("h_eta_all_corr");
    TH1D* h_pp_u  = mk("h_eta_pp_uncorr");
    TH1D* h_pp_c  = mk("h_eta_pp_corr");
    TH1D* h_mm_u  = mk("h_eta_mm_uncorr");
    TH1D* h_mm_c  = mk("h_eta_mm_corr");

    std::cout << "ECAL η validation, flavour=" << fl << "\n";
    fillEtaChannel("../output_epem_exp.root", look, h_all_u, h_all_c, "all");
    fillEtaChannel("../output_epep_exp.root", look, h_pp_u,  h_pp_c,  "pp ");
    fillEtaChannel("../output_emem_exp.root", look, h_mm_u,  h_mm_c,  "mm ");

    TH1D* h_sig_u = makeSignal(h_all_u, h_pp_u, h_mm_u, "u_" + fl);
    TH1D* h_sig_c = makeSignal(h_all_c, h_pp_c, h_mm_c, "c_" + fl);

    EtaFit f_u = fitEta(h_sig_u, "uncorr_" + fl);
    EtaFit f_c = fitEta(h_sig_c, "corr_"   + fl);

    std::cout << "\n=== η peak validation (" << fl << ") ===\n";
    std::cout << std::fixed << std::setprecision(5);
    std::cout << "  PDG η mass:     " << kEtaPDG << "\n";
    std::cout << "  μ uncorrected:  " << f_u.mu << " ± " << f_u.mu_err
              << "  Δ=" << (f_u.mu - kEtaPDG)
              << "  σ=" << f_u.sigma
              << "  ok=" << (f_u.ok ? "Y" : "N") << "\n";
    std::cout << "  μ corrected:    " << f_c.mu << " ± " << f_c.mu_err
              << "  Δ=" << (f_c.mu - kEtaPDG)
              << "  σ=" << f_c.sigma
              << "  ok=" << (f_c.ok ? "Y" : "N") << "\n";
    const double improv = std::fabs(f_u.mu - kEtaPDG) - std::fabs(f_c.mu - kEtaPDG);
    std::cout << "  Improvement |Δμ| reduced by: " << improv * 1000.0 << " MeV\n";

    // Comparison plot: uncorr vs corr signal histograms with fits.
    auto draw = [&](TH1D* h, EtaFit& f, int color, const std::string& label) {
        h->SetMarkerStyle(20);
        h->SetMarkerSize(0.7);
        h->SetMarkerColor(color);
        h->SetLineColor(color);
        return h;
    };

    draw(h_sig_u, f_u, kBlue,  "uncorrected");
    draw(h_sig_c, f_c, kRed,   "corrected");

    const TString ctitle = TString::Format(
        "#eta validation (%s);M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Signal counts",
        fl.c_str());
    h_sig_u->SetTitle(ctitle);
    h_sig_u->GetXaxis()->SetRangeUser(0.40, 0.70);

    auto* c = new TCanvas(("c_eta_val_" + fl).c_str(),
                          "eta validation", 1100, 700);
    c->SetMargin(0.13, 0.05, 0.12, 0.08);
    c->SetGrid();
    h_sig_u->Draw("E");
    h_sig_c->Draw("E SAME");

    auto* fdraw_u = new TF1(("fd_u_" + fl).c_str(),
                            FitFnEta, kFitMin, kFitMax, 10);
    for (int p = 0; p < 10; ++p) fdraw_u->SetParameter(p, f_u.pars[p]);
    fdraw_u->SetLineColor(kBlue);   fdraw_u->SetLineStyle(2);
    fdraw_u->SetLineWidth(2);       fdraw_u->Draw("SAME");

    auto* fdraw_c = new TF1(("fd_c_" + fl).c_str(),
                            FitFnEta, kFitMin, kFitMax, 10);
    for (int p = 0; p < 10; ++p) fdraw_c->SetParameter(p, f_c.pars[p]);
    fdraw_c->SetLineColor(kRed);    fdraw_c->SetLineStyle(2);
    fdraw_c->SetLineWidth(2);       fdraw_c->Draw("SAME");

    auto* lpdg = new TLine(kEtaPDG, 0, kEtaPDG, h_sig_u->GetMaximum() * 1.05);
    lpdg->SetLineColor(kBlack);
    lpdg->SetLineStyle(2);
    lpdg->Draw();

    auto* leg = new TLegend(0.55, 0.65, 0.94, 0.92);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.030);
    leg->AddEntry(h_sig_u, TString::Format("uncorr: #mu=%.4f #pm %.4f",
                                           f_u.mu, f_u.mu_err), "lpe");
    leg->AddEntry(h_sig_c, TString::Format("corr:   #mu=%.4f #pm %.4f",
                                           f_c.mu, f_c.mu_err), "lpe");
    leg->AddEntry(lpdg,    TString::Format("PDG #eta = %.4f", kEtaPDG), "l");
    leg->Draw();

    c->Update();
    const std::string base = "plots/output/ecal_eta_validate_corr_" + fl + map_tag;
    c->SaveAs((base + ".pdf").c_str());
    c->SaveAs((base + ".png").c_str());

    // Save per-stream hists for later inspection.
    TFile* fout = TFile::Open(("ecal_eta_validate_corr_" + fl + map_tag + ".root").c_str(),
                              "RECREATE");
    h_sig_u->Write();
    h_sig_c->Write();
    h_all_u->Write(); h_pp_u->Write(); h_mm_u->Write();
    h_all_c->Write(); h_pp_c->Write(); h_mm_c->Write();
    fout->Close();

    std::cout << "\nWrote: " << base << ".{pdf,png}\n";
}
