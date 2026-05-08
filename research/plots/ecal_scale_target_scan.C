// research/plots/ecal_scale_target_scan.C — Phase A of constraint-based
// ECAL energy-scale calibration.
//
// IDEA: every e+e-γ event near a known meson peak gives a per-event
// measurement of the ECAL energy-scale correction. Using the invariant-mass
// constraint
//
//     M² = m_ee² + 2·E_γ_true·(E_ee − p_ee·cos α)
//                 ≡ m_ee² + 2·E_γ_true·D
//
// with M = m(π⁰) for events near the π⁰ peak (or M = m(η) near the η peak),
// and D recomputed from raw lepton kinematics, we get the "ideal" photon
// energy:
//
//     E_γ_constraint = (M² − m_ee²) / (2·D)
//
// and the per-event scale:
//
//     scale_target = E_γ_constraint / E_γ_measured
//
// Each event provides one independent measurement.  Histogramming
// scale_target − 1 in (E_γ, θ_γ, φ_γ, OA) gives the SHAPE of the
// calibration error vs each kinematic axis — input for Phase B (smooth
// model fit). π⁰ events dominate the low-E_γ range, η events extend
// the coverage to higher E_γ where π⁰ has no leverage.
//
// Background under the mass peaks contributes wrong scale_target values,
// but they form a broad distribution. Per-cell MEDIAN is robust against
// these tails — no explicit CB subtraction needed at this diagnostic
// stage. (For final calibration we can switch to sideband-subtracted
// or CB-subtracted means.)
//
// Usage (from research/):
//   root -l -b -q plots/ecal_scale_target_scan.C            # exp REC
//   root -l -b -q 'plots/ecal_scale_target_scan.C("sim")'   # sim
//
// Output:
//   ecal_scale_target_<dataset>.root  with:
//     - TTree scale_events   (one row per accepted event)
//     - TH1D  h_scale_vs_logE, h_scale_vs_theta, h_scale_vs_phi, h_scale_vs_oa
//   plots/output/ecal_scale_target_proj_<dataset>.{pdf,png}

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TVector3.h>
#include <TLorentzVector.h>
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
#include <algorithm>

namespace {
    // PDG masses
    constexpr double kPi0PDG = 0.13498;
    constexpr double kEtaPDG = 0.5478;
    constexpr double kMe     = 0.000511;

    // Tight mass windows around observed peak positions (≈ ±1σ): wider
    // windows let too much bg in and median(scale−1) gets pulled toward 0.
    // For exp π⁰ peak ≈ 0.140 ± 0.013 → window [0.130, 0.150].
    // For exp η peak ≈ 0.553 ± 0.030 → window [0.530, 0.575].
    // These also work for sim (peaks at lower mass) — the windows still
    // cover ±1σ around the observed peaks in sim.
    constexpr double kPi0WinLo = 0.128;
    constexpr double kPi0WinHi = 0.150;
    constexpr double kEtaWinLo = 0.520;
    constexpr double kEtaWinHi = 0.575;

    // Quality cuts
    constexpr double kOAMin    = 0.5;     // tight pairs are bg-dominated
    constexpr double kOAMax    = 15.0;
    constexpr double kEgMin    = 0.10;    // very low-E photons unreliable
    constexpr double kThetaMin = 15.0;
    constexpr double kThetaMax = 55.0;

    // Sanity bounds for scale_target. Outside this range the event is
    // almost certainly bg or a kinematic outlier.
    constexpr double kScaleMin = 0.6;
    constexpr double kScaleMax = 1.6;

    // Compute D = E_ee − p_ee · n_γ  (returned in GeV).
    // Inputs: lepton momenta in MeV (FAT ntuple convention) + angles in deg;
    //         γ direction angles in deg.
    double computeD(double ep_p_MeV, double ep_theta_deg, double ep_phi_deg,
                    double em_p_MeV, double em_theta_deg, double em_phi_deg,
                    double g_theta_deg, double g_phi_deg)
    {
        const double d2r = M_PI / 180.0;
        const double ep_p_GeV = ep_p_MeV * 1e-3;
        const double em_p_GeV = em_p_MeV * 1e-3;
        TVector3 p_ep, p_em, n_g;
        p_ep.SetMagThetaPhi(ep_p_GeV, ep_theta_deg * d2r, ep_phi_deg * d2r);
        p_em.SetMagThetaPhi(em_p_GeV, em_theta_deg * d2r, em_phi_deg * d2r);
        n_g .SetMagThetaPhi(1.0,      g_theta_deg  * d2r, g_phi_deg  * d2r);

        TLorentzVector lv_ep(p_ep, std::sqrt(p_ep.Mag2() + kMe * kMe));
        TLorentzVector lv_em(p_em, std::sqrt(p_em.Mag2() + kMe * kMe));
        TLorentzVector lv_ee = lv_ep + lv_em;

        const double E_ee = lv_ee.E();
        const double pdotn = lv_ee.Vect().Dot(n_g);   // p_ee · n_γ  (GeV)
        return E_ee - pdotn;
    }

    // Robust median of a vector. Returns NaN if empty.
    double median(std::vector<double>& v) {
        if (v.empty()) return std::nan("");
        std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
        return v[v.size() / 2];
    }
}

void ecal_scale_target_scan(const char* dataset = "exp") {

    const std::string ds = dataset;
    const bool is_sim = (ds == "sim" || ds == "SIM");

    const std::string in_path = is_sim
        ? "../output_epem_sim.root"
        : "../output_epem_exp.root";

    TFile* fin = TFile::Open(in_path.c_str(), "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "Cannot open " << in_path << "\n"; return;
    }
    auto* t = (TTree*)fin->Get("meson_dalitz_nt");
    if (!t) { std::cerr << "meson_dalitz_nt missing\n"; return; }

    t->SetBranchStatus("*", 0);
    for (const char* b : {
            "m_ee", "m_epemg",
            "ecal_quality_pass",
            "neutr_cluster_energy", "neutr_cluster_theta", "neutr_cluster_phi",
            "oa_epem",
            "ep_p_rec", "ep_theta", "ep_phi",
            "em_p_rec", "em_theta", "em_phi",
            "sim_genweight"
         }) {
        if (t->FindBranch(b)) t->SetBranchStatus(b, 1);
    }

    float m_ee=0, m_eg=0, ecal_q=0, ne_E=0, ne_th=0, ne_ph=0, oa=0;
    float ep_p=0, ep_th=0, ep_phi=0, em_p=0, em_th=0, em_phi=0;
    float w = 1.0f;
    t->SetBranchAddress("m_ee",                 &m_ee);
    t->SetBranchAddress("m_epemg",              &m_eg);
    t->SetBranchAddress("ecal_quality_pass",    &ecal_q);
    t->SetBranchAddress("neutr_cluster_energy", &ne_E);
    t->SetBranchAddress("neutr_cluster_theta",  &ne_th);
    t->SetBranchAddress("neutr_cluster_phi",    &ne_ph);
    t->SetBranchAddress("oa_epem",              &oa);
    t->SetBranchAddress("ep_p_rec",             &ep_p);
    t->SetBranchAddress("ep_theta",             &ep_th);
    t->SetBranchAddress("ep_phi",               &ep_phi);
    t->SetBranchAddress("em_p_rec",             &em_p);
    t->SetBranchAddress("em_theta",             &em_th);
    t->SetBranchAddress("em_phi",               &em_phi);
    if (is_sim && t->FindBranch("sim_genweight"))
        t->SetBranchAddress("sim_genweight", &w);

    gSystem->mkdir("plots/output", kTRUE);

    // Output file.
    const std::string out_root = "ecal_scale_target_" + ds + ".root";
    TFile* fout = TFile::Open(out_root.c_str(), "RECREATE");

    // Per-event TTree (small subset, just what we need for further fitting).
    TTree* tev = new TTree("scale_events",
        "Per-event scale_target with kinematics for π⁰+η ECAL calibration");
    int   meson;             // 0 = π⁰, 1 = η
    float ev_logE, ev_E, ev_th, ev_ph, ev_oa;
    float ev_scale, ev_w;
    tev->Branch("meson",  &meson);
    tev->Branch("logE",   &ev_logE);
    tev->Branch("E",      &ev_E);
    tev->Branch("theta",  &ev_th);
    tev->Branch("phi",    &ev_ph);
    tev->Branch("oa",     &ev_oa);
    tev->Branch("scale",  &ev_scale);
    tev->Branch("weight", &ev_w);

    // 1D profile-style histos: per-axis MEDIAN of (scale − 1).
    // Built post-loop from the per-event TTree to use the median operator.
    // Here we just keep accumulators.
    auto makeBins = [](double lo, double hi, int n) {
        std::vector<double> e(n + 1);
        for (int i = 0; i <= n; ++i) e[i] = lo + (hi - lo) * i / n;
        return e;
    };
    std::vector<double> binsLogE  = makeBins(std::log(0.20), std::log(2.50), 14);
    std::vector<double> binsTheta = makeBins(15.0, 55.0, 8);
    std::vector<double> binsPhi   = makeBins(0.0, 60.0, 6);   // φ_local
    std::vector<double> binsOA    = makeBins(0.5, 15.0, 14);

    // Two separate sets of bags — one for π⁰ events (mflag=0), one for η.
    // We'll plot both projections on the same canvas to see if they agree.
    std::vector<std::vector<double>> bagLogE_pi0 (binsLogE .size() - 1);
    std::vector<std::vector<double>> bagTheta_pi0(binsTheta.size() - 1);
    std::vector<std::vector<double>> bagPhi_pi0  (binsPhi  .size() - 1);
    std::vector<std::vector<double>> bagOA_pi0   (binsOA   .size() - 1);
    std::vector<std::vector<double>> bagLogE_eta (binsLogE .size() - 1);
    std::vector<std::vector<double>> bagTheta_eta(binsTheta.size() - 1);
    std::vector<std::vector<double>> bagPhi_eta  (binsPhi  .size() - 1);
    std::vector<std::vector<double>> bagOA_eta   (binsOA   .size() - 1);

    auto pickBin = [](double v, const std::vector<double>& edges) {
        if (v < edges.front() || v >= edges.back()) return -1;
        for (size_t i = 0; i + 1 < edges.size(); ++i)
            if (v >= edges[i] && v < edges[i + 1]) return (int)i;
        return -1;
    };

    const Long64_t N = t->GetEntries();
    Long64_t n_pi0 = 0, n_eta = 0, n_pass = 0;
    std::cout << "Phase A: scanning scale_target on " << N
              << " entries (" << ds << ")\n";

    for (Long64_t ev = 0; ev < N; ++ev) {
        t->GetEntry(ev);
        if (ecal_q != 1.0f) continue;
        if (oa < kOAMin || oa >= kOAMax) continue;
        if (ne_th < kThetaMin || ne_th >= kThetaMax) continue;
        const double E_GeV = ne_E / 1000.0;
        if (E_GeV < kEgMin) continue;

        // Determine meson hypothesis.
        double M = -1.0;
        int    mflag = -1;
        if (m_eg >= kPi0WinLo && m_eg < kPi0WinHi) { M = kPi0PDG; mflag = 0; ++n_pi0; }
        else if (m_eg >= kEtaWinLo && m_eg < kEtaWinHi) { M = kEtaPDG; mflag = 1; ++n_eta; }
        else continue;
        ++n_pass;

        // Compute D and scale_target.
        const double D = computeD(ep_p, ep_th, ep_phi,
                                  em_p, em_th, em_phi,
                                  ne_th, ne_ph);
        if (D <= 1e-6) continue;
        const double E_constr = (M * M - m_ee * m_ee) / (2.0 * D);
        if (E_constr <= 0) continue;
        const double scale = E_constr / E_GeV;
        if (scale < kScaleMin || scale > kScaleMax) continue;

        // φ_local = φ mod 60°
        double phi_loc = std::fmod(ne_ph + 360.0, 60.0);

        meson    = mflag;
        ev_E     = E_GeV;
        ev_logE  = std::log(E_GeV);
        ev_th    = ne_th;
        ev_ph    = ne_ph;     // full φ kept; projection histos use φ_local
        ev_oa    = oa;
        ev_scale = scale;
        ev_w     = is_sim ? w : 1.0f;
        tev->Fill();

        // Accumulate for median histos — separate by meson.
        const double sm1 = scale - 1.0;
        auto& bL = (mflag == 0) ? bagLogE_pi0  : bagLogE_eta;
        auto& bT = (mflag == 0) ? bagTheta_pi0 : bagTheta_eta;
        auto& bP = (mflag == 0) ? bagPhi_pi0   : bagPhi_eta;
        auto& bO = (mflag == 0) ? bagOA_pi0    : bagOA_eta;
        if (int b = pickBin(ev_logE, binsLogE);  b >= 0) bL[b].push_back(sm1);
        if (int b = pickBin(ev_th,   binsTheta); b >= 0) bT[b].push_back(sm1);
        if (int b = pickBin(phi_loc, binsPhi);   b >= 0) bP[b].push_back(sm1);
        if (int b = pickBin(ev_oa,   binsOA);    b >= 0) bO[b].push_back(sm1);
    }
    std::cout << "  selected: " << n_pass
              << "  (π⁰: " << n_pi0 << ",  η: " << n_eta << ")\n";

    // Build median-vs-axis histograms.
    auto buildMedianHist = [](std::vector<std::vector<double>>& bag,
                              const std::vector<double>& edges,
                              const char* name, const char* title) {
        TH1D* h = new TH1D(name, title,
                           (int)edges.size() - 1, edges.data());
        h->SetDirectory(nullptr);
        for (int b = 0; b < (int)bag.size(); ++b) {
            if (bag[b].size() < 20) continue;   // need a few events for stable median
            const double med = median(bag[b]);
            h->SetBinContent(b + 1, med);
            // 1.4826·MAD/√N as error proxy (robust)
            std::vector<double> dev;
            dev.reserve(bag[b].size());
            for (double v : bag[b]) dev.push_back(std::fabs(v - med));
            const double mad = median(dev);
            h->SetBinError(b + 1, 1.4826 * mad / std::sqrt((double)bag[b].size()));
        }
        return h;
    };

    TH1D* h_logE_pi0  = buildMedianHist(bagLogE_pi0,  binsLogE,
        "h_scale_vs_logE_pi0",  ";log(E_{#gamma}/GeV);median(scale - 1)");
    TH1D* h_theta_pi0 = buildMedianHist(bagTheta_pi0, binsTheta,
        "h_scale_vs_theta_pi0", ";#theta_{#gamma} [deg];median(scale - 1)");
    TH1D* h_phi_pi0   = buildMedianHist(bagPhi_pi0,   binsPhi,
        "h_scale_vs_philoc_pi0",";#phi_{#gamma} mod 60° [deg];median(scale - 1)");
    TH1D* h_oa_pi0    = buildMedianHist(bagOA_pi0,    binsOA,
        "h_scale_vs_oa_pi0",    ";OA(e^{+}e^{-}) [deg];median(scale - 1)");

    TH1D* h_logE_eta  = buildMedianHist(bagLogE_eta,  binsLogE,
        "h_scale_vs_logE_eta",  ";log(E_{#gamma}/GeV);median(scale - 1)");
    TH1D* h_theta_eta = buildMedianHist(bagTheta_eta, binsTheta,
        "h_scale_vs_theta_eta", ";#theta_{#gamma} [deg];median(scale - 1)");
    TH1D* h_phi_eta   = buildMedianHist(bagPhi_eta,   binsPhi,
        "h_scale_vs_philoc_eta",";#phi_{#gamma} mod 60° [deg];median(scale - 1)");
    TH1D* h_oa_eta    = buildMedianHist(bagOA_eta,    binsOA,
        "h_scale_vs_oa_eta",    ";OA(e^{+}e^{-}) [deg];median(scale - 1)");

    // Diagnostic 4-pad canvas — π⁰ (black filled circles) + η (red squares).
    gStyle->SetOptStat(0);
    auto* c = new TCanvas(("c_scale_" + ds).c_str(),
                          ("scale_target diag " + ds).c_str(), 1500, 1000);
    c->Divide(2, 2);
    auto styled = [](TH1D* h, int color, int marker) {
        h->SetMarkerStyle(marker); h->SetMarkerSize(0.9);
        h->SetMarkerColor(color); h->SetLineColor(color);
        h->SetLineWidth(2);
    };
    styled(h_logE_pi0,  kBlack, 20);
    styled(h_theta_pi0, kBlack, 20);
    styled(h_phi_pi0,   kBlack, 20);
    styled(h_oa_pi0,    kBlack, 20);
    styled(h_logE_eta,  kRed,   24);
    styled(h_theta_eta, kRed,   24);
    styled(h_phi_eta,   kRed,   24);
    styled(h_oa_eta,    kRed,   24);

    auto drawPad = [&](int pad, TH1D* h_p, TH1D* h_e, double ylo, double yhi) {
        c->cd(pad);
        gPad->SetGrid();
        gPad->SetMargin(0.13, 0.05, 0.13, 0.10);
        h_p->SetMinimum(ylo);
        h_p->SetMaximum(yhi);
        h_p->Draw("P");
        h_e->Draw("P SAME");
        const double xlo = h_p->GetXaxis()->GetXmin();
        const double xhi = h_p->GetXaxis()->GetXmax();
        auto* l0 = new TLine(xlo, 0.0, xhi, 0.0);
        l0->SetLineStyle(2); l0->SetLineColor(kBlue); l0->SetLineWidth(2);
        l0->Draw();
        if (pad == 1) {
            auto* leg = new TLegend(0.55, 0.78, 0.94, 0.92);
            leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.040);
            leg->AddEntry(h_p, "#pi^{0} events",  "lp");
            leg->AddEntry(h_e, "#eta events",     "lp");
            leg->AddEntry(l0,  "scale = 1",       "l");
            leg->Draw();
        }
    };
    const double ylo = -0.15;
    const double yhi =  0.20;
    drawPad(1, h_logE_pi0,  h_logE_eta,  ylo, yhi);
    drawPad(2, h_theta_pi0, h_theta_eta, ylo, yhi);
    drawPad(3, h_phi_pi0,   h_phi_eta,   ylo, yhi);
    drawPad(4, h_oa_pi0,    h_oa_eta,    ylo, yhi);

    c->Update();
    const std::string out_base = "plots/output/ecal_scale_target_proj_" + ds;
    c->SaveAs((out_base + ".pdf").c_str());
    c->SaveAs((out_base + ".png").c_str());

    // Save TTree + histograms.
    fout->cd();
    tev->Write();
    h_logE_pi0->Write(); h_theta_pi0->Write();
    h_phi_pi0 ->Write(); h_oa_pi0   ->Write();
    h_logE_eta->Write(); h_theta_eta->Write();
    h_phi_eta ->Write(); h_oa_eta   ->Write();
    fout->Close();
    fin->Close();

    std::cout << "\nWrote:\n"
              << "  " << out_root << "  (TTree scale_events + 4 histos)\n"
              << "  " << out_base << ".{pdf,png}\n";
}
