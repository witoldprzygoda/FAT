// research/plots/ecal_scale_target_3dmap.C — differential ECAL energy-scale
// correction map s(E_γ, θ_γ, φ_γ) built from per-event scale_target values.
//
// Method: every e+e-γ event near a known meson peak (π⁰ or η) gives one
// independent measurement of the scale_target via the constraint
//
//     M² = m_ee² + 2·E_γ_true·D    where D = E_ee − p_ee·cos α
//     scale_target = (M² − m_ee²) / (2·D·E_γ_meas)
//
// The cell map is then the median of scale_target over events whose
// (E_γ, θ_γ, φ_γ) falls inside the cell. Median is robust against
// background contamination in the mass window. φ is binned over the FULL
// [0°, 360°) range — every sector and every position within sector gets
// its own value. No folding, no cross-sector averaging, no accumulation.
//
// Usage (from research/):
//   root -l -b -q 'plots/ecal_scale_target_3dmap.C("exp")'  // exp data
//   root -l -b -q 'plots/ecal_scale_target_3dmap.C("sim")'  // sim
//
// Output:
//   ecal_scale_target_3dmap_<dataset>.root  with TH3D h_s_scaletgt_3d
//      (calibration map; loadable by EcalLookup::load3D)
//   plots/output/ecal_scale_target_3dmap_<dataset>_phi_pattern.png
//      (averaged-over-(E,θ) view; visible dead-sector gap).

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TCanvas.h>
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
    constexpr double kPi0PDG = 0.13498;
    constexpr double kEtaPDG = 0.5478;
    constexpr double kMe     = 0.000511;

    // Tight mass windows (≈ ±1σ around observed peaks).
    constexpr double kPi0WinLo = 0.128;
    constexpr double kPi0WinHi = 0.150;
    constexpr double kEtaWinLo = 0.520;
    constexpr double kEtaWinHi = 0.575;

    constexpr double kOAMin = 0.5;
    constexpr double kOAMax = 15.0;
    constexpr double kEgMin = 0.10;

    // Sanity bounds for scale_target.
    constexpr double kScaleMin = 0.6;
    constexpr double kScaleMax = 1.6;

    // 3D grid for the map. Same E_γ × θ_γ binning as the previous 4×3,
    // and 36 phi-bins of 10° over the FULL [0, 360°) range.
    const std::vector<double> kE_edges     = {0.20, 0.45, 0.75, 1.20, 2.50};
    const std::vector<double> kTheta_edges = {15.0, 28.0, 42.0, 55.0};
    std::vector<double> makePhiEdges() {
        std::vector<double> v; v.reserve(37);
        for (int i = 0; i <= 36; ++i) v.push_back(i * 10.0);
        return v;
    }
    const std::vector<double> kPhi_edges = makePhiEdges();

    // Minimum events in a cell to compute a stable median.
    constexpr int kMinCellEvents = 20;

    int findBin(const std::vector<double>& edges, double v) {
        if (v < edges.front() || v >= edges.back()) return -1;
        for (size_t i = 0; i + 1 < edges.size(); ++i)
            if (v >= edges[i] && v < edges[i + 1]) return (int)i;
        return -1;
    }

    double median(std::vector<double>& v) {
        if (v.empty()) return std::nan("");
        std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
        return v[v.size() / 2];
    }
}

void ecal_scale_target_3dmap(const char* dataset = "exp") {

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
            "m_ee", "m_epemg", "ecal_quality_pass",
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

    const int nE  = (int)kE_edges.size()     - 1;
    const int nTh = (int)kTheta_edges.size() - 1;
    const int nPh = (int)kPhi_edges.size()   - 1;
    const int n_cells = nE * nTh * nPh;

    // Per-cell containers for scale_target values.
    std::vector<std::vector<std::vector<std::vector<double>>>> bag(
        nE, std::vector<std::vector<std::vector<double>>>(
            nTh, std::vector<std::vector<double>>(nPh)));

    const Long64_t N = t->GetEntries();
    std::cout << "Building 3D scale_target map on " << N
              << " entries (" << ds << ")\n";

    const double d2r = M_PI / 180.0;
    Long64_t n_pi0 = 0, n_eta = 0;
    for (Long64_t ev = 0; ev < N; ++ev) {
        t->GetEntry(ev);
        if (ecal_q != 1.0f) continue;
        if (oa < kOAMin || oa >= kOAMax) continue;
        const double E_GeV = ne_E / 1000.0;
        if (E_GeV < kEgMin) continue;

        double M = -1.0;
        if      (m_eg >= kPi0WinLo && m_eg < kPi0WinHi) { M = kPi0PDG; ++n_pi0; }
        else if (m_eg >= kEtaWinLo && m_eg < kEtaWinHi) { M = kEtaPDG; ++n_eta; }
        else continue;

        // D from raw lepton kinematics. Lepton momenta are in MeV.
        TVector3 p_ep, p_em, n_g;
        p_ep.SetMagThetaPhi(ep_p * 1e-3, ep_th * d2r, ep_phi * d2r);
        p_em.SetMagThetaPhi(em_p * 1e-3, em_th * d2r, em_phi * d2r);
        n_g .SetMagThetaPhi(1.0,         ne_th * d2r, ne_ph  * d2r);
        TLorentzVector lv_ep(p_ep, std::sqrt(p_ep.Mag2() + kMe * kMe));
        TLorentzVector lv_em(p_em, std::sqrt(p_em.Mag2() + kMe * kMe));
        TLorentzVector lv_ee = lv_ep + lv_em;
        const double D = lv_ee.E() - lv_ee.Vect().Dot(n_g);
        if (D <= 1e-6) continue;

        const double E_constr = (M * M - m_ee * m_ee) / (2.0 * D);
        if (E_constr <= 0) continue;
        const double scale = E_constr / E_GeV;
        if (scale < kScaleMin || scale > kScaleMax) continue;

        const int ei = findBin(kE_edges,     E_GeV);
        const int ti = findBin(kTheta_edges, ne_th);
        const int pi = findBin(kPhi_edges,   ne_ph);   // FULL phi
        if (ei < 0 || ti < 0 || pi < 0) continue;

        bag[ei][ti][pi].push_back(scale);
    }
    std::cout << "  used π⁰ events: " << n_pi0
              << "  η events: " << n_eta << "\n";

    // -- Build TH3D with median per cell -----------------------------------
    auto* h_s = new TH3D(("h_s_scaletgt_3d"),
        TString::Format("scale_target median per (E, #theta, #phi) — %s;"
                        "E_{#gamma} [GeV];#theta_{#gamma} [deg];#phi_{#gamma} [deg]",
                        ds.c_str()),
        nE,  kE_edges.data(),
        nTh, kTheta_edges.data(),
        nPh, kPhi_edges.data());
    h_s->SetDirectory(nullptr);
    auto* h_n = new TH3D(("h_n_scaletgt_3d"),
        TString::Format("Cell event count — %s;"
                        "E_{#gamma} [GeV];#theta_{#gamma} [deg];#phi_{#gamma} [deg]",
                        ds.c_str()),
        nE,  kE_edges.data(),
        nTh, kTheta_edges.data(),
        nPh, kPhi_edges.data());
    h_n->SetDirectory(nullptr);

    int n_filled = 0;
    for (int ei = 0; ei < nE; ++ei)
        for (int ti = 0; ti < nTh; ++ti)
            for (int pi = 0; pi < nPh; ++pi) {
                auto& v = bag[ei][ti][pi];
                h_n->SetBinContent(ei+1, ti+1, pi+1, (double)v.size());
                if ((int)v.size() < kMinCellEvents) continue;
                const double med = median(v);
                h_s->SetBinContent(ei+1, ti+1, pi+1, med);
                ++n_filled;
            }
    std::cout << "  cells filled: " << n_filled << " / " << n_cells
              << " (" << (100.0 * n_filled / n_cells) << "%)\n";

    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);
    gSystem->mkdir("plots/output", kTRUE);

    // -- Quick projection: median(s) vs phi, averaged over (E, θ) cells with
    //    enough stats. Sectors are 60° apart; dead sector should show as a gap.
    {
        TH1D h_phi("h_scale_vs_phi",
                   ";#phi_{#gamma} [deg];median scale_target",
                   nPh, kPhi_edges.data());
        for (int pi = 0; pi < nPh; ++pi) {
            std::vector<double> vs;
            for (int ei = 0; ei < nE; ++ei)
                for (int ti = 0; ti < nTh; ++ti) {
                    const double v = h_s->GetBinContent(ei+1, ti+1, pi+1);
                    if (v > 0) vs.push_back(v);
                }
            if (vs.size() >= 3)
                h_phi.SetBinContent(pi + 1, median(vs));
        }
        TCanvas c("c_phi", "phi pattern", 1300, 600);
        c.SetGrid();
        c.SetMargin(0.10, 0.05, 0.13, 0.08);
        h_phi.SetMinimum(0.85); h_phi.SetMaximum(1.15);
        h_phi.SetMarkerStyle(20); h_phi.SetMarkerSize(1.0);
        h_phi.SetMarkerColor(kBlack); h_phi.SetLineColor(kBlack);
        h_phi.Draw("P");
        // Sector boundaries
        for (int s = 1; s < 6; ++s) {
            auto* l = new TLine(s * 60.0, 0.85, s * 60.0, 1.15);
            l->SetLineStyle(3); l->SetLineColor(kGray + 2);
            l->Draw();
        }
        auto* l1 = new TLine(0, 1.0, 360, 1.0);
        l1->SetLineStyle(2); l1->SetLineColor(kBlue); l1->SetLineWidth(2);
        l1->Draw();
        const std::string out = "plots/output/ecal_scale_target_3dmap_"
                              + ds + "_phi_pattern";
        c.SaveAs((out + ".pdf").c_str());
        c.SaveAs((out + ".png").c_str());
        std::cout << "  phi-pattern: " << out << ".{pdf,png}\n";
    }

    // Save
    const std::string out_root = "ecal_scale_target_3dmap_" + ds + ".root";
    TFile* fout = TFile::Open(out_root.c_str(), "RECREATE");
    h_s->Write();
    h_n->Write();
    fout->Close();
    fin->Close();

    std::cout << "Wrote: " << out_root << "  (TH3D h_s_scaletgt_3d, h_n_scaletgt_3d)\n";
}
