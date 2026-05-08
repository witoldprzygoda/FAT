// research/plots/apply_ecal_correction.C — produce research_sim_ecalcor.root.
//
// Goal (SIM): bring the post-hoc ECAL energy-scale correction into the
// standard sim pipeline. Output mirrors `meson_research`'s SIM layout
// (single file, full + OA-sliced histograms, three mass flavours
// REC/COR/TRU, all weighted by sim_genweight), so the existing downstream
// tooling (fit_pi0.C, yields_pi0.C, fit_eta.C, yields_eta.C, …) can read
// the corrected mass with a single path swap.
//
// Closed-form correction used here:
//     m²_corr = (1 − s) · m_ee² + s · m_epemg²
// derived from m² = m_ee² + 2·E_γ·D and m²_corr = m_ee² + 2·s·E_γ·D.
// `m_ee` is the REC dilepton mass — same value used for all three flavours
// (REC/COR/TRU) since the SIM ntuple does not store m_ee_cor / m_ee_sim
// separately at the meson_dalitz_nt level. TRU is shipped uncorrected
// (m_epemg_tru = m_epemg_sim) — calorimeter response is irrelevant for
// MC truth.
//
// Per event passing ecal_quality_pass==1, fills:
//   - m_epemg_full       (REC corrected, integrated over OA)
//   - m_epemg_cor_full   (COR corrected, integrated over OA)
//   - m_epemg_tru_full   (TRU uncorrected, integrated over OA)
//   - m_epemg_oa_<lo>_<hi>      (REC, fine 0.2° → 75 slices, π⁰)
//   - m_epemg_oa_<lo>_<hi>      (REC, coarse 0.5° → 30 slices, η)
//   - same for m_epemg_cor_*  and  m_epemg_tru_*
// All fills weighted by sim_genweight.
//
// Usage (from research/):
//   root -l -b -q plots/apply_ecal_correction.C
//
// Pre-requisite: ecal_pi0_2dscan_rec_sim.root must exist (run
// ecal_pi0_2dscan.C first to build the sim calibration map).

#include "ecal_apply_lookup.h"

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TString.h>
#include <TSystem.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <tuple>
#include <vector>

namespace {
    constexpr double kSliceMin     = 0.0;
    constexpr double kSliceMax     = 15.0;
    constexpr double kSliceStep    = 0.2;    // π⁰
    constexpr double kSliceStepEta = 0.5;    // η

    constexpr int    kHistNBins = 160;       // 5 MeV/bin
    constexpr double kHistMin   = 0.0;
    constexpr double kHistMax   = 0.8;

    std::string fmtEdge(double x) {
        char buf[32];
        snprintf(buf, sizeof(buf), "%.1f", x);
        std::string s(buf);
        for (auto& c : s) if (c == '.') c = 'p';
        return s;
    }
}

// -----------------------------------------------------------------------------
// Process the single sim channel: loop events, fill REC/COR/TRU histograms
// (full + per-slice), write a single file mimicking meson_research's
// sim output layout.
// -----------------------------------------------------------------------------
void processSim(const std::string& in_path,
                const std::string& out_path,
                const EcalLookup&  look)
{
    TFile* fin = TFile::Open(in_path.c_str(), "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "Cannot open " << in_path << "\n"; return;
    }
    auto* t = (TTree*)fin->Get("meson_dalitz_nt");
    if (!t) { std::cerr << "meson_dalitz_nt missing in " << in_path << "\n"; return; }

    // Activate only the branches we read.
    t->SetBranchStatus("*", 0);
    for (const char* b : {
            "m_ee",
            "m_epemg",
            "m_epemg_cor",
            "m_epemg_sim",
            "ecal_quality_pass",
            "neutr_cluster_energy",
            "neutr_cluster_theta",
            "neutr_cluster_phi",
            "oa_epem",
            "sim_genweight"
         }) t->SetBranchStatus(b, 1);

    float m_ee=0, m_eg=0, m_eg_cor=0, m_eg_sim=0,
          ecal_q=0, ne_E=0, ne_th=0, ne_ph=0, oa=0, w=1.0f;
    t->SetBranchAddress("m_ee",                 &m_ee);
    t->SetBranchAddress("m_epemg",              &m_eg);
    t->SetBranchAddress("m_epemg_cor",          &m_eg_cor);
    t->SetBranchAddress("m_epemg_sim",          &m_eg_sim);
    t->SetBranchAddress("ecal_quality_pass",    &ecal_q);
    t->SetBranchAddress("neutr_cluster_energy", &ne_E);
    t->SetBranchAddress("neutr_cluster_theta",  &ne_th);
    t->SetBranchAddress("neutr_cluster_phi",    &ne_ph);
    t->SetBranchAddress("oa_epem",              &oa);
    t->SetBranchAddress("sim_genweight",        &w);

    TFile* fout = TFile::Open(out_path.c_str(), "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "Cannot create " << out_path << "\n";
        fin->Close();  return;
    }
    fout->cd();

    // -- Allocate full histograms for the three flavours (REC/COR/TRU). ------
    auto makeFull = [&](const char* name, const char* tag) {
        TH1D* h = new TH1D(name,
            TString::Format(
                "M(e^{+}e^{-}#gamma) %s ECALCOR, OA #in [%.1f, %.1f] deg (full);"
                "M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Counts (sim_genweight)",
                tag, kSliceMin, kSliceMax),
            kHistNBins, kHistMin, kHistMax);
        h->Sumw2();
        return h;
    };
    TH1D* h_full_rec = makeFull("m_epemg_full",     "REC");
    TH1D* h_full_cor = makeFull("m_epemg_cor_full", "COR");
    TH1D* h_full_tru = makeFull("m_epemg_tru_full", "TRU");

    // -- Allocate per-slice histograms for two grids × three flavours. -------
    auto buildSliceHists = [&](double step, const char* prefix, const char* tag) {
        const int n = static_cast<int>(std::round((kSliceMax - kSliceMin) / step));
        std::vector<TH1D*>  hs(n, nullptr);
        std::vector<double> los(n), his(n);
        for (int i = 0; i < n; ++i) {
            const double lo = kSliceMin + i * step;
            const double hi = kSliceMin + (i + 1) * step;
            los[i] = lo; his[i] = hi;
            const std::string suff = "oa_" + fmtEdge(lo) + "_" + fmtEdge(hi);
            const std::string hname = std::string(prefix) + "_" + suff;
            // Avoid duplicate ROOT names if fine and coarse grids produce
            // identical edge tags (only at exact step coincidence — does not
            // happen for 0.2° vs 0.5°, but the guard is cheap).
            if (gDirectory->FindObject(hname.c_str())) continue;
            TH1D* h = new TH1D(hname.c_str(),
                TString::Format(
                    "M(e^{+}e^{-}#gamma) %s ECALCOR, OA #in [%.1f, %.1f] deg;"
                    "M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Counts (sim_genweight)",
                    tag, lo, hi),
                kHistNBins, kHistMin, kHistMax);
            h->Sumw2();
            hs[i] = h;
        }
        return std::make_tuple(hs, los, his);
    };

    // Fine grid (π⁰, 0.2°) and coarse grid (η, 0.5°), each for all flavours.
    auto [hs_pi0_rec, los_pi0, his_pi0] = buildSliceHists(kSliceStep, "m_epemg",     "REC");
    auto [hs_pi0_cor, _2,      _3]      = buildSliceHists(kSliceStep, "m_epemg_cor", "COR");
    auto [hs_pi0_tru, _4,      _5]      = buildSliceHists(kSliceStep, "m_epemg_tru", "TRU");
    (void)_2; (void)_3; (void)_4; (void)_5;
    auto [hs_eta_rec, los_eta, his_eta] = buildSliceHists(kSliceStepEta, "m_epemg",     "REC");
    auto [hs_eta_cor, _6,      _7]      = buildSliceHists(kSliceStepEta, "m_epemg_cor", "COR");
    auto [hs_eta_tru, _8,      _9]      = buildSliceHists(kSliceStepEta, "m_epemg_tru", "TRU");
    (void)los_pi0; (void)his_pi0; (void)los_eta; (void)his_eta;
    (void)_6; (void)_7; (void)_8; (void)_9;

    auto pickSlice = [&](double v, double step,
                         const std::vector<TH1D*>& hs) -> TH1D* {
        if (v < kSliceMin || v >= kSliceMax) return nullptr;
        const int i = static_cast<int>((v - kSliceMin) / step);
        if (i < 0 || i >= (int)hs.size()) return nullptr;
        return hs[i];
    };

    // -- Event loop -----------------------------------------------------------
    const Long64_t N = t->GetEntries();
    std::cout << "  " << in_path << ": " << N << " entries\n";

    Long64_t n_pass = 0, n_filled = 0;
    for (Long64_t ev = 0; ev < N; ++ev) {
        t->GetEntry(ev);
        if (ecal_q != 1.0f) continue;
        ++n_pass;

        // Lookup s and compute corrected masses for each flavour.
        // Same s for REC and COR (single calibration map). TRU bypasses
        // calorimeter and is taken straight from the simulation truth.
        const double E_GeV = ne_E / 1000.0;
        const double s     = look.s(E_GeV, ne_th, ne_ph);
        const double m_ee2 = m_ee * m_ee;

        const double m2_rec = (1.0 - s) * m_ee2 + s * (m_eg     * m_eg);
        const double m2_cor = (1.0 - s) * m_ee2 + s * (m_eg_cor * m_eg_cor);
        if (m2_rec <= 0 || m2_cor <= 0) continue;

        const double m_corr_rec = std::sqrt(m2_rec);
        const double m_corr_cor = std::sqrt(m2_cor);
        const double m_tru      = m_eg_sim;   // truth left as-is

        if (oa < kSliceMin || oa >= kSliceMax) continue;

        h_full_rec->Fill(m_corr_rec, w);
        h_full_cor->Fill(m_corr_cor, w);
        h_full_tru->Fill(m_tru,      w);

        if (TH1D* h = pickSlice(oa, kSliceStep,    hs_pi0_rec)) h->Fill(m_corr_rec, w);
        if (TH1D* h = pickSlice(oa, kSliceStep,    hs_pi0_cor)) h->Fill(m_corr_cor, w);
        if (TH1D* h = pickSlice(oa, kSliceStep,    hs_pi0_tru)) h->Fill(m_tru,      w);
        if (TH1D* h = pickSlice(oa, kSliceStepEta, hs_eta_rec)) h->Fill(m_corr_rec, w);
        if (TH1D* h = pickSlice(oa, kSliceStepEta, hs_eta_cor)) h->Fill(m_corr_cor, w);
        if (TH1D* h = pickSlice(oa, kSliceStepEta, hs_eta_tru)) h->Fill(m_tru,      w);

        ++n_filled;
    }
    std::cout << "  ecal_quality pass: " << n_pass
              << "  filled: "             << n_filled
              << "\n  full REC entries: " << (Long64_t)h_full_rec->GetEntries()
              << "  COR: "                << (Long64_t)h_full_cor->GetEntries()
              << "  TRU: "                << (Long64_t)h_full_tru->GetEntries() << "\n";

    fout->Write();
    fout->Close();
    fin->Close();

    std::cout << "  → " << out_path << "\n";
}

void apply_ecal_correction(const char* flavour = "rec",
                           const char* map_dim = "2d") {

    std::string fl = flavour;
    for (auto& c : fl) c = std::tolower(c);
    if (fl != "rec") {
        std::cerr << "Only 'rec' is supported as the calibration flavour here\n";
        return;
    }
    std::string md = map_dim ? map_dim : "2d";
    for (auto& c : md) c = std::tolower(c);

    EcalLookup look;
    std::string map_path, map_hist;
    bool loaded = false;
    if (md == "combined") {
        map_path = "ecal_combined_3dscan_sim.root";
        map_hist = "h_s_combined_3d";
        loaded = look.load3D(map_path, map_hist);
    } else if (md == "3d") {
        map_path = "ecal_pi0_3dscan_" + fl + "_sim.root";
        map_hist = "h_s_" + fl + "_3d";
        loaded = look.load3D(map_path, map_hist);
    } else {
        map_path = "ecal_pi0_2dscan_" + fl + "_sim.root";
        map_hist = "h_s_" + fl;
        loaded = look.load(map_path, map_hist);
    }
    if (!loaded) {
        std::cerr << "Cannot open " << map_path << "\n";
        return;
    }
    std::cout << "ECAL correction map loaded from " << map_path
              << " (" << md << ")\n";

    gSystem->mkdir("plots/output", kTRUE);

    const std::string sfx = (md == "combined") ? "_combined"
                          : (md == "3d")        ? "_3dmap"
                                                : "";
    processSim("../output_epem_sim.root",
               "research_sim_ecalcor" + sfx + ".root", look);

    std::cout << "\nDone. research_sim_ecalcor" << sfx << ".root produced.\n";
}
