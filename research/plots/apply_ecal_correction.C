// research/plots/apply_ecal_correction.C — produce research_*_ecalcor.root.
//
// Goal: bring the post-hoc ECAL energy-scale correction into the standard
// pipeline. Output mirrors `meson_research`'s layout (full + OA-sliced
// histograms in two slice grids), so the existing downstream tooling
// (slices_pi0.C, fit_pi0.C, yields_pi0.C, fit_eta.C, yields_eta.C, …) can
// read the corrected mass with a single path swap.
//
// What it does, per channel ∈ {epem, epep, emem}:
//   1) read meson_dalitz_nt (event by event)
//   2) for each event passing ecal_quality_pass==1, look up
//        s = s(E_γ, θ_γ)   from the 2D calibration map h_s_<fl>
//        m²_corr = (1 − s) · m_ee² + s · m_epemg²        (REC flavour)
//      The closed-form above is equivalent to
//        m²_corr = m_ee² + 2·s·E_γ·D
//      and avoids needing the per-event gamma_D branch.
//   3) fill OA-sliced histograms with m_corr using the same naming as
//      meson_research:
//        m_epemg_full
//        m_epemg_oa_<lo>_<hi>      (fine step 0.2° → 75 slices, for π⁰)
//        m_epemg_oa_<lo>_<hi>      (coarse step 0.5° → 30 slices, for η)
//   4) write research_<channel>_ecalcor.root.
//
// Usage (from research/):
//   root -l -b -q plots/apply_ecal_correction.C
//
// Pre-requisite: ecal_pi0_2dscan_rec.root must exist (run ecal_pi0_2dscan.C
// first to build the calibration map).

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
#include <vector>

namespace {
    // OA slice grids — exact match to meson_research.cc / SliceConfig.
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
// Process a single channel: loop events, fill all histograms, write file.
// -----------------------------------------------------------------------------
void processChannel(const std::string& in_path,
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
            "ecal_quality_pass",
            "neutr_cluster_energy",
            "neutr_cluster_theta",
            "oa_epem"
         }) t->SetBranchStatus(b, 1);

    float m_ee=0, m_eg=0, ecal_q=0, ne_E=0, ne_th=0, oa=0;
    t->SetBranchAddress("m_ee",                 &m_ee);
    t->SetBranchAddress("m_epemg",              &m_eg);
    t->SetBranchAddress("ecal_quality_pass",    &ecal_q);
    t->SetBranchAddress("neutr_cluster_energy", &ne_E);
    t->SetBranchAddress("neutr_cluster_theta",  &ne_th);
    t->SetBranchAddress("oa_epem",              &oa);

    TFile* fout = TFile::Open(out_path.c_str(), "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "Cannot create " << out_path << "\n";
        fin->Close();  return;
    }
    fout->cd();

    // -- Allocate output histograms in the exact meson_research naming. -------
    TH1D* h_full = new TH1D("m_epemg_full",
        TString::Format(
            "M(e^{+}e^{-}#gamma) ECALCOR, OA(e^{+}e^{-}) #in [%.1f, %.1f] deg (full);"
            "M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Counts",
            kSliceMin, kSliceMax),
        kHistNBins, kHistMin, kHistMax);
    h_full->Sumw2();

    auto buildSliceHists = [&](double step) {
        const int n = static_cast<int>(std::round((kSliceMax - kSliceMin) / step));
        std::vector<TH1D*> hs(n, nullptr);
        std::vector<double> los(n), his(n);
        for (int i = 0; i < n; ++i) {
            const double lo = kSliceMin + i * step;
            const double hi = kSliceMin + (i + 1) * step;
            los[i] = lo; his[i] = hi;
            const std::string suff = "oa_" + fmtEdge(lo) + "_" + fmtEdge(hi);
            const std::string hname = "m_epemg_" + suff;
            const TString title = TString::Format(
                "M(e^{+}e^{-}#gamma) ECALCOR, OA(e^{+}e^{-}) #in [%.1f, %.1f] deg;"
                "M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Counts",
                lo, hi);
            // ROOT throws on duplicate names — avoid by checking if a hist of
            // this name already lives in gDirectory (fine vs coarse grids
            // overlap on slice 0.0–0.5 only when the step coincides).
            if (gDirectory->FindObject(hname.c_str())) continue;
            TH1D* h = new TH1D(hname.c_str(), title,
                               kHistNBins, kHistMin, kHistMax);
            h->Sumw2();
            hs[i] = h;
        }
        return std::make_tuple(hs, los, his);
    };

    auto [hs_pi0,  los_pi0,  his_pi0]  = buildSliceHists(kSliceStep);
    auto [hs_eta,  los_eta,  his_eta]  = buildSliceHists(kSliceStepEta);

    // Helper: bin lookup for variable-width grid via integer math.
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

        // Lookup s and compute corrected mass squared:
        //   m²_corr = (1 − s) · m_ee² + s · m²
        const double E_GeV = ne_E / 1000.0;
        const double s     = look.s(E_GeV, ne_th);
        const double m2    = m_eg * m_eg;
        const double m_ee2 = m_ee * m_ee;
        const double m2c   = (1.0 - s) * m_ee2 + s * m2;
        if (m2c <= 0.0) continue;
        const double m_corr = std::sqrt(m2c);

        if (oa >= kSliceMin && oa < kSliceMax) {
            h_full->Fill(m_corr);
            if (TH1D* h = pickSlice(oa, kSliceStep,    hs_pi0)) h->Fill(m_corr);
            if (TH1D* h = pickSlice(oa, kSliceStepEta, hs_eta)) h->Fill(m_corr);
            ++n_filled;
        }
    }
    std::cout << "  ecal_quality pass: " << n_pass
              << "  filled: "             << n_filled
              << "  full hist entries: "  << (Long64_t)h_full->GetEntries() << "\n";

    fout->Write();
    fout->Close();
    fin->Close();

    std::cout << "  → " << out_path << "\n";
}

void apply_ecal_correction(const char* flavour = "rec") {

    std::string fl = flavour;
    for (auto& c : fl) c = std::tolower(c);
    if (fl != "rec") {
        std::cerr << "Only 'rec' is supported here. The COR flavour would need\n"
                  << "m_ee_cor in the ntuple to apply the same closed-form\n"
                  << "correction; consider extending main.cc if needed.\n";
        return;
    }

    EcalLookup look;
    const std::string map_path = "ecal_pi0_2dscan_" + fl + ".root";
    if (!look.load(map_path, "h_s_" + fl)) {
        std::cerr << "Run ecal_pi0_2dscan.C first to produce " << map_path << "\n";
        return;
    }
    std::cout << "ECAL correction map loaded from " << map_path << "\n";

    gSystem->mkdir("plots/output", kTRUE);

    processChannel("../output_epem_exp.root", "research_epem_ecalcor.root", look);
    processChannel("../output_epep_exp.root", "research_epep_ecalcor.root", look);
    processChannel("../output_emem_exp.root", "research_emem_ecalcor.root", look);

    std::cout << "\nDone. Three files produced — same layout as research_<ch>.root,\n"
                 "but with the ECAL energy-scale correction applied event-by-event.\n";
}
