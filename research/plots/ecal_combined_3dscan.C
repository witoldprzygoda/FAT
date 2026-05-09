// research/plots/ecal_combined_3dscan.C — 3D ECAL energy-scale calibration
// map built by COMBINING per-cell π⁰ and η peak-position fits.
//
// Per (E_γ, θ_γ, φ_γ) cell we measure TWO independent peak positions in
// the m_epemg distribution: μ_π⁰ ≈ 0.135 and μ_η ≈ 0.548. Each yields a
// scale estimate
//
//     s_X = (m_X_PDG / μ_X_fit)²       X ∈ {π⁰, η}
//
// The two are combined as an inverse-variance-weighted average:
//
//     s = (w_π⁰ s_π⁰ + w_η s_η) / (w_π⁰ + w_η),    w_X = 1/σ_X²
//
// where σ_X is propagated from the μ-fit error. This gives a single,
// physically-well-defined per-cell scale. Cells where only one of the
// two fits succeeds use that one alone. Cells where neither succeeds
// stay empty (handled exactly as before — sector-aware neighbour fill,
// dead sectors skipped).
//
// Why combine: π⁰ Dalitz photons sample lower E_γ (peak ~0.4 GeV);
// η Dalitz photons sample higher E_γ (peak ~0.8 GeV). Within a shared
// cell both peaks are present (with different statistical weights), but
// each meson constrains the ECAL response at a slightly different
// energy regime — together they pin down the calibration over a wider
// dynamic range than either alone. In particular η helps in high-E_γ
// cells where π⁰ stats are sparse.
//
// φ axis: full [0°, 360°), per-sector, dead sector kept empty.
// Pass-2 / pass-3 (neighbour fill / median smoothing) restricted to
// same sector — same logic as ecal_pi0_3dscan.C.
//
// Usage (from research/):
//   root -l -b -q 'plots/ecal_combined_3dscan.C("exp")'   // exp data
//   root -l -b -q 'plots/ecal_combined_3dscan.C("sim")'   // sim
//
// Output:
//   ecal_combined_3dscan_<dataset>.root  with TH3D h_s_combined_3d
//   plots/output/ecal_combined_3dscan_<dataset>_phi_pattern.{pdf,png}

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
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

    // m_epemg histogram binning (for both peak fits).
    constexpr int    kHistNBins = 200;
    constexpr double kHistMin   = 0.0;
    constexpr double kHistMax   = 0.8;

    // Apex windows for the two peaks. Tight enough that a Gauss + linear
    // baseline captures only the peak position. The exact window center
    // does not matter — we'll find the local maximum first.
    constexpr double kPi0ApexLo = 0.110, kPi0ApexHi = 0.160;
    constexpr double kEtaApexLo = 0.510, kEtaApexHi = 0.590;

    // 3D grid — same as ecal_pi0_3dscan with full φ.
    // E_γ edges chosen denser at low E where π⁰ Dalitz photons dominate
    // (~0.2–0.6 GeV is ~95% of statistics) and stop at 1.5 GeV (above
    // that the sample is essentially empty after ncells==1 / OA>4 cuts).
    // 10 bins, narrower below 0.5 GeV.
    const std::vector<double> kE_edges     = {0.20, 0.26, 0.33, 0.40, 0.47,
                                              0.55, 0.65, 0.78, 0.95, 1.20,
                                              1.50};
    // θ_γ edges: 4° steps from 15° to 47°. The HADES ECAL acceptance
    // ends below ~47° in this sample (no events between 47° and 55°),
    // so the previous {15, 28, 42, 55} grid wasted the upper bin. 8 bins.
    const std::vector<double> kTheta_edges = {15.0, 19.0, 23.0, 27.0, 31.0,
                                              35.0, 39.0, 43.0, 47.0};
    std::vector<double> makePhiEdges() {
        std::vector<double> v; v.reserve(37);
        for (int i = 0; i <= 36; ++i) v.push_back(i * 10.0);
        return v;
    }
    const std::vector<double> kPhi_edges = makePhiEdges();

    constexpr int    kPhiPerSector  = 6;     // 36 phi bins / 6 sectors
    constexpr double kMinSignalApex = 50.0;  // events in apex window for fit

    const char* tag_name(int ei, int ti, int pi) {
        static thread_local char buf[32];
        snprintf(buf, sizeof(buf), "e%02d_t%02d_p%02d", ei, ti, pi);
        return buf;
    }
}

// OA(e+e-) > 4° cut — calibration uses only well-resolved dilepton pairs:
// removes γ→e+e- conversions and bremsstrahlung pairs (both populate
// tight OA), and makes π⁰/η kinematics more comparable inside each cell.
constexpr double kCalibOAMin = 4.0;

// Photon-quality cuts on the cell sample (applied to BOTH exp & sim):
//   neutr_mult            == 1   — exactly one ECAL cluster (pure
//                                  one-photon hypothesis)
//   neutr_cluster_ncells  == 1   — single-cell cluster → tightest
//                                  shower-shape bin, minimal cluster
//                                  reconstruction systematics. Removes
//                                  ~63% of events but keeps ~2.8 M (exp)
//                                  events with the cleanest photon
//                                  energy/angle estimate.

// Build CB-subtracted (or single-channel weighted) signal for one cell.
// EXP: 3 channels, CB = 2√(N₊₊·N₋₋), signal = all − CB.
// SIM: single channel, weighted by sim_genweight.
TH1D* buildSignal3D(TTree* t_all, TTree* t_pp, TTree* t_mm,
                    bool is_sim,
                    const std::string& mass_var,
                    double e_lo, double e_hi,
                    double th_lo, double th_hi,
                    double ph_lo, double ph_hi,
                    const std::string& tag)
{
    if (is_sim) {
        TH1D* h = new TH1D(("h_sig_" + tag).c_str(),
                           ";M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Counts (sim_genweight)",
                           kHistNBins, kHistMin, kHistMax);
        h->Sumw2();
        TString cut = TString::Format(
            "(ecal_quality_pass==1 && oa_epem>%g"
            " && neutr_mult==1 && neutr_cluster_ncells==1"
            " && (neutr_cluster_energy/1000.0)>=%g"
            " && (neutr_cluster_energy/1000.0)<%g"
            " && neutr_cluster_theta>=%g && neutr_cluster_theta<%g"
            " && neutr_cluster_phi>=%g && neutr_cluster_phi<%g)*sim_genweight",
            kCalibOAMin, e_lo, e_hi, th_lo, th_hi, ph_lo, ph_hi);
        TString cmd = TString::Format("%s>>h_sig_%s",
                                      mass_var.c_str(), tag.c_str());
        t_all->Draw(cmd, cut, "goff");
        h->SetDirectory(nullptr);
        return h;
    }

    auto draw = [&](TTree* t, const std::string& hname) -> TH1D* {
        TH1D* h = new TH1D(hname.c_str(),
                           ";M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Counts",
                           kHistNBins, kHistMin, kHistMax);
        h->Sumw2();
        TString cut = TString::Format(
            "ecal_quality_pass==1 && oa_epem>%g"
            " && neutr_mult==1 && neutr_cluster_ncells==1"
            " && (neutr_cluster_energy/1000.0)>=%g"
            " && (neutr_cluster_energy/1000.0)<%g"
            " && neutr_cluster_theta>=%g && neutr_cluster_theta<%g"
            " && neutr_cluster_phi>=%g && neutr_cluster_phi<%g",
            kCalibOAMin, e_lo, e_hi, th_lo, th_hi, ph_lo, ph_hi);
        TString cmd = TString::Format("%s>>%s",
                                      mass_var.c_str(), hname.c_str());
        t->Draw(cmd, cut, "goff");
        h->SetDirectory(nullptr);
        return h;
    };
    TH1D* h_all = draw(t_all, "h_all_" + tag);
    TH1D* h_pp  = draw(t_pp,  "h_pp_"  + tag);
    TH1D* h_mm  = draw(t_mm,  "h_mm_"  + tag);
    TH1D* h_cb = (TH1D*)h_all->Clone(("h_cb_" + tag).c_str());
    h_cb->Reset();  h_cb->SetDirectory(nullptr);
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
    delete h_all; delete h_pp; delete h_mm; delete h_cb;
    return h_sig;
}

// Local Gauss + linear baseline fit on a tight apex window. Same recipe
// as in ecal_pi0_3dscan.C — finds peak μ robustly without modelling
// background details.
struct PeakFit {
    double mu = 0, mu_err = 0;
    double sigma = 0;
    double yield = 0;     // apex window integral (proxy for stats)
    bool   ok = false;
};

PeakFit findPeak(TH1D* h, double apex_lo, double apex_hi,
                 const std::string& base)
{
    PeakFit r;
    if (!h || h->GetEntries() == 0) return r;
    const int b_lo = h->FindBin(apex_lo);
    const int b_hi = h->FindBin(apex_hi);
    if (b_hi <= b_lo) return r;
    const double apex_sum = h->Integral(b_lo, b_hi);
    r.yield = apex_sum;
    if (apex_sum < kMinSignalApex) return r;

    int    pb = b_lo;
    double pv = h->GetBinContent(b_lo);
    for (int b = b_lo; b <= b_hi; ++b)
        if (h->GetBinContent(b) > pv) { pv = h->GetBinContent(b); pb = b; }
    if (pv <= 0) return r;
    const double mu_seed = h->GetBinCenter(pb);

    // Tight ±halfwidth window — narrow for π⁰, wider for η (broader peak).
    const double halfw = (apex_lo > 0.4) ? 0.025 : 0.012;
    const double fit_lo = std::max(apex_lo, mu_seed - halfw);
    const double fit_hi = std::min(apex_hi, mu_seed + halfw);

    auto* fit = new TF1(("apex_" + base).c_str(),
                        "[0]*exp(-0.5*((x-[1])/[2])^2) + [3] + [4]*x",
                        fit_lo, fit_hi);
    fit->SetParameter(0, pv);
    fit->SetParameter(1, mu_seed);
    fit->SetParameter(2, halfw * 0.5);
    fit->SetParameter(3, 0.0);
    fit->SetParameter(4, 0.0);
    fit->SetParLimits(0, 0.0, 1e12);
    fit->SetParLimits(1, fit_lo, fit_hi);
    fit->SetParLimits(2, 0.003, halfw * 1.2);

    int status = h->Fit(fit, "RQN");
    if (status != 0) {
        fit->SetParameter(2, halfw * 0.7);
        h->Fit(fit, "RQN");
    }
    r.mu     = fit->GetParameter(1);
    r.mu_err = fit->GetParError(1);
    r.sigma  = fit->GetParameter(2);
    delete fit;

    r.ok = (r.mu > fit_lo + 1e-4 && r.mu < fit_hi - 1e-4
            && r.mu_err > 0.0 && r.mu_err < 0.005
            && r.sigma > 0.003);
    return r;
}

// Combine π⁰ and η scale measurements via inverse-variance weighting.
// Returns NaN if neither fit is OK.
struct CombFit {
    double s = 0;        // combined scale
    double s_err = 0;    // combined uncertainty
    bool   has_pi0 = false, has_eta = false;
    double s_pi0 = 0, s_eta = 0;
};

CombFit combineScales(const PeakFit& fp, const PeakFit& fe) {
    CombFit r;
    auto sigma_s = [](const PeakFit& f, double m_pdg) {
        // s = (m_pdg/μ)²  →  σ_s = 2·m_pdg²/μ³ · σ_μ = 2 s · σ_μ / μ
        const double s = (m_pdg / f.mu) * (m_pdg / f.mu);
        return 2.0 * s * f.mu_err / f.mu;
    };
    double sum_w = 0, sum_ws = 0;
    if (fp.ok) {
        const double s = (kPi0PDG / fp.mu) * (kPi0PDG / fp.mu);
        const double sig = sigma_s(fp, kPi0PDG);
        if (sig > 0) {
            const double w = 1.0 / (sig * sig);
            sum_w  += w;
            sum_ws += w * s;
            r.has_pi0 = true;
            r.s_pi0   = s;
        }
    }
    if (fe.ok) {
        const double s = (kEtaPDG / fe.mu) * (kEtaPDG / fe.mu);
        const double sig = sigma_s(fe, kEtaPDG);
        if (sig > 0) {
            const double w = 1.0 / (sig * sig);
            sum_w  += w;
            sum_ws += w * s;
            r.has_eta = true;
            r.s_eta   = s;
        }
    }
    if (sum_w <= 0) return r;
    r.s     = sum_ws / sum_w;
    r.s_err = 1.0 / std::sqrt(sum_w);
    return r;
}

void ecal_combined_3dscan(const char* dataset = "exp") {

    const std::string ds = dataset;
    const bool is_sim = (ds == "sim" || ds == "SIM");

    TFile *f_all = nullptr, *f_pp = nullptr, *f_mm = nullptr;
    if (is_sim) {
        f_all = TFile::Open("../output_epem_sim.root", "READ");
        if (!f_all || f_all->IsZombie()) { std::cerr << "Cannot open sim\n"; return; }
    } else {
        f_all = TFile::Open("../output_epem_exp.root", "READ");
        f_pp  = TFile::Open("../output_epep_exp.root", "READ");
        f_mm  = TFile::Open("../output_emem_exp.root", "READ");
        if (!f_all || f_all->IsZombie()) { std::cerr << "Cannot open epem exp\n"; return; }
        if (!f_pp  || f_pp ->IsZombie()) { std::cerr << "Cannot open epep exp\n"; return; }
        if (!f_mm  || f_mm ->IsZombie()) { std::cerr << "Cannot open emem exp\n"; return; }
    }
    auto* t_all = (TTree*)f_all->Get("meson_dalitz_nt");
    auto* t_pp  = is_sim ? nullptr : (TTree*)f_pp ->Get("meson_dalitz_nt");
    auto* t_mm  = is_sim ? nullptr : (TTree*)f_mm ->Get("meson_dalitz_nt");
    if (!t_all || (!is_sim && (!t_pp || !t_mm))) {
        std::cerr << "meson_dalitz_nt missing\n"; return;
    }

    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);
    gSystem->mkdir("plots/output", kTRUE);

    const int nE  = (int)kE_edges.size()     - 1;
    const int nTh = (int)kTheta_edges.size() - 1;
    const int nPh = (int)kPhi_edges.size()   - 1;
    const int n_total_cells = nE * nTh * nPh;

    // Per-cell results: combined s, error, plus separate π⁰/η for diagnostics.
    struct Cell {
        double s = 0, s_err = 0;
        bool   ok = false;
        bool   has_pi0 = false, has_eta = false;
        double s_pi0 = 0, s_eta = 0;
        double apex_pi0 = 0, apex_eta = 0;
    };
    std::vector<std::vector<std::vector<Cell>>> cell(
        nE, std::vector<std::vector<Cell>>(nTh, std::vector<Cell>(nPh)));

    std::cout << "Combined π⁰+η 3D scan — " << ds
              << "  cells=" << n_total_cells
              << " (" << nE << "×" << nTh << "×" << nPh << ")\n";

    int n_cells_done = 0;
    for (int ei = 0; ei < nE; ++ei) {
        for (int ti = 0; ti < nTh; ++ti) {
            for (int pi = 0; pi < nPh; ++pi) {
                const double el = kE_edges[ei],     eh = kE_edges[ei + 1];
                const double tl = kTheta_edges[ti], th = kTheta_edges[ti + 1];
                const double pl = kPhi_edges[pi],   ph = kPhi_edges[pi + 1];
                const std::string tag = tag_name(ei, ti, pi);

                TH1D* h_sig = buildSignal3D(t_all, t_pp, t_mm, is_sim, "m_epemg",
                                            el, eh, tl, th, pl, ph, tag);

                PeakFit fp = findPeak(h_sig, kPi0ApexLo, kPi0ApexHi, "p0_" + tag);
                PeakFit fe = findPeak(h_sig, kEtaApexLo, kEtaApexHi, "et_" + tag);

                CombFit cf = combineScales(fp, fe);
                Cell& c = cell[ei][ti][pi];
                c.has_pi0 = cf.has_pi0;
                c.has_eta = cf.has_eta;
                c.s_pi0   = cf.s_pi0;
                c.s_eta   = cf.s_eta;
                c.s       = cf.s;
                c.s_err   = cf.s_err;
                c.ok      = (c.has_pi0 || c.has_eta);
                c.apex_pi0 = fp.yield;
                c.apex_eta = fe.yield;

                delete h_sig;
                ++n_cells_done;
                if (n_cells_done % 50 == 0)
                    std::cout << "  done " << n_cells_done << " / "
                              << n_total_cells << "\n";
            }
        }
    }

    // Stats summary.
    int n_pi0_only = 0, n_eta_only = 0, n_both = 0, n_neither = 0;
    for (int ei = 0; ei < nE; ++ei)
      for (int ti = 0; ti < nTh; ++ti)
        for (int pi = 0; pi < nPh; ++pi) {
            const Cell& c = cell[ei][ti][pi];
            if      ( c.has_pi0 &&  c.has_eta) ++n_both;
            else if ( c.has_pi0 && !c.has_eta) ++n_pi0_only;
            else if (!c.has_pi0 &&  c.has_eta) ++n_eta_only;
            else                                ++n_neither;
        }
    std::cout << "\n  cells with both fits: " << n_both
              << "  π⁰-only: "                << n_pi0_only
              << "  η-only: "                 << n_eta_only
              << "  neither: "                << n_neither << "\n";

    // ---- Pass 2 / Pass 3 (sector-aware fill + median smoothing) ----------
    // Identify dead sectors first.
    std::vector<bool> dead_sector(nPh / kPhiPerSector, false);
    for (int sec = 0; sec < (int)dead_sector.size(); ++sec) {
        int n_cells = 0, n_low = 0;
        for (int ei = 0; ei < nE; ++ei)
          for (int ti = 0; ti < nTh; ++ti)
            for (int dp = 0; dp < kPhiPerSector; ++dp) {
                const int pi = sec * kPhiPerSector + dp;
                ++n_cells;
                const Cell& c = cell[ei][ti][pi];
                if (!c.ok && c.apex_pi0 < kMinSignalApex
                          && c.apex_eta < kMinSignalApex) ++n_low;
            }
        if (n_cells > 0 && (double)n_low / n_cells > 0.90) {
            dead_sector[sec] = true;
            std::cout << "  DEAD SECTOR " << sec
                      << " (φ ∈ [" << sec * 60 << ", " << (sec + 1) * 60
                      << ")°) — skipping neighbour fill.\n";
        }
    }

    // Pass-2 fill: weighted average of OK same-sector neighbours.
    auto neighbourS = [&](int ei, int ti, int pi) -> double {
        const int my_sec = pi / kPhiPerSector;
        double sum_w = 0, sum_ws = 0;
        for (int de = -1; de <= 1; ++de)
          for (int dt = -1; dt <= 1; ++dt)
            for (int dp = -1; dp <= 1; ++dp) {
                if (de == 0 && dt == 0 && dp == 0) continue;
                const int ne = ei + de;
                const int nt = ti + dt;
                const int np = pi + dp;
                if (ne < 0 || ne >= nE)  continue;
                if (nt < 0 || nt >= nTh) continue;
                if (np < 0 || np >= nPh) continue;
                if ((np / kPhiPerSector) != my_sec) continue;
                const Cell& nc = cell[ne][nt][np];
                if (!nc.ok || nc.s_err <= 0) continue;
                const double w = 1.0 / (nc.s_err * nc.s_err);
                sum_w  += w;
                sum_ws += w * nc.s;
            }
        return (sum_w > 0) ? sum_ws / sum_w : -1.0;
    };

    int n_recovered = 0;
    for (int ei = 0; ei < nE; ++ei)
      for (int ti = 0; ti < nTh; ++ti)
        for (int pi = 0; pi < nPh; ++pi) {
            Cell& c = cell[ei][ti][pi];
            if (c.ok) continue;
            const int sec = pi / kPhiPerSector;
            if (dead_sector[sec]) continue;
            const double s_est = neighbourS(ei, ti, pi);
            if (s_est <= 0) continue;
            c.s     = s_est;
            c.s_err = 0.05;     // conservative
            c.ok    = true;
            ++n_recovered;
        }
    std::cout << "Pass-2 filled " << n_recovered
              << " cells from same-sector neighbours.\n";

    // Pass-3: median smoothing (same-sector ring), only outliers.
    auto angularMedian = [&](int ei, int ti, int pi) -> double {
        const int my_sec = pi / kPhiPerSector;
        std::vector<double> vs;
        for (int dt = -1; dt <= 1; ++dt)
          for (int dp = -1; dp <= 1; ++dp) {
            if (dt == 0 && dp == 0) continue;
            const int nt = ti + dt;
            const int np = pi + dp;
            if (nt < 0 || nt >= nTh) continue;
            if (np < 0 || np >= nPh) continue;
            if ((np / kPhiPerSector) != my_sec) continue;
            const Cell& nc = cell[ei][nt][np];
            if (nc.ok) vs.push_back(nc.s);
          }
        if (vs.size() < 3) return -1.0;
        std::nth_element(vs.begin(), vs.begin() + vs.size() / 2, vs.end());
        return vs[vs.size() / 2];
    };

    int n_smoothed = 0;
    for (int ei = 0; ei < nE; ++ei)
      for (int ti = 0; ti < nTh; ++ti)
        for (int pi = 0; pi < nPh; ++pi) {
            Cell& c = cell[ei][ti][pi];
            if (!c.ok) continue;
            const double med = angularMedian(ei, ti, pi);
            if (med <= 0) continue;
            // Outlier if |s − median| > 0.04 AND own error large
            if (std::fabs(c.s - med) > 0.04 && c.s_err > 0.02) {
                c.s = med;
                ++n_smoothed;
            }
        }
    std::cout << "Pass-3 smoothed " << n_smoothed
              << " outlier cells.\n";

    // ---- Output: TH3D + diagnostic phi-pattern projection ---------------
    auto* h_s = new TH3D("h_s_combined_3d",
        TString::Format("Combined π⁰+η scale s(E, θ, φ) — %s;"
                        "E_{#gamma} [GeV];#theta_{#gamma} [deg];#phi_{#gamma} [deg]",
                        ds.c_str()),
        nE,  kE_edges.data(),
        nTh, kTheta_edges.data(),
        nPh, kPhi_edges.data());
    h_s->SetDirectory(nullptr);

    int n_filled = 0;
    for (int ei = 0; ei < nE; ++ei)
      for (int ti = 0; ti < nTh; ++ti)
        for (int pi = 0; pi < nPh; ++pi) {
            const Cell& c = cell[ei][ti][pi];
            if (c.ok && c.s > 0) {
                h_s->SetBinContent(ei + 1, ti + 1, pi + 1, c.s);
                ++n_filled;
            }
        }
    std::cout << "  cells filled in TH3D: " << n_filled << " / "
              << n_total_cells << " (" << (100.0 * n_filled / n_total_cells)
              << "%)\n";

    // Phi-pattern projection (median over (E, θ) per phi-bin).
    TH1D h_phi("h_s_vs_phi", ";#phi_{#gamma} [deg];median s(combined)",
               nPh, kPhi_edges.data());
    for (int pi = 0; pi < nPh; ++pi) {
        std::vector<double> vs;
        for (int ei = 0; ei < nE; ++ei)
          for (int ti = 0; ti < nTh; ++ti) {
              const Cell& c = cell[ei][ti][pi];
              if (c.ok && c.s > 0) vs.push_back(c.s);
          }
        if (vs.size() >= 2) {
            std::nth_element(vs.begin(), vs.begin() + vs.size() / 2, vs.end());
            h_phi.SetBinContent(pi + 1, vs[vs.size() / 2]);
        }
    }
    TCanvas c("c_phi", "phi pattern", 1300, 600);
    c.SetGrid();  c.SetMargin(0.10, 0.05, 0.13, 0.08);
    h_phi.SetMinimum(0.85); h_phi.SetMaximum(1.15);
    h_phi.SetMarkerStyle(20); h_phi.SetMarkerSize(1.0);
    h_phi.SetMarkerColor(kBlack); h_phi.SetLineColor(kBlack);
    h_phi.Draw("P");
    for (int s = 1; s < 6; ++s) {
        auto* l = new TLine(s * 60.0, 0.85, s * 60.0, 1.15);
        l->SetLineStyle(3); l->SetLineColor(kGray + 2); l->Draw();
    }
    auto* l1 = new TLine(0, 1.0, 360, 1.0);
    l1->SetLineStyle(2); l1->SetLineColor(kBlue); l1->SetLineWidth(2);
    l1->Draw();
    const std::string out_phi = "plots/output/ecal_combined_3dscan_"
                              + ds + "_phi_pattern";
    c.SaveAs((out_phi + ".pdf").c_str());
    c.SaveAs((out_phi + ".png").c_str());

    const std::string out_root = "ecal_combined_3dscan_" + ds + ".root";
    TFile* fout = TFile::Open(out_root.c_str(), "RECREATE");
    h_s->Write();
    fout->Close();
    f_all->Close();
    if (f_pp) f_pp->Close();
    if (f_mm) f_mm->Close();

    std::cout << "\nWrote:\n"
              << "  " << out_root << "  (TH3D h_s_combined_3d)\n"
              << "  " << out_phi  << ".{pdf,png}\n";
}
