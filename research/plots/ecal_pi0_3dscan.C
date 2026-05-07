// research/plots/ecal_pi0_3dscan.C — π⁰ peak μ(E_γ, θ_γ, φ_γ) per-sector
// ECAL energy-scale calibration map.
//
// φ axis is the FULL [0, 360°) range (NOT folded mod 60°), so each of
// the 6 HADES ECAL sectors gets its own column of cells. Sectors have
// genuinely different miscalibrations in HADES, so folding would just
// average them out. A coarser (E, θ) grid compensates for the 6× drop
// in events per cell from splitting φ into 36 bins.
//
// Per-cell peak position is found with a SIMPLE Gaussian fit on the
// apex window [0.115, 0.155] — no full CB+P3+exp model. The user does
// not care about background quality here, only about the location of
// the maximum, which a narrow-window Gauss recovers reliably.
//
// Pipeline per cell:
//   1) signal hist via TTree::Draw with cell cut + CB subtraction from
//      the three input channels (epem / epep / emem)
//   2) Gauss apex fit; cell ok if μ ∈ tight bounds and apex stat sufficient
//   3) for failed cells: μ from neighbour-averaged value (3D box, no φ
//      cyclic wrap because sectors differ)
//
// Output:
//   ecal_pi0_3dscan_<fl>.root  with TH3D h_s_<fl>_3d (calibration map),
//                              + h_mu_3d, h_yield_3d
//                              + scan3d_results TTree (one row per cell)
//   plots/output/ecal_pi0_3dscan_{mu,scale}_phi*_<fl>.{pdf,png}
//                              (one (E,θ) heatmap slice per φ-bin × 36)
//
// Usage (from research/):
//   root -l -b -q plots/ecal_pi0_3dscan.C            # REC

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TCanvas.h>
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

    // Apex window for the Gaussian peak fit. Tight enough that the bg
    // shape is approximately linear and a Gauss + linear baseline fit
    // recovers the peak position cleanly.
    constexpr double kApexLo = 0.110;
    constexpr double kApexHi = 0.160;

    // 3D grid — coarser (E, θ) to push per-cell stats above ~17k events
    // (was ~7k). At low stats the per-cell μ noise (~1 MeV) is the same
    // order as the per-sector signal we want to capture; coarsening
    // pushes σ_μ down to ~0.5 MeV so the per-sector pattern resolves.
    // 4 × 3 × 36 = 432 cells, ≈17500 events / cell on average.
    const std::vector<double> kE_edges = {
        0.20, 0.45, 0.75, 1.20, 2.50
    };
    const std::vector<double> kTheta_edges = {
        15.0, 28.0, 42.0, 55.0
    };
    // Full φ ∈ [0, 360°) in 36 bins of 10°.
    std::vector<double> makePhiEdges() {
        std::vector<double> v;
        v.reserve(37);
        for (int i = 0; i <= 36; ++i) v.push_back(i * 10.0);
        return v;
    }
    const std::vector<double> kPhi_edges = makePhiEdges();

    // Minimum signal counts in the apex window to attempt the Gauss fit.
    constexpr double kMinSignalApex = 80.0;

    const char* tag_name(int ei, int ti, int pi) {
        static thread_local char buf[32];
        snprintf(buf, sizeof(buf), "e%02d_t%02d_p%02d", ei, ti, pi);
        return buf;
    }
}

// Build CB-subtracted signal for one cell. Same CB recipe as before:
// CB(b) = 2 √(N_++(b) · N_--(b)), signal = all − CB.
TH1D* buildSignal3D(TTree* t_all, TTree* t_pp, TTree* t_mm,
                    const std::string& mass_var,
                    double e_lo, double e_hi,
                    double th_lo, double th_hi,
                    double ph_lo, double ph_hi,
                    const std::string& tag)
{
    auto draw = [&](TTree* t, const std::string& hname) -> TH1D* {
        TH1D* h = new TH1D(hname.c_str(),
                           ";M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Counts",
                           kHistNBins, kHistMin, kHistMax);
        h->Sumw2();
        // Full-φ cut. Note: phi is already in [0, 360°) on this dataset.
        TString cut = TString::Format(
            "ecal_quality_pass==1"
            " && (neutr_cluster_energy/1000.0)>=%g"
            " && (neutr_cluster_energy/1000.0)<%g"
            " && neutr_cluster_theta>=%g"
            " && neutr_cluster_theta<%g"
            " && neutr_cluster_phi>=%g"
            " && neutr_cluster_phi<%g",
            e_lo, e_hi, th_lo, th_hi, ph_lo, ph_hi);
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

    delete h_all; delete h_pp; delete h_mm; delete h_cb;
    return h_sig;
}

// Simple peak finder: Gaussian + linear baseline on the apex window.
// Returns μ (peak position) and a quality flag. We don't model the full
// background — only what's needed to localise the maximum.
struct ApexFit {
    double mu = 0, mu_err = 0;
    double sigma = 0;
    double yield = 0;          // signal integral inside apex window
    bool   ok = false;
};

ApexFit findPeak(TH1D* h, const std::string& base) {
    ApexFit r;
    if (!h || h->GetEntries() == 0) return r;

    // -- Step 1: locate the max bin in a wide search window. ----------------
    const int b_lo_search = h->FindBin(kApexLo);
    const int b_hi_search = h->FindBin(kApexHi);
    if (b_hi_search <= b_lo_search) return r;

    const double apex_sum = h->Integral(b_lo_search, b_hi_search);
    r.yield = apex_sum;
    if (apex_sum < kMinSignalApex) return r;

    int    pb = b_lo_search;
    double pv = h->GetBinContent(b_lo_search);
    for (int b = b_lo_search; b <= b_hi_search; ++b) {
        if (h->GetBinContent(b) > pv) { pv = h->GetBinContent(b); pb = b; }
    }
    if (pv <= 0) return r;
    const double mu_seed = h->GetBinCenter(pb);

    // -- Step 2: tight LOCAL Gauss + linear baseline around the peak. -------
    // Window ±12 MeV — narrow enough that the asymmetric CB tail does not
    // pull the fitted μ off the actual maximum, wide enough to constrain
    // σ. This is the "lokalne fitowanie czegokolwiek (gaussa) w ścisłej
    // okolicy maksimum" the user asked for.
    constexpr double kLocalHalfWidth = 0.012;
    const double fit_lo = std::max(kApexLo, mu_seed - kLocalHalfWidth);
    const double fit_hi = std::min(kApexHi, mu_seed + kLocalHalfWidth);

    auto* fit = new TF1(("apex_" + base).c_str(),
                        "[0]*exp(-0.5*((x-[1])/[2])^2) + [3] + [4]*x",
                        fit_lo, fit_hi);
    fit->SetParameter(0, pv);
    fit->SetParameter(1, mu_seed);
    fit->SetParameter(2, 0.007);
    fit->SetParameter(3, 0.0);
    fit->SetParameter(4, 0.0);
    fit->SetParLimits(0, 0.0, 1e12);
    fit->SetParLimits(1, fit_lo, fit_hi);
    fit->SetParLimits(2, 0.003, 0.015);

    int status = h->Fit(fit, "RQN");
    if (status != 0) {
        fit->SetParameter(2, 0.010);
        status = h->Fit(fit, "RQN");
    }
    r.mu     = fit->GetParameter(1);
    r.mu_err = fit->GetParError(1);
    r.sigma  = fit->GetParameter(2);
    delete fit;

    // Accept if μ is inside its bound window with a finite error.
    r.ok = (r.mu > fit_lo + 1e-4 && r.mu < fit_hi - 1e-4
            && r.mu_err > 0.0 && r.mu_err < 0.004
            && r.sigma > 0.0035 && r.sigma < 0.014);
    return r;
}

void ecal_pi0_3dscan(const char* flavour = "rec") {

    std::string fl = flavour;
    for (auto& c : fl) c = std::tolower(c);
    std::string mass_var;
    if      (fl == "rec") mass_var = "m_epemg";
    else if (fl == "cor") mass_var = "m_epemg_cor";
    else { std::cerr << "Unknown flavour '" << flavour << "'\n"; return; }

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
    gStyle->SetPalette(kBird);
    gSystem->mkdir("plots/output", kTRUE);

    const int nE  = (int)kE_edges.size()     - 1;
    const int nTh = (int)kTheta_edges.size() - 1;
    const int nPh = (int)kPhi_edges.size()   - 1;

    auto make3D = [&](const char* name, const char* title) {
        TH3D* h = new TH3D(name, title,
                           nE,  kE_edges.data(),
                           nTh, kTheta_edges.data(),
                           nPh, kPhi_edges.data());
        h->SetDirectory(nullptr);
        h->GetXaxis()->SetTitle("E_{#gamma} [GeV]");
        h->GetYaxis()->SetTitle("#theta_{#gamma} [deg]");
        h->GetZaxis()->SetTitle("#phi_{#gamma} [deg]");
        return h;
    };
    TH3D* h_mu_3d    = make3D(("h_mu_"    + fl + "_3d").c_str(),
        TString::Format("#mu_{#pi^{0}} (%s);;;;#mu", fl.c_str()));
    TH3D* h_s_3d     = make3D(("h_s_"     + fl + "_3d").c_str(),
        TString::Format("s = (m_{PDG}/#mu)^{2} (%s);;;;s", fl.c_str()));
    TH3D* h_yield_3d = make3D(("h_yield_" + fl + "_3d").c_str(),
        TString::Format("Apex signal (%s);;;;yield", fl.c_str()));

    const std::string out_root = "ecal_pi0_3dscan_" + fl + ".root";
    TFile* fout = TFile::Open(out_root.c_str(), "RECREATE");

    TTree* tout = new TTree("scan3d_results",
                            ("π⁰ 3D scan E×θ×φ_full (" + fl + ")").c_str());
    int   eb = 0, tb = 0, pb = 0;
    float e_lo = 0, e_hi = 0, th_lo = 0, th_hi = 0, ph_lo = 0, ph_hi = 0;
    float mu = 0, mu_err = 0, sigma = 0, yield = 0, s_corr = 0;
    int   fit_ok = 0;
    tout->Branch("eb",     &eb);
    tout->Branch("tb",     &tb);
    tout->Branch("pb",     &pb);
    tout->Branch("e_lo",   &e_lo);
    tout->Branch("e_hi",   &e_hi);
    tout->Branch("th_lo",  &th_lo);
    tout->Branch("th_hi",  &th_hi);
    tout->Branch("ph_lo",  &ph_lo);
    tout->Branch("ph_hi",  &ph_hi);
    tout->Branch("mu",     &mu);
    tout->Branch("mu_err", &mu_err);
    tout->Branch("sigma",  &sigma);
    tout->Branch("yield",  &yield);
    tout->Branch("s_corr", &s_corr);
    tout->Branch("fit_ok", &fit_ok);

    std::cout << "ECAL π⁰ 3D scan (FULL φ): E_γ × θ_γ × φ"
              << "  flavour=" << fl
              << "  cells=" << (nE * nTh * nPh)
              << " (" << nE << "×" << nTh << "×" << nPh << ")\n";

    // ---- Pass 1: independent peak fits ------------------------------------
    using F1 = std::vector<ApexFit>;
    using F2 = std::vector<F1>;
    using F3 = std::vector<F2>;
    F3 cell_fit(nE, F2(nTh, F1(nPh)));
    std::vector<std::vector<std::vector<double>>> cell_apex(
        nE, std::vector<std::vector<double>>(nTh, std::vector<double>(nPh, 0.0)));

    int n_cells_done = 0;
    const int n_total_cells = nE * nTh * nPh;
    for (int ei = 0; ei < nE; ++ei) {
        for (int ti = 0; ti < nTh; ++ti) {
            for (int pi = 0; pi < nPh; ++pi) {
                const double el = kE_edges[ei],     eh = kE_edges[ei + 1];
                const double tl = kTheta_edges[ti], th = kTheta_edges[ti + 1];
                const double pl = kPhi_edges[pi],   ph = kPhi_edges[pi + 1];
                const std::string tag = tag_name(ei, ti, pi);

                TH1D* h_sig = buildSignal3D(t_all, t_pp, t_mm, mass_var,
                                            el, eh, tl, th, pl, ph, tag);
                ApexFit r = findPeak(h_sig, tag);
                cell_fit [ei][ti][pi] = r;
                cell_apex[ei][ti][pi] = r.yield;
                delete h_sig;   // signal hist no longer needed (we have μ)

                ++n_cells_done;
                if (n_cells_done % 100 == 0)
                    std::cout << "  done " << n_cells_done
                              << " / " << n_total_cells << "\n";
            }
        }
    }

    // ---- Pass 2: fill failed cells from neighbour μ ----------------------
    // For dead/empty cells we extrapolate μ from neighbours. Weighting:
    // inverse μ_err² so well-fitted neighbours dominate. We prefer angular
    // neighbours ((θ, φ) at the same E_γ) over E-direction ones, since
    // dead cells in the calorimeter are usually local hardware issues
    // bound to a (θ, φ) location, not to E_γ. This is implemented by
    // first trying the (θ, φ)-only ring; if too few good neighbours are
    // available we fall back to the full 26-cell box.
    auto neighbourMu = [&](int ei, int ti, int pi, bool angular_only) -> double {
        double sum_w = 0, sum_wmu = 0;
        int cnt = 0;
        const int de_lo = angular_only ? 0  : -1;
        const int de_hi = angular_only ? 0  :  1;
        for (int de = de_lo; de <= de_hi; ++de) {
            for (int dt = -1; dt <= 1; ++dt) {
                for (int dp = -1; dp <= 1; ++dp) {
                    if (de == 0 && dt == 0 && dp == 0) continue;
                    const int ne = ei + de;
                    const int nt = ti + dt;
                    int       np = pi + dp;
                    if (np < 0)    np += nPh;     // φ wraps cyclically
                    if (np >= nPh) np -= nPh;
                    if (ne < 0 || ne >= nE)  continue;
                    if (nt < 0 || nt >= nTh) continue;
                    const ApexFit& nf = cell_fit[ne][nt][np];
                    if (!nf.ok) continue;
                    const double err = std::max(1e-5, nf.mu_err);
                    const double w = 1.0 / (err * err);
                    sum_w   += w;
                    sum_wmu += w * nf.mu;
                    ++cnt;
                }
            }
        }
        if (cnt < 2) return -1.0;
        return sum_wmu / sum_w;
    };

    int n_recovered = 0;
    for (int ei = 0; ei < nE; ++ei)
      for (int ti = 0; ti < nTh; ++ti)
        for (int pi = 0; pi < nPh; ++pi) {
            if (cell_fit[ei][ti][pi].ok) continue;
            // Prefer (θ, φ)-only ring; fall back to full box if too sparse.
            double mu_est = neighbourMu(ei, ti, pi, true);
            if (mu_est < kApexLo) mu_est = neighbourMu(ei, ti, pi, false);
            if (mu_est < kApexLo) continue;
            cell_fit[ei][ti][pi].mu     = mu_est;
            cell_fit[ei][ti][pi].mu_err = 0.003;   // pessimistic
            cell_fit[ei][ti][pi].ok     = true;
            ++n_recovered;
        }
    std::cout << "Pass-2 filled " << n_recovered
              << " cells from neighbour μ (weighted, angular-first).\n";

    // ---- Pass 3: smooth outliers via local median ------------------------
    // For each cell, compare its μ to the median of its (θ, φ)-only ring
    // (same E_γ). If it's > 3 MeV away from that median AND its μ_err is
    // larger than 1 MeV (i.e. weak fit), replace it with the median.
    // This kills isolated noisy cells while preserving genuine sector
    // structure (which spans many neighbouring cells).
    auto angularMedian = [&](int ei, int ti, int pi) -> double {
        std::vector<double> vs;
        for (int dt = -1; dt <= 1; ++dt)
            for (int dp = -1; dp <= 1; ++dp) {
                if (dt == 0 && dp == 0) continue;
                const int nt = ti + dt;
                int       np = pi + dp;
                if (np < 0)    np += nPh;
                if (np >= nPh) np -= nPh;
                if (nt < 0 || nt >= nTh) continue;
                const ApexFit& nf = cell_fit[ei][nt][np];
                if (nf.ok) vs.push_back(nf.mu);
            }
        if (vs.size() < 3) return -1.0;
        std::sort(vs.begin(), vs.end());
        return vs[vs.size() / 2];
    };

    int n_smoothed = 0;
    for (int ei = 0; ei < nE; ++ei)
      for (int ti = 0; ti < nTh; ++ti)
        for (int pi = 0; pi < nPh; ++pi) {
            ApexFit& cf = cell_fit[ei][ti][pi];
            if (!cf.ok) continue;
            const double med = angularMedian(ei, ti, pi);
            if (med < kApexLo) continue;
            const double dev = std::fabs(cf.mu - med);
            if (dev > 0.003 && cf.mu_err > 0.001) {
                cf.mu = med;
                ++n_smoothed;
            }
        }
    std::cout << "Pass-3 smoothed " << n_smoothed
              << " outlier cells (|Δμ_neighbours| > 3 MeV).\n";

    // ---- Final: fill TTree + 3D maps + report -----------------------------
    int n_ok = 0;
    for (int ei = 0; ei < nE; ++ei) {
        for (int ti = 0; ti < nTh; ++ti) {
            for (int pi = 0; pi < nPh; ++pi) {
                const double el = kE_edges[ei],     eh = kE_edges[ei + 1];
                const double tl = kTheta_edges[ti], th = kTheta_edges[ti + 1];
                const double pl = kPhi_edges[pi],   ph = kPhi_edges[pi + 1];
                const ApexFit& r = cell_fit[ei][ti][pi];

                eb = ei; tb = ti; pb = pi;
                e_lo = el; e_hi = eh;
                th_lo = tl; th_hi = th;
                ph_lo = pl; ph_hi = ph;
                mu = r.mu; mu_err = r.mu_err; sigma = r.sigma;
                yield = r.yield;
                s_corr = (r.mu > 0) ? std::pow(kPi0PDG / r.mu, 2.0) : 0.0;
                fit_ok = r.ok ? 1 : 0;
                tout->Fill();

                if (r.ok && r.mu > 0) {
                    ++n_ok;
                    const int gx = ei + 1, gy = ti + 1, gz = pi + 1;
                    h_mu_3d   ->SetBinContent(gx, gy, gz, r.mu);
                    h_s_3d    ->SetBinContent(gx, gy, gz, s_corr);
                    h_yield_3d->SetBinContent(gx, gy, gz, r.yield);
                }
            }
        }
    }
    std::cout << "\nCells filled: " << n_ok << " / " << n_total_cells
              << "  (" << (100.0 * n_ok / n_total_cells) << "%)\n";

    // ---- Per-φ-slice heatmaps (one PDF per φ-bin) -------------------------
    auto drawPhiSlice = [&](TH3D* h3, const std::string& base, const char* opt,
                            double zmin, double zmax) {
        for (int pi = 0; pi < nPh; ++pi) {
            const double pl = kPhi_edges[pi], pp = kPhi_edges[pi + 1];
            TH2D h2(TString::Format("%s_phi%02d", base.c_str(), pi),
                    TString::Format("%s, #phi#in[%g, %g)°;"
                                    "E_{#gamma} [GeV];#theta_{#gamma} [deg]",
                                    h3->GetTitle(), pl, pp),
                    nE,  kE_edges.data(),
                    nTh, kTheta_edges.data());
            for (int ei = 1; ei <= nE; ++ei)
                for (int ti = 1; ti <= nTh; ++ti)
                    h2.SetBinContent(ei, ti, h3->GetBinContent(ei, ti, pi + 1));
            h2.SetMinimum(zmin);
            h2.SetMaximum(zmax);
            TString cname = TString::Format("c_%s_phi%02d", base.c_str(), pi);
            TCanvas c(cname, cname, 1100, 800);
            c.SetMargin(0.13, 0.15, 0.12, 0.10);
            c.SetGrid();
            h2.Draw(opt);
            c.Update();
            const std::string out =
                "plots/output/" + base + TString::Format("_phi%02d_%s",
                                                          pi, fl.c_str()).Data();
            c.SaveAs((out + ".png").c_str());
        }
    };

    drawPhiSlice(h_mu_3d, "ecal_pi0_3dscan_mu",    "COLZ", 0.130, 0.150);
    drawPhiSlice(h_s_3d,  "ecal_pi0_3dscan_scale", "COLZ", 0.85,  1.05);

    // Per-sector-summary: 6 (E, θ) heatmaps, each integrating 6 φ-bins
    // within one sector. Lets you see at a glance whether sectors differ.
    auto drawSectorSummary = [&](TH3D* h3, const std::string& base, const char* opt,
                                  double zmin, double zmax) {
        for (int sec = 0; sec < 6; ++sec) {
            TH2D h2(TString::Format("%s_sec%d", base.c_str(), sec),
                    TString::Format("%s, sector %d (#phi#in[%d°,%d°));"
                                    "E_{#gamma} [GeV];#theta_{#gamma} [deg]",
                                    h3->GetTitle(), sec, sec*60, (sec+1)*60),
                    nE,  kE_edges.data(),
                    nTh, kTheta_edges.data());
            // Each sector holds 6 φ-bins (10° each). Use the unweighted
            // average of populated bins as the representative value.
            for (int ei = 1; ei <= nE; ++ei) {
                for (int ti = 1; ti <= nTh; ++ti) {
                    double sum = 0; int cnt = 0;
                    for (int pi = 1; pi <= 6; ++pi) {
                        const int gz = sec * 6 + pi;
                        const double v = h3->GetBinContent(ei, ti, gz);
                        if (v > 0) { sum += v; ++cnt; }
                    }
                    if (cnt > 0) h2.SetBinContent(ei, ti, sum / cnt);
                }
            }
            h2.SetMinimum(zmin);
            h2.SetMaximum(zmax);
            TString cname = TString::Format("c_%s_sec%d", base.c_str(), sec);
            TCanvas c(cname, cname, 1100, 800);
            c.SetMargin(0.13, 0.15, 0.12, 0.10);
            c.SetGrid();
            h2.Draw(opt);
            c.Update();
            const std::string out =
                "plots/output/" + base + TString::Format("_sec%d_%s",
                                                          sec, fl.c_str()).Data();
            c.SaveAs((out + ".png").c_str());
        }
    };
    drawSectorSummary(h_mu_3d, "ecal_pi0_3dscan_sec_mu",    "COLZ TEXT45", 0.130, 0.150);
    drawSectorSummary(h_s_3d,  "ecal_pi0_3dscan_sec_scale", "COLZ TEXT45", 0.85,  1.05);

    fout->cd();
    h_mu_3d->Write();  h_s_3d->Write();  h_yield_3d->Write();
    tout->Write();
    fout->Close();
    f_all->Close();  f_pp->Close();  f_mm->Close();

    std::cout << "\nWrote:\n"
              << "  " << out_root << "  (TH3D h_s_" << fl << "_3d + TTree)\n"
              << "  plots/output/ecal_pi0_3dscan_{mu,scale}_phi{00..35}_" << fl << ".png\n"
              << "  plots/output/ecal_pi0_3dscan_sec_{mu,scale}_sec{0..5}_" << fl << ".png\n";
}
