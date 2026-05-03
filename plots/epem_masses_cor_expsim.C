// epem_masses_cor_expsim.C — Joint EXP + SIM mass spectra (CORRECTED).
//
// Mirrors the FULL plot list of epem_masses_cor_exp.C, with sim
// (output_epem_sim.root) overlaid as a brown step line on top of a
// light-gray stat-error band, rescaled to the exp signal via per-plot
// findBestScale (with integral clamp).
//
// Reads from `*_nt_cor` (mirror ntuples where compound observables are
// computed from KinematicType::CORRECTED) on both exp and sim sides.
// Photon (gamma) has no measured momentum correction, so pure-gamma
// observables (REC≡COR) are NOT duplicated here; see epem_masses_rec_expsim.C.
//
// Y-axis log/linear matches the corresponding plot in epem_masses_cor_exp.C.
//
// Reads in CWD:
//   output_epem_exp.root, output_epep_exp.root, output_emem_exp.root  (exp)
//   output_epem_sim.root                                              (sim)
//
// Usage:  root -l -b -q plots/epem_masses_cor_expsim.C
//
// Output: plots/output/joint_m_*_cor_expsim.{pdf,png}

#include "PlotUtils.h"
#include "JointPlotter.h"
#include <TFile.h>
#include <TH1.h>
#include <TNtuple.h>
#include <iostream>
#include <string>
#include <tuple>

namespace {

void printIntegrals(const char* label, TH1D* all, TH1D* cb, TH1D* sig, TH1D* sim) {
    std::cout << "\n=== " << label << " ===\n";
    std::cout << "  exp: all=" << (all ? all->Integral() : 0)
              << "   CB="    << (cb  ? cb->Integral()  : 0)
              << "   sig="   << (sig ? sig->Integral() : 0) << "\n";
    if (sim) std::cout << "  sim: " << sim->Integral() << "  (auto-scaled to exp signal)\n";
    else     std::cout << "  sim: (no contribution)\n";
}

// Best-chi2 picker reading from epemggg_nt_cor (CORRECTED), unweighted (exp).
TH1D* bestChi2ExpCor(const std::string& file_path,
                     int nbins, double xmin, double xmax,
                     const std::string& hname,
                     double m_pi0, double sigma_pi0,
                     double m_eta, double sigma_eta,
                     double pi0_lo, double pi0_hi,
                     double eta_lo, double eta_hi)
{
    TFile* f = TFile::Open(file_path.c_str(), "READ");
    if (!f || f->IsZombie()) { std::cerr << "bestChi2ExpCor: cannot open " << file_path << "\n"; return nullptr; }
    auto t = (TNtuple*)f->Get("epemggg_nt_cor");
    if (!t) { std::cerr << "bestChi2ExpCor: no epemggg_nt_cor in " << file_path << "\n"; return nullptr; }

    TH1D* h = new TH1D(hname.c_str(), "", nbins, xmin, xmax);
    h->SetDirectory(nullptr);
    h->Sumw2();

    Float_t m_epemg, m_gg, m_epemggg;
    t->SetBranchAddress("m_epemg",   &m_epemg);
    t->SetBranchAddress("m_gg",      &m_gg);
    t->SetBranchAddress("m_epemggg", &m_epemggg);

    Long64_t n = t->GetEntries();
    Long64_t n_kept = 0;
    for (Long64_t i = 0; i + 2 < n; i += 3) {
        double best_chi2 = 1e30, best_m_full = 0, best_m_epemg = 0, best_m_gg = 0;
        for (int r = 0; r < 3; ++r) {
            t->GetEntry(i + r);
            double dx_eta = (m_epemg - m_eta) / sigma_eta;
            double dx_pi0 = (m_gg    - m_pi0) / sigma_pi0;
            double chi2   = dx_eta*dx_eta + dx_pi0*dx_pi0;
            if (chi2 < best_chi2) {
                best_chi2 = chi2; best_m_full = m_epemggg;
                best_m_epemg = m_epemg; best_m_gg = m_gg;
            }
        }
        if (best_m_epemg >= eta_lo && best_m_epemg <= eta_hi &&
            best_m_gg    >= pi0_lo && best_m_gg    <= pi0_hi) {
            h->Fill(best_m_full); ++n_kept;
        }
    }
    std::cout << "  bestChi2ExpCor(" << file_path << "): kept=" << n_kept << "\n";
    return h;
}

// Best-chi2 picker on epemggg_nt_cor for sim — applies sim_genweight at fill.
TH1D* bestChi2SimCor(const std::string& file_path,
                     int nbins, double xmin, double xmax,
                     const std::string& hname,
                     double m_pi0, double sigma_pi0,
                     double m_eta, double sigma_eta,
                     double pi0_lo, double pi0_hi,
                     double eta_lo, double eta_hi)
{
    TFile* f = TFile::Open(file_path.c_str(), "READ");
    if (!f || f->IsZombie()) { std::cerr << "bestChi2SimCor: cannot open " << file_path << "\n"; return nullptr; }
    auto t = (TNtuple*)f->Get("epemggg_nt_cor");
    if (!t) { std::cerr << "bestChi2SimCor: no epemggg_nt_cor in " << file_path << "\n"; return nullptr; }

    TH1D* h = new TH1D(hname.c_str(), "", nbins, xmin, xmax);
    h->SetDirectory(nullptr);
    h->Sumw2();

    Float_t m_epemg, m_gg, m_epemggg, w;
    t->SetBranchAddress("m_epemg",       &m_epemg);
    t->SetBranchAddress("m_gg",          &m_gg);
    t->SetBranchAddress("m_epemggg",     &m_epemggg);
    t->SetBranchAddress("sim_genweight", &w);

    Long64_t n = t->GetEntries();
    Long64_t n_kept = 0;
    for (Long64_t i = 0; i + 2 < n; i += 3) {
        double best_chi2 = 1e30, best_m_full = 0, best_m_epemg = 0, best_m_gg = 0, best_w = 0;
        for (int r = 0; r < 3; ++r) {
            t->GetEntry(i + r);
            double dx_eta = (m_epemg - m_eta) / sigma_eta;
            double dx_pi0 = (m_gg    - m_pi0) / sigma_pi0;
            double chi2   = dx_eta*dx_eta + dx_pi0*dx_pi0;
            if (chi2 < best_chi2) {
                best_chi2 = chi2; best_m_full = m_epemggg;
                best_m_epemg = m_epemg; best_m_gg = m_gg; best_w = w;
            }
        }
        if (best_m_epemg >= eta_lo && best_m_epemg <= eta_hi &&
            best_m_gg    >= pi0_lo && best_m_gg    <= pi0_hi) {
            h->Fill(best_m_full, best_w); ++n_kept;
        }
    }
    std::cout << "  bestChi2SimCor(" << file_path << "): kept=" << n_kept << "\n";
    return h;
}

}  // namespace

void epem_masses_cor_expsim() {

    PlotUtils exp("output_epem_exp.root",
                  "output_epep_exp.root",
                  "output_emem_exp.root");
    JointPlotter::SimSource sim("output_epem_sim.root");

    auto plot = [&](const char* nt, const char* var,
                    int nbins, double xmin, double xmax,
                    const char* cut, const char* title,
                    const char* basename, bool logy,
                    bool also_linear = false, bool linear_only = false) {

        TH1D *a, *c, *s;
        std::tie(a, c, s) = exp.drawSignal(nt, var, nbins, xmin, xmax, cut, title);
        if (!a) {
            std::cerr << "epem_masses_cor_expsim: drawSignal failed for " << basename << "\n";
            return;
        }

        TH1D* h_sim = sim.draw(nt, var, nbins, xmin, xmax, cut);
        if (h_sim) JointPlotter::styleSimLine(h_sim);
        double scale = h_sim ? JointPlotter::rescaleSimToData(h_sim, s) : 1.0;

        std::string base = std::string(basename) + "_cor_expsim";

        if (linear_only) {
            auto* cv = JointPlotter::drawJoint(
                a, c, s, h_sim, title, ("c_" + base).c_str(), /*logy=*/false, scale);
            JointPlotter::save(cv, base);
        } else {
            auto* cv_log = JointPlotter::drawJoint(
                a, c, s, h_sim, title,
                ("c_" + base + "_log").c_str(), /*logy=*/true, scale);
            JointPlotter::save(cv_log, also_linear ? base + "_log" : base);

            if (also_linear) {
                auto* cv_lin = JointPlotter::drawJoint(
                    a, c, s, h_sim, title,
                    ("c_" + base + "_lin").c_str(), /*logy=*/false, scale);
                JointPlotter::save(cv_lin, base + "_lin");
            }
        }
        printIntegrals(base.c_str(), a, c, s, h_sim);
    };

    std::cout << "\n=== Section 1: dilepton + mult==1 / mult==2 (cor) ===\n";

    // -- M(e+e-) deliberately NOT here — covered by mass_spectra_cor_expsim.C
    //    (joint_mass_ee_no_oa_cor_expsim / joint_mass_ee_oa4_cor_expsim).

    // -- M(e+e-#gamma) ECAL N_gamma==1, log + linear --------------------
    plot("epemg_nt_cor", "epemg_mass",
         200, 0.0, 1.0, "",
         "M(e^{+}e^{-}#gamma) (cor)  (ECAL N_{#gamma}=1, #pi^{0}/#eta Dalitz);"
         "M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Counts",
         "m_epemg", /*logy=*/true, /*also_linear=*/true);

    plot("epemgg_nt_cor", "m_epemgg",
         280, 0.0, 1.4, "",
         "M(e^{+}e^{-}#gamma#gamma) (cor)  (ECAL N_{#gamma}=2);"
         "M_{e^{+}e^{-}#gamma#gamma} [GeV/c^{2}];Counts",
         "m_epemgg", /*logy=*/true);

    plot("epemgg_nt_cor", "m_epemgg",
         280, 0.0, 1.4, "pi0_pass_narrow==1",
         "#omega(782) #rightarrow #pi^{0}(#gamma#gamma) e^{+}e^{-} (cor)  "
         "[M(#gamma#gamma) #in #pi^{0} narrow];"
         "M_{e^{+}e^{-}#gamma#gamma} [GeV/c^{2}];Counts",
         "m_epemgg_pi0narrow", /*logy=*/true);

    std::cout << "\n=== Section 2: mult==3 rotational combinatorics (cor) ===\n";

    plot("epemggg_nt_cor", "m_epemg",
         280, 0.0, 1.4, "",
         "M(e^{+}e^{-}#gamma) (cor)  (ECAL N_{#gamma}=3, all rotations);"
         "M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Counts",
         "m_epemg_mult3", /*logy=*/true);

    plot("epemggg_nt_cor", "m_epemg",
         280, 0.0, 1.4, "pi0_pass_narrow==1",
         "M(e^{+}e^{-}#gamma) (cor)  (ECAL N_{#gamma}=3, "
         "M(#gamma#gamma)_{other} #in #pi^{0} narrow);"
         "M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Counts",
         "m_epemg_mult3_pi0narrow", /*logy=*/true);

    // -- M(#gamma#gamma) under M(epemg)_cor in eta — linear only --------
    plot("epemggg_nt_cor", "m_gg",
         200, 0.0, 1.0, "eta_pass==1",
         "M(#gamma#gamma)  (ECAL N_{#gamma}=3, "
         "M(e^{+}e^{-}#gamma)_{cor} #in #eta window);"
         "M_{#gamma#gamma} [GeV/c^{2}];Counts",
         "m_gg_mult3_eta", /*logy=*/false, /*also_linear=*/false, /*linear_only=*/true);

    plot("epemggg_nt_cor", "m_epemggg",
         360, 0.0, 1.8, "",
         "M(e^{+}e^{-}#gamma#gamma#gamma) (cor)  (ECAL N_{#gamma}=3);"
         "M_{e^{+}e^{-}#gamma#gamma#gamma} [GeV/c^{2}];Counts",
         "m_epemggg", /*logy=*/true);

    plot("epemggg_nt_cor", "m_epemggg",
         240, 0.4, 1.6, "pi0_pass_narrow==1 && eta_pass==1",
         "a_{0}(980) #rightarrow #eta(e^{+}e^{-}#gamma) + #pi^{0}(#gamma#gamma) (cor)  "
         "[M(e^{+}e^{-}#gamma) #in #eta, M(#gamma#gamma) #in #pi^{0} narrow];"
         "M_{e^{+}e^{-}#gamma#gamma#gamma} [GeV/c^{2}];Counts",
         "m_epemggg_both", /*logy=*/true);

    // -- a0 candidate: best-chi2 rotation per event + window cuts (cor) -
    {
        const double m_pi0 = 0.135, m_eta = 0.547;
        const double sigma_pi0 = 0.010, sigma_eta = 0.050;
        const double pi0_lo = 0.125, pi0_hi = 0.145;
        const double eta_lo = 0.500, eta_hi = 0.600;

        const int nbins = 240;
        const double xmin = 0.4, xmax = 1.6;
        const std::string base = "m_epemggg_bestchi2_cor_expsim";
        const char* title =
            "a_{0}(980) #rightarrow #eta(e^{+}e^{-}#gamma) + #pi^{0}(#gamma#gamma) (cor)  "
            "[best #chi^{2} rotation, M(#gamma#gamma) #in #pi^{0}, M(e^{+}e^{-}#gamma) #in #eta];"
            "M_{e^{+}e^{-}#gamma#gamma#gamma} [GeV/c^{2}];Counts";

        TH1D* h_all  = bestChi2ExpCor("output_epem_exp.root", nbins, xmin, xmax,
                                      "h_bestchi2_cor_all",
                                      m_pi0, sigma_pi0, m_eta, sigma_eta,
                                      pi0_lo, pi0_hi, eta_lo, eta_hi);
        TH1D* h_pp   = bestChi2ExpCor("output_epep_exp.root", nbins, xmin, xmax,
                                      "h_bestchi2_cor_pp",
                                      m_pi0, sigma_pi0, m_eta, sigma_eta,
                                      pi0_lo, pi0_hi, eta_lo, eta_hi);
        TH1D* h_mm   = bestChi2ExpCor("output_emem_exp.root", nbins, xmin, xmax,
                                      "h_bestchi2_cor_mm",
                                      m_pi0, sigma_pi0, m_eta, sigma_eta,
                                      pi0_lo, pi0_hi, eta_lo, eta_hi);
        TH1D* h_simb = bestChi2SimCor("output_epem_sim.root", nbins, xmin, xmax,
                                      "h_bestchi2_cor_sim",
                                      m_pi0, sigma_pi0, m_eta, sigma_eta,
                                      pi0_lo, pi0_hi, eta_lo, eta_hi);

        if (h_all && h_pp && h_mm) {
            TH1D* h_cb  = exp.makeCB(h_pp, h_mm, "h_bestchi2_cor_cb");
            TH1D* h_sig = exp.makeSignal(h_all, h_cb, "h_bestchi2_cor_sig");
            exp.styleAll(h_all);
            exp.styleCB(h_cb);
            exp.styleSignal(h_sig);

            if (h_simb) JointPlotter::styleSimLine(h_simb);
            double scale = h_simb ? JointPlotter::rescaleSimToData(h_simb, h_sig) : 1.0;

            auto* cv = JointPlotter::drawJoint(
                h_all, h_cb, h_sig, h_simb, title,
                "c_m_epemggg_bestchi2_cor_expsim_log", /*logy=*/true, scale);
            JointPlotter::save(cv, base);
            printIntegrals(base.c_str(), h_all, h_cb, h_sig, h_simb);
        }
    }

    std::cout << "\nDone. Joint cor plots in plots/output/joint_m_*_cor_expsim.{pdf,png}\n";
}
