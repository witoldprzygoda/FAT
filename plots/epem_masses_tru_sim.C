// epem_masses_tru_sim.C — Mass spectra suite for the e+e- + ECAL analysis.
//                          SIMULATED-truth (Geant) kinematics, sim-weighted.
//
// Reads `_sim`-suffixed fields from the *_cor ntuples (dilepton_nt_cor,
// epemg_nt_cor, epemgg_nt_cor, epemggg_nt_cor). Photon kinematics are mirrored
// REC≡COR≡SIM at fill time, but lepton/pion four-momenta carry their Geant
// truth. All fills weighted by per-event sim_genweight.
//
// Reads single output in CWD: output_epem_sim.root
//
// Usage:  root -l -b -q plots/epem_masses_tru_sim.C
//
// Output: plots/output/m_*_tru_sim.{pdf,png}

#include "PlotUtils.h"
#include <TFile.h>
#include <TH1.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TSystem.h>
#include <iostream>
#include <string>

namespace {

std::string wcut(const std::string& filter = "") {
    return filter.empty() ? std::string("sim_genweight")
                          : "(" + filter + ")*sim_genweight";
}

void printIntegrals(const char* label, TH1D* h) {
    std::cout << "  " << label
              << "  integral=" << (h ? h->Integral() : 0)
              << "  entries=" << (h ? h->GetEntries() : 0) << "\n";
}

void forceYCap(TH1D* h, double data_max, TCanvas* cv, bool logy) {
    h->SetMaximum(data_max * (logy ? 3.0 : 1.2));
    h->SetMinimum(logy ? 0.5 : 0.0);
    cv->Update();
}

void plotSingle(PlotUtils& pu,
                const char* nt_name, const char* varexpr,
                int nbins, double xmin, double xmax,
                const char* filter,
                const char* title_with_axes,
                const char* basename,
                bool also_linear = false,
                bool linear_only = false)
{
    TH1D* h = pu.drawNtupleSingle(nt_name, varexpr, nbins, xmin, xmax,
                                  wcut(filter), title_with_axes);
    if (!h) {
        std::cerr << "plotSingle: drawNtupleSingle failed for " << basename << "\n";
        return;
    }
    double data_max = h->GetMaximum();

    if (linear_only) {
        auto* cv = pu.drawSingle(h, title_with_axes,
                                 (std::string("c_") + basename).c_str(), false);
        forceYCap(h, data_max, cv, false);
        pu.save(cv, basename);
    } else {
        auto* cv_log = pu.drawSingle(h, title_with_axes,
                                     (std::string("c_") + basename + "_log").c_str(), true);
        forceYCap(h, data_max, cv_log, true);
        pu.save(cv_log, also_linear ? (std::string(basename) + "_log").c_str() : basename);

        if (also_linear) {
            auto* cv_lin = pu.drawSingle(h, title_with_axes,
                                         (std::string("c_") + basename + "_lin").c_str(), false);
            forceYCap(h, data_max, cv_lin, false);
            pu.save(cv_lin, (std::string(basename) + "_lin").c_str());
        }
    }
    printIntegrals(basename, h);
}

// Best-chi2 picker on SIM-truth M(epemg) and M(epemggg). M(gg) is REC≡COR≡SIM
// so the gg branch is shared with rec/cor versions of this rotation analysis.
TH1D* bestChi2HistogramTruSim(const char* file_path,
                              int nbins, double xmin, double xmax,
                              const char* hname,
                              double m_pi0, double sigma_pi0,
                              double m_eta, double sigma_eta,
                              double pi0_lo, double pi0_hi,
                              double eta_lo, double eta_hi)
{
    TFile* f = TFile::Open(file_path, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "bestChi2HistogramTruSim: cannot open " << file_path << "\n";
        return nullptr;
    }
    auto t = (TNtuple*)f->Get("epemggg_nt_cor");
    if (!t) {
        std::cerr << "bestChi2HistogramTruSim: no epemggg_nt_cor in " << file_path << "\n";
        return nullptr;
    }

    TH1D* h = new TH1D(hname, "", nbins, xmin, xmax);
    h->SetDirectory(nullptr);
    h->Sumw2();

    Float_t m_epemg_sim, m_gg, m_epemggg_sim, w;
    t->SetBranchAddress("m_epemg_sim",   &m_epemg_sim);
    t->SetBranchAddress("m_gg",          &m_gg);
    t->SetBranchAddress("m_epemggg_sim", &m_epemggg_sim);
    t->SetBranchAddress("sim_genweight", &w);

    Long64_t n = t->GetEntries();
    Long64_t n_events = 0, n_kept = 0;
    for (Long64_t i = 0; i + 2 < n; i += 3) {
        double best_chi2 = 1e30, best_m_full = 0, best_m_epemg = 0, best_m_gg = 0, best_w = 0;
        for (int r = 0; r < 3; ++r) {
            t->GetEntry(i + r);
            double dx_eta = (m_epemg_sim - m_eta) / sigma_eta;
            double dx_pi0 = (m_gg        - m_pi0) / sigma_pi0;
            double chi2   = dx_eta*dx_eta + dx_pi0*dx_pi0;
            if (chi2 < best_chi2) {
                best_chi2    = chi2;
                best_m_full  = m_epemggg_sim;
                best_m_epemg = m_epemg_sim;
                best_m_gg    = m_gg;
                best_w       = w;
            }
        }
        ++n_events;
        if (best_m_epemg >= eta_lo && best_m_epemg <= eta_hi &&
            best_m_gg    >= pi0_lo && best_m_gg    <= pi0_hi) {
            h->Fill(best_m_full, best_w);
            ++n_kept;
        }
    }
    std::cout << "  bestChi2TruSim: " << n_kept << " / " << n_events
              << " events kept from " << file_path << "\n";
    return h;
}

}  // namespace

void epem_masses_tru_sim() {

    PlotUtils pu("output_epem_sim.root");

    std::cout << "\n=== Section 1: dilepton + mult==1 / mult==2 (truth) ===\n";

    plotSingle(pu, "dilepton_nt_cor", "m_ee_sim",
               240, 0.0, 1.2, "",
               "M(e^{+}e^{-}) (truth);M_{e^{+}e^{-}} [GeV/c^{2}];a.u.",
               "m_ee_tru_sim");

    plotSingle(pu, "epemg_nt_cor", "epemg_mass_sim",
               200, 0.0, 1.0, "",
               "M(e^{+}e^{-}#gamma) (truth)  (ECAL N_{#gamma}=1);"
               "M_{e^{+}e^{-}#gamma} [GeV/c^{2}];a.u.",
               "m_epemg_tru_sim",
               /*also_linear=*/true);

    plotSingle(pu, "epemgg_nt_cor", "m_epemgg_sim",
               280, 0.0, 1.4, "",
               "M(e^{+}e^{-}#gamma#gamma) (truth)  (ECAL N_{#gamma}=2);"
               "M_{e^{+}e^{-}#gamma#gamma} [GeV/c^{2}];a.u.",
               "m_epemgg_tru_sim");

    plotSingle(pu, "epemgg_nt_cor", "m_epemgg_sim",
               280, 0.0, 1.4, "pi0_pass_narrow==1",
               "#omega(782) #rightarrow #pi^{0}(#gamma#gamma) e^{+}e^{-} (truth)  "
               "[M(#gamma#gamma) #in #pi^{0} narrow];"
               "M_{e^{+}e^{-}#gamma#gamma} [GeV/c^{2}];a.u.",
               "m_epemgg_pi0narrow_tru_sim");

    std::cout << "\n=== Section 2: mult==3 rotational combinatorics (truth) ===\n";

    plotSingle(pu, "epemggg_nt_cor", "m_epemg_sim",
               280, 0.0, 1.4, "",
               "M(e^{+}e^{-}#gamma) (truth)  (ECAL N_{#gamma}=3, all rotations);"
               "M_{e^{+}e^{-}#gamma} [GeV/c^{2}];a.u.",
               "m_epemg_mult3_tru_sim");

    // pi0 narrow cut from M(gg) — REC≡COR≡SIM, same flag works
    plotSingle(pu, "epemggg_nt_cor", "m_epemg_sim",
               280, 0.0, 1.4, "pi0_pass_narrow==1",
               "M(e^{+}e^{-}#gamma) (truth)  (ECAL N_{#gamma}=3, "
               "M(#gamma#gamma)_{other} #in #pi^{0} narrow);"
               "M_{e^{+}e^{-}#gamma} [GeV/c^{2}];a.u.",
               "m_epemg_mult3_pi0narrow_tru_sim");

    plotSingle(pu, "epemggg_nt_cor", "m_epemggg_sim",
               360, 0.0, 1.8, "",
               "M(e^{+}e^{-}#gamma#gamma#gamma) (truth)  (ECAL N_{#gamma}=3);"
               "M_{e^{+}e^{-}#gamma#gamma#gamma} [GeV/c^{2}];a.u.",
               "m_epemggg_tru_sim");

    // a0 candidate using truth M(epemg) for the eta window — eta_pass flag is
    // set on REC; for truth selection we'd need a separate tru-window flag.
    // Use the COR eta_pass as a proxy; the a0 best-chi2 picker below applies
    // the cuts directly on truth values.
    plotSingle(pu, "epemggg_nt_cor", "m_epemggg_sim",
               240, 0.4, 1.6, "pi0_pass_narrow==1 && eta_pass==1",
               "a_{0}(980) candidate (truth)  "
               "[M(#gamma#gamma) #in #pi^{0} narrow, M(e^{+}e^{-}#gamma)_{cor} #in #eta];"
               "M_{e^{+}e^{-}#gamma#gamma#gamma} [GeV/c^{2}];a.u.",
               "m_epemggg_both_tru_sim");

    {
        const double m_pi0     = 0.135;
        const double m_eta     = 0.547;
        const double sigma_pi0 = 0.010;
        const double sigma_eta = 0.050;
        const double pi0_lo    = 0.125, pi0_hi = 0.145;
        const double eta_lo    = 0.500, eta_hi = 0.600;

        const int    nbins = 240;
        const double xmin  = 0.4;
        const double xmax  = 1.6;
        const char*  base  = "m_epemggg_bestchi2_tru_sim";
        const char*  title =
            "a_{0}(980) #rightarrow #eta(e^{+}e^{-}#gamma) + #pi^{0}(#gamma#gamma) (truth)  "
            "[best #chi^{2} rotation on truth, windows applied to truth];"
            "M_{e^{+}e^{-}#gamma#gamma#gamma} [GeV/c^{2}];a.u.";

        TH1D* h = bestChi2HistogramTruSim("output_epem_sim.root", nbins, xmin, xmax,
                                          "h_bestchi2_tru_sim",
                                          m_pi0, sigma_pi0, m_eta, sigma_eta,
                                          pi0_lo, pi0_hi, eta_lo, eta_hi);
        if (h) {
            pu.styleAll(h);
            double data_max = h->GetMaximum();
            auto* cv_log = pu.drawSingle(h, title,
                                         "c_m_epemggg_bestchi2_tru_sim_log", true);
            forceYCap(h, data_max, cv_log, true);
            pu.save(cv_log, base);
            printIntegrals(base, h);
        }
    }

    std::cout << "\nDone. All truth plots in plots/output/\n";
}
