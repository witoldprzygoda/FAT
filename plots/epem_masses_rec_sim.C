// epem_masses_rec_sim.C — Mass spectra suite for the e+e- + ECAL analysis.
//                          RECONSTRUCTED kinematics, simulation (single-file).
//
// Reads from `dilepton_nt`, `epemg_nt`, `epemgg_nt`, `epemggg_nt` (REC).
// All fills weighted by the per-event sim_genweight stored on the ntuples.
//
// Reads single output in CWD:
//   output_epem_sim.root
//
// Usage:  root -l -b -q plots/epem_masses_rec_sim.C
//
// Output: plots/output/m_*_rec_sim.{pdf,png}

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

// Helper: compose a weighted TTree::Draw cut. Per-event sim_genweight from
// the ntuple is multiplied by an optional boolean filter.
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

// Single-curve plot helper for sim — wraps drawNtupleSingle + drawSingle with
// optional dual-Y (log + linear) output and Y-cap forcing.
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
                                 (std::string("c_") + basename).c_str(),
                                 /*logy=*/false);
        forceYCap(h, data_max, cv, /*logy=*/false);
        pu.save(cv, basename);
    } else {
        auto* cv_log = pu.drawSingle(h, title_with_axes,
                                     (std::string("c_") + basename + "_log").c_str(),
                                     /*logy=*/true);
        forceYCap(h, data_max, cv_log, /*logy=*/true);
        pu.save(cv_log, also_linear ? (std::string(basename) + "_log").c_str() : basename);

        if (also_linear) {
            auto* cv_lin = pu.drawSingle(h, title_with_axes,
                                         (std::string("c_") + basename + "_lin").c_str(),
                                         /*logy=*/false);
            forceYCap(h, data_max, cv_lin, /*logy=*/false);
            pu.save(cv_lin, (std::string(basename) + "_lin").c_str());
        }
    }

    printIntegrals(basename, h);
}

// Best-chi2 rotation picker for sim — same logic as exp version, but reads from
// epemggg_nt (REC) and applies sim_genweight per-event during histogram fill.
TH1D* bestChi2HistogramSim(const char* file_path,
                           int nbins, double xmin, double xmax,
                           const char* hname,
                           double m_pi0, double sigma_pi0,
                           double m_eta, double sigma_eta,
                           double pi0_lo, double pi0_hi,
                           double eta_lo, double eta_hi)
{
    TFile* f = TFile::Open(file_path, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "bestChi2HistogramSim: cannot open " << file_path << "\n";
        return nullptr;
    }
    auto t = (TNtuple*)f->Get("epemggg_nt");
    if (!t) {
        std::cerr << "bestChi2HistogramSim: no epemggg_nt in " << file_path << "\n";
        return nullptr;
    }

    TH1D* h = new TH1D(hname, "", nbins, xmin, xmax);
    h->SetDirectory(nullptr);
    h->Sumw2();

    Float_t m_epemg, m_gg, m_epemggg, w;
    t->SetBranchAddress("m_epemg",       &m_epemg);
    t->SetBranchAddress("m_gg",          &m_gg);
    t->SetBranchAddress("m_epemggg",     &m_epemggg);
    t->SetBranchAddress("sim_genweight", &w);

    Long64_t n = t->GetEntries();
    Long64_t n_events = 0, n_kept = 0;
    for (Long64_t i = 0; i + 2 < n; i += 3) {
        double best_chi2   = 1e30;
        double best_m_full = 0.0, best_m_epemg = 0.0, best_m_gg = 0.0;
        double best_w = 0.0;
        for (int r = 0; r < 3; ++r) {
            t->GetEntry(i + r);
            double dx_eta = (m_epemg - m_eta) / sigma_eta;
            double dx_pi0 = (m_gg    - m_pi0) / sigma_pi0;
            double chi2   = dx_eta*dx_eta + dx_pi0*dx_pi0;
            if (chi2 < best_chi2) {
                best_chi2    = chi2;
                best_m_full  = m_epemggg;
                best_m_epemg = m_epemg;
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
    std::cout << "  bestChi2Sim: " << n_kept << " / " << n_events
              << " events kept from " << file_path << "\n";
    return h;
}

// Direct ECAL histogram (already weighted at fill time in main.cc), no CB.
void drawDirectHisto(const char* file_path, const char* hist_path,
                     const char* title_with_axes,
                     const char* basename,
                     bool logy = true,
                     double xmin = 0.0, double xmax = 0.0)
{
    TFile* f = TFile::Open(file_path, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "drawDirectHisto: cannot open " << file_path << "\n";
        return;
    }
    auto h = (TH1*)f->Get(hist_path);
    if (!h) {
        std::cerr << "drawDirectHisto: missing " << hist_path << "\n";
        return;
    }

    auto* c = new TCanvas((std::string("c_") + basename).c_str(), basename, 800, 600);
    c->SetMargin(0.12, 0.05, 0.12, 0.08);
    if (logy) c->SetLogy();

    h->SetTitle(title_with_axes);
    h->SetMarkerStyle(20);
    h->SetMarkerSize(0.6);
    h->SetMarkerColor(kBlack);
    h->SetLineColor(kBlack);
    h->GetYaxis()->SetMaxDigits(3);
    if (xmin < xmax) h->GetXaxis()->SetRangeUser(xmin, xmax);

    double data_max = h->GetMaximum();
    h->SetMaximum(data_max * (logy ? 3.0 : 1.2));
    h->SetMinimum(logy ? 0.5 : 0.0);

    h->Draw("E");
    c->Update();

    gSystem->mkdir("plots/output", kTRUE);
    c->SaveAs((std::string("plots/output/") + basename + ".pdf").c_str());
    c->SaveAs((std::string("plots/output/") + basename + ".png").c_str());
    std::cout << "  " << basename << ": integral=" << h->Integral()
              << "  entries=" << h->GetEntries() << "\n";
}

}  // namespace

void epem_masses_rec_sim() {

    PlotUtils pu("output_epem_sim.root");

    std::cout << "\n=== Section 1: dilepton + mult==1 / mult==2 ===\n";

    // -- M(e+e-) ---------------------------------------------------------
    plotSingle(pu, "dilepton_nt", "m_ee",
               240, 0.0, 1.2, "",
               "M(e^{+}e^{-});M_{e^{+}e^{-}} [GeV/c^{2}];a.u.",
               "m_ee_rec_sim");

    // -- M(e+e- gamma)  ECAL N_gamma==1  (log + linear) ------------------
    plotSingle(pu, "epemg_nt", "epemg_mass",
               200, 0.0, 1.0, "",
               "M(e^{+}e^{-}#gamma)  (ECAL N_{#gamma}=1, #pi^{0}/#eta Dalitz);"
               "M_{e^{+}e^{-}#gamma} [GeV/c^{2}];a.u.",
               "m_epemg_rec_sim",
               /*also_linear=*/true);

    // -- M(gg) mult==2 — direct ECAL histogram, no CB --------------------
    drawDirectHisto("output_epem_sim.root", "ecal_only/mass_gg",
                    "M(#gamma#gamma)  (ECAL N_{#gamma}=2);"
                    "M_{#gamma#gamma} [GeV/c^{2}];a.u.",
                    "m_gg_sim",
                    /*logy=*/false);

    // -- M(epem + gg) ECAL N_gamma==2 ------------------------------------
    plotSingle(pu, "epemgg_nt", "m_epemgg",
               280, 0.0, 1.4, "",
               "M(e^{+}e^{-}#gamma#gamma)  (ECAL N_{#gamma}=2);"
               "M_{e^{+}e^{-}#gamma#gamma} [GeV/c^{2}];a.u.",
               "m_epemgg_rec_sim");

    // -- M(epem + gg) under narrow pi0 — omega candidate -----------------
    plotSingle(pu, "epemgg_nt", "m_epemgg",
               280, 0.0, 1.4, "pi0_pass_narrow==1",
               "#omega(782) #rightarrow #pi^{0}(#gamma#gamma) e^{+}e^{-}  "
               "[M(#gamma#gamma) #in #pi^{0} narrow];"
               "M_{e^{+}e^{-}#gamma#gamma} [GeV/c^{2}];a.u.",
               "m_epemgg_pi0narrow_rec_sim");

    std::cout << "\n=== Section 2: mult==3 rotational combinatorics ===\n";

    plotSingle(pu, "epemggg_nt", "m_epemg",
               280, 0.0, 1.4, "",
               "M(e^{+}e^{-}#gamma)  (ECAL N_{#gamma}=3, all rotations);"
               "M_{e^{+}e^{-}#gamma} [GeV/c^{2}];a.u.",
               "m_epemg_mult3_rec_sim");

    plotSingle(pu, "epemggg_nt", "m_gg",
               200, 0.0, 1.0, "",
               "M(#gamma#gamma)  (ECAL N_{#gamma}=3, all rotations);"
               "M_{#gamma#gamma} [GeV/c^{2}];a.u.",
               "m_gg_mult3_sim",
               /*also_linear=*/false,
               /*linear_only=*/true);

    plotSingle(pu, "epemggg_nt", "m_epemg",
               280, 0.0, 1.4, "pi0_pass_narrow==1",
               "M(e^{+}e^{-}#gamma)  (ECAL N_{#gamma}=3, "
               "M(#gamma#gamma)_{other} #in #pi^{0} narrow);"
               "M_{e^{+}e^{-}#gamma} [GeV/c^{2}];a.u.",
               "m_epemg_mult3_pi0narrow_rec_sim");

    plotSingle(pu, "epemggg_nt", "m_gg",
               200, 0.0, 1.0, "eta_pass==1",
               "M(#gamma#gamma)  (ECAL N_{#gamma}=3, "
               "M(e^{+}e^{-}#gamma) #in #eta window);"
               "M_{#gamma#gamma} [GeV/c^{2}];a.u.",
               "m_gg_mult3_eta_rec_sim",
               /*also_linear=*/false,
               /*linear_only=*/true);

    plotSingle(pu, "epemggg_nt", "m_epemggg",
               360, 0.0, 1.8, "",
               "M(e^{+}e^{-}#gamma#gamma#gamma)  (ECAL N_{#gamma}=3);"
               "M_{e^{+}e^{-}#gamma#gamma#gamma} [GeV/c^{2}];a.u.",
               "m_epemggg_rec_sim");

    plotSingle(pu, "epemggg_nt", "m_epemggg",
               240, 0.4, 1.6, "pi0_pass_narrow==1 && eta_pass==1",
               "a_{0}(980) #rightarrow #eta(e^{+}e^{-}#gamma) + #pi^{0}(#gamma#gamma)  "
               "[M(e^{+}e^{-}#gamma) #in #eta, M(#gamma#gamma) #in #pi^{0} narrow];"
               "M_{e^{+}e^{-}#gamma#gamma#gamma} [GeV/c^{2}];a.u.",
               "m_epemggg_both_rec_sim");

    // -- a0 candidate: best-chi2 rotation per event + window cuts --------
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
        const char*  base  = "m_epemggg_bestchi2_rec_sim";
        const char*  title =
            "a_{0}(980) #rightarrow #eta(e^{+}e^{-}#gamma) + #pi^{0}(#gamma#gamma)  "
            "[best #chi^{2} rotation, M(#gamma#gamma) #in #pi^{0}, M(e^{+}e^{-}#gamma) #in #eta];"
            "M_{e^{+}e^{-}#gamma#gamma#gamma} [GeV/c^{2}];a.u.";

        TH1D* h = bestChi2HistogramSim("output_epem_sim.root", nbins, xmin, xmax,
                                       "h_bestchi2_rec_sim",
                                       m_pi0, sigma_pi0, m_eta, sigma_eta,
                                       pi0_lo, pi0_hi, eta_lo, eta_hi);
        if (h) {
            pu.styleAll(h);
            double data_max = h->GetMaximum();
            auto* cv_log = pu.drawSingle(h, title,
                                         "c_m_epemggg_bestchi2_rec_sim_log", true);
            forceYCap(h, data_max, cv_log, true);
            pu.save(cv_log, base);
            printIntegrals(base, h);
        }
    }

    std::cout << "\n=== Section 3: mult==4 — gggg under (pi0, pi0) ===\n";

    drawDirectHisto("output_epem_sim.root", "ecal_only/mass_gggg_pi0pi0",
                    "K^{0}_{S}(498) #rightarrow #pi^{0}(#gamma#gamma) + #pi^{0}(#gamma#gamma)  "
                    "[both M(#gamma#gamma) #in #pi^{0} narrow];"
                    "M_{#gamma#gamma#gamma#gamma} [GeV/c^{2}];a.u.",
                    "m_gggg_pi0pi0_sim",
                    /*logy=*/false,
                    /*xmin=*/0.0, /*xmax=*/1.0);

    std::cout << "\nDone. All plots in plots/output/\n";
}
