// epem_masses_cor_exp.C — CORRECTED-kinematics counterpart of epem_masses_rec_exp.C.
//
// Reads from `dilepton_nt_cor`, `epemg_nt_cor`, `epemgg_nt_cor`, `epemggg_nt_cor`
// (mirror ntuples where compound observables are computed from
// KinematicType::CORRECTED). Field names match the RECONSTRUCTED ntuples, so
// cut expressions stay identical; only the underlying values differ.
//
// Photon (gamma) has no measured momentum correction — REC≡COR for ECAL clusters,
// so the pure-gamma plots (M(gg) mult==2 direct, M(gggg) mult==4) are NOT
// duplicated here; see epem_masses_rec_exp.C for those.
//
// Reads three outputs in CWD:
//   output_epem_exp.root, output_epep_exp.root, output_emem_exp.root
//
// Usage:  root -l -b -q plots/epem_masses_cor_exp.C
//
// Output: plots/output/m_*_cor_exp.{pdf,png}

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

void printIntegrals(const char* label, TH1D* a, TH1D* c, TH1D* s) {
    std::cout << "  " << label
              << "  all=" << (a ? a->Integral() : 0)
              << "  CB=" << (c ? c->Integral() : 0)
              << "  sig=" << (s ? s->Integral() : 0) << "\n";
}

void forceYCap(TH1D* h_all, double data_max, TCanvas* cv, bool logy) {
    h_all->SetMaximum(data_max * (logy ? 3.0 : 1.2));
    h_all->SetMinimum(logy ? 0.5 : 0.0);
    cv->Update();
}

void plotTriple(PlotUtils& pu,
                const char* nt_name, const char* varexpr,
                int nbins, double xmin, double xmax,
                const char* cut,
                const char* title_with_axes,
                const char* basename,
                bool also_linear = false,
                bool linear_only = false)
{
    TH1D *a, *c, *s;
    std::tie(a, c, s) = pu.drawSignal(nt_name, varexpr, nbins, xmin, xmax, cut, "");
    if (!a) {
        std::cerr << "plotTriple: drawSignal failed for " << basename << "\n";
        return;
    }

    double data_max = a->GetMaximum();

    if (linear_only) {
        auto* cv = pu.drawTriple(a, c, s, title_with_axes,
                                 (std::string("c_") + basename).c_str(),
                                 /*logy=*/false);
        forceYCap(a, data_max, cv, /*logy=*/false);
        pu.save(cv, basename);
    } else {
        auto* cv_log = pu.drawTriple(a, c, s, title_with_axes,
                                     (std::string("c_") + basename + "_log").c_str(),
                                     /*logy=*/true);
        forceYCap(a, data_max, cv_log, /*logy=*/true);
        pu.save(cv_log, also_linear ? (std::string(basename) + "_log").c_str() : basename);

        if (also_linear) {
            auto* cv_lin = pu.drawTriple(a, c, s, title_with_axes,
                                         (std::string("c_") + basename + "_lin").c_str(),
                                         /*logy=*/false);
            forceYCap(a, data_max, cv_lin, /*logy=*/false);
            pu.save(cv_lin, (std::string(basename) + "_lin").c_str());
        }
    }

    printIntegrals(basename, a, c, s);
}

// Best-chi2 rotation picker on CORRECTED kinematics.
//
// Reads from `epemggg_nt_cor` — same field layout as the REC ntuple, but
// m_epemg / m_epemggg are computed from CORRECTED four-momenta. The picked
// rotation must minimise chi2 against (eta, pi0) AND fall inside both narrow
// mass windows (chi2 alone could pick a rotation kinematically far from the
// resonances).
TH1D* bestChi2HistogramCor(const char* file_path,
                           int nbins, double xmin, double xmax,
                           const char* hname,
                           double m_pi0, double sigma_pi0,
                           double m_eta, double sigma_eta,
                           double pi0_lo, double pi0_hi,
                           double eta_lo, double eta_hi)
{
    TFile* f = TFile::Open(file_path, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "bestChi2HistogramCor: cannot open " << file_path << "\n";
        return nullptr;
    }
    auto t = (TNtuple*)f->Get("epemggg_nt_cor");
    if (!t) {
        std::cerr << "bestChi2HistogramCor: no epemggg_nt_cor in " << file_path << "\n";
        return nullptr;
    }

    TH1D* h = new TH1D(hname, "", nbins, xmin, xmax);
    h->SetDirectory(nullptr);
    h->Sumw2();

    Float_t m_epemg, m_gg, m_epemggg;
    t->SetBranchAddress("m_epemg",   &m_epemg);
    t->SetBranchAddress("m_gg",      &m_gg);
    t->SetBranchAddress("m_epemggg", &m_epemggg);

    Long64_t n = t->GetEntries();
    Long64_t n_events = 0, n_kept = 0;
    for (Long64_t i = 0; i + 2 < n; i += 3) {
        double best_chi2   = 1e30;
        double best_m_full = 0.0, best_m_epemg = 0.0, best_m_gg = 0.0;
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
            }
        }
        ++n_events;
        if (best_m_epemg >= eta_lo && best_m_epemg <= eta_hi &&
            best_m_gg    >= pi0_lo && best_m_gg    <= pi0_hi) {
            h->Fill(best_m_full);
            ++n_kept;
        }
    }
    std::cout << "  bestChi2Cor: " << n_kept << " / " << n_events
              << " events kept from " << file_path << "\n";
    return h;
}

}  // namespace

void epem_masses_cor_exp() {

    PlotUtils pu("output_epem_exp.root",
                 "output_epep_exp.root",
                 "output_emem_exp.root");

    std::cout << "\n=== Section 1: dilepton + mult==1 / mult==2 (cor) ===\n";

    // -- M(e+e-) ---------------------------------------------------------
    plotTriple(pu, "dilepton_nt_cor", "m_ee",
               240, 0.0, 1.2, "",
               "M(e^{+}e^{-}) (cor);M_{e^{+}e^{-}} [GeV/c^{2}];Counts",
               "m_ee_cor_exp");

    // -- M(e+e- gamma)  ECAL N_gamma==1 ----------------------------------
    plotTriple(pu, "epemg_nt_cor", "epemg_mass",
               200, 0.0, 1.0, "",
               "M(e^{+}e^{-}#gamma) (cor)  (ECAL N_{#gamma}=1, #pi^{0}/#eta Dalitz);"
               "M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Counts",
               "m_epemg_cor_exp",
               /*also_linear=*/true);

    // -- M(epem + gg) ECAL N_gamma==2 ------------------------------------
    plotTriple(pu, "epemgg_nt_cor", "m_epemgg",
               280, 0.0, 1.4, "",
               "M(e^{+}e^{-}#gamma#gamma) (cor)  (ECAL N_{#gamma}=2);"
               "M_{e^{+}e^{-}#gamma#gamma} [GeV/c^{2}];Counts",
               "m_epemgg_cor_exp");

    // -- M(epem + gg) under narrow pi0 — omega candidate -----------------
    plotTriple(pu, "epemgg_nt_cor", "m_epemgg",
               280, 0.0, 1.4, "pi0_pass_narrow==1",
               "#omega(782) #rightarrow #pi^{0}(#gamma#gamma) e^{+}e^{-} (cor)  "
               "[M(#gamma#gamma) #in #pi^{0} narrow];"
               "M_{e^{+}e^{-}#gamma#gamma} [GeV/c^{2}];Counts",
               "m_epemgg_pi0narrow_cor_exp");

    std::cout << "\n=== Section 2: mult==3 rotational combinatorics (cor) ===\n";

    // -- M(epemg) all rotations ------------------------------------------
    plotTriple(pu, "epemggg_nt_cor", "m_epemg",
               280, 0.0, 1.4, "",
               "M(e^{+}e^{-}#gamma) (cor)  (ECAL N_{#gamma}=3, all rotations);"
               "M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Counts",
               "m_epemg_mult3_cor_exp");

    // -- M(epemg) under M(gg-other) #in pi0 narrow -----------------------
    plotTriple(pu, "epemggg_nt_cor", "m_epemg",
               280, 0.0, 1.4, "pi0_pass_narrow==1",
               "M(e^{+}e^{-}#gamma) (cor)  (ECAL N_{#gamma}=3, "
               "M(#gamma#gamma)_{other} #in #pi^{0} narrow);"
               "M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Counts",
               "m_epemg_mult3_pi0narrow_cor_exp");

    // -- M(gg) under M(epemg) in eta window — eta_pass uses CORRECTED ----
    plotTriple(pu, "epemggg_nt_cor", "m_gg",
               200, 0.0, 1.0, "eta_pass==1",
               "M(#gamma#gamma)  (ECAL N_{#gamma}=3, "
               "M(e^{+}e^{-}#gamma)_{cor} #in #eta window);"
               "M_{#gamma#gamma} [GeV/c^{2}];Counts",
               "m_gg_mult3_eta_cor_exp",
               /*also_linear=*/false,
               /*linear_only=*/true);

    // -- M(epemggg) full compound ----------------------------------------
    plotTriple(pu, "epemggg_nt_cor", "m_epemggg",
               360, 0.0, 1.8, "",
               "M(e^{+}e^{-}#gamma#gamma#gamma) (cor)  (ECAL N_{#gamma}=3);"
               "M_{e^{+}e^{-}#gamma#gamma#gamma} [GeV/c^{2}];Counts",
               "m_epemggg_cor_exp");

    // -- M(epemggg) under both cuts — a0(980) candidate ------------------
    plotTriple(pu, "epemggg_nt_cor", "m_epemggg",
               240, 0.4, 1.6, "pi0_pass_narrow==1 && eta_pass==1",
               "a_{0}(980) #rightarrow #eta(e^{+}e^{-}#gamma) + #pi^{0}(#gamma#gamma) (cor)  "
               "[M(e^{+}e^{-}#gamma) #in #eta, M(#gamma#gamma) #in #pi^{0} narrow];"
               "M_{e^{+}e^{-}#gamma#gamma#gamma} [GeV/c^{2}];Counts",
               "m_epemggg_both_cor_exp");

    // -- a0 candidate: best-chi2 rotation per event + window cuts (cor) --
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
        const char*  base  = "m_epemggg_bestchi2_cor_exp";
        const char*  title =
            "a_{0}(980) #rightarrow #eta(e^{+}e^{-}#gamma) + #pi^{0}(#gamma#gamma) (cor)  "
            "[best #chi^{2} rotation, M(#gamma#gamma) #in #pi^{0}, M(e^{+}e^{-}#gamma) #in #eta];"
            "M_{e^{+}e^{-}#gamma#gamma#gamma} [GeV/c^{2}];Counts";

        TH1D* h_all = bestChi2HistogramCor("output_epem_exp.root", nbins, xmin, xmax,
                                           "h_bestchi2_cor_all",
                                           m_pi0, sigma_pi0, m_eta, sigma_eta,
                                           pi0_lo, pi0_hi, eta_lo, eta_hi);
        TH1D* h_pp  = bestChi2HistogramCor("output_epep_exp.root", nbins, xmin, xmax,
                                           "h_bestchi2_cor_pp",
                                           m_pi0, sigma_pi0, m_eta, sigma_eta,
                                           pi0_lo, pi0_hi, eta_lo, eta_hi);
        TH1D* h_mm  = bestChi2HistogramCor("output_emem_exp.root", nbins, xmin, xmax,
                                           "h_bestchi2_cor_mm",
                                           m_pi0, sigma_pi0, m_eta, sigma_eta,
                                           pi0_lo, pi0_hi, eta_lo, eta_hi);

        if (h_all && h_pp && h_mm) {
            TH1D* h_cb  = pu.makeCB(h_pp, h_mm, "h_bestchi2_cor_cb");
            TH1D* h_sig = pu.makeSignal(h_all, h_cb, "h_bestchi2_cor_sig");
            pu.styleAll(h_all);
            pu.styleCB(h_cb);
            pu.styleSignal(h_sig);

            double data_max = h_all->GetMaximum();

            auto* cv_log = pu.drawTriple(h_all, h_cb, h_sig, title,
                                         "c_m_epemggg_bestchi2_cor_exp_log", true);
            forceYCap(h_all, data_max, cv_log, true);
            pu.save(cv_log, base);
            printIntegrals(base, h_all, h_cb, h_sig);
        }
    }

    std::cout << "\nDone. All cor plots in plots/output/\n";
}
