// epem_masses_rec_expsim.C — Joint EXP + SIM mass spectra (RECONSTRUCTED).
//
// Mirrors the FULL plot list of epem_masses_rec_exp.C, with sim
// (output_epem_sim.root) overlaid as a brown step line on top of a
// light-gray stat-error band, rescaled to the exp signal via per-plot
// findBestScale (with integral clamp).
//
// Two drawing modes:
//   - drawJoint        — for CB-extractable observables (3-curve all/CB/sig + sim)
//   - drawJointSingle  — for pure-ECAL observables (single data curve + sim),
//                        used at mult==2 M(γγ) and mult==4 M(γγγγ).
//
// Y-axis log/linear matches the corresponding plot in epem_masses_rec_exp.C.
//   - also_linear = true   → save both _log and _lin variants
//   - linear_only = true   → save only one linear variant (no _log/_lin suffix)
//
// Reads in CWD:
//   output_epem_exp.root, output_epep_exp.root, output_emem_exp.root  (exp)
//   output_epem_sim.root                                              (sim)
//
// Usage:  root -l -b -q plots/epem_masses_rec_expsim.C
//
// Output: plots/output/joint_m_*_rec_expsim.{pdf,png}

#include "PlotUtils.h"
#include "JointPlotter.h"
#include <TFile.h>
#include <TH1.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TSystem.h>
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

void printIntegralsSingle(const char* label, TH1* data, TH1D* sim) {
    std::cout << "\n=== " << label << " ===\n";
    std::cout << "  exp: integral=" << (data ? data->Integral() : 0) << "\n";
    if (sim) std::cout << "  sim: " << sim->Integral() << "  (auto-scaled to exp)\n";
    else     std::cout << "  sim: (no contribution)\n";
}

// Read the same `histpath` from a ROOT file. Sim histograms in main.cc
// are filled with sim_genweight at fill time — no further per-bin
// weighting is needed here.
TH1D* readHistFromFile(const std::string& fname, const std::string& histpath,
                       const std::string& target_name)
{
    TFile* f = TFile::Open(fname.c_str(), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "readHistFromFile: cannot open " << fname << "\n";
        return nullptr;
    }
    auto* h_raw = (TH1*)f->Get(histpath.c_str());
    if (!h_raw) {
        std::cerr << "readHistFromFile: missing " << histpath << " in " << fname << "\n";
        f->Close();
        return nullptr;
    }
    auto* h = (TH1D*)h_raw->Clone(target_name.c_str());
    h->SetDirectory(nullptr);
    if (!h->GetSumw2N()) h->Sumw2();
    f->Close();
    return h;
}

// Best-chi2 rotation picker (exp-style, unweighted) for one of the three
// exp inputs. Mirrors bestChi2Histogram from epem_masses_rec_exp.C.
TH1D* bestChi2Exp(const std::string& file_path,
                  int nbins, double xmin, double xmax,
                  const std::string& hname,
                  double m_pi0, double sigma_pi0,
                  double m_eta, double sigma_eta,
                  double pi0_lo, double pi0_hi,
                  double eta_lo, double eta_hi)
{
    TFile* f = TFile::Open(file_path.c_str(), "READ");
    if (!f || f->IsZombie()) { std::cerr << "bestChi2Exp: cannot open " << file_path << "\n"; return nullptr; }
    auto t = (TNtuple*)f->Get("epemggg_nt");
    if (!t) { std::cerr << "bestChi2Exp: no epemggg_nt in " << file_path << "\n"; return nullptr; }

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
    std::cout << "  bestChi2Exp(" << file_path << "): kept=" << n_kept << "\n";
    return h;
}

// Best-chi2 picker for sim — applies sim_genweight at fill time.
TH1D* bestChi2Sim(const std::string& file_path,
                  int nbins, double xmin, double xmax,
                  const std::string& hname,
                  double m_pi0, double sigma_pi0,
                  double m_eta, double sigma_eta,
                  double pi0_lo, double pi0_hi,
                  double eta_lo, double eta_hi)
{
    TFile* f = TFile::Open(file_path.c_str(), "READ");
    if (!f || f->IsZombie()) { std::cerr << "bestChi2Sim: cannot open " << file_path << "\n"; return nullptr; }
    auto t = (TNtuple*)f->Get("epemggg_nt");
    if (!t) { std::cerr << "bestChi2Sim: no epemggg_nt in " << file_path << "\n"; return nullptr; }

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
    std::cout << "  bestChi2Sim(" << file_path << "): kept=" << n_kept << "\n";
    return h;
}

}  // namespace

void epem_masses_rec_expsim() {

    PlotUtils exp("output_epem_exp.root",
                  "output_epep_exp.root",
                  "output_emem_exp.root");
    JointPlotter::SimSource sim("output_epem_sim.root");

    // CB-extractable joint plot helper.
    //   logy        — log Y for the primary plot
    //   also_linear — also emit a _lin (linear) sibling next to _log
    //   linear_only — single linear plot, no log variant
    auto plot = [&](const char* nt, const char* var,
                    int nbins, double xmin, double xmax,
                    const char* cut, const char* title,
                    const char* basename, bool logy,
                    bool also_linear = false, bool linear_only = false) {

        TH1D *a, *c, *s;
        std::tie(a, c, s) = exp.drawSignal(nt, var, nbins, xmin, xmax, cut, title);
        if (!a) {
            std::cerr << "epem_masses_rec_expsim: drawSignal failed for " << basename << "\n";
            return;
        }

        TH1D* h_sim = sim.draw(nt, var, nbins, xmin, xmax, cut);
        if (h_sim) JointPlotter::styleSimLine(h_sim);
        double scale = h_sim ? JointPlotter::rescaleSimToData(h_sim, s) : 1.0;

        std::string base = std::string(basename) + "_rec_expsim";

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

    // No-CB joint plot helper (read pre-filled histogram from exp & sim file).
    auto plotDirect = [&](const std::string& histpath,
                          const std::string& title,
                          const std::string& basename, bool logy,
                          double xmin = 0.0, double xmax = 0.0) {

        TH1D* h_data = readHistFromFile("output_epem_exp.root", histpath,
                                        std::string("h_data_") + basename);
        TH1D* h_sim  = readHistFromFile("output_epem_sim.root", histpath,
                                        std::string("h_sim_")  + basename);
        if (!h_data) return;

        // Style data with PlotUtils convention (black markers, error bars).
        h_data->SetMarkerStyle(20);
        h_data->SetMarkerSize(0.6);
        h_data->SetMarkerColor(kBlack);
        h_data->SetLineColor(kBlack);

        double scale = h_sim ? JointPlotter::rescaleSimToData(h_sim, h_data) : 1.0;

        std::string base = basename + "_rec_expsim";
        auto* cv = JointPlotter::drawJointSingle(
            h_data, h_sim, title, ("c_" + base).c_str(), logy, scale,
            xmin, xmax);
        JointPlotter::save(cv, base);
        printIntegralsSingle(base.c_str(), h_data, h_sim);
    };

    std::cout << "\n=== Section 1: dilepton + mult==1 / mult==2 ===\n";

    // -- M(e+e-) deliberately NOT here — covered by mass_spectra_rec_expsim.C
    //    (joint_mass_ee_no_oa_rec_expsim / joint_mass_ee_oa4_rec_expsim).

    // -- M(e+e-#gamma) ECAL N_gamma==1, log + linear --------------------
    plot("epemg_nt", "epemg_mass",
         200, 0.0, 1.0, "",
         "M(e^{+}e^{-}#gamma)  (ECAL N_{#gamma}=1, #pi^{0}/#eta Dalitz);"
         "M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Counts",
         "m_epemg", /*logy=*/true, /*also_linear=*/true);

    // -- M(#gamma#gamma) mult==2 — direct histogram, NO CB, linear -----
    plotDirect("ecal_only/mass_gg",
               "M(#gamma#gamma)  (ECAL N_{#gamma}=2);"
               "M_{#gamma#gamma} [GeV/c^{2}];Counts",
               "m_gg", /*logy=*/false);

    // -- M(epem + gg) mult==2 -------------------------------------------
    plot("epemgg_nt", "m_epemgg",
         280, 0.0, 1.4, "",
         "M(e^{+}e^{-}#gamma#gamma)  (ECAL N_{#gamma}=2);"
         "M_{e^{+}e^{-}#gamma#gamma} [GeV/c^{2}];Counts",
         "m_epemgg", /*logy=*/true);

    // -- M(epem + gg) under narrow pi0 — omega candidate ----------------
    plot("epemgg_nt", "m_epemgg",
         280, 0.0, 1.4, "pi0_pass_narrow==1",
         "#omega(782) #rightarrow #pi^{0}(#gamma#gamma) e^{+}e^{-}  "
         "[M(#gamma#gamma) #in #pi^{0} narrow];"
         "M_{e^{+}e^{-}#gamma#gamma} [GeV/c^{2}];Counts",
         "m_epemgg_pi0narrow", /*logy=*/true);

    std::cout << "\n=== Section 2: mult==3 rotational combinatorics ===\n";

    // -- M(epemg) all rotations ------------------------------------------
    plot("epemggg_nt", "m_epemg",
         280, 0.0, 1.4, "",
         "M(e^{+}e^{-}#gamma)  (ECAL N_{#gamma}=3, all rotations);"
         "M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Counts",
         "m_epemg_mult3", /*logy=*/true);

    // -- M(#gamma#gamma) all rotations — linear only --------------------
    plot("epemggg_nt", "m_gg",
         200, 0.0, 1.0, "",
         "M(#gamma#gamma)  (ECAL N_{#gamma}=3, all rotations);"
         "M_{#gamma#gamma} [GeV/c^{2}];Counts",
         "m_gg_mult3", /*logy=*/false, /*also_linear=*/false, /*linear_only=*/true);

    // -- M(epemg) under M(gg-other) #in pi0 narrow -----------------------
    plot("epemggg_nt", "m_epemg",
         280, 0.0, 1.4, "pi0_pass_narrow==1",
         "M(e^{+}e^{-}#gamma)  (ECAL N_{#gamma}=3, "
         "M(#gamma#gamma)_{other} #in #pi^{0} narrow);"
         "M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Counts",
         "m_epemg_mult3_pi0narrow", /*logy=*/true);

    // -- M(#gamma#gamma) under M(epemg) in eta — linear only ------------
    plot("epemggg_nt", "m_gg",
         200, 0.0, 1.0, "eta_pass==1",
         "M(#gamma#gamma)  (ECAL N_{#gamma}=3, "
         "M(e^{+}e^{-}#gamma) #in #eta window);"
         "M_{#gamma#gamma} [GeV/c^{2}];Counts",
         "m_gg_mult3_eta", /*logy=*/false, /*also_linear=*/false, /*linear_only=*/true);

    // -- M(epemggg) full compound ----------------------------------------
    plot("epemggg_nt", "m_epemggg",
         360, 0.0, 1.8, "",
         "M(e^{+}e^{-}#gamma#gamma#gamma)  (ECAL N_{#gamma}=3);"
         "M_{e^{+}e^{-}#gamma#gamma#gamma} [GeV/c^{2}];Counts",
         "m_epemggg", /*logy=*/true);

    // -- M(epemggg) under both cuts — a0(980) candidate ------------------
    plot("epemggg_nt", "m_epemggg",
         240, 0.4, 1.6, "pi0_pass_narrow==1 && eta_pass==1",
         "a_{0}(980) #rightarrow #eta(e^{+}e^{-}#gamma) + #pi^{0}(#gamma#gamma)  "
         "[M(e^{+}e^{-}#gamma) #in #eta, M(#gamma#gamma) #in #pi^{0} narrow];"
         "M_{e^{+}e^{-}#gamma#gamma#gamma} [GeV/c^{2}];Counts",
         "m_epemggg_both", /*logy=*/true);

    // -- a0 candidate: best-chi2 rotation per event + window cuts --------
    {
        const double m_pi0 = 0.135, m_eta = 0.547;
        const double sigma_pi0 = 0.010, sigma_eta = 0.050;
        const double pi0_lo = 0.125, pi0_hi = 0.145;
        const double eta_lo = 0.500, eta_hi = 0.600;

        const int nbins = 240;
        const double xmin = 0.4, xmax = 1.6;
        const std::string base = "m_epemggg_bestchi2_rec_expsim";
        const char* title =
            "a_{0}(980) #rightarrow #eta(e^{+}e^{-}#gamma) + #pi^{0}(#gamma#gamma)  "
            "[best #chi^{2} rotation, M(#gamma#gamma) #in #pi^{0}, M(e^{+}e^{-}#gamma) #in #eta];"
            "M_{e^{+}e^{-}#gamma#gamma#gamma} [GeV/c^{2}];Counts";

        TH1D* h_all  = bestChi2Exp("output_epem_exp.root", nbins, xmin, xmax,
                                   "h_bestchi2_rec_all",
                                   m_pi0, sigma_pi0, m_eta, sigma_eta,
                                   pi0_lo, pi0_hi, eta_lo, eta_hi);
        TH1D* h_pp   = bestChi2Exp("output_epep_exp.root", nbins, xmin, xmax,
                                   "h_bestchi2_rec_pp",
                                   m_pi0, sigma_pi0, m_eta, sigma_eta,
                                   pi0_lo, pi0_hi, eta_lo, eta_hi);
        TH1D* h_mm   = bestChi2Exp("output_emem_exp.root", nbins, xmin, xmax,
                                   "h_bestchi2_rec_mm",
                                   m_pi0, sigma_pi0, m_eta, sigma_eta,
                                   pi0_lo, pi0_hi, eta_lo, eta_hi);
        TH1D* h_simb = bestChi2Sim("output_epem_sim.root", nbins, xmin, xmax,
                                   "h_bestchi2_rec_sim",
                                   m_pi0, sigma_pi0, m_eta, sigma_eta,
                                   pi0_lo, pi0_hi, eta_lo, eta_hi);

        if (h_all && h_pp && h_mm) {
            TH1D* h_cb  = exp.makeCB(h_pp, h_mm, "h_bestchi2_rec_cb");
            TH1D* h_sig = exp.makeSignal(h_all, h_cb, "h_bestchi2_rec_sig");
            exp.styleAll(h_all);
            exp.styleCB(h_cb);
            exp.styleSignal(h_sig);

            if (h_simb) JointPlotter::styleSimLine(h_simb);
            double scale = h_simb ? JointPlotter::rescaleSimToData(h_simb, h_sig) : 1.0;

            auto* cv = JointPlotter::drawJoint(
                h_all, h_cb, h_sig, h_simb, title,
                "c_m_epemggg_bestchi2_rec_expsim_log", /*logy=*/true, scale);
            JointPlotter::save(cv, base);
            printIntegrals(base.c_str(), h_all, h_cb, h_sig, h_simb);
        }
    }

    std::cout << "\n=== Section 3: mult==4 — gggg under (pi0, pi0) ===\n";

    // -- M(gggg) ECAL N_gamma==4 (no CB, linear) -------------------------
    plotDirect("ecal_only/mass_gggg_pi0pi0",
               "K^{0}_{S}(498) #rightarrow #pi^{0}(#gamma#gamma) + #pi^{0}(#gamma#gamma)  "
               "[both M(#gamma#gamma) #in #pi^{0} narrow];"
               "M_{#gamma#gamma#gamma#gamma} [GeV/c^{2}];Counts",
               "m_gggg_pi0pi0", /*logy=*/false,
               /*xmin=*/0.0, /*xmax=*/1.0);

    std::cout << "\nDone. Joint plots in plots/output/joint_m_*_rec_expsim.{pdf,png}\n";
}
