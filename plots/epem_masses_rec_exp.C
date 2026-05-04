// epem_masses_rec_exp.C — Mass spectra suite (RECONSTRUCTED, exp data).
//
// Reads from `dilepton_nt`, `epemg_nt`, `epemgg_nt`, `epemggg_nt` (REC).
// CB extraction uses the like-sign exp samples.
//
// Reads three outputs in CWD:
//   output_epem_exp.root, output_epep_exp.root, output_emem_exp.root
//
// Usage:  root -l -b -q plots/epem_masses_rec_exp.C
//
// Output: plots/output/m_*_rec_exp.{pdf,png}

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

// Force Y-axis cap based on the raw data max:
//   log    -> 3.0 * data_max  (log needs lots of headroom — legend + decade margin)
//   linear -> 1.2 * data_max  (tight, almost touches the top bin)
// Overrides drawTriple's defaults and works around a bug where a second
// drawTriple call (linear after log) would multiply the already-modified max.
void forceYCap(TH1D* h_all, double data_max, TCanvas* cv, bool logy) {
    h_all->SetMaximum(data_max * (logy ? 3.0 : 1.2));
    h_all->SetMinimum(logy ? 0.5 : 0.0);
    cv->Update();
}

// drawSignal + drawTriple for log Y, optionally also linear Y.
// title_with_axes must be in ROOT format "annotation;xtitle;ytitle" — the
// semicolons make ROOT split the string into title + axis labels. Without them
// SetTitle wipes the axis labels set by drawSignal.
void plotTriple(PlotUtils& pu,
                const char* nt_name, const char* varexpr,
                int nbins, double xmin, double xmax,
                const char* cut,
                const char* title_with_axes,
                const char* basename,
                bool also_linear = false,
                bool linear_only = false)   // skip log, do single linear plot
{
    TH1D *a, *c, *s;
    std::tie(a, c, s) = pu.drawSignal(nt_name, varexpr, nbins, xmin, xmax, cut, "");
    if (!a) {
        std::cerr << "plotTriple: drawSignal failed for " << basename << "\n";
        return;
    }

    // Cache raw data max BEFORE any drawTriple modifies SetMaximum.
    double data_max = a->GetMaximum();

    if (linear_only) {
        // Single linear plot, no _log/_lin suffix.
        auto* cv = pu.drawTriple(a, c, s, title_with_axes,
                                 (std::string("c_") + basename).c_str(),
                                 /*logy=*/false);
        forceYCap(a, data_max, cv, /*logy=*/false);
        pu.save(cv, basename);
    } else {
        // -- log Y --
        auto* cv_log = pu.drawTriple(a, c, s, title_with_axes,
                                     (std::string("c_") + basename + "_log").c_str(),
                                     /*logy=*/true);
        forceYCap(a, data_max, cv_log, /*logy=*/true);
        pu.save(cv_log, also_linear ? (std::string(basename) + "_log").c_str() : basename);

        // -- linear Y --
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

// === best-chi2 rotation picker for mult==3 a0 candidate =====================
// epemggg_nt has 3 consecutive rows per event (one per gamma rotation). For
// each event:
//   1) compute chi2 ~ ((M_epemg-m_eta)/s_eta)^2 + ((M_gg-m_pi0)/s_pi0)^2 for
//      each rotation,
//   2) pick the rotation with smallest chi2,
//   3) ALSO require that picked rotation falls inside both mass windows
//      (M_epemg in [eta_lo, eta_hi] AND M_gg in [pi0_lo, pi0_hi]).
// If the best rotation fails the windows the event is dropped — chi2 alone
// could pick a rotation that minimises distance to (eta, pi0) but is still
// kinematically far from those resonances.
//
// Done entirely in the macro from existing ntuples — no analysis rerun.
TH1D* bestChi2Histogram(const char* file_path,
                        int nbins, double xmin, double xmax,
                        const char* hname,
                        double m_pi0, double sigma_pi0,
                        double m_eta, double sigma_eta,
                        double pi0_lo, double pi0_hi,
                        double eta_lo, double eta_hi)
{
    TFile* f = TFile::Open(file_path, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "bestChi2Histogram: cannot open " << file_path << "\n";
        return nullptr;
    }
    auto t = (TNtuple*)f->Get("epemggg_nt");
    if (!t) {
        std::cerr << "bestChi2Histogram: no epemggg_nt in " << file_path << "\n";
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
        // Best rotation must also satisfy both mass windows
        if (best_m_epemg >= eta_lo && best_m_epemg <= eta_hi &&
            best_m_gg    >= pi0_lo && best_m_gg    <= pi0_hi) {
            h->Fill(best_m_full);
            ++n_kept;
        }
    }
    std::cout << "  bestChi2: " << n_kept << " / " << n_events
              << " events kept from " << file_path << "\n";
    return h;
}

// Single-file, single-curve plot (no CB extraction). Used when the observable
// is purely from ECAL (gamma-gamma in mult==2, gggg in mult==4) and CB
// subtraction would just remove signal that doesn't depend on lepton sign.
// Optional X range (0,0) means "use histogram's own binning range".
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

    h->SetTitle(title_with_axes);                 // sets main + axis labels at once
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
    std::cout << "  " << basename << ": " << h->GetEntries() << " entries\n";
}

}  // namespace

void epem_masses_rec_exp() {

    PlotUtils pu("output_epem_exp.root",
                 "output_epep_exp.root",
                 "output_emem_exp.root");

    std::cout << "\n=== Section 1: dilepton + mult==1 / mult==2 ===\n";

    // -- M(e+e-) deliberately NOT here — covered by mass_spectra_rec_exp.C
    //    (mass_ee_no_oa_rec_exp_{log,lin} / mass_ee_oa4_rec_exp_{log,lin}).

    // -- M(e+e- gamma)  ECAL N_gamma==1  (log + linear)
    //    Both #pi^{0} (135 MeV) AND #eta (547 MeV) Dalitz decays produce e+e-#gamma.
    plotTriple(pu, "epemg_nt", "epemg_mass",
               200, 0.0, 1.0, "",
               "M(e^{+}e^{-}#gamma)  (ECAL N_{#gamma}=1, #pi^{0}/#eta Dalitz);"
               "M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Counts",
               "m_epemg_rec_exp",
               /*also_linear=*/true);

    // -- 1a) M(gg) mult==2 — direct ECAL, no CB --------------------------
    drawDirectHisto("output_epem_exp.root", "ecal_only/mass_gg",
                    "M(#gamma#gamma)  (ECAL N_{#gamma}=2);"
                    "M_{#gamma#gamma} [GeV/c^{2}];Counts",
                    "m_gg_exp",
                    /*logy=*/false);

    // -- 1d) M(epem + gg) ECAL N_gamma==2 --------------------------------
    plotTriple(pu, "epemgg_nt", "m_epemgg",
               280, 0.0, 1.4, "",
               "M(e^{+}e^{-}#gamma#gamma)  (ECAL N_{#gamma}=2);"
               "M_{e^{+}e^{-}#gamma#gamma} [GeV/c^{2}];Counts",
               "m_epemgg_rec_exp");

    // -- 1e) M(epem + gg) under narrow pi0 — omega candidate -------------
    plotTriple(pu, "epemgg_nt", "m_epemgg",
               280, 0.0, 1.4, "pi0_pass_narrow==1",
               "#omega(782) #rightarrow #pi^{0}(#gamma#gamma) e^{+}e^{-}  "
               "[M(#gamma#gamma) #in #pi^{0} narrow];"
               "M_{e^{+}e^{-}#gamma#gamma} [GeV/c^{2}];Counts",
               "m_epemgg_pi0narrow_rec_exp");

    std::cout << "\n=== Section 2: mult==3 rotational combinatorics ===\n";

    // -- 2a) M(epemg) all rotations --------------------------------------
    plotTriple(pu, "epemggg_nt", "m_epemg",
               280, 0.0, 1.4, "",
               "M(e^{+}e^{-}#gamma)  (ECAL N_{#gamma}=3, all rotations);"
               "M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Counts",
               "m_epemg_mult3_rec_exp");

    // -- 2b) M(gg) all rotations — linear Y for clarity of pi0 peak ------
    plotTriple(pu, "epemggg_nt", "m_gg",
               200, 0.0, 1.0, "",
               "M(#gamma#gamma)  (ECAL N_{#gamma}=3, all rotations);"
               "M_{#gamma#gamma} [GeV/c^{2}];Counts",
               "m_gg_mult3_exp",
               /*also_linear=*/false,
               /*linear_only=*/true);

    // -- 2c) M(epemg) under M(gg-other) #in pi0 narrow -------------------
    plotTriple(pu, "epemggg_nt", "m_epemg",
               280, 0.0, 1.4, "pi0_pass_narrow==1",
               "M(e^{+}e^{-}#gamma)  (ECAL N_{#gamma}=3, "
               "M(#gamma#gamma)_{other} #in #pi^{0} narrow);"
               "M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Counts",
               "m_epemg_mult3_pi0narrow_rec_exp");

    // -- 2d) M(gg) under M(epemg) in eta window — linear Y for clarity ---
    plotTriple(pu, "epemggg_nt", "m_gg",
               200, 0.0, 1.0, "eta_pass==1",
               "M(#gamma#gamma)  (ECAL N_{#gamma}=3, "
               "M(e^{+}e^{-}#gamma) #in #eta window);"
               "M_{#gamma#gamma} [GeV/c^{2}];Counts",
               "m_gg_mult3_eta_rec_exp",
               /*also_linear=*/false,
               /*linear_only=*/true);

    // -- 2e) M(epemggg) full compound ------------------------------------
    plotTriple(pu, "epemggg_nt", "m_epemggg",
               360, 0.0, 1.8, "",
               "M(e^{+}e^{-}#gamma#gamma#gamma)  (ECAL N_{#gamma}=3);"
               "M_{e^{+}e^{-}#gamma#gamma#gamma} [GeV/c^{2}];Counts",
               "m_epemggg_rec_exp");

    // -- 2f) M(epemggg) under both cuts — a0(980) candidate --------------
    plotTriple(pu, "epemggg_nt", "m_epemggg",
               240, 0.4, 1.6, "pi0_pass_narrow==1 && eta_pass==1",
               "a_{0}(980) #rightarrow #eta(e^{+}e^{-}#gamma) + #pi^{0}(#gamma#gamma)  "
               "[M(e^{+}e^{-}#gamma) #in #eta, M(#gamma#gamma) #in #pi^{0} narrow];"
               "M_{e^{+}e^{-}#gamma#gamma#gamma} [GeV/c^{2}];Counts",
               "m_epemggg_both_rec_exp");

    // -- 2g) a0 candidate: best-chi2 rotation per event + window cuts ----
    //   For each ECAL N_gamma==3 event we have 3 ways to assign one gamma to
    //   (epem + gamma) = eta-candidate and the other two to (gamma+gamma) =
    //   pi0-candidate. Compute chi2 for each rotation against the (eta, pi0)
    //   nominal masses and pick the rotation with smallest chi2. THEN require
    //   the picked rotation to also fall inside the eta and pi0 narrow mass
    //   windows — chi2 alone could prefer a rotation that's still kinematically
    //   too far from the resonances. Window cuts use the same ranges as the
    //   defineRangeCut() entries in setup_cuts.h.
    {
        const double m_pi0     = 0.135;   // nominal pi0 mass [GeV/c^2]
        const double m_eta     = 0.547;   // nominal eta mass [GeV/c^2]
        const double sigma_pi0 = 0.010;   // ~ half-width of pi0_mass_window_narrow
        const double sigma_eta = 0.050;   // ~ half-width of eta_mass_window
        const double pi0_lo    = 0.125, pi0_hi = 0.145;   // = pi0_mass_window_narrow
        const double eta_lo    = 0.500, eta_hi = 0.600;   // = eta_mass_window

        const int    nbins = 240;
        const double xmin  = 0.4;
        const double xmax  = 1.6;
        const char*  base  = "m_epemggg_bestchi2_rec_exp";
        const char*  title =
            "a_{0}(980) #rightarrow #eta(e^{+}e^{-}#gamma) + #pi^{0}(#gamma#gamma)  "
            "[best #chi^{2} rotation, M(#gamma#gamma) #in #pi^{0}, M(e^{+}e^{-}#gamma) #in #eta];"
            "M_{e^{+}e^{-}#gamma#gamma#gamma} [GeV/c^{2}];Counts";

        TH1D* h_all = bestChi2Histogram("output_epem_exp.root", nbins, xmin, xmax,
                                        "h_bestchi2_all",
                                        m_pi0, sigma_pi0, m_eta, sigma_eta,
                                        pi0_lo, pi0_hi, eta_lo, eta_hi);
        TH1D* h_pp  = bestChi2Histogram("output_epep_exp.root", nbins, xmin, xmax,
                                        "h_bestchi2_pp",
                                        m_pi0, sigma_pi0, m_eta, sigma_eta,
                                        pi0_lo, pi0_hi, eta_lo, eta_hi);
        TH1D* h_mm  = bestChi2Histogram("output_emem_exp.root", nbins, xmin, xmax,
                                        "h_bestchi2_mm",
                                        m_pi0, sigma_pi0, m_eta, sigma_eta,
                                        pi0_lo, pi0_hi, eta_lo, eta_hi);

        if (h_all && h_pp && h_mm) {
            TH1D* h_cb  = pu.makeCB(h_pp, h_mm, "h_bestchi2_cb");
            TH1D* h_sig = pu.makeSignal(h_all, h_cb, "h_bestchi2_sig");
            pu.styleAll(h_all);
            pu.styleCB(h_cb);
            pu.styleSignal(h_sig);

            double data_max = h_all->GetMaximum();

            auto* cv_log = pu.drawTriple(h_all, h_cb, h_sig, title,
                                         "c_m_epemggg_bestchi2_rec_exp_log", true);
            forceYCap(h_all, data_max, cv_log, true);
            pu.save(cv_log, base);
            printIntegrals(base, h_all, h_cb, h_sig);
        }
    }

    std::cout << "\n=== Section 3: mult==4 — gggg under (pi0, pi0) ===\n";

    // -- 3) M(gggg) ECAL N_gamma==4 with both gg pairs in narrow pi0
    //       both M(#gamma#gamma) cuts applied at fill time in ana.
    drawDirectHisto("output_epem_exp.root", "ecal_only/mass_gggg_pi0pi0",
                    "K^{0}_{S}(498) #rightarrow #pi^{0}(#gamma#gamma) + #pi^{0}(#gamma#gamma)  "
                    "[both M(#gamma#gamma) #in #pi^{0} narrow];"
                    "M_{#gamma#gamma#gamma#gamma} [GeV/c^{2}];Counts",
                    "m_gggg_pi0pi0_exp",
                    /*logy=*/false,
                    /*xmin=*/0.0, /*xmax=*/1.0);

    std::cout << "\nDone. All plots in plots/output/\n";
}
