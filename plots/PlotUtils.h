/**
 * @file PlotUtils.h
 * @brief Signal extraction and plotting utilities for dilepton analysis
 *
 * Handles: file loading, combinatorial background (CB = 2*sqrt(N++*N--)),
 * signal extraction (sig = all - CB), and consistent styling.
 *
 * Usage in ROOT macros:
 *   #include "PlotUtils.h"
 *   PlotUtils pu("output_epem.root", "output_epep.root", "output_emem.root");
 *   auto [all, cb, sig] = pu.drawSignal("dilepton_nt", "m_ee", 200, 0, 1.0);
 *
 * @author Witold Przygoda (witold.przygoda@uj.edu.pl)
 * @date 2025
 */

#ifndef PLOTUTILS_H
#define PLOTUTILS_H

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TMath.h>
#include <iostream>
#include <tuple>
#include <string>

class PlotUtils {
public:
    /// 3-file constructor: signal + like-sign pairs for CB extraction (data analysis).
    PlotUtils(const std::string& f_all,
              const std::string& f_pp,
              const std::string& f_mm)
    {
        f_all_ = TFile::Open(f_all.c_str(), "READ");
        f_pp_  = TFile::Open(f_pp.c_str(), "READ");
        f_mm_  = TFile::Open(f_mm.c_str(), "READ");

        if (!f_all_ || f_all_->IsZombie()) {
            std::cerr << "PlotUtils: Cannot open " << f_all << "\n";
        }
        if (!f_pp_ || f_pp_->IsZombie()) {
            std::cerr << "PlotUtils: Cannot open " << f_pp << "\n";
        }
        if (!f_mm_ || f_mm_->IsZombie()) {
            std::cerr << "PlotUtils: Cannot open " << f_mm << "\n";
        }

        applyStyle();
        counter_ = 0;
    }

    /// Single-file constructor: one input only (simulation, no CB extraction).
    /// Use drawNtupleSingle / getHistSingle / drawSingle for plotting.
    explicit PlotUtils(const std::string& f_all)
    {
        f_all_ = TFile::Open(f_all.c_str(), "READ");
        f_pp_  = nullptr;
        f_mm_  = nullptr;

        if (!f_all_ || f_all_->IsZombie()) {
            std::cerr << "PlotUtils: Cannot open " << f_all << "\n";
        }

        applyStyle();
        counter_ = 0;
    }

    ~PlotUtils() {
        if (f_all_) f_all_->Close();
        if (f_pp_)  f_pp_->Close();
        if (f_mm_)  f_mm_->Close();
    }

    // ========================================================================
    // Signal from existing histograms
    // ========================================================================

    /**
     * @brief Extract signal from an existing histogram in all 3 files
     * @param histpath Path inside ROOT file (e.g. "dilepton/mass_ee")
     * @return {all, cb, signal} histograms
     */
    std::tuple<TH1D*, TH1D*, TH1D*> getSignal(const std::string& histpath) {
        TH1* h_all_raw = getHist(f_all_, histpath);
        TH1* h_pp_raw  = getHist(f_pp_, histpath);
        TH1* h_mm_raw  = getHist(f_mm_, histpath);

        if (!h_all_raw || !h_pp_raw || !h_mm_raw) {
            std::cerr << "PlotUtils::getSignal: histogram '" << histpath
                      << "' not found in all files\n";
            return {nullptr, nullptr, nullptr};
        }

        std::string base = uniqueName(histpath);

        TH1D* h_all = toTH1D(h_all_raw, base + "_all");
        TH1D* h_pp  = toTH1D(h_pp_raw,  base + "_pp");
        TH1D* h_mm  = toTH1D(h_mm_raw,  base + "_mm");

        TH1D* h_cb  = makeCB(h_pp, h_mm, base + "_cb");
        TH1D* h_sig = makeSignal(h_all, h_cb, base + "_sig");

        styleAll(h_all);
        styleCB(h_cb);
        styleSignal(h_sig);

        return {h_all, h_cb, h_sig};
    }

    // ========================================================================
    // Signal from ntuples
    // ========================================================================

    /**
     * @brief Draw variable from ntuple in all 3 files, extract signal
     * @param ntuple_name Ntuple name (e.g. "dilepton_nt")
     * @param varexpr Variable expression (e.g. "m_ee", "mm_epem_mass2")
     * @param nbins Number of bins
     * @param xmin X-axis minimum
     * @param xmax X-axis maximum
     * @param cut Optional cut string (e.g. "oa_pass==1")
     * @param title Histogram title (e.g. ";M_{ee} [GeV/c^{2}];Counts")
     * @return {all, cb, signal} histograms
     */
    std::tuple<TH1D*, TH1D*, TH1D*> drawSignal(
        const std::string& ntuple_name,
        const std::string& varexpr,
        int nbins, double xmin, double xmax,
        const std::string& cut = "",
        const std::string& title = "")
    {
        std::string base = uniqueName(varexpr);

        TH1D* h_all = drawFromNtuple(f_all_, ntuple_name, varexpr,
                                      nbins, xmin, xmax, cut, base + "_all");
        TH1D* h_pp  = drawFromNtuple(f_pp_, ntuple_name, varexpr,
                                      nbins, xmin, xmax, cut, base + "_pp");
        TH1D* h_mm  = drawFromNtuple(f_mm_, ntuple_name, varexpr,
                                      nbins, xmin, xmax, cut, base + "_mm");

        if (!h_all || !h_pp || !h_mm) {
            std::cerr << "PlotUtils::drawSignal: ntuple drawing failed for '"
                      << varexpr << "'\n";
            return {nullptr, nullptr, nullptr};
        }

        if (!title.empty()) {
            h_all->SetTitle(title.c_str());
        }

        TH1D* h_cb  = makeCB(h_pp, h_mm, base + "_cb");
        TH1D* h_sig = makeSignal(h_all, h_cb, base + "_sig");

        styleAll(h_all);
        styleCB(h_cb);
        styleSignal(h_sig);

        return {h_all, h_cb, h_sig};
    }

    // ========================================================================
    // CB and Signal calculation
    // ========================================================================

    /**
     * @brief Combinatorial background: CB_i = 2 * sqrt(N++_i * N--_i)
     * Error: sigma_CB_i = sqrt(N++_i + N--_i)
     */
    TH1D* makeCB(TH1* h_pp, TH1* h_mm, const std::string& name) {
        TH1D* h_cb = (TH1D*)h_pp->Clone(name.c_str());
        h_cb->Reset();

        int nbins = h_pp->GetNbinsX();
        for (int i = 0; i <= nbins + 1; ++i) {  // include under/overflow
            double npp = h_pp->GetBinContent(i);
            double nmm = h_mm->GetBinContent(i);

            double cb = 2.0 * TMath::Sqrt(npp * nmm);
            double err = TMath::Sqrt(npp + nmm);

            h_cb->SetBinContent(i, cb);
            h_cb->SetBinError(i, err);
        }
        return h_cb;
    }

    /**
     * @brief Signal: sig_i = all_i - CB_i
     * Error: sigma_sig_i = sqrt(N_all_i + N++_i + N--_i)
     * (uses the raw counts, not the CB value, for proper error propagation)
     */
    TH1D* makeSignal(TH1* h_all, TH1* h_cb, const std::string& name) {
        TH1D* h_sig = (TH1D*)h_all->Clone(name.c_str());

        int nbins = h_all->GetNbinsX();
        for (int i = 0; i <= nbins + 1; ++i) {
            double n_all = h_all->GetBinContent(i);
            double n_cb  = h_cb->GetBinContent(i);
            double err_all = h_all->GetBinError(i);
            double err_cb  = h_cb->GetBinError(i);

            h_sig->SetBinContent(i, n_all - n_cb);
            h_sig->SetBinError(i, TMath::Sqrt(err_all * err_all + err_cb * err_cb));
        }
        return h_sig;
    }

    // ========================================================================
    // Drawing
    // ========================================================================

    /**
     * @brief Draw all, CB, and signal on a canvas with legend
     * @param canvasTitle Title shown on top of pad
     * @param canvasName Canvas object name (used for saving)
     * @param logy       true for log-Y, false (default) for linear-Y.
     *                   Y-range defaults are tuned to the chosen scale.
     */
    TCanvas* drawTriple(TH1* h_all, TH1* h_cb, TH1* h_sig,
                        const std::string& canvasTitle = "",
                        const std::string& canvasName = "c1",
                        bool logy = false)
    {
        TCanvas* c = new TCanvas(canvasName.c_str(), canvasTitle.c_str(), 800, 600);
        c->SetMargin(0.12, 0.05, 0.12, 0.08);
        if (logy) c->SetLogy();

        if (!canvasTitle.empty()) {
            h_all->SetTitle(canvasTitle.c_str());
        }

        if (logy) {
            h_all->SetMaximum(h_all->GetMaximum() * 3.0);
            h_all->SetMinimum(0.5);
        } else {
            h_all->SetMaximum(h_all->GetMaximum() * 1.2);
            h_all->SetMinimum(0.0);
        }

        h_all->Draw("E");
        h_cb->Draw("E SAME");
        h_sig->Draw("E SAME");

        TLegend* leg = new TLegend(0.65, 0.72, 0.92, 0.90);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.035);
        leg->AddEntry(h_all, "e^{+}e^{-} (all)", "lpe");
        leg->AddEntry(h_cb, "CB (2#sqrt{N_{++}N_{--}})", "lpe");
        leg->AddEntry(h_sig, "Signal (all - CB)", "lpe");
        leg->Draw();

        c->Update();
        return c;
    }

    // ========================================================================
    // Single-file helpers (no CB extraction)
    // ========================================================================
    // For simulation or any analysis where like-sign samples don't exist.
    // Pair drawNtupleSingle / getHistSingle (data source) with drawSingle
    // (canvas) the same way drawSignal pairs with drawTriple.

    /// Draw a variable from a single ntuple in f_all_, return styled TH1D.
    TH1D* drawNtupleSingle(const std::string& ntuple_name,
                           const std::string& varexpr,
                           int nbins, double xmin, double xmax,
                           const std::string& cut = "",
                           const std::string& title = "")
    {
        TH1D* h = drawFromNtuple(f_all_, ntuple_name, varexpr,
                                 nbins, xmin, xmax, cut,
                                 uniqueName(varexpr));
        if (!h) {
            std::cerr << "PlotUtils::drawNtupleSingle: failed for '" << varexpr << "'\n";
            return nullptr;
        }
        if (!title.empty()) h->SetTitle(title.c_str());
        styleAll(h);
        return h;
    }

    /// Get an existing histogram from f_all_, return styled TH1D copy.
    TH1D* getHistSingle(const std::string& histpath) {
        TH1* raw = getHist(f_all_, histpath);
        if (!raw) {
            std::cerr << "PlotUtils::getHistSingle: histogram '" << histpath
                      << "' not found in f_all_\n";
            return nullptr;
        }
        TH1D* h = toTH1D(raw, uniqueName(histpath));
        styleAll(h);
        return h;
    }

    /// Draw a single histogram on a canvas with consistent styling + log/lin Y.
    TCanvas* drawSingle(TH1* h,
                        const std::string& canvasTitle = "",
                        const std::string& canvasName  = "c1",
                        bool logy = false,
                        const std::string& legendLabel = "e^{+}e^{-} (sim)")
    {
        TCanvas* c = new TCanvas(canvasName.c_str(), canvasTitle.c_str(), 800, 600);
        c->SetMargin(0.12, 0.05, 0.12, 0.08);
        if (logy) c->SetLogy();

        if (!canvasTitle.empty()) {
            h->SetTitle(canvasTitle.c_str());
        }

        if (logy) {
            h->SetMaximum(h->GetMaximum() * 3.0);
            h->SetMinimum(0.5);
        } else {
            h->SetMaximum(h->GetMaximum() * 1.2);
            h->SetMinimum(0.0);
        }

        h->Draw("E");

        TLegend* leg = new TLegend(0.65, 0.82, 0.92, 0.90);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.035);
        leg->AddEntry(h, legendLabel.c_str(), "lpe");
        leg->Draw();

        c->Update();
        return c;
    }

    /**
     * @brief Save canvas to plots/output/{basename}.pdf and .png
     */
    void save(TCanvas* c, const std::string& basename) {
        gSystem->mkdir("plots/output", kTRUE);
        std::string pdf = "plots/output/" + basename + ".pdf";
        std::string png = "plots/output/" + basename + ".png";
        c->SaveAs(pdf.c_str());
        c->SaveAs(png.c_str());
        std::cout << "Saved: " << pdf << "\n";
    }

    // ========================================================================
    // Styling
    // ========================================================================

    void styleAll(TH1* h) {
        h->SetMarkerStyle(20);
        h->SetMarkerSize(0.7);
        h->SetMarkerColor(kBlack);
        h->SetLineColor(kBlack);
        h->SetLineWidth(1);
    }

    void styleCB(TH1* h) {
        h->SetMarkerStyle(20);
        h->SetMarkerSize(0.7);
        h->SetMarkerColor(kRed);
        h->SetLineColor(kRed);
        h->SetLineWidth(1);
    }

    void styleSignal(TH1* h) {
        h->SetMarkerStyle(20);
        h->SetMarkerSize(0.7);
        h->SetMarkerColor(kBlue);
        h->SetLineColor(kBlue);
        h->SetLineWidth(1);
    }

private:
    TFile* f_all_ = nullptr;
    TFile* f_pp_  = nullptr;
    TFile* f_mm_  = nullptr;
    int counter_;

    void applyStyle() {
        gStyle->SetOptStat(0);
        gStyle->SetOptTitle(1);
        gStyle->SetTitleSize(0.05, "XYZ");
        gStyle->SetLabelSize(0.045, "XYZ");
        gStyle->SetTitleOffset(1.0, "X");
        gStyle->SetTitleOffset(1.1, "Y");
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);
    }

    std::string uniqueName(const std::string& base) {
        // Strip path separators for ROOT name safety
        std::string clean = base;
        for (auto& c : clean) {
            if (c == '/' || c == ' ' || c == '(' || c == ')') c = '_';
        }
        return clean + "_" + std::to_string(counter_++);
    }

    TH1* getHist(TFile* f, const std::string& path) {
        if (!f) return nullptr;
        return dynamic_cast<TH1*>(f->Get(path.c_str()));
    }

    /**
     * @brief Convert TH1F (from TNtuple) to TH1D for precise CB calculation
     */
    TH1D* toTH1D(TH1* h, const std::string& name) {
        if (auto* hd = dynamic_cast<TH1D*>(h)) {
            return (TH1D*)hd->Clone(name.c_str());
        }
        // Convert TH1F -> TH1D
        int nb = h->GetNbinsX();
        TH1D* hd = new TH1D(name.c_str(), h->GetTitle(),
                             nb, h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
        hd->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle());
        hd->GetYaxis()->SetTitle(h->GetYaxis()->GetTitle());
        for (int i = 0; i <= nb + 1; ++i) {
            hd->SetBinContent(i, h->GetBinContent(i));
            hd->SetBinError(i, h->GetBinError(i));
        }
        return hd;
    }

    /**
     * @brief Draw variable from ntuple in a single file
     */
    TH1D* drawFromNtuple(TFile* f, const std::string& ntuple_name,
                          const std::string& varexpr,
                          int nbins, double xmin, double xmax,
                          const std::string& cut,
                          const std::string& hname)
    {
        if (!f) return nullptr;

        TTree* tree = dynamic_cast<TTree*>(f->Get(ntuple_name.c_str()));
        if (!tree) {
            std::cerr << "PlotUtils: ntuple '" << ntuple_name
                      << "' not found in " << f->GetName() << "\n";
            return nullptr;
        }

        // Create histogram with desired binning
        TH1D* h = new TH1D(hname.c_str(), "", nbins, xmin, xmax);
        h->Sumw2();

        // Use TTree::Draw to fill it
        std::string drawcmd = varexpr + ">>" + hname;
        tree->Draw(drawcmd.c_str(), cut.c_str(), "goff");

        return h;
    }
};

#endif // PLOTUTILS_H
