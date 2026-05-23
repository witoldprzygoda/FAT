// AdaptiveSigPlot.h — shared logic for adaptive-binning CB-subtracted plots.
//
// A per-variable macro fills an `AdaptiveSigPlot::Config` and calls
// `AdaptiveSigPlot::run(cfg)`. The run function produces TWO canvases:
//
//   1) <base>_adaptive_sig.{pdf,png}             — three panels:
//      • PT3 spectra (all / CB / signal)          adaptive bins, log Y
//      • PT2 spectra (all / CB / signal)          adaptive bins, log Y
//      • 63 · N_PT2 / N_PT3 (signal only)        same bins, linear Y
//
//   2) <base>_adaptive_sig_ratios_coarse.{pdf,png} — three panels (signal only):
//      • 63 · N_PT2 / N_PT3   (correction factor)
//      • N_PT3 / (63 · N_PT2) (trigger efficiency, ≡ 1/<w>)
//      • N_PT3 / N_PT2        (raw firing-rate ratio, ref. line at 63)
//
// Adaptive binning: walk a fine-grained PT2 epem histogram, close each
// bin when accumulated PT2 count reaches `N_min_pt2_spectra` (or
// `N_min_pt2_ratio` for the coarser canvas) OR width hits
// `max_bin_width`. Histogram contents are then divided by bin width
// (counts / unit) so a variable-width binning shows no spurious steps.
//
// To extend to a new observable, write a thin wrapper macro:
//
//     #include "AdaptiveSigPlot.h"
//     void pt_cms_adaptive_sig_exp(double N1 = 200, double mw = 150,
//                                  double N2 = 2000) {
//         AdaptiveSigPlot::Config c;
//         c.var_expr           = "pt";
//         c.x_axis_label       = "p_{T} [MeV/c]";
//         c.y_axis_unit        = "MeV/c";
//         c.xmin = 0.0;  c.xmax = 1500.0;  c.fine_nbins = 1500;
//         c.max_bin_width      = mw;
//         c.N_min_pt2_spectra  = N1;
//         c.N_min_pt2_ratio    = N2;
//         c.base_filename      = "pt_cms";
//         c.canvas_title       = "p_{T} (CMS)";
//         AdaptiveSigPlot::run(c);
//     }
//
// @author Witold Przygoda (witold.przygoda@uj.edu.pl)
// @date 2026

#ifndef ADAPTIVE_SIG_PLOT_H
#define ADAPTIVE_SIG_PLOT_H

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace AdaptiveSigPlot {

// ---------------------------------------------------------------------------
// Configuration for one variable
// ---------------------------------------------------------------------------
struct Config {
    std::string var_expr           = "m_ee";
    std::string x_axis_label       = "M_{e^{+}e^{-}} [GeV/c^{2}]";
    std::string y_axis_unit        = "GeV/c^{2}";  // "" → just "Counts/bin"
    double      xmin               = 0.0;
    double      xmax               = 1.4;
    int         fine_nbins         = 1400;
    double      max_bin_width      = 0.2;
    double      N_min_pt2_spectra  = 200.0;
    double      N_min_pt2_ratio    = 2000.0;
    std::string base_filename      = "mass_ee";
    std::string canvas_title       = "M_{e^{+}e^{-}}";
    std::string nt_name            = "dilepton_nt";
    bool        log_y_spectra      = true;
};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------
inline TH1D* drawFromNt(TTree* t, const std::string& expr,
                        const std::string& cut, const std::string& name,
                        int nb, const double* edges) {
    TH1D* h = new TH1D(name.c_str(), "", nb, edges);
    h->Sumw2();
    t->Draw((expr + ">>" + name).c_str(), cut.c_str(), "goff");
    h->SetDirectory(nullptr);
    return h;
}

inline TH1D* makeCB(TH1D* pp, TH1D* mm, const std::string& name) {
    TH1D* cb = static_cast<TH1D*>(pp->Clone(name.c_str()));
    cb->SetDirectory(nullptr);
    cb->Reset();
    const int nb = pp->GetNbinsX();
    for (int b = 1; b <= nb; ++b) {
        const double np = pp->GetBinContent(b);
        const double ep = pp->GetBinError(b);
        const double nm = mm->GetBinContent(b);
        const double em = mm->GetBinError(b);
        if (np > 0 && nm > 0) {
            const double cval = 2.0 * std::sqrt(np * nm);
            const double rel  = std::sqrt(std::pow(ep / np, 2) +
                                           std::pow(em / nm, 2));
            cb->SetBinContent(b, cval);
            cb->SetBinError  (b, cval * 0.5 * rel);
        }
    }
    return cb;
}

inline TH1D* makeSig(TH1D* all, TH1D* cb, const std::string& name) {
    TH1D* sig = static_cast<TH1D*>(all->Clone(name.c_str()));
    sig->SetDirectory(nullptr);
    sig->Add(cb, -1.0);
    return sig;
}

/// Per-bin ratio scale · num/den with independent-variable error propagation.
inline TH1D* makeRatio(TH1D* num, TH1D* den, const std::string& name,
                       double scale = 1.0) {
    TH1D* r = static_cast<TH1D*>(num->Clone(name.c_str()));
    r->SetDirectory(nullptr);
    r->Reset();
    const int nb = num->GetNbinsX();
    for (int b = 1; b <= nb; ++b) {
        const double n  = num->GetBinContent(b);
        const double en = num->GetBinError(b);
        const double d  = den->GetBinContent(b);
        const double ed = den->GetBinError(b);
        if (d != 0.0 && std::isfinite(d) && std::isfinite(n)) {
            const double val = scale * n / d;
            const double rel = std::sqrt(
                (n != 0.0 ? std::pow(en / n, 2) : 0.0) +
                std::pow(ed / d, 2));
            r->SetBinContent(b, val);
            r->SetBinError  (b, std::abs(val) * rel);
        }
    }
    return r;
}

inline double safeMaxAbs(TH1D* h, bool include_err = true) {
    double mx = 0.0;
    const int nb = h->GetNbinsX();
    for (int b = 1; b <= nb; ++b) {
        const double v = std::abs(h->GetBinContent(b))
                       + (include_err ? h->GetBinError(b) : 0.0);
        if (std::isfinite(v) && v > mx) mx = v;
    }
    return mx;
}

/// Robust max — ignore bins with relative error > `rel_err_cap` (outliers
/// in low-stat tails that otherwise blow up the auto Y range).
inline double safeMaxRobust(TH1D* h, double rel_err_cap = 0.50) {
    double mx = 0.0;
    const int nb = h->GetNbinsX();
    for (int b = 1; b <= nb; ++b) {
        const double v = h->GetBinContent(b);
        const double e = h->GetBinError(b);
        if (!std::isfinite(v) || v <= 0.0) continue;
        if (e > 0.0 && e / std::abs(v) > rel_err_cap) continue;
        if (v + e > mx) mx = v + e;
    }
    return mx;
}

inline double safeMinPositive(TH1D* h) {
    double mn = std::numeric_limits<double>::infinity();
    const int nb = h->GetNbinsX();
    for (int b = 1; b <= nb; ++b) {
        const double v = h->GetBinContent(b);
        if (std::isfinite(v) && v > 0.0 && v < mn) mn = v;
    }
    return std::isfinite(mn) ? mn : 1.0;
}

inline void styleAll(TH1D* h) {
    h->SetMarkerStyle(20); h->SetMarkerSize(0.65);
    h->SetMarkerColor(kBlack); h->SetLineColor(kBlack);
}
inline void styleCB(TH1D* h) {
    h->SetMarkerStyle(21); h->SetMarkerSize(0.65);
    h->SetMarkerColor(kRed + 1); h->SetLineColor(kRed + 1);
}
inline void styleSig(TH1D* h) {
    h->SetMarkerStyle(22); h->SetMarkerSize(0.7);
    h->SetMarkerColor(kBlue + 1); h->SetLineColor(kBlue + 1);
}
inline void styleAxes(TH1D* h) {
    h->GetXaxis()->SetTitleSize(0.050);
    h->GetYaxis()->SetTitleSize(0.050);
    h->GetXaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetLabelSize(0.045);
    h->GetXaxis()->SetTitleOffset(1.05);
    h->GetYaxis()->SetTitleOffset(1.20);
}

/// Derive adaptive bin edges from epem PT2 distribution.
inline std::vector<double> deriveEdges(TTree* nt_em, const Config& cfg,
                                       double N_min,
                                       const std::string& seed_name) {
    const double fine_w = (cfg.xmax - cfg.xmin) / cfg.fine_nbins;
    TH1D* h_seed = new TH1D(seed_name.c_str(), "",
                            cfg.fine_nbins, cfg.xmin, cfg.xmax);
    nt_em->Draw((cfg.var_expr + ">>" + seed_name).c_str(),
                "trigbit==4096", "goff");
    h_seed->SetDirectory(nullptr);

    std::vector<double> edges;
    edges.reserve(cfg.fine_nbins + 1);
    edges.push_back(cfg.xmin);
    double acc = 0.0;
    for (int b = 1; b <= cfg.fine_nbins; ++b) {
        acc += h_seed->GetBinContent(b);
        const double upper = h_seed->GetBinLowEdge(b + 1);
        const double width = upper - edges.back();
        if (acc >= N_min || width >= cfg.max_bin_width - 0.5 * fine_w) {
            edges.push_back(upper);
            acc = 0.0;
        }
    }
    if (edges.back() < cfg.xmax - 0.5 * fine_w) edges.push_back(cfg.xmax);
    delete h_seed;
    return edges;
}

// ---------------------------------------------------------------------------
// Main entry point
// ---------------------------------------------------------------------------
inline void run(const Config& cfg) {
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    TFile* fe = TFile::Open("output_epem_exp.root", "READ");
    TFile* fp = TFile::Open("output_epep_exp.root", "READ");
    TFile* fm = TFile::Open("output_emem_exp.root", "READ");
    if (!fe || !fp || !fm) {
        std::cerr << "AdaptiveSigPlot::run — cannot open one of the input files\n";
        return;
    }
    TTree* nt_em = dynamic_cast<TTree*>(fe->Get(cfg.nt_name.c_str()));
    TTree* nt_pp = dynamic_cast<TTree*>(fp->Get(cfg.nt_name.c_str()));
    TTree* nt_mm = dynamic_cast<TTree*>(fm->Get(cfg.nt_name.c_str()));
    if (!nt_em || !nt_pp || !nt_mm) {
        std::cerr << "AdaptiveSigPlot::run — tree '" << cfg.nt_name << "' missing\n";
        return;
    }

    const std::string seed_fine_name   = "h_seed_pt2_" + cfg.base_filename + "_f";
    const std::string seed_coarse_name = "h_seed_pt2_" + cfg.base_filename + "_c";

    // -------- FINE binning for spectra --------
    const std::vector<double> edges = deriveEdges(
        nt_em, cfg, cfg.N_min_pt2_spectra, seed_fine_name);
    const int nb = static_cast<int>(edges.size()) - 1;
    std::cout << "[" << cfg.base_filename << "] spectra: "
              << nb << " bins, N_min_pt2_spectra=" << cfg.N_min_pt2_spectra
              << ", max_width=" << cfg.max_bin_width << "\n";

    // -------- Per-channel histograms (raw counts) -----------------
    const std::string n0 = cfg.base_filename + "_";
    TH1D* em_p3 = drawFromNt(nt_em, cfg.var_expr, "trigbit==8192", n0 + "em_p3", nb, edges.data());
    TH1D* pp_p3 = drawFromNt(nt_pp, cfg.var_expr, "trigbit==8192", n0 + "pp_p3", nb, edges.data());
    TH1D* mm_p3 = drawFromNt(nt_mm, cfg.var_expr, "trigbit==8192", n0 + "mm_p3", nb, edges.data());
    TH1D* em_p2 = drawFromNt(nt_em, cfg.var_expr, "trigbit==4096", n0 + "em_p2", nb, edges.data());
    TH1D* pp_p2 = drawFromNt(nt_pp, cfg.var_expr, "trigbit==4096", n0 + "pp_p2", nb, edges.data());
    TH1D* mm_p2 = drawFromNt(nt_mm, cfg.var_expr, "trigbit==4096", n0 + "mm_p2", nb, edges.data());

    TH1D* cb_p3  = makeCB (pp_p3, mm_p3, n0 + "cb_p3");
    TH1D* cb_p2  = makeCB (pp_p2, mm_p2, n0 + "cb_p2");
    TH1D* sig_p3 = makeSig(em_p3, cb_p3, n0 + "sig_p3");
    TH1D* sig_p2 = makeSig(em_p2, cb_p2, n0 + "sig_p2");

    TH1D* r_sig  = makeRatio(sig_p2, sig_p3, n0 + "r_sig", 63.0);

    // Scale spectra by bin width — counts → counts / unit
    for (TH1D* h : {em_p3, cb_p3, sig_p3, em_p2, cb_p2, sig_p2}) {
        h->Scale(1.0, "width");
    }

    // Styling
    styleAll(em_p3); styleCB(cb_p3); styleSig(sig_p3);
    styleAll(em_p2); styleCB(cb_p2); styleSig(sig_p2);
    styleSig(r_sig);
    for (TH1D* h : {em_p3, cb_p3, sig_p3, em_p2, cb_p2, sig_p2, r_sig}) styleAxes(h);

    // Axis titles
    const std::string ylabel = cfg.y_axis_unit.empty()
        ? "Counts / bin"
        : ("Counts / (" + cfg.y_axis_unit + ")");
    em_p3->SetTitle(("PT3: all / CB / signal  (adaptive bins);"
                     + cfg.x_axis_label + ";" + ylabel).c_str());
    em_p2->SetTitle(("PT2: all / CB / signal  (adaptive bins);"
                     + cfg.x_axis_label + ";" + ylabel).c_str());
    r_sig->SetTitle(("63 #upoint N_{PT2}/N_{PT3}  (signal);"
                     + cfg.x_axis_label + ";"
                     "63 #upoint N_{PT2}/N_{PT3}  (sig)").c_str());

    // -------- Canvas 1 (spectra + fine ratio) --------
    TCanvas* c = new TCanvas(("c_" + cfg.base_filename + "_adaptive_sig").c_str(),
                             ("Adaptive " + cfg.canvas_title + ": PT3 / PT2 / ratio").c_str(),
                             1800, 600);
    c->Divide(3, 1, 0.001, 0.001);

    auto drawSpectraPad = [&](int pad_idx, TH1D* all, TH1D* cb, TH1D* sig) {
        c->cd(pad_idx);
        if (cfg.log_y_spectra) gPad->SetLogy();
        gPad->SetMargin(0.14, 0.04, 0.13, 0.10);
        const double ymax = std::max({safeMaxAbs(all), safeMaxAbs(cb), safeMaxAbs(sig)});
        const double ymin = std::max(0.5, std::min({safeMinPositive(all),
                                                    safeMinPositive(cb),
                                                    safeMinPositive(sig)}));
        if (cfg.log_y_spectra) {
            all->GetYaxis()->SetRangeUser(ymin * 0.5, ymax * 4.0);
        } else {
            all->GetYaxis()->SetRangeUser(0.0, ymax * 1.25);
        }
        all->Draw("E1");
        cb ->Draw("E1 SAME");
        sig->Draw("E1 SAME");
        TLegend* leg = new TLegend(0.55, 0.74, 0.95, 0.90);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.040);
        leg->AddEntry(all, "e^{+}e^{-} (all)",                 "lpe");
        leg->AddEntry(cb,  "CB = 2#sqrt{N_{++}#upointN_{--}}", "lpe");
        leg->AddEntry(sig, "Signal (all - CB)",                "lpe");
        leg->Draw();
    };
    drawSpectraPad(1, em_p3, cb_p3, sig_p3);
    drawSpectraPad(2, em_p2, cb_p2, sig_p2);

    // Pad 3: fine ratio — signal only
    c->cd(3);
    gPad->SetMargin(0.14, 0.04, 0.13, 0.10);
    {
        const double ymax = safeMaxRobust(r_sig);
        const double y_hi = (ymax > 0.0 && std::isfinite(ymax)) ? ymax * 1.25 : 5.0;
        r_sig->GetYaxis()->SetRangeUser(0.0, y_hi);
    }
    r_sig->Draw("E1");
    TLine* l_one = new TLine(cfg.xmin, 1.0, cfg.xmax, 1.0);
    l_one->SetLineStyle(3); l_one->SetLineColor(kGray + 2); l_one->Draw();
    {
        TLegend* leg = new TLegend(0.50, 0.78, 0.95, 0.90);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.038);
        leg->AddEntry(r_sig, "63 #upoint N_{PT2}/N_{PT3}  (signal)", "lpe");
        leg->AddEntry(l_one, "1.0 (no correction)",                  "l");
        leg->Draw();
    }

    gSystem->mkdir("plots/output", true);
    const std::string out1 = "plots/output/" + cfg.base_filename + "_adaptive_sig";
    c->SaveAs((out1 + ".pdf").c_str());
    c->SaveAs((out1 + ".png").c_str());
    std::cout << "Saved: " << out1 << ".{pdf,png}\n";

    // ===================================================================
    // COARSE binning for ratios — signal only
    // ===================================================================
    const std::vector<double> edges_c = deriveEdges(
        nt_em, cfg, cfg.N_min_pt2_ratio, seed_coarse_name);
    const int nb_c = static_cast<int>(edges_c.size()) - 1;
    std::cout << "[" << cfg.base_filename << "] coarse: "
              << nb_c << " bins, N_min_pt2_ratio=" << cfg.N_min_pt2_ratio << "\n";

    TH1D* em_p3c = drawFromNt(nt_em, cfg.var_expr, "trigbit==8192", n0 + "em_p3c", nb_c, edges_c.data());
    TH1D* pp_p3c = drawFromNt(nt_pp, cfg.var_expr, "trigbit==8192", n0 + "pp_p3c", nb_c, edges_c.data());
    TH1D* mm_p3c = drawFromNt(nt_mm, cfg.var_expr, "trigbit==8192", n0 + "mm_p3c", nb_c, edges_c.data());
    TH1D* em_p2c = drawFromNt(nt_em, cfg.var_expr, "trigbit==4096", n0 + "em_p2c", nb_c, edges_c.data());
    TH1D* pp_p2c = drawFromNt(nt_pp, cfg.var_expr, "trigbit==4096", n0 + "pp_p2c", nb_c, edges_c.data());
    TH1D* mm_p2c = drawFromNt(nt_mm, cfg.var_expr, "trigbit==4096", n0 + "mm_p2c", nb_c, edges_c.data());
    TH1D* cb_p3c  = makeCB (pp_p3c, mm_p3c, n0 + "cb_p3c");
    TH1D* cb_p2c  = makeCB (pp_p2c, mm_p2c, n0 + "cb_p2c");
    TH1D* sig_p3c = makeSig(em_p3c, cb_p3c, n0 + "sig_p3c");
    TH1D* sig_p2c = makeSig(em_p2c, cb_p2c, n0 + "sig_p2c");

    // Three ratio variants — all on the signal only.
    TH1D* r63_sig_c  = makeRatio(sig_p2c, sig_p3c, n0 + "r63_sig_c",  63.0);
    TH1D* reff_sig_c = makeRatio(sig_p3c, sig_p2c, n0 + "reff_sig_c", 1.0 / 63.0);
    TH1D* rinv_sig_c = makeRatio(sig_p3c, sig_p2c, n0 + "rinv_sig_c", 1.0);

    styleSig(r63_sig_c); styleSig(reff_sig_c); styleSig(rinv_sig_c);
    styleAxes(r63_sig_c); styleAxes(reff_sig_c); styleAxes(rinv_sig_c);

    TCanvas* c2 = new TCanvas(
        ("c_" + cfg.base_filename + "_adaptive_sig_ratios_coarse").c_str(),
        ("Coarse-binning ratios: 63·PT2/PT3, PT3/(63·PT2), PT3/PT2 — " + cfg.canvas_title).c_str(),
        2100, 600);
    c2->Divide(3, 1, 0.001, 0.001);

    auto drawRatioPad = [&](int idx, TH1D* h, const std::string& title,
                            const std::string& ylab, const std::string& leg_label,
                            double ref_y, const std::string& ref_label,
                            int ref_style, int ref_color, double y_max_floor) {
        c2->cd(idx);
        gPad->SetMargin(0.14, 0.04, 0.13, 0.10);
        const double ymax = safeMaxRobust(h);
        const double y_hi = (ymax > 0.0 && std::isfinite(ymax)) ? ymax * 1.25 : y_max_floor;
        h->SetTitle((title + ";" + cfg.x_axis_label + ";" + ylab).c_str());
        h->GetYaxis()->SetRangeUser(0.0, std::max(y_hi, y_max_floor));
        h->Draw("E1");
        TLine* lref = new TLine(cfg.xmin, ref_y, cfg.xmax, ref_y);
        lref->SetLineStyle(ref_style); lref->SetLineColor(ref_color);
        if (ref_style == 2) lref->SetLineWidth(2);
        lref->Draw();
        TLegend* leg = new TLegend(0.40, 0.78, 0.95, 0.90);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.038);
        leg->AddEntry(h,    leg_label.c_str(), "lpe");
        leg->AddEntry(lref, ref_label.c_str(), "l");
        leg->Draw();
    };

    drawRatioPad(1, r63_sig_c,
        "63 #upoint N_{PT2}/N_{PT3}  signal (coarse bins)",
        "63 #upoint N_{PT2}/N_{PT3}  (sig)",
        "63 #upoint N_{PT2}/N_{PT3}  (signal)",
        1.0, "1.0 (no correction)", 3, kGray + 2, 5.0);

    drawRatioPad(2, reff_sig_c,
        "N_{PT3} / (63 #upoint N_{PT2})  signal (coarse bins)",
        "N_{PT3} / (63 #upoint N_{PT2})  (sig)",
        "N_{PT3} / (63 #upoint N_{PT2})  (signal)",
        1.0, "1.0 (PT3 catches all PT2)", 3, kGray + 2, 1.2);

    drawRatioPad(3, rinv_sig_c,
        "N_{PT3} / N_{PT2}  signal (coarse bins)",
        "N_{PT3} / N_{PT2}  (sig)",
        "N_{PT3} / N_{PT2}  (signal)",
        63.0, "63 (= no-bias, PT2 downscale 64#minus1)", 2, kGray + 3, 70.0);

    const std::string out2 = "plots/output/" + cfg.base_filename + "_adaptive_sig_ratios_coarse";
    c2->SaveAs((out2 + ".pdf").c_str());
    c2->SaveAs((out2 + ".png").c_str());
    std::cout << "Saved: " << out2 << ".{pdf,png}\n";

    // -------- Integrals (raw counts) --------
    auto rawIntegral = [&](TH1D* h) {
        double s = 0.0;
        const int n = h->GetNbinsX();
        for (int b = 1; b <= n; ++b) s += h->GetBinContent(b) * h->GetBinWidth(b);
        return s;
    };
    std::cout << "[" << cfg.base_filename << "] PT3 integrals (raw): "
              << "all=" << rawIntegral(em_p3)
              << "  CB=" << rawIntegral(cb_p3)
              << "  sig=" << rawIntegral(sig_p3) << "\n";
    std::cout << "[" << cfg.base_filename << "] PT2 integrals (raw): "
              << "all=" << rawIntegral(em_p2)
              << "  CB=" << rawIntegral(cb_p2)
              << "  sig=" << rawIntegral(sig_p2) << "\n";

    fe->Close(); fp->Close(); fm->Close();
}

}  // namespace AdaptiveSigPlot

#endif  // ADAPTIVE_SIG_PLOT_H
