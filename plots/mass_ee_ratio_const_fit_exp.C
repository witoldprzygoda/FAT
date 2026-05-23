// mass_ee_ratio_const_fit_exp.C — three ratio panels for M_ee SIGNAL with
// FIXED 20-bin binning over [0, 1.4] GeV/c², round-dot markers, vertical
// error bars only, and a constant fit y = a drawn on each panel.
//
// Standalone, self-contained — does NOT depend on AdaptiveSigPlot.h or
// PlotUtils.h. Reads ntuples from output_{epem,epep,emem}_exp.root.
//
// Three panels (all on signal = all − CB, all-channel files combined):
//   1) 63 · N_PT2 / N_PT3   (correction factor)
//   2) N_PT3 / (63 · N_PT2) (trigger efficiency)
//   3) N_PT3 / N_PT2        (raw firing-rate ratio, ref at 63)
//
// On each panel: red horizontal line is the fit  y = a  (constant), with
// the fitted value ± error and χ²/ndf annotated in red just above.
//
// Arguments (all optional):
//   1) oa_cut_deg  — OA threshold; events with oa > oa_cut_deg pass
//                    (default 0 ⇒ no cut, all events used)
//   2) use_weights — true → weighted χ² fit (1/σ² per bin, ROOT default)
//                    false → unweighted (all bins weight 1, simple mean)
//                    default: true
//   3) fit_xmin    — lower bound of the fit range in GeV/c²
//                    default −1 ⇒ use full histogram range
//   4) fit_xmax    — upper bound of the fit range in GeV/c²
//                    default −1 ⇒ use full histogram range
//
// The histogram is always drawn over [0, 1.4] GeV/c²; only the fit is
// restricted. Useful to exclude the π⁰ Dalitz peak bin (which dominates
// the weighted fit) or the high-mass tail (where bins fluctuate wildly).
//
// Usage:
//   root -l -b -q plots/mass_ee_ratio_const_fit_exp.C                       # default
//   root -l -b -q 'plots/mass_ee_ratio_const_fit_exp.C(4)'                  # OA > 4°
//   root -l -b -q 'plots/mass_ee_ratio_const_fit_exp.C(0, false)'           # unweighted
//   root -l -b -q 'plots/mass_ee_ratio_const_fit_exp.C(0, true, 0.07, 0.8)' # skip first bin and high-mass tail
//   root -l -b -q 'plots/mass_ee_ratio_const_fit_exp.C(4, true, 0.1, 0.7)'  # OA cut + restricted fit
//
// Output filename picks up suffixes "_oa{N}", "_unw", "_fit{lo}-{hi}"
// only when each option differs from default, so successive runs do not
// overwrite each other.
//
// Fit is weighted least squares (ROOT default): each bin contributes
// (y - a)² / σ² to χ², so points with larger error bars influence the
// fit less. Empty bins (σ = 0) are excluded by ROOT automatically.
//
// @author Witold Przygoda (witold.przygoda@uj.edu.pl)
// @date 2026

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

namespace {

TH1D* drawFromNt(TTree* t, const std::string& var, const std::string& cut,
                 const std::string& name, int nb, double xmin, double xmax) {
    TH1D* h = new TH1D(name.c_str(), "", nb, xmin, xmax);
    h->Sumw2();
    t->Draw((var + ">>" + name).c_str(), cut.c_str(), "goff");
    h->SetDirectory(nullptr);
    return h;
}

TH1D* makeCB(TH1D* pp, TH1D* mm, const std::string& name) {
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

TH1D* makeSig(TH1D* all, TH1D* cb, const std::string& name) {
    TH1D* sig = static_cast<TH1D*>(all->Clone(name.c_str()));
    sig->SetDirectory(nullptr);
    sig->Add(cb, -1.0);
    return sig;
}

TH1D* makeRatio(TH1D* num, TH1D* den, const std::string& name,
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

void styleDot(TH1D* h, Color_t color = kBlue + 1) {
    h->SetMarkerStyle(20);   // full circle
    h->SetMarkerSize(0.9);
    h->SetMarkerColor(color);
    h->SetLineColor(color);
    h->GetXaxis()->SetTitleSize(0.050);
    h->GetYaxis()->SetTitleSize(0.050);
    h->GetXaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetLabelSize(0.045);
    h->GetXaxis()->SetTitleOffset(1.05);
    h->GetYaxis()->SetTitleOffset(1.25);
}

struct FitResult { double a; double e; double chi2_ndf; int ndf; };

/// Fit y = a constant. With `use_weights = true` (ROOT default) the fit
/// minimises χ² = Σ (y_i − a)² / σ_i² — points with smaller error bars
/// dominate. With `use_weights = false` ("W" option) every non-empty
/// bin is set to weight 1, so the fit becomes a simple unweighted mean
/// that visually centres on the bulk of points.
FitResult fitConst(TH1D* h, double fit_xmin, double fit_xmax,
                   bool use_weights, bool draw, Color_t color) {
    const std::string fname = std::string(h->GetName()) +
                              (use_weights ? "_fconst_w" : "_fconst_u");
    // TF1 range == fit range; "R" option below makes Fit honour it.
    TF1* f = new TF1(fname.c_str(), "[0]", fit_xmin, fit_xmax);
    f->SetParName(0, "a");
    f->SetLineColor(color);
    f->SetLineWidth(2);
    f->SetLineStyle(1);
    const std::string opt =
        std::string("RQ") + (use_weights ? "" : "W") + (draw ? "" : "N0");
    h->Fit(f, opt.c_str());
    FitResult r;
    r.a        = f->GetParameter(0);
    r.e        = f->GetParError(0);
    r.ndf      = f->GetNDF();
    r.chi2_ndf = (r.ndf > 0) ? f->GetChisquare() / r.ndf : 0.0;
    return r;
}

}  // anonymous namespace

void mass_ee_ratio_const_fit_exp(double oa_cut_deg = 0.0,
                                 bool   use_weights = true,
                                 double fit_xmin = -1.0,
                                 double fit_xmax = -1.0) {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetErrorX(0.0);   // no horizontal error bars

    const int    nb   = 20;
    const double xmin = 0.0;
    const double xmax = 1.4;

    // Resolve fit range: sentinel −1 (or out-of-bounds) ⇒ full histogram range.
    const bool   custom_fit_range =
        (fit_xmin >= 0.0 && fit_xmin < xmax) ||
        (fit_xmax >  0.0 && fit_xmax <= xmax);
    const double fxlo = (fit_xmin >= 0.0 && fit_xmin < xmax) ? fit_xmin : xmin;
    const double fxhi = (fit_xmax >  0.0 && fit_xmax <= xmax) ? fit_xmax : xmax;

    TFile* fe = TFile::Open("output_epem_exp.root", "READ");
    TFile* fp = TFile::Open("output_epep_exp.root", "READ");
    TFile* fm = TFile::Open("output_emem_exp.root", "READ");
    if (!fe || !fp || !fm) {
        std::cerr << "Cannot open one of the input files\n";
        return;
    }
    TTree* nt_em = dynamic_cast<TTree*>(fe->Get("dilepton_nt"));
    TTree* nt_pp = dynamic_cast<TTree*>(fp->Get("dilepton_nt"));
    TTree* nt_mm = dynamic_cast<TTree*>(fm->Get("dilepton_nt"));
    if (!nt_em || !nt_pp || !nt_mm) {
        std::cerr << "dilepton_nt missing in one of the files\n";
        return;
    }

    // --- Compose cut strings. OA cut is "oa > oa_cut_deg"; default 0 ⇒
    //     no OA filtering (every event with the matching trigbit passes).
    const std::string oa_clause =
        (oa_cut_deg > 0.0) ? std::string(" && oa>") + std::to_string(oa_cut_deg)
                           : std::string{};
    const std::string cut_pt3 = "trigbit==8192" + oa_clause;
    const std::string cut_pt2 = "trigbit==4096" + oa_clause;
    std::cout << "Cuts:\n  PT3 = " << cut_pt3
              << "\n  PT2 = " << cut_pt2
              << "\nFit mode (drawn): " << (use_weights ? "WEIGHTED (uses bin errors as 1/sigma^2)"
                                                        : "UNWEIGHTED (all bins weight 1, simple mean)")
              << "\nFit range:        [" << fxlo << ", " << fxhi << "] GeV/c^2"
              << (custom_fit_range ? "  (custom)" : "  (full)")
              << "\n";

    // --- Build per-channel histograms with fixed 20 bins -------------
    TH1D* em_p3 = drawFromNt(nt_em, "m_ee", cut_pt3, "em_p3_fit", nb, xmin, xmax);
    TH1D* pp_p3 = drawFromNt(nt_pp, "m_ee", cut_pt3, "pp_p3_fit", nb, xmin, xmax);
    TH1D* mm_p3 = drawFromNt(nt_mm, "m_ee", cut_pt3, "mm_p3_fit", nb, xmin, xmax);
    TH1D* em_p2 = drawFromNt(nt_em, "m_ee", cut_pt2, "em_p2_fit", nb, xmin, xmax);
    TH1D* pp_p2 = drawFromNt(nt_pp, "m_ee", cut_pt2, "pp_p2_fit", nb, xmin, xmax);
    TH1D* mm_p2 = drawFromNt(nt_mm, "m_ee", cut_pt2, "mm_p2_fit", nb, xmin, xmax);

    TH1D* cb_p3  = makeCB (pp_p3, mm_p3, "cb_p3_fit");
    TH1D* cb_p2  = makeCB (pp_p2, mm_p2, "cb_p2_fit");
    TH1D* sig_p3 = makeSig(em_p3, cb_p3, "sig_p3_fit");
    TH1D* sig_p2 = makeSig(em_p2, cb_p2, "sig_p2_fit");

    // --- Three ratios on signal only ---------------------------------
    TH1D* r63  = makeRatio(sig_p2, sig_p3, "r63_fit",  63.0);
    TH1D* reff = makeRatio(sig_p3, sig_p2, "reff_fit", 1.0 / 63.0);
    TH1D* rinv = makeRatio(sig_p3, sig_p2, "rinv_fit", 1.0);

    styleDot(r63);  styleDot(reff);  styleDot(rinv);

    // OA tag injected into titles only when a cut is applied.
    const std::string oa_tag =
        (oa_cut_deg > 0.0)
            ? Form(", OA > %g#circ", oa_cut_deg)
            : std::string{};
    r63 ->SetTitle((std::string(
        "63 #upoint N_{PT2}/N_{PT3}  signal (20 fixed bins") + oa_tag + ");"
        "M_{e^{+}e^{-}} [GeV/c^{2}];63 #upoint N_{PT2}/N_{PT3} (sig)").c_str());
    reff->SetTitle((std::string(
        "N_{PT3} / (63 #upoint N_{PT2})  signal (20 fixed bins") + oa_tag + ");"
        "M_{e^{+}e^{-}} [GeV/c^{2}];N_{PT3} / (63 #upoint N_{PT2}) (sig)").c_str());
    rinv->SetTitle((std::string(
        "N_{PT3} / N_{PT2}  signal (20 fixed bins") + oa_tag + ");"
        "M_{e^{+}e^{-}} [GeV/c^{2}];N_{PT3} / N_{PT2} (sig)").c_str());

    // --- Canvas: 3 panels --------------------------------------------
    TCanvas* c = new TCanvas("c_mass_ee_ratio_const_fit",
                             "M_{ee} ratio constant fit (20 fixed bins)",
                             1800, 600);
    c->Divide(3, 1, 0.001, 0.001);

    auto drawPad = [&](int idx, TH1D* h, double ref_y,
                       const std::string& /*ref_label*/,
                       int ref_style, int ref_color, double y_max_floor) {
        c->cd(idx);
        const double pad_lm = 0.14;   // left margin (NDC)
        const double pad_rm = 0.04;
        const double pad_bm = 0.13;
        const double pad_tm = 0.10;
        gPad->SetMargin(pad_lm, pad_rm, pad_bm, pad_tm);

        // Auto Y range from data (+25% margin, never below y_max_floor).
        double ymax = 0.0;
        for (int b = 1; b <= nb; ++b) {
            const double v = h->GetBinContent(b) + h->GetBinError(b);
            if (std::isfinite(v) && v > ymax) ymax = v;
        }
        const double y_hi = std::max(ymax * 1.25, y_max_floor);
        h->GetYaxis()->SetRangeUser(0.0, y_hi);

        h->Draw("E1");

        // Reference line (drawn under fit).
        TLine* lref = new TLine(xmin, ref_y, xmax, ref_y);
        lref->SetLineStyle(ref_style);
        lref->SetLineColor(ref_color);
        if (ref_style == 2) lref->SetLineWidth(2);
        lref->Draw();

        // Constant fit — drawn version respects user choice; we also
        // compute the OTHER option silently and print both to stdout.
        const FitResult fit  = fitConst(h, fxlo, fxhi, use_weights,
                                        /*draw=*/true,  kRed + 1);
        const FitResult fit2 = fitConst(h, fxlo, fxhi, !use_weights,
                                        /*draw=*/false, kRed + 1);

        // Annotate "a = ... ± ..." JUST ABOVE the fit line itself. Compute
        // NDC y from the fit value; flip to "below the line" if the line
        // is so high in the pad that "above" would clip the top margin.
        const double y_frac = (y_hi > 0) ? (fit.a / y_hi) : 0.0;
        const double y_line_ndc =
            pad_bm + (1.0 - pad_bm - pad_tm) * y_frac;
        const double y_text =
            (y_line_ndc < 0.70) ? y_line_ndc + 0.045 : y_line_ndc - 0.075;

        TLatex tex;
        tex.SetNDC();
        tex.SetTextAlign(11);   // left-bottom
        tex.SetTextSize(0.050);
        tex.SetTextColor(kRed + 1);
        tex.DrawLatex(0.52, y_text,
            Form("a = %.4f #pm %.4f", fit.a, fit.e));
        tex.SetTextSize(0.035);
        tex.SetTextColor(kGray + 3);
        tex.DrawLatex(0.52, y_text - 0.045,
            Form("#chi^{2}/ndf = %.1f / %d",
                 fit.chi2_ndf * fit.ndf, fit.ndf));

        const char* tag_drawn = use_weights ? "weighted  " : "unweighted";
        const char* tag_other = use_weights ? "unweighted" : "weighted  ";
        std::cout << "  pad " << idx
                  << "  [" << tag_drawn << " drawn] a = " << fit.a
                  << " ± " << fit.e
                  << "   chi2/ndf = " << (fit.chi2_ndf * fit.ndf)
                  << "/" << fit.ndf << "\n"
                  << "          [" << tag_other << "]       a = " << fit2.a
                  << " ± " << fit2.e
                  << "   chi2/ndf = " << (fit2.chi2_ndf * fit2.ndf)
                  << "/" << fit2.ndf << "\n";
    };

    std::cout << "Constant fit y = a:\n";
    drawPad(1, r63,   1.0, "1.0 (no correction)",        3, kGray + 2,  5.0);
    drawPad(2, reff,  1.0, "1.0 (PT3 catches all PT2)",  3, kGray + 2,  1.2);
    drawPad(3, rinv, 63.0, "63 (no-bias reference)",     2, kGray + 3, 70.0);

    gSystem->mkdir("plots/output", true);
    const std::string oa_suffix =
        (oa_cut_deg > 0.0) ? Form("_oa%g", oa_cut_deg) : std::string{};
    const std::string wt_suffix = use_weights ? std::string{} : "_unw";
    const std::string fr_suffix =
        custom_fit_range ? Form("_fit%g-%g", fxlo, fxhi) : std::string{};
    const std::string base =
        "plots/output/mass_ee_ratio_const_fit" + oa_suffix + wt_suffix + fr_suffix;
    c->SaveAs((base + ".pdf").c_str());
    c->SaveAs((base + ".png").c_str());
    std::cout << "Saved: " << base << ".{pdf,png}\n";

    fe->Close(); fp->Close(); fm->Close();
}
