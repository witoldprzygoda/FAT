// mass_ee_ratio_const_fit_exp_vs_sim.C — overlay of EXP data and SMASH SIM
// ratios in the same canvas, with separate constant fits per sample using
// SHARED fit parameters (range, OA cut, weighted/unweighted flag).
//
// Two panels (signal only, 20 fixed bins on [0, 1.4] GeV/c²):
//
//   Pad 1 — "Trigger correction factor":
//     exp: 63 · N_PT2 / N_PT3    (factor 63 = PT2 downscale 64−1)
//     sim:      N_PT2 / N_PT3    (no downscale in simulation)
//     Both represent the multiplier to apply to PT3-observed counts to
//     recover the unbiased event count. Exp ≈ 2.5, sim ≈ 1.6 for pp45.
//
//   Pad 2 — "Trigger efficiency":
//     exp: N_PT3 / (63 · N_PT2)  (= 1/pad1 for exp)
//     sim: N_PT3 /       N_PT2   (= 1/pad1 for sim)
//     Both represent the fraction of unbiased events that PT3 catches.
//     Exp ≈ 0.40, sim ≈ 0.61 for pp45 (sim is higher because the SMASH
//     sample is filtered to true dileptons, while exp data has CB and
//     mixed lepton/hadron rejection).
//
// On each panel both data points and their constant fits are drawn:
//   exp = blue filled circles + blue solid fit line
//   sim = red  filled squares + red  solid fit line
// Fit values ± error and χ²/ndf annotated in matching colors at the top
// of the pad (exp on the left, sim on the right).
//
// EXP signal: CB-subtracted from three input files (output_{epem,epep,emem}
//             _exp.root) using like-sign 2·√(N_++ · N_--).
// SIM signal: dilepton_nt from output_epem_sim.root, weighted by
//             sim_genweight (SMASH luminosity weight); no CB needed because
//             the MC purity gate already selects true e+e- pairs.
//
// Arguments (all optional, identical semantics to the per-sample macros):
//   1) oa_cut_deg   — OA threshold; default 0 ⇒ no cut
//   2) use_weights  — true → weighted χ² fit; false → unweighted; default true
//   3) fit_xmin     — fit range lower bound [GeV/c²]; -1 ⇒ full
//   4) fit_xmax     — fit range upper bound [GeV/c²]; -1 ⇒ full
//
// Usage:
//   root -l -b -q plots/mass_ee_ratio_const_fit_exp_vs_sim.C
//   root -l -b -q 'plots/mass_ee_ratio_const_fit_exp_vs_sim.C(4)'
//   root -l -b -q 'plots/mass_ee_ratio_const_fit_exp_vs_sim.C(0, true, 0.1, 0.8)'
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
#include <limits>
#include <string>

namespace {

constexpr int    kNb   = 20;
constexpr double kXmin = 0.0;
constexpr double kXmax = 1.4;

constexpr Color_t kColExp     = kBlue + 1;
constexpr Color_t kColExpFit  = kBlue + 2;
constexpr Color_t kColSim     = kRed  + 1;
constexpr Color_t kColSimFit  = kRed  + 2;

TH1D* drawFromNt(TTree* t, const std::string& cut_or_weight,
                 const std::string& name) {
    TH1D* h = new TH1D(name.c_str(), "", kNb, kXmin, kXmax);
    h->Sumw2();
    t->Draw(("m_ee>>" + name).c_str(), cut_or_weight.c_str(), "goff");
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

TH1D* makeRatio(TH1D* num, TH1D* den, const std::string& name, double scale) {
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

/// Scaled sim / exp ratio with full error propagation. Scale uncertainty
/// (from the two constant fits) propagates through as well — added in
/// quadrature with the per-bin sim and exp relative errors.
TH1D* makeScaledSimOverExp(TH1D* sim, TH1D* exp_h, double scale,
                           double scale_err, const std::string& name) {
    TH1D* r = static_cast<TH1D*>(sim->Clone(name.c_str()));
    r->SetDirectory(nullptr);
    r->Reset();
    const int nb = sim->GetNbinsX();
    const double rel_scale =
        (scale != 0.0) ? std::abs(scale_err / scale) : 0.0;
    for (int b = 1; b <= nb; ++b) {
        const double s  = sim->GetBinContent(b);
        const double es = sim->GetBinError(b);
        const double e  = exp_h->GetBinContent(b);
        const double ee = exp_h->GetBinError(b);
        if (e != 0.0 && std::isfinite(e) && std::isfinite(s)) {
            const double val = scale * s / e;
            const double rel = std::sqrt(
                (s != 0.0 ? std::pow(es / s, 2) : 0.0) +
                std::pow(ee / e, 2) +
                rel_scale * rel_scale);
            r->SetBinContent(b, val);
            r->SetBinError  (b, std::abs(val) * rel);
        }
    }
    return r;
}

void styleDot(TH1D* h, Color_t color, Style_t marker) {
    h->SetMarkerStyle(marker);
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

FitResult fitConst(TH1D* h, double fit_xmin, double fit_xmax,
                   bool use_weights, Color_t color) {
    const std::string fname = std::string(h->GetName()) + "_fconst";
    TF1* f = new TF1(fname.c_str(), "[0]", fit_xmin, fit_xmax);
    f->SetParName(0, "a");
    f->SetLineColor(color);
    f->SetLineWidth(2);
    f->SetLineStyle(1);
    const std::string opt = std::string("RQ") + (use_weights ? "" : "W");
    h->Fit(f, opt.c_str());
    FitResult r;
    r.a        = f->GetParameter(0);
    r.e        = f->GetParError(0);
    r.ndf      = f->GetNDF();
    r.chi2_ndf = (r.ndf > 0) ? f->GetChisquare() / r.ndf : 0.0;
    return r;
}

double safeMaxAbs(TH1D* h) {
    double mx = 0.0;
    const int nb = h->GetNbinsX();
    for (int b = 1; b <= nb; ++b) {
        const double v = std::abs(h->GetBinContent(b)) + h->GetBinError(b);
        if (std::isfinite(v) && v > mx) mx = v;
    }
    return mx;
}

}  // anonymous namespace

void mass_ee_ratio_const_fit_exp_vs_sim(double oa_cut_deg  = 0.0,
                                        bool   use_weights = true,
                                        double fit_xmin    = -1.0,
                                        double fit_xmax    = -1.0) {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetErrorX(0.0);

    // Resolve shared fit range.
    const bool   custom_range =
        (fit_xmin >= 0.0 && fit_xmin < kXmax) ||
        (fit_xmax >  0.0 && fit_xmax <= kXmax);
    const double fxlo = (fit_xmin >= 0.0 && fit_xmin < kXmax) ? fit_xmin : kXmin;
    const double fxhi = (fit_xmax >  0.0 && fit_xmax <= kXmax) ? fit_xmax : kXmax;

    // -------- Open four files --------------------------------------
    TFile* fe_em = TFile::Open("output_epem_exp.root", "READ");
    TFile* fe_pp = TFile::Open("output_epep_exp.root", "READ");
    TFile* fe_mm = TFile::Open("output_emem_exp.root", "READ");
    //--- TFile* fs    = TFile::Open("output_epem_sim.root", "READ"); SMASH LEPTON
    TFile* fs    = TFile::Open("output_epem_sim.root", "READ"); 
    //TFile* fs    = TFile::Open("output_epem_sim_std.root", "READ");
    if (!fe_em || !fe_pp || !fe_mm || !fs) {
        std::cerr << "Cannot open one of the input files\n";
        return;
    }
    TTree* nt_em  = dynamic_cast<TTree*>(fe_em->Get("dilepton_nt"));
    TTree* nt_pp  = dynamic_cast<TTree*>(fe_pp->Get("dilepton_nt"));
    TTree* nt_mm  = dynamic_cast<TTree*>(fe_mm->Get("dilepton_nt"));
    TTree* nt_sim = dynamic_cast<TTree*>(fs   ->Get("dilepton_nt"));
    if (!nt_em || !nt_pp || !nt_mm || !nt_sim) {
        std::cerr << "dilepton_nt missing in one of the files\n";
        return;
    }

    // -------- EXP cut strings (PT3 / PT2 + shared OA) --------------
    const std::string oa_exp =
        (oa_cut_deg > 0.0)
            ? std::string(" && oa>") + std::to_string(oa_cut_deg)
            : std::string{};
    const std::string cut_pt3_exp = "trigbit==8192" + oa_exp;
    const std::string cut_pt2_exp = "trigbit==4096" + oa_exp;

    // -------- SIM weight expressions (pt3/pt2 flag · weight · OA) ---
    // Macro-level autodetect (NOT per-event ternary — that would
    // over-weight stray zero-weight events in lepton mode). Lepton-mode
    // SMASH stores the dilepton-process weight in `sim_genweight`; std-
    // mode stores 0 — in that case fall back to unweighted (1.0).
    std::string sim_w;
    {
        TH1D* h = new TH1D("h_wsum_detect", "", 1, 0, 2);
        nt_sim->Draw("1>>h_wsum_detect", "sim_genweight", "goff");
        const double sum = h->Integral();
        delete h;
        sim_w = (sum > 0.0) ? "sim_genweight" : "1.0";
        std::cout << "SIM weight detection: sim_genweight integral = " << sum
                  << " → using " << sim_w << "\n";
    }
    const std::string oa_sim =
        (oa_cut_deg > 0.0)
            ? std::string(" * (oa>") + std::to_string(oa_cut_deg) + ")"
            : std::string{};
    const std::string w_pt3_sim = "(pt3==1) * " + sim_w + oa_sim;
    const std::string w_pt2_sim = "(pt2==1) * " + sim_w + oa_sim;

    std::cout << "EXP cuts:\n  PT3 = " << cut_pt3_exp
              << "\n  PT2 = " << cut_pt2_exp
              << "\nSIM weights:\n  PT3 = " << w_pt3_sim
              << "\n  PT2 = " << w_pt2_sim
              << "\nFit mode: "
              << (use_weights ? "WEIGHTED" : "UNWEIGHTED")
              << "\nFit range: [" << fxlo << ", " << fxhi << "] GeV/c^2"
              << (custom_range ? "  (custom)" : "  (full)") << "\n";

    // ============================================================================
    // EXP: build all/CB/sig for PT3 and PT2
    // ============================================================================
    TH1D* em_p3 = drawFromNt(nt_em, cut_pt3_exp, "em_p3_e");
    TH1D* pp_p3 = drawFromNt(nt_pp, cut_pt3_exp, "pp_p3_e");
    TH1D* mm_p3 = drawFromNt(nt_mm, cut_pt3_exp, "mm_p3_e");
    TH1D* em_p2 = drawFromNt(nt_em, cut_pt2_exp, "em_p2_e");
    TH1D* pp_p2 = drawFromNt(nt_pp, cut_pt2_exp, "pp_p2_e");
    TH1D* mm_p2 = drawFromNt(nt_mm, cut_pt2_exp, "mm_p2_e");
    TH1D* cb_p3  = makeCB (pp_p3, mm_p3, "cb_p3_e");
    TH1D* cb_p2  = makeCB (pp_p2, mm_p2, "cb_p2_e");
    TH1D* sig_p3 = makeSig(em_p3, cb_p3, "sig_p3_e");
    TH1D* sig_p2 = makeSig(em_p2, cb_p2, "sig_p2_e");

    TH1D* r_corr_exp = makeRatio(sig_p2, sig_p3, "r_corr_exp", 63.0);
    TH1D* r_eff_exp  = makeRatio(sig_p3, sig_p2, "r_eff_exp",  1.0 / 63.0);

    // ============================================================================
    // SIM: build pt3/pt2 weighted spectra (no CB, no factor 63)
    // ============================================================================
    TH1D* h_p3_sim = drawFromNt(nt_sim, w_pt3_sim, "h_p3_sim");
    TH1D* h_p2_sim = drawFromNt(nt_sim, w_pt2_sim, "h_p2_sim");

    TH1D* r_corr_sim = makeRatio(h_p2_sim, h_p3_sim, "r_corr_sim", 1.0);
    TH1D* r_eff_sim  = makeRatio(h_p3_sim, h_p2_sim, "r_eff_sim",  1.0);

    // ============================================================================
    // Style and titles
    // ============================================================================
    styleDot(r_corr_exp, kColExp, 20);  // filled circle
    styleDot(r_eff_exp,  kColExp, 20);
    styleDot(r_corr_sim, kColSim, 21);  // filled square
    styleDot(r_eff_sim,  kColSim, 21);

    const std::string oa_tag =
        (oa_cut_deg > 0.0)
            ? Form(", OA > %g#circ", oa_cut_deg)
            : std::string{};
    r_corr_exp->SetTitle((std::string(
        "Trigger correction factor  (exp vs sim, 20 fixed bins") + oa_tag + ");"
        "M_{e^{+}e^{-}} [GeV/c^{2}];"
        "63 #upoint N_{PT2}/N_{PT3} (exp)   N_{PT2}/N_{PT3} (sim)").c_str());
    r_eff_exp->SetTitle((std::string(
        "Trigger efficiency  (exp vs sim, 20 fixed bins") + oa_tag + ");"
        "M_{e^{+}e^{-}} [GeV/c^{2}];"
        "N_{PT3}/(63#upointN_{PT2}) (exp)   N_{PT3}/N_{PT2} (sim)").c_str());

    // ============================================================================
    // Fits — run BEFORE drawing so canvas 2 (sim/exp ratio) can reuse them.
    // TF1::Fit attaches the function to the histogram, so when we later draw
    // the histograms on canvas 1 the fit lines appear automatically.
    // ============================================================================
    std::cout << "\nConstant fit y = a (both samples, shared range):\n";
    const FitResult f_corr_exp = fitConst(r_corr_exp, fxlo, fxhi, use_weights, kColExpFit);
    const FitResult f_corr_sim = fitConst(r_corr_sim, fxlo, fxhi, use_weights, kColSimFit);
    const FitResult f_eff_exp  = fitConst(r_eff_exp,  fxlo, fxhi, use_weights, kColExpFit);
    const FitResult f_eff_sim  = fitConst(r_eff_sim,  fxlo, fxhi, use_weights, kColSimFit);
    std::cout << "  pad 1 (correction factor):\n"
              << "    exp  a = " << f_corr_exp.a << " ± " << f_corr_exp.e
              << "   chi2/ndf = " << (f_corr_exp.chi2_ndf * f_corr_exp.ndf)
              << "/" << f_corr_exp.ndf << "\n"
              << "    sim  a = " << f_corr_sim.a << " ± " << f_corr_sim.e
              << "   chi2/ndf = " << (f_corr_sim.chi2_ndf * f_corr_sim.ndf)
              << "/" << f_corr_sim.ndf << "\n"
              << "  pad 2 (efficiency):\n"
              << "    exp  a = " << f_eff_exp.a << " ± " << f_eff_exp.e
              << "   chi2/ndf = " << (f_eff_exp.chi2_ndf * f_eff_exp.ndf)
              << "/" << f_eff_exp.ndf << "\n"
              << "    sim  a = " << f_eff_sim.a << " ± " << f_eff_sim.e
              << "   chi2/ndf = " << (f_eff_sim.chi2_ndf * f_eff_sim.ndf)
              << "/" << f_eff_sim.ndf << "\n";

    // ============================================================================
    // Canvas 1: 2 panels (correction factor, efficiency) — exp + sim overlay
    // ============================================================================
    TCanvas* c = new TCanvas("c_mass_ee_ratio_exp_vs_sim",
                             "Exp vs Sim constant fit (M_{ee}, 20 fixed bins)",
                             1500, 600);
    c->Divide(2, 1, 0.001, 0.001);

    auto drawComparePad = [&](int idx, TH1D* h_exp, TH1D* h_sim,
                              const FitResult& f_exp, const FitResult& f_sim,
                              double ref_y, double y_max_floor) {
        c->cd(idx);
        const double pad_lm = 0.13, pad_rm = 0.04, pad_bm = 0.13, pad_tm = 0.10;
        gPad->SetMargin(pad_lm, pad_rm, pad_bm, pad_tm);

        // Y range covering both samples with margin.
        const double ymax = std::max(safeMaxAbs(h_exp), safeMaxAbs(h_sim));
        const double y_hi = std::max(ymax * 1.30, y_max_floor);
        h_exp->GetYaxis()->SetRangeUser(0.0, y_hi);
        h_exp->Draw("E1");          // fit line drawn with it (TF1 attached)
        h_sim->Draw("E1 SAME");

        // Reference line at ref_y (typically 1.0 = no correction / full eff).
        TLine* lref = new TLine(kXmin, ref_y, kXmax, ref_y);
        lref->SetLineStyle(3); lref->SetLineColor(kGray + 2); lref->Draw();

        // Annotation: exp top-left, sim top-right. Stack value + χ²/ndf.
        TLatex tex;
        tex.SetNDC();
        tex.SetTextAlign(11);   // left-bottom
        tex.SetTextSize(0.040);
        tex.SetTextColor(kColExpFit);
        tex.DrawLatex(0.16, 0.85,
            Form("exp:  a = %.4f #pm %.4f", f_exp.a, f_exp.e));
        tex.SetTextSize(0.030);
        tex.DrawLatex(0.16, 0.81,
            Form("       #chi^{2}/ndf = %.1f / %d",
                 f_exp.chi2_ndf * f_exp.ndf, f_exp.ndf));

        tex.SetTextAlign(31);   // right-bottom
        tex.SetTextSize(0.040);
        tex.SetTextColor(kColSimFit);
        tex.DrawLatex(0.95, 0.85,
            Form("sim:  a = %.4f #pm %.4f", f_sim.a, f_sim.e));
        tex.SetTextSize(0.030);
        tex.DrawLatex(0.95, 0.81,
            Form("#chi^{2}/ndf = %.1f / %d       ",
                 f_sim.chi2_ndf * f_sim.ndf, f_sim.ndf));

        // Marker legend (top-center, between annotations).
        TLegend* leg = new TLegend(0.40, 0.74, 0.65, 0.83);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.034);
        leg->AddEntry(h_exp, "exp",  "lpe");
        leg->AddEntry(h_sim, "sim",  "lpe");
        leg->Draw();
    };

    drawComparePad(1, r_corr_exp, r_corr_sim, f_corr_exp, f_corr_sim, 1.0, 3.5);
    drawComparePad(2, r_eff_exp,  r_eff_sim,  f_eff_exp,  f_eff_sim,  1.0, 1.2);

    // ============================================================================
    // Save with descriptive suffixes
    // ============================================================================
    gSystem->mkdir("plots/output", true);
    const std::string oa_suffix =
        (oa_cut_deg > 0.0) ? Form("_oa%g", oa_cut_deg) : std::string{};
    const std::string wt_suffix = use_weights ? std::string{} : "_unw";
    const std::string fr_suffix =
        custom_range ? Form("_fit%g-%g", fxlo, fxhi) : std::string{};
    const std::string base =
        "plots/output/mass_ee_ratio_const_fit_exp_vs_sim" +
        oa_suffix + wt_suffix + fr_suffix;
    c->SaveAs((base + ".pdf").c_str());
    c->SaveAs((base + ".png").c_str());
    std::cout << "Saved: " << base << ".{pdf,png}\n";

    // ============================================================================
    // Canvas 2: (sim · scale) / exp ratios — shape comparison after
    // rescaling sim to match exp's overall magnitude. Scale per pad:
    //   pad 1  scale_corr = a_corr_exp / a_corr_sim   (≈ 1.56 for pp45)
    //   pad 2  scale_eff  = a_eff_exp  / a_eff_sim    (≈ 0.66 for pp45)
    // After rescaling, the curve should cluster around 1 if exp and sim
    // have the SAME mass-shape (only different overall magnitude).
    // Deviations from 1.0 expose residual shape differences.
    // ============================================================================
    const double scale_corr =
        (f_corr_sim.a != 0.0) ? f_corr_exp.a / f_corr_sim.a : 1.0;
    const double scale_eff  =
        (f_eff_sim.a  != 0.0) ? f_eff_exp.a  / f_eff_sim.a  : 1.0;
    // σ(scale) = scale · √((σ_exp/a_exp)² + (σ_sim/a_sim)²)
    const double scale_corr_err = std::abs(scale_corr) * std::sqrt(
        (f_corr_exp.a != 0.0 ? std::pow(f_corr_exp.e / f_corr_exp.a, 2) : 0.0) +
        (f_corr_sim.a != 0.0 ? std::pow(f_corr_sim.e / f_corr_sim.a, 2) : 0.0));
    const double scale_eff_err  = std::abs(scale_eff)  * std::sqrt(
        (f_eff_exp.a  != 0.0 ? std::pow(f_eff_exp.e  / f_eff_exp.a,  2) : 0.0) +
        (f_eff_sim.a  != 0.0 ? std::pow(f_eff_sim.e  / f_eff_sim.a,  2) : 0.0));

    std::cout << "\nRescaling factors (a_exp / a_sim):\n"
              << "  pad 1 (correction):  scale = " << scale_corr
              << " ± " << scale_corr_err << "\n"
              << "  pad 2 (efficiency):  scale = " << scale_eff
              << " ± " << scale_eff_err  << "\n";

    TH1D* h_ratio_corr = makeScaledSimOverExp(
        r_corr_sim, r_corr_exp, scale_corr, scale_corr_err, "h_ratio_corr");
    TH1D* h_ratio_eff  = makeScaledSimOverExp(
        r_eff_sim,  r_eff_exp,  scale_eff,  scale_eff_err,  "h_ratio_eff");

    // Style: violet (between blue exp and red sim) to mark these as
    // "exp/sim comparison" rather than pure-exp or pure-sim curves.
    const Color_t kColRatio = kViolet + 1;
    styleDot(h_ratio_corr, kColRatio, 20);
    styleDot(h_ratio_eff,  kColRatio, 20);

    const std::string oa_tag_c2 =
        (oa_cut_deg > 0.0) ? Form(", OA > %g#circ", oa_cut_deg) : std::string{};
    h_ratio_corr->SetTitle((std::string(
        "(sim #upoint scale) / exp  —  correction factor (20 fixed bins") +
        oa_tag_c2 + ");"
        "M_{e^{+}e^{-}} [GeV/c^{2}];"
        "(N_{PT2}/N_{PT3})_{sim} #upoint scale / (63 N_{PT2}/N_{PT3})_{exp}").c_str());
    h_ratio_eff->SetTitle((std::string(
        "(sim #upoint scale) / exp  —  efficiency (20 fixed bins") +
        oa_tag_c2 + ");"
        "M_{e^{+}e^{-}} [GeV/c^{2}];"
        "(N_{PT3}/N_{PT2})_{sim} #upoint scale / (N_{PT3}/(63N_{PT2}))_{exp}").c_str());

    TCanvas* c2 = new TCanvas("c_mass_ee_ratio_exp_vs_sim_ratio",
                              "(sim rescaled) / exp shape ratio",
                              1500, 600);
    c2->Divide(2, 1, 0.001, 0.001);

    auto drawRatioPad = [&](int idx, TH1D* h, double scale, double scale_err,
                            const std::string& scale_label) {
        c2->cd(idx);
        const double pad_lm = 0.13, pad_rm = 0.04, pad_bm = 0.13, pad_tm = 0.10;
        gPad->SetMargin(pad_lm, pad_rm, pad_bm, pad_tm);

        // Y range centred on 1: cover bulk of data with margin.
        double ymax = 0.0, ymin = std::numeric_limits<double>::infinity();
        for (int b = 1; b <= kNb; ++b) {
            const double v = h->GetBinContent(b);
            const double e = h->GetBinError(b);
            if (!std::isfinite(v) || v == 0.0) continue;
            // robust: ignore points with rel err > 50% from range calc
            if (e > 0.0 && std::abs(e / v) > 0.5) continue;
            ymax = std::max(ymax, v + e);
            ymin = std::min(ymin, std::max(0.0, v - e));
        }
        if (!std::isfinite(ymin)) ymin = 0.0;
        if (ymax <= 0.0)          ymax = 2.0;
        // Symmetric padding around 1 if data clusters there
        const double y_lo = std::max(0.0, std::min(ymin, 0.5) - 0.15);
        const double y_hi = std::max(ymax * 1.20, 1.5);
        h->GetYaxis()->SetRangeUser(y_lo, y_hi);

        h->Draw("E1");

        // Reference line at 1.0 (perfect agreement after rescaling).
        TLine* lref = new TLine(kXmin, 1.0, kXmax, 1.0);
        lref->SetLineStyle(2); lref->SetLineColor(kGray + 3);
        lref->SetLineWidth(2);
        lref->Draw();

        // Annotation: scale factor used (top-left).
        TLatex tex;
        tex.SetNDC();
        tex.SetTextAlign(11);
        tex.SetTextSize(0.040);
        tex.SetTextColor(kColRatio);
        tex.DrawLatex(0.16, 0.85,
            Form("scale = a_{exp} / a_{sim} = %.4f #pm %.4f",
                 scale, scale_err));
        tex.SetTextSize(0.033);
        tex.SetTextColor(kGray + 3);
        tex.DrawLatex(0.16, 0.80,
            Form("(%s)", scale_label.c_str()));
        tex.DrawLatex(0.16, 0.76,
            "= 1 #Rightarrow sim and exp have identical mass shape");
    };

    drawRatioPad(1, h_ratio_corr, scale_corr, scale_corr_err,
                 Form("= %.4f / %.4f", f_corr_exp.a, f_corr_sim.a));
    drawRatioPad(2, h_ratio_eff,  scale_eff,  scale_eff_err,
                 Form("= %.4f / %.4f", f_eff_exp.a, f_eff_sim.a));

    const std::string base2 =
        "plots/output/mass_ee_ratio_const_fit_exp_vs_sim_ratio" +
        oa_suffix + wt_suffix + fr_suffix;
    c2->SaveAs((base2 + ".pdf").c_str());
    c2->SaveAs((base2 + ".png").c_str());
    std::cout << "Saved: " << base2 << ".{pdf,png}\n";

    fe_em->Close(); fe_pp->Close(); fe_mm->Close(); fs->Close();
}
