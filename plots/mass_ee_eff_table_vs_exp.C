// mass_ee_eff_table_vs_exp.C — overlay of RAYANE trigger-efficiency data
// (provided externally as a HARDCODED table of (M, y, err_y) tuples) with
// experimental data (CB-subtracted), plus a (exp·scale)/RAYANE shape-ratio
// panel.
//
// Two panels:
//
//   Pad 1 — Trigger efficiency overlay:
//     RAYANE: 13 fixed bins on [0, 1.04] GeV/c² with values + per-bin
//             errors taken directly from kRayaneTable below (PT3/PT2 per
//             bin from external analysis).
//     EXP:    SAME 13 bins on [0, 1.04] from output_{epem,epep,emem}_exp
//             .root with CB subtraction (signal = N_ep_em - 2·√(N_++·N_--)),
//             computed as N_PT3 / (63·N_PT2) (factor 63 = PT2 downscale).
//     Both fitted with constant y = a (separate fits), values annotated.
//
//   Pad 2 — (exp · scale) / RAYANE shape ratio:
//     scale = a_RAYANE / a_EXP from the two fits in pad 1.
//     A third constant fit is overlaid on the ratio; values near 1.0 mean
//     EXP and RAYANE have the same mass shape, deviations expose residual
//     differences after normalisation.
//
// Bin layout: 13 bins of width 80 MeV, edges 0.00, 0.08, ..., 1.04.
// Bin i (1-based) center = 0.04 + (i-1)·0.08, matching kRayaneTable order.
//
// Arguments (all optional, semantics identical to the other macros):
//   1) oa_cut_deg   — OA threshold applied to EXP only (table is fixed);
//                     default 0 ⇒ no cut
//   2) use_weights  — true → weighted χ² fit; false → unweighted; default true
//   3) fit_xmin     — fit range lower bound [GeV/c²]; -1 ⇒ full
//   4) fit_xmax     — fit range upper bound [GeV/c²]; -1 ⇒ full
//
// Usage:
//   root -l -b -q plots/mass_ee_eff_table_vs_exp.C
//   root -l -b -q 'plots/mass_ee_eff_table_vs_exp.C(4)'
//   root -l -b -q 'plots/mass_ee_eff_table_vs_exp.C(0, true, 0.1, 0.8)'
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
#include <vector>

namespace {

// -------- Binning (table-driven) --------
constexpr int    kNb   = 13;
constexpr double kXmin = 0.0;
constexpr double kXmax = 1.04;     // 13 bins × 0.08 GeV
constexpr double kBinW = 0.08;

// -------- RAYANE table — bin centers, PT3/PT2 values and errors --------
struct RayanePoint { double m; double y; double ey; };
const std::vector<RayanePoint> kRayaneTable = {
    {0.04, 0.532, 0.025},
    {0.12, 0.620, 0.022},
    {0.20, 0.627, 0.024},
    {0.28, 0.664, 0.030},
    {0.36, 0.613, 0.032},
    {0.44, 0.607, 0.038},
    {0.52, 0.532, 0.043},
    {0.60, 0.636, 0.065},
    {0.68, 0.592, 0.060},
    {0.76, 0.558, 0.050},
    {0.84, 0.590, 0.295},
    {0.92, 0.085, 0.085},
    {1.00, 0.535, 0.310},
};

// -------- Colors --------
constexpr Color_t kColRay    = kGreen  + 2;
constexpr Color_t kColRayFit = kGreen  + 3;
constexpr Color_t kColExp    = kBlue   + 1;
constexpr Color_t kColExpFit = kBlue   + 2;
constexpr Color_t kColRatio  = kViolet + 1;

// -------- Helpers --------

TH1D* buildRayaneTable(const std::string& name) {
    TH1D* h = new TH1D(name.c_str(), "", kNb, kXmin, kXmax);
    h->SetDirectory(nullptr);
    h->Sumw2();
    for (size_t i = 0; i < kRayaneTable.size(); ++i) {
        const int b = static_cast<int>(i) + 1;
        if (b > kNb) break;
        h->SetBinContent(b, kRayaneTable[i].y);
        h->SetBinError  (b, kRayaneTable[i].ey);
    }
    return h;
}

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

/// Scaled-numerator / denominator ratio with full error propagation.
/// Scale uncertainty (from the two constant fits) propagates through too —
/// added in quadrature with per-bin relative errors of num and den.
TH1D* makeScaledNumOverDen(TH1D* num, TH1D* den, double scale,
                           double scale_err, const std::string& name) {
    TH1D* r = static_cast<TH1D*>(num->Clone(name.c_str()));
    r->SetDirectory(nullptr);
    r->Reset();
    const int nb = num->GetNbinsX();
    const double rel_scale =
        (scale != 0.0) ? std::abs(scale_err / scale) : 0.0;
    for (int b = 1; b <= nb; ++b) {
        const double s  = num->GetBinContent(b);
        const double es = num->GetBinError(b);
        const double e  = den->GetBinContent(b);
        const double ee = den->GetBinError(b);
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

void mass_ee_eff_table_vs_exp(double oa_cut_deg  = 0.0,
                              bool   use_weights = true,
                              double fit_xmin    = -1.0,
                              double fit_xmax    = -1.0) {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetErrorX(0.0);

    // Resolve fit range: sentinel −1 ⇒ use the histogram's full range.
    const bool   custom_range =
        (fit_xmin >= 0.0 && fit_xmin < kXmax) ||
        (fit_xmax >  0.0 && fit_xmax <= kXmax);
    const double fxlo = (fit_xmin >= 0.0 && fit_xmin < kXmax) ? fit_xmin : kXmin;
    const double fxhi = (fit_xmax >  0.0 && fit_xmax <= kXmax) ? fit_xmax : kXmax;

    // -------- RAYANE from table --------
    TH1D* h_ray = buildRayaneTable("h_rayane_eff_table");

    // -------- EXP (CB-subtracted) from three files --------
    TFile* fe_em = TFile::Open("output_epem_exp.root", "READ");
    TFile* fe_pp = TFile::Open("output_epep_exp.root", "READ");
    TFile* fe_mm = TFile::Open("output_emem_exp.root", "READ");
    if (!fe_em || !fe_pp || !fe_mm ||
        fe_em->IsZombie() || fe_pp->IsZombie() || fe_mm->IsZombie()) {
        std::cerr << "Cannot open one of the EXP input files\n";
        return;
    }
    TTree* nt_em = dynamic_cast<TTree*>(fe_em->Get("dilepton_nt"));
    TTree* nt_pp = dynamic_cast<TTree*>(fe_pp->Get("dilepton_nt"));
    TTree* nt_mm = dynamic_cast<TTree*>(fe_mm->Get("dilepton_nt"));
    if (!nt_em || !nt_pp || !nt_mm) {
        std::cerr << "dilepton_nt missing in one of the EXP files\n";
        return;
    }

    const std::string oa_exp =
        (oa_cut_deg > 0.0)
            ? std::string(" && oa>") + std::to_string(oa_cut_deg)
            : std::string{};
    const std::string cut_pt3_exp = "trigbit==8192" + oa_exp;
    const std::string cut_pt2_exp = "trigbit==4096" + oa_exp;

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

    // PT3 / (63 · PT2)  →  divide PT2 by 63 effectively via scale = 1/63.
    TH1D* h_exp = makeRatio(sig_p3, sig_p2, "h_exp_eff_table", 1.0 / 63.0);

    std::cout << "Fit mode: "
              << (use_weights ? "WEIGHTED" : "UNWEIGHTED")
              << "\nFit range: [" << fxlo << ", " << fxhi << "] GeV/c^2"
              << (custom_range ? "  (custom)" : "  (full)") << "\n";

    // -------- Styling and titles --------
    styleDot(h_ray, kColRay, 20);
    styleDot(h_exp, kColExp, 21);

    const std::string oa_tag =
        (oa_cut_deg > 0.0)
            ? Form(", EXP OA > %g#circ", oa_cut_deg)
            : std::string{};
    h_ray->SetTitle((std::string(
        "Trigger efficiency  RAYANE(table) vs exp — 13 fixed bins") +
        oa_tag + ");"
        "M_{e^{+}e^{-}} [GeV/c^{2}];"
        "PT3/PT2 (RAYANE)   PT3/(63 #upoint PT2) (exp)").c_str());

    // -------- Fits (BEFORE drawing so canvas reuses results) --------
    const FitResult f_ray = fitConst(h_ray, fxlo, fxhi, use_weights, kColRayFit);
    const FitResult f_exp = fitConst(h_exp, fxlo, fxhi, use_weights, kColExpFit);
    std::cout << "\nPad 1 fits (constant y = a):\n"
              << "  RAYANE(table)  a = " << f_ray.a << " ± " << f_ray.e
              << "   chi2/ndf = " << (f_ray.chi2_ndf * f_ray.ndf)
              << "/" << f_ray.ndf << "\n"
              << "  exp            a = " << f_exp.a << " ± " << f_exp.e
              << "   chi2/ndf = " << (f_exp.chi2_ndf * f_exp.ndf)
              << "/" << f_exp.ndf << "\n";

    // -------- Scale + (exp · scale) / RAYANE ratio --------
    const double scale =
        (f_exp.a != 0.0) ? f_ray.a / f_exp.a : 1.0;
    const double scale_err = std::abs(scale) * std::sqrt(
        (f_ray.a != 0.0 ? std::pow(f_ray.e / f_ray.a, 2) : 0.0) +
        (f_exp.a != 0.0 ? std::pow(f_exp.e / f_exp.a, 2) : 0.0));

    TH1D* h_ratio = makeScaledNumOverDen(
        h_exp, h_ray, scale, scale_err, "h_ratio_eff_table");
    styleDot(h_ratio, kColRatio, 20);
    h_ratio->SetTitle((std::string(
        "(exp #upoint scale) / RAYANE — shape ratio (13 fixed bins") + oa_tag + ");"
        "M_{e^{+}e^{-}} [GeV/c^{2}];"
        "(PT3/(63#upointPT2))_{exp} #upoint scale / (PT3/PT2)_{RAYANE}").c_str());

    const FitResult f_ratio =
        fitConst(h_ratio, fxlo, fxhi, use_weights, kColRatio);
    std::cout << "\nPad 2 fit (constant y = a on the ratio):\n"
              << "  scale = a_RAYANE / a_exp = " << scale << " ± " << scale_err << "\n"
              << "  ratio fit a = " << f_ratio.a << " ± " << f_ratio.e
              << "   chi2/ndf = " << (f_ratio.chi2_ndf * f_ratio.ndf)
              << "/" << f_ratio.ndf << "\n";

    // -------- Canvas: 2 panels --------
    TCanvas* c = new TCanvas(
        "c_mass_ee_eff_table_vs_exp",
        "Trigger eff (RAYANE table vs exp) + (exp #upoint scale) / RAYANE ratio",
        1500, 600);
    c->Divide(2, 1, 0.001, 0.001);

    // Pad 1 — Efficiency overlay
    c->cd(1);
    {
        const double pad_lm = 0.13, pad_rm = 0.04, pad_bm = 0.13, pad_tm = 0.10;
        gPad->SetMargin(pad_lm, pad_rm, pad_bm, pad_tm);

        const double ymax = std::max(safeMaxAbs(h_ray), safeMaxAbs(h_exp));
        const double y_hi = std::max(ymax * 1.30, 1.2);
        h_ray->GetYaxis()->SetRangeUser(0.0, y_hi);
        h_ray->Draw("E1");          // fit attached
        h_exp->Draw("E1 SAME");

        TLine* l1 = new TLine(kXmin, 1.0, kXmax, 1.0);
        l1->SetLineStyle(3); l1->SetLineColor(kGray + 2); l1->Draw();

        TLatex tex;
        tex.SetNDC();
        tex.SetTextAlign(11);
        tex.SetTextSize(0.040);
        tex.SetTextColor(kColRayFit);
        tex.DrawLatex(0.16, 0.85,
            Form("RAYANE:  a = %.4f #pm %.4f", f_ray.a, f_ray.e));
        tex.SetTextSize(0.030);
        tex.DrawLatex(0.16, 0.81,
            Form("           #chi^{2}/ndf = %.1f / %d",
                 f_ray.chi2_ndf * f_ray.ndf, f_ray.ndf));

        tex.SetTextAlign(31);
        tex.SetTextSize(0.040);
        tex.SetTextColor(kColExpFit);
        tex.DrawLatex(0.95, 0.85,
            Form("exp:  a = %.4f #pm %.4f", f_exp.a, f_exp.e));
        tex.SetTextSize(0.030);
        tex.DrawLatex(0.95, 0.81,
            Form("#chi^{2}/ndf = %.1f / %d       ",
                 f_exp.chi2_ndf * f_exp.ndf, f_exp.ndf));

        TLegend* leg = new TLegend(0.40, 0.74, 0.65, 0.83);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.034);
        leg->AddEntry(h_ray, "RAYANE (table)", "lpe");
        leg->AddEntry(h_exp, "exp",            "lpe");
        leg->Draw();
    }

    // Pad 2 — (exp · scale) / RAYANE ratio with a constant fit
    c->cd(2);
    {
        const double pad_lm = 0.13, pad_rm = 0.04, pad_bm = 0.13, pad_tm = 0.10;
        gPad->SetMargin(pad_lm, pad_rm, pad_bm, pad_tm);

        // Y range: centred on 1, robust ignoring big-error bins.
        double ymax = 0.0, ymin = std::numeric_limits<double>::infinity();
        for (int b = 1; b <= kNb; ++b) {
            const double v = h_ratio->GetBinContent(b);
            const double e = h_ratio->GetBinError(b);
            if (!std::isfinite(v) || v == 0.0) continue;
            if (e > 0.0 && std::abs(e / v) > 0.5) continue;
            ymax = std::max(ymax, v + e);
            ymin = std::min(ymin, std::max(0.0, v - e));
        }
        if (!std::isfinite(ymin)) ymin = 0.0;
        if (ymax <= 0.0)          ymax = 2.0;
        const double y_lo = std::max(0.0, std::min(ymin, 0.5) - 0.15);
        const double y_hi = std::max(ymax * 1.20, 1.5);
        h_ratio->GetYaxis()->SetRangeUser(y_lo, y_hi);

        h_ratio->Draw("E1");

        // Reference line at 1.0 (perfect agreement after rescaling).
        TLine* l1 = new TLine(kXmin, 1.0, kXmax, 1.0);
        l1->SetLineStyle(2); l1->SetLineColor(kGray + 3); l1->SetLineWidth(2);
        l1->Draw();

        TLatex tex;
        tex.SetNDC();
        tex.SetTextAlign(11);
        tex.SetTextSize(0.040);
        tex.SetTextColor(kColRatio);
        tex.DrawLatex(0.16, 0.85,
            Form("scale = a_{RAYANE}/a_{exp} = %.4f #pm %.4f",
                 scale, scale_err));
        tex.SetTextSize(0.035);
        tex.DrawLatex(0.16, 0.80,
            Form("ratio fit  a = %.4f #pm %.4f", f_ratio.a, f_ratio.e));
        tex.SetTextSize(0.030);
        tex.SetTextColor(kGray + 3);
        tex.DrawLatex(0.16, 0.76,
            Form("#chi^{2}/ndf = %.1f / %d",
                 f_ratio.chi2_ndf * f_ratio.ndf, f_ratio.ndf));
        tex.DrawLatex(0.16, 0.72,
            "fit = 1 #Rightarrow exp and RAYANE identical mass shape");
    }

    // -------- Save --------
    gSystem->mkdir("plots/output", true);
    const std::string oa_suffix =
        (oa_cut_deg > 0.0) ? Form("_oa%g", oa_cut_deg) : std::string{};
    const std::string wt_suffix = use_weights ? std::string{} : "_unw";
    const std::string fr_suffix =
        custom_range ? Form("_fit%g-%g", fxlo, fxhi) : std::string{};
    const std::string base =
        "plots/output/mass_ee_eff_table_vs_exp" +
        oa_suffix + wt_suffix + fr_suffix;
    c->SaveAs((base + ".pdf").c_str());
    c->SaveAs((base + ".png").c_str());
    std::cout << "Saved: " << base << ".{pdf,png}\n";

    fe_em->Close();
    fe_pp->Close();
    fe_mm->Close();
}
