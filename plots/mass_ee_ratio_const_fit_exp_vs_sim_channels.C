// mass_ee_ratio_const_fit_exp_vs_sim_channels.C — like
// mass_ee_ratio_const_fit_exp_vs_sim.C but with three extra SIM curves on each
// panel, restricted to single physics channels via the new `ep_sim_geninfo1`
// branch added to dilepton_nt:
//   • sim_geninfo1 == 7051  → π⁰ Dalitz (π⁰ → e⁺e⁻γ), clean (m_ee < m_π⁰)
//   • sim_geninfo1 == 17051 → η  Dalitz (η  → e⁺e⁻γ), clean (m_ee < m_η)
//   • sim_geninfo1 == 41    → ω all-decays MIX (Dalitz + direct + conv.)
//                             ω → π⁰ e⁺e⁻ Dalitz (BR≈7.7e-4, m_ee≤0.65 GeV) +
//                             ω → e⁺e⁻ direct (BR≈7.4e-5, m_ee≈0.78 GeV) are
//                             both encoded as the same geninfo (=410203) in
//                             SMASH — cannot be separated by geninfo alone.
//
// Two panels (20 fixed bins on [0, 1.4] GeV/c²):
//
//   Pad 1 — Trigger correction factor:
//     exp:        63 · N_PT2 / N_PT3                 (CB-subtracted)
//     sim (all):       N_PT2 / N_PT3                 (no downscale)
//     sim (π⁰):        N_PT2 / N_PT3   filtered to geninfo1 == 7051
//     sim (η ):        N_PT2 / N_PT3   filtered to geninfo1 == 17051
//
//   Pad 2 — Trigger efficiency:
//     exp:        N_PT3 / (63 · N_PT2)
//     sim (all):  N_PT3 / N_PT2
//     sim (π⁰):   N_PT3 / N_PT2 filtered
//     sim (η ):   N_PT3 / N_PT2 filtered
//
// Each curve gets its own constant fit y = a, drawn in matching color.
//
// Arguments (identical semantics to the reference macro):
//   1) oa_cut_deg   — OA threshold; default 0 ⇒ no cut
//   2) use_weights  — true → weighted χ² fit; false → unweighted; default true
//   3) fit_xmin     — fit range lower bound [GeV/c²]; -1 ⇒ full
//   4) fit_xmax     — fit range upper bound [GeV/c²]; -1 ⇒ full
//
// Usage:
//   root -l -b -q plots/mass_ee_ratio_const_fit_exp_vs_sim_channels.C
//   root -l -b -q 'plots/mass_ee_ratio_const_fit_exp_vs_sim_channels.C(0,true,0.1,0.8)'
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
#include <vector>

namespace {

constexpr int    kNb   = 20;
constexpr double kXmin = 0.0;
constexpr double kXmax = 1.4;

constexpr Color_t kColExp     = kBlue    + 1;
constexpr Color_t kColExpFit  = kBlue    + 2;
constexpr Color_t kColSim     = kRed     + 1;
constexpr Color_t kColSimFit  = kRed     + 2;
constexpr Color_t kColPi0     = kGreen   + 2;
constexpr Color_t kColPi0Fit  = kGreen   + 3;
constexpr Color_t kColEta     = kMagenta + 1;
constexpr Color_t kColEtaFit  = kMagenta + 2;
constexpr Color_t kColOmega   = kBlack;
constexpr Color_t kColOmegaFit= kBlack;

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

void styleDot(TH1D* h, Color_t color, Style_t marker, double size = 0.9) {
    h->SetMarkerStyle(marker);
    h->SetMarkerSize(size);
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

std::string detectSimWeight(TTree* nt) {
    TH1D* h = new TH1D("h_wsum_detect_chs", "", 1, 0, 2);
    nt->Draw("1>>h_wsum_detect_chs", "sim_genweight", "goff");
    const double sum = h->Integral();
    delete h;
    if (sum > 0.0) {
        std::cout << "SIM weight: sim_genweight integral = " << sum
                  << "  → using sim_genweight\n";
        return "sim_genweight";
    }
    std::cout << "SIM weight: sim_genweight integral = 0  → unweighted (1.0)\n";
    return "1.0";
}

}  // anonymous namespace

void mass_ee_ratio_const_fit_exp_vs_sim_channels(double oa_cut_deg  = 0.0,
                                                 bool   use_weights = true,
                                                 double fit_xmin    = -1.0,
                                                 double fit_xmax    = -1.0) {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetErrorX(0.0);

    const bool   custom_range =
        (fit_xmin >= 0.0 && fit_xmin < kXmax) ||
        (fit_xmax >  0.0 && fit_xmax <= kXmax);
    const double fxlo = (fit_xmin >= 0.0 && fit_xmin < kXmax) ? fit_xmin : kXmin;
    const double fxhi = (fit_xmax >  0.0 && fit_xmax <= kXmax) ? fit_xmax : kXmax;

    // -------- Open files --------
    TFile* fe_em = TFile::Open("output_epem_exp.root", "READ");
    TFile* fe_pp = TFile::Open("output_epep_exp.root", "READ");
    TFile* fe_mm = TFile::Open("output_emem_exp.root", "READ");
    TFile* fs    = TFile::Open("output_epem_sim.root", "READ");
    if (!fe_em || !fe_pp || !fe_mm || !fs ||
        fe_em->IsZombie() || fe_pp->IsZombie() ||
        fe_mm->IsZombie() || fs->IsZombie()) {
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

    // -------- EXP cuts (PT3 / PT2 + shared OA) --------
    const std::string oa_exp =
        (oa_cut_deg > 0.0)
            ? std::string(" && oa>") + std::to_string(oa_cut_deg)
            : std::string{};
    const std::string cut_pt3_exp = "trigbit==8192" + oa_exp;
    const std::string cut_pt2_exp = "trigbit==4096" + oa_exp;

    // -------- SIM weight (autodetect) + cuts ----------
    const std::string sim_w = detectSimWeight(nt_sim);
    const std::string oa_sim =
        (oa_cut_deg > 0.0)
            ? std::string(" * (oa>") + std::to_string(oa_cut_deg) + ")"
            : std::string{};
    const std::string w_pt3_sim = "(pt3==1) * " + sim_w + oa_sim;
    const std::string w_pt2_sim = "(pt2==1) * " + sim_w + oa_sim;

    // Channel-filtered SIM weights — multiply in (geninfo1 == X)
    const std::string w_pt3_pi0 =
        "(pt3==1) * (ep_sim_geninfo1==7051) * " + sim_w + oa_sim;
    const std::string w_pt2_pi0 =
        "(pt2==1) * (ep_sim_geninfo1==7051) * " + sim_w + oa_sim;
    const std::string w_pt3_eta =
        "(pt3==1) * (ep_sim_geninfo1==17051) * " + sim_w + oa_sim;
    const std::string w_pt2_eta =
        "(pt2==1) * (ep_sim_geninfo1==17051) * " + sim_w + oa_sim;
    const std::string w_pt3_omg =
        "(pt3==1) * (ep_sim_geninfo1==41) * " + sim_w + oa_sim;
    const std::string w_pt2_omg =
        "(pt2==1) * (ep_sim_geninfo1==41) * " + sim_w + oa_sim;

    std::cout << "Fit mode: " << (use_weights ? "WEIGHTED" : "UNWEIGHTED")
              << "   Fit range: [" << fxlo << ", " << fxhi << "] GeV/c^2"
              << (custom_range ? "  (custom)" : "  (full)") << "\n\n";

    // ============================================================================
    // EXP: build all/CB/sig for PT3 and PT2
    // ============================================================================
    TH1D* em_p3 = drawFromNt(nt_em, cut_pt3_exp, "em_p3_chs");
    TH1D* pp_p3 = drawFromNt(nt_pp, cut_pt3_exp, "pp_p3_chs");
    TH1D* mm_p3 = drawFromNt(nt_mm, cut_pt3_exp, "mm_p3_chs");
    TH1D* em_p2 = drawFromNt(nt_em, cut_pt2_exp, "em_p2_chs");
    TH1D* pp_p2 = drawFromNt(nt_pp, cut_pt2_exp, "pp_p2_chs");
    TH1D* mm_p2 = drawFromNt(nt_mm, cut_pt2_exp, "mm_p2_chs");
    TH1D* cb_p3  = makeCB (pp_p3, mm_p3, "cb_p3_chs");
    TH1D* cb_p2  = makeCB (pp_p2, mm_p2, "cb_p2_chs");
    TH1D* sig_p3 = makeSig(em_p3, cb_p3, "sig_p3_chs");
    TH1D* sig_p2 = makeSig(em_p2, cb_p2, "sig_p2_chs");

    TH1D* r_corr_exp = makeRatio(sig_p2, sig_p3, "r_corr_exp_chs", 63.0);
    TH1D* r_eff_exp  = makeRatio(sig_p3, sig_p2, "r_eff_exp_chs",  1.0 / 63.0);

    // ============================================================================
    // SIM total + channel-filtered spectra (no CB, no factor 63)
    // ============================================================================
    TH1D* h_p3_sim_all = drawFromNt(nt_sim, w_pt3_sim, "h_p3_sim_all_chs");
    TH1D* h_p2_sim_all = drawFromNt(nt_sim, w_pt2_sim, "h_p2_sim_all_chs");
    TH1D* h_p3_sim_pi0 = drawFromNt(nt_sim, w_pt3_pi0, "h_p3_sim_pi0_chs");
    TH1D* h_p2_sim_pi0 = drawFromNt(nt_sim, w_pt2_pi0, "h_p2_sim_pi0_chs");
    TH1D* h_p3_sim_eta = drawFromNt(nt_sim, w_pt3_eta, "h_p3_sim_eta_chs");
    TH1D* h_p2_sim_eta = drawFromNt(nt_sim, w_pt2_eta, "h_p2_sim_eta_chs");
    TH1D* h_p3_sim_omg = drawFromNt(nt_sim, w_pt3_omg, "h_p3_sim_omg_chs");
    TH1D* h_p2_sim_omg = drawFromNt(nt_sim, w_pt2_omg, "h_p2_sim_omg_chs");

    TH1D* r_corr_sim_all = makeRatio(h_p2_sim_all, h_p3_sim_all, "r_corr_sim_all_chs", 1.0);
    TH1D* r_eff_sim_all  = makeRatio(h_p3_sim_all, h_p2_sim_all, "r_eff_sim_all_chs",  1.0);
    TH1D* r_corr_sim_pi0 = makeRatio(h_p2_sim_pi0, h_p3_sim_pi0, "r_corr_sim_pi0_chs", 1.0);
    TH1D* r_eff_sim_pi0  = makeRatio(h_p3_sim_pi0, h_p2_sim_pi0, "r_eff_sim_pi0_chs",  1.0);
    TH1D* r_corr_sim_eta = makeRatio(h_p2_sim_eta, h_p3_sim_eta, "r_corr_sim_eta_chs", 1.0);
    TH1D* r_eff_sim_eta  = makeRatio(h_p3_sim_eta, h_p2_sim_eta, "r_eff_sim_eta_chs",  1.0);
    TH1D* r_corr_sim_omg = makeRatio(h_p2_sim_omg, h_p3_sim_omg, "r_corr_sim_omg_chs", 1.0);
    TH1D* r_eff_sim_omg  = makeRatio(h_p3_sim_omg, h_p2_sim_omg, "r_eff_sim_omg_chs",  1.0);

    // ============================================================================
    // Style and titles
    // ============================================================================
    // Markers: exp/sim-total keep filled circle/square (default size 0.9);
    // channel-filtered curves use open shapes per user spec.
    styleDot(r_corr_exp,     kColExp,   20);            // filled circle
    styleDot(r_eff_exp,      kColExp,   20);
    styleDot(r_corr_sim_all, kColSim,   21);            // filled square
    styleDot(r_eff_sim_all,  kColSim,   21);
    styleDot(r_corr_sim_pi0, kColPi0,   25, 1.1);       // open square, green, 1.1
    styleDot(r_eff_sim_pi0,  kColPi0,   25, 1.1);
    styleDot(r_corr_sim_eta, kColEta,   24, 1.1);       // open circle, magenta, 1.1
    styleDot(r_eff_sim_eta,  kColEta,   24, 1.1);
    styleDot(r_corr_sim_omg, kColOmega, 26, 1.3);       // open triangle-up, black, 1.3
    styleDot(r_eff_sim_omg,  kColOmega, 26, 1.3);

    const std::string oa_tag =
        (oa_cut_deg > 0.0)
            ? Form(", OA > %g#circ", oa_cut_deg)
            : std::string{};
    r_corr_exp->SetTitle((std::string(
        "Trigger correction factor  (exp + sim total / #pi^{0} Dal. / #eta Dal. / #omega all, 20 fixed bins") + oa_tag + ");"
        "M_{e^{+}e^{-}} [GeV/c^{2}];"
        "63 #upoint N_{PT2}/N_{PT3} (exp)   N_{PT2}/N_{PT3} (sim)").c_str());
    r_eff_exp->SetTitle((std::string(
        "Trigger efficiency  (exp + sim total / #pi^{0} Dal. / #eta Dal. / #omega all, 20 fixed bins") + oa_tag + ");"
        "M_{e^{+}e^{-}} [GeV/c^{2}];"
        "N_{PT3}/(63#upointN_{PT2}) (exp)   N_{PT3}/N_{PT2} (sim)").c_str());

    // ============================================================================
    // Fits (drawn-with-data; ran BEFORE the canvas)
    // ============================================================================
    std::cout << "Constant fit y = a:\n";
    const FitResult f_corr_exp = fitConst(r_corr_exp,     fxlo, fxhi, use_weights, kColExpFit);
    const FitResult f_corr_all = fitConst(r_corr_sim_all, fxlo, fxhi, use_weights, kColSimFit);
    const FitResult f_corr_pi0 = fitConst(r_corr_sim_pi0, fxlo, fxhi, use_weights, kColPi0Fit);
    const FitResult f_corr_eta = fitConst(r_corr_sim_eta, fxlo, fxhi, use_weights, kColEtaFit);
    const FitResult f_corr_omg = fitConst(r_corr_sim_omg, fxlo, fxhi, use_weights, kColOmegaFit);
    const FitResult f_eff_exp  = fitConst(r_eff_exp,      fxlo, fxhi, use_weights, kColExpFit);
    const FitResult f_eff_all  = fitConst(r_eff_sim_all,  fxlo, fxhi, use_weights, kColSimFit);
    const FitResult f_eff_pi0  = fitConst(r_eff_sim_pi0,  fxlo, fxhi, use_weights, kColPi0Fit);
    const FitResult f_eff_eta  = fitConst(r_eff_sim_eta,  fxlo, fxhi, use_weights, kColEtaFit);
    const FitResult f_eff_omg  = fitConst(r_eff_sim_omg,  fxlo, fxhi, use_weights, kColOmegaFit);

    auto pf = [](const char* tag, const FitResult& f) {
        printf("  %-22s a = %.4f #pm %.4f   chi2/ndf = %8.1f / %d\n",
               tag, f.a, f.e, f.chi2_ndf * f.ndf, f.ndf);
    };
    std::cout << " Pad 1 (correction factor):\n";
    pf("exp",              f_corr_exp);
    pf("sim total",        f_corr_all);
    pf("sim #pi^0 Dalitz", f_corr_pi0);
    pf("sim #eta Dalitz",  f_corr_eta);
    pf("sim #omega all (41)",  f_corr_omg);
    std::cout << " Pad 2 (efficiency):\n";
    pf("exp",              f_eff_exp);
    pf("sim total",        f_eff_all);
    pf("sim #pi^0 Dalitz", f_eff_pi0);
    pf("sim #eta Dalitz",  f_eff_eta);
    pf("sim #omega all (41)",  f_eff_omg);

    // ============================================================================
    // Canvas
    // ============================================================================
    TCanvas* c = new TCanvas("c_mass_ee_ratio_exp_vs_sim_channels",
                             "Exp vs Sim (total + pi0/eta Dalitz) constant fit",
                             1600, 700);
    c->Divide(2, 1, 0.001, 0.001);

    auto drawComparePad = [&](int idx,
                              TH1D* h_exp, TH1D* h_all,
                              TH1D* h_pi0, TH1D* h_eta, TH1D* h_omg,
                              const FitResult& f_exp, const FitResult& f_all,
                              const FitResult& f_pi0, const FitResult& f_eta,
                              const FitResult& f_omg,
                              double ref_y, double y_max_floor) {
        c->cd(idx);
        const double pad_lm = 0.13, pad_rm = 0.04, pad_bm = 0.13, pad_tm = 0.10;
        gPad->SetMargin(pad_lm, pad_rm, pad_bm, pad_tm);

        const double ymax = std::max({safeMaxAbs(h_exp), safeMaxAbs(h_all),
                                      safeMaxAbs(h_pi0), safeMaxAbs(h_eta),
                                      safeMaxAbs(h_omg)});
        const double y_hi = std::max(ymax * 1.35, y_max_floor);
        h_exp->GetYaxis()->SetRangeUser(0.0, y_hi);
        h_exp->Draw("E1");
        h_all->Draw("E1 SAME");
        h_pi0->Draw("E1 SAME");
        h_eta->Draw("E1 SAME");
        h_omg->Draw("E1 SAME");

        TLine* lref = new TLine(kXmin, ref_y, kXmax, ref_y);
        lref->SetLineStyle(3); lref->SetLineColor(kGray + 2); lref->Draw();

        // Annotation: five lines, color-coded.
        TLatex tex;
        tex.SetNDC();
        tex.SetTextAlign(11);
        tex.SetTextSize(0.032);
        double yT = 0.87;
        tex.SetTextColor(kColExpFit);
        tex.DrawLatex(0.16, yT,
            Form("exp:           a = %.4f #pm %.4f", f_exp.a, f_exp.e));
        yT -= 0.042;
        tex.SetTextColor(kColSimFit);
        tex.DrawLatex(0.16, yT,
            Form("sim total:   a = %.4f #pm %.4f", f_all.a, f_all.e));
        yT -= 0.042;
        tex.SetTextColor(kColPi0Fit);
        tex.DrawLatex(0.16, yT,
            Form("sim #pi^{0} Dal.: a = %.4f #pm %.4f", f_pi0.a, f_pi0.e));
        yT -= 0.042;
        tex.SetTextColor(kColEtaFit);
        tex.DrawLatex(0.16, yT,
            Form("sim #eta Dal.: a = %.4f #pm %.4f", f_eta.a, f_eta.e));
        yT -= 0.042;
        tex.SetTextColor(kColOmegaFit);
        tex.DrawLatex(0.16, yT,
            Form("sim #omega all (41): a = %.4f #pm %.4f", f_omg.a, f_omg.e));

        // Legend (markers) — top-right corner.
        TLegend* leg = new TLegend(0.62, 0.62, 0.95, 0.88);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.030);
        leg->AddEntry(h_exp, "exp",                          "lpe");
        leg->AddEntry(h_all, "sim total",                    "lpe");
        leg->AddEntry(h_pi0, "sim #pi^{0} Dalitz (7051)",     "lpe");
        leg->AddEntry(h_eta, "sim #eta Dalitz (17051)",       "lpe");
        leg->AddEntry(h_omg, "sim #omega all (41) Dalitz + direct", "lpe");
        leg->Draw();
    };

    drawComparePad(1, r_corr_exp, r_corr_sim_all,
                   r_corr_sim_pi0, r_corr_sim_eta, r_corr_sim_omg,
                   f_corr_exp, f_corr_all, f_corr_pi0, f_corr_eta, f_corr_omg,
                   1.0, 3.5);
    drawComparePad(2, r_eff_exp,  r_eff_sim_all,
                   r_eff_sim_pi0, r_eff_sim_eta, r_eff_sim_omg,
                   f_eff_exp,  f_eff_all,  f_eff_pi0,  f_eff_eta, f_eff_omg,
                   1.0, 1.2);

    // -------- Save --------
    gSystem->mkdir("plots/output", true);
    const std::string oa_suffix =
        (oa_cut_deg > 0.0) ? Form("_oa%g", oa_cut_deg) : std::string{};
    const std::string wt_suffix = use_weights ? std::string{} : "_unw";
    const std::string fr_suffix =
        custom_range ? Form("_fit%g-%g", fxlo, fxhi) : std::string{};
    const std::string base =
        "plots/output/mass_ee_ratio_const_fit_exp_vs_sim_channels" +
        oa_suffix + wt_suffix + fr_suffix;
    c->SaveAs((base + ".pdf").c_str());
    c->SaveAs((base + ".png").c_str());
    std::cout << "\nSaved: " << base << ".{pdf,png}\n";

    fe_em->Close();
    fe_pp->Close();
    fe_mm->Close();
    fs->Close();
}
