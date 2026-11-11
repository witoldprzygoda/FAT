// isbest_emulation_analyze_betanew.C
// =========================================================================
// Analysis of the isBest ranking emulation (plots/isbest_emulation_betanew.C),
// per multiplicity class c (total identified leptons n = 2, 3, 4, >= 5):
//   * eps_data_X,c(m) = N_best / N_all     (X = pm, pp, mm; stored isBest)
//   * eps_mix_X,c(m)  = N_crosswin / N_cross (emulated win probability of a
//                       pure combinatorial pair in the real event environment)
//   * closure of the emulation on the like-sign patterns per class:
//       rho_X,c = eps_mix / eps_data
//   * like-sign identity per class: CB_c(m) = 2 sqrt(N++_c N--_c) (holds for
//     independent legs at fixed total multiplicity); summed over classes it is
//     compared with the usual k 2 sqrt(N++ N--) (k from the no-class mixing).
//   * background of the isBest scheme:
//       emulated : CB_best^emu(m)   = sum_c eps_B,c(m) CB_c(m)   (eps_B,c = eps_mix_pm,c)
//       corrected: eps_B,c -> eps_B,c / sqrt(rho_pp,c rho_mm,c)  (like-sign closure)
//       naive    : CB_best^naive(m) = 2 sqrt(N++_best N--_best)
//       per class naive: sum_c 2 sqrt(N++_best,c N--_best,c)
//     f(m) = CB_best^emu / CB_best^naive is the correction of the naive
//     isBest like-sign estimator; the per-class condition
//     eps_B,c = sqrt(eps_pp,c eps_mm,c) is printed for MID.
//   * trigger correction r(m) = 63 S_PT2 / S_PT3, pol0 on [0.14, 0.80] -> C_MID
//     for: all (k = 1), all + k, all per class (+ k), best naive, best naive
//     + k, best emulated (raw, closure-corrected, closure-corrected + k).
//   * per class: N+-_best,c versus eps_B,c CB_c.  In MID the classes c >= 3
//     are ~99% combinatorial (S_c/N+-_c ~ 1%), so this is a direct test of
//     the emulated weights; the inclusive "implied signal win probability"
//     eps_S = S_best/S_all is dominated by class 2 (eps = 1 by construction)
//     and is printed only as a bound, not as a validation.
//   * systematic variations of the emulated scheme (schemes 7-9 of the
//     trigger table): arithmetic mean N++_c + N--_c instead of the geometric
//     mean for c >= 3 (exact for "correlated pair + one extra leg"), and the
//     data weights eps_pm,c (all e+e- pairs, valid where S_c ~ 0) instead of
//     the emulated eps_B,c.
// Errors: binomial for eps (floored at 1/N); sqrt(N) for counts; the eps
// errors are propagated into CB_best^emu; k errors from the mixed counts.
// The cross pairs of one event share legs and competitors (n legs x R donors
// x partners), so the eps_mix errors are optimistic by a factor of a few;
// this is immaterial next to sqrt(N++ + N--).  The C_MID errors are
// statistical only; the spread between the schemes 6-10 (and the two class
// definitions, see the emulation macro) is the systematic of the emulated
// scheme.
//
// Usage:
//   root -l -b -q 'plots/isbest_emulation_analyze_betanew.C+("plots/output/isbest_emu_vtx_%s.root","vtx")'
//   (first argument: printf pattern with %s = pm / pp / mm)
//
// @author Witold Przygoda (witold.przygoda@uj.edu.pl)
// @date 2026

#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TString.h>
#include <TMath.h>
#include <TROOT.h>
#include <TPad.h>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace ibea {

constexpr int    kNb58   = 58;
constexpr int    kNCl    = 4;
constexpr double kMidLo  = 0.14, kMidHi = 0.80;
constexpr double kTrigCorr = 63.0;
const char* kTrig[2] = {"PT3", "PT2"};
const char* kPat[3]  = {"pm", "pp", "mm"};
const char* kCl[kNCl] = {"c2", "c3", "c4", "c5"};
const char* kClTitle[kNCl] = {"n = 2", "n = 3", "n = 4", "n #geq 5"};

TH1D* getH(TFile* f, const std::string& name, bool required = true) {
    TH1D* h = dynamic_cast<TH1D*>(f->Get(name.c_str()));
    if (!h && required) { printf("ERROR: histogram %s missing in %s\n", name.c_str(), f->GetName()); throw 1; }
    if (h) h->SetDirectory(nullptr);
    return h;
}
TH1D* cloneEmpty(TH1D* a, const char* name) { TH1D* r = (TH1D*)a->Clone(name); r->SetDirectory(nullptr); r->Reset(); return r; }
TH1D* sumClasses(TH1D* const* hc, const char* name) {
    TH1D* s = (TH1D*)hc[0]->Clone(name); s->SetDirectory(nullptr);
    for (int c = 1; c < kNCl; ++c) s->Add(hc[c]);
    return s;
}
// binomial efficiency num/den (num subset of den), error floored at 1/den
TH1D* effBinomial(TH1D* num, TH1D* den, const char* name) {
    TH1D* e = cloneEmpty(num, name);
    for (int b = 1; b <= num->GetNbinsX(); ++b) {
        const double n = num->GetBinContent(b), d = den->GetBinContent(b);
        if (d <= 0.) continue;
        const double p = n / d;
        e->SetBinContent(b, p);
        e->SetBinError(b, std::max(std::sqrt(p * (1. - p) / d), 1. / d));
    }
    return e;
}
TH1D* ratioH(TH1D* a, TH1D* b, const char* name) {
    TH1D* r = cloneEmpty(a, name);
    for (int i = 1; i <= a->GetNbinsX(); ++i) {
        const double x = a->GetBinContent(i), y = b->GetBinContent(i);
        if (x <= 0. || y <= 0.) continue;
        const double ex = a->GetBinError(i), ey = b->GetBinError(i);
        r->SetBinContent(i, x / y);
        r->SetBinError(i, x / y * std::sqrt(ex * ex / (x * x) + ey * ey / (y * y)));
    }
    return r;
}
TH1D* geoMean(TH1D* a, TH1D* b, const char* name, double scale = 1.) {   // scale * sqrt(a b)
    TH1D* r = cloneEmpty(a, name);
    for (int i = 1; i <= a->GetNbinsX(); ++i) {
        const double x = a->GetBinContent(i), y = b->GetBinContent(i);
        if (x <= 0. || y <= 0.) continue;
        const double ex = a->GetBinError(i), ey = b->GetBinError(i);
        const double g = scale * std::sqrt(x * y);
        r->SetBinContent(i, g);
        r->SetBinError(i, 0.5 * g * std::sqrt(ex * ex / (x * x) + ey * ey / (y * y)));
    }
    return r;
}
TH1D* productH(TH1D* a, TH1D* b, const char* name) {
    TH1D* r = cloneEmpty(a, name);
    for (int i = 1; i <= a->GetNbinsX(); ++i) {
        const double x = a->GetBinContent(i), y = b->GetBinContent(i);
        const double ex = a->GetBinError(i), ey = b->GetBinError(i);
        r->SetBinContent(i, x * y);
        r->SetBinError(i, std::sqrt(ex * ex * y * y + ey * ey * x * x));
    }
    return r;
}
double pol0(TH1D* h, double lo, double hi, double& err, double& chi2, int& ndf) {
    TF1 f("f_pol0_tmp", "[0]", lo, hi);
    h->Fit(&f, "RQ0N");
    err = f.GetParError(0); chi2 = f.GetChisquare(); ndf = f.GetNDF();
    return f.GetParameter(0);
}
double avgEff(TH1D* num, TH1D* den, double lo, double hi, double& err) {
    double n = 0., d = 0.;
    for (int b = 1; b <= num->GetNbinsX(); ++b) {
        const double c = num->GetBinCenter(b);
        if (c < lo || c > hi) continue;
        n += num->GetBinContent(b); d += den->GetBinContent(b);
    }
    if (d <= 0.) { err = 0.; return 0.; }
    const double p = n / d; err = std::max(std::sqrt(p * (1. - p) / d), 1. / d);
    return p;
}
double integ(TH1D* h, double lo, double hi) { return h->Integral(h->FindBin(lo + 1e-6), h->FindBin(hi - 1e-6)); }
// S = N+- - CB with independent errors
TH1D* signalH(TH1D* pm, TH1D* cb, const char* name) {
    TH1D* s = cloneEmpty(pm, name);
    for (int b = 1; b <= pm->GetNbinsX(); ++b) {
        s->SetBinContent(b, pm->GetBinContent(b) - cb->GetBinContent(b));
        s->SetBinError(b, std::sqrt(pm->GetBinContent(b) + cb->GetBinError(b) * cb->GetBinError(b)));
    }
    return s;
}
// CB = g 2 sqrt(N++ N--), error from counts and from g
TH1D* cbH(TH1D* pp, TH1D* mm, TH1D* g, const char* name) {
    TH1D* r = cloneEmpty(pp, name);
    for (int b = 1; b <= pp->GetNbinsX(); ++b) {
        const double npp = pp->GetBinContent(b), nmm = mm->GetBinContent(b);
        double gg = 1., eg = 0.;
        if (g) { gg = g->GetBinContent(b); eg = g->GetBinError(b); if (gg <= 0.) { gg = 1.; eg = 0.; } }
        const double cb = gg * 2. * std::sqrt(npp * nmm);
        double var = gg * gg * (npp + nmm);
        if (gg > 0.) var += (cb * eg / gg) * (cb * eg / gg);
        r->SetBinContent(b, cb); r->SetBinError(b, std::sqrt(var));
    }
    return r;
}
TH1D* makeRatio63(TH1D* s2, TH1D* s3, const char* name) {
    TH1D* r = cloneEmpty(s2, name);
    for (int b = 1; b <= s2->GetNbinsX(); ++b) {
        const double a = s2->GetBinContent(b), c = s3->GetBinContent(b);
        if (a <= 0. || c <= 0.) continue;
        const double ea = s2->GetBinError(b), ec = s3->GetBinError(b);
        const double v = kTrigCorr * a / c;
        r->SetBinContent(b, v);
        r->SetBinError(b, v * std::sqrt(ea * ea / (a * a) + ec * ec / (c * c)));
    }
    return r;
}
void style(TH1D* h, int color, int marker, double size = 0.9) {
    h->SetLineColor(color); h->SetMarkerColor(color); h->SetMarkerStyle(marker); h->SetMarkerSize(size);
    h->SetStats(0); h->GetXaxis()->SetTitleSize(0.05); h->GetYaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetLabelSize(0.045); h->GetYaxis()->SetLabelSize(0.045);
}
void saveC(TCanvas* c, const std::string& base) { c->SaveAs((base + ".pdf").c_str()); c->SaveAs((base + ".png").c_str()); }

} // namespace ibea

void isbest_emulation_analyze_betanew(const char* file_pattern = "plots/output/isbest_emu_vtx_%s.root",
                                      const char* tag          = "vtx",
                                      const char* mixing_file  = "plots/output/mixing_betanew.root",
                                      const char* out_dir      = "plots/output")
{
    using namespace ibea;
    gStyle->SetOptStat(0); gStyle->SetOptFit(0); gROOT->SetBatch(true);
    TFile* f[3];
    for (int x = 0; x < 3; ++x) {
        f[x] = TFile::Open(Form(file_pattern, kPat[x]));
        if (!f[x] || f[x]->IsZombie()) { printf("ERROR: cannot open %s\n", Form(file_pattern, kPat[x])); return; }
    }
    TFile* fmix = TFile::Open(mixing_file);
    if (!fmix || fmix->IsZombie()) { printf("WARNING: mixing file %s not found, k = 1\n", mixing_file); fmix = nullptr; }
    const std::string od = out_dir;

    // ---- load per class ---------------------------------------------------
    TH1D *allc[3][2][kNCl], *bestc[3][2][kNCl], *crossc[3][2][kNCl], *winc[3][2][kNCl], *nopidc[3][2][kNCl];
    TH1D *eDc[3][2][kNCl], *eMc[3][2][kNCl], *rhoc[3][2][kNCl];
    TH1D *all[3][2], *best[3][2], *cross[3][2], *win[3][2], *nopid[3][2], *eData[3][2], *eMix[3][2];
    try {
        for (int x = 0; x < 3; ++x)
            for (int t = 0; t < 2; ++t) {
                for (int c = 0; c < kNCl; ++c) {
                    allc[x][t][c]   = getH(f[x], Form("h_all_%s_%s_%s", kPat[x], kTrig[t], kCl[c]));
                    bestc[x][t][c]  = getH(f[x], Form("h_best_%s_%s_%s", kPat[x], kTrig[t], kCl[c]));
                    crossc[x][t][c] = getH(f[x], Form("h_cross_%s_%s_%s", kPat[x], kTrig[t], kCl[c]));
                    winc[x][t][c]   = getH(f[x], Form("h_crosswin_%s_%s_%s", kPat[x], kTrig[t], kCl[c]));
                    nopidc[x][t][c] = getH(f[x], Form("h_nopid_%s_%s_%s", kPat[x], kTrig[t], kCl[c]));
                    eDc[x][t][c]    = effBinomial(bestc[x][t][c], allc[x][t][c], Form("eData_%s_%s_%s", kPat[x], kTrig[t], kCl[c]));
                    eMc[x][t][c]    = effBinomial(winc[x][t][c], crossc[x][t][c], Form("eMix_%s_%s_%s", kPat[x], kTrig[t], kCl[c]));
                    rhoc[x][t][c]   = ratioH(eMc[x][t][c], eDc[x][t][c], Form("rho_%s_%s_%s", kPat[x], kTrig[t], kCl[c]));
                }
                all[x][t]   = sumClasses(allc[x][t],   Form("all_%s_%s", kPat[x], kTrig[t]));
                best[x][t]  = sumClasses(bestc[x][t],  Form("best_%s_%s", kPat[x], kTrig[t]));
                cross[x][t] = sumClasses(crossc[x][t], Form("cross_%s_%s", kPat[x], kTrig[t]));
                win[x][t]   = sumClasses(winc[x][t],   Form("win_%s_%s", kPat[x], kTrig[t]));
                nopid[x][t] = sumClasses(nopidc[x][t], Form("nopid_%s_%s", kPat[x], kTrig[t]));
                eData[x][t] = effBinomial(best[x][t], all[x][t], Form("eData_%s_%s", kPat[x], kTrig[t]));
                eMix[x][t]  = effBinomial(win[x][t], cross[x][t], Form("eMix_%s_%s", kPat[x], kTrig[t]));
            }
    } catch (...) { return; }
    for (int x = 0; x < 3; ++x) {
        TH1D* c = getH(f[x], Form("h_counters_%s", kPat[x]), false);
        if (!c) continue;
        printf("[%s] ", kPat[x]);
        for (int b = 1; b <= c->GetNbinsX(); ++b) if (c->GetXaxis()->GetBinLabel(b)[0]) printf("%s=%.0f  ", c->GetXaxis()->GetBinLabel(b), c->GetBinContent(b));
        printf("\n");
    }
    // k(m) from the no-class mixing pool
    TH1D* kfac[2] = {nullptr, nullptr};
    if (fmix) {
        for (int t = 0; t < 2; ++t) {
            TH1D* mpm = dynamic_cast<TH1D*>(fmix->Get(Form("h_mixnc_pm_%s", kTrig[t])));
            TH1D* mpp = dynamic_cast<TH1D*>(fmix->Get(Form("h_mixnc_pp_%s", kTrig[t])));
            TH1D* mmm = dynamic_cast<TH1D*>(fmix->Get(Form("h_mixnc_mm_%s", kTrig[t])));
            if (!mpm || !mpp || !mmm || mpm->GetNbinsX() != kNb58) { printf("WARNING: h_mixnc_* (%s) missing or wrong binning, k = 1\n", kTrig[t]); continue; }
            kfac[t] = cloneEmpty(all[0][t], Form("k_%s", kTrig[t]));
            for (int b = 1; b <= kNb58; ++b) {
                const double a = mpm->GetBinContent(b), p = mpp->GetBinContent(b), m = mmm->GetBinContent(b);
                if (a <= 0. || p <= 0. || m <= 0.) continue;
                const double k = a / (2. * std::sqrt(p * m));
                kfac[t]->SetBinContent(b, k);
                kfac[t]->SetBinError(b, k * std::sqrt(1. / a + 0.25 * (1. / p + 1. / m)));
            }
        }
    }

    // ---- per-class background and the isBest estimators -------------------
    TH1D *cbAllSum[2], *cbAllK[2], *cbAllClass[2], *cbBestNaive[2], *cbBestNaiveK[2], *cbBestClassNaive[2], *cbBestEmu[2], *cbBestEmuC[2], *fEmu[2], *fEmuC[2], *epsS[2];
    TH1D *cbBestAM[2], *cbBestData[2];
    TH1D *cbc[2][kNCl];
    printf("\n================ isBest emulation (%s) ================\n", tag);
    for (int t = 0; t < 2; ++t) {
        cbAllSum[t]      = cbH(all[1][t], all[2][t], nullptr, Form("cb_all_%s", kTrig[t]));
        cbAllK[t]        = cbH(all[1][t], all[2][t], kfac[t], Form("cb_all_k_%s", kTrig[t]));
        cbBestNaive[t]   = cbH(best[1][t], best[2][t], nullptr, Form("cb_best_naive_%s", kTrig[t]));
        cbBestNaiveK[t]  = cbH(best[1][t], best[2][t], kfac[t], Form("cb_best_naive_k_%s", kTrig[t]));
        cbAllClass[t]    = cloneEmpty(all[0][t], Form("cb_all_class_%s", kTrig[t]));
        cbBestClassNaive[t] = cloneEmpty(all[0][t], Form("cb_best_class_naive_%s", kTrig[t]));
        cbBestEmu[t]     = cloneEmpty(all[0][t], Form("cb_best_emu_%s", kTrig[t]));
        cbBestEmuC[t]    = cloneEmpty(all[0][t], Form("cb_best_emu_corr_%s", kTrig[t]));
        cbBestAM[t]   = cloneEmpty(all[0][t], Form("cb_best_emu_am_%s", kTrig[t]));
        cbBestData[t] = cloneEmpty(all[0][t], Form("cb_best_dataw_%s", kTrig[t]));
        for (int c = 0; c < kNCl; ++c) {
            cbc[t][c] = cbH(allc[1][t][c], allc[2][t][c], nullptr, Form("cb_%s_%s", kTrig[t], kCl[c]));
            cbAllClass[t]->Add(cbc[t][c]);
            // arithmetic-mean variant for c >= 3 (c2 keeps the geometric mean)
            TH1D* cbam = c == 0 ? (TH1D*)cbc[t][c]->Clone(Form("cbam_%s_%s", kTrig[t], kCl[c])) : (TH1D*)allc[1][t][c]->Clone(Form("cbam_%s_%s", kTrig[t], kCl[c]));
            cbam->SetDirectory(nullptr);
            if (c > 0) cbam->Add(allc[2][t][c]);
            cbBestAM[t]->Add(productH(eMc[0][t][c], cbam, Form("cbbest_am_%s_%s", kTrig[t], kCl[c])));
            // data-weight variant: eps_pm,c (all e+e- pairs of the class)
            cbBestData[t]->Add(productH(eDc[0][t][c], cbc[t][c], Form("cbbest_dataw_%s_%s", kTrig[t], kCl[c])));
            TH1D* nb = cbH(bestc[1][t][c], bestc[2][t][c], nullptr, Form("cbbest_naive_%s_%s", kTrig[t], kCl[c]));
            cbBestClassNaive[t]->Add(nb);
            TH1D* e = productH(eMc[0][t][c], cbc[t][c], Form("cbbest_emu_%s_%s", kTrig[t], kCl[c]));
            cbBestEmu[t]->Add(e);
            TH1D* sr = geoMean(rhoc[1][t][c], rhoc[2][t][c], Form("sqrtrho_%s_%s", kTrig[t], kCl[c]));
            TH1D* ec = ratioH(eMc[0][t][c], sr, Form("epsBc_%s_%s", kTrig[t], kCl[c]));   // closure-corrected eps_B,c
            // bins where the closure is undefined (no like-sign pairs): fall back to the raw eps
            for (int b = 1; b <= kNb58; ++b) if (ec->GetBinContent(b) <= 0.) { ec->SetBinContent(b, eMc[0][t][c]->GetBinContent(b)); ec->SetBinError(b, eMc[0][t][c]->GetBinError(b)); }
            TH1D* e2 = productH(ec, cbc[t][c], Form("cbbest_emuc_%s_%s", kTrig[t], kCl[c]));
            cbBestEmuC[t]->Add(e2);
        }
        fEmu[t]  = ratioH(cbBestEmu[t],  cbBestNaive[t], Form("f_emu_%s", kTrig[t]));
        fEmuC[t] = ratioH(cbBestEmuC[t], cbBestNaive[t], Form("f_emu_corr_%s", kTrig[t]));
        epsS[t] = cloneEmpty(all[0][t], Form("epsS_%s", kTrig[t]));
        for (int b = 1; b <= kNb58; ++b) {
            const double npm = all[0][t]->GetBinContent(b), cb = cbAllClass[t]->GetBinContent(b), s = npm - cb;
            if (s <= 0. || npm <= 0.) continue;
            epsS[t]->SetBinContent(b, (eData[0][t]->GetBinContent(b) * npm - cbBestEmuC[t]->GetBinContent(b)) / s);
        }

        // ---- printout ----
        double e1, e2, ch; int nd;
        printf("\n--- %s ---\n", kTrig[t]);
        printf("  class shares of MID pairs [%.2f,%.2f]: ", kMidLo, kMidHi);
        for (int c = 0; c < kNCl; ++c) printf("%s: e+e- %.1f%% (cross %.1f%%)  e+e+ %.1f%%  e-e- %.1f%% | ", kCl[c],
            100. * integ(allc[0][t][c], kMidLo, kMidHi) / std::max(1., integ(all[0][t], kMidLo, kMidHi)),
            100. * integ(crossc[0][t][c], kMidLo, kMidHi) / std::max(1., integ(cross[0][t], kMidLo, kMidHi)),
            100. * integ(allc[1][t][c], kMidLo, kMidHi) / std::max(1., integ(all[1][t], kMidLo, kMidHi)),
            100. * integ(allc[2][t][c], kMidLo, kMidHi) / std::max(1., integ(all[2][t], kMidLo, kMidHi)));
        printf("\n  pair-weighted MID averages per class (eps_data / eps_mix):\n");
        for (int c = 0; c < kNCl; ++c) {
            printf("   %s:", kCl[c]);
            for (int x = 0; x < 3; ++x) {
                const double ed = avgEff(bestc[x][t][c], allc[x][t][c], kMidLo, kMidHi, e1);
                const double em = avgEff(winc[x][t][c], crossc[x][t][c], kMidLo, kMidHi, e2);
                printf("  %s %.4f+-%.4f / %.4f+-%.4f", kPat[x], ed, e1, em, e2);
            }
            double epp, emm, eB, eB2;
            const double epsPP = avgEff(bestc[1][t][c], allc[1][t][c], kMidLo, kMidHi, epp);
            const double epsMM = avgEff(bestc[2][t][c], allc[2][t][c], kMidLo, kMidHi, emm);
            const double epsB  = avgEff(winc[0][t][c], crossc[0][t][c], kMidLo, kMidHi, eB);
            const double epsPPm = avgEff(winc[1][t][c], crossc[1][t][c], kMidLo, kMidHi, eB2);
            const double epsMMm = avgEff(winc[2][t][c], crossc[2][t][c], kMidLo, kMidHi, eB2);
            const double rr = (epsPPm > 0 && epsMMm > 0 && epsPP > 0 && epsMM > 0) ? std::sqrt(epsPPm / epsPP * epsMMm / epsMM) : 0.;
            printf("\n        condition: eps_B,c = %.4f  vs sqrt(eps_pp,c eps_mm,c) data = %.4f  (ratio %.4f; closure sqrt(rho_pp rho_mm) = %.4f -> corrected ratio %.4f)\n",
                   epsB, std::sqrt(epsPP * epsMM), epsPP * epsMM > 0 ? epsB / std::sqrt(epsPP * epsMM) : 0., rr, (rr > 0 && epsPP * epsMM > 0) ? epsB / std::sqrt(epsPP * epsMM) / rr : 0.);
        }
        const double mAll = integ(all[0][t], kMidLo, kMidHi), mBest = integ(best[0][t], kMidLo, kMidHi);
        printf("  MID sums: N+-_all %.0f  N+-_best %.0f | CB_all 2sqrt %.0f, with k %.0f, per class %.0f | CB_best naive %.0f, naive+k %.0f, per-class naive %.0f, emulated %.0f, emulated+closure %.0f\n",
               mAll, mBest, integ(cbAllSum[t], kMidLo, kMidHi), integ(cbAllK[t], kMidLo, kMidHi), integ(cbAllClass[t], kMidLo, kMidHi),
               integ(cbBestNaive[t], kMidLo, kMidHi), integ(cbBestNaiveK[t], kMidLo, kMidHi), integ(cbBestClassNaive[t], kMidLo, kMidHi), integ(cbBestEmu[t], kMidLo, kMidHi), integ(cbBestEmuC[t], kMidLo, kMidHi));
        const double fe = pol0(fEmu[t], kMidLo, kMidHi, e1, ch, nd);
        printf("  f = CB_best^emu / CB_best^naive (pol0 MID) = %.4f +- %.4f (chi2/ndf %.1f/%d)", fe, e1, ch, nd);
        const double fc = pol0(fEmuC[t], kMidLo, kMidHi, e2, ch, nd);
        printf(";  closure-corrected = %.4f +- %.4f (chi2/ndf %.1f/%d)\n", fc, e2, ch, nd);
        if (kfac[t]) { double ek; const double kk = pol0(kfac[t], kMidLo, kMidHi, ek, ch, nd); printf("  k (no-class mixing, pol0 MID) = %.4f +- %.4f;  CB_all per class / CB_all 2sqrt (MID sums) = %.4f\n", kk, ek, integ(cbAllClass[t], kMidLo, kMidHi) / integ(cbAllSum[t], kMidLo, kMidHi)); }
        printf("  cross pairs of active legs failing the recomputed PID: %.3f%% (all masses)\n", 100. * nopid[0][t]->Integral() / std::max(1., nopid[0][t]->Integral() + cross[0][t]->Integral()));
        printf("  per-class test in MID (classes >= 3 are ~99%% combinatorial): N+-_all,c | CB_c (GM) | CB_c (AM) | S_all,c | N+-_best,c | eps_B,c CB_c (emu) | eps_pm,c CB_c (data weights) | naive 2sqrt(best)\n");
        for (int c = 0; c < kNCl; ++c) {
            const double na = integ(allc[0][t][c], kMidLo, kMidHi), cb = integ(cbc[t][c], kMidLo, kMidHi);
            const double am = integ(allc[1][t][c], kMidLo, kMidHi) + integ(allc[2][t][c], kMidLo, kMidHi);
            const double nb = integ(bestc[0][t][c], kMidLo, kMidHi);
            double e1_, e2_;
            const double eB = avgEff(winc[0][t][c], crossc[0][t][c], kMidLo, kMidHi, e1_);
            const double eD = avgEff(bestc[0][t][c], allc[0][t][c], kMidLo, kMidHi, e2_);
            const double nv = 2. * std::sqrt(integ(bestc[1][t][c], kMidLo, kMidHi) * integ(bestc[2][t][c], kMidLo, kMidHi));
            printf("     %s: %9.0f | %9.0f | %9.0f | %8.0f | %9.0f | %9.0f | %9.0f | %9.0f\n", kCl[c], na, cb, am, na - cb, nb, eB * cb, eD * cb, nv);
        }
        printf("  inclusive S_best/S_all (class-2 dominated, a bound only) in MID bins: ");
        for (int b = 1; b <= kNb58; ++b) { const double c = epsS[t]->GetBinCenter(b); if (c > kMidLo && c < kMidHi) printf("%.2f ", epsS[t]->GetBinContent(b)); }
        printf("\n");
    }

    // ---- trigger correction, several CB schemes ---------------------------
    printf("\n--- trigger correction C_MID (pol0 on [%.2f, %.2f] of 63 S_PT2/S_PT3) ---\n", kMidLo, kMidHi);
    // the acceptance factor k (e+ and e- single-leg spectra differ, so the
    // geometric mean of the like-sign spectra is not the e+e- combinatorial
    // spectrum even for independent legs) applies to the per-class sums as well
    TH1D *cbAllClassK[2], *cbBestEmuCK[2];
    for (int t = 0; t < 2; ++t) {
        cbAllClassK[t] = kfac[t] ? productH(cbAllClass[t], kfac[t], Form("cb_all_class_k_%s", kTrig[t])) : cbAllClass[t];
        cbBestEmuCK[t] = kfac[t] ? productH(cbBestEmuC[t], kfac[t], Form("cb_best_emu_corr_k_%s", kTrig[t])) : cbBestEmuC[t];
    }
    TH1D *cbBestAMK[2], *cbBestDataK[2];
    for (int t = 0; t < 2; ++t) {
        cbBestAMK[t]   = kfac[t] ? productH(cbBestAM[t], kfac[t], Form("cb_best_emu_am_k_%s", kTrig[t])) : cbBestAM[t];
        cbBestDataK[t] = kfac[t] ? productH(cbBestData[t], kfac[t], Form("cb_best_dataw_k_%s", kTrig[t])) : cbBestData[t];
    }
    const int nSch = 11;
    const char* schName[nSch] = {"all, k=1", "all + k", "all per class", "all per class + k", "best naive", "best naive + k", "best emulated", "best emulated+closure", "best emulated+closure+k", "best emulated AM(c>=3)+k", "best data-weights eps_pm,c+k"};
    TH1D* rSch[nSch];
    for (int s = 0; s < nSch; ++s) {
        TH1D* sig[2];
        for (int t = 0; t < 2; ++t) {
            TH1D* cb = s == 0 ? cbAllSum[t] : s == 1 ? cbAllK[t] : s == 2 ? cbAllClass[t] : s == 3 ? cbAllClassK[t] : s == 4 ? cbBestNaive[t] : s == 5 ? cbBestNaiveK[t] : s == 6 ? cbBestEmu[t] : s == 7 ? cbBestEmuC[t] : s == 8 ? cbBestEmuCK[t] : s == 9 ? cbBestAMK[t] : cbBestDataK[t];
            TH1D* pm_ = s <= 3 ? all[0][t] : best[0][t];
            sig[t] = signalH(pm_, cb, Form("sig_%d_%s", s, kTrig[t]));
        }
        rSch[s] = makeRatio63(sig[1], sig[0], Form("ratio_%d", s));
        double err, ch; int nd;
        const double c = pol0(rSch[s], kMidLo, kMidHi, err, ch, nd);
        printf("  %-24s : C_MID = %.4f +- %.4f  (chi2/ndf %.1f/%d)   S_MID PT3 = %.0f  PT2 = %.0f\n", schName[s], c, err, ch, nd, integ(sig[0], kMidLo, kMidHi), integ(sig[1], kMidLo, kMidHi));
    }

    // ---- figures ---------------------------------------------------------
    const int clCol[kNCl] = {kGray + 2, kRed + 1, kBlue + 1, kGreen + 2};
    // Fig 1: like-sign closure per class (c3..c5), both patterns, one canvas per trigger
    for (int t = 0; t < 2; ++t) {
        TCanvas* c1 = new TCanvas(Form("c_closure_%s", kTrig[t]), "", 1400, 600); c1->Divide(2, 1);
        for (int x = 1; x <= 2; ++x) {
            c1->cd(x);
            TH1D* fr = cloneEmpty(eDc[x][t][0], Form("fr_cl_%d_%d", t, x));
            fr->SetTitle(Form("%s %s: isBest win probability per multiplicity class;M_{ee} [GeV/c^{2}];#varepsilon", x == 1 ? "e^{+}e^{+}" : "e^{-}e^{-}", kTrig[t]));
            fr->GetYaxis()->SetRangeUser(0., 1.1); fr->GetXaxis()->SetRangeUser(0., 1.0); style(fr, kBlack, 1); fr->Draw();
            TLegend* l = new TLegend(0.5, 0.12, 0.9, 0.45); l->SetBorderSize(0); l->SetFillStyle(0);
            for (int c = 1; c < kNCl; ++c) {
                style(eDc[x][t][c], clCol[c], 20); style(eMc[x][t][c], clCol[c], 24);
                eDc[x][t][c]->Draw("same e"); eMc[x][t][c]->Draw("same e");
                l->AddEntry(eDc[x][t][c], Form("%s data (isBest)", kClTitle[c]), "lpe");
                l->AddEntry(eMc[x][t][c], Form("%s emulation", kClTitle[c]), "lpe");
            }
            l->Draw();
        }
        saveC(c1, od + "/isbest_emu_closure_" + kTrig[t] + "_" + tag);
    }
    // Fig 2: e+e- per class: eps_B,c (emulation) vs sqrt(eps_pp,c eps_mm,c) (data) and eps_data_pm,c
    for (int t = 0; t < 2; ++t) {
        TCanvas* c2 = new TCanvas(Form("c_epsB_%s", kTrig[t]), "", 1400, 600); c2->Divide(2, 1);
        for (int pane = 1; pane <= 2; ++pane) {
            c2->cd(pane);
            TH1D* fr = cloneEmpty(eDc[0][t][0], Form("fr_eb_%d_%d", t, pane));
            fr->SetTitle(Form("e^{+}e^{-} %s, %s;M_{ee} [GeV/c^{2}];#varepsilon", kTrig[t], pane == 1 ? "n = 3" : "n = 4, n #geq 5"));
            fr->GetYaxis()->SetRangeUser(0., 1.1); fr->GetXaxis()->SetRangeUser(0., 1.0); style(fr, kBlack, 1); fr->Draw();
            TLegend* l = new TLegend(0.45, 0.12, 0.9, 0.42); l->SetBorderSize(0); l->SetFillStyle(0);
            for (int c = (pane == 1 ? 1 : 2); c < (pane == 1 ? 2 : kNCl); ++c) {
                TH1D* sq = geoMean(eDc[1][t][c], eDc[2][t][c], Form("sqrtData_%s_%s", kTrig[t], kCl[c]));
                style(eMc[0][t][c], clCol[c], 24); style(sq, clCol[c], 22); style(eDc[0][t][c], clCol[c], 20, 0.7);
                eMc[0][t][c]->Draw("same e"); sq->Draw("same e"); eDc[0][t][c]->Draw("same e");
                l->AddEntry(eMc[0][t][c], Form("%s #varepsilon_{B} emulation (combinatorial)", kClTitle[c]), "lpe");
                l->AddEntry(sq, Form("%s #sqrt{#varepsilon_{++}#varepsilon_{--}} data", kClTitle[c]), "lpe");
                l->AddEntry(eDc[0][t][c], Form("%s #varepsilon_{+-} data (all pairs)", kClTitle[c]), "lpe");
            }
            l->Draw();
        }
        saveC(c2, od + "/isbest_emu_epsB_" + kTrig[t] + "_" + tag);
    }
    // Fig 3: f(m) for both triggers
    {
        TCanvas* c3 = new TCanvas("c_f", "", 1000, 700);
        TH1D* fr = cloneEmpty(fEmu[0], "fr_f"); fr->SetTitle("isBest scheme: emulated / naive like-sign background;M_{ee} [GeV/c^{2}];f = CB_{best}^{emu} / 2#sqrt{N_{++}^{best}N_{--}^{best}}");
        fr->GetYaxis()->SetRangeUser(0.6, 1.6); fr->GetXaxis()->SetRangeUser(0., 1.0); style(fr, kBlack, 1); fr->Draw();
        style(fEmu[0], kBlue + 1, 24); style(fEmuC[0], kBlue + 1, 20); style(fEmu[1], kRed + 1, 25); style(fEmuC[1], kRed + 1, 21);
        fEmu[0]->Draw("same e"); fEmuC[0]->Draw("same e"); fEmu[1]->Draw("same e"); fEmuC[1]->Draw("same e");
        TLine* ln = new TLine(0., 1., 1.0, 1.); ln->SetLineStyle(2); ln->Draw();
        TLegend* l = new TLegend(0.15, 0.65, 0.6, 0.9); l->SetBorderSize(0); l->SetFillStyle(0);
        l->AddEntry(fEmu[0], "PT3 raw", "lpe"); l->AddEntry(fEmuC[0], "PT3 closure-corrected", "lpe");
        l->AddEntry(fEmu[1], "PT2 raw", "lpe"); l->AddEntry(fEmuC[1], "PT2 closure-corrected", "lpe"); l->Draw();
        saveC(c3, od + "/isbest_emu_f_" + tag);
    }
    // Fig 4: trigger ratio for the schemes
    {
        TCanvas* c4 = new TCanvas("c_ratio", "", 1100, 750);
        TH1D* fr = cloneEmpty(rSch[0], "fr5"); fr->SetTitle("63 S_{PT2}/S_{PT3};M_{ee} [GeV/c^{2}];63 S_{PT2}/S_{PT3}");
        fr->GetYaxis()->SetRangeUser(0., 4.); fr->GetXaxis()->SetRangeUser(0.1, 1.0); style(fr, kBlack, 1); fr->Draw();
        const int cs[nSch] = {kGray + 2, kBlack, kGray + 1, kGray + 3, kRed + 1, kOrange + 7, kBlue + 1, kGreen + 2, kMagenta + 2, kCyan + 2, kViolet - 1};
        const int ms[nSch] = {24, 20, 26, 32, 25, 27, 21, 22, 23, 33, 34};
        TLegend* l = new TLegend(0.13, 0.58, 0.66, 0.9); l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(0.025);
        for (int s = 0; s < nSch; ++s) {
            style(rSch[s], cs[s], ms[s]); rSch[s]->Draw("same e");
            double err, ch; int nd; const double c = pol0(rSch[s], kMidLo, kMidHi, err, ch, nd);
            TLine* ln = new TLine(kMidLo, c, kMidHi, c); ln->SetLineColor(cs[s]); ln->SetLineWidth(2); ln->Draw();
            l->AddEntry(rSch[s], Form("%s: C_{MID} = %.3f #pm %.3f", schName[s], c, err), "lpe");
        }
        l->Draw();
        saveC(c4, od + "/isbest_emu_trigger_ratio_" + tag);
    }
    // Fig 5: diagnostics D and chi2, e-e- PT3 native vs cross
    {
        TCanvas* c5 = new TCanvas("c_diag", "", 1400, 600); c5->Divide(2, 1);
        TH1D* dn = getH(f[2], "h_D_native_mm_PT3"); TH1D* dc = getH(f[2], "h_D_cross_mm_PT3");
        TH1D* cn = getH(f[2], "h_chi2_native_mm_PT3"); TH1D* cc = getH(f[2], "h_chi2_cross_mm_PT3");
        c5->cd(1); gPad->SetLogy();
        dn->Scale(1. / std::max(1., dn->Integral())); dc->Scale(1. / std::max(1., dc->Integral()));
        style(dn, kBlack, 20, 0.6); style(dc, kRed + 1, 24, 0.6); dn->SetTitle("e^{-}e^{-} PT3, 0.14<M<0.80: D = (tof_{1}-tof_{2}) - (tof_{mom,1}-tof_{mom,2});D [ns];normalised");
        dn->Draw("e"); dc->Draw("same e");
        TLegend* l = new TLegend(0.6, 0.75, 0.9, 0.9); l->SetBorderSize(0); l->SetFillStyle(0); l->AddEntry(dn, "native pairs", "lpe"); l->AddEntry(dc, "cross pairs (emulated)", "lpe"); l->Draw();
        c5->cd(2); gPad->SetLogy();
        cn->Scale(1. / std::max(1., cn->Integral())); cc->Scale(1. / std::max(1., cc->Integral()));
        style(cn, kBlack, 20, 0.6); style(cc, kRed + 1, 24, 0.6); cn->SetTitle("e^{-}e^{-} PT3, 0.14<M<0.80: ranking #chi^{2};log_{10}#chi^{2};normalised");
        cn->Draw("e"); cc->Draw("same e");
        saveC(c5, od + "/isbest_emu_diag_mm_PT3_" + tag);
    }
    // ---- store ------------------------------------------------------------
    TFile* fo = new TFile((od + "/isbest_emu_factors_" + tag + ".root").c_str(), "recreate");
    for (int t = 0; t < 2; ++t) {
        for (int x = 0; x < 3; ++x) { eData[x][t]->Write(); eMix[x][t]->Write(); for (int c = 0; c < kNCl; ++c) { eDc[x][t][c]->Write(); eMc[x][t][c]->Write(); rhoc[x][t][c]->Write(); } }
        for (int c = 0; c < kNCl; ++c) cbc[t][c]->Write();
        cbAllSum[t]->Write(); cbAllK[t]->Write(); cbAllClass[t]->Write(); cbBestNaive[t]->Write(); cbBestNaiveK[t]->Write(); cbBestClassNaive[t]->Write();
        cbBestEmu[t]->Write(); cbBestEmuC[t]->Write(); cbBestAM[t]->Write(); cbBestData[t]->Write(); fEmu[t]->Write(); fEmuC[t]->Write(); epsS[t]->Write(); if (kfac[t]) { kfac[t]->Write(); cbAllClassK[t]->Write(); cbBestEmuCK[t]->Write(); }
    }
    for (int s = 0; s < nSch; ++s) rSch[s]->Write();
    fo->Close();
    printf("\nwritten: %s/isbest_emu_{closure,epsB}_{PT3,PT2}_%s.{pdf,png}, isbest_emu_f_%s, isbest_emu_trigger_ratio_%s, isbest_emu_diag_mm_PT3_%s, isbest_emu_factors_%s.root\n", od.c_str(), tag, tag, tag, tag, tag);
}
