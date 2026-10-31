// mixing_analyze_betanew.C
// =========================================================================
// Analysis of the event-mixing output of mixing_build_betanew.C:
//
//   for each trigger (PT3, PT2) and scheme (all = all pair combinations,
//   best = isBest==1 pairs / "best legs" mixing):
//   (a) k-factor k(m) = N_mix+-(m) / (2 sqrt(N_mix++(m) N_mix--(m))) on the
//       58-bin recipe grid and on 5 MeV bins, with errors; global k in the
//       windows LOW [0, 0.14], MID [0.14, 0.80], HIGH [0.80, 1.40];
//   (b) shape: same-event CB_ls = 2 sqrt(N_se++ N_se--) vs N_mix+-
//       normalised to it in [norm_lo, norm_hi] (default [0.30, 0.80]:
//       above the pi0/eta Dalitz region, inside the window used for the
//       trigger ratio, where S/CB is smallest); ratio vs m; also the
//       like-sign spectra ++ and -- individually vs their mixed shapes;
//   (c) signals: S_ls = N+- - CB_ls, S_k = N+- - k(m) CB_ls,
//       S_mix = N+- - alpha N_mix+-;
//   (d) r = 63 S_PT2 / S_PT3 on the 58-bin grid for S_ls, S_k, S_mix,
//       weighted pol0 in MID [0.14, 0.80] (C_MID) and width-weighted means
//       LOW/MID/HIGH (same code as plots/mass_ee_pt_ratio_betanew.C);
//   (e) plots to <out_dir>/mixing_<tag>_*.{pdf,png}.
//
// Internal check (full data, scheme all): S_ls PT3 MID must reproduce
// 768 003 and C_MID = 1.607 of mass_ee_pt_ratio_betanew.C.
//
// mix_kind: "mix"   = pools in classes trigger x n+ x n- (specification);
//           "mixnc" = pools per trigger only (no multiplicity classes).
// Output names: <out_dir>/mixing_<tag>_<mix_kind>_{k,shape,ratio}_{all,best}.{pdf,png}
//
// Usage (ACLiC):
//   root -l -b -q 'plots/mixing_analyze_betanew.C+("plots/output/mixing_betanew.root","betanew")'
//   root -l -b -q 'plots/mixing_analyze_betanew.C+("plots/output/mixing_betanew.root","betanew","plots/output","mixnc")'
//
// @author Witold Przygoda (witold.przygoda@uj.edu.pl)
// @date 2026

#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TBox.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TParameter.h>
#include <TString.h>

#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

namespace mixa {

constexpr double kTrigCorr = 63.0;
constexpr double kLowHi  = 0.14;
constexpr double kHighLo = 0.80;
constexpr double kXmax   = 1.40;
const char* trigName(int t) { return t == 0 ? "PT3" : "PT2"; }
const char* schemeName(int s) { return s == 0 ? "all" : "best"; }

TH1D* getH(TFile* f, const std::string& name) {
    TH1D* h = dynamic_cast<TH1D*>(f->Get(name.c_str()));
    if (!h) { std::cerr << "missing histogram " << name << "\n"; return nullptr; }
    h->SetDirectory(nullptr);
    return h;
}

long long getL(TFile* f, const char* name) {
    TParameter<Long64_t>* p = dynamic_cast<TParameter<Long64_t>*>(f->Get(name));
    return p ? static_cast<long long>(p->GetVal()) : -1;
}

TH1D* cloneEmpty(TH1D* h, const std::string& name) {
    TH1D* c = static_cast<TH1D*>(h->Clone(name.c_str()));
    c->SetDirectory(nullptr);
    c->Reset();
    return c;
}

// CB = 2 sqrt(N++ N--) (as in the recipe macro).
TH1D* makeCB(TH1D* h_pp, TH1D* h_mm, const std::string& name) {
    TH1D* cb = cloneEmpty(h_pp, name);
    for (int b = 1; b <= h_pp->GetNbinsX(); ++b) {
        const double n_pp = h_pp->GetBinContent(b), e_pp = h_pp->GetBinError(b);
        const double n_mm = h_mm->GetBinContent(b), e_mm = h_mm->GetBinError(b);
        if (n_pp <= 0.0 || n_mm <= 0.0) { cb->SetBinContent(b, 0.0); cb->SetBinError(b, 0.0); continue; }
        const double val = 2.0 * std::sqrt(n_pp * n_mm);
        const double err = val * 0.5 * std::sqrt(std::pow(e_pp / n_pp, 2) + std::pow(e_mm / n_mm, 2));
        cb->SetBinContent(b, val);
        cb->SetBinError(b, err);
    }
    return cb;
}

// k(m) = N+- / (2 sqrt(N++ N--)) from three independent Poisson spectra.
TH1D* makeK(TH1D* pm, TH1D* pp, TH1D* mm, const std::string& name) {
    TH1D* k = cloneEmpty(pm, name);
    for (int b = 1; b <= pm->GetNbinsX(); ++b) {
        const double a = pm->GetBinContent(b), p = pp->GetBinContent(b), q = mm->GetBinContent(b);
        if (a <= 0.0 || p <= 0.0 || q <= 0.0) { k->SetBinContent(b, 0.0); k->SetBinError(b, 0.0); continue; }
        const double val = a / (2.0 * std::sqrt(p * q));
        const double rel = std::sqrt(1.0 / a + 0.25 / p + 0.25 / q);
        k->SetBinContent(b, val);
        k->SetBinError(b, val * rel);
    }
    return k;
}

// Integral and its error (sqrt of sum of squared bin errors) in [xlo, xhi].
void integ(TH1D* h, double xlo, double xhi, double& sum, double& err) {
    sum = 0.0; err = 0.0;
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
        const double lo = h->GetXaxis()->GetBinLowEdge(b), hi = h->GetXaxis()->GetBinUpEdge(b);
        if (lo < xlo - 1e-9 || hi > xhi + 1e-9) continue;
        sum += h->GetBinContent(b);
        err += std::pow(h->GetBinError(b), 2);
    }
    err = std::sqrt(err);
}

// global k in a window from integrated counts
void kWindow(TH1D* pm, TH1D* pp, TH1D* mm, double xlo, double xhi, double& k, double& ek) {
    double a, ea, p, ep, q, eq;
    integ(pm, xlo, xhi, a, ea); integ(pp, xlo, xhi, p, ep); integ(mm, xlo, xhi, q, eq);
    if (a <= 0 || p <= 0 || q <= 0) { k = 0; ek = 0; return; }
    k = a / (2.0 * std::sqrt(p * q));
    ek = k * std::sqrt(1.0 / a + 0.25 / p + 0.25 / q);
}

// scaled copy: alpha * h, with alpha +- ealpha (fully correlated across bins)
TH1D* scaled(TH1D* h, double alpha, double ealpha, const std::string& name) {
    TH1D* s = cloneEmpty(h, name);
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
        const double n = h->GetBinContent(b), e = h->GetBinError(b);
        s->SetBinContent(b, alpha * n);
        s->SetBinError(b, std::sqrt(std::pow(alpha * e, 2) + std::pow(n * ealpha, 2)));
    }
    return s;
}

// product k(m) * CB(m), independent errors
TH1D* product(TH1D* k, TH1D* cb, const std::string& name) {
    TH1D* p = cloneEmpty(cb, name);
    for (int b = 1; b <= cb->GetNbinsX(); ++b) {
        const double kv = k->GetBinContent(b), ke = k->GetBinError(b);
        const double cv = cb->GetBinContent(b), ce = cb->GetBinError(b);
        p->SetBinContent(b, kv * cv);
        p->SetBinError(b, std::sqrt(std::pow(kv * ce, 2) + std::pow(cv * ke, 2)));
    }
    return p;
}

TH1D* makeSignal(TH1D* h_em, TH1D* h_cb, const std::string& name) {
    TH1D* sig = cloneEmpty(h_em, name);
    for (int b = 1; b <= h_em->GetNbinsX(); ++b) {
        sig->SetBinContent(b, h_em->GetBinContent(b) - h_cb->GetBinContent(b));
        sig->SetBinError(b, std::sqrt(std::pow(h_em->GetBinError(b), 2) + std::pow(h_cb->GetBinError(b), 2)));
    }
    return sig;
}

TH1D* makeRatio(TH1D* num, TH1D* den, double factor, const std::string& name) {
    TH1D* r = cloneEmpty(num, name);
    for (int b = 1; b <= num->GetNbinsX(); ++b) {
        const double n = num->GetBinContent(b), en = num->GetBinError(b);
        const double d = den->GetBinContent(b), ed = den->GetBinError(b);
        if (d <= 0.0 || !std::isfinite(d) || !std::isfinite(n)) { r->SetBinContent(b, 0.0); r->SetBinError(b, 0.0); continue; }
        const double val = factor * n / d;
        const double rel = std::sqrt((n != 0.0 ? std::pow(en / n, 2) : 0.0) + std::pow(ed / d, 2));
        r->SetBinContent(b, val);
        r->SetBinError(b, std::abs(val) * rel);
    }
    return r;
}

// width-weighted mean as in the recipe macro
void regionMean(TH1D* h, double xlo, double xhi, double& mean, double& sigma, int& nb) {
    double sum_vw = 0.0, sum_w = 0.0, sum_e2w2 = 0.0;
    nb = 0;
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
        const double lo = h->GetXaxis()->GetBinLowEdge(b), hi = h->GetXaxis()->GetBinUpEdge(b);
        if (lo < xlo - 1e-9 || hi > xhi + 1e-9) continue;
        const double w = h->GetBinWidth(b), v = h->GetBinContent(b), e = h->GetBinError(b);
        if (v == 0.0 && e == 0.0) continue;
        sum_vw += v * w; sum_w += w; sum_e2w2 += (e * w) * (e * w);
        ++nb;
    }
    mean = sum_w > 0 ? sum_vw / sum_w : 0.0;
    sigma = sum_w > 0 ? std::sqrt(sum_e2w2) / sum_w : 0.0;
}

void fitMid(TH1D* r, double& C, double& eC, double& chi2, int& ndf) {
    TF1 f("pol0_mid", "[0]", kLowHi, kHighLo);
    r->Fit(&f, "RQ0");
    C = f.GetParameter(0); eC = f.GetParError(0); chi2 = f.GetChisquare(); ndf = f.GetNDF();
    if (r->GetListOfFunctions()) r->GetListOfFunctions()->Clear();
}

void style(TH1D* h, Color_t c, Style_t m, double msize = 0.8) {
    h->SetLineColor(c); h->SetMarkerColor(c); h->SetMarkerStyle(m); h->SetMarkerSize(msize);
    h->SetLineWidth(1);
    h->GetXaxis()->SetTitleSize(0.045); h->GetYaxis()->SetTitleSize(0.045);
    h->GetXaxis()->SetLabelSize(0.040); h->GetYaxis()->SetLabelSize(0.040);
    h->GetXaxis()->SetTitleOffset(1.05); h->GetYaxis()->SetTitleOffset(1.25);
}

struct Set {
    // 58-bin grid
    TH1D *se_pm, *se_pp, *se_mm, *mix_pm, *mix_pp, *mix_mm;
    // fine grid
    TH1D *se_pmF, *se_ppF, *se_mmF, *mix_pmF, *mix_ppF, *mix_mmF;
    TH1D *k, *kF, *cb_ls, *cb_mix, *cb_k, *s_ls, *s_k, *s_mix;
    TH1D *shape_pm, *shape_pp, *shape_mm;       // same-event / normalised mixed
    double alpha, ealpha, k_win[3], ek_win[3];
};

}  // namespace mixa

using namespace mixa;

void mixing_analyze_betanew(const char* in_file = "plots/output/mixing_betanew.root",
                            const char* tag = "betanew",
                            const char* out_dir = "plots/output",
                            const char* mix_kind = "mix",
                            double norm_lo = 0.30,
                            double norm_hi = 0.80) {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetErrorX(0.0);
    gSystem->mkdir(out_dir, true);

    TFile* f = TFile::Open(in_file, "READ");
    if (!f || f->IsZombie()) { std::cerr << "cannot open " << in_file << "\n"; return; }
    std::printf("Input: %s\n", in_file);
    std::printf("  events=%lld (e+e- %lld, e+e+ only %lld, e-e- only %lld)  PT3=%lld PT2=%lld  runs=%lld [%lld, %lld]  pool_depth=%lld  max_rows=%lld\n",
                getL(f, "n_events"), getL(f, "n_events_epem"), getL(f, "n_events_epep_only"),
                getL(f, "n_events_emem_only"), getL(f, "n_events_PT3"), getL(f, "n_events_PT2"),
                getL(f, "n_runs"), getL(f, "first_run"), getL(f, "last_run"),
                getL(f, "pool_depth"), getL(f, "max_rows"));
    std::printf("  rows: epem=%lld epep=%lld emem=%lld; sanity: n11=%lld bad=%lld; rows!=n+n-: %lld; both-charge legs: %lld\n",
                getL(f, "rows_epem"), getL(f, "rows_epep"), getL(f, "rows_emem"),
                getL(f, "n_events_11"), getL(f, "n_events_11_bad"), getL(f, "n_rows_pm_ne_npnm"),
                getL(f, "n_legs_bothcharge"));
    std::printf("Mixed-event kind: %s (%s)\n", mix_kind,
                std::string(mix_kind) == "mixnc" ? "pools per trigger only, no multiplicity classes"
                                                 : "pools in classes trigger x n+ x n-");
    std::printf("Normalisation window for mixed +-: [%.2f, %.2f]\n", norm_lo, norm_hi);
    const std::string mk = std::string("h_") + mix_kind + "_";

    const double win_lo[3] = {0.0, kLowHi, kHighLo};
    const double win_hi[3] = {kLowHi, kHighLo, kXmax};
    const char*  win_nm[3] = {"LOW ", "MID ", "HIGH"};

    Set S[2][2];   // [scheme][trigger]
    for (int s = 0; s < 2; ++s)
        for (int t = 0; t < 2; ++t) {
            Set& X = S[s][t];
            const std::string suf = std::string(s == 1 ? "_best" : "") + "_" + trigName(t);
            X.se_pm  = getH(f, "h_se_pm" + suf);   X.se_pp  = getH(f, "h_se_pp" + suf);   X.se_mm  = getH(f, "h_se_mm" + suf);
            X.mix_pm = getH(f, mk + "pm" + suf);   X.mix_pp = getH(f, mk + "pp" + suf);   X.mix_mm = getH(f, mk + "mm" + suf);
            X.se_pmF  = getH(f, "h_se_pm" + suf + "_fine");  X.se_ppF  = getH(f, "h_se_pp" + suf + "_fine");  X.se_mmF  = getH(f, "h_se_mm" + suf + "_fine");
            X.mix_pmF = getH(f, mk + "pm" + suf + "_fine"); X.mix_ppF = getH(f, mk + "pp" + suf + "_fine"); X.mix_mmF = getH(f, mk + "mm" + suf + "_fine");
            if (!X.se_pm || !X.se_pp || !X.se_mm || !X.mix_pm || !X.mix_pp || !X.mix_mm ||
                !X.se_pmF || !X.se_ppF || !X.se_mmF || !X.mix_pmF || !X.mix_ppF || !X.mix_mmF) return;
            const std::string id = std::string("_") + schemeName(s) + "_" + trigName(t);
            X.k  = makeK(X.mix_pm,  X.mix_pp,  X.mix_mm,  "k" + id);
            X.kF = makeK(X.mix_pmF, X.mix_ppF, X.mix_mmF, "kF" + id);
            for (int w = 0; w < 3; ++w) kWindow(X.mix_pm, X.mix_pp, X.mix_mm, win_lo[w], win_hi[w], X.k_win[w], X.ek_win[w]);
            X.cb_ls = makeCB(X.se_pp, X.se_mm, "cb_ls" + id);
            // alpha: normalise N_mix+- to CB_ls in [norm_lo, norm_hi]
            double A, eA, B, eB;
            integ(X.cb_ls, norm_lo, norm_hi, A, eA);
            integ(X.mix_pm, norm_lo, norm_hi, B, eB);
            X.alpha  = (B > 0) ? A / B : 0.0;
            X.ealpha = (A > 0 && B > 0) ? X.alpha * std::sqrt(std::pow(eA / A, 2) + std::pow(eB / B, 2)) : 0.0;
            X.cb_mix = scaled(X.mix_pm, X.alpha, X.ealpha, "cb_mix" + id);
            X.cb_k   = product(X.k, X.cb_ls, "cb_k" + id);
            X.s_ls   = makeSignal(X.se_pm, X.cb_ls,  "s_ls" + id);
            X.s_k    = makeSignal(X.se_pm, X.cb_k,   "s_k" + id);
            X.s_mix  = makeSignal(X.se_pm, X.cb_mix, "s_mix" + id);
            X.shape_pm = makeRatio(X.cb_ls, X.cb_mix, 1.0, "shape_pm" + id);
            // like-sign individually: N_se++ / (beta N_mix++), beta from the same window
            double a2, ea2, b2, eb2;
            integ(X.se_pp, norm_lo, norm_hi, a2, ea2); integ(X.mix_pp, norm_lo, norm_hi, b2, eb2);
            TH1D* mpp = scaled(X.mix_pp, b2 > 0 ? a2 / b2 : 0.0, 0.0, "mixpp_n" + id);
            X.shape_pp = makeRatio(X.se_pp, mpp, 1.0, "shape_pp" + id);
            integ(X.se_mm, norm_lo, norm_hi, a2, ea2); integ(X.mix_mm, norm_lo, norm_hi, b2, eb2);
            TH1D* mmm = scaled(X.mix_mm, b2 > 0 ? a2 / b2 : 0.0, 0.0, "mixmm_n" + id);
            X.shape_mm = makeRatio(X.se_mm, mmm, 1.0, "shape_mm" + id);
        }

    // ------------------------------------------------------------ tables
    for (int s = 0; s < 2; ++s) {
        std::printf("\n################ SCHEME: %s ################\n", schemeName(s));
        for (int t = 0; t < 2; ++t) {
            Set& X = S[s][t];
            std::printf("\n--- %s %s: window integrals\n", schemeName(s), trigName(t));
            std::printf("  mixed pairs: +-=%.0f ++=%.0f --=%.0f  alpha(norm [%.2f,%.2f]) = %.5g +- %.2g\n",
                        X.mix_pm->Integral(), X.mix_pp->Integral(), X.mix_mm->Integral(), norm_lo, norm_hi, X.alpha, X.ealpha);
            std::printf("  %-5s %12s %10s %10s %12s %12s | %8s %12s %12s | %12s %12s %10s\n",
                        "win", "N+-", "N++", "N--", "CB_ls", "S_ls", "k_win", "S_k(bin k)", "CB_k",
                        "CB_mix", "S_mix", "CBls/CBmix");
            for (int w = 0; w < 3; ++w) {
                double npm, e1, npp, e2, nmm, e3, cb, ecb, sls, esls, sk, esk, cbk, ecbk, cbm, ecbm, smx, esmx;
                integ(X.se_pm, win_lo[w], win_hi[w], npm, e1);
                integ(X.se_pp, win_lo[w], win_hi[w], npp, e2);
                integ(X.se_mm, win_lo[w], win_hi[w], nmm, e3);
                integ(X.cb_ls, win_lo[w], win_hi[w], cb, ecb);
                integ(X.s_ls, win_lo[w], win_hi[w], sls, esls);
                integ(X.s_k, win_lo[w], win_hi[w], sk, esk);
                integ(X.cb_k, win_lo[w], win_hi[w], cbk, ecbk);
                integ(X.cb_mix, win_lo[w], win_hi[w], cbm, ecbm);
                integ(X.s_mix, win_lo[w], win_hi[w], smx, esmx);
                std::printf("  %-5s %12.0f %10.0f %10.0f %12.1f %12.1f | %8.4f %12.1f %12.1f | %12.1f %12.1f %10.4f\n",
                            win_nm[w], npm, npp, nmm, cb, sls, X.k_win[w], sk, cbk, cbm, smx, cbm > 0 ? cb / cbm : 0.0);
                std::printf("  %-5s %12s %10s %10s %12.1f %12.1f | %8.4f %12.1f %12.1f | %12.1f %12.1f\n",
                            "+-err", "", "", "", ecb, esls, X.ek_win[w], esk, ecbk, ecbm, esmx);
            }
            if (s == 0 && t == 0) {
                double sls, esls;
                integ(X.s_ls, kLowHi, kHighLo, sls, esls);
                std::printf("  CHECK (full data only): S_ls PT3 MID = %.1f  (recipe: 768003; diff %.1f)\n", sls, sls - 768003.0);
            }
        }

        // ratios 63 S_PT2 / S_PT3
        Set& P3 = S[s][0];
        Set& P2 = S[s][1];
        TH1D* r_ls  = makeRatio(P2.s_ls,  P3.s_ls,  kTrigCorr, std::string("r_ls_")  + schemeName(s));
        TH1D* r_k   = makeRatio(P2.s_k,   P3.s_k,   kTrigCorr, std::string("r_k_")   + schemeName(s));
        TH1D* r_mix = makeRatio(P2.s_mix, P3.s_mix, kTrigCorr, std::string("r_mix_") + schemeName(s));
        TH1D* rr[3] = {r_ls, r_k, r_mix};
        const char* rn[3] = {"S_ls ", "S_k  ", "S_mix"};
        double C[3], eC[3], chi2[3]; int ndf[3];
        std::printf("\n--- %s: 63 * S_PT2 / S_PT3\n", schemeName(s));
        std::printf("  %-6s %18s %12s | %18s %18s %18s\n", "var", "C_MID pol0 [0.14,0.80]", "chi2/ndf", "<r> LOW", "<r> MID", "<r> HIGH");
        for (int v = 0; v < 3; ++v) {
            fitMid(rr[v], C[v], eC[v], chi2[v], ndf[v]);
            double m[3], e[3]; int nb[3];
            for (int w = 0; w < 3; ++w) regionMean(rr[v], win_lo[w], win_hi[w], m[w], e[w], nb[w]);
            std::printf("  %-6s %9.4f +- %6.4f %6.1f/%-3d  | %8.4f +- %6.4f %8.4f +- %6.4f %8.4f +- %6.4f\n",
                        rn[v], C[v], eC[v], chi2[v], ndf[v], m[0], e[0], m[1], e[1], m[2], e[2]);
        }
        if (s == 0)
            std::printf("  CHECK (full data only): C_MID(S_ls) = %.4f  (recipe: 1.607; diff %+.4f)\n", C[0], C[0] - 1.607);
        std::printf("  k_win %s: PT3 LOW %.4f MID %.4f HIGH %.4f | PT2 LOW %.4f MID %.4f HIGH %.4f\n", schemeName(s),
                    P3.k_win[0], P3.k_win[1], P3.k_win[2], P2.k_win[0], P2.k_win[1], P2.k_win[2]);

        // per-bin table
        std::printf("\n--- %s: per-bin table (58-bin grid)\n", schemeName(s));
        std::printf("  %3s %6s %6s | %7s %7s | %7s %7s | %8s %8s %8s | %8s %8s %8s\n",
                    "bin", "lo", "hi", "k_PT3", "k_PT2", "shp_PT3", "shp_PT2", "S_ls3", "S_k3", "S_mix3", "r_ls", "r_k", "r_mix");
        for (int b = 1; b <= r_ls->GetNbinsX(); ++b) {
            std::printf("  %3d %6.4f %6.4f | %7.4f %7.4f | %7.4f %7.4f | %8.0f %8.0f %8.0f | %8.4f %8.4f %8.4f\n",
                        b, r_ls->GetXaxis()->GetBinLowEdge(b), r_ls->GetXaxis()->GetBinUpEdge(b),
                        P3.k->GetBinContent(b), P2.k->GetBinContent(b),
                        P3.shape_pm->GetBinContent(b), P2.shape_pm->GetBinContent(b),
                        P3.s_ls->GetBinContent(b), P3.s_k->GetBinContent(b), P3.s_mix->GetBinContent(b),
                        r_ls->GetBinContent(b), r_k->GetBinContent(b), r_mix->GetBinContent(b));
        }

        // -------------------------------------------------------- plots
        const std::string base = std::string(out_dir) + "/mixing_" + tag + "_" + mix_kind + "_";
        // (e1) k(m), fine grid, PT3 and PT2
        {
            TCanvas* c = new TCanvas(Form("c_k_%s", schemeName(s)), "k-factor", 1000, 650);
            c->SetMargin(0.11, 0.04, 0.13, 0.08);
            TH1D* k3 = P3.kF; TH1D* k2 = P2.kF;
            style(k3, kBlack, 20, 0.6); style(k2, kRed, 24, 0.6);
            k3->SetTitle(Form("k(m) = N_{mix}^{+-} / 2#sqrt{N_{mix}^{++} N_{mix}^{--}}  (%s, %s);M_{ee} [GeV/c^{2}];k", schemeName(s), mix_kind));
            double kmax = 1.4;
            for (int w = 0; w < 3; ++w) kmax = std::max(kmax, 1.2 * std::max(P3.k_win[w], P2.k_win[w]));
            k3->GetYaxis()->SetRangeUser(kmax > 1.4 ? 0.0 : 0.6, kmax);
            k3->GetXaxis()->SetRangeUser(0.0, kXmax);
            k3->Draw("E1");
            k2->Draw("E1 SAME");
            TLine l1(0.0, 1.0, kXmax, 1.0); l1.SetLineStyle(2); l1.SetLineColor(kGray + 2); l1.DrawClone();
            for (int w = 0; w < 3; ++w) {
                TLine* lw = new TLine(win_lo[w], P3.k_win[w], win_hi[w], P3.k_win[w]);
                lw->SetLineColor(kBlue); lw->SetLineWidth(2); lw->Draw();
                TLine* lw2 = new TLine(win_lo[w], P2.k_win[w], win_hi[w], P2.k_win[w]);
                lw2->SetLineColor(kMagenta); lw2->SetLineWidth(2); lw2->SetLineStyle(2); lw2->Draw();
            }
            TLegend* leg = new TLegend(0.14, 0.70, 0.60, 0.90);
            leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.028);
            leg->AddEntry(k3, "PT3, 5 MeV bins", "lpe");
            leg->AddEntry(k2, "PT2, 5 MeV bins", "lpe");
            leg->AddEntry((TObject*)nullptr, Form("PT3 windows: LOW %.3f  MID %.3f  HIGH %.3f", P3.k_win[0], P3.k_win[1], P3.k_win[2]), "");
            leg->AddEntry((TObject*)nullptr, Form("PT2 windows: LOW %.3f  MID %.3f  HIGH %.3f", P2.k_win[0], P2.k_win[1], P2.k_win[2]), "");
            leg->Draw();
            c->SaveAs((base + "k_" + schemeName(s) + ".pdf").c_str());
            c->SaveAs((base + "k_" + schemeName(s) + ".png").c_str());
        }
        // (e2) same-event CB_ls vs normalised mixed +-, and ratios, PT3 / PT2
        {
            TCanvas* c = new TCanvas(Form("c_shape_%s", schemeName(s)), "same-event vs mixed", 1400, 900);
            c->Divide(2, 2);
            for (int t = 0; t < 2; ++t) {
                Set& X = S[s][t];
                c->cd(1 + t); gPad->SetLogy(); gPad->SetMargin(0.12, 0.04, 0.13, 0.08);
                TH1D* a = static_cast<TH1D*>(X.cb_ls->Clone()); a->SetDirectory(nullptr);
                TH1D* b = static_cast<TH1D*>(X.cb_mix->Clone()); b->SetDirectory(nullptr);
                TH1D* d = static_cast<TH1D*>(X.se_pm->Clone()); d->SetDirectory(nullptr);
                style(a, kRed, 20); style(b, kBlue, 24); style(d, kBlack, 21);
                a->SetTitle(Form("%s %s: 2#sqrt{N^{++}N^{--}} (same event) vs #alpha N_{mix}^{+-};M_{ee} [GeV/c^{2}];pairs / bin",
                                 schemeName(s), trigName(t)));
                a->GetXaxis()->SetRangeUser(0.0, kXmax);
                const double ymax = std::max(d->GetMaximum(), a->GetMaximum()) * 3.0;
                a->GetYaxis()->SetRangeUser(0.5, ymax);
                a->Draw("E1"); b->Draw("E1 SAME"); d->Draw("E1 SAME");
                TLegend* leg = new TLegend(0.50, 0.70, 0.95, 0.90);
                leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.03);
                leg->AddEntry(d, "N^{+-} same event", "lpe");
                leg->AddEntry(a, "2#sqrt{N^{++}N^{--}} same event", "lpe");
                leg->AddEntry(b, Form("#alpha N_{mix}^{+-}, norm [%.2f,%.2f]", norm_lo, norm_hi), "lpe");
                leg->Draw();
                c->cd(3 + t); gPad->SetMargin(0.12, 0.04, 0.13, 0.08);
                TH1D* r = static_cast<TH1D*>(X.shape_pm->Clone()); r->SetDirectory(nullptr);
                TH1D* rp = static_cast<TH1D*>(X.shape_pp->Clone()); rp->SetDirectory(nullptr);
                TH1D* rm = static_cast<TH1D*>(X.shape_mm->Clone()); rm->SetDirectory(nullptr);
                style(r, kBlack, 20); style(rp, kRed, 24); style(rm, kBlue, 25);
                r->SetTitle(Form("%s %s: same-event / normalised mixed;M_{ee} [GeV/c^{2}];ratio", schemeName(s), trigName(t)));
                r->GetYaxis()->SetRangeUser(0.0, 3.0);
                r->GetXaxis()->SetRangeUser(0.0, kXmax);
                r->Draw("E1"); rp->Draw("E1 SAME"); rm->Draw("E1 SAME");
                TLine* l1 = new TLine(0.0, 1.0, kXmax, 1.0); l1->SetLineStyle(2); l1->SetLineColor(kGray + 2); l1->Draw();
                TBox* bx = new TBox(norm_lo, 0.0, norm_hi, 3.0); bx->SetFillColorAlpha(kGreen + 2, 0.10); bx->Draw();
                TLegend* leg2 = new TLegend(0.50, 0.72, 0.95, 0.90);
                leg2->SetBorderSize(0); leg2->SetFillStyle(0); leg2->SetTextSize(0.03);
                leg2->AddEntry(r, "2#sqrt{N^{++}N^{--}} / #alpha N_{mix}^{+-}", "lpe");
                leg2->AddEntry(rp, "N^{++} / norm. N_{mix}^{++}", "lpe");
                leg2->AddEntry(rm, "N^{--} / norm. N_{mix}^{--}", "lpe");
                leg2->Draw();
            }
            c->SaveAs((base + "shape_" + schemeName(s) + ".pdf").c_str());
            c->SaveAs((base + "shape_" + schemeName(s) + ".png").c_str());
        }
        // (e3) 63 S_PT2 / S_PT3 in three variants
        {
            TCanvas* c = new TCanvas(Form("c_ratio_%s", schemeName(s)), "trigger ratio", 1100, 700);
            c->SetMargin(0.11, 0.04, 0.13, 0.09);
            style(r_ls, kBlack, 20); style(r_k, kBlue, 24); style(r_mix, kRed, 25);
            r_ls->SetTitle(Form("63 S_{PT2} / S_{PT3} (%s, %s): like-sign, k-corrected, mixed;M_{ee} [GeV/c^{2}];63 #upoint S_{PT2} / S_{PT3}", schemeName(s), mix_kind));
            r_ls->GetXaxis()->SetRangeUser(0.0, kXmax);
            r_ls->GetYaxis()->SetRangeUser(0.5, 4.0);
            r_ls->Draw("E1"); r_k->Draw("E1 SAME"); r_mix->Draw("E1 SAME");
            Color_t cc[3] = {kBlack, kBlue, kRed};
            for (int v = 0; v < 3; ++v) {
                TLine* lf = new TLine(kLowHi, C[v], kHighLo, C[v]);
                lf->SetLineColor(cc[v]); lf->SetLineWidth(2); lf->SetLineStyle(v == 0 ? 1 : 2); lf->Draw();
            }
            TLine* l1 = new TLine(0.0, 1.0, kXmax, 1.0); l1->SetLineStyle(2); l1->SetLineColor(kGray + 2); l1->Draw();
            TLegend* leg = new TLegend(0.40, 0.68, 0.96, 0.90);
            leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.026);
            leg->AddEntry(r_ls,  Form("S_{ls} = N^{+-} - 2#sqrt{N^{++}N^{--}}:  C_{MID} = %.3f #pm %.3f, #chi^{2}/ndf = %.1f/%d", C[0], eC[0], chi2[0], ndf[0]), "lpe");
            leg->AddEntry(r_k,   Form("S_{k} = N^{+-} - k(m) 2#sqrt{N^{++}N^{--}}:  C_{MID} = %.3f #pm %.3f, #chi^{2}/ndf = %.1f/%d", C[1], eC[1], chi2[1], ndf[1]), "lpe");
            leg->AddEntry(r_mix, Form("S_{mix} = N^{+-} - #alpha N_{mix}^{+-}:  C_{MID} = %.3f #pm %.3f, #chi^{2}/ndf = %.1f/%d", C[2], eC[2], chi2[2], ndf[2]), "lpe");
            leg->Draw();
            c->SaveAs((base + "ratio_" + schemeName(s) + ".pdf").c_str());
            c->SaveAs((base + "ratio_" + schemeName(s) + ".png").c_str());
        }
        std::printf("  plots: %sk_%s, %sshape_%s, %sratio_%s .{pdf,png}\n", base.c_str(), schemeName(s),
                    base.c_str(), schemeName(s), base.c_str(), schemeName(s));
    }
    f->Close();
}
