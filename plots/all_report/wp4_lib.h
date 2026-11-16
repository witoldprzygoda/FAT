// ============================================================================
// wp4_lib.h  --  model, fits and estimators of the WP4 stress closure (prefix wp4)
// ============================================================================
// SPEC report_all, WP4. This header is the implementation that WP4 validates:
//   * Smooth    : exp(pol3) smoothing of the non-eta templates (D13), effective-count
//                 Poisson likelihood with one common scale sum(var)/sum(content) per histogram,
//                 bin-integrated (8-point Gauss-Legendre), one fit per
//                 F10 histogram (slice x E1/E2) and per F5g slice. Sparse histograms:
//                 N_eff 20-50 exp(pol2); N_eff < 20 the shape of the parent histogram
//                 (F5g -> F10 E2 of the slice; omega/rho0 split -> its total) scaled to the
//                 own content (keeps the kinematic thresholds of S4/S5)
//   * CB        : 2 sqrt(N++ N--) in 10 MeV cells, interpolated to the 2.5 MeV fine
//                 bins, spread by f(cell | 50 MeV band) from N++ + N-- (= WP0 rule,
//                 reproduced to 1e-9 from the like-sign G3 counts)
//   * NLFit     : Poisson likelihood fit of a model that is non-linear in its
//                 parameters (products c_j b, c_j kappa_B, a_eta g), by Fisher scoring:
//                 every step is one call of bla::LinFit (s5_blind_ana_fit.h) on the
//                 linearised model with the variance V = max(mu, 1) + k_e^2 var(CB)
//                 + sum_t coef_t^2 var(T_t) (Barlow-Beeston-lite on the raw signal
//                 templates); the fixed point is the IRLS / Poisson ML solution of the
//                 fitter. Gaussian constraints g = 1 +- 0.30, k_E2 = 1 +- 0.05.
//   * F-main    : mu = k_e CB + c_j (1 + b x)[G_D + kappa_B G_C] + a [C2 + kappa (C3 + C4)
//                 + rho C5 + g G_eta], F10 (5 slices x E1/E2 x 70 bins), internally with
//                 the linear amplitudes beta2 = a, beta34 = a kappa, beta5 = a rho
//                 (identical model; kappa and rho by the delta method)
//   * F-free    : alpha2 C2 + alpha3 C3 + alpha4 C4 + alpha5 C5 + alpha2 g G_eta
//   * stage A   : F10 restricted to CR-U (0.565-0.700, all slices) and CR-K (S5,
//                 0.400-0.5225), background free, eta terms fixed to F-main
//   * F-gold    : F5g (M(ee) >= 20 MeV), alpha2, alpha5, g, c_j, b free, kappa_B and
//                 a kappa (C3 + C4) fixed to F-main
//   * B4        : Y = (N - CB)_B - [(N - CB)_U - eta_U] Bk_A(B)/Bk_A(U) - eta_other(B)
//   * peak SB   : Y = [(N - CB)_P - R' (N - CB)_SB] / (1 - R' f_leak),
//                 R' = [Bk_A + bump tail + G_eta](P) / [same](SB), f_leak = (C4 + C5)(SB) / (C4 + C5)(P)
// All counts are RAW data rows at L = 6079.2 nb^-1 (templates: nb x 2 (half) x kNbToRaw);
// yields Y are raw rows x 64/63.
// ============================================================================
#ifndef WP4_LIB_H
#define WP4_LIB_H

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <array>
#include <vector>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TRandom3.h"
#include "TVectorD.h"

#include "s5_blind_ana_fit.h"
#include "hists_io.h"
#include "wp4_grid.h"

namespace wp4 {

typedef std::vector<double> Vec;
const double kHalfToRaw = 2.0 * kNbToRaw;   // half-sample template [nb] -> expected raw rows at L_data
const double kYf = kPT3;                    // yields: raw rows x 64/63

// ===========================================================================
// 1. smoothing exp(pol n) in x = (M4 - 0.55)/0.1
// ===========================================================================
const double kGLx[8] = {-0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
                        0.1834346424956498,  0.5255324099163290,  0.7966664774136267,  0.9602898564975363};
const double kGLw[8] = {0.1012285362903763, 0.2223810344533745, 0.3137066860094368, 0.3626837833783620,
                        0.3626837833783620, 0.3137066860094368, 0.2223810344533745, 0.1012285362903763};
struct Smooth {
    int order = -1;                 // -1 = zero function
    double p[4] = {0, 0, 0, 0};
    double neff = 0;
    bool conv = false;
    bool scaled = false;            // shape of the parent histogram, scaled
    inline double Dens(double m) const {
        if (order < 0) return 0.;
        const double x = TiltX(m);
        return std::exp(p[0] + x * (p[1] + x * (p[2] + x * p[3])));
    }
    // integral over [lo, hi] of f, and of x f (for the exact tilt of the truth)
    void Int2(double lo, double hi, double& i0, double& i1) const {
        i0 = i1 = 0.;
        if (order < 0 || !(hi > lo)) return;
        const double h = 0.5 * (hi - lo), c = 0.5 * (hi + lo);
        for (int q = 0; q < 8; ++q) {
            const double m = c + h * kGLx[q];
            const double f = Dens(m) * kGLw[q] * h;
            i0 += f; i1 += f * TiltX(m);
        }
    }
    double Int(double lo, double hi) const { double a, b; Int2(lo, hi, a, b); return a; }
};
// fit to the bin contents c (bins 1..n of edges e) with event-group variances v
inline Smooth SmoothFit(const Vec& c, const Vec& v, const std::vector<double>& e, const Smooth* parent = nullptr, int forceOrder = -1) {
    Smooth S;
    const int n = (int)e.size() - 1;
    double C = 0, V = 0;
    for (int i = 1; i <= n; ++i) { C += c[i]; V += v[i]; }
    if (!(C > 0)) return S;
    S.neff = V > 0 ? C * C / V : 0;
    if (S.neff < 20 && parent && parent->order >= 0 && forceOrder < 0) {
        const double neff = S.neff;
        S = *parent;
        S.p[0] += std::log(C / parent->Int(e[0], e[n]));
        S.neff = neff; S.scaled = true;
        return S;
    }
    int ord = S.neff >= 50 ? 3 : S.neff >= 20 ? 2 : S.neff >= 8 ? 1 : 0;
    if (forceOrder >= 0) ord = forceOrder;
    S.order = ord;
    const int np = ord + 1;
    const double sbar = V > 0 ? V / C : 1.0;
    Vec s(n + 1), nn(n + 1);
    // one common scale per histogram (scaled Poisson): per-bin scales v_i/c_i correlate the weight with the
    // fluctuation (a bin hit by a heavy event group gets a small weight) and bias the normalisation low
    for (int i = 1; i <= n; ++i) { s[i] = sbar; nn[i] = c[i] / s[i]; }
    for (int k = 0; k < 4; ++k) S.p[k] = 0;
    S.p[0] = std::log(C / (e[n] - e[0]));
    // precompute nodes
    std::vector<std::vector<double>> xm(n + 1, Vec(8)), wq(n + 1, Vec(8));
    for (int i = 1; i <= n; ++i) {
        const double h = 0.5 * (e[i] - e[i - 1]), cc = 0.5 * (e[i] + e[i - 1]);
        for (int q = 0; q < 8; ++q) { xm[i][q] = TiltX(cc + h * kGLx[q]); wq[i][q] = kGLw[q] * h; }
    }
    auto evalL = [&](const double* p, Vec* F, std::vector<Vec>* dF) {
        double L = 0;
        for (int i = 1; i <= n; ++i) {
            double f = 0; double d[4] = {0, 0, 0, 0};
            for (int q = 0; q < 8; ++q) {
                const double x = xm[i][q];
                const double ex = wq[i][q] * std::exp(p[0] + x * (p[1] + x * (p[2] + x * p[3])));
                f += ex; double xk = 1; for (int k = 0; k < np; ++k) { d[k] += xk * ex; xk *= x; }
            }
            const double m = std::max(f / s[i], 1e-300);
            L += (nn[i] > 0 ? nn[i] * std::log(m) : 0.) - m;
            if (F) (*F)[i] = f;
            if (dF) for (int k = 0; k < np; ++k) (*dF)[i][k] = d[k];
        }
        return L;
    };
    Vec F(n + 1); std::vector<Vec> dF(n + 1, Vec(4));
    double L = evalL(S.p, &F, &dF);
    for (int it = 0; it < 200; ++it) {
        bla::Mat I = bla::Zero(np, np); Vec g(np, 0.);
        for (int i = 1; i <= n; ++i) {
            const double m = std::max(F[i] / s[i], 1e-300);
            for (int k = 0; k < np; ++k) {
                const double dk = dF[i][k] / s[i];
                g[k] += (nn[i] / m - 1.0) * dk;
                for (int l = 0; l < np; ++l) I[k][l] += dk * dF[i][l] / s[i] / m;
            }
        }
        bla::Mat Ii = I;
        if (!bla::CholInv(Ii)) break;
        Vec d(np, 0.); double q = 0;
        for (int k = 0; k < np; ++k) for (int l = 0; l < np; ++l) d[k] += Ii[k][l] * g[l];
        for (int k = 0; k < np; ++k) q += d[k] * g[k];
        if (q < 1e-12) { S.conv = true; break; }
        double t = 1.0, pn[4];
        bool acc = false;
        for (int h = 0; h < 30; ++h) {
            for (int k = 0; k < 4; ++k) pn[k] = S.p[k] + (k < np ? t * d[k] : 0.);
            const double Ln = evalL(pn, nullptr, nullptr);
            if (Ln >= L - 1e-12) { acc = true; break; }
            t *= 0.5;
        }
        if (!acc) { S.conv = true; break; }
        for (int k = 0; k < 4; ++k) S.p[k] = pn[k];
        L = evalL(S.p, &F, &dF);
    }
    return S;
}

// ===========================================================================
// 2. templates of one half and one mass variant (raw rows at L_data)
// ===========================================================================
struct Tmpl {
    Vec R[kNTComp][kNh15];                     // refined cells, index 0..nR+1
    Vec F15[kNTComp][kNh15], F15v[kNTComp][kNh15];
    Vec F10[kNTComp][10], F10v[kNTComp][10];
    double tot[kNTComp];
};
inline bool LoadTmpl(TFile* f, int mv, int half, const Grid& G, Tmpl& T) {
    const std::string dir = Form("%s_h%s/", MVarTag(mv), half == 0 ? "A" : "B");
    for (int c = 0; c < kNTComp; ++c) {
        for (int h15 = 0; h15 < kNh15; ++h15) {
            const int j = h15 / kNe3, e3 = h15 % kNe3;
            TH1* r = (TH1*)f->Get((dir + Form("%s_R_S%de%d", TCompTag(c), j + 1, e3)).c_str());
            TH1* x = (TH1*)f->Get((dir + Form("%s_F15_S%de%d", TCompTag(c), j + 1, e3)).c_str());
            if (!r || !x) { printf("wp4: template %s missing\n", dir.c_str()); return false; }
            T.R[c][h15].assign(G.nR + 2, 0.); T.F15[c][h15].assign(G.nF + 2, 0.); T.F15v[c][h15].assign(G.nF + 2, 0.);
            for (int b = 0; b <= G.nR + 1; ++b) T.R[c][h15][b] = r->GetBinContent(b) * kHalfToRaw;
            for (int b = 0; b <= G.nF + 1; ++b) {
                T.F15[c][h15][b] = x->GetBinContent(b) * kHalfToRaw;
                T.F15v[c][h15][b] = std::pow(x->GetBinError(b) * kHalfToRaw, 2);
            }
        }
        for (int k = 0; k < 10; ++k) {
            TH1* x = (TH1*)f->Get((dir + Form("%s_F10_S%dE%d", TCompTag(c), k / 2 + 1, k % 2 + 1)).c_str());
            if (!x) return false;
            T.F10[c][k].assign(G.nF + 2, 0.); T.F10v[c][k].assign(G.nF + 2, 0.);
            for (int b = 0; b <= G.nF + 1; ++b) {
                T.F10[c][k][b] = x->GetBinContent(b) * kHalfToRaw;
                T.F10v[c][k][b] = std::pow(x->GetBinError(b) * kHalfToRaw, 2);
            }
        }
        TH1* t = (TH1*)f->Get((dir + Form("%s_tot", TCompTag(c))).c_str());
        T.tot[c] = t ? t->GetBinContent(1) * kHalfToRaw : 0.;
    }
    return true;
}
// full sample = half A + half B (contents and event-group variances add; nb x kNbToRaw = (A + B) x kHalfToRaw / 2)
inline bool LoadTmplFull(TFile* f, int mv, const Grid& G, Tmpl& T) {
    Tmpl* B = new Tmpl();
    if (!LoadTmpl(f, mv, 0, G, T) || !LoadTmpl(f, mv, 1, G, *B)) { delete B; return false; }
    for (int c = 0; c < kNTComp; ++c) {
        for (int h = 0; h < kNh15; ++h) {
            for (size_t b = 0; b < T.R[c][h].size(); ++b) T.R[c][h][b] = 0.5 * (T.R[c][h][b] + B->R[c][h][b]);
            for (size_t b = 0; b < T.F15[c][h].size(); ++b) {
                T.F15[c][h][b] = 0.5 * (T.F15[c][h][b] + B->F15[c][h][b]); T.F15v[c][h][b] = 0.25 * (T.F15v[c][h][b] + B->F15v[c][h][b]);
            }
        }
        for (int k = 0; k < 10; ++k)
            for (size_t b = 0; b < T.F10[c][k].size(); ++b) {
                T.F10[c][k][b] = 0.5 * (T.F10[c][k][b] + B->F10[c][k][b]); T.F10v[c][k][b] = 0.25 * (T.F10v[c][k][b] + B->F10v[c][k][b]);
            }
        T.tot[c] = 0.5 * (T.tot[c] + B->tot[c]);
    }
    delete B;
    return true;
}
// window sum of a component over refined cells (all 15 cells, or only slice js)
inline double WinSum(const Tmpl& T, const Grid& G, int c, double lo, double hi, int js = -1) {
    double s = 0;
    for (int r : G.RefIn(lo, hi))
        for (int h15 = 0; h15 < kNh15; ++h15) { if (js >= 0 && h15 / kNe3 != js) continue; s += T.R[c][h15][r]; }
    return s;
}

// smoothed non-eta templates used by the fit (one per F10 histogram and per F5g slice)
struct BkgSmooth {
    Smooth GD10[10], GC10[10], GD5[5], GC5[5];
};
inline void MakeFitSmooth(const Tmpl& T, const Grid& G, BkgSmooth& B) {
    for (int k = 0; k < 10; ++k) {
        Vec c(G.nF + 2, 0.), v(G.nF + 2, 0.), cc(G.nF + 2, 0.), vc(G.nF + 2, 0.);
        for (int b = 1; b <= G.nF; ++b) {
            c[b] = T.F10[kTGDom][k][b] + T.F10[kTGDnom][k][b]; v[b] = T.F10v[kTGDom][k][b] + T.F10v[kTGDnom][k][b];
            cc[b] = T.F10[kTGCom][k][b] + T.F10[kTGCnom][k][b]; vc[b] = T.F10v[kTGCom][k][b] + T.F10v[kTGCnom][k][b];
        }
        B.GD10[k] = SmoothFit(c, v, G.f10); B.GC10[k] = SmoothFit(cc, vc, G.f10);
    }
    for (int j = 0; j < kNj; ++j) {
        const int h = H15(j, 2);
        Vec c(G.nF + 2, 0.), v(G.nF + 2, 0.), cc(G.nF + 2, 0.), vc(G.nF + 2, 0.);
        for (int b = 1; b <= G.nF; ++b) {
            c[b] = T.F15[kTGDom][h][b] + T.F15[kTGDnom][h][b]; v[b] = T.F15v[kTGDom][h][b] + T.F15v[kTGDnom][h][b];
            cc[b] = T.F15[kTGCom][h][b] + T.F15[kTGCnom][h][b]; vc[b] = T.F15v[kTGCom][h][b] + T.F15v[kTGCnom][h][b];
        }
        B.GD5[j] = SmoothFit(c, v, G.f10, &B.GD10[2 * j + 1]); B.GC5[j] = SmoothFit(cc, vc, G.f10, &B.GC10[2 * j + 1]);
    }
}
// smoothed truth: the same smoothing as the fit (per F10 histogram and per F5g slice) on the truth half,
// plus the omega/rho0 split per F10 histogram (scenario (e)). Toy cell (j, e3): e3 = 0 -> F10 E1; e3 = 2 -> F5g;
// e3 = 1 -> F10 E2 - F5g (clipped at 0).
struct TruthSmooth {
    BkgSmooth B;
    Smooth GDom[10], GDnom[10], GCom[10], GCnom[10];
};
inline void MakeTruthSmooth(const Tmpl& T, const Grid& G, TruthSmooth& S) {
    MakeFitSmooth(T, G, S.B);
    for (int k = 0; k < 10; ++k) {
        auto one = [&](int comp, const Smooth* par) {
            Vec c(G.nF + 2, 0.), v(G.nF + 2, 0.);
            for (int b = 1; b <= G.nF; ++b) { c[b] = T.F10[comp][k][b]; v[b] = T.F10v[comp][k][b]; }
            return SmoothFit(c, v, G.f10, par);
        };
        S.GDom[k] = one(kTGDom, &S.B.GD10[k]); S.GDnom[k] = one(kTGDnom, &S.B.GD10[k]);
        S.GCom[k] = one(kTGCom, &S.B.GC10[k]); S.GCnom[k] = one(kTGCnom, &S.B.GC10[k]);
    }
}
// integral of a truth background over [lo, hi) of cell (j, e3): i0 = int f, i1 = int x f
// mode 0: G_D + G_C (each x 1); mode 1: 0.3 x omega/rho0 part + 1.6 x the rest, split into e3 like mode 0
inline void TruthBkgCell(const TruthSmooth& S, int j, int e3, double lo, double hi, int mode, double& i0, double& i1) {
    auto pair = [&](const Smooth& a, const Smooth& b, double& x0, double& x1) {
        double a0, a1, b0, b1; a.Int2(lo, hi, a0, a1); b.Int2(lo, hi, b0, b1); x0 = a0 + b0; x1 = a1 + b1;
    };
    double s0, s1;   // standard (mode 0) content of the cell
    const int kE1 = 2 * j, kE2 = 2 * j + 1;
    if (e3 == 0) pair(S.B.GD10[kE1], S.B.GC10[kE1], s0, s1);
    else if (e3 == 2) pair(S.B.GD5[j], S.B.GC5[j], s0, s1);
    else {
        double e0, e1, f0, f1;
        pair(S.B.GD10[kE2], S.B.GC10[kE2], e0, e1); pair(S.B.GD5[j], S.B.GC5[j], f0, f1);
        s0 = e0 - f0; s1 = e1 - f1;
        if (!(s0 > 0)) { s0 = 0; s1 = 0; }
    }
    if (mode == 0) { i0 = s0; i1 = s1; return; }
    const int k = e3 == 0 ? kE1 : kE2;
    double p0, p1, o0, o1, q0, q1, r0, r1, t0, t1;
    pair(S.B.GD10[k], S.B.GC10[k], p0, p1);                // standard F10 histogram content
    S.GDom[k].Int2(lo, hi, o0, o1); S.GCom[k].Int2(lo, hi, q0, q1);
    S.GDnom[k].Int2(lo, hi, r0, r1); S.GCnom[k].Int2(lo, hi, t0, t1);
    const double fr = p0 > 0 ? s0 / p0 : 0.;                // share of the class e3 in the F10 histogram
    i0 = fr * (0.3 * (o0 + q0) + 1.6 * (r0 + t0));
    i1 = fr * (0.3 * (o1 + q1) + 1.6 * (r1 + t1));
}

// ===========================================================================
// 3. CB from like-sign counts (WP0 rule)
// ===========================================================================
const int kNh16 = kNh15 + 1;   // + "other" (no slice / class)
struct LS { Vec pp, mm; };     // [fine i * 16 + h]
struct CBres {
    Vec cb;                    // fine bins [kNfine]
    double f15[kNBand][kNh15];
};
inline double FracFine(int i, double a, double b) {
    double fa = (a - kFineLo) / kFineW, fb = (b - kFineLo) / kFineW;
    if (std::fabs(fa - std::round(fa)) < 1e-6) fa = std::round(fa);
    if (std::fabs(fb - std::round(fb)) < 1e-6) fb = std::round(fb);
    const double o = std::min((double)(i + 1), fb) - std::max((double)i, fa);
    return o > 0 ? std::min(1., o) : 0.;
}
inline void BuildCB(const LS& L, CBres& out) {
    const int reb = kCBreb, nw = kNfine / reb;
    Vec wp(nw, 0.), wm(nw, 0.);
    for (int i = 0; i < kNfine; ++i)
        for (int h = 0; h < kNh16; ++h) { wp[i / reb] += L.pp[i * kNh16 + h]; wm[i / reb] += L.mm[i * kNh16 + h]; }
    Vec cbW(nw);
    for (int j = 0; j < nw; ++j) cbW[j] = 2. * std::sqrt(wp[j] * wm[j]);
    out.cb.assign(kNfine, 0.);
    for (int i = 0; i < kNfine; ++i) {
        const double t = (i + 0.5) / reb - 0.5;
        int k0, k1; double f;
        if (t <= 0) { k0 = k1 = 0; f = 0; }
        else if (t >= nw - 1) { k0 = k1 = nw - 1; f = 0; }
        else { k0 = (int)t; k1 = k0 + 1; f = t - k0; }
        out.cb[i] = ((1 - f) * cbW[k0] + f * cbW[k1]) / reb;
    }
    // fractions f(h | band) over the 15 valid cells; empty band -> pooled fraction
    double pool[kNh15] = {0}, ptot = 0;
    double cnt[kNBand][kNh15];
    for (int b = 0; b < kNBand; ++b) for (int h = 0; h < kNh15; ++h) cnt[b][h] = 0;
    for (int i = 0; i < kNfine; ++i) {
        const int b = BandOfFine(i);
        for (int h = 0; h < kNh15; ++h) { const double x = L.pp[i * kNh16 + h] + L.mm[i * kNh16 + h]; cnt[b][h] += x; pool[h] += x; ptot += x; }
    }
    for (int b = 0; b < kNBand; ++b) {
        double t = 0; for (int h = 0; h < kNh15; ++h) t += cnt[b][h];
        for (int h = 0; h < kNh15; ++h) out.f15[b][h] = t > 0 ? cnt[b][h] / t : (ptot > 0 ? pool[h] / ptot : 0.);
    }
}
// CB over [lo, hi) (fractional fine bins), all cells
inline double CBWin(const CBres& C, double lo, double hi) {
    double s = 0;
    for (int i = 0; i < kNfine; ++i) { const double fr = FracFine(i, lo, hi); if (fr > 0) s += fr * C.cb[i]; }
    return s;
}
// CB of a set of toy cells h15 over [lo, hi) (inside one band)
inline double CBCell(const CBres& C, double lo, double hi, int h15) {
    const int band = Band(0.5 * (lo + hi));
    return CBWin(C, lo, hi) * C.f15[band][h15];
}

// ===========================================================================
// 4. data inputs (like-sign counts, CB variances)
// ===========================================================================
struct DataIn {
    LS ls;
    CBres cb;
    Vec vCB10[10], vCB5[5];          // replica variance per F10 / F5g bin
    double covW[5][5];               // replica covariance of the window CBs: B, P, SBlo, SBhi, CRU
    Vec N10[10];                     // data N (for information only; never used in the toys)
    std::vector<double> otherPP, otherMM;
};
const int kWinIdx[5] = {kWB, kWP, kWSBlo, kWSBhi, kWCRU};
inline bool LoadData(TFile* f, const Grid& G, DataIn& D, FILE* log = nullptr) {
    D.ls.pp.assign(kNfine * kNh16, 0.); D.ls.mm.assign(kNfine * kNh16, 0.);
    TH3* g[2] = {(TH3*)f->Get("data/S0/all/Npp_G3"), (TH3*)f->Get("data/S0/all/Nmm_G3")};
    TH1* fine[2] = {(TH1*)f->Get("data/S0/all/Npp_fine_M4"), (TH1*)f->Get("data/S0/all/Nmm_fine_M4")};
    if (!g[0] || !g[1] || !fine[0] || !fine[1]) { printf("wp4: like-sign histograms missing\n"); return false; }
    const int i0 = (int)std::lround((kG3xlo - kFineLo) / kFineW);
    double maxdiff = 0, other = 0;
    for (int s = 0; s < 2; ++s) {
        Vec& out = s == 0 ? D.ls.pp : D.ls.mm;
        TH3* h = g[s];
        const int nx = h->GetNbinsX(), ny = h->GetNbinsY(), nz = h->GetNbinsZ();
        for (int x = 1; x <= nx; ++x) {
            const int i = i0 + x - 1;
            double sum = 0;
            for (int y = 0; y <= ny + 1; ++y)
                for (int z = 0; z <= nz + 1; ++z) {
                    const double c = h->GetBinContent(x, y, z);
                    if (c == 0) continue;
                    sum += c;
                    int cell = kNh15;
                    if (y >= 1 && y <= ny && z >= 1 && z <= nz) {
                        const int j = Slice(h->GetYaxis()->GetBinCenter(y));
                        const int e3 = E3class(h->GetZaxis()->GetBinCenter(z));
                        if (j >= 0 && e3 >= 0) cell = H15(j, e3);
                    }
                    if (cell == kNh15) other += c;
                    out[i * kNh16 + cell] += c;
                }
            maxdiff = std::max(maxdiff, std::fabs(sum - fine[s]->GetBinContent(i + 1)));
        }
        for (int i = 0; i < kNfine; ++i)
            if (i < i0 || i >= i0 + nx) out[i * kNh16 + kNh15] = fine[s]->GetBinContent(i + 1);
    }
    if (log) fprintf(log, "like-sign G3 vs fine_M4: max |difference| %.3g rows; like-sign rows without a toy cell in 0.40-0.85: %.0f\n", maxdiff, other);
    BuildCB(D.ls, D.cb);
    // replica variances
    for (int k = 0; k < 10; ++k) {
        TH2* r = (TH2*)f->Get(Form("data/S0/all/CBrep_%s", allio::F10(k / 2, k % 2).c_str()));
        D.vCB10[k].assign(G.nF + 2, 0.);
        if (!r) return false;
        const int nr = r->GetNbinsY();
        for (int x = 1; x <= G.nF; ++x) {
            double s = 0, q = 0;
            for (int y = 1; y <= nr; ++y) { const double c = r->GetBinContent(x, y); s += c; q += c * c; }
            D.vCB10[k][x] = std::max(0., q / nr - (s / nr) * (s / nr));
        }
        TH1* n = (TH1*)f->Get(Form("data/S0/all/N_%s", allio::F10(k / 2, k % 2).c_str()));
        D.N10[k].assign(G.nF + 2, 0.);
        for (int x = 1; x <= G.nF; ++x) D.N10[k][x] = n->GetBinContent(x);
    }
    for (int j = 0; j < kNj; ++j) {
        TH2* r = (TH2*)f->Get(Form("data/S0/all/CBrep_%s", allio::F5g(j).c_str()));
        D.vCB5[j].assign(G.nF + 2, 0.);
        if (!r) return false;
        const int nr = r->GetNbinsY();
        for (int x = 1; x <= G.nF; ++x) {
            double s = 0, q = 0;
            for (int y = 1; y <= nr; ++y) { const double c = r->GetBinContent(x, y); s += c; q += c * c; }
            D.vCB5[j][x] = std::max(0., q / nr - (s / nr) * (s / nr));
        }
    }
    {
        TH2* r = (TH2*)f->Get("data/S0/all/CBrep_win");
        if (!r) return false;
        const int nr = r->GetNbinsY();
        double m[5] = {0};
        for (int a = 0; a < 5; ++a) { for (int y = 1; y <= nr; ++y) m[a] += r->GetBinContent(kWinIdx[a] + 1, y); m[a] /= nr; }
        for (int a = 0; a < 5; ++a)
            for (int b = 0; b < 5; ++b) {
                double s = 0;
                for (int y = 1; y <= nr; ++y) s += (r->GetBinContent(kWinIdx[a] + 1, y) - m[a]) * (r->GetBinContent(kWinIdx[b] + 1, y) - m[b]);
                D.covW[a][b] = s / nr;
            }
    }
    return true;
}

// ===========================================================================
// 5. observed (toy) data on the toy cells
// ===========================================================================
struct Obs {
    Vec n;          // [h15 * (nR + 2) + r]
    LS ls;          // toy like-sign
    CBres cb;
    // derived
    Vec y10[10], cb10[10], y5[5], cb5[5];
    double NB = 0, NP = 0, NSB = 0, NU = 0, CBB = 0, CBP = 0, CBSB = 0, CBU = 0;
};
inline void Derive(const Grid& G, Obs& O) {
    BuildCB(O.ls, O.cb);
    const int st = G.nR + 2;
    for (int k = 0; k < 10; ++k) { O.y10[k].assign(G.nF + 2, 0.); O.cb10[k].assign(G.nF + 2, 0.); }
    for (int j = 0; j < kNj; ++j) { O.y5[j].assign(G.nF + 2, 0.); O.cb5[j].assign(G.nF + 2, 0.); }
    for (int h15 = 0; h15 < kNh15; ++h15) {
        const int j = h15 / kNe3, e3 = h15 % kNe3, k = 2 * j + E2of(e3);
        for (int r = 1; r <= G.nR; ++r) {
            const int x = G.refToF10[r];
            O.y10[k][x] += O.n[h15 * st + r];
            if (e3 == 2) O.y5[j][x] += O.n[h15 * st + r];
        }
    }
    // CB per F10 / F5g cell: CB(F10 bin) x f(cell | band)
    for (int x = 1; x <= G.nF; ++x) {
        const double lo = G.f10[x - 1], hi = G.f10[x];
        const double c1 = CBWin(O.cb, lo, hi);
        const int band = Band(0.5 * (lo + hi));
        for (int j = 0; j < kNj; ++j) {
            O.cb10[2 * j][x] = c1 * O.cb.f15[band][H15(j, 0)];
            O.cb10[2 * j + 1][x] = c1 * (O.cb.f15[band][H15(j, 1)] + O.cb.f15[band][H15(j, 2)]);
            O.cb5[j][x] = c1 * O.cb.f15[band][H15(j, 2)];
        }
    }
    auto nwin = [&](double lo, double hi) {
        double s = 0;
        for (int r : G.RefIn(lo, hi)) for (int h15 = 0; h15 < kNh15; ++h15) s += O.n[h15 * st + r];
        return s;
    };
    O.NB = nwin(kWin[kWB].lo, kWin[kWB].hi); O.NP = nwin(kWin[kWP].lo, kWin[kWP].hi);
    O.NSB = nwin(kWin[kWSBlo].lo, kWin[kWSBlo].hi) + nwin(kWin[kWSBhi].lo, kWin[kWSBhi].hi);
    O.NU = nwin(kWin[kWCRU].lo, kWin[kWCRU].hi);
    O.CBB = CBWin(O.cb, kWin[kWB].lo, kWin[kWB].hi); O.CBP = CBWin(O.cb, kWin[kWP].lo, kWin[kWP].hi);
    O.CBSB = CBWin(O.cb, kWin[kWSBlo].lo, kWin[kWSBlo].hi) + CBWin(O.cb, kWin[kWSBhi].lo, kWin[kWSBhi].hi);
    O.CBU = CBWin(O.cb, kWin[kWCRU].lo, kWin[kWCRU].hi);
}

// ===========================================================================
// 6. non-linear Poisson fit (Fisher scoring around bla::LinFit)
// ===========================================================================
enum SigT { kS2 = 0, kS3, kS4, kS5, kSeta, kNSig };
struct FitData {
    int n = 0;
    Vec y, cb, vcb, x, gd, gc, fix;
    std::vector<int> j, e;
    Vec sig[kNSig], vsig[kNSig];
    void Clear() {
        n = 0; y.clear(); cb.clear(); vcb.clear(); x.clear(); gd.clear(); gc.clear(); fix.clear(); j.clear(); e.clear();
        for (int t = 0; t < kNSig; ++t) { sig[t].clear(); vsig[t].clear(); }
    }
    void Push(double yy, double c, double vc, double xx, int jj, int ee, double d, double g, double fx, const double* s, const double* vs) {
        y.push_back(yy); cb.push_back(c); vcb.push_back(vc); x.push_back(xx); j.push_back(jj); e.push_back(ee);
        gd.push_back(d); gc.push_back(g); fix.push_back(fx);
        for (int t = 0; t < kNSig; ++t) { sig[t].push_back(s[t]); vsig[t].push_back(vs[t]); }
        ++n;
    }
};
struct Term { int t, p1, p2; };   // coefficient theta[p1] (x theta[p2] if p2 >= 0)
struct Spec {
    int np = 0;
    int iK[2] = {-1, -1};     // k_E1, k_E2 (-1: fixed to kFix)
    double kFix[2] = {1, 1};
    int iC[kNj] = {-1, -1, -1, -1, -1};
    int iB = -1, iKB = -1;
    double kBfix = 1;
    std::vector<Term> terms;
    std::vector<bla::Constraint> cons;
    Vec start;
    std::vector<std::string> names;
};
struct FitOut {
    bool ok = false;
    int nit = 0;
    Vec th;
    bla::Mat C, Cd, Cs;       // covariance with the template variance (BB-lite); Fisher with the data variance only;
                              // sandwich C M C, M = J^T V^-1 Vd V^-1 J (+ constraints): stat error of the BB-weighted estimator
    double chi2 = 0, chi2P = 0;
    int ndf = 0;
};
inline void Model(const FitData& D, const Spec& S, const Vec& th, Vec& mu, std::vector<Vec>* J, Vec* V, Vec* Vd) {
    const int n = D.n;
    mu.assign(n, 0.);
    if (J) J->assign(n, Vec(S.np, 0.));
    if (V) V->assign(n, 0.);
    if (Vd) Vd->assign(n, 0.);
    const double kB = S.iKB >= 0 ? th[S.iKB] : S.kBfix;
    const double b = S.iB >= 0 ? th[S.iB] : 0.;
    for (int i = 0; i < n; ++i) {
        const int e = D.e[i], j = D.j[i];
        const double k = S.iK[e] >= 0 ? th[S.iK[e]] : S.kFix[e];
        const double c = S.iC[j] >= 0 ? th[S.iC[j]] : 1.;
        const double bg = D.gd[i] + kB * D.gc[i];
        const double tl = 1. + b * D.x[i];
        double m = k * D.cb[i] + c * tl * bg + D.fix[i];
        double vs = 0;
        for (const auto& T : S.terms) {
            const double c2 = T.p2 >= 0 ? th[T.p2] : 1.;
            const double coef = th[T.p1] * c2;
            m += coef * D.sig[T.t][i];
            vs += coef * coef * D.vsig[T.t][i];
            if (J) {
                (*J)[i][T.p1] += c2 * D.sig[T.t][i];
                if (T.p2 >= 0) (*J)[i][T.p2] += th[T.p1] * D.sig[T.t][i];
            }
        }
        mu[i] = m;
        if (J) {
            if (S.iK[e] >= 0) (*J)[i][S.iK[e]] = D.cb[i];
            if (S.iC[j] >= 0) (*J)[i][S.iC[j]] = tl * bg;
            if (S.iB >= 0) (*J)[i][S.iB] = c * D.x[i] * bg;
            if (S.iKB >= 0) (*J)[i][S.iKB] = c * tl * D.gc[i];
        }
        const double vbase = std::max(m, 1.0) + k * k * D.vcb[i];
        if (V) (*V)[i] = vbase + vs;
        if (Vd) (*Vd)[i] = vbase;
    }
}
inline bool CovAt(const FitData& D, const Spec& S, const std::vector<Vec>& J, const Vec& V, bla::Mat& C) {
    const int np = S.np;
    bla::Mat A = bla::Zero(np, np);
    for (int i = 0; i < D.n; ++i) {
        const double w = 1.0 / V[i];
        for (int k = 0; k < np; ++k) {
            const double fk = J[i][k] * w;
            if (fk == 0) continue;
            for (int l = 0; l <= k; ++l) A[k][l] += fk * J[i][l];
        }
    }
    for (const auto& c : S.cons) A[c.k][c.k] += 1.0 / (c.sig * c.sig);
    for (int k = 0; k < np; ++k) for (int l = k + 1; l < np; ++l) A[k][l] = A[l][k];
    Vec sc(np, 1.0);
    for (int k = 0; k < np; ++k) sc[k] = A[k][k] > 0 ? 1.0 / std::sqrt(A[k][k]) : 1.0;
    for (int k = 0; k < np; ++k) for (int l = 0; l < np; ++l) A[k][l] *= sc[k] * sc[l];
    if (!bla::CholInv(A)) return false;
    for (int k = 0; k < np; ++k) for (int l = 0; l < np; ++l) A[k][l] *= sc[k] * sc[l];
    C = A;
    return true;
}
inline FitOut NLFit(const FitData& D, const Spec& S, const Vec* start = nullptr, int maxIt = 80) {
    FitOut out;
    const int n = D.n, np = S.np;
    Vec th = start ? *start : S.start;
    Vec mu, V;
    std::vector<Vec> J;
    bool conv = false;
    for (int it = 0; it < maxIt; ++it) {
        Model(D, S, th, mu, &J, &V, nullptr);
        bla::LinProblem P;
        P.n = n; P.np = np; P.poissonVar = false; P.maxIt = 1; P.cons = S.cons;
        P.y.resize(n); P.F = J; P.var0 = V;
        for (int i = 0; i < n; ++i) { double s = 0; for (int k = 0; k < np; ++k) s += J[i][k] * th[k]; P.y[i] = D.y[i] - mu[i] + s; }
        const bla::LinResult r = bla::LinFit(P, &th);
        if (!r.ok) { out.nit = it; return out; }
        double dmax = 0, lim = 1.0;
        for (int k = 0; k < np; ++k) {
            const double sk = std::sqrt(std::max(r.C[k][k], 1e-300));
            const double d = std::fabs(r.a[k] - th[k]) / sk;
            dmax = std::max(dmax, d);
            if (d > 3.0) lim = std::min(lim, 3.0 / d);
        }
        for (int k = 0; k < np; ++k) th[k] += lim * (r.a[k] - th[k]);
        out.nit = it + 1;
        if (dmax < 1e-5) { conv = true; break; }
    }
    if (!conv) return out;
    Vec Vd;
    Model(D, S, th, mu, &J, &V, &Vd);
    if (!CovAt(D, S, J, V, out.C) || !CovAt(D, S, J, Vd, out.Cd)) return out;
    {   // sandwich
        bla::Mat M = bla::Zero(np, np);
        for (int i = 0; i < n; ++i) {
            const double w = Vd[i] / (V[i] * V[i]);
            for (int k = 0; k < np; ++k) { const double fk = J[i][k] * w; if (fk == 0) continue; for (int l = 0; l < np; ++l) M[k][l] += fk * J[i][l]; }
        }
        for (const auto& c : S.cons) M[c.k][c.k] += 1.0 / (c.sig * c.sig);
        out.Cs = bla::Zero(np, np);
        bla::Mat T = bla::Zero(np, np);
        for (int k = 0; k < np; ++k) for (int l = 0; l < np; ++l) { double x = 0; for (int m = 0; m < np; ++m) x += out.C[k][m] * M[m][l]; T[k][l] = x; }
        for (int k = 0; k < np; ++k) for (int l = 0; l < np; ++l) { double x = 0; for (int m = 0; m < np; ++m) x += T[k][m] * out.C[m][l]; out.Cs[k][l] = x; }
    }
    out.th = th;
    out.chi2 = 0; out.chi2P = 0;
    for (int i = 0; i < n; ++i) {
        const double m = std::max(mu[i], 1e-6);
        out.chi2 += 2.0 * (m - D.y[i] + (D.y[i] > 0 ? D.y[i] * std::log(D.y[i] / m) : 0.0));
        out.chi2P += (D.y[i] - mu[i]) * (D.y[i] - mu[i]) / V[i];
    }
    out.ndf = n - np + (int)S.cons.size();
    out.ok = true;
    return out;
}
inline double VarOf(const bla::Mat& C, const Vec& g) {
    double s = 0;
    for (size_t k = 0; k < g.size(); ++k) for (size_t l = 0; l < g.size(); ++l) s += g[k] * C[k][l] * g[l];
    return s;
}

// ===========================================================================
// 7. the fit inputs of one toy
// ===========================================================================
struct FitTmpl {        // fit templates (half A or B, nominal mass) and derived integrals
    const Tmpl* T = nullptr;
    BkgSmooth B;
    // F10 bin integrals of the smoothed background (per F10 histogram) and F5g
    Vec gd10[10], gc10[10], gd5[5], gc5[5];
    // refined-cell integrals per F10 histogram (for window sums of stage A)
    Vec gdR[10], gcR[10];
    // window sums of the signal templates (all cells): [sig][win] win 0 B, 1 P, 2 SB, 3 U
    double W[kNSig][4];
};
inline void MakeFitTmpl(const Tmpl& T, const Grid& G, FitTmpl& F) {
    F.T = &T;
    MakeFitSmooth(T, G, F.B);
    for (int k = 0; k < 10; ++k) {
        F.gd10[k].assign(G.nF + 2, 0.); F.gc10[k].assign(G.nF + 2, 0.);
        for (int x = 1; x <= G.nF; ++x) { F.gd10[k][x] = F.B.GD10[k].Int(G.f10[x - 1], G.f10[x]); F.gc10[k][x] = F.B.GC10[k].Int(G.f10[x - 1], G.f10[x]); }
        F.gdR[k].assign(G.nR + 2, 0.); F.gcR[k].assign(G.nR + 2, 0.);
        for (int r = 1; r <= G.nR; ++r) { F.gdR[k][r] = F.B.GD10[k].Int(G.ref[r - 1], G.ref[r]); F.gcR[k][r] = F.B.GC10[k].Int(G.ref[r - 1], G.ref[r]); }
    }
    for (int j = 0; j < kNj; ++j) {
        F.gd5[j].assign(G.nF + 2, 0.); F.gc5[j].assign(G.nF + 2, 0.);
        for (int x = 1; x <= G.nF; ++x) { F.gd5[j][x] = F.B.GD5[j].Int(G.f10[x - 1], G.f10[x]); F.gc5[j][x] = F.B.GC5[j].Int(G.f10[x - 1], G.f10[x]); }
    }
    const int sc[kNSig] = {kTC2, kTC3, kTC4, kTC5, kTGeta};
    for (int t = 0; t < kNSig; ++t) {
        F.W[t][0] = WinSum(T, G, sc[t], kWin[kWB].lo, kWin[kWB].hi);
        F.W[t][1] = WinSum(T, G, sc[t], kWin[kWP].lo, kWin[kWP].hi);
        F.W[t][2] = WinSum(T, G, sc[t], kWin[kWSBlo].lo, kWin[kWSBlo].hi) + WinSum(T, G, sc[t], kWin[kWSBhi].lo, kWin[kWSBhi].hi);
        F.W[t][3] = WinSum(T, G, sc[t], kWin[kWCRU].lo, kWin[kWCRU].hi);
    }
}
const int kSigComp[kNSig] = {kTC2, kTC3, kTC4, kTC5, kTGeta};
// bootstrap replicas of the full sample (hists_v1.root sim/S0__nom/<comp>_<layout>_rep, 100 Poisson event-group
// replicas made by WP0): fit templates only (F10, F5g = class e3 2 of F15, window sums); the variances are those of
// the full sample (Tfull); G_D and G_C enter as totals (stored in the "_nom" slots, "_om" slots empty)
inline bool LoadReplicas(TFile* fh, const Grid& G, const Tmpl& Tfull, int nrep, std::vector<FitTmpl*>& out, std::vector<Tmpl*>& keep) {
    const int nc = 7;
    const int tc[nc] = {kTC2, kTC3, kTC4, kTC5, kTGeta, kTGDnom, kTGCnom};
    const int hc[nc] = {kC2, kC3, kC4, kC5, kGeta, kGD, kGconv};
    keep.assign(nrep, nullptr);
    for (int r = 0; r < nrep; ++r) {
        keep[r] = new Tmpl();
        Tmpl& T = *keep[r];
        for (int c = 0; c < kNTComp; ++c) {
            for (int h = 0; h < kNh15; ++h) { T.R[c][h].assign(G.nR + 2, 0.); T.F15[c][h].assign(G.nF + 2, 0.); T.F15v[c][h].assign(G.nF + 2, 0.); }
            for (int k = 0; k < 10; ++k) { T.F10[c][k].assign(G.nF + 2, 0.); T.F10v[c][k].assign(G.nF + 2, 0.); }
            T.tot[c] = 0;
        }
    }
    std::vector<std::vector<std::array<double, 4>>> W(nrep, std::vector<std::array<double, 4>>(kNSig));
    for (int q = 0; q < nc; ++q) {
        const int c = tc[q];
        const bool gd = c == kTGDnom, gc = c == kTGCnom;
        for (int k = 0; k < 10; ++k) {
            TH2* h = (TH2*)fh->Get(Form("sim/S0__nom/%s_%s_rep", CompTag(hc[q]), F10Name(k / 2, k % 2).c_str()));
            if (!h || h->GetNbinsY() < nrep) { printf("wp4: replica F10 %s missing\n", CompTag(hc[q])); return false; }
            for (int r = 0; r < nrep; ++r)
                for (int x = 1; x <= G.nF; ++x) {
                    keep[r]->F10[c][k][x] = h->GetBinContent(x, r + 1) * kNbToRaw;
                    keep[r]->F10v[c][k][x] = gd ? Tfull.F10v[kTGDom][k][x] + Tfull.F10v[kTGDnom][k][x]
                                           : gc ? Tfull.F10v[kTGCom][k][x] + Tfull.F10v[kTGCnom][k][x] : Tfull.F10v[c][k][x];
                }
        }
        for (int j = 0; j < kNj; ++j) {
            TH2* h = (TH2*)fh->Get(Form("sim/S0__nom/%s_%s_rep", CompTag(hc[q]), F5gName(j).c_str()));
            if (!h || h->GetNbinsY() < nrep) { printf("wp4: replica F5g %s missing\n", CompTag(hc[q])); return false; }
            const int h15 = H15(j, 2);
            for (int r = 0; r < nrep; ++r)
                for (int x = 1; x <= G.nF; ++x) {
                    keep[r]->F15[c][h15][x] = h->GetBinContent(x, r + 1) * kNbToRaw;
                    keep[r]->F15v[c][h15][x] = gd ? Tfull.F15v[kTGDom][h15][x] + Tfull.F15v[kTGDnom][h15][x]
                                             : gc ? Tfull.F15v[kTGCom][h15][x] + Tfull.F15v[kTGCnom][h15][x] : Tfull.F15v[c][h15][x];
                }
        }
        if (q < kNSig) {
            TH2* h = (TH2*)fh->Get(Form("sim/S0__nom/%s_win_rep", CompTag(hc[q])));
            if (!h) return false;
            for (int r = 0; r < nrep; ++r) {
                W[r][q][0] = h->GetBinContent(kWB + 1, r + 1) * kNbToRaw;
                W[r][q][1] = h->GetBinContent(kWP + 1, r + 1) * kNbToRaw;
                W[r][q][2] = (h->GetBinContent(kWSBlo + 1, r + 1) + h->GetBinContent(kWSBhi + 1, r + 1)) * kNbToRaw;
                W[r][q][3] = h->GetBinContent(kWCRU + 1, r + 1) * kNbToRaw;
            }
        }
    }
    out.assign(nrep, nullptr);
    for (int r = 0; r < nrep; ++r) {
        out[r] = new FitTmpl();
        MakeFitTmpl(*keep[r], G, *out[r]);
        for (int t = 0; t < kNSig; ++t) for (int w = 0; w < 4; ++w) out[r]->W[t][w] = W[r][t][w];
    }
    return true;
}
// F10 fit data; sel(k, x) chooses the bins (all: F-main / F-free; CR: stage A); fixEta: per-bin fixed eta term
inline void MakeF10Data(const Obs& O, const DataIn& Din, const FitTmpl& F, const Grid& G, FitData& D,
                        const std::function<bool(int, int)>& sel, const std::function<double(int, int)>& fixEta) {
    D.Clear();
    const Tmpl& T = *F.T;
    for (int k = 0; k < 10; ++k)
        for (int x = 1; x <= G.nF; ++x) {
            if (!sel(k, x)) continue;
            double s[kNSig], vs[kNSig];
            for (int t = 0; t < kNSig; ++t) { s[t] = T.F10[kSigComp[t]][k][x]; vs[t] = T.F10v[kSigComp[t]][k][x]; }
            D.Push(O.y10[k][x], O.cb10[k][x], Din.vCB10[k][x], G.xF[x], k / 2, k % 2, F.gd10[k][x], F.gc10[k][x],
                   fixEta ? fixEta(k, x) : 0., s, vs);
        }
}
inline void MakeF5gData(const Obs& O, const DataIn& Din, const FitTmpl& F, const Grid& G, FitData& D,
                        const std::function<double(int, int)>& fixEta) {
    D.Clear();
    const Tmpl& T = *F.T;
    for (int j = 0; j < kNj; ++j)
        for (int x = 1; x <= G.nF; ++x) {
            const int h = H15(j, 2);
            double s[kNSig], vs[kNSig];
            for (int t = 0; t < kNSig; ++t) { s[t] = T.F15[kSigComp[t]][h][x]; vs[t] = T.F15v[kSigComp[t]][h][x]; }
            D.Push(O.y5[j][x], O.cb5[j][x], Din.vCB5[j][x], G.xF[x], j, 1, F.gd5[j][x], F.gc5[j][x], fixEta ? fixEta(j, x) : 0., s, vs);
        }
}

// parameter layouts
// F-main: 0 beta2, 1 beta34, 2 beta5, 3 g, 4-8 c_j, 9 b, 10 kB, 11 kE1, 12 kE2
inline Spec SpecFMain() {
    Spec S; S.np = 13;
    S.iC[0] = 4; S.iC[1] = 5; S.iC[2] = 6; S.iC[3] = 7; S.iC[4] = 8; S.iB = 9; S.iKB = 10; S.iK[0] = 11; S.iK[1] = 12;
    S.terms = {{kS2, 0, -1}, {kS3, 1, -1}, {kS4, 1, -1}, {kS5, 2, -1}, {kSeta, 0, 3}};
    S.cons = {{3, 1.0, 0.30}, {12, 1.0, 0.05}};
    S.start = {0.5, 0.5, 0.5, 1.0, 1, 1, 1, 1, 1, 0.0, 1.0, 1.0, 1.0};
    S.names = {"beta2", "beta34", "beta5", "g", "c1", "c2", "c3", "c4", "c5", "b", "kB", "kE1", "kE2"};
    return S;
}
// F-free: 0 a2, 1 a3, 2 a4, 3 a5, 4 g, 5-9 c_j, 10 b, 11 kB, 12 kE1, 13 kE2
inline Spec SpecFFree() {
    Spec S; S.np = 14;
    for (int j = 0; j < kNj; ++j) S.iC[j] = 5 + j;
    S.iB = 10; S.iKB = 11; S.iK[0] = 12; S.iK[1] = 13;
    S.terms = {{kS2, 0, -1}, {kS3, 1, -1}, {kS4, 2, -1}, {kS5, 3, -1}, {kSeta, 0, 4}};
    S.cons = {{4, 1.0, 0.30}, {13, 1.0, 0.05}};
    S.start = {0.5, 0.5, 0.5, 0.5, 1.0, 1, 1, 1, 1, 1, 0.0, 1.0, 1.0, 1.0};
    S.names = {"a2", "a3", "a4", "a5", "g", "c1", "c2", "c3", "c4", "c5", "b", "kB", "kE1", "kE2"};
    return S;
}
// stage A: 0-4 c_j, 5 b, 6 kB, 7 kE1, 8 kE2
inline Spec SpecStageA() {
    Spec S; S.np = 9;
    for (int j = 0; j < kNj; ++j) S.iC[j] = j;
    S.iB = 5; S.iKB = 6; S.iK[0] = 7; S.iK[1] = 8;
    S.cons = {{8, 1.0, 0.05}};
    S.start = {1, 1, 1, 1, 1, 0.0, 1.0, 1.0, 1.0};
    S.names = {"c1", "c2", "c3", "c4", "c5", "b", "kB", "kE1", "kE2"};
    return S;
}
// F-gold (F5g, class E2 only): 0 a2, 1 a5, 2 g, 3-7 c_j, 8 b, 9 kE2; kappa_B fixed
inline Spec SpecFGold(double kB) {
    Spec S; S.np = 10;
    for (int j = 0; j < kNj; ++j) S.iC[j] = 3 + j;
    S.iB = 8; S.iKB = -1; S.kBfix = kB; S.iK[0] = -1; S.iK[1] = 9;
    S.terms = {{kS2, 0, -1}, {kS5, 1, -1}, {kSeta, 0, 2}};
    S.cons = {{2, 1.0, 0.30}, {9, 1.0, 0.05}};
    S.start = {0.5, 0.5, 1.0, 1, 1, 1, 1, 1, 0.0, 1.0};
    S.names = {"a2", "a5", "g", "c1", "c2", "c3", "c4", "c5", "b", "kE2"};
    return S;
}
inline bool InCR(const Grid& G, int k, int x) {
    const double lo = G.f10[x - 1], hi = G.f10[x];
    return InCRU(lo, hi, k / 2) || InCRK(lo, hi, k / 2);
}

// stage-A background in a set of refined cells [lo, hi): sum over F10 histograms of c_j (1 + b x_r)(GD + kB GC)
inline double BkA(const FitTmpl& F, const Grid& G, const Vec& th, double lo, double hi) {
    double s = 0;
    for (int r : G.RefIn(lo, hi))
        for (int k = 0; k < 10; ++k) s += th[k / 2] * (1. + th[5] * G.xc[r]) * (F.gdR[k][r] + th[6] * F.gcR[k][r]);
    return s;
}

// ===========================================================================
// 8. estimators of one toy
// ===========================================================================
// quantity ids
enum Qid {
    qFMa = 0, qFMkap, qFMrho, qFMg, qFMYb, qFMYp, qFMY5, qFMY4,
    qFFa2, qFFa3, qFFa4, qFFa5, qFFYb, qFFYp, qFFr5, qFFR34, qFFBRC,
    qFGa2, qFGa5, qFGBRG,
    qB4Yb, qPSBYp, kNQ };
inline const char* QName(int q) {
    static const char* n[kNQ] = {"FM_a", "FM_kappa", "FM_rho", "FM_g", "FM_Ybump", "FM_Ypeak", "FM_Y5", "FM_Y4",
                                 "FF_a2", "FF_a3", "FF_a4", "FF_a5", "FF_Ybump", "FF_Ypeak", "FF_r5", "FF_R34", "FF_BRC",
                                 "FG_a2", "FG_a5", "FG_BRG", "B4_Ybump", "PSB_Ypeak"};
    return n[q];
}
struct Result {
    double est[kNQ], err[kNQ], errd[kNQ], errs[kNQ];   // fit (BB-lite), Fisher data-only, sandwich (stat of the BB-weighted fit)
    bool ok[kNQ];
    double chi2FM = 0, chi2FF = 0, chi2FG = 0, chi2SA = 0;
    int ndfFM = 0, ndfFF = 0, ndfFG = 0, ndfSA = 0;
    int nitFM = 0;
    double fmTh[13];
    double saTh[9];
    double rhoAbE = 0, rhoAbB = 0;   // corr(a, b) of F-main and corr(rho_d, kappa)
};

// ratio helpers (delta method) of linear functionals u^T th / v^T th
inline void Ratio(const Vec& th, const bla::Mat& C, const bla::Mat& Cd, const bla::Mat& Cs, const Vec& u, const Vec& v, double& r, double& e,
                  double& ed, double& es) {
    double U = 0, Vv = 0;
    for (size_t k = 0; k < th.size(); ++k) { U += u[k] * th[k]; Vv += v[k] * th[k]; }
    r = U / Vv;
    Vec g(th.size(), 0.);
    for (size_t k = 0; k < th.size(); ++k) g[k] = (u[k] - r * v[k]) / Vv;
    e = std::sqrt(std::max(0., VarOf(C, g))); ed = std::sqrt(std::max(0., VarOf(Cd, g))); es = std::sqrt(std::max(0., VarOf(Cs, g)));
}

inline void RunEstimators(const Obs& O, const DataIn& Din, const FitTmpl& F, const Grid& G, Result& R) {
    for (int q = 0; q < kNQ; ++q) { R.est[q] = R.err[q] = R.errd[q] = R.errs[q] = 0; R.ok[q] = false; }
    FitData D;
    // ---------------- F-main
    MakeF10Data(O, Din, F, G, D, [](int, int) { return true; }, nullptr);
    const Spec SM = SpecFMain();
    const FitOut fm = NLFit(D, SM);
    R.nitFM = fm.nit;
    if (!fm.ok) return;
    R.chi2FM = fm.chi2; R.ndfFM = fm.ndf;
    for (int k = 0; k < 13; ++k) R.fmTh[k] = fm.th[k];
    const Vec& t = fm.th;
    auto lin = [&](int q, const Vec& g, double scale) {
        double s = 0; for (int k = 0; k < (int)g.size(); ++k) s += g[k] * t[k];
        R.est[q] = s * scale; R.err[q] = std::sqrt(std::max(0., VarOf(fm.C, g))) * scale;
        R.errd[q] = std::sqrt(std::max(0., VarOf(fm.Cd, g))) * scale; R.errs[q] = std::sqrt(std::max(0., VarOf(fm.Cs, g))) * scale; R.ok[q] = true;
    };
    {
        Vec g(13, 0.);
        g[0] = 1; lin(qFMa, g, 1.); g[0] = 0;
        g[3] = 1; lin(qFMg, g, 1.); g[3] = 0;
        g[0] = F.W[kS2][0]; g[1] = F.W[kS3][0]; lin(qFMYb, g, kYf); g[0] = g[1] = 0;
        g[1] = F.W[kS4][1]; g[2] = F.W[kS5][1]; lin(qFMYp, g, kYf); g[1] = g[2] = 0;
        g[2] = F.W[kS5][1]; lin(qFMY5, g, kYf); g[2] = 0;
        g[1] = F.W[kS4][1]; lin(qFMY4, g, kYf); g[1] = 0;
        Vec u(13, 0.), v(13, 0.);
        v[0] = 1;
        u[1] = 1; Ratio(t, fm.C, fm.Cd, fm.Cs, u, v, R.est[qFMkap], R.err[qFMkap], R.errd[qFMkap], R.errs[qFMkap]); R.ok[qFMkap] = true; u[1] = 0;
        u[2] = 1; Ratio(t, fm.C, fm.Cd, fm.Cs, u, v, R.est[qFMrho], R.err[qFMrho], R.errd[qFMrho], R.errs[qFMrho]); R.ok[qFMrho] = true;
        // correlations: corr(a, b), corr(rho, kappa)
        R.rhoAbE = fm.C[0][9] / std::sqrt(fm.C[0][0] * fm.C[9][9]);
        // rho = beta5/beta2, kappa = beta34/beta2: gradients
        Vec gr(13, 0.), gk(13, 0.);
        gr[2] = 1 / t[0]; gr[0] = -t[2] / (t[0] * t[0]);
        gk[1] = 1 / t[0]; gk[0] = -t[1] / (t[0] * t[0]);
        double crk = 0;
        for (int a = 0; a < 13; ++a) for (int b = 0; b < 13; ++b) crk += gr[a] * fm.C[a][b] * gk[b];
        R.rhoAbB = crk / std::sqrt(std::max(1e-300, VarOf(fm.C, gr) * VarOf(fm.C, gk)));
    }
    // eta term of F-main in an F10 cell / F5g cell (fixed terms of stage A and F-gold)
    const Tmpl& T = *F.T;
    auto etaF10 = [&](int k, int x) {
        return t[0] * T.F10[kTC2][k][x] + t[1] * (T.F10[kTC3][k][x] + T.F10[kTC4][k][x]) + t[2] * T.F10[kTC5][k][x] +
               t[0] * t[3] * T.F10[kTGeta][k][x];
    };
    // ---------------- F-free
    {
        const Spec SF = SpecFFree();
        Vec st = {t[0], t[1], t[1], t[2], t[3], t[4], t[5], t[6], t[7], t[8], t[9], t[10], t[11], t[12]};
        const FitOut ff = NLFit(D, SF, &st);
        if (ff.ok) {
            R.chi2FF = ff.chi2; R.ndfFF = ff.ndf;
            const Vec& a = ff.th;
            auto linF = [&](int q, const Vec& g, double scale) {
                double s = 0; for (int k = 0; k < (int)g.size(); ++k) s += g[k] * a[k];
                R.est[q] = s * scale; R.err[q] = std::sqrt(std::max(0., VarOf(ff.C, g))) * scale;
                R.errd[q] = std::sqrt(std::max(0., VarOf(ff.Cd, g))) * scale; R.errs[q] = std::sqrt(std::max(0., VarOf(ff.Cs, g))) * scale;
                R.ok[q] = true;
            };
            Vec g(14, 0.);
            for (int k = 0; k < 4; ++k) { g[k] = 1; linF(qFFa2 + k, g, 1.); g[k] = 0; }
            g[0] = F.W[kS2][0]; g[1] = F.W[kS3][0]; linF(qFFYb, g, kYf); g[0] = g[1] = 0;
            g[2] = F.W[kS4][1]; g[3] = F.W[kS5][1]; linF(qFFYp, g, kYf); g[2] = g[3] = 0;
            Vec u(14, 0.), v(14, 0.);
            u[3] = 1; v[0] = 1; Ratio(a, ff.C, ff.Cd, ff.Cs, u, v, R.est[qFFr5], R.err[qFFr5], R.errd[qFFr5], R.errs[qFFr5]); R.ok[qFFr5] = true;
            u.assign(14, 0.); v.assign(14, 0.); u[2] = 1; v[1] = 1;
            Ratio(a, ff.C, ff.Cd, ff.Cs, u, v, R.est[qFFR34], R.err[qFFR34], R.errd[qFFR34], R.errs[qFFR34]); R.ok[qFFR34] = true;
            // BR_C = 2.68e-4 a5 (a3/a2)/a4: log-derivative
            const double brc = kBRdir * a[3] * a[1] / (a[0] * a[2]);
            Vec gl(14, 0.); gl[3] = brc / a[3]; gl[1] = brc / a[1]; gl[0] = -brc / a[0]; gl[2] = -brc / a[2];
            R.est[qFFBRC] = brc; R.err[qFFBRC] = std::sqrt(std::max(0., VarOf(ff.C, gl)));
            R.errd[qFFBRC] = std::sqrt(std::max(0., VarOf(ff.Cd, gl))); R.errs[qFFBRC] = std::sqrt(std::max(0., VarOf(ff.Cs, gl)));
            R.ok[qFFBRC] = true;
        }
    }
    // ---------------- stage A (CR-U + CR-K), eta fixed to F-main
    FitData DA;
    MakeF10Data(O, Din, F, G, DA, [&](int k, int x) { return InCR(G, k, x); }, etaF10);
    const Spec SA = SpecStageA();
    Vec sta = {t[4], t[5], t[6], t[7], t[8], t[9], t[10], t[11], t[12]};
    const FitOut sa = NLFit(DA, SA, &sta);
    if (sa.ok) {
        R.chi2SA = sa.chi2; R.ndfSA = sa.ndf;
        for (int k = 0; k < 9; ++k) R.saTh[k] = sa.th[k];
        // eta template sums (F-main parameters) in windows: idx 0 B, 1 P, 2 SB, 3 U
        auto etaW = [&](const Vec& p, int w, bool c2, bool c3, bool c4, bool c5, bool ge) {
            return (c2 ? p[0] * F.W[kS2][w] : 0.) + (c3 ? p[1] * F.W[kS3][w] : 0.) + (c4 ? p[1] * F.W[kS4][w] : 0.) +
                   (c5 ? p[2] * F.W[kS5][w] : 0.) + (ge ? p[0] * p[3] * F.W[kSeta][w] : 0.);
        };
        const double lo[4] = {kWin[kWB].lo, kWin[kWP].lo, 0, kWin[kWCRU].lo}, hi[4] = {kWin[kWB].hi, kWin[kWP].hi, 0, kWin[kWCRU].hi};
        auto bk = [&](const Vec& s, int w) {
            if (w == 2) return BkA(F, G, s, kWin[kWSBlo].lo, kWin[kWSBlo].hi) + BkA(F, G, s, kWin[kWSBhi].lo, kWin[kWSBhi].hi);
            return BkA(F, G, s, lo[w], hi[w]);
        };
        const double sB = O.NB - O.CBB, sU = O.NU - O.CBU, sP = O.NP - O.CBP, sSB = O.NSB - O.CBSB;
        // B4
        auto b4 = [&](const Vec& s, const Vec& p) {
            const double Rr = bk(s, 0) / bk(s, 3);
            return kYf * (sB - (sU - etaW(p, 3, true, true, true, true, true)) * Rr - etaW(p, 0, false, false, true, true, true));
        };
        // peak side bands; leakage of the peak into the side bands with the unweighted C4 + C5 shape (SPEC WP6: "iterated
        // out with C4 + C5"; with the fitted kappa, rho the ratio is undefined when the fitted peak amplitude is near 0)
        const double fleak = (F.W[kS4][2] + F.W[kS5][2]) / (F.W[kS4][1] + F.W[kS5][1]);
        auto psb = [&](const Vec& s, const Vec& p) {
            const double bP = bk(s, 1) + etaW(p, 1, true, true, false, false, true);
            const double bS = bk(s, 2) + etaW(p, 2, true, true, false, false, true);
            const double Rp = bP / bS;
            return kYf * (sP - Rp * sSB) / (1. - Rp * fleak);
        };
        Vec pF(t.begin(), t.begin() + 4);
        const Vec& s0 = sa.th;
        R.est[qB4Yb] = b4(s0, pF); R.est[qPSBYp] = psb(s0, pF);
        // errors: data counts + window CB (replica covariance) + stage-A covariance + F-main covariance (numerical gradients)
        auto numgrad = [&](const std::function<double(const Vec&)>& fn, const Vec& x0, const bla::Mat& C) {
            Vec g(x0.size(), 0.);
            for (size_t k = 0; k < x0.size(); ++k) {
                const double h = 1e-4 * std::sqrt(std::max(C[k][k], 1e-20));
                Vec xp = x0, xm = x0; xp[k] += h; xm[k] -= h;
                g[k] = (fn(xp) - fn(xm)) / (2 * h);
            }
            return g;
        };
        // F-main parameter block (beta2, beta34, beta5, g) covariance
        bla::Mat CF = bla::Zero(4, 4), CFd = bla::Zero(4, 4), CFs = bla::Zero(4, 4);
        for (int a = 0; a < 4; ++a) for (int b = 0; b < 4; ++b) { CF[a][b] = fm.C[a][b]; CFd[a][b] = fm.Cd[a][b]; CFs[a][b] = fm.Cs[a][b]; }
        {
            const double Rr = bk(s0, 0) / bk(s0, 3);
            // d Y / d(N_B, CB_B, N_U, CB_U): kYf (1, -1, -R, +R)
            const double vN = O.NB + Rr * Rr * O.NU;
            const double vCB = Din.covW[0][0] + Rr * Rr * Din.covW[4][4] - 2 * Rr * Din.covW[0][4];
            const Vec gs = numgrad([&](const Vec& s) { return b4(s, pF); }, s0, sa.C);
            const Vec gp = numgrad([&](const Vec& p) { return b4(s0, p); }, pF, CF);
            const double v = kYf * kYf * (vN + vCB) + VarOf(sa.C, gs) + VarOf(CF, gp);
            const double vd = kYf * kYf * (vN + vCB) + VarOf(sa.Cd, gs) + VarOf(CFd, gp);
            const double vs = kYf * kYf * (vN + vCB) + VarOf(sa.Cs, gs) + VarOf(CFs, gp);
            R.err[qB4Yb] = std::sqrt(std::max(0., v)); R.errd[qB4Yb] = std::sqrt(std::max(0., vd)); R.errs[qB4Yb] = std::sqrt(std::max(0., vs));
            R.ok[qB4Yb] = true;
        }
        {
            const double bP = bk(s0, 1) + etaW(pF, 1, true, true, false, false, true);
            const double bS = bk(s0, 2) + etaW(pF, 2, true, true, false, false, true);
            const double Rp = bP / bS;
            const double Dn = 1. - Rp * fleak;
            const double vN = (O.NP + Rp * Rp * O.NSB) / (Dn * Dn);
            const double cSB = Din.covW[2][2] + Din.covW[3][3] + 2 * Din.covW[2][3];
            const double cPSB = Din.covW[1][2] + Din.covW[1][3];
            const double vCB = (Din.covW[1][1] + Rp * Rp * cSB - 2 * Rp * cPSB) / (Dn * Dn);
            const Vec gs = numgrad([&](const Vec& s) { return psb(s, pF); }, s0, sa.C);
            const Vec gp = numgrad([&](const Vec& p) { return psb(s0, p); }, pF, CF);
            const double v = kYf * kYf * (vN + vCB) + VarOf(sa.C, gs) + VarOf(CF, gp);
            const double vd = kYf * kYf * (vN + vCB) + VarOf(sa.Cd, gs) + VarOf(CFd, gp);
            const double vs = kYf * kYf * (vN + vCB) + VarOf(sa.Cs, gs) + VarOf(CFs, gp);
            R.err[qPSBYp] = std::sqrt(std::max(0., v)); R.errd[qPSBYp] = std::sqrt(std::max(0., vd)); R.errs[qPSBYp] = std::sqrt(std::max(0., vs));
            R.ok[qPSBYp] = true;
        }
    }
    // ---------------- F-gold (F5g), C3 + C4 fixed from F-main, kappa_B fixed
    {
        FitData DG;
        MakeF5gData(O, Din, F, G, DG, [&](int j, int x) {
            const int h = H15(j, 2);
            return t[1] * (T.F15[kTC3][h][x] + T.F15[kTC4][h][x]);
        });
        const Spec SG = SpecFGold(t[10]);
        Vec stg = {t[0], t[2], t[3], t[4], t[5], t[6], t[7], t[8], t[9], t[12]};
        const FitOut fg = NLFit(DG, SG, &stg);
        if (fg.ok) {
            R.chi2FG = fg.chi2; R.ndfFG = fg.ndf;
            const Vec& a = fg.th;
            R.est[qFGa2] = a[0]; R.err[qFGa2] = std::sqrt(fg.C[0][0]); R.errd[qFGa2] = std::sqrt(fg.Cd[0][0]); R.errs[qFGa2] = std::sqrt(fg.Cs[0][0]);
            R.ok[qFGa2] = true;
            R.est[qFGa5] = a[1]; R.err[qFGa5] = std::sqrt(fg.C[1][1]); R.errd[qFGa5] = std::sqrt(fg.Cd[1][1]); R.errs[qFGa5] = std::sqrt(fg.Cs[1][1]);
            R.ok[qFGa5] = true;
            Vec u(10, 0.), v(10, 0.); u[1] = kBRdir; v[0] = 1;
            Ratio(a, fg.C, fg.Cd, fg.Cs, u, v, R.est[qFGBRG], R.err[qFGBRG], R.errd[qFGBRG], R.errs[qFGBRG]); R.ok[qFGBRG] = true;
        }
    }
}

// double-bootstrap replica sets of the full sample (wp4_fill_boot.C): index p (parent) and P + p C + c (child c of p)
inline bool LoadBoot(const char* file, const Grid& G, const Tmpl& Tfull, int& P, int& C, std::vector<FitTmpl*>& out, std::vector<Tmpl*>& keep) {
    TFile* f = TFile::Open(file);
    if (!f || f->IsZombie()) return false;
    TVectorD* meta = (TVectorD*)f->Get("meta");
    if (!meta) return false;
    P = (int)(*meta)[0]; C = (int)(*meta)[1];
    const int nF = (int)(*meta)[2];
    const int nRep = P * (1 + C);
    const size_t sF10 = 10 * nF, sF5 = 5 * nF, stride = sF10 + sF5 + 5;
    const int tc[7] = {kTC2, kTC3, kTC4, kTC5, kTGeta, kTGDnom, kTGCnom};
    const int hc[7] = {kC2, kC3, kC4, kC5, kGeta, kGD, kGconv};
    keep.assign(nRep, nullptr); out.assign(nRep, nullptr);
    for (int r = 0; r < nRep; ++r) {
        keep[r] = new Tmpl();
        Tmpl& T = *keep[r];
        for (int c = 0; c < kNTComp; ++c) {
            for (int h = 0; h < kNh15; ++h) { T.R[c][h].assign(G.nR + 2, 0.); T.F15[c][h].assign(G.nF + 2, 0.); T.F15v[c][h].assign(G.nF + 2, 0.); }
            for (int k = 0; k < 10; ++k) { T.F10[c][k].assign(G.nF + 2, 0.); T.F10v[c][k].assign(G.nF + 2, 0.); }
            T.tot[c] = 0;
        }
    }
    double W[7][5];
    std::vector<std::array<std::array<double, 4>, kNSig>> Wr(nRep);
    for (int q = 0; q < 7; ++q) {
        TVectorD* v = (TVectorD*)f->Get(CompTag(hc[q]));
        if (!v) return false;
        const int c = tc[q];
        const bool gd = c == kTGDnom, gc = c == kTGCnom;
        for (int r = 0; r < nRep; ++r) {
            const double* a = v->GetMatrixArray() + (size_t)r * stride;
            Tmpl& T = *keep[r];
            for (int k = 0; k < 10; ++k)
                for (int x = 1; x <= G.nF; ++x) {
                    T.F10[c][k][x] = a[k * nF + x] * kNbToRaw;
                    T.F10v[c][k][x] = gd ? Tfull.F10v[kTGDom][k][x] + Tfull.F10v[kTGDnom][k][x]
                                     : gc ? Tfull.F10v[kTGCom][k][x] + Tfull.F10v[kTGCnom][k][x] : Tfull.F10v[c][k][x];
                }
            for (int j = 0; j < kNj; ++j) {
                const int h15 = H15(j, 2);
                for (int x = 1; x <= G.nF; ++x) {
                    T.F15[c][h15][x] = a[sF10 + j * nF + x] * kNbToRaw;
                    T.F15v[c][h15][x] = gd ? Tfull.F15v[kTGDom][h15][x] + Tfull.F15v[kTGDnom][h15][x]
                                       : gc ? Tfull.F15v[kTGCom][h15][x] + Tfull.F15v[kTGCnom][h15][x] : Tfull.F15v[c][h15][x];
                }
            }
            if (q < kNSig) {
                for (int w = 0; w < 5; ++w) W[q][w] = a[sF10 + sF5 + w] * kNbToRaw;
                Wr[r][q][0] = W[q][0]; Wr[r][q][1] = W[q][1]; Wr[r][q][2] = W[q][2] + W[q][3]; Wr[r][q][3] = W[q][4];
            }
        }
    }
    for (int r = 0; r < nRep; ++r) {
        out[r] = new FitTmpl();
        MakeFitTmpl(*keep[r], G, *out[r]);
        for (int t = 0; t < kNSig; ++t) for (int w = 0; w < 4; ++w) out[r]->W[t][w] = Wr[r][t][w];
    }
    f->Close();
    return true;
}

// ===========================================================================
// 9. scenarios and truth
// ===========================================================================
struct Scen {
    const char* tag; const char* text;
    double a, kap, rho, g;
    int bkg;            // 0: uniform x cbk; 1: omega/rho0 x 0.3, rest x 1.6
    double cbk, b;
    int mv;             // mass variant of the truth
    int cbTruth;        // 0: CB spread rule (= estimate); 1: per-cell geometric mean of the band (S1 variant)
};
const int kNScen = 12;
const Scen kScen[kNScen] = {
    {"a", "background x 1.35, tilt -0.15, a_eta 0.4, rho_d 1, kappa 1, g 1", 0.4, 1.0, 1.0, 1.0, 0, 1.35, -0.15, kMnom, 0},
    {"b-", "(a) with rho_d = 0.5", 0.4, 1.0, 0.5, 1.0, 0, 1.35, -0.15, kMnom, 0},
    {"b+", "(a) with rho_d = 1.5", 0.4, 1.0, 1.5, 1.0, 0, 1.35, -0.15, kMnom, 0},
    {"c-", "(a) with kappa = 0.8", 0.4, 0.8, 1.0, 1.0, 0, 1.35, -0.15, kMnom, 0},
    {"c+", "(a) with kappa = 1.2", 0.4, 1.2, 1.0, 1.0, 0, 1.35, -0.15, kMnom, 0},
    {"d-", "(a) with g = 0.5", 0.4, 1.0, 1.0, 0.5, 0, 1.35, -0.15, kMnom, 0},
    {"d+", "(a) with g = 1.5", 0.4, 1.0, 1.0, 1.5, 0, 1.35, -0.15, kMnom, 0},
    {"e", "(a) with omega/rho0-pion background x 0.3, the rest x 1.6 (tilt -0.15 kept)", 0.4, 1.0, 1.0, 1.0, 1, 1.35, -0.15, kMnom, 0},
    {"f", "(a) with truth mass transfer s + 0.0009 (0.99705) and r = 0.89", 0.4, 1.0, 1.0, 1.0, 0, 1.35, -0.15, kMf, 0},
    {"g", "(a) with a_eta = 0 (null)", 0.0, 1.0, 1.0, 1.0, 0, 1.35, -0.15, kMnom, 0},
    {"h", "(a) with a_eta = 1", 1.0, 1.0, 1.0, 1.0, 0, 1.35, -0.15, kMnom, 0},
    {"i", "(a) with the true CB distributed as the per-cell geometric mean (informational, S1 variant)", 0.4, 1.0, 1.0, 1.0, 0, 1.35, -0.15, kMnom, 1}};

struct Truth {
    Vec mu;            // expected toy counts [h15 * (nR + 2) + r]
    Vec cbTrue;        // the CB-like part of mu
    double val[kNQ];   // true value of every quantity
    double cbTot = 0, bkgTot = 0, etaTot = 0;
};
// truth for scenario s from the truth templates TT (half and mass variant chosen by the caller), truth smoothing TS
inline void MakeTruth(const Scen& s, const Tmpl& TT, const TruthSmooth& TS, const DataIn& Din, const Grid& G, Truth& U) {
    const int st = G.nR + 2;
    U.mu.assign(kNh15 * st, 0.); U.cbTrue.assign(kNh15 * st, 0.);
    U.cbTot = U.bkgTot = U.etaTot = 0;
    // geometric-mean fractions (scenario i): g(h | band) = 2 sqrt(N++_h N--_h) / 2 sqrt(N++ N--) of the band
    double gm[kNBand][kNh15];
    for (int b = 0; b < kNBand; ++b) {
        double P = 0, M = 0, pp[kNh15] = {0}, mm[kNh15] = {0};
        for (int i = b * 20; i < (b + 1) * 20 && i < kNfine; ++i)
            for (int h = 0; h < kNh15; ++h) { pp[h] += Din.ls.pp[i * kNh16 + h]; mm[h] += Din.ls.mm[i * kNh16 + h]; }
        for (int h = 0; h < kNh15; ++h) { P += pp[h]; M += mm[h]; }
        for (int h = 0; h < kNh15; ++h) gm[b][h] = (P > 0 && M > 0) ? std::sqrt(pp[h] * mm[h]) / std::sqrt(P * M) : 0.;
    }
    for (int h = 0; h < kNh15; ++h)
        for (int r = 1; r <= G.nR; ++r) {
            const double lo = G.ref[r - 1], hi = G.ref[r];
            const int band = Band(0.5 * (lo + hi));
            const double cb1 = CBWin(Din.cb, lo, hi);
            const double cb = cb1 * (s.cbTruth == 1 ? gm[band][h] : Din.cb.f15[band][h]);
            double i0, i1;
            TruthBkgCell(TS, h / kNe3, h % kNe3, lo, hi, s.bkg, i0, i1);
            const double bg = (s.bkg == 0 ? s.cbk : 1.0) * (i0 + s.b * i1);
            const double et = s.a * (TT.R[kTC2][h][r] + s.kap * (TT.R[kTC3][h][r] + TT.R[kTC4][h][r]) + s.rho * TT.R[kTC5][h][r] +
                                     s.g * TT.R[kTGeta][h][r]);
            U.mu[h * st + r] = cb + bg + et; U.cbTrue[h * st + r] = cb;
            U.cbTot += cb; U.bkgTot += bg; U.etaTot += et;
        }
    // true values
    auto W = [&](int c, int w) {
        if (w == 2) return WinSum(TT, G, c, kWin[kWSBlo].lo, kWin[kWSBlo].hi) + WinSum(TT, G, c, kWin[kWSBhi].lo, kWin[kWSBhi].hi);
        const int id = w == 0 ? kWB : w == 1 ? kWP : kWCRU;
        return WinSum(TT, G, c, kWin[id].lo, kWin[id].hi);
    };
    const double Yb = kYf * s.a * (W(kTC2, 0) + s.kap * W(kTC3, 0));
    const double Yp = kYf * s.a * (s.kap * W(kTC4, 1) + s.rho * W(kTC5, 1));
    for (int q = 0; q < kNQ; ++q) U.val[q] = 0;
    U.val[qFMa] = s.a; U.val[qFMkap] = s.kap; U.val[qFMrho] = s.rho; U.val[qFMg] = s.g;
    U.val[qFMYb] = Yb; U.val[qFMYp] = Yp;
    U.val[qFMY5] = kYf * s.a * s.rho * W(kTC5, 1); U.val[qFMY4] = kYf * s.a * s.kap * W(kTC4, 1);
    U.val[qFFa2] = s.a; U.val[qFFa3] = s.a * s.kap; U.val[qFFa4] = s.a * s.kap; U.val[qFFa5] = s.a * s.rho;
    U.val[qFFYb] = Yb; U.val[qFFYp] = Yp; U.val[qFFr5] = s.rho; U.val[qFFR34] = 1.0; U.val[qFFBRC] = kBRdir * s.rho;
    U.val[qFGa2] = s.a; U.val[qFGa5] = s.a * s.rho; U.val[qFGBRG] = kBRdir * s.rho;
    U.val[qB4Yb] = Yb; U.val[qPSBYp] = Yp;
}
// one toy (toy < 0: Asimov, the expectation itself; like-sign = data counts)
inline void MakeToy(const Truth& U, const DataIn& Din, const Grid& G, TRandom3* rnd, Obs& O) {
    const int st = G.nR + 2;
    O.n.assign(kNh15 * st, 0.);
    O.ls = Din.ls;
    if (rnd) {
        for (size_t i = 0; i < O.n.size(); ++i) O.n[i] = U.mu[i] > 0 ? rnd->Poisson(U.mu[i]) : 0.;
        for (size_t i = 0; i < O.ls.pp.size(); ++i) {
            O.ls.pp[i] = Din.ls.pp[i] > 0 ? rnd->Poisson(Din.ls.pp[i]) : 0.;
            O.ls.mm[i] = Din.ls.mm[i] > 0 ? rnd->Poisson(Din.ls.mm[i]) : 0.;
        }
    } else {
        O.n = U.mu;
    }
    Derive(G, O);
}

}  // namespace wp4

#endif
