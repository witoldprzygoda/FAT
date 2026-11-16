// ============================================================================
// wp5_lib.h  --  stage-B extraction of SPEC report_all (WP5, prefix wp5)
// ============================================================================
// Inputs: hists_v1.root (WP0), the event-group covariances (wp3_rows.root, wp5_egcov.root).
// The estimator is the one WP4 validated (wp4_lib.h): the same smoothing of the non-eta templates
// (wp4::SmoothFit, exp(pol3) in x = (M4 - 0.55)/0.1 with one common scale per histogram, parent
// shape for sparse F5g histograms), the same Fisher-scoring iteration (Poisson / IRLS fixed point,
// variance V = max(mu, 1) + k_e^2 var(CB) + sum coef^2 var(raw signal template) = BB-lite),
// the same tilt variable at the F10 bin centre. It is written here in a general form (terms per
// slice, fixed parameters, quadratic and per-slice tilts) for the stage-B configurations, and it is
// checked against wp4::NLFit on the configuration that WP4 used (F-main, E1 + E2).
//
// Model, bin i of histogram h = (slice j, M(ee) class e), RAW rows:
//   mu_i = k_e CB_i + c_j (1 + b x_i [+ q x_i^2]) [G_D + kappa_B G_C]_i (smoothed) + sum_t coef_t(theta, j) T_t,i + fix_i
//   T_t = C2, C3, C4, C5, G_eta (raw simulated rows x L x 63/64); coef from "terms" (theta[p1] x theta[p2], per slice).
// Errors (all on the same fixed point):
//   * C     : (J^T V^-1 J + P)^-1 (the fitter's own covariance, BB-lite);
//   * stat  : sandwich C (J^T V^-1 S V^-1 J + P) C, split into its parts:
//             N  = event-group covariance of the observed counts (scaled to the model, as WP3),
//             CB = covariance of the CB from the 200 like-sign replicas (k_e k_e' S_CB),
//             one part per Gaussian constraint (k_E2, g, kappa);
//   * MC    : covariance over 100 fits with the WP0 bootstrap replicas of the simulation (smoothing redone);
//   * total : stat + MC.
// ============================================================================
#ifndef WP5_LIB_H
#define WP5_LIB_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <functional>
#include <map>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TDecompChol.h"

#include "s5_blind_ana_fit.h"
#include "hists_io.h"
#include "wp4_lib.h"

namespace wp5 {
using namespace alldefs;
typedef std::vector<double> Vec;
typedef std::vector<Vec> Mat;

const char* const kOut = "/hdd3/przygoda/HADES/FAT_out/all_report/wp5";
const char* const kFigDir = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_all";
const char* const kTabDir = "/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all";
const char* const kEg3 = "/hdd3/przygoda/HADES/FAT_out/all_report/wp3/wp3_rows.root";
const char* const kEg5 = "/hdd3/przygoda/HADES/FAT_out/all_report/wp5/wp5_egcov.root";
const double kYf = kPT3;   // yields = raw rows x 64/63
const int kNB = kNF10;     // 70 bins per histogram
const int kNH = 15;        // 10 F10 histograms (h = 2 j + e) + 5 F5g (h = 10 + j)

enum SigT { kS2 = 0, kS3, kS4, kS5, kSeta, kNSig };
const int kSigComp[kNSig] = {kC2, kC3, kC4, kC5, kGeta};
inline const char* SigTag(int t) { static const char* n[kNSig] = {"C2", "C3", "C4", "C5", "Geta"}; return n[t]; }
// yield windows: bump, peak, side bands (both), CR-U, fit range, all M4
enum YWin { kYB = 0, kYP, kYSB, kYU, kYFit, kYTot, kNYW };
inline const char* YWTag(int w) { static const char* n[kNYW] = {"B", "P", "SB", "CRU", "Fit", "tot"}; return n[w]; }

inline Mat Zero(int n, int m) { return Mat(n, Vec(m, 0.)); }
inline double Quad(const Mat& C, const Vec& g) {
    double s = 0;
    for (size_t k = 0; k < g.size(); ++k) { if (g[k] == 0) continue; for (size_t l = 0; l < g.size(); ++l) s += g[k] * C[k][l] * g[l]; }
    return s;
}
inline Mat Add(const Mat& a, const Mat& b) { Mat c = a; for (size_t i = 0; i < a.size(); ++i) for (size_t j = 0; j < a[i].size(); ++j) c[i][j] += b[i][j]; return c; }
inline Mat Scale(const Mat& a, double f) { Mat c = a; for (auto& r : c) for (auto& x : r) x *= f; return c; }

// ===========================================================================
// 1. templates (simulation, expected raw data rows)
// ===========================================================================
struct Templ {
    Vec s[kNSig][kNH], vs[kNSig][kNH];       // raw signal templates per histogram, bins 1..70
    Vec gd[kNH], vgd[kNH], gc[kNH], vgc[kNH]; // raw non-eta templates (G_D, G_C)
    wp4::Smooth sgd[kNH], sgc[kNH];          // smoothing
    Vec gdI[kNH], gcI[kNH];                  // smoothed bin integrals
    double W[kNSig][kNYW];                   // window sums of the signal templates (all slices, all M(ee))
    double Wbj[kNSig][kNSlice];              // bump window per slice (F10 E1 + E2, M4 0.455-0.5225)
    double Wgd[kNYW], Wgc[kNYW];             // window sums of the raw non-eta templates
};
inline std::string HName(int h) { return h < 10 ? F10Name(h / 2, h % 2) : F5gName(h - 10); }
inline void Alloc(Templ& T) {
    for (int h = 0; h < kNH; ++h) {
        for (int t = 0; t < kNSig; ++t) { T.s[t][h].assign(kNB + 2, 0.); T.vs[t][h].assign(kNB + 2, 0.); }
        T.gd[h].assign(kNB + 2, 0.); T.vgd[h].assign(kNB + 2, 0.); T.gc[h].assign(kNB + 2, 0.); T.vgc[h].assign(kNB + 2, 0.);
        T.gdI[h].assign(kNB + 2, 0.); T.gcI[h].assign(kNB + 2, 0.);
    }
    for (int t = 0; t < kNSig; ++t) { for (int w = 0; w < kNYW; ++w) T.W[t][w] = 0; for (int j = 0; j < kNSlice; ++j) T.Wbj[t][j] = 0; }
    for (int w = 0; w < kNYW; ++w) T.Wgd[w] = T.Wgc[w] = 0;
}
inline void SmoothAll(Templ& T, int forceOrder = -1) {
    const std::vector<double> e = F10Edges();
    for (int h = 0; h < 10; ++h) { T.sgd[h] = wp4::SmoothFit(T.gd[h], T.vgd[h], e, nullptr, forceOrder); T.sgc[h] = wp4::SmoothFit(T.gc[h], T.vgc[h], e, nullptr, forceOrder); }
    for (int j = 0; j < kNSlice; ++j) {
        const int h = 10 + j, par = 2 * j + 1;
        T.sgd[h] = wp4::SmoothFit(T.gd[h], T.vgd[h], e, &T.sgd[par], forceOrder);
        T.sgc[h] = wp4::SmoothFit(T.gc[h], T.vgc[h], e, &T.sgc[par], forceOrder);
    }
    for (int h = 0; h < kNH; ++h)
        for (int x = 1; x <= kNB; ++x) { T.gdI[h][x] = T.sgd[h].Int(e[x - 1], e[x]); T.gcI[h][x] = T.sgc[h].Int(e[x - 1], e[x]); }
}
// raw (unsmoothed) non-eta templates as the "smoothed" ones (S4 variant "raw BB-lite")
inline void RawAsSmooth(Templ& T) { for (int h = 0; h < kNH; ++h) { T.gdI[h] = T.gd[h]; T.gcI[h] = T.gc[h]; } }
inline void WinFromF10(Templ& T) {
    const std::vector<double> e = F10Edges();
    for (int t = 0; t < kNSig; ++t)
        for (int j = 0; j < kNSlice; ++j) {
            double s = 0;
            for (int ec = 0; ec < 2; ++ec) for (int x = 1; x <= kNB; ++x) if (e[x - 1] >= kBlo - 1e-9 && e[x] <= kBhi + 1e-9) s += T.s[t][2 * j + ec][x];
            T.Wbj[t][j] = s;
        }
}
inline int WinBin(int w) {   // bin of the "win" histogram (WinId + 1); -1 for SB (two bins) and tot
    switch (w) { case kYB: return kWB + 1; case kYP: return kWP + 1; case kYU: return kWCRU + 1; case kYFit: return kWFit + 1; }
    return -1;
}
// full-sample templates of sel__var, scaled by f (luminosity fraction of a data subset; 1 for all data)
inline bool LoadTempl(TFile* fh, const char* sel, const char* var, double f, Templ& T, bool smooth = true) {
    Alloc(T);
    const double k = kNbToRaw * f;
    const std::string d = std::string("sim/") + sel + "__" + var + "/";
    for (int h = 0; h < kNH; ++h) {
        for (int t = 0; t < kNSig + 2; ++t) {
            const int comp = t < kNSig ? kSigComp[t] : (t == kNSig ? kGD : kGconv);
            TH1* x = (TH1*)fh->Get((d + CompTag(comp) + "_" + HName(h)).c_str());
            if (!x) { printf("wp5: %s%s_%s missing\n", d.c_str(), CompTag(comp), HName(h).c_str()); return false; }
            Vec& c = t < kNSig ? T.s[t][h] : (t == kNSig ? T.gd[h] : T.gc[h]);
            Vec& v = t < kNSig ? T.vs[t][h] : (t == kNSig ? T.vgd[h] : T.vgc[h]);
            for (int b = 1; b <= kNB; ++b) { c[b] = x->GetBinContent(b) * k; v[b] = std::pow(x->GetBinError(b) * k, 2); }
        }
    }
    for (int t = 0; t < kNSig + 2; ++t) {
        const int comp = t < kNSig ? kSigComp[t] : (t == kNSig ? kGD : kGconv);
        TH1* w = (TH1*)fh->Get((d + CompTag(comp) + "_win").c_str());
        TH1* tt = (TH1*)fh->Get((d + CompTag(comp) + "_tot").c_str());
        if (!w || !tt) return false;
        double* out = t < kNSig ? T.W[t] : (t == kNSig ? T.Wgd : T.Wgc);
        for (int q = 0; q < kNYW; ++q) {
            if (q == kYSB) out[q] = (w->GetBinContent(kWSBlo + 1) + w->GetBinContent(kWSBhi + 1)) * k;
            else if (q == kYTot) out[q] = tt->GetBinContent(1) * k;
            else out[q] = w->GetBinContent(WinBin(q)) * k;
        }
    }
    WinFromF10(T);
    if (smooth) SmoothAll(T);
    return true;
}
// the 100 WP0 bootstrap replicas of sel__nom; variances (BB-lite weights and smoothing scale) = those of the full sample
inline bool LoadReplicas(TFile* fh, const char* sel, double f, const Templ& Tfull, int nrep, std::vector<Templ>& R, int forceOrder = -1, bool raw = false) {
    R.assign(nrep, Templ());
    for (auto& T : R) Alloc(T);
    const double k = kNbToRaw * f;
    const std::string d = std::string("sim/") + sel + "__nom/";
    for (int t = 0; t < kNSig + 2; ++t) {
        const int comp = t < kNSig ? kSigComp[t] : (t == kNSig ? kGD : kGconv);
        for (int h = 0; h < kNH; ++h) {
            TH2* x = (TH2*)fh->Get((d + CompTag(comp) + "_" + HName(h) + "_rep").c_str());
            if (!x || x->GetNbinsY() < nrep) { printf("wp5: replica %s %s missing\n", CompTag(comp), HName(h).c_str()); return false; }
            for (int r = 0; r < nrep; ++r)
                for (int b = 1; b <= kNB; ++b) {
                    const double c = x->GetBinContent(b, r + 1) * k;
                    if (t < kNSig) { R[r].s[t][h][b] = c; R[r].vs[t][h][b] = Tfull.vs[t][h][b]; }
                    else if (t == kNSig) { R[r].gd[h][b] = c; R[r].vgd[h][b] = Tfull.vgd[h][b]; }
                    else { R[r].gc[h][b] = c; R[r].vgc[h][b] = Tfull.vgc[h][b]; }
                }
        }
        TH2* w = (TH2*)fh->Get((d + CompTag(comp) + "_win_rep").c_str());
        TH2* tt = (TH2*)fh->Get((d + CompTag(comp) + "_tot_rep").c_str());
        if (!w || !tt) return false;
        for (int r = 0; r < nrep; ++r) {
            double* out = t < kNSig ? R[r].W[t] : (t == kNSig ? R[r].Wgd : R[r].Wgc);
            for (int q = 0; q < kNYW; ++q) {
                if (q == kYSB) out[q] = (w->GetBinContent(kWSBlo + 1, r + 1) + w->GetBinContent(kWSBhi + 1, r + 1)) * k;
                else if (q == kYTot) out[q] = tt->GetBinContent(1, r + 1) * k;
                else out[q] = w->GetBinContent(WinBin(q), r + 1) * k;
            }
        }
    }
    for (auto& T : R) { WinFromF10(T); if (raw) RawAsSmooth(T); else SmoothAll(T, forceOrder); }
    return true;
}

// ===========================================================================
// 2. data (raw rows): N, CB, CB covariance (200 replicas), event-group covariance of N
// ===========================================================================
struct DataSet {
    std::string sel, sub;
    Vec y[kNH], cb[kNH], vcb[kNH];
    std::vector<double> covCB10, covCB5;      // dense 700 x 700 and 350 x 350 (cell = h x 70 + x - 1)
    std::vector<double> eg10, eg5;            // event-group covariance of N
    bool egFull = false;
    double Nwin[kNYW], CBwin[kNYW];
};
inline int Cell(int h, int x) { return (h < 10 ? h : h - 10) * kNB + (x - 1); }
// CB covariance from the replicas of sel/all; for a subset, cov_ab x sqrt(r_a r_b), r = CB_sub/CB_all per bin
inline bool LoadData(TFile* fh, const char* sel, const char* sub, DataSet& D, FILE* log = nullptr) {
    D.sel = sel; D.sub = sub;
    const bool all = std::string(sub) == "all";
    Vec cbAll[kNH];
    for (int h = 0; h < kNH; ++h) {
        TH1* n = allio::Data(fh, sel, sub, "N", HName(h));
        TH1* c = allio::Data(fh, sel, sub, "CB", HName(h));
        TH1* ca = allio::Data(fh, sel, "all", "CB", HName(h));
        if (!n || !c || !ca) return false;
        D.y[h].assign(kNB + 2, 0.); D.cb[h].assign(kNB + 2, 0.); D.vcb[h].assign(kNB + 2, 0.); cbAll[h].assign(kNB + 2, 0.);
        for (int x = 1; x <= kNB; ++x) { D.y[h][x] = n->GetBinContent(x); D.cb[h][x] = c->GetBinContent(x); cbAll[h][x] = ca->GetBinContent(x); }
        delete n; delete c; delete ca;
    }
    // CB replica covariance (dense) per layout
    for (int lay = 0; lay < 2; ++lay) {
        const int h0 = lay == 0 ? 0 : 10, nh = lay == 0 ? 10 : 5, nc = nh * kNB;
        std::vector<double>& C = lay == 0 ? D.covCB10 : D.covCB5;
        C.assign((size_t)nc * nc, 0.);
        std::vector<Vec> rep(nc);
        int nr = 0;
        for (int q = 0; q < nh; ++q) {
            TH2* r = (TH2*)fh->Get(Form("data/%s/all/CBrep_%s", sel, HName(h0 + q).c_str()));
            if (!r) { printf("wp5: CBrep %s missing\n", HName(h0 + q).c_str()); return false; }
            nr = r->GetNbinsY();
            for (int x = 1; x <= kNB; ++x) { Vec& v = rep[q * kNB + x - 1]; v.resize(nr); for (int k = 0; k < nr; ++k) v[k] = r->GetBinContent(x, k + 1); }
        }
        Vec m(nc, 0.);
        for (int a = 0; a < nc; ++a) { for (int k = 0; k < nr; ++k) m[a] += rep[a][k]; m[a] /= nr; }
        for (int a = 0; a < nc; ++a)
            for (int b = a; b < nc; ++b) {
                double s = 0;
                for (int k = 0; k < nr; ++k) s += (rep[a][k] - m[a]) * (rep[b][k] - m[b]);
                s /= nr;
                C[(size_t)a * nc + b] = C[(size_t)b * nc + a] = s;
            }
        if (!all) {
            Vec rr(nc, 0.);
            for (int q = 0; q < nh; ++q) for (int x = 1; x <= kNB; ++x) { const double A = cbAll[h0 + q][x]; rr[q * kNB + x - 1] = A > 0 ? D.cb[h0 + q][x] / A : 0.; }
            for (int a = 0; a < nc; ++a) for (int b = 0; b < nc; ++b) C[(size_t)a * nc + b] *= std::sqrt(std::max(0., rr[a] * rr[b]));
        }
        for (int q = 0; q < nh; ++q) for (int x = 1; x <= kNB; ++x) { const int a = q * kNB + x - 1; D.vcb[h0 + q][x] = C[(size_t)a * nc + a]; }
    }
    // event-group covariance of N
    D.egFull = false;
    {
        TFile* f5 = TFile::Open(kEg5);
        TMatrixDSym* m10 = f5 ? (TMatrixDSym*)f5->Get(Form("egcov_%s_%s_F10", sel, sub)) : nullptr;
        TMatrixDSym* m5 = f5 ? (TMatrixDSym*)f5->Get(Form("egcov_%s_%s_F5g", sel, sub)) : nullptr;
        if (m10 && m5) {
            D.eg10.assign(700 * 700, 0.); D.eg5.assign(350 * 350, 0.);
            for (int a = 0; a < 700; ++a) for (int b = 0; b < 700; ++b) D.eg10[(size_t)a * 700 + b] = (*m10)(a, b);
            for (int a = 0; a < 350; ++a) for (int b = 0; b < 350; ++b) D.eg5[(size_t)a * 350 + b] = (*m5)(a, b);
            D.egFull = true;
        }
        if (f5) f5->Close();
        if (!D.egFull) {   // diagonal from the bin errors (event-group errors of hists_v1)
            D.eg10.assign(700 * 700, 0.); D.eg5.assign(350 * 350, 0.);
            for (int h = 0; h < kNH; ++h) {
                TH1* n = allio::Data(fh, sel, sub, "N", HName(h));
                for (int x = 1; x <= kNB; ++x) {
                    const int a = Cell(h, x);
                    if (h < 10) D.eg10[(size_t)a * 700 + a] = std::pow(n->GetBinError(x), 2);
                    else D.eg5[(size_t)a * 350 + a] = std::pow(n->GetBinError(x), 2);
                }
                delete n;
            }
            if (log) fprintf(log, "  [data %s/%s: no event-group covariance file entry -> diagonal event-group variances]\n", sel, sub);
        }
    }
    TH1* nw = allio::Data(fh, sel, sub, "N", "win");
    TH1* cw = allio::Data(fh, sel, sub, "CB", "win");
    TH1* nt = allio::Data(fh, sel, sub, "N", "tot");
    TH1* ct = allio::Data(fh, sel, sub, "CB", "tot");
    for (int q = 0; q < kNYW; ++q) {
        if (q == kYSB) { D.Nwin[q] = nw->GetBinContent(kWSBlo + 1) + nw->GetBinContent(kWSBhi + 1); D.CBwin[q] = cw->GetBinContent(kWSBlo + 1) + cw->GetBinContent(kWSBhi + 1); }
        else if (q == kYTot) { D.Nwin[q] = nt->GetBinContent(1); D.CBwin[q] = ct->GetBinContent(1); }
        else { D.Nwin[q] = nw->GetBinContent(WinBin(q)); D.CBwin[q] = cw->GetBinContent(WinBin(q)); }
    }
    delete nw; delete cw; delete nt; delete ct;
    return true;
}

// ===========================================================================
// 3. fit bins and the model
// ===========================================================================
struct Bins {
    int n = 0, lay = 0;             // lay 0: F10 (cells of 700), 1: F5g (cells of 350)
    Vec y, cb, vcb, x, gd, gc, fix, lo, hi;
    std::vector<int> j, e, h, xb, cell;
    Vec s[kNSig], vs[kNSig];
};
// histograms used: F10 with the classes in useE (bit 1 = E1, bit 2 = E2), or F5g; range [mlo, mhi)
inline void MakeBins(const DataSet& D, const Templ& T, int lay, int useE, Bins& B, double mlo = kF10Lo, double mhi = kF10Hi,
                     const std::function<double(int, int)>& fixfn = nullptr) {
    B = Bins(); B.lay = lay;
    const std::vector<double> e = F10Edges();
    const int h0 = lay == 0 ? 0 : 10, nh = lay == 0 ? 10 : 5;
    for (int q = 0; q < nh; ++q) {
        const int h = h0 + q;
        const int jj = lay == 0 ? q / 2 : q, ec = lay == 0 ? q % 2 : 1;
        if (lay == 0 && !((useE >> ec) & 1)) continue;
        for (int x = 1; x <= kNB; ++x) {
            if (e[x - 1] < mlo - 1e-9 || e[x] > mhi + 1e-9) continue;
            B.y.push_back(D.y[h][x]); B.cb.push_back(D.cb[h][x]); B.vcb.push_back(D.vcb[h][x]);
            B.x.push_back(TiltX(0.5 * (e[x - 1] + e[x]))); B.lo.push_back(e[x - 1]); B.hi.push_back(e[x]);
            B.gd.push_back(T.gdI[h][x]); B.gc.push_back(T.gcI[h][x]);
            B.fix.push_back(fixfn ? fixfn(h, x) : 0.);
            B.j.push_back(jj); B.e.push_back(ec); B.h.push_back(h); B.xb.push_back(x); B.cell.push_back(Cell(h, x));
            for (int t = 0; t < kNSig; ++t) { B.s[t].push_back(T.s[t][h][x]); B.vs[t].push_back(T.vs[t][h][x]); }
            ++B.n;
        }
    }
}
// replace the template parts of B by those of T (same bins)
inline void SwapTempl(Bins& B, const Templ& T) {
    for (int i = 0; i < B.n; ++i) {
        const int h = B.h[i], x = B.xb[i];
        B.gd[i] = T.gdI[h][x]; B.gc[i] = T.gcI[h][x];
        for (int t = 0; t < kNSig; ++t) { B.s[t][i] = T.s[t][h][x]; B.vs[t][i] = T.vs[t][h][x]; }
    }
}

struct Term { int t, p1, p2, slice; };   // coef = theta[p1] (x theta[p2]) on bins of slice (all if < 0)
struct Spec {
    std::vector<std::string> names;
    Vec start;
    std::vector<int> fixed;              // 1 = held at start
    std::vector<bla::Constraint> cons;
    int iK[2] = {-1, -1}; double kFix[2] = {1, 1};
    int iC[kNSlice] = {-1, -1, -1, -1, -1};
    int iB = -1; double bFix = 0;
    int iBj[kNSlice] = {-1, -1, -1, -1, -1};
    int iQ = -1;
    int iKB = -1; double kBfix = 1;
    std::vector<Term> terms;
    int np() const { return (int)names.size(); }
    int Idx(const std::string& n) const { for (int k = 0; k < np(); ++k) if (names[k] == n) return k; return -1; }
    int Add(const std::string& n, double st) { names.push_back(n); start.push_back(st); fixed.push_back(0); return np() - 1; }
};
inline double Coef(const Spec& S, const Term& T, const Vec& th) { return th[T.p1] * (T.p2 >= 0 ? th[T.p2] : 1.); }
// mu, Jacobian, V (BB-lite) and Vd (data part) at theta
inline void Eval(const Bins& B, const Spec& S, const Vec& th, Vec& mu, Mat* J, Vec* V, Vec* Vd) {
    const int n = B.n, np = S.np();
    mu.assign(n, 0.);
    if (J) J->assign(n, Vec(np, 0.));
    if (V) V->assign(n, 0.);
    if (Vd) Vd->assign(n, 0.);
    const double kB = S.iKB >= 0 ? th[S.iKB] : S.kBfix;
    for (int i = 0; i < n; ++i) {
        const int e = B.e[i], j = B.j[i];
        const double k = S.iK[e] >= 0 ? th[S.iK[e]] : S.kFix[e];
        const double c = S.iC[j] >= 0 ? th[S.iC[j]] : 1.;
        const double b = S.iBj[j] >= 0 ? th[S.iBj[j]] : (S.iB >= 0 ? th[S.iB] : S.bFix);
        const double q = S.iQ >= 0 ? th[S.iQ] : 0.;
        const double bg = B.gd[i] + kB * B.gc[i];
        const double x = B.x[i];
        const double tl = 1. + b * x + q * x * x;
        double m = k * B.cb[i] + c * tl * bg + B.fix[i], vs = 0;
        for (const auto& T : S.terms) {
            if (T.slice >= 0 && T.slice != j) continue;
            const double c2 = T.p2 >= 0 ? th[T.p2] : 1.;
            const double co = th[T.p1] * c2;
            m += co * B.s[T.t][i];
            vs += co * co * B.vs[T.t][i];
            if (J) { (*J)[i][T.p1] += c2 * B.s[T.t][i]; if (T.p2 >= 0) (*J)[i][T.p2] += th[T.p1] * B.s[T.t][i]; }
        }
        mu[i] = m;
        if (J) {
            if (S.iK[e] >= 0) (*J)[i][S.iK[e]] = B.cb[i];
            if (S.iC[j] >= 0) (*J)[i][S.iC[j]] = tl * bg;
            if (S.iBj[j] >= 0) (*J)[i][S.iBj[j]] = c * x * bg;
            else if (S.iB >= 0) (*J)[i][S.iB] = c * x * bg;
            if (S.iQ >= 0) (*J)[i][S.iQ] = c * x * x * bg;
            if (S.iKB >= 0) (*J)[i][S.iKB] = c * tl * B.gc[i];
        }
        const double vb = std::max(m, 1.0) + k * k * B.vcb[i];
        if (V) (*V)[i] = vb + vs;
        if (Vd) (*Vd)[i] = vb;
    }
}
// parts of a model at theta (for figures): CB, background, each signal class
inline void Parts(const Bins& B, const Spec& S, const Vec& th, Vec& pcb, Vec& pbg, Vec* psig /*[kNSig]*/) {
    const int n = B.n;
    pcb.assign(n, 0.); pbg.assign(n, 0.);
    for (int t = 0; t < kNSig; ++t) psig[t].assign(n, 0.);
    const double kB = S.iKB >= 0 ? th[S.iKB] : S.kBfix;
    for (int i = 0; i < n; ++i) {
        const int e = B.e[i], j = B.j[i];
        const double k = S.iK[e] >= 0 ? th[S.iK[e]] : S.kFix[e];
        const double c = S.iC[j] >= 0 ? th[S.iC[j]] : 1.;
        const double b = S.iBj[j] >= 0 ? th[S.iBj[j]] : (S.iB >= 0 ? th[S.iB] : S.bFix);
        const double q = S.iQ >= 0 ? th[S.iQ] : 0.;
        pcb[i] = k * B.cb[i];
        pbg[i] = c * (1. + b * B.x[i] + q * B.x[i] * B.x[i]) * (B.gd[i] + kB * B.gc[i]);
        for (const auto& T : S.terms) { if (T.slice >= 0 && T.slice != j) continue; psig[T.t][i] += Coef(S, T, th) * B.s[T.t][i]; }
        // the fixed term (F-gold: C3 + C4 from F-main) is counted with C3 for display
        psig[kS3][i] += B.fix[i];
    }
}

struct FitRes {
    bool ok = false;
    int nit = 0, n = 0, nfree = 0, ndf = 0;
    Vec th;
    Mat C;                    // fitter covariance (BB-lite)
    Mat Cn, Ccb;              // sandwich parts: event-group N, CB replicas
    std::vector<Mat> Ccons;   // sandwich part of each constraint
    Mat Cstat;                // Cn + Ccb + sum Ccons
    Vec mu, V, Vd;
    double chi2BC = 0, chi2P = 0, chi2E = 0;   // Baker-Cousins deviance, Pearson with V, Pearson with V x event-group factor
};
inline Mat CovStat(const FitRes& r) { return r.Cstat; }

// Fisher scoring (Gauss-Newton IRLS); covariances; sandwich with the full covariances of N and CB
// Vfix: if given, the variance is frozen to it (fixed-weight chi2 = sum (y - mu)^2 / Vfix: its minimum is the fitter's fixed point
// when Vfix = V of that fit, so it is the consistent metric for profile scans)
inline FitRes Fit(const Bins& B, const Spec& S, const DataSet& D, const Vec* start = nullptr, bool sandwich = true, int maxIt = 300, const Vec* Vfix = nullptr) {
    FitRes R;
    const int n = B.n, np = S.np();
    R.n = n;
    std::vector<int> fr;
    for (int k = 0; k < np; ++k) if (!S.fixed[k]) fr.push_back(k);
    const int nf = (int)fr.size();
    R.nfree = nf;
    Vec th = start ? *start : S.start;
    for (int k = 0; k < np; ++k) if (S.fixed[k]) th[k] = S.start[k];
    Vec mu, V;
    Mat J;
    // constraint precision per free index
    Vec Pd(nf, 0.), P0(nf, 0.);
    for (const auto& c : S.cons) for (int a = 0; a < nf; ++a) if (fr[a] == c.k) { Pd[a] += 1. / (c.sig * c.sig); P0[a] = c.val; }
    bool conv = false;
    Mat Ainv;
    for (int it = 0; it < maxIt; ++it) {
        Eval(B, S, th, mu, &J, &V, nullptr);
        if (Vfix) V = *Vfix;
        Mat A = Zero(nf, nf); Vec b(nf, 0.);
        for (int i = 0; i < n; ++i) {
            const double w = 1. / V[i], r = B.y[i] - mu[i];
            for (int a = 0; a < nf; ++a) {
                const double ja = J[i][fr[a]] * w;
                if (ja == 0) continue;
                b[a] += ja * r;
                for (int c = 0; c <= a; ++c) A[a][c] += ja * J[i][fr[c]];
            }
        }
        for (int a = 0; a < nf; ++a) { A[a][a] += Pd[a]; b[a] -= Pd[a] * (th[fr[a]] - P0[a]); }
        for (int a = 0; a < nf; ++a) for (int c = a + 1; c < nf; ++c) A[a][c] = A[c][a];
        Vec sc(nf, 1.);
        for (int a = 0; a < nf; ++a) sc[a] = A[a][a] > 0 ? 1. / std::sqrt(A[a][a]) : 1.;
        Mat As = A;
        for (int a = 0; a < nf; ++a) for (int c = 0; c < nf; ++c) As[a][c] *= sc[a] * sc[c];
        if (!bla::CholInv(As)) { R.nit = it; return R; }
        for (int a = 0; a < nf; ++a) for (int c = 0; c < nf; ++c) As[a][c] *= sc[a] * sc[c];
        Ainv = As;
        Vec d(nf, 0.);
        for (int a = 0; a < nf; ++a) for (int c = 0; c < nf; ++c) d[a] += Ainv[a][c] * b[c];
        double dmax = 0, lim = 1.;
        for (int a = 0; a < nf; ++a) {
            const double z = std::fabs(d[a]) / std::sqrt(std::max(Ainv[a][a], 1e-300));
            dmax = std::max(dmax, z);
            if (z > 3.) lim = std::min(lim, 3. / z);
        }
        for (int a = 0; a < nf; ++a) th[fr[a]] += lim * d[a];
        R.nit = it + 1;
        if (it > 0 && dmax < 1e-5) { conv = true; break; }
    }
    if (!conv) return R;
    Vec Vd;
    Eval(B, S, th, mu, &J, &V, &Vd);
    if (Vfix) V = *Vfix;
    // fitter covariance
    Mat A = Zero(nf, nf);
    for (int i = 0; i < n; ++i) for (int a = 0; a < nf; ++a) { const double ja = J[i][fr[a]] / V[i]; if (ja == 0) continue; for (int c = 0; c < nf; ++c) A[a][c] += ja * J[i][fr[c]]; }
    for (int a = 0; a < nf; ++a) A[a][a] += Pd[a];
    Vec sc(nf, 1.);
    for (int a = 0; a < nf; ++a) sc[a] = A[a][a] > 0 ? 1. / std::sqrt(A[a][a]) : 1.;
    for (int a = 0; a < nf; ++a) for (int c = 0; c < nf; ++c) A[a][c] *= sc[a] * sc[c];
    if (!bla::CholInv(A)) return R;
    for (int a = 0; a < nf; ++a) for (int c = 0; c < nf; ++c) A[a][c] *= sc[a] * sc[c];
    auto embed = [&](const Mat& m) { Mat o = Zero(np, np); for (int a = 0; a < nf; ++a) for (int c = 0; c < nf; ++c) o[fr[a]][fr[c]] = m[a][c]; return o; };
    R.C = embed(A);
    // sandwich parts
    auto sand = [&](const Mat& M) {
        Mat T = Zero(nf, nf), O = Zero(nf, nf);
        for (int a = 0; a < nf; ++a) for (int c = 0; c < nf; ++c) { double s = 0; for (int k = 0; k < nf; ++k) s += A[a][k] * M[k][c]; T[a][c] = s; }
        for (int a = 0; a < nf; ++a) for (int c = 0; c < nf; ++c) { double s = 0; for (int k = 0; k < nf; ++k) s += T[a][k] * A[k][c]; O[a][c] = s; }
        return embed(O);
    };
    if (sandwich) {
        // u_i = J_i / V_i ; meat_N = sum_ab u_a S^N_ab u_b, S^N_ab = eg_ab sqrt(mu_a mu_b / (y_a y_b)); meat_CB = sum_ab u_a k_a k_b S^CB_ab u_b
        const int nc = B.lay == 0 ? 700 : 350;
        const std::vector<double>& eg = B.lay == 0 ? D.eg10 : D.eg5;
        const std::vector<double>& cc = B.lay == 0 ? D.covCB10 : D.covCB5;
        Mat U = Zero(n, nf);
        Vec sf(n, 0.), kk(n, 0.);
        for (int i = 0; i < n; ++i) {
            for (int a = 0; a < nf; ++a) U[i][a] = J[i][fr[a]] / V[i];
            sf[i] = B.y[i] > 0 ? std::sqrt(std::max(mu[i], 1e-9) / B.y[i]) : 0.;
            const int e = B.e[i];
            kk[i] = S.iK[e] >= 0 ? th[S.iK[e]] : S.kFix[e];
        }
        Mat Mn = Zero(nf, nf), Mc = Zero(nf, nf);
        Vec tn(nf), tc(nf);
        for (int i = 0; i < n; ++i) {
            std::fill(tn.begin(), tn.end(), 0.); std::fill(tc.begin(), tc.end(), 0.);
            const size_t ri = (size_t)B.cell[i] * nc;
            for (int l = 0; l < n; ++l) {
                double sn = eg[ri + B.cell[l]];
                if (B.y[i] > 0 && B.y[l] > 0) sn *= sf[i] * sf[l];
                else sn = (i == l) ? std::max(mu[i], 0.) : 0.;
                const double scb = cc[ri + B.cell[l]] * kk[i] * kk[l];
                if (sn == 0 && scb == 0) continue;
                for (int a = 0; a < nf; ++a) { tn[a] += sn * U[l][a]; tc[a] += scb * U[l][a]; }
            }
            for (int a = 0; a < nf; ++a) for (int c = 0; c < nf; ++c) { Mn[a][c] += U[i][a] * tn[c]; Mc[a][c] += U[i][a] * tc[c]; }
        }
        R.Cn = sand(Mn); R.Ccb = sand(Mc);
        R.Cstat = Add(R.Cn, R.Ccb);
        R.Ccons.clear();
        for (const auto& c : S.cons) {
            Mat Mp = Zero(nf, nf);
            for (int a = 0; a < nf; ++a) if (fr[a] == c.k) Mp[a][a] = 1. / (c.sig * c.sig);
            R.Ccons.push_back(sand(Mp));
            R.Cstat = Add(R.Cstat, R.Ccons.back());
        }
    } else {
        R.Cn = R.C; R.Ccb = Zero(np, np); R.Cstat = R.C;
        for (size_t q = 0; q < S.cons.size(); ++q) R.Ccons.push_back(Zero(np, np));
    }
    R.th = th; R.mu = mu; R.V = V; R.Vd = Vd;
    R.chi2BC = R.chi2P = R.chi2E = 0;
    const std::vector<double>& eg = B.lay == 0 ? D.eg10 : D.eg5;
    const int nc = B.lay == 0 ? 700 : 350;
    for (int i = 0; i < n; ++i) {
        const double m = std::max(mu[i], 1e-9), y = B.y[i];
        R.chi2BC += 2. * (m - y + (y > 0 ? y * std::log(y / m) : 0.));
        R.chi2P += (y - mu[i]) * (y - mu[i]) / V[i];
        const double se = y > 0 ? eg[(size_t)B.cell[i] * nc + B.cell[i]] / y : 1.;
        R.chi2E += (y - mu[i]) * (y - mu[i]) / (V[i] + (se - 1.) * std::max(mu[i], 1.));
    }
    for (const auto& c : S.cons) { const double z = (th[c.k] - c.val) / c.sig; R.chi2BC += z * z; R.chi2P += z * z; R.chi2E += z * z; }
    R.ndf = n - nf + (int)S.cons.size();
    R.ok = true;
    return R;
}

// ===========================================================================
// 4. specs
// ===========================================================================
struct Opt {
    int useE = 2;                // bit mask: 1 = E1, 2 = E2 (F10 fits)
    bool kappaCons = true;       // kappa = 1 +- 0.2 (fall-back "drop E1")
    double kappaSig = 0.2;
    bool quad = false, perSliceTilt = false;
    bool fixB = false; double bVal = 0, bSig = 0;   // b fixed (bSig = 0) or constrained (bSig > 0)
    bool fixKB = false; double kBval = 1;
    bool fixKE2 = false; double kE2val = 1;
    double gSig = 0.30;
    bool gFree = false;
};
// common background part; returns nothing, adds parameters
inline void AddBkg(Spec& S, const Opt& o, bool useE1, bool useE2) {
    for (int j = 0; j < kNSlice; ++j) S.iC[j] = S.Add(Form("c_%d", j + 1), 1.0);
    if (o.perSliceTilt) { for (int j = 0; j < kNSlice; ++j) S.iBj[j] = S.Add(Form("b_%d", j + 1), -0.07); }
    else {
        S.iB = S.Add("b", o.fixB ? o.bVal : -0.07);
        if (o.fixB && o.bSig <= 0) S.fixed[S.iB] = 1;
        if (o.fixB && o.bSig > 0) S.cons.push_back({S.iB, o.bVal, o.bSig});
    }
    if (o.quad) S.iQ = S.Add("q", 0.0);
    S.iKB = S.Add("kappa_B", o.fixKB ? o.kBval : 2.0);
    if (o.fixKB) S.fixed[S.iKB] = 1;
    if (useE1) S.iK[0] = S.Add("k_E1", 1.0);
    if (useE2) {
        S.iK[1] = S.Add("k_E2", o.fixKE2 ? o.kE2val : 1.0);
        if (o.fixKE2) S.fixed[S.iK[1]] = 1; else S.cons.push_back({S.iK[1], 1.0, 0.05});
    }
}
// F-main: a_eta, kappa, rho_d, g (+ background)
inline Spec SpecMain(const Opt& o, bool perSliceA = false) {
    Spec S;
    std::vector<int> ia;
    if (perSliceA) for (int j = 0; j < kNSlice; ++j) ia.push_back(S.Add(Form("a_eta_%d", j + 1), 0.4));
    else ia.push_back(S.Add("a_eta", 0.4));
    const int ik = S.Add("kappa", 1.0), ir = S.Add("rho_d", 1.0), ig = S.Add("g", 1.0);
    if (o.kappaCons) S.cons.push_back({ik, 1.0, o.kappaSig});
    if (!o.gFree) S.cons.push_back({ig, 1.0, o.gSig});
    for (size_t q = 0; q < ia.size(); ++q) {
        const int sl = perSliceA ? (int)q : -1;
        S.terms.push_back({kS2, ia[q], -1, sl});
        S.terms.push_back({kS3, ia[q], ik, sl});
        S.terms.push_back({kS4, ia[q], ik, sl});
        S.terms.push_back({kS5, ia[q], ir, sl});
        S.terms.push_back({kSeta, ia[q], ig, sl});
    }
    AddBkg(S, o, (o.useE & 1) != 0, (o.useE & 2) != 0);
    return S;
}
// F-free: alpha_2..alpha_5, g (G_eta x alpha_2 g)
inline Spec SpecFree(const Opt& o) {
    Spec S;
    const int i2 = S.Add("alpha_2", 0.4), i3 = S.Add("alpha_3", 0.4), i4 = S.Add("alpha_4", 0.4), i5 = S.Add("alpha_5", 0.4), ig = S.Add("g", 1.0);
    if (!o.gFree) S.cons.push_back({ig, 1.0, o.gSig});
    S.terms = {{kS2, i2, -1, -1}, {kS3, i3, -1, -1}, {kS4, i4, -1, -1}, {kS5, i5, -1, -1}, {kSeta, i2, ig, -1}};
    AddBkg(S, o, (o.useE & 1) != 0, (o.useE & 2) != 0);
    return S;
}
// F-gold (F5g): alpha_2, alpha_5, g; kappa_B fixed (from F-main); C3 + C4 in the fixed term
inline Spec SpecGold(const Opt& o, double kB) {
    Spec S;
    const int i2 = S.Add("alpha_2", 0.4), i5 = S.Add("alpha_5", 0.4), ig = S.Add("g", 1.0);
    if (!o.gFree) S.cons.push_back({ig, 1.0, o.gSig});
    S.terms = {{kS2, i2, -1, -1}, {kS5, i5, -1, -1}, {kSeta, i2, ig, -1}};
    Opt oo = o; oo.fixKB = true; oo.kBval = kB;
    AddBkg(S, oo, false, true);
    return S;
}

// ===========================================================================
// 5. derived quantities, replicas, errors
// ===========================================================================
struct Deriv {
    std::string name, text;
    std::function<double(const Vec&, const Templ&)> f;
};
struct DVal {
    std::string name, text;
    double v = 0, eN = 0, eCB = 0, eStat = 0, eMC = 0, eTot = 0, shift = 0, bc = 0;
    double eMC68 = 0, median = 0;   // robust MC spread (half of the central 68 % interval) and replica median
    Vec eCons;   // per constraint
    Vec rep;     // replica values
};
inline Vec NumGrad(const std::function<double(const Vec&)>& f, const Vec& th, const Mat& C) {
    Vec g(th.size(), 0.);
    for (size_t k = 0; k < th.size(); ++k) {
        const double h = std::max(1e-7, 1e-4 * std::sqrt(std::max(C[k][k], 0.)));
        if (C[k][k] <= 0) continue;
        Vec a = th, b = th; a[k] += h; b[k] -= h;
        g[k] = (f(a) - f(b)) / (2 * h);
    }
    return g;
}
struct Replicas {
    std::vector<FitRes> fit;
    Mat cov;          // parameter covariance over replicas
    Vec mean;
    int nok = 0;
};
// fit every replica (same data, same spec; start = nominal)
inline Replicas RunReplicas(const Bins& B0, const Spec& S, const DataSet& D, const std::vector<Templ>& R, const Vec& th0,
                            const std::function<void(Bins&, const Templ&)>& rebuild = nullptr) {
    Replicas out;
    const int np = S.np();
    out.fit.resize(R.size());
    Bins B = B0;
    for (size_t r = 0; r < R.size(); ++r) {
        if (rebuild) rebuild(B, R[r]); else SwapTempl(B, R[r]);
        out.fit[r] = Fit(B, S, D, &th0, false);
    }
    out.mean.assign(np, 0.); out.cov = Zero(np, np);
    for (auto& f : out.fit) if (f.ok) { ++out.nok; for (int k = 0; k < np; ++k) out.mean[k] += f.th[k]; }
    for (int k = 0; k < np; ++k) out.mean[k] /= std::max(1, out.nok);
    for (auto& f : out.fit) if (f.ok) for (int k = 0; k < np; ++k) for (int l = 0; l < np; ++l) out.cov[k][l] += (f.th[k] - out.mean[k]) * (f.th[l] - out.mean[l]);
    for (int k = 0; k < np; ++k) for (int l = 0; l < np; ++l) out.cov[k][l] /= std::max(1, out.nok - 1);
    return out;
}
inline DVal Evaluate(const Deriv& d, const FitRes& F, const Templ& T, const Replicas* rep, const std::vector<Templ>* R) {
    DVal o; o.name = d.name; o.text = d.text;
    auto fth = [&](const Vec& th) { return d.f(th, T); };
    o.v = fth(F.th);
    const Vec g = NumGrad(fth, F.th, F.C);
    o.eN = std::sqrt(std::max(0., Quad(F.Cn, g)));
    o.eCB = std::sqrt(std::max(0., Quad(F.Ccb, g)));
    for (const auto& C : F.Ccons) o.eCons.push_back(std::sqrt(std::max(0., Quad(C, g))));
    o.eStat = std::sqrt(std::max(0., Quad(F.Cstat, g)));
    if (rep && R) {
        double s = 0, q = 0; int n = 0;
        for (size_t r = 0; r < R->size(); ++r) {
            if (!rep->fit[r].ok) continue;
            const double v = d.f(rep->fit[r].th, (*R)[r]);
            o.rep.push_back(v); s += v; q += v * v; ++n;
        }
        if (n > 1) {
            const double m = s / n;
            o.eMC = std::sqrt(std::max(0., (q - n * m * m) / (n - 1)));
            o.shift = m - o.v;
            o.bc = 2 * o.v - m;
            Vec srt = o.rep; std::sort(srt.begin(), srt.end());
            auto qt = [&](double f) { const double x = f * (n - 1); const int i = (int)x; const double w = x - i; return i + 1 < n ? (1 - w) * srt[i] + w * srt[i + 1] : srt[n - 1]; };
            o.eMC68 = 0.5 * (qt(0.8413) - qt(0.1587));
            o.median = qt(0.5);
        }
    }
    o.eTot = std::sqrt(o.eStat * o.eStat + o.eMC * o.eMC);
    return o;
}

// standard derived quantities of a F-main type spec (a, kappa, rho, g; per-slice a if present)
inline std::vector<Deriv> DerivMain(const Spec& S) {
    std::vector<Deriv> v;
    const int ia = S.Idx("a_eta"), ik = S.Idx("kappa"), ir = S.Idx("rho_d"), ig = S.Idx("g");
    if (ia < 0) return v;
    auto Y = [=](int t, int w, bool kap, bool rho) {
        return [=](const Vec& th, const Templ& T) { return kYf * th[ia] * (kap ? th[ik] : 1.) * (rho ? th[ir] : 1.) * T.W[t][w]; };
    };
    v.push_back({"a_eta", "eta scale relative to SMASH x PDG BR", [=](const Vec& th, const Templ&) { return th[ia]; }});
    v.push_back({"kappa", "conversion scale of the eta classes C3, C4", [=](const Vec& th, const Templ&) { return th[ik]; }});
    v.push_back({"rho_d", "direct-decay scale C5 relative to a_eta", [=](const Vec& th, const Templ&) { return th[ir]; }});
    v.push_back({"g", "G_eta factor", [=](const Vec& th, const Templ&) { return th[ig]; }});
    v.push_back({"kappa_B", "conversion scale of the non-eta background G_C", [=](const Vec& th, const Templ&) { return th[S.iKB]; }});
    v.push_back({"Y2_B", "C2 (pi0 Dalitz) in B", Y(kS2, kYB, false, false)});
    v.push_back({"Y3_B", "C3 (conversion) in B", Y(kS3, kYB, true, false)});
    v.push_back({"Y_bump", "Y2 + Y3 in B", [=](const Vec& th, const Templ& T) { return kYf * th[ia] * (T.W[kS2][kYB] + th[ik] * T.W[kS3][kYB]); }});
    v.push_back({"Y4_P", "C4 (eta -> pi+pi-gamma conv.) in P", Y(kS4, kYP, true, false)});
    v.push_back({"Y5_P", "C5 (direct) in P", Y(kS5, kYP, false, true)});
    v.push_back({"Y_peak", "Y4 + Y5 in P", [=](const Vec& th, const Templ& T) { return kYf * th[ia] * (th[ik] * T.W[kS4][kYP] + th[ir] * T.W[kS5][kYP]); }});
    v.push_back({"Y2_tot", "C2, all M4", Y(kS2, kYTot, false, false)});
    v.push_back({"Y3_tot", "C3, all M4", Y(kS3, kYTot, true, false)});
    v.push_back({"Y4_tot", "C4, all M4", Y(kS4, kYTot, true, false)});
    v.push_back({"Y5_tot", "C5, all M4", Y(kS5, kYTot, false, true)});
    v.push_back({"Yeta_B", "G_eta in B", [=](const Vec& th, const Templ& T) { return kYf * th[ia] * th[ig] * T.W[kSeta][kYB]; }});
    v.push_back({"Yeta_P", "G_eta in P", [=](const Vec& th, const Templ& T) { return kYf * th[ia] * th[ig] * T.W[kSeta][kYP]; }});
    v.push_back({"DR", "(Y_peak/Y_bump)_data / (Y_peak/Y_bump)_sim",
                 [=](const Vec& th, const Templ& T) {
                     const double d = (th[ik] * T.W[kS4][kYP] + th[ir] * T.W[kS5][kYP]) / (T.W[kS2][kYB] + th[ik] * T.W[kS3][kYB]);
                     const double s = (T.W[kS4][kYP] + T.W[kS5][kYP]) / (T.W[kS2][kYB] + T.W[kS3][kYB]);
                     return d / s; }});
    v.push_back({"BR_A", "2.68e-4 rho_d", [=](const Vec& th, const Templ&) { return kBRdir * th[ir]; }});
    v.push_back({"sigvis_bump", "Y_bump / L [nb]", [=](const Vec& th, const Templ& T) { return kYf * th[ia] * (T.W[kS2][kYB] + th[ik] * T.W[kS3][kYB]) / kL; }});
    v.push_back({"sigvis_peak", "Y_peak / L [nb]", [=](const Vec& th, const Templ& T) { return kYf * th[ia] * (th[ik] * T.W[kS4][kYP] + th[ir] * T.W[kS5][kYP]) / kL; }});
    return v;
}
inline std::vector<Deriv> DerivFree(const Spec& S) {
    std::vector<Deriv> v;
    const int i2 = S.Idx("alpha_2"), i3 = S.Idx("alpha_3"), i4 = S.Idx("alpha_4"), i5 = S.Idx("alpha_5");
    for (int q : {i2, i3, i4, i5}) { const std::string nm = S.names[q]; v.push_back({nm, "free class scale", [=](const Vec& th, const Templ&) { return th[q]; }}); }
    v.push_back({"Y_bump", "alpha_2 C2 + alpha_3 C3 in B", [=](const Vec& th, const Templ& T) { return kYf * (th[i2] * T.W[kS2][kYB] + th[i3] * T.W[kS3][kYB]); }});
    v.push_back({"Y_peak", "alpha_4 C4 + alpha_5 C5 in P", [=](const Vec& th, const Templ& T) { return kYf * (th[i4] * T.W[kS4][kYP] + th[i5] * T.W[kS5][kYP]); }});
    v.push_back({"Y5_P", "alpha_5 C5 in P", [=](const Vec& th, const Templ& T) { return kYf * th[i5] * T.W[kS5][kYP]; }});
    v.push_back({"Y4_P", "alpha_4 C4 in P", [=](const Vec& th, const Templ& T) { return kYf * th[i4] * T.W[kS4][kYP]; }});
    v.push_back({"R_34", "alpha_4/alpha_3", [=](const Vec& th, const Templ&) { return th[i4] / th[i3]; }});
    v.push_back({"r5", "alpha_5/alpha_2", [=](const Vec& th, const Templ&) { return th[i5] / th[i2]; }});
    v.push_back({"BR_C", "2.68e-4 alpha_5 (alpha_3/alpha_2)/alpha_4", [=](const Vec& th, const Templ&) { return kBRdir * th[i5] * th[i3] / (th[i2] * th[i4]); }});
    v.push_back({"kappa_B", "conversion scale of G_C", [=](const Vec& th, const Templ&) { return th[S.iKB]; }});
    return v;
}
inline std::vector<Deriv> DerivGold(const Spec& S) {
    std::vector<Deriv> v;
    const int i2 = S.Idx("alpha_2"), i5 = S.Idx("alpha_5"), ig = S.Idx("g");
    v.push_back({"alpha_2", "C2 scale (M(ee) >= 20 MeV)", [=](const Vec& th, const Templ&) { return th[i2]; }});
    v.push_back({"alpha_5", "C5 scale (M(ee) >= 20 MeV)", [=](const Vec& th, const Templ&) { return th[i5]; }});
    v.push_back({"g", "G_eta factor", [=](const Vec& th, const Templ&) { return th[ig]; }});
    v.push_back({"BR_G", "2.68e-4 alpha_5/alpha_2", [=](const Vec& th, const Templ&) { return kBRdir * th[i5] / th[i2]; }});
    return v;
}

// ===========================================================================
// 6. GLS chi2 with the MC covariance of the prediction (100 replicas)
// ===========================================================================
struct GLSRes { double chi2 = 0, chi2diag = 0, chi2min = 0; int ndf = 0; Vec th; bool ok = false; };
inline bool InvSym(std::vector<double>& M, int n) {
    TMatrixDSym m(n);
    for (int a = 0; a < n; ++a) for (int b = 0; b < n; ++b) m(a, b) = M[(size_t)a * n + b];
    TDecompChol ch(m);
    if (!ch.Decompose()) return false;
    TMatrixDSym inv(n);
    ch.Invert(inv);
    for (int a = 0; a < n; ++a) for (int b = 0; b < n; ++b) M[(size_t)a * n + b] = inv(a, b);
    return true;
}
// prediction MC covariance at theta: mu_r(theta) over replica templates
inline std::vector<double> PredCovMC(const Bins& B0, const Spec& S, const Vec& th, const std::vector<Templ>& R) {
    const int n = B0.n;
    std::vector<Vec> mus;
    Bins B = B0;
    for (const auto& T : R) { SwapTempl(B, T); Vec mu; Eval(B, S, th, mu, nullptr, nullptr, nullptr); mus.push_back(mu); }
    Vec m(n, 0.);
    for (auto& v : mus) for (int i = 0; i < n; ++i) m[i] += v[i];
    for (int i = 0; i < n; ++i) m[i] /= mus.size();
    std::vector<double> C((size_t)n * n, 0.);
    for (auto& v : mus) for (int a = 0; a < n; ++a) { const double da = v[a] - m[a]; for (int b = 0; b < n; ++b) C[(size_t)a * n + b] += da * (v[b] - m[b]); }
    for (auto& c : C) c /= (mus.size() - 1);
    return C;
}
// data covariance of the counts (event groups scaled to the model + CB), full or diagonal
inline std::vector<double> DataCov(const Bins& B, const DataSet& D, const Spec& S, const Vec& th, const Vec& mu, bool full) {
    const int n = B.n, nc = B.lay == 0 ? 700 : 350;
    const std::vector<double>& eg = B.lay == 0 ? D.eg10 : D.eg5;
    const std::vector<double>& cc = B.lay == 0 ? D.covCB10 : D.covCB5;
    std::vector<double> C((size_t)n * n, 0.);
    Vec sf(n), kk(n);
    for (int i = 0; i < n; ++i) { sf[i] = B.y[i] > 0 ? std::sqrt(std::max(mu[i], 1e-9) / B.y[i]) : 0.; kk[i] = S.iK[B.e[i]] >= 0 ? th[S.iK[B.e[i]]] : S.kFix[B.e[i]]; }
    for (int a = 0; a < n; ++a)
        for (int b = 0; b < n; ++b) {
            if (!full && a != b) continue;
            double s = eg[(size_t)B.cell[a] * nc + B.cell[b]];
            if (B.y[a] > 0 && B.y[b] > 0) s *= sf[a] * sf[b]; else s = (a == b) ? std::max(mu[a], 1.) : 0.;
            s += kk[a] * kk[b] * cc[(size_t)B.cell[a] * nc + B.cell[b]];
            C[(size_t)a * n + b] = s;
        }
    return C;
}
// GLS: chi2 at theta_hat (full data covariance + MC), with diagonal data covariance, and minimised (Gauss-Newton, 3 updates of MC)
inline GLSRes GLS(const Bins& B, const Spec& S, const DataSet& D, const FitRes& F, const std::vector<Templ>& R) {
    GLSRes g;
    const int n = B.n, np = S.np();
    auto chi2At = [&](const Vec& th, const std::vector<double>& Winv, Vec* r_out) {
        Vec mu; Eval(B, S, th, mu, nullptr, nullptr, nullptr);
        Vec r(n); for (int i = 0; i < n; ++i) r[i] = B.y[i] - mu[i];
        double c = 0;
        for (int a = 0; a < n; ++a) { double s = 0; for (int b = 0; b < n; ++b) s += Winv[(size_t)a * n + b] * r[b]; c += r[a] * s; }
        for (const auto& q : S.cons) { const double z = (th[q.k] - q.val) / q.sig; c += z * z; }
        if (r_out) *r_out = r;
        return c;
    };
    std::vector<double> Cmc = PredCovMC(B, S, F.th, R);
    {
        std::vector<double> W = DataCov(B, D, S, F.th, F.mu, true);
        for (size_t k = 0; k < W.size(); ++k) W[k] += Cmc[k];
        if (!InvSym(W, n)) return g;
        g.chi2 = chi2At(F.th, W, nullptr);
        std::vector<double> Wd = DataCov(B, D, S, F.th, F.mu, false);
        for (size_t k = 0; k < Wd.size(); ++k) Wd[k] += Cmc[k];
        if (!InvSym(Wd, n)) return g;
        g.chi2diag = chi2At(F.th, Wd, nullptr);
    }
    // minimisation
    std::vector<int> fr; for (int k = 0; k < np; ++k) if (!S.fixed[k]) fr.push_back(k);
    const int nf = (int)fr.size();
    Vec th = F.th;
    for (int outer = 0; outer < 3; ++outer) {
        Vec mu0; Eval(B, S, th, mu0, nullptr, nullptr, nullptr);
        std::vector<double> W = DataCov(B, D, S, th, mu0, true);
        std::vector<double> Cm = PredCovMC(B, S, th, R);
        for (size_t k = 0; k < W.size(); ++k) W[k] += Cm[k];
        if (!InvSym(W, n)) return g;
        for (int it = 0; it < 30; ++it) {
            Vec mu; Mat J; Eval(B, S, th, mu, &J, nullptr, nullptr);
            Vec r(n); for (int i = 0; i < n; ++i) r[i] = B.y[i] - mu[i];
            Vec Wr(n, 0.); Mat WJ = Zero(n, nf);
            for (int a = 0; a < n; ++a) for (int b = 0; b < n; ++b) { const double w = W[(size_t)a * n + b]; if (w == 0) continue; Wr[a] += w * r[b]; for (int k = 0; k < nf; ++k) WJ[a][k] += w * J[b][fr[k]]; }
            Mat A = Zero(nf, nf); Vec bb(nf, 0.);
            for (int a = 0; a < n; ++a) for (int k = 0; k < nf; ++k) { bb[k] += J[a][fr[k]] * Wr[a]; for (int l = 0; l < nf; ++l) A[k][l] += J[a][fr[k]] * WJ[a][l]; }
            for (const auto& q : S.cons) for (int k = 0; k < nf; ++k) if (fr[k] == q.k) { A[k][k] += 1. / (q.sig * q.sig); bb[k] -= (th[q.k] - q.val) / (q.sig * q.sig); }
            Mat Ai = A;
            Vec sc(nf, 1.); for (int k = 0; k < nf; ++k) sc[k] = Ai[k][k] > 0 ? 1. / std::sqrt(Ai[k][k]) : 1.;
            for (int k = 0; k < nf; ++k) for (int l = 0; l < nf; ++l) Ai[k][l] *= sc[k] * sc[l];
            if (!bla::CholInv(Ai)) return g;
            for (int k = 0; k < nf; ++k) for (int l = 0; l < nf; ++l) Ai[k][l] *= sc[k] * sc[l];
            double dmax = 0;
            for (int k = 0; k < nf; ++k) { double d = 0; for (int l = 0; l < nf; ++l) d += Ai[k][l] * bb[l]; th[fr[k]] += d; dmax = std::max(dmax, std::fabs(d) / std::sqrt(std::max(Ai[k][k], 1e-300))); }
            if (dmax < 1e-6) break;
        }
        if (outer == 2) g.chi2min = chi2At(th, W, nullptr);
    }
    g.th = th;
    g.ndf = F.ndf;
    g.ok = true;
    return g;
}

// longest run of same-sign pulls in a set of bins (ordered); also the longest run with |pull| > 1
inline void Runs(const Vec& p, int& longest, int& longest1) {
    longest = longest1 = 0;
    int cur = 0, cur1 = 0, sg = 0, sg1 = 0;
    for (double v : p) {
        const int s = v > 0 ? 1 : (v < 0 ? -1 : 0);
        if (s != 0 && s == sg) ++cur; else { cur = s != 0 ? 1 : 0; sg = s; }
        longest = std::max(longest, cur);
        const int s1 = std::fabs(v) > 1 ? s : 0;
        if (s1 != 0 && s1 == sg1) ++cur1; else { cur1 = s1 != 0 ? 1 : 0; sg1 = s1; }
        longest1 = std::max(longest1, cur1);
    }
}
// probability that the longest same-sign run among m random signs is >= L (exact, dynamic programming)
inline double ProbRunGE(int m, int L) {
    // states: current run length 1..L-1 (no run >= L yet); start with run 1 after the first bin
    std::vector<double> p(L, 0.);
    if (L <= 1) return 1.;
    p[1] = 1.;
    for (int k = 2; k <= m; ++k) {
        std::vector<double> q(L, 0.);
        for (int r = 1; r < L; ++r) { q[1] += 0.5 * p[r]; if (r + 1 < L) q[r + 1] += 0.5 * p[r]; }
        p = q;
    }
    double s = 0; for (int r = 1; r < L; ++r) s += p[r];
    return 1. - s;
}

}  // namespace wp5

#endif
