// wp3_stageA_cnt_plain2.C -- COPY of wp3_stageA_cnt_plain.C (2026-11-16, writer of report_all_short.tex). Only two label strings of the
// map change: "control (dropped)" -> "control (not used)" and the legend "(dashed: 0.565, first choice)" -> "(dashed: 0.565, free-scale fits)"
// (0.565-0.70 is the region of the fits with a free scale per SMASH channel). Writes docs/latex/figs_all/wp3_map_cnt_plain2.{pdf,png};
// everything else goes to the scratch folder kCntOut. Every number unchanged.
// Run: root -l -b -q -e 'gSystem->SetBuildDir("<scratch>/aclic_wp3p2", kTRUE)' 'plots/all_report/wp3_stageA_cnt_plain2.C+O'
// wp3_stageA_cnt_plain.C -- COPY of wp3_stageA_cnt.C (2026-11-16, short-report figure agent for report_all_short.tex). Runs the
// same stage A on the same inputs; writes to docs/latex only figs_all/wp3_map_cnt_plain.{pdf,png}; its log, stageA.txt/.root,
// tables, the EPS files and the other figures go to the scratch folder kCntOut. The figures are saved with SaveNoRm (EPS ->
// ps2pdf -dEPSCrop -> PDF, then PNG, as bla::SaveFig, but the EPS is kept: nothing is removed). Label strings of the map in plain
// words (instruction I2): S0, isBest, spread rule -> "both pions in vertex fit"; CR-U -> normalisation region (0.575-0.700; dashed
// 0.565, the first choice); CR-K -> "control (dropped)"; CR-L -> "test"; B -> "bump"; P -> "peak" (written upright along the narrow
// window, as "test"); edge E(M(pipi), q) -> "edge, M_ee = q" (the kinematic edge of the bump). Every number unchanged.
// Run: root -l -b -q -e 'gSystem->SetBuildDir("<scratch>/aclic/wp3_plain", kTRUE)' 'plots/all_report/wp3_stageA_cnt_plain.C+O'
// wp3_stageA_cnt.C -- COPY of wp3_stageA.C (figure agent B, 2026-11-15): "counts" instead of "rows" in the y title of
// Fig. wp3_smoothing and the z title of Fig. wp3_map (-> docs/latex/figs_all/wp3_smoothing_cnt, wp3_map_cnt .pdf/.png). It reads
// the same inputs (hists_v1.root, wp3/wp3_rows.root, wp3_rows_summary.txt); its log, stageA.txt/.root, tables and the other
// figures go to /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figB/out/wp3_stageA/ (the WP3 files are not touched). Nothing else changed.
// Run: root -l -b -q -e 'gSystem->SetBuildDir("<scratch>/aclic_wp3s_cnt", kTRUE)' 'plots/all_report/wp3_stageA_cnt.C+O'
// ============================================================================
// wp3_stageA.C  --  WP3 control regions, stage A (SPEC report_all, D9-D11, WP3)
// ============================================================================
// Written by WP3 (prefix wp3), 2026-11-14. Reads hists_v1.root (WP0) and wp3_rows.root (wp3_rows.C).
//
// What it does, in plain words
//   1. Inputs. The 10 F10 histograms (5 M(pi+pi-) slices x 2 M(ee) classes, M4 0.40-0.70) of the data
//      (N, CB = 2 sqrt(N++ N--) spread rule, CB replica variance, event-group covariance of N) and of the
//      simulation (SMASH, mode A, w x V3, at SMASH sigma_prod x PDG BR) in expected RAW data rows
//      (nb x kNbToRaw = 6079.2 x 63/64).
//   2. Smoothing (D13). The two non-eta background templates G_D = G_pi0D + G_vec + G_oth and G_C = G_conv
//      are replaced, per histogram, by a smooth function exp(p0 + p1 t + p2 t^2 + p3 t^3), t = (M4 - 0.55)/0.15,
//      integrated over each bin and fitted to the simulated histogram by a Poisson-type likelihood (this keeps
//      the integral). The fit starts at the first bin with simulated content (kinematic threshold), and is
//      zero below it. The eta templates (C2-C5, G_eta) stay raw; their bin variance enters the fit (BB-lite).
//   3. Stage-A model in every CR bin i of histogram (j = slice, e = M(ee) class):
//        mu = k_e CB + c_j (1 + b x + q x^2) [G_D + kappa_B G_C] + a_eta [C2 + kappa (C3 + C4) + rho_d C5 + g G_eta]
//      with x = (M4 - 0.55)/0.1. Free: c_1..c_5, b, kappa_B, k_E1; constrained k_E2 = 1 +- 0.05; the eta part
//      FIXED at a_eta(0) = 0.4, kappa = rho_d = g = 1. Poisson likelihood of the unsubtracted N (iteratively
//      re-weighted least squares with variance mu + k^2 var(CB) + var(eta templates)); errors: sandwich with
//      the event-group covariance of N; MC error separately from 100 Poisson replicas of the simulation
//      (each one smoothed again and refitted).
//   4. The D11 rule: a quadratic term q, then per-slice tilts b_j, are added only if they lower chi2 by more
//      than 4. Nested models M0 (scale only) .. M4 are all fitted and listed.
//   5. Tests A1-A5, the fall-backs of Sect. 6, variants (pre, a_eta = 0/1/2, CB per-cell, k_E1 = 1, CR edges).
//   6. Outputs: tables C1-C3 (docs/latex/tables_all/wp3_*), figures C1-C3 (docs/latex/figs_all/wp3_*),
//      the frozen parameters /hdd3/przygoda/HADES/FAT_out/all_report/wp3/stageA.txt (+ stageA.root).
//
// Run (from /home/przygoda/HADES/FAT/FAT; about 2 GB, a few minutes):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_wp3", kTRUE)' \
//        'plots/all_report/wp3_stageA.C+O'
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <map>
#include <string>
#include <vector>

#include "TBox.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "THStack.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TMatrixDSym.h"
#include "TPad.h"
#include "TStyle.h"
#include "TSystem.h"

#include "hists_io.h"
#include "s5_blind_ana_fit.h"

namespace wp3 {
using namespace alldefs;
typedef std::vector<double> Vec;
typedef std::vector<Vec> Mat;

const char* const kOutDir = "/hdd3/przygoda/HADES/FAT_out/all_report/wp3";
const char* const kFigDir = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_all";
const char* const kTabDir = "/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ6/allshort_final/wp3out";   // tables to my scratch
const char* const kCntOut = "/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ6/allshort_final/wp3out";
// _plain: as bla::SaveFig (EPS -> ps2pdf -dEPSCrop -> PDF, then PNG) but the EPS is written to kCntOut and kept
void SaveNoRm(TCanvas* c, const std::string& base) {
    const std::string eps = std::string(kCntOut) + "/" + gSystem->BaseName(base.c_str()) + ".eps";
    c->SaveAs(eps.c_str());
    gSystem->Exec(Form("ps2pdf -dEPSCrop %s %s.pdf", eps.c_str(), base.c_str()));
    c->SaveAs((base + ".png").c_str());
}   // _cnt copy: log, stageA.txt/.root to scratch
double kAeta0 = 0.4;                // a_eta(0) of stage A (argument of wp3_stageA; 0.4 nominal, rerun at the WP5 value)
const double kKE2sig = 0.05;        // k_E2 = 1 +- 0.05
const double kStepS = 0.0054;       // edge smearing (= peak width)
const int kNH = 10, kNB = kNF10;

FILE* gOut = nullptr;
double gChi2MC = 0;     // chi2 of the adopted fit including the MC template covariance
bool gSimLimited = false;
#define P(...) do { printf(__VA_ARGS__); if (gOut) fprintf(gOut, __VA_ARGS__); } while (0)

// ------------------------------------------------------------------------------------------------
// simulation components used here (index into the local arrays)
// ------------------------------------------------------------------------------------------------
enum LC { lGD = 0, lGC, lGeta, lC2, lC3, lC4, lC5, lT12, lGDom, lGCom, kNLC };
const int kLCcomp[kNLC] = {kGD, kGconv, kGeta, kC2, kC3, kC4, kC5, kT12, -1, kGconvom};   // lGDom = Gpi0D_om + Gvec_om + Goth_om
const char* const kLCname[kNLC] = {"GD", "GC", "Geta", "C2", "C3", "C4", "C5", "T12", "GDom", "GCom"};

// ------------------------------------------------------------------------------------------------
// one histogram of the fit (normally an F10 histogram; variants may be built from G3)
// ------------------------------------------------------------------------------------------------
struct HIn {
    std::string name;
    int j = 0, e = 0;                        // slice (c_j, tilt b_j), M(ee) class (k_e)
    bool cruOK = true, crkOK = false;        // may its bins belong to CR-U / CR-K
    Vec lo, hi, x;                           // M4 bin edges and tilt variable
    Vec N, CB, vCB, Npp, Nmm, CBgm;          // data (raw rows); CBgm = per-cell geometric-mean CB variant
    std::vector<int> eg;                     // index into the event-group covariance (-1: none)
    Vec s[kNLC], v[kNLC];                    // simulation raw rows and variance
    Vec D, C, Dom, Com;                      // smoothed background (raw rows)
    int domD = 0, domC = 0, domDom = 0, domCom = 0;
    double pD[4] = {0, 0, 0, 0}, pC[4] = {0, 0, 0, 0};
    double chiD = 0, chiC = 0; int ndfD = 0, ndfC = 0;
};

// ------------------------------------------------------------------------------------------------
// smoothing: exp(pol3) integrated over the bins, Poisson-type likelihood (Newton with step halving)
// ------------------------------------------------------------------------------------------------
const double kGLx[4] = {-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526};
const double kGLw[4] = {0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538};
inline double Tsm(double m) { return (m - 0.55) / 0.15; }
struct SmRes { bool ok = false; int dom = 0; double p[4] = {0, 0, 0, 0}; Vec f; double chi2 = 0; int ndf = 0; };
// bin integral of exp(P(t)) and of t^k exp(P(t)), k = 0..6
inline void BinMoments(double lo, double hi, const double* p, double* mom) {
    for (int k = 0; k < 7; ++k) mom[k] = 0;
    const double h = 0.5 * (hi - lo), c = 0.5 * (hi + lo);
    for (int g = 0; g < 4; ++g) {
        const double m = c + h * kGLx[g], t = Tsm(m);
        const double e = std::exp(std::min(700., p[0] + t * (p[1] + t * (p[2] + t * p[3])))) * h * kGLw[g];
        double tk = 1;
        for (int k = 0; k < 7; ++k) { mom[k] += tk * e; tk *= t; }
    }
}
SmRes Smooth(const Vec& y, const Vec& e2, const Vec& lo, const Vec& hi, int dom = -1, const double* pstart = nullptr) {
    SmRes r;
    const int n = (int)y.size();
    if (dom < 0) { dom = n; for (int i = 0; i < n; ++i) if (y[i] > 0) { dom = i; break; } }
    r.dom = dom;
    r.f.assign(n, 0.);
    double tot = 0; for (int i = dom; i < n; ++i) tot += y[i];
    if (dom >= n - 4 || !(tot > 0)) return r;
    double p[4] = {0, 0, 0, 0};
    if (pstart) for (int k = 0; k < 4; ++k) p[k] = pstart[k];
    else {
        // start: weighted least squares of ln(density) on bins with content
        Mat A = bla::Zero(4, 4); Vec b(4, 0.);
        for (int i = dom; i < n; ++i) {
            if (!(y[i] > 0)) continue;
            const double t = Tsm(0.5 * (lo[i] + hi[i])), w = y[i], z = std::log(y[i] / (hi[i] - lo[i]));
            double tk[4] = {1, t, t * t, t * t * t};
            for (int k = 0; k < 4; ++k) { b[k] += w * tk[k] * z; for (int l = 0; l < 4; ++l) A[k][l] += w * tk[k] * tk[l]; }
        }
        for (int k = 0; k < 4; ++k) A[k][k] += 1e-9;
        if (bla::CholInv(A)) for (int k = 0; k < 4; ++k) { p[k] = 0; for (int l = 0; l < 4; ++l) p[k] += A[k][l] * b[l]; }
        else { p[0] = std::log(tot / (hi[n - 1] - lo[dom])); }
    }
    auto F = [&](const double* q) {
        double s = 0, mom[7];
        for (int i = dom; i < n; ++i) {
            BinMoments(lo[i], hi[i], q, mom);
            const double f = std::max(mom[0], 1e-300);
            s += f - (y[i] > 0 ? y[i] * std::log(f) : 0.);
        }
        return s;
    };
    double Fc = F(p);
    for (int it = 0; it < 100; ++it) {
        Mat H = bla::Zero(4, 4); Vec g(4, 0.);
        double mom[7];
        for (int i = dom; i < n; ++i) {
            BinMoments(lo[i], hi[i], p, mom);
            const double f = std::max(mom[0], 1e-300), a = 1.0 - y[i] / f, c2 = y[i] / (f * f);
            for (int k = 0; k < 4; ++k) {
                g[k] += a * mom[k];
                for (int l = 0; l < 4; ++l) H[k][l] += a * mom[k + l] + c2 * mom[k] * mom[l];
            }
        }
        // make H positive definite if needed (Levenberg)
        Mat Hi = H;
        double lam = 0;
        bool inv = false;
        for (int tr = 0; tr < 20 && !inv; ++tr) {
            Hi = H;
            for (int k = 0; k < 4; ++k) Hi[k][k] += lam * std::max(1e-12, std::fabs(H[k][k]));
            inv = bla::CholInv(Hi);
            lam = lam == 0 ? 1e-6 : lam * 10;
        }
        if (!inv) break;
        double d[4];
        for (int k = 0; k < 4; ++k) { d[k] = 0; for (int l = 0; l < 4; ++l) d[k] -= Hi[k][l] * g[l]; }
        double st = 1.0, Fn = 0, pn[4];
        for (int h = 0; h < 40; ++h) {
            for (int k = 0; k < 4; ++k) pn[k] = p[k] + st * d[k];
            Fn = F(pn);
            if (Fn <= Fc + 1e-12 * std::fabs(Fc)) break;
            st *= 0.5;
        }
        const double dF = Fc - Fn;
        if (Fn <= Fc) { for (int k = 0; k < 4; ++k) p[k] = pn[k]; Fc = Fn; }
        if (dF >= 0 && dF < 1e-10 * std::max(1.0, std::fabs(Fc))) break;
        if (st < 1e-10) break;
    }
    for (int k = 0; k < 4; ++k) r.p[k] = p[k];
    double sy = 0, se = 0;
    for (int i = dom; i < n; ++i) { double mom[7]; BinMoments(lo[i], hi[i], p, mom); r.f[i] = mom[0]; sy += y[i]; se += e2[i]; }
    const double wbar = sy > 0 ? se / sy : 1.0;
    r.chi2 = 0; r.ndf = n - dom - 4;
    for (int i = dom; i < n; ++i) if (r.f[i] > 0) r.chi2 += (y[i] - r.f[i]) * (y[i] - r.f[i]) / (wbar * r.f[i]);
    r.ok = true;
    return r;
}

// ------------------------------------------------------------------------------------------------
// generic fitter: Poisson ML by iteratively re-weighted least squares (Gauss-Newton), BB-lite variance,
// Gaussian constraints, fixed parameters, sandwich covariance with the event-group covariance of N.
// ------------------------------------------------------------------------------------------------
struct GProb {
    int np = 0;
    Vec th0;
    std::vector<bool> fixed;
    std::vector<int> consK; Vec consV, consS;
    Vec y;                                  // data per fit bin
    std::vector<int> eg;                    // event-group index per fit bin (-1: none)
    const TMatrixDSym* egcov = nullptr;
    // th, bin -> mu, gradient (np), extra variance
    std::function<void(const Vec&, int, double&, Vec&, double&)> eval;
    bool sandwich = true;
    double vfloor = 0.25;
};
struct GRes {
    bool ok = false;
    Vec th; Mat C, Cs; Vec mu, V, vex;
    double chi2P = 0, chi2BC = 0, chi2E = 0; int ndf = 0, n = 0, nfree = 0, it = 0;
    double Err(int k) const { return (k < (int)Cs.size() && Cs[k][k] > 0) ? std::sqrt(Cs[k][k]) : 0.; }
    double ErrM(int k) const { return (k < (int)C.size() && C[k][k] > 0) ? std::sqrt(C[k][k]) : 0.; }
};
GRes GFit(const GProb& P_) {
    GRes r;
    const int n = (int)P_.y.size(), np = P_.np;
    r.n = n;
    std::vector<int> fr;
    for (int k = 0; k < np; ++k) if (P_.fixed.empty() || !P_.fixed[k]) fr.push_back(k);
    const int nf = (int)fr.size();
    r.nfree = nf;
    Vec th = P_.th0;
    Vec mu(n), V(n), vex(n), g(np);
    Mat J(n, Vec(np, 0.));
    auto evalAll = [&](const Vec& t) {
        for (int i = 0; i < n; ++i) { P_.eval(t, i, mu[i], g, vex[i]); J[i] = g; V[i] = std::max(mu[i], P_.vfloor) + vex[i]; }
    };
    Mat Hi;
    bool conv = false;
    for (int it = 0; it < 200 && !conv; ++it) {
        evalAll(th);
        Mat H = bla::Zero(nf, nf); Vec b(nf, 0.);
        for (int i = 0; i < n; ++i) {
            const double w = 1.0 / V[i], res = P_.y[i] - mu[i];
            for (int a = 0; a < nf; ++a) {
                const double ja = J[i][fr[a]] * w;
                if (ja == 0) continue;
                b[a] += ja * res;
                for (int c = 0; c <= a; ++c) H[a][c] += ja * J[i][fr[c]];
            }
        }
        for (size_t q = 0; q < P_.consK.size(); ++q) {
            for (int a = 0; a < nf; ++a) if (fr[a] == P_.consK[q]) {
                const double s2 = P_.consS[q] * P_.consS[q];
                H[a][a] += 1.0 / s2; b[a] += (P_.consV[q] - th[fr[a]]) / s2;
            }
        }
        for (int a = 0; a < nf; ++a) for (int c = a + 1; c < nf; ++c) H[a][c] = H[c][a];
        Vec sc(nf, 1.);
        for (int a = 0; a < nf; ++a) sc[a] = H[a][a] > 0 ? 1.0 / std::sqrt(H[a][a]) : 1.0;
        Mat Hs = H;
        for (int a = 0; a < nf; ++a) for (int c = 0; c < nf; ++c) Hs[a][c] *= sc[a] * sc[c];
        if (!bla::CholInv(Hs)) { r.ok = false; return r; }
        for (int a = 0; a < nf; ++a) for (int c = 0; c < nf; ++c) Hs[a][c] *= sc[a] * sc[c];
        Hi = Hs;
        Vec d(nf, 0.);
        for (int a = 0; a < nf; ++a) for (int c = 0; c < nf; ++c) d[a] += Hi[a][c] * b[c];
        // damping: keep the model positive in the fit bins
        double st = 1.0;
        for (int h = 0; h < 30; ++h) {
            Vec tn = th;
            for (int a = 0; a < nf; ++a) tn[fr[a]] += st * d[a];
            bool pos = true;
            for (int i = 0; i < n && pos; ++i) { double m, ve; Vec gg(np); P_.eval(tn, i, m, gg, ve); if (m < 0 && P_.y[i] > 0) pos = false; }
            if (pos) break;
            st *= 0.5;
        }
        double dmax = 0;
        for (int a = 0; a < nf; ++a) { th[fr[a]] += st * d[a]; dmax = std::max(dmax, std::fabs(st * d[a]) / std::sqrt(std::max(Hi[a][a], 1e-300))); }
        r.it = it + 1;
        if (it > 1 && dmax < 1e-7) conv = true;
    }
    evalAll(th);
    // final covariance (full np x np, zero for fixed parameters)
    {
        Mat H = bla::Zero(nf, nf);
        for (int i = 0; i < n; ++i) {
            const double w = 1.0 / V[i];
            for (int a = 0; a < nf; ++a) for (int c = 0; c < nf; ++c) H[a][c] += J[i][fr[a]] * w * J[i][fr[c]];
        }
        Mat Hc = H;
        for (size_t q = 0; q < P_.consK.size(); ++q)
            for (int a = 0; a < nf; ++a) if (fr[a] == P_.consK[q]) Hc[a][a] += 1.0 / (P_.consS[q] * P_.consS[q]);
        Mat Hinv = Hc;
        if (!bla::CholInv(Hinv)) { r.ok = false; return r; }
        r.C = bla::Zero(np, np);
        for (int a = 0; a < nf; ++a) for (int c = 0; c < nf; ++c) r.C[fr[a]][fr[c]] = Hinv[a][c];
        r.Cs = r.C;
        if (P_.sandwich) {
            // meat = sum_ab u_a Sigma_ab u_b^T + constraint precision; u_a = J_a / V_a
            Mat M = bla::Zero(nf, nf);
            Vec sfac(n, 1.0);
            for (int i = 0; i < n; ++i) sfac[i] = P_.y[i] > 0 ? std::sqrt(std::max(mu[i], 1e-9) / P_.y[i]) : 0.;
            for (int a = 0; a < n; ++a) {
                double sa = 0;   // diagonal: model-scaled event-group variance (+ template variance)
                if (P_.egcov && P_.eg[a] >= 0 && P_.y[a] > 0) sa = (*P_.egcov)(P_.eg[a], P_.eg[a]) * sfac[a] * sfac[a];
                else sa = std::max(mu[a], 0.);
                sa += vex[a];
                for (int k = 0; k < nf; ++k) for (int l = 0; l < nf; ++l) M[k][l] += J[a][fr[k]] * J[a][fr[l]] * sa / (V[a] * V[a]);
                if (!(P_.egcov && P_.eg[a] >= 0 && P_.y[a] > 0)) continue;
                for (int b2 = 0; b2 < n; ++b2) {
                    if (b2 == a || P_.eg[b2] < 0 || !(P_.y[b2] > 0)) continue;
                    const double s = (*P_.egcov)(P_.eg[a], P_.eg[b2]);
                    if (s == 0) continue;
                    const double sab = s * sfac[a] * sfac[b2] / (V[a] * V[b2]);
                    for (int k = 0; k < nf; ++k) for (int l = 0; l < nf; ++l) M[k][l] += J[a][fr[k]] * J[b2][fr[l]] * sab;
                }
            }
            for (size_t q = 0; q < P_.consK.size(); ++q)
                for (int a = 0; a < nf; ++a) if (fr[a] == P_.consK[q]) M[a][a] += 1.0 / (P_.consS[q] * P_.consS[q]);
            Mat Cs = bla::Zero(nf, nf);
            for (int a = 0; a < nf; ++a) for (int c = 0; c < nf; ++c) {
                double s = 0;
                for (int k = 0; k < nf; ++k) for (int l = 0; l < nf; ++l) s += Hinv[a][k] * M[k][l] * Hinv[l][c];
                Cs[a][c] = s;
            }
            r.Cs = bla::Zero(np, np);
            for (int a = 0; a < nf; ++a) for (int c = 0; c < nf; ++c) r.Cs[fr[a]][fr[c]] = Cs[a][c];
        }
    }
    r.th = th; r.mu = mu; r.V = V; r.vex = vex;
    r.chi2P = 0; r.chi2BC = 0;
    for (int i = 0; i < n; ++i) {
        const double m = std::max(mu[i], 1e-9), y = P_.y[i];
        r.chi2P += (y - mu[i]) * (y - mu[i]) / V[i];
        r.chi2BC += 2.0 * (m - y + (y > 0 ? y * std::log(y / m) : 0.));
        double se = (P_.egcov && P_.eg[i] >= 0 && y > 0) ? (*P_.egcov)(P_.eg[i], P_.eg[i]) / y : 1.0;
        r.chi2E += (y - mu[i]) * (y - mu[i]) / (std::max(mu[i], P_.vfloor) * se + vex[i]);
    }
    for (size_t q = 0; q < P_.consK.size(); ++q) {
        const double z = (th[P_.consK[q]] - P_.consV[q]) / P_.consS[q];
        r.chi2P += z * z; r.chi2BC += z * z; r.chi2E += z * z;
    }
    r.ndf = n - nf + (int)P_.consK.size();
    r.ok = conv;
    return r;
}

// ------------------------------------------------------------------------------------------------
// stage-A model
// ------------------------------------------------------------------------------------------------
struct Model {
    bool tilt = true, quad = false, perSlice = false, omega = false;
    bool kE1free = true;
    bool useCBvar = true;
    bool cbGM = false;                       // CB variant: per-cell geometric mean (10 MeV cells of the histogram)
    int eOnly = -1;                          // -1: E1 and E2; 1: E2 only (fall-back "drop E1"; k_E1 then fixed to 1)
    double kBsig = 0;                        // > 0: kappa_B = 1 +- kBsig (Gaussian constraint); 0: free
    double aEta = kAeta0, kappa = 1, rhod = 1, g = 1;
    // parameter layout
    int ic = 0, ib = -1, iq = -1, ikB = -1, ikE1 = -1, ikE2 = -1, iom = -1, np = 0;
    std::vector<std::string> names;
    void Layout() {
        names.clear();
        ic = 0; for (int j = 0; j < 5; ++j) names.push_back(Form("c_%d", j + 1));
        ib = -1;
        if (tilt) { ib = (int)names.size(); if (perSlice) for (int j = 0; j < 5; ++j) names.push_back(Form("b_%d", j + 1)); else names.push_back("b"); }
        iq = -1; if (quad) { iq = (int)names.size(); names.push_back("q"); }
        ikB = (int)names.size(); names.push_back("kappa_B");
        ikE1 = (int)names.size(); names.push_back("k_E1");
        ikE2 = (int)names.size(); names.push_back("k_E2");
        iom = -1; if (omega) { iom = (int)names.size(); names.push_back("c_om"); }
        np = (int)names.size();
    }
    std::string Tag() const {
        std::string s = !tilt ? "M0 scale only" : (perSlice ? (quad ? "M4 per-slice tilt + quadratic" : "M3 per-slice tilt") : (quad ? "M2 tilt + quadratic" : "M1 tilt"));
        if (omega) s += " + omega split";
        if (eOnly == 1) s += ", E2 only";
        if (kBsig > 0) s += Form(", kappa_B = 1 +- %.2f", kBsig);
        return s;
    }
};
// eta part (fixed) and its BB variance
inline double EtaSum(const HIn& h, int i, const Model& m, double* var = nullptr, double injS = 0, double injG = 0) {
    const double a = m.aEta;
    const double s = a * (h.s[lC2][i] + m.kappa * (h.s[lC3][i] + h.s[lC4][i]) + m.rhod * h.s[lC5][i] + m.g * h.s[lGeta][i]);
    if (var) *var = a * a * (h.v[lC2][i] + m.kappa * m.kappa * (h.v[lC3][i] + h.v[lC4][i]) + m.rhod * m.rhod * h.v[lC5][i] + m.g * m.g * h.v[lGeta][i]);
    (void)injS; (void)injG;
    return s;
}
// background (renormalised) and its gradient
inline double Bkg(const HIn& h, int i, const Model& m, const Vec& th, Vec* grad) {
    const int j = h.j;
    const double x = h.x[i];
    const double b = m.tilt ? th[m.ib + (m.perSlice ? j : 0)] : 0.;
    const double q = m.quad ? th[m.iq] : 0.;
    const double T = 1 + b * x + q * x * x;
    const double kB = th[m.ikB];
    double Dn = h.D[i], Cn = h.C[i], Do = 0, Co = 0, com = 0;
    if (m.omega) { Do = h.Dom[i]; Co = h.Com[i]; Dn = h.D[i] - Do; Cn = h.C[i] - Co; com = th[m.iom]; }
    const double Bk = (Dn + kB * Cn) + com * (Do + kB * Co);
    const double c = th[m.ic + j];
    if (grad) {
        Vec& g = *grad;
        g[m.ic + j] += T * Bk;
        if (m.tilt) g[m.ib + (m.perSlice ? j : 0)] += c * x * Bk;
        if (m.quad) g[m.iq] += c * x * x * Bk;
        g[m.ikB] += c * T * (Cn + com * Co);
        if (m.omega) g[m.iom] += c * T * (Do + kB * Co);
    }
    return c * T * Bk;
}
inline double CBof(const HIn& h, int i, const Model& m) { return m.cbGM ? h.CBgm[i] : h.CB[i]; }

// control-region definition
struct CRdef {
    double uLo = 0.565, uHi = 0.700;
    bool useK = true;
    double kLo = 0.400, kHi = 0.5225;
    std::string tag = "CR-U 0.565-0.700 + CR-K (S5, 0.400-0.5225)";
};
inline bool InCR(const HIn& h, int i, const CRdef& c) {
    const double lo = h.lo[i], hi = h.hi[i];
    if (h.cruOK && lo >= c.uLo - 1e-9 && hi <= c.uHi + 1e-9) return true;
    if (c.useK && h.crkOK && lo >= c.kLo - 1e-9 && hi <= c.kHi + 1e-9) return true;
    return false;
}

struct StageA {
    Model m; CRdef cr; GRes r;
    std::vector<std::pair<int, int>> bins;    // (hist, bin) of the fit bins
    bool ok() const { return r.ok; }
};

// injection: pseudo-data = N + dS x (C2..C5) + dG x G_eta  (dS, dG in units of SMASH)
StageA FitStageA(const std::vector<HIn>& H, const Model& m0, const CRdef& cr, const TMatrixDSym* eg, double dS = 0, double dG = 0,
                 bool sandwich = true, const Vec* start = nullptr) {
    StageA A; A.m = m0; A.m.Layout(); A.cr = cr;
    const Model& m = A.m;
    GProb P;
    P.np = m.np;
    P.th0.assign(m.np, 0.);
    for (int j = 0; j < 5; ++j) P.th0[m.ic + j] = 1.3;
    P.th0[m.ikB] = 1.0; P.th0[m.ikE1] = 1.0; P.th0[m.ikE2] = 1.0;
    if (m.omega) P.th0[m.iom] = 1.0;
    if (start && (int)start->size() == m.np) P.th0 = *start;
    P.fixed.assign(m.np, false);
    if (!m.kE1free || m.eOnly == 1) { P.fixed[m.ikE1] = true; P.th0[m.ikE1] = 1.0; }
    P.consK.push_back(m.ikE2); P.consV.push_back(1.0); P.consS.push_back(kKE2sig);
    if (m.kBsig > 0) { P.consK.push_back(m.ikB); P.consV.push_back(1.0); P.consS.push_back(m.kBsig); }
    for (size_t h = 0; h < H.size(); ++h)
        for (int i = 0; i < (int)H[h].lo.size(); ++i) {
            if (m.eOnly >= 0 && H[h].e != m.eOnly) continue;
            if (!InCR(H[h], i, cr)) continue;
            const HIn& x = H[h];
            const bool any = x.N[i] > 0 || CBof(x, i, m) > 0 || x.D[i] > 0 || x.C[i] > 0 || x.s[lGeta][i] > 0 || x.s[lC5][i] > 0;
            if (!any) continue;
            A.bins.push_back(std::make_pair((int)h, i));
            P.y.push_back(x.N[i] + dS * (x.s[lC2][i] + x.s[lC3][i] + x.s[lC4][i] + x.s[lC5][i]) + dG * x.s[lGeta][i]);
            P.eg.push_back(x.eg[i]);
        }
    P.egcov = eg;
    P.sandwich = sandwich;
    const std::vector<std::pair<int, int>> bins = A.bins;
    P.eval = [&H, m, bins](const Vec& th, int k, double& mu, Vec& g, double& vex) {
        const HIn& h = H[bins[k].first];
        const int i = bins[k].second;
        std::fill(g.begin(), g.end(), 0.);
        double ve = 0;
        const double eta = EtaSum(h, i, m, &ve);
        const double kE = th[h.e == 0 ? m.ikE1 : m.ikE2];
        const double cb = CBof(h, i, m);
        mu = kE * cb + Bkg(h, i, m, th, &g) + eta;
        g[h.e == 0 ? m.ikE1 : m.ikE2] += cb;
        vex = ve + (m.useCBvar ? kE * kE * h.vCB[i] : 0.);
    };
    A.r = GFit(P);
    return A;
}

// ------------------------------------------------------------------------------------------------
// inputs
// ------------------------------------------------------------------------------------------------
// per-cell geometric-mean CB in 10 MeV cells of one histogram (spread uniformly over the F10 bins of the cell)
void MakeCBgm(HIn& h) {
    const int n = (int)h.lo.size();
    h.CBgm.assign(n, 0.);
    std::map<int, double> pp, mm, w;
    for (int i = 0; i < n; ++i) { const int c = (int)std::floor((h.lo[i] - 0.30) / 0.010 + 1e-6); pp[c] += h.Npp[i]; mm[c] += h.Nmm[i]; w[c] += h.hi[i] - h.lo[i]; }
    for (int i = 0; i < n; ++i) {
        const int c = (int)std::floor((h.lo[i] - 0.30) / 0.010 + 1e-6);
        h.CBgm[i] = 2.0 * std::sqrt(pp[c] * mm[c]) * (h.hi[i] - h.lo[i]) / w[c];
    }
}
std::vector<HIn> LoadF10(TFile* f, const char* sel, bool withRep, std::vector<std::vector<Vec>>* rep = nullptr, const char* var = "nom") {
    std::vector<HIn> H(kNH);
    const std::vector<double> E = F10Edges();
    for (int h = 0; h < kNH; ++h) {
        HIn& x = H[h];
        x.j = h / 2; x.e = h % 2;
        x.name = F10Name(x.j, x.e);
        x.cruOK = true; x.crkOK = (x.j == 4);
        x.lo.resize(kNB); x.hi.resize(kNB); x.x.resize(kNB); x.eg.resize(kNB);
        for (int i = 0; i < kNB; ++i) { x.lo[i] = E[i]; x.hi[i] = E[i + 1]; x.x[i] = TiltX(0.5 * (E[i] + E[i + 1])); x.eg[i] = h * kNB + i; }
        TH1* N = allio::Data(f, sel, "all", "N", x.name);
        TH1* CB = allio::Data(f, sel, "all", "CB", x.name);
        TH1* Npp = allio::Data(f, sel, "all", "Npp", x.name);
        TH1* Nmm = allio::Data(f, sel, "all", "Nmm", x.name);
        TH2* R = (TH2*)allio::Data(f, sel, "all", "CBrep", x.name);
        x.N.resize(kNB); x.CB.resize(kNB); x.vCB.resize(kNB); x.Npp.resize(kNB); x.Nmm.resize(kNB);
        for (int i = 0; i < kNB; ++i) {
            x.N[i] = N->GetBinContent(i + 1); x.CB[i] = CB->GetBinContent(i + 1);
            x.Npp[i] = Npp->GetBinContent(i + 1); x.Nmm[i] = Nmm->GetBinContent(i + 1);
            double s = 0, s2 = 0; const int nr = R->GetNbinsY();
            for (int k = 1; k <= nr; ++k) { const double v = R->GetBinContent(i + 1, k); s += v; s2 += v * v; }
            s /= nr; s2 /= nr;
            x.vCB[i] = std::max(0., s2 - s * s);
        }
        delete N; delete CB; delete Npp; delete Nmm; delete R;
        MakeCBgm(x);
        for (int c = 0; c < kNLC; ++c) {
            x.s[c].assign(kNB, 0.); x.v[c].assign(kNB, 0.);
            std::vector<int> comps;
            if (c == lGDom) comps = {kGpi0Dom, kGvecom, kGothom};
            else comps = {kLCcomp[c]};
            for (int cc : comps) {
                TH1* s = allio::Sim(f, sel, var, cc, x.name);
                for (int i = 0; i < kNB; ++i) { x.s[c][i] += s->GetBinContent(i + 1) * kNbToRaw; x.v[c][i] += std::pow(s->GetBinError(i + 1) * kNbToRaw, 2); }
                delete s;
            }
        }
        if (withRep && rep) {
            // replicas of the components used in the fit: GD, GC, Geta, C2..C5   [comp][hist][rep*kNB + bin]
            for (int c : {lGD, lGC, lGeta, lC2, lC3, lC4, lC5}) {
                TH2* s = (TH2*)allio::Sim(f, sel, "nom", kLCcomp[c], x.name, "_rep");
                const int nr = s->GetNbinsY();
                Vec& out = (*rep)[c][h];
                out.assign(nr * kNB, 0.);
                for (int k = 0; k < nr; ++k) for (int i = 0; i < kNB; ++i) out[k * kNB + i] = s->GetBinContent(i + 1, k + 1) * kNbToRaw;
                delete s;
            }
        }
    }
    return H;
}
// a histogram built from G3: M(pipi) in [pLo, pHi) minus the veto [vLo, vHi), M(ee) class e; F10 binning
HIn FromG3(TFile* f, const char* sel, const char* name, int j, int e, double pLo, double pHi, double vLo, double vHi,
           const HIn& ref, bool cruOK, bool crkOK) {
    HIn x;
    x.name = name; x.j = j; x.e = e; x.cruOK = cruOK; x.crkOK = crkOK;
    x.lo = ref.lo; x.hi = ref.hi; x.x = ref.x;
    x.eg.assign(kNB, -1);
    TH3* N = (TH3*)allio::Data(f, sel, "all", "N", "G3");
    TH3* CB = (TH3*)allio::Data(f, sel, "all", "CB", "G3");
    TH3* Npp = (TH3*)allio::Data(f, sel, "all", "Npp", "G3");
    TH3* Nmm = (TH3*)allio::Data(f, sel, "all", "Nmm", "G3");
    auto ysel = [&](int yb) {
        const double a = N->GetYaxis()->GetBinLowEdge(yb), b = N->GetYaxis()->GetBinUpEdge(yb);
        if (a < pLo - 1e-9 || b > pHi + 1e-9) return false;
        if (vHi > vLo && a >= vLo - 1e-9 && b <= vHi + 1e-9) return false;
        return true;
    };
    auto zsel = [&](int zb) { const double a = N->GetZaxis()->GetBinLowEdge(zb); return e == 0 ? a < kEsplit - 1e-9 : a >= kEsplit - 1e-9; };
    auto xmap = [&](int xb) { const double c = N->GetXaxis()->GetBinCenter(xb); for (int i = 0; i < kNB; ++i) if (c >= x.lo[i] && c < x.hi[i]) return i; return -1; };
    auto fill = [&](TH3* h3, Vec& out, Vec* var, double scale) {
        out.assign(kNB, 0.); if (var) var->assign(kNB, 0.);
        for (int xb = 1; xb <= h3->GetNbinsX(); ++xb) {
            const int i = xmap(xb); if (i < 0) continue;
            for (int yb = 1; yb <= h3->GetNbinsY(); ++yb) { if (!ysel(yb)) continue;
                for (int zb = 1; zb <= h3->GetNbinsZ(); ++zb) { if (!zsel(zb)) continue;
                    out[i] += h3->GetBinContent(xb, yb, zb) * scale;
                    if (var) (*var)[i] += std::pow(h3->GetBinError(xb, yb, zb) * scale, 2);
                } }
        }
    };
    Vec dummy;
    fill(N, x.N, nullptr, 1.); fill(CB, x.CB, nullptr, 1.); fill(Npp, x.Npp, nullptr, 1.); fill(Nmm, x.Nmm, nullptr, 1.);
    // G3 has no CB replicas: CB variance = replica variance of the reference F10 histogram x (CB / CB_ref)^2 per bin
    x.vCB.assign(kNB, 0.);
    for (int i = 0; i < kNB; ++i) x.vCB[i] = ref.CB[i] > 0 ? ref.vCB[i] * std::pow(x.CB[i] / ref.CB[i], 2) : 0.;
    MakeCBgm(x);
    for (int c = 0; c < kNLC; ++c) {
        x.s[c].assign(kNB, 0.); x.v[c].assign(kNB, 0.);
        std::vector<int> comps;
        if (c == lGDom) comps = {kGpi0Dom, kGvecom, kGothom}; else comps = {kLCcomp[c]};
        for (int cc : comps) {
            TH3* s = (TH3*)allio::Sim(f, sel, "nom", cc, "G3");
            Vec a, va; fill(s, a, &va, kNbToRaw);
            for (int i = 0; i < kNB; ++i) { x.s[c][i] += a[i]; x.v[c][i] += va[i]; }
            delete s;
        }
    }
    delete N; delete CB; delete Npp; delete Nmm;   // the axis lambdas use N: delete at the end
    return x;
}
void SmoothAll(HIn& x, bool omega = false) {
    SmRes d = Smooth(x.s[lGD], x.v[lGD], x.lo, x.hi);
    SmRes c = Smooth(x.s[lGC], x.v[lGC], x.lo, x.hi);
    x.D = d.f; x.C = c.f; x.domD = d.dom; x.domC = c.dom;
    for (int k = 0; k < 4; ++k) { x.pD[k] = d.p[k]; x.pC[k] = c.p[k]; }
    x.chiD = d.chi2; x.chiC = c.chi2; x.ndfD = d.ndf; x.ndfC = c.ndf;
    if (omega) {
        SmRes a = Smooth(x.s[lGDom], x.v[lGDom], x.lo, x.hi);
        SmRes b = Smooth(x.s[lGCom], x.v[lGCom], x.lo, x.hi);
        x.Dom = a.f; x.Com = b.f; x.domDom = a.dom; x.domCom = b.dom;
    } else { x.Dom.assign(x.lo.size(), 0.); x.Com.assign(x.lo.size(), 0.); }
}

// ------------------------------------------------------------------------------------------------
// helpers for windows (fractional F10 bins)
// ------------------------------------------------------------------------------------------------
inline double Frac(double lo, double hi, double a, double b) { const double o = std::min(hi, b) - std::max(lo, a); return o > 0 ? o / (hi - lo) : 0.; }
struct WinVal { double N = 0, CB = 0, B = 0, eta1 = 0, eta04 = 0, D = 0, C = 0, Geta = 0, C2 = 0, C3 = 0, C4 = 0, C5 = 0, T12 = 0; Vec gB; };
WinVal WinSum(const std::vector<HIn>& H, const StageA& A, double a, double b, int jsel = -1, int esel = -1) {
    WinVal w; w.gB.assign(A.m.np, 0.);
    Model m1 = A.m; m1.aEta = 1.0;
    for (const auto& h : H) {
        if (jsel >= 0 && h.j != jsel) continue;
        if (esel >= 0 && h.e != esel) continue;
        const double kE = A.r.th[h.e == 0 ? A.m.ikE1 : A.m.ikE2];
        for (int i = 0; i < (int)h.lo.size(); ++i) {
            const double f = Frac(h.lo[i], h.hi[i], a, b);
            if (f <= 0) continue;
            Vec g(A.m.np, 0.);
            const double B = Bkg(h, i, A.m, A.r.th, &g);
            w.B += f * B; for (int k = 0; k < A.m.np; ++k) w.gB[k] += f * g[k];
            w.N += f * h.N[i]; w.CB += f * kE * CBof(h, i, A.m);
            w.eta1 += f * EtaSum(h, i, m1); w.eta04 += f * EtaSum(h, i, A.m);
            w.D += f * h.D[i]; w.C += f * h.C[i];
            w.Geta += f * h.s[lGeta][i]; w.C2 += f * h.s[lC2][i]; w.C3 += f * h.s[lC3][i]; w.C4 += f * h.s[lC4][i]; w.C5 += f * h.s[lC5][i];
            w.T12 += f * h.s[lT12][i];
        }
    }
    return w;
}
inline double QF(const Vec& g, const Mat& C) { double s = 0; for (size_t k = 0; k < g.size(); ++k) for (size_t l = 0; l < g.size(); ++l) s += g[k] * C[k][l] * g[l]; return s; }

std::string Fmt(double v, double e, int prec = 4) { return std::string(Form("%.*f +- %.*f", prec, v, prec, e)); }

// ------------------------------------------------------------------------------------------------
// helpers of the main part
// ------------------------------------------------------------------------------------------------
inline bool Used(const HIn& h, const Model& m) { return m.eOnly < 0 || h.e == m.eOnly; }
// GLS fit with the MC covariance of the (smoothed) templates: V = V_MC(theta) + diag(mu + k^2 var CB); V_MC from the simulation
// replicas HR at the current parameters (recomputed in 3 outer iterations). Used for the D11 decisions and the chi2 criterion.
struct GLSRes { bool ok = false; Vec th; double chi2 = 0; int ndf = 0; Mat C; };
inline double MuStageA(const HIn& h, int i, const Model& m, const Vec& th, Vec* g) {
    const double kE = th[h.e == 0 ? m.ikE1 : m.ikE2], cb = CBof(h, i, m);
    if (g) { std::fill(g->begin(), g->end(), 0.); (*g)[h.e == 0 ? m.ikE1 : m.ikE2] += cb; }
    return kE * cb + Bkg(h, i, m, th, g) + EtaSum(h, i, m);
}
GLSRes GLS(const std::vector<HIn>& H, const StageA& A, const std::vector<std::vector<HIn>>& HR) {
    GLSRes o;
    const Model& m = A.m;
    const int n = (int)A.bins.size(), np = m.np;
    if ((int)A.r.th.size() != np || HR.size() < 10) return o;
    std::vector<int> fr;
    for (int k = 0; k < np; ++k) if (!(k == m.ikE1 && (m.eOnly == 1 || !m.kE1free))) fr.push_back(k);
    const int nf = (int)fr.size();
    std::vector<int> cK; Vec cV, cS;
    cK.push_back(m.ikE2); cV.push_back(1.0); cS.push_back(kKE2sig);
    if (m.kBsig > 0) { cK.push_back(m.ikB); cV.push_back(1.0); cS.push_back(m.kBsig); }
    Vec th = A.r.th;
    Vec y(n);
    for (int k = 0; k < n; ++k) y[k] = H[A.bins[k].first].N[A.bins[k].second];
    TMatrixDSym Vi(n);
    Mat Hinv;
    for (int outer = 0; outer < 3; ++outer) {
        std::vector<Vec> mus(HR.size(), Vec(n, 0.));
        for (size_t r = 0; r < HR.size(); ++r) for (int k = 0; k < n; ++k) mus[r][k] = MuStageA(HR[r][A.bins[k].first], A.bins[k].second, m, th, nullptr);
        Vec mean(n, 0.);
        for (auto& v : mus) for (int k = 0; k < n; ++k) mean[k] += v[k] / mus.size();
        TMatrixDSym Vmc(n);
        for (int a = 0; a < n; ++a) for (int b = 0; b <= a; ++b) {
            double s2 = 0; for (auto& v : mus) s2 += (v[a] - mean[a]) * (v[b] - mean[b]);
            Vmc(a, b) = Vmc(b, a) = s2 / (mus.size() - 1);
        }
        for (int it = 0; it < 40; ++it) {
            Vec mu(n); Mat J(n, Vec(np, 0.)); Vec g(np);
            TMatrixDSym V = Vmc;
            for (int k = 0; k < n; ++k) {
                const HIn& h = H[A.bins[k].first]; const int i = A.bins[k].second;
                mu[k] = MuStageA(h, i, m, th, &g); J[k] = g;
                const double kE = th[h.e == 0 ? m.ikE1 : m.ikE2];
                V(k, k) += std::max(mu[k], 0.25) + kE * kE * h.vCB[i];
            }
            Vi = V; Vi.Invert();
            Mat Hm = bla::Zero(nf, nf); Vec b(nf, 0.);
            Mat VJ(n, Vec(nf, 0.));
            for (int a = 0; a < n; ++a) for (int q = 0; q < nf; ++q) { double s2 = 0; for (int c = 0; c < n; ++c) s2 += Vi(a, c) * J[c][fr[q]]; VJ[a][q] = s2; }
            for (int q = 0; q < nf; ++q) {
                for (int a = 0; a < n; ++a) b[q] += VJ[a][q] * (y[a] - mu[a]);
                for (int w = 0; w < nf; ++w) { double s2 = 0; for (int a = 0; a < n; ++a) s2 += J[a][fr[q]] * VJ[a][w]; Hm[q][w] = s2; }
            }
            for (size_t c = 0; c < cK.size(); ++c) for (int q = 0; q < nf; ++q) if (fr[q] == cK[c]) { Hm[q][q] += 1 / (cS[c] * cS[c]); b[q] += (cV[c] - th[cK[c]]) / (cS[c] * cS[c]); }
            Hinv = Hm;
            if (!bla::CholInv(Hinv)) return o;
            double dmax = 0;
            for (int q = 0; q < nf; ++q) { double d = 0; for (int w = 0; w < nf; ++w) d += Hinv[q][w] * b[w]; th[fr[q]] += d; dmax = std::max(dmax, std::fabs(d) / std::sqrt(Hinv[q][q])); }
            if (dmax < 1e-6) break;
        }
    }
    // final chi2 with the last covariance
    Vec r(n);
    for (int k = 0; k < n; ++k) r[k] = y[k] - MuStageA(H[A.bins[k].first], A.bins[k].second, m, th, nullptr);
    double c2 = 0;
    for (int a = 0; a < n; ++a) for (int b = 0; b < n; ++b) c2 += r[a] * Vi(a, b) * r[b];
    for (size_t c = 0; c < cK.size(); ++c) c2 += std::pow((th[cK[c]] - cV[c]) / cS[c], 2);
    o.ok = true; o.th = th; o.chi2 = c2; o.ndf = n - nf + (int)cK.size();
    o.C = bla::Zero(np, np);
    for (int q = 0; q < nf; ++q) for (int w = 0; w < nf; ++w) o.C[fr[q]][fr[w]] = Hinv[q][w];
    return o;
}
struct NestOut { std::vector<StageA> M; std::vector<GLSRes> G; int adopt = 1, adoptData = 1; double d01 = 0, d12 = 0, dps = 0, g01 = 0, g12 = 0, gps = 0; };
NestOut RunNested(const std::vector<HIn>& H, const CRdef& cr, const Model& base, const TMatrixDSym* eg, const std::vector<std::vector<HIn>>* HR = nullptr) {
    NestOut o;
    for (int k = 0; k < 5; ++k) {
        Model m = base;
        m.tilt = (k != 0); m.quad = (k == 2 || k == 4); m.perSlice = (k >= 3);
        o.M.push_back(FitStageA(H, m, cr, eg));
        if (HR) o.G.push_back(GLS(H, o.M.back(), *HR));
    }
    // D11 with the data-only chi2 (literal)
    o.d01 = o.M[0].r.chi2P - o.M[1].r.chi2P;
    o.d12 = o.M[1].r.chi2P - o.M[2].r.chi2P;
    o.adoptData = o.d12 > 4 ? 2 : 1;
    o.dps = o.M[o.adoptData].r.chi2P - o.M[o.adoptData == 1 ? 3 : 4].r.chi2P;
    if (o.dps > 4) o.adoptData = o.adoptData == 1 ? 3 : 4;
    o.adopt = o.adoptData;
    // D11 with the chi2 including the MC covariance of the templates (used when the replicas are given)
    if (HR && o.G.size() == 5 && o.G[0].ok && o.G[1].ok && o.G[2].ok && o.G[3].ok && o.G[4].ok) {
        o.g01 = o.G[0].chi2 - o.G[1].chi2;
        int a = o.g01 > 4 ? 1 : 0;
        if (a == 1) {
            o.g12 = o.G[1].chi2 - o.G[2].chi2;
            if (o.g12 > 4) a = 2;
            o.gps = o.G[a].chi2 - o.G[a == 1 ? 3 : 4].chi2;
            if (o.gps > 4) a = a == 1 ? 3 : 4;
        }
        o.adopt = a;
    }
    return o;
}
// A1 in the CR bins of a configuration
struct A1Bin { std::string h; double lo, hi, sraw, ssm, sfit, neff, sig, geta; bool k; };
struct A1Out { double maxRaw = 0, maxSm = 0, maxFit = 0, sigU = 0, getaU = 0, totU = 0, sigK = 0, getaK = 0, totK = 0, maxGeta = 0;
               int nRaw = 0, nSm = 0, nFit = 0, nSmEff = 0, nBins = 0; std::vector<A1Bin> bins; };
A1Out EvalA1(const std::vector<HIn>& H, const StageA& A) {
    A1Out o;
    for (const auto& h : H) {
        if (!Used(h, A.m)) continue;
        for (int i = 0; i < (int)h.lo.size(); ++i) {
            if (!InCR(h, i, A.cr)) continue;
            const double sig = h.s[lC2][i] + h.s[lC3][i] + h.s[lC4][i] + h.s[lC5][i];
            const double vs = h.v[lC2][i] + h.v[lC3][i] + h.v[lC4][i] + h.v[lC5][i];
            const double t12 = h.s[lT12][i], tsm = h.D[i] + h.C[i] + h.s[lGeta][i] + sig;
            if (!(t12 > 0) && !(tsm > 0)) continue;
            const bool isK = !(h.cruOK && h.lo[i] >= A.cr.uLo - 1e-9 && h.hi[i] <= A.cr.uHi + 1e-9);
            const double mu = A.r.th[h.e == 0 ? A.m.ikE1 : A.m.ikE2] * CBof(h, i, A.m) + Bkg(h, i, A.m, A.r.th, nullptr) + EtaSum(h, i, A.m);
            A1Bin b{h.name, h.lo[i], h.hi[i], t12 > 0 ? sig / t12 : 0., tsm > 0 ? sig / tsm : 0., mu > 0 ? A.m.aEta * sig / mu : 0.,
                    vs > 0 ? sig * sig / vs : 0., sig, tsm > 0 ? h.s[lGeta][i] / tsm : 0., isK};
            o.bins.push_back(b); ++o.nBins;
            o.maxRaw = std::max(o.maxRaw, b.sraw); o.maxSm = std::max(o.maxSm, b.ssm); o.maxFit = std::max(o.maxFit, b.sfit);
            o.maxGeta = std::max(o.maxGeta, b.geta);
            if (b.sraw > 0.01) ++o.nRaw;
            if (b.ssm > 0.01) ++o.nSm;
            if (b.sfit > 0.01) ++o.nFit;
            if (b.ssm > 0.01 && b.neff > 2) ++o.nSmEff;
            if (isK) { o.sigK += sig; o.getaK += h.s[lGeta][i]; o.totK += tsm; } else { o.sigU += sig; o.getaU += h.s[lGeta][i]; o.totU += tsm; }
        }
    }
    return o;
}
// A2, template version: CR-K bins (S5) of the used M(ee) classes, N = k CB + c' B_A + a_eta(0) eta + h step
struct A2T { double h = 0, eh = 0, c = 0, ec = 0, chi2 = 0; int ndf = 0; bool ok = false; };
A2T A2Template(const std::vector<HIn>& HH, const StageA& A, bool gm, double kLo, double kHi) {
    A2T out;
    const double eK = Edge(kEdgeLimit, 0.);
    std::vector<std::pair<int, int>> bins;
    GProb Pp; Pp.np = 2; Pp.th0 = {1.0, 0.0};
    for (size_t h = 0; h < HH.size(); ++h) {
        if (!HH[h].crkOK || !Used(HH[h], A.m)) continue;
        for (int i = 0; i < kNB; ++i) {
            if (HH[h].lo[i] < kLo - 1e-9 || HH[h].hi[i] > kHi + 1e-9) continue;
            const double cb = gm ? HH[h].CBgm[i] : HH[h].CB[i];
            if (!(HH[h].N[i] > 0 || cb > 0 || HH[h].D[i] > 0 || HH[h].C[i] > 0)) continue;
            bins.push_back(std::make_pair((int)h, i)); Pp.y.push_back(HH[h].N[i]); Pp.eg.push_back(-1);
        }
    }
    if (bins.size() < 4) return out;
    const Model m = A.m; const Vec th = A.r.th;
    Pp.eval = [&HH, bins, m, th, gm, eK](const Vec& p, int k, double& mu, Vec& g, double& vex) {
        const HIn& h = HH[bins[k].first]; const int i = bins[k].second;
        double ve = 0; const double eta = EtaSum(h, i, m, &ve);
        const double kE = th[h.e == 0 ? m.ikE1 : m.ikE2];
        const double cb = gm ? h.CBgm[i] : h.CB[i];
        const double B = Bkg(h, i, m, th, nullptr);
        const double st = bla::StepBin(h.lo[i], h.hi[i], eK, kStepS) * (h.hi[i] - h.lo[i]) / 0.0025;
        mu = kE * cb + p[0] * B + eta + p[1] * st;
        g[0] = B; g[1] = st;
        vex = ve + kE * kE * h.vCB[i];
    };
    Pp.sandwich = false;
    GRes r = GFit(Pp);
    out.ok = r.ok; out.c = r.th[0]; out.ec = r.ErrM(0); out.h = r.th[1]; out.eh = r.ErrM(1); out.chi2 = r.chi2P; out.ndf = r.ndf;
    return out;
}
// prediction in a set of bins (closure): sum of N, model and the variance of the model from the stage-A covariance
struct Pred { double N = 0, mu = 0, vfit = 0, vcb = 0, B = 0, eta = 0, cb = 0; };
Pred Predict(const std::vector<HIn>& H, const StageA& A, std::function<bool(const HIn&, int)> sel, const Mat* C = nullptr) {
    Pred p; Vec g(A.m.np, 0.);
    for (const auto& h : H) {
        if (!Used(h, A.m)) continue;
        for (int i = 0; i < (int)h.lo.size(); ++i) {
            if (!sel(h, i)) continue;
            Vec gg(A.m.np, 0.);
            const double b = Bkg(h, i, A.m, A.r.th, &gg);
            const double kE = A.r.th[h.e == 0 ? A.m.ikE1 : A.m.ikE2];
            p.N += h.N[i]; p.B += b; p.eta += EtaSum(h, i, A.m); p.cb += kE * CBof(h, i, A.m); p.vcb += kE * kE * h.vCB[i];
            for (int k = 0; k < A.m.np; ++k) g[k] += gg[k];
        }
    }
    p.mu = p.B + p.eta + p.cb;
    p.vfit = QF(g, C ? *C : A.r.Cs);
    return p;
}
std::map<std::string, std::vector<double>> ReadRowsSummary(const char* path) {
    std::map<std::string, std::vector<double>> m;
    FILE* f = fopen(path, "r");
    if (!f) return m;
    char line[1024];
    while (fgets(line, sizeof line, f)) {
        char sel[32]; double lo, hi, e, ho, eho, co, hs, ehs, cs, nb; int no, ns;
        if (sscanf(line, "A2row %31s %lf %lf %lf own %lf %lf %lf %d spread %lf %lf %lf %d nbelow %lf", sel, &lo, &hi, &e, &ho, &eho, &co, &no, &hs, &ehs, &cs, &ns, &nb) == 13)
            m[Form("A2 %s %.2f-%.2f", sel, lo, hi)] = {e, ho, eho, co, (double)no, hs, ehs, cs, (double)ns, nb};
        int k; double nr, T, eT, sob, c2; int nd;
        if (sscanf(line, "TAG %31s %d %lf %lf %lf %lf %lf %d", sel, &k, &nr, &T, &eT, &sob, &c2, &nd) == 8) m[Form("TAG %s %d", sel, k)] = {nr, T, eT, sob, c2, (double)nd};
    }
    fclose(f);
    return m;
}

}  // namespace wp3

// ================================================================================================
void wp3_stageA_cnt_plain2(double aEta0 = 0.4, const char* tag = "") {
    using namespace alldefs;
    using namespace wp3;
    kAeta0 = aEta0;
    const bool full = (tag == nullptr || tag[0] == 0);   // full run: tables and figures; tagged run: only stageA_<tag>.txt and its log
    const std::string sfx = full ? std::string("") : std::string("_") + tag;
    gSystem->mkdir(kOutDir, kTRUE);
    gOut = fopen(Form("%s/wp3_stageA_log%s.txt", kCntOut, sfx.c_str()), "w");
    TFile* f = TFile::Open(kHistFile);
    TFile* fr = TFile::Open(Form("%s/wp3_rows.root", kOutDir));
    TMatrixDSym* egS0 = (TMatrixDSym*)fr->Get("egcov_S0");
    TMatrixDSym* egPre = (TMatrixDSym*)fr->Get("egcov_pre");
    std::map<std::string, std::vector<double>> RS = ReadRowsSummary(Form("%s/wp3_rows_summary.txt", kOutDir));
    P("wp3_stageA.C  WP3 control regions, stage A. hists: %s\n", kHistFile);
    P("L = %.1f nb^-1, sim templates in expected raw rows = nb x %.2f; a_eta(0) = %.2f; k_E2 = 1 +- %.2f\n\n", kL, kNbToRaw, kAeta0, kKE2sig);

    // -------------------------------------------------------------------------------------------
    // 1. inputs and smoothing
    // -------------------------------------------------------------------------------------------
    std::vector<std::vector<Vec>> rep(kNLC, std::vector<Vec>(kNH));
    std::vector<HIn> H = LoadF10(f, "S0", true, &rep);
    std::vector<HIn> Hpre = LoadF10(f, "pre", false);
    std::vector<HIn> Hgh = LoadF10(f, "S0", false, nullptr, "ghostA");
    for (auto& h : H) SmoothAll(h, true);
    for (auto& h : Hpre) SmoothAll(h, true);
    for (auto& h : Hgh) SmoothAll(h, true);
    P("1. Smoothing of the non-eta templates, S0 (exp(pol3), t = (M4 - 0.55)/0.15; chi2 with the mean weight of the histogram)\n");
    P("   %-9s | %-50s | %-50s\n", "hist", "G_D: start, sum raw / smooth, chi2/ndf (p)", "G_C: start, sum raw / smooth, chi2/ndf (p)");
    for (const auto& h : H) {
        double sd = 0, fd = 0, sc = 0, fc = 0;
        for (int i = 0; i < kNB; ++i) { sd += h.s[lGD][i]; fd += h.D[i]; sc += h.s[lGC][i]; fc += h.C[i]; }
        P("   %-9s | %.4f  %9.1f / %9.1f  %6.1f/%-3d (p %.3f) | %.4f  %9.1f / %9.1f  %6.1f/%-3d (p %.3f)\n", h.name.c_str(), h.lo[h.domD], sd, fd,
          h.chiD, h.ndfD, TMath::Prob(h.chiD, h.ndfD), h.lo[h.domC], sc, fc, h.chiC, h.ndfC, TMath::Prob(h.chiC, h.ndfC));
    }
    auto nestPrint = [&](const char* title, const NestOut& n) {
        P("%s\n", title);
        P("   %-48s | %8s %8s %4s %4s %7s | %8s %7s |\n", "model", "chi2_P", "chi2_BC", "nbin", "ndf", "/ndf", "chi2_MC", "/ndf");
        for (size_t k = 0; k < n.M.size(); ++k) {
            const StageA& A = n.M[k];
            const bool g = k < n.G.size() && n.G[k].ok;
            P("   %-48s | %8.1f %8.1f %4d %4d %7.3f | %8.1f %7.3f | %s", (A.m.Tag() + ((int)k == n.adopt ? " [D11]" : "") + ((int)k == n.adoptData && n.adoptData != n.adopt ? " [D11 data-only]" : "")).c_str(),
              A.r.chi2P, A.r.chi2BC, A.r.n, A.r.ndf, A.r.chi2P / A.r.ndf, g ? n.G[k].chi2 : 0., g ? n.G[k].chi2 / n.G[k].ndf : 0., A.r.ok ? "" : "[NOT CONVERGED] ");
            for (int p = 0; p < A.m.np; ++p) P("%s %.4f(%.4f) ", A.m.names[p].c_str(), A.r.th[p], A.r.Err(p));
            P("\n");
        }
        P("   D11 with the data-only chi2: dchi2(+ tilt) = %.1f, (+ quadratic) = %.1f, (+ per-slice tilts) = %.1f -> %s\n", n.d01, n.d12, n.dps, n.M[n.adoptData].m.Tag().c_str());
        if (!n.G.empty()) P("   D11 with the MC-inclusive chi2 (GLS): dchi2(+ tilt) = %.1f, (+ quadratic) = %.1f, (+ per-slice tilts) = %.1f -> %s\n", n.g01, n.g12, n.gps, n.M[n.adopt].m.Tag().c_str());
    };

    // -------------------------------------------------------------------------------------------
    // 2. step 0: the specification as written (E1 + E2, CR-U 0.565-0.700 + CR-K)
    // -------------------------------------------------------------------------------------------
    CRdef crN;
    CRdef crFB1; crFB1.uLo = 0.575; crFB1.tag = "CR-U 0.575-0.700 + CR-K (S5, M(pipi) >= 0.44, 0.400-0.5225)";
    CRdef crFB2; crFB2.uLo = 0.575; crFB2.useK = false; crFB2.tag = "CR-U 0.575-0.700 (no CR-K)";
    CRdef crU565; crU565.useK = false; crU565.tag = "CR-U 0.565-0.700 (no CR-K)";
    const NestOut n0 = RunNested(H, crN, Model(), egS0);
    nestPrint("\n2. Step 0, the specification as written: S0, E1 + E2, CR-U 0.565-0.700 + CR-K, a_eta(0) = 0.4", n0);
    const StageA& s0 = n0.M[n0.adopt];
    const double kE1v = s0.r.th[s0.m.ikE1], kE1e = s0.r.Err(s0.m.ikE1);
    const bool fbKE1 = kE1v < 0.8 || kE1v > 3.0;
    P("   k_E1 = %.2f +- %.2f: outside 0.8-3 -> %s\n", kE1v, kE1e, fbKE1 ? "fall-back 'k_E1 outside 0.8-3: drop E1; kappa = 1 +- 0.2' FIRES" : "no fall-back");
    // diagnosis: the same model with the ghost variant (a) templates
    const StageA sGh = FitStageA(Hgh, s0.m, crN, egS0);
    P("   diagnosis: same model with the ghost-variant (a) templates (ghost-flagged true pairs in T2): k_E1 = %.2f +- %.2f, chi2/ndf %.1f/%d\n",
      sGh.r.th[sGh.m.ikE1], sGh.r.Err(sGh.m.ikE1), sGh.r.chi2P, sGh.r.ndf);
    {
        double t12[2] = {0, 0}, t3[2] = {0, 0};
        for (int h = 0; h < kNH; ++h) {
            TH1* a = allio::Sim(f, "S0", "nom", kT3, H[h].name);
            for (int i = 0; i < kNB; ++i) if (H[h].lo[i] >= 0.565 - 1e-9) { t12[H[h].e] += H[h].s[lT12][i]; t3[H[h].e] += a->GetBinContent(i + 1) * kNbToRaw; }
            delete a;
        }
        P("   in CR-U the fake-lepton template T3 (not in the model, D6) is %.1f %% of T12 in E1 and %.1f %% in E2\n", 100 * t3[0] / t12[0], 100 * t3[1] / t12[1]);
    }

    // -------------------------------------------------------------------------------------------
    // 3. after the k_E1 fall-back: E2 only, kappa_B = 1 +- 0.2; the CR chain of A1/A2
    // -------------------------------------------------------------------------------------------
    Model base;
    if (fbKE1) { base.eOnly = 1; base.kBsig = 0; }   // spec: 'drop E1; kappa = 1 +- 0.2' (kappa = signal conversion factor, fixed in stage A); kappa_B stays free
    // histogram sets: H (F10), H44 (S5 split at M(pipi) 0.44 from G3 for the CR-K >= 0.44 fall-back)
    // fall-back 1 set: the ten F10 histograms (S5 used for CR-U only) + S5 with M(pipi) >= 0.44 from G3 used for CR-K only (same c_5, b_5)
    std::vector<HIn> H44 = H;
    {
        H44[8].crkOK = false; H44[9].crkOK = false;
        for (int e = 0; e < 2; ++e) {
            HIn hi = FromG3(f, "S0", Form("S5[0.44,inf)E%d", e + 1), 4, e, 0.44, 2.0, 0, 0, H[8 + e], false, true);
            SmoothAll(hi, true);
            H44.push_back(hi);
        }
    }
    // simulation replicas of the histogram sets (templates replaced by replica k, smoothing redone; the omega parts keep their nominal share)
    auto replicaH = [&](const std::vector<HIn>& H0, int k) {
        std::vector<HIn> Hr = H0;
        for (size_t h = 0; h < Hr.size(); ++h) {
            if (Hr[h].eg.empty() || Hr[h].eg[0] < 0) continue;   // G3-built histograms keep the nominal templates
            const int hf = Hr[h].eg[0] / kNB;
            HIn& x = Hr[h];
            for (int c : {lGD, lGC, lGeta, lC2, lC3, lC4, lC5}) for (int i = 0; i < kNB; ++i) x.s[c][i] = rep[c][hf][k * kNB + i];
            x.D = Smooth(x.s[lGD], x.v[lGD], x.lo, x.hi, H0[h].domD, H0[h].pD).f;
            x.C = Smooth(x.s[lGC], x.v[lGC], x.lo, x.hi, H0[h].domC, H0[h].pC).f;
            for (int i = 0; i < kNB; ++i) {
                x.Dom[i] = H0[h].D[i] > 0 ? H0[h].Dom[i] * x.D[i] / H0[h].D[i] : 0.;
                x.Com[i] = H0[h].C[i] > 0 ? H0[h].Com[i] * x.C[i] / H0[h].C[i] : 0.;
            }
        }
        return Hr;
    };
    std::vector<std::vector<HIn>> HRH, HR44;
    for (int k = 0; k < kNRepSim; ++k) { HRH.push_back(replicaH(H, k)); HR44.push_back(replicaH(H44, k)); }
    struct CfgRes { std::string tag; CRdef cr; NestOut n; A1Out a1; bool a1fail = false, a2fail = false, hasK = false; A2T a2t, a2tg;
                    double a2h[2] = {0, 0}, a2e[2] = {0, 0}, a2c[2] = {0, 0}; int a2n[2] = {0, 0}; std::string a2sl[2]; double tagT = 0, tagE = 0; };
    std::vector<CfgRes> chain(3);
    chain[0].tag = "spec CRs"; chain[0].cr = crN;
    chain[1].tag = "fall-back 1"; chain[1].cr = crFB1;
    chain[2].tag = "fall-back 2"; chain[2].cr = crFB2;
    for (int c = 0; c < 3; ++c) {
        CfgRes& R = chain[c];
        const std::vector<HIn>& HH = c == 1 ? H44 : H;
        R.n = RunNested(HH, R.cr, base, c == 1 ? nullptr : egS0, c == 1 ? &HR44 : &HRH);
        const StageA& A = R.n.M[R.n.adopt];
        R.a1 = EvalA1(HH, A);
        R.a1fail = R.a1.maxRaw >= 0.01;
        R.hasK = R.cr.useK;
        if (R.hasK) {
            const char* sl0 = c == 0 ? "0.43-0.46" : "0.44-0.46";
            const char* sls[2] = {sl0, "0.46-0.55"};
            for (int k = 0; k < 2; ++k) {
                R.a2sl[k] = sls[k];
                const std::vector<double>& v = RS[Form("A2 S0 %s", sls[k])];
                if (v.size() >= 10) { R.a2h[k] = v[5]; R.a2e[k] = v[6]; R.a2c[k] = v[7]; R.a2n[k] = (int)v[8]; }
                if (std::fabs(R.a2h[k] / R.a2e[k]) >= 2) R.a2fail = true;
            }
            const std::vector<double>& t = RS["TAG S0 0"];
            if (t.size() >= 3) { R.tagT = t[1]; R.tagE = t[2]; if (std::fabs(R.tagT / R.tagE) >= 2) R.a2fail = true; }
            R.a2t = A2Template(HH, A, false, 0.400, 0.5225);
            R.a2tg = A2Template(HH, A, true, 0.400, 0.5225);
        }
        char t1[256];
        snprintf(t1, sizeof t1, "\n3.%d %s: %s, %s", c + 1, R.tag.c_str(), R.cr.tag.c_str(), fbKE1 ? "E2 only (k_E1 fall-back)" : "E1 + E2");
        nestPrint(t1, R.n);
        P("   A1 (sim C2-C5 / raw T12 at SMASH < 1 %% in every CR bin): max %.2f %% (%d of %d bins above 1 %%); smoothed denominator: max %.2f %% (%d bins, %d with N_eff > 2); "
          "relative to the fitted model at a_eta(0): max %.2f %% (%d bins) -> %s\n", 100 * R.a1.maxRaw, R.a1.nRaw, R.a1.nBins, 100 * R.a1.maxSm, R.a1.nSm, R.a1.nSmEff,
          100 * R.a1.maxFit, R.a1.nFit, R.a1fail ? "FAIL" : "PASS");
        for (const auto& b : R.a1.bins)
            if (b.sraw > 0.01 || b.ssm > 0.01 || b.sfit > 0.01)
                P("      %-18s %.4f-%.4f %s raw %.2f %% smooth %.2f %% fit %.2f %% N_eff(signal) %.1f signal %.1f rows\n", b.h.c_str(), b.lo, b.hi, b.k ? "CR-K" : "CR-U",
                  100 * b.sraw, 100 * b.ssm, 100 * b.sfit, b.neff, b.sig);
        P("      integrated: CR-U C2-C5 %.3f %%, G_eta %.2f %%", 100 * R.a1.sigU / std::max(R.a1.totU, 1e-9), 100 * R.a1.getaU / std::max(R.a1.totU, 1e-9));
        if (R.hasK) P("; CR-K C2-C5 %.3f %%, G_eta %.2f %%", 100 * R.a1.sigK / std::max(R.a1.totK, 1e-9), 100 * R.a1.getaK / std::max(R.a1.totK, 1e-9));
        P("; max G_eta share in a bin %.1f %%\n", 100 * R.a1.maxGeta);
        if (R.hasK) {
            P("   A2 (row-level, all M(ee), CB spread): h(%s) = %.2f +- %.2f (%.1f sigma, chi2/ndf %.1f/%d), h(%s) = %.2f +- %.2f (%.1f sigma, chi2/ndf %.1f/%d); "
              "ECAL tag in CR-K T = %.1f +- %.1f -> %s\n", R.a2sl[0].c_str(), R.a2h[0], R.a2e[0], R.a2h[0] / R.a2e[0], R.a2c[0], R.a2n[0], R.a2sl[1].c_str(), R.a2h[1], R.a2e[1],
              R.a2h[1] / R.a2e[1], R.a2c[1], R.a2n[1], R.tagT, R.tagE, R.a2fail ? "FAIL" : "PASS");
            P("   A2 template version (CR-K bins of the fit, stage-A shape + step): CB spread h = %.2f +- %.2f (%.1f sigma), c' = %.3f +- %.3f, chi2/ndf %.1f/%d;"
              " CB per-cell h = %.2f +- %.2f (%.1f sigma), c' = %.3f, chi2/ndf %.1f/%d\n", R.a2t.h, R.a2t.eh, R.a2t.h / R.a2t.eh, R.a2t.c, R.a2t.ec, R.a2t.chi2, R.a2t.ndf,
              R.a2tg.h, R.a2tg.eh, R.a2tg.h / R.a2tg.eh, R.a2tg.c, R.a2tg.chi2, R.a2tg.ndf);
        }
    }
    // the chain of Sect. 6: A1 or A2 fails -> fall-back 1; still fails -> fall-back 2 (drop CR-K, b from CR-U only)
    int adoptCfg = 0;
    if (chain[0].a1fail || chain[0].a2fail) adoptCfg = 1;
    if (adoptCfg == 1 && (chain[1].a1fail || chain[1].a2fail)) adoptCfg = 2;
    P("\n   Fall-back chain (A1/A2): spec CRs %s; fall-back 1 %s; -> adopted: %s (%s)\n", (chain[0].a1fail || chain[0].a2fail) ? "FAIL" : "PASS",
      (chain[1].a1fail || chain[1].a2fail) ? "FAIL" : "PASS", chain[adoptCfg].tag.c_str(), chain[adoptCfg].cr.tag.c_str());
    if (adoptCfg == 2)
        P("   -> CR-K dropped, the tilt comes from CR-U only: the bump background is EXTRAPOLATED; the CR-L closure term enters S10.\n");

    // -------------------------------------------------------------------------------------------
    // 4. adopted configuration, chi2 fall-back (omega split, then error scaling)
    // -------------------------------------------------------------------------------------------
    const CRdef crA = chain[adoptCfg].cr;
    std::vector<HIn>& HA = adoptCfg == 1 ? H44 : H;
    const TMatrixDSym* egA = adoptCfg == 1 ? nullptr : egS0;
    const std::vector<std::vector<HIn>>& HRA = adoptCfg == 1 ? HR44 : HRH;
    StageA nom = chain[adoptCfg].n.M[chain[adoptCfg].n.adopt];
    GLSRes nomG = chain[adoptCfg].n.G[chain[adoptCfg].n.adopt];
    StageA nomOmega;
    GLSRes nomOmegaG;
    bool omegaAdopted = false;
    double errScale = 1.0;
    {
        Model mo = nom.m; mo.omega = true;
        nomOmega = FitStageA(HA, mo, crA, egA);
        nomOmegaG = GLS(HA, nomOmega, HRA);
        const double r0 = nomG.chi2 / nomG.ndf;
        P("\n4. Adopted %s, %s: chi2/ndf with the MC template covariance (GLS) = %.1f/%d = %.3f; data-only Pearson %.1f/%d = %.3f\n", nom.m.Tag().c_str(),
          crA.tag.c_str(), nomG.chi2, nomG.ndf, r0, nom.r.chi2P, nom.r.ndf, nom.r.chi2P / nom.r.ndf);
        P("   omega split on top (information unless needed): GLS dchi2 = %.1f (data-only %.1f), c_om = %.3f +- %.3f (stat)\n", nomG.chi2 - nomOmegaG.chi2,
          nom.r.chi2P - nomOmega.r.chi2P, nomOmega.r.th[nomOmega.m.iom], nomOmega.r.Err(nomOmega.m.iom));
        if (r0 > 1.5) {
            P("   chi2/ndf > 1.5 after the D11 terms: fall-back 'add the omega split'\n");
            if (nomG.chi2 - nomOmegaG.chi2 > 4) { nom = nomOmega; nomG = nomOmegaG; omegaAdopted = true; P("   omega split adopted (dchi2 > 4)\n"); }
            const double r1 = nomG.chi2 / nomG.ndf;
            if (r1 > 1.5) { errScale = std::sqrt(r1); P("   still > 1.5: stat errors x sqrt(chi2/ndf) = x %.3f, 'background-model limited'\n", errScale); }
        }
    }
    const Model mN = nom.m;
    // per-histogram chi2 and residual pattern
    {
        P("   chi2 per histogram (data-only Pearson) and mean pull in M4 thirds of the CR (%zu fit bins, %zu mu):\n", nom.bins.size(), nom.r.mu.size());
        std::map<int, double> c2; std::map<int, int> nb; std::map<int, std::vector<double>> pt;
        for (size_t k = 0; k < nom.bins.size(); ++k) {
            const int h = nom.bins[k].first, i = nom.bins[k].second;
            const double y = HA[h].N[i], mu = nom.r.mu[k], V = nom.r.V[k];
            c2[h] += (y - mu) * (y - mu) / V; nb[h] += 1;
            const double x = 0.5 * (HA[h].lo[i] + HA[h].hi[i]);
            const int t = x < 0.615 ? 0 : (x < 0.655 ? 1 : 2);
            pt[h * 3 + t].push_back((y - mu) / std::sqrt(V));
        }
        for (auto& kv : c2) {
            P("      %-18s chi2 %7.1f / %3d bins", HA[kv.first].name.c_str(), kv.second, nb[kv.first]);
            for (int t = 0; t < 3; ++t) { const auto& v = pt[kv.first * 3 + t]; double m = 0; for (double z : v) m += z; if (!v.empty()) m /= v.size(); P("  <pull>%s %+.2f (%zu)", t == 0 ? "[<0.615]" : (t == 1 ? "[0.615-0.655]" : "[>0.655]"), m, v.size()); }
            P("\n");
        }
    }

    // -------------------------------------------------------------------------------------------
    // 5. MC error (100 simulation replicas, smoothing redone) and CB replica check
    // -------------------------------------------------------------------------------------------
    std::vector<Vec> thMC;
    const std::vector<std::vector<HIn>>& HR = adoptCfg == 1 ? HR44 : HRH;   // replica histogram sets of the adopted configuration
    for (int k = 0; k < kNRepSim; ++k) thMC.push_back(FitStageA(HR[k], mN, crA, egA, 0, 0, false, &nom.r.th).r.th);
    gChi2MC = nomG.chi2;
    P("   GLS (MC-inclusive) estimates of the adopted model, for comparison:");
    for (int p = 0; p < mN.np; ++p) P(" %s %.4f(%.4f)", mN.names[p].c_str(), nomG.th.size() == (size_t)mN.np ? nomG.th[p] : 0., (nomG.C.size() == (size_t)mN.np && nomG.C[p][p] > 0) ? std::sqrt(nomG.C[p][p]) : 0.);
    P("\n");
    Mat Cmc = bla::Zero(mN.np, mN.np);
    {
        Vec mean(mN.np, 0.);
        for (auto& t : thMC) for (int p = 0; p < mN.np; ++p) mean[p] += t[p] / thMC.size();
        for (auto& t : thMC) for (int p = 0; p < mN.np; ++p) for (int q = 0; q < mN.np; ++q) Cmc[p][q] += (t[p] - mean[p]) * (t[q] - mean[q]) / (thMC.size() - 1);
    }
    Mat Ccb = bla::Zero(mN.np, mN.np);
    StageA nomNoCB;
    {
        Model m = mN; m.useCBvar = false;
        nomNoCB = FitStageA(HA, m, crA, egA);
        std::vector<Vec> th;
        std::vector<TH2*> R(kNH);
        for (int h = 0; h < kNH; ++h) R[h] = (TH2*)allio::Data(f, "S0", "all", "CBrep", H[h].name);
        for (int k = 0; k < kNRepCB; ++k) {
            std::vector<HIn> Hr = HA;
            for (size_t h = 0; h < Hr.size(); ++h) {
                if (Hr[h].eg.empty() || Hr[h].eg[0] < 0) continue;
                const int hf = Hr[h].eg[0] / kNB;
                for (int i = 0; i < kNB; ++i) Hr[h].CB[i] = R[hf]->GetBinContent(i + 1, k + 1);
            }
            th.push_back(FitStageA(Hr, m, crA, egA, 0, 0, false, &nom.r.th).r.th);
        }
        for (auto r : R) delete r;
        Vec mean(mN.np, 0.);
        for (auto& t : th) for (int p = 0; p < mN.np; ++p) mean[p] += t[p] / th.size();
        for (auto& t : th) for (int p = 0; p < mN.np; ++p) for (int q = 0; q < mN.np; ++q) Ccb[p][q] += (t[p] - mean[p]) * (t[q] - mean[q]) / (th.size() - 1);
    }
    Mat Cstat = nom.r.Cs;
    for (int p = 0; p < mN.np; ++p) for (int q = 0; q < mN.np; ++q) Cstat[p][q] *= errScale * errScale;
    Mat Ctot = Cstat;
    for (int p = 0; p < mN.np; ++p) for (int q = 0; q < mN.np; ++q) Ctot[p][q] += Cmc[p][q];
    auto sdv = [](const Mat& C, int p) { return C[p][p] > 0 ? std::sqrt(C[p][p]) : 0.; };
    P("\n5. Adopted model %s, %s: parameters and errors\n", mN.Tag().c_str(), crA.tag.c_str());
    P("   %-8s %10s | %9s %9s %7s | %9s | %9s | %9s %9s | %9s\n", "param", "value", "stat(mod)", "stat(sw)", "sw/mod", "x scale", "MC", "CBrep", "CB(BB)", "total");
    for (int p = 0; p < mN.np; ++p) {
        const double cbBB = std::sqrt(std::max(0., nom.r.C[p][p] - nomNoCB.r.C[p][p]));
        P("   %-8s %10.5f | %9.5f %9.5f %7.3f | %9.5f | %9.5f | %9.5f %9.5f | %9.5f\n", mN.names[p].c_str(), nom.r.th[p], nom.r.ErrM(p), nom.r.Err(p),
          nom.r.ErrM(p) > 0 ? nom.r.Err(p) / nom.r.ErrM(p) : 0., sdv(Cstat, p), sdv(Cmc, p), sdv(Ccb, p), cbBB, sdv(Ctot, p));
    }
    P("   stat(mod): model covariance (CB and eta-template variances in the weights); stat(sw): sandwich with the event-group covariance of N;\n"
      "   x scale: stat(sw) x %.3f (chi2 fall-back); MC: 100 simulation replicas (smoothing redone); CBrep: 200 CB replicas, fit without the CB variance\n"
      "   (check of CB(BB) = the CB part of stat(mod)); total = stat(sw) x scale (+) MC\n", errScale);
    P("   correlation matrix (total):\n");
    for (int p = 0; p < mN.np; ++p) {
        P("   %-8s", mN.names[p].c_str());
        for (int q = 0; q < mN.np; ++q) P(" %6.3f", (Ctot[p][p] > 0 && Ctot[q][q] > 0) ? Ctot[p][q] / std::sqrt(Ctot[p][p] * Ctot[q][q]) : 0.);
        P("\n");
    }
    auto isBkgPar = [&](int p) { return p < 5 || (mN.tilt && p >= mN.ib && p < mN.ib + (mN.perSlice ? 5 : 1)) || (mN.quad && p == mN.iq); };

    // -------------------------------------------------------------------------------------------
    // 6. a_eta dependence and injection (A3)
    // -------------------------------------------------------------------------------------------
    std::map<std::string, StageA> V;
    for (double a : {0.0, 1.0, 2.0}) { Model m = mN; m.aEta = a; V[Form("aeta%.1f", a)] = FitStageA(HA, m, crA, egA); }
    Vec dth(mN.np, 0.);
    for (int p = 0; p < mN.np; ++p) dth[p] = V["aeta1.0"].r.th[p] - V["aeta0.0"].r.th[p];
    P("\n6. a_eta(0) dependence and injection test A3 (shift / total error of the adopted fit)\n   %-24s |", "case");
    for (int p = 0; p < mN.np; ++p) P(" %8s", mN.names[p].c_str());
    P(" | chi2/ndf | max (c_j, b, q)\n");
    double maxInjSig = 0, maxInjAll = 0;
    auto shiftLine = [&](const std::string& tag, const StageA& A, double* mx) {
        P("   %-24s |", tag.c_str());
        double m = 0;
        for (int p = 0; p < mN.np; ++p) {
            const double e = sdv(Ctot, p);
            const double z = e > 0 ? (A.r.th[p] - nom.r.th[p]) / e : 0.;
            P(" %8.3f", z);
            if (isBkgPar(p)) m = std::max(m, std::fabs(z));
        }
        P(" | %8.3f | %.3f\n", A.r.chi2P / A.r.ndf, m);
        if (mx) *mx = std::max(*mx, m);
    };
    for (double a : {0.0, 1.0, 2.0}) shiftLine(Form("fit at a_eta = %.1f", a), V[Form("aeta%.1f", a)], nullptr);
    for (double xx : {0.0, 1.0, 2.0}) {
        StageA As = FitStageA(HA, mN, crA, egA, xx - kAeta0, 0.);
        StageA Aa = FitStageA(HA, mN, crA, egA, xx - kAeta0, xx - kAeta0);
        shiftLine(Form("inject C2-C5 x %.1f", xx), As, &maxInjSig);
        shiftLine(Form("inject all eta x %.1f", xx), Aa, &maxInjAll);
    }
    P("   A3 (signal classes C2-C5 injected at 0, 1, 2 x SMASH, fit at a_eta(0) = %.2f): max shift %.3f sigma -> %s\n", kAeta0, maxInjSig, maxInjSig < 0.3 ? "PASS" : "FAIL");
    P("      with G_eta injected as well: max shift %.3f sigma (the c_j depend on the assumed a_eta through G_eta, D10: rerun at the WP5 value)\n", maxInjAll);

    // -------------------------------------------------------------------------------------------
    // 7. variants
    // -------------------------------------------------------------------------------------------
    {
        Model mg = mN; mg.cbGM = true; V["CBgm"] = FitStageA(HA, mg, crA, egA);
        Model mk = mN; mk.kBsig = 0.2; V["kBcons"] = FitStageA(HA, mk, crA, egA);
        V["pre"] = FitStageA(Hpre, mN, crA, egPre);
        V["ghostA"] = FitStageA(Hgh, mN, crA, egS0);
        Model me = mN; me.eOnly = -1; me.kBsig = 0; V["E1E2"] = FitStageA(H, me, crA, egS0);
        V["CRU565"] = FitStageA(H, mN, crU565, egS0);
        V["specCR"] = FitStageA(H, mN, crN, egS0);
        V["fb1"] = FitStageA(H44, mN, crFB1, nullptr);
        Model mt = mN; mt.perSlice = false; mt.quad = false; mt.tilt = true; mt.omega = false; V["M1"] = FitStageA(HA, mt, crA, egA);
        V["D11data"] = chain[adoptCfg].n.M[chain[adoptCfg].n.adoptData];
        V["D11dataOmega"] = nomOmega.m.perSlice ? nomOmega : V["D11data"];
        { Model mo = chain[adoptCfg].n.M[chain[adoptCfg].n.adoptData].m; mo.omega = true; V["D11dataOmega"] = FitStageA(HA, mo, crA, egA); }
        // K0s veto (A5): S5 without M(pipi) 0.485-0.507, and the same S5 rebuilt from G3 without veto
        std::vector<HIn> HK0 = H, HG3 = H;
        for (int e = 0; e < 2; ++e) {
            HK0[8 + e] = FromG3(f, "S0", Form("S5vetoE%d", e + 1), 4, e, 0.43, 2.0, 0.485, 0.507, H[8 + e], true, true); SmoothAll(HK0[8 + e], mN.omega);
            HG3[8 + e] = FromG3(f, "S0", Form("S5g3E%d", e + 1), 4, e, 0.43, 2.0, 0, 0, H[8 + e], true, true); SmoothAll(HG3[8 + e], mN.omega);
        }
        V["K0veto"] = FitStageA(HK0, mN, crA, nullptr);
        V["G3ref"] = FitStageA(HG3, mN, crA, nullptr);
        // CR-U split at M(pipi) 0.42 (A4): S4 -> [0.39, 0.42) and [0.42, 0.43) (both c_4)
    }
    P("\n7. Variants (adopted model and CR unless stated): chi2/ndf and parameters (stat errors, sandwich)\n");
    auto varLine = [&](const std::string& tag, const StageA& A) {
        if ((int)A.r.th.size() != A.m.np) { P("   %-46s | FIT FAILED\n", tag.c_str()); return; }
        P("   %-46s | %8.1f %4d %6.3f |", tag.c_str(), A.r.chi2P, A.r.ndf, A.r.chi2P / A.r.ndf);
        for (int p = 0; p < A.m.np; ++p) P(" %s=%.4f(%.4f)", A.m.names[p].c_str(), A.r.th[p], A.r.Err(p));
        P("%s\n", A.r.ok ? "" : " [NOT CONVERGED]");
    };
    varLine("adopted", nom);
    varLine("adopted, CB variance not in the weights", nomNoCB);
    if (!omegaAdopted) varLine("+ omega split (c_om)", nomOmega);
    for (double a : {0.0, 1.0, 2.0}) varLine(Form("a_eta(0) = %.1f", a), V[Form("aeta%.1f", a)]);
    varLine("common tilt only (M1)", V["M1"]);
    varLine("D11 with the data-only chi2 (literal)", V["D11data"]);
    varLine("D11 data-only + omega split (literal chain)", V["D11dataOmega"]);
    varLine("kappa_B = 1 +- 0.2 (constrained)", V["kBcons"]);
    varLine("CB per-cell geometric mean (S1)", V["CBgm"]);
    varLine("ghost variant (a) templates", V["ghostA"]);
    varLine("E1 + E2, k_E1 free, kappa_B free", V["E1E2"]);
    varLine("CR-U 0.565-0.700", V["CRU565"]);
    varLine("spec CRs (CR-U 0.565 + CR-K >= 0.43)", V["specCR"]);
    varLine("fall-back 1 (CR-U 0.575 + CR-K >= 0.44)", V["fb1"]);
    varLine("S5 rebuilt from G3, no veto (A5 reference)", V["G3ref"]);
    varLine("K0s veto 0.485-0.507 in S5 (A5)", V["K0veto"]);
    varLine("pre (parallel column)", V["pre"]);

    // -------------------------------------------------------------------------------------------
    // 8. acceptance tests A4, A5, chi2; closure in CR-K and CR-L
    // -------------------------------------------------------------------------------------------
    P("\n8. Acceptance tests on the adopted configuration\n");
    const double chi2ndf = nom.r.chi2P / nom.r.ndf;
    P("   chi2/ndf <= 1.5, with the MC template covariance (GLS): %.1f/%d = %.3f -> %s; data only: Pearson %.1f/%d = %.3f (Baker-Cousins %.1f, event-group-scaled %.1f) -> %s\n",
      nomG.chi2, nomG.ndf, nomG.chi2 / nomG.ndf, nomG.chi2 / nomG.ndf <= 1.5 ? "PASS" : "FAIL", nom.r.chi2P, nom.r.ndf, chi2ndf, nom.r.chi2BC, nom.r.chi2E, chi2ndf <= 1.5 ? "PASS" : "FAIL");
    {
        double sMC = 0, sSt = 0;
        for (int p = 0; p < 5; ++p) { sMC += std::sqrt(Cmc[mN.ic + p][mN.ic + p]); sSt += std::sqrt(Cstat[mN.ic + p][mN.ic + p]); }
        P("   MC error vs data statistical error (mean over c_j): %.4f vs %.4f -> fall-back 'MC error > data statistical error: simulation-limited' %s\n", sMC / 5, sSt / 5,
          sMC > sSt ? "FIRES" : "not fired");
        gSimLimited = sMC > sSt;
    }
    const NestOut& nA = chain[adoptCfg].n;
    double a4p = 1, a4dc = 0;
    {
        const int i1 = nA.adopt == 2 || nA.adopt == 4 ? 2 : 1, i3 = i1 == 2 ? 4 : 3;
        a4dc = nA.G[i1].chi2 - nA.G[i3].chi2;
        a4p = TMath::Prob(std::max(0., a4dc), 4);
        const double a4dd = nA.M[i1].r.chi2P - nA.M[i3].r.chi2P;
        P("   A4 per-slice tilts vs a common b (%s vs %s): dchi2 (MC-inclusive) = %.2f for 4 dof, p = %.3g -> %s (data-only dchi2 %.1f, p = %.2e)\n",
          nA.M[i1].m.Tag().c_str(), nA.M[i3].m.Tag().c_str(), a4dc, a4p, a4p > 0.01 ? "PASS" : "FAIL", a4dd, TMath::Prob(std::max(0., a4dd), 4));
        const StageA& A3 = nA.M[i3];
        P("      per-slice tilts:");
        for (int j = 0; j < 5; ++j) P(" b_%d = %.3f +- %.3f", j + 1, A3.r.th[A3.m.ib + j], A3.r.Err(A3.m.ib + j));
        P("  (common b = %.3f +- %.3f)\n", nA.M[i1].r.th[nA.M[i1].m.ib], nA.M[i1].r.Err(nA.M[i1].m.ib));
    }
    double splitR[2] = {0, 0}, splitE[2] = {0, 0}, splitZ = 0, splitS4Z = 0;
    {
        // data / model in CR-U for M(pipi) < 0.42 and >= 0.42; S4 halves from G3 with raw templates; model = adopted fit
        std::vector<HIn> H42;
        for (int h = 0; h < kNH; ++h) if (H[h].j != 3) H42.push_back(H[h]);
        for (int e = 0; e < 2; ++e) {
            H42.push_back(FromG3(f, "S0", Form("S4[0.39,0.42)E%d", e + 1), 3, e, 0.39, 0.42, 0, 0, H[6 + e], true, false));
            H42.push_back(FromG3(f, "S0", Form("S4[0.42,0.43)E%d", e + 1), 3, e, 0.42, 0.43, 0, 0, H[6 + e], true, false));
        }
        for (auto& x : H42) { x.D = x.s[lGD]; x.C = x.s[lGC]; x.Dom = x.s[lGDom]; x.Com = x.s[lGCom]; }   // raw templates for the split sums
        double n[2] = {0, 0}, m[2] = {0, 0}, vm[2] = {0, 0}, nS4[2] = {0, 0}, mS4[2] = {0, 0}, vS4[2] = {0, 0};
        for (const auto& x : H42) {
            if (!Used(x, mN)) continue;
            const bool hiPart = (x.j == 4) || (x.name.find("[0.42") != std::string::npos);
            const bool isS4 = x.name.find("S4[") != std::string::npos;
            for (int i = 0; i < kNB; ++i) {
                if (!(x.lo[i] >= crA.uLo - 1e-9 && x.hi[i] <= crA.uHi + 1e-9)) continue;
                const double kE = nom.r.th[x.e == 0 ? mN.ikE1 : mN.ikE2];
                const double T = 1 + (mN.tilt ? nom.r.th[mN.ib + (mN.perSlice ? x.j : 0)] * x.x[i] : 0.) + (mN.quad ? nom.r.th[mN.iq] * x.x[i] * x.x[i] : 0.);
                const double c = nom.r.th[mN.ic + x.j], kB = nom.r.th[mN.ikB];
                const double B = Bkg(x, i, mN, nom.r.th, nullptr);   // RAW templates (set below), omega parts included
                const double vB = std::pow(c * T, 2) * (x.v[lGD][i] + kB * kB * x.v[lGC][i]);
                const double mu = kE * x.CB[i] + B + EtaSum(x, i, mN);
                const int k = hiPart ? 1 : 0;
                n[k] += x.N[i]; m[k] += mu; vm[k] += vB + x.vCB[i];
                if (isS4) { nS4[k] += x.N[i]; mS4[k] += mu; vS4[k] += vB + x.vCB[i]; }
            }
        }
        for (int k = 0; k < 2; ++k) { splitR[k] = n[k] / m[k]; splitE[k] = splitR[k] * std::sqrt(1.0 / n[k] + vm[k] / (m[k] * m[k])); }
        splitZ = (splitR[0] - splitR[1]) / std::hypot(splitE[0], splitE[1]);
        P("   A4 CR-U split at M(pipi) 0.42: data/model %.4f +- %.4f (< 0.42, N %.0f) vs %.4f +- %.4f (>= 0.42, N %.0f): %.2f sigma -> %s\n", splitR[0], splitE[0], n[0],
          splitR[1], splitE[1], n[1], splitZ, std::fabs(splitZ) < 2 ? "PASS" : "FAIL");
        const double r0 = nS4[0] / mS4[0], r1 = nS4[1] / mS4[1];
        const double e0 = r0 * std::sqrt(1 / nS4[0] + vS4[0] / (mS4[0] * mS4[0])), e1 = r1 * std::sqrt(1 / nS4[1] + vS4[1] / (mS4[1] * mS4[1]));
        splitS4Z = (r0 - r1) / std::hypot(e0, e1);
        P("      inside S4: %.4f +- %.4f (0.39-0.42) vs %.4f +- %.4f (0.42-0.43): %.2f sigma\n", r0, e0, r1, e1, splitS4Z);
    }
    double a5max = 0;
    {
        const StageA& K = V["K0veto"]; const StageA& G = V["G3ref"];
        P("   A5 K0s veto (M(pipi) 0.485-0.507 removed from S5 in data and simulation; reference: S5 rebuilt from G3 without veto): shift / total error\n     ");
        for (int p = 0; p < mN.np; ++p) { const double e = sdv(Ctot, p); const double z = e > 0 ? (K.r.th[p] - G.r.th[p]) / e : 0.; if (e > 0) a5max = std::max(a5max, std::fabs(z)); P(" %s %.2f", mN.names[p].c_str(), z); }
        P(" -> max %.2f sigma -> %s\n", a5max, a5max < 1 ? "PASS" : "FAIL");
        double mm = 0; for (int p = 0; p < mN.np; ++p) { const double e = sdv(Ctot, p); if (e > 0) mm = std::max(mm, std::fabs(G.r.th[p] - nom.r.th[p]) / e); }
        P("      (G3 rebuild without veto vs the F10 fit: max shift %.3f sigma; it lacks only the event-group covariance)\n", mm);
    }
    // closure: CR-K (E-classes of the model) and CR-L, prediction from the adopted fit
    Pred pK = Predict(HA, nom, [](const HIn& h, int i) { return h.crkOK && h.hi[i] <= 0.5225 + 1e-9; }, &Ctot);
    Pred pL = Predict(HA, nom, [](const HIn& h, int i) { return h.lo[i] >= 0.5225 - 1e-9 && h.hi[i] <= 0.5300 + 1e-9; }, &Ctot);
    Mat CtotN = Ctot;
    auto predLine = [&](const char* tag, const Pred& p) {
        const double sp = std::sqrt(p.vfit + p.vcb);
        P("   closure %s: N = %.0f, predicted %.1f (B_A %.1f + k CB %.1f + eta(0.4) %.1f) +- %.1f (stage A) +- %.1f (CB): pull %.2f; S10 term max(|pred - obs|, sigma_pred) = %.1f rows (%.2f %% of N - k CB)\n",
          tag, p.N, p.mu, p.B, p.cb, p.eta, std::sqrt(p.vfit), std::sqrt(p.vcb), (p.N - p.mu) / std::sqrt(p.N + p.vfit + p.vcb), std::max(std::fabs(p.mu - p.N), sp),
          100 * std::max(std::fabs(p.mu - p.N), sp) / (p.N - p.cb));
    };
    predLine("CR-K (S5, 0.400-0.5225)", pK);
    predLine("CR-L (0.5225-0.5300)", pL);

    // -------------------------------------------------------------------------------------------
    // 9. Table C3 (counted "after" display) and per-bin background
    // -------------------------------------------------------------------------------------------
    std::vector<Vec> binB(H.size(), Vec(kNB, 0.)), binEs(H.size(), Vec(kNB, 0.)), binEmc(H.size(), Vec(kNB, 0.));
    for (int h = 0; h < kNH; ++h) {
        std::vector<double> mc;
        for (int i = 0; i < kNB; ++i) {
            Vec g(mN.np, 0.);
            binB[h][i] = Bkg(H[h], i, mN, nom.r.th, &g);
            binEs[h][i] = std::sqrt(QF(g, Cstat));
            double m = 0, m2 = 0;
            for (int k = 0; k < (int)HR.size(); ++k) {
                // replica k: its templates and its parameters (F10 histogram h of the replica set; HA = H unless fall-back 1 is adopted)
                const HIn& xr = HR[k][std::min<size_t>(h, HR[k].size() - 1)];
                const double v = (xr.name == H[h].name) ? Bkg(xr, i, mN, thMC[k], nullptr) : binB[h][i];
                m += v; m2 += v * v;
            }
            m /= HR.size(); m2 /= HR.size();
            binEmc[h][i] = std::sqrt(std::max(0., m2 - m * m));
        }
    }
    struct Reg { const char* tag; double lo, hi; int win; bool crk; };
    const std::vector<Reg> regs = {{"R0' 0.400-0.455", 0.400, 0.455, -1, false}, {"R1 = B 0.4550-0.5225", 0.455, 0.5225, kWB, false},
                                   {"R2 0.5225-0.5336", 0.5225, 0.53364, kWR2, false}, {"R3 = P 0.5336-0.5542", 0.53364, 0.55416, kWP, false},
                                   {"R4 0.5542-0.565", 0.55416, 0.565, kWR4, false}, {"R5 0.565-0.700", 0.565, 0.700, kWCRU, false},
                                   {"CR-U 0.575-0.700", 0.575, 0.700, -1, false}, {"CR-L 0.5225-0.5300", 0.5225, 0.5300, kWCRL, false},
                                   {"CR-K S5 0.400-0.5225", 0.400, 0.5225, kWCRK, true}};
    struct RegOut { double N = 0, CB = 0, eCB = 0, B = 0, eBs = 0, eBmc = 0, eta1 = 0, eta04 = 0, Geta1 = 0, sig1 = 0; };
    auto regCalc = [&](const Reg& R, int esel) {
        RegOut o;
        Vec g(mN.np, 0.);
        std::vector<double> mcs(HR.size(), 0.);
        Model m1 = mN; m1.aEta = 1.0;
        for (int h = 0; h < kNH; ++h) {
            const HIn& x = H[h];
            if (esel >= 0 && x.e != esel) continue;
            if (R.crk && x.j != 4) continue;
            const double kE = nom.r.th[x.e == 0 ? mN.ikE1 : mN.ikE2];
            for (int i = 0; i < kNB; ++i) {
                const double fr2 = Frac(x.lo[i], x.hi[i], R.lo, R.hi);
                if (fr2 <= 0) continue;
                Vec gg(mN.np, 0.);
                const double b = Bkg(x, i, mN, nom.r.th, &gg);
                for (int k = 0; k < mN.np; ++k) g[k] += fr2 * gg[k];
                o.B += fr2 * b; o.N += fr2 * x.N[i]; o.CB += fr2 * kE * x.CB[i]; o.eCB += fr2 * fr2 * kE * kE * x.vCB[i];
                o.eta1 += fr2 * EtaSum(x, i, m1); o.eta04 += fr2 * EtaSum(x, i, mN); o.Geta1 += fr2 * x.s[lGeta][i];
                o.sig1 += fr2 * (x.s[lC2][i] + x.s[lC3][i] + x.s[lC4][i] + x.s[lC5][i]);
                for (size_t k = 0; k < HR.size(); ++k) {
                    const HIn& xr = HR[k][std::min<size_t>(h, HR[k].size() - 1)];
                    mcs[k] += fr2 * ((xr.name == x.name) ? Bkg(xr, i, mN, thMC[k], nullptr) : b);
                }
            }
        }
        // all M(ee): exact window counts of N (fractional bins otherwise)
        if (esel < 0 && R.win >= 0) { TH1* w = allio::Data(f, "S0", "all", "N", "win"); o.N = allio::Win(w, R.win); delete w; }
        o.eCB = std::sqrt(o.eCB);
        o.eBs = std::sqrt(QF(g, Cstat));
        double m = 0, m2 = 0; for (double v : mcs) { m += v; m2 += v * v; } m /= mcs.size(); m2 /= mcs.size();
        o.eBmc = std::sqrt(std::max(0., m2 - m * m));
        return o;
    };
    std::vector<RegOut> roE2, roAll;
    for (const auto& R : regs) { roE2.push_back(regCalc(R, 1)); roAll.push_back(regCalc(R, -1)); }
    auto c3Print = [&](const char* title, const std::vector<RegOut>& ro) {
        P("%s\n", title);
        P("   %-22s | %8s %8s %6s | %9s %7s %7s | %8s %8s %8s | %9s %7s | %7s\n", "region", "N", "k CB", "err", "B_A", "stat", "MC", "eta(1)", "eta(0.4)", "Geta(1)",
          "N-kCB-B_A", "err", "ratio");
        for (size_t q = 0; q < regs.size(); ++q) {
            const RegOut& o = ro[q];
            const double ex = o.N - o.CB - o.B, eex = std::sqrt(o.N + o.eCB * o.eCB + o.eBs * o.eBs + o.eBmc * o.eBmc);
            P("   %-22s | %8.0f %8.1f %6.1f | %9.1f %7.1f %7.1f | %8.1f %8.1f %8.1f | %9.1f %7.1f | %7.4f\n", regs[q].tag, o.N, o.CB, o.eCB, o.B, o.eBs, o.eBmc, o.eta1,
              o.eta04, o.Geta1, ex, eex, (o.N - o.CB) / (o.B + o.eta04));
        }
        P("   ratio = (N - k CB) / (B_A + 0.4 eta). DISPLAY COUNTS, NOT YIELDS (D20).\n");
    };
    c3Print("\n9. Table C3, M(ee) >= 10 MeV (E2, fitted), raw rows", roE2);
    c3Print("   Table C3, all M(ee) (information: E1 predicted with k_E1 = 1 and the E2 parameters), raw rows", roAll);

    // -------------------------------------------------------------------------------------------
    // 10. Table C1: candidate control regions
    // -------------------------------------------------------------------------------------------
    P("\n10. Table C1: candidate control regions (S0, isBest >= 0; simulation at SMASH x PDG BR, mode A, w x V3)\n");
    struct C1Row { std::string reg, n, cbn, sig, geta, ratio, verdict; };
    std::vector<C1Row> c1;
    {
        // helper over F10 cells
        auto cells = [&](std::function<bool(const HIn&, int)> sel, double& N, double& CB, double& T12, double& sig, double& geta) {
            N = CB = T12 = sig = geta = 0;
            for (const auto& h : H) for (int i = 0; i < kNB; ++i) {
                if (!sel(h, i)) continue;
                N += h.N[i]; CB += h.CB[i]; T12 += h.s[lT12][i]; geta += h.s[lGeta][i];
                sig += h.s[lC2][i] + h.s[lC3][i] + h.s[lC4][i] + h.s[lC5][i];
            }
        };
        auto add = [&](const std::string& reg, std::function<bool(const HIn&, int)> sel, const std::string& verdict) {
            double N, CB, T, s, g; cells(sel, N, CB, T, s, g);
            c1.push_back({reg, Form("%.0f", N), Form("%.3f", CB / N), Form("%.2f", 100 * s / T), Form("%.1f", 100 * g / T), Form("%.3f", (N - CB) / T), verdict});
        };
        add("CR-U 0.565-0.700, all slices, all M(ee)", [](const HIn& h, int i) { return h.lo[i] >= 0.565 - 1e-9; },
            "spec; E1 part dropped (k_E1 fall-back), first 10 MeV dropped (A1 fall-back)");
        add("CR-U 0.565-0.700, E1 (M(ee) < 10 MeV)", [](const HIn& h, int i) { return h.e == 0 && h.lo[i] >= 0.565 - 1e-9; },
            Form("rejected (fall-back): k_E1 = %.1f, %.1f with ghost variant (a); outside 0.8-3", kE1v, sGh.r.th[sGh.m.ikE1]));
        add("CR-U 0.575-0.700, E2 (adopted)", [](const HIn& h, int i) { return h.e == 1 && h.lo[i] >= 0.575 - 1e-9; }, "ADOPTED (normalisation and tilt)");
        add("CR-K S5 (M(pipi) >= 0.43), 0.400-0.5225", [](const HIn& h, int i) { return h.j == 4 && h.hi[i] <= 0.5225 + 1e-9; },
            "rejected by the A1/A2 fall-back (C5 radiative tail, S5 threshold, CB N++ = 2)");
        add("CR-L 0.5225-0.5300", [](const HIn& h, int i) { return h.lo[i] >= 0.5225 - 1e-9 && h.hi[i] <= 0.5300 + 1e-9; }, "test only (closure term S10)");
        add("below the bump 0.400-0.455", [](const HIn& h, int i) { return h.hi[i] <= 0.455 + 1e-9; }, "rejected: bump tail (C2/C3 at q > 0)");
        // 0.60-0.80 from fine_M4 (no slices)
        TH1* Nf = allio::Data(f, "S0", "all", "N", "fine_M4"); TH1* CBf = allio::Data(f, "S0", "all", "CB", "fine_M4");
        TH1* Tf = allio::Sim(f, "S0", "nom", kT12, "fine_M4"); TH1* Gf = allio::Sim(f, "S0", "nom", kGeta, "fine_M4");
        TH1* Sf[4] = {allio::Sim(f, "S0", "nom", kC2, "fine_M4"), allio::Sim(f, "S0", "nom", kC3, "fine_M4"), allio::Sim(f, "S0", "nom", kC4, "fine_M4"), allio::Sim(f, "S0", "nom", kC5, "fine_M4")};
        auto fineRow = [&](const char* reg, double a, double b, const std::string& verdict) {
            double N = 0, CB = 0, T = 0, G = 0, S = 0;
            for (int i = 1; i <= Nf->GetNbinsX(); ++i) {
                const double c = Nf->GetXaxis()->GetBinCenter(i);
                if (c < a || c >= b) continue;
                N += Nf->GetBinContent(i); CB += CBf->GetBinContent(i); T += Tf->GetBinContent(i) * kNbToRaw; G += Gf->GetBinContent(i) * kNbToRaw;
                for (int k = 0; k < 4; ++k) S += Sf[k]->GetBinContent(i) * kNbToRaw;
            }
            c1.push_back({reg, Form("%.0f", N), Form("%.3f", CB / N), Form("%.2f", 100 * S / T), Form("%.1f", 100 * G / T), Form("%.3f", (N - CB) / T), verdict});
        };
        fineRow("0.60-0.65 (part of 0.60-0.80)", 0.60, 0.65, "");
        fineRow("0.65-0.70", 0.65, 0.70, "");
        fineRow("0.70-0.75", 0.70, 0.75, "");
        fineRow("0.75-0.80", 0.75, 0.80, "rejected 0.60-0.80: data/sim drifts, long extrapolation");
        // vz < -160 mm: none - pre (data only; the simulation has 0.15 % of its rows there)
        TH1* wn = allio::Data(f, "none", "all", "N", "win"); TH1* wp = allio::Data(f, "pre", "all", "N", "win");
        TH1* cn = allio::Data(f, "none", "all", "CB", "win"); TH1* cp = allio::Data(f, "pre", "all", "CB", "win");
        const double nU = allio::Win(wn, kWCRU) - allio::Win(wp, kWCRU), cU = allio::Win(cn, kWCRU) - allio::Win(cp, kWCRU);
        const double nP = allio::Win(wn, kWP) - allio::Win(wp, kWP), cP = allio::Win(cn, kWP) - allio::Win(cp, kWP);
        c1.push_back({"vz < -160 mm (none - pre), CR-U range", Form("%.0f", nU), Form("%.3f", cU / nU), "n/a", "n/a", "no sim",
                      Form("rejected: no simulation (0.15 %% of sim rows), CB/N %.2f in P; no eta (A1 notes: S = -15 +- 118)", cP / nP)});
        // like-sign
        double npp = 0, nmm = 0, nU2 = 0;
        for (const auto& h : H) for (int i = 0; i < kNB; ++i) if (h.lo[i] >= 0.565 - 1e-9) { npp += h.Npp[i]; nmm += h.Nmm[i]; nU2 += h.N[i]; }
        c1.push_back({"like-sign e+e+, e-e- (CR-U range)", Form("%.0f + %.0f", npp, nmm), "--", "0", "0", "--", "not a background CR: it is the CB estimate (k_E2 constrained)"});
        const std::vector<double>& tg = RS["TAG S0 2"];
        const std::vector<double>& tk = RS["TAG S0 0"];
        if (tg.size() >= 4 && tk.size() >= 3)
            c1.push_back({"ECAL-tagged eta (all M4)", Form("%.0f", tg[0]), "--", Form("T = %.0f +- %.0f", tg[1], tg[2]), "--", Form("S/B %.2f", tg[3]),
                          Form("signal-enriched: cross-check of a_eta only (D14); in CR-K T = %.1f +- %.1f", tk[1], tk[2])});
        delete Nf; delete CBf; delete Tf; delete Gf; for (auto x : Sf) delete x; delete wn; delete wp; delete cn; delete cp;
    }
    P("   %-44s | %10s | %6s | %8s | %6s | %8s | %s\n", "region", "N (rows)", "CB/N", "C2-C5 %", "Geta %", "(N-CB)/T12", "verdict");
    for (const auto& r : c1) P("   %-44s | %10s | %6s | %8s | %6s | %8s | %s\n", r.reg.c_str(), r.n.c_str(), r.cbn.c_str(), r.sig.c_str(), r.geta.c_str(), r.ratio.c_str(), r.verdict.c_str());

    // -------------------------------------------------------------------------------------------
    // 11. frozen output: stageA.txt, stageA.root
    // -------------------------------------------------------------------------------------------
    {
        FILE* o = fopen(Form("%s/stageA%s.txt", kCntOut, sfx.c_str()), "w");
        fprintf(o, "# =====================================================================================================\n");
        fprintf(o, "# WP3 stage A (SPEC report_all, D9-D11): FROZEN non-eta background normalisation from the control regions\n");
        fprintf(o, "# written by plots/all_report/wp3_stageA.C, 2026-11-14; inputs hists_v1.root (WP0), wp3_rows.root (wp3_rows.C)\n");
        fprintf(o, "# Fall-backs that fired (SPEC Sect. 6): %s%s%s%s\n", fbKE1 ? "[k_E1 outside 0.8-3 -> E1 dropped (kappa = 1 +- 0.2 applies to stage B)] " : "",
                adoptCfg >= 1 ? "[A1/A2 -> CR-U from 0.575, CR-K at 0.44] " : "", adoptCfg == 2 ? "[still failing -> CR-K dropped, tilt from CR-U only, bump EXTRAPOLATED, CR-L term in S10] " : "",
                errScale > 1 ? Form("[chi2/ndf > 1.5 -> %sstat errors x %.3f, background-model limited]", omegaAdopted ? "omega split, " : "", errScale) : "");
        fprintf(o, "# =====================================================================================================\n");
        fprintf(o, "# Model, bin i of F10 histogram (j = M(pipi) slice 1..5, e = M(ee) class 1 = E1 < 10 MeV, 2 = E2 >= 10 MeV), RAW data rows:\n");
        fprintf(o, "#   mu = k_e CB + c_j (1 + b_j x + q x^2) [G_D + kappa_B G_C]%s + a_eta [C2 + kappa (C3 + C4) + rho_d C5 + g G_eta]\n", mN.omega ? " (omega parts x c_om)" : "");
        fprintf(o, "#   (common b if model M1/M2; b_j per slice if M3/M4; q only for M2/M4)  x = (M4 - 0.55)/0.1 at the F10 bin centre\n");
        fprintf(o, "#   G_D, G_C = SMOOTHED templates (exp(pol3) per histogram, parameters below); sim in expected raw rows = nb x %.4f\n", kNbToRaw);
        fprintf(o, "#   (L = %.1f nb^-1, x 63/64); eta part FIXED at a_eta = %.2f, kappa = rho_d = g = 1; CB = 2 sqrt(N++ N--), spread rule\n", kL, kAeta0);
        fprintf(o, "#   (hists_v1 CB_F10_*); k_E2 = 1 +- 0.05 constrained; k_E1 %s.\n", mN.eOnly == 1 ? "not used (E1 dropped; fixed to 1)" : "free");
        fprintf(o, "# Errors: stat = sandwich with the event-group covariance of N (CB replica variance and eta-template variance in the weights)\n");
        fprintf(o, "#         x error_scale; mc = 100 simulation replicas (smoothing redone); total = stat + mc (covariances added).\n");
        fprintf(o, "# Model selection (D11) and the chi2 criterion use the chi2 WITH the MC covariance of the smoothed templates (GLS, 100 replicas):\n");
        fprintf(o, "# the templates are simulation-limited (MC error of c_j %.3f vs stat %.3f: fall-back 'simulation-limited' fired); with the data-only\n",
                std::sqrt(Cmc[mN.ic][mN.ic]), sdv(Cstat, mN.ic));
        fprintf(o, "# chi2 D11 would add a quadratic term and per-slice tilts (then the omega split) that follow MC fluctuations (variants D11data, D11dataOmega).\n");
        fprintf(o, "# Parameter values: Poisson-ML fit with the data variance (spec D12/D13); values_gls = the GLS estimates (cross-check).\n");
        fprintf(o, "simulation_limited %d\nmodel_selection D11_with_MC_covariance\n", gSimLimited ? 1 : 0);
        fprintf(o, "selection S0\nlayout F10\nmeec_classes_used %s\ncontrol_region %s\nmodel %s\n", mN.eOnly == 1 ? "E2" : "E1 E2", crA.tag.c_str(), mN.Tag().c_str());
        fprintf(o, "cr_u %.4f %.4f\ncr_k %s\n", crA.uLo, crA.uHi, crA.useK ? Form("%.4f %.4f S5", crA.kLo, crA.kHi) : "none");
        fprintf(o, "bump_background_label %s\n", crA.useK ? "interpolated (CR-K + CR-U)" : "EXTRAPOLATED (tilt from CR-U only)");
        fprintf(o, "a_eta0 %.3f\nchi2_with_mc_cov %.3f\nndf_with_mc_cov %d\nchi2_pearson_data_only %.3f\nchi2_bakercousins_data_only %.3f\nnbins %d\nndf %d\nerror_scale %.5f\n", kAeta0, nomG.chi2, nomG.ndf, nom.r.chi2P, nom.r.chi2BC, nom.r.n, nom.r.ndf, errScale);
        fprintf(o, "values_gls"); for (int p = 0; p < mN.np; ++p) fprintf(o, " %.8g", nomG.th.size() == (size_t)mN.np ? nomG.th[p] : 0.); fprintf(o, "\n");
        fprintf(o, "npar %d\nnames", mN.np); for (auto& s : mN.names) fprintf(o, " %s", s.c_str()); fprintf(o, "\n");
        fprintf(o, "fixed"); for (int p = 0; p < mN.np; ++p) fprintf(o, " %d", (p == mN.ikE1 && (mN.eOnly == 1 || !mN.kE1free)) ? 1 : 0); fprintf(o, "\n");
        fprintf(o, "constraints k_E2 1 %.3f", kKE2sig); if (mN.kBsig > 0) fprintf(o, " kappa_B 1 %.3f", mN.kBsig); fprintf(o, "\n");
        fprintf(o, "values"); for (int p = 0; p < mN.np; ++p) fprintf(o, " %.8g", nom.r.th[p]); fprintf(o, "\n");
        fprintf(o, "err_stat"); for (int p = 0; p < mN.np; ++p) fprintf(o, " %.6g", sdv(Cstat, p)); fprintf(o, "\n");
        fprintf(o, "err_mc"); for (int p = 0; p < mN.np; ++p) fprintf(o, " %.6g", sdv(Cmc, p)); fprintf(o, "\n");
        fprintf(o, "err_total"); for (int p = 0; p < mN.np; ++p) fprintf(o, " %.6g", sdv(Ctot, p)); fprintf(o, "\n");
        auto wm = [&](const char* tag, const Mat& C) { fprintf(o, "%s %d\n", tag, mN.np); for (int p = 0; p < mN.np; ++p) { for (int q = 0; q < mN.np; ++q) fprintf(o, " %.8e", C[p][q]); fprintf(o, "\n"); } };
        wm("cov_stat", Cstat); wm("cov_mc", Cmc); wm("cov_total", Ctot);
        fprintf(o, "# a_eta dependence (same model and CR): dtheta_daeta = theta(a_eta = 1) - theta(a_eta = 0); linear interpolation is good to\n");
        fprintf(o, "# the level shown by values_aeta2.0; exact: rerun wp3_stageA.C with kAeta0 = the WP5 value\n");
        fprintf(o, "dtheta_daeta"); for (int p = 0; p < mN.np; ++p) fprintf(o, " %.8g", dth[p]); fprintf(o, "\n");
        for (double a : {0.0, 1.0, 2.0}) { const StageA& A = V[Form("aeta%.1f", a)]; fprintf(o, "values_aeta%.1f", a); for (int p = 0; p < mN.np; ++p) fprintf(o, " %.8g", A.r.th[p]); fprintf(o, "\n"); }
        fprintf(o, "# smoothing per F10 histogram: name start_GD p0 p1 p2 p3 start_GC p0 p1 p2 p3; density [raw rows per GeV/c2] =\n");
        fprintf(o, "#   exp(p0 + p1 t + p2 t^2 + p3 t^3), t = (M4 - 0.55)/0.15, for M4 >= start (zero below); bin content = integral over the bin\n");
        for (int h = 0; h < kNH; ++h) fprintf(o, "smooth %s %.4f %.10g %.10g %.10g %.10g %.4f %.10g %.10g %.10g %.10g\n", H[h].name.c_str(), H[h].lo[H[h].domD], H[h].pD[0], H[h].pD[1],
                                           H[h].pD[2], H[h].pD[3], H[h].lo[H[h].domC], H[h].pC[0], H[h].pC[1], H[h].pC[2], H[h].pC[3]);
        fprintf(o, "# per-bin stage-A background B_A = c_j (...) [G_D + kappa_B G_C] in raw rows, all 10 F10 histograms (E1 rows = prediction\n");
        fprintf(o, "# with the E2-derived parameters if E1 is dropped); columns: hist bin M4lo M4hi B_A err_stat err_mc G_D_smooth G_C_smooth in_CR\n");
        for (int h = 0; h < kNH; ++h) for (int i = 0; i < kNB; ++i)
            fprintf(o, "bin %s %2d %.5f %.5f %.6g %.6g %.6g %.6g %.6g %d\n", H[h].name.c_str(), i + 1, H[h].lo[i], H[h].hi[i], binB[h][i], binEs[h][i], binEmc[h][i], H[h].D[i],
                    H[h].C[i], (Used(H[h], mN) && InCR(H[h], i, crA)) ? 1 : 0);
        fprintf(o, "# variants: name chi2P ndf name=value(stat error)...\n");
        auto vline = [&](const char* nm, const StageA& A) { if ((int)A.r.th.size() != A.m.np) { fprintf(o, "variant %s FAILED\n", nm); return; } fprintf(o, "variant %s %.3f %d", nm, A.r.chi2P, A.r.ndf); for (int p = 0; p < A.m.np; ++p) fprintf(o, " %s=%.6g(%.3g)", A.m.names[p].c_str(), A.r.th[p], A.r.Err(p)); fprintf(o, "\n"); };
        vline("adopted", nom);
        for (const char* nm : {"pre", "CRU565", "specCR", "fb1", "M1", "D11data", "D11dataOmega", "kBcons", "CBgm", "ghostA", "E1E2", "K0veto", "G3ref", "aeta0.0", "aeta1.0", "aeta2.0"}) vline(nm, V[nm]);
        if (!omegaAdopted) vline("omega", nomOmega);
        vline("spec_step0_E1E2_kE1free", s0);
        for (int c = 0; c < 3; ++c) for (size_t k = 0; k < chain[c].n.M.size(); ++k) vline(Form("chain%d_M%zu", c, k), chain[c].n.M[k]);
        fprintf(o, "# closure (adopted fit): region N predicted err_stageA err_CB S10_term\n");
        fprintf(o, "closure CR-K %.1f %.2f %.2f %.2f %.2f\n", pK.N, pK.mu, std::sqrt(pK.vfit), std::sqrt(pK.vcb), std::max(std::fabs(pK.mu - pK.N), std::sqrt(pK.vfit + pK.vcb)));
        fprintf(o, "closure CR-L %.1f %.2f %.2f %.2f %.2f\n", pL.N, pL.mu, std::sqrt(pL.vfit), std::sqrt(pL.vcb), std::max(std::fabs(pL.mu - pL.N), std::sqrt(pL.vfit + pL.vcb)));
        fclose(o);
    }
    {
        TFile* fo = TFile::Open(Form("%s/stageA%s.root", kCntOut, sfx.c_str()), "RECREATE");
        std::vector<double> e = F10Edges();
        for (int h = 0; h < kNH; ++h) {
            const HIn& x = H[h];
            TH1D* hd = new TH1D(Form("GDsmooth_%s", x.name.c_str()), "smoothed G_D, raw rows", kNB, e.data());
            TH1D* hc = new TH1D(Form("GCsmooth_%s", x.name.c_str()), "smoothed G_C, raw rows", kNB, e.data());
            TH1D* hb = new TH1D(Form("BA_%s", x.name.c_str()), "stage-A background, raw rows; error = stat (+) MC", kNB, e.data());
            TH1D* hbs = new TH1D(Form("BAstat_%s", x.name.c_str()), "stage-A background, error = stat", kNB, e.data());
            TH1D* hbm = new TH1D(Form("BAmc_%s", x.name.c_str()), "stage-A background, error = MC", kNB, e.data());
            for (int i = 0; i < kNB; ++i) {
                hd->SetBinContent(i + 1, x.D[i]); hc->SetBinContent(i + 1, x.C[i]);
                hb->SetBinContent(i + 1, binB[h][i]); hb->SetBinError(i + 1, std::hypot(binEs[h][i], binEmc[h][i]));
                hbs->SetBinContent(i + 1, binB[h][i]); hbs->SetBinError(i + 1, binEs[h][i]);
                hbm->SetBinContent(i + 1, binB[h][i]); hbm->SetBinError(i + 1, binEmc[h][i]);
            }
            hd->Write(); hc->Write(); hb->Write(); hbs->Write(); hbm->Write();
        }
        TMatrixDSym cs(mN.np), cm(mN.np), ct(mN.np);
        for (int p = 0; p < mN.np; ++p) for (int q = 0; q < mN.np; ++q) { cs(p, q) = Cstat[p][q]; cm(p, q) = Cmc[p][q]; ct(p, q) = Ctot[p][q]; }
        cs.Write("cov_stat"); cm.Write("cov_mc"); ct.Write("cov_total");
        fo->Close();
    }

    // -------------------------------------------------------------------------------------------
    // 12. tables
    // -------------------------------------------------------------------------------------------
    if (!full) { fclose(gOut); gOut = nullptr; f->Close(); fr->Close(); printf("done (tag %s)\n", tag); return; }
    {
        FILE* t = fopen(Form("%s/wp3_tableC1.txt", kTabDir), "w");
        FILE* x = fopen(Form("%s/wp3_tableC1.tex", kTabDir), "w");
        fprintf(t, "# Table C1: candidate control regions, data S0 (isBest >= 0), raw rows; simulation SMASH (mode A, w x V3, sigma_prod x PDG BR)\n");
        fprintf(t, "# C2-C5 %%, Geta %% = share in the raw simulated sum T12; (N-CB)/T12 = data/sim at SMASH scale (L = %.1f nb^-1, x 63/64 in T12)\n", kL);
        fprintf(x, "%% Table C1 (wp3_stageA.C): candidate control regions\n\\begin{tabular}{p{5.2cm}rrrrrp{5.6cm}}\n\\toprule\n");
        fprintf(x, "region & $N$ & CB/$N$ & C2--C5 [\\%%] & $G_\\eta$ [\\%%] & $\\frac{N-\\mathrm{CB}}{T_{12}}$ & verdict \\\\\n\\midrule\n");
        for (const auto& r : c1) {
            fprintf(t, "%-44s | %10s | %6s | %8s | %6s | %8s | %s\n", r.reg.c_str(), r.n.c_str(), r.cbn.c_str(), r.sig.c_str(), r.geta.c_str(), r.ratio.c_str(), r.verdict.c_str());
            auto esc = [](std::string s) { std::string o; for (char ch : s) { if (ch == '%') o += "\\%"; else if (ch == '_') o += "\\_"; else if (ch == '>') o += "$>$"; else if (ch == '<') o += "$<$"; else o += ch; } return o; };
            fprintf(x, "%s & %s & %s & %s & %s & %s & %s \\\\\n", esc(r.reg).c_str(), esc(r.n).c_str(), r.cbn.c_str(), esc(r.sig).c_str(), r.geta.c_str(), esc(r.ratio).c_str(), esc(r.verdict).c_str());
        }
        fprintf(x, "\\bottomrule\n\\end{tabular}\n");
        fclose(t); fclose(x);
    }
    {
        FILE* t = fopen(Form("%s/wp3_tableC2.txt", kTabDir), "w");
        FILE* x = fopen(Form("%s/wp3_tableC2.tex", kTabDir), "w");
        fprintf(t, "# Table C2: stage-A fits, data S0, F10, eta part fixed at a_eta(0) = %.2f; errors = stat (sandwich, event groups); dchi2 = vs the next simpler model\n", kAeta0);
        auto blockT = [&](const char* title, const NestOut& n) {
            fprintf(t, "\n# %s\n", title);
            for (size_t k = 0; k < n.M.size(); ++k) {
                const StageA& A = n.M[k];
                const double dc = k == 0 ? 0. : (k == 1 ? n.d01 : (k == 2 ? n.d12 : (k == 3 ? n.M[1].r.chi2P - A.r.chi2P : n.M[2].r.chi2P - A.r.chi2P)));
                const bool g = k < n.G.size() && n.G[k].ok;
                fprintf(t, "%-50s | %8.1f %4d %6.3f dchi2 %7.1f | MC-incl. %7.1f/%d |", (A.m.Tag() + ((int)k == n.adopt ? " [D11]" : "") + ((int)k == n.adoptData && n.adoptData != n.adopt ? " [D11 data-only]" : "")).c_str(),
                        A.r.chi2P, A.r.ndf, A.r.chi2P / A.r.ndf, dc, g ? n.G[k].chi2 : 0., g ? n.G[k].ndf : 0);
                for (int p = 0; p < A.m.np; ++p) fprintf(t, " %s=%.4f(%.4f)", A.m.names[p].c_str(), A.r.th[p], A.r.Err(p));
                fprintf(t, "\n");
            }
        };
        blockT("step 0, specification as written: E1 + E2, CR-U 0.565-0.700 + CR-K", n0);
        for (int c = 0; c < 3; ++c) blockT((chain[c].tag + ": " + chain[c].cr.tag + (fbKE1 ? ", E2 only" : "")).c_str(), chain[c].n);
        fprintf(t, "\n# ADOPTED %s, %s: value, stat (x %.3f), MC, total\n", mN.Tag().c_str(), crA.tag.c_str(), errScale);
        for (int p = 0; p < mN.np; ++p) fprintf(t, "%-8s %10.5f %9.5f %9.5f %9.5f\n", mN.names[p].c_str(), nom.r.th[p], sdv(Cstat, p), sdv(Cmc, p), sdv(Ctot, p));
        fprintf(t, "\n# correlation (total)\n");
        for (int p = 0; p < mN.np; ++p) { fprintf(t, "%-8s", mN.names[p].c_str()); for (int q = 0; q < mN.np; ++q) fprintf(t, " %6.3f", (Ctot[p][p] > 0 && Ctot[q][q] > 0) ? Ctot[p][q] / std::sqrt(Ctot[p][p] * Ctot[q][q]) : 0.); fprintf(t, "\n"); }
        fclose(t);
        // tex
        fprintf(x, "%% Table C2 (wp3_stageA.C): nested stage-A models (chi2_P/ndf, dchi2 against the next simpler model)\n");
        fprintf(x, "%% entries: chi2/ndf with the data variance only (in brackets: with the MC covariance of the templates, GLS); bold = adopted by D11\n");
        fprintf(x, "\\begin{tabular}{lrrrr}\n\\toprule\nmodel & spec CRs, E1+E2 & spec CRs, E2 & fall-back 1, E2 & fall-back 2, E2 \\\\\n\\midrule\n");
        const char* nm[5] = {"M0 scale $c_j$ only", "M1 + common tilt $b$", "M2 + quadratic $q$", "M3 per-slice tilts $b_j$", "M4 per-slice tilts + $q$"};
        for (int k = 0; k < 5; ++k) {
            fprintf(x, "%s", nm[k]);
            const NestOut* ns[4] = {&n0, &chain[0].n, &chain[1].n, &chain[2].n};
            for (int c = 0; c < 4; ++c) {
                const StageA& A = ns[c]->M[k];
                const bool g = k < (int)ns[c]->G.size() && ns[c]->G[k].ok;
                fprintf(x, " & %s%.0f/%d%s%s", k == ns[c]->adopt ? "\\textbf{" : "", A.r.chi2P, A.r.ndf, g ? Form(" (%.0f)", ns[c]->G[k].chi2) : "", k == ns[c]->adopt ? "}" : "");
            }
            fprintf(x, " \\\\\n");
        }
        fprintf(x, "\\bottomrule\n\\end{tabular}\n\n");
        fprintf(x, "\\begin{tabular}{lrrrrr}\n\\toprule\nparameter & value & stat & MC & total & $\\partial/\\partial a_\\eta$ \\\\\n\\midrule\n");
        for (int p = 0; p < mN.np; ++p) {
            if (p == mN.ikE1 && mN.eOnly == 1) continue;
            std::string n2 = mN.names[p];
            std::string tx = n2 == "kappa_B" ? "$\\kappa_B$" : (n2 == "k_E2" ? "$k_{E2}$" : (n2 == "k_E1" ? "$k_{E1}$" : (n2 == "c_om" ? "$c_\\omega$" : "$" + n2 + "$")));
            fprintf(x, "%s & %.4f & %.4f & %.4f & %.4f & %+.4f \\\\\n", tx.c_str(), nom.r.th[p], sdv(Cstat, p), sdv(Cmc, p), sdv(Ctot, p), dth[p]);
        }
        fprintf(x, "\\bottomrule\n\\end{tabular}\n");
        fclose(x);
    }
    {
        FILE* t = fopen(Form("%s/wp3_tableC3.txt", kTabDir), "w");
        FILE* x = fopen(Form("%s/wp3_tableC3.tex", kTabDir), "w");
        fprintf(t, "# Table C3: counted 'after' display, data S0 (isBest >= 0), RAW rows (x 64/63 / %.1f nb^-1 for nb). B_A = stage-A background (%s, %s);\n", kL, mN.Tag().c_str(), crA.tag.c_str());
        fprintf(t, "# eta(1) = all eta components at SMASH; eta(0.4) = x 0.4. DISPLAY COUNTS, NOT YIELDS (D20: yields only after WP4). Windows with fractional F10 bins\n");
        fprintf(t, "# (N of 'all M(ee)' from the exact window histograms where available). ratio = (N - k CB)/(B_A + 0.4 eta)\n");
        fprintf(x, "%% Table C3 (wp3_stageA.C): counted after display, data S0, raw rows; M(ee) >= 10 MeV (fitted) and all M(ee)\n");
        fprintf(x, "\\begin{tabular}{lrrrrrr}\n\\toprule\nregion [GeV/$c^2$] & $N$ & $k\\,$CB & $B_A$ & $\\eta$ (SMASH) & $N - k\\,\\mathrm{CB} - B_A$ & ratio \\\\\n\\midrule\n");
        auto blk = [&](const char* title, const std::vector<RegOut>& ro) {
            fprintf(t, "\n# %s\n# region | N | kCB err | B_A err_stat err_MC | eta(1) eta(0.4) Geta(1) C2-C5(1) | N-kCB-B_A err | ratio err\n", title);
            fprintf(x, "\\multicolumn{7}{l}{\\emph{%s}} \\\\\n", title);
            for (size_t q = 0; q < regs.size(); ++q) {
                const RegOut& o = ro[q];
                const double ex = o.N - o.CB - o.B, eex = std::sqrt(o.N + o.eCB * o.eCB + o.eBs * o.eBs + o.eBmc * o.eBmc);
                const double rr = (o.N - o.CB) / (o.B + o.eta04);
                const double err = rr * std::sqrt((o.N + o.eCB * o.eCB) / std::pow(o.N - o.CB, 2) + (o.eBs * o.eBs + o.eBmc * o.eBmc) / std::pow(o.B + o.eta04, 2));
                fprintf(t, "%-22s | %8.0f | %8.1f %6.1f | %9.1f %7.1f %7.1f | %8.1f %8.1f %8.1f %8.1f | %9.1f %7.1f | %.4f +- %.4f\n", regs[q].tag, o.N, o.CB, o.eCB, o.B, o.eBs,
                        o.eBmc, o.eta1, o.eta04, o.Geta1, o.sig1, ex, eex, rr, err);
                fprintf(x, "%s & %.0f & $%.0f \\pm %.0f$ & $%.0f \\pm %.0f$ & %.0f & $%.0f \\pm %.0f$ & $%.3f \\pm %.3f$ \\\\\n", regs[q].tag, o.N, o.CB, o.eCB, o.B,
                        std::hypot(o.eBs, o.eBmc), o.eta1, ex, eex, rr, err);
            }
        };
        blk("M(ee) >= 10 MeV/c2 (E2, the fitted class)", roE2);
        if (!fbKE1) blk("all M(ee)", roAll);
        fprintf(x, "\\bottomrule\n\\end{tabular}\n");
        fclose(t); fclose(x);
    }
    {
        FILE* t = fopen(Form("%s/wp3_A1_bins.txt", kTabDir), "w");
        fprintf(t, "# A1: share of the simulated signal classes C2-C5 in every CR bin (S0, SMASH scale), per configuration of the fall-back chain\n");
        fprintf(t, "# config hist M4lo M4hi region share_rawT12 share_smoothT12 share_fit(a_eta(0)) Neff(signal) signal_rows Geta_share_smooth\n");
        for (int c = 0; c < 3; ++c)
            for (const auto& b : chain[c].a1.bins)
                fprintf(t, "%s %-18s %.4f %.4f %s %.5f %.5f %.5f %.2f %.3f %.4f\n", c == 0 ? "spec" : (c == 1 ? "fb1" : "fb2"), b.h.c_str(), b.lo, b.hi, b.k ? "CRK" : "CRU",
                        b.sraw, b.ssm, b.sfit, b.neff, b.sig, b.geta);
        fclose(t);
    }
    {
        FILE* t = fopen(Form("%s/wp3_acceptance.txt", kTabDir), "w");
        fprintf(t, "# WP3 stage A acceptance tests (SPEC WP3) on the configuration they refer to\n");
        fprintf(t, "chi2/ndf<=1.5 adopted (MC template covariance included): %.1f/%d = %.3f -> %s (error scale %.3f); data-only Pearson %.1f/%d = %.3f -> %s\n", nomG.chi2, nomG.ndf,
                nomG.chi2 / nomG.ndf, nomG.chi2 / nomG.ndf <= 1.5 ? "PASS" : "FAIL", errScale, nom.r.chi2P, nom.r.ndf, chi2ndf, chi2ndf <= 1.5 ? "PASS" : "FAIL");
        for (int c = 0; c < 3; ++c) {
            const CfgRes& R = chain[c];
            fprintf(t, "A1 %s: max raw share %.4f (%d/%d bins > 1%%), smoothed %.4f (%d; %d with Neff>2), fitted-model %.4f (%d) -> %s\n", R.tag.c_str(), R.a1.maxRaw,
                    R.a1.nRaw, R.a1.nBins, R.a1.maxSm, R.a1.nSm, R.a1.nSmEff, R.a1.maxFit, R.a1.nFit, R.a1fail ? "FAIL" : "PASS");
            if (R.hasK) fprintf(t, "A2 %s: h(%s) %.3f +- %.3f, h(%s) %.3f +- %.3f, tag T %.2f +- %.2f, template h %.3f +- %.3f (gm %.3f +- %.3f) -> %s\n", R.tag.c_str(),
                                R.a2sl[0].c_str(), R.a2h[0], R.a2e[0], R.a2sl[1].c_str(), R.a2h[1], R.a2e[1], R.tagT, R.tagE, R.a2t.h, R.a2t.eh, R.a2tg.h, R.a2tg.eh,
                                R.a2fail ? "FAIL" : "PASS");
        }
        fprintf(t, "A3 signal C2-C5 injected x{0,1,2}: max shift %.3f sigma -> %s; with G_eta: %.3f sigma (information)\n", maxInjSig, maxInjSig < 0.3 ? "PASS" : "FAIL", maxInjAll);
        fprintf(t, "A4 per-slice tilts: dchi2 %.2f / 4 dof, p = %.3g -> %s; CR-U split 0.42: %.2f sigma -> %s (S4 inside: %.2f sigma)\n", a4dc, a4p, a4p > 0.01 ? "PASS" : "FAIL",
                splitZ, std::fabs(splitZ) < 2 ? "PASS" : "FAIL", splitS4Z);
        fprintf(t, "A5 K0s veto: max shift %.3f sigma -> %s\n", a5max, a5max < 1 ? "PASS" : "FAIL");
        fprintf(t, "k_E1 step 0: %.3f +- %.3f -> %s\n", kE1v, kE1e, fbKE1 ? "fall-back drop E1" : "ok");
        fclose(t);
    }

    // -------------------------------------------------------------------------------------------
    // 13. figures
    // -------------------------------------------------------------------------------------------
    bla::Style();
    // smoothing QA
    {
        TCanvas* c = new TCanvas("c_sm", "", 1800, 760);
        c->Divide(5, 2, 0.001, 0.001);
        std::vector<double> e = F10Edges();
        for (int h = 0; h < kNH; ++h) {
            c->cd(1 + (h % 2) * 5 + h / 2);
            gPad->SetLogy(); gPad->SetLeftMargin(0.17); gPad->SetRightMargin(0.02);
            const HIn& x = H[h];
            TH1D* rd = new TH1D(Form("rd%d", h), "", kNB, e.data()), *rc = new TH1D(Form("rc%d", h), "", kNB, e.data());
            TH1D* sd = new TH1D(Form("sd%d", h), "", kNB, e.data()), *sc = new TH1D(Form("sc%d", h), "", kNB, e.data());
            for (int i = 0; i < kNB; ++i) {
                const double w = 1e3 * (x.hi[i] - x.lo[i]);
                rd->SetBinContent(i + 1, x.s[lGD][i] / w); rd->SetBinError(i + 1, std::sqrt(x.v[lGD][i]) / w);
                rc->SetBinContent(i + 1, x.s[lGC][i] / w); rc->SetBinError(i + 1, std::sqrt(x.v[lGC][i]) / w);
                sd->SetBinContent(i + 1, x.D[i] / w); sc->SetBinContent(i + 1, x.C[i] / w);
            }
            rd->SetMinimum(0.05); rd->SetMaximum(2000);
            rd->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]"); rd->GetYaxis()->SetTitle("expected counts / (MeV/c^{2})");
            rd->GetXaxis()->SetNdivisions(505);
            rd->SetMarkerStyle(20); rd->SetMarkerSize(0.6); rd->SetMarkerColor(kAzure + 2); rd->SetLineColor(kAzure + 2);
            rc->SetMarkerStyle(24); rc->SetMarkerSize(0.6); rc->SetMarkerColor(kOrange + 7); rc->SetLineColor(kOrange + 7);
            sd->SetLineColor(kAzure + 2); sd->SetLineWidth(3); sc->SetLineColor(kOrange + 7); sc->SetLineWidth(3);
            rd->Draw("E1"); rc->Draw("E1 same"); sd->Draw("hist same ]["); sc->Draw("hist same ][");
            TLatex t; t.SetNDC(); t.SetTextSize(0.06);
            t.DrawLatex(0.22, 0.86, Form("S%d E%d", x.j + 1, x.e + 1));
            if (h == 0) {
                TLegend* l = new TLegend(0.40, 0.17, 0.97, 0.42);
                l->AddEntry(rd, "G_{D} simulated", "pe"); l->AddEntry(sd, "G_{D} exp(pol3)", "l");
                l->AddEntry(rc, "G_{C} simulated", "pe"); l->AddEntry(sc, "G_{C} exp(pol3)", "l");
                l->SetTextSize(0.05); l->Draw();
            }
        }
        SaveNoRm(c, Form("%s/wp3_smoothing_cnt_plainrun", kCntOut));
    }
    // Fig. C1 map
    {
        TH2* N = (TH2*)allio::Data(f, "S0", "all", "N", "B_map");
        TH2* CB = (TH2*)allio::Data(f, "S0", "all", "CB", "B_map");
        N->Add(CB, -1);
        TCanvas* c = new TCanvas("c_map", "", 1150, 960);
        gPad->SetRightMargin(0.17); gPad->SetLeftMargin(0.13); gPad->SetTopMargin(0.17);
        N->GetXaxis()->SetRangeUser(0.40, 0.70); N->GetYaxis()->SetRangeUser(0.28, 0.58);
        N->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]"); N->GetYaxis()->SetTitle("M(#pi^{+}#pi^{-}) [GeV/c^{2}]");
        N->GetZaxis()->SetTitle("N - CB [counts per 10 #times 10 (MeV/c^{2})^{2}]"); N->GetZaxis()->SetTitleOffset(1.35);
        N->SetMinimum(0);
        gStyle->SetPalette(kBird);
        N->Draw("COLZ");
        auto box = [&](double x1, double y1, double x2, double y2, int col, int sty, int w) {
            TBox* b = new TBox(x1, y1, x2, y2); b->SetFillStyle(0); b->SetLineColor(col); b->SetLineWidth(w); b->SetLineStyle(sty); b->Draw(); };
        box(0.565, 0.28, 0.700, 0.58, kGreen + 2, 2, 3);           // spec CR-U
        box(crA.uLo, 0.28, 0.700, 0.58, kGreen + 2, 1, 5);         // adopted CR-U
        box(0.400, 0.43, 0.5225, 0.58, kMagenta + 1, crA.useK ? 1 : 2, 3);   // CR-K
        box(0.5225, 0.28, 0.530, 0.58, kOrange + 1, 2, 3);         // CR-L
        box(kBlo, 0.28, kBhi, 0.43, kBlue + 1, 1, 3);              // bump window (S1-S4)
        box(kPlo, 0.28, kPhi, 0.58, kRed + 1, 1, 3);               // peak window
        for (double y : {0.34, 0.37, 0.39, 0.43}) { TLine* l = new TLine(0.40, y, 0.70, y); l->SetLineStyle(3); l->SetLineColor(kGray + 2); l->SetLineWidth(2); l->Draw(); }
        TGraph* g0 = new TGraph(), *g6 = new TGraph();
        int n0p = 0;
        for (double m = 0.28; m <= kEdgeLimit + 1e-6; m += 0.0005) { const double e0 = Edge(m, 0.), e6 = Edge(m, 0.060); if (e0 > 0 && e6 > 0) { g0->SetPoint(n0p, e0, m); g6->SetPoint(n0p, e6, m); ++n0p; } }
        g0->SetLineColor(kBlack); g0->SetLineWidth(4); g0->Draw("L");
        g6->SetLineColor(kBlack); g6->SetLineWidth(4); g6->SetLineStyle(2); g6->Draw("L");
        TLatex t; t.SetTextSize(0.034); t.SetTextFont(62);
        t.SetTextColor(kGreen + 2); t.DrawLatex(0.615, 0.548, "normalisation");
        t.SetTextColor(kMagenta + 1); t.DrawLatex(0.405, 0.437, crA.useK ? "control" : "control (not used)");
        t.SetTextColor(kBlue + 1); t.DrawLatex(0.476, 0.290, "bump");
        t.SetTextColor(kRed + 1); t.SetTextAngle(90); t.DrawLatex(0.5505, 0.29, "peak"); t.SetTextAngle(0);
        t.SetTextColor(kOrange + 1); t.SetTextAngle(90); t.DrawLatex(0.5285, 0.29, "test"); t.SetTextAngle(0);
        TLegend* l = new TLegend(0.08, 0.835, 0.95, 0.995);
        l->SetNColumns(2);
        l->SetHeader("data, both pions in the vertex fit, all M(e^{+}e^{-}): N - CB");
        l->AddEntry(g0, "edge, M_{ee} = 0", "l"); l->AddEntry(g6, "edge, M_{ee} = 60 MeV/c^{2}", "l");
        TLine* lu = new TLine(); lu->SetLineColor(kGreen + 2); lu->SetLineWidth(5);
        l->AddEntry(lu, Form("normalisation region %.3f-0.700 (dashed: 0.565, free-scale fits)", crA.uLo), "l");
        l->SetTextSize(0.028); l->SetFillStyle(0); l->Draw();
        SaveNoRm(c, Form("%s/wp3_map_cnt_plain2", kFigDir));
    }
    // Fig. C2 ratio per slice and class
    {
        TCanvas* c = new TCanvas("c_rat", "", 2000, 560);
        c->Divide(5, 1, 0.001, 0.001);
        Model m1 = mN; m1.aEta = 1.0;
        std::vector<double> e = F10Edges();
        for (int h = 0; h < kNH; ++h) {
            if (!Used(H[h], mN)) continue;
            c->cd(1 + h / 2);
            gPad->SetBottomMargin(0.15);
            gPad->SetLeftMargin(0.17); gPad->SetRightMargin(0.02); gPad->SetTopMargin(0.08);
            const HIn& x = H[h];
            const double ylo = 0.6, yhi = 1.8;
            TH1D* fr0 = new TH1D(Form("fr%d", h), "", kNB, e.data());
            fr0->SetMinimum(ylo); fr0->SetMaximum(yhi);
            fr0->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]");
            fr0->GetYaxis()->SetTitle("(N - k CB) / (B_{A} + 0.4 #eta)");
            fr0->GetXaxis()->SetNdivisions(505); fr0->GetYaxis()->SetTitleSize(0.055);
            fr0->Draw("AXIS");
            if (Used(x, mN)) { TBox* bu = new TBox(crA.uLo, ylo, crA.uHi, yhi); bu->SetFillColorAlpha(kGreen - 9, 0.40); bu->Draw(); }
            if (x.crkOK && crA.useK && Used(x, mN)) { TBox* bk = new TBox(0.400, ylo, 0.5225, yhi); bk->SetFillColorAlpha(kMagenta - 9, 0.40); bk->Draw(); }
            TGraphErrors* gin = new TGraphErrors(), *gout = new TGraphErrors();
            TGraph* gexp = new TGraph(), *gu = new TGraph(), *gl = new TGraph();
            int ni = 0, no = 0, ne = 0;
            const double kE = nom.r.th[x.e == 0 ? mN.ikE1 : mN.ikE2];
            for (int i = 0; i < kNB; ++i) {
                const double B = binB[h][i], eta = EtaSum(x, i, mN), eta1 = EtaSum(x, i, m1);
                const double den = B + eta;
                if (!(den > 0) || x.N[i] < 5) continue;
                const double r = (x.N[i] - kE * x.CB[i]) / den;
                const double er = std::sqrt(x.N[i] + kE * kE * x.vCB[i]) / den;
                const double xc = 0.5 * (x.lo[i] + x.hi[i]);
                if (Used(x, mN) && InCR(x, i, crA)) { gin->SetPoint(ni, xc, r); gin->SetPointError(ni++, 0, er); }
                else { gout->SetPoint(no, xc, r); gout->SetPointError(no++, 0, er); }
                gexp->SetPoint(ne, xc, (B + eta1) / den);
                const double eb = std::hypot(binEmc[h][i], binEs[h][i]) / den;
                gu->SetPoint(ne, xc, 1 + eb); gl->SetPoint(ne++, xc, 1 - eb);
            }
            gu->SetLineColor(kGray + 2); gl->SetLineColor(kGray + 2); gu->SetLineWidth(2); gl->SetLineWidth(2);
            gu->Draw("L"); gl->Draw("L");
            TLine* one = new TLine(0.40, 1, 0.70, 1); one->SetLineStyle(2); one->Draw();
            gexp->SetLineColor(kRed + 1); gexp->SetLineWidth(3); gexp->Draw("L");
            gout->SetMarkerStyle(24); gout->SetMarkerSize(0.7); gout->SetMarkerColor(kBlack); gout->SetLineColor(kBlack); gout->Draw("P");
            gin->SetMarkerStyle(20); gin->SetMarkerSize(0.7); gin->SetMarkerColor(kBlack); gin->SetLineColor(kBlack); gin->Draw("P");
            TLatex t; t.SetNDC(); t.SetTextSize(0.060);
            t.DrawLatex(0.21, 0.84, Form("S%d E%d%s", x.j + 1, x.e + 1, Used(x, mN) ? "" : " (not fitted)"));
            if (h == 1) {
                TLegend* l = new TLegend(0.19, 0.58, 0.98, 0.82);
                l->AddEntry(gin, "data, bin of the fit", "pe"); l->AddEntry(gout, "data, outside the fit", "pe");
                l->AddEntry(gexp, "#eta at SMASH instead of x 0.4", "l"); l->AddEntry(gu, "#pm stage-A error", "l");
                l->SetTextSize(0.050); l->SetFillStyle(0); l->Draw();
            }
        }
        SaveNoRm(c, Form("%s/wp3_ratio", kCntOut));   // _cnt copy: scratch only
    }
    // Fig. C3 "after" spectra
    auto afterFig = [&](int esel, const char* tag) {
        std::vector<double> e = F10Edges();
        TH1D* hD = new TH1D(Form("aD%s", tag), "", kNB, e.data()), *hC = new TH1D(Form("aC%s", tag), "", kNB, e.data());
        TH1D* hG = new TH1D(Form("aG%s", tag), "", kNB, e.data()), *h2 = new TH1D(Form("a2%s", tag), "", kNB, e.data());
        TH1D* h3 = new TH1D(Form("a3%s", tag), "", kNB, e.data()), *h4 = new TH1D(Form("a4%s", tag), "", kNB, e.data());
        TH1D* h5 = new TH1D(Form("a5%s", tag), "", kNB, e.data()), *hd = new TH1D(Form("ad%s", tag), "", kNB, e.data());
        TH1D* h04 = new TH1D(Form("a04%s", tag), "", kNB, e.data()), *hB = new TH1D(Form("aB%s", tag), "", kNB, e.data());
        std::vector<double> vd(kNB, 0.), vb(kNB, 0.);
        for (int h = 0; h < kNH; ++h) {
            const HIn& x = H[h];
            if (esel >= 0 && x.e != esel) continue;
            const double kE = nom.r.th[x.e == 0 ? mN.ikE1 : mN.ikE2];
            for (int i = 0; i < kNB; ++i) {
                Vec g(mN.np, 0.);
                // split of B_A into the G_D and G_C parts
                const double c = nom.r.th[mN.ic + x.j];
                const double T = 1 + (mN.tilt ? nom.r.th[mN.ib + (mN.perSlice ? x.j : 0)] * x.x[i] : 0.) + (mN.quad ? nom.r.th[mN.iq] * x.x[i] * x.x[i] : 0.);
                const double fD = binB[h][i] > 0 ? c * T * (x.D[i] + (mN.omega ? (nom.r.th[mN.iom] - 1) * x.Dom[i] : 0.)) / binB[h][i] : 0.;
                hD->AddBinContent(i + 1, binB[h][i] * fD); hC->AddBinContent(i + 1, binB[h][i] * (1 - fD));
                hG->AddBinContent(i + 1, x.s[lGeta][i]); h2->AddBinContent(i + 1, x.s[lC2][i]); h3->AddBinContent(i + 1, x.s[lC3][i]);
                h4->AddBinContent(i + 1, x.s[lC4][i]); h5->AddBinContent(i + 1, x.s[lC5][i]);
                hd->AddBinContent(i + 1, x.N[i] - kE * x.CB[i]); vd[i] += x.N[i] + kE * kE * x.vCB[i];
                h04->AddBinContent(i + 1, binB[h][i] + EtaSum(x, i, mN)); hB->AddBinContent(i + 1, binB[h][i]);
            }
        }
        for (TH1D* hh : {hD, hC, hG, h2, h3, h4, h5, hd, h04, hB}) for (int i = 0; i < kNB; ++i) {
            const double w = e[i + 1] - e[i];
            const double conv = (hh == hd) ? kRawToNb / w : 1.0 / (kNbToRaw * w);
            if (hh == hd) hh->SetBinError(i + 1, std::sqrt(vd[i]) * conv); else hh->SetBinError(i + 1, 0.);
            hh->SetBinContent(i + 1, hh->GetBinContent(i + 1) * conv);
        }
        TCanvas* c = new TCanvas(Form("c_after%s", tag), "", 1100, 1000);
        TPad* p1 = new TPad(Form("p1%s", tag), "", 0, 0.30, 1, 1); p1->SetBottomMargin(0.015); p1->SetLeftMargin(0.15); p1->SetRightMargin(0.03); p1->Draw();
        TPad* p2 = new TPad(Form("p2%s", tag), "", 0, 0, 1, 0.30); p2->SetTopMargin(0.02); p2->SetBottomMargin(0.36); p2->SetLeftMargin(0.15); p2->SetRightMargin(0.03); p2->Draw();
        p1->cd();
        THStack* st = new THStack(Form("st%s", tag), "");
        const int cols[7] = {kAzure - 4, kOrange - 3, kViolet - 4, kGreen + 1, kCyan + 1, kPink + 1, kRed + 1};
        TH1D* comps[7] = {hD, hC, hG, h2, h3, h4, h5};
        for (int k = 0; k < 7; ++k) { comps[k]->SetFillColor(cols[k]); comps[k]->SetLineColor(cols[k]); st->Add(comps[k], "hist"); }
        hd->SetMarkerStyle(20); hd->SetMarkerColor(kBlue + 1); hd->SetLineColor(kBlue + 1); hd->SetMarkerSize(1.0);
        h04->SetLineColor(kBlack); h04->SetLineWidth(3); h04->SetLineStyle(2);
        const double ymax = 1.30 * std::max(hd->GetMaximum(), st->GetMaximum());
        st->SetMaximum(ymax); st->SetMinimum(0);
        st->Draw();
        st->GetYaxis()->SetTitle("d#sigma/dM [nb/(GeV/c^{2})]"); st->GetYaxis()->SetTitleSize(0.050); st->GetYaxis()->SetLabelSize(0.045);
        st->GetXaxis()->SetLabelSize(0); st->GetYaxis()->SetTitleOffset(1.35);
        h04->Draw("hist same"); hd->Draw("E1 same");
        TBox* bu = new TBox(crA.uLo, 0, crA.uHi, ymax * 0.025); bu->SetFillColor(kGreen + 2); bu->Draw();
        TLegend* l = new TLegend(0.17, 0.36, 0.56, 0.79);
        l->AddEntry(hd, "data N - k CB", "pe");
        l->AddEntry(hD, "G_{D} renormalised (stage A)", "f"); l->AddEntry(hC, "G_{C} renormalised (stage A)", "f");
        l->AddEntry(hG, "G_{#eta} (SMASH)", "f"); l->AddEntry(h2, "C2 bump, #pi^{0} Dalitz (SMASH)", "f"); l->AddEntry(h3, "C3 bump, conversion (SMASH)", "f");
        l->AddEntry(h4, "C4 peak, conversion (SMASH)", "f"); l->AddEntry(h5, "C5 #eta #rightarrow #pi^{+}#pi^{-}e^{+}e^{-} (SMASH)", "f");
        l->AddEntry(h04, "total with all #eta x 0.4", "l");
        l->SetTextSize(0.031); l->SetFillStyle(0); l->Draw();
        TLatex t; t.SetNDC(); t.SetTextSize(0.038);
        t.DrawLatex(0.18, 0.87, Form("S0, isBest #geq 0, %s", esel < 0 ? "all M(e^{+}e^{-})" : (esel == 0 ? "M(e^{+}e^{-}) < 10 MeV/c^{2}" : "M(e^{+}e^{-}) #geq 10 MeV/c^{2}")));
        t.DrawLatex(0.18, 0.81, Form("L = %.1f nb^{-1}, data x 64/63; green bar: CR-U", kL));
        p2->cd();
        TH1D* r1 = (TH1D*)hd->Clone(Form("r1%s", tag)), *r04 = (TH1D*)hd->Clone(Form("r04%s", tag));
        TH1D* den1 = (TH1D*)hB->Clone(Form("den1%s", tag));
        for (TH1D* hh : {hG, h2, h3, h4, h5}) den1->Add(hh);
        for (int i = 1; i <= kNB; ++i) {
            const double d1 = den1->GetBinContent(i), d04 = h04->GetBinContent(i);
            r1->SetBinContent(i, d1 > 0 ? hd->GetBinContent(i) / d1 : 0); r1->SetBinError(i, d1 > 0 ? hd->GetBinError(i) / d1 : 0);
            r04->SetBinContent(i, d04 > 0 ? hd->GetBinContent(i) / d04 : 0); r04->SetBinError(i, d04 > 0 ? hd->GetBinError(i) / d04 : 0);
        }
        r1->SetMinimum(0.7); r1->SetMaximum(1.3);
        r1->GetYaxis()->SetTitle("data / sim"); r1->GetYaxis()->SetNdivisions(505);
        r1->GetYaxis()->SetTitleSize(0.11); r1->GetYaxis()->SetLabelSize(0.10); r1->GetYaxis()->SetTitleOffset(0.60);
        r1->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]"); r1->GetXaxis()->SetTitleSize(0.12); r1->GetXaxis()->SetLabelSize(0.10);
        r1->SetMarkerColor(kRed + 1); r1->SetLineColor(kRed + 1); r1->SetMarkerStyle(24);
        r04->SetMarkerColor(kBlack); r04->SetLineColor(kBlack); r04->SetMarkerStyle(20);
        r1->Draw("E1"); r04->Draw("E1 same");
        TLine* one = new TLine(0.40, 1, 0.70, 1); one->SetLineStyle(2); one->Draw();
        TLegend* l2 = new TLegend(0.17, 0.74, 0.75, 0.97); l2->SetNColumns(2);
        l2->AddEntry(r04, "all #eta x 0.4", "pe"); l2->AddEntry(r1, "#eta at SMASH", "pe"); l2->SetTextSize(0.09); l2->SetFillStyle(0); l2->Draw();
        SaveNoRm(c, Form("%s/wp3_after%s", kCntOut, tag));   // _cnt copy: scratch only
    };
    afterFig(fbKE1 ? 1 : -1, fbKE1 ? "_E2" : "");

    fclose(gOut); gOut = nullptr;
    f->Close(); fr->Close();
    printf("done\n");
}
