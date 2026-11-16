// wp6_ana_cnt.C -- COPY of wp6_ana.C (figure agent B, 2026-11-15): "counts" instead of "rows" in the axis titles and the
// sigma_tr annotation of Figs. wp6_E4, wp6_edge, wp6_tag (-> docs/latex/figs_all/wp6_E4_cnt, wp6_edge_cnt, wp6_tag_cnt .pdf/.png;
// fix_wp6_E4_unrot_cnt.pdf = pdftocairo -pdf of wp6_E4_cnt.pdf, as fix_wp6_E4_unrot.pdf of wp6_E4.pdf). Its log, results.txt and
// Table E3 go to /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figB/out/wp6/ (WP6 files untouched); same inputs. Nothing else changed.
// Run: root -l -b -q -e 'gSystem->SetBuildDir("<scratch>/aclic_wp6_cnt", kTRUE)' 'plots/all_report/wp6_ana_cnt.C+O'
// ============================================================================
// wp6_ana.C  --  WP6 side bands and cross-checks (SPEC report_all), prefix wp6
// ============================================================================
// Written by WP6, 2026-11-14. Inputs (all read-only):
//   hists_v1.root (WP0: data F10 and CB replicas, simulation templates and 100 replicas),
//   wp5/wp5_egcov.root (event-group covariances), wp5/results.txt (F-main values, for the check),
//   wp3/stageA.txt (frozen stage A, used as a variant), wp6/wp6_rows.root and wp6/wp6_fill.root (wp6_fill.C).
// Library: wp5_lib.h (= the fitter WP4 validated, general form) and wp4_lib.h (smoothing).
//
// What is computed (all on S0, isBest >= 0, raw data rows; yields x 64/63):
//   0. F-main (WP5 nominal: F10 E2, kappa = 1 +- 0.2) refitted and its 100 replica fits (the "other eta" inputs);
//      the values must equal WP5 results.txt.
//   1. Stage A re-run at the WP5 eta values (spec WP3: "rerun at the WP5 value"): frozen configuration of WP3
//      (model M1 c_j (1 + b x), E2 only, CR-U 0.575-0.700, no CR-K), eta terms fixed to F-main, smoothing and
//      fitter of WP4. The background B_A(window) = c_j * integral of (1 + b x)[G_D + kappa_B G_C] (smoothed).
//   2. B4 (bump, E2):  Y = (64/63) {(N - CB)_B - [(N - CB)_U - eta_U] B_A(B)/B_A(U) - [a kappa C4 + a rho C5 + a g G_eta]_B}
//      U = 0.575-0.700 (the adopted CR-U); eta_U = all eta classes in U at the F-main values (WP4 implementation).
//      All M(ee): Y_E2 / f_E2, f_E2 = (C2 + kappa C3)_E2 / (C2 + kappa C3)_all in B (simulation, kappa of F-main).
//   3. Peak side bands (all M(ee)): Y = (64/63) [(N - CB)_P - R' (N - CB)_SB] / (1 - R' f),
//      R' = [B_A + a (C2 + kappa C3) + a g G_eta](P) / [same](SB), f = (C4 + C5)(SB) / (C4 + C5)(P): the leakage of
//      the peak into the side bands, iterated out (the closed form is the limit of the iteration, shown).
//   4. Closure: 13 displaced windows (centres 0.575-0.695), each removed from CR-U before the stage-A refit;
//      pulls of the side-band estimator; sigma_tr = sigma_stat sqrt(max(0, rms^2 - 1)). CR-K and CR-L prediction pulls.
//   5. Errors: stat = data N (event-group covariance) + CB (200 replicas) + stage A (sandwich, numerical gradient)
//      + F-main (sandwich, total derivative including the stage-A refit); MC = RMS over the 100 replica chains
//      (F-main replica r -> stage A with replica r -> estimator with the replica-r window sums).
//   6. Table E3, Fig. E4 (+ wp6_edge, wp6_tag figures from wp6_fill.root).
//
// Run (from /home/przygoda/HADES/FAT/FAT; about 1-2 min, 3 GB):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_wp6", kTRUE)' \
//        'plots/all_report/wp6_ana.C+O'
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <functional>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "TBox.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TNamed.h"
#include "TPad.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include "wp5_lib.h"

namespace wp6 {
using namespace alldefs;
typedef std::vector<double> Vec;
using wp5::Templ;
using wp5::DataSet;
using wp5::Bins;
using wp5::Spec;
using wp5::FitRes;
using wp5::Opt;

const char* const kOut = "/hdd3/przygoda/HADES/FAT_out/all_report/wp6";
const char* const kFig = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_all";
const char* const kTab = "/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figB/out/wp6";   // _cnt copy: tables to scratch
const char* const kCntOut = "/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figB/out/wp6";   // _cnt copy: log and results.txt to scratch
const char* const kStageA = "/hdd3/przygoda/HADES/FAT_out/all_report/wp3/stageA.txt";
const char* const kWP5 = "/hdd3/przygoda/HADES/FAT_out/all_report/wp5/results.txt";
const int kNRep = 100, kNRepCB = 200, kNCell = 10, kNFine = 320;
const double kYf = kPT3;
FILE* gLog = nullptr;
FILE* gRes = nullptr;
#define P(...) do { printf(__VA_ARGS__); if (gLog) fprintf(gLog, __VA_ARGS__); } while (0)
#define RS(...) do { if (gRes) fprintf(gRes, __VA_ARGS__); } while (0)

// ---------------------------------------------------------------------------------------------------------------
// windows (identical list to wp6_fill.C, read back from wp6_fill.root)
// ---------------------------------------------------------------------------------------------------------------
struct Win { std::string tag; double lo, hi; };
std::vector<Win> gW;
int nW = 0;
int WI(const std::string& t) {
    for (int i = 0; i < nW; ++i) if (gW[i].tag == t) return i;
    printf("wp6: window %s not found\n", t.c_str());
    return -1;
}
inline int DI(int cell, int w) { return cell * nW + w; }
enum { wC2 = 0, wC3, wC4, wC5, wGeta, wGD, wGC, wT12, kNWc };
inline double SW(const Vec& s, int comp, int cell, int w) { return s[((size_t)comp * kNCell + cell) * nW + w]; }
std::vector<int> Cells(int e, int slice = -1) {   // e: 0 E1, 1 E2, 2 both; slice -1 = all
    std::vector<int> c;
    for (int j = 0; j < kNSlice; ++j) {
        if (slice >= 0 && j != slice) continue;
        if (e == 0 || e == 2) c.push_back(2 * j);
        if (e == 1 || e == 2) c.push_back(2 * j + 1);
    }
    return c;
}
const double kPc = 0.5439;
const double kDSl1 = 0.5250 - kPc, kDSh2 = 0.5628 - kPc;
const int kNDisp = 13;
inline double DispCentre(int k) { return 0.575 + 0.010 * k; }

// ---------------------------------------------------------------------------------------------------------------
// replica generator (copy of wp0_fill.C) for the CB replicas of the data windows
// ---------------------------------------------------------------------------------------------------------------
inline uint64_t Mix(uint64_t x) {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}
inline int Pois1(uint64_t h) {
    const double u = (double)(h >> 11) * (1.0 / 9007199254740992.0);
    double p = 0.36787944117144233, c = p;
    int k = 0;
    while (u >= c && k < 30) { ++k; p /= k; c += p; }
    return k;
}
inline void GroupPois(uint64_t seed, int tree, int grp, int n, double* k) {
    const uint64_t base = Mix(seed ^ Mix(((uint64_t)(uint32_t)tree << 32) | (uint64_t)(uint32_t)grp));
    for (int r = 0; r < n; ++r) k[r] = Pois1(Mix(base + 0x632BE59BD9B4E019ULL * (uint64_t)(r + 1)));
}
const uint64_t kSeedCB6 = 0xCB2026ULL;

// ---------------------------------------------------------------------------------------------------------------
// data window sums: N (event-group covariance) and CB (spread rule, 200 replicas)
// ---------------------------------------------------------------------------------------------------------------
inline double FracFine(int i, double a, double b) {
    if (a < -1e8) return 1.0;
    const double x1 = kFineLo + i * kFineW, x2 = x1 + kFineW;
    const double o = std::min(x2, b) - std::max(x1, a);
    return o > 0 ? o / kFineW : 0.;
}
// pp, mm: [cell * kNFine + i]
Vec CBWin(const double* pp, const double* mm) {
    Vec out((size_t)kNCell * nW, 0.);
    std::vector<double> tp(kNFine, 0.), tm(kNFine, 0.);
    for (int c = 0; c < kNCell; ++c) for (int i = 0; i < kNFine; ++i) { tp[i] += pp[c * kNFine + i]; tm[i] += mm[c * kNFine + i]; }
    const int nw = kNFine / kCBreb;
    std::vector<double> cbW(nw, 0.);
    for (int j = 0; j < nw; ++j) {
        double a = 0, b = 0;
        for (int i = j * kCBreb; i < (j + 1) * kCBreb; ++i) { a += tp[i]; b += tm[i]; }
        cbW[j] = 2. * std::sqrt(a * b);
    }
    std::vector<double> cbF(kNFine, 0.);
    for (int i = 0; i < kNFine; ++i) {
        const double x = kFineLo + (i + 0.5) * kFineW;
        const double t = (x - (kFineLo + 0.5 * kCBreb * kFineW)) / (kCBreb * kFineW);
        double v;
        if (t <= 0) v = cbW[0];
        else if (t >= nw - 1) v = cbW[nw - 1];
        else { const int j0 = (int)t; const double f = t - j0; v = (1 - f) * cbW[j0] + f * cbW[j0 + 1]; }
        cbF[i] = v / kCBreb;
    }
    std::vector<double> La(kNBand, 0.), Lc((size_t)kNCell * kNBand, 0.);
    for (int i = 0; i < kNFine; ++i) La[BandOfFine(i)] += tp[i] + tm[i];
    for (int c = 0; c < kNCell; ++c) for (int i = 0; i < kNFine; ++i) Lc[(size_t)c * kNBand + BandOfFine(i)] += pp[c * kNFine + i] + mm[c * kNFine + i];
    for (int c = 0; c < kNCell; ++c)
        for (int w = 0; w < nW; ++w) {
            double s = 0;
            for (int i = 0; i < kNFine; ++i) {
                const double u = FracFine(i, gW[w].lo, gW[w].hi);
                if (u <= 0) continue;
                const int b = BandOfFine(i);
                if (La[b] > 0) s += cbF[i] * u * Lc[(size_t)c * kNBand + b] / La[b];
            }
            out[DI(c, w)] = s;
        }
    return out;
}
struct DataW {
    Vec N, CB;
    std::vector<double> covN, covCB;   // dense (kNCell nW)^2
    Vec cbRepMean;
    long nUS = 0, nLS = 0;
};
bool LoadDataW(DataW& D) {
    TFile* f = TFile::Open(Form("%s/wp6_rows.root", kOut));
    if (!f) return false;
    TTree* t = (TTree*)f->Get("rows");
    Int_t k, grp; Float_t m4, mpipi, mee;
    t->SetBranchAddress("k", &k); t->SetBranchAddress("grp", &grp);
    t->SetBranchAddress("m4", &m4); t->SetBranchAddress("mpipi", &mpipi); t->SetBranchAddress("mee", &mee);
    const int nI = kNCell * nW;
    D.N.assign(nI, 0.); D.covN.assign((size_t)nI * nI, 0.);
    std::vector<double> pp((size_t)kNCell * kNFine, 0.), mm((size_t)kNCell * kNFine, 0.);
    std::vector<double> rpp((size_t)kNRepCB * kNCell * kNFine, 0.), rmm((size_t)kNRepCB * kNCell * kNFine, 0.);
    std::map<int, int> grpCnt;
    int curK = -1, curG = -1;
    std::vector<double> kr(kNRepCB, 1.);
    auto flushN = [&]() {
        for (auto& a : grpCnt) for (auto& b : grpCnt) D.covN[(size_t)a.first * nI + b.first] += (double)a.second * b.second;
        grpCnt.clear();
    };
    const Long64_t n = t->GetEntries();
    for (Long64_t i = 0; i < n; ++i) {
        t->GetEntry(i);
        if (k != curK || grp != curG) {
            if (curK == 0) flushN();
            curK = k; curG = grp;
            if (k > 0) GroupPois(kSeedCB6, k, grp, kNRepCB, kr.data());
        }
        const int s = Slice(mpipi), e = Eclass(mee);
        if (s < 0 || e < 0) continue;
        const int cell = 2 * s + e;
        if (k == 0) {
            ++D.nUS;
            for (int w = 0; w < nW; ++w)
                if (m4 >= gW[w].lo && m4 < gW[w].hi) { D.N[DI(cell, w)] += 1; grpCnt[DI(cell, w)] += 1; }
        } else {
            ++D.nLS;
            const int fi = FineBin(m4);
            if (fi < 0) continue;
            std::vector<double>& A = k == 1 ? pp : mm;
            std::vector<double>& RA = k == 1 ? rpp : rmm;
            A[(size_t)cell * kNFine + fi] += 1;
            for (int r = 0; r < kNRepCB; ++r) RA[((size_t)r * kNCell + cell) * kNFine + fi] += kr[r];
        }
    }
    if (curK == 0) flushN();
    f->Close();
    D.CB = CBWin(pp.data(), mm.data());
    std::vector<Vec> rep(kNRepCB);
    D.cbRepMean.assign(nI, 0.);
    for (int r = 0; r < kNRepCB; ++r) {
        rep[r] = CBWin(&rpp[(size_t)r * kNCell * kNFine], &rmm[(size_t)r * kNCell * kNFine]);
        for (int a = 0; a < nI; ++a) D.cbRepMean[a] += rep[r][a] / kNRepCB;
    }
    D.covCB.assign((size_t)nI * nI, 0.);
    for (int r = 0; r < kNRepCB; ++r)
        for (int a = 0; a < nI; ++a) {
            const double da = rep[r][a] - D.cbRepMean[a];
            if (da == 0) continue;
            for (int b = 0; b < nI; ++b) D.covCB[(size_t)a * nI + b] += da * (rep[r][b] - D.cbRepMean[b]) / (kNRepCB - 1);
        }
    return true;
}

// ---------------------------------------------------------------------------------------------------------------
// WP3 smoothing (copied from plots/all_report/wp3_stageA.C, Smooth/BinMoments, unchanged): exp(pol3) in
// t = (M4 - 0.55)/0.15, Poisson likelihood of the expected rows, fitted from the first bin with simulated content
// (the kinematic threshold of the histogram); replicas start from the nominal domain and parameters.
// (The WP4 smoothing has the same likelihood maximum per domain but always starts at 0.400.)
// ---------------------------------------------------------------------------------------------------------------
namespace w3 {
const double kGLx[4] = {-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526};
const double kGLw[4] = {0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538};
inline double Tsm(double m) { return (m - 0.55) / 0.15; }
struct SmRes { bool ok = false; int dom = 0; double p[4] = {0, 0, 0, 0}; Vec f; double chi2 = 0; int ndf = 0; };
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
        bla::Mat A = bla::Zero(4, 4); Vec b(4, 0.);
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
        bla::Mat H = bla::Zero(4, 4); Vec g(4, 0.);
        double mom[7];
        for (int i = dom; i < n; ++i) {
            BinMoments(lo[i], hi[i], p, mom);
            const double f = std::max(mom[0], 1e-300), a = 1.0 - y[i] / f, c2 = y[i] / (f * f);
            for (int k = 0; k < 4; ++k) {
                g[k] += a * mom[k];
                for (int l = 0; l < 4; ++l) H[k][l] += a * mom[k + l] + c2 * mom[k] * mom[l];
            }
        }
        bla::Mat Hi = H;
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
}  // namespace w3
// smoothing parameters of the 10 F10 histograms (G_D = 0, G_C = 1) in the WP3 form
struct SmSet { double start[10][2]; int dom[10][2]; double p[10][2][4]; };
inline double DensW3(const SmSet& s, int h, int c, double m) {
    if (m < s.start[h][c]) return 0.;
    const double t = w3::Tsm(m);
    const double* p = s.p[h][c];
    return std::exp(std::min(700., p[0] + t * (p[1] + t * (p[2] + t * p[3]))));
}
// replace the smoothed bin integrals of T (F10 histograms) by the WP3 smoothing; nominal (ref == nullptr) or replica
void SmoothTemplW3(Templ& T, SmSet& out, const SmSet* ref, double* chi2 = nullptr) {
    const std::vector<double> e = F10Edges();
    Vec lo(kNF10), hi(kNF10);
    for (int i = 0; i < kNF10; ++i) { lo[i] = e[i]; hi[i] = e[i + 1]; }
    for (int h = 0; h < 10; ++h)
        for (int c = 0; c < 2; ++c) {
            const Vec& src = c == 0 ? T.gd[h] : T.gc[h];
            const Vec& var = c == 0 ? T.vgd[h] : T.vgc[h];
            Vec y(kNF10), v(kNF10);
            for (int i = 0; i < kNF10; ++i) { y[i] = src[i + 1]; v[i] = var[i + 1]; }
            const w3::SmRes r = ref ? w3::Smooth(y, v, lo, hi, ref->dom[h][c], ref->p[h][c]) : w3::Smooth(y, v, lo, hi);
            out.dom[h][c] = r.dom;
            out.start[h][c] = r.dom < kNF10 ? lo[r.dom] : 1e9;
            for (int k = 0; k < 4; ++k) out.p[h][c][k] = r.p[k];
            Vec& dst = c == 0 ? T.gdI[h] : T.gcI[h];
            for (int i = 0; i < kNF10; ++i) dst[i + 1] = r.f[i];
            if (chi2) chi2[2 * h + c] = r.ndf > 0 ? r.chi2 / r.ndf : 0.;
        }
}
inline double IntTiltW3(const SmSet& s, int h, int c, double lo, double hi, double b, double q) {
    if (!(hi > lo)) return 0.;
    lo = std::max(lo, s.start[h][c]);
    if (!(hi > lo)) return 0.;
    const int n = std::max(1, (int)std::ceil((hi - lo) / 0.0025 - 1e-9));
    const double w = (hi - lo) / n;
    double tot = 0;
    for (int k = 0; k < n; ++k) {
        const double cc = lo + (k + 0.5) * w;
        for (int g = 0; g < 8; ++g) {
            const double m = cc + 0.5 * w * wp4::kGLx[g];
            const double x = TiltX(m);
            tot += DensW3(s, h, c, m) * (1. + b * x + q * x * x) * wp4::kGLw[g] * 0.5 * w;
        }
    }
    return tot;
}

// ---------------------------------------------------------------------------------------------------------------
// smoothed background integrals
// ---------------------------------------------------------------------------------------------------------------
inline double IntTilt(const wp4::Smooth& sm, double lo, double hi, double b, double q) {
    if (sm.order < 0 || !(hi > lo)) return 0.;
    const int n = std::max(1, (int)std::ceil((hi - lo) / 0.0025 - 1e-9));
    const double h = (hi - lo) / n;
    double tot = 0;
    for (int k = 0; k < n; ++k) {
        const double c = lo + (k + 0.5) * h;
        for (int g = 0; g < 8; ++g) {
            const double m = c + 0.5 * h * wp4::kGLx[g];
            const double x = TiltX(m);
            tot += sm.Dens(m) * (1. + b * x + q * x * x) * wp4::kGLw[g] * 0.5 * h;
        }
    }
    return tot;
}
// stage A background of histogram h over [lo, hi): c_j int (1 + b x + q x^2) [G_D + kappa_B G_C]
inline double BAint(const Templ& T, const SmSet* sm, const Spec& S, const Vec& th, int h, double lo, double hi) {
    const double kB = S.iKB >= 0 ? th[S.iKB] : S.kBfix;
    const double b = S.iB >= 0 ? th[S.iB] : S.bFix;
    const double q = S.iQ >= 0 ? th[S.iQ] : 0.;
    const double c = th[S.iC[h / 2]];
    if (sm) return c * (IntTiltW3(*sm, h, 0, lo, hi, b, q) + kB * IntTiltW3(*sm, h, 1, lo, hi, b, q));
    return c * (IntTilt(T.sgd[h], lo, hi, b, q) + kB * IntTilt(T.sgc[h], lo, hi, b, q));
}
// ... in every (cell, window); the F10 histogram index h = 2 slice + class = cell
Vec BATable(const Templ& T, const SmSet* sm, const Spec& S, const Vec& th) {
    Vec t((size_t)kNCell * nW, 0.);
    for (int cell = 0; cell < kNCell; ++cell)
        for (int w = 0; w < nW; ++w) {
            if (gW[w].lo < -1e8) continue;
            t[DI(cell, w)] = BAint(T, sm, S, th, cell, gW[w].lo, gW[w].hi);
        }
    return t;
}
// WP3 frozen stage A (stageA.txt): smoothing exp(pol3) in t = (M4 - 0.55)/0.15 from 'start'; values c_j, b, kappa_B
struct WP3Smooth { double sD = 0, pD[4] = {0, 0, 0, 0}, sC = 0, pC[4] = {0, 0, 0, 0}; };
struct WP3Stage { Vec v; WP3Smooth sm[10]; Vec binBA[10]; bool ok = false; };
bool LoadWP3(WP3Stage& W) {
    std::ifstream in(kStageA);
    if (!in) return false;
    std::string line;
    for (int h = 0; h < 10; ++h) W.binBA[h].assign(kNF10 + 2, 0.);
    while (std::getline(in, line)) {
        std::istringstream is(line);
        std::string key; is >> key;
        if (key == "values") { W.v.clear(); double x; while (is >> x) W.v.push_back(x); }
        else if (key == "smooth") {
            std::string hn; is >> hn;
            int j = hn[5] - '1', e = hn[7] - '1';
            WP3Smooth& s = W.sm[2 * j + e];
            is >> s.sD >> s.pD[0] >> s.pD[1] >> s.pD[2] >> s.pD[3] >> s.sC >> s.pC[0] >> s.pC[1] >> s.pC[2] >> s.pC[3];
        } else if (key == "bin") {
            std::string hn; int x; double lo, hi, ba;
            is >> hn >> x >> lo >> hi >> ba;
            const int j = hn[5] - '1', e = hn[7] - '1';
            W.binBA[2 * j + e][x] = ba;
        }
    }
    W.ok = W.v.size() >= 9;
    return W.ok;
}
inline double WP3Dens(double m, double start, const double* p) {
    if (m < start) return 0.;
    const double t = (m - 0.55) / 0.15;
    return std::exp(p[0] + t * (p[1] + t * (p[2] + t * p[3])));
}
Vec BATableWP3(const WP3Stage& W, const Vec& v) {   // v = c_1..c_5, b, kappa_B, k_E1, k_E2
    Vec t((size_t)kNCell * nW, 0.);
    for (int cell = 0; cell < kNCell; ++cell) {
        const WP3Smooth& s = W.sm[cell];
        const double c = v[cell / 2], b = v[5], kB = v[6];
        for (int w = 0; w < nW; ++w) {
            if (gW[w].lo < -1e8) continue;
            const double lo = gW[w].lo, hi = gW[w].hi;
            const int n = std::max(1, (int)std::ceil((hi - lo) / 0.0025 - 1e-9));
            const double h = (hi - lo) / n;
            double tot = 0;
            for (int k = 0; k < n; ++k) {
                const double cc = lo + (k + 0.5) * h;
                for (int g = 0; g < 8; ++g) {
                    const double m = cc + 0.5 * h * wp4::kGLx[g];
                    tot += (1. + b * TiltX(m)) * (WP3Dens(m, s.sD, s.pD) + kB * WP3Dens(m, s.sC, s.pC)) * wp4::kGLw[g] * 0.5 * h;
                }
            }
            t[DI(cell, w)] = c * tot;
        }
    }
    return t;
}

// ---------------------------------------------------------------------------------------------------------------
// stage A (WP3 configuration) with the eta terms fixed
// ---------------------------------------------------------------------------------------------------------------
Spec SpecStageA(bool quad) { Spec S; Opt o; o.quad = quad; wp5::AddBkg(S, o, false, true); return S; }
std::function<double(int, int)> EtaFix(const Templ* T, const Vec& f) {
    return [T, f](int h, int x) {
        const double a = f[0], k = f[1], r = f[2], g = f[3];
        return a * T->s[wp5::kS2][h][x] + a * k * (T->s[wp5::kS3][h][x] + T->s[wp5::kS4][h][x]) + a * r * T->s[wp5::kS5][h][x] +
               a * g * T->s[wp5::kSeta][h][x];
    };
}
Bins FilterBins(const Bins& B, const std::function<bool(double, double)>& drop) {
    if (!drop) return B;
    Bins o; o.lay = B.lay;
    for (int i = 0; i < B.n; ++i) {
        if (drop(B.lo[i], B.hi[i])) continue;
        o.y.push_back(B.y[i]); o.cb.push_back(B.cb[i]); o.vcb.push_back(B.vcb[i]); o.x.push_back(B.x[i]);
        o.gd.push_back(B.gd[i]); o.gc.push_back(B.gc[i]); o.fix.push_back(B.fix[i]); o.lo.push_back(B.lo[i]); o.hi.push_back(B.hi[i]);
        o.j.push_back(B.j[i]); o.e.push_back(B.e[i]); o.h.push_back(B.h[i]); o.xb.push_back(B.xb[i]); o.cell.push_back(B.cell[i]);
        for (int t = 0; t < wp5::kNSig; ++t) { o.s[t].push_back(B.s[t][i]); o.vs[t].push_back(B.vs[t][i]); }
        ++o.n;
    }
    return o;
}
struct StageA { Spec S; Bins B; FitRes F; };
StageA FitStageA(const DataSet& D, const Templ& T, const Vec& thF, double lo, double hi, const std::function<bool(double, double)>& drop,
                 bool quad, const Vec* start, bool sandwich) {
    StageA A;
    A.S = SpecStageA(quad);
    Bins B;
    wp5::MakeBins(D, T, 0, 2, B, lo, hi, EtaFix(&T, thF));
    A.B = FilterBins(B, drop);
    A.F = wp5::Fit(A.B, A.S, D, start, sandwich);
    return A;
}

// ---------------------------------------------------------------------------------------------------------------
// estimator state and contexts
// ---------------------------------------------------------------------------------------------------------------
struct St { const Vec* N; const Vec* CB; const Vec* BA; const Vec* S; Vec thF; double kE2; };
typedef std::function<double(const St&)> Est;
struct Ctx {
    std::string tag;
    bool ok = false;
    StageA A;
    Vec BA0, thF0;
    double kE2 = 1;
    // F-main total derivative: stage A refitted at thF +- h
    Vec hF;
    std::vector<Vec> BAFp, BAFm;
    std::vector<double> kE2p, kE2m;
    // replica chain
    std::vector<Vec> BAr, thFr;
    std::vector<double> kE2r;
    std::vector<int> repIdx;
    // stage-A parameter gradient helper; smoothing: sm != nullptr = WP3 form (mode 1), else WP4 (mode 0)
    const Templ* T = nullptr;
    const SmSet* sm = nullptr;
    int mode = 1;
    bool wp3 = false;
    const WP3Stage* W3 = nullptr;
};
struct Glob {
    DataSet D0;
    Templ T0;
    std::vector<Templ> R0;
    FitRes FM;
    Spec SM;
    wp5::Replicas RepFM;
    DataW DW;
    Vec S0, Svar;
    std::vector<Vec> Srep;
    Vec thF0;
    Templ T3;                      // templates with the WP3 smoothing (nominal stage A of this package)
    std::vector<Templ> R3;
    SmSet sm3;
    std::vector<SmSet> sm3R;
    wp5::Mat CF;   // F-main stat covariance of (a, kappa, rho, g)
    wp5::Mat CFmc; // replica covariance of (a, kappa, rho, g)
};
Glob* G = nullptr;
inline Vec ThF(const Vec& th) { return Vec(th.begin(), th.begin() + 4); }

Ctx BuildCtx(const std::string& tag, double lo, double hi, const std::function<bool(double, double)>& drop, bool quad, bool withF, bool withRep,
             int mode = 1) {
    Ctx C; C.tag = tag; C.mode = mode; C.thF0 = G->thF0;
    const Templ& TT = mode == 1 ? G->T3 : G->T0;
    const std::vector<Templ>& RR = mode == 1 ? G->R3 : G->R0;
    C.T = &TT; C.sm = mode == 1 ? &G->sm3 : nullptr;
    C.A = FitStageA(G->D0, TT, C.thF0, lo, hi, drop, quad, nullptr, true);
    if (!C.A.F.ok) { P("  [%s] stage A FAILED\n", tag.c_str()); return C; }
    C.BA0 = BATable(TT, C.sm, C.A.S, C.A.F.th);
    C.kE2 = C.A.F.th[C.A.S.iK[1]];
    if (withF) {
        C.hF.assign(4, 0.);
        for (int k = 0; k < 4; ++k) {
            C.hF[k] = 0.25 * std::sqrt(std::max(G->CF[k][k], 1e-12));
            Vec tp = C.thF0, tm = C.thF0;
            tp[k] += C.hF[k]; tm[k] -= C.hF[k];
            const StageA ap = FitStageA(G->D0, TT, tp, lo, hi, drop, quad, &C.A.F.th, false);
            const StageA am = FitStageA(G->D0, TT, tm, lo, hi, drop, quad, &C.A.F.th, false);
            C.BAFp.push_back(BATable(TT, C.sm, ap.S, ap.F.th)); C.BAFm.push_back(BATable(TT, C.sm, am.S, am.F.th));
            C.kE2p.push_back(ap.F.th[ap.S.iK[1]]); C.kE2m.push_back(am.F.th[am.S.iK[1]]);
        }
    }
    if (withRep) {
        for (int r = 0; r < (int)RR.size(); ++r) {
            if (!G->RepFM.fit[r].ok) continue;
            const Vec tf = ThF(G->RepFM.fit[r].th);
            const StageA ar = FitStageA(G->D0, RR[r], tf, lo, hi, drop, quad, &C.A.F.th, false);
            if (!ar.F.ok) continue;
            C.BAr.push_back(BATable(RR[r], mode == 1 ? &G->sm3R[r] : nullptr, ar.S, ar.F.th));
            C.thFr.push_back(tf);
            C.kE2r.push_back(ar.F.th[ar.S.iK[1]]);
            C.repIdx.push_back(r);
        }
    }
    C.ok = true;
    return C;
}
// context from the frozen WP3 stage A (no refit; stat and MC from the WP3 parameter covariances)
Ctx BuildCtxWP3(const WP3Stage& W) {
    Ctx C; C.tag = "WP3frozen"; C.thF0 = G->thF0; C.wp3 = true; C.W3 = &W;
    C.BA0 = BATableWP3(W, W.v);
    C.kE2 = W.v[8];
    C.ok = true;
    return C;
}

struct ERes {
    std::string name;
    double v = 0, eN = 0, eCB = 0, eA = 0, eF = 0, stat = 0, mc = 0, shift = 0, bc = 0, mc68 = 0, tot = 0;
    int nrep = 0;
    Vec rep;
};
ERes Evaluate(const std::string& name, const Est& f, const Ctx& C, bool withMC, const wp5::Mat* wp3CovStat = nullptr,
              const wp5::Mat* wp3CovMC = nullptr) {
    ERes o; o.name = name;
    const DataW& D = G->DW;
    St s0{&D.N, &D.CB, &C.BA0, &G->S0, C.thF0, C.kE2};
    o.v = f(s0);
    const int nI = kNCell * nW;
    // data N and CB (the estimators are linear in N and CB: exact finite differences)
    {
        Vec gN(nI, 0.), gC(nI, 0.);
        Vec N = D.N, CB = D.CB;
        for (int i = 0; i < nI; ++i) {
            N[i] += 1.; St s{&N, &D.CB, &C.BA0, &G->S0, C.thF0, C.kE2}; gN[i] = f(s) - o.v; N[i] -= 1.;
            CB[i] += 1.; St t{&D.N, &CB, &C.BA0, &G->S0, C.thF0, C.kE2}; gC[i] = f(t) - o.v; CB[i] -= 1.;
        }
        double vN = 0, vC = 0;
        for (int a = 0; a < nI; ++a) {
            if (gN[a] != 0) for (int b = 0; b < nI; ++b) if (gN[b] != 0) vN += gN[a] * D.covN[(size_t)a * nI + b] * gN[b];
            if (gC[a] != 0) for (int b = 0; b < nI; ++b) if (gC[b] != 0) vC += gC[a] * D.covCB[(size_t)a * nI + b] * gC[b];
        }
        o.eN = std::sqrt(std::max(0., vN)); o.eCB = std::sqrt(std::max(0., vC));
    }
    // stage-A parameters
    if (!C.wp3) {
        const Spec& S = C.A.S;
        const Vec& th = C.A.F.th;
        const int np = S.np();
        Vec g(np, 0.);
        for (int k = 0; k < np; ++k) {
            if (S.fixed[k] || C.A.F.Cstat[k][k] <= 0) continue;
            const double h = 0.05 * std::sqrt(C.A.F.Cstat[k][k]);
            Vec tp = th, tm = th; tp[k] += h; tm[k] -= h;
            const Vec Bp = BATable(*C.T, C.sm, S, tp), Bm = BATable(*C.T, C.sm, S, tm);
            const double kp = tp[S.iK[1]], km = tm[S.iK[1]];
            St a{&D.N, &D.CB, &Bp, &G->S0, C.thF0, kp}, b{&D.N, &D.CB, &Bm, &G->S0, C.thF0, km};
            g[k] = (f(a) - f(b)) / (2 * h);
        }
        o.eA = std::sqrt(std::max(0., wp5::Quad(C.A.F.Cstat, g)));
    } else if (wp3CovStat) {
        const Vec& v = C.W3->v;
        Vec g(v.size(), 0.);
        for (size_t k = 0; k < v.size(); ++k) {
            if ((*wp3CovStat)[k][k] <= 0) continue;
            const double h = 0.05 * std::sqrt((*wp3CovStat)[k][k]);
            Vec tp = v, tm = v; tp[k] += h; tm[k] -= h;
            const Vec Bp = BATableWP3(*C.W3, tp), Bm = BATableWP3(*C.W3, tm);
            St a{&D.N, &D.CB, &Bp, &G->S0, C.thF0, tp[8]}, b{&D.N, &D.CB, &Bm, &G->S0, C.thF0, tm[8]};
            g[k] = (f(a) - f(b)) / (2 * h);
        }
        o.eA = std::sqrt(std::max(0., wp5::Quad(*wp3CovStat, g)));
        if (wp3CovMC) o.mc = std::sqrt(std::max(0., wp5::Quad(*wp3CovMC, g)));
    }
    // F-main (total derivative with the stage-A refit when available)
    {
        Vec g(4, 0.);
        for (int k = 0; k < 4; ++k) {
            const double h = C.hF.empty() ? 0.25 * std::sqrt(std::max(G->CF[k][k], 1e-12)) : C.hF[k];
            Vec tp = C.thF0, tm = C.thF0; tp[k] += h; tm[k] -= h;
            const Vec& Bp = C.hF.empty() ? C.BA0 : C.BAFp[k];
            const Vec& Bm = C.hF.empty() ? C.BA0 : C.BAFm[k];
            const double kp = C.hF.empty() ? C.kE2 : C.kE2p[k], km = C.hF.empty() ? C.kE2 : C.kE2m[k];
            St a{&D.N, &D.CB, &Bp, &G->S0, tp, kp}, b{&D.N, &D.CB, &Bm, &G->S0, tm, km};
            g[k] = (f(a) - f(b)) / (2 * h);
        }
        o.eF = std::sqrt(std::max(0., wp5::Quad(G->CF, g)));
    }
    o.stat = std::sqrt(o.eN * o.eN + o.eCB * o.eCB + o.eA * o.eA + o.eF * o.eF);
    if (withMC && !C.BAr.empty()) {
        double s1 = 0, s2 = 0; int n = 0;
        for (size_t q = 0; q < C.BAr.size(); ++q) {
            const int r = C.repIdx[q];
            St s{&D.N, &D.CB, &C.BAr[q], &G->Srep[r], C.thFr[q], C.kE2r[q]};
            const double v = f(s);
            o.rep.push_back(v); s1 += v; s2 += v * v; ++n;
        }
        if (n > 1) {
            const double m = s1 / n;
            o.mc = std::sqrt(std::max(0., (s2 - n * m * m) / (n - 1)));
            o.shift = m - o.v; o.bc = 2 * o.v - m; o.nrep = n;
            Vec srt = o.rep; std::sort(srt.begin(), srt.end());
            auto qt = [&](double x) { const double y = x * (n - 1); const int i = (int)y; const double w = y - i; return i + 1 < n ? (1 - w) * srt[i] + w * srt[i + 1] : srt[n - 1]; };
            o.mc68 = 0.5 * (qt(0.8413) - qt(0.1587));
        }
    }
    o.tot = std::sqrt(o.stat * o.stat + o.mc * o.mc);
    return o;
}
void PrintE(const ERes& e, const char* unit = "") {
    P("    %-46s %11.2f +- %8.2f stat (N %.2f, CB %.2f, stage A %.2f, F-main %.2f) +- %8.2f MC  [total %.2f; MC shift %+.2f, bias-corr. %.2f, mc68 %.2f] %s\n",
      e.name.c_str(), e.v, e.stat, e.eN, e.eCB, e.eA, e.eF, e.mc, e.tot, e.shift, e.bc, e.mc68, unit);
    RS("est %s %.6g stat %.6g N %.6g CB %.6g A %.6g F %.6g mc %.6g tot %.6g shift %.6g bc %.6g mc68 %.6g nrep %d\n", e.name.c_str(), e.v, e.stat, e.eN,
       e.eCB, e.eA, e.eF, e.mc, e.tot, e.shift, e.bc, e.mc68, e.nrep);
}

// ---------------------------------------------------------------------------------------------------------------
// the estimators
// ---------------------------------------------------------------------------------------------------------------
inline double SumC(const Vec& v, const std::vector<int>& cells, int w) { double s = 0; for (int c : cells) s += v[DI(c, w)]; return s; }
inline double SumS(const Vec& S, int comp, const std::vector<int>& cells, int w) { double s = 0; for (int c : cells) s += SW(S, comp, c, w); return s; }
inline double EtaAll(const Vec& S, const Vec& f, const std::vector<int>& cells, int w) {
    return f[0] * SumS(S, wC2, cells, w) + f[0] * f[1] * (SumS(S, wC3, cells, w) + SumS(S, wC4, cells, w)) + f[0] * f[2] * SumS(S, wC5, cells, w) +
           f[0] * f[3] * SumS(S, wGeta, cells, w);
}
inline double EtaOth(const Vec& S, const Vec& f, const std::vector<int>& cells, int w) {
    return f[0] * f[1] * SumS(S, wC4, cells, w) + f[0] * f[2] * SumS(S, wC5, cells, w) + f[0] * f[3] * SumS(S, wGeta, cells, w);
}
inline double Tail(const Vec& S, const Vec& f, const std::vector<int>& cells, int w) {   // bump tail + G_eta (in R')
    return f[0] * SumS(S, wC2, cells, w) + f[0] * f[1] * SumS(S, wC3, cells, w) + f[0] * f[3] * SumS(S, wGeta, cells, w);
}
inline double E45(const Vec& S, const Vec& f, const std::vector<int>& cells, int w) {
    return f[0] * f[1] * SumS(S, wC4, cells, w) + f[0] * f[2] * SumS(S, wC5, cells, w);
}
inline double U45(const Vec& S, const std::vector<int>& cells, int w) { return SumS(S, wC4, cells, w) + SumS(S, wC5, cells, w); }

// B4: subU = subtract the eta content of U; kcb = factor on the CB (1 nominal)
Est MakeB4(const std::vector<int>& cells, int wB, int wU, bool subU = true, bool kcbFromA = false) {
    return [=](const St& s) {
        const double k = kcbFromA ? s.kE2 : 1.;
        const double sB = SumC(*s.N, cells, wB) - k * SumC(*s.CB, cells, wB);
        const double sU = SumC(*s.N, cells, wU) - k * SumC(*s.CB, cells, wU);
        const double R = SumC(*s.BA, cells, wB) / SumC(*s.BA, cells, wU);
        const double eU = subU ? EtaAll(*s.S, s.thF, cells, wU) : 0.;
        return kYf * (sB - (sU - eU) * R - EtaOth(*s.S, s.thF, cells, wB));
    };
}
// B4 with the transfer done slice by slice and summed
Est MakeB4PerSlice(int e, int wB, int wU) {
    return [=](const St& s) {
        double y = 0;
        for (int j = 0; j < kNSlice; ++j) y += MakeB4(Cells(e, j), wB, wU)(s);
        return y;
    };
}
// fraction of the bump (C2 + kappa C3) in B that falls into the cells (simulation)
inline double FracE2(const Vec& S, const Vec& f, int wB) {
    const std::vector<int> e2 = Cells(1), all = Cells(2);
    return (SumS(S, wC2, e2, wB) + f[1] * SumS(S, wC3, e2, wB)) / (SumS(S, wC2, all, wB) + f[1] * SumS(S, wC3, all, wB));
}
Est MakeB4All(int wB, int wU) {
    return [=](const St& s) { return MakeB4(Cells(1), wB, wU)(s) / FracE2(*s.S, s.thF, wB); };
}
// peak side bands. mode 0: R' scales background + bump tail + G_eta (WP4 implementation); 1: only the background is
// normalised by the side bands, bump tail and G_eta fixed at the F-main values. fFix >= 0: leakage fraction given.
// leakW 1: leakage shape kappa C4 + rho C5 instead of C4 + C5. nIter > 0: return the n-th iterate instead of the limit.
Est MakePSB(const std::vector<int>& cells, int wP, int wL, int wH, int mode = 0, double fFix = -1, int leakW = 0, int nIter = 0) {
    return [=](const St& s) {
        const double sP = SumC(*s.N, cells, wP) - SumC(*s.CB, cells, wP);
        const double sS = SumC(*s.N, cells, wL) - SumC(*s.CB, cells, wL) + SumC(*s.N, cells, wH) - SumC(*s.CB, cells, wH);
        const double bP = SumC(*s.BA, cells, wP), bS = SumC(*s.BA, cells, wL) + SumC(*s.BA, cells, wH);
        const double tP = Tail(*s.S, s.thF, cells, wP), tS = Tail(*s.S, s.thF, cells, wL) + Tail(*s.S, s.thF, cells, wH);
        double f = fFix;
        if (f < 0) {
            if (leakW == 0) f = (U45(*s.S, cells, wL) + U45(*s.S, cells, wH)) / U45(*s.S, cells, wP);
            else f = (E45(*s.S, s.thF, cells, wL) + E45(*s.S, s.thF, cells, wH)) / E45(*s.S, s.thF, cells, wP);
        }
        double R, a, b;
        if (mode == 0) { R = (bP + tP) / (bS + tS); a = sP; b = sS; }
        else { R = bP / bS; a = sP - tP; b = sS - tS; }
        if (nIter > 0) {   // Y_0 = a - R b; Y_{n+1} = a - R (b - f Y_n)
            double y = a - R * b;
            for (int it = 1; it < nIter; ++it) y = a - R * (b - f * y);
            return kYf * y;
        }
        return kYf * (a - R * b) / (1. - R * f);
    };
}
// CR prediction: k_E2 CB + B_A + eta(F-main), and the observed N (E2 cells)
Est MakePred(const std::vector<int>& cells, int w) {
    return [=](const St& s) { return s.kE2 * SumC(*s.CB, cells, w) + SumC(*s.BA, cells, w) + EtaAll(*s.S, s.thF, cells, w); };
}

// ---------------------------------------------------------------------------------------------------------------
// figures
// ---------------------------------------------------------------------------------------------------------------
void Style() {
    bla::Style();
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadBottomMargin(0.14);
}
void Save(TCanvas* c, const std::string& base) {
    c->SaveAs((base + ".pdf").c_str());
    c->SaveAs((base + ".png").c_str());
}
TLatex* Txt(double x, double y, const char* s, double size = 0.045, int col = kBlack, int align = 11) {
    TLatex* t = new TLatex(x, y, s);
    t->SetNDC(); t->SetTextSize(size); t->SetTextColor(col); t->SetTextAlign(align); t->SetTextFont(42);
    t->Draw();
    return t;
}

}  // namespace wp6

void wp6_ana_cnt() {
    using namespace alldefs;
    using namespace wp6;
    TStopwatch total;
    gSystem->mkdir(kOut, kTRUE);
    gLog = fopen(Form("%s/wp6_ana.log", kCntOut), "w");
    gRes = fopen(Form("%s/results.txt", kCntOut), "w");
    RS("# =====================================================================================================\n");
    RS("# WP6 side bands and cross-checks (SPEC report_all) -- plots/all_report/wp6_ana.C (+ wp6_fill.C), 2026-11-14\n");
    RS("# S0, isBest >= 0, raw data rows (yields x 64/63), L = %.1f nb^-1; simulation SMASH mode A, w x V3, M4' = %.5f M4\n", kL, kS);
    RS("# est <name> value stat N CB A(stage A) F(F-main) mc tot shift(replica mean - value) bc(2v - mean) mc68 nrep\n");
    RS("# =====================================================================================================\n");
    G = new Glob();
    Glob& g = *G;

    // ---------------------------------------------------------------- inputs
    {
        TFile* ff = TFile::Open(Form("%s/wp6_fill.root", kOut));
        TNamed* wn = (TNamed*)ff->Get("sim/windows");
        std::stringstream ss(wn->GetTitle());
        std::string item;
        while (std::getline(ss, item, ';')) {
            if (item.empty()) continue;
            std::istringstream is(item); Win w; is >> w.tag >> w.lo >> w.hi; gW.push_back(w);
        }
        nW = (int)gW.size();
        TH1* hs = (TH1*)ff->Get("sim/win_sum");
        TH2* hr = (TH2*)ff->Get("sim/win_rep");
        const int n = hs->GetNbinsX();
        g.S0.assign(n, 0.); g.Svar.assign(n, 0.);
        for (int i = 0; i < n; ++i) { g.S0[i] = hs->GetBinContent(i + 1); g.Svar[i] = std::pow(hs->GetBinError(i + 1), 2); }
        g.Srep.assign(kNRep, Vec(n, 0.));
        for (int r = 0; r < kNRep; ++r) for (int i = 0; i < n; ++i) g.Srep[r][i] = hr->GetBinContent(i + 1, r + 1);
        ff->Close();
        P("wp6_ana: %d windows, %d simulated window sums, %d replicas\n", nW, n, kNRep);
    }
    TFile* fh = TFile::Open(kHistFile);
    if (!wp5::LoadTempl(fh, "S0", "nom", 1.0, g.T0)) return;
    if (!wp5::LoadReplicas(fh, "S0", 1.0, g.T0, kNRep, g.R0)) return;
    // the WP3 smoothing (nominal of this package; the WP4 smoothing of T0/R0 is a variant)
    {
        double chi[20];
        g.T3 = g.T0;
        SmoothTemplW3(g.T3, g.sm3, nullptr, chi);
        g.R3 = g.R0;
        g.sm3R.resize(kNRep);
        for (int r = 0; r < kNRep; ++r) SmoothTemplW3(g.R3[r], g.sm3R[r], &g.sm3);
        P("WP3 smoothing of G_D, G_C (F10): start of the domain and chi2/ndf (mean weight):\n  ");
        for (int h = 0; h < 10; ++h) P(" %s %.4f/%.4f (%.2f, %.2f)", F10Name(h / 2, h % 2).c_str(), g.sm3.start[h][0], g.sm3.start[h][1], chi[2 * h], chi[2 * h + 1]);
        P("\n");
    }
    if (!wp5::LoadData(fh, "S0", "all", g.D0, gLog)) return;
    if (!LoadDataW(g.DW)) return;
    P("data rows: e+e- %ld, like-sign %ld (S0)\n", g.DW.nUS, g.DW.nLS);
    const int iB = WI("B"), iP = WI("P"), iSl = WI("SBlo"), iSh = WI("SBhi"), iU = WI("U575"), iU565 = WI("U565"), iL = WI("CRL"),
              iK = WI("K400"), iR0 = WI("R0"), iSlN = WI("SBloN"), iShN = WI("SBhiN"), iSlW = WI("SBloW"), iShW = WI("SBhiW"), iTot = WI("tot");
    const std::vector<int> cE1 = Cells(0), cE2 = Cells(1), cAll = Cells(2);

    // ---------------------------------------------------------------- 0. checks of the inputs
    P("\n0. CHECKS OF THE INPUTS\n");
    {
        TH1* nw = allio::Data(fh, "S0", "all", "N", "win");
        TH1* cw = allio::Data(fh, "S0", "all", "CB", "win");
        TH2* cr = (TH2*)fh->Get("data/S0/all/CBrep_win");
        const int wid[6] = {kWB, kWP, kWSBlo, kWSBhi, kWCRL, kWCRU};
        const int mine[6] = {iB, iP, iSl, iSh, iL, iU565};
        for (int q = 0; q < 6; ++q) {
            const double n = SumC(g.DW.N, cAll, mine[q]), c = SumC(g.DW.CB, cAll, mine[q]);
            double vv = 0;
            const int nI = kNCell * nW;
            for (int a : cAll) for (int b : cAll) vv += g.DW.covCB[(size_t)DI(a, mine[q]) * nI + DI(b, mine[q])];
            double m = 0, m2 = 0;
            if (cr) { for (int r = 1; r <= cr->GetNbinsY(); ++r) { const double x = cr->GetBinContent(wid[q] + 1, r); m += x; m2 += x * x; } m /= cr->GetNbinsY(); m2 = std::sqrt(std::max(0., m2 / cr->GetNbinsY() - m * m)); }
            P("   window %-5s data N %8.0f (hists_v1 %8.0f), CB %10.3f (hists_v1 %10.3f, diff %.2e); CB replica RMS %.2f (hists_v1 %.2f)\n", gW[mine[q]].tag.c_str(), n,
              nw->GetBinContent(wid[q] + 1), c, cw->GetBinContent(wid[q] + 1), c - cw->GetBinContent(wid[q] + 1), std::sqrt(vv), m2);
        }
        // E2 sums of B and U575 against the F10 histograms
        const std::vector<double> e = F10Edges();
        double nB = 0, cbB = 0, nU = 0, cbU = 0, sB = 0, sBr0 = 0;
        for (int j = 0; j < kNSlice; ++j) {
            const int h = 2 * j + 1;
            for (int x = 1; x <= kNF10; ++x) {
                const bool inB = e[x - 1] >= kBlo - 1e-9 && e[x] <= kBhi + 1e-9, inU = e[x - 1] >= 0.575 - 1e-9 && e[x] <= 0.700 + 1e-9;
                if (inB) { nB += g.D0.y[h][x]; cbB += g.D0.cb[h][x]; sB += g.T0.s[wp5::kS2][h][x]; sBr0 += g.R0[0].s[wp5::kS2][h][x]; }
                if (inU) { nU += g.D0.y[h][x]; cbU += g.D0.cb[h][x]; }
            }
        }
        P("   E2: N(B) %.0f vs F10 %.0f; CB(B) %.4f vs %.4f; N(U575) %.0f vs %.0f; CB(U575) %.4f vs %.4f\n", SumC(g.DW.N, cE2, iB), nB,
          SumC(g.DW.CB, cE2, iB), cbB, SumC(g.DW.N, cE2, iU), nU, SumC(g.DW.CB, cE2, iU), cbU);
        P("   sim C2 in B, E2: window sums %.4f vs F10 templates %.4f; replica 0: %.4f vs %.4f (same replica generator)\n", SumS(g.S0, wC2, cE2, iB), sB,
          SumS(g.Srep[0], wC2, cE2, iB), sBr0);
        TH1* c2w = allio::Sim(fh, "S0", "nom", kC2, "win");
        TH1* g5w = allio::Sim(fh, "S0", "nom", kC5, "win");
        P("   sim C2 in P (all cells) %.4f vs hists_v1 win %.4f; C5 in P %.4f vs %.4f\n", SumS(g.S0, wC2, cAll, iP), c2w->GetBinContent(kWP + 1) * kNbToRaw,
          SumS(g.S0, wC5, cAll, iP), g5w->GetBinContent(kWP + 1) * kNbToRaw);
        RS("check CB_P %.6f %.6f\n", SumC(g.DW.CB, cAll, iP), cw->GetBinContent(kWP + 1));
        delete nw; delete cw; delete c2w; delete g5w;
    }

    // ---------------------------------------------------------------- 1. F-main (WP5 nominal) and its replicas
    P("\n1. F-MAIN (WP5 nominal: F10 E2, kappa = 1 +- 0.2, g = 1 +- 0.3, k_E2 = 1 +- 0.05) refitted with wp5_lib.h\n");
    Opt oNom; oNom.useE = 2; oNom.kappaCons = true;
    g.SM = wp5::SpecMain(oNom);
    Bins BM; wp5::MakeBins(g.D0, g.T0, 0, 2, BM);
    g.FM = wp5::Fit(BM, g.SM, g.D0);
    if (!g.FM.ok) { P("F-main failed\n"); return; }
    g.RepFM = wp5::RunReplicas(BM, g.SM, g.D0, g.R0, g.FM.th);
    g.thF0 = ThF(g.FM.th);
    g.CF = wp5::Zero(4, 4); g.CFmc = wp5::Zero(4, 4);
    for (int a = 0; a < 4; ++a) for (int b = 0; b < 4; ++b) { g.CF[a][b] = g.FM.Cstat[a][b]; g.CFmc[a][b] = g.RepFM.cov[a][b]; }
    {
        // compare with WP5 results.txt
        std::ifstream in(kWP5);
        std::string line; Vec v5;
        while (std::getline(in, line)) if (line.rfind("values FMAIN ", 0) == 0) { std::istringstream is(line.substr(13)); double x; while (is >> x) v5.push_back(x); break; }
        double md = 0;
        for (size_t k = 0; k < v5.size() && k < g.FM.th.size(); ++k) md = std::max(md, std::fabs(v5[k] - g.FM.th[k]));
        P("   a_eta %.5f +- %.5f (stat) +- %.5f (MC), kappa %.4f, rho_d %.4f, g %.4f; replicas %d of %zu; max |this - WP5| over 12 parameters = %.2e\n",
          g.thF0[0], std::sqrt(g.CF[0][0]), std::sqrt(g.CFmc[0][0]), g.thF0[1], g.thF0[2], g.thF0[3], g.RepFM.nok, g.R0.size(), md);
        RS("fmain a_eta %.6f kappa %.6f rho_d %.6f g %.6f stat %.6f %.6f %.6f %.6f mc %.6f %.6f %.6f %.6f maxdiff_wp5 %.3g\n", g.thF0[0], g.thF0[1], g.thF0[2], g.thF0[3],
           std::sqrt(g.CF[0][0]), std::sqrt(g.CF[1][1]), std::sqrt(g.CF[2][2]), std::sqrt(g.CF[3][3]), std::sqrt(g.CFmc[0][0]), std::sqrt(g.CFmc[1][1]),
           std::sqrt(g.CFmc[2][2]), std::sqrt(g.CFmc[3][3]), md);
    }

    // ---------------------------------------------------------------- 2. stage A
    P("\n2. STAGE A re-run at the WP5 eta values (WP3 frozen configuration: M1, E2, CR-U 0.575-0.700, no CR-K, WP3 smoothing; wp5_lib fitter)\n");
    Ctx C0 = BuildCtx("nominal", 0.575, 0.700, nullptr, false, true, true);
    if (!C0.ok) return;
    WP3Stage W3;
    const bool haveW3 = LoadWP3(W3);
    {
        const FitRes& F = C0.A.F;
        P("   %-8s", "");
        for (int k = 0; k < C0.A.S.np(); ++k) P(" %10s", C0.A.S.names[k].c_str());
        P("\n   %-8s", "this");
        for (int k = 0; k < C0.A.S.np(); ++k) P(" %10.4f", F.th[k]);
        P("\n   %-8s", "stat");
        for (int k = 0; k < C0.A.S.np(); ++k) P(" %10.4f", std::sqrt(std::max(0., F.Cstat[k][k])));
        // replica spread of the stage-A parameters (re-fit in the chain)
        P("\n   chi2 Pearson %.1f / %d (WP3 data-only Pearson 326.4/118); %zu replica chains\n", F.chi2P, F.ndf, C0.BAr.size());
        // the same configuration with the WP3 eta values (a 0.4, kappa = rho = g = 1): comparison with the frozen WP3 numbers
        const Vec f04 = {0.4, 1., 1., 1.};
        const StageA a04 = FitStageA(g.D0, g.T3, f04, 0.575, 0.700, nullptr, false, nullptr, false);
        P("   chi2 Pearson with the WP3 eta: %.1f / %d\n", a04.F.chi2P, a04.F.ndf);
        P("   same fit with the WP3 eta (a = 0.4, kappa = rho = g = 1):");
        for (int k = 0; k < a04.S.np(); ++k) P(" %s %.4f", a04.S.names[k].c_str(), a04.F.th[k]);
        P("\n   WP3 frozen (stageA.txt, its own smoothing):              c_1..c_5 %.4f %.4f %.4f %.4f %.4f b %.4f kappa_B %.3f k_E2 %.4f\n", W3.v[0], W3.v[1],
          W3.v[2], W3.v[3], W3.v[4], W3.v[5], W3.v[6], W3.v[8]);
        RS("stageA names");
        for (int k = 0; k < C0.A.S.np(); ++k) RS(" %s", C0.A.S.names[k].c_str());
        RS("\nstageA values");
        for (int k = 0; k < C0.A.S.np(); ++k) RS(" %.6g", F.th[k]);
        RS("\nstageA stat");
        for (int k = 0; k < C0.A.S.np(); ++k) RS(" %.6g", std::sqrt(std::max(0., F.Cstat[k][k])));
        RS("\nstageA chi2P %.3f ndf %d\nstageA_aeta04 ", F.chi2P, F.ndf);
        for (int k = 0; k < a04.S.np(); ++k) RS(" %.6g", a04.F.th[k]);
        RS("\n");
        // GLS chi2 with the MC covariance of the templates (as WP3)
        const wp5::GLSRes gl = wp5::GLS(C0.A.B, C0.A.S, g.D0, F, g.R3);
        if (gl.ok) { P("   GLS chi2 with the template MC covariance: %.1f / %d = %.2f (WP3: 129.8/118 = 1.10)\n", gl.chi2, gl.ndf, gl.chi2 / gl.ndf); RS("stageA_gls %.3f %d\n", gl.chi2, gl.ndf); }
        // background in B and U: exact integral vs the fitted bin model
        const std::vector<double> e = F10Edges();
        double bB = 0, bU = 0;
        const Spec& S = C0.A.S;
        for (int j = 0; j < kNSlice; ++j) {
            const int h = 2 * j + 1;
            for (int x = 1; x <= kNF10; ++x) {
                const double xc = TiltX(0.5 * (e[x - 1] + e[x]));
                const double v = F.th[S.iC[j]] * (1. + F.th[S.iB] * xc) * (g.T3.gdI[h][x] + F.th[S.iKB] * g.T3.gcI[h][x]);
                if (e[x - 1] >= kBlo - 1e-9 && e[x] <= kBhi + 1e-9) bB += v;
                if (e[x - 1] >= 0.575 - 1e-9 && e[x] <= 0.700 + 1e-9) bU += v;
            }
        }
        P("   B_A(B, E2) = %.1f (bin model %.1f), B_A(U, E2) = %.1f (bin model %.1f)\n", SumC(C0.BA0, cE2, iB), bB, SumC(C0.BA0, cE2, iU), bU);
        if (haveW3) {
            double wB = 0, wU = 0;
            for (int j = 0; j < kNSlice; ++j)
                for (int x = 1; x <= kNF10; ++x) {
                    if (e[x - 1] >= kBlo - 1e-9 && e[x] <= kBhi + 1e-9) wB += W3.binBA[2 * j + 1][x];
                    if (e[x - 1] >= 0.575 - 1e-9 && e[x] <= 0.700 + 1e-9) wU += W3.binBA[2 * j + 1][x];
                }
            const Vec T3 = BATableWP3(W3, W3.v);
            P("   WP3 frozen: B_A(B, E2) = %.1f (stageA.txt bins %.1f), B_A(U, E2) = %.1f (bins %.1f)\n", SumC(T3, cE2, iB), wB, SumC(T3, cE2, iU), wU);
        }
    }

    // ---------------------------------------------------------------- 3. estimators, nominal
    P("\n3. SIDE-BAND ESTIMATORS (nominal stage A of section 2; errors: stat = N + CB + stage A + F-main, MC = 100 replica chains)\n");
    std::vector<ERes> E;
    auto add = [&](const std::string& n, const Est& f, const Ctx& C, bool mc = true) { E.push_back(Evaluate(n, f, C, mc)); PrintE(E.back()); return E.back(); };
    const ERes b4 = add("B4_Ybump_E2", MakeB4(cE2, iB, iU), C0);
    const ERes b4all = add("B4_Ybump_all_extrap", MakeB4All(iB, iU), C0);
    const ERes psb = add("PSB_Ypeak", MakePSB(cAll, iP, iSl, iSh), C0);
    const ERes psbE1 = add("PSB_Ypeak_E1", MakePSB(cE1, iP, iSl, iSh), C0);
    const ERes psbE2 = add("PSB_Ypeak_E2", MakePSB(cE2, iP, iSl, iSh), C0);
    // ingredients
    const St s0{&g.DW.N, &g.DW.CB, &C0.BA0, &g.S0, g.thF0, C0.kE2};
    const double fE2 = FracE2(g.S0, g.thF0, iB);
    {
        const double sB = SumC(g.DW.N, cE2, iB) - SumC(g.DW.CB, cE2, iB), sU = SumC(g.DW.N, cE2, iU) - SumC(g.DW.CB, cE2, iU);
        const double R = SumC(C0.BA0, cE2, iB) / SumC(C0.BA0, cE2, iU);
        const double eU = EtaAll(g.S0, g.thF0, cE2, iU), eB = EtaOth(g.S0, g.thF0, cE2, iB);
        P("   B4 ingredients (E2): N_B %.0f, CB_B %.1f, N_U %.0f, CB_U %.1f, eta_U %.1f (G_eta %.1f), B_A(B)/B_A(U) = %.5f, other eta in B %.1f (C4 %.1f, C5 %.1f, G_eta %.1f)\n",
          SumC(g.DW.N, cE2, iB), SumC(g.DW.CB, cE2, iB), SumC(g.DW.N, cE2, iU), SumC(g.DW.CB, cE2, iU), eU, g.thF0[0] * g.thF0[3] * SumS(g.S0, wGeta, cE2, iU), R, eB,
          g.thF0[0] * g.thF0[1] * SumS(g.S0, wC4, cE2, iB), g.thF0[0] * g.thF0[2] * SumS(g.S0, wC5, cE2, iB), g.thF0[0] * g.thF0[3] * SumS(g.S0, wGeta, cE2, iB));
        P("   background in B used by B4 = (N - CB - eta)_U x R = %.1f rows; f_E2 (bump in B with M(ee) >= 10 MeV, sim, kappa %.3f) = %.4f\n", (sU - eU) * R,
          g.thF0[1], fE2);
        RS("b4parts NB %.1f CBB %.4f NU %.1f CBU %.4f etaU %.4f R %.6f etaB %.4f bkgB %.4f fE2 %.6f\n", SumC(g.DW.N, cE2, iB), SumC(g.DW.CB, cE2, iB), SumC(g.DW.N, cE2, iU),
           SumC(g.DW.CB, cE2, iU), eU, R, eB, (sU - eU) * R, fE2);
        // simulation at SMASH (kappa = rho = 1) and at the F-main values, for the eta-scale comparison
        const double smB = kYf * (SumS(g.S0, wC2, cE2, iB) + SumS(g.S0, wC3, cE2, iB));
        const double smBall = kYf * (SumS(g.S0, wC2, cAll, iB) + SumS(g.S0, wC3, cAll, iB));
        const double smP = kYf * U45(g.S0, cAll, iP), smPE2 = kYf * U45(g.S0, cE2, iP);
        P("   simulation at SMASH x PDG BR (x 64/63): bump in B: E2 %.1f, all %.1f; peak C4 + C5 in P: all %.1f, E2 %.1f\n", smB, smBall, smP, smPE2);
        RS("simsmash bumpB_E2 %.4f bumpB_all %.4f peakP_all %.4f peakP_E2 %.4f C2B_E2 %.4f C3B_E2 %.4f C4P %.4f C5P %.4f\n", smB, smBall, smP, smPE2, kYf * SumS(g.S0, wC2, cE2, iB),
           kYf * SumS(g.S0, wC3, cE2, iB), kYf * SumS(g.S0, wC4, cAll, iP), kYf * SumS(g.S0, wC5, cAll, iP));
    }
    // side bands: the ingredients and the iteration
    {
        const double sP = SumC(g.DW.N, cAll, iP) - SumC(g.DW.CB, cAll, iP);
        const double sS = SumC(g.DW.N, cAll, iSl) - SumC(g.DW.CB, cAll, iSl) + SumC(g.DW.N, cAll, iSh) - SumC(g.DW.CB, cAll, iSh);
        const double bP = SumC(C0.BA0, cAll, iP), bS = SumC(C0.BA0, cAll, iSl) + SumC(C0.BA0, cAll, iSh);
        const double tP = Tail(g.S0, g.thF0, cAll, iP), tS = Tail(g.S0, g.thF0, cAll, iSl) + Tail(g.S0, g.thF0, cAll, iSh);
        const double f = (U45(g.S0, cAll, iSl) + U45(g.S0, cAll, iSh)) / U45(g.S0, cAll, iP);
        const double R = (bP + tP) / (bS + tS);
        P("   side bands (all M(ee)): N_P %.0f, CB_P %.1f, N_SB %.0f, CB_SB %.1f; B_A(P) %.1f, B_A(SB) %.1f; bump tail + G_eta: P %.1f, SB %.1f\n",
          SumC(g.DW.N, cAll, iP), SumC(g.DW.CB, cAll, iP), SumC(g.DW.N, cAll, iSl) + SumC(g.DW.N, cAll, iSh), SumC(g.DW.CB, cAll, iSl) + SumC(g.DW.CB, cAll, iSh), bP, bS, tP, tS);
        P("   R' = %.5f (background alone %.5f; flat background: width ratio %.5f); leakage f = (C4+C5)(SB)/(C4+C5)(P) = %.5f\n", R, bP / bS,
          (kPhi - kPlo) / (0.5336 - 0.5250 + 0.5628 - 0.5542), f);
        P("   iteration Y_{n+1} = (N - CB)_P - R' [(N - CB)_SB - f Y_n], raw rows:");
        double y = sP - R * sS;
        for (int it = 0; it < 6; ++it) { P(" %.1f", y); y = sP - R * (sS - f * y); }
        P(" -> limit (N - CB)_P - R' (N - CB)_SB)/(1 - R' f) = %.1f; x 64/63 = %.1f\n", (sP - R * sS) / (1 - R * f), kYf * (sP - R * sS) / (1 - R * f));
        RS("psbparts NP %.1f CBP %.4f NSB %.1f CBSB %.4f BAP %.4f BASB %.4f tailP %.4f tailSB %.4f Rp %.6f f %.6f\n", SumC(g.DW.N, cAll, iP), SumC(g.DW.CB, cAll, iP),
           SumC(g.DW.N, cAll, iSl) + SumC(g.DW.N, cAll, iSh), SumC(g.DW.CB, cAll, iSl) + SumC(g.DW.CB, cAll, iSh), bP, bS, tP, tS, R, f);
    }
    // derived: visible cross sections, eta scales, double ratios
    const double smBE2 = kYf * (SumS(g.S0, wC2, cE2, iB) + SumS(g.S0, wC3, cE2, iB));
    const double smP = kYf * U45(g.S0, cAll, iP);
    const ERes aB4 = add("a_bump_B4_E2 (Y_B4,E2 / SMASH)", [&](const St& s) { return MakeB4(cE2, iB, iU)(s) / (kYf * (SumS(*s.S, wC2, cE2, iB) + SumS(*s.S, wC3, cE2, iB))); }, C0);
    const ERes aPSB = add("a_peak_PSB (Y_PSB / SMASH C4+C5)", [&](const St& s) { return MakePSB(cAll, iP, iSl, iSh)(s) / (kYf * U45(*s.S, cAll, iP)); }, C0);
    const ERes DRall = add("DR_sideband (all M(ee), bump extrapolated)", [&](const St& s) {
        const double yb = MakeB4All(iB, iU)(s), yp = MakePSB(cAll, iP, iSl, iSh)(s);
        const double sb = SumS(*s.S, wC2, Cells(2), iB) + SumS(*s.S, wC3, Cells(2), iB), sp = U45(*s.S, Cells(2), iP);
        return (yp / yb) / (sp / sb); }, C0);
    const ERes DRE2 = add("DR_sideband_E2 (both in M(ee) >= 10 MeV)", [&](const St& s) {
        const double yb = MakeB4(Cells(1), iB, iU)(s), yp = MakePSB(Cells(1), iP, iSl, iSh)(s);
        const double sb = SumS(*s.S, wC2, Cells(1), iB) + SumS(*s.S, wC3, Cells(1), iB), sp = U45(*s.S, Cells(1), iP);
        return (yp / yb) / (sp / sb); }, C0);
    const ERes shE2 = add("peak share M(ee) >= 10 MeV (PSB_E2 / PSB)", [&](const St& s) { return MakePSB(Cells(1), iP, iSl, iSh)(s) / MakePSB(Cells(2), iP, iSl, iSh)(s); }, C0);
    P("   visible cross sections: bump (B4, all M(ee) extrapolated) %.3f +- %.3f (stat) +- %.3f (MC) nb; peak (side bands) %.4f +- %.4f +- %.4f nb [L = %.1f nb^-1]\n",
      b4all.v / kL, b4all.stat / kL, b4all.mc / kL, psb.v / kL, psb.stat / kL, psb.mc / kL, kL);
    RS("sigvis bump %.6g %.6g %.6g peak %.6g %.6g %.6g\n", b4all.v / kL, b4all.stat / kL, b4all.mc / kL, psb.v / kL, psb.stat / kL, psb.mc / kL);
    // extrapolation error of the all-M(ee) bump from kappa (F-main total error)
    {
        const double k0 = g.thF0[1], ek = std::sqrt(g.CF[1][1] + g.CFmc[1][1]);
        Vec fp = g.thF0, fm = g.thF0; fp[1] = k0 + ek; fm[1] = std::max(0., k0 - ek);
        const double y0 = b4.v / FracE2(g.S0, g.thF0, iB), yp = b4.v / FracE2(g.S0, fp, iB), ym = b4.v / FracE2(g.S0, fm, iB);
        P("   extrapolation E2 -> all M(ee): kappa = %.3f +- %.3f gives f_E2 %.4f (%.4f / %.4f): Y_bump(all) %.0f, %+.0f / %+.0f\n", k0, ek, FracE2(g.S0, g.thF0, iB),
          FracE2(g.S0, fp, iB), FracE2(g.S0, fm, iB), y0, yp - y0, ym - y0);
        RS("extrap kappa %.5f %.5f fE2 %.6f %.6f %.6f Yall %.4f dYp %.4f dYm %.4f\n", k0, ek, FracE2(g.S0, g.thF0, iB), FracE2(g.S0, fp, iB), FracE2(g.S0, fm, iB), y0, yp - y0, ym - y0);
    }

    // ---------------------------------------------------------------- 4. variants
    P("\n4. VARIANTS (stat and MC as above unless stated)\n");
    std::vector<ERes> V;
    auto addV = [&](const std::string& n, const Est& f, const Ctx& C, bool mc = true) { V.push_back(Evaluate(n, f, C, mc)); PrintE(V.back()); return V.back(); };
    // 4a. side-band definitions
    addV("PSB per class (E1 + E2 separately)", [&](const St& s) { return MakePSB(cE1, iP, iSl, iSh)(s) + MakePSB(cE2, iP, iSl, iSh)(s); }, C0);
    addV("PSB background-only normalisation", MakePSB(cAll, iP, iSl, iSh, 1), C0);
    addV("PSB leakage shape kappa C4 + rho C5", MakePSB(cAll, iP, iSl, iSh, 0, -1, 1), C0);
    addV("PSB no leakage correction (f = 0)", MakePSB(cAll, iP, iSl, iSh, 0, 0.0), C0);
    addV("PSB narrow side bands (5.6 MeV)", MakePSB(cAll, iP, iSlN, iShN), C0);
    addV("PSB wide side bands (11.6 MeV)", MakePSB(cAll, iP, iSlW, iShW), C0);
    addV("PSB flat background (R' = width ratio)", [&](const St& s) {
        const double sP = SumC(*s.N, cAll, iP) - SumC(*s.CB, cAll, iP);
        const double sS = SumC(*s.N, cAll, iSl) - SumC(*s.CB, cAll, iSl) + SumC(*s.N, cAll, iSh) - SumC(*s.CB, cAll, iSh);
        const double R = (kPhi - kPlo) / (0.5336 - 0.5250 + 0.5628 - 0.5542);
        const double f = (U45(*s.S, cAll, iSl) + U45(*s.S, cAll, iSh)) / U45(*s.S, cAll, iP);
        return kYf * (sP - R * sS) / (1 - R * f); }, C0, false);
    // 4b. B4 definitions
    addV("B4 E2 without the eta_U subtraction (spec literal)", MakeB4(cE2, iB, iU, false), C0);
    addV("B4 E2 with CB x k_E2 (stage A)", MakeB4(cE2, iB, iU, true, true), C0);
    addV("B4 E2 per-slice transfer", MakeB4PerSlice(1, iB, iU), C0);
    addV("B4 E2 normalised on U 0.565-0.700 (same stage A)", MakeB4(cE2, iB, iU565), C0);
    addV("B4 E1 + E2 direct (E1 shape predicted)", MakeB4(cAll, iB, iU), C0);
    for (int j = 0; j < kNSlice; ++j) addV(Form("B4 E2 slice S%d alone", j + 1), MakeB4(Cells(1, j), iB, iU), C0);
    // 4c. stage-A variants
    Ctx Cq, C4s;
    {
        Cq = BuildCtx("quad", 0.575, 0.700, nullptr, true, true, true);
        if (Cq.ok) {
            P("   stage A + quadratic term: q = %.4f +- %.4f\n", Cq.A.F.th[Cq.A.S.iQ], std::sqrt(Cq.A.F.Cstat[Cq.A.S.iQ][Cq.A.S.iQ]));
            addV("B4 E2, stage A + quadratic", MakeB4(cE2, iB, iU), Cq);
            addV("PSB, stage A + quadratic", MakePSB(cAll, iP, iSl, iSh), Cq);
        }
        C4s = BuildCtx("WP4smooth", 0.575, 0.700, nullptr, false, true, true, 0);
        if (C4s.ok) {
            P("   stage A with the WP4 smoothing (as WP5 F-main): c_1..c_5 %.4f %.4f %.4f %.4f %.4f b %.4f kappa_B %.3f k_E2 %.4f; chi2 Pearson %.1f/%d\n", C4s.A.F.th[0],
              C4s.A.F.th[1], C4s.A.F.th[2], C4s.A.F.th[3], C4s.A.F.th[4], C4s.A.F.th[5], C4s.A.F.th[6], C4s.A.F.th[7], C4s.A.F.chi2P, C4s.A.F.ndf);
            const wp5::GLSRes gl = wp5::GLS(C4s.A.B, C4s.A.S, g.D0, C4s.A.F, g.R0);
            if (gl.ok) P("   ... GLS chi2 with the template MC covariance %.1f / %d = %.2f\n", gl.chi2, gl.ndf, gl.chi2 / gl.ndf);
            addV("B4 E2, stage A with the WP4 smoothing", MakeB4(cE2, iB, iU), C4s);
            addV("PSB, stage A with the WP4 smoothing", MakePSB(cAll, iP, iSl, iSh), C4s);
        }
        Ctx C565 = BuildCtx("CRU565", 0.565, 0.700, nullptr, false, true, true);
        if (C565.ok) {
            addV("B4 E2, stage A and U on 0.565-0.700", MakeB4(cE2, iB, iU565), C565);
            addV("PSB, stage A on 0.565-0.700", MakePSB(cAll, iP, iSl, iSh), C565);
        }
        if (haveW3) {
            // WP3 frozen stage A: stat and MC from its parameter covariances (cov_stat, cov_mc of stageA.txt)
            wp5::Mat cs = wp5::Zero(9, 9), cm = wp5::Zero(9, 9);
            std::ifstream in(kStageA);
            std::string line; int mode = 0, row = 0;
            while (std::getline(in, line)) {
                if (line.rfind("cov_stat", 0) == 0) { mode = 1; row = 0; continue; }
                if (line.rfind("cov_mc", 0) == 0) { mode = 2; row = 0; continue; }
                if (line.rfind("cov_total", 0) == 0) { mode = 0; continue; }
                if (mode > 0 && row < 9) { std::istringstream is(line); for (int c = 0; c < 9; ++c) is >> (mode == 1 ? cs : cm)[row][c]; ++row; }
            }
            Ctx C3 = BuildCtxWP3(W3);
            V.push_back(Evaluate("B4 E2, WP3 frozen stage A (eta 0.4, g = 1)", MakeB4(cE2, iB, iU), C3, false, &cs, &cm)); PrintE(V.back());
            V.push_back(Evaluate("PSB, WP3 frozen stage A", MakePSB(cAll, iP, iSl, iSh), C3, false, &cs, &cm)); PrintE(V.back());
        }
    }

    // ---------------------------------------------------------------- 5. closure
    P("\n5. CLOSURE\n");
    // 5a. displaced windows
    const double fNom = (U45(g.S0, cAll, iSl) + U45(g.S0, cAll, iSh)) / U45(g.S0, cAll, iP);
    const double fNomE2 = (U45(g.S0, cE2, iSl) + U45(g.S0, cE2, iSh)) / U45(g.S0, cE2, iP);
    struct DispRes { double c, y, s, sMC, exp, pull, yNR, sNR, pullNR, yE2, sE2, pullE2; };
    std::vector<DispRes> disp;
    P("   5a. displaced windows: the side-band estimator (all M(ee), leakage f of the nominal peak = %.4f) at 13 positions; the window and its side bands\n", fNom);
    P("       are removed from CR-U before the stage-A refit. pull = (Y - Y_expected)/sigma_stat; Y_expected = the simulated C4 + C5 content (F-main scale).\n");
    P("       %-7s | %9s %8s %8s %7s | %9s %8s %7s (no refit) | %9s %8s %7s (E2 only)\n", "centre", "Y", "stat", "Y_exp", "pull", "Y", "stat", "pull", "Y", "stat", "pull");
    double sp2 = 0, spNR2 = 0, spE2 = 0, spm = 0;
    for (int k = 0; k < kNDisp; ++k) {
        const double c = DispCentre(k);
        const int wPk = WI(Form("P_d%d", k)), wLk = WI(Form("SBlo_d%d", k)), wHk = WI(Form("SBhi_d%d", k));
        auto drop = [c](double lo, double hi) { return hi > c + kDSl1 + 1e-9 && lo < c + kDSh2 - 1e-9; };
        Ctx Ck = BuildCtx(Form("disp%d", k), 0.575, 0.700, drop, false, false, false);
        DispRes d; d.c = c;
        const Est est = MakePSB(cAll, wPk, wLk, wHk, 0, fNom);
        const ERes r = Evaluate(Form("disp%d", k), est, Ck, false);
        const ERes rNR = Evaluate(Form("disp%d_norefit", k), est, C0, false);
        const ERes rE2 = Evaluate(Form("disp%d_E2", k), MakePSB(cE2, wPk, wLk, wHk, 0, fNomE2), Ck, false);
        // expected: the simulated C4 + C5 content through the same estimator
        {
            const St sk{&g.DW.N, &g.DW.CB, &Ck.BA0, &g.S0, g.thF0, Ck.kE2};
            const double bP = SumC(Ck.BA0, cAll, wPk) + Tail(g.S0, g.thF0, cAll, wPk);
            const double bS = SumC(Ck.BA0, cAll, wLk) + SumC(Ck.BA0, cAll, wHk) + Tail(g.S0, g.thF0, cAll, wLk) + Tail(g.S0, g.thF0, cAll, wHk);
            const double R = bP / bS;
            const double e45P = E45(g.S0, g.thF0, cAll, wPk), e45S = E45(g.S0, g.thF0, cAll, wLk) + E45(g.S0, g.thF0, cAll, wHk);
            d.exp = kYf * (e45P - R * e45S) / (1 - R * fNom);
            (void)sk;
        }
        d.y = r.v; d.s = r.stat; d.pull = (r.v - d.exp) / r.stat;
        d.yNR = rNR.v; d.sNR = rNR.stat; d.pullNR = (rNR.v - d.exp) / rNR.stat;
        d.yE2 = rE2.v; d.sE2 = rE2.stat; d.pullE2 = rE2.v / rE2.stat;
        d.sMC = 0;
        disp.push_back(d);
        sp2 += d.pull * d.pull; spNR2 += d.pullNR * d.pullNR; spE2 += d.pullE2 * d.pullE2; spm += d.pull;
        P("       %.4f | %9.1f %8.1f %8.1f %7.2f | %9.1f %8.1f %7.2f            | %9.1f %8.1f %7.2f\n", c, d.y, d.s, d.exp, d.pull, d.yNR, d.sNR, d.pullNR, d.yE2, d.sE2, d.pullE2);
        RS("disp %d %.4f Y %.4f stat %.4f exp %.4f pull %.4f Ynorefit %.4f stat %.4f pull %.4f YE2 %.4f stat %.4f pull %.4f\n", k, c, d.y, d.s, d.exp, d.pull, d.yNR, d.sNR,
           d.pullNR, d.yE2, d.sE2, d.pullE2);
    }
    const double rms = std::sqrt(sp2 / kNDisp), rmsNR = std::sqrt(spNR2 / kNDisp), rmsE2 = std::sqrt(spE2 / kNDisp), meanP = spm / kNDisp;
    double sdev = 0; for (const auto& d : disp) sdev += (d.pull - meanP) * (d.pull - meanP); sdev = std::sqrt(sdev / (kNDisp - 1));
    const double sigTr = psb.stat * std::sqrt(std::max(0., rms * rms - 1.));
    P("       rms of the pulls (about 0) %.3f (mean %.3f, std. dev. %.3f); without the refit %.3f; E2 only %.3f\n", rms, meanP, sdev, rmsNR, rmsE2);
    P("       sigma_tr = sigma_stat(peak SB) sqrt(max(0, rms^2 - 1)) = %.1f x sqrt(%.3f) = %.1f rows (%.1f %% of Y_peak)\n", psb.stat, std::max(0., rms * rms - 1.), sigTr,
      100. * sigTr / psb.v);
    RS("closure_disp rms %.5f mean %.5f sdev %.5f rms_norefit %.5f rms_E2 %.5f sigma_tr %.4f sigma_stat %.4f rel %.5f\n", rms, meanP, sdev, rmsNR, rmsE2, sigTr, psb.stat,
       sigTr / psb.v);
    // 5b. CR-L and CR-K prediction pulls (E2)
    struct CRp { std::string name; double N, vN, pred, sStat, sMC, pull, s10, rel, dYb; };
    std::vector<CRp> crp;
    {
        const int nI = kNCell * nW;
        auto cr = [&](const std::string& name, const std::vector<int>& cells, int w, const Ctx& CC = Ctx(), bool keep = true) {
            const Ctx& Cx = CC.ok ? CC : C0;
            CRp o; o.name = name;
            const ERes pr = Evaluate("pred_" + name, MakePred(cells, w), Cx, true);
            o.N = SumC(g.DW.N, cells, w);
            o.vN = 0; for (int a : cells) for (int b : cells) o.vN += g.DW.covN[(size_t)DI(a, w) * nI + DI(b, w)];
            // the prediction's stat error: CB, stage A and F-main only (its N part is the CR data used in the fit: none here)
            o.pred = pr.v; o.sStat = std::sqrt(pr.eCB * pr.eCB + pr.eA * pr.eA + pr.eF * pr.eF); o.sMC = pr.mc;
            const double sp = std::sqrt(o.sStat * o.sStat + o.sMC * o.sMC);
            o.pull = (o.N - o.pred) / std::sqrt(sp * sp + o.vN);
            o.s10 = std::max(std::fabs(o.N - o.pred), sp);
            o.rel = o.s10 / SumC(Cx.BA0, cells, w);
            const double sU = SumC(g.DW.N, cE2, iU) - SumC(g.DW.CB, cE2, iU) - EtaAll(g.S0, g.thF0, cE2, iU);
            o.dYb = kYf * o.rel * sU * SumC(Cx.BA0, cE2, iB) / SumC(Cx.BA0, cE2, iU);
            if (!keep) {
                P("       [%s] %-28s N = %8.0f; predicted %9.1f +- %6.1f (stat) +- %6.1f (MC); pull %+.2f\n", Cx.tag.c_str(), name.c_str(), o.N, o.pred, o.sStat, o.sMC, o.pull);
                RS("crpull_%s %s N %.1f pred %.4f stat %.4f mc %.4f pull %.4f\n", Cx.tag.c_str(), name.c_str(), o.N, o.pred, o.sStat, o.sMC, o.pull);
                return;
            }
            P("   5b. %-28s N = %8.0f +- %6.1f; predicted %9.1f +- %6.1f (stat) +- %6.1f (MC); pull %+.2f; S10 term max(|N - pred|, sigma_pred) = %.1f rows = %.1f %% of B_A\n",
              name.c_str(), o.N, std::sqrt(o.vN), o.pred, o.sStat, o.sMC, o.pull, o.s10, 100 * o.rel);
            P("       -> applied to the B4 background in B (E2): %.0f rows = %.1f %% of Y_bump(E2)\n", o.dYb, 100. * o.dYb / b4.v);
            RS("crpull %s N %.1f eN %.4f pred %.4f stat %.4f mc %.4f pull %.4f s10 %.4f rel %.6f dYbump %.4f relY %.6f\n", name.c_str(), o.N, std::sqrt(o.vN), o.pred,
               o.sStat, o.sMC, o.pull, o.s10, o.rel, o.dYb, o.dYb / b4.v);
            crp.push_back(o);
        };
        cr("CR-L (0.5225-0.5300, E2)", cE2, iL);
        cr("CR-K (S5, 0.400-0.5225, E2)", Cells(1, 4), iK);
        cr("0.400-0.455, E2 (below B)", cE2, iR0);
        for (const Ctx* x : {&Cq, &C4s}) {
            if (!x->ok) continue;
            cr("CR-L (0.5225-0.5300, E2)", cE2, iL, *x, false);
            cr("CR-K (S5, 0.400-0.5225, E2)", Cells(1, 4), iK, *x, false);
            cr("0.400-0.455, E2 (below B)", cE2, iR0, *x, false);
        }
        // the same with the WP3 eta (stage A at a 0.4), as WP3 quoted
        const Vec f04 = {0.4, 1., 1., 1.};
        const StageA a04 = FitStageA(g.D0, g.T3, f04, 0.575, 0.700, nullptr, false, nullptr, false);
        const Vec B04 = BATable(g.T3, &g.sm3, a04.S, a04.F.th);
        const double k04 = a04.F.th[a04.S.iK[1]];
        auto pr04 = [&](const std::vector<int>& cells, int w) { return k04 * SumC(g.DW.CB, cells, w) + SumC(B04, cells, w) + EtaAll(g.S0, f04, cells, w); };
        P("       with the WP3 eta values (a = 0.4, g = 1): predicted CR-L %.1f, CR-K %.1f (WP3 stageA.txt: 8078.0, 440.4)\n", pr04(cE2, iL), pr04(Cells(1, 4), iK));
    }

    // ---------------------------------------------------------------- 6. cross-checks from wp6_fill (edge, tag, P_main)
    P("\n6. CROSS-CHECKS (from wp6_fill.txt)\n");
    std::map<std::string, std::vector<std::string>> key;
    {
        std::ifstream in(Form("%s/wp6_fill.txt", kOut));
        std::string line;
        while (std::getline(in, line)) {
            if (line.rfind("KEY ", 0) != 0) continue;
            std::istringstream is(line.substr(4));
            std::string k; is >> k;
            std::vector<std::string> v; std::string x; while (is >> x) v.push_back(x);
            std::string kk = k;
            if (k == "edge" || k == "tagmc" || k == "tagmatched") kk = k + " " + v[0];
            if (k == "tag") kk = k + " " + v[0] + " " + v[1];
            if (k == "pmain") kk = k + " " + v[0];
            if (k == "edgeslice") kk = k + " " + v[0] + " " + v[1];
            key[kk] = v;
        }
    }
    auto num = [&](const std::string& k, int i) { auto it = key.find(k); return (it != key.end() && i < (int)it->second.size()) ? std::atof(it->second[i].c_str()) : NAN; };
    // edge: "edge all E_data v e E_sim v mc v f v stat v mc v model v nvar n"
    const double fEd = num("edge all", 9), fEdSt = num("edge all", 11), fEdMC = num("edge all", 13), fEdMod = num("edge all", 15);
    const double Edat = num("edge all", 2), eEdat = num("edge all", 3), Esim = num("edge all", 5), eEsim = num("edge all", 7);
    P("   f_eta(edge) = E_data/E_sim = %.1f/%.1f = %.4f +- %.4f (stat) +- %.4f (MC) +- %.4f (model); E1 %.4f, E2 %.4f\n", Edat, Esim, fEd, fEdSt, fEdMC, fEdMod,
      num("edge E1", 9), num("edge E2", 9));
    // tag: "tag g e data S eS chi2 ndf free mean sigma S sim S free S mean sigma truth tr trw"
    const double Td = num("tag 0 0", 3), eTd = num("tag 0 0", 4), Ts = num("tag 0 0", 12), Tsf = num("tag 0 0", 14), TdF = num("tag 0 0", 10);
    const double Tmc = num("tagmc 0", 1), TmcF = num("tagmc 0", 2);
    const double TsM = num("tagmatched 0", 3), TsMf = num("tagmatched 0", 4);
    const double fTag = Td / Ts, fTagSt = eTd / Ts, fTagMC = fTag * Tmc / Ts;
    P("   f_eta(tag) = T_data/T_sim = %.1f/%.1f = %.4f +- %.4f (stat) +- %.4f (MC) [frozen data shape]; free shapes %.1f/%.1f = %.4f; sim pi0 window matched: %.4f (frozen), %.4f (free)\n",
      Td, Ts, fTag, fTagSt, fTagMC, TdF, Tsf, TdF / Tsf, Td / TsM, TdF / TsMf);
    const double eShD = num("tag 0 2", 3) / Td, eShS = num("tag 0 2", 12) / Ts, eShT = num("tag 0 2", 18) / num("tag 0 0", 18);
    P("   tag share with M(ee) >= 10 MeV: data %.3f +- %.3f, simulation %.3f (frozen fit), %.3f (truth C2+C3); bump extrapolation factor f_E2 (sim, B) = %.3f\n", eShD,
      std::sqrt(eShD * (1 - eShD) / Td * 1.6), eShS, eShT, fE2);
    // P_main
    const double Spm = num("pmain data", 1), eSpm = num("pmain data", 2), sprPm = num("pmain data", 3), SwinPm = num("pmain data", 9);
    const double Spm1 = num("pmain sim0", 1), Spm089 = num("pmain sim1", 1), truthP = num("pmain sim0", 7);
    P("   P_main on S0: S = %.1f +- %.1f (stat) +- %.1f (spread); in P x 64/63: %.1f; sim T1+T2 with the same frozen fit: %.1f (r = 1), %.1f (r = 0.89); S_data/S_sim = %.4f / %.4f\n",
      Spm, eSpm, sprPm, kYf * SwinPm, Spm1, Spm089, Spm / Spm1, Spm / Spm089);
    RS("xcheck f_edge %.6f %.6f %.6f %.6f fE1 %.6f fE2 %.6f\n", fEd, fEdSt, fEdMC, fEdMod, num("edge E1", 9), num("edge E2", 9));
    RS("xcheck f_tag %.6f %.6f %.6f free %.6f matched %.6f %.6f shareE2 data %.5f sim %.5f truth %.5f\n", fTag, fTagSt, fTagMC, TdF / Tsf, Td / TsM, TdF / TsMf, eShD, eShS, eShT);
    RS("xcheck pmain S %.4f %.4f %.4f Swin64 %.4f simS %.4f %.4f ratio %.6f %.6f\n", Spm, eSpm, sprPm, kYf * SwinPm, Spm1, Spm089, Spm / Spm1, Spm / Spm089);
    // F-main values of WP5 (cross-check)
    const double fmYb = 19426.131, fmYbT = 2478.0, fmYp = 5044.7742, fmYpT = 804.37;
    {
        std::ifstream in(kWP5);
        std::string line;
        // (read Y_bump and Y_peak of FMAIN to stay exact)
        while (std::getline(in, line)) {
            if (line.rfind("derived FMAIN Y_bump ", 0) == 0 || line.rfind("derived FMAIN Y_peak ", 0) == 0) P("   WP5 %s\n", line.substr(8, 140).c_str());
        }
    }

    // ---------------------------------------------------------------- 7. acceptance and summary
    P("\n7. SUMMARY (S0, isBest >= 0; raw rows x 64/63; L = %.1f nb^-1)\n", kL);
    P("   Y_bump (B4, M(ee) >= 10 MeV)      = %8.0f +- %5.0f (stat) +- %5.0f (MC); bias-corrected %8.0f [WP4: FAIL marginal; PASS with bias correction]\n", b4.v, b4.stat, b4.mc, b4.bc);
    P("   Y_bump (B4, all M(ee), extrapol.) = %8.0f +- %5.0f (stat) +- %5.0f (MC)\n", b4all.v, b4all.stat, b4all.mc);
    P("   Y_peak (side bands, all M(ee))    = %8.0f +- %5.0f (stat) +- %5.0f (MC) +- %5.0f (sigma_tr) [WP4: PASS]\n", psb.v, psb.stat, psb.mc, sigTr);
    P("   F-main (WP5, cross-check)         : Y_bump %.0f +- %.0f, Y_peak %.0f +- %.0f (total)\n", fmYb, fmYbT, fmYp, fmYpT);

    // ---------------------------------------------------------------- 8. Table E3
    {
        FILE* t = fopen(Form("%s/wp6_tableE3.txt", kTab), "w");
        fprintf(t, "# Table E3 -- side-band estimators and cross-checks (WP6; plots/all_report/wp6_ana.C, wp6_fill.C), 2026-11-14\n");
        fprintf(t, "# data S0 (pre + both pions in the vertex fit), isBest >= 0, L = %.1f nb^-1; yields = raw rows x 64/63, not efficiency corrected\n", kL);
        fprintf(t, "# errors: stat = data N (event groups) + CB (200 like-sign replicas) + stage A (sandwich) + F-main 'other eta' (sandwich);\n");
        fprintf(t, "#         MC = RMS over 100 simulation replica chains (F-main -> stage A -> estimator); bias-corr = 2 v - replica mean (WP4 'RB')\n");
        fprintf(t, "%-52s %11s %9s %9s %9s %11s %s\n", "quantity", "value", "stat", "MC", "total", "bias-corr", "status / note");
        auto row = [&](const std::string& n, double v, double s, double m, double bc, const std::string& note, const char* fmt = "%11.1f") {
            fprintf(t, "%-52s ", n.c_str());
            fprintf(t, fmt, v);
            fprintf(t, " %9.3g %9.3g %9.3g ", s, m, std::sqrt(s * s + m * m));
            if (bc == bc) fprintf(t, fmt, bc); else fprintf(t, "%11s", "-");
            fprintf(t, " %s\n", note.c_str());
        };
        row("Y_bump B4, M(ee) >= 10 MeV [rows]", b4.v, b4.stat, b4.mc, b4.bc, "MAIN (bump); WP4 FAIL marginal, PASS with bias correction; EXTRAPOLATED (CR-U only)");
        row("Y_bump B4, all M(ee) (E2 / f_E2) [rows]", b4all.v, b4all.stat, b4all.mc, b4all.bc, Form("f_E2 = %.4f from the simulation (kappa %.2f)", fE2, g.thF0[1]));
        row("Y_peak side bands, all M(ee) [rows]", psb.v, psb.stat, psb.mc, psb.bc, Form("MAIN (peak); WP4 PASS; + sigma_tr %.0f", sigTr));
        row("Y_peak side bands, M(ee) < 10 MeV [rows]", psbE1.v, psbE1.stat, psbE1.mc, psbE1.bc, "per class");
        row("Y_peak side bands, M(ee) >= 10 MeV [rows]", psbE2.v, psbE2.stat, psbE2.mc, psbE2.bc, "per class");
        row("sigma_vis bump (B4 all M(ee)) [nb]", b4all.v / kL, b4all.stat / kL, b4all.mc / kL, b4all.bc / kL, "Y/L", "%11.4f");
        row("sigma_vis peak (side bands) [nb]", psb.v / kL, psb.stat / kL, psb.mc / kL, psb.bc / kL, "Y/L", "%11.4f");
        row("bump/SMASH, M(ee) >= 10 MeV (B4)", aB4.v, aB4.stat, aB4.mc, aB4.bc, "Y_B4 / (C2 + C3 at SMASH x PDG BR)", "%11.4f");
        row("peak/SMASH, all M(ee) (side bands)", aPSB.v, aPSB.stat, aPSB.mc, aPSB.bc, "Y_SB / (C4 + C5 at SMASH x PDG BR)", "%11.4f");
        row("DR (side bands, all M(ee))", DRall.v, DRall.stat, DRall.mc, DRall.bc, "(Y_peak/Y_bump)_data / (..)_sim; bump extrapolated", "%11.4f");
        row("DR (side bands, both M(ee) >= 10 MeV)", DRE2.v, DRE2.stat, DRE2.mc, DRE2.bc, "no M(ee) extrapolation", "%11.4f");
        row("peak share M(ee) >= 10 MeV (side bands)", shE2.v, shE2.stat, shE2.mc, shE2.bc, "", "%11.4f");
        fprintf(t, "# closure\n");
        fprintf(t, "%-52s %11.3f   rms of 13 pulls (mean %.3f); sigma_tr = %.1f rows (%.1f %% of Y_peak) -> S10 (peak)\n", "displaced windows 0.575-0.695", rms, meanP, sigTr,
                100 * sigTr / psb.v);
        for (const auto& c : crp)
            fprintf(t, "%-52s %11.2f   pull; N %.0f, pred %.1f +- %.1f; S10 term %.0f rows (%.1f %% of B_A) -> %.0f rows = %.1f %% of Y_bump(E2)\n", c.name.c_str(), c.pull, c.N,
                    c.pred, std::sqrt(c.sStat * c.sStat + c.sMC * c.sMC), c.s10, 100 * c.rel, c.dYb, 100 * c.dYb / b4.v);
        fprintf(t, "# cross-checks (not combined, D14)\n");
        fprintf(t, "%-52s %11.4f %9.4f %9.4f %9.4f %11s model %.4f\n", "f_eta(edge) = E_data/E_sim (T1+T2)", fEd, fEdSt, fEdMC, std::sqrt(fEdSt * fEdSt + fEdMC * fEdMC + fEdMod * fEdMod), "-", fEdMod);
        fprintf(t, "%-52s %11.4f %9.4f %9.4f %9.4f %11s sim tag emulated (UNVERIFIED); matched pi0 window %.4f\n", "f_eta(tag) = T_data/T_sim", fTag, fTagSt, fTagMC,
                std::sqrt(fTagSt * fTagSt + fTagMC * fTagMC), "-", Td / TsM);
        fprintf(t, "%-52s %11.4f %9.4f %9.4f %9.4f %11s P_main %.0f +- %.0f +- %.0f (spread) vs sim %.0f (r = 1)\n", "P_main S_data/S_sim (same fit)", Spm / Spm1, eSpm / Spm1, 0.,
                std::sqrt(eSpm * eSpm + sprPm * sprPm) / Spm1, "-", Spm, eSpm, sprPm, Spm1);
        fprintf(t, "%-52s %11.4f %9.4f %9.4f %9.4f %11s WP5, F10 E2 (CROSS-CHECK, WP4 FAIL)\n", "a_eta (F-main)", g.thF0[0], std::sqrt(g.CF[0][0]), std::sqrt(g.CFmc[0][0]),
                std::sqrt(g.CF[0][0] + g.CFmc[0][0]), "-");
        fprintf(t, "%-52s %11.1f %9s %9s %9.1f %11s WP5 (CROSS-CHECK)\n", "Y_bump F-main [rows]", fmYb, "-", "-", fmYbT, "-");
        fprintf(t, "%-52s %11.1f %9s %9s %9.1f %11s WP5 (CROSS-CHECK)\n", "Y_peak F-main [rows]", fmYp, "-", "-", fmYpT, "-");
        fprintf(t, "%-52s %11.1f %9.1f %9s %9.1f %11s whole Gaussian %.0f; model spread %.0f\n", "P_main in P x 64/63 [rows]", kYf * SwinPm, kYf * eSpm * SwinPm / Spm, "-",
                kYf * std::sqrt(eSpm * eSpm + sprPm * sprPm) * SwinPm / Spm, "-", Spm, sprPm);
        fprintf(t, "# variants\n");
        for (const auto& v : V) {
            fprintf(t, "%-52s %11.1f %9.3g %9.3g %9.3g ", v.name.substr(0, 52).c_str(), v.v, v.stat, v.mc, v.tot);
            if (v.nrep > 0) fprintf(t, "%11.1f\n", v.bc); else fprintf(t, "%11s %s\n", "-", v.mc > 0 ? "(MC from the WP3 parameter covariance)" : "(no MC)");
        }
        fclose(t);
        // LaTeX snippet
        FILE* x = fopen(Form("%s/wp6_tableE3.tex", kTab), "w");
        fprintf(x, "%% Table E3 (WP6) -- plots/all_report/wp6_ana.C; numbers in tables_all/wp6_tableE3.txt\n");
        fprintf(x, "\\begin{tabular}{lrrrrl}\n\\toprule\nquantity & value & stat & MC & bias-corr. & note \\\\\n\\midrule\n");
        auto tr = [&](const char* n, const ERes& e, const char* fmt, const char* note) {
            fprintf(x, "%s & ", n); fprintf(x, fmt, e.v); fprintf(x, " & "); fprintf(x, fmt, e.stat); fprintf(x, " & "); fprintf(x, fmt, e.mc); fprintf(x, " & ");
            fprintf(x, fmt, e.bc); fprintf(x, " & %s \\\\\n", note);
        };
        tr("$Y_\\mathrm{bump}$ (B4), $M_{ee} \\ge 10$\\,MeV", b4, "%.0f", "main; WP4 fail (m.)");
        tr("$Y_\\mathrm{bump}$ (B4), all $M_{ee}$", b4all, "%.0f", "extrapolated in $M_{ee}$");
        tr("$Y_\\mathrm{peak}$ (side bands), all $M_{ee}$", psb, "%.0f", "main; WP4 pass");
        tr("$Y_\\mathrm{peak}$, $M_{ee} < 10$\\,MeV", psbE1, "%.0f", "");
        tr("$Y_\\mathrm{peak}$, $M_{ee} \\ge 10$\\,MeV", psbE2, "%.0f", "");
        tr("bump / SMASH ($M_{ee} \\ge 10$\\,MeV)", aB4, "%.3f", "");
        tr("peak / SMASH", aPSB, "%.3f", "");
        tr("DR (side bands)", DRall, "%.3f", "");
        tr("DR ($M_{ee} \\ge 10$\\,MeV)", DRE2, "%.3f", "");
        fprintf(x, "\\midrule\n");
        fprintf(x, "$\\sigma_\\mathrm{tr}$ (13 displaced windows) & %.0f & \\multicolumn{3}{l}{rms of pulls %.2f} & S10 (peak) \\\\\n", sigTr, rms);
        for (const auto& c : crp) {
            std::string nm = c.name.substr(0, c.name.find(','));
            fprintf(x, "pull %s & %+.2f & \\multicolumn{3}{l}{S10 %.0f rows} & %.1f\\,\\%% of $Y_\\mathrm{bump}(E2)$ \\\\\n", nm.c_str(), c.pull, c.dYb, 100 * c.dYb / b4.v);
        }
        fprintf(x, "\\midrule\n");
        fprintf(x, "$f_\\eta$(edge) & %.3f & %.3f & %.3f & -- & model %.3f \\\\\n", fEd, fEdSt, fEdMC, fEdMod);
        fprintf(x, "$f_\\eta$(tag) & %.3f & %.3f & %.3f & -- & emulated tag, UNVERIFIED \\\\\n", fTag, fTagSt, fTagMC);
        fprintf(x, "$P_\\mathrm{main}$ $S$ & %.0f & %.0f & -- & -- & spread %.0f \\\\\n", Spm, eSpm, sprPm);
        fprintf(x, "$a_\\eta$ (F-main, WP5) & %.3f & %.3f & %.3f & -- & cross-check \\\\\n", g.thF0[0], std::sqrt(g.CF[0][0]), std::sqrt(g.CFmc[0][0]));
        fprintf(x, "\\bottomrule\n\\end{tabular}\n");
        fclose(x);
    }

    // ---------------------------------------------------------------- 9. Fig. E4
    Style();
    {
        TCanvas* c = new TCanvas("cE4", "", 1500, 1200);
        c->Divide(2, 2, 0.004, 0.004);
        const std::vector<double> e = F10Edges();
        // ---- (a) B4, E2
        c->cd(1);
        gPad->SetLeftMargin(0.16); gPad->SetRightMargin(0.03);
        {
            TH1D* hd = new TH1D("a_d", "", kNF10, e.data());
            TH1D* hb = new TH1D("a_b", "", kNF10, e.data());
            TH1D* ho = new TH1D("a_o", "", kNF10, e.data());
            TH1D* hp = new TH1D("a_p", "", kNF10, e.data());
            const Spec& S = C0.A.S; const Vec& th = C0.A.F.th;
            const double sU = SumC(g.DW.N, cE2, iU) - SumC(g.DW.CB, cE2, iU) - EtaAll(g.S0, g.thF0, cE2, iU);
            const double norm = sU / SumC(C0.BA0, cE2, iU);
            const double a = g.thF0[0], k = g.thF0[1], r = g.thF0[2], gg = g.thF0[3];
            for (int x = 1; x <= kNF10; ++x) {
                const double bwm = 1e3 * (e[x] - e[x - 1]);
                double y = 0, v = 0, b = 0, o = 0, p = 0;
                for (int j = 0; j < kNSlice; ++j) {
                    const int h = 2 * j + 1;
                    y += g.D0.y[h][x] - g.D0.cb[h][x];
                    v += g.D0.y[h][x] + g.D0.vcb[h][x];
                    b += BAint(g.T3, &g.sm3, S, th, h, e[x - 1], e[x]);
                    o += a * k * g.T0.s[wp5::kS4][h][x] + a * r * g.T0.s[wp5::kS5][h][x] + a * gg * g.T0.s[wp5::kSeta][h][x];
                    p += a * g.T0.s[wp5::kS2][h][x] + a * k * g.T0.s[wp5::kS3][h][x];
                }
                hd->SetBinContent(x, y / bwm); hd->SetBinError(x, std::sqrt(v) / bwm);
                hb->SetBinContent(x, norm * b / bwm); ho->SetBinContent(x, o / bwm); hp->SetBinContent(x, p / bwm);
            }
            hb->SetFillColor(kGray); hb->SetLineColor(kGray + 2);
            ho->SetFillColor(kViolet - 9); ho->SetLineColor(kViolet + 2);
            hp->SetFillColor(kOrange - 4); hp->SetLineColor(kOrange + 7);
            THStack* st = new THStack("a_st", "");
            st->Add(hb); st->Add(ho); st->Add(hp);
            hd->SetMarkerStyle(20); hd->SetMarkerSize(0.8); hd->SetMarkerColor(kBlue + 1); hd->SetLineColor(kBlue + 1);
            hd->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]");
            hd->GetYaxis()->SetTitle("N - CB [counts / (MeV/c^{2})]");
            hd->SetMinimum(0); hd->SetMaximum(1.75 * hd->GetMaximum());
            hd->Draw("E1");
            st->Draw("hist same");
            hd->Draw("E1 same");
            const double ymax = hd->GetMaximum();
            TBox* bB = new TBox(kBlo, 0, kBhi, 0.04 * ymax); bB->SetFillColor(kAzure - 4); bB->Draw();
            TBox* bU = new TBox(0.575, 0, 0.700, 0.04 * ymax); bU->SetFillColor(kGreen - 6); bU->Draw();
            TLegend* lg = new TLegend(0.19, 0.45, 0.75, 0.70);
            lg->SetTextSize(0.036);
            lg->AddEntry(hd, "data N - CB, M(ee) #geq 10 MeV/c^{2}", "pe");
            lg->AddEntry(hb, "stage-A background, normalised on U", "f");
            lg->AddEntry(ho, "other #eta (C4, C5, G_{#eta}; F-main values)", "f");
            lg->AddEntry(hp, "bump C2 + #kappaC3 (F-main values)", "f");
            lg->Draw();
            Txt(0.19, 0.87, "(a) B4, S0, isBest #geq 0, L = 6079 nb^{-1}", 0.042);
            Txt(0.19, 0.81, Form("Y_{bump}(E2) = %.0f #pm %.0f (stat) #pm %.0f (MC)", b4.v, b4.stat, b4.mc), 0.040, kBlue + 2);
            Txt(0.19, 0.75, "coloured bars at the bottom: B (blue), U (green)", 0.034);
        }
        // ---- (b) peak side bands, all M(ee)
        c->cd(2);
        gPad->SetLeftMargin(0.16); gPad->SetRightMargin(0.03);
        {
            const double lo = 0.515, hi = 0.575;
            std::vector<double> ee;
            for (double x = lo; x < hi + 1e-9; x += 0.0025) ee.push_back(x);
            const int nb = (int)ee.size() - 1;
            TH1D* hd = new TH1D("b_d", "", nb, ee.data());
            for (int x = 1; x <= kNF10; ++x) {
                const double xc = 0.5 * (e[x - 1] + e[x]);
                if (xc < lo || xc > hi) continue;
                double y = 0, v = 0;
                for (int h = 0; h < 10; ++h) { y += g.D0.y[h][x] - g.D0.cb[h][x]; v += g.D0.y[h][x] + g.D0.vcb[h][x]; }
                const int b = hd->FindBin(xc);
                const double bwm = 1e3 * (e[x] - e[x - 1]);
                hd->SetBinContent(b, hd->GetBinContent(b) + y / bwm); hd->SetBinError(b, std::sqrt(std::pow(hd->GetBinError(b), 2) + v / (bwm * bwm)));
            }
            // background model scaled to the side bands (as R' does): (B_A + tail + G_eta) x [(N - CB)_SB - f Y] / [same](SB)
            const Spec& S = C0.A.S; const Vec& th = C0.A.F.th;
            const double sS = SumC(g.DW.N, cAll, iSl) - SumC(g.DW.CB, cAll, iSl) + SumC(g.DW.N, cAll, iSh) - SumC(g.DW.CB, cAll, iSh);
            const double bS = SumC(C0.BA0, cAll, iSl) + SumC(C0.BA0, cAll, iSh) + Tail(g.S0, g.thF0, cAll, iSl) + Tail(g.S0, g.thF0, cAll, iSh);
            const double scl = (sS - fNom * psb.v / kYf) / bS;
            // per 2.5 MeV: B_A (smooth) + bump tail + G_eta (simulation, 1 MeV templates with fractional bins), and + the peak
            // shape C4 + C5 scaled to the extracted yield
            TH1D* gb = new TH1D("b_bk", "", nb, ee.data());
            TH1D* gt = new TH1D("b_tot", "", nb, ee.data());
            TFile* ff = TFile::Open(Form("%s/wp6_fill.root", kOut));
            auto tint = [&](TH1* t, double a, double b) {
                double v = 0;
                for (int q = 1; q <= t->GetNbinsX(); ++q) {
                    const double o = std::min(b, t->GetXaxis()->GetBinUpEdge(q)) - std::max(a, t->GetXaxis()->GetBinLowEdge(q));
                    if (o > 0) v += t->GetBinContent(q) * o / t->GetXaxis()->GetBinWidth(q);
                }
                return v;
            };
            double u45P = U45(g.S0, cAll, iP);
            for (int q = 1; q <= nb; ++q) {
                const double m1 = ee[q - 1], m2 = ee[q];
                double b = 0, pk = 0;
                for (int h = 0; h < 10; ++h) {
                    const int j = h / 2;
                    b += BAint(g.T3, &g.sm3, S, th, h, m1, m2);
                    const std::string sfx = Form("_S%dE%d", j + 1, h % 2 + 1);
                    b += g.thF0[0] * tint((TH1*)ff->Get(("sim/tpl_C2" + sfx).c_str()), m1, m2) +
                         g.thF0[0] * g.thF0[1] * tint((TH1*)ff->Get(("sim/tpl_C3" + sfx).c_str()), m1, m2) +
                         g.thF0[0] * g.thF0[3] * tint((TH1*)ff->Get(("sim/tpl_Geta" + sfx).c_str()), m1, m2);
                    pk += tint((TH1*)ff->Get(("sim/tpl_C4" + sfx).c_str()), m1, m2) + tint((TH1*)ff->Get(("sim/tpl_C5" + sfx).c_str()), m1, m2);
                }
                const double bwm = 1e3 * (m2 - m1);
                gb->SetBinContent(q, scl * b / bwm);
                gt->SetBinContent(q, (scl * b + (psb.v / kYf) * pk / u45P) / bwm);
            }
            ff->Close();
            hd->SetMarkerStyle(20); hd->SetMarkerSize(0.9); hd->SetMarkerColor(kBlue + 1); hd->SetLineColor(kBlue + 1);
            hd->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]");
            hd->GetYaxis()->SetTitle("N - CB [counts / (MeV/c^{2})]");
            hd->SetMinimum(0); hd->SetMaximum(1.30 * hd->GetMaximum());
            hd->Draw("E1");
            const double ymax = hd->GetMaximum();
            TBox* bp = new TBox(kPlo, 0, kPhi, ymax); bp->SetFillColorAlpha(kRed - 9, 0.35); bp->SetLineWidth(0); bp->Draw();
            TBox* bl = new TBox(0.5250, 0, 0.5336, ymax); bl->SetFillColorAlpha(kGreen - 9, 0.45); bl->SetLineWidth(0); bl->Draw();
            TBox* bh = new TBox(0.5542, 0, 0.5628, ymax); bh->SetFillColorAlpha(kGreen - 9, 0.45); bh->SetLineWidth(0); bh->Draw();
            hd->Draw("E1 same");
            gb->SetLineColor(kGreen + 3); gb->SetLineWidth(3); gb->Draw("hist same");
            gt->SetLineColor(kRed + 1); gt->SetLineWidth(2); gt->SetLineStyle(2); gt->Draw("hist same");
            hd->Draw("E1 same");
            TLegend* lg = new TLegend(0.19, 0.16, 0.80, 0.36);
            lg->SetTextSize(0.036);
            lg->AddEntry(hd, "data N - CB, all M(ee)", "pe");
            lg->AddEntry(gb, "B_{A} + bump tail + G_{#eta}, scaled to the side bands", "l");
            lg->AddEntry(gt, "+ C4 + C5 shape (sim) at the extracted yield", "l");
            lg->Draw();
            Txt(0.19, 0.86, "(b) peak side bands (green), window P (red)", 0.045);
            Txt(0.19, 0.80, Form("Y_{peak} = %.0f #pm %.0f (stat) #pm %.0f (MC)", psb.v, psb.stat, psb.mc), 0.040, kBlue + 2);
            Txt(0.19, 0.75, Form("R' = %.3f, leakage f = %.3f", (SumC(C0.BA0, cAll, iP) + Tail(g.S0, g.thF0, cAll, iP)) / bS, fNom), 0.036);
        }
        // ---- (c) closure pulls
        c->cd(3);
        gPad->SetLeftMargin(0.16); gPad->SetRightMargin(0.03);
        {
            TH1D* fr = new TH1D("c_fr", "", 10, 0.56, 0.71);
            fr->SetMinimum(-6.5); fr->SetMaximum(4.5);
            fr->GetXaxis()->SetTitle("centre of the displaced window [GeV/c^{2}]");
            fr->GetYaxis()->SetTitle("pull = (Y - Y_{exp}) / #sigma_{stat}");
            fr->Draw("AXIS");
            TBox* b1 = new TBox(0.56, -1, 0.71, 1); b1->SetFillColor(kGray); b1->Draw();
            TLine* l0 = new TLine(0.56, 0, 0.71, 0); l0->SetLineStyle(2); l0->Draw();
            TGraph* gp = new TGraph(), *gn = new TGraph(), *g2 = new TGraph();
            for (int k = 0; k < kNDisp; ++k) { gp->SetPoint(k, disp[k].c, disp[k].pull); gn->SetPoint(k, disp[k].c + 0.0015, disp[k].pullNR); g2->SetPoint(k, disp[k].c - 0.0015, disp[k].pullE2); }
            gp->SetMarkerStyle(20); gp->SetMarkerSize(1.3); gp->Draw("P");
            gn->SetMarkerStyle(24); gn->SetMarkerSize(1.1); gn->SetMarkerColor(kGray + 2); gn->Draw("P");
            g2->SetMarkerStyle(25); g2->SetMarkerSize(1.1); g2->SetMarkerColor(kBlue + 1); g2->Draw("P");
            fr->Draw("AXIS same");
            TLegend* lg = new TLegend(0.19, 0.16, 0.95, 0.32);
            lg->SetTextSize(0.033); lg->SetNColumns(1);
            lg->AddEntry(gp, Form("all M(ee), window removed from CR-U, refit (rms %.2f)", rms), "p");
            lg->AddEntry(gn, Form("all M(ee), nominal stage A (rms %.2f)", rmsNR), "p");
            lg->AddEntry(g2, Form("M(ee) #geq 10 MeV/c^{2} only (rms %.2f)", rmsE2), "p");
            lg->Draw();
            Txt(0.19, 0.86, "(c) side-band estimator at 13 displaced windows", 0.045);
            Txt(0.19, 0.80, Form("#sigma_{tr} = #sigma_{stat} #sqrt{rms^{2} - 1} = %.0f counts (%.1f%% of Y_{peak})", sigTr, 100 * sigTr / psb.v), 0.038);
            Txt(0.19, 0.75, Form("CR-L pull %+.2f, CR-K pull %+.2f (E2, stage A + F-main)", crp[0].pull, crp[1].pull), 0.038);
        }
        // ---- (d) eta scales
        c->cd(4);
        gPad->SetLeftMargin(0.40); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.10);
        {
            struct Pt { std::string n; double v, s, t; int col; };
            std::vector<Pt> pts = {
                {"a_{#eta} F-main (WP5, cross-check)", g.thF0[0], std::sqrt(g.CF[0][0]), std::sqrt(g.CF[0][0] + g.CFmc[0][0]), kGray + 2},
                {"bump B4 / SMASH, M_{ee} #geq 10 MeV", aB4.v, aB4.stat, aB4.tot, kOrange + 7},
                {"peak side bands / SMASH", aPSB.v, aPSB.stat, aPSB.tot, kRed + 1},
                {"f_{#eta}(edge) = E_{data}/E_{sim}", fEd, fEdSt, std::sqrt(fEdSt * fEdSt + fEdMC * fEdMC + fEdMod * fEdMod), kGreen + 3},
                {"f_{#eta}(tag) = T_{data}/T_{sim} (UNVERIFIED)", fTag, fTagSt, std::sqrt(fTagSt * fTagSt + fTagMC * fTagMC), kAzure + 2},
                {"P_{main}: S_{data}/S_{sim}", Spm / Spm1, eSpm / Spm1, std::sqrt(eSpm * eSpm + sprPm * sprPm) / Spm1, kViolet + 1}};
            const int n = (int)pts.size();
            TH2D* fr = new TH2D("d_fr", "", 10, 0.0, 0.9, n, 0, n);
            for (int i = 0; i < n; ++i) fr->GetYaxis()->SetBinLabel(n - i, pts[i].n.c_str());
            fr->GetYaxis()->SetLabelSize(0.042);
            fr->GetXaxis()->SetTitle("data / simulation (SMASH #sigma #times PDG BR)");
            fr->Draw("AXIS");
            for (int i = 0; i < n; ++i) {
                const double y = n - i - 0.5;
                TGraphErrors* gt = new TGraphErrors(1); gt->SetPoint(0, pts[i].v, y); gt->SetPointError(0, pts[i].t, 0);
                gt->SetLineColor(pts[i].col); gt->SetLineWidth(1); gt->SetMarkerStyle(1); gt->Draw("P");
                TGraphErrors* gs = new TGraphErrors(1); gs->SetPoint(0, pts[i].v, y); gs->SetPointError(0, pts[i].s, 0);
                gs->SetLineColor(pts[i].col); gs->SetLineWidth(4); gs->SetMarkerStyle(20); gs->SetMarkerColor(pts[i].col); gs->SetMarkerSize(1.3); gs->Draw("P");
            }
            TLine* l1 = new TLine(g.thF0[0], 0, g.thF0[0], n); l1->SetLineStyle(2); l1->SetLineColor(kGray + 1); l1->Draw();
            Txt(0.02, 0.975, "(d) #eta scale estimates relative to SMASH #times PDG BR (thick: stat, thin: total)", 0.038, kBlack, 13);
        }
        Save(c, Form("%s/wp6_E4_cnt", kFig));
    }
    // ---- edge-step figure (data vs simulation)
    {
        TFile* ff = TFile::Open(Form("%s/wp6_fill.root", kOut));
        TCanvas* c = new TCanvas("cEdge", "", 1800, 900);
        c->Divide(4, 2, 0.003, 0.003);
        const char* sl[4] = {"0.30-0.34", "0.34-0.37", "0.37-0.39", "0.39-0.40"};
        for (int s = 0; s < 4; ++s)
            for (int row = 0; row < 2; ++row) {
                c->cd(1 + s + 4 * row);
                gPad->SetLeftMargin(0.17); gPad->SetRightMargin(0.03);
                TH1* h = (TH1*)ff->Get(Form("%s/edge_%s_sl%d_e0", row == 0 ? "data" : "sim", row == 0 ? "sig" : "T12", s));
                TGraph* gf = (TGraph*)ff->Get(Form("%s/edgefit_sl%d_e0", row == 0 ? "data" : "sim", s));
                TGraph* g0 = (TGraph*)ff->Get(Form("%s/edgefit0_sl%d_e0", row == 0 ? "data" : "sim", s));
                if (!h || !gf) continue;
                double lo = gf->GetX()[0] - 0.005, hi = gf->GetX()[gf->GetN() - 1] + 0.005;
                h->GetXaxis()->SetRangeUser(lo, hi);
                h->SetMarkerStyle(20); h->SetMarkerSize(0.7);
                h->SetMarkerColor(row == 0 ? kBlue + 1 : kBlack); h->SetLineColor(row == 0 ? kBlue + 1 : kBlack);
                h->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]");
                h->GetYaxis()->SetTitle(row == 0 ? "data N - CB [counts / 2.5 MeV/c^{2}]" : "sim T1+T2 [exp. counts / 2.5 MeV/c^{2}]");
                h->GetXaxis()->SetNdivisions(505);
                double mx = 0, mn = 1e9;
                for (int b = h->FindBin(lo); b <= h->FindBin(hi); ++b) { mx = std::max(mx, h->GetBinContent(b) + h->GetBinError(b)); mn = std::min(mn, h->GetBinContent(b) - h->GetBinError(b)); }
                h->SetMaximum(1.25 * mx); h->SetMinimum(std::min(0., mn));
                h->Draw("E1");
                TH1* hc = nullptr;
                if (row == 1) {
                    hc = (TH1*)ff->Get(Form("sim/edge_C23_sl%d_e0", s));
                    hc->SetLineColor(kOrange + 7); hc->SetLineWidth(2); hc->SetFillColorAlpha(kOrange - 4, 0.5); hc->Draw("hist same"); h->Draw("E1 same");
                }
                gf->SetLineColor(kGreen + 2); gf->SetLineWidth(3); gf->Draw("L");
                g0->SetLineColor(kGreen + 2); g0->SetLineStyle(2); g0->SetLineWidth(2); g0->Draw("L");
                const double hh = num(Form("edgeslice all %d", s), row == 0 ? 5 : 10), eh = num(Form("edgeslice all %d", s), row == 0 ? 6 : 11);
                Txt(0.20, 0.86, Form("%s, M(#pi^{+}#pi^{-}) %s", row == 0 ? "data S0" : "simulation", sl[s]), 0.050);
                Txt(0.20, 0.79, Form("h = %.0f #pm %.0f", hh, eh), 0.050, row == 0 ? kBlue + 2 : kBlack);
                if (row == 1 && s == 0) {
                    TLegend* lg = new TLegend(0.20, 0.55, 0.70, 0.75); lg->SetTextSize(0.045);
                    lg->AddEntry(h, "T1 + T2 at SMASH", "pe"); lg->AddEntry(hc, "of it: bump C2 + C3", "f"); lg->AddEntry(gf, "pol2 + step", "l"); lg->Draw();
                }
            }
        Save(c, Form("%s/wp6_edge_cnt", kFig));
        // ---- tag figure
        TCanvas* t = new TCanvas("cTag", "", 1800, 600);
        t->Divide(3, 1, 0.003, 0.003);
        for (int p = 0; p < 2; ++p) {
            t->cd(p + 1);
            gPad->SetLeftMargin(0.17); gPad->SetRightMargin(0.03);
            TH1* h = (TH1*)ff->Get(Form("%s/tag_%s_r0_e0", p == 0 ? "data" : "sim", p == 0 ? "sig" : "T12"));
            TGraph* gf = (TGraph*)ff->Get(Form("%s/tagfit_r0_e0", p == 0 ? "data" : "sim"));
            TGraph* g0 = (TGraph*)ff->Get(Form("%s/tagfit0_r0_e0", p == 0 ? "data" : "sim"));
            h->GetXaxis()->SetRangeUser(0.45, 0.66);
            h->SetMarkerStyle(20); h->SetMarkerSize(0.8);
            h->SetMarkerColor(p == 0 ? kBlue + 1 : kBlack); h->SetLineColor(p == 0 ? kBlue + 1 : kBlack);
            h->GetXaxis()->SetTitle("m_{5c} = M(#pi^{+}#pi^{-}e^{+}e^{-}#gamma), #gamma rescaled to m_{#pi^{0}} [GeV/c^{2}]");
            h->GetYaxis()->SetTitle(p == 0 ? "data N - CB [counts / 2.5 MeV/c^{2}]" : "sim T1+T2 [exp. counts / 2.5 MeV/c^{2}]");
            h->SetMinimum(0); h->SetMaximum(1.3 * h->GetMaximum());
            h->Draw("E1");
            if (p == 1) { TH1* hc = (TH1*)ff->Get("sim/tag_C23_r0_e0"); hc->SetLineColor(kOrange + 7); hc->SetFillColorAlpha(kOrange - 4, 0.5); hc->Draw("hist same"); h->Draw("E1 same"); }
            if (gf) { gf->SetLineColor(kGreen + 2); gf->SetLineWidth(3); gf->Draw("L"); }
            if (g0) { g0->SetLineColor(kGreen + 2); g0->SetLineStyle(2); g0->Draw("L"); }
            Txt(0.20, 0.86, p == 0 ? "data S0, ECAL-tagged" : "simulation, emulated tag (UNVERIFIED)", 0.045);
            Txt(0.20, 0.80, Form("T = %.0f%s", p == 0 ? Td : Ts, p == 0 ? Form(" #pm %.0f", eTd) : Form(" #pm %.0f (MC)", Tmc)), 0.045, p == 0 ? kBlue + 2 : kBlack);
            if (p == 1) Txt(0.20, 0.74, Form("T_{data}/T_{sim} = %.3f", fTag), 0.045);
        }
        t->cd(3);
        gPad->SetLeftMargin(0.17); gPad->SetRightMargin(0.03);
        {
            TH1* b = (TH1*)ff->Get("sim/meeg_C23");
            TH1* a = b;
            b->GetXaxis()->SetRangeUser(0.05, 0.25);
            b->SetLineColor(kOrange + 7); b->SetFillColorAlpha(kOrange - 4, 0.6);
            b->GetXaxis()->SetTitle("M(e^{+}e^{-}#gamma) [GeV/c^{2}]");
            b->GetYaxis()->SetTitle("sim bump C2 + C3 [exp. counts / MeV/c^{2}]");
            b->SetMaximum(1.6 * b->GetMaximum());
            b->Draw("hist");
            const double ym = a->GetMaximum();
            TLine* l1 = new TLine(0.119, 0, 0.119, ym); l1->SetLineColor(kBlue + 1); l1->SetLineWidth(2); l1->Draw();
            TLine* l2 = new TLine(0.173, 0, 0.173, ym); l2->SetLineColor(kBlue + 1); l2->SetLineWidth(2); l2->Draw();
            TLine* l3 = new TLine(num("tagmatched 0", 1), 0, num("tagmatched 0", 1), ym); l3->SetLineColor(kRed + 1); l3->SetLineStyle(2); l3->SetLineWidth(2); l3->Draw();
            TLine* l4 = new TLine(num("tagmatched 0", 2), 0, num("tagmatched 0", 2), ym); l4->SetLineColor(kRed + 1); l4->SetLineStyle(2); l4->SetLineWidth(2); l4->Draw();
            TLine* l5 = new TLine(0.1459, 0, 0.1459, 0.5 * ym); l5->SetLineColor(kBlack); l5->SetLineStyle(3); l5->SetLineWidth(3); l5->Draw();
            TLegend* lg = new TLegend(0.48, 0.55, 0.97, 0.88); lg->SetTextSize(0.038);
            lg->AddEntry(b, "sim bump C2 + C3, good photon", "f");
            lg->AddEntry(l1, "data #pi^{0} window 0.119-0.173", "l"); lg->AddEntry(l3, "matched sim window", "l"); lg->AddEntry(l5, "data #pi^{0} position 0.1459", "l");
            lg->Draw();
            Txt(0.20, 0.86, Form("sim #pi^{0}: %.4f, #sigma %.1f MeV/c^{2}", num("pi0sim", 0), 1e3 * num("pi0sim", 1)), 0.042);
        }
        Save(t, Form("%s/wp6_tag_cnt", kFig));
        ff->Close();
    }
    fclose(gRes); gRes = nullptr;
    P("\nwritten %s/results.txt, wp6_ana.log; %s/wp6_tableE3.{txt,tex}; %s/wp6_{E4,edge,tag}.{pdf,png}; %.0f s\n", kOut, kTab, kFig, total.RealTime());
    fclose(gLog); gLog = nullptr;
    fh->Close();
}
