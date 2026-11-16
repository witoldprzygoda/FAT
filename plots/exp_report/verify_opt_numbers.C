// ============================================================================
// verify_opt_numbers.C -- independent re-derivation of the key numbers of the
//                         chapter "Optimisation of the selection" (prefix opt)
// ============================================================================
// HADES p+p, T = 4.53 GeV (feb22), channel pi+ pi- e+ e-. EXPERIMENTAL DATA ONLY.
// Written by the verifier of the opt chapter; it does NOT use opt_common.h.
//
// What it does
//   Part A (no ScanEngine code at all): reads the compact rows with plain
//     TTree branches, applies the selections written out here from their
//     definitions, histograms M4 (2.5 MeV, 0.35-0.95) for e+e-, e+e+, e-e-,
//     builds CB = 2 sqrt(N++ N--) in 10 MeV bins, interpolates it linearly to
//     the 2.5 MeV bins (own code), counts N and CB in the windows (exact row
//     masses, event-group variance), and fits N - CB with a Gaussian (0.5439,
//     5.4 MeV, content over each bin) + pol2 in 0.5225-0.6100 by generalised
//     least squares (own code: expected-count variances iterated + covariance
//     of the interpolated CB), and also with a plain Neyman chi2 (observed
//     variances, no CB covariance) as a method cross-check.
//     Also: recomputes dz_pi = max(|pip_z - vz|, |pim_z - vz|) and compares it
//     with the stored column.
//   Part B (ScanEngine fits, but with selections written here as lambdas, not
//     opt::Cuts): FillMany + FitPeakMain + PeakModelSpread on the same
//     selections, the edge-step sum E (four M(pi+pi-) slices) for pre, PEAK and
//     bump nominal + vertex, and the claimed best-F grid point.
//
// Selections (all: cor_ok == 1, isBest >= 0; FAT already applied vz > -500, PT3)
//   pre        vz > -160
//   DIRECT     pre + dz_pi < 40
//   PEAK       pre + dz_pi < 40 + MM4 < 2.85
//   PNV        pre + oa_lab < 55 + M(pipi) > 0.30 + pip/pim_isUsedVertex
//   X17        pre + pip/pim_isUsedVertex + 15.5 <= 1000 M(ee) < 24.5
//   + isBest == 1 and OA > 2 variants of pre / PEAK, the vz < -160 rows, the
//     best-F point (oa_lab < 50, M(pipi) > 0.30, dz_pi < 30)
//
// Input : /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/scan_v2/scan_compact.root
// Output: stdout and
//         /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/opt/verify/verify_opt_numbers.txt
// Run (from /home/przygoda/HADES/FAT/FAT):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/aclic_verify_opt", kTRUE)' 'plots/exp_report/verify_opt_numbers.C+'
// ============================================================================
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

#include "TFile.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TDecompChol.h"
#include "TMath.h"
#include "TSystem.h"
#include "TTree.h"
#include "TVectorD.h"

#include "../ScanEngine.h"

namespace vopt {

const char* const kIn = "/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/scan_v2/scan_compact.root";
const char* const kOutDir = "/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/opt/verify";

struct R {   // the few columns needed
    int cor_ok, isBest, grp, fold, pipV, pimV;
    float m4, mee, mpipi, oa, oa_lab, oa_eta, mm4, vz, dz_pi, pip_z, pim_z;
};
typedef std::function<bool(const R&)> VSel;

const int NB = 240; const double LO = 0.35, HI = 0.95, BW = 0.0025; const int REB = 4, NW = NB / REB;
const double MU = 0.5439, SG = 0.0054, FLO = 0.5225, FHI = 0.6100;
const double PLO = 0.53364, PHI = 0.55416;          // data peak window
const double SLO = 0.5333, SHI = 0.5587;            // simulation peak window
const double XLO = 0.5333 * 0.5439 / 0.5460, XHI = 0.5587 * 0.5439 / 0.5460;  // shifted

struct Win { double lo, hi, n = 0, var = 0, nls[2] = {0, 0}; int cur = 0; };
struct Acc {
    std::string name; VSel sel;
    std::vector<double> h[3];
    std::vector<Win> w;      // windows (us) ; like-sign counts kept too
    long long nsel = 0;
};

double GFrac(double a, double b) {
    const double s = 1.0 / (std::sqrt(2.0) * SG);
    return 0.5 * (std::erf((b - MU) * s) - std::erf((a - MU) * s));
}

struct MyFit { double S = 0, eS = 0, chi2 = 0; int ndf = 0; double B = 0, Swin = 0, N = 0, eN = 0, CB = 0; double par[4]; };

// CB: wide bins, interpolation weights
void CBbuild(const Acc& a, std::vector<double>& cbW, std::vector<double>& cbWV, std::vector<double>& cb,
             std::vector<double>& cbV, std::vector<int>& j0, std::vector<int>& j1, std::vector<double>& a0,
             std::vector<double>& a1) {
    cbW.assign(NW, 0); cbWV.assign(NW, 0);
    std::vector<double> p(NW, 0), m(NW, 0);
    for (int i = 0; i < NB; ++i) { p[i / REB] += a.h[1][i]; m[i / REB] += a.h[2][i]; }
    for (int j = 0; j < NW; ++j) { cbW[j] = 2 * std::sqrt(p[j] * m[j]); cbWV[j] = p[j] + m[j]; }
    cb.assign(NB, 0); cbV.assign(NB, 0); j0.assign(NB, 0); j1.assign(NB, 0); a0.assign(NB, 0); a1.assign(NB, 0);
    for (int i = 0; i < NB; ++i) {
        const double x = LO + (i + 0.5) * BW;
        // centres of wide bins: LO + (j + 0.5) * REB * BW
        double t = (x - LO) / (REB * BW) - 0.5;
        int k0, k1; double f;
        if (t <= 0) { k0 = k1 = 0; f = 0; }
        else if (t >= NW - 1) { k0 = k1 = NW - 1; f = 0; }
        else { k0 = (int)std::floor(t); k1 = k0 + 1; f = t - k0; }
        j0[i] = k0; j1[i] = k1; a0[i] = (1 - f) / REB; a1[i] = f / REB;
        cb[i] = a0[i] * cbW[k0] + a1[i] * cbW[k1];
        cbV[i] = a0[i] * a0[i] * cbWV[k0] + a1[i] * a1[i] * cbWV[k1];
        if (k0 == k1) cbV[i] = (a0[i] + a1[i]) * (a0[i] + a1[i]) * cbWV[k0];
    }
}
double CBwin(const std::vector<double>& cb, double a, double b) {
    double s = 0;
    for (int i = 0; i < NB; ++i) {
        const double x1 = LO + i * BW, x2 = x1 + BW;
        const double o = std::min(x2, b) - std::max(x1, a);
        if (o > 0) s += cb[i] * o / BW;
    }
    return s;
}

// GLS fit of N - CB; mode 0: expected variances + CB covariance (iterated); 1: Neyman diag
MyFit Fit(const Acc& a, int mode, double kcb = 1.0) {
    std::vector<double> cbW, cbWV, cb, cbV, A0, A1; std::vector<int> J0, J1;
    CBbuild(a, cbW, cbWV, cb, cbV, J0, J1, A0, A1);
    for (auto& v : cb) v *= kcb;
    std::vector<int> idx;
    for (int i = 0; i < NB; ++i) { const double x = LO + (i + 0.5) * BW; if (x >= FLO && x < FHI) idx.push_back(i); }
    const int n = idx.size(), np = 4;
    const double x0 = 0.5 * (FLO + FHI), dx = 0.5 * (FHI - FLO);
    TMatrixD F(n, np); std::vector<double> y(n), c(n), var(n);
    for (int q = 0; q < n; ++q) {
        const int i = idx[q]; const double x1 = LO + i * BW, xc = x1 + 0.5 * BW, t = (xc - x0) / dx;
        F(q, 0) = GFrac(x1, x1 + BW); F(q, 1) = 1; F(q, 2) = t; F(q, 3) = t * t;
        y[q] = a.h[0][i] - cb[i]; c[q] = cb[i];
        var[q] = mode == 1 ? std::max(a.h[0][i] + kcb * kcb * cbV[i], 1.0) : std::max(a.h[0][i], 1.0);
    }
    TMatrixDSym C0(n);
    if (mode == 0)
        for (int q = 0; q < n; ++q) for (int p = 0; p < n; ++p) {
            const int i = idx[q], k = idx[p];
            double s = 0;
            const int ji[2] = {J0[i], J1[i]}, jk[2] = {J0[k], J1[k]};
            const double wi[2] = {A0[i], A1[i]}, wk[2] = {A0[k], A1[k]};
            for (int u = 0; u < 2; ++u) for (int v = 0; v < 2; ++v) if (ji[u] == jk[v]) s += wi[u] * wk[v] * cbWV[ji[u]];
            C0(q, p) = kcb * kcb * s;
        }
    TVectorD beta(np); TMatrixD covb(np, np); double chi2 = 0;
    for (int it = 0; it < 12; ++it) {
        TMatrixDSym V(n);
        if (mode == 0) V = C0; else V.Zero();
        for (int q = 0; q < n; ++q) V(q, q) += var[q];
        TMatrixDSym Vi(V); Vi.Invert();
        TMatrixD Ft(TMatrixD::kTransposed, F);
        TMatrixD A = Ft * Vi * F; A.Invert();
        TVectorD yy(n); for (int q = 0; q < n; ++q) yy[q] = y[q];
        beta = A * (Ft * (Vi * yy));
        covb = A;
        TVectorD r = yy - F * beta;
        chi2 = r * (Vi * r);
        if (mode == 1) break;
        for (int q = 0; q < n; ++q) {
            double f = 0; for (int k = 0; k < np; ++k) f += F(q, k) * beta[k];
            var[q] = std::max(f + c[q], 1.0);
        }
    }
    MyFit r; r.S = beta[0]; r.eS = std::sqrt(covb(0, 0)); r.chi2 = chi2; r.ndf = n - np;
    for (int k = 0; k < 4; ++k) r.par[k] = beta[k];
    r.Swin = r.S * GFrac(PLO, PHI);
    const int nst = 400; const double hh = (PHI - PLO) / nst; double s = 0;
    for (int q = 0; q < nst; ++q) { const double t = (PLO + (q + 0.5) * hh - x0) / dx; s += beta[1] + beta[2] * t + beta[3] * t * t; }
    r.B = s * hh / BW;
    r.CB = CBwin(cb, PLO, PHI);
    return r;
}

FILE* gF = nullptr;
void P(const char* fmt, ...) {
    va_list a, b; va_start(a, fmt); va_copy(b, a);
    vprintf(fmt, a); if (gF) vfprintf(gF, fmt, b);
    va_end(a); va_end(b);
}

} // namespace vopt

void verify_opt_numbers() {
    using namespace vopt;
    gSystem->mkdir(kOutDir, kTRUE);
    gF = fopen(Form("%s/verify_opt_numbers.txt", kOutDir), "w");
    P("# verify_opt_numbers.C: independent check of the opt chapter (input %s)\n", kIn);

    // ---------------- selections (written from the definitions) ----------------
    auto pre  = [](const R& r) { return r.cor_ok == 1 && r.isBest >= 0 && r.vz > -160.f; };
    auto vtx2 = [](const R& r) { return r.pipV == 1 && r.pimV == 1; };
    std::vector<Acc> acc;
    auto add = [&](const char* n, VSel s) { Acc a; a.name = n; a.sel = s; acc.push_back(a); };
    add("pre",                    [=](const R& r) { return pre(r); });
    add("DIRECT",                 [=](const R& r) { return pre(r) && r.dz_pi < 40.f; });
    add("PEAK",                   [=](const R& r) { return pre(r) && r.dz_pi < 40.f && r.mm4 < 2.85f; });
    add("peak nominal + vertex",  [=](const R& r) { return pre(r) && r.oa_lab < 55.f && r.mpipi > 0.30f && vtx2(r); });
    add("X17 band",               [=](const R& r) { return pre(r) && vtx2(r) && r.mee * 1000.f >= 15.5f && r.mee * 1000.f < 24.5f; });
    add("pre isBest==1",          [=](const R& r) { return pre(r) && r.isBest == 1; });
    add("PEAK isBest==1",         [=](const R& r) { return pre(r) && r.isBest == 1 && r.dz_pi < 40.f && r.mm4 < 2.85f; });
    add("PEAK OA>2",              [=](const R& r) { return pre(r) && r.oa > 2.f && r.dz_pi < 40.f && r.mm4 < 2.85f; });
    add("pre OA>2",               [=](const R& r) { return pre(r) && r.oa > 2.f; });
    add("pre without vz",         [=](const R& r) { return r.cor_ok == 1 && r.isBest >= 0; });
    add("vz <= -160 rows",        [=](const R& r) { return r.cor_ok == 1 && r.isBest >= 0 && !(r.vz > -160.f); });
    add("best-F point",           [=](const R& r) { return pre(r) && r.oa_lab < 50.f && r.mpipi > 0.30f && r.dz_pi < 30.f; });
    for (auto& a : acc) {
        for (int k = 0; k < 3; ++k) a.h[k].assign(NB, 0.);
        double wl[3][2] = {{PLO, PHI}, {SLO, SHI}, {XLO, XHI}};
        for (int w = 0; w < 3; ++w) { Win x; x.lo = wl[w][0]; x.hi = wl[w][1]; a.w.push_back(x); }
    }

    // ---------------- Part A: read the rows ----------------
    TFile* f = TFile::Open(kIn, "READ");
    const char* tn[3] = {"us", "pp", "mm"};
    long long nDz = 0, nDzBad = 0, nDzNaN = 0;
    for (int k = 0; k < 3; ++k) {
        TTree* t = (TTree*)f->Get(tn[k]);
        R r;
        t->SetBranchStatus("*", 0);
        const char* ib[] = {"cor_ok", "isBest", "grp", "fold", "pip_isUsedVertex", "pim_isUsedVertex"};
        int* ip[] = {&r.cor_ok, &r.isBest, &r.grp, &r.fold, &r.pipV, &r.pimV};
        for (int q = 0; q < 6; ++q) { t->SetBranchStatus(ib[q], 1); t->SetBranchAddress(ib[q], ip[q]); }
        const char* fb[] = {"m4", "mee", "mpipi", "oa", "oa_lab", "oa_eta", "mm4", "vz", "dz_pi", "pip_z", "pim_z"};
        float* fp[] = {&r.m4, &r.mee, &r.mpipi, &r.oa, &r.oa_lab, &r.oa_eta, &r.mm4, &r.vz, &r.dz_pi, &r.pip_z, &r.pim_z};
        for (int q = 0; q < 11; ++q) { t->SetBranchStatus(fb[q], 1); t->SetBranchAddress(fb[q], fp[q]); }
        const Long64_t n = t->GetEntries();
        int cur = -1;
        std::vector<std::vector<int>> kc(acc.size(), std::vector<int>(3, 0));
        auto flush = [&]() {
            for (size_t s = 0; s < acc.size(); ++s)
                for (int w = 0; w < 3; ++w) { if (k == 0) acc[s].w[w].var += (double)kc[s][w] * kc[s][w]; kc[s][w] = 0; }
        };
        for (Long64_t i = 0; i < n; ++i) {
            t->GetEntry(i);
            if (r.grp != cur) { flush(); cur = r.grp; }
            // dz_pi recomputation (float arithmetic as documented)
            if (k == 0) {
                const float d = std::max(std::fabs(r.pip_z - r.vz), std::fabs(r.pim_z - r.vz));
                ++nDz;
                if (!(d == d) || !(r.dz_pi == r.dz_pi)) { if (!((d == d) == (r.dz_pi == r.dz_pi))) ++nDzBad; else ++nDzNaN; }
                else if (d != r.dz_pi) ++nDzBad;
            }
            const double m = r.m4;
            if (!(m == m)) continue;
            const int b = (m >= LO && m < HI) ? (int)((m - LO) / BW) : -1;
            for (size_t s = 0; s < acc.size(); ++s) {
                if (!acc[s].sel(r)) continue;
                if (k == 0) acc[s].nsel++;
                if (b >= 0) acc[s].h[k][b] += 1;
                for (int w = 0; w < 3; ++w)
                    if (m >= acc[s].w[w].lo && m < acc[s].w[w].hi) {
                        if (k == 0) { acc[s].w[w].n += 1; kc[s][w]++; }
                        else acc[s].w[w].nls[k - 1] += 1;
                    }
            }
        }
        flush();
    }
    f->Close();
    P("\n[dz_pi] e+e- rows %lld: recomputed max(|pip_z-vz|,|pim_z-vz|) differs from the stored dz_pi in %lld rows; both NaN in %lld rows\n",
      nDz, nDzBad, nDzNaN);

    P("\n[Part A] own CB + own GLS fit (mode 0: expected variances + CB covariance; mode 1: Neyman diag)\n");
    P("%-24s | %7s %9s %6s | %8s %6s %6s %5s | %7s %5s | %8s %6s | %6s %6s | %s\n", "selection", "N_win", "+-grp", "CB/N",
      "S(0)", "eS", "F", "chi2", "S/B", "", "S(1)", "eS", "S(CB-3%)", "S(+3%)", "sim win N, CB ; shifted win N, CB");
    for (auto& a : acc) {
        const MyFit f0 = Fit(a, 0), f1 = Fit(a, 1), fm = Fit(a, 0, 0.97), fp = Fit(a, 0, 1.03);
        std::vector<double> cbW, cbWV, cb, cbV, A0, A1; std::vector<int> J0, J1;
        CBbuild(a, cbW, cbWV, cb, cbV, J0, J1, A0, A1);
        const double cbS = CBwin(cb, SLO, SHI), cbX = CBwin(cb, XLO, XHI);
        P("%-24s | %7.0f %9.1f %6.3f | %8.1f %6.1f %6.2f %5.1f | %7.4f %5s | %8.1f %6.1f | %6.0f %6.0f | %6.0f +- %4.1f, %5.1f ; %6.0f +- %4.1f, %5.1f\n",
          a.name.c_str(), a.w[0].n, std::sqrt(a.w[0].var), f0.CB / a.w[0].n, f0.S, f0.eS, f0.S / f0.eS, f0.chi2,
          f0.Swin / (f0.B + f0.CB), "", f1.S, f1.eS, fm.S, fp.S, a.w[1].n, std::sqrt(a.w[1].var), cbS, a.w[2].n,
          std::sqrt(a.w[2].var), cbX);
    }

    // ---------------- Part B: engine fits on the same (lambda) selections ----------------
    scan::Data d;
    scan::Load(d);
    using scan::Row;
    auto preE = [](const Row& r) { return r.cor_ok == 1 && r.isBest >= 0 && r.vz > -160.f; };
    auto v2E = [](const Row& r) { return r.pip_isUsedVertex == 1 && r.pim_isUsedVertex == 1; };
    std::vector<std::string> names;
    std::vector<scan::Sel> sels;
    auto addE = [&](const char* n, scan::Sel s) { names.push_back(n); sels.push_back(s); };
    addE("pre",                   [=](const Row& r) { return preE(r); });
    addE("DIRECT",                [=](const Row& r) { return preE(r) && r.dz_pi < 40.f; });
    addE("PEAK",                  [=](const Row& r) { return preE(r) && r.dz_pi < 40.f && r.mm4 < 2.85f; });
    addE("peak nominal + vertex", [=](const Row& r) { return preE(r) && r.oa_lab < 55.f && r.mpipi > 0.30f && v2E(r); });
    addE("X17 band",              [=](const Row& r) { return preE(r) && v2E(r) && r.mee * 1000.f >= 15.5f && r.mee * 1000.f < 24.5f; });
    addE("pre isBest==1",         [=](const Row& r) { return preE(r) && r.isBest == 1; });
    addE("PEAK isBest==1",        [=](const Row& r) { return preE(r) && r.isBest == 1 && r.dz_pi < 40.f && r.mm4 < 2.85f; });
    addE("PEAK OA>2",             [=](const Row& r) { return preE(r) && r.oa > 2.f && r.dz_pi < 40.f && r.mm4 < 2.85f; });
    addE("pre OA>2",              [=](const Row& r) { return preE(r) && r.oa > 2.f; });
    addE("best-F point",          [=](const Row& r) { return preE(r) && r.oa_lab < 50.f && r.mpipi > 0.30f && r.dz_pi < 30.f; });
    addE("vz <= -160 rows",       [=](const Row& r) { return r.cor_ok == 1 && r.isBest >= 0 && !(r.vz > -160.f); });
    addE("bump nominal + vertex", [=](const Row& r) { return preE(r) && r.oa_lab < 40.f && r.mpipi > 0.30f && r.mpipi < 0.41f && v2E(r); });
    const std::vector<scan::Spectra> sp = scan::FillMany(d, sels);
    P("\n[Part B] ScanEngine FitPeakMain + PeakModelSpread on lambda selections\n");
    P("%-24s | %7s %6s | %8s %6s %7s(%2s) | %6s %7s | %s\n", "selection", "N_win", "CB/N", "S", "stat", "spread", "n", "F", "S/B",
      "S x 64/63 / 6079.2 [nb]");
    for (size_t i = 0; i < sels.size(); ++i) {
        const scan::PeakFit pf = scan::FitPeakMain(sp[i]);
        const scan::Spread spr = scan::PeakModelSpread(sp[i]);
        P("%-24s | %7.0f %6.3f | %8.1f %6.1f %7.1f(%2d) | %6.2f %7.4f | %.4f +- %.4f +- %.4f\n", names[i].c_str(), pf.N, pf.CB / pf.N,
          pf.S, pf.eS, spr.half, spr.npass, pf.F, pf.SoB, pf.S * 64. / 63. / 6079.2, pf.eS * 64. / 63. / 6079.2,
          spr.half * 64. / 63. / 6079.2);
    }

    // edge-step sum E (procedure of scan_bump.C; slices 0.30-0.34, 0.34-0.37, 0.37-0.39, 0.39-0.40)
    const double sl[5] = {0.30, 0.34, 0.37, 0.39, 0.40};
    std::vector<scan::Sel> es;
    const int ie[3] = {0, 2, 11};
    for (int q = 0; q < 3; ++q)
        for (int i = 0; i < 4; ++i) {
            const scan::Sel s = sels[ie[q]]; const float a = sl[i], b = sl[i + 1];
            es.push_back([=](const Row& r) { return r.mpipi >= a && r.mpipi < b && s(r); });
        }
    const std::vector<scan::Spectra> spe = scan::FillMany(d, es);
    P("\n[Part B] edge-step sum E (rows per 2.5 MeV)\n");
    for (int q = 0; q < 3; ++q) {
        double E = 0, v = 0;
        for (int i = 0; i < 4; ++i) {
            scan::StepCfg c; c.mpLo = sl[i]; c.mpHi = sl[i + 1];
            if (i == 3) c.halfRange = 0.035;
            const scan::StepFit ft = scan::FitEdgeStep(spe[q * 4 + i], c);
            if (!ft.ok) continue;
            E += ft.h; v += ft.eh * ft.eh;
        }
        P("%-24s E = %.1f +- %.1f ; N_b (bump window) = %.0f, W = N - CB = %.0f\n", names[ie[q]].c_str(), E, std::sqrt(v),
          sp[ie[q]].win[1].n[0], sp[ie[q]].win[1].sig());
    }
    if (gF) fclose(gF);
}
