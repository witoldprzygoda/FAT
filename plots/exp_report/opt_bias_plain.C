// opt_bias_plain.C -- COPY of plots/exp_report/opt_bias.C made by the short-report figure agent (integ6, 2026-11-16):
// internal catalogue labels in the figure text replaced by plain names (glossary_v2.md); outputs redirected:
// figures of the scope -> docs/latex/<dir>/<base>_plain.{pdf,png}; everything else -> /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ6/shortfig/run/opt
// Built by /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ6/shortfig/mk_plain.py; the physics code is unchanged.
// ============================================================================
// opt_bias.C  --  optimisation bias of the data-driven peak selection, measured
//                 with independent halves of the events (fold and 40 random halvings)
// ============================================================================
// HADES p+p, T = 4.53 GeV (feb22), channel pi+ pi- e+ e-. EXPERIMENTAL DATA ONLY.
//
// Question: when a selection is chosen on the data by maximising F = S/sigma_S,
// how much of its F is a statistical up-fluctuation that the choice picked up?
// And does the chosen selection beat "pre" on events it was NOT chosen on?
//
// What it does
//   1. grid of working points on top of pre (cor_ok && isBest >= 0 && vz > -160):
//        oa_lab   < 40, 45, 50, 55, 60, 65, no cut
//        M(pipi)  > 0.34, 0.32, 0.30, no cut
//        pions in the event-vertex fit: both, at least one, no cut
//        dz_pi    < 20, 30, 40, 60 mm, no cut
//        MM4      < 2.85 GeV/c2, no cut  (MM4 < 2.85 is flagged by the sculpting
//                 test, opt_sculpt.txt: never chosen by the rules, only used for the
//                 fixed point PEAK)
//      = 840 points. Every row gets one "type" (the tightest value of each
//      variable it passes); the M4 spectra of a grid point are prefix sums over
//      the types, so all points of all halvings come from one pass over the rows.
//   2. halvings of the event groups: the parity of the group index (fold, as in
//      scan_peak.C) and 40 random halvings (hash of the group index, the same hash
//      as scan_peak.C, seeds 1..40);
//   3. on each half A the working point is chosen by four fixed rules
//        best-F  : largest F_A
//        loose   : most rows in the window among points with F_A >= F_A,best - 1
//        pure    : largest S/B among points with F_A >= F_A,best - 1
//        nominal : the self-consistent plateau centre of scan_peak.C (N-1 scans
//                  around the current point; each variable goes to the centre of
//                  its values with F >= F_best - 1, or stays open; start oa_lab < 50;
//                  repeated until stable, at most 8 iterations)
//      and evaluated on the other half B;
//   4. bias estimator (scan_peak.C): b = 1/2 [ (F_A(c_A) - F_B(c_A)) + (F_B(c_B) - F_A(c_B)) ];
//      the plain difference of the halves cancels; mean +- rms/sqrt(n) over the
//      halvings (the halvings reuse the same events: the error is too small by an
//      unknown factor); fixed points (pre, peak nominal, peak nominal + vertex,
//      DIRECT, PEAK) give b = 0 by construction (check);
//   5. independent gain: F_B(c_A) - F_B(pre), averaged (the improvement over pre on
//      events that did not take part in the choice), compared with the same for the
//      frozen DIRECT and PEAK (chosen on simulation: no data bias).
//   F of one half is about F(full)/sqrt(2); numbers are given per half and scaled
//   by sqrt(2) to the full sample.
//   P_main as everywhere (N - CB, frozen Gaussian 0.5439 GeV/c2 / 5.4 MeV + pol2,
//   0.5225-0.6100, CB = 2 sqrt(N++ N--) per 10 MeV, k = 1); rows of 0.50-0.63 only.
//
// Inputs : /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/scan_v2/scan_compact.root (scan::kCompact)
//          first-study numbers quoted from /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/scan/peak_bias.txt
// Outputs: /home/przygoda/HADES/FAT/FAT/docs/latex/tables_exp/opt_bias.txt, opt_bias.tex
//          /home/przygoda/HADES/FAT/FAT/docs/latex/figs_exp/opt_bias.pdf (+ .png)
//
// Verifier changes (2026-11-14): legend of panel (b) moved above the data (y range
// to 3.2); the LaTeX caption names the MM4 grid dimension and explains that b = 0
// of the fixed points does not make the points chosen on these data unbiased.
//
// Run (from the repository root /home/przygoda/HADES/FAT/FAT), about 2 min:
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/aclic_opt_bias", kTRUE)' 'plots/exp_report/opt_bias.C+'
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TColor.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TSystem.h"

#include "opt_common_plain.h"

namespace optbias {

using namespace opt;

// ---- grid (ordered from the tightest value, index 0, to "no cut") ---------------
const int kNA = 7, kNM = 4, kNV = 3, kND = 5, kNX = 2;
const float kOaLabV[kNA] = {40, 45, 50, 55, 60, 65, 1e30f};
const float kMpLoV[kNM]  = {0.34f, 0.32f, 0.30f, -1e30f};
const int   kVtxV[kNV]   = {2, 1, 0};
const float kDzV[kND]    = {20, 30, 40, 60, 1e30f};
const float kMm4V[kNX]   = {2.85f, 1e30f};
const int kNT = kNA * kNM * kNV * kND * kNX;
inline int TIdx(int a, int m, int v, int d, int x) { return (((a * kNM + m) * kNV + v) * kND + d) * kNX + x; }
inline void TDec(int t, int& a, int& m, int& v, int& d, int& x) {
    x = t % kNX; t /= kNX; d = t % kND; t /= kND; v = t % kNV; t /= kNV; m = t % kNM; a = t / kNM;
}
std::string PtName(int g) {
    int a, m, v, d, x; TDec(g, a, m, v, d, x);
    std::string s;
    auto add = [&](const std::string& t) { if (!s.empty()) s += ", "; s += t; };
    if (a < kNA - 1) add(Form("oa_lab<%.0f", kOaLabV[a]));
    if (m < kNM - 1) add(Form("Mpipi>%.2f", kMpLoV[m]));
    if (v == 0) add("both pi in vtx"); if (v == 1) add(">=1 pi in vtx");
    if (d < kND - 1) add(Form("dz_pi<%.0f", kDzV[d]));
    if (x == 0) add("MM4<2.85");
    return s.empty() ? "pre" : s;
}
// type of a row: tightest value passed per variable
inline int RowType(const Row& r) {
    int a = kNA - 1; for (int i = 0; i < kNA - 1; ++i) if (r.oa_lab < kOaLabV[i]) { a = i; break; }
    int m = kNM - 1; for (int i = 0; i < kNM - 1; ++i) if (r.mpipi > kMpLoV[i]) { m = i; break; }
    const int nv = (r.pip_isUsedVertex == 1) + (r.pim_isUsedVertex == 1);
    const int v = nv == 2 ? 0 : (nv == 1 ? 1 : 2);
    int d = kND - 1; for (int i = 0; i < kND - 1; ++i) if (r.dz_pi < kDzV[i]) { d = i; break; }
    const int x = r.mm4 < kMm4V[0] ? 0 : 1;
    return TIdx(a, m, v, d, x);
}
// the same hash as scan_peak.C InSamp(mode 2)
inline int RndHalf(int grp, unsigned seed) {
    unsigned x = (unsigned)grp * 2654435761u + seed * 0x9E3779B9u;
    x ^= x >> 16; x *= 0x85ebca6bu; x ^= x >> 13; x *= 0xc2b2ae35u; x ^= x >> 16;
    return (int)(x & 1u);
}

// ---- spectra from bin counts (the CB part of scan::FillMany) -------------------
const double kRLo = 0.50, kRHi = 0.63;
const int kB0 = (int)std::lround((kRLo - scan::kMlo) / scan::kBinW);   // first engine bin used
const int kB1 = (int)std::lround((kRHi - scan::kMlo) / scan::kBinW);   // one past the last
const int kNB = kB1 - kB0;
Spectra MakeSpectra(const double* n0, const double* n1, const double* n2) {
    Spectra sp;
    sp.nbin = scan::kNbin; sp.reb = scan::kCBreb; sp.nw = scan::kNbin / scan::kCBreb;
    sp.lo = scan::kMlo; sp.hi = scan::kMhi; sp.bw = scan::kBinW;
    for (int k = 0; k < 3; ++k) { sp.n[k].assign(sp.nbin, 0.); sp.w[k].assign(sp.nw, 0.); }
    for (int i = 0; i < kNB; ++i) { sp.n[0][kB0 + i] = n0[i]; sp.n[1][kB0 + i] = n1[i]; sp.n[2][kB0 + i] = n2[i]; }
    for (int k = 0; k < 3; ++k) for (int i = 0; i < sp.nbin; ++i) sp.w[k][i / sp.reb] += sp.n[k][i];
    sp.cbW.assign(sp.nw, 0.); sp.cbWVar.assign(sp.nw, 0.);
    for (int j = 0; j < sp.nw; ++j) { sp.cbW[j] = 2.0 * std::sqrt(sp.w[1][j] * sp.w[2][j]); sp.cbWVar[j] = sp.w[1][j] + sp.w[2][j]; }
    sp.cb.assign(sp.nbin, 0.); sp.cbVar.assign(sp.nbin, 0.); sp.sig.assign(sp.nbin, 0.); sp.sigVar.assign(sp.nbin, 0.);
    for (int i = 0; i < sp.nbin; ++i) {
        int j0, j1; double f; sp.Interp(i, j0, j1, f);
        const double a0 = (1 - f) / sp.reb, a1 = f / sp.reb;
        sp.cb[i] = a0 * sp.cbW[j0] + a1 * sp.cbW[j1];
        sp.cbVar[i] = a0 * a0 * sp.cbWVar[j0] + a1 * a1 * sp.cbWVar[j1];
        sp.sig[i] = sp.n[0][i] - sp.cb[i];
        sp.sigVar[i] = sp.n[0][i] + sp.cbVar[i];
    }
    return sp;
}

struct PtVal { double F = 0, S = 0, eS = 0, N = 0, SoB = 0; };

// ---- plateau of a 1D scan (copy of scan_peak.C Plateau1D, loose = high index) ------
struct Plat { int lo = 0, hi = 0, ctr = 0, best = 0; double fmax = 0; };
Plat Plateau1D(const std::vector<double>& F, const std::vector<int>& flag) {
    const int n = (int)F.size();
    Plat p; p.best = -1;
    for (int i = 0; i < n; ++i) if (!flag[i] && (p.best < 0 || F[i] > F[p.best])) p.best = i;
    if (p.best < 0) { p.best = 0; return p; }
    p.fmax = F[p.best]; p.lo = p.hi = p.best;
    while (p.lo > 0 && !flag[p.lo - 1] && F[p.lo - 1] >= p.fmax - 1.0) p.lo--;
    while (p.hi < n - 1 && !flag[p.hi + 1] && F[p.hi + 1] >= p.fmax - 1.0) p.hi++;
    while (p.hi - p.lo < 2) {
        const bool cl = p.lo > 0 && !flag[p.lo - 1], ch = p.hi < n - 1 && !flag[p.hi + 1];
        if (!cl && !ch) break;
        if (cl && (!ch || F[p.lo - 1] >= F[p.hi + 1])) p.lo--; else p.hi++;
    }
    const int s = p.lo + p.hi;
    p.ctr = (s % 2 == 0) ? s / 2 : (s + 1) / 2;     // half-integer middle -> the looser value
    return p;
}

} // namespace optbias

// ============================================================================
void opt_bias_plain(int nSplit = 40) {
    using namespace opt;
    using namespace optbias;
    TStopwatch sw;
    gSystem->mkdir(kFigOut, kTRUE); gSystem->mkdir(kTabOut, kTRUE); gSystem->mkdir(kPlainScratch, kTRUE);
    Style();
    scan::Data d;
    if (!scan::Load(d, scan::kCompact, kRLo - 0.001, kRHi + 0.001)) return;
    const Cuts pre = Pre();
    const int nH = nSplit + 1;                       // 0 = fold, 1..nSplit = random

    // ---- counts per (halving, half, type, sample, bin) ----
    std::vector<double> cnt((size_t)nH * 2 * kNT * 3 * kNB, 0.);
    auto at = [&](int h, int s, int t, int k, int b) -> double& { return cnt[((((size_t)h * 2 + s) * kNT + t) * 3 + k) * kNB + b]; };
    long nrow = 0;
    for (int k = 0; k < 3; ++k)
        for (const Row& r : d.s[k]) {
            if (!Pass(pre, r)) continue;
            const double m = r.m4;
            if (!(m >= kRLo && m < kRHi)) continue;
            const int ib = (int)((m - scan::kMlo) / scan::kBinW) - kB0;
            if (ib < 0 || ib >= kNB) continue;
            const int t = RowType(r);
            at(0, r.fold, t, k, ib) += 1;
            for (int h = 1; h < nH; ++h) at(h, RndHalf(r.grp, (unsigned)h), t, k, ib) += 1;
            nrow++;
        }
    printf("rows filled: %ld at %.1f s\n", nrow, sw.RealTime()); sw.Continue();
    // prefix sums over the five type indices (in place): count of grid point g = sum of types t <= g
    for (int h = 0; h < nH; ++h) for (int s = 0; s < 2; ++s) for (int k = 0; k < 3; ++k) for (int b = 0; b < kNB; ++b) {
        for (int dim = 0; dim < 5; ++dim)
            for (int t = 0; t < kNT; ++t) {
                int idx[5]; TDec(t, idx[0], idx[1], idx[2], idx[3], idx[4]);
                if (idx[dim] == 0) continue;
                idx[dim]--;
                at(h, s, t, k, b) += at(h, s, TIdx(idx[0], idx[1], idx[2], idx[3], idx[4]), k, b);
            }
    }
    printf("prefix sums at %.1f s\n", sw.RealTime()); sw.Continue();
    // NOTE: the in-place prefix sum above walks t in increasing order, so the predecessor
    // along 'dim' has already been summed along the same dim -> a correct cumulative sum.

    // ---- F of every grid point on every half ----
    std::vector<PtVal> V((size_t)nH * 2 * kNT);
    auto val = [&](int h, int s, int g) -> PtVal& { return V[((size_t)h * 2 + s) * kNT + g]; };
    for (int h = 0; h < nH; ++h) for (int s = 0; s < 2; ++s) for (int g = 0; g < kNT; ++g) {
        const Spectra sp = MakeSpectra(&at(h, s, g, 0, 0), &at(h, s, g, 1, 0), &at(h, s, g, 2, 0));
        const PeakFit f = scan::FitPeakMain(sp);
        PtVal& v = val(h, s, g); v.F = f.F; v.S = f.S; v.eS = f.eS; v.N = f.N; v.SoB = f.SoB;
    }
    printf("fits at %.1f s\n", sw.RealTime()); sw.Continue();
    // full sample (both halves of halving 1) for the reference values
    std::vector<PtVal> Vfull(kNT);
    for (int g = 0; g < kNT; ++g) {
        std::vector<double> n[3];
        for (int k = 0; k < 3; ++k) { n[k].assign(kNB, 0.); for (int b = 0; b < kNB; ++b) n[k][b] = at(0, 0, g, k, b) + at(0, 1, g, k, b); }
        const PeakFit f = scan::FitPeakMain(MakeSpectra(n[0].data(), n[1].data(), n[2].data()));
        Vfull[g].F = f.F; Vfull[g].S = f.S; Vfull[g].eS = f.eS; Vfull[g].N = f.N; Vfull[g].SoB = f.SoB;
    }

    // ---- rules ----
    auto flagged = [](int g) { int a, m, v, dd, x; TDec(g, a, m, v, dd, x); return x == 0; };   // MM4 < 2.85 not chosen
    const int gPre = TIdx(kNA - 1, kNM - 1, kNV - 1, kND - 1, kNX - 1);
    const int gPeakNom = TIdx(3, 2, 2, kND - 1, kNX - 1);
    const int gPeakNomV = TIdx(3, 2, 0, kND - 1, kNX - 1);
    const int gDirect = TIdx(kNA - 1, kNM - 1, kNV - 1, 2, kNX - 1);
    const int gPeakF = TIdx(kNA - 1, kNM - 1, kNV - 1, 2, 0);
    const int gPreV = TIdx(kNA - 1, kNM - 1, 0, kND - 1, kNX - 1);
    const int gNomDz = TIdx(3, 2, kNV - 1, 2, kNX - 1);
    auto ruleBest = [&](int h, int s) { int b = -1; for (int g = 0; g < kNT; ++g) if (!flagged(g) && (b < 0 || val(h, s, g).F > val(h, s, b).F)) b = g; return b; };
    auto ruleLoose = [&](int h, int s) {
        const double fb = val(h, s, ruleBest(h, s)).F; int b = -1;
        for (int g = 0; g < kNT; ++g) if (!flagged(g) && val(h, s, g).F >= fb - 1.0 && (b < 0 || val(h, s, g).N > val(h, s, b).N)) b = g;
        return b; };
    auto rulePure = [&](int h, int s) {
        const double fb = val(h, s, ruleBest(h, s)).F; int b = -1;
        for (int g = 0; g < kNT; ++g) if (!flagged(g) && val(h, s, g).F >= fb - 1.0 && (b < 0 || val(h, s, g).SoB > val(h, s, b).SoB)) b = g;
        return b; };
    auto ruleNominal = [&](int h, int s) {
        int cur[5] = {2, kNM - 1, kNV - 1, kND - 1, kNX - 1};         // start: oa_lab < 50
        const int nvals[4] = {kNA, kNM, kNV, kND};
        for (int it = 0; it < 8; ++it) {
            int old[5]; for (int q = 0; q < 5; ++q) old[q] = cur[q];
            for (int dim = 0; dim < 4; ++dim) {
                std::vector<double> F; std::vector<int> fl;
                for (int i = 0; i < nvals[dim]; ++i) { int c[5]; for (int q = 0; q < 5; ++q) c[q] = cur[q]; c[dim] = i;
                    F.push_back(val(h, s, TIdx(c[0], c[1], c[2], c[3], c[4])).F); fl.push_back(0); }
                const Plat p = Plateau1D(F, fl);
                // "stays open" when no cut has the largest F
                cur[dim] = (p.best == nvals[dim] - 1) ? nvals[dim] - 1 : p.ctr;
            }
            bool same = true; for (int q = 0; q < 5; ++q) if (cur[q] != old[q]) same = false;
            if (same) break;
        }
        return TIdx(cur[0], cur[1], cur[2], cur[3], cur[4]);
    };
    const int nR = 4;
    const char* rName[nR] = {"nominal", "loose", "pure", "best-F"};
    const char* rPlain[nR] = {"plateau centre", "most counts", "largest S/B", "largest F"};
    const int nFix = 7;
    const int gFix[nFix] = {gPre, gPeakNom, gPeakNomV, gPreV, gNomDz, gDirect, gPeakF};
    const char* fName[nFix] = {"pre (fixed)", "peak nominal (fixed)", "peak nominal + vertex (fixed)", "pre + vertex (fixed)",
                               "peak nominal + dz_pi<40 (fixed)", "DIRECT (frozen)", "PEAK (frozen)"};
    // per halving: chosen points, bias, independent gain
    struct HRes { int c[2][nR]; double b[nR], gain[nR], bFix[nFix], gainFix[nFix]; };
    std::vector<HRes> H(nH);
    for (int h = 0; h < nH; ++h) {
        HRes& R = H[h];
        for (int s = 0; s < 2; ++s) { R.c[s][0] = ruleNominal(h, s); R.c[s][1] = ruleLoose(h, s); R.c[s][2] = rulePure(h, s); R.c[s][3] = ruleBest(h, s); }
        for (int q = 0; q < nR; ++q) {
            const int cA = R.c[0][q], cB = R.c[1][q];
            R.b[q] = 0.5 * ((val(h, 0, cA).F - val(h, 1, cA).F) + (val(h, 1, cB).F - val(h, 0, cB).F));
            // independent gain over pre, averaged over the two directions
            R.gain[q] = 0.5 * ((val(h, 1, cA).F - val(h, 1, gPre).F) + (val(h, 0, cB).F - val(h, 0, gPre).F));
        }
        for (int q = 0; q < nFix; ++q) {
            const int g = gFix[q];
            R.bFix[q] = 0.5 * ((val(h, 0, g).F - val(h, 1, g).F) + (val(h, 1, g).F - val(h, 0, g).F));
            R.gainFix[q] = 0.5 * ((val(h, 1, g).F - val(h, 1, gPre).F) + (val(h, 0, g).F - val(h, 0, gPre).F));
        }
    }
    auto meanErr = [&](std::function<double(const HRes&)> f, double& m, double& e) {
        double s = 0, s2 = 0; int n = 0;
        for (int h = 1; h < nH; ++h) { const double x = f(H[h]); s += x; s2 += x * x; n++; }
        m = s / n; e = std::sqrt(std::max(0., s2 / n - m * m) / n);
    };

    // ---- table ----
    Out o;
    o.Open(std::string(kTabOut) + "/opt_bias.txt");
    o.P("# opt_bias.txt  --  written by plots/exp_report/opt_bias.C\n");
    o.P("# EXPERIMENTAL DATA. Base pre = cor_ok && isBest >= 0 && vz > -160 mm (FAT: vz > -500, PT3); rows with 0.50 <= M4 < 0.63 (%ld rows).\n", nrow);
    o.P("# P_main (N - CB, frozen Gaussian 0.5439 / 5.4 MeV + pol2, 0.5225-0.6100), F = S/sigma_S statistical only.\n");
    o.P("# Grid: oa_lab < 40..65/none x M(pipi) > 0.34/0.32/0.30/none x pions in vertex fit both/>=1/none x dz_pi < 20/30/40/60/none x MM4 < 2.85/none\n");
    o.P("#   = %d points; MM4 < 2.85 is never chosen by the rules (sculpting flag), only used for the fixed PEAK point.\n", kNT);
    o.P("# Halvings: fold (parity of the event-group index) and %d random halvings (hash of the group index, seeds 1..%d, as scan_peak.C).\n", nSplit, nSplit);
    o.P("# Bias b = 1/2 [(F_A(c_A) - F_B(c_A)) + (F_B(c_B) - F_A(c_B))], c_X = point chosen on half X; F of one half ~ F(full)/sqrt(2).\n");
    o.P("# Independent gain = 1/2 [(F_B(c_A) - F_B(pre)) + (F_A(c_B) - F_A(pre))]: improvement over pre on events not used for the choice.\n");
    o.P("# mean +- rms/sqrt(n) over the random halvings (the halvings reuse the same events: errors too small by an unknown factor).\n#\n");
    o.P("FULL SAMPLE, reference points: F (S +- stat)\n");
    for (int q = 0; q < nFix; ++q) o.P("   %-34s F = %6.2f  (S = %5.0f +- %4.0f)   %s\n", fName[q], Vfull[gFix[q]].F, Vfull[gFix[q]].S, Vfull[gFix[q]].eS, PtName(gFix[q]).c_str());
    { int b = -1; for (int g = 0; g < kNT; ++g) if (!flagged(g) && (b < 0 || Vfull[g].F > Vfull[b].F)) b = g;
      o.P("   %-34s F = %6.2f  (S = %5.0f +- %4.0f)   %s\n", "best-F on the full sample", Vfull[b].F, Vfull[b].S, Vfull[b].eS, PtName(b).c_str()); }
    o.P("\nFOLD (parity) halving:\n");
    for (int q = 0; q < nR; ++q)
        o.P("   %-8s chosen on fold 0: %-55s F0 %5.2f F1 %5.2f | on fold 1: %-55s F1 %5.2f F0 %5.2f | b = %5.2f, gain = %5.2f\n", rName[q],
            PtName(H[0].c[0][q]).c_str(), val(0, 0, H[0].c[0][q]).F, val(0, 1, H[0].c[0][q]).F, PtName(H[0].c[1][q]).c_str(),
            val(0, 1, H[0].c[1][q]).F, val(0, 0, H[0].c[1][q]).F, H[0].b[q], H[0].gain[q]);
    o.P("\n%d RANDOM HALVINGS (per half; x sqrt(2) = full-sample scale):\n", nSplit);
    o.P("   %-34s | %14s | %14s | %18s\n", "rule / point", "bias b(F)", "x sqrt(2)", "independent gain");
    double bm[nR], be[nR], gm[nR], ge[nR], bfm[nFix], bfe[nFix], gfm[nFix], gfe[nFix];
    for (int q = 0; q < nR; ++q) {
        meanErr([q](const HRes& r) { return r.b[q]; }, bm[q], be[q]);
        meanErr([q](const HRes& r) { return r.gain[q]; }, gm[q], ge[q]);
        o.P("   %-34s | %6.2f +- %4.2f | %6.2f +- %4.2f | %6.2f +- %4.2f (x1.41 = %5.2f)\n", rName[q], bm[q], be[q], bm[q] * std::sqrt(2.), be[q] * std::sqrt(2.), gm[q], ge[q], gm[q] * std::sqrt(2.));
    }
    for (int q = 0; q < nFix; ++q) {
        meanErr([q](const HRes& r) { return r.bFix[q]; }, bfm[q], bfe[q]);
        meanErr([q](const HRes& r) { return r.gainFix[q]; }, gfm[q], gfe[q]);
        o.P("   %-34s | %6.2f +- %4.2f | %6.2f +- %4.2f | %6.2f +- %4.2f (x1.41 = %5.2f)\n", fName[q], bfm[q], bfe[q], bfm[q] * std::sqrt(2.), bfe[q] * std::sqrt(2.), gfm[q], gfe[q], gfm[q] * std::sqrt(2.));
    }
    o.P("\nPOINTS CHOSEN on the %d random halves (2 per halving):\n", 2 * nSplit);
    for (int q = 0; q < nR; ++q) {
        std::map<std::string, int> cnt2;
        for (int h = 1; h < nH; ++h) for (int s = 0; s < 2; ++s) cnt2[PtName(H[h].c[s][q])]++;
        std::vector<std::pair<int, std::string>> v;
        for (auto& e : cnt2) v.push_back(std::make_pair(e.second, e.first));
        std::sort(v.rbegin(), v.rend());
        o.P("   %s (%zu distinct):", rName[q], v.size());
        for (size_t i = 0; i < std::min<size_t>(v.size(), 6); ++i) o.P("  %dx [%s]", v[i].first, v[i].second.c_str());
        o.P("\n");
    }
    o.P("\nFIRST STUDY (scan/peak_bias.txt; base WITHOUT vz > -160, grid OA lo/hi x oa_lab x M(pipi) lo, no vertex cuts), per half:\n");
    o.P("   two folds: nominal 0.07, loose 0.29, pure 1.25, best-F 1.34 ; 40 random halvings: nominal 0.47 +- 0.06, loose 0.16 +- 0.03,\n");
    o.P("   pure 0.96 +- 0.06, best-F 1.13 +- 0.05 ; gain of the kinematic nominal over 'none' dF = 0.48 +- 0.47, minus bias 0.47 = 0.01.\n");
    o.P("\nrun time %.1f s\n", sw.RealTime());
    o.Close();

    // ---- LaTeX ----
    {
        FILE* t = fopen((std::string(kTabOut) + "/opt_bias.tex").c_str(), "w");
        fprintf(t, "%% opt_bias.tex -- written by plots/exp_report/opt_bias.C (numbers: tables_exp/opt_bias.txt)\n");
        fprintf(t, "\\begin{table}[H]\n\\centering\\small\n\\setlength{\\tabcolsep}{4pt}\n");
        fprintf(t, "\\caption{Optimisation bias of data-driven choices of the $\\eta$-peak selection and the gain over pre on independent events. "
                   "%d random halvings of the event groups; the working point is chosen on one half by a fixed rule (grid of %d points on top of pre: "
                   "$\\alpha_\\mathrm{lab}$, $M_{\\pi\\pi}$ lower cut, pions in the vertex fit, $\\Delta z_\\pi$, MM$_4<2.85$ or none; MM$_4<2.85$ is never chosen by the rules because of its sculpting flag) and evaluated on the other. "
                   "$b$ = mean of $\\frac12[(F_A(c_A)-F_B(c_A))+(F_B(c_B)-F_A(c_B))]$; gain = $F$ on the other half minus $F$(pre) there; per half and scaled by $\\sqrt2$ to the full sample. "
                   "Errors: rms$/\\sqrt{n}$ (too small, the halvings share events). $^*$For a fixed point the halving test gives $b=0$ by construction: it cannot see a choice made on all data, "
                   "so the points chosen on these data (peak nominal, + vertex, pre + vertex, + $\\Delta z_\\pi<40$) still carry the bias of the rule that found them "
                   "(nominal rule: %.2f on the full-sample scale); only the frozen DIRECT and PEAK are free of it. Data, isBest $\\geq0$, $P_\\mathrm{main}$.}\n", nSplit, kNT, bm[0] * std::sqrt(2.));
        fprintf(t, "\\label{tab:opt_bias}\n\\begin{tabular}{lcccc}\n\\toprule\n");
        fprintf(t, "Rule or point & $b(F)$ per half & $b(F)\\times\\sqrt2$ & gain per half & gain $\\times\\sqrt2$ \\\\\n\\midrule\n");
        for (int q = 0; q < nR; ++q)
            fprintf(t, "%s (data-driven) & $%.2f\\pm%.2f$ & %.2f & $%.2f\\pm%.2f$ & %.2f \\\\\n", rName[q], bm[q], be[q], bm[q] * std::sqrt(2.), gm[q], ge[q], gm[q] * std::sqrt(2.));
        fprintf(t, "\\midrule\n");
        const char* texName[nFix] = {"pre", "peak nominal (chosen on data)", "peak nominal + vertex (chosen on data)",
                                     "pre + vertex (chosen on data)", "peak nominal + $\\Delta z_\\pi<40$ (chosen on data)", "DIRECT (frozen on simulation)", "PEAK (frozen on simulation)"};
        for (int q = 1; q < nFix; ++q)
            fprintf(t, "%s & 0$^*$ & -- & $%.2f$ & %.2f \\\\\n", texName[q], gfm[q], gfm[q] * std::sqrt(2.));
        fprintf(t, "\\bottomrule\n\\end{tabular}\n\\end{table}\n");
        fclose(t);
    }

    // ---- figure: (a) bias per rule (this study and the first study), (b) independent gain per rule/point ----
    {
        TCanvas* c = new TCanvas("c_bias", "", 1500, 620);
        c->Divide(2, 1, 0.001, 0.001);
        // (a)
        TPad* pa = (TPad*)c->cd(1);
        pa->SetLeftMargin(0.13); pa->SetBottomMargin(0.16); pa->SetTopMargin(0.10);
        TH1D* fa = new TH1D("fba", "", nR, 0, nR); fa->SetDirectory(nullptr);
        for (int q = 0; q < nR; ++q) fa->GetXaxis()->SetBinLabel(q + 1, rPlain[q]);
        fa->SetMinimum(-0.3); fa->SetMaximum(2.2); fa->GetYaxis()->SetTitle("optimisation bias b(F), per half");
        fa->GetXaxis()->SetLabelSize(0.06); fa->Draw("AXIS");
        HLine(0, nR, 0, kGray + 1, 2, 1);
        TGraphErrors* g1 = GE(kBlue + 1, 20, 1.4); TGraphErrors* g2 = GE(kGray + 2, 24, 1.4); TGraphErrors* g3 = GE(kRed + 1, 21, 1.2);
        const double old40[nR] = {0.47, 0.16, 0.96, 1.13}, old40e[nR] = {0.06, 0.03, 0.06, 0.05};
        for (int q = 0; q < nR; ++q) {
            AddPt(g1, q + 0.35, bm[q], be[q]); AddPt(g3, q + 0.5, H[0].b[q], 0); AddPt(g2, q + 0.65, old40[q], old40e[q]);
        }
        g2->Draw("PZ SAME"); g3->Draw("P SAME"); g1->Draw("PZ SAME");
        Txt(0.13, 0.93, "(a) bias of the choice: chosen on one half, tested on the other", 0.042);
        TLegend* l = Leg(0.16, 0.66, 0.70, 0.88, 0.040);
        l->AddEntry(g1, Form("this grid (with vertex cuts), %d random halvings", nSplit), "pe");
        l->AddEntry(g3, "this grid, odd/even halves (fold)", "p");
        l->AddEntry(g2, "earlier study (no vz cut, kinematic grid)", "pe");
        l->Draw();
        // (b)
        TPad* pb = (TPad*)c->cd(2);
        pb->SetLeftMargin(0.13); pb->SetBottomMargin(0.22); pb->SetTopMargin(0.10);
        const int nb = nR + nFix - 1;
        const char* shortName[nFix] = {"base", "#splitline{#alpha_{lab}<55#circ,}{M_{#pi#pi}>0.3}", "#splitline{#alpha_{lab}<55#circ,}{M_{#pi#pi}>0.3, #pi vtx}", "base + #pi vtx", "#splitline{#alpha_{lab}<55#circ,}{#splitline{M_{#pi#pi}>0.3,}{#Deltaz_{#pi}<40 mm}}", "#Deltaz_{#pi}<40 mm", "#splitline{#Deltaz_{#pi}<40 mm,}{MM_{4}<2.85}"};
        TH1D* fb = new TH1D("fbb", "", nb, 0, nb); fb->SetDirectory(nullptr);
        for (int q = 0; q < nR; ++q) fb->GetXaxis()->SetBinLabel(q + 1, rPlain[q]);
        for (int q = 1; q < nFix; ++q) fb->GetXaxis()->SetBinLabel(nR + q, shortName[q]);
        fb->LabelsOption("v", "X"); fb->GetXaxis()->SetLabelSize(0.045);
        // verifier change 2026-11-14: y range up to 3.2 so that the legend marker sits above the data
        fb->SetMinimum(-0.3); fb->SetMaximum(3.2); fb->GetYaxis()->SetTitle("gain in F over base (full-sample scale)");
        fb->GetYaxis()->SetTitleSize(0.050);
        fb->Draw("AXIS");
        HLine(0, nb, 0, kGray + 1, 2, 1);
        TBox* bx0 = new TBox(0, -0.3, nR, 3.2); bx0->SetFillColor(TColor::GetColor("#eef0fb")); bx0->SetLineWidth(0); bx0->Draw();
        TBox* bx = new TBox(nR + 4, -0.3, nb, 3.2); bx->SetFillColor(TColor::GetColor("#fdf0d2")); bx->SetLineWidth(0); bx->Draw();
        TGraphErrors* gg = GE(kBlue + 1, 20, 1.4); TGraphErrors* gc = GE(kRed + 1, 24, 1.4);
        for (int q = 0; q < nR; ++q) AddPt(gg, q + 0.5, gm[q] * std::sqrt(2.), ge[q] * std::sqrt(2.));
        for (int q = 1; q < nFix; ++q) {
            AddPt(gg, nR + q - 0.5, gfm[q] * std::sqrt(2.), 0);
            if (q <= 4) AddPt(gc, nR + q - 0.5, (gfm[q] - bm[0]) * std::sqrt(2.), 0);   // minus the bias of the nominal rule
        }
        gc->Draw("P SAME"); gg->Draw("PZ SAME"); fb->Draw("AXIS SAME");
        Txt(0.13, 0.93, "(b) gain over base on events not used for the choice", 0.042);
        Txt(0.15, 0.84, "blue band: chosen on one half, measured on the other", 0.034, kBlue + 2);
        Txt(0.15, 0.79, "white: fixed points, chosen with all data", 0.034, kBlack);
        Txt(0.15, 0.74, "orange: frozen on simulation (no data bias)", 0.034, kOrange + 7);
        TLegend* l2 = Leg(0.15, 0.655, 0.70, 0.725, 0.034);
        l2->AddEntry(gc, "fixed point minus bias of the plateau-centre rule", "p"); l2->Draw();
        Save(c, "opt_bias");
    }
    printf("opt_bias: done in %.1f s\n", sw.RealTime());
}
