// ============================================================================
// prep_dzpi.C  --  distance of the pion tracks from the event vertex (dz_pi), data
// ============================================================================
// HADES p+p, T = 4.53 GeV (feb22), channel pi+ pi- e+ e-. EXPERIMENTAL DATA ONLY.
// dz_pi = max(|pip_z - vz|, |pim_z - vz|) [mm]; pip_z, pim_z = z of the pion track
// at its closest approach to the beam axis (PAT X_z), vz = event vertex
// (eVertReco_z). It is the variable of the frozen simulation selections PEAK
// (dz_pi < 40 mm, MM4 < 2.85 GeV/c2) and DIRECT (dz_pi < 40 mm).
//
// What it does
//   1. row pass fractions of dz_pi < X (X = 10 ... 150 mm) of the unlike-sign
//      rows and of the like-sign rows (e+e+ plus e-e-) in six M4 regions:
//      all 0.35-0.95, bump window 0.4550-0.5225, lower side band 0.5225-0.53364,
//      peak window 0.53364-0.55416, upper side band 0.55416-0.56530 (side bands
//      as in scan_extra.C), side bands + window 0.5225-0.56530;
//      error of a fraction with the event-group correlation (ratio estimator);
//   2. eta-peak efficiency eps_S = S_pass / (S_pass + S_fail) from P_main fits
//      of the pass and the fail spectra (N - CB, frozen Gaussian 0.5439 GeV/c2 /
//      5.4 MeV + pol2, 0.5225-0.6100, CB = 2 sqrt(N++ N--) in 10 MeV bins,
//      k = 1); stat error from the two fits (taken as independent); model error
//      = half range of eps_S over the 48 P_main variants in which BOTH fits have
//      p > 0.05; eps_B = kept part of (N - S) in the peak window, G = eps_S/sqrt(eps_B);
//   3. dz_pi shape of the eta-peak signal from P_main fits in dz_pi bins;
//   4. figure: (a) normalised dz_pi shapes, peak window against the side bands,
//      like-sign rows and the fitted signal; (b) cumulative pass fractions.
//
// Selections (rows = all stored combinations, isBest >= 0; corrected momenta,
// cor_ok == 1; already in the FAT production: vertex z > -500 mm, PT3 trigbit 8192)
//   pre          = cor_ok && isBest >= 0 && vz > -160 mm
//   pre0         = cor_ok && isBest >= 0 (no vz > -160; only for the cross-check
//                  with scan_extra.C, whose start selection had no vz cut)
//   peak nominal = pre + oa_lab < 55 deg + M(pi+pi-) > 0.30 GeV/c2
//   bump nominal = pre + oa_lab < 40 deg + 0.30 < M(pi+pi-) < 0.41 GeV/c2
//   variants     : isBest == 1 (comparison) ; OA(e+e-) > 2 deg (not a fixed cut)
// Windows: data peak window 0.53364-0.55416 and data bump window 0.4550-0.5225
// (the frozen simulation windows 0.5333-0.5587 / 0.4576-0.5255 are sim-centred).
//
// Inputs : /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/scan_v2/scan_compact.root
//          (scan::kCompact, written by plots/scan_prepare.C) through ../ScanEngine.h
// Outputs: /home/przygoda/HADES/FAT/FAT/docs/latex/figs_exp/prep_dzpi.pdf (+ .png; PDF via EPS + ps2pdf -dEPSCrop)
//          /home/przygoda/HADES/FAT/FAT/docs/latex/tables_exp/prep_dzpi.txt, prep_dzpi.tex
//          /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/prep/prep_dzpi.root
//
// Run (from the repository root /home/przygoda/HADES/FAT/FAT):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/aclic_prep", kTRUE)' 'plots/exp_report/prep_dzpi.C+'
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TStyle.h"
#include "TSystem.h"

#include "../ScanEngine.h"

namespace prepdz {

using scan::PeakCfg;
using scan::PeakFit;
using scan::Row;
using scan::Sel;
using scan::Spectra;

const char* const kFigOut = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_exp";
const char* const kTabOut = "/home/przygoda/HADES/FAT/FAT/docs/latex/tables_exp";
const char* const kWork   = "/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/prep";

// ---- M4 regions [GeV/c2] ------------------------------------------------------
const double kSbLoA = 0.5225, kSbHiB = 0.56530;            // side bands of scan_extra.C
const int kNReg = 6;
const double kRegLo[kNReg] = {scan::kMlo, scan::kBumpLo, kSbLoA, scan::kPeakLo, scan::kPeakHi, kSbLoA};
const double kRegHi[kNReg] = {scan::kMhi, scan::kBumpHi, scan::kPeakLo, scan::kPeakHi, kSbHiB, kSbHiB};
const char* const kRegName[kNReg] = {"all M4 0.35-0.95", "bump window 0.4550-0.5225",
                                     "lower SB 0.5225-0.5336", "peak window 0.5336-0.5542",
                                     "upper SB 0.5542-0.5653", "SB+window 0.5225-0.5653"};
enum { kAll = 0, kBumpW = 1, kSbLo = 2, kPeakW = 3, kSbHi = 4, kSbWin = 5 };

// ---- thresholds [mm] ------------------------------------------------------------
const int kNThr = 8;
const float kThr[kNThr] = {10.f, 15.f, 20.f, 30.f, 40.f, 60.f, 100.f, 150.f};
const int kI40 = 4;

// ---- dz_pi bins for the signal shape [mm] ---------------------------------------
const int kNShape = 10;
const double kShapeEdge[kNShape + 1] = {0, 5, 10, 15, 20, 30, 40, 60, 100, 150, 1e9};

// ---- selections ---------------------------------------------------------------------
inline bool Pre(const Row& r)     { return r.cor_ok == 1 && r.isBest >= 0 && r.vz > -160.f; }
inline bool Pre0(const Row& r)    { return r.cor_ok == 1 && r.isBest >= 0; }
inline bool PeakKin(const Row& r) { return r.oa_lab < 55.f && r.mpipi > 0.30f; }
inline bool BumpKin(const Row& r) { return r.oa_lab < 40.f && r.mpipi > 0.30f && r.mpipi < 0.41f; }

struct Base { std::string name, tag; Sel sel; };
std::vector<Base> Bases() {
    std::vector<Base> b;
    b.push_back({"pre", "pre", [](const Row& r) { return Pre(r); }});
    b.push_back({"pre + OA>2", "preOA2", [](const Row& r) { return Pre(r) && r.oa > 2.f; }});
    b.push_back({"pre, isBest==1", "preBest", [](const Row& r) { return Pre(r) && r.isBest == 1; }});
    b.push_back({"peak nominal", "peakNom", [](const Row& r) { return Pre(r) && PeakKin(r); }});
    b.push_back({"bump nominal", "bumpNom", [](const Row& r) { return Pre(r) && BumpKin(r); }});
    b.push_back({"pre0 + peak kin.", "pre0PeakKin", [](const Row& r) { return Pre0(r) && PeakKin(r); }});
    return b;
}

// ---- pass fraction with the event-group correlation (ratio estimator) ----------------
struct Acc {
    double N = 0, P = 0, snn = 0, spp = 0, snp = 0;
    void Add(double n, double p) { N += n; P += p; snn += n * n; spp += p * p; snp += n * p; }
    double Eps() const { return N > 0 ? P / N : 0.; }
    double Err() const {
        if (!(N > 0)) return 0.;
        const double e = Eps();
        return std::sqrt(std::max(0., spp - 2 * e * snp + e * e * snn)) / N;
    }
};

struct Eff { double e = 0, err = 0; };
// a / (a + b) for independent a, b
Eff Ratio(double a, double ea, double b, double eb) {
    Eff r; const double t = a + b;
    if (!(t > 0)) return r;
    r.e = a / t; r.err = std::sqrt(b * b * ea * ea + a * a * eb * eb) / (t * t);
    return r;
}

FILE* gOut = nullptr;
FILE* gTex = nullptr;
void P(const char* fmt, ...) {
    va_list ap; va_start(ap, fmt); vprintf(fmt, ap); va_end(ap);
    if (gOut) { va_start(ap, fmt); vfprintf(gOut, fmt, ap); va_end(ap); }
}

void Style() {
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    gStyle->SetTextFont(42); gStyle->SetLabelFont(42, "xyz"); gStyle->SetTitleFont(42, "xyz");
    gStyle->SetLabelSize(0.050, "xyz"); gStyle->SetTitleSize(0.056, "xyz");
    gStyle->SetTitleOffset(1.05, "x"); gStyle->SetTitleOffset(1.25, "y");
    gStyle->SetPadLeftMargin(0.15); gStyle->SetPadRightMargin(0.03);
    gStyle->SetPadTopMargin(0.08); gStyle->SetPadBottomMargin(0.14);
    gStyle->SetEndErrorSize(0);
}

struct PF { Eff eps; double spread = 0; int nvar = 0; double epsB = 0, G = 0, eG = 0;
            PeakFit pass, fail, base; };

PF PassFail(const Spectra& sb, const Spectra& sp, const Spectra& sf) {
    PF o;
    o.base = scan::FitPeakMain(sb); o.pass = scan::FitPeakMain(sp); o.fail = scan::FitPeakMain(sf);
    o.eps = Ratio(o.pass.S, o.pass.eS, o.fail.S, o.fail.eS);
    const double ba = o.pass.N - o.pass.Swin, bb = o.fail.N - o.fail.Swin;
    o.epsB = (ba + bb) > 0 ? ba / (ba + bb) : 0;
    if (o.epsB > 0) { o.G = o.eps.e / std::sqrt(o.epsB); o.eG = o.eps.err / std::sqrt(o.epsB); }
    double emin = 1e9, emax = -1e9;
    for (int ir = 0; ir < scan::kNVarRange; ++ir)
        for (int ord = scan::kVarOrdMin; ord <= scan::kVarOrdMax; ++ord)
            for (int iw = -1; iw <= 1; ++iw) {
                PeakCfg c;
                c.lo = scan::kVarLo[ir]; c.hi = scan::kVarHi[ir]; c.order = ord;
                c.sigma = scan::kPeakSigma + iw * scan::kVarDSigma;
                const PeakFit a = scan::FitPeakMain(sp, c), b = scan::FitPeakMain(sf, c);
                if (!a.ok || !b.ok || !(a.prob > scan::kVarPmin) || !(b.prob > scan::kVarPmin)) continue;
                if (!(a.S + b.S > 0)) continue;
                const double e = a.S / (a.S + b.S);
                emin = std::min(emin, e); emax = std::max(emax, e);
                o.nvar++;
            }
    o.spread = o.nvar > 0 ? 0.5 * (emax - emin) : 0.;
    return o;
}

} // namespace prepdz

void prep_dzpi() {
    using namespace prepdz;
    TStopwatch swAll;
    Style();
    gSystem->mkdir(kFigOut, kTRUE); gSystem->mkdir(kTabOut, kTRUE); gSystem->mkdir(kWork, kTRUE);
    gOut = fopen((std::string(kTabOut) + "/prep_dzpi.txt").c_str(), "w");

    scan::Data d;
    if (!scan::Load(d)) return;
    {
        long nnan = 0;
        for (int k = 0; k < 3; ++k) for (const auto& r : d.s[k]) if (std::isnan(r.dz_pi)) ++nnan;
        if (nnan > 0) { printf("prep_dzpi: %ld rows without dz_pi: wrong compact file\n", nnan); return; }
    }
    const std::vector<Base> bases = Bases();
    const int nb = (int)bases.size();

    // ------------------------------------------------------------------------
    // 1. row pass fractions + histograms (one pass, event-group accumulation)
    // ------------------------------------------------------------------------
    // sample index: 0 us, 1 pp, 2 mm, 3 like-sign pp + mm
    std::vector<Acc> acc(nb * 4 * kNReg * kNThr);
    auto A = [&](int b, int k, int q, int t) -> Acc& { return acc[((b * 4 + k) * kNReg + q) * kNThr + t]; };
    std::vector<TH1D*> hUS(nb * kNReg), hLS(nb * kNReg);
    for (int b = 0; b < nb; ++b)
        for (int q = 0; q < kNReg; ++q) {
            hUS[b * kNReg + q] = new TH1D(Form("h_us_%s_%d", bases[b].tag.c_str(), q), "", 100, 0, 200);
            hLS[b * kNReg + q] = new TH1D(Form("h_ls_%s_%d", bases[b].tag.c_str(), q), "", 100, 0, 200);
            hUS[b * kNReg + q]->SetDirectory(nullptr); hLS[b * kNReg + q]->SetDirectory(nullptr);
            hUS[b * kNReg + q]->Sumw2(); hLS[b * kNReg + q]->Sumw2();
        }
    for (int k = 0; k < 3; ++k) {
        const std::vector<Row>& v = d.s[k];
        std::vector<int> cn(nb * kNReg, 0), cp(nb * kNReg * kNThr, 0);
        int cur = -1; bool dirty = false;
        auto flush = [&]() {
            for (int b = 0; b < nb; ++b)
                for (int q = 0; q < kNReg; ++q) {
                    int& n = cn[b * kNReg + q];
                    if (!n) continue;
                    for (int t = 0; t < kNThr; ++t) {
                        int& p = cp[(b * kNReg + q) * kNThr + t];
                        A(b, k, q, t).Add(n, p);
                        if (k > 0) A(b, 3, q, t).Add(n, p);
                        p = 0;
                    }
                    n = 0;
                }
            dirty = false;
        };
        for (size_t i = 0; i < v.size(); ++i) {
            const Row& r = v[i];
            if (r.grp != cur) { if (dirty) flush(); cur = r.grp; }
            if (r.cor_ok != 1) continue;
            const double m = r.m4;
            if (!(m == m)) continue;
            for (int b = 0; b < nb; ++b) {
                if (!bases[b].sel(r)) continue;
                for (int q = 0; q < kNReg; ++q) {
                    if (!(m >= kRegLo[q] && m < kRegHi[q])) continue;
                    cn[b * kNReg + q]++;
                    for (int t = 0; t < kNThr; ++t)
                        if (r.dz_pi < kThr[t]) cp[(b * kNReg + q) * kNThr + t]++;
                    dirty = true;
                    (k == 0 ? hUS : hLS)[b * kNReg + q]->Fill(r.dz_pi);
                }
            }
        }
        if (dirty) flush();
    }

    // ------------------------------------------------------------------------
    // 2./3. spectra: bases, pass/fail per threshold, signal shape bins (pre)
    // ------------------------------------------------------------------------
    std::vector<Sel> sels;
    for (int b = 0; b < nb; ++b) sels.push_back(bases[b].sel);
    for (int b = 0; b < nb; ++b)
        for (int t = 0; t < kNThr; ++t) {
            const Sel s = bases[b].sel; const float c = kThr[t];
            sels.push_back([s, c](const Row& r) { return s(r) && r.dz_pi < c; });
            sels.push_back([s, c](const Row& r) { return s(r) && !(r.dz_pi < c); });
        }
    const int iShape0 = (int)sels.size();
    for (int j = 0; j < kNShape; ++j) {
        const float a = (float)kShapeEdge[j], c = (float)kShapeEdge[j + 1];
        sels.push_back([a, c](const Row& r) { return Pre(r) && r.dz_pi >= a && r.dz_pi < c; });
    }
    scan::FillOpt fo;
    fo.windows.push_back(std::make_pair(kSbLoA, scan::kPeakLo));
    fo.windows.push_back(std::make_pair(scan::kPeakHi, kSbHiB));
    TStopwatch sw;
    const std::vector<Spectra> sp = scan::FillMany(d, sels, fo);
    printf("FillMany: %zu selections in %.1f s\n", sels.size(), sw.RealTime());
    auto iPass = [&](int b, int t) { return nb + 2 * (b * kNThr + t); };

    std::vector<PF> pf(nb * kNThr);
    for (int b = 0; b < nb; ++b)
        for (int t = 0; t < kNThr; ++t)
            pf[b * kNThr + t] = PassFail(sp[b], sp[iPass(b, t)], sp[iPass(b, t) + 1]);

    // ------------------------------------------------------------------------
    // tables
    // ------------------------------------------------------------------------
    P("# prep_dzpi.txt  --  written by plots/exp_report/prep_dzpi.C\n");
    P("# data: %s\n", scan::kCompact);
    P("# EXPERIMENTAL DATA. rows = stored particle combinations (isBest >= 0 unless stated), corrected momenta, cor_ok == 1;\n");
    P("#   production: vertex z > -500 mm, PT3 (trigbit 8192). pre = + vz > -160 mm. pre0 = without vz > -160 (cross-check only).\n");
    P("#   peak nominal = pre + oa_lab < 55 + M(pipi) > 0.30 ; bump nominal = pre + oa_lab < 40 + 0.30 < M(pipi) < 0.41\n");
    P("# dz_pi = max(|pip_z - vz|, |pim_z - vz|) [mm]; pip_z, pim_z = z of the pion track at the closest approach to the beam axis\n");
    P("#   (same definition as the simulation study, FAT_simwt plots/eta_sim/s5_blind_prep.C DzPi). Frozen sim cut: dz_pi < 40 mm.\n");
    P("# Windows: data peak window 0.53364-0.55416, side bands 0.5225-0.53364 and 0.55416-0.56530 (as scan_extra.C), data bump window 0.4550-0.5225.\n\n");

    // ---- A. row pass fractions at 40 mm ----
    P("A. ROW PASS FRACTION of dz_pi < 40 mm (rows passing / rows; error with the event-group correlation; samples treated separately)\n");
    P("%-17s %-27s | %9s %15s | %8s %15s | %8s %15s | %8s %15s\n", "selection", "M4 region", "N e+e-", "f(e+e-)", "N e+e+",
      "f(e+e+)", "N e-e-", "f(e-e-)", "N LS", "f(LS = ++ + --)");
    for (int b = 0; b < nb; ++b) {
        for (int q = 0; q < kNReg; ++q) {
            P("%-17s %-27s |", bases[b].name.c_str(), kRegName[q]);
            for (int k = 0; k < 4; ++k) {
                const Acc& a = A(b, k, q, kI40);
                P(" %8.0f  %6.4f+-%6.4f |", a.N, a.Eps(), a.Err());
            }
            P("\n");
        }
        P("\n");
    }

    // ---- A2. thresholds, pre ----
    P("A2. ROW PASS FRACTION vs threshold X (dz_pi < X), selection pre (unlike-sign e+e- unless stated)\n");
    P("%-32s |", "region / sample");
    for (int t = 0; t < kNThr; ++t) P("  X=%-5.0f       ", kThr[t]);
    P("\n");
    const int qList[5] = {kAll, kBumpW, kSbLo, kPeakW, kSbHi};
    for (int qi = 0; qi < 5; ++qi) {
        P("%-32s |", kRegName[qList[qi]]);
        for (int t = 0; t < kNThr; ++t) P(" %6.4f+-%6.4f", A(0, 0, qList[qi], t).Eps(), A(0, 0, qList[qi], t).Err());
        P("\n");
    }
    P("%-32s |", "LS rows, SB+window 0.5225-0.5653");
    for (int t = 0; t < kNThr; ++t) P(" %6.4f+-%6.4f", A(0, 3, kSbWin, t).Eps(), A(0, 3, kSbWin, t).Err());
    P("\n%-32s |", "LS rows, all M4");
    for (int t = 0; t < kNThr; ++t) P(" %6.4f+-%6.4f", A(0, 3, kAll, t).Eps(), A(0, 3, kAll, t).Err());
    P("\n\n");

    // ---- A3. CB and N - CB fractions from the spectra (40 mm) ----
    P("A3. CB = 2 sqrt(N++ N--) (10 MeV bins, interpolated, k = 1) and N - CB kept by dz_pi < 40 mm, in the windows (no errors)\n");
    P("%-17s | %8s %8s %8s %8s | %10s %10s %8s\n", "selection", "CB peak", "CB lowSB", "CB upSB", "CB bump", "N-CB bump", "kept", "frac");
    for (int b = 0; b < nb; ++b) {
        const Spectra& s0 = sp[b]; const Spectra& s1 = sp[iPass(b, kI40)];
        auto fr = [](double a, double c) { return c > 0 ? a / c : 0.; };
        P("%-17s | %8.4f %8.4f %8.4f %8.4f | %10.0f %10.0f %8.4f\n", bases[b].name.c_str(), fr(s1.win[0].cb, s0.win[0].cb),
          fr(s1.win[2].cb, s0.win[2].cb), fr(s1.win[3].cb, s0.win[3].cb), fr(s1.win[1].cb, s0.win[1].cb), s0.win[1].sig(),
          s1.win[1].sig(), fr(s1.win[1].sig(), s0.win[1].sig()));
    }
    P("\n");

    // ---- B. signal efficiency ----
    P("B. ETA-PEAK SIGNAL EFFICIENCY eps_S = S_pass / (S_pass + S_fail), P_main fits (frozen 0.5439 / 5.4 MeV, pol2, 0.5225-0.6100, N - CB)\n");
    P("#  stat: from the two fits, taken as independent; model: half range over the variants (of 48) with p > 0.05 in BOTH fits (nvar)\n");
    P("#  eps_B = kept part of (N - S_window) in the peak window; G = eps_S / sqrt(eps_B); S, F = S/stat of the base and after the cut;\n");
    P("#  clos = S_pass + S_fail - S_base\n");
    P("%-17s %5s | %6s %6s %6s %4s | %6s %6s %5s | %7s %5s %6s | %7s %5s %6s | %7s %5s | %5s | %9s\n", "selection", "X", "eps_S", "stat",
      "model", "nvar", "eps_B", "G", "err", "S_base", "err", "F_base", "S_pass", "err", "F_pass", "S_fail", "err", "clos",
      "chi2/ndf p");
    for (int b = 0; b < nb; ++b) {
        for (int t = 0; t < kNThr; ++t) {
            const PF& o = pf[b * kNThr + t];
            P("%-17s %5.0f | %6.3f %6.3f %6.3f %4d | %6.3f %6.3f %5.3f | %7.0f %5.0f %6.2f | %7.0f %5.0f %6.2f | %7.0f %5.0f | %5.0f | %4.1f/%-2d %4.1f/%-2d\n",
              bases[b].name.c_str(), kThr[t], o.eps.e, o.eps.err, o.spread, o.nvar, o.epsB, o.G, o.eG, o.base.S, o.base.eS,
              o.base.F, o.pass.S, o.pass.eS, o.pass.F, o.fail.S, o.fail.eS, o.pass.S + o.fail.S - o.base.S, o.pass.chi2,
              o.pass.ndf, o.fail.chi2, o.fail.ndf);
        }
        P("\n");
    }

    // ---- C. signal shape ----
    P("C. dz_pi SHAPE OF THE ETA-PEAK SIGNAL (selection pre): P_main fit in each dz_pi bin\n");
    double sumS = 0, sumV = 0;
    std::vector<PeakFit> shp(kNShape);
    for (int j = 0; j < kNShape; ++j) {
        shp[j] = scan::FitPeakMain(sp[iShape0 + j]);
        sumS += shp[j].S; sumV += shp[j].eS * shp[j].eS;
    }
    P("%-14s | %7s %6s | %7s %6s | %8s | %9s\n", "dz_pi [mm]", "S", "err", "frac", "err", "N window", "chi2/ndf");
    for (int j = 0; j < kNShape; ++j)
        P("%5.0f - %-6.0f | %7.0f %6.0f | %7.4f %6.4f | %8.0f | %5.1f/%d\n", kShapeEdge[j], std::min(kShapeEdge[j + 1], 9999.),
          shp[j].S, shp[j].eS, shp[j].S / sumS, shp[j].eS / sumS, shp[j].N, shp[j].chi2, shp[j].ndf);
    P("sum of bins S = %.0f +- %.0f (base fit S = %.0f +- %.0f)\n\n", sumS, std::sqrt(sumV), pf[0].base.S, pf[0].base.eS);

    // ---- D. simulation values (read from the simulation tables, not recomputed) ----
    P("D. SIMULATION VALUES FOR COMPARISON (copied from FAT_simwt/plots/eta_sim/tables; NOT recomputed here)\n");
    P("   s5_opt_eta_frozen.txt [DIRECT]: efficiency of dz_pi < 40 for class 5 (eta -> pi+pi-e+e-) in the sim window 0.5333-0.5587,\n");
    P("     base pre, physics weights V3: 0.9625 +- 0.0005 ; [PEAK] (dz_pi < 40 and MM4 < 2.85, classes 4+5): 0.9310 +- 0.0072\n");
    P("   eta_peak_eff.txt (all M4, mode A, PDG weights), dz_pi < X at base pre: direct / conversion / total\n");
    P("     X=20 0.809/0.847/0.827 ; 30 0.917/0.932/0.924 ; 40 0.955/0.954+-0.012/0.955+-0.006 ; 60 0.982/0.977/0.980 ;\n");
    P("     100 0.995/0.987/0.991 ; 150 0.998/0.993/0.996\n");
    P("   eta_peak_eff.txt, base 'start' (pre0 + oa_lab<55 + M(pipi)>0.30), X=40: 0.957/0.952+-0.013/0.955+-0.006 vs data(old) 0.968+-0.036\n");
    P("   sim_eta_peak.C / eta_peak_eff.txt, base 'md' (pre0 + oa_lab<50 + M(pipi)<0.46): sim total 0.956+-0.006 vs data 0.974+-0.038\n");
    P("     (the 0.974 +- 0.038 is a DATA value of the first scan_extra.C version, hard-coded in sim_eta_peak.C)\n");
    P("   cuts_eff_single.txt (base pre, std + lep, w PDG), dz_pi < 40: bump Dalitz 0.9448+-0.0056, bump conversion 0.9537+-0.0061,\n");
    P("     peak conversion 0.9544+-0.0119, S_bump (bump in 0.4550-0.5225) 0.9455+-0.0049, physical background in the peak window\n");
    P("     B_peak 0.8265+-0.0069, in the bump window B_bump 0.8075+-0.0061 (simulation has NO combinatorial background)\n\n");

    // ---- summary lines ----
    const PF& m0 = pf[0 * kNThr + kI40];
    P("SUMMARY (dz_pi < 40 mm)\n");
    P("  pre: rows e+e- all M4 %.4f, peak window %.4f +- %.4f, side bands %.4f / %.4f, like-sign (SB+window) %.4f +- %.4f\n",
      A(0, 0, kAll, kI40).Eps(), A(0, 0, kPeakW, kI40).Eps(), A(0, 0, kPeakW, kI40).Err(), A(0, 0, kSbLo, kI40).Eps(),
      A(0, 0, kSbHi, kI40).Eps(), A(0, 3, kSbWin, kI40).Eps(), A(0, 3, kSbWin, kI40).Err());
    for (int b = 0; b < nb; ++b) {
        const PF& o = pf[b * kNThr + kI40];
        P("  %-17s eps_S = %.3f +- %.3f (stat) +- %.3f (model, %d var.), eps_B %.3f, G %.3f +- %.3f, F %.2f -> %.2f\n",
          bases[b].name.c_str(), o.eps.e, o.eps.err, o.spread, o.nvar, o.epsB, o.G, o.eG, o.base.F, o.pass.F);
    }
    P("  pull (sim DIRECT 0.9625 - data pre) / err_data = %+.2f ; (sim total 0.955 - data pre) / sqrt(...) = %+.2f\n",
      (0.9625 - m0.eps.e) / std::hypot(m0.eps.err, m0.spread),
      (0.955 - m0.eps.e) / std::sqrt(m0.eps.err * m0.eps.err + m0.spread * m0.spread + 0.006 * 0.006));
    fclose(gOut); gOut = nullptr;

    // ------------------------------------------------------------------------
    // LaTeX snippet
    // ------------------------------------------------------------------------
    gTex = fopen((std::string(kTabOut) + "/prep_dzpi.tex").c_str(), "w");
    if (gTex) {
        FILE* o = gTex;
        fprintf(o, "%% prep_dzpi.tex -- written by plots/exp_report/prep_dzpi.C (numbers: tables_exp/prep_dzpi.txt)\n");
        fprintf(o, "\\begin{table}[H]\n\\centering\\small\n\\setlength{\\tabcolsep}{4pt}\n");
        fprintf(o, "\\caption{Pion-vertex distance $\\Delta z_\\pi=\\max(|z_{\\pi^+}-v_z|,|z_{\\pi^-}-v_z|)$ in the data: pass fraction of "
                   "$\\Delta z_\\pi<40$~mm. Rows: all stored combinations (isBest $\\geq 0$) unless stated, corrected momenta, "
                   "PT3, $v_z>-160$~mm. Row fractions: data peak window 0.53364--0.55416, side bands 0.5225--0.53364 and "
                   "0.55416--0.56530~GeV/$c^2$, errors with the event-group correlation. Signal: $\\varepsilon_S=S_\\mathrm{pass}/"
                   "(S_\\mathrm{pass}+S_\\mathrm{fail})$ from $P_\\mathrm{main}$ fits (stat $\\pm$ model). Not normalised to luminosity. "
                   "Simulation: values copied from the simulation study (no combinatorial background there).}\n");
        fprintf(o, "\\label{tab:prep_dzpi}\n\\begin{tabular}{llc}\n\\toprule\n");
        fprintf(o, "Sample & Selection & pass fraction $\\Delta z_\\pi<40$~mm \\\\\n\\midrule\n");
        auto rowf = [&](const char* s, const char* sel, const Acc& a) {
            fprintf(o, "%s & %s & $%.3f\\pm%.3f$ \\\\\n", s, sel, a.Eps(), a.Err()); };
        rowf("$e^+e^-$ rows, all $M_4$", "pre", A(0, 0, kAll, kI40));
        rowf("$e^+e^-$ rows, peak window", "pre", A(0, 0, kPeakW, kI40));
        rowf("$e^+e^-$ rows, lower side band", "pre", A(0, 0, kSbLo, kI40));
        rowf("$e^+e^-$ rows, upper side band", "pre", A(0, 0, kSbHi, kI40));
        rowf("$e^\\pm e^\\pm$ rows, 0.5225--0.5653", "pre", A(0, 3, kSbWin, kI40));
        rowf("$e^+e^-$ rows, bump window", "pre", A(0, 0, kBumpW, kI40));
        fprintf(o, "\\midrule\n");
        const char* texName[6] = {"pre", "pre + OA$>2^\\circ$", "pre, isBest $=1$", "peak nominal", "bump nominal",
                                  "no $v_z$ cut + $\\alpha_\\mathrm{lab}<55^\\circ$, $M_{\\pi\\pi}>0.30$"};
        for (int b = 0; b < nb; ++b) {
            const PF& q = pf[b * kNThr + kI40];
            fprintf(o, "$\\eta$ peak signal $\\varepsilon_S$ & %s & $%.3f\\pm%.3f\\pm%.3f$ \\\\\n", texName[b], q.eps.e, q.eps.err,
                    q.spread);
        }
        fprintf(o, "\\midrule\n");
        fprintf(o, "Simulation, $\\eta\\to\\pi^+\\pi^-e^+e^-$ (sim window) & pre (DIRECT) & $0.9625\\pm0.0005$ \\\\\n");
        fprintf(o, "Simulation, direct / conversion / total (all $M_4$) & pre & $0.955$ / $0.954$ / $0.955\\pm0.006$ \\\\\n");
        fprintf(o, "Simulation, physical background in the peak window & pre & $0.827\\pm0.007$ \\\\\n");
        fprintf(o, "\\bottomrule\n\\end{tabular}\n\\end{table}\n");
        fclose(gTex);
    }

    // ------------------------------------------------------------------------
    // 4. figure
    // ------------------------------------------------------------------------
    TFile fr((std::string(kWork) + "/prep_dzpi.root").c_str(), "RECREATE");
    for (auto h : hUS) h->Write();
    for (auto h : hLS) h->Write();

    auto norm = [](TH1D* h) {
        TH1D* c = (TH1D*)h->Clone(Form("%s_n", h->GetName())); c->SetDirectory(nullptr);
        const double tot = h->Integral(0, h->GetNbinsX() + 1);
        if (tot > 0) c->Scale(1.0 / tot, "width");
        return c;
    };
    TH1D* hW  = norm(hUS[0 * kNReg + kPeakW]);
    TH1D* hL  = norm(hUS[0 * kNReg + kSbLo]);
    TH1D* hU  = norm(hUS[0 * kNReg + kSbHi]);
    TH1D* hLs = norm(hLS[0 * kNReg + kSbWin]);
    TGraphErrors* gS = new TGraphErrors();
    for (int j = 0; j < kNShape; ++j) {
        if (kShapeEdge[j + 1] > 1e3) continue;
        if (!(shp[j].S > 0)) continue;           // log scale: bins with S <= 0 are not drawn (see table C)
        const double w = kShapeEdge[j + 1] - kShapeEdge[j];
        const double y = shp[j].S / sumS / w, ey = shp[j].eS / sumS / w;
        const int n = gS->GetN();
        gS->SetPoint(n, 0.5 * (kShapeEdge[j] + kShapeEdge[j + 1]), y);
        gS->SetPointError(n, 0.5 * w, ey);
    }

    TCanvas* c = new TCanvas("c_dzpi", "", 1600, 720);
    c->Divide(2, 1, 0.004, 0.004);
    // (a) shapes
    c->cd(1);
    gPad->SetLogy();
    TH1D* fa = new TH1D("fa", "", 100, 0, 200);
    fa->SetDirectory(nullptr);
    fa->GetXaxis()->SetNdivisions(505);
    fa->GetXaxis()->SetTitle("#Deltaz_{#pi} = max(|z_{#pi^{+}} - v_{z}|, |z_{#pi^{-}} - v_{z}|) [mm]");
    fa->GetYaxis()->SetTitle("1/N dN/d#Deltaz_{#pi} [1/mm]");
    fa->SetMinimum(2e-5); fa->SetMaximum(0.6);
    fa->GetXaxis()->SetTitleSize(0.050);
    fa->Draw("axis");
    hL->SetLineColor(kOrange + 7); hL->SetLineWidth(3); hL->Draw("hist same");
    hU->SetLineColor(kGreen + 2); hU->SetLineWidth(3); hU->SetLineStyle(2); hU->Draw("hist same");
    hLs->SetMarkerStyle(25); hLs->SetMarkerColor(kRed + 1); hLs->SetLineColor(kRed + 1); hLs->SetMarkerSize(1.0);
    hLs->Draw("e1 same");
    hW->SetMarkerStyle(20); hW->SetMarkerColor(kBlack); hW->SetLineColor(kBlack); hW->SetMarkerSize(0.9);
    hW->Draw("e1 same");
    gS->SetMarkerStyle(21); gS->SetMarkerColor(kBlue + 1); gS->SetLineColor(kBlue + 1); gS->SetLineWidth(2);
    gS->SetMarkerSize(1.3);
    gS->Draw("pz same");
    TLine* l40 = new TLine(40, 2e-5, 40, 0.6); l40->SetLineStyle(7); l40->SetLineColor(kGray + 2); l40->SetLineWidth(2);
    l40->Draw();
    TLegend* la = new TLegend(0.44, 0.60, 0.97, 0.91);
    la->SetBorderSize(0); la->SetFillStyle(0); la->SetTextSize(0.040);
    la->AddEntry(hW, "e^{+}e^{-} rows, peak window", "pe");
    la->AddEntry(hL, "e^{+}e^{-} rows, lower side band", "l");
    la->AddEntry(hU, "e^{+}e^{-} rows, upper side band", "l");
    la->AddEntry(hLs, "e^{#pm}e^{#pm} rows, 0.5225-0.5653", "pe");
    la->AddEntry(gS, "#eta-peak signal (fit per bin)", "pe");
    la->Draw();
    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.042);
    tx.DrawLatex(0.16, 0.94, "(a) data, selection pre: normalised shapes");

    // (b) cumulative pass fractions
    c->cd(2);
    TH1D* fb = new TH1D("fb", "", 100, 0, 160);
    fb->SetDirectory(nullptr);
    fb->GetXaxis()->SetNdivisions(505);
    fb->GetXaxis()->SetTitle("threshold X of #Deltaz_{#pi} < X [mm]");
    fb->GetYaxis()->SetTitle("pass fraction");
    fb->SetMinimum(0.0); fb->SetMaximum(1.08);
    fb->Draw("axis");
    auto mk = [&](int k, int q, int col, int ms, int ls) {
        TGraphErrors* g = new TGraphErrors();
        for (int t = 0; t < kNThr; ++t) { g->SetPoint(t, kThr[t], A(0, k, q, t).Eps()); g->SetPointError(t, 0, A(0, k, q, t).Err()); }
        g->SetMarkerStyle(ms); g->SetMarkerColor(col); g->SetLineColor(col); g->SetLineStyle(ls); g->SetLineWidth(2);
        g->SetMarkerSize(1.1);
        return g;
    };
    TGraphErrors* cW = mk(0, kPeakW, kBlack, 20, 1);
    TGraphErrors* cL = mk(0, kSbLo, kOrange + 7, 24, 1);
    TGraphErrors* cU = mk(0, kSbHi, kGreen + 2, 26, 2);
    TGraphErrors* cLs = mk(3, kSbWin, kRed + 1, 25, 1);
    TGraphErrors* cS = new TGraphErrors();
    for (int t = 0; t < kNThr; ++t) {
        const PF& o = pf[0 * kNThr + t];
        cS->SetPoint(t, kThr[t] + 1.5, o.eps.e);
        cS->SetPointError(t, 0, std::hypot(o.eps.err, o.spread));
    }
    cS->SetMarkerStyle(21); cS->SetMarkerColor(kBlue + 1); cS->SetLineColor(kBlue + 1); cS->SetLineWidth(2); cS->SetMarkerSize(1.3);
    cL->Draw("lp same"); cU->Draw("lp same"); cLs->Draw("lp same"); cW->Draw("lp same"); cS->Draw("pz same");
    TLine* m40 = new TLine(40, 0.0, 40, 1.08); m40->SetLineStyle(7); m40->SetLineColor(kGray + 2); m40->SetLineWidth(2);
    m40->Draw();
    TLegend* lb = new TLegend(0.40, 0.17, 0.97, 0.47);
    lb->SetBorderSize(0); lb->SetFillStyle(0); lb->SetTextSize(0.040);
    lb->AddEntry(cW, "e^{+}e^{-} rows, peak window", "lp");
    lb->AddEntry(cL, "e^{+}e^{-} rows, lower side band", "lp");
    lb->AddEntry(cU, "e^{+}e^{-} rows, upper side band", "lp");
    lb->AddEntry(cLs, "e^{#pm}e^{#pm} rows, 0.5225-0.5653", "lp");
    lb->AddEntry(cS, "#eta-peak signal #varepsilon_{S} (fits)", "pe");
    lb->Draw();
    tx.DrawLatex(0.16, 0.94, "(b) data, selection pre: pass fractions");

    // vector PDF cropped to the canvas: EPS -> ps2pdf -dEPSCrop (as in the simulation report)
    const std::string fbase = std::string(kFigOut) + "/prep_dzpi";
    c->SaveAs((fbase + ".eps").c_str());
    gSystem->Exec(Form("ps2pdf -dEPSCrop %s.eps %s.pdf && rm -f %s.eps", fbase.c_str(), fbase.c_str(), fbase.c_str()));
    c->SaveAs((fbase + ".png").c_str());
    fr.cd();
    c->Write();
    fr.Close();
    printf("prep_dzpi: done in %.1f s\n", swAll.RealTime());
}
