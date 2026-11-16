// res_mee_cnt.C -- COPY of plots/exp_report/res_mee.C made by figure agent A (integ4, 2026-11-15) with
// the figure strings 'rows' -> 'counts' (analyst's terminology) and the outputs of the copy redirected:
// figures of the scope -> docs/latex/figs_exp/<base>_cnt.{pdf,png}; everything else -> /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figA/work/res_mee
// Built by /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figA/work/mk_cnt.py; the physics code is unchanged.
// ============================================================================
// res_mee.C  --  content of the eta peak versus M(e+e-) and OA(e+e-); the X17 window (data only)
// ============================================================================
// HADES p+p, T = 4.53 GeV (feb22, 29 days), channel pi+ pi- e+ e-. EXPERIMENTAL
// DATA ONLY. Chapter "Results without simulation" (prefix res), tasks 3 and 4.
//
// Method (side-band subtraction in M4, res_common.h section D), per class c
// (a bin of M(ee) or OA, optionally split by the RICH-ring flag):
//   D_c = (N - CB)_W,c - r (N - CB)_SB,c          W  = data peak window 0.53364-0.55416
//                                                 SB = 0.5225-0.53364 and 0.55416-0.56530
//   r   = integral of the P_main pol2 of the base selection (all classes) over W / over SB
//   S_c = D_c / (fW - r fSB)  (whole-Gaussian content; fW, fSB = shares of the frozen
//         Gaussian 0.5439 / 5.4 MeV in W and in SB)
//   CB  = 2 sqrt(N++ N--) of the same class (10 MeV bins, interpolated), k = 1.
//   Errors per event group (unlike-sign weights +1 / -r) and per like-sign group.
//   Model error: r from the P_main fits with pol1 and pol3 instead of pol2 (half range).
//   Check: P_main fit of each class separately (method of plots/scan_peak.C, peak_classes).
// Task 3: classes M(ee) 0-5, 5-10, 10-20, 20-50, 50-100, 100-270 MeV/c2 and
//   OA 0-1, 1-2, 2-4, 4-6, 6-10, 10-20, > 20 deg, each split into "shared RICH ring"
//   (shared_ring == 1) and "two rings"; shares and the share above 20 MeV/c2.
// Task 4: X17 window. Base = pre + both pions in the vertex fit (the frozen X17
//   selection without its M(ee) window; data peak window instead of the simulation
//   window 0.5333-0.5587). M(ee) windows of 9 MeV: 6.5-15.5, 15.5-24.5 (X17),
//   24.5-33.5, 33.5-42.5, 42.5-51.5, 51.5-60.5: N, CB, N - CB in W, side-band
//   background, D and S. No limit (needs the simulated efficiency; next report).
// Selections (rows = stored combinations; cor_ok == 1; FAT: vertex z > -500 mm, PT3):
//   pre = cor_ok && isBest >= 0 && vz > -160 mm; pre isBest == 1; peak nominal = pre +
//   oa_lab < 55 + M(pipi) > 0.30; peak nominal + vertex = + both pions in the vertex fit;
//   "nominal, no vz" = oa_lab < 55 + M(pipi) > 0.30 WITHOUT vz > -160 (definition of the
//   earlier number 30.6 +- 3.8 % above 20 MeV/c2, scan/peak_classes.txt).
//
// Inputs : /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/scan_v2/scan_compact.root (scan::kCompact)
// Outputs: docs/latex/tables_exp/res_mee.txt, res_mee.tex, res_x17.tex
//          docs/latex/figs_exp/res_mee_content.pdf/.png, res_x17_window.pdf/.png
//          /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/res/res_mee.root
//
// Run (from /home/przygoda/HADES/FAT/FAT, about 1 min):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/aclic_res", kTRUE)' 'plots/exp_report/res_mee.C+'
// ============================================================================
#include "TBox.h"
#include "TDatime.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TPad.h"
#include "TStopwatch.h"

#include "res_common_cnt.h"

namespace resme {

using namespace res;

const int kNMee = 6;
const double kMeeE[kNMee + 1] = {0.0, 0.005, 0.010, 0.020, 0.050, 0.100, 0.270};
const char* const kMeeName[kNMee] = {"0-5", "5-10", "10-20", "20-50", "50-100", "100-270"};
const int kNOa = 7;
const double kOaE[kNOa + 1] = {0, 1, 2, 4, 6, 10, 20, 180};
const char* const kOaName[kNOa] = {"0-1", "1-2", "2-4", "4-6", "6-10", "10-20", ">20"};
// X17 windows (9 MeV)
const int kNX = 6;
const double kXlo[kNX] = {0.0065, 0.0155, 0.0245, 0.0335, 0.0425, 0.0515};
const double kXhi[kNX] = {0.0155, 0.0245, 0.0335, 0.0425, 0.0515, 0.0605};
// fine M(ee) bins for the X17 figure: 3 MeV from 0.5 to 60.5 MeV
const int kNF = 20;
inline double FineLo(int k) { return 0.0005 + 0.003 * k; }

struct Base { std::string key, tex; Sel f; };

// classes of one base: 0 = all; 1..6 M(ee); 7..13 OA; then ring / two rings for the 13 bins
const int kNCls = 1 + kNMee + kNOa + 2 * (kNMee + kNOa);
inline Sel ClassSel(int c) {
    if (c == 0) return [](const Row&) { return true; };
    int b = c - 1;
    int ring = -1;
    if (b >= kNMee + kNOa) { const int q = b - (kNMee + kNOa); ring = q % 2; b = q / 2; }
    Sel s;
    if (b < kNMee) { const float a = kMeeE[b], e = kMeeE[b + 1]; s = [a, e](const Row& r) { return r.mee >= a && r.mee < e; }; }
    else { const float a = kOaE[b - kNMee], e = kOaE[b - kNMee + 1]; s = [a, e](const Row& r) { return r.oa >= a && r.oa < e; }; }
    if (ring == 0) return [s](const Row& r) { return r.shared_ring == 1 && s(r); };
    if (ring == 1) return [s](const Row& r) { return r.shared_ring != 1 && s(r); };
    return s;
}
inline int ClsBin(int bin, int ring) { return ring < 0 ? 1 + bin : 1 + kNMee + kNOa + 2 * bin + ring; }

struct BaseRes {
    PeakFit pm[3];          // P_main of the whole base with pol1, pol2, pol3
    double r[3];
    std::vector<SBRes> sb;  // per class, r from pol2
    std::vector<PeakFit> fit;   // per class P_main (check)
    std::vector<double> rOwn;   // r from the P_main pol2 of the class itself (variant)
};
// S of a class for another r (D linear in r)
inline double SAt(const SBRes& s, double r) { return (s.DW - r * s.DSB) / (s.fW - r * s.fSB); }
// change of S for CB x 1.03
inline double SCB3(const SBRes& s) { return -0.03 * (s.cbW - s.r * s.cbSB) / (s.fW - s.r * s.fSB); }
inline double SModel(const SBRes& s, const BaseRes& b) {
    const double a = SAt(s, b.r[0]), c = SAt(s, b.r[2]), m = SAt(s, b.r[1]);
    return 0.5 * (std::max(a, std::max(c, m)) - std::min(a, std::min(c, m)));
}

void Txt(double x, double y, const char* t, double size = 0.050, int col = kBlack, int align = 11) {
    TLatex l; l.SetNDC(); l.SetTextFont(42); l.SetTextSize(size); l.SetTextColor(col); l.SetTextAlign(align);
    l.DrawLatex(x, y, t);
}

} // namespace resme

void res_mee_cnt() {
    using namespace resme;
    TStopwatch swAll;
    gSystem->mkdir(kFigOut, kTRUE); gSystem->mkdir(kTabOut, kTRUE); gSystem->mkdir(kWork, kTRUE); gSystem->mkdir(kCntScratch, kTRUE);
    Style();

    Data d;
    {
        Data a;
        if (!scan::Load(a)) return;
        Filter(a, d, [](const Row& r) { return r.m4_fb >= 0.45f && r.m4_fb < 0.70f; });
    }
    printf("rows 0.45-0.70: %zu / %zu / %zu\n", d.s[0].size(), d.s[1].size(), d.s[2].size());

    std::vector<Base> bases = {
        {"pre", "pre", [](const Row& r) { return Pre(r); }},
        {"pre, isBest==1", "pre, isBest $=1$", [](const Row& r) { return Pre(r) && r.isBest == 1; }},
        {"peak nominal", "peak nominal", [](const Row& r) { return Pre(r) && PeakNomKin(r); }},
        {"peak nominal+vertex", "peak nominal + vertex", [](const Row& r) { return Pre(r) && PeakNomKin(r) && PiVtx(r); }},
        {"nominal, no vz cut", "nominal, no $v_z$ cut", [](const Row& r) { return PeakNomKin(r); }},
        {"pre+vertex (X17 base)", "pre + vertex", [](const Row& r) { return Pre(r) && PiVtx(r); }},
        {"pre+vertex, isBest==1", "pre + vertex, isBest $=1$", [](const Row& r) { return Pre(r) && PiVtx(r) && r.isBest == 1; }},
    };
    const int nB = (int)bases.size();
    const int iX0 = 5;   // the two X17 bases

    // ---- base fits and r ----
    std::vector<BaseRes> B(nB);
    {
        std::vector<Sel> s;
        for (const auto& b : bases) s.push_back(b.f);
        const std::vector<Spectra> sp = scan::FillMany(d, s);
        for (int b = 0; b < nB; ++b)
            for (int o = 1; o <= 3; ++o) {
                PeakCfg c; c.order = o;
                B[b].pm[o - 1] = scan::FitPeakMain(sp[b], c);
                B[b].r[o - 1] = RatioWSB(B[b].pm[o - 1]);
            }
    }
    // ---- classes: side-band subtraction and per-class fits ----
    std::vector<Sel> cs;
    for (int b = 0; b < iX0; ++b)
        for (int c = 0; c < kNCls; ++c) cs.push_back(And(bases[b].f, ClassSel(c)));
    {
        TStopwatch sw;
        const std::vector<Spectra> sp = scan::FillMany(d, cs);
        for (int b = 0; b < iX0; ++b) {
            std::vector<Sel> one(cs.begin() + b * kNCls, cs.begin() + (b + 1) * kNCls);
            B[b].sb = SideBandMany(d, one, B[b].r[1]);
            for (int c = 0; c < kNCls; ++c) {
                B[b].fit.push_back(scan::FitPeakMain(sp[b * kNCls + c]));
                B[b].rOwn.push_back(B[b].fit.back().ok ? RatioWSB(B[b].fit.back()) : B[b].r[1]);
            }
        }
        printf("classes: %.1f s\n", sw.RealTime());
    }
    // ---- X17 windows (+ OA > 2 variant) and fine bins ----
    std::vector<SBRes> xw[3], xf;    // [0] pre+vertex, [1] isBest == 1, [2] pre+vertex + OA > 2
    std::vector<PeakFit> xfit[3];
    {
        for (int v = 0; v < 3; ++v) {
            const Sel base = v == 1 ? bases[iX0 + 1].f : (v == 2 ? OA2(bases[iX0].f) : bases[iX0].f);
            std::vector<Sel> s;
            for (int k = 0; k < kNX; ++k) { const float a = kXlo[k], e = kXhi[k]; s.push_back([base, a, e](const Row& r) { return r.mee >= a && r.mee < e && base(r); }); }
            xw[v] = SideBandMany(d, s, B[v == 1 ? iX0 + 1 : iX0].r[1]);
            const std::vector<Spectra> sp = scan::FillMany(d, s);
            for (int k = 0; k < kNX; ++k) xfit[v].push_back(scan::FitPeakMain(sp[k]));
        }
        std::vector<Sel> s;
        for (int k = 0; k < kNF; ++k) {
            const float a = FineLo(k), e = FineLo(k + 1);
            const Sel base = bases[iX0].f;
            s.push_back([base, a, e](const Row& r) { return r.mee >= a && r.mee < e && base(r); });
        }
        xf = SideBandMany(d, s, B[iX0].r[1]);
    }

    // ------------------------------------------------------------------------
    // tables
    // ------------------------------------------------------------------------
    OpenTxt("res_mee.txt");
    P("# res_mee.txt -- written by plots/exp_report/res_mee.C (%s)\n", TDatime().AsSQLString());
    P("# EXPERIMENTAL DATA, p+p 4.53 GeV (feb22, 29 days). Rows = stored combinations; cor_ok == 1; FAT: vertex z > -500 mm, PT3 (trigbit 8192).\n");
    P("# Side-band subtraction in M4 of N - CB: W = %.5f-%.5f, SB = %.4f-%.5f and %.5f-%.4f GeV/c2;\n", scan::kPeakLo, scan::kPeakHi, kSb1Lo, kSb1Hi, kSb2Lo, kSb2Hi);
    P("#   D = (N-CB)_W - r (N-CB)_SB, r = integral of the P_main pol2 of the base (all classes) over W / over SB; S = D / (fW - r fSB).\n");
    P("#   CB = 2 sqrt(N++ N--) of the same class, k = 1 (at small OA / M(ee) like-sign pairs are almost absent: CB there is a lower limit).\n");
    P("#   model = half range of S for r from pol1, pol2, pol3. S_fit = P_main fit of the class alone (frozen 0.5439 / 5.4 MeV + pol2, 0.5225-0.6100).\n");
    P("#   S(r_own) = the same subtraction with r from the P_main pol2 of the class itself (variant); dS(CB+3%%) = change of S for CB x 1.03.\n");
    P("#   shares: S of a class / sum over the classes of its group; errors treat the classes as independent.\n\n");
    P("A. BASE FITS\n");
    for (int b = 0; b < nB; ++b)
        P("   %-24s P_main pol1/pol2/pol3: S = %5.0f / %5.0f +- %3.0f / %5.0f (chi2 %.1f/%.1f/%.1f), r = %.4f / %.4f / %.4f; fW %.4f, fSB %.4f\n",
          bases[b].key.c_str(), B[b].pm[0].S, B[b].pm[1].S, B[b].pm[1].eS, B[b].pm[2].S, B[b].pm[0].chi2, B[b].pm[1].chi2, B[b].pm[2].chi2,
          B[b].r[0], B[b].r[1], B[b].r[2], scan::GaussFrac(scan::kPeakLo, scan::kPeakHi, scan::kPeakMean, scan::kPeakSigma),
          scan::GaussFrac(kSb1Lo, kSb1Hi, scan::kPeakMean, scan::kPeakSigma) + scan::GaussFrac(kSb2Lo, kSb2Hi, scan::kPeakMean, scan::kPeakSigma));

    // per base: table of classes
    double shareHi[iX0], eShareHi[iX0], mShareHi[iX0], sHi[iX0], eSHi[iX0], shareHiFit[iX0], eShareHiFit[iX0];
    double oa2[iX0], eOa2[iX0], ring[iX0], eRing[iX0], shareHiOwn[iX0], shareHiCB[iX0];
    for (int b = 0; b < iX0; ++b) {
        const BaseRes& R = B[b];
        P("\nB. %s (%s)\n", bases[b].key.c_str(), b == 4 ? "earlier definition, for the comparison with 30.6 +- 3.8 %" : "isBest as stated");
        P("   %-16s | %7s %6s %6s %5s %5s | %8s %6s | %6s %5s %5s | %6s %5s | %6s %5s | %8s %6s | %-15s | %6s %5s\n", "class", "N_W", "CB_W", "CB/N", "N++", "N--",
          "(N-CB)SB", "err", "S", "stat", "model", "share", "err", "S_fit", "err", "S(r_own)", "dS(CB)", "shared ring", "S/MeV", "err");
        auto line = [&](const char* name, int c, double share, double eshare, double width) {
            const SBRes& s = R.sb[c];
            P("   %-16s | %7.0f %6.0f %6.3f %5.0f %5.0f | %8.0f %6.0f | %6.0f %5.0f %5.0f | %6.3f %5.3f | %6.0f %5.0f | %8.0f %+6.1f |", name, s.nW, s.cbW,
              s.nW > 0 ? s.cbW / s.nW : 0., s.ppW, s.mmW, s.DSB, s.eDSB, s.S, s.eS, SModel(s, R), share, eshare, R.fit[c].S, R.fit[c].eS,
              SAt(s, R.rOwn[c]), SCB3(s));
            return width;
        };
        // M(ee)
        double sum = 0, sum2 = 0;
        for (int k = 0; k < kNMee; ++k) { sum += R.sb[ClsBin(k, -1)].S; sum2 += std::pow(R.sb[ClsBin(k, -1)].eS, 2); }
        for (int k = 0; k < kNMee; ++k) {
            const SBRes& s = R.sb[ClsBin(k, -1)];
            double f, ef; Share(s.S, s.eS, sum - s.S, std::sqrt(std::max(0., sum2 - s.eS * s.eS)), f, ef);
            line((std::string("M(ee) ") + kMeeName[k]).c_str(), ClsBin(k, -1), f, ef, 0);
            double fr, efr; Share(R.sb[ClsBin(k, 0)].S, R.sb[ClsBin(k, 0)].eS, R.sb[ClsBin(k, 1)].S, R.sb[ClsBin(k, 1)].eS, fr, efr);
            const double w = 1e3 * (kMeeE[k + 1] - kMeeE[k]);
            P(" %6.3f +- %5.3f | %6.1f %5.1f\n", fr, efr, s.S / w, s.eS / w);
        }
        // share above 20 MeV (bins 3..5)
        {
            double a = 0, a2 = 0, c = 0, c2 = 0, af = 0, af2 = 0, cf = 0, cf2 = 0, am[3] = {0, 0, 0}, cm[3] = {0, 0, 0}, ao = 0, co = 0, ak = 0, ck = 0;
            for (int k = 0; k < kNMee; ++k) {
                const SBRes& s = R.sb[ClsBin(k, -1)];
                const PeakFit& f = R.fit[ClsBin(k, -1)];
                const double so = SAt(s, R.rOwn[ClsBin(k, -1)]);
                if (k >= 3) { a += s.S; a2 += s.eS * s.eS; af += f.S; af2 += f.eS * f.eS; ao += so; ak += s.S + SCB3(s); for (int o = 0; o < 3; ++o) am[o] += SAt(s, R.r[o]); }
                else { c += s.S; c2 += s.eS * s.eS; cf += f.S; cf2 += f.eS * f.eS; co += so; ck += s.S + SCB3(s); for (int o = 0; o < 3; ++o) cm[o] += SAt(s, R.r[o]); }
            }
            Share(a, std::sqrt(a2), c, std::sqrt(c2), shareHi[b], eShareHi[b]);
            Share(af, std::sqrt(af2), cf, std::sqrt(cf2), shareHiFit[b], eShareHiFit[b]);
            double mn = 1, mx = 0;
            for (int o = 0; o < 3; ++o) { const double f = am[o] / (am[o] + cm[o]); mn = std::min(mn, f); mx = std::max(mx, f); }
            mShareHi[b] = 0.5 * (mx - mn); sHi[b] = a; eSHi[b] = std::sqrt(a2);
            shareHiOwn[b] = ao / (ao + co); shareHiCB[b] = ak / (ak + ck);
            P("   share of the content above 20 MeV/c2: %.3f +- %.3f (stat) +- %.3f (model), S(>20) = %.0f +- %.0f; per-class fits: %.3f +- %.3f;\n"
              "   with r of each class: %.3f; with CB x 1.03: %.3f\n",
              shareHi[b], eShareHi[b], mShareHi[b], a, std::sqrt(a2), shareHiFit[b], eShareHiFit[b], shareHiOwn[b], shareHiCB[b]);
            P("   sum of the M(ee) classes %.0f +- %.0f against the whole base (side band) %.0f +- %.0f and P_main %.0f +- %.0f\n", a + c,
              std::sqrt(a2 + c2), R.sb[0].S, R.sb[0].eS, R.pm[1].S, R.pm[1].eS);
        }
        // OA
        sum = 0; sum2 = 0;
        for (int k = 0; k < kNOa; ++k) { sum += R.sb[ClsBin(kNMee + k, -1)].S; sum2 += std::pow(R.sb[ClsBin(kNMee + k, -1)].eS, 2); }
        double aO = 0, aO2 = 0, cO = 0, cO2 = 0;
        for (int k = 0; k < kNOa; ++k) {
            const SBRes& s = R.sb[ClsBin(kNMee + k, -1)];
            double f, ef; Share(s.S, s.eS, sum - s.S, std::sqrt(std::max(0., sum2 - s.eS * s.eS)), f, ef);
            line((std::string("OA ") + kOaName[k]).c_str(), ClsBin(kNMee + k, -1), f, ef, 0);
            double fr, efr; Share(R.sb[ClsBin(kNMee + k, 0)].S, R.sb[ClsBin(kNMee + k, 0)].eS, R.sb[ClsBin(kNMee + k, 1)].S, R.sb[ClsBin(kNMee + k, 1)].eS, fr, efr);
            const double w = kOaE[k + 1] - kOaE[k];
            P(" %6.3f +- %5.3f | %6.1f %5.1f (per deg)\n", fr, efr, s.S / w, s.eS / w);
            if (k >= 2) { aO += s.S; aO2 += s.eS * s.eS; } else { cO += s.S; cO2 += s.eS * s.eS; }
        }
        Share(aO, std::sqrt(aO2), cO, std::sqrt(cO2), oa2[b], eOa2[b]);
        P("   OA > 2 deg share of the content: %.3f +- %.3f\n", oa2[b], eOa2[b]);
        // shared ring, all
        double r1 = 0, r12 = 0, r0 = 0, r02 = 0;
        for (int k = 0; k < kNMee; ++k) { r1 += R.sb[ClsBin(k, 0)].S; r12 += std::pow(R.sb[ClsBin(k, 0)].eS, 2); r0 += R.sb[ClsBin(k, 1)].S; r02 += std::pow(R.sb[ClsBin(k, 1)].eS, 2); }
        Share(r1, std::sqrt(r12), r0, std::sqrt(r02), ring[b], eRing[b]);
        P("   shared RICH ring share of the content: %.3f +- %.3f\n", ring[b], eRing[b]);
    }
    P("\nREADING. Small M(ee) and small OA are filled by two sources that the data alone cannot separate: eta -> pi+pi-gamma with the photon\n");
    P("converted in the material (target, beam pipe, RICH radiator/mirror; M(ee) of a few MeV/c2 from the opening angle and resolution) and\n");
    P("the direct decay eta -> pi+pi-e+e- whose M(ee) rises like 1/M(ee) towards 2 m_e. Such pairs open by less than the ring size and share\n");
    P("one RICH ring (shared_ring == 1). Above 20 MeV/c2 (OA of several degrees) conversions are absent: this content is direct decays.\n");
    P("The split into decay modes needs the simulated acceptance per mode (next report).\n");

    P("\nC. X17 WINDOW (base: pre + both pions in the vertex fit; data peak window; side-band subtraction with r of the base)\n");
    const char* xv[3] = {"isBest>=0", "isBest==1", "isBest>=0, OA>2"};
    for (int v = 0; v < 3; ++v) {
        P("   %s (r = %.4f)\n", xv[v], v == 1 ? B[iX0 + 1].r[1] : B[iX0].r[1]);
        P("   %-14s | %6s %5s %5s %5s | %8s %5s | %8s %5s | %7s %5s | %6s %5s | %6s %5s\n", "M(ee) [MeV]", "N_W", "CB_W", "N++", "N--", "(N-CB)_W", "err",
          "bkg r*SB", "err", "D", "err", "S", "err", "S_fit", "err");
        for (int k = 0; k < kNX; ++k) {
            const SBRes& s = xw[v][k];
            P("   %5.1f-%5.1f%s | %6.0f %5.1f %5.0f %5.0f | %8.1f %5.1f | %8.1f %5.1f | %7.1f %5.1f | %6.1f %5.1f | %6.1f %5.1f\n", 1e3 * kXlo[k], 1e3 * kXhi[k],
              k == 1 ? " *" : "  ", s.nW, s.cbW, s.ppW, s.mmW, s.DW, s.eDW, s.bkgW(), s.r * s.eDSB, s.D, s.eD, s.S, s.eS, xfit[v][k].S, xfit[v][k].eS);
        }
    }
    P("   (* = X17 window 15.5-24.5 MeV/c2). D = rows of the eta peak inside the window after the side-band subtraction; S = D / (fW - r fSB).\n");
    {
        const double g = std::sqrt(xw[0][0].D * xw[0][2].D);
        const double eg = 0.5 * g * std::hypot(xw[0][0].eD / xw[0][0].D, xw[0][2].eD / xw[0][2].D);
        P("   Log-linear interpolation from the two neighbours (geometric mean of D(6.5-15.5) and D(24.5-33.5); the spectrum falls steeply):\n");
        P("   %.0f +- %.0f expected against %.1f +- %.1f observed in 15.5-24.5 (difference %.0f +- %.0f). This is a look, not a test.\n",
          g, eg, xw[0][1].D, xw[0][1].eD, xw[0][1].D - g, std::hypot(eg, xw[0][1].eD));
    }
    P("\n   fine M(ee) bins (3 MeV), base pre + vertex, isBest >= 0:\n");
    for (int k = 0; k < kNF; ++k)
        P("   %5.1f-%5.1f | N_W %6.0f CB_W %6.1f (N-CB)_W %7.1f +- %5.1f bkg %7.1f +- %5.1f D %7.1f +- %5.1f\n", 1e3 * FineLo(k), 1e3 * FineLo(k + 1), xf[k].nW,
          xf[k].cbW, xf[k].DW, xf[k].eDW, xf[k].bkgW(), xf[k].r * xf[k].eDSB, xf[k].D, xf[k].eD);
    P("\nrun time %.0f s\n", swAll.RealTime());
    CloseTxt();

    // ---- LaTeX: content by class ----
    {
        Tex t("res_mee.tex");
        t("%% res_mee.tex -- written by plots/exp_report/res_mee.C (numbers: tables_exp/res_mee.txt)\n");
        t("\\begin{table}[H]\n\\centering\\footnotesize\n\\setlength{\\tabcolsep}{3pt}\n");
        t("\\caption{Content of the $\\eta$ peak in classes of $M(e^+e^-)$ and OA$(e^+e^-)$, experimental data only. "
          "Side-band subtraction in $M_4$: $S_\\mathrm{SB}=[(N-\\mathrm{CB})_{R_\\mathrm{pk}}-r\\,(N-\\mathrm{CB})_{R_\\mathrm{SB}}]/(f_\\mathrm{pk}-r f_{SB})$ with the data peak window $R_\\mathrm{pk}$ = 0.53364--0.55416, "
          "side bands $R_\\mathrm{SB}$ = 0.5225--0.53364 and 0.55416--0.56530~GeV/$c^2$, $r$ from the $P_\\mathrm{main}$ pol2 of the whole selection; CB $=2\\sqrt{N_{++}N_{--}}$ of the class ($k=1$). "
          "Errors: statistical (event groups) $\\pm$ model ($r$ from pol1/pol3). Ring: share of the class content with both leptons on one RICH ring (it follows the opening angle and does not separate direct decays from conversions). "
          "$S_\\mathrm{fit}$: $P_\\mathrm{main}$ fit of the class alone. Selection pre, all stored combinations (isBest $\\geq0$); not normalised to luminosity.}\n");
        t("\\label{tab:res_mee}\n\\begin{tabular}{lrrcccc}\n\\toprule\n");
        t("Class & $N_\\mathrm{pk}$ & CB$/N_\\mathrm{pk}$ & $S_\\mathrm{SB}\\pm$stat$\\pm$model & share & ring share & $S_\\mathrm{fit}$ \\\\\n\\midrule\n");
        const BaseRes& R = B[0];
        double sum = 0, sum2 = 0;
        for (int k = 0; k < kNMee; ++k) { sum += R.sb[ClsBin(k, -1)].S; sum2 += std::pow(R.sb[ClsBin(k, -1)].eS, 2); }
        for (int k = 0; k < kNMee; ++k) {
            const SBRes& s = R.sb[ClsBin(k, -1)];
            double f, ef, fr, efr;
            Share(s.S, s.eS, sum - s.S, std::sqrt(std::max(0., sum2 - s.eS * s.eS)), f, ef);
            Share(R.sb[ClsBin(k, 0)].S, R.sb[ClsBin(k, 0)].eS, R.sb[ClsBin(k, 1)].S, R.sb[ClsBin(k, 1)].eS, fr, efr);
            t("$M_{ee}$ %s MeV/$c^2$ & %.0f & %.3f & $%.0f\\pm%.0f\\pm%.0f$ & $%.3f\\pm%.3f$ & $%.2f\\pm%.2f$ & $%.0f\\pm%.0f$ \\\\\n", kMeeName[k], s.nW,
              s.nW > 0 ? s.cbW / s.nW : 0., s.S, s.eS, SModel(s, R), f, ef, fr, efr, R.fit[ClsBin(k, -1)].S, R.fit[ClsBin(k, -1)].eS);
        }
        t("\\midrule\n");
        sum = 0; sum2 = 0;
        for (int k = 0; k < kNOa; ++k) { sum += R.sb[ClsBin(kNMee + k, -1)].S; sum2 += std::pow(R.sb[ClsBin(kNMee + k, -1)].eS, 2); }
        for (int k = 0; k < kNOa; ++k) {
            const SBRes& s = R.sb[ClsBin(kNMee + k, -1)];
            double f, ef, fr, efr;
            Share(s.S, s.eS, sum - s.S, std::sqrt(std::max(0., sum2 - s.eS * s.eS)), f, ef);
            Share(R.sb[ClsBin(kNMee + k, 0)].S, R.sb[ClsBin(kNMee + k, 0)].eS, R.sb[ClsBin(kNMee + k, 1)].S, R.sb[ClsBin(kNMee + k, 1)].eS, fr, efr);
            t("OA %s$^\\circ$ & %.0f & %.3f & $%.0f\\pm%.0f\\pm%.0f$ & $%.3f\\pm%.3f$ & $%.2f\\pm%.2f$ & $%.0f\\pm%.0f$ \\\\\n", kOaName[k], s.nW,
              s.nW > 0 ? s.cbW / s.nW : 0., s.S, s.eS, SModel(s, R), f, ef, fr, efr, R.fit[ClsBin(kNMee + k, -1)].S, R.fit[ClsBin(kNMee + k, -1)].eS);
        }
        t("\\midrule\n");
        t("\\multicolumn{7}{p{0.93\\textwidth}}{Share of the content with $M_{ee}>20$~MeV/$c^2$ (stat $\\pm$ model; in brackets: per-class fits, and with $r$ of each class):} \\\\\n");
        for (int b = 0; b < iX0; ++b)
            t("\\multicolumn{7}{p{0.93\\textwidth}}{\\quad %s: $%.3f\\pm%.3f\\pm%.3f$ ($%.3f\\pm%.3f$; own $r$: %.3f); OA$>2^\\circ$ share $%.3f\\pm%.3f$; shared ring $%.3f\\pm%.3f$} \\\\\n",
              bases[b].tex.c_str(), shareHi[b], eShareHi[b], mShareHi[b], shareHiFit[b], eShareHiFit[b], shareHiOwn[b], oa2[b], eOa2[b], ring[b], eRing[b]);
        t("\\bottomrule\n\\end{tabular}\n\\end{table}\n");
    }
    // ---- LaTeX: X17 window ----
    {
        Tex t("res_x17.tex");
        t("%% res_x17.tex -- written by plots/exp_report/res_mee.C (numbers: tables_exp/res_mee.txt)\n");
        t("\\begin{table}[H]\n\\centering\\small\n\\setlength{\\tabcolsep}{4pt}\n");
        t("\\caption{The X17 band, data only: rows in the data $\\eta$ peak window (0.53364--0.55416~GeV/$c^2$) for the selection pre + both pions in the vertex fit "
          "(frozen X17 selection with the data peak window as X17 $M_4$ window), in $M(e^+e^-)$ bands of 9~MeV/$c^2$ (bold: the X17 band). "
          "CB $=2\\sqrt{N_{++}N_{--}}$ ($k=1$); background $=r\\,(N-\\mathrm{CB})_{R_\\mathrm{SB}}$ from the $M_4$ side bands (0.5225--0.53364, 0.55416--0.56530) with $r=%.3f$; "
          "$D$ = $\\eta$-peak rows inside $R_\\mathrm{pk}$ after the subtraction. All stored combinations (isBest $\\geq0$); in brackets isBest $=1$. "
          "No limit is derived here: it needs the simulated X17 efficiency.}\n", B[iX0].r[1]);
        t("\\label{tab:res_x17}\n\\begin{tabular}{lrrrrr}\n\\toprule\n");
        t("$M_{ee}$ [MeV/$c^2$] & $N_\\mathrm{pk}$ & CB$_\\mathrm{pk}$ & $(N-\\mathrm{CB})_\\mathrm{pk}$ & background & $D$ \\\\\n\\midrule\n");
        for (int k = 0; k < kNX; ++k) {
            const SBRes& s = xw[0][k]; const SBRes& q = xw[1][k];
            t("%s%.1f--%.1f%s & %.0f (%.0f) & %.1f & $%.0f\\pm%.0f$ & $%.0f\\pm%.0f$ & $%.0f\\pm%.0f$ ($%.0f\\pm%.0f$) \\\\\n", k == 1 ? "\\textbf{" : "",
              1e3 * kXlo[k], 1e3 * kXhi[k], k == 1 ? "}" : "", s.nW, q.nW, s.cbW, s.DW, s.eDW, s.bkgW(), s.r * s.eDSB, s.D, s.eD, q.D, q.eD);
        }
        t("\\bottomrule\n\\end{tabular}\n\\end{table}\n");
    }

    // ------------------------------------------------------------------------
    // figures
    // ------------------------------------------------------------------------
    TFile fr(Form("%s/res_mee.root", kWork), "RECREATE");
    {
        TCanvas* c = new TCanvas("c_cont", "", 2000, 820);
        for (int g = 0; g < 2; ++g) {
            c->cd();
            TPad* p = new TPad(Form("pc_%d", g), "", g / 2.0, 0, (g + 1) / 2.0, 1);
            p->SetLeftMargin(0.15); p->SetRightMargin(0.03); p->SetTopMargin(0.10); p->SetBottomMargin(0.16);
            p->Draw(); p->cd();
            const int n = g == 0 ? kNMee : kNOa, off = g == 0 ? 0 : kNMee;
            const BaseRes& R = B[0];
            double sum = 0;
            for (int k = 0; k < n; ++k) sum += R.sb[ClsBin(off + k, -1)].S;
            TH1D* f = new TH1D(Form("fc_%d", g), "", n, 0, n); f->SetDirectory(nullptr); f->SetStats(0);
            for (int k = 0; k < n; ++k) f->GetXaxis()->SetBinLabel(k + 1, g == 0 ? kMeeName[k] : kOaName[k]);
            f->GetXaxis()->SetLabelSize(0.065);
            f->GetXaxis()->SetTitle(g == 0 ? "M(e^{+}e^{#minus}) [MeV/c^{2}]" : "OA(e^{+}e^{#minus}) [deg]");
            f->GetYaxis()->SetTitle("share of the #eta-peak content");
            f->SetMinimum(-0.05); f->SetMaximum(0.75);
            f->Draw("axis");
            TLine l0; l0.SetLineStyle(3); l0.SetLineColor(kGray + 1); l0.DrawLine(0, 0, n, 0);
            TGraphErrors* ga = new TGraphErrors(); TGraphErrors* gr = new TGraphErrors(); TGraphErrors* g2 = new TGraphErrors(); TGraphErrors* gf = new TGraphErrors();
            for (int k = 0; k < n; ++k) {
                const SBRes& s = R.sb[ClsBin(off + k, -1)];
                ga->SetPoint(k, k + 0.5, s.S / sum); ga->SetPointError(k, 0, std::sqrt(s.eS * s.eS + std::pow(SModel(s, R), 2)) / sum);
                gr->SetPoint(k, k + 0.35, R.sb[ClsBin(off + k, 0)].S / sum); gr->SetPointError(k, 0, R.sb[ClsBin(off + k, 0)].eS / sum);
                g2->SetPoint(k, k + 0.65, R.sb[ClsBin(off + k, 1)].S / sum); g2->SetPointError(k, 0, R.sb[ClsBin(off + k, 1)].eS / sum);
                gf->SetPoint(k, k + 0.5, R.fit[ClsBin(off + k, -1)].S / sum);
            }
            ga->SetMarkerStyle(20); ga->SetMarkerSize(1.5); ga->SetMarkerColor(kBlue + 1); ga->SetLineColor(kBlue + 1); ga->SetLineWidth(2);
            gr->SetMarkerStyle(25); gr->SetMarkerSize(1.3); gr->SetMarkerColor(kMagenta + 1); gr->SetLineColor(kMagenta + 1); gr->SetLineWidth(2);
            g2->SetMarkerStyle(26); g2->SetMarkerSize(1.3); g2->SetMarkerColor(kGreen + 2); g2->SetLineColor(kGreen + 2); g2->SetLineWidth(2);
            gf->SetMarkerStyle(5); gf->SetMarkerSize(2.0); gf->SetMarkerColor(kBlack);
            gr->Draw("pz same"); g2->Draw("pz same"); ga->Draw("pz same"); gf->Draw("p same");
            TLegend* lg = new TLegend(0.47, 0.62, 0.97, 0.88); lg->SetTextSize(0.040);
            lg->AddEntry(ga, "all (side-band subtracted)", "pe");
            lg->AddEntry(gr, "pairs on one RICH ring", "pe");
            lg->AddEntry(g2, "pairs on two RICH rings", "pe");
            lg->AddEntry(gf, "all, P_{main} fit per class", "p");
            lg->Draw();
            Txt(0.15, 0.93, g == 0 ? "(a) #eta-peak content vs M(e^{+}e^{#minus}), pre, isBest #geq 0" : "(b) #eta-peak content vs OA(e^{+}e^{#minus}), pre, isBest #geq 0", 0.046);
            Txt(g == 0 ? 0.58 : 0.55, 0.55, g == 0 ? Form("share > 20 MeV/c^{2}: %.2f #pm %.2f", shareHi[0], std::hypot(eShareHi[0], mShareHi[0]))
                                   : Form("share above 2#circ: %.2f #pm %.2f", oa2[0], eOa2[0]), 0.038, kBlue + 1);
            Txt(g == 0 ? 0.58 : 0.55, 0.49, Form("shared RICH ring: %.2f #pm %.2f", ring[0], eRing[0]), 0.038, kMagenta + 1);
        }
        Save(c, "res_mee_content");
    }
    {
        TCanvas* c = new TCanvas("c_x17", "", 2000, 820);
        for (int g = 0; g < 2; ++g) {
            c->cd();
            TPad* p = new TPad(Form("px_%d", g), "", g / 2.0, 0, (g + 1) / 2.0, 1);
            p->SetLeftMargin(0.15); p->SetRightMargin(0.03); p->SetTopMargin(0.10); p->SetBottomMargin(0.15);
            p->Draw(); p->cd();
            if (g == 0) p->SetLogy();
            TH1D* f = new TH1D(Form("fx_%d", g), "", 10, 0, 61); f->SetDirectory(nullptr); f->SetStats(0);
            f->GetXaxis()->SetTitle("M(e^{+}e^{#minus}) [MeV/c^{2}]");
            f->GetYaxis()->SetTitle("counts / 3 MeV/c^{2}");
            double mx = 0, mnD = 0, mxD = 0;
            for (int k = 0; k < kNF; ++k) { mx = std::max(mx, xf[k].nW); mnD = std::min(mnD, xf[k].D - xf[k].eD); mxD = std::max(mxD, xf[k].D + xf[k].eD); }
            if (g == 0) { f->SetMinimum(0.5); f->SetMaximum(mx * 30); } else { f->SetMinimum(mnD - 0.1 * (mxD - mnD)); f->SetMaximum(mxD + 0.5 * (mxD - mnD)); }
            f->Draw("axis");
            const double y0 = f->GetMinimum(), y1 = g == 0 ? mx * 1.5 : f->GetMinimum() + 0.62 * (f->GetMaximum() - f->GetMinimum());
            TBox* bx = new TBox(15.5, y0, 24.5, y1); bx->SetFillColorAlpha(kOrange, 0.25); bx->Draw();
            TGraphErrors* gN = new TGraphErrors(); TGraphErrors* gC = new TGraphErrors(); TGraphErrors* gB = new TGraphErrors(); TGraphErrors* gD = new TGraphErrors();
            for (int k = 0; k < kNF; ++k) {
                const double x = 1e3 * 0.5 * (FineLo(k) + FineLo(k + 1));
                const SBRes& s = xf[k];
                gN->SetPoint(k, x, s.nW); gN->SetPointError(k, 0, std::sqrt(s.nW));
                gC->SetPoint(k, x, std::max(s.cbW, 1e-3)); gC->SetPointError(k, 0, std::sqrt(s.ppW + s.mmW) * (s.cbW > 0 ? s.cbW / std::max(1., s.ppW + s.mmW) : 0));
                gB->SetPoint(k, x, s.cbW + s.bkgW()); gB->SetPointError(k, 0, s.r * s.eDSB);
                gD->SetPoint(k, x, s.D); gD->SetPointError(k, 0, s.eD);
            }
            for (auto* q : {gN, gC, gB, gD}) { q->SetMarkerSize(1.2); q->SetLineWidth(2); }
            gN->SetMarkerStyle(20); gN->SetMarkerColor(kBlack); gN->SetLineColor(kBlack);
            gC->SetMarkerStyle(20); gC->SetMarkerColor(kRed + 1); gC->SetLineColor(kRed + 1);
            gB->SetMarkerStyle(24); gB->SetMarkerColor(kGray + 2); gB->SetLineColor(kGray + 2);
            gD->SetMarkerStyle(20); gD->SetMarkerColor(kBlue + 1); gD->SetLineColor(kBlue + 1);
            if (g == 0) {
                gN->Draw("pz same"); gB->Draw("pz same"); gC->Draw("pz same");
                TLegend* lg = new TLegend(0.40, 0.70, 0.97, 0.88); lg->SetTextSize(0.038); lg->SetNColumns(1); lg->SetFillStyle(0);
                lg->AddEntry(gN, "N in the peak window", "pe");
                lg->AddEntry(gB, "CB + r (N#minusCB)_{side bands}", "pe");
                lg->AddEntry(gC, "CB = 2#sqrt{N_{++}N_{#minus#minus}}", "pe");
                lg->Draw();
                Txt(0.15, 0.93, "(a) pre + both pions in the vertex fit, isBest #geq 0", 0.044);
            } else {
                TLine l0; l0.SetLineStyle(3); l0.SetLineColor(kGray + 1); l0.DrawLine(0, 0, 61, 0);
                gD->Draw("pz same");
                Txt(0.15, 0.93, "(b) #eta-peak counts D after the side-band subtraction", 0.044);
                Txt(0.36, 0.82, Form("X17 band 15.5#minus24.5 MeV/c^{2}: D = %.0f #pm %.0f", xw[0][1].D, xw[0][1].eD), 0.038, kOrange + 7);
                Txt(0.36, 0.76, Form("6.5#minus15.5: %.0f #pm %.0f;  24.5#minus33.5: %.0f #pm %.0f", xw[0][0].D, xw[0][0].eD, xw[0][2].D, xw[0][2].eD), 0.036);
            }
        }
        Save(c, "res_x17_window");
    }
    fr.Close();
    printf("res_mee: done in %.0f s\n", swAll.RealTime());
}
