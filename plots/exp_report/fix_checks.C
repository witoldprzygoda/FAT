// ============================================================================
// fix_checks.C  --  checks requested by the review of report_exp (prefix fix)
// ============================================================================
// HADES p+p, T = 4.53 GeV (feb22, 29 days), channel pi+ pi- e+ e-. EXPERIMENTAL
// DATA ONLY. Written by the review fixer of 2026-11-14; it changes no number of
// the chapter macros, it adds the checks that the reviewers asked for.
//
// What it does
//   A. Model uncertainty of the eta-peak content S (P_main: chi2 fit of N - CB,
//      frozen Gaussian 0.5439 GeV/c2 / 5.4 MeV + pol2 in 0.5225-0.6100) for the
//      seven final candidates and four variants:
//        - the old envelope of the chapters (ScanEngine PeakModelSpread: order
//          1-4 x ranges 0.520-0.590, 0.515-0.610, 0.510-0.620, 0.505-0.650 x
//          width +-0.5 MeV, half range over the fits with p > 0.05), with its
//          rms, its background-only and width-only parts and the share of the
//          variations below the nominal S;
//        - an edge-safe model uncertainty: every lower bound >= 0.5225 (the
//          highest bump edge is 0.5202 measured, smeared by 5.4 MeV), split into
//          components: order (pol1..pol4 in the nominal range), range (lower
//          bound 0.5225/0.5250/0.5275 x upper bound 0.590/0.610/0.630/0.650,
//          pol2), width (+-0.5 MeV); each component = half range over its fits
//          with p > 0.05 (the nominal fit included), and their quadrature sum;
//          the full edge-safe grid (4 x 12 x 3 = 144 fits): its rms is the
//          quoted model uncertainty (the one-at-a-time components miss the
//          interplay of order and range);
//        - the event-group factor var/N (sum over event groups of k^2 / rows)
//          of the unlike-sign rows in the data peak window and in the fit range
//          (the fit error of P_main assumes var/N = 1).
//   B. Calibration of the statistical errors on the data: the fold halving and
//      40 random halvings of the event groups (the hash of opt_bias.C, seeds
//      1..40); for each halving the two halves are fitted separately and the
//      pull z = (X_A - X_B)/sqrt(err_A^2 + err_B^2) is formed for S (P_main)
//      and for the bump edge-step sum E (four M(pipi) slices, fixed edge, as
//      res_bump.C); rms of z over the 40 random halvings (1 if the errors are right).
//   C. P_main on N - CB of "pre" at 13 displaced positions 0.575-0.695 (the
//      test of cb_checks.C) with the upper end of each fit range: which ranges
//      cross the structure near 0.745 GeV/c2.
//   D. Share of the eta-peak content above M(ee) = 20 MeV/c2 with the side-band
//      method of res_mee.C (data peak window W, side bands, r from the P_main
//      pol2 of the base) for different side bands: both (as in the chapter),
//      upper only, lower only, upper with double width; and the same for the
//      base pre + M(pipi) > 0.33 (bump edge below the lower side band for every
//      M(ee)); per-class P_main fits for both bases.
//   E. Efficiency of MM4 < 2.85 given dz_pi < 40 mm on the data peak (pass/fail
//      P_main fits of DIRECT rows).
//
// Selections: as res_common.h (pre = cor_ok && isBest >= 0 && vz > -160 mm; in
// FAT already vz > -500 mm and PT3). CB = 2 sqrt(N++ N--) per 10 MeV, k = 1.
//
// Inputs : /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/scan_v2/scan_compact.root
// Outputs: /home/przygoda/HADES/FAT/FAT/docs/latex/tables_exp/fix_checks.txt, fix_spread.tex, fix_stat.tex
//          /home/przygoda/HADES/FAT/FAT/docs/latex/figs_exp/fix_spread.pdf (+ .png)
//
// Run (from /home/przygoda/HADES/FAT/FAT, about 2 min):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/aclic_fix", kTRUE)' 'plots/exp_report/fix_checks.C+'
// ============================================================================
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TBox.h"
#include "TMath.h"
#include "TPad.h"
#include "TStopwatch.h"

#include "res_common.h"

namespace fixc {
using namespace res;
using scan::PeakCfg;
using scan::PeakFit;
using scan::StepCfg;
using scan::StepFit;

inline int RndHalf(int grp, unsigned seed) {   // as opt_bias.C / scan_peak.C
    unsigned x = (unsigned)grp * 2654435761u + seed * 0x9E3779B9u;
    x ^= x >> 16; x *= 0x85ebca6bu; x ^= x >> 13; x *= 0xc2b2ae35u; x ^= x >> 16;
    return (int)(x & 1u);
}

const int kNLo = 3, kNHi = 4;
const double kLoS[kNLo] = {0.5225, 0.5250, 0.5275};
const double kHiS[kNHi] = {0.590, 0.610, 0.630, 0.650};

struct Model {
    PeakFit nom;
    scan::Spread old;
    double oldBelow = 0;          // share of the passing old variations with S < S_nom
    double ordHalf = 0; int ordPass = 0;
    double rngHalf = 0; int rngPass = 0;
    double widHalf = 0;
    double quad = 0;
    double gHalf = 0, gRms = 0, gMean = 0, gBelow = 0; int gPass = 0, gN = 0;
    double vnWin = 1, vnFit = 1;  // var/N of the unlike-sign rows
    std::vector<double> oLo, oS, gLo, gS;   // for the figure: lower bound and S of passing fits
};

inline Model Evaluate(const Spectra& sp) {
    Model m;
    m.nom = scan::FitPeakMain(sp);
    m.old = scan::PeakModelSpread(sp, PeakCfg(), true);
    {
        int nb = 0, np = 0;
        for (const auto& f : m.old.fits) {
            if (!f.ok || !(f.prob > scan::kVarPmin)) continue;
            np++; if (f.S < m.nom.S) nb++;
            m.oLo.push_back(f.cfg.lo); m.oS.push_back(f.S);
        }
        m.oldBelow = np > 0 ? (double)nb / np : 0;
    }
    // components (nominal fit always included)
    {
        double mn = m.nom.S, mx = m.nom.S;
        for (int o = 1; o <= 4; ++o) {
            PeakCfg c; c.order = o;
            const PeakFit f = scan::FitPeakMain(sp, c);
            if (!f.ok || !(f.prob > scan::kVarPmin)) continue;
            m.ordPass++; mn = std::min(mn, f.S); mx = std::max(mx, f.S);
        }
        m.ordHalf = 0.5 * (mx - mn);
    }
    {
        double mn = m.nom.S, mx = m.nom.S;
        for (int a = 0; a < kNLo; ++a)
            for (int b = 0; b < kNHi; ++b) {
                PeakCfg c; c.lo = kLoS[a]; c.hi = kHiS[b];
                const PeakFit f = scan::FitPeakMain(sp, c);
                if (!f.ok || !(f.prob > scan::kVarPmin)) continue;
                m.rngPass++; mn = std::min(mn, f.S); mx = std::max(mx, f.S);
            }
        m.rngHalf = 0.5 * (mx - mn);
    }
    {
        PeakCfg a, b; a.sigma -= scan::kVarDSigma; b.sigma += scan::kVarDSigma;
        m.widHalf = 0.5 * std::fabs(scan::FitPeakMain(sp, b).S - scan::FitPeakMain(sp, a).S);
    }
    m.quad = std::sqrt(m.ordHalf * m.ordHalf + m.rngHalf * m.rngHalf + m.widHalf * m.widHalf);
    {
        double mn = 0, mx = 0, s1 = 0, s2 = 0; int nb = 0;
        for (int o = 1; o <= 4; ++o)
            for (int a = 0; a < kNLo; ++a)
                for (int b = 0; b < kNHi; ++b)
                    for (int w = -1; w <= 1; ++w) {
                        PeakCfg c; c.order = o; c.lo = kLoS[a]; c.hi = kHiS[b]; c.sigma += w * scan::kVarDSigma;
                        const PeakFit f = scan::FitPeakMain(sp, c);
                        m.gN++;
                        if (!f.ok || !(f.prob > scan::kVarPmin)) continue;
                        if (m.gPass == 0) mn = mx = f.S;
                        mn = std::min(mn, f.S); mx = std::max(mx, f.S);
                        s1 += f.S; s2 += f.S * f.S; m.gPass++;
                        if (f.S < m.nom.S) nb++;
                        m.gLo.push_back(c.lo); m.gS.push_back(f.S);
                    }
        if (m.gPass > 0) {
            m.gHalf = 0.5 * (mx - mn); m.gMean = s1 / m.gPass;
            m.gRms = std::sqrt(std::max(0., s2 / m.gPass - m.gMean * m.gMean));
            m.gBelow = (double)nb / m.gPass;
        }
    }
    const scan::WinSum& w0 = sp.win[0];
    m.vnWin = w0.n[0] > 0 ? w0.var[0] / w0.n[0] : 1;
    const scan::WinSum& w2 = sp.win[2];
    m.vnFit = w2.n[0] > 0 ? w2.var[0] / w2.n[0] : 1;
    return m;
}

// ---- edge steps (copy of res_bump.C, slices 0..3) ---------------------------------
const int kNSig = 4;
const double kSlLo[kNSig] = {0.30, 0.34, 0.37, 0.39};
const double kSlHi[kNSig] = {0.34, 0.37, 0.39, 0.40};
inline StepCfg SliceCfg(int i) {
    StepCfg c; c.mpLo = kSlLo[i]; c.mpHi = kSlHi[i];
    if (i == 3) c.halfRange = 0.035;
    return c;
}

// ---- generalised side-band subtraction (res_common.h section D with free side bands)
struct SB2 { double S = 0, eS = 0, D = 0, eD = 0; };
inline std::vector<SB2> SideBandGen(const Data& d, const std::vector<Sel>& sels, const PeakFit& base,
                                     const std::vector<std::pair<double, double>>& sbs) {
    const double wLo = scan::kPeakLo, wHi = scan::kPeakHi;
    auto inSB = [&](double m) { for (const auto& p : sbs) if (m >= p.first && m < p.second) return true; return false; };
    double pw = PolyInt(base, wLo, wHi), ps = 0, fSB = 0;
    for (const auto& p : sbs) { ps += PolyInt(base, p.first, p.second); fSB += scan::GaussFrac(p.first, p.second, scan::kPeakMean, scan::kPeakSigma); }
    const double r = pw / ps;
    const double fW = scan::GaussFrac(wLo, wHi, scan::kPeakMean, scan::kPeakSigma);
    const int ns = (int)sels.size();
    Spectra g; g.nbin = scan::kNbin; g.lo = scan::kMlo; g.hi = scan::kMhi; g.bw = scan::kBinW; g.reb = scan::kCBreb; g.nw = scan::kNbin / scan::kCBreb;
    const int nw = g.nw;
    std::vector<double> cW(nw, 0.), cS(nw, 0.);
    for (int i = 0; i < g.nbin; ++i) {
        const double uW = g.Frac(i, wLo, wHi);
        double uS = 0; for (const auto& p : sbs) uS += g.Frac(i, p.first, p.second);
        if (uW <= 0 && uS <= 0) continue;
        int j0, j1; double f; g.Interp(i, j0, j1, f);
        cW[j0] += uW * (1 - f) / g.reb; cW[j1] += uW * f / g.reb;
        cS[j0] += uS * (1 - f) / g.reb; cS[j1] += uS * f / g.reb;
    }
    std::vector<double> aD(ns, 0.), vD(ns, 0.), nW(ns, 0.), nS(ns, 0.);
    std::vector<std::vector<double>> wl[2];
    std::vector<std::vector<std::pair<int, int>>> lsr[2];
    for (int q = 0; q < 2; ++q) { wl[q].assign(ns, std::vector<double>(nw, 0.)); lsr[q].assign(ns, {}); }
    for (int k = 0; k < 3; ++k) {
        const std::vector<Row>& v = d.s[k];
        int cur = -1;
        for (size_t i = 0; i < v.size(); ++i) {
            const Row& row = v[i];
            if (!row.cor_ok) continue;
            if (k == 0 && row.grp != cur) {
                for (int s = 0; s < ns; ++s) { vD[s] += aD[s] * aD[s]; aD[s] = 0; }
                cur = row.grp;
            }
            const double m = row.m4;
            if (!(m >= 0.50 && m < 0.60)) continue;
            const bool inW = m >= wLo && m < wHi;
            const bool inS = inSB(m);
            const int jb = (int)((m - scan::kMlo) / (g.reb * g.bw));
            for (int s = 0; s < ns; ++s) {
                if (!sels[s](row)) continue;
                if (k == 0) {
                    if (inW) { nW[s] += 1; aD[s] += 1; }
                    else if (inS) { nS[s] += 1; aD[s] -= r; }
                } else if (jb >= 0 && jb < nw) { wl[k - 1][s][jb] += 1; lsr[k - 1][s].push_back(std::make_pair(row.grp, jb)); }
            }
        }
        if (k == 0) for (int s = 0; s < ns; ++s) { vD[s] += aD[s] * aD[s]; aD[s] = 0; }
    }
    std::vector<SB2> out(ns);
    for (int s = 0; s < ns; ++s) {
        std::vector<double> cb(nw, 0.), cD(nw, 0.);
        double cbW = 0, cbS = 0;
        for (int j = 0; j < nw; ++j) { cb[j] = 2.0 * std::sqrt(wl[0][s][j] * wl[1][s][j]); cbW += cW[j] * cb[j]; cbS += cS[j] * cb[j]; cD[j] = cW[j] - r * cS[j]; }
        double var = 0;
        for (int q = 0; q < 2; ++q) {
            std::vector<double> wt(nw, 0.);
            for (int j = 0; j < nw; ++j) {
                const double own = wl[q][s][j], oth = wl[1 - q][s][j];
                wt[j] = cD[j] * ((own > 0 && oth > 0) ? std::sqrt(oth / own) : 1.0);
            }
            int cur = -1; double acc = 0;
            for (const auto& e : lsr[q][s]) { if (e.first != cur) { var += acc * acc; acc = 0; cur = e.first; } acc += wt[e.second]; }
            var += acc * acc;
        }
        out[s].D = (nW[s] - cbW) - r * (nS[s] - cbS);
        out[s].eD = std::sqrt(vD[s] + var);
        out[s].S = out[s].D / (fW - r * fSB); out[s].eS = out[s].eD / (fW - r * fSB);
    }
    return out;
}

void Txt(double x, double y, const char* t, double size = 0.05, int col = kBlack, int align = 11) {
    TLatex l; l.SetNDC(); l.SetTextFont(42); l.SetTextSize(size); l.SetTextColor(col); l.SetTextAlign(align);
    l.DrawLatex(x, y, t);
}
} // namespace fixc

void fix_checks() {
    using namespace fixc;
    TStopwatch sw;
    gSystem->mkdir(kFigOut, kTRUE); gSystem->mkdir(kTabOut, kTRUE);
    gSystem->mkdir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/fix", kTRUE);
    Style();

    Data d;
    {
        Data all;
        if (!scan::Load(all)) return;
        Filter(all, d, [](const Row& r) { return Pre(r) && r.m4_fb >= 0.40f && r.m4_fb < 0.80f; });
    }
    printf("rows (pre, 0.40-0.80): %zu / %zu / %zu\n", d.s[0].size(), d.s[1].size(), d.s[2].size());
    OpenTxt("fix_checks.txt");
    P("# fix_checks.txt -- written by plots/exp_report/fix_checks.C (review fixes, 2026-11-14). EXPERIMENTAL DATA ONLY.\n");
    P("# Base pre = cor_ok && isBest >= 0 && vz > -160 mm (FAT: vz > -500 mm, PT3). CB = 2 sqrt(N++ N--) per 10 MeV, interpolated, k = 1.\n");
    P("# P_main = chi2 fit of N - CB, frozen Gaussian 0.5439 GeV/c2 / 5.4 MeV + pol2 in 0.5225-0.6100, S = whole Gaussian.\n\n");

    // ------------------------------------------------------------------ A
    std::vector<NSel> cand = {SelPre(), SelPreVtx(), SelPeakNom(), SelPeakNomVtx(), SelDIRECT(), SelPEAK(), SelX17()};
    const int nc = (int)cand.size();
    std::vector<Sel> sels;
    for (int v = 0; v < 4; ++v) for (int i = 0; i < nc; ++i) sels.push_back(Variant(cand[i].f, v));
    scan::FillOpt fo; fo.windows.push_back(std::make_pair(scan::kFitLo, scan::kFitHi));
    std::vector<Spectra> sp = scan::FillMany(d, sels, fo);
    std::vector<Model> M(sels.size());
    for (size_t k = 0; k < sels.size(); ++k) M[k] = Evaluate(sp[k]);

    P("A. MODEL UNCERTAINTY OF S. old = envelope of the chapters (48 variations, lower bounds 0.505-0.520 cross the bump edge 0.5202);\n");
    P("   new = edge-safe components (all lower bounds >= 0.5225): order pol1..4 | range lo 0.5225/0.5250/0.5275 x hi 0.590/0.610/0.630/0.650 |\n");
    P("   width +-0.5 MeV; each = half range over its fits with p > 0.05 (nominal included); quad = quadrature sum; quoted model uncertainty = rms of the full edge-safe grid.\n");
    P("   grid = all 144 edge-safe fits (4 orders x 12 ranges x 3 widths). below = share of the passing variations with S < S(nominal).\n");
    P("   var/N = sum over event groups of k^2 / rows (unlike-sign) in the peak window and in the fit range 0.5225-0.6100; stat x sqrt(var/N)\n");
    P("   is the fit error with the event-group term. F = S/stat, Fg = S/(stat sqrt(var/N_win)), Ftot_o (old envelope), Ftot_n (quadrature sum) = S/sqrt(stat_g^2 + model^2).\n\n");
    P("%-24s %-16s | %6s %5s %5s | old: %5s %4s %5s %5s %5s %5s | new: %4s(%2s) %4s(%2s) %4s %5s | grid: %5s %5s %3s/%3s %5s | %6s %6s | %6s %6s %6s %6s\n",
      "selection", "variant", "S", "stat", "chi2", "half", "n", "rms", "bkg", "width", "below", "ord", "n", "rng", "n", "wid", "quad",
      "half", "rms", "np", "nf", "below", "vN_win", "vN_fit", "F", "Fg", "Ftot_o", "Ftot_n");
    for (int v = 0; v < 4; ++v)
        for (int i = 0; i < nc; ++i) {
            const Model& m = M[v * nc + i];
            const double eg = m.nom.eS * std::sqrt(m.vnWin);
            P("%-24s %-16s | %6.0f %5.0f %5.1f | old: %5.0f %4d %5.0f %5.0f %5.0f %5.2f | new: %4.0f(%2d) %4.0f(%2d) %4.0f %5.0f | grid: %5.0f %5.0f %3d/%3d %5.2f | %6.3f %6.3f | %6.2f %6.2f %6.2f %6.2f\n",
              cand[i].key.c_str(), VariantName(v), m.nom.S, m.nom.eS, m.nom.chi2, m.old.half, m.old.npass, m.old.rms, m.old.halfBkg, m.old.halfWidth,
              m.oldBelow, m.ordHalf, m.ordPass, m.rngHalf, m.rngPass, m.widHalf, m.quad, m.gHalf, m.gRms, m.gPass, m.gN, m.gBelow,
              m.vnWin, m.vnFit, m.nom.F, m.nom.S / eg, m.nom.S / std::sqrt(eg * eg + m.old.half * m.old.half),
              m.nom.S / std::sqrt(eg * eg + m.quad * m.quad));
        }
    P("\n");

    // quoted numbers: model = rms of the edge-safe grid; stat_g = fit error x sqrt(var/N) of the peak window
    P("A2. QUOTED NUMBERS (model = rms of the edge-safe grid; stat_g = stat x sqrt(var/N_win)); Y = S x 64/63 / 6079.2 [rows per nb^-1]\n");
    P("%-24s %-16s | %6s %5s %5s %5s | %6s %6s | %7s %7s %7s\n", "selection", "variant", "S", "stat", "stat_g", "model", "F_g", "Ftot", "Y", "Y stat", "Y model");
    for (int v = 0; v < 4; ++v)
        for (int i = 0; i < nc; ++i) {
            const Model& m = M[v * nc + i];
            const double eg = m.nom.eS * std::sqrt(m.vnWin);
            P("%-24s %-16s | %6.0f %5.0f %5.0f %5.0f | %6.2f %6.2f | %7.4f %7.4f %7.4f\n", cand[i].key.c_str(), VariantName(v), m.nom.S, m.nom.eS, eg, m.gRms,
              m.nom.S / eg, m.nom.S / std::sqrt(eg * eg + m.gRms * m.gRms), Yield(m.nom.S), Yield(eg), Yield(m.gRms));
        }
    P("\n");

    // X17 band in three M4 windows: N, CB, pol2 B and the Gaussian part inside (P_main of the band)
    {
        const Model& m = M[6];
        const double wl[3][2] = {{0.5333, 0.5587}, {0.5312, 0.5566}, {scan::kPeakLo, scan::kPeakHi}};
        const char* wn[3] = {"simulation window", "moved to the data scale", "data peak window"};
        P("A3. X17 BAND (pre + both pions in vertex fit + 15.5 <= M(ee) < 24.5 MeV/c2, isBest >= 0) in three M4 windows; P_main S = %.0f +- %.0f\n", m.nom.S, m.nom.eS);
        for (int w = 0; w < 3; ++w) {
            PeakCfg c; c.winLo = wl[w][0]; c.winHi = wl[w][1];
            const PeakFit f = scan::FitPeakMain(sp[6], c);
            P("   %-24s %.5f-%.5f: N = %.0f +- %.0f, CB = %.1f +- %.1f, N-CB = %.0f, pol2 B = %.0f, S inside = %.0f, S/(B+CB) = %.4f\n",
              wn[w], wl[w][0], wl[w][1], f.N, f.eN, f.CB, f.eCB, f.N - f.CB, f.B, f.Swin, f.SoB);
        }
        P("\n");
    }

    // LaTeX table of the model uncertainty (isBest >= 0 main, plus isBest == 1)
    {
        Tex t("fix_spread.tex");
        t("%% written by plots/exp_report/fix_checks.C (review fixes)\n");
        t("\\begin{table}[H]\\centering\\footnotesize\\setlength{\\tabcolsep}{3pt}\n");
        t("\\caption{Model uncertainty of the $\\eta$-peak content $S$ (\\Pmain), experimental data, all combinations (isBest $\\geq0$, top) and isBest $=1$ (bottom), no OA cut. "
          "\\emph{Old envelope} (chapter tables): half range of $S$ over the 48 variations (order 1--4, four ranges with lower bounds 0.505--0.520\\,\\GeVcc, width $\\pm0.5$\\,\\MeVcc) "
          "with $p>0.05$ ($n$ fits), its rms and the share of these fits below the nominal $S$. \\emph{Edge-safe} (every lower bound $\\geq0.5225$\\,\\GeVcc): half range over the fits with $p>0.05$ "
          "(nominal included) of the polynomial order alone (pol1--pol4, nominal range), of the range alone (lower bound 0.5225/0.5250/0.5275 $\\times$ upper bound 0.590/0.610/0.630/0.650, pol2) "
          "and of the width alone ($\\pm0.5$\\,\\MeVcc), their quadrature sum, and the half range and rms of the full grid of 144 combinations ($n$ fits with $p>0.05$). "
          "\\textbf{Quoted model uncertainty $=$ rms of the edge-safe grid.} stat$_g$ $=$ fit error $\\times\\sqrt{\\mathrm{var}/N}$ (event-group factor of the unlike-sign rows in the data peak window); "
          "$F_{\\mathrm{tot}}=S/\\sqrt{\\mathrm{stat}_g^2+\\mathrm{model}^2}$. CB with $k=1$; rows, not normalised to luminosity.}\n");
        t("\\label{tab:fix_spread}\n\\resizebox{\\textwidth}{!}{%%\n\\begin{tabular}{@{}l r r r | r r r | r r r r r r | r@{}}\n\\toprule\n");
        t(" & & & & \\multicolumn{3}{c|}{old envelope} & \\multicolumn{6}{c|}{edge-safe} & \\\\\n");
        t("selection & $S$ & stat & stat$_g$ & half ($n$) & rms & below & order & range & width & quad.\\ & grid half ($n$) & \\textbf{model} (rms) & $F_{\\mathrm{tot}}$ \\\\\n\\midrule\n");
        const char* tn[7] = {"pre", "pre + vertex", "peak nominal", "peak nominal + vertex", "DIRECT", "PEAK", "X17 band"};
        for (int v = 0; v < 2; ++v) {
            if (v == 1) t("\\midrule\n\\multicolumn{14}{@{}l}{isBest $=1$} \\\\\n");
            for (int i = 0; i < nc; ++i) {
                const Model& m = M[v * nc + i];
                const double eg = m.nom.eS * std::sqrt(m.vnWin);
                t("%s & %.0f & %.0f & %.0f & %.0f (%d) & %.0f & %.2f & %.0f & %.0f & %.0f & %.0f & %.0f (%d) & \\textbf{%.0f} & %.2f \\\\\n",
                  tn[i], m.nom.S, m.nom.eS, eg, m.old.half, m.old.npass, m.old.rms, m.oldBelow,
                  m.ordHalf, m.rngHalf, m.widHalf, m.quad, m.gHalf, m.gPass, m.gRms, m.nom.S / std::sqrt(eg * eg + m.gRms * m.gRms));
            }
        }
        t("\\bottomrule\n\\end{tabular}}\n\\end{table}\n");
    }

    // ------------------------------------------------------------------ B
    {
        std::vector<NSel> hs = {SelPre(), SelPreVtx(), SelPeakNomVtx(), SelDIRECT(), SelPEAK()};
        std::vector<NSel> he = {SelPre(), SelBumpNomVtx()};
        const int nH = 41;   // 0 = fold, 1..40 random seeds
        auto half = [](int h, int a) -> Sel {
            if (h == 0) return [a](const Row& r) { return (int)r.fold == a; };
            const unsigned seed = (unsigned)h;
            return [seed, a](const Row& r) { return RndHalf(r.grp, seed) == a; };
        };
        std::vector<Sel> s2;
        for (int h = 0; h < nH; ++h)
            for (int a = 0; a < 2; ++a) {
                const Sel hf = half(h, a);
                for (const auto& x : hs) s2.push_back(And(x.f, hf));
                for (const auto& x : he)
                    for (int i = 0; i < kNSig; ++i) {
                        const double lo = kSlLo[i], hi = kSlHi[i];
                        s2.push_back(And(And(x.f, hf), [lo, hi](const Row& r) { return r.mpipi >= lo && r.mpipi < hi; }));
                    }
            }
        printf("halvings: %zu spectra\n", s2.size());
        std::vector<Spectra> h2 = scan::FillMany(d, s2);
        const int per = (int)hs.size() + (int)he.size() * kNSig;
        std::vector<std::vector<double>> zS(hs.size()), zE(he.size());
        std::vector<double> zSf(hs.size()), zEf(he.size());
        std::vector<double> sA(hs.size()), sB(hs.size()), eA(he.size()), eB(he.size());
        for (int h = 0; h < nH; ++h) {
            double S[2][8], eS[2][8], E[2][4], eE[2][4];
            for (int a = 0; a < 2; ++a) {
                const int base = (h * 2 + a) * per;
                for (size_t j = 0; j < hs.size(); ++j) { const PeakFit f = scan::FitPeakMain(h2[base + j]); S[a][j] = f.S; eS[a][j] = f.eS; }
                for (size_t j = 0; j < he.size(); ++j) {
                    double e = 0, v = 0;
                    for (int i = 0; i < kNSig; ++i) {
                        const StepFit f = scan::FitEdgeStep(h2[base + hs.size() + j * kNSig + i], SliceCfg(i));
                        if (!f.ok) continue;
                        e += f.h; v += f.eh * f.eh;
                    }
                    E[a][j] = e; eE[a][j] = std::sqrt(v);
                }
            }
            for (size_t j = 0; j < hs.size(); ++j) {
                const double z = (S[0][j] - S[1][j]) / std::sqrt(eS[0][j] * eS[0][j] + eS[1][j] * eS[1][j]);
                if (h == 0) zSf[j] = z; else zS[j].push_back(z);
                if (h == 0) { sA[j] = S[0][j]; sB[j] = S[1][j]; }
            }
            for (size_t j = 0; j < he.size(); ++j) {
                const double z = (E[0][j] - E[1][j]) / std::sqrt(eE[0][j] * eE[0][j] + eE[1][j] * eE[1][j]);
                if (h == 0) zEf[j] = z; else zE[j].push_back(z);
                if (h == 0) { eA[j] = E[0][j]; eB[j] = E[1][j]; }
            }
        }
        auto stat = [](const std::vector<double>& z, double& mean, double& rms, double& rms0) {
            double s1 = 0, s2 = 0; for (double x : z) { s1 += x; s2 += x * x; }
            mean = s1 / z.size(); rms0 = std::sqrt(s2 / z.size()); rms = std::sqrt(std::max(0., s2 / z.size() - mean * mean));
        };
        P("B. CALIBRATION OF THE STATISTICAL ERRORS WITH INDEPENDENT HALVES (pull z = (X_A - X_B)/sqrt(err_A^2 + err_B^2), fit errors as quoted).\n");
        P("   fold = parity of the event-group index; 40 random halvings (hash of opt_bias.C, seeds 1..40). rms0 = sqrt(mean z^2) (expected 1),\n");
        P("   its statistical uncertainty about 1/sqrt(2 x 40) = 0.11 if the halvings were independent.\n");
        P("%-24s | %-5s | %8s %8s %6s | 40 random: %6s %6s %6s | n(|z|>2)\n", "selection", "X", "fold A", "fold B", "z", "mean", "rms", "rms0");
        Tex t("fix_stat.tex");
        t("%% written by plots/exp_report/fix_checks.C (review fixes)\n");
        t("\\begin{table}[H]\\centering\\small\n");
        t("\\caption{Calibration of the statistical errors on the data with independent halves of the event groups, all combinations (isBest $\\geq0$). "
          "For each halving both halves are fitted alone; pull $z=(X_A-X_B)/\\sqrt{\\sigma_A^2+\\sigma_B^2}$ with the fit errors as quoted ($S$: \\Pmain; $E$: sum of the four "
          "edge steps). Fold: parity of the group index; 40 random halvings with the hash of the bias study (seeds 1--40). If the errors are right, "
          "$\\sqrt{\\langle z^2\\rangle}=1$ (uncertainty about 0.11). CB with $k=1$; rows.}\n\\label{tab:fix_stat}\n");
        t("\\begin{tabular}{@{}l l r r r r@{}}\n\\toprule\nselection & $X$ & $z$ (fold) & $\\langle z\\rangle$ (40) & $\\sqrt{\\langle z^2\\rangle}$ (40) & $|z|>2$ \\\\\n\\midrule\n");
        const char* hn[5] = {"pre", "pre + vertex", "peak nominal + vertex", "DIRECT", "PEAK"};
        for (size_t j = 0; j < hs.size(); ++j) {
            double mean, rms, rms0; stat(zS[j], mean, rms, rms0);
            int n2 = 0; for (double x : zS[j]) if (std::fabs(x) > 2) n2++;
            P("%-24s | %-5s | %8.0f %8.0f %6.2f | %17.2f %6.2f %6.2f | %d\n", hs[j].key.c_str(), "S", sA[j], sB[j], zSf[j], mean, rms, rms0, n2);
            t("%s & $S$ & $%+.2f$ & $%+.2f$ & %.2f & %d \\\\\n", hn[j], zSf[j], mean, rms0, n2);
        }
        const char* en[2] = {"pre", "bump nominal + vertex"};
        for (size_t j = 0; j < he.size(); ++j) {
            double mean, rms, rms0; stat(zE[j], mean, rms, rms0);
            int n2 = 0; for (double x : zE[j]) if (std::fabs(x) > 2) n2++;
            P("%-24s | %-5s | %8.0f %8.0f %6.2f | %17.2f %6.2f %6.2f | %d\n", he[j].key.c_str(), "E", eA[j], eB[j], zEf[j], mean, rms, rms0, n2);
            t("%s & $E$ & $%+.2f$ & $%+.2f$ & %.2f & %d \\\\\n", en[j], zEf[j], mean, rms0, n2);
        }
        t("\\bottomrule\n\\end{tabular}\n\\end{table}\n");
        P("   per-halving pulls, pre S : "); for (double x : zS[0]) P("%.2f ", x); P("\n");
        P("   per-halving pulls, pre E : "); for (double x : zE[0]) P("%.2f ", x); P("\n\n");
    }

    // ------------------------------------------------------------------ C
    {
        const Spectra& s0 = sp[0];
        P("C. P_main ON N - CB OF pre AT DISPLACED POSITIONS (fit range and window moved with the Gaussian; as cb_checks.C).\n");
        P("   position  range_lo  range_hi    S    stat   pull  chi2/ndf\n");
        double a1 = 0, a2 = 0, b1 = 0, b2 = 0; int na = 0, nb = 0;
        for (int i = 0; i < 13; ++i) {
            const double pos = 0.575 + 0.010 * i, sh = pos - scan::kPeakMean;
            PeakCfg c; c.mean = pos; c.lo = scan::kFitLo + sh; c.hi = scan::kFitHi + sh; c.winLo = scan::kPeakLo + sh; c.winHi = scan::kPeakHi + sh;
            const PeakFit f = scan::FitPeakMain(s0, c);
            const double p = f.eS > 0 ? f.S / f.eS : 0;
            P("   %.3f   %.4f   %.4f  %6.0f %5.0f %6.2f  %5.1f/%d%s\n", pos, c.lo, c.hi, f.S, f.eS, p, f.chi2, f.ndf, c.hi > 0.745 ? "   range crosses 0.745" : "");
            if (c.hi <= 0.745) { a1 += p; a2 += p * p; na++; } else { b1 += p; b2 += p * p; nb++; }
        }
        P("   ranges below 0.745 (%d positions): pull mean %.2f, rms %.2f; ranges crossing 0.745 (%d): mean %.2f, rms %.2f\n\n",
          na, a1 / na, std::sqrt(a2 / na), nb, nb ? b1 / nb : 0., nb ? std::sqrt(b2 / nb) : 0.);
        // the like-sign / N - CB spectrum around the structure: local slope change
        P("   N - CB of pre in 10 MeV bins 0.70-0.80 (rows): ");
        for (int j = 0; j < 10; ++j) {
            double s = 0; for (int i = 0; i < s0.nbin; ++i) { const double x = s0.x(i); if (x >= 0.70 + 0.01 * j && x < 0.71 + 0.01 * j) s += s0.sig[i]; }
            P("%.0f ", s);
        }
        P("\n\n");
    }

    // ------------------------------------------------------------------ D
    {
        const int nMee = 6;
        const double mE[nMee + 1] = {0.0, 0.005, 0.010, 0.020, 0.050, 0.100, 0.270};
        P("D. SHARE OF THE ETA-PEAK CONTENT WITH M(ee) > 20 MeV/c2, side-band method (res_mee.C) with different side bands; W = 0.53364-0.55416.\n");
        P("   share = sum S(classes 20-50, 50-100, 100-270) / sum S(all six classes); errors treat the classes as independent (as res_mee.C).\n");
        for (int b = 0; b < 2; ++b) {
            const Sel base = b == 0 ? Sel([](const Row& r) { return Pre(r); }) : Sel([](const Row& r) { return Pre(r) && r.mpipi > 0.33f; });
            std::vector<Sel> cl;
            for (int k = 0; k < nMee; ++k) { const float lo = mE[k], hi = mE[k + 1]; cl.push_back(And(base, [lo, hi](const Row& r) { return r.mee >= lo && r.mee < hi; })); }
            std::vector<Sel> all1(1, base);
            const Spectra bs = scan::FillMany(d, all1)[0];
            const PeakFit bf = scan::FitPeakMain(bs);
            std::vector<Spectra> cs = scan::FillMany(d, cl);
            // per-class P_main fits
            double fa = 0, fc = 0, fva = 0, fvc = 0;
            for (int k = 0; k < nMee; ++k) {
                const PeakFit f = scan::FitPeakMain(cs[k]);
                if (k >= 3) { fa += f.S; fva += f.eS * f.eS; } else { fc += f.S; fvc += f.eS * f.eS; }
            }
            auto share = [](double a, double va, double c, double vc, double& sh, double& e) {
                sh = a / (a + c); e = std::sqrt(c * c * va + a * a * vc) / std::pow(a + c, 2);
            };
            double shf, eshf; share(fa, fva, fc, fvc, shf, eshf);
            P("   base %-22s: P_main S = %.0f +- %.0f; per-class fits: share %.3f +- %.3f (sum of classes %.0f)\n",
              b == 0 ? "pre" : "pre + M(pipi) > 0.33", bf.S, bf.eS, shf, eshf, fa + fc);
            struct SBdef { const char* name; std::vector<std::pair<double, double>> sb; };
            std::vector<SBdef> defs = {
                {"both (chapter)", {{0.5225, scan::kPeakLo}, {scan::kPeakHi, 0.56530}}},
                {"upper only", {{scan::kPeakHi, 0.56530}}},
                {"lower only", {{0.5225, scan::kPeakLo}}},
                {"upper, double width", {{scan::kPeakHi, 0.57644}}}};
            for (const auto& df : defs) {
                const std::vector<SB2> r = SideBandGen(d, cl, bf, df.sb);
                double a = 0, c = 0, va = 0, vc = 0;
                for (int k = 0; k < nMee; ++k) { if (k >= 3) { a += r[k].S; va += r[k].eS * r[k].eS; } else { c += r[k].S; vc += r[k].eS * r[k].eS; } }
                double sh, esh; share(a, va, c, vc, sh, esh);
                P("      side bands %-20s: share %.3f +- %.3f; class contents", df.name, sh, esh);
                for (int k = 0; k < nMee; ++k) P(" %.0f", r[k].S);
                P(" (sum %.0f)\n", a + c);
            }
        }
        P("\n");
    }

    // ------------------------------------------------------------------ E
    {
        std::vector<Sel> e2 = {SelPEAK().f, [](const Row& r) { return Pre(r) && r.dz_pi < 40.f && !(r.mm4 < 2.85f); }};
        std::vector<Spectra> es = scan::FillMany(d, e2);
        const PeakFit p = scan::FitPeakMain(es[0]), f = scan::FitPeakMain(es[1]);
        const double eps = p.S / (p.S + f.S);
        const double err = std::sqrt(f.S * f.S * p.eS * p.eS + p.S * p.S * f.eS * f.eS) / std::pow(p.S + f.S, 2);
        P("E. EFFICIENCY OF MM4 < 2.85 GIVEN dz_pi < 40 mm (pass/fail P_main fits of the DIRECT rows): pass %.0f +- %.0f, fail %.0f +- %.0f, eps = %.3f +- %.3f\n\n",
          p.S, p.eS, f.S, f.eS, eps, err);
    }

    // ------------------------------------------------------------------ figure: S against the lower fit bound
    {
        TCanvas* c = new TCanvas("cfs", "", 1500, 620);
        c->Divide(2, 1, 0.004, 0.004);
        const int pick[2] = {5, 0};   // PEAK, pre (isBest >= 0)
        const char* lab[2] = {"(a) PEAK", "(b) pre"};
        for (int q = 0; q < 2; ++q) {
            c->cd(q + 1);
            gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.15); gPad->SetTopMargin(0.10); gPad->SetRightMargin(0.03);
            const Model& m = M[pick[q]];
            double ymin = m.nom.S, ymax = m.nom.S;
            for (double s : m.oS) { ymin = std::min(ymin, s); ymax = std::max(ymax, s); }
            for (double s : m.gS) { ymin = std::min(ymin, s); ymax = std::max(ymax, s); }
            TH1D* fr = new TH1D(Form("fr%d", q), "", 10, 0.5025, 0.5300);
            fr->SetDirectory(nullptr);
            fr->SetMinimum(ymin - 0.10 * (ymax - ymin) - 100); fr->SetMaximum(ymax + 0.85 * (ymax - ymin) + 100);
            fr->GetXaxis()->SetTitle("lower bound of the fit range [GeV/c^{2}]");
            fr->GetYaxis()->SetTitle("#eta-peak content S [rows]");
            fr->GetXaxis()->SetNdivisions(506);
            fr->GetXaxis()->SetTitleSize(0.058); fr->GetYaxis()->SetTitleSize(0.058);
            fr->GetXaxis()->SetLabelSize(0.052); fr->GetYaxis()->SetLabelSize(0.052);
            fr->GetYaxis()->SetTitleOffset(1.30);
            fr->Draw("axis");
            TBox* bx = new TBox(0.5025, m.nom.S - m.nom.eS, 0.5300, m.nom.S + m.nom.eS);
            bx->SetFillColorAlpha(kMagenta, 0.12); bx->Draw();
            TLine* le = new TLine(0.5202, fr->GetMinimum(), 0.5202, ymax + 0.08 * (ymax - ymin));
            le->SetLineColor(kGreen + 2); le->SetLineStyle(2); le->SetLineWidth(2); le->Draw();
            TLine* ln = new TLine(0.5025, m.nom.S, 0.5300, m.nom.S);
            ln->SetLineColor(kMagenta + 1); ln->SetLineWidth(2); ln->Draw();
            TGraph* go = new TGraph();
            for (size_t k = 0; k < m.oS.size(); ++k) go->SetPoint(go->GetN(), m.oLo[k] + 0.00025 * ((int)(k % 5) - 2) * 0.4, m.oS[k]);
            go->SetMarkerStyle(24); go->SetMarkerColor(kGray + 2); go->SetMarkerSize(1.1); go->Draw("P");
            TGraph* gn = new TGraph();
            for (size_t k = 0; k < m.gS.size(); ++k) gn->SetPoint(gn->GetN(), m.gLo[k] + 0.00025 * ((int)(k % 5) - 2) * 0.4, m.gS[k]);
            gn->SetMarkerStyle(20); gn->SetMarkerColor(kBlue + 1); gn->SetMarkerSize(0.9); gn->Draw("P");
            Txt(0.19, 0.82, lab[q], 0.062);
            Txt(0.19, 0.75, Form("S = %.0f #pm %.0f (stat)", m.nom.S, m.nom.eS), 0.050, kMagenta + 1);
            Txt(0.19, 0.69, Form("old envelope #pm%.0f (rms %.0f)", m.old.half, m.old.rms), 0.050, kGray + 2);
            Txt(0.19, 0.63, Form("edge-safe grid rms %.0f", m.gRms), 0.050, kBlue + 1);
            TLegend* lg = new TLegend(0.60, 0.66, 0.96, 0.88);
            lg->SetTextSize(0.046); lg->SetFillStyle(0);
            lg->AddEntry(go, "old, p > 0.05", "p");
            lg->AddEntry(gn, "edge-safe, p > 0.05", "p");
            lg->AddEntry(le, "bump edge 0.5202", "l");
            lg->Draw();
        }
        Save(c, "fix_spread");
    }
    CloseTxt();
    printf("fix_checks: done in %.0f s\n", sw.RealTime());
}
