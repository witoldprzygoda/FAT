// spec_m4_cnt.C -- COPY of plots/exp_report/spec_m4.C made by figure agent A (integ4, 2026-11-15) with
// the figure strings 'rows' -> 'counts' (analyst's terminology) and the outputs of the copy redirected:
// figures of the scope -> docs/latex/figs_exp/<base>_cnt.{pdf,png}; everything else -> /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figA/work/spec
// Built by /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figA/work/mk_cnt.py; the physics code is unchanged.
// ============================================================================
// spec_m4.C  --  M(pi+pi-e+e-) spectra, like-sign composition, vertex-cut effect
//                and counts in the standard M4 regions (current production)
// ============================================================================
// HADES p+p, T = 4.53 GeV (feb22), channel pi+ pi- e+ e-. EXPERIMENTAL DATA ONLY.
// Production pippimepem_isBestAll_T4530 (FAT: event vertex z > -500 mm, PT3
// trigbit == 8192, all stored combinations), corrected momenta (cor_ok == 1).
//
// What it does
//   1. fills M4 = M(pi+pi-e+e-) of the unlike-sign rows (N) and of the two
//      like-sign hypotheses (N++, N--) in 2.5 MeV bins (0.25-3.25 GeV/c2), for
//      8 base selections x {isBest >= 0, isBest == 1} x {no OA cut, OA > 2 deg};
//      CB = 2 sqrt(N++ N--) in 10 MeV bins, linearly interpolated (ScanEngine,
//      k = 1); N - CB;
//   2. sums N, N++, N--, CB, N - CB in the standard M4 regions with the
//      event-group errors of ScanEngine (variance = sum over groups of adjacent
//      rows with equal (file_id, run, event) of k^2; CB error propagated per
//      group); visible yields (N - CB) x 64/63 / L, L = 6079.2 nb^-1;
//   3. reference eta-peak fit P_main (frozen Gaussian 0.5439 GeV/c2 / 5.4 MeV +
//      pol2 on N - CB, 0.5225-0.6100) and its model spread (48 variations,
//      p > 0.05) per selection, only to annotate the figures;
//   4. like-sign composition (N++, N--, N++/N--, CB/N versus M4) and the
//      effect of the vertex cuts (vz > -160 mm, both pions in the vertex fit,
//      dz_pi < 40 mm) on unlike-sign versus like-sign rows, as row pass
//      fractions versus M4 with the event-group (ratio-estimator) errors.
//
// Selections (see spec_common.h; all rows of the FAT production, which already
// has vertex z > -500 mm and PT3):
//   pre = cor_ok && isBest >= 0 && vz > -160 mm   (= frozen BUMP)
//   pre + vertex = pre + both pions in the event-vertex fit
//   peak nominal (+ vertex) = pre + oa_lab < 55 + M(pipi) > 0.30 (+ vertex)
//   bump nominal (+ vertex) = pre + oa_lab < 40 + 0.30 < M(pipi) < 0.41 (+ vertex)
//   frozen PEAK = pre + dz_pi < 40 mm + MM4 < 2.85; frozen DIRECT = pre + dz_pi < 40 mm
//   variants: isBest == 1, OA(e+e-) > 2 deg
// Regions [GeV/c2]: < 0.455, bump window 0.4550-0.5225, gap 0.5225-0.53364,
//   peak window 0.53364-0.55416, upper gap 0.55416-0.56, 0.56-0.70, 0.70-0.90,
//   0.90-1.10, > 1.10, total; and the simulation-centred windows of the frozen
//   selections 0.4576-0.5255, 0.5333-0.5587 (reference only).
//
// Inputs : /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/scan_v2/scan_compact.root
//          (scan::kCompact, written by plots/scan_prepare.C) through ../ScanEngine.h, spec_common.h
// Outputs: /home/przygoda/HADES/FAT/FAT/docs/latex/figs_exp/spec_m4_overview, spec_m4_zoom,
//          spec_m4_isbest, spec_likesign, spec_vertexcut  (.pdf vector + .png)
//          /home/przygoda/HADES/FAT/FAT/docs/latex/tables_exp/spec_counts.txt, spec_likesign.txt,
//          spec_counts_pre.tex, spec_counts_sel.tex, spec_counts_sel_oa2.tex, spec_yields.tex, spec_vertexcut.tex
//          /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/spec/spec_m4.root
//
// Run (from the repository root /home/przygoda/HADES/FAT/FAT; about 1-2 min, 4 GB RAM):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/aclic_spec", kTRUE)' 'plots/exp_report/spec_m4.C+'
// ============================================================================
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TSystem.h"

#include "spec_common_cnt.h"

namespace specm4 {

using namespace spec;
using scan::PeakFit;
using scan::Spectra;

// ---- regions --------------------------------------------------------------------
struct Region { const char* tag; const char* name; const char* tex; double lo, hi; };
const int kNReg = 12;
const Region kReg[kNReg] = {
    {"below",   "M4 < 0.455",            "$<0.455$",               0.0,    0.455},
    {"bump",    "bump window",           "bump window 0.455--0.5225", 0.455, 0.5225},
    {"gap",     "gap",                   "gap 0.5225--0.53364",    0.5225, scan::kPeakLo},
    {"peak",    "peak window",           "peak window 0.53364--0.55416", scan::kPeakLo, scan::kPeakHi},
    {"ugap",    "upper gap",             "0.55416--0.56",          scan::kPeakHi, 0.56},
    {"r56_70",  "0.56-0.70",             "0.56--0.70",             0.56,   0.70},
    {"r70_90",  "0.70-0.90",             "0.70--0.90",             0.70,   0.90},
    {"r90_110", "0.90-1.10",             "0.90--1.10",             0.90,   1.10},
    {"above",   "M4 > 1.10",             "$>1.10$",                1.10,   3.25},
    {"total",   "total",                 "total",                  0.0,    3.25},
    {"simbump", "sim bump win 0.4576-0.5255", "sim. bump window 0.4576--0.5255", kSimBumpLo, kSimBumpHi},
    {"simpeak", "sim peak win 0.5333-0.5587", "sim. peak window 0.5333--0.5587", kSimPeakLo, kSimPeakHi},
};
const int kNStd = 9;   // regions 0..8 partition the M4 axis

// ---- variants -----------------------------------------------------------------------
const int kNVar = 4;
const char* const kVarTag[kNVar]  = {"best>=0", "best>=0,OA>2", "best==1", "best==1,OA>2"};
const char* const kVarName[kNVar] = {"isBest >= 0", "isBest >= 0, OA > 2 deg", "isBest == 1", "isBest == 1, OA > 2 deg"};

// ---- text output ------------------------------------------------------------------
FILE* gOut = nullptr;
void P(const char* fmt, ...) {
    va_list a; va_start(a, fmt); vprintf(fmt, a); va_end(a);
    if (gOut) { va_start(a, fmt); vfprintf(gOut, fmt, a); va_end(a); }
}

// ---- row pass fractions with the event-group (ratio estimator) error -----------------
// f = sum kp / sum kt; var(f) = sum_g (kp_g - f kt_g)^2 / (sum kt)^2
struct FracAcc {
    std::vector<double> p, t, pp, tt, pt;
    void Init(int n) { p.assign(n, 0); t.assign(n, 0); pp.assign(n, 0); tt.assign(n, 0); pt.assign(n, 0); }
    double F(int b) const { return t[b] > 0 ? p[b] / t[b] : 0; }
    double E(int b) const {
        if (t[b] <= 0) return 0;
        const double f = F(b);
        return std::sqrt(std::max(pp[b] - 2 * f * pt[b] + f * f * tt[b], 0.)) / t[b];
    }
};
typedef std::function<bool(const Row&)> Pred;
typedef std::function<int(const Row&)> BinOf;
void FracFill(const std::vector<Row>& v, const Pred& base, const Pred& cut, const BinOf& binOf, FracAcc& a) {
    int cur = -1;
    std::vector<std::array<int, 3>> g;   // bin, kp, kt
    auto flush = [&]() {
        for (const auto& e : g) {
            a.p[e[0]] += e[1]; a.t[e[0]] += e[2];
            a.pp[e[0]] += (double)e[1] * e[1]; a.tt[e[0]] += (double)e[2] * e[2]; a.pt[e[0]] += (double)e[1] * e[2];
        }
        g.clear();
    };
    for (const Row& r : v) {
        if (r.grp != cur) { flush(); cur = r.grp; }
        if (!base(r)) continue;
        const int b = binOf(r);
        if (b < 0) continue;
        const int c = cut(r) ? 1 : 0;
        bool found = false;
        for (auto& e : g) if (e[0] == b) { e[1] += c; e[2] += 1; found = true; break; }
        if (!found) g.push_back({b, c, 1});
    }
    flush();
}

// ---- helpers for the figures --------------------------------------------------------
TGraphErrors* Gr(int col, int mk, double ms = 0.9) {
    TGraphErrors* g = new TGraphErrors();
    g->SetLineColor(col); g->SetMarkerColor(col); g->SetMarkerStyle(mk); g->SetMarkerSize(ms);
    return g;
}
TGraph* Curve(const PeakFit& f, bool total, double x1, double x2, int col, int sty, double yscale = 1.0) {
    TGraph* g = new TGraph();
    const int n = 300;
    for (int i = 0; i <= n; ++i) {
        const double m = x1 + (x2 - x1) * i / n;
        g->SetPoint(i, m, yscale * (total ? f.Total(m) : f.Bkg(m)));
    }
    g->SetLineColor(col); g->SetLineStyle(sty); g->SetLineWidth(2);
    return g;
}
TH1D* Frame(const char* name, double x1, double x2, double y1, double y2, const char* xt, const char* yt) {
    TH1D* h = new TH1D(name, "", 100, x1, x2);
    h->SetDirectory(nullptr);
    h->SetMinimum(y1); h->SetMaximum(y2);
    h->GetXaxis()->SetTitle(Fx(xt).c_str()); h->GetYaxis()->SetTitle(Fx(yt).c_str());
    h->GetXaxis()->SetNdivisions(506);
    h->Draw("axis");
    return h;
}
// histogram scaled by f (content and error), for axes in units of 10^3 rows
TH1D* Scaled(TH1D* h, double f) {
    for (int i = 0; i <= h->GetNbinsX() + 1; ++i) {
        h->SetBinContent(i, f * h->GetBinContent(i));
        h->SetBinError(i, f * h->GetBinError(i));
    }
    return h;
}
// copy of a spectrum with the CB scaled by f (for the +-3 % CB sensitivity of P_main)
Spectra ScaleCB(const Spectra& in, double f) {
    Spectra o = in;
    for (auto& x : o.cbW) x *= f;
    for (auto& x : o.cbWVar) x *= f * f;
    for (int i = 0; i < o.nbin; ++i) {
        o.cb[i] *= f; o.cbVar[i] *= f * f;
        o.sig[i] = o.n[0][i] - o.cb[i]; o.sigVar[i] = o.n[0][i] + o.cbVar[i];
    }
    for (auto& w : o.win) { w.cb *= f; w.cbVar *= f * f; }
    return o;
}
std::string Pm(double v, double e, int prec = 0) {
    return std::string(Form("%.*f +- %.*f", prec, v, prec, e));
}
std::string TexPm(double v, double e, int prec = 0) {
    return std::string(Form("$%.*f\\pm%.*f$", prec, v, prec, e));
}

} // namespace specm4

void spec_m4_cnt() {
    using namespace specm4;
    TStopwatch sw;
    Style();
    gSystem->mkdir(kWork, kTRUE); gSystem->mkdir(kCntScratch, kTRUE);
    gSystem->mkdir(kFigOut, kTRUE);
    gSystem->mkdir(kTabOut, kTRUE);

    scan::Data d;
    if (!scan::Load(d)) return;

    // ---- selections: 8 base x 4 variants -----------------------------------------
    const std::vector<NSel> base = BaseSelections();
    const int nb = (int)base.size();
    std::vector<Sel> sels;
    for (int b = 0; b < nb; ++b) {
        sels.push_back(base[b].sel);
        sels.push_back(WithOA(base[b].sel));
        sels.push_back(scan::SelBest(base[b].sel));
        sels.push_back(scan::SelBest(WithOA(base[b].sel)));
    }
    auto idx = [](int b, int v) { return 4 * b + v; };

    scan::FillOpt opt;
    opt.lo = 0.25; opt.hi = 3.25; opt.nbin = 1200; opt.reb = 4;     // 2.5 MeV, CB in 10 MeV (same grid as ScanEngine)
    for (int r = 0; r < kNReg; ++r) opt.windows.push_back(std::make_pair(kReg[r].lo, kReg[r].hi));
    TStopwatch sf;
    const std::vector<Spectra> sp = scan::FillMany(d, sels, opt);
    printf("FillMany of %zu selections: %.1f s\n", sels.size(), sf.RealTime());

    // ---- P_main + model spread -----------------------------------------------------
    std::vector<PeakFit> fit(sels.size());
    std::vector<scan::Spread> spr(sels.size());
    std::vector<double> sUp(sels.size()), sDn(sels.size());
    for (size_t s = 0; s < sels.size(); ++s) {
        fit[s] = scan::FitPeakMain(sp[s]);
        spr[s] = scan::PeakModelSpread(sp[s]);
        sUp[s] = scan::FitPeakMain(ScaleCB(sp[s], 1.0 + kCBsyst)).S;
        sDn[s] = scan::FitPeakMain(ScaleCB(sp[s], 1.0 - kCBsyst)).S;
    }

    // ---- rows of the no_hf runs (L = 0; the luminosity rule drops them) -----------
    const int kNoHf[4] = {1244739, 2109357, 2109532, 3274681};
    double nNoHf[3] = {0, 0, 0}, nPre[3] = {0, 0, 0};
    for (int k = 0; k < 3; ++k)
        for (const Row& r : d.s[k]) {
            if (!Pre(r)) continue;
            nPre[k] += 1;
            for (int q = 0; q < 4; ++q) if (r.run == kNoHf[q]) nNoHf[k] += 1;
        }

    // ================================================================================
    // A. tables
    // ================================================================================
    gOut = fopen(Form("%s/spec_counts.txt", kTabOut), "w");
    P("# spec_counts.txt -- counts of the current production in the standard M(pi+pi-e+e-) regions\n");
    P("# written by plots/exp_report/spec_m4.C ; data: %s\n", scan::kCompact);
    P("# HADES p+p 4.53 GeV feb22, PT3 (trigbit 8192), event vertex z > -500 mm (FAT), corrected momenta (cor_ok == 1)\n");
    P("# Selections: pre = cor_ok && isBest >= 0 && vz > -160 mm (= frozen BUMP); +vertex = both pions in the vertex fit;\n");
    P("#   peak nominal = pre + oa_lab < 55 + M(pipi) > 0.30; bump nominal = pre + oa_lab < 40 + 0.30 < M(pipi) < 0.41;\n");
    P("#   frozen PEAK = pre + dz_pi < 40 mm + MM4 < 2.85 GeV/c2; frozen DIRECT = pre + dz_pi < 40 mm.\n");
    P("# Variants: isBest >= 0 (main) / isBest == 1 (comparison); without / with OA(e+e-) > 2 deg (OA > 2 is not a fixed cut).\n");
    P("# N = unlike-sign rows; N++, N-- = like-sign rows; CB = 2 sqrt(N++ N--) in 10 MeV bins interpolated to 2.5 MeV, k = 1.\n");
    P("# Errors: event-group errors (groups = adjacent rows with equal (file_id, run, event); var = sum k^2), CB propagated per group.\n");
    P("# dCB3 = 3 %% of CB (sensitivity of N - CB to a +-3 %% CB change). Y = (N - CB) x 64/63 / L, L = %.1f nb^-1 (runs present;\n", kLumi);
    P("#   all runs %.1f); Y in rows per nb^-1 (= nb), stat. error only; L has %.1f %% total uncertainty; PT3 live fraction: one-sided up to -6.8 %%.\n", kLumiAll, 100 * kLumiRel);
    P("# Not corrected for acceptance or efficiency. Regions in GeV/c2. The data windows are used; the sim-centred windows only for reference.\n");
    P("# Rows of the four no_hf runs (L = 0) inside pre: us %.0f of %.0f (%.4f %%), pp %.0f, mm %.0f -- kept (negligible).\n",
      nNoHf[0], nPre[0], 100 * nNoHf[0] / std::max(nPre[0], 1.), nNoHf[1], nNoHf[2]);
    P("#\n");
    for (int b = 0; b < nb; ++b) {
        for (int v = 0; v < kNVar; ++v) {
            const int s = idx(b, v);
            const Spectra& q = sp[s];
            P("== %s | %s | selected rows (cor_ok, any M4): us %.0f pp %.0f mm %.0f\n", base[b].name.c_str(), kVarName[v],
              q.nsel[0], q.nsel[1], q.nsel[2]);
            P("   %-27s %17s %8s %8s %17s %6s %19s %7s %15s %9s\n", "region", "N", "N++", "N--", "CB", "CB/N", "N - CB", "dCB3",
              "Y=(N-CB)/L [nb]", "N/L [nb]");
            for (int r = 0; r < kNReg; ++r) {
                const scan::WinSum& w = q.win[2 + r];
                P("   %-27s %17s %8.0f %8.0f %17s %6.3f %19s %7.0f %15s %9.3f\n", kReg[r].name, Pm(w.n[0], std::sqrt(w.var[0])).c_str(),
                  w.n[1], w.n[2], Pm(w.cb, std::sqrt(w.cbVar)).c_str(), w.cbShare(), Pm(w.sig(), w.sigErr()).c_str(), kCBsyst * w.cb,
                  Pm(w.sig() * kPT3 / kLumi, w.sigErr() * kPT3 / kLumi, 3).c_str(), w.n[0] * kPT3 / kLumi);
            }
            const PeakFit& f = fit[s];
            P("   P_main (frozen 0.5439/5.4 MeV, pol2, 0.5225-0.6100): S = %.0f +- %.0f (stat) +- %.0f (model, %d of %d fits p > 0.05), "
              "F = %.2f, chi2/ndf = %.1f/%d (p %.3f), S/B = %.4f; S x 64/63 / L = %.3f +- %.3f (stat) +- %.3f (model) nb;\n"
              "   S with CB x 1.03 / x 0.97: %.0f / %.0f (change %+.0f / %+.0f)\n\n",
              f.S, f.eS, spr[s].half, spr[s].npass, spr[s].nfit, f.F, f.chi2, f.ndf, f.prob, f.SoB, f.S * kPT3 / kLumi, f.eS * kPT3 / kLumi,
              spr[s].half * kPT3 / kLumi, sUp[s], sDn[s], sUp[s] - f.S, sDn[s] - f.S);
        }
    }
    // consistency: the nine standard regions add up to the total
    {
        P("# closure: sum of the 9 standard regions - total (pre, isBest >= 0): N %.3f, CB %.3f\n",
          [&]() { double s = 0; for (int r = 0; r < kNStd; ++r) s += sp[0].win[2 + r].n[0]; return s - sp[0].win[2 + 9].n[0]; }(),
          [&]() { double s = 0; for (int r = 0; r < kNStd; ++r) s += sp[0].win[2 + r].cb; return s - sp[0].win[2 + 9].cb; }());
        P("# check: win[0] (ScanEngine peak window) = region 'peak': N %.0f / %.0f, CB %.2f / %.2f\n", sp[0].win[0].n[0],
          sp[0].win[2 + 3].n[0], sp[0].win[0].cb, sp[0].win[2 + 3].cb);
    }
    fclose(gOut); gOut = nullptr;

    // ---- LaTeX: pre in detail ------------------------------------------------------------
    {
        FILE* o = fopen(Form("%s/spec_counts_pre.tex", kTabOut), "w");
        fprintf(o, "%% written by plots/exp_report/spec_m4.C -- selection pre, counts in the standard M4 regions\n");
        fprintf(o, "\\begin{table}[H]\n\\centering\n\\scriptsize\n");
        fprintf(o, "\\caption{Rows of the selection \\emph{pre} (corrected momenta, $\\mathrm{cor\\_ok}=1$, event vertex $z>-160$~mm, "
                   "PT3, all stored combinations $\\mathrm{isBest}\\ge0$) in the standard regions of $M(\\pi^+\\pi^-e^+e^-)$. "
                   "$N$: unlike-sign rows; CB $=2\\sqrt{N_{++}N_{--}}$ (10~MeV/$c^2$ bins interpolated, $k=1$). Errors: statistical, "
                   "with the event-group correlation. $\\Delta_{3\\%%}$: change of $N-\\mathrm{CB}$ for a 3~\\%% change of CB. "
                   "$Y=(N-\\mathrm{CB})\\times\\frac{64}{63}/L$ with $L=%.1f$~nb$^{-1}$ (rows per nb$^{-1}$; statistical error only; "
                   "$L$ has a 6.9~\\%% uncertainty). Last two columns: $N-\\mathrm{CB}$ with $\\mathrm{OA}(e^+e^-)>2^\\circ$ and with "
                   "$\\mathrm{isBest}=1$. Not corrected for acceptance or efficiency.}\n", kLumi);
        fprintf(o, "\\label{tab:spec_counts_pre}\n");
        fprintf(o, "\\resizebox{\\textwidth}{!}{\\begin{tabular}{lrrrrrrrr}\n\\toprule\n");
        fprintf(o, "region [GeV/$c^2$] & $N$ & CB & CB/$N$ & $N-\\mathrm{CB}$ & $\\Delta_{3\\%%}$ & $Y$ [nb] & OA$>2^\\circ$: $N-\\mathrm{CB}$ & "
                   "isBest$=1$: $N-\\mathrm{CB}$ \\\\\n\\midrule\n");
        for (int r = 0; r < kNReg; ++r) {
            if (r == kNStd || r == kNStd + 1) fprintf(o, "\\midrule\n");
            const scan::WinSum& w = sp[idx(0, 0)].win[2 + r];
            const scan::WinSum& wo = sp[idx(0, 1)].win[2 + r];
            const scan::WinSum& w1 = sp[idx(0, 2)].win[2 + r];
            fprintf(o, "%s & %s & %s & %.3f & %s & %.0f & %s & %s & %s \\\\\n", kReg[r].tex, TexPm(w.n[0], std::sqrt(w.var[0])).c_str(),
                    TexPm(w.cb, std::sqrt(w.cbVar)).c_str(), w.cbShare(), TexPm(w.sig(), w.sigErr()).c_str(), kCBsyst * w.cb,
                    TexPm(w.sig() * kPT3 / kLumi, w.sigErr() * kPT3 / kLumi, 2).c_str(), TexPm(wo.sig(), wo.sigErr()).c_str(),
                    TexPm(w1.sig(), w1.sigErr()).c_str());
        }
        fprintf(o, "\\bottomrule\n\\end{tabular}}\n\\end{table}\n");
        fclose(o);
    }
    // ---- LaTeX: N - CB of all base selections (isBest >= 0, without and with OA > 2) ----
    for (int oa = 0; oa < 2; ++oa) {
        FILE* o = fopen(Form("%s/%s", kTabOut, oa == 0 ? "spec_counts_sel.tex" : "spec_counts_sel_oa2.tex"), "w");
        fprintf(o, "%% written by plots/exp_report/spec_m4.C -- N - CB per selection and region, isBest >= 0%s\n", oa ? ", OA > 2 deg" : "");
        fprintf(o, "\\begin{table}[H]\n\\centering\n\\scriptsize\n\\setlength{\\tabcolsep}{3pt}\n");
        fprintf(o, "\\caption{$N-\\mathrm{CB}$ (unlike-sign rows minus CB $=2\\sqrt{N_{++}N_{--}}$, $k=1$) in the standard regions of "
                   "$M(\\pi^+\\pi^-e^+e^-)$ for the selections of the report, all stored combinations ($\\mathrm{isBest}\\ge0$)%s. "
                   "Statistical errors with the event-group correlation. The last row is the CB share CB/$N$ in the peak window. "
                   "The frozen PEAK and DIRECT selections were defined on simulation; here they are evaluated in the data windows.}\n",
                oa ? ", with $\\mathrm{OA}(e^+e^-)>2^\\circ$" : ", no OA cut");
        fprintf(o, "\\label{tab:spec_counts_sel%s}\n", oa ? "_oa2" : "");
        fprintf(o, "\\resizebox{\\textwidth}{!}{\\begin{tabular}{l");
        for (int b = 0; b < nb; ++b) fprintf(o, "r");
        fprintf(o, "}\n\\toprule\nregion [GeV/$c^2$]");
        const char* tn[8] = {"pre", "pre+vtx", "peak nom.", "peak nom.+vtx", "bump nom.", "bump nom.+vtx", "PEAK", "DIRECT"};
        for (int b = 0; b < nb; ++b) fprintf(o, " & %s", tn[b]);
        fprintf(o, " \\\\\n\\midrule\n");
        for (int r = 0; r < kNReg; ++r) {
            if (r == kNStd || r == kNStd + 1) fprintf(o, "\\midrule\n");
            fprintf(o, "%s", kReg[r].tex);
            for (int b = 0; b < nb; ++b) {
                const scan::WinSum& w = sp[idx(b, oa)].win[2 + r];
                fprintf(o, " & %s", TexPm(w.sig(), w.sigErr()).c_str());
            }
            fprintf(o, " \\\\\n");
        }
        fprintf(o, "\\midrule\nCB/$N$ in the peak window");
        for (int b = 0; b < nb; ++b) fprintf(o, " & %.3f", sp[idx(b, oa)].win[2 + 3].cbShare());
        fprintf(o, " \\\\\n\\bottomrule\n\\end{tabular}}\n\\end{table}\n");
        fclose(o);
    }
    // ---- LaTeX: visible yields per nb^-1 -------------------------------------------------
    {
        FILE* o = fopen(Form("%s/spec_yields.tex", kTabOut), "w");
        fprintf(o, "%% written by plots/exp_report/spec_m4.C -- visible yields (N - CB) x 64/63 / L\n");
        fprintf(o, "\\begin{table}[H]\n\\centering\n\\scriptsize\n\\setlength{\\tabcolsep}{3pt}\n");
        fprintf(o, "\\caption{Visible yields $Y=(N-\\mathrm{CB})\\times\\frac{64}{63}/L$ in rows per nb$^{-1}$ (i.e.\\ in nb), "
                   "$L=%.1f$~nb$^{-1}$ (START-free luminosity of the analysed runs), all stored combinations ($\\mathrm{isBest}\\ge0$); "
                   "first line of each region without, second line with $\\mathrm{OA}(e^+e^-)>2^\\circ$. Statistical errors only; "
                   "$L$ carries a 6.9~\\%% uncertainty and a one-sided $-6.8$~\\%% PT3 live-fraction uncertainty. $Y$ is not corrected "
                   "for acceptance or efficiency and contains the physical background under the peak and the bump.}\n", kLumi);
        fprintf(o, "\\label{tab:spec_yields}\n");
        fprintf(o, "\\resizebox{\\textwidth}{!}{\\begin{tabular}{ll");
        for (int b = 0; b < nb; ++b) fprintf(o, "r");
        fprintf(o, "}\n\\toprule\nregion [GeV/$c^2$] & OA cut");
        const char* tn[8] = {"pre", "pre+vtx", "peak nom.", "peak nom.+vtx", "bump nom.", "bump nom.+vtx", "PEAK", "DIRECT"};
        for (int b = 0; b < nb; ++b) fprintf(o, " & %s", tn[b]);
        fprintf(o, " \\\\\n\\midrule\n");
        for (int r = 0; r < kNStd + 1; ++r) {
            for (int oa = 0; oa < 2; ++oa) {
                fprintf(o, "%s & %s", oa == 0 ? kReg[r].tex : "", oa == 0 ? "none" : "$>2^\\circ$");
                for (int b = 0; b < nb; ++b) {
                    const scan::WinSum& w = sp[idx(b, oa)].win[2 + r];
                    fprintf(o, " & %s", TexPm(w.sig() * kPT3 / kLumi, w.sigErr() * kPT3 / kLumi, 2).c_str());
                }
                fprintf(o, " \\\\\n");
            }
        }
        fprintf(o, "\\bottomrule\n\\end{tabular}}\n\\end{table}\n");
        fclose(o);
    }

    // ================================================================================
    // B. like-sign composition and vertex-cut pass fractions
    // ================================================================================
    // pass fractions versus M4 (10 MeV, 0.35-1.10) and in the standard regions
    const double fLo = 0.35, fHi = 1.10; const int fNb = 75;
    BinOf binM4 = [=](const Row& r) { const double m = r.m4; return (m >= fLo && m < fHi) ? (int)((m - fLo) / (fHi - fLo) * fNb) : -1; };
    BinOf binReg = [](const Row& r) {
        const double m = r.m4;
        if (!(m == m)) return -1;
        for (int q = 0; q < kNStd; ++q) if (m >= kReg[q].lo && m < kReg[q].hi) return q;
        return -1;
    };
    const int nCut = 3;
    const char* cutName[nCut] = {"vz > -160 mm (base: cor_ok, isBest >= 0)", "both pions in vertex fit (base: pre)", "dz_pi < 40 mm (base: pre)"};
    const char* cutRoot[nCut] = {"v_{z} > #minus160 mm", "both #pi in the vertex fit", "#Deltaz_{#pi} < 40 mm"};
    Pred cutBase[nCut] = {[](const Row& r) { return PreNoVz(r); }, [](const Row& r) { return Pre(r); }, [](const Row& r) { return Pre(r); }};
    Pred cutPass[nCut] = {[](const Row& r) { return r.vz > kVz; }, [](const Row& r) { return PiVtx(r); }, [](const Row& r) { return r.dz_pi < kDzPi; }};
    FracAcc fM4[nCut][3], fRg[nCut][3];
    for (int c = 0; c < nCut; ++c)
        for (int k = 0; k < 3; ++k) {
            fM4[c][k].Init(fNb); fRg[c][k].Init(kNStd);
            FracFill(d.s[k], cutBase[c], cutPass[c], binM4, fM4[c][k]);
            FracFill(d.s[k], cutBase[c], cutPass[c], binReg, fRg[c][k]);
        }
    // the same for isBest == 1 (one row per event; vz > -160 only)
    FracAcc fRgBest[3];
    for (int k = 0; k < 3; ++k) {
        fRgBest[k].Init(kNStd);
        FracFill(d.s[k], [](const Row& r) { return PreNoVz(r) && r.isBest == 1; }, [](const Row& r) { return r.vz > kVz; }, binReg, fRgBest[k]);
    }
    // vz shapes of the three samples (cor_ok, isBest >= 0, no vz > -160)
    TH1D* hVz[3];
    const char* smpRoot[3] = {"e^{+}e^{-} (unlike-sign)", "e^{+}e^{+}", "e^{-}e^{-}"};
    for (int k = 0; k < 3; ++k) {
        hVz[k] = new TH1D(Form("hvz_%s", scan::SampleName(k)), "", 120, -500, 100);
        hVz[k]->SetDirectory(nullptr);
        for (const Row& r : d.s[k]) if (PreNoVz(r)) hVz[k]->Fill(r.vz);
    }

    gOut = fopen(Form("%s/spec_likesign.txt", kTabOut), "w");
    P("# spec_likesign.txt -- like-sign composition and the effect of the vertex cuts on like-sign vs unlike-sign rows\n");
    P("# written by plots/exp_report/spec_m4.C ; data: %s ; FAT: vertex z > -500 mm, PT3, corrected momenta (cor_ok == 1)\n", scan::kCompact);
    P("# A. N++, N--, N++/N--, CB/N per region and selection (isBest >= 0 and isBest == 1)\n");
    for (int b = 0; b < nb; ++b)
        for (int v = 0; v < kNVar; v += 2) {
            P("== %s | %s\n", base[b].name.c_str(), kVarName[v]);
            P("   %-27s %9s %8s %8s %16s %16s %7s\n", "region", "N", "N++", "N--", "N++/N--", "CB/N", "N--/N++");
            for (int r = 0; r < kNStd + 1; ++r) {
                const scan::WinSum& w = sp[idx(b, v)].win[2 + r];
                const double rt = w.n[2] > 0 ? w.n[1] / w.n[2] : 0;
                const double ert = (w.n[1] > 0 && w.n[2] > 0) ? rt * std::sqrt(w.var[1] / (w.n[1] * w.n[1]) + w.var[2] / (w.n[2] * w.n[2])) : 0;
                const double cs = w.cbShare();
                const double ecs = w.n[0] > 0 ? cs * std::sqrt(w.cbVar / std::max(w.cb * w.cb, 1e-9) + w.var[0] / (w.n[0] * w.n[0])) : 0;
                P("   %-27s %9.0f %8.0f %8.0f %16s %16s %7.2f\n", kReg[r].name, w.n[0], w.n[1], w.n[2], Pm(rt, ert, 4).c_str(),
                  Pm(cs, ecs, 4).c_str(), w.n[1] > 0 ? w.n[2] / w.n[1] : 0.);
            }
        }
    P("\n# B. Row pass fractions of the vertex cuts, unlike-sign (us) vs e+e+ (pp) vs e-e- (mm), per region; errors with the\n");
    P("#    event-group correlation (ratio estimator over groups). CBkept = sqrt(f++ f--) = kept fraction of 2 sqrt(N++ N--) of the\n");
    P("#    region (the M4-binned region CB after the cut is in section A); f(us)/CBkept = factor by which the cut lowers CB/N.\n");
    for (int c = 0; c < nCut; ++c) {
        P("== %s\n", cutName[c]);
        P("   %-27s %20s %20s %20s %8s %8s\n", "region", "f(us)", "f(pp)", "f(mm)", "CBkept", "f(us)/CBkept");
        for (int r = 0; r < kNStd; ++r) {
            const double cbk = std::sqrt(fRg[c][1].F(r) * fRg[c][2].F(r));
            P("   %-27s %20s %20s %20s %8.4f %8.3f\n", kReg[r].name, Pm(fRg[c][0].F(r), fRg[c][0].E(r), 4).c_str(),
              Pm(fRg[c][1].F(r), fRg[c][1].E(r), 4).c_str(), Pm(fRg[c][2].F(r), fRg[c][2].E(r), 4).c_str(), cbk,
              cbk > 0 ? fRg[c][0].F(r) / cbk : 0.);
        }
        // all regions together
        double p[3] = {0, 0, 0}, t[3] = {0, 0, 0};
        for (int k = 0; k < 3; ++k) for (int r = 0; r < kNStd; ++r) { p[k] += fRg[c][k].p[r]; t[k] += fRg[c][k].t[r]; }
        P("   %-27s %20.4f %20.4f %20.4f   (rows: us %.0f, pp %.0f, mm %.0f)\n", "all M4", p[0] / t[0], p[1] / t[1], p[2] / t[2], t[0], t[1], t[2]);
    }
    P("== vz > -160 mm, isBest == 1 rows only (base: cor_ok, isBest == 1)\n");
    for (int r = 0; r < kNStd; ++r)
        P("   %-27s %20s %20s %20s\n", kReg[r].name, Pm(fRgBest[0].F(r), fRgBest[0].E(r), 4).c_str(),
          Pm(fRgBest[1].F(r), fRgBest[1].E(r), 4).c_str(), Pm(fRgBest[2].F(r), fRgBest[2].E(r), 4).c_str());
    P("\n# C. The mandatory FAT cut vertex z > -500 mm is applied before the compact file; its effect on like-sign vs unlike-sign\n");
    P("#    rows is measured from the PAT input by plots/exp_report/spec_vz500.C (table spec_vz500.txt).\n");
    // denominators include the overflow (vz > 100 mm; verifier fix, the in-range integral made the fractions 0.1-0.4 % too high)
    P("\n# D. vz shapes (cor_ok, isBest >= 0, no vz > -160): fraction of rows with vz <= -160 mm: us %.4f, pp %.4f, mm %.4f\n",
      hVz[0]->Integral(1, hVz[0]->FindBin(-160.001)) / hVz[0]->Integral(0, hVz[0]->GetNbinsX() + 1),
      hVz[1]->Integral(1, hVz[1]->FindBin(-160.001)) / hVz[1]->Integral(0, hVz[1]->GetNbinsX() + 1),
      hVz[2]->Integral(1, hVz[2]->FindBin(-160.001)) / hVz[2]->Integral(0, hVz[2]->GetNbinsX() + 1));
    fclose(gOut); gOut = nullptr;
    // LaTeX: pass fractions of the three vertex cuts in the two windows and for all M4
    {
        FILE* o = fopen(Form("%s/spec_vertexcut.tex", kTabOut), "w");
        fprintf(o, "%% written by plots/exp_report/spec_m4.C -- row pass fractions of the vertex cuts, unlike-sign vs like-sign\n");
        fprintf(o, "\\begin{table}[H]\n\\centering\n\\scriptsize\n");
        fprintf(o, "\\caption{Row pass fractions of the vertex cuts for the unlike-sign rows ($e^+e^-$) and the two like-sign hypotheses "
                   "($e^+e^+$, $e^-e^-$), all stored combinations ($\\mathrm{isBest}\\ge0$), corrected momenta, PT3, event vertex "
                   "$z>-500$~mm (FAT). Base of $v_z>-160$~mm: $\\mathrm{cor\\_ok}=1$; base of the other two cuts: \\emph{pre}. "
                   "CB kept $=\\sqrt{f_{++}f_{--}}$ (kept part of $2\\sqrt{N_{++}N_{--}}$); the last column is the factor by which "
                   "the cut lowers CB/$N$. Errors: statistical, with the event-group correlation (ratio estimator).}\n");
        fprintf(o, "\\label{tab:spec_vertexcut}\n\\begin{tabular}{llrrrrr}\n\\toprule\n");
        fprintf(o, "cut & $M_4$ region & $f(e^+e^-)$ & $f(e^+e^+)$ & $f(e^-e^-)$ & CB kept & CB/$N$ factor \\\\\n\\midrule\n");
        const char* tc[nCut] = {"$v_z>-160$~mm", "both $\\pi$ in vertex fit", "$\\Delta z_\\pi<40$~mm"};
        for (int c = 0; c < nCut; ++c) {
            const int rr[2] = {1, 3};
            for (int q = 0; q < 3; ++q) {
                double f[3], e[3];
                if (q < 2) { for (int k = 0; k < 3; ++k) { f[k] = fRg[c][k].F(rr[q]); e[k] = fRg[c][k].E(rr[q]); } }
                else {
                    for (int k = 0; k < 3; ++k) {
                        double pp = 0, tt = 0; for (int r = 0; r < kNStd; ++r) { pp += fRg[c][k].p[r]; tt += fRg[c][k].t[r]; }
                        f[k] = pp / tt; e[k] = 0;
                    }
                }
                const double cbk = std::sqrt(f[1] * f[2]);
                if (q < 2)
                    fprintf(o, "%s & %s & $%.3f\\pm%.3f$ & $%.3f\\pm%.3f$ & $%.3f\\pm%.3f$ & %.3f & %.3f \\\\\n", q == 0 ? tc[c] : "",
                            q == 0 ? "bump window" : "peak window", f[0], e[0], f[1], e[1], f[2], e[2], cbk, cbk / f[0]);
                else
                    fprintf(o, " & all & $%.3f$ & $%.3f$ & $%.3f$ & %.3f & %.3f \\\\\n", f[0], f[1], f[2], cbk, cbk / f[0]);
            }
            if (c + 1 < nCut) fprintf(o, "\\midrule\n");
        }
        fprintf(o, "\\bottomrule\n\\end{tabular}\n\\end{table}\n");
        fclose(o);
    }

    // ================================================================================
    // C. figures
    // ================================================================================
    const char* xM4 = "M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]";
    const int cN = kBlack, cCB = kRed + 1, cS = kBlue + 1;
    const int colPeak = ColPeak(), colBump = ColBump();     // window shading (solid, light)
    // the six panels: pre (= BUMP), pre+vtx, peak nom.+vtx, bump nom.+vtx, PEAK, DIRECT
    const int pan[6] = {0, 1, 3, 5, 6, 7};
    const bool frozenPan[6] = {true, false, false, false, true, true};
    const char* panTitle[6] = {"pre (= frozen BUMP)", "pre + vertex", "peak nominal + vertex", "bump nominal + vertex",
                               "frozen PEAK", "frozen DIRECT"};
    TFile* fout = TFile::Open(Form("%s/spec_m4.root", kWork), "RECREATE");
    TBox* boxB = new TBox(); boxB->SetFillColor(colBump); boxB->SetLineColor(colBump);
    TBox* boxP = new TBox(); boxP->SetFillColor(colPeak); boxP->SetLineColor(colPeak);

    // ---- C1. log overview 0.35-1.10, 5 MeV bins, common y range ---------------------
    {
        TH1D *hn[6], *hc[6], *hs[6];
        double ymax = 0, ymin = 1e30;
        for (int p = 0; p < 6; ++p) {
            const Spectra& q = sp[idx(pan[p], 0)];
            hn[p] = Merge(q, "us", 2, Form("ov_us_%s", base[pan[p]].tag.c_str()));
            hc[p] = Merge(q, "cb", 2, Form("ov_cb_%s", base[pan[p]].tag.c_str()));
            hs[p] = Merge(q, "sig", 2, Form("ov_sig_%s", base[pan[p]].tag.c_str()));
            Deco(hn[p], cN, 20, 0.55); Deco(hc[p], cCB, 20, 0.55); Deco(hs[p], cS, 20, 0.55);
            for (int i = hn[p]->FindBin(0.3501); i <= hn[p]->FindBin(1.0999); ++i) {
                ymax = std::max(ymax, hn[p]->GetBinContent(i));
                if (hc[p]->GetBinContent(i) > 0) ymin = std::min(ymin, hc[p]->GetBinContent(i));
            }
            fout->cd(); hn[p]->Write(); hc[p]->Write(); hs[p]->Write();
        }
        const double y1 = std::max(0.5, 0.5 * ymin), y2 = 40 * ymax;
        printf("overview: common y range %.2f - %.0f\n", y1, y2);
        TCanvas* c = new TCanvas("c_ov", "", 2100, 1300);
        c->Divide(3, 2, 0.001, 0.001);
        for (int p = 0; p < 6; ++p) {
            c->cd(p + 1); gPad->SetLogy();
            Frame(Form("fr_ov_%d", p), 0.35, 1.10, y1, y2, xM4, "rows / 5 MeV/c^{2}");
            Band(scan::kBumpLo, scan::kBumpHi, y1, y2 / 30, colBump);
            Band(scan::kPeakLo, scan::kPeakHi, y1, y2 / 30, colPeak);
            hc[p]->Draw("e same"); hn[p]->Draw("e same"); hs[p]->Draw("e same");
            gPad->RedrawAxis();
            Txt(0.16, 0.935, Form("(%c) %s", 'a' + p, panTitle[p]), 0.054);
            const scan::WinSum& wp = sp[idx(pan[p], 0)].win[0];
            const scan::WinSum& wb = sp[idx(pan[p], 0)].win[1];
            Txt(0.95, 0.865, Form("CB/N: bump window %.3f, peak window %.3f", wb.cbShare(), wp.cbShare()), 0.042, 31);
            if (p == 0) {
                TLegend* lg = Leg(0.50, 0.17, 0.96, 0.40, 0.046);
                LE(lg, hn[p], "N (e^{+}e^{-} rows)", "pe");
                LE(lg, hc[p], "CB = 2#sqrt{N_{++}N_{#minus#scale[0.3]{ }#minus}}", "pe");
                LE(lg, hs[p], "N #minus CB", "pe");
                lg->Draw();
            }
            if (p == 1) {
                TLegend* lg = Leg(0.50, 0.20, 0.96, 0.36, 0.044);
                LE(lg, boxB, "bump window (data)", "f");
                LE(lg, boxP, "peak window (data)", "f");
                lg->Draw();
            }
        }
        Save(c, "spec_m4_overview");
    }

    // ---- C2. linear zoom 0.42-0.62 of N - CB, 2.5 MeV, windows shaded, P_main --------
    {
        TCanvas* c = new TCanvas("c_zoom", "", 2100, 1300);
        c->Divide(3, 2, 0.001, 0.001);
        for (int p = 0; p < 6; ++p) {
            const int s = idx(pan[p], 0);
            const Spectra& q = sp[s];
            TH1D* hs = Scaled(q.Hist("sig", Form("zoom_sig_%s", base[pan[p]].tag.c_str())), 1e-3);
            Deco(hs, cS, 20, 0.8);
            fout->cd(); hs->Write();
            double ymax = 0;
            for (int i = hs->FindBin(0.4201); i <= hs->FindBin(0.6199); ++i) ymax = std::max(ymax, hs->GetBinContent(i) + hs->GetBinError(i));
            c->cd(p + 1);
            const double y2 = 1.32 * ymax;
            Frame(Form("fr_zoom_%d", p), 0.42, 0.62, 0, y2, xM4, "N #minus CB [10^{3} / 2.5 MeV/c^{2}]");
            Band(scan::kBumpLo, scan::kBumpHi, 0, y2, colBump);
            Band(scan::kPeakLo, scan::kPeakHi, 0, y2, colPeak);
            if (frozenPan[p]) {
                VLine(kSimBumpLo, 0, 0.74 * y2, kGreen + 3, 2); VLine(kSimBumpHi, 0, 0.74 * y2, kGreen + 3, 2);
                VLine(kSimPeakLo, 0, 0.74 * y2, kAzure + 3, 2); VLine(kSimPeakHi, 0, 0.74 * y2, kAzure + 3, 2);
            }
            hs->Draw("e same");
            const PeakFit& f = fit[s];
            Curve(f, false, f.cfg.lo, f.cfg.hi, kMagenta + 2, 2, 1e-3)->Draw("l same");
            Curve(f, true, f.cfg.lo, f.cfg.hi, kMagenta + 2, 1, 1e-3)->Draw("l same");
            gPad->RedrawAxis();
            Txt(0.16, 0.935, Form("(%c) %s", 'a' + p, panTitle[p]), 0.054);
            Txt(0.20, 0.845, Form("P_{main}: S = %.0f #pm %.0f #pm %.0f", f.S, f.eS, spr[s].half), 0.046, 11, kMagenta + 2);
            Txt(0.20, 0.775, Form("bump window: N #minus CB = %.0f #pm %.0f", q.win[1].sig(), q.win[1].sigErr()), 0.046, 11, kGreen + 3);
            if (p == 0) {
                TLegend* lg = Leg(0.64, 0.17, 0.97, 0.36, 0.038);
                LE(lg, hs, "N #minus CB", "pe");
                TGraph* gl = new TGraph(); gl->SetLineColor(kMagenta + 2); gl->SetLineWidth(2);
                TGraph* gl2 = new TGraph(); gl2->SetLineColor(kMagenta + 2); gl2->SetLineWidth(2); gl2->SetLineStyle(2);
                LE(lg, gl, "P_{main} fit", "l");
                LE(lg, gl2, "its pol2 part", "l");
                TLine* ld = new TLine(); ld->SetLineStyle(2); ld->SetLineWidth(2); ld->SetLineColor(kGray + 2);
                LE(lg, ld, "sim. windows", "l");
                lg->Draw();
            }
        }
        Save(c, "spec_m4_zoom");
    }

    // ---- C3. isBest >= 0 versus isBest == 1 ---------------------------------------------
    {
        TCanvas* c = new TCanvas("c_best", "", 2400, 1250);
        c->Divide(4, 2, 0.001, 0.001);
        for (int p = 0; p < 6; ++p) {
            const int s0 = idx(pan[p], 0), s1 = idx(pan[p], 2);
            TH1D* h0 = Scaled(sp[s0].Hist("sig", Form("best_sig0_%s", base[pan[p]].tag.c_str())), 1e-3);
            TH1D* h1 = Scaled(sp[s1].Hist("sig", Form("best_sig1_%s", base[pan[p]].tag.c_str())), 1e-3);
            Deco(h0, cS, 20, 0.75); Deco(h1, kCyan + 3, 24, 0.75);
            fout->cd(); h0->Write(); h1->Write();
            double ymax = 0;
            for (int i = h0->FindBin(0.4201); i <= h0->FindBin(0.6199); ++i) ymax = std::max(ymax, h0->GetBinContent(i) + h0->GetBinError(i));
            c->cd(p + 1);
            const double y2 = 1.38 * ymax;
            Frame(Form("fr_best_%d", p), 0.42, 0.62, 0, y2, xM4, "N #minus CB [10^{3} / 2.5 MeV/c^{2}]");
            Band(scan::kBumpLo, scan::kBumpHi, 0, y2, colBump);
            Band(scan::kPeakLo, scan::kPeakHi, 0, y2, colPeak);
            h0->Draw("e same"); h1->Draw("e same");
            gPad->RedrawAxis();
            Txt(0.16, 0.935, Form("(%c) %s", 'a' + p, panTitle[p]), 0.052);
            Txt(0.20, 0.85, Form("isBest #geq 0: S = %.0f #pm %.0f", fit[s0].S, fit[s0].eS), 0.044, 11, cS);
            Txt(0.20, 0.79, Form("isBest = 1: S = %.0f #pm %.0f", fit[s1].S, fit[s1].eS), 0.044, 11, kCyan + 3);
            if (p == 0) {
                TLegend* lg = Leg(0.52, 0.17, 0.96, 0.33, 0.044);
                LE(lg, h0, "isBest #geq 0", "pe");
                LE(lg, h1, "isBest = 1", "pe");
                lg->Draw();
            }
        }
        // panel 7: kept fractions of isBest == 1 within pre, versus M4 (10 MeV)
        {
            const Spectra& a = sp[idx(0, 0)]; const Spectra& b = sp[idx(0, 2)];
            TGraphErrors* gN = Gr(cN, 20); TGraphErrors* gC = Gr(cCB, 21);
            for (int j = 0; j < a.nw; ++j) {
                const double m = a.lo + (j + 0.5) * a.reb * a.bw;
                if (m < 0.35 || m > 1.10) continue;
                if (a.w[0][j] > 0) {
                    const double f = b.w[0][j] / a.w[0][j];
                    gN->SetPoint(gN->GetN(), m, f); gN->SetPointError(gN->GetN() - 1, 0, std::sqrt(f * (1 - f) / a.w[0][j]));
                }
                if (a.cbW[j] > 0 && b.cbW[j] > 0) {
                    const double fp = b.w[1][j] / a.w[1][j], fm = b.w[2][j] / a.w[2][j];
                    const double f = b.cbW[j] / a.cbW[j];
                    const double er = 0.5 * std::sqrt((1 - fp) / (fp * a.w[1][j]) + (1 - fm) / (fm * a.w[2][j]));
                    gC->SetPoint(gC->GetN(), m, f); gC->SetPointError(gC->GetN() - 1, 0, f * er);
                }
            }
            fout->cd(); gN->Write("best1_share_N"); gC->Write("best1_share_CB");
            c->cd(7);
            Frame("fr_best_frac", 0.35, 1.10, 0.4, 1.15, xM4, "isBest = 1 / isBest #geq 0");
            gN->Draw("p same"); gC->Draw("p same");
            Txt(0.16, 0.935, "(g) pre: share of isBest = 1 (10 MeV/c^{2})", 0.050);
            TLegend* lg = Leg(0.20, 0.20, 0.60, 0.36, 0.046);
            LE(lg, gN, "N", "pe"); LE(lg, gC, "CB", "pe"); lg->Draw();
        }
        // panel 8: log overview of pre, isBest == 1
        {
            const Spectra& q = sp[idx(0, 2)];
            TH1D* hn = Merge(q, "us", 2, "ov1_us"); TH1D* hc = Merge(q, "cb", 2, "ov1_cb"); TH1D* hs = Merge(q, "sig", 2, "ov1_sig");
            Deco(hn, cN, 20, 0.5); Deco(hc, cCB, 20, 0.5); Deco(hs, cS, 20, 0.5);
            fout->cd(); hn->Write(); hc->Write(); hs->Write();
            double ymax = 0;
            for (int i = hn->FindBin(0.3501); i <= hn->FindBin(1.0999); ++i) ymax = std::max(ymax, hn->GetBinContent(i));
            c->cd(8); gPad->SetLogy();
            Frame("fr_best_ov", 0.35, 1.10, 5, 40 * ymax, xM4, "rows / 5 MeV/c^{2}");
            Band(scan::kBumpLo, scan::kBumpHi, 5, 40 * ymax / 30, colBump);
            Band(scan::kPeakLo, scan::kPeakHi, 5, 40 * ymax / 30, colPeak);
            hc->Draw("e same"); hn->Draw("e same"); hs->Draw("e same");
            gPad->RedrawAxis();
            Txt(0.16, 0.935, "(h) pre, isBest = 1", 0.050);
            TLegend* lg = Leg(0.54, 0.17, 0.96, 0.38, 0.044);
            LE(lg, hn, "N", "pe"); LE(lg, hc, "CB", "pe"); LE(lg, hs, "N #minus CB", "pe"); lg->Draw();
        }
        Save(c, "spec_m4_isbest");
    }

    // ---- C4. like-sign composition ----------------------------------------------------------
    {
        TCanvas* c = new TCanvas("c_ls", "", 2100, 700);
        c->Divide(3, 1, 0.001, 0.001);
        // (a) N, N++, N--, CB for pre, 10 MeV, log
        {
            const Spectra& q = sp[idx(0, 0)];
            TH1D* hN = Merge(q, "us", 4, "ls_us"); TH1D* hP = Merge(q, "pp", 4, "ls_pp"); TH1D* hM = Merge(q, "mm", 4, "ls_mm");
            TH1D* hC = Merge(q, "cb", 4, "ls_cb");
            Deco(hN, cN, 20, 0.7); Deco(hP, kGreen + 2, 25, 0.8); Deco(hM, kOrange + 7, 26, 0.9); Deco(hC, cCB, 20, 0.7);   // [review fix] e+e+ green, e-e- orange
            fout->cd(); hN->Write(); hP->Write(); hM->Write(); hC->Write();
            double ymax = 0;
            for (int i = hN->FindBin(0.3501); i <= hN->FindBin(1.0999); ++i) ymax = std::max(ymax, hN->GetBinContent(i));
            c->cd(1); gPad->SetLogy();
            Frame("fr_ls_a", 0.35, 1.10, 1, 30 * ymax, xM4, "counts / 10 MeV/c^{2}");
            hN->Draw("e same"); hC->Draw("e same"); hM->Draw("e same"); hP->Draw("e same");
            Txt(0.16, 0.935, "(a) pre: unlike-sign and like-sign counts", 0.050);
            TLegend* lg = Leg(0.54, 0.17, 0.96, 0.43, 0.044);
            LE(lg, hN, "N (e^{+}e^{-})", "pe"); LE(lg, hM, Form("%s (e^{-}e^{-})", kNmm), "pe");
            LE(lg, hP, "N_{++} (e^{+}e^{+})", "pe"); LE(lg, hC, Form("CB = 2#sqrt{N_{++}%s}", kNmm), "pe");
            lg->Draw();
        }
        // (b) N++/N-- and (c) CB/N for the selections, 20 MeV bins
        const int nSel = 6;
        // [review fix] one selection key with cb_share: pre black, pre + vertex blue, peak nom. + vertex violet,
        // bump nom. + vertex olive, PEAK orange, DIRECT green (red is reserved for the CB)
        const int sCol[nSel] = {kBlack, kBlue + 1, kViolet + 1, kYellow + 3, kOrange + 7, kGreen + 2};
        const int sMk[nSel] = {20, 21, 33, 34, 23, 22};
        TGraphErrors* gR[nSel]; TGraphErrors* gF[nSel];
        for (int p = 0; p < nSel; ++p) {
            const Spectra& q = sp[idx(pan[p], 0)];
            gR[p] = Gr(sCol[p], sMk[p], 0.9); gF[p] = Gr(sCol[p], sMk[p], 0.9);
            for (int j = 0; j + 1 < q.nw; j += 2) {           // 2 x 10 MeV
                const double m = q.lo + (j + 1) * q.reb * q.bw;
                if (m < 0.35 || m > 1.10) continue;
                const double np = q.w[1][j] + q.w[1][j + 1], nm = q.w[2][j] + q.w[2][j + 1], nn = q.w[0][j] + q.w[0][j + 1];
                const double cb = q.cbW[j] + q.cbW[j + 1], cbv = q.cbWVar[j] + q.cbWVar[j + 1];
                const double dx = 0.003 * (p - 2.5);
                if (np > 0 && nm > 0) {
                    const double r = nm / np;   // [review fix] N-- / N++ as everywhere in the report
                    gR[p]->SetPoint(gR[p]->GetN(), m + dx, r); gR[p]->SetPointError(gR[p]->GetN() - 1, 0, r * std::sqrt(1 / np + 1 / nm));
                }
                if (nn > 0 && cb > 0) {
                    const double f = cb / nn;
                    gF[p]->SetPoint(gF[p]->GetN(), m + dx, f); gF[p]->SetPointError(gF[p]->GetN() - 1, 0, f * std::sqrt(cbv / (cb * cb) + 1 / nn));
                }
            }
            fout->cd(); gR[p]->Write(Form("ratio_pp_mm_%s", base[pan[p]].tag.c_str())); gF[p]->Write(Form("cb_over_n_%s", base[pan[p]].tag.c_str()));
        }
        c->cd(2);
        Frame("fr_ls_b", 0.35, 1.10, 0.0, 14.0, xM4, Form("%s / N_{++}", kNmm));
        Band(scan::kBumpLo, scan::kBumpHi, 0.0, 14.0, colBump); Band(scan::kPeakLo, scan::kPeakHi, 0.0, 14.0, colPeak);
        for (int p = 0; p < nSel; ++p) gR[p]->Draw("p same");
        gPad->RedrawAxis();
        Txt(0.16, 0.935, "(b) like-sign ratio (20 MeV/c^{2} bins)", 0.050);
        c->cd(3);
        Frame("fr_ls_c", 0.35, 1.10, 0.0, 0.45, xM4, "CB / N");
        Band(scan::kBumpLo, scan::kBumpHi, 0.0, 0.45, colBump); Band(scan::kPeakLo, scan::kPeakHi, 0.0, 0.45, colPeak);
        for (int p = 0; p < nSel; ++p) gF[p]->Draw("p same");
        gPad->RedrawAxis();
        Txt(0.16, 0.935, "(c) CB share (20 MeV/c^{2} bins)", 0.050);
        TLegend* lg = Leg(0.40, 0.62, 0.96, 0.90, 0.040);
        for (int p = 0; p < nSel; ++p) LE(lg, gF[p], panTitle[p], "pe");
        lg->Draw();
        Save(c, "spec_likesign");
    }

    // ---- C5. vertex cuts: unlike-sign vs like-sign ---------------------------------------
    {
        TCanvas* c = new TCanvas("c_vtx", "", 1700, 1300);
        c->Divide(2, 2, 0.001, 0.001);
        const int kCol[3] = {kBlack, kGreen + 2, kOrange + 7};   // [review fix] e+e+ green squares, e-e- orange triangles
        const int kMk[3] = {20, 25, 26};
        // (a) vz shapes
        c->cd(1); gPad->SetLogy();
        double ym = 0;
        for (int k = 0; k < 3; ++k) {
            hVz[k]->Scale(1.0 / hVz[k]->Integral(), "width");     // per mm
            Deco(hVz[k], kCol[k], kMk[k], 0.6);
            ym = std::max(ym, hVz[k]->GetMaximum());
            fout->cd(); hVz[k]->Write();
        }
        Frame("fr_vz", -500, 100, 1e-6, 20 * ym, "event vertex v_{z} [mm]", "fraction of rows / mm");
        for (int k = 0; k < 3; ++k) hVz[k]->Draw("e same");
        VLine(-160, 1e-6, 2 * ym, kRed + 1, 2);
        Txt(0.16, 0.935, "(a) v_{z} of the rows (cor_ok, isBest #geq 0, no v_{z} cut at #minus160 mm)", 0.044);
        TLegend* lg = Leg(0.20, 0.66, 0.60, 0.88, 0.044);
        for (int k = 0; k < 3; ++k) LE(lg, hVz[k], smpRoot[k], "pe");
        lg->Draw();
        Txt(0.42, 0.20, "#leftarrow removed by v_{z} > #minus160 mm", 0.042, 11, kRed + 1);
        // (b)-(d) pass fractions versus M4
        const char* lab[nCut] = {"(b) ", "(c) ", "(d) "};
        const char* baseTxt[nCut] = {"base: cor_ok, isBest #geq 0", "base: pre", "base: pre"};
        for (int cc = 0; cc < nCut; ++cc) {
            c->cd(2 + cc);
            Frame(Form("fr_pass_%d", cc), 0.35, 1.10, 0.0, 1.25, xM4, "row pass fraction");
            Band(scan::kBumpLo, scan::kBumpHi, 0.0, 1.0, colBump); Band(scan::kPeakLo, scan::kPeakHi, 0.0, 1.0, colPeak);
            TGraphErrors* gk[3];
            for (int k = 0; k < 3; ++k) {
                TGraphErrors* g = Gr(kCol[k], kMk[k], 0.75);
                for (int b = 0; b < fNb; ++b) {
                    if (fM4[cc][k].t[b] < 5) continue;
                    const double m = fLo + (b + 0.5) * (fHi - fLo) / fNb;
                    g->SetPoint(g->GetN(), m, fM4[cc][k].F(b)); g->SetPointError(g->GetN() - 1, 0, fM4[cc][k].E(b));
                }
                g->Draw("p same");
                gk[k] = g;
                fout->cd(); g->Write(Form("pass_%d_%s", cc, scan::SampleName(k)));
            }
            gPad->RedrawAxis();
            Txt(0.16, 0.935, Form("%s%s (%s)", lab[cc], cutRoot[cc], baseTxt[cc]), 0.046);
            TLegend* l3 = Leg(0.18, 0.80, 0.96, 0.89, 0.042);
            l3->SetNColumns(3);
            for (int k = 0; k < 3; ++k) LE(l3, gk[k], smpRoot[k], "pe");
            l3->Draw();
        }
        Save(c, "spec_vertexcut");
    }
    fout->Close();
    printf("spec_m4 done in %.1f s\n", sw.RealTime());
}
