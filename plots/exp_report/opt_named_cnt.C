// opt_named_cnt.C -- COPY of plots/exp_report/opt_named.C made by figure agent A (integ4, 2026-11-15) with
// the figure strings 'rows' -> 'counts' (analyst's terminology) and the outputs of the copy redirected:
// figures of the scope -> docs/latex/figs_exp/<base>_cnt.{pdf,png}; everything else -> /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figA/work/opt_named
// Built by /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figA/work/mk_cnt.py; the physics code is unchanged.
// ============================================================================
// opt_named.C  --  all named selections on the data: eta peak, bump, tagged eta,
//                  X17 band; frozen simulation selections against the
//                  data-optimised ones
// ============================================================================
// HADES p+p, T = 4.53 GeV (feb22), channel pi+ pi- e+ e-. EXPERIMENTAL DATA ONLY
// (the simulation enters only through the FROZEN cut values PEAK, DIRECT, BUMP
// and X17, which were chosen on the simulation alone by a freeze rule declared
// before the search; these data were not used in that choice).
// Edge-step sums E whose summed edge fit has p < 1e-3 (band-pass cuts: slides,
// old code) are marked as invalid in the tables and in the summary figure
// (verifier change, 2026-11-14).
//
// What it does
//   For every named selection (see opt_common.h, section B) and the four variants
//   {isBest >= 0, isBest == 1} x {no OA cut, OA(e+e-) > 2 deg}:
//   1. eta peak (P_main: N - CB, frozen Gaussian 0.5439 GeV/c2 / 5.4 MeV + pol2,
//      0.5225-0.6100): S +- stat +- model spread (48 variations, p > 0.05),
//      S/B in the DATA window 0.53364-0.55416, F = S/stat, F_tot = S/sqrt(stat^2 +
//      spread^2), S with CB x 0.97 and x 1.03, eps_S = S/S(pre);
//   2. bump: edge-step visibility E +- stat +- model spread (four M(pi+pi-) slices,
//      scan_bump.C procedure), N and N - CB in the DATA bump window 0.4550-0.5225,
//      Q = E/sqrt(N), E with CB x 0.97 / 1.03;
//   3. ECAL-tagged eta T (pi0 window 0.119-0.173 in M(e+e-gamma), m5c fit) +- stat
//      +- spread and eff_T = T/T(pre): an independent check of the bump selections
//      (ECAL is used for verification only);
//   4. the two halves of the event groups (fold 0 / 1) for the main variant;
//   5. the X17 band (frozen: pre + both pions in the vertex fit + 15.5 <= M(ee) <
//      24.5 MeV/c2): rows, CB and N - CB in the simulation window 0.5333-0.5587 and
//      in the data window, and the eta content of the band (P_main);
//   6. visible yields per luminosity N x 64/63 / L (L = 6079.2 nb^-1, runs in the
//      data; not corrected for acceptance or efficiency).
//   Figures: M4 spectra of four selections (pre, peak nominal + vertex, DIRECT,
//   PEAK), the X17 band, the tagged eta for four bump selections, and a summary
//   chart of S, F, E and T of all named selections.
//
// Selections (all rows: cor_ok == 1, isBest >= 0; FAT: vertex z > -500 mm, PT3)
//   pre                   = vz > -160 mm                       (= frozen BUMP)
//   peak nominal          = pre + oa_lab < 55 + M(pipi) > 0.30
//   peak nominal + vertex = + both pions in the vertex fit (pip/pim_isUsedVertex)
//   bump nominal          = pre + oa_lab < 40 + 0.30 < M(pipi) < 0.41
//   bump nominal + vertex = + both pions in the vertex fit
//   DIRECT (frozen)       = pre + dz_pi < 40 mm
//   PEAK (frozen)         = pre + dz_pi < 40 mm + MM4 < 2.85 GeV/c2
//   X17 (frozen)          = pre + both pions in the vertex fit + 15.5 <= M(ee) < 24.5 MeV/c2
//   + references: pre without vz > -160, pre + vertex, peak/bump nominal + dz_pi < 40,
//     slides, old code
//
// Inputs : /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/scan_v2/scan_compact.root
//          (scan::kCompact, with dz_pi) through ../ScanEngine.h and opt_common.h
// Outputs: /home/przygoda/HADES/FAT/FAT/docs/latex/tables_exp/opt_named.txt (all numbers),
//          opt_named.tex, opt_named_variants.tex (LaTeX tables)
//          /home/przygoda/HADES/FAT/FAT/docs/latex/figs_exp/opt_named_spectra.pdf, opt_named_x17.pdf,
//          opt_named_tagged.pdf, opt_named_summary.pdf (+ .png)
//          /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/opt/opt_named.root (spectra)
//
// Run (from the repository root /home/przygoda/HADES/FAT/FAT), about 3 min:
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/aclic_opt", kTRUE)' 'plots/exp_report/opt_named.C+'
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "TBox.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TSystem.h"

#include "opt_common_cnt.h"

namespace optnamed {

using namespace opt;

const int kNVar = 4;   // 0: isBest >= 0, 1: isBest >= 0 + OA > 2, 2: isBest == 1, 3: isBest == 1 + OA > 2
const char* const kVarName[kNVar] = {"isBest>=0", "isBest>=0, OA>2", "isBest==1", "isBest==1, OA>2"};
inline Cuts Variant(const Cuts& c, int v) {
    Cuts o = c;
    if (v == 1 || v == 3) o = WithOA2(o);
    if (v >= 2) o = WithBest(o);
    return o;
}

struct Entry {
    Cuts c;                  // base definition (variant 0)
    std::string tex;         // LaTeX name
    int kind = 0;            // 0 data-optimised / reference, 1 frozen simulation selection
    PeakRes pk[kNVar];
    BumpRes bp[kNVar];
    PeakFit fold[2]; BumpRes bfold[2];
    double epsS[kNVar], eEpsS[kNVar], effT[kNVar], eEffT[kNVar];
};

// the summed edge fit of E is bad (p < 1e-3): band-pass cut in M4, E is not a valid measure
// (the same criterion as BumpFlag in opt_scans.C)
inline bool EdgeBad(const BumpRes& b) { return b.ndf > 0 && TMath::Prob(b.chi2, b.ndf) < 1e-3; }

// axis fonts in pixels (font 43): the same size in every sub-pad
void AxPx(TH1* h, int lab = 20, int tit = 22, double xoff = 3.2, double yoff = 2.0) {
    for (TAxis* ax : {h->GetXaxis(), h->GetYaxis()}) {
        ax->SetLabelFont(43); ax->SetTitleFont(43); ax->SetLabelSize(lab); ax->SetTitleSize(tit);
    }
    h->GetXaxis()->SetTitleOffset(xoff); h->GetYaxis()->SetTitleOffset(yoff);
}
TLatex* TxtPx(double x, double y, const char* t, int px = 20, int col = kBlack, int align = 11) {
    TLatex* l = new TLatex(x, y, t); l->SetNDC(); l->SetTextFont(43); l->SetTextSize(px);
    l->SetTextColor(col); l->SetTextAlign(align); l->Draw(); return l;
}

// draw one M4 spectrum into the cell pad p: top = N black, CB red, N - CB blue, fit;
// bottom = residual N - CB - pol2 with the fitted Gaussian (rows / 2.5 MeV)
const double kX1 = 0.49, kX2 = 0.62;
void DrawSpectrum(TPad* p, const Spectra& sp, const PeakRes& r, const char* title, double ymax, double rlo, double rhi,
                  bool legend, const char* tag, bool simWin = false) {
    p->cd();
    TPad* pt = new TPad(Form("pt_%s", tag), "", 0, 0.34, 1, 1);
    TPad* pb = new TPad(Form("pb_%s", tag), "", 0, 0, 1, 0.34);
    pt->SetLeftMargin(0.16); pt->SetRightMargin(0.03); pt->SetTopMargin(0.10); pt->SetBottomMargin(0.02);
    pb->SetLeftMargin(0.16); pb->SetRightMargin(0.03); pb->SetTopMargin(0.03); pb->SetBottomMargin(0.36);
    pt->Draw(); pb->Draw();
    // ---- top ----
    pt->cd();
    TH1D* fr = Frame(Form("fr_%s", tag), kX1, kX2, 0, ymax, "", "counts / 2.5 MeV/c^{2}");
    AxPx(fr); fr->GetXaxis()->SetLabelSize(0); fr->GetYaxis()->SetMaxDigits(6);
    TBox* bw = new TBox(scan::kPeakLo, 0, scan::kPeakHi, ymax); bw->SetFillColor(TColor::GetColor("#e6e6e6")); bw->SetLineWidth(0); bw->Draw();
    if (simWin) { TBox* b = new TBox(opt::kSimPeakLo, 0, opt::kSimPeakHi, ymax); b->SetFillStyle(3354); b->SetFillColor(TColor::GetColor("#f7b58c")); b->SetLineWidth(0); b->Draw(); }
    TH1D* hn = sp.Hist("us", Form("hn_%s", tag)); TH1D* hc = sp.Hist("cb", Form("hc_%s", tag)); TH1D* hs = sp.Hist("sig", Form("hs_%s", tag));
    for (TH1D* h : {hn, hc, hs}) h->GetXaxis()->SetRangeUser(kX1, kX2 - 1e-6);
    hn->SetMarkerStyle(20); hn->SetMarkerSize(0.6); hn->SetMarkerColor(kBlack); hn->SetLineColor(kBlack);
    hc->SetMarkerStyle(21); hc->SetMarkerSize(0.55); hc->SetMarkerColor(kRed + 1); hc->SetLineColor(kRed + 1);
    hs->SetMarkerStyle(20); hs->SetMarkerSize(0.6); hs->SetMarkerColor(kBlue + 1); hs->SetLineColor(kBlue + 1);
    TGraph* gt = new TGraph(); TGraph* gb = new TGraph(); TGraph* gg = new TGraph();
    for (int i = 0; i <= 200; ++i) {
        const double m = scan::kFitLo + i * (scan::kFitHi - scan::kFitLo) / 200.;
        gt->SetPoint(i, m, r.f.Total(m)); gb->SetPoint(i, m, r.f.Bkg(m)); gg->SetPoint(i, m, r.f.Total(m) - r.f.Bkg(m));
    }
    // [review fix 2026-11-14] fit magenta, pol2 grey dashed (the key of the other chapters; was blue on blue points)
    gt->SetLineColor(kMagenta + 1); gt->SetLineWidth(3); gb->SetLineColor(kGray + 2); gb->SetLineStyle(2); gb->SetLineWidth(2);
    gg->SetLineColor(kMagenta + 1); gg->SetLineWidth(3);
    hn->Draw("PE SAME"); hc->Draw("PE SAME"); hs->Draw("PE SAME"); gb->Draw("L SAME"); gt->Draw("L SAME");
    fr->Draw("AXIS SAME");
    TxtPx(0.16, 0.93, title, 22);
    TxtPx(0.19, 0.80, Form("S = %.0f #pm %.0f (stat) #pm %.0f (model, old envelope)", r.f.S, r.f.eS, r.spr.half), 19, kBlue + 2);
    TxtPx(0.19, 0.71, Form("S/B = %.3f,  F = S/#sigma_{S} = %.1f", r.f.SoB, r.f.F), 19, kBlue + 2);
    if (legend) {
        TLegend* l = new TLegend(0.19, 0.30, 0.62, 0.68);
        l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextFont(43); l->SetTextSize(18);
        l->AddEntry(hn, "all counts N", "pe"); l->AddEntry(hc, "CB = 2#sqrt{N_{++}N_{--}}", "pe");
        l->AddEntry(hs, "N - CB", "pe"); l->AddEntry(gt, "fit: Gaussian + pol2", "l"); l->AddEntry(gb, "pol2 part", "l");
        l->Draw();
    }
    // ---- bottom: residual ----
    pb->cd();
    TH1D* fb = Frame(Form("fb_%s", tag), kX1, kX2, rlo, rhi, "M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]", "N - CB - pol2");
    AxPx(fb, 20, 22, 1.45, 2.0); fb->GetYaxis()->SetNdivisions(504);
    TBox* bw2 = new TBox(scan::kPeakLo, rlo, scan::kPeakHi, rhi); bw2->SetFillColor(TColor::GetColor("#e6e6e6")); bw2->SetLineWidth(0); bw2->Draw();
    TGraphErrors* gr = GE(kBlue + 1, 20, 0.6);
    for (int i = 0; i < sp.nbin; ++i) {
        const double x = sp.x(i);
        if (x < scan::kFitLo || x >= scan::kFitHi) continue;
        AddPt(gr, x, sp.sig[i] - r.f.Bkg(x), std::sqrt(sp.sigVar[i]));
    }
    HLine(kX1, kX2, 0, kGray + 2, 2, 1);
    gr->Draw("PZ SAME"); gg->Draw("L SAME");
    fb->Draw("AXIS SAME");
}

} // namespace optnamed

// ============================================================================
void opt_named_cnt() {
    using namespace opt;
    using namespace optnamed;
    TStopwatch sw;
    gSystem->mkdir(kFigOut, kTRUE); gSystem->mkdir(kTabOut, kTRUE); gSystem->mkdir(kWork, kTRUE); gSystem->mkdir(kCntScratch, kTRUE);
    Style();

    scan::Data d;
    if (!scan::Load(d, scan::kCompact, scan::kMlo, scan::kMhi)) return;
    scan::Data dT;
    Filter(d, dT, [](const Row& r) { return HasGam(r) && r.m_eeg >= 0.05f && r.m_eeg < 0.25f; });
    printf("photon sample: %zu / %zu / %zu rows\n", dT.s[0].size(), dT.s[1].size(), dT.s[2].size());

    // ---- the named selections ----
    std::vector<Entry> E;
    auto add = [&](const Cuts& c, const char* tex, int kind) { Entry e; e.c = c; e.tex = tex; e.kind = kind; E.push_back(e); };
    add(Pre0(),      "pre without $v_z>-160$", 0);
    add(Pre(),       "pre = BUMP (frozen)", 1);
    add(PeakNom(),   "peak nominal", 0);
    add(PeakNomV(),  "peak nominal + vertex", 0);
    add(PeakNomDz(), "peak nominal + $\\Delta z_\\pi<40$", 0);
    add(Direct(),    "DIRECT (frozen)", 1);
    add(PeakF(),     "PEAK (frozen)", 1);
    add(PreV(),      "pre + vertex", 0);
    add(BumpNom(),   "bump nominal", 0);
    add(BumpNomV(),  "bump nominal + vertex", 0);
    add(BumpNomDz(), "bump nominal + $\\Delta z_\\pi<40$", 0);
    add(Slides(),    "slides", 0);
    add(OldCode(),   "old code", 0);
    add(X17(),       "X17 band (frozen)", 1);
    const int ne = (int)E.size();
    const int iPre = 1, iX17 = ne - 1;

    // ---- peak spectra: all variants + folds, one pass ----
    std::vector<Sel> sels;
    for (int i = 0; i < ne; ++i) {
        for (int v = 0; v < kNVar; ++v) sels.push_back(MkSel(Variant(E[i].c, v)));
        for (int f = 0; f < 2; ++f) sels.push_back(MkSel(WithFold(E[i].c, f)));
    }
    const std::vector<Spectra> sp = scan::FillMany(d, sels, OptWin());
    printf("peak spectra filled (%zu) at %.1f s\n", sp.size(), sw.RealTime()); sw.Continue();
    for (int i = 0; i < ne; ++i) {
        for (int v = 0; v < kNVar; ++v) E[i].pk[v] = EvalPeak(sp[i * (kNVar + 2) + v]);
        for (int f = 0; f < 2; ++f) E[i].fold[f] = scan::FitPeakMain(sp[i * (kNVar + 2) + kNVar + f]);
    }
    printf("peak fits done at %.1f s\n", sw.RealTime()); sw.Continue();

    // ---- bump: E (+ spread, CB +-3 %) and T for all variants; folds without spread ----
    {
        std::vector<Cuts> cs;
        for (int i = 0; i < ne; ++i) for (int v = 0; v < kNVar; ++v) cs.push_back(Variant(E[i].c, v));
        const std::vector<BumpRes> br = EvalBump(d, cs, true, true, &dT, true);
        for (int i = 0; i < ne; ++i) for (int v = 0; v < kNVar; ++v) E[i].bp[v] = br[i * kNVar + v];
        std::vector<Cuts> cf;
        for (int i = 0; i < ne; ++i) for (int f = 0; f < 2; ++f) cf.push_back(WithFold(E[i].c, f));
        const std::vector<BumpRes> bf = EvalBump(d, cf, false, false);
        for (int i = 0; i < ne; ++i) for (int f = 0; f < 2; ++f) E[i].bfold[f] = bf[i * 2 + f];
    }
    printf("bump fits done at %.1f s\n", sw.RealTime()); sw.Continue();

    // ---- efficiencies relative to pre (same isBest / OA variant) ----
    for (int i = 0; i < ne; ++i)
        for (int v = 0; v < kNVar; ++v) {
            const PeakRes& p0 = E[iPre].pk[v];
            SubRatio(E[i].pk[v].f.S, E[i].pk[v].f.eS, p0.f.S, p0.f.eS, E[i].epsS[v], E[i].eEpsS[v]);
            const BumpRes& b0 = E[iPre].bp[v];
            SubRatio(E[i].bp[v].T.S, E[i].bp[v].T.eS, b0.T.S, b0.T.eS, E[i].effT[v], E[i].eEffT[v]);
        }

    // ========================================================================
    // tables
    // ========================================================================
    Out o;
    o.Open(std::string(kTabOut) + "/opt_named.txt");
    o.P("# opt_named.txt  --  written by plots/exp_report/opt_named.C\n");
    o.P("# EXPERIMENTAL DATA, p+p 4.53 GeV feb22. Rows = stored particle combinations, cor_ok == 1, corrected momenta (COR).\n");
    o.P("# FAT production: vertex z > -500 mm, PT3 (trigbit 8192). pre = isBest >= 0 && vz > -160 mm. CB = 2 sqrt(N++ N--), k = 1.\n");
    o.P("# PEAK: P_main = chi2 fit of N - CB, frozen Gaussian %.4f GeV/c2 / %.1f MeV + pol%d in %.4f-%.4f; S = whole Gaussian;\n",
        scan::kPeakMean, 1e3 * scan::kPeakSigma, scan::kFitOrder, scan::kFitLo, scan::kFitHi);
    o.P("#   spread = half range of S over the 48 variations (order 1-4 x 4 ranges x width +-0.5 MeV) with p > 0.05 (n kept in brackets);\n");
    o.P("#   S/B = S in the DATA window %.5f-%.5f / (pol2 + CB) there; F = S/stat (compares selections, NOT a significance);\n", scan::kPeakLo, scan::kPeakHi);
    o.P("#   Ftot = S/sqrt(stat^2 + spread^2); S(CB-3%%), S(CB+3%%) = S with the CB scaled by 0.97 / 1.03; eps_S = S / S(pre, same variant).\n");
    o.P("# BUMP: E = sum of the edge-step heights (rows per 2.5 MeV) in M(pipi) 0.30-0.34, 0.34-0.37, 0.37-0.39, 0.39-0.40 (scan_bump.C),\n");
    o.P("#   a visibility measure, NOT a yield; spread = slices added linearly (order 2,3 x half range 35,45,55 MeV x smearing 4.4,5.4,6.4 MeV);\n");
    o.P("#   N_b, W = N - CB in the DATA bump window %.4f-%.4f; Q = E/sqrt(N_b); effN = N_b/N_b(pre); gainE = (E/E(pre))/sqrt(effN),\n"
        "#   gainT = effT/sqrt(effN) (> 1: the cut improves the statistical visibility of the bump; same variant).\n", scan::kBumpLo, scan::kBumpHi);
    o.P("# TAGGED ETA (ECAL, verification only): rows with a good photon, %.3f <= M(e+e-gamma) < %.3f, M(ee) < 0.1; m5c fit %.4f / %.1f MeV, pol2,\n",
        kPi0Lo, kPi0Hi, kTagMean, 1e3 * kTagSigma);
    o.P("#   %.2f-%.2f (P_main); T = Gaussian content +- stat +- spread; eff_T = T / T(pre, same variant).\n", kTagLo, kTagHi);
    o.P("# Errors: stat from the fits (event-group variance in the window sums). Visible yield per luminosity = S x 64/63 / L, L = %.1f nb^-1\n", kLumi);
    o.P("#   (runs present in the data; START-free elastic-scattering luminosity, total syst 6.9 %%, PT3 live fraction: up to -6.8 %% on L).\n");
    o.P("# Frozen selections (PEAK, DIRECT, BUMP, X17) were chosen on the simulation alone (freeze rule declared before the search);\n"
        "#   these data were not used in the choice, so they carry no data-optimisation bias. (The data had been studied before, in the\n"
        "#   first cut-scan study, which had already named |z_pi - vz| < 40 mm as the geometric alternative to the vertex-fit flag.)\n");
    o.P("# E marked '!' in the DETAILS section: summed edge fit p < 1e-3 (band-pass cut in M4), E is NOT a valid visibility measure.\n\n");

    for (int v = 0; v < kNVar; ++v) {
        o.P("==== VARIANT %s ====\n", kVarName[v]);
        o.P("%-28s | %7s %6s | %6s %5s %5s(%2s) | %6s %6s | %6s %5s | %6s %6s | %5s %5s | %6s | %5s %4s %4s(%2s) %6s %7s | %5s %4s %4s %5s %5s | %5s %5s %5s\n",
            "selection", "N_win", "CB/N", "S", "stat", "sprd", "n", "S/B", "F", "Ftot", "chi2", "S(-3%)", "S(+3%)", "epsS", "err", "S/nb",
            "E", "stat", "sprd", "n", "N_b", "W", "T", "stat", "sprd", "effT", "err", "effN", "gainE", "gainT");
        for (int i = 0; i < ne; ++i) {
            const PeakRes& p = E[i].pk[v]; const BumpRes& b = E[i].bp[v];
            const double effN = E[iPre].bp[v].win.n[0] > 0 ? b.win.n[0] / E[iPre].bp[v].win.n[0] : 0;
            const double gE = effN > 0 ? (b.E / E[iPre].bp[v].E) / std::sqrt(effN) : 0, gT = effN > 0 ? E[i].effT[v] / std::sqrt(effN) : 0;
            o.P("%-28s | %7.0f %6.3f | %6.0f %5.0f %5.0f(%2d) | %6.4f %6.2f | %6.2f %5.1f | %6.0f %6.0f | %5.3f %5.3f | %6.3f | %5.0f %4.0f %4.0f(--) %6.0f %7.0f | %5.0f %4.0f %4.0f %5.3f %5.3f | %5.3f %5.3f %5.3f\n",
                E[i].c.name.c_str(), p.f.N, p.cbs, p.f.S, p.f.eS, p.spr.half, p.spr.npass, p.f.SoB, p.f.F, p.Ftot, p.f.chi2,
                p.Sm3, p.Sp3, E[i].epsS[v], E[i].eEpsS[v], p.f.S * kPT3fac / kLumi,
                b.E, b.eE, b.spreadE, b.win.n[0], b.win.sig(), b.T.S, b.T.eS, b.spreadT, E[i].effT[v], E[i].eEffT[v], effN, gE, gT);
        }
        o.P("\n");
    }
    o.P("==== DEFINITIONS ====\n");
    for (int i = 0; i < ne; ++i) o.P("%-28s : %s%s\n", E[i].c.name.c_str(), Describe(E[i].c).c_str(), E[i].kind ? "   [frozen on simulation]" : "");

    o.P("\n==== DETAILS, variant isBest>=0 (no OA cut) ====\n");
    o.P("%-28s | %7s %7s %7s | %8s %6s | %7s %7s | %6s %6s | %s\n", "selection", "N_win", "CB", "B", "Swin", "chi2/n", "E(-3%)", "E(+3%)",
        "Q", "N_simW", "edge-step heights h per slice (0.30-0.34, 0.34-0.37, 0.37-0.39, 0.39-0.40 | 0.40-0.41 | control 0.42-0.46), dchi2");
    for (int i = 0; i < ne; ++i) {
        const PeakRes& p = E[i].pk[0]; const BumpRes& b = E[i].bp[0];
        std::string hs;
        for (int k = 0; k < kNSl; ++k) hs += b.okSl[k] ? Form(" %.0f+-%.0f(%.0f)", b.h[k], b.eh[k], b.dchi2[k]) : " --";
        hs += Form(" | summed edge fit chi2 %.0f/%d%s", b.chi2, b.ndf, EdgeBad(b) ? " ! INVALID (p < 1e-3)" : "");
        o.P("%-28s | %7.0f %7.0f %7.0f | %8.0f %3.0f/%-2d | %7.0f %7.0f | %6.2f %6.0f |%s\n", E[i].c.name.c_str(), p.f.N, p.f.CB, p.f.B,
            p.f.Swin, p.f.chi2, p.f.ndf, b.Em3, b.Ep3, b.Q, p.wSim.n[0], hs.c_str());
    }

    o.P("\n==== HALVES OF THE EVENT GROUPS (fold = parity of the group index), variant isBest>=0 ====\n");
    o.P("%-28s | %6s %5s %6s | %6s %5s %6s | %7s %5s | %6s %4s | %6s %4s | %7s\n", "selection", "S f0", "err", "F f0", "S f1", "err", "F f1",
        "S0-S1", "pull", "E f0", "err", "E f1", "err", "pull E");
    for (int i = 0; i < ne; ++i) {
        const PeakFit& a = E[i].fold[0]; const PeakFit& b = E[i].fold[1];
        const BumpRes& x = E[i].bfold[0]; const BumpRes& y = E[i].bfold[1];
        const double dS = a.S - b.S, eS = std::hypot(a.eS, b.eS), dE = x.E - y.E, eE = std::hypot(x.eE, y.eE);
        o.P("%-28s | %6.0f %5.0f %6.2f | %6.0f %5.0f %6.2f | %7.0f %5.2f | %6.0f %4.0f | %6.0f %4.0f | %7.2f\n", E[i].c.name.c_str(),
            a.S, a.eS, a.F, b.S, b.eS, b.F, dS, eS > 0 ? dS / eS : 0, x.E, x.eE, y.E, y.eE, eE > 0 ? dE / eE : 0);
    }

    // ---- X17 band ----
    {
        o.P("\n==== X17 BAND (frozen: pre + both pions in vertex fit + %.1f <= M(ee) < %.1f MeV/c2) ====\n", kX17Lo, kX17Hi);
        for (int v = 0; v < kNVar; ++v) {
            const PeakRes& p = E[iX17].pk[v];
            const scan::WinSum& ws = p.wSim;
            const double frSim = scan::GaussFrac(kSimPeakLo, kSimPeakHi, scan::kPeakMean, scan::kPeakSigma);
            // polynomial (non-peak, non-CB background) integrated over the simulation window
            double bsim = 0; { const int n = 400; const double h = (kSimPeakHi - kSimPeakLo) / n;
                for (int q = 0; q < n; ++q) bsim += p.f.Bkg(kSimPeakLo + (q + 0.5) * h); bsim *= h / scan::kBinW; }
            o.P("%-16s simulation window %.4f-%.4f: N = %.0f +- %.0f, CB = %.0f +- %.0f, N - CB = %.0f +- %.0f ; eta S inside = %.0f, pol2 B = %.0f\n",
                kVarName[v], kSimPeakLo, kSimPeakHi, ws.n[0], std::sqrt(ws.var[0]), ws.cb, std::sqrt(ws.cbVar), ws.sig(), ws.sigErr(),
                p.f.S * frSim, bsim);
            {
                const scan::WinSum& wh = p.wSh;
                const double fr = scan::GaussFrac(kShPeakLo, kShPeakHi, scan::kPeakMean, scan::kPeakSigma);
                double bsh = 0; { const int n = 400; const double h = (kShPeakHi - kShPeakLo) / n;
                    for (int q = 0; q < n; ++q) bsh += p.f.Bkg(kShPeakLo + (q + 0.5) * h); bsh *= h / scan::kBinW; }
                o.P("%-16s simulation window on the data mass scale %.4f-%.4f (x 0.5439/0.5460): N = %.0f +- %.0f, CB = %.0f +- %.0f, N - CB = %.0f +- %.0f ; eta S inside = %.0f, pol2 B = %.0f\n",
                    "", kShPeakLo, kShPeakHi, wh.n[0], std::sqrt(wh.var[0]), wh.cb, std::sqrt(wh.cbVar), wh.sig(), wh.sigErr(), p.f.S * fr, bsh);
            }
            o.P("%-16s data window %.5f-%.5f: N = %.0f, CB = %.0f (CB/N %.3f), S/B %.4f ; eta content of the band S = %.0f +- %.0f +- %.0f (%d var.), chi2 %.1f/%d\n",
                "", scan::kPeakLo, scan::kPeakHi, p.f.N, p.f.CB, p.cbs, p.f.SoB, p.f.S, p.f.eS, p.spr.half, p.spr.npass, p.f.chi2, p.f.ndf);
        }
        o.P("   (the simulation expected B = 5342 +- 239 rows of physical background at L = 6080.5 nb^-1 in the simulation window; the simulation has no CB:\n");
        o.P("    a direct comparison needs per-channel cross-section scaling and is NOT made here)\n");
    }
    o.P("\nrun time %.1f s\n", sw.RealTime());
    o.Close();
    sw.Continue();

    // ---- LaTeX: main table (isBest >= 0, no OA cut) ----
    {
        FILE* t = fopen((std::string(kTabOut) + "/opt_named.tex").c_str(), "w");
        fprintf(t, "%% opt_named.tex -- written by plots/exp_report/opt_named.C (numbers: tables_exp/opt_named.txt)\n");
        fprintf(t, "\\begin{table}[H]\n\\centering\\footnotesize\n\\setlength{\\tabcolsep}{3pt}\n");
        fprintf(t, "\\caption{Named selections on the data, all stored combinations (isBest $\\geq 0$), no OA cut. pre = $v_z>-160$~mm, corrected momenta (cor\\_ok $=1$), PT3, FAT vertex $z>-500$~mm. "
                   "$\\eta$ peak: $P_\\mathrm{main}$ fit of $N-\\mathrm{CB}$ (frozen Gaussian 0.5439~GeV/$c^2$, $\\sigma=5.4$~MeV, pol2, 0.5225--0.6100), $S$ = whole Gaussian $\\pm$ stat $\\pm$ model spread (old envelope; revised model uncertainty in Table~\\ref{tab:fix_spread}); "
                   "stat and the $E$ errors are per-row Poisson fit errors (no event-group term; the halving test of Table~\\ref{tab:fix_stat} finds those of $E$ too small by a factor 1.2); "
                   "$N$, CB/$N$ and $S/B$ in the data window 0.53364--0.55416~GeV/$c^2$; $F=S/\\sigma_S^\\mathrm{stat}$ compares selections and is not a significance. "
                   "Bump: edge-step visibility $E$ (rows per 2.5~MeV, not a yield) $\\pm$ stat $\\pm$ spread. "
                   "$T$: ECAL-tagged $\\eta\\to\\pi^+\\pi^-\\pi^0$ (verification only), $\\varepsilon_T=T/T(\\mathrm{pre})$. CB $=2\\sqrt{N_{++}N_{--}}$, $k=1$. "
                   "Frozen = chosen on the simulation alone (freeze rule declared before the search); these data were not used in the choice. "
                   "$E$ in brackets with $^\\dagger$: the summed edge fit fails ($p<10^{-3}$, the $\\eta$-rest angle cut is a band-pass in $M_4$), $E$ is not valid. "
                   "Rows, not normalised to luminosity.}\n");
        fprintf(t, "\\label{tab:opt_named}\n\\begin{tabular}{lrrrrrrrr}\n\\toprule\n");
        fprintf(t, "Selection & $N$ & CB/$N$ & $S\\pm$stat$\\pm$spread & $S/B$ & $F$ & $E\\pm$stat$\\pm$spread & $T\\pm$stat & $\\varepsilon_T$ \\\\\n\\midrule\n");
        for (int i = 0; i < ne; ++i) {
            if (i == iX17) continue;
            const PeakRes& p = E[i].pk[0]; const BumpRes& b = E[i].bp[0];
            const std::string es = EdgeBad(b) ? std::string(Form("(%.0f)$^\\dagger$", b.E))
                                              : std::string(Form("$%.0f\\pm%.0f\\pm%.0f$", b.E, b.eE, b.spreadE));
            fprintf(t, "%s & %.0f & %.3f & $%.0f\\pm%.0f\\pm%.0f$ & %.3f & %.1f & %s & $%.0f\\pm%.0f$ & %.3f \\\\\n",
                    E[i].tex.c_str(), p.f.N, p.cbs, p.f.S, p.f.eS, p.spr.half, p.f.SoB, p.f.F, es.c_str(), b.T.S, b.T.eS, E[i].effT[0]);
            if (i == 0 || i == 7 || i == 10) fprintf(t, "\\midrule\n");
        }
        fprintf(t, "\\bottomrule\n\\end{tabular}\n\\end{table}\n");
        fclose(t);
    }
    // ---- LaTeX: variants (S, F for isBest and OA) ----
    {
        FILE* t = fopen((std::string(kTabOut) + "/opt_named_variants.tex").c_str(), "w");
        fprintf(t, "%% opt_named_variants.tex -- written by plots/exp_report/opt_named.C (numbers: tables_exp/opt_named.txt)\n");
        fprintf(t, "\\begin{table}[H]\n\\centering\\footnotesize\n\\setlength{\\tabcolsep}{3pt}\n");
        fprintf(t, "\\caption{$\\eta$ peak content $S\\pm$stat$\\pm$model spread (old envelope; revised: Table~\\ref{tab:fix_spread}) and $F=S/\\sigma_S^\\mathrm{stat}$ of the named selections in the four variants: all stored combinations "
                   "(isBest $\\geq 0$, main) or isBest $=1$ (comparison), without and with OA$(e^+e^-)>2^\\circ$ (OA $>2^\\circ$ is not a fixed cut). Last column: edge-step visibility $E$ of the bump "
                   "(isBest $\\geq0$ / with OA $>2^\\circ$); $^\\dagger$: summed edge fit fails ($p<10^{-3}$), $E$ not valid. Same fit and conventions as Table~\\ref{tab:opt_named}. Rows, not normalised to luminosity.}\n");
        fprintf(t, "\\label{tab:opt_named_variants}\n\\resizebox{\\textwidth}{!}{%%\n\\begin{tabular}{lcccccc}\n\\toprule\n");
        fprintf(t, " & \\multicolumn{2}{c}{isBest $\\geq 0$} & \\multicolumn{2}{c}{isBest $=1$} & \\\\\n");
        fprintf(t, "Selection & no OA cut & OA $>2^\\circ$ & no OA cut & OA $>2^\\circ$ & $E$ / $E$(OA$>2^\\circ$) \\\\\n\\midrule\n");
        for (int i = 0; i < ne; ++i) {
            fprintf(t, "%s", E[i].tex.c_str());
            const int ord[4] = {0, 1, 2, 3};
            for (int k = 0; k < 4; ++k) {
                const PeakRes& p = E[i].pk[ord[k]];
                fprintf(t, " & $%.0f\\pm%.0f\\pm%.0f$ (%.1f)", p.f.S, p.f.eS, p.spr.half, p.f.F);
            }
            auto ef = [&](const BumpRes& b) { return EdgeBad(b) ? std::string(Form("(%.0f)$^\\dagger$", b.E)) : std::string(Form("%.0f", b.E)); };
            fprintf(t, " & %s / %s \\\\\n", ef(E[i].bp[0]).c_str(), ef(E[i].bp[1]).c_str());
            if (i == 0 || i == 7 || i == 10 || i == 12) fprintf(t, "\\midrule\n");
        }
        fprintf(t, "\\bottomrule\n\\end{tabular}}\n\\end{table}\n");
        fclose(t);
    }

    // ========================================================================
    // figures
    // ========================================================================
    TFile* fo = TFile::Open(Form("%s/opt_named.root", kWork), "RECREATE");
    // ---- (1) M4 spectra of four selections, equal y ranges ----
    auto yMax = [](const Spectra& s, double a, double b) { double m = 0; for (int i = 0; i < s.nbin; ++i) if (s.x(i) > a && s.x(i) < b) m = std::max(m, s.n[0][i]); return m; };
    {
        const int idx[4] = {iPre, 3, 5, 6};
        const char* ttl[4] = {"(a) pre (= BUMP, frozen)", "(b) peak nominal + vertex", "(c) DIRECT (frozen): pre + #Deltaz_{#pi} < 40 mm",
                              "(d) PEAK (frozen): DIRECT + MM_{4} < 2.85 GeV/c^{2}"};
        double ymax = 0;
        for (int k = 0; k < 4; ++k) ymax = std::max(ymax, yMax(sp[idx[k] * (kNVar + 2)], kX1, kX2));
        ymax *= 1.06;
        TCanvas* c = new TCanvas("c_spec", "", 1400, 1200);
        c->Divide(2, 2, 0.001, 0.001);
        for (int k = 0; k < 4; ++k) {
            DrawSpectrum((TPad*)c->cd(k + 1), sp[idx[k] * (kNVar + 2)], E[idx[k]].pk[0], ttl[k], ymax, -450, 1150, k == 1, Form("s%d", k));
            fo->cd(); sp[idx[k] * (kNVar + 2)].Hist("us", Form("m4_us_%d", idx[k]))->Write();
            sp[idx[k] * (kNVar + 2)].Hist("cb", Form("m4_cb_%d", idx[k]))->Write();
        }
        Save(c, "opt_named_spectra");
    }
    // ---- (2) X17 band ----
    {
        TCanvas* c = new TCanvas("c_x17", "", 800, 700);
        const Spectra& s = sp[iX17 * (kNVar + 2)];
        DrawSpectrum((TPad*)c->cd(), s, E[iX17].pk[0], "X17 band (frozen): pre + vertex + 15.5 #leq M(e^{+}e^{-}) < 24.5 MeV/c^{2}",
                     yMax(s, kX1, kX2) * 1.75, -110, 230, false, "x17", true);
        c->cd(); TPad* pt = (TPad*)c->GetPrimitive("pt_x17"); pt->cd();
        // [review fix 2026-11-14] the recommended X17 M4 window (simulation window moved to the data scale) is drawn too
        const double ytop = yMax(s, kX1, kX2) * 1.75;
        VLine(kShPeakLo, 0, 0.42 * ytop, kOrange + 7, 2, 3); VLine(kShPeakHi, 0, 0.42 * ytop, kOrange + 7, 2, 3);
        TxtPx(0.19, 0.63, "hatched: frozen simulation window 0.5333-0.5587", 17, kOrange + 7);
        TxtPx(0.19, 0.57, Form("orange dashed: moved to the data scale %.4f-%.4f (recommended)", kShPeakLo, kShPeakHi), 17, kOrange + 7);
        TxtPx(0.19, 0.51, "grey: data peak window 0.5336-0.5542", 17, kGray + 2);
        TxtPx(0.19, 0.45, "black N, red CB, blue N - CB; magenta: Gaussian + pol2", 17, kBlack);
        Save(c, "opt_named_x17");
    }
    // ---- (3) tagged eta for the bump selections ----
    {
        const int idx[4] = {iPre, 8, 9, 6};
        const char* ttl[4] = {"(a) pre (= BUMP, frozen)", "(b) bump nominal", "(c) bump nominal + vertex", "(d) PEAK (frozen)"};
        std::vector<Sel> ts;
        for (int k = 0; k < 4; ++k) { const Cuts w = E[idx[k]].c; ts.push_back([=](const Row& r) { return Tagged(r) && Pass(w, r); }); }
        const std::vector<Spectra> st = scan::FillMany(dT, ts, TagOpt());
        const double tx1 = 0.44, tx2 = 0.68;
        double ymax = 0;
        for (int k = 0; k < 4; ++k) ymax = std::max(ymax, yMax(st[k], tx1, tx2));
        ymax *= 1.10;
        TCanvas* c = new TCanvas("c_tag", "", 1400, 1050);
        c->Divide(2, 2, 0.001, 0.001);
        for (int k = 0; k < 4; ++k) {
            TPad* p = (TPad*)c->cd(k + 1);
            p->SetLeftMargin(0.15); p->SetRightMargin(0.03); p->SetTopMargin(0.09); p->SetBottomMargin(0.16);
            const PeakFit f = scan::FitPeakMain(st[k], TagCfg());
            TH1D* fr = Frame(Form("fr_t%d", k), tx1, tx2, 0, ymax, "m_{5c} = M(#pi^{+}#pi^{-}e^{+}e^{-}#gamma) with M(e^{+}e^{-}#gamma) #rightarrow m_{#pi^{0}} [GeV/c^{2}]", "counts / 2.5 MeV/c^{2}");
            AxPx(fr, 20, 22, 1.6, 2.2);
            TH1D* hn = st[k].Hist("us", Form("tn%d", k)); TH1D* hc = st[k].Hist("cb", Form("tc%d", k)); TH1D* hs = st[k].Hist("sig", Form("ts%d", k));
            for (TH1D* h : {hn, hc, hs}) h->GetXaxis()->SetRangeUser(tx1, tx2 - 1e-6);
            hn->SetMarkerStyle(20); hn->SetMarkerSize(0.6); hn->SetMarkerColor(kBlack); hn->SetLineColor(kBlack);
            hc->SetMarkerStyle(21); hc->SetMarkerSize(0.55); hc->SetMarkerColor(kRed + 1); hc->SetLineColor(kRed + 1);
            hs->SetMarkerStyle(24); hs->SetMarkerSize(0.6); hs->SetMarkerColor(kBlue + 1); hs->SetLineColor(kBlue + 1);
            TGraph* gt = new TGraph(); TGraph* gb = new TGraph();
            for (int i = 0; i <= 200; ++i) { const double m = kTagLo + i * (kTagHi - kTagLo) / 200.; gt->SetPoint(i, m, f.Total(m)); gb->SetPoint(i, m, f.Bkg(m)); }
            gt->SetLineColor(kMagenta + 1); gt->SetLineWidth(3); gb->SetLineColor(kGray + 2); gb->SetLineStyle(2); gb->SetLineWidth(2);
            hn->Draw("PE SAME"); hc->Draw("PE SAME"); hs->Draw("PE SAME"); gb->Draw("L SAME"); gt->Draw("L SAME");
            TxtPx(0.15, 0.93, ttl[k], 22);
            TxtPx(0.19, 0.82, Form("T = %.0f #pm %.0f (stat),  #varepsilon_{T} = %.3f #pm %.3f", f.S, f.eS, E[idx[k]].effT[0], E[idx[k]].eEffT[0]), 19, kBlue + 2);
            if (k == 0) {
                TLegend* l = new TLegend(0.19, 0.48, 0.55, 0.78);
                l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextFont(43); l->SetTextSize(18);
                l->AddEntry(hn, "tagged counts N", "pe"); l->AddEntry(hc, "CB = 2#sqrt{N_{++}N_{--}}", "pe"); l->AddEntry(hs, "N - CB", "pe");
                l->AddEntry(gt, "fit: Gaussian + pol2", "l"); l->Draw();
            }
        }
        Save(c, "opt_named_tagged");
    }
    // ---- (4) summary chart: S, F, E and eff_T per selection and variant ----
    {
        std::vector<int> rows;
        for (int i = 0; i < ne; ++i) if (i != iX17) rows.push_back(i);
        const int nr = (int)rows.size();
        TCanvas* c = new TCanvas("c_sum", "", 1600, 900);
        const double xl = 0.25;
        const double dx = (1.0 - xl) / 4.0;
        const double xs[5] = {xl, xl + dx, xl + 2 * dx, xl + 3 * dx, 1.0};
        const char* ttl[4] = {"(a) #eta peak S [counts]", "(b) F = S/#sigma_{S}^{stat}", "(c) bump: edge step E", "(d) tagged #eta: #varepsilon_{T}"};
        const double lo[4] = {0, 7, 0, 0.45}, hi[4] = {6400, 16, 1750, 1.08};
        const int col[kNVar] = {kBlue + 1, kAzure + 7, kCyan + 3, kGreen + 2};   // [review fix] isBest == 1 no longer red (red = CB)
        const int mk[kNVar] = {20, 24, 21, 25};
        for (int q = 0; q < 4; ++q) {
            c->cd();
            TPad* p = new TPad(Form("p%d", q), "", q == 0 ? 0.0 : xs[q], 0.07, xs[q + 1], 0.95);
            const double w = (q == 0 ? xs[1] : xs[q + 1] - xs[q]);
            p->SetLeftMargin(q == 0 ? (xl / w) : 0.03); p->SetRightMargin(0.04); p->SetTopMargin(0.07); p->SetBottomMargin(0.08);
            p->Draw(); p->cd();
            TH1D* fr = new TH1D(Form("fs%d", q), "", 10, lo[q], hi[q]); fr->SetDirectory(nullptr);
            fr->SetMinimum(0); fr->SetMaximum(nr); fr->GetYaxis()->SetNdivisions(0); fr->GetYaxis()->SetLabelSize(0);
            fr->GetXaxis()->SetNdivisions(504);
            fr->GetXaxis()->SetLabelFont(43); fr->GetXaxis()->SetLabelSize(26);
            fr->Draw("AXIS");
            for (int r = 0; r < nr; ++r) {
                const int i = rows[r];
                const double y0 = nr - r - 0.5;
                if (E[i].kind == 1) { TBox* b = new TBox(lo[q], y0 - 0.5, hi[q], y0 + 0.5); b->SetFillColor(TColor::GetColor("#fdf0d2")); b->SetLineWidth(0); b->Draw(); }
                if (q == 2 && EdgeBad(E[i].bp[0])) {   // E not valid: hatched row + note
                    TBox* b = new TBox(lo[q], y0 - 0.5, hi[q], y0 + 0.5); b->SetFillColor(kGray + 1); b->SetFillStyle(3354);
                    b->SetLineWidth(0); b->Draw();
                    TLatex* t = new TLatex(lo[q] + 0.03 * (hi[q] - lo[q]), y0, "edge fit invalid");
                    t->SetTextAlign(12); t->SetTextSize(22); t->SetTextFont(43); t->SetTextColor(kGray + 2); t->Draw();
                }
                if (q == 0) {
                    TLatex* t = new TLatex(lo[q] - 0.03 * (hi[q] - lo[q]), y0, E[i].c.name.c_str());
                    t->SetTextAlign(32); t->SetTextSize(25); t->SetTextFont(43); t->Draw();
                }
            }
            if (q == 0) VLine(E[iPre].pk[0].f.S, 0, nr, kGray + 1, 2, 1);
            if (q == 1) VLine(E[iPre].pk[0].f.F, 0, nr, kGray + 1, 2, 1);
            if (q == 2) VLine(E[iPre].bp[0].E, 0, nr, kGray + 1, 2, 1);
            if (q == 3) VLine(1.0, 0, nr, kGray + 1, 2, 1);
            for (int v = 0; v < kNVar; ++v) {
                TGraphErrors* g = GE(col[v], mk[v], 1.1);
                TGraphErrors* gtot = GE(col[v], 1, 0.1); gtot->SetLineWidth(1);
                for (int r = 0; r < nr; ++r) {
                    const int i = rows[r];
                    const double y = nr - r - 0.5 + 0.30 - 0.20 * v;
                    double x = 0, e = 0, et = 0;
                    if (q == 0) { x = E[i].pk[v].f.S; e = E[i].pk[v].f.eS; et = std::hypot(e, E[i].pk[v].spr.half); }
                    if (q == 1) { x = E[i].pk[v].f.F; }
                    if (q == 2) { x = E[i].bp[v].E; e = E[i].bp[v].eE; et = std::hypot(e, E[i].bp[v].spreadE); }
                    if (q == 3) { x = E[i].effT[v]; e = E[i].eEffT[v]; et = e; }
                    g->SetPoint(r, x, y); g->SetPointError(r, e, 0);
                    gtot->SetPoint(r, x, y); gtot->SetPointError(r, et, 0);
                }
                if (q == 0 || q == 2) gtot->Draw("Z SAME");
                g->Draw("P SAME");
            }
            c->cd();
            TLatex* t = new TLatex(q == 0 ? xl + 0.5 * (xs[1] - xl) : 0.5 * (xs[q] + xs[q + 1]), 0.955, ttl[q]);
            t->SetNDC(); t->SetTextAlign(21); t->SetTextFont(43); t->SetTextSize(27); t->Draw();
        }
        c->cd();
        TLegend* l2 = new TLegend(0.27, 0.000, 0.99, 0.055);
        l2->SetBorderSize(0); l2->SetFillStyle(0); l2->SetTextFont(43); l2->SetTextSize(25);
        l2->SetNColumns(4);
        for (int v = 0; v < kNVar; ++v) { TGraphErrors* g = GE(col[v], mk[v], 1.1); l2->AddEntry(g, kVarName[v], "p"); }
        l2->Draw();
        TxtPx(0.005, 0.015, "thick: stat; thin: stat #oplus model", 22, kGray + 2, 11);
        Save(c, "opt_named_summary");
    }
    fo->Close();
    printf("opt_named: done in %.1f s\n", sw.RealTime());
}
