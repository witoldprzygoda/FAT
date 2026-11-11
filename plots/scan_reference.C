// ============================================================================
// scan_reference.C  --  frozen numbers, reference table and reference figures
// ============================================================================
// Uses scan_compact.root (plots/scan_prepare.C) through plots/ScanEngine.h.
//
//   A. how the frozen numbers of ScanEngine.h were obtained: free peak fits on
//      the reference selection (oa_lab < 50, cor_ok, isBest >= 0) as a function
//      of the fit range; K0s -> pi+pi- as a cross-check of the mass scale
//   B. reference table: N, CB, CB share, S +- stat, model spread, S/B, F for
//      the selections none, oa_lab<50, slides, old code, old code + OA 2-6;
//      isBest >= 0 and isBest == 1; fallback variant; like-sign null test
//   C. edge-step fit in three M(pi+pi-) slices (check of the engine)
//   D. figures in plots/output/eta_scan/ref_*.png (.pdf for the main ones)
//
// Output: <prod>/scan/reference_numbers.txt, <prod>/scan/scan_reference.root
//
// Run:
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/scan/aclic_engine", kTRUE)' 'plots/scan_reference.C+'
// ============================================================================
#include <cstdarg>
#include <cstdio>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TStyle.h"
#include "TSystem.h"

#include "ScanEngine.h"

namespace sref {

FILE* gOut = nullptr;
// print to the screen and to reference_numbers.txt
void P(const char* fmt, ...) {
    va_list a, b;
    va_start(a, fmt);
    va_copy(b, a);
    vprintf(fmt, a);
    if (gOut) vfprintf(gOut, fmt, b);
    va_end(a);
    va_end(b);
}

const char* kXt = "M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]";
const char* kYt = "rows / 2.5 MeV/c^{2}";

void Style() {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetPadLeftMargin(0.14);
    gStyle->SetPadRightMargin(0.04);
    gStyle->SetPadBottomMargin(0.14);
    gStyle->SetPadTopMargin(0.08);
    gStyle->SetLabelSize(0.05, "xy");
    gStyle->SetTitleSize(0.055, "xy");
    gStyle->SetTitleOffset(1.25, "y");
    gStyle->SetTitleOffset(1.10, "x");
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);
    gStyle->SetEndErrorSize(0);
}
void Deco(TH1* h, int col, int mk, double ms = 0.7) {
    h->SetLineColor(col); h->SetMarkerColor(col); h->SetMarkerStyle(mk); h->SetMarkerSize(ms);
    h->GetXaxis()->SetTitle(kXt); h->GetYaxis()->SetTitle(kYt);
    h->GetXaxis()->SetNdivisions(506);
}
void Head(const char* txt, double x = 0.14, double y = 0.94, double size = 0.042) {
    TLatex l; l.SetNDC(); l.SetTextSize(size); l.SetTextFont(42);
    l.DrawLatex(x, y, txt);
}
void VLine(double x, double y1, double y2, int col, int sty) {
    TLine* l = new TLine(x, y1, x, y2);
    l->SetLineColor(col); l->SetLineStyle(sty); l->SetLineWidth(2);
    l->Draw();
}
void Save(TCanvas* c, const char* name, bool pdf) {
    c->SaveAs(Form("%s/%s.png", scan::kFigDir, name));
    if (pdf) c->SaveAs(Form("%s/%s.pdf", scan::kFigDir, name));
}
// fitted curve of a peak fit, rows per fine bin
TGraph* Curve(const scan::PeakFit& f, bool total, int col, int sty) {
    TGraph* g = new TGraph();
    const int n = 400;
    for (int i = 0; i <= n; ++i) {
        const double m = f.cfg.lo + (f.cfg.hi - f.cfg.lo) * i / n;
        g->SetPoint(i, m, total ? f.Total(m) : f.Bkg(m));
    }
    g->SetLineColor(col); g->SetLineStyle(sty); g->SetLineWidth(2);
    return g;
}
// spectrum minus the fitted polynomial
TH1D* Resid(const scan::Spectra& sp, const std::vector<double>& y, const std::vector<double>& var,
            const scan::PeakFit& f, const char* name) {
    TH1D* h = new TH1D(name, "", sp.nbin, sp.lo, sp.hi);
    h->SetDirectory(nullptr);
    for (int i = 0; i < sp.nbin; ++i) {
        h->SetBinContent(i + 1, y[i] - f.Bkg(sp.x(i)));
        h->SetBinError(i + 1, std::sqrt(var[i]));
    }
    return h;
}
TGraph* GaussCurve(const scan::PeakFit& f, int col) {
    TGraph* g = new TGraph();
    const int n = 400;
    for (int i = 0; i <= n; ++i) {
        const double m = f.cfg.lo + (f.cfg.hi - f.cfg.lo) * i / n;
        g->SetPoint(i, m, f.Total(m) - f.Bkg(m));
    }
    g->SetLineColor(col); g->SetLineWidth(2);
    return g;
}

struct Res {
    std::string name;
    scan::Spectra sp;
    scan::PeakFit pm, pp, lsP, lsM;
    scan::Spread sd;
};
Res Analyse(const std::string& name, const scan::Spectra& sp) {
    Res r; r.name = name; r.sp = sp;
    r.pm = scan::FitPeakMain(sp);
    r.pp = scan::FitPeakPois(sp);
    r.lsP = scan::FitPeakLikeSign(sp, 1);
    r.lsM = scan::FitPeakLikeSign(sp, 2);
    r.sd = scan::PeakModelSpread(sp);
    return r;
}
void TableHead() {
    P("%-20s %-9s %8s %7s %6s | %6s %5s %6s %5s %5s %6s | %7s %5s %6s | %6s %5s %6s | %6s %5s | %6s %5s\n",
      "selection", "isBest", "N", "CB", "CB/N", "S", "stat", "spread", "S/B", "F", "chi2", "B", "prob", "rms",
      "S_pois", "stat", "chi2", "S(++)", "err", "S(--)", "err");
}
void TableRow(const Res& r, const char* best) {
    P("%-20s %-9s %8.0f %7.0f %6.3f | %6.0f %5.0f %6.0f %5.3f %5.1f %4.1f/%-2d| %7.0f %5.2f %6.0f | %6.0f %5.0f %4.1f/%-2d| %6.0f %5.0f | %6.0f %5.0f\n",
      r.name.c_str(), best, r.pm.N, r.pm.CB, r.pm.N > 0 ? r.pm.CB / r.pm.N : 0.,
      r.pm.S, r.pm.eS, r.sd.half, r.pm.SoB, r.pm.F, r.pm.chi2, r.pm.ndf,
      r.pm.B, r.pm.prob, r.sd.rms,
      r.pp.S, r.pp.eS, r.pp.chi2, r.pp.ndf, r.lsP.S, r.lsP.eS, r.lsM.S, r.lsM.eS);
}

} // namespace sref

void scan_reference() {
    using namespace sref;
    TStopwatch swAll;
    Style();
    gSystem->mkdir(scan::kFigDir, kTRUE);
    const std::string txt = std::string(scan::kScanDir) + "/reference_numbers.txt";
    gOut = fopen(txt.c_str(), "w");

    scan::Data d;
    if (!scan::Load(d)) return;

    // ------------------------------------------------------------------------
    // spectra: 5 selections x (all rows, isBest == 1, the rows added by isBest >= 0)
    // ------------------------------------------------------------------------
    const int nSel = 5;
    const char* selName[nSel] = {"none", "oa_lab<50", "slides", "old code", "old code + OA 2-6"};
    const char* selTag[nSel] = {"none", "oalab50", "slides", "oldcode", "oldcodeOA"};
    const scan::Sel base[nSel] = {scan::SelNone(), scan::SelOaLab50(), scan::SelSlides(),
                                  scan::SelOldCode(), scan::SelOldCodeOA()};
    std::vector<scan::Sel> sels;
    for (int i = 0; i < nSel; ++i) sels.push_back(base[i]);
    for (int i = 0; i < nSel; ++i) sels.push_back(scan::SelBest(base[i]));
    for (int i = 0; i < nSel; ++i) {
        const scan::Sel b = base[i];
        sels.push_back([b](const scan::Row& r) { return r.isBest != 1 && b(r); });
    }
    TStopwatch sw;
    const std::vector<scan::Spectra> spv = scan::FillMany(d, sels);
    const double tFill = sw.RealTime();
    // fallback variant: reconstructed pi- momentum for the cor_ok == 0 rows
    d.SetFallback(true);
    std::vector<scan::Sel> selsFb(sels.begin(), sels.begin() + nSel);
    const std::vector<scan::Spectra> spFb = scan::FillMany(d, selsFb);
    d.SetFallback(false);

    sw.Start();
    std::vector<Res> all, best, fb;
    for (int i = 0; i < nSel; ++i) all.push_back(Analyse(selName[i], spv[i]));
    for (int i = 0; i < nSel; ++i) best.push_back(Analyse(selName[i], spv[nSel + i]));
    for (int i = 0; i < nSel; ++i) fb.push_back(Analyse(selName[i], spFb[i]));
    const double tFit = sw.RealTime();
    const scan::Spectra& ref = spv[1];           // reference selection: oa_lab < 50

    // ------------------------------------------------------------------------
    // A. frozen numbers
    // ------------------------------------------------------------------------
    P("==========================================================================================\n");
    P("reference_numbers.txt  --  written by plots/scan_reference.C\n");
    P("data: %s\n", scan::kCompact);
    P("rows = stored particle combinations (isBest >= 0), corrected momenta, cor_ok == 1\n");
    P("==========================================================================================\n\n");
    P("A. FROZEN NUMBERS (plots/ScanEngine.h)\n");
    P("   binning            %.1f MeV, %.2f-%.2f GeV/c2 (%d bins); CB in %.0f MeV bins, linear interpolation\n",
      1e3 * scan::kBinW, scan::kMlo, scan::kMhi, scan::kNbin, 1e3 * scan::kBinW * scan::kCBreb);
    P("   peak position      %.4f GeV/c2\n", scan::kPeakMean);
    P("   peak width         %.1f MeV (Gaussian sigma)\n", 1e3 * scan::kPeakSigma);
    P("   fit range          %.4f-%.4f GeV/c2, polynomial order %d\n", scan::kFitLo, scan::kFitHi, scan::kFitOrder);
    P("   peak window        %.5f-%.5f GeV/c2 (+-%.1f sigma, holds %.4f of the Gaussian)\n", scan::kPeakLo, scan::kPeakHi,
      scan::kWinNsig, scan::GaussFrac(scan::kPeakLo, scan::kPeakHi, scan::kPeakMean, scan::kPeakSigma));
    P("   bump window        %.4f-%.4f GeV/c2\n", scan::kBumpLo, scan::kBumpHi);
    P("   mass scale         %.5f for M4 (peak position / %.6f), %.4f for M(pipi) (K0s, see A2)\n", scan::kScale, scan::kMeta, scan::kScalePiPi);
    P("   edge smearing      %.1f MeV\n", 1e3 * scan::kEdgeSigma);
    P("   bump edge, measured scale, for M(pipi) = 0.2791 (true), 0.30, 0.34, 0.37, 0.39, 0.405 (measured):\n     ");
    P(" %.4f", scan::kScale * scan::M4max(2 * scan::kMpi));
    for (double m : {0.30, 0.34, 0.37, 0.39, 0.405}) P(" %.4f", scan::EdgeMeasured(m));
    P("\n\n");

    P("A1. Free peak fit on the reference selection (oa_lab < 50), pol2, as a function of the range.\n");
    P("    main = chi2 fit of N - CB, pois = Poisson likelihood fit of N. 'frozen' = fit with the frozen shape.\n");
    P("    %-15s | %-44s | %-44s | %-27s | %-18s\n", "range", "main, free position and width",
      "pois, free position and width", "main, frozen shape", "pois, frozen shape");
    P("    %-15s | %8s %7s %6s %5s %6s %5s %4s    | %8s %7s %6s %5s %6s %5s %4s    | %6s %5s %9s %5s | %6s %5s %6s\n", "",
      "position", "err", "width", "err", "S", "err", "chi2", "position", "err", "width", "err", "S", "err", "chi2",
      "S", "err", "chi2/ndf", "prob", "S", "err", "chi2");
    struct RangePt { double lo, hi; scan::FreeFit fm, fp; scan::PeakFit pm, pp; };
    std::vector<RangePt> rlo, rhi;
    auto rangeLine = [&](double lo, double hi, std::vector<RangePt>& store) {
        scan::PeakCfg c; c.lo = lo; c.hi = hi;
        RangePt q; q.lo = lo; q.hi = hi;
        q.fm = scan::FitPeakFree(ref, c, false);
        q.fp = scan::FitPeakFree(ref, c, true);
        q.pm = scan::FitPeakMain(ref, c);
        q.pp = scan::FitPeakPois(ref, c);
        store.push_back(q);
        P("    %.4f-%.4f%s | %8.5f %7.5f %6.2f %5.2f %6.0f %5.0f %4.1f/%-2d | %8.5f %7.5f %6.2f %5.2f %6.0f %5.0f %4.1f/%-2d | %6.0f %5.0f %6.1f/%-2d %5.3f | %6.0f %5.0f %6.1f\n",
          lo, hi, (lo == scan::kFitLo && hi == scan::kFitHi) ? "*" : " ",
          q.fm.mean, q.fm.eMean, 1e3 * q.fm.sigma, 1e3 * q.fm.eSigma, q.fm.S, q.fm.eS, q.fm.chi2, q.fm.ndf,
          q.fp.mean, q.fp.eMean, 1e3 * q.fp.sigma, 1e3 * q.fp.eSigma, q.fp.S, q.fp.eS, q.fp.chi2, q.fp.ndf,
          q.pm.S, q.pm.eS, q.pm.chi2, q.pm.ndf, q.pm.prob, q.pp.S, q.pp.eS, q.pp.chi2);
    };
    for (double lo : {0.5050, 0.5100, 0.5150, 0.5175, 0.5200, 0.5225, 0.5250, 0.5275, 0.5300})
        rangeLine(lo, scan::kFitHi, rlo);
    P("    upper bound:\n");
    for (double hi : {0.5800, 0.5900, 0.6000, 0.6100, 0.6200, 0.6300, 0.6500}) rangeLine(scan::kFitLo, hi, rhi);
    P("    (* = frozen range; widths in MeV)\n");
    {
        const scan::FreeFit f0 = scan::FitPeakFree(spv[0], scan::PeakCfg(), false);
        P("    selection 'none', frozen range, main free: position %.5f +- %.5f, width %.2f +- %.2f MeV, S %.0f +- %.0f, chi2 %.1f/%d\n",
          f0.mean, f0.eMean, 1e3 * f0.sigma, 1e3 * f0.eSigma, f0.S, f0.eS, f0.chi2, f0.ndf);
        const scan::FreeFit f1 = scan::FitPeakFree(spv[nSel + 1], scan::PeakCfg(), false);
        P("    oa_lab<50, isBest == 1, frozen range, main free: position %.5f +- %.5f, width %.2f +- %.2f MeV, S %.0f +- %.0f, chi2 %.1f/%d\n",
          f1.mean, f1.eMean, 1e3 * f1.sigma, 1e3 * f1.eSigma, f1.S, f1.eS, f1.chi2, f1.ndf);
    }

    // ---- K0s -> pi+pi- in M(pi+pi-) as a cross-check of the mass scale ----
    // hK: pion pairs with a displaced crossing point; hK0: all rows with M4 > 0.60
    TH1D* hK = new TH1D("hK0s", "", 100, 0.42, 0.57);
    TH1D* hK0 = new TH1D("hK0s_incl", "", 100, 0.42, 0.57);
    hK->SetDirectory(nullptr); hK0->SetDirectory(nullptr);
    for (const auto& r : d.us()) {
        if (!r.cor_ok) continue;
        if (r.dz_pipi > 20.f && r.doca_pipi < 10.f) hK->Fill(r.mpipi);
        if (r.m4 > 0.60f) hK0->Fill(r.mpipi);
    }
    TF1* fK = new TF1("fK", "gausn(0)*0.0015+pol2(3)", 0.455, 0.54);
    fK->SetParameters(0.03 * hK->Integral(), 0.496, 0.006, hK->GetBinContent(hK->FindBin(0.47)), 0, 0);
    fK->SetParLimits(1, 0.485, 0.505); fK->SetParLimits(2, 0.003, 0.012);
    TFitResultPtr rK = hK->Fit(fK, "RQSN0L");
    const double kPos = rK->Parameter(1), kPosE = rK->ParError(1);
    TF1* fK0 = new TF1("fK0", "gausn(0)*0.0015+pol2(3)", 0.46, 0.535);
    fK0->SetParameters(15000, 0.495, 0.0065, hK0->GetBinContent(hK0->FindBin(0.495)), 0, 0);
    fK0->SetParLimits(1, 0.485, 0.505); fK0->SetParLimits(2, 0.003, 0.012);
    TFitResultPtr rK0 = hK0->Fit(fK0, "RQSN0L");
    P("\nA2. K0s -> pi+pi- in M(pi+pi-) (isBest >= 0: one pion pair can enter several times, errors slightly too small).\n");
    P("    pion crossing point displaced (dz_pipi > 20 mm, doca_pipi < 10 mm): position %.5f +- %.5f, width %.2f +- %.2f MeV,\n",
      kPos, kPosE, 1e3 * rK->Parameter(2), 1e3 * rK->ParError(2));
    P("    content %.0f +- %.0f, chi2/ndf %.1f/%d. Scale of M(pipi): %.5f +- %.5f (frozen: %.4f).\n",
      rK->Parameter(0), rK->ParError(0), rK->Chi2(), rK->Ndf(), kPos / scan::kMK0, kPosE / scan::kMK0, scan::kScalePiPi);
    P("    all rows with M4 > 0.60 (peak on a kink of the background, less reliable): position %.5f +- %.5f, width %.2f +- %.2f MeV, scale %.5f.\n",
      rK0->Parameter(1), rK0->ParError(1), 1e3 * rK0->Parameter(2), 1e3 * rK0->ParError(2), rK0->Parameter(1) / scan::kMK0);
    P("    Scale of M4 from the eta peak (free fit, frozen range): %.5f +- %.5f.\n",
      rlo[5].fm.mean / scan::kMeta, rlo[5].fm.eMean / scan::kMeta);

    // ------------------------------------------------------------------------
    // B. reference table
    // ------------------------------------------------------------------------
    P("\nB. REFERENCE TABLE. Peak window %.5f-%.5f. N = unlike-sign rows, CB = 2 sqrt(N++ N--) (interpolated),\n", scan::kPeakLo, scan::kPeakHi);
    P("   S = content of the whole Gaussian from P_main (chi2 fit of N - CB, frozen shape, pol%d, %.4f-%.4f), stat = its fit error,\n",
      scan::kFitOrder, scan::kFitLo, scan::kFitHi);
    P("   spread = half of (max - min) of S over the model variations with p > 0.05 (orders 1-4, four ranges, width +-0.5 MeV; 48 fits),\n");
    P("   rms = their rms, B = fitted polynomial in the window, S/B = S in the window / (B + CB), F = S/stat,\n");
    P("   S_pois = Poisson fit of the unsubtracted N, S(++), S(--) = the same Poisson fit on the like-sign spectra (null test).\n\n");
    TableHead();
    for (int i = 0; i < nSel; ++i) TableRow(all[i], ">= 0");
    for (int i = 0; i < nSel; ++i) TableRow(best[i], "== 1");
    P("fallback variant (reconstructed pi- momentum where the corrected one is missing; cor_ok not required):\n");
    for (int i = 0; i < nSel; ++i) TableRow(fb[i], ">= 0, fb");
    P("\n   window sums with the event-group errors (sqrt of the sum over groups of k^2), isBest >= 0:\n");
    P("   %-20s | %-39s | %-48s\n", "selection", "peak window: N, error, sqrt(N), CB, error", "bump window 0.4550-0.5225: N, error, CB, error, N-CB, error");
    for (int i = 0; i < nSel; ++i) {
        const scan::WinSum& a = all[i].sp.win[0];
        const scan::WinSum& b = all[i].sp.win[1];
        P("   %-20s | %8.0f %6.1f %6.1f %7.1f %5.1f      | %8.0f %6.1f %7.1f %5.1f %8.0f %6.1f\n", selName[i],
          a.n[0], std::sqrt(a.var[0]), std::sqrt(a.n[0]), a.cb, std::sqrt(a.cbVar),
          b.n[0], std::sqrt(b.var[0]), b.cb, std::sqrt(b.cbVar), b.sig(), b.sigErr());
    }
    P("\n   split of the model spread, isBest >= 0 (fits kept of 48; half-range with the nominal width only; effect of the width alone):\n");
    for (int i = 0; i < nSel; ++i)
        P("   %-20s kept %2d, S from %5.0f to %5.0f, mean %5.0f; background model only (kept %2d of 16): %4.0f; width +-0.5 MeV only: %4.0f\n",
          selName[i], all[i].sd.npass, all[i].sd.smin, all[i].sd.smax, all[i].sd.mean, all[i].sd.npassBkg,
          all[i].sd.halfBkg, all[i].sd.halfWidth);

    // ---- the first-look configuration of the old production, on the new files ----
    {
        scan::PeakCfg c; c.mean = 0.5442; c.sigma = 0.0052; c.lo = 0.515; c.hi = 0.610; c.order = 2;
        c.winLo = 0.530; c.winHi = 0.550;
        const scan::PeakFit a = scan::FitPeakPois(spv[nSel + 1], c), b = scan::FitPeakMain(spv[nSel + 1], c);
        const scan::PeakFit a0 = scan::FitPeakPois(spv[1], c);
        P("\n   check against the old production: its first-look configuration (position 0.5442, width 5.2 MeV, range 0.515-0.610, pol2,\n");
        P("   Poisson fit of N) on oa_lab<50, isBest == 1 gives S = %.0f +- %.0f (old: 3554 +- 281; chi2 fit of N - CB: %.0f +- %.0f, old 3681 +- 304);\n",
          a.S, a.eS, b.S, b.eS);
        P("   the same configuration with isBest >= 0: %.0f +- %.0f.\n", a0.S, a0.eS);
    }

    // ---- isBest >= 0 against isBest == 1 ----
    P("\nB2. isBest >= 0 against isBest == 1 (P_main). The error of the ratio uses the fit of the rows added by isBest >= 0\n");
    P("    (Poisson fit, S_add), which are independent of the isBest == 1 rows.\n");
    P("    %-20s %8s %8s %6s | %6s %5s | %6s %5s | %6s %5s | %7s %6s | %5s %5s %6s\n", "selection", "N(>=0)", "N(==1)", "ratio",
      "S(>=0)", "stat", "S(==1)", "stat", "S_add", "stat", "S ratio", "err", "F>=0", "F==1", "ratio");
    for (int i = 0; i < nSel; ++i) {
        const scan::PeakFit add = scan::FitPeakPois(spv[2 * nSel + i]);
        const double s1 = best[i].pm.S, e1 = best[i].pm.eS;
        const double R = all[i].pm.S / s1;
        const double eR = std::sqrt(std::pow(add.eS / s1, 2) + std::pow(add.S * e1 / (s1 * s1), 2));
        P("    %-20s %8.0f %8.0f %6.3f | %6.0f %5.0f | %6.0f %5.0f | %6.0f %5.0f | %7.3f %6.3f | %5.2f %5.2f %6.3f\n",
          selName[i], all[i].pm.N, best[i].pm.N, all[i].pm.N / best[i].pm.N, all[i].pm.S, all[i].pm.eS, s1, e1,
          add.S, add.eS, R, eR, all[i].pm.F, best[i].pm.F, all[i].pm.F / best[i].pm.F);
    }

    // ------------------------------------------------------------------------
    // C. edge-step fit in M(pi+pi-) slices (oa_lab < 50)
    // ------------------------------------------------------------------------
    const int nSl = 4;
    const double slLo[nSl] = {0.30, 0.34, 0.37, 0.42}, slHi[nSl] = {0.34, 0.37, 0.39, 0.46};
    std::vector<scan::Sel> slSel;
    for (int i = 0; i < nSl; ++i) {
        const float a = slLo[i], b = slHi[i];
        slSel.push_back([a, b](const scan::Row& r) { return r.oa_lab < 50.f && r.mpipi >= a && r.mpipi < b; });
    }
    const std::vector<scan::Spectra> slSp = scan::FillMany(d, slSel);
    std::vector<scan::StepFit> slFit;
    P("\nC. EDGE-STEP FIT E_step (oa_lab < 50; polynomial order 2 + eta Gaussian + smeared step at the FIXED predicted edge).\n");
    P("   h = step height in rows per 2.5 MeV. The last line is a control slice above the kinematic limit, edge set by hand.\n");
    P("   %-12s %7s %15s %15s | %6s %5s %6s | %9s %7s %6s | %6s %5s\n", "M(pipi)", "edge", "edge range", "fit range",
      "h", "err", "bkg", "chi2/ndf", "no step", "dchi2", "S_eta", "err");
    for (int i = 0; i < nSl; ++i) {
        scan::StepCfg c; c.mpLo = slLo[i]; c.mpHi = slHi[i];
        if (i == 3) { c.useEdge = true; c.edge = 0.5150; }
        const scan::StepFit f = scan::FitEdgeStep(slSp[i], c);
        slFit.push_back(f);
        P("   %.2f-%.2f    %7.4f %7.4f-%7.4f %7.4f-%7.4f | %6.1f %5.1f %6.0f | %6.1f/%-2d %7.1f %6.1f | %6.0f %5.0f\n", slLo[i], slHi[i],
          f.edge, f.edgeLo, f.edgeHi, f.lo, f.hi, f.h, f.eh, f.bkgAtEdge, f.chi2, f.ndf, f.chi2NoStep, f.dchi2, f.Seta, f.eSeta);
    }

    P("\nrun times: load see above; FillMany of %zu selections %.2f s; 15 x (P_main + P_pois + 2 like-sign fits + 48 model fits) %.3f s\n",
      sels.size(), tFill, tFit);

    // ------------------------------------------------------------------------
    // ROOT file with the reference spectra
    // ------------------------------------------------------------------------
    {
        TFile fo((std::string(scan::kScanDir) + "/scan_reference.root").c_str(), "RECREATE");
        for (int v = 0; v < 3; ++v)
            for (int i = 0; i < nSel; ++i) {
                const scan::Spectra& sp = v == 0 ? spv[i] : (v == 1 ? spv[nSel + i] : spFb[i]);
                const char* vt = v == 0 ? "all" : (v == 1 ? "best" : "fb");
                for (const char* w : {"us", "pp", "mm", "cb", "sig"}) {
                    TH1D* h = sp.Hist(w, Form("%s_%s_%s", w, selTag[i], vt));
                    h->SetTitle(Form("%s, %s, %s;%s;%s", w, selName[i], vt, kXt, kYt));
                    h->Write();
                    delete h;
                }
            }
        hK->Write(); hK0->Write();
        fo.Close();
    }

    // ------------------------------------------------------------------------
    // D. figures
    // ------------------------------------------------------------------------
    const char* hdr = "p+p 4.53 GeV, all combinations (isBest #geq 0), corrected momenta";

    // ---- D1. overview spectra: log and linear zoom ----
    for (int i = 0; i < 2; ++i) {
        const scan::Spectra& sp = spv[i];
        TH1D* hn = sp.Hist("us", Form("ov_us_%d", i));
        TH1D* hc = sp.Hist("cb", Form("ov_cb_%d", i));
        TH1D* hs = sp.Hist("sig", Form("ov_sig_%d", i));
        Deco(hn, kBlack, 20); Deco(hc, kRed + 1, 20); Deco(hs, kBlue + 1, 20);
        TCanvas* c = new TCanvas(Form("c_ov_%d", i), "", 2000, 900);
        c->Divide(2, 1);
        c->cd(1); gPad->SetLogy(); gPad->SetLeftMargin(0.16);
        hn->GetYaxis()->SetTitleOffset(1.35);
        hn->SetMinimum(30); hn->SetMaximum(3.0 * hn->GetMaximum());
        hn->Draw("e"); hc->Draw("e same"); hs->Draw("e same");
        Head(hdr, 0.16, 0.94, 0.034);
        TLegend* lg = new TLegend(0.42, 0.18, 0.93, 0.42);
        lg->SetTextSize(0.042);
        lg->SetHeader(Form("selection: %s", i == 0 ? "no cut" : "oa_lab((#pi#pi),(ee)) < 50#circ"));
        lg->AddEntry(hn, "all pairs e^{+}e^{-} (N)", "pe");
        lg->AddEntry(hc, "CB = 2#sqrt{N_{++}N_{--}}", "pe");
        lg->AddEntry(hs, "N #minus CB", "pe");
        lg->Draw();
        c->cd(2); gPad->SetLeftMargin(0.21);
        TH1D* zn = (TH1D*)hn->Clone(Form("ovz_us_%d", i));
        zn->GetXaxis()->SetRangeUser(0.4001, 0.6999);
        zn->GetYaxis()->SetTitleOffset(1.95);
        zn->SetMinimum(0); zn->SetMaximum(1.12 * zn->GetBinContent(zn->FindBin(0.699)));
        zn->Draw("e"); hc->Draw("e same"); hs->Draw("e same");
        const double ym = zn->GetMaximum();
        VLine(scan::kPeakLo, 0, 0.62 * ym, kGray + 2, 2); VLine(scan::kPeakHi, 0, 0.62 * ym, kGray + 2, 2);
        VLine(scan::kBumpLo, 0, 0.42 * ym, kGreen + 2, 7); VLine(scan::kBumpHi, 0, 0.42 * ym, kGreen + 2, 7);
        Head(Form("the same, linear zoom (%s)", i == 0 ? "no cut" : "oa_lab < 50#circ"), 0.21, 0.94, 0.034);
        TLine* lp = new TLine(); lp->SetLineColor(kGray + 2); lp->SetLineStyle(2); lp->SetLineWidth(2);
        TLine* lb = new TLine(); lb->SetLineColor(kGreen + 2); lb->SetLineStyle(7); lb->SetLineWidth(2);
        TLegend* l2 = new TLegend(0.23, 0.66, 0.76, 0.89);
        l2->SetTextSize(0.038);
        l2->AddEntry(lp, Form("#eta peak window %.4f#minus%.4f", scan::kPeakLo, scan::kPeakHi), "l");
        l2->AddEntry(lb, Form("bump window %.3f#minus%.4f", scan::kBumpLo, scan::kBumpHi), "l");
        l2->AddEntry((TObject*)nullptr, "(bump: #eta#rightarrow#pi^{+}#pi^{-}#pi^{0}, #pi^{0}#rightarrowe^{+}e^{-}#gamma)", "");
        l2->Draw();
        Save(c, Form("ref_overview_%s", selTag[i]), true);
    }

    // ---- D2. reference peak fit ----
    {
        const Res& R = all[1];
        TCanvas* c = new TCanvas("c_fit", "", 1300, 1200);
        TPad* p1 = new TPad("p1", "", 0, 0.42, 1, 1); TPad* p2 = new TPad("p2", "", 0, 0, 1, 0.42);
        p1->SetBottomMargin(0.02); p2->SetTopMargin(0.03); p2->SetBottomMargin(0.24);
        p1->Draw(); p2->Draw();
        p1->cd();
        TH1D* hs = R.sp.Hist("sig", "fit_sig"); Deco(hs, kBlue + 1, 20, 1.0);
        hs->GetXaxis()->SetRangeUser(0.4951, 0.6349);
        hs->SetMinimum(0.80 * hs->GetBinContent(hs->FindBin(0.4975)));
        hs->SetMaximum(1.04 * hs->GetBinContent(hs->FindBin(0.634)));
        hs->GetXaxis()->SetLabelSize(0);
        hs->Draw("e");
        TGraph* gt = Curve(R.pm, true, kBlack, 1); TGraph* gb = Curve(R.pm, false, kBlack, 2);
        gb->Draw("l same"); gt->Draw("l same");
        const double y0 = hs->GetMinimum(), y1 = hs->GetMaximum();
        VLine(scan::kFitLo, y0, y0 + 0.35 * (y1 - y0), kGray + 2, 3); VLine(scan::kFitHi, y0, y0 + 0.9 * (y1 - y0), kGray + 2, 3);
        VLine(scan::kScale * scan::M4max(2 * scan::kMpi), y0, y0 + 0.25 * (y1 - y0), kGreen + 2, 7);
        Head(hdr, 0.14, 0.94, 0.040);
        TLegend* lg = new TLegend(0.16, 0.62, 0.62, 0.90);
        lg->SetTextSize(0.043);
        lg->SetHeader("oa_lab < 50#circ: reference peak fit P_main");
        lg->AddEntry(hs, "N #minus CB", "pe");
        lg->AddEntry(gt, Form("Gaussian (%.4f, #sigma = %.1f MeV, frozen) + pol2", scan::kPeakMean, 1e3 * scan::kPeakSigma), "l");
        lg->AddEntry(gb, "pol2 alone", "l");
        lg->Draw();
        TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.043);
        t.DrawLatex(0.56, 0.34, Form("S = %.0f #pm %.0f (stat) #pm %.0f (model)", R.pm.S, R.pm.eS, R.sd.half));
        t.DrawLatex(0.56, 0.27, Form("#chi^{2}/ndf = %.1f/%d, fit range %.4f#minus%.4f", R.pm.chi2, R.pm.ndf, scan::kFitLo, scan::kFitHi));
        t.DrawLatex(0.56, 0.20, Form("window #pm1.9#sigma: N = %.0f, CB = %.0f", R.pm.N, R.pm.CB));
        t.DrawLatex(0.56, 0.13, Form("S/B = %.3f, S/#sigma_{S} = %.1f", R.pm.SoB, R.pm.F));
        t.DrawLatex(0.56, 0.06, Form("cross-check P_pois: S = %.0f #pm %.0f", R.pp.S, R.pp.eS));
        t.SetTextColor(kGreen + 2); t.SetTextSize(0.036);
        t.DrawLatex(0.17, 0.07, "green: highest bump edge");
        p2->cd();
        TH1D* hr = Resid(R.sp, R.sp.sig, R.sp.sigVar, R.pm, "fit_res"); Deco(hr, kBlue + 1, 20, 1.0);
        hr->GetXaxis()->SetRangeUser(0.4951, 0.6349);
        hr->GetYaxis()->SetTitle("N #minus CB #minus pol2");
        hr->GetXaxis()->SetLabelSize(0.07); hr->GetXaxis()->SetTitleSize(0.08);
        hr->GetYaxis()->SetLabelSize(0.065); hr->GetYaxis()->SetTitleSize(0.075); hr->GetYaxis()->SetTitleOffset(0.9);
        hr->SetMinimum(-450); hr->SetMaximum(950);
        hr->Draw("e");
        TLine* l0 = new TLine(0.495, 0, 0.635, 0); l0->SetLineStyle(2); l0->Draw();
        GaussCurve(R.pm, kBlack)->Draw("l same");
        VLine(scan::kPeakLo, -450, 950, kGray + 2, 2); VLine(scan::kPeakHi, -450, 950, kGray + 2, 2);
        VLine(scan::kFitLo, -450, 300, kGray + 2, 3); VLine(scan::kFitHi, -450, 300, kGray + 2, 3);
        TLatex t2; t2.SetNDC(); t2.SetTextFont(42); t2.SetTextSize(0.055);
        t2.DrawLatex(0.56, 0.88, "dashed: peak window, dotted: fit range");
        Save(c, "ref_peakfit_oalab50", true);
    }

    // ---- D3. isBest >= 0 against isBest == 1 ----
    {
        const Res& A = all[1]; const Res& B = best[1];
        TCanvas* c = new TCanvas("c_best", "", 2000, 1300);
        c->Divide(2, 2);
        c->cd(1);
        TH1D* ha = A.sp.Hist("sig", "b_sig_all"); TH1D* hb = B.sp.Hist("sig", "b_sig_best");
        Deco(ha, kBlue + 1, 20, 0.9); Deco(hb, kBlue + 1, 24, 0.9);
        ha->GetXaxis()->SetRangeUser(0.4951, 0.6349);
        ha->SetMinimum(0.75 * hb->GetBinContent(hb->FindBin(0.4975)));
        ha->SetMaximum(1.05 * ha->GetBinContent(ha->FindBin(0.634)));
        ha->Draw("e"); hb->Draw("e same");
        Curve(A.pm, true, kBlack, 1)->Draw("l same"); Curve(A.pm, false, kBlack, 2)->Draw("l same");
        Curve(B.pm, true, kBlack, 1)->Draw("l same"); Curve(B.pm, false, kBlack, 2)->Draw("l same");
        Head("N #minus CB, oa_lab < 50#circ, with the reference fit P_main", 0.14, 0.94, 0.040);
        TLegend* lg = new TLegend(0.16, 0.70, 0.60, 0.89); lg->SetTextSize(0.042);
        lg->AddEntry(ha, "all combinations (isBest #geq 0)", "pe");
        lg->AddEntry(hb, "best combination only (isBest = 1)", "pe");
        lg->Draw();
        c->cd(2);
        TH1D* ra = Resid(A.sp, A.sp.sig, A.sp.sigVar, A.pm, "b_res_all");
        TH1D* rb = Resid(B.sp, B.sp.sig, B.sp.sigVar, B.pm, "b_res_best");
        Deco(ra, kBlue + 1, 20, 0.9); Deco(rb, kBlue + 1, 24, 0.9);
        ra->GetXaxis()->SetRangeUser(0.5101, 0.6149);
        ra->GetYaxis()->SetTitle("N #minus CB #minus pol2  [rows / 2.5 MeV/c^{2}]");
        ra->SetMinimum(-450); ra->SetMaximum(1550);
        ra->Draw("e"); rb->Draw("e same");
        GaussCurve(A.pm, kBlack)->Draw("l same");
        TGraph* gB = GaussCurve(B.pm, kBlack); gB->SetLineStyle(2); gB->Draw("l same");
        TLine* l0 = new TLine(0.51, 0, 0.615, 0); l0->SetLineStyle(3); l0->Draw();
        Head("the same after subtracting the fitted pol2", 0.14, 0.94, 0.040);
        TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.042);
        t.DrawLatex(0.30, 0.86, Form("isBest #geq 0: S = %.0f #pm %.0f, S/#sigma_{S} = %.1f", A.pm.S, A.pm.eS, A.pm.F));
        t.DrawLatex(0.30, 0.80, Form("isBest = 1: S = %.0f #pm %.0f, S/#sigma_{S} = %.1f", B.pm.S, B.pm.eS, B.pm.F));
        t.DrawLatex(0.30, 0.74, Form("S/B: %.3f and %.3f", A.pm.SoB, B.pm.SoB));
        c->cd(3);
        // ratio of the row counts, all / best, in 10 MeV bins
        TH1D* hr = new TH1D("b_ratio", "", 60, 0.35, 0.95); hr->SetDirectory(nullptr);
        TH1D* hrc = new TH1D("b_ratio_cb", "", 60, 0.35, 0.95); hrc->SetDirectory(nullptr);
        const scan::Spectra& sa = spv[1]; const scan::Spectra& sb = spv[nSel + 1]; const scan::Spectra& sx = spv[2 * nSel + 1];
        for (int j = 0; j < 60; ++j) {
            const double nb = sb.w[0][j], nx = sx.w[0][j];
            if (nb > 0 && nx > 0) { hr->SetBinContent(j + 1, 1 + nx / nb); hr->SetBinError(j + 1, nx / nb * std::sqrt(1 / nx + 1 / nb)); }
            if (sb.cbW[j] > 0) {
                hrc->SetBinContent(j + 1, sa.cbW[j] / sb.cbW[j]);
                hrc->SetBinError(j + 1, sa.cbW[j] / sb.cbW[j] * std::sqrt(sx.cbWVar[j]) / std::max(sa.cbW[j], 1.));
            }
        }
        Deco(hr, kBlack, 20, 0.9); Deco(hrc, kRed + 1, 20, 0.9);
        hr->GetYaxis()->SetTitle("rows(isBest #geq 0) / rows(isBest = 1)");
        hr->SetMinimum(1.0); hr->SetMaximum(1.6);
        hr->Draw("e"); hrc->Draw("e same");
        VLine(scan::kPeakLo, 1.0, 1.45, kGray + 2, 2); VLine(scan::kPeakHi, 1.0, 1.45, kGray + 2, 2);
        Head("oa_lab < 50#circ: gain in rows from using all combinations, 10 MeV bins", 0.14, 0.94, 0.040);
        TLegend* l3 = new TLegend(0.45, 0.72, 0.93, 0.89); l3->SetTextSize(0.042);
        l3->AddEntry(hr, "all pairs e^{+}e^{-} (N)", "pe");
        l3->AddEntry(hrc, "CB", "pe");
        l3->Draw();
        c->cd(4);
        TLatex q; q.SetNDC(); q.SetTextFont(42); q.SetTextSize(0.040);
        q.DrawLatex(0.05, 0.92, "Peak content S (whole Gaussian) and S/#sigma_{S}, fit P_main");
        q.SetTextSize(0.036);
        q.DrawLatex(0.05, 0.84, "selection                    S (isBest #geq 0)      S (isBest = 1)      ratio of S       S/#sigma_{S}: #geq 0 / = 1");
        for (int i = 0; i < nSel; ++i) {
            const scan::PeakFit add = scan::FitPeakPois(spv[2 * nSel + i]);
            const double s1 = best[i].pm.S, e1 = best[i].pm.eS;
            const double Rr = all[i].pm.S / s1;
            const double eR = std::sqrt(std::pow(add.eS / s1, 2) + std::pow(add.S * e1 / (s1 * s1), 2));
            q.DrawLatex(0.05, 0.76 - 0.08 * i, selName[i]);
            q.DrawLatex(0.30, 0.76 - 0.08 * i, Form("%.0f #pm %.0f", all[i].pm.S, all[i].pm.eS));
            q.DrawLatex(0.48, 0.76 - 0.08 * i, Form("%.0f #pm %.0f", s1, e1));
            q.DrawLatex(0.65, 0.76 - 0.08 * i, Form("%.2f #pm %.2f", Rr, eR));
            q.DrawLatex(0.82, 0.76 - 0.08 * i, Form("%.1f / %.1f", all[i].pm.F, best[i].pm.F));
        }
        q.DrawLatex(0.05, 0.30, "The two samples overlap: the error of the ratio comes from the");
        q.DrawLatex(0.05, 0.24, "fit of the rows that isBest #geq 0 adds (independent of isBest = 1).");
        q.DrawLatex(0.05, 0.16, "slides: oa_lab<50, M(#pi#pi)<0.46, eta_rest>140, OA>2");
        q.DrawLatex(0.05, 0.10, "old code: oa_lab<50, M(#pi#pi)<0.42, eta_rest>140 (and 2<OA<6)");
        Save(c, "ref_isbest_compare", true);
    }

    // ---- D4. range study ----
    {
        TCanvas* c = new TCanvas("c_range", "", 1300, 1500);
        c->Divide(1, 3);
        const int n = (int)rlo.size();
        TGraphErrors* gS = new TGraphErrors(); TGraphErrors* gSp = new TGraphErrors();
        TGraph* gC = new TGraph(); TGraphErrors* gW = new TGraphErrors(); TGraphErrors* gWp = new TGraphErrors();
        for (int i = 0; i < n; ++i) {
            gS->SetPoint(i, rlo[i].lo - 0.0002, rlo[i].pm.S); gS->SetPointError(i, 0, rlo[i].pm.eS);
            gSp->SetPoint(i, rlo[i].lo + 0.0002, rlo[i].pp.S); gSp->SetPointError(i, 0, rlo[i].pp.eS);
            gC->SetPoint(i, rlo[i].lo, rlo[i].pm.chi2 / rlo[i].pm.ndf);
            gW->SetPoint(i, rlo[i].lo - 0.0002, 1e3 * rlo[i].fm.sigma); gW->SetPointError(i, 0, 1e3 * rlo[i].fm.eSigma);
            gWp->SetPoint(i, rlo[i].lo + 0.0002, 1e3 * rlo[i].fp.sigma); gWp->SetPointError(i, 0, 1e3 * rlo[i].fp.eSigma);
        }
        auto frame = [&](int pad, double y0, double y1, const char* yt) {
            c->cd(pad); gPad->SetGridy();
            TH1D* f = new TH1D(Form("fr%d", pad), "", 10, 0.5030, 0.5320); f->SetDirectory(nullptr);
            f->SetMinimum(y0); f->SetMaximum(y1);
            f->GetXaxis()->SetTitle("lower bound of the fit range [GeV/c^{2}]  (upper bound 0.610, pol2)");
            f->GetYaxis()->SetTitle(yt); f->GetYaxis()->SetTitleOffset(0.9);
            f->Draw();
            VLine(scan::kFitLo, y0, y1, kGray + 2, 2);
            VLine(scan::kScale * scan::M4max(2 * scan::kMpi), y0, y1, kGreen + 2, 7);
        };
        frame(1, 3000, 5200, "peak content S, frozen shape");
        gS->SetMarkerStyle(20); gS->SetMarkerColor(kBlue + 1); gS->SetLineColor(kBlue + 1); gS->SetMarkerSize(1.3);
        gSp->SetMarkerStyle(24); gSp->SetMarkerColor(kBlack); gSp->SetMarkerSize(1.3);
        gS->Draw("p same"); gSp->Draw("p same");
        Head("oa_lab < 50#circ, isBest #geq 0: peak fit against the lower bound of the fit range", 0.14, 0.94, 0.050);
        TLegend* lg = new TLegend(0.16, 0.16, 0.62, 0.36); lg->SetTextSize(0.050);
        lg->AddEntry(gS, "P_main: #chi^{2} fit of N #minus CB", "pe");
        lg->AddEntry(gSp, "P_pois: likelihood fit of N", "pe");
        lg->Draw();
        frame(2, 0, 3.2, "#chi^{2}/ndf of P_main, frozen shape");
        gC->SetMarkerStyle(20); gC->SetMarkerColor(kBlue + 1); gC->SetMarkerSize(1.3); gC->Draw("p same");
        Head("dashed grey: chosen bound 0.5225; green: highest bump edge (measured scale)", 0.14, 0.94, 0.050);
        frame(3, 3.5, 7.5, "peak width from the free fit [MeV]");
        gW->SetMarkerStyle(20); gW->SetMarkerColor(kBlue + 1); gW->SetLineColor(kBlue + 1); gW->SetMarkerSize(1.3);
        gWp->SetMarkerStyle(24); gWp->SetMarkerColor(kBlack); gWp->SetMarkerSize(1.3);
        gW->Draw("p same"); gWp->Draw("p same");
        TLine* lw = new TLine(0.503, 1e3 * scan::kPeakSigma, 0.532, 1e3 * scan::kPeakSigma); lw->SetLineStyle(2); lw->Draw();
        Head(Form("free position and width; dashed line: frozen width %.1f MeV", 1e3 * scan::kPeakSigma), 0.14, 0.94, 0.050);
        Save(c, "ref_range_study", false);
    }

    // ---- D5. model variations ----
    {
        const scan::Spread sd = scan::PeakModelSpread(ref, scan::PeakCfg(), true);
        TCanvas* c = new TCanvas("c_var", "", 1700, 900);
        gPad->SetBottomMargin(0.16);
        TH1D* f = new TH1D("fr_var", "", 48, 0, 48); f->SetDirectory(nullptr);
        f->SetMinimum(2200); f->SetMaximum(6400);
        f->GetYaxis()->SetTitle("peak content S");
        f->GetXaxis()->SetTitle("model variation: range (blocks of 12), polynomial order 1-4 (groups of 3), width #minus0.5 / 0 / +0.5 MeV");
        f->GetXaxis()->SetTitleSize(0.040); f->GetXaxis()->SetLabelSize(0);
        f->Draw();
        TGraphErrors* g[3]; TGraphErrors* gx = new TGraphErrors();
        const int mk[3] = {23, 20, 22};
        for (int w = 0; w < 3; ++w) { g[w] = new TGraphErrors(); g[w]->SetMarkerStyle(mk[w]); g[w]->SetMarkerColor(kBlue + 1); g[w]->SetLineColor(kBlue + 1); g[w]->SetMarkerSize(1.4); }
        gx->SetMarkerStyle(5); gx->SetMarkerColor(kGray + 1); gx->SetLineColor(kGray + 1); gx->SetMarkerSize(1.4);
        for (int i = 0; i < (int)sd.fits.size(); ++i) {
            const scan::PeakFit& q = sd.fits[i];
            if (!q.ok) continue;
            TGraphErrors* t = q.prob > scan::kVarPmin ? g[i % 3] : gx;
            const int k = t->GetN();
            t->SetPoint(k, i + 0.5, q.S); t->SetPointError(k, 0, q.eS);
        }
        gx->Draw("p same"); for (int w = 0; w < 3; ++w) g[w]->Draw("p same");
        const Res& R = all[1];
        TLine* a = new TLine(0, R.pm.S, 48, R.pm.S); a->SetLineWidth(2); a->Draw();
        TLine* b1 = new TLine(0, R.pm.S + R.pm.eS, 48, R.pm.S + R.pm.eS); b1->SetLineStyle(2); b1->Draw();
        TLine* b2 = new TLine(0, R.pm.S - R.pm.eS, 48, R.pm.S - R.pm.eS); b2->SetLineStyle(2); b2->Draw();
        TLatex t; t.SetTextFont(42); t.SetTextSize(0.035);
        for (int ir = 0; ir < scan::kNVarRange; ++ir) {
            VLine(12 * ir, 2200, 6400, kGray + 1, 3);
            t.DrawLatex(12 * ir + 0.5, 2330, Form("%.3f#minus%.3f", scan::kVarLo[ir], scan::kVarHi[ir]));
            for (int o = 0; o < 4; ++o) t.DrawLatex(12 * ir + 3 * o + 0.6, 2530, Form("pol%d", o + 1));
        }
        Head("oa_lab < 50#circ, isBest #geq 0: peak content S for the 48 model variations of P_main", 0.14, 0.94, 0.038);
        TLegend* lg = new TLegend(0.16, 0.72, 0.62, 0.90); lg->SetTextSize(0.034);
        lg->AddEntry(g[1], Form("width %.1f MeV; triangles: #minus0.5 and +0.5 MeV", 1e3 * scan::kPeakSigma), "pe");
        lg->AddEntry(gx, "fit rejected (p < 0.05)", "p");
        lg->AddEntry(a, Form("reference fit %.0f #pm %.0f (stat); half-range of kept fits: %.0f", R.pm.S, R.pm.eS, sd.half), "l");
        lg->Draw();
        Save(c, "ref_model_variations", false);
    }

    // ---- D6. edge steps ----
    {
        TCanvas* c = new TCanvas("c_step", "", 2000, 1400);
        c->Divide(2, 2);
        for (int i = 0; i < nSl; ++i) {
            c->cd(i + 1);
            const scan::StepFit& f = slFit[i];
            TH1D* h = slSp[i].Hist("sig", Form("st_sig_%d", i)); Deco(h, kBlue + 1, 20, 0.9);
            h->GetXaxis()->SetRangeUser(f.lo - 0.0149, f.hi + 0.0149);
            double ymax = 0, ymin = 1e30;
            for (int b = h->FindBin(f.lo + 1e-6); b <= h->FindBin(f.hi - 1e-6); ++b) {
                ymax = std::max(ymax, h->GetBinContent(b)); ymin = std::min(ymin, h->GetBinContent(b));
            }
            const double yA = std::max(0., ymin - 0.30 * (ymax - ymin)), yB = ymax + 0.55 * (ymax - ymin);
            h->SetMinimum(yA); h->SetMaximum(yB);
            h->Draw("e");
            TGraph* g1 = new TGraph(); TGraph* g0 = new TGraph();
            for (size_t q = 0; q < f.tx.size(); ++q) { g1->SetPoint(q, f.tx[q], f.fy[q]); g0->SetPoint(q, f.tx[q], f.fy0[q]); }
            g1->SetLineWidth(2); g0->SetLineWidth(2); g0->SetLineStyle(2);
            g0->Draw("l same"); g1->Draw("l same");
            VLine(f.edge, yA, yA + 0.45 * (yB - yA), kGreen + 2, 7);
            Head(Form("N #minus CB, oa_lab < 50#circ, %.2f < M(#pi^{+}#pi^{-}) < %.2f GeV/c^{2}%s", slLo[i], slHi[i],
                      i == 3 ? " (control)" : ""), 0.14, 0.94, 0.038);
            TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.042);
            t.DrawLatex(0.17, 0.84, Form("step at the %s edge %.4f:", i == 3 ? "hand-set" : "predicted", f.edge));
            t.DrawLatex(0.17, 0.77, Form("h = %.0f #pm %.0f rows / 2.5 MeV,  #Delta#chi^{2} = %.1f", f.h, f.eh, f.dchi2));
            t.DrawLatex(0.17, 0.70, Form("#chi^{2}/ndf = %.1f/%d", f.chi2, f.ndf));
            t.SetTextSize(0.036);
            t.DrawLatex(0.50, 0.19, "full line: fit; dashed: the fit without its step term");
            t.SetTextColor(kGreen + 2);
            t.DrawLatex(0.50, 0.24, i == 3 ? "green: hand-set edge (no edge expected)" : "green: predicted edge (slice centre)");
        }
        Save(c, "ref_edge_steps", false);
    }

    // ---- D7. K0s ----
    {
        TCanvas* c = new TCanvas("c_k0", "", 1200, 900);
        gPad->SetLeftMargin(0.17);
        Deco(hK, kBlack, 20, 0.9);
        hK->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-}) [GeV/c^{2}]"); hK->GetYaxis()->SetTitle("rows / 1.5 MeV/c^{2}");
        hK->GetYaxis()->SetTitleOffset(1.5);
        hK->SetMinimum(0); hK->SetMaximum(1.45 * hK->GetMaximum());
        hK->Draw("e");
        fK->SetLineColor(kBlack); fK->SetLineWidth(2); fK->Draw("same");
        VLine(scan::kMK0, 0, 0.72 * hK->GetMaximum(), kGray + 2, 2);
        Head("all combinations, corrected momenta, #pi^{+}#pi^{-} crossing point > 20 mm downstream of the vertex, doca < 10 mm", 0.17, 0.94, 0.027);
        TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.034);
        t.DrawLatex(0.20, 0.86, Form("K^{0}_{S} #rightarrow #pi^{+}#pi^{-}: position %.4f #pm %.4f GeV/c^{2} (PDG 0.4976, dashed)", kPos, kPosE));
        t.DrawLatex(0.20, 0.80, Form("width %.1f #pm %.1f MeV, mass scale of M(#pi#pi): %.4f #pm %.4f", 1e3 * rK->Parameter(2), 1e3 * rK->ParError(2), kPos / scan::kMK0, kPosE / scan::kMK0));
        t.DrawLatex(0.20, 0.74, Form("for comparison, scale of M(#pi#piee) from the #eta peak: %.4f", scan::kScale));
        Save(c, "ref_k0s_scale", false);
    }

    // ---- D8. like-sign null test ----
    {
        const Res& R = all[1];
        TCanvas* c = new TCanvas("c_null", "", 2000, 850);
        c->Divide(2, 1);
        for (int k = 1; k <= 2; ++k) {
            c->cd(k); gPad->SetLeftMargin(0.17);
            TH1D* h = R.sp.Hist(k == 1 ? "pp" : "mm", Form("null_%d", k)); Deco(h, kRed + 1, 20, 0.9);
            h->GetYaxis()->SetTitleOffset(1.5);
            h->GetXaxis()->SetTitle(k == 1 ? "M(#pi^{+}#pi^{-}e^{+}e^{+}) [GeV/c^{2}]" : "M(#pi^{+}#pi^{-}e^{-}e^{-}) [GeV/c^{2}]");
            h->GetXaxis()->SetRangeUser(0.4951, 0.6349);
            const double a = h->GetBinContent(h->FindBin(0.4975)), b = h->GetBinContent(h->FindBin(0.63));
            h->SetMinimum(0.6 * std::min(a, b)); h->SetMaximum(1.45 * std::max(a, b));
            h->Draw("e");
            const scan::PeakFit& f = k == 1 ? R.lsP : R.lsM;
            Curve(f, true, kBlack, 1)->Draw("l same"); Curve(f, false, kBlack, 2)->Draw("l same");
            VLine(scan::kPeakLo, h->GetMinimum(), 0.8 * h->GetMaximum(), kGray + 2, 2);
            VLine(scan::kPeakHi, h->GetMinimum(), 0.8 * h->GetMaximum(), kGray + 2, 2);
            Head(Form("null test: like-sign pairs %s, oa_lab < 50#circ, isBest #geq 0", k == 1 ? "e^{+}e^{+}" : "e^{#minus}e^{#minus}"), 0.17, 0.94, 0.040);
            TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.042);
            t.DrawLatex(0.20, 0.84, "same peak shape and range, Poisson fit");
            t.DrawLatex(0.20, 0.77, Form("fitted 'peak' content: %.0f #pm %.0f,  #chi^{2}/ndf = %.1f/%d", f.S, f.eS, f.chi2, f.ndf));
        }
        Save(c, "ref_null_likesign", false);
    }

    P("\ntotal run time %.1f s\n", swAll.RealTime());
    if (gOut) fclose(gOut);
    gOut = nullptr;
    printf("written: %s\n", txt.c_str());
}
