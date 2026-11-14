// ============================================================================
// sim_summary.C  --  key figures of the simulation study, one consistent style
// ============================================================================
// HADES p+p, T = 4.53 GeV (feb22), pi+ pi- e+ e-, SMASH std + lep sets.
// Everything is recomputed here from sim_compact.root (SimEngine.h) and the
// data compact file (../ScanEngine.h); the template fit reuses the machinery
// of sim_templatefit.C (namespace stf, included below, not modified). The
// data-derived cut efficiencies of the data study are read from
// tables/cuts_compare_data.txt (written by sim_cuts.C).
//
// Rules: simulated spectrum = signal (T1) + physical background (T2), rows
// owned in partition mode A, std + lep, weight w (PDG BR); fake lepton (T3)
// and different sources (T4) only as dashed information curves;
// sigma = sum w / L_sim, L_sim = 221.05 nb^-1; data (N - CB) / L_data,
// L_data = 6080.5 nb^-1 (UNCONFIRMED: every figure says so); CB = 2 sqrt(N++ N--)
// in 10 MeV M4 bins without a k factor; all weighted sums in double; errors
// of the simulation from event groups. eta -> pi+pi-e+e- (direct decay) is
// NOT in the simulation (eta set missing): every figure says so.
//
// Figures (plots/eta_sim/figs/, png + pdf):
//   summary_a_sim_m4_pre       simulated dsigma/dM4 per component, "pre": log, stacked linear zoom, signals alone
//   summary_b_datasim_abs      data N, CB, N - CB against the simulation, absolute scale, ratio pads (pre, pre+OA>2)
//   summary_c_datasim_R        the same with the simulation x R (R = data/sim in 0.60-0.80 of each selection)
//   summary_d_oa_classes       data/sim per OA(e+e-) class in three M4 regions, and against OA in 1 deg bins
//   summary_e_cut_eff          efficiency of the main cuts: simulated signal and background next to the data-derived values
//   summary_f_templatefit_pre  PRELIMINARY template fit of data N - CB ("pre"): models M1 and M6, pulls
// Table: plots/eta_sim/tables/summary_recomputed.txt (every number this macro computes)
//
// parts (bit mask): 1 (a), 2 (b)+(c), 4 (d), 8 (e), 16 (f); default 31 = all (about 2-3 min, 9 GB)
// Run from /home/przygoda/HADES/FAT/FAT:
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/aclic_summary", kTRUE)' 'plots/eta_sim/sim_summary.C+'
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TStopwatch.h"
#include "TStyle.h"

#include "SimEngine.h"
#include "sim_ghostvar.h"     // review finding 1: ghost-flagged true pairs, variant (a)
#include "sim_templatefit.C"   // stf::FillTpl, FillData, MakeIn, DoFit, Implied, ModelFine, Models (reused, not modified)

namespace ssum {
using namespace simeng;

const char* const kXt = "M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}]";
const char* const kYt = "d#sigma/dM [nb/(GeV/c^{2})]";
const int kColN = 1, kColCB = 633 /* kRed+1 */, kColSig = 601 /* kBlue+1 */;
const double kRLo = 0.60, kRHi = 0.80;               // window of the renormalisation factor R
const double kSimPkLo = 0.5332, kSimPkHi = 0.5588;   // simulated peak window (sim_cuts.C: 0.5460 +- 1.9 sigma)
const std::string kNoteL = "L_{data} = 6080.5 nb^{-1}, UNCONFIRMED: every data value scales with 1/L_{data}";
const std::string kNoteD = "#eta #rightarrow #pi^{+}#pi^{-}e^{+}e^{-} (direct decay) is NOT in the simulation (eta set not available)";
const std::string kNoteS = "simulation: SMASH std + lep, owned rows (mode A), w (PDG BR); #sigma = #Sigma w / L_{sim}, L_{sim} = 221.05 nb^{-1}; no efficiency correction";
const std::string kNoteC = "data: CB = 2#sqrt{N_{++}N_{--}} in 10 MeV bins (no k factor); fake-lepton and different-source rows of the simulation are NOT in the simulated spectrum";

FILE* gOut = nullptr;
void P(const char* fmt, ...) {
    char buf[8192];
    va_list ap; va_start(ap, fmt); vsnprintf(buf, sizeof(buf), fmt, ap); va_end(ap);
    printf("%s", buf); fflush(stdout);
    if (gOut) { fputs(buf, gOut); fflush(gOut); }
}

// ---------------------------------------------------------------- plotting helpers
void Txt(double x, double y, const char* s, double size, int col = 1, int align = 11, int font = 42) {
    TLatex t; t.SetNDC(); t.SetTextFont(font); t.SetTextSize(size); t.SetTextColor(col); t.SetTextAlign(align); t.DrawLatex(x, y, s);
}
TPad* MkPad(const char* name, double x1, double y1, double x2, double y2, bool logy, double lm, double rm, double tm, double bm) {
    TPad* p = new TPad(name, "", x1, y1, x2, y2);
    p->SetLeftMargin(lm); p->SetRightMargin(rm); p->SetTopMargin(tm); p->SetBottomMargin(bm);
    p->SetLogy(logy); p->SetTicks(1, 1); p->SetFillStyle(0);
    p->Draw(); p->cd();
    return p;
}
// label size ls and title size ts in NDC of a pad of height h (canvas fraction): absolute sizes a, t (canvas fraction)
void Ax(TH1* h, const char* xt, const char* yt, double hpad, double a = 0.0150, double t = 0.0170, double yoff = 1.0, double xoff = 1.0) {
    h->GetXaxis()->SetTitle(xt); h->GetYaxis()->SetTitle(yt);
    h->GetXaxis()->SetLabelSize(a / hpad); h->GetYaxis()->SetLabelSize(a / hpad);
    h->GetXaxis()->SetTitleSize(t / hpad); h->GetYaxis()->SetTitleSize(t / hpad);
    h->GetYaxis()->SetTitleOffset(yoff); h->GetXaxis()->SetTitleOffset(xoff);
    h->GetXaxis()->SetTickLength(std::min(0.06, 0.012 / hpad)); h->GetYaxis()->SetNdivisions(505);
    h->GetYaxis()->SetMaxDigits(4);
    h->SetStats(0);
}
// y-title offset that puts the title at 78 % of the left margin lm (ROOT scales the y offset with the pad width;
// t = absolute title size as a fraction of the canvas height, hpad = pad height as a fraction of the canvas)
double YO(double lm, double t, double hpad) { return 0.78 * lm / (1.6 * t / hpad); }
void Ln(TH1D* h, int col, int lw, int ls, int fill = 0) {
    h->SetLineColor(col); h->SetLineWidth(lw); h->SetLineStyle(ls); h->SetMarkerSize(0); h->SetMarkerColor(col);
    if (fill) { h->SetFillColor(fill); h->SetFillStyle(1001); } else h->SetFillStyle(0);
}
void Mk(TH1* h, int col, int st, double sz = 1.2) {
    h->SetMarkerStyle(st); h->SetMarkerSize(sz); h->SetMarkerColor(col); h->SetLineColor(col); h->SetLineWidth(2);
}
void MkG(TGraph* g, int col, int st, double sz = 1.6) {
    g->SetMarkerStyle(st); g->SetMarkerSize(sz); g->SetMarkerColor(col); g->SetLineColor(col); g->SetLineWidth(2);
}
TH1D* Frame(const char* name, double lo, double hi, double ymin, double ymax) {
    TH1D* f = new TH1D(name, "", 1, lo, hi);
    f->SetDirectory(nullptr); f->SetMinimum(ymin); f->SetMaximum(ymax); f->SetStats(0);
    return f;
}
struct HL { std::string s; int col; bool bold; };
// header over the full width of a canvas (y1..1): one text line per entry
void Header(TCanvas* c, const std::vector<HL>& L, double y1) {
    c->cd();
    TPad* p = new TPad(Form("hdr_%s", c->GetName()), "", 0, y1, 1, 1);
    p->SetMargin(0, 0, 0, 0); p->SetFillStyle(0); p->Draw(); p->cd();
    const int n = (int)L.size();
    const double h = 1.0 / (n + 0.3);
    for (int i = 0; i < n; ++i) Txt(0.010, 1 - (i + 0.88) * h, L[i].s.c_str(), 0.70 * h, L[i].col, 11, L[i].bold ? 62 : 42);
    c->cd();
}
void Save(TCanvas* c, const std::string& base) {
    c->SaveAs((base + ".png").c_str());
    c->SaveAs((base + ".pdf").c_str());
}
TLine* HLine(double x1, double x2, double y, int col, int ls, int lw = 2) {
    TLine* l = new TLine(x1, y, x2, y); l->SetLineColor(col); l->SetLineStyle(ls); l->SetLineWidth(lw); l->Draw(); return l;
}
TLine* VLine(double x, double y1, double y2, int col, int ls, int lw = 1) {
    TLine* l = new TLine(x, y1, x, y2); l->SetLineColor(col); l->SetLineStyle(ls); l->SetLineWidth(lw); l->Draw(); return l;
}
TLegend* Leg(double x1, double y1, double x2, double y2, double ts, int ncol = 1) {
    TLegend* l = new TLegend(x1, y1, x2, y2);
    l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextFont(42); l->SetTextSize(ts); l->SetNColumns(ncol);
    return l;
}

// ---------------------------------------------------------------- data and simulation helpers
// data histogram in nb/(GeV/c2): what = us, pp, mm, cb, sig (errors: Poisson for N, ScanEngine CB propagation)
TH1D* DH(const scan::Spectra& sp, const char* what, const char* name, double Ld = kLdata) {
    TH1D* h = sp.Hist(what, name);
    const double f = 1.0 / (Ld * sp.bw);
    for (int i = 1; i <= h->GetNbinsX(); ++i) { h->SetBinContent(i, h->GetBinContent(i) * f); h->SetBinError(i, h->GetBinError(i) * f); }
    return h;
}
std::vector<scan::Spectra> DFillMany(const scan::Data& dd, const std::vector<scan::Sel>& sels, int nbin, double lo, double hi, int reb,
                                     const std::vector<std::pair<double, double>>& win) {
    scan::FillOpt fo; fo.nbin = nbin; fo.lo = lo; fo.hi = hi; fo.reb = reb; fo.windows = win;
    return scan::FillMany(dd, sels, fo);
}
// ratio num / (f * den) per bin; bins with den <= 0, num <= 0 or a relative error of den above maxRel are left empty
TH1D* Ratio(const TH1D* num, const TH1D* den, const char* name, double f = 1.0, double maxRel = 0.35) {
    TH1D* r = (TH1D*)num->Clone(name);
    r->Reset(); r->SetDirectory(nullptr);
    for (int i = 1; i <= r->GetNbinsX(); ++i) {
        const double n = num->GetBinContent(i), en = num->GetBinError(i);
        const double dv = den->GetBinContent(i) * f, ed = den->GetBinError(i) * f;
        if (dv <= 0 || n <= 0 || ed / dv > maxRel) continue;
        const double v = n / dv;
        r->SetBinContent(i, v); r->SetBinError(i, v * std::sqrt(en * en / (n * n) + ed * ed / (dv * dv)));
    }
    return r;
}
// simulated sum in [a, b) of one component with the exact event-group variance (one-bin fill)
struct XS { double s = 0, e = 0, neff = 0, rows = 0; };
XS SimWin(const SimData& d, const Sel& sel, const Comp& c, double a, double b, const SimOpt& o = SimOpt()) {
    Binning bb; bb.nbin = 1; bb.lo = a; bb.hi = b;
    const WSpec w = Fill(d, sel, c, o, bb);
    XS x; x.s = w.S / kLsim; x.e = std::sqrt(w.V) / kLsim; x.neff = w.NeffIn(); x.rows = w.N;
    return x;
}
struct RatioV { double r = 0, e = 0; };
RatioV Div(double a, double ea, double b, double eb) {
    RatioV q;
    if (a > 0 && b > 0) { q.r = a / b; q.e = q.r * std::sqrt(ea * ea / (a * a) + eb * eb / (b * b)); }
    return q;
}
struct Region { std::string name; double a, b; };
std::vector<Region> Regions() {
    return {{"below the bump", 0.30, 0.45},           {"bump window", scan::kBumpLo, scan::kBumpHi},
            {"gap bump-peak", 0.525, 0.530},           {"peak window", scan::kPeakLo, scan::kPeakHi},
            {"above the peak 0.56-0.70", 0.56, 0.70},  {"0.70-0.90", 0.70, 0.90},
            {"high mass 0.90-1.10", 0.90, 1.10},       {"R window 0.60-0.80", kRLo, kRHi},
            {"all 0.30-1.10", 0.30, 1.10}};
}
const int kIdxR = 7;    // index of the R window in Regions()

// ============================================================================
// (a) simulated dsigma/dM4 per component, "pre"
// ============================================================================
void PartA(const SimData& d, const std::string& fig) {
    P("\n## (a) simulated cross section per component, pre (isBest >= 0, PT3, vz > -160, cor_ok), mode A, w (PDG), std + lep\n");
    const std::vector<Comp> comps = {CompBumpDalitz(), CompBumpConv(), CompPeakConv(), CompPhysBkg(), CompSim(), CompFake(), CompDiff(), CompSignal()};
    const std::vector<WSum> t = SumMany(d, SelPreset(kPre), comps, SimOpt());
    P("%-34s %10s %12s %20s %10s %14s\n", "component (all M4)", "rows", "sum w", "sigma [nb]", "N_eff", "data rows*");
    for (size_t i = 0; i < comps.size(); ++i)
        P("%-34s %10.0f %12.2f %11.3f +- %6.3f %10.0f %14.0f\n", comps[i].name.c_str(), t[i].N, t[i].S, t[i].Xsec(), t[i].XsecErr(), t[i].Neff(), t[i].DataRows());
    P("* expected data rows = sum w * L_data / L_sim with L_data = 6080.5 nb^-1 (UNCONFIRMED)\n");
    // windows
    P("%-34s %22s %22s\n", "component in a window", "bump 0.4550-0.5225 [nb]", "peak 0.5336-0.5542 [nb]");
    for (int i : {0, 1, 2, 7, 3, 4, 5, 6}) {
        const XS b = SimWin(d, SelPreset(kPre), comps[i], scan::kBumpLo, scan::kBumpHi);
        const XS p = SimWin(d, SelPreset(kPre), comps[i], scan::kPeakLo, scan::kPeakHi);
        P("%-34s %10.3f +- %6.3f (N_eff %6.0f) %8.3f +- %6.3f (N_eff %6.0f)\n", comps[i].name.c_str(), b.s, b.e, b.neff, p.s, p.e, p.neff);
    }

    Binning b5; b5.nbin = 160; b5.lo = 0.30; b5.hi = 1.10;
    const std::vector<WSpec> w = FillMany(d, SelPreset(kPre), comps, SimOpt(), b5);
    Binning b25; b25.nbin = 88; b25.lo = 0.40; b25.hi = 0.62;
    const std::vector<WSpec> wz = FillMany(d, SelPreset(kPre), {CompBumpDalitz(), CompBumpConv(), CompPeakConv(), CompSignal()}, SimOpt(), b25);
    std::vector<TH1D*> h(comps.size());
    for (size_t i = 0; i < comps.size(); ++i) h[i] = w[i].Hist(Form("hA_%zu", i), 2);

    TCanvas* c = new TCanvas("cSumA", "", 2000, 1700);
    Header(c, {{"(a) Simulated cross section within the HADES acceptance after the selection: d#sigma/dM_{#pi^{+}#pi^{-}e^{+}e^{-}} per component", 1, true},
               {"selection \"pre\": isBest #geq 0, PT3, v_{z} > -160 mm, corrected momenta (e^{+}e^{-}); 5 MeV bins (panel c: 2.5 MeV)", 1, false},
               {kNoteS, 1, false},
               {kNoteD, kRed + 1, false}},
           0.905);
    // ---- (a) log, all components separately
    c->cd();
    MkPad("pA1", 0.0, 0.475, 1.0, 0.905, true, 0.075, 0.015, 0.03, 0.15);
    const double H1 = 0.43;
    TH1D* f1 = Frame("fA1", 0.30, 1.10, 0.05, 3e6);
    Ax(f1, kXt, kYt, H1, 0.0185, 0.0200, YO(0.075, 0.0200, H1));
    f1->Draw("AXIS");
    Ln(h[3], kColPhysBkg, 1, 1, kColPhysFill); h[3]->Draw("HIST SAME");
    Ln(h[0], kColBumpDalitz, 3, 1); h[0]->Draw("HIST SAME");
    Ln(h[1], kColBumpConv, 3, 1); h[1]->Draw("HIST SAME");
    Ln(h[2], kColPeakConv, 3, 1); h[2]->Draw("HIST SAME");
    Ln(h[5], kColFake, 3, 2); h[5]->Draw("HIST SAME");
    Ln(h[6], kColDiff, 3, 2); h[6]->Draw("HIST SAME");
    Ln(h[4], kColSim, 4, 1); h[4]->Draw("HIST SAME");
    gPad->RedrawAxis();
    Txt(0.085, 0.95, "(a) all components, logarithmic", 0.046, 1, 13, 62);
    TLegend* lg = Leg(0.085, 0.67, 0.985, 0.89, 0.040, 2);
    lg->AddEntry(h[4], "simulated spectrum = signal + physical background", "l");
    lg->AddEntry(h[0], "signal: bump, #eta#rightarrow#pi^{+}#pi^{-}#pi^{0}, #pi^{0} Dalitz", "l");
    lg->AddEntry(h[3], "physical background (one correlated e^{+}e^{-} source)", "f");
    lg->AddEntry(h[1], "signal: bump, #eta#rightarrow#pi^{+}#pi^{-}#pi^{0}, converted #gamma", "l");
    lg->AddEntry(h[5], "fake lepton (information only, not in the spectrum)", "l");
    lg->AddEntry(h[2], "signal: peak, #eta#rightarrow#pi^{+}#pi^{-}#gamma, converted #gamma", "l");
    lg->AddEntry(h[6], "different sources (information only, not in the spectrum)", "l");
    lg->Draw();
    // ---- (b) linear zoom, stacked
    c->cd();
    MkPad("pA2", 0.0, 0.0, 0.5, 0.475, false, 0.16, 0.02, 0.10, 0.15);
    const double H2 = 0.475;
    TH1D* s1 = (TH1D*)h[3]->Clone("hA_s1"); TH1D* s2 = (TH1D*)s1->Clone("hA_s2"); s2->Add(h[0]);
    TH1D* s3 = (TH1D*)s2->Clone("hA_s3"); s3->Add(h[1]); TH1D* s4 = (TH1D*)s3->Clone("hA_s4"); s4->Add(h[2]);
    double ymax = 0;
    for (int i = 1; i <= s4->GetNbinsX(); ++i) if (s4->GetBinCenter(i) > 0.40 && s4->GetBinCenter(i) < 0.70) ymax = std::max(ymax, s4->GetBinContent(i));
    TH1D* f2 = Frame("fA2", 0.40, 0.70, 0, 1.30 * ymax);
    Ax(f2, kXt, kYt, H2, 0.0185, 0.0200, YO(0.16, 0.0200, H2));
    f2->Draw("AXIS");
    Ln(s4, kColPeakConv, 1, 1, kColPeakConv); s4->Draw("HIST SAME");
    Ln(s3, kColBumpConv, 1, 1, kColBumpConv); s3->Draw("HIST SAME");
    Ln(s2, kColBumpDalitz, 1, 1, kColBumpDalitz); s2->Draw("HIST SAME");
    Ln(s1, kColPhysBkg, 1, 1, kColPhysFill); s1->Draw("HIST SAME");
    TH1D* fk = (TH1D*)h[5]->Clone("hA_fk"); fk->Draw("HIST SAME");
    gPad->RedrawAxis();
    Txt(0.17, 0.975, "(b) linear zoom, components stacked", 0.046, 1, 13, 62);
    TLegend* lg2 = Leg(0.18, 0.56, 0.62, 0.88, 0.040);
    lg2->AddEntry(s4, "+ peak, conversion", "f");
    lg2->AddEntry(s3, "+ bump, conversion", "f");
    lg2->AddEntry(s2, "+ bump, #pi^{0} Dalitz", "f");
    lg2->AddEntry(s1, "physical background", "f");
    lg2->AddEntry(fk, "fake lepton (not stacked, info)", "l");
    lg2->Draw();
    // ---- (c) the signal channels alone, 2.5 MeV, with errors
    c->cd();
    MkPad("pA3", 0.5, 0.0, 1.0, 0.475, false, 0.16, 0.02, 0.10, 0.15);
    std::vector<TH1D*> hz(4);
    double ym3 = 0;
    for (int i = 0; i < 4; ++i) {
        hz[i] = wz[i].Hist(Form("hA_z%d", i), 2);
        for (int k = 1; k <= hz[i]->GetNbinsX(); ++k) ym3 = std::max(ym3, hz[i]->GetBinContent(k) + hz[i]->GetBinError(k));
    }
    TH1D* f3 = Frame("fA3", 0.40, 0.62, 0, 1.75 * ym3);
    Ax(f3, kXt, kYt, H2, 0.0185, 0.0200, YO(0.16, 0.0200, H2));
    f3->Draw("AXIS");
    const int col3[4] = {kColBumpDalitz, kColBumpConv, kColPeakConv, 1};
    for (int i = 0; i < 4; ++i) {
        Ln(hz[i], col3[i], i == 3 ? 2 : 3, i == 3 ? 7 : 1); hz[i]->SetMarkerSize(0); hz[i]->Draw("HIST SAME");
        if (i < 3) hz[i]->Draw("E1 SAME");
    }
    Txt(0.17, 0.975, "(c) the three signal channels alone (2.5 MeV bins)", 0.046, 1, 13, 62);
    TLegend* lg3 = Leg(0.18, 0.64, 0.80, 0.88, 0.038);
    lg3->AddEntry(hz[0], Form("bump, #pi^{0} Dalitz: %.2f #pm %.2f nb", t[0].Xsec(), t[0].XsecErr()), "l");
    lg3->AddEntry(hz[1], Form("bump, conversion: %.2f #pm %.2f nb", t[1].Xsec(), t[1].XsecErr()), "l");
    lg3->AddEntry(hz[2], Form("peak, conversion: %.2f #pm %.2f nb", t[2].Xsec(), t[2].XsecErr()), "l");
    lg3->AddEntry(hz[3], Form("sum of the signals: %.2f #pm %.2f nb", t[7].Xsec(), t[7].XsecErr()), "l");
    lg3->Draw();
    Txt(0.18, 0.595, "#sigma: whole M_{4} range; sum of signals without errors", 0.034, kRed + 1);
    Save(c, fig + "summary_a_sim_m4_pre");
}

// ============================================================================
// (b), (c) data N - CB against the simulated sum: absolute and x R
// ============================================================================
struct ColBC {
    int preset;
    std::vector<WSpec> w5;        // 0 sum, 1 phys bkg, 2 bump D, 3 bump C, 4 peak C, 5 fake, 6 diff
    scan::Spectra sp;             // data, 5 MeV, CB 10 MeV, region windows
    std::vector<XS> simReg;       // T1 + T2 per region
    std::vector<XS> simRegA;      // variant (a): T1 + T2 + ghost-flagged true pairs, per region
    double R = 1, eR = 0, Ra = 1, eRa = 0;
};
void DrawBC(std::vector<ColBC>& cols, bool useR, const std::string& fig) {
    const char* tag = useR ? "R" : "abs";
    TCanvas* c = new TCanvas(Form("cSum%s", tag), "", 2000, 2000);
    std::vector<HL> hdr;
    if (!useR) hdr.push_back({"(b) Data (N - CB) / L_{data} against the simulated spectrum (signal + physical background) in absolute units, NOT renormalised", 1, true});
    else hdr.push_back({"(c) The same with the simulation multiplied by R = data (N - CB) / simulation in M_{4} = 0.60 - 0.80 GeV/c^{2} of each selection", 1, true});
    hdr.push_back({"data: black N (all e^{+}e^{-} combinations), red CB, blue N - CB; " + kNoteL, kColSig, false});
    hdr.push_back({kNoteS, 1, false});
    hdr.push_back({kNoteC, 1, false});
    hdr.push_back({kNoteD, kRed + 1, false});
    hdr.push_back({std::string("ratio pads: data (N - CB) / simulation") + (useR ? " #times R" : "") + "; thick black lines: averages over the regions below the bump, bump window, peak window, 0.56-0.70, 0.70-0.90, 0.90-1.10", 1, false});
    hdr.push_back({"dashed teal / open green: variant (a) = T1 + T2 + ghost-flagged true e^{+}e^{-} pairs (91 % of the 'fake lepton' class; treatment to be confirmed by the user)", kGreen + 3, false});
    hdr.push_back({"R depends on the M_{4} window, the v_{z} slice and the ghost treatment (0.93 - 1.45 for v_{z} > -160 and -140 < v_{z} < -100): absolute scale not determined better than 15 - 20 % (tables/review_R_variants.txt)", kRed + 1, false});
    Header(c, hdr, 0.905);
    // common y ranges
    double zmax = 0;
    for (auto& col : cols) {
        const double f = useR ? col.R : 1.0;
        TH1D* sig = DH(col.sp, "sig", "tmpz");
        for (int i = 1; i <= sig->GetNbinsX(); ++i) {
            const double x = sig->GetBinCenter(i);
            if (x > 0.40 && x < 0.70) zmax = std::max({zmax, sig->GetBinContent(i), f * col.w5[0].s[i - 1] / (kLsim * col.w5[0].bw)});
        }
        delete sig;
    }
    const std::vector<Region> reg = Regions();
    for (size_t k = 0; k < cols.size(); ++k) {
        ColBC& col = cols[k];
        const double f = useR ? col.R : 1.0;
        const double x0 = 0.5 * k, x1 = x0 + 0.5;
        const std::string pn = Form("%s_%zu", tag, k);
        const double fA = useR ? col.Ra : 1.0;
        std::vector<TH1D*> hs(8);
        for (int i = 0; i < 8; ++i) { hs[i] = col.w5[i].Hist(Form("hS_%s_%d", pn.c_str(), i), 2); hs[i]->Scale(i == 7 ? fA : f); }
        TH1D* hN = DH(col.sp, "us", Form("hN_%s", pn.c_str()));
        TH1D* hC = DH(col.sp, "cb", Form("hC_%s", pn.c_str()));
        TH1D* hD = DH(col.sp, "sig", Form("hD_%s", pn.c_str()));
        Mk(hN, kColN, 20, 0.9); Mk(hC, kColCB, 21, 0.8); Mk(hD, kColSig, 20, 0.9);
        TH1D* hr = Ratio(hD, hs[0], Form("hR_%s", pn.c_str()));
        Mk(hr, kColSig, 20, 0.9);
        TH1D* hrA = Ratio(hD, hs[7], Form("hRA_%s", pn.c_str()));
        Mk(hrA, kGreen + 2, 24, 0.8);
        const std::string selName = PresetName(col.preset);
        // ---- log pad
        c->cd();
        MkPad(Form("pL_%s", pn.c_str()), x0, 0.60, x1, 0.905, true, 0.16, 0.02, 0.05, 0.0);
        const double HL1 = 0.305;
        TH1D* fl = Frame(Form("fL_%s", pn.c_str()), 0.30, 1.10, 0.3, 3e7);
        Ax(fl, kXt, kYt, HL1, 0.0175, 0.0190, YO(0.16, 0.0190, HL1));
        fl->Draw("AXIS");
        Ln(hs[1], kColPhysBkg, 1, 1, kColPhysFill); hs[1]->Draw("HIST SAME");
        Ln(hs[2], kColBumpDalitz, 2, 1); hs[2]->Draw("HIST SAME");
        Ln(hs[3], kColBumpConv, 2, 1); hs[3]->Draw("HIST SAME");
        Ln(hs[4], kColPeakConv, 2, 1); hs[4]->Draw("HIST SAME");
        Ln(hs[5], kColFake, 2, 2); hs[5]->Draw("HIST SAME");
        Ln(hs[6], kColDiff, 2, 2); hs[6]->Draw("HIST SAME");
        Ln(hs[0], kColSim, 3, 1); hs[0]->Draw("HIST SAME");
        Ln(hs[7], kColSim, 3, 2); hs[7]->Draw("HIST SAME");
        hN->Draw("P SAME"); hC->Draw("P SAME"); hD->Draw("P SAME");
        gPad->RedrawAxis();
        const std::string scl = useR ? std::string(Form("sim #times R = %.3f #pm %.3f; (a) #times %.3f", col.R, col.eR, col.Ra)) : std::string("absolute scale");
        Txt(0.17, 0.90, (selName + " (" + scl + ")").c_str(), 0.050, 1, 13, 62);
        TLegend* lg = Leg(0.18, 0.66, 0.42, 0.84, 0.042);
        lg->AddEntry(hN, "data: N", "p");
        lg->AddEntry(hC, "data: CB", "p");
        lg->AddEntry(hD, "data: N - CB", "p");
        lg->Draw();
        TLegend* lgs = Leg(0.43, 0.48, 0.98, 0.84, 0.036);
        lgs->AddEntry(hs[0], useR ? "sim: signal + phys. bkg (#times R)" : "sim: signal + phys. bkg", "l");
        lgs->AddEntry(hs[7], useR ? "sim (a): + ghost-flagged true pairs (#times R_{a})" : "sim (a): + ghost-flagged true pairs", "l");
        lgs->AddEntry(hs[1], "sim: physical background", "f");
        lgs->AddEntry(hs[2], "sim: bump, #pi^{0} Dalitz", "l");
        lgs->AddEntry(hs[3], "sim: bump, conversion", "l");
        lgs->AddEntry(hs[4], "sim: peak, conversion", "l");
        lgs->AddEntry(hs[5], "sim: fake lepton (info only)", "l");
        lgs->AddEntry(hs[6], "sim: different sources (info only)", "l");
        lgs->Draw();
        // ---- ratio pad
        c->cd();
        MkPad(Form("pLr_%s", pn.c_str()), x0, 0.47, x1, 0.60, false, 0.16, 0.02, 0.0, 0.30);
        const double HL2 = 0.13;
        TH1D* fr = Frame(Form("fLr_%s", pn.c_str()), 0.30, 1.10, useR ? 0.3 : 0.3, useR ? 1.6 : 1.9);
        Ax(fr, kXt, useR ? "data/(R sim)" : "data/sim", HL2, 0.0175, 0.0160, YO(0.16, 0.0160, HL2), 0.95);
        fr->GetXaxis()->SetTitleSize(0.0190 / HL2);
        fr->Draw("AXIS");
        HLine(0.30, 1.10, 1.0, kGray + 2, 2, 1);
        for (double xb : {scan::kBumpLo, scan::kBumpHi, scan::kPeakLo, scan::kPeakHi}) VLine(xb, useR ? 0.3 : 0.3, useR ? 1.6 : 1.9, kGray + 1, 3);
        hrA->Draw("P SAME");
        hr->Draw("P SAME");
        for (int ir : {0, 1, 3, 4, 5, 6}) {
            const scan::WinSum& ws = col.sp.win[2 + ir];
            const RatioV q = Div(ws.sig() / kLdata, ws.sigErr() / kLdata, f * col.simReg[ir].s, f * col.simReg[ir].e);
            HLine(reg[ir].a, reg[ir].b, q.r, 1, 1, 4);
        }
        // ---- zoom pad
        c->cd();
        MkPad(Form("pZ_%s", pn.c_str()), x0, 0.165, x1, 0.455, false, 0.16, 0.02, 0.05, 0.0);
        const double HZ = 0.29;
        TH1D* s1 = (TH1D*)hs[1]->Clone(Form("hZ1_%s", pn.c_str()));
        TH1D* s2 = (TH1D*)s1->Clone(Form("hZ2_%s", pn.c_str())); s2->Add(hs[2]);
        TH1D* s3 = (TH1D*)s2->Clone(Form("hZ3_%s", pn.c_str())); s3->Add(hs[3]);
        TH1D* s4 = (TH1D*)s3->Clone(Form("hZ4_%s", pn.c_str())); s4->Add(hs[4]);
        TH1D* fz = Frame(Form("fZ_%s", pn.c_str()), 0.40, 0.70, 0, 1.12 * zmax);
        Ax(fz, kXt, kYt, HZ, 0.0175, 0.0190, YO(0.16, 0.0190, HZ));
        fz->Draw("AXIS");
        Ln(s4, kColPeakConv, 1, 1, kColPeakConv); s4->Draw("HIST SAME");
        Ln(s3, kColBumpConv, 1, 1, kColBumpConv); s3->Draw("HIST SAME");
        Ln(s2, kColBumpDalitz, 1, 1, kColBumpDalitz); s2->Draw("HIST SAME");
        Ln(s1, kColPhysBkg, 1, 1, kColPhysFill); s1->Draw("HIST SAME");
        TH1D* hD2 = (TH1D*)hD->Clone(Form("hDz_%s", pn.c_str())); hD2->Draw("P SAME");
        TH1D* hC2 = (TH1D*)hC->Clone(Form("hCz_%s", pn.c_str())); hC2->Draw("P SAME");
        gPad->RedrawAxis();
        Txt(0.17, 0.92, "linear zoom, simulation stacked", 0.050, 1, 13, 62);
        TLegend* lz = Leg(0.17, 0.50, 0.62, 0.85, 0.045);
        lz->AddEntry(hD2, "data: N - CB", "p");
        lz->AddEntry(hC2, "data: CB", "p");
        lz->AddEntry(s4, "sim: + peak, conversion", "f");
        lz->AddEntry(s3, "sim: + bump, conversion", "f");
        lz->AddEntry(s2, "sim: + bump, #pi^{0} Dalitz", "f");
        lz->AddEntry(s1, "sim: physical background", "f");
        lz->Draw();
        // ---- zoom ratio
        c->cd();
        MkPad(Form("pZr_%s", pn.c_str()), x0, 0.0, x1, 0.165, false, 0.16, 0.02, 0.0, 0.33);
        const double HZ2 = 0.165;
        TH1D* fzr = Frame(Form("fZr_%s", pn.c_str()), 0.40, 0.70, useR ? 0.5 : 0.6, useR ? 1.5 : 1.8);
        Ax(fzr, kXt, useR ? "data/(R sim)" : "data/sim", HZ2, 0.0175, 0.0160, YO(0.16, 0.0160, HZ2), 0.95);
        fzr->GetXaxis()->SetTitleSize(0.0190 / HZ2);
        fzr->Draw("AXIS");
        HLine(0.40, 0.70, 1.0, kGray + 2, 2, 1);
        for (double xb : {scan::kBumpLo, scan::kBumpHi, scan::kPeakLo, scan::kPeakHi}) VLine(xb, useR ? 0.5 : 0.6, useR ? 1.5 : 1.8, kGray + 1, 3);
        TH1D* hrA2 = (TH1D*)hrA->Clone(Form("hRAz_%s", pn.c_str())); hrA2->Draw("P SAME");
        TH1D* hr2 = (TH1D*)hr->Clone(Form("hRz_%s", pn.c_str())); hr2->Draw("P SAME");
        for (int ir : {0, 1, 3, 4}) {
            const scan::WinSum& ws = col.sp.win[2 + ir];
            const RatioV q = Div(ws.sig() / kLdata, ws.sigErr() / kLdata, f * col.simReg[ir].s, f * col.simReg[ir].e);
            HLine(std::max(0.40, reg[ir].a), std::min(0.70, reg[ir].b), q.r, 1, 1, 4);
        }
    }
    Save(c, fig + (useR ? "summary_c_datasim_R" : "summary_b_datasim_abs"));
}
void PartBC(const SimData& d, const scan::Data& dd, const std::string& fig) {
    const std::vector<Region> reg = Regions();
    std::vector<std::pair<double, double>> win;
    for (const auto& r : reg) win.push_back(std::make_pair(r.a, r.b));
    std::vector<ColBC> cols;
    Binning b5; b5.nbin = 160; b5.lo = 0.30; b5.hi = 1.10;
    const std::vector<Comp> comps = {CompSim(), CompPhysBkg(), CompBumpDalitz(), CompBumpConv(), CompPeakConv(), CompFake(), CompDiff(), CompSimVar(1)};
    const std::vector<scan::Spectra> sps = DFillMany(dd, {DataSelPreset(kPre), DataSelPreset(kPreOA2)}, 160, 0.30, 1.10, 2, win);
    int k = 0;
    for (int p : {kPre, kPreOA2}) {
        ColBC col; col.preset = p;
        col.w5 = FillMany(d, SelPreset(p), comps, SimOpt(), b5);
        col.sp = sps[k++];
        for (const auto& r : reg) col.simReg.push_back(SimWin(d, SelPreset(p), CompSim(), r.a, r.b));
        for (const auto& r : reg) col.simRegA.push_back(SimWin(d, SelPreset(p), CompSimVar(1), r.a, r.b));
        const scan::WinSum& wr = col.sp.win[2 + kIdxR];
        const RatioV q = Div(wr.sig() / kLdata, wr.sigErr() / kLdata, col.simReg[kIdxR].s, col.simReg[kIdxR].e);
        col.R = q.r; col.eR = q.e;
        const RatioV qa = Div(wr.sig() / kLdata, wr.sigErr() / kLdata, col.simRegA[kIdxR].s, col.simRegA[kIdxR].e);
        col.Ra = qa.r; col.eRa = qa.e;
        P("\n## (b,c) %s: data (N - CB) / L_data against the simulation T1 + T2 (mode A, w); R = %.4f +- %.4f (0.60-0.80); variant (a) R = %.4f +- %.4f\n", PresetName(p), col.R, col.eR, col.Ra, col.eRa);
        P("%-34s %9s %8s %18s %18s %8s %16s %16s %16s\n", "region", "N", "CB", "data [nb]", "sim [nb]", "N_eff", "data/sim", "data/(R sim)", "data/sim (a)");
        for (size_t i = 0; i < reg.size(); ++i) {
            const scan::WinSum& ws = col.sp.win[2 + i];
            const XS& s = col.simReg[i];
            const RatioV q1 = Div(ws.sig() / kLdata, ws.sigErr() / kLdata, s.s, s.e);
            const RatioV qa1 = Div(ws.sig() / kLdata, ws.sigErr() / kLdata, col.simRegA[i].s, col.simRegA[i].e);
            P("%-34s %9.0f %8.0f %9.4f +- %6.4f %9.3f +- %6.3f %8.0f %7.3f +- %6.3f %7.3f +- %6.3f %7.3f +- %6.3f\n", reg[i].name.c_str(), ws.n[0], ws.cb, ws.sig() / kLdata,
              ws.sigErr() / kLdata, s.s, s.e, s.neff, q1.r, q1.e, q1.r / col.R, q1.e / col.R, qa1.r, qa1.e);
        }
        // R variants
        P("R variants (0.60-0.80):");
        const scan::WinSum& wr2 = col.sp.win[2 + kIdxR];
        for (int v = 0; v < 3; ++v) {
            SimOpt o; if (v == 0) o.mode = 1; if (v == 1) o.db = true; if (v == 2) { o.mode = 1; o.db = true; }
            const XS s = SimWin(d, SelPreset(p), CompSim(), kRLo, kRHi, o);
            const RatioV q = Div(wr2.sig() / kLdata, wr2.sigErr() / kLdata, s.s, s.e);
            P("  %s %.3f +- %.3f;", v == 0 ? "mode B, w" : v == 1 ? "mode A, w_db" : "mode B, w_db (not physical)", q.r, q.e);
        }
        P("\n");
        cols.push_back(col);
    }
    DrawBC(cols, false, fig);
    DrawBC(cols, true, fig);
}

// ============================================================================
// (d) data / simulation per OA(e+e-) class
// ============================================================================
void PartD(const SimData& d, const scan::Data& dd, const std::string& fig, double Rpre, double eRpre) {
    const int nC = 5;
    const double cLo[nC] = {0, 2, 4, 6, 10}, cHi[nC] = {2, 4, 6, 10, 1e9};
    const char* cLab[nC] = {"OA < 2#circ", "2 - 4#circ", "4 - 6#circ", "6 - 10#circ", "OA > 10#circ"};
    const char* cTxt[nC] = {"OA < 2", "2-4", "4-6", "6-10", "> 10"};
    const int nR = 3;
    const double rLo[nR] = {0.30, scan::kBumpLo, 0.56}, rHi[nR] = {0.45, scan::kBumpHi, 0.80};
    const char* rLab[nR] = {"below the bump, 0.30 - 0.45", "bump window, 0.455 - 0.5225 (#eta signals)", "above the peak, 0.56 - 0.80"};
    const int rCol[nR] = {kOrange + 7, kMagenta + 1, kColSig};
    const int rSt[nR] = {21, 22, 20};
    std::vector<std::pair<double, double>> win;
    for (int r = 0; r < nR; ++r) win.push_back(std::make_pair(rLo[r], rHi[r]));
    // data: one selection per class and per 1 deg OA bin (0-20)
    const int nB = 20;
    std::vector<scan::Sel> sels;
    for (int i = 0; i < nC; ++i) {
        const float a = cLo[i], b = cHi[i];
        sels.push_back([a, b](const scan::Row& r) { return PassPre(r) && r.oa >= a && r.oa < b; });
    }
    for (int i = 0; i < nB; ++i) {
        const float a = i, b = i + 1;
        sels.push_back([a, b](const scan::Row& r) { return PassPre(r) && r.oa >= a && r.oa < b; });
    }
    const std::vector<scan::Spectra> sp = DFillMany(dd, sels, 160, 0.30, 1.10, 2, win);
    P("\n## (d) data (N - CB) / simulation (T1 + T2, mode A, w) per OA(e+e-) class, pre; CB per class in 10 MeV M4 bins; R(pre) = %.4f\n", Rpre);
    std::vector<TGraphErrors*> g(nR), gB(nR), gA(nR);
    for (int r = 0; r < nR; ++r) {
        g[r] = new TGraphErrors(); gB[r] = new TGraphErrors(); gA[r] = new TGraphErrors();
        P("%s\n%-8s %9s %8s %18s %18s %8s %16s %16s %16s %16s\n", rLab[r], "class", "N", "CB", "data [nb]", "sim A [nb]", "N_eff", "data/sim A", "data/(R sim)", "data/sim B", "data/sim (a)");
        for (int i = 0; i < nC; ++i) {
            const float a = cLo[i], b = cHi[i];
            const Sel s = SelAnd(SelPreset(kPre), [a, b](const Row& rr) { return rr.oa >= a && rr.oa < b; });
            const XS x = SimWin(d, s, CompSim(), rLo[r], rHi[r]);
            SimOpt oB; oB.mode = 1;
            const XS xB = SimWin(d, s, CompSim(), rLo[r], rHi[r], oB);
            const XS xA = SimWin(d, s, CompSimVar(1), rLo[r], rHi[r]);
            const scan::WinSum& ws = sp[i].win[2 + r];
            const RatioV q = Div(ws.sig() / kLdata, ws.sigErr() / kLdata, x.s, x.e);
            const RatioV qB = Div(ws.sig() / kLdata, ws.sigErr() / kLdata, xB.s, xB.e);
            const RatioV qA = Div(ws.sig() / kLdata, ws.sigErr() / kLdata, xA.s, xA.e);
            P("%-8s %9.0f %8.0f %9.4f +- %6.4f %9.4f +- %6.4f %8.0f %7.3f +- %6.3f %7.3f +- %6.3f %7.3f +- %6.3f %7.3f +- %6.3f\n", cTxt[i], ws.n[0], ws.cb, ws.sig() / kLdata,
              ws.sigErr() / kLdata, x.s, x.e, x.neff, q.r, q.e, q.r / Rpre, q.e / Rpre, qB.r, qB.e, qA.r, qA.e);
            const double xo = i + 0.5 + (r - 1) * 0.18;
            g[r]->SetPoint(i, xo, q.r); g[r]->SetPointError(i, 0, q.e);
            gB[r]->SetPoint(i, xo + 0.06, qB.r); gB[r]->SetPointError(i, 0, qB.e);
            gA[r]->SetPoint(i, xo + 0.12, qA.r); gA[r]->SetPointError(i, 0, qA.e);
        }
    }
    // against OA in 1 deg bins
    std::vector<TGraphErrors*> gO(nR);
    for (int r = 0; r < nR; ++r) {
        gO[r] = new TGraphErrors();
        Binning bo; bo.nbin = nB; bo.lo = 0; bo.hi = nB;
        const double a = rLo[r], b = rHi[r];
        bo.x = [](const Row& rr) { return (double)rr.oa; };
        const WSpec w = Fill(d, SelAnd(SelPreset(kPre), [a, b](const Row& rr) { return rr.m4 >= a && rr.m4 < b; }), CompSim(), SimOpt(), bo);
        int n = 0;
        for (int i = 0; i < nB; ++i) {
            const scan::WinSum& ws = sp[nC + i].win[2 + r];
            const RatioV q = Div(ws.sig() / kLdata, ws.sigErr() / kLdata, w.s[i] / kLsim, std::sqrt(w.v[i]) / kLsim);
            if (q.r <= 0 || q.e > 0.5 * q.r) continue;
            gO[r]->SetPoint(n, i + 0.5 + (r - 1) * 0.12, q.r); gO[r]->SetPointError(n, 0, q.e); ++n;
        }
    }
    TCanvas* c = new TCanvas("cSumD", "", 2000, 1100);
    Header(c, {{"(d) Data (N - CB) / simulation (signal + physical background) per e^{+}e^{-} opening-angle class, selection \"pre\", absolute scale", 1, true},
               {"data: CB = 2#sqrt{N_{++}N_{--}} of each class in 10 MeV M_{4} bins (no k factor), event-group errors; " + kNoteL, kColSig, false},
               {kNoteS, 1, false},
               {"filled markers: mode A (default); open markers: mode B (HGeant #pi^{0} Dalitz kinematics); diamonds: mode A, variant (a) with ghost-flagged true pairs; dashed line: R(pre)", 1, false},
               {"the small-OA question stays open: the pattern depends on the ghost treatment, the partition mode and the v_{z} slice (tables/review_oa_classes.txt)", kRed + 1, false},
               {kNoteD, kRed + 1, false}},
           0.85);
    c->cd();
    MkPad("pD1", 0.0, 0.0, 0.5, 0.85, false, 0.12, 0.02, 0.04, 0.10);
    TH1D* f1 = new TH1D("fD1", "", nC, 0, nC);
    f1->SetDirectory(nullptr); f1->SetMinimum(0.5); f1->SetMaximum(1.85); f1->SetStats(0);
    for (int i = 0; i < nC; ++i) f1->GetXaxis()->SetBinLabel(i + 1, cLab[i]);
    Ax(f1, "", "data (N - CB) / simulation", 0.85, 0.021, 0.023, YO(0.12, 0.023, 0.85));
    f1->GetXaxis()->SetLabelSize(0.050); f1->GetXaxis()->SetLabelOffset(0.012);
    f1->Draw("AXIS");
    HLine(0, nC, 1.0, kGray + 2, 1, 1);
    HLine(0, nC, Rpre, kGray + 2, 2, 2);
    for (int r = 0; r < nR; ++r) {
        MkG(g[r], rCol[r], rSt[r], 2.0); g[r]->Draw("P SAME");
        MkG(gB[r], rCol[r], rSt[r] + 4, 1.4); gB[r]->Draw("P SAME");
        MkG(gA[r], rCol[r], 33, 1.8); gA[r]->Draw("P SAME");
    }
    TLegend* lg = Leg(0.14, 0.76, 0.97, 0.95, 0.030);
    lg->SetHeader("M_{4} region [GeV/c^{2}]: (filled mode A, open mode B, diamonds (a))");
    for (int r = 0; r < nR; ++r) lg->AddEntry(g[r], rLab[r], "p");
    lg->Draw();
    Txt(0.14, 0.715, Form("dashed: R = %.3f #pm %.3f", Rpre, eRpre), 0.028, kGray + 2);
    c->cd();
    MkPad("pD2", 0.5, 0.0, 1.0, 0.85, false, 0.12, 0.02, 0.04, 0.10);
    TH1D* f2 = Frame("fD2", 0, nB, 0.5, 1.85);
    Ax(f2, "e^{+}e^{-} opening angle OA [deg]", "data (N - CB) / simulation", 0.85, 0.021, 0.023, YO(0.12, 0.023, 0.85));
    f2->Draw("AXIS");
    HLine(0, nB, 1.0, kGray + 2, 1, 1);
    HLine(0, nB, Rpre, kGray + 2, 2, 2);
    VLine(2, 0.5, 1.85, kGray + 1, 3);
    for (int r = 0; r < nR; ++r) { MkG(gO[r], rCol[r], rSt[r], 1.5); gO[r]->Draw("P SAME"); }
    TLegend* lg2 = Leg(0.14, 0.80, 0.97, 0.95, 0.030);
    lg2->SetHeader("1#circ bins, mode A; points with error > 50 % dropped");
    for (int r = 0; r < nR; ++r) lg2->AddEntry(gO[r], rLab[r], "p");
    lg2->Draw();
    Txt(0.14, 0.13, "dotted: OA = 2#circ (usual minimal conversion cut)", 0.030, kGray + 2);
    Save(c, fig + "summary_d_oa_classes");
}

// ============================================================================
// (e) efficiencies of the main cuts: simulation (signal, background) and data
// ============================================================================
template <class R> inline bool PassCut(int id, const R& r) {
    switch (id) {
        case 0: return r.oa > 2.f;
        case 1: return r.oa < 6.f;
        case 2: return r.oa_lab < 50.f;
        case 3: return r.oa_lab < 40.f;
        case 4: return r.mpipi > 0.30f;
        case 5: return r.mpipi > 0.34f;
        case 6: return r.mpipi < 0.41f;
        case 7: return r.oa_eta > 140.f;
        case 8: return r.mm4 < 2.80f;
        case 9: return r.pt4 > 300.f;
        case 10: return r.pip_isUsedVertex == 1 && r.pim_isUsedVertex == 1;
        case 11: return r.pip_isUsedVertex == 1 && r.pim_isUsedVertex == 1 && r.ep_isUsedVertex == 1 && r.em_isUsedVertex == 1;
    }
    return false;
}
const int kNCut = 12;
const char* const kCutLab[kNCut] = {"OA > 2#circ", "OA < 6#circ", "#alpha_{lab} < 50#circ", "#alpha_{lab} < 40#circ", "M_{#pi#pi} > 0.30",
                                    "M_{#pi#pi} > 0.34", "M_{#pi#pi} < 0.41", "#alpha_{#eta} > 140#circ", "MM_{4} < 2.80", "p_{T,4} > 300 *",
                                    "#pi^{+}#pi^{-} in vtx *", "4 tracks in vtx *"};
const char* const kCutTxt[kNCut] = {"OA > 2", "OA < 6", "oa_lab < 50", "oa_lab < 40", "M(pipi) > 0.30", "M(pipi) > 0.34", "M(pipi) < 0.41",
                                    "eta_rest > 140", "MM4 < 2.80", "pT4 > 300 *", "both pions in vtx *", "all four in vtx *"};
const char* const kCutKey[kNCut] = {"oa_lo 2 no cut", "oa_hi 6 no cut", "oa_lab 50 no cut", "oa_lab 40 no cut", "mpipi_lo 0.30 no cut",
                                    "mpipi_lo 0.34 no cut", "mpipi_hi 0.41 no cut", "eta_rest 140 no cut", "mm4_hi 2.80 no cut", "pt4_lo 300 start",
                                    "vtx both pions start", "vtx all four start"};
inline bool CutStart(int id) { return id >= 9; }
struct EV { bool ok = false; double v = 0, e = 0; };
struct TabLine { EV dPeak, sPeak, dTag, sTag, dEdge, sEdge; };
// parse "0.779 +- 0.033    0.564 +- 0.028   -5.0" (data, sim, pull) or "-  -"
void ParseSeg(const std::string& s, EV& dat, EV& sim) {
    std::istringstream is(s);
    std::vector<std::string> t; std::string x;
    while (is >> x) t.push_back(x);
    if (t.size() >= 6 && t[0] != "-") {
        dat.ok = true; dat.v = atof(t[0].c_str()); dat.e = atof(t[2].c_str());
        sim.ok = true; sim.v = atof(t[3].c_str()); sim.e = atof(t[5].c_str());
    }
}
std::string Squeeze(const std::string& s) {
    std::istringstream is(s); std::string x, o;
    while (is >> x) { if (!o.empty()) o += " "; o += x; }
    return o;
}
bool ReadDataEff(const std::string& path, std::vector<TabLine>& out) {
    out.assign(kNCut, TabLine());
    std::ifstream f(path.c_str());
    if (!f) { P("cannot read %s\n", path.c_str()); return false; }
    std::string line;
    int found = 0;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::vector<std::string> seg;
        size_t a = 0, b;
        while ((b = line.find('|', a)) != std::string::npos) { seg.push_back(line.substr(a, b - a)); a = b + 1; }
        seg.push_back(line.substr(a));
        if (seg.size() < 4) continue;
        const std::string key = Squeeze(seg[0]);
        for (int i = 0; i < kNCut; ++i) {
            if (key != kCutKey[i]) continue;
            ParseSeg(seg[1], out[i].dPeak, out[i].sPeak);
            ParseSeg(seg[2], out[i].dTag, out[i].sTag);
            ParseSeg(seg[3], out[i].dEdge, out[i].sEdge);
            ++found;
        }
    }
    P("data-derived efficiencies: %d of %d cut lines found in %s\n", found, kNCut, path.c_str());
    return found == kNCut;
}
void PartE(const SimData& d, const scan::Data& dd, const std::string& fig, const std::string& tab) {
    std::vector<TabLine> T;
    ReadDataEff(tab + "cuts_compare_data.txt", T);
    // simulation bases
    const Sel pre = SelPreset(kPre);
    const Sel pre0 = [](const Row& r) { return r.isBest >= 0 && r.vz > -500.f && r.cor_ok == 1 && PT3(r); };
    const Sel peakStart = [pre0](const Row& r) { return pre0(r) && r.oa_lab < 55.f && r.mpipi > 0.30f; };
    const Sel bumpStart = [pre0](const Row& r) { return pre0(r) && r.oa_lab < 40.f && r.mpipi > 0.30f && r.mpipi < 0.41f; };
    const Comp cPeak = CompPeakConv();
    const Comp cTag = {"bump D+C, photon-tagged", [](const Row& r) { return r.tmpl == 1 && (r.cls == 2 || r.cls == 3) && r.g_idx > 0 && r.m_eeg >= 0.119f && r.m_eeg < 0.173f; }};
    const Comp cSbu = {"bump D+C in the bump window", [](const Row& r) { return r.tmpl == 1 && (r.cls == 2 || r.cls == 3) && r.m4 >= scan::kBumpLo && r.m4 < scan::kBumpHi; }};
    const Comp cBbu = {"T2 in the bump window", [](const Row& r) { return r.tmpl == 2 && r.m4 >= scan::kBumpLo && r.m4 < scan::kBumpHi; }};
    const Comp cBpk = {"T2 in the sim peak window", [](const Row& r) { return r.tmpl == 2 && r.m4 >= kSimPkLo && r.m4 < kSimPkHi; }};
    const Comp cBR = {"T1+T2 in 0.60-0.80", [](const Row& r) { return (r.tmpl == 1 || r.tmpl == 2) && r.m4 >= kRLo && r.m4 < kRHi; }};
    // data background: N - CB in 0.60-0.80 at pre and pre + cut (approximate error, see below)
    std::vector<scan::Sel> ds;
    ds.push_back(DataSelPreset(kPre));
    for (int i = 0; i < kNCut; ++i) ds.push_back([i](const scan::Row& r) { return PassPre(r) && PassCut(i, r); });
    const std::vector<scan::Spectra> dsp = DFillMany(dd, ds, 160, 0.30, 1.10, 2, {{kRLo, kRHi}});
    const double S0 = dsp[0].win[2].sig(), V0 = std::pow(dsp[0].win[2].sigErr(), 2);
    struct Row6 { EffRes peak, tag, sbu, bbu, bpk, bR; EV dR; };
    std::vector<Row6> E(kNCut);
    P("\n## (e) efficiency of each cut alone. Signal: base \"pre\" (cuts marked * : the scan_extra start selections on pre0, as the data);\n");
    P("#     background: base \"pre\" for every cut. Simulation: mode A, w, event-group errors (pass/fail correlated).\n");
    P("#     data: peak eps_S, tagged eta eps_T, edge-step eps_E read from tables/cuts_compare_data.txt (data study; data peak contains the direct decay);\n");
    P("#     data background: (N - CB) kept in 0.60-0.80 (no eta signal there), error approximate: [(1-e)^2 V_pass + e^2 (V_base - V_pass)]^0.5 / S_base\n");
    P("%-20s | %-13s %-13s %-8s | %-13s %-13s %-8s | %-13s %-13s %-8s | %-13s %-13s %-13s %-13s\n", "cut", "data peak", "sim peak", "check",
      "data tagged", "sim tagged", "check", "data edge E", "sim S_bump", "check", "sim B_bump", "sim B_peak", "sim T12 R-win", "data N-CB R-win");
    double worst = 0;
    for (int i = 0; i < kNCut; ++i) {
        const Sel pass = [i](const Row& r) { return PassCut(i, r); };
        const Sel bPk = CutStart(i) ? peakStart : pre;
        const Sel bBu = CutStart(i) ? bumpStart : pre;
        E[i].peak = Efficiency(d, bPk, pass, cPeak);
        E[i].tag = Efficiency(d, bBu, pass, cTag);
        E[i].sbu = Efficiency(d, pre, pass, cSbu);
        E[i].bbu = Efficiency(d, pre, pass, cBbu);
        E[i].bpk = Efficiency(d, pre, pass, cBpk);
        E[i].bR = Efficiency(d, pre, pass, cBR);
        const double S1 = dsp[1 + i].win[2].sig(), V1 = std::pow(dsp[1 + i].win[2].sigErr(), 2);
        const double e = S1 / S0;
        E[i].dR.ok = true; E[i].dR.v = e; E[i].dR.e = std::sqrt(std::max(0., (1 - e) * (1 - e) * V1 + e * e * std::max(0., V0 - V1))) / S0;
        auto ck = [&](const EV& tabSim, const EffRes& mine, bool use) -> std::string {
            if (!tabSim.ok || !use) return "-";
            const double dlt = mine.eps - tabSim.v;
            worst = std::max(worst, std::fabs(dlt));
            return Form("%+.3f", dlt);
        };
        auto fe = [](const EV& x) -> std::string { return x.ok ? Form("%.3f+-%.3f", x.v, x.e) : std::string("-"); };
        auto fe4 = [](const EV& x) -> std::string { return x.ok ? Form("%.4f+-%.4f", x.v, x.e) : std::string("-"); };
        auto fr = [](const EffRes& x) -> std::string { return Form("%.3f+-%.3f", x.eps, x.err); };
        P("%-20s | %-13s %-13s %-8s | %-13s %-13s %-8s | %-13s %-13s %-8s | %-13s %-13s %-13s %-13s\n", kCutTxt[i], fe(T[i].dPeak).c_str(), fr(E[i].peak).c_str(),
          ck(T[i].sPeak, E[i].peak, true).c_str(), fe(T[i].dTag).c_str(), fr(E[i].tag).c_str(), ck(T[i].sTag, E[i].tag, true).c_str(), fe(T[i].dEdge).c_str(),
          fr(E[i].sbu).c_str(), ck(T[i].sEdge, E[i].sbu, !CutStart(i)).c_str(), fr(E[i].bbu).c_str(), fr(E[i].bpk).c_str(), fr(E[i].bR).c_str(), fe4(E[i].dR).c_str());
    }
    P("check = this macro minus the simulated value of tables/cuts_compare_data.txt; largest |difference| %.4f\n", worst);
    P("N_eff(base, pre): peak conv %.0f, tagged %.0f, S_bump %.0f, B_bump %.0f, B_peak %.0f, T12 R-window %.0f\n",
      Efficiency(d, pre, [](const Row&) { return true; }, cPeak).NeffAll, Efficiency(d, pre, [](const Row&) { return true; }, cTag).NeffAll,
      Efficiency(d, pre, [](const Row&) { return true; }, cSbu).NeffAll, Efficiency(d, pre, [](const Row&) { return true; }, cBbu).NeffAll,
      Efficiency(d, pre, [](const Row&) { return true; }, cBpk).NeffAll, Efficiency(d, pre, [](const Row&) { return true; }, cBR).NeffAll);

    // ---- figure
    TCanvas* c = new TCanvas("cSumE", "", 2200, 1900);
    Header(c, {{"(e) Efficiency of each cut alone: simulation (open markers) next to the data-derived values (blue, full markers)", 1, true},
               {"base: selection \"pre\" (isBest #geq 0, PT3, v_{z} > -160 mm, corrected momenta); efficiencies do not depend on L_{data}", 1, false},
               {"* signal base = start selection of the data study: pre without v_{z} cut + #alpha_{lab} < 55#circ, M_{#pi#pi} > 0.30 (peak) or #alpha_{lab} < 40#circ, 0.30 < M_{#pi#pi} < 0.41 (bump)", 1, false},
               {"data: fits of the data study (tables/cuts_compare_data.txt); the data #eta peak contains #eta #rightarrow #pi^{+}#pi^{-}e^{+}e^{-}, which the simulation lacks", kRed + 1, false},
               {"simulation: SMASH std + lep, mode A, w (PDG BR), event-group errors; simulated photon tag 0.119 #leq M_{ee#gamma} < 0.173 emulates the ECAL tag (UNVERIFIED)", 1, false},
               {"background: physical background (sim) in the bump and peak windows; sim T1 + T2 and data N - CB in M_{4} 0.60 - 0.80 (no #eta there; data error approximate)", 1, false}},
           0.885);
    const double pY[4] = {0.885, 0.63, 0.375, 0.0};
    const char* pT[3] = {"#eta peak (#eta #rightarrow #pi^{+}#pi^{-}#gamma, converted photon)", "#eta bump (#eta #rightarrow #pi^{+}#pi^{-}#pi^{0})",
                         "background (no #eta)"};
    for (int pnl = 0; pnl < 3; ++pnl) {
        c->cd();
        const double bm = pnl == 2 ? 0.23 : 0.03;
        MkPad(Form("pE%d", pnl), 0.0, pY[pnl + 1], 1.0, pY[pnl], false, 0.07, 0.20, 0.07, bm);
        const double H = pY[pnl] - pY[pnl + 1];
        TH1D* f = new TH1D(Form("fE%d", pnl), "", kNCut, 0, kNCut);
        f->SetDirectory(nullptr); f->SetMinimum(pnl == 2 ? 0.0 : 0.0); f->SetMaximum(1.12); f->SetStats(0);
        for (int i = 0; i < kNCut; ++i) f->GetXaxis()->SetBinLabel(i + 1, pnl == 2 ? kCutLab[i] : "");
        Ax(f, "", "efficiency", H, 0.0170, 0.0190, YO(0.07, 0.0190, H));
        f->GetXaxis()->SetLabelSize(0.0155 / H); f->GetXaxis()->SetLabelOffset(0.012); f->GetXaxis()->LabelsOption("h");
        f->Draw("AXIS");
        HLine(0, kNCut, 1.0, kGray + 1, 2, 1);
        for (int i = 1; i < kNCut; ++i) VLine(i, 0, 1.12, kGray, 3);
        Txt(0.075, 0.90, pT[pnl], 0.0165 / H, 1, 13, 62);
        TLegend* lg = Leg(0.805, pnl == 2 ? 0.25 : 0.06, 0.995, 0.93, 0.0118 / H);
        auto series = [&](std::function<EV(int)> get, double dx, int col, int st, double sz, const char* lab) {
            TGraphErrors* g = new TGraphErrors();
            int n = 0;
            for (int i = 0; i < kNCut; ++i) {
                const EV v = get(i);
                if (!v.ok) continue;
                g->SetPoint(n, i + 0.5 + dx, v.v); g->SetPointError(n, 0, v.e); ++n;
            }
            MkG(g, col, st, sz); g->Draw("P SAME");
            lg->AddEntry(g, lab, "p");
        };
        auto fromE = [](const EffRes& x) { EV v; v.ok = true; v.v = x.eps; v.e = x.err; return v; };
        if (pnl == 0) {
            series([&](int i) { return T[i].dPeak; }, -0.12, kColSig, 20, 2.0, "data: #eta peak #varepsilon_{S} (fit)");
            series([&](int i) { return fromE(E[i].peak); }, 0.12, kColPeakConv, 24, 2.0, "sim: peak conversion");
        } else if (pnl == 1) {
            series([&](int i) { return T[i].dTag; }, -0.27, kColSig, 21, 1.9, "data: tagged #eta #varepsilon_{T}");
            series([&](int i) { return fromE(E[i].tag); }, -0.09, kColBumpDalitz, 25, 1.9, "sim: tagged bump (D + C)");
            series([&](int i) { return T[i].dEdge; }, 0.09, kColSig, 22, 2.1, "data: edge step #varepsilon_{E}");
            series([&](int i) { EV v = fromE(E[i].sbu); v.ok = !CutStart(i); return v; }, 0.27, kColBumpConv, 26, 2.1, "sim: bump D+C, bump window");
        } else {
            series([&](int i) { return E[i].dR; }, -0.27, kColSig, 33, 2.5, "data: N-CB, M_{4} 0.60-0.80");
            series([&](int i) { return fromE(E[i].bR); }, -0.09, kColSim, 27, 2.5, "sim: T1+T2, M_{4} 0.60-0.80");
            series([&](int i) { return fromE(E[i].bbu); }, 0.09, kGray + 2, 25, 1.9, "sim: T2, bump window");
            series([&](int i) { return fromE(E[i].bpk); }, 0.27, kGray + 2, 24, 1.9, "sim: T2, peak window");
        }
        lg->Draw();
    }
    Save(c, fig + "summary_e_cut_eff");
}

// ============================================================================
// (f) PRELIMINARY template fit of "pre" (machinery of sim_templatefit.C)
// ============================================================================
void PartF(const SimData& d, const scan::Data& dd, const std::string& fig) {
    stf::TplSet T;
    stf::FillTpl(d, SelPreset(kPre), stf::OptOf(0), 1.0, T);
    stf::DSet D;
    stf::FillData(dd, DataSelPreset(kPre), D);
    // width of the peak-conversion template (iterated +-2 sigma Gaussian fit, as sim_templatefit.C)
    double gs = 0.006, gmSim = 0;
    {
        TH1D* h = new TH1D("hF_peak", "", stf::kNB, stf::kLo, stf::kHi);
        h->SetDirectory(nullptr);
        for (int i = 0; i < stf::kNB; ++i) { h->SetBinContent(i + 1, T.s[stf::cPeak][i]); h->SetBinError(i + 1, std::sqrt(T.C[stf::cPeak][stf::cPeak][i * stf::kNB + i])); }
        TF1* g = new TF1("gF_peak", "gaus", 0.53, 0.56);
        g->SetParameters(h->GetMaximum(), 0.5445, 0.006);
        h->Fit(g, "QNR0");
        for (int it = 0; it < 3; ++it) {
            g->SetRange(g->GetParameter(1) - 2 * g->GetParameter(2), g->GetParameter(1) + 2 * g->GetParameter(2));
            h->Fit(g, "QNR0");
        }
        gmSim = g->GetParameter(1); gs = std::fabs(g->GetParameter(2));
        P("\n## (f) PRELIMINARY template fit, pre, 2.5 MeV, 0.43-0.75, mode A, w\npeak-conversion template: Gaussian mean %.5f +- %.5f, sigma %.5f +- %.5f GeV/c2 (data 0.5439, 0.0054)\n",
          gmSim, g->GetParError(1), gs, g->GetParError(2));
    }
    const double gm = scan::kPeakMean;
    const stf::FitIn in = stf::MakeIn(D, T, 1, 0.43, 0.75, gm, gs);
    const std::vector<stf::Model> M = stf::Models();
    const int use[3] = {0, 5, 6};   // M1, M6, M7
    stf::FitOut o[3];
    for (int j = 0; j < 3; ++j) {
        o[j] = stf::DoFit(in, M[use[j]]);
        stf::Implied(o[j], M[use[j]], T, D, gm, gs);
        P("%s (%s): %s chi2/ndf %.1f/%d, template share of V %.2f; implied S = %.0f +- %.0f +- %.0f, E = %.0f +- %.0f +- %.0f (data rows)\n", M[use[j]].tag.c_str(),
          M[use[j]].label.c_str(), stf::ParStr(M[use[j]], o[j], false).c_str(), o[j].chi2, o[j].ndf, o[j].tvShare, o[j].S, o[j].eS, o[j].tS, o[j].E, o[j].eE, o[j].tE);
    }
    P("data: S (P_main) = %.0f +- %.0f +- %.0f (model spread), E = %.0f +- %.0f (data rows)\n", D.pf.S, D.pf.eS, D.spr.half, D.E, D.eE);

    TCanvas* c = new TCanvas("cSumF", "", 2000, 1500);
    Header(c, {{"(f) PRELIMINARY template fit of the data N - CB in M_{4}, selection \"pre\", 2.5 MeV bins, fit range 0.43 - 0.75 GeV/c^{2}", 1, true},
               {"model: #alpha_{bump} #times bump (#pi^{0} Dalitz + conversion) + #alpha_{peak} #times peak (conversion) + #alpha_{bkg} #times physical background", 1, false},
               {"M1: as written; M6: background coefficient #alpha_{bkg} + #beta (M - 0.6)/0.1 (linear tilt); #alpha = 1: SMASH as simulated (std + lep, mode A, w)", 1, false},
               {"#chi^{2}: data covariance (event groups, CB) + template covariance (event groups); data: blue N - CB, red CB; " + kNoteL, 1, false},
               {kNoteD + "; #alpha_{peak} absorbs it", kRed + 1, false}},
           0.895);
    std::vector<double> yD = D.sp.sig;
    TH1D* hD = stf::HFine("hF_data", yD, &D.cov, 1);
    Mk(hD, kColSig, 20, 1.0);
    TH1D* hCB = new TH1D("hF_cb", "", stf::kNB, stf::kLo, stf::kHi);
    hCB->SetDirectory(nullptr);
    for (int i = 0; i < stf::kNB; ++i) { hCB->SetBinContent(i + 1, D.sp.cb[i] / (kLdata * stf::kBW)); hCB->SetBinError(i + 1, std::sqrt(D.sp.cbVar[i]) / (kLdata * stf::kBW)); }
    Mk(hCB, kColCB, 21, 0.8);
    stf::Coef one; for (int k = 0; k < stf::kNC; ++k) one.c[k] = 1.;
    TH1D* hOne = stf::HFine("hF_one", stf::ModelFine(T, one, gm, gs), nullptr, 1);
    Ln(hOne, kColSim, 3, 2);
    double ymax = 0;
    for (int i = 1; i <= hD->GetNbinsX(); ++i) if (hD->GetBinCenter(i) > 0.43 && hD->GetBinCenter(i) < 0.75) ymax = std::max(ymax, hD->GetBinContent(i));
    for (int j = 0; j < 2; ++j) {
        const stf::Model& m = M[use[j]];
        const stf::FitOut& r = o[j];
        const double x0 = 0.5 * j, x1 = x0 + 0.5;
        stf::Coef cB = r.q; cB.c[stf::cBump] = 0; cB.c[stf::cPeak] = 0; cB.AG = 0;
        stf::Coef cBb = r.q; cBb.c[stf::cPeak] = 0; cBb.AG = 0;
        stf::Coef cAll = r.q;
        TH1D* h1 = stf::HFine(Form("hF_b%d", j), stf::ModelFine(T, cB, gm, gs), nullptr, 1);
        TH1D* h2 = stf::HFine(Form("hF_bb%d", j), stf::ModelFine(T, cBb, gm, gs), nullptr, 1);
        TH1D* h3 = stf::HFine(Form("hF_all%d", j), stf::ModelFine(T, cAll, gm, gs), nullptr, 1);
        for (TH1D* h : {h1, h2, h3})
            for (int i = 1; i <= h->GetNbinsX(); ++i) if (h->GetBinCenter(i) < 0.43 || h->GetBinCenter(i) > 0.75) h->SetBinContent(i, 0);
        c->cd();
        MkPad(Form("pF%d", j), x0, 0.30, x1, 0.895, false, 0.16, 0.02, 0.04, 0.0);
        const double H = 0.595;
        TH1D* f = Frame(Form("fF%d", j), 0.43, 0.75, 0, 1.25 * ymax);
        Ax(f, kXt, kYt, H, 0.0175, 0.0190, YO(0.16, 0.0190, H));
        f->Draw("AXIS");
        Ln(h3, kColPeakConv, 1, 1, kColPeakConv); h3->Draw("HIST SAME");
        Ln(h2, kColBumpDalitz, 1, 1, kColBumpDalitz); h2->Draw("HIST SAME");
        Ln(h1, kColPhysBkg, 1, 1, kColPhysFill); h1->Draw("HIST SAME");
        TH1D* h3l = (TH1D*)h3->Clone(Form("hF_alll%d", j)); Ln(h3l, 1, 2, 1); h3l->Draw("HIST SAME");
        TH1D* hO = (TH1D*)hOne->Clone(Form("hF_one%d", j)); hO->Draw("HIST SAME");
        hD->Draw("P SAME"); hCB->Draw("P SAME");
        gPad->RedrawAxis();
        Txt(0.18, 0.93, Form("%s: %s", m.tag.c_str(), j == 0 ? "three templates" : "three templates, background tilt"), 0.034, 1, 13, 62);
        std::string ps = stf::ParStr(m, r, true);
        Txt(0.18, 0.875, ps.c_str(), 0.026, 1, 13);
        Txt(0.18, 0.835, Form("#chi^{2}/ndf = %.0f/%d; implied peak S = %.0f #pm %.0f #pm %.0f (data %.0f #pm %.0f)", r.chi2, r.ndf, r.S, r.eS, r.tS, D.pf.S, D.pf.eS), 0.026, 1, 13);
        TLegend* lg = Leg(0.18, 0.50, 0.70, 0.80, 0.028);
        lg->AddEntry(hD, "data: N - CB", "p");
        lg->AddEntry(hCB, "data: CB", "p");
        lg->AddEntry(h3l, "fit: sum", "l");
        lg->AddEntry(h3, "fit: + #alpha_{peak} #times peak conversion", "f");
        lg->AddEntry(h2, "fit: + #alpha_{bump} #times bump", "f");
        lg->AddEntry(h1, "fit: #alpha_{bkg} #times physical background", "f");
        lg->AddEntry(hO, "SMASH alone (all #alpha = 1)", "l");
        lg->Draw();
        // pulls
        c->cd();
        MkPad(Form("pFp%d", j), x0, 0.0, x1, 0.30, false, 0.16, 0.02, 0.0, 0.25);
        const double H2 = 0.30;
        TH1D* fp = Frame(Form("fFp%d", j), 0.43, 0.75, -5.5, 5.5);
        Ax(fp, kXt, "pull", H2, 0.0175, 0.0190, YO(0.16, 0.0190, H2));
        fp->Draw("AXIS");
        HLine(0.43, 0.75, 0, kGray + 2, 1, 1); HLine(0.43, 0.75, 2, kGray + 1, 2, 1); HLine(0.43, 0.75, -2, kGray + 1, 2, 1);
        for (double xb : {scan::kBumpLo, scan::kBumpHi, scan::kPeakLo, scan::kPeakHi}) VLine(xb, -5.5, 5.5, kGray + 1, 3);
        TGraph* gp = new TGraph();
        for (int q = 0; q < in.npts; ++q) gp->SetPoint(q, in.x[q], r.pull[q]);
        MkG(gp, kColSig, 20, 1.0); gp->Draw("P SAME");
        Txt(0.18, 0.95, "pull = (data - fit) / #sigma, #sigma from data and template covariance; dotted: bump and peak windows", 0.050, 1, 13);
    }
    Save(c, fig + "summary_f_templatefit_pre");
}

}  // namespace ssum

// ============================================================================
void sim_summary(int parts = 31) {
    using namespace ssum;
    TStopwatch sw;
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0); gStyle->SetEndErrorSize(0); gStyle->SetLegendFont(42);
    gStyle->SetTextFont(42); gStyle->SetLabelFont(42, "XYZ"); gStyle->SetTitleFont(42, "XYZ");
    const std::string fig = std::string(simeng::kFigDir) + "/";
    const std::string tab = std::string(simeng::kTabDir) + "/";
    gOut = fopen((tab + (parts == 31 ? "summary_recomputed.txt" : Form("summary_recomputed_part%d.txt", parts))).c_str(), "w");
    P("# summary_recomputed.txt -- plots/eta_sim/sim_summary.C (every number below is recomputed by this macro)\n");
    P("# simulation: SMASH std + lep, rows owned in mode A, weight w (PDG BR), L_sim = 221.05 nb^-1, event-group errors, double precision\n");
    P("# data: (N - CB) / L_data, L_data = 6080.5 nb^-1 (UNCONFIRMED), CB = 2 sqrt(N++ N--) in 10 MeV M4 bins, no k; ScanEngine event-group errors\n");
    P("# eta -> pi+pi-e+e- (direct decay) is NOT in the simulation\n");

    simeng::SimData d;
    if (!simeng::Load(d, simeng::kCompact, (1 << simeng::kStdUS) | (1 << simeng::kLepUS))) return;
    // event groups on the full trees (all rows, owned or not): one isBest == 1 row per group (adjacent rows with equal grp)
    {
        TFile* f = TFile::Open(simeng::kCompact, "READ");
        for (int k : {simeng::kStdUS, simeng::kLepUS}) {
            TTree* t = (TTree*)f->Get(simeng::TreeName(k));
            Int_t grp = 0, isBest = 0;
            t->SetBranchStatus("*", 0); t->SetBranchStatus("grp", 1); t->SetBranchStatus("isBest", 1);
            t->SetBranchAddress("grp", &grp); t->SetBranchAddress("isBest", &isBest);
            long ng = 0, bad = 0; int cur = -1, nb = 0;
            const Long64_t n = t->GetEntries();
            for (Long64_t i = 0; i < n; ++i) {
                t->GetEntry(i);
                if (grp != cur) { if (cur >= 0 && nb != 1) ++bad; cur = grp; nb = 0; ++ng; }
                if (isBest == 1) ++nb;
            }
            if (nb != 1) ++bad;
            P("event groups %s (all %lld rows): %ld groups, %ld without exactly one isBest == 1 row\n", simeng::TreeName(k), n, ng, bad);
            t->ResetBranchAddresses();
        }
        f->Close(); delete f;
    }
    scan::Data dd;
    if (parts & 30) if (!scan::Load(dd, scan::kCompact, 0.30, 1.10)) return;

    if (parts & 1) PartA(d, fig);
    double Rpre = 1.2514, eRpre = 0.0051;
    if (parts & 2) {
        PartBC(d, dd, fig);
    }
    if (parts & 4) {
        // R(pre) recomputed (the same as in PartBC)
        const std::vector<scan::Spectra> s = DFillMany(dd, {simeng::DataSelPreset(simeng::kPre)}, 160, 0.30, 1.10, 2, {{kRLo, kRHi}});
        const XS x = SimWin(d, simeng::SelPreset(simeng::kPre), simeng::CompSim(), kRLo, kRHi);
        const RatioV q = Div(s[0].win[2].sig() / simeng::kLdata, s[0].win[2].sigErr() / simeng::kLdata, x.s, x.e);
        Rpre = q.r; eRpre = q.e;
        PartD(d, dd, fig, Rpre, eRpre);
    }
    if (parts & 8) PartE(d, dd, fig, tab);
    if (parts & 16) PartF(d, dd, fig);
    P("\n# run time %.0f s\n", sw.RealTime());
    if (gOut) fclose(gOut);
    gOut = nullptr;
}
