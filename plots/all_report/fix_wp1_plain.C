// ============================================================================
// fix_wp1_plain.C -- COPY of fix_wp1.C (2026-11-16, short-report figure agent for report_all_short.tex). Writes only
// docs/latex/figs_all/fix_B1_plain.{pdf,png} (= fix_B1 of fix_wp1.C) and wp1_B2_plain.{pdf,png} (= wp1_B2 of wp1_before.C, whose
// FigM4 is identical to this one except the CB legend, here the corrected 2 (N++ N--)^{1/2}); the map fix_B3 and its statistics
// (needed first: FigM4 alone crashes) and the intermediate EPS files go to the scratch directory kScr; nothing is removed.
// Only label strings in plain words (readability instruction I2): S0 / pre / isBest / cor_ok -> nominal / base selection / all
// combinations / corr. momenta; SMASH T1+T2 -> SMASH total; G_pi0D, G_conv, G_vec + G_oth, G_eta, C2-C5 -> physics names; the
// window letters B, P -> "bump", "peak"; w(V3) -> eta-decay weights (centred on the window; in the log figure "peak" is drawn a little higher so that the two
// words do not touch); "(a)" ghost variant -> "SMASH + ghost-flagged pairs" / "(SMASH + ghosts)". Same input hists_v1.root, same drawing code, numbers unchanged.
// Run: root -l -b -q -e 'gSystem->SetBuildDir("<scratch>/aclic_wp1_plain", kTRUE)' 'plots/all_report/fix_wp1_plain.C+O'
// ============================================================================
// fix_wp1.C -- review fix (prefix fix): Figs. B1 and B3 of wp1_before.C (WP1) redrawn with (i) the CB legend written as
// 2 (N++ N--)^{1/2} (the TLatex square root was broken) and (ii) the window boxes of B3 in the common palette (P red, CR-U green,
// as Figs. 1 and C1). Everything else is wp1_before.C unchanged (namespace renamed). Reads hists_v1.root only.
// Output: docs/latex/figs_all/fix_B1.{pdf,png}, fix_B3.{pdf,png}; the map statistics go to /hdd3/.../all_report/fix/fix_wp1_map.txt
// Run: root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_fix", kTRUE)' 'plots/all_report/fix_wp1.C+O'
// ============================================================================
// ============================================================================
// wp1_before.C  --  WP1 "before": data against SMASH in absolute units, every
//                   channel at its own cross section, NO renormalisation (prefix wp1)
// ============================================================================
// SPEC report_all, WP1. Input: hists_v1.root (WP0) only.
//   data  : (N - CB) x 64/63 / L, L = 6079.2 nb^-1, CB = 2 sqrt(N++ N--) (k = 1, spread rule of WP0),
//           errors: N event-group, CB 200 replicas (where stored; B_map and B_<var> only linear)
//   sim   : SMASH std + lep + eta, mode A, w x V3 (code 108), T1 + T2, M4' = 0.99615 M4, M(pipi) x 0.9952,
//           sigma = sum w / 221.05 nb^-1, event-group errors; ghost variant (a) dashed; weights off as a column
// Outputs
//   docs/latex/figs_all/wp1_B1.{pdf,png}            M4 0.30-1.10, log, S0 | pre, stacked groups, ratio pads
//   docs/latex/figs_all/wp1_B2.{pdf,png}            M4 0.40-0.70, linear, 2.5 MeV, S0 | pre
//   docs/latex/figs_all/wp1_B3.{pdf,png}            data/sim map M4 x M(pipi) with E(m,0), E(m,60 MeV), peak band, CR boxes
//   docs/latex/figs_all/wp1_B4_<var>[_pre].{pdf,png} small multiples of M(pipi), M(ee), OA, vz, MM4, pT4, y4 in the 8 regions
//   docs/latex/figs_all/wp1_K0S.{pdf,png}           K0S peak in M(pipi), data and simulation (S0, pre)
//   docs/latex/tables_all/wp1_B1.{txt,tex}          Table B1 (+ sigma per group, isBest == 1, K0S, shape chi2)
//   docs/latex/tables_all/wp1_A2.{txt,tex}          Table A2: corrections in plain words with their procedure
// Run (from /home/przygoda/HADES/FAT/FAT; < 1 min, < 2 GB):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_wp1", kTRUE)' 'plots/all_report/wp1_before.C+O'
// Blinding (D19): data N of B_mee_R2..R4 in 12-28 MeV is masked in hists_v1.root and not drawn; in addition this macro does
// not draw the data OA of R2..R4 in 4-30 deg (an X17 near 17 MeV sits there; v0 recommendation). D20: the window contents
// below are N - CB including all background, not signal yields.
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "TBox.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"

#include "defs.h"
#include "hists_io.h"
#include "s5_blind_ana_fit.h"

namespace wp1fixplain {
using namespace alldefs;
const char* const kFig = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_all";
const char* const kTab = "/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all";
TFile* gF = nullptr;
int gUid = 0;
const char* const kScr = "/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ6/figplain_allshort/out/wp1";
// _plain: as the SaveFig of s5_blind_ana_fit.h (EPS -> ps2pdf -dEPSCrop -> PDF, then PNG), but the intermediate EPS is written
// to kScr and kept (the library function removes it with rm, which this copy must not do)
void SaveNoRm(TCanvas* c, const std::string& base) {
    const std::string eps = std::string(kScr) + "/" + gSystem->BaseName(base.c_str()) + ".eps";
    c->SaveAs(eps.c_str());
    gSystem->Exec(Form("ps2pdf -dEPSCrop %s %s.pdf", eps.c_str(), base.c_str()));
    c->SaveAs((base + ".png").c_str());
}
std::string U(const char* b) { return Form("%s_%d", b, gUid++); }

struct VE { double v = 0, e = 0; };
VE Div(VE a, VE b) {
    VE r;
    if (b.v == 0) return r;
    r.v = a.v / b.v;
    r.e = std::sqrt(std::pow(a.e / b.v, 2) + std::pow(a.v * b.e / (b.v * b.v), 2));
    return r;
}
double RepRms(const TH1* rep, int bx) {
    const int n = rep->GetNbinsY();
    double m = 0, m2 = 0;
    for (int r = 1; r <= n; ++r) { const double x = rep->GetBinContent(bx, r); m += x; m2 += x * x; }
    m /= n;
    return std::sqrt(std::max(0., m2 / n - m * m));
}

// ---------------------------------------------------------------- data and simulation access (nb units)
// what: 0 = N, 1 = CB, 2 = N - CB. CB error = replica RMS when the layout has replicas (1D only), else the stored linear error.
TH1* DataNb(const char* sel, const std::string& lay, int what) {
    TH1* n = allio::Data(gF, sel, "all", "N", lay);
    TH1* c = allio::Data(gF, sel, "all", "CB", lay);
    if (!n || !c) return nullptr;
    if (n->GetDimension() == 1) {
        TH1* r = (TH1*)gF->Get(Form("data/%s/all/CBrep_%s", sel, lay.c_str()));
        if (r) for (int b = 1; b <= c->GetNbinsX(); ++b) c->SetBinError(b, RepRms(r, b));
    }
    TH1* o = (TH1*)n->Clone(U("dnb").c_str());
    o->SetDirectory(nullptr);
    for (int b = 0; b < n->GetNcells(); ++b) {
        double v = 0, e = 0;
        if (what == 0) { v = n->GetBinContent(b); e = n->GetBinError(b); }
        else if (what == 1) { v = c->GetBinContent(b); e = c->GetBinError(b); }
        else { v = n->GetBinContent(b) - c->GetBinContent(b); e = std::hypot(n->GetBinError(b), c->GetBinError(b)); }
        o->SetBinContent(b, v * kRawToNb);
        o->SetBinError(b, e * kRawToNb);
    }
    delete n; delete c;
    return o;
}
TH1* SimSum(const char* sel, const char* var, const std::vector<int>& comps, const std::string& lay) {
    TH1* s = nullptr;
    for (int c : comps) {
        TH1* h = allio::Sim(gF, sel, var, c, lay);
        if (!h) continue;
        if (!s) { s = (TH1*)h->Clone(U("ss").c_str()); s->SetDirectory(nullptr); }
        else s->Add(h);
        delete h;
    }
    return s;
}
VE DataWin(const char* sel, int wid, int what = 2) {
    TH1* n = allio::Data(gF, sel, "all", "N", "win");
    TH1* c = allio::Data(gF, sel, "all", "CB", "win");
    TH1* r = allio::Data(gF, sel, "all", "CBrep", "win");
    const double ce = r ? RepRms(r, wid + 1) : allio::WinErr(c, wid);
    VE o;
    if (what == 0) { o.v = allio::Win(n, wid); o.e = allio::WinErr(n, wid); }
    else if (what == 1) { o.v = allio::Win(c, wid); o.e = ce; }
    else { o.v = allio::Win(n, wid) - allio::Win(c, wid); o.e = std::hypot(allio::WinErr(n, wid), ce); }
    delete n; delete c; delete r;
    return o;
}
VE Nb(VE raw) { return {raw.v * kRawToNb, raw.e * kRawToNb}; }
VE SimWin(const char* sel, const char* var, int comp, int wid) {
    TH1* h = allio::Sim(gF, sel, var, comp, "win");
    VE o;
    if (h) { o.v = allio::Win(h, wid); o.e = allio::WinErr(h, wid); }
    delete h;
    return o;
}
VE SimWinSum(const char* sel, const char* var, const std::vector<int>& comps, int wid) {
    VE o;
    for (int c : comps) { const VE x = SimWin(sel, var, c, wid); o.v += x.v; o.e = std::hypot(o.e, x.e); }
    return o;
}

// ---------------------------------------------------------------- shape test: one free scale a, chi2 = sum (d - a s)^2 / (vd + a^2 vs)
struct Shape { double a = 0, ea = 0, chi2 = 0; int ndf = 0; };
Shape ShapeFit(const char* sel, const char* var, const std::vector<std::string>& lays, double lo, double hi) {
    std::vector<double> d, vd, s, vs;
    for (const auto& l : lays) {
        TH1* D = DataNb(sel, l, 2);
        TH1* S = allio::Sim(gF, sel, var, kT12, l);
        if (!D || !S) { delete D; delete S; continue; }
        for (int b = 1; b <= D->GetNbinsX(); ++b) {
            const double x = D->GetXaxis()->GetBinCenter(b);
            if (x < lo || x > hi) continue;
            if (S->GetBinContent(b) <= 0) continue;
            d.push_back(D->GetBinContent(b)); vd.push_back(std::pow(D->GetBinError(b), 2));
            s.push_back(S->GetBinContent(b)); vs.push_back(std::pow(S->GetBinError(b), 2));
        }
        delete D; delete S;
    }
    Shape r;
    double sd = 0, ss = 0;
    for (size_t i = 0; i < d.size(); ++i) { sd += d[i]; ss += s[i]; }
    double a = ss > 0 ? sd / ss : 1.;
    double den = 0;
    for (int it = 0; it < 100; ++it) {
        double num = 0; den = 0;
        for (size_t i = 0; i < d.size(); ++i) { const double w = vd[i] + a * a * vs[i]; num += d[i] * s[i] / w; den += s[i] * s[i] / w; }
        const double an = num / den;
        if (std::fabs(an - a) < 1e-12) { a = an; break; }
        a = an;
    }
    r.a = a; r.ea = 1. / std::sqrt(den);
    for (size_t i = 0; i < d.size(); ++i) r.chi2 += std::pow(d[i] - a * s[i], 2) / (vd[i] + a * a * vs[i]);
    r.ndf = (int)d.size() - 1;
    return r;
}

// ---------------------------------------------------------------- simulated stacks
struct Stk { std::vector<int> comps; const char* label; int col; };
std::vector<Stk> StackDef(bool merged) {
    if (merged)
        return {{{kGpi0D, kGconv, kGvec, kGoth}, "multi-pion G_{#pi^{0}D}+G_{conv}+G_{vec}+G_{oth}", kGray},
                {{kGeta}, "G_{#eta} (other #eta decays)", kViolet - 9},
                {{kC2, kC3}, "bump C_{2}+C_{3}", kOrange + 7},
                {{kC4, kC5}, "peak C_{4}+C_{5}", kGreen + 2}};
    return {{{kGpi0D}, "non-#eta #pi^{0} Dalitz pair + pions", kGray},
            {{kGconv}, "non-#eta conversion pair + pions", kGray + 2},
            {{kGvec, kGoth}, "vector mesons + other pairs", kSpring - 9},
            {{kGeta}, "other #eta decays", kViolet - 9},
            {{kC2}, "bump, #pi^{0} Dalitz", kOrange + 7},
            {{kC3}, "bump, conversion", kMagenta + 1},
            {{kC4}, "peak, #eta#rightarrow#pi^{+}#pi^{-}#gamma, conversion", kGreen + 2},
            {{kC5}, "peak, #eta#rightarrow#pi^{+}#pi^{-}e^{+}e^{-}", kCyan + 1}};
}
const int kColA = kTeal + 3;   // ghost variant (a)

void Ax(TH1* h, int lab, int tit, double xoff, double yoff) {
    TAxis* ax[3] = {h->GetXaxis(), h->GetYaxis(), h->GetZaxis()};
    for (TAxis* a : ax) { a->SetLabelFont(43); a->SetLabelSize(lab); a->SetTitleFont(43); a->SetTitleSize(tit); }
    h->GetXaxis()->SetTitleOffset(xoff);
    h->GetYaxis()->SetTitleOffset(yoff);
}
void Mark(TH1* h, int col, int sty, double size = 1.0) {
    h->SetMarkerColor(col); h->SetLineColor(col); h->SetMarkerStyle(sty); h->SetMarkerSize(size); h->SetLineWidth(2);
}
// scale contents and errors by f
void ScaleH(TH1* h, double f) { for (int b = 0; b < h->GetNcells(); ++b) { h->SetBinContent(b, h->GetBinContent(b) * f); h->SetBinError(b, h->GetBinError(b) * f); } }
TH1* Rebinned(TH1* h, int k) { TH1* r = h->Rebin(k, U("rb").c_str()); r->SetDirectory(nullptr); return r; }

// ratio pad: data (N - CB) / sim with data errors; band 1 +- (MC (+) 6.9 %); live time 1 .. 1.073 (hatched); (a) open markers
void RatioPad(TH1* sig, TH1* sim, TH1* simA, double xlo, double xhi, double ylo, double yhi, const char* xt, int lab, int tit,
              double xoff, double yoff, const std::vector<std::pair<double, double>>& masked = {}) {
    TH1* fr = new TH1D(U("frr").c_str(), "", 100, xlo, xhi);
    fr->SetDirectory(nullptr);
    fr->SetMinimum(ylo); fr->SetMaximum(yhi);
    fr->GetXaxis()->SetTitle(xt);
    fr->GetYaxis()->SetTitle("data/sim");
    fr->GetYaxis()->SetNdivisions(505);
    fr->GetXaxis()->SetNdivisions(505);
    Ax(fr, lab, tit, xoff, yoff);
    fr->Draw("AXIS");
    TH1* band = (TH1*)sim->Clone(U("band").c_str());
    band->SetDirectory(nullptr);
    for (int b = 1; b <= band->GetNbinsX(); ++b) {
        const double s = sim->GetBinContent(b);
        const double rel = s > 0 ? std::hypot(sim->GetBinError(b) / s, kLrelErr) : 0;
        band->SetBinContent(b, s > 0 ? 1. : -10.);
        band->SetBinError(b, rel);
    }
    band->SetFillColor(kGray); band->SetLineColor(kGray); band->SetMarkerSize(0); band->SetFillStyle(1001);
    band->Draw("E2 SAME");
    TBox* lt = new TBox(xlo, 1.0, xhi, 1.0 + kLiveUp);
    lt->SetFillStyle(3354); lt->SetFillColor(kOrange + 7); lt->SetLineColor(kOrange + 7);
    lt->Draw();
    TLine* one = new TLine(xlo, 1, xhi, 1); one->SetLineStyle(2); one->Draw();
    auto mk = [&](TH1* den, int col, int sty) {
        if (!den) return;
        TH1* r = (TH1*)sig->Clone(U("rat").c_str());
        r->SetDirectory(nullptr);
        for (int b = 1; b <= r->GetNbinsX(); ++b) {
            const double s = den->GetBinContent(b), x = r->GetXaxis()->GetBinCenter(b);
            bool m = false;
            for (auto& p : masked) if (x > p.first && x < p.second) m = true;
            if (s <= 0 || m || x < xlo || x > xhi) { r->SetBinContent(b, -100); r->SetBinError(b, 0); continue; }
            r->SetBinContent(b, sig->GetBinContent(b) / s);
            r->SetBinError(b, sig->GetBinError(b) / s);
        }
        Mark(r, col, sty, 0.9);
        r->Draw("E1 X0 SAME");
    };
    mk(simA, kColA, 24);
    mk(sim, kBlue + 1, 20);
    fr->Draw("AXIS SAME");
}

// legend entries shared by the figures
TH1* LegH(int col, int mst, int fill = -1, int fst = 1001, int lst = 1, int lw = 2) {
    TH1* h = new TH1D(U("lg").c_str(), "", 1, 0, 1);
    h->SetDirectory(nullptr);
    h->SetLineColor(col); h->SetMarkerColor(col); h->SetMarkerStyle(mst); h->SetLineStyle(lst); h->SetLineWidth(lw);
    if (fill >= 0) { h->SetFillColor(fill); h->SetFillStyle(fst); h->SetLineColor(fill); }
    return h;
}
void RatioLegend(TLegend* lg) {
    TH1* l1 = LegH(kBlue + 1, 20); l1->SetMarkerSize(1.4);
    TH1* l2 = LegH(kColA, 24); l2->SetMarkerSize(1.4);
    lg->AddEntry(l1, "(N - CB) / SMASH total", "p");
    lg->AddEntry(l2, "(N - CB) / (SMASH + ghosts)", "p");
    lg->AddEntry(LegH(kGray, 1, kGray), "band: MC #oplus L (6.9 %)", "f");
    lg->AddEntry(LegH(kOrange + 7, 1, kOrange + 7, 3354), "live time +7.3/-0 % (not applied)", "f");
}

// ---------------------------------------------------------------- B1 / B2: M4 spectra, S0 | pre | legend column
void FigM4(const char* lay, double xlo, double xhi, bool logy, const char* base) {
    const char* sels[2] = {"S0", "pre"};
    const char* selTxt[2] = {"(a) nominal: base selection + both pions in the vertex fit", "(b) base selection: PT3, v_{z} > -160 mm, corr. momenta"};
    TCanvas* c = new TCanvas(U("cm4").c_str(), "", 2100, 1150);
    double ymax = 0, ymin = 1e30;
    TH1 *hN[2], *hC[2], *hS[2], *hT[2], *hA[2];
    std::vector<TH1*> st[2];
    const std::vector<Stk> sd = StackDef(false);
    for (int i = 0; i < 2; ++i) {
        hN[i] = DataNb(sels[i], lay, 0); hC[i] = DataNb(sels[i], lay, 1); hS[i] = DataNb(sels[i], lay, 2);
        hT[i] = allio::Sim(gF, sels[i], "nom", kT12, lay); hA[i] = allio::Sim(gF, sels[i], "ghostA", kT12, lay);
        const double bw = hN[i]->GetXaxis()->GetBinWidth(1);
        for (TH1* h : {hN[i], hC[i], hS[i], hT[i], hA[i]}) ScaleH(h, 1. / bw);
        for (const Stk& s : sd) { TH1* h = SimSum(sels[i], "nom", s.comps, lay); ScaleH(h, 1. / bw); st[i].push_back(h); }
        for (int b = 1; b <= hN[i]->GetNbinsX(); ++b) {
            const double x = hN[i]->GetXaxis()->GetBinCenter(b);
            if (x < xlo || x > xhi) continue;
            ymax = std::max(ymax, std::max(hN[i]->GetBinContent(b), hA[i]->GetBinContent(b)));
            if (hS[i]->GetBinContent(b) > 0) ymin = std::min(ymin, hS[i]->GetBinContent(b));
        }
    }
    const double yhi = logy ? ymax * 3 : ymax * 1.08, ylo = logy ? std::max(ymin * 0.5, ymax * 1e-4) : 0.;
    const double W = 0.40;
    for (int i = 0; i < 2; ++i) {
        const double x0 = i * W;
        c->cd();
        TPad* pt = new TPad(U("pt").c_str(), "", x0, 0.31, x0 + W, 1.0);
        pt->SetLeftMargin(0.17); pt->SetRightMargin(0.025); pt->SetTopMargin(0.075); pt->SetBottomMargin(0.012);
        pt->SetLogy(logy);
        pt->Draw(); pt->cd();
        TH1* fr = new TH1D(U("fr").c_str(), "", 100, xlo, xhi);
        fr->SetDirectory(nullptr);
        fr->SetMinimum(ylo); fr->SetMaximum(yhi);
        fr->GetYaxis()->SetTitle("d#sigma/dM_{4} [nb/(GeV/c^{2})]");
        Ax(fr, 30, 34, 1.0, 1.75);
        fr->GetXaxis()->SetLabelSize(0);
        fr->Draw("AXIS");
        THStack* hs = new THStack(U("hs").c_str(), "");
        for (size_t k = 0; k < sd.size(); ++k) {
            st[i][k]->SetFillColor(sd[k].col); st[i][k]->SetLineColor(sd[k].col); st[i][k]->SetLineWidth(1);
            hs->Add(st[i][k], "HIST");
        }
        hs->Draw("SAME");
        hA[i]->SetLineColor(kColA); hA[i]->SetLineStyle(2); hA[i]->SetLineWidth(2);
        hA[i]->Draw("HIST SAME");
        for (int w : {kWB, kWP}) for (double x : {kWin[w].lo, kWin[w].hi}) {
            if (x < xlo || x > xhi) continue;
            TLine* l = new TLine(x, ylo, x, logy ? ylo * std::pow(yhi / ylo, 0.80) : 0.80 * yhi);
            l->SetLineStyle(3); l->SetLineColor(w == kWB ? kOrange + 7 : kGreen + 3); l->SetLineWidth(2); l->Draw();
        }
        Mark(hN[i], kBlack, 20, 0.85); Mark(hC[i], kRed + 1, 21, 0.8); Mark(hS[i], kBlue + 1, 20, 0.85);
        hC[i]->Draw("E1 X0 SAME"); hN[i]->Draw("E1 X0 SAME"); hS[i]->Draw("E1 X0 SAME");
        fr->Draw("AXIS SAME");
        TLatex tx; tx.SetNDC(); tx.SetTextFont(43); tx.SetTextSize(28);
        tx.DrawLatex(0.17, 0.945, selTxt[i]);
        TLatex tw; tw.SetTextFont(43); tw.SetTextSize(24);
        const double yl = logy ? ylo * std::pow(yhi / ylo, 0.81) : 0.81 * yhi;
        tw.SetTextAlign(21);   // _plain: words centred on the windows
        tw.SetTextColor(kOrange + 7); tw.DrawLatex(0.5 * (kWin[kWB].lo + kWin[kWB].hi), yl, "bump");
        tw.SetTextColor(kGreen + 3); tw.DrawLatex(0.5 * (kWin[kWP].lo + kWin[kWP].hi), logy ? ylo * std::pow(yhi / ylo, 0.87) : yl, "peak");
        c->cd();
        TPad* pb = new TPad(U("pb").c_str(), "", x0, 0.0, x0 + W, 0.31);
        pb->SetLeftMargin(0.17); pb->SetRightMargin(0.025); pb->SetTopMargin(0.02); pb->SetBottomMargin(0.36);
        pb->Draw(); pb->cd();
        RatioPad(hS[i], hT[i], hA[i], xlo, xhi, logy ? 0.3 : 0.5, logy ? 2.1 : 1.9, "M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]", 30, 34, 1.05, 1.75);
    }
    // legend column
    c->cd();
    TPad* pl = new TPad(U("pl").c_str(), "", 2 * W, 0.0, 1.0, 1.0);
    pl->Draw(); pl->cd();
    {
        TLatex tl; tl.SetNDC(); tl.SetTextFont(43); tl.SetTextSize(23);
        tl.DrawLatex(0.03, 0.965, "data: L = 6079.2 nb^{-1}, #times 64/63");
        tl.DrawLatex(0.03, 0.935, "SMASH: #sigma_{prod} #times BR #times #eta-decay weights");
        tl.DrawLatex(0.03, 0.905, "no renormalisation");
    }
    TLegend* lg = new TLegend(0.02, 0.34, 0.99, 0.89);
    lg->SetTextFont(43); lg->SetTextSize(24); lg->SetFillStyle(0); lg->SetBorderSize(0);
    lg->SetHeader("data (e^{+}e^{-}, all combinations)");
    lg->AddEntry(hN[0], "N", "p");
    lg->AddEntry(hC[0], "CB = 2 (N_{++}N_{--})^{1/2}", "p");
    lg->AddEntry(hS[0], "N - CB", "p");
    lg->AddEntry((TObject*)nullptr, "SMASH total (stacked):", "");
    for (int k = (int)sd.size() - 1; k >= 0; --k) lg->AddEntry(st[0][k], sd[k].label, "f");
    lg->AddEntry(hA[0], "SMASH + ghost-flagged pairs", "l");
    lg->Draw();
    TLegend* lr = new TLegend(0.02, 0.06, 0.99, 0.29);
    lr->SetTextFont(43); lr->SetTextSize(24); lr->SetFillStyle(0); lr->SetBorderSize(0);
    lr->SetHeader("ratio panels");
    RatioLegend(lr);
    lr->Draw();
    SaveNoRm(c, Form("%s/%s", kFig, base));
    delete c;
}

// ---------------------------------------------------------------- B3: map data/sim in M4 x M(pipi)
void FigMap(FILE* o) {
    TCanvas* c = new TCanvas(U("cmap").c_str(), "", 1800, 1500);
    c->Divide(2, 2, 0.002, 0.002);
    const char* sels[2] = {"S0", "pre"};
    const double zlo = 0.4, zhi = 1.6;
    {
        double st[3] = {0., 0.5, 1.};
        double r[3] = {0.15, 1., 0.80}, g[3] = {0.30, 1., 0.10}, b[3] = {0.75, 1., 0.10};
        TColor::CreateGradientColorTable(3, st, r, g, b, 99);
    }
    fprintf(o, "\n## B3 map (data/sim per cell, shown when the relative error < 0.20 and N - CB > 0): cells shown / cells with data\n");
    for (int i = 0; i < 2; ++i) {
        TH2* N = (TH2*)allio::Data(gF, sels[i], "all", "N", "B_map");
        TH2* CB = (TH2*)allio::Data(gF, sels[i], "all", "CB", "B_map");
        TH2* T = (TH2*)allio::Sim(gF, sels[i], "nom", kT12, "B_map");
        for (int z = 0; z < 2; ++z) {
            const int rb = z == 0 ? 2 : 1;   // full range: 20 MeV cells; zoom: 10 MeV
            TH2* n = (TH2*)N->Clone(U("mn").c_str()); TH2* cb = (TH2*)CB->Clone(U("mc").c_str()); TH2* t = (TH2*)T->Clone(U("mt").c_str());
            if (rb > 1) { n->Rebin2D(rb, rb); cb->Rebin2D(rb, rb); t->Rebin2D(rb, rb); }
            TH2* r = (TH2*)t->Clone(U("mr").c_str());
            r->Reset();
            int shown = 0, withData = 0;
            for (int bx = 1; bx <= r->GetNbinsX(); ++bx)
                for (int by = 1; by <= r->GetNbinsY(); ++by) {
                    const double nn = n->GetBinContent(bx, by), cc = cb->GetBinContent(bx, by), ss = t->GetBinContent(bx, by);
                    if (nn > 0) ++withData;
                    if (nn <= 0 || ss <= 0) { r->SetBinContent(bx, by, -10); continue; }
                    const double d = (nn - cc) * kRawToNb, de = std::hypot(n->GetBinError(bx, by), cb->GetBinError(bx, by)) * kRawToNb;
                    const double q = d / ss, qe = std::fabs(q) * std::hypot(de / std::max(1e-30, std::fabs(d)), t->GetBinError(bx, by) / ss);
                    if (!(qe < 0.20 * std::fabs(q)) || d <= 0) { r->SetBinContent(bx, by, -10); continue; }
                    ++shown;
                    r->SetBinContent(bx, by, std::min(zhi - 1e-3, std::max(zlo + 1e-3, q)));
                }
            fprintf(o, "%-4s %-5s cells shown %d of %d with data (cell %d MeV)\n", sels[i], z ? "zoom" : "full", shown, withData, 10 * rb);
            c->cd(1 + z + 2 * i);
            gPad->SetLeftMargin(0.13); gPad->SetRightMargin(0.16); gPad->SetTopMargin(0.07); gPad->SetBottomMargin(0.12);
            r->SetMinimum(zlo); r->SetMaximum(zhi);
            r->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]");
            r->GetYaxis()->SetTitle("M(#pi^{+}#pi^{-}) [GeV/c^{2}]");
            r->GetZaxis()->SetTitle("data (N - CB) / sim T1+T2");
            Ax(r, 26, 28, 1.3, 1.6);
            r->GetZaxis()->SetTitleOffset(1.5);
            const double x1 = z ? 0.40 : 0.30, x2 = z ? 0.70 : 1.10, y1 = z ? 0.27 : 0.25, y2 = z ? 0.60 : 1.00;
            r->GetXaxis()->SetRangeUser(x1, x2 - 1e-6);
            r->GetYaxis()->SetRangeUser(y1, y2 - 1e-6);
            r->SetContour(99);
            r->Draw("COLZ");
            TLine* kl = new TLine(std::max(x1, y1), std::max(x1, y1), std::min(x2, y2), std::min(x2, y2));
            kl->SetLineStyle(3); kl->SetLineColor(kGray + 2); kl->SetLineWidth(2); kl->Draw();
            for (int k = 0; k < 2; ++k) {
                TGraph* g = new TGraph();
                for (int j = 0; j <= 200; ++j) {
                    const double m = 0.27 + (kEdgeLimit - 0.27) * j / 200.;
                    const double e = Edge(m, k == 0 ? 0. : 0.060);
                    if (e > 0 && m >= y1) g->SetPoint(g->GetN(), e, m);
                }
                g->SetLineColor(kBlack); g->SetLineWidth(3); g->SetLineStyle(k == 0 ? 1 : 7);
                g->Draw("L");
            }
            auto box = [&](double a, double b, double cc, double d, int col, int sty) {
                TBox* bx = new TBox(std::max(a, x1), std::max(b, y1), std::min(cc, x2), std::min(d, y2));
                bx->SetFillStyle(0); bx->SetLineColor(col); bx->SetLineWidth(3); bx->SetLineStyle(sty); bx->Draw();
            };
            box(kWin[kWP].lo, y1, kWin[kWP].hi, std::min(y2, kWin[kWP].hi), kRed + 1, 1);   // peak band (M(pipi) < M4) [fix: red, as Figs. 1, C1]
            box(kWin[kWCRU].lo, y1, kWin[kWCRU].hi, std::min(y2, kWin[kWCRU].hi), kGreen + 2, 1); // CR-U [fix: green]
            box(kWin[kWCRK].lo, 0.43, kWin[kWCRK].hi, kWin[kWCRK].hi, kMagenta + 2, 1);         // CR-K (S5)
            box(kWin[kWCRL].lo, y1, kWin[kWCRL].hi, kWin[kWCRL].hi, kOrange + 7, 2);            // CR-L (test)
            TLatex tx; tx.SetTextFont(43); tx.SetTextSize(24);
            tx.SetTextColor(kRed + 1); tx.DrawLatex(kWin[kWP].lo + 0.001, std::min(y2, kWin[kWP].hi) + (z ? 0.004 : 0.012), "P");
            tx.SetTextColor(kGreen + 2); tx.DrawLatex(z ? 0.62 : 0.60, std::min(y2, 0.70) + (z ? -0.025 : 0.012), "CR-U");
            tx.SetTextColor(kMagenta + 2); tx.DrawLatex(z ? 0.47 : 0.405, z ? 0.50 : 0.5345, "CR-K");
            if (z) { tx.SetTextColor(kOrange + 7); tx.SetTextAngle(90); tx.DrawLatex(0.5285, 0.30, "CR-L"); tx.SetTextAngle(0); }
            TLatex tn; tn.SetNDC(); tn.SetTextFont(43); tn.SetTextSize(26);
            tn.DrawLatex(0.13, 0.945, Form("(%c) %s, %s, %d MeV cells", 'a' + z + 2 * i, sels[i], z ? "zoom" : "full range", 10 * rb));
            if (i == 0) {
                TLegend* lg = z ? new TLegend(0.14, 0.775, 0.41, 0.925) : new TLegend(0.15, 0.68, 0.55, 0.91);
                lg->SetTextFont(43); lg->SetTextSize(z ? 18 : 20); lg->SetFillStyle(1001); lg->SetFillColor(kWhite); lg->SetBorderSize(0);
                TGraph* g0 = new TGraph(); g0->SetLineWidth(3);
                TGraph* g1 = new TGraph(); g1->SetLineWidth(3); g1->SetLineStyle(7);
                lg->AddEntry(g0, "bump edge E(m, 0)", "l");
                lg->AddEntry(g1, "bump edge E(m, 60 MeV/c^{2})", "l");
                lg->AddEntry(kl, "M_{#pi#pi} = M_{4}", "l");
                if (!z) {
                    TGraph* gb = new TGraph(); gb->SetLineWidth(3); gb->SetLineColor(kRed + 1);
                    lg->AddEntry(gb, "peak band P", "l");
                    TGraph* gu = new TGraph(); gu->SetLineWidth(3); gu->SetLineColor(kGreen + 2);
                    lg->AddEntry(gu, "CR-U; CR-K (magenta); CR-L (orange)", "l");
                }
                lg->Draw();
            }
            delete n; delete cb; delete t;
        }
        delete N; delete CB; delete T;
    }
    SaveNoRm(c, Form("%s/fix_B3_plainrun", kScr));
    delete c;
}

// ---------------------------------------------------------------- B4: small multiples of the other variables
struct VarCfg { int var; int reb; double xlo, xhi; bool logy; const char* unit; double xconv; };
const char* RegText(int k) {
    static const char* t[kNReg] = {"R0: M_{4} < 0.455", "R1 = B: 0.455-0.5225", "R2: 0.5225-0.5336", "R3 = P: 0.5336-0.5542",
                                   "R4: 0.5542-0.565", "R5: 0.565-0.70", "R6: 0.70-0.90", "R7: 0.90-1.10"};
    return t[k];
}
void FigVar(const char* sel, const VarCfg& vc, FILE* o) {
    TCanvas* c = new TCanvas(U("cv").c_str(), "", 2000, 1250);
    const std::vector<Stk> sd = StackDef(true);
    fprintf(o, "%-4s %-6s", sel, BVarTag(vc.var));
    const double Htop = 0.935;   // legend strip above
    std::vector<TH1*> legH;
    for (int k = 0; k < kNReg; ++k) {
        const std::string lay = BVarName(vc.var, k);
        TH1* rawN = allio::Data(gF, sel, "all", "N", lay);
        std::vector<char> mask(rawN->GetNbinsX() + 2, 0);
        for (int b = 1; b <= rawN->GetNbinsX(); ++b) {
            if (rawN->GetBinContent(b) < 0) mask[b] = 1;
            if (vc.var == kBoa && k >= 2 && k <= 4) {
                const double a = rawN->GetXaxis()->GetBinLowEdge(b), bb = rawN->GetXaxis()->GetBinUpEdge(b);
                if (a >= 4. - 1e-9 && bb <= 30. + 1e-9) mask[b] = 1;
            }
        }
        TH1* S = DataNb(sel, lay, 2);
        for (int b = 1; b <= S->GetNbinsX(); ++b) if (mask[b]) { S->SetBinContent(b, 0); S->SetBinError(b, 0); }
        TH1* T = allio::Sim(gF, sel, "nom", kT12, lay);
        TH1* A = allio::Sim(gF, sel, "ghostA", kT12, lay);
        std::vector<TH1*> st;
        for (const Stk& s : sd) st.push_back(SimSum(sel, "nom", s.comps, lay));
        double ds = 0, ts = 0;
        for (int b = 1; b <= S->GetNbinsX(); ++b) {
            if (mask[b]) continue;
            ds += S->GetBinContent(b); ts += T->GetBinContent(b);
        }
        fprintf(o, " | %6.3f", ts > 0 ? ds / ts : 0.);
        TH1* Sr = Rebinned(S, vc.reb); TH1* Tr = Rebinned(T, vc.reb); TH1* Ar = Rebinned(A, vc.reb);
        std::vector<TH1*> str;
        for (TH1* h : st) str.push_back(Rebinned(h, vc.reb));
        std::vector<std::pair<double, double>> mranges;
        for (int b = 1; b <= Sr->GetNbinsX(); ++b) {
            bool m = false;
            for (int j = 0; j < vc.reb; ++j) if (mask[(b - 1) * vc.reb + j + 1]) m = true;
            if (m) {
                mranges.push_back({Sr->GetXaxis()->GetBinLowEdge(b), Sr->GetXaxis()->GetBinUpEdge(b)});
                Sr->SetBinContent(b, -1e30); Sr->SetBinError(b, 0);
            }
        }
        const double bw = Sr->GetXaxis()->GetBinWidth(1) * vc.xconv;
        for (TH1* h : {Sr, Tr, Ar}) ScaleH(h, 1. / bw);
        for (TH1* h : str) ScaleH(h, 1. / bw);
        double xhi = vc.xhi;
        if (vc.var == kBmpipi) xhi = std::min(vc.xhi, kWin[kWR0 + k].hi + 0.0001);
        if (vc.var == kBmee) xhi = std::min(vc.xhi, std::max(0.15, kWin[kWR0 + k].hi - 0.27));
        double ymax = 0, ymin = 1e30;
        for (int b = 1; b <= Sr->GetNbinsX(); ++b) {
            const double x = Sr->GetXaxis()->GetBinCenter(b);
            if (x < vc.xlo || x > xhi) continue;
            ymax = std::max(ymax, std::max(Sr->GetBinContent(b) + Sr->GetBinError(b), std::max(Tr->GetBinContent(b), Ar->GetBinContent(b))));
            if (Tr->GetBinContent(b) > 0) ymin = std::min(ymin, Tr->GetBinContent(b));
        }
        const double yhi = vc.logy ? ymax * 6 : ymax * 1.28, ylo = vc.logy ? std::max(ymin * 0.5, ymax * 1e-4) : 0.;
        const int col = k % 4, row = k / 4;
        const double x0 = col * 0.25, y0 = row == 0 ? 0.5 * Htop : 0.0, hh = 0.5 * Htop;
        c->cd();
        TPad* pt = new TPad(U("vt").c_str(), "", x0, y0 + 0.36 * hh, x0 + 0.25, y0 + hh);
        pt->SetLeftMargin(0.22); pt->SetRightMargin(0.03); pt->SetTopMargin(0.11); pt->SetBottomMargin(0.015);
        pt->SetLogy(vc.logy);
        pt->Draw(); pt->cd();
        TH1* fr = new TH1D(U("fv").c_str(), "", 100, vc.xlo, xhi);
        fr->SetDirectory(nullptr);
        fr->SetMinimum(ylo); fr->SetMaximum(yhi);
        fr->GetYaxis()->SetTitle(vc.var == kBy4 ? "d#sigma/dy_{cm} [nb]" : Form("d#sigma/dx [nb/%s]", vc.unit));
        Ax(fr, 22, 24, 1.0, 2.5);
        fr->GetXaxis()->SetLabelSize(0);
        fr->GetYaxis()->SetNdivisions(505);
        fr->Draw("AXIS");
        THStack* hs = new THStack(U("hsv").c_str(), "");
        for (size_t j = 0; j < sd.size(); ++j) {
            str[j]->SetFillColor(sd[j].col); str[j]->SetLineColor(sd[j].col); str[j]->SetLineWidth(1);
            hs->Add(str[j], "HIST");
        }
        hs->Draw("SAME");
        Ar->SetLineColor(kColA); Ar->SetLineStyle(2); Ar->SetLineWidth(2); Ar->Draw("HIST SAME");
        for (auto& p : mranges) {
            if (p.second <= vc.xlo || p.first >= xhi) continue;
            TBox* bx = new TBox(std::max(p.first, vc.xlo), ylo, std::min(p.second, xhi), yhi);
            bx->SetFillStyle(3345); bx->SetFillColor(kGray + 2); bx->SetLineColor(kGray + 2); bx->Draw();
        }
        Mark(Sr, kBlue + 1, 20, 0.7);
        Sr->Draw("E1 X0 SAME");
        fr->Draw("AXIS SAME");
        TLatex tx; tx.SetNDC(); tx.SetTextFont(43); tx.SetTextSize(22);
        tx.DrawLatex(0.22, 0.92, Form("%s", RegText(k)));
        tx.SetTextAlign(31); tx.DrawLatex(0.96, 0.80, Form("data/sim = %.3f", ts > 0 ? ds / ts : 0.)); tx.SetTextAlign(11);
        if (k == 0) { legH.push_back(Sr); for (TH1* h : str) legH.push_back(h); legH.push_back(Ar); }
        c->cd();
        TPad* pb = new TPad(U("vb").c_str(), "", x0, y0, x0 + 0.25, y0 + 0.36 * hh);
        pb->SetLeftMargin(0.22); pb->SetRightMargin(0.03); pb->SetTopMargin(0.02); pb->SetBottomMargin(0.46);
        pb->Draw(); pb->cd();
        RatioPad(Sr, Tr, Ar, vc.xlo, xhi, 0.2, 2.2, BVarTitle(vc.var), 22, 24, 0.95, 2.5, mranges);
        delete rawN; delete S; delete T; delete A;
        for (TH1* h : st) delete h;
    }
    fprintf(o, "\n");
    // legend strip
    c->cd();
    TPad* pl = new TPad(U("vl").c_str(), "", 0.0, Htop, 1.0, 1.0);
    pl->Draw(); pl->cd();
    TLegend* lg = new TLegend(0.005, 0.02, 0.995, 0.98);
    lg->SetNColumns(4); lg->SetTextFont(43); lg->SetTextSize(21); lg->SetFillStyle(0); lg->SetBorderSize(0);
    legH[0]->SetMarkerSize(1.2);
    lg->AddEntry(legH[0], Form("data N - CB (%s)", sel), "p");
    for (int j = (int)sd.size() - 1; j >= 0; --j) lg->AddEntry(legH[1 + j], sd[j].label, "f");
    lg->AddEntry(legH.back(), "sim (a): + ghost-flagged true pairs", "l");
    lg->AddEntry(LegH(kGray + 2, 1, kGray + 2, 3345), "data not shown (blinding)", "f");
    lg->AddEntry(LegH(kGray, 1, kGray), "ratio band: MC #oplus L 6.9 %", "f");
    lg->Draw();
    SaveNoRm(c, Form("%s/wp1_B4_%s%s", kScr, BVarTag(vc.var), std::string(sel) == "S0" ? "" : "_pre"));
    delete c;
}

// ---------------------------------------------------------------- K0S peak in M(pipi), all regions R0-R7
struct K0Res { double y = 0, ey = 0, mu = 0, emu = 0, sg = 0, esg = 0, chi2 = 0; int ndf = 0; bool ok = false; };
K0Res FitK0(TH1* h, const char* name) {
    K0Res r;
    TF1* f = new TF1(name, "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]+[4]*(x-0.4961)+[5]*(x-0.4961)^2", 0.455, 0.545);
    const int b1 = h->FindBin(0.462), b2 = h->FindBin(0.530), bp = h->FindBin(0.4961);
    const double bl = 0.5 * (h->GetBinContent(b1) + h->GetBinContent(b2));
    f->SetParameters(std::max(1e-9, h->GetBinContent(bp) - bl), 0.4961, 0.0058, bl, (h->GetBinContent(b2) - h->GetBinContent(b1)) / 0.068, 0.);
    f->SetParLimits(1, 0.488, 0.504);
    f->SetParLimits(2, 0.0025, 0.012);
    TFitResultPtr fp = h->Fit(f, "QRNS0");
    if (fp.Get() && (int)fp == 0) {
        const double bw = h->GetXaxis()->GetBinWidth(1), k = std::sqrt(2 * M_PI) / bw;
        r.y = f->GetParameter(0) * f->GetParameter(2) * k;
        const double cov02 = fp->CovMatrix(0, 2);
        r.ey = k * std::sqrt(std::pow(f->GetParameter(2) * f->GetParError(0), 2) + std::pow(f->GetParameter(0) * f->GetParError(2), 2) +
                             2 * f->GetParameter(0) * f->GetParameter(2) * cov02);
        r.mu = f->GetParameter(1); r.emu = f->GetParError(1); r.sg = f->GetParameter(2); r.esg = f->GetParError(2);
        r.chi2 = f->GetChisquare(); r.ndf = f->GetNDF(); r.ok = true;
    }
    return r;
}
}  // namespace wp1

void fix_wp1_plain() {
    using namespace wp1fixplain;
    bla::Style();
    gStyle->SetHatchesLineWidth(1);
    gStyle->SetHatchesSpacing(1.2);
    gF = TFile::Open(kHistFile);
    if (!gF || gF->IsZombie()) { printf("fix_wp1: cannot open %s\n", kHistFile); return; }
    // _plain: the call order of fix_wp1.C is kept (FigM4 alone crashes in THStack painting when FigMap has not run first);
    // the map and its statistics go to the scratch directory kScr, the two M4 figures of the short report to figs_all.
    FILE* o = fopen(Form("%s/fix_wp1_map_plainrun.txt", kScr), "w");
    FigMap(o);
    FigM4("B_M4", 0.30, 1.10, true, "fix_B1_plain");
    FigM4("fine_M4", 0.40, 0.70, false, "wp1_B2_plain");
    fclose(o);
}
