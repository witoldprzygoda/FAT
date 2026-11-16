// ============================================================================
// opt_scans.C  --  cut scans on the data for the report: single-cut scans,
//                  N-1 scans, cut flows and the vertex cuts
// ============================================================================
// HADES p+p, T = 4.53 GeV (feb22), channel pi+ pi- e+ e-. EXPERIMENTAL DATA ONLY.
// Re-presentation of the cut-scan study (plots/scan_peak.C, scan_bump.C,
// scan_extra.C, scan_ribbons.C) on the base selection of the report,
// pre = cor_ok && isBest >= 0 && vz > -160 mm (the first study started without
// the vz cut; its "none" = "pre without vz > -160" here).
//
// What it does
//   1. eta peak, single-cut scans from pre: every cut alone, P_main fit of the
//      passing rows (S, sigma_S, F = S/sigma_S, S/B in the data window) and of
//      the failing rows; eps_S = S_pass/(S_pass + S_fail), eps_B = kept share of
//      (pol2 + CB) in the window; the cut improves the statistical precision only
//      if eps_S > sqrt(eps_B); the two halves of the event groups (fold 0, 1);
//      vz threshold scanned from "pre without vz > -160";
//   2. eta peak, N-1 scans around "peak nominal + vertex" (each cut varied with
//      the others at their nominal values);
//   3. bump, single-cut scans from pre and N-1 scans around "bump nominal +
//      vertex": edge-step visibility E (scan_bump.C procedure), N_b = rows in the
//      bump window 0.4550-0.5225, gain_E = (E/E_base)/sqrt(N_b/N_b,base), and the
//      ECAL-tagged eta T with gain_T = (T/T_base)/sqrt(N_b/N_b,base) (ECAL:
//      verification only);
//   4. cut flows of the peak and bump selections (and of the frozen PEAK) with
//      bands: at each step the cut of this step moved one step down and up;
//      not-applied cuts shown after the last step;
//   5. vertex: signal per vz slice (peak S and tagged T, CB share) and the vertex
//      cuts (vz threshold, pions in the vertex fit, dz_pi).
//   Sculpting flags (hatched values) are read from tables_exp/opt_sculpt.txt
//   (written by opt_sculpt.C; run that first).
//
// Fixed conventions: corrected momenta, cor_ok == 1, isBest >= 0 (isBest == 1
// only as a comparison elsewhere), FAT vertex z > -500 mm and PT3; CB = 2
// sqrt(N++ N--) per 10 MeV, k = 1; P_main with the frozen Gaussian 0.5439 GeV/c2,
// sigma 5.4 MeV + pol2 in 0.5225-0.6100; F compares selections, it is NOT a
// significance; E is NOT a yield.
//
// Inputs : /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/scan_v2/scan_compact.root (scan::kCompact)
//          /home/przygoda/HADES/FAT/FAT/docs/latex/tables_exp/opt_sculpt.txt (flags)
// Outputs: /home/przygoda/HADES/FAT/FAT/docs/latex/tables_exp/opt_scans.txt, opt_cutflow.tex
//          /home/przygoda/HADES/FAT/FAT/docs/latex/figs_exp/opt_peak_single.pdf, opt_peak_single_eff.pdf, opt_peak_nm1.pdf,
//          opt_bump_single.pdf, opt_bump_nm1.pdf, opt_cutflow.pdf, opt_vertex.pdf (+ .png)
//
// Verifier changes (2026-11-14), layout only (numbers unchanged): scan figures in
// 3 columns instead of 5, crowded category axes labelled at every second value,
// cut-flow canvas 1300 x 1000 (was 1800 x 1250) with 100 px left margins, vertex
// canvas 1150 x 800 (was 1500 x 1000) with a one-column legend in panel (a): the
// labels were 3-5 pt at the text width.
//
// Run (from the repository root /home/przygoda/HADES/FAT/FAT), about 5 min:
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/aclic_opt", kTRUE)' 'plots/exp_report/opt_scans.C+'
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "TBox.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TSystem.h"

#include "opt_common.h"

namespace optscan {

using namespace opt;

// ============================================================================
// variables
// ============================================================================
enum VarId { kOaLo = 0, kOaHi, kOaLab, kMpLo, kMpHi, kMm4, kEta, kDz, kVtx, kVz, kNVar };
struct VarDef {
    int id = 0;
    std::string fam;              // family name as in opt_sculpt.txt
    std::string axis;             // axis title
    std::vector<double> val;      // NaN = no cut
    std::vector<std::string> key; // value key as in opt_sculpt.txt
    std::vector<std::string> lab; // label on the axis
    std::vector<int> flag;        // 1 = flagged by the sculpting test
};
const double kNaN = std::numeric_limits<double>::quiet_NaN();
inline void SetVar(Cuts& c, int id, double v) {
    const bool none = !(v == v);
    switch (id) {
        case kOaLo:  c.oaLo = none ? kNoCut : (float)v; break;
        case kOaHi:  c.oaHi = none ? kNoCut : (float)v; break;
        case kOaLab: c.oaLab = none ? kNoCut : (float)v; break;
        case kMpLo:  c.mpLo = none ? kNoCut : (float)v; break;
        case kMpHi:  c.mpHi = none ? kNoCut : (float)v; break;
        case kMm4:   c.mm4Hi = none ? kNoCut : (float)v; break;
        case kEta:   c.etaLo = none ? kNoCut : (float)v; break;
        case kDz:    c.dzPi = none ? kNoCut : (float)v; break;
        case kVtx:   c.vtx = none ? 0 : (int)v; break;
        case kVz:    c.vzMin = none ? kNoCut : (float)v; break;
    }
}
std::vector<VarDef> gV;
void BuildVars() {
    gV.assign(kNVar, VarDef());
    auto mk = [&](int id, const char* fam, const char* axis, const std::vector<double>& v, const char* fmt) {
        VarDef& d = gV[id]; d.id = id; d.fam = fam; d.axis = axis; d.val = v;
        for (double x : v) {
            if (!(x == x)) { d.key.push_back("none"); d.lab.push_back("none"); continue; }
            d.key.push_back(Form(fmt, x)); d.lab.push_back(Form(fmt, x));
        }
        d.flag.assign(v.size(), 0);
    };
    mk(kOaLo, "OA >", "lower cut on OA(e^{+}e^{-}) [deg]", {kNaN, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0}, "%.1f");
    mk(kOaHi, "OA <", "upper cut on OA(e^{+}e^{-}) [deg]", {4, 6, 8, 10, 15, 20, 30, kNaN}, "%.0f");
    { std::vector<double> v; for (int i = 25; i <= 70; i += 5) v.push_back(i); v.push_back(kNaN);
      mk(kOaLab, "oa_lab <", "upper cut on oa_{lab}((#pi#pi),(ee)) [deg]", v, "%.0f"); }
    mk(kMpLo, "M(pipi) >", "lower cut on M(#pi^{+}#pi^{-}) [GeV/c^{2}]", {kNaN, 0.28, 0.30, 0.32, 0.34, 0.36, 0.38}, "%.2f");
    mk(kMpHi, "M(pipi) <", "upper cut on M(#pi^{+}#pi^{-}) [GeV/c^{2}]", {0.38, 0.39, 0.40, 0.41, 0.42, 0.43, 0.44, 0.46, 0.48, 0.50, kNaN}, "%.2f");
    mk(kMm4, "MM4 <", "upper cut on MM_{4} [GeV/c^{2}]", {2.60, 2.65, 2.70, 2.75, 2.80, 2.85, 2.90, 2.95, kNaN}, "%.2f");
    mk(kEta, "eta_rest >", "lower cut on the #eta-rest angle [deg]", {kNaN, 100, 120, 140, 150, 160}, "%.0f");
    mk(kDz, "dz_pi <", "upper cut on #Deltaz_{#pi} [mm]", {10, 15, 20, 30, 40, 60, 100, kNaN}, "%.0f");
    mk(kVtx, "vertex fit", "tracks in the event-vertex fit", {kNaN, 1, 2, 4}, "%.0f");
    gV[kVtx].key = {"none", ">=1 pi", "both pi", "all 4"};
    gV[kVtx].lab = {"none", "#geq1 #pi", "both #pi", "all 4"};
    mk(kVz, "vz >", "lower cut on the vertex z [mm]", {kNaN, -200, -180, -170, -160, -150, -140, -130, -120}, "%.0f");
    // flags from the sculpting test
    std::ifstream in(std::string(kTabOut) + "/opt_sculpt.txt");
    int nflag = 0;
    if (!in) printf("WARNING: %s/opt_sculpt.txt not found: no sculpting flags\n", kTabOut);
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#' || line.size() < 40 || line.find('|') == std::string::npos) continue;
        std::string fam = line.substr(0, 11), key = line.substr(12, 22);
        auto trim = [](std::string s) { size_t a = s.find_first_not_of(' '), b = s.find_last_not_of(' ');
                                        return a == std::string::npos ? std::string() : s.substr(a, b - a + 1); };
        fam = trim(fam); key = trim(key);
        const size_t p = line.rfind('|');
        int G = 0, M = 0;
        if (sscanf(line.c_str() + p + 1, "%d %d", &G, &M) != 2) continue;
        for (auto& d : gV)
            if (d.fam == fam)
                for (size_t i = 0; i < d.key.size(); ++i)
                    if (d.key[i] == key && (G || M)) { d.flag[i] = 1; nflag++; }
    }
    printf("sculpting flags read: %d\n", nflag);
}
int NoneIdx(int v) { for (size_t i = 0; i < gV[v].val.size(); ++i) if (!(gV[v].val[i] == gV[v].val[i])) return (int)i; return -1; }

// ============================================================================
// peak scan points
// ============================================================================
struct PPt { PeakFit pass, fail, f0, f1; double epsS = 1, eEpsS = 0, epsB = 1, z = 0; bool hasFail = false; };
void PassFail(PPt& q) {
    if (!q.hasFail) return;
    Ratio(q.pass.S, q.pass.eS, q.fail.S, q.fail.eS, q.epsS, q.eEpsS);
    const double bp = q.pass.B + q.pass.CB, bf = q.fail.B + q.fail.CB;
    q.epsB = (bp + bf) > 0 ? bp / (bp + bf) : 1;
    q.z = q.eEpsS > 0 ? (q.epsS - std::sqrt(std::max(q.epsB, 0.))) / q.eEpsS : 0;
}
// scan of the variables 'vars' around 'base' (each variable alone at all values, the rest as in base)
struct PScan { Cuts base; std::vector<int> vars; std::vector<std::vector<PPt>> pt; };
PScan PeakScan(const scan::Data& d, const Cuts& base, const std::vector<int>& vars, const Cuts* vzBase = nullptr) {
    PScan s; s.base = base; s.vars = vars;
    std::vector<Sel> sels;
    for (int v : vars) {
        const Cuts b0 = (v == kVz && vzBase) ? *vzBase : base;
        Cuts open = b0; SetVar(open, v, kNaN);
        for (double x : gV[v].val) {
            Cuts w = b0; SetVar(w, v, x);
            sels.push_back(MkSel(w));
            sels.push_back([=](const Row& r) { return Pass(open, r) && !Pass(w, r); });
            sels.push_back(MkSel(WithFold(w, 0)));
            sels.push_back(MkSel(WithFold(w, 1)));
        }
    }
    const std::vector<Spectra> sp = scan::FillMany(d, sels);
    size_t k = 0;
    for (int v : vars) {
        std::vector<PPt> row;
        for (size_t i = 0; i < gV[v].val.size(); ++i) {
            PPt q;
            q.pass = scan::FitPeakMain(sp[k]);
            q.hasFail = sp[k + 1].nsel[0] > 0;
            q.fail = q.hasFail ? scan::FitPeakMain(sp[k + 1]) : PeakFit();
            if (!q.hasFail) { q.fail.S = 0; q.fail.eS = 0; }
            q.hasFail = true;
            q.f0 = scan::FitPeakMain(sp[k + 2]); q.f1 = scan::FitPeakMain(sp[k + 3]);
            PassFail(q);
            row.push_back(q);
            k += 4;
        }
        s.pt.push_back(row);
    }
    return s;
}

// ============================================================================
// bump scan points
// ============================================================================
struct BPt { BumpRes r; double effE = 0, eEffE = 0, effN = 0, gE = 0, egE = 0, effT = 0, eEffT = 0, gT = 0, egT = 0; };
struct BScan { Cuts base; std::vector<int> vars; std::vector<std::vector<BPt>> pt; };
BScan BumpScan(const scan::Data& d, const scan::Data& dT, const Cuts& base, const std::vector<int>& vars, const Cuts* vzBase = nullptr) {
    BScan s; s.base = base; s.vars = vars;
    std::vector<Cuts> cs, refs;
    for (int v : vars) {
        const Cuts b0 = (v == kVz && vzBase) ? *vzBase : base;
        Cuts open = b0; SetVar(open, v, kNaN);
        refs.push_back(open);
        for (double x : gV[v].val) { Cuts w = b0; SetVar(w, v, x); cs.push_back(w); }
    }
    const std::vector<BumpRes> br = EvalBump(d, cs, false, true, &dT);
    const std::vector<BumpRes> b0 = EvalBump(d, refs, false, true, &dT);
    size_t k = 0;
    for (size_t iv = 0; iv < vars.size(); ++iv) {
        const int v = vars[iv];
        const BumpRes& R0 = b0[iv];
        std::vector<BPt> row;
        for (size_t i = 0; i < gV[v].val.size(); ++i) {
            BPt q; q.r = br[k++];
            SubRatio(q.r.E, q.r.eE, R0.E, R0.eE, q.effE, q.eEffE);
            q.effN = R0.win.n[0] > 0 ? q.r.win.n[0] / R0.win.n[0] : 0;
            SubRatio(q.r.T.S, q.r.T.eS, R0.T.S, R0.T.eS, q.effT, q.eEffT);
            if (q.effN > 0) {
                const double sq = std::sqrt(q.effN);
                q.gE = q.effE / sq; q.egE = q.eEffE / sq; q.gT = q.effT / sq; q.egT = q.eEffT / sq;
            }
            row.push_back(q);
        }
        s.pt.push_back(row);
    }
    return s;
}
// bump flag of a scan point: the summed edge fit is bad (p < 1e-3), e.g. a band-pass in M4
inline bool BumpFlag(const BPt& q) { return q.r.ndf > 0 && TMath::Prob(q.r.chi2, q.r.ndf) < 1e-3; }

// ============================================================================
// drawing helpers
// ============================================================================
void AxPx(TH1* h, int lab = 22, int tit = 24, double xoff = 1.0, double yoff = 1.6) {   // [review fix] larger text
    for (TAxis* ax : {h->GetXaxis(), h->GetYaxis()}) { ax->SetLabelFont(43); ax->SetTitleFont(43); ax->SetLabelSize(lab); ax->SetTitleSize(tit); }
    h->GetXaxis()->SetTitleOffset(xoff); h->GetYaxis()->SetTitleOffset(yoff);
}
TLatex* TxtPx(double x, double y, const char* t, int px = 18, int col = kBlack, int align = 11) {
    TLatex* l = new TLatex(x, y, t); l->SetNDC(); l->SetTextFont(43); l->SetTextSize(px);
    l->SetTextColor(col); l->SetTextAlign(align); l->Draw(); return l;
}
// categorical frame of one variable
TH1D* VarFrame(const char* name, int v, double y1, double y2, const char* yt) {
    const int n = (int)gV[v].val.size();
    TH1D* h = new TH1D(name, "", n, 0, n); h->SetDirectory(nullptr);
    // crowded axes (9 or more values): label every second value (both ends kept; verifier change 2026-11-14)
    for (int i = 0; i < n; ++i) h->GetXaxis()->SetBinLabel(i + 1, (n < 9 || i % 2 == 0 || i == n - 1) ? gV[v].lab[i].c_str() : "");
    h->SetMinimum(y1); h->SetMaximum(y2);
    h->GetXaxis()->SetTitle(gV[v].axis.c_str()); h->GetYaxis()->SetTitle(yt);
    AxPx(h, n > 5 ? 20 : 22, 24, 1.10, 1.45);   // [review fix] larger text (was 16/18, 20)
    h->GetXaxis()->SetLabelOffset(0.012);
    h->LabelsOption("h", "X");
    h->Draw("AXIS");
    return h;
}
void Hatch(int v, double y1, double y2) {
    for (size_t i = 0; i < gV[v].flag.size(); ++i)
        if (gV[v].flag[i]) { TBox* b = new TBox(i, y1, i + 1, y2); b->SetFillStyle(3354); b->SetFillColor(kGray + 1); b->SetLineWidth(0); b->Draw(); }
}
// panel grid inside a canvas: n panels in nc columns
TPad* Cell(TCanvas* c, int k, int nc, int nr, double top = 0.93, double bot = 0.0) {
    c->cd();
    const int ix = k % nc, iy = k / nc;
    const double h = (top - bot) / nr;
    TPad* p = new TPad(Form("%s_c%d", c->GetName(), k), "", (double)ix / nc, top - (iy + 1) * h, (double)(ix + 1) / nc, top - iy * h);
    p->SetLeftMargin(0.20); p->SetRightMargin(0.03); p->SetTopMargin(0.06); p->SetBottomMargin(0.21);
    p->Draw(); p->cd();
    return p;
}

// Scan figures: 3 columns (verifier change 2026-11-14: was 5 columns, which gave
// 3-4 pt labels at the text width); the canvas has a header of kHeadPx pixels
// with up to three text lines.
const int kScanNc = 3, kHeadPx = 115;
TCanvas* ScanCanvas(const char* fname, int nv, int& nc, int& nr, double& top) {
    nc = kScanNc; nr = (nv + nc - 1) / nc;
    const int H = 330 * nr + kHeadPx;
    TCanvas* c = new TCanvas(Form("c_%s", fname), "", 400 * nc, H);
    top = 1.0 - (double)kHeadPx / H;
    return c;
}
void ScanHead(TCanvas* c, const char* l1, const char* l2, const char* l3) {
    const double H = c->GetWh();
    c->cd();
    TxtPx(0.01, 1.0 - 8.0 / H, l1, 23, kBlack, 13);
    if (l2 && *l2) TxtPx(0.01, 1.0 - 42.0 / H, l2, 17, kGray + 3, 13);
    if (l3 && *l3) TxtPx(0.01, 1.0 - 68.0 / H, l3, 17, kGray + 3, 13);
}

// F (and folds) of a peak scan: one panel per variable
void DrawPeakF(const char* fname, const char* head, const PScan& s, double Fref, const char* refLab) {
    const int nv = (int)s.vars.size();
    int nc, nr; double top;
    TCanvas* c = ScanCanvas(fname, nv, nc, nr, top);
    ScanHead(c, head, Form("blue: all rows; grey / orange open: the two halves of the event groups (F #times #sqrt{2}); dashed: %s;", refLab),
             "green: values with F #geq F_{best} - 1; hatched: flagged by the sculpting test");
    for (int iv = 0; iv < nv; ++iv) {
        const int v = s.vars[iv];
        Cell(c, iv, nc, nr, top);
        double fmax = 0, fmin = 1e9;
        for (const PPt& q : s.pt[iv]) { fmax = std::max(fmax, q.pass.F); fmin = std::min(fmin, q.pass.F);
            fmax = std::max(fmax, std::sqrt(2.) * std::max(q.f0.F, q.f1.F)); fmin = std::min(fmin, std::sqrt(2.) * std::min(q.f0.F, q.f1.F)); }
        const double y1 = std::max(0., fmin - 1.0), y2 = fmax + 1.2;
        TH1D* fr = VarFrame(Form("f_%s_%d", fname, iv), v, y1, y2, "F = S/#sigma_{S}");
        // plateau
        double fb = -1; for (size_t i = 0; i < s.pt[iv].size(); ++i) if (!gV[v].flag[i]) fb = std::max(fb, s.pt[iv][i].pass.F);
        for (size_t i = 0; i < s.pt[iv].size(); ++i)
            if (!gV[v].flag[i] && s.pt[iv][i].pass.F >= fb - 1.0) { TBox* b = new TBox(i, fb - 1.0, i + 1, fb); b->SetFillColor(TColor::GetColor("#cdf2cd")); b->SetLineWidth(0); b->Draw(); }
        Hatch(v, y1, y2);
        HLine(0, gV[v].val.size(), Fref, kGray + 2, 2, 1);
        TGraphErrors* g = GE(kBlue + 1, 20, 1.0); TGraphErrors* g0 = GE(kGray + 2, 24, 0.8); TGraphErrors* g1 = GE(kOrange + 7, 25, 0.8);
        for (size_t i = 0; i < s.pt[iv].size(); ++i) {
            AddPt(g, i + 0.5, s.pt[iv][i].pass.F, 0);
            AddPt(g0, i + 0.35, std::sqrt(2.) * s.pt[iv][i].f0.F, 0);
            AddPt(g1, i + 0.65, std::sqrt(2.) * s.pt[iv][i].f1.F, 0);
        }
        g->SetLineWidth(2); g0->Draw("P SAME"); g1->Draw("P SAME"); g->Draw("LP SAME");
        fr->Draw("AXIS SAME");
    }
    Save(c, fname);
}
// eps_S against sqrt(eps_B) of a peak scan
void DrawPeakEff(const char* fname, const char* head, const PScan& s) {
    const int nv = (int)s.vars.size();
    int nc, nr; double top;
    TCanvas* c = ScanCanvas(fname, nv, nc, nr, top);
    ScanHead(c, head, "blue: #varepsilon_{S} = S_{pass}/(S_{pass} + S_{fail}) (stat. error); black line: #sqrt{#varepsilon_{B}} (kept share of pol2 + CB in the window);",
             "the cut helps only where blue lies above black");
    for (int iv = 0; iv < nv; ++iv) {
        const int v = s.vars[iv];
        Cell(c, iv, nc, nr, top);
        TH1D* fr = VarFrame(Form("e_%s_%d", fname, iv), v, 0, 1.25, "fraction kept");
        Hatch(v, 0, 1.25);
        HLine(0, gV[v].val.size(), 1, kGray + 1, 2, 1);
        TGraphErrors* g = GE(kBlue + 1, 20, 1.0); TGraphErrors* gb = GE(kBlack, 1, 0.1);
        for (size_t i = 0; i < s.pt[iv].size(); ++i) {
            AddPt(g, i + 0.5, s.pt[iv][i].epsS, s.pt[iv][i].eEpsS);
            AddPt(gb, i + 0.5, std::sqrt(std::max(s.pt[iv][i].epsB, 0.)), 0);
        }
        gb->SetLineWidth(2); gb->Draw("L SAME"); g->Draw("PZ SAME");
        fr->Draw("AXIS SAME");
    }
    Save(c, fname);
}
// gain_E and gain_T of a bump scan
void DrawBump(const char* fname, const char* head, const BScan& s) {
    const int nv = (int)s.vars.size();
    int nc, nr; double top;
    TCanvas* c = ScanCanvas(fname, nv, nc, nr, top);
    ScanHead(c, head, "gain = (kept share of the bump measure) / #sqrt{kept share of the rows N_{b} in the bump window}; > 1: the cut improves the visibility;",
             "violet: edge step E (stat. error); blue: ECAL-tagged #eta T (verification); hatched: edge fit fails (p < 10^{-3})");
    for (int iv = 0; iv < nv; ++iv) {
        const int v = s.vars[iv];
        Cell(c, iv, nc, nr, top);
        TH1D* fr = VarFrame(Form("b_%s_%d", fname, iv), v, 0.55, 1.45, "gain");
        for (size_t i = 0; i < s.pt[iv].size(); ++i)
            if (BumpFlag(s.pt[iv][i])) { TBox* b = new TBox(i, 0.55, i + 1, 1.45); b->SetFillStyle(3354); b->SetFillColor(kGray + 1); b->SetLineWidth(0); b->Draw(); }
        HLine(0, gV[v].val.size(), 1, kGray + 2, 2, 1);
        TGraphErrors* gE = GE(kViolet + 1, 21, 1.0); TGraphErrors* gT = GE(kBlue + 1, 20, 1.1);   // [review fix] E not red (red = CB)
        for (size_t i = 0; i < s.pt[iv].size(); ++i) {
            AddPt(gE, i + 0.38, s.pt[iv][i].gE, s.pt[iv][i].egE);
            AddPt(gT, i + 0.62, s.pt[iv][i].gT, s.pt[iv][i].egT);
        }
        gE->Draw("PZ SAME"); gT->Draw("PZ SAME");
        fr->Draw("AXIS SAME");
    }
    Save(c, fname);
}

// ============================================================================
// cut flow
// ============================================================================
struct Step { std::string lab, sub; Cuts c, lo, hi; bool applied = true; bool band = false; };
struct StepRes { PeakRes p, pLo, pHi; BumpRes b, bLo, bHi; };

} // namespace optscan

// ============================================================================
void opt_scans() {
    using namespace opt;
    using namespace optscan;
    TStopwatch sw;
    gSystem->mkdir(kFigOut, kTRUE); gSystem->mkdir(kTabOut, kTRUE);
    Style();
    BuildVars();
    scan::Data d;
    if (!scan::Load(d, scan::kCompact, scan::kMlo, scan::kMhi)) return;
    scan::Data dT;
    Filter(d, dT, [](const Row& r) { return HasGam(r) && r.m_eeg >= 0.05f && r.m_eeg < 0.25f; });
    const Cuts pre = Pre(), pre0 = Pre0();

    // ========================================================================
    // 1-2. peak scans
    // ========================================================================
    const std::vector<int> vSingle = {kOaLo, kOaHi, kOaLab, kMpLo, kMpHi, kMm4, kEta, kDz, kVtx, kVz};
    const PScan ps = PeakScan(d, pre, vSingle, &pre0);
    printf("peak single at %.1f s\n", sw.RealTime()); sw.Continue();
    const Cuts pnv = PeakNomV();
    const std::vector<int> vNm1 = {kOaLab, kMpLo, kVtx, kOaLo, kOaHi, kDz, kMm4, kMpHi, kVz};
    const PScan pn = PeakScan(d, pnv, vNm1);
    printf("peak N-1 at %.1f s\n", sw.RealTime()); sw.Continue();
    const PeakFit fPre = scan::FitPeakMain(scan::Fill(d, MkSel(pre)));
    const PeakFit fPnv = scan::FitPeakMain(scan::Fill(d, MkSel(pnv)));

    // ========================================================================
    // 3. bump scans
    // ========================================================================
    const std::vector<int> bSingle = {kOaLo, kOaHi, kOaLab, kMpLo, kMpHi, kMm4, kEta, kDz, kVtx, kVz};
    const BScan bs = BumpScan(d, dT, pre, bSingle, &pre0);
    printf("bump single at %.1f s\n", sw.RealTime()); sw.Continue();
    const Cuts bnv = BumpNomV();
    const std::vector<int> bNm1 = {kOaLab, kMpLo, kMpHi, kVtx, kOaLo, kOaHi, kDz, kMm4, kVz};
    const BScan bn = BumpScan(d, dT, bnv, bNm1);
    printf("bump N-1 at %.1f s\n", sw.RealTime()); sw.Continue();

    // ========================================================================
    // 4. cut flows
    // ========================================================================
    auto mkStep = [](const char* lab, const char* sub, const Cuts& c, bool applied) { Step s; s.lab = lab; s.sub = sub; s.c = c; s.applied = applied; return s; };
    std::vector<Step> fp, fb, ff;
    {
        Cuts c = pre0;
        fp.push_back(mkStep("pre without vz cut", "isBest #geq 0, PT3, z>-500", c, true));
        c.vzMin = -160; { Step s = mkStep("vz > -160 mm", "= pre; band -170, -150", c, true); s.band = true; s.lo = c; s.lo.vzMin = -170; s.hi = c; s.hi.vzMin = -150; fp.push_back(s); }
        c.oaLab = 55; { Step s = mkStep("oa_{lab} < 55#circ", "band 50, 60", c, true); s.band = true; s.lo = c; s.lo.oaLab = 50; s.hi = c; s.hi.oaLab = 60; fp.push_back(s); }
        c.mpLo = 0.30f; { Step s = mkStep("M(#pi#pi) > 0.30", "= peak nominal; band 0.28, 0.32", c, true); s.band = true; s.lo = c; s.lo.mpLo = 0.28f; s.hi = c; s.hi.mpLo = 0.32f; fp.push_back(s); }
        c.vtx = 2; { Step s = mkStep("both #pi in vertex fit", "= peak nom. + vertex; band #geq1 #pi, all 4", c, true); s.band = true; s.lo = c; s.lo.vtx = 1; s.hi = c; s.hi.vtx = 4; fp.push_back(s); }
        { Cuts x = c; x.oaLo = 2; fp.push_back(mkStep("+ OA > 2#circ", "not applied", x, false)); }
        { Cuts x = c; x.dzPi = 40; x.vtx = 0; fp.push_back(mkStep("#Deltaz_{#pi}<40 instead of vtx", "not applied", x, false)); }
        { Cuts x = c; x.mm4Hi = 2.85f; fp.push_back(mkStep("+ MM_{4} < 2.85", "not applied (flag)", x, false)); }
        { Cuts x = c; x.mpHi = 0.46f; fp.push_back(mkStep("+ M(#pi#pi) < 0.46", "not applied (sculpts)", x, false)); }
    }
    {
        Cuts c = pre0;
        fb.push_back(mkStep("pre without vz cut", "isBest #geq 0, PT3, z>-500", c, true));
        c.vzMin = -160; { Step s = mkStep("vz > -160 mm", "= pre = BUMP (frozen)", c, true); s.band = true; s.lo = c; s.lo.vzMin = -170; s.hi = c; s.hi.vzMin = -150; fb.push_back(s); }
        c.oaLab = 40; { Step s = mkStep("oa_{lab} < 40#circ", "band 35, 45", c, true); s.band = true; s.lo = c; s.lo.oaLab = 35; s.hi = c; s.hi.oaLab = 45; fb.push_back(s); }
        c.mpLo = 0.30f; { Step s = mkStep("M(#pi#pi) > 0.30", "band 0.28, 0.32", c, true); s.band = true; s.lo = c; s.lo.mpLo = 0.28f; s.hi = c; s.hi.mpLo = 0.32f; fb.push_back(s); }
        c.mpHi = 0.41f; { Step s = mkStep("M(#pi#pi) < 0.41", "= bump nominal; band 0.40, 0.42", c, true); s.band = true; s.lo = c; s.lo.mpHi = 0.40f; s.hi = c; s.hi.mpHi = 0.42f; fb.push_back(s); }
        c.vtx = 2; { Step s = mkStep("both #pi in vertex fit", "= bump nom. + vertex; band #geq1 #pi, all 4", c, true); s.band = true; s.lo = c; s.lo.vtx = 1; s.hi = c; s.hi.vtx = 4; fb.push_back(s); }
        { Cuts x = c; x.oaLo = 2; fb.push_back(mkStep("+ OA > 2#circ", "not applied", x, false)); }
        { Cuts x = c; x.dzPi = 40; x.vtx = 0; fb.push_back(mkStep("#Deltaz_{#pi}<40 instead of vtx", "not applied", x, false)); }
        { Cuts x = c; x.oaHi = 20; fb.push_back(mkStep("+ OA < 20#circ", "not applied", x, false)); }
    }
    {
        Cuts c = pre;
        ff.push_back(mkStep("pre", "= BUMP (frozen)", c, true));
        c.dzPi = 40; { Step s = mkStep("#Deltaz_{#pi} < 40 mm", "= DIRECT (frozen); band 30, 60", c, true); s.band = true; s.lo = c; s.lo.dzPi = 30; s.hi = c; s.hi.dzPi = 60; ff.push_back(s); }
        c.mm4Hi = 2.85f; { Step s = mkStep("MM_{4} < 2.85", "= PEAK (frozen); band 2.80, 2.90", c, true); s.band = true; s.lo = c; s.lo.mm4Hi = 2.80f; s.hi = c; s.hi.mm4Hi = 2.90f; ff.push_back(s); }
        { Cuts x = c; x.oaLo = 2; ff.push_back(mkStep("+ OA > 2#circ", "not applied", x, false)); }
        { Cuts x = c; x.best = 1; ff.push_back(mkStep("isBest == 1", "comparison", x, false)); }
    }
    auto runFlow = [&](const std::vector<Step>& st) {
        std::vector<Sel> sels; std::vector<Cuts> cs;
        for (const Step& s : st) {
            sels.push_back(MkSel(s.c)); cs.push_back(s.c);
            sels.push_back(MkSel(s.band ? s.lo : s.c)); cs.push_back(s.band ? s.lo : s.c);
            sels.push_back(MkSel(s.band ? s.hi : s.c)); cs.push_back(s.band ? s.hi : s.c);
        }
        const std::vector<Spectra> sp = scan::FillMany(d, sels);
        const std::vector<BumpRes> br = EvalBump(d, cs, true, true, &dT);
        std::vector<StepRes> out(st.size());
        for (size_t i = 0; i < st.size(); ++i) {
            out[i].p = EvalPeak(sp[3 * i], true, false); out[i].pLo = EvalPeak(sp[3 * i + 1], false, false); out[i].pHi = EvalPeak(sp[3 * i + 2], false, false);
            out[i].b = br[3 * i]; out[i].bLo = br[3 * i + 1]; out[i].bHi = br[3 * i + 2];
        }
        return out;
    };
    const std::vector<StepRes> rp = runFlow(fp), rb = runFlow(fb), rf = runFlow(ff);
    printf("cut flows at %.1f s\n", sw.RealTime()); sw.Continue();

    // ========================================================================
    // 5. vertex: signal per vz slice (base: pre without vz cut)
    // ========================================================================
    const int nz = 11;
    const double zl[nz] = {-500, -200, -170, -160, -150, -140, -130, -120, -110, -100, -90};
    const double zh[nz] = {-200, -170, -160, -150, -140, -130, -120, -110, -100, -90, 0};
    std::vector<PeakFit> zS(nz); std::vector<PeakFit> zT(nz); std::vector<scan::WinSum> zB(nz);
    {
        std::vector<Sel> sels, ts;
        for (int i = 0; i < nz; ++i) {
            const float a = zl[i], b = zh[i];
            sels.push_back([=](const Row& r) { return Pass(pre0, r) && r.vz > a && r.vz <= b; });
            ts.push_back([=](const Row& r) { return Tagged(r) && Pass(pre0, r) && r.vz > a && r.vz <= b; });
        }
        const std::vector<Spectra> sp = scan::FillMany(d, sels);
        const std::vector<Spectra> st = scan::FillMany(dT, ts, TagOpt());
        for (int i = 0; i < nz; ++i) { zS[i] = scan::FitPeakMain(sp[i]); zT[i] = scan::FitPeakMain(st[i], TagCfg()); zB[i] = sp[i].win[1]; }
    }
    const PeakFit fAll0 = scan::FitPeakMain(scan::Fill(d, MkSel(pre0)));
    double tAll0 = 0, nbAll0 = 0;
    {
        std::vector<Cuts> one(1, pre0);
        const std::vector<BumpRes> b = EvalBump(d, one, false, true, &dT);
        tAll0 = b[0].T.S; nbAll0 = b[0].win.n[0];
    }

    // ========================================================================
    // tables
    // ========================================================================
    Out o;
    o.Open(std::string(kTabOut) + "/opt_scans.txt");
    o.P("# opt_scans.txt  --  written by plots/exp_report/opt_scans.C\n");
    o.P("# EXPERIMENTAL DATA. Rows = stored particle combinations, cor_ok == 1, corrected momenta; FAT: vertex z > -500 mm, PT3.\n");
    o.P("# pre = isBest >= 0 && vz > -160 mm. CB = 2 sqrt(N++ N--), k = 1. P_main: N - CB, frozen Gaussian 0.5439 / 5.4 MeV + pol2, 0.5225-0.6100.\n");
    o.P("# S = whole Gaussian; F = S/stat (compares selections, NOT a significance); S/B in the data window 0.53364-0.55416.\n");
    o.P("# eps_S = S_pass/(S_pass + S_fail), eps_B = kept share of (pol2 + CB) in the window; z = (eps_S - sqrt(eps_B))/err(eps_S).\n");
    o.P("# f0/f1 = the two halves of the event groups. Peak flag = 1: flagged by the sculpting test (opt_sculpt.txt); bump flag = 1: summed edge fit p < 1e-3.\n");
    o.P("# BUMP: E = edge-step visibility (rows / 2.5 MeV, NOT a yield), N_b = rows in the bump window 0.4550-0.5225, T = ECAL-tagged eta;\n");
    o.P("#   effE = E/E_ref, effT = T/T_ref, effN = N_b/N_b,ref (ref = the same point with this cut open), gain = eff/sqrt(effN).\n#\n");
    auto peakTable = [&](const char* title, const PScan& s) {
        o.P("\n%s\n", title);
        for (size_t iv = 0; iv < s.vars.size(); ++iv) {
            const int v = s.vars[iv];
            o.P("%-11s %6s %4s | %7s %6s | %6s %5s %6s %6s %8s | %6s %5s | %6s %5s %6s %6s %6s | %6s %5s %6s | %6s %5s %6s\n", gV[v].fam.c_str(), "value", "flag",
                "N", "CB/N", "S", "err", "F", "S/B", "chi2/ndf", "S_fail", "err", "eps_S", "err", "eps_B", "sqrt", "z", "S f0", "err", "F f0", "S f1", "err", "F f1");
            for (size_t i = 0; i < gV[v].val.size(); ++i) {
                const PPt& q = s.pt[iv][i];
                o.P("%-11s %6s %4d | %7.0f %6.3f | %6.0f %5.0f %6.2f %6.4f %5.1f/%-2d | %6.0f %5.0f | %6.3f %5.3f %6.3f %6.3f %6.2f | %6.0f %5.0f %6.2f | %6.0f %5.0f %6.2f\n",
                    gV[v].fam.c_str(), gV[v].key[i].c_str(), gV[v].flag[i], q.pass.N, q.pass.N > 0 ? q.pass.CB / q.pass.N : 0, q.pass.S, q.pass.eS, q.pass.F,
                    q.pass.SoB, q.pass.chi2, q.pass.ndf, q.fail.S, q.fail.eS, q.epsS, q.eEpsS, q.epsB, std::sqrt(std::max(q.epsB, 0.)), q.z,
                    q.f0.S, q.f0.eS, q.f0.F, q.f1.S, q.f1.eS, q.f1.F);
            }
            o.P("\n");
        }
    };
    auto bumpTable = [&](const char* title, const BScan& s) {
        o.P("\n%s\n", title);
        for (size_t iv = 0; iv < s.vars.size(); ++iv) {
            const int v = s.vars[iv];
            o.P("%-11s %6s %4s | %6s %4s %8s | %7s | %5s %5s %5s | %6s %4s | %5s %5s %5s %5s %5s | %5s %5s\n", gV[v].fam.c_str(), "value", "flag",
                "E", "err", "chi2/ndf", "N_b", "effE", "err", "effN", "T", "err", "effT", "err", "gainT", "err", "gainE", "err", "Q");
            for (size_t i = 0; i < gV[v].val.size(); ++i) {
                const BPt& q = s.pt[iv][i];
                o.P("%-11s %6s %4d | %6.0f %4.0f %4.0f/%-3d | %7.0f | %5.3f %5.3f %5.3f | %6.0f %4.0f | %5.3f %5.3f %5.3f %5.3f %5.3f | %5.3f %5.2f\n",
                    gV[v].fam.c_str(), gV[v].key[i].c_str(), (int)BumpFlag(q), q.r.E, q.r.eE, q.r.chi2, q.r.ndf, q.r.win.n[0], q.effE, q.eEffE, q.effN,
                    q.r.T.S, q.r.T.eS, q.effT, q.eEffT, q.gT, q.egT, q.gE, q.egE, q.r.Q);
            }
            o.P("\n");
        }
    };
    o.P("PRE: S = %.0f +- %.0f, F = %.2f ; PEAK NOMINAL + VERTEX: S = %.0f +- %.0f, F = %.2f\n", fPre.S, fPre.eS, fPre.F, fPnv.S, fPnv.eS, fPnv.F);
    peakTable("==== A. ETA PEAK, SINGLE-CUT SCANS FROM pre (vz threshold from pre without vz > -160) ====", ps);
    peakTable("==== B. ETA PEAK, N-1 SCANS AROUND peak nominal + vertex (oa_lab < 55, M(pipi) > 0.30, both pions in vertex fit) ====", pn);
    bumpTable("==== C. BUMP, SINGLE-CUT SCANS FROM pre (vz threshold from pre without vz > -160) ====", bs);
    bumpTable("==== D. BUMP, N-1 SCANS AROUND bump nominal + vertex (oa_lab < 40, 0.30 < M(pipi) < 0.41, both pions in vertex fit) ====", bn);
    auto flowTable = [&](const char* title, const std::vector<Step>& st, const std::vector<StepRes>& r) {
        o.P("\n%s\n", title);
        o.P("%-34s %-8s | %7s %7s | %6s %5s %5s | %6s %6s | %5s-%5s | %6s %4s %4s | %7s | %5s %5s | %5s-%5s\n", "step", "applied", "N_win", "CB",
            "S", "stat", "sprd", "S/B", "F", "F lo", "F hi", "E", "stat", "sprd", "N_b", "T", "err", "E lo", "E hi");
        for (size_t i = 0; i < st.size(); ++i) {
            const StepRes& q = r[i];
            o.P("%-34s %-8s | %7.0f %7.0f | %6.0f %5.0f %5.0f | %6.4f %6.2f | %5.2f-%5.2f | %6.0f %4.0f %4.0f | %7.0f | %5.0f %5.0f | %5.0f-%5.0f   %s\n",
                st[i].lab.c_str(), st[i].applied ? "yes" : "no", q.p.f.N, q.p.f.CB, q.p.f.S, q.p.f.eS, q.p.spr.half, q.p.f.SoB, q.p.f.F,
                std::min(q.pLo.f.F, q.pHi.f.F), std::max(q.pLo.f.F, q.pHi.f.F), q.b.E, q.b.eE, q.b.spreadE, q.b.win.n[0], q.b.T.S, q.b.T.eS,
                std::min(q.bLo.E, q.bHi.E), std::max(q.bLo.E, q.bHi.E), Describe(st[i].c).c_str());
        }
    };
    flowTable("==== E. CUT FLOW, ETA PEAK (data-optimised) ====", fp, rp);
    flowTable("==== F. CUT FLOW, BUMP (data-optimised) ====", fb, rb);
    flowTable("==== G. CUT FLOW, FROZEN SIMULATION SELECTIONS ====", ff, rf);
    o.P("\n==== H. SIGNAL PER VERTEX-Z SLICE (base: pre without vz > -160; no other cut) ====\n");
    o.P("# N = rows in the peak window, S = P_main; N_b = rows in the bump window, T = tagged eta; ratio = (S/N)/(S/N)_all, (T/N_b)/(T/N_b)_all\n");
    o.P("%-12s | %7s %6s | %6s %5s %6s %6s | %7s %5s %4s %6s %6s\n", "vz slice", "N", "CB/N", "S", "err", "ratio", "err", "N_b", "T", "err", "ratio", "err");
    const double sn0 = fAll0.S / fAll0.N, tn0 = tAll0 / nbAll0;
    for (int i = 0; i < nz; ++i) {
        const double r1 = zS[i].N > 0 ? (zS[i].S / zS[i].N) / sn0 : 0, e1 = zS[i].N > 0 ? (zS[i].eS / zS[i].N) / sn0 : 0;
        const double r2 = zB[i].n[0] > 0 ? (zT[i].S / zB[i].n[0]) / tn0 : 0, e2 = zB[i].n[0] > 0 ? (zT[i].eS / zB[i].n[0]) / tn0 : 0;
        o.P("%5.0f..%-5.0f | %7.0f %6.3f | %6.0f %5.0f %6.2f %6.2f | %7.0f %5.0f %4.0f %6.2f %6.2f\n", zl[i], zh[i], zS[i].N, zS[i].N > 0 ? zS[i].CB / zS[i].N : 0,
            zS[i].S, zS[i].eS, r1, e1, zB[i].n[0], zT[i].S, zT[i].eS, r2, e2);
    }
    o.P("\nrun time %.1f s\n", sw.RealTime());
    o.Close();

    // ---- LaTeX cut-flow table ----
    {
        FILE* t = fopen((std::string(kTabOut) + "/opt_cutflow.tex").c_str(), "w");
        fprintf(t, "%% opt_cutflow.tex -- written by plots/exp_report/opt_scans.C (numbers: tables_exp/opt_scans.txt)\n");
        fprintf(t, "\\begin{table}[H]\n\\centering\\footnotesize\n\\setlength{\\tabcolsep}{3pt}\n");
        fprintf(t, "\\caption{Cut flows on the data, all stored combinations (isBest $\\geq 0$), corrected momenta, PT3, FAT vertex $z>-500$~mm; cuts added from top to bottom inside each block, "
                   "rows marked ``no'' are not applied (each added to the last applied step). $N$, CB in the data peak window 0.53364--0.55416~GeV/$c^2$; $S\\pm$stat$\\pm$model spread from $P_\\mathrm{main}$; "
                   "$F=S/\\sigma_S^\\mathrm{stat}$; $E\\pm$stat$\\pm$spread = edge-step visibility of the bump (not a yield); $T$ = ECAL-tagged $\\eta$. CB $=2\\sqrt{N_{++}N_{--}}$, $k=1$. Rows, not normalised to luminosity.}\n");
        fprintf(t, "\\label{tab:opt_cutflow}\n\\begin{tabular}{llrrrrrrr}\n\\toprule\n");
        fprintf(t, "Step & applied & $N$ & CB & $S\\pm$stat$\\pm$spread & $S/B$ & $F$ & $E\\pm$stat$\\pm$spread & $T\\pm$stat \\\\\n\\midrule\n");
        auto block = [&](const char* name, const std::vector<Step>& st, const std::vector<StepRes>& r) {
            fprintf(t, "\\multicolumn{9}{l}{\\textbf{%s}} \\\\\n", name);
            for (size_t i = 0; i < st.size(); ++i) {
                std::string lab = st[i].lab;
                // TLatex -> LaTeX for the few symbols used
                auto rep = [&](const std::string& a, const std::string& b) { size_t p; while ((p = lab.find(a)) != std::string::npos) lab.replace(p, a.size(), b); };
                rep("#circ", "$^\\circ$"); rep("#Deltaz_{#pi}", "$\\Delta z_\\pi$"); rep("#pi#pi", "$\\pi\\pi$"); rep("#pi", "$\\pi$"); rep("oa_{lab}", "$\\alpha_\\mathrm{lab}$");
                rep("MM_{4}", "MM$_4$"); rep("#geq", "$\\geq$");
                const StepRes& q = r[i];
                fprintf(t, "%s & %s & %.0f & %.0f & $%.0f\\pm%.0f\\pm%.0f$ & %.3f & %.1f & $%.0f\\pm%.0f\\pm%.0f$ & $%.0f\\pm%.0f$ \\\\\n", lab.c_str(),
                        st[i].applied ? "yes" : "no", q.p.f.N, q.p.f.CB, q.p.f.S, q.p.f.eS, q.p.spr.half, q.p.f.SoB, q.p.f.F, q.b.E, q.b.eE, q.b.spreadE, q.b.T.S, q.b.T.eS);
            }
        };
        block("$\\eta$ peak, data-optimised", fp, rp); fprintf(t, "\\midrule\n");
        block("bump, data-optimised", fb, rb); fprintf(t, "\\midrule\n");
        block("frozen simulation selections", ff, rf);
        fprintf(t, "\\bottomrule\n\\end{tabular}\n\\end{table}\n");
        fclose(t);
    }

    // ========================================================================
    // figures
    // ========================================================================
    DrawPeakF("opt_peak_single", "#eta peak, single-cut scans from pre (data, all combinations isBest #geq 0, P_{main}, F = S/#sigma_{S} statistical only)", ps, fPre.F, Form("F(pre) = %.1f", fPre.F));
    DrawPeakEff("opt_peak_single_eff", "#eta peak, single-cut scans from pre: signal kept against #sqrt{background kept} (data, isBest #geq 0)", ps);
    DrawPeakF("opt_peak_nm1", "#eta peak, N-1 scans around peak nominal + vertex (data, isBest #geq 0, P_{main})", pn, fPnv.F, Form("F(peak nominal + vertex) = %.1f", fPnv.F));
    DrawBump("opt_bump_single", "bump, single-cut scans from pre (data, isBest #geq 0): edge-step visibility E and tagged #eta T", bs);
    DrawBump("opt_bump_nm1", "bump, N-1 scans around bump nominal + vertex (data, isBest #geq 0)", bn);

    // ---- cut-flow figure: three columns (peak, bump, frozen), three rows ----
    {
        // verifier change 2026-11-14: canvas 1800 x 1250 -> 1300 x 1000 (labels were 4 pt at the text width)
        TCanvas* c = new TCanvas("c_flow", "", 1300, 1000);
        const std::vector<Step>* st[3] = {&fp, &fb, &ff};
        const std::vector<StepRes>* rr[3] = {&rp, &rb, &rf};
        const char* ttl[3] = {"(a) #eta peak, data-optimised", "(b) bump, data-optimised", "(c) frozen on simulation"};
        const double xw[4] = {0.0, 0.40, 0.80, 1.0};
        double nmax = 0;
        for (int k = 0; k < 3; ++k) for (const StepRes& q : *rr[k]) nmax = std::max(nmax, std::max(q.p.f.N, q.b.win.n[0]));
        for (int k = 0; k < 3; ++k) {
            const int ns = (int)st[k]->size();
            const double x0 = k == 0 ? 0.0 : (k == 1 ? 0.40 : 0.75), x1 = k == 0 ? 0.40 : (k == 1 ? 0.75 : 1.0);
            (void)xw;
            for (int row = 0; row < 3; ++row) {
                c->cd();
                const double ytop = 0.95 - row * 0.255, ybot = ytop - 0.255 - (row == 2 ? 0.10 : 0.0);
                TPad* p = new TPad(Form("pf%d%d", k, row), "", x0, row == 2 ? 0.0 : ybot, x1, ytop);
                p->SetLeftMargin(std::min(0.35, 100.0 / ((x1 - x0) * c->GetWw()))); p->SetRightMargin(0.02);   // 100 px for the y labels and title
                p->SetTopMargin(row == 0 ? 0.10 : 0.03); p->SetBottomMargin(row == 2 ? 0.42 : 0.03);
                if (row == 0) p->SetLogy();
                p->Draw(); p->cd();
                TH1D* fr = new TH1D(Form("ff%d%d", k, row), "", ns, 0, ns); fr->SetDirectory(nullptr);
                for (int i = 0; i < ns; ++i) fr->GetXaxis()->SetBinLabel(i + 1, (*st[k])[i].lab.c_str());
                AxPx(fr, 19, 21, 1.0, 2.4);
                if (row < 2) fr->GetXaxis()->SetLabelSize(0);
                else { fr->LabelsOption("v", "X"); fr->GetXaxis()->SetLabelSize(18); }
                double y1 = 0, y2 = 1;
                const char* yt = "";
                if (row == 0) { y1 = 30; y2 = 3 * nmax; yt = "rows in the window"; }
                if (row == 1) { y1 = 0; y2 = 7600; yt = "S  and  4 #times E [rows]"; }
                if (row == 2) { y1 = 0.0; y2 = 18; yt = "F  and  10 #times gain"; }
                fr->SetMinimum(y1); fr->SetMaximum(y2); fr->GetYaxis()->SetTitle(yt);
                fr->Draw("AXIS");
                for (int i = 0; i < ns; ++i)
                    if (!(*st[k])[i].applied) { TBox* b = new TBox(i, row == 0 ? y1 : y1, i + 1, y2); b->SetFillColor(TColor::GetColor("#fde9d3")); b->SetLineWidth(0); b->Draw(); }
                const bool bump = (k == 1);
                if (row == 0) {
                    TGraphErrors* gN = GE(kBlack, 20, 1.0); TGraphErrors* gC = GE(kRed + 1, 21, 1.0); TGraphErrors* gNb = GE(kBlack, 24, 1.0); TGraphErrors* gCb = GE(kRed + 1, 25, 1.0);
                    for (int i = 0; i < ns; ++i) {
                        const StepRes& q = (*rr[k])[i];
                        AddPt(gN, i + 0.4, q.p.f.N, 0); AddPt(gC, i + 0.4, q.p.f.CB, 0);
                        AddPt(gNb, i + 0.6, q.b.win.n[0], 0); AddPt(gCb, i + 0.6, q.b.win.cb, 0);
                    }
                    gN->Draw("P SAME"); gC->Draw("P SAME"); gNb->Draw("P SAME"); gCb->Draw("P SAME");
                    TxtPx(0.5 * (1 + p->GetLeftMargin()), 0.92, ttl[k], 22, kBlack, 21);
                    if (k == 0) {
                        TLegend* l = new TLegend(p->GetLeftMargin() + 0.01, 0.04, 0.99, 0.30); l->SetNColumns(2);
                        l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextFont(43); l->SetTextSize(18);
                        l->AddEntry(gN, "N, peak window", "p"); l->AddEntry(gNb, "N_{b}, bump window", "p");
                        l->AddEntry(gC, "CB, peak window", "p"); l->AddEntry(gCb, "CB, bump window", "p"); l->Draw();
                    }
                } else if (row == 1) {
                    TGraphErrors* gS = GE(kBlue + 1, 20, 1.1); TGraphErrors* gSt = GE(kBlue + 1, 1, 0.1); gSt->SetLineWidth(1);
                    TGraphErrors* gE = GE(kViolet + 1, 21, 1.0); TGraphErrors* gEt = GE(kViolet + 1, 1, 0.1); gEt->SetLineWidth(1);
                    for (int i = 0; i < ns; ++i) {
                        const StepRes& q = (*rr[k])[i];
                        if ((*st[k])[i].band) {
                            TBox* b = new TBox(i + 0.15, std::min(q.pLo.f.S, q.pHi.f.S), i + 0.45, std::max(q.pLo.f.S, q.pHi.f.S)); b->SetFillColor(TColor::GetColor("#c9d0f2")); b->SetLineWidth(0); b->Draw();
                            TBox* e = new TBox(i + 0.55, 4 * std::min(q.bLo.E, q.bHi.E), i + 0.85, 4 * std::max(q.bLo.E, q.bHi.E)); e->SetFillColor(TColor::GetColor("#f6caca")); e->SetLineWidth(0); e->Draw();
                        }
                        AddPt(gS, i + 0.3, q.p.f.S, q.p.f.eS); AddPt(gSt, i + 0.3, q.p.f.S, std::hypot(q.p.f.eS, q.p.spr.half));
                        AddPt(gE, i + 0.7, 4 * q.b.E, 4 * q.b.eE); AddPt(gEt, i + 0.7, 4 * q.b.E, 4 * std::hypot(q.b.eE, q.b.spreadE));
                    }
                    gSt->Draw("Z SAME"); gEt->Draw("Z SAME"); gS->Draw("PZ SAME"); gE->Draw("PZ SAME");
                    if (k == 0) {
                        TLegend* l = new TLegend(p->GetLeftMargin() + 0.01, 0.04, 0.99, 0.26); l->SetNColumns(2);
                        l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextFont(43); l->SetTextSize(18);
                        l->AddEntry(gS, "#eta peak S", "pe"); l->AddEntry(gE, "4 #times bump step E", "pe"); l->Draw();
                        TxtPx(p->GetLeftMargin() + 0.02, 0.90, "thick: stat; thin: #oplus model; boxes: cut #pm one step", 17, kGray + 2);
                    }
                } else {
                    TGraphErrors* gF = GE(kBlue + 1, 20, 1.1); TGraphErrors* gG = GE(kViolet + 1, 21, 1.0);
                    TGraphErrors* gGT = GE(kAzure + 7, 25, 1.0);
                    const double gT0 = 1.0;
                    for (int i = 0; i < ns; ++i) {
                        const StepRes& q = (*rr[k])[i];
                        if ((*st[k])[i].band) {
                            TBox* b = new TBox(i + 0.15, std::min(q.pLo.f.F, q.pHi.f.F), i + 0.45, std::max(q.pLo.f.F, q.pHi.f.F)); b->SetFillColor(TColor::GetColor("#c9d0f2")); b->SetLineWidth(0); b->Draw();
                        }
                        AddPt(gF, i + 0.3, q.p.f.F, 0);
                        // gain of the tagged eta and of E relative to the first step of the block
                        const StepRes& q0 = (*rr[k])[k == 2 ? 0 : 1];
                        const double effN = q0.b.win.n[0] > 0 ? q.b.win.n[0] / q0.b.win.n[0] : 0;
                        if (effN > 0 && q0.b.T.S > 0) AddPt(gGT, i + 0.7, 10 * gT0 * (q.b.T.S / q0.b.T.S) / std::sqrt(effN), 0);
                        if (effN > 0 && q0.b.E > 0) AddPt(gG, i + 0.7, 10 * (q.b.E / q0.b.E) / std::sqrt(effN), 0);
                    }
                    HLine(0, ns, 10, kGray + 1, 2, 1);
                    gF->Draw("P SAME"); gG->Draw("P SAME"); gGT->Draw("P SAME");
                    if (k == 0) {
                        TLegend* l = new TLegend(p->GetLeftMargin() + 0.01, 0.45, 0.99, 0.60); l->SetNColumns(3);
                        l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextFont(43); l->SetTextSize(18);
                        l->AddEntry(gF, "F = S/#sigma_{S}", "p"); l->AddEntry(gG, "10 #times gain_{E}", "p"); l->AddEntry(gGT, "10 #times gain_{T}", "p"); l->Draw();
                    }
                }
                fr->Draw("AXIS SAME");
            }
        }
        c->cd();
        TxtPx(0.01, 0.995, "Cut flows on the data (all combinations isBest #geq 0, corrected momenta). Shaded columns: cut not applied (added to the last applied step).", 17, kGray + 3, 13);
        TxtPx(0.01, 0.973, "Gains relative to pre (bump: N_{b} = rows in the bump window).", 17, kGray + 3, 13);
        Save(c, "opt_cutflow");
    }

    // ---- vertex figure ----
    {
        // verifier change 2026-11-14: canvas 1500 x 1000 -> 1150 x 800 (larger labels at the text width)
        TCanvas* c = new TCanvas("c_vtx", "", 1150, 800);
        // (a) signal share per vz slice
        TPad* pa = Cell(c, 0, 2, 2, 0.97);
        pa->SetBottomMargin(0.16); pa->SetTopMargin(0.10);
        TH1D* fa = Frame("fva", -260, 0, -0.3, 3.0, "vertex z [mm]", "(S/N) / (S/N)_{all z}");
        AxPx(fa, 18, 20, 1.2, 1.6);
        TGraphErrors* g1 = GE(kBlue + 1, 20, 1.1); TGraphErrors* g2 = GE(kAzure + 7, 25, 1.1); TGraphErrors* g3 = GE(kRed + 1, 1, 0.1);
        for (int i = 0; i < nz; ++i) {
            const double xl = std::max(zl[i], -260.), xh = std::min(zh[i], 0.), xc = 0.5 * (xl + xh);
            if (zS[i].N > 0) AddPt(g1, xc - 1.5, (zS[i].S / zS[i].N) / sn0, (zS[i].eS / zS[i].N) / sn0, 0.5 * (xh - xl));
            if (zB[i].n[0] > 0) AddPt(g2, xc + 1.5, (zT[i].S / zB[i].n[0]) / tn0, (zT[i].eS / zB[i].n[0]) / tn0, 0);
            if (zS[i].N > 0) AddPt(g3, xc, zS[i].CB / zS[i].N * 2, 0, 0.5 * (xh - xl));
        }
        HLine(-260, 0, 1, kGray + 1, 2, 1); HLine(-260, 0, 0, kGray + 1, 1, 1);
        VLine(-160, -0.3, 3.0, kGreen + 2, 2, 2);
        g3->SetLineWidth(3); g3->Draw("Z SAME"); g1->Draw("PZ SAME"); g2->Draw("PZ SAME");
        TxtPx(0.17, 0.95, "(a) signal per vertex-z slice (pre without vz cut)", 19, kBlack, 13);
        { TLegend* l = new TLegend(0.19, 0.70, 0.80, 0.885); l->SetMargin(0.12); l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextFont(43); l->SetTextSize(15);
          l->AddEntry(g1, "#eta peak S per peak-window row", "pe"); l->AddEntry(g2, "tagged #eta T per bump-window row", "pe");
          l->AddEntry(g3, "2 #times CB/N, peak window", "l"); l->Draw(); }
        // (b) vz threshold, (c) vertex-fit categories, (d) dz_pi threshold: eps_S and sqrt(eps_B), gain_T
        const int vv[3] = {kVz, kVtx, kDz};
        const char* tt[3] = {"(b) lower cut on vz (from pre without vz cut)", "(c) tracks in the vertex fit (from pre)", "(d) #Deltaz_{#pi} = max|z_{#pi} - v_{z}| (from pre)"};
        for (int q = 0; q < 3; ++q) {
            const int v = vv[q];
            int iv = -1; for (size_t k = 0; k < ps.vars.size(); ++k) if (ps.vars[k] == v) iv = (int)k;
            int jb = -1; for (size_t k = 0; k < bs.vars.size(); ++k) if (bs.vars[k] == v) jb = (int)k;
            TPad* p = Cell(c, q + 1, 2, 2, 0.97);
            p->SetBottomMargin(0.21); p->SetTopMargin(0.10);
            TH1D* fr = VarFrame(Form("fvx%d", q), v, 0.0, 1.35, "fraction kept / gain");
            fr->GetXaxis()->SetTitleOffset(1.45);
            HLine(0, gV[v].val.size(), 1, kGray + 1, 2, 1);
            TGraphErrors* gS = GE(kBlue + 1, 20, 1.1); TGraphErrors* gB = GE(kBlack, 1, 0.1); TGraphErrors* gT = GE(kAzure + 7, 25, 1.1);
            TGraphErrors* gF = GE(kMagenta + 2, 22, 1.2);
            for (size_t i = 0; i < gV[v].val.size(); ++i) {
                AddPt(gS, i + 0.4, ps.pt[iv][i].epsS, ps.pt[iv][i].eEpsS);
                AddPt(gB, i + 0.5, std::sqrt(std::max(ps.pt[iv][i].epsB, 0.)), 0);
                AddPt(gT, i + 0.6, bs.pt[jb][i].gT, bs.pt[jb][i].egT);
                AddPt(gF, i + 0.5, ps.pt[iv][i].pass.F / ps.pt[iv][NoneIdx(v)].pass.F, 0);
            }
            gB->SetLineWidth(2); gB->Draw("L SAME"); gS->Draw("PZ SAME"); gT->Draw("PZ SAME"); gF->Draw("P SAME");
            fr->Draw("AXIS SAME");
            TxtPx(0.17, 0.95, tt[q], 19, kBlack, 13);
            if (q == 0) {
                TLegend* l = new TLegend(0.20, 0.20, 0.95, 0.42); l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextFont(43); l->SetTextSize(15); l->SetNColumns(2);
                l->AddEntry(gS, "#varepsilon_{S}, #eta peak (pass/fail fits)", "pe"); l->AddEntry(gB, "#sqrt{#varepsilon_{B}}, peak window", "l");
                l->AddEntry(gT, "gain_{T}, tagged #eta (bump)", "pe"); l->AddEntry(gF, "F / F(no cut), #eta peak", "p"); l->Draw();
            }
        }
        c->cd();
        TxtPx(0.01, 0.995, "Vertex cuts on the data (all combinations isBest #geq 0, corrected momenta, PT3, FAT z > -500 mm)", 18, kGray + 3, 13);
        Save(c, "opt_vertex");
    }
    printf("opt_scans: done in %.1f s\n", sw.RealTime());
}
