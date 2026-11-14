// ============================================================================
// s5_acc.C  --  steps 1 and 2 of the report's section 5: the eta signals at generator level and what HADES keeps
// ============================================================================
// HADES p+p, T = 4.53 GeV, pi+ pi- e+ e-; SMASH simulation ONLY (no measured data is read, no CB exists or is added).
// Signals (src/sim_truth.h four-particle classes): 2 bump eta -> pi+pi-pi0, pi0 -> e+e-gamma (std); 3 bump with a converted
// pi0 photon (std); 4 peak eta -> pi+pi-gamma, photon converted (std); 5 peak eta -> pi+pi-e+e- (eta set); 6 X17
// (eta -> pi+pi- X17, X17 -> e+e-, eta set; hypothesis, W = 2 BR per row, never summed).
// Weights: rows owned in mode A, std + lep + eta (SimOpt::sets = 7), w (PDG), physics weights of the eta decays ON with the
// main variant (simphys::kPhysMain = V3 + KLOE bump + conversion F_V + background rules). Base selection "pre"
// (SimEngine kPre: isBest >= 0, vz > -500, PT3 bit, vz > -160 mm, corrected momenta, cor_ok).
//
// Part bits (default 15):
//   1 generator level: the PLUTO sample of the direct decay (/hdd3/.../generator/eta_pipiee_gen_kin_1000files.txt.gz,
//     1 181 801 decays from 1000 of the 100 000 files; columns M(ee), M(pipi), cos th*_e+, cos th*_pi+, M(eta)) with the
//     V3 weight; analytic production densities of the bump (HGeant phase space x KLOE-2 Dalitz plot), of the conversion
//     peak (phase space x |F_V|^2 p*^2 k^2), of the pi0 Dalitz pair (Kroll-Wada) and of X17 (3-body phase space ASSUMED);
//     absolute production cross sections from N_eta = 2 x 100 x N_sample (every eta of the eta set is forced, q = 1/2)
//   2 acceptance x efficiency of "pre" per produced decay, total and against true M(ee), M(pipi), OA, p(e+e-)
//   4 resolution and scale (M4, M(ee), M(pipi): reconstructed corrected minus true), position and width of the
//     simulated peak with the physics weights and the check of the fixed PEAK window 0.5333-0.5587
//   8 reconstructed OA against M(ee) and the OA separation of direct decay, conversion peak, X17 and bump Dalitz pairs
// Tables plots/eta_sim/tables/s5_acc_*.txt, figures plots/eta_sim/figs/s5_acc_*.{png,pdf} (pdf via EPS + ps2pdf -dEPSCrop).
//
// Run from /home/przygoda/HADES/FAT/FAT (about 3 min, about 8 GB while loading, one process):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/aclic_s5_acc", kTRUE)' \
//        'plots/eta_sim/s5_acc.C+(15)' 2>&1 | tee /hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/s5_acc.log
// All sums in double; histograms TH1D/TH2D; errors of accepted sums from event groups (SimEngine).
// ============================================================================
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <functional>
#include <map>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TColor.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMarker.h"
#include "TMath.h"
#include "TPad.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TSystem.h"

#include "SimEngine.h"
#include "sim_eta_physics.h"

namespace s5acc {
using namespace simeng;
using namespace simphys;

// ---------------------------------------------------------------- output helpers
FILE* gF = nullptr;
void P(const char* fmt, ...) {
    char buf[32768];
    va_list ap; va_start(ap, fmt); vsnprintf(buf, sizeof buf, fmt, ap); va_end(ap);
    fputs(buf, stdout);
    if (gF) fputs(buf, gF);
}
void OpenTab(const char* name) {
    const std::string p = std::string(kTabDir) + "/" + name;
    gF = fopen(p.c_str(), "w");
    if (!gF) printf("cannot write %s\n", p.c_str());
}
void CloseTab() { if (gF) fclose(gF); gF = nullptr; }
int gCnt = 0;
std::string Nm(const char* s) { return Form("%s_%d", s, gCnt++); }
// PNG + vector PDF with a tight page (EPS, then ps2pdf -dEPSCrop, EPS deleted)
void SaveFig(TCanvas* c, const char* name) {
    const std::string b = std::string(kFigDir) + "/" + name;
    c->SaveAs((b + ".png").c_str());
    c->SaveAs((b + ".eps").c_str());
    gSystem->Exec(Form("ps2pdf -dEPSCrop %s.eps %s.pdf && rm -f %s.eps", b.c_str(), b.c_str(), b.c_str()));
}
void Txt(double x, double y, const std::string& s, double size = 0.05, int col = 1, int align = 11) {
    TLatex t; t.SetNDC(); t.SetTextSize(size); t.SetTextColor(col); t.SetTextAlign(align); t.SetTextFont(42); t.DrawLatex(x, y, s.c_str());
}
void Style() {
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
    gStyle->SetLabelSize(0.05, "XYZ"); gStyle->SetTitleSize(0.055, "XYZ");
    gStyle->SetTitleOffset(1.10, "X"); gStyle->SetTitleOffset(1.25, "Y");
    gStyle->SetEndErrorSize(0); gStyle->SetTextFont(42); gStyle->SetLabelFont(42, "XYZ"); gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLegendFont(42);
    gStyle->SetPalette(kBird);
}
void PadM(double l = 0.15, double r = 0.04, double b = 0.15, double t = 0.09) {
    gPad->SetLeftMargin(l); gPad->SetRightMargin(r); gPad->SetBottomMargin(b); gPad->SetTopMargin(t);
}
TLegend* Leg(double x1, double y1, double x2, double y2, double ts = 0.045) {
    TLegend* l = new TLegend(x1, y1, x2, y2); l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(ts); l->SetTextFont(42);
    return l;
}
TH1D* H1v(const char* n, const std::vector<double>& e) {
    TH1D* h = new TH1D(Nm(n).c_str(), "", (int)e.size() - 1, e.data()); h->Sumw2(); h->SetDirectory(nullptr); return h;
}
TH1D* H1(const char* n, int nb, double lo, double hi) {
    TH1D* h = new TH1D(Nm(n).c_str(), "", nb, lo, hi); h->Sumw2(); h->SetDirectory(nullptr); return h;
}
void Sty(TH1* h, int col, int lw = 3, int ls = 1, int ms = 0, double msz = 1.2) {
    h->SetLineColor(col); h->SetLineWidth(lw); h->SetLineStyle(ls); h->SetMarkerColor(col);
    if (ms) { h->SetMarkerStyle(ms); h->SetMarkerSize(msz); }
}
void Sty(TGraph* g, int col, int ms, double msz = 1.3, int lw = 2, int ls = 1) {
    g->SetLineColor(col); g->SetMarkerColor(col); g->SetMarkerStyle(ms); g->SetMarkerSize(msz); g->SetLineWidth(lw);
    g->SetLineStyle(ls);
}
void AxisT(TH1* h, const char* xt, const char* yt) {
    h->GetXaxis()->SetTitle(xt); h->GetYaxis()->SetTitle(yt);
    h->GetXaxis()->SetLabelSize(0.05); h->GetYaxis()->SetLabelSize(0.05);
    h->GetXaxis()->SetTitleSize(0.055); h->GetYaxis()->SetTitleSize(0.055);
    h->GetXaxis()->SetMoreLogLabels(false);
}
std::vector<double> LogEdges(double a, double b, int n) {
    std::vector<double> e; for (int i = 0; i <= n; ++i) e.push_back(a * std::pow(b / a, (double)i / n)); return e;
}
std::vector<double> LinEdges(double a, double b, int n) {
    std::vector<double> e; for (int i = 0; i <= n; ++i) e.push_back(a + (b - a) * i / n); return e;
}
int EdgeBin(const std::vector<double>& e, double x) {
    if (!(x >= e.front() && x < e.back())) return -1;
    return (int)(std::upper_bound(e.begin(), e.end(), x) - e.begin()) - 1;
}

// ---------------------------------------------------------------- constants of this study
const double kPeakLo = 0.5333, kPeakHi = 0.5587;   // PEAK window (fixed for all agents)
const double kBumpLo = 0.4576, kBumpHi = 0.5255;   // BUMP window
const double kFileScale = 100.;                    // the generator sample covers 1000 of the 100 000 files
const double kMX17 = 0.017;                        // GeV
const double kBrX17Ill = 1e-5;                     // illustrative BR(eta -> pi+pi- X17, X17 -> e+e-)
const double kBRpi0D = simtruth::kPDG_pi0D, kBRppp0 = simtruth::kPDG_etaPPP0, kBRppg = simtruth::kPDG_etaPPG;
const double kQthr = 1.552e-3;                     // threshold bin of the generator density G (stored precision)

// signal classes of this step
struct CDef { int cls; const char* key; const char* label; int col; int ms; };
const int kNC = 5;
const CDef kC[kNC] = {
    {5, "direct (5)", "#eta#rightarrow#pi^{+}#pi^{-}e^{+}e^{-} (5)", 633, 20},      // kRed+1
    {6, "X17 (6)", "X17 (6)", 635, 29},                                               // kRed+3
    {2, "bump Dalitz (2)", "bump, #pi^{0} Dalitz (2)", 807, 21},                     // kOrange+7
    {3, "bump conv. (3)", "bump, #gamma conversion (3)", 902, 25},                   // kPink+2
    {4, "peak conv. (4)", "peak, #eta#rightarrow#pi^{+}#pi^{-}#gamma conv. (4)", 618, 22}  // kMagenta+2
};
int CIdx(int cls) { for (int i = 0; i < kNC; ++i) if (kC[i].cls == cls) return i; return -1; }

// component of one class: X17 rows weighted 1 per row (they carry W = 0 in the production)
Comp CompOf(int cls) {
    if (cls == 6) {
        Comp c{"X17 (1 per row)", [](const Row& r) { return r.tmpl == 1 && r.cls == 6; }};
        c.wf = [](const Row&) { return 1.; };
        return c;
    }
    return CompSigCls(cls);
}
// weight of a row of class cls as the engine uses it (owned, mode, phys); X17: 1 per owned row
double RowW(const Row& r, int cls, const SimOpt& o) {
    if (!Owned(r, o)) return 0.;
    if (cls == 6) return 1.;
    return Weight(r, o);
}

// ---------------------------------------------------------------- generator sample
struct Gen { std::vector<double> q, m, ce, cp, me, w3; };
bool ReadGen(Gen& g) {
    TStopwatch sw;
    const std::string cmd = std::string("zcat ") + kGenText;
    FILE* p = gSystem->OpenPipe(cmd.c_str(), "r");
    if (!p) { printf("cannot read %s\n", kGenText); return false; }
    double a, b, c, d, e;
    while (fscanf(p, "%lf %lf %lf %lf %lf", &a, &b, &c, &d, &e) == 5) {
        g.q.push_back(a); g.m.push_back(b); g.ce.push_back(c); g.cp.push_back(d); g.me.push_back(e);
    }
    gSystem->ClosePipe(p);
    g.w3.resize(g.q.size());
    for (size_t i = 0; i < g.q.size(); ++i) g.w3[i] = W5(g.q[i], g.m[i], kV3, 0, g.ce[i], g.cp[i], true, 0.);
    printf("ReadGen: %zu decays in %.1f s\n", g.q.size(), sw.RealTime());
    return !g.q.empty();
}

// ---------------------------------------------------------------- analytic production densities (unnormalised)
double KWpi0(double q) {   // Kroll-Wada pi0 -> e+e-gamma, per dq (no form factor; HGeant pi0 Dalitz ASSUMED QED)
    if (q <= 2 * kMe || q >= kMpi0) return 0.;
    const double x = q * q / (kMpi0 * kMpi0), r = kMe * kMe / (q * q);
    return (1. / q) * std::pow(1 - x, 3) * (1 + 2 * r) * std::sqrt(1 - 4 * r);
}
double RhoBump(double M, bool kloe) { return Rho3(M) * (kloe ? WBump(M, kBumpKLOE) : 1.); }     // eta -> pi+pi-pi0
double RhoConv(double M, bool fv) { return RhoG(M) * (fv ? WConv(M, 0) : 1.); }                // eta -> pi+pi-gamma
double RhoX17(double M) {   // eta -> pi+pi- X17: 3-body phase space (ASSUMED; PLUTO's X17 M(pipi) not checked)
    return (M > 2 * kMpi && M < kMetaGen - kMX17) ? PStar(M, kMpi) * K3(kMetaGen, M, kMX17) : 0.;
}
// integral of f over [a, b]; log substitution when logsub (densities ~ 1/q)
double Integ(const std::function<double(double)>& f, double a, double b, bool logsub = false, int n = 2000) {
    if (b <= a) return 0.;
    double s = 0;
    if (logsub && a > 0) {
        const double la = std::log(a), lb = std::log(b);
        for (int i = 0; i < n; ++i) { const double u = la + (lb - la) * (i + 0.5) / n, x = std::exp(u); s += f(x) * x * (lb - la) / n; }
    } else {
        for (int i = 0; i < n; ++i) { const double x = a + (b - a) * (i + 0.5) / n; s += f(x) * (b - a) / n; }
    }
    return s;
}

// ---------------------------------------------------------------- weighted quantiles
struct Quant { double med = 0, q16 = 0, q84 = 0, q025 = 0, q975 = 0, q05 = 0, q95 = 0, sw = 0, sw2 = 0, mean = 0; };
Quant Quantiles(std::vector<std::pair<double, double>> v) {   // (value, weight)
    Quant q;
    if (v.empty()) return q;
    std::sort(v.begin(), v.end());
    double sw = 0, sw2 = 0, sx = 0;
    for (auto& p : v) { sw += p.second; sw2 += p.second * p.second; sx += p.first * p.second; }
    q.sw = sw; q.sw2 = sw2; q.mean = sx / sw;
    auto at = [&](double f) {
        double c = 0; const double t = f * sw;
        for (auto& p : v) { c += p.second; if (c >= t) return p.first; }
        return v.back().first;
    };
    q.med = at(0.5); q.q16 = at(0.15866); q.q84 = at(0.84134); q.q025 = at(0.025); q.q975 = at(0.975); q.q05 = at(0.05); q.q95 = at(0.95);
    return q;
}
double NeffOf(const Quant& q) { return q.sw2 > 0 ? q.sw * q.sw / q.sw2 : 0.; }

// ---------------------------------------------------------------- Gaussian core fit (iterated mean +- n sigma)
struct GFit { double mu = 0, emu = 0, sig = 0, esig = 0, chi2 = 0; int ndf = 0; bool ok = false; };
GFit FitCore(TH1D* h, double mu0, double s0, double nsig = 2., int nit = 8) {
    GFit g; double mu = mu0, s = std::max(s0, 1e-9);
    for (int it = 0; it < nit; ++it) {
        TF1 f(Nm("fg").c_str(), "gaus", mu - nsig * s, mu + nsig * s);
        const int b0 = h->FindBin(mu);
        f.SetParameters(std::max(h->GetBinContent(b0), 1e-30), mu, s);
        const int st = h->Fit(&f, "QRN0");
        if (st != 0 && it == nit - 1) { g.ok = false; return g; }
        const double m1 = f.GetParameter(1), s1 = std::fabs(f.GetParameter(2));
        if (!(s1 > 0) || !(std::fabs(m1 - mu) < 5 * s)) break;
        mu = m1; s = s1;
        g.mu = mu; g.emu = f.GetParError(1); g.sig = s; g.esig = f.GetParError(2); g.chi2 = f.GetChisquare(); g.ndf = f.GetNDF();
        g.ok = true;
    }
    return g;
}

// ---------------------------------------------------------------- data held across the parts
struct Ctx {
    Gen gen;
    double Ngen = 0, NgenW3 = 0;   // sample decays, sum of V3 weights over the sample
    double Ndir = 0, Neta = 0;     // full production: direct decays, etas
    double Wprod[7] = {0};         // produced weight units per class (denominator of eps, same units as the accepted sums)
    double sigProd[7] = {0};       // production cross section per class [nb] (X17 for BR = kBrX17Ill)
    double relN = 0;               // relative statistical error of the production scale
    SimData sig;                   // signal rows (template 1, classes 2-6) of std, lep, eta e+e- trees
    bool loaded = false;
    double x17lo = 0, x17hi = 0;        // shortest reconstructed M(ee) interval with 90 % of the X17 rows [MeV]
    std::vector<std::string> summary;   // key lines for the log end
};
Ctx& C() { static Ctx c; return c; }

SimOpt OptMain() { SimOpt o; o.mode = 0; o.db = false; o.sets = kSetsAll; o.hyp = kUS; o.phys = kPhysMain; return o; }
SimOpt OptOff() { SimOpt o = OptMain(); o.phys = 0; return o; }

void SetupProduction() {
    Ctx& c = C();
    if (c.Ngen > 0) return;
    ReadGen(c.gen);
    c.Ngen = (double)c.gen.q.size();
    for (double w : c.gen.w3) c.NgenW3 += w;
    c.Ndir = kFileScale * c.Ngen;
    c.Neta = 2. * c.Ndir;                         // q = 1/2: every eta forced into pi+pi-e+e- or X17 pi+pi-
    c.relN = 1. / std::sqrt(c.Ngen);
    // produced weight units: class 5 rows carry W = 2 BR(ppee) (x V3 factor); X17 rows 1 per row (wf); std rows W = PDG/HGeant
    // factors, so their accepted sum counts physical (PDG) decays
    c.Wprod[5] = c.Ndir * 2. * kBRppee;
    c.Wprod[6] = c.Ndir;                          // N(X17 decays) = N_eta / 2
    c.Wprod[2] = c.Neta * kBRppp0 * kBRpi0D;
    c.Wprod[3] = c.Neta * kBRppp0;                // per eta -> pi+pi-pi0 decay (any pi0 decay; conversion is a detector process)
    c.Wprod[4] = c.Neta * kBRppg;                 // per eta -> pi+pi-gamma decay
    c.sigProd[5] = c.Neta * kBRppee / kLsim;
    c.sigProd[6] = c.Neta * kBrX17Ill / kLsim;
    c.sigProd[2] = c.Neta * kBRppp0 * kBRpi0D / kLsim;
    c.sigProd[3] = c.Neta * kBRppp0 / kLsim;
    c.sigProd[4] = c.Neta * kBRppg / kLsim;
}

bool LoadSim() {
    Ctx& c = C();
    if (c.loaded) return true;
    TStopwatch sw;
    if (!Mdl().ok && !Init(kRwGenFile)) return false;
    Install();
    LoadAngles(kRwAngFile);
    SimData d;
    if (!Load(d, kCompact, kMaskUS)) return false;
    // keep only the signal rows (template 1, classes 2-6), in file order (event groups stay adjacent)
    long long nsig[kNTreeAll] = {0};
    for (int k = 0; k < kNTreeAll; ++k) {
        for (const Row& r : d.t[k])
            if (r.tmpl == 1 && r.cls >= 2 && r.cls <= 6) { c.sig.t[k].push_back(r); ++nsig[k]; }
        std::vector<Row>().swap(d.t[k]);
    }
    c.sig.loaded = d.loaded;
    c.loaded = true;
    printf("LoadSim: signal rows std_us %lld, lep_us %lld, eta_us %lld (%.0f s)\n", nsig[kStdUS], nsig[kLepUS], nsig[kEtaUS], sw.RealTime());
    return true;
}

// rows of one class passing sel: (x, y, w) with x, y from functions
std::vector<std::array<double, 3>> Collect(int cls, const Sel& sel, const SimOpt& o, const std::function<double(const Row&)>& fx,
                                           const std::function<double(const Row&)>& fy = std::function<double(const Row&)>()) {
    std::vector<std::array<double, 3>> out;
    const SimData& d = C().sig;
    for (int set = 1; set <= 3; ++set) {
        if (!(o.sets & SetBit(set))) continue;
        for (const Row& r : d.tree(set, o.hyp)) {
            if (!(r.tmpl == 1 && r.cls == cls)) continue;
            if (!sel(r)) continue;
            const double w = RowW(r, cls, o);
            if (w == 0) continue;
            out.push_back({fx(r), fy ? fy(r) : 0., w});
        }
    }
    return out;
}
// binning of SimEngine over arbitrary edges: x = bin index + 0.5
Binning IdxBinning(const std::vector<double>& e, const std::function<double(const Row&)>& val) {
    Binning b; b.nbin = (int)e.size() - 1; b.lo = 0; b.hi = b.nbin;
    b.x = [e, val](const Row& r) { const int k = EdgeBin(e, val(r)); return k < 0 ? -1. : k + 0.5; };
    return b;
}
// true e+e- pair momentum [MeV/c] from the true momenta and the true opening angle
double PeeTrue(const Row& r) {
    if (!(r.oa_true >= 0)) return -1.;
    const double a = r.ep_sim_p, b = r.em_sim_p, c = std::cos(r.oa_true * kD2R);
    return std::sqrt(std::max(0., a * a + b * b + 2 * a * b * c));
}

// graph of y per bin with x errors = bin extent (x centre: geometric for log bins)
TGraphAsymmErrors* BinGraph(const std::vector<double>& e, const std::vector<double>& y, const std::vector<double>& ey, bool logx,
                            double xscale = 1., bool skipZero = true) {
    TGraphAsymmErrors* g = new TGraphAsymmErrors();
    for (size_t i = 0; i + 1 < e.size(); ++i) {
        if (skipZero && !(y[i] > 0)) continue;
        const double xc = logx ? std::sqrt(e[i] * e[i + 1]) : 0.5 * (e[i] + e[i + 1]);
        const int n = g->GetN();
        g->SetPoint(n, xc * xscale, y[i]);
        g->SetPointError(n, (xc - e[i]) * xscale, (e[i + 1] - xc) * xscale, ey[i], ey[i]);
    }
    return g;
}

// =====================================================================================================
// Part 1: generator level
// =====================================================================================================
void PartGen() {
    TStopwatch sw;
    Style();
    if (!Mdl().ok && !Init(kRwGenFile)) return;
    SetupProduction();
    Ctx& c = C();
    const Gen& g = c.gen;
    OpenTab("s5_acc_generator.txt");
    P("# s5_acc_generator: the eta signals at generator level (s5_acc.C part 1). SIMULATION ONLY.\n");
    P("# Already shown elsewhere (cited, not redone): generator M(ee), M(pipi) and decay angles of the direct decay against the\n");
    P("#   target models V0-V2 (figs/rw_weights_generator, tables/rw_weights_generator.txt, rw_weights_model.txt); the V3-reweighted\n");
    P("#   generator sample against T in M(ee), M(pipi) and both angles (figs/rw_weights_closure, tables/rw_weights_closure.txt);\n");
    P("#   the bump weight wB(M(pipi)), <X^2>_Y check, conversion weight wC (figs/rw_weights_bump, tables/rw_weights_bump.txt).\n");
    P("# Added here: absolute production per channel, produced M(ee) / M(pipi) spectra of all signal channels in one scale,\n");
    P("#   M(ee) fractions relevant for OA and X17, and the generator-level direct decays under the X17 line.\n\n");
    P("## 1. production scale\n");
    P("generator sample (direct eta -> pi+pi-e+e-, 1000 of 100 000 files): N = %.0f decays; sum of V3 weights %.6f (mean %.6f)\n", c.Ngen, c.NgenW3,
      c.NgenW3 / c.Ngen);
    P("full production: N_direct = 100 x N = %.4e; N_eta = 2 N_direct = %.4e (eta set: every eta forced, q = 1/2 into pi+pi-e+e- and X17 pi+pi-)\n",
      c.Ndir, c.Neta);
    P("relative statistical error of the scale 1/sqrt(N) = %.2e (+ the ASSUMPTION that the 1000 files are representative, UNVERIFIED)\n", c.relN);
    P("sigma(p p -> eta X) in SMASH = N_eta / L_sim = %.4f mb (L_sim = %.2f nb^-1)\n", c.Neta / kLsim * 1e-6, kLsim);
    P("PDG BR used: pi+pi-e+e- %.3e, pi+pi-pi0 %.4f, pi0 -> e+e-gamma %.5f, pi+pi-gamma %.4f; X17 illustrative %.0e\n", kBRppee, kBRppp0, kBRpi0D,
      kBRppg, kBrX17Ill);
    P("\nchannel                                      produced decays (10^10 events)    sigma_prod [nb]\n");
    const char* pn[7] = {"", "", "class 2 source: eta->pi+pi-pi0, pi0->e+e-gamma", "class 3 source: eta->pi+pi-pi0 (all pi0)",
                         "class 4 source: eta->pi+pi-gamma", "class 5: eta->pi+pi-e+e-", "class 6 source: X17 (BR 1e-5)"};
    for (int k = 2; k <= 6; ++k) {
        const double nprod = c.sigProd[k] * kLsim;
        P("%-46s %14.4e                  %12.4f\n", pn[k], nprod, c.sigProd[k]);
    }
    // ---- histograms (per GeV) of the produced spectra in nb/(GeV/c2)
    // threshold bin 2 m_e - 1.552 MeV (generator stored precision: 4.9 % of the values below 2 m_e), then log bins
    std::vector<double> qe = {2 * kMe};
    for (double x : LogEdges(kQthr, QMax(), 56)) qe.push_back(x);
    const std::vector<double> me = LinEdges(2 * kMpi, kMetaGen, 54);
    TH1D* hqG = H1v("hqG", qe); TH1D* hqV = H1v("hqV", qe);
    TH1D* hmG = H1v("hmG", me); TH1D* hmV = H1v("hmV", me);
    for (size_t i = 0; i < g.q.size(); ++i) {
        const double q = std::max(g.q[i], 2 * kMe * 1.0000001);
        hqG->Fill(q); hqV->Fill(q, g.w3[i]); hmG->Fill(g.m[i]); hmV->Fill(g.m[i], g.w3[i]);
    }
    const double s5 = c.sigProd[5];
    for (TH1D* h : {hqG, hmG}) h->Scale(s5 / c.Ngen, "width");
    for (TH1D* h : {hqV, hmV}) h->Scale(s5 / c.NgenW3, "width");
    // bump pi0 Dalitz pair (class 2 source) in M(ee): Kroll-Wada; X17 a line at 17 MeV
    const double zKW = Integ(KWpi0, 2 * kMe * (1 + 1e-9), kMpi0, true, 20000);
    TH1D* hqB = H1v("hqB", qe);
    for (int b = 1; b <= hqB->GetNbinsX(); ++b)
        hqB->SetBinContent(b, c.sigProd[2] * Integ(KWpi0, qe[b - 1], qe[b], true, 200) / zKW / hqB->GetBinWidth(b));
    // M(pipi): bump source (all pi0 decays: class 3; class 2 = x BR(pi0 Dalitz)), conversion source, X17
    TH1D* hmB0 = H1v("hmB0", me); TH1D* hmB1 = H1v("hmB1", me); TH1D* hmC0 = H1v("hmC0", me); TH1D* hmC1 = H1v("hmC1", me);
    TH1D* hmX = H1v("hmX", me);
    const double zB0 = Integ([](double M) { return RhoBump(M, false); }, 2 * kMpi, kMetaHG - kMpi0, false, 20000);
    const double zB1 = Integ([](double M) { return RhoBump(M, true); }, 2 * kMpi, kMetaHG - kMpi0, false, 20000);
    const double zC0 = Integ([](double M) { return RhoConv(M, false); }, 2 * kMpi, kMetaHG, false, 20000);
    const double zC1 = Integ([](double M) { return RhoConv(M, true); }, 2 * kMpi, kMetaHG, false, 20000);
    const double zX = Integ(RhoX17, 2 * kMpi, kMetaGen - kMX17, false, 20000);
    for (int b = 1; b <= hmB0->GetNbinsX(); ++b) {
        const double a = me[b - 1], bb = me[b], w = bb - a;
        hmB0->SetBinContent(b, c.sigProd[3] * Integ([](double M) { return RhoBump(M, false); }, a, bb) / zB0 / w);
        hmB1->SetBinContent(b, c.sigProd[3] * Integ([](double M) { return RhoBump(M, true); }, a, bb) / zB1 / w);
        hmC0->SetBinContent(b, c.sigProd[4] * Integ([](double M) { return RhoConv(M, false); }, a, bb) / zC0 / w);
        hmC1->SetBinContent(b, c.sigProd[4] * Integ([](double M) { return RhoConv(M, true); }, a, bb) / zC1 / w);
        hmX->SetBinContent(b, c.sigProd[6] * Integ(RhoX17, a, bb) / zX / w);
    }
    P("\n## 2. normalisation checks of the analytic densities (should be 1): <wB>_PS(KLOE) = %.5f, <wC>_PS = %.5f\n", zB1 / zB0, zC1 / zC0);
    // ---- M(ee) fractions
    auto fracAbove = [&](double q0, bool v3) {
        double s = 0, t = 0;
        for (size_t i = 0; i < g.q.size(); ++i) { const double w = v3 ? g.w3[i] : 1.; t += w; if (g.q[i] > q0) s += w; }
        return s / t;
    };
    auto fracKW = [&](double q0) { return Integ(KWpi0, std::max(q0, 2 * kMe * (1 + 1e-9)), kMpi0, true, 20000) / zKW; };
    P("\n## 3. produced M(ee): fraction of decays above a threshold (generator = PLUTO sample unweighted; V3 = physical model;\n");
    P("#   pi0 Dalitz = Kroll-Wada of the class-2 source)\n");
    P("M(ee) >      direct gen    direct V3    pi0 Dalitz (KW)\n");
    for (double q0 : {0.002, 0.005, 0.010, 0.017, 0.020, 0.050, 0.100}) P("%5.0f MeV    %8.4f     %8.4f     %8.4f\n", q0 * 1e3, fracAbove(q0, false), fracAbove(q0, true), fracKW(q0));
    P("mean M(pipi) of the produced spectra [GeV/c2]: direct gen %.4f, direct V3 %.4f, bump PS %.4f, bump KLOE %.4f, conv PS %.4f, conv wC %.4f, X17 PS %.4f\n",
      hmG->GetMean(), hmV->GetMean(), hmB0->GetMean(), hmB1->GetMean(), hmC0->GetMean(), hmC1->GetMean(), hmX->GetMean());
    // ---- direct decays under the X17 line at generator level
    P("\n## 4. generator-level direct decays in M(ee) windows around 17 MeV (V3) against X17 at BR = 1e-5 (%.4f nb, all of it at 17.000 MeV)\n",
      c.sigProd[6]);
    P("half width [MeV]   sigma_direct(|M(ee) - 17| < d) [nb]   (gen. unweighted)   ratio direct/X17(1e-5)\n");
    std::vector<double> hw = {0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0};
    std::vector<double> sdV, sdG;
    for (double d : hw) {
        double sv = 0, sg = 0;
        for (size_t i = 0; i < g.q.size(); ++i) if (std::fabs(g.q[i] - kMX17) < d * 1e-3) { sv += g.w3[i]; sg += 1; }
        sv *= s5 / c.NgenW3; sg *= s5 / c.Ngen;
        sdV.push_back(sv); sdG.push_back(sg);
        P("%8.2f             %10.4f                            %10.4f            %8.2f\n", d, sv, sg, sv / c.sigProd[6]);
    }
    CloseTab();

    // ---- figure
    TCanvas* cv = new TCanvas(Nm("cgen").c_str(), "", 1700, 1400);
    cv->Divide(2, 2, 0.003, 0.003);
    // (a) M(ee)
    cv->cd(1); PadM(0.16, 0.03, 0.14, 0.09); gPad->SetLogx(); gPad->SetLogy();
    TH1D* fa = H1("fa", 10, 2 * kMe, QMax());
    AxisT(fa, "M(e^{+}e^{-})_{true} [GeV/c^{2}]", "d#sigma_{prod}/dM(e^{+}e^{-}) [nb/(GeV/c^{2})]");
    fa->SetMinimum(2e1); fa->SetMaximum(3e9); fa->Draw("axis");
    Sty(hqB, 807, 3, 1); hqB->Draw("hist same");
    Sty(hqG, 1, 3, 2); hqG->Draw("hist same");
    Sty(hqV, 633, 4, 1); hqV->Draw("hist same");
    // X17: a line of sigma 0.107 nb (BR 1e-5) at 17 MeV, drawn as a 0.2 MeV wide bin for visibility
    {
        const double dX = 0.0002;
        TLine lx; lx.SetLineColor(635); lx.SetLineWidth(4); lx.SetLineStyle(2);
        lx.DrawLine(kMX17, 2e1, kMX17, c.sigProd[6] / dX);
        TMarker mk(kMX17, c.sigProd[6] / dX, 29); mk.SetMarkerColor(635); mk.SetMarkerSize(2.5); mk.DrawClone();
    }
    TLegend* la = Leg(0.18, 0.62, 0.97, 0.89, 0.042);
    la->AddEntry(hqV, Form("#eta#rightarrow#pi^{+}#pi^{-}e^{+}e^{-}, model V3 (%.0f nb)", c.sigProd[5]), "l");
    la->AddEntry(hqG, "#eta#rightarrow#pi^{+}#pi^{-}e^{+}e^{-}, PLUTO generator (not QED)", "l");
    la->AddEntry(hqB, Form("#eta#rightarrow#pi^{+}#pi^{-}#pi^{0}, #pi^{0}#rightarrowe^{+}e^{-}#gamma (KW, %.0f nb)", c.sigProd[2]), "l");
    TH1D* dumX = H1("dumX", 1, 0, 1); Sty(dumX, 635, 4, 2, 29, 2.0);
    la->AddEntry(dumX, Form("X17 (BR 10^{-5}, %.1f nb), per 0.2 MeV/c^{2}", c.sigProd[6]), "lp");
    la->Draw();
    Txt(0.16, 0.93, "(a) produced M(e^{+}e^{-}), generator level", 0.055);
    // (b) M(pipi)
    cv->cd(2); PadM(0.16, 0.03, 0.14, 0.09); gPad->SetLogy();
    TH1D* fb = H1("fb", 10, 0.27, 0.56);
    AxisT(fb, "M(#pi^{+}#pi^{-})_{true} [GeV/c^{2}]", "d#sigma_{prod}/dM(#pi^{+}#pi^{-}) [nb/(GeV/c^{2})]");
    fb->SetMinimum(1); fb->SetMaximum(3e12); fb->Draw("axis");
    Sty(hmB0, 807, 3, 3); Sty(hmB1, 807, 4, 1); Sty(hmC0, 618, 3, 3); Sty(hmC1, 618, 4, 1);
    Sty(hmG, 1, 3, 2); Sty(hmV, 633, 4, 1); Sty(hmX, 635, 3, 2);
    for (TH1D* h : {hmB0, hmB1, hmC0, hmC1, hmG, hmV, hmX}) h->Draw("hist same");
    TLegend* lb = Leg(0.18, 0.60, 0.97, 0.90, 0.040);
    lb->AddEntry(hmB1, "#eta#rightarrow#pi^{+}#pi^{-}#pi^{0} (all #pi^{0}), KLOE-2 Dalitz plot", "l");
    lb->AddEntry(hmB0, "  the same, phase space (HGeant)", "l");
    lb->AddEntry(hmC1, "#eta#rightarrow#pi^{+}#pi^{-}#gamma, |F_{V}|^{2}p*^{2}k^{2}", "l");
    lb->AddEntry(hmC0, "  the same, phase space (HGeant)", "l");
    lb->AddEntry(hmV, "#eta#rightarrow#pi^{+}#pi^{-}e^{+}e^{-}, V3 (generator dashed black)", "l");
    lb->AddEntry(hmX, "X17 (BR 10^{-5}), 3-body phase space ASSUMED", "l");
    lb->Draw();
    Txt(0.16, 0.93, "(b) produced M(#pi^{+}#pi^{-}), generator level", 0.055);
    // (c) fraction above M(ee)
    cv->cd(3); PadM(0.16, 0.03, 0.14, 0.09); gPad->SetLogx();
    const std::vector<double> qq = LogEdges(2 * kMe * 1.001, QMax(), 120);
    TGraph* gG = new TGraph(); TGraph* gV = new TGraph(); TGraph* gK = new TGraph();
    {
        // cumulative from sorted sample
        std::vector<std::pair<double, double>> s;
        s.reserve(g.q.size());
        for (size_t i = 0; i < g.q.size(); ++i) s.push_back(std::make_pair(g.q[i], g.w3[i]));
        std::sort(s.begin(), s.end());
        size_t j = 0; double cw = 0, cn = 0;
        for (double q0 : qq) {
            while (j < s.size() && s[j].first <= q0) { cw += s[j].second; cn += 1; ++j; }
            gG->SetPoint(gG->GetN(), q0, 1 - cn / c.Ngen);
            gV->SetPoint(gV->GetN(), q0, 1 - cw / c.NgenW3);
            if (q0 < kMpi0) gK->SetPoint(gK->GetN(), q0, fracKW(q0));
        }
    }
    TH1D* fc = H1("fc", 10, 2 * kMe, QMax());
    AxisT(fc, "M(e^{+}e^{-})_{true} threshold [GeV/c^{2}]", "fraction of decays above");
    fc->SetMinimum(0); fc->SetMaximum(1.05); fc->Draw("axis");
    Sty(gG, 1, 1, 0.1, 4, 2); Sty(gV, 633, 1, 0.1, 4, 1); Sty(gK, 807, 1, 0.1, 4, 1);
    gG->Draw("l same"); gV->Draw("l same"); gK->Draw("l same");
    TLine l17; l17.SetLineColor(635); l17.SetLineStyle(2); l17.SetLineWidth(3); l17.DrawLine(kMX17, 0, kMX17, 1.05);
    TLegend* lc = Leg(0.36, 0.62, 0.97, 0.89, 0.042);
    lc->AddEntry(gV, Form("direct, V3 (> 20 MeV: %.3f)", fracAbove(0.020, true)), "l");
    lc->AddEntry(gG, Form("direct, generator (> 20 MeV: %.3f)", fracAbove(0.020, false)), "l");
    lc->AddEntry(gK, Form("#pi^{0} Dalitz, KW (> 20 MeV: %.3f)", fracKW(0.020)), "l");
    lc->AddEntry(&l17, "X17 mass, 17 MeV/c^{2}", "l");
    lc->Draw();
    Txt(0.16, 0.93, "(c) produced decays above a given M(e^{+}e^{-})", 0.055);
    // (d) direct decays under the X17 line
    cv->cd(4); PadM(0.16, 0.03, 0.14, 0.09);
    TGraph* gdV = new TGraph(); TGraph* gdG = new TGraph();
    for (size_t i = 0; i < hw.size(); ++i) { gdV->SetPoint(i, hw[i], sdV[i]); gdG->SetPoint(i, hw[i], sdG[i]); }
    TH1D* fd = H1("fd", 10, 0, 6.5);
    AxisT(fd, "half width d of the M(e^{+}e^{-}) window around 17 MeV/c^{2} [MeV/c^{2}]", "#sigma_{prod} in |M_{ee} - 17| < d [nb]");
    fd->GetXaxis()->SetTitleSize(0.048);
    fd->SetMinimum(0); fd->SetMaximum(1.25 * std::max(sdV.back(), c.sigProd[6])); fd->Draw("axis");
    Sty(gdV, 633, 20, 1.6, 3); Sty(gdG, 1, 24, 1.6, 3, 2);
    gdV->Draw("pl same"); gdG->Draw("pl same");
    TLine lX; lX.SetLineColor(635); lX.SetLineWidth(4); lX.SetLineStyle(2); lX.DrawLine(0, c.sigProd[6], 6.5, c.sigProd[6]);
    TLegend* ld = Leg(0.19, 0.66, 0.95, 0.88, 0.045);
    ld->AddEntry(gdV, "#eta#rightarrow#pi^{+}#pi^{-}e^{+}e^{-}, V3", "pl");
    ld->AddEntry(gdG, "#eta#rightarrow#pi^{+}#pi^{-}e^{+}e^{-}, generator", "pl");
    ld->AddEntry(&lX, Form("X17, BR = 10^{-5} (%.1f nb)", c.sigProd[6]), "l");
    ld->Draw();
    Txt(0.16, 0.93, "(d) direct decays near M(e^{+}e^{-}) = 17 MeV/c^{2}", 0.055);
    cv->cd(0);
    SaveFig(cv, "s5_acc_gen");
    c.summary.push_back(Form("sigma_prod: eta %.4f mb; direct %.2f nb; X17(1e-5) %.3f nb; class-2 source %.1f nb", c.Neta / kLsim * 1e-6, c.sigProd[5],
                             c.sigProd[6], c.sigProd[2]));
    printf("PartGen done in %.0f s\n", sw.RealTime());
}

// =====================================================================================================
// Part 2: acceptance x efficiency of "pre" per produced decay
// =====================================================================================================
struct EffBins {
    std::vector<double> e, acc, accV, prod, prodE, eps, epsE;   // eps = acc / prod, epsE from groups (+) production stat
    int n() const { return (int)e.size() - 1; }
};
// accepted sums per bin of val(row) (SimEngine FillMany, group variance)
void FillAcc(EffBins& b, int cls, const Sel& sel, const SimOpt& o, const std::function<double(const Row&)>& val) {
    const WSpec sp = Fill(C().sig, sel, CompOf(cls), o, IdxBinning(b.e, val));
    b.acc.assign(b.n(), 0.); b.accV.assign(b.n(), 0.);
    for (int i = 0; i < b.n(); ++i) { b.acc[i] = sp.s[i]; b.accV[i] = sp.v[i]; }
}
void Divide(EffBins& b) {
    b.eps.assign(b.n(), 0.); b.epsE.assign(b.n(), 0.);
    for (int i = 0; i < b.n(); ++i) {
        if (!(b.prod[i] > 0)) continue;
        b.eps[i] = b.acc[i] / b.prod[i];
        const double r1 = b.acc[i] > 0 ? std::sqrt(b.accV[i]) / b.acc[i] : 0., r2 = b.prodE[i] / b.prod[i];
        b.epsE[i] = b.eps[i] * std::sqrt(r1 * r1 + r2 * r2);
        if (!(b.acc[i] > 0)) b.epsE[i] = std::sqrt(b.accV[i]) / b.prod[i];
    }
}
// produced class-5 decays per bin from the generator sample (weights w5 V3 or 1), in produced weight units
void ProdFromGen(EffBins& b, bool v3, bool useQ) {
    Ctx& c = C();
    b.prod.assign(b.n(), 0.); b.prodE.assign(b.n(), 0.);
    std::vector<double> s2(b.n(), 0.);
    const double norm = v3 ? c.Ngen / c.NgenW3 : 1.;
    for (size_t i = 0; i < c.gen.q.size(); ++i) {
        double x = useQ ? std::max(c.gen.q[i], 2 * kMe * 1.0000001) : c.gen.m[i];
        const int k = EdgeBin(b.e, x);
        if (k < 0) continue;
        const double w = (v3 ? c.gen.w3[i] : 1.) * norm;
        b.prod[k] += w; s2[k] += w * w;
    }
    const double f = c.Wprod[5] / c.Ngen;
    for (int k = 0; k < b.n(); ++k) { b.prodE[k] = std::sqrt(s2[k]) * f; b.prod[k] *= f; }
}
void ProdAnalytic(EffBins& b, double total, const std::function<double(double)>& rho, double lo, double hi, bool logsub) {
    b.prod.assign(b.n(), 0.); b.prodE.assign(b.n(), 0.);
    const double z = Integ(rho, lo, hi, logsub, 20000);
    for (int k = 0; k < b.n(); ++k) {
        const double a = std::max(b.e[k], lo), bb = std::min(b.e[k + 1], hi);
        b.prod[k] = bb > a ? total * Integ(rho, a, bb, logsub, 400) / z : 0.;
    }
}

void PartAcc() {
    TStopwatch sw;
    Style();
    SetupProduction();
    if (!LoadSim()) return;
    Ctx& c = C();
    const Sel pre = SelPreset(kPre);
    const SimOpt oM = OptMain(), o0 = OptOff();
    OpenTab("s5_acc_efficiency.txt");
    P("# s5_acc_efficiency: acceptance x efficiency of 'pre' per PRODUCED decay (s5_acc.C part 2). SIMULATION ONLY, no CB.\n");
    P("# 'pre' = isBest >= 0, vz > -500, PT3 bit, vz > -160 mm, corrected momenta, cor_ok; rows owned in mode A, std + lep + eta (sets = 7),\n");
    P("#   w (PDG), physics weights V3 main (simphys::kPhysMain = %d); 'off' = without the physics weights.\n", kPhysMain);
    P("# eps = (accepted sum of row weights) / (produced decays in the same weight units): ROWS per produced decay (a decay can give\n");
    P("#   more than one row; the rows/decay ratio is listed). Produced decays:\n");
    P("#   class 5: N_direct = 100 x %.0f (generator sample of 1000 of the 100 000 files), each row W = 2 BR(ppee) (x V3 factor)\n", c.Ngen);
    P("#   class 6: N_X17 = N_eta / 2 = N_direct (q = 1/2), one per row (production W = 0)\n");
    P("#   class 2: N_eta BR(pi+pi-pi0) BR(pi0 Dalitz), class 3: N_eta BR(pi+pi-pi0) (any pi0 decay), class 4: N_eta BR(pi+pi-gamma) (PDG;\n");
    P("#     std rows carry W = PDG/HGeant factors, so their accepted sum counts PDG decays); N_eta = 2 N_direct = %.4e\n", c.Neta);
    P("# errors: accepted sums from event groups (SimEngine) (+) production scale %.1e (relative, statistical)\n\n", c.relN);

    // ---- totals
    P("## 1. totals ('pre', all M4 and in the M4 windows PEAK %.4f-%.4f, BUMP %.4f-%.4f GeV/c2)\n", kPeakLo, kPeakHi, kBumpLo, kBumpHi);
    P("class               rows    groups  rows/grp   N_eff       sigma_acc [nb]        eps (V3) [%%]           eps (off) [%%]     eps in PEAK [%%]    eps in BUMP [%%]\n");
    for (int i = 0; i < kNC; ++i) {
        const int cls = kC[i].cls;
        const std::vector<Comp> cm = {CompOf(cls)};
        const WSum sM = SumMany(c.sig, pre, cm, oM)[0];
        const WSum s0 = SumMany(c.sig, pre, cm, o0)[0];
        const WSum sP = SumMany(c.sig, [pre](const Row& r) { return pre(r) && r.m4 >= kPeakLo && r.m4 < kPeakHi; }, cm, oM)[0];
        const WSum sB = SumMany(c.sig, [pre](const Row& r) { return pre(r) && r.m4 >= kBumpLo && r.m4 < kBumpHi; }, cm, oM)[0];
        // groups with >= 1 row
        long long ng = 0;
        for (int set = 1; set <= 3; ++set) {
            int cur = -1; bool had = false;
            for (const Row& r : c.sig.tree(set, kUS)) {
                if (r.grp != cur) { if (had) ++ng; had = false; cur = r.grp; }
                if (r.tmpl == 1 && r.cls == cls && Owned(r, oM) && pre(r)) had = true;
            }
            if (had) ++ng;
        }
        auto ep = [&](const WSum& s, double& e) { const double v = s.S / c.Wprod[cls]; const double rs = s.S > 0 ? std::sqrt(s.V) / s.S : 0.;
                                                  e = v * std::sqrt(rs * rs + c.relN * c.relN); return v; };
        double eM, e0, eP, eB;
        const double vM = ep(sM, eM), v0 = ep(s0, e0), vP = ep(sP, eP), vB = ep(sB, eB);
        const double sig = cls == 6 ? sM.S * 2 * kBrX17Ill / kLsim : sM.S / kLsim;
        const double sigE = cls == 6 ? std::sqrt(sM.V) * 2 * kBrX17Ill / kLsim : std::sqrt(sM.V) / kLsim;
        P("%-16s %8.0f %8lld  %7.4f  %9.0f   %9.5f +- %7.5f   %8.5f +- %7.5f   %8.5f +- %7.5f   %8.5f +- %7.5f   %8.5f +- %7.5f\n", kC[i].key, sM.N, ng,
          ng > 0 ? sM.N / ng : 0., sM.Neff(), sig, sigE, 100 * vM, 100 * eM, 100 * v0, 100 * e0, 100 * vP, 100 * eP, 100 * vB, 100 * eB);
        c.summary.push_back(Form("eps(pre) %-16s %.4f +- %.4f %% (off %.4f), N_eff %.0f; in PEAK %.4f %%, BUMP %.4f %%", kC[i].key, 100 * vM, 100 * eM,
                                 100 * v0, sM.Neff(), 100 * vP, 100 * vB));
    }
    P("# sigma_acc of class 6 for BR = 1e-5 (W = 2 BR per row); the ratio X17 / direct (generator, off) of eps is the 'factor 1.379' of\n");
    P("#   docs/SYMULACJA_ETA_GENERATOR.md 4.3. lep set: %zu signal rows (expected 0: the lep set forces every eta into e+e-gamma).\n", c.sig.t[kLepUS].size());

    // ---- binned
    struct VarDef { const char* name; const char* unit; std::function<double(const Row&)> f; };
    auto fMee = [](const Row& r) { return (double)r.mee_true; };
    auto fMpp = [](const Row& r) { return (double)r.mpipi_true; };
    auto fOA = [](const Row& r) { return (double)r.oa_true; };
    auto fPee = [](const Row& r) { return PeeTrue(r); };
    // M(ee) edges: threshold bin 2 m_e - 1.552 MeV (generator stored precision), then log bins
    std::vector<double> eQ5 = {2 * kMe};
    for (double x : LogEdges(kQthr, QMax(), 20)) eQ5.push_back(x);
    const std::vector<double> eQ2 = {2 * kMe, kQthr, 0.0025, 0.004, 0.0065, 0.010, 0.016, 0.025, 0.040, 0.065, kMpi0};
    const std::vector<double> eM5 = LinEdges(2 * kMpi, kMetaGen, 27);
    const std::vector<double> eM23 = LinEdges(2 * kMpi, kMetaHG - kMpi0, 7);
    const std::vector<double> eM4 = LinEdges(2 * kMpi, kMetaHG, 6);
    const std::vector<double> eM6 = LinEdges(2 * kMpi, kMetaGen - kMX17, 13);

    // class 5 against true M(ee) and M(pipi), V3 and off
    EffBins q5V, q5G, m5V, m5G;
    q5V.e = q5G.e = eQ5; m5V.e = m5G.e = eM5;
    FillAcc(q5V, 5, pre, oM, fMee); ProdFromGen(q5V, true, true); Divide(q5V);
    FillAcc(q5G, 5, pre, o0, fMee); ProdFromGen(q5G, false, true); Divide(q5G);
    FillAcc(m5V, 5, pre, oM, fMpp); ProdFromGen(m5V, true, false); Divide(m5V);
    FillAcc(m5G, 5, pre, o0, fMpp); ProdFromGen(m5G, false, false); Divide(m5G);
    // class 5 restricted to true M(ee) 14-20 MeV (comparison with X17 in M(pipi))
    EffBins m5X; m5X.e = eM6;
    {
        const Sel s1720 = [pre](const Row& r) { return pre(r) && r.mee_true >= 0.014f && r.mee_true < 0.020f; };
        FillAcc(m5X, 5, s1720, oM, fMpp);
        // produced: generator decays with 14 <= q < 20 MeV per M(pipi) bin (V3)
        m5X.prod.assign(m5X.n(), 0.); m5X.prodE.assign(m5X.n(), 0.);
        std::vector<double> s2(m5X.n(), 0.);
        const double norm = c.Ngen / c.NgenW3, f = c.Wprod[5] / c.Ngen;
        for (size_t i = 0; i < c.gen.q.size(); ++i) {
            if (!(c.gen.q[i] >= 0.014 && c.gen.q[i] < 0.020)) continue;
            const int k = EdgeBin(m5X.e, c.gen.m[i]); if (k < 0) continue;
            const double w = c.gen.w3[i] * norm; m5X.prod[k] += w; s2[k] += w * w;
        }
        for (int k = 0; k < m5X.n(); ++k) { m5X.prodE[k] = std::sqrt(s2[k]) * f; m5X.prod[k] *= f; }
        Divide(m5X);
    }
    EffBins m5Xo; m5Xo.e = eM6;
    {
        const Sel s1720 = [pre](const Row& r) { return pre(r) && r.mee_true >= 0.014f && r.mee_true < 0.020f; };
        FillAcc(m5Xo, 5, s1720, o0, fMpp);
        m5Xo.prod.assign(m5Xo.n(), 0.); m5Xo.prodE.assign(m5Xo.n(), 0.);
        const double f = c.Wprod[5] / c.Ngen;
        for (size_t i = 0; i < c.gen.q.size(); ++i) {
            if (!(c.gen.q[i] >= 0.014 && c.gen.q[i] < 0.020)) continue;
            const int k = EdgeBin(m5Xo.e, c.gen.m[i]); if (k < 0) continue;
            m5Xo.prod[k] += 1;
        }
        for (int k = 0; k < m5Xo.n(); ++k) { m5Xo.prodE[k] = std::sqrt(m5Xo.prod[k]) * f; m5Xo.prod[k] *= f; }
        Divide(m5Xo);
    }
    // class 2 against M(ee) (KW) and M(pipi) (PS x KLOE); class 3 M(pipi); class 4 M(pipi) (PS x wC); class 6 M(pipi) (PS ASSUMED)
    EffBins q2, m2, m3, m4, m6;
    q2.e = eQ2; m2.e = eM23; m3.e = eM23; m4.e = eM4; m6.e = eM6;
    FillAcc(q2, 2, pre, oM, fMee); ProdAnalytic(q2, c.Wprod[2], KWpi0, 2 * kMe * (1 + 1e-9), kMpi0, true); Divide(q2);
    auto rB = [](double M) { return RhoBump(M, true); };
    auto rC = [](double M) { return RhoConv(M, true); };
    FillAcc(m2, 2, pre, oM, fMpp); ProdAnalytic(m2, c.Wprod[2], rB, 2 * kMpi, kMetaHG - kMpi0, false); Divide(m2);
    FillAcc(m3, 3, pre, oM, fMpp); ProdAnalytic(m3, c.Wprod[3], rB, 2 * kMpi, kMetaHG - kMpi0, false); Divide(m3);
    FillAcc(m4, 4, pre, oM, fMpp); ProdAnalytic(m4, c.Wprod[4], rC, 2 * kMpi, kMetaHG, false); Divide(m4);
    FillAcc(m6, 6, pre, oM, fMpp); ProdAnalytic(m6, c.Wprod[6], RhoX17, 2 * kMpi, kMetaGen - kMX17, false); Divide(m6);
    // X17 at its M(ee): one point
    const WSum s6 = SumMany(c.sig, pre, {CompOf(6)}, oM)[0];
    const double eps6 = s6.S / c.Wprod[6], eps6E = eps6 * std::sqrt(s6.V / (s6.S * s6.S) + c.relN * c.relN);

    auto PrintBins = [&](const char* title, const EffBins& b, double xs, const char* xu) {
        P("\n%s\n", title);
        P("  bin [%s]                accepted (W units)         produced (W units)        eps [%%]\n", xu);
        for (int i = 0; i < b.n(); ++i)
            P("  %9.4f - %9.4f   %14.6g +- %10.4g   %14.6g +- %10.4g   %9.5f +- %8.5f\n", b.e[i] * xs, b.e[i + 1] * xs, b.acc[i],
              std::sqrt(b.accV[i]), b.prod[i], b.prodE[i], 100 * b.eps[i], 100 * b.epsE[i]);
    };
    P("\n## 2. against the TRUE M(ee) (produced: class 5 from the generator sample x 100; class 2 Kroll-Wada; X17 one line)\n");
    P("# the threshold bin 2m_e - 1.552 MeV holds the generator values below 2 m_e (4.9 %% of the sample, stored precision): the first two bins\n");
    P("#   of class 5 are affected by the generator's stored precision (UNRELIABLE at the 10 %% level)\n");
    PrintBins("class 5, V3 weights (accepted and produced weighted)", q5V, 1e3, "MeV");
    PrintBins("class 5, no physics weights (generator and accepted unweighted)", q5G, 1e3, "MeV");
    PrintBins("class 2 (bump pi0 Dalitz), produced = Kroll-Wada (HGeant pi0 Dalitz ASSUMED QED)", q2, 1e3, "MeV");
    P("\nclass 6 (X17, M(ee) = 17.000 MeV): eps = %.5f +- %.5f %%   (class 5 V3 at true M(ee) 15.6-18.7 MeV in the table above for comparison)\n",
      100 * eps6, 100 * eps6E);
    P("\n## 3. against the TRUE M(pipi)\n");
    PrintBins("class 5, V3", m5V, 1, "GeV/c2");
    PrintBins("class 5, no physics weights", m5G, 1, "GeV/c2");
    PrintBins("class 5, V3, true M(ee) 14-20 MeV only (for the comparison with X17)", m5X, 1, "GeV/c2");
    PrintBins("class 5, no physics weights (flat e+ angle as the X17 of PLUTO, ASSUMED isotropic), true M(ee) 14-20 MeV only", m5Xo, 1, "GeV/c2");
    PrintBins("class 6 (X17), produced M(pipi) = 3-body phase space ASSUMED (PLUTO X17 M(pipi) not available)", m6, 1, "GeV/c2");
    PrintBins("class 2 (bump pi0 Dalitz), produced = phase space x KLOE-2 wB", m2, 1, "GeV/c2");
    PrintBins("class 3 (bump conversion), produced = per eta -> pi+pi-pi0 decay, phase space x KLOE-2 wB", m3, 1, "GeV/c2");
    PrintBins("class 4 (peak conversion), produced = per eta -> pi+pi-gamma decay, phase space x |F_V|^2 p*^2 k^2", m4, 1, "GeV/c2");

    // true OA and true p(ee): the produced distribution is not available -> d eps / dx
    const std::vector<double> eOA = LogEdges(0.003, 180., 26);
    const std::vector<double> eP = LinEdges(0., 3000., 30);
    const std::vector<double> ePc = LinEdges(0., 3000., 15);
    std::vector<EffBins> bOA(kNC), bP(kNC);
    for (int i = 0; i < kNC; ++i) {
        const int cls = kC[i].cls;
        bOA[i].e = eOA; bP[i].e = (cls == 5 || cls == 6) ? eP : ePc;
        FillAcc(bOA[i], cls, pre, oM, fOA); FillAcc(bP[i], cls, pre, oM, fPee);
    }
    P("\n## 4. against the TRUE OA and the TRUE e+e- pair momentum: d eps / dx = accepted per produced decay per unit of x\n");
    P("# The produced distributions in OA and p(e+e-) need the eta lab momentum, which the generator sample does not hold (only rest-frame\n");
    P("#   M(ee), M(pipi), cos th*_e+, cos th*_pi+, M(eta)) and the SMASH output is not on this machine: the ratio accepted / produced\n");
    P("#   per bin cannot be formed. d eps / d log10(OA/deg) and d eps / dp [1/(GeV/c)] integrate to eps of section 1.\n");
    P("class          | d eps/d log10(OA) per OA bin [%%]  (bins: log10 OA from %.3f to %.0f deg, %d bins)\n", eOA.front(), eOA.back(), (int)eOA.size() - 1);
    for (int i = 0; i < kNC; ++i) {
        P("%-14s |", kC[i].key);
        for (int k = 0; k < bOA[i].n(); ++k) {
            const double dl = std::log10(bOA[i].e[k + 1] / bOA[i].e[k]);
            P(" %.4g", 100 * bOA[i].acc[k] / c.Wprod[kC[i].cls] / dl);
        }
        P("\n");
    }
    // quantiles of true OA and p(ee) of the accepted rows
    P("\naccepted rows ('pre', V3): weighted quantiles of the true OA [deg] and the true p(e+e-) [MeV/c]\n");
    P("class          | OA: q05    median   q95     | p(ee): q05   median   q95   | N_eff\n");
    for (int i = 0; i < kNC; ++i) {
        std::vector<std::pair<double, double>> vo, vp;
        for (auto& a : Collect(kC[i].cls, pre, oM, fOA, fPee)) { vo.push_back({a[0], a[2]}); vp.push_back({a[1], a[2]}); }
        const Quant qo = Quantiles(vo), qp = Quantiles(vp);
        P("%-14s | %8.3f %8.3f %8.3f | %8.0f %8.0f %8.0f | %8.0f\n", kC[i].key, qo.q05, qo.med, qo.q95, qp.q05, qp.med, qp.q95, NeffOf(qo));
    }

    // ---------------- figures
    // (1) eps against the true masses
    {
        TCanvas* cv = new TCanvas(Nm("ceff").c_str(), "", 1700, 1400);
        cv->Divide(2, 2, 0.003, 0.003);
        // (a) M(ee)
        cv->cd(1); PadM(0.15, 0.03, 0.14, 0.09); gPad->SetLogx();
        TH1D* fa = H1("fa", 10, 2 * kMe, QMax());
        AxisT(fa, "M(e^{+}e^{-})_{true} [GeV/c^{2}]", "acceptance #times efficiency [%]");
        fa->SetMinimum(0); fa->SetMaximum(2.0); fa->Draw("axis");
        auto pc = [](std::vector<double> v) { for (double& x : v) x *= 100; return v; };
        TGraphAsymmErrors* g5V = BinGraph(q5V.e, pc(q5V.eps), pc(q5V.epsE), true);
        TGraphAsymmErrors* g5G = BinGraph(q5G.e, pc(q5G.eps), pc(q5G.epsE), true);
        TGraphAsymmErrors* g2 = BinGraph(q2.e, pc(q2.eps), pc(q2.epsE), true);
        Sty(g5V, 633, 20, 1.4); Sty(g5G, 1, 24, 1.3); Sty(g2, 807, 21, 1.5, 3);
        g5G->Draw("p same"); g2->Draw("p same"); g5V->Draw("p same");
        TGraphAsymmErrors* g6 = new TGraphAsymmErrors(1); g6->SetPoint(0, kMX17, 100 * eps6); g6->SetPointError(0, 0, 0, 100 * eps6E, 100 * eps6E);
        Sty(g6, 635, 29, 3.0); g6->Draw("p same");
        TLegend* l = Leg(0.18, 0.64, 0.95, 0.89, 0.045);
        l->AddEntry(g5V, "#eta#rightarrow#pi^{+}#pi^{-}e^{+}e^{-} (5), V3", "pe");
        l->AddEntry(g5G, "the same, no physics weights", "pe");
        l->AddEntry(g6, "X17 (6), M(e^{+}e^{-}) = 17 MeV/c^{2}", "p");
        l->AddEntry(g2, "bump, #pi^{0} Dalitz pair (2)", "pe");
        l->Draw();
        Txt(0.15, 0.93, "(a) per produced decay vs true M(e^{+}e^{-}), pre", 0.055);
        Txt(0.18, 0.17, "first two bins: generator stored precision", 0.042);
        // (b) M(pipi), log
        cv->cd(2); PadM(0.15, 0.03, 0.14, 0.09); gPad->SetLogy();
        TH1D* fb = H1("fb", 10, 0.27, 0.56);
        AxisT(fb, "M(#pi^{+}#pi^{-})_{true} [GeV/c^{2}]", "acceptance #times efficiency [%]");
        fb->SetMinimum(3e-4); fb->SetMaximum(1000); fb->Draw("axis");
        TGraphAsymmErrors* h5 = BinGraph(m5V.e, pc(m5V.eps), pc(m5V.epsE), false);
        TGraphAsymmErrors* h6 = BinGraph(m6.e, pc(m6.eps), pc(m6.epsE), false);
        TGraphAsymmErrors* h2 = BinGraph(m2.e, pc(m2.eps), pc(m2.epsE), false);
        TGraphAsymmErrors* h3 = BinGraph(m3.e, pc(m3.eps), pc(m3.epsE), false);
        TGraphAsymmErrors* h4 = BinGraph(m4.e, pc(m4.eps), pc(m4.epsE), false);
        Sty(h5, 633, 20, 1.3); Sty(h6, 635, 30, 1.8); Sty(h2, 807, 21, 1.5, 3); Sty(h3, 902, 25, 1.5, 3); Sty(h4, 618, 22, 1.7, 3);
        for (TGraph* gg : {(TGraph*)h2, (TGraph*)h3, (TGraph*)h4, (TGraph*)h6, (TGraph*)h5}) gg->Draw("p same");
        TLegend* lb = Leg(0.36, 0.63, 0.97, 0.90, 0.042);
        lb->AddEntry(h5, "direct (5), V3", "pe"); lb->AddEntry(h6, "X17 (6), PS ASSUMED", "pe");
        lb->AddEntry(h2, "bump Dalitz (2)", "pe");
        lb->AddEntry(h3, "bump conv. (3), per #eta#rightarrow#pi^{+}#pi^{-}#pi^{0}", "pe");
        lb->AddEntry(h4, "peak conv. (4), per #eta#rightarrow#pi^{+}#pi^{-}#gamma", "pe");
        lb->Draw();
        Txt(0.15, 0.93, "(b) per produced decay vs true M(#pi^{+}#pi^{-}), pre", 0.055);
        // (c) ratio V3 / off for class 5
        cv->cd(3); PadM(0.15, 0.03, 0.14, 0.09); gPad->SetLogx();
        // ratio = <f>_accepted / <w5 V3>_generated per bin (f = V3 factor of the row): errors from the spread of f and w5
        // ratio = <f>_accepted / <w5>_generated per bin (f = weight factor of the row): errors from the spread of f and w5.
        // V3 (main, with the decay angles) and V2 (same M(ee), M(pipi) model, generator's flat e+ angle): V2/off isolates the
        // in-bin mass effects, V3/V2 the decay angles
        std::vector<double> rq(q5V.n(), 0.), rqe(q5V.n(), 0.), r2(q5V.n(), 0.), r2e(q5V.n(), 0.);
        {
            P("\nclass 5: eps(V)/eps(off) per true M(ee) bin = <f>_accepted / <w5>_generated, V = V3 (main) and V2 (no decay-angle factor)\n");
            P("  bin [MeV]                 V3/off              V2/off              V3/V2\n");
            for (int iv = 0; iv < 2; ++iv) {
                const int code = iv == 0 ? kPhysMain : kPhysV2;
                std::vector<double> af(q5V.n(), 0.), af2(q5V.n(), 0.), an(q5V.n(), 0.), gw(q5V.n(), 0.), gw2(q5V.n(), 0.), gn(q5V.n(), 0.);
                for (auto& a : Collect(5, pre, o0, fMee, [code](const Row& r) { return RowFactor(r, code); })) {
                    const int k = EdgeBin(q5V.e, a[0]); if (k < 0) continue;
                    af[k] += a[1]; af2[k] += a[1] * a[1]; an[k] += 1;
                }
                double sgw = 0;
                std::vector<double> wv(c.gen.q.size());
                for (size_t j = 0; j < c.gen.q.size(); ++j) {
                    wv[j] = iv == 0 ? c.gen.w3[j] : W5(c.gen.q[j], c.gen.m[j], kV2, 0);
                    sgw += wv[j];
                }
                const double norm = c.Ngen / sgw;
                for (size_t j = 0; j < c.gen.q.size(); ++j) {
                    const int k = EdgeBin(q5V.e, std::max(c.gen.q[j], 2 * kMe * 1.0000001)); if (k < 0) continue;
                    const double w = wv[j] * norm; gw[k] += w; gw2[k] += w * w; gn[k] += 1;
                }
                for (int k = 0; k < q5V.n(); ++k) {
                    if (an[k] < 5 || gn[k] < 5) continue;
                    const double mf = af[k] / an[k], sf = std::sqrt(std::max(0., af2[k] / an[k] - mf * mf) / an[k]);
                    const double mg = gw[k] / gn[k], sg = std::sqrt(std::max(0., gw2[k] / gn[k] - mg * mg) / gn[k]);
                    (iv == 0 ? rq : r2)[k] = mf / mg; (iv == 0 ? rqe : r2e)[k] = (mf / mg) * std::hypot(sf / mf, sg / mg);
                }
            }
            for (int k = 0; k < q5V.n(); ++k)
                if (rq[k] > 0 && r2[k] > 0)
                    P("  %8.3f - %8.3f    %.4f +- %.4f    %.4f +- %.4f    %.4f\n", 1e3 * q5V.e[k], 1e3 * q5V.e[k + 1], rq[k], rqe[k], r2[k], r2e[k], rq[k] / r2[k]);
        }
        TH1D* fcc = H1("fcc", 10, 2 * kMe, QMax());
        AxisT(fcc, "M(e^{+}e^{-})_{true} [GeV/c^{2}]", "#varepsilon(V3) / #varepsilon(no weights)");
        fcc->SetMinimum(0.6); fcc->SetMaximum(1.3); fcc->Draw("axis");
        TGraphAsymmErrors* gr = BinGraph(q5V.e, rq, rqe, true); Sty(gr, 633, 20, 1.4); gr->Draw("p same");
        TGraphAsymmErrors* gr2 = BinGraph(q5V.e, r2, r2e, true); Sty(gr2, 1, 24, 1.3); gr2->Draw("p same");
        TLine one; one.SetLineStyle(2); one.SetLineWidth(2); one.DrawLine(2 * kMe, 1, QMax(), 1);
        TLegend* lr = Leg(0.18, 0.74, 0.95, 0.89, 0.045);
        lr->AddEntry(gr, "V3 (main) / no weights", "pe");
        lr->AddEntry(gr2, "V2 (no angle factor) / no weights", "pe");
        lr->Draw();
        Txt(0.15, 0.93, "(c) direct: #varepsilon ratio at fixed M(e^{+}e^{-})", 0.055);
        Txt(0.18, 0.19, "V3/V2: QED 1 + cos^{2}#theta*_{e} lowers the acceptance", 0.042);
        // (d) X17 against the direct decay at the same M(ee), in M(pipi)
        cv->cd(4); PadM(0.15, 0.03, 0.14, 0.09);
        TH1D* fd = H1("fd", 10, 0.27, 0.54);
        AxisT(fd, "M(#pi^{+}#pi^{-})_{true} [GeV/c^{2}]", "acceptance #times efficiency [%]");
        fd->SetMinimum(0); fd->SetMaximum(2.2); fd->Draw("axis");
        TGraphAsymmErrors* x5 = BinGraph(m5X.e, pc(m5X.eps), pc(m5X.epsE), false);
        TGraphAsymmErrors* x6 = BinGraph(m6.e, pc(m6.eps), pc(m6.epsE), false);
        TGraphAsymmErrors* x5o = BinGraph(m5Xo.e, pc(m5Xo.eps), pc(m5Xo.epsE), false);
        Sty(x5, 633, 20, 1.4); Sty(x6, 635, 30, 2.0, 2); Sty(x5o, 1, 24, 1.4);
        x5o->Draw("p same"); x5->Draw("p same"); x6->Draw("p same");
        TLegend* ldd = Leg(0.50, 0.66, 0.97, 0.89, 0.045);
        ldd->AddEntry(x5, "direct (5), V3", "pe");
        ldd->AddEntry(x5o, "direct (5), no weights", "pe");
        ldd->AddEntry(x6, "X17 (6), PS ASSUMED", "pe");
        ldd->Draw();
        Txt(0.15, 0.93, "(d) X17 vs direct, true M(e^{+}e^{-}) 14-20 MeV/c^{2}", 0.055);
        Txt(0.18, 0.20, Form("X17 total: %.3f #pm %.3f %%", 100 * eps6, 100 * eps6E), 0.045, 635);
        cv->cd(0);
        SaveFig(cv, "s5_acc_eff_mass");
    }
    // (2) d eps / dx in true OA and true p(ee)
    {
        TCanvas* cv = new TCanvas(Nm("ckin").c_str(), "", 1700, 750);
        cv->Divide(2, 1, 0.003, 0.003);
        cv->cd(1); PadM(0.15, 0.03, 0.15, 0.10); gPad->SetLogx(); gPad->SetLogy();
        TH1D* fa = H1("fka", 10, eOA.front(), eOA.back());
        AxisT(fa, "OA_{true} [deg]", "d#varepsilon / dlog_{10}(OA_{true}/deg) [%]");
        fa->SetMinimum(1e-5); fa->SetMaximum(20); fa->Draw("axis");
        TLegend* l = Leg(0.17, 0.66, 0.60, 0.89, 0.042);
        for (int i = 0; i < kNC; ++i) {
            std::vector<double> y(bOA[i].n()), ey(bOA[i].n());
            for (int k = 0; k < bOA[i].n(); ++k) {
                const double dl = std::log10(bOA[i].e[k + 1] / bOA[i].e[k]);
                y[k] = 100 * bOA[i].acc[k] / C().Wprod[kC[i].cls] / dl; ey[k] = 100 * std::sqrt(bOA[i].accV[k]) / C().Wprod[kC[i].cls] / dl;
            }
            TGraphAsymmErrors* gg = BinGraph(bOA[i].e, y, ey, true);
            Sty(gg, kC[i].col, kC[i].ms, kC[i].cls == 6 ? 2.0 : 1.3); gg->Draw("p same");
            l->AddEntry(gg, kC[i].key, "pe");
        }
        l->Draw();
        Txt(0.15, 0.93, "(a) d#varepsilon/dx vs true OA, pre", 0.055);
        Txt(0.17, 0.20, "produced OA distribution not available: d#varepsilon/dx", 0.042);
        cv->cd(2); PadM(0.15, 0.03, 0.15, 0.10); gPad->SetLogy();
        TH1D* fb = H1("fkb", 10, 0, 3.0);
        AxisT(fb, "p(e^{+}e^{-})_{true} [GeV/c]", "d#varepsilon / dp(e^{+}e^{-})_{true} [%/(GeV/c)]");
        fb->SetMinimum(1e-5); fb->SetMaximum(20); fb->Draw("axis");
        TLegend* lb = Leg(0.55, 0.66, 0.95, 0.89, 0.042);
        for (int i = 0; i < kNC; ++i) {
            std::vector<double> y(bP[i].n()), ey(bP[i].n());
            for (int k = 0; k < bP[i].n(); ++k) {
                const double dp = (bP[i].e[k + 1] - bP[i].e[k]) * 1e-3;
                y[k] = 100 * bP[i].acc[k] / C().Wprod[kC[i].cls] / dp; ey[k] = 100 * std::sqrt(bP[i].accV[k]) / C().Wprod[kC[i].cls] / dp;
            }
            TGraphAsymmErrors* gg = BinGraph(bP[i].e, y, ey, false, 1e-3);
            Sty(gg, kC[i].col, kC[i].ms, kC[i].cls == 6 ? 2.0 : 1.3); gg->Draw("p same");
            lb->AddEntry(gg, kC[i].key, "pe");
        }
        lb->Draw();
        Txt(0.15, 0.93, "(b) d#varepsilon/dx vs true p(e^{+}e^{-}), pre", 0.055);
        cv->cd(0);
        SaveFig(cv, "s5_acc_eff_kin");
    }
    CloseTab();
    printf("PartAcc done in %.0f s\n", sw.RealTime());
}

// =====================================================================================================
// Part 3: resolution, scale, the simulated peak and the PEAK window
// =====================================================================================================
void PartRes() {
    TStopwatch sw;
    Style();
    if (!LoadSim()) return;
    Ctx& c = C();
    const Sel pre = SelPreset(kPre);
    const SimOpt oM = OptMain(), o0 = OptOff();
    OpenTab("s5_acc_resolution.txt");
    P("# s5_acc_resolution: resolution and scale of the reconstructed (energy-loss corrected) M4, M(ee), M(pipi) against the GEANT truth,\n");
    P("#   'pre', owned rows mode A, std + lep + eta, w (PDG) x physics weights V3 (X17: 1 per row). s5_acc.C part 4. SIMULATION ONLY.\n");
    P("# core: Gaussian fit iterated in mean +- 2 sigma (bins with event-group errors); quantiles weighted; tails = weighted fractions outside\n");
    P("#   core mean +- 3 sigma_core and |Delta| > 20 MeV; error of a median ~ 1.25 (q84-q16)/2 / sqrt(N_eff).\n\n");
    // ---- Delta M4
    P("## 1. Delta M4 = M4 - M4_true [MeV]\n");
    P("class            N_eff     median     (q84-q16)/2   core mean        core sigma       low tail   high tail  |D|>20MeV   median M4/M4_true\n");
    struct DR { TH1D* h = nullptr; GFit g; Quant q; };
    std::vector<DR> dm4(kNC);
    for (int i = 0; i < kNC; ++i) {
        const int cls = kC[i].cls;
        Binning b; b.nbin = (cls == 5 || cls == 6) ? 240 : 60; b.lo = -60; b.hi = 60;   // 0.5 MeV bins, 2 MeV for the std classes
        b.x = [](const Row& r) { return 1e3 * (r.m4 - r.m4_true); };
        const WSpec sp = Fill(c.sig, pre, CompOf(cls), oM, b);
        dm4[i].h = sp.Hist(Nm("dm4").c_str(), 0);
        std::vector<std::pair<double, double>> v, vr;
        for (auto& a : Collect(cls, pre, oM, [](const Row& r) { return 1e3 * (r.m4 - r.m4_true); },
                               [](const Row& r) { return r.m4_true > 0 ? r.m4 / r.m4_true : -1.; })) {
            v.push_back({a[0], a[2]}); vr.push_back({a[1], a[2]});
        }
        dm4[i].q = Quantiles(v);
        const Quant qr = Quantiles(vr);
        const double s0 = 0.5 * (dm4[i].q.q84 - dm4[i].q.q16);
        dm4[i].g = FitCore(dm4[i].h, dm4[i].q.med, std::min(s0, 8.));
        double lo = 0, hi = 0, t20 = 0, tot = 0;
        for (auto& p : v) {
            tot += p.second;
            if (p.first < dm4[i].g.mu - 3 * dm4[i].g.sig) lo += p.second;
            if (p.first > dm4[i].g.mu + 3 * dm4[i].g.sig) hi += p.second;
            if (std::fabs(p.first) > 20) t20 += p.second;
        }
        const double ne = NeffOf(dm4[i].q);
        P("%-15s %8.0f   %7.2f +- %4.2f   %6.2f    %6.2f +- %4.2f    %6.2f +- %4.2f    %6.4f     %6.4f     %6.4f     %.5f +- %.5f\n", kC[i].key, ne,
          dm4[i].q.med, 1.25 * s0 / std::sqrt(std::max(ne, 1.)), s0, dm4[i].g.mu, dm4[i].g.emu, dm4[i].g.sig, dm4[i].g.esig, lo / tot, hi / tot, t20 / tot, qr.med,
          1.25 * 0.5 * (qr.q84 - qr.q16) / std::sqrt(std::max(ne, 1.)));
        if (cls == 5 || cls == 4 || cls == 2)
            c.summary.push_back(Form("Delta M4 %-15s median %.2f MeV, core sigma %.2f +- %.2f MeV, |D|>20 MeV %.3f, M4/M4true %.5f", kC[i].key, dm4[i].q.med,
                                     dm4[i].g.sig, dm4[i].g.esig, t20 / tot, qr.med));
    }
    // ---- Delta M(ee) against the true M(ee)
    P("\n## 2. Delta M(ee) = M(ee) - M(ee)_true [MeV] in bins of the true M(ee) (class 5 V3; class 2; X17)\n");
    P("class / true M(ee) bin [MeV]      N_eff     median Delta    (q84-q16)/2    core sigma          core sigma / M_true\n");
    const std::vector<double> eQr = {0.001022, 0.002, 0.004, 0.007, 0.010, 0.014, 0.0155, 0.0185, 0.022, 0.030, 0.045, 0.070, 0.110, 0.270};
    std::vector<double> xsC, sC, sCe, qC, xsCe;
    for (int k = 0; k + 1 < (int)eQr.size(); ++k) {
        const double a = eQr[k], bq = eQr[k + 1];
        const Sel s = [pre, a, bq](const Row& r) { return pre(r) && r.mee_true >= a && r.mee_true < bq; };
        std::vector<std::pair<double, double>> v;
        for (auto& x : Collect(5, s, oM, [](const Row& r) { return 1e3 * (r.mee - r.mee_true); })) v.push_back({x[0], x[2]});
        if (v.size() < 50) continue;
        const Quant q = Quantiles(v);
        const double s0 = std::max(0.5 * (q.q84 - q.q16), 0.05);
        Binning b; b.nbin = 120; b.lo = q.med - 6 * s0; b.hi = q.med + 6 * s0;
        b.x = [](const Row& r) { return 1e3 * (r.mee - r.mee_true); };
        TH1D* h = Fill(c.sig, s, CompOf(5), oM, b).Hist(Nm("dmee").c_str(), 0);
        const GFit g = FitCore(h, q.med, s0);
        const double mt = std::sqrt(a * bq) * 1e3;
        P("direct  %7.2f - %7.2f        %9.0f   %8.3f        %8.3f       %6.3f +- %5.3f     %.4f\n", a * 1e3, bq * 1e3, NeffOf(q), q.med, s0, g.sig, g.esig,
          g.sig / mt);
        xsC.push_back(mt); sC.push_back(g.sig); sCe.push_back(g.esig); qC.push_back(s0); xsCe.push_back(0);
        delete h;
    }
    // class 2 in a few bins
    std::vector<double> x2, s2v, s2e;
    {
        const std::vector<double> e2 = {0.001022, 0.004, 0.010, 0.025, 0.135};
        for (int k = 0; k + 1 < (int)e2.size(); ++k) {
            const double a = e2[k], bq = e2[k + 1];
            const Sel s = [pre, a, bq](const Row& r) { return pre(r) && r.mee_true >= a && r.mee_true < bq; };
            std::vector<std::pair<double, double>> v;
            for (auto& x : Collect(2, s, oM, [](const Row& r) { return 1e3 * (r.mee - r.mee_true); })) v.push_back({x[0], x[2]});
            if (v.size() < 30) continue;
            const Quant q = Quantiles(v);
            const double s0 = std::max(0.5 * (q.q84 - q.q16), 0.05);
            P("bump(2) %7.2f - %7.2f        %9.0f   %8.3f        %8.3f       (quantiles only, N_eff small)\n", a * 1e3, bq * 1e3, NeffOf(q), q.med, s0);
            x2.push_back(std::sqrt(a * bq) * 1e3); s2v.push_back(s0); s2e.push_back(s0 / std::sqrt(2 * std::max(NeffOf(q), 2.)));
        }
    }
    // X17: Delta M(ee) and the window
    TH1D* hX = nullptr; GFit gX; Quant qX;
    {
        std::vector<std::pair<double, double>> v, vm;
        for (auto& x : Collect(6, pre, oM, [](const Row& r) { return 1e3 * (r.mee - r.mee_true); }, [](const Row& r) { return 1e3 * r.mee; })) {
            v.push_back({x[0], x[2]}); vm.push_back({x[1], x[2]});
        }
        qX = Quantiles(v);
        const Quant qm = Quantiles(vm);
        Binning b; b.nbin = 100; b.lo = -25; b.hi = 25;
        b.x = [](const Row& r) { return 1e3 * (r.mee - r.mee_true); };
        hX = Fill(c.sig, pre, CompOf(6), oM, b).Hist(Nm("dmX").c_str(), 0);
        gX = FitCore(hX, qX.med, 0.5 * (qX.q84 - qX.q16));
        // smallest windows around the median of the reconstructed M(ee) containing 68 / 90 / 95 %
        std::vector<std::pair<double, double>> s = vm; std::sort(s.begin(), s.end());
        double tot = 0; for (auto& p : s) tot += p.second;
        auto shortest = [&](double frac, double& lo, double& hi) {
            // two-pointer over the sorted values
            size_t j = 0; double acc = 0, best = 1e9; lo = hi = 0;
            for (size_t i = 0; i < s.size(); ++i) {
                while (j < s.size() && acc < frac * tot) { acc += s[j].second; ++j; }
                if (acc < frac * tot) break;
                if (s[j - 1].first - s[i].first < best) { best = s[j - 1].first - s[i].first; lo = s[i].first; hi = s[j - 1].first; }
                acc -= s[i].second;
            }
        };
        double l68, h68, l90, h90, l95, h95;
        shortest(0.68, l68, h68); shortest(0.90, l90, h90); shortest(0.95, l95, h95);
        c.x17lo = l90; c.x17hi = h90;
        P("X17     true 17.000 (all rows)   %9.0f   %8.3f        %8.3f       %6.3f +- %5.3f     %.4f   (core mean %.3f MeV)\n", NeffOf(qX), qX.med,
          0.5 * (qX.q84 - qX.q16), gX.sig, gX.esig, gX.sig / 17., gX.mu);
        P("X17 reconstructed M(ee): median %.3f MeV; shortest intervals containing 68 %% [%.2f, %.2f], 90 %% [%.2f, %.2f], 95 %% [%.2f, %.2f] MeV/c2\n", qm.med,
          l68, h68, l90, h90, l95, h95);
        P("   -> X17 window suggestion (90 %% of the X17 rows): %.1f - %.1f MeV/c2 (half width %.2f MeV); core sigma %.2f MeV at 17 MeV\n", l90, h90,
          0.5 * (h90 - l90), gX.sig);
        c.summary.push_back(Form("X17 M(ee): core sigma %.3f +- %.3f MeV, median shift %.3f MeV, 90%% interval %.2f-%.2f MeV/c2", gX.sig, gX.esig, qX.med, l90, h90));
    }
    // ---- Delta M(pipi)
    P("\n## 3. Delta M(pipi) = M(pipi) - M(pipi)_true [MeV]\n");
    P("class            N_eff     median     (q84-q16)/2   core sigma\n");
    std::vector<TH1D*> hpp(kNC, nullptr);
    for (int i = 0; i < kNC; ++i) {
        const int cls = kC[i].cls;
        Binning b; b.nbin = (cls == 5 || cls == 6) ? 160 : 40; b.lo = -40; b.hi = 40;
        b.x = [](const Row& r) { return 1e3 * (r.mpipi - r.mpipi_true); };
        hpp[i] = Fill(c.sig, pre, CompOf(cls), oM, b).Hist(Nm("dpp").c_str(), 0);
        std::vector<std::pair<double, double>> v;
        for (auto& x : Collect(cls, pre, oM, [](const Row& r) { return 1e3 * (r.mpipi - r.mpipi_true); })) v.push_back({x[0], x[2]});
        const Quant q = Quantiles(v);
        const GFit g = FitCore(hpp[i], q.med, 0.5 * (q.q84 - q.q16));
        P("%-15s %8.0f   %7.3f      %7.3f       %6.3f +- %5.3f\n", kC[i].key, NeffOf(q), q.med, 0.5 * (q.q84 - q.q16), g.sig, g.esig);
    }
    P("class 5, against the true M(pipi):\n");
    for (double a = 0.28; a < 0.54; a += 0.04) {
        const double bq = a + 0.04;
        const Sel s = [pre, a, bq](const Row& r) { return pre(r) && r.mpipi_true >= a && r.mpipi_true < bq; };
        std::vector<std::pair<double, double>> v;
        for (auto& x : Collect(5, s, oM, [](const Row& r) { return 1e3 * (r.mpipi - r.mpipi_true); })) v.push_back({x[0], x[2]});
        if (v.size() < 50) continue;
        const Quant q = Quantiles(v);
        P("  true M(pipi) %.2f-%.2f: median %6.3f MeV, (q84-q16)/2 %6.3f MeV, N_eff %.0f\n", a, bq, q.med, 0.5 * (q.q84 - q.q16), NeffOf(q));
    }

    // ---- the simulated peak with the physics weights: position, width, PEAK window
    P("\n## 4. simulated eta peak (pre, reconstructed corrected M4) with the physics weights V3; window = mean +- 1.9 sigma\n");
    P("# fixed PEAK window of this study: %.4f - %.4f GeV/c2 (0.5460 +- 1.9 x 0.0067)\n", kPeakLo, kPeakHi);
    P("sample     opt   bin [MeV]  range   mean [GeV/c2]          sigma [MeV]       chi2/ndf     window (mean +- 1.9 sigma)   d(lo) d(hi) [MeV]\n");
    struct PF { std::string name; double mu, emu, sig, esig, lo, hi; std::string lab; };
    std::vector<PF> pf;
    TH1D* hpk[3] = {nullptr, nullptr, nullptr};
    const char* pname[3] = {"total", "direct", "conversion"};
    const Comp pcomp[3] = {CompPeakTotal(), CompPeakDirect(), CompPeakConv()};
    double mainMu = 0, mainSig = 0, mainEmu = 0, mainEsig = 0;
    std::vector<double> spreadMu, spreadSig;
    for (int ip = 0; ip < 3; ++ip)
        for (int io = 0; io < 2; ++io)
            for (double bw : {0.5, 1.0, 2.0})
                for (double ns : {1.5, 2.0, 2.5}) {
                    if (io == 1 && !(bw == 1.0 && ns == 2.0)) continue;
                    Binning b; b.lo = 0.49; b.hi = 0.61; b.nbin = (int)std::lround((b.hi - b.lo) / (bw * 1e-3));
                    const WSpec sp = Fill(c.sig, pre, pcomp[ip], io == 0 ? oM : o0, b);
                    TH1D* h = sp.Hist(Nm("hpk").c_str(), 0);
                    const GFit g = FitCore(h, 0.5462, 0.0065, ns);
                    const double lo = g.mu - 1.9 * g.sig, hi = g.mu + 1.9 * g.sig;
                    P("%-10s %-4s  %5.1f     +-%.1f   %.5f +- %.5f    %5.2f +- %4.2f    %6.1f/%-3d    %.4f - %.4f            %+5.1f %+5.1f\n", pname[ip],
                      io ? "off" : "V3", bw, ns, g.mu, g.emu, 1e3 * g.sig, 1e3 * g.esig, g.chi2, g.ndf, lo, hi, 1e3 * (lo - kPeakLo), 1e3 * (hi - kPeakHi));
                    if (io == 0 && bw == 1.0 && ns == 2.0) {
                        if (!hpk[ip]) hpk[ip] = (TH1D*)h->Clone(Nm("hpkc").c_str());
                        if (ip == 0) { mainMu = g.mu; mainSig = g.sig; mainEmu = g.emu; mainEsig = g.esig; }
                    }
                    if (io == 0 && ip == 0) { spreadMu.push_back(g.mu); spreadSig.push_back(g.sig); }
                    if (io == 0) {
                        pf.push_back({Form("%s %.1f MeV +-%.1f", pname[ip], bw, ns), g.mu, g.emu, g.sig, g.esig, lo, hi, ""});
                        pf.back().lab = Form("%s %.1f, %.1f#sigma", ip == 0 ? "total" : ip == 1 ? "direct" : "conv.", bw, ns);
                    }
                    delete h;
                }
    const double muLo = *std::min_element(spreadMu.begin(), spreadMu.end()), muHi = *std::max_element(spreadMu.begin(), spreadMu.end());
    const double sgLo = *std::min_element(spreadSig.begin(), spreadSig.end()), sgHi = *std::max_element(spreadSig.begin(), spreadSig.end());
    const double wLo = mainMu - 1.9 * mainSig, wHi = mainMu + 1.9 * mainSig;
    P("\nMAIN (total = direct + conversion, V3, 1 MeV bins, +-2 sigma): mean %.5f +- %.5f (stat) [fit-choice range %.5f - %.5f] GeV/c2,\n", mainMu,
      mainEmu, muLo, muHi);
    P("   sigma %.2f +- %.2f (stat) [fit-choice range %.2f - %.2f] MeV -> window %.4f - %.4f GeV/c2\n", 1e3 * mainSig, 1e3 * mainEsig, 1e3 * sgLo,
      1e3 * sgHi, wLo, wHi);
    P("   difference to the fixed window: lower edge %+.1f MeV, upper edge %+.1f MeV -> %s\n", 1e3 * (wLo - kPeakLo), 1e3 * (wHi - kPeakHi),
      (std::fabs(wLo - kPeakLo) > 0.001 || std::fabs(wHi - kPeakHi) > 0.001) ? "DIFFERS by more than 1 MeV/c2 (see the report)"
                                                                              : "CONFIRMED within 1 MeV/c2");
    // window content of the simulated peak
    {
        const std::vector<Comp> cm = {CompPeakTotal(), CompPeakDirect(), CompPeakConv()};
        const std::vector<WSum> all = SumMany(c.sig, pre, cm, oM);
        const std::vector<WSum> inF = SumMany(c.sig, [pre](const Row& r) { return pre(r) && r.m4 >= kPeakLo && r.m4 < kPeakHi; }, cm, oM);
        const std::vector<WSum> inN = SumMany(c.sig, [pre, wLo, wHi](const Row& r) { return pre(r) && r.m4 >= wLo && r.m4 < wHi; }, cm, oM);
        P("fraction of the simulated peak (all M4) inside the fixed window / the recomputed window: total %.4f / %.4f, direct %.4f / %.4f, conversion %.4f / %.4f\n",
          inF[0].S / all[0].S, inN[0].S / all[0].S, inF[1].S / all[1].S, inN[1].S / all[1].S, inF[2].S / all[2].S, inN[2].S / all[2].S);
        P("sigma in the fixed window [nb]: total %.4f +- %.4f, direct %.4f +- %.4f, conversion %.4f +- %.4f\n", inF[0].S / kLsim, std::sqrt(inF[0].V) / kLsim,
          inF[1].S / kLsim, std::sqrt(inF[1].V) / kLsim, inF[2].S / kLsim, std::sqrt(inF[2].V) / kLsim);
        c.summary.push_back(Form("peak (V3, total): mean %.5f +- %.5f, sigma %.2f +- %.2f MeV (fit range mean %.5f-%.5f, sigma %.2f-%.2f); window %.4f-%.4f; "
                                 "edges %+.1f/%+.1f MeV vs fixed; fixed-window containment %.3f",
                                 mainMu, mainEmu, 1e3 * mainSig, 1e3 * mainEsig, muLo, muHi, 1e3 * sgLo, 1e3 * sgHi, wLo, wHi, 1e3 * (wLo - kPeakLo),
                                 1e3 * (wHi - kPeakHi), inF[0].S / all[0].S));
    }
    CloseTab();

    // ---------------- figures
    {
        TCanvas* cv = new TCanvas(Nm("cres").c_str(), "", 1700, 1400);
        cv->Divide(2, 2, 0.003, 0.003);
        // (a) Delta M4
        cv->cd(1); PadM(0.15, 0.03, 0.14, 0.09); gPad->SetLogy();
        TH1D* fa = H1("fra", 10, -60, 40);
        AxisT(fa, "M4 - M4_{true} [MeV/c^{2}]", "probability density [1/(MeV/c^{2})]");
        fa->SetMinimum(2e-5); fa->SetMaximum(20); fa->Draw("axis");
        TLegend* l = Leg(0.17, 0.62, 0.70, 0.89, 0.042);
        for (int i = 0; i < kNC; ++i) {
            TH1D* h = dm4[i].h; const double s = h->Integral();
            if (s > 0) h->Scale(1. / (s * h->GetBinWidth(1)));
            Sty(h, kC[i].col, kC[i].cls == 5 ? 4 : 3, kC[i].cls == 6 ? 2 : 1);
            h->Draw("hist same");
            l->AddEntry(h, Form("%s: #sigma_{core} %.1f MeV", kC[i].key, dm4[i].g.sig), "l");
        }
        l->Draw();
        Txt(0.15, 0.93, "(a) M4 - M4_{true}, pre", 0.055);
        Txt(0.17, 0.20, "bins: 0.5 MeV/c^{2} (5, 6), 2 MeV/c^{2} (2, 3, 4)", 0.040);
        // (b) M(ee) resolution against the true M(ee)
        cv->cd(2); PadM(0.15, 0.03, 0.14, 0.09); gPad->SetLogx();
        TH1D* fb = H1("frb", 10, 1, 300);
        AxisT(fb, "M(e^{+}e^{-})_{true} [MeV/c^{2}]", "width of M(e^{+}e^{-}) - M(e^{+}e^{-})_{true} [MeV/c^{2}]");
        fb->SetMinimum(0); fb->SetMaximum(10); fb->Draw("axis");
        TGraphErrors* gs = new TGraphErrors((int)xsC.size(), xsC.data(), sC.data(), xsCe.data(), sCe.data());
        TGraph* gq = new TGraph((int)xsC.size(), xsC.data(), qC.data());
        TGraphErrors* g2 = new TGraphErrors((int)x2.size(), x2.data(), s2v.data(), nullptr, s2e.data());
        Sty(gs, 633, 20, 1.5); Sty(gq, 633, 24, 1.5, 2, 2); Sty(g2, 807, 21, 1.6);
        gs->Draw("p same"); gq->Draw("p same"); g2->Draw("p same");
        TGraphErrors* gx = new TGraphErrors(1); gx->SetPoint(0, 17., gX.sig); gx->SetPointError(0, 0, gX.esig); Sty(gx, 635, 29, 3.0); gx->Draw("p same");
        TLegend* lb = Leg(0.17, 0.66, 0.95, 0.89, 0.042);
        lb->AddEntry(gs, "direct (5): Gaussian core #sigma", "pe"); lb->AddEntry(gq, "direct (5): (q_{84} - q_{16})/2", "p");
        lb->AddEntry(gx, Form("X17 (6): core #sigma = %.2f MeV/c^{2}", gX.sig), "p"); lb->AddEntry(g2, "bump #pi^{0} Dalitz (2): (q_{84} - q_{16})/2", "p");
        lb->Draw();
        Txt(0.15, 0.93, "(b) M(e^{+}e^{-}) resolution, pre", 0.055);
        // (c) X17 Delta M(ee)
        cv->cd(3); PadM(0.15, 0.03, 0.14, 0.09); gPad->SetLogy();
        hX->Scale(1. / hX->Integral());
        Sty(hX, 635, 3, 1);
        AxisT(hX, "M(e^{+}e^{-}) - 17 MeV/c^{2} (X17 rows) [MeV/c^{2}]", "fraction per 0.5 MeV/c^{2}");
        hX->SetMinimum(1e-5); hX->SetMaximum(1.0); hX->Draw("hist");
        TF1* fg = new TF1(Nm("fgx").c_str(), "gaus", gX.mu - 2.5 * gX.sig, gX.mu + 2.5 * gX.sig);
        fg->SetParameters(hX->GetBinContent(hX->FindBin(gX.mu)), gX.mu, gX.sig); fg->SetLineColor(1); fg->SetLineStyle(2); fg->SetLineWidth(3);
        fg->Draw("same");
        Txt(0.18, 0.82, Form("core #sigma = %.2f #pm %.2f MeV/c^{2}", gX.sig, gX.esig), 0.048);
        Txt(0.18, 0.75, Form("median shift %.2f MeV/c^{2}", qX.med), 0.048);
        Txt(0.18, 0.68, Form("(q_{84} - q_{16})/2 = %.2f MeV/c^{2}", 0.5 * (qX.q84 - qX.q16)), 0.048);
        Txt(0.15, 0.93, "(c) X17: M(e^{+}e^{-}) - 17 MeV/c^{2}, pre", 0.055);
        // (d) Delta M(pipi)
        cv->cd(4); PadM(0.15, 0.03, 0.14, 0.09); gPad->SetLogy();
        TH1D* fd = H1("frd", 10, -40, 40);
        AxisT(fd, "M(#pi^{+}#pi^{-}) - M(#pi^{+}#pi^{-})_{true} [MeV/c^{2}]", "probability density [1/(MeV/c^{2})]");
        fd->SetMinimum(2e-5); fd->SetMaximum(20); fd->Draw("axis");
        TLegend* ld = Leg(0.17, 0.64, 0.60, 0.89, 0.042);
        for (int i = 0; i < kNC; ++i) {
            TH1D* h = hpp[i]; const double s = h->Integral(); if (s > 0) h->Scale(1. / (s * h->GetBinWidth(1)));
            Sty(h, kC[i].col, kC[i].cls == 5 ? 4 : 3, kC[i].cls == 6 ? 2 : 1); h->Draw("hist same");
            ld->AddEntry(h, kC[i].key, "l");
        }
        ld->Draw();
        Txt(0.15, 0.93, "(d) M(#pi^{+}#pi^{-}) resolution, pre", 0.055);
        cv->cd(0);
        SaveFig(cv, "s5_acc_res");
    }
    {
        TCanvas* cv = new TCanvas(Nm("cpk").c_str(), "", 1700, 750);
        cv->Divide(2, 1, 0.003, 0.003);
        cv->cd(1); PadM(0.15, 0.03, 0.15, 0.10);
        TH1D* ht = hpk[0]; TH1D* hd = hpk[1]; TH1D* hc = hpk[2];
        for (TH1D* h : {ht, hd, hc}) { h->Scale(1. / (kLsim * h->GetBinWidth(1))); }
        AxisT(ht, "M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]", "d#sigma/dM4 [nb/(GeV/c^{2})]");
        ht->GetXaxis()->SetRangeUser(0.50, 0.59); ht->GetXaxis()->SetNdivisions(505);
        Sty(ht, 1, 3, 1, 20, 1.0); Sty(hd, 633, 3, 1); Sty(hc, 618, 3, 1);
        ht->SetMaximum(1.75 * ht->GetMaximum()); ht->SetMinimum(0);
        ht->Draw("pe"); hd->Draw("hist same"); hc->Draw("hist same");
        TF1* fg = new TF1(Nm("fgp").c_str(), "gaus", mainMu - 2 * mainSig, mainMu + 2 * mainSig);
        fg->SetParameters(ht->GetBinContent(ht->FindBin(mainMu)), mainMu, mainSig); fg->SetLineColor(601); fg->SetLineWidth(3);
        ht->Fit(fg, "QR0"); fg->Draw("same");
        const double ymax = ht->GetMaximum();
        TLine lf; lf.SetLineColor(922); lf.SetLineWidth(3); lf.SetLineStyle(2);
        lf.DrawLine(kPeakLo, 0, kPeakLo, ymax); lf.DrawLine(kPeakHi, 0, kPeakHi, ymax);
        TLine ln; ln.SetLineColor(601); ln.SetLineWidth(3); ln.SetLineStyle(1);
        ln.DrawLine(wLo, 0, wLo, 0.8 * ymax); ln.DrawLine(wHi, 0, wHi, 0.8 * ymax);
        TLegend* l = Leg(0.17, 0.66, 0.55, 0.89, 0.042);
        l->AddEntry(ht, "total (4 + 5), V3", "pe"); l->AddEntry(hd, "direct (5)", "l"); l->AddEntry(hc, "conversion (4)", "l");
        l->AddEntry(fg, "Gaussian #pm2#sigma", "l"); l->Draw();
        TLegend* l2 = Leg(0.62, 0.78, 0.97, 0.89, 0.042);
        l2->AddEntry(&lf, "fixed window", "l"); l2->AddEntry(&ln, "#mu #pm 1.9#sigma", "l"); l2->Draw();
        Txt(0.17, 0.54, Form("#mu = %.4f GeV/c^{2}", mainMu), 0.045, 601);
        Txt(0.17, 0.48, Form("#sigma = %.2f MeV/c^{2}", 1e3 * mainSig), 0.045, 601);
        Txt(0.15, 0.93, "(a) simulated #eta peak, weights V3, pre", 0.055);
        // (b) window edges per fit choice (total V3)
        cv->cd(2); PadM(0.40, 0.03, 0.15, 0.10);
        std::vector<PF> sel;
        for (auto& q : pf) if (q.name.find("total") == 0 || q.name.find("1.0 MeV +-2.0") != std::string::npos) sel.push_back(q);
        pf.swap(sel);
        const int n = (int)pf.size();
        TH1D* fb = H1("fpb", 10, 0.531, 0.561);
        AxisT(fb, "window #mu #pm 1.9#sigma [GeV/c^{2}]", "");
        fb->GetXaxis()->SetNdivisions(506);
        fb->SetMinimum(-0.5); fb->SetMaximum(n - 0.5); fb->GetYaxis()->SetLabelSize(0); fb->GetYaxis()->SetTickLength(0); fb->Draw("axis");
        for (int k = 0; k < n; ++k) {
            const int col = pf[k].name.find("total") == 0 ? 1 : pf[k].name.find("direct") == 0 ? 633 : 618;
            TLine lk; lk.SetLineColor(col); lk.SetLineWidth(4); lk.DrawLine(pf[k].lo, k, pf[k].hi, k);
            TMarker mk(pf[k].mu, k, 20); mk.SetMarkerColor(col); mk.SetMarkerSize(1.2); mk.DrawClone();
            TLatex t; t.SetTextSize(0.040); t.SetTextAlign(32); t.SetTextColor(col); t.SetTextFont(42);
            t.DrawLatex(0.5305, k, pf[k].lab.c_str());
        }
        TLine lf2; lf2.SetLineColor(922); lf2.SetLineWidth(3); lf2.SetLineStyle(2);
        lf2.DrawLine(kPeakLo, -0.5, kPeakLo, n - 0.5); lf2.DrawLine(kPeakHi, -0.5, kPeakHi, n - 0.5);
        Txt(0.03, 0.93, "(b) #mu #pm 1.9#sigma per fit choice (bin [MeV], range)", 0.050);
        Txt(0.03, 0.04, "dashed: fixed window", 0.045, 922);
        cv->cd(0);
        SaveFig(cv, "s5_acc_peak");
    }
    printf("PartRes done in %.0f s\n", sw.RealTime());
}

// =====================================================================================================
// Part 4: reconstructed OA against M(ee) and the OA separation
// =====================================================================================================
void PartOA() {
    TStopwatch sw;
    Style();
    if (!LoadSim()) return;
    Ctx& c = C();
    const Sel pre = SelPreset(kPre);
    const SimOpt oM = OptMain();
    OpenTab("s5_acc_oa.txt");
    P("# s5_acc_oa: reconstructed e+e- opening angle (OA) and M(ee) of the signal classes, 'pre', V3 (s5_acc.C part 8). SIMULATION ONLY.\n");
    P("# fractions are weighted (rows owned in mode A, w x V3; X17 1 per row); errors ~ sqrt(f(1-f)/N_eff)\n\n");
    std::vector<std::vector<std::array<double, 3>>> rows(kNC);
    for (int i = 0; i < kNC; ++i) rows[i] = Collect(kC[i].cls, pre, oM, [](const Row& r) { return (double)r.oa; }, [](const Row& r) { return 1e3 * r.mee; });
    auto frac = [&](int i, const std::function<bool(double, double)>& f, double& err) {
        double s = 0, t = 0, t2 = 0;
        for (auto& a : rows[i]) { t += a[2]; t2 += a[2] * a[2]; if (f(a[0], a[1])) s += a[2]; }
        const double p = t > 0 ? s / t : 0., ne = t2 > 0 ? t * t / t2 : 0.;
        err = ne > 0 ? std::sqrt(p * (1 - p) / ne) : 0.;
        return p;
    };
    P("## 1. fraction passing OA > c\n");
    P("c [deg]  ");
    for (int i = 0; i < kNC; ++i) P("| %-22s", kC[i].key);
    P("\n");
    const std::vector<double> cuts = {0.5, 1, 1.5, 2, 3, 4, 6, 9, 12, 20};
    for (double cc : cuts) {
        P("%6.1f   ", cc);
        for (int i = 0; i < kNC; ++i) { double e; const double f = frac(i, [cc](double oa, double) { return oa > cc; }, e); P("| %7.4f +- %6.4f      ", f, e); }
        P("\n");
    }
    P("\n## 2. fraction passing M(ee)_reco > c\n");
    P("c [MeV]  ");
    for (int i = 0; i < kNC; ++i) P("| %-22s", kC[i].key);
    P("\n");
    for (double cc : {2., 3., 5., 8., 10., 15., 20., 30.}) {
        P("%6.1f   ", cc);
        for (int i = 0; i < kNC; ++i) { double e; const double f = frac(i, [cc](double, double m) { return m > cc; }, e); P("| %7.4f +- %6.4f      ", f, e); }
        P("\n");
    }
    // OA quantiles
    P("\n## 3. weighted quantiles of the reconstructed OA [deg] and M(ee) [MeV/c2]\n");
    P("class            N_eff     OA: q05   q16   median   q84   q95   | M(ee): q05   median   q95\n");
    for (int i = 0; i < kNC; ++i) {
        std::vector<std::pair<double, double>> vo, vm;
        for (auto& a : rows[i]) { vo.push_back({a[0], a[2]}); vm.push_back({a[1], a[2]}); }
        const Quant qo = Quantiles(vo), qm = Quantiles(vm);
        P("%-15s %8.0f   %6.2f %6.2f %6.2f %6.2f %6.2f   | %7.2f %7.2f %7.2f\n", kC[i].key, NeffOf(qo), qo.q05, qo.q16, qo.med, qo.q84, qo.q95, qm.q05, qm.med,
          qm.q95);
    }
    // OA resolution: reconstructed minus true OA
    P("\n## 3b. OA resolution: Delta OA = OA - OA_true [deg] (weighted quantiles); the true OA of conversion pairs is small, their\n");
    P("#     reconstructed OA is set by the resolution and by the two-track separation\n");
    P("class            N_eff     true OA median   Delta OA: q16    median    q84    (q84-q16)/2\n");
    for (int i = 0; i < kNC; ++i) {
        std::vector<std::pair<double, double>> vt, vd;
        for (auto& a : Collect(kC[i].cls, pre, oM, [](const Row& r) { return (double)r.oa_true; }, [](const Row& r) { return (double)(r.oa - r.oa_true); })) {
            vt.push_back({a[0], a[2]}); vd.push_back({a[1], a[2]});
        }
        const Quant qt = Quantiles(vt), qd = Quantiles(vd);
        P("%-15s %8.0f   %8.3f          %7.3f   %7.3f   %7.3f   %7.3f\n", kC[i].key, NeffOf(qd), qt.med, qd.q16, qd.med, qd.q84, 0.5 * (qd.q84 - qd.q16));
    }
    // separation direct vs conversion: OA cut and M(ee) cut, and the 2D cut OA > a or M(ee) > b
    P("\n## 4. direct (5) against conversion peak (4): kept fractions for OA > c and for M(ee) > c (ROC points); bump Dalitz (2) against bump\n");
    P("#    conversion (3) with the same cuts\n");
    P("cut              eps(5)    eps(4)    eps(4)/eps(5)   | eps(2)    eps(3)\n");
    const int i5 = CIdx(5), i4 = CIdx(4), i2 = CIdx(2), i3 = CIdx(3), i6 = CIdx(6);
    for (double cc : cuts) {
        double e;
        const double a5 = frac(i5, [cc](double oa, double) { return oa > cc; }, e), a4 = frac(i4, [cc](double oa, double) { return oa > cc; }, e);
        const double a2 = frac(i2, [cc](double oa, double) { return oa > cc; }, e), a3 = frac(i3, [cc](double oa, double) { return oa > cc; }, e);
        P("OA > %4.1f deg     %7.4f   %7.4f   %7.4f         | %7.4f   %7.4f\n", cc, a5, a4, a4 / a5, a2, a3);
    }
    for (double cc : {3., 5., 8., 10., 15., 20.}) {
        double e;
        const double a5 = frac(i5, [cc](double, double m) { return m > cc; }, e), a4 = frac(i4, [cc](double, double m) { return m > cc; }, e);
        const double a2 = frac(i2, [cc](double, double m) { return m > cc; }, e), a3 = frac(i3, [cc](double, double m) { return m > cc; }, e);
        P("M(ee) > %4.1f MeV  %7.4f   %7.4f   %7.4f         | %7.4f   %7.4f\n", cc, a5, a4, a4 / a5, a2, a3);
    }
    // X17 against the direct decay in the X17 M(ee) band
    P("\n## 5. X17 (6) against the direct decay (5) at the same reconstructed M(ee) (15.5-18.5 MeV/c2): OA does not separate them\n");
    P("cut          eps(5 | band)   eps(6 | band)\n");
    for (double cc : {2., 3., 4., 5., 6., 8.}) {
        double e;
        const double a5 = frac(i5, [cc](double oa, double m) { return m > 15.5 && m < 18.5 ? oa > cc : false; }, e) /
                          std::max(1e-12, frac(i5, [](double, double m) { return m > 15.5 && m < 18.5; }, e));
        const double a6 = frac(i6, [cc](double oa, double m) { return m > 15.5 && m < 18.5 ? oa > cc : false; }, e) /
                          std::max(1e-12, frac(i6, [](double, double m) { return m > 15.5 && m < 18.5; }, e));
        P("OA > %4.1f     %7.4f          %7.4f\n", cc, a5, a6);
    }
    {
        double e;
        const double b5 = frac(i5, [](double, double m) { return m > 15.5 && m < 18.5; }, e), b6 = frac(i6, [](double, double m) { return m > 15.5 && m < 18.5; }, e);
        P("share of the class inside the band 15.5-18.5 MeV/c2: direct %.4f, X17 %.4f\n", b5, b6);
        if (c.x17hi > c.x17lo) {
            const double lo = c.x17lo * 1e-3, hi = c.x17hi * 1e-3;
            const Sel inW = [pre, lo, hi](const Row& r) { return pre(r) && r.mee >= lo && r.mee < hi; };
            const std::vector<Comp> cm = {CompOf(5), CompX17(kBrX17Ill), CompOf(4), CompOf(2)};
            const std::vector<WSum> sw = SumMany(c.sig, inW, cm, oM);
            P("accepted cross sections ('pre', all M4) with the reconstructed M(ee) in the X17 90 %% window %.2f-%.2f MeV/c2 [nb]:\n", c.x17lo, c.x17hi);
            P("  direct (5) %.4f +- %.4f; X17 (6, BR = 1e-5) %.5f +- %.5f; peak conv. (4) %.4f +- %.4f; bump Dalitz (2) %.4f +- %.4f\n", sw[0].S / kLsim,
              std::sqrt(sw[0].V) / kLsim, sw[1].S / kLsim, std::sqrt(sw[1].V) / kLsim, sw[2].S / kLsim, std::sqrt(sw[2].V) / kLsim, sw[3].S / kLsim,
              std::sqrt(sw[3].V) / kLsim);
            const Sel inWP = [inW](const Row& r) { return inW(r) && r.m4 >= kPeakLo && r.m4 < kPeakHi; };
            const std::vector<WSum> sp = SumMany(c.sig, inWP, cm, oM);
            P("  the same inside the PEAK window: direct %.4f +- %.4f; X17 (1e-5) %.5f +- %.5f; peak conv. %.4f +- %.4f\n", sp[0].S / kLsim,
              std::sqrt(sp[0].V) / kLsim, sp[1].S / kLsim, std::sqrt(sp[1].V) / kLsim, sp[2].S / kLsim, std::sqrt(sp[2].V) / kLsim);
            c.summary.push_back(Form("X17 90%% window %.1f-%.1f MeV, PEAK window: direct %.4f nb, X17(1e-5) %.5f nb, conv. %.4f nb", c.x17lo, c.x17hi,
                                     sp[0].S / kLsim, sp[1].S / kLsim, sp[2].S / kLsim));
        }
        double e2;
        const double f5 = frac(i5, [](double oa, double) { return oa > 2; }, e2), f4 = frac(i4, [](double oa, double) { return oa > 2; }, e2),
                     f6 = frac(i6, [](double oa, double) { return oa > 2; }, e2), f2 = frac(i2, [](double oa, double) { return oa > 2; }, e2),
                     f3 = frac(i3, [](double oa, double) { return oa > 2; }, e2);
        c.summary.push_back(Form("OA > 2 deg keeps: direct %.3f, conv. peak %.3f, X17 %.3f, bump Dalitz %.3f, bump conv. %.3f", f5, f4, f6, f2, f3));
    }
    CloseTab();

    // ---------------- figures: 2D OA vs M(ee)
    {
        TCanvas* cv = new TCanvas(Nm("coa2").c_str(), "", 1700, 1400);
        cv->Divide(2, 2, 0.003, 0.003);
        const std::vector<double> xe = LogEdges(1.0, 300., 60), ye = LogEdges(0.1, 120., 60);
        const int order[4] = {i5, i4, i6, i2};
        for (int k = 0; k < 4; ++k) {
            const int i = order[k];
            cv->cd(k + 1); PadM(0.14, 0.19, 0.14, 0.09); gPad->SetLogx(); gPad->SetLogy(); gPad->SetLogz();
            const bool coarse = kC[i].cls == 2 || kC[i].cls == 4 || kC[i].cls == 3;
            const std::vector<double> xe2 = coarse ? LogEdges(1.0, 300., 25) : xe, ye2 = coarse ? LogEdges(0.1, 120., 25) : ye;
            TH2D* h = new TH2D(Nm("h2").c_str(), "", (int)xe2.size() - 1, xe2.data(), (int)ye2.size() - 1, ye2.data());
            h->SetDirectory(nullptr); h->Sumw2();
            double t = 0;
            for (auto& a : rows[i]) { h->Fill(std::max(1.0001, a[1]), std::max(0.10001, a[0]), a[2]); t += a[2]; }
            if (t > 0) h->Scale(1. / t);
            AxisT(h, "M(e^{+}e^{-}) reconstructed [MeV/c^{2}]", "OA [deg]");
            h->GetZaxis()->SetTitle("fraction per cell"); h->GetZaxis()->SetTitleSize(0.05); h->GetZaxis()->SetLabelSize(0.045);
            h->GetZaxis()->SetTitleOffset(1.35);
            h->SetMinimum(coarse ? 1e-4 : 1e-5); h->SetMaximum(coarse ? 0.3 : 0.1);
            h->Draw("colz");
            TLine l2; l2.SetLineColor(632); l2.SetLineWidth(3); l2.SetLineStyle(2); l2.DrawLine(1, 2, 300, 2);
            TLine l17; l17.SetLineColor(635); l17.SetLineWidth(2); l17.SetLineStyle(3); l17.DrawLine(17, 0.1, 17, 120);
            double e;
            const double f = frac(i, [](double oa, double) { return oa > 2; }, e);
            Txt(0.14, 0.93, Form("(%c) %s, pre%s", 'a' + k, kC[i].key, coarse ? " (coarse bins)" : ""), 0.055, kC[i].col);
            Txt(0.17, 0.82, Form("OA > 2#circ keeps %.3f", f), 0.048);
            Txt(0.17, 0.17, "dashed: OA = 2#circ, dotted: 17 MeV/c^{2}", 0.042);
        }
        cv->cd(0);
        SaveFig(cv, "s5_acc_oa2d");
    }
    // 1D OA, cumulative, ROC
    {
        TCanvas* cv = new TCanvas(Nm("coa1").c_str(), "", 1800, 700);
        cv->Divide(3, 1, 0.003, 0.003);
        const std::vector<double> oe = LogEdges(0.1, 120., 45), oe1 = LogEdges(0.1, 120., 25);
        cv->cd(1); PadM(0.16, 0.03, 0.15, 0.10); gPad->SetLogx(); gPad->SetLogy();
        TH1D* fa = H1v("foa", oe1);
        AxisT(fa, "OA [deg]", "fraction per bin (0.12 decades)");
        fa->SetMinimum(1e-4); fa->SetMaximum(100); fa->Draw("axis");
        TLegend* l = Leg(0.18, 0.64, 0.70, 0.89, 0.042);
        for (int i = 0; i < kNC; ++i) {
            TH1D* h = H1v("hoa", oe1); double t = 0;
            for (auto& a : rows[i]) { h->Fill(std::max(0.10001, a[0]), a[2]); t += a[2]; }
            if (t > 0) h->Scale(1. / t);
            Sty(h, kC[i].col, kC[i].cls == 5 ? 4 : 3, kC[i].cls == 6 ? 2 : 1); h->Draw("hist same");
            l->AddEntry(h, kC[i].key, "l");
        }
        l->Draw();
        TLine l2; l2.SetLineColor(922); l2.SetLineWidth(2); l2.SetLineStyle(2); l2.DrawLine(2, 1e-4, 2, 0.6);
        Txt(0.16, 0.93, "(a) reconstructed OA, pre", 0.055);
        cv->cd(2); PadM(0.16, 0.03, 0.15, 0.10); gPad->SetLogx();
        TH1D* fb = H1v("foab", oe);
        AxisT(fb, "OA cut c [deg]", "fraction with OA > c");
        fb->SetMinimum(0); fb->SetMaximum(1.05); fb->Draw("axis");
        for (int i = 0; i < kNC; ++i) {
            TGraph* g = new TGraph();
            for (int k = 0; k < (int)oe.size(); ++k) {
                double e; const double cc = oe[k];
                g->SetPoint(k, cc, frac(i, [cc](double oa, double) { return oa > cc; }, e));
            }
            Sty(g, kC[i].col, 1, 0.1, kC[i].cls == 5 ? 4 : 3, kC[i].cls == 6 ? 2 : 1); g->Draw("l same");
        }
        TLine l2b; l2b.SetLineColor(922); l2b.SetLineWidth(2); l2b.SetLineStyle(2); l2b.DrawLine(2, 0, 2, 1.05);
        Txt(0.16, 0.93, "(b) fraction passing OA > c", 0.055);
        cv->cd(3); PadM(0.16, 0.03, 0.15, 0.10);
        TH1D* fc = H1("foac", 10, 0, 1.02);
        AxisT(fc, "kept fraction of the direct decay (5)", "kept fraction of the background class");
        fc->SetMinimum(0); fc->SetMaximum(1.02); fc->Draw("axis");
        TGraph* r4 = new TGraph(); TGraph* r4m = new TGraph(); TGraph* r3 = new TGraph();
        for (int k = 0; k < (int)oe.size(); ++k) {
            double e; const double cc = oe[k];
            r4->SetPoint(k, frac(i5, [cc](double oa, double) { return oa > cc; }, e), frac(i4, [cc](double oa, double) { return oa > cc; }, e));
            r3->SetPoint(k, frac(i2, [cc](double oa, double) { return oa > cc; }, e), frac(i3, [cc](double oa, double) { return oa > cc; }, e));
        }
        const std::vector<double> me = LogEdges(1.0, 100., 40);
        for (int k = 0; k < (int)me.size(); ++k) {
            double e; const double cc = me[k];
            r4m->SetPoint(k, frac(i5, [cc](double, double m) { return m > cc; }, e), frac(i4, [cc](double, double m) { return m > cc; }, e));
        }
        Sty(r4, 618, 1, 0.1, 4, 1); Sty(r4m, 618, 1, 0.1, 4, 2); Sty(r3, 902, 1, 0.1, 4, 1);
        r4->Draw("l same"); r4m->Draw("l same"); r3->Draw("l same");
        // OA > 2 points
        double e;
        TMarker m4(frac(i5, [](double oa, double) { return oa > 2; }, e), frac(i4, [](double oa, double) { return oa > 2; }, e), 22);
        m4.SetMarkerColor(618); m4.SetMarkerSize(2.2); m4.DrawClone();
        TMarker m3(frac(i2, [](double oa, double) { return oa > 2; }, e), frac(i3, [](double oa, double) { return oa > 2; }, e), 25);
        m3.SetMarkerColor(902); m3.SetMarkerSize(2.2); m3.DrawClone();
        TLine dg; dg.SetLineStyle(3); dg.SetLineColor(922); dg.DrawLine(0, 0, 1.02, 1.02);
        TLegend* lc = Leg(0.18, 0.62, 0.70, 0.89, 0.040);
        lc->AddEntry(r4, "conv. peak (4) vs direct (5), OA cut", "l");
        lc->AddEntry(r4m, "the same, M(e^{+}e^{-}) cut", "l");
        lc->AddEntry(r3, "bump conv. (3) vs bump Dalitz (2), OA", "l");
        lc->Draw();
        Txt(0.50, 0.20, "markers: OA > 2#circ", 0.045);
        Txt(0.16, 0.93, "(c) separation: kept fractions", 0.055);
        cv->cd(0);
        SaveFig(cv, "s5_acc_oa");
    }
    printf("PartOA done in %.0f s\n", sw.RealTime());
}

}  // namespace s5acc

void s5_acc(int parts = 15) {
    using namespace s5acc;
    TStopwatch sw;
    gROOT->SetBatch(kTRUE);
    if (!simphys::Mdl().ok && !simphys::Init(simphys::kRwGenFile)) return;
    if (parts & 1) PartGen();
    if (parts & 14) SetupProduction();
    if (parts & 2) PartAcc();
    if (parts & 4) PartRes();
    if (parts & 8) PartOA();
    printf("\n==== s5_acc summary ====\n");
    for (auto& s : C().summary) printf("%s\n", s.c_str());
    printf("angles: hit %lld, miss %lld\n", simphys::AngHit(), simphys::AngMiss());
    printf("s5_acc done in %.0f s\n", sw.RealTime());
}
