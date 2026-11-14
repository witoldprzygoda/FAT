// ============================================================================
// sim_datacomp.C  --  data against simulation: absolute scale and shape
// ============================================================================
// HADES p+p, T = 4.53 GeV, pi+ pi- e+ e-. Simulation (SMASH std + lep, owned
// rows in partition mode A, weight w = PDG BR) as a cross section within the
// HADES acceptance after the selection, without efficiency correction:
//     dsigma/dx = sum w / L_sim / bin width,  L_sim = 221.05 nb^-1.
// Simulated spectrum = signal (template 1) + physical background (template 2);
// fake lepton (3) and different sources (4) only as separate curves.
// Data in the same units: (N - CB) / L_data / bin width, CB = 2 sqrt(N++ N--)
// (M4: 10 MeV CB bins as in the data study; other variables: CB per bin of the
// variable), L_data = 6080.5 nb^-1 (producers' note, UNCONFIRMED).
// No trigger-bias or efficiency correction on either side.
//
// Statistics
//   simulation: per bin sum w and event-group variance (SimEngine FillMany)
//   data:       per bin rows and EVENT-GROUP variance (sum over groups of the
//               squared number of the group's rows in the bin; this file,
//               dcmp::FillData). ScanEngine's per-bin errors are Poisson; the
//               group variance matters for the pair variables of one lepton
//               pair combined with several pions. CB variance from the
//               like-sign group variances (linear propagation). Window sums of
//               data: ScanEngine windows (group variance, 10 MeV CB bins).
//   single-lepton variables of like-sign rows: both leptons enter with weight
//               1/2 (the like-sign "ep"/"em" slots are not charge slots)
//
// Parts (argument of the macro; 0 = all):
//   1  M4 spectra for the presets pre, pre+OA>2, pre+pions in vertex fit,
//      slides: absolute comparison, renormalisation R (0.60-0.80; variants
//      0.56-0.62, 0.65-0.85; mode B and w_db), shape after R with chi2/ndf in
//      the regions outside the eta structures; components of the physical
//      background (pion-pair classes); check of the data filler against
//      ScanEngine
//   2  OA(e+e-) and M(e+e-) classes in M4 0.56-0.80 and 0.30-0.45: ratios
//      (absolute, after R), distributions, shared-RICH-ring share; content of
//      the data CB against the simulated information classes fake lepton and
//      different sources
//   3  subsystems of the pi+pi-e+e- hypothesis (pi+pi- pair, e+e- pair, single
//      tracks, oa_lab, mm4, pt4) for all M4, M4 0.56-0.80 and the bump window;
//      presets pre and pre+OA>2
//
// Figures: plots/eta_sim/figs/datacomp_*.png (.pdf for the main ones)
// Tables:  plots/eta_sim/tables/datacomp_*.txt
//
// Run from /home/przygoda/HADES/FAT/FAT (compile once, then the parts may run
// in parallel):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/aclic_datacomp", kTRUE)' 'plots/eta_sim/sim_datacomp.C+(1)'
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

#include "TBox.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TPad.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TSystem.h"

#include "SimEngine.h"

namespace dcmp {
using namespace simeng;

// ============================================================================
// 1. Data histograms with event-group variance
// ============================================================================
// entries of one row: up to two (value, weight) pairs; k = sample (0 us, 1 pp, 2 mm)
typedef std::function<int(const scan::Row&, int, double*, double*)> DFill;
inline DFill DVar(const std::function<double(const scan::Row&)>& v) {
    return [v](const scan::Row& r, int, double* x, double* w) { x[0] = v(r); w[0] = 1; return 1; };
}
// single-lepton variable: unlike-sign rows use slot lep (0 = e+ "ep", 1 = e- "em");
// like-sign rows give both leptons with weight 1/2
inline DFill DLep(const std::function<double(const scan::Row&, int)>& v, int lep) {
    return [v, lep](const scan::Row& r, int k, double* x, double* w) {
        if (k == 0) { x[0] = v(r, lep); w[0] = 1; return 1; }
        x[0] = v(r, 0); x[1] = v(r, 1); w[0] = w[1] = 0.5; return 2;
    };
}

struct DH {
    int nbin = 0, reb = 1, nw = 0;
    double lo = 0, hi = 0, bw = 0, Ldata = kLdata;
    std::vector<double> n[3], v[3], wn[3], wv[3];
    std::vector<double> cbW, cbWVar, cb, cbVar, sig, sigVar;
    double x(int i) const { return lo + (i + 0.5) * bw; }
    void Interp(int i, int& j0, int& j1, double& f) const {   // as scan::Spectra::Interp
        const double t = (x(i) - (lo + 0.5 * reb * bw)) / (reb * bw);
        if (t <= 0) { j0 = j1 = 0; f = 0; return; }
        if (t >= nw - 1) { j0 = j1 = nw - 1; f = 0; return; }
        j0 = (int)t; j1 = j0 + 1; f = t - j0;
    }
    // what: 0 N, 1 N++, 2 N--, 3 CB, 4 N - CB
    double val(int what, int i) const { return what < 3 ? n[what][i] : what == 3 ? cb[i] : sig[i]; }
    double var(int what, int i) const { return what < 3 ? v[what][i] : what == 3 ? cbVar[i] : sigVar[i]; }
    // unit 0: rows; 2: nb per unit of x (rows / L_data / bin width)
    TH1D* Hist(int what, const char* name, int unit = 2) const {
        TH1D* h = new TH1D(name, "", nbin, lo, hi);
        h->SetDirectory(nullptr);
        h->Sumw2();
        const double f = unit == 0 ? 1. : 1. / (Ldata * bw);
        for (int i = 0; i < nbin; ++i) {
            h->SetBinContent(i + 1, val(what, i) * f);
            h->SetBinError(i + 1, std::sqrt(std::max(0., var(what, i))) * f);
        }
        return h;
    }
};

inline DH FillData(const scan::Data& d, const scan::Sel& sel, int nbin, double lo, double hi, int reb,
                   const DFill& fx = DFill(), double Ldata = kLdata) {
    DH h;
    h.nbin = nbin; h.lo = lo; h.hi = hi; h.bw = (hi - lo) / nbin; h.Ldata = Ldata;
    h.reb = (reb > 0 && nbin % reb == 0) ? reb : 1;
    h.nw = nbin / h.reb;
    for (int k = 0; k < 3; ++k) {
        h.n[k].assign(nbin, 0.); h.v[k].assign(nbin, 0.);
        h.wn[k].assign(h.nw, 0.); h.wv[k].assign(h.nw, 0.);
    }
    std::vector<std::pair<int, double>> buf;
    for (int k = 0; k < 3; ++k) {
        auto flush = [&]() {
            if (buf.empty()) return;
            std::sort(buf.begin(), buf.end());
            size_t i = 0;
            while (i < buf.size()) {
                const int b = buf[i].first; double a = 0;
                while (i < buf.size() && buf[i].first == b) { a += buf[i].second; ++i; }
                h.v[k][b] += a * a;
            }
            i = 0;
            while (i < buf.size()) {
                const int j = buf[i].first / h.reb; double a = 0;
                while (i < buf.size() && buf[i].first / h.reb == j) { a += buf[i].second; ++i; }
                h.wv[k][j] += a * a;
            }
            buf.clear();
        };
        int cur = -1;
        for (const scan::Row& r : d.s[k]) {
            if (!r.cor_ok && !d.fallback) continue;
            if (r.grp != cur) { flush(); cur = r.grp; }
            if (!sel(r)) continue;
            double xs[2] = {0, 0}, ws[2] = {0, 0};
            int ne = 1;
            if (fx) ne = fx(r, k, xs, ws);
            else { xs[0] = r.m4; ws[0] = 1; }
            for (int e = 0; e < ne; ++e) {
                const double m = xs[e];
                if (!(m == m) || m < lo || m >= hi) continue;
                const int ib = (int)((m - lo) / h.bw);
                if (ib < 0 || ib >= nbin) continue;
                h.n[k][ib] += ws[e];
                buf.push_back(std::make_pair(ib, ws[e]));
            }
        }
        flush();
    }
    for (int k = 0; k < 3; ++k)
        for (int i = 0; i < nbin; ++i) h.wn[k][i / h.reb] += h.n[k][i];
    h.cbW.assign(h.nw, 0.); h.cbWVar.assign(h.nw, 0.);
    for (int j = 0; j < h.nw; ++j) {
        const double a = h.wn[1][j], b = h.wn[2][j];
        h.cbW[j] = 2.0 * std::sqrt(a * b);
        // dCB/dN++ = sqrt(N--/N++), dCB/dN-- = sqrt(N++/N--); group variances of N++, N--
        h.cbWVar[j] = (a > 0 && b > 0) ? (b / a) * h.wv[1][j] + (a / b) * h.wv[2][j] : h.wv[1][j] + h.wv[2][j];
    }
    h.cb.assign(nbin, 0.); h.cbVar.assign(nbin, 0.); h.sig.assign(nbin, 0.); h.sigVar.assign(nbin, 0.);
    for (int i = 0; i < nbin; ++i) {
        int j0, j1; double f; h.Interp(i, j0, j1, f);
        const double a0 = (1 - f) / h.reb, a1 = f / h.reb;
        h.cb[i] = a0 * h.cbW[j0] + a1 * h.cbW[j1];
        h.cbVar[i] = a0 * a0 * h.cbWVar[j0] + a1 * a1 * h.cbWVar[j1];
        h.sig[i] = h.n[0][i] - h.cb[i];
        h.sigVar[i] = h.v[0][i] + h.cbVar[i];
    }
    return h;
}

// ============================================================================
// 2. Ratios, chi2, window sums
// ============================================================================
inline TH1D* Ratio(const TH1D* d, const TH1D* s, const char* name, double scaleSim = 1.0) {
    TH1D* r = (TH1D*)d->Clone(name);
    r->SetDirectory(nullptr);
    r->Reset();
    for (int i = 1; i <= d->GetNbinsX(); ++i) {
        const double sv = s->GetBinContent(i) * scaleSim, se = s->GetBinError(i) * scaleSim;
        const double dv = d->GetBinContent(i), de = d->GetBinError(i);
        if (sv > 0) {
            const double q = dv / sv;
            r->SetBinContent(i, q);
            r->SetBinError(i, std::sqrt(de * de + q * q * se * se) / sv);
        }
    }
    return r;
}

// chi2 of the shape: bins whose centre is in [a, b), with data rows >= minN and
// simulated rows >= minS. Fixed factor R: ndf = n. Free factor: ndf = n - 1.
struct Chi {
    double chiR = 0, chiF = 0, rF = 0, rFe = 0;
    int nR = 0, nF = 0;
};
inline Chi ShapeChi(const TH1D* d, const TH1D* s, const std::vector<double>& dRows, const std::vector<double>& sRows,
                    double a, double b, double R, double minN = 10, double minS = 5) {
    Chi c;
    std::vector<double> dv, de2, sv, se2;
    for (int i = 1; i <= d->GetNbinsX(); ++i) {
        const double xc = d->GetBinCenter(i);
        if (xc < a || xc >= b) continue;
        if (dRows[i - 1] < minN || sRows[i - 1] < minS) continue;
        dv.push_back(d->GetBinContent(i)); de2.push_back(std::pow(d->GetBinError(i), 2));
        sv.push_back(s->GetBinContent(i)); se2.push_back(std::pow(s->GetBinError(i), 2));
    }
    const int n = (int)dv.size();
    for (int i = 0; i < n; ++i) {
        const double v = de2[i] + R * R * se2[i];
        if (v > 0) { c.chiR += std::pow(dv[i] - R * sv[i], 2) / v; ++c.nR; }
    }
    double r = R, den = 0;
    for (int it = 0; it < 30; ++it) {
        double num = 0; den = 0;
        for (int i = 0; i < n; ++i) {
            const double v = de2[i] + r * r * se2[i];
            if (v <= 0) continue;
            num += dv[i] * sv[i] / v; den += sv[i] * sv[i] / v;
        }
        if (den <= 0) break;
        const double rn = num / den;
        if (std::fabs(rn - r) < 1e-9) { r = rn; break; }
        r = rn;
    }
    c.rF = r; c.rFe = den > 0 ? 1. / std::sqrt(den) : 0.;
    int m = 0;
    for (int i = 0; i < n; ++i) {
        const double v = de2[i] + r * r * se2[i];
        if (v > 0) { c.chiF += std::pow(dv[i] - r * sv[i], 2) / v; ++m; }
    }
    c.nF = std::max(0, m - 1);
    return c;
}

// data window sum (N - CB) in nb and the simulated sum, with errors
struct WinCmp {
    double a = 0, b = 0, N = 0, CB = 0, ds = 0, de = 0, ss = 0, se = 0, sNeff = 0;
    double r() const { return ss > 0 ? ds / ss : 0; }
    double re() const { return (ss > 0 && ds != 0) ? r() * std::sqrt(std::pow(de / ds, 2) + std::pow(se / ss, 2)) : 0; }
};
inline Sel SelM4(double a, double b) { return [a, b](const Row& r) { return r.m4 >= a && r.m4 < b; }; }
inline scan::Sel DSelM4(double a, double b) { return [a, b](const scan::Row& r) { return r.m4 >= a && r.m4 < b; }; }

// ============================================================================
// 3. Plot helpers
// ============================================================================
int gCnt = 0;
inline const char* Nm(const char* base) { return Form("%s_%d", base, gCnt++); }
inline void Txt(double x, double y, const char* s, double size, int col = 1, int align = 11) {
    TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(size); t.SetTextColor(col); t.SetTextAlign(align);
    t.DrawLatex(x, y, s);
}
inline void StyleL(TH1D* h, int col, int lw, int ls, int fill = 0) {
    h->SetLineColor(col); h->SetLineWidth(lw); h->SetLineStyle(ls); h->SetMarkerSize(0); h->SetMarkerColor(col);
    if (fill) { h->SetFillColor(fill); h->SetFillStyle(1001); } else h->SetFillStyle(0);
}
inline void StyleM(TH1D* h, int col, int st, double sz) {
    h->SetMarkerStyle(st); h->SetMarkerColor(col); h->SetMarkerSize(sz); h->SetLineColor(col); h->SetLineWidth(2);
}
struct PP { TPad* up; TPad* lo; double hu, hl; };
// a cell with an upper pad (logy) and a ratio pad below; hu, hl: heights in canvas units
inline PP MkPair(TCanvas* c, double x1, double y1, double x2, double y2, double frac, bool logy, double lm = 0.16,
                 double tm = 0.10) {
    PP p;
    const double ym = y1 + frac * (y2 - y1);
    c->cd();
    p.up = new TPad(Nm("pu"), "", x1, ym, x2, y2);
    p.up->SetLeftMargin(lm); p.up->SetRightMargin(0.03); p.up->SetTopMargin(tm); p.up->SetBottomMargin(0.015);
    p.up->SetLogy(logy); p.up->Draw();
    c->cd();
    p.lo = new TPad(Nm("pl"), "", x1, y1, x2, ym);
    p.lo->SetLeftMargin(lm); p.lo->SetRightMargin(0.03); p.lo->SetTopMargin(0.02); p.lo->SetBottomMargin(0.36);
    p.lo->Draw();
    p.hu = y2 - ym; p.hl = ym - y1;
    return p;
}
// frame histogram for a pad; sz = label size in canvas units
inline TH1D* Frame(const TH1D* like, double ymin, double ymax, const char* xt, const char* yt, double sz, double h,
                   bool xlab, double yoff = 1.3) {
    TH1D* f = (TH1D*)like->Clone(Nm("fr"));
    f->SetDirectory(nullptr); f->Reset();
    f->SetMinimum(ymin); f->SetMaximum(ymax);
    f->GetXaxis()->SetTitle(xlab ? xt : "");
    f->GetYaxis()->SetTitle(yt);
    f->GetXaxis()->SetLabelSize(xlab ? sz / h : 0); f->GetXaxis()->SetTitleSize(xlab ? 1.1 * sz / h : 0);
    f->GetYaxis()->SetLabelSize(sz / h); f->GetYaxis()->SetTitleSize(1.1 * sz / h);
    f->GetXaxis()->SetTitleOffset(1.05); f->GetYaxis()->SetTitleOffset(yoff);
    f->GetYaxis()->SetNdivisions(505); f->GetXaxis()->SetTickLength(0.03 * 0.5 / h);
    f->GetYaxis()->SetMaxDigits(3);
    return f;
}
inline double MaxOf(const std::vector<TH1D*>& hs, bool withErr) {
    double m = 0;
    for (TH1D* h : hs)
        for (int i = 1; i <= h->GetNbinsX(); ++i) m = std::max(m, h->GetBinContent(i) + (withErr ? h->GetBinError(i) : 0.));
    return m;
}
inline double MinPos(const std::vector<TH1D*>& hs) {
    double m = 1e300;
    for (TH1D* h : hs)
        for (int i = 1; i <= h->GetNbinsX(); ++i) if (h->GetBinContent(i) > 0) m = std::min(m, h->GetBinContent(i));
    return m;
}
inline void HLine(double x1, double x2, double y, int col, int ls, int lw = 2) {
    TLine* l = new TLine(x1, y, x2, y); l->SetLineColor(col); l->SetLineStyle(ls); l->SetLineWidth(lw); l->Draw();
}
inline void VLine(double x, double y1, double y2, int col, int ls, int lw = 1) {
    TLine* l = new TLine(x, y1, x, y2); l->SetLineColor(col); l->SetLineStyle(ls); l->SetLineWidth(lw); l->Draw();
}
inline void Band(double x1, double x2, double y1, double y2, int col) {
    TBox* b = new TBox(x1, y1, x2, y2); b->SetFillColorAlpha(col, 0.25); b->SetLineWidth(0); b->Draw();
}
inline void Save(TCanvas* c, const std::string& base, bool pdf) {
    c->SaveAs((std::string(kFigDir) + "/" + base + ".png").c_str());
    if (pdf) c->SaveAs((std::string(kFigDir) + "/" + base + ".pdf").c_str());
}
inline FILE* OpenTab(const char* name) {
    FILE* f = fopen((std::string(kTabDir) + "/" + name).c_str(), "w");
    if (!f) printf("cannot write %s\n", name);
    return f;
}
// print to stdout and to a table file
inline void P2(FILE* f, const char* fmt, ...) __attribute__((format(printf, 2, 3)));
inline void P2(FILE* f, const char* fmt, ...) {
    va_list a; va_start(a, fmt); vprintf(fmt, a); va_end(a);
    if (f) { va_list b; va_start(b, fmt); vfprintf(f, fmt, b); va_end(b); }
}

const char* kYM4 = "d#sigma/dM_{#pi^{+}#pi^{-}e^{+}e^{-}} [nb/(GeV/c^{2})]";
const char* kXM4 = "M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}]";
const char* kDirect = "#eta #rightarrow #pi^{+}#pi^{-}e^{+}e^{-} (direct decay) is NOT in the simulation (eta set not available)";
inline std::string DataLine() {
    return Form("data: (N - CB) / L_{data} / bin width, %s; CB = 2#sqrt{N_{++}N_{--}} (no k factor); errors from event groups", kLdataNote);
}
inline std::string SimLine(const SimOpt& o = SimOpt()) {
    return Form("simulation: signal + physical background (templates 1+2), SMASH %s; #Sigma w / L_{sim}, L_{sim} = %.2f nb^{-1}; no efficiency correction",
                OptLabel(o).c_str(), kLsim);
}
// presets of this study and file tags
const int kNP = 4;
const int kPresets[kNP] = {kPre, kPreOA2, kPrePiVtx, kSlides};
const char* kTag[kNP] = {"pre", "preOA2", "prePiVtx", "slides"};

// renormalisation windows (0 = default)
const int kNRw = 3;
const double kRlo[kNRw] = {0.60, 0.56, 0.65};
const double kRhi[kNRw] = {0.80, 0.62, 0.85};

// data window sums for one selection: N, CB, N - CB (ScanEngine windows: 2.5 MeV
// bins, CB in 10 MeV bins interpolated, group variances), windows 2.. of the list
inline std::vector<scan::WinSum> DataWindows(const scan::Data& dd, const scan::Sel& s,
                                             const std::vector<std::pair<double, double>>& w) {
    scan::FillOpt fo; fo.nbin = 320; fo.lo = 0.30; fo.hi = 1.10; fo.reb = 4; fo.windows = w;
    const scan::Spectra sp = scan::Fill(dd, s, fo);
    return std::vector<scan::WinSum>(sp.win.begin() + 2, sp.win.end());
}
// R in the window [a, b) for a preset and a simulation option
inline WinCmp CompareWindow(const SimData& sd, const scan::WinSum& w, int preset, double a, double b,
                            const SimOpt& o = SimOpt()) {
    WinCmp c; c.a = a; c.b = b; c.N = w.n[0]; c.CB = w.cb;
    c.ds = w.sig() / kLdata; c.de = w.sigErr() / kLdata;
    const WSum s = SumMany(sd, SelAnd(SelPreset(preset), SelM4(a, b)), std::vector<Comp>(1, CompSim()), o)[0];
    c.ss = s.Xsec(); c.se = s.XsecErr(); c.sNeff = s.Neff();
    return c;
}
inline double GetR(const SimData& sd, const scan::Data& dd, int preset, double* err = nullptr, int iw = 0) {
    std::vector<std::pair<double, double>> w(1, std::make_pair(kRlo[iw], kRhi[iw]));
    const std::vector<scan::WinSum> ws = DataWindows(dd, DataSelPreset(preset), w);
    const WinCmp c = CompareWindow(sd, ws[0], preset, kRlo[iw], kRhi[iw]);
    if (err) *err = c.re();
    return c.r();
}

}  // namespace dcmp

// ============================================================================
// Part 1: M4 spectra, renormalisation, shape
// ============================================================================
namespace dcmp {

void Part1(const SimData& sd, const scan::Data& dd) {
    TStopwatch sw;
    // ---- check of the data filler against ScanEngine ----
    {
        FILE* f = OpenTab("datacomp_check.txt");
        P2(f, "# check: dcmp::FillData (event-group variance) against ScanEngine (Poisson per bin), data, selection pre\n");
        struct Cfg { const char* name; int nbin; double lo, hi; int reb; scan::Var x; };
        std::vector<Cfg> cfgs = {
            {"M4 5 MeV, CB 10 MeV", 160, 0.30, 1.10, 2, scan::Var()},
            {"M4 2.5 MeV, CB 10 MeV", 120, 0.40, 0.70, 4, scan::Var()},
            {"M(ee) 0-0.15, 2.5 MeV", 60, 0.0, 0.15, 1, [](const scan::Row& r) { return (double)r.mee; }},
            {"OA 0-30 deg, 0.5 deg", 60, 0.0, 30.0, 1, [](const scan::Row& r) { return (double)r.oa; }},
            {"M(pipi) 0.27-1.07", 40, 0.27, 1.07, 1, [](const scan::Row& r) { return (double)r.mpipi; }}};
        P2(f, "%-24s %12s %12s %12s %14s %14s %14s\n", "binning", "max|dN|", "max|dCB|", "max|dN-CB|", "<V_N grp/N>",
           "<V_CB grp/P>", "<V_sig grp/P>");
        for (const Cfg& c : cfgs) {
            scan::FillOpt fo; fo.nbin = c.nbin; fo.lo = c.lo; fo.hi = c.hi; fo.reb = c.reb; fo.x = c.x;
            const scan::Spectra sp = scan::Fill(dd, DataSelPreset(kPre), fo);
            const DH h = FillData(dd, DataSelPreset(kPre), c.nbin, c.lo, c.hi, c.reb,
                                  c.x ? DVar(c.x) : DFill());
            double dN = 0, dC = 0, dS = 0, rN = 0, rC = 0, rS = 0, wN = 0, wC = 0, wS = 0;
            for (int i = 0; i < c.nbin; ++i) {
                dN = std::max(dN, std::fabs(h.n[0][i] - sp.n[0][i]));
                dC = std::max(dC, std::fabs(h.cb[i] - sp.cb[i]));
                dS = std::max(dS, std::fabs(h.sig[i] - sp.sig[i]));
                if (sp.n[0][i] > 0) { rN += h.v[0][i]; wN += sp.n[0][i]; }
                if (sp.cbVar[i] > 0) { rC += h.cbVar[i]; wC += sp.cbVar[i]; }
                if (sp.sigVar[i] > 0) { rS += h.sigVar[i]; wS += sp.sigVar[i]; }
            }
            P2(f, "%-24s %12.3g %12.3g %12.3g %14.4f %14.4f %14.4f\n", c.name, dN, dC, dS, wN > 0 ? rN / wN : 0,
               wC > 0 ? rC / wC : 0, wS > 0 ? rS / wS : 0);
        }
        P2(f, "# <V grp / Poisson> = summed group variance / summed Poisson variance over the bins (1 = no correlation)\n");
        if (f) fclose(f);
    }

    // ---- components ----
    std::vector<Comp> comps = {CompSim(),       CompPhysBkg(), CompBumpDalitz(), CompBumpConv(),
                               CompPeakConv(),  CompFake(),    CompDiff(),       CompSignal()};
    enum { cSim = 0, cT2, cBD, cBC, cPC, cT3, cT4, cT1, cN };
    const int col[cN] = {kColSim, kColPhysBkg, kColBumpDalitz, kColBumpConv, kColPeakConv, kColFake, kColDiff, 1};

    Binning b5;  b5.nbin = 160;  b5.lo = 0.30; b5.hi = 1.10;
    Binning b10; b10.nbin = 80;  b10.lo = 0.30; b10.hi = 1.10;
    Binning b25; b25.nbin = 320; b25.lo = 0.30; b25.hi = 1.10;
    Binning bz;  bz.nbin = 120;  bz.lo = 0.40;  bz.hi = 0.70;
    Binning bz5; bz5.nbin = 60;  bz5.lo = 0.40; bz5.hi = 0.70;

    // regions of the tables
    struct Reg { const char* name; double a, b; };
    const std::vector<Reg> regs = {{"below the bump", 0.30, 0.45},       {"bump window", scan::kBumpLo, scan::kBumpHi},
                                   {"gap bump-peak", 0.525, 0.530},       {"peak window", scan::kPeakLo, scan::kPeakHi},
                                   {"above the peak", 0.56, 0.90},        {"  0.56-0.70", 0.56, 0.70},
                                   {"  0.70-0.90", 0.70, 0.90},           {"high mass", 0.90, 1.10},
                                   {"R window 0.60-0.80", 0.60, 0.80},    {"R window 0.56-0.62", 0.56, 0.62},
                                   {"R window 0.65-0.85", 0.65, 0.85},    {"all 0.30-1.10", 0.30, 1.10}};
    std::vector<std::pair<double, double>> wl;
    for (const Reg& r : regs) wl.push_back(std::make_pair(r.a, r.b));

    FILE* tReg = OpenTab("datacomp_regions.txt");
    FILE* tChi = OpenTab("datacomp_chi2.txt");
    FILE* tR = OpenTab("datacomp_R.txt");
    P2(tReg, "# data (N - CB) / L_data against the simulation (signal + physical background) in M4 regions\n");
    P2(tReg, "# data: L_data = %.1f nb^-1 (UNCONFIRMED); CB = 2 sqrt(N++ N--) in 10 MeV bins, no k; group variances (ScanEngine windows)\n", kLdata);
    P2(tReg, "# sim: std + lep owned rows, mode A, w (PDG), L_sim = %.2f nb^-1; eta -> pi+pi-e+e- (direct) NOT in the simulation\n", kLsim);
    P2(tReg, "# after R = data / (R sim), R = data/sim in 0.60-0.80 of the same preset; the error of R is not included in that column\n");
    P2(tR, "# renormalisation factor R = data (N - CB) / sim (signal + physical background) in an M4 window\n");
    P2(tR, "# L_data = %.1f nb^-1 (UNCONFIRMED: R scales with 1/L_data); sim std + lep owned rows, L_sim = %.2f nb^-1\n", kLdata, kLsim);
    P2(tR, "%-26s %-12s %14s %14s %14s %14s\n", "preset", "window", "A, w (PDG)", "A, w_db", "B, w", "B, w_db");
    P2(tChi, "# shape after renormalisation: chi2 = sum (d - R s)^2 / (var d + R^2 var s); bins with >= 10 data rows and >= 5 sim rows\n");
    P2(tChi, "# fixed: R from 0.60-0.80 (ndf = bins); free: best factor r inside the region (ndf = bins - 1; pure shape)\n");
    P2(tChi, "# bins: 10 MeV (gap 0.525-0.530: 2.5 MeV); data group variances, CB variance included; sim group variances\n");
    P2(tChi, "%-26s %-16s %12s %10s %10s %16s %10s\n", "preset", "region", "M4", "chi2/ndf R", "p(R)", "free r", "chi2/ndf r");

    std::vector<TH1D*> ratAbs(kNP), ratShape(kNP);
    std::vector<double> Rv(kNP), Re(kNP);
    double yLoLog = 0, yHiLog = 0;
    for (int ip = 0; ip < kNP; ++ip) {
        const int p = kPresets[ip];
        const Sel sel = SelPreset(p);
        const scan::Sel dsel = DataSelPreset(p);
        SimOpt opt;
        const std::vector<WSpec> s5 = FillMany(sd, sel, comps, opt, b5);
        const std::vector<WSpec> s10 = FillMany(sd, sel, comps, opt, b10);
        const std::vector<WSpec> s25 = FillMany(sd, sel, comps, opt, b25);
        const std::vector<WSpec> sz = FillMany(sd, sel, comps, opt, bz);
        const std::vector<WSpec> sz5 = FillMany(sd, sel, comps, opt, bz5);
        const DH d5 = FillData(dd, dsel, 160, 0.30, 1.10, 2);
        const DH d10 = FillData(dd, dsel, 80, 0.30, 1.10, 1);
        const DH d25 = FillData(dd, dsel, 320, 0.30, 1.10, 4);
        const DH dz = FillData(dd, dsel, 120, 0.40, 0.70, 4);
        const DH dz5 = FillData(dd, dsel, 60, 0.40, 0.70, 2);

        // ---- windows: data and simulation (variants) ----
        const std::vector<scan::WinSum> ws = DataWindows(dd, dsel, wl);
        std::vector<WinCmp> wc(regs.size());
        for (size_t r = 0; r < regs.size(); ++r) wc[r] = CompareWindow(sd, ws[r], p, regs[r].a, regs[r].b);
        const double R = wc[8].r(), Rerr = wc[8].re();
        Rv[ip] = R; Re[ip] = Rerr;
        // R variants
        for (int iw = 0; iw < kNRw; ++iw) {
            const size_t ir = 8 + iw;
            std::string line = Form("%-26s %4.2f-%-7.2f", PresetName(p), kRlo[iw], kRhi[iw]);
            for (int v = 0; v < 4; ++v) {
                SimOpt o; o.mode = v / 2; o.db = (v % 2) == 1;
                const WinCmp c = CompareWindow(sd, ws[ir], p, regs[ir].a, regs[ir].b, o);
                line += Form(" %6.3f+-%-6.3f", c.r(), c.re());
            }
            P2(tR, "%s\n", line.c_str());
        }
        // region table with N_eff of every template
        P2(tReg, "\n## %s   (R = %.3f +- %.3f from 0.60-0.80)\n", PresetName(p), R, Rerr);
        P2(tReg, "%-20s %13s %9s %8s %17s %17s %15s %13s\n", "region", "M4", "N", "CB", "data [nb]", "sim [nb]",
           "data/sim", "data/(R sim)");
        for (size_t r = 0; r < regs.size(); ++r) {
            const WinCmp& c = wc[r];
            const double ra = c.r(), re = c.re();
            P2(tReg, "%-20s %6.4f-%-6.4f %9.0f %8.0f %8.4g +- %-6.2g %8.4g +- %-6.2g %6.3f +- %-6.3f %6.3f +- %.3f\n",
               regs[r].name, regs[r].a, regs[r].b, c.N, c.CB, c.ds, c.de, c.ss, c.se, ra, re, ra / R, re / R);
        }
        P2(tReg, "# simulated templates per region: sigma [nb] (N_eff)\n");
        P2(tReg, "%-20s %22s %22s %22s %22s\n", "region", "T1 signal", "T2 phys. bkg", "T3 fake lepton", "T4 diff. sources");
        for (size_t r = 0; r < regs.size(); ++r) {
            if (r == 9 || r == 10) continue;
            const std::vector<WSum> t = SumMany(sd, SelAnd(sel, SelM4(regs[r].a, regs[r].b)),
                                                {CompSignal(), CompPhysBkg(), CompFake(), CompDiff()}, opt);
            std::string line = Form("%-20s", regs[r].name);
            for (const WSum& x : t) line += Form(" %8.4g (%9.1f)   ", x.Xsec(), x.Neff());
            P2(tReg, "%s\n", line.c_str());
        }

        // ---- histograms ----
        TH1D *hs5[cN], *hs10[cN], *hsz[cN], *hsz5[cN];
        for (int c = 0; c < cN; ++c) {
            hs5[c] = s5[c].Hist(Nm("s5"), 2); hs10[c] = s10[c].Hist(Nm("s10"), 2);
            hsz[c] = sz[c].Hist(Nm("sz"), 2); hsz5[c] = sz5[c].Hist(Nm("sz5"), 2);
        }
        TH1D* hN = d5.Hist(0, Nm("dN")); TH1D* hC = d5.Hist(3, Nm("dC")); TH1D* hS = d5.Hist(4, Nm("dS"));
        TH1D* hS10 = d10.Hist(4, Nm("dS10")); TH1D* hSz = dz.Hist(4, Nm("dSz")); TH1D* hSz5 = dz5.Hist(4, Nm("dSz5"));
        TH1D* hS25 = d25.Hist(4, Nm("dS25"));
        TH1D* hs25 = s25[cSim].Hist(Nm("s25"), 2);
        StyleM(hN, 1, 20, 1.0); StyleM(hC, 2, 21, 0.9); StyleM(hS, 4, 20, 1.0);
        StyleM(hS10, 4, 20, 1.2); StyleM(hSz, 4, 20, 1.0); StyleM(hSz5, 4, 20, 1.2);
        ratAbs[ip] = Ratio(hS10, hs10[cSim], Nm("rA"));
        ratShape[ip] = Ratio(hS10, hs10[cSim], Nm("rS"), R);

        // ---- chi2 of the shape ----
        std::vector<double> dRows10 = d10.n[0], sRows10 = s10[cSim].n, dRows25 = d25.n[0], sRows25 = s25[cSim].n;
        struct CR { const char* name; double a, b; int fine; };
        const std::vector<CR> creg = {{"below the bump", 0.30, 0.45, 0}, {"gap bump-peak", 0.525, 0.530, 1},
                                      {"above the peak", 0.56, 0.90, 0}, {"  0.56-0.70", 0.56, 0.70, 0},
                                      {"  0.70-0.90", 0.70, 0.90, 0},    {"  0.56-0.73", 0.56, 0.73, 0},
                                      {"  0.76-0.90", 0.76, 0.90, 0},    {"high 0.90-1.10", 0.90, 1.10, 0}};
        std::vector<Chi> chis;
        for (const CR& c : creg) {
            const Chi x = c.fine ? ShapeChi(hS25, hs25, dRows25, sRows25, c.a, c.b, R)
                                 : ShapeChi(hS10, hs10[cSim], dRows10, sRows10, c.a, c.b, R);
            chis.push_back(x);
            P2(tChi, "%-26s %-16s %5.3f-%-5.3f %5.1f/%-4d %10.2g %7.3f+-%-6.3f %5.1f/%-4d\n", PresetName(p), c.name, c.a,
               c.b, x.chiR, x.nR, TMath::Prob(x.chiR, x.nR), x.rF, x.rFe, x.chiF, x.nF);
        }

        // ---- figures: absolute (shape = 0) and after R (shape = 1) ----
        // shared y range of the log panels (from "pre", the largest selection)
        if (ip == 0) {
            const double ymax = MaxOf({hN, hs5[cSim]}, false);
            const double ymin = MinPos({hs5[cBD], hs5[cBC], hs5[cPC]});
            yLoLog = std::max(ymin * 0.3, ymax * 1e-5); yHiLog = ymax * 300;
        }
        for (int shape = 0; shape < 2; ++shape) {
            const double f = shape ? R : 1.0;
            TCanvas* c = new TCanvas(Nm("cM4"), "", 2600, 1900);
            c->cd();
            // header
            const std::string rtxt = shape ? std::string(Form("simulation #times R = %.3f #pm %.3f (R = data/sim in 0.60-0.80)", R, Rerr))
                                           : std::string("simulation NOT renormalised");
            Txt(0.01, 0.975, Form("%s M_{#pi^{+}#pi^{-}e^{+}e^{-}}: data against simulation, %s; selection %s",
                                  shape ? "SHAPE" : "ABSOLUTE", rtxt.c_str(), PresetName(p)), 0.017);
            Txt(0.01, 0.952, PresetLabel(p), 0.014);
            Txt(0.01, 0.931, DataLine().c_str(), 0.014, 4);
            Txt(0.01, 0.910, SimLine().c_str(), 0.014, kColSim);
            Txt(0.01, 0.889, kDirect, 0.014, 2);
            // left: log full range + ratio (10 MeV)
            PP L = MkPair(c, 0.0, 0.0, 0.55, 0.875, 0.30, true, 0.13, 0.03);
            L.up->cd();
            TH1D* fr = Frame(hN, yLoLog, yHiLog, kXM4, kYM4, 0.018, L.hu, false);
            fr->GetYaxis()->SetTitleOffset(0.95);
            fr->Draw("axis");
            TH1D* t2 = (TH1D*)hs5[cT2]->Clone(Nm("t2")); t2->Scale(f); StyleL(t2, kColPhysBkg, 1, 1, kColPhysFill); t2->Draw("hist same");
            TH1D* sm = (TH1D*)hs5[cSim]->Clone(Nm("sm")); sm->Scale(f); StyleL(sm, kColSim, 4, 1); sm->Draw("hist same");
            TH1D* k3[3];
            for (int k = 0; k < 3; ++k) { k3[k] = (TH1D*)hs5[cBD + k]->Clone(Nm("k3")); k3[k]->Scale(f); StyleL(k3[k], col[cBD + k], 3, 1); k3[k]->Draw("hist same"); }
            TH1D* t3 = (TH1D*)hs5[cT3]->Clone(Nm("t3")); t3->Scale(f); StyleL(t3, kColFake, 3, 2); t3->Draw("hist same");
            TH1D* t4 = (TH1D*)hs5[cT4]->Clone(Nm("t4")); t4->Scale(f); StyleL(t4, kColDiff, 3, 2); t4->Draw("hist same");
            hN->Draw("E1 P same"); hC->Draw("E1 P same"); hS->Draw("E1 P same");
            fr->Draw("axis same");
            TLegend* lg = new TLegend(0.14, 0.70, 0.97, 0.97);
            lg->SetNColumns(2); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.013 / L.hu);
            lg->AddEntry(hN, "data: N (all e^{+}e^{-} combinations)", "p");
            lg->AddEntry(sm, shape ? "sim: signal + phys. bkg, #times R" : "sim: signal + physical background", "l");
            lg->AddEntry(hC, "data: CB", "p");
            lg->AddEntry(t2, "sim: physical background", "f");
            lg->AddEntry(hS, "data: N - CB", "p");
            lg->AddEntry(k3[0], "sim signal: bump, #pi^{0} Dalitz", "l");
            lg->AddEntry(t3, "sim: fake lepton (not in the sum)", "l");
            lg->AddEntry(k3[1], "sim signal: bump, conversion", "l");
            lg->AddEntry(t4, "sim: different sources (not in the sum)", "l");
            lg->AddEntry(k3[2], "sim signal: peak, conversion", "l");
            lg->Draw();
            Txt(0.14, 0.64, "5 MeV bins", 0.013 / L.hu);
            L.lo->cd();
            TH1D* rr = shape ? ratShape[ip] : ratAbs[ip];
            StyleM(rr, 4, 20, 1.2);
            TH1D* fl = Frame(hS10, 0.3, 2.0, kXM4, shape ? "data / (R #times sim)" : "data / sim", 0.018, L.hl, true);
            fl->GetYaxis()->SetTitleOffset(0.95 * L.hl / L.hu); fl->GetXaxis()->SetTitleOffset(1.0);
            fl->Draw("axis");
            Band(scan::kBumpLo, scan::kBumpHi, 0.3, 2.0, kOrange);
            Band(scan::kPeakLo, scan::kPeakHi, 0.3, 2.0, kGreen);
            HLine(0.30, 1.10, 1.0, kGray + 2, 1);
            if (!shape) HLine(0.30, 1.10, R, kColSim, 2, 3);
            rr->Draw("E1 P same");
            Txt(0.14, 0.88, shape ? Form("(N - CB) / (R #times sim), 10 MeV bins; R = %.3f", R)
                                  : Form("(N - CB) / sim, 10 MeV bins; dashed: R = %.3f", R), 0.013 / L.hl);
            if (shape)
                Txt(0.14, 0.76, Form("#chi^{2}/ndf: < 0.45: %.0f/%d;  0.525-0.53: %.1f/%d;  0.56-0.90: %.0f/%d (0.56-0.70: %.0f/%d, 0.70-0.90: %.0f/%d)",
                                     chis[0].chiR, chis[0].nR, chis[1].chiR, chis[1].nR, chis[2].chiR, chis[2].nR,
                                     chis[3].chiR, chis[3].nR, chis[4].chiR, chis[4].nR), 0.012 / L.hl);
            // right: linear zoom 0.40-0.70 (2.5 MeV) + ratio (5 MeV)
            PP Rz = MkPair(c, 0.55, 0.0, 1.0, 0.875, 0.30, false, 0.15, 0.03);
            Rz.up->cd();
            TH1D* hb = (TH1D*)hsz[cT2]->Clone(Nm("hb")); hb->Scale(f); StyleL(hb, kColPhysBkg, 1, 1, kColPhysFill);
            const int fillc[3] = {kColBumpDalitz, kColBumpConv, kColPeakConv};
            TH1D* hc[3];
            for (int k = 0; k < 3; ++k) {
                hc[k] = (TH1D*)(k == 0 ? hb : hc[k - 1])->Clone(Nm("hc"));
                TH1D* add = (TH1D*)hsz[cBD + k]->Clone(Nm("ad")); add->Scale(f);
                hc[k]->Add(add);
                StyleL(hc[k], col[cBD + k], 1, 1, fillc[k]);
            }
            TH1D* hsum = (TH1D*)hsz[cSim]->Clone(Nm("hsum")); hsum->Scale(f); StyleL(hsum, kColSim, 3, 1);
            // y range: the same for the absolute and the shape figure of this preset
            double ym = std::max(MaxOf({hSz}, true), MaxOf({hsz[cSim]}, false) * std::max(1.0, R));
            TH1D* fz = Frame(hSz, 0, ym * 1.30, kXM4, kYM4, 0.018, Rz.hu, false, 1.0);
            fz->GetYaxis()->SetTitleOffset(1.25);
            fz->Draw("axis");
            for (int k = 2; k >= 0; --k) hc[k]->Draw("hist same");
            hb->Draw("hist same");
            hsum->Draw("hist same");
            hSz->Draw("E1 P same");
            fz->Draw("axis same");
            TLegend* lz = new TLegend(0.17, 0.66, 0.70, 0.97);
            lz->SetBorderSize(0); lz->SetFillStyle(0); lz->SetTextSize(0.013 / Rz.hu);
            lz->AddEntry(hSz, "data: N - CB", "p");
            lz->AddEntry(hsum, shape ? "sim: signal + phys. bkg, #times R" : "sim: signal + phys. bkg", "l");
            lz->AddEntry(hc[2], "sim: + peak, conversion (stacked)", "f");
            lz->AddEntry(hc[1], "sim: + bump, conversion (stacked)", "f");
            lz->AddEntry(hc[0], "sim: + bump, #pi^{0} Dalitz (stacked)", "f");
            lz->AddEntry(hb, "sim: physical background", "f");
            lz->Draw();
            Txt(0.17, 0.62, "linear zoom, 2.5 MeV bins", 0.013 / Rz.hu);
            Rz.lo->cd();
            TH1D* rz = Ratio(hSz5, hsz5[cSim], Nm("rz"), f); StyleM(rz, 4, 20, 1.2);
            TH1D* flz = Frame(hSz5, 0.3, 2.0, kXM4, shape ? "data / (R #times sim)" : "data / sim", 0.018, Rz.hl, true);
            flz->GetYaxis()->SetTitleOffset(1.25 * Rz.hl / Rz.hu);
            flz->Draw("axis");
            Band(scan::kBumpLo, scan::kBumpHi, 0.3, 2.0, kOrange);
            Band(scan::kPeakLo, scan::kPeakHi, 0.3, 2.0, kGreen);
            HLine(0.40, 0.70, 1.0, kGray + 2, 1);
            if (!shape) HLine(0.40, 0.70, R, kColSim, 2, 3);
            rz->Draw("E1 P same");
            Txt(0.17, 0.88, "5 MeV bins; bands: bump (orange) and peak (green) windows of the data", 0.013 / Rz.hl);
            Save(c, std::string(shape ? "datacomp_shape_" : "datacomp_m4_") + kTag[ip], true);
            delete c;
        }

        // ---- pion-pair classes of the physical background (pre only) ----
        if (p == kPre) {
            std::vector<Comp> pc;
            for (int k = 0; k <= 8; ++k) pc.push_back(CompBkgPcls(k));
            pc.push_back(CompSignal());
            const std::vector<WSpec> sp5 = FillMany(sd, sel, pc, opt, b5);
            const std::vector<WSpec> sp10 = FillMany(sd, sel, pc, opt, b10);
            const int pcol[10] = {kGray + 1, kOrange + 7, kMagenta + 1, kRed + 1, kBlue - 4, kCyan + 2, kSpring + 5, kViolet + 1, kGray + 3, kBlack};
            const int pst[10] = {2, 2, 1, 1, 1, 1, 1, 1, 9, 3};
            TCanvas* c = new TCanvas(Nm("cPc"), "", 2400, 1800);
            c->cd();
            Txt(0.01, 0.975, "Physical background of the simulation split by the origin of the pion pair; selection pre; simulation NOT renormalised", 0.017);
            Txt(0.01, 0.952, DataLine().c_str(), 0.014, 4);
            Txt(0.01, 0.931, SimLine().c_str(), 0.014, kColSim);
            Txt(0.01, 0.910, kDirect, 0.014, 2);
            PP L = MkPair(c, 0.0, 0.0, 1.0, 0.90, 0.38, true, 0.09, 0.03);
            L.up->cd();
            TH1D* fr = Frame(hN, yLoLog, MaxOf({hS, hs5[cSim]}, false) * 40, kXM4, kYM4, 0.018, L.hu, false);
            fr->GetYaxis()->SetTitleOffset(0.8);
            fr->Draw("axis");
            TH1D* sm = (TH1D*)hs5[cSim]->Clone(Nm("sm")); StyleL(sm, kColSim, 4, 1); sm->Draw("hist same");
            TLegend* lg = new TLegend(0.10, 0.76, 0.97, 0.97);
            lg->SetNColumns(4); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.013 / L.hu);
            lg->AddEntry(hS, "data: N - CB", "p");
            lg->AddEntry(sm, "sim: signal + phys. bkg", "l");
            TH1D* hp[10];
            for (int k = 0; k < 10; ++k) {
                hp[k] = sp5[k].Hist(Nm("hp"), 2); StyleL(hp[k], pcol[k], 3, pst[k]); hp[k]->Draw("hist same");
                lg->AddEntry(hp[k], Form("%s: %.3g nb", k < 9 ? PclsName(k) : "signal (T1)", sp5[k].S / kLsim), "l");
            }
            hS->Draw("E1 P same");
            fr->Draw("axis same");
            lg->Draw();
            L.lo->cd();
            TH1D* fl = Frame(hS10, 0.0, 1.0, kXM4, "fraction of sim", 0.018, L.hl, true);
            fl->GetYaxis()->SetTitleOffset(0.8 * L.hl / L.hu);
            fl->Draw("axis");
            for (int k = 0; k < 10; ++k) {
                TH1D* fr10 = sp10[k].Hist(Nm("f10"), 2);
                fr10->Divide(hs10[cSim]);
                for (int i = 1; i <= fr10->GetNbinsX(); ++i) fr10->SetBinError(i, 0);
                StyleL(fr10, pcol[k], 3, pst[k]); fr10->Draw("hist same");
            }
            Txt(0.10, 0.88, "fraction of each pion-pair class in the simulated spectrum (signal + physical background), 10 MeV bins", 0.013 / L.hl);
            Save(c, "datacomp_m4_pcls_pre", true);
            delete c;
        }
        printf("Part1: preset %s done (%.0f s)\n", PresetName(p), sw.RealTime()); sw.Continue();
    }
    if (tReg) fclose(tReg);
    if (tChi) fclose(tChi);
    if (tR) fclose(tR);

    // ---- overlay of the ratios of the four presets ----
    {
        TCanvas* c = new TCanvas(Nm("cOv"), "", 2400, 1600);
        c->cd();
        Txt(0.01, 0.97, "Ratio data (N - CB) / simulation (signal + physical background) for four selections, 10 MeV bins", 0.020);
        Txt(0.01, 0.94, DataLine().c_str(), 0.016, 4);
        Txt(0.01, 0.91, (SimLine() + "; " + kDirect).c_str(), 0.014, kColSim);
        const int pcol[kNP] = {kBlue, kRed + 1, kGreen + 2, kMagenta + 1};
        const int pst[kNP] = {20, 21, 22, 33};
        for (int sh = 0; sh < 2; ++sh) {
            c->cd();
            TPad* pd = new TPad(Nm("ov"), "", 0, sh ? 0.0 : 0.45, 1, sh ? 0.45 : 0.90);
            pd->SetLeftMargin(0.08); pd->SetRightMargin(0.02); pd->SetTopMargin(0.04); pd->SetBottomMargin(0.16);
            pd->Draw(); pd->cd();
            TH1D* fr = Frame(ratAbs[0], 0.3, 2.0, kXM4, sh ? "data / (R #times sim)" : "data / sim", 0.020, 0.45, true);
            fr->GetYaxis()->SetTitleOffset(0.7);
            fr->Draw("axis");
            Band(scan::kBumpLo, scan::kBumpHi, 0.3, 2.0, kOrange);
            Band(scan::kPeakLo, scan::kPeakHi, 0.3, 2.0, kGreen);
            HLine(0.30, 1.10, 1.0, kGray + 2, 1);
            TLegend* lg = new TLegend(0.36, 0.76, 0.97, 0.95);
            lg->SetNColumns(2); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.015 / 0.45);
            for (int ip = 0; ip < kNP; ++ip) {
                TH1D* h = (TH1D*)(sh ? ratShape[ip] : ratAbs[ip])->Clone(Nm("ovh"));
                // shift the markers a little for readability
                StyleM(h, pcol[ip], pst[ip], 1.3);
                h->Draw("E1 P same");
                std::string lab = PresetName(kPresets[ip]);
                if (sh) lab += Form(" (R = %.3f #pm %.3f)", Rv[ip], Re[ip]);
                lg->AddEntry(h, lab.c_str(), "p");
            }
            lg->Draw();
            Txt(0.09, 0.19, sh ? "(b) after R (R = data/sim in 0.60-0.80 of each selection)" : "(a) absolute (L_{data} = 6080.5 nb^{-1}, unconfirmed)", 0.016 / 0.45);
        }
        Save(c, "datacomp_shape_ratio_presets", true);
        delete c;
    }
}

}  // namespace dcmp

// ============================================================================
// Part 2: OA and M(ee) classes; content of the data CB
// ============================================================================
namespace dcmp {

inline int OAcls(double oa) { return oa < 2 ? 0 : oa < 4 ? 1 : oa < 6 ? 2 : oa < 10 ? 3 : 4; }
inline int MEEcls(double m) { return m < 0.005 ? 0 : m < 0.010 ? 1 : m < 0.020 ? 2 : m < 0.050 ? 3 : 4; }
const char* kOAname[5] = {"OA < 2", "2-4", "4-6", "6-10", "> 10"};
const char* kMEEname[5] = {"M < 5", "5-10", "10-20", "20-50", "> 50"};
inline bool T3Same(const Row& r) { return r.tmpl == 3 && (r.flags & simtruth::kFlagSameParticle); }
inline bool T3Ghost(const Row& r) { return r.tmpl == 3 && !(r.flags & simtruth::kFlagSameParticle) && (r.flags & simtruth::kFlagGhostLepton); }
inline bool T3Misid(const Row& r) { return r.tmpl == 3 && !(r.flags & simtruth::kFlagSameParticle) && !(r.flags & simtruth::kFlagGhostLepton); }

void Part2(const SimData& sd, const scan::Data& dd) {
    const int p = kPre;
    const Sel sel = SelPreset(p);
    const scan::Sel dsel = DataSelPreset(p);
    double Rerr = 0;
    const double R = GetR(sd, dd, p, &Rerr);
    printf("Part2: R(pre) = %.4f +- %.4f\n", R, Rerr);
    struct Reg { const char* name; double a, b; };
    const std::vector<Reg> regs = {{"M4 0.56-0.80", 0.56, 0.80}, {"M4 0.30-0.45", 0.30, 0.45}, {"M4 0.30-1.10", 0.30, 1.10}};
    std::vector<std::pair<double, double>> wl;
    for (const Reg& r : regs) wl.push_back(std::make_pair(r.a, r.b));

    // ---- data: class selections, window sums ----
    // 0 all; 1-5 OA classes; 6-10 M(ee) classes; 11-15 OA classes with shared ring; 16 all with shared ring
    std::vector<scan::Sel> ds;
    ds.push_back(dsel);
    for (int k = 0; k < 5; ++k) ds.push_back([dsel, k](const scan::Row& r) { return dsel(r) && OAcls(r.oa) == k; });
    for (int k = 0; k < 5; ++k) ds.push_back([dsel, k](const scan::Row& r) { return dsel(r) && MEEcls(r.mee) == k; });
    for (int k = 0; k < 5; ++k)
        ds.push_back([dsel, k](const scan::Row& r) { return dsel(r) && OAcls(r.oa) == k && r.shared_ring == 1; });
    ds.push_back([dsel](const scan::Row& r) { return dsel(r) && r.shared_ring == 1; });
    scan::FillOpt fo; fo.nbin = 320; fo.lo = 0.30; fo.hi = 1.10; fo.reb = 4; fo.windows = wl;
    const std::vector<scan::Spectra> dsp = scan::FillMany(dd, ds, fo);
    auto DW = [&](int s, int r) -> const scan::WinSum& { return dsp[s].win[2 + r]; };

    // ---- simulation: per class (binning = class index), per region ----
    std::vector<Comp> comps = {CompSim(), CompSignal(), CompPhysBkg(), CompFake(), CompDiff(),
                               {"T3 same particle", T3Same}, {"T3 ghost lepton", T3Ghost}, {"T3 misidentified", T3Misid},
                               {"T1+T2 shared ring", [](const Row& r) { return (r.tmpl == 1 || r.tmpl == 2) && r.shared_ring == 1; }},
                               {"T3 shared ring", [](const Row& r) { return r.tmpl == 3 && r.shared_ring == 1; }},
                               {"T4 shared ring", [](const Row& r) { return r.tmpl == 4 && r.shared_ring == 1; }}};
    enum { kSim = 0, kT1, kT2, kT3, kT4, kT3s, kT3g, kT3m, kSimSh, kT3Sh, kT4Sh, kNC };
    Binning bOA; bOA.nbin = 5; bOA.lo = 0; bOA.hi = 5; bOA.x = [](const Row& r) { return OAcls(r.oa) + 0.5; };
    Binning bME; bME.nbin = 5; bME.lo = 0; bME.hi = 5; bME.x = [](const Row& r) { return MEEcls(r.mee) + 0.5; };
    std::vector<std::vector<WSpec>> sOA(regs.size()), sME(regs.size()), sOAB(regs.size()), sMEB(regs.size());
    SimOpt oB; oB.mode = 1;
    for (size_t r = 0; r < regs.size(); ++r) {
        const Sel s = SelAnd(sel, SelM4(regs[r].a, regs[r].b));
        sOA[r] = FillMany(sd, s, comps, SimOpt(), bOA);
        sME[r] = FillMany(sd, s, comps, SimOpt(), bME);
        sOAB[r] = FillMany(sd, s, std::vector<Comp>(1, CompSim()), oB, bOA);
        sMEB[r] = FillMany(sd, s, std::vector<Comp>(1, CompSim()), oB, bME);
    }
    auto Xs = [](const WSpec& w, int k) { return w.s[k] / kLsim; };
    auto Xe = [](const WSpec& w, int k) { return std::sqrt(w.v[k]) / kLsim; };
    auto RatioE = [](double a, double ae, double b, double be, double& r, double& re) {
        r = b > 0 ? a / b : 0;
        re = b > 0 ? std::sqrt(ae * ae + r * r * be * be) / b : 0;
    };

    // ---- table: classes ----
    FILE* t = OpenTab("datacomp_classes.txt");
    P2(t, "# OA(e+e-) and M(e+e-) classes: data (N - CB) / L_data against the simulation (signal + physical background)\n");
    P2(t, "# selection pre; data L_data = %.1f nb^-1 (UNCONFIRMED), CB in 10 MeV M4 bins per class, group variances\n", kLdata);
    P2(t, "# sim std + lep owned rows, w (PDG), L_sim = %.2f nb^-1; mode A default, mode B as a variant (HGeant pi0 Dalitz kinematics)\n", kLsim);
    P2(t, "# after R: R = %.4f +- %.4f (pre, 0.60-0.80); the error of R is not included in the column 'data/(R sim)'\n", R, Rerr);
    double rat[2][3][5], ratE[2][3][5];   // [OA/MEE][region][class]: absolute ratio, mode A
    double ratB[2][3][5], ratBE[2][3][5];
    for (int kind = 0; kind < 2; ++kind) {
        for (size_t r = 0; r < regs.size(); ++r) {
            P2(t, "\n## %s classes, %s\n", kind ? "M(e+e-) [MeV/c2]" : "OA(e+e-) [deg]", regs[r].name);
            P2(t, "%-8s %9s %8s %6s %17s %22s %15s %13s %15s %14s\n", "class", "N", "CB", "CB/N", "data [nb]",
               "sim A [nb] (N_eff)", "data/sim A", "data/(R sim)", "data/sim B", "sim T1 [nb]");
            const std::vector<WSpec>& S = kind ? sME[r] : sOA[r];
            const std::vector<WSpec>& SB = kind ? sMEB[r] : sOAB[r];
            double tot[6] = {0, 0, 0, 0, 0, 0};
            for (int k = 0; k < 5; ++k) {
                const scan::WinSum& w = DW(1 + 5 * kind + k, r);
                const double ds = w.sig() / kLdata, de = w.sigErr() / kLdata;
                double ra, re, rb, rbe;
                RatioE(ds, de, Xs(S[kSim], k), Xe(S[kSim], k), ra, re);
                RatioE(ds, de, Xs(SB[0], k), Xe(SB[0], k), rb, rbe);
                rat[kind][r][k] = ra; ratE[kind][r][k] = re; ratB[kind][r][k] = rb; ratBE[kind][r][k] = rbe;
                P2(t, "%-8s %9.0f %8.0f %6.3f %8.4g +- %-6.2g %8.4g +- %-5.2g (%7.0f) %6.3f +- %-5.3f %6.3f +- %-5.3f %6.3f +- %-5.3f %8.4g\n",
                   kind ? kMEEname[k] : kOAname[k], w.n[0], w.cb, w.n[0] > 0 ? w.cb / w.n[0] : 0, ds, de, Xs(S[kSim], k),
                   Xe(S[kSim], k), S[kSim].Neff(k), ra, re, ra / R, re / R, rb, rbe, Xs(S[kT1], k));
                tot[0] += w.n[0]; tot[1] += w.cb;
            }
            const scan::WinSum& wa = DW(0, r);
            P2(t, "%-8s %9.0f %8.0f %6.3f %8.4g +- %-6.2g %8.4g +- %-5.2g (%7.0f) %6.3f\n", "all", wa.n[0], wa.cb,
               wa.cb / wa.n[0], wa.sig() / kLdata, wa.sigErr() / kLdata, S[kSim].S / kLsim, std::sqrt(S[kSim].V) / kLsim,
               S[kSim].NeffIn(), (wa.sig() / kLdata) / (S[kSim].S / kLsim));
            P2(t, "# sum of the class CBs %.0f vs CB of all %.0f (CB is not additive: geometric mean per class)\n", tot[1], wa.cb);
        }
    }
    // shared ring
    P2(t, "\n## share of pairs whose two leptons share one RICH ring, per OA class\n");
    P2(t, "# data N: N_shared / N; data N - CB: (N - CB)_shared / (N - CB); sim: sum w shared / sum w (mode A); errors approximate (subset formula)\n");
    double shD[3][5], shDe[3][5], shS[3][5], shSe[3][5], shN[3][5], shT3[3][5];
    for (size_t r = 0; r < regs.size(); ++r) {
        P2(t, "### %s\n%-8s %16s %16s %16s %16s %16s\n", regs[r].name, "class", "data N", "data N - CB", "sim T1+T2",
           "sim T3", "sim T4");
        for (int k = 0; k < 5; ++k) {
            const scan::WinSum& wa = DW(1 + k, r);
            const scan::WinSum& wsh = DW(11 + k, r);
            const double pN = wa.n[0] > 0 ? wsh.n[0] / wa.n[0] : 0, pNe = wa.n[0] > 0 ? std::sqrt(pN * (1 - pN) / wa.n[0]) : 0;
            const double num = wsh.sig(), den = wa.sig(), vn = std::pow(wsh.sigErr(), 2), vd = std::pow(wa.sigErr(), 2);
            const double pS = den > 0 ? num / den : 0;
            const double pSe = den > 0 ? std::sqrt(std::max(0., vn * (1 - 2 * pS) + pS * pS * vd)) / den : 0;
            const WSpec& S = sOA[r][kSim]; const WSpec& Sh = sOA[r][kSimSh];
            const double q = S.s[k] > 0 ? Sh.s[k] / S.s[k] : 0;
            const double qe = S.s[k] > 0 ? std::sqrt(std::max(0., Sh.v[k] * (1 - 2 * q) + q * q * S.v[k])) / S.s[k] : 0;
            const double q3 = sOA[r][kT3].s[k] > 0 ? sOA[r][kT3Sh].s[k] / sOA[r][kT3].s[k] : 0;
            const double q4 = sOA[r][kT4].s[k] > 0 ? sOA[r][kT4Sh].s[k] / sOA[r][kT4].s[k] : 0;
            shN[r][k] = pN; shD[r][k] = pS; shDe[r][k] = pSe; shS[r][k] = q; shSe[r][k] = qe; shT3[r][k] = q3;
            P2(t, "%-8s %7.4f +- %-6.4f %7.4f +- %-6.4f %7.4f +- %-6.4f %16.4f %16.4f\n", kOAname[k], pN, pNe, pS, pSe, q, qe, q3, q4);
        }
    }
    if (t) fclose(t);

    // ---- table: content of the data CB against the simulated information classes ----
    FILE* tc = OpenTab("datacomp_cb_content.txt");
    P2(tc, "# SIDE NOTE, NOT A RESULT: the simulation has no combinatorial background, so background levels of data and simulation are not compared\n");
    P2(tc, "# (user rule). No conclusion is drawn from the ratios below. About 90 %% of T3 are ghost-flagged true pairs (tables/review_ghost.txt).\n");
    P2(tc, "# What the CB subtraction of the data removes, against the simulated information classes (not part of the simulated spectrum)\n");
    P2(tc, "# T3 fake lepton = misidentified lepton, ghost track in a lepton slot, or the same particle in both lepton slots;\n");
    P2(tc, "# T4 different sources = true leptons from two different decays/conversions (L7-L9). Selection pre, mode A, w (PDG)\n");
    P2(tc, "# data CB = 2 sqrt(N++ N--) in 10 MeV M4 bins, no k; sigma [nb] with L_data = %.1f nb^-1 (UNCONFIRMED), L_sim = %.2f nb^-1\n", kLdata, kLsim);
    P2(tc, "# columns: data CB, data N - CB [nb]; sim T1+T2, T3 (split: same particle / ghost / misidentified), T4 [nb] (N_eff); ratios to the data CB\n");
    double cbT3[3][5], cbT4[3][5], cbT34e[3][5], cbT34[3][5];
    for (size_t r = 0; r < regs.size(); ++r) {
        P2(tc, "\n## %s\n%-8s %16s %10s %10s %18s %8s %8s %8s %16s %7s %7s %13s\n", regs[r].name, "OA class", "data CB [nb]",
           "data N-CB", "sim T1+T2", "sim T3 (N_eff)", "same", "ghost", "misid", "sim T4 (N_eff)", "T3/CB", "T4/CB", "(T3+T4)/CB");
        for (int k = 0; k <= 5; ++k) {   // 5 = all classes
            const scan::WinSum& w = k < 5 ? DW(1 + k, r) : DW(0, r);
            const double cb = w.cb / kLdata, cbe = std::sqrt(w.cbVar) / kLdata, sg = w.sig() / kLdata;
            auto G = [&](int c, double& s, double& v, double& ne) {
                const WSpec& x = sOA[r][c];
                if (k < 5) { s = x.s[k]; v = x.v[k]; } else { s = x.S; v = x.V; }
                ne = v > 0 ? s * s / v : 0; s /= kLsim; v /= kLsim * kLsim;
            };
            double s0, v0, n0, s3, v3, n3, s4, v4, n4, ss, vs, ns, sg3, vg, ng, sm, vm, nm;
            G(kSim, s0, v0, n0); G(kT3, s3, v3, n3); G(kT4, s4, v4, n4); G(kT3s, ss, vs, ns); G(kT3g, sg3, vg, ng); G(kT3m, sm, vm, nm);
            double r3, r3e, r4, r4e, r34, r34e;
            RatioE(s3, std::sqrt(v3), cb, cbe, r3, r3e);   // ratio sim / data: same propagation
            RatioE(s4, std::sqrt(v4), cb, cbe, r4, r4e);
            RatioE(s3 + s4, std::sqrt(v3 + v4), cb, cbe, r34, r34e);
            if (k < 5) { cbT3[r][k] = r3; cbT4[r][k] = r4; cbT34[r][k] = r34; cbT34e[r][k] = r34e; }
            P2(tc, "%-8s %8.4g +- %-6.2g %10.4g %10.4g %9.4g (%6.0f) %8.4g %8.4g %8.4g %8.4g (%6.0f) %7.2f %7.3f %6.2f +- %.2f\n",
               k < 5 ? kOAname[k] : "all", cb, cbe, sg, s0, s3, n3, ss, sg3, sm, s4, n4, r3, r4, r34, r34e);
        }
    }
    if (tc) fclose(tc);

    // ======================= figure: OA and M(ee) distributions =======================
    {
        struct DV { const char* name; const char* xt; const char* yt; int nb; double lo, hi; };
        const DV dv[2] = {{"oa", "OA(e^{+}e^{-}) [deg]", "d#sigma/dOA [nb/deg]", 60, 0.0, 30.0},
                          {"mee", "M_{e^{+}e^{-}} [GeV/c^{2}]", "d#sigma/dM_{ee} [nb/(GeV/c^{2})]", 60, 0.0, 0.15}};
        TCanvas* c = new TCanvas(Nm("cOM"), "", 2600, 2200);
        c->cd();
        Txt(0.01, 0.978, "OA(e^{+}e^{-}) and M(e^{+}e^{-}) of data and simulation in two M4 regions; selection pre (isBest #geq 0, PT3, v_{z} > -160 mm, no OA cut)", 0.016);
        Txt(0.01, 0.958, DataLine().c_str(), 0.013, 4);
        Txt(0.01, 0.940, SimLine().c_str(), 0.013, kColSim);
        Txt(0.01, 0.922, Form("solid: simulation absolute; dashed: simulation #times R, R = %.3f #pm %.3f (pre, 0.60-0.80); dotted lines: class boundaries; %s", R, Rerr, "#eta #rightarrow #pi^{+}#pi^{-}e^{+}e^{-} NOT in sim"), 0.013, 1);
        const double rlo[2] = {0.56, 0.30}, rhi[2] = {0.80, 0.45};
        for (int v = 0; v < 2; ++v)
            for (int r = 0; r < 2; ++r) {
                const double x1 = r * 0.5, x2 = x1 + 0.5, y2 = v == 0 ? 0.91 : 0.455, y1 = v == 0 ? 0.455 : 0.0;
                const DV& D = dv[v];
                const Sel s = SelAnd(sel, SelM4(rlo[r], rhi[r]));
                Binning b; b.nbin = D.nb; b.lo = D.lo; b.hi = D.hi;
                if (v == 0) b.x = [](const Row& x) { return (double)x.oa; };
                else b.x = [](const Row& x) { return (double)x.mee; };
                const std::vector<WSpec> S = FillMany(sd, s, {CompSim(), CompFake(), CompDiff()}, SimOpt(), b);
                const double a = rlo[r], bb = rhi[r];
                const scan::Sel dsr = [dsel, a, bb](const scan::Row& x) { return dsel(x) && x.m4 >= a && x.m4 < bb; };
                const DH H = FillData(dd, dsr, D.nb, D.lo, D.hi, 1,
                                      v == 0 ? DVar([](const scan::Row& x) { return (double)x.oa; })
                                             : DVar([](const scan::Row& x) { return (double)x.mee; }));
                TH1D* hN = H.Hist(0, Nm("oN")); TH1D* hC = H.Hist(3, Nm("oC")); TH1D* hS = H.Hist(4, Nm("oS"));
                StyleM(hN, 1, 20, 0.9); StyleM(hC, 2, 21, 0.8); StyleM(hS, 4, 20, 0.9);
                TH1D* sS = S[0].Hist(Nm("oSim"), 2); StyleL(sS, kColSim, 4, 1);
                TH1D* sR = (TH1D*)sS->Clone(Nm("oSR")); sR->Scale(R); StyleL(sR, kColSim, 3, 2);
                TH1D* s3 = S[1].Hist(Nm("o3"), 2); StyleL(s3, kColFake, 3, 2);
                TH1D* s4 = S[2].Hist(Nm("o4"), 2); StyleL(s4, kColDiff, 3, 2);
                PP P = MkPair(c, x1, y1, x2, y2, 0.32, true, 0.15, 0.10);
                P.up->cd();
                const double ymax = MaxOf({hN, sS, sR}, false), ymin = std::max(MinPos({hS, sS}) * 0.3, ymax * 1e-4);
                TH1D* fr = Frame(hN, ymin, ymax * 30, D.xt, D.yt, 0.016, P.hu, false);
                fr->GetYaxis()->SetTitleOffset(1.2);
                fr->Draw("axis");
                sS->Draw("hist same"); sR->Draw("hist same"); s3->Draw("hist same"); s4->Draw("hist same");
                hN->Draw("E1 P same"); hC->Draw("E1 P same"); hS->Draw("E1 P same");
                const double edges[4] = {v == 0 ? 2. : 0.005, v == 0 ? 4. : 0.010, v == 0 ? 6. : 0.020, v == 0 ? 10. : 0.050};
                for (double e : edges) VLine(e, ymin, ymax * 30, kGray + 1, 3);
                fr->Draw("axis same");
                TLegend* lg = new TLegend(0.40, 0.55, 0.97, 0.89);
                lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.012 / P.hu);
                lg->AddEntry(hN, "data N", "p"); lg->AddEntry(hC, "data CB", "p"); lg->AddEntry(hS, "data N - CB", "p");
                lg->AddEntry(sS, "sim signal + phys. bkg", "l"); lg->AddEntry(sR, "sim #times R", "l");
                lg->AddEntry(s3, "sim fake lepton (info)", "l"); lg->AddEntry(s4, "sim different sources (info)", "l");
                lg->Draw();
                Txt(0.16, 0.92, Form("(%c) %s, M_{#pi^{+}#pi^{-}e^{+}e^{-}} %.2f-%.2f GeV/c^{2}", 'a' + 2 * v + r,
                                     v == 0 ? "OA(e^{+}e^{-})" : "M(e^{+}e^{-})", rlo[r], rhi[r]), 0.014 / P.hu);
                P.lo->cd();
                TH1D* rr = Ratio(hS, sS, Nm("oR")); StyleM(rr, 4, 20, 1.0);
                TH1D* fl = Frame(hS, 0.0, 2.5, D.xt, "data / sim", 0.016, P.hl, true);
                fl->GetYaxis()->SetTitleOffset(1.2 * P.hl / P.hu);
                fl->Draw("axis");
                HLine(D.lo, D.hi, 1.0, kGray + 2, 1);
                HLine(D.lo, D.hi, R, kColSim, 2, 3);
                for (double e : edges) VLine(e, 0.0, 2.5, kGray + 1, 3);
                rr->Draw("E1 P same");
            }
        Save(c, "datacomp_oa_mee_dist", true);
        delete c;
    }

    // ======================= figure: class ratios and shared ring =======================
    {
        TCanvas* c = new TCanvas(Nm("cCl"), "", 2400, 1900);
        c->cd();
        Txt(0.01, 0.975, "Data (N - CB) / simulation (signal + physical background) per OA(e^{+}e^{-}) and M(e^{+}e^{-}) class; share of pairs with a shared RICH ring; selection pre", 0.016);
        Txt(0.01, 0.953, DataLine().c_str(), 0.013, 4);
        Txt(0.01, 0.934, (SimLine() + Form("; R = %.3f #pm %.3f (0.60-0.80)", R, Rerr)).c_str(), 0.013, kColSim);
        Txt(0.01, 0.915, kDirect, 0.013, 2);
        const int rcol[2] = {kBlue, kRed + 1};
        for (int pnl = 0; pnl < 4; ++pnl) {
            c->cd();
            const double x1 = (pnl % 2) * 0.5, y2 = pnl < 2 ? 0.90 : 0.45;
            TPad* pd = new TPad(Nm("cl"), "", x1, y2 - 0.45, x1 + 0.5, y2);
            pd->SetLeftMargin(0.14); pd->SetRightMargin(0.03); pd->SetTopMargin(0.10); pd->SetBottomMargin(0.14);
            pd->Draw(); pd->cd();
            TH1D* fr = new TH1D(Nm("clf"), "", 5, 0, 5);
            fr->SetDirectory(nullptr);
            for (int k = 0; k < 5; ++k) fr->GetXaxis()->SetBinLabel(k + 1, pnl == 1 ? kMEEname[k] : kOAname[k]);
            fr->GetXaxis()->SetLabelSize(0.065); fr->GetYaxis()->SetLabelSize(0.050); fr->GetYaxis()->SetTitleSize(0.055);
            fr->GetXaxis()->SetTitleSize(0.055); fr->GetYaxis()->SetTitleOffset(1.15);
            fr->GetXaxis()->SetTitle(pnl == 1 ? "M(e^{+}e^{-}) class [MeV/c^{2}]" : "OA(e^{+}e^{-}) class [deg]");
            if (pnl < 2) {
                fr->SetMinimum(0.0); fr->SetMaximum(2.4); fr->GetYaxis()->SetTitle("data (N - CB) / sim");
                fr->Draw("axis");
                HLine(0, 5, 1.0, kGray + 2, 1); HLine(0, 5, R, kColSim, 2, 3);
                TLegend* lg = new TLegend(0.16, 0.66, 0.97, 0.89);
                lg->SetNColumns(2); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.040);
                for (int r = 0; r < 2; ++r) {
                    TH1D* hA = new TH1D(Nm("ha"), "", 5, 0, 5); hA->SetDirectory(nullptr);
                    TH1D* hB = new TH1D(Nm("hb"), "", 5, 0, 5); hB->SetDirectory(nullptr);
                    for (int k = 0; k < 5; ++k) {
                        hA->SetBinContent(k + 1, rat[pnl][r][k]); hA->SetBinError(k + 1, ratE[pnl][r][k]);
                        hB->SetBinContent(k + 1, ratB[pnl][r][k]); hB->SetBinError(k + 1, ratBE[pnl][r][k]);
                    }
                    StyleM(hA, rcol[r], 20, 2.0); StyleM(hB, rcol[r], 24, 2.0);
                    hA->Draw("E1 P X0 same"); hB->Draw("E1 P X0 same");
                    lg->AddEntry(hA, Form("%s, mode A", regs[r].name), "p");
                    lg->AddEntry(hB, Form("%s, mode B", regs[r].name), "p");
                }
                lg->Draw();
                Txt(0.14, 0.92, Form("(%c) %s; dashed: R = %.3f", 'a' + pnl,
                                     pnl ? "M(e^{+}e^{-}) classes" : "OA classes", R), 0.045);
            } else {
                const int r = pnl - 2;
                fr->SetMinimum(0.0); fr->SetMaximum(1.15); fr->GetYaxis()->SetTitle("share with a shared RICH ring");
                fr->Draw("axis");
                TH1D* h[4];
                const char* ln[4] = {"data N", "data N - CB", "sim signal + phys. bkg", "sim fake lepton (info)"};
                const int cc[4] = {1, 4, kColSim, kColFake}, st[4] = {24, 20, 21, 25};
                TLegend* lg = new TLegend(0.52, 0.55, 0.97, 0.88);
                lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.042);
                for (int q = 0; q < 4; ++q) {
                    h[q] = new TH1D(Nm("hs"), "", 5, 0, 5); h[q]->SetDirectory(nullptr);
                    for (int k = 0; k < 5; ++k) {
                        const double v = q == 0 ? shN[r][k] : q == 1 ? shD[r][k] : q == 2 ? shS[r][k] : shT3[r][k];
                        const double e = q == 1 ? shDe[r][k] : q == 2 ? shSe[r][k] : 0;
                        h[q]->SetBinContent(k + 1, v); h[q]->SetBinError(k + 1, e);
                    }
                    StyleM(h[q], cc[q], st[q], 2.0);
                    h[q]->Draw("E1 P X0 same");
                    lg->AddEntry(h[q], ln[q], "p");
                }
                lg->Draw();
                Txt(0.14, 0.92, Form("(%c) shared-RICH-ring share per OA class, %s", 'a' + pnl, regs[r].name), 0.045);
            }
        }
        Save(c, "datacomp_classes", true);
        delete c;
    }

    // ======================= figure: content of the data CB =======================
    {
        TCanvas* c = new TCanvas(Nm("cCB"), "", 2400, 1900);
        c->cd();
        Txt(0.01, 0.975, "SIDE NOTE, NOT A RESULT (the simulation has no CB; background levels of data and simulation are not compared): data CB next to the simulated classes T3 and T4; pre", 0.016, 2);
        Txt(0.01, 0.953, Form("data: CB = 2#sqrt{N_{++}N_{--}} / L_{data} / bin width (no k), %s; simulation: SMASH std + lep, mode A, w (PDG), #Sigma w / L_{sim}, L_{sim} = 221.05 nb^{-1}", kLdataNote), 0.013, 4);
        Txt(0.01, 0.934, "T3 and T4 are NOT part of the simulated spectrum used for the comparisons; about 90 % of T3 are ghost-flagged TRUE pairs (review finding 1, tables/review_ghost.txt)", 0.013, kColSim);
        // (a) M4
        {
            c->cd();
            TPad* pd = new TPad(Nm("cb"), "", 0, 0.465, 0.5, 0.925);
            pd->SetLeftMargin(0.14); pd->SetRightMargin(0.03); pd->SetTopMargin(0.10); pd->SetBottomMargin(0.14); pd->SetLogy();
            pd->Draw(); pd->cd();
            Binning b; b.nbin = 80; b.lo = 0.30; b.hi = 1.10;
            const std::vector<WSpec> S = FillMany(sd, sel, {CompSim(), CompFake(), CompDiff(), {"T3+T4", [](const Row& r) { return r.tmpl == 3 || r.tmpl == 4; }}}, SimOpt(), b);
            const DH H = FillData(dd, dsel, 80, 0.30, 1.10, 1);
            TH1D* hC = H.Hist(3, Nm("cC")); TH1D* hS = H.Hist(4, Nm("cS"));
            StyleM(hC, 2, 21, 1.1); StyleM(hS, 4, 20, 1.1);
            TH1D* s0 = S[0].Hist(Nm("c0"), 2); StyleL(s0, kColSim, 3, 1);
            TH1D* s3 = S[1].Hist(Nm("c3"), 2); StyleL(s3, kColFake, 3, 2);
            TH1D* s4 = S[2].Hist(Nm("c4"), 2); StyleL(s4, kColDiff, 3, 2);
            TH1D* s34 = S[3].Hist(Nm("c34"), 2); StyleL(s34, 1, 3, 1);
            const double ymax = MaxOf({s0, hS}, false);
            TH1D* fr = Frame(hC, ymax * 2e-5, ymax * 2e3, kXM4, kYM4, 0.016, 0.46, true);
            fr->GetYaxis()->SetTitleOffset(1.2);
            fr->Draw("axis");
            s0->Draw("hist same"); s3->Draw("hist same"); s4->Draw("hist same"); s34->Draw("hist same");
            hS->Draw("E1 P same"); hC->Draw("E1 P same");
            TLegend* lg = new TLegend(0.16, 0.66, 0.97, 0.89);
            lg->SetNColumns(2); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.038);
            lg->AddEntry(hC, "data CB", "p"); lg->AddEntry(s34, "sim T3 + T4", "l");
            lg->AddEntry(hS, "data N - CB", "p"); lg->AddEntry(s3, "sim T3 fake lepton", "l");
            lg->AddEntry(s0, "sim signal + phys. bkg", "l"); lg->AddEntry(s4, "sim T4 different sources", "l");
            lg->Draw();
            Txt(0.14, 0.92, "(a) M_{#pi^{+}#pi^{-}e^{+}e^{-}}, 10 MeV bins, simulation not renormalised", 0.045);
        }
        // (b) OA and (d) M(ee) in M4 0.30-1.10
        for (int v = 0; v < 2; ++v) {
            c->cd();
            TPad* pd = new TPad(Nm("cb"), "", 0.5, v ? 0.0 : 0.465, 1.0, v ? 0.465 : 0.925);
            pd->SetLeftMargin(0.14); pd->SetRightMargin(0.03); pd->SetTopMargin(0.10); pd->SetBottomMargin(0.14); pd->SetLogy();
            pd->Draw(); pd->cd();
            Binning b; b.nbin = 60; b.lo = 0; b.hi = v ? 0.15 : 30.0;
            if (v == 0) b.x = [](const Row& x) { return (double)x.oa; };
            else b.x = [](const Row& x) { return (double)x.mee; };
            const Sel s = SelAnd(sel, SelM4(0.30, 1.10));
            const std::vector<WSpec> S = FillMany(sd, s, {CompFake(), CompDiff(), {"same", T3Same}, {"ghost", T3Ghost}, {"misid", T3Misid}}, SimOpt(), b);
            const scan::Sel dsr = [dsel](const scan::Row& x) { return dsel(x) && x.m4 >= 0.30 && x.m4 < 1.10; };
            const DH H = FillData(dd, dsr, 60, 0, v ? 0.15 : 30.0, 1,
                                  v == 0 ? DVar([](const scan::Row& x) { return (double)x.oa; })
                                         : DVar([](const scan::Row& x) { return (double)x.mee; }));
            TH1D* hC = H.Hist(3, Nm("bC")); StyleM(hC, 2, 21, 1.0);
            TH1D* h[5];
            const int cc[5] = {kColFake, kColDiff, kMagenta + 2, kOrange + 7, kGray + 2}, ls[5] = {1, 1, 2, 2, 2};
            const char* ln[5] = {"sim T3 fake lepton", "sim T4 different sources", "T3: same particle twice", "T3: ghost lepton", "T3: misidentified"};
            for (int q = 0; q < 5; ++q) { h[q] = S[q].Hist(Nm("bh"), 2); StyleL(h[q], cc[q], 3, ls[q]); }
            const double ymax = MaxOf({hC, h[0]}, false);
            TH1D* fr = Frame(hC, std::max(MinPos({h[1]}) * 0.3, ymax * 1e-4), ymax * 30,
                             v ? "M_{e^{+}e^{-}} [GeV/c^{2}]" : "OA(e^{+}e^{-}) [deg]",
                             v ? "d#sigma/dM_{ee} [nb/(GeV/c^{2})]" : "d#sigma/dOA [nb/deg]", 0.016, 0.46, true);
            fr->GetYaxis()->SetTitleOffset(1.2);
            fr->Draw("axis");
            for (int q = 0; q < 5; ++q) h[q]->Draw("hist same");
            hC->Draw("E1 P same");
            TLegend* lg = new TLegend(0.40, 0.58, 0.97, 0.89);
            lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.036);
            lg->AddEntry(hC, "data CB", "p");
            for (int q = 0; q < 5; ++q) lg->AddEntry(h[q], ln[q], "l");
            lg->Draw();
            Txt(0.14, 0.92, Form("(%c) %s, M_{#pi^{+}#pi^{-}e^{+}e^{-}} 0.30-1.10 GeV/c^{2}", v ? 'd' : 'b', v ? "M(e^{+}e^{-})" : "OA(e^{+}e^{-})"), 0.045);
        }
        // (c) ratio (T3 + T4) / CB per OA class
        {
            c->cd();
            TPad* pd = new TPad(Nm("cb"), "", 0, 0.0, 0.5, 0.465);
            pd->SetLeftMargin(0.14); pd->SetRightMargin(0.03); pd->SetTopMargin(0.10); pd->SetBottomMargin(0.14); pd->SetLogy();
            pd->Draw(); pd->cd();
            TH1D* fr = new TH1D(Nm("cf"), "", 5, 0, 5); fr->SetDirectory(nullptr);
            for (int k = 0; k < 5; ++k) fr->GetXaxis()->SetBinLabel(k + 1, kOAname[k]);
            fr->GetXaxis()->SetLabelSize(0.065); fr->GetYaxis()->SetLabelSize(0.050); fr->GetYaxis()->SetTitleSize(0.055);
            fr->GetXaxis()->SetTitleSize(0.055); fr->GetYaxis()->SetTitleOffset(1.15);
            fr->GetXaxis()->SetTitle("OA(e^{+}e^{-}) class [deg]"); fr->GetYaxis()->SetTitle("sim (L_{sim}) / data CB (L_{data})");
            fr->SetMinimum(0.005); fr->SetMaximum(3000); fr->Draw("axis");
            HLine(0, 5, 1.0, kGray + 2, 1);
            TLegend* lg = new TLegend(0.40, 0.68, 0.97, 0.89);
            lg->SetNColumns(3); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.034);
            const int rc[3] = {kBlue, kRed + 1, kBlack};
            for (int r = 0; r < 3; ++r) {
                TH1D* a = new TH1D(Nm("ca"), "", 5, 0, 5); a->SetDirectory(nullptr);
                TH1D* b3 = new TH1D(Nm("cb3"), "", 5, 0, 5); b3->SetDirectory(nullptr);
                TH1D* b4 = new TH1D(Nm("cb4"), "", 5, 0, 5); b4->SetDirectory(nullptr);
                for (int k = 0; k < 5; ++k) {
                    a->SetBinContent(k + 1, cbT34[r][k]); a->SetBinError(k + 1, cbT34e[r][k]);
                    b3->SetBinContent(k + 1, cbT3[r][k]); b4->SetBinContent(k + 1, cbT4[r][k]);
                }
                StyleM(a, rc[r], 20, 2.0); StyleM(b3, rc[r], 25, 1.8); StyleM(b4, rc[r], 26, 1.8);
                a->Draw("E1 P X0 same"); b3->Draw("P X0 same"); b4->Draw("P X0 same");
                lg->AddEntry(a, Form("(T3+T4)/CB %s", regs[r].name + 3), "p");
                lg->AddEntry(b3, "T3/CB", "p");
                lg->AddEntry(b4, "T4/CB", "p");
            }
            lg->Draw();
            Txt(0.14, 0.92, "(c) simulated information classes relative to the data CB, per OA class", 0.045);
        }
        Save(c, "datacomp_cb_content", true);
        delete c;
    }
}

}  // namespace dcmp

// ============================================================================
// Part 3: subsystems
// ============================================================================
namespace dcmp {

struct VD {
    std::string key, title, xt, yt;
    int nb; double lo, hi;
    Var s;      // simulation
    DFill d;    // data
};
inline double LepP(const scan::Row& r, int l) { return l == 0 ? r.ep_p : r.em_p; }
inline double LepTh(const scan::Row& r, int l) { return l == 0 ? r.ep_theta : r.em_theta; }

// variables; region 0 all M4, 1 M4 0.56-0.80, 2 bump window
inline std::vector<VD> PionVars(int reg) {
    const char* GM = "[GeV/c^{2}]"; (void)GM;
    std::vector<VD> v;
    const int h = reg == 2 ? 2 : 1;   // fewer bins in the bump window
    double mlo = 0.27, mhi = 1.07; int mnb = 40;
    if (reg == 1) { mhi = 0.79; mnb = 26; }
    if (reg == 2) { mhi = 0.53; mnb = 26; }
    v.push_back({"mpipi", "M(#pi^{+}#pi^{-})", "M_{#pi^{+}#pi^{-}} [GeV/c^{2}]", "nb/(GeV/c^{2})", mnb, mlo, mhi,
                 [](const Row& r) { return (double)r.mpipi; }, DVar([](const scan::Row& r) { return (double)r.mpipi; })});
    v.push_back({"p_pipi", "p(#pi^{+}#pi^{-})", "p_{#pi^{+}#pi^{-}} [MeV/c]", "nb/(MeV/c)", 50 / h, 0, 2500,
                 [](const Row& r) { return (double)r.p_pipi; }, DVar([](const scan::Row& r) { return (double)PairKinOf(r).p_pipi; })});
    v.push_back({"th_pipi", "#theta(#pi^{+}#pi^{-})", "#theta_{#pi^{+}#pi^{-}} [deg]", "nb/deg", 40 / h, 0, 80,
                 [](const Row& r) { return (double)r.th_pipi; }, DVar([](const scan::Row& r) { return (double)PairKinOf(r).th_pipi; })});
    v.push_back({"pt_pipi", "p_{T}(#pi^{+}#pi^{-})", "p_{T,#pi^{+}#pi^{-}} [MeV/c]", "nb/(MeV/c)", 40 / h, 0, 1000,
                 [](const Row& r) { return (double)r.pt_pipi; }, DVar([](const scan::Row& r) { return (double)PairKinOf(r).pt_pipi; })});
    v.push_back({"mm_pipi", "MM(#pi^{+}#pi^{-})", "MM_{#pi^{+}#pi^{-}} [GeV/c^{2}]", "nb/(GeV/c^{2})", 36 / h, 1.80, 3.24,
                 [](const Row& r) { return (double)r.mm_pipi; }, DVar([](const scan::Row& r) { return (double)r.mm_pipi; })});
    v.push_back({"pip_p", "p(#pi^{+})", "p_{#pi^{+}} [MeV/c]", "nb/(MeV/c)", 40 / h, 0, 2000,
                 [](const Row& r) { return (double)r.pip_p; }, DVar([](const scan::Row& r) { return (double)r.pip_p; })});
    v.push_back({"pip_theta", "#theta(#pi^{+})", "#theta_{#pi^{+}} [deg]", "nb/deg", 44 / h, 0, 88,
                 [](const Row& r) { return (double)r.pip_theta; }, DVar([](const scan::Row& r) { return (double)r.pip_theta; })});
    v.push_back({"pim_p", "p(#pi^{-})", "p_{#pi^{-}} [MeV/c]", "nb/(MeV/c)", 40 / h, 0, 2000,
                 [](const Row& r) { return (double)r.pim_p; }, DVar([](const scan::Row& r) { return (double)r.pim_p; })});
    v.push_back({"pim_theta", "#theta(#pi^{-})", "#theta_{#pi^{-}} [deg]", "nb/deg", 44 / h, 0, 88,
                 [](const Row& r) { return (double)r.pim_theta; }, DVar([](const scan::Row& r) { return (double)r.pim_theta; })});
    v.push_back({"oa_lab", "#alpha_{lab}(#pi^{+}#pi^{-}, e^{+}e^{-})", "#alpha_{lab} [deg]", "nb/deg", 50 / h, 0, 150,
                 [](const Row& r) { return (double)r.oa_lab; }, DVar([](const scan::Row& r) { return (double)r.oa_lab; })});
    v.push_back({"mm4", "MM(#pi^{+}#pi^{-}e^{+}e^{-})", "MM_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}]", "nb/(GeV/c^{2})", 40 / h, 1.40, 3.00,
                 [](const Row& r) { return (double)r.mm4; }, DVar([](const scan::Row& r) { return (double)r.mm4; })});
    v.push_back({"pt4", "p_{T}(#pi^{+}#pi^{-}e^{+}e^{-})", "p_{T,4} [MeV/c]", "nb/(MeV/c)", 40 / h, 0, 1000,
                 [](const Row& r) { return (double)r.pt4; }, DVar([](const scan::Row& r) { return (double)r.pt4; })});
    return v;
}
inline std::vector<VD> LeptonVars(int reg) {
    std::vector<VD> v;
    const int h = reg == 2 ? 2 : 1;
    v.push_back({"mee", "M(e^{+}e^{-})", "M_{e^{+}e^{-}} [GeV/c^{2}]", "nb/(GeV/c^{2})", 50 / h, 0, 0.20,
                 [](const Row& r) { return (double)r.mee; }, DVar([](const scan::Row& r) { return (double)r.mee; })});
    v.push_back({"oa", "OA(e^{+}e^{-})", "OA(e^{+}e^{-}) [deg]", "nb/deg", 40 / h, 0, 40,
                 [](const Row& r) { return (double)r.oa; }, DVar([](const scan::Row& r) { return (double)r.oa; })});
    v.push_back({"p_ee", "p(e^{+}e^{-})", "p_{e^{+}e^{-}} [MeV/c]", "nb/(MeV/c)", 50 / h, 0, 1500,
                 [](const Row& r) { return (double)r.p_ee; }, DVar([](const scan::Row& r) { return (double)PairKinOf(r).p_ee; })});
    v.push_back({"th_ee", "#theta(e^{+}e^{-})", "#theta_{e^{+}e^{-}} [deg]", "nb/deg", 44 / h, 0, 88,
                 [](const Row& r) { return (double)r.th_ee; }, DVar([](const scan::Row& r) { return (double)PairKinOf(r).th_ee; })});
    v.push_back({"pt_ee", "p_{T}(e^{+}e^{-})", "p_{T,e^{+}e^{-}} [MeV/c]", "nb/(MeV/c)", 40 / h, 0, 800,
                 [](const Row& r) { return (double)r.pt_ee; }, DVar([](const scan::Row& r) { return (double)PairKinOf(r).pt_ee; })});
    v.push_back({"mm_ee", "MM(e^{+}e^{-})", "MM_{e^{+}e^{-}} [GeV/c^{2}]", "nb/(GeV/c^{2})", 40 / h, 2.40, 3.40,
                 [](const Row& r) { return (double)r.mm_ee; }, DVar([](const scan::Row& r) { return (double)r.mm_ee; })});
    v.push_back({"ep_p", "p(e^{+})", "p_{e^{+}} [MeV/c]", "nb/(MeV/c)", 40 / h, 0, 1000,
                 [](const Row& r) { return (double)r.ep_p; }, DLep(LepP, 0)});
    v.push_back({"ep_theta", "#theta(e^{+})", "#theta_{e^{+}} [deg]", "nb/deg", 44 / h, 0, 88,
                 [](const Row& r) { return (double)r.ep_theta; }, DLep(LepTh, 0)});
    v.push_back({"em_p", "p(e^{-})", "p_{e^{-}} [MeV/c]", "nb/(MeV/c)", 40 / h, 0, 1000,
                 [](const Row& r) { return (double)r.em_p; }, DLep(LepP, 1)});
    v.push_back({"em_theta", "#theta(e^{-})", "#theta_{e^{-}} [deg]", "nb/deg", 44 / h, 0, 88,
                 [](const Row& r) { return (double)r.em_theta; }, DLep(LepTh, 1)});
    return v;
}

void Part3(const SimData& sd, const scan::Data& dd) {
    const int pres[2] = {kPre, kPreOA2};
    const char* ptag[2] = {"pre", "preOA2"};
    struct Reg { const char* name; const char* tag; double a, b; };
    const Reg regs[3] = {{"all M4", "allM4", -1e9, 1e9}, {"M4 0.56-0.80", "m4_056_080", 0.56, 0.80},
                         {"bump window M4 0.455-0.5225", "bump", scan::kBumpLo, scan::kBumpHi}};
    FILE* t = OpenTab("datacomp_subsys.txt");
    P2(t, "# subsystems of the pi+pi-e+e- hypothesis: data (N - CB) against the simulation (signal + physical background)\n");
    P2(t, "# data: L_data = %.1f nb^-1 (UNCONFIRMED), CB = 2 sqrt(N++ N--) per bin of the variable (integrated over the M4 region), group variances;\n", kLdata);
    P2(t, "#       single-lepton variables: like-sign rows enter with both leptons, weight 1/2 each\n");
    P2(t, "# sim: std + lep owned rows, mode A, w (PDG), L_sim = %.2f nb^-1. R = data/sim in M4 0.60-0.80 of the preset\n", kLsim);
    P2(t, "# columns: data/sim of the integrals in the histogram range; chi2/ndf with R fixed; free factor r and chi2/ndf (pure shape)\n");
    P2(t, "# chi2 bins: >= 10 data rows and >= 5 sim rows per bin\n");
    for (int ip = 0; ip < 2; ++ip) {
        const int p = pres[ip];
        double Rerr = 0;
        const double R = GetR(sd, dd, p, &Rerr);
        printf("Part3: %s R = %.4f +- %.4f\n", PresetName(p), R, Rerr);
        for (int ir = 0; ir < 3; ++ir) {
            const Reg& rg = regs[ir];
            const double a = rg.a, b = rg.b;
            const Sel sel = ir == 0 ? SelPreset(p) : SelAnd(SelPreset(p), SelM4(a, b));
            const scan::Sel ds0 = DataSelPreset(p);
            const scan::Sel dsel = ir == 0 ? ds0 : scan::Sel([ds0, a, b](const scan::Row& r) { return ds0(r) && r.m4 >= a && r.m4 < b; });
            P2(t, "\n## %s, %s (R = %.4f +- %.4f)\n", PresetName(p), rg.name, R, Rerr);
            P2(t, "%-10s %16s %14s %16s %14s\n", "variable", "data/sim (int.)", "chi2/ndf R", "free r", "chi2/ndf r");
            for (int grp = 0; grp < 2; ++grp) {
                const std::vector<VD> vars = grp == 0 ? PionVars(ir) : LeptonVars(ir);
                const int ncol = grp == 0 ? 4 : 5, nrow = grp == 0 ? 3 : 2;
                TCanvas* c = new TCanvas(Nm("cSub"), "", 2900, grp == 0 ? 2300 : 1750);
                c->cd();
                const double top = grp == 0 ? 0.925 : 0.90;
                const double ts = grp == 0 ? 0.0122 : 0.0150;      // header text size (canvas units)
                const double dy = grp == 0 ? 0.0178 : 0.0235;      // header line spacing
                const double sz = grp == 0 ? 0.0105 : 0.0130;      // axis label size (canvas units)
                Txt(0.005, 1 - 1.0 * dy, Form("%s subsystem: data against simulation; selection %s; %s; d#sigma/dx in nb per unit of x",
                                            grp == 0 ? "#pi^{+}#pi^{-} (and four-particle)" : "e^{+}e^{-}", PresetName(p), rg.name), ts * 1.15);
                Txt(0.005, 1 - 2.0 * dy, (DataLine() + "; like-sign single leptons: both, weight 1/2").c_str(), ts, 4);
                Txt(0.005, 1 - 3.0 * dy, SimLine().c_str(), ts, kColSim);
                Txt(0.005, 1 - 4.0 * dy,
                    Form("data: N black, CB red, N - CB blue; sim: solid absolute, dashed #times R = %.3f #pm %.3f (0.60-0.80), dotted: signal (T1) alone; ratio pad: (N - CB)/sim, grey line 1, dashed R; %s",
                         R, Rerr, "#eta #rightarrow #pi^{+}#pi^{-}e^{+}e^{-} NOT in sim"), ts, 1);
                for (size_t iv = 0; iv < vars.size(); ++iv) {
                    const VD& V = vars[iv];
                    Binning bb; bb.nbin = V.nb; bb.lo = V.lo; bb.hi = V.hi; bb.x = V.s;
                    const std::vector<WSpec> S = FillMany(sd, sel, {CompSim(), CompSignal()}, SimOpt(), bb);
                    const DH H = FillData(dd, dsel, V.nb, V.lo, V.hi, 1, V.d);
                    TH1D* hN = H.Hist(0, Nm("vN")); TH1D* hC = H.Hist(3, Nm("vC")); TH1D* hS = H.Hist(4, Nm("vS"));
                    StyleM(hN, 1, 20, 0.6); StyleM(hC, 2, 21, 0.55); StyleM(hS, 4, 20, 0.65);
                    TH1D* sS = S[0].Hist(Nm("vSim"), 2); StyleL(sS, kColSim, 3, 1);
                    TH1D* sR = (TH1D*)sS->Clone(Nm("vSR")); sR->Scale(R); StyleL(sR, kColSim, 2, 2);
                    TH1D* s1 = S[1].Hist(Nm("vT1"), 2); StyleL(s1, kColBumpConv, 2, 3);
                    // integrals and chi2
                    double di = 0, dv = 0;
                    for (int i = 0; i < V.nb; ++i) { di += H.sig[i]; dv += H.sigVar[i]; }
                    const double dI = di / kLdata, dE = std::sqrt(dv) / kLdata, sI = S[0].S / kLsim, sE = std::sqrt(S[0].V) / kLsim;
                    const double rI = sI > 0 ? dI / sI : 0, rIe = sI > 0 ? std::sqrt(dE * dE + rI * rI * sE * sE) / sI : 0;
                    std::vector<double> dRows = H.n[0];
                    const Chi ch = ShapeChi(hS, sS, dRows, S[0].n, V.lo, V.hi + 1e-9, R);
                    P2(t, "%-10s %7.3f +- %-6.3f %8.1f/%-4d %8.3f+-%-6.3f %8.1f/%-4d\n", V.key.c_str(), rI, rIe, ch.chiR, ch.nR, ch.rF,
                       ch.rFe, ch.chiF, ch.nF);
                    // pads
                    const int cx = iv % ncol, cy = (int)iv / ncol;
                    const double w = 1.0 / ncol, hcell = top / nrow;
                    const double x1 = cx * w, y2 = top - cy * hcell, y1 = y2 - hcell;
                    PP P = MkPair(c, x1, y1, x1 + w, y2, 0.33, true, 0.20, 0.09);
                    P.up->cd();
                    const double ymax = MaxOf({hN, sS, sR}, false);
                    const double ymin = std::max(MinPos({hS, sS}) * 0.3, ymax * 2e-4);
                    TH1D* fr = Frame(hN, ymin, ymax * 8, V.xt.c_str(), Form("d#sigma/dx [%s]", V.yt.c_str()), sz, P.hu, false);
                    fr->GetYaxis()->SetTitleOffset(1.45);
                    fr->Draw("axis");
                    sS->Draw("hist same"); sR->Draw("hist same"); s1->Draw("hist same");
                    hN->Draw("E1 P same"); hC->Draw("E1 P same"); hS->Draw("E1 P same");
                    fr->Draw("axis same");
                    Txt(0.21, 0.93, Form("%s", V.title.c_str()), 1.15 * sz / P.hu);
                    Txt(0.96, 0.93, Form("#chi^{2}/ndf (R): %.0f/%d; free r %.2f: %.0f/%d", ch.chiR, ch.nR, ch.rF, ch.chiF, ch.nF), 0.9 * sz / P.hu, 1, 31);
                    P.lo->cd();
                    TH1D* rr = Ratio(hS, sS, Nm("vR")); StyleM(rr, 4, 20, 0.7);
                    TH1D* fl = Frame(hS, 0.0, 2.5, V.xt.c_str(), "data/sim", sz, P.hl, true);
                    fl->GetYaxis()->SetTitleOffset(1.45 * P.hl / P.hu);
                    fl->GetXaxis()->SetTitleOffset(1.0);
                    fl->Draw("axis");
                    HLine(V.lo, V.hi, 1.0, kGray + 2, 1);
                    HLine(V.lo, V.hi, R, kColSim, 2, 2);
                    rr->Draw("E1 P same");
                }
                Save(c, Form("datacomp_subsys_%s_%s_%s", grp == 0 ? "pipi" : "ee", ptag[ip], rg.tag), ip == 0 && ir == 1);
                delete c;
            }
        }
    }
    if (t) fclose(t);
}

}  // namespace dcmp

// ============================================================================
void sim_datacomp(int part = 0) {
    using namespace dcmp;
    TStopwatch sw;
    TH1::AddDirectory(kFALSE);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetEndErrorSize(0);
    gStyle->SetLineScalePS(1.5);
    if (part < 0) { printf("sim_datacomp: compiled only\n"); return; }
    SimData sd;
    if (!Load(sd, kCompact, (1 << kStdUS) | (1 << kLepUS))) return;
    scan::Data dd;
    if (!scan::Load(dd, scan::kCompact)) return;
    if (part == 0 || part == 1) Part1(sd, dd);
    if (part == 0 || part == 2) Part2(sd, dd);
    if (part == 0 || part == 3) Part3(sd, dd);
    printf("sim_datacomp(%d): total %.0f s\n", part, sw.RealTime());
}
