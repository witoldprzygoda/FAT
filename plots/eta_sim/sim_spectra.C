// ============================================================================
// sim_spectra.C  --  simulated spectra in absolute scale and the eta structures
// ============================================================================
// HADES p+p, T = 4.53 GeV, pi+ pi- e+ e- (e+e- hypothesis). SMASH std + lep,
// rows owned in partition mode A (mode B as a variant), weight w (PDG BR;
// w_db as a variant). Cross section within the HADES acceptance after the
// selection, no efficiency correction:
//     sigma per bin = sum w / L_sim (L_sim = 221.05 nb^-1 per set; std and lep
//     owned rows added with the same L_sim), dsigma/dx = that / bin width.
// Data in the same units: (N - CB) / L_data / bin width, L_data = 6080.5 nb^-1
// (producers' note, UNCONFIRMED), CB = 2 sqrt(N++ N--) (no k factor).
// Simulated spectrum = signal (template 1) + physical background (template 2);
// fake lepton (3) and different sources (4) only as dashed information curves.
// eta -> pi+pi-e+e- (direct) is NOT in the simulation (eta set missing).
//
// Parts (argument "parts", bit mask; default all):
//   1  M4 per selection: components, stacked and separate, log 0.3-1.1 and
//      linear 0.40-0.70, mode A/B and w/w_db of the sum, data in the same
//      units (figs spectra_m4_<sel>, spectra_m4data_<sel>, spectra_m4_presets)
//   2  subsystems per selection (figs spectra_sub_<sel>_{pipi,ee,tracks,four})
//   4  eta structures: reco and true M4, peak fit, bump edge per M(pipi)
//      slice, resolutions, mass scales (figs spectra_eta_*, table spectra_eta.txt)
//   8  window and region integrals, templates with N_eff (tables
//      spectra_windows.txt, spectra_regions.txt, spectra_templates.txt)
//
// Statistics: per bin S = sum w, V = sum over event groups of (group sum)^2
// (simeng::FillMany), N_eff = S^2 / V; all sums in double (TH1D, double
// accumulators). Fits: weighted least squares with variances q * model
// (q = V / S of the component, iterated), never a Poisson likelihood.
//
// Run from /home/przygoda/HADES/FAT/FAT:
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/aclic_spectra", kTRUE)' 'plots/eta_sim/sim_spectra.C+'
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <map>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "TBox.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TError.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMarker.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TSystem.h"

#include "SimEngine.h"

namespace ssp {
using namespace simeng;

// ============================================================================
// selections, components, styles
// ============================================================================
const int kNP = 6;
const int kPres[kNP] = {kPre, kPreOA2, kPrePiVtx, kSlides, kPeakNom, kBumpNom};
const char* const kPTag[kNP] = {"pre", "preOA2", "prePiVtx", "slides", "peakNom", "bumpNom"};

enum { cBD = 0, cBC, cPC, cP7, cP8, cP1, cP0, cP3, cP4, cP6, cP25, cT2, cT1, cSum, cFake, cDiff, kNC };
struct CDef { const char* tag; const char* leg; const char* txt; int line, fill; };
const CDef kCD[kNC] = {
    {"bump_dalitz", "signal: bump, #pi^{0} Dalitz", "signal: bump, pi0 Dalitz", 807, 807},
    {"bump_conv", "signal: bump, conversion", "signal: bump, conversion", 617, 617},
    {"peak_conv", "signal: peak, conversion", "signal: peak, conversion", 418, 418},
    {"bkg_primary", "phys. bkg: both #pi primary", "phys. bkg: both pions primary (P7)", 921, 920},
    {"bkg_mixed", "phys. bkg: #pi mixed / secondary", "phys. bkg: pions mixed/secondary (P8)", 803, 795},
    {"bkg_pimu", "phys. bkg: #pi #rightarrow #mu", "phys. bkg: pi -> mu (P1)", 402, 390},
    {"bkg_misid", "phys. bkg: misidentified #pi", "phys. bkg: misidentified pion (P0)", 434, 423},
    {"bkg_omega", "phys. bkg: #pi^{+}#pi^{-} of one #omega", "phys. bkg: pions of one omega (P3)", 881, 871},
    {"bkg_rho", "phys. bkg: #pi^{+}#pi^{-} of one #rho^{0}", "phys. bkg: pions of one rho0 (P4)", 902, 893},
    {"bkg_k0s", "phys. bkg: #pi^{+}#pi^{-} of one K^{0}_{S}", "phys. bkg: pions of one K0S (P6)", 824, 813},
    {"bkg_eta_phi", "phys. bkg: #pi^{+}#pi^{-} of one #eta or #phi", "phys. bkg: pions of one eta (ee not of it) or phi (P2,P5)", 801, 797},
    {"bkg_total", "physical background, total", "physical background, total (T2)", 923, 920},
    {"signal_total", "signal, total", "signal, total (T1)", 807, 807},
    {"sim_sum", "SUM = signal + phys. bkg", "SUM: signal + physical background (T1+T2)", 843, 843},
    {"fake_lepton", "fake lepton (not in the sum)", "fake lepton (T3, not in the sum)", 861, 861},
    {"diff_sources", "different sources (not in the sum)", "different sources (T4, not in the sum)", 403, 403}};
// colours (set by InitColors): signals keep the engine colours; background lines dark, fills light
int gLine[kNC], gFill[kNC];
void InitColors() {
    const char* L[kNC] = {nullptr, nullptr, nullptr, "#808080", "#8c510a", "#4d6a8a", "#00838f", "#6a3d9a", "#c51b7d", "#1b7837",
                          "#c08a00", "#404040", nullptr, nullptr, nullptr, nullptr};
    const char* F[kNC] = {nullptr, nullptr, nullptr, "#d9d9d9", "#dfc27d", "#c6d3e3", "#a8e0e8", "#cab2d6", "#f4b6d8", "#a6dba0",
                          "#fee391", "#d9d9d9", nullptr, nullptr, nullptr, nullptr};
    for (int k = 0; k < kNC; ++k) {
        gLine[k] = L[k] ? TColor::GetColor(L[k]) : kCD[k].line;
        gFill[k] = F[k] ? TColor::GetColor(F[k]) : kCD[k].fill;
    }
}
const int kNBkg = 8;
const int kBkgOrder[kNBkg] = {cP7, cP8, cP1, cP0, cP3, cP4, cP6, cP25};
const int kSigOrder[3] = {cBD, cBC, cPC};

inline std::vector<Comp> Comps() {
    std::vector<Comp> c(kNC);
    c[cBD] = CompBumpDalitz(); c[cBC] = CompBumpConv(); c[cPC] = CompPeakConv();
    c[cP7] = CompBkgPcls(7); c[cP8] = CompBkgPcls(8); c[cP1] = CompBkgPcls(1); c[cP0] = CompBkgPcls(0);
    c[cP3] = CompBkgPcls(3); c[cP4] = CompBkgPcls(4); c[cP6] = CompBkgPcls(6);
    c[cP25] = {"phys. bkg, P2 same eta / P5 same phi / P9",
               [](const Row& r) { return r.tmpl == 2 && (r.pcls == 2 || r.pcls == 5 || r.pcls == 9); }};
    c[cT2] = CompPhysBkg(); c[cT1] = CompSignal(); c[cSum] = CompSim(); c[cFake] = CompFake(); c[cDiff] = CompDiff();
    return c;
}

// data values of the data study (reports engine.md, scan_bump.md; not recomputed here)
const double kDataPeak = 0.5439, kDataPeakErr = 0.00049, kDataSigma = 0.0054, kDataSigmaErr = 0.0006;
const double kDataScale = 0.9928, kDataScaleErr = 0.0009, kDataScalePiPi = 0.9970;
const int kNDataEdge = 4;   // slices 0.30-0.34, 0.34-0.37, 0.37-0.39, 0.39-0.40 ("no cut")
const double kDataEdgeLo[kNDataEdge] = {0.30, 0.34, 0.37, 0.39}, kDataEdgeHi[kNDataEdge] = {0.34, 0.37, 0.39, 0.40};
const double kDataEdge[kNDataEdge] = {0.5158, 0.5116, 0.5048, 0.5004};
const double kDataEdgeF[kNDataEdge] = {0.9912, 0.9932, 0.9927, 0.9972};
const double kDataEdgeFm[kNDataEdge] = {0.0030, 0.0025, 0.0020, 0.0047}, kDataEdgeFp[kNDataEdge] = {0.0027, 0.0022, 0.0022, 0.0040};

// windows and regions in M4 [GeV/c2]
const int kNReg = 7;
const double kRegLo[kNReg] = {0.30, scan::kBumpLo, scan::kBumpHi, scan::kPeakLo, scan::kPeakHi, 0.58, 0.70};
const double kRegHi[kNReg] = {scan::kBumpLo, scan::kBumpHi, scan::kPeakLo, scan::kPeakHi, 0.58, 0.70, 1.10};
const char* const kRegName[kNReg] = {"below the bump", "bump window", "between bump and peak", "peak window",
                                     "just above the peak", "above the peak", "high mass"};
const int kRegNorm = 5;   // renormalisation region of the figures: 0.58-0.70 ("above the peak")

std::string FigPath(const std::string& n) { return std::string(kFigDir) + "/" + n; }
std::string TabPath(const std::string& n) { return std::string(kTabDir) + "/" + n; }

// ---------------- drawing helpers ----------------
void StyleLine(TH1D* h, int col, int lw, int ls) {
    h->SetLineColor(col); h->SetLineWidth(lw); h->SetLineStyle(ls); h->SetFillStyle(0);
    h->SetMarkerColor(col); h->SetMarkerSize(0);
}
void StyleFill(TH1D* h, int line, int fill) {
    h->SetLineColor(line); h->SetLineWidth(1); h->SetLineStyle(1); h->SetFillColor(fill); h->SetFillStyle(1001);
    h->SetMarkerSize(0);
}
void StyleComp(TH1D* h, int k, double scale = 1.0) {
    if (k == cSum) StyleLine(h, gLine[k], (int)std::lround(4 * scale), 1);
    else if (k == cFake || k == cDiff) StyleLine(h, gLine[k], (int)std::lround(3 * scale), 2);
    else if (k == cBD || k == cBC || k == cPC) StyleLine(h, gLine[k], (int)std::lround(3 * scale), 1);
    else if (k == cT2) StyleLine(h, gLine[k], (int)std::lround(3 * scale), 1);
    else StyleLine(h, gLine[k], std::max(1, (int)std::lround(2 * scale)), 1);
}
void StyleData(TH1D* h, int col, int mst, double ms = 1.0) {
    h->SetMarkerStyle(mst); h->SetMarkerSize(ms); h->SetMarkerColor(col); h->SetLineColor(col); h->SetLineWidth(1);
}
void Axes(TH1* h, const char* xt, const char* yt, double ts = 0.055, double yoff = 1.25) {
    h->GetXaxis()->SetTitle(xt); h->GetYaxis()->SetTitle(yt);
    h->GetXaxis()->SetTitleSize(ts); h->GetYaxis()->SetTitleSize(ts);
    h->GetXaxis()->SetLabelSize(ts * 0.87); h->GetYaxis()->SetLabelSize(ts * 0.87);
    h->GetXaxis()->SetTitleOffset(1.0); h->GetYaxis()->SetTitleOffset(yoff);
    h->GetYaxis()->SetMaxDigits(4);
}
void Txt(double x, double y, const char* s, double size, int col = 1, int align = 11) {
    TLatex t; t.SetNDC(); t.SetTextSize(size); t.SetTextColor(col); t.SetTextAlign(align); t.SetTextFont(42);
    t.DrawLatex(x, y, s);
}
TPad* MkPad(const char* name, double x1, double y1, double x2, double y2, bool logy, double lm, double rm, double tm,
            double bm) {
    TPad* p = new TPad(name, "", x1, y1, x2, y2);
    p->SetLeftMargin(lm); p->SetRightMargin(rm); p->SetTopMargin(tm); p->SetBottomMargin(bm);
    p->SetLogy(logy);
    p->SetTicks(1, 1);
    p->Draw();
    p->cd();
    return p;
}
TH1D* Frame(const TH1D* like, const char* name, double ymin, double ymax) {
    TH1D* f = (TH1D*)like->Clone(name);
    f->Reset();
    f->SetMinimum(ymin); f->SetMaximum(ymax);
    f->SetLineColor(0);
    return f;
}
// cumulative histograms of the components in `order` (first = bottom); returned bottom first
std::vector<TH1D*> Cumul(TH1D* const h[], const std::vector<int>& order, const std::string& pre, double scale = 1.0) {
    std::vector<TH1D*> out;
    TH1D* prev = nullptr;
    for (int k : order) {
        TH1D* c = (TH1D*)h[k]->Clone(Form("%s_cum%d", pre.c_str(), k));
        c->Scale(scale);
        if (prev) c->Add(prev);
        StyleFill(c, gLine[k], gFill[k]);
        out.push_back(c);
        prev = c;
    }
    return out;
}
void DrawStack(const std::vector<TH1D*>& st) { for (int i = (int)st.size() - 1; i >= 0; --i) st[i]->Draw("hist same"); }
double MaxOf(const TH1D* h, bool withErr = false) {
    double m = 0;
    for (int i = 1; i <= h->GetNbinsX(); ++i) m = std::max(m, h->GetBinContent(i) + (withErr ? h->GetBinError(i) : 0.));
    return m;
}
double MinPos(const TH1D* h) {
    double m = 1e300;
    for (int i = 1; i <= h->GetNbinsX(); ++i) if (h->GetBinContent(i) > 0) m = std::min(m, h->GetBinContent(i));
    return m;
}
void VLine(double x, double y1, double y2, int col, int ls, int lw = 2) {
    TLine* l = new TLine(x, y1, x, y2); l->SetLineColor(col); l->SetLineStyle(ls); l->SetLineWidth(lw); l->Draw();
}
void HLine(double x1, double x2, double y, int col, int ls, int lw = 2) {
    TLine* l = new TLine(x1, y, x2, y); l->SetLineColor(col); l->SetLineStyle(ls); l->SetLineWidth(lw); l->Draw();
}
// a model given per bin (expected sum w per bin, x0 = bin low edge) drawn as a smooth line in [lo, hi), in dsigma/dx units
TGraph* Curve(const std::function<double(double, double)>& perBin, double bw, double lo, double hi, double conv, int col, int lw = 3) {
    TGraph* g = new TGraph();
    const int n = 200;
    for (int i = 0; i <= n; ++i) {
        const double x = lo + (hi - lo) * i / n;
        g->SetPoint(i, x, perBin(x - 0.5 * bw, bw) * conv);
    }
    g->SetLineColor(col); g->SetLineWidth(lw);
    return g;
}
void Mark(double x, double y, int col, int st, double sz = 1.6) {
    TMarker* m = new TMarker(x, y, st); m->SetMarkerColor(col); m->SetMarkerSize(sz); m->Draw();
}
TH1D* Ratio(const TH1D* num, const TH1D* den, const char* name) {
    TH1D* r = (TH1D*)num->Clone(name);
    for (int i = 1; i <= r->GetNbinsX(); ++i) {
        const double a = num->GetBinContent(i), ea = num->GetBinError(i), b = den->GetBinContent(i), eb = den->GetBinError(i);
        if (b > 0) {
            r->SetBinContent(i, a / b);
            r->SetBinError(i, std::sqrt(ea * ea / (b * b) + a * a * eb * eb / (b * b * b * b)));
        } else { r->SetBinContent(i, -10); r->SetBinError(i, 0); }
    }
    return r;
}

const std::string kSimLine = Form("simulation: SMASH std + lep, owned rows, mode A, w (PDG BR); #sigma = #Sigma w / L_{sim}, L_{sim} = %.2f nb^{-1}; HADES acceptance, no efficiency correction", kLsim);
const std::string kDataLine = Form("data: (N - CB) / L_{data} / bin width, %s; CB = 2#sqrt{N_{++}N_{--}} (no k factor); no efficiency correction", kLdataNote);
const char* const kDirect = "#eta #rightarrow #pi^{+}#pi^{-}e^{+}e^{-} (direct decay) is NOT in the simulation (eta set not available)";
const char* const kM4t = "M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}]";
const char* const kM4y = "d#sigma/dM_{#pi^{+}#pi^{-}e^{+}e^{-}} [nb/(GeV/c^{2})]";

// ============================================================================
// weighted least-squares fits (Minuit2)
// ============================================================================
struct MinRes { bool ok = false; std::vector<double> p, e; double fval = 0; };
MinRes Minimize(const std::function<double(const double*)>& f, const std::vector<double>& p0, const std::vector<double>& st,
                const std::vector<double>& lo, const std::vector<double>& hi, const std::vector<bool>& fix) {
    MinRes r;
    std::unique_ptr<ROOT::Math::Minimizer> m(ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
    if (!m) return r;
    m->SetMaxFunctionCalls(200000); m->SetMaxIterations(200000); m->SetTolerance(0.01);
    m->SetPrintLevel(-1); m->SetErrorDef(1.0); m->SetStrategy(1);
    ROOT::Math::Functor fn(f, (unsigned)p0.size());
    m->SetFunction(fn);
    for (size_t i = 0; i < p0.size(); ++i) {
        const std::string n = Form("p%zu", i);
        // a double-sided limit with a huge range freezes the parameter (sin transformation at its edge):
        // upper limits >= 1e11 mean "lower limit only"
        if (!fix.empty() && fix[i]) m->SetFixedVariable(i, n, p0[i]);
        else if (lo[i] < hi[i] && hi[i] < 1e11) m->SetLimitedVariable(i, n, p0[i], st[i], lo[i], hi[i]);
        else if (lo[i] < hi[i]) m->SetLowerLimitedVariable(i, n, p0[i], st[i], lo[i]);
        else m->SetVariable(i, n, p0[i], st[i]);
    }
    r.ok = m->Minimize();
    if (r.ok) r.ok = m->Hesse();
    const double* x = m->X();
    const double* e = m->Errors();
    r.p.assign(x, x + p0.size());
    r.e.assign(e, e + p0.size());
    r.fval = m->MinValue();
    return r;
}

// bins of a weighted spectrum inside [lo, hi): sum w, group variance; q = V / S of the whole spectrum
struct BinData { std::vector<double> x, y, v; double bw = 0, q = 1; };
BinData Bins(const WSpec& s, double lo, double hi) {
    BinData b; b.bw = s.bw;
    b.q = (s.S > 0 && s.V > 0) ? s.V / s.S : 1.0;
    for (int i = 0; i < s.nbin; ++i) {
        if (s.x(i) < lo || s.x(i) >= hi) continue;
        b.x.push_back(s.x(i)); b.y.push_back(s.s[i]); b.v.push_back(s.v[i]);
    }
    return b;
}
typedef std::function<double(double, double, const double*)> Model;   // (bin low edge, bin width, params) -> expected sum w
struct WFit { bool ok = false; std::vector<double> p, e; double chi2 = 0; int ndf = 0, npts = 0; };
// iterated weighted least squares: variance of bin i = q * max(model_i, q/2) (first pass: observed variance, floor q^2)
WFit FitW(const BinData& b, const Model& model, std::vector<double> p, const std::vector<double>& st,
          const std::vector<double>& lo, const std::vector<double>& hi, const std::vector<bool>& fix = {}) {
    WFit r;
    const int n = (int)b.x.size();
    r.npts = n;
    int nfree = 0;
    for (size_t i = 0; i < p.size(); ++i) if (fix.empty() || !fix[i]) ++nfree;
    if (n <= nfree) return r;
    std::vector<double> var(n);
    for (int i = 0; i < n; ++i) var[i] = std::max(b.v[i], b.q * b.q);
    for (int it = 0; it < 5; ++it) {
        auto chi2 = [&](const double* pp) {
            double c = 0;
            for (int i = 0; i < n; ++i) {
                const double d = b.y[i] - model(b.x[i] - 0.5 * b.bw, b.bw, pp);
                c += d * d / var[i];
            }
            return c;
        };
        const MinRes m = Minimize(chi2, p, st, lo, hi, fix);
        r.ok = m.ok;
        p = m.p; r.e = m.e;
        for (int i = 0; i < n; ++i) var[i] = b.q * std::max(model(b.x[i] - 0.5 * b.bw, b.bw, p.data()), 0.5 * b.q);
    }
    r.p = p;
    double c = 0;
    for (int i = 0; i < n; ++i) {
        const double d = b.y[i] - model(b.x[i] - 0.5 * b.bw, b.bw, p.data());
        c += d * d / var[i];
    }
    r.chi2 = c; r.ndf = n - nfree;
    return r;
}
// Gaussian (bin integrated): p = {A (sum w under the Gaussian), mu, sigma}
double GaussBin(double x0, double bw, const double* p) {
    const double s = std::max(std::fabs(p[2]), 1e-7) * std::sqrt(2.0);
    return p[0] * 0.5 * (std::erfc((x0 - p[1]) / s) - std::erfc((x0 + bw - p[1]) / s));
}
struct GRes { bool ok = false; double A = 0, mu = 0, sig = 0, eA = 0, emu = 0, esig = 0, chi2 = 0, lo = 0, hi = 0; int ndf = 0; };
// Gaussian core: range mu -+ nsig sigma, iterated
GRes FitGaussCore(const WSpec& s, double mu0, double sig0, double nsig = 2.0) {
    GRes g; g.mu = mu0; g.sig = sig0;
    for (int it = 0; it < 6; ++it) {
        const double lo = std::min(g.mu - nsig * g.sig, g.mu - 3 * s.bw), hi = std::max(g.mu + nsig * g.sig, g.mu + 3 * s.bw);
        const BinData b = Bins(s, lo, hi);
        double a0 = 0; for (double y : b.y) a0 += y;
        const WFit f = FitW(b, GaussBin, {std::max(a0, 1e-9) * 1.2, g.mu, g.sig}, {std::max(a0, 1e-9) * 0.1, g.sig * 0.1, g.sig * 0.1},
                            {0, mu0 - 0.05, 1e-5}, {1e12, mu0 + 0.05, 0.05});
        if (f.p.size() < 3 || !(f.p[2] == f.p[2])) return g;
        g.ok = f.ok; g.A = f.p[0]; g.mu = f.p[1]; g.sig = std::fabs(f.p[2]);
        g.eA = f.e[0]; g.emu = f.e[1]; g.esig = f.e[2]; g.chi2 = f.chi2; g.ndf = f.ndf; g.lo = lo; g.hi = hi;
    }
    return g;
}
// weighted mean and rms of a spectrum inside [lo, hi) (bin centres)
void MeanRms(const WSpec& s, double lo, double hi, double& mean, double& rms) {
    double sw = 0, sx = 0, sx2 = 0;
    for (int i = 0; i < s.nbin; ++i) {
        if (s.x(i) < lo || s.x(i) >= hi) continue;
        sw += s.s[i]; sx += s.s[i] * s.x(i); sx2 += s.s[i] * s.x(i) * s.x(i);
    }
    mean = sw > 0 ? sx / sw : 0;
    rms = sw > 0 ? std::sqrt(std::max(0., sx2 / sw - mean * mean)) : 0;
}

// ---- bump edge: step averaged over a M(pipi) slice (as the data study, scan_bump.C) ----
// expected sum w per bin = A (1 + B (m - ec) / 0.04) T(m), T = mean over 21 points of the slice of
// 0.5 erfc((m - k M4max(mpipi / sPiPi)) / (sqrt2 s)); p = {A, B, k, s}; ec = M4max(centre / sPiPi)
struct EdgeModel {
    std::vector<double> E;   // M4max at the 21 slice points
    double ec = 0;
    EdgeModel() {}
    // aligned variable u = M4 / M4max(M(pipi)): one edge at u = k
    static EdgeModel Aligned() { EdgeModel m; m.E.push_back(1.0); m.ec = 1.0; return m; }
    EdgeModel(double mpLo, double mpHi, double sPiPi) {
        for (int j = 0; j < 21; ++j) {
            const double mp = mpLo + (j + 0.5) * (mpHi - mpLo) / 21.0;
            const double e = scan::M4max(mp / sPiPi);
            if (e > 0) E.push_back(e);
        }
        ec = scan::M4max(0.5 * (mpLo + mpHi) / sPiPi);
    }
    double T(double m, double k, double s) const {
        double t = 0;
        const double q = 1.0 / (std::sqrt(2.0) * std::max(s, 1e-7));
        for (double e : E) t += 0.5 * std::erfc((m - k * e) * q);
        return E.empty() ? 0. : t / E.size();
    }
    double operator()(double x0, double bw, const double* p) const {   // 5-point bin average
        double f = 0;
        for (int i = 0; i < 5; ++i) {
            const double m = x0 + (i + 0.5) * bw / 5;
            f += (1.0 + p[1] * (m - ec) / 0.04) * T(m, p[2], p[3]);
        }
        return p[0] * f / 5;
    }
};
struct ERes { bool ok = false; double k = 0, ek = 0, s = 0, es = 0, A = 0, B = 0, chi2 = 0, ec = 0, edge = 0, eedge = 0, S = 0, Neff = 0; int ndf = 0; double lo = 0, hi = 0; };
// The amplitudes A, A B enter linearly: for fixed (k, s) they follow from the normal equations. k and s are
// scanned on a grid (coarse, then fine around the minimum), as the data study profiles chi2 over the scale.
// Variance of a bin: first pass observed group variance (floor q^2), then q * model of the best point (3 passes).
// Errors of k and s: delta chi2 = 1 of the profile (minimum over the other grid variable).
ERes FitEdge(const WSpec& sp, const EdgeModel& em, double k0, double sMin, double sMax, double below,
             double above, bool fixS = false, double sFix = 0) {
    ERes r;
    r.ec = em.ec;
    r.lo = k0 * em.ec - below; r.hi = k0 * em.ec + above;
    const BinData b = Bins(sp, r.lo, r.hi);
    const int n = (int)b.x.size();
    if (n < 6) return r;
    std::vector<double> var(n);
    for (int i = 0; i < n; ++i) var[i] = std::max(b.v[i], b.q * b.q);
    // chi2 at (k, s) with the linear amplitudes; also returns A and A*B
    auto eval = [&](double k, double s, double& a1, double& a2, std::vector<double>* model) {
        double s11 = 0, s12 = 0, s22 = 0, y1 = 0, y2 = 0;
        std::vector<double> g1(n), g2(n);
        for (int i = 0; i < n; ++i) {
            double t1 = 0, t2 = 0;
            const double x0 = b.x[i] - 0.5 * b.bw;
            for (int j = 0; j < 5; ++j) {
                const double m = x0 + (j + 0.5) * b.bw / 5;
                const double t = em.T(m, k, s);
                t1 += t; t2 += t * (m - em.ec) / 0.04;
            }
            g1[i] = t1 / 5; g2[i] = t2 / 5;
            const double w = 1.0 / var[i];
            s11 += w * g1[i] * g1[i]; s12 += w * g1[i] * g2[i]; s22 += w * g2[i] * g2[i];
            y1 += w * b.y[i] * g1[i]; y2 += w * b.y[i] * g2[i];
        }
        const double det = s11 * s22 - s12 * s12;
        if (!(det > 0)) { a1 = a2 = 0; return 1e30; }
        a1 = (y1 * s22 - y2 * s12) / det; a2 = (y2 * s11 - y1 * s12) / det;
        double c = 0;
        for (int i = 0; i < n; ++i) {
            const double f = a1 * g1[i] + a2 * g2[i];
            if (model) (*model)[i] = f;
            c += (b.y[i] - f) * (b.y[i] - f) / var[i];
        }
        return a1 > 0 ? c : c + 1e6;
    };
    double kb = k0, sb = fixS ? sFix : 0.5 * (sMin + sMax), cb = 1e300, a1b = 0, a2b = 0;
    std::vector<double> ks, ss, prof_k, prof_s;
    for (int pass = 0; pass < 3; ++pass) {
        // coarse grid
        cb = 1e300;
        const int nkC = 61, nsC = fixS ? 1 : 30;
        for (int ik = 0; ik < nkC; ++ik) {
            const double k = k0 - 0.03 + 0.06 * ik / (nkC - 1);
            for (int is = 0; is < nsC; ++is) {
                const double s = fixS ? sFix : sMin * std::pow(sMax / sMin, is / (double)(nsC - 1));
                double a1, a2;
                const double c = eval(k, s, a1, a2, nullptr);
                if (c < cb) { cb = c; kb = k; sb = s; }
            }
        }
        // fine grid around the coarse minimum
        const int nkF = 121, nsF = fixS ? 1 : 41;
        const double kLo = kb - 0.012, kHi = kb + 0.012;
        const double sLo = fixS ? sFix : std::max(sMin, sb / 2.5), sHi = fixS ? sFix : std::min(sMax, sb * 2.5);
        ks.assign(nkF, 0); prof_k.assign(nkF, 1e300); ss.assign(nsF, 0); prof_s.assign(nsF, 1e300);
        for (int ik = 0; ik < nkF; ++ik) {
            const double k = kLo + (kHi - kLo) * ik / (nkF - 1);
            ks[ik] = k;
            for (int is = 0; is < nsF; ++is) {
                const double s = nsF == 1 ? sLo : sLo * std::pow(sHi / sLo, is / (double)(nsF - 1));
                ss[is] = s;
                double a1, a2;
                const double c = eval(k, s, a1, a2, nullptr);
                prof_k[ik] = std::min(prof_k[ik], c);
                prof_s[is] = std::min(prof_s[is], c);
                if (c < cb) { cb = c; kb = k; sb = s; a1b = a1; a2b = a2; }
            }
        }
        // new variances from the best model
        std::vector<double> model(n);
        eval(kb, sb, a1b, a2b, &model);
        for (int i = 0; i < n; ++i) var[i] = b.q * std::max(model[i], 0.5 * b.q);
    }
    // final chi2 and profiles with the final variances
    {
        double a1, a2;
        cb = eval(kb, sb, a1, a2, nullptr); a1b = a1; a2b = a2;
        for (size_t ik = 0; ik < ks.size(); ++ik) { prof_k[ik] = 1e300; }
        for (size_t is = 0; is < ss.size(); ++is) { prof_s[is] = 1e300; }
        for (size_t ik = 0; ik < ks.size(); ++ik)
            for (size_t is = 0; is < ss.size(); ++is) {
                const double c = eval(ks[ik], ss[is], a1, a2, nullptr);
                prof_k[ik] = std::min(prof_k[ik], c);
                prof_s[is] = std::min(prof_s[is], c);
                if (c < cb) { cb = c; kb = ks[ik]; sb = ss[is]; a1b = a1; a2b = a2; }
            }
    }
    auto err = [&](const std::vector<double>& x, const std::vector<double>& pr, double xb) {
        double lo = xb, hi = xb;
        for (size_t i = 0; i < x.size(); ++i) if (pr[i] < cb + 1.0) { lo = std::min(lo, x[i]); hi = std::max(hi, x[i]); }
        const bool atEdge = (lo <= x.front() + 1e-12) || (hi >= x.back() - 1e-12);
        return atEdge ? -0.5 * (hi - lo) : 0.5 * (hi - lo);   // negative: interval reaches the grid edge
    };
    r.k = kb; r.ek = err(ks, prof_k, kb); r.s = sb; r.es = fixS ? 0 : err(ss, prof_s, sb);
    r.A = a1b; r.B = a1b > 0 ? a2b / a1b : 0;
    r.chi2 = cb; r.ndf = n - (fixS ? 3 : 4);
    r.edge = r.k * em.ec; r.eedge = std::fabs(r.ek) * em.ec;
    r.ok = a1b > 0 && r.k > k0 - 0.0299 && r.k < k0 + 0.0299;
    return r;
}

// ============================================================================
// data: pair quantities per data row (scan::Row has no pair columns)
// ============================================================================
struct DataPK {
    const scan::Row* b[3] = {nullptr, nullptr, nullptr};
    size_t n[3] = {0, 0, 0};
    std::vector<PairKin> pk[3];
    PairKin nan;
    void Build(const scan::Data& d) {
        const float q = std::numeric_limits<float>::quiet_NaN();
        nan = PairKin{q, q, q, q, q, q, q, q};
        for (int k = 0; k < 3; ++k) {
            b[k] = d.s[k].data(); n[k] = d.s[k].size();
            pk[k].resize(n[k]);
            for (size_t i = 0; i < n[k]; ++i) pk[k][i] = d.s[k][i].cor_ok ? PairKinOf(d.s[k][i]) : nan;
        }
    }
    const PairKin& Of(const scan::Row& r) const {
        for (int k = 0; k < 3; ++k)
            if (n[k] && &r >= b[k] && &r < b[k] + n[k]) return pk[k][&r - b[k]];
        return nan;
    }
};

// ============================================================================
// Part 1: M4 per selection
// ============================================================================
struct M4Set {
    std::vector<WSpec> log, lin;          // components, mode A, w
    std::vector<WSpec> var;               // sum T1+T2 in A/w, A/w_db, B/w, B/w_db (log binning)
    std::vector<WSpec> varT1;             // signal T1 in the same variants
    DSpec dLog, dLin;
};
const char* const kVarName[4] = {"mode A, w (PDG)", "mode A, w_{db}", "mode B, w (PDG)", "mode B, w_{db}"};

void DrawM4Figure(int ip, const M4Set& m, double yLo, double yHi) {
    const std::string tag = kPTag[ip];
    const int p = kPres[ip];
    TH1D *hL[kNC], *hZ[kNC];
    for (int k = 0; k < kNC; ++k) {
        hL[k] = m.log[k].Hist(Form("m4L_%s_%d", tag.c_str(), k), 2);
        hZ[k] = m.lin[k].Hist(Form("m4Z_%s_%d", tag.c_str(), k), 2);
    }
    TCanvas* c = new TCanvas(Form("cM4_%s", tag.c_str()), "", 2000, 2400);
    // ---- (a) separate curves, log ----
    c->cd();
    MkPad(Form("a_%s", tag.c_str()), 0.0, 0.60, 1.0, 1.0, true, 0.07, 0.015, 0.17, 0.11);
    TH1D* fr = Frame(hL[cSum], Form("frA_%s", tag.c_str()), yLo, yHi);
    Axes(fr, kM4t, kM4y, 0.050, 0.62);
    fr->Draw("axis");
    for (int k : kBkgOrder) { StyleComp(hL[k], k); hL[k]->Draw("hist same"); }
    StyleComp(hL[cT2], cT2); hL[cT2]->Draw("hist same");
    for (int k : {cFake, cDiff}) { StyleComp(hL[k], k); hL[k]->Draw("hist same"); }
    for (int k : kSigOrder) { StyleComp(hL[k], k); hL[k]->Draw("hist same"); }
    StyleComp(hL[cSum], cSum); hL[cSum]->Draw("hist same");
    fr->Draw("axis same");
    TLegend* lg = new TLegend(0.075, 0.585, 0.985, 0.825);
    lg->SetNColumns(3); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.027);
    const int ordL[16] = {cSum, cT2, cFake, cBD, cP7, cDiff, cBC, cP8, cP3, cPC, cP1, cP4, cP0, cP6, cP25, -1};
    for (int k : ordL) {
        if (k < 0) continue;
        lg->AddEntry(hL[k], Form("%s: %.3g nb", kCD[k].leg, m.log[k].S / kLsim), "l");
    }
    lg->Draw();
    Txt(0.07, 0.94, Form("(a) simulated d#sigma/dM_{#pi^{+}#pi^{-}e^{+}e^{-}} per component, selection \"%s\": %s; 5 MeV bins; legend: #sigma in 0.30-1.10 GeV/c^{2}",
                         PresetName(p), PresetLabel(p)), 0.034);
    Txt(0.07, 0.895, kSimLine.c_str(), 0.030);
    Txt(0.07, 0.855, kDirect, 0.030, 2);
    // ---- (b) stacked, log ----
    c->cd();
    MkPad(Form("b_%s", tag.c_str()), 0.0, 0.30, 0.5, 0.60, true, 0.14, 0.02, 0.19, 0.14);
    std::vector<int> ord(kBkgOrder, kBkgOrder + kNBkg);
    for (int k : kSigOrder) ord.push_back(k);
    std::vector<TH1D*> stL = Cumul(hL, ord, "stL_" + tag);
    TH1D* frb = Frame(hL[cSum], Form("frB_%s", tag.c_str()), yLo, yHi / 60.);
    Axes(frb, kM4t, kM4y, 0.055, 1.2);
    frb->Draw("axis");
    DrawStack(stL);
    for (int k : {cFake, cDiff}) hL[k]->Draw("hist same");
    frb->Draw("axis same");
    Txt(0.14, 0.945, "(b) stacked: background by pion origin, signals on top (colours as in (a)), log", 0.042);
    Txt(0.14, 0.905, "dashed: fake lepton, different sources (not stacked, not in the sum)", 0.036);
    // ---- (c) stacked, linear zoom ----
    c->cd();
    MkPad(Form("c_%s", tag.c_str()), 0.5, 0.30, 1.0, 0.60, false, 0.14, 0.02, 0.19, 0.14);
    std::vector<TH1D*> stZ = Cumul(hZ, ord, "stZ_" + tag);
    const double ymZ = MaxOf(stZ.back()) * 1.12;
    TH1D* frc = Frame(hZ[cSum], Form("frC_%s", tag.c_str()), 0, ymZ);
    Axes(frc, kM4t, kM4y, 0.055, 1.2);
    frc->Draw("axis");
    DrawStack(stZ);
    TH1D* hsZ = (TH1D*)hZ[cSum]->Clone(Form("hsZ_%s", tag.c_str()));
    StyleLine(hsZ, gLine[cSum], 2, 1); hsZ->Draw("E1 same");
    for (double xl : {scan::kBumpLo, scan::kBumpHi, scan::kPeakLo, scan::kPeakHi}) VLine(xl, 0, ymZ, kGray + 2, 3, 1);
    frc->Draw("axis same");
    Txt(0.14, 0.945, "(c) linear zoom, 2.5 MeV bins, stacked (colours as in (a))", 0.042);
    Txt(0.14, 0.905, "error bars: sum, event-group errors; dotted: data windows (bump, peak)", 0.036);
    // ---- (d) separate curves, linear zoom ----
    c->cd();
    MkPad(Form("d_%s", tag.c_str()), 0.0, 0.0, 0.5, 0.30, false, 0.14, 0.02, 0.19, 0.14);
    double ymd = 0;
    for (int k = 0; k < kNC; ++k) if (k != cSum && k != cT1) ymd = std::max(ymd, MaxOf(hZ[k]));
    TH1D* frd = Frame(hZ[cSum], Form("frD_%s", tag.c_str()), 0, ymd * 1.12);
    Axes(frd, kM4t, kM4y, 0.055, 1.2);
    frd->Draw("axis");
    for (int k : kBkgOrder) { StyleComp(hZ[k], k); hZ[k]->Draw("hist same"); }
    StyleComp(hZ[cT2], cT2); hZ[cT2]->Draw("hist same");
    for (int k : {cFake, cDiff}) { StyleComp(hZ[k], k); hZ[k]->Draw("hist same"); }
    for (int k : kSigOrder) {
        StyleComp(hZ[k], k); hZ[k]->Draw("hist same");
        TH1D* he = (TH1D*)hZ[k]->Clone(Form("heZ_%s_%d", tag.c_str(), k)); StyleLine(he, gLine[k], 1, 1); he->Draw("E1 same");
    }
    frd->Draw("axis same");
    Txt(0.14, 0.945, "(d) linear zoom, separate curves (colours as in (a)), 2.5 MeV bins", 0.042);
    Txt(0.14, 0.905, Form("signals: bump #pi^{0} Dal. %.3g, bump conv. %.3g, peak conv. %.3g nb (0.40-0.70)",
                          m.lin[cBD].S / kLsim, m.lin[cBC].S / kLsim, m.lin[cPC].S / kLsim), 0.036);
    // ---- (e) variants of the sum ----
    c->cd();
    MkPad(Form("e1_%s", tag.c_str()), 0.5, 0.11, 1.0, 0.30, true, 0.14, 0.02, 0.17, 0.0);
    TH1D* hv[4];
    const int vc[4] = {gLine[cSum], kAzure + 2, kOrange + 7, kMagenta + 2};
    const int vs[4] = {1, 2, 1, 2};
    for (int v = 0; v < 4; ++v) {
        hv[v] = m.var[v].Hist(Form("var_%s_%d", tag.c_str(), v), 2);
        StyleLine(hv[v], vc[v], v == 0 ? 4 : 3, vs[v]);
    }
    double ymv = 0;
    for (int v = 0; v < 4; ++v) ymv = std::max(ymv, MaxOf(hv[v]));
    TH1D* fre = Frame(hv[0], Form("frE_%s", tag.c_str()), ymv * 2e-3, ymv * 30);
    Axes(fre, "", kM4y, 0.075, 0.85);
    fre->GetXaxis()->SetLabelSize(0);
    fre->Draw("axis");
    for (int v = 3; v >= 0; --v) hv[v]->Draw("hist same");
    TLegend* lv = new TLegend(0.16, 0.50, 0.98, 0.80);
    lv->SetNColumns(2); lv->SetBorderSize(0); lv->SetFillStyle(0); lv->SetTextSize(0.060);
    for (int v = 0; v < 4; ++v) lv->AddEntry(hv[v], Form("%s: %.4g nb", kVarName[v], m.var[v].S / kLsim), "l");
    lv->Draw();
    Txt(0.14, 0.905, "(e) sum (signal + phys. bkg): partition mode A / B, weight w (PDG) / w_{db}", 0.062);
    c->cd();
    MkPad(Form("e2_%s", tag.c_str()), 0.5, 0.0, 1.0, 0.11, false, 0.14, 0.02, 0.0, 0.33);
    TH1D* fre2 = Frame(hv[0], Form("frE2_%s", tag.c_str()), 0.85, 1.65);
    Axes(fre2, kM4t, "/ mode A, w", 0.13, 0.45);
    fre2->GetYaxis()->SetNdivisions(504);
    fre2->Draw("axis");
    HLine(0.3, 1.1, 1.0, kGray + 2, 2, 1);
    for (int v = 1; v < 4; ++v) {
        TH1D* r = Ratio(hv[v], hv[0], Form("rv_%s_%d", tag.c_str(), v));
        StyleLine(r, vc[v], 3, vs[v]);
        r->Draw("hist same");
    }
    c->SaveAs(FigPath("spectra_m4_" + tag + ".png").c_str());
    c->SaveAs(FigPath("spectra_m4_" + tag + ".pdf").c_str());
}

// data against simulation, M4
struct RegRes { double dS = 0, dE = 0, sS = 0, sE = 0, r = 0, rE = 0, N = 0, CB = 0; };
void DrawM4DataFigure(int ip, const M4Set& m, double yLo, double yHi, const std::vector<RegRes>& reg) {
    const std::string tag = kPTag[ip];
    const int p = kPres[ip];
    TH1D *hL[kNC], *hZ[kNC];
    for (int k = 0; k < kNC; ++k) {
        hL[k] = m.log[k].Hist(Form("dm4L_%s_%d", tag.c_str(), k), 2);
        hZ[k] = m.lin[k].Hist(Form("dm4Z_%s_%d", tag.c_str(), k), 2);
    }
    TH1D* hN = m.dLog.Hist("us", Form("dN_%s", tag.c_str()), 2);
    TH1D* hC = m.dLog.Hist("cb", Form("dC_%s", tag.c_str()), 2);
    TH1D* hS = m.dLog.Hist("sig", Form("dS_%s", tag.c_str()), 2);
    TH1D* hSZ = m.dLin.Hist("sig", Form("dSZ_%s", tag.c_str()), 2);
    StyleData(hN, 1, 20, 1.0); StyleData(hC, 2, 21, 1.0); StyleData(hS, 4, 20, 1.0); StyleData(hSZ, 4, 20, 1.1);
    const double f = reg[kRegNorm].r;   // data / sim in 0.58-0.70
    TCanvas* c = new TCanvas(Form("cD_%s", tag.c_str()), "", 2000, 2400);
    // ---- (a) log ----
    c->cd();
    MkPad(Form("da_%s", tag.c_str()), 0.0, 0.60, 1.0, 1.0, true, 0.07, 0.015, 0.17, 0.11);
    TH1D* fr = Frame(hN, Form("frDa_%s", tag.c_str()), yLo, yHi);
    Axes(fr, kM4t, kM4y, 0.050, 0.62);
    fr->Draw("axis");
    TH1D* b2 = (TH1D*)hL[cT2]->Clone(Form("b2_%s", tag.c_str())); StyleFill(b2, gLine[cP7], gFill[cP7]); b2->Draw("hist same");
    for (int k : {cFake, cDiff}) { StyleComp(hL[k], k); hL[k]->Draw("hist same"); }
    for (int k : kSigOrder) { StyleComp(hL[k], k); hL[k]->Draw("hist same"); }
    StyleComp(hL[cSum], cSum); hL[cSum]->Draw("hist same");
    hN->Draw("E1 P same"); hC->Draw("E1 P same"); hS->Draw("E1 P same");
    fr->Draw("axis same");
    TLegend* lg = new TLegend(0.075, 0.62, 0.985, 0.825);
    lg->SetNColumns(3); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.028);
    lg->AddEntry(hN, "data: all e^{+}e^{-} combinations N", "p");
    lg->AddEntry(hL[cSum], "sim: signal + phys. bkg (sum)", "l");
    lg->AddEntry(hL[cBD], "sim: bump, #pi^{0} Dalitz", "l");
    lg->AddEntry(hC, "data: CB = 2#sqrt{N_{++}N_{--}} (no k)", "p");
    lg->AddEntry(b2, "sim: physical background", "f");
    lg->AddEntry(hL[cBC], "sim: bump, conversion", "l");
    lg->AddEntry(hS, "data: N - CB", "p");
    lg->AddEntry(hL[cFake], "sim: fake lepton, diff. sources (dashed, not in sum)", "l");
    lg->AddEntry(hL[cPC], "sim: peak, conversion", "l");
    lg->Draw();
    Txt(0.07, 0.94, Form("(a) d#sigma/dM_{#pi^{+}#pi^{-}e^{+}e^{-}}, data and sim, same units, NOT renormalised, 5 MeV; \"%s\"", PresetName(p)), 0.031);
    Txt(0.985, 0.835, PresetLabel(p), 0.026, 1, 31);
    Txt(0.07, 0.895, kDataLine.c_str(), 0.030, 4);
    Txt(0.07, 0.855, kSimLine.c_str(), 0.026);
    Txt(0.60, 0.94, kDirect, 0.024, 2);
    // ---- (b) linear zoom, not renormalised ----
    std::vector<int> ord(kBkgOrder, kBkgOrder + kNBkg);
    for (int k : kSigOrder) ord.push_back(k);
    for (int pass = 0; pass < 2; ++pass) {
        const double sc = pass == 0 ? 1.0 : f;
        c->cd();
        MkPad(Form("db%d_%s", pass, tag.c_str()), pass == 0 ? 0.0 : 0.5, 0.30, pass == 0 ? 0.5 : 1.0, 0.60, false, 0.14, 0.02, 0.19, 0.14);
        std::vector<TH1D*> st = Cumul(hZ, ord, Form("dst%d_%s", pass, tag.c_str()), sc);
        const double ym = std::max(MaxOf(st.back()), MaxOf(hSZ, true)) * 1.12;
        TH1D* frb = Frame(hZ[cSum], Form("frDb%d_%s", pass, tag.c_str()), 0, ym);
        Axes(frb, kM4t, kM4y, 0.055, 1.2);
        frb->Draw("axis");
        DrawStack(st);
        hSZ->Draw("E1 P same");
        for (double xl : {scan::kBumpLo, scan::kBumpHi, scan::kPeakLo, scan::kPeakHi}) VLine(xl, 0, ym, kGray + 2, 3, 1);
        frb->Draw("axis same");
        if (pass == 0) {
            Txt(0.14, 0.945, "(b) linear zoom, 2.5 MeV: data N - CB (blue) and the simulation", 0.042);
            Txt(0.14, 0.905, "stacked (colours as in spectra_m4), NOT renormalised; dotted: windows", 0.036);
        } else {
            Txt(0.14, 0.945, Form("(c) the same, simulation #times %.3f (data/sim in 0.58-0.70)", f), 0.042);
            Txt(0.14, 0.905, "renormalisation only for the shape comparison; factor scales with 1/L_{data}", 0.036);
        }
        Txt(0.14, 0.865, "#eta #rightarrow #pi^{+}#pi^{-}e^{+}e^{-} (direct) NOT in the simulation", 0.036, 2);
    }
    // ---- (d) ratio, full range ----
    c->cd();
    MkPad(Form("dd_%s", tag.c_str()), 0.0, 0.0, 0.5, 0.30, false, 0.14, 0.02, 0.19, 0.14);
    TH1D* hr = Ratio(hS, hL[cSum], Form("dr_%s", tag.c_str()));
    StyleData(hr, 4, 20, 0.9);
    TH1D* frr = Frame(hr, Form("frDd_%s", tag.c_str()), 0.0, 2.5);
    Axes(frr, kM4t, "data (N - CB) / simulation", 0.055, 1.1);
    frr->Draw("axis");
    HLine(0.3, 1.1, 1.0, kGray + 2, 2, 1);
    for (double xl : {scan::kBumpLo, scan::kBumpHi, scan::kPeakLo, scan::kPeakHi}) VLine(xl, 0, 2.5, kGray + 1, 3, 1);
    hr->Draw("E1 P same");
    for (int r = 0; r < kNReg; ++r) {
        if (reg[r].sS <= 0) continue;
        HLine(kRegLo[r], kRegHi[r], reg[r].r, kOrange + 7, 1, 4);
    }
    Txt(0.14, 0.945, "(d) data (N - CB) / simulation (signal + phys. bkg), 5 MeV", 0.042);
    Txt(0.14, 0.905, "orange bars: ratio integrated in the regions (table spectra_regions)", 0.036);
    Txt(0.14, 0.865, "ratio scales with 1/L_{data} (unconfirmed); no direct #eta in sim", 0.036, 2);
    // ---- (e) residual: data - renormalised sim background vs renormalised sim signals ----
    c->cd();
    MkPad(Form("de_%s", tag.c_str()), 0.5, 0.0, 1.0, 0.30, false, 0.14, 0.02, 0.19, 0.14);
    TH1D* res = (TH1D*)hSZ->Clone(Form("res_%s", tag.c_str()));
    for (int i = 1; i <= res->GetNbinsX(); ++i) {
        res->SetBinContent(i, hSZ->GetBinContent(i) - f * hZ[cT2]->GetBinContent(i));
        res->SetBinError(i, std::sqrt(std::pow(hSZ->GetBinError(i), 2) + std::pow(f * hZ[cT2]->GetBinError(i), 2)));
    }
    std::vector<TH1D*> sig = Cumul(hZ, {cBD, cBC, cPC}, Form("dsig_%s", tag.c_str()), f);
    double ym = std::max(MaxOf(res, true), MaxOf(sig.back())) * 1.2;
    double yl = 0;
    for (int i = 1; i <= res->GetNbinsX(); ++i) yl = std::min(yl, res->GetBinContent(i) - res->GetBinError(i));
    TH1D* fre = Frame(hZ[cSum], Form("frDe_%s", tag.c_str()), yl * 1.1, ym);
    Axes(fre, kM4t, kM4y, 0.055, 1.2);
    fre->Draw("axis");
    DrawStack(sig);
    HLine(0.40, 0.70, 0, kGray + 2, 2, 1);
    res->Draw("E1 P same");
    for (double xl : {scan::kBumpLo, scan::kBumpHi, scan::kPeakLo, scan::kPeakHi}) VLine(xl, yl * 1.1, ym, kGray + 2, 3, 1);
    fre->Draw("axis same");
    Txt(0.14, 0.945, Form("(e) data (N - CB) - %.3f #times sim phys. bkg (blue)", f), 0.042);
    Txt(0.14, 0.905, Form("against %.3f #times sim signals (stacked); bkg shape not verified", f), 0.036);
    Txt(0.14, 0.865, "#eta #rightarrow #pi^{+}#pi^{-}e^{+}e^{-} (direct) NOT in the simulation", 0.036, 2);
    c->SaveAs(FigPath("spectra_m4data_" + tag + ".png").c_str());
    c->SaveAs(FigPath("spectra_m4data_" + tag + ".pdf").c_str());
}

// ============================================================================
// Part 2: subsystems
// ============================================================================
struct SubVar {
    std::string tag;
    std::string xt, ysym, unit;
    int nbin; double lo, hi;
    std::function<double(const Row&)> fs;
    std::function<double(const scan::Row&)> fd;
    int win;      // 0 all M4, 1 bump window, 2 peak window
    bool logy;
};
std::vector<SubVar> SubVars(const DataPK& pk, const std::string& group) {
    std::vector<SubVar> v;
    const DataPK* P = &pk;
    auto F = [](float Row::*m) { return [m](const Row& r) { return (double)(r.*m); }; };
    auto D = [](float scan::Row::*m) { return [m](const scan::Row& r) { return (double)(r.*m); }; };
    auto K = [P](float PairKin::*m) { return [P, m](const scan::Row& r) { return (double)(P->Of(r).*m); }; };
    const char* GeV = "GeV/c^{2}";
    if (group == "pipi") {
        v.push_back({"mpipi", "M_{#pi^{+}#pi^{-}} [GeV/c^{2}]", "M_{#pi#pi}", GeV, 100, 0.27, 1.27, F(&Row::mpipi), D(&scan::Row::mpipi), 0, true});
        v.push_back({"p_pipi", "p_{#pi^{+}#pi^{-}} [MeV/c]", "p_{#pi#pi}", "MeV/c", 70, 0, 3500, F(&Row::p_pipi), K(&PairKin::p_pipi), 0, true});
        v.push_back({"th_pipi", "#theta_{#pi^{+}#pi^{-}} [deg]", "#theta_{#pi#pi}", "deg", 45, 0, 90, F(&Row::th_pipi), K(&PairKin::th_pipi), 0, true});
        v.push_back({"pt_pipi", "p_{T,#pi^{+}#pi^{-}} [MeV/c]", "p_{T,#pi#pi}", "MeV/c", 60, 0, 1200, F(&Row::pt_pipi), K(&PairKin::pt_pipi), 0, true});
        v.push_back({"mm_pipi", "MM(#pi^{+}#pi^{-}) [GeV/c^{2}]", "MM_{#pi#pi}", GeV, 85, 1.5, 3.2, F(&Row::mm_pipi), D(&scan::Row::mm_pipi), 0, true});
    } else if (group == "ee") {
        v.push_back({"mee", "M_{e^{+}e^{-}} [GeV/c^{2}]", "M_{ee}", GeV, 150, 0, 0.6, F(&Row::mee), D(&scan::Row::mee), 0, true});
        v.push_back({"oa", "OA(e^{+}e^{-}) [deg]", "OA", "deg", 90, 0, 90, F(&Row::oa), D(&scan::Row::oa), 0, true});
        v.push_back({"p_ee", "p_{e^{+}e^{-}} [MeV/c]", "p_{ee}", "MeV/c", 50, 0, 2500, F(&Row::p_ee), K(&PairKin::p_ee), 0, true});
        v.push_back({"th_ee", "#theta_{e^{+}e^{-}} [deg]", "#theta_{ee}", "deg", 45, 0, 90, F(&Row::th_ee), K(&PairKin::th_ee), 0, true});
        v.push_back({"pt_ee", "p_{T,e^{+}e^{-}} [MeV/c]", "p_{T,ee}", "MeV/c", 50, 0, 1000, F(&Row::pt_ee), K(&PairKin::pt_ee), 0, true});
        v.push_back({"mm_ee", "MM(e^{+}e^{-}) [GeV/c^{2}]", "MM_{ee}", GeV, 60, 2.2, 3.4, F(&Row::mm_ee), D(&scan::Row::mm_ee), 0, true});
    } else if (group == "tracks") {
        v.push_back({"pip_p", "p_{#pi^{+}} [MeV/c]", "p_{#pi^{+}}", "MeV/c", 52, 0, 2600, F(&Row::pip_p), D(&scan::Row::pip_p), 0, true});
        v.push_back({"pim_p", "p_{#pi^{-}} [MeV/c]", "p_{#pi^{-}}", "MeV/c", 52, 0, 2600, F(&Row::pim_p), D(&scan::Row::pim_p), 0, true});
        v.push_back({"ep_p", "p_{e^{+}} [MeV/c]", "p_{e^{+}}", "MeV/c", 40, 0, 2000, F(&Row::ep_p), D(&scan::Row::ep_p), 0, true});
        v.push_back({"em_p", "p_{e^{-}} [MeV/c]", "p_{e^{-}}", "MeV/c", 40, 0, 2000, F(&Row::em_p), D(&scan::Row::em_p), 0, true});
        v.push_back({"pip_theta", "#theta_{#pi^{+}} [deg]", "#theta_{#pi^{+}}", "deg", 45, 0, 90, F(&Row::pip_theta), D(&scan::Row::pip_theta), 0, true});
        v.push_back({"pim_theta", "#theta_{#pi^{-}} [deg]", "#theta_{#pi^{-}}", "deg", 45, 0, 90, F(&Row::pim_theta), D(&scan::Row::pim_theta), 0, true});
        v.push_back({"ep_theta", "#theta_{e^{+}} [deg]", "#theta_{e^{+}}", "deg", 45, 0, 90, F(&Row::ep_theta), D(&scan::Row::ep_theta), 0, true});
        v.push_back({"em_theta", "#theta_{e^{-}} [deg]", "#theta_{e^{-}}", "deg", 45, 0, 90, F(&Row::em_theta), D(&scan::Row::em_theta), 0, true});
    } else {
        v.push_back({"oa_lab", "#alpha_{lab}(#pi^{+}#pi^{-}, e^{+}e^{-}) [deg]", "#alpha_{lab}", "deg", 90, 0, 180, F(&Row::oa_lab), D(&scan::Row::oa_lab), 0, true});
        v.push_back({"mm4", "MM(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]", "MM_{4}", GeV, 105, 1.0, 3.1, F(&Row::mm4), D(&scan::Row::mm4), 0, true});
        v.push_back({"pt4", "p_{T,#pi^{+}#pi^{-}e^{+}e^{-}} [MeV/c]", "p_{T,4}", "MeV/c", 60, 0, 1200, F(&Row::pt4), D(&scan::Row::pt4), 0, true});
        v.push_back({"mpipi_bump", "M_{#pi^{+}#pi^{-}} [GeV/c^{2}], M_{4} in the bump window", "M_{#pi#pi}", GeV, 56, 0.27, 0.55, F(&Row::mpipi), D(&scan::Row::mpipi), 1, false});
        v.push_back({"mpipi_peak", "M_{#pi^{+}#pi^{-}} [GeV/c^{2}], M_{4} in the peak window", "M_{#pi#pi}", GeV, 56, 0.27, 0.55, F(&Row::mpipi), D(&scan::Row::mpipi), 2, false});
        v.push_back({"mee_bump", "M_{e^{+}e^{-}} [GeV/c^{2}], M_{4} in the bump window", "M_{ee}", GeV, 100, 0, 0.25, F(&Row::mee), D(&scan::Row::mee), 1, true});
        v.push_back({"mee_peak", "M_{e^{+}e^{-}} [GeV/c^{2}], M_{4} in the peak window", "M_{ee}", GeV, 100, 0, 0.25, F(&Row::mee), D(&scan::Row::mee), 2, true});
    }
    return v;
}

void DrawSubFigure(const SimData& d, const scan::Data& dd, const DataPK& pk, int ip, const std::string& group,
                   int ncol, int nrow, std::map<std::string, std::pair<double, double>>& yr, FILE* tab) {
    const std::string tag = kPTag[ip];
    const int p = kPres[ip];
    const std::vector<SubVar> vars = SubVars(pk, group);
    const std::vector<Comp> comps = Comps();
    SimOpt opt;
    const double W = 640.0 * ncol, H = 600.0 * nrow / 0.92;
    TCanvas* c = new TCanvas(Form("cS_%s_%s", tag.c_str(), group.c_str()), "", (int)W, (int)H);
    const char* gname = group == "pipi" ? "#pi^{+}#pi^{-} pair" : group == "ee" ? "e^{+}e^{-} pair"
                      : group == "tracks" ? "single tracks" : "four-body quantities and the windows";
    c->cd();
    const double fsT = W / 95.0 / H, fsL = W / 140.0 / H;
    Txt(0.01, 0.975, Form("Simulated d#sigma/dx per component, %s; selection \"%s\": %s; panels: top = cross section, bottom = data (N - CB) / sim (sum)",
                          gname, PresetName(p), PresetLabel(p)), fsT);
    Txt(0.01, 0.957, kSimLine.c_str(), fsL, 1);
    Txt(0.01, 0.942, kDataLine.c_str(), fsL, 4);
    Txt(0.01, 0.927, (std::string(kDirect) + "; log panels: same y range in every selection (that of \"pre\"); linear panels scaled per selection").c_str(), fsL, 2);
    const double top = 0.92;
    const double cw = 1.0 / ncol, ch = top / nrow;
    TH1D* legH[kNC] = {nullptr};
    TH1D* legD = nullptr;
    for (size_t iv = 0; iv < vars.size(); ++iv) {
        const SubVar& v = vars[iv];
        const int col = (int)iv % ncol, row = (int)iv / ncol;
        const double x1 = col * cw, x2 = x1 + cw, y2 = top - row * ch, y1 = y2 - ch;
        const double ylim = y1 + 0.30 * ch;
        // simulation
        const int win = v.win;
        Sel sel = SelPreset(p);
        if (win) {
            const double a = win == 1 ? scan::kBumpLo : scan::kPeakLo, b = win == 1 ? scan::kBumpHi : scan::kPeakHi;
            sel = SelAnd(sel, [a, b](const Row& r) { return r.m4 >= a && r.m4 < b; });
        }
        Binning bn; bn.nbin = v.nbin; bn.lo = v.lo; bn.hi = v.hi; bn.x = v.fs;
        const std::vector<WSpec> sp = FillMany(d, sel, comps, opt, bn);
        // data
        scan::Sel dsel = DataSelPreset(p);
        if (win) {
            const float a = win == 1 ? scan::kBumpLo : scan::kPeakLo, b = win == 1 ? scan::kBumpHi : scan::kPeakHi;
            dsel = [p, a, b](const scan::Row& r) { return PassPreset(p, r) && r.m4 >= a && r.m4 < b; };
        }
        const DSpec ds = DataSpectrum(dd, dsel, v.nbin, v.lo, v.hi, v.fd, kLdata, 1);
        TH1D* h[kNC];
        for (int k = 0; k < kNC; ++k) h[k] = sp[k].Hist(Form("sub_%s_%s_%d", tag.c_str(), v.tag.c_str(), k), 2);
        TH1D* hd = ds.Hist("sig", Form("subd_%s_%s", tag.c_str(), v.tag.c_str()), 2);
        StyleData(hd, 4, 20, 0.8);
        // table: fractions outside the range, integrals
        if (tab) {
            double dIn = 0, dAll = ds.sp.nsel[0];
            for (int i = 0; i < ds.nbin; ++i) dIn += ds.sp.n[0][i];
            fprintf(tab, "%-9s %-11s sim sum: in range %8.4g nb (all %8.4g nb, N_eff %8.0f); signal T1 %7.4g nb; data N - CB in range %8.4g nb; data rows outside %.3f%%\n",
                    tag.c_str(), v.tag.c_str(), sp[cSum].S / kLsim, sp[cSum].Sall / kLsim, sp[cSum].NeffIn(),
                    sp[cT1].S / kLsim, (std::accumulate(ds.sp.sig.begin(), ds.sp.sig.end(), 0.0)) / kLdata,
                    dAll > 0 ? 100.0 * (1 - dIn / dAll) : 0.);
        }
        // y range
        double ymax = std::max(MaxOf(h[cSum]), MaxOf(hd)), ymin;
        if (v.logy) {
            auto it = yr.find(v.tag);
            if (it == yr.end()) {
                double mp = 1e300;
                for (int k : kSigOrder) mp = std::min(mp, MinPos(h[k]));
                ymin = std::max(ymax * 3e-6, mp * 0.3);
                if (mp > 1e299) ymin = ymax * 1e-4;
                ymax *= 8;
                yr[v.tag] = std::make_pair(ymin, ymax);
            } else { ymin = it->second.first; ymax = it->second.second; }
        } else { ymin = 0; ymax *= 1.25; }
        // main pad
        c->cd();
        TPad* pm = MkPad(Form("pm_%s_%s", tag.c_str(), v.tag.c_str()), x1, ylim, x2, y2, v.logy, 0.17, 0.03, 0.06, 0.0);
        TH1D* fr = Frame(h[cSum], Form("frS_%s_%s", tag.c_str(), v.tag.c_str()), ymin, ymax);
        Axes(fr, "", Form("d#sigma/d%s [nb/(%s)]", v.ysym.c_str(), v.unit.c_str()), 0.060, 1.35);
        fr->GetXaxis()->SetLabelSize(0);
        fr->Draw("axis");
        std::vector<int> ord(kBkgOrder, kBkgOrder + kNBkg);
        if (!v.logy) for (int k : kSigOrder) ord.push_back(k);
        std::vector<TH1D*> st = Cumul(h, ord, Form("sst_%s_%s", tag.c_str(), v.tag.c_str()));
        DrawStack(st);
        for (int k : {cFake, cDiff}) { StyleComp(h[k], k, 0.7); h[k]->Draw("hist same"); }
        if (v.logy) for (int k : kSigOrder) { StyleComp(h[k], k, 0.8); h[k]->Draw("hist same"); }
        StyleComp(h[cSum], cSum, 0.6); h[cSum]->Draw("hist same");
        hd->Draw("E1 P same");
        fr->Draw("axis same");
        if (win) Txt(0.20, 0.86, win == 1 ? "M_{4} in 0.455-0.5225 (bump)" : "M_{4} in 0.5336-0.5542 (peak)", 0.065, 1);
        if (!legH[0]) {
            for (int k = 0; k < kNC; ++k) legH[k] = h[k];
            for (size_t j = 0; j < st.size(); ++j) legH[ord[j]] = st[j];
            legD = hd;
        }
        pm->RedrawAxis();
        // ratio pad
        c->cd();
        MkPad(Form("pr_%s_%s", tag.c_str(), v.tag.c_str()), x1, y1, x2, ylim, false, 0.17, 0.03, 0.0, 0.36);
        TH1D* r = Ratio(hd, h[cSum], Form("sr_%s_%s", tag.c_str(), v.tag.c_str()));
        StyleData(r, 4, 20, 0.6);
        TH1D* frr = Frame(r, Form("frSr_%s_%s", tag.c_str(), v.tag.c_str()), 0.0, 2.5);
        Axes(frr, v.xt.c_str(), "data/sim", 0.13, 0.62);
        frr->GetYaxis()->SetNdivisions(503);
        frr->GetXaxis()->SetTitleOffset(1.15);
        frr->Draw("axis");
        HLine(v.lo, v.hi, 1.0, kGray + 2, 2, 1);
        r->Draw("E1 P same");
    }
    // legend cell
    {
        const int iv = (int)vars.size();
        const int col = iv % ncol, row = iv / ncol;
        const double x1 = col * cw, x2 = x1 + cw, y2 = top - row * ch, y1 = y2 - ch;
        c->cd();
        TPad* pl = MkPad(Form("pl_%s_%s", tag.c_str(), group.c_str()), x1, y1, x2, y2, false, 0.02, 0.02, 0.02, 0.02);
        TLegend* lg = new TLegend(0.02, 0.02, 0.98, 0.98);
        lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.046);
        lg->AddEntry(legD, "data: N - CB (blue)", "p");
        lg->AddEntry(legH[cSum], "sim: signal + phys. bkg (sum)", "l");
        for (int k : kSigOrder) lg->AddEntry(legH[k], kCD[k].leg, "l");
        for (int j = kNBkg - 1; j >= 0; --j) lg->AddEntry(legH[kBkgOrder[j]], kCD[kBkgOrder[j]].leg, "f");
        lg->AddEntry(legH[cFake], "fake lepton (not in the sum)", "l");
        lg->AddEntry(legH[cDiff], "different sources (not in the sum)", "l");
        lg->AddEntry((TObject*)nullptr, "signals: lines (log), stacked (linear)", "");
        lg->Draw();
        (void)pl;
    }
    c->SaveAs(FigPath("spectra_sub_" + tag + "_" + group + ".png").c_str());
    if (ip == 0) c->SaveAs(FigPath("spectra_sub_" + tag + "_" + group + ".pdf").c_str());
}

}  // namespace ssp

// ============================================================================
// main
// ============================================================================
void sim_spectra(int parts = 15) {
    using namespace ssp;
    InitColors();
    TStopwatch swAll;
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetEndErrorSize(0);
    gErrorIgnoreLevel = kError;

    SimData d;
    // all rows of std_us and lep_us (owned = 0: the isBest check needs every row; spectra use Owned())
    if (!Load(d, kCompact, (1 << kStdUS) | (1 << kLepUS), -1e9, 1e9, 0)) return;
    // event groups: exactly one isBest == 1 row per group (adjacent rows, grp column)
    for (int set = 1; set <= 2; ++set) {
        const std::vector<Row>& v = d.tree(set, kUS);
        long ng = 0, bad = 0; int cur = -1, nb = 0;
        for (size_t i = 0; i <= v.size(); ++i) {
            if (i == v.size() || v[i].grp != cur) {
                if (i > 0) { ++ng; if (nb != 1) ++bad; }
                if (i == v.size()) break;
                cur = v[i].grp; nb = 0;
            }
            if (v[i].isBest == 1) ++nb;
        }
        printf("sim_spectra: %s: %zu rows, %ld event groups, groups without exactly one isBest == 1 row: %ld\n",
               TreeName(TreeIndex(set, kUS)), v.size(), ng, bad);
    }
    scan::Data dd;
    if (!scan::Load(dd, scan::kCompact)) return;
    const std::vector<Comp> comps = Comps();
    SimOpt optA;
    Binning bLog; bLog.nbin = 160; bLog.lo = 0.30; bLog.hi = 1.10;
    Binning bLin; bLin.nbin = 120; bLin.lo = 0.40; bLin.hi = 0.70;
    TFile* fo = (parts & 1) ? TFile::Open((std::string(kScanDir) + "/sim_spectra.root").c_str(), "RECREATE") : nullptr;

    // ======================= Part 1 and 8: M4 =======================
    std::vector<M4Set> M(kNP);
    std::vector<std::vector<RegRes>> REG(kNP, std::vector<RegRes>(kNReg + 1));
    if (parts & (1 | 8)) {
        for (int ip = 0; ip < kNP; ++ip) {
            const int p = kPres[ip];
            const Sel sel = SelPreset(p);
            M4Set& m = M[ip];
            m.log = FillMany(d, sel, comps, optA, bLog);
            m.lin = FillMany(d, sel, comps, optA, bLin);
            for (int v = 0; v < 4; ++v) {
                SimOpt o; o.mode = v / 2; o.db = v % 2;
                const std::vector<WSpec> s = FillMany(d, sel, {comps[cSum], comps[cT1]}, o, bLog);
                m.var.push_back(s[0]); m.varT1.push_back(s[1]);
            }
            m.dLog = DataSpectrumPreset(dd, p, bLog.nbin, bLog.lo, bLog.hi);
            m.dLin = DataSpectrumPreset(dd, p, bLin.nbin, bLin.lo, bLin.hi);
            // regions: simulation (one pass, exact group variance), data (ScanEngine windows)
            Binning br; br.nbin = kNReg; br.lo = 0; br.hi = kNReg;
            br.x = [](const Row& r) {
                for (int k = 0; k < kNReg; ++k) if (r.m4 >= kRegLo[k] && r.m4 < kRegHi[k]) return k + 0.5;
                return -1.0;
            };
            const WSpec sr = Fill(d, sel, comps[cSum], optA, br);
            Binning bt; bt.nbin = 1; bt.lo = 0.30; bt.hi = 1.10;
            const WSpec stt = Fill(d, sel, comps[cSum], optA, bt);
            scan::FillOpt fo2; fo2.nbin = 320; fo2.lo = 0.30; fo2.hi = 1.10; fo2.reb = 4;
            for (int r = 0; r < kNReg; ++r) fo2.windows.push_back(std::make_pair(kRegLo[r], kRegHi[r]));
            fo2.windows.push_back(std::make_pair(0.30, 1.10));
            const scan::Spectra dsp = scan::Fill(dd, DataSelPreset(p), fo2);
            for (int r = 0; r <= kNReg; ++r) {
                RegRes& q = REG[ip][r];
                const scan::WinSum& w = dsp.win[2 + r];
                q.N = w.n[0]; q.CB = w.cb;
                q.dS = w.sig() / kLdata; q.dE = w.sigErr() / kLdata;
                q.sS = (r < kNReg ? sr.s[r] : stt.S) / kLsim;
                q.sE = std::sqrt(r < kNReg ? sr.v[r] : stt.V) / kLsim;
                q.r = q.sS > 0 ? q.dS / q.sS : 0;
                q.rE = q.sS > 0 && q.dS != 0 ? q.r * std::sqrt(std::pow(q.dE / q.dS, 2) + std::pow(q.sE / q.sS, 2)) : 0;
            }
            printf("sim_spectra: M4 fills for \"%s\" done (%.1f s)\n", PresetName(p), swAll.RealTime());
            swAll.Continue();
        }
    }
    if (parts & 1) {
        // common log range (from "pre", the loosest selection; every other is a subset)
        double ymax = std::max(MaxOf(M[0].log[cSum].Hist("tmp1", 2)), MaxOf(M[0].dLog.Hist("us", "tmp2", 2)));
        double ymin = 1e300;
        for (int ip = 0; ip < kNP; ++ip)
            for (int k : kSigOrder) {
                TH1D* h = M[ip].log[k].Hist(Form("tmp3_%d_%d", ip, k), 2);
                ymin = std::min(ymin, MinPos(h));
            }
        const double yLo = std::max(ymin * 0.5, ymax * 1e-6), yHi = ymax * 600;
        for (int ip = 0; ip < kNP; ++ip) {
            DrawM4Figure(ip, M[ip], yLo, yHi);
            DrawM4DataFigure(ip, M[ip], yLo, yHi, REG[ip]);
            // store histograms
            fo->cd();
            TDirectory* dir = fo->mkdir(kPTag[ip], "", true);
            dir->cd();
            for (int k = 0; k < kNC; ++k) {
                M[ip].log[k].Hist(Form("m4_%s", kCD[k].tag), 2)->Write();
                M[ip].lin[k].Hist(Form("m4zoom_%s", kCD[k].tag), 2)->Write();
            }
            for (int v = 0; v < 4; ++v) M[ip].var[v].Hist(Form("m4_sum_variant%d", v), 2)->Write();
            for (const char* w : {"us", "cb", "sig"}) {
                M[ip].dLog.Hist(w, Form("data_m4_%s", w), 2)->Write();
                M[ip].dLin.Hist(w, Form("data_m4zoom_%s", w), 2)->Write();
            }
        }
        // ---- overview of the selections ----
        {
            TCanvas* c = new TCanvas("cPres", "", 2000, 1500);
            const int pc[kNP] = {kTeal + 3, kAzure + 2, kOrange + 7, kMagenta + 1, kGray + 2, kViolet + 1};
            const int pls[kNP] = {1, 2, 1, 2, 1, 2};
            c->cd();
            MkPad("pp1", 0.0, 0.45, 0.5, 1.0, true, 0.14, 0.02, 0.17, 0.12);
            TH1D* f1 = nullptr;
            TLegend* lg = new TLegend(0.55, 0.15, 0.97, 0.47);
            lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.034);
            for (int ip = 0; ip < kNP; ++ip) {
                TH1D* h = M[ip].log[cSum].Hist(Form("ps_%d", ip), 2);
                StyleLine(h, pc[ip], 3, pls[ip]);
                if (!f1) { f1 = Frame(h, "fp1", yLo * 10, MaxOf(h) * 5); Axes(f1, kM4t, kM4y, 0.050, 1.3); f1->Draw("axis"); }
                h->Draw("hist same");
                lg->AddEntry(h, Form("%s: %.4g nb", PresetName(kPres[ip]), M[ip].log[cSum].S / kLsim), "l");
            }
            lg->Draw();
            Txt(0.14, 0.93, "(a) simulated sum (signal + phys. bkg) per selection, 5 MeV", 0.042);
            Txt(0.14, 0.885, "std + lep, mode A, w (PDG), #sigma = #Sigma w / L_{sim}; legend: 0.30-1.10", 0.036);
            c->cd();
            MkPad("pp2", 0.5, 0.45, 1.0, 1.0, false, 0.14, 0.02, 0.17, 0.12);
            TH1D* f2 = nullptr;
            double ym = 0;
            std::vector<TH1D*> hs;
            for (int ip = 0; ip < kNP; ++ip) {
                TH1D* h = M[ip].lin[cT1].Hist(Form("pt1_%d", ip), 2);
                h->Rebin(2); h->Scale(0.5);
                StyleLine(h, pc[ip], 3, pls[ip]);
                ym = std::max(ym, MaxOf(h));
                hs.push_back(h);
            }
            f2 = Frame(hs[0], "fp2", 0, ym * 1.3); Axes(f2, kM4t, kM4y, 0.050, 1.3); f2->GetXaxis()->SetRangeUser(0.40, 0.60);
            f2->Draw("axis");
            TLegend* lg2 = new TLegend(0.16, 0.55, 0.60, 0.80);
            lg2->SetBorderSize(0); lg2->SetFillStyle(0); lg2->SetTextSize(0.034);
            for (int ip = 0; ip < kNP; ++ip) {
                hs[ip]->Draw("hist same");
                lg2->AddEntry(hs[ip], Form("%s: bump %.3g, peak conv. %.3g nb", kPTag[ip], (M[ip].lin[cBD].S + M[ip].lin[cBC].S) / kLsim,
                                           M[ip].lin[cPC].S / kLsim), "l");
            }
            lg2->Draw();
            Txt(0.14, 0.93, "(b) simulated signal (bump + peak conversion) per selection, 5 MeV", 0.042);
            Txt(0.14, 0.885, "#eta #rightarrow #pi^{+}#pi^{-}e^{+}e^{-} (direct) NOT in the simulation", 0.036, 2);
            c->cd();
            MkPad("pp3", 0.0, 0.0, 1.0, 0.45, false, 0.07, 0.015, 0.14, 0.14);
            TH1D* f3 = nullptr;
            TLegend* lg3 = new TLegend(0.20, 0.70, 0.98, 0.85);
            lg3->SetNColumns(3); lg3->SetBorderSize(0); lg3->SetFillStyle(0); lg3->SetTextSize(0.045);
            for (int ip = 0; ip < kNP; ++ip) {
                TH1D* hd = M[ip].dLog.Hist("sig", Form("pr_d_%d", ip), 2);
                TH1D* hsim = M[ip].log[cSum].Hist(Form("pr_s_%d", ip), 2);
                hd->Rebin(2); hsim->Rebin(2);
                TH1D* r = Ratio(hd, hsim, Form("prr_%d", ip));
                StyleData(r, pc[ip], 20 + ip, 1.0);
                if (!f3) { f3 = Frame(r, "fp3", 0.4, 2.3); Axes(f3, kM4t, "data (N - CB) / sim", 0.055, 0.6); f3->Draw("axis"); HLine(0.3, 1.1, 1, kGray + 2, 2, 1); }
                r->Draw("E1 P same");
                lg3->AddEntry(r, PresetName(kPres[ip]), "p");
            }
            for (double xl : {scan::kBumpLo, scan::kBumpHi, scan::kPeakLo, scan::kPeakHi}) VLine(xl, 0.4, 1.95, kGray + 1, 3, 1);
            lg3->Draw();
            Txt(0.07, 0.93, Form("(c) data (N - CB) / simulation (signal + phys. bkg) per selection, 10 MeV bins; %s; no direct #eta in sim", kLdataNote), 0.050);
            c->SaveAs(FigPath("spectra_m4_presets.png").c_str());
            c->SaveAs(FigPath("spectra_m4_presets.pdf").c_str());
        }
    }

    // ======================= Part 8: tables =======================
    if (parts & 8) {
        // ---- windows ----
        FILE* o = fopen(TabPath("spectra_windows.txt").c_str(), "w");
        fprintf(o, "# Cross sections integrated in the data windows, simulation (std + lep owned rows), e+e- hypothesis\n");
        fprintf(o, "# bump window %.4f-%.4f, peak window %.5f-%.5f GeV/c2 (measured, corrected M4)\n", scan::kBumpLo, scan::kBumpHi, scan::kPeakLo, scan::kPeakHi);
        fprintf(o, "# sigma = sum w / L_sim, L_sim = %.2f nb^-1; errors: event-group variance; N_eff = (sum w)^2 / V\n", kLsim);
        fprintf(o, "# expected data rows = sum w * L_data / L_sim with L_data = %.1f nb^-1 (UNCONFIRMED); no efficiency correction\n", kLdata);
        fprintf(o, "# eta -> pi+pi-e+e- (direct) is NOT in the simulation (eta set missing)\n");
        Binning bw; bw.nbin = 2; bw.lo = 0; bw.hi = 2;
        bw.x = [](const Row& r) {
            if (r.m4 >= scan::kBumpLo && r.m4 < scan::kBumpHi) return 0.5;
            if (r.m4 >= scan::kPeakLo && r.m4 < scan::kPeakHi) return 1.5;
            return -1.0;
        };
        for (int ip = 0; ip < kNP; ++ip) {
            const int p = kPres[ip];
            const Sel sel = SelPreset(p);
            fprintf(o, "\n## selection \"%s\": %s\n", PresetName(p), PresetLabel(p));
            fprintf(o, "%-58s | %9s %-8s %7s %10s %-8s | %9s %-8s %7s %10s %-8s\n", "component (mode A, w)", "bump [nb]", "+- err", "N_eff",
                    "data rows", "+- err", "peak [nb]", "+- err", "N_eff", "data rows", "+- err");
            const std::vector<WSpec> s = FillMany(d, sel, comps, optA, bw);
            for (int k = 0; k < kNC; ++k) {
                fprintf(o, "%-58s |", kCD[k].txt);
                for (int b = 0; b < 2; ++b) {
                    const double S = s[k].s[b], V = s[k].v[b];
                    fprintf(o, " %9.4f %-8.4f %7.0f %10.0f %-8.0f |", S / kLsim, std::sqrt(V) / kLsim, V > 0 ? S * S / V : 0.,
                            S * kLdata / kLsim, std::sqrt(V) * kLdata / kLsim);
                }
                fprintf(o, "\n");
            }
            // variants of the sum and of the signal
            for (int v = 1; v < 4; ++v) {
                SimOpt ov; ov.mode = v / 2; ov.db = v % 2;
                const std::vector<WSpec> sv = FillMany(d, sel, {comps[cSum], comps[cT1], comps[cT2]}, ov, bw);
                for (int j = 0; j < 3; ++j) {
                    const char* nm = j == 0 ? "SUM" : j == 1 ? "signal T1" : "phys. bkg T2";
                    fprintf(o, "%-58s |", Form("variant %s: %s", v == 1 ? "A, w_db" : v == 2 ? "B, w" : "B, w_db", nm));
                    for (int b = 0; b < 2; ++b) {
                        const double S = sv[j].s[b], V = sv[j].v[b];
                        fprintf(o, " %9.4f %-8.4f %7.0f %10.0f %-8.0f |", S / kLsim, std::sqrt(V) / kLsim, V > 0 ? S * S / V : 0.,
                                S * kLdata / kLsim, std::sqrt(V) * kLdata / kLsim);
                    }
                    fprintf(o, "\n");
                }
            }
            // data in the same windows
            for (int b = 0; b < 2; ++b) {
                const RegRes& q = REG[ip][b == 0 ? 1 : 3];
                fprintf(o, "data %s window: N %.0f, CB %.0f, N - CB = %.0f +- %.0f rows -> %.4f +- %.4f nb (L_data %.1f, unconfirmed); data / sim SUM = %.3f +- %.3f\n",
                        b == 0 ? "bump" : "peak", q.N, q.CB, q.N - q.CB, q.dE * kLdata, q.dS, q.dE, kLdata, q.r, q.rE);
            }
        }
        fclose(o);
        // ---- regions ----
        o = fopen(TabPath("spectra_regions.txt").c_str(), "w");
        fprintf(o, "# data (N - CB) / L_data against the simulation (signal + physical background, mode A, w) in M4 regions\n");
        fprintf(o, "# data: L_data = %.1f nb^-1 (UNCONFIRMED), CB = 2 sqrt(N++ N--) in 10 MeV bins, no k factor, window sums with event-group variance\n", kLdata);
        fprintf(o, "# sim: L_sim = %.2f nb^-1, event-group variance; eta -> pi+pi-e+e- (direct) missing in the simulation\n", kLsim);
        for (int ip = 0; ip < kNP; ++ip) {
            fprintf(o, "\n## selection \"%s\"\n", PresetName(kPres[ip]));
            fprintf(o, "%-22s %15s %10s %10s %20s %20s %16s\n", "region", "M4 [GeV/c2]", "N", "CB", "data sigma [nb]", "sim sigma [nb]", "data / sim");
            for (int r = 0; r <= kNReg; ++r) {
                const RegRes& q = REG[ip][r];
                const double lo = r < kNReg ? kRegLo[r] : 0.30, hi = r < kNReg ? kRegHi[r] : 1.10;
                fprintf(o, "%-22s %7.4f-%-7.4f %10.0f %10.0f %10.4g +- %-7.2g %10.4g +- %-7.2g %7.3f +- %.3f\n", r < kNReg ? kRegName[r] : "all 0.30-1.10",
                        lo, hi, q.N, q.CB, q.dS, q.dE, q.sS, q.sE, q.r, q.rE);
            }
        }
        fclose(o);
        // ---- templates (all M4) ----
        o = fopen(TabPath("spectra_templates.txt").c_str(), "w");
        fprintf(o, "# templates and components at any M4 (std + lep owned rows), e+e- hypothesis; sigma = sum w / %.2f nb^-1; N_eff = (sum w)^2 / V\n", kLsim);
        for (int ip = 0; ip < kNP; ++ip) {
            fprintf(o, "\n## selection \"%s\"\n%-58s %9s %12s %10s %18s | %18s %18s %18s\n", PresetName(kPres[ip]), "component", "rows", "sum w", "N_eff",
                    "sigma A,w [nb]", "A,w_db", "B,w", "B,w_db");
            std::vector<WSum> sv[4];
            for (int v = 0; v < 4; ++v) { SimOpt ov; ov.mode = v / 2; ov.db = v % 2; sv[v] = SumMany(d, SelPreset(kPres[ip]), comps, ov); }
            for (int k = 0; k < kNC; ++k) {
                fprintf(o, "%-58s %9.0f %12.6g %10.0f %9.4g +- %-6.2g |", kCD[k].txt, sv[0][k].N, sv[0][k].S, sv[0][k].Neff(), sv[0][k].Xsec(), sv[0][k].XsecErr());
                for (int v = 1; v < 4; ++v) fprintf(o, " %9.4g +- %-6.2g", sv[v][k].Xsec(), sv[v][k].XsecErr());
                fprintf(o, "\n");
            }
        }
        fclose(o);
        printf("sim_spectra: tables written (%.1f s)\n", swAll.RealTime());
        swAll.Continue();
    }

    // ======================= Part 2: subsystems =======================
    if (parts & 2) {
        DataPK pk;
        pk.Build(dd);
        std::map<std::string, std::pair<double, double>> yr;
        FILE* tab = fopen(TabPath("spectra_subsystems.txt").c_str(), "w");
        fprintf(tab, "# subsystems: integrals inside the plotted range (sim sum T1+T2 mode A w; data N - CB / L_data, L_data = %.1f nb^-1 UNCONFIRMED)\n", kLdata);
        for (int ip = 0; ip < kNP; ++ip) {
            DrawSubFigure(d, dd, pk, ip, "pipi", 3, 2, yr, tab);
            DrawSubFigure(d, dd, pk, ip, "ee", 4, 2, yr, tab);
            DrawSubFigure(d, dd, pk, ip, "tracks", 3, 3, yr, tab);
            DrawSubFigure(d, dd, pk, ip, "four", 4, 2, yr, tab);
            printf("sim_spectra: subsystems of \"%s\" done (%.1f s)\n", PresetName(kPres[ip]), swAll.RealTime());
            swAll.Continue();
        }
        fclose(tab);
    }

    // ======================= Part 4: eta structures =======================
    if (parts & 4) {
        FILE* o = fopen(TabPath("spectra_eta.txt").c_str(), "w");
#define P(...) do { printf(__VA_ARGS__); fprintf(o, __VA_ARGS__); } while (0)
        P("# eta structures in the simulation (std set only: the signal classes exist only there), std + lep owned rows, mode A, w (PDG)\n");
        P("# fits: weighted least squares, variance of a bin = q * model (q = V/S of the class, event-group variance), iterated; no Poisson likelihood\n");
        P("# data values from the data study (reports engine.md, scan_bump.md), NOT recomputed: peak %.4f +- %.5f, sigma %.1f +- %.1f MeV, scale %.4f +- %.4f, M(pipi) scale %.4f\n",
          kDataPeak, kDataPeakErr, kDataSigma * 1e3, kDataSigmaErr * 1e3, kDataScale, kDataScaleErr, kDataScalePiPi);
        const Sel selPre = SelPreset(kPre);
        const Comp cBoth = {"bump (Dalitz + conversion)", [](const Row& r) { return r.tmpl == 1 && (r.cls == 2 || r.cls == 3); }};
        const Comp cPCng = {"peak conversion, no ghost lepton", [](const Row& r) { return r.tmpl == 1 && r.cls == 4 && !GhostLepton(r); }};
        std::vector<Comp> sc = {comps[cBD], comps[cBC], comps[cPC], cBoth, cPCng};

        // ---- 1. peak: reco and true M4 of the peak-conversion class, per selection ----
        P("\nA. Peak conversion class (eta -> pi+pi-gamma, photon converted): Gaussian core fit (range mu -+ 2 sigma, iterated), 1 MeV bins\n");
        P("   %-24s %-22s %6s %8s | %-20s %-16s %-14s %9s | mean, rms in 0.50-0.59\n", "selection", "sample", "rows", "N_eff", "mu [GeV/c2]",
          "sigma [MeV]", "chi2/ndf", "mu/m_eta");
        Binning bp; bp.nbin = 90; bp.lo = 0.50; bp.hi = 0.59;
        Binning bpt; bpt.nbin = 250; bpt.lo = 0.535; bpt.hi = 0.560; bpt.x = [](const Row& r) { return (double)r.m4_true; };
        GRes gPre, gPreNG, gTrue, gRatio;
        WSpec spPeakRatio;
        std::vector<GRes> gVar;
        WSpec spPeak, spPeakTrue, spPeakNG;
        for (int ip = 0; ip < kNP; ++ip) {
            const std::vector<WSpec> s = FillMany(d, SelPreset(kPres[ip]), {comps[cPC], cPCng}, optA, bp);
            for (int j = 0; j < 2; ++j) {
                const GRes g = FitGaussCore(s[j], 0.5445, 0.006);
                double mean, rms; MeanRms(s[j], 0.50, 0.59, mean, rms);
                P("   %-24s %-22s %6.0f %8.0f | %.5f +- %.5f   %5.2f +- %-5.2f  %5.1f/%-3d      %.4f +- %.4f | %.5f, %.2f MeV | Gauss area %.3f +- %.3f nb = %.0f +- %.0f rows; 0.50-0.59: %.3f nb = %.0f rows\n",
                  PresetName(kPres[ip]), j == 0 ? "reco, all" : "reco, no ghost lepton", s[j].N, s[j].NeffIn(), g.mu, g.emu, g.sig * 1e3, g.esig * 1e3,
                  g.chi2, g.ndf, g.mu / scan::kMeta, g.emu / scan::kMeta, mean, rms * 1e3, g.A / kLsim, g.eA / kLsim, g.A * kLdata / kLsim,
                  g.eA * kLdata / kLsim, s[j].S / kLsim, s[j].S * kLdata / kLsim);
                if (ip == 0 && j == 0) { gPre = g; spPeak = s[0]; }
                if (ip == 0 && j == 1) { gPreNG = g; spPeakNG = s[1]; }
            }
        }
        {
            spPeakTrue = Fill(d, selPre, comps[cPC], optA, bpt);
            // true M4 of the class: weighted median, mean, share within 1 MeV of m_eta (a Gaussian makes no sense: the core is < 0.1 MeV)
            std::vector<std::pair<double, double>> tv;
            double sw = 0, sx = 0, sNear = 0;
            for (const Row& r : d.tree(1, kUS)) {
                if (!Owned(r, optA) || !selPre(r) || r.tmpl != 1 || r.cls != 4) continue;
                tv.push_back(std::make_pair((double)r.m4_true, r.w)); sw += r.w; sx += r.w * r.m4_true;
                if (std::fabs(r.m4_true - scan::kMeta) < 0.001) sNear += r.w;
            }
            std::sort(tv.begin(), tv.end());
            double c = 0, med = 0;
            for (auto& q : tv) { c += q.second; if (c >= 0.5 * sw) { med = q.first; break; } }
            gTrue.mu = med;
            P("   %-24s %-22s %6zu          | TRUE M4: weighted median %.5f, mean %.5f; share within 1 MeV of m_eta = %.3f (a Gaussian core is < 0.1 MeV wide)\n",
              "pre", "TRUE M4", tv.size(), med, sx / sw, sNear / sw);
            Binning brat; brat.nbin = 80; brat.lo = 0.95; brat.hi = 1.03;
            brat.x = [](const Row& r) { return r.m4_true > 0 ? (double)r.m4 / r.m4_true : -1.0; };
            spPeakRatio = Fill(d, selPre, comps[cPC], optA, brat);
            gRatio = FitGaussCore(spPeakRatio, 0.994, 0.010);
            P("   %-24s %-22s %6.0f %8.0f | per-row m4/m4_true, Gaussian core: %.5f +- %.5f, sigma %.5f +- %.5f, chi2/ndf %.1f/%d\n", "pre", "reco / true",
              spPeakRatio.N, spPeakRatio.NeffIn(), gRatio.mu, gRatio.emu, gRatio.sig, gRatio.esig, gRatio.chi2, gRatio.ndf);
            P("   fit-range variation (pre, reco, all): ");
            for (double ns : {1.5, 2.0, 2.5, 3.0}) {
                const GRes g = FitGaussCore(spPeak, 0.5445, 0.006, ns);
                gVar.push_back(g);
                P(" +-%.1f sigma: %.5f / %.2f MeV;", ns, g.mu, g.sig * 1e3);
            }
            P("\n   mass scale reco/true (pre): fitted mu_reco / median true = %.5f +- %.5f; mu_reco / m_eta = %.5f +- %.5f; per-row core %.5f +- %.5f; data: %.4f +- %.4f\n",
              gPre.mu / gTrue.mu, gPre.emu / gTrue.mu, gPre.mu / scan::kMeta, gPre.emu / scan::kMeta, gRatio.mu, gRatio.emu, kDataScale, kDataScaleErr);
            P("   sim - data: position %+.2f +- %.2f MeV, width %+.2f +- %.2f MeV\n", (gPre.mu - kDataPeak) * 1e3,
              std::hypot(gPre.emu, kDataPeakErr) * 1e3, (gPre.sig - kDataSigma) * 1e3, std::hypot(gPre.esig, kDataSigmaErr) * 1e3);
        }
        // ---- 2. resolutions and per-row ratios ----
        P("\nB. Resolution (pre): Gaussian core of reco - true (range -+2 sigma, iterated) and weighted rms in -+50 MeV (M4) / -+20 MeV (M(pipi))\n");
        P("   %-28s | %-17s %-15s %9s %9s | %-17s %-15s %9s %9s | %9s %9s\n", "class", "M4: core mean", "core sigma", "mean", "rms", "M(pipi): core mean",
          "core sigma", "mean", "rms", "<m4/true>", "<mpp/true>");
        Binning bres; bres.nbin = 100; bres.lo = -0.05; bres.hi = 0.05;
        bres.x = [](const Row& r) { return r.m4_true > 0 ? (double)r.m4 - r.m4_true : -999.; };
        Binning bresp; bresp.nbin = 80; bresp.lo = -0.02; bresp.hi = 0.02;
        bresp.x = [](const Row& r) { return r.mpipi_true > 0 ? (double)r.mpipi - r.mpipi_true : -999.; };
        const std::vector<WSpec> rM4 = FillMany(d, selPre, sc, optA, bres);
        const std::vector<WSpec> rPP = FillMany(d, selPre, sc, optA, bresp);
        std::vector<GRes> gM4(sc.size()), gPP(sc.size());
        double ratM4[8], ratPP[8];
        for (size_t j = 0; j < sc.size(); ++j) {
            gM4[j] = FitGaussCore(rM4[j], -0.003, 0.005);
            gPP[j] = FitGaussCore(rPP[j], 0.0, 0.002);
            double m1, r1, m2, r2;
            MeanRms(rM4[j], -0.05, 0.05, m1, r1); MeanRms(rPP[j], -0.02, 0.02, m2, r2);
            double sw = 0, sr = 0, swp = 0, srp = 0;
            for (int set = 1; set <= 2; ++set)
                for (const Row& r : d.tree(set, kUS)) {
                    if (!Owned(r, optA) || !selPre(r) || !sc[j].f(r)) continue;
                    if (r.m4_true > 0) { sw += r.w; sr += r.w * r.m4 / r.m4_true; }
                    if (r.mpipi_true > 0) { swp += r.w; srp += r.w * r.mpipi / r.mpipi_true; }
                }
            ratM4[j] = sr / sw; ratPP[j] = srp / swp;
            P("   %-28s | %+7.2f +- %-5.2f %5.2f +- %-5.2f %+9.2f %9.2f | %+7.2f +- %-5.2f %5.2f +- %-5.2f %+9.2f %9.2f | %9.5f %9.5f   [MeV]\n", sc[j].name.c_str(),
              gM4[j].mu * 1e3, gM4[j].emu * 1e3, gM4[j].sig * 1e3, gM4[j].esig * 1e3, m1 * 1e3, r1 * 1e3, gPP[j].mu * 1e3, gPP[j].emu * 1e3,
              gPP[j].sig * 1e3, gPP[j].esig * 1e3, m2 * 1e3, r2 * 1e3, ratM4[j], ratPP[j]);
        }
        // M(pipi) scale: signal rows and the K0S analogue of the data method
        double sPiPiSim = 0;
        {
            double sw = 0, sr = 0;
            for (const Row& r : d.tree(1, kUS)) {
                if (!Owned(r, optA) || !selPre(r) || r.tmpl != 1 || r.mpipi_true <= 0) continue;
                sw += r.w; sr += r.w * r.mpipi / r.mpipi_true;
            }
            sPiPiSim = sr / sw;
            const Comp ck = {"K0S pions", [](const Row& r) { return r.pcls == 6; }};
            const Sel sk = SelAnd(selPre, [](const Row& r) { return r.dz_pipi > 20.f && r.doca_pipi < 10.f; });
            Binning bk; bk.nbin = 120; bk.lo = 0.47; bk.hi = 0.53; bk.x = [](const Row& r) { return (double)r.mpipi; };
            const WSpec k0 = Fill(d, sk, ck, optA, bk);
            const GRes gk = FitGaussCore(k0, 0.4975, 0.003);
            double sw2 = 0, sr2 = 0;
            for (int set = 1; set <= 2; ++set)
                for (const Row& r : d.tree(set, kUS)) {
                    if (!Owned(r, optA) || !sk(r) || r.pcls != 6 || r.mpipi_true <= 0) continue;
                    sw2 += r.w; sr2 += r.w * r.mpipi / r.mpipi_true;
                }
            P("\n   M(pipi) scale of the simulation: signal rows <mpipi/mpipi_true> = %.5f; K0S pions (pre, dz_pipi > 20, doca_pipi < 10, any template, N_eff %.0f):\n",
              sPiPiSim, k0.NeffIn());
            P("   Gaussian %.5f +- %.5f, sigma %.2f MeV -> scale %.5f +- %.5f (<mpipi/true> %.5f); data (K0S): %.4f +- 0.0004\n", gk.mu, gk.emu, gk.sig * 1e3,
              gk.mu / scan::kMK0, gk.emu / scan::kMK0, sr2 / sw2, kDataScalePiPi);
        }
        const double kSimScale = gPre.mu / scan::kMeta;

        // ---- 3. bump edge per M(pipi) slice ----
        const int nsl = 5;
        const double slLo[nsl] = {0.30, 0.34, 0.37, 0.39, 0.39}, slHi[nsl] = {0.34, 0.37, 0.39, 0.41, 0.40};
        std::vector<WSpec> eReco[nsl], eTrue[nsl];
        ERes fR[nsl][3], fRfix[nsl], fT[nsl][3];
        P("\nC. Falling edge of the bump per M(pipi) slice (pre). Model: A (1 + B (m - e_c)/0.04) x step, the step averaged over the slice as in the data\n");
        P("   study (scan_bump.C): 0.5 erfc((m - k M4max(M(pipi)/s_pipi)) / (sqrt2 s)), M4max^2 = m_eta^2 - m_eta (E0 - p0), k and s free.\n");
        P("   reco: measured M4 in slices of measured M(pipi), s_pipi = %.5f (simulation), 2 MeV bins, range k0 e_c - 40 .. + 25 MeV (k0 = %.4f);\n", sPiPiSim, kSimScale);
        P("   true: m4_true in slices of mpipi_true, s_pipi = 1, 1 MeV bins, range e_c - 30 .. + 12 MeV. Edge = k M4max(slice centre / s_pipi).\n");
        P("   A and A B solved linearly; k, s on a grid (coarse, then fine); errors: delta chi2 = 1 of the profile; a NEGATIVE error means the interval reaches the grid edge.\n");
        P("   predictions at the slice centre c: true M4max(c); data scales %.4f M4max(c/%.4f); sim scales %.5f M4max(c/%.5f)\n", kDataScale, kDataScalePiPi,
          kSimScale, sPiPiSim);
        const std::vector<Comp> ec = {cBoth, comps[cBD], comps[cBC]};
        for (int i = 0; i < nsl; ++i) {
            const double a = slLo[i], b = slHi[i];
            const EdgeModel emR(a, b, sPiPiSim), emT(a, b, 1.0);
            Binning br; br.nbin = 50; br.lo = kSimScale * emR.ec - 0.06; br.hi = br.lo + 0.10;
            Binning bt; bt.nbin = 70; bt.lo = emT.ec - 0.05; bt.hi = bt.lo + 0.07; bt.x = [](const Row& r) { return (double)r.m4_true; };
            eReco[i] = FillMany(d, SelAnd(selPre, [a, b](const Row& r) { return r.mpipi >= a && r.mpipi < b; }), ec, optA, br);
            eTrue[i] = FillMany(d, SelAnd(selPre, [a, b](const Row& r) { return r.mpipi_true >= a && r.mpipi_true < b; }), ec, optA, bt);
            for (int j = 0; j < 3; ++j) {
                fR[i][j] = FitEdge(eReco[i][j], EdgeModel(a, b, sPiPiSim), kSimScale, 0.0015, 0.020, 0.040, 0.025);
                fT[i][j] = FitEdge(eTrue[i][j], EdgeModel(a, b, 1.0), 1.0, 0.0001, 0.006, 0.030, 0.012);
                fR[i][j].S = eReco[i][j].S; fR[i][j].Neff = eReco[i][j].NeffIn();
                fT[i][j].S = eTrue[i][j].S; fT[i][j].Neff = eTrue[i][j].NeffIn();
            }
            fRfix[i] = FitEdge(eReco[i][0], EdgeModel(a, b, sPiPiSim), kSimScale, 0.0015, 0.020, 0.040, 0.025, true, kDataSigma);
        }
        P("   %-10s %-30s | %8s | %-18s %-8s %-16s %-11s | %-18s | %-18s %-8s %-9s | %8s %8s %8s | %-18s\n", "M(pipi)", "class", "N_eff", "reco edge (s free)",
          "k", "", "s [MeV]", "reco edge (s=5.4)", "true edge", "k_true", "s [MeV]", "pr.true", "pr.data", "pr.sim", "data fitted edge");
        for (int i = 0; i < nsl; ++i) {
            const double cen = 0.5 * (slLo[i] + slHi[i]);
            const double pT = scan::M4max(cen), pD = kDataScale * scan::M4max(cen / kDataScalePiPi), pS = kSimScale * scan::M4max(cen / sPiPiSim);
            int di = -1;
            for (int k = 0; k < kNDataEdge; ++k) if (std::fabs(kDataEdgeLo[k] - slLo[i]) < 1e-6 && std::fabs(kDataEdgeHi[k] - slHi[i]) < 1e-6) di = k;
            for (int j = 0; j < 3; ++j) {
                const ERes& r = fR[i][j];
                const ERes& t = fT[i][j];
                std::string dstr = "-";
                if (di >= 0 && j == 0) {
                    const double e = scan::M4max(cen / kDataScalePiPi);
                    dstr = Form("%.4f -%.4f +%.4f", kDataEdge[di], kDataEdgeFm[di] * e, kDataEdgeFp[di] * e);
                } else if (j == 0) dstr = "no data fit";
                P("   %.2f-%.2f  %-30s | %8.0f | %.4f +- %.4f  %.4f +- %.4f    %5.2f +- %-4.2f | %s | %.4f +- %.4f  %.4f   %5.2f    | %8.4f %8.4f %8.4f | %s   (chi2/ndf reco %.1f/%d, true %.1f/%d)\n",
                  slLo[i], slHi[i], ec[j].name.c_str(), r.Neff, r.edge, r.eedge, r.k, r.ek, r.s * 1e3, r.es * 1e3,
                  j == 0 ? Form("%.4f +- %.4f", fRfix[i].edge, fRfix[i].eedge) : "        -         ", t.edge, t.eedge, t.k, t.s * 1e3, pT, pD, pS,
                  dstr.c_str(), r.chi2, r.ndf, t.chi2, t.ndf);
            }
        }
        {
            double sw = 0, sx = 0;
            for (int i = 0; i < 3; ++i) { const double w = 1.0 / std::pow(std::max(std::fabs(fR[i][0].ek), 1e-5), 2); sw += w; sx += w * fR[i][0].k; }
            const double wi = 1.0 / std::pow(std::max(std::fabs(fR[3][0].ek), 1e-5), 2);
            P("   weighted mean edge factor k (reco, bump both classes): slices 1-3: %.5f +- %.5f; slices 1-4 (0.39-0.41): %.5f +- %.5f; sim peak scale %.5f +- %.5f\n",
              sx / sw, 1 / std::sqrt(sw), (sx + wi * fR[3][0].k) / (sw + wi), 1 / std::sqrt(sw + wi), kSimScale, gPre.emu / scan::kMeta);
            double tw = 0, tx = 0;
            for (int i = 0; i < 4; ++i) { const double w = 1.0 / std::pow(std::max(std::fabs(fT[i][0].ek), 1e-5), 2); tw += w; tx += w * fT[i][0].k; }
            P("   weighted mean k of the TRUE spectra (slices 1-4): %.5f +- %.5f (1 = the formula holds)\n", tx / tw, 1 / std::sqrt(tw));
            P("   data (scan_bump.md, no cut): factors 0.9912, 0.9932, 0.9927, 0.9972 (0.39-0.40), weighted mean 0.9930 +- 0.0013; data peak scale 0.9928 +- 0.0009\n");
        }
        // ---- row-by-row aligned edge: u = M4 / M4max(M(pipi)/s_pipi), all bump rows with 0.30 <= M(pipi) < 0.41 ----
        const double sPP = sPiPiSim;
        Binning bu; bu.nbin = 100; bu.lo = 0.90; bu.hi = 1.05;
        bu.x = [sPP](const Row& r) { const double e = scan::M4max(r.mpipi / sPP); return e > 0 ? r.m4 / e : -1.0; };
        Binning but; but.nbin = 120; but.lo = 0.92; but.hi = 1.04;
        but.x = [](const Row& r) { const double e = scan::M4max(r.mpipi_true); return (e > 0 && r.m4_true > 0) ? r.m4_true / e : -1.0; };
        const std::vector<WSpec> uR = FillMany(d, SelAnd(selPre, [](const Row& r) { return r.mpipi >= 0.30f && r.mpipi < 0.41f; }), ec, optA, bu);
        const std::vector<WSpec> uT = FillMany(d, SelAnd(selPre, [](const Row& r) { return r.mpipi_true >= 0.30f && r.mpipi_true < 0.41f; }), ec, optA, but);
        ERes aR[3], aT[3], aRfix;
        for (int j = 0; j < 3; ++j) {
            aR[j] = FitEdge(uR[j], EdgeModel::Aligned(), kSimScale, 0.003, 0.04, 0.07, 0.04);
            aT[j] = FitEdge(uT[j], EdgeModel::Aligned(), 1.0, 0.0002, 0.012, 0.05, 0.02);
        }
        aRfix = FitEdge(uR[0], EdgeModel::Aligned(), kSimScale, 0.003, 0.04, 0.07, 0.04, true, kDataSigma / 0.515);
        P("\nD. Row-by-row aligned edge (all bump rows, 0.30 <= M(pipi) < 0.41): u = M4 / M4max(M(pipi)/%.5f) (reco), u = m4_true / M4max(mpipi_true) (true);\n", sPiPiSim);
        P("   step at u = k, smearing s in u; same model and grid as C (slope over 0.04 in u)\n");
        P("   %-30s | %8s | %-22s %-18s %-10s | %-22s %-18s %-10s\n", "class", "N_eff", "reco k", "s (u units)", "chi2/ndf", "true k", "s (u units)", "chi2/ndf");
        for (int j = 0; j < 3; ++j)
            P("   %-30s | %8.0f | %.5f +- %-10.5f   %.4f +- %-8.4f %5.1f/%-3d | %.5f +- %-10.5f   %.4f +- %-8.4f %5.1f/%-3d\n", ec[j].name.c_str(),
              uR[j].NeffIn(), aR[j].k, aR[j].ek, aR[j].s, aR[j].es, aR[j].chi2, aR[j].ndf, aT[j].k, aT[j].ek, aT[j].s, aT[j].es, aT[j].chi2, aT[j].ndf);
        P("   bump both, reco, s fixed to 5.4 MeV (0.0105 in u, the data choice): k = %.5f +- %.5f, chi2/ndf %.1f/%d\n", aRfix.k, aRfix.ek, aRfix.chi2, aRfix.ndf);
        P("   compare: sim peak scale %.5f +- %.5f; data edge factor 0.9930 +- 0.0013 and data peak scale 0.9928 +- 0.0009 (data use s_pipi = 0.9970)\n",
          kSimScale, gPre.emu / scan::kMeta);
        fclose(o);
#undef P

        // ======== figures of the eta structures ========
        // ---- (i) reco and true M4 of the signal classes, resolutions ----
        {
            TCanvas* c = new TCanvas("cEtaRT", "", 2100, 1500);
            Binning bb; bb.nbin = 100; bb.lo = 0.35; bb.hi = 0.60;
            Binning bbt = bb; bbt.x = [](const Row& r) { return (double)r.m4_true; };
            Binning bpk; bpk.nbin = 80; bpk.lo = 0.50; bpk.hi = 0.58;
            Binning bpkt = bpk; bpkt.x = [](const Row& r) { return (double)r.m4_true; };
            const std::vector<WSpec> rb = FillMany(d, selPre, {comps[cBD], comps[cBC]}, optA, bb);
            const std::vector<WSpec> tb = FillMany(d, selPre, {comps[cBD], comps[cBC]}, optA, bbt);
            const WSpec rp = Fill(d, selPre, comps[cPC], optA, bpk), tp = Fill(d, selPre, comps[cPC], optA, bpkt);
            const char* cn[3] = {"bump, #pi^{0} Dalitz", "bump, conversion", "peak, conversion"};
            for (int j = 0; j < 3; ++j) {
                c->cd();
                MkPad(Form("rt%d", j), j / 3.0, 0.5, (j + 1) / 3.0, 0.98, j == 2, 0.17, 0.03, 0.14, 0.13);
                TH1D* hr = j < 2 ? rb[j].Hist(Form("rtR%d", j), 2) : rp.Hist("rtR2", 2);
                TH1D* ht = j < 2 ? tb[j].Hist(Form("rtT%d", j), 2) : tp.Hist("rtT2", 2);
                StyleLine(hr, gLine[kSigOrder[j]], 3, 1); StyleLine(ht, kGray + 3, 2, 2);
                TH1D* fr = Frame(hr, Form("frRT%d", j), j == 2 ? 0.3 : 0, std::max(MaxOf(hr, true), MaxOf(ht)) * (j == 2 ? 20 : 1.25));
                Axes(fr, kM4t, kM4y, 0.050, 1.6);
                fr->Draw("axis");
                ht->Draw("hist same"); hr->Draw("hist same");
                TH1D* he = (TH1D*)hr->Clone(Form("rtE%d", j)); StyleLine(he, gLine[kSigOrder[j]], 1, 1); he->Draw("E1 same");
                const double ym = fr->GetMaximum();
                if (j < 2) {
                    VLine(scan::kBumpLo, 0, ym, kGray + 1, 3, 1); VLine(scan::kBumpHi, 0, ym, kGray + 1, 3, 1);
                } else {
                    VLine(kDataPeak, 0, ym, 4, 2, 2); VLine(scan::kMeta, 0, ym, kGray + 2, 3, 1);
                }
                TLegend* lg = new TLegend(0.19, 0.68, 0.97, 0.80);
                lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.036);
                lg->AddEntry(hr, Form("reco (corrected): %.3g nb", (j < 2 ? rb[j].S : rp.S) / kLsim), "l");
                lg->AddEntry(ht, Form("true (generator momenta): %.3g nb", (j < 2 ? tb[j].S : tp.S) / kLsim), "l");
                lg->Draw();
                Txt(0.17, 0.93, Form("(%c) %s, pre, %s bins", 'a' + j, cn[j], j < 2 ? "2.5 MeV" : "1 MeV"), 0.050);
                Txt(0.17, 0.885, j < 2 ? "dotted: bump window of the data" : "log y; blue dashed: data peak 0.5439; grey: m_{#eta}", 0.036, j < 2 ? 1 : 4);
            }
            // resolutions
            for (int q = 0; q < 2; ++q) {
                c->cd();
                MkPad(Form("rs%d", q), q / 3.0, 0.0, (q + 1) / 3.0, 0.48, true, 0.17, 0.03, 0.14, 0.13);
                const std::vector<WSpec>& R = q == 0 ? rM4 : rPP;
                TH1D* f = nullptr;
                TLegend* lg = new TLegend(0.19, 0.64, 0.97, 0.80);
                lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.031);
                for (int j = 0; j < 3; ++j) {
                    TH1D* h = R[j].Hist(Form("res%d_%d", q, j), 0);
                    if (R[j].S > 0) h->Scale(1.0 / R[j].S);
                    StyleLine(h, gLine[kSigOrder[j]], 3, 1);
                    if (!f) {
                        f = Frame(h, Form("frRes%d", q), 2e-4, 1.5);
                        Axes(f, q == 0 ? "M_{4}^{reco} - M_{4}^{true} [GeV/c^{2}]" : "M_{#pi#pi}^{reco} - M_{#pi#pi}^{true} [GeV/c^{2}]", "fraction per bin", 0.050, 1.5);
                        f->GetXaxis()->SetNdivisions(505);
                        f->Draw("axis");
                    }
                    h->Draw("hist same");
                    const GRes& g = q == 0 ? gM4[j] : gPP[j];
                    lg->AddEntry(h, Form("%s: core %+.2f, #sigma %.2f MeV", cn[j], g.mu * 1e3, g.sig * 1e3), "l");
                }
                lg->Draw();
                Txt(0.17, 0.93, q == 0 ? "(d) M_{4} resolution, signal classes, pre" : "(e) M_{#pi#pi} resolution, signal classes, pre", 0.046);
                Txt(0.17, 0.885, "area 1 per class; core: Gaussian fit #pm2#sigma (iterated)", 0.036);
            }
            // M(pipi) reco and true of the signal classes
            {
                c->cd();
                MkPad("rs2", 2 / 3.0, 0.0, 1.0, 0.48, false, 0.17, 0.03, 0.14, 0.13);
                Binning bm; bm.nbin = 56; bm.lo = 0.27; bm.hi = 0.55; bm.x = [](const Row& r) { return (double)r.mpipi; };
                Binning bmt = bm; bmt.x = [](const Row& r) { return (double)r.mpipi_true; };
                const std::vector<WSpec> mr = FillMany(d, selPre, {comps[cBD], comps[cBC], comps[cPC]}, optA, bm);
                const std::vector<WSpec> mt = FillMany(d, selPre, {comps[cBD], comps[cBC], comps[cPC]}, optA, bmt);
                TH1D* f = nullptr;
                double ym = 0;
                std::vector<TH1D*> hh;
                for (int j = 0; j < 3; ++j) {
                    TH1D* h = mr[j].Hist(Form("mpr%d", j), 2); TH1D* ht = mt[j].Hist(Form("mpt%d", j), 2);
                    StyleLine(h, gLine[kSigOrder[j]], 3, 1); StyleLine(ht, gLine[kSigOrder[j]], 2, 2);
                    ym = std::max(ym, std::max(MaxOf(h), MaxOf(ht)));
                    hh.push_back(h); hh.push_back(ht);
                }
                f = Frame(hh[0], "frMP", 0, ym * 1.35);
                Axes(f, "M_{#pi^{+}#pi^{-}} [GeV/c^{2}]", "d#sigma/dM_{#pi#pi} [nb/(GeV/c^{2})]", 0.050, 1.6);
                f->Draw("axis");
                for (TH1D* h : hh) h->Draw("hist same");
                TLegend* lg = new TLegend(0.19, 0.64, 0.97, 0.80);
                lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.040);
                lg->AddEntry(hh[0], "reco (solid), per class colour", "l");
                lg->AddEntry(hh[1], "true (dashed)", "l");
                lg->Draw();
                Txt(0.17, 0.93, "(f) M_{#pi^{+}#pi^{-}} of the signal classes, pre, 5 MeV", 0.050);
                Txt(0.17, 0.885, "orange: bump #pi^{0} Dal., magenta: bump conv., green: peak conv.", 0.032);
            }
            c->cd();
            Txt(0.005, 0.985, (kSimLine + "; " + kDirect).c_str(), 0.0115);
            c->SaveAs(FigPath("spectra_eta_reco_true.png").c_str());
            c->SaveAs(FigPath("spectra_eta_reco_true.pdf").c_str());
        }
        // ---- (ii) peak fit ----
        {
            TCanvas* c = new TCanvas("cEtaPk", "", 2100, 800);
            for (int q = 0; q < 3; ++q) {
                c->cd();
                MkPad(Form("pk%d", q), q / 3.0, 0.0, (q + 1) / 3.0, 0.95, false, 0.17, 0.03, 0.20, 0.13);
                const WSpec& s = q == 0 ? spPeak : q == 1 ? spPeakNG : spPeakRatio;
                const GRes& g = q == 0 ? gPre : q == 1 ? gPreNG : gRatio;
                TH1D* h = s.Hist(Form("pkh%d", q), 2);
                StyleData(h, gLine[cPC], 20, 0.9);
                h->SetLineWidth(2);
                TH1D* f = Frame(h, Form("frPk%d", q), 0, MaxOf(h, true) * 1.45);
                Axes(f, q < 2 ? kM4t : "M_{4}^{reco} / M_{4}^{true} (per row)", q < 2 ? kM4y : "d#sigma/d(ratio) [nb]", 0.050, 1.6);
                f->Draw("axis");
                // fitted Gaussian (dsigma/dx units) and the data Gaussian with the same area
                const double conv = 1.0 / (kLsim * s.bw);
                const GRes gg = g;
                TGraph* gf = Curve([gg](double x0, double bw) { const double p[3] = {gg.A, gg.mu, gg.sig}; return GaussBin(x0, bw, p); },
                                   s.bw, gg.lo, gg.hi, conv, kRed + 1, 3);
                if (q < 2) {
                    TGraph* gd = Curve([gg](double x0, double bw) { const double p[3] = {gg.A, kDataPeak, kDataSigma}; return GaussBin(x0, bw, p); },
                                       s.bw, kDataPeak - 0.02, kDataPeak + 0.02, conv, 4, 3);
                    gd->SetLineStyle(2);
                    gd->Draw("L same");
                }
                h->Draw("E1 P same");
                gf->Draw("L same");
                if (q < 2) VLine(scan::kMeta, 0, f->GetMaximum(), kGray + 2, 3, 1);
                else { VLine(kDataScale, 0, f->GetMaximum(), 4, 2, 2); VLine(1.0, 0, f->GetMaximum(), kGray + 2, 3, 1); }
                const char* tt[3] = {"(a) peak conversion, reco M_{4}, pre, 1 MeV", "(b) the same without ghost-lepton rows",
                                     "(c) peak conversion: reco / true M_{4} per row"};
                Txt(0.17, 0.93, tt[q], 0.044);
                Txt(0.19, 0.75, Form("Gaussian core: #mu = %.5f #pm %.5f", g.mu, g.emu), 0.034, kRed + 1);
                if (q < 2) Txt(0.19, 0.71, Form("#sigma = %.2f #pm %.2f MeV, #chi^{2}/ndf = %.1f/%d", g.sig * 1e3, g.esig * 1e3, g.chi2, g.ndf), 0.034, kRed + 1);
                else Txt(0.19, 0.71, Form("#sigma = %.4f #pm %.4f, #chi^{2}/ndf = %.1f/%d", g.sig, g.esig, g.chi2, g.ndf), 0.034, kRed + 1);
                if (q < 2) {
                    Txt(0.19, 0.67, Form("#mu/m_{#eta} = %.4f #pm %.4f (data %.4f #pm %.4f)", g.mu / scan::kMeta, g.emu / scan::kMeta, kDataScale, kDataScaleErr), 0.034, kRed + 1);
                    Txt(0.19, 0.63, Form("blue dashed: data %.4f, #sigma %.1f MeV (same area)", kDataPeak, kDataSigma * 1e3), 0.034, 4);
                    Txt(0.19, 0.59, Form("rows %.0f, N_{eff} %.0f, %.3g #pm %.2g nb in 0.50-0.59", s.N, s.NeffIn(), s.S / kLsim, std::sqrt(s.V) / kLsim), 0.034);
                } else {
                    Txt(0.19, 0.67, Form("true M_{4}: median %.5f (m_{#eta} = 0.547862)", gTrue.mu), 0.034);
                    Txt(0.19, 0.63, Form("blue dashed: data scale %.4f", kDataScale), 0.034, 4);
                }
                Txt(0.17, 0.875, "weighted LSQ, variance q #times model (event groups)", 0.036);
            }
            c->cd();
            Txt(0.005, 0.975, (kSimLine + "; " + kDirect).c_str(), 0.019);
            c->SaveAs(FigPath("spectra_eta_peak.png").c_str());
            c->SaveAs(FigPath("spectra_eta_peak.pdf").c_str());
        }
        // ---- (iii) edge fits per slice ----
        {
            TCanvas* c = new TCanvas("cEtaEdge", "", 2500, 1300);
            for (int row = 0; row < 2; ++row)
                for (int i = 0; i < nsl; ++i) {
                    c->cd();
                    MkPad(Form("ed%d_%d", row, i), i / 5.0, row == 0 ? 0.47 : 0.0, (i + 1) / 5.0, row == 0 ? 0.94 : 0.47, false, 0.17, 0.03, 0.17, 0.13);
                    const WSpec& s = row == 0 ? eReco[i][0] : eTrue[i][0];
                    const ERes& r = row == 0 ? fR[i][0] : fT[i][0];
                    TH1D* h = s.Hist(Form("edh%d_%d", row, i), 2);
                    StyleData(h, gLine[cBD], 20, 0.8);
                    TH1D* hD = (row == 0 ? eReco[i][1] : eTrue[i][1]).Hist(Form("edD%d_%d", row, i), 2);
                    TH1D* hC = (row == 0 ? eReco[i][2] : eTrue[i][2]).Hist(Form("edC%d_%d", row, i), 2);
                    StyleLine(hD, gLine[cBD], 2, 2); StyleLine(hC, gLine[cBC], 2, 2);
                    const double ym = MaxOf(h, true) * 1.55;
                    TH1D* f = Frame(h, Form("frEd%d_%d", row, i), 0, ym);
                    Axes(f, row == 0 ? kM4t : "M_{#pi^{+}#pi^{-}e^{+}e^{-}}^{true} [GeV/c^{2}]", kM4y, 0.055, 1.45);
                    f->GetXaxis()->SetNdivisions(505);
                    f->Draw("axis");
                    hD->Draw("hist same"); hC->Draw("hist same");
                    h->Draw("E1 P same");
                    // fitted model
                    const EdgeModel em(slLo[i], slHi[i], row == 0 ? sPiPiSim : 1.0);
                    const ERes rr = r;
                    TGraph* gf = Curve([em, rr](double x0, double bw) { const double p[4] = {rr.A, rr.B, rr.k, rr.s}; return em(x0, bw, p); },
                                       s.bw, rr.lo, rr.hi, 1.0 / (kLsim * s.bw), kRed + 1, 3);
                    gf->Draw("L same");
                    const double cen = 0.5 * (slLo[i] + slHi[i]);
                    VLine(r.edge, 0, ym * 0.80, kRed + 1, 1, 2);
                    VLine(scan::M4max(cen), 0, ym * 0.80, kGray + 2, 3, 2);
                    if (row == 0) {
                        VLine(kDataScale * scan::M4max(cen / kDataScalePiPi), 0, ym * 0.80, 4, 2, 2);
                        for (int k = 0; k < kNDataEdge; ++k)
                            if (std::fabs(kDataEdgeLo[k] - slLo[i]) < 1e-6 && std::fabs(kDataEdgeHi[k] - slHi[i]) < 1e-6) VLine(kDataEdge[k], 0, ym * 0.80, 4, 1, 3);
                    }
                    Txt(0.17, 0.93, Form("%s, M_{#pi#pi}%s %.2f-%.2f", row == 0 ? "reco" : "TRUE", row == 0 ? "" : "^{true}", slLo[i], slHi[i]), 0.055);
                    Txt(0.19, 0.76, Form("edge %.4f #pm %.4f", r.edge, r.eedge), 0.044, kRed + 1);
                    Txt(0.19, 0.71, Form("k = %.4f #pm %.4f, s = %.1f MeV", r.k, std::fabs(r.ek), r.s * 1e3), 0.044, kRed + 1);
                    Txt(0.19, 0.66, Form("#chi^{2}/ndf %.1f/%d, N_{eff} %.0f", r.chi2, r.ndf, s.NeffIn()), 0.044);
                }
            c->cd();
            Txt(0.005, 0.975, "Bump edge in the simulation (pre): markers = bump Dalitz + conversion, dashed = each class (orange Dalitz, magenta conversion), red = fitted step model (k, s free) and its edge;", 0.0125);
            Txt(0.005, 0.955, Form("grey dotted = M4max(centre) (true formula), blue dashed = prediction with the data scales (%.4f, %.4f), blue solid = edge fitted in the data (scan_bump, no cut; only 0.30-0.34, 0.34-0.37, 0.37-0.39, 0.39-0.40)",
                                   kDataScale, kDataScalePiPi), 0.0125, 4);
            c->SaveAs(FigPath("spectra_eta_edge.png").c_str());
            c->SaveAs(FigPath("spectra_eta_edge.pdf").c_str());
        }
        // ---- (iv) edge summary ----
        {
            TCanvas* c = new TCanvas("cEtaEdgeS", "", 2700, 850);
            c->cd();
            MkPad("es1", 0.0, 0.0, 1 / 3.0, 0.92, false, 0.16, 0.03, 0.06, 0.13);
            TH1D* f = new TH1D("frEs1", "", 100, 0.29, 0.42);
            f->SetMinimum(0.493); f->SetMaximum(0.530);
            Axes(f, "M_{#pi^{+}#pi^{-}} slice centre [GeV/c^{2}]", "edge position [GeV/c^{2}]", 0.050, 1.55);
            f->Draw("axis");
            TH1D *cT = new TH1D("cT", "", 260, 0.29, 0.42), *cD = new TH1D("cDp", "", 260, 0.29, 0.42), *cS = new TH1D("cSp", "", 260, 0.29, 0.42);
            for (int b = 1; b <= 260; ++b) {
                const double x = cT->GetBinCenter(b);
                cT->SetBinContent(b, scan::M4max(x)); cD->SetBinContent(b, kDataScale * scan::M4max(x / kDataScalePiPi));
                cS->SetBinContent(b, kSimScale * scan::M4max(x / sPiPiSim));
            }
            StyleLine(cT, kGray + 2, 2, 3); StyleLine(cD, 4, 2, 2); StyleLine(cS, kRed + 1, 2, 2);
            cT->Draw("hist l same"); cD->Draw("hist l same"); cS->Draw("hist l same");
            TH1D *pR = new TH1D("pR", "", 1, 0, 1), *pT = new TH1D("pT", "", 1, 0, 1), *pD = new TH1D("pD", "", 1, 0, 1);
            StyleData(pR, kRed + 1, 20, 1.6); StyleData(pT, kGray + 2, 24, 1.6); StyleData(pD, 4, 21, 1.6);
            for (int i = 0; i < 4; ++i) {
                const double cen = 0.5 * (slLo[i] + slHi[i]), hw = 0.5 * (slHi[i] - slLo[i]);
                TLine* l = new TLine(cen - hw, fR[i][0].edge, cen + hw, fR[i][0].edge); l->SetLineColor(kRed + 1); l->SetLineWidth(2); l->Draw();
                TLine* le = new TLine(cen, fR[i][0].edge - fR[i][0].eedge, cen, fR[i][0].edge + fR[i][0].eedge); le->SetLineColor(kRed + 1); le->SetLineWidth(2); le->Draw();
                Mark(cen, fR[i][0].edge, kRed + 1, 20);
                Mark(cen, fT[i][0].edge, kGray + 2, 24);
            }
            for (int k = 0; k < kNDataEdge; ++k) {
                const double cen = 0.5 * (kDataEdgeLo[k] + kDataEdgeHi[k]) + 0.002, e = scan::M4max((cen - 0.002) / kDataScalePiPi);
                TLine* le = new TLine(cen, kDataEdge[k] - kDataEdgeFm[k] * e, cen, kDataEdge[k] + kDataEdgeFp[k] * e); le->SetLineColor(4); le->SetLineWidth(2); le->Draw();
                Mark(cen, kDataEdge[k], 4, 21);
            }
            TLegend* lg = new TLegend(0.17, 0.14, 0.62, 0.40);
            lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.030);
            lg->AddEntry(pR, "sim reco: fitted edge (bar = slice width)", "p");
            lg->AddEntry(pT, "sim TRUE: fitted edge (true M_{#pi#pi} slices)", "p");
            lg->AddEntry(pD, "data: fitted edge (no cut; 4th: 0.39-0.40)", "p");
            lg->AddEntry(cT, "M4max(M_{#pi#pi}) (true formula)", "l");
            lg->AddEntry(cS, Form("%.4f M4max(M_{#pi#pi}/%.4f) (sim scales)", kSimScale, sPiPiSim), "l");
            lg->AddEntry(cD, Form("%.4f M4max(M_{#pi#pi}/%.4f) (data scales)", kDataScale, kDataScalePiPi), "l");
            lg->Draw();
            Txt(0.16, 0.955, "(a) bump edge against M_{#pi#pi} (bump Dalitz + conv., pre)", 0.040);
            // factors
            c->cd();
            MkPad("es2", 1 / 3.0, 0.0, 2 / 3.0, 0.92, false, 0.16, 0.03, 0.06, 0.13);
            TH1D* f2 = new TH1D("frEs2", "", 100, 0.29, 0.42);
            f2->SetMinimum(0.975); f2->SetMaximum(1.030);
            Axes(f2, "M_{#pi^{+}#pi^{-}} slice centre [GeV/c^{2}]", "edge factor k = edge / M4max(M_{#pi#pi}/s_{#pi#pi})", 0.050, 1.55);
            f2->Draw("axis");
            TBox* bs = new TBox(0.29, kSimScale - gPre.emu / scan::kMeta, 0.42, kSimScale + gPre.emu / scan::kMeta);
            bs->SetFillColorAlpha(kRed - 9, 0.5); bs->Draw();
            TBox* bd = new TBox(0.29, kDataScale - kDataScaleErr, 0.42, kDataScale + kDataScaleErr);
            bd->SetFillColorAlpha(kAzure - 9, 0.5); bd->Draw();
            HLine(0.29, 0.42, 1.0, kGray + 2, 3, 2);
            for (int i = 0; i < 4; ++i) {
                const double cen = 0.5 * (slLo[i] + slHi[i]);
                TLine* le = new TLine(cen, fR[i][0].k - std::fabs(fR[i][0].ek), cen, fR[i][0].k + std::fabs(fR[i][0].ek)); le->SetLineColor(kRed + 1); le->SetLineWidth(2); le->Draw();
                Mark(cen, fR[i][0].k, kRed + 1, 20);
                TLine* lt = new TLine(cen - 0.002, fT[i][0].k - std::fabs(fT[i][0].ek), cen - 0.002, fT[i][0].k + std::fabs(fT[i][0].ek)); lt->SetLineColor(kGray + 2); lt->SetLineWidth(2); lt->Draw();
                Mark(cen - 0.002, fT[i][0].k, kGray + 2, 24);
            }
            for (int k = 0; k < kNDataEdge; ++k) {
                const double cen = 0.5 * (kDataEdgeLo[k] + kDataEdgeHi[k]) + 0.002;
                TLine* le = new TLine(cen, kDataEdgeF[k] - kDataEdgeFm[k], cen, kDataEdgeF[k] + kDataEdgeFp[k]); le->SetLineColor(4); le->SetLineWidth(2); le->Draw();
                Mark(cen, kDataEdgeF[k], 4, 21);
            }
            TBox* ba = new TBox(0.29, aR[0].k - std::fabs(aR[0].ek), 0.42, aR[0].k + std::fabs(aR[0].ek));
            ba->SetFillStyle(3354); ba->SetFillColor(kRed + 1); ba->Draw();
            Txt(0.16, 0.95, "(b) edge factor k per slice", 0.042);
            Txt(0.18, 0.75, Form("hatched red: aligned fit, all slices: k = %.4f #pm %.4f", aR[0].k, std::fabs(aR[0].ek)), 0.032, kRed + 1);
            Txt(0.18, 0.87, Form("bands: peak scale, sim %.4f (red), data %.4f (blue)", kSimScale, kDataScale), 0.032);
            Txt(0.18, 0.83, Form("s_{#pi#pi}: sim %.4f, data %.4f, TRUE 1", sPiPiSim, kDataScalePiPi), 0.032);
            Txt(0.18, 0.79, "full red: sim reco; open grey: sim TRUE; blue square: data", 0.032);
            // (c) aligned spectrum
            c->cd();
            MkPad("es3", 2 / 3.0, 0.0, 1.0, 0.92, false, 0.16, 0.03, 0.06, 0.13);
            TH1D* hu = uR[0].Hist("huR", 1);
            TH1D* ht = uT[0].Hist("huT", 1);
            hu->Scale(1.0 / hu->GetBinWidth(1)); ht->Scale(1.0 / ht->GetBinWidth(1));
            StyleData(hu, kRed + 1, 20, 0.9); StyleLine(ht, kGray + 2, 2, 2);
            TH1D* f3 = Frame(hu, "frEs3", 0, std::max(MaxOf(hu, true), MaxOf(ht)) * 1.35);
            f3->GetXaxis()->SetLimits(0.90, 1.05);
            Axes(f3, "u = M_{4} / M4max(M_{#pi#pi}/s_{#pi#pi})", "d#sigma/du [nb]", 0.050, 1.55);
            f3->Draw("axis");
            ht->Draw("hist same");
            hu->Draw("E1 P same");
            const EdgeModel ea = EdgeModel::Aligned();
            const ERes ra = aR[0];
            TGraph* gfu = Curve([ea, ra](double x0, double bw) { const double p[4] = {ra.A, ra.B, ra.k, ra.s}; return ea(x0, bw, p); },
                                uR[0].bw, ra.lo, ra.hi, 1.0 / (kLsim * uR[0].bw), kRed + 1, 3);
            gfu->Draw("L same");
            VLine(aR[0].k, 0, f3->GetMaximum() * 0.75, kRed + 1, 1, 2);
            VLine(kDataScale, 0, f3->GetMaximum() * 0.75, 4, 2, 2);
            VLine(1.0, 0, f3->GetMaximum() * 0.75, kGray + 2, 3, 2);
            Txt(0.16, 0.955, "(c) all bump rows aligned row by row, 0.30 #leq M_{#pi#pi} < 0.41", 0.040);
            Txt(0.18, 0.87, Form("reco (red): k = %.4f #pm %.4f, s = %.4f, #chi^{2}/ndf %.0f/%d", aR[0].k, std::fabs(aR[0].ek), aR[0].s, aR[0].chi2, aR[0].ndf), 0.032, kRed + 1);
            Txt(0.18, 0.83, Form("TRUE (grey dashed): k = %.4f #pm %.4f", aT[0].k, std::fabs(aT[0].ek)), 0.032, kGray + 2);
            Txt(0.18, 0.79, Form("blue dashed: data scale %.4f", kDataScale), 0.032, 4);
            c->cd();
            Txt(0.005, 0.965, (kSimLine + "; data values from the data study (not recomputed)").c_str(), 0.0175);
            c->SaveAs(FigPath("spectra_eta_edge_summary.png").c_str());
            c->SaveAs(FigPath("spectra_eta_edge_summary.pdf").c_str());
        }
    }
    if (fo) { fo->Close(); delete fo; }
    printf("sim_spectra: total %.1f s\n", swAll.RealTime());
}
