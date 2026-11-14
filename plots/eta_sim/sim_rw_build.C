// ============================================================================
// sim_rw_build.C  --  physics reweighting of the simulated eta decays: build and validation
// ============================================================================
// HADES p+p, T = 4.53 GeV, pi+ pi- e+ e-. Model and formulas: sim_eta_physics.h.
//
// Part 1 (bit 1) generator and model, no simulation rows needed (about 30 s):
//   - reads the PLUTO generator sample of the eta set's eta -> pi+pi-e+e-
//     (/hdd3/.../generator/eta_pipiee_gen_kin_1000files.txt.gz: M(ee), M(pipi),
//     cos th*_e+, cos th*_pi+, M(eta); 1 181 801 decays from 1000 files)
//   - fits the factorised generator density G(q, M) = g(q) h(M) lambda^{1/2} by
//     maximum likelihood (q: log bins from 2 m_e with a threshold bin up to 1.553 MeV,
//     M: 2 MeV bins), tests the factorisation and the angular distributions (flat?
//     independent of the masses?), measures the generator angular densities
//   - eta -> pi+pi-pi0 phase-space Monte Carlo (TGenPhaseSpace, 4e6 decays): <f>_PS of
//     the KLOE Dalitz-plot function, the slides' 1 + aY, check of <X^2>_Y = Xmax^2/3;
//     eta -> pi+pi-gamma phase space: check of the conversion normalisation
//   - writes /hdd3/.../scan/rw_generator.root (g_q, h_m, consts and diagnostic
//     histograms), measures the generator-level mean weights (renormalisation to 1)
//   - target model T: PDG ratio check Gamma(pi pi e e)/Gamma(pi pi gamma), shape
//     against the QED reference of tables/eta_peak_mee_true.txt; closure of the
//     reweighted generator sample against T (V0-V3)
//   tables rw_weights_model.txt, rw_weights_generator.txt, rw_weights_closure.txt,
//   rw_weights_bump.txt; figures rw_weights_generator, rw_weights_model, rw_weights_closure, rw_weights_bump
// Part 2 (bit 2) true angles of class-5 rows (about 30 s): reads the true momenta of the
//   four tracks from scan_nt of output_eta_pippimepem.root, aligned entry by entry with
//   the compact tree eta_us (same order, sim_prepare.C keeps every row), checks the
//   alignment with the true masses, writes /hdd3/.../scan/rw_angles_eta.root (tree ang)
// Part 3 (bit 4) simulation (about 3 min, 8 GB): weights per class, cross sections per
//   component, selection and window before / after, peak classes of OA and M(ee) and the
//   OA > 2 efficiency against data, accepted direct-decay shapes, acceptance in true M(ee)
//   tables rw_weights_classes.txt, rw_weights_xsec.txt, rw_weights_peak.txt;
//   figures rw_weights_classes, rw_weights_xsec, rw_weights_peak, rw_weights_acceptance
//
// Run from /home/przygoda/HADES/FAT/FAT:
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/aclic_rw_build", kTRUE)' 'plots/eta_sim/sim_rw_build.C+(7)'
// L_data = 6080.5 nb^-1 UNCONFIRMED; simulation without CB; TH1D and double only.
// ============================================================================
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TF1.h"
#include "TGenPhaseSpace.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TPad.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TSystem.h"

#include "SimEngine.h"
#include "sim_eta_physics.h"

namespace srw {
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
void Save(TCanvas* c, const char* name, bool pdf) {
    c->SaveAs((std::string(kFigDir) + "/" + name + ".png").c_str());
    if (pdf) c->SaveAs((std::string(kFigDir) + "/" + name + ".pdf").c_str());
}
void Txt(double x, double y, const std::string& s, double size = 0.04, int col = 1, int align = 11) {
    TLatex t; t.SetNDC(); t.SetTextSize(size); t.SetTextColor(col); t.SetTextAlign(align); t.DrawLatex(x, y, s.c_str());
}
void Style() {
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
    gStyle->SetLabelSize(0.045, "XY"); gStyle->SetTitleSize(0.05, "XY");
    gStyle->SetTitleOffset(1.15, "Y"); gStyle->SetEndErrorSize(0);
}
TH1D* H1(const char* n, int nb, double lo, double hi) { TH1D* h = new TH1D(Nm(n).c_str(), "", nb, lo, hi); h->Sumw2(); h->SetDirectory(nullptr); return h; }
TH1D* H1v(const char* n, const std::vector<double>& e) { TH1D* h = new TH1D(Nm(n).c_str(), "", (int)e.size() - 1, e.data()); h->Sumw2(); h->SetDirectory(nullptr); return h; }
void Sty(TH1* h, int col, int lw = 2, int ls = 1, int ms = 0) {
    h->SetLineColor(col); h->SetLineWidth(lw); h->SetLineStyle(ls); h->SetMarkerColor(col);
    if (ms) { h->SetMarkerStyle(ms); h->SetMarkerSize(0.8); }
}
// divide contents by bin width
void PerWidth(TH1D* h, double scale = 1.) {
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        const double w = h->GetBinWidth(i);
        h->SetBinContent(i, h->GetBinContent(i) * scale / w); h->SetBinError(i, h->GetBinError(i) * scale / w);
    }
}
std::vector<double> LogEdges(double a, double b, int n) { std::vector<double> e; for (int i = 0; i <= n; ++i) e.push_back(a * std::pow(b / a, (double)i / n)); return e; }

// colours of the variants
const int kColGen = 1, kColV0 = 618 /*kMagenta+2*/, kColV1 = 807 /*kOrange+7*/, kColV2 = 601 /*kBlue+1*/, kColV3 = 418 /*kGreen+2*/,
          kColGS = 921 /*kGray+1*/, kColRef = 632 /*kRed*/;
const char* const kLdNote = "L_{data} = 6080.5 nb^{-1} UNCONFIRMED";

// ---------------------------------------------------------------- generator sample
struct Gen { std::vector<double> q, m, ce, cp, me; };
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
    printf("ReadGen: %zu decays in %.1f s\n", g.q.size(), sw.RealTime());
    return !g.q.empty();
}
// G binning: q log bins (40 log steps from 2 m_e to q_max, the first three merged into the threshold bin);
// M: 2 MeV bins from 2 m_pi, the last bin merged up to m_eta - 2 m_e
std::vector<double> QEdgesG() {
    const double a = 2 * kMe, b = QMax(); const int n = 40;
    std::vector<double> e; e.push_back(a);
    for (int k = 3; k <= n; ++k) e.push_back(a * std::pow(b / a, (double)k / n));
    return e;
}
std::vector<double> MEdgesG() {
    std::vector<double> e;
    for (int k = 0; k <= 132; ++k) e.push_back(2 * kMpi + 0.002 * k);
    e.push_back(kMetaGen - 2 * kMe);
    return e;
}
inline int Bin(const std::vector<double>& e, double x) { return BinOf(e, x); }

// cell integrals of lambda^{1/2} (q bin a, M bin b) inside q + M < m_eta
std::vector<double> CellLam(const std::vector<double>& qE, const std::vector<double>& mE) {
    const int nq = (int)qE.size() - 1, nm = (int)mE.size() - 1;
    std::vector<double> L(nq * nm, 0.), x, w;
    GL(16, x, w);
    for (int a = 0; a < nq; ++a)
        for (int i = 0; i < 16; ++i) {
            const double q = 0.5 * (qE[a] + qE[a + 1]) + 0.5 * (qE[a + 1] - qE[a]) * x[i], wq = 0.5 * (qE[a + 1] - qE[a]) * w[i];
            for (int b = 0; b < nm; ++b) {
                const double m0 = mE[b], m1 = std::min(mE[b + 1], kMetaGen - q);
                if (m1 <= m0) continue;
                double s = 0;
                for (int j = 0; j < 16; ++j) {
                    const double M = 0.5 * (m0 + m1) + 0.5 * (m1 - m0) * x[j];
                    const double l = Lam(kMetaGen * kMetaGen, M * M, q * q);
                    s += 0.5 * (m1 - m0) * w[j] * (l > 0 ? std::sqrt(l) : 0.);
                }
                L[a * nm + b] += wq * s;
            }
        }
    return L;
}
// method of moments for a(c) ~ 1 + b c^2: <c^2> = (2/3 + 2b/5) / (2 + 2b/3)
struct AngFit { double n = 0, b = 0, eb = 0, c2 = 0, ec2 = 0, c1 = 0, ec1 = 0; };
AngFit FitAng(const std::vector<double>& c, const std::vector<char>& use) {
    double n = 0, s1 = 0, s2 = 0, s4 = 0;
    for (size_t i = 0; i < c.size(); ++i) if (use[i]) { n += 1; s1 += c[i]; s2 += c[i] * c[i]; s4 += std::pow(c[i], 4); }
    AngFit f; f.n = n;
    if (n < 10) return f;
    f.c1 = s1 / n; f.ec1 = std::sqrt(std::max(0., s2 / n - f.c1 * f.c1) / n);
    f.c2 = s2 / n; f.ec2 = std::sqrt(std::max(0., s4 / n - f.c2 * f.c2) / n);
    auto bOf = [](double m2) { return (2 * m2 - 2. / 3.) / (0.4 - 2 * m2 / 3.); };
    f.b = bOf(f.c2);
    f.eb = 0.5 * std::fabs(bOf(f.c2 + f.ec2) - bOf(f.c2 - f.ec2));
    return f;
}

// expected fraction of T_v (dir variant dv) in [q0, q1] (q-marginal) and in [m0, m1] (M-marginal)
double TeeInt(double q0, double q1, int dv, int gs) {
    const int v = dv == kV3 ? 2 : dv - 1;
    const Model& m = Mdl();
    q0 = std::max(q0, 2 * kMe); q1 = std::min(q1, QMax());
    if (q1 <= q0) return 0.;
    const int n = 4000;
    double s = 0;
    // split at the threshold-bin edge (KWeff is a step there)
    const double cut = m.qE[1];
    auto seg = [&](double a, double b) {
        if (b <= a) return;
        const double smax = std::sqrt(b - a);
        for (int i = 0; i < n; ++i) {
            const double sv = smax * (i + 0.5) / n, q = a + sv * sv, dq = 2 * sv * smax / n;
            s += KWeff(q) * Interp(m.Iee[v][gs], m.lq0, m.lq1, std::log(q)) * dq;
        }
    };
    seg(q0, std::min(q1, cut)); seg(std::max(q0, cut), q1);
    return s / m.ZT[v][gs];
}
double TppInt(double m0, double m1, int dv, int gs) {
    const int v = dv == kV3 ? 2 : dv - 1;
    const Model& m = Mdl();
    const int n = 200; double s = 0;
    for (int i = 0; i < n; ++i) { const double x = m0 + (m1 - m0) * (i + 0.5) / n; s += Interp(m.Tpp[v][gs], m.m0, m.m1, x) * (m1 - m0) / n; }
    return s / m.ZT[v][gs];
}
// PDG-ratio integral: int dq KW(q) int dM R_v(q, M) / int dM R_v(0, M)
double PdgRatio(int v, int gs) {
    const Model& m = Mdl();
    const int n = 40000; const double a = 2 * kMe, smax = std::sqrt(QMax() - a);
    double num = 0;
    for (int i = 0; i < n; ++i) {
        const double sv = smax * (i + 0.5) / n, q = a + sv * sv, dq = 2 * sv * smax / n;
        num += KW(q) * Interp(m.Iee[v][gs], m.lq0, m.lq1, std::log(q)) * dq;
    }
    std::vector<double> x, w; GL(64, x, w);
    double den = 0;
    const double lo = 2 * kMpi, hi = kMetaGen;
    for (int k = 0; k < 64; ++k) { const double M = 0.5 * (lo + hi) + 0.5 * (hi - lo) * x[k]; den += 0.5 * (hi - lo) * w[k] * Rmass(0., M, v, gs); }
    return num / den;
}

// =====================================================================================================
// Part 1
// =====================================================================================================
void Part1() {
    TStopwatch sw;
    Style();
    Gen G;
    if (!ReadGen(G)) return;
    const size_t N = G.q.size();
    const std::vector<double> qE = QEdgesG(), mE = MEdgesG();
    const int nq = (int)qE.size() - 1, nm = (int)mE.size() - 1;
    // ---------------- counts
    std::vector<double> nQ(nq, 0.), nMb(nm, 0.), O(nq * nm, 0.);
    long long nBelow = 0, nZero = 0;
    double metaS = 0, metaS2 = 0;
    for (size_t i = 0; i < N; ++i) {
        if (G.q[i] < 2 * kMe) ++nBelow;
        if (G.q[i] == 0) ++nZero;
        const int a = Bin(qE, G.q[i]), b = Bin(mE, G.m[i]);
        nQ[a] += 1; nMb[b] += 1; O[a * nm + b] += 1;
        metaS += G.me[i]; metaS2 += G.me[i] * G.me[i];
    }
    const double metaMean = metaS / N, metaRms = std::sqrt(std::max(0., metaS2 / N - metaMean * metaMean));
    // ---------------- ML fit of g, h
    const std::vector<double> L = CellLam(qE, mE);
    std::vector<double> g(nq, 1.), h(nm, 1.);
    for (int b = 0; b < nm; ++b) h[b] = nMb[b] / (mE[b + 1] - mE[b]) + 1e-12;
    int it = 0; double change = 1;
    for (it = 0; it < 2000 && change > 1e-13; ++it) {
        change = 0;
        for (int a = 0; a < nq; ++a) {
            double den = 0; for (int b = 0; b < nm; ++b) den += h[b] * L[a * nm + b];
            const double gn = den > 0 ? nQ[a] / den : 0.;
            change = std::max(change, std::fabs(gn - g[a]) / std::max(gn, 1e-300)); g[a] = gn;
        }
        for (int b = 0; b < nm; ++b) {
            double den = 0; for (int a = 0; a < nq; ++a) den += g[a] * L[a * nm + b];
            h[b] = den > 0 ? nMb[b] / den : 0.;
        }
    }
    double Z = 0; for (int a = 0; a < nq; ++a) for (int b = 0; b < nm; ++b) Z += g[a] * h[b] * L[a * nm + b];
    for (int a = 0; a < nq; ++a) g[a] /= Z;
    // ---------------- factorisation test: observed against predicted cell counts (M bins merged by 4)
    double chiTot = 0; int ndfTot = 0;
    std::vector<double> chiQ(nq, 0.); std::vector<int> nQc(nq, 0);
    const int mg = 4;
    for (int a = 0; a < nq; ++a)
        for (int b0 = 0; b0 < nm; b0 += mg) {
            double o = 0, p = 0;
            for (int b = b0; b < std::min(nm, b0 + mg); ++b) { o += O[a * nm + b]; p += N * g[a] * h[b] * L[a * nm + b]; }
            if (p < 10) continue;
            const double c = (o - p) * (o - p) / p;
            chiTot += c; ndfTot += 1; chiQ[a] += c; nQc[a] += 1;
        }
    const int nPar = nq + nm - 1;
    // slices of q for the report (8 groups of q bins)
    // ---------------- angular distributions
    std::vector<char> useAll(N, 1), useHi(N, 0);
    for (size_t i = 0; i < N; ++i) useHi[i] = G.q[i] > 0.002;
    const AngFit fe = FitAng(G.ce, useHi), fe0 = FitAng(G.ce, useAll), fp = FitAng(G.cp, useAll);
    // slices: q (8 log groups above 2 MeV + threshold), M (4)
    const std::vector<double> qS = {0., 0.002, 0.004, 0.008, 0.016, 0.032, 0.064, 0.128, 0.27};
    const std::vector<double> mS = {0.27, 0.36, 0.40, 0.44, 0.55};
    std::vector<AngFit> feQ, fpQ, feM, fpM;
    for (size_t s = 0; s + 1 < qS.size(); ++s) {
        std::vector<char> u(N, 0); for (size_t i = 0; i < N; ++i) u[i] = G.q[i] >= qS[s] && G.q[i] < qS[s + 1];
        feQ.push_back(FitAng(G.ce, u)); fpQ.push_back(FitAng(G.cp, u));
    }
    for (size_t s = 0; s + 1 < mS.size(); ++s) {
        std::vector<char> u(N, 0); for (size_t i = 0; i < N; ++i) u[i] = G.m[i] >= mS[s] && G.m[i] < mS[s + 1] && G.q[i] > 0.002;
        feM.push_back(FitAng(G.ce, u));
        for (size_t i = 0; i < N; ++i) u[i] = G.m[i] >= mS[s] && G.m[i] < mS[s + 1];
        fpM.push_back(FitAng(G.cp, u));
    }
    // decision on the generator angular densities
    const bool piIsSin2 = std::fabs(fp.b + 1) < 3 * fp.eb + 0.003;
    const double bEuse = 0.;   // e+: flat above 2 MeV (decided below from fe)
    const double bPiUse = piIsSin2 ? -1. : fp.b;
    // angular histograms
    TH1D* hCe = H1("hCe", 40, -1, 1); TH1D* hCeLo = H1("hCeLo", 40, -1, 1); TH1D* hCp = H1("hCp", 40, -1, 1);
    for (size_t i = 0; i < N; ++i) { (G.q[i] > 0.002 ? hCe : hCeLo)->Fill(G.ce[i]); hCp->Fill(G.cp[i]); }
    // correlation of the two angles (independence)
    double sxy = 0, sx = 0, sy = 0, sxx = 0, syy = 0;
    for (size_t i = 0; i < N; ++i) { const double x = G.ce[i] * G.ce[i], y = G.cp[i] * G.cp[i]; sx += x; sy += y; sxy += x * y; sxx += x * x; syy += y * y; }
    const double corr = (sxy / N - sx / N * sy / N) / std::sqrt((sxx / N - sx * sx / N / N) * (syy / N - sy * sy / N / N));

    // ---------------- eta -> pi+pi-pi0 phase-space Monte Carlo
    TRandom3 rnd(4530); gRandom = &rnd;
    double sW = 0, sF = 0, sFbar = 0, sS = 0, sFx = 0;
    const int nYB = 10; std::vector<double> yX2(nYB, 0.), yXf(nYB, 0.), yW(nYB, 0.);
    double xOverMax = 0;
    TH1D* hPSm = H1("hPSm", 70, 0.27, 0.42);
    {
        TGenPhaseSpace ps; TLorentzVector eta(0, 0, 0, kMetaHG); double ms[3] = {kMpi, kMpi, kMpi0};
        ps.SetDecay(eta, 3, ms);
        const double Q = DalQ(kMetaHG);
        for (int i = 0; i < 4000000; ++i) {
            const double wt = ps.Generate();
            const TLorentzVector& a = *ps.GetDecay(0); const TLorentzVector& b = *ps.GetDecay(1); const TLorentzVector& c = *ps.GetDecay(2);
            const double Tp = a.E() - kMpi, Tm = b.E() - kMpi, T0 = c.E() - kMpi0;
            const double X = std::sqrt(3.) * (Tp - Tm) / Q, Y = 3 * T0 / Q - 1;
            const double Mpp = (a + b).M();
            sW += wt; sF += wt * DalF(X, Y, kBumpKLOE); sFbar += wt * DalFbar(Mpp, kBumpKLOE); sS += wt * DalF(X, Y, kBumpSlides);
            sFx += wt * (1 + kDalA * Y + kDalB * Y * Y + kDalF * Y * Y * Y);
            const double xm = DalXmax(Mpp);
            if (xm > 0) xOverMax = std::max(xOverMax, std::fabs(X) / xm);
            const int k = std::min(nYB - 1, std::max(0, (int)((Y + 1) / 2. * nYB)));
            yX2[k] += wt * X * X; yXf[k] += wt * xm * xm / 3.; yW[k] += wt;
            hPSm->Fill(Mpp, wt);
        }
    }
    const double normKLOE = sF / sW, normKLOEbar = sFbar / sW, normSlides = sS / sW, normNoX = sFx / sW;
    // eta -> pi+pi-gamma phase space: check of normC
    double sWg = 0, sCg[2] = {0, 0};
    {
        TGenPhaseSpace ps; TLorentzVector eta(0, 0, 0, kMetaHG); double ms[3] = {kMpi, kMpi, 0.};
        ps.SetDecay(eta, 3, ms);
        for (int i = 0; i < 2000000; ++i) {
            const double wt = ps.Generate();
            const double Mpp = (*ps.GetDecay(0) + *ps.GetDecay(1)).M();
            const double p = PStar(Mpp, kMpi), k = K3(kMetaHG, Mpp, 0.);
            sWg += wt; for (int gs = 0; gs < 2; ++gs) sCg[gs] += wt * p * p * k * k * FV2(Mpp * Mpp, gs);
        }
    }

    // ---------------- write rw_generator.root (renormalisation 1 at first)
    auto writeFile = [&](const double ren5[4][2], const double renEE[3][2], const double renPP[3][2]) {
        TFile f(kRwGenFile, "RECREATE");
        TH1D hg("g_q", "generator density factor g(q) [per GeV, product normalisation]", nq, qE.data());
        for (int a = 0; a < nq; ++a) hg.SetBinContent(a + 1, g[a]);
        TH1D hh("h_m", "generator density factor h(M(pipi))", nm, mE.data());
        for (int b = 0; b < nm; ++b) hh.SetBinContent(b + 1, h[b]);
        TH1D hc("consts", "1-8 ren5, 9-14 renEE, 15-20 renPP, 21 bE, 22 bPi, 23 <f>_PS KLOE, 24 <f>_PS slides", 30, 0.5, 30.5);
        for (int v = 0; v < 4; ++v) for (int fv = 0; fv < 2; ++fv) hc.SetBinContent(1 + 2 * v + fv, ren5[v][fv]);
        for (int v = 0; v < 3; ++v) for (int fv = 0; fv < 2; ++fv) { hc.SetBinContent(9 + 2 * v + fv, renEE[v][fv]); hc.SetBinContent(15 + 2 * v + fv, renPP[v][fv]); }
        hc.SetBinContent(21, bEuse); hc.SetBinContent(22, bPiUse); hc.SetBinContent(23, normKLOEbar); hc.SetBinContent(24, normSlides);
        TH1D hn("gen_q_counts", "generator decays per q bin of G", nq, qE.data());
        for (int a = 0; a < nq; ++a) hn.SetBinContent(a + 1, nQ[a]);
        TH1D hm("gen_m_counts", "generator decays per M(pipi) bin of G", nm, mE.data());
        for (int b = 0; b < nm; ++b) hm.SetBinContent(b + 1, nMb[b]);
        TH2D h2("gen_qm_counts", "generator decays, q bin x M(pipi) bin", nq, qE.data(), nm, mE.data());
        for (int a = 0; a < nq; ++a) for (int b = 0; b < nm; ++b) h2.SetBinContent(a + 1, b + 1, O[a * nm + b]);
        TH2D h2p("model_qm_counts", "factorised model, expected decays per cell", nq, qE.data(), nm, mE.data());
        for (int a = 0; a < nq; ++a) for (int b = 0; b < nm; ++b) h2p.SetBinContent(a + 1, b + 1, N * g[a] * h[b] * L[a * nm + b]);
        TH1D* c1 = (TH1D*)hCe->Clone("gen_cos_e_q_gt_2MeV"); TH1D* c2 = (TH1D*)hCeLo->Clone("gen_cos_e_q_lt_2MeV"); TH1D* c3 = (TH1D*)hCp->Clone("gen_cos_pi");
        hg.Write(); hh.Write(); hc.Write(); hn.Write(); hm.Write(); h2.Write(); h2p.Write(); c1->Write(); c2->Write(); c3->Write();
        f.Close();
    };
    double one5[4][2], one3[3][2];
    for (int v = 0; v < 4; ++v) for (int fv = 0; fv < 2; ++fv) one5[v][fv] = 1.;
    for (int v = 0; v < 3; ++v) for (int fv = 0; fv < 2; ++fv) one3[v][fv] = 1.;
    writeFile(one5, one3, one3);
    if (!Init(kRwGenFile)) return;
    // raw generator-level means
    double raw5[4][2], rawEE[3][2], rawPP[3][2];
    for (int gs = 0; gs < 2; ++gs) {
        double s5[4] = {0, 0, 0, 0}, se[3] = {0, 0, 0}, sp[3] = {0, 0, 0};
        for (size_t i = 0; i < N; ++i) {
            for (int dv = kV0; dv <= kV3; ++dv) s5[dv - 1] += W5(G.q[i], G.m[i], dv, gs, G.ce[i], G.cp[i], true);
            for (int dv = kV0; dv <= kV2; ++dv) { se[dv - 1] += Wee(G.q[i], dv, gs); sp[dv - 1] += Wpp(G.m[i], dv, gs); }
        }
        for (int v = 0; v < 4; ++v) raw5[v][gs] = s5[v] / N;
        for (int v = 0; v < 3; ++v) { rawEE[v][gs] = se[v] / N; rawPP[v][gs] = sp[v] / N; }
    }
    double r5[4][2], rE[3][2], rP[3][2];
    for (int v = 0; v < 4; ++v) for (int gs = 0; gs < 2; ++gs) r5[v][gs] = 1. / raw5[v][gs];
    // V3 is stored relative to the raw V3 mean (W5 applies ren5[3] for V3)
    for (int v = 0; v < 3; ++v) for (int gs = 0; gs < 2; ++gs) { rE[v][gs] = 1. / rawEE[v][gs]; rP[v][gs] = 1. / rawPP[v][gs]; }
    writeFile(r5, rE, rP);
    if (!Init(kRwGenFile)) return;

    // ================= tables
    OpenTab("rw_weights_generator.txt");
    P("# rw_weights_generator: generator density G of the eta set's eta -> pi+pi-e+e- (PLUTO, makeStaticData()->AddDecay, eta_physics plugin off)\n");
    P("# input %s\n", kGenText);
    P("decays %zu (1000 files); M(eta) column %.7f +- %.7f GeV (rms)\n", N, metaMean, metaRms);
    P("M(ee) below 2 m_e = 1.0220 MeV: %lld (%.2f %%), exactly 0: %lld -> limited precision of the stored momenta (UNVERIFIED cause);\n", nBelow, 100. * nBelow / N, nZero);
    P("   these values are counted in the threshold bin [2 m_e, %.3f MeV] of G; above about 2 MeV the effect is negligible\n", qE[1] * 1e3);
    double f152 = 0, f20 = 0; for (size_t i = 0; i < N; ++i) { if (G.q[i] < 0.00152) f152 += 1; if (G.q[i] > 0.020) f20 += 1; }
    P("fraction M(ee) < 1.52 MeV (incl. below 2 m_e): %.4f; M(ee) > 20 MeV: %.4f\n\n", f152 / N, f20 / N);
    P("## factorised model G(q, M) = g(q) h(M) lambda^{1/2}(m_eta^2, M^2, q^2), m_eta = %.3f GeV, maximum likelihood (%d iterations, last relative change %.1e)\n", kMetaGen, it, change);
    P("binning: q %d bins (threshold bin %.4f-%.4f MeV, then log steps x %.4f up to %.2f MeV); M(pipi) %d bins (2 MeV, last %.4f-%.4f GeV)\n", nq,
      qE[0] * 1e3, qE[1] * 1e3, qE[2] / qE[1], QMax() * 1e3, nm, mE[nm - 1], mE[nm]);
    const int ndfApprox = std::max(1, ndfTot - (nq - 1) - ((nm + mg - 1) / mg - 1));
    P("goodness of fit (cells q bin x 8 MeV in M(pipi) with >= 10 expected): chi2 = %.1f for %d cells, approximate ndf %d (p = %.3g); per q bin chi2/cells:\n",
      chiTot, ndfTot, ndfApprox, TMath::Prob(chiTot, ndfApprox));
    for (int a = 0; a < nq; ++a) P("  q %7.3f-%7.3f MeV  decays %7.0f  chi2 %6.1f / %d cells\n", qE[a] * 1e3, qE[a + 1] * 1e3, nQ[a], chiQ[a], nQc[a]);
    P("# (the fit has %d free shape parameters on the fine bins; approximate ndf = cells - (%d - 1) q - (%d - 1) merged-M parameters)\n", nPar, nq, (nm + mg - 1) / mg);
    P("\n## generator angular distributions (densities in cos th*; fit a(c) ~ 1 + b c^2 by moments; flat: b = 0, <c^2> = 1/3; sin^2: b = -1, <c^2> = 0.2)\n");
    P("cos th*_e+ (e+e- rest frame, w.r.t. the e+e- direction in the eta frame), M(ee) > 2 MeV: N %.0f, <c> %.5f +- %.5f, <c^2> %.5f +- %.5f, b = %.4f +- %.4f\n", fe.n, fe.c1, fe.ec1, fe.c2, fe.ec2, fe.b, fe.eb);
    P("cos th*_e+, all decays (incl. threshold): <c^2> %.5f +- %.5f, b = %.4f +- %.4f\n", fe0.c2, fe0.ec2, fe0.b, fe0.eb);
    P("cos th*_pi+ (pi+pi- rest frame, w.r.t. the pi+pi- direction): N %.0f, <c> %.5f +- %.5f, <c^2> %.5f +- %.5f, b = %.4f +- %.4f  -> %s\n", fp.n, fp.c1, fp.ec1, fp.c2,
      fp.ec2, fp.b, fp.eb, piIsSin2 ? "consistent with sin^2 (P wave) already in the generator" : "NOT sin^2");
    P("slices of M(ee) [GeV]:            b(e+)               b(pi+)          N\n");
    for (size_t s = 0; s < feQ.size(); ++s)
        P("  %6.3f-%6.3f            %7.4f +- %.4f    %7.4f +- %.4f   %8.0f\n", qS[s], qS[s + 1], feQ[s].b, feQ[s].eb, fpQ[s].b, fpQ[s].eb, feQ[s].n);
    P("slices of M(pipi) [GeV] (e+: M(ee) > 2 MeV):\n");
    for (size_t s = 0; s < feM.size(); ++s)
        P("  %6.3f-%6.3f            %7.4f +- %.4f    %7.4f +- %.4f   %8.0f\n", mS[s], mS[s + 1], feM[s].b, feM[s].eb, fpM[s].b, fpM[s].eb, fpM[s].n);
    P("correlation of cos^2 th*_e and cos^2 th*_pi: %.5f (N = %zu, 1/sqrt(N) = %.5f)\n", corr, N, 1 / std::sqrt((double)N));
    P("chi2 of the histograms (40 bins): e+ (M(ee) > 2 MeV) against flat %.1f / 39; pi+ against sin^2 %.1f / 39\n",
      [&]() { double c = 0; const double e = hCe->Integral() / 40; for (int i = 1; i <= 40; ++i) c += std::pow(hCe->GetBinContent(i) - e, 2) / e; return c; }(),
      [&]() { double c = 0; const double n = hCp->Integral(); for (int i = 1; i <= 40; ++i) { const double a = hCp->GetBinLowEdge(i), b = a + 0.05;
              const double e = n * 0.75 * ((b - b * b * b / 3) - (a - a * a * a / 3)); c += std::pow(hCp->GetBinContent(i) - e, 2) / e; } return c; }());
    P("decision: generator e+ angle FLAT (b_e = 0; near threshold the stored precision distorts it, where the QED factor is ~1 anyway); generator pi+ angle %s\n",
      piIsSin2 ? "sin^2 (the V3 pion factor is exactly 1)" : Form("1 + %.4f c^2", fp.b));
    CloseTab();

    OpenTab("rw_weights_bump.txt");
    P("# rw_weights_bump: eta -> pi+pi-pi0 Dalitz-plot weight of the bump classes (2, 3) and the eta -> pi+pi-gamma weight of the conversion class (4)\n");
    P("# f(X,Y) = 1 + aY + bY^2 + dX^2 + fY^3 + gX^2Y, KLOE-2 2016 (A. Anastasi et al., JHEP 05 (2016) 019, arXiv:1601.06985): a %.3f b %.3f d %.3f f %.3f g %.3f\n",
      kDalA, kDalB, kDalD, kDalF, kDalG);
    P("# variant 'slides' (X17 analysis, K. Proscinski): 1 + aY, a = %.3f\n", kDalA);
    P("# HGeant (GEANT3) eta -> pi+pi-pi0 ASSUMED phase space (flat Dalitz plot); m_eta = %.6f, m_pi = %.6f, m_pi0 = %.6f GeV, Q = %.5f GeV\n", kMetaHG, kMpi, kMpi0, DalQ(kMetaHG));
    P("phase-space MC (TGenPhaseSpace, 4e6 weighted decays, seed 4530):\n");
    P("  <f>_PS full KLOE f(X,Y)                        = %.5f\n", normKLOE);
    P("  <f>_PS with X^2 -> Xmax(Y)^2/3 (used)          = %.5f   (difference to full: %.2e, MC statistics)\n", normKLOEbar, normKLOEbar - normKLOE);
    P("  <1 + aY + bY^2 + fY^3>_PS (X terms dropped, the earlier sim_cuts.C variant) = %.5f\n", normNoX);
    P("  <1 + aY>_PS (slides)                            = %.5f\n", normSlides);
    P("  largest |X| / Xmax(Y) over the MC events: %.6f (Dalitz boundary formula Xmax = 2 sqrt(3) p* q0* / (m_eta Q))\n", xOverMax);
    P("  <X^2> per Y bin, MC against Xmax^2/3:\n");
    for (int k = 0; k < nYB; ++k) if (yW[k] > 0) P("    Y %5.2f-%5.2f   MC %.5f   formula %.5f\n", -1 + 2. * k / nYB, -1 + 2. * (k + 1) / nYB, yX2[k] / yW[k], yXf[k] / yW[k]);
    P("weight wB(M(pipi)_true) = f_bar(Y) / <f>_PS (KLOE) and (1 + aY) / <1 + aY>_PS (slides):\n  M(pipi)    Y      Xmax    wB KLOE   wB slides   p3 (pi+pi-pi0 share of pion class 2)   wC (GP)   wC (GS)\n");
    for (double x = 0.285; x < 0.545; x += 0.01)
        P("  %.3f  %6.3f  %6.3f   %7.4f    %7.4f     %6.4f    %7.4f   %7.4f\n", x, x < kMetaHG - kMpi0 ? DalY(x) : -1., x < kMetaHG - kMpi0 ? DalXmax(x) : 0.,
          x < kMetaHG - kMpi0 ? WBump(x, kBumpKLOE) : 0., x < kMetaHG - kMpi0 ? WBump(x, kBumpSlides) : 0., P3(x), WConv(x, 0), WConv(x, 1));
    P("conversion: wC = |F_V|^2 p*^2 k^2 / <|F_V|^2 p*^2 k^2>_PS (HGeant eta -> pi+pi-gamma ASSUMED phase space, dPhi/dM ~ p* k)\n");
    P("  normalisation numeric (GP / GS): %.6e / %.6e; TGenPhaseSpace 2e6: %.6e / %.6e\n", Mdl().normC[0], Mdl().normC[1], sCg[0] / sWg, sCg[1] / sWg);
    P("pion class 2 mixture: p3(M) = B3 rho3 / (B3 rho3 + Bg rhog), HGeant B3 = %.4f, Bg = %.4f, rho = phase-space M(pipi) densities\n", kHGppp0, kHGppg);
    CloseTab();

    // ---------------- model table
    OpenTab("rw_weights_model.txt");
    P("# rw_weights_model: target model T of eta -> pi+pi-e+e- (formulas and sources in plots/eta_sim/sim_eta_physics.h)\n");
    P("# T(q, M) = KW(q) x p*^3 lambda^{3/2}(m_eta^2, M^2, q^2) x |F_V(M^2)|^2 x |F_VMD(q^2)|^2 [x angular factors]\n");
    P("#  KW(q) = (2 alpha/(3 pi q)) (1 + 2 m_e^2/q^2) sqrt(1 - 4 m_e^2/q^2)                 Kroll-Wada, Phys. Rev. 98 (1955) 1355\n");
    P("#  p*^3 lambda^{3/2}: P-wave pi pi gamma* rate per dM(pipi) (|A|^2 ~ p*^2 lambda sin^2 th_pi, phase space ~ lambda^{1/2} p* dM^2 / M)\n");
    P("#      = (beta_pi^3 M^3 / 8) lambda^{3/2}; q -> 0: p*^3 (2 m_eta k)^3, the eta -> pi+pi-gamma rate (Stollenwerk et al., PLB 707 (2012) 184)\n");
    P("#  F_V: Guerrero-Pich, PLB 412 (1997) 382, hep-ph/9707347 (pi + K loops, f_pi = 92.4 MeV); variant Gounaris-Sakurai, PRL 21 (1968) 244\n");
    P("#       m_rho = %.4f GeV, Gamma_rho = %.4f GeV (as in the X17 slides)\n", kMrho, kGrho);
    P("#  F_VMD = m_V^2 / (m_V^2 - q^2 - i m_V Gamma(q^2)), m_V = m_rho; Gamma(q^2) = 0 below 4 m_pi^2 (whole q range)\n");
    P("#  angles (V3): e+ (1 + cos^2 + 4 m_e^2/q^2 sin^2) [transverse gamma*], pi+ sin^2 [P wave]: the phi-average of the full density\n");
    P("#      sin^2 th_pi (1 - beta^2 sin^2 th_e sin^2 phi) of the same amplitude; the phi term (V4) is checked in rw_weights_planes.txt (below 0.5 %%)\n");
    P("#  V0 = KW x P-wave phase space; V1 = V0 x |F_V|^2; V2 = V1 x |F_VMD|^2; V3 = V2 x angular factors (MAIN since the review of 2026-11-12)\n");
    P("#  the M(ee) and M(pipi) marginals of V3 equal those of V2 (the angular factors are normalised at each q): rows 'V2' below hold for V3\n");
    P("# m_eta = %.3f GeV (PLUTO value of the generator), m_e = %.8f, m_pi = %.8f GeV, alpha = 1/137.036\n\n", kMetaGen, kMe, kMpi);
    P("## validation 1: Gamma(pi+pi-e+e-) / Gamma(pi+pi-gamma) = int dq KW int dM R(q,M) / int dM R(0,M), same F_V in both (PDG: 2.68e-4 / 4.22e-2 = %.3e +- %.2e)\n",
      kBRppee / kBRppg, kBRppee / kBRppg * std::sqrt(std::pow(0.11 / 2.68, 2) + std::pow(0.08 / 4.22, 2)));
    const double pr[5] = {PdgRatio(0, 0), PdgRatio(1, 0), PdgRatio(2, 0), PdgRatio(1, 1), PdgRatio(2, 1)};
    P("  V0 (no form factors)          %.4e   (ratio to PDG %.3f)\n", pr[0], pr[0] / (kBRppee / kBRppg));
    P("  V1 (F_V Guerrero-Pich)        %.4e   (%.3f)\n", pr[1], pr[1] / (kBRppee / kBRppg));
    P("  V2 = V3 (F_V GP + F_VMD, main) %.4e  (%.3f)   -> BR(eta -> pi+pi-e+e-) = %.3e x BR(pi+pi-gamma) 4.22e-2 = %.3e\n", pr[2], pr[2] / (kBRppee / kBRppg), pr[2], pr[2] * kBRppg);
    P("  V1 (F_V Gounaris-Sakurai)     %.4e   (%.3f)\n", pr[3], pr[3] / (kBRppee / kBRppg));
    P("  V2 (F_V GS + F_VMD)           %.4e   (%.3f)\n", pr[4], pr[4] / (kBRppee / kBRppg));
    P("  (the weights keep the generated number of decays, i.e. BR = 2.68e-4 of the production; this ratio only validates the model)\n");
    P("\n## |F_V(M^2)|^2 and |F_VMD(q^2)|^2\n  M [GeV]   |F_V|^2 GP   |F_V|^2 GS   GS/GP      q [GeV]  |F_VMD|^2\n");
    for (int i = 0; i < 9; ++i) {
        const double x = 0.28 + 0.033 * i, qq = 0.03 * i;
        P("  %.3f      %.4f       %.4f     %.4f     %.3f    %.4f\n", x, FV2(x * x, 0), FV2(x * x, 1), FV2(x * x, 1) / FV2(x * x, 0), qq, FVMD2(qq * qq));
    }
    // shape against the QED reference of eta_peak_mee_true.txt
    P("\n## validation 2: M(ee) shape against the QED reference of tables/eta_peak_mee_true.txt (model column; KW x P-wave pi pi gamma*, its M(pipi)\n");
    P("#  weight p*^3/M per dM instead of p*^3, no form factors): fractions per log bin of the true M(ee)\n");
    std::vector<double> refE = LogEdges(0.001022, 0.2685, 14), refM;
    {
        std::ifstream in(std::string(kTabDir) + "/eta_peak_mee_true.txt");
        std::string line;
        while (std::getline(in, line)) {
            double a, b, c, d, e;
            if (line.empty() || line[0] == '#') { if (refM.size() == 14) break; continue; }
            if (sscanf(line.c_str(), " %lf - %lf %lf %lf %lf", &a, &b, &c, &d, &e) == 5 && refM.size() < 14) refM.push_back(c);
        }
    }
    P("  M(ee) [MeV]          ref.model   V0       V0/ref   V1       V2       V2/ref    generator\n");
    double fg[14];
    for (int k = 0; k < 14; ++k) { fg[k] = 0; }
    for (size_t i = 0; i < N; ++i) { const double qq = std::max(G.q[i], 0.001022); for (int k = 0; k < 14; ++k) if (qq >= refE[k] && (qq < refE[k + 1] || k == 13)) { fg[k] += 1; break; } }
    for (int k = 0; k < 14; ++k) {
        const double v0 = TeeInt(refE[k], refE[k + 1], kV0, 0), v1 = TeeInt(refE[k], refE[k + 1], kV1, 0), v2 = TeeInt(refE[k], refE[k + 1], kV2, 0);
        const double rm = k < (int)refM.size() ? refM[k] : NAN;
        P("  %7.2f - %-7.2f     %.4f    %.4f   %.3f    %.4f   %.4f   %.3f     %.4f\n", refE[k] * 1e3, refE[k + 1] * 1e3, rm, v0, v0 / rm, v1, v2, v2 / rm, fg[k] / N);
    }
    double a20[3] = {TeeInt(0.020, 1, kV0, 0), TeeInt(0.020, 1, kV1, 0), TeeInt(0.020, 1, kV2, 0)};
    P("  fraction above 20 MeV: V0 %.4f, V1 %.4f, V2 %.4f, generator %.4f (QED reference quoted: 0.28)\n", a20[0], a20[1], a20[2], f20 / N);
    P("  fraction 1.02-1.52 MeV: V2 %.4f, generator %.4f\n", TeeInt(0, 0.00152, kV2, 0), f152 / N);
    P("\n## generator-level normalisation (mean weight over the %zu generated decays before renormalisation; T and G each normalised to 1)\n", N);
    const char* vn[4] = {"V0", "V1", "V2", "V3"};
    for (int gs = 0; gs < 2; ++gs) {
        P("  F_V %s:", gs ? "GS" : "GP");
        for (int v = 0; v < 4; ++v) P("  w5 %s %.5f", vn[v], raw5[v][gs]);
        for (int v = 0; v < 3; ++v) P("  wee %s %.5f", vn[v], rawEE[v][gs]);
        for (int v = 0; v < 3; ++v) P("  wpp %s %.5f", vn[v], rawPP[v][gs]);
        P("\n");
    }
    P("  -> the weights are divided by these means, so the generator-level mean is exactly 1 (number of decays and BR unchanged)\n");
    P("  T normalisations ZT: V0 %.5e, V1 GP %.5e, V2 GP %.5e, V1 GS %.5e, V2 GS %.5e; KW bin average in the threshold bin %.4f /GeV\n", Mdl().ZT[0][0],
      Mdl().ZT[1][0], Mdl().ZT[2][0], Mdl().ZT[1][1], Mdl().ZT[2][1], Mdl().KWbar0);
    CloseTab();

    // ---------------- closure: reweighted generator against T
    OpenTab("rw_weights_closure.txt");
    P("# rw_weights_closure: generator sample (%zu decays) reweighted with w5 against the target T (fractions per bin; errors from sum w^2)\n", N);
    std::vector<double> cqE = LogEdges(Mdl().qE[1], QMax(), 29);
    cqE.insert(cqE.begin(), 2 * kMe);   // first closure bin = threshold bin of G
    std::vector<double> cmE; for (int k = 0; k <= 27; ++k) cmE.push_back(2 * kMpi + (kMetaGen - 2 * kMe - 2 * kMpi) * k / 27.);
    struct Clo { TH1D *q, *m, *ce, *cp, *ceHi; double sw = 0, sw2 = 0, swHi = 0, wmin = 1e30, wmax = 0, chiQ = 0, chiM = 0, chiE = 0, chiP = 0, chiEHi = 0; int nQ = 0, nM = 0; };
    const int dvs[5] = {0, kV0, kV1, kV2, kV3};
    std::vector<Clo> clo(5);
    for (int k = 0; k < 5; ++k) {
        clo[k].q = H1v("cq", cqE); clo[k].m = H1v("cm", cmE); clo[k].ce = H1("cce", 20, -1, 1); clo[k].cp = H1("ccp", 20, -1, 1);
        clo[k].ceHi = H1("cceh", 20, -1, 1);
    }
    for (size_t i = 0; i < N; ++i) {
        const double qq = std::max(G.q[i], 2 * kMe * 1.0000001);
        for (int k = 0; k < 5; ++k) {
            const double w = k == 0 ? 1. : W5(G.q[i], G.m[i], dvs[k], 0, G.ce[i], G.cp[i], true);
            Clo& c = clo[k];
            c.q->Fill(qq, w); c.m->Fill(G.m[i], w); c.ce->Fill(G.ce[i], w); c.cp->Fill(G.cp[i], w);
            if (G.q[i] > 0.002) { c.ceHi->Fill(G.ce[i], w); c.swHi += w; }
            c.sw += w; c.sw2 += w * w; c.wmin = std::min(c.wmin, w); c.wmax = std::max(c.wmax, w);
        }
    }
    P("variant  mean w    rms w     min       max        N_eff/N  | chi2/bins: M(ee) (threshold bin + 29 log bins)  M(pipi) (27)  cos e+ (20)  cos pi+ (20) | cos e+, M(ee) > 2 MeV\n");
    std::vector<std::vector<double>> expQ(5), expM(5), expE(5), expP(5), expEHi(5);
    for (int k = 0; k < 5; ++k) {
        Clo& c = clo[k];
        const double mean = c.sw / N, rms = std::sqrt(std::max(0., c.sw2 / N - mean * mean));
        // expected
        for (int b = 0; b < (int)cqE.size() - 1; ++b) expQ[k].push_back(k == 0 ? NAN : TeeInt(cqE[b], cqE[b + 1], dvs[k], 0));
        for (int b = 0; b < (int)cmE.size() - 1; ++b) expM[k].push_back(k == 0 ? NAN : TppInt(cmE[b], cmE[b + 1], dvs[k], 0));
        for (int b = 0; b < 20; ++b) {
            const double a = -1 + 0.1 * b, bb = a + 0.1;
            expP[k].push_back(0.75 * ((bb - bb * bb * bb / 3) - (a - a * a * a / 3)));
            if (dvs[k] == kV3) {   // e+ angle: int dq T_ee(q) int A_e dc (all q, and q > 2 MeV normalised to that range)
                double s = 0, sh = 0, nh = 0;
                for (int j = 0; j < 400; ++j) {
                    const double q0 = std::exp(std::log(2 * kMe) + (std::log(QMax()) - std::log(2 * kMe)) * j / 400.);
                    const double q1 = std::exp(std::log(2 * kMe) + (std::log(QMax()) - std::log(2 * kMe)) * (j + 1) / 400.);
                    const double qm = std::sqrt(q0 * q1), r = kMe * kMe / (qm * qm);
                    auto I = [&](double c) { return 0.375 * ((c + c * c * c / 3) + 4 * r * (c - c * c * c / 3)) / (1 + 2 * r); };
                    const double t = TeeInt(q0, q1, kV2, 0);
                    s += t * (I(bb) - I(a));
                    const double th = TeeInt(std::max(q0, 0.002), std::max(q1, 0.002), kV2, 0);
                    sh += th * (I(bb) - I(a)); nh += th;
                }
                expE[k].push_back(s); expEHi[k].push_back(sh / nh);
            } else { expE[k].push_back(0.05); expEHi[k].push_back(0.05); }
        }
        auto chi = [&](TH1D* h, const std::vector<double>& e, double sw, int& n) {
            double cs = 0; n = 0;
            for (int b = 1; b <= h->GetNbinsX(); ++b) {
                const double o = h->GetBinContent(b) / sw, oe = h->GetBinError(b) / sw, x = e[b - 1];
                if (!(x == x) || oe <= 0) continue;
                cs += (o - x) * (o - x) / (oe * oe); ++n;
            }
            return cs;
        };
        int n1 = 0, n2 = 0, n3 = 0, n4 = 0;
        if (k > 0) { c.chiQ = chi(c.q, expQ[k], c.sw, n1); c.chiM = chi(c.m, expM[k], c.sw, n2); }
        c.chiE = chi(c.ce, expE[k], c.sw, n3); c.chiP = chi(c.cp, expP[k], c.sw, n4);
        int n5 = 0; c.chiEHi = chi(c.ceHi, expEHi[k], c.swHi, n5);
        P("%-7s  %.5f  %.4f  %.3e  %.3e  %.4f   | %8.1f / %d   %8.1f / %d   %8.1f / %d   %8.1f / %d | %8.1f / 20\n", k == 0 ? "gen" : vn[k - 1], mean, rms, c.wmin, c.wmax,
          c.sw * c.sw / c.sw2 / N, c.chiQ, n1, c.chiM, n2, c.chiE, n3, c.chiP, n4, c.chiEHi);
    }
    P("# 'gen' = unweighted generator; its chi2 columns compare it with flat (e+) and sin^2 (pi+) angular densities\n");
    P("# expected e+ angle for V0-V2: flat (no angular factor); for V3: int dq T_ee(q) A_e(c; q); last column: decays with M(ee) > 2 MeV only (below, the\n");
    P("#   stored precision distorts the generator e+ angle: b = 0.26 there, see rw_weights_generator.txt)\n");
    P("\nper bin, M(ee) [MeV]: generator fraction, V2 reweighted, T_V2, (V2 - T)/err\n");
    for (int b = 0; b < (int)cqE.size() - 1; ++b) {
        const double g0 = clo[0].q->GetBinContent(b + 1) / clo[0].sw, o = clo[3].q->GetBinContent(b + 1) / clo[3].sw, oe = clo[3].q->GetBinError(b + 1) / clo[3].sw;
        P("  %8.3f-%-8.3f  %.5f  %.5f +- %.5f  %.5f  %6.2f\n", cqE[b] * 1e3, cqE[b + 1] * 1e3, g0, o, oe, expQ[3][b], oe > 0 ? (o - expQ[3][b]) / oe : 0.);
    }
    CloseTab();

    // ================= figures
    // (1) generator: G_ee against T_ee, ratio W_ee, M(pipi) marginals, factorisation residuals, angles
    {
        TCanvas* c = new TCanvas(Nm("cg").c_str(), "", 1800, 1100);
        c->Divide(3, 2, 0.004, 0.006);
        const std::vector<double> fq = LogEdges(2 * kMe, QMax(), 48);
        TH1D* hq = H1v("gq", fq);
        for (size_t i = 0; i < N; ++i) hq->Fill(std::max(G.q[i], 2 * kMe * 1.0000001));
        hq->Scale(1. / N); PerWidth(hq, 1e-3);   // per MeV
        TH1D* tq[4];
        const int dv4[4] = {kV0, kV1, kV2, -1};
        for (int k = 0; k < 3; ++k) {
            tq[k] = H1v("tq", fq);
            for (int b = 0; b < (int)fq.size() - 1; ++b) tq[k]->SetBinContent(b + 1, TeeInt(fq[b], fq[b + 1], dv4[k], 0));
            PerWidth(tq[k], 1e-3);
        }
        c->cd(1); gPad->SetLogx(); gPad->SetLogy(); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
        Sty(hq, kColGen, 2, 1, 20); hq->GetXaxis()->SetTitle("M(e^{+}e^{-}) [GeV/c^{2}]"); hq->GetYaxis()->SetTitle("probability density [1/(MeV/c^{2})]");
        hq->SetMinimum(1e-6); hq->SetMaximum(3); hq->Draw("pe");
        Sty(tq[0], kColV0, 2, 2); Sty(tq[1], kColV1, 2, 3); Sty(tq[2], kColV2, 3, 1);
        for (int k = 0; k < 3; ++k) tq[k]->Draw("hist same");
        TLegend* lg = new TLegend(0.16, 0.17, 0.62, 0.42); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.037);
        lg->AddEntry(hq, "PLUTO generator (eta set)", "pe"); lg->AddEntry(tq[0], "V0: QED x P-wave phase space", "l");
        lg->AddEntry(tq[1], "V1: V0 x |F_{V}(M_{#pi#pi}^{2})|^{2}", "l"); lg->AddEntry(tq[2], "V2 = V3 (main): V1 x |F_{VMD}(M_{ee}^{2})|^{2}", "l"); lg->Draw();
        Txt(0.15, 0.94, "(a) M(e^{+}e^{-}) of #eta #rightarrow #pi^{+}#pi^{-}e^{+}e^{-}, generator level", 0.045);
        Txt(0.17, 0.44, Form("generator: %.1f %% below 2m_{e} (stored precision)", 100. * nBelow / N), 0.033);
        c->cd(2); gPad->SetLogx(); gPad->SetLogy(); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
        TH1D* rq[3];
        for (int k = 0; k < 3; ++k) { rq[k] = (TH1D*)tq[k]->Clone(Nm("rq").c_str()); rq[k]->Divide(hq); }
        Sty(rq[0], kColV0, 2, 2); Sty(rq[1], kColV1, 2, 3); Sty(rq[2], kColV2, 3, 1);
        rq[2]->GetXaxis()->SetTitle("M(e^{+}e^{-}) [GeV/c^{2}]"); rq[2]->GetYaxis()->SetTitle("T_{ee} / G_{ee} (M(ee)-only weight)");
        rq[2]->SetMinimum(0.05); rq[2]->SetMaximum(30); rq[2]->Draw("hist"); rq[0]->Draw("hist same"); rq[1]->Draw("hist same");
        TLine l1; l1.SetLineStyle(2); l1.DrawLine(2 * kMe, 1, QMax(), 1);
        Txt(0.15, 0.94, "(b) ratio target / generator in M(e^{+}e^{-})", 0.045);
        // (c) M(pipi): generator, T variants
        c->cd(3); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
        TH1D* hm = H1v("gm", cmE); for (size_t i = 0; i < N; ++i) hm->Fill(G.m[i]);
        hm->Scale(1. / N); PerWidth(hm);
        TH1D* tm[3];
        for (int k = 0; k < 3; ++k) {
            tm[k] = H1v("tm", cmE);
            for (int b = 0; b < (int)cmE.size() - 1; ++b) tm[k]->SetBinContent(b + 1, TppInt(cmE[b], cmE[b + 1], dv4[k], 0));
            PerWidth(tm[k]);
        }
        Sty(hm, kColGen, 2, 1, 20); hm->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-}) [GeV/c^{2}]"); hm->GetYaxis()->SetTitle("probability density [1/(GeV/c^{2})]");
        hm->SetMaximum(1.35 * std::max(hm->GetMaximum(), tm[2]->GetMaximum())); hm->SetMinimum(0); hm->Draw("pe");
        Sty(tm[0], kColV0, 2, 2); Sty(tm[1], kColV1, 2, 3); Sty(tm[2], kColV2, 3, 1);
        for (int k = 0; k < 3; ++k) tm[k]->Draw("hist same");
        Txt(0.15, 0.94, "(c) M(#pi^{+}#pi^{-}), generator level (linear)", 0.045);
        Txt(0.17, 0.86, Form("mean: gen %.4f, V0 %.4f, V2 %.4f GeV/c^{2}", hm->GetMean(), tm[0]->GetMean(), tm[2]->GetMean()), 0.035);
        // (d) factorisation: observed / model per q slice in M(pipi)
        c->cd(4); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
        const int sl[5][2] = {{0, 4}, {5, 12}, {13, 20}, {21, 28}, {29, nq - 1}};
        const int scol[5] = {1, kColV0, kColV1, kColV2, kColV3};
        TH1D* fr[5];
        TLegend* lf = new TLegend(0.16, 0.15, 0.93, 0.36); lf->SetBorderSize(0); lf->SetFillStyle(0); lf->SetTextSize(0.034); lf->SetNColumns(2);
        for (int s = 0; s < 5; ++s) {
            std::vector<double> me2; for (int b = 0; b <= nm; b += 6) me2.push_back(mE[b]); if (me2.back() < mE[nm]) me2.push_back(mE[nm]);
            fr[s] = H1v("fr", me2); TH1D* fp2 = H1v("fp", me2);
            for (int a = sl[s][0]; a <= sl[s][1]; ++a) for (int b = 0; b < nm; ++b) {
                const double x = 0.5 * (mE[b] + mE[b + 1]);
                fr[s]->Fill(x, O[a * nm + b]); fp2->Fill(x, N * g[a] * h[b] * L[a * nm + b]);
            }
            for (int b = 1; b <= fr[s]->GetNbinsX(); ++b) {
                const double o = fr[s]->GetBinContent(b), p = fp2->GetBinContent(b);
                fr[s]->SetBinContent(b, p > 20 ? o / p : 0); fr[s]->SetBinError(b, p > 20 ? std::sqrt(o) / p : 0);
            }
            Sty(fr[s], scol[s], 1, 1, 20 + s);
            fr[s]->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-}) [GeV/c^{2}]"); fr[s]->GetYaxis()->SetTitle("generator / factorised model");
            fr[s]->SetMinimum(0.75); fr[s]->SetMaximum(1.25);
            fr[s]->Draw(s ? "pe same" : "pe");
            lf->AddEntry(fr[s], Form("M_{ee} %.1f-%.1f MeV", qE[sl[s][0]] * 1e3, qE[sl[s][1] + 1] * 1e3), "pe");
        }
        lf->Draw(); l1.DrawLine(2 * kMpi, 1, kMetaGen, 1);
        Txt(0.15, 0.94, Form("(d) G = g(q) h(M) #lambda^{1/2}: #chi^{2} = %.0f for %d cells", chiTot, ndfTot), 0.045);
        // (e) cos e+
        c->cd(5); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
        TH1D* e1 = (TH1D*)hCe->Clone(Nm("e1").c_str()); e1->Scale(1. / e1->Integral() / 0.05);
        TH1D* e2 = (TH1D*)hCeLo->Clone(Nm("e2").c_str()); e2->Scale(1. / e2->Integral() / 0.05);
        Sty(e1, kColGen, 2, 1, 20); Sty(e2, kColRef, 1, 1, 24);
        e1->GetXaxis()->SetTitle("cos #theta*_{e^{+}} (e^{+}e^{-} rest frame)"); e1->GetYaxis()->SetTitle("probability density");
        e1->SetMinimum(0.3); e1->SetMaximum(0.85); e1->Draw("pe"); e2->Draw("pe same");
        TF1* fl = new TF1(Nm("fl").c_str(), "0.5", -1, 1); fl->SetLineColor(kColV2); fl->SetLineStyle(2); fl->Draw("same");
        TF1* fk = new TF1(Nm("fk").c_str(), "0.375*(1+x*x)", -1, 1); fk->SetLineColor(kColV3); fk->Draw("same");
        TLegend* le = new TLegend(0.16, 0.66, 0.93, 0.89); le->SetBorderSize(0); le->SetFillStyle(0); le->SetTextSize(0.035);
        le->AddEntry(e1, Form("generator, M_{ee} > 2 MeV: b = %.4f #pm %.4f", fe.b, fe.eb), "pe");
        le->AddEntry(e2, Form("generator, M_{ee} < 2 MeV (precision artefact)"), "pe");
        le->AddEntry(fl, "flat", "l"); le->AddEntry(fk, "QED (1 + cos^{2}#theta)/(8/3), M_{ee} >> 2m_{e}", "l"); le->Draw();
        Txt(0.15, 0.94, "(e) e^{+} decay angle, generator level", 0.045);
        // (f) cos pi+
        c->cd(6); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
        TH1D* p1 = (TH1D*)hCp->Clone(Nm("p1").c_str()); p1->Scale(1. / p1->Integral() / 0.05);
        Sty(p1, kColGen, 2, 1, 20);
        p1->GetXaxis()->SetTitle("cos #theta*_{#pi^{+}} (#pi^{+}#pi^{-} rest frame)"); p1->GetYaxis()->SetTitle("probability density");
        p1->SetMinimum(0); p1->SetMaximum(1.0); p1->Draw("pe");
        TF1* fs = new TF1(Nm("fs").c_str(), "0.75*(1-x*x)", -1, 1); fs->SetLineColor(kColV3); fs->Draw("same");
        fl->Draw("same");
        TLegend* lp = new TLegend(0.16, 0.74, 0.93, 0.89); lp->SetBorderSize(0); lp->SetFillStyle(0); lp->SetTextSize(0.035);
        lp->AddEntry(p1, Form("generator: b = %.4f #pm %.4f", fp.b, fp.eb), "pe"); lp->AddEntry(fs, "P wave (3/4) sin^{2}#theta", "l"); lp->AddEntry(fl, "flat", "l"); lp->Draw();
        Txt(0.15, 0.94, "(f) #pi^{+} decay angle, generator level", 0.045);
        c->cd(0);
        Save(c, "rw_weights_generator", true);
    }
    // (2) model: form factors, T_ee vs QED reference, PDG ratio
    {
        TCanvas* c = new TCanvas(Nm("cm").c_str(), "", 1800, 600);
        c->Divide(3, 1, 0.004, 0.006);
        c->cd(1); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
        TH1D* f1 = H1("f1", 280, 0.0, 0.56); TH1D* f2 = H1("f2", 280, 0.0, 0.56); TH1D* f3 = H1("f3", 280, 0.0, 0.56);
        for (int i = 1; i <= 280; ++i) {
            const double x = f1->GetBinCenter(i);
            f1->SetBinContent(i, x > 2 * kMpi ? FV2(x * x, 0) : 0); f2->SetBinContent(i, x > 2 * kMpi ? FV2(x * x, 1) : 0);
            f3->SetBinContent(i, x < QMax() ? FVMD2(x * x) : 0);
        }
        Sty(f1, kColV2, 3); Sty(f2, kColGS, 2, 2); Sty(f3, kColV1, 2, 3);
        f1->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-}) for F_{V}, M(e^{+}e^{-}) for F_{VMD} [GeV/c^{2}]"); f1->GetYaxis()->SetTitle("|F|^{2}");
        f1->SetMinimum(0); f1->SetMaximum(1.15 * std::max(f1->GetMaximum(), f2->GetMaximum())); f1->Draw("hist l"); f2->Draw("hist l same"); f3->Draw("hist l same");
        TLine ln; ln.SetLineStyle(3); ln.DrawLine(kMetaGen, 0, kMetaGen, f1->GetMaximum()); ln.DrawLine(2 * kMpi, 0, 2 * kMpi, f1->GetMaximum());
        TLegend* lg = new TLegend(0.15, 0.62, 0.70, 0.88); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.037);
        lg->AddEntry(f1, "|F_{V}|^{2} Guerrero-Pich (main)", "l"); lg->AddEntry(f2, "|F_{V}|^{2} Gounaris-Sakurai", "l");
        lg->AddEntry(f3, "|F_{VMD}(M_{ee}^{2})|^{2}, #rho pole", "l"); lg->Draw();
        Txt(0.13, 0.94, "(a) form factors; dotted: 2m_{#pi} and m_{#eta}", 0.045);
        c->cd(2); gPad->SetLogx(); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
        TH1D* rr[3]; TH1D* gg = H1v("gg", refE);
        for (int k = 0; k < 3; ++k) rr[k] = H1v("rr", refE);
        for (int b = 0; b < 14 && b < (int)refM.size(); ++b) {
            rr[0]->SetBinContent(b + 1, TeeInt(refE[b], refE[b + 1], kV0, 0) / refM[b]);
            rr[1]->SetBinContent(b + 1, TeeInt(refE[b], refE[b + 1], kV1, 0) / refM[b]);
            rr[2]->SetBinContent(b + 1, TeeInt(refE[b], refE[b + 1], kV2, 0) / refM[b]);
            gg->SetBinContent(b + 1, fg[b] / N / refM[b]);
        }
        Sty(rr[0], kColV0, 2, 2); Sty(rr[1], kColV1, 2, 3); Sty(rr[2], kColV2, 3); Sty(gg, kColGen, 2, 1, 20);
        gPad->SetLogy();
        gg->GetXaxis()->SetTitle("M(e^{+}e^{-}) [GeV/c^{2}]"); gg->GetYaxis()->SetTitle("fraction / earlier QED reference");
        gg->SetMinimum(0.05); gg->SetMaximum(20); gg->Draw("p"); for (int k = 0; k < 3; ++k) rr[k]->Draw("hist same");
        ln.SetLineStyle(2); ln.DrawLine(refE[0], 1, refE[14], 1);
        TLegend* l2 = new TLegend(0.15, 0.66, 0.70, 0.88); l2->SetBorderSize(0); l2->SetFillStyle(0); l2->SetTextSize(0.037);
        l2->AddEntry(gg, "generator", "p"); l2->AddEntry(rr[0], "V0", "l"); l2->AddEntry(rr[1], "V1", "l"); l2->AddEntry(rr[2], "V2 = V3 (main)", "l"); l2->Draw();
        Txt(0.13, 0.94, "(b) M(e^{+}e^{-}) fractions / earlier QED reference", 0.045);
        c->cd(3);
        Txt(0.05, 0.90, "(c) #Gamma(#pi^{+}#pi^{-}e^{+}e^{-}) / #Gamma(#pi^{+}#pi^{-}#gamma)", 0.055);
        Txt(0.05, 0.80, Form("PDG: 2.68#times10^{-4} / 4.22#times10^{-2} = %.2f#times10^{-3}", 1e3 * kBRppee / kBRppg), 0.045);
        const char* lab[5] = {"V0 (no form factor)", "V1 (F_{V} GP)", "V2 (F_{V} GP, F_{VMD})", "V1 (F_{V} GS)", "V2 (F_{V} GS, F_{VMD})"};
        for (int k = 0; k < 5; ++k) Txt(0.05, 0.70 - 0.08 * k, Form("%s: %.2f#times10^{-3} (%.3f of PDG)", lab[k], 1e3 * pr[k], pr[k] / (kBRppee / kBRppg)), 0.045);
        Txt(0.05, 0.24, "same F_{V} in numerator and denominator; the weights", 0.04);
        Txt(0.05, 0.18, "keep the produced BR = 2.68#times10^{-4} (mean weight 1)", 0.04);
        c->cd(0);
        Save(c, "rw_weights_model", false);
    }
    // (3) closure
    {
        TCanvas* c = new TCanvas(Nm("cc").c_str(), "", 1800, 1000);
        c->Divide(3, 2, 0.004, 0.006);
        const int ks[3] = {1, 3, 4};
        const int kc[5] = {kColGen, kColV0, kColV1, kColV2, kColV3};
        // (a) M(ee) density, (b) ratio reweighted/T, (c) M(pipi), (d) cos e+ (V3), (e) cos pi+, (f) weight distribution
        c->cd(1); gPad->SetLogx(); gPad->SetLogy(); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
        TH1D* hq[5]; TH1D* eq[5];
        for (int k = 0; k < 5; ++k) {
            hq[k] = (TH1D*)clo[k].q->Clone(Nm("hq").c_str()); hq[k]->Scale(1. / clo[k].sw); PerWidth(hq[k], 1e-3);
            Sty(hq[k], kc[k], 1, 1, k == 0 ? 24 : 20);
            if (k) { eq[k] = H1v("eq", cqE); for (int b = 0; b < (int)cqE.size() - 1; ++b) eq[k]->SetBinContent(b + 1, expQ[k][b]); PerWidth(eq[k], 1e-3); Sty(eq[k], kc[k], 2); }
        }
        hq[0]->GetXaxis()->SetTitle("M(e^{+}e^{-}) [GeV/c^{2}]"); hq[0]->GetYaxis()->SetTitle("density [1/(MeV/c^{2})]");
        hq[0]->SetMinimum(1e-6); hq[0]->SetMaximum(3); hq[0]->Draw("pe");
        for (int j = 0; j < 3; ++j) { hq[ks[j]]->Draw("pe same"); eq[ks[j]]->Draw("hist same"); }
        TLegend* lg = new TLegend(0.16, 0.17, 0.62, 0.40); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.036);
        lg->AddEntry(hq[0], "generator", "p"); lg->AddEntry(hq[1], "reweighted V0 (points), T_{V0} (line)", "pl");
        lg->AddEntry(hq[3], "reweighted V2, T_{V2}", "pl"); lg->AddEntry(hq[4], "reweighted V3, T_{V3}", "pl"); lg->Draw();
        Txt(0.15, 0.94, "(a) closure in M(e^{+}e^{-})", 0.045);
        c->cd(2); gPad->SetLogx(); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
        TH1D* rq[5];
        for (int j = 0; j < 3; ++j) {
            const int k = ks[j];
            rq[k] = (TH1D*)clo[k].q->Clone(Nm("rq").c_str()); rq[k]->Scale(1. / clo[k].sw);
            for (int b = 1; b <= rq[k]->GetNbinsX(); ++b) { const double e = expQ[k][b - 1]; rq[k]->SetBinContent(b, rq[k]->GetBinContent(b) / e); rq[k]->SetBinError(b, rq[k]->GetBinError(b) / e); }
            Sty(rq[k], kc[k], 1, 1, 20 + j); rq[k]->SetMinimum(0.9); rq[k]->SetMaximum(1.1);
            rq[k]->GetXaxis()->SetTitle("M(e^{+}e^{-}) [GeV/c^{2}]"); rq[k]->GetYaxis()->SetTitle("reweighted generator / T");
            rq[k]->Draw(j ? "pe same" : "pe");
        }
        TLine l1; l1.SetLineStyle(2); l1.DrawLine(2 * kMe, 1, QMax(), 1);
        Txt(0.15, 0.94, Form("(b) ratio; #chi^{2}: V0 %.0f, V2 %.0f, V3 %.0f / 30 bins", clo[1].chiQ, clo[3].chiQ, clo[4].chiQ), 0.045);
        c->cd(3); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
        TH1D* hm[5]; TH1D* em[5];
        for (int k = 0; k < 5; ++k) {
            hm[k] = (TH1D*)clo[k].m->Clone(Nm("hm").c_str()); hm[k]->Scale(1. / clo[k].sw); PerWidth(hm[k]); Sty(hm[k], kc[k], 1, 1, k == 0 ? 24 : 20);
            if (k) { em[k] = H1v("em", cmE); for (int b = 0; b < (int)cmE.size() - 1; ++b) em[k]->SetBinContent(b + 1, expM[k][b]); PerWidth(em[k]); Sty(em[k], kc[k], 2); }
        }
        hm[0]->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-}) [GeV/c^{2}]"); hm[0]->GetYaxis()->SetTitle("density [1/(GeV/c^{2})]");
        hm[0]->SetMinimum(0); hm[0]->SetMaximum(1.3 * std::max(hm[0]->GetMaximum(), em[3]->GetMaximum())); hm[0]->Draw("pe");
        for (int j = 0; j < 3; ++j) { hm[ks[j]]->Draw("pe same"); em[ks[j]]->Draw("hist same"); }
        Txt(0.15, 0.94, Form("(c) M(#pi^{+}#pi^{-}); #chi^{2}: V0 %.0f, V2 %.0f, V3 %.0f / 27", clo[1].chiM, clo[3].chiM, clo[4].chiM), 0.045);
        c->cd(4); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
        TH1D* he[5];
        for (int k : {0, 3, 4}) { he[k] = (TH1D*)clo[k].ce->Clone(Nm("he").c_str()); he[k]->Scale(1. / clo[k].sw / 0.1); Sty(he[k], kc[k], 1, 1, k == 0 ? 24 : 20); }
        TH1D* ee3 = H1("ee3", 20, -1, 1); for (int b = 0; b < 20; ++b) ee3->SetBinContent(b + 1, expE[4][b] / 0.1); Sty(ee3, kColV3, 2);
        he[0]->GetXaxis()->SetTitle("cos #theta*_{e^{+}}"); he[0]->GetYaxis()->SetTitle("probability density"); he[0]->SetMinimum(0.3); he[0]->SetMaximum(0.8);
        he[0]->Draw("pe"); he[3]->Draw("pe same"); he[4]->Draw("pe same"); ee3->Draw("hist same");
        Txt(0.15, 0.94, Form("(d) e^{+} angle; V3 vs T_{V3}: #chi^{2} %.0f / 20 (%.0f for M_{ee} > 2 MeV)", clo[4].chiE, clo[4].chiEHi), 0.042);
        c->cd(5); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
        TH1D* hp[5];
        for (int k : {0, 3, 4}) { hp[k] = (TH1D*)clo[k].cp->Clone(Nm("hp").c_str()); hp[k]->Scale(1. / clo[k].sw / 0.1); Sty(hp[k], kc[k], 1, 1, k == 0 ? 24 : 20); }
        TH1D* ep3 = H1("ep3", 20, -1, 1); for (int b = 0; b < 20; ++b) ep3->SetBinContent(b + 1, expP[4][b] / 0.1); Sty(ep3, kColV3, 2);
        hp[0]->GetXaxis()->SetTitle("cos #theta*_{#pi^{+}}"); hp[0]->GetYaxis()->SetTitle("probability density"); hp[0]->SetMinimum(0); hp[0]->SetMaximum(1.0);
        hp[0]->Draw("pe"); hp[3]->Draw("pe same"); hp[4]->Draw("pe same"); ep3->Draw("hist same");
        Txt(0.15, 0.94, Form("(e) #pi^{+} angle; V3 against sin^{2}: #chi^{2} %.0f / 20", clo[4].chiP), 0.045);
        c->cd(6); gPad->SetLogy(); gPad->SetLogx(); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
        const std::vector<double> we = LogEdges(1e-3, 1e3, 60);
        TH1D* hw[5];
        for (int k = 1; k < 5; ++k) hw[k] = H1v("hw", we);
        for (size_t i = 0; i < N; ++i) for (int k = 1; k < 5; ++k) hw[k]->Fill(std::max(1.0001e-3, W5(G.q[i], G.m[i], dvs[k], 0, G.ce[i], G.cp[i], true)));
        for (int k = 1; k < 5; ++k) { hw[k]->Scale(1. / N); Sty(hw[k], kc[k], 2); }
        hw[1]->GetXaxis()->SetTitle("w_{5} = T / G"); hw[1]->GetYaxis()->SetTitle("fraction of generated decays"); hw[1]->SetMinimum(1e-6); hw[1]->SetMaximum(1);
        hw[1]->Draw("hist"); for (int k = 2; k < 5; ++k) hw[k]->Draw("hist same");
        TLegend* l6 = new TLegend(0.52, 0.70, 0.93, 0.90); l6->SetBorderSize(0); l6->SetFillStyle(0); l6->SetTextSize(0.037);
        for (int k = 1; k < 5; ++k) l6->AddEntry(hw[k], Form("%s: mean %.3f, rms %.2f", vn[k - 1], clo[k].sw / N, std::sqrt(clo[k].sw2 / N - std::pow(clo[k].sw / N, 2))), "l");
        l6->Draw();
        Txt(0.15, 0.94, "(f) weights of the generated decays", 0.045);
        c->cd(0);
        Save(c, "rw_weights_closure", true);
    }
    // (4) bump and conversion weights
    {
        TCanvas* c = new TCanvas(Nm("cb").c_str(), "", 1800, 600);
        c->Divide(3, 1, 0.004, 0.006);
        c->cd(1); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
        TH1D* wk = H1("wk", 300, 2 * kMpi, kMetaHG - kMpi0); TH1D* ws = H1("ws", 300, 2 * kMpi, kMetaHG - kMpi0); TH1D* wn = H1("wn", 300, 2 * kMpi, kMetaHG - kMpi0);
        for (int i = 1; i <= 300; ++i) {
            const double x = wk->GetBinCenter(i);
            wk->SetBinContent(i, WBump(x, kBumpKLOE)); ws->SetBinContent(i, WBump(x, kBumpSlides));
            const double Y = DalY(x); wn->SetBinContent(i, (1 + kDalA * Y + kDalB * Y * Y + kDalF * Y * Y * Y) / normNoX);
        }
        Sty(wk, kColV2, 3); Sty(ws, kColV0, 2, 2); Sty(wn, kColGS, 2, 3);
        wk->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-})_{true} [GeV/c^{2}]"); wk->GetYaxis()->SetTitle("w_{B}"); wk->SetMinimum(0); wk->SetMaximum(3.3);
        wk->Draw("hist l"); ws->Draw("hist l same"); wn->Draw("hist l same");
        TH1D* ps = (TH1D*)hPSm->Clone(Nm("ps").c_str()); ps->Scale(1.5 / ps->GetMaximum()); Sty(ps, kColGen, 1, 3); ps->Draw("hist same");
        TLegend* lg = new TLegend(0.15, 0.65, 0.92, 0.89); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.035);
        lg->AddEntry(wk, Form("KLOE 2016, X^{2} #rightarrow X_{max}^{2}/3 (main), <f>_{PS} = %.4f", normKLOEbar), "l");
        lg->AddEntry(ws, Form("slides: 1 + aY, <f>_{PS} = %.4f", normSlides), "l");
        lg->AddEntry(wn, "X terms dropped (earlier sim_cuts.C variant)", "l");
        lg->AddEntry(ps, "phase-space M(#pi#pi) of #eta#rightarrow#pi^{+}#pi^{-}#pi^{0} (arb.)", "l"); lg->Draw();
        Txt(0.13, 0.94, "(a) bump weight #eta #rightarrow #pi^{+}#pi^{-}#pi^{0} (classes 2, 3)", 0.045);
        c->cd(2); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
        TH1D* xm = H1("xm", nYB, -1, 1); TH1D* xf = H1("xf", nYB, -1, 1);
        for (int k = 0; k < nYB; ++k) if (yW[k] > 0) { xm->SetBinContent(k + 1, yX2[k] / yW[k]); xf->SetBinContent(k + 1, yXf[k] / yW[k]); }
        Sty(xm, kColGen, 1, 1, 20); Sty(xf, kColV2, 2);
        xm->GetXaxis()->SetTitle("Y"); xm->GetYaxis()->SetTitle("<X^{2}> at fixed Y"); xm->SetMinimum(0); xm->SetMaximum(1.4 * xm->GetMaximum());
        xm->Draw("p"); xf->Draw("hist same");
        TLegend* l2 = new TLegend(0.15, 0.75, 0.92, 0.89); l2->SetBorderSize(0); l2->SetFillStyle(0); l2->SetTextSize(0.037);
        l2->AddEntry(xm, "TGenPhaseSpace, 4#times10^{6} decays", "p"); l2->AddEntry(xf, "X_{max}(Y)^{2}/3 (uniform X at fixed Y)", "l"); l2->Draw();
        Txt(0.13, 0.94, "(b) replacement of the X^{2} terms", 0.045);
        c->cd(3); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
        TH1D* c1 = H1("c1", 300, 2 * kMpi, kMetaHG); TH1D* c2 = H1("c2", 300, 2 * kMpi, kMetaHG); TH1D* p3 = H1("p3", 300, 2 * kMpi, kMetaHG);
        for (int i = 1; i <= 300; ++i) { const double x = c1->GetBinCenter(i); c1->SetBinContent(i, WConv(x, 0)); c2->SetBinContent(i, WConv(x, 1)); p3->SetBinContent(i, P3(x)); }
        Sty(c1, kColV2, 3); Sty(c2, kColGS, 2, 2); Sty(p3, kColV1, 2, 3);
        c1->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-})_{true} [GeV/c^{2}]"); c1->GetYaxis()->SetTitle("w_{C}; p_{3}"); c1->SetMinimum(0); c1->SetMaximum(1.15 * std::max(c1->GetMaximum(), c2->GetMaximum()));
        c1->Draw("hist l"); c2->Draw("hist l same"); p3->Draw("hist l same");
        TLegend* l3 = new TLegend(0.15, 0.70, 0.75, 0.89); l3->SetBorderSize(0); l3->SetFillStyle(0); l3->SetTextSize(0.037);
        l3->AddEntry(c1, "w_{C}: #eta #rightarrow #pi^{+}#pi^{-}#gamma, F_{V} GP (class 4)", "l"); l3->AddEntry(c2, "w_{C}, F_{V} GS", "l");
        l3->AddEntry(p3, "p_{3}: #pi^{+}#pi^{-}#pi^{0} share of pion class 2 (HGeant)", "l"); l3->Draw();
        Txt(0.13, 0.94, "(c) conversion weight, pion-class-2 mixture", 0.042);
        c->cd(0);
        Save(c, "rw_weights_bump", false);
    }
    printf("Part1 done in %.1f s\n", sw.RealTime());
}

// =====================================================================================================
// Part 2: true angles of the class-5 rows
// =====================================================================================================
void Part2() {
    TStopwatch sw;
    TFile* fc = TFile::Open(kCompactEta, "READ");
    TTree* tc = fc ? (TTree*)fc->Get("eta_us") : nullptr;
    const std::string outName = std::string(kProdDir) + "/output_eta_pippimepem.root";
    TFile* fs = TFile::Open(outName.c_str(), "READ");
    TTree* ts = fs ? (TTree*)fs->Get("scan_nt") : nullptr;
    if (!tc || !ts) { printf("Part2: input missing\n"); return; }
    const Long64_t n = tc->GetEntries();
    if (ts->GetEntries() != n) { printf("Part2: entries differ: compact %lld, scan_nt %lld\n", n, ts->GetEntries()); return; }
    Row r; BindRow(tc, r, false);
    tc->SetBranchStatus("*", 0);
    for (const char* b : {"grp", "cls", "set", "ep_sim_p", "em_sim_p", "pip_sim_p", "pim_sim_p", "mee_true", "mpipi_true", "m4_true", "ownA"}) tc->SetBranchStatus(b, 1);
    Float_t px[4], py[4], pz[4], sp[4], scls;
    const char* tn[4] = {"ep", "em", "pip", "pim"};
    ts->SetBranchStatus("*", 0);
    for (int k = 0; k < 4; ++k) {
        for (const char* c : {"px", "py", "pz", "p"}) ts->SetBranchStatus(Form("%s_sim_%s", tn[k], c), 1);
        ts->SetBranchAddress(Form("%s_sim_px", tn[k]), &px[k]); ts->SetBranchAddress(Form("%s_sim_py", tn[k]), &py[k]);
        ts->SetBranchAddress(Form("%s_sim_pz", tn[k]), &pz[k]); ts->SetBranchAddress(Form("%s_sim_p", tn[k]), &sp[k]);
    }
    ts->SetBranchStatus("sim_cls", 1); ts->SetBranchAddress("sim_cls", &scls);
    TFile* fo = new TFile(kRwAngFile, "RECREATE");
    TTree* tp = new TTree("ang", "true decay angles of the class-5 rows of eta_us (key: grp + true momenta)");
    TTree& to = *tp;
    Int_t grp; Float_t a, b, c, d; Double_t ce, cp, c2phi, mee, mpp, m4;
    to.Branch("grp", &grp, "grp/I"); to.Branch("ep_sim_p", &a, "ep_sim_p/F"); to.Branch("em_sim_p", &b, "em_sim_p/F");
    to.Branch("pip_sim_p", &c, "pip_sim_p/F"); to.Branch("pim_sim_p", &d, "pim_sim_p/F");
    to.Branch("ce", &ce, "ce/D"); to.Branch("cp", &cp, "cp/D"); to.Branch("c2phi", &c2phi, "c2phi/D");
    to.Branch("mee", &mee, "mee/D"); to.Branch("mpipi", &mpp, "mpipi/D"); to.Branch("m4", &m4, "m4/D");
    long long n5 = 0, bad = 0, badP = 0, badCls = 0;
    double dMax[3] = {0, 0, 0};
    TH1D hce("acc_cos_e", "class 5 rows (all, owned): cos th*_e+", 40, -1, 1), hcp("acc_cos_pi", "class 5 rows (all, owned): cos th*_pi+", 40, -1, 1);
    TH1D hph("acc_c2phi", "class 5 rows (all, owned): cos 2 phi (angle between the decay planes)", 40, -1, 1);
    const double me = simtruth::kMassElectron, mp = simtruth::kMassPion;   // MeV, as sim_truth.h
    for (Long64_t i = 0; i < n; ++i) {
        tc->GetEntry(i);
        if (r.cls != 5) continue;
        ts->GetEntry(i);
        ++n5;
        if (std::lround(scls) != 5) ++badCls;
        if (sp[0] != r.ep_sim_p || sp[1] != r.em_sim_p || sp[2] != r.pip_sim_p || sp[3] != r.pim_sim_p) ++badP;
        TLorentzVector v[4];
        for (int k = 0; k < 4; ++k) v[k].SetXYZM(px[k], py[k], pz[k], k < 2 ? me : mp);
        const TLorentzVector eta = v[0] + v[1] + v[2] + v[3];
        const double rm[3] = {(v[0] + v[1]).M() / 1000., (v[2] + v[3]).M() / 1000., eta.M() / 1000.};
        const double cm[3] = {r.mee_true, r.mpipi_true, r.m4_true};
        for (int k = 0; k < 3; ++k) dMax[k] = std::max(dMax[k], std::fabs(rm[k] - cm[k]));
        if (std::fabs(rm[0] - cm[0]) > 1e-5 || std::fabs(rm[1] - cm[1]) > 1e-5) ++bad;
        const TVector3 be = -eta.BoostVector();
        TLorentzVector w[4]; for (int k = 0; k < 4; ++k) { w[k] = v[k]; w[k].Boost(be); }
        const TLorentzVector ee = w[0] + w[1], pp = w[2] + w[3];
        TLorentzVector e1 = w[0]; e1.Boost(-ee.BoostVector());
        TLorentzVector p1 = w[2]; p1.Boost(-pp.BoostVector());
        ce = e1.Vect().Unit().Dot(ee.Vect().Unit());
        cp = p1.Vect().Unit().Dot(pp.Vect().Unit());
        // angle phi between the e+e- plane (z, e+) and the pi+pi- plane (z, pi+), z = e+e- direction in the eta frame;
        // the components transverse to z are boost invariant along z, so the eta frame is used for both
        {
            const TVector3 z = ee.Vect().Unit();
            const TVector3 ne = z.Cross(w[0].Vect()), np = z.Cross(w[2].Vect());
            const double cphi = (ne.Mag() > 0 && np.Mag() > 0) ? ne.Unit().Dot(np.Unit()) : 1.;
            c2phi = 2 * cphi * cphi - 1;
        }
        grp = r.grp; a = r.ep_sim_p; b = r.em_sim_p; c = r.pip_sim_p; d = r.pim_sim_p; mee = rm[0]; mpp = rm[1]; m4 = rm[2];
        to.Fill();
        if (r.ownA) { hce.Fill(ce); hcp.Fill(cp); hph.Fill(c2phi); }
    }
    fo->cd(); to.Write(); hce.Write(); hcp.Write(); hph.Write();
    fo->Close(); delete fo;
    fc->Close(); fs->Close();
    printf("Part2: eta_us %lld rows, class 5: %lld; scan_nt class != 5: %lld; true-momentum magnitude mismatches: %lld; recomputed masses off by > 10 keV: %lld; "
           "max |diff| mee %.2e, mpipi %.2e, m4 %.2e GeV; %.1f s\n", n, n5, badCls, badP, bad, dMax[0], dMax[1], dMax[2], sw.RealTime());
    FILE* f = fopen((std::string(kScanDir) + "/rw_angles_eta.log").c_str(), "w");
    if (f) {
        fprintf(f, "rw_angles_eta: eta_us %lld rows, class 5: %lld; scan_nt sim_cls != 5 at those entries: %lld; true |p| mismatches compact/scan_nt: %lld;\n"
                   "recomputed masses off by > 10 keV: %lld; max |diff| mee %.2e, mpipi %.2e, m4 %.2e GeV\n", n, n5, badCls, badP, bad, dMax[0], dMax[1], dMax[2]);
        fclose(f);
    }
}

// =====================================================================================================
// Part 3: simulation
// =====================================================================================================
struct DataCls { std::vector<std::string> name; std::vector<double> d, de, simDir, simConv; };
// parse tables/eta_summary.txt section 4 ("# pre, OA [deg] classes" etc.)
bool ParseClasses(const std::string& head, DataCls& out) {
    std::ifstream in(std::string(kTabDir) + "/eta_summary.txt");
    std::string line; bool on = false;
    while (std::getline(in, line)) {
        if (line.find(head) == 0) { on = true; continue; }
        if (!on) continue;
        if (line.rfind("class", 0) == 0) continue;
        if (line.rfind("sum", 0) == 0) break;
        const size_t bar = line.find('|');
        if (bar == std::string::npos) continue;
        std::istringstream ss(line.substr(0, bar));
        std::vector<std::string> t; std::string x;
        while (ss >> x) t.push_back(x);
        if (t.size() < 14) continue;
        out.name.push_back(t[0]);
        out.simDir.push_back(std::atof(t[1].c_str())); out.simConv.push_back(std::atof(t[4].c_str()));
        out.d.push_back(std::atof(t[t.size() - 3].c_str())); out.de.push_back(std::atof(t[t.size() - 1].c_str()));
    }
    return out.d.size() == 5;
}
// data efficiencies (section 5): returns the n-th "x+-y" token of the line starting with cut
bool ParseEff(const std::string& cut, int idx, double& v, double& e) {
    std::ifstream in(std::string(kTabDir) + "/eta_summary.txt");
    std::string line;
    while (std::getline(in, line)) {
        if (line.rfind(cut, 0) != 0) continue;
        std::vector<std::pair<double, double>> vals;
        size_t p = 0;
        while ((p = line.find("+-", p)) != std::string::npos) {
            size_t s = p; while (s > 0 && (std::isdigit(line[s - 1]) || line[s - 1] == '.')) --s;
            size_t q = p + 2; while (q < line.size() && (std::isdigit(line[q]) || line[q] == '.')) ++q;
            vals.push_back(std::make_pair(std::atof(line.substr(s, p - s).c_str()), std::atof(line.substr(p + 2, q - p - 2).c_str())));
            p = q;
        }
        if (idx < (int)vals.size()) { v = vals[idx].first; e = vals[idx].second; return true; }
    }
    return false;
}
struct Fit1 { double a = 0, ea = 0, chi = 0; int ndf = 0; };
Fit1 OneFactor(const std::vector<double>& d, const std::vector<double>& de, const std::vector<double>& s, const std::vector<double>& se) {
    Fit1 f; double a = 0.3;
    for (int it = 0; it < 5; ++it) {
        double sxy = 0, sxx = 0;
        for (size_t i = 0; i < d.size(); ++i) { const double v = de[i] * de[i] + a * a * se[i] * se[i]; sxy += d[i] * s[i] / v; sxx += s[i] * s[i] / v; }
        a = sxy / sxx; f.ea = 1 / std::sqrt(sxx);
    }
    f.a = a;
    for (size_t i = 0; i < d.size(); ++i) { const double v = de[i] * de[i] + a * a * se[i] * se[i]; f.chi += std::pow(d[i] - a * s[i], 2) / v; }
    f.ndf = (int)d.size() - 1;
    return f;
}
struct Fit2 { double a = 0, b = 0, ea = 0, eb = 0, rho = 0, chi = 0; int ndf = 0; };
Fit2 TwoFactor(const std::vector<double>& d, const std::vector<double>& de, const std::vector<double>& s1, const std::vector<double>& s2) {
    double A11 = 0, A12 = 0, A22 = 0, B1 = 0, B2 = 0;
    for (size_t i = 0; i < d.size(); ++i) {
        const double v = de[i] * de[i];
        A11 += s1[i] * s1[i] / v; A12 += s1[i] * s2[i] / v; A22 += s2[i] * s2[i] / v; B1 += d[i] * s1[i] / v; B2 += d[i] * s2[i] / v;
    }
    const double det = A11 * A22 - A12 * A12;
    Fit2 f;
    f.a = (A22 * B1 - A12 * B2) / det; f.b = (A11 * B2 - A12 * B1) / det;
    f.ea = std::sqrt(A22 / det); f.eb = std::sqrt(A11 / det); f.rho = -A12 / std::sqrt(A11 * A22);
    for (size_t i = 0; i < d.size(); ++i) f.chi += std::pow(d[i] - f.a * s1[i] - f.b * s2[i], 2) / (de[i] * de[i]);
    f.ndf = (int)d.size() - 2;
    return f;
}

void Part3() {
    TStopwatch sw;
    Style();
    if (!Init(kRwGenFile)) return;
    Install();
    LoadAngles(kRwAngFile);
    SimData d;
    if (!Load(d, kCompact, kMaskUS)) return;
    printf("Part3: loaded in %.0f s\n", sw.RealTime()); sw.Continue();
    const Sel pre = SelPreset(kPre);
    const Sel selE2 = [](const Row& r) { return PassPre(r) && r.oa_lab < 50.f && r.mpipi < 0.46f && r.oa_eta > 140.f; };
    struct PS { const char* name; Sel s; };
    const std::vector<PS> presets = {{"pre", pre}, {"pre+OA>2", SelPreset(kPreOA2)}, {"E1 (slides)", SelPreset(kSlides)}, {"E2 (slides w/o OA)", selE2}};
    // variants
    struct Var { const char* name; int code; };
    // main = V3 (decay angles) since the review of 2026-11-12; V2 (no angles) was the main variant before
    const std::vector<Var> vars = {{"off", 0},
                                   {"main (V3, KLOE, conv, bkg)", kPhysMain},
                                   {"direct V0", Code(kV0, kBumpKLOE, 1, 1)},
                                   {"direct V1", Code(kV1, kBumpKLOE, 1, 1)},
                                   {"direct V2 (no angles)", kPhysV2},
                                   {"direct V4 (V3 x plane term)", Code(kV4, kBumpKLOE, 1, 1)},
                                   {"bump 1+aY (slides)", Code(kV3, kBumpSlides, 1, 1)},
                                   {"F_V Gounaris-Sakurai", Code(kV3, kBumpKLOE, 1, 1, 1)}};
    const int nV = (int)vars.size();

    // ---------------- weights per class (pre, owned, mode A)
    OpenTab("rw_weights_classes.txt");
    P("# rw_weights_classes: physics-weight factor per class, rows of 'pre' (isBest >= 0, PT3, vz > -160, cor_ok), owned in mode A, std + lep + eta (e+e- trees)\n");
    P("# factor f multiplies the stored weight w (PDG BR); mean_w = sum(w f) / sum(w) = change of the cross section of the class in pre (all M4)\n");
    struct Grp { const char* name; std::function<bool(const Row&)> f; };
    const std::vector<Grp> grps = {
        {"class 2 bump pi0 Dalitz", [](const Row& r) { return r.tmpl == 1 && r.cls == 2; }},
        {"class 3 bump conversion", [](const Row& r) { return r.tmpl == 1 && r.cls == 3; }},
        {"class 4 peak conversion", [](const Row& r) { return r.tmpl == 1 && r.cls == 4; }},
        {"class 5 peak direct (eta set)", [](const Row& r) { return r.tmpl == 1 && r.cls == 5; }},
        {"class 6 X17 (unweighted)", [](const Row& r) { return r.cls == 6; }},
        {"eta set, class 7, lcls 10 (forced e+e-)", [](const Row& r) { return r.set == 3 && r.cls == 7 && r.lcls == 10; }},
        {"eta set, other lcls 10 (cls 0/1/8)", [](const Row& r) { return r.set == 3 && r.cls != 7 && r.cls != 5 && r.lcls == 10; }},
        {"eta set, pcls 9 not cls 5, lcls != 10", [](const Row& r) { return r.set == 3 && r.cls != 5 && r.pcls == 9 && r.lcls != 10; }},
        {"phys. bkg, pcls 2 (std + lep)", [](const Row& r) { return r.tmpl == 2 && r.pcls == 2 && r.set != 3; }},
        {"phys. bkg, pcls 2 (eta set)", [](const Row& r) { return r.tmpl == 2 && r.pcls == 2 && r.set == 3; }}};
    for (int iv = 1; iv < nV; ++iv) {
        P("\n## variant %s (code %d: %s)\n", vars[iv].name, vars[iv].code, PhysLabel(vars[iv].code).c_str());
        P("%-42s %9s %12s %9s %9s %9s %10s %10s\n", "group", "rows", "sum w", "mean_w", "mean", "rms", "min", "max");
        const long long m0 = AngMiss(), h0 = AngHit();
        for (const Grp& gr : grps) {
            double n = 0, sw0 = 0, swf = 0, sf = 0, sf2 = 0, fmin = 1e30, fmax = -1e30;
            for (int set = 1; set <= 3; ++set)
                for (const Row& r : d.tree(set, kUS)) {
                    if (!r.ownA || !pre(r) || !gr.f(r)) continue;
                    const double f = RowFactor(r, vars[iv].code);
                    n += 1; sw0 += r.w; swf += r.w * f; sf += f; sf2 += f * f; fmin = std::min(fmin, f); fmax = std::max(fmax, f);
                }
            if (n == 0) { P("%-42s %9d\n", gr.name, 0); continue; }
            const double mn = sf / n;
            P("%-42s %9.0f %12.4f %9.4f %9.4f %9.4f %10.4g %10.4g\n", gr.name, n, sw0, sw0 > 0 ? swf / sw0 : 0., mn, std::sqrt(std::max(0., sf2 / n - mn * mn)), fmin, fmax);
        }
        if (CodeDir(vars[iv].code) >= kV3) P("# V3/V4: true angles found for %lld class-5 rows, missing for %lld (fall back to V2)\n", AngHit() - h0, AngMiss() - m0);
    }
    CloseTab();
    printf("Part3: classes done %.0f s\n", sw.RealTime()); sw.Continue();

    // ---------------- cross sections
    OpenTab("rw_weights_xsec.txt");
    P("# rw_weights_xsec: simulated cross sections [nb] per component, selection and M4 window, without and with the physics weights\n");
    P("# sigma = sum w f / L_sim (221.05 nb^-1), std + lep (mode A) + eta owned rows, w (PDG); errors from event groups; no efficiency correction\n");
    P("# windows: all M4; bump %.4f-%.4f; peak %.4f-%.4f GeV/c2. Ratios = with / without weights (same rows, correlated: no error quoted)\n", scan::kBumpLo,
      scan::kBumpHi, scan::kPeakLo, scan::kPeakHi);
    struct CompD { const char* name; std::function<bool(const Row&)> f; };
    const std::vector<CompD> cds = {
        {"bump pi0 Dalitz (2)", [](const Row& r) { return r.tmpl == 1 && r.cls == 2; }},
        {"bump conversion (3)", [](const Row& r) { return r.tmpl == 1 && r.cls == 3; }},
        {"peak conversion (4)", [](const Row& r) { return r.tmpl == 1 && r.cls == 4; }},
        {"peak direct (5)", [](const Row& r) { return r.tmpl == 1 && r.cls == 5; }},
        {"peak total (4+5)", [](const Row& r) { return r.tmpl == 1 && (r.cls == 4 || r.cls == 5); }},
        {"bump total (2+3)", [](const Row& r) { return r.tmpl == 1 && (r.cls == 2 || r.cls == 3); }},
        {"eta-set phys. bkg (T2)", [](const Row& r) { return r.tmpl == 2 && r.set == 3; }},
        {"pcls-2 phys. bkg (std+lep)", [](const Row& r) { return r.tmpl == 2 && r.pcls == 2 && r.set != 3; }},
        {"T1 + T2 (simulated sum)", [](const Row& r) { return (r.tmpl == 1 && r.cls != 6) || r.tmpl == 2; }}};
    const int nC = (int)cds.size();
    const double win[3][2] = {{-1e9, 1e9}, {scan::kBumpLo, scan::kBumpHi}, {scan::kPeakLo, scan::kPeakHi}};
    const char* wn[3] = {"all M4", "bump win", "peak win"};
    std::vector<Comp> comps;
    for (int w = 0; w < 3; ++w)
        for (int k = 0; k < nC; ++k) {
            const double lo = win[w][0], hi = win[w][1]; const auto f = cds[k].f;
            comps.push_back(Comp(std::string(cds[k].name) + " " + wn[w], [f, lo, hi](const Row& r) { return f(r) && r.m4 >= lo && r.m4 < hi; }));
        }
    // xs[preset][variant][comp]
    std::vector<std::vector<std::vector<WSum>>> xs(presets.size(), std::vector<std::vector<WSum>>(nV));
    for (size_t ip = 0; ip < presets.size(); ++ip)
        for (int iv = 0; iv < nV; ++iv) { SimOpt o; o.phys = vars[iv].code; xs[ip][iv] = SumMany(d, presets[ip].s, comps, o); }
    for (size_t ip = 0; ip < presets.size(); ++ip) {
        P("\n## %s\n", presets[ip].name);
        P("%-28s %-9s | %18s %18s  ratio | ratios: %s", "component", "window", "sigma off", "sigma main", "");
        for (int iv = 2; iv < nV; ++iv) P(" %s |", vars[iv].name);
        P("\n");
        for (int w = 0; w < 3; ++w)
            for (int k = 0; k < nC; ++k) {
                const int ic = w * nC + k;
                const WSum& a = xs[ip][0][ic]; const WSum& b = xs[ip][1][ic];
                P("%-28s %-9s | %9.4f +- %6.4f %9.4f +- %6.4f  %6.4f |", cds[k].name, wn[w], a.Xsec(), a.XsecErr(), b.Xsec(), b.XsecErr(), a.S > 0 ? b.S / a.S : 0.);
                for (int iv = 2; iv < nV; ++iv) P(" %6.4f", a.S > 0 ? xs[ip][iv][ic].S / a.S : 0.);
                P("\n");
            }
    }
    CloseTab();
    printf("Part3: xsec done %.0f s\n", sw.RealTime()); sw.Continue();

    // ---------------- plane term (V4 / V3) on the accepted direct decay
    OpenTab("rw_weights_planes.txt");
    P("# rw_weights_planes: check of the term in the angle phi between the e+e- and pi+pi- decay planes (variant V4 = V3 x R_phi)\n");
    P("# full angular density of eta -> pi+pi- gamma* -> pi+pi- e+e- (one magnetic amplitude): sin^2 th_pi [1 - beta^2 sin^2 th_e sin^2 phi];\n");
    P("# V3 = its phi average; R_phi = 1 + beta^2 sin^2 th_e cos 2phi / (2 - beta^2 sin^2 th_e); generator phi ASSUMED flat (no phi in the generator file),\n");
    P("# so <R_phi> = 1 at generator level. True phi from the true momenta in scan_nt (rw_angles_eta.root, branch c2phi).\n");
    P("# class-5 rows (direct decay), owned, w x V3 weight; <x>_V3 = sum w w3 x / sum w w3; ratio = sum w w4 / sum w w3 (V4 / V3 cross section)\n");
    P("%-20s %-9s | %9s | %16s %16s | %10s | %18s\n", "selection", "window", "rows", "<cos 2phi>_V3", "<cos 2phi>_unw", "V4/V3", "OA>2 eff V3 / V4");
    {
        const int c3 = Code(kV3, 0, 0, 0), c4 = Code(kV4, 0, 0, 0);
        for (size_t ip = 0; ip < presets.size(); ++ip)
            for (int w = 0; w < 3; ++w) {
                if (w == 1) continue;
                double n = 0, s3 = 0, s4 = 0, sc = 0, sc2 = 0, su = 0, su2 = 0, p3 = 0, p4 = 0, sw3sq = 0;
                for (const Row& r : d.tree(3, kUS)) {
                    if (!r.ownA || r.tmpl != 1 || r.cls != 5 || !presets[ip].s(r)) continue;
                    if (!(r.m4 >= win[w][0] && r.m4 < win[w][1])) continue;
                    double ce = 0, cp = 0, c2p = 0;
                    if (!FindAngles(r, ce, cp, c2p)) continue;
                    const double w3 = r.w * RowFactor(r, c3), w4 = r.w * RowFactor(r, c4);
                    n += 1; s3 += w3; s4 += w4; sc += w3 * c2p; sc2 += w3 * c2p * c2p; su += c2p; su2 += c2p * c2p; sw3sq += w3 * w3;
                    if (r.oa > 2.f) { p3 += w3; p4 += w4; }
                }
                const double m3 = sc / s3, neff = s3 * s3 / sw3sq, e3 = std::sqrt(std::max(0., sc2 / s3 - m3 * m3) / neff);
                const double mu = su / n, eu = std::sqrt(std::max(0., su2 / n - mu * mu) / n);
                P("%-20s %-9s | %9.0f | %+7.4f +- %-6.4f %+7.4f +- %-6.4f | %10.5f | %.4f / %.4f\n", presets[ip].name, wn[w], n, m3, e3, mu, eu, s4 / s3, p3 / s3, p4 / s4);
            }
    }
    P("# |V4/V3 - 1| is the size of the neglected plane term on the accepted rows; the OA > 2 efficiency of the direct decay moves by the same order\n");
    CloseTab();
    printf("Part3: planes done %.0f s\n", sw.RealTime()); sw.Continue();

    // ---------------- peak classes, efficiencies, shapes
    OpenTab("rw_weights_peak.txt");
    P("# rw_weights_peak: the simulated eta peak with the physics weights against the data peak (P_main fits of N - CB per class, stat (+) model errors,\n");
    P("# read from tables/eta_summary.txt section 4). Simulation in expected data rows (sum w f x L_data / L_sim, %s), all M4.\n", "L_data = 6080.5 nb^-1 UNCONFIRMED");
    P("# 'off' reproduces the earlier columns sim dir / sim conv exactly (checked below). Fits: data_i = a1 (sim total_i) (one factor, sim errors\n");
    P("# added in quadrature) and data_i = a (conv_i) + b (direct_i) (data errors only).\n");
    const double kToRows = kLdata / kLsim;
    struct CSet { const char* tag; const char* head; int var; std::vector<double> lo, hi; Sel base; const char* pname; };
    const std::vector<CSet> csets = {
        {"OA5", "# pre, OA [deg] classes", 0, {0, 2, 4, 6, 10}, {2, 4, 6, 10, 1e9}, pre, "pre"},
        {"MEE5", "# pre, M(ee) [MeV/c2] classes", 1, {0, 0.008, 0.020, 0.050, 0.135}, {0.008, 0.020, 0.050, 0.135, 1e9}, pre, "pre"},
        {"OA5", "# E2, OA [deg] classes", 0, {0, 2, 4, 6, 10}, {2, 4, 6, 10, 1e9}, selE2, "E2"},
        {"MEE5", "# E2, M(ee) [MeV/c2] classes", 1, {0, 0.008, 0.020, 0.050, 0.135}, {0.008, 0.020, 0.050, 0.135, 1e9}, selE2, "E2"}};
    // direct variants (codes with only the direct part, plus conversion on/off)
    // index 4 = main (V3); index 3 = V2 (the main variant before the review of 2026-11-12)
    const std::vector<Var> dvars = {{"off", 0}, {"V0", Code(kV0, 0, 1, 0)}, {"V1", Code(kV1, 0, 1, 0)}, {"V2", Code(kV2, 0, 1, 0)},
                                    {"V3 main", Code(kV3, 0, 1, 0)}, {"V3 F_V GS", Code(kV3, 0, 1, 0, 1)}};
    struct ClsRes { std::vector<double> dir[6], dirE[6], conv[3], convE[3]; DataCls data; };
    std::vector<ClsRes> cres(csets.size());
    for (size_t is = 0; is < csets.size(); ++is) {
        const CSet& cs = csets[is];
        ClsRes& cr = cres[is];
        ParseClasses(cs.head, cr.data);
        std::vector<Comp> cc;
        for (size_t i = 0; i < cs.lo.size(); ++i) {
            const double lo = cs.lo[i], hi = cs.hi[i]; const int var = cs.var;
            cc.push_back(Comp("dir", [lo, hi, var](const Row& r) { const double x = var ? r.mee : r.oa; return r.tmpl == 1 && r.cls == 5 && x >= lo && x < hi; }));
            cc.push_back(Comp("conv", [lo, hi, var](const Row& r) { const double x = var ? r.mee : r.oa; return r.tmpl == 1 && r.cls == 4 && x >= lo && x < hi; }));
        }
        for (size_t iv = 0; iv < dvars.size(); ++iv) {
            SimOpt o; o.phys = dvars[iv].code;
            const std::vector<WSum> s = SumMany(d, cs.base, cc, o);
            for (size_t i = 0; i < cs.lo.size(); ++i) {
                cr.dir[iv].push_back(s[2 * i].S * kToRows); cr.dirE[iv].push_back(std::sqrt(s[2 * i].V) * kToRows);
                if (iv == 0 || iv == 3 || iv == 5) {
                    const int k = iv == 0 ? 0 : (iv == 3 ? 1 : 2);
                    cr.conv[k].push_back(s[2 * i + 1].S * kToRows); cr.convE[k].push_back(std::sqrt(s[2 * i + 1].V) * kToRows);
                }
            }
        }
        P("\n## %s, %s classes (%s)\n", cs.pname, cs.tag, cs.var ? "M(ee) [GeV/c2]" : "OA [deg]");
        P("check 'off' against tables/eta_summary.txt: sim dir");
        for (size_t i = 0; i < cs.lo.size(); ++i) P(" %.0f/%.0f", cr.dir[0][i], i < cr.data.simDir.size() ? cr.data.simDir[i] : NAN);
        P("; sim conv"); for (size_t i = 0; i < cs.lo.size(); ++i) P(" %.0f/%.0f", cr.conv[0][i], i < cr.data.simConv.size() ? cr.data.simConv[i] : NAN);
        P("\n%-8s | %8s %8s %8s %8s %8s %8s | %8s %8s | %14s | fractions: data  dir gen  dir V3  conv  conv wC  tot off  tot main(V3)\n", "class", "dir off", "dir V0", "dir V1",
          "dir V2", "dir V3", "dir GS", "conv off", "conv wC", "data");
        // sc[2] = conversion with F_V GS (used by the GS row; fixed 2026-11-12: the array had 2 elements, sc[2] was read out of bounds)
        double sd[6] = {0}, sc[3] = {0}, sdat = 0;
        for (size_t i = 0; i < cs.lo.size(); ++i) { for (int v = 0; v < 6; ++v) sd[v] += cr.dir[v][i]; for (int k = 0; k < 3; ++k) sc[k] += cr.conv[k][i]; sdat += cr.data.d[i]; }
        for (size_t i = 0; i < cs.lo.size(); ++i) {
            P("%-8s | %8.0f %8.0f %8.0f %8.0f %8.0f %8.0f | %8.0f %8.0f | %6.0f +- %5.0f | %.3f+-%.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f\n", cr.data.name[i].c_str(), cr.dir[0][i],
              cr.dir[1][i], cr.dir[2][i], cr.dir[3][i], cr.dir[4][i], cr.dir[5][i], cr.conv[0][i], cr.conv[1][i], cr.data.d[i], cr.data.de[i], cr.data.d[i] / sdat,
              cr.data.de[i] / sdat, cr.dir[0][i] / sd[0], cr.dir[4][i] / sd[4], cr.conv[0][i] / sc[0], cr.conv[1][i] / sc[1],
              (cr.dir[0][i] + cr.conv[0][i]) / (sd[0] + sc[0]), (cr.dir[4][i] + cr.conv[1][i]) / (sd[4] + sc[1]));
        }
        P("%-8s | %8.0f %8.0f %8.0f %8.0f %8.0f %8.0f | %8.0f %8.0f | %6.0f\n", "sum", sd[0], sd[1], sd[2], sd[3], sd[4], sd[5], sc[0], sc[1], sdat);
        // fits
        for (int iv = 0; iv < 6; ++iv) {
            const int kc = (iv == 0) ? 0 : (iv == 5 ? 2 : 1);
            std::vector<double> tot, totE;
            for (size_t i = 0; i < cs.lo.size(); ++i) { tot.push_back(cr.dir[iv][i] + cr.conv[kc][i]); totE.push_back(std::hypot(cr.dirE[iv][i], cr.convE[kc][i])); }
            const Fit1 f1 = OneFactor(cr.data.d, cr.data.de, tot, totE);
            const Fit2 f2 = TwoFactor(cr.data.d, cr.data.de, cr.conv[kc], cr.dir[iv]);
            // (fixed 2026-11-12: the p value of the two-factor fit was printed under the label 'direct share')
            P("  direct %-10s + conv %-5s: one factor a1 = %.3f +- %.3f, chi2/ndf %.1f/%d (p %.3f) | two factors: conv %.3f +- %.3f, direct %.3f +- %.3f (rho %.2f), chi2/ndf %.1f/%d (p %.3f); direct share of the sim peak %.3f\n",
              dvars[iv].name, kc == 0 ? "off" : (kc == 1 ? "wC" : "wC GS"), f1.a, f1.ea, f1.chi, f1.ndf, TMath::Prob(f1.chi, f1.ndf), f2.a, f2.ea, f2.b, f2.eb, f2.rho, f2.chi, f2.ndf,
              TMath::Prob(f2.chi, std::max(1, f2.ndf)), sd[iv] / (sd[iv] + sc[kc]));
        }
    }
    // efficiencies, base pre
    P("\n## single-cut efficiencies, base pre, all M4 (event-group errors); data from tables/eta_summary.txt section 5 (peak: P_main pass/fail, stat (+) model)\n");
    struct CutD { const char* name; Sel pass; const char* key; };
    const std::vector<CutD> cuts = {{"OA > 2", [](const Row& r) { return r.oa > 2.f; }, "OA > 2"},
                                    {"OA < 4", [](const Row& r) { return r.oa < 4.f; }, "OA < 4"},
                                    {"M(pipi) > 0.34", [](const Row& r) { return r.mpipi > 0.34f; }, "M(pipi) > 0.34"},
                                    {"M(pipi) < 0.40", [](const Row& r) { return r.mpipi < 0.40f; }, "M(pipi) < 0.40 (s)"}};
    const Comp cDir = CompPeakDirect(), cConv = CompPeakConv(), cTot = CompPeakTotal(),
               cBump("bump D+C", [](const Row& r) { return r.tmpl == 1 && (r.cls == 2 || r.cls == 3); });
    P("%-16s | %-30s | %-26s | %-26s | %-14s || %-26s | %s\n", "cut", "direct: gen / V0 / V2 / V3", "conversion: off / wC", "total: off / main", "data peak",
      "bump D+C: off / KLOE / 1+aY", "data bump tagged / edge");
    for (const CutD& ct : cuts) {
        const Sel ps = SelAnd(pre, ct.pass);
        auto E = [&](const Comp& c, int code) { SimOpt o; o.phys = code; return Efficiency(d, pre, ps, c, o); };
        const EffRes d0 = E(cDir, 0), d1 = E(cDir, Code(kV0, 0, 0, 0)), d2 = E(cDir, Code(kV2, 0, 0, 0)), d3 = E(cDir, Code(kV3, 0, 0, 0));
        const EffRes c0 = E(cConv, 0), c1 = E(cConv, Code(0, 0, 1, 0));
        const EffRes t0 = E(cTot, 0), t1 = E(cTot, kPhysMain);
        const EffRes b0 = E(cBump, 0), b1 = E(cBump, Code(0, kBumpKLOE, 0, 0)), b2 = E(cBump, Code(0, kBumpSlides, 0, 0));
        double dv = NAN, de = NAN, tv = NAN, te = NAN, ev = NAN, ee = NAN;
        ParseEff(ct.key, 4, dv, de); ParseEff(ct.key, 8, tv, te); ParseEff(ct.key, 9, ev, ee);
        P("%-16s | %.3f / %.3f / %.3f / %.3f    | %.3f+-%.3f / %.3f+-%.3f | %.3f+-%.3f / %.3f+-%.3f | %.3f +- %.3f || %.3f / %.3f / %.3f (+-%.3f) | %.3f+-%.3f / %.3f+-%.3f\n", ct.name,
          d0.eps, d1.eps, d2.eps, d3.eps, c0.eps, c0.err, c1.eps, c1.err, t0.eps, t0.err, t1.eps, t1.err, dv, de, b0.eps, b1.eps, b2.eps, b1.err, tv, te, ev, ee);
    }
    P("# direct-decay efficiency errors are about 0.001 (event groups); data bump columns: tagged eta and edge-step (data study, no vz cut)\n");
    CloseTab();
    printf("Part3: peak done %.0f s\n", sw.RealTime()); sw.Continue();

    // ================= figures
    // (a) weight distributions per class
    {
        TCanvas* c = new TCanvas(Nm("cw").c_str(), "", 1800, 600);
        c->Divide(3, 1, 0.004, 0.006);
        const std::vector<double> we = LogEdges(1e-3, 1e3, 60);
        c->cd(1); gPad->SetLogx(); gPad->SetLogy(); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
        TH1D* hd[4]; const int dvc[4] = {Code(kV0, 0, 0, 0), Code(kV1, 0, 0, 0), Code(kV2, 0, 0, 0), Code(kV3, 0, 0, 0)};
        const int col[4] = {kColV0, kColV1, kColV2, kColV3};
        for (int k = 0; k < 4; ++k) hd[k] = H1v("hd", we);
        double sw0 = 0;
        for (const Row& r : d.tree(3, kUS)) {
            if (!r.ownA || r.tmpl != 1 || r.cls != 5 || !pre(r)) continue;
            sw0 += r.w;
            for (int k = 0; k < 4; ++k) hd[k]->Fill(std::max(1.0001e-3, RowFactor(r, dvc[k])), r.w);
        }
        for (int k = 0; k < 4; ++k) { hd[k]->Scale(1. / sw0); Sty(hd[k], col[k], 2); }
        hd[0]->GetXaxis()->SetTitle("factor w_{5}"); hd[0]->GetYaxis()->SetTitle("fraction of the accepted direct decay (pre)"); hd[0]->SetMinimum(1e-6); hd[0]->SetMaximum(1);
        hd[0]->Draw("hist"); for (int k = 1; k < 4; ++k) hd[k]->Draw("hist same");
        TLegend* lg = new TLegend(0.15, 0.72, 0.6, 0.9); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.037);
        const char* vl[4] = {"V0", "V1", "V2", "V3 (main)"};
        for (int k = 0; k < 4; ++k) lg->AddEntry(hd[k], vl[k], "l");
        lg->Draw();
        Txt(0.13, 0.94, "(a) direct decay (class 5), accepted rows in pre", 0.045);
        c->cd(2); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
        TH1D* hb[3]; const int bc[3] = {Code(0, kBumpKLOE, 0, 0), Code(0, kBumpSlides, 0, 0), Code(0, 0, 1, 0)};
        const int bcol[3] = {kColBumpDalitz, kColV0, kColPeakConv};
        for (int k = 0; k < 3; ++k) hb[k] = H1("hb", 80, 0, 2.5);
        double swb = 0, swc = 0;
        for (int set = 1; set <= 3; ++set)
            for (const Row& r : d.tree(set, kUS)) {
                if (!r.ownA || r.tmpl != 1 || !pre(r)) continue;
                if (r.cls == 2 || r.cls == 3) { swb += r.w; hb[0]->Fill(RowFactor(r, bc[0]), r.w); hb[1]->Fill(RowFactor(r, bc[1]), r.w); }
                if (r.cls == 4) { swc += r.w; hb[2]->Fill(RowFactor(r, bc[2]), r.w); }
            }
        hb[0]->Scale(1. / swb); hb[1]->Scale(1. / swb); hb[2]->Scale(1. / swc);
        for (int k = 0; k < 3; ++k) Sty(hb[k], bcol[k], 2, k == 1 ? 2 : 1);
        hb[0]->GetXaxis()->SetTitle("factor"); hb[0]->GetYaxis()->SetTitle("fraction of the accepted rows (pre, weighted)");
        hb[0]->SetMaximum(1.3 * std::max(hb[0]->GetMaximum(), std::max(hb[1]->GetMaximum(), hb[2]->GetMaximum()))); hb[0]->Draw("hist"); hb[1]->Draw("hist same"); hb[2]->Draw("hist same");
        TLegend* l2 = new TLegend(0.40, 0.72, 0.93, 0.9); l2->SetBorderSize(0); l2->SetFillStyle(0); l2->SetTextSize(0.037);
        l2->AddEntry(hb[0], "bump (2+3), w_{B} KLOE", "l"); l2->AddEntry(hb[1], "bump (2+3), w_{B} = 1 + aY", "l"); l2->AddEntry(hb[2], "peak conversion (4), w_{C}", "l"); l2->Draw();
        Txt(0.13, 0.94, "(b) bump and conversion classes, pre", 0.045);
        // (c) accepted true M(pipi) of the bump and of the conversion, before/after
        c->cd(3); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
        TH1D* mb0 = H1("mb0", 54, 0.275, 0.545); TH1D* mb1 = H1("mb1", 54, 0.275, 0.545); TH1D* mc0 = H1("mc0", 54, 0.275, 0.545); TH1D* mc1 = H1("mc1", 54, 0.275, 0.545);
        for (int set = 1; set <= 3; ++set)
            for (const Row& r : d.tree(set, kUS)) {
                if (!r.ownA || r.tmpl != 1 || !pre(r)) continue;
                if (r.cls == 2 || r.cls == 3) { mb0->Fill(r.mpipi_true, r.w / kLsim); mb1->Fill(r.mpipi_true, r.w * RowFactor(r, bc[0]) / kLsim); }
                if (r.cls == 4) { mc0->Fill(r.mpipi_true, r.w / kLsim); mc1->Fill(r.mpipi_true, r.w * RowFactor(r, bc[2]) / kLsim); }
            }
        for (TH1D* h : {mb0, mb1, mc0, mc1}) PerWidth(h);
        Sty(mb0, kColBumpDalitz, 2, 2); Sty(mb1, kColBumpDalitz, 3); Sty(mc0, kColPeakConv, 2, 2); Sty(mc1, kColPeakConv, 3);
        gPad->SetLogy();
        mb0->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-})_{true} [GeV/c^{2}]"); mb0->GetYaxis()->SetTitle("d#sigma/dM [nb/(GeV/c^{2})]");
        mb0->SetMinimum(0.5); mb0->SetMaximum(5 * std::max(mb0->GetMaximum(), mb1->GetMaximum())); mb0->Draw("hist"); mb1->Draw("hist same"); mc0->Draw("hist same"); mc1->Draw("hist same");
        TLegend* l3 = new TLegend(0.15, 0.72, 0.93, 0.9); l3->SetBorderSize(0); l3->SetFillStyle(0); l3->SetTextSize(0.035);
        l3->AddEntry(mb0, "bump 2+3: HGeant phase space (dashed)", "l"); l3->AddEntry(mb1, "bump 2+3: w_{B} KLOE", "l");
        l3->AddEntry(mc0, "peak conversion (4): phase space (dashed)", "l"); l3->AddEntry(mc1, "peak conversion (4): w_{C}", "l"); l3->Draw();
        Txt(0.13, 0.94, "(c) accepted true M(#pi^{+}#pi^{-}), pre, all M_{4}", 0.045);
        c->cd(0);
        Save(c, "rw_weights_classes", false);
    }
    // (b) cross-section ratios
    {
        TCanvas* c = new TCanvas(Nm("cx").c_str(), "", 1800, 650);
        c->Divide(3, 1, 0.004, 0.006);
        const int show[6] = {0, 1, 2, 3, 4, 5};   // comps: bumpD, bumpC, peakC, peakD, peakT, bumpT
        const char* sn[6] = {"bump Dal.", "bump conv", "peak conv", "peak direct", "peak total", "bump total"};
        const int vcol[8] = {1, kColV3, kColV0, kColV1, kColV2, kColRef, kColBumpDalitz, kColGS};
        for (int w = 0; w < 3; ++w) {
            c->cd(w + 1); gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.16); gPad->SetGridy();
            TH1D* fr = H1("fr", 6 * 4, 0, 6 * 4);
            for (int k = 0; k < 6; ++k) for (size_t ip = 0; ip < presets.size(); ++ip) fr->GetXaxis()->SetBinLabel(k * 4 + ip + 1, ip == 1 ? sn[k] : " ");
            fr->GetXaxis()->LabelsOption("h");
            fr->SetMinimum(0.3); fr->SetMaximum(2.4); fr->GetYaxis()->SetTitle("#sigma with / without physics weights"); fr->GetXaxis()->SetLabelSize(0.042);
            fr->Draw("axis");
            { TLine sep; sep.SetLineStyle(3); sep.SetLineColor(kGray + 1); for (int k = 1; k < 6; ++k) sep.DrawLine(4 * k, 0.3, 4 * k, 2.4); }
            TLegend* lg = new TLegend(0.14, 0.64, 0.62, 0.89); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.033);
            for (int iv = 1; iv < nV; ++iv) {
                TGraph* g = new TGraph();
                for (int k = 0; k < 6; ++k)
                    for (size_t ip = 0; ip < presets.size(); ++ip) {
                        const int ic = w * nC + show[k];
                        const double a = xs[ip][0][ic].S, b = xs[ip][iv][ic].S;
                        if (a > 0) g->SetPoint(g->GetN(), k * 4 + ip + 0.2 + 0.1 * iv, b / a);
                    }
                g->SetMarkerStyle(iv == 1 ? 20 : 24 + iv); g->SetMarkerColor(vcol[iv]); g->SetMarkerSize(iv == 1 ? 1.3 : 1.0);
                g->Draw("p");
                lg->AddEntry(g, vars[iv].name, "p");
            }
            lg->Draw();
            TLine l; l.SetLineStyle(2); l.DrawLine(0, 1, 24, 1);
            Txt(0.12, 0.94, Form("(%c) %s; points per component: pre, pre+OA>2, E1, E2", 'a' + w, wn[w]), 0.042);
        }
        c->cd(0);
        Save(c, "rw_weights_xsec", false);
    }
    // (c) peak: accepted direct-decay M(ee) and OA shapes, and classes against data
    {
        TCanvas* c = new TCanvas(Nm("cp").c_str(), "", 1800, 1100);
        c->Divide(3, 2, 0.004, 0.006);
        const int dvc[5] = {0, Code(kV0, 0, 0, 0), Code(kV1, 0, 0, 0), Code(kV2, 0, 0, 0), Code(kV3, 0, 0, 0)};
        const int col[5] = {kColGen, kColV0, kColV1, kColV2, kColV3};
        const char* vl[5] = {"generator (no weight)", "V0", "V1", "V2", "V3 (main)"};
        const std::vector<double> me = LogEdges(0.001, 0.30, 36);
        TH1D* hm[5]; TH1D* ho[5];
        for (int k = 0; k < 5; ++k) { hm[k] = H1v("hme", me); ho[k] = H1("hoa", 60, 0, 60); }
        for (const Row& r : d.tree(3, kUS)) {
            if (!r.ownA || r.tmpl != 1 || r.cls != 5 || !pre(r)) continue;
            for (int k = 0; k < 5; ++k) {
                const double w = r.w * kToRows * (k ? RowFactor(r, dvc[k]) : 1.);
                hm[k]->Fill(std::max((double)r.mee, 0.001001), w); ho[k]->Fill((double)r.oa, w);
            }
        }
        c->cd(1); gPad->SetLogx(); gPad->SetLogy(); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
        for (int k = 0; k < 5; ++k) { PerWidth(hm[k], 1e-3); Sty(hm[k], col[k], k == 4 ? 3 : 2, k == 0 ? 2 : 1); }
        hm[0]->GetXaxis()->SetTitle("M(e^{+}e^{-}) reconstructed [GeV/c^{2}]"); hm[0]->GetYaxis()->SetTitle("expected rows / (MeV/c^{2})");
        hm[0]->SetMinimum(0.05); hm[0]->SetMaximum(60 * hm[0]->GetMaximum()); hm[0]->Draw("hist"); for (int k = 1; k < 5; ++k) hm[k]->Draw("hist same");
        TLegend* lg = new TLegend(0.40, 0.68, 0.93, 0.90); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.034);
        for (int k = 0; k < 5; ++k) lg->AddEntry(hm[k], Form("%s: %.0f rows", vl[k], hm[k]->Integral("width") * 1e3), "l");
        lg->Draw();
        Txt(0.15, 0.94, "(a) accepted direct decay, M(e^{+}e^{-}), pre, all M_{4}", 0.045);
        Txt(0.16, 0.17, Form("expected data rows, %s", kLdNote), 0.033);
        c->cd(2); gPad->SetLogy(); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
        for (int k = 0; k < 5; ++k) Sty(ho[k], col[k], k == 4 ? 3 : 2, k == 0 ? 2 : 1);
        ho[0]->GetXaxis()->SetTitle("OA(e^{+}e^{-}) [deg]"); ho[0]->GetYaxis()->SetTitle("expected rows / deg");
        ho[0]->SetMinimum(1); ho[0]->SetMaximum(20 * std::max(ho[0]->GetMaximum(), ho[4]->GetMaximum())); ho[0]->Draw("hist"); for (int k = 1; k < 5; ++k) ho[k]->Draw("hist same");
        TLine l; l.SetLineStyle(3); l.DrawLine(2, 1, 2, ho[0]->GetMaximum());
        auto eff2 = [&](TH1D* h) { const int b = h->FindBin(2.0001); return h->Integral(b, h->GetNbinsX()) / h->Integral(); };
        Txt(0.30, 0.85, "fraction with OA > 2#circ (0-60#circ):", 0.035);
        Txt(0.30, 0.80, Form("gen %.3f, V0 %.3f, V2 %.3f, V3 (main) %.3f", eff2(ho[0]), eff2(ho[1]), eff2(ho[3]), eff2(ho[4])), 0.035);
        Txt(0.15, 0.94, "(b) accepted direct decay, OA(e^{+}e^{-}), pre", 0.045);
        // (c)-(f): class fractions, data against simulated total (off / main) and the components
        for (int is = 0; is < 4; ++is) {
            c->cd(3 + is); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
            const ClsRes& cr = cres[is];
            const int n = (int)cr.data.d.size();
            TH1D* hd = H1("hdat", n, 0, n); TH1D* h0 = H1("hs0", n, 0, n); TH1D* h1 = H1("hs1", n, 0, n); TH1D* hD = H1("hD", n, 0, n); TH1D* hC = H1("hC", n, 0, n);
            double sdat = 0, s0 = 0, s1 = 0, sD = 0, sC = 0;
            for (int i = 0; i < n; ++i) { sdat += cr.data.d[i]; s0 += cr.dir[0][i] + cr.conv[0][i]; s1 += cr.dir[4][i] + cr.conv[1][i]; sD += cr.dir[4][i]; sC += cr.conv[1][i]; }
            for (int i = 0; i < n; ++i) {
                hd->SetBinContent(i + 1, cr.data.d[i] / sdat); hd->SetBinError(i + 1, cr.data.de[i] / sdat);
                h0->SetBinContent(i + 1, (cr.dir[0][i] + cr.conv[0][i]) / s0); h1->SetBinContent(i + 1, (cr.dir[4][i] + cr.conv[1][i]) / s1);
                hD->SetBinContent(i + 1, cr.dir[4][i] / sD); hC->SetBinContent(i + 1, cr.conv[1][i] / sC);
                hd->GetXaxis()->SetBinLabel(i + 1, cr.data.name[i].c_str());
            }
            Sty(hd, kBlue, 2, 1, 20); Sty(h0, kColGS, 2, 3); Sty(h1, kColPeakTotal, 3); Sty(hD, kColPeakDirect, 2, 2); Sty(hC, kColPeakConv, 2, 2);
            hd->GetXaxis()->SetLabelSize(0.06); hd->GetXaxis()->SetTitle(csets[is].var ? "M(e^{+}e^{-}) class [MeV/c^{2}]" : "OA(e^{+}e^{-}) class [deg]");
            hd->GetYaxis()->SetTitle("fraction of the peak"); hd->SetMinimum(0); hd->SetMaximum(0.85);
            hd->Draw("pe"); h0->Draw("hist same"); h1->Draw("hist same"); hD->Draw("hist same"); hC->Draw("hist same"); hd->Draw("pe same");
            std::vector<double> tot, totE, tot0, tot0E;
            for (int i = 0; i < n; ++i) {
                tot.push_back(cr.dir[4][i] + cr.conv[1][i]); totE.push_back(std::hypot(cr.dirE[4][i], cr.convE[1][i]));
                tot0.push_back(cr.dir[0][i] + cr.conv[0][i]); tot0E.push_back(std::hypot(cr.dirE[0][i], cr.convE[0][i]));
            }
            const Fit1 f0 = OneFactor(cr.data.d, cr.data.de, tot0, tot0E), f1 = OneFactor(cr.data.d, cr.data.de, tot, totE);
            TLegend* lc = new TLegend(0.33, 0.58, 0.93, 0.90); lc->SetBorderSize(0); lc->SetFillStyle(0); lc->SetTextSize(0.031);
            lc->AddEntry(hd, "data N - CB (P_{main}, stat #oplus model)", "pe");
            lc->AddEntry(h0, Form("sim, no weights: #chi^{2}/ndf %.1f/%d", f0.chi, f0.ndf), "l");
            lc->AddEntry(h1, Form("sim, main weights (V3): #chi^{2}/ndf %.1f/%d", f1.chi, f1.ndf), "l");
            lc->AddEntry(hD, "direct V3 alone", "l"); lc->AddEntry(hC, "conversion w_{C} alone", "l"); lc->Draw();
            Txt(0.15, 0.94, Form("(%c) eta peak per %s class, %s", 'c' + is, csets[is].var ? "M(ee)" : "OA", csets[is].pname), 0.045);
            Txt(0.35, 0.53, Form("one factor data/sim: %.2f #pm %.2f #rightarrow %.2f #pm %.2f", f0.a, f0.ea, f1.a, f1.ea), 0.031);
        }
        c->cd(0);
        Save(c, "rw_weights_peak", true);
    }
    // (d) acceptance in the true M(ee): accepted (pre) / generated, normalised
    {
        TCanvas* c = new TCanvas(Nm("ca").c_str(), "", 1200, 600);
        c->Divide(2, 1, 0.004, 0.006);
        const std::vector<double> qe = LogEdges(2 * kMe, QMax(), 30);
        Gen G; ReadGen(G);
        TH1D* hg = H1v("hgen", qe); for (size_t i = 0; i < G.q.size(); ++i) hg->Fill(std::max(G.q[i], 2 * kMe * 1.0000001));
        hg->Scale(1. / G.q.size());
        TH1D* ha = H1v("hacc", qe); TH1D* hx = H1v("hx17", qe);
        double sa = 0;
        for (const Row& r : d.tree(3, kUS)) {
            if (!r.ownA || r.tmpl != 1 || r.cls != 5 || !pre(r)) continue;
            ha->Fill(std::max((double)r.mee_true, 2 * kMe * 1.0000001), r.w); sa += r.w;
        }
        ha->Scale(1. / sa);
        TH1D* rat = (TH1D*)ha->Clone(Nm("rat").c_str()); rat->Divide(hg);
        // absolute acceptance: accepted direct decays per generated decay = sum w (pre) / (2 BR x N_eta); use the mean of the ratio as scale
        c->cd(1); gPad->SetLogx(); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
        Sty(rat, kColV2, 2, 1, 20); rat->GetXaxis()->SetTitle("M(e^{+}e^{-})_{true} [GeV/c^{2}]");
        rat->GetYaxis()->SetTitle("accepted fraction / generated fraction"); rat->SetMinimum(0); rat->SetMaximum(2.5); rat->Draw("pe");
        TLine l; l.SetLineStyle(2); l.DrawLine(2 * kMe, 1, QMax(), 1);
        Txt(0.13, 0.94, "(a) relative acceptance, direct decay, pre", 0.045);
        Txt(0.15, 0.85, "accepted (class 5, pre, w) / generated (PLUTO sample)", 0.033);
        Txt(0.15, 0.80, "both normalised to 1; first bin: generator values < 2m_{e} included", 0.033);
        c->cd(2); gPad->SetLogx(); gPad->SetLogy(); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
        TH1D* a1 = (TH1D*)ha->Clone(Nm("a1").c_str()); PerWidth(a1, 1e-3); TH1D* g1 = (TH1D*)hg->Clone(Nm("g1").c_str()); PerWidth(g1, 1e-3);
        TH1D* a2 = H1v("a2", qe);
        double s2 = 0;
        for (const Row& r : d.tree(3, kUS)) {
            if (!r.ownA || r.tmpl != 1 || r.cls != 5 || !pre(r)) continue;
            const double w = r.w * RowFactor(r, kPhysMain); a2->Fill(std::max((double)r.mee_true, 2 * kMe * 1.0000001), w); s2 += w;
        }
        a2->Scale(1. / s2); PerWidth(a2, 1e-3);
        TH1D* t2 = H1v("t2", qe); for (int b = 0; b < (int)qe.size() - 1; ++b) t2->SetBinContent(b + 1, TeeInt(qe[b], qe[b + 1], kV2, 0)); PerWidth(t2, 1e-3);
        Sty(g1, kColGen, 2, 2); Sty(a1, kColGen, 2, 1, 24); Sty(a2, kColV2, 2, 1, 20); Sty(t2, kColV2, 2, 2);
        g1->GetXaxis()->SetTitle("M(e^{+}e^{-})_{true} [GeV/c^{2}]"); g1->GetYaxis()->SetTitle("density [1/(MeV/c^{2})]"); g1->SetMinimum(1e-6); g1->SetMaximum(5);
        g1->Draw("hist"); a1->Draw("pe same"); t2->Draw("hist same"); a2->Draw("pe same");
        TLegend* lg = new TLegend(0.40, 0.68, 0.93, 0.90); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.037);
        lg->AddEntry(g1, "generator (dashed)", "l"); lg->AddEntry(a1, "accepted, no weight", "p"); lg->AddEntry(t2, "target T, M(ee) marginal (V2 = V3, dashed)", "l");
        lg->AddEntry(a2, "accepted, main weight (V3)", "p"); lg->Draw();
        Txt(0.13, 0.94, "(b) true M(e^{+}e^{-}) shapes, each normalised", 0.042);
        c->cd(0);
        Save(c, "rw_weights_acceptance", false);
    }
    printf("Part3 done in %.0f s\n", sw.RealTime());
}

}  // namespace srw

void sim_rw_build(int parts = 7) {
    TH1::AddDirectory(kFALSE);
    if (parts & 1) srw::Part1();
    if (parts & 2) srw::Part2();
    if (parts & 4) srw::Part3();
}
