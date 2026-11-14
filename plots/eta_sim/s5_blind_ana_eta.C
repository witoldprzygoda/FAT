// ============================================================================
// s5_blind_ana_eta.C  --  blind analysis of the eta pseudo-data (label blind_ana)
// ============================================================================
// Applies the FROZEN procedure tables/s5_blind_ana_procedure.txt to every toy of
// blind/pseudo/eta_toys.root, treated like measured data. Inputs: the toys and the
// half-A templates (blind/templates/templates.root) only; nothing from blind/sealed.
//
// Per toy and selection (PEAK, DIRECT, BUMP, X17SEL):
//   PEAK yield : M1 (Gaussian frozen + pol2), M2 / M2b (template fits), M2d,
//                M3 (linear side bands) and M3c (calibrated side bands)
//   BUMP yield : B1 (bifurcated Gaussian + pol2), M2 / M2b, B3 (edge step), B4
//                (side band + background shape)
//   DIRECT     : M2d (conv4 constrained), D5 (M(ee) split, DIRECT selection only)
// Errors: event-group (sandwich) stat, Poisson fit, half-A template replicas,
// expected parent (half-B) MC noise; see the procedure file.
//
// Output: /hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/blind/results/eta_results.root
//         (tree "eta", one entry per toy x selection x method x quantity; example
//         spectra of toy 0) and eta_results.txt (per toy, plain text).
// Run (about 1-2 min):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/s5_blind_ana_build", kTRUE)' plots/eta_sim/s5_blind_ana_eta.C+
// ============================================================================
#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TNamed.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TTree.h"

#include "s5_blind_ana_fit.h"
#include "s5_blind_ana_tmpl.h"

using namespace bla;

namespace beta {

const char* const kToyFile = "/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/blind/pseudo/eta_toys.root";
const char* const kOutRoot = "/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/blind/results/eta_results.root";
const char* const kOutTxt  = "/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/blind/results/eta_results.txt";
const double kSc = 0.5;           // L_toy / L_sim
const int kNRep = 50;             // template replicas
const int kNB = 450;              // 2 MeV grid 0.30-1.20
const double kG0 = 0.30, kBW = 0.002;
inline double Xc(int i) { return kG0 + (i + 0.5) * kBW; }
inline int Bin2(float m) { int b = (int)std::floor((m - kG0) / kBW); return (b >= 0 && b < kNB) ? b : -1; }
// D5 grid: 2.5 MeV from 0.483, 50 bins (= m4mee_X17 templates)
const int kNB5 = 50;
const double kG5 = 0.483, kBW5 = 0.0025;
inline int Bin5(float m) { int b = (int)std::floor((m - kG5) / kBW5); return (b >= 0 && b < kNB5) ? b : -1; }

// frozen shapes
const double kBifM[4] = {0.4971, 0.4983, 0.4982, 0.5010}, kBifL[4] = {0.0314, 0.0319, 0.0317, 0.0409}, kBifR[4] = {0.0128, 0.0123, 0.0122, 0.0135};
const double kEdgeE[4] = {0.5072, 0.5090, 0.5090, 0.5060}, kEdgeS[4] = {0.0108, 0.0101, 0.0100, 0.0120};
const double kSBw = 0.0078;                       // M3 side-band width
const double kSBlo = kPeakLo - kSBw, kSBhi = kPeakHi + kSBw;
const double kR3 = (kPeakHi - kPeakLo) / kSBw;    // 25.4 / 7.8
const double kB4lo = 0.5665, kB4hi = 0.7000;      // B4 side band
const double kMeeSplit = 0.010f;                  // D5 split [GeV/c2]

enum { mM1 = 1, mM2, mM2b, mM2d, mM3, mM3c, mB1, mB3, mB4, mD5 };
enum { qPeak = 0, qBump, qDirect, qConv };

struct Row { int evt; float m4, mee, mm4, dz; int pivtx; };

inline bool PassSel(int sel, const Row& r) {
    switch (sel) {
        case 0: return r.dz < 40.f && r.mm4 < 2.85f;
        case 1: return r.dz < 40.f;
        case 2: return true;
        case 3: return r.pivtx == 1 && r.mee * 1000.f >= 13.0 && r.mee * 1000.f < 26.5;
    }
    return false;
}

// ---------------------------------------------------------------------------
// template set of one selection (scaled to L_toy), nominal or replica
// ---------------------------------------------------------------------------
struct TSet {
    std::vector<double> peak, bump, bkg, geta, bkgne, dir, conv, nonpeak, nonbump;   // 2 MeV arrays
    std::vector<double> vpeak, vbump, vbkg, vgeta, vbkgne, vdir, vconv, vnonpeak, vnonbump;
    double peakWin = 0, bumpWin = 0, dirWin = 0, convWin = 0, convWinRel = 0;
    double lamM3 = 0, dB = 0;            // M3 leakage and interpolation bias
    double rhoB4 = 0, pkBump = 0;        // B4 ratio and peak tail in the bump window
    double vNP_win = 0, vNP_L = 0, vNP_U = 0;     // parent variances (non-peak) in M3 regions
    double vNB_win = 0, vNB_SB = 0;               // parent variances (non-bump) in B4 regions
};

// pol2 fit of a 1 MeV array (centres in [lo, hi)), return integrals helper for dB
double InterpBias(const std::vector<double>& c1, const std::vector<double>& e21, double scale) {
    LinProblem P; P.poissonVar = false;
    const double lo = 0.526, hi = 0.566, x0 = 0.546, dx = 0.02;
    for (int b = 0; b < 900; ++b) {
        const double xc = 0.30 + 0.001 * (b + 0.5);
        if (xc < lo || xc >= hi) continue;
        P.y.push_back(scale * c1[b]);
        P.var0.push_back(std::max(scale * scale * e21[b], 1e-6));
        const double t = (xc - x0) / dx;
        P.F.push_back({1.0, t, t * t});
    }
    P.n = P.y.size(); P.np = 3;
    LinResult r = LinFit(P);
    // integral of the polynomial per 1 MeV bin over [x1, x2]
    auto I = [&](double x1, double x2) {
        auto prim = [&](double x) { const double t = (x - x0) / dx; return dx * (r.a[0] * t + r.a[1] * t * t / 2 + r.a[2] * t * t * t / 3); };
        return (prim(x2) - prim(x1)) / 0.001;
    };
    const double bl = I(kSBlo, kPeakLo), bu = I(kPeakHi, kSBhi), bw = I(kPeakLo, kPeakHi);
    return kR3 * 0.5 * (bl + bu) - bw;
}

TSet MakeSet(const M4Templ& t, const std::vector<std::vector<double>>* rep) {
    TSet s;
    s.peak = Sum(t, AtPeak(), 2, kSc, &s.vpeak, rep);
    s.bump = Sum(t, AtBump(), 2, kSc, &s.vbump, rep);
    s.bkg = Sum(t, AtBkg(), 2, kSc, &s.vbkg, rep);
    s.geta = Sum(t, {kGeta}, 2, kSc, &s.vgeta, rep);
    s.bkgne = Sum(t, AtBkgNoEta(), 2, kSc, &s.vbkgne, rep);
    s.dir = Sum(t, {kDirect5}, 2, kSc, &s.vdir, rep);
    s.conv = Sum(t, {kConv4}, 2, kSc, &s.vconv, rep);
    s.nonpeak = Sum(t, {kBump2, kBump3, kGeta, kGpi0D, kGconv, kGvec, kGoth}, 2, kSc, &s.vnonpeak, rep);
    s.nonbump = Sum(t, {kDirect5, kConv4, kGeta, kGpi0D, kGconv, kGvec, kGoth}, 2, kSc, &s.vnonbump, rep);
    if (!rep) {   // exact window sums (reg_ histograms)
        s.peakWin = kSc * (t.regPeak[kDirect5] + t.regPeak[kConv4]);
        s.bumpWin = kSc * (t.regBump[kBump2] + t.regBump[kBump3]);
        s.dirWin = kSc * t.regPeak[kDirect5];
        s.convWin = kSc * t.regPeak[kConv4];
        s.pkBump = kSc * (t.regBump[kDirect5] + t.regBump[kConv4]);
    } else {
        s.peakWin = SumWin(t, AtPeak(), kPeakLo, kPeakHi, kSc, nullptr, rep);
        s.bumpWin = SumWin(t, AtBump(), kBumpLo, kBumpHi, kSc, nullptr, rep);
        s.dirWin = SumWin(t, {kDirect5}, kPeakLo, kPeakHi, kSc, nullptr, rep);
        s.convWin = SumWin(t, {kConv4}, kPeakLo, kPeakHi, kSc, nullptr, rep);
        s.pkBump = SumWin(t, AtPeak(), kBumpLo, kBumpHi, kSc, nullptr, rep);
    }
    s.convWinRel = t.regPeakE[kConv4] / std::max(t.regPeak[kConv4], 1e-9);
    {
        const double Pw = SumWin(t, AtPeak(), kPeakLo, kPeakHi, kSc, nullptr, rep);
        const double PL = SumWin(t, AtPeak(), kSBlo, kPeakLo, kSc, nullptr, rep);
        const double PU = SumWin(t, AtPeak(), kPeakHi, kSBhi, kSc, nullptr, rep);
        s.lamM3 = kR3 * 0.5 * (PL + PU) / Pw;
        // non-peak 1 MeV array
        std::vector<double> c1(900, 0.), e1(900, 0.);
        for (int a : {kBump2, kBump3, kGeta, kGpi0D, kGconv, kGvec, kGoth})
            for (int b = 0; b < 900; ++b) { c1[b] += rep ? (*rep)[a][b] : t.c[a][b]; e1[b] += t.e2[a][b]; }
        s.dB = InterpBias(c1, e1, kSc);
        double e;
        SumWin(t, {kBump2, kBump3, kGeta, kGpi0D, kGconv, kGvec, kGoth}, kPeakLo, kPeakHi, kSc, &e); s.vNP_win = e * e;
        SumWin(t, {kBump2, kBump3, kGeta, kGpi0D, kGconv, kGvec, kGoth}, kSBlo, kPeakLo, kSc, &e); s.vNP_L = e * e;
        SumWin(t, {kBump2, kBump3, kGeta, kGpi0D, kGconv, kGvec, kGoth}, kPeakHi, kSBhi, kSc, &e); s.vNP_U = e * e;
    }
    {
        const double bw = SumWin(t, AtBkg(), kBumpLo, kBumpHi, kSc, nullptr, rep), bs = SumWin(t, AtBkg(), kB4lo, kB4hi, kSc, nullptr, rep);
        s.rhoB4 = bw / bs;
        double e;
        SumWin(t, {kDirect5, kConv4, kGeta, kGpi0D, kGconv, kGvec, kGoth}, kBumpLo, kBumpHi, kSc, &e); s.vNB_win = e * e;
        SumWin(t, {kDirect5, kConv4, kGeta, kGpi0D, kGconv, kGvec, kGoth}, kB4lo, kB4hi, kSc, &e); s.vNB_SB = e * e;
    }
    return s;
}

// ---------------------------------------------------------------------------
// D5 templates: M4 (2.5 MeV, 0.483) projections of m4mee_X17_* for M(ee) < / >= 10 MeV
// ---------------------------------------------------------------------------
struct D5Set { std::vector<double> dL, dH, cL, cH, oL, oH, vdL, vdH, vcL, vcH, voL, voH; };
D5Set LoadD5(TFile& f) {
    D5Set d;
    auto proj = [&](const char* comp, std::vector<double>& L, std::vector<double>& H, std::vector<double>& vL, std::vector<double>& vH) {
        TH2D* h = (TH2D*)f.Get(Form("m4mee_X17_%s", comp));
        L.assign(kNB5, 0.); H.assign(kNB5, 0.); vL.assign(kNB5, 0.); vH.assign(kNB5, 0.);
        for (int i = 1; i <= kNB5; ++i)
            for (int j = 1; j <= h->GetNbinsY(); ++j) {
                const double c = kSc * h->GetBinContent(i, j), e2 = kSc * kSc * std::pow(h->GetBinError(i, j), 2);
                if (h->GetYaxis()->GetBinUpEdge(j) <= 10.0 + 1e-9) { L[i - 1] += c; vL[i - 1] += e2; }
                else { H[i - 1] += c; vH[i - 1] += e2; }
            }
    };
    proj("direct5", d.dL, d.dH, d.vdL, d.vdH);
    proj("conv4", d.cL, d.cH, d.vcL, d.vcH);
    proj("otherX17", d.oL, d.oH, d.voL, d.voH);
    return d;
}

// ---------------------------------------------------------------------------
// model builders (2 MeV grid)
// ---------------------------------------------------------------------------
struct Fit { LinResult r; std::vector<int> idx; bool ok = false; };

Fit FitM1(const std::vector<double>& y, double lo, double hi, int ord, double mean, double sig, double g0, double bw, int nb) {
    Fit f; LinProblem P;
    const double x0 = 0.5 * (lo + hi), dx = 0.5 * (hi - lo);
    for (int i = 0; i < nb; ++i) {
        const double xc = g0 + (i + 0.5) * bw;
        if (xc < lo || xc >= hi) continue;
        f.idx.push_back(i); P.y.push_back(y[i]);
        Vec v(ord + 2);
        v[0] = GaussFrac(g0 + i * bw, g0 + (i + 1) * bw, mean, sig);
        double t = (xc - x0) / dx, tp = 1;
        for (int j = 0; j <= ord; ++j) { v[1 + j] = tp; tp *= t; }
        P.F.push_back(v);
    }
    P.n = f.idx.size(); P.np = ord + 2;
    f.r = LinFit(P); f.ok = f.r.ok;
    return f;
}
// template fit: columns = arrays (2 MeV), with variances, range [lo, hi), optional constraints
Fit FitTmpl(const std::vector<double>& y, const std::vector<const std::vector<double>*>& cols,
            const std::vector<const std::vector<double>*>& vars, double lo, double hi, const std::vector<Constraint>& cons = {}) {
    Fit f; LinProblem P;
    const int np = cols.size();
    for (int i = 0; i < kNB; ++i) {
        if (Xc(i) < lo || Xc(i) >= hi) continue;
        f.idx.push_back(i); P.y.push_back(y[i]);
        Vec v(np), tv(np);
        for (int k = 0; k < np; ++k) { v[k] = (*cols[k])[i]; tv[k] = (*vars[k])[i]; }
        P.F.push_back(v); P.T.push_back(tv);
    }
    P.n = f.idx.size(); P.np = np; P.cons = cons;
    f.r = LinFit(P); f.ok = f.r.ok;
    return f;
}
Fit FitB1(const std::vector<double>& y, int sel) {
    Fit f; LinProblem P;
    const double lo = 0.400, hi = 0.620, x0 = 0.51, dx = 0.11;
    for (int i = 0; i < kNB; ++i) {
        if (Xc(i) < lo || Xc(i) >= hi) continue;
        f.idx.push_back(i); P.y.push_back(y[i]);
        const double x1 = kG0 + i * kBW, x2 = x1 + kBW, t = (Xc(i) - x0) / dx;
        P.F.push_back({BifurFrac(x1, x2, kBifM[sel], kBifL[sel], kBifR[sel]), GaussFrac(x1, x2, kPeakMean, kPeakSigma), 1.0, t, t * t});
    }
    P.n = f.idx.size(); P.np = 5;
    f.r = LinFit(P); f.ok = f.r.ok;
    return f;
}
Fit FitB3(const std::vector<double>& y, int sel) {
    Fit f; LinProblem P;
    const double lo = 0.470, hi = 0.566, x0 = 0.518, dx = 0.048;
    for (int i = 0; i < kNB; ++i) {
        if (Xc(i) < lo || Xc(i) >= hi) continue;
        f.idx.push_back(i); P.y.push_back(y[i]);
        const double x1 = kG0 + i * kBW, x2 = x1 + kBW, t = (Xc(i) - x0) / dx;
        P.F.push_back({StepBin(x1, x2, kEdgeE[sel], kEdgeS[sel]), GaussFrac(x1, x2, kPeakMean, kPeakSigma), 1.0, t, t * t});
    }
    P.n = f.idx.size(); P.np = 5;
    f.r = LinFit(P); f.ok = f.r.ok;
    return f;
}
// weights of g^T a on the full grid (zero outside the fit)
std::vector<double> GridW(const Fit& f, const Vec& g, int nb) {
    std::vector<double> w(nb, 0.);
    const Vec wi = f.r.Weights(g);
    for (size_t q = 0; q < f.idx.size(); ++q) w[f.idx[q]] = wi[q];
    return w;
}
double Dot(const std::vector<double>& w, const std::vector<double>& t) { double s = 0; for (size_t i = 0; i < w.size(); ++i) s += w[i] * t[i]; return s; }
double Dot2(const std::vector<double>& w, const std::vector<double>& v) { double s = 0; for (size_t i = 0; i < w.size(); ++i) s += w[i] * w[i] * v[i]; return s; }

// ---------------------------------------------------------------------------
// one result
// ---------------------------------------------------------------------------
struct Res {
    int sel, meth, qty;
    double Y = 0, raw = 0, cal = 1, eFit = 0, eTmpl = 0, eParent = 0, spread = 0, chi2 = 0; int ndf = 0;
    std::function<double(const Row&)> wrow;   // row weight for the sandwich
    double var = 0;                           // sandwich accumulator
    std::vector<double> rep;                  // replica values
};

void Run() {
    TStopwatch sw; sw.Start();
    TFile ft(kTemplFile);
    M4Templ T[kNSel];
    TSet S[kNSel];
    std::vector<std::vector<TSet>> SR(kNSel);
    TRandom3 rnd(4242);
    for (int s = 0; s < kNSel; ++s) {
        if (!LoadM4(ft, s, T[s])) return;
        S[s] = MakeSet(T[s], nullptr);
        for (int k = 0; k < kNRep; ++k) { auto rp = Replica(T[s], rnd); SR[s].push_back(MakeSet(T[s], &rp)); }
        printf("selection %s: T_peak45(win) %.2f T_bump23(win) %.2f T_direct5(win) %.2f T_conv4(win) %.2f (rel %.3f); lambda_M3 %.4f dB %.2f; rho_B4 %.4f peak in bump win %.2f\n",
               kSel[s], S[s].peakWin, S[s].bumpWin, S[s].dirWin, S[s].convWin, S[s].convWinRel, S[s].lamM3, S[s].dB, S[s].rhoB4, S[s].pkBump);
    }
    D5Set D5 = LoadD5(ft);
    std::vector<D5Set> D5R;
    for (int k = 0; k < kNRep; ++k) {
        D5Set r = D5;
        auto fl = [&](std::vector<double>& a, const std::vector<double>& v) { for (size_t i = 0; i < a.size(); ++i) a[i] += v[i] > 0 ? rnd.Gaus(0, std::sqrt(v[i])) : 0.; };
        fl(r.dL, D5.vdL); fl(r.dH, D5.vdH); fl(r.cL, D5.vcL); fl(r.cH, D5.vcH); fl(r.oL, D5.voL); fl(r.oH, D5.voH);
        D5R.push_back(r);
    }
    std::vector<double> d5DirWinRep(kNRep), d5ConvWinRep(kNRep);
    for (int k = 0; k < kNRep; ++k) {
        d5DirWinRep[k] = S[1].dirWin + rnd.Gaus(0, kSc * T[1].regPeakE[kDirect5]);
        d5ConvWinRep[k] = S[1].convWin + rnd.Gaus(0, kSc * T[1].regPeakE[kConv4]);
    }

    // output
    TFile fo(kOutRoot, "RECREATE");
    TTree tr("eta", "blind analysis of the eta toys (label blind_ana): one entry per toy x selection x method x quantity");
    Int_t oToy, oSel, oMeth, oQty, oNdf;
    Double_t oY, oRaw, oCal, oStat, oFit, oTmpl, oPar, oTot, oSpread, oChi2;
    tr.Branch("toy", &oToy, "toy/I"); tr.Branch("sel", &oSel, "sel/I"); tr.Branch("meth", &oMeth, "meth/I"); tr.Branch("qty", &oQty, "qty/I");
    tr.Branch("Y", &oY, "Y/D"); tr.Branch("raw", &oRaw, "raw/D"); tr.Branch("cal", &oCal, "cal/D");
    tr.Branch("eStat", &oStat, "eStat/D"); tr.Branch("eFit", &oFit, "eFit/D"); tr.Branch("eTmpl", &oTmpl, "eTmpl/D");
    tr.Branch("eParent", &oPar, "eParent/D"); tr.Branch("eTot", &oTot, "eTot/D"); tr.Branch("spread", &oSpread, "spread/D");
    tr.Branch("chi2", &oChi2, "chi2/D"); tr.Branch("ndf", &oNdf, "ndf/I");
    FILE* ftxt = fopen(kOutTxt, "w");
    fprintf(ftxt, "# eta_results.txt -- blind analysis (label blind_ana), frozen procedure plots/eta_sim/tables/s5_blind_ana_procedure.txt\n");
    fprintf(ftxt, "# sel: 0 PEAK 1 DIRECT 2 BUMP 3 X17SEL; meth: 1 M1 2 M2 3 M2b 4 M2d 5 M3 6 M3c 7 B1 8 B3 9 B4 10 D5; qty: 0 peak(4+5 in peak window) 1 bump(2+3 in bump window) 2 direct(5 in peak window) 3 conv(4 in peak window)\n");
    fprintf(ftxt, "# Y = estimate (rows at L_toy = 110.525 nb^-1); eStat event-group; eFit Poisson fit; eTmpl half-A template MC; eParent expected half-B parent MC noise (common to all toys); eTot = sqrt(eStat^2+eTmpl^2)\n");
    fprintf(ftxt, "# toy sel meth qty Y raw cal eStat eFit eTmpl eParent eTot spread chi2 ndf\n");

    // read the toys
    TFile fin(kToyFile);
    TTree* t = (TTree*)fin.Get("eta_toys");
    Int_t toy, evt, pivtx; Float_t m4, mee, mm4, dz;
    t->SetBranchAddress("toy", &toy); t->SetBranchAddress("evt", &evt); t->SetBranchAddress("m4", &m4);
    t->SetBranchAddress("mee", &mee); t->SetBranchAddress("mm4", &mm4); t->SetBranchAddress("dz_pi", &dz); t->SetBranchAddress("pivtx", &pivtx);
    const Long64_t nent = t->GetEntries();
    std::vector<Row> rows; rows.reserve(80000);
    int curToy = -1;
    Long64_t ient = 0;
    const double fwG = GaussFrac(kPeakLo, kPeakHi, kPeakMean, kPeakSigma);

    auto processToy = [&](int itoy) {
        std::vector<Res> R;
        // histograms per selection
        std::vector<std::vector<double>> y(kNSel, std::vector<double>(kNB, 0.));
        for (const Row& r : rows) {
            const int b = Bin2(r.m4); if (b < 0) continue;
            for (int s = 0; s < kNSel; ++s) if (PassSel(s, r)) y[s][b] += 1;
        }
        for (int s = 0; s < kNSel; ++s) {
            const TSet& ts = S[s];
            const std::vector<double>& ys = y[s];
            // ---------------- M1
            {
                Fit f = FitM1(ys, 0.526, 0.614, 2, kPeakMean, kPeakSigma, kG0, kBW, kNB);
                Vec g(4, 0.); g[0] = 1;
                std::vector<double> w = GridW(f, g, kNB);
                const double Rs = Dot(w, ts.peak), c = ts.peakWin / Rs;
                Res r; r.sel = s; r.meth = mM1; r.qty = qPeak;
                r.raw = f.r.a[0] * fwG; r.cal = c; r.Y = f.r.a[0] * c;
                r.eFit = std::sqrt(f.r.C[0][0]) * c; r.chi2 = f.r.chi2P; r.ndf = f.r.ndf;
                r.eParent = std::sqrt(Dot2(w, ts.vnonpeak)) * c;
                for (int k = 0; k < kNRep; ++k) r.rep.push_back(f.r.a[0] * SR[s][k].peakWin / Dot(w, SR[s][k].peak));
                std::vector<double> wc = w; for (auto& x : wc) x *= c;
                r.wrow = [wc, s](const Row& q) { if (!PassSel(s, q)) return 0.; const int b = Bin2(q.m4); return b < 0 ? 0. : wc[b]; };
                // model spread (as the data): 4 ranges x orders 1-4 x widths -0.5/0/+0.5 MeV, p > 0.05
                const double los[4] = {0.522, 0.518, 0.512, 0.508}, his[4] = {0.594, 0.614, 0.624, 0.654};
                double smin = 1e30, smax = -1e30; int npass = 0;
                for (int ir = 0; ir < 4; ++ir) for (int o = 1; o <= 4; ++o) for (int iw = -1; iw <= 1; ++iw) {
                    const double sg = kPeakSigma + iw * 0.0005;
                    Fit q = FitM1(ys, los[ir], his[ir], o, kPeakMean, sg, kG0, kBW, kNB);
                    if (!q.ok || !(q.r.prob > 0.05)) continue;
                    const double v = q.r.a[0] * GaussFrac(kPeakLo, kPeakHi, kPeakMean, sg);
                    smin = std::min(smin, v); smax = std::max(smax, v); npass++;
                }
                r.spread = npass > 1 ? 0.5 * (smax - smin) : -1;
                R.push_back(r);
            }
            // ---------------- M2 and M2b and M2d
            for (int variant = 0; variant < 3; ++variant) {
                auto doFit = [&](const TSet& q) -> Fit {
                    if (variant == 0) return FitTmpl(ys, {&q.peak, &q.bump, &q.bkg}, {&q.vpeak, &q.vbump, &q.vbkg}, 0.400, 0.700);
                    if (variant == 1) return FitTmpl(ys, {&q.peak, &q.bump, &q.geta, &q.bkgne}, {&q.vpeak, &q.vbump, &q.vgeta, &q.vbkgne}, 0.400, 0.700);
                    return FitTmpl(ys, {&q.dir, &q.conv, &q.bump, &q.bkg}, {&q.vdir, &q.vconv, &q.vbump, &q.vbkg}, 0.400, 0.700,
                                   {Constraint{1, 1.0, ts.convWinRel}});
                };
                Fit f = doFit(ts);
                std::vector<Fit> fr; for (int k = 0; k < kNRep; ++k) fr.push_back(doFit(SR[s][k]));
                const int np = f.r.np;
                const int meth = variant == 0 ? mM2 : variant == 1 ? mM2b : mM2d;
                // quantities: peak and bump (M2, M2b); direct, peak, bump (M2d)
                std::vector<int> qs = variant < 2 ? std::vector<int>{qPeak, qBump} : std::vector<int>{qDirect, qPeak, qBump, qConv};
                for (int qq : qs) {
                    Vec g(np, 0.);
                    auto fillG = [&](Vec& gg, const TSet& q) {
                        gg.assign(np, 0.);
                        if (variant < 2) { if (qq == qPeak) gg[0] = q.peakWin; else gg[1] = q.bumpWin; }
                        else { if (qq == qDirect) gg[0] = q.dirWin; else if (qq == qConv) gg[1] = q.convWin; else if (qq == qPeak) { gg[0] = q.dirWin; gg[1] = q.convWin; } else gg[2] = q.bumpWin; }
                    };
                    fillG(g, ts);
                    Res r; r.sel = s; r.meth = meth; r.qty = qq;
                    double Y = 0; for (int k = 0; k < np; ++k) Y += g[k] * f.r.a[k];
                    r.Y = Y; r.raw = Y; r.cal = 1; r.eFit = std::sqrt(f.r.Var(g)); r.chi2 = f.r.chi2P; r.ndf = f.r.ndf;
                    std::vector<double> w = GridW(f, g, kNB);
                    const std::vector<double>& vnt = (qq == qBump) ? ts.vnonbump : ts.vnonpeak;
                    r.eParent = std::sqrt(Dot2(w, vnt));
                    for (int k = 0; k < kNRep; ++k) { Vec gk; fillG(gk, SR[s][k]); double v = 0; for (int j = 0; j < np; ++j) v += gk[j] * fr[k].r.a[j]; r.rep.push_back(v); }
                    r.wrow = [w, s](const Row& q) { if (!PassSel(s, q)) return 0.; const int b = Bin2(q.m4); return b < 0 ? 0. : w[b]; };
                    R.push_back(r);
                }
            }
            // ---------------- M3 and M3c (exact row counts)
            {
                double Nw = 0, NL = 0, NU = 0;
                for (const Row& q : rows) {
                    if (!PassSel(s, q)) continue;
                    if (q.m4 >= kPeakLo && q.m4 < kPeakHi) Nw++;
                    else if (q.m4 >= kSBlo && q.m4 < kPeakLo) NL++;
                    else if (q.m4 >= kPeakHi && q.m4 < kSBhi) NU++;
                }
                const double lam = ts.lamM3, a = 1.0 / (1 - lam);
                const double raw = (Nw - kR3 * 0.5 * (NL + NU)) * a;
                for (int c = 0; c < 2; ++c) {
                    Res r; r.sel = s; r.meth = c == 0 ? mM3 : mM3c; r.qty = qPeak;
                    r.raw = raw; r.Y = c == 0 ? raw : raw - ts.dB * a; r.cal = a;
                    r.eFit = std::sqrt(Nw + kR3 * kR3 * 0.25 * (NL + NU)) * a;
                    r.eParent = std::sqrt(ts.vNP_win + kR3 * kR3 * 0.25 * (ts.vNP_L + ts.vNP_U)) * a;
                    for (int k = 0; k < kNRep; ++k) {
                        const double ak = 1.0 / (1 - SR[s][k].lamM3);
                        const double rk = (Nw - kR3 * 0.5 * (NL + NU)) * ak;
                        r.rep.push_back(c == 0 ? rk : rk - SR[s][k].dB * ak);
                    }
                    r.wrow = [s, a](const Row& q) {
                        if (!PassSel(s, q)) return 0.;
                        if (q.m4 >= kPeakLo && q.m4 < kPeakHi) return a;
                        if ((q.m4 >= kSBlo && q.m4 < kPeakLo) || (q.m4 >= kPeakHi && q.m4 < kSBhi)) return -0.5 * kR3 * a;
                        return 0.;
                    };
                    R.push_back(r);
                }
            }
            // ---------------- B1
            {
                Fit f = FitB1(ys, s);
                Vec g(5, 0.); g[0] = 1;
                std::vector<double> w = GridW(f, g, kNB);
                const double c = ts.bumpWin / Dot(w, ts.bump);
                Res r; r.sel = s; r.meth = mB1; r.qty = qBump;
                r.raw = f.r.a[0] * BifurFrac(kBumpLo, kBumpHi, kBifM[s], kBifL[s], kBifR[s]); r.cal = c; r.Y = f.r.a[0] * c;
                r.eFit = std::sqrt(f.r.C[0][0]) * c; r.chi2 = f.r.chi2P; r.ndf = f.r.ndf;
                r.eParent = std::sqrt(Dot2(w, ts.vnonbump)) * c;
                for (int k = 0; k < kNRep; ++k) r.rep.push_back(f.r.a[0] * SR[s][k].bumpWin / Dot(w, SR[s][k].bump));
                std::vector<double> wc = w; for (auto& x : wc) x *= c;
                r.wrow = [wc, s](const Row& q) { if (!PassSel(s, q)) return 0.; const int b = Bin2(q.m4); return b < 0 ? 0. : wc[b]; };
                R.push_back(r);
                // M1-like peak from the B1 fit is not used (frozen procedure)
            }
            // ---------------- B3 edge step
            {
                Fit f = FitB3(ys, s);
                Vec g(5, 0.); g[0] = 1;
                std::vector<double> w = GridW(f, g, kNB);
                const double c = ts.bumpWin / Dot(w, ts.bump);
                Res r; r.sel = s; r.meth = mB3; r.qty = qBump;
                r.raw = f.r.a[0]; r.cal = c; r.Y = f.r.a[0] * c;
                r.eFit = std::sqrt(f.r.C[0][0]) * std::fabs(c); r.chi2 = f.r.chi2P; r.ndf = f.r.ndf;
                r.eParent = std::sqrt(Dot2(w, ts.vnonbump)) * std::fabs(c);
                for (int k = 0; k < kNRep; ++k) r.rep.push_back(f.r.a[0] * SR[s][k].bumpWin / Dot(w, SR[s][k].bump));
                std::vector<double> wc = w; for (auto& x : wc) x *= c;
                r.wrow = [wc, s](const Row& q) { if (!PassSel(s, q)) return 0.; const int b = Bin2(q.m4); return b < 0 ? 0. : wc[b]; };
                R.push_back(r);
            }
            // ---------------- B4
            {
                double Nw = 0, Ns = 0;
                for (const Row& q : rows) {
                    if (!PassSel(s, q)) continue;
                    if (q.m4 >= kBumpLo && q.m4 < kBumpHi) Nw++;
                    else if (q.m4 >= kB4lo && q.m4 < kB4hi) Ns++;
                }
                const double rho = ts.rhoB4;
                Res r; r.sel = s; r.meth = mB4; r.qty = qBump;
                r.Y = Nw - rho * Ns - ts.pkBump; r.raw = r.Y; r.cal = rho;
                r.eFit = std::sqrt(Nw + rho * rho * Ns);
                r.eParent = std::sqrt(ts.vNB_win + rho * rho * ts.vNB_SB);
                for (int k = 0; k < kNRep; ++k) r.rep.push_back(Nw - SR[s][k].rhoB4 * Ns - SR[s][k].pkBump);
                r.wrow = [s, rho](const Row& q) {
                    if (!PassSel(s, q)) return 0.;
                    if (q.m4 >= kBumpLo && q.m4 < kBumpHi) return 1.;
                    if (q.m4 >= kB4lo && q.m4 < kB4hi) return -rho;
                    return 0.;
                };
                R.push_back(r);
            }
        }
        // ---------------- D5 (DIRECT selection quoted; rows pivtx == 1)
        {
            std::vector<double> yL(kNB5, 0.), yH(kNB5, 0.);
            for (const Row& q : rows) {
                if (q.pivtx != 1) continue;
                const int b = Bin5(q.m4); if (b < 0) continue;
                if (q.mee < kMeeSplit) yL[b] += 1; else yH[b] += 1;
            }
            Fit fL = FitM1(yL, 0.5255, 0.6080, 2, kPeakMean, kPeakSigma, kG5, kBW5, kNB5);
            Fit fH = FitM1(yH, 0.5255, 0.6080, 2, kPeakMean, kPeakSigma, kG5, kBW5, kNB5);
            Vec g(4, 0.); g[0] = 1;
            std::vector<double> wL = GridW(fL, g, kNB5), wH = GridW(fH, g, kNB5);
            auto solve = [&](const D5Set& d, double dirWin, double convWin, double SL, double SH, double& D, double& C, double& aL, double& aH) {
                const double dL = Dot(wL, d.dL) / dirWin, dH = Dot(wH, d.dH) / dirWin;
                const double cL = Dot(wL, d.cL) / convWin, cH = Dot(wH, d.cH) / convWin;
                const double det = dL * cH - cL * dH;
                D = (cH * SL - cL * SH) / det; C = (dL * SH - dH * SL) / det;
                aL = cH / det; aH = -cL / det;   // D = aL SL + aH SH
            };
            const double SL = fL.r.a[0], SH = fH.r.a[0];
            double D, C, aL, aH;
            solve(D5, S[1].dirWin, S[1].convWin, SL, SH, D, C, aL, aH);
            const double dL = Dot(wL, D5.dL) / S[1].dirWin, dH = Dot(wH, D5.dH) / S[1].dirWin;
            const double cL = Dot(wL, D5.cL) / S[1].convWin, cH = Dot(wH, D5.cH) / S[1].convWin;
            const double det = dL * cH - cL * dH;
            const double bL = -dH / det, bH = dL / det;   // C = bL SL + bH SH
            for (int qq : {qDirect, qConv}) {
                Res r; r.sel = 1; r.meth = mD5; r.qty = qq;
                const double xL = qq == qDirect ? aL : bL, xH = qq == qDirect ? aH : bH;
                r.Y = qq == qDirect ? D : C; r.raw = qq == qDirect ? SH : SL; r.cal = det;
                r.eFit = std::sqrt(xL * xL * fL.r.C[0][0] + xH * xH * fH.r.C[0][0]);
                r.chi2 = fL.r.chi2P + fH.r.chi2P; r.ndf = fL.r.ndf + fH.r.ndf;
                r.eParent = std::sqrt(xL * xL * (Dot2(wL, D5.voL)) + xH * xH * (Dot2(wH, D5.voH)));
                for (int k = 0; k < kNRep; ++k) {
                    double Dk, Ck, a1, a2;
                    solve(D5R[k], d5DirWinRep[k], d5ConvWinRep[k], SL, SH, Dk, Ck, a1, a2);
                    r.rep.push_back(qq == qDirect ? Dk : Ck);
                }
                std::vector<double> vL = wL, vH = wH; for (auto& x : vL) x *= xL; for (auto& x : vH) x *= xH;
                r.wrow = [vL, vH](const Row& q) {
                    if (q.pivtx != 1) return 0.;
                    const int b = Bin5(q.m4); if (b < 0) return 0.;
                    return q.mee < kMeeSplit ? vL[b] : vH[b];
                };
                R.push_back(r);
            }
        }
        // ---------------- sandwich over events (rows of one event are adjacent)
        {
            const int nr = R.size();
            std::vector<double> acc(nr, 0.);
            int last = -1;
            for (size_t i = 0; i <= rows.size(); ++i) {
                const bool end = (i == rows.size()) || (rows[i].evt != last);
                if (end && i > 0) { for (int j = 0; j < nr; ++j) { R[j].var += acc[j] * acc[j]; acc[j] = 0; } }
                if (i == rows.size()) break;
                last = rows[i].evt;
                for (int j = 0; j < nr; ++j) acc[j] += R[j].wrow(rows[i]);
            }
        }
        // ---------------- output
        for (Res& r : R) {
            double m = 0, m2 = 0; for (double v : r.rep) { m += v; m2 += v * v; }
            m /= r.rep.size(); r.eTmpl = std::sqrt(std::max(0., m2 / r.rep.size() - m * m));
            oToy = itoy; oSel = r.sel; oMeth = r.meth; oQty = r.qty; oY = r.Y; oRaw = r.raw; oCal = r.cal;
            oStat = std::sqrt(r.var); oFit = r.eFit; oTmpl = r.eTmpl; oPar = r.eParent; oTot = std::sqrt(r.var + r.eTmpl * r.eTmpl);
            oSpread = r.spread; oChi2 = r.chi2; oNdf = r.ndf;
            tr.Fill();
            fprintf(ftxt, "%d %d %d %d %.3f %.3f %.5f %.3f %.3f %.3f %.3f %.3f %.3f %.2f %d\n", oToy, oSel, oMeth, oQty, oY, oRaw, oCal, oStat, oFit, oTmpl, oPar, oTot, oSpread, oChi2, oNdf);
        }
        // example spectra of toy 0
        if (itoy == 0) {
            fo.cd();
            for (int s = 0; s < kNSel; ++s) {
                TH1D h(Form("toy0_m4_%s", kSel[s]), "toy 0: M(#pi^{+}#pi^{-}e^{+}e^{-}), rows per 2 MeV/c^{2}", kNB, 0.30, 1.20);
                for (int i = 0; i < kNB; ++i) { h.SetBinContent(i + 1, y[s][i]); h.SetBinError(i + 1, std::sqrt(y[s][i])); }
                h.Write();
                // fitted models: M1, M2 (and components), B1
                Fit f1 = FitM1(y[s], 0.526, 0.614, 2, kPeakMean, kPeakSigma, kG0, kBW, kNB);
                TH1D hm1(Form("toy0_M1_%s", kSel[s]), "M1 model", kNB, 0.30, 1.20), hm1b(Form("toy0_M1bkg_%s", kSel[s]), "M1 background", kNB, 0.30, 1.20);
                for (size_t q = 0; q < f1.idx.size(); ++q) {
                    hm1.SetBinContent(f1.idx[q] + 1, f1.r.mu[q]);
                    hm1b.SetBinContent(f1.idx[q] + 1, f1.r.mu[q] - f1.r.a[0] * GaussFrac(kG0 + f1.idx[q] * kBW, kG0 + (f1.idx[q] + 1) * kBW, kPeakMean, kPeakSigma));
                }
                hm1.Write(); hm1b.Write();
                Fit f2 = FitTmpl(y[s], {&S[s].peak, &S[s].bump, &S[s].bkg}, {&S[s].vpeak, &S[s].vbump, &S[s].vbkg}, 0.400, 0.700);
                const char* cn[3] = {"peak45", "bump23", "bkgT2"};
                const std::vector<double>* cv[3] = {&S[s].peak, &S[s].bump, &S[s].bkg};
                for (int k = 0; k < 3; ++k) {
                    TH1D hc(Form("toy0_M2_%s_%s", cn[k], kSel[s]), "M2 component", kNB, 0.30, 1.20);
                    for (size_t q = 0; q < f2.idx.size(); ++q) hc.SetBinContent(f2.idx[q] + 1, f2.r.a[k] * (*cv[k])[f2.idx[q]]);
                    hc.Write();
                }
                Fit fb = FitB1(y[s], s);
                TH1D hb(Form("toy0_B1_%s", kSel[s]), "B1 model", kNB, 0.30, 1.20), hbb(Form("toy0_B1bump_%s", kSel[s]), "B1 bump", kNB, 0.30, 1.20);
                for (size_t q = 0; q < fb.idx.size(); ++q) {
                    hb.SetBinContent(fb.idx[q] + 1, fb.r.mu[q]);
                    const double x1 = kG0 + fb.idx[q] * kBW;
                    hbb.SetBinContent(fb.idx[q] + 1, fb.r.a[0] * BifurFrac(x1, x1 + kBW, kBifM[s], kBifL[s], kBifR[s]));
                }
                hb.Write(); hbb.Write();
            }
        }
    };

    for (ient = 0; ient < nent; ++ient) {
        t->GetEntry(ient);
        if (toy != curToy) {
            if (curToy >= 0) { processToy(curToy); if (curToy % 20 == 0) printf("toy %d done (%.0f s)\n", curToy, sw.RealTime()), sw.Continue(); }
            rows.clear(); curToy = toy;
        }
        rows.push_back(Row{evt, m4, mee, mm4, dz, pivtx});
    }
    processToy(curToy);
    fclose(ftxt);
    fo.cd();
    tr.Write();
    TNamed info("info", "blind analysis of the eta toys (label blind_ana, s5_blind_ana_eta.C); frozen procedure plots/eta_sim/tables/s5_blind_ana_procedure.txt; "
                        "sel 0 PEAK 1 DIRECT 2 BUMP 3 X17SEL; meth 1 M1 2 M2 3 M2b 4 M2d 5 M3 6 M3c 7 B1 8 B3 9 B4 10 D5; qty 0 peak 1 bump 2 direct 3 conv; "
                        "Y in rows at L_toy = 110.525 nb^-1; no sealed file read");
    info.Write();
    fo.Close();
    printf("done in %.0f s\n", sw.RealTime());
}

} // namespace beta

void s5_blind_ana_eta() { beta::Run(); }
