// ============================================================================
// wp6_fill.C  --  WP6 (side bands and cross-checks): the parts that need the ROWS
// ============================================================================
// SPEC report_all, WP6. Written by WP6 (prefix wp6), 2026-11-14.
//
// One pass over the data rows (scan_v2 compact) and one over the simulated rows (SimEngine copy, mode A,
// w x V3 = code 108, owned rows, std + lep + eta), selection S0 (isBest >= 0 nominal).
//
//   1. Data rows of S0 (e+e-, e+e+, e-e-) are written to a small tree (grp, M4, M(pipi), M(ee)) for the
//      window counts, event-group covariances and CB replicas of wp6_ana.C.
//   2. Edge step E_step (procedure of plots/scan_bump.C, frozen): N - CB in the M(pi+pi-) slices 0.30-0.34,
//      0.34-0.37, 0.37-0.39, 0.39-0.40 (+ 0.40-0.41 and the control 0.42-0.46), pol2 [+ eta Gaussian] + h x step
//      at the fixed edge 0.99277 M4max(M(pipi)/0.9970), smearing 5.4 MeV; E = sum of h over the four slices.
//      The SAME fit is applied to the simulated sum T1 + T2 at SMASH sigma x PDG BR (expected raw data rows at
//      L = 6079.2 nb^-1 x 63/64; no CB in the simulation). f_eta(edge) = E_data / E_sim.
//      Simulation on the data scale: M4' = 0.99615 M4 (r = 1), M(pipi)' = 0.9952 M(pipi).
//      MC error: 100 Poisson event-group replicas (the WP0 generator and seed, so replica r = hists_v1 replica r).
//   3. ECAL tag (frozen definition of scan_bump.C): good photon, M(e+e-gamma) 0.119-0.173, M(ee) < 0.10;
//      m5c = M(pi+pi-e+e-gamma) with the photon energy rescaled to M(e+e-gamma) = m_pi0. T = Gaussian content
//      of the m5c peak (P_main: frozen 0.5463 / 6.2 MeV, pol2, 0.49-0.62). The simulated ("emulated", UNVERIFIED)
//      tag uses the same row definition; the same fit is applied to the simulated T1 + T2 tagged rows.
//   4. P_main on S0 (frozen 0.5439 / 5.4 MeV, pol2, 0.5225-0.61) with its 48-fit model spread.
//   5. Simulated template sums in exact windows (B, P, side bands, CR-U, CR-L, CR-K range, displaced windows)
//      per (M(pipi) slice x M(ee) class E1/E2) for C2, C3, C4, C5, G_eta, G_D, G_C, T12, with event-group
//      variances and the 100 replicas; 1 MeV template histograms for the figures.
//
// Outputs: /hdd3/przygoda/HADES/FAT_out/all_report/wp6/{wp6_rows.root, wp6_fill.root, wp6_fill.txt}
// Run (from /home/przygoda/HADES/FAT/FAT; about 10-12 GB, 3-5 min):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_wp6", kTRUE)' \
//        'plots/all_report/wp6_fill.C+O'
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <map>
#include <string>
#include <vector>

#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TTree.h"

#include "hists_io.h"
#include "sim_eta_physics.h"

namespace wp6f {
using namespace alldefs;

const char* const kOutDir = "/hdd3/przygoda/HADES/FAT_out/all_report/wp6";
FILE* gOut = nullptr;
#define P(...) do { printf(__VA_ARGS__); if (gOut) fprintf(gOut, __VA_ARGS__); } while (0)

// ---- counter-based Poisson(1) replica weights: copy of wp0_fill.C (same seed -> same replicas) ----------------
inline uint64_t Mix(uint64_t x) {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}
inline int Pois1(uint64_t h) {
    const double u = (double)(h >> 11) * (1.0 / 9007199254740992.0);
    double p = 0.36787944117144233, c = p;
    int k = 0;
    while (u >= c && k < 30) { ++k; p /= k; c += p; }
    return k;
}
inline void GroupPois(uint64_t seed, int tree, int grp, int n, double* k) {
    const uint64_t base = Mix(seed ^ Mix(((uint64_t)(uint32_t)tree << 32) | (uint64_t)(uint32_t)grp));
    for (int r = 0; r < n; ++r) k[r] = Pois1(Mix(base + 0x632BE59BD9B4E019ULL * (uint64_t)(r + 1)));
}
const uint64_t kSeedSim = 0x5157A11ULL;
const int kNRep = 100;

// ---- ECAL tag (frozen, scan_bump.C; identical to wp3_rows.C) ---------------------------------------------------
const double kPi0Lo = 0.1190, kPi0Hi = 0.1730;
const double kMeeMax = 0.100;
const double kTagMean = 0.5463, kTagSigma = 0.0062, kTagLo = 0.49, kTagHi = 0.62;
const int kTagOrder = 2, kTagNbin = 160, kTagReb = 8;
const double kTagHLo = 0.40, kTagHHi = 0.80;
template <class R> inline double M5c(const R& r) {
    const double a = (double)r.m_eeg * r.m_eeg - (double)r.mee * r.mee;
    if (!(a > 0)) return -1.0;
    const double f = (scan::kMpi0 * scan::kMpi0 - (double)r.mee * r.mee) / a;
    const double m2 = (double)r.m4 * r.m4 + f * ((double)r.m5 * r.m5 - (double)r.m4 * r.m4);
    return m2 > 0 ? std::sqrt(m2) : -1.0;
}
template <class R> inline bool TaggedW(const R& r, double lo, double hi) {
    return r.g_idx > 0 && r.mee < kMeeMax && r.m_eeg >= lo && r.m_eeg < hi;
}
inline scan::PeakCfg TagCfg() {
    scan::PeakCfg c;
    c.mean = kTagMean; c.sigma = kTagSigma; c.lo = kTagLo; c.hi = kTagHi; c.order = kTagOrder;
    c.winLo = c.mean - scan::kWinNsig * c.sigma; c.winHi = c.mean + scan::kWinNsig * c.sigma;
    return c;
}
// tag regions: 0 all M4, 1 bump region (M(pipi) < 0.43, M4 < 0.5225), 2 bump window B, 3 CR-K (M(pipi) >= 0.43, M4 0.400-0.5225)
const int kNTagReg = 4;
const char* const kTagRegName[kNTagReg] = {"all M4", "M(pipi) < 0.43, M4 < 0.5225", "bump window B", "CR-K"};
inline bool InTagReg(int g, double m4, double mpp) {
    switch (g) {
        case 0: return true;
        case 1: return mpp < 0.43 && m4 < 0.5225;
        case 2: return m4 >= kBlo && m4 < kBhi;
        case 3: return mpp >= 0.43 && m4 >= 0.400 && m4 < 0.5225;
    }
    return false;
}

// ---- edge-step slices (scan_bump.C) ------------------------------------------------------------------------
const int kNSl = 6, kNSigSl = 4;
const double kSlLo[kNSl] = {0.30, 0.34, 0.37, 0.39, 0.40, 0.42};
const double kSlHi[kNSl] = {0.34, 0.37, 0.39, 0.40, 0.41, 0.46};
const double kCtrlEdge = 0.5150;
inline scan::StepCfg SliceCfg(int i) {
    scan::StepCfg c;
    c.mpLo = kSlLo[i]; c.mpHi = kSlHi[i];
    if (i == 3 || i == 4) c.halfRange = 0.035;
    if (i == 5) { c.useEdge = true; c.edge = kCtrlEdge; }
    return c;
}
// M(ee) classes of the edge and tag spectra: 0 all, 1 E1 (< 10 MeV), 2 E2 (>= 10 MeV)
const int kNEc = 3;
const char* const kEcName[kNEc] = {"all M(ee)", "E1 M(ee) < 10 MeV", "E2 M(ee) >= 10 MeV"};

// ---- simulation groups for the edge spectra --------------------------------------------------------------------
// 0 T12 (simulated sum), 1 C23 (bump classes), 2 other eta (C4 + C5 + G_eta), 3 non-eta background (G_D + G_C)
const int kNEg = 4;
const char* const kEgName[kNEg] = {"T12", "C23", "Oeta", "Bkg"};

// ---- template components of the window sums ---------------------------------------------------------------------
const int kNWc = 8;
const int kWcComp[kNWc] = {kC2, kC3, kC4, kC5, kGeta, kGD, kGconv, kT12};
const char* const kWcTag[kNWc] = {"C2", "C3", "C4", "C5", "Geta", "GD", "GC", "T12"};
inline int WcOf(int comp) { for (int c = 0; c < kNWc; ++c) if (kWcComp[c] == comp) return c; return -1; }

// ---- exact windows of the template sums ------------------------------------------------------------------------
struct Win { std::string tag; double lo, hi; };
std::vector<Win> gW;
// peak-relative offsets of P and the side bands (P = 0.5439 -+ 1.9 x 5.4 MeV; SB 0.5250-0.5336, 0.5542-0.5628)
const double kPc = 0.5439;
const double kDPlo = kPlo - kPc, kDPhi = kPhi - kPc;
const double kDSl1 = 0.5250 - kPc, kDSl2 = 0.5336 - kPc, kDSh1 = 0.5542 - kPc, kDSh2 = 0.5628 - kPc;
const int kNDisp = 13;
inline double DispCentre(int k) { return 0.575 + 0.010 * k; }
void InitWindows() {
    gW.clear();
    gW.push_back({"B", kBlo, kBhi});                 // 0
    gW.push_back({"P", kPlo, kPhi});                 // 1
    gW.push_back({"SBlo", 0.5250, 0.5336});          // 2
    gW.push_back({"SBhi", 0.5542, 0.5628});          // 3
    gW.push_back({"U575", 0.575, 0.700});            // 4
    gW.push_back({"U565", 0.565, 0.700});            // 5
    gW.push_back({"CRL", 0.5225, 0.5300});           // 6
    gW.push_back({"K400", 0.400, 0.5225});           // 7  (CR-K = slice S5 of this range)
    gW.push_back({"R0", 0.400, 0.455});              // 8
    gW.push_back({"Fit", 0.400, 0.700});             // 9
    gW.push_back({"SBloN", 0.5280, 0.5336});         // 10 narrow side bands (5.6 MeV)
    gW.push_back({"SBhiN", 0.5542, 0.5598});         // 11
    gW.push_back({"SBloW", 0.5220, 0.5336});         // 12 wide side bands (11.6 MeV)
    gW.push_back({"SBhiW", 0.5542, 0.5658});         // 13
    gW.push_back({"tot", -1e9, 1e9});                // 14 all M4
    for (int k = 0; k < kNDisp; ++k) {               // 15 + 3k: P_k, SBlo_k, SBhi_k
        const double c = DispCentre(k);
        gW.push_back({Form("P_d%d", k), c + kDPlo, c + kDPhi});
        gW.push_back({Form("SBlo_d%d", k), c + kDSl1, c + kDSl2});
        gW.push_back({Form("SBhi_d%d", k), c + kDSh1, c + kDSh2});
    }
}
const int kNCell = 10;   // slice (5) x M(ee) class E1/E2 (2): cell = 2 slice + class

// ---- accumulators with event-group variance (sum over groups of the squared group sums) --------------------
struct Acc {
    std::vector<double> sum, var, tmp;
    std::vector<int> touched;
    void Init(size_t n) { sum.assign(n, 0.); var.assign(n, 0.); tmp.assign(n, 0.); touched.clear(); }
    inline void Add(size_t i, double w) {
        if (tmp[i] == 0.) touched.push_back((int)i);
        tmp[i] += w; sum[i] += w;
        if (tmp[i] == 0.) tmp[i] = 1e-300;   // keep it marked
    }
    void Flush() {
        for (int i : touched) { var[i] += tmp[i] * tmp[i]; tmp[i] = 0.; }
        touched.clear();
    }
};
struct RepAcc {
    std::vector<double> r;   // r[i * kNRep + q]
    void Init(size_t n) { r.assign(n * kNRep, 0.); }
    inline void Add(size_t i, double w, const double* k) { double* p = &r[i * kNRep]; for (int q = 0; q < kNRep; ++q) p[q] += w * k[q]; }
};

// ---- a scan::Spectra holding a simulated spectrum (no CB) ---------------------------------------------------
scan::Spectra SimSpectra(const double* y, int nbin, double lo, double hi, int reb) {
    scan::Spectra s;
    s.nbin = nbin; s.reb = reb; s.nw = nbin / reb; s.lo = lo; s.hi = hi; s.bw = (hi - lo) / nbin;
    for (int k = 0; k < 3; ++k) { s.n[k].assign(nbin, 0.); s.w[k].assign(s.nw, 0.); }
    for (int i = 0; i < nbin; ++i) s.n[0][i] = y[i];
    for (int i = 0; i < nbin; ++i) s.w[0][i / reb] += y[i];
    s.cbW.assign(s.nw, 0.); s.cbWVar.assign(s.nw, 0.);
    s.cb.assign(nbin, 0.); s.cbVar.assign(nbin, 0.);
    s.sig.assign(y, y + nbin); s.sigVar.assign(y, y + nbin);
    for (int k = 0; k < 3; ++k) s.nsel[k] = 0;
    for (int i = 0; i < nbin; ++i) s.nsel[0] += y[i];
    return s;
}
// spread-rule CB of a class spectrum: CB of the whole slice x f(class | 50 MeV band) from N++ + N-- (WP3 SpreadCB)
inline int BandOfX(double x) { return (int)std::floor((x - 0.30) / 0.05 + 1e-9); }
scan::Spectra SpreadCB(const scan::Spectra& all, const scan::Spectra& sl) {
    scan::Spectra s = sl;
    std::map<int, double> la, ls;
    for (int i = 0; i < sl.nbin; ++i) {
        const int b = BandOfX(sl.x(i));
        la[b] += all.n[1][i] + all.n[2][i];
        ls[b] += sl.n[1][i] + sl.n[2][i];
    }
    auto f = [&](double x) { const int b = BandOfX(x); return la[b] > 0 ? ls[b] / la[b] : 0.; };
    for (int i = 0; i < s.nbin; ++i) {
        const double fi = f(s.x(i));
        s.cb[i] = all.cb[i] * fi;
        s.cbVar[i] = all.cbVar[i] * fi * fi;
        s.sig[i] = s.n[0][i] - s.cb[i];
        s.sigVar[i] = s.n[0][i] + s.cbVar[i];
    }
    for (int j = 0; j < s.nw; ++j) {
        const double xc = s.lo + (j + 0.5) * s.reb * s.bw;
        const double fj = f(xc);
        s.cbW[j] = all.cbW[j] * fj;
        s.cbWVar[j] = all.cbWVar[j] * fj * fj;
    }
    return s;
}

// ---- histogram helpers ---------------------------------------------------------------------------------------
TH1D* H1(const char* name, const char* title, int n, double lo, double hi, const double* y, const double* v = nullptr) {
    TH1D* h = new TH1D(name, title, n, lo, hi);
    h->SetDirectory(nullptr);
    for (int i = 0; i < n; ++i) { h->SetBinContent(i + 1, y[i]); if (v) h->SetBinError(i + 1, std::sqrt(std::max(0., v[i]))); }
    return h;
}
void Write(TDirectory* d, TObject* o) { d->WriteTObject(o, o->GetName(), "Overwrite"); delete o; }
TGraph* Curve(const std::vector<double>& x, const std::vector<double>& y, const char* name) {
    TGraph* g = new TGraph((int)x.size());
    for (size_t i = 0; i < x.size(); ++i) g->SetPoint((int)i, x[i], y[i]);
    g->SetName(name);
    return g;
}
TGraph* PeakCurve(const scan::PeakFit& f, double bw, const char* name, bool bkgOnly) {
    TGraph* g = new TGraph();
    int n = 0;
    for (double m = f.cfg.lo + 0.5 * bw; m < f.cfg.hi; m += 0.25 * bw) g->SetPoint(n++, m, bkgOnly ? f.Bkg(m) : f.Total(m));
    g->SetName(name);
    return g;
}

// ---- edge-step sums ------------------------------------------------------------------------------------------
struct ERes { double E = 0, eE = 0; double h[kNSl], eh[kNSl], chi2[kNSl]; int ndf[kNSl]; bool ok[kNSl]; double chi2s = 0; int ndfs = 0; };
ERes EdgeSum(const std::vector<scan::Spectra>& sp, const std::function<void(int, scan::StepCfg&)>& mod = nullptr,
             std::vector<scan::StepFit>* keep = nullptr) {
    ERes r;
    double v = 0;
    for (int i = 0; i < kNSl; ++i) {
        scan::StepCfg c = SliceCfg(i);
        if (mod) mod(i, c);
        const scan::StepFit f = scan::FitEdgeStep(sp[i], c);
        if (keep) keep->push_back(f);
        r.ok[i] = f.ok; r.h[i] = f.h; r.eh[i] = f.eh; r.chi2[i] = f.chi2; r.ndf[i] = f.ndf;
        if (i < kNSigSl && f.ok) { r.E += f.h; v += f.eh * f.eh; r.chi2s += f.chi2; r.ndfs += f.ndf; }
    }
    r.eE = std::sqrt(v);
    return r;
}

}  // namespace wp6f

void wp6_fill() {
    using namespace alldefs;
    using namespace wp6f;
    TStopwatch total;
    gSystem->mkdir(kOutDir, kTRUE);
    gOut = fopen(Form("%s/wp6_fill.txt", kOutDir), "w");
    InitWindows();
    const int nW = (int)gW.size();
    P("# wp6_fill.C (WP6 row-level part), 2026-11-14. Selection S0 = pre + both pions in the vertex fit (isBest >= 0).\n");
    P("# data: %s\n# sim : %s + %s (mode A, w x V3 code %d, owned rows); M4' = %.5f M4 (r = 1), M(pipi)' = %.4f M(pipi)\n",
      scan::kCompact, simeng::kCompact, simeng::kCompactEta, kPhysNom, kS, kSPiPi);
    P("# simulation in expected RAW data rows: sum w x %.4f (L = %.1f nb^-1 x 63/64, L_sim = %.2f nb^-1)\n\n", kSimToRaw, kL, kLsim);

    TFile* fo = TFile::Open(Form("%s/wp6_fill.root", kOutDir), "RECREATE");
    TDirectory* dData = fo->mkdir("data");
    TDirectory* dSim = fo->mkdir("sim");

    // =========================================================================================================
    // A. DATA
    // =========================================================================================================
    std::vector<scan::Spectra> dEdge;      // [sl * 3 + ec]
    std::vector<scan::Spectra> dTag;       // [reg * 3 + ec]
    scan::Spectra dS0;
    {
        scan::Data d;
        if (!scan::Load(d)) return;
        auto S0 = [](const scan::Row& r) { return PassSel(kSelS0, r, r.cor_ok == 1, r.dz_pi, r.mm4); };
        // ---- 1. rows of S0 to a small tree
        {
            TFile* fr = TFile::Open(Form("%s/wp6_rows.root", kOutDir), "RECREATE");
            TTree* t = new TTree("rows", "S0 rows (cor_ok == 1): k = 0 e+e-, 1 e+e+, 2 e-e-; file order (event groups adjacent)");
            Int_t k, grp; Float_t m4, mpipi, mee;
            t->Branch("k", &k, "k/I"); t->Branch("grp", &grp, "grp/I");
            t->Branch("m4", &m4, "m4/F"); t->Branch("mpipi", &mpipi, "mpipi/F"); t->Branch("mee", &mee, "mee/F");
            long nk[3] = {0, 0, 0};
            for (int s = 0; s < 3; ++s)
                for (const auto& r : d.s[s]) {
                    if (r.cor_ok != 1 || !S0(r)) continue;
                    k = s; grp = r.grp; m4 = r.m4; mpipi = r.mpipi; mee = r.mee;
                    t->Fill(); ++nk[s];
                }
            fr->cd(); t->Write(); fr->Close();
            P("A1. data S0 rows written: e+e- %ld, e+e+ %ld, e-e- %ld (wp6_rows.root)\n", nk[0], nk[1], nk[2]);
        }
        // ---- 2. edge spectra: slice x class (own CB of the slice; classes with the spread rule)
        {
            std::vector<scan::Sel> vs;
            for (int i = 0; i < kNSl; ++i)
                for (int e = 0; e < kNEc; ++e) {
                    const float a = kSlLo[i], b = kSlHi[i];
                    vs.push_back([=](const scan::Row& r) {
                        if (!S0(r) || !(r.mpipi >= a && r.mpipi < b)) return false;
                        if (e == 1) return r.mee < 0.010f;
                        if (e == 2) return r.mee >= 0.010f;
                        return true;
                    });
                }
            vs.push_back(S0);
            std::vector<scan::Spectra> sp = scan::FillMany(d, vs);
            for (int i = 0; i < kNSl; ++i)
                for (int e = 0; e < kNEc; ++e) {
                    const scan::Spectra& all = sp[i * kNEc];
                    dEdge.push_back(e == 0 ? all : SpreadCB(all, sp[i * kNEc + e]));
                }
            dS0 = sp.back();
        }
        // ---- 3. tag spectra (m5c): region x class
        {
            std::vector<scan::Sel> vt;
            for (int g = 0; g < kNTagReg; ++g)
                for (int e = 0; e < kNEc; ++e)
                    vt.push_back([=](const scan::Row& r) {
                        if (!TaggedW(r, kPi0Lo, kPi0Hi) || !S0(r)) return false;
                        if (!InTagReg(g, r.m4, r.mpipi)) return false;
                        if (e == 1) return r.mee < 0.010f;
                        if (e == 2) return r.mee >= 0.010f;
                        return true;
                    });
            scan::FillOpt o;
            o.nbin = kTagNbin; o.lo = kTagHLo; o.hi = kTagHHi; o.reb = kTagReb;
            o.x = [](const scan::Row& r) { return M5c(r); };
            dTag = scan::FillMany(d, vt, o);
        }
    }

    // =========================================================================================================
    // B. SIMULATION
    // =========================================================================================================
    const int nbE = scan::kNbin;                    // 240 x 2.5 MeV, 0.35-0.95
    const size_t nEdge = (size_t)kNSl * kNEc * kNEg * nbE;
    const size_t nTag = (size_t)kNTagReg * kNEc * 2 * kTagNbin;          // groups: 0 T12, 1 C23 (truth)
    const size_t nTagM = (size_t)2 * 2 * kTagNbin;                        // matched pi0 window: region 0/1 x T12/C23
    const int nbEeg = 300;                                                // M(e+e-gamma) 0-0.30, 1 MeV
    const size_t nWin = (size_t)kNWc * kNCell * nW;
    const int nbT = 350;                                                  // 1 MeV templates 0.40-0.75
    const size_t nTpl = (size_t)kNWc * kNCell * nbT;
    const size_t nM4 = (size_t)4 * nbE;                                    // P_main inputs: T12 r=1, T12 r=0.89, C45 r=1, C45 r=0.89
    Acc aEdge, aTag, aTagM, aEeg, aWin, aTpl, aM4;
    aEdge.Init(nEdge); aTag.Init(nTag); aTagM.Init(nTagM); aEeg.Init(2 * nbEeg); aWin.Init(nWin); aTpl.Init(nTpl); aM4.Init(nM4);
    RepAcc rEdge, rTag, rWin;
    rEdge.Init((size_t)kNSl * kNEc * 2 * nbE);        // groups T12, C23
    rTag.Init((size_t)2 * kTagNbin);                  // T12, regions 0 and 1, all M(ee)
    rWin.Init(nWin);
    double pi0Mean = 0, pi0Sig = 0, pi0MatchLo = 0, pi0MatchHi = 0;
    long nSimRows = 0, nSimS0 = 0;
    {
        simphys::Init();
        simphys::Install();
        simphys::LoadAngles();
        simeng::SimData* D = new simeng::SimData();
        if (!simeng::Load(*D, simeng::kCompact, simeng::kMaskUS, -1e9, 1e9, 1)) return;
        const int trees[3] = {simeng::kStdUS, simeng::kLepUS, simeng::kEtaUS};
        // ---- pass 1: M(e+e-gamma) of the bump classes with a good photon (pi0 position of the simulation)
        {
            std::vector<double> h(nbEeg, 0.);
            for (int t = 0; t < 3; ++t)
                for (const auto& r : D->t[trees[t]]) {
                    if (!(r.tmpl == 1 && (r.cls == 2 || r.cls == 3))) continue;
                    if (!(r.g_idx > 0 && r.mee < kMeeMax)) continue;
                    if (!PassSel(kSelS0, r, r.cor_ok == 1, 0.f, r.mm4) || r.ownA != 1) continue;
                    const double w = r.w * simphys::RowFactor(r, kPhysNom) * kSimToRaw;
                    const int b = (int)std::floor(r.m_eeg / 0.001);
                    if (b >= 0 && b < nbEeg) h[b] += w;
                }
            // Gaussian fit of the peak region 0.10-0.17 with a linear background (TH1 fit)
            TH1D* hh = H1("tmpEeg", "", nbEeg, 0., 0.30, h.data());
            TF1 fg("fg", "gaus(0) + pol1(3)", 0.105, 0.165);
            fg.SetParameters(hh->GetMaximum(), 0.135, 0.010, 0., 0.);
            hh->Fit(&fg, "QRN0");
            pi0Mean = fg.GetParameter(1); pi0Sig = std::fabs(fg.GetParameter(2));
            pi0MatchLo = std::floor((pi0Mean - 2 * pi0Sig) * 1000. + 0.5) / 1000.;
            pi0MatchHi = std::floor((pi0Mean + 2 * pi0Sig) * 1000. + 0.5) / 1000.;
            P("B1. simulated pi0 in M(e+e-gamma) (bump classes C2+C3, good photon, M(ee) < 0.10, S0): position %.4f, width %.1f MeV\n",
              pi0Mean, 1e3 * pi0Sig);
            P("    data (A1 notes, scan_bump.C): position 0.1459 +- 0.0005, width 13.25 +- 0.53 MeV (photon energy scale about 8 %% high).\n");
            P("    frozen data window %.3f-%.3f; matched simulation window (+-2 sigma of the simulation) %.3f-%.3f\n", kPi0Lo, kPi0Hi,
              pi0MatchLo, pi0MatchHi);
            delete hh;
        }
        // ---- pass 2: all spectra
        std::vector<double> kr(kNRep, 0.);
        TStopwatch sw;
        for (int t = 0; t < 3; ++t) {
            const std::vector<simeng::Row>& v = D->t[trees[t]];
            int cur = -1;
            bool haveK = false;
            auto flush = [&]() { aEdge.Flush(); aTag.Flush(); aTagM.Flush(); aEeg.Flush(); aWin.Flush(); aTpl.Flush(); aM4.Flush(); };
            for (size_t i = 0; i < v.size(); ++i) {
                const simeng::Row& r = v[i];
                if (r.grp != cur) { flush(); cur = r.grp; haveK = false; }
                ++nSimRows;
                if (r.ownA != 1) continue;
                if (!PassSel(kSelS0, r, r.cor_ok == 1, 0.f, r.mm4)) continue;
                int cn[6];
                const int ncn = CompsOf(r, false, cn);
                if (ncn == 0) continue;
                const double w = r.w * simphys::RowFactor(r, kPhysNom) * kSimToRaw;
                if (w == 0.) continue;
                ++nSimS0;
                bool isT12 = false, isC23 = false, isOeta = false, isBkg = false;
                int wc[6]; int nwc = 0;
                for (int q = 0; q < ncn; ++q) {
                    const int c = cn[q];
                    if (c == kT12) isT12 = true;
                    if (c == kC2 || c == kC3) isC23 = true;
                    if (c == kC4 || c == kC5 || c == kGeta) isOeta = true;
                    if (c == kGD || c == kGconv) isBkg = true;
                    const int x = WcOf(c);
                    if (x >= 0) wc[nwc++] = x;
                }
                if (!isT12 && nwc == 0) continue;   // C6 (X17) and T3/T4 rows: not used here
                auto needK = [&]() { if (!haveK) { GroupPois(kSeedSim, trees[t], r.grp, kNRep, kr.data()); haveK = true; } };
                const double m4p = M4sim(r.m4, r.m4_true, kS, kR);
                const double m4r = M4sim(r.m4, r.m4_true, kS, kRlo);
                const double mpp = kSPiPi * r.mpipi;
                const double mee = r.mee;
                const int slice = Slice(mpp), ecl = Eclass(mee);
                // (a) edge spectra
                if (isT12) {
                    const int ib = (m4p >= scan::kMlo && m4p < scan::kMhi) ? (int)((m4p - scan::kMlo) / scan::kBinW) : -1;
                    if (ib >= 0 && ib < nbE)
                        for (int s = 0; s < kNSl; ++s) {
                            if (!(mpp >= kSlLo[s] && mpp < kSlHi[s])) continue;
                            for (int e = 0; e < kNEc; ++e) {
                                if (e == 1 && !(mee < 0.010)) continue;
                                if (e == 2 && !(mee >= 0.010)) continue;
                                const bool gs[kNEg] = {true, isC23, isOeta, isBkg};
                                for (int g = 0; g < kNEg; ++g) {
                                    if (!gs[g]) continue;
                                    aEdge.Add((((size_t)s * kNEc + e) * kNEg + g) * nbE + ib, w);
                                    if (g < 2) { needK(); rEdge.Add((((size_t)s * kNEc + e) * 2 + g) * nbE + ib, w, kr.data()); }
                                }
                            }
                        }
                    // (f) M4 spectra for P_main (all slices, all M(ee))
                    if (ib >= 0 && ib < nbE) {
                        aM4.Add(0 * nbE + ib, w);
                        if (r.cls == 4 || r.cls == 5) if (r.tmpl == 1) aM4.Add(2 * nbE + ib, w);
                    }
                    const int ibr = (m4r >= scan::kMlo && m4r < scan::kMhi) ? (int)((m4r - scan::kMlo) / scan::kBinW) : -1;
                    if (ibr >= 0 && ibr < nbE) {
                        aM4.Add(1 * nbE + ibr, w);
                        if (r.tmpl == 1 && (r.cls == 4 || r.cls == 5)) aM4.Add(3 * nbE + ibr, w);
                    }
                    // (c) tag
                    const double m5c = (r.g_idx > 0 && r.mee < kMeeMax) ? M5c(r) : -1.;
                    const int it = (m5c >= kTagHLo && m5c < kTagHHi) ? (int)((m5c - kTagHLo) / ((kTagHHi - kTagHLo) / kTagNbin)) : -1;
                    if (it >= 0 && it < kTagNbin) {
                        if (TaggedW(r, kPi0Lo, kPi0Hi)) {
                            for (int g = 0; g < kNTagReg; ++g) {
                                if (!InTagReg(g, m4p, mpp)) continue;
                                for (int e = 0; e < kNEc; ++e) {
                                    if (e == 1 && !(mee < 0.010)) continue;
                                    if (e == 2 && !(mee >= 0.010)) continue;
                                    aTag.Add((((size_t)g * kNEc + e) * 2 + 0) * kTagNbin + it, w);
                                    if (isC23) aTag.Add((((size_t)g * kNEc + e) * 2 + 1) * kTagNbin + it, w);
                                }
                                if (g < 2) { needK(); rTag.Add((size_t)g * kTagNbin + it, w, kr.data()); }
                            }
                        }
                        if (TaggedW(r, pi0MatchLo, pi0MatchHi))
                            for (int g = 0; g < 2; ++g) {
                                if (!InTagReg(g, m4p, mpp)) continue;
                                aTagM.Add(((size_t)g * 2 + 0) * kTagNbin + it, w);
                                if (isC23) aTagM.Add(((size_t)g * 2 + 1) * kTagNbin + it, w);
                            }
                    }
                    // (d) M(e+e-gamma) of all rows with a good photon (T12 and C23)
                    if (r.g_idx > 0 && r.mee < kMeeMax) {
                        const int b = (int)std::floor(r.m_eeg / 0.001);
                        if (b >= 0 && b < nbEeg) { aEeg.Add(b, w); if (isC23) aEeg.Add(nbEeg + b, w); }
                    }
                }
                // (e) windows and 1 MeV templates per (slice, class) and component
                if (slice >= 0 && ecl >= 0 && nwc > 0) {
                    const int cell = 2 * slice + ecl;
                    for (int x = 0; x < nW; ++x) {
                        if (!(m4p >= gW[x].lo && m4p < gW[x].hi)) continue;
                        needK();
                        for (int q = 0; q < nwc; ++q) {
                            const size_t idx = ((size_t)wc[q] * kNCell + cell) * nW + x;
                            aWin.Add(idx, w);
                            rWin.Add(idx, w, kr.data());
                        }
                    }
                    if (m4p >= 0.40 && m4p < 0.75) {
                        const int b = std::min(nbT - 1, (int)((m4p - 0.40) / 0.001));
                        for (int q = 0; q < nwc; ++q) aTpl.Add(((size_t)wc[q] * kNCell + cell) * nbT + b, w);
                    }
                }
            }
            flush();
            P("B2. sim tree %s: %zu rows, %.0f s\n", simeng::TreeName(trees[t]), v.size(), sw.RealTime());
            sw.Continue();
        }
        delete D;
    }
    P("    simulated rows read %ld, used (owned, S0, with a component) %ld\n", nSimRows, nSimS0);

    // =========================================================================================================
    // C. FITS: edge step (data and simulation), tag, P_main
    // =========================================================================================================
    auto simEdge = [&](int s, int e, int g, const std::vector<double>* sub = nullptr, double subf = 0.) {
        std::vector<double> y(nbE);
        for (int i = 0; i < nbE; ++i) {
            y[i] = aEdge.sum[(((size_t)s * kNEc + e) * kNEg + g) * nbE + i];
            if (sub) y[i] -= subf * (*sub)[i];
        }
        return y;
    };
    auto simEdgeRep = [&](int s, int e, int g, int q) {
        std::vector<double> y(nbE);
        for (int i = 0; i < nbE; ++i) y[i] = rEdge.r[((((size_t)s * kNEc + e) * 2 + g) * nbE + i) * kNRep + q];
        return y;
    };
    // ---- C1. edge step
    P("\nC1. EDGE STEP E (rows per 2.5 MeV/c2 summed over M(pipi) 0.30-0.34, 0.34-0.37, 0.37-0.39, 0.39-0.40; scan_bump.C procedure:\n");
    P("    pol2 [+ eta Gaussian 0.5439 / 5.4 MeV] + h x step at the fixed edge 0.99277 M4max(M(pipi)/0.9970), smearing 5.4 MeV,\n");
    P("    range edge -+ 45 MeV (-+ 35 MeV from 0.39 on); data: N - CB, CB = own 2 sqrt(N++ N--) of the slice (classes: spread rule);\n");
    P("    simulation: T1 + T2 at SMASH x PDG BR, expected raw data rows, the SAME fit (Poisson weights of the expected rows))\n");
    std::vector<scan::Spectra> sEdge[kNEc];           // simulation spectra T12 per class [slice]
    std::vector<std::vector<double>> sEdgeY[kNEc];
    ERes Ed[kNEc], Es[kNEc];
    std::vector<scan::StepFit> keepD[kNEc], keepS[kNEc];
    double fEdge[kNEc], fEdgeMC[kNEc], fEdgeStat[kNEc], fEdgeSpread[kNEc], EsRms[kNEc], EsMean[kNEc];
    FILE* mr = gOut;
    for (int e = 0; e < kNEc; ++e) {
        std::vector<scan::Spectra> dd;
        for (int s = 0; s < kNSl; ++s) dd.push_back(dEdge[s * kNEc + e]);
        for (int s = 0; s < kNSl; ++s) { sEdgeY[e].push_back(simEdge(s, e, 0)); sEdge[e].push_back(SimSpectra(sEdgeY[e][s].data(), nbE, scan::kMlo, scan::kMhi, scan::kCBreb)); }
        Ed[e] = EdgeSum(dd, nullptr, &keepD[e]);
        Es[e] = EdgeSum(sEdge[e], nullptr, &keepS[e]);
        // replicas
        double s1 = 0, s2 = 0; int n = 0;
        for (int q = 0; q < kNRep; ++q) {
            std::vector<scan::Spectra> rr;
            std::vector<std::vector<double>> yy;
            for (int s = 0; s < kNSl; ++s) yy.push_back(simEdgeRep(s, e, 0, q));
            for (int s = 0; s < kNSl; ++s) rr.push_back(SimSpectra(yy[s].data(), nbE, scan::kMlo, scan::kMhi, scan::kCBreb));
            const ERes x = EdgeSum(rr);
            s1 += x.E; s2 += x.E * x.E; ++n;
        }
        EsMean[e] = s1 / n; EsRms[e] = std::sqrt(std::max(0., (s2 - n * EsMean[e] * EsMean[e]) / (n - 1)));
        // model variations applied to data and simulation together (order 2, 3; half range 35/45/55 MeV; smearing 4.4/5.4/6.4 MeV)
        double rmin = 1e9, rmax = -1e9, emin = 1e9, emax = -1e9; int kept = 0;
        for (int ord = 2; ord <= 3; ++ord)
            for (int ih = 0; ih < 3; ++ih)
                for (int is = -1; is <= 1; ++is) {
                    auto mod = [&](int, scan::StepCfg& c) { c.order = ord; c.halfRange = 0.035 + 0.010 * ih; c.s = scan::kEdgeSigma + 0.001 * is; };
                    const ERes a = EdgeSum(dd, mod), b = EdgeSum(sEdge[e], mod);
                    // per slice (scan_bump.C EdgeSpread rule): a varied data fit with p <= 0.05 keeps the nominal
                    // value of that slice, in the data AND in the simulation
                    double Ea = 0, Eb = 0; int nv = 0;
                    for (int s = 0; s < kNSigSl; ++s) {
                        const bool use = a.ok[s] && b.ok[s] && TMath::Prob(a.chi2[s], a.ndf[s]) > 0.05;
                        Ea += use ? a.h[s] : Ed[e].h[s];
                        Eb += use ? b.h[s] : Es[e].h[s];
                        if (use) ++nv;
                    }
                    if (nv == 0 || !(Eb > 0)) continue;
                    ++kept;
                    rmin = std::min(rmin, Ea / Eb); rmax = std::max(rmax, Ea / Eb);
                    emin = std::min(emin, Ea); emax = std::max(emax, Ea);
                }
        fEdge[e] = Ed[e].E / Es[e].E;
        fEdgeStat[e] = Ed[e].eE / Es[e].E;
        fEdgeMC[e] = fEdge[e] * EsRms[e] / Es[e].E;
        fEdgeSpread[e] = kept > 0 ? 0.5 * (std::max(rmax, fEdge[e]) - std::min(rmin, fEdge[e])) : 0.;
        P("\n  [%s]\n", kEcName[e]);
        P("    %-10s | %-34s | %-34s | %s\n", "M(pipi)", "data: h +- err (chi2/ndf)", "sim T1+T2: h +- Poisson err (chi2/ndf)", "data/sim");
        for (int s = 0; s < kNSl; ++s)
            P("    %.2f-%.2f%s | %8.1f +- %6.1f (%5.1f/%-3d)%s | %8.1f +- %6.1f (%7.1f/%-3d)    | %6.3f\n", kSlLo[s], kSlHi[s], s == 5 ? "*" : " ",
              Ed[e].h[s], Ed[e].eh[s], Ed[e].chi2[s], Ed[e].ndf[s], s >= kNSigSl ? " (not in E)" : "           ", Es[e].h[s], Es[e].eh[s],
              Es[e].chi2[s], Es[e].ndf[s], Es[e].h[s] != 0 ? Ed[e].h[s] / Es[e].h[s] : 0.);
        P("    E_data = %.1f +- %.1f (stat);  E_sim = %.1f +- %.1f (MC, 100 replicas; replica mean %.1f)\n", Ed[e].E, Ed[e].eE, Es[e].E, EsRms[e],
          EsMean[e]);
        P("    f_eta(edge) = E_data/E_sim = %.4f +- %.4f (stat) +- %.4f (MC) +- %.4f (model: half range of the ratio over %d joint variations; a varied slice fit with p <= 0.05 keeps its nominal h)\n",
          fEdge[e], fEdgeStat[e], fEdgeMC[e], fEdgeSpread[e], kept);
        P("    data E over the same variations: %.0f-%.0f\n", emin, emax);
        if (mr) fprintf(mr, "KEY edge %s E_data %.4f %.4f E_sim %.4f mc %.4f f %.5f stat %.5f mc %.5f model %.5f nvar %d\n", e == 0 ? "all" : e == 1 ? "E1" : "E2",
                        Ed[e].E, Ed[e].eE, Es[e].E, EsRms[e], fEdge[e], fEdgeStat[e], fEdgeMC[e], fEdgeSpread[e], kept);
        for (int s = 0; s < kNSl; ++s)
            if (mr) fprintf(mr, "KEY edgeslice %s %d %.2f %.2f data %.4f %.4f %.3f %d sim %.4f %.4f %.3f %d\n", e == 0 ? "all" : e == 1 ? "E1" : "E2", s,
                            kSlLo[s], kSlHi[s], Ed[e].h[s], Ed[e].eh[s], Ed[e].chi2[s], Ed[e].ndf[s], Es[e].h[s], Es[e].eh[s], Es[e].chi2[s], Es[e].ndf[s]);
    }
    // data-only model spread of E in the scan_bump way (per slice half range over 18 fits with p > 0.05, added linearly)
    {
        double sum = 0;
        for (int s = 0; s < kNSigSl; ++s) {
            const scan::StepFit nom = scan::FitEdgeStep(dEdge[s * kNEc], SliceCfg(s));
            double mn = nom.h, mx = nom.h;
            for (int ord = 2; ord <= 3; ++ord)
                for (int ih = 0; ih < 3; ++ih)
                    for (int is = -1; is <= 1; ++is) {
                        scan::StepCfg c = SliceCfg(s);
                        c.order = ord; c.halfRange = 0.035 + 0.010 * ih; c.s = scan::kEdgeSigma + 0.001 * is;
                        const scan::StepFit f = scan::FitEdgeStep(dEdge[s * kNEc], c);
                        if (!f.ok || !(TMath::Prob(f.chi2, f.ndf) > 0.05)) continue;
                        mn = std::min(mn, f.h); mx = std::max(mx, f.h);
                    }
            sum += 0.5 * (mx - mn);
        }
        P("\n    data E (all M(ee)) with the scan_bump model spread: %.0f +- %.0f (stat) +- %.0f (spread)  [data team, pre + vertex: 1142 +- 59 +- 236]\n",
          Ed[0].E, Ed[0].eE, sum);
        if (mr) fprintf(mr, "KEY edge_data_spread %.4f\n", sum);
    }
    // decomposition of the simulated step: without the bump classes, and without the non-eta background
    {
        P("\n    decomposition of the simulated step (all M(ee)): E(T12) = %.1f; E(T12 - C23) = ", Es[0].E);
        std::vector<scan::Spectra> a, b, c;
        std::vector<std::vector<double>> ya, yb, yc;
        for (int s = 0; s < kNSl; ++s) {
            std::vector<double> t = simEdge(s, 0, 0), c23 = simEdge(s, 0, 1), bk = simEdge(s, 0, 3), oe = simEdge(s, 0, 2);
            std::vector<double> u(nbE), v(nbE), w(nbE);
            for (int i = 0; i < nbE; ++i) { u[i] = t[i] - c23[i]; v[i] = t[i] - bk[i]; w[i] = c23[i]; }
            ya.push_back(u); yb.push_back(v); yc.push_back(w);
        }
        for (int s = 0; s < kNSl; ++s) { a.push_back(SimSpectra(ya[s].data(), nbE, scan::kMlo, scan::kMhi, 4)); b.push_back(SimSpectra(yb[s].data(), nbE, scan::kMlo, scan::kMhi, 4)); c.push_back(SimSpectra(yc[s].data(), nbE, scan::kMlo, scan::kMhi, 4)); }
        const ERes ra = EdgeSum(a), rb = EdgeSum(b), rc = EdgeSum(c);
        P("%.1f (non-bump part), E(T12 - G_D - G_C) = %.1f (eta part), E(C23 alone) = %.1f\n", ra.E, rb.E, rc.E);
        if (mr) fprintf(mr, "KEY edge_decomp T12 %.4f noC23 %.4f etaonly %.4f C23 %.4f\n", Es[0].E, ra.E, rb.E, rc.E);
    }

    // ---- C2. tag
    P("\nC2. ECAL TAG: m5c peak content T (P_main fit: frozen data shape %.4f / %.1f MeV, pol%d, %.2f-%.2f; free = position and width free)\n",
      kTagMean, 1e3 * kTagSigma, kTagOrder, kTagLo, kTagHi);
    P("    data: N - CB of the tagged rows (S0); simulation (emulated tag, UNVERIFIED): T1 + T2 tagged rows, the same fits; truth = C2 + C3 rows\n");
    auto simTag = [&](const Acc& A, size_t off) { std::vector<double> y(kTagNbin); for (int i = 0; i < kTagNbin; ++i) y[i] = A.sum[off + i]; return y; };
    double Td[kNTagReg][kNEc], eTd[kNTagReg][kNEc], Ts[kNTagReg][kNEc], Tsf[kNTagReg][kNEc], Ttr[kNTagReg][kNEc], TtrW[kNTagReg][kNEc];
    std::vector<scan::PeakFit> tfD, tfS;
    for (int g = 0; g < kNTagReg; ++g)
        for (int e = 0; e < kNEc; ++e) {
            const scan::Spectra& sd = dTag[g * kNEc + e];
            const scan::PeakFit fd = scan::FitPeakMain(sd, TagCfg());
            const scan::FreeFit ffd = scan::FitPeakFree(sd, TagCfg());
            const std::vector<double> ys = simTag(aTag, (((size_t)g * kNEc + e) * 2 + 0) * kTagNbin);
            const std::vector<double> yc = simTag(aTag, (((size_t)g * kNEc + e) * 2 + 1) * kTagNbin);
            const scan::Spectra ss = SimSpectra(ys.data(), kTagNbin, kTagHLo, kTagHHi, kTagReb);
            const scan::PeakFit fs = scan::FitPeakMain(ss, TagCfg());
            const scan::FreeFit ffs = scan::FitPeakFree(ss, TagCfg());
            double tr = 0, trw = 0;
            for (int i = 0; i < kTagNbin; ++i) {
                const double x = kTagHLo + (i + 0.5) * (kTagHHi - kTagHLo) / kTagNbin;
                tr += yc[i];
                if (x >= kTagLo && x < kTagHi) trw += yc[i];
            }
            Td[g][e] = fd.S; eTd[g][e] = fd.eS; Ts[g][e] = fs.S; Tsf[g][e] = ffs.ok ? ffs.S : 0.; Ttr[g][e] = tr; TtrW[g][e] = trw;
            tfD.push_back(fd); tfS.push_back(fs);
            P("  %-28s %-18s | data T = %7.1f +- %5.1f (chi2 %5.1f/%d; free: %.4f, %.1f MeV, T %.1f) | sim T = %7.1f (frozen), %7.1f (free: %.4f, %.1f MeV) | truth C2+C3 %7.1f (all m5c), %7.1f (0.49-0.62)\n",
              kTagRegName[g], kEcName[e], fd.S, fd.eS, fd.chi2, fd.ndf, ffd.mean, 1e3 * ffd.sigma, ffd.S, fs.S, ffs.S, ffs.mean, 1e3 * ffs.sigma, tr, trw);
            if (mr) fprintf(mr, "KEY tag %d %d data %.4f %.4f %.3f %d free %.5f %.5f %.4f sim %.4f free %.4f %.5f %.5f truth %.4f %.4f\n", g, e, fd.S, fd.eS, fd.chi2, fd.ndf,
                            ffd.mean, ffd.sigma, ffd.S, fs.S, ffs.S, ffs.mean, ffs.sigma, tr, trw);
        }
    // tag replicas (frozen fit and free fit, regions 0 and 1, all M(ee))
    double tagRms[2] = {0, 0}, tagRmsF[2] = {0, 0};
    for (int g = 0; g < 2; ++g) {
        double s1 = 0, s2 = 0, f1 = 0, f2 = 0; int n = 0, nf = 0;
        for (int q = 0; q < kNRep; ++q) {
            std::vector<double> y(kTagNbin);
            for (int i = 0; i < kTagNbin; ++i) y[i] = rTag.r[((size_t)g * kTagNbin + i) * kNRep + q];
            const scan::Spectra ss = SimSpectra(y.data(), kTagNbin, kTagHLo, kTagHHi, kTagReb);
            const scan::PeakFit fs = scan::FitPeakMain(ss, TagCfg());
            if (fs.ok) { s1 += fs.S; s2 += fs.S * fs.S; ++n; }
            const scan::FreeFit ff = scan::FitPeakFree(ss, TagCfg());
            if (ff.ok) { f1 += ff.S; f2 += ff.S * ff.S; ++nf; }
        }
        if (n > 1) tagRms[g] = std::sqrt(std::max(0., (s2 - s1 * s1 / n) / (n - 1)));
        if (nf > 1) tagRmsF[g] = std::sqrt(std::max(0., (f2 - f1 * f1 / nf) / (nf - 1)));
        P("    MC error of the simulated T (%s): %.1f (frozen fit), %.1f (free fit; %d of %d replicas converged)\n", kTagRegName[g], tagRms[g], tagRmsF[g], nf, kNRep);
        if (mr) fprintf(mr, "KEY tagmc %d %.4f %.4f %d\n", g, tagRms[g], tagRmsF[g], nf);
    }
    // matched pi0 window in the simulation
    for (int g = 0; g < 2; ++g) {
        const std::vector<double> ys = simTag(aTagM, ((size_t)g * 2 + 0) * kTagNbin), yc = simTag(aTagM, ((size_t)g * 2 + 1) * kTagNbin);
        const scan::Spectra ss = SimSpectra(ys.data(), kTagNbin, kTagHLo, kTagHHi, kTagReb);
        const scan::PeakFit fs = scan::FitPeakMain(ss, TagCfg());
        const scan::FreeFit ff = scan::FitPeakFree(ss, TagCfg());
        double tr = 0; for (double x : yc) tr += x;
        P("    matched window %.3f-%.3f, %s: sim T = %.1f (frozen), %.1f (free); truth C2+C3 %.1f\n", pi0MatchLo, pi0MatchHi, kTagRegName[g], fs.S, ff.S, tr);
        if (mr) fprintf(mr, "KEY tagmatched %d %.4f %.4f %.4f %.4f %.4f\n", g, pi0MatchLo, pi0MatchHi, fs.S, ff.S, tr);
    }
    if (mr) fprintf(mr, "KEY pi0sim %.5f %.5f\n", pi0Mean, pi0Sig);

    // ---- C3. P_main
    P("\nC3. P_main on S0 (N - CB, frozen Gaussian %.4f / %.1f MeV + pol%d, %.4f-%.4f; spread = 48 model fits, kept p > 0.05)\n", scan::kPeakMean,
      1e3 * scan::kPeakSigma, scan::kFitOrder, scan::kFitLo, scan::kFitHi);
    const scan::PeakFit pm = scan::FitPeakMain(dS0);
    const scan::Spread psp = scan::PeakModelSpread(dS0);
    const scan::FreeFit pmf = scan::FitPeakFree(dS0);
    P("    data: S = %.1f +- %.1f (stat) +- %.1f (spread, %d of %d fits kept); chi2/ndf %.1f/%d; window N %.0f, CB %.1f, S_win %.1f, S/B %.4f\n",
      pm.S, pm.eS, psp.half, psp.npass, psp.nfit, pm.chi2, pm.ndf, pm.N, pm.CB, pm.Swin, pm.SoB);
    P("          free shape: position %.5f +- %.5f, width %.2f +- %.2f MeV, S %.1f +- %.1f\n", pmf.mean, pmf.eMean, 1e3 * pmf.sigma, 1e3 * pmf.eSigma, pmf.S, pmf.eS);
    P("          [A1 / data team: 4388 +- 307 +- 1078]\n");
    if (mr) fprintf(mr, "KEY pmain data %.4f %.4f %.4f %d %d %.4f %d Swin %.4f free %.6f %.6f %.6f %.6f %.4f %.4f\n", pm.S, pm.eS, psp.half, psp.npass, psp.nfit, pm.chi2, pm.ndf,
                    pm.Swin, pmf.mean, pmf.eMean, pmf.sigma, pmf.eSigma, pmf.S, pmf.eS);
    for (int v = 0; v < 2; ++v) {
        std::vector<double> y(nbE), yc(nbE);
        for (int i = 0; i < nbE; ++i) { y[i] = aM4.sum[(size_t)v * nbE + i]; yc[i] = aM4.sum[(size_t)(2 + v) * nbE + i]; }
        const scan::Spectra ss = SimSpectra(y.data(), nbE, scan::kMlo, scan::kMhi, 4);
        const scan::PeakFit fs = scan::FitPeakMain(ss);
        const scan::FreeFit ff = scan::FitPeakFree(ss);
        double c45 = 0, c45w = 0;
        for (int i = 0; i < nbE; ++i) { c45 += yc[i]; c45w += yc[i] * ss.Frac(i, scan::kPeakLo, scan::kPeakHi); }
        P("    sim T1+T2 (r = %.2f): S = %.1f (frozen shape), free: position %.5f, width %.2f MeV, S %.1f; truth C4+C5: %.1f in P, %.1f in 0.35-0.95\n",
          v == 0 ? 1.0 : kRlo, fs.S, ff.mean, 1e3 * ff.sigma, ff.S, c45w, c45);
        if (mr) fprintf(mr, "KEY pmain sim%d %.4f free %.6f %.6f %.4f truthP %.4f truth %.4f\n", v, fs.S, ff.mean, ff.sigma, ff.S, c45w, c45);
    }

    // =========================================================================================================
    // D. OUTPUT
    // =========================================================================================================
    {
        // windows (flat index ((comp * 10 + cell) * nW + win)) + variance + replicas
        TH1D* hw = new TH1D("win_sum", "sim window sums [raw rows]; index ((comp*10 + cell)*nW + win)", (int)nWin, 0, (double)nWin);
        for (size_t i = 0; i < nWin; ++i) { hw->SetBinContent((int)i + 1, aWin.sum[i]); hw->SetBinError((int)i + 1, std::sqrt(aWin.var[i])); }
        Write(dSim, hw);
        TH2D* hr = new TH2D("win_rep", "sim window sums, replicas", (int)nWin, 0, (double)nWin, kNRep, -0.5, kNRep - 0.5);
        for (size_t i = 0; i < nWin; ++i) for (int q = 0; q < kNRep; ++q) hr->SetBinContent((int)i + 1, q + 1, rWin.r[i * kNRep + q]);
        Write(dSim, hr);
        std::string wl;
        for (int x = 0; x < nW; ++x) wl += Form("%s %.6f %.6f;", gW[x].tag.c_str(), gW[x].lo, gW[x].hi);
        TNamed* nm = new TNamed("windows", wl.c_str()); dSim->WriteTObject(nm, "windows", "Overwrite"); delete nm;
        std::string cl; for (int c = 0; c < kNWc; ++c) cl += std::string(kWcTag[c]) + " ";
        TNamed* nc = new TNamed("comps", cl.c_str()); dSim->WriteTObject(nc, "comps", "Overwrite"); delete nc;
        // 1 MeV templates
        for (int c = 0; c < kNWc; ++c)
            for (int cell = 0; cell < kNCell; ++cell) {
                const size_t o = ((size_t)c * kNCell + cell) * nbT;
                Write(dSim, H1(Form("tpl_%s_S%dE%d", kWcTag[c], cell / 2 + 1, cell % 2 + 1), "raw rows per 1 MeV", nbT, 0.40, 0.75, &aTpl.sum[o], &aTpl.var[o]));
            }
        // edge spectra (sim, groups) and data spectra
        for (int s = 0; s < kNSl; ++s)
            for (int e = 0; e < kNEc; ++e) {
                for (int g = 0; g < kNEg; ++g) {
                    const size_t o = (((size_t)s * kNEc + e) * kNEg + g) * nbE;
                    Write(dSim, H1(Form("edge_%s_sl%d_e%d", kEgName[g], s, e), "raw rows per 2.5 MeV", nbE, scan::kMlo, scan::kMhi, &aEdge.sum[o], &aEdge.var[o]));
                }
                const scan::Spectra& sd = dEdge[s * kNEc + e];
                Write(dData, H1(Form("edge_N_sl%d_e%d", s, e), "N", nbE, scan::kMlo, scan::kMhi, sd.n[0].data(), sd.n[0].data()));
                Write(dData, H1(Form("edge_CB_sl%d_e%d", s, e), "CB", nbE, scan::kMlo, scan::kMhi, sd.cb.data(), sd.cbVar.data()));
                Write(dData, H1(Form("edge_sig_sl%d_e%d", s, e), "N - CB", nbE, scan::kMlo, scan::kMhi, sd.sig.data(), sd.sigVar.data()));
                if (s < (int)keepD[e].size() && keepD[e][s].ok) {
                    Write(dData, Curve(keepD[e][s].tx, keepD[e][s].fy, Form("edgefit_sl%d_e%d", s, e)));
                    Write(dData, Curve(keepD[e][s].tx, keepD[e][s].fy0, Form("edgefit0_sl%d_e%d", s, e)));
                }
                if (s < (int)keepS[e].size() && keepS[e][s].ok) {
                    Write(dSim, Curve(keepS[e][s].tx, keepS[e][s].fy, Form("edgefit_sl%d_e%d", s, e)));
                    Write(dSim, Curve(keepS[e][s].tx, keepS[e][s].fy0, Form("edgefit0_sl%d_e%d", s, e)));
                }
            }
        // tag spectra
        const double tbw = (kTagHHi - kTagHLo) / kTagNbin;
        for (int g = 0; g < kNTagReg; ++g)
            for (int e = 0; e < kNEc; ++e) {
                const scan::Spectra& sd = dTag[g * kNEc + e];
                Write(dData, H1(Form("tag_sig_r%d_e%d", g, e), "N - CB", kTagNbin, kTagHLo, kTagHHi, sd.sig.data(), sd.sigVar.data()));
                Write(dData, H1(Form("tag_N_r%d_e%d", g, e), "N", kTagNbin, kTagHLo, kTagHHi, sd.n[0].data(), sd.n[0].data()));
                Write(dData, H1(Form("tag_CB_r%d_e%d", g, e), "CB", kTagNbin, kTagHLo, kTagHHi, sd.cb.data(), sd.cbVar.data()));
                for (int k = 0; k < 2; ++k) {
                    const size_t o = (((size_t)g * kNEc + e) * 2 + k) * kTagNbin;
                    Write(dSim, H1(Form("tag_%s_r%d_e%d", k == 0 ? "T12" : "C23", g, e), "raw rows per 2.5 MeV", kTagNbin, kTagHLo, kTagHHi, &aTag.sum[o], &aTag.var[o]));
                }
                const scan::PeakFit& fd = tfD[g * kNEc + e];
                const scan::PeakFit& fs = tfS[g * kNEc + e];
                if (fd.ok) { Write(dData, PeakCurve(fd, tbw, Form("tagfit_r%d_e%d", g, e), false)); Write(dData, PeakCurve(fd, tbw, Form("tagfit0_r%d_e%d", g, e), true)); }
                if (fs.ok) { Write(dSim, PeakCurve(fs, tbw, Form("tagfit_r%d_e%d", g, e), false)); Write(dSim, PeakCurve(fs, tbw, Form("tagfit0_r%d_e%d", g, e), true)); }
            }
        Write(dSim, H1("meeg_T12", "raw rows per 1 MeV", nbEeg, 0., 0.30, &aEeg.sum[0], &aEeg.var[0]));
        Write(dSim, H1("meeg_C23", "raw rows per 1 MeV", nbEeg, 0., 0.30, &aEeg.sum[nbEeg], &aEeg.var[nbEeg]));
        // P_main inputs
        const char* m4n[4] = {"m4_T12_r1", "m4_T12_r089", "m4_C45_r1", "m4_C45_r089"};
        for (int v = 0; v < 4; ++v) Write(dSim, H1(m4n[v], "raw rows per 2.5 MeV", nbE, scan::kMlo, scan::kMhi, &aM4.sum[(size_t)v * nbE], &aM4.var[(size_t)v * nbE]));
        Write(dData, H1("m4_N", "N", nbE, scan::kMlo, scan::kMhi, dS0.n[0].data(), dS0.n[0].data()));
        Write(dData, H1("m4_CB", "CB", nbE, scan::kMlo, scan::kMhi, dS0.cb.data(), dS0.cbVar.data()));
        Write(dData, H1("m4_sig", "N - CB", nbE, scan::kMlo, scan::kMhi, dS0.sig.data(), dS0.sigVar.data()));
        Write(dData, PeakCurve(pm, scan::kBinW, "pmain_fit", false));
        Write(dData, PeakCurve(pm, scan::kBinW, "pmain_fit0", true));
    }
    fo->Close();
    P("\nwritten %s/wp6_fill.root, wp6_rows.root, wp6_fill.txt; total %.0f s\n", kOutDir, total.RealTime());
    fclose(gOut); gOut = nullptr;
}
