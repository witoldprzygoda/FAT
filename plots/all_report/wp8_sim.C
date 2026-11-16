// ============================================================================
// wp8_sim.C  --  WP8 (physics quantities): the parts that need the ROWS
// ============================================================================
// SPEC report_all, WP8. Written by WP8 (prefix wp8), 2026-11-14.
//
// One pass over the simulated rows (SimEngine copy, mode A, owned rows, std + lep + eta, e+e- trees) and one over
// the data rows (scan_v2 compact). Nothing here changes a frozen definition: selections, windows, weights, mass
// transfer and the ECAL tag are those of defs.h / wp0_fill.C / wp6_fill.C (copied helpers are marked).
//
//   1. Table P1 inputs: simulated visible cross sections (sum w / L_sim) of the classes C2-C6 for
//      selection S0 and pre, physics weights V3 (code 108) and off (code 0), windows all M4 / B / P / X17c-M4,
//      M(ee) classes all / E1 < 10 / E2 >= 10 / >= 20 / 15.5-24.5 MeV (reconstructed) and >= 20 / >= 10 MeV (true M(ee));
//      event-group variances. Mass transfer M4' = 0.99615 [m4_true + (m4 - m4_true)], M(pipi)' = 0.9952 M(pipi).
//   2. eps2 QED check, part 2: the ECAL tag (frozen definition of scan_bump.C / wp6_fill.C) in M(ee) classes
//      < 8, 8-20, 20-30, 30-100 MeV (+ < 20, >= 20, < 10, >= 10): data m5c fits (N - CB, frozen 0.5463 / 6.2 MeV, pol2,
//      0.49-0.62) per class, the emulated simulated tag (T1 + T2 tagged rows, the same fit) and the truth sums of
//      C2 and C3. Regions: S0 all M4; pre all M4; and the bump start selection of the data scan (scan_extra.C:
//      base 'none', oa_lab < 40, 0.30 < M(pipi) < 0.41, M(e+e-gamma) 0.110-0.170) to reproduce the quoted 0.077 / 0.144.
//   3. (pT4, y4) of the peak: data S0 rows (e+e-, e+e+, e-e-) in P and in the side bands, and the simulated C4, C5 in P
//      and C2, C3 in B (V3, S0), on the grid of the generator histogram (pT 20 MeV/c, y_cm 0.05). y4 of the simulation is
//      recomputed at T = 4500 MeV (wp0_fill.C Kin4), the data keep 4530 MeV (stored y4).
//
// Outputs: /hdd3/przygoda/HADES/FAT_out/all_report/wp8/{wp8_sim.root, wp8_sim.txt}
// Run (from /home/przygoda/HADES/FAT/FAT; about 9-10 GB, 3-4 min):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_wp8", kTRUE)' \
//        'plots/all_report/wp8_sim.C+O'
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TStopwatch.h"
#include "TSystem.h"

#include "hists_io.h"
#include "sim_eta_physics.h"

namespace wp8s {
using namespace alldefs;

const char* const kOutDir = "/hdd3/przygoda/HADES/FAT_out/all_report/wp8";
FILE* gOut = nullptr;
#define P(...) do { printf(__VA_ARGS__); if (gOut) fprintf(gOut, __VA_ARGS__); } while (0)

// ---- accumulators with event-group variance (copy of wp6_fill.C) ------------------------------------------------
struct Acc {
    std::vector<double> sum, var, tmp;
    std::vector<int> touched;
    void Init(size_t n) { sum.assign(n, 0.); var.assign(n, 0.); tmp.assign(n, 0.); touched.clear(); }
    inline void Add(size_t i, double w) {
        if (tmp[i] == 0.) touched.push_back((int)i);
        tmp[i] += w; sum[i] += w;
        if (tmp[i] == 0.) tmp[i] = 1e-300;
    }
    void Flush() {
        for (int i : touched) { var[i] += tmp[i] * tmp[i]; tmp[i] = 0.; }
        touched.clear();
    }
};

// ---- ECAL tag (frozen, copy of wp6_fill.C) --------------------------------------------------------------------
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
scan::Spectra SimSpectra(const double* y, int nbin, double lo, double hi, int reb) {   // copy of wp6_fill.C
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
// tag M(ee) classes
const int kNTc = 9;
const char* const kTcName[kNTc] = {"all", "<8", "8-20", "20-30", "30-100", "<20", ">=20", "<10", ">=10"};
inline bool InTc(int k, double mee) {
    switch (k) {
        case 0: return true;
        case 1: return mee < 0.008;
        case 2: return mee >= 0.008 && mee < 0.020;
        case 3: return mee >= 0.020 && mee < 0.030;
        case 4: return mee >= 0.030;              // the tag itself requires M(ee) < 0.100
        case 5: return mee < 0.020;
        case 6: return mee >= 0.020;
        case 7: return mee < 0.010;
        case 8: return mee >= 0.010;
    }
    return false;
}
// tag regions: 0 S0 all M4; 1 pre all M4; 2 bump start selection of scan_extra.C (base 'none', oa_lab < 40,
// 0.30 < M(pipi) < 0.41, pi0 window 0.110-0.170)
const int kNTr = 3;
const char* const kTrName[kNTr] = {"S0, all M4 (frozen tag)", "pre, all M4 (frozen tag)", "A1 start selection (none, oa_lab<40, 0.30<M(pipi)<0.41, 0.110-0.170)"};
const double kTrLo[kNTr] = {kPi0Lo, kPi0Lo, 0.110}, kTrHi[kNTr] = {kPi0Hi, kPi0Hi, 0.170};

// ---- beam-dependent kinematics of a simulated row (copy of wp0_fill.C Kin4) ------------------------------------
inline void Kin4(const simeng::Row& r, double T, double& mm4, double& y4) {
    const double mp = simeng::kMpMeV;
    const double E = T + mp, pz = std::sqrt(E * E - mp * mp);
    TLorentzVector ini(0, 0, pz, E + mp);
    const TLorentzVector f = simeng::P4(r.ep_p, r.ep_theta, r.ep_phi, simeng::kMeMeV) + simeng::P4(r.em_p, r.em_theta, r.em_phi, simeng::kMeMeV) +
                             simeng::P4(r.pip_p, r.pip_theta, r.pip_phi, simeng::kMpiMeV) +
                             simeng::P4(r.pim_p, r.pim_theta, r.pim_phi, simeng::kMpiMeV);
    const TLorentzVector miss = ini - f;
    const double m2 = miss.M2();
    mm4 = (m2 >= 0 ? std::sqrt(m2) : -std::sqrt(-m2)) * 1e-3;
    y4 = f.Rapidity() - ini.Rapidity();
}

// ---- efficiency layout ------------------------------------------------------------------------------------------
const int kNEs = 2;   // selection: 0 S0, 1 pre
const int kNEw = 2;   // weights: 0 V3 (108), 1 off (0)
const int kNEc = 5;   // C2..C6
const int kNEwin = 4; // 0 all M4, 1 B, 2 P, 3 X17c M4 range 0.5312-0.5566
const int kNEm = 7;   // 0 all, 1 E1 < 10, 2 E2 >= 10, 3 >= 20, 4 15.5-24.5 (reco), 5 true >= 20, 6 true >= 10
const char* const kEsName[kNEs] = {"S0", "pre"};
const char* const kEwName[kNEw] = {"V3", "off"};
const char* const kEcName[kNEc] = {"C2", "C3", "C4", "C5", "C6"};
const char* const kEwinName[kNEwin] = {"all", "B", "P", "X17cM4"};
const char* const kEmName[kNEm] = {"all", "E1", "E2", "ge20", "x17", "tge20", "tge10"};
inline size_t EIdx(int s, int w, int c, int win, int m) { return ((((size_t)s * kNEw + w) * kNEc + c) * kNEwin + win) * kNEm + m; }
inline bool InEwin(int k, double m4) {
    switch (k) {
        case 0: return true;
        case 1: return m4 >= kBlo && m4 < kBhi;
        case 2: return m4 >= kPlo && m4 < kPhi;
        case 3: return m4 >= 0.5312 && m4 < 0.5566;
    }
    return false;
}
inline bool InEm(int k, double mee, double meet) {
    switch (k) {
        case 0: return true;
        case 1: return mee < 0.010;
        case 2: return mee >= 0.010;
        case 3: return mee >= 0.020;
        case 4: return mee >= 0.0155 && mee < 0.0245;
        case 5: return meet >= 0.020;
        case 6: return meet >= 0.010;
    }
    return false;
}

// ---- (pT4, y4) grid (= the generator histogram of the SMASH standard set) ---------------------------------------
const int kNPt = 75;  const double kPtLo = 0., kPtHi = 1500.;   // MeV/c, 20 MeV/c
const int kNY = 100;  const double kYLo = -2.5, kYHi = 2.5;     // y_cm, 0.05
inline int PtBin(double pt) { return (pt >= kPtLo && pt < kPtHi) ? (int)((pt - kPtLo) / 20.) : (pt >= kPtHi ? kNPt - 1 : -1); }
inline int YBin(double y) { if (!(y == y)) return -1; int b = (int)std::floor((y - kYLo) / 0.05); return std::max(0, std::min(kNY - 1, b)); }
const int kNKg = 6;   // sim groups: 0 C4 in P, 1 C5 in P, 2 C2 in B, 3 C3 in B, 4 G_eta in P, 5 C4+C5 in the side bands
const char* const kKgName[kNKg] = {"C4_P", "C5_P", "C2_B", "C3_B", "Geta_P", "C45_SB"};

void WriteH2(TFile* f, const char* name, const std::vector<double>& s, const std::vector<double>* v, size_t off) {
    TH2D* h = new TH2D(name, ";p_{T}(4) [MeV/c];y_{cm}(4)", kNPt, kPtLo, kPtHi, kNY, kYLo, kYHi);
    h->SetDirectory(nullptr);
    for (int i = 0; i < kNPt; ++i)
        for (int j = 0; j < kNY; ++j) {
            const size_t k = off + (size_t)i * kNY + j;
            h->SetBinContent(i + 1, j + 1, s[k]);
            if (v) h->SetBinError(i + 1, j + 1, std::sqrt(std::max(0., (*v)[k])));
        }
    f->WriteTObject(h, name, "Overwrite");
    delete h;
}

}  // namespace wp8s

void wp8_sim() {
    using namespace alldefs;
    using namespace wp8s;
    TStopwatch total;
    gSystem->mkdir(kOutDir, kTRUE);
    gOut = fopen(Form("%s/wp8_sim.txt", kOutDir), "w");
    P("# wp8_sim.C (WP8 row-level part), 2026-11-14. isBest >= 0 nominal. Simulation in sum w / L_sim [nb] (L_sim = %.2f nb^-1),\n", kLsim);
    P("# data in raw rows. sim: %s + %s (mode A, owned rows); M4' = %.5f [m4_true + (m4 - m4_true)], M(pipi)' = %.4f M(pipi)\n",
      simeng::kCompact, simeng::kCompactEta, kS, kSPiPi);
    P("# data: %s\n\n", scan::kCompact);
    TFile* fo = TFile::Open(Form("%s/wp8_sim.root", kOutDir), "RECREATE");

    // =========================================================================================================
    // A. SIMULATION
    // =========================================================================================================
    const size_t nEff = (size_t)kNEs * kNEw * kNEc * kNEwin * kNEm;
    const size_t nTagS = (size_t)kNTr * kNTc * kTagNbin;            // T12 m5c spectra (V3)
    const size_t nTagT = (size_t)kNTr * kNTc * 2 * 2;               // truth: [reg][class][C2/C3][all m5c / 0.49-0.62]
    const size_t nGrid = (size_t)kNPt * kNY;
    Acc aEff, aTagS, aTagT, aK;
    aEff.Init(nEff); aTagS.Init(nTagS); aTagT.Init(nTagT); aK.Init((size_t)kNKg * nGrid);
    long nRead = 0, nUsed = 0;
    {
        simphys::Init();
        simphys::Install();
        simphys::LoadAngles();
        simeng::SimData* D = new simeng::SimData();
        if (!simeng::Load(*D, simeng::kCompact, simeng::kMaskUS, -1e9, 1e9, 1)) return;
        const int trees[3] = {simeng::kStdUS, simeng::kLepUS, simeng::kEtaUS};
        TStopwatch sw;
        for (int t = 0; t < 3; ++t) {
            const std::vector<simeng::Row>& v = D->t[trees[t]];
            int cur = -1;
            auto flush = [&]() { aEff.Flush(); aTagS.Flush(); aTagT.Flush(); aK.Flush(); };
            for (size_t i = 0; i < v.size(); ++i) {
                const simeng::Row& r = v[i];
                if (r.grp != cur) { flush(); cur = r.grp; }
                ++nRead;
                if (r.ownA != 1) continue;
                if (!PassBase(r, r.cor_ok == 1)) continue;
                const bool pS0 = PassSel(kSelS0, r, r.cor_ok == 1, 0.f, r.mm4);
                const bool pPre = PassSel(kSelPre, r, r.cor_ok == 1, 0.f, r.mm4);
                const bool pNone = true;   // PassBase
                int cn[6];
                const int ncn = CompsOf(r, false, cn);
                if (ncn == 0) continue;
                const double wV3 = r.w * simphys::RowFactor(r, kPhysNom);
                const double wOff = r.w;
                const double m4p = M4sim(r.m4, r.m4_true, kS, kR);
                const double mpp = kSPiPi * r.mpipi;
                const double mee = r.mee, meet = r.mee_true;
                bool isT12 = false, isC2 = false, isC3 = false;
                for (int q = 0; q < ncn; ++q) {
                    if (cn[q] == kT12) isT12 = true;
                    if (cn[q] == kC2) isC2 = true;
                    if (cn[q] == kC3) isC3 = true;
                }
                ++nUsed;
                // (1) efficiencies
                for (int q = 0; q < ncn; ++q) {
                    const int c = cn[q];
                    if (c < kC2 || c > kC6) continue;
                    const int ec = c - kC2;
                    for (int s = 0; s < kNEs; ++s) {
                        if (!(s == 0 ? pS0 : pPre)) continue;
                        for (int w = 0; w < kNEw; ++w) {
                            const double ww = (c == kC6) ? 2.0 : (w == 0 ? wV3 : wOff);
                            if (ww == 0.) continue;
                            for (int win = 0; win < kNEwin; ++win) {
                                if (!InEwin(win, m4p)) continue;
                                for (int m = 0; m < kNEm; ++m)
                                    if (InEm(m, mee, meet)) aEff.Add(EIdx(s, w, ec, win, m), ww);
                            }
                        }
                    }
                }
                // (2) tag (T12 spectra, truth C2 / C3), V3
                if (isT12 && r.g_idx > 0 && mee < kMeeMax) {
                    const double m5c = M5c(r);
                    const int it = (m5c >= kTagHLo && m5c < kTagHHi) ? (int)((m5c - kTagHLo) / ((kTagHHi - kTagHLo) / kTagNbin)) : -1;
                    for (int g = 0; g < kNTr; ++g) {
                        bool sel = false;
                        if (g == 0) sel = pS0;
                        if (g == 1) sel = pPre;
                        if (g == 2) sel = pNone && r.oa_lab < 40.f && mpp > 0.30 && mpp < 0.41;
                        if (!sel || !TaggedW(r, kTrLo[g], kTrHi[g])) continue;
                        for (int k = 0; k < kNTc; ++k) {
                            if (!InTc(k, mee)) continue;
                            if (it >= 0 && it < kTagNbin) aTagS.Add(((size_t)g * kNTc + k) * kTagNbin + it, wV3);
                            if (isC2 || isC3) {
                                const int cc = isC2 ? 0 : 1;
                                aTagT.Add((((size_t)g * kNTc + k) * 2 + cc) * 2 + 0, wV3);
                                if (m5c >= kTagLo && m5c < kTagHi) aTagT.Add((((size_t)g * kNTc + k) * 2 + cc) * 2 + 1, wV3);
                            }
                        }
                    }
                }
                // (3) (pT4, y4) grid, S0, V3
                if (pS0 && wV3 != 0.) {
                    int kg[3]; int nkg = 0;
                    const bool inP = m4p >= kPlo && m4p < kPhi, inB = m4p >= kBlo && m4p < kBhi;
                    const bool inSB = (m4p >= 0.5250 && m4p < 0.5336) || (m4p >= 0.5542 && m4p < 0.5628);
                    for (int q = 0; q < ncn; ++q) {
                        if (cn[q] == kC4 && inP) kg[nkg++] = 0;
                        if (cn[q] == kC5 && inP) kg[nkg++] = 1;
                        if (cn[q] == kC2 && inB) kg[nkg++] = 2;
                        if (cn[q] == kC3 && inB) kg[nkg++] = 3;
                        if (cn[q] == kGeta && inP) kg[nkg++] = 4;
                        if ((cn[q] == kC4 || cn[q] == kC5) && inSB) kg[nkg++] = 5;
                    }
                    if (nkg > 0) {
                        double mm4s, y4s;
                        Kin4(r, kTsim, mm4s, y4s);
                        const int ip = PtBin(r.pt4), iy = YBin(y4s);
                        if (ip >= 0 && iy >= 0)
                            for (int q = 0; q < nkg; ++q) aK.Add((size_t)kg[q] * nGrid + (size_t)ip * kNY + iy, wV3);
                    }
                }
            }
            flush();
            P("A1. sim tree %s: %zu rows, %.0f s\n", simeng::TreeName(trees[t]), v.size(), sw.RealTime());
            sw.Continue();
        }
        delete D;
    }
    P("    simulated rows read %ld, used (owned, base, with a component) %ld\n\n", nRead, nUsed);

    // ---- efficiency sums (nb)
    P("A2. SIMULATED VISIBLE CROSS SECTIONS sum w / L_sim [nb] +- event-group error (KEY eff lines)\n");
    for (int s = 0; s < kNEs; ++s)
        for (int w = 0; w < kNEw; ++w)
            for (int c = 0; c < kNEc; ++c)
                for (int win = 0; win < kNEwin; ++win) {
                    std::string line = Form("KEY eff %s %s %s %s", kEsName[s], kEwName[w], kEcName[c], kEwinName[win]);
                    for (int m = 0; m < kNEm; ++m) {
                        const size_t k = EIdx(s, w, c, win, m);
                        line += Form(" %s %.8g %.4g", kEmName[m], aEff.sum[k] / kLsim, std::sqrt(aEff.var[k]) / kLsim);
                    }
                    P("%s\n", line.c_str());
                }

    // =========================================================================================================
    // B. DATA
    // =========================================================================================================
    std::vector<scan::Spectra> dTag;   // [reg * kNTc + class]
    std::vector<double> dK((size_t)2 * 3 * nGrid, 0.);   // [win P/SB][sample us/pp/mm][grid]
    double dWin[2][3] = {{0, 0, 0}, {0, 0, 0}};
    long dNoBin = 0;
    {
        scan::Data d;
        if (!scan::Load(d)) return;
        auto S0 = [](const scan::Row& r) { return PassSel(kSelS0, r, r.cor_ok == 1, r.dz_pi, r.mm4); };
        auto PRE = [](const scan::Row& r) { return PassSel(kSelPre, r, r.cor_ok == 1, r.dz_pi, r.mm4); };
        auto NONE = [](const scan::Row& r) { return PassSel(kSelNone, r, r.cor_ok == 1, r.dz_pi, r.mm4); };
        // ---- tag spectra
        std::vector<scan::Sel> vt;
        for (int g = 0; g < kNTr; ++g)
            for (int k = 0; k < kNTc; ++k) {
                const double lo = kTrLo[g], hi = kTrHi[g];
                vt.push_back([=](const scan::Row& r) {
                    if (!TaggedW(r, lo, hi)) return false;
                    if (g == 0 && !S0(r)) return false;
                    if (g == 1 && !PRE(r)) return false;
                    if (g == 2 && !(NONE(r) && r.oa_lab < 40.f && r.mpipi > 0.30f && r.mpipi < 0.41f)) return false;
                    return InTc(k, r.mee);
                });
            }
        scan::FillOpt o;
        o.nbin = kTagNbin; o.lo = kTagHLo; o.hi = kTagHHi; o.reb = kTagReb;
        o.x = [](const scan::Row& r) { return M5c(r); };
        dTag = scan::FillMany(d, vt, o);
        // ---- (pT4, y4) of the S0 rows in P and in the side bands
        for (int smp = 0; smp < 3; ++smp)
            for (const auto& r : d.s[smp]) {
                if (r.cor_ok != 1 || !S0(r)) continue;
                const double m4 = r.m4;
                int w = -1;
                if (m4 >= kPlo && m4 < kPhi) w = 0;
                else if ((m4 >= 0.5250 && m4 < 0.5336) || (m4 >= 0.5542 && m4 < 0.5628)) w = 1;
                if (w < 0) continue;
                dWin[w][smp] += 1;
                const int ip = PtBin(r.pt4), iy = YBin(r.y4);
                if (ip < 0 || iy < 0) { ++dNoBin; continue; }
                dK[((size_t)w * 3 + smp) * nGrid + (size_t)ip * kNY + iy] += 1;
            }
    }
    P("\nB1. data S0 rows: P us/pp/mm %.0f / %.0f / %.0f; side bands %.0f / %.0f / %.0f; rows without a grid bin %ld\n",
      dWin[0][0], dWin[0][1], dWin[0][2], dWin[1][0], dWin[1][1], dWin[1][2], dNoBin);
    P("KEY dwin P %.0f %.0f %.0f SB %.0f %.0f %.0f\n", dWin[0][0], dWin[0][1], dWin[0][2], dWin[1][0], dWin[1][1], dWin[1][2]);

    // =========================================================================================================
    // C. TAG FITS per class: data, emulated simulation, truth
    // =========================================================================================================
    P("\nC1. ECAL TAG per M(ee) class: frozen m5c fit (%.4f / %.1f MeV, pol%d, %.2f-%.2f); data N - CB of the tagged rows,\n", kTagMean,
      1e3 * kTagSigma, kTagOrder, kTagLo, kTagHi);
    P("    simulation = T1 + T2 tagged rows (emulated tag, UNVERIFIED) with the same fit; truth = C2 and C3 rows (V3), all m5c\n");
    for (int g = 0; g < kNTr; ++g) {
        P("  [%s]\n", kTrName[g]);
        for (int k = 0; k < kNTc; ++k) {
            const scan::Spectra& sd = dTag[g * kNTc + k];
            const scan::PeakFit fd = scan::FitPeakMain(sd, TagCfg());
            std::vector<double> ys(kTagNbin);
            for (int i = 0; i < kTagNbin; ++i) ys[i] = aTagS.sum[((size_t)g * kNTc + k) * kTagNbin + i] * kSimToRaw;
            const scan::Spectra ss = SimSpectra(ys.data(), kTagNbin, kTagHLo, kTagHHi, kTagReb);
            const scan::PeakFit fs = scan::FitPeakMain(ss, TagCfg());
            double tr[2][2], trv[2][2];
            for (int cc = 0; cc < 2; ++cc)
                for (int a = 0; a < 2; ++a) {
                    const size_t idx = (((size_t)g * kNTc + k) * 2 + cc) * 2 + a;
                    tr[cc][a] = aTagT.sum[idx] * kSimToRaw;
                    trv[cc][a] = aTagT.var[idx] * kSimToRaw * kSimToRaw;
                }
            P("    %-7s data T = %7.1f +- %5.1f (chi2 %5.1f/%d, N %6.0f) | sim T = %8.1f (frozen fit) | truth C2 %8.1f +- %6.1f, C3 %8.1f +- %6.1f (all m5c); "
              "C2 %8.1f, C3 %8.1f (0.49-0.62)\n",
              kTcName[k], fd.S, fd.eS, fd.chi2, fd.ndf, sd.nsel[0], fs.S, tr[0][0], std::sqrt(trv[0][0]), tr[1][0], std::sqrt(trv[1][0]), tr[0][1], tr[1][1]);
            if (gOut) fprintf(gOut, "KEY tag %d %d %s data %.4f %.4f %.3f %d ok %d sim %.4f ok %d C2 %.4f %.4f C3 %.4f %.4f C2w %.4f %.4f C3w %.4f %.4f\n", g, k,
                              kTcName[k], fd.S, fd.eS, fd.chi2, fd.ndf, (int)fd.ok, fs.S, (int)fs.ok, tr[0][0], std::sqrt(trv[0][0]), tr[1][0],
                              std::sqrt(trv[1][0]), tr[0][1], std::sqrt(trv[0][1]), tr[1][1], std::sqrt(trv[1][1]));
            // spectra for the figure
            TH1D* hd = new TH1D(Form("tag_data_r%d_c%d", g, k), "N - CB;m5c [GeV/c^{2}];rows", kTagNbin, kTagHLo, kTagHHi);
            TH1D* hs = new TH1D(Form("tag_sim_r%d_c%d", g, k), "T1+T2 tagged;m5c [GeV/c^{2}];expected raw rows", kTagNbin, kTagHLo, kTagHHi);
            hd->SetDirectory(nullptr); hs->SetDirectory(nullptr);
            for (int i = 0; i < kTagNbin; ++i) {
                hd->SetBinContent(i + 1, sd.sig[i]); hd->SetBinError(i + 1, std::sqrt(std::max(0., sd.sigVar[i])));
                hs->SetBinContent(i + 1, ys[i]);
                hs->SetBinError(i + 1, std::sqrt(aTagS.var[((size_t)g * kNTc + k) * kTagNbin + i]) * kSimToRaw);
            }
            fo->WriteTObject(hd, hd->GetName(), "Overwrite"); fo->WriteTObject(hs, hs->GetName(), "Overwrite");
            delete hd; delete hs;
        }
    }

    // =========================================================================================================
    // D. (pT4, y4) grids
    // =========================================================================================================
    for (int g = 0; g < kNKg; ++g) {
        std::vector<double> s(nGrid), v(nGrid);
        for (size_t i = 0; i < nGrid; ++i) { s[i] = aK.sum[(size_t)g * nGrid + i] / kLsim; v[i] = aK.var[(size_t)g * nGrid + i] / (kLsim * kLsim); }
        WriteH2(fo, Form("sim_%s", kKgName[g]), s, &v, 0);
    }
    const char* smpn[3] = {"us", "pp", "mm"};
    for (int w = 0; w < 2; ++w)
        for (int smp = 0; smp < 3; ++smp)
            WriteH2(fo, Form("data_%s_%s", w == 0 ? "P" : "SB", smpn[smp]), dK, nullptr, ((size_t)w * 3 + smp) * nGrid);
    // text dump of the grids (sparse) for the python step
    {
        FILE* fg = fopen(Form("%s/wp8_grid.txt", kOutDir), "w");
        fprintf(fg, "# wp8_sim.C grids: pT bin (20 MeV/c from 0), y bin (0.05 from -2.5); sim in nb (sum w / L_sim, var), data in raw rows\n");
        for (int g = 0; g < kNKg; ++g)
            for (size_t i = 0; i < nGrid; ++i)
                if (aK.sum[(size_t)g * nGrid + i] != 0.)
                    fprintf(fg, "sim %s %zu %zu %.8g %.8g\n", kKgName[g], i / kNY, i % kNY, aK.sum[(size_t)g * nGrid + i] / kLsim,
                            aK.var[(size_t)g * nGrid + i] / (kLsim * kLsim));
        for (int w = 0; w < 2; ++w)
            for (int smp = 0; smp < 3; ++smp)
                for (size_t i = 0; i < nGrid; ++i) {
                    const double c = dK[((size_t)w * 3 + smp) * nGrid + i];
                    if (c != 0.) fprintf(fg, "data %s %s %zu %zu %.0f\n", w == 0 ? "P" : "SB", smpn[smp], i / kNY, i % kNY, c);
                }
        fclose(fg);
    }
    fo->Close();
    P("\nwritten %s/wp8_sim.root, wp8_sim.txt, wp8_grid.txt; total %.0f s\n", kOutDir, total.RealTime());
    fclose(gOut); gOut = nullptr;
}
