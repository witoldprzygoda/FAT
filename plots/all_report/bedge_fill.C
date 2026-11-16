// ============================================================================
// bedge_fill.C  --  package bedge (the eta bump measured through its kinematic edge): the row pass
// ============================================================================
// Written 2026-11-14 (prefix bedge). Reads the data compact rows (scan_v2) and the simulated e+e- rows (std + lep +
// eta, SimEngine copy) ONCE and writes flat double arrays (layout in bedge_common.h and <out>/layout.txt):
//
//   data   N, N++, N-- per (configuration x M(ee) class x 5 MeV M(pipi) bin x 2.5 MeV M4 bin), configurations S0,
//          pre, S0 & isBest == 1, the four periods, the two folds; cor_ok == 1 rows only (as wp6_fill.C).
//          200 bootstrap replicas (Poisson(1) weight per event group, own seed; unlike-sign and like-sign groups
//          independent, as the WP0 CB replicas) of S0, pre and S0best in the sub-window M4 0.40-0.62,
//          M(pipi) 0.28-0.47.
//   sim    expected RAW data rows (sum w x RowFactor x 27.0718) per (configuration x group x class x M(pipi) bin x
//          M4 bin), with event-group variances; groups C2, C3, C4, C5, G_eta, G_D (omega/rho0 pions / other),
//          G_C (omega/rho0 / other); 15 configurations (nominal, s, r, M(pipi) scale, physics weights, pre,
//          isBest == 1, ownership mode B). 100 replicas (WP0 generator and seed) of T12, C2, C3 in the sub-window,
//          for every configuration (same events -> correlated variant differences). Window sums (bump window B,
//          all M4) with replicas; M(ee) of C2 + C3 near the edge (nominal) for the q-averaged step template.
//
// Selections and mass transfer exactly as wp6_fill.C: data S0 = PassSel(kSelS0, r, cor_ok == 1, dz_pi, mm4) and
// cor_ok == 1; simulation owned rows (mode A), PassSel(kSelS0, r, cor_ok == 1, 0, mm4), CompsOf(r, false),
// M4' = M4sim(m4, m4_true, s, r), M(pipi)' = 0.9952 M(pipi).
//
// Run (from /home/przygoda/HADES/FAT/FAT; about 12 GB, a few minutes):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_bedge", kTRUE)' \
//        'plots/all_report/bedge_fill.C+O'
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "TStopwatch.h"
#include "TSystem.h"

#include "hists_io.h"
#include "sim_eta_physics.h"
#include "bedge_common.h"

namespace bedgef {
using namespace bedge;
FILE* gLog = nullptr;
#define PL(...) do { printf(__VA_ARGS__); if (gLog) fprintf(gLog, __VA_ARGS__); fflush(stdout); } while (0)

// accumulator with event-group variance (copy of wp6_fill.C Acc)
struct Acc {
    std::vector<double> sum, var, tmp;
    std::vector<size_t> touched;
    void Init(size_t n) { sum.assign(n, 0.); var.assign(n, 0.); tmp.assign(n, 0.); touched.clear(); }
    inline void Add(size_t i, double w) {
        if (tmp[i] == 0.) touched.push_back(i);
        tmp[i] += w; sum[i] += w;
        if (tmp[i] == 0.) tmp[i] = 1e-300;
    }
    void Flush() {
        for (size_t i : touched) { var[i] += tmp[i] * tmp[i]; tmp[i] = 0.; }
        touched.clear();
    }
};
inline int GroupOf(const int* cn, int ncn, bool& isT12) {
    isT12 = false; int g = -1;
    for (int q = 0; q < ncn; ++q) {
        const int c = cn[q];
        if (c == kT12) isT12 = true;
        if (g >= 0) continue;
        switch (c) {
            case kC2: g = gC2; break;
            case kC3: g = gC3; break;
            case kC4: g = gC4; break;
            case kC5: g = gC5; break;
            case kGeta: g = gGeta; break;
            case kGpi0Dom: case kGvecom: case kGothom: g = gGDom; break;
            case kGpi0Dnom: case kGvecnom: case kGothnom: g = gGDnom; break;
            case kGconvom: g = gGCom; break;
            case kGconvnom: g = gGCnom; break;
        }
    }
    return g;
}
}  // namespace bedgef

void bedge_fill() {
    using namespace bedgef;
    TStopwatch total;
    gSystem->mkdir(kOut, kTRUE);
    gLog = fopen(Form("%s/bedge_fill.txt", kOut), "w");
    PL("# bedge_fill.C, 2026-11-14: row pass of the package bedge (edge step as a measurement)\n");
    PL("# data %s\n# sim  %s + %s\n", scan::kCompact, simeng::kCompact, simeng::kCompactEta);

    // ---------------------------------------------------------------------------------------------------- data
    std::vector<double> dN(SizeDN(), 0.), dR(SizeDR(), 0.);
    long nData[kNDCfg][3] = {{0}};
    long nNaN = 0;
    {
        scan::Data d;
        if (!scan::Load(d)) return;
        float eF[kNPP + 1];
        for (int k = 0; k <= kNPP; ++k) eF[k] = (float)PPEdge(k);
        const double bwD = (scan::kMhi - scan::kMlo) / scan::kNbin;
        std::vector<double> kr(kNRepD, 0.);
        TStopwatch sw;
        for (int k = 0; k < 3; ++k) {
            int cur = -1; bool haveK = false;
            for (const auto& r : d.s[k]) {
                if (r.grp != cur) { cur = r.grp; haveK = false; }
                if (r.cor_ok != 1) continue;
                const float mp = r.mpipi;
                if (!(mp >= eF[0] && mp < eF[kNPP])) continue;
                const int p = (int)(std::upper_bound(eF, eF + kNPP + 1, mp) - eF) - 1;
                if (p < 0 || p >= kNPP) continue;
                const double m = r.m4;
                if (!(m >= scan::kMlo && m < scan::kMhi)) continue;
                const int ib = (int)((m - scan::kMlo) / bwD);
                if (ib < 0 || ib >= kNM4) continue;
                const int e = (r.mee < 0.010f) ? 0 : (r.mee >= 0.010f ? 1 : 2);
                const bool s0 = PassSel(kSelS0, r, true, r.dz_pi, r.mm4);
                const bool pre = PassSel(kSelPre, r, true, r.dz_pi, r.mm4);
                const int per = PeriodOfDay(r.file_id);
                bool pass[kNDCfg];
                pass[dS0] = s0; pass[dPre] = pre; pass[dBest] = s0 && r.isBest == 1;
                pass[dPerI] = s0 && per == kPerI; pass[dPerIIA] = s0 && per == kPerIIA;
                pass[dPerIIB] = s0 && per == kPerIIB; pass[dPerIII] = s0 && per == kPerIII;
                pass[dFold0] = s0 && r.fold == 0; pass[dFold1] = s0 && r.fold == 1;
                for (int c = 0; c < kNDCfg; ++c)
                    if (pass[c]) { dN[IdxDN(c, e, p, k, ib)] += 1; ++nData[c][k]; }
                if (e == 2 && s0) ++nNaN;
                if (p >= kRPPlo && ib >= kRM4lo && ib < kRM4lo + kNRM4 && (pass[dS0] || pass[dPre] || pass[dBest])) {
                    if (!haveK) { GroupPois(kSeedData, k, r.grp, kNRepD, kr.data()); haveK = true; }
                    for (int c = 0; c < kNDRep; ++c) {
                        if (!pass[c]) continue;
                        double* q = &dR[IdxDR(c, e, p - kRPPlo, k, ib - kRM4lo)];
                        for (int rr = 0; rr < kNRepD; ++rr) q[rr] += kr[rr];
                    }
                }
            }
        }
        PL("A. data rows in 0.27 <= M(pipi) < 0.47, 0.35 <= M4 < 0.95 (us / pp / mm):\n");
        for (int c = 0; c < kNDCfg; ++c) PL("   %-7s %9ld %9ld %9ld\n", kDTag[c], nData[c][0], nData[c][1], nData[c][2]);
        PL("   S0 rows without an M(ee) value (class 2): %ld; %.1f s\n", nNaN, sw.RealTime());
    }
    WriteArr(Form("%s/data_N.bin", kOut), dN);
    WriteArr(Form("%s/data_R.bin", kOut), dR);
    { std::vector<double>().swap(dR); }

    // ---------------------------------------------------------------------------------------------------- simulation
    std::vector<Acc> acc(kNSCfg), accW(kNSCfg);
    for (int c = 0; c < kNSCfg; ++c) { acc[c].Init((size_t)kNG * kNCls * kNPP * kNM4); accW[c].Init((size_t)kNG * kNCls * kNWS); }
    std::vector<double> sR(SizeSR(), 0.), sWR(SizeSWR(), 0.), hq(SizeHq(), 0.);
    long nUsed[kNSCfg] = {0};
    {
        simphys::Init();
        simphys::Install();
        simphys::LoadAngles();
        simeng::SimData* D = new simeng::SimData();
        if (!simeng::Load(*D, simeng::kCompact, simeng::kMaskUS, -1e9, 1e9, 1)) return;
        const int trees[3] = {simeng::kStdUS, simeng::kLepUS, simeng::kEtaUS};
        double eD[kNPP + 1];
        for (int k = 0; k <= kNPP; ++k) eD[k] = PPEdge(k);
        std::vector<double> kr(kNRepS, 0.);
        const int nPhys = 4; const int physCodes[nPhys] = {kPhysNom, kPhysOff, kPhysBump1aY, kPhysBumpPS};
        TStopwatch sw;
        for (int t = 0; t < 3; ++t) {
            const std::vector<simeng::Row>& v = D->t[trees[t]];
            int cur = -1; bool haveK = false;
            for (size_t i = 0; i < v.size(); ++i) {
                const simeng::Row& r = v[i];
                if (r.grp != cur) {
                    for (int c = 0; c < kNSCfg; ++c) { acc[c].Flush(); accW[c].Flush(); }
                    cur = r.grp; haveK = false;
                }
                if (r.ownA != 1 && r.ownB != 1) continue;
                int cn[6];
                const int ncn = CompsOf(r, false, cn);
                if (ncn == 0) continue;
                bool isT12 = false;
                const int g = GroupOf(cn, ncn, isT12);
                if (g < 0 || !isT12) continue;
                const bool s0 = PassSel(kSelS0, r, r.cor_ok == 1, 0.f, r.mm4);
                const bool pre = PassSel(kSelPre, r, r.cor_ok == 1, 0.f, r.mm4);
                if (!s0 && !pre) continue;
                double pf[nPhys]; bool havePf[nPhys] = {false, false, false, false};
                auto phys = [&](int code) {
                    for (int q = 0; q < nPhys; ++q)
                        if (physCodes[q] == code) { if (!havePf[q]) { pf[q] = simphys::RowFactor(r, code); havePf[q] = true; } return pf[q]; }
                    return simphys::RowFactor(r, code);
                };
                const int e = (r.mee < 0.010) ? 0 : (r.mee >= 0.010 ? 1 : 2);
                for (int c = 0; c < kNSCfg; ++c) {
                    const SCfgDef& C = kSC[c];
                    if ((C.mode == 0 ? r.ownA : r.ownB) != 1) continue;
                    if (!(C.sel == kSelS0 ? s0 : pre)) continue;
                    if (C.best && r.isBest != 1) continue;
                    const double w = r.w * phys(C.phys) * kSimToRaw;
                    if (w == 0.) continue;
                    ++nUsed[c];
                    const double m4p = M4sim(r.m4, r.m4_true, C.s, C.r);
                    const double mpp = C.spp * r.mpipi;
                    // window sums (all M(pipi))
                    const bool inB = (m4p >= kBlo && m4p < kBhi);
                    accW[c].Add(((size_t)g * kNCls + e) * kNWS + 1, w);
                    if (inB) accW[c].Add(((size_t)g * kNCls + e) * kNWS + 0, w);
                    if (g == gC2 || g == gC3) {
                        if (!haveK) { GroupPois(kSeedSim, trees[t], r.grp, kNRepS, kr.data()); haveK = true; }
                        for (int ww = 0; ww < 2; ++ww) {
                            if (ww == 0 && !inB) continue;
                            double* q = &sWR[IdxSWR(c, g == gC2 ? 0 : 1, e, ww)];
                            for (int rr = 0; rr < kNRepS; ++rr) q[rr] += w * kr[rr];
                        }
                    }
                    // spectra
                    const int ib = (m4p >= scan::kMlo && m4p < scan::kMhi) ? (int)((m4p - scan::kMlo) / scan::kBinW) : -1;
                    if (ib < 0 || ib >= kNM4) continue;
                    if (!(mpp >= eD[0] && mpp < eD[kNPP])) continue;
                    const int p = (int)(std::upper_bound(eD, eD + kNPP + 1, mpp) - eD) - 1;
                    if (p < 0 || p >= kNPP) continue;
                    acc[c].Add((((size_t)g * kNCls + e) * kNPP + p) * kNM4 + ib, w);
                    if (p >= kRPPlo && ib >= kRM4lo && ib < kRM4lo + kNRM4) {
                        if (!haveK) { GroupPois(kSeedSim, trees[t], r.grp, kNRepS, kr.data()); haveK = true; }
                        const int rg[3] = {rT12, g == gC2 ? rC2 : -1, g == gC3 ? rC3 : -1};
                        for (int q3 = 0; q3 < 3; ++q3) {
                            if (rg[q3] < 0) continue;
                            double* q = &sR[IdxSR(c, rg[q3], e, p - kRPPlo, ib - kRM4lo)];
                            for (int rr = 0; rr < kNRepS; ++rr) q[rr] += w * kr[rr];
                        }
                    }
                    if (c == sNom && (g == gC2 || g == gC3)) {
                        const double ed = Edge(mpp, 0.);
                        if (ed > 0 && m4p >= ed - 0.030 && m4p < ed + 0.005) {
                            const int b = (int)(r.mee / 0.001);
                            if (b >= 0 && b < kNHq) hq[IdxHq(e, p, b)] += w;
                        }
                    }
                }
            }
            for (int c = 0; c < kNSCfg; ++c) { acc[c].Flush(); accW[c].Flush(); }
            PL("B. sim tree %s: %zu rows, %.0f s\n", simeng::TreeName(trees[t]), v.size(), sw.RealTime());
            sw.Continue();
        }
        delete D;
    }
    PL("   used rows per configuration:");
    for (int c = 0; c < kNSCfg; ++c) PL(" %s %ld", kSC[c].tag, nUsed[c]);
    PL("\n");
    // flatten
    std::vector<double> sS(SizeSS(), 0.), sV(SizeSS(), 0.), sW(SizeSW(), 0.), sWV(SizeSW(), 0.);
    for (int c = 0; c < kNSCfg; ++c) {
        const size_t o = IdxSS(c, 0, 0, 0, 0), n = (size_t)kNG * kNCls * kNPP * kNM4;
        std::copy(acc[c].sum.begin(), acc[c].sum.end(), sS.begin() + o);
        std::copy(acc[c].var.begin(), acc[c].var.end(), sV.begin() + o);
        (void)n;
        const size_t ow = IdxSW(c, 0, 0, 0);
        std::copy(accW[c].sum.begin(), accW[c].sum.end(), sW.begin() + ow);
        std::copy(accW[c].var.begin(), accW[c].var.end(), sWV.begin() + ow);
    }
    WriteArr(Form("%s/sim_S.bin", kOut), sS);
    WriteArr(Form("%s/sim_V.bin", kOut), sV);
    WriteArr(Form("%s/sim_R.bin", kOut), sR);
    WriteArr(Form("%s/sim_W.bin", kOut), sW);
    WriteArr(Form("%s/sim_WV.bin", kOut), sWV);
    WriteArr(Form("%s/sim_WR.bin", kOut), sWR);
    WriteArr(Form("%s/sim_hq.bin", kOut), hq);
    // quick checks: nominal C2 + C3 in B, all M(ee), x 64/63 (WP6/v8: 54758.9 = SMASH (C2 + C3) in B in yield units)
    double b23 = 0, b23e2 = 0;
    for (int g = gC2; g <= gC3; ++g)
        for (int e = 0; e < kNCls; ++e) { b23 += sW[IdxSW(sNom, g, e, 0)]; if (e == 1) b23e2 += sW[IdxSW(sNom, g, e, 0)]; }
    PL("C. nominal (C2 + C3) in B: all M(ee) %.1f raw rows (x 64/63 = %.1f); E2 %.1f (x 64/63 = %.1f)\n", b23, b23 * kPT3, b23e2, b23e2 * kPT3);
    FILE* fl = fopen(Form("%s/layout.txt", kOut), "w");
    if (fl) {
        fprintf(fl, "# bedge arrays (flat double, little endian, first 8 bytes = uint64 count); index functions in bedge_common.h\n");
        fprintf(fl, "data_N.bin  [dcfg %d][cls %d][pp %d][k 3: us pp mm][m4 %d]\n", kNDCfg, kNCls, kNPP, kNM4);
        fprintf(fl, "data_R.bin  [rcfg %d: S0 pre S0best][cls][pp %d from bin %d][k 3][m4 %d from bin %d][rep %d]\n", kNDRep, kNRPP, kRPPlo, kNRM4, kRM4lo, kNRepD);
        fprintf(fl, "sim_S.bin, sim_V.bin  [scfg %d][group %d][cls][pp %d][m4 %d]  (raw data rows; event-group variance)\n", kNSCfg, kNG, kNPP, kNM4);
        fprintf(fl, "sim_R.bin   [scfg][rgroup 3: T12 C2 C3][cls][pp %d from %d][m4 %d from %d][rep %d]\n", kNRPP, kRPPlo, kNRM4, kRM4lo, kNRepS);
        fprintf(fl, "sim_W.bin, sim_WV.bin  [scfg][group][cls][w 2: B, all M4] (all M(pipi))\n");
        fprintf(fl, "sim_WR.bin  [scfg][2: C2 C3][cls][w 2][rep %d]\n", kNRepS);
        fprintf(fl, "sim_hq.bin  [cls][pp][150 x 1 MeV M(ee)]  C2 + C3, nominal, M4' in [E(m,0) - 30 MeV, E(m,0) + 5 MeV]\n");
        fprintf(fl, "M4 bins: 0.35 + 0.0025 i; M(pipi) bins: 0.27 + 0.005 p (data scale); cls 0 E1 < 10 MeV, 1 E2, 2 no M(ee)\n");
        fprintf(fl, "dcfg: "); for (int c = 0; c < kNDCfg; ++c) fprintf(fl, "%d %s; ", c, kDTag[c]); fprintf(fl, "\n");
        fprintf(fl, "scfg: "); for (int c = 0; c < kNSCfg; ++c) fprintf(fl, "%d %s (%s); ", c, kSC[c].tag, kSC[c].text); fprintf(fl, "\n");
        fprintf(fl, "groups: "); for (int g = 0; g < kNG; ++g) fprintf(fl, "%d %s; ", g, kGTag[g]); fprintf(fl, "\n");
        fclose(fl);
    }
    PL("done in %.0f s\n", total.RealTime());
    if (gLog) fclose(gLog);
}
