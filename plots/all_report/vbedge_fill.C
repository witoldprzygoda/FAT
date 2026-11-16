// ============================================================================
// vbedge_fill.C  --  INDEPENDENT VERIFIER of the package bedge: own row pass
// ============================================================================
// Written 2026-11-15 (prefix vbedge). Does NOT use bedge_common.h, bedge_fill.C or the bedge arrays.
// Own tree loops (data: own SetBranchAddress on the compact data trees; simulation: own loop over the
// compact simulated trees with simeng::BindRow), own binning, own selection code (written from the
// definitions in defs.h Sect. D), own bootstrap generator (std::mt19937_64 + std::poisson_distribution, own seeds).
// Frozen library pieces used on purpose (they are definitions, not the package's logic):
//   simphys::RowFactor(row, 108) (physics weights), alldefs::CompsOf (class of a row), kSimToRaw = 27.0718.
//
// Output (/hdd3/przygoda/HADES/FAT_out/all_report/vbedge/), flat float64 little endian, no header:
//   dN.bin   [sample 3: us pp mm][cls 3: E1 E2 NaN][pp 40: 0.27 + 5 MeV][m4 240: 0.35 + 2.5 MeV]      S0 data
//   dNpre.bin  same, selection pre
//   dB.bin   [rep 200][sample 3][cls 3][pp 40][m4 240]   S0 data bootstrap (Poisson(1) per event group)
//   sS.bin   [cfg 4: nominal, r = 0.89, bump phase space (100), weights off (0)][grp 6: C2 C3 T12 OtherEta G_D G_C][cls 3][pp 40][m4 240]  expected raw rows
//   sR.bin   [rep 100][cfg 4][grp 3: C2 C3 T12][cls 3][pp 40][m4 240]   simulation replicas (own generator)
//   sW.bin   [grp 2: C2 C3][cls 3][win 2: B, any M4]   nominal window sums (all M(pipi))
//   sWR.bin  [rep 100][grp 2][cls 3][win 2]
//
// Run (from /home/przygoda/HADES/FAT/FAT):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/vbedge/aclic", kTRUE)' \
//        'plots/all_report/vbedge_fill.C+O'
// ============================================================================
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

#include "TFile.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TTree.h"

#include "sim_eta_physics.h"
#include "defs.h"

namespace vbf {
const char* const kOutDir = "/hdd3/przygoda/HADES/FAT_out/all_report/vbedge";
const int NM4 = 240, NPP = 40, NCL = 3, NS = 3;
const double M4LO = 0.35, M4W = 0.0025, PPLO = 0.27, PPW = 0.005;
const int NRD = 200, NRS = 100;
inline size_t iD(int s, int c, int p, int m) { return (((size_t)s * NCL + c) * NPP + p) * NM4 + m; }
const size_t SZD = (size_t)NS * NCL * NPP * NM4;
const int NCFG = 4, NGS = 6, NGR = 3;   // cfg: 0 nominal, 1 r = 0.89, 2 bump phase space (code 100), 3 physics weights off (0)
inline size_t iS(int cf, int g, int c, int p, int m) { return ((((size_t)cf * NGS + g) * NCL + c) * NPP + p) * NM4 + m; }
const size_t SZS = (size_t)NCFG * NGS * NCL * NPP * NM4;
inline size_t iR(int cf, int g, int c, int p, int m) { return ((((size_t)cf * NGR + g) * NCL + c) * NPP + p) * NM4 + m; }
const size_t SZR1 = (size_t)NCFG * NGR * NCL * NPP * NM4;
inline int M4Bin(double m) { if (!(m >= M4LO)) return -1; const int i = (int)std::floor((m - M4LO) / M4W); return i < NM4 ? i : -1; }
inline int PPBin(double m) { if (!(m >= PPLO)) return -1; const int i = (int)std::floor((m - PPLO) / PPW + 1e-12); return i < NPP ? i : -1; }
inline int Cls(float mee) { if (!(mee == mee)) return 2; return mee < 0.010f ? 0 : 1; }
// own 64-bit hash (Murmur3 finaliser) for the per-group seeds
inline uint64_t H64(uint64_t k) {
    k ^= k >> 33; k *= 0xff51afd7ed558ccdULL; k ^= k >> 33; k *= 0xc4ceb9fe1a85ec53ULL; k ^= k >> 33; return k;
}
void Poisson1(uint64_t key, int n, std::vector<int>& out) {
    std::mt19937_64 g(H64(key));
    std::poisson_distribution<int> P(1.0);
    out.resize(n);
    for (int i = 0; i < n; ++i) out[i] = P(g);
}
bool Dump(const std::string& fn, const std::vector<double>& v) {
    FILE* f = fopen(fn.c_str(), "wb");
    if (!f) { printf("cannot write %s\n", fn.c_str()); return false; }
    fwrite(v.data(), sizeof(double), v.size(), f); fclose(f);
    printf("wrote %s (%zu values)\n", fn.c_str(), v.size());
    return true;
}
}  // namespace vbf

void vbedge_fill(int dataOnly = 0) {
    using namespace vbf;
    TStopwatch tot;
    gSystem->mkdir(kOutDir, kTRUE);
    // ------------------------------------------------------------------------------------------- data
    {
        std::vector<double> N(SZD, 0.), Np(SZD, 0.), B((size_t)NRD * SZD, 0.), Sub(6 * SZD, 0.);   // Sub: fold0 fold1 perI perIIA perIIB perIII
        TFile* f = TFile::Open(scan::kCompact, "READ");
        if (!f || f->IsZombie()) { printf("cannot open data\n"); return; }
        const char* tn[3] = {"us", "pp", "mm"};
        long nsel[3] = {0, 0, 0}, npre[3] = {0, 0, 0};
        std::vector<int> kw;
        for (int s = 0; s < 3; ++s) {
            TTree* t = (TTree*)f->Get(tn[s]);
            Int_t isBest, cor_ok, pipv, pimv, grp, fid;
            Float_t vz, m4, mpipi, mee;
            t->SetBranchStatus("*", 0);
            for (const char* b : {"isBest", "cor_ok", "pip_isUsedVertex", "pim_isUsedVertex", "grp", "vz", "m4", "mpipi", "mee", "file_id"}) t->SetBranchStatus(b, 1);
            t->SetBranchAddress("file_id", &fid);
            t->SetBranchAddress("isBest", &isBest); t->SetBranchAddress("cor_ok", &cor_ok);
            t->SetBranchAddress("pip_isUsedVertex", &pipv); t->SetBranchAddress("pim_isUsedVertex", &pimv);
            t->SetBranchAddress("grp", &grp); t->SetBranchAddress("vz", &vz); t->SetBranchAddress("m4", &m4);
            t->SetBranchAddress("mpipi", &mpipi); t->SetBranchAddress("mee", &mee);
            const Long64_t n = t->GetEntries();
            int lastG = -1;
            for (Long64_t i = 0; i < n; ++i) {
                t->GetEntry(i);
                // pre = isBest >= 0, vz > -500, PT3 (all data rows), vz > -160, cor_ok == 1; S0 = pre and both pions in the vertex fit
                const bool pre = isBest >= 0 && vz > -500.f && vz > -160.f && cor_ok == 1;
                if (!pre) continue;
                const bool s0 = pipv == 1 && pimv == 1;
                const int p = PPBin(mpipi), m = M4Bin(m4);
                if (p < 0 || m < 0) continue;
                const int c = Cls(mee);
                Np[iD(s, c, p, m)] += 1; ++npre[s];
                if (!s0) continue;
                N[iD(s, c, p, m)] += 1; ++nsel[s];
                // subsets: event-group parity; period by day (I 034-046, IIA 047-048, IIB 049-053, III 059-067)
                Sub[(size_t)(grp & 1) * SZD + iD(s, c, p, m)] += 1;
                const int per = (fid >= 34 && fid <= 46) ? 0 : (fid >= 47 && fid <= 48) ? 1 : (fid >= 49 && fid <= 53) ? 2 : (fid >= 59 && fid <= 67) ? 3 : -1;
                if (per >= 0) Sub[(size_t)(2 + per) * SZD + iD(s, c, p, m)] += 1;
                if (grp != lastG) { Poisson1(((uint64_t)(s + 1) << 40) ^ (uint64_t)(uint32_t)grp ^ 0x7E57ULL, NRD, kw); lastG = grp; }
                for (int r = 0; r < NRD; ++r) if (kw[r]) B[(size_t)r * SZD + iD(s, c, p, m)] += kw[r];
            }
        }
        f->Close();
        printf("data S0 rows (0.27 <= M(pipi) < 0.47, 0.35 <= M4 < 0.95): us %ld pp %ld mm %ld; pre: %ld %ld %ld\n", nsel[0], nsel[1], nsel[2], npre[0], npre[1], npre[2]);
        Dump(Form("%s/dN.bin", kOutDir), N); Dump(Form("%s/dNpre.bin", kOutDir), Np); Dump(Form("%s/dB.bin", kOutDir), B); Dump(Form("%s/dSub.bin", kOutDir), Sub);
    }
    if (dataOnly) { printf("data only: done in %.0f s\n", tot.RealTime()); return; }
    // ------------------------------------------------------------------------------------------- simulation
    {
        simphys::Init(); simphys::Install(); simphys::LoadAngles();
        std::vector<double> S(SZS, 0.), R((size_t)NRS * SZR1, 0.), W(2 * NCL * 2, 0.), WR((size_t)NRS * 2 * NCL * 2, 0.);
        const double kS = 0.99615, kSPP = 0.9952, toRaw = alldefs::kL / alldefs::kLsim / (64.0 / 63.0);
        printf("sum w -> raw rows: %.5f\n", toRaw);
        const double rr[NCFG] = {1.0, 0.89, 1.0, 1.0};
        const int phc[NCFG] = {108, 108, 100, 0};
        struct Src { const char* file; const char* tree; int id; };
        const Src src[3] = {{simeng::kCompact, "std_us", 0}, {simeng::kCompact, "lep_us", 3}, {simeng::kCompactEta, "eta_us", 6}};
        std::vector<int> kw;
        long nUsed = 0;
        double sumC23B = 0;
        for (int q = 0; q < 3; ++q) {
            TFile* f = TFile::Open(src[q].file, "READ");
            TTree* t = (TTree*)f->Get(src[q].tree);
            simeng::Row r;
            simeng::BindRow(t, r, false);
            const Long64_t n = t->GetEntries();
            int lastG = -1;
            for (Long64_t i = 0; i < n; ++i) {
                t->GetEntry(i);
                if (r.ownA != 1) continue;                       // ownership mode A
                // S0 for simulated rows: isBest >= 0, vz > -500, PT3 = trigbit & 8192, vz > -160, cor_ok == 1, both pions in the vertex
                if (!(r.isBest >= 0 && r.vz > -500.f && (r.trigbit & 8192) != 0 && r.vz > -160.f && r.cor_ok == 1 &&
                      r.pip_isUsedVertex == 1 && r.pim_isUsedVertex == 1)) continue;
                int cn[6];
                const int ncn = alldefs::CompsOf(r, false, cn);
                bool t12 = false; int cl = -1;   // 0 C2, 1 C3, 2 other eta (C4, C5, G_eta), 3 G_D, 4 G_C
                for (int k = 0; k < ncn; ++k) {
                    if (cn[k] == alldefs::kT12) t12 = true;
                    if (cn[k] == alldefs::kC2) cl = 0;
                    else if (cn[k] == alldefs::kC3) cl = 1;
                    else if (cn[k] == alldefs::kC4 || cn[k] == alldefs::kC5 || cn[k] == alldefs::kGeta) cl = 2;
                    else if (cl < 0 && cn[k] == alldefs::kGD) cl = 3;
                    else if (cl < 0 && cn[k] == alldefs::kGconv) cl = 4;
                }
                if (!t12 || cl < 0) continue;
                const double w0 = r.w * toRaw;
                double wc[NCFG];
                for (int cf = 0; cf < NCFG; ++cf) wc[cf] = w0 * simphys::RowFactor(r, phc[cf]);
                if (wc[0] == 0. && wc[2] == 0. && wc[3] == 0.) continue;
                ++nUsed;
                if (r.grp != lastG) { Poisson1(((uint64_t)(src[q].id + 11) << 40) ^ (uint64_t)(uint32_t)r.grp ^ 0x51DE5ULL, NRS, kw); lastG = r.grp; }
                const int c = Cls(r.mee);
                const double mpp = kSPP * r.mpipi;
                const int p = PPBin(mpp);
                for (int cf = 0; cf < NCFG; ++cf) {
                    const double w = wc[cf];
                    if (w == 0.) continue;
                    const double m4p = (r.m4_true > 0.f) ? kS * (r.m4_true + rr[cf] * (r.m4 - r.m4_true)) : kS * r.m4;
                    if (cf == 0 && cl <= 1) {
                        const bool inB = m4p >= 0.4550 && m4p < 0.5225;
                        for (int ww = 0; ww < 2; ++ww) {
                            if (ww == 0 && !inB) continue;
                            const size_t j = ((size_t)cl * NCL + c) * 2 + ww;
                            W[j] += w;
                            for (int k = 0; k < NRS; ++k) if (kw[k]) WR[(size_t)k * (2 * NCL * 2) + j] += w * kw[k];
                        }
                        if (inB) sumC23B += w;
                    }
                    const int m = M4Bin(m4p);
                    if (m < 0 || p < 0) continue;
                    // groups of S: 0 C2, 1 C3, 2 T12 (everything), 3 other eta, 4 G_D, 5 G_C
                    const int gS = cl <= 1 ? cl : cl + 1;
                    S[iS(cf, gS, c, p, m)] += w;
                    S[iS(cf, 2, c, p, m)] += w;
                    // replicas: C2, C3, T12
                    for (int k = 0; k < NRS; ++k) {
                        if (!kw[k]) continue;
                        const double x = w * kw[k];
                        double* base = &R[(size_t)k * SZR1];
                        if (cl <= 1) base[iR(cf, cl, c, p, m)] += x;
                        base[iR(cf, 2, c, p, m)] += x;
                    }
                }
            }
            t->ResetBranchAddresses();
            f->Close();
            printf("sim tree %s done, %.0f s\n", src[q].tree, tot.RealTime()); tot.Continue();
        }
        // note: group index 2 of S holds T12 = C2 + C3 + other eta + G_D + G_C; index 3 other eta, 4 G_D, 5 G_C
        printf("sim rows used: %ld; C2 + C3 in B (nominal, all M(pipi)): %.2f raw rows, x 64/63 = %.2f\n", nUsed, sumC23B, sumC23B * 64. / 63.);
        Dump(Form("%s/sS.bin", kOutDir), S); Dump(Form("%s/sR.bin", kOutDir), R); Dump(Form("%s/sW.bin", kOutDir), W); Dump(Form("%s/sWR.bin", kOutDir), WR);
    }
    printf("done in %.0f s\n", tot.RealTime());
}
