// ============================================================================
// e1v_verify.C  --  independent verifier of package e1 (prefix e1v)
// ============================================================================
// Own code (does not use e1_fill.C / e1_ana.py, nor the wp0 accumulators): streams
// the data compact rows (scan_v2) and the SMASH e+e- rows (std_us, lep_us, eta_us)
// and writes raw sums per selection to /hdd3/przygoda/HADES/FAT_out/all_report/e1v/:
//   data: fine M4 counts (2.5 MeV, 0.30-1.10) of us, pp, mm with the event-group
//         variance of us; exact-M4 window counts; like-sign counts per 50 MeV band
//         inside / outside the X17 M(ee) window (for the CB of the X17 window)
//   sim : sum w x RowFactor(108), own ownA rows, per component (eta = C2-C5 + G_eta,
//         non-eta = rest of template 2, C5, C6 (weight 2), C2+C3, G_eta, G_conv),
//         fine M4' bins with event-group variance, exact-M4' windows, C6 opening
//         angle and pair momentum in S0
// Selections (own implementation of pre / S0 and of the cuts of sim_eta_enriched.C):
//   pre = isBest >= 0, vz > -500, vz > -160, cor_ok == 1, PT3 (sim: trigbit & 8192)
//   S0  = pre + pip_isUsedVertex == 1 + pim_isUsedVertex == 1
//   E1 = oa_lab < 50, oa_eta > 140, M(pipi) < 0.46, OA > 2;  E2 = E1 without OA > 2
//   X1 = E2 + OA < 6;  OA6 = OA < 6 (X3 on S0)
// SMASH: M4' = 0.99615 M4 (r = 1), M(pipi)' = 0.9952 M(pipi) in the cuts (Part III);
//   extra variants with the raw M(pipi) and (old convention) the raw M4.
// Run (from /home/przygoda/HADES/FAT/FAT):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_e1v", kTRUE)' 'plots/all_report/e1v_verify.C+O(1)'  data
//   ... 'plots/all_report/e1v_verify.C+O(2)'  SMASH
// ============================================================================
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "TFile.h"
#include "TStopwatch.h"
#include "TTree.h"

#include "defs.h"
#include "sim_eta_physics.h"

namespace e1v {

const char* const kOut = "/hdd3/przygoda/HADES/FAT_out/all_report/e1v";
const int NF = 320;
const double FLO = 0.30, FW = 0.0025;
inline int Fine(double m) {
    if (!(m >= FLO && m < 1.10)) return -1;
    int i = (int)((m - FLO) / FW);
    return i >= NF ? NF - 1 : i;
}
// windows by exact M4 (data) / M4' (sim)
enum { wB = 0, wP, w56, wCR, wPB, wX17, wX2D, wX2Da, wFit, wR5, NW };
const char* const kWn[NW] = {"B", "P", "r56", "CR", "PB", "X17", "X2D", "X2Dall", "Fit", "R5"};
const double kWlo[NW] = {0.4550, 0.53364, 0.56, 0.575, 0.5325, 0.5325, 0.5225, 0.5225, 0.40, 0.565};
const double kWhi[NW] = {0.5225, 0.55416, 0.70, 0.70, 0.5575, 0.5575, 0.6100, 0.6100, 0.70, 0.70};
inline bool InW(int w, double m4, double mee) {
    if (!(m4 >= kWlo[w] && m4 < kWhi[w])) return false;
    if (w == wX17) return mee >= 0.0155 && mee < 0.0245;
    if (w == wX2D) return mee >= 0.002 && mee < 0.062;
    return true;
}
const int NOA = 41;   // 1 deg bins 0-40, last = overflow
// cells of a block: fine NF | windows NW | oa in PB NOA | oa in X17 window NOA
const int oFine = 0, oWin = NF, oOaPB = NF + NW, oOaX = NF + NW + NOA, NCELL = NF + NW + 2 * NOA;

// cut sets on top of the base
enum { cNone = 0, cE1, cE2, cX1, cOA6, NCUT };
inline bool Cut(int c, double mpp, double oa, double oal, double oae) {
    switch (c) {
        case cNone: return true;
        case cOA6: return oa < 6.;
        case cE1: return oal < 50. && oae > 140. && mpp < 0.46 && oa > 2.;
        case cE2: return oal < 50. && oae > 140. && mpp < 0.46;
        case cX1: return oal < 50. && oae > 140. && mpp < 0.46 && oa < 6.;
    }
    return false;
}
// selections: base (0 pre, 1 S0) x cut
struct SelD { const char* name; int base, cut; };
const int NS = 8;
const SelD kS[NS] = {{"pre", 0, cNone}, {"E1", 0, cE1}, {"E2", 0, cE2}, {"X1", 0, cX1},
                     {"S0", 1, cNone}, {"S0+E1", 1, cE1}, {"X3", 1, cOA6}, {"S0+X1", 1, cX1}};

struct GAcc {
    std::vector<double> s, v, cur;
    std::vector<int> touched;
    std::vector<char> mark;
    void Init(size_t n) { s.assign(n, 0.); v.assign(n, 0.); cur.assign(n, 0.); mark.assign(n, 0); }
    inline void Add(size_t i, double w) {
        s[i] += w;
        if (!mark[i]) { mark[i] = 1; touched.push_back((int)i); }
        cur[i] += w;
    }
    void Flush() {
        for (int i : touched) { v[i] += cur[i] * cur[i]; cur[i] = 0.; mark[i] = 0; }
        touched.clear();
    }
};

void Arr(FILE* f, const char* key, const double* a, size_t n, bool comma = true) {
    fprintf(f, "%s\"%s\": [", comma ? ", " : "", key);
    for (size_t i = 0; i < n; ++i) fprintf(f, i ? ",%.12g" : "%.12g", std::isfinite(a[i]) ? a[i] : 0.);
    fprintf(f, "]");
}

inline void FillCells(GAcc& A, size_t base, double m4, double mee, double oa, double w) {
    const int fi = Fine(m4);
    if (fi >= 0) A.Add(base + oFine + fi, w);
    for (int q = 0; q < NW; ++q) if (InW(q, m4, mee)) A.Add(base + oWin + q, w);
    if (m4 >= kWlo[wPB] && m4 < kWhi[wPB] && oa == oa && oa >= 0) {
        const int j = std::min(NOA - 1, (int)oa);
        A.Add(base + oOaPB + j, w);
        if (mee >= 0.0155 && mee < 0.0245) A.Add(base + oOaX + j, w);
    }
}

// =========================================================================== data
void RunData() {
    TStopwatch sw;
    TFile* f = TFile::Open(scan::kCompact, "READ");
    if (!f || f->IsZombie()) { printf("e1v: cannot open data\n"); return; }
    GAcc A[3];
    for (int k = 0; k < 3; ++k) A[k].Init((size_t)NS * NCELL);
    // like-sign rows per 50 MeV band inside / outside the X17 M(ee) window [sel][band][in/out] (pp + mm)
    std::vector<double> lsx((size_t)NS * 16 * 2, 0.);
    long long nread[3] = {0, 0, 0}, npass[NS][3] = {};
    for (int k = 0; k < 3; ++k) {
        TTree* t = (TTree*)f->Get(scan::SampleName(k));
        scan::Row r;
        scan::BindRow(t, r, false);
        t->SetBranchStatus("*", 0);
        for (const char* b : {"grp", "isBest", "vz", "cor_ok", "pip_isUsedVertex", "pim_isUsedVertex", "m4", "mee", "mpipi", "oa",
                              "oa_lab", "oa_eta"})
            t->SetBranchStatus(b, 1);
        const Long64_t n = t->GetEntries();
        int cur = -1;
        for (Long64_t i = 0; i < n; ++i) {
            t->GetEntry(i);
            nread[k]++;
            if (r.grp != cur) { A[k].Flush(); cur = r.grp; }
            if (!(r.isBest >= 0 && r.vz > -500.f && r.vz > -160.f && r.cor_ok == 1)) continue;
            if (!(r.m4 == r.m4)) continue;
            const bool s0 = r.pip_isUsedVertex == 1 && r.pim_isUsedVertex == 1;
            for (int s = 0; s < NS; ++s) {
                if (kS[s].base == 1 && !s0) continue;
                if (!Cut(kS[s].cut, r.mpipi, r.oa, r.oa_lab, r.oa_eta)) continue;
                npass[s][k]++;
                FillCells(A[k], (size_t)s * NCELL, r.m4, r.mee, r.oa, 1.);
                if (k > 0) {
                    const int fi = Fine(r.m4);
                    if (fi >= 0) {
                        const int band = std::min(15, fi / 20);
                        const int in = (r.mee >= 0.0155f && r.mee < 0.0245f) ? 0 : 1;
                        lsx[((size_t)s * 16 + band) * 2 + in] += 1.;
                    }
                }
            }
        }
        A[k].Flush();
        t->ResetBranchAddresses();
        printf("e1v: data %s %lld rows (%.1f s)\n", scan::SampleName(k), nread[k], sw.RealTime());
        sw.Continue();
    }
    f->Close();
    FILE* o = fopen(Form("%s/e1v_data.json", kOut), "w");
    fprintf(o, "{\"ncell\": %d, \"oWin\": %d, \"oOaPB\": %d, \"oOaX\": %d, \"wins\": [", NCELL, oWin, oOaPB, oOaX);
    for (int q = 0; q < NW; ++q) fprintf(o, "%s\"%s\"", q ? ", " : "", kWn[q]);
    fprintf(o, "], \"nread\": [%lld, %lld, %lld], \"sel\": {", nread[0], nread[1], nread[2]);
    for (int s = 0; s < NS; ++s) {
        fprintf(o, "%s\"%s\": {\"npass\": [%lld, %lld, %lld]", s ? ", " : "", kS[s].name, npass[s][0], npass[s][1], npass[s][2]);
        const size_t b = (size_t)s * NCELL;
        Arr(o, "us", &A[0].s[b], NCELL); Arr(o, "usV", &A[0].v[b], NCELL);
        Arr(o, "pp", &A[1].s[b], NCELL); Arr(o, "mm", &A[2].s[b], NCELL);
        Arr(o, "lsx", &lsx[(size_t)s * 32], 32);
        fprintf(o, "}");
    }
    fprintf(o, "}}\n");
    fclose(o);
    printf("e1v: data done (%.1f s)\n", sw.RealTime());
}

// =========================================================================== simulation
enum { kEta = 0, kNon, kC5, kC6, kC23, kGeta, kGconv, NC };
const char* const kCn[NC] = {"eta", "non", "C5", "C6", "C23", "Geta", "Gconv"};
// sim selections: the 8 Part III ones (M4', scaled M(pipi)), then variants
const int NSS = NS + 4;   // + E1raw, X1raw (M4', raw M(pipi)); pre_old, E1_old (raw M4, raw M(pipi))
const char* const kSSn[NSS] = {"pre", "E1", "E2", "X1", "S0", "S0+E1", "X3", "S0+X1", "E1raw", "X1raw", "pre_old", "E1_old"};

void RunSim() {
    TStopwatch sw;
    simphys::Init();
    simphys::Install();
    simphys::LoadAngles();
    GAcc A;
    A.Init((size_t)NSS * NC * NCELL);
    // C6 pair momentum in S0, peak band (all M(ee)) and X17 window: weighted lists
    std::vector<double> pPB, pX, wPB, wX;
    double nOther1 = 0, nOtherW = 0;   // template-1 rows with cls outside 2..6 (should be none)
    long long nread = 0, nown = 0;
    const char* files[3] = {simeng::kCompact, simeng::kCompact, simeng::kCompactEta};
    const char* trees[3] = {"std_us", "lep_us", "eta_us"};
    for (int tt = 0; tt < 3; ++tt) {
        TFile* f = TFile::Open(files[tt], "READ");
        TTree* t = (TTree*)f->Get(trees[tt]);
        simeng::Row r;
        simeng::BindRow(t, r, false);
        const Long64_t n = t->GetEntries();
        int cur = -1;
        for (Long64_t i = 0; i < n; ++i) {
            t->GetEntry(i);
            nread++;
            if (r.grp != cur) { A.Flush(); cur = r.grp; }
            if (r.ownA != 1) continue;
            nown++;
            if (!(r.isBest >= 0 && r.vz > -500.f && r.vz > -160.f && r.cor_ok == 1 && (r.trigbit & 8192) != 0)) continue;
            // components
            int comps[4], nc = 0;
            if (r.tmpl == 1) {
                if (r.cls >= 2 && r.cls <= 5) {
                    comps[nc++] = kEta;
                    if (r.cls == 5) comps[nc++] = kC5;
                    if (r.cls == 2 || r.cls == 3) comps[nc++] = kC23;
                } else if (r.cls == 6) comps[nc++] = kC6;
                else { nOther1 += 1; nOtherW += r.w; }
            } else if (r.tmpl == 2) {
                if (r.pcls == 2 || r.lcls == 2) { comps[nc++] = kEta; comps[nc++] = kGeta; }
                else { comps[nc++] = kNon; if (r.lcls == 3) comps[nc++] = kGconv; }
            }
            if (nc == 0) continue;
            const double w = r.w * simphys::RowFactor(r, 108);
            const bool s0 = r.pip_isUsedVertex == 1 && r.pim_isUsedVertex == 1;
            const double m4p = 0.99615 * r.m4, mpps = 0.9952 * r.mpipi;
            for (int s = 0; s < NSS; ++s) {
                int base = 0, cut = 0;
                double m4 = m4p, mpp = mpps;
                if (s < NS) { base = kS[s].base; cut = kS[s].cut; }
                else if (s == NS) { cut = cE1; mpp = r.mpipi; }
                else if (s == NS + 1) { cut = cX1; mpp = r.mpipi; }
                else if (s == NS + 2) { cut = cNone; m4 = r.m4; }
                else { cut = cE1; m4 = r.m4; mpp = r.mpipi; }
                if (base == 1 && !s0) continue;
                if (!Cut(cut, mpp, r.oa, r.oa_lab, r.oa_eta)) continue;
                for (int q = 0; q < nc; ++q) {
                    const double wc = comps[q] == kC6 ? 2. : w;
                    FillCells(A, ((size_t)s * NC + comps[q]) * NCELL, m4, r.mee, r.oa, wc);
                }
                if (s == 4 && comps[0] == kC6 && m4 >= 0.5325 && m4 < 0.5575) {
                    pPB.push_back(r.p_ee); wPB.push_back(2.);
                    if (r.mee >= 0.0155f && r.mee < 0.0245f) { pX.push_back(r.p_ee); wX.push_back(2.); }
                }
            }
        }
        A.Flush();
        t->ResetBranchAddresses();
        f->Close();
        printf("e1v: sim %s done, %lld rows read, %lld owned A (%.1f s)\n", trees[tt], nread, nown, sw.RealTime());
        sw.Continue();
    }
    printf("e1v: template-1 rows with cls outside 2..6: %.0f (sum w %.4g)\n", nOther1, nOtherW);
    FILE* o = fopen(Form("%s/e1v_sim.json", kOut), "w");
    fprintf(o, "{\"ncell\": %d, \"oWin\": %d, \"oOaPB\": %d, \"oOaX\": %d, \"Lsim\": 221.05, \"nread\": %lld, \"nown\": %lld, \"nOther1\": %.0f, \"comps\": [",
            NCELL, oWin, oOaPB, oOaX, nread, nown, nOther1);
    for (int q = 0; q < NC; ++q) fprintf(o, "%s\"%s\"", q ? ", " : "", kCn[q]);
    fprintf(o, "], \"pPB\": [");
    for (size_t i = 0; i < pPB.size(); ++i) fprintf(o, i ? ",%.6g" : "%.6g", pPB[i]);
    fprintf(o, "], \"pX\": [");
    for (size_t i = 0; i < pX.size(); ++i) fprintf(o, i ? ",%.6g" : "%.6g", pX[i]);
    fprintf(o, "], \"sel\": {");
    for (int s = 0; s < NSS; ++s) {
        fprintf(o, "%s\"%s\": {\"_\": 0", s ? ", " : "", kSSn[s]);
        for (int q = 0; q < NC; ++q) {
            const size_t b = ((size_t)s * NC + q) * NCELL;
            fprintf(o, ", \"%s\": {\"_\": 0", kCn[q]);
            Arr(o, "s", &A.s[b], NCELL); Arr(o, "v", &A.v[b], NCELL);
            fprintf(o, "}");
        }
        fprintf(o, "}");
    }
    fprintf(o, "}}\n");
    fclose(o);
    printf("e1v: sim done (%.1f s)\n", sw.RealTime());
}

}  // namespace e1v

void e1v_verify(int mode = 0) {
    if (mode & 1) e1v::RunData();
    if (mode & 2) e1v::RunSim();
}
