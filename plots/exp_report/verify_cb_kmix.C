// ============================================================================
// verify_cb_kmix.C  --  independent check of the event-mixing k of the chapter
//                       "Combinatorial background" (prefix cb, cb_kmix.C variant A)
// ============================================================================
// HADES p+p, T = 4.53 GeV (feb22), channel pi+ pi- e+ e-. EXPERIMENTAL DATA ONLY.
// Written by the verifier of the cb chapter. It does NOT use the compact file,
// ScanEngine.h or cb_common.h: it reads the FAT output of the e+e- hypothesis
// (tree scan_nt of output_pippimepem_exp.root, corrected momenta p_cor) and the
// like-sign outputs (for the data CB used as weight), and recomputes
//   k(M4) = N+-_mix / (2 sqrt(N++_mix N--_mix))      in 10 MeV bins of M4,
// mixed combination = (pi+ pi- l1) of event E + l2 of an earlier event B of the
// same pool, both events passing pre (cor_ok == 1, isBest >= 0, vz > -160 mm;
// FAT: vz > -500 mm, PT3). Event = adjacent rows with equal (file_id, runnumber,
// event, vertex z); tracks = distinct (p, theta, phi); pion pairs as they occur
// in the rows; a lepton of E equal to a pion of the pair (same theta, phi) is
// skipped. Selections on the mixed combination: pre; peak nominal (oa_lab < 55,
// M(pipi) > 0.30); bump nominal (oa_lab < 40, 0.30 < M(pipi) < 0.41).
// Pools: P1 = 5 vertex-z classes of 20 mm from -160 mm (last open), last 10
// events (the definition of cb_kmix.C variant A); P2 = 10 mm classes, last 30
// events (pool dependence).
// Output: CB-weighted mean of k per M4 region (weights = data CB of the same
// selection, 2 sqrt(N++ N--) per 10 MeV bin) and the plain mixed count ratio.
//
// Inputs : /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/output_pippim{epem,epep,emem}_exp.root
// Outputs: printout (redirect to a log), e.g.
//          /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/cb_verify/verify_cb_kmix.log
// Run (from /home/przygoda/HADES/FAT/FAT), about 2 min:
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/aclic_cb_verify", kTRUE)' 'plots/exp_report/verify_cb_kmix.C+' > /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/cb_verify/verify_cb_kmix.log 2>&1
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <deque>
#include <vector>

#include "TFile.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TVector3.h"

namespace vkm {

const char* const kDir = "/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530";
const double kMe = 0.51099895e-3, kMpi = 0.13956995;   // GeV
const int kNw = 60; const double kLo = 0.35, kW = 0.010;
const int kNS = 3;                                       // pre, peak nominal, bump nominal

struct Trk { TLorentzVector v; float p, th, ph; };
struct PPair { TLorentzVector s; double mpp; float thp, php, thm, phm; };
struct Ev { std::vector<PPair> pr; std::vector<Trk> l[2]; };

inline TLorentzVector V(float p, float th, float ph, double m) {
    TVector3 u; u.SetMagThetaPhi(p * 1e-3, th * TMath::DegToRad(), ph * TMath::DegToRad());
    TLorentzVector v; v.SetVectM(u, m); return v;
}

struct Pool {
    int nz; double zw; int depth;
    std::vector<std::deque<Ev>> buf;
    double h[kNS][3][kNw];
    void Init(int n, double w, int d) { nz = n; zw = w; depth = d; buf.assign(nz, std::deque<Ev>()); for (int s = 0; s < kNS; ++s) for (int c = 0; c < 3; ++c) for (int j = 0; j < kNw; ++j) h[s][c][j] = 0; }
    int Cls(float vz) const { return std::max(0, std::min(nz - 1, (int)std::floor((vz + 160.) / zw))); }
    void Mix(const Ev& E, float vz) {
        std::deque<Ev>& b = buf[Cls(vz)];
        for (const Ev& B : b)
            for (const PPair& p : E.pr)
                for (int qa = 0; qa < 2; ++qa)
                    for (const Trk& a : E.l[qa]) {
                        if ((a.th == p.thp && a.ph == p.php) || (a.th == p.thm && a.ph == p.phm)) continue;
                        for (int qb = 0; qb < 2; ++qb)
                            for (const Trk& o : B.l[qb]) {
                                const int c = qa != qb ? 0 : (qa == 0 ? 1 : 2);
                                const TLorentzVector ee = a.v + o.v, f = p.s + ee;
                                const double m = f.M();
                                if (!(m >= kLo && m < kLo + kNw * kW)) continue;
                                const int j = (int)((m - kLo) / kW);
                                const double oal = p.s.Vect().Angle(ee.Vect()) * TMath::RadToDeg();
                                h[0][c][j] += 1;
                                if (oal < 55 && p.mpp > 0.30) h[1][c][j] += 1;
                                if (oal < 40 && p.mpp > 0.30 && p.mpp < 0.41) h[2][c][j] += 1;
                            }
                    }
        b.push_back(E);
        if ((int)b.size() > depth) b.pop_front();
    }
    double K(int s, int j) const { return (h[s][1][j] > 0 && h[s][2][j] > 0) ? h[s][0][j] / (2 * std::sqrt(h[s][1][j] * h[s][2][j])) : 1.0; }
};

} // namespace vkm

void verify_cb_kmix() {
    using namespace vkm;
    // ---- data CB per 10 MeV bin for the three selections (weights) ----
    double cbw[kNS][kNw];
    {
        double n[3][kNS][kNw];
        for (int k = 0; k < 3; ++k) for (int s = 0; s < kNS; ++s) for (int j = 0; j < kNw; ++j) n[k][s][j] = 0;
        const char* tag[3] = {"epem", "epep", "emem"};
        for (int k = 1; k < 3; ++k) {
            TFile* f = TFile::Open(Form("%s/output_pippim%s_exp.root", kDir, tag[k]));
            TTree* t = (TTree*)f->Get("scan_nt");
            float m4, mpp, oal, vz, isb, cor;
            t->SetBranchStatus("*", 0);
            const char* nm[] = {"m_pippimepem_cor", "m_pippim_cor", "oa_pippim_epem_lab_cor", "eVertReco_z", "isBest", "cor_ok"};
            float* ad[] = {&m4, &mpp, &oal, &vz, &isb, &cor};
            for (int i = 0; i < 6; ++i) { t->SetBranchStatus(nm[i], 1); t->SetBranchAddress(nm[i], ad[i]); }
            for (Long64_t i = 0; i < t->GetEntries(); ++i) {
                t->GetEntry(i);
                if (!(cor == 1 && isb >= 0 && vz > -160.f)) continue;
                if (!(m4 >= kLo && m4 < kLo + kNw * kW)) continue;
                const int j = (int)((m4 - kLo) / kW);
                n[k][0][j] += 1;
                if (oal < 55.f && mpp > 0.30f) n[k][1][j] += 1;
                if (oal < 40.f && mpp > 0.30f && mpp < 0.41f) n[k][2][j] += 1;
            }
            f->Close();
        }
        for (int s = 0; s < kNS; ++s) for (int j = 0; j < kNw; ++j) cbw[s][j] = 2 * std::sqrt(n[1][s][j] * n[2][s][j]);
    }

    // ---- stream the e+e- events ----
    TFile* f = TFile::Open(Form("%s/output_pippimepem_exp.root", kDir));
    TTree* t = (TTree*)f->Get("scan_nt");
    float file, run, ev, vz, isb, cor;
    float pp[4], th[4], ph[4];
    t->SetBranchStatus("*", 0);
    const char* key[] = {"file_id", "runnumber", "event", "eVertReco_z", "isBest", "cor_ok"};
    float* ka[] = {&file, &run, &ev, &vz, &isb, &cor};
    for (int i = 0; i < 6; ++i) { t->SetBranchStatus(key[i], 1); t->SetBranchAddress(key[i], ka[i]); }
    const char* pfx[4] = {"pip", "pim", "ep", "em"};
    for (int q = 0; q < 4; ++q) {
        t->SetBranchStatus(Form("%s_p_cor", pfx[q]), 1); t->SetBranchAddress(Form("%s_p_cor", pfx[q]), &pp[q]);
        t->SetBranchStatus(Form("%s_theta", pfx[q]), 1); t->SetBranchAddress(Form("%s_theta", pfx[q]), &th[q]);
        t->SetBranchStatus(Form("%s_phi", pfx[q]), 1); t->SetBranchAddress(Form("%s_phi", pfx[q]), &ph[q]);
    }
    Pool P1, P2; P1.Init(5, 20., 10); P2.Init(16, 10., 30);
    const Long64_t n = t->GetEntries();
    Ev E; bool have = false; float evz = 0;
    float cf = -1, cr = -1, ce = -1, cz = 1e9;
    long nEv = 0;
    std::vector<std::pair<int, int>> seen;
    std::vector<Trk> pips, pims;
    auto findAdd = [](std::vector<Trk>& v, float p, float t_, float f_, double m) {
        for (size_t k = 0; k < v.size(); ++k) if (v[k].p == p && v[k].th == t_ && v[k].ph == f_) return (int)k;
        Trk x; x.v = V(p, t_, f_, m); x.p = p; x.th = t_; x.ph = f_; v.push_back(x); return (int)v.size() - 1;
    };
    auto flush = [&]() {
        if (have) { P1.Mix(E, evz); P2.Mix(E, evz); ++nEv; }
        E = Ev(); have = false; seen.clear(); pips.clear(); pims.clear();
    };
    for (Long64_t i = 0; i < n; ++i) {
        t->GetEntry(i);
        if (file != cf || run != cr || ev != ce || vz != cz) { flush(); cf = file; cr = run; ce = ev; cz = vz; }
        if (!(cor == 1 && isb >= 0 && vz > -160.f)) continue;
        have = true; evz = vz;
        const int a = findAdd(pips, pp[0], th[0], ph[0], kMpi), b = findAdd(pims, pp[1], th[1], ph[1], kMpi);
        findAdd(E.l[0], pp[2], th[2], ph[2], kMe);
        findAdd(E.l[1], pp[3], th[3], ph[3], kMe);
        if (std::find(seen.begin(), seen.end(), std::make_pair(a, b)) == seen.end()) {
            seen.push_back(std::make_pair(a, b));
            PPair x; x.s = pips[a].v + pims[b].v; x.mpp = x.s.M();
            x.thp = th[0]; x.php = ph[0]; x.thm = th[1]; x.phm = ph[1];
            E.pr.push_back(x);
        }
    }
    flush();
    f->Close();
    printf("events in the pool: %ld\n", nEv);

    const double ra[] = {0.35, 0.455, 0.5225, 0.56, 0.60, 0.70, 0.35}, rb[] = {0.455, 0.5225, 0.56, 0.60, 0.70, 0.95, 0.95};
    const char* sn[kNS] = {"pre", "peak nominal", "bump nominal"};
    for (int pl = 0; pl < 2; ++pl) {
        const Pool& P = pl == 0 ? P1 : P2;
        printf("\nPool %s: CB-weighted mean k per region (bins whose centre lies in the region)\n", pl == 0 ? "P1 (5 x 20 mm, depth 10)" : "P2 (16 x 10 mm, depth 30)");
        printf("%-14s | 0.35-0.455 | 0.455-0.5225 | 0.5225-0.56 | 0.56-0.60 | 0.60-0.70 | 0.70-0.95 | 0.35-0.95 | plain count ratio 0.35-0.95\n", "selection");
        for (int s = 0; s < kNS; ++s) {
            printf("%-14s", sn[s]);
            for (int q = 0; q < 7; ++q) {
                double num = 0, den = 0;
                for (int j = 0; j < kNw; ++j) {
                    const double x = kLo + (j + 0.5) * kW;
                    if (x < ra[q] || x >= rb[q]) continue;
                    num += P.K(s, j) * cbw[s][j]; den += cbw[s][j];
                }
                printf(" | %9.4f", den > 0 ? num / den : 0.);
            }
            double c[3] = {0, 0, 0};
            for (int j = 0; j < kNw; ++j) for (int q = 0; q < 3; ++q) c[q] += P.h[s][q][j];
            printf(" | %.4f (mixed +- %.0f ++ %.0f -- %.0f, --/++ %.3f)\n", c[0] / (2 * std::sqrt(c[1] * c[2])), c[0], c[1], c[2], c[2] / c[1]);
        }
    }
}
