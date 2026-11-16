// ============================================================================
// c2recon_rows.C  --  package c2recon: the simulated C2 rows with their true momenta (row pass)
// ============================================================================
// Written 2026-11-15 (prefix c2recon). New file; nothing else is changed. Reads the SimEngine compact (std + lep + eta,
// the same loader and selection as ecalacc_fill.C / bedge_fill.C) and, for every simulated bump row of class C2
// (eta -> pi+pi-pi0, pi0 -> e+e-gamma from HGeant, photon lost; r.tmpl == 1, r.cls == 2, owned in mode A, "pre"),
// the PAT entry of the standard set (true momenta of e+, e-, pi+, pi-, the vertex, the ECAL clusters with their GEANT
// truth). One text line per row:
//   tree idx set pat_entry grp s0 best cor_ok | W = w x RowFactor(108) x 27.0718, and W for the codes 0, 116, 100
//   | m4 m4_true mpipi mpipi_true mee mee_true | bedge bins (cls e, pp bin p, M4 bin ib) of the nominal configuration
//   | tagged (compact) m5c g_e | true momenta [MeV/c] e+ e- pi+ pi- (px py pz) | vertex (mm)
//   | found photon (cluster truth: GEANT id 1, parent pi0, the codes and the vertex of the pions): flag px py pz Ecl,
//     tagged-photon-is-this-one flag
// The QED gamma* weights are computed from these momenta in c2recon_weights.py (the photon is integrated out on its
// constraint circle, as in ecalacc; exact where the photon is known).
//
// Output: /hdd3/przygoda/HADES/FAT_out/all_report/c2recon/c2recon_rows.txt (+ c2recon_rows.log)
// Run (from /home/przygoda/HADES/FAT/FAT; about 10 GB, about 1 minute):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/c2recon/aclic", kTRUE)' \
//        'plots/all_report/c2recon_rows.C+O'
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "TFile.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TTree.h"

#include "sim_eta_physics.h"
#include "defs.h"

namespace c2rr {
using namespace alldefs;
const char* const kOutDir = "/hdd3/przygoda/HADES/FAT_out/all_report/c2recon";
const char* const kPatFile[4] = {"", "/hdd3/przygoda/HADES/files/sim/ppee_std_all.root",
                                 "/hdd3/przygoda/HADES/files/sim/ppee_lep_all.root",
                                 "/hdd3/przygoda/HADES/files/sim/ppee_eta_all.root"};
FILE* gLog = nullptr;
#define PL(...) do { printf(__VA_ARGS__); if (gLog) fprintf(gLog, __VA_ARGS__); fflush(stdout); } while (0)

// ---- tag (frozen; copy of ecalacc_fill.C) ----
const double kPi0Lo = 0.1190, kPi0Hi = 0.1730, kMeeMax = 0.100;
template <class R> inline double M5c(const R& r) {
    const double a = (double)r.m_eeg * r.m_eeg - (double)r.mee * r.mee;
    if (!(a > 0)) return -1.0;
    const double f = (scan::kMpi0 * scan::kMpi0 - (double)r.mee * r.mee) / a;
    const double m2 = (double)r.m4 * r.m4 + f * ((double)r.m5 * r.m5 - (double)r.m4 * r.m4);
    return m2 > 0 ? std::sqrt(m2) : -1.0;
}
template <class R> inline bool Tagged(const R& r) {
    return r.g_idx > 0 && r.mee < kMeeMax && r.m_eeg >= kPi0Lo && r.m_eeg < kPi0Hi;
}
// bedge binning (bedge_common.h): M4 240 x 2.5 MeV from scan::kMlo; M(pipi) 40 x 5 MeV from 0.27 (edges rounded to 1e-3)
const int kNM4 = 240, kNPP = 40;
inline double PPEdge(int k) { return std::round((0.27 + 0.005 * k) * 1000.) / 1000.; }

// ---- PAT truth (copy of the ecalacc_fill.C reader) ----
const int kNN = 5;
struct PatRow {
    float ep_px, ep_py, ep_pz, em_px, em_py, em_pz, pip_px, pip_py, pip_pz, pim_px, pim_py, pim_pz;
    float ep_p, pip_p, pip_vx, pip_vy, pip_vz;
    float pip_g0, pip_g1, pip_g2;
    float nmult;
    float pid[kNN], beta[kNN], ecl[kNN];
    float sid[kNN], sppid[kNN], sg0[kNN], sg1[kNN], sg2[kNN], spx[kNN], spy[kNN], spz[kNN], svx[kNN], svy[kNN], svz[kNN];
};
struct PatReader {
    TFile* f = nullptr; TTree* t = nullptr; PatRow r;
    bool Open(const char* fn) {
        f = TFile::Open(fn, "READ");
        if (!f || f->IsZombie()) { PL("cannot open %s\n", fn); return false; }
        t = (TTree*)f->Get("PipPimEpEm_ID");
        if (!t) { PL("no PipPimEpEm_ID in %s\n", fn); return false; }
        t->SetBranchStatus("*", 0);
        auto B = [&](const char* n, float* a) { t->SetBranchStatus(n, 1); t->SetBranchAddress(n, a); };
        B("ep_sim_px", &r.ep_px); B("ep_sim_py", &r.ep_py); B("ep_sim_pz", &r.ep_pz);
        B("em_sim_px", &r.em_px); B("em_sim_py", &r.em_py); B("em_sim_pz", &r.em_pz);
        B("pip_sim_px", &r.pip_px); B("pip_sim_py", &r.pip_py); B("pip_sim_pz", &r.pip_pz);
        B("pim_sim_px", &r.pim_px); B("pim_sim_py", &r.pim_py); B("pim_sim_pz", &r.pim_pz);
        B("ep_sim_p", &r.ep_p); B("pip_sim_p", &r.pip_p);
        B("pip_sim_vertexx", &r.pip_vx); B("pip_sim_vertexy", &r.pip_vy); B("pip_sim_vertexz", &r.pip_vz);
        B("pip_sim_geninfo", &r.pip_g0); B("pip_sim_geninfo1", &r.pip_g1); B("pip_sim_geninfo2", &r.pip_g2);
        B("neutr_mult", &r.nmult);
        for (int k = 0; k < kNN; ++k) {
            const int i = k + 1;
            B(Form("neutr_pid_%d", i), &r.pid[k]); B(Form("neutr_beta_%d", i), &r.beta[k]);
            B(Form("neutr_cluster_energy_%d", i), &r.ecl[k]);
            B(Form("neutr_sim_id_%d", i), &r.sid[k]); B(Form("neutr_sim_parentpid_%d", i), &r.sppid[k]);
            B(Form("neutr_sim_geninfo_%d", i), &r.sg0[k]); B(Form("neutr_sim_geninfo1_%d", i), &r.sg1[k]);
            B(Form("neutr_sim_geninfo2_%d", i), &r.sg2[k]);
            B(Form("neutr_sim_px_%d", i), &r.spx[k]); B(Form("neutr_sim_py_%d", i), &r.spy[k]); B(Form("neutr_sim_pz_%d", i), &r.spz[k]);
            B(Form("neutr_sim_vx_%d", i), &r.svx[k]); B(Form("neutr_sim_vy_%d", i), &r.svy[k]); B(Form("neutr_sim_vz_%d", i), &r.svz[k]);
        }
        return true;
    }
    bool Get(Long64_t e) { return t && e >= 0 && e < t->GetEntries() && t->GetEntry(e) > 0; }
};
inline int MatchPhoton(const PatRow& p) {   // copy of ecalacc_fill.C matchPhoton
    int best = -1;
    const int n = std::min(kNN, (int)std::lround(p.nmult));
    for (int k = 0; k < n; ++k) {
        if (std::lround(p.sid[k]) != 1 || std::lround(p.sppid[k]) != 7) continue;
        if (p.sg0[k] != p.pip_g0 || p.sg1[k] != p.pip_g1 || p.sg2[k] != p.pip_g2) continue;
        const double dv = std::sqrt(std::pow(p.svx[k] - p.pip_vx, 2) + std::pow(p.svy[k] - p.pip_vy, 2) + std::pow(p.svz[k] - p.pip_vz, 2));
        if (dv > 0.01) continue;
        if (best < 0 || p.ecl[k] > p.ecl[best]) best = k;
    }
    return best;
}

void Run() {
    TStopwatch sw;
    gSystem->mkdir(kOutDir, kTRUE);
    gLog = fopen(Form("%s/c2recon_rows.log", kOutDir), "w");
    PL("c2recon_rows.C -- simulated C2 rows with their true momenta, %s\n", TDatime().AsString());
    simphys::Init();
    simphys::Install();
    simphys::LoadAngles();
    simeng::SimData* D = new simeng::SimData();
    if (!simeng::Load(*D, simeng::kCompact, simeng::kMaskUS, -1e9, 1e9, 1)) return;
    const int trees[3] = {simeng::kStdUS, simeng::kLepUS, simeng::kEtaUS};
    double eD[kNPP + 1];
    for (int k = 0; k <= kNPP; ++k) eD[k] = PPEdge(k);
    struct Want { int tree; size_t idx; };
    std::vector<Want> want;
    long nSet[4] = {0, 0, 0, 0};
    for (int t = 0; t < 3; ++t) {
        const std::vector<simeng::Row>& v = D->t[trees[t]];
        for (size_t i = 0; i < v.size(); ++i) {
            const simeng::Row& r = v[i];
            if (r.ownA != 1 || r.tmpl != 1 || r.cls != 2) continue;
            if (!PassSel(kSelPre, r, r.cor_ok == 1, 0.f, r.mm4)) continue;
            int cn[6]; const int ncn = CompsOf(r, false, cn);
            bool isT12 = false; for (int q = 0; q < ncn; ++q) if (cn[q] == kT12) isT12 = true;
            if (!isT12) continue;
            want.push_back({trees[t], i}); ++nSet[r.set];
        }
    }
    PL("C2 rows (ownA, pre, T12): %zu (std %ld, lep %ld, eta %ld)\n", want.size(), nSet[1], nSet[2], nSet[3]);
    std::sort(want.begin(), want.end(), [&](const Want& a, const Want& b) {
        const simeng::Row& ra = D->t[a.tree][a.idx]; const simeng::Row& rb = D->t[b.tree][b.idx];
        if (ra.set != rb.set) return ra.set < rb.set;
        return ra.pat_entry < rb.pat_entry; });
    PatReader pr[4];
    FILE* fo = fopen(Form("%s/c2recon_rows.txt", kOutDir), "w");
    fprintf(fo, "# tree idx set pat_entry grp s0 best cor_ok W108 W0 W116 W100 m4 m4_true mpipi mpipi_true mee mee_true e p ib tag m5c g_e "
                "ep_px ep_py ep_pz em_px em_py em_pz pip_px pip_py pip_pz pim_px pim_py pim_pz vx vy vz gfound g_px g_py g_pz g_ecl g_istag\n");
    long nBad = 0, nOut = 0, nFound = 0;
    double sW = 0;
    for (const Want& q : want) {
        const simeng::Row& r = D->t[q.tree][q.idx];
        if (!pr[r.set].t && !pr[r.set].Open(kPatFile[r.set])) return;
        if (!pr[r.set].Get(r.pat_entry)) { ++nBad; continue; }
        const PatRow& p = pr[r.set].r;
        if (std::fabs(p.ep_p - r.ep_sim_p) > 1e-3f * std::max(1.f, std::fabs(r.ep_sim_p)) ||
            std::fabs(p.pip_p - r.pip_sim_p) > 1e-3f * std::max(1.f, std::fabs(r.pip_sim_p))) { ++nBad; continue; }
        const bool s0 = PassSel(kSelS0, r, r.cor_ok == 1, 0.f, r.mm4);
        const double W108 = r.w * simphys::RowFactor(r, kPhysNom) * kSimToRaw;
        const double W0 = r.w * kSimToRaw;
        const double W116 = r.w * simphys::RowFactor(r, kPhysBump1aY) * kSimToRaw;
        const double W100 = r.w * simphys::RowFactor(r, kPhysBumpPS) * kSimToRaw;
        // bedge bins of the nominal configuration (bedge_fill.C)
        const double m4p = M4sim(r.m4, r.m4_true, kS, kR);
        const double mpp = kSPiPi * r.mpipi;
        const int e = (r.mee < 0.010) ? 0 : (r.mee >= 0.010 ? 1 : 2);
        int ib = (m4p >= scan::kMlo && m4p < scan::kMhi) ? (int)((m4p - scan::kMlo) / scan::kBinW) : -1;
        if (ib >= kNM4) ib = -1;
        int pb = -1;
        if (mpp >= eD[0] && mpp < eD[kNPP]) { pb = (int)(std::upper_bound(eD, eD + kNPP + 1, mpp) - eD) - 1; if (pb >= kNPP) pb = -1; }
        const int km = MatchPhoton(p);
        const bool tag = Tagged(r);
        int isTag = 0;
        if (km >= 0 && tag) isTag = std::fabs(r.g_e - p.ecl[km]) < 1e-3f * std::max(1.f, p.ecl[km]) ? 1 : 0;
        if (km >= 0) ++nFound;
        fprintf(fo, "%d %zu %d %d %d %d %d %d %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g %d %d %d %d %.9g %.9g "
                    "%.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g %d %.9g %.9g %.9g %.9g %d\n",
                q.tree, q.idx, r.set, r.pat_entry, r.grp, s0 ? 1 : 0, r.isBest == 1 ? 1 : 0, r.cor_ok, W108, W0, W116, W100,
                r.m4, r.m4_true, r.mpipi, r.mpipi_true, r.mee, r.mee_true, e, pb, ib, tag ? 1 : 0, tag ? M5c(r) : -1., r.g_e,
                p.ep_px, p.ep_py, p.ep_pz, p.em_px, p.em_py, p.em_pz, p.pip_px, p.pip_py, p.pip_pz, p.pim_px, p.pim_py, p.pim_pz,
                p.pip_vx, p.pip_vy, p.pip_vz, km >= 0 ? 1 : 0, km >= 0 ? p.spx[km] : 0.f, km >= 0 ? p.spy[km] : 0.f,
                km >= 0 ? p.spz[km] : 0.f, km >= 0 ? p.ecl[km] : 0.f, isTag);
        ++nOut; sW += W108;
    }
    fclose(fo);
    PL("rows written: %ld (PAT read or identity failures %ld); found photons %ld; sum W108 = %.3f expected raw rows\n", nOut, nBad, nFound, sW);
    PL("done in %.1f s\n", sw.RealTime());
    delete D;
    if (gLog) fclose(gLog);
}
}  // namespace c2rr

void c2recon_rows() { c2rr::Run(); }
