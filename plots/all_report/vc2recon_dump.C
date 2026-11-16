// ============================================================================
// vc2recon_dump.C  --  independent verification of package c2recon: the simulated C2 rows with their PAT truth
// ============================================================================
// Written 2026-11-15 (prefix vc2recon, verifier). New file; nothing else is changed; no c2recon file is read.
// Loads the SimEngine compact with the SAME filter and the SAME order as edgeu_dump.C (owned in mode A or B, template
// T12, S0 or pre, 0.26 <= 0.9952 M(pipi) < 0.46, 0.29 <= M4 < 0.86), so that row j of this dump is row j of the
// edgeu column files. For every row of the bedge/edgeu group C2 it reads the PAT entry (standard set) with an own
// reader and writes one text line:
//   seq (edgeu row index) tree idx set pat_entry grp bits isBest ownA ownB cls lcls pcls tmpl | W0 = w x 27.0718,
//   f108 | m4 m4_true mpipi mpipi_true mee mee_true | ep_sim_p pip_sim_p (compact, identity check) |
//   PAT truth: ep em pip pim (px py pz, MeV/c), ep/em parent id | first primary pi0 (geant_sim slot 3: id px py pz) |
//   ECAL tag fields of the compact (g_idx g_e m_eeg m5 mm4 vz) | neutral clusters k = 1..5: id parentpid px py pz
//   ecl (cluster energy) vx vy vz | pion vertex
// The QED weights, the photon solution and every fit are done in vc2recon_ana.py (own code).
//
// Output: /hdd3/przygoda/HADES/FAT_out/all_report/vc2recon/vc2recon_rows.txt (+ .log)
// Run (from /home/przygoda/HADES/FAT/FAT; about 10 GB, about 1 minute):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/vc2recon/aclic", kTRUE)' \
//        'plots/all_report/vc2recon_dump.C+O'
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
#include "bedge_common.h"

namespace vc2d {
using namespace bedge;
const char* const kOutV = "/hdd3/przygoda/HADES/FAT_out/all_report/vc2recon";
const char* const kPat = "/hdd3/przygoda/HADES/files/sim/ppee_std_all.root";
FILE* gLog = nullptr;
#define VL(...) do { printf(__VA_ARGS__); if (gLog) fprintf(gLog, __VA_ARGS__); fflush(stdout); } while (0)

// group of a row (the same switch as edgeu_dump.C / bedge_fill.C)
inline int GroupOfV(const int* cn, int ncn, bool& isT12) {
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

struct Pat {
    float l[12];          // ep em pip pim px py pz
    float epp, emp_, eppar, empar, pipp;
    float gid3, gpx3, gpy3, gpz3;
    float nmult;
    float sid[5], spp[5], spx[5], spy[5], spz[5], ecl[5], svx[5], svy[5], svz[5];
    float pvx, pvy, pvz;
};
}  // namespace vc2d

void vc2recon_dump() {
    using namespace vc2d;
    TStopwatch total;
    gSystem->mkdir(kOutV, kTRUE);
    gLog = fopen(Form("%s/vc2recon_rows.log", kOutV), "w");
    VL("# vc2recon_dump.C, %s: C2 rows of the edgeu dump order with PAT truth (verifier)\n", TDatime().AsString());
    simphys::Init();
    simphys::Install();
    simphys::LoadAngles();
    simeng::SimData* D = new simeng::SimData();
    if (!simeng::Load(*D, simeng::kCompact, simeng::kMaskUS, -1e9, 1e9, 1)) return;

    TFile* fp = TFile::Open(kPat, "READ");
    if (!fp || fp->IsZombie()) { VL("cannot open %s\n", kPat); return; }
    TTree* t = (TTree*)fp->Get("PipPimEpEm_ID");
    Pat P;
    t->SetBranchStatus("*", 0);
    auto B = [&](const char* n, float* a) {
        if (!t->GetBranch(n)) { VL("missing branch %s\n", n); return; }
        t->SetBranchStatus(n, 1); t->SetBranchAddress(n, a); };
    const char* ln[4] = {"ep", "em", "pip", "pim"};
    for (int k = 0; k < 4; ++k) {
        B(Form("%s_sim_px", ln[k]), &P.l[3 * k]); B(Form("%s_sim_py", ln[k]), &P.l[3 * k + 1]); B(Form("%s_sim_pz", ln[k]), &P.l[3 * k + 2]);
    }
    B("ep_sim_p", &P.epp); B("em_sim_p", &P.emp_); B("pip_sim_p", &P.pipp);
    B("ep_sim_parentid", &P.eppar); B("em_sim_parentid", &P.empar);
    B("geant_sim_id_3", &P.gid3); B("geant_sim_px_3", &P.gpx3); B("geant_sim_py_3", &P.gpy3); B("geant_sim_pz_3", &P.gpz3);
    B("neutr_mult", &P.nmult);
    for (int k = 0; k < 5; ++k) {
        const int i = k + 1;
        B(Form("neutr_sim_id_%d", i), &P.sid[k]); B(Form("neutr_sim_parentpid_%d", i), &P.spp[k]);
        B(Form("neutr_sim_px_%d", i), &P.spx[k]); B(Form("neutr_sim_py_%d", i), &P.spy[k]); B(Form("neutr_sim_pz_%d", i), &P.spz[k]);
        B(Form("neutr_cluster_energy_%d", i), &P.ecl[k]);
        B(Form("neutr_sim_vx_%d", i), &P.svx[k]); B(Form("neutr_sim_vy_%d", i), &P.svy[k]); B(Form("neutr_sim_vz_%d", i), &P.svz[k]);
    }
    B("pip_sim_vertexx", &P.pvx); B("pip_sim_vertexy", &P.pvy); B("pip_sim_vertexz", &P.pvz);
    const Long64_t nPat = t->GetEntries();

    FILE* fo = fopen(Form("%s/vc2recon_rows.txt", kOutV), "w");
    fprintf(fo, "# seq tree idx set pat_entry grp bits isBest ownA ownB cls lcls pcls tmpl W0 f108 m4 m4_true mpipi mpipi_true mee mee_true "
                "ep_sim_p pip_sim_p | ep_px ep_py ep_pz em_px em_py em_pz pip_px pip_py pip_pz pim_px pim_py pim_pz ep_par em_par | "
                "g3id g3px g3py g3pz | g_idx g_e m_eeg m5 mm4 vz | nmult (id ppid px py pz ecl vx vy vz) x5 | pvx pvy pvz\n");
    const int trees[3] = {simeng::kStdUS, simeng::kLepUS, simeng::kEtaUS};
    long seq = 0, nC2 = 0, nBad = 0, ng[kNG] = {0};
    for (int ti = 0; ti < 3; ++ti) {
        const std::vector<simeng::Row>& v = D->t[trees[ti]];
        for (size_t i = 0; i < v.size(); ++i) {
            const simeng::Row& r = v[i];
            if (r.ownA != 1 && r.ownB != 1) continue;
            int cn[6];
            const int ncn = CompsOf(r, false, cn);
            if (ncn == 0) continue;
            bool isT12 = false;
            const int g = GroupOfV(cn, ncn, isT12);
            if (g < 0 || !isT12) continue;
            const bool s0 = PassSel(kSelS0, r, r.cor_ok == 1, 0.f, r.mm4);
            const bool pre = PassSel(kSelPre, r, r.cor_ok == 1, 0.f, r.mm4);
            if (!s0 && !pre) continue;
            const double mpp = kSPiPi * r.mpipi;
            if (!(mpp >= 0.26 && mpp < 0.46)) continue;
            if (!(r.m4 >= 0.29f && r.m4 < 0.86f)) continue;
            const long mySeq = seq++;
            ++ng[g];
            if (g != gC2) continue;
            ++nC2;
            const int bits = (s0 ? 1 : 0) | (pre ? 2 : 0) | (r.isBest == 1 ? 4 : 0) | (r.ownA == 1 ? 8 : 0) | (r.ownB == 1 ? 16 : 0);
            if (r.set != 1 || r.pat_entry < 0 || r.pat_entry >= nPat || t->GetEntry(r.pat_entry) <= 0) {
                ++nBad; VL("row seq %ld: set %d pat_entry %d not readable\n", mySeq, r.set, r.pat_entry); continue;
            }
            if (std::fabs(P.epp - r.ep_sim_p) > 1e-3f * std::max(1.f, std::fabs(r.ep_sim_p)) ||
                std::fabs(P.pipp - r.pip_sim_p) > 1e-3f * std::max(1.f, std::fabs(r.pip_sim_p))) ++nBad;
            fprintf(fo, "%ld %d %zu %d %d %d %d %d %d %d %d %d %d %d %.10g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g |",
                    mySeq, trees[ti], i, r.set, r.pat_entry, r.grp, bits, r.isBest, r.ownA, r.ownB, r.cls, r.lcls, r.pcls, r.tmpl,
                    r.w * kSimToRaw, simphys::RowFactor(r, kPhysNom), r.m4, r.m4_true, r.mpipi, r.mpipi_true, r.mee, r.mee_true,
                    r.ep_sim_p, r.pip_sim_p);
            for (int k = 0; k < 12; ++k) fprintf(fo, " %.8g", P.l[k]);
            fprintf(fo, " %.0f %.0f | %.0f %.8g %.8g %.8g | %d %.8g %.8g %.8g %.8g %.8g | %.0f", P.eppar, P.empar, P.gid3, P.gpx3, P.gpy3, P.gpz3,
                    r.g_idx, r.g_e, r.m_eeg, r.m5, r.mm4, r.vz, P.nmult);
            for (int k = 0; k < 5; ++k)
                fprintf(fo, " %.0f %.0f %.8g %.8g %.8g %.8g %.6g %.6g %.6g", P.sid[k], P.spp[k], P.spx[k], P.spy[k], P.spz[k], P.ecl[k], P.svx[k], P.svy[k], P.svz[k]);
            fprintf(fo, " | %.6g %.6g %.6g\n", P.pvx, P.pvy, P.pvz);
        }
    }
    fclose(fo);
    VL("rows in the edgeu order: %ld (edgeu dump: 2995797); groups:", seq);
    for (int g = 0; g < kNG; ++g) VL(" %s %ld", kGTag[g], ng[g]);
    VL("\nC2 rows written %ld, PAT read/identity failures %ld; %.0f s\n", nC2, nBad, total.RealTime());
    delete D;
    if (gLog) fclose(gLog);
}
