// ============================================================================
// vecalacc_rows.C -- INDEPENDENT VERIFIER of package ecalacc (2026-11-15): own row dump
// ============================================================================
// Selects, with the frozen definitions (defs.h, SimEngine.h, sim_eta_physics.h; NOT ecalacc_fill.C), the owned
// simulated pi+pi-e+e- rows of the classes C2, C3 (bump) and C4, C5 (peak) in "pre" (template 1), and writes one text
// line per row with the row keys, the weight w x RowFactor(108) x kSimToRaw, the reconstructed masses, the compact's
// ECAL-photon columns, and, read from the PAT file of the set (entry pat_entry), the TRUE momenta of the four charged
// tracks, the pion vertex and the five neutral candidates with their cluster and GEANT truth.
// A second output lists every neutral candidate of every entry of the standard-set PAT file (GEANT id, parent, true
// momentum and creation vertex, cluster theta/phi/sector/cell/energy), for an empirical map of the ECAL acceptance.
// Outputs: /hdd3/przygoda/HADES/FAT_out/all_report/vecalacc/vrows.txt, vclusters_std.txt
// Run (from /home/przygoda/HADES/FAT/FAT):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_vecalacc", kTRUE)' \
//        'plots/all_report/vecalacc_rows.C+O'
// ============================================================================
#include <cstdio>
#include <cmath>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "sim_eta_physics.h"
#include "defs.h"

namespace vecalacc {
using namespace alldefs;
const char* const kOut = "/hdd3/przygoda/HADES/FAT_out/all_report/vecalacc";
const char* const kPat[4] = {"", "/hdd3/przygoda/HADES/files/sim/ppee_std_all.root", "/hdd3/przygoda/HADES/files/sim/ppee_lep_all.root",
                             "/hdd3/przygoda/HADES/files/sim/ppee_eta_all.root"};
const int NN = 5;
const char* const kScal[] = {"ep_sim_px", "ep_sim_py", "ep_sim_pz", "em_sim_px", "em_sim_py", "em_sim_pz", "pip_sim_px", "pip_sim_py",
                             "pip_sim_pz", "pim_sim_px", "pim_sim_py", "pim_sim_pz", "pip_sim_vertexx", "pip_sim_vertexy",
                             "pip_sim_vertexz", "ep_sim_vertexx", "ep_sim_vertexy", "ep_sim_vertexz", "pip_sim_geninfo",
                             "pip_sim_geninfo1", "pip_sim_geninfo2", "ep_sim_p", "pip_sim_p", "neutr_mult"};
const int NS = sizeof(kScal) / sizeof(kScal[0]);
const char* const kNeu[] = {"neutr_pid_%d", "neutr_beta_%d", "neutr_cluster_energy_%d", "neutr_cluster_theta_%d", "neutr_cluster_phi_%d",
                            "neutr_cluster_sector_%d", "neutr_cluster_cell_%d", "neutr_sim_id_%d", "neutr_sim_parentpid_%d",
                            "neutr_sim_geninfo_%d", "neutr_sim_geninfo1_%d", "neutr_sim_geninfo2_%d", "neutr_sim_px_%d",
                            "neutr_sim_py_%d", "neutr_sim_pz_%d", "neutr_sim_vx_%d", "neutr_sim_vy_%d", "neutr_sim_vz_%d"};
const int NV = sizeof(kNeu) / sizeof(kNeu[0]);
struct Pat {
    TFile* f = nullptr; TTree* t = nullptr; float s[NS]; float n[NV][NN];
    bool Open(const char* fn) {
        f = TFile::Open(fn, "READ"); if (!f || f->IsZombie()) return false;
        t = (TTree*)f->Get("PipPimEpEm_ID"); if (!t) return false;
        t->SetBranchStatus("*", 0);
        for (int i = 0; i < NS; ++i) { t->SetBranchStatus(kScal[i], 1); t->SetBranchAddress(kScal[i], &s[i]); }
        for (int j = 0; j < NV; ++j) for (int k = 0; k < NN; ++k) {
            TString b = Form(kNeu[j], k + 1); t->SetBranchStatus(b, 1); t->SetBranchAddress(b, &n[j][k]); }
        return true;
    }
};

void Run(int doClusters) {
    gSystem->mkdir(kOut, kTRUE);
    simphys::Init(); simphys::Install(); simphys::LoadAngles();
    simeng::SimData* D = new simeng::SimData();
    if (!simeng::Load(*D, simeng::kCompact, simeng::kMaskUS, -1e9, 1e9, 1)) return;
    const int trees[3] = {simeng::kStdUS, simeng::kLepUS, simeng::kEtaUS};
    Pat pat[4];
    FILE* fo = fopen(Form("%s/vrows.txt", kOut), "w");
    fprintf(fo, "# set pat_entry cls s0 w m4 m4true mpipi mee | g_idx g_e g_theta g_phi m_eeg m5 | %d PAT scalars | 5 x %d neutral columns (column-major per variable) | reconstructed p theta phi of e+ e- pi+ pi-\n", NS, NV);
    long n[6] = {0, 0, 0, 0, 0, 0}, nbad = 0; double sw[6] = {0, 0, 0, 0, 0, 0};
    for (int it = 0; it < 3; ++it) {
        for (const simeng::Row& r : D->t[trees[it]]) {
            if (r.ownA != 1 || r.tmpl != 1 || r.cls < 2 || r.cls > 5) continue;
            if (!PassSel(kSelPre, r, r.cor_ok == 1, 0.f, r.mm4)) continue;
            const double w = r.w * simphys::RowFactor(r, kPhysNom) * kSimToRaw;
            if (w == 0.) continue;
            const int s0 = PassSel(kSelS0, r, r.cor_ok == 1, 0.f, r.mm4) ? 1 : 0;
            Pat& P = pat[r.set];
            if (!P.t && !P.Open(kPat[r.set])) { printf("cannot open %s\n", kPat[r.set]); return; }
            if (P.t->GetEntry(r.pat_entry) <= 0) { ++nbad; continue; }
            if (std::fabs(P.s[21] - r.ep_sim_p) > 1e-3 * std::max(1.f, std::fabs(r.ep_sim_p)) ||
                std::fabs(P.s[22] - r.pip_sim_p) > 1e-3 * std::max(1.f, std::fabs(r.pip_sim_p))) { ++nbad; continue; }
            ++n[r.cls]; sw[r.cls] += w;
            fprintf(fo, "%d %d %d %d %.8g %.7g %.7g %.7g %.7g | %d %.6g %.6g %.6g %.7g %.7g |", r.set, r.pat_entry, r.cls, s0, w, r.m4,
                    r.m4_true, r.mpipi, r.mee, r.g_idx, r.g_e, r.g_theta, r.g_phi, r.m_eeg, r.m5);
            for (int i = 0; i < NS; ++i) fprintf(fo, " %.7g", P.s[i]);
            fprintf(fo, " |");
            for (int j = 0; j < NV; ++j) for (int k = 0; k < NN; ++k) fprintf(fo, " %.7g", P.n[j][k]);
            // reconstructed (corrected) track kinematics of the compact row: p [MeV/c], theta, phi [deg] of e+ e- pi+ pi-
            fprintf(fo, " | %.7g %.7g %.7g %.7g %.7g %.7g %.7g %.7g %.7g %.7g %.7g %.7g", r.ep_p, r.ep_theta, r.ep_phi, r.em_p, r.em_theta,
                    r.em_phi, r.pip_p, r.pip_theta, r.pip_phi, r.pim_p, r.pim_theta, r.pim_phi);
            fprintf(fo, "\n");
        }
    }
    fclose(fo);
    for (int c = 2; c <= 5; ++c) printf("class C%d: %ld rows, sum of weights x kSimToRaw %.2f\n", c, n[c], sw[c]);
    printf("rows dropped (PAT read or identity check): %ld\n", nbad);
    if (!doClusters) return;
    // all neutral candidates of the standard-set PAT file
    Pat Q; if (!Q.Open(kPat[1])) return;
    FILE* fc = fopen(Form("%s/vclusters_std.txt", kOut), "w");
    fprintf(fc, "# entry k pid beta E_cl theta_cl phi_cl sector cell sim_id parentpid px py pz vx vy vz (one line per candidate)\n");
    long nc = 0; float last[6] = {-1, -1, -1, -1, -1, -1};
    for (Long64_t e = 0; e < Q.t->GetEntries(); ++e) {
        Q.t->GetEntry(e);
        // the rows of one event repeat the same clusters: write an event once (same pion-vertex + neutral multiplicity + first cluster)
        const float key[6] = {Q.s[12], Q.s[13], Q.s[14], Q.s[23], Q.n[2][0], Q.n[3][0]};
        bool same = true; for (int i = 0; i < 6; ++i) same = same && key[i] == last[i];
        if (same) continue;
        for (int i = 0; i < 6; ++i) last[i] = key[i];
        const int m = std::min(NN, (int)std::lround(Q.s[23]));
        for (int k = 0; k < m; ++k) {
            fprintf(fc, "%lld %d", e, k);
            const int cols[16] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 15, 16, 17, -1};
            for (int q = 0; q < 15; ++q) fprintf(fc, " %.6g", Q.n[cols[q]][k]);
            fprintf(fc, "\n"); ++nc;
        }
    }
    fclose(fc);
    printf("neutral candidates written: %ld\n", nc);
}
}  // namespace vecalacc

void vecalacc_rows(int doClusters = 0) { vecalacc::Run(doClusters); }
