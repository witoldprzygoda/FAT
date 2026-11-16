// ============================================================================
// veestep_pat2.C  --  INDEPENDENT verification of eestep (prefix veestep): ring and energy-loss columns of the like-sign rows
// ============================================================================
// Written 2026-11-15 by the independent verifier. Test of a sub-hypothesis that would connect the like-sign step with the
// close-pair question of D32: a conversion pair reconstructed as ONE lepton (two electrons in one ring and one track) carries
// the full photon momentum and would sit at the edge with a sharp step. Such a lepton has about twice the ring signal and
// twice the MDC energy loss. Reads, for every like-sign row of veestep_dump.C (same row order), from its PAT entry:
// rich_padnum, rich_amp, rich_avg_ringcharge, dedx_mdc, dedx_tof, mdcinnerchi2 of both lepton slots.
// Run (from /home/przygoda/HADES/FAT/FAT):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_veestep", kTRUE)' \
//        'plots/all_report/veestep_pat2.C+O'
// ============================================================================
#include <cstdint>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include "TFile.h"
#include "TNtuple.h"
#include "TStopwatch.h"
#include "TString.h"

void veestep_pat2() {
    const char* out = "/hdd3/przygoda/HADES/FAT_out/all_report/veestep";
    TStopwatch sw;
    long n = 0;
    { FILE* f = fopen(Form("%s/veestep_layout.txt", out), "r"); char s[32]; if (f) { if (fscanf(f, "%31s %ld", s, &n) != 2) n = 0; fclose(f); } }
    std::vector<int8_t> k(n); std::vector<int16_t> day(n); std::vector<int32_t> ent(n);
    { FILE* f = fopen(Form("%s/d_k.i8", out), "rb"); if (fread(k.data(), 1, n, f) != (size_t)n) printf("short k\n"); fclose(f); }
    { FILE* f = fopen(Form("%s/d_day.i16", out), "rb"); if (fread(day.data(), 2, n, f) != (size_t)n) printf("short day\n"); fclose(f); }
    { FILE* f = fopen(Form("%s/d_ent.i32", out), "rb"); if (fread(ent.data(), 4, n, f) != (size_t)n) printf("short ent\n"); fclose(f); }
    const int NV = 6;
    const char* vn[NV] = {"rich_padnum", "rich_amp", "rich_avg_ringcharge", "dedx_mdc", "dedx_tof", "mdcinnerchi2"};
    std::vector<float> col[2][NV];
    for (int s = 0; s < 2; ++s) for (int j = 0; j < NV; ++j) col[s][j].assign(n, -999.f);
    std::map<int, std::vector<long>> byday;
    for (long i = 0; i < n; ++i) if (k[i] > 0) byday[day[i]].push_back(i);
    for (auto& kv : byday) {
        TFile* f = TFile::Open(Form("/hdd3/przygoda/HADES/files/exp/gen4/lep_%03d.root", kv.first), "READ");
        for (int kk = 1; kk <= 2; ++kk) {
            TNtuple* t = (TNtuple*)f->Get(kk == 1 ? "PipPimEpEp_ID" : "PipPimEmEm_ID");
            const char* p[2] = {kk == 1 ? "ep1" : "em1", kk == 1 ? "ep2" : "em2"};
            t->SetBranchStatus("*", 0);
            float v[2][NV];
            for (int s = 0; s < 2; ++s)
                for (int j = 0; j < NV; ++j) {
                    TString b = Form("%s_%s", p[s], vn[j]);
                    t->SetBranchStatus(b, 1); t->SetBranchAddress(b, &v[s][j]);
                }
            for (long i : kv.second) {
                if (k[i] != kk) continue;
                t->GetEntry(ent[i]);
                for (int s = 0; s < 2; ++s) for (int j = 0; j < NV; ++j) col[s][j][i] = v[s][j];
            }
        }
        f->Close(); delete f;
    }
    for (int s = 0; s < 2; ++s)
        for (int j = 0; j < NV; ++j) {
            FILE* fo = fopen(Form("%s/d2_%s_%d.f32", out, vn[j], s + 1), "wb");
            fwrite(col[s][j].data(), 4, n, fo); fclose(fo);
        }
    printf("veestep_pat2: %ld rows, like-sign columns written, %.0f s\n", n, sw.RealTime());
}
