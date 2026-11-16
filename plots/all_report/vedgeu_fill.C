// ============================================================================
// vedgeu_fill.C  --  INDEPENDENT VERIFIER of the package edgeu: own row pass
// ============================================================================
// Written 2026-11-15 (prefix vedgeu). Does NOT use edgeu_dump.C, edgeu_core.py or any edgeu column file.
// Own tree loops: data with own SetBranchAddress on the compact data trees us/pp/mm; simulation with
// simeng::BindRow over std_us, lep_us, eta_us. Own selection code, written from defs.h Sect. D:
//   pre = isBest >= 0 && vz > -500 && vz > -160 && cor_ok == 1 (&& PT3: all data rows; sim (trigbit & 8192) != 0)
//   S0  = pre && pip_isUsedVertex == 1 && pim_isUsedVertex == 1
// Frozen definitions used on purpose: simphys::RowFactor(row, 108) (physics weights V3), alldefs::CompsOf (class of
// a row), the conversion sum w -> expected raw data rows L / L_sim x 63/64.
//
// Output (/hdd3/przygoda/HADES/FAT_out/all_report/vedgeu/), flat little-endian columns (no header):
//   d_k.i1 (0 us, 1 pp, 2 mm) d_grp.i4 d_fid.i2 d_m4.f4 d_pp.f4 d_ee.f4 d_bits.u1 (1 S0, 2 pre, 4 isBest == 1)
//     rows: pre, 0.25 <= M(pipi) < 0.47, 0.30 <= M4 < 0.90 (data scale, raw columns)
//   s_tree.i1 (0 std, 1 lep, 2 eta) s_grp.i4 s_c.i1 (0 C2 1 C3 2 C4 3 C5 4 G_eta 5 G_D 6 G_C) s_om.i1 s_bits.u1 (1 S0, 2 pre)
//   s_m4.f4 s_m4t.f4 s_pp.f4 s_ppt.f4 s_ee.f4 s_eet.f4 s_w.f8 (= w x RowFactor(108) x 27.0718: expected raw rows)
//   s_w0.f8 (= w x 27.0718, physics weights off)
//     rows: ownership mode A, pre, class in C2-C5 / G_eta / G_D / G_C, 0.30 <= M4 < 0.90, 0.25 <= 0.9952 M(pipi) < 0.47
//
// Run (from /home/przygoda/HADES/FAT/FAT):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/vedgeu/aclic", kTRUE)' \
//        'plots/all_report/vedgeu_fill.C+O'
// ============================================================================
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "TFile.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TTree.h"

#include "sim_eta_physics.h"
#include "defs.h"

namespace vef {
const char* const kOut = "/hdd3/przygoda/HADES/FAT_out/all_report/vedgeu";
template <class T> void W(const char* name, const std::vector<T>& v) {
    std::string fn = std::string(kOut) + "/" + name;
    FILE* f = fopen(fn.c_str(), "wb");
    if (!f) { printf("cannot write %s\n", fn.c_str()); return; }
    fwrite(v.data(), sizeof(T), v.size(), f);
    fclose(f);
}
}  // namespace vef

void vedgeu_fill() {
    using namespace vef;
    TStopwatch tot;
    gSystem->mkdir(kOut, kTRUE);
    FILE* lg = fopen(Form("%s/vedgeu_fill.log", kOut), "w");
    // ------------------------------------------------------------------------------------------ data
    {
        std::vector<int8_t> ck; std::vector<int32_t> cg; std::vector<int16_t> cf; std::vector<float> c4, cp, ce; std::vector<uint8_t> cb;
        TFile* f = TFile::Open(scan::kCompact, "READ");
        if (!f || f->IsZombie()) { printf("cannot open data\n"); return; }
        const char* tn[3] = {"us", "pp", "mm"};
        long nk[3] = {0, 0, 0}, ns0[3] = {0, 0, 0};
        for (int s = 0; s < 3; ++s) {
            TTree* t = (TTree*)f->Get(tn[s]);
            Int_t isBest, cor_ok, pipv, pimv, grp, fid;
            Float_t vz, m4, mpipi, mee;
            t->SetBranchStatus("*", 0);
            for (const char* b : {"isBest", "cor_ok", "pip_isUsedVertex", "pim_isUsedVertex", "grp", "vz", "m4", "mpipi", "mee", "file_id"}) t->SetBranchStatus(b, 1);
            t->SetBranchAddress("isBest", &isBest); t->SetBranchAddress("cor_ok", &cor_ok);
            t->SetBranchAddress("pip_isUsedVertex", &pipv); t->SetBranchAddress("pim_isUsedVertex", &pimv);
            t->SetBranchAddress("grp", &grp); t->SetBranchAddress("vz", &vz); t->SetBranchAddress("m4", &m4);
            t->SetBranchAddress("mpipi", &mpipi); t->SetBranchAddress("mee", &mee); t->SetBranchAddress("file_id", &fid);
            const Long64_t n = t->GetEntries();
            for (Long64_t i = 0; i < n; ++i) {
                t->GetEntry(i);
                const bool pre = isBest >= 0 && vz > -500.f && vz > -160.f && cor_ok == 1;
                if (!pre) continue;
                if (!(mpipi >= 0.25f && mpipi < 0.47f && m4 >= 0.30f && m4 < 0.90f)) continue;
                const bool s0 = pipv == 1 && pimv == 1;
                ck.push_back((int8_t)s); cg.push_back(grp); cf.push_back((int16_t)fid); c4.push_back(m4); cp.push_back(mpipi); ce.push_back(mee);
                cb.push_back((uint8_t)((s0 ? 1 : 0) | 2 | (isBest == 1 ? 4 : 0)));
                ++nk[s]; if (s0) ++ns0[s];
            }
            t->ResetBranchAddresses();
        }
        f->Close();
        fprintf(lg, "data rows (pre, 0.25 <= M(pipi) < 0.47, 0.30 <= M4 < 0.90): us %ld pp %ld mm %ld; S0: %ld %ld %ld\n", nk[0], nk[1], nk[2], ns0[0], ns0[1], ns0[2]);
        printf("data rows: us %ld pp %ld mm %ld; S0 %ld %ld %ld (%.0f s)\n", nk[0], nk[1], nk[2], ns0[0], ns0[1], ns0[2], tot.RealTime()); tot.Continue();
        W("d_k.i1", ck); W("d_grp.i4", cg); W("d_fid.i2", cf); W("d_m4.f4", c4); W("d_pp.f4", cp); W("d_ee.f4", ce); W("d_bits.u1", cb);
        fprintf(lg, "data %zu\n", ck.size());
    }
    // ------------------------------------------------------------------------------------------ simulation
    {
        simphys::Init(); simphys::Install(); simphys::LoadAngles();
        const double toRaw = alldefs::kL / alldefs::kLsim / (64.0 / 63.0);
        std::vector<int8_t> ct, cc, co; std::vector<int32_t> cg; std::vector<uint8_t> cb;
        std::vector<float> c4, c4t, cp, cpt, ce, cet; std::vector<double> cw, cw0;
        struct Src { const char* file; const char* tree; int id; };
        const Src src[3] = {{simeng::kCompact, "std_us", 0}, {simeng::kCompact, "lep_us", 1}, {simeng::kCompactEta, "eta_us", 2}};
        double sumC[7] = {0, 0, 0, 0, 0, 0, 0};
        long nC[7] = {0, 0, 0, 0, 0, 0, 0};
        for (int q = 0; q < 3; ++q) {
            TFile* f = TFile::Open(src[q].file, "READ");
            TTree* t = (TTree*)f->Get(src[q].tree);
            simeng::Row r;
            simeng::BindRow(t, r, false);
            const Long64_t n = t->GetEntries();
            for (Long64_t i = 0; i < n; ++i) {
                t->GetEntry(i);
                if (r.ownA != 1) continue;
                const bool pre = r.isBest >= 0 && r.vz > -500.f && (r.trigbit & 8192) != 0 && r.vz > -160.f && r.cor_ok == 1;
                if (!pre) continue;
                const bool s0 = r.pip_isUsedVertex == 1 && r.pim_isUsedVertex == 1;
                if (!(r.m4 >= 0.30f && r.m4 < 0.90f)) continue;
                const double mppp = 0.9952 * r.mpipi;
                if (!(mppp >= 0.25 && mppp < 0.47)) continue;
                int cn[6];
                const int ncn = alldefs::CompsOf(r, false, cn);
                int cl = -1; bool t12 = false;
                for (int k = 0; k < ncn; ++k) {
                    const int c = cn[k];
                    if (c == alldefs::kT12) t12 = true;
                    if (c == alldefs::kC2) cl = 0;
                    else if (c == alldefs::kC3) cl = 1;
                    else if (c == alldefs::kC4) cl = 2;
                    else if (c == alldefs::kC5) cl = 3;
                    else if (c == alldefs::kGeta) cl = 4;
                    else if (cl < 0 && c == alldefs::kGD) cl = 5;
                    else if (cl < 0 && c == alldefs::kGconv) cl = 6;
                }
                if (!t12 || cl < 0) continue;
                const double w0 = r.w * toRaw;
                const double w = w0 * simphys::RowFactor(r, 108);
                if (w == 0. && w0 == 0.) continue;
                ct.push_back((int8_t)src[q].id); cg.push_back(r.grp); cc.push_back((int8_t)cl);
                co.push_back((int8_t)((r.pcls == 3 || r.pcls == 4) ? 1 : 0));
                cb.push_back((uint8_t)((s0 ? 1 : 0) | 2));
                c4.push_back(r.m4); c4t.push_back(r.m4_true); cp.push_back(r.mpipi); cpt.push_back(r.mpipi_true);
                ce.push_back(r.mee); cet.push_back(r.mee_true); cw.push_back(w); cw0.push_back(w0);
                if (s0) { sumC[cl] += w; ++nC[cl]; }
            }
            t->ResetBranchAddresses();
            f->Close();
            printf("sim tree %s done, rows kept so far %zu (%.0f s)\n", src[q].tree, ct.size(), tot.RealTime()); tot.Continue();
        }
        const char* nm[7] = {"C2", "C3", "C4", "C5", "Geta", "GD", "GC"};
        for (int k = 0; k < 7; ++k) fprintf(lg, "sim S0 %s: rows %ld, expected raw rows %.2f\n", nm[k], nC[k], sumC[k]);
        W("s_tree.i1", ct); W("s_grp.i4", cg); W("s_c.i1", cc); W("s_om.i1", co); W("s_bits.u1", cb);
        W("s_m4.f4", c4); W("s_m4t.f4", c4t); W("s_pp.f4", cp); W("s_ppt.f4", cpt); W("s_ee.f4", ce); W("s_eet.f4", cet);
        W("s_w.f8", cw); W("s_w0.f8", cw0);
        fprintf(lg, "sim %zu\n", ct.size());
    }
    fprintf(lg, "done in %.0f s\n", tot.RealTime());
    fclose(lg);
    printf("done in %.0f s\n", tot.RealTime());
}
