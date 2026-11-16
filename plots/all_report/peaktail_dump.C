// ============================================================================
// peaktail_dump.C  --  package peaktail (the low-mass tail of the eta peak and the bump edge): row dump
// ============================================================================
// Written 2026-11-15 (prefix peaktail). A copy of the logic of edgeu_dump.C (same selections, same simulated groups and
// weights) with a WIDER M(pi+pi-) range (0.26-0.62 instead of 0.26-0.46), so that the eta peak can be studied where the
// bump cannot contribute (M(pi+pi-) > 0.4116 GeV/c2: no eta -> pi+pi-pi0 there), and with extra columns:
//   data   cor_ok == 1 rows passing "pre", 0.26 <= M(pipi) < 0.62, 0.30 <= M4 < 0.85. Columns: k (0 us, 1 pp, 2 mm), grp,
//          m4, mpipi, mee, oa, mm4, bits (1 S0, 2 pre, 4 isBest == 1, 8 PEAK = pre + dz_pi < 40 + MM4 < 2.85), day, fold,
//          ep_p, em_p.
//   sim    owned rows (mode A or B) of template T12 (classes C2-C5 and the groups G_eta, G_D, G_C), passing S0 or pre,
//          0.26 <= 0.9952 M(pipi) < 0.62, 0.29 <= M4 < 0.86. Columns as edgeu_dump.C (tree, grp, g, cls, lcls, pcls, bits,
//          m4, m4_true, mpipi, mpipi_true, mee, mee_true, w = sim w x 27.0718, f108) and the track truth: reconstructed
//          (corrected) and true momenta, polar angles and the azimuths (reco) of the four particles, the creation vertex
//          (r, z) of the two leptons, oa, oa_true.
// Nothing else is read or written; every existing file is only read.
//
// Run (from /home/przygoda/HADES/FAT/FAT; about 12 GB, about a minute):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/peaktail/aclic", kTRUE)' \
//        'plots/all_report/peaktail_dump.C+O'
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "TStopwatch.h"
#include "TSystem.h"

#include "sim_eta_physics.h"
#include "bedge_common.h"

namespace peaktaild {
using namespace bedge;
const char* const kOutP = "/hdd3/przygoda/HADES/FAT_out/all_report/peaktail";
FILE* gLog = nullptr;
#define PL(...) do { printf(__VA_ARGS__); if (gLog) fprintf(gLog, __VA_ARGS__); fflush(stdout); } while (0)

template <class T> bool WCol(const std::string& name, const std::vector<T>& v) {
    const std::string f = std::string(kOutP) + "/" + name;
    FILE* fp = fopen(f.c_str(), "wb");
    if (!fp) { PL("peaktail_dump: cannot write %s\n", f.c_str()); return false; }
    fwrite(v.data(), sizeof(T), v.size(), fp);
    fclose(fp);
    return true;
}
// copy of edgeu_dump.C GroupOf
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
}  // namespace peaktaild

void peaktail_dump() {
    using namespace peaktaild;
    TStopwatch total;
    gSystem->mkdir(kOutP, kTRUE);
    gLog = fopen(Form("%s/peaktail_dump.log", kOutP), "w");
    PL("# peaktail_dump.C, 2026-11-15: row dump of the package peaktail (eta peak tail), M(pipi) 0.26-0.62\n");
    // ---------------------------------------------------------------------------------------------------- data
    {
        std::vector<int8_t> ck, cfold; std::vector<int32_t> cgrp; std::vector<float> cm4, cpp, cee, coa, cmm4, cep, cem;
        std::vector<uint8_t> cbits; std::vector<int16_t> cday;
        scan::Data d;
        if (!scan::Load(d)) return;
        long nk[3] = {0, 0, 0};
        for (int k = 0; k < 3; ++k) {
            for (const auto& r : d.s[k]) {
                if (r.cor_ok != 1) continue;
                const bool pre = PassSel(kSelPre, r, true, r.dz_pi, r.mm4);
                if (!pre) continue;
                const bool s0 = PassSel(kSelS0, r, true, r.dz_pi, r.mm4);
                const bool pk = PassSel(kSelPeak, r, true, r.dz_pi, r.mm4);
                if (!(r.mpipi >= 0.26f && r.mpipi < 0.62f)) continue;
                if (!(r.m4 >= 0.30f && r.m4 < 0.85f)) continue;
                ck.push_back((int8_t)k); cgrp.push_back(r.grp); cm4.push_back(r.m4); cpp.push_back(r.mpipi); cee.push_back(r.mee);
                coa.push_back(r.oa); cmm4.push_back(r.mm4); cep.push_back(r.ep_p); cem.push_back(r.em_p);
                cbits.push_back((uint8_t)((s0 ? 1 : 0) | 2 | (r.isBest == 1 ? 4 : 0) | (pk ? 8 : 0)));
                cday.push_back((int16_t)r.file_id); cfold.push_back((int8_t)r.fold);
                ++nk[k];
            }
        }
        PL("A. data rows kept (pre, cor_ok == 1, 0.26 <= M(pipi) < 0.62, 0.30 <= M4 < 0.85): us %ld pp %ld mm %ld\n", nk[0], nk[1], nk[2]);
        WCol("data_k.i8", ck); WCol("data_grp.i32", cgrp); WCol("data_m4.f32", cm4); WCol("data_mpipi.f32", cpp); WCol("data_mee.f32", cee);
        WCol("data_oa.f32", coa); WCol("data_mm4.f32", cmm4); WCol("data_ep_p.f32", cep); WCol("data_em_p.f32", cem);
        WCol("data_bits.u8", cbits); WCol("data_day.i16", cday); WCol("data_fold.i8", cfold);
        FILE* fl = fopen(Form("%s/dump_layout.txt", kOutP), "w");
        if (fl) { fprintf(fl, "data %zu\n", ck.size()); fclose(fl); }
    }
    // ---------------------------------------------------------------------------------------------------- simulation
    {
        simphys::Init();
        simphys::Install();
        simphys::LoadAngles();
        simeng::SimData* D = new simeng::SimData();
        if (!simeng::Load(*D, simeng::kCompact, simeng::kMaskUS, -1e9, 1e9, 1)) return;
        const int trees[3] = {simeng::kStdUS, simeng::kLepUS, simeng::kEtaUS};
        std::vector<int8_t> ct, cg, ccls, clcls, cpcls; std::vector<int32_t> cgrp; std::vector<uint8_t> cbits;
        std::vector<float> cm4, cm4t, cpp, cppt, cee, ceet, f108, coa, coat;
        std::vector<float> tp[4][2], tth[4][2], tph[4], lvr[2], lvz[2];   // [particle ep em pip pim][reco true]
        std::vector<double> cw;
        long ng[kNG] = {0};
        for (int t = 0; t < 3; ++t) {
            const std::vector<simeng::Row>& v = D->t[trees[t]];
            for (size_t i = 0; i < v.size(); ++i) {
                const simeng::Row& r = v[i];
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
                const double mpp = kSPiPi * r.mpipi;
                if (!(mpp >= 0.26 && mpp < 0.62)) continue;
                if (!(r.m4 >= 0.29f && r.m4 < 0.86f)) continue;
                ct.push_back((int8_t)trees[t]); cgrp.push_back(r.grp); cg.push_back((int8_t)g);
                ccls.push_back((int8_t)r.cls); clcls.push_back((int8_t)r.lcls); cpcls.push_back((int8_t)r.pcls);
                cbits.push_back((uint8_t)((s0 ? 1 : 0) | (pre ? 2 : 0) | (r.isBest == 1 ? 4 : 0) | (r.ownA == 1 ? 8 : 0) | (r.ownB == 1 ? 16 : 0)));
                cm4.push_back(r.m4); cm4t.push_back(r.m4_true); cpp.push_back(r.mpipi); cppt.push_back(r.mpipi_true);
                cee.push_back(r.mee); ceet.push_back(r.mee_true); coa.push_back(r.oa); coat.push_back(r.oa_true);
                cw.push_back(r.w * kSimToRaw);
                f108.push_back((float)simphys::RowFactor(r, kPhysNom));
                const float P[4][2] = {{r.ep_p, r.ep_sim_p}, {r.em_p, r.em_sim_p}, {r.pip_p, r.pip_sim_p}, {r.pim_p, r.pim_sim_p}};
                const float TH[4][2] = {{r.ep_theta, r.ep_sim_theta}, {r.em_theta, r.em_sim_theta}, {r.pip_theta, r.pip_sim_theta}, {r.pim_theta, r.pim_sim_theta}};
                const float PH[4] = {r.ep_phi, r.em_phi, r.pip_phi, r.pim_phi};
                for (int a = 0; a < 4; ++a) {
                    for (int b = 0; b < 2; ++b) { tp[a][b].push_back(P[a][b]); tth[a][b].push_back(TH[a][b]); }
                    tph[a].push_back(PH[a]);
                }
                lvr[0].push_back(r.ep_sim_vr); lvr[1].push_back(r.em_sim_vr); lvz[0].push_back(r.ep_sim_vz); lvz[1].push_back(r.em_sim_vz);
                ++ng[g];
            }
            PL("B. sim tree %s: %zu rows read, %zu kept so far\n", simeng::TreeName(trees[t]), v.size(), ct.size());
        }
        delete D;
        PL("   kept rows per group:");
        for (int g = 0; g < kNG; ++g) PL(" %s %ld", kGTag[g], ng[g]);
        PL("\n");
        WCol("sim_tree.i8", ct); WCol("sim_grp.i32", cgrp); WCol("sim_g.i8", cg); WCol("sim_cls.i8", ccls); WCol("sim_lcls.i8", clcls);
        WCol("sim_pcls.i8", cpcls); WCol("sim_bits.u8", cbits);
        WCol("sim_m4.f32", cm4); WCol("sim_m4true.f32", cm4t); WCol("sim_mpipi.f32", cpp); WCol("sim_mpipitrue.f32", cppt);
        WCol("sim_mee.f32", cee); WCol("sim_meetrue.f32", ceet); WCol("sim_w.f64", cw); WCol("sim_f108.f32", f108);
        WCol("sim_oa.f32", coa); WCol("sim_oatrue.f32", coat);
        const char* pn[4] = {"ep", "em", "pip", "pim"};
        for (int a = 0; a < 4; ++a) {
            WCol(Form("sim_%s_p.f32", pn[a]), tp[a][0]); WCol(Form("sim_%s_ptrue.f32", pn[a]), tp[a][1]);
            WCol(Form("sim_%s_th.f32", pn[a]), tth[a][0]); WCol(Form("sim_%s_thtrue.f32", pn[a]), tth[a][1]);
            WCol(Form("sim_%s_ph.f32", pn[a]), tph[a]);
        }
        WCol("sim_ep_vr.f32", lvr[0]); WCol("sim_em_vr.f32", lvr[1]); WCol("sim_ep_vz.f32", lvz[0]); WCol("sim_em_vz.f32", lvz[1]);
        FILE* fl = fopen(Form("%s/dump_layout.txt", kOutP), "a");
        if (fl) {
            fprintf(fl, "sim %zu\n", ct.size());
            fprintf(fl, "# groups: 0 C2, 1 C3, 2 C4, 3 C5, 4 Geta, 5 GDom, 6 GDnom, 7 GCom, 8 GCnom; bits data 1 S0 2 pre 4 best 8 PEAK; sim 1 S0 2 pre 4 best 8 ownA 16 ownB\n");
            fclose(fl);
        }
    }
    PL("done in %.0f s\n", total.RealTime());
    if (gLog) fclose(gLog);
}
