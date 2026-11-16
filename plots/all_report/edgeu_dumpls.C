// ============================================================================
// edgeu_dumpls.C  --  package edgeu: dump of the SIMULATED like-sign rows (e+e+ and e-e- trees) for the like-sign null test
// ============================================================================
// Written 2026-11-15 (prefix edgeu). The data like-sign rows show a step at u = 0 in e-e- (not in e+e+). This macro dumps
// the simulated like-sign rows of the three sets (std, lep, eta; owned rows, S0 or pre) with their truth (template, lepton
// class, GEANT ids and parents of the two lepton slots, ghost bits, flags) so that edgeu_lsnull.py can find which true
// particle combinations make an edge-like structure in the like-sign sample. The simulation has no combinatorial
// background (rule of the simulation study): these rows are used ONLY to identify the origin of the like-sign step.
// Columns (flat files in /hdd3/przygoda/HADES/FAT_out/all_report/edgeu/, prefix simls_): tree, grp, tmpl, cls, lcls, pcls,
// flags, ghost, bits (1 S0, 2 pre, 8 ownA, 16 ownB), id1, id2, par1, par2, gpar1, gpar2, proc1, proc2, m4, m4t, mpipi,
// mpipit, mee, meet, w (x 27.0718), f108.
// Run (from /home/przygoda/HADES/FAT/FAT; about 8 GB, 1 min):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_edgeu", kTRUE)' \
//        'plots/all_report/edgeu_dumpls.C+O'
// ============================================================================
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "TStopwatch.h"
#include "TSystem.h"

#include "sim_eta_physics.h"
#include "bedge_common.h"

namespace edgeuls {
using namespace bedge;
const char* const kOutU = "/hdd3/przygoda/HADES/FAT_out/all_report/edgeu";
template <class T> void WCol(const std::string& name, const std::vector<T>& v) {
    FILE* fp = fopen((std::string(kOutU) + "/" + name).c_str(), "wb");
    if (!fp) { printf("cannot write %s\n", name.c_str()); return; }
    fwrite(v.data(), sizeof(T), v.size(), fp);
    fclose(fp);
}
}  // namespace edgeuls

void edgeu_dumpls() {
    using namespace edgeuls;
    TStopwatch sw;
    simphys::Init();
    simphys::Install();
    simphys::LoadAngles();
    simeng::SimData* D = new simeng::SimData();
    const int mask = (1 << simeng::kStdPP) | (1 << simeng::kStdMM) | (1 << simeng::kLepPP) | (1 << simeng::kLepMM) |
                     (1 << simeng::kEtaPP) | (1 << simeng::kEtaMM);
    if (!simeng::Load(*D, simeng::kCompact, mask, -1e9, 1e9, 1)) return;
    const int trees[6] = {simeng::kStdPP, simeng::kStdMM, simeng::kLepPP, simeng::kLepMM, simeng::kEtaPP, simeng::kEtaMM};
    std::vector<int8_t> ct, ctm, ccls, clcls, cpcls; std::vector<int32_t> cgrp, cfl, cgh, cid1, cid2, cpa1, cpa2, cgp1, cgp2, cpr1, cpr2;
    std::vector<uint8_t> cbits; std::vector<float> cm4, cm4t, cpp, cppt, cee, ceet, cf; std::vector<double> cw;
    for (int t = 0; t < 6; ++t) {
        const std::vector<simeng::Row>& v = D->t[trees[t]];
        for (const auto& r : v) {
            if (r.ownA != 1 && r.ownB != 1) continue;
            const bool s0 = PassSel(kSelS0, r, r.cor_ok == 1, 0.f, r.mm4);
            const bool pre = PassSel(kSelPre, r, r.cor_ok == 1, 0.f, r.mm4);
            if (!s0 && !pre) continue;
            const double mpp = kSPiPi * r.mpipi;
            if (!(mpp >= 0.26 && mpp < 0.46)) continue;
            if (!(r.m4 >= 0.29f && r.m4 < 0.86f)) continue;
            ct.push_back((int8_t)trees[t]); cgrp.push_back(r.grp); ctm.push_back((int8_t)r.tmpl); ccls.push_back((int8_t)r.cls);
            clcls.push_back((int8_t)r.lcls); cpcls.push_back((int8_t)r.pcls); cfl.push_back(r.flags); cgh.push_back(r.ghost);
            cbits.push_back((uint8_t)((s0 ? 1 : 0) | (pre ? 2 : 0) | (r.isBest == 1 ? 4 : 0) | (r.ownA == 1 ? 8 : 0) | (r.ownB == 1 ? 16 : 0)));
            cid1.push_back(r.ep_sim_id); cid2.push_back(r.em_sim_id); cpa1.push_back(r.ep_sim_parentid); cpa2.push_back(r.em_sim_parentid);
            cgp1.push_back(r.ep_sim_grandparentid); cgp2.push_back(r.em_sim_grandparentid); cpr1.push_back(r.ep_sim_processid); cpr2.push_back(r.em_sim_processid);
            cm4.push_back(r.m4); cm4t.push_back(r.m4_true); cpp.push_back(r.mpipi); cppt.push_back(r.mpipi_true); cee.push_back(r.mee); ceet.push_back(r.mee_true);
            cw.push_back(r.w * kSimToRaw); cf.push_back((float)simphys::RowFactor(r, kPhysNom));
        }
        printf("tree %s: %zu rows read, %zu kept so far\n", simeng::TreeName(trees[t]), v.size(), ct.size());
    }
    delete D;
    WCol("simls_tree.i8", ct); WCol("simls_grp.i32", cgrp); WCol("simls_tmpl.i8", ctm); WCol("simls_cls.i8", ccls); WCol("simls_lcls.i8", clcls);
    WCol("simls_pcls.i8", cpcls); WCol("simls_flags.i32", cfl); WCol("simls_ghost.i32", cgh); WCol("simls_bits.u8", cbits);
    WCol("simls_id1.i32", cid1); WCol("simls_id2.i32", cid2); WCol("simls_par1.i32", cpa1); WCol("simls_par2.i32", cpa2);
    WCol("simls_gpar1.i32", cgp1); WCol("simls_gpar2.i32", cgp2); WCol("simls_proc1.i32", cpr1); WCol("simls_proc2.i32", cpr2);
    WCol("simls_m4.f32", cm4); WCol("simls_m4true.f32", cm4t); WCol("simls_mpipi.f32", cpp); WCol("simls_mpipitrue.f32", cppt);
    WCol("simls_mee.f32", cee); WCol("simls_meetrue.f32", ceet); WCol("simls_w.f64", cw); WCol("simls_f108.f32", cf);
    FILE* fl = fopen(Form("%s/dumpls_layout.txt", kOutU), "w");
    if (fl) { fprintf(fl, "simls %zu\n", ct.size()); fclose(fl); }
    printf("done: %zu rows, %.0f s\n", ct.size(), sw.RealTime());
}
