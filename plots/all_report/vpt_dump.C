// ============================================================================
// vpt_dump.C  --  INDEPENDENT VERIFIER of the package peaktail (low-mass tail of the eta peak): own row dump
// ============================================================================
// Written 2026-11-15 (prefix vpt). No peaktail_* code or output is used. Frozen headers only (defs.h selections,
// CompsOf, simphys::RowFactor). Writes flat column files to /hdd3/przygoda/HADES/FAT_out/all_report/vpt/.
//   data  cor_ok == 1 rows passing "pre", 0.26 <= M(pipi) < 0.70, 0.30 <= M4 < 0.85: k, grp, m4, mpipi, mee, bits (1 S0, 2 pre)
//   sim   owned rows (mode A or B) of T12 (groups C2-C5, G_eta, G_D, G_C), passing S0 or pre, all M(pipi) (no cut),
//         0.29 <= M4 < 0.86: tree, grp, g, bits (1 S0, 2 pre, 8 ownA, 16 ownB), m4, m4_true, mpipi, mpipi_true, mee,
//         mee_true, w (x kSimToRaw), f108 (V3), momenta reco (corrected and uncorrected) and true, polar angles reco and
//         true, azimuth reco, lepton creation vertex r, z, lepton GEANT process id and parent id.
// Run (from /home/przygoda/HADES/FAT/FAT):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/vpt/aclic", kTRUE)' 'plots/all_report/vpt_dump.C+O'
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

namespace vptd {
using namespace bedge;
const char* const kOutV = "/hdd3/przygoda/HADES/FAT_out/all_report/vpt";
template <class T> void W(const std::string& name, const std::vector<T>& v) {
    const std::string f = std::string(kOutV) + "/" + name;
    FILE* fp = fopen(f.c_str(), "wb");
    if (!fp) { printf("cannot write %s\n", f.c_str()); return; }
    fwrite(v.data(), sizeof(T), v.size(), fp);
    fclose(fp);
}
inline int Grp(const int* cn, int ncn, bool& t12) {
    t12 = false; int g = -1;
    for (int q = 0; q < ncn; ++q) {
        const int c = cn[q];
        if (c == kT12) t12 = true;
        if (g >= 0) continue;
        if (c == kC2) g = 0; else if (c == kC3) g = 1; else if (c == kC4) g = 2; else if (c == kC5) g = 3;
        else if (c == kGeta) g = 4;
        else if (c == kGpi0Dom || c == kGvecom || c == kGothom) g = 5;
        else if (c == kGpi0Dnom || c == kGvecnom || c == kGothnom) g = 6;
        else if (c == kGconvom) g = 7;
        else if (c == kGconvnom) g = 8;
    }
    return g;
}
}  // namespace vptd

void vpt_dump() {
    using namespace vptd;
    TStopwatch sw;
    gSystem->mkdir(kOutV, kTRUE);
    FILE* lay = fopen(Form("%s/layout.txt", kOutV), "w");
    {
        std::vector<int8_t> ck; std::vector<int32_t> cg; std::vector<float> cm4, cpp, cee; std::vector<uint8_t> cb;
        scan::Data d;
        if (!scan::Load(d)) return;
        for (int k = 0; k < 3; ++k)
            for (const auto& r : d.s[k]) {
                if (r.cor_ok != 1) continue;
                if (!PassSel(kSelPre, r, true, r.dz_pi, r.mm4)) continue;
                const bool s0 = PassSel(kSelS0, r, true, r.dz_pi, r.mm4);
                if (!(r.mpipi >= 0.26f && r.mpipi < 0.70f)) continue;
                if (!(r.m4 >= 0.30f && r.m4 < 0.85f)) continue;
                ck.push_back((int8_t)k); cg.push_back(r.grp); cm4.push_back(r.m4); cpp.push_back(r.mpipi); cee.push_back(r.mee);
                cb.push_back((uint8_t)((s0 ? 1 : 0) | 2));
            }
        W("d_k.i8", ck); W("d_grp.i32", cg); W("d_m4.f32", cm4); W("d_mpipi.f32", cpp); W("d_mee.f32", cee); W("d_bits.u8", cb);
        fprintf(lay, "data %zu\n", ck.size());
        printf("data rows %zu (%.0f s)\n", ck.size(), sw.RealTime()); sw.Continue();
    }
    {
        simphys::Init(); simphys::Install(); simphys::LoadAngles();
        simeng::SimData* D = new simeng::SimData();
        if (!simeng::Load(*D, simeng::kCompact, simeng::kMaskUS, -1e9, 1e9, 1)) return;
        const int trees[3] = {simeng::kStdUS, simeng::kLepUS, simeng::kEtaUS};
        std::vector<int8_t> ct, cgg; std::vector<int32_t> cgrp, proc[2], par[2]; std::vector<uint8_t> cb;
        std::vector<float> m4, m4t, pp, ppt, ee, eet, f108, P[4], PT[4], PR[4], TH[4], THT[4], PH[4], VR[2], VZ[2];
        std::vector<double> w;
        for (int t = 0; t < 3; ++t) {
            for (const simeng::Row& r : D->t[trees[t]]) {
                if (r.ownA != 1 && r.ownB != 1) continue;
                int cn[6]; const int ncn = CompsOf(r, false, cn);
                if (ncn == 0) continue;
                bool t12 = false; const int g = Grp(cn, ncn, t12);
                if (g < 0 || !t12) continue;
                const bool s0 = PassSel(kSelS0, r, r.cor_ok == 1, 0.f, r.mm4);
                const bool pre = PassSel(kSelPre, r, r.cor_ok == 1, 0.f, r.mm4);
                if (!s0 && !pre) continue;
                if (!(r.m4 >= 0.29f && r.m4 < 0.86f)) continue;
                ct.push_back((int8_t)trees[t]); cgrp.push_back(r.grp); cgg.push_back((int8_t)g);
                cb.push_back((uint8_t)((s0 ? 1 : 0) | (pre ? 2 : 0) | (r.ownA == 1 ? 8 : 0) | (r.ownB == 1 ? 16 : 0)));
                m4.push_back(r.m4); m4t.push_back(r.m4_true); pp.push_back(r.mpipi); ppt.push_back(r.mpipi_true);
                ee.push_back(r.mee); eet.push_back(r.mee_true); w.push_back(r.w * kSimToRaw);
                f108.push_back((float)simphys::RowFactor(r, kPhysNom));
                const float a[4] = {r.ep_p, r.em_p, r.pip_p, r.pim_p}, at[4] = {r.ep_sim_p, r.em_sim_p, r.pip_sim_p, r.pim_sim_p};
                const float ar[4] = {r.ep_p_rec, r.em_p_rec, r.pip_p_rec, r.pim_p_rec};
                const float th[4] = {r.ep_theta, r.em_theta, r.pip_theta, r.pim_theta};
                const float tht[4] = {r.ep_sim_theta, r.em_sim_theta, r.pip_sim_theta, r.pim_sim_theta};
                const float ph[4] = {r.ep_phi, r.em_phi, r.pip_phi, r.pim_phi};
                for (int i = 0; i < 4; ++i) { P[i].push_back(a[i]); PT[i].push_back(at[i]); PR[i].push_back(ar[i]); TH[i].push_back(th[i]); THT[i].push_back(tht[i]); PH[i].push_back(ph[i]); }
                VR[0].push_back(r.ep_sim_vr); VR[1].push_back(r.em_sim_vr); VZ[0].push_back(r.ep_sim_vz); VZ[1].push_back(r.em_sim_vz);
                proc[0].push_back(r.ep_sim_processid); proc[1].push_back(r.em_sim_processid);
                par[0].push_back(r.ep_sim_parentid); par[1].push_back(r.em_sim_parentid);
            }
        }
        delete D;
        W("s_tree.i8", ct); W("s_grp.i32", cgrp); W("s_g.i8", cgg); W("s_bits.u8", cb); W("s_m4.f32", m4); W("s_m4t.f32", m4t);
        W("s_pp.f32", pp); W("s_ppt.f32", ppt); W("s_ee.f32", ee); W("s_eet.f32", eet); W("s_w.f64", w); W("s_f108.f32", f108);
        const char* n[4] = {"ep", "em", "pip", "pim"};
        for (int i = 0; i < 4; ++i) {
            W(Form("s_%s_p.f32", n[i]), P[i]); W(Form("s_%s_pt.f32", n[i]), PT[i]); W(Form("s_%s_pr.f32", n[i]), PR[i]);
            W(Form("s_%s_th.f32", n[i]), TH[i]); W(Form("s_%s_tht.f32", n[i]), THT[i]); W(Form("s_%s_ph.f32", n[i]), PH[i]);
        }
        W("s_ep_vr.f32", VR[0]); W("s_em_vr.f32", VR[1]); W("s_ep_vz.f32", VZ[0]); W("s_em_vz.f32", VZ[1]);
        W("s_ep_proc.i32", proc[0]); W("s_em_proc.i32", proc[1]); W("s_ep_par.i32", par[0]); W("s_em_par.i32", par[1]);
        fprintf(lay, "sim %zu\n", ct.size());
        printf("sim rows %zu (%.0f s)\n", ct.size(), sw.RealTime());
    }
    fclose(lay);
}
