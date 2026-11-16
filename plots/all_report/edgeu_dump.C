// ============================================================================
// edgeu_dump.C  --  package edgeu (the eta bump in the distance to the kinematic edge, u = M4 - E(M(pipi), M(ee))): row dump
// ============================================================================
// Written 2026-11-15 (prefix edgeu). Reads the data compact rows (scan_v2, ScanEngine.h) and the simulated e+e- rows
// (std + lep + eta, SimEngine.h) ONCE and writes the few columns that the u analysis needs as flat column files
// (little endian, no header; the row counts and the column types are in <out>/dump_layout.txt). The analysis itself
// (edgeu_core.py, edgeu_ana.py, edgeu_toys.py) is done in numpy on these columns, so that every u definition (edge
// function, mass scales, M(ee) dependence, resolution transfer) can be recomputed row by row without a new pass.
//
// data   cor_ok == 1 rows passing "pre" (S0 is a subset), 0.26 <= M(pipi) < 0.46, 0.30 <= M4 < 0.85 (as bedge_fill.C:
//        the cor_ok fallback rows are not used). Columns: k (0 us, 1 pp, 2 mm), grp, m4, mpipi, mee, bits (1 S0, 2 pre,
//        4 isBest == 1), file_id (day), fold.
// sim    owned rows (mode A or B) of template T12 (classes C2-C5 and the template-2 groups G_eta, G_D, G_C; never T3, T4,
//        X17), passing S0 or pre (PassSel with cor_ok == 1, dz = 0, as bedge_fill.C), 0.26 <= 0.9952 M(pipi) < 0.46,
//        0.29 <= M4 < 0.86. Columns: tree (simeng TreeIdx 0 std, 3 lep, 6 eta), grp, g (bedge group: 0 C2, 1 C3, 2 C4,
//        3 C5, 4 G_eta, 5 G_D omega/rho0 pions, 6 G_D other, 7 G_C omega/rho0, 8 G_C other), cls, lcls, pcls, bits (1 S0,
//        2 pre, 4 isBest == 1, 8 ownA, 16 ownB), m4, m4_true, mpipi, mpipi_true, mee, mee_true, w (= sim w x 27.0718:
//        expected RAW data rows at L = 6079.2 nb^-1 before the physics factor), and the physics factors
//        simphys::RowFactor for the codes 108 (V3, nominal), 116 (bump 1 + aY), 100 (bump phase space), 107 (direct V2),
//        76 (conversion F_V off); code 0 (weights off) is 1 by definition.
// checks the frozen edge function alldefs::Edge at the report's check points; the WP0 replica generator
//        (bedge::GroupPois, seed 0x5157A11) for the first groups of each tree (edgeu_core.py reproduces it).
//
// Run (from /home/przygoda/HADES/FAT/FAT; about 10 GB, a few minutes):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_edgeu", kTRUE)' \
//        'plots/all_report/edgeu_dump.C+O'
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

namespace edgeud {
using namespace bedge;
const char* const kOutU = "/hdd3/przygoda/HADES/FAT_out/all_report/edgeu";
FILE* gLog = nullptr;
#define PL(...) do { printf(__VA_ARGS__); if (gLog) fprintf(gLog, __VA_ARGS__); fflush(stdout); } while (0)

template <class T> bool WCol(const std::string& name, const std::vector<T>& v) {
    const std::string f = std::string(kOutU) + "/" + name;
    FILE* fp = fopen(f.c_str(), "wb");
    if (!fp) { PL("edgeu_dump: cannot write %s\n", f.c_str()); return false; }
    fwrite(v.data(), sizeof(T), v.size(), fp);
    fclose(fp);
    return true;
}
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
}  // namespace edgeud

void edgeu_dump() {
    using namespace edgeud;
    TStopwatch total;
    gSystem->mkdir(kOutU, kTRUE);
    gLog = fopen(Form("%s/edgeu_dump.log", kOutU), "w");
    PL("# edgeu_dump.C, 2026-11-15: row dump of the package edgeu (bump in u = M4 - E(M(pipi), M(ee)))\n");
    PL("# data %s\n# sim  %s + %s\n", scan::kCompact, simeng::kCompact, simeng::kCompactEta);
    // ---- frozen edge function at the check points of the report (Sect. 2): 0.51859, 0.50077, 0.50958, 0.47390
    PL("A. alldefs::Edge: E(0.30,0) %.5f  E(0.39,0) %.5f  E(0.39,0.060) %.5f  E(0.4116,0) %.5f  E(0.30,0.100) %.5f  E(0.4120,0) %.5f\n",
       Edge(0.30, 0.), Edge(0.39, 0.), Edge(0.39, 0.060), Edge(0.4116, 0.), Edge(0.30, 0.100), Edge(0.4120, 0.));
    PL("   scan::EdgeMeasured (bedge/wp6 edge, scale %.6f): E(0.30) %.5f  E(0.39) %.5f\n", scan::kScale, scan::EdgeMeasured(0.30), scan::EdgeMeasured(0.39));
    // ---- WP0 replica generator: reference values for edgeu_core.py
    {
        double k[8];
        const int trees[3] = {simeng::kStdUS, simeng::kLepUS, simeng::kEtaUS};
        for (int t = 0; t < 3; ++t)
            for (int g = 0; g < 3; ++g) {
                GroupPois(kSeedSim, trees[t], 1000 * g + 7, 8, k);
                PL("   GroupPois(seed 0x5157A11, tree %d, grp %d): %g %g %g %g %g %g %g %g\n", trees[t], 1000 * g + 7, k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7]);
            }
    }

    // ---------------------------------------------------------------------------------------------------- data
    {
        std::vector<int8_t> ck, cfold; std::vector<int32_t> cgrp; std::vector<float> cm4, cpp, cee; std::vector<uint8_t> cbits; std::vector<int16_t> cday;
        scan::Data d;
        if (!scan::Load(d)) return;
        long nNaN = 0, nk[3] = {0, 0, 0};
        for (int k = 0; k < 3; ++k) {
            for (const auto& r : d.s[k]) {
                if (r.cor_ok != 1) continue;
                const bool pre = PassSel(kSelPre, r, true, r.dz_pi, r.mm4);
                if (!pre) continue;
                const bool s0 = PassSel(kSelS0, r, true, r.dz_pi, r.mm4);
                if (!(r.mpipi >= 0.26f && r.mpipi < 0.46f)) continue;
                if (!(r.m4 >= 0.30f && r.m4 < 0.85f)) continue;
                if (!(r.mee == r.mee)) ++nNaN;
                ck.push_back((int8_t)k); cgrp.push_back(r.grp); cm4.push_back(r.m4); cpp.push_back(r.mpipi); cee.push_back(r.mee);
                cbits.push_back((uint8_t)((s0 ? 1 : 0) | 2 | (r.isBest == 1 ? 4 : 0)));
                cday.push_back((int16_t)r.file_id); cfold.push_back((int8_t)r.fold);
                ++nk[k];
            }
        }
        PL("B. data rows kept (pre, cor_ok == 1, 0.26 <= M(pipi) < 0.46, 0.30 <= M4 < 0.85): us %ld pp %ld mm %ld; without M(ee): %ld\n", nk[0], nk[1], nk[2], nNaN);
        WCol("data_k.i8", ck); WCol("data_grp.i32", cgrp); WCol("data_m4.f32", cm4); WCol("data_mpipi.f32", cpp); WCol("data_mee.f32", cee);
        WCol("data_bits.u8", cbits); WCol("data_day.i16", cday); WCol("data_fold.i8", cfold);
        FILE* fl = fopen(Form("%s/dump_layout.txt", kOutU), "w");
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
        std::vector<float> cm4, cm4t, cpp, cppt, cee, ceet, f108, f116, f100, f107, f76; std::vector<double> cw;
        const int codes[5] = {kPhysNom, kPhysBump1aY, kPhysBumpPS, kPhysV2, kPhysConvOff};
        std::vector<float>* fc[5] = {&f108, &f116, &f100, &f107, &f76};
        long ng[kNG] = {0};
        double sw108[kNG] = {0};
        TStopwatch sw;
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
                if (!(mpp >= 0.26 && mpp < 0.46)) continue;
                if (!(r.m4 >= 0.29f && r.m4 < 0.86f)) continue;
                ct.push_back((int8_t)trees[t]); cgrp.push_back(r.grp); cg.push_back((int8_t)g);
                ccls.push_back((int8_t)r.cls); clcls.push_back((int8_t)r.lcls); cpcls.push_back((int8_t)r.pcls);
                cbits.push_back((uint8_t)((s0 ? 1 : 0) | (pre ? 2 : 0) | (r.isBest == 1 ? 4 : 0) | (r.ownA == 1 ? 8 : 0) | (r.ownB == 1 ? 16 : 0)));
                cm4.push_back(r.m4); cm4t.push_back(r.m4_true); cpp.push_back(r.mpipi); cppt.push_back(r.mpipi_true);
                cee.push_back(r.mee); ceet.push_back(r.mee_true);
                cw.push_back(r.w * kSimToRaw);
                for (int q = 0; q < 5; ++q) fc[q]->push_back((float)simphys::RowFactor(r, codes[q]));
                ++ng[g];
                if (s0 && r.ownA == 1) sw108[g] += r.w * kSimToRaw * f108.back();
            }
            PL("C. sim tree %s: %zu rows read, %zu kept so far, %.0f s\n", simeng::TreeName(trees[t]), v.size(), ct.size(), sw.RealTime());
            sw.Continue();
        }
        delete D;
        PL("   kept rows per group:");
        for (int g = 0; g < kNG; ++g) PL(" %s %ld", kGTag[g], ng[g]);
        PL("\n   S0, mode A, V3: expected raw rows per group in the dump range:");
        for (int g = 0; g < kNG; ++g) PL(" %s %.1f", kGTag[g], sw108[g]);
        PL("\n");
        WCol("sim_tree.i8", ct); WCol("sim_grp.i32", cgrp); WCol("sim_g.i8", cg); WCol("sim_cls.i8", ccls); WCol("sim_lcls.i8", clcls);
        WCol("sim_pcls.i8", cpcls); WCol("sim_bits.u8", cbits);
        WCol("sim_m4.f32", cm4); WCol("sim_m4true.f32", cm4t); WCol("sim_mpipi.f32", cpp); WCol("sim_mpipitrue.f32", cppt);
        WCol("sim_mee.f32", cee); WCol("sim_meetrue.f32", ceet); WCol("sim_w.f64", cw);
        WCol("sim_f108.f32", f108); WCol("sim_f116.f32", f116); WCol("sim_f100.f32", f100); WCol("sim_f107.f32", f107); WCol("sim_f76.f32", f76);
        FILE* fl = fopen(Form("%s/dump_layout.txt", kOutU), "a");
        if (fl) {
            fprintf(fl, "sim %zu\n", ct.size());
            fprintf(fl, "# data columns: data_k.i8 data_grp.i32 data_m4.f32 data_mpipi.f32 data_mee.f32 data_bits.u8 data_day.i16 data_fold.i8\n");
            fprintf(fl, "# sim columns: sim_tree.i8 sim_grp.i32 sim_g.i8 sim_cls.i8 sim_lcls.i8 sim_pcls.i8 sim_bits.u8 sim_m4.f32 sim_m4true.f32 sim_mpipi.f32\n");
            fprintf(fl, "#   sim_mpipitrue.f32 sim_mee.f32 sim_meetrue.f32 sim_w.f64 sim_f108.f32 sim_f116.f32 sim_f100.f32 sim_f107.f32 sim_f76.f32\n");
            fprintf(fl, "# groups: 0 C2, 1 C3, 2 C4, 3 C5, 4 Geta, 5 GDom, 6 GDnom, 7 GCom, 8 GCnom; bits data 1 S0 2 pre 4 best; sim 1 S0 2 pre 4 best 8 ownA 16 ownB\n");
            fclose(fl);
        }
    }
    PL("done in %.0f s\n", total.RealTime());
    if (gLog) fclose(gLog);
}
