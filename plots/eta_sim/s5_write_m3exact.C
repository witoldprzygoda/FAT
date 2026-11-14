// s5_write_m3exact.C -- exact (unbinned) half-A and half-B sums for the linear side-band method M3 of the blind test
// (review finding L5: s5_blind_eval.C used the binned dev.log sums of the analyst for its "second realisation").
// SIMULATION ONLY. std + lep + eta e+e- trees, rows owned in mode A, w (PDG) x physics weights V3 (kPhysMain), base pre.
// Halves: A = even grp, B = odd grp (as s5_blind_prep.C). Sums of w of one half = expected rows at L_toy = L_sim / 2.
// M3: N_win [0.5333,0.5587), N_L [0.5255,0.5333), N_U [0.5587,0.5665) (T1 + T2 rows); B_lin = (N_L + N_U)/2 x 25.4/7.8;
//     lambda = leakage of the half-A peak45 (classes 4+5) into the side bands, (T_L + T_U)/2 x 25.4/7.8 / T_win;
//     S_raw = (N_win - B_lin) / (1 - lambda); relative bias = S_raw / T_win(same half) - 1.
// Selections (frozen eta cuts): PEAK = dz_pi < 40 mm && MM4 < 2.85; DIRECT = dz_pi < 40 mm; BUMP = pre only.
// Run from /home/przygoda/HADES/FAT/FAT:
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/s5_write/aclic_m3", kTRUE)' 'plots/eta_sim/s5_write_m3exact.C+'
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>
#include "SimEngine.h"
#include "sim_eta_physics.h"

using namespace simeng;
namespace s5m3 {
struct R { float m4, mm4, dz; unsigned char half, pk45, t12; double w; };
std::vector<R> gR;
}  // namespace s5m3

void s5_write_m3exact() {
    using namespace s5m3;
    simphys::Init();
    simphys::Install();
    simphys::LoadAngles();
    {
        SimData d;
        if (!Load(d, kCompact, kMaskUS)) return;
        const int trees[3] = {kStdUS, kLepUS, kEtaUS};
        for (int it = 0; it < 3; ++it)
            for (const Row& r : d.t[trees[it]]) {
                if (r.ownA != 1 || !PassPre(r)) continue;
                if (!(r.m4 >= 0.52f && r.m4 < 0.57f)) continue;
                if (r.tmpl != 1 && r.tmpl != 2) continue;           // T1 + T2 only (the simulated sum)
                if (r.tmpl == 1 && !(r.cls >= 2 && r.cls <= 5)) continue;   // X17 (class 6) never summed
                R x;
                x.m4 = r.m4; x.mm4 = r.mm4;
                x.dz = std::max(std::fabs(r.pip_z - r.vz), std::fabs(r.pim_z - r.vz));
                x.half = r.grp & 1;
                x.pk45 = (r.tmpl == 1 && (r.cls == 4 || r.cls == 5));
                x.t12 = 1;
                x.w = r.w * simphys::RowFactor(r, simphys::kPhysMain);
                gR.push_back(x);
            }
        printf("rows (pre, owned A, T1 classes 2-5 + T2, M4 0.52-0.57): %zu; angle hits %lld misses %lld\n", gR.size(), simphys::AngHit(),
               simphys::AngMiss());
    }
    const char* sn[3] = {"PEAK", "DIRECT", "BUMP"};
    printf("\n# exact M3 sums, rows at L_toy = 110.525 nb^-1 (sum w of one half); simulation only\n");
    printf("%-7s %-5s %9s %9s %9s %9s %9s %9s %8s %9s %9s\n", "sel", "half", "N_win", "N_L", "N_U", "T_win", "T_L", "T_U", "lambda", "S_raw", "rel%");
    for (int s = 0; s < 3; ++s) {
        auto sel = [s](const R& r) { return s == 0 ? (r.dz < 40.f && r.mm4 < 2.85f) : s == 1 ? (r.dz < 40.f) : true; };
        double lamA = 0;
        for (int h = 0; h < 2; ++h) {
            double N[3] = {0, 0, 0}, T[3] = {0, 0, 0};
            for (const R& r : gR) {
                if (r.half != h || !sel(r)) continue;
                const int k = (r.m4 >= 0.5333f && r.m4 < 0.5587f) ? 0 : (r.m4 >= 0.5255f && r.m4 < 0.5333f) ? 1 : (r.m4 >= 0.5587f && r.m4 < 0.5665f) ? 2 : -1;
                if (k < 0) continue;
                N[k] += r.w;
                if (r.pk45) T[k] += r.w;
            }
            const double lam = (T[1] + T[2]) / 2 * 25.4 / 7.8 / T[0];
            if (h == 0) lamA = lam;
            // the analyst's lambda comes from half A for both halves
            const double Sraw = (N[0] - (N[1] + N[2]) / 2 * 25.4 / 7.8) / (1 - lamA);
            printf("%-7s %-5s %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %8.4f %9.2f %+9.1f\n", sn[s], h ? "B" : "A", N[0], N[1], N[2], T[0], T[1], T[2], lam,
                   Sraw, 100 * (Sraw / T[0] - 1));
        }
    }
}
