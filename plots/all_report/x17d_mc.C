// ============================================================================
// x17d_mc.C  --  WP9 X-d (prefix x17d): template statistics ('MC', quoted separately, D13; x17_frozen.txt Sect. 4, 10 step 5)
// ============================================================================
// The observed limit and r_hat are recomputed with the 200 CB replicas (CBrep_X2D: Poisson event-group replicas of the
// like-sign rows) and with the 100 simulation replicas of the conversion template C4 (C4_X2D_rep); everything else is the
// frozen nominal. The RMS of the UL over the replicas is the 'MC' uncertainty of the limit.
// One call = one mass, all three normalisations. Output: /hdd3/.../x17d/mc_m<100 m>.txt (one line per replica)
// Run: root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_x17d", kTRUE)' \
//        'plots/all_report/x17d_mc.C+O(17.0)'
// ============================================================================
#include <cmath>
#include <cstdio>
#include <vector>

#include "TFile.h"
#include "TStopwatch.h"

#include "x17d_lib.h"

void x17d_mc(double mX = 17.0, int nCB = 200, int nC4 = 100) {
    using namespace x17d;
    TStopwatch sw;
    Ctx c;
    if (!Init(c)) { printf("x17d_mc: init failed\n"); return; }
    FILE* fo = fopen(Form("%s/mc_m%04d.txt", kOutD, (int)std::lround(100 * mX)), "w");
    fprintf(fo, "# x17d_mc: m %.2f; columns: kind(0 nominal, 1 CB replica, 2 C4 replica) index UL1 UL2 UL3 rhat1 rhat2 rhat3 ok1 ok2 ok3 (BR units)\n", mX);
    auto one = [&](int kind, int idx) {
        BuildOpt b;
        if (kind == 1) b.cbRep = idx;
        if (kind == 2) b.simRep = idx;
        x17c::Inputs Iv;
        if (!Build(c.f, Iv, c.w, b, c.z.a0)) return;
        if (!SetLine(Iv, mX)) return;
        double ul[3], rh[3]; int ok[3];
        for (int a = 0; a < 3; ++a) {
            Lim l = LimitBR(Iv, MakeConfig(c.z, a + 1), c.D.n, false);
            ul[a] = l.ul; rh[a] = l.rhat; ok[a] = l.ok;
        }
        fprintf(fo, "%d %d %.5e %.5e %.5e %+.5e %+.5e %+.5e %d %d %d\n", kind, idx, ul[0], ul[1], ul[2], rh[0], rh[1], rh[2], ok[0], ok[1], ok[2]);
        fflush(fo);
    };
    one(0, -1);
    for (int k = 0; k < nCB; ++k) one(1, k);
    for (int k = 0; k < nC4; ++k) one(2, k);
    fclose(fo);
    printf("x17d_mc: m %.2f done, %.0f s\n", mX, sw.RealTime());
    c.f->Close();
}
