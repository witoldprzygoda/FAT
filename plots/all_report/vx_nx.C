// ============================================================================
// vx_nx.C  --  independent verification (prefix vx, 2026-11-14): a direct 90 % CLs limit on the NUMBER of X17 rows N_X
// (raw rows, X2D M4 range, all M(ee)) with the frozen model, data and test statistic. POI = N_X itself: normalisation-(2)
// form X17 rows = r xA N5ref R_X with the X17 yield nuisance switched off (xSig = 0), because a limit on rows needs no
// efficiency, BR_dir or eta-yield factor. Line-shape nuisances (width, M(ee) scale) stay as frozen. READ-ONLY.
// Run: root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_vx", kTRUE)' 'plots/all_report/vx_nx.C+O'
// ============================================================================
#include <cmath>
#include <cstdio>

#include "x17d_lib.h"

void vx_nx() {
    using namespace x17d;
    Ctx c;
    if (!Init(c)) { printf("vx_nx: init failed\n"); return; }
    FILE* fo = fopen("/hdd3/przygoda/HADES/FAT_out/all_report/vx/vx_nx.txt", "w");
    for (double m : {16.7, 16.85, 17.0}) {
        x17c::Inputs I = c.I;
        if (!SetLine(I, m)) continue;
        x17c::Config C = MakeConfig(c.z, 2);
        C.xSig = 0.;                 // no yield nuisance: the POI is the number of rows
        const double perR = C.xA * I.N5ref * I.RX;   // X17 rows (all M(ee)) per unit r
        Lim l = LimitBR(I, C, c.D.n, true);
        const double ulRows = l.ul / kBR * perR, rh = l.rhat / kBR * perR, re = l.rerr / kBR * perR, med = l.exp[2] / kBR * perR;
        // the same in BR units for the three normalisations (information): rows / (N5_hat(r=0) R_X) x BR_dir, rows / (a s_X2D)
        const double sX2D = I.N5ref * I.RX / kBR;
        char s[1024];
        snprintf(s, sizeof s,
                 "vx_nx: m %.2f: 90%% CLs UL on N_X = %.1f raw rows (%.1f after x 64/63); N_X_hat %+.1f +- %.1f; Asimov median (observed r = 0 fit) %.1f rows;"
                 " N_X^UL / (a_eta s_X2D) = %.3e, N_X^UL / (s_X2D N5_hat/N5ref) = %.3e (N5_hat of the free fit %.1f)\n",
                 m, ulRows, ulRows * 64. / 63., rh, re, med, ulRows / (c.z.a0 * sX2D), ulRows / (sX2D * l.N5 / I.N5ref), l.N5);
        printf("%s", s); fputs(s, fo);
    }
    fclose(fo);
    c.f->Close();
}
