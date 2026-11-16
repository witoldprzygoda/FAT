// ============================================================================
// x17d_fc.C  --  WP9 X-d (prefix x17d): Feldman-Cousins cross-check at m_X = 17 MeV/c2 (x17_frozen.txt Sect. 6)
// ============================================================================
// Unified interval in r from toys of the profile-likelihood ratio: t_r = (F(r) - F(max(r_hat, 0))) / qScale (two-sided,
// physical boundary r >= 0). For each r of the list: t_obs(r) on the data; toys from the conditional fit of the DATA at r
// (profile construction; x17c::ToyData, over-dispersed by qScale, constraint centres drawn); t_c(r) = 90 % quantile of
// t_r of the toys. r belongs to the 90 % CL interval if t_obs(r) <= t_c(r). Frozen model and configuration.
// One call = one value of r. Output line appended to /hdd3/.../x17d/fc_m<100 m>_n<norm>_r<1e5 r>.txt
// Run: root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_x17d", kTRUE)' \
//        'plots/all_report/x17d_fc.C+O(0.10, 2000)'
// ============================================================================
#include <cmath>
#include <cstdio>
#include <vector>

#include "TFile.h"
#include "TStopwatch.h"

#include "x17d_lib.h"

void x17d_fc(double r = 0.10, int ntoys = 2000, double mX = 17.0, int norm = 1, unsigned long long seed = 20261116ULL) {
    using namespace x17d;
    TStopwatch sw;
    Ctx c;
    if (!Init(c)) { printf("x17d_fc: init failed\n"); return; }
    x17c::Inputs I = c.I;
    if (!SetLine(I, mX)) return;
    const x17c::Config C = MakeConfig(c.z, norm);
    const x17c::Model M(I, C);
    const x17c::Glob g = x17c::GlobOf(C);
    const std::vector<double>& n = c.D.n;
    const x17c::Profile P = x17c::ProfileFits(M, n, g);
    const unsigned long long sd = seed + 1000003ULL * (unsigned long long)std::lround(1e5 * r) + 7919ULL * norm;
    FCPoint T = FCToys(M, n, g, P, r, ntoys, sd);
    FILE* fo = fopen(Form("%s/fc_m%04d_n%d_r%06d.txt", kOutD, (int)std::lround(100 * mX), norm, (int)std::lround(1e5 * r)), "w");
    fprintf(fo, "# x17d_fc: m %.2f norm %d; r_hat %.6g +- %.6g; columns: r BR t_obs t68 t90 t95 p(t >= t_obs) ntoys nok\n", mX, norm, P.rhat, P.rhatErr);
    fprintf(fo, "%.6g %.5e %.5f %.5f %.5f %.5f %.5f %d %d\n", T.r, T.r * kBR, T.tobs, T.t68, T.t90, T.t95, T.pval, T.nt, T.nok);
    fclose(fo);
    printf("x17d_fc: m %.2f norm %d r %.5f BR %.4e t_obs %.4f t68 %.4f t90 %.4f t95 %.4f p %.4f (%d/%d toys) %.0f s\n", mX, norm, T.r, T.r * kBR, T.tobs, T.t68, T.t90,
           T.t95, T.pval, T.nok, T.nt, sw.RealTime());
    c.f->Close();
}
