// ============================================================================
// x17d_toycls.C  --  WP9 X-d (prefix x17d): toy-based CLs on the UNBLINDED DATA (x17_frozen.txt Sect. 6 and 10 step 3)
// ============================================================================
// Frozen procedure: x17c::ToyCLs (1000 toys per point and hypothesis; pseudo-data from the conditional fits of the data at
// r and at 0, over-dispersed by qScale, constraint centres drawn; q~_r / qScale of every toy); CLs scanned on a grid of r
// around the observed asymptotic UL; toy UL = x17c::ToyUL (weighted quadratic fit of ln CLs, CLs = 0.10). If toy and
// asymptotic UL differ by more than 10 %, the toy UL is primary. If CL_b = 0, it is reported and the asymptotic CLs is used.
// One call = one mass and one normalisation. Output: /hdd3/.../x17d/toycls_m<100 m>_n<norm>.txt
// Run: root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_x17d", kTRUE)' \
//        'plots/all_report/x17d_toycls.C+O(17.0, 1, 1000)'
// ============================================================================
#include <cmath>
#include <cstdio>
#include <vector>

#include "TFile.h"
#include "TStopwatch.h"
#include "TSystem.h"

#include "x17d_lib.h"

void x17d_toycls(double mX = 17.0, int norm = 1, int ntoys = 1000, unsigned long long seed = 20261115ULL) {
    using namespace x17d;
    TStopwatch sw;
    Ctx c;
    if (!Init(c)) { printf("x17d_toycls: init failed\n"); return; }
    x17c::Inputs I = c.I;
    if (!SetLine(I, mX)) return;
    const x17c::Config C = MakeConfig(c.z, norm);
    const x17c::Model M(I, C);
    const x17c::Glob g = x17c::GlobOf(C);
    const std::vector<double>& n = c.D.n;
    // observed asymptotic limit and its Asimov set (exactly x17c::LimitAsym)
    const x17c::Limit La = x17c::LimitAsym(M, n, g, 0.1, false);
    const x17c::Profile& P = La.prof;
    x17c::Profile PA;
    PA.zero = M.Fit(La.asimov, 0., false, g, &P.zero.p); PA.hat = PA.zero; PA.rhat = 0; PA.Fref = PA.zero.F;
    const double ulr = La.ul;
    FILE* fo = fopen(Form("%s/toycls_m%04d_n%d.txt", kOutD, (int)std::lround(100 * mX), norm), "w");
    fprintf(fo, "# x17d_toycls: m %.2f norm %d: observed asymptotic UL r %.6g (BR %.4e), r_hat %.5g +- %.5g; %d toys per point and hypothesis; qScale %.4f\n", mX, norm,
            ulr, ulr * kBR, P.rhat, P.rhatErr, ntoys, C.qScale);
    fprintf(fo, "# r  BR  CLs_toy  err  CLsb  CLb  nt  CLs_asym  q_obs\n");
    printf("x17d_toycls: m %.2f norm %d: asymptotic UL %.4e (r %.5f)\n", mX, norm, ulr * kBR, ulr);
    const double fac[9] = {0.55, 0.70, 0.85, 0.95, 1.05, 1.15, 1.30, 1.50, 1.80};
    std::vector<x17c::ToyCLsPoint> pts;
    bool clb0 = false;
    for (int k = 0; k < 9; ++k) {
        const double r = ulr * fac[k];
        const double qo = x17c::QTilde(M, n, g, P, r), qa = x17c::QTilde(M, La.asimov, g, PA, r);
        const double clsA = x17c::CLsAsym(qo, qa);
        const unsigned long long sd = seed + 1000003ULL * (unsigned long long)std::lround(100 * mX) + 7919ULL * norm + 104729ULL * k;
        x17c::ToyCLsPoint T = x17c::ToyCLs(M, n, g, P, r, ntoys, sd);
        if (T.clb <= 0) clb0 = true;
        pts.push_back(T);
        fprintf(fo, "%.6g %.5e %.5f %.5f %.5f %.5f %d %.5f %.4f\n", r, r * kBR, T.cls, T.clsErr, T.clsb, T.clb, T.nt, clsA, T.qobs);
        fflush(fo);
        printf("x17d_toycls: m %.2f norm %d r %.5f BR %.4e CLs_toy %.4f +- %.4f (CLsb %.4f CLb %.4f) CLs_asym %.4f q_obs %.3f (%.0f s)\n", mX, norm, r, r * kBR, T.cls,
               T.clsErr, T.clsb, T.clb, clsA, T.qobs, sw.RealTime());
        sw.Continue();
    }
    const double ulT = x17c::ToyUL(pts);
    fprintf(fo, "RESULT m %.2f norm %d UL_asym %.5e UL_toy %.5e ratio %.4f clb0 %d\n", mX, norm, ulr * kBR, ulT * kBR, ulT / ulr, clb0 ? 1 : 0);
    fclose(fo);
    printf("x17d_toycls: RESULT m %.2f norm %d UL asym %.4e toy %.4e ratio %.4f%s, %.0f s\n", mX, norm, ulr * kBR, ulT * kBR, ulT / ulr, clb0 ? " (CL_b = 0 at some point)" : "",
           sw.RealTime());
    c.f->Close();
}
