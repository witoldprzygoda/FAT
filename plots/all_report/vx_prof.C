// ============================================================================
// vx_prof.C  --  independent verification of the X17 limit (prefix vx, 2026-11-14). READ-ONLY use of the frozen model.
// ============================================================================
// With the frozen inputs (x17d::Init = frozen build + full data) at a mass m_X and normalisation 1:
//   1. profile scan F(r) = -2 ln L_profile(r) (all nuisances refitted) on a grid of r;
//   2. Bayesian 90 % upper limit with a flat prior in r >= 0 on the profile likelihood exp(-(F(r) - F_min) / (2 qScale));
//   3. the asymptotic CLs recomputed by hand from q~_obs(r) and q_A(r) of the same scan (Cowan et al. eqs. 65-66), to check
//      the frozen LimitAsym arithmetic;
//   4. the conditional fit at the observed UL: N5, theta_X and the X17 rows actually in the model at the limit;
//   5. a background-only Asimov check of sigma(r) at small and large r (nonlinearity of normalisation 1).
// Output: stdout and /hdd3/przygoda/HADES/FAT_out/all_report/vx/vx_prof_m<100 m>.txt
// Run: root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_vx", kTRUE)' \
//        'plots/all_report/vx_prof.C+O(16.85)'
// ============================================================================
#include <cmath>
#include <cstdio>
#include <vector>

#include "TFile.h"
#include "TMath.h"

#include "x17d_lib.h"

void vx_prof(double mX = 16.85, int norm = 1) {
    using namespace x17d;
    Ctx c;
    if (!Init(c)) { printf("vx_prof: init failed\n"); return; }
    x17c::Inputs I = c.I;
    if (!SetLine(I, mX)) return;
    const x17c::Config C = MakeConfig(c.z, norm);
    const x17c::Model M(I, C);
    const x17c::Glob g = x17c::GlobOf(C);
    const std::vector<double>& n = c.D.n;
    const double qs = C.qScale;
    FILE* fo = fopen(Form("/hdd3/przygoda/HADES/FAT_out/all_report/vx/vx_prof_m%04d_n%d.txt", (int)std::lround(100 * mX), norm), "w");
    auto P = [&](const char* s) { printf("%s", s); fputs(s, fo); };
    // frozen limit (reference)
    const x17c::Limit L = x17c::LimitAsym(M, n, g, 0.1, true);
    const x17c::Profile& Pr = L.prof;
    P(Form("vx_prof: m %.2f norm %d: frozen LimitAsym UL r %.6f (BR %.4e), r_hat %.5f +- %.5f, Asimov exp median %.4e\n", mX, norm, L.ul, L.ul * kBR, Pr.rhat,
           Pr.rhatErr, L.exp[2] * kBR));
    // Asimov profile
    x17c::Profile PA;
    PA.zero = M.Fit(L.asimov, 0., false, g, &Pr.zero.p); PA.hat = PA.zero; PA.rhat = 0; PA.Fref = PA.zero.F;
    // 1. profile scan
    std::vector<double> rr, FF, FA;
    const double rmax = 0.45, dr = 0.0025;
    std::vector<double> st = Pr.hat.p;
    double Fmin = 1e300;
    for (double r = -0.10; r <= rmax + 1e-9; r += dr) {
        x17c::FitResult f = M.Fit(n, r, false, g, &st);
        if (!f.ok) { x17c::FitResult f2 = M.Fit(n, r, false, g); if (f2.F < f.F) f = f2; }
        st = f.p;
        x17c::FitResult fa = M.Fit(L.asimov, r, false, g, &PA.zero.p);
        rr.push_back(r); FF.push_back(f.F); FA.push_back(fa.F);
        Fmin = std::min(Fmin, f.F);
    }
    // 2. Bayesian flat prior r >= 0
    double tot = 0; std::vector<double> cum(rr.size(), 0.);
    for (size_t k = 1; k < rr.size(); ++k) {
        if (rr[k] <= 0) { cum[k] = 0; continue; }
        const double a = rr[k - 1] < 0 ? 0. : std::exp(-0.5 * (FF[k - 1] - Fmin) / qs), b = std::exp(-0.5 * (FF[k] - Fmin) / qs);
        tot += 0.5 * (a + b) * dr; cum[k] = tot;
    }
    double ulB = NAN;
    for (size_t k = 1; k < rr.size(); ++k) if (cum[k - 1] < 0.9 * tot && cum[k] >= 0.9 * tot) ulB = rr[k - 1] + dr * (0.9 * tot - cum[k - 1]) / (cum[k] - cum[k - 1]);
    P(Form("vx_prof: Bayesian 90%% UL (flat prior r >= 0, profile likelihood / qScale): r %.5f -> BR %.4e (posterior tail beyond r = %.2f: %.2e of total)\n", ulB,
           ulB * kBR, rmax, std::exp(-0.5 * (FF.back() - Fmin) / qs)));
    // 3. CLs by hand from the scan
    const double Fref = Pr.Fref, FA0 = PA.zero.F;
    double ulH = NAN, prev = 1.;
    for (size_t k = 0; k < rr.size(); ++k) {
        const double r = rr[k];
        if (r <= 0) continue;
        const double q = (r < Pr.rhat) ? 0. : std::max(0., FF[k] - Fref) / qs;
        const double qa = std::max(1e-12, FA[k] - FA0) / qs;
        const double sq = std::sqrt(q), sa = std::sqrt(qa);
        double ps, pb;
        if (q <= qa) { ps = 1 - TMath::Freq(sq); pb = TMath::Freq(sa - sq); }
        else { ps = 1 - TMath::Freq((q + qa) / (2 * sa)); pb = 1 - TMath::Freq((q - qa) / (2 * sa)); }
        const double cls = ps / pb;
        if (prev > 0.10 && cls <= 0.10 && k > 0) ulH = rr[k - 1] + dr * (prev - 0.10) / (prev - cls);
        if (std::fabs(r - 0.10) < 1e-6 || std::fabs(r - 0.115) < 1e-6)
            P(Form("vx_prof:   r %.4f: q~ %.4f qA %.4f CLs+b %.4f CLb %.4f CLs %.4f\n", r, q, qa, ps, pb, cls));
        prev = cls;
    }
    P(Form("vx_prof: CLs recomputed by hand from the scan (linear interpolation in r): UL r %.5f -> BR %.4e (frozen LimitAsym %.4e)\n", ulH, ulH * kBR, L.ul * kBR));
    // one-sided profile-likelihood limit without CLs (q = 2.706)
    P(Form("vx_prof: one-sided PL limit (q~ = 2.706, no CLs): BR %.4e\n", L.ulPL * kBR));
    // 4. conditional fit at the UL
    {
        x17c::FitResult f = M.Fit(n, L.ul, false, g, &Pr.hat.p);
        const x17c::Layout& Ly = M.L;
        const double lxs = C.xSig > 0 ? std::log(1. + C.xSig) : 0.;
        const double xf = std::exp(lxs * f.p[Ly.iTX]);
        double sx = 0; for (int j = 0; j < I.G.ny; ++j) sx += M.Xeff(f.p, j);
        const double rows = L.ul * M.Nnorm(f.p) * I.RX * xf;
        P(Form("vx_prof: conditional fit at the UL: N5 %.1f (r = 0 fit %.1f), theta_X %+.3f (yield factor %.3f), thW %+.3f thS %+.3f, a_eta %.3f;"
               " X17 rows in the model %.1f (all M(ee)), %.1f inside 2-62 MeV/c2; table D1 quotes r_UL N5(r=0) R_X = %.1f\n",
               f.p[Ly.iN5], Pr.zero.p[Ly.iN5], f.p[Ly.iTX], xf, f.p[Ly.iTW], f.p[Ly.iTS], f.p[Ly.iA], rows, rows * sx, L.ul * Pr.zero.p[Ly.iN5] * I.RX));
    }
    // 5. sigma(r) on the Asimov set
    for (double r : {0.02, 0.05, 0.10, 0.15}) {
        size_t k = (size_t)std::lround((r + 0.10) / dr);
        const double qa = std::max(1e-12, FA[k] - FA0) / qs;
        P(Form("vx_prof: Asimov sigma_A(r = %.3f) = r / sqrt(qA) = %.5f (BR %.3e)\n", rr[k], rr[k] / std::sqrt(qa), rr[k] / std::sqrt(qa) * kBR));
    }
    fprintf(fo, "# scan: r BR F_obs-Fmin F_asimov-FA0\n");
    for (size_t k = 0; k < rr.size(); ++k) fprintf(fo, "%.4f %.4e %.4f %.4f\n", rr[k], rr[k] * kBR, FF[k] - Fmin, FA[k] - FA0);
    fclose(fo);
    c.f->Close();
}
