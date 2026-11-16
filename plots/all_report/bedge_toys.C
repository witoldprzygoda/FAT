// ============================================================================
// bedge_toys.C  --  package bedge: closure of the edge estimator on pseudo-data with a known truth
// ============================================================================
// Written 2026-11-14 (prefix bedge). The WP4 machinery (wp4_lib.h) applied to the edge spectra:
//   truth (expected raw rows per 2.5 MeV M4 bin, per M(pipi) slice of the edge estimator and per M(ee) class E1/E2):
//     mu = CB (the data CB of the slice: 2 sqrt(N++ N--) of the data like-sign rows, class share by the spread rule)
//        + c x S(m) x (1 + b x + q x^2) [x (1 + A exp(-(m - m0)^2 / 2 w^2))]    non-eta background
//        + a x [c2 C2 + kappa C3 + kappa C4 + C5 + G_eta]                        eta classes, raw simulation (full sample)
//     S = the smoothed simulated non-eta background G_D + G_C of the slice and class (wp4::SmoothFit, exp(pol3) on
//     0.40-0.70); x = (M4 - 0.55)/0.1 as in D11; WP4 scenario (a): c = 1.35, b = -0.15, a = 0.4, kappa = c2 = 1.
//   toy: N = Poisson(mu) per bin and class; N++ and N-- = Poisson of the data like-sign counts; the CB is rebuilt by
//   the rule (2 sqrt in 10 MeV cells, class spread rule), exactly as for the data.
//   estimator: the frozen edge fit (FitStep = scan::FitEdgeStep) on the toy; f = E_toy / E_sim with
//     AA: E_sim of the full simulated sample (the same events as the truth: the estimator alone);
//     RF: E_sim of bootstrap replica r = toy mod 100 (template noise included, as WP4 "RF").
//   The corrected quantities (C2 angle, kappa fixed) use the linear decomposition of E_sim (data covariance), as
//   bedge_ana.C. Pull = (q - truth)/sqrt(stat_fit^2 + MC^2), MC = the relative replica spread of the nominal.
// Scenarios: background with a quadratic term (q = -0.173 is the stage-A value that moves B4 by +150 %, WP6), x2,
// opposite sign; WP4 (e) omega/rho0-pion background x 0.3 and the rest x 1.6; a broad (15 MeV) and a narrow (8 MeV)
// local background structure under the bump; mass transfer (s + 0.0009, r = 0.89) of the truth; a = 0 and a = 1;
// the gamma* angle (C2 x 0.8855) and a conversion scale kappa = 0.48 in the truth.
//
// Run (from /home/przygoda/HADES/FAT/FAT; about 2 GB):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_bedge", kTRUE)' \
//        'plots/all_report/bedge_toys.C+O(1000)'
// ============================================================================
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "TRandom3.h"
#include "TStopwatch.h"
#include "TSystem.h"

#include "wp4_lib.h"
#include "bedge_common.h"

namespace bedget {
using namespace bedge;
FILE* gL = nullptr;
#define PL(...) do { printf(__VA_ARGS__); if (gL) fprintf(gL, __VA_ARGS__); fflush(stdout); } while (0)

std::vector<double> DN, SS, SV, SR;
const double kBwD = (scan::kMhi - scan::kMlo) / scan::kNbin;
const double kRLo = scan::kMlo + kRM4lo * kBwD;
const int kNSl = 4;
const double kSlLo[kNSl] = {0.30, 0.34, 0.37, 0.39}, kSlHi[kNSl] = {0.34, 0.37, 0.39, 0.40};
inline int PB(double x) { return (int)std::lround((x - kPPlo) / kPPw); }
scan::StepCfg Cfg(int i) { scan::StepCfg c; c.mpLo = kSlLo[i]; c.mpHi = kSlHi[i]; c.halfRange = kSlLo[i] >= 0.39 - 1e-9 ? 0.035 : 0.045; return c; }

std::vector<double> SYg(int sc, int g, int cls, int pa, int pb, const std::vector<double>& A) {
    std::vector<double> y(kNM4, 0.);
    for (int q = pa; q < pb; ++q) { const double* a = &A[IdxSS(sc, g, cls, q, 0)]; for (int i = 0; i < kNM4; ++i) y[i] += a[i]; }
    return y;
}
std::vector<double> DC(int cls, int k, int pa, int pb) {
    std::vector<double> y(kNM4, 0.);
    for (int q = pa; q < pb; ++q) { const double* a = &DN[IdxDN(dS0, cls, q, k, 0)]; for (int i = 0; i < kNM4; ++i) y[i] += a[i]; }
    return y;
}

// quantities (copy of bedge_ana.C Quant, subset)
double CF[3][3] = {{0.880, 0.889, 0.868}, {0.891, 0.869, 0.896}, {0, 0, 0}};
const double kKapPeak = 0.48;
const int kNQ = 8;
const char* const kQN[kNQ] = {"f_all", "f_E1", "f_E2", "fc_rec", "a_k48", "a_D", "kappa", "a_k48_unc"};
struct Inp { double Ed[3], Es[3], lin[3][3]; };
void Quant(const Inp& x, double* q) {
    for (int k = 0; k < 3; ++k) q[k] = x.Ed[k] / x.Es[k];
    double sC2[3], sC3[3];
    for (int k = 0; k < 3; ++k) { sC2[k] = x.lin[k][1] / x.lin[k][0]; sC3[k] = x.lin[k][2] / x.lin[k][0]; }
    const double red = ((1 - CF[2][1]) * x.lin[1][1] + (1 - CF[2][2]) * x.lin[2][1]) / (x.lin[1][0] + x.lin[2][0]);
    q[3] = x.Ed[0] / (x.Es[0] * (1 - red));
    const double t = x.lin[1][0] + x.lin[2][0];
    const double c2 = (CF[2][1] * x.lin[1][1] + CF[2][2] * x.lin[2][1]) / (x.lin[1][1] + x.lin[2][1]);
    const double s2 = (x.lin[1][1] + x.lin[2][1]) / t, s3 = (x.lin[1][2] + x.lin[2][2]) / t;
    q[4] = x.Ed[0] / (x.Es[0] * (c2 * s2 + kKapPeak * s3 + 1 - s2 - s3));
    q[7] = x.Ed[0] / (x.Es[0] * (s2 + kKapPeak * s3 + 1 - s2 - s3));
    const double A1 = x.Es[1] * (CF[2][1] * sC2[1] + 1 - sC2[1] - sC3[1]), B1 = x.Es[1] * sC3[1];
    const double A2 = x.Es[2] * (CF[2][2] * sC2[2] + 1 - sC2[2] - sC3[2]), B2 = x.Es[2] * sC3[2];
    const double rho = x.Ed[1] / x.Ed[2];
    const double kap = (rho * A2 - A1) / (B1 - rho * B2);
    q[5] = x.Ed[2] / (A2 + kap * B2); q[6] = kap;
}

struct Scen { const char* tag; const char* text; double a, kap, c2, rho, g, c, b, q; int bkg; double kB, dipA, dipM, dipW; int etaCfg; int raw; };
const int kNScen = 16;
const Scen kScen[kNScen] = {
    {"a", "WP4 (a) with the raw simulated background: x 1.35, tilt -0.15, a = 0.4", 0.4, 1, 1, 1, 1, 1.35, -0.15, 0, 0, 1, 0, 0, 0, sNom, 1},
    {"c1", "raw SMASH background x 1, a = 0.4", 0.4, 1, 1, 1, 1, 1.0, 0, 0, 0, 1, 0, 0, 0, sNom, 1},
    {"a1", "raw SMASH background x 1, a = 1 (exact closure expected)", 1.0, 1, 1, 1, 1, 1.0, 0, 0, 0, 1, 0, 0, 0, sNom, 1},
    {"null", "(a) with a = 0", 0.0, 1, 1, 1, 1, 1.35, -0.15, 0, 0, 1, 0, 0, 0, sNom, 1},
    {"sA", "stage-A background at the F-main eta (c 0.98, b -0.064, kappa_B 2.93; a 0.382, kappa 0.832, rho 1.315, g 1.717)", 0.382, 0.832, 1, 1.315, 1.717, 0.98, -0.0639, 0, 2, 2.932, 0, 0, 0, sNom, 1},
    {"q-", "(a) + quadratic term q = -0.173 (stage-A value that moves B4 by +150 %)", 0.4, 1, 1, 1, 1, 1.35, -0.15, -0.173, 0, 1, 0, 0, 0, sNom, 1},
    {"q--", "(a) + q = -0.346", 0.4, 1, 1, 1, 1, 1.35, -0.15, -0.346, 0, 1, 0, 0, 0, sNom, 1},
    {"q+", "(a) + q = +0.173", 0.4, 1, 1, 1, 1, 1.35, -0.15, 0.173, 0, 1, 0, 0, 0, sNom, 1},
    {"om", "WP4 (e): omega/rho0-pion background x 0.3, the rest x 1.6", 0.4, 1, 1, 1, 1, 1.0, -0.15, 0, 1, 1, 0, 0, 0, sNom, 1},
    {"wide", "(a) x (1 + 0.10 Gauss(0.490, 15 MeV)) under the bump", 0.4, 1, 1, 1, 1, 1.35, -0.15, 0, 0, 1, 0.10, 0.490, 0.015, sNom, 1},
    {"narrow", "(a) x (1 + 0.05 Gauss(0.505, 8 MeV)) at the edges", 0.4, 1, 1, 1, 1, 1.35, -0.15, 0, 0, 1, 0.05, 0.505, 0.008, sNom, 1},
    {"smooth", "(a) with the SMOOTHED background (exp(pol3), WP4 style) instead of the raw one", 0.4, 1, 1, 1, 1, 1.35, -0.15, 0, 0, 1, 0, 0, 0, sNom, 0},
    {"mass", "(a) with the eta truth at s + 0.0009, r = 0.89 (templates nominal)", 0.4, 1, 1, 1, 1, 1.35, -0.15, 0, 0, 1, 0, 0, 0, sF, 1},
    {"qed", "(a) with C2 x 0.8855 in the truth (QED gamma* angle)", 0.4, 1, 0.8855, 1, 1, 1.35, -0.15, 0, 0, 1, 0, 0, 0, sNom, 1},
    {"kap", "(a) with conversions (C3, C4) x 0.48 in the truth", 0.4, 0.48, 1, 1, 1, 1.35, -0.15, 0, 0, 1, 0, 0, 0, sNom, 1},
    {"datalike", "a = 0.55, C2 x 0.8855, conversions x 0.48, raw SMASH background x 1", 0.55, 0.48, 0.8855, 1, 1, 1.0, 0, 0, 0, 1, 0, 0, 0, sNom, 1}};
}  // namespace bedget

void bedge_toys(int nToy = 1000) {
    using namespace bedget;
    TStopwatch total;
    gL = fopen(Form("%s/bedge_toys.txt", kOut), "w");
    if (!ReadArr(Form("%s/data_N.bin", kOut), DN, SizeDN()) || !ReadArr(Form("%s/sim_S.bin", kOut), SS, SizeSS()) ||
        !ReadArr(Form("%s/sim_V.bin", kOut), SV, SizeSS()) || !ReadArr(Form("%s/sim_R.bin", kOut), SR, SizeSR())) return;
    for (int k = 0; k < 3; ++k) CF[2][k] = 0.5 * (CF[0][k] + CF[1][k]);
    PL("# bedge_toys.C, 2026-11-14: closure of the edge estimator on pseudo-data, %d toys per scenario\n", nToy);

    // ---------------------------------------------------------------- inputs per slice and class (E1 = 0, E2 = 1)
    struct In { std::vector<double> cb, pp, mm, Som, Snom, SGD, SGC, Rom, Rnom, RGD, RGC, eta[2][5]; };
    std::vector<std::array<In, 2>> in(kNSl);
    std::vector<double> edges; for (int i = 0; i <= 120; ++i) edges.push_back(0.40 + 0.0025 * i);
    const int iSm0 = (int)std::lround((0.40 - scan::kMlo) / 0.0025);
    for (int s = 0; s < kNSl; ++s) {
        const int pa = PB(kSlLo[s]), pb = PB(kSlHi[s]);
        std::vector<double> au(kNM4, 0.), ap(kNM4, 0.), am(kNM4, 0.);
        for (int e = 0; e < 2; ++e) {
            const std::vector<double> u = DC(e, 0, pa, pb), p = DC(e, 1, pa, pb), m = DC(e, 2, pa, pb);
            for (int i = 0; i < kNM4; ++i) { au[i] += u[i]; ap[i] += p[i]; am[i] += m[i]; }
        }
        const scan::Spectra all = MakeSpec(au.data(), ap.data(), am.data(), kNM4, scan::kMlo, kBwD);
        for (int e = 0; e < 2; ++e) {
            In& I = in[s][e];
            const std::vector<double> u = DC(e, 0, pa, pb);
            I.pp = DC(e, 1, pa, pb); I.mm = DC(e, 2, pa, pb);
            const scan::Spectra c = SpreadCB(all, MakeSpec(u.data(), I.pp.data(), I.mm.data(), kNM4, scan::kMlo, kBwD));
            I.cb = c.cb;
            // smoothed non-eta background (omega/rho0 pions; the rest | G_D; G_C)
            for (int part = 0; part < 4; ++part) {
                const int g1 = part == 0 ? gGDom : part == 1 ? gGDnom : part == 2 ? gGDom : gGCom;
                const int g2 = part == 0 ? gGCom : part == 1 ? gGCnom : part == 2 ? gGDnom : gGCnom;
                const std::vector<double> y1 = SYg(sNom, g1, e, pa, pb, SS), y2 = SYg(sNom, g2, e, pa, pb, SS);
                const std::vector<double> v1 = SYg(sNom, g1, e, pa, pb, SV), v2 = SYg(sNom, g2, e, pa, pb, SV);
                wp4::Vec cc(122, 0.), vv(122, 0.);
                for (int i = 1; i <= 120; ++i) { cc[i] = y1[iSm0 + i - 1] + y2[iSm0 + i - 1]; vv[i] = v1[iSm0 + i - 1] + v2[iSm0 + i - 1]; }
                const wp4::Smooth S = wp4::SmoothFit(cc, vv, edges);
                std::vector<double>& out = part == 0 ? I.Som : part == 1 ? I.Snom : part == 2 ? I.SGD : I.SGC;
                std::vector<double>& rawv = part == 0 ? I.Rom : part == 1 ? I.Rnom : part == 2 ? I.RGD : I.RGC;
                rawv.assign(kNM4, 0.);
                for (int i = 0; i < kNM4; ++i) rawv[i] = y1[i] + y2[i];
                out.assign(kNM4, 0.);
                double raw = 0, sm = 0;
                for (int i = 0; i < 120; ++i) { out[iSm0 + i] = S.Int(edges[i], edges[i + 1]); raw += cc[i + 1]; sm += out[iSm0 + i]; }
                if (s == 0 || s == 3) PL("   slice %d class E%d %s: smoothing order %d, N_eff %.0f, raw %.1f, smooth %.1f\n", s, e + 1, part == 0 ? "om" : part == 1 ? "rest" : part == 2 ? "G_D" : "G_C", S.order, S.neff, raw, sm);
            }
            // eta classes, raw: [cfg 0 nominal, 1 sF][C2, C3, C4, C5, G_eta]
            for (int t = 0; t < 2; ++t) {
                const int sc = t == 0 ? sNom : sF;
                const int gg[5] = {gC2, gC3, gC4, gC5, gGeta};
                for (int j = 0; j < 5; ++j) I.eta[t][j] = SYg(sc, gg[j], e, pa, pb, SS);
            }
        }
    }
    // ---------------------------------------------------------------- E_sim (full) and replicas, linear decomposition with the data covariance
    Inp base;
    std::vector<std::vector<double>> Vd(3 * kNSl);
    for (int k = 0; k < 3; ++k) {
        base.Ed[k] = base.Es[k] = 0; for (int g = 0; g < 3; ++g) base.lin[k][g] = 0;
        for (int s = 0; s < kNSl; ++s) {
            const int pa = PB(kSlLo[s]), pb = PB(kSlHi[s]);
            std::vector<double> au(kNM4, 0.), ap(kNM4, 0.), am(kNM4, 0.);
            for (int e = 0; e < 3; ++e) {
                const std::vector<double> u = DC(e, 0, pa, pb), p = DC(e, 1, pa, pb), m = DC(e, 2, pa, pb);
                for (int i = 0; i < kNM4; ++i) { au[i] += u[i]; ap[i] += p[i]; am[i] += m[i]; }
            }
            scan::Spectra d = MakeSpec(au.data(), ap.data(), am.data(), kNM4, scan::kMlo, kBwD);
            if (k > 0) {
                const std::vector<double> u = DC(k - 1, 0, pa, pb), p = DC(k - 1, 1, pa, pb), m = DC(k - 1, 2, pa, pb);
                d = SpreadCB(d, MakeSpec(u.data(), p.data(), m.data(), kNM4, scan::kMlo, kBwD));
            }
            const scan::StepFit fd = FitStep(d, Cfg(s), nullptr, nullptr, &Vd[k * kNSl + s]);
            base.Ed[k] += fd.h;
            std::vector<double> yT(kNM4, 0.), y2(kNM4, 0.), y3(kNM4, 0.);
            for (int e = 0; e < 3; ++e) {
                if (k > 0 && e != k - 1) continue;
                for (int g = 0; g < kNG; ++g) {
                    const std::vector<double> y = SYg(sNom, g, e, pa, pb, SS);
                    for (int i = 0; i < kNM4; ++i) { yT[i] += y[i]; if (g == gC2) y2[i] += y[i]; if (g == gC3) y3[i] += y[i]; }
                }
            }
            base.Es[k] += FitStep(SimSpec(yT.data(), kNM4, scan::kMlo, kBwD), Cfg(s)).h;
            base.lin[k][0] += FitStep(SimSpec(yT.data(), kNM4, scan::kMlo, kBwD), Cfg(s), nullptr, &Vd[k * kNSl + s]).h;
            base.lin[k][1] += FitStep(SimSpec(y2.data(), kNM4, scan::kMlo, kBwD), Cfg(s), nullptr, &Vd[k * kNSl + s]).h;
            base.lin[k][2] += FitStep(SimSpec(y3.data(), kNM4, scan::kMlo, kBwD), Cfg(s), nullptr, &Vd[k * kNSl + s]).h;
        }
    }
    PL("   data: E %.2f / %.2f / %.2f; E_sim %.2f / %.2f / %.2f (all / E1 / E2)\n", base.Ed[0], base.Ed[1], base.Ed[2], base.Es[0], base.Es[1], base.Es[2]);
    std::vector<Inp> rep(kNRepS, base);
    for (int r = 0; r < kNRepS; ++r)
        for (int k = 0; k < 3; ++k) {
            Inp& x = rep[r];
            x.Es[k] = 0; for (int g = 0; g < 3; ++g) x.lin[k][g] = 0;
            for (int s = 0; s < kNSl; ++s) {
                const int pa = PB(kSlLo[s]), pb = PB(kSlHi[s]);
                for (int g = 0; g < 3; ++g) {
                    std::vector<double> y(kNRM4, 0.);
                    for (int e = 0; e < kNCls; ++e) {
                        if (k > 0 && e != k - 1) continue;
                        for (int q = pa; q < pb; ++q) for (int i = 0; i < kNRM4; ++i) y[i] += SR[IdxSR(sNom, g, e, q - kRPPlo, i) + r];
                    }
                    const scan::Spectra sp = SimSpec(y.data(), kNRM4, kRLo, kBwD);
                    if (g == 0) x.Es[k] += FitStep(sp, Cfg(s)).h;
                    x.lin[k][g] += FitStep(sp, Cfg(s), nullptr, &Vd[k * kNSl + s]).h;
                }
            }
        }
    // relative MC spread of every quantity (data fixed)
    double q0[kNQ], mcRel[kNQ];
    Quant(base, q0);
    for (int j = 0; j < kNQ; ++j) {
        std::vector<double> v;
        for (int r = 0; r < kNRepS; ++r) { double q[kNQ]; Quant(rep[r], q); v.push_back(q[j]); }
        std::sort(v.begin(), v.end());
        mcRel[j] = 0.5 * (v[84] - v[15]) / std::fabs(q0[j]);
    }
    PL("   data quantities:"); for (int j = 0; j < kNQ; ++j) PL(" %s %.4f (MC68 %.3f)", kQN[j], q0[j], mcRel[j]); PL("\n");

    // ---------------------------------------------------------------- toys
    FILE* fo = fopen(Form("%s/toys_results.txt", kOut), "w");
    fprintf(fo, "# scenario mode quantity truth mean bias relbias rms meanpull pullsd cover1 n\n");
    TRandom3 rnd(20261114);
    for (int sc = 0; sc < kNScen; ++sc) {
        const Scen& S = kScen[sc];
        TStopwatch sw;
        // truth per slice and class
        std::vector<std::array<std::vector<double>, 2>> mu(kNSl);
        for (int s = 0; s < kNSl; ++s)
            for (int e = 0; e < 2; ++e) {
                const In& I = in[s][e];
                std::vector<double>& m = mu[s][e];
                m.assign(kNM4, 0.);
                const int t = S.etaCfg == sF ? 1 : 0;
                for (int i = 0; i < kNM4; ++i) {
                    const double x = scan::kMlo + (i + 0.5) * kBwD, tx = TiltX(x);
                    const std::vector<double>& Bom = S.raw ? I.Rom : I.Som; const std::vector<double>& Bnom = S.raw ? I.Rnom : I.Snom;
                    const std::vector<double>& BGD = S.raw ? I.RGD : I.SGD; const std::vector<double>& BGC = S.raw ? I.RGC : I.SGC;
                    double bg = S.bkg == 0 ? S.c * (Bom[i] + Bnom[i]) : S.bkg == 1 ? (0.3 * Bom[i] + 1.6 * Bnom[i]) : S.c * (BGD[i] + S.kB * BGC[i]);
                    bg *= (1 + S.b * tx + S.q * tx * tx);
                    if (S.dipA != 0) bg *= 1 + S.dipA * std::exp(-0.5 * (x - S.dipM) * (x - S.dipM) / (S.dipW * S.dipW));
                    const double et = S.a * (S.c2 * I.eta[t][0][i] + S.kap * (I.eta[t][1][i] + I.eta[t][2][i]) + S.rho * I.eta[t][3][i] + S.g * I.eta[t][4][i]);
                    m[i] = std::max(0., I.cb[i] + bg + et);
                }
            }
        // Asimov diagnostic (scenario a only): fit of the expectation and its linear decomposition with the Asimov covariance
        if (sc == 0) {
            double Ea = 0, Eeta = 0, Ebkg = 0, Ecb = 0, Esv = 0;
            for (int s = 0; s < kNSl; ++s) {
                std::vector<double> au(kNM4, 0.), ap(kNM4, 0.), am(kNM4, 0.), ye(kNM4, 0.), yb(kNM4, 0.), yt(kNM4, 0.);
                for (int e = 0; e < 2; ++e) {
                    const In& I = in[s][e];
                    for (int i = 0; i < kNM4; ++i) {
                        const double x = scan::kMlo + (i + 0.5) * kBwD, tx = TiltX(x);
                        const double bg = S.c * (I.Som[i] + I.Snom[i]) * (1 + S.b * tx);
                        const double et = S.a * (I.eta[0][0][i] + I.eta[0][1][i] + I.eta[0][2][i] + I.eta[0][3][i] + I.eta[0][4][i]);
                        au[i] += mu[s][e][i]; ap[i] += I.pp[i]; am[i] += I.mm[i]; ye[i] += et; yb[i] += bg;
                    }
                }
                for (int g = 0; g < kNG; ++g) for (int e = 0; e < 3; ++e) { const std::vector<double> y = SYg(sNom, g, e, PB(kSlLo[s]), PB(kSlHi[s]), SS); for (int i = 0; i < kNM4; ++i) yt[i] += y[i]; }
                const scan::Spectra d = MakeSpec(au.data(), ap.data(), am.data(), kNM4, scan::kMlo, kBwD);
                std::vector<double> V;
                const scan::StepFit f = FitStep(d, Cfg(s), nullptr, nullptr, &V);
                Ea += f.h;
                Eeta += FitStep(SimSpec(ye.data(), kNM4, scan::kMlo, kBwD), Cfg(s), nullptr, &V).h;
                Ebkg += FitStep(SimSpec(yb.data(), kNM4, scan::kMlo, kBwD), Cfg(s), nullptr, &V).h;
                std::vector<double> r(kNM4); for (int i = 0; i < kNM4; ++i) r[i] = d.sig[i] - ye[i] - yb[i];
                Ecb += FitStep(SimSpec(r.data(), kNM4, scan::kMlo, kBwD), Cfg(s), nullptr, &V).h;
                Esv += FitStep(SimSpec(yt.data(), kNM4, scan::kMlo, kBwD), Cfg(s), nullptr, &V).h;
            }
            PL("   ASIMOV (scenario a: raw background x 1.35, tilt -0.15, a = 0.4): E = %.2f = eta %.2f + smoothed background %.2f + (raw - smoothed background, CB) %.2f; E_sim with the Asimov covariance %.2f -> f %.4f (WP6 weights: %.4f)\n",
               Ea, Eeta, Ebkg, Ecb, Esv, Ea / Esv, Ea / base.Es[0]);
        }
        std::vector<std::array<double, kNQ>> qa, qr, ea;
        for (int it = 0; it < nToy; ++it) {
            Inp x = base;
            double vEd[3] = {0, 0, 0};
            for (int k = 0; k < 3; ++k) x.Ed[k] = 0;
            for (int s = 0; s < kNSl; ++s) {
                std::vector<double> u[2], p[2], m[2];
                for (int e = 0; e < 2; ++e) {
                    u[e].resize(kNM4); p[e].resize(kNM4); m[e].resize(kNM4);
                    for (int i = 0; i < kNM4; ++i) {
                        u[e][i] = mu[s][e][i] > 0 ? rnd.Poisson(mu[s][e][i]) : 0.;
                        p[e][i] = in[s][e].pp[i] > 0 ? rnd.Poisson(in[s][e].pp[i]) : 0.;
                        m[e][i] = in[s][e].mm[i] > 0 ? rnd.Poisson(in[s][e].mm[i]) : 0.;
                    }
                }
                std::vector<double> au(kNM4), ap(kNM4), am(kNM4);
                for (int i = 0; i < kNM4; ++i) { au[i] = u[0][i] + u[1][i]; ap[i] = p[0][i] + p[1][i]; am[i] = m[0][i] + m[1][i]; }
                const scan::Spectra all = MakeSpec(au.data(), ap.data(), am.data(), kNM4, scan::kMlo, kBwD);
                for (int k = 0; k < 3; ++k) {
                    const scan::Spectra d = k == 0 ? all : SpreadCB(all, MakeSpec(u[k - 1].data(), p[k - 1].data(), m[k - 1].data(), kNM4, scan::kMlo, kBwD));
                    const scan::StepFit f = FitStep(d, Cfg(s));
                    x.Ed[k] += f.h; vEd[k] += f.eh * f.eh;
                }
            }
            std::array<double, kNQ> q, qq, e;
            Quant(x, q.data());
            Inp xr = rep[it % kNRepS]; for (int k = 0; k < 3; ++k) xr.Ed[k] = x.Ed[k];
            Quant(xr, qq.data());
            // fit error propagated
            for (int j = 0; j < kNQ; ++j) e[j] = 0;
            for (int k = 0; k < 3; ++k) {
                Inp y = x; y.Ed[k] += std::sqrt(vEd[k]);
                double qd[kNQ]; Quant(y, qd);
                for (int j = 0; j < kNQ; ++j) e[j] += (qd[j] - q[j]) * (qd[j] - q[j]);
            }
            for (int j = 0; j < kNQ; ++j) e[j] = std::sqrt(e[j]);
            qa.push_back(q); qr.push_back(qq); ea.push_back(e);
        }
        PL("\nscenario %-8s %s (%.0f s)\n", S.tag, S.text, sw.RealTime());
        PL("   %-7s %-4s %8s %9s %8s %8s %8s %8s %6s\n", "q", "mode", "truth", "mean", "relbias", "rms/tr", "pullm", "pullsd", "cov");
        for (int mode = 0; mode < 2; ++mode)
            for (int j = 0; j < kNQ; ++j) {
                const double tr = S.a;
                double m = 0, m2 = 0, pm = 0, p2 = 0; int cov = 0, n = 0;
                for (int it = 0; it < nToy; ++it) {
                    const double v = mode == 0 ? qa[it][j] : qr[it][j];
                    if (!std::isfinite(v)) continue;
                    const double sig = std::sqrt(ea[it][j] * ea[it][j] + (mode == 1 ? std::pow(mcRel[j] * v, 2) : 0.));
                    const double pull = (v - tr) / sig;
                    m += v; m2 += v * v; pm += pull; p2 += pull * pull; if (std::fabs(pull) < 1) ++cov; ++n;
                }
                m /= n; const double rms = std::sqrt(std::max(0., m2 / n - m * m)); pm /= n; const double psd = std::sqrt(std::max(0., p2 / n - pm * pm));
                const double rb = tr != 0 ? (m - tr) / tr : m - tr;
                if (j == 0 || j == 3 || j == 4 || j == 7 || mode == 0)
                    PL("   %-7s %-4s %8.4f %9.5f %+8.4f %8.4f %+8.3f %8.3f %6.3f\n", kQN[j], mode == 0 ? "AA" : "RF", tr, m, rb, tr != 0 ? rms / tr : rms, pm, psd, (double)cov / n);
                fprintf(fo, "%s %s %s %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %d\n", S.tag, mode == 0 ? "AA" : "RF", kQN[j], tr, m, m - tr, rb, rms, pm, psd, (double)cov / n, n);
            }
    }
    fclose(fo);
    PL("\ndone in %.0f s\n", total.RealTime());
    fclose(gL);
}
