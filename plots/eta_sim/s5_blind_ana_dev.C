// ============================================================================
// s5_blind_ana_dev.C  --  method development of the blind analysis on the
// half-A TEMPLATES only (Asimov = template sum, no pseudo-data are read)
// ============================================================================
// Label blind_ana. Reads only blind/templates/templates.root. For every frozen
// selection (PEAK, DIRECT, BUMP, X17SEL) and the eta toys' luminosity
// (scale 0.5 = L_sim / 2):
//   A. free Gaussian + pol2 fit (position, width) on the Asimov sum
//   B. M1 (Gaussian frozen 0.5460 / 6.7 MeV + pol2, 0.526-0.614, 2 MeV bins):
//      Asimov closure and expected error
//   C. M3 (linear side bands 0.5255-0.5333 and 0.5587-0.5665): closure
//   D. bump shape of the half-A bump23 template (Gaussian, bifurcated Gaussian,
//      erfc edge) and the Asimov closure of the bump fits B1 and B3
// Output: printed numbers (log in the scan directory).
// Run: root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/s5_blind_ana_build", kTRUE)' plots/eta_sim/s5_blind_ana_dev.C+
// ============================================================================
#include <cmath>
#include <cstdio>
#include <vector>

#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"

#include "s5_blind_ana_fit.h"
#include "s5_blind_ana_tmpl.h"

using namespace bla;

namespace dev {

const int kRb = 2;             // 2 MeV bins
const double kGridLo = 0.30, kBw = 0.002;
inline double X(int i) { return kGridLo + (i + 0.5) * kBw; }

// peak fit problem: Gaussian (bin fraction) + polynomial of order `ord`, range [lo, hi) (bin centres)
LinProblem PeakProblem(const std::vector<double>& y, double lo, double hi, int ord, double mean, double sig, std::vector<int>& idx) {
    LinProblem P; idx.clear();
    for (int i = 0; i < (int)y.size(); ++i) if (X(i) >= lo && X(i) < hi) idx.push_back(i);
    P.n = idx.size(); P.np = ord + 2;
    const double x0 = 0.5 * (lo + hi), dx = 0.5 * (hi - lo);
    for (int q = 0; q < P.n; ++q) {
        const int i = idx[q];
        P.y.push_back(y[i]);
        Vec f(P.np);
        f[0] = GaussFrac(kGridLo + i * kBw, kGridLo + (i + 1) * kBw, mean, sig);
        double t = (X(i) - x0) / dx, tp = 1;
        for (int j = 0; j <= ord; ++j) { f[1 + j] = tp; tp *= t; }
        P.F.push_back(f);
    }
    return P;
}

void FreeFit(const std::vector<double>& y, double lo, double hi, int ord, double& bm, double& bs, double& bS, double& chi2) {
    double best = 1e30;
    for (int im = 0; im <= 60; ++im)
        for (int is = 0; is <= 50; ++is) {
            const double m = 0.5430 + 0.0001 * im, s = 0.0045 + 0.0001 * is;
            std::vector<int> idx;
            LinProblem P = PeakProblem(y, lo, hi, ord, m, s, idx);
            LinResult r = LinFit(P);
            if (r.ok && r.chi2 < best) { best = r.chi2; bm = m; bs = s; bS = r.a[0]; }
        }
    chi2 = best;
}

void Run() {
    TFile f(kTemplFile);
    for (int sel = 0; sel < kNSel; ++sel) {
        M4Templ t;
        if (!LoadM4(f, sel, t)) return;
        const double sc = 0.5;
        std::vector<double> vs;
        const std::vector<double> sum = Sum(t, AtAll(), kRb, sc, &vs);
        const double truePeak = sc * (t.regPeak[kDirect5] + t.regPeak[kConv4]);
        const double trueDir = sc * t.regPeak[kDirect5];
        const double trueBump = sc * (t.regBump[kBump2] + t.regBump[kBump3]);
        printf("\n==================== selection %s (L_toy = L_sim/2) ====================\n", kSel[sel]);
        printf("truth (half A x 0.5): peak45 in peak window %.2f, direct5 %.2f, bump23 in bump window %.2f\n", truePeak, trueDir, trueBump);
        // A. free fit
        double bm, bs, bS, c2;
        FreeFit(sum, 0.526, 0.614, 2, bm, bs, bS, c2);
        printf("A. free Gauss+pol2 0.526-0.614 on Asimov sum: mean %.4f sigma %.2f MeV S %.1f (chi2 %.2f)\n", bm, 1e3 * bs, bS, c2);
        {
            const std::vector<double> pk = Sum(t, AtPeak(), kRb, sc);
            FreeFit(pk, 0.526, 0.614, 0, bm, bs, bS, c2);
            printf("   free Gauss+pol0 on peak45 alone: mean %.4f sigma %.2f MeV S %.1f\n", bm, 1e3 * bs, bS);
        }
        // B. M1 nominal and variations
        {
            std::vector<int> idx;
            LinProblem P = PeakProblem(sum, 0.526, 0.614, 2, kPeakMean, kPeakSigma, idx);
            LinResult r = LinFit(P);
            const double fw = GaussFrac(kPeakLo, kPeakHi, kPeakMean, kPeakSigma);
            printf("B. M1 nominal: S = %.1f +- %.1f (whole Gaussian), S_win = %.1f (truth %.1f, ratio %.3f), chi2 %.2f/%d\n",
                   r.a[0], std::sqrt(r.C[0][0]), r.a[0] * fw, truePeak, r.a[0] * fw / truePeak, r.chi2, r.ndf);
            for (int ord = 1; ord <= 4; ++ord) {
                LinProblem Q = PeakProblem(sum, 0.526, 0.614, ord, kPeakMean, kPeakSigma, idx);
                LinResult q = LinFit(Q);
                printf("   order %d: S_win %.1f ratio %.3f chi2 %.2f\n", ord, q.a[0] * fw, q.a[0] * fw / truePeak, q.chi2);
            }
            const double los[4] = {0.522, 0.518, 0.512, 0.508}, his[4] = {0.594, 0.614, 0.624, 0.654};
            for (int k = 0; k < 4; ++k) {
                LinProblem Q = PeakProblem(sum, los[k], his[k], 2, kPeakMean, kPeakSigma, idx);
                LinResult q = LinFit(Q);
                printf("   range %.3f-%.3f: S_win %.1f ratio %.3f chi2 %.2f/%d\n", los[k], his[k], q.a[0] * fw, q.a[0] * fw / truePeak, q.chi2, q.ndf);
            }
            for (int k = -1; k <= 1; k += 2) {
                LinProblem Q = PeakProblem(sum, 0.526, 0.614, 2, kPeakMean, kPeakSigma + k * 0.0005, idx);
                LinResult q = LinFit(Q);
                const double fw2 = GaussFrac(kPeakLo, kPeakHi, kPeakMean, kPeakSigma + k * 0.0005);
                printf("   width %+.1f MeV: S_win %.1f ratio %.3f\n", k * 0.5, q.a[0] * fw2, q.a[0] * fw2 / truePeak);
            }
        }
        // C. M3 side bands
        {
            const double w = kPeakHi - kPeakLo, sbw = 0.0078;
            const double lL = kPeakLo - sbw, uU = kPeakHi + sbw;
            const double Nw = SumWin(t, AtAll(), kPeakLo, kPeakHi, sc), NL = SumWin(t, AtAll(), lL, kPeakLo, sc), NU = SumWin(t, AtAll(), kPeakHi, uU, sc);
            const double Pw = SumWin(t, AtPeak(), kPeakLo, kPeakHi, sc), PL = SumWin(t, AtPeak(), lL, kPeakLo, sc), PU = SumWin(t, AtPeak(), kPeakHi, uU, sc);
            const double Bw = SumWin(t, {kBump2, kBump3, kGeta, kGpi0D, kGconv, kGvec, kGoth}, kPeakLo, kPeakHi, sc);
            const double BL = SumWin(t, {kBump2, kBump3, kGeta, kGpi0D, kGconv, kGvec, kGoth}, lL, kPeakLo, sc);
            const double BU = SumWin(t, {kBump2, kBump3, kGeta, kGpi0D, kGconv, kGvec, kGoth}, kPeakHi, uU, sc);
            const double rr = w / sbw;
            const double lam = rr * 0.5 * (PL + PU) / Pw;
            const double raw = Nw - rr * 0.5 * (NL + NU);
            printf("C. M3: N_win %.1f N_L %.1f N_U %.1f; peak45 win %.1f L %.2f U %.2f; leakage lambda %.4f\n", Nw, NL, NU, Pw, PL, PU, lam);
            printf("   raw S = %.1f, corrected S = %.1f (truth %.1f, ratio %.3f); background win %.1f vs interpolated %.1f (bias of B %.1f)\n",
                   raw, raw / (1 - lam), Pw, raw / (1 - lam) / Pw, Bw, rr * 0.5 * (BL + BU), rr * 0.5 * (BL + BU) - Bw);
            printf("   expected stat error (Poisson, Asimov): %.1f\n", std::sqrt(Nw + rr * rr * 0.25 * (NL + NU)) / (1 - lam));
        }
        // D. bump shape
        {
            std::vector<double> vb;
            const std::vector<double> bp = Sum(t, AtBump(), kRb, sc, &vb);
            TH1D hb(Form("hb%d", sel), "", 450, 0.30, 1.20);
            for (int i = 0; i < 450; ++i) { hb.SetBinContent(i + 1, bp[i]); hb.SetBinError(i + 1, std::sqrt(vb[i])); }
            TF1 g("g", "[0]*exp(-0.5*((x-[1])/[2])^2)", 0.40, 0.54);
            g.SetParameters(20, 0.495, 0.02);
            hb.Fit(&g, "QNR0");
            TF1 bf("bf", "[0]*(x<[1]?exp(-0.5*((x-[1])/[2])^2):exp(-0.5*((x-[1])/[3])^2))", 0.40, 0.54);
            bf.SetParameters(20, 0.50, 0.03, 0.012);
            hb.Fit(&bf, "QNR0");
            TF1 ed("ed", "[0]*(1+[1]*(x-0.50))*0.5*TMath::Erfc((x-[2])/(sqrt(2.)*[3]))", 0.47, 0.546);
            ed.SetParameters(20, 0, 0.512, 0.006);
            hb.Fit(&ed, "QNR0");
            printf("D. bump23 shape: Gauss mean %.4f sigma %.4f chi2/ndf %.1f/%d | bifur mean %.4f sL %.4f sR %.4f chi2/ndf %.1f/%d | edge e %.4f s %.4f slope %.2f chi2/ndf %.1f/%d\n",
                   g.GetParameter(1), g.GetParameter(2), g.GetChisquare(), g.GetNDF(),
                   bf.GetParameter(1), bf.GetParameter(2), bf.GetParameter(3), bf.GetChisquare(), bf.GetNDF(),
                   ed.GetParameter(2), ed.GetParameter(3), ed.GetParameter(1), ed.GetChisquare(), ed.GetNDF());
            // B1 closure: bump shape + eta Gauss + pol(ord) in [0.400, 0.620)
            for (int shape = 0; shape < 2; ++shape)
                for (int ord = 2; ord <= 4; ++ord) {
                    LinProblem P; std::vector<int> idx;
                    const double lo = 0.400, hi = 0.620, x0 = 0.5 * (lo + hi), dx = 0.5 * (hi - lo);
                    for (int i = 0; i < 450; ++i) if (X(i) >= lo && X(i) < hi) idx.push_back(i);
                    P.n = idx.size(); P.np = ord + 3;
                    for (int i : idx) {
                        P.y.push_back(sum[i]);
                        Vec fv(P.np);
                        const double x1 = kGridLo + i * kBw, x2 = x1 + kBw;
                        fv[0] = shape == 0 ? GaussFrac(x1, x2, g.GetParameter(1), g.GetParameter(2))
                                           : BifurFrac(x1, x2, bf.GetParameter(1), bf.GetParameter(2), bf.GetParameter(3));
                        fv[1] = GaussFrac(x1, x2, kPeakMean, kPeakSigma);
                        double tt = (X(i) - x0) / dx, tp = 1;
                        for (int j = 0; j <= ord; ++j) { fv[2 + j] = tp; tp *= tt; }
                        P.F.push_back(fv);
                    }
                    LinResult r = LinFit(P);
                    const double fwin = shape == 0 ? GaussFrac(kBumpLo, kBumpHi, g.GetParameter(1), g.GetParameter(2))
                                                   : BifurFrac(kBumpLo, kBumpHi, bf.GetParameter(1), bf.GetParameter(2), bf.GetParameter(3));
                    const double Y = r.a[0] * fwin;
                    const double fw = GaussFrac(kPeakLo, kPeakHi, kPeakMean, kPeakSigma);
                    printf("   B1 shape %s pol%d: bump Y_win %.1f +- %.1f (truth %.1f, ratio %.3f); peak S_win %.1f (ratio %.3f); chi2 %.1f/%d\n",
                           shape == 0 ? "gauss" : "bifur", ord, Y, std::sqrt(r.C[0][0]) * fwin, trueBump, Y / trueBump, r.a[1] * fw, r.a[1] * fw / truePeak, r.chi2, r.ndf);
                }
            // B3 edge step: pol2 + eta Gauss + h * step in [0.470, 0.566)
            for (int ord = 1; ord <= 3; ++ord) {
                LinProblem P; std::vector<int> idx;
                const double lo = 0.470, hi = 0.566, x0 = 0.5 * (lo + hi), dx = 0.5 * (hi - lo);
                for (int i = 0; i < 450; ++i) if (X(i) >= lo && X(i) < hi) idx.push_back(i);
                P.n = idx.size(); P.np = ord + 3;
                for (int i : idx) {
                    P.y.push_back(sum[i]);
                    Vec fv(P.np);
                    const double x1 = kGridLo + i * kBw, x2 = x1 + kBw;
                    fv[0] = StepBin(x1, x2, ed.GetParameter(2), ed.GetParameter(3));
                    fv[1] = GaussFrac(x1, x2, kPeakMean, kPeakSigma);
                    double tt = (X(i) - x0) / dx, tp = 1;
                    for (int j = 0; j <= ord; ++j) { fv[2 + j] = tp; tp *= tt; }
                    P.F.push_back(fv);
                }
                LinResult r = LinFit(P);
                printf("   B3 edge step pol%d: h %.2f +- %.2f per 2 MeV, kappa = truth/h = %.1f, rel. error %.3f; chi2 %.1f/%d\n",
                       ord, r.a[0], std::sqrt(r.C[0][0]), trueBump / r.a[0], std::sqrt(r.C[0][0]) / r.a[0], r.chi2, r.ndf);
            }
            // B4 side-band normalised background template
            {
                const double SBlo = 0.5665, SBhi = 0.700;
                const double Tb_w = SumWin(t, AtBkg(), kBumpLo, kBumpHi, sc), Tb_sb = SumWin(t, AtBkg(), SBlo, SBhi, sc);
                const double N_sb = SumWin(t, AtAll(), SBlo, SBhi, sc), N_w = SumWin(t, AtAll(), kBumpLo, kBumpHi, sc);
                const double Pk_w = SumWin(t, AtPeak(), kBumpLo, kBumpHi, sc);
                const double Y = N_w - N_sb * Tb_w / Tb_sb - Pk_w;
                printf("   B4 SB-normalised background: N_win %.1f N_SB %.1f ratio %.4f, Y %.1f (truth %.1f), stat %.1f\n", N_w, N_sb, Tb_w / Tb_sb, Y, trueBump,
                       std::sqrt(N_w + N_sb * std::pow(Tb_w / Tb_sb, 2)));
            }
        }
    }
}

} // namespace dev

void s5_blind_ana_dev() { dev::Run(); }
