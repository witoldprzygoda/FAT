// ============================================================================
// wp2_vzfit.C  --  WP2: cross-check of the flux ratio b (full/empty) by a fit of the full-target vertex-z profile
// ============================================================================
// Plain words: the full-target vertex distribution of the channel is the sum of the LH2 cell (shape from the
// simulation, which contains only the LH2 cell) and of the material (shape = the measured empty-target
// distribution of the same day). A fit F(z) = a x S(z; delta, sigma) + b x E(z) gives the material scale b without
// any normaliser window. S is the simulated reconstructed vz of T12 (S0 without the vz cut, wp2_simvz.root), shifted
// by delta and smeared by an extra Gaussian sigma (data resolution may differ). Rows: unlike-sign, base selection +
// both pions in the vertex fit, all M4 (wp2_win.root, 2 mm bins). Fit range [-250,-60) mm (variants: [-220,-70),
// sigma fixed to 0). chi2 with the variance F + b^2 E (+ a^2 var S).
// Then f_win(fit) = b_fit E(vz > -160) / F(vz > -160) for these rows, to compare with b of the normaliser window.
// Output: docs/latex/tables_all/wp2_vzfit.txt, /hdd3/.../wp2/wp2_vzfit.root (data, fitted components per day)
// Run: root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_wp2", kTRUE)' \
//           'plots/all_report/wp2_vzfit.C+O'
// ============================================================================
#include <cmath>
#include <cstdio>
#include <functional>
#include <vector>

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "TFile.h"
#include "TH1D.h"

namespace wp2vf {
const int kDays[3] = {34, 47, 65};
TH1D* Smear(TH1D* s, double delta, double sigma, const char* name) {
    TH1D* o = (TH1D*)s->Clone(name);
    o->Reset();
    o->SetDirectory(nullptr);
    const int n = s->GetNbinsX();
    const double w = s->GetBinWidth(1);
    for (int j = 1; j <= n; ++j) {
        const double c = s->GetBinContent(j), v = std::pow(s->GetBinError(j), 2);
        if (c == 0) continue;
        const double z0 = s->GetBinCenter(j) + delta;
        if (sigma < 0.05) {
            // shift only: split between the two neighbouring bins
            const double x = (z0 - s->GetBinCenter(1)) / w + 1;
            const int i0 = (int)std::floor(x);
            const double fr = x - i0;
            if (i0 >= 1 && i0 <= n) { o->AddBinContent(i0, c * (1 - fr)); o->SetBinError(i0, std::sqrt(std::pow(o->GetBinError(i0), 2) + v * (1 - fr) * (1 - fr))); }
            if (i0 + 1 >= 1 && i0 + 1 <= n) { o->AddBinContent(i0 + 1, c * fr); o->SetBinError(i0 + 1, std::sqrt(std::pow(o->GetBinError(i0 + 1), 2) + v * fr * fr)); }
            continue;
        }
        const int k1 = std::max(1, s->FindBin(z0 - 5 * sigma)), k2 = std::min(n, s->FindBin(z0 + 5 * sigma));
        for (int i = k1; i <= k2; ++i) {
            const double lo = s->GetBinLowEdge(i), hi = lo + w;
            const double p = 0.5 * (std::erf((hi - z0) / (std::sqrt(2.) * sigma)) - std::erf((lo - z0) / (std::sqrt(2.) * sigma)));
            o->AddBinContent(i, c * p);
            o->SetBinError(i, std::sqrt(std::pow(o->GetBinError(i), 2) + v * p * p));
        }
    }
    return o;
}
}  // namespace wp2vf

void wp2_vzfit() {
    using namespace wp2vf;
    TFile* fw = TFile::Open("/hdd3/przygoda/HADES/FAT_out/all_report/wp2/wp2_win.root", "READ");
    TFile* fs = TFile::Open("/hdd3/przygoda/HADES/FAT_out/all_report/wp2/wp2_simvz.root", "READ");
    TH1D* S = (TH1D*)fs->Get("simvz_T12");
    S->SetDirectory(nullptr);
    TFile* fo = TFile::Open("/hdd3/przygoda/HADES/FAT_out/all_report/wp2/wp2_vzfit.root", "RECREATE");
    FILE* o = fopen("/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all/wp2_vzfit.txt", "w");
    fprintf(o, "# wp2_vzfit.C: F(z) = a S(z; delta, sigma) + b E(z), unlike-sign rows, base + pion vertex, all M4, 2 mm bins\n");
    fprintf(o, "# S = simulated LH2 (T12) reconstructed vz; E = empty-target vz of the same day; b = material scale (flux ratio full/empty)\n");
    fprintf(o, "# f_fit = b E(vz>-160)/F(vz>-160) (these rows); b_U = normaliser-window value ch[-200,-170) from wp2_win.txt for comparison\n");
    fprintf(o, "# variant  day  range  b +- err  a  delta[mm]  sigma[mm]  chi2/ndf  f_fit[%%]  F(vz>-160)  E(vz>-160)\n");
    struct Var { const char* name; double lo, hi; bool fixSig; };
    const Var vars[3] = {{"nominal", -250, -60, false}, {"narrow", -220, -70, false}, {"sigma0", -250, -60, true}};
    for (int d = 0; d < 3; ++d) {
        TH1D* F = (TH1D*)fw->Get(Form("vz_full_us_%03d", kDays[d]));
        TH1D* E = (TH1D*)fw->Get(Form("vz_empty_us_%03d", kDays[d]));
        if (!F || !E) { printf("wp2_vzfit: missing histograms for day %03d\n", kDays[d]); continue; }
        const int c1 = F->FindBin(-159.999);
        const double Fin = F->Integral(c1, F->GetNbinsX()), Ein = E->Integral(c1, E->GetNbinsX());
        const double sscale = F->Integral(F->FindBin(-139.999), F->FindBin(-110.001)) / S->Integral(S->FindBin(-139.999), S->FindBin(-110.001));
        for (const Var& V : vars) {
            const int i1 = F->FindBin(V.lo + 1e-3), i2 = F->FindBin(V.hi - 1e-3);
            int ndf = 0;
            auto chi2 = [&](const double* p) -> double {
                TH1D* Sp = Smear(S, p[2], V.fixSig ? 0. : std::fabs(p[3]), "tmpS");
                double c = 0;
                ndf = 0;
                for (int i = i1; i <= i2; ++i) {
                    const double f = F->GetBinContent(i), e = E->GetBinContent(i);
                    const double m = p[0] * sscale * Sp->GetBinContent(i) + p[1] * e;
                    const double var = std::max(1.0, f) + p[1] * p[1] * e + std::pow(p[0] * sscale * Sp->GetBinError(i), 2);
                    c += (f - m) * (f - m) / var;
                    ++ndf;
                }
                delete Sp;
                return c;
            };
            ROOT::Math::Minimizer* mn = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
            ROOT::Math::Functor fun(chi2, 4);
            mn->SetFunction(fun);
            mn->SetPrintLevel(0);
            mn->SetMaxFunctionCalls(20000);
            const double b0 = d == 2 ? 0.9 : 12.;
            mn->SetVariable(0, "a", 1.0, 0.01);
            mn->SetVariable(1, "b", b0, 0.1 * b0);
            mn->SetVariable(2, "delta", 0.0, 0.2);
            if (V.fixSig) mn->SetFixedVariable(3, "sigma", 0.0);
            else mn->SetLimitedVariable(3, "sigma", 2.0, 0.2, 0.0, 15.0);
            mn->Minimize();
            mn->Hesse();
            const double* p = mn->X();
            const double* e = mn->Errors();
            const double c2 = chi2(p);
            const int nfree = V.fixSig ? 3 : 4;
            const double ff = p[1] * Ein / Fin;
            fprintf(o, "  %-8s %03d [%4.0f,%4.0f)  b %7.3f +- %6.3f  a %.4f  delta %6.2f  sigma %5.2f  chi2/ndf %7.1f/%d  f_fit %6.2f +- %5.2f  %8.0f %7.0f\n",
                    V.name, kDays[d], V.lo, V.hi, p[1], e[1], p[0], p[2], std::fabs(p[3]), c2, ndf - nfree, 100 * ff,
                    100 * ff * std::sqrt(std::pow(e[1] / p[1], 2) + 1 / std::max(1., Ein)), Fin, Ein);
            if (V.name[0] == 'n' && V.name[1] == 'o') {
                fo->cd();
                TH1D* Sp = Smear(S, p[2], std::fabs(p[3]), Form("fit_lh2_%03d", kDays[d]));
                Sp->Scale(p[0] * sscale);
                Sp->Write();
                TH1D* Eb = (TH1D*)E->Clone(Form("fit_mat_%03d", kDays[d]));
                Eb->Scale(p[1]);
                Eb->Write();
                F->Clone(Form("data_full_%03d", kDays[d]))->Write();
            }
            delete mn;
        }
    }
    fclose(o);
    fo->Close();
    printf("wp2_vzfit: done; table docs/latex/tables_all/wp2_vzfit.txt\n");
    FILE* in = fopen("/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all/wp2_vzfit.txt", "r");
    char line[512];
    while (fgets(line, sizeof line, in)) printf("%s", line);
    fclose(in);
}
