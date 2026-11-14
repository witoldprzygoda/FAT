// ============================================================================
// s5_blind_ana_xdev.C  --  development of the X17 fit on the half-A TEMPLATES
// only (no pseudo-data are read). Label blind_ana.
// ============================================================================
// M(ee) templates in the X17 M4 window (pre + both pions in the vertex fit),
// scaled to L_data = 6080.5 nb^-1 (UNCONFIRMED): x 27.5074.
//   1. smooth description of 'other' (= otherX17 + conv4) on 6-60 MeV/c2
//   2. Asimov fits in the frozen window 13.0-26.5 MeV/c2: closure, expected
//      statistical error of BR, effect of the template MC statistics
// Run: root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/s5_blind_ana_build", kTRUE)' plots/eta_sim/s5_blind_ana_xdev.C+
// ============================================================================
#include <cmath>
#include <cstdio>
#include <vector>

#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TH1D.h"
#include "TRandom3.h"

#include "s5_blind_ana_fit.h"

using namespace bla;

namespace xdev {
const double kSc = 6080.5 / 221.05;

void Run() {
    TFile f("/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/blind/templates/templates.root");
    TH1D* D = (TH1D*)f.Get("mee_X17_direct5")->Clone("D"); D->Scale(kSc);
    TH1D* C = (TH1D*)f.Get("mee_X17_conv4")->Clone("C"); C->Scale(kSc);
    TH1D* O = (TH1D*)f.Get("mee_X17_otherX17")->Clone("O"); O->Scale(kSc);
    TH1D* X = (TH1D*)f.Get("mee_X17_x17perBR")->Clone("X"); X->Scale(kSc);
    TH1D* OC = (TH1D*)O->Clone("OC"); OC->Add(C);
    // 1. smooth fits of other + conv
    const char* forms[4] = {"exp([0]+[1]*x+[2]*x*x)", "exp([0]+[1]*x+[2]*x*x+[3]*x*x*x)", "[0]*pow(x,[1])*exp(-[2]*x)", "[0]*exp(-x/[1])+[2]*exp(-x/[3])"};
    for (double lo : {6.0, 8.0, 10.0}) {
        for (int k = 0; k < 4; ++k) {
            TF1 fn(Form("fn%d", k), forms[k], lo, 60.);
            if (k == 0) fn.SetParameters(7.5, -0.08, 0.0);
            if (k == 1) fn.SetParameters(7.5, -0.08, 0.0, 0.0);
            if (k == 2) fn.SetParameters(3000, 0.0, 0.1);
            if (k == 3) fn.SetParameters(3000, 5, 300, 30);
            TFitResultPtr r = OC->Fit(&fn, "QNRS0");
            // integral in window from function (bin centres x 1 bin) and from template
            double sf = 0, st = 0, et2 = 0;
            for (int b = 27; b <= 53; ++b) { sf += fn.Eval(OC->GetBinCenter(b)); st += OC->GetBinContent(b); et2 += std::pow(OC->GetBinError(b), 2); }
            printf("lo %.0f form %d: chi2/ndf %.1f/%d; window: function %.1f, template %.1f +- %.1f\n", lo, k, fn.GetChisquare(), fn.GetNDF(), sf, st, std::sqrt(et2));
        }
    }
    // the X17 shape inside the window
    double xs = 0; for (int b = 27; b <= 53; ++b) xs += X->GetBinContent(b);
    printf("X17 rows per unit BR in window: %.4g; direct5 in window %.1f, conv4 %.1f, other %.1f\n", xs, D->Integral(27, 53), C->Integral(27, 53), O->Integral(27, 53));

    // 2. Asimov fits in the window with smooth other (exp(pol2) from 8-60 MeV)
    TF1 sm("sm", forms[0], 8., 60.); sm.SetParameters(7.5, -0.08, 0.0);
    OC->Fit(&sm, "QNR0");
    std::vector<double> d, o, x, oraw, ovar, y;
    for (int b = 27; b <= 53; ++b) {
        d.push_back(D->GetBinContent(b));
        // smooth: average of the function over the bin
        double s = 0; for (int q = 0; q < 10; ++q) s += sm.Eval(OC->GetBinLowEdge(b) + (q + 0.5) * 0.05); o.push_back(s / 10.);
        oraw.push_back(OC->GetBinContent(b)); ovar.push_back(std::pow(OC->GetBinError(b), 2));
        x.push_back(X->GetBinContent(b));
    }
    auto fitBR = [&](const std::vector<double>& yy, bool withD, LinResult& res) {
        LinProblem P; P.n = yy.size(); P.np = withD ? 3 : 2; P.y = yy;
        for (int i = 0; i < P.n; ++i) { Vec fv; fv.push_back(x[i]); if (withD) fv.push_back(d[i]); fv.push_back(o[i]); P.F.push_back(fv); }
        res = LinFit(P);
        return res.a[0];
    };
    for (double br : {0.0, 5e-6, 1e-5}) {
        std::vector<double> yy(d.size());
        for (size_t i = 0; i < d.size(); ++i) yy[i] = d[i] + o[i] + br * x[i];
        LinResult r; double b = fitBR(yy, true, r);
        printf("Asimov smooth, injected %.1e: BR-hat %.3e +- %.3e (D, O free) chi2 %.2f\n", br, b, std::sqrt(r.C[0][0]), r.chi2);
        LinResult r2; double b2 = fitBR(yy, false, r2);
        printf("   only O free (D in O): BR-hat %.3e +- %.3e\n", b2, std::sqrt(r2.C[0][0]));
        // weights and parent-noise term
        Vec g(3, 0.); g[0] = 1;
        Vec w = r.Weights(g);
        double vp = 0; for (size_t i = 0; i < w.size(); ++i) vp += w[i] * w[i] * ovar[i];
        printf("   sigma from half-A-size MC noise of 'other' (per bin): %.3e\n", std::sqrt(vp));
    }
    {
        std::vector<double> yy(d.size());
        for (size_t i = 0; i < d.size(); ++i) yy[i] = d[i] + oraw[i];
        LinResult r; double b = fitBR(yy, true, r);
        printf("Asimov with RAW half-A other (spikes): BR-hat %.3e +- %.3e chi2 %.1f/%d (this is one realisation of the parent-noise offset)\n", b, std::sqrt(r.C[0][0]), r.chi2P, r.ndf);
    }
    // replicas of half A noise: distribution of the offset
    TRandom3 rnd(4242);
    double s1 = 0, s2 = 0; int n = 0;
    for (int k = 0; k < 2000; ++k) {
        std::vector<double> yy(d.size());
        for (size_t i = 0; i < d.size(); ++i) yy[i] = d[i] + o[i] + rnd.Gaus(0, std::sqrt(ovar[i]));
        LinResult r; double b = fitBR(yy, true, r);
        s1 += b; s2 += b * b; n++;
    }
    printf("parent-noise replicas: mean offset %.3e rms %.3e\n", s1 / n, std::sqrt(s2 / n - s1 * s1 / n / n));
}
}
void s5_blind_ana_xdev() { xdev::Run(); }
