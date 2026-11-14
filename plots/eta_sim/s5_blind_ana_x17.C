// ============================================================================
// s5_blind_ana_x17.C  --  blind analysis of the X17 pseudo-data (label blind_ana)
// ============================================================================
// Applies the FROZEN procedure tables/s5_blind_ana_procedure.txt (X1 primary,
// X2, X3 variants) to every toy of blind/pseudo/x17_toys.root, treated like
// measured data at L_data = 6080.5 nb^-1 (UNCONFIRMED). Inputs: the toys and the
// half-A templates only; nothing from blind/sealed.
//
//   X1: Poisson ML fit of M(ee) in 13.0-26.5 MeV/c2: BR x X17 line (template)
//       + a_D x direct5 (template) + a_O x smooth 'other' (double exponential fitted
//       to the half-A otherX17 + conv4 on 8-60 MeV/c2); over-dispersion scale
//       S = max(1, chi2/ndf); UL90 asymptotic CLs (q~_mu) with sigma_tot.
//   X2: as X1 with the raw half-A 'other' template + template variance.
//   X3: as X1 on 8-40 MeV/c2.
//   POST-FREEZE (declared in the procedure file): X4 = X1 with a_D constrained to 1 +- 0.10;
//   X5 = counting in the window with B from the half-A sumT12 template.
// Output: /hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/blind/results/x17_results.root
//         (tree "x17", one entry per toy id) and x17_results.txt
// Run: root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/s5_blind_ana_build", kTRUE)' plots/eta_sim/s5_blind_ana_x17.C+
// ============================================================================
#include <cmath>
#include <cstdio>
#include <vector>

#include "Math/ProbFuncMathCore.h"
#include "Math/QuantFuncMathCore.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TNamed.h"
#include "TTree.h"

#include "s5_blind_ana_fit.h"

using namespace bla;

namespace bx17 {

const char* const kToyFile = "/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/blind/pseudo/x17_toys.root";
const char* const kTplFile = "/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/blind/templates/templates.root";
const char* const kOutRoot = "/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/blind/results/x17_results.root";
const char* const kOutTxt  = "/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/blind/results/x17_results.txt";
const double kSc = 6080.5 / 221.05;   // L_data (UNCONFIRMED) / L_sim

// asymptotic CLs upper limit for a Gaussian estimator b +- s (q~_mu, Cowan et al. 2011)
double ULCLs(double b, double s, double cl = 0.90) {
    const double a = 1 - cl;
    const double pb = ROOT::Math::normal_cdf(b / s);
    return b - s * ROOT::Math::normal_quantile(a * pb, 1.0);
}
// one-sided profile-likelihood limit without CLs (q~_mu = z^2, z = 1.2816 for 90 %)
double ULPL(double b, double s, double cl = 0.90) {
    const double z = ROOT::Math::normal_quantile(cl, 1.0);
    if (b >= 0) return b + z * s;
    return b + std::sqrt(b * b + z * z * s * s);
}

struct XFit { LinResult r; double BR = 0, eStat = 0, S = 1, eTot = 0; std::vector<double> w; };

XFit FitX(const TH1D* h, int b1, int b2, const std::vector<double>& X, const std::vector<double>& D, const std::vector<double>& O,
          const std::vector<double>* vO, bool withX = true, double consD = -1) {
    LinProblem P;
    if (consD > 0 && withX) P.cons.push_back(Constraint{1, 1.0, consD});
    for (int b = b1; b <= b2; ++b) {
        P.y.push_back(h->GetBinContent(b));
        Vec f; if (withX) f.push_back(X[b]); f.push_back(D[b]); f.push_back(O[b]);
        P.F.push_back(f);
        if (vO) { Vec tv; if (withX) tv.push_back(0.); tv.push_back(0.); tv.push_back((*vO)[b]); P.T.push_back(tv); }
    }
    P.n = P.y.size(); P.np = withX ? 3 : 2;
    XFit x; x.r = LinFit(P);
    if (!withX) return x;
    x.BR = x.r.a[0]; x.eStat = std::sqrt(x.r.C[0][0]);
    x.S = std::max(1.0, x.r.chi2P / x.r.ndf);
    x.eTot = x.eStat * std::sqrt(x.S);
    Vec g(3, 0.); g[0] = 1; x.w = x.r.Weights(g);
    return x;
}

void Run() {
    // ---- templates (bins 1..200, index = bin number)
    TFile ft(kTplFile);
    TH1D* hD = (TH1D*)ft.Get("mee_X17_direct5");
    TH1D* hC = (TH1D*)ft.Get("mee_X17_conv4");
    TH1D* hO = (TH1D*)ft.Get("mee_X17_otherX17");
    TH1D* hX = (TH1D*)ft.Get("mee_X17_x17perBR");
    std::vector<double> X(201, 0.), D(201, 0.), Oraw(201, 0.), vO(201, 0.), Os(201, 0.), vX(201, 0.);
    TH1D* hOC = (TH1D*)hO->Clone("hOC"); hOC->Add(hC); hOC->Scale(kSc);
    for (int b = 1; b <= 200; ++b) {
        X[b] = kSc * hX->GetBinContent(b); vX[b] = std::pow(kSc * hX->GetBinError(b), 2);
        D[b] = kSc * hD->GetBinContent(b);
        Oraw[b] = hOC->GetBinContent(b); vO[b] = std::pow(hOC->GetBinError(b), 2);
    }
    TF1 sm("sm", "[0]*exp(-x/[1])+[2]*exp(-x/[3])", 8., 60.);
    sm.SetParameters(3000, 5, 300, 30);
    hOC->Fit(&sm, "QNR0");
    printf("smooth 'other': %.4g exp(-m/%.3f) + %.4g exp(-m/%.3f), chi2/ndf %.1f/%d\n", sm.GetParameter(0), sm.GetParameter(1), sm.GetParameter(2), sm.GetParameter(3), sm.GetChisquare(), sm.GetNDF());
    for (int b = 1; b <= 200; ++b) { double s = 0; for (int q = 0; q < 10; ++q) s += sm.Eval(hOC->GetBinLowEdge(b) + (q + 0.5) * 0.05); Os[b] = s / 10.; }
    double xw = 0, dw = 0, ow = 0, osw = 0;
    for (int b = 27; b <= 53; ++b) { xw += X[b]; dw += D[b]; ow += Oraw[b]; osw += Os[b]; }
    printf("window 13.0-26.5: X17 rows per unit BR %.4g; direct5 %.1f; other raw %.1f, smooth %.1f (L_data UNCONFIRMED)\n", xw, dw, ow, osw);

    // ---- toys
    TFile fin(kToyFile);
    TTree* tt = (TTree*)fin.Get("toys");
    Int_t id, index; tt->SetBranchAddress("id", &id); tt->SetBranchAddress("index", &index);
    TFile fo(kOutRoot, "RECREATE");
    TTree tr("x17", "blind analysis of the X17 toys (label blind_ana): one entry per toy");
    Int_t oId, oIndex, oNdf;
    Double_t oBR, oStat, oTot, oS, oSsb, oULcls, oULclsStat, oULpl, oAD, oAO, oChi2, oPar, oNwin;
    Double_t o2BR, o2Err, o2UL, o3BR, o3Stat, o3S, o3Tot, o3UL;
    Double_t o4BR, o4Stat, o4S, o4Tot, o4UL, o4AD, o5BR, o5Stat, o5Tot, o5UL;
    tr.Branch("id", &oId, "id/I"); tr.Branch("index", &oIndex, "index/I");
    tr.Branch("BR", &oBR, "BR/D"); tr.Branch("eBRstat", &oStat, "eBRstat/D"); tr.Branch("eBRtot", &oTot, "eBRtot/D");
    tr.Branch("S", &oS, "S/D"); tr.Branch("S_sideband", &oSsb, "S_sideband/D");
    tr.Branch("UL90_CLs", &oULcls, "UL90_CLs/D"); tr.Branch("UL90_CLs_stat", &oULclsStat, "UL90_CLs_stat/D"); tr.Branch("UL90_PL", &oULpl, "UL90_PL/D");
    tr.Branch("aD", &oAD, "aD/D"); tr.Branch("aO", &oAO, "aO/D"); tr.Branch("chi2", &oChi2, "chi2/D"); tr.Branch("ndf", &oNdf, "ndf/I");
    tr.Branch("eBRparent", &oPar, "eBRparent/D"); tr.Branch("Nwin", &oNwin, "Nwin/D");
    tr.Branch("X2_BR", &o2BR, "X2_BR/D"); tr.Branch("X2_eBR", &o2Err, "X2_eBR/D"); tr.Branch("X2_UL90_CLs", &o2UL, "X2_UL90_CLs/D");
    tr.Branch("X3_BR", &o3BR, "X3_BR/D"); tr.Branch("X3_eBRstat", &o3Stat, "X3_eBRstat/D"); tr.Branch("X3_S", &o3S, "X3_S/D");
    tr.Branch("X3_eBRtot", &o3Tot, "X3_eBRtot/D"); tr.Branch("X3_UL90_CLs", &o3UL, "X3_UL90_CLs/D");
    // post-freeze variants (declared in the procedure file)
    tr.Branch("X4_BR", &o4BR, "X4_BR/D"); tr.Branch("X4_eBRstat", &o4Stat, "X4_eBRstat/D"); tr.Branch("X4_S", &o4S, "X4_S/D");
    tr.Branch("X4_eBRtot", &o4Tot, "X4_eBRtot/D"); tr.Branch("X4_UL90_CLs", &o4UL, "X4_UL90_CLs/D"); tr.Branch("X4_aD", &o4AD, "X4_aD/D");
    tr.Branch("X5_BR", &o5BR, "X5_BR/D"); tr.Branch("X5_eBRstat", &o5Stat, "X5_eBRstat/D"); tr.Branch("X5_eBRtot", &o5Tot, "X5_eBRtot/D");
    tr.Branch("X5_UL90_CLs", &o5UL, "X5_UL90_CLs/D");
    double Bsum = 0, Bsum2 = 0;
    { TH1D* hs = (TH1D*)ft.Get("mee_X17_sumT12"); for (int b = 27; b <= 53; ++b) { Bsum += kSc * hs->GetBinContent(b); Bsum2 += std::pow(kSc * hs->GetBinError(b), 2); } }
    printf("X5: B (half-A sumT12 in window) = %.1f +- %.1f rows (L_data UNCONFIRMED)\n", Bsum, std::sqrt(Bsum2));
    FILE* ftxt = fopen(kOutTxt, "w");
    fprintf(ftxt, "# x17_results.txt -- blind analysis (label blind_ana), frozen procedure plots/eta_sim/tables/s5_blind_ana_procedure.txt\n");
    fprintf(ftxt, "# L_data = 6080.5 nb^-1 (UNCONFIRMED). X1 primary: BR +- eBRtot (= eBRstat x sqrt(S)), UL90 = asymptotic CLs (q~_mu) with eBRtot\n");
    fprintf(ftxt, "# id index BR eBRstat eBRtot S S_sideband UL90_CLs UL90_CLs_stat UL90_PL aD aO chi2 ndf | X2: BR eBR UL90 | X3: BR eBRstat S eBRtot UL90 | X4 (post-freeze): BR eBRstat S eBRtot UL90 aD | X5 (post-freeze): BR eBRstat eBRtot UL90\n");
    for (Long64_t i = 0; i < tt->GetEntries(); ++i) {
        tt->GetEntry(i);
        TH1D* h = (TH1D*)fin.Get(Form("toy_%d_mee", id));
        if (!h) { printf("missing toy %d\n", id); continue; }
        XFit x1 = FitX(h, 27, 53, X, D, Os, nullptr);
        XFit x2 = FitX(h, 27, 53, X, D, Oraw, &vO);
        XFit x3 = FitX(h, 17, 80, X, D, Os, nullptr);
        oNwin = h->Integral(27, 53);
        // side-band over-dispersion: D + O_s (free norms) on 26.5-60 MeV/c2
        XFit sb = FitX(h, 54, 120, X, D, Os, nullptr, false);
        oId = id; oIndex = index;
        oBR = x1.BR; oStat = x1.eStat; oTot = x1.eTot; oS = x1.S; oSsb = sb.r.chi2P / sb.r.ndf;
        oULcls = ULCLs(x1.BR, x1.eTot); oULclsStat = ULCLs(x1.BR, x1.eStat); oULpl = ULPL(x1.BR, x1.eTot);
        oAD = x1.r.a[1]; oAO = x1.r.a[2]; oChi2 = x1.r.chi2P; oNdf = x1.r.ndf;
        double vp = 0; for (int b = 27; b <= 53; ++b) vp += x1.w[b - 27] * x1.w[b - 27] * vO[b];
        oPar = std::sqrt(vp);
        o2BR = x2.BR; o2Err = x2.eStat; o2UL = ULCLs(x2.BR, x2.eStat);
        o3BR = x3.BR; o3Stat = x3.eStat; o3S = x3.S; o3Tot = x3.eTot; o3UL = ULCLs(x3.BR, x3.eTot);
        XFit x4 = FitX(h, 27, 53, X, D, Os, nullptr, true, 0.10);
        o4BR = x4.BR; o4Stat = x4.eStat; o4S = x4.S; o4Tot = x4.eTot; o4UL = ULCLs(x4.BR, x4.eTot); o4AD = x4.r.a[1];
        o5BR = (oNwin - Bsum) / xw; o5Stat = std::sqrt(oNwin) / xw; o5Tot = std::sqrt(oNwin + Bsum2) / xw; o5UL = ULCLs(o5BR, o5Tot);
        tr.Fill();
        fprintf(ftxt, "%d %d %.4e %.3e %.3e %.2f %.2f %.4e %.4e %.4e %.4f %.4f %.1f %d | %.4e %.3e %.4e | %.4e %.3e %.2f %.3e %.4e | %.4e %.3e %.2f %.3e %.4e %.3f | %.4e %.3e %.3e %.4e\n",
                oId, oIndex, oBR, oStat, oTot, oS, oSsb, oULcls, oULclsStat, oULpl, oAD, oAO, oChi2, oNdf, o2BR, o2Err, o2UL, o3BR, o3Stat, o3S, o3Tot, o3UL,
                o4BR, o4Stat, o4S, o4Tot, o4UL, o4AD, o5BR, o5Stat, o5Tot, o5UL);
        if (i < 3) {   // example fits
            fo.cd();
            TH1D hd(Form("ex%lld_data", i), Form("toy id %d", id), 200, 0., 100.);
            TH1D hm(Form("ex%lld_model", i), "X1 model", 200, 0., 100.), hb(Form("ex%lld_bkg", i), "X1 background", 200, 0., 100.), hx(Form("ex%lld_x17", i), "X1 X17", 200, 0., 100.);
            for (int b = 1; b <= 200; ++b) { hd.SetBinContent(b, h->GetBinContent(b)); hd.SetBinError(b, std::sqrt(h->GetBinContent(b))); }
            for (int b = 27; b <= 53; ++b) {
                hm.SetBinContent(b, x1.r.mu[b - 27]);
                hx.SetBinContent(b, x1.BR * X[b]);
                hb.SetBinContent(b, x1.r.mu[b - 27] - x1.BR * X[b]);
            }
            hd.Write(); hm.Write(); hb.Write(); hx.Write();
        }
        delete h;
    }
    fclose(ftxt);
    fo.cd();
    tr.Write();
    // templates used (L_data scale) for the figures
    TH1D hXs("tpl_X17_perBR", "X17 line per unit BR (L_data)", 200, 0., 100.), hDs("tpl_direct5", "direct5 (L_data)", 200, 0., 100.),
         hOs("tpl_other_smooth", "smooth other (L_data)", 200, 0., 100.), hOr("tpl_other_raw", "raw half-A other + conv (L_data)", 200, 0., 100.);
    for (int b = 1; b <= 200; ++b) { hXs.SetBinContent(b, X[b]); hDs.SetBinContent(b, D[b]); hOs.SetBinContent(b, Os[b]); hOr.SetBinContent(b, Oraw[b]); hOr.SetBinError(b, std::sqrt(vO[b])); }
    hXs.Write(); hDs.Write(); hOs.Write(); hOr.Write();
    TNamed info("info", "blind analysis of the X17 toys (label blind_ana, s5_blind_ana_x17.C), frozen procedure plots/eta_sim/tables/s5_blind_ana_procedure.txt; "
                        "X1 primary: BR +- eBRtot, UL90_CLs = asymptotic CLs with the one-sided profile-likelihood statistic q~_mu; L_data = 6080.5 nb^-1 UNCONFIRMED; no sealed file read");
    info.Write();
    fo.Close();
    printf("done: %lld toys\n", tt->GetEntries());
}
}
void s5_blind_ana_x17() { bx17::Run(); }
