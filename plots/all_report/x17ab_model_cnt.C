// x17ab_model_cnt.C -- COPY of x17ab_model.C (figure agent B, 2026-11-15): "counts" instead of "rows" in the annotation of
// Fig. x17ab_model_mee; writes ONLY docs/latex/figs_all/x17ab_model_mee_cnt.{pdf,png}; its tables, root file and the other two
// figures go to /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figB/out/x17ab/. Nothing else changed.
// Run: root -l -b -q -e 'gSystem->SetBuildDir("<scratch>/aclic_x17ab_cnt", kTRUE)' 'plots/all_report/x17ab_model_cnt.C+O'
// ============================================================================
// x17ab_model.C  --  WP9 X-a: the X17 model on the X2D layout (data scale) and its simulation-based projection
// ============================================================================
// SPEC report_all. NO DATA ROW (N) IS READ: the model uses the data CB (like-sign, not blinded) and the simulation
// templates of hists_v1.root (sim/S0__*); the "projection" data set is the model itself (simulation Asimov set):
//   background per M(ee) bin n_j = 1.35 x the simulated (G_D + G_C) of that bin (data/SMASH indication of the earlier
//   studies), peak N5, N4 and the eta-related terms at a_eta = 0.40 (indication 0.3-0.5), CB = data CB (k = 1).
// It is a PRE-DATA projection, superseded by the data-driven Asimov limit of X-c (side-band fit, D19).
// Steps:
//   1. Model inputs (x17_model.h BuildFromHists): X2D 35 x 30 (M4 0.5225-0.6100 x M(ee) 2-62 MeV/c2), classes E1/E2,
//      B = smoothed (exp(pol2)) simulated G_D + G_C per class, peak shape = C5 per class, D = C5 (V3), C = C4,
//      X = C6 (17 MeV/c2), R_X = eps_X/eps_5; C2+C3 and G_eta factorised and smoothed; B smoothed with exp(pol2), corrected in the fit
//      by (1 + (b + b'u) t + (c + c'u) t^2) per class (bkgMode 4).
//   2. Uncertainty of the conversion tail (M(ee) >= 12 MeV/c2) from the 100 simulation replicas of C4.
//   3. Projection: Asimov limit (CLs 90 %, asymptotic) with bands; the effect of every model choice (nuisance list of WP9)
//      on the Asimov fit: refit the nominal Asimov set with the changed model, report r_hat (bias) and the limit.
//   4. Projection versus m_X (10-40 MeV/c2): lines from x17_lines.root (morphed C6 line, R_X(m)).
// Output: docs/latex/tables_all/x17ab_model.txt, x17ab_nuis.txt, x17ab_proj_mass.txt; figures x17ab_model_mee, x17ab_model_m4,
//         x17ab_proj_mass (.pdf/.png); /hdd3/.../x17ab/x17ab_model.root
// Run: root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_x17ab", kTRUE)' 'plots/all_report/x17ab_model.C+O'
// ============================================================================
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TStyle.h"
#include "TSystem.h"

#include "defs.h"
#include "hists_io.h"
#include "x17_model.h"

namespace x17mod {
const char* const kOut = "/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figB/out/x17ab";   // _cnt copy: writes to scratch
const char* const kFig = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_all";
const char* const kTab = "/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figB/out/x17ab";   // _cnt copy: tables to scratch
const double kAproj = 0.40;     // a_eta of the projection
const double kBgScale = 1.35;   // data/SMASH background of the projection
const double kRscale = 1e-5 / x17m::kBRdir;

struct Proj {
    x17m::Inputs I;
    x17m::Config C;
    std::vector<double> asimov;
    std::vector<double> ptrue;
};
// conversion-tail relative uncertainty from the C4 replicas
double CtailSigma(TFile* f, const x17m::Inputs& I, double& nom, double& rms) {
    TH3* h = (TH3*)f->Get("sim/S0__nom/C4_X2D_rep");
    TH2* h0 = (TH2*)f->Get("sim/S0__nom/C4_X2D");
    if (!h || !h0) return 0.3;
    auto tailFrac = [&](std::function<double(int, int)> get) {
        double t = 0, a = 0;
        for (int j = 0; j <= I.G.ny + 1; ++j) for (int i = 1; i <= I.G.nx; ++i) {
            const double v = get(i, j);
            a += v;
            if (j >= 1 && j <= I.G.ny && I.Ctail[j - 1] > 0) t += v;
            if (j == I.G.ny + 1) t += v;   // M(ee) above 62 MeV/c2 belongs to the tail as well
        }
        return a > 0 ? t / a : 0.;
    };
    nom = tailFrac([&](int i, int j) { return h0->GetBinContent(i, j); });
    double s = 0, s2 = 0; int n = 0;
    for (int r = 1; r <= h->GetNbinsZ(); ++r) {
        const double v = tailFrac([&](int i, int j) { return h->GetBinContent(i, j, r); });
        s += v; s2 += v * v; ++n;
    }
    const double m = s / n;
    rms = std::sqrt(std::max(0., s2 / n - m * m));
    return nom > 0 ? rms / nom : 0.3;
}
// projection parameters: n_j from the simulation, peak at a_eta
std::vector<double> ProjParams(TFile* f, const x17m::Model& M, const x17m::Inputs& I) {
    std::vector<double> p(M.L.np, 0.);
    TH2* gd = (TH2*)f->Get("sim/S0__nom/GD_X2D");
    TH2* gc = (TH2*)f->Get("sim/S0__nom/Gconv_X2D");
    for (int j = 0; j < I.G.ny; ++j) {
        double s = 0;
        for (int i = 1; i <= I.G.nx; ++i) s += gd->GetBinContent(i, j + 1) + gc->GetBinContent(i, j + 1);
        p[j] = kBgScale * s * alldefs::kNbToRaw;
    }
    p[M.L.iN5] = kAproj * I.N5start; p[M.L.iN4] = kAproj * I.N4start;
    p[M.L.iK] = 1; p[M.L.iK + 1] = 1; p[M.L.iA] = kAproj; p[M.L.iG] = 1; p[M.L.iR] = 0;
    return p;
}
}  // namespace x17mod

void x17ab_model_cnt() {
    using namespace x17mod;
    gSystem->mkdir(kOut, kTRUE);
    gStyle->SetOptStat(0);
    TFile* f = TFile::Open(alldefs::kHistFile, "READ");
    // ---------------------------------------------------------------- 1. nominal inputs
    x17m::HistOpt o;
    o.aEta = 1.;
    x17m::Inputs I;
    if (!x17m::BuildFromHists(f, I, o, alldefs::kNbToRaw)) return;
    double tailNom = 0, tailRms = 0;
    const double sC = CtailSigma(f, I, tailNom, tailRms);
    x17m::Config C;
    C.kSig[0] = 0.5; C.kSig[1] = 0.05;
    C.a0 = kAproj; C.aSig = 0.10; C.gSig = 0.30;
    C.bkgMode = 4;   // B (1 + (b + b'u) t + (c + c'u) t^2), u = (M(ee) - 30)/30 MeV, per class (see the results note)
    C.cTailSig = sC;
    x17m::Model M(I, C);
    const x17m::Glob g0 = x17m::GlobOf(C);
    std::vector<double> ptrue = ProjParams(f, M, I);
    std::vector<double> asimov = M.Expected(ptrue);
    double sumD = 0, sumC = 0, sumX = 0, dWin = 0, cWin = 0, xWin = 0;
    for (int j = 0; j < I.G.ny; ++j) {
        sumD += I.D[j]; sumC += I.C[j]; sumX += I.X[j];
        if (I.G.ye[j] >= 0.0139 && I.G.ye[j + 1] <= 0.0261) { dWin += I.D[j]; cWin += I.C[j]; xWin += I.X[j]; }
    }
    double nAs = 0, nCB = 0;
    for (size_t b = 0; b < asimov.size(); ++b) { nAs += asimov[b]; nCB += I.CB[b]; }
    printf("x17ab_model: R_X = %.4f; N5(sim, a=1) %.1f, N4 %.1f rows; conversion tail fraction %.4f +- %.4f (rel %.3f)\n", I.RX, I.N5start, I.N4start, tailNom, tailRms, sC);
    printf("x17ab_model: fractions in 2-62: D %.4f C %.4f X %.4f; in 14-26: D %.4f C %.4f X %.4f; projection rows %.0f (CB %.0f)\n", sumD, sumC, sumX, dWin, cWin, xWin, nAs, nCB);
    const double s1 = ptrue[M.L.iN5] * I.RX / x17m::kBRdir;   // X17 rows per unit BR at a_eta = 0.4 in this model
    // ---------------------------------------------------------------- 2. projection limit
    x17m::Limit L0 = x17m::LimitAsym(M, asimov, g0, kRscale, true);
    printf("x17ab_model: projection: rhat %.3e +- %.3e; UL(CLs) r %.4e = BR %.3e; expected -2..+2: %.3e %.3e %.3e %.3e %.3e (BR); UL PL %.3e\n", L0.rhat,
           L0.rhatErr, L0.ul, L0.ul * x17m::kBRdir, L0.exp[0] * x17m::kBRdir, L0.exp[1] * x17m::kBRdir, L0.exp[2] * x17m::kBRdir, L0.exp[3] * x17m::kBRdir,
           L0.exp[4] * x17m::kBRdir, L0.ulPL * x17m::kBRdir);
    // ---------------------------------------------------------------- 3. model choices on the projection
    struct Var { std::string name, text; double rhat, rerr, ul, N5; bool ok; };
    std::vector<Var> vars;
    auto runVar = [&](const std::string& name, const std::string& text, const x17m::Inputs& Iv, const x17m::Config& Cv) {
        x17m::Model Mv(Iv, Cv);
        x17m::Limit Lv = x17m::LimitAsym(Mv, asimov, x17m::GlobOf(Cv), kRscale, false);
        Var v{name, text, Lv.rhat, Lv.rhatErr, Lv.ul, Lv.prof.hat.p.empty() ? 0. : Lv.prof.hat.p[Mv.L.iN5], Lv.ok};
        vars.push_back(v);
        printf("x17ab_model: variant %-12s rhat %+.3e (BR %+.2e) +- %.2e, UL BR %.3e, N5 %.1f (%s)\n", name.c_str(), v.rhat, v.rhat * x17m::kBRdir,
               v.rerr * x17m::kBRdir, v.ul * x17m::kBRdir, v.N5, text.c_str());
    };
    runVar("nominal", "nominal model", I, C);
    // X line: M(ee) scale and resolution
    for (int k = 0; k < 4; ++k) {
        const double sc[4] = {1.02, 0.98, 1., 1.}, rs[4] = {0., 0., 0.10, -0.10};
        const char* nm[4] = {"meeScale+2%", "meeScale-2%", "lineRes+10%", "lineRes-10%"};
        x17m::Inputs Iv = I;
        double RX;
        if (!x17m::LineFor(17., Iv.G, Iv.X, RX, 0, sc[k], rs[k])) continue;
        Iv.RX = I.RX;
        runVar(nm[k], "X17 line with M(ee) scale x " + std::to_string(sc[k]) + ", width x " + std::to_string(1 + rs[k]), Iv, C);
    }
    // line from the C5 slice instead of C6
    {
        x17m::Inputs Iv = I; double RX;
        if (x17m::LineFor(17., Iv.G, Iv.X, RX, 0, 1., 0., true)) { Iv.RX = I.RX; runVar("lineSlice", "X17 line from the reweighted C5 slice at 17 MeV/c2", Iv, C); }
    }
    // eps_X/eps_5 variants: R_X scaled
    {
        const char* nm[3] = {"epsDL", "epsC2", "epsNP"};
        const char* tx[3] = {"R_X: direct-like M(pipi)", "R_X: e+ decay angle 1+cos^2", "R_X: generator pion angle kept"};
        for (int v = 1; v <= 3; ++v) {
            std::vector<double> X; double RXv, RX0;
            if (!x17m::LineFor(17., I.G, X, RX0, 0) || !x17m::LineFor(17., I.G, X, RXv, v)) continue;
            x17m::Inputs Iv = I; Iv.RX = I.RX * RXv / RX0;
            runVar(nm[v - 1], tx[v - 1], Iv, C);
        }
    }
    // direct-decay shape (D and peak shape from the variant directories)
    for (const char* var : {"V2", "GS", "r089", "r110", "sUp", "sDn"}) {
        x17m::HistOpt ov = o; ov.var = var;
        x17m::Inputs Iv;
        if (!x17m::BuildFromHists(f, Iv, ov, alldefs::kNbToRaw)) continue;
        x17m::Inputs Iw = I;   // keep everything but the peak shape, D and R_X of the variant
        Iw.D = Iv.D; Iw.gP = Iv.gP; Iw.gX = Iv.gX; Iw.RX = Iv.RX;
        runVar(std::string("sim_") + var, std::string("peak/direct templates of S0__") + var, Iw, C);
    }
    // conversion template: 20 replicas
    std::vector<double> repR;
    for (int r = 0; r < 20; ++r) {
        x17m::HistOpt ov = o; ov.simRep = r;
        x17m::Inputs Iv;
        if (!x17m::BuildFromHists(f, Iv, ov, alldefs::kNbToRaw)) continue;
        x17m::Inputs Iw = I; Iw.C = Iv.C;
        x17m::Model Mv(Iw, C);
        x17m::Profile P = x17m::ProfileFits(Mv, asimov, g0);
        repR.push_back(P.rhat);
    }
    double rm = 0, rr = 0;
    for (double v : repR) rm += v;
    rm /= std::max<size_t>(1, repR.size());
    for (double v : repR) rr += (v - rm) * (v - rm);
    rr = std::sqrt(rr / std::max<size_t>(1, repR.size() - 1));
    printf("x17ab_model: C4 replicas (%zu): rhat mean %.3e, RMS %.3e (BR %.2e)\n", repR.size(), rm, rr, rr * x17m::kBRdir);
    // background model and CB
    {
        x17m::Config Cv = C; Cv.bkgMode = 2; runVar("bkgPol2", "background flat x (1 + b t + c t^2) per class", I, Cv);
        Cv = C; Cv.bkgMode = 0; runVar("bkgTilt", "background B x (1 + b t) per class (tilt only)", I, Cv);
        Cv = C; Cv.bkgMode = 3; runVar("bkgCurv", "background B x (1 + b t + c t^2) per class (no M(ee) dependence)", I, Cv);
        Cv = C; Cv.bkgMode = 1; runVar("bkgFixed", "background shape fixed (no tilt)", I, Cv);
        x17m::HistOpt ov = o; ov.nbcl = 3;
        x17m::Inputs Iv;
        if (x17m::BuildFromHists(f, Iv, ov, alldefs::kNbToRaw)) {
            // the projection truth uses two classes; refit with three
            runVar("bkg3cl", "three background classes (E1, 10-30, >= 30 MeV), bkgMode 4", Iv, C);
            x17m::Config C3 = C; C3.bkgMode = 3;
            runVar("bkg3clCurv", "three background classes, bkgMode 3", Iv, C3);
        }
        for (int ordv : {1, 3}) {
            ov = o; ov.smoothOrder = ordv;
            if (x17m::BuildFromHists(f, Iv, ov, alldefs::kNbToRaw)) runVar(ordv == 1 ? "Bexppol1" : "Bexppol3", ordv == 1 ? "background template smoothed with exp(pol1)" : "background template smoothed with exp(pol3)", Iv, C);
        }
        ov = o; ov.smoothB = false;
        if (x17m::BuildFromHists(f, Iv, ov, alldefs::kNbToRaw)) runVar("bkgRaw", "background shape raw (not smoothed)", Iv, C);
        ov = o; ov.kappaB = 1.3;
        if (x17m::BuildFromHists(f, Iv, ov, alldefs::kNbToRaw)) runVar("kappaB1.3", "B = G_D + 1.3 G_C", Iv, C);
        for (int cbm = 1; cbm <= 2; ++cbm) {
            ov = o; ov.cbMode = cbm;
            if (x17m::BuildFromHists(f, Iv, ov, alldefs::kNbToRaw)) runVar(cbm == 1 ? "CBfact" : "CBgeo", cbm == 1 ? "CB factorised (sum rule)" : "CB per-row geometric mean", Iv, C);
        }
        Cv = C; Cv.kSig[0] = -1.; runVar("kE1free", "k_E1 free", I, Cv);
        Cv = C; Cv.kSig[1] = 0.; runVar("kE2fixed", "k_E2 fixed at 1", I, Cv);
        Cv = C; Cv.cTailSig = 0.; runVar("noCtail", "no conversion-tail nuisance", I, Cv);
        Cv = C; Cv.gSig = 0.; Cv.aSig = 0.; runVar("etaFixed", "a_eta and g fixed", I, Cv);
        Cv = C; Cv.xSig = 0.10; runVar("xSig10", "X17 yield nuisance 10 %", I, Cv);
    }
    // N5 bias over the background-shape variants
    double n5min = 1e30, n5max = -1e30;
    for (const auto& v : vars)
        if ((v.name.rfind("bkg", 0) == 0 && v.name != "bkgRaw") || v.name == "nominal" || v.name.rfind("CB", 0) == 0 || v.name == "kappaB1.3") {
            n5min = std::min(n5min, v.N5); n5max = std::max(n5max, v.N5);
        }
    printf("x17ab_model: N5 over background/CB variants: %.1f - %.1f (truth %.1f): half spread %.2f %%\n", n5min, n5max, ptrue[M.L.iN5],
           50. * (n5max - n5min) / ptrue[M.L.iN5]);
    // ---------------------------------------------------------------- 4. projection versus mass
    std::vector<double> mgrid, ulm[5], rxm;
    for (double m = 10.; m <= 40.01; m += 1.0) {
        x17m::Inputs Iv = I;
        double RX;
        if (!x17m::LineFor(m, Iv.G, Iv.X, RX)) continue;
        Iv.RX = RX;
        x17m::Model Mv(Iv, C);
        x17m::Limit Lv = x17m::LimitAsym(Mv, asimov, g0, kRscale, true);
        mgrid.push_back(m); rxm.push_back(RX);
        for (int k = 0; k < 5; ++k) ulm[k].push_back(Lv.exp[k] * x17m::kBRdir);
        printf("x17ab_model: m %.1f R_X %.4f expected BR UL median %.3e (-1s %.3e, +1s %.3e)\n", m, RX, Lv.exp[2] * x17m::kBRdir, Lv.exp[1] * x17m::kBRdir, Lv.exp[3] * x17m::kBRdir);
    }

    // ---------------------------------------------------------------- tables
    FILE* t = fopen(Form("%s/x17ab_model.txt", kTab), "w");
    fprintf(t, "# x17ab_model.txt -- X17 model inputs on X2D (x17ab_model.C, x17_model.h). Simulation templates sim/S0__nom of hists_v1.root, data CB of\n");
    fprintf(t, "# data/S0/all (like-sign, not blinded). NO data N used. Units: expected RAW rows at L = 6079.2 nb^-1 x 63/64 (sim nb x %.4f).\n", alldefs::kNbToRaw);
    fprintf(t, "grid: M4 %.4f-%.4f GeV/c2 in %d bins x M(ee) %.3f-%.3f GeV/c2 in %d bins; classes E1 (< 10 MeV), E2\n", I.G.xe.front(), I.G.xe.back(), I.G.nx,
            I.G.ye.front(), I.G.ye.back(), I.G.ny);
    fprintf(t, "N5 (direct, C5 V3, all M(ee), M4 range) at a_eta = 1: %.1f rows; N4 (conversion C4): %.1f rows\n", I.N5start, I.N4start);
    fprintf(t, "R_X = eps_X/eps_5 (C6 per unit BR x 2.68e-4 / C5): %.4f\n", I.RX);
    fprintf(t, "X17 rows per unit BR_X BR_ee at a_eta = 1: s = N5 R_X / 2.68e-4 = %.4e (at a_eta = %.2f: %.4e); at BR 1e-5 and a_eta = 1: %.1f rows\n",
            I.N5start * I.RX / x17m::kBRdir, kAproj, s1, I.N5start * I.RX / x17m::kBRdir * 1e-5);
    fprintf(t, "fractions in M(ee) 2-62 MeV/c2: D %.4f, C %.4f, X %.4f; in 14-26 MeV/c2: D %.4f, C %.4f, X %.4f\n", sumD, sumC, sumX, dWin, cWin, xWin);
    fprintf(t, "conversion tail (M(ee) >= 12 MeV/c2) fraction of C4: %.4f +- %.4f (100 replicas) -> relative uncertainty sC = %.3f\n", tailNom, tailRms, sC);
    fprintf(t, "data CB (spread rule) in X2D: %.1f rows\n", nCB);
    fprintf(t, "\n# per M(ee) bin: D, C, X fractions; data CB (row sum); simulated G_D + G_C (a_eta-free background, raw rows at SMASH scale)\n");
    fprintf(t, "# Mee_lo  Mee_hi        D          C          X      CB_row   bkg_sim_row\n");
    for (int j = 0; j < I.G.ny; ++j) {
        double cb = 0; for (int i = 0; i < I.G.nx; ++i) cb += I.CB[I.G.idx(i, j)];
        fprintf(t, "%6.1f %6.1f  %9.5f  %9.5f  %9.5f  %8.1f  %10.1f\n", 1e3 * I.G.ye[j], 1e3 * I.G.ye[j + 1], I.D[j], I.C[j], I.X[j], cb, ptrue[j] / kBgScale);
    }
    fclose(t);
    t = fopen(Form("%s/x17ab_nuis.txt", kTab), "w");
    fprintf(t, "# x17ab_nuis.txt -- PRE-DATA PROJECTION (simulation Asimov set, no data N): background n_j = %.2f x simulated G_D + G_C per M(ee) bin,\n", kBgScale);
    fprintf(t, "# a_eta = %.2f, data CB. Nominal model: k_E1 1 +- 0.5, k_E2 1 +- 0.05, a_eta %.2f +- 0.10, g 1 +- 0.30, B (1 + (b + b'u) t + (c + c'u) t^2) per class, u = (M(ee) - 30 MeV)/30 MeV, conversion tail +- %.3f.\n", kAproj, kAproj, sC);
    fprintf(t, "# Each line: the nominal Asimov set refitted with the changed model: r_hat (the bias the change would cause if the nominal\n");
    fprintf(t, "# were true), its error, the asymptotic CLs 90 %% limit; BR = r x 2.68e-4.\n");
    fprintf(t, "projection (nominal): UL BR %.3e; expected bands -2s %.3e -1s %.3e median %.3e +1s %.3e +2s %.3e; profile-likelihood UL (no CLs) %.3e\n",
            L0.ul * x17m::kBRdir, L0.exp[0] * x17m::kBRdir, L0.exp[1] * x17m::kBRdir, L0.exp[2] * x17m::kBRdir, L0.exp[3] * x17m::kBRdir, L0.exp[4] * x17m::kBRdir,
            L0.ulPL * x17m::kBRdir);
    fprintf(t, "X17 rows per unit BR at a_eta %.2f: %.4e; UL in rows: %.1f\n", kAproj, s1, L0.ul * ptrue[M.L.iN5] * I.RX);
    fprintf(t, "%-12s %12s %10s %10s %10s %9s  %s\n", "variant", "rhat(BR)", "err(BR)", "UL(BR)", "UL/UL0", "N5", "description");
    for (const auto& v : vars)
        fprintf(t, "%-12s %+12.3e %10.3e %10.3e %10.3f %9.1f  %s%s\n", v.name.c_str(), v.rhat * x17m::kBRdir, v.rerr * x17m::kBRdir, v.ul * x17m::kBRdir,
                vars[0].ul > 0 ? v.ul / vars[0].ul : 0., v.N5, v.text.c_str(), v.ok ? "" : " [fit not converged]");
    fprintf(t, "C4 replicas (%zu): r_hat RMS = %.3e in BR (mean %+.3e)\n", repR.size(), rr * x17m::kBRdir, rm * x17m::kBRdir);
    fprintf(t, "N5 over the background and CB variants (bkgRaw excluded: MC-spike artefact): %.1f - %.1f (truth %.1f), half spread %.2f %% = N5-bias nuisance\n", n5min, n5max, ptrue[M.L.iN5],
            50. * (n5max - n5min) / ptrue[M.L.iN5]);
    fclose(t);
    t = fopen(Form("%s/x17ab_proj_mass.txt", kTab), "w");
    fprintf(t, "# x17ab_proj_mass.txt -- PRE-DATA projection of the expected 90 %% CLs limit on BR(eta->pi+pi-X17) BR(X17->ee) versus m_X (simulation Asimov,\n");
    fprintf(t, "# a_eta = %.2f, background 1.35 x SMASH, data CB; X17 line = morphed C6 line, R_X(m) from the reweighted C5 slices). Superseded by X-c/X-d.\n", kAproj);
    fprintf(t, "# m[MeV]   R_X    -2s        -1s        median     +1s        +2s\n");
    for (size_t k = 0; k < mgrid.size(); ++k)
        fprintf(t, "%6.1f %7.4f %10.3e %10.3e %10.3e %10.3e %10.3e\n", mgrid[k], rxm[k], ulm[0][k], ulm[1][k], ulm[2][k], ulm[3][k], ulm[4][k]);
    fclose(t);

    // ---------------------------------------------------------------- figures
    // (a) M(ee) shapes of D, C (with replica band), X(17), X at 12 and 30 MeV
    {
        TCanvas c("cm", "", 900, 700);
        gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.12); gPad->SetLogy();
        TH1D hD("hD", ";M(e^{+}e^{-}) [MeV/c^{2}];fraction of the component per 2 MeV/c^{2}", I.G.ny, 1e3 * I.G.ye.front(), 1e3 * I.G.ye.back());
        TH1D hC("hC", "", I.G.ny, 1e3 * I.G.ye.front(), 1e3 * I.G.ye.back());
        TH1D hX("hX", "", I.G.ny, 1e3 * I.G.ye.front(), 1e3 * I.G.ye.back());
        TH1D hX12("hX12", "", I.G.ny, 1e3 * I.G.ye.front(), 1e3 * I.G.ye.back());
        TH1D hX30("hX30", "", I.G.ny, 1e3 * I.G.ye.front(), 1e3 * I.G.ye.back());
        std::vector<double> X12, X30; double rx;
        x17m::LineFor(12., I.G, X12, rx); x17m::LineFor(30., I.G, X30, rx);
        for (int j = 0; j < I.G.ny; ++j) {
            hD.SetBinContent(j + 1, I.D[j]); hC.SetBinContent(j + 1, I.C[j]); hX.SetBinContent(j + 1, I.X[j]);
            hX12.SetBinContent(j + 1, X12[j]); hX30.SetBinContent(j + 1, X30[j]);
        }
        hD.SetLineColor(kBlue + 1); hD.SetLineWidth(3);
        hC.SetLineColor(kOrange + 7); hC.SetLineWidth(3); hC.SetLineStyle(2);
        hX.SetLineColor(kRed + 1); hX.SetLineWidth(3);
        hX12.SetLineColor(kMagenta + 1); hX12.SetLineWidth(2); hX12.SetLineStyle(3);
        hX30.SetLineColor(kGreen + 2); hX30.SetLineWidth(2); hX30.SetLineStyle(3);
        hD.SetMinimum(1e-5); hD.SetMaximum(1.);
        hD.GetXaxis()->SetTitleSize(0.05); hD.GetYaxis()->SetTitleSize(0.045);
        hD.Draw("HIST"); hC.Draw("HIST SAME"); hX.Draw("HIST SAME"); hX12.Draw("HIST SAME"); hX30.Draw("HIST SAME");
        TLegend lg(0.45, 0.66, 0.89, 0.89); lg.SetBorderSize(0); lg.SetTextSize(0.033);
        lg.AddEntry(&hD, "D: direct #eta#rightarrow#pi^{+}#pi^{-}e^{+}e^{-} (C5, V3)", "l");
        lg.AddEntry(&hC, Form("C: conversion (C4); tail #pm%.0f%%", 100 * sC), "l");
        lg.AddEntry(&hX, Form("X: X17, 17 MeV/c^{2} (C6), R_{X} = %.3f", I.RX), "l");
        lg.AddEntry(&hX12, "X at 12 MeV/c^{2} (morphed)", "l");
        lg.AddEntry(&hX30, "X at 30 MeV/c^{2} (morphed)", "l");
        lg.Draw();
        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.033);
        tx.DrawLatex(0.15, 0.17, "S0, M4' 0.5225-0.610 GeV/c^{2}; fractions of all counts (flows included)");
        c.SaveAs(Form("%s/x17ab_model_mee_cnt.pdf", kFig)); c.SaveAs(Form("%s/x17ab_model_mee_cnt.png", kFig));
    }
    // (b) M4 shapes: peak E1/E2, X17, background E1/E2 (smoothed), CB per class
    {
        TCanvas c("c4", "", 900, 700);
        gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.12);
        const int nx = I.G.nx;
        TH1D* h[7];
        const char* nm[7] = {"gP_E1", "gP_E2", "gX", "B_E1", "B_E2", "CB_E1", "CB_E2"};
        for (int k = 0; k < 7; ++k) { h[k] = new TH1D(nm[k], ";M'(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}];fraction per 2.5 MeV/c^{2}", nx, I.G.xe.data()); h[k]->SetDirectory(nullptr); }
        std::vector<double> cb[2] = {std::vector<double>(nx, 0.), std::vector<double>(nx, 0.)};
        for (int j = 0; j < I.G.ny; ++j) for (int i = 0; i < nx; ++i) cb[I.kcl[j]][i] += I.CB[I.G.idx(i, j)];
        for (int e = 0; e < 2; ++e) x17m::Normalise(cb[e]);
        for (int i = 0; i < nx; ++i) {
            h[0]->SetBinContent(i + 1, I.gP[0][i]); h[1]->SetBinContent(i + 1, I.gP[1][i]); h[2]->SetBinContent(i + 1, I.gX[i]);
            h[3]->SetBinContent(i + 1, I.B[0][i]); h[4]->SetBinContent(i + 1, I.B[1][i]);
            h[5]->SetBinContent(i + 1, cb[0][i]); h[6]->SetBinContent(i + 1, cb[1][i]);
        }
        const int col[7] = {kBlue + 1, kAzure + 7, kRed + 1, kGray + 2, kBlack, kOrange + 7, kRed - 7};
        const int sty[7] = {1, 2, 1, 2, 1, 2, 1};
        for (int k = 0; k < 7; ++k) { h[k]->SetLineColor(col[k]); h[k]->SetLineStyle(sty[k]); h[k]->SetLineWidth(3); }
        h[0]->SetMaximum(0.2); h[0]->SetMinimum(0);
        h[0]->GetXaxis()->SetTitleSize(0.05); h[0]->GetYaxis()->SetTitleSize(0.045); h[0]->GetXaxis()->SetNdivisions(505);
        h[0]->Draw("HIST");
        for (int k = 1; k < 7; ++k) h[k]->Draw("HIST SAME");
        TLegend lg(0.55, 0.55, 0.89, 0.89); lg.SetBorderSize(0); lg.SetTextSize(0.032);
        lg.AddEntry(h[0], "peak shape E1 (C5)", "l"); lg.AddEntry(h[1], "peak shape E2 (C5)", "l"); lg.AddEntry(h[2], "X17 (C6)", "l");
        lg.AddEntry(h[3], "background E1 (G_{D}+G_{C}, smoothed)", "l"); lg.AddEntry(h[4], "background E2", "l");
        lg.AddEntry(h[5], "data CB E1", "l"); lg.AddEntry(h[6], "data CB E2", "l");
        lg.Draw();
        c.SaveAs(Form("%s/x17ab_model_m4.pdf", kOut)); c.SaveAs(Form("%s/x17ab_model_m4.png", kOut));   // _cnt copy: scratch only
        for (int k = 0; k < 7; ++k) delete h[k];
    }
    // (c) projection vs mass
    if (!mgrid.empty()) {
        TCanvas c("cp", "", 900, 700);
        gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.12); gPad->SetLogy();
        const int n = (int)mgrid.size();
        TGraphAsymmErrors b2(n), b1(n);
        TGraph med(n);
        for (int k = 0; k < n; ++k) {
            b2.SetPoint(k, mgrid[k], ulm[2][k]); b2.SetPointError(k, 0.5, 0.5, ulm[2][k] - ulm[0][k], ulm[4][k] - ulm[2][k]);
            b1.SetPoint(k, mgrid[k], ulm[2][k]); b1.SetPointError(k, 0.5, 0.5, ulm[2][k] - ulm[1][k], ulm[3][k] - ulm[2][k]);
            med.SetPoint(k, mgrid[k], ulm[2][k]);
        }
        b2.SetFillColor(kOrange); b1.SetFillColor(kGreen + 1);
        b2.SetTitle(";m_{X} [MeV/c^{2}];90% CL_{s} on BR(#eta#rightarrow#pi^{+}#pi^{-}X17) BR(X17#rightarrowe^{+}e^{-})");
        b2.GetXaxis()->SetTitleSize(0.045); b2.GetYaxis()->SetTitleSize(0.04);
        b2.SetMinimum(3e-6); b2.SetMaximum(5e-4);
        b2.Draw("A3"); b1.Draw("3 SAME");
        med.SetLineStyle(2); med.SetLineWidth(2); med.Draw("L SAME");
        TLegend lg(0.42, 0.74, 0.89, 0.89); lg.SetBorderSize(0); lg.SetTextSize(0.033); lg.SetFillStyle(0);
        lg.AddEntry(&med, "expected (simulation Asimov, pre-data)", "l"); lg.AddEntry(&b1, "#pm 1#sigma", "f"); lg.AddEntry(&b2, "#pm 2#sigma", "f");
        lg.Draw();
        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.03);
        tx.DrawLatex(0.16, 0.17, Form("S0, X2D; a_{#eta} = %.2f, background 1.35 #times SMASH, data CB; L = 6079.2 nb^{-1}", kAproj));
        c.SaveAs(Form("%s/x17ab_proj_mass.pdf", kOut)); c.SaveAs(Form("%s/x17ab_proj_mass.png", kOut));   // _cnt copy: scratch only
    }
    // file with the inputs
    TFile* fo = TFile::Open(Form("%s/x17ab_model.root", kOut), "RECREATE");
    TH2D ha("asimov_proj", "projection Asimov set (simulation, r = 0)", I.G.nx, I.G.xe.data(), I.G.ny, I.G.ye.data());
    for (int j = 0; j < I.G.ny; ++j) for (int i = 0; i < I.G.nx; ++i) ha.SetBinContent(i + 1, j + 1, asimov[I.G.idx(i, j)]);
    ha.Write();
    fo->Close();
    f->Close();
    printf("x17ab_model: done\n");
}
