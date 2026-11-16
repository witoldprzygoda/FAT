// ============================================================================
// x17d_variants.C  --  WP9 X-d (prefix x17d): INFORMATION ONLY -- the observed limit at 17 MeV/c2 under the model variants
//                      that X-c used to size the frozen nuisances (x17c_freeze.C, runVar list without the masked-row fill
//                      variants, which have no meaning once all rows are fitted).
// ============================================================================
// The frozen limit (x17d_unblind.C) already contains these choices through the N5-bias and the X17 yield nuisance sX;
// this table is NOT a second limit. It is input for the UL_X17 column of Table S1 (WP7) and shows the size of each choice
// on the observed data. Every variant: its own inputs (x17cl::Build options as in x17c_freeze.C), the frozen configuration
// (MakeConfig, normalisation 1) with the variant's bkgMode / k_E1 setting, all 1050 cells.
// Output: docs/latex/tables_all/x17d_variants.txt
// Run: root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_x17d", kTRUE)' \
//        'plots/all_report/x17d_variants.C+O'
// ============================================================================
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "TFile.h"

#include "x17d_lib.h"

void x17d_variants(double mX = 17.0) {
    using namespace x17d;
    Ctx c;
    if (!Init(c)) { printf("x17d_variants: init failed\n"); return; }
    double rxHistNom = 0;
    { x17c::Inputs It; Build(c.f, It, c.w, BuildOpt(), c.z.a0, &rxHistNom); }
    FILE* tb = fopen(Form("%s/x17d_variants.txt", kTabD), "w");
    fprintf(tb, "# x17d_variants.txt -- INFORMATION ONLY (x17d_variants.C): observed 90 %% CLs UL at m_X = %.2f MeV/c2, normalisation (1), under the model\n", mX);
    fprintf(tb, "# variants that X-c used to size the frozen nuisances (N5-bias, sX). The frozen limit already contains them; this is not a second limit.\n");
    fprintf(tb, "# Data S0, isBest >= 0, all 1050 X2D cells, L = 6079.2 nb^-1. Pearson chi2/ndf of the r = 0 fit; N5 of the r = 0 fit.\n");
    fprintf(tb, "%-10s %10s %9s %11s %8s %11s %10s  %s\n", "variant", "chi2/ndf", "N5", "UL(1)", "UL/nom", "r_hat(BR)", "err", "description");
    double ulNom = 0;
    auto run = [&](const std::string& name, const std::string& text, const BuildOpt& b, int bkgMode, double kE1sig) {
        x17c::Inputs Iv;
        double rxh = 0;
        if (!Build(c.f, Iv, c.w, b, c.z.a0, &rxh)) return;
        if (!SetLine(Iv, mX)) return;
        Iv.RX *= rxh / rxHistNom;   // as x17c_freeze.C: R_X change of a simulation variant carried onto the line
        x17c::Config Cv = MakeConfig(c.z, 1);
        if (bkgMode >= 0) Cv.bkgMode = bkgMode;
        if (kE1sig != 0) Cv.kSig[0] = kE1sig;
        Lim l = LimitBR(Iv, Cv, c.D.n, false);
        x17c::Model Mv(Iv, Cv);
        const x17c::Profile P = x17c::ProfileFits(Mv, c.D.n, x17c::GlobOf(Cv));
        const Gof G = GoodnessOfFit(Mv, P.zero, c.D.n);
        if (name == "nominal") ulNom = l.ul;
        fprintf(tb, "%-10s %10.4f %9.1f %11.4e %8.3f %+11.3e %10.3e  %s\n", name.c_str(), G.chi2 / (G.ncell - G.nfree), P.zero.p[Mv.L.iN5], l.ul, ulNom > 0 ? l.ul / ulNom : 1.,
                l.rhat, l.rerr, text.c_str());
        printf("x17d_variants: %-10s UL %.4e (x %.3f) r_hat %+.3e chi2/ndf %.4f N5 %.1f\n", name.c_str(), l.ul, ulNom > 0 ? l.ul / ulNom : 1., l.rhat, G.chi2 / (G.ncell - G.nfree),
               P.zero.p[Mv.L.iN5]);
    };
    BuildOpt b;
    run("nominal", "frozen nominal (WP5 background shapes, bkgMode 4, 2 classes, spread-rule CB)", b, -1, 0);
    b = BuildOpt(); b.bshape = 1; run("Bsim2", "background shape: simulated G_D + G_C, exp(pol2)", b, -1, 0);
    b = BuildOpt(); b.bshape = 2; run("Bsim1", "background shape: simulated, exp(pol1)", b, -1, 0);
    b = BuildOpt(); b.bshape = 3; run("Bsim3", "background shape: simulated, exp(pol3)", b, -1, 0);
    b = BuildOpt(); b.bshape = 4; run("BkappaB", "background shape: simulated G_D + kappa_B G_C, exp(pol2)", b, -1, 0);
    b = BuildOpt(); run("flatPol2", "background flat x (1 + b t + c t^2) per class", b, 2, 0);
    b = BuildOpt(); b.bshape = 1; b.nbcl = 3; run("3class", "three background classes, simulated exp(pol2)", b, -1, 0);
    b = BuildOpt(); run("tiltOnly", "WP5 shapes x (1 + b t) per class", b, 0, 0);
    b = BuildOpt(); run("noMeeDep", "WP5 shapes x (1 + b t + c t^2), no M(ee) dependence", b, 3, 0);
    b = BuildOpt(); b.cbMode = 1; run("CBfact", "CB factorised (sum rule)", b, -1, 0);
    b = BuildOpt(); b.cbMode = 2; run("CBgeo", "CB per-row geometric mean", b, -1, 0);
    b = BuildOpt(); run("kE1free", "k_E1 free", b, -1, -1);
    for (const char* var : {"V2", "GS", "r089", "r110", "sUp", "sDn"}) {
        b = BuildOpt(); b.var = var;
        run(std::string("sim_") + var, std::string("peak/direct templates of S0__") + var, b, -1, 0);
    }
    fclose(tb);
    c.f->Close();
}
