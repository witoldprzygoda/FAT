// ============================================================================
// x17p26_limit.C  --  X17 limit with PDG 2026 (analyst ruling D25, 15 November 2026). Prefix x17p26. Written 2026-11-15.
// ============================================================================
// Question: where does BR_dir = BR(eta -> pi+pi-e+e-) enter the frozen X17 procedure, and what changes with PDG 2026
// (2.76 +- 0.16 e-4, S = 1.5) instead of PDG 2020 (2.68 +- 0.11 e-4)?
//   - Central value: NOT in the likelihood. The POI is r = BR_X BR_ee / BR_dir, X17 rows = r N5 R_X with N5 fitted;
//     R_X = eps_X / eps_5 = C6 x 2.68e-4 / C5 in which 2.68e-4 only undoes the simulation weight of C5 (2 BR per forced
//     decay), so R_X is independent of the PDG edition. The BR limit is r_UL x BR_dir: a pure factor 2.76/2.68.
//   - Relative uncertainty: INSIDE the likelihood, as one term of the frozen X17 yield nuisance of normalisation (1):
//     sX1 = 0.246705 contains BR_dir 0.11/2.68 = 4.10 % (x17c_freeze.C). With PDG 2026 that term is 0.16/2.76 = 5.80 %.
// This macro reruns the FROZEN code (x17c_model.h, x17c_lib.h, x17d_lib.h, unchanged, only included) with ONLY that
// term changed: sX1' = sqrt(sX1^2 - (0.11/2.68)^2 + (0.16/2.76)^2). Everything else (data, configuration, qScale, grid,
// side-band fit, Asimov set) is the frozen one. Normalisations (2) and (3) contain no BR_dir and are not rerun.
// Normalisations (2) and (3) contain no BR_dir; their fixed eta scales were fitted with PDG 2020 BRs (a_eta: C2/C3 bump,
// x 1/1.00436 with PDG 2026; a_P: SMASH peak C4 + C5, x 1/1.0224; verifier v8). They are recomputed here with xA rescaled
// (this only checks the exact linearity UL ~ 1/xA). The 1 + cos^2 variant of (1) is rerun like (1).
// Output (new files only): /hdd3/przygoda/HADES/FAT_out/all_report/x17p26/x17p26_limit.log, x17p26_scan.txt, x17p26_state.root
//   (tree "obs26": PDG 2026 BR units: ul[3], rhat[3], rerr[3] (norm 1 rerun; 2, 3 with rescaled eta scales), ulC2, frozen-Asimov
//    bands e1[5] (norm 1 rerun), e2[5], e3[5]; and the frozen-procedure values in PDG 2020 units for the check)
// Run: root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_x17p26", kTRUE)' \
//        'plots/all_report/x17p26_limit.C+O'
// ============================================================================
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <string>
#include <vector>

#include "TFile.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TTree.h"

#include "x17d_lib.h"

namespace x17p26 {
const char* const kOut = "/hdd3/przygoda/HADES/FAT_out/all_report/x17p26";
const double kBR20 = 2.68e-4, kE20 = 0.11e-4;   // PDG 2020 (frozen)
const double kBR26 = 2.76e-4, kE26 = 0.16e-4;   // PDG 2026 (F. Takahashi et al., IJMPA 41 (2026) 2630011), S = 1.5
FILE* gLog = nullptr;
void P(const char* fmt, ...) __attribute__((format(printf, 1, 2)));
void P(const char* fmt, ...) {
    va_list a;
    va_start(a, fmt); vprintf(fmt, a); va_end(a);
    if (gLog) { va_start(a, fmt); vfprintf(gLog, fmt, a); va_end(a); fflush(gLog); }
}
}  // namespace x17p26

void x17p26_limit() {
    using namespace x17d;
    using x17p26::P;
    TStopwatch sw;
    gSystem->mkdir(x17p26::kOut, kTRUE);
    x17p26::gLog = fopen(Form("%s/x17p26_limit.log", x17p26::kOut), "w");
    Ctx c;
    if (!Init(c)) { P("x17p26: init failed\n"); return; }
    x17c::Inputs& I = c.I;
    const Frozen z = c.z;
    Frozen z26 = z;
    const double d20 = x17p26::kE20 / x17p26::kBR20, d26 = x17p26::kE26 / x17p26::kBR26;
    z26.xSig1 = std::sqrt(z.xSig1 * z.xSig1 - d20 * d20 + d26 * d26);
    const double fac = x17p26::kBR26 / x17p26::kBR20;
    P("x17p26: frozen sX1 %.6f (BR_dir part %.4f); PDG 2026 sX1' %.6f (BR_dir part %.4f); BR_dir factor 2.76/2.68 = %.6f; qScale %.5f\n", z.xSig1, d20,
      z26.xSig1, d26, fac, z.qScale);
    P("x17p26: code constant kBR (simulation-weight convention, used in R_X and in the BR printout) = %.3e; R_X(17) = %.4f (frozen 1.3455)\n", kBR, I.RX);
    const std::vector<double>& n = c.D.n;
    // frozen background-only Asimov set (side-band fit with the masked rows predicted), as in x17d_unblind.C step 0
    const std::vector<double> simRow = SimRow(c.f, I), simRel = SimRowRelErr(c.f, I);
    x17c::Config C1 = MakeConfig(z, 1), C1n = MakeConfig(z26, 1);
    x17c::Config Csb = C1; Csb.use = c.D.useSB;
    SideBand SB = FitSideBands(I, Csb, c.D.nSBonly, simRow, z.fillOrder, &simRel);
    x17c::Model M1(I, C1);
    const std::vector<double> asimovFrozen = M1.Expected(SB.pFill);
    {
        Lim lFz = LimitBR(I, C1, asimovFrozen, true);
        P("x17p26: reproduction: B_pred(X17c) %.1f (frozen %.1f); frozen Asimov median 17 MeV/c2 %.4e (frozen %.4e)\n", WindowX17c(I.G, asimovFrozen), kBpredX17c,
          lFz.exp[2], kFrozenExpMed17);
    }
    FILE* tb = fopen(Form("%s/x17p26_scan.txt", x17p26::kOut), "w");
    fprintf(tb, "# x17p26_scan.txt -- X17, normalisation (1), frozen procedure rerun with ONLY the BR_dir part of the X17 yield nuisance changed\n");
    fprintf(tb, "# (PDG 2020 0.11/2.68 = 4.10 %% -> PDG 2026 0.16/2.76 = 5.80 %%; sX1 %.6f -> %.6f). r = BR_X BR_ee / BR_dir (dimensionless).\n", z.xSig1, z26.xSig1);
    fprintf(tb, "# obs = observed asymptotic 90 %% CLs UL on r (all 1050 cells); expA = Asimov median and bands of the FROZEN background-only set (side bands).\n");
    fprintf(tb, "# BR20 = r x 2.68e-4 (frozen), BR26a = r(frozen sX) x 2.76e-4 (pure conversion), BR26b = r(PDG 2026 sX) x 2.76e-4 (rerun).\n");
    fprintf(tb, "# m      rUL_frozen  rUL_sX26   ratio  | BR20       BR26a      BR26b     | expA med r: frozen  sX26    | expA -2s -1s med +1s +2s (BR26b)\n");
    const double f2 = 1.00436, f3 = 1.0224;   // PDG 2026 rescaling of the WP5 eta scales (v8): a_eta / 1.00436, a_P / 1.0224
    x17c::Config C2n = MakeConfig(z, 2), C3n = MakeConfig(z, 3), C2o = MakeConfig(z, 2), C3o = MakeConfig(z, 3);
    C2n.xA = z.a0 / f2; C3n.xA = z.aP / f3;
    TFile* fo = TFile::Open(Form("%s/x17p26_state.root", x17p26::kOut), "RECREATE");
    TTree tr("obs26", "X17 limits in PDG 2026 BR units (x17p26_limit.C)");
    double m26, ul26[3], rh26[3], re26[3], c226, e126[5], e226[5], e326[5], ul20[3], c220, e120[5], rUL20, rUL26;
    tr.Branch("m", &m26); tr.Branch("ul", ul26, "ul[3]/D"); tr.Branch("rhat", rh26, "rhat[3]/D"); tr.Branch("rerr", re26, "rerr[3]/D");
    tr.Branch("ulC2", &c226); tr.Branch("e1", e126, "e1[5]/D"); tr.Branch("e2", e226, "e2[5]/D"); tr.Branch("e3", e326, "e3[5]/D");
    tr.Branch("ul20", ul20, "ul20[3]/D"); tr.Branch("ulC2_20", &c220); tr.Branch("e1_20", e120, "e1_20[5]/D"); tr.Branch("rUL20", &rUL20); tr.Branch("rUL26", &rUL26);
    for (double mm : MassGrid()) {
        x17c::Inputs Im = I;
        if (!SetLine(Im, mm)) continue;
        Lim a = LimitBR(Im, C1, n, false), b = LimitBR(Im, C1n, n, false);
        Lim ea = LimitBR(Im, C1, asimovFrozen, true), eb = LimitBR(Im, C1n, asimovFrozen, true);
        const double ra = a.ul / kBR, rb = b.ul / kBR;
        {
            const double g26 = x17p26::kBR26 / kBR;   // r -> BR with PDG 2026 (normalisation 1 only)
            Lim l2o = LimitBR(Im, C2o, n, false), l3o = LimitBR(Im, C3o, n, false);
            Lim l2 = LimitBR(Im, C2n, n, false), l3 = LimitBR(Im, C3n, n, false);
            Lim e2 = LimitBR(Im, C2n, asimovFrozen, true), e3 = LimitBR(Im, C3n, asimovFrozen, true);
            x17c::Inputs Ic = I;
            SetLine(Ic, mm, 2);
            Lim lc = LimitBR(Ic, C1, n, false), lcn = LimitBR(Ic, C1n, n, false);
            m26 = mm;
            ul26[0] = b.ul * g26; ul26[1] = l2.ul; ul26[2] = l3.ul;
            rh26[0] = b.rhat * g26; rh26[1] = l2.rhat; rh26[2] = l3.rhat;
            re26[0] = b.rerr * g26; re26[1] = l2.rerr; re26[2] = l3.rerr;
            c226 = lcn.ul * g26; c220 = lc.ul;
            for (int k = 0; k < 5; ++k) { e126[k] = eb.exp[k] * g26; e226[k] = e2.exp[k]; e326[k] = e3.exp[k]; e120[k] = ea.exp[k]; }
            ul20[0] = a.ul; ul20[1] = l2o.ul; ul20[2] = l3o.ul;
            rUL20 = ra; rUL26 = rb;
            tr.Fill();
            P("x17p26: m %5.2f  (2) %.4e -> %.4e (x %.5f)  (3) %.4e -> %.4e (x %.5f)  cos2 %.4e -> %.4e (x %.5f incl. 2.76/2.68)  r_hat(1) %+.4e -> %+.4e +- %.4e (BR26)\n",
              mm, l2o.ul, l2.ul, l2.ul / l2o.ul, l3o.ul, l3.ul, l3.ul / l3o.ul, lc.ul, c226, c226 / lc.ul, a.rhat, rh26[0], re26[0]);
        }
        P("x17p26: m %5.2f  r_UL frozen %.5f  sX26 %.5f  (x %.5f) | BR %.4e -> %.4e (conv) / %.4e (rerun) | Asimov med r %.5f / %.5f (x %.5f) | ok %d%d%d%d\n", mm, ra, rb,
          rb / ra, a.ul, ra * x17p26::kBR26, rb * x17p26::kBR26, ea.exp[2] / kBR, eb.exp[2] / kBR, eb.exp[2] / ea.exp[2], a.ok, b.ok, ea.ok, eb.ok);
        fprintf(tb, "%6.2f %10.5f %10.5f %7.5f | %10.4e %10.4e %10.4e | %9.5f %9.5f | %10.4e %10.4e %10.4e %10.4e %10.4e\n", mm, ra, rb, rb / ra, a.ul,
                ra * x17p26::kBR26, rb * x17p26::kBR26, ea.exp[2] / kBR, eb.exp[2] / kBR, eb.exp[0] / kBR * x17p26::kBR26, eb.exp[1] / kBR * x17p26::kBR26,
                eb.exp[2] / kBR * x17p26::kBR26, eb.exp[3] / kBR * x17p26::kBR26, eb.exp[4] / kBR * x17p26::kBR26);
    }
    fclose(tb);
    fo->cd(); tr.Write(); fo->Close();
    c.f->Close();
    P("x17p26_limit: done, %.0f s\n", sw.RealTime());
    fclose(x17p26::gLog);
}
