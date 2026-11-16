// ============================================================================
// x17c_toycls.C  --  WP9 X-c (prefix x17c): toy-based CLs against the asymptotic formula on the background-only Asimov
//                    set of the data side-band fit (m_X given, normalisation (1)), one value of r per call.
// x17c::ToyCLs: pseudo-data from the conditional fits at r and at 0 (over-dispersed by qScale, constraint centres drawn),
// q~_r of every toy; CLs = CL_{s+b}/CL_b. x17c_summary-like collection: x17c_toycls_sum.C is not needed, the text lines of
// the logs are collected by x17c_frozen.py (one line per r: "x17c_toycls: m r CLs_toy err CLs_asym").
// Run: root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_x17c", kTRUE)' \
//        'plots/all_report/x17c_toycls.C+O(17.0, 0.10, 600)'
// ============================================================================
#include <cmath>
#include <cstdio>
#include <vector>

#include "TFile.h"

#include "x17c_lib.h"

void x17c_toycls(double mX = 17.0, double r = 0.10, int ntoys = 600, unsigned long long seed = 777ULL) {
    using namespace x17cl;
    TFile* f = TFile::Open(alldefs::kHistFile, "READ");
    Wp5 w = ReadWp5();
    Frozen z;
    if (!w.ok || !ReadConfig(z)) return;
    x17c::Inputs I;
    if (!Build(f, I, w, BuildOpt(), z.a0)) return;
    if (!SetLine(I, mX)) return;
    std::vector<double> n;
    std::vector<char> use;
    if (!LoadSideBands(f, n, use)) return;
    const std::vector<double> simRow = SimRow(f, I), simRel = SimRowRelErr(f, I);
    x17c::Config C1 = MakeConfig(z, 1);
    x17c::Config Csb = C1; Csb.use = use;
    SideBand SB = FitSideBands(I, Csb, n, simRow, z.fillOrder, &simRel);
    x17c::Model M(I, C1);
    const std::vector<double> asimov = M.Expected(SB.pFill);
    const x17c::Glob g = x17c::GlobOf(C1);
    const x17c::Profile P = x17c::ProfileFits(M, asimov, g);
    // asymptotic CLs at r for the same data set
    x17c::Profile PA;
    PA.zero = M.Fit(asimov, 0., false, g, &P.zero.p); PA.hat = PA.zero; PA.rhat = 0; PA.Fref = PA.zero.F;
    const double qo = x17c::QTilde(M, asimov, g, P, r), qa = x17c::QTilde(M, asimov, g, PA, r);
    const double clsA = x17c::CLsAsym(qo, qa);
    x17c::ToyCLsPoint T = x17c::ToyCLs(M, asimov, g, P, r, ntoys, seed + (unsigned long long)std::lround(1e4 * r));
    printf("x17c_toycls: m %.2f r %.5f BR %.4e CLs_toy %.5f +- %.5f (CLsb %.4f CLb %.4f, %d toys) CLs_asym %.5f q_obs %.4f\n", mX, r, r * kBR, T.cls, T.clsErr, T.clsb,
           T.clb, T.nt, clsA, T.qobs);
    f->Close();
}
