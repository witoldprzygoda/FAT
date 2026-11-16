// x17p26_p0.C -- X17 PDG 2026 (ruling D25), prefix x17p26, 2026-11-15: local p0 / Z of normalisation (1) at 17.0 and 22.0 MeV/c2 with the
// frozen sX1 and with the PDG 2026 BR_dir uncertainty in sX1 (check that Z does not depend on the PDG edition). Reads only.
// Run: root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_x17p26", kTRUE)' 'plots/all_report/x17p26_p0.C+O'
#include <cmath>
#include <cstdio>
#include "x17d_lib.h"
void x17p26_p0() {
    using namespace x17d;
    Ctx c;
    if (!Init(c)) { printf("x17p26_p0: init failed\n"); return; }
    Frozen z26 = c.z;
    z26.xSig1 = std::sqrt(c.z.xSig1 * c.z.xSig1 - std::pow(0.11 / 2.68, 2) + std::pow(0.16 / 2.76, 2));
    for (double m : {17.0, 22.0}) {
        x17c::Inputs Im = c.I;
        SetLine(Im, m);
        for (int v = 0; v < 2; ++v) {
            x17c::Config C = MakeConfig(v == 0 ? c.z : z26, 1);
            x17c::Model M(Im, C);
            P0 p = LocalP0(M, c.D.n, x17c::GlobOf(C));
            printf("x17p26_p0: m %.2f sX1 %.6f: r_hat %.5f +- %.5f, q0 %.4f, Z %.4f, p0 %.4f\n", m, C.xSig, p.rhat, p.rerr, p.q0, p.Z, p.p0);
        }
    }
    c.f->Close();
}
