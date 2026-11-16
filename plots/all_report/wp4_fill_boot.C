// ============================================================================
// wp4_fill_boot.C  --  WP4: double bootstrap of the simulation templates (prefix wp4)
// ============================================================================
// SPEC report_all, WP4, extra test "RB" (realistic MC + bootstrap bias correction).
// One pass over the simulated e+e- rows (std_us + lep_us + eta_us, owned rows, mode A,
// w x V3 = physics code 108, selection S0, nominal mass transfer s = 0.99615, r = 1).
// For every event group g (tree, grp):
//     parent weight  k_p(g)  ~ Poisson(1),        p = 0 .. P-1
//     child  weight  k_pc(g) ~ Poisson(k_p(g)),   c = 0 .. C-1   (a bootstrap replica of the parent replica)
// (counter-based generator, reproducible). Filled per replica set: the F10 layout (5 slices x E1/E2 x 70 bins),
// F5g (5 slices, M(ee) >= 20 MeV) and the windows B, P, SBlo, SBhi, CR-U, for the components C2, C3, C4, C5,
// G_eta, G_D, G_C (= G_conv). Units nb (sum w / L_sim). Index of a replica set: parent p -> p, child (p, c) ->
// P + p C + c; set nRep = P (1 + C). Also the unweighted full sample (check against hists_v1 S0__nom).
// Output: /hdd3/przygoda/HADES/FAT_out/all_report/wp4/wp4_boot.root (TVectorD per component and layout)
//
// Run (from /home/przygoda/HADES/FAT/FAT):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_wp4", kTRUE)' \
//        'plots/all_report/wp4_fill_boot.C+O(20, 20)'
// ============================================================================
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "TFile.h"
#include "TStopwatch.h"
#include "TVectorD.h"

#include "wp4_grid.h"
#include "sim_eta_physics.h"
#include "hists_io.h"

namespace wp4boot {
using namespace wp4;
inline uint64_t Mix(uint64_t x) {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}
inline double U01(uint64_t h) { return ((Mix(h) >> 11) + 0.5) * (1.0 / 9007199254740992.0); }
// Poisson(lambda) by inversion (lambda <= ~20)
inline int PoisInv(double lam, uint64_t h) {
    if (!(lam > 0)) return 0;
    const double u = U01(h);
    double p = std::exp(-lam), c = p;
    int k = 0;
    while (u > c && k < 200) { ++k; p *= lam / k; c += p; }
    return k;
}
const int kNBC = 7;
const int kBC[kNBC] = {kC2, kC3, kC4, kC5, kGeta, kGD, kGconv};
const int kNBW = 5;
const int kBW[kNBW] = {kWB, kWP, kWSBlo, kWSBhi, kWCRU};
}  // namespace wp4boot

void wp4_fill_boot(int P = 20, int C = 20) {
    using namespace wp4boot;
    TStopwatch sw;
    Grid G; G.Init();
    std::vector<int> f10x(kNfine);
    for (int i = 0; i < kNfine; ++i) f10x[i] = FindBinE(G.f10, kFineLo + (i + 0.5) * kFineW);
    const int nRep = P * (1 + C);
    const int nF = G.nF + 2;
    const size_t sF10 = 10 * nF, sF5 = 5 * nF, sW = kNBW;
    const size_t stride = sF10 + sF5 + sW;
    // arrays [comp][rep][cell] (+ one extra replica slot nRep = the unweighted full sample)
    std::vector<std::vector<double>> A(kNBC, std::vector<double>((size_t)(nRep + 1) * stride, 0.));
    simphys::Init();
    simphys::Install();
    simphys::LoadAngles();
    simeng::SimData* D = new simeng::SimData();
    if (!simeng::Load(*D, simeng::kCompact, simeng::kMaskUS, -1e9, 1e9, 1)) return;
    const int trees[3] = {simeng::kStdUS, simeng::kLepUS, simeng::kEtaUS};
    std::vector<double> k(nRep + 1, 1.0);
    const uint64_t seedP = 0xB0075EEDULL, seedC = 0xC41DB007ULL;
    for (int t = 0; t < 3; ++t) {
        const std::vector<simeng::Row>& vr = D->t[trees[t]];
        int cur = -1;
        for (size_t i = 0; i < vr.size(); ++i) {
            const simeng::Row& r = vr[i];
            if (r.grp != cur) {
                cur = r.grp;
                const uint64_t base = ((uint64_t)trees[t] << 40) ^ (uint64_t)(uint32_t)r.grp;
                for (int p = 0; p < P; ++p) {
                    const int kp = PoisInv(1.0, Mix(seedP ^ (base * 1000003ULL + (uint64_t)p)));
                    k[p] = kp;
                    for (int c = 0; c < C; ++c)
                        k[P + p * C + c] = PoisInv((double)kp, Mix(seedC ^ (base * 1000003ULL + (uint64_t)(p * 1000 + c))));
                }
                k[nRep] = 1.0;
            }
            if (!(r.ownA == 1)) continue;
            const float dzpi = DzPi(r.pip_z, r.pim_z, r.vz);
            if (!PassSel(kSelS0, r, r.cor_ok == 1, dzpi, (float)r.mm4)) continue;
            int cl[6];
            const int nc = CompsOf(r, false, cl);
            int qi[6], nq = 0;
            for (int q = 0; q < nc; ++q) for (int b = 0; b < kNBC; ++b) if (cl[q] == kBC[b]) qi[nq++] = b;
            if (nq == 0) continue;
            const double w = r.w * simphys::RowFactor(r, kPhysNom);
            if (w == 0.) continue;
            const double m4 = M4sim(r.m4, r.m4_true, kS, kR);
            const double mpp = kSPiPi * r.mpipi;
            const int j = Slice(mpp), E = Eclass(r.mee);
            const int fi = FineBin(m4);
            const int fb = fi >= 0 ? f10x[fi] : (m4 < kFineLo ? 0 : G.nF + 1);
            int cells[4], ncell = 0;
            if (j >= 0 && E >= 0) {
                cells[ncell++] = (2 * j + E) * nF + fb;
                if (r.mee >= kF5gMee) cells[ncell++] = (int)sF10 + j * nF + fb;
            }
            for (int wv = 0; wv < kNBW; ++wv) if (InWindow(kBW[wv], m4, r.mee, mpp)) cells[ncell++] = (int)(sF10 + sF5) + wv;
            if (ncell == 0) continue;
            for (int q = 0; q < nq; ++q) {
                double* a = A[qi[q]].data();
                for (int rr = 0; rr <= nRep; ++rr) {
                    const double ww = w * k[rr];
                    if (ww == 0.) continue;
                    double* ar = a + (size_t)rr * stride;
                    for (int cc = 0; cc < ncell; ++cc) ar[cells[cc]] += ww;
                }
            }
        }
        printf("wp4_fill_boot: tree %s done, %.1f s\n", simeng::TreeName(trees[t]), sw.RealTime()); sw.Continue();
    }
    delete D;
    // check the unweighted full sample against hists_v1 (S0__nom = half A + half B)
    TFile* fh = TFile::Open(kHistFile);
    double worst = 0;
    for (int b = 0; b < kNBC; ++b) {
        const double* a = A[b].data() + (size_t)nRep * stride;
        for (int kk = 0; kk < 10; ++kk) {
            TH1* h = allio::Sim(fh, "S0", "nom", kBC[b], allio::F10(kk / 2, kk % 2));
            for (int x = 1; x <= G.nF; ++x) {
                const double m = a[kk * nF + x] / kLsim, f = h->GetBinContent(x);
                if (f != 0 || m != 0) worst = std::max(worst, std::fabs(m - f) / std::max(std::fabs(f), 1e-12));
            }
            delete h;
        }
        TH1* hw = allio::Sim(fh, "S0", "nom", kBC[b], "win");
        for (int wv = 0; wv < kNBW; ++wv) {
            const double m = a[sF10 + sF5 + wv] / kLsim, f = allio::Win(hw, kBW[wv]);
            if (f != 0 || m != 0) worst = std::max(worst, std::fabs(m - f) / std::max(std::fabs(f), 1e-12));
        }
        delete hw;
    }
    printf("wp4_fill_boot: full sample vs hists_v1 S0__nom (F10, win): worst relative difference %.3g\n", worst);
    TFile* fo = TFile::Open(Form("%s/wp4_boot.root", kOutDir), "RECREATE");
    TVectorD meta(4); meta[0] = P; meta[1] = C; meta[2] = nF; meta[3] = worst;
    fo->WriteTObject(&meta, "meta");
    for (int b = 0; b < kNBC; ++b) {
        TVectorD v((int)A[b].size());
        for (size_t i = 0; i < A[b].size(); ++i) v[(int)i] = A[b][i] / kLsim;
        fo->WriteTObject(&v, Form("%s", CompTag(kBC[b])));
    }
    fo->Close();
    printf("wp4_fill_boot: %d parents x %d children written in %.1f s\n", P, C, sw.RealTime());
}
