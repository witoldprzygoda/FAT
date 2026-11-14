// ============================================================================
// s5_blind_prep_varcheck.C  --  validation of the Poisson bootstrap of the eta pseudo-data with many toys
// ============================================================================
// Includes s5_blind_prep.C (unchanged) and reuses its candidate rows of half B. Generates N toys with an INDEPENDENT
// random seed (std::random_device; nothing is written to the blind directories, no toy is stored) and compares, per
// frozen selection and M4 region, the mean and the variance of the TOTAL number of rows per toy (all groups together:
// observable quantities, no composition) with the bootstrap prediction: mean = half-B sum of w, variance V_coh =
// sum over groups of sum_rs min(w_r, w_s) (coherent thinning). The error of the toy variance uses the fourth central
// moment (the per-toy counts are not Gaussian when events contribute several rows): Var(s^2) = (m4 - (n-3)/(n-1) m2^2)/n.
// Output: plots/eta_sim/tables/s5_blind_prep_varcheck.txt. SIMULATION ONLY.
// Run from /home/przygoda/HADES/FAT/FAT (about 2 min for 2000 toys, 10 GB):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/aclic_s5_blind_prep", kTRUE)' \
//        'plots/eta_sim/s5_blind_prep_varcheck.C+(2000)'
// ============================================================================
#include "s5_blind_prep.C"

namespace s5bv {
using namespace s5b;

void Run(int nToy) {
    TStopwatch sw;
    if (!SetupCuts()) return;
    if (!simphys::Init()) return;
    simphys::Install();
    if (!simphys::LoadAngles()) return;
    gO = SimOpt();
    gO.phys = simphys::kPhysMain;
    gD = new SimData;
    if (!Load(*gD, kCompact, kMaskUS, kM4Lo, kM4Hi, 1)) return;
    BuildCandidates();
    Expectations();
    std::random_device rd;
    std::mt19937_64 rng(((ULong64_t)rd() << 32) ^ (ULong64_t)rd());
    // per toy: totals per selection x region
    std::vector<double> tot((size_t)nToy * NSW * NRX, 0.);
    auto T = [&](int t, int s, int k) -> double& { return tot[((size_t)t * NSW + s) * NRX + k]; };
    for (int t = 0; t < nToy; ++t) {
        for (const CGrp& G : gG) {
            const int n = Pois(rng, G.wg, G.L);
            for (int c = 0; c < n; ++c) {
                const double u = U01(rng);
                for (int i = G.first; i < G.first + G.n; ++i) {
                    const Cand& x = gC[i];
                    if (!(x.ratio > u)) continue;
                    for (int s = 0; s < NSW; ++s) {
                        if (!(x.bits & (1 << s))) continue;
                        if (x.reg >= 0) T(t, s, x.reg) += 1;
                        if (x.xwin) T(t, s, NR) += 1;
                    }
                }
            }
        }
        if (t % 250 == 0) { printf("toy %d (%.0f s)\n", t, sw.RealTime()); sw.Continue(); }
    }
    FILE* f = fopen((std::string(kTabDir) + "/s5_blind_prep_varcheck.txt").c_str(), "w");
    auto P = [&](const char* fmt, auto... a) { printf(fmt, a...); fprintf(f, fmt, a...); };
    P("# s5_blind_prep_varcheck: validation of the Poisson bootstrap (coherent thinning) of the eta pseudo-data, %d toys\n", nToy);
    P("# independent random seed, nothing stored; totals of ALL rows per toy (observable), per frozen selection and M4 region\n");
    P("# prediction: mean = half-B sum of w; V_coh = sum_groups sum_rs min(w_r, w_s); V_ind = independent lepton pairs;\n");
    P("# err(var) from the 4th central moment; pull_var = (var - V_coh) / err(var)\n");
    P("%-7s %-8s %10s %10s %7s %11s %10s %9s %8s %8s %8s\n", "sel", "region", "toy mean", "predicted", "pull", "toy var", "V_coh",
      "err(var)", "var/Vc", "pull_var", "Vind/Vc");
    double chi2m = 0, chi2v = 0; int nd = 0;
    for (int s = 0; s < NSW; ++s)
        for (int k = 0; k < NRX; ++k) {
            double m = 0;
            for (int t = 0; t < nToy; ++t) m += T(t, s, k);
            m /= nToy;
            double m2 = 0, m4 = 0;
            for (int t = 0; t < nToy; ++t) { const double d = T(t, s, k) - m; m2 += d * d; m4 += d * d * d * d; }
            const double var = m2 / (nToy - 1);
            m2 /= nToy; m4 /= nToy;
            const double evar = std::sqrt(std::max(0., (m4 - (nToy - 3.) / (nToy - 1.) * m2 * m2) / nToy));
            const double vc = gVcoh[s][k];
            const double pm = (m - gMuTot[s][k]) / std::sqrt(vc / nToy), pv = evar > 0 ? (var - vc) / evar : 0.;
            P("%-7s %-8s %10.2f %10.2f %7.2f %11.1f %10.1f %9.1f %8.3f %8.2f %8.4f\n", kSN[s], RegShort(k).c_str(), m, gMuTot[s][k], pm, var,
              vc, evar, var / vc, pv, gVind[s][k] / vc);
            if (k < NR && s != sDIRECT) { chi2m += pm * pm; chi2v += pv * pv; ++nd; }   // DIRECT shares most rows with PEAK
        }
    P("\nsum of pull^2 over the disjoint regions of PEAK, BUMP, X17SEL (%d entries; regions of one selection nearly independent,\n"
      "the selections share rows): means %.1f, variances %.1f\n", nd, chi2m, chi2v);
    P("time %.0f s\n", sw.RealTime());
    fclose(f);
}
}  // namespace s5bv

void s5_blind_prep_varcheck(int nToy = 2000) { s5bv::Run(nToy); }
