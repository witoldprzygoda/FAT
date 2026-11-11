// isbest_emu_vs_truth_sim.C
//
// Closure test of the data-only isBest emulation in simulation.
// The emulation (plots/isbest_emulation_betanew.C, run on the raw SMASH PAT
// files with sim_trigger_bits = true) measures eps_B,c: the probability that a
// combinatorial e+e- pair wins the TIME-chi2 ranking, built from cross pairs
// (a leg replaced by a donor leg from another event). In simulation the same
// probability is known from MC truth (plots/isbest_truth_sim.C, ranking
// "time": argmin of the stored relative-time chi2). Both use PT3 events, the
// FAT event selection and the same lepton-multiplicity classes.
//
// Run: root -l -b -q 'plots/isbest_emu_vs_truth_sim.C+("emu.root","truth.root")'

#include <TFile.h>
#include <TH1.h>
#include <cmath>
#include <cstdio>

void isbest_emu_vs_truth_sim(const char* emu_file = "plots/output/isbest_emu_simstd_pm.root",
                             const char* truth_file = "plots/output/isbest_truth_sim_std_pt3.root")
{
    TFile fe(emu_file), ft(truth_file);
    const double rng[4][2] = {{0.0, 0.14}, {0.14, 0.80}, {0.14, 0.30}, {0.30, 0.80}};
    printf("eps_B,c = probability that a combinatorial e+e- pair wins the time-chi2 ranking (PT3)\n");
    printf("emulation: winning / active cross pairs;  truth: chosen / all combinatorial pairs (MC truth)\n");
    for (auto& r : rng) {
        printf("\nM %.2f-%.2f GeV/c2\n   class | emulation              | truth (time ranking)    | ratio emu/truth\n", r[0], r[1]);
        for (int c = 2; c <= 5; ++c) {
            TH1* hc = (TH1*)fe.Get(Form("h_cross_pm_PT3_c%d", c));
            TH1* hw = (TH1*)fe.Get(Form("h_crosswin_pm_PT3_c%d", c));
            TH1* ha = (TH1*)ft.Get(Form("h_all_comb_c%d", c));
            TH1* hs = (TH1*)ft.Get(Form("h_time_comb_c%d", c));
            if (!hc || !hw || !ha || !hs) { printf("   missing histograms for class %d\n", c); continue; }
            auto integ = [&](TH1* h) { return h->Integral(h->FindBin(r[0] + 1e-6), h->FindBin(r[1] - 1e-6)); };
            const double nc = integ(hc), nw = integ(hw), na = integ(ha), ns = integ(hs);
            const double ee = nc > 0 ? nw / nc : 0., et = na > 0 ? ns / na : 0.;
            const double dee = nc > 0 ? std::sqrt(std::max(ee * (1 - ee), 1. / nc) / nc) : 0.;
            const double det = na > 0 ? std::sqrt(std::max(et * (1 - et), 1. / na) / na) : 0.;
            printf("   n=%d%s | %.3f +- %.3f (%9.0f) | %.3f +- %.3f (%8.0f) | %.3f\n", c, c == 5 ? "+" : " ",
                   ee, dee, nc, et, det, na, et > 0 ? ee / et : 0.);
        }
    }
}
