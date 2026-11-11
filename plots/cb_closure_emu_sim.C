// cb_closure_emu_sim.C
//
// Closure of the data-driven correction for the best-pair CB in simulation.
// Corrected estimate (as used on exp data, k = 1 here):
//    CB_best^emu(M) = sum_c eps_B,c^emu(M) * 2 sqrt(N++_c(M) N--_c(M))
// eps_B,c^emu = winning / active cross pairs (time-chi2 ranking) from
// plots/isbest_emulation_betanew.C run on the same PAT file (patterns pm, pp,
// mm; sim_trigger_bits and vertex_key on), N++_c / N--_c = native like-sign
// pairs per lepton-multiplicity class from the same runs.
// Compared with the true CB of the best pairs (time ranking) and with the
// naive estimate 2 sqrt(N++_best N--_best), both from plots/cb_closure_sim.C.
//
// Run: root -l -b -q 'plots/cb_closure_emu_sim.C+("plots/output/isbest_emu_simjan/%s_jan.root","plots/output/cb_closure_sim_jan26.root")'

#include <TFile.h>
#include <TH1.h>
#include <cmath>
#include <cstdio>

void cb_closure_emu_sim(const char* emu_pattern = "plots/output/isbest_emu_simjan/%s_jan.root",
                        const char* closure_file = "plots/output/cb_closure_sim_jan26.root")
{
    TFile fpm(Form(emu_pattern, "pm")), fpp(Form(emu_pattern, "pp")), fmm(Form(emu_pattern, "mm")), fc(closure_file);
    TH1* Bt   = (TH1*)fc.Get("Btruth_time");
    TH1* Bp   = (TH1*)fc.Get("Btruth_prod");
    TH1* Ba   = (TH1*)fc.Get("Btruth_all");
    TH1* CBnT = (TH1*)fc.Get("CBhat_time");
    TH1* CBa  = (TH1*)fc.Get("CBhat_all");
    TH1* Npm  = (TH1*)fc.Get("Npm_all");
    if (!Bt || !CBnT || !CBa) { printf("closure histograms missing\n"); return; }
    const int nb = Bt->GetNbinsX();
    TH1* emu = (TH1*)Bt->Clone("CBbest_emu"); emu->Reset();
    TH1* cls = (TH1*)Bt->Clone("CBall_perclass"); cls->Reset();
    TH1* npmEmu = (TH1*)Bt->Clone("Npm_emu_native"); npmEmu->Reset();
    double e2[64] = {0};
    for (int c = 2; c <= 5; ++c) {
        TH1* hc = (TH1*)fpm.Get(Form("h_cross_pm_PT3_c%d", c));
        TH1* hw = (TH1*)fpm.Get(Form("h_crosswin_pm_PT3_c%d", c));
        TH1* hp = (TH1*)fpp.Get(Form("h_all_pp_PT3_c%d", c));
        TH1* hm = (TH1*)fmm.Get(Form("h_all_mm_PT3_c%d", c));
        TH1* hn = (TH1*)fpm.Get(Form("h_all_pm_PT3_c%d", c));
        if (!hc || !hw || !hp || !hm) { printf("missing class %d histograms\n", c); continue; }
        if (hn) npmEmu->Add(hn);
        for (int b = 1; b <= nb; ++b) {
            const double a = hp->GetBinContent(b), m = hm->GetBinContent(b);
            if (a <= 0 || m <= 0) continue;
            const double cb = 2 * std::sqrt(a * m);
            const double ncr = hc->GetBinContent(b);
            // empty cross-pair bins: class-2 has eps = 1 by construction; otherwise use the class average
            double eps = ncr > 0 ? hw->GetBinContent(b) / ncr : (c == 2 ? 1. : hw->Integral() / std::max(1., hc->Integral()));
            cls->AddBinContent(b, cb);
            emu->AddBinContent(b, eps * cb);
            e2[b] += std::pow(eps * 0.5 * cb, 2) * (1. / a + 1. / m);
        }
    }
    const double rng[3][2] = {{0.0, 0.14}, {0.14, 0.80}, {0.0, 1.40}};
    auto I = [](TH1* h, double lo, double hi) { return h->Integral(h->FindBin(lo + 1e-6), h->FindBin(hi - 1e-6)); };
    printf("consistency: native e+e- pairs, emulation vs closure macro (all masses): %.0f vs %.0f\n", I(npmEmu, 0, 1.4), Npm ? I(Npm, 0, 1.4) : -1.);
    printf("\n   M range   | B_truth best(time) | naive 2sqrt(N++best N--best) | corrected sum_c eps_B,c CB_c | naive/truth | corrected/truth | (all comb.: per-class CB / B_truth_all)\n");
    for (auto& r : rng) {
        const double bt = I(Bt, r[0], r[1]), nv = I(CBnT, r[0], r[1]), em = I(emu, r[0], r[1]);
        double ee = 0; for (int b = Bt->FindBin(r[0] + 1e-6); b <= Bt->FindBin(r[1] - 1e-6); ++b) ee += e2[b];
        const double ba = I(Ba, r[0], r[1]), cl = I(cls, r[0], r[1]);
        printf("  %.2f-%.2f  | %10.0f +- %4.0f | %16.0f             | %14.0f +- %4.0f        | %.3f      | %.3f +- %.3f   | %.3f\n",
               r[0], r[1], bt, std::sqrt(bt), nv, em, std::sqrt(ee), bt > 0 ? nv / bt : 0., bt > 0 ? em / bt : 0.,
               bt > 0 ? (em / bt) * std::sqrt(ee / (em * em) + 1. / bt) : 0., ba > 0 ? cl / ba : 0.);
    }
    if (Bp) printf("\n(prod ranking, for reference) B_truth best(prod) MID = %.0f\n", I(Bp, 0.14, 0.80));
}
