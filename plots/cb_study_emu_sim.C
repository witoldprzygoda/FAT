// cb_study_emu_sim.C
//
// Best-pair (isBest == 1) background in the full SMASH simulation, PT3 events,
// no opening-angle cut: true background against
//   naive     : 2 sqrt(N++ N--) of isBest == 1 pairs, per mass bin,
//   truth eps : sum_c eps_B,c^truth CB_c with CB_c = 2 sqrt(N++_c N--_c) of all pairs,
//   emu eps   : the same with eps_B,c from the data-driven emulation
//               (plots/isbest_emulation_betanew.C: winning / active cross pairs),
// each without and with the k factor from event mixing done on this simulation
// (hmix_* of plots/cb_study_sim.C). Inputs: cb_study_sim_<tag>.root and the
// emulation outputs for patterns pm, pp, mm.
//
// Run: root -l -b -q 'plots/cb_study_emu_sim.C+'

#include <TFile.h>
#include <TH1.h>
#include <cmath>
#include <cstdio>

void cb_study_emu_sim(const char* study = "plots/output/cb_study_sim_std_full.root",
                      const char* emu_pattern = "plots/output/isbest_emu_simfull_v5/%s_full.root")
{
    TFile fs(study), fpm(Form(emu_pattern, "pm")), fpp(Form(emu_pattern, "pp")), fmm(Form(emu_pattern, "mm"));
    TH1* Bb = (TH1*)fs.Get("Btrue_PT3_best_oa0"); TH1* Bg = (TH1*)fs.Get("Sgen_PT3_best_oa0");
    TH1* CBn = (TH1*)fs.Get("CBhat_PT3_best_oa0"); TH1* Npm = (TH1*)fs.Get("Npm_PT3_best_oa0");
    TH1* mpm = (TH1*)fs.Get("hmix_PT3_oa0_pm"); TH1* mpp = (TH1*)fs.Get("hmix_PT3_oa0_pp"); TH1* mmm = (TH1*)fs.Get("hmix_PT3_oa0_mm");
    if (!Bb || !CBn || !mpm) { printf("missing study histograms\n"); return; }
    const int nb = Bb->GetNbinsX();
    TH1* corrT = (TH1*)Bb->Clone("corrT"); corrT->Reset(); TH1* corrE = (TH1*)Bb->Clone("corrE"); corrE->Reset();
    double chk[3] = {0, 0, 0}, chk2[3] = {0, 0, 0};
    for (int c = 2; c <= 5; ++c) {
        TH1* cr = (TH1*)fpm.Get(Form("h_cross_pm_PT3_c%d", c)); TH1* wn = (TH1*)fpm.Get(Form("h_crosswin_pm_PT3_c%d", c));
        TH1* pe = (TH1*)fpp.Get(Form("h_all_pp_PT3_c%d", c)); TH1* me = (TH1*)fmm.Get(Form("h_all_mm_PT3_c%d", c)); TH1* ne = (TH1*)fpm.Get(Form("h_all_pm_PT3_c%d", c));
        TH1* ps = (TH1*)fs.Get(Form("hc_PT3_all_oa0_pp_c%d", c)); TH1* ms = (TH1*)fs.Get(Form("hc_PT3_all_oa0_mm_c%d", c)); TH1* ns = (TH1*)fs.Get(Form("hc_PT3_all_oa0_pm_c%d", c));
        TH1* ba = (TH1*)fs.Get(Form("hcB_PT3_all_oa0_c%d", c)); TH1* bb = (TH1*)fs.Get(Form("hcB_PT3_best_oa0_c%d", c));
        if (!cr || !wn || !pe || !me || !ps || !ms || !ba || !bb) { printf("missing class %d\n", c); continue; }
        chk[0] += ne->Integral(); chk[1] += pe->Integral(); chk[2] += me->Integral(); chk2[0] += ns->Integral(); chk2[1] += ps->Integral(); chk2[2] += ms->Integral();
        for (int b = 1; b <= nb; ++b) {
            const double a = ps->GetBinContent(b), m = ms->GetBinContent(b);
            if (a <= 0 || m <= 0) continue;
            const double g = 2 * std::sqrt(a * m);
            const double nc = cr->GetBinContent(b);
            const double eE = nc > 0 ? wn->GetBinContent(b) / nc : (c == 2 ? 1. : wn->Integral() / std::max(1., cr->Integral()));
            const double bab = ba->GetBinContent(b);
            const double eT = bab > 0 ? bb->GetBinContent(b) / bab : (c == 2 ? 1. : bb->Integral() / std::max(1., ba->Integral()));
            corrE->AddBinContent(b, eE * g); corrT->AddBinContent(b, eT * g);
        }
    }
    printf("consistency of the two codes (PT3, classes summed): N+- emulation %.0f / study %.0f ; N++ %.0f / %.0f ; N-- %.0f / %.0f\n", chk[0], chk2[0], chk[1], chk2[1], chk[2], chk2[2]);
    const double R[6][2] = {{0, 0.05}, {0.05, 0.14}, {0.14, 0.30}, {0.30, 0.56}, {0.0, 0.14}, {0.14, 0.80}};
    auto sum = [&](TH1* h, double lo, double hi, bool withK) { double s = 0; for (int b = h->FindBin(lo + 1e-6); b <= h->FindBin(hi - 1e-6); ++b) { double k = 1.;
            if (withK) { const double a = mpp->GetBinContent(b), m = mmm->GetBinContent(b); if (a > 0 && m > 0) k = mpm->GetBinContent(b) / (2 * std::sqrt(a * m)); } s += k * h->GetBinContent(b); } return s; };
    printf("\nRatio to the true best-pair background (PT3, OA>0), and signal N+- - CB over the true signal:\n");
    printf("  M range   | B_true(best) | naive | eps truth | eps emu | eps truth x k_mix | eps emu x k_mix || S_hat/S_true: naive | eps truth | eps emu | eps emu x k_mix\n");
    for (auto& r : R) {
        const double B = sum(Bb, r[0], r[1], false), N = sum(Npm, r[0], r[1], false), S = sum(Bg, r[0], r[1], false);
        const double n = sum(CBn, r[0], r[1], false), t = sum(corrT, r[0], r[1], false), e = sum(corrE, r[0], r[1], false), tk = sum(corrT, r[0], r[1], true), ek = sum(corrE, r[0], r[1], true);
        printf("  %.2f-%.2f | %8.0f | %.3f | %.3f | %.3f | %.3f | %.3f || %.4f | %.4f | %.4f | %.4f\n", r[0], r[1], B, n / B, t / B, e / B, tk / B, ek / B, (N - n) / S, (N - t) / S, (N - e) / S, (N - ek) / S);
    }
    printf("\neps_B,c in 0.14-0.80 (PT3): class | truth | emulation\n");
    for (int c = 3; c <= 5; ++c) { TH1* cr = (TH1*)fpm.Get(Form("h_cross_pm_PT3_c%d", c)); TH1* wn = (TH1*)fpm.Get(Form("h_crosswin_pm_PT3_c%d", c)); TH1* ba = (TH1*)fs.Get(Form("hcB_PT3_all_oa0_c%d", c)); TH1* bb = (TH1*)fs.Get(Form("hcB_PT3_best_oa0_c%d", c));
        auto s = [&](TH1* h) { return h->Integral(h->FindBin(0.14 + 1e-6), h->FindBin(0.80 - 1e-6)); };
        printf("  n=%d%s | %.4f (%.0f pairs) | %.4f\n", c, c == 5 ? "+" : " ", s(bb) / s(ba), s(ba), s(wn) / s(cr)); }
}
