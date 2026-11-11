// cb_study_sim_plot.C
//
// Figures for plots/cb_study_sim.C (full SMASH simulation with EpEm_ID, EpEp_ID,
// EmEm_ID): (1) all / CB / signal for isBest >= 0 and isBest == 1 with the true
// CB and true signal; (2) k factors; (3) production channels and CB composition;
// (4) best-pair background estimates against the truth.
//
// Run: root -l -b -q 'plots/cb_study_sim_plot.C+'

#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TString.h>
#include <TStyle.h>
#include <cmath>
#include <cstdio>
#include <vector>

namespace csp {
void axes(TH1* h, const char* y, double ymin, double ymax, double xmax) {
    h->SetTitle(Form(";M_{ee} [GeV/c^{2}];%s", y)); h->SetMinimum(ymin); h->SetMaximum(ymax); h->GetXaxis()->SetRangeUser(0, xmax);
    h->GetXaxis()->SetTitleSize(0.055); h->GetYaxis()->SetTitleSize(0.055); h->GetXaxis()->SetLabelSize(0.05); h->GetYaxis()->SetLabelSize(0.05); h->GetYaxis()->SetTitleOffset(1.15); h->GetXaxis()->SetTitleOffset(1.05);
}
TH1D* frame(const char* n, const char* y, double ymin, double ymax, double xmax) { TH1D* f = new TH1D(n, "", 10, 0, xmax); axes(f, y, ymin, ymax, xmax); return f; }
void mk(TH1* h, int c, int m, double s = 0.8) { h->SetMarkerColor(c); h->SetLineColor(c); h->SetMarkerStyle(m); h->SetMarkerSize(s); }
void pad(bool logy) { gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.08); if (logy) gPad->SetLogy(); else gPad->SetGridy(); }
void title(const char* t) { TLatex l; l.SetNDC(); l.SetTextSize(0.05); l.DrawLatex(0.14, 0.935, t); }
void unity(double xmax) { TLine* l = new TLine(0, 1, xmax, 1); l->SetLineStyle(2); l->SetLineWidth(2); l->Draw(); }
const std::vector<double> kCoarse = {0, 0.005, 0.01, 0.02, 0.03, 0.05, 0.08, 0.11, 0.14, 0.16, 0.18, 0.20, 0.22, 0.24, 0.26, 0.28, 0.30, 0.34, 0.38, 0.42, 0.46, 0.50, 0.56, 0.60, 0.70, 0.80};
// ratio a/b on coarse bins with sqrt(N) errors of a and b
TGraphErrors* ratioG(TH1* a, TH1* b, int c, int m, double minN = 30) {
    TGraphErrors* g = new TGraphErrors();
    for (size_t i = 0; i + 1 < kCoarse.size(); ++i) { const int b1 = a->FindBin(kCoarse[i] + 1e-6), b2 = a->FindBin(kCoarse[i + 1] - 1e-6);
        double A = 0, B = 0, eA = 0, eB = 0; for (int k = b1; k <= b2; ++k) { A += a->GetBinContent(k); B += b->GetBinContent(k); eA += std::pow(a->GetBinError(k), 2); eB += std::pow(b->GetBinError(k), 2); }
        if (B < minN || A <= 0) continue;
        const int n = g->GetN(); const double r = A / B;
        g->SetPoint(n, 0.5 * (kCoarse[i] + kCoarse[i + 1]), r); g->SetPointError(n, 0.5 * (kCoarse[i + 1] - kCoarse[i]), r * std::sqrt(eA / (A * A) + eB / (B * B)));
    }
    g->SetMarkerColor(c); g->SetLineColor(c); g->SetMarkerStyle(m); g->SetMarkerSize(1.0);
    return g;
}
// k = N / 2 sqrt(P M) on coarse bins (sum of spectra inside a coarse bin)
TGraphErrors* kG(TH1* n, TH1* p, TH1* m, int c, int mk, double minN = 30) {
    TGraphErrors* g = new TGraphErrors();
    for (size_t i = 0; i + 1 < kCoarse.size(); ++i) { const int b1 = n->FindBin(kCoarse[i] + 1e-6), b2 = n->FindBin(kCoarse[i + 1] - 1e-6);
        const double N = n->Integral(b1, b2), P = p->Integral(b1, b2), M = m->Integral(b1, b2); if (N < minN || P <= 0 || M <= 0) continue;
        const int ip = g->GetN(); const double k = N / (2 * std::sqrt(P * M));
        g->SetPoint(ip, 0.5 * (kCoarse[i] + kCoarse[i + 1]), k); g->SetPointError(ip, 0.5 * (kCoarse[i + 1] - kCoarse[i]), k * std::sqrt(1 / N + 0.25 / P + 0.25 / M)); }
    g->SetMarkerColor(c); g->SetLineColor(c); g->SetMarkerStyle(mk); g->SetMarkerSize(1.0);
    return g;
}
TH1D* get(TFile& f, const char* n) { TH1D* h = (TH1D*)f.Get(n); if (!h) printf("missing %s\n", n); else { h = (TH1D*)h->Clone(Form("%s_c", n)); h->SetDirectory(nullptr); } return h; }
TH1D* sumCats(TFile& f, const char* set, const char* v, double oa, const char* p, std::vector<int> cats) { TH1D* s = nullptr; for (int c : cats) { TH1D* h = get(f, Form("h_%s_%s_oa%g_%s_cat%d", set, v, oa, p, c)); if (!h) continue; if (!s) s = h; else s->Add(h); } return s; }
}

void cb_study_sim_plot(const char* study = "plots/output/cb_study_sim_std_full.root", const char* mixdata = "plots/output/mixing_betanew.root",
                       const char* emu_pattern = "plots/output/isbest_emu_simfull_v5/%s_full.root", const char* out = "plots/output/cb_study_sim_std_full")
{
    using namespace csp;
    gStyle->SetOptStat(0);
    TFile fs(study);
    const char* set = "all"; const double OA = 2.;
    const char* vs[2] = {"all", "best"}; const char* vt[2] = {"isBest #geq 0, OA > 2#circ", "isBest = 1, OA > 2#circ"};

    // ---------- (1) all / CB / signal ----------
    TCanvas* c1 = new TCanvas("c1", "", 1500, 1000); c1->Divide(2, 2);
    for (int v = 0; v < 2; ++v) {
        TH1D* N = get(fs, Form("Npm_%s_%s_oa%g", set, vs[v], OA)), *CB = get(fs, Form("CBhat_%s_%s_oa%g", set, vs[v], OA)), *B = get(fs, Form("Btrue_%s_%s_oa%g", set, vs[v], OA)), *S = get(fs, Form("Sgen_%s_%s_oa%g", set, vs[v], OA));
        TH1D* Sh = (TH1D*)N->Clone(Form("Sh%d", v)); Sh->Add(CB, -1.); for (int b = 1; b <= Sh->GetNbinsX(); ++b) Sh->SetBinError(b, std::sqrt(N->GetBinContent(b) + std::pow(CB->GetBinError(b), 2)));
        for (TH1D* h : {B, S}) for (int b = 1; b <= h->GetNbinsX(); ++b) h->SetBinError(b, std::sqrt(h->GetBinContent(b)));
        c1->cd(v + 1); pad(true); mk(N, kBlack, 20); mk(CB, kRed + 1, 21); mk(B, kRed + 1, 25, 1.0); mk(Sh, kBlue + 1, 22); mk(S, kBlue + 1, 26, 1.0);
        axes(N, "pairs / bin", 5., 3e6, 1.0); N->Draw("E"); B->Draw("E same"); CB->Draw("E same"); S->Draw("E same"); Sh->Draw("E same");
        TLegend* lg = new TLegend(0.44, 0.62, 0.97, 0.90); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.042);
        lg->AddEntry(N, "all: N_{+-}", "p"); lg->AddEntry(CB, "CB = 2#sqrt{N_{++}N_{--}}", "p"); lg->AddEntry(B, "CB true (particle history)", "p"); lg->AddEntry(Sh, "signal = all - CB", "p"); lg->AddEntry(S, "signal true", "p"); lg->Draw(); title(vt[v]);
        c1->cd(v + 3); pad(false); frame(Form("f1_%d", v), "ratio", 0.6, 1.6, 1.0)->Draw(); unity(1.0);
        TGraphErrors* r1 = ratioG(CB, B, kRed + 1, 21), *r2 = ratioG(Sh, S, kBlue + 1, 22); r1->Draw("P same"); r2->Draw("P same");
        TLegend* lr = new TLegend(0.42, 0.74, 0.97, 0.90); lr->SetBorderSize(0); lr->SetFillStyle(0); lr->SetTextSize(0.045); lr->AddEntry(r1, "CB formula / CB true", "p"); lr->AddEntry(r2, "signal (all - CB) / signal true", "p"); lr->Draw(); title(vt[v]);
    }
    c1->SaveAs(Form("%s_signal_cb.pdf", out)); c1->SaveAs(Form("%s_signal_cb.png", out));

    // ---------- (2) k factors, PT3 ----------
    TFile fd(mixdata);
    TH1D* dpm = get(fd, "h_mixnc_pm_PT3"), *dpp = get(fd, "h_mixnc_pp_PT3"), *dmm = get(fd, "h_mixnc_mm_PT3");
    TH1D* spm = get(fs, "hmix_PT3_oa0_pm"), *spp = get(fs, "hmix_PT3_oa0_pp"), *smm = get(fs, "hmix_PT3_oa0_mm");
    TH1D* B = get(fs, "Btrue_PT3_all_oa0"), *Np = get(fs, "Npp_PT3_all_oa0"), *Nm = get(fs, "Nmm_PT3_all_oa0");
    TH1D* Bl = get(fs, "Blep_PT3_all_oa0"), *Nlp = get(fs, "NppLep_PT3_all_oa0"), *Nlm = get(fs, "NmmLep_PT3_all_oa0");
    TCanvas* c2 = new TCanvas("c2", "", 1500, 1000); c2->Divide(2, 2);
    TGraphErrors* gT = kG(B, Np, Nm, kRed + 1, 25), *gL = kG(Bl, Nlp, Nlm, kBlue + 1, 21), *gS = kG(spm, spp, smm, kGreen + 2, 24), *gD = kG(dpm, dpp, dmm, kBlack, 20);
    const double xm[2] = {0.8, 0.2};
    for (int p = 0; p < 2; ++p) { c2->cd(p + 1); pad(false); frame(Form("f2_%d", p), "k = N_{+-} / 2#sqrt{N_{++}N_{--}}", 0.7, 1.4, xm[p])->Draw(); unity(xm[p]); gT->Draw("P same"); gS->Draw("P same"); gD->Draw("P same"); gL->Draw("P same");
        if (p == 0) { TLegend* lg = new TLegend(0.30, 0.68, 0.97, 0.90); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.042); lg->AddEntry(gL, "k_{true}: simulation, true leptons", "p"); lg->AddEntry(gT, "k_{true}: simulation, with fake legs", "p"); lg->AddEntry(gS, "k_{mix}: event mixing, simulation", "p"); lg->AddEntry(gD, "k_{mix}: event mixing, data", "p"); lg->Draw(); }
        title(p == 0 ? "PT3, all combinations, no OA cut" : "low-mass zoom");
        c2->cd(p + 3); pad(false); frame(Form("f3_%d", p), "r = k_{true} / k_{mix,sim}", 0.7, 1.4, xm[p])->Draw(); unity(xm[p]);
        TGraphErrors* gr = new TGraphErrors(); for (int i = 0; i < gL->GetN(); ++i) { double x = gL->GetX()[i]; for (int j = 0; j < gS->GetN(); ++j) if (std::fabs(gS->GetX()[j] - x) < 1e-6) { const double r = gL->GetY()[i] / gS->GetY()[j]; int n = gr->GetN(); gr->SetPoint(n, x, r); gr->SetPointError(n, gL->GetEX()[i], r * std::sqrt(std::pow(gL->GetEY()[i] / gL->GetY()[i], 2) + std::pow(gS->GetEY()[j] / gS->GetY()[j], 2))); } }
        gr->SetMarkerColor(kBlue + 1); gr->SetLineColor(kBlue + 1); gr->SetMarkerStyle(21); gr->Draw("P same"); title("true leptons"); }
    c2->SaveAs(Form("%s_k.pdf", out)); c2->SaveAs(Form("%s_k.png", out));

    // ---------- (3) channels and CB composition ----------
    TCanvas* c3 = new TCanvas("c3", "", 1800, 1000); c3->Divide(3, 2);
    const int chC[4] = {kBlue + 1, kMagenta + 1, kOrange + 1, kCyan + 2}, chM[4] = {20, 21, 22, 23}; std::vector<std::vector<int>> chCats = {{0}, {1}, {2}, {3, 4}}; const char* chN[4] = {"#pi^{0} Dalitz", "#eta Dalitz", "conversion, #gamma from #pi^{0}", "conversion, #gamma from #eta and other"};
    const double oas[3] = {0., 2., 9.}; const char* oat[3] = {"no OA cut, isBest #geq 0", "OA > 2#circ, isBest #geq 0", "OA > 9#circ, isBest #geq 0"};
    for (int i = 0; i < 3; ++i) { c3->cd(i + 1); pad(true); TH1D* f = frame(Form("f4_%d", i), "genuine pairs / bin", 5, 2e6, 0.7); f->Draw(); TLegend* lg = new TLegend(0.38, 0.66, 0.97, 0.90); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.04);
        for (int k = 0; k < 4; ++k) { TH1D* h = sumCats(fs, set, "all", oas[i], "pm", chCats[k]); if (!h) continue; mk(h, chC[k], chM[k]); h->Draw("E same"); lg->AddEntry(h, chN[k], "p"); } lg->Draw(); title(oat[i]); }
    c3->cd(4); pad(false); frame("f5", "isBest = 1 / all, OA > 2#circ", 0.0, 1.1, 0.7)->Draw(); unity(0.7); { TLegend* lg = new TLegend(0.18, 0.18, 0.75, 0.46); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.04);
        for (int k = 0; k < 4; ++k) { TH1D* a = sumCats(fs, set, "all", 2., "pm", chCats[k]), *b = sumCats(fs, set, "best", 2., "pm", chCats[k]); TGraphErrors* g = ratioG(b, a, chC[k], chM[k]); g->Draw("P same"); lg->AddEntry(g, chN[k], "p"); } lg->Draw(); title("genuine pairs kept by isBest = 1"); }
    std::vector<std::vector<int>> cbCats = {{6}, {7}, {8}, {9, 10}, {11, 12, 13}}; const char* cbN[5] = {"Dalitz + Dalitz", "Dalitz + conversion", "conversion + conversion", "with secondary e, other", "with fake leg"}; const int cbC[5] = {kRed + 1, kOrange + 7, kMagenta + 1, kGreen + 2, kGray + 2};
    c3->cd(5); pad(true); frame("f6", "true CB e^{+}e^{-} pairs / bin", 5, 2e5, 0.7)->Draw(); { TLegend* lg = new TLegend(0.42, 0.62, 0.97, 0.90); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.04);
        for (int k = 0; k < 5; ++k) { TH1D* h = sumCats(fs, set, "all", 2., "pm", cbCats[k]); mk(h, cbC[k], 20 + k); h->Draw("E same"); lg->AddEntry(h, cbN[k], "p"); } lg->Draw(); title("true CB by source, OA > 2#circ, isBest #geq 0"); }
    c3->cd(6); pad(false); frame("f7", "isBest = 1 / all, OA > 2#circ", 0.0, 1.1, 0.7)->Draw(); unity(0.7); { TLegend* lg = new TLegend(0.18, 0.18, 0.80, 0.50); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.04);
        for (int k = 0; k < 3; ++k) { TH1D* a = sumCats(fs, set, "all", 2., "pm", cbCats[k]), *b = sumCats(fs, set, "best", 2., "pm", cbCats[k]); TGraphErrors* g = ratioG(b, a, cbC[k], 20 + k); g->Draw("P same"); lg->AddEntry(g, Form("e^{+}e^{-} %s", cbN[k]), "p"); }
        for (int k = 0; k < 2; ++k) { TH1D* a = get(fs, Form("%s_%s_all_oa2", k == 0 ? "Npp" : "Nmm", set)), *b = get(fs, Form("%s_%s_best_oa2", k == 0 ? "Npp" : "Nmm", set)); TGraphErrors* g = ratioG(b, a, k == 0 ? kBlack : kGray + 1, 24 + k); g->Draw("P same"); lg->AddEntry(g, k == 0 ? "e^{+}e^{+}, all origins" : "e^{-}e^{-}, all origins", "p"); } lg->Draw(); title("CB and like-sign pairs kept by isBest = 1"); }
    c3->SaveAs(Form("%s_channels.pdf", out)); c3->SaveAs(Form("%s_channels.png", out));

    // ---------- (4) best-pair background estimates ----------
    TFile fpm(Form(emu_pattern, "pm")), fpp(Form(emu_pattern, "pp")), fmm(Form(emu_pattern, "mm"));
    TH1D* Bb = get(fs, "Btrue_PT3_best_oa0"), *CBn = get(fs, "CBhat_PT3_best_oa0"); TH1D* corrT = (TH1D*)Bb->Clone("corrT4"); corrT->Reset(); TH1D* corrE = (TH1D*)Bb->Clone("corrE4"); corrE->Reset(); TH1D* corrK = (TH1D*)Bb->Clone("corrK4"); corrK->Reset();
    for (int b = 1; b <= Bb->GetNbinsX(); ++b) { double kmx = 1.; { const double a = spp->GetBinContent(b), m = smm->GetBinContent(b); if (a > 0 && m > 0) kmx = spm->GetBinContent(b) / (2 * std::sqrt(a * m)); }
        for (int cl = 2; cl <= 5; ++cl) { TH1* cr = (TH1*)fpm.Get(Form("h_cross_pm_PT3_c%d", cl)); TH1* wn = (TH1*)fpm.Get(Form("h_crosswin_pm_PT3_c%d", cl)); TH1* ps = (TH1*)fs.Get(Form("hc_PT3_all_oa0_pp_c%d", cl)); TH1* ms = (TH1*)fs.Get(Form("hc_PT3_all_oa0_mm_c%d", cl)); TH1* ba = (TH1*)fs.Get(Form("hcB_PT3_all_oa0_c%d", cl)); TH1* bb = (TH1*)fs.Get(Form("hcB_PT3_best_oa0_c%d", cl)); if (!cr || !wn || !ps || !ms || !ba || !bb) continue;
            const double a = ps->GetBinContent(b), m = ms->GetBinContent(b); if (a <= 0 || m <= 0) continue; const double g = 2 * std::sqrt(a * m); const double nc = cr->GetBinContent(b), bab = ba->GetBinContent(b);
            const double eE = nc > 0 ? wn->GetBinContent(b) / nc : (cl == 2 ? 1. : wn->Integral() / std::max(1., cr->Integral())); const double eT = bab > 0 ? bb->GetBinContent(b) / bab : (cl == 2 ? 1. : bb->Integral() / std::max(1., ba->Integral()));
            corrE->AddBinContent(b, eE * g); corrT->AddBinContent(b, eT * g); corrK->AddBinContent(b, eE * g * kmx); } }
    for (TH1D* h : {corrT, corrE, corrK, CBn, Bb}) for (int b = 1; b <= h->GetNbinsX(); ++b) h->SetBinError(b, std::sqrt(std::max(0., h->GetBinContent(b))));
    TCanvas* c4 = new TCanvas("c4", "", 1000, 650); pad(false); gPad->SetLeftMargin(0.11); frame("f8", "best-pair CB estimate / true CB", 0.7, 1.6, 0.8)->Draw(); unity(0.8);
    TGraphErrors* g1 = ratioG(CBn, Bb, kRed + 1, 21), *g2 = ratioG(corrT, Bb, kBlack, 24), *g3 = ratioG(corrE, Bb, kBlue + 1, 22), *g4 = ratioG(corrK, Bb, kGreen + 2, 20); g1->Draw("P same"); g2->Draw("P same"); g3->Draw("P same"); g4->Draw("P same");
    TLegend* lg = new TLegend(0.40, 0.68, 0.97, 0.90); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.04); lg->AddEntry(g1, "naive 2#sqrt{N_{++}N_{--}} of isBest = 1 pairs", "p"); lg->AddEntry(g3, "#Sigma_{c} #varepsilon_{B,c}^{emu} 2#sqrt{N_{++,c}N_{--,c}}", "p"); lg->AddEntry(g2, "#Sigma_{c} #varepsilon_{B,c}^{truth} 2#sqrt{N_{++,c}N_{--,c}}", "p"); lg->AddEntry(g4, "#varepsilon_{B,c}^{emu} formula #times k_{mix,sim}", "p"); lg->Draw(); title("PT3, isBest = 1, no OA cut");
    c4->SaveAs(Form("%s_best_cb.pdf", out)); c4->SaveAs(Form("%s_best_cb.png", out));
}
