// k_factor_compare_note.C
//
// k from event mixing in the p+p 4.53 GeV data (plots/k_mix_pools.C, PT3, pool of all events, raw pair counts) with and
// without OA > 9 deg, against Fig. 18 of the analysis note of K. Scharmann (p+p 1.58 GeV, 15 Oct 2026). The note gives no
// table: its points below were read from the figure by eye (about +- 0.01). Differences of the note: 1.58 GeV, OA > 9 deg,
// acceptance- and efficiency-corrected mixed spectra, five vertex-z slices in the mixer, theta-p acceptance window.
// k = sum N+-_mix / sum 2 sqrt(N++_mix N--_mix) over 5 MeV bins.
//
// Run: root -l -b -q 'plots/k_factor_compare_note.C+'
#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <cmath>
#include <cstdio>
#include <vector>

namespace kcn {
struct K { double k, e; };
K kOf(TH1* pm, TH1* pp, TH1* mm, double lo, double hi) {
    const int b1 = pm->FindBin(lo + 1e-6), b2 = pm->FindBin(hi - 1e-6); double n = 0, g = 0, sp = 0, sm = 0;
    for (int b = b1; b <= b2; ++b) { const double p = pp->GetBinContent(b), m = mm->GetBinContent(b); n += pm->GetBinContent(b); sp += p; sm += m; if (p > 0 && m > 0) g += 2 * std::sqrt(p * m); }
    if (g <= 0 || n <= 0) return {0, 0};
    const double k = n / g; return {k, k * std::sqrt(1 / n + 0.25 / std::max(sp, 1.) + 0.25 / std::max(sm, 1.))};
}
TGraphErrors* graph(TH1* pm, TH1* pp, TH1* mm, const std::vector<double>& e) { TGraphErrors* g = new TGraphErrors(); for (size_t i = 0; i + 1 < e.size(); ++i) { K k = kOf(pm, pp, mm, e[i], e[i + 1]); if (k.k <= 0) continue; const int n = g->GetN(); g->SetPoint(n, 500 * (e[i] + e[i + 1]), k.k); g->SetPointError(n, 500 * (e[i + 1] - e[i]), k.e); } return g; }
}

void k_factor_compare_note(const char* f9 = "plots/output/k_mix_pools_data_oa9.root", const char* f0 = "plots/output/k_mix_pools_data_oa0.root", const char* out = "plots/output/k_factor_compare_note")
{
    using namespace kcn; gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    TFile F9(f9), F0(f0);
    TH1* a[3], *b[3]; const char* pn[3] = {"pm", "pp", "mm"}; for (int i = 0; i < 3; ++i) { a[i] = (TH1*)F9.Get(Form("hf_all_%s", pn[i])); b[i] = (TH1*)F0.Get(Form("hf_all_%s", pn[i])); if (!a[i] || !b[i]) { printf("missing histograms\n"); return; } }
    // ---- points read from Fig. 18 of the note: M [MeV], k
    const std::vector<double> nx = {5, 18, 30, 42, 55, 65, 75, 88, 100, 112, 140, 180, 225, 260, 290, 320, 350, 380, 450, 550};
    const std::vector<double> ny = {0.685, 0.77, 0.825, 0.87, 0.925, 0.955, 0.98, 1.007, 1.047, 1.075, 1.08, 1.05, 1.022, 1.01, 1.005, 1.0, 1.0, 0.998, 1.0, 1.0};
    const std::vector<double> e = {0, .01, .02, .03, .04, .05, .06, .07, .08, .09, .10, .11, .12, .16, .20, .24, .28, .32, .40, .60};
    TGraphErrors* g9 = graph(a[0], a[1], a[2], e), *g0 = graph(b[0], b[1], b[2], e);
    printf(" M [MeV]    | note (eye) | 4.53 GeV OA>9 | 4.53 GeV no OA\n");
    const double W[11][3] = {{10, 30, 18}, {30, 50, 42}, {50, 80, 65}, {80, 100, 88}, {100, 120, 112}, {120, 160, 140}, {160, 200, 180}, {200, 250, 225}, {250, 280, 260}, {280, 340, 320}, {340, 600, 450}};
    for (auto& w : W) { double nv = 0; for (size_t i = 0; i < nx.size(); ++i) if (std::fabs(nx[i] - w[2]) < 0.5) nv = ny[i]; K k9 = kOf(a[0], a[1], a[2], w[0] / 1000, w[1] / 1000), k0 = kOf(b[0], b[1], b[2], w[0] / 1000, w[1] / 1000); printf(" %4.0f-%4.0f  |   %.3f    |  %.3f |  %.3f\n", w[0], w[1], nv, k9.k, k0.k); }
    TCanvas c("c", "", 950, 700); c.SetLeftMargin(0.12); c.SetRightMargin(0.03); c.SetTopMargin(0.03); c.SetBottomMargin(0.13); c.SetTicks(1, 1);
    TH1D fr("fr", ";M_{ee} [MeV/c^{2}];k = N^{mix}_{+-} / 2#sqrt{N^{mix}_{++}N^{mix}_{--}}", 10, 0, 600); fr.SetMinimum(0.65); fr.SetMaximum(1.15); fr.GetXaxis()->SetTitleSize(0.05); fr.GetYaxis()->SetTitleSize(0.05); fr.GetXaxis()->SetLabelSize(0.045); fr.GetYaxis()->SetLabelSize(0.045); fr.GetYaxis()->SetTitleOffset(1.15); fr.Draw("AXIS");
    TLine l(0, 1, 600, 1); l.SetLineStyle(2); l.Draw();
    TGraphErrors* gn = new TGraphErrors(); for (size_t i = 0; i < nx.size(); ++i) { gn->SetPoint(i, nx[i], ny[i]); gn->SetPointError(i, 0, 0.01); }
    gn->SetMarkerStyle(24); gn->SetMarkerSize(1.3); gn->SetMarkerColor(kBlack); gn->SetLineColor(kBlack); gn->Draw("P SAME");
    g9->SetMarkerStyle(20); g9->SetMarkerColor(kRed + 1); g9->SetLineColor(kRed + 1); g9->SetMarkerSize(1.1); g9->Draw("P SAME");
    g0->SetMarkerStyle(21); g0->SetMarkerColor(kBlue + 1); g0->SetLineColor(kBlue + 1); g0->SetMarkerSize(1.0); g0->Draw("P SAME");
    TLegend lg(0.33, 0.18, 0.95, 0.40); lg.SetBorderSize(0); lg.SetFillStyle(0); lg.SetTextSize(0.036); lg.AddEntry(gn, "p+p 1.58 GeV, note Fig. 18 (read by eye), OA > 9#circ", "pe"); lg.AddEntry(g9, "p+p 4.53 GeV, data, OA > 9#circ", "pe"); lg.AddEntry(g0, "p+p 4.53 GeV, data, no OA cut", "pe"); lg.Draw();
    c.SaveAs(Form("%s.pdf", out)); c.SaveAs(Form("%s.png", out));
}
