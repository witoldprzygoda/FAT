// wp2_fig_cnt.C -- COPY of wp2_fig.C (figure agent B, 2026-11-15): "counts" instead of "rows" in three axis/annotation
// strings; writes ONLY docs/latex/figs_all/wp2_N1v2_cnt.{pdf,png}. Nothing else changed.
// Run: root -l -b -q -e 'gSystem->SetBuildDir("<scratch>/aclic_wp2_cnt", kTRUE)' 'plots/all_report/wp2_fig_cnt.C+'
// ============================================================================
// wp2_fig.C  --  WP2: Fig. N1 (normalisation items: sector-2 deficit and target-window contribution)
// ============================================================================
// Panels
//   (a) early/late track rate per 10 deg of azimuth (S0, all four tracks), normalised to the bins outside sector 2
//   (b) eps per time block of ~25 nb^-1 (S0 rows with a sector-2 track, double ratio against rows without one)
//   (c) vertex z of the S0 rows without the vz cut, days 034 + 047 + 065: full target, empty target scaled to the same
//       flux (normaliser window), their difference and the fitted LH2 simulation
//   (d) f_win per day and combined for N - CB in six M4 windows, with the normaliser range and the elastic w6 of L
// Inputs: wp2_sector.root, wp2_blocks_S0.csv, wp2_win.root, wp2_vzfit.root, docs/latex/tables_all/wp2_win.txt
// Output: docs/latex/figs_all/wp2_N1v2.pdf/.png (wp2_N1.pdf/.png = the first attempt, panels (a), (b) only) and
//         wp2_N1e.pdf/.png (M4 spectra of the full target and of the material, with f_win(M4))
// Run: root -l -b -q 'plots/all_report/wp2_fig.C+'  (ACLiC; build dir aclic_wp2)
// ============================================================================
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <cmath>

#include "TBox.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"

namespace wp2fig {
const char* kOut = "/hdd3/przygoda/HADES/FAT_out/all_report/wp2";
const char* kFig = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_all";
const char* kTab = "/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all";
void Style(TH1* h, double ls = 0.050, double ts = 0.055) {
    h->GetXaxis()->SetLabelSize(ls); h->GetYaxis()->SetLabelSize(ls);
    h->GetXaxis()->SetTitleSize(ts); h->GetYaxis()->SetTitleSize(ts);
    h->GetXaxis()->SetTitleOffset(1.0); h->GetYaxis()->SetTitleOffset(1.25);
}
void PanelA(TPad* p) {
    p->cd(); p->SetLeftMargin(0.15); p->SetBottomMargin(0.14); p->SetTopMargin(0.08); p->SetRightMargin(0.04);
    TFile* f = TFile::Open(Form("%s/wp2_sector.root", kOut));
    TH1D* h = (TH1D*)f->Get("hPhiRatio")->Clone("hPhiRatioC");
    h->SetDirectory(nullptr);
    f->Close();
    h->SetTitle(";#varphi [deg];track rate (034#minus045) / (046#minus067), normalised");
    // the profile of wp2_sector.root uses early = 034-044, late = 045-067
    h->GetYaxis()->SetTitle("early / late track rate (norm.)");
    h->SetMinimum(0.86); h->SetMaximum(1.14);
    h->SetMarkerStyle(20); h->SetMarkerSize(0.9); h->SetMarkerColor(kBlack); h->SetLineColor(kBlack);
    Style(h);
    h->Draw("E1");
    TBox* b = new TBox(180, 0.86, 240, 1.14); b->SetFillColorAlpha(kRed, 0.12); b->Draw();
    h->Draw("E1 same");
    TLine l; l.SetLineStyle(2); l.SetLineColor(kGray + 2);
    for (int k = 1; k < 6; ++k) l.DrawLine(60 * k, 0.86, 60 * k, 1.14);
    l.SetLineStyle(1); l.DrawLine(0, 1, 360, 1);
    TLatex t; t.SetTextSize(0.040); t.SetTextAlign(22);
    const int secOf[6] = {5, 0, 1, 2, 3, 4};
    for (int k = 0; k < 6; ++k) t.DrawLatex(60 * k + 30, 1.125, Form("S%d", secOf[k]));
    t.SetNDC(); t.SetTextAlign(11); t.SetTextSize(0.042);
    t.DrawLatex(0.17, 0.17, "S0 e^{+}e^{-} counts, tracks per nb^{-1}: days 034#minus045 / 046#minus067");
    t.SetTextSize(0.055); t.DrawLatex(0.17, 0.79, "(a)");
}
void PanelB(TPad* p) {
    p->cd(); p->SetLeftMargin(0.15); p->SetBottomMargin(0.14); p->SetTopMargin(0.08); p->SetRightMargin(0.04);
    std::ifstream in(Form("%s/wp2_blocks_S0.csv", kOut));
    std::string line;
    std::getline(in, line);
    std::vector<double> x, y, ey, ex;
    std::vector<int> day;
    double cum = 0;
    while (std::getline(in, line)) {
        std::stringstream ss(line);
        std::string tok;
        std::vector<std::string> t;
        while (std::getline(ss, tok, ',')) t.push_back(tok);
        const double L = std::stod(t[4]);
        x.push_back(cum + 0.5 * L); ex.push_back(0.5 * L);
        y.push_back(std::stod(t[5])); ey.push_back(std::stod(t[6]));
        day.push_back(std::stoi(t[1]));
        cum += L;
    }
    TH1D* fr = new TH1D("frB", ";integrated luminosity of the analysed runs [nb^{-1}];#varepsilon_{s2} (counts with a sector-2 track)", 10, 0, cum);
    fr->SetMinimum(0.45); fr->SetMaximum(1.15);
    Style(fr);
    fr->Draw("AXIS");
    // deficit span 034-045
    double xe = 0;
    for (size_t i = 0; i < x.size(); ++i) if (day[i] <= 45) xe = x[i] + ex[i];
    TBox* b = new TBox(0, 0.45, xe, 1.15); b->SetFillColorAlpha(kRed, 0.10); b->Draw();
    TGraphErrors* g = new TGraphErrors((int)x.size(), x.data(), y.data(), ex.data(), ey.data());
    g->SetMarkerStyle(20); g->SetMarkerSize(0.5); g->SetLineColor(kBlack);
    g->Draw("P same");
    TLine l; l.SetLineStyle(2); l.DrawLine(0, 1, cum, 1);
    // day labels at the first block of 034, 040, 045, 047, 053, 059, 065
    TLatex t; t.SetTextSize(0.036); t.SetTextAlign(21);
    int last = -1;
    for (size_t i = 0; i < x.size(); ++i) {
        if (day[i] != last && (day[i] == 34 || day[i] == 40 || day[i] == 45 || day[i] == 49 || day[i] == 59 || day[i] == 65)) {
            TLine dl; dl.SetLineColor(kGray); dl.SetLineStyle(3); dl.DrawLine(x[i] - ex[i], 0.45, x[i] - ex[i], 1.10);
            t.DrawLatex(x[i] - ex[i] + 110, 1.115, Form("%03d", day[i]));
        }
        last = day[i];
    }
    t.SetNDC(); t.SetTextAlign(11); t.SetTextSize(0.045);
    t.DrawLatex(0.17, 0.17, "blocks of 25 nb^{-1}; reference: days 047#minus067");
    t.SetTextSize(0.055); t.DrawLatex(0.17, 0.30, "(b)");
}
// (c) vertex z of the unlike-sign rows (base + both pions in the vertex fit, all M4), days 034 + 047 + 065: full target,
//     empty target x b (each day with its own flux ratio b = U_full/U_empty, channel rows in [-200,-170) mm), the
//     difference, and the LH2 simulation fitted to the full-target profile (wp2_vzfit.root)
void PanelC(TPad* p) {
    p->cd(); p->SetLeftMargin(0.15); p->SetBottomMargin(0.14); p->SetTopMargin(0.08); p->SetRightMargin(0.04);
    p->SetLogy();
    TFile* fw = TFile::Open(Form("%s/wp2_win.root", kOut));
    TFile* ff = TFile::Open(Form("%s/wp2_vzfit.root", kOut));
    const int days[3] = {34, 47, 65};
    TH1D *hF = nullptr, *hM = nullptr, *hL = nullptr;
    for (int d = 0; d < 3; ++d) {
        double U[2] = {0, 0};
        for (int t = 0; t < 2; ++t)
            for (const char* smp : {"us", "pp", "mm"}) {
                TH1D* h = (TH1D*)fw->Get(Form("vzbase_%s_%s_%03d", t ? "empty" : "full", smp, days[d]));
                U[t] += h->Integral(h->FindBin(-199.999), h->FindBin(-170.001));
            }
        const double b = U[0] / U[1];
        TH1D* f = (TH1D*)fw->Get(Form("vz_full_us_%03d", days[d]));
        TH1D* e = (TH1D*)fw->Get(Form("vz_empty_us_%03d", days[d]));
        TH1D* l = (TH1D*)ff->Get(Form("fit_lh2_%03d", days[d]));
        if (!hF) {
            hF = (TH1D*)f->Clone("pcF"); hF->SetDirectory(nullptr);
            hM = (TH1D*)e->Clone("pcM"); hM->SetDirectory(nullptr); hM->Scale(b);
            hL = (TH1D*)l->Clone("pcL"); hL->SetDirectory(nullptr);
        } else { hF->Add(f); hM->Add(e, b); hL->Add(l); }
    }
    hF->Rebin(2); hM->Rebin(2); hL->Rebin(2);
    TH1D* hD = (TH1D*)hF->Clone("pcD"); hD->SetDirectory(nullptr); hD->Add(hM, -1);
    hF->SetTitle(";v_{z} [mm];counts / 4 mm");
    hF->GetXaxis()->SetRangeUser(-260, -40);
    hF->SetMinimum(20); hF->SetMaximum(2e5);
    Style(hF);
    hF->SetMarkerStyle(20); hF->SetMarkerSize(0.7); hF->SetMarkerColor(kBlack); hF->SetLineColor(kBlack);
    hF->Draw("E1");
    TBox* nb = new TBox(-200, 20, -170, 2e5); nb->SetFillColorAlpha(kOrange, 0.25); nb->Draw();
    hM->SetFillColorAlpha(kRed, 0.35); hM->SetLineColor(kRed + 1);
    hM->Draw("HIST same");
    hL->SetLineColor(kGreen + 2); hL->SetLineWidth(2); hL->SetLineStyle(2);
    hL->Draw("HIST same");
    hD->SetMarkerStyle(24); hD->SetMarkerSize(0.6); hD->SetMarkerColor(kBlue + 1); hD->SetLineColor(kBlue + 1);
    hD->Draw("E1 same");
    hF->Draw("E1 same");
    TLine l; l.SetLineStyle(2); l.SetLineColor(kGray + 2); l.DrawLine(-160, 20, -160, 2e5);
    TLegend* lg = new TLegend(0.17, 0.62, 0.60, 0.90);
    lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.040);
    lg->AddEntry(hF, "full target (days 034, 047, 065)", "pe");
    lg->AddEntry(hM, "empty target #times b (same day)", "f");
    lg->AddEntry(hD, "full #minus empty #times b", "pe");
    lg->AddEntry(hL, "LH_{2} simulation (fit)", "l");
    lg->AddEntry(nb, "flux normaliser window", "f");
    lg->Draw();
    TLatex t; t.SetNDC(); t.SetTextSize(0.040); t.SetTextAlign(31);
    t.DrawLatex(0.95, 0.17, "unlike-sign, isBest #geq 0, PT3, cor_ok, #pi^{#pm} in vertex");
    t.SetTextAlign(11); t.SetTextSize(0.055); t.DrawLatex(0.88, 0.84, "(c)");
    fw->Close(); ff->Close();
}
// (d) f_win (N - CB, S0) per M4 window: the three days and their L-weighted combination; grey box = normaliser range
//     (vz-profile fit = lower, positivity bound = upper); dotted band = elastic window fraction inside L (w6)
void PanelD(TPad* p) {
    p->cd(); p->SetLeftMargin(0.15); p->SetBottomMargin(0.14); p->SetTopMargin(0.08); p->SetRightMargin(0.04);
    const int nW = 6;
    const char* wins[nW] = {"all M4", "B", "P", "B+P", "0.40-0.70", "CR-U"};
    const char* lab[nW] = {"all M_{4}", "bump B", "peak P", "B + P", "0.40#minus0.70", "CR-U"};
    double f[4][nW] = {}, ef[4][nW] = {};   // [day 0..2, 3 = combined]
    std::ifstream in(Form("%s/wp2_win.txt", kTab));
    std::string line;
    while (std::getline(in, line)) {
        if (line.size() < 40 || line.compare(0, 2, "  ") != 0) continue;
        const std::string sel = line.substr(2, 5), win = line.substr(8, 13), qn = line.substr(22, 5), day = line.substr(28, 3);
        auto trim = [](std::string x) { while (!x.empty() && x.back() == ' ') x.pop_back(); return x; };
        if (trim(sel) != "S0" || trim(qn) != "N-CB") continue;
        int iw = -1;
        for (int k = 0; k < nW; ++k) if (trim(win) == wins[k]) iw = k;
        if (iw < 0) continue;
        const int id = day == "034" ? 0 : day == "047" ? 1 : day == "065" ? 2 : day == "all" ? 3 : -1;
        if (id < 0) continue;
        const size_t bar = line.find('|');
        double a, e;
        if (sscanf(line.c_str() + bar + 1, "%lf +- %lf", &a, &e) == 2) { f[id][iw] = a; ef[id][iw] = e; }
    }
    TH1D* fr = new TH1D("frD", ";;f_{win} (N #minus CB, S0) [%]", nW, 0, nW);
    for (int k = 0; k < nW; ++k) fr->GetXaxis()->SetBinLabel(k + 1, lab[k]);
    fr->SetMinimum(0); fr->SetMaximum(36);
    Style(fr);
    fr->GetXaxis()->SetLabelSize(0.060);
    fr->Draw("AXIS");
    TBox* w6 = new TBox(0, 0.70, nW, 0.88); w6->SetFillColor(kGray + 1); w6->SetFillStyle(3002); w6->Draw();
    const double lo = 0.867, hi = 1.159;   // normaliser range (L-weighted b_fit/b_ch and b_max/b_ch)
    for (int k = 0; k < nW; ++k) {
        TBox* bx = new TBox(k + 0.62, f[3][k] * lo, k + 0.90, f[3][k] * hi);
        bx->SetFillColorAlpha(kGray, 0.6); bx->Draw();
    }
    const int col[4] = {kRed + 1, kOrange + 7, kAzure + 2, kBlack};
    const int mst[4] = {21, 22, 23, 20};
    const double off[4] = {0.18, 0.32, 0.46, 0.76};
    TLegend* lg = new TLegend(0.17, 0.70, 0.93, 0.90);
    lg->SetNColumns(2); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.040);
    const char* nm[4] = {"day 034 (period I)", "day 047 (IIA, IIB)", "day 065 (III)", "L-weighted, stat"};
    for (int d = 0; d < 4; ++d) {
        TGraphErrors* g = new TGraphErrors(nW);
        for (int k = 0; k < nW; ++k) { g->SetPoint(k, k + off[d], f[d][k]); g->SetPointError(k, 0, ef[d][k]); }
        g->SetMarkerStyle(mst[d]); g->SetMarkerColor(col[d]); g->SetLineColor(col[d]); g->SetMarkerSize(d == 3 ? 1.3 : 1.0);
        g->Draw("P same");
        lg->AddEntry(g, nm[d], "pe");
    }
    TBox* dummy = new TBox(); dummy->SetFillColorAlpha(kGray, 0.6);
    lg->AddEntry(dummy, "normaliser range", "f");
    lg->AddEntry(w6, "elastic w_{6} in L (0.70#minus0.88 %)", "f");
    lg->Draw();
    TLatex t; t.SetNDC(); t.SetTextSize(0.055); t.DrawLatex(0.88, 0.60, "(d)");
}
// (e) M4 spectra of S0 N - CB, days 034 + 047 + 065: full target (black) and material = sum_d b_d x empty_d (red);
//     lower pad: the ratio f_win(M4) (20 MeV bins)
void PanelE(TCanvas* c) {
    TFile* fw = TFile::Open(Form("%s/wp2_win.root", kOut));
    const int days[3] = {34, 47, 65};
    TH1D *F[3] = {}, *M[3] = {};   // [sample us, pp, mm]
    const char* smp[3] = {"us", "pp", "mm"};
    for (int d = 0; d < 3; ++d) {
        double U[2] = {0, 0};
        for (int t = 0; t < 2; ++t)
            for (int k = 0; k < 3; ++k) {
                TH1D* h = (TH1D*)fw->Get(Form("vzbase_%s_%s_%03d", t ? "empty" : "full", smp[k], days[d]));
                U[t] += h->Integral(h->FindBin(-199.999), h->FindBin(-170.001));
            }
        const double b = U[0] / U[1];
        for (int k = 0; k < 3; ++k) {
            TH1D* f = (TH1D*)fw->Get(Form("m4_full_%s_%03d", smp[k], days[d]));
            TH1D* e = (TH1D*)fw->Get(Form("m4_empty_%s_%03d", smp[k], days[d]));
            if (!F[k]) {
                F[k] = (TH1D*)f->Clone(Form("peF%d", k)); F[k]->SetDirectory(nullptr);
                M[k] = (TH1D*)e->Clone(Form("peM%d", k)); M[k]->SetDirectory(nullptr); M[k]->Scale(b);
            } else { F[k]->Add(f); M[k]->Add(e, b); }
        }
    }
    fw->Close();
    for (int k = 0; k < 3; ++k) { F[k]->Rebin(2); M[k]->Rebin(2); }
    // N - CB per bin, CB = 2 sqrt(N++ N--) (scaled counts for the material: CB of the scaled like-sign)
    auto nmcb = [](TH1D* us, TH1D* pp, TH1D* mm, const char* name) {
        TH1D* h = (TH1D*)us->Clone(name); h->SetDirectory(nullptr);
        for (int i = 1; i <= h->GetNbinsX(); ++i) {
            const double np = pp->GetBinContent(i), nm = mm->GetBinContent(i);
            const double cb = (np > 0 && nm > 0) ? 2 * std::sqrt(np * nm) : 0;
            const double vcb = (np > 0 && nm > 0) ? (nm / np) * std::pow(pp->GetBinError(i), 2) + (np / nm) * std::pow(mm->GetBinError(i), 2) : 0;
            h->SetBinContent(i, us->GetBinContent(i) - cb);
            h->SetBinError(i, std::sqrt(std::pow(us->GetBinError(i), 2) + vcb));
        }
        return h;
    };
    TH1D* hF = nmcb(F[0], F[1], F[2], "peFull");
    TH1D* hM = nmcb(M[0], M[1], M[2], "peMat");
    TH1D* hR = (TH1D*)hM->Clone("peRatio"); hR->SetDirectory(nullptr);
    hR->Divide(hF);
    c->cd();
    TPad* p1 = new TPad("pe1", "", 0, 0.32, 1, 1); p1->SetBottomMargin(0.02); p1->SetLeftMargin(0.13); p1->SetRightMargin(0.04); p1->SetLogy();
    TPad* p2 = new TPad("pe2", "", 0, 0, 1, 0.32); p2->SetTopMargin(0.02); p2->SetBottomMargin(0.32); p2->SetLeftMargin(0.13); p2->SetRightMargin(0.04);
    p1->Draw(); p2->Draw();
    p1->cd();
    hF->SetTitle(";M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}];N #minus CB / 20 MeV");
    hF->SetMinimum(5); hF->SetMaximum(3e5);
    hF->GetXaxis()->SetLabelSize(0); hF->GetYaxis()->SetLabelSize(0.05); hF->GetYaxis()->SetTitleSize(0.06); hF->GetYaxis()->SetTitleOffset(0.95);
    hF->SetMarkerStyle(20); hF->SetMarkerSize(0.8); hF->SetMarkerColor(kBlack); hF->SetLineColor(kBlack);
    hF->Draw("E1");
    hM->SetMarkerStyle(21); hM->SetMarkerSize(0.8); hM->SetMarkerColor(kRed + 1); hM->SetLineColor(kRed + 1);
    hM->Draw("E1 same");
    TBox* bB = new TBox(0.4550, 5, 0.5225, 3e5); bB->SetFillColorAlpha(kAzure, 0.10); bB->Draw();
    TBox* bP = new TBox(0.53364, 5, 0.55416, 3e5); bP->SetFillColorAlpha(kGreen + 2, 0.15); bP->Draw();
    hF->Draw("E1 same"); hM->Draw("E1 same");
    TLegend* lg = new TLegend(0.50, 0.66, 0.94, 0.92);
    lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.050);
    lg->AddEntry(hF, "full target, days 034, 047, 065", "pe");
    lg->AddEntry(hM, "material: #Sigma_{d} b_{d} #times empty target", "pe");
    lg->AddEntry(bB, "bump window B", "f");
    lg->AddEntry(bP, "peak window P", "f");
    lg->Draw();
    TLatex t; t.SetNDC(); t.SetTextSize(0.050);
    t.DrawLatex(0.16, 0.08, "S0, isBest #geq 0, unlike-sign minus 2#sqrt{N_{++}N_{--}}");
    p2->cd();
    hR->SetTitle(";M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}];f_{win}");
    hR->SetMinimum(0); hR->SetMaximum(0.45);
    hR->GetXaxis()->SetLabelSize(0.11); hR->GetXaxis()->SetTitleSize(0.12); hR->GetXaxis()->SetTitleOffset(1.1);
    hR->GetYaxis()->SetLabelSize(0.10); hR->GetYaxis()->SetTitleSize(0.12); hR->GetYaxis()->SetTitleOffset(0.45); hR->GetYaxis()->SetNdivisions(505);
    hR->SetMarkerStyle(20); hR->SetMarkerSize(0.8); hR->SetMarkerColor(kRed + 1); hR->SetLineColor(kRed + 1);
    hR->Draw("E1");
    TLine l; l.SetLineStyle(2); l.SetLineColor(kGray + 2);
    l.DrawLine(0.30, 0.1186, 1.10, 0.1186);
}
}  // namespace wp2fig

void wp2_fig_cnt() {
    using namespace wp2fig;
    gStyle->SetOptStat(0);
    gStyle->SetEndErrorSize(0);
    TCanvas* c = new TCanvas("cN1", "Fig. N1", 1400, 1100);
    c->Divide(2, 2);
    PanelA((TPad*)c->cd(1));
    PanelB((TPad*)c->cd(2));
    PanelC((TPad*)c->cd(3));
    PanelD((TPad*)c->cd(4));
    c->SaveAs(Form("%s/wp2_N1v2_cnt.pdf", kFig));
    c->SaveAs(Form("%s/wp2_N1v2_cnt.png", kFig));
    return;   // _cnt copy: Fig. N1e is not redrawn
    TCanvas* ce = new TCanvas("cN1e", "Fig. N1e", 900, 750);
    PanelE(ce);
    ce->SaveAs(Form("%s/wp2_N1e.pdf", kFig));
    ce->SaveAs(Form("%s/wp2_N1e.png", kFig));
}
