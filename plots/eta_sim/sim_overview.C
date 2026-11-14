// ============================================================================
// sim_overview.C  --  overview of the simulated M(pi+pi-e+e-) cross section
// ============================================================================
// Selection "pre" (isBest >= 0, PT3, vz > -160 mm), e+e- hypothesis, std + lep
// owned rows in partition mode A, weight w (PDG BR). Cross section within the
// HADES acceptance after the selection, no efficiency correction:
//   dsigma/dM = sum w / L_sim / bin width, L_sim = 221.05 nb^-1.
//
// Figures (plots/eta_sim/figs/):
//   engine_overview_pre.png/.pdf    per template: signal split into bump Dalitz,
//                                   bump conversion, peak conversion; physical
//                                   background; fake lepton and different
//                                   sources dashed (not part of the simulated
//                                   spectrum). (a) log 0.3-1.1, (b) linear zoom
//                                   0.40-0.70 with the signals stacked on the
//                                   physical background, (c) the signals alone
//   engine_data_vs_sim_pre.png/.pdf the same simulation with the data of the
//                                   same selection in the same units: N (black),
//                                   CB = 2 sqrt(N++ N--) (red), N - CB (blue),
//                                   L_data = 6080.5 nb^-1 (unconfirmed);
//                                   (a) log, (b) linear zoom, (c) ratio
//                                   (N - CB) / simulation per bin
// Table: plots/eta_sim/tables/engine_data_vs_sim_pre.txt (region integrals)
//
// Run from /home/przygoda/HADES/FAT/FAT:
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/aclic_engine", kTRUE)' 'plots/eta_sim/sim_overview.C+'
// ============================================================================
#include <cstdio>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TStopwatch.h"
#include "TStyle.h"

#include "SimEngine.h"

namespace sov {

void StyleH(TH1D* h, int col, int lw, int ls, int fill = 0) {
    h->SetLineColor(col); h->SetLineWidth(lw); h->SetLineStyle(ls);
    h->SetMarkerColor(col); h->SetMarkerSize(0);
    if (fill) { h->SetFillColor(fill); h->SetFillStyle(1001); }
    else h->SetFillStyle(0);
}
void Axes(TH1D* h, const char* ytitle, double scale = 1.0) {
    h->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}]");
    h->GetYaxis()->SetTitle(ytitle);
    h->GetXaxis()->SetTitleSize(0.055 * scale); h->GetYaxis()->SetTitleSize(0.055 * scale);
    h->GetXaxis()->SetLabelSize(0.048 * scale); h->GetYaxis()->SetLabelSize(0.048 * scale);
    h->GetXaxis()->SetTitleOffset(1.0); h->GetYaxis()->SetTitleOffset(1.25 / scale);
    h->GetYaxis()->SetMaxDigits(4);
}
void Txt(double x, double y, const char* s, double size, int col = 1) {
    TLatex t; t.SetNDC(); t.SetTextSize(size); t.SetTextColor(col); t.SetTextAlign(11); t.DrawLatex(x, y, s);
}
TPad* MkPad(const char* name, double x1, double y1, double x2, double y2, bool logy, double lm, double tm) {
    TPad* p = new TPad(name, "", x1, y1, x2, y2);
    p->SetLeftMargin(lm); p->SetRightMargin(0.02); p->SetTopMargin(tm); p->SetBottomMargin(0.13);
    p->SetLogy(logy);
    p->Draw();
    p->cd();
    return p;
}
const char* kYtit = "d#sigma/dM [nb/(GeV/c^{2})]";

}  // namespace sov

void sim_overview() {
    using namespace simeng;
    using namespace sov;
    TStopwatch swAll;
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetEndErrorSize(0);

    SimData d;
    if (!Load(d, kCompact, (1 << kStdUS) | (1 << kLepUS))) return;
    const Sel sel = SelPreset(kPre);
    SimOpt opt;   // mode A, w (PDG), std + lep, e+e-

    // components: 0 bump Dalitz, 1 bump conversion, 2 peak conversion, 3 physical background,
    //             4 fake lepton, 5 different sources, 6 signal + physical background
    std::vector<Comp> comps = {CompBumpDalitz(), CompBumpConv(), CompPeakConv(), CompPhysBkg(), CompFake(),
                               CompDiff(), CompSim()};
    const int col[7] = {kColBumpDalitz, kColBumpConv, kColPeakConv, kColPhysBkg, kColFake, kColDiff, kColSim};
    const char* leg[7] = {"signal: bump, #pi^{0} Dalitz", "signal: bump, conversion", "signal: peak, conversion",
                          "physical background", "fake lepton (not in the sum)", "different sources (not in the sum)",
                          "simulated spectrum = signal + physical background"};
    const char* shortName[3] = {"bump, #pi^{0} Dalitz", "bump, conversion", "peak, conversion"};
    Binning bLog; bLog.nbin = 160; bLog.lo = 0.30; bLog.hi = 1.10;     // 5 MeV
    Binning bLin; bLin.nbin = 120; bLin.lo = 0.40; bLin.hi = 0.70;     // 2.5 MeV
    TStopwatch swFill;
    const std::vector<WSpec> sLog = FillMany(d, sel, comps, opt, bLog);
    const std::vector<WSpec> sLin = FillMany(d, sel, comps, opt, bLin);
    printf("sim_overview: fills %.2f s\n", swFill.RealTime());
    TH1D *hLog[7], *hLin[7];
    for (int c = 0; c < 7; ++c) {
        hLog[c] = sLog[c].Hist(Form("hLog%d", c), 2);
        hLin[c] = sLin[c].Hist(Form("hLin%d", c), 2);
        printf("  %-40s sigma in 0.3-1.1: %.4g +- %.2g nb (all M4: %.4g nb), rows %.0f, N_eff %.0f\n", comps[c].name.c_str(),
               sLog[c].S / kLsim, std::sqrt(sLog[c].V) / kLsim, sLog[c].Sall / kLsim, sLog[c].N, sLog[c].NeffIn());
    }

    // ---- data (same selection, same binning, same units) ----
    scan::Data dd;
    if (!scan::Load(dd, scan::kCompact, 0.30, 1.10)) return;
    const DSpec dLog = DataSpectrumPreset(dd, kPre, bLog.nbin, bLog.lo, bLog.hi);
    const DSpec dLin = DataSpectrumPreset(dd, kPre, bLin.nbin, bLin.lo, bLin.hi);
    TH1D* hN = dLog.Hist("us", "dN", 2);
    TH1D* hC = dLog.Hist("cb", "dC", 2);
    TH1D* hS = dLog.Hist("sig", "dS", 2);
    TH1D* hS2 = dLin.Hist("sig", "dS2", 2);
    auto mk = [](TH1D* h, int c, int st) { h->SetMarkerStyle(st); h->SetMarkerSize(1.0); h->SetMarkerColor(c); h->SetLineColor(c); };
    mk(hN, 1, 20); mk(hC, 2, 21); mk(hS, 4, 20); mk(hS2, 4, 20);

    // y range of the log panels, shared by both figures (room for the legend on top)
    double ymax = 0, ymin = 1e30;
    for (int i = 1; i <= hLog[6]->GetNbinsX(); ++i) ymax = std::max(ymax, std::max(hLog[6]->GetBinContent(i), hN->GetBinContent(i)));
    for (int c = 0; c < 3; ++c)
        for (int i = 1; i <= hLog[c]->GetNbinsX(); ++i)
            if (hLog[c]->GetBinContent(i) > 0) ymin = std::min(ymin, hLog[c]->GetBinContent(i));
    const double yLo = std::max(ymin * 0.5, ymax * 1e-6), yHi = ymax * 400;

    const std::string simLine = Form("simulation: SMASH std + lep, owned rows (mode A), w (PDG BR); #sigma = #Sigma w / L_{sim}, L_{sim} = %.2f nb^{-1}", kLsim);
    const char* selLine = "selection \"pre\": isBest #geq 0, PT3, v_{z} > -160 mm, e^{+}e^{-}; within the HADES acceptance, no efficiency correction";
    const char* direct = "#eta #rightarrow #pi^{+}#pi^{-}e^{+}e^{-} (direct decay) is NOT in the simulation (eta set not available)";
    const std::string dataLine = Form("data: (N - CB) / L_{data} / bin width, %s; CB = 2#sqrt{N_{++}N_{--}} (10 MeV bins, no k factor); cor_ok = 1", kLdataNote);

    // ======================= figure 1: overview =======================
    {
        TCanvas* c = new TCanvas("cOv", "", 2000, 1700);
        // (a) log, full width
        c->cd(); MkPad("pA", 0.0, 0.5, 1.0, 1.0, true, 0.08, 0.24);
        TH1D* fr = (TH1D*)hLog[6]->Clone("frLog"); fr->Reset(); Axes(fr, kYtit);
        fr->GetYaxis()->SetTitleOffset(0.65);
        fr->SetMinimum(yLo); fr->SetMaximum(yHi); fr->Draw("axis");
        StyleH(hLog[3], kColPhysBkg, 1, 1, kColPhysFill); hLog[3]->Draw("hist same");
        StyleH(hLog[6], col[6], 4, 1); hLog[6]->Draw("hist same");
        for (int k = 0; k < 3; ++k) { StyleH(hLog[k], col[k], 3, 1); hLog[k]->Draw("hist same"); }
        StyleH(hLog[4], col[4], 3, 2); hLog[4]->Draw("hist same");
        StyleH(hLog[5], col[5], 3, 2); hLog[5]->Draw("hist same");
        fr->Draw("axis same");
        TLegend* lg = new TLegend(0.09, 0.57, 0.97, 0.75);
        lg->SetNColumns(2); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.040);
        const int ord[7] = {6, 0, 3, 1, 4, 2, 5};
        for (int k : ord)
            lg->AddEntry(hLog[k], Form("%s: %.3g nb", leg[k], sLog[k].S / kLsim), k == 3 ? "f" : "l");
        lg->Draw();
        Txt(0.08, 0.94, "(a) simulated d#sigma/dM_{#pi^{+}#pi^{-}e^{+}e^{-}} per template, 5 MeV bins (cross sections in the legend: 0.30-1.10 GeV/c^{2})", 0.045);
        Txt(0.08, 0.89, selLine, 0.036);
        Txt(0.08, 0.85, simLine.c_str(), 0.036);
        Txt(0.08, 0.805, direct, 0.036, 2);
        // (b) linear zoom, signals stacked on the physical background
        c->cd(); MkPad("pB", 0.0, 0.0, 0.5, 0.5, false, 0.16, 0.20);
        TH1D* hb = (TH1D*)hLin[3]->Clone("hbLin"); StyleH(hb, kColPhysBkg, 1, 1, kColPhysFill);
        TH1D* hs[3];
        const int fillc[3] = {kColBumpDalitz, kColBumpConv, kColPeakConv};
        for (int k = 0; k < 3; ++k) {
            hs[k] = (TH1D*)(k == 0 ? hb : hs[k - 1])->Clone(Form("hsLin%d", k));
            hs[k]->Add(hLin[k]);
            StyleH(hs[k], col[k], 1, 1, fillc[k]);
        }
        double ym = 0;
        for (int i = 1; i <= hLin[6]->GetNbinsX(); ++i) ym = std::max(ym, hLin[6]->GetBinContent(i) + hLin[6]->GetBinError(i));
        TH1D* fr2 = (TH1D*)hLin[6]->Clone("frLin"); fr2->Reset(); Axes(fr2, kYtit);
        fr2->SetMinimum(0); fr2->SetMaximum(ym * 1.05); fr2->Draw("axis");
        for (int k = 2; k >= 0; --k) hs[k]->Draw("hist same");
        hb->Draw("hist same");
        TH1D* hsum = (TH1D*)hLin[6]->Clone("hsumLin"); StyleH(hsum, col[6], 2, 1); hsum->Draw("E1 same");
        StyleH(hLin[4], col[4], 3, 2); hLin[4]->Draw("hist same");
        StyleH(hLin[5], col[5], 3, 2); hLin[5]->Draw("hist same");
        fr2->Draw("axis same");
        TLegend* lg2 = new TLegend(0.18, 0.47, 0.60, 0.79);
        lg2->SetBorderSize(0); lg2->SetFillStyle(0); lg2->SetTextSize(0.038);
        lg2->AddEntry(hsum, "signal + physical background", "le");
        lg2->AddEntry(hs[2], "+ peak, conversion", "f");
        lg2->AddEntry(hs[1], "+ bump, conversion", "f");
        lg2->AddEntry(hs[0], "+ bump, #pi^{0} Dalitz", "f");
        lg2->AddEntry(hb, "physical background", "f");
        lg2->AddEntry(hLin[4], "fake lepton (not stacked)", "l");
        lg2->AddEntry(hLin[5], "different sources (not stacked)", "l");
        lg2->Draw();
        Txt(0.16, 0.93, "(b) linear zoom, 2.5 MeV bins, signals stacked", 0.050);
        Txt(0.16, 0.875, "pre, e^{+}e^{-}; std + lep, mode A, w (PDG), L_{sim} = 221.05 nb^{-1}", 0.036);
        Txt(0.16, 0.825, "#eta #rightarrow #pi^{+}#pi^{-}e^{+}e^{-} (direct) NOT in the simulation", 0.038, 2);
        // (c) the signal channels alone
        c->cd(); MkPad("pC", 0.5, 0.0, 1.0, 0.5, false, 0.16, 0.20);
        double ys = 0;
        for (int k = 0; k < 3; ++k)
            for (int i = 1; i <= hLin[k]->GetNbinsX(); ++i) ys = std::max(ys, hLin[k]->GetBinContent(i) + hLin[k]->GetBinError(i));
        TH1D* fr3 = (TH1D*)hLin[0]->Clone("frSig"); fr3->Reset(); Axes(fr3, kYtit);
        fr3->SetMinimum(0); fr3->SetMaximum(ys * 1.45); fr3->Draw("axis");
        TLegend* lg3 = new TLegend(0.18, 0.62, 0.98, 0.79);
        lg3->SetBorderSize(0); lg3->SetFillStyle(0); lg3->SetTextSize(0.038);
        for (int k = 0; k < 3; ++k) {
            TH1D* h = (TH1D*)hLin[k]->Clone(Form("hsig%d", k)); StyleH(h, col[k], 3, 1);
            h->Draw("hist same");
            TH1D* he = (TH1D*)hLin[k]->Clone(Form("hsige%d", k)); StyleH(he, col[k], 1, 1); he->Draw("E1 same");
            lg3->AddEntry(h, Form("%s: %.3g #pm %.2g nb (N_{eff} %.0f)", shortName[k], sLin[k].S / kLsim,
                                  std::sqrt(sLin[k].V) / kLsim, sLin[k].NeffIn()), "l");
        }
        lg3->Draw();
        Txt(0.16, 0.93, "(c) signals alone (std set), 2.5 MeV bins", 0.050);
        Txt(0.16, 0.875, "pre, e^{+}e^{-}; #sigma in 0.40-0.70 GeV/c^{2}; errors from event groups", 0.038);
        Txt(0.16, 0.825, "#eta #rightarrow #pi^{+}#pi^{-}e^{+}e^{-} (direct) NOT in the simulation", 0.038, 2);
        c->SaveAs((std::string(kFigDir) + "/engine_overview_pre.png").c_str());
        c->SaveAs((std::string(kFigDir) + "/engine_overview_pre.pdf").c_str());
    }

    // ======================= figure 2: data in the same units =======================
    {
        TCanvas* c = new TCanvas("cDS", "", 2000, 1700);
        // (a) log, full width
        c->cd(); MkPad("dA", 0.0, 0.5, 1.0, 1.0, true, 0.08, 0.24);
        TH1D* fr = (TH1D*)hN->Clone("frD"); fr->Reset(); Axes(fr, kYtit);
        fr->GetYaxis()->SetTitleOffset(0.65);
        fr->SetMinimum(yLo); fr->SetMaximum(yHi); fr->Draw("axis");
        TH1D* b = (TH1D*)hLog[3]->Clone("bD"); StyleH(b, kColPhysBkg, 1, 1, kColPhysFill); b->Draw("hist same");
        TH1D* s = (TH1D*)hLog[6]->Clone("sD"); StyleH(s, col[6], 4, 1); s->Draw("hist same");
        TH1D* ks[3];
        for (int k = 0; k < 3; ++k) { ks[k] = (TH1D*)hLog[k]->Clone(Form("kD%d", k)); StyleH(ks[k], col[k], 3, 1); ks[k]->Draw("hist same"); }
        hN->Draw("E1 P same"); hC->Draw("E1 P same"); hS->Draw("E1 P same");
        fr->Draw("axis same");
        TLegend* lg = new TLegend(0.09, 0.57, 0.97, 0.75);
        lg->SetNColumns(2); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.040);
        lg->AddEntry(hN, "data: all e^{+}e^{-} combinations N", "p");
        lg->AddEntry(s, "sim: signal + physical background", "l");
        lg->AddEntry(hC, "data: CB = 2#sqrt{N_{++}N_{--}} (no k factor)", "p");
        lg->AddEntry(b, "sim: physical background", "f");
        lg->AddEntry(hS, "data: N - CB", "p");
        lg->AddEntry(ks[0], "sim: bump #pi^{0} Dalitz / bump conv. / peak conv.", "l");
        lg->Draw();
        Txt(0.08, 0.94, "(a) d#sigma/dM_{#pi^{+}#pi^{-}e^{+}e^{-}}, data and simulation in the same units, 5 MeV bins, simulation NOT renormalised", 0.045);
        Txt(0.08, 0.89, dataLine.c_str(), 0.036, 4);
        Txt(0.08, 0.85, (simLine + "; pre: isBest #geq 0, PT3, v_{z} > -160 mm").c_str(), 0.036);
        Txt(0.08, 0.805, direct, 0.036, 2);
        // (b) linear zoom: N - CB against the stacked simulation
        c->cd(); MkPad("dB", 0.0, 0.0, 0.5, 0.5, false, 0.16, 0.20);
        TH1D* hb = (TH1D*)hLin[3]->Clone("hbD"); StyleH(hb, kColPhysBkg, 1, 1, kColPhysFill);
        const int fillc[3] = {kColBumpDalitz, kColBumpConv, kColPeakConv};
        TH1D* hc[3];
        for (int k = 0; k < 3; ++k) {
            hc[k] = (TH1D*)(k == 0 ? hb : hc[k - 1])->Clone(Form("hsD%d", k));
            hc[k]->Add(hLin[k]);
            StyleH(hc[k], col[k], 1, 1, fillc[k]);
        }
        double ym = 0;
        for (int i = 1; i <= hS2->GetNbinsX(); ++i) ym = std::max(ym, hS2->GetBinContent(i) + hS2->GetBinError(i));
        for (int i = 1; i <= hLin[6]->GetNbinsX(); ++i) ym = std::max(ym, hLin[6]->GetBinContent(i));
        TH1D* fr2 = (TH1D*)hS2->Clone("frD2"); fr2->Reset(); Axes(fr2, kYtit);
        fr2->SetMinimum(0); fr2->SetMaximum(ym * 1.05); fr2->Draw("axis");
        for (int k = 2; k >= 0; --k) hc[k]->Draw("hist same");
        hb->Draw("hist same");
        hS2->Draw("E1 P same");
        fr2->Draw("axis same");
        TLegend* lg2 = new TLegend(0.18, 0.52, 0.60, 0.79);
        lg2->SetBorderSize(0); lg2->SetFillStyle(0); lg2->SetTextSize(0.038);
        lg2->AddEntry(hS2, "data: N - CB", "p");
        lg2->AddEntry(hc[2], "sim: + peak, conversion", "f");
        lg2->AddEntry(hc[1], "sim: + bump, conversion", "f");
        lg2->AddEntry(hc[0], "sim: + bump, #pi^{0} Dalitz", "f");
        lg2->AddEntry(hb, "sim: physical background", "f");
        lg2->Draw();
        Txt(0.16, 0.93, "(b) linear zoom, 2.5 MeV bins, sim stacked", 0.050);
        Txt(0.16, 0.875, "data: L_{data} = 6080.5 nb^{-1} (unconfirmed), CB without k", 0.036, 4);
        Txt(0.16, 0.825, "pre, e^{+}e^{-}; sim not renormalised, no direct #eta decay", 0.036);
        // (c) ratio (N - CB) / simulation, 5 MeV
        c->cd(); MkPad("dC", 0.5, 0.0, 1.0, 0.5, false, 0.16, 0.20);
        TH1D* hr = (TH1D*)hS->Clone("ratio");
        hr->Divide(hLog[6]);    // errors of both added in quadrature (independent samples)
        mk(hr, 4, 20);
        TH1D* fr3 = (TH1D*)hr->Clone("frR"); fr3->Reset(); Axes(fr3, "data (N - CB) / simulation");
        fr3->SetMinimum(0); fr3->SetMaximum(2.5); fr3->Draw("axis");
        TLine one(bLog.lo, 1, bLog.hi, 1); one.SetLineStyle(2); one.SetLineColor(kGray + 2); one.Draw();
        for (double xl : {scan::kBumpLo, scan::kBumpHi, scan::kPeakLo, scan::kPeakHi}) {
            TLine* l = new TLine(xl, 0, xl, 2.5); l->SetLineStyle(3); l->SetLineColor(kGray + 1); l->Draw();
        }
        hr->Draw("E1 P same");
        Txt(0.16, 0.93, "(c) data (N - CB) / sim (signal + phys. bkg)", 0.050);
        Txt(0.16, 0.875, "5 MeV bins; dotted: bump and peak windows of the data", 0.036);
        Txt(0.16, 0.825, "ratio scales with 1/L_{data} (unconfirmed); the peak-window ratio follows the", 0.036, 2);
        Txt(0.16, 0.780, "background trend of its neighbours: it does not measure the missing direct decay", 0.036, 2);
        c->SaveAs((std::string(kFigDir) + "/engine_data_vs_sim_pre.png").c_str());
        c->SaveAs((std::string(kFigDir) + "/engine_data_vs_sim_pre.pdf").c_str());
    }

    // ======================= region integrals =======================
    {
        const int nreg = 6;
        const double rlo[nreg] = {0.30, scan::kBumpLo, scan::kPeakLo, 0.58, 0.70, 0.30};
        const double rhi[nreg] = {scan::kBumpLo, scan::kBumpHi, scan::kPeakHi, 0.70, 1.10, 1.10};
        const char* rname[nreg] = {"below the bump", "bump window", "peak window", "above the peak", "high mass", "all 0.30-1.10"};
        // simulation: exact group variance (one sum per region)
        std::vector<WSum> sim(nreg);
        for (int r = 0; r < nreg; ++r) {
            const double a = rlo[r], bb = rhi[r];
            const Sel sr = SelAnd(sel, [a, bb](const Row& x) { return x.m4 >= a && x.m4 < bb; });
            sim[r] = SumMany(d, sr, std::vector<Comp>(1, CompSim()), opt)[0];
        }
        // data: windows of ScanEngine (per-group variance of N, CB variance per group)
        scan::FillOpt fo; fo.nbin = 320; fo.lo = 0.30; fo.hi = 1.10; fo.reb = 4;
        for (int r = 0; r < nreg; ++r) fo.windows.push_back(std::make_pair(rlo[r], rhi[r]));
        const scan::Spectra sp = scan::Fill(dd, DataSelPreset(kPre), fo);
        FILE* o = fopen((std::string(kTabDir) + "/engine_data_vs_sim_pre.txt").c_str(), "w");
        for (int pass = 0; pass < 2; ++pass) {
            FILE* f = pass ? o : stdout;
            if (!f) continue;
            fprintf(f, "# selection \"pre\" (e+e-): data (N - CB) / L_data against the simulation (signal + physical background)\n");
            fprintf(f, "# data: L_data = %.1f nb^-1 (UNCONFIRMED), CB = 2 sqrt(N++ N--) in 10 MeV bins, no k factor\n", kLdata);
            fprintf(f, "# sim: std + lep owned rows (mode A), w (PDG), L_sim = %.2f nb^-1; eta -> pi+pi-e+e- (direct) missing in sim\n", kLsim);
            fprintf(f, "# the peak-window ratio equals the background-shape trend around it (gap 0.5225-0.5336: 1.29, 0.554-0.58: 1.46, review finding 4);\n");
            fprintf(f, "# it does NOT measure the missing eta -> pi+pi-e+e-. Absolute ratios depend on the vz slice and the ghost treatment (tables/review_regions.txt)\n");
            fprintf(f, "%-16s %15s %10s %10s %18s %18s %16s\n", "region", "M4 [GeV/c2]", "N", "CB", "data sigma [nb]", "sim sigma [nb]", "data / sim");
            for (int r = 0; r < nreg; ++r) {
                const scan::WinSum& w = sp.win[2 + r];
                const double ds = w.sig() / kLdata, de = w.sigErr() / kLdata;
                const double ss = sim[r].Xsec(), se = sim[r].XsecErr();
                const double ra = ss > 0 ? ds / ss : 0, re = ss > 0 ? ra * std::sqrt(std::pow(de / ds, 2) + std::pow(se / ss, 2)) : 0;
                fprintf(f, "%-16s %7.4f-%-7.4f %10.0f %10.0f %9.4g +- %-6.2g %9.4g +- %-6.2g %7.3f +- %.3f\n", rname[r], rlo[r], rhi[r],
                        w.n[0], w.cb, ds, de, ss, se, ra, re);
            }
        }
        // position and width of the simulated peak (conversion) and of the bump edge region, M4 scale
        {
            double sw = 0, sx = 0, sx2 = 0, sr = 0, swr = 0;
            for (int set = 1; set <= 2; ++set)
                for (const Row& x : d.tree(set, kUS)) {
                    if (!Owned(x, opt) || x.tmpl != 1 || x.cls != 4 || !sel(x)) continue;
                    if (x.m4 > 0.52f && x.m4 < 0.58f) { sw += x.w; sx += x.w * x.m4; sx2 += x.w * x.m4 * x.m4; }
                    if (x.m4_true > 0) { swr += x.w; sr += x.w * x.m4 / x.m4_true; }
                }
            const double mean = sx / sw, rms = std::sqrt(std::max(0., sx2 / sw - mean * mean));
            for (int pass = 0; pass < 2; ++pass) {
                FILE* f = pass ? o : stdout;
                if (!f) continue;
                fprintf(f, "# sim peak (conversion, pre): weighted mean M4 in 0.52-0.58 = %.5f GeV/c2, rms %.5f (%.0f rows); mean m4/m4_true = %.5f\n",
                        mean, rms, sw / simtruth::fEtaPPG, sr / swr);   // class-4 rows carry w = fEtaPPG (rules v1: 0.849 hard-coded)
                fprintf(f, "# data peak (data study): 0.5439 GeV/c2, sigma 5.4 MeV; m_eta = 0.547862\n");
            }
        }
        if (o) fclose(o);
    }
    printf("sim_overview: total %.1f s\n", swAll.RealTime());
}
