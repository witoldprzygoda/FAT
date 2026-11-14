// ============================================================================
// sim_eta_overview.C  --  simulated M(pi+pi-e+e-) cross section per component, with the eta set
// ============================================================================
// e+e- hypothesis, std (owned, mode A) + lep (owned, mode A) + eta (owned)
// rows, weight w (PDG BR). Cross section within the HADES acceptance after
// the selection, no efficiency correction:
//   dsigma/dM = sum w / L_sim / bin width, L_sim = 221.05 nb^-1 per set.
// Selections "pre" (top row) and "slides" (bottom row).
//   column 1: log, 0.30-1.00 GeV/c2, 5 MeV bins, every component as a line,
//             the physical background of std + lep filled grey; the same y
//             range in both rows
//   column 2: linear zoom 0.40-0.70, 2.5 MeV bins, components stacked
//             (physical background std + lep, physical background eta, bump
//             pi0 Dalitz, bump conversion, peak conversion, peak direct) with
//             the sum and its errors; own y range per row (labelled)
//   column 3: linear zoom 0.40-0.70, the signals alone (bump Dalitz, bump
//             conversion, peak conversion, peak direct, peak total) and the
//             X17 shape for BR = 1e-5 (dashed, information only); the same y
//             range in both rows
// Fake lepton and different sources are not drawn (information only; see
// tables/eta_engine_reference.txt and figs/engine_overview_pre.png).
//
// Output: figs/eta_overview.png, .pdf; tables/eta_overview.txt (cross sections
// per component in 0.30-1.00 and in the windows, for both selections)
// Run from /home/przygoda/HADES/FAT/FAT:
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/aclic_eta", kTRUE)' 'plots/eta_sim/sim_eta_overview.C+'
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

namespace seov {

void StyleH(TH1D* h, int col, int lw, int ls, int fill = 0) {
    h->SetLineColor(col); h->SetLineWidth(lw); h->SetLineStyle(ls);
    h->SetMarkerColor(col); h->SetMarkerSize(0);
    if (fill) { h->SetFillColor(fill); h->SetFillStyle(1001); }
    else h->SetFillStyle(0);
}
void Axes(TH1D* h, double ts = 0.055) {
    h->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}]");
    h->GetYaxis()->SetTitle("d#sigma/dM [nb/(GeV/c^{2})]");
    h->GetXaxis()->SetTitleSize(ts); h->GetYaxis()->SetTitleSize(ts);
    h->GetXaxis()->SetLabelSize(ts * 0.87); h->GetYaxis()->SetLabelSize(ts * 0.87);
    h->GetXaxis()->SetTitleOffset(1.05); h->GetYaxis()->SetTitleOffset(1.35);
    h->GetYaxis()->SetMaxDigits(4);
    h->GetXaxis()->SetNdivisions(506);
}
void Txt(double x, double y, const char* s, double size, int col = 1) {
    TLatex t; t.SetNDC(); t.SetTextSize(size); t.SetTextColor(col); t.SetTextAlign(11); t.DrawLatex(x, y, s);
}
TPad* MkPad(const char* name, double x1, double y1, double x2, double y2, bool logy) {
    TPad* p = new TPad(name, "", x1, y1, x2, y2);
    p->SetLeftMargin(0.15); p->SetRightMargin(0.03); p->SetTopMargin(0.25); p->SetBottomMargin(0.13);
    p->SetLogy(logy);
    p->Draw();
    p->cd();
    return p;
}
void Windows(double ylo, double yhi) {
    for (double x : {scan::kBumpLo, scan::kBumpHi, scan::kPeakLo, scan::kPeakHi}) {
        TLine* l = new TLine(x, ylo, x, yhi); l->SetLineStyle(3); l->SetLineColor(kGray + 1); l->Draw();
    }
}

}  // namespace seov

void sim_eta_overview() {
    using namespace simeng;
    using namespace seov;
    TStopwatch swAll;
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetEndErrorSize(0);

    SimData d;
    if (!Load(d, kCompact, kMaskUS)) return;
    SimOpt opt;   // std + lep + eta, mode A, w (PDG), e+e-

    // components
    enum { cBD = 0, cBC, cPC, cPD, cPT, cX17, cBkgSL, cBkgEta, cSum, cN };
    std::vector<Comp> comps = {CompBumpDalitz(), CompBumpConv(), CompPeakConv(), CompPeakDirect(), CompPeakTotal(),
                               CompX17(1e-5), CompPhysBkgSets(kSetStd | kSetLep), CompEtaBkg(), CompSim()};
    const int col[cN] = {kColBumpDalitz, kColBumpConv, kColPeakConv, kColPeakDirect, kColPeakTotal, kColX17,
                         kColPhysBkg, kColEtaBkg, kColSim};
    const char* leg[cN] = {"bump: #eta#rightarrow#pi^{+}#pi^{-}#pi^{0}, #pi^{0} Dalitz",
                           "bump: #eta#rightarrow#pi^{+}#pi^{-}#pi^{0}, conversion",
                           "peak: #eta#rightarrow#pi^{+}#pi^{-}#gamma, conversion",
                           "peak: #eta#rightarrow#pi^{+}#pi^{-}e^{+}e^{-} direct (eta set)",
                           "peak total (conversion + direct)",
                           "X17, BR = 10^{-5} (illustration, not in the sum)",
                           "physical background, std + lep",
                           "physical background, eta set",
                           "simulated spectrum (signal + physical background)"};
    const int presets[2] = {kPre, kSlides};
    Binning bLog; bLog.nbin = 140; bLog.lo = 0.30; bLog.hi = 1.00;     // 5 MeV
    Binning bLin; bLin.nbin = 120; bLin.lo = 0.40; bLin.hi = 0.70;     // 2.5 MeV
    std::vector<WSpec> sLog[2], sLin[2];
    TH1D *hLog[2][cN], *hLin[2][cN];
    TStopwatch swFill;
    for (int ip = 0; ip < 2; ++ip) {
        const Sel sel = SelPreset(presets[ip]);
        sLog[ip] = FillMany(d, sel, comps, opt, bLog);
        sLin[ip] = FillMany(d, sel, comps, opt, bLin);
        for (int c = 0; c < cN; ++c) {
            hLog[ip][c] = sLog[ip][c].Hist(Form("hLog%d_%d", ip, c), 2);
            hLin[ip][c] = sLin[ip][c].Hist(Form("hLin%d_%d", ip, c), 2);
        }
    }
    printf("sim_eta_overview: fills %.2f s\n", swFill.RealTime());

    // ---- table: cross sections per component ----
    {
        FILE* o = fopen((std::string(kTabDir) + "/eta_overview.txt").c_str(), "w");
        for (int pass = 0; pass < 2; ++pass) {
            FILE* f = pass ? o : stdout;
            if (!f) continue;
            fprintf(f, "# eta_overview: simulated sigma [nb] per component (e+e-), std + lep + eta owned rows, mode A, w (PDG)\n");
            fprintf(f, "# sigma = sum w / L_sim, L_sim = %.2f nb^-1; within the acceptance, no efficiency correction; errors from event groups\n", kLsim);
            fprintf(f, "# windows from 2.5 MeV bins whose centre is inside (exact window sums: tables/eta_engine_reference.txt)\n");
            fprintf(f, "# X17: BR(eta -> pi+pi- X17, X17 -> e+e-) = 1e-5, illustration only, not in the sum\n");
            for (int ip = 0; ip < 2; ++ip) {
                fprintf(f, "\n## selection \"%s\"\n", PresetName(presets[ip]));
                fprintf(f, "%-62s %22s %22s %22s %10s\n", "component", "0.30-1.00", "bump 0.455-0.5225", "peak 0.5336-0.5542", "N_eff");
                for (int c = 0; c < cN; ++c) {
                    const WSpec& L = sLog[ip][c];
                    const auto b = sLin[ip][c].Window(scan::kBumpLo, scan::kBumpHi);
                    const auto p = sLin[ip][c].Window(scan::kPeakLo, scan::kPeakHi);
                    fprintf(f, "%-62s %10.4g +- %-9.2g %10.4g +- %-9.2g %10.4g +- %-9.2g %10.0f\n", comps[c].name.c_str(), L.S / kLsim,
                            std::sqrt(L.V) / kLsim, b.first / kLsim, std::sqrt(b.second) / kLsim, p.first / kLsim,
                            std::sqrt(p.second) / kLsim, L.NeffIn());
                }
            }
        }
        if (o) fclose(o);
    }

    // ---- y ranges ----
    double yMax = 0, yMin = 1e30;
    for (int ip = 0; ip < 2; ++ip)
        for (int c = 0; c < cN; ++c)
            for (int i = 1; i <= hLog[ip][c]->GetNbinsX(); ++i) {
                const double v = hLog[ip][c]->GetBinContent(i);
                if (v > yMax) yMax = v;
                if (v > 0 && v < yMin) yMin = v;
            }
    const double yLo = std::max(yMin * 0.5, yMax * 1e-6), yHi = yMax * 400;
    double ySig = 0;
    for (int ip = 0; ip < 2; ++ip)
        for (int c : {cBD, cBC, cPC, cPD, cPT, cX17})
            for (int i = 1; i <= hLin[ip][c]->GetNbinsX(); ++i)
                ySig = std::max(ySig, hLin[ip][c]->GetBinContent(i) + hLin[ip][c]->GetBinError(i));
    const double ySigTop = ySig * 2.3;

    const char* presetLine[2] = {"pre: isBest #geq 0, PT3, v_{z} > -160 mm, corrected momenta, e^{+}e^{-}",
                                 "slides: pre + #alpha_{lab} < 50#circ, M_{#pi#pi} < 0.46 GeV/c^{2}, #alpha_{#eta} > 140#circ, OA > 2#circ"};

    TCanvas* cv = new TCanvas("cEtaOv", "", 2700, 2250);
    // ---------------- header strip: quantity, normalisation, shared legend ----------------
    cv->cd();
    {
        TPad* ph = new TPad("pHead", "", 0.0, 0.855, 1.0, 1.0);
        ph->SetMargin(0, 0, 0, 0);
        ph->Draw(); ph->cd();
        Txt(0.012, 0.80, "Simulated d#sigma/dM_{#pi^{+}#pi^{-}e^{+}e^{-}} per component, with the SMASH eta set (direct decay #eta#rightarrow#pi^{+}#pi^{-}e^{+}e^{-} included)", 0.150);
        Txt(0.012, 0.60, Form("SMASH std + lep + eta, owned rows (mode A), w (PDG BR); d#sigma/dM = #Sigma w / L_{sim} / bin width, L_{sim} = %.2f nb^{-1} per set; in the acceptance, no efficiency corr.", kLsim), 0.105);
        Txt(0.012, 0.44, "#eta#rightarrow#pi^{+}#pi^{-}e^{+}e^{-}: BR = 2.68#times10^{-4} (PDG), W = 2 BR per row (q = 1/2); fake lepton and different sources not drawn (information only, tables/eta_engine_reference.txt)", 0.105);
        Txt(0.012, 0.28, NoteX17(1e-5).c_str(), 0.105, kColX17);
        // dummy histograms for the shared legend (line and fill of each component)
        TLegend* lg = new TLegend(0.005, 0.0, 0.995, 0.25);
        lg->SetNColumns(5); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.095);
        auto dummy = [&](const char* n, int lc, int lw, int ls, int fc) {
            TH1D* h = new TH1D(n, "", 1, 0, 1); h->SetDirectory(nullptr);
            h->SetLineColor(lc); h->SetLineWidth(lw); h->SetLineStyle(ls);
            if (fc) { h->SetFillColor(fc); h->SetFillStyle(1001); } else h->SetFillStyle(0);
            return h;
        };
        lg->AddEntry(dummy("lgSum", col[cSum], 4, 1, 0), "sum: signal + phys. background", "l");
        lg->AddEntry(dummy("lgB", kColPhysBkg, 1, 1, kColPhysFill), "phys. background, std + lep", "f");
        lg->AddEntry(dummy("lgE", col[cBkgEta], 3, 1, kColEtaFill), "phys. background, eta set", "lf");
        lg->AddEntry(dummy("lgBD", col[cBD], 3, 1, col[cBD]), "bump: #pi^{0} Dalitz", "lf");
        lg->AddEntry(dummy("lgBC", col[cBC], 3, 1, col[cBC]), "bump: conversion", "lf");
        lg->AddEntry(dummy("lgPC", col[cPC], 3, 1, col[cPC]), "peak: conversion (#pi^{+}#pi^{-}#gamma)", "lf");
        lg->AddEntry(dummy("lgPD", col[cPD], 3, 1, col[cPD]), "peak: direct #pi^{+}#pi^{-}e^{+}e^{-} (eta set; gen. M_{ee} not QED)", "lf");
        lg->AddEntry(dummy("lgPT", col[cPT], 2, 7, 0), "peak total (conv. + direct)", "l");
        lg->AddEntry(dummy("lgX", col[cX17], 3, 2, 0), "X17, BR 10^{-5} (not in the sum)", "l");
        lg->AddEntry(dummy("lgW", kGray + 1, 1, 3, 0), "bump / peak windows", "l");
        lg->Draw();
    }
    for (int ip = 0; ip < 2; ++ip) {
        const double y1 = ip == 0 ? 0.4275 : 0.0, y2 = ip == 0 ? 0.855 : 0.4275;
        const char* tag = ip == 0 ? "pre" : "slides";
        // ---------------- column 1: log ----------------
        cv->cd();
        MkPad(Form("pL%d", ip), 0.0, y1, 0.34, y2, true);
        gPad->SetTopMargin(0.13);
        TH1D* fr = (TH1D*)hLog[ip][cSum]->Clone(Form("frL%d", ip)); fr->Reset(); Axes(fr, 0.050);
        fr->SetMinimum(yLo); fr->SetMaximum(yHi); fr->Draw("axis");
        TH1D* hb = (TH1D*)hLog[ip][cBkgSL]->Clone(Form("bL%d", ip)); StyleH(hb, kColPhysBkg, 1, 1, kColPhysFill); hb->Draw("hist same");
        StyleH(hLog[ip][cSum], col[cSum], 4, 1); hLog[ip][cSum]->Draw("hist same");
        StyleH(hLog[ip][cBkgEta], col[cBkgEta], 3, 1); hLog[ip][cBkgEta]->Draw("hist same");
        for (int c : {cBD, cBC, cPC, cPD}) { StyleH(hLog[ip][c], col[c], 3, 1); hLog[ip][c]->Draw("hist same"); }
        StyleH(hLog[ip][cX17], col[cX17], 3, 2); hLog[ip][cX17]->Draw("hist same");
        Windows(yLo, yHi);
        fr->Draw("axis same");
        Txt(0.15, 0.935, Form("(%s) %s: log scale, 5 MeV bins", ip == 0 ? "a" : "d", tag), 0.046);
        Txt(0.15, 0.893, presetLine[ip], ip == 0 ? 0.030 : 0.027);
        Txt(0.18, 0.82, Form("#sigma(0.30-1.00): sum %.4g nb, direct %.3g nb", sLog[ip][cSum].S / kLsim, sLog[ip][cPD].S / kLsim), 0.034);
        Txt(0.18, 0.775, Form("X17 (BR 10^{-5}) %.3g nb, eta-set bkg %.3g nb", sLog[ip][cX17].S / kLsim, sLog[ip][cBkgEta].S / kLsim), 0.034);
        Txt(0.18, 0.73, "same y range in (a) and (d)", 0.034);
        // ---------------- column 2: linear zoom, stacked ----------------
        cv->cd();
        MkPad(Form("pZ%d", ip), 0.34, y1, 0.67, y2, false);
        gPad->SetTopMargin(0.13);
        const int stackOrder[6] = {cBkgSL, cBkgEta, cBD, cBC, cPC, cPD};
        const int fillc[6] = {kColPhysFill, kColEtaFill, kColBumpDalitz, kColBumpConv, kColPeakConv, kColPeakDirect};
        TH1D* hs[6];
        for (int k = 0; k < 6; ++k) {
            hs[k] = (TH1D*)hLin[ip][stackOrder[k]]->Clone(Form("hs%d_%d", ip, k));
            if (k > 0) hs[k]->Add(hs[k - 1]);
            StyleH(hs[k], k == 0 ? kColPhysBkg : col[stackOrder[k]], 1, 1, fillc[k]);
        }
        double ym = 0;
        for (int i = 1; i <= hLin[ip][cSum]->GetNbinsX(); ++i)
            ym = std::max(ym, hLin[ip][cSum]->GetBinContent(i) + hLin[ip][cSum]->GetBinError(i));
        TH1D* fz = (TH1D*)hLin[ip][cSum]->Clone(Form("frZ%d", ip)); fz->Reset(); Axes(fz, 0.050);
        fz->SetMinimum(0); fz->SetMaximum(ym * 1.30); fz->Draw("axis");
        for (int k = 5; k >= 0; --k) hs[k]->Draw("hist same");
        TH1D* hsum = (TH1D*)hLin[ip][cSum]->Clone(Form("hsum%d", ip)); StyleH(hsum, col[cSum], 2, 1); hsum->Draw("E1 same");
        TH1D* hx = (TH1D*)hLin[ip][cX17]->Clone(Form("hx%d", ip)); StyleH(hx, col[cX17], 3, 2); hx->Draw("hist same");
        Windows(0, ym * 1.30);
        fz->Draw("axis same");
        Txt(0.15, 0.935, Form("(%s) %s: linear, 2.5 MeV bins, stacked", ip == 0 ? "b" : "e", tag), 0.044);
        Txt(0.15, 0.893, "stack: bkg std+lep, bkg eta, bump Dalitz, bump conv., peak conv., peak direct", 0.0225);
        Txt(0.18, 0.82, "error bars: event-group variance of the sum", 0.034);
        Txt(0.18, 0.775, "y range of this panel only (the selections differ a lot)", 0.034);
        Txt(0.18, 0.73, "X17 (BR 10^{-5}) dashed from zero, not stacked", 0.034, kColX17);
        // ---------------- column 3: signals alone ----------------
        cv->cd();
        MkPad(Form("pS%d", ip), 0.67, y1, 1.0, y2, false);
        gPad->SetTopMargin(0.13);
        TH1D* fs = (TH1D*)hLin[ip][cBD]->Clone(Form("frS%d", ip)); fs->Reset(); Axes(fs, 0.050);
        fs->SetMinimum(0); fs->SetMaximum(ySigTop); fs->Draw("axis");
        TLegend* ls = new TLegend(0.17, 0.535, 0.97, 0.845);
        ls->SetBorderSize(0); ls->SetFillStyle(0); ls->SetTextSize(0.034);
        const char* sname[cN] = {"bump, #pi^{0} Dalitz", "bump, conversion", "peak, conversion", "peak, direct (eta set)",
                                 "peak total", "X17, BR 10^{-5} (illustr.)", "", "", ""};
        for (int c : {cBD, cBC, cPC, cPD, cPT, cX17}) {
            TH1D* h = (TH1D*)hLin[ip][c]->Clone(Form("hsg%d_%d", ip, c));
            StyleH(h, col[c], c == cPT ? 2 : 3, c == cX17 ? 2 : (c == cPT ? 7 : 1));
            h->Draw("hist same");
            if (c != cPT && c != cX17) {
                TH1D* he = (TH1D*)hLin[ip][c]->Clone(Form("hsge%d_%d", ip, c)); StyleH(he, col[c], 1, 1); he->Draw("E1 same");
            }
            const auto w = sLin[ip][c].Window(0.40, 0.70);
            ls->AddEntry(h, Form("%s: %.3g #pm %.1g nb", sname[c], w.first / kLsim, std::sqrt(w.second) / kLsim), "l");
        }
        Windows(0, ySig * 1.05);
        fs->Draw("axis same");
        ls->Draw();
        Txt(0.15, 0.935, Form("(%s) %s: signals alone, 2.5 MeV bins", ip == 0 ? "c" : "f", tag), 0.046);
        Txt(0.15, 0.893, "#sigma in 0.40-0.70 GeV/c^{2}; same y range in (c) and (f)", 0.030);
    }
    cv->SaveAs((std::string(kFigDir) + "/eta_overview.png").c_str());
    cv->SaveAs((std::string(kFigDir) + "/eta_overview.pdf").c_str());
    printf("sim_eta_overview: total %.1f s\n", swAll.RealTime());
}
