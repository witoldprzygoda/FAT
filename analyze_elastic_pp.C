////////////////////////////////////////////////////////////////////////////////
// analyze_elastic_pp.C
//
// Two independent fit series for elastic pp scattering:
//   Series A: slice Y around 0.54, project X -> fit dphi peak
//   Series B: slice X around 180,  project Y -> fit tantan peak
//
// For each slice: find actual histogram maximum, then find where it drops
// by 6% on each side -> that defines the allowed range for Gauss mean.
// Polynomial: only degree 2 or 3 (higher eats signal).
// Sigma: hard bounds (dphi: 0.3-3 deg, tantan: 0.005-0.05).
//
// Control: 8 canvases (2 directions x 4 widths), slices from center outward.
//
// Usage: root -l 'analyze_elastic_pp.C("output_pp.root")'
////////////////////////////////////////////////////////////////////////////////

#include <TFile.h>
#include <TH2.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCutG.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TMath.h>
#include <TLatex.h>
#include <TGraph.h>
#include <TLine.h>
#include <TPad.h>
#include <TROOT.h>
#include <TMarker.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>

struct FitResult {
    double sliceCenter, sliceLo, sliceHi;
    int    sliceWidth;
    double amp, mean, sigma;
    double ampErr, meanErr, sigmaErr;
    double chi2ndf;
    int    polyDeg;
    bool   valid;
    std::vector<double> pars;

    FitResult() : valid(false), chi2ndf(1e9), polyDeg(-1),
                  amp(0), mean(0), sigma(0), ampErr(0), meanErr(0), sigmaErr(0),
                  sliceCenter(0), sliceLo(0), sliceHi(0), sliceWidth(0) {}
};

// ============================================================================
// Find peak and 6%-drop range in a histogram near expectedPeak
// Returns: peakPos, rangeLo, rangeHi
// ============================================================================
void findPeakRange(TH1D* h, double expectedPeak, double searchWindow,
                   double& peakPos, double& rangeLo, double& rangeHi)
{
    // Find maximum bin within search window around expectedPeak
    int bLo = h->FindBin(expectedPeak - searchWindow);
    int bHi = h->FindBin(expectedPeak + searchWindow);
    double maxVal = -1;
    int maxBin = -1;
    for (int b = bLo; b <= bHi; b++) {
        double c = h->GetBinContent(b);
        if (c > maxVal) { maxVal = c; maxBin = b; }
    }
    peakPos = h->GetBinCenter(maxBin);

    // Find where content drops to 94% of max (6% drop) on each side
    double threshold = maxVal * 0.94;

    // Scan left
    rangeLo = peakPos;
    for (int b = maxBin; b >= bLo; b--) {
        if (h->GetBinContent(b) < threshold) {
            rangeLo = h->GetBinCenter(b);
            break;
        }
    }

    // Scan right
    rangeHi = peakPos;
    for (int b = maxBin; b <= bHi; b++) {
        if (h->GetBinContent(b) < threshold) {
            rangeHi = h->GetBinCenter(b);
            break;
        }
    }

    // Safety: ensure minimum range (at least 5 bins or reasonable physics range)
    double minHalfRange = std::max(h->GetBinWidth(1) * 5, searchWindow * 0.05);
    if (rangeHi - rangeLo < 2 * minHalfRange) {
        rangeLo = peakPos - minHalfRange;
        rangeHi = peakPos + minHalfRange;
    }
}

// ============================================================================
// Fit one slice: Gauss + pol2 or pol3
// Mean range determined dynamically from peak finding
// ============================================================================
FitResult fitSlice(TH1D* h,
                   double fitMin, double fitMax,
                   double expectedPeak, double searchWindow,
                   double sigGuess, double sigMin, double sigMax,
                   double ampGuess)
{
    static int debugCount = 0;  // print first 10 slice details
    bool debug = (debugCount < 10);

    FitResult best;
    if (!h || h->GetEntries() < 50) return best;

    // Dynamic peak finding
    double peakPos, meanLo, meanHi;
    findPeakRange(h, expectedPeak, searchWindow, peakPos, meanLo, meanHi);

    if (ampGuess <= 0) {
        ampGuess = h->GetBinContent(h->FindBin(peakPos));
        if (ampGuess <= 0) ampGuess = h->GetMaximum();
    }

    // Background from edges
    double bgL = 0, bgR = 0;
    int cL = 0, cR = 0;
    for (int b = h->FindBin(fitMin); b <= h->FindBin(fitMin + (fitMax-fitMin)*0.1); b++) { bgL += h->GetBinContent(b); cL++; }
    for (int b = h->FindBin(fitMax - (fitMax-fitMin)*0.1); b <= h->FindBin(fitMax); b++) { bgR += h->GetBinContent(b); cR++; }
    if (cL > 0) bgL /= cL;
    if (cR > 0) bgR /= cR;
    double bgEst = (bgL + bgR) / 2.0;

    // Clamp initial guesses to valid ranges
    double ampInit = std::max(2.0, ampGuess - bgEst);
    double sigInit = std::max(sigMin, std::min(sigMax, sigGuess));

    double bestChi2 = 1e9;

    for (int deg = 2; deg <= 3; deg++) {
        TString fn = Form("ftmp_%s_p%d", h->GetName(), deg);
        TF1* f = new TF1(fn, Form("gaus(0)+pol%d(3)", deg), fitMin, fitMax);
        f->SetNpx(500);

        // Set limits FIRST
        f->SetParLimits(0, 1.0, ampGuess * 20);
        f->SetParLimits(1, meanLo, meanHi);
        f->SetParLimits(2, sigMin, sigMax);

        // Then set initial values clamped within limits
        f->SetParameter(0, std::min(ampInit, ampGuess * 20));
        f->SetParameter(1, std::max(meanLo, std::min(meanHi, peakPos)));
        f->SetParameter(2, sigInit);

        f->SetParameter(3, bgEst);
        for (int i = 1; i <= deg; i++) f->SetParameter(3+i, 0);

        h->Fit(f, "RQN", "", fitMin, fitMax);
        int st = h->Fit(f, "RQN", "", fitMin, fitMax);

        double chi2 = (f->GetNDF() > 0) ? f->GetChisquare()/f->GetNDF() : 1e9;
        double fA = f->GetParameter(0), fS = fabs(f->GetParameter(2));

        // Accept fit if Gauss component is meaningful
        // chi2 cut relaxed: high-stats bins can give large chi2/ndf even for good fits
        bool ok = chi2 > 0 && chi2 < 1000 &&
                  fA > 1.0 && fS >= sigMin && fS <= sigMax;

        if (debug) {
            std::cout << Form("      [DBG] %s pol%d: st=%d chi2/ndf=%.1f amp=%.1f sig=%.4f mean=%.4f meanRange=[%.4f,%.4f] -> %s",
                h->GetName(), deg, st, chi2, fA, fS, f->GetParameter(1), meanLo, meanHi, ok?"OK":"REJECT") << std::endl;
        }

        if (ok && chi2 < bestChi2) {
            bestChi2       = chi2;
            best.amp       = fA;
            best.mean      = f->GetParameter(1);
            best.sigma     = fS;
            best.ampErr    = f->GetParError(0);
            best.meanErr   = f->GetParError(1);
            best.sigmaErr  = f->GetParError(2);
            best.chi2ndf   = chi2;
            best.polyDeg   = deg;
            best.valid     = true;
            best.pars.clear();
            for (int i = 0; i < f->GetNpar(); i++) best.pars.push_back(f->GetParameter(i));
        }
        delete f;
    }
    if (debug) debugCount++;
    return best;
}

// ============================================================================
// Run fit series: dir=0 (slice Y, proj X), dir=1 (slice X, proj Y)
// Returns ordered: center first, then alternating up/down
// ============================================================================
std::vector<FitResult> runSeries(TH2* h2, int dir, int sliceWidth)
{
    TAxis* slAx = (dir == 0) ? h2->GetYaxis() : h2->GetXaxis();

    double fitMin, fitMax, expPeak, searchWin, sigG, sigMin, sigMax;
    double slMin, slMax, slCen;

    if (dir == 0) {
        fitMin = 175; fitMax = 185;
        expPeak = 180.0; searchWin = 3.0;
        sigG = 0.8; sigMin = 0.3; sigMax = 3.0;
        slMin = 0.45; slMax = 0.62; slCen = 0.54;
    } else {
        fitMin = 0.45; fitMax = 0.62;
        expPeak = 0.54; searchWin = 0.05;
        sigG = 0.012; sigMin = 0.005; sigMax = 0.05;
        slMin = 175; slMax = 185; slCen = 180.0;
    }

    int cenBin = slAx->FindBin(slCen);
    int bLo = slAx->FindBin(slMin), bHi = slAx->FindBin(slMax);
    int cStart = cenBin - sliceWidth / 2;

    // Build ordered slice list: center, then up1, down1, up2, down2, ...
    std::vector<int> starts_up, starts_down;
    for (int s = cStart; s + sliceWidth - 1 <= bHi; s += sliceWidth) starts_up.push_back(s);
    for (int s = cStart - sliceWidth; s >= bLo; s -= sliceWidth)     starts_down.push_back(s);

    // Interleave: center, up1, down1, up2, down2, ...
    std::vector<int> ordered;
    int nu = starts_up.size(), nd = starts_down.size();
    for (int i = 0; i < std::max(nu, nd); i++) {
        if (i < nu) ordered.push_back(starts_up[i]);
        if (i < nd) ordered.push_back(starts_down[i]);
    }

    std::vector<FitResult> results;
    double pMean = expPeak, pSig = sigG, pAmp = 0;
    TString lbl = (dir == 0) ? "pX" : "pY";

    for (int sb : ordered) {
        int eb = std::min(sb + sliceWidth - 1, bHi);
        sb = std::max(sb, bLo);
        double sLo = slAx->GetBinLowEdge(sb), sHi = slAx->GetBinUpEdge(eb);

        TH1D* hp = (dir == 0)
            ? h2->ProjectionX(Form("h_%s_w%d_b%d", lbl.Data(), sliceWidth, sb), sb, eb)
            : h2->ProjectionY(Form("h_%s_w%d_b%d", lbl.Data(), sliceWidth, sb), sb, eb);

        // Use previous good fit as expected peak hint (but findPeakRange will find actual max)
        FitResult r = fitSlice(hp, fitMin, fitMax, pMean, searchWin, pSig, sigMin, sigMax, pAmp);
        r.sliceCenter = (sLo + sHi) / 2.0;
        r.sliceLo = sLo; r.sliceHi = sHi; r.sliceWidth = sliceWidth;

        if (r.valid) { pMean = r.mean; pSig = r.sigma; pAmp = r.amp; }
        results.push_back(r);
    }
    return results;
}

// ============================================================================
// Draw one fit pad
// ============================================================================
void drawPad(TH1D* h, const FitResult& r, double fMin, double fMax)
{
    h->SetLineColor(kBlack); h->SetMarkerStyle(20); h->SetMarkerSize(0.4);
    h->GetXaxis()->SetRangeUser(fMin, fMax);
    h->Draw("E");
    if (!r.valid || r.pars.empty()) {
        TLatex* l = new TLatex(); l->SetNDC(); l->SetTextSize(0.05); l->SetTextColor(kRed);
        l->DrawLatex(0.3, 0.5, "FIT FAILED"); return;
    }
    TF1* ft = new TF1(Form("ft_%s",h->GetName()), Form("gaus(0)+pol%d(3)",r.polyDeg), fMin, fMax);
    for (int i=0;i<(int)r.pars.size();i++) ft->SetParameter(i, r.pars[i]);
    ft->SetLineColor(kRed); ft->SetLineWidth(2); ft->SetNpx(500); ft->Draw("same");

    TF1* fg = new TF1(Form("fg_%s",h->GetName()), "gaus", fMin, fMax);
    fg->SetParameters(r.amp, r.mean, r.sigma);
    fg->SetLineColor(kBlue); fg->SetLineStyle(2); fg->SetLineWidth(2); fg->SetNpx(500); fg->Draw("same");

    TF1* fp = new TF1(Form("fp_%s",h->GetName()), Form("pol%d",r.polyDeg), fMin, fMax);
    for (int i=0;i<=r.polyDeg;i++) fp->SetParameter(i, r.pars[3+i]);
    fp->SetLineColor(kGreen+2); fp->SetLineStyle(3); fp->SetLineWidth(2); fp->SetNpx(500); fp->Draw("same");

    TLatex* l = new TLatex(); l->SetNDC(); l->SetTextSize(0.040);
    l->DrawLatex(0.45, 0.85, Form("#mu=%.4f #sigma=%.4f", r.mean, r.sigma));
    l->DrawLatex(0.45, 0.79, Form("#chi^{2}/ndf=%.2f pol%d", r.chi2ndf, r.polyDeg));
    l->DrawLatex(0.45, 0.73, Form("slice[%.4f,%.4f]", r.sliceLo, r.sliceHi));
}

// ============================================================================
// Weighted average
// ============================================================================
struct WAvg {
    double sw, swx;
    WAvg() : sw(0), swx(0) {}
    void add(double v, double e, double w) { if(e>0&&w>0){double wt=w/(e*e); swx+=v*wt; sw+=wt;} }
    double val() const { return sw>0 ? swx/sw : 0; }
    double err() const { return sw>0 ? 1.0/sqrt(sw) : 0; }
    bool   ok()  const { return sw>0; }
};

// ============================================================================
// MAIN
// ============================================================================
void analyze_elastic_pp(const char* filename = "output_pp.root")
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    TFile* fin = TFile::Open(filename, "READ");
    if (!fin || fin->IsZombie()) { std::cerr << "Cannot open " << filename << std::endl; return; }
    fin->cd("debug");
    TH2* h2 = (TH2*)gDirectory->Get("dphi_vs_tantan");
    if (!h2) { std::cerr << "Cannot find dphi_vs_tantan" << std::endl; return; }
    std::cout << "Loaded: " << h2->GetEntries() << " entries, "
              << h2->GetNbinsX() << "x" << h2->GetNbinsY() << " bins" << std::endl;

    const int W[] = {10, 20, 30, 40};
    const int nW = 4;
    std::vector<FitResult> resX[nW], resY[nW];

    // ---- SERIES A ----
    std::cout << "\n===== SERIES A: dphi peak (slice Y, project X) =====" << std::endl;
    std::cout << "  Peak found dynamically, mean range = 6%-drop zone" << std::endl;
    for (int i = 0; i < nW; i++) {
        resX[i] = runSeries(h2, 0, W[i]);
        int nv = 0;
        for (int j = 0; j < (int)resX[i].size(); j++) {
            auto& r = resX[i][j];
            if (r.valid) {
                nv++;
                if (i == 0)
                    std::cout << Form("    Y=[%.4f,%.4f] mu=%.3f sig=%.3f amp=%.0f chi2=%.2f p%d",
                        r.sliceLo, r.sliceHi, r.mean, r.sigma, r.amp, r.chi2ndf, r.polyDeg) << std::endl;
            } else if (i == 0 && j < 4) {
                std::cout << Form("    Y=[%.4f,%.4f] ** FAILED **", r.sliceLo, r.sliceHi) << std::endl;
            }
        }
        std::cout << "  w=" << W[i] << ": " << nv << "/" << resX[i].size() << " valid" << std::endl;
    }

    // ---- SERIES B ----
    std::cout << "\n===== SERIES B: tantan peak (slice X, project Y) =====" << std::endl;
    for (int i = 0; i < nW; i++) {
        resY[i] = runSeries(h2, 1, W[i]);
        int nv = 0;
        for (auto& r : resY[i]) {
            if (r.valid) {
                nv++;
                if (i == 0)
                    std::cout << Form("    X=[%.2f,%.2f] mu=%.5f sig=%.5f amp=%.0f chi2=%.2f p%d",
                        r.sliceLo, r.sliceHi, r.mean, r.sigma, r.amp, r.chi2ndf, r.polyDeg) << std::endl;
            }
        }
        std::cout << "  w=" << W[i] << ": " << nv << "/" << resY[i].size() << " valid" << std::endl;
    }

    // ---- Combine ----
    WAvg aMx, aSx, aMy, aSy;
    for (int i = 0; i < nW; i++) {
        for (auto& r : resX[i]) { if (r.valid) { aMx.add(r.mean, r.meanErr, r.amp); aSx.add(r.sigma, r.sigmaErr, r.amp); } }
        for (auto& r : resY[i]) { if (r.valid) { aMy.add(r.mean, r.meanErr, r.amp); aSy.add(r.sigma, r.sigmaErr, r.amp); } }
    }

    double cx = aMx.ok() ? aMx.val() : 180.2;
    double sx = aSx.ok() ? aSx.val() : 1.0;
    double cy = aMy.ok() ? aMy.val() : 0.54;
    double sy = aSy.ok() ? aSy.val() : 0.015;

    if (!aMx.ok()) std::cout << "\n  *** WARNING: X-direction failed ***" << std::endl;
    if (!aMy.ok()) std::cout << "\n  *** WARNING: Y-direction failed ***" << std::endl;

    std::cout << "\n========== RESULTS ==========" << std::endl;
    std::cout << "  dphi:   center=" << cx << " +/- " << aMx.err() << "  sigma=" << sx << " +/- " << aSx.err() << std::endl;
    std::cout << "  tantan: center=" << cy << " +/- " << aMy.err() << "  sigma=" << sy << " +/- " << aSy.err() << std::endl;
    std::cout << "  1sig: dphi[" << cx-sx << "," << cx+sx << "] tantan[" << cy-sy << "," << cy+sy << "]" << std::endl;
    std::cout << "  3sig: dphi[" << cx-3*sx << "," << cx+3*sx << "] tantan[" << cy-3*sy << "," << cy+3*sy << "]" << std::endl;
    std::cout << "  5sig: dphi[" << cx-5*sx << "," << cx+5*sx << "] tantan[" << cy-5*sy << "," << cy+5*sy << "]" << std::endl;

    for (int i = 0; i < nW; i++) {
        WAvg mx, sigx, my, sigy;
        for (auto& r : resX[i]) { if (r.valid) { mx.add(r.mean, r.meanErr, r.amp); sigx.add(r.sigma, r.sigmaErr, r.amp); } }
        for (auto& r : resY[i]) { if (r.valid) { my.add(r.mean, r.meanErr, r.amp); sigy.add(r.sigma, r.sigmaErr, r.amp); } }
        std::cout << "  w=" << W[i] << ": X mu=" << mx.val() << " sig=" << sigx.val()
                  << " | Y mu=" << my.val() << " sig=" << sigy.val() << std::endl;
    }

    // ---- TGraph for drawing, TCutG for file ----
    const int nP = 200;
    auto mkG = [&](double nS, int col, int sty, int lw) -> TGraph* {
        TGraph* g = new TGraph(nP+1);
        for (int i=0;i<=nP;i++) { double t=2*TMath::Pi()*i/nP; g->SetPoint(i, cx+nS*sx*cos(t), cy+nS*sy*sin(t)); }
        g->SetLineColor(col); g->SetLineWidth(lw); g->SetLineStyle(sty); g->SetFillStyle(0);
        return g;
    };
    auto mkC = [&](const char* nm, double nS) -> TCutG* {
        TCutG* c = new TCutG(nm, nP+1);
        c->SetVarX("dphi"); c->SetVarY("tantan");
        for (int i=0;i<=nP;i++) { double t=2*TMath::Pi()*i/nP; c->SetPoint(i, cx+nS*sx*cos(t), cy+nS*sy*sin(t)); }
        return c;
    };

    TGraph* g1=mkG(1,kRed,2,2); TGraph* g3=mkG(3,kGreen+2,9,2); TGraph* g5=mkG(5,kMagenta,1,3);
    TCutG* c1=mkC("cutg_1sigma",1); TCutG* c3=mkC("cutg_3sigma",3);
    TCutG* c5=mkC("cutg_5sigma",5); TCutG* cm=mkC("cutg",5);

    TFile* fout = new TFile("CUT_dphi_tantan_fitted.root", "RECREATE");
    c1->Write(); c3->Write(); c5->Write(); cm->Write(); fout->Close();
    std::cout << "\nSaved: CUT_dphi_tantan_fitted.root (cutg_1sigma, cutg_3sigma, cutg_5sigma, cutg)" << std::endl;

    // ========================================================================
    // CONTROL PLOTS: 8 canvases = 2 directions x 4 widths
    // Each shows slices ordered from center outward (already ordered that way)
    // ========================================================================
    const char* dirName[] = {"X-dir (dphi)", "Y-dir (tantan)"};
    double fRange[2][2] = {{175, 185}, {0.45, 0.62}};

    for (int dir = 0; dir < 2; dir++) {
        for (int iw = 0; iw < nW; iw++) {
            std::vector<FitResult>& res = (dir == 0) ? resX[iw] : resY[iw];
            int nSlices = res.size();
            if (nSlices == 0) continue;

            // Pick up to 8 evenly spaced, starting from index 0 (=center)
            std::vector<int> pick;
            int nShow = std::min(nSlices, 8);
            for (int j = 0; j < nShow; j++)
                pick.push_back(j * nSlices / nShow);

            int nCols = (nShow <= 4) ? nShow : 4;
            int nRows = (nShow + nCols - 1) / nCols;

            TCanvas* cc = new TCanvas(
                Form("c_%s_w%d", (dir==0?"X":"Y"), W[iw]),
                Form("%s  w=%d bins", dirName[dir], W[iw]),
                nCols * 420, nRows * 350);
            cc->Divide(nCols, nRows);

            for (int j = 0; j < nShow; j++) {
                cc->cd(j + 1);
                gPad->SetLeftMargin(0.14);
                gPad->SetBottomMargin(0.12);
                gPad->SetTopMargin(0.08);
                FitResult& r = res[pick[j]];

                TH1D* hp = nullptr;
                if (dir == 0) {
                    int b1 = h2->GetYaxis()->FindBin(r.sliceLo + 1e-4);
                    int b2 = h2->GetYaxis()->FindBin(r.sliceHi - 1e-4);
                    hp = h2->ProjectionX(Form("hc_%s_w%d_%d", (dir==0?"X":"Y"), W[iw], j), b1, b2);
                    hp->SetTitle(Form("Y[%.3f,%.3f] w=%d", r.sliceLo, r.sliceHi, W[iw]));
                    hp->GetXaxis()->SetTitle("#Delta#phi [deg]");
                } else {
                    int b1 = h2->GetXaxis()->FindBin(r.sliceLo + 1e-4);
                    int b2 = h2->GetXaxis()->FindBin(r.sliceHi - 1e-4);
                    hp = h2->ProjectionY(Form("hc_%s_w%d_%d", (dir==0?"X":"Y"), W[iw], j), b1, b2);
                    hp->SetTitle(Form("X[%.1f,%.1f] w=%d", r.sliceLo, r.sliceHi, W[iw]));
                    hp->GetXaxis()->SetTitle("tan#theta_{1}#upointtan#theta_{2}");
                }
                hp->GetYaxis()->SetTitle("Counts");
                drawPad(hp, r, fRange[dir][0], fRange[dir][1]);
            }
            cc->Update();
        }
    }

    // ---- Summary canvas ----
    {
        int col[] = {kBlack, kRed, kBlue, kGreen+2};
        int mrk[] = {20, 21, 22, 23};
        TCanvas* cs = new TCanvas("c_summ", "Fit summary", 1400, 900);
        cs->Divide(2, 2);

        cs->cd(1); gPad->SetLeftMargin(0.15);
        bool first = true;
        for (int i=0;i<nW;i++) {
            std::vector<double> vx,vy;
            for (auto& r:resX[i]) if(r.valid){vx.push_back(r.sliceCenter);vy.push_back(r.mean);}
            if(vx.empty()) continue;
            TGraph* g=new TGraph(vx.size(),vx.data(),vy.data());
            g->SetMarkerColor(col[i]);g->SetMarkerStyle(mrk[i]);g->SetMarkerSize(0.8);
            if(first){g->SetTitle("dphi mean vs Y;Y;#mu_{#Delta#phi}");g->Draw("AP");first=false;}
            else g->Draw("P same");
        }
        TLine* la=new TLine(0.45,cx,0.62,cx);la->SetLineColor(kRed);la->SetLineStyle(2);la->Draw();

        cs->cd(2); gPad->SetLeftMargin(0.15);
        first=true;
        for (int i=0;i<nW;i++) {
            std::vector<double> vx,vy;
            for (auto& r:resX[i]) if(r.valid){vx.push_back(r.sliceCenter);vy.push_back(r.sigma);}
            if(vx.empty()) continue;
            TGraph* g=new TGraph(vx.size(),vx.data(),vy.data());
            g->SetMarkerColor(col[i]);g->SetMarkerStyle(mrk[i]);g->SetMarkerSize(0.8);
            if(first){g->SetTitle("dphi #sigma vs Y;Y;#sigma_{#Delta#phi}");g->Draw("AP");first=false;}
            else g->Draw("P same");
        }
        TLine* lb=new TLine(0.45,sx,0.62,sx);lb->SetLineColor(kRed);lb->SetLineStyle(2);lb->Draw();

        cs->cd(3); gPad->SetLeftMargin(0.15);
        first=true;
        for (int i=0;i<nW;i++) {
            std::vector<double> vx,vy;
            for (auto& r:resY[i]) if(r.valid){vx.push_back(r.sliceCenter);vy.push_back(r.mean);}
            if(vx.empty()) continue;
            TGraph* g=new TGraph(vx.size(),vx.data(),vy.data());
            g->SetMarkerColor(col[i]);g->SetMarkerStyle(mrk[i]);g->SetMarkerSize(0.8);
            if(first){g->SetTitle("tantan mean vs X;X [deg];#mu_{tan#theta}");g->Draw("AP");first=false;}
            else g->Draw("P same");
        }
        TLine* lc=new TLine(175,cy,185,cy);lc->SetLineColor(kRed);lc->SetLineStyle(2);lc->Draw();

        cs->cd(4); gPad->SetLeftMargin(0.15);
        first=true;
        for (int i=0;i<nW;i++) {
            std::vector<double> vx,vy;
            for (auto& r:resY[i]) if(r.valid){vx.push_back(r.sliceCenter);vy.push_back(r.sigma);}
            if(vx.empty()) continue;
            TGraph* g=new TGraph(vx.size(),vx.data(),vy.data());
            g->SetMarkerColor(col[i]);g->SetMarkerStyle(mrk[i]);g->SetMarkerSize(0.8);
            if(first){g->SetTitle("tantan #sigma vs X;X [deg];#sigma_{tan#theta}");g->Draw("AP");first=false;}
            else g->Draw("P same");
        }
        TLine* ld=new TLine(175,sy,185,sy);ld->SetLineColor(kRed);ld->SetLineStyle(2);ld->Draw();

        cs->cd(1);
        TLegend* leg=new TLegend(0.55,0.65,0.88,0.88);
        for(int i=0;i<nW;i++){TMarker* m=new TMarker();m->SetMarkerColor(col[i]);m->SetMarkerStyle(mrk[i]);leg->AddEntry(m,Form("w=%d",W[i]),"p");}
        leg->AddEntry(la,"Avg","l"); leg->Draw();
        cs->Update();
    }

    // ---- Final 2D ----
    {
        TCanvas* cf = new TCanvas("c_final", "2D + ellipses", 1100, 800);
        gPad->SetRightMargin(0.14);
        h2->GetXaxis()->SetRangeUser(160, 200);
        h2->GetYaxis()->SetRangeUser(0.0, 1.0);
        h2->SetTitle("#Delta#phi vs tan#theta_{1}#upointtan#theta_{2}");
        h2->Draw("colz");
        g1->Draw("L same"); g3->Draw("L same"); g5->Draw("L same");
        TLegend* lf=new TLegend(0.65,0.78,0.85,0.90);
        lf->SetFillStyle(0);lf->SetBorderSize(0);lf->SetTextColor(kWhite);
        lf->AddEntry(g1,"1#sigma","l");lf->AddEntry(g3,"3#sigma","l");lf->AddEntry(g5,"5#sigma","l");
        lf->Draw();
        TLatex* tx=new TLatex();tx->SetNDC();tx->SetTextColor(kWhite);tx->SetTextSize(0.028);
        tx->DrawLatex(0.15,0.85,Form("Center: (%.3f, %.5f)",cx,cy));
        tx->DrawLatex(0.15,0.81,Form("#sigma_{x}=%.3f  #sigma_{y}=%.5f",sx,sy));
        cf->Update();
    }
    {
        TCanvas* cz = new TCanvas("c_zoom", "2D zoom", 1100, 800);
        gPad->SetRightMargin(0.14);
        h2->GetXaxis()->SetRangeUser(170, 190);
        h2->GetYaxis()->SetRangeUser(0.3, 0.8);
        h2->SetTitle("Zoomed");
        h2->Draw("colz");
        g1->Draw("L same"); g3->Draw("L same"); g5->Draw("L same");
        TLegend* lz=new TLegend(0.65,0.78,0.85,0.90);
        lz->SetFillStyle(0);lz->SetBorderSize(0);
        lz->AddEntry(g1,"1#sigma","l");lz->AddEntry(g3,"3#sigma","l");lz->AddEntry(g5,"5#sigma","l");
        lz->Draw(); cz->Update();
    }

    // ---- Table ----
    std::cout << "\n====== FIT TABLE ======" << std::endl;
    for (int i=0;i<nW;i++) {
        for (auto& r:resX[i]) if(r.valid)
            std::cout << Form(" X w=%2d Y=[%.4f,%.4f] mu=%.4f sig=%.4f amp=%.0f chi2=%.2f p%d",W[i],r.sliceLo,r.sliceHi,r.mean,r.sigma,r.amp,r.chi2ndf,r.polyDeg) << std::endl;
        for (auto& r:resY[i]) if(r.valid)
            std::cout << Form(" Y w=%2d X=[%.2f,%.2f] mu=%.6f sig=%.6f amp=%.0f chi2=%.2f p%d",W[i],r.sliceLo,r.sliceHi,r.mean,r.sigma,r.amp,r.chi2ndf,r.polyDeg) << std::endl;
    }
    std::cout << "\nDONE." << std::endl;
}
