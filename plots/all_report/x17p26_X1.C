// ============================================================================
// x17p26_X1.C -- X17 limit with PDG 2026 (analyst ruling D25, 15 November 2026), prefix x17p26, written 2026-11-15: copy of
// fix_x17X1.C whose ONLY change is the printed limit of the right panel, now the PDG 2026 value 3.1589e-5 (normalisation (1),
// 17 MeV/c2; x17p26_limit.C: frozen code rerun with BR_dir = 2.76e-4 and its PDG 2026 uncertainty in s_X) and the output name
// fix_X1_pdg26 (the old fix_X1 is not touched). Everything drawn is in rows and does not depend on the PDG edition.
// Run: root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_x17p26", kTRUE)' 'plots/all_report/x17p26_X1.C+O'
// --- header of fix_x17X1.C follows ---
// ============================================================================
// fix_x17X1.C -- review fix (prefix fix): Fig. X1 of x17d_figs.C (WP9 X-d, as corrected by the verifier vx) redrawn in the
// common palette of report_all: CB red, other background grey, direct decays D (C5) cyan, conversions C (C4) green (was blue,
// the colour of N - CB), the eta-related term (C2 + C3 + G_eta, ONE smoothed term of the frozen model) light violet, X17 magenta.
// The X1 block is copied unchanged except the colours, the legend text of the eta-related term and the output names; N_X^UL is
// taken from vx_nx.txt (391.9 raw rows at 17 MeV/c2) instead of being recomputed. Reads only (hists_v1.root, x17d/x17c files).
// Run: root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_fix", kTRUE)' 'plots/all_report/fix_x17X1.C+O'
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <map>
#include <string>
#include <vector>

#include "TBox.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "THStack.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TVectorD.h"

#include "x17d_lib.h"

namespace x17dfg {
using namespace x17d;
const double kQ[5] = {0.02275, 0.15866, 0.5, 0.84134, 0.97725};
// provisional sigma(pp -> eta X) for sigma x BR x BR (WP8 not available): a_eta(WP5) x 1.0665 mb x C_sec (1 - f_eta) (WP2)
// [vx verification 2026-11-14: was C_sec (1 - f_win) = 0.842; SPEC AMENDMENT 1 (08:20) makes C_sec x (1 - f_eta) binding for
//  every eta yield: C_sec(eta) 1.0110 x (1 - f_eta) 0.97 = 0.9807 (+0.03/-0.04); (1 - f_win) would remove the material twice]
const double kCsecFwin = 1.0110 * 0.97;
struct Obs { double m, RX, ul[3], rhat[3], rerr[3], ea[3][5], q0[3], Z[3], p0[3], nX[3], N5, ulC2, s1x; };
struct Frz { double m, RX, e[3][5], nX; };
struct ToyBO { double m; std::vector<double> ul[3], pull; };
struct TC { double r, br, cls, err, clsb, clb, clsA, qobs; int nt; };
struct TCRes { double m = 0; int norm = 0; double ulA = 0, ulT = 0, ratio = 0; int clb0 = 0; std::vector<TC> pts; bool ok = false; };
struct FC { double r, br, tobs, t68, t90, t95, p; int nt, nok; };
struct MCr { int kind, idx; double ul[3], rh[3]; };

bool ReadToysBO(double m, const std::vector<int>& firsts, ToyBO& T) {
    T.m = m;
    for (int fst : firsts) {
        TFile* f = TFile::Open(Form("/hdd3/przygoda/HADES/FAT_out/all_report/x17c/toys_m%04d_i0000_%04d.root", (int)std::lround(100 * m), fst), "READ");
        if (!f || f->IsZombie()) return false;
        TTree* t = (TTree*)f->Get("toys");
        int ok[3]; double ul[3], rh[3], re[3];
        t->SetBranchAddress("ok", ok); t->SetBranchAddress("ul", ul); t->SetBranchAddress("rhat", rh); t->SetBranchAddress("rerr", re);
        for (Long64_t k = 0; k < t->GetEntries(); ++k) {
            t->GetEntry(k);
            for (int a = 0; a < 3; ++a) if (ok[a] && std::isfinite(ul[a])) T.ul[a].push_back(ul[a]);
            if (ok[0] && re[0] > 0) T.pull.push_back(rh[0] / re[0]);
        }
        f->Close();
    }
    return true;
}
TCRes ReadTC(double m, int norm) {
    TCRes R; R.m = m; R.norm = norm;
    FILE* f = fopen(Form("%s/toycls_m%04d_n%d.txt", kOutD, (int)std::lround(100 * m), norm), "r");
    if (!f) return R;
    char line[1024];
    while (fgets(line, sizeof line, f)) {
        if (line[0] == '#') continue;
        if (strncmp(line, "RESULT", 6) == 0) {
            double mm; int nn;
            if (sscanf(line, "RESULT m %lf norm %d UL_asym %lf UL_toy %lf ratio %lf clb0 %d", &mm, &nn, &R.ulA, &R.ulT, &R.ratio, &R.clb0) == 6) R.ok = true;
            continue;
        }
        TC p;
        if (sscanf(line, "%lf %lf %lf %lf %lf %lf %d %lf %lf", &p.r, &p.br, &p.cls, &p.err, &p.clsb, &p.clb, &p.nt, &p.clsA, &p.qobs) == 9) R.pts.push_back(p);
    }
    fclose(f);
    return R;
}
std::vector<FC> ReadFC(double m, int norm) {
    std::vector<FC> v;
    void* d = gSystem->OpenDirectory(kOutD);
    const char* e;
    std::vector<std::string> files;
    const std::string pre = Form("fc_m%04d_n%d_r", (int)std::lround(100 * m), norm);
    while ((e = gSystem->GetDirEntry(d))) { std::string s(e); if (s.rfind(pre, 0) == 0) files.push_back(s); }
    gSystem->FreeDirectory(d);
    for (const auto& s : files) {
        FILE* f = fopen(Form("%s/%s", kOutD, s.c_str()), "r");
        if (!f) continue;
        char line[1024];
        while (fgets(line, sizeof line, f)) {
            if (line[0] == '#') continue;
            FC p;
            if (sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %d %d", &p.r, &p.br, &p.tobs, &p.t68, &p.t90, &p.t95, &p.p, &p.nt, &p.nok) == 9) v.push_back(p);
        }
        fclose(f);
    }
    std::sort(v.begin(), v.end(), [](const FC& a, const FC& b) { return a.r < b.r; });
    return v;
}
std::vector<MCr> ReadMC(double m) {
    std::vector<MCr> v;
    FILE* f = fopen(Form("%s/mc_m%04d.txt", kOutD, (int)std::lround(100 * m)), "r");
    if (!f) return v;
    char line[1024];
    while (fgets(line, sizeof line, f)) {
        if (line[0] == '#') continue;
        MCr r; int o[3];
        if (sscanf(line, "%d %d %lf %lf %lf %lf %lf %lf %d %d %d", &r.kind, &r.idx, &r.ul[0], &r.ul[1], &r.ul[2], &r.rh[0], &r.rh[1], &r.rh[2], &o[0], &o[1], &o[2]) == 11) v.push_back(r);
    }
    fclose(f);
    return v;
}
double Rms(const std::vector<double>& v) {
    if (v.size() < 2) return 0;
    double m = 0; for (double x : v) m += x; m /= v.size();
    double s = 0; for (double x : v) s += (x - m) * (x - m);
    return std::sqrt(s / (v.size() - 1));
}
double Mean(const std::vector<double>& v) { double m = 0; for (double x : v) m += x; return v.empty() ? 0 : m / v.size(); }
void StyleAxis(TH1* h, double ts = 0.05, double ls = 0.045) {
    h->GetXaxis()->SetTitleSize(ts); h->GetYaxis()->SetTitleSize(ts); h->GetXaxis()->SetLabelSize(ls); h->GetYaxis()->SetLabelSize(ls);
}
}  // namespace x17dfg
void x17p26_X1() {
    using namespace x17dfg;
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    Ctx c;
    if (!Init(c)) { printf("x17p26_X1: init failed\n"); return; }
    x17c::Inputs& I = c.I;
    const std::vector<double>& n = c.D.n;
    const int nx = I.G.nx, ny = I.G.ny;
    const x17c::Config C1 = MakeConfig(c.z, 1);
    x17c::Inputs I17 = I;
    SetLine(I17, 17.0);
    const x17c::Model M1(I17, C1);
    const x17c::Profile Pf = x17c::ProfileFits(M1, n, x17c::GlobOf(C1));
    const std::vector<double>& p0 = Pf.zero.p;
    const double nXul = 391.9;            // vx_nx.txt, 17.00 MeV/c2
    struct { double ul[3]; } o17v = {{3.1589e-5, 0, 0}};   // x17p26: PDG 2026 (was 3.0636e-5, PDG 2020 frozen)
    const auto* o17 = &o17v;
    // ------------------------------------------------------------ Fig. X1: M(ee) spectra of the peak band
    {
        const double pkLo = 0.5325, pkHi = 0.5575;
        TH1D* hN = new TH1D("x1N", "", ny, 1e3 * I.G.ye.front(), 1e3 * I.G.ye.back());
        TH1D* hc[6];
        const char* cn[6] = {"x1CB", "x1Bg", "x1D", "x1C", "x1Et", "x1X"};
        for (int k = 0; k < 6; ++k) { hc[k] = new TH1D(cn[k], "", ny, 1e3 * I.G.ye.front(), 1e3 * I.G.ye.back()); }
        TH1D* hX = new TH1D("x1Xul", "", ny, 1e3 * I.G.ye.front(), 1e3 * I.G.ye.back());
        TH1D* hTot = new TH1D("x1Tot", "", ny, 1e3 * I.G.ye.front(), 1e3 * I.G.ye.back());
        TH1D* hSig = new TH1D("x1Sig", "", ny, 1e3 * I.G.ye.front(), 1e3 * I.G.ye.back());   // N - (CB + other + eta-related)
        TH1D* hDC = new TH1D("x1DC", "", ny, 1e3 * I.G.ye.front(), 1e3 * I.G.ye.back());
        TH1D* hD = new TH1D("x1Donly", "", ny, 1e3 * I.G.ye.front(), 1e3 * I.G.ye.back());
        TH1D* hDCX = new TH1D("x1DCX", "", ny, 1e3 * I.G.ye.front(), 1e3 * I.G.ye.back());
        TH1D* hPull = new TH1D("x1Pull", "", ny, 1e3 * I.G.ye.front(), 1e3 * I.G.ye.back());
        TH1D* hPullX = new TH1D("x1PullX", "", ny, 1e3 * I.G.ye.front(), 1e3 * I.G.ye.back());
        for (int j = 0; j < ny; ++j) {
            double sn = 0, s[6] = {0, 0, 0, 0, 0, 0}, sx = 0;
            for (int i = 0; i < nx; ++i) {
                if (!(I.G.xc(i) > pkLo && I.G.xc(i) < pkHi)) continue;
                double cc[6]; Comps(M1, p0, i, j, cc);
                for (int k = 0; k < 6; ++k) s[k] += cc[k];
                sn += n[I.G.idx(i, j)];
                sx += nXul * I17.gX[i] * I17.X[j];
            }
            const double w = 0.5;   // per MeV/c2 (2 MeV/c2 bins)
            hN->SetBinContent(j + 1, sn * w); hN->SetBinError(j + 1, std::sqrt(sn) * w);
            for (int k = 0; k < 6; ++k) hc[k]->SetBinContent(j + 1, s[k] * w);
            const double tot = s[0] + s[1] + s[2] + s[3] + s[4];
            hTot->SetBinContent(j + 1, tot * w);
            hX->SetBinContent(j + 1, (tot + sx) * w);
            hSig->SetBinContent(j + 1, (sn - s[0] - s[1] - s[4]) * w); hSig->SetBinError(j + 1, std::sqrt(sn) * w);
            hD->SetBinContent(j + 1, s[2] * w);
            hDC->SetBinContent(j + 1, (s[2] + s[3]) * w);
            hDCX->SetBinContent(j + 1, (s[2] + s[3] + sx) * w);
            hPull->SetBinContent(j + 1, tot > 0 ? (sn - tot) / std::sqrt(tot) : 0.); hPull->SetBinError(j + 1, 1.);
            hPullX->SetBinContent(j + 1, tot > 0 ? sx / std::sqrt(tot) : 0.);
        }
        TCanvas cv("cx1", "", 1600, 800);
        // left: full composition (log)
        TPad* l1 = new TPad("l1", "", 0.0, 0.30, 0.5, 1.0), *l2 = new TPad("l2", "", 0.0, 0.0, 0.5, 0.30);
        TPad* r1 = new TPad("r1", "", 0.5, 0.30, 1.0, 1.0), *r2 = new TPad("r2", "", 0.5, 0.0, 1.0, 0.30);
        for (TPad* p : {l1, r1}) { p->SetBottomMargin(0.02); p->SetLeftMargin(0.15); p->SetRightMargin(0.03); p->SetTopMargin(0.05); }
        for (TPad* p : {l2, r2}) { p->SetTopMargin(0.02); p->SetBottomMargin(0.34); p->SetLeftMargin(0.15); p->SetRightMargin(0.03); }
        l1->SetLogy();
        cv.cd(); l1->Draw(); l2->Draw(); r1->Draw(); r2->Draw();
        l1->cd();
        THStack* st = new THStack("x1st", ";M(e^{+}e^{-}) [MeV/c^{2}];rows per MeV/c^{2}");
        const int col[6] = {kRed + 1, kGray + 1, kCyan + 1, kGreen + 2, kViolet - 9, kMagenta};
        for (int k : {0, 4, 3, 2, 1}) { hc[k]->SetFillColor(col[k]); hc[k]->SetLineColor(col[k]); st->Add(hc[k]); }
        st->SetMinimum(5); st->SetMaximum(4e5);
        st->Draw("HIST");
        st->GetYaxis()->SetTitleSize(0.055); st->GetYaxis()->SetLabelSize(0.05); st->GetXaxis()->SetLabelSize(0); st->GetYaxis()->SetTitleOffset(1.2);
        hX->SetLineColor(kMagenta + 1); hX->SetLineWidth(3); hX->SetLineStyle(2); hX->Draw("HIST SAME");
        hN->SetMarkerStyle(20); hN->SetMarkerSize(1.1); hN->SetMarkerColor(kBlack); hN->SetLineColor(kBlack); hN->Draw("E1 SAME");
        TLegend* lg = new TLegend(0.47, 0.50, 0.97, 0.94); lg->SetBorderSize(0); lg->SetTextSize(0.037); lg->SetFillStyle(0);
        lg->AddEntry(hN, "data N (S0, isBest #geq 0)", "pe");
        lg->AddEntry(hc[1], "other background n_{j} B(M4)", "f");
        lg->AddEntry(hc[2], "#eta#rightarrow#pi^{+}#pi^{-}e^{+}e^{-} (D, fitted)", "f");
        lg->AddEntry(hc[3], "#eta#rightarrow#pi^{+}#pi^{-}#gamma, conv. (C, fitted)", "f");
        lg->AddEntry(hc[4], "#eta-related: C2+C3+G_{#eta} (one term)", "f");
        lg->AddEntry(hc[0], "CB (like-sign data)", "f");
        lg->AddEntry(hX, Form("+ X17 (17 MeV/c^{2}) at UL: N_{X} = %.0f", nXul), "l");
        lg->Draw();
        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.042);
        tx.SetTextSize(0.036);
        tx.DrawLatex(0.18, 0.89, "M4 0.5325-0.5575 GeV/c^{2}");
        tx.DrawLatex(0.18, 0.84, "(peak band)");
        tx.DrawLatex(0.18, 0.79, "fit r = 0, 1050 cells");
        l2->cd();
        hPull->SetTitle(";M(e^{+}e^{-}) [MeV/c^{2}];(N - #mu_{0})/#sqrt{#mu_{0}}");
        StyleAxis(hPull, 0.12, 0.10);
        hPull->GetYaxis()->SetTitleSize(0.10); hPull->GetYaxis()->SetTitleOffset(0.6); hPull->GetYaxis()->SetNdivisions(505);
        hPull->SetMinimum(-3.5); hPull->SetMaximum(3.5); hPull->SetMarkerStyle(20); hPull->SetMarkerColor(kBlack); hPull->SetLineColor(kBlack);
        hPull->Draw("E1");
        hPullX->SetLineColor(kMagenta + 1); hPullX->SetLineWidth(3); hPullX->SetLineStyle(2); hPullX->Draw("HIST SAME");
        TLine* z0 = new TLine(2, 0, 62, 0); z0->SetLineStyle(2); z0->Draw();
        TBox* bl = new TBox(12, -3.5, 28, 3.5); bl->SetFillColorAlpha(kYellow - 9, 0.35); bl->SetLineWidth(0); bl->Draw();
        hPull->Draw("E1 SAME");
        // right: peak-band data minus the non-peak model terms, with D, C and the X17 line at the limit
        r1->cd();
        double ymax = 0, ymin = 0;
        for (int j = 1; j <= ny; ++j) { ymax = std::max(ymax, hSig->GetBinContent(j) + hSig->GetBinError(j)); ymin = std::min(ymin, hSig->GetBinContent(j) - hSig->GetBinError(j)); }
        hSig->SetTitle(";M(e^{+}e^{-}) [MeV/c^{2}];N - (CB + other + #eta-rel.) per MeV/c^{2}");
        StyleAxis(hSig, 0.05, 0.05);
        hSig->GetXaxis()->SetLabelSize(0); hSig->GetYaxis()->SetTitleOffset(1.25); hSig->GetYaxis()->SetTitleSize(0.050);
        hSig->SetMaximum(1.25 * ymax); hSig->SetMinimum(std::min(-30., 1.2 * ymin));
        hSig->SetMarkerStyle(20); hSig->SetMarkerColor(kBlack); hSig->SetLineColor(kBlack);
        hSig->Draw("E1");
        THStack* st2 = new THStack("x1st2", "");
        TH1D* hD2 = (TH1D*)hc[2]->Clone("x1D2"); TH1D* hC2 = (TH1D*)hc[3]->Clone("x1C2");
        st2->Add(hD2); st2->Add(hC2);
        st2->Draw("HIST SAME");
        hDCX->SetLineColor(kMagenta + 1); hDCX->SetLineWidth(3); hDCX->SetLineStyle(2); hDCX->Draw("HIST SAME");
        hSig->Draw("E1 SAME");
        TLine* z1 = new TLine(2, 0, 62, 0); z1->SetLineStyle(3); z1->Draw();
        TLegend* lg2 = new TLegend(0.43, 0.48, 0.97, 0.78); lg2->SetBorderSize(0); lg2->SetTextSize(0.040); lg2->SetFillStyle(0);
        lg2->AddEntry(hSig, "data N minus non-peak model terms", "pe");
        lg2->AddEntry(hD2, "direct decays D (fitted N5)", "f");
        lg2->AddEntry(hC2, "conversions C (fitted N4)", "f");
        lg2->AddEntry(hDCX, Form("D + C + X17 at the UL (%.0f rows)", nXul), "l");
        lg2->Draw();
        tx.SetTextSize(0.038);
        tx.SetTextSize(0.036); tx.DrawLatex(0.27, 0.88, Form("90%% CL_{s} UL (1), 17 MeV/c^{2}: BR_{X}BR_{ee} < %.2f#times10^{-5} (PDG 2026)", o17 ? 1e5 * o17->ul[0] : 0.));   // [vx: x 0.44 -> 0.40, exponent was clipped]
        tx.DrawLatex(0.44, 0.82, "M4 0.5325-0.5575 GeV/c^{2}; S0; L = 6079.2 nb^{-1}");
        r2->cd();
        TH1D* hRes = (TH1D*)hSig->Clone("x1Res");
        for (int j = 1; j <= ny; ++j) { hRes->SetBinContent(j, hSig->GetBinContent(j) - hDC->GetBinContent(j)); }
        hRes->SetTitle(";M(e^{+}e^{-}) [MeV/c^{2}];minus D+C");
        StyleAxis(hRes, 0.12, 0.10);
        hRes->GetXaxis()->SetLabelSize(0.10); hRes->GetYaxis()->SetTitleSize(0.10); hRes->GetYaxis()->SetTitleOffset(0.6); hRes->GetYaxis()->SetNdivisions(505);
        double rr = 0; for (int j = 1; j <= ny; ++j) rr = std::max(rr, std::fabs(hRes->GetBinContent(j)) + hRes->GetBinError(j));
        hRes->SetMaximum(1.1 * rr); hRes->SetMinimum(-1.1 * rr);
        hRes->Draw("E1");
        TH1D* hXo = (TH1D*)hDCX->Clone("x1Xo"); hXo->Add(hDC, -1); hXo->Draw("HIST SAME");
        TLine* z2 = new TLine(2, 0, 62, 0); z2->SetLineStyle(2); z2->Draw();
        TBox* br2 = new TBox(12, -1.1 * rr, 28, 1.1 * rr); br2->SetFillColorAlpha(kYellow - 9, 0.35); br2->SetLineWidth(0); br2->Draw();
        hRes->Draw("E1 SAME");
        cv.SaveAs(Form("%s/fix_X1_pdg26.pdf", kFigD)); cv.SaveAs(Form("%s/fix_X1_pdg26.png", kFigD));
    }

}
