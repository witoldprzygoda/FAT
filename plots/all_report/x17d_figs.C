// ============================================================================
// x17d_figs.C  --  WP9 X-d (prefix x17d): figures X1-X4 and Table D1 of the X17 unblinding
// ============================================================================
// Inputs: x17d_state.root (x17d_unblind.C), toycls_*.txt (x17d_toycls.C), fc_*.txt (x17d_fc.C), mc_*.txt (x17d_mc.C);
//         the FROZEN expectation: x17c/x17c_state.root (Asimov bands vs m_X) and x17c/toys_*.root (2000 background-only toys
//         at the quoted masses, 500 at 10/25/40 MeV/c2), read only.
// Output: docs/latex/figs_all/x17d_X1_mee, x17d_X2_limit, x17d_X3_p0, x17d_X4_valid (.pdf + .png);
//         docs/latex/tables_all/x17d_D1.txt, x17d_D1.tex, x17d_checks.txt
// Run: root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_x17d", kTRUE)' \
//        'plots/all_report/x17d_figs.C+O'
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

void x17d_figs() {
    using namespace x17dfg;
    gStyle->SetOptStat(0);
    gStyle->SetEndErrorSize(0);
    // ------------------------------------------------------------ inputs
    Ctx c;
    if (!Init(c)) { printf("x17d_figs: init failed\n"); return; }
    x17c::Inputs& I = c.I;
    const std::vector<double>& n = c.D.n;
    const int nx = I.G.nx, ny = I.G.ny;
    TFile* fs = TFile::Open(Form("%s/x17d_state.root", kOutD), "READ");
    TTree* to = (TTree*)fs->Get("obs");
    Obs o;
    to->SetBranchAddress("m", &o.m); to->SetBranchAddress("RX", &o.RX); to->SetBranchAddress("ul", o.ul); to->SetBranchAddress("rhat", o.rhat);
    to->SetBranchAddress("rerr", o.rerr); to->SetBranchAddress("expA", o.ea); to->SetBranchAddress("q0", o.q0); to->SetBranchAddress("Z", o.Z);
    to->SetBranchAddress("p0", o.p0); to->SetBranchAddress("nX", o.nX); to->SetBranchAddress("N5", &o.N5); to->SetBranchAddress("ulCos2", &o.ulC2);
    to->SetBranchAddress("s1x", &o.s1x);
    std::vector<Obs> obs;
    for (Long64_t k = 0; k < to->GetEntries(); ++k) { to->GetEntry(k); obs.push_back(o); }
    TVectorD* misc = (TVectorD*)fs->Get("misc");
    TFile* fz = TFile::Open("/hdd3/przygoda/HADES/FAT_out/all_report/x17c/x17c_state.root", "READ");
    TTree* te = (TTree*)fz->Get("exp");
    Frz fr;
    double e1[5], e2[5], e3[5];
    te->SetBranchAddress("m", &fr.m); te->SetBranchAddress("RX", &fr.RX); te->SetBranchAddress("e1", e1); te->SetBranchAddress("e2", e2); te->SetBranchAddress("e3", e3);
    te->SetBranchAddress("nX", &fr.nX);
    std::vector<Frz> frz;
    for (Long64_t k = 0; k < te->GetEntries(); ++k) {
        te->GetEntry(k);
        for (int q = 0; q < 5; ++q) { fr.e[0][q] = e1[q]; fr.e[1][q] = e2[q]; fr.e[2][q] = e3[q]; }
        frz.push_back(fr);
    }
    auto obsAt = [&](double m) -> const Obs* { for (const auto& x : obs) if (std::fabs(x.m - m) < 1e-6) return &x; return nullptr; };
    auto frzAt = [&](double m) -> const Frz* { for (const auto& x : frz) if (std::fabs(x.m - m) < 1e-6) return &x; return nullptr; };
    std::vector<ToyBO> bo;
    for (double m : {16.7, 16.85, 17.0}) { ToyBO T; if (ReadToysBO(m, {0, 500, 1000, 1500}, T)) bo.push_back(T); }
    for (double m : {10.0, 25.0, 40.0}) { ToyBO T; if (ReadToysBO(m, {0}, T)) bo.push_back(T); }
    auto boAt = [&](double m) -> const ToyBO* { for (const auto& x : bo) if (std::fabs(x.m - m) < 1e-6) return &x; return nullptr; };
    std::map<std::pair<int, int>, TCRes> tcs;
    for (double m : kQuoted) for (int a = 1; a <= 3; ++a) tcs[{(int)std::lround(100 * m), a}] = ReadTC(m, a);
    const std::vector<FC> fc = ReadFC(17.0, 1);
    std::map<int, std::vector<MCr>> mcs;
    for (double m : kQuoted) mcs[(int)std::lround(100 * m)] = ReadMC(m);

    // ------------------------------------------------------------ the full fits (r = 0) for the components
    const x17c::Config C1 = MakeConfig(c.z, 1);
    x17c::Inputs I17 = I;
    SetLine(I17, 17.0);
    const x17c::Model M1(I17, C1);
    const x17c::Profile Pf = x17c::ProfileFits(M1, n, x17c::GlobOf(C1));
    const std::vector<double>& p0 = Pf.zero.p;
    const Obs* o17 = obsAt(17.0);
    // [vx verification 2026-11-14] The previous N_X (o17->nX = r_UL x N5(r = 0 fit) x R_X = 563 rows) is a conversion of the BR
    // limit, not a limit on the number of X17 rows: at the UL the conditional fit has N5 2687 and theta_X -0.46, i.e. 374 X17
    // rows in the model. The spec asks for a limit on N_X: a direct 90 % CLs limit with POI = rows (normalisation-(2) form,
    // X17 rows = r xA N5ref R_X, yield nuisance off because rows need no efficiency or BR factor; everything else frozen).
    std::map<int, double> nXrows;
    for (double m : kQuoted) {
        x17c::Inputs Im = I;
        if (!SetLine(Im, m)) continue;
        x17c::Config Cr = MakeConfig(c.z, 2);
        Cr.xSig = 0.;
        const Lim l = LimitBR(Im, Cr, n, false);
        nXrows[(int)std::lround(100 * m)] = l.ok ? l.ul / kBR * Cr.xA * Im.N5ref * Im.RX : NAN;
        printf("x17d_figs: [vx] direct 90%% CLs UL on N_X at %.2f MeV/c2: %.1f raw rows (r_UL N5(r=0) R_X of normalisation 1: %.1f)\n", m,
               nXrows[(int)std::lround(100 * m)], obsAt(m) ? obsAt(m)->nX[0] : NAN);
    }
    const double nXul = nXrows[1700];   // 90 % CLs UL on the X17 rows at 17 MeV/c2 (raw rows, X2D M4 range, all M(ee))

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
        const int col[6] = {kRed + 1, kGray + 1, kCyan + 1, kAzure + 2, kViolet + 2, kMagenta};
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
        lg->AddEntry(hc[4], "#eta-related (C2+C3, G_{#eta})", "f");
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
        tx.DrawLatex(0.40, 0.88, Form("90%% CL_{s} UL (1), 17 MeV/c^{2}: BR_{X}BR_{ee} < %.2f#times10^{-5}", o17 ? 1e5 * o17->ul[0] : 0.));   // [vx: x 0.44 -> 0.40, exponent was clipped]
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
        cv.SaveAs(Form("%s/x17d_X1_mee.pdf", kFigD)); cv.SaveAs(Form("%s/x17d_X1_mee.png", kFigD));
    }

    // ------------------------------------------------------------ Fig. X2: limit vs m_X
    {
        TCanvas cv("cx2", "", 1000, 780);
        gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.12); gPad->SetLogy(); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.04);
        const int nf = (int)frz.size();
        TGraphAsymmErrors* b2 = new TGraphAsymmErrors(nf), *b1 = new TGraphAsymmErrors(nf);
        TGraph* fmed = new TGraph(nf);
        for (int k = 0; k < nf; ++k) {
            const Frz& f = frz[k];
            b2->SetPoint(k, f.m, f.e[0][2]); b2->SetPointError(k, 0, 0, f.e[0][2] - f.e[0][0], f.e[0][4] - f.e[0][2]);
            b1->SetPoint(k, f.m, f.e[0][2]); b1->SetPointError(k, 0, 0, f.e[0][2] - f.e[0][1], f.e[0][3] - f.e[0][2]);
            fmed->SetPoint(k, f.m, f.e[0][2]);
        }
        TGraph* g1 = new TGraph(), *g2 = new TGraph(), *g3 = new TGraph(), *gA = new TGraph(), *gC = new TGraph();
        for (const auto& x : obs) {
            g1->SetPoint(g1->GetN(), x.m, x.ul[0]); g2->SetPoint(g2->GetN(), x.m, x.ul[1]); g3->SetPoint(g3->GetN(), x.m, x.ul[2]);
            gA->SetPoint(gA->GetN(), x.m, x.ea[0][2]); gC->SetPoint(gC->GetN(), x.m, x.ulC2);
        }
        b2->SetFillColor(kOrange - 2); b2->SetLineColor(kOrange - 2); b1->SetFillColor(kGreen + 1); b1->SetLineColor(kGreen + 1);
        b2->SetTitle(";m_{X} [MeV/c^{2}];90% CL_{s} UL on BR(#eta#rightarrow#pi^{+}#pi^{-}X) BR(X#rightarrowe^{+}e^{-})");
        b2->GetXaxis()->SetTitleSize(0.048); b2->GetYaxis()->SetTitleSize(0.042); b2->GetXaxis()->SetLabelSize(0.042); b2->GetYaxis()->SetLabelSize(0.042);
        b2->GetYaxis()->SetTitleOffset(1.5);
        b2->SetMinimum(3e-6); b2->SetMaximum(4e-3); b2->GetXaxis()->SetLimits(9.5, 40.5);
        b2->Draw("A3"); b1->Draw("3 SAME");
        fmed->SetLineStyle(2); fmed->SetLineWidth(2); fmed->Draw("L SAME");
        gA->SetLineColor(kGray + 2); gA->SetLineStyle(3); gA->SetLineWidth(2); gA->Draw("L SAME");
        g2->SetLineColor(kBlue + 1); g2->SetLineWidth(2); g2->SetLineStyle(7); g2->Draw("L SAME");
        g3->SetLineColor(kMagenta + 2); g3->SetLineWidth(2); g3->SetLineStyle(9); g3->Draw("L SAME");
        gC->SetLineColor(kRed + 1); gC->SetLineWidth(2); gC->SetLineStyle(4); gC->Draw("L SAME");
        g1->SetLineColor(kBlack); g1->SetLineWidth(3); g1->SetMarkerStyle(20); g1->SetMarkerSize(0.9); g1->Draw("LP SAME");
        // frozen toy expectation at the quoted masses (and 10/25/40)
        TGraphAsymmErrors* tq = new TGraphAsymmErrors();
        for (const auto& T : bo) {
            const double md = x17d::Quantile(T.ul[0], 0.5), lo = x17d::Quantile(T.ul[0], kQ[1]), hi = x17d::Quantile(T.ul[0], kQ[3]);
            const int k = tq->GetN(); tq->SetPoint(k, T.m + (T.m < 16.9 && T.m > 16.6 ? 0 : 0.0), md); tq->SetPointError(k, 0, 0, md - lo, hi - md);
        }
        tq->SetMarkerStyle(24); tq->SetMarkerSize(1.3); tq->SetMarkerColor(kGreen + 3); tq->SetLineColor(kGreen + 3); tq->SetLineWidth(2); tq->Draw("P SAME");
        // toy CLs observed UL at the quoted masses (normalisation 1)
        TGraph* gt = new TGraph();
        for (double m : kQuoted) { const TCRes& R = tcs[{(int)std::lround(100 * m), 1}]; if (R.ok && std::isfinite(R.ulT)) gt->SetPoint(gt->GetN(), m, R.ulT); }
        gt->SetMarkerStyle(29); gt->SetMarkerSize(2.0); gt->SetMarkerColor(kRed + 1); gt->Draw("P SAME");
        TLegend* lg = new TLegend(0.36, 0.63, 0.97, 0.955); lg->SetBorderSize(0); lg->SetTextSize(0.026); lg->SetFillStyle(0);
        lg->AddEntry(g1, "observed (1) N5 fitted, asymptotic CL_{s} (PRIMARY)", "lp");
        if (gt->GetN() > 0) lg->AddEntry(gt, "observed (1), toy CL_{s} (16.7, 16.85, 17.0)", "p");
        lg->AddEntry(fmed, "expected (1), FROZEN Asimov median", "l");
        lg->AddEntry(b1, "expected #pm 1#sigma (frozen)", "f"); lg->AddEntry(b2, "expected #pm 2#sigma (frozen)", "f");
        lg->AddEntry(tq, "expected (1), FROZEN toys: median, 16-84%", "pe");
        lg->AddEntry(g2, Form("observed (2) #eta yield a_{#eta} = %.2f (WP5)", c.z.a0), "l");
        lg->AddEntry(g3, Form("observed (3) #eta yield from the peak, a_{P} = %.2f", c.z.aP), "l");
        lg->AddEntry(gC, "observed (1), e^{+} decay angle 1 + cos^{2}#theta", "l");
        lg->AddEntry(gA, "Asimov median of the observed r = 0 fit (information)", "l");
        lg->Draw();
        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.030);
        tx.DrawLatex(0.18, 0.20, "S0, isBest #geq 0, X2D: M4 0.5225-0.610 GeV/c^{2} #times M(e^{+}e^{-}) 2-62 MeV/c^{2}");
        tx.DrawLatex(0.18, 0.155, "L = 6079.2 nb^{-1}; 2D Poisson likelihood, q~_{r}, frozen procedure");
        cv.SaveAs(Form("%s/x17d_X2_limit.pdf", kFigD)); cv.SaveAs(Form("%s/x17d_X2_limit.png", kFigD));
    }

    // ------------------------------------------------------------ Fig. X3: local p0 and r_hat vs m_X
    {
        TCanvas cv("cx3", "", 1000, 900);
        TPad* p1 = new TPad("p31", "", 0, 0.42, 1, 1), *p2 = new TPad("p32", "", 0, 0, 1, 0.42);
        p1->SetBottomMargin(0.02); p1->SetLeftMargin(0.13); p1->SetRightMargin(0.03); p1->SetTopMargin(0.04); p1->SetLogy();
        p2->SetTopMargin(0.02); p2->SetBottomMargin(0.22); p2->SetLeftMargin(0.13); p2->SetRightMargin(0.03);
        p1->Draw(); p2->Draw();
        p1->cd();
        TGraph* q1 = new TGraph(), *q2 = new TGraph(), *q3 = new TGraph();
        for (const auto& x : obs) { q1->SetPoint(q1->GetN(), x.m, x.p0[0]); q2->SetPoint(q2->GetN(), x.m, x.p0[1]); q3->SetPoint(q3->GetN(), x.m, x.p0[2]); }
        TH1D* fr1 = new TH1D("x3fr", ";m_{X} [MeV/c^{2}];local p_{0}", 10, 9.5, 40.5);
        fr1->SetMinimum(5e-4); fr1->SetMaximum(1.5); StyleAxis(fr1, 0.06, 0.055); fr1->GetXaxis()->SetLabelSize(0); fr1->GetYaxis()->SetTitleOffset(0.9);
        fr1->Draw("AXIS");
        for (int s = 1; s <= 3; ++s) {
            const double pz = 1 - x17c::Phi(s);
            TLine* l = new TLine(9.5, pz, 40.5, pz); l->SetLineColor(kRed + 1); l->SetLineStyle(2); l->Draw();
            TLatex* t = new TLatex(40.0, pz * 1.15, Form("%d#sigma", s)); t->SetTextSize(0.05); t->SetTextColor(kRed + 1); t->SetTextAlign(31); t->Draw();
        }
        q2->SetLineColor(kBlue + 1); q2->SetLineStyle(7); q2->SetLineWidth(2); q2->Draw("L SAME");
        q3->SetLineColor(kMagenta + 2); q3->SetLineStyle(9); q3->SetLineWidth(2); q3->Draw("L SAME");
        q1->SetLineWidth(3); q1->SetMarkerStyle(20); q1->Draw("LP SAME");
        TLegend* lg = new TLegend(0.15, 0.06, 0.62, 0.30); lg->SetBorderSize(0); lg->SetTextSize(0.045); lg->SetFillStyle(0);
        lg->AddEntry(q1, "(1) N5 fitted (PRIMARY)", "lp"); lg->AddEntry(q2, "(2) a_{#eta} normalisation", "l"); lg->AddEntry(q3, "(3) peak normalisation", "l");
        lg->Draw();
        double zmax = 0, mz = 0; for (const auto& x : obs) if (x.Z[0] > zmax) { zmax = x.Z[0]; mz = x.m; }
        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.048);
        tx.DrawLatex(0.52, 0.22, Form("largest local Z = %.2f at m_{X} = %.1f MeV/c^{2}", zmax, mz));
        p2->cd();
        TGraphErrors* gr = new TGraphErrors();
        for (const auto& x : obs) { const int k = gr->GetN(); gr->SetPoint(k, x.m, 1e5 * x.rhat[0]); gr->SetPointError(k, 0, 1e5 * x.rerr[0]); }
        TH1D* fr2 = new TH1D("x3fr2", ";m_{X} [MeV/c^{2}];fitted BR_{X}BR_{ee} [10^{-5}]", 10, 9.5, 40.5);
        fr2->SetMinimum(-6.5); fr2->SetMaximum(6.5); StyleAxis(fr2, 0.075, 0.07); fr2->GetYaxis()->SetTitleOffset(0.7); fr2->GetYaxis()->SetNdivisions(505);
        fr2->Draw("AXIS");
        gr->SetFillColorAlpha(kAzure - 4, 0.5); gr->Draw("3 SAME");
        TGraph* gl = new TGraph(); for (const auto& x : obs) gl->SetPoint(gl->GetN(), x.m, 1e5 * x.rhat[0]);
        gl->SetLineWidth(2); gl->Draw("L SAME");
        TLine* z = new TLine(9.5, 0, 40.5, 0); z->SetLineStyle(2); z->Draw();
        TLatex t2; t2.SetNDC(); t2.SetTextSize(0.065);
        t2.DrawLatex(0.16, 0.86, "normalisation (1): #hat{r} #times 2.68#times10^{-4} #pm Hessian error (band)");
        cv.SaveAs(Form("%s/x17d_X3_p0.pdf", kFigD)); cv.SaveAs(Form("%s/x17d_X3_p0.png", kFigD));
    }

    // ------------------------------------------------------------ Fig. X4: validation panels
    double pctUL = -1, pEmp = -1, ulFCup = NAN, ulFClo = NAN;
    {
        TCanvas cv("cx4", "", 1800, 1100);
        cv.Divide(3, 2, 0.004, 0.004);
        // (a) CLs vs BR at 17 MeV/c2: toys vs asymptotic
        cv.cd(1); gPad->SetLogy(); gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.13); gPad->SetRightMargin(0.03);
        {
            const TCRes& R = tcs[{1700, 1}];
            TGraphErrors* gT = new TGraphErrors(); TGraph* gA = new TGraph();
            for (const auto& p : R.pts) { const int k = gT->GetN(); gT->SetPoint(k, 1e5 * p.br, p.cls); gT->SetPointError(k, 0, p.err); gA->SetPoint(k, 1e5 * p.br, p.clsA); }
            TH1D* fr = new TH1D("x4a", "(a) m_{X} = 17 MeV/c^{2}, normalisation (1);BR_{X}BR_{ee} [10^{-5}];CL_{s} on the data", 10, 0, gT->GetN() ? 1e5 * R.pts.back().br * 1.05 : 6);
            fr->SetMinimum(0.01); fr->SetMaximum(1.0); StyleAxis(fr); fr->Draw("AXIS");
            { TLatex tt; tt.SetNDC(); tt.SetTextSize(0.045); tt.DrawLatex(0.15, 0.94, "(a) toy CL_{s} on the data, 17 MeV/c^{2}, (1)"); }
            gA->SetLineColor(kRed + 1); gA->SetLineWidth(2); gA->Draw("L SAME");
            gT->SetMarkerStyle(20); gT->Draw("P SAME");
            TLine* l = new TLine(0, 0.1, fr->GetXaxis()->GetXmax(), 0.1); l->SetLineStyle(2); l->Draw();
            TLegend* lg = new TLegend(0.35, 0.74, 0.97, 0.90); lg->SetBorderSize(0); lg->SetTextSize(0.040); lg->SetFillStyle(0);
            lg->AddEntry(gT, "toys (1000 per point)", "pe"); lg->AddEntry(gA, "asymptotic", "l"); lg->Draw();
            TLatex tx; tx.SetNDC(); tx.SetTextSize(0.042);
            if (R.ok) tx.DrawLatex(0.20, 0.24, Form("UL: asymptotic %.2f, toys %.2f #times10^{-5}", 1e5 * R.ulA, 1e5 * R.ulT));
            if (R.ok) tx.DrawLatex(0.20, 0.18, Form("toys / asymptotic = %.3f", R.ratio));
        }
        // (b) FC construction at 17 MeV/c2
        cv.cd(2); gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.13); gPad->SetRightMargin(0.03);
        {
            TGraph* go = new TGraph(), *g90 = new TGraph(), *g68 = new TGraph();
            for (const auto& p : fc) { go->SetPoint(go->GetN(), 1e5 * p.br, p.tobs); g90->SetPoint(g90->GetN(), 1e5 * p.br, p.t90); g68->SetPoint(g68->GetN(), 1e5 * p.br, p.t68); }
            // interval: t_obs <= t90 (linear interpolation of the difference)
            double lo = NAN, hi = NAN;
            for (size_t k = 0; k < fc.size(); ++k) {
                const double d = fc[k].tobs - fc[k].t90;
                if (d <= 0 && std::isnan(lo)) {
                    if (k == 0) lo = fc[k].br;
                    else { const double d0 = fc[k - 1].tobs - fc[k - 1].t90; lo = fc[k - 1].br + (fc[k].br - fc[k - 1].br) * d0 / (d0 - d); }
                }
                if (d <= 0) hi = fc[k].br;
                if (d > 0 && !std::isnan(lo) && k > 0 && fc[k - 1].tobs - fc[k - 1].t90 <= 0) {
                    const double d0 = fc[k - 1].tobs - fc[k - 1].t90;
                    hi = fc[k - 1].br + (fc[k].br - fc[k - 1].br) * d0 / (d0 - d);
                    break;
                }
            }
            ulFClo = lo; ulFCup = hi;
            TH1D* fr = new TH1D("x4b", "(b) Feldman-Cousins, 17 MeV/c^{2}, (1);BR_{X}BR_{ee} [10^{-5}];t_{r} = -2 ln #lambda(r) / qScale", 10, 0,
                                fc.empty() ? 5 : 1e5 * fc.back().br * 1.03);
            double ym = 0; for (const auto& p : fc) ym = std::max(ym, std::max(p.tobs, p.t90));
            fr->SetMinimum(0); fr->SetMaximum(1.6 * ym); StyleAxis(fr); fr->Draw("AXIS");
            { TLatex tt; tt.SetNDC(); tt.SetTextSize(0.045); tt.DrawLatex(0.15, 0.94, "(b) Feldman-Cousins, 17 MeV/c^{2}, (1)"); }
            if (std::isfinite(lo) && std::isfinite(hi)) { TBox* b = new TBox(1e5 * lo, 0, 1e5 * hi, 1.6 * ym); b->SetFillColorAlpha(kGreen - 9, 0.5); b->Draw(); }
            g90->SetLineColor(kRed + 1); g90->SetLineWidth(2); g90->SetMarkerStyle(21); g90->SetMarkerColor(kRed + 1); g90->Draw("LP SAME");
            g68->SetLineColor(kOrange + 1); g68->SetLineWidth(2); g68->SetLineStyle(2); g68->Draw("L SAME");
            go->SetLineWidth(3); go->SetMarkerStyle(20); go->Draw("LP SAME");
            TLegend* lg = new TLegend(0.18, 0.70, 0.70, 0.90); lg->SetBorderSize(0); lg->SetTextSize(0.040); lg->SetFillStyle(0);
            lg->AddEntry(go, "t_{obs}(r), data", "lp"); lg->AddEntry(g90, "t_{c}(r): 90% quantile of toys", "lp"); lg->AddEntry(g68, "68% quantile", "l");
            lg->Draw();
            TLatex tx; tx.SetNDC(); tx.SetTextSize(0.042);
            if (std::isfinite(hi)) tx.DrawLatex(0.20, 0.66, Form("90%% CL interval: [%.2f, %.2f] #times10^{-5}", 1e5 * lo, 1e5 * hi));
        }
        // (c) observed UL in the frozen background-only toy distribution (17 MeV/c2)
        cv.cd(3); gPad->SetLogx(); gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.13); gPad->SetRightMargin(0.03);
        {
            const ToyBO* T = boAt(17.0);
            const int nb = 50;
            std::vector<double> be(nb + 1);
            for (int k = 0; k <= nb; ++k) be[k] = 5e-6 * std::pow(40., (double)k / nb);
            TH1D* h = new TH1D("x4c", "(c) 17 MeV/c^{2}: frozen background-only toys;90% CL_{s} UL on BR_{X}BR_{ee} (1);toys per bin", nb, be.data());
            if (T) for (double v : T->ul[0]) h->Fill(v);
            StyleAxis(h); h->SetLineWidth(2); h->SetLineColor(kGray + 3); h->SetFillColor(kGray); h->SetMaximum(h->GetMaximum() * 1.3);
            h->Draw("HIST");
            if (T && o17) { int nb2 = 0; for (double v : T->ul[0]) if (v <= o17->ul[0]) ++nb2; pctUL = (double)nb2 / T->ul[0].size(); }
            if (o17) { TLine* l = new TLine(o17->ul[0], 0, o17->ul[0], h->GetMaximum()); l->SetLineColor(kRed + 1); l->SetLineWidth(3); l->Draw(); }
            TLatex tx; tx.SetNDC(); tx.SetTextSize(0.042);
            if (T) tx.DrawLatex(0.18, 0.84, Form("%zu toys; median %.2f#times10^{-5}", T->ul[0].size(), 1e5 * x17d::Quantile(T->ul[0], 0.5)));
            tx.SetTextColor(kRed + 1);
            if (o17) tx.DrawLatex(0.18, 0.78, Form("observed %.2f#times10^{-5}: %.0f%% of toys lower", 1e5 * o17->ul[0], 100 * pctUL));
        }
        // (d) r_hat / sigma: observed vs background-only toys (17 MeV/c2) -> empirical p0
        cv.cd(4); gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.13); gPad->SetRightMargin(0.03);
        {
            const ToyBO* T = boAt(17.0);
            TH1D* h = new TH1D("x4d", "(d) 17 MeV/c^{2}: #hat{r}/#sigma_{#hat{r}}, background-only toys;#hat{r}/#sigma_{#hat{r}} (normalisation 1);toys per 0.2", 40, -4, 4);
            if (T) for (double v : T->pull) h->Fill(v);
            StyleAxis(h); h->SetLineWidth(2); h->SetLineColor(kGray + 3); h->SetFillColor(kGray); h->SetMaximum(h->GetMaximum() * 1.3); h->Draw("HIST");
            const double zo = o17 ? o17->rhat[0] / o17->rerr[0] : 0;
            if (T) { int k = 0; for (double v : T->pull) if (v >= zo) ++k; pEmp = (double)k / T->pull.size(); }
            TLine* l = new TLine(zo, 0, zo, h->GetMaximum()); l->SetLineColor(kRed + 1); l->SetLineWidth(3); l->Draw();
            TLatex tx; tx.SetNDC(); tx.SetTextSize(0.040); tx.SetTextColor(kRed + 1);
            if (o17) tx.DrawLatex(0.18, 0.84, Form("observed %.2f", zo));
            if (o17) tx.DrawLatex(0.18, 0.78, Form("p (toys) %.3f, p_{0} (asymptotic) %.3f", pEmp, o17->p0[0]));
        }
        // (e) template statistics: UL over CB and C4 replicas (17 MeV/c2)
        cv.cd(5); gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.13); gPad->SetRightMargin(0.03);
        {
            const auto& v = mcs[1700];
            double nom = 0; std::vector<double> a, b;
            for (const auto& r : v) { if (r.kind == 0) nom = r.ul[0]; else if (r.kind == 1) a.push_back(r.ul[0]); else b.push_back(r.ul[0]); }
            const double lo = nom > 0 ? nom * 0.70 : 2e-5, hi = nom > 0 ? nom * 1.20 : 4e-5;
            TH1D* ha = new TH1D("x4e1", "(e) 17 MeV/c^{2}: template replicas;observed 90% CL_{s} UL (1) [10^{-5}];replicas per bin", 40, 1e5 * lo, 1e5 * hi);
            TH1D* hb = new TH1D("x4e2", "", 40, 1e5 * lo, 1e5 * hi);
            for (double x : a) ha->Fill(1e5 * x);
            for (double x : b) hb->Fill(1e5 * x);
            StyleAxis(ha); ha->SetLineColor(kRed + 1); ha->SetLineWidth(2); hb->SetLineColor(kAzure + 2); hb->SetLineWidth(2); hb->SetLineStyle(2);
            ha->SetMaximum(1.4 * std::max(ha->GetMaximum(), hb->GetMaximum()));
            ha->Draw("HIST"); hb->Draw("HIST SAME");
            if (nom > 0) { TLine* l = new TLine(1e5 * nom, 0, 1e5 * nom, ha->GetMaximum()); l->SetLineWidth(3); l->Draw(); }
            TLegend* lg = new TLegend(0.18, 0.70, 0.97, 0.90); lg->SetBorderSize(0); lg->SetTextSize(0.040); lg->SetFillStyle(0);
            lg->AddEntry(ha, Form("CB replicas (%zu): RMS %.2f#times10^{-5}", a.size(), 1e5 * Rms(a)), "l");
            lg->AddEntry(hb, Form("C4 replicas (%zu): RMS %.2f#times10^{-5}", b.size(), 1e5 * Rms(b)), "l");
            if (nom > 0) lg->AddEntry((TObject*)nullptr, Form("black line: nominal %.2f#times10^{-5}", 1e5 * nom), "");
            lg->Draw();
        }
        // (f) formerly masked rows: fitted n_j against the frozen prediction
        cv.cd(6); gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.13); gPad->SetRightMargin(0.03);
        {
            TVectorD* pz = (TVectorD*)fs->Get("pZero");
            TVectorD* ez = (TVectorD*)fs->Get("eZero");
            TVectorD* pf = (TVectorD*)fs->Get("pSBfill");
            TGraphErrors* gfit = new TGraphErrors(), *gmask = new TGraphErrors();
            for (int j = 0; j < ny; ++j) {
                const double y = 1e3 * I.G.yc(j), rat = (*pz)(j) / (*pf)(j), er = (*ez)(j) / (*pf)(j);
                TGraphErrors* g = (j >= kMaskLo && j <= kMaskHi) ? gmask : gfit;
                const int k = g->GetN(); g->SetPoint(k, y, rat); g->SetPointError(k, 0, er);
            }
            TH1D* fr = new TH1D("x4f", "(f) other background per M(ee) row;M(e^{+}e^{-}) [MeV/c^{2}];n_{j}(full fit) / n_{j}(frozen side-band fit)", 10, 0, 64);
            fr->SetMinimum(0.7); fr->SetMaximum(1.3); StyleAxis(fr); fr->GetYaxis()->SetTitleSize(0.042); fr->Draw("AXIS");
            { TLatex tt; tt.SetNDC(); tt.SetTextSize(0.045); tt.DrawLatex(0.15, 0.94, "(f) other background n_{j} per M(e^{+}e^{-}) row"); }
            TBox* b = new TBox(12, 0.7, 28, 1.3); b->SetFillColorAlpha(kYellow - 9, 0.4); b->Draw();
            TLine* l = new TLine(0, 1, 64, 1); l->SetLineStyle(2); l->Draw();
            gfit->SetMarkerStyle(24); gfit->SetMarkerColor(kGray + 2); gfit->SetLineColor(kGray + 2); gfit->Draw("P SAME");
            gmask->SetMarkerStyle(20); gmask->SetMarkerColor(kBlue + 1); gmask->SetLineColor(kBlue + 1); gmask->Draw("P SAME");
            // [vx verification: legend moved to the empty lower right; its markers sat at y 1.22-1.27 inside the 12-28 MeV band and
            //  looked like data points]
            TLegend* lg = new TLegend(0.52, 0.15, 0.97, 0.29); lg->SetBorderSize(0); lg->SetTextSize(0.034); lg->SetFillStyle(0);
            lg->AddEntry(gmask, "12-28 MeV: fit / prediction", "pe");
            lg->AddEntry(gfit, "side bands: full / side-band fit", "pe");
            lg->Draw();
        }
        cv.SaveAs(Form("%s/x17d_X4_valid.pdf", kFigD)); cv.SaveAs(Form("%s/x17d_X4_valid.png", kFigD));
    }

    // ------------------------------------------------------------ Table D1 and the checks table
    const double sigEtaProv = c.z.a0 * alldefs::kSigEtaSMASH * 1e6 * kCsecFwin;   // nb
    FILE* tb = fopen(Form("%s/x17d_D1.txt", kTabD), "w");
    fprintf(tb, "# x17d_D1.txt -- Table D1: OBSERVED 90 %% CLs upper limits on BR(eta -> pi+pi-X17) BR(X17 -> e+e-) (x17d_unblind.C, x17d_toycls.C, x17d_mc.C,\n");
    fprintf(tb, "# x17d_figs.C). Frozen procedure of x17_frozen.txt; S0, isBest >= 0, X2D (M4 0.5225-0.6100 GeV/c2 x M(ee) 2-62 MeV/c2), L = 6079.2 nb^-1; isotropic e+ angle.\n");
    fprintf(tb, "# (1) r x 2.68e-4, N5 fitted (PRIMARY); (2) N_X / (a_eta s_X2D), a_eta = %.4f (WP5); (3) relative to Y_peak: a_P = %.4f.\n", c.z.a0, c.z.aP);
    fprintf(tb, "# exp = FROZEN expectation (2000 background-only toys: median [16 %%, 84 %%]); UL_toy = toy CLs on the data (1000 toys per point and hypothesis);\n");
    fprintf(tb, "# N_X_UL = direct 90 %% CLs UL on the number of X17 rows (raw rows, X2D M4 range, all M(ee); POI = rows, no yield nuisance; the same for all\n");
    fprintf(tb, "#   normalisations) [vx verification]; rUL_N5_RX = the BR limit converted back with the r = 0 fit (r_UL N5 R_X, resp. r_UL a N5ref R_X): NOT a limit on rows;\n");
    fprintf(tb, "# MC = RMS of the UL over 200 CB + 100 C4 template replicas (quoted separately);\n");
    fprintf(tb, "# cos2 = UL (1) with e+ decay angle 1 + cos^2 (model dependence, quoted separately); sigma x BR x BR = UL(1) x sigma(pp -> eta X) PROVISIONAL\n");
    fprintf(tb, "# (WP8 not available): sigma(pp -> eta X) = a_eta x 1.0665 mb x C_sec (1 - f_eta) = %.4f x 1.0665 x %.3f = %.4f mb (spec amendment 1; +3/-4 %%).\n", c.z.a0, kCsecFwin, sigEtaProv * 1e-6);
    fprintf(tb, "# MC_RMS = RMS over the 200 CB replicas (+) RMS over the 100 C4 replicas, in quadrature [vx verification: was the RMS of the pooled 300 values].\n");
    fprintf(tb, "%-6s %-4s %10s %10s %7s %26s %13s %6s %6s %8s %9s %10s %10s %9s\n", "m", "norm", "UL_asym", "UL_toy", "toy/as", "exp median [16%, 84%]", "r_hat(BR)", "Z",
            "p0", "N_X_UL", "MC_RMS", "UL_cos2", "sBRBR[nb]", "rUL_N5_RX");
    struct D1row { double m; int norm; double ulA, ulT, ratio, em, el, eh, rh, re, Z, p0, nX, mc, c2, sb, nXconv; };
    std::vector<D1row> d1;
    for (double m : kQuoted) {
        const Obs* x = obsAt(m);
        const ToyBO* T = boAt(m);
        const auto& v = mcs[(int)std::lround(100 * m)];
        for (int a = 0; a < 3; ++a) {
            D1row r{};
            r.m = m; r.norm = a + 1;
            r.ulA = x ? x->ul[a] : NAN;
            const TCRes& R = tcs[{(int)std::lround(100 * m), a + 1}];
            r.ulT = R.ok ? R.ulT : NAN; r.ratio = R.ok ? R.ratio : NAN;
            r.em = T ? x17d::Quantile(T->ul[a], 0.5) : NAN; r.el = T ? x17d::Quantile(T->ul[a], kQ[1]) : NAN; r.eh = T ? x17d::Quantile(T->ul[a], kQ[3]) : NAN;
            r.rh = x ? x->rhat[a] : NAN; r.re = x ? x->rerr[a] : NAN; r.Z = x ? x->Z[a] : NAN; r.p0 = x ? x->p0[a] : NAN; r.nXconv = x ? x->nX[a] : NAN; r.nX = nXrows.count((int)std::lround(100 * m)) ? nXrows[(int)std::lround(100 * m)] : NAN;
            // [vx verification: CB and C4 replicas are independent sources -> RMS in quadrature (the pooled RMS of the 300 values,
            //  used before, mixes the two samples and gave 2.4e-6 instead of 3.4e-6 at 17 MeV/c2, normalisation 1)]
            std::vector<double> u1, u2; for (const auto& q : v) { if (q.kind == 1) u1.push_back(q.ul[a]); if (q.kind == 2) u2.push_back(q.ul[a]); }
            r.mc = std::hypot(Rms(u1), Rms(u2));
            r.c2 = (a == 0 && x) ? x->ulC2 : NAN;
            r.sb = (a == 0 && x) ? x->ul[0] * sigEtaProv : NAN;
            d1.push_back(r);
            fprintf(tb, "%6.2f %4d %10.3e %10.3e %7.3f %10.3e [%.2e, %.2e] %+.2e+-%.2e %6.2f %6.3f %8.1f %9.2e %10.3e %10.2f %9.1f\n", m, a + 1, r.ulA, r.ulT, r.ratio, r.em, r.el,
                    r.eh, r.rh, r.re, r.Z, r.p0, r.nX, r.mc, r.c2, r.sb, r.nXconv);
        }
    }
    const double nXc = (*misc)(5), exO = (*misc)(8), erO = (*misc)(9), brO = (*misc)(10), bseO = (*misc)(11);
    fprintf(tb, "counting X17c (M4 0.5312-0.5566, M(ee) 15.5-24.5): N %.1f (overlap; exact %.0f), B_pred %.1f +- %.1f (frozen), excess %+.1f +- %.1f rows, BR = %+.2e +- %.2e\n",
            nXc, (*misc)(6), kBpredX17c, kSigBX17c, exO, erO, brO, bseO);
    if (std::isfinite(ulFCup)) fprintf(tb, "Feldman-Cousins 90 %% CL interval at 17 MeV/c2 (normalisation 1): [%.3e, %.3e]\n", ulFClo, ulFCup);
    fprintf(tb, "observed UL at 17 MeV/c2 in the frozen background-only toy distribution: %.1f %% of the toys have a lower UL; r_hat/sigma empirical p %.3f\n", 100 * pctUL, pEmp);
    fprintf(tb, "simulation expectation of the simulation study (comparison only): 1.9e-5 / 5.7e-6\n");
    fclose(tb);
    // LaTeX snippet
    tb = fopen(Form("%s/x17d_D1.tex", kTabD), "w");
    fprintf(tb, "%% x17d_D1.tex -- Table D1 (x17d_figs.C). Observed 90%% CLs upper limits on BR(eta->pi+pi-X17)BR(X17->e+e-), frozen procedure.\n");
    fprintf(tb, "\\begin{tabular}{llcccccc}\n\\toprule\n");
    fprintf(tb, "$m_X$ [MeV/$c^2$] & norm. & observed UL & toy-CL$_s$ UL & expected (median, 68\\%%) & $\\hat{r}\\times2.68\\cdot10^{-4}$ & local $Z$ & $N_X^{\\rm UL}$ \\\\\n\\midrule\n");
    for (const auto& r : d1) {
        fprintf(tb, "%.2f & (%d) & $%.2f\\cdot10^{-5}$ & $%.2f\\cdot10^{-5}$ & $%.2f\\,[%.2f, %.2f]\\cdot10^{-5}$ & $(%+.2f\\pm%.2f)\\cdot10^{-5}$ & %.2f & %.0f \\\\\n", r.m, r.norm,
                1e5 * r.ulA, 1e5 * r.ulT, 1e5 * r.em, 1e5 * r.el, 1e5 * r.eh, 1e5 * r.rh, 1e5 * r.re, r.Z, r.nX);
        if (r.norm == 3 && r.m < 16.99) fprintf(tb, "\\addlinespace\n");
    }
    fprintf(tb, "\\bottomrule\n\\end{tabular}\n");
    fclose(tb);
    // checks
    tb = fopen(Form("%s/x17d_checks.txt", kTabD), "w");
    fprintf(tb, "# x17d_checks.txt -- X17 unblinding checks (x17d_figs.C)\n");
    fprintf(tb, "toy CLs vs asymptotic (data):\n");
    for (double m : kQuoted) for (int a = 1; a <= 3; ++a) {
        const TCRes& R = tcs[{(int)std::lround(100 * m), a}];
        fprintf(tb, "  m %.2f norm %d: UL asym %.4e toy %.4e ratio %.4f %s%s\n", m, a, R.ulA, R.ulT, R.ratio, R.ok ? "" : "[missing]", R.clb0 ? " CL_b = 0 at some point" : "");
    }
    fprintf(tb, "FC at 17 MeV/c2 (norm 1): r, BR, t_obs, t68, t90, t95, p(t >= t_obs), toys\n");
    for (const auto& p : fc) fprintf(tb, "  %.4f %.4e %.4f %.4f %.4f %.4f %.4f %d/%d\n", p.r, p.br, p.tobs, p.t68, p.t90, p.t95, p.p, p.nok, p.nt);
    fprintf(tb, "FC 90 %% CL interval: [%.4e, %.4e]\n", ulFClo, ulFCup);
    for (double m : kQuoted) {
        const auto& v = mcs[(int)std::lround(100 * m)];
        for (int a = 0; a < 3; ++a) {
            std::vector<double> u1, u2, r1, r2; double nom = 0;
            for (const auto& q : v) { if (q.kind == 0) nom = q.ul[a]; if (q.kind == 1) { u1.push_back(q.ul[a]); r1.push_back(q.rh[a]); } if (q.kind == 2) { u2.push_back(q.ul[a]); r2.push_back(q.rh[a]); } }
            std::vector<double> u = u1; u.insert(u.end(), u2.begin(), u2.end());
            fprintf(tb, "MC m %.2f norm %d: nominal UL %.4e; CB replicas %zu: UL mean %.4e RMS %.3e, r_hat RMS %.3e; C4 replicas %zu: UL mean %.4e RMS %.3e, r_hat RMS %.3e; combined RMS %.3e (%.1f %%)\n",
                    m, a + 1, nom, u1.size(), Mean(u1), Rms(u1), Rms(r1), u2.size(), Mean(u2), Rms(u2), Rms(r2), std::hypot(Rms(u1), Rms(u2)), nom > 0 ? 100 * std::hypot(Rms(u1), Rms(u2)) / nom : 0.);
        }
    }
    fprintf(tb, "observed UL(1) at 17 MeV/c2 percentile in the frozen toys: %.3f; empirical p of r_hat/sigma: %.3f\n", pctUL, pEmp);
    for (const auto& T : bo) {
        const Obs* x = obsAt(T.m);
        if (!x) continue;
        int k = 0; for (double v : T.ul[0]) if (v <= x->ul[0]) ++k;
        int kp = 0; const double zo = x->rhat[0] / x->rerr[0]; for (double v : T.pull) if (v >= zo) ++kp;
        fprintf(tb, "frozen toys m %.2f: %zu toys; observed UL(1) %.3e at percentile %.3f (median %.3e); r_hat/sigma %.3f, empirical p %.3f, asymptotic p0 %.3f\n", T.m,
                T.ul[0].size(), x->ul[0], (double)k / T.ul[0].size(), x17d::Quantile(T.ul[0], 0.5), zo, (double)kp / T.pull.size(), x->p0[0]);
    }
    fclose(tb);
    printf("x17d_figs: done (FC interval [%.3e, %.3e], UL percentile %.3f, empirical p %.3f)\n", ulFClo, ulFCup, pctUL, pEmp);
    fs->Close(); fz->Close(); c.f->Close();
}
