// ============================================================================
// x17p26_figs.C  --  X17 limit with PDG 2026 (analyst ruling D25, 15 November 2026), prefix x17p26, written 2026-11-15.
// NEW figures (suffix _pdg26; the frozen ones are not touched), drawn with the code blocks of x17d_figs.C (Figs. X3, X4) and
// x17c_summary.C (x17c_expected, x17c_toys_ul, x17c_inject), copied verbatim except for the PDG 2026 factors and the output names:
//   normalisation (1): asymptotic values from the rerun of the frozen code with ONLY BR_dir changed (x17p26_limit.C ->
//     x17p26_state.root); toys, FC and template replicas: frozen x 2.76/2.68 x k (k = rerun/frozen at that mass: observed ratio for
//     the data toys, FC and replicas, Asimov ratio for the expected toys); r_hat x 2.76/2.68; injected BR x 2.76/2.68;
//   normalisation (2) x 1.00436, (3) x 1.0224 (WP5 eta scales fitted with PDG 2020 BRs; the rerun confirms the linearity);
//   the pre-data projection of x17ab (simulation Asimov) x 2.76/2.68.
// Assembled by scratchpad/x17p26/assemble_figs.py from x17d_figs.C and x17c_summary.C (blocks copied, replacements asserted).
// Run: root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_x17p26", kTRUE)' 'plots/all_report/x17p26_figs.C+O'
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

namespace x17p26f {
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
}  // namespace x17p26f
namespace x17cs {
using namespace x17cl;
struct ToySet { double m = 0, inj = 0; std::vector<double> ul[3], rhat[3], rerr[3]; int nTot = 0, nOk[3] = {0, 0, 0}; };
bool ReadToys(double m, double inj, const std::vector<int>& firsts, ToySet& T) {
    T.m = m; T.inj = inj;
    for (int fst : firsts) {
        TFile* f = TFile::Open(Form("%s/toys_m%04d_i%04d_%04d.root", kOutDir, (int)std::lround(100 * m), (int)std::lround(inj * 1e7), fst), "READ");
        if (!f || f->IsZombie()) { printf("x17c_summary: toys m %.2f inj %.1e first %d missing\n", m, inj, fst); return false; }
        TTree* t = (TTree*)f->Get("toys");
        int ok[3]; double ul[3], rh[3], re[3];
        t->SetBranchAddress("ok", ok); t->SetBranchAddress("ul", ul); t->SetBranchAddress("rhat", rh); t->SetBranchAddress("rerr", re);
        for (Long64_t k = 0; k < t->GetEntries(); ++k) {
            t->GetEntry(k);
            ++T.nTot;
            for (int a = 0; a < 3; ++a) {
                if (!ok[a] || !std::isfinite(ul[a])) continue;
                ++T.nOk[a];
                T.ul[a].push_back(ul[a]); T.rhat[a].push_back(rh[a]); T.rerr[a].push_back(re[a]);
            }
        }
        f->Close();
    }
    return true;
}
double Quant(std::vector<double> v, double q) {
    if (v.empty()) return (double)NAN;
    std::sort(v.begin(), v.end());
    const double x = q * (v.size() - 1);
    const size_t i = (size_t)std::floor(x);
    const double w = x - i;
    return i + 1 < v.size() ? (1 - w) * v[i] + w * v[i + 1] : v.back();
}
const double kQ[5] = {0.02275, 0.15866, 0.5, 0.84134, 0.97725};
}  // namespace x17cs

void x17p26_figs() {
    using namespace x17p26f;
    gStyle->SetOptStat(0);
    gStyle->SetEndErrorSize(0);
    const double FBR = 2.76 / 2.68, F2 = 1.00436, F3 = 1.0224;
    // ------------------------------------------------------------ PDG 2026 rerun (x17p26_limit.C)
    TFile* f26 = TFile::Open("/hdd3/przygoda/HADES/FAT_out/all_report/x17p26/x17p26_state.root", "READ");
    TTree* t26 = (TTree*)f26->Get("obs26");
    double m26, ul26[3], e126[5], e226[5], e326[5], e120[5], r20, r26;
    t26->SetBranchAddress("m", &m26); t26->SetBranchAddress("ul", ul26); t26->SetBranchAddress("e1", e126); t26->SetBranchAddress("e2", e226);
    t26->SetBranchAddress("e3", e326); t26->SetBranchAddress("e1_20", e120); t26->SetBranchAddress("rUL20", &r20); t26->SetBranchAddress("rUL26", &r26);
    std::map<int, double> kObs, kAsi, ulP26;
    std::vector<double> am, ae[3][5];
    for (Long64_t k = 0; k < t26->GetEntries(); ++k) {
        t26->GetEntry(k);
        const int key = (int)std::lround(100 * m26);
        kObs[key] = r26 / r20; kAsi[key] = e126[2] / (e120[2] * FBR); ulP26[key] = ul26[0];
        am.push_back(m26);
        for (int q = 0; q < 5; ++q) { ae[0][q].push_back(e126[q]); ae[1][q].push_back(e226[q]); ae[2][q].push_back(e326[q]); }
    }
    const double F1o17 = FBR * kObs[1700], F1a17 = FBR * kAsi[1700], ul26_17 = ulP26[1700];
    printf("x17p26_figs: F1(obs, 17) %.5f, F1(Asimov, 17) %.5f, UL(1) 17 %.4e\n", F1o17, F1a17, ul26_17);
    auto asimovAt = [&](double mm, int norm, int q) -> double { for (size_t k = 0; k < am.size(); ++k) if (std::fabs(am[k] - mm) < 1e-6) return ae[norm][q][k]; return (double)NAN; };
    // ------------------------------------------------------------ frozen inputs (as x17d_figs.C)
    Ctx c;
    if (!Init(c)) { printf("x17p26_figs: init failed\n"); return; }
    x17c::Inputs& I = c.I;
    const int ny = I.G.ny;
    TFile* fs = TFile::Open(Form("%s/x17d_state.root", kOutD), "READ");
    TTree* to = (TTree*)fs->Get("obs");
    Obs o;
    to->SetBranchAddress("m", &o.m); to->SetBranchAddress("RX", &o.RX); to->SetBranchAddress("ul", o.ul); to->SetBranchAddress("rhat", o.rhat);
    to->SetBranchAddress("rerr", o.rerr); to->SetBranchAddress("expA", o.ea); to->SetBranchAddress("q0", o.q0); to->SetBranchAddress("Z", o.Z);
    to->SetBranchAddress("p0", o.p0); to->SetBranchAddress("nX", o.nX); to->SetBranchAddress("N5", &o.N5); to->SetBranchAddress("ulCos2", &o.ulC2);
    to->SetBranchAddress("s1x", &o.s1x);
    std::vector<Obs> obs;
    for (Long64_t k = 0; k < to->GetEntries(); ++k) { to->GetEntry(k); obs.push_back(o); }
    auto obsAt = [&](double m) -> const Obs* { for (const auto& x : obs) if (std::fabs(x.m - m) < 1e-6) return &x; return nullptr; };
    std::vector<ToyBO> bo;
    for (double m : {16.7, 16.85, 17.0}) { ToyBO T; if (ReadToysBO(m, {0, 500, 1000, 1500}, T)) bo.push_back(T); }
    auto boAt = [&](double m) -> const ToyBO* { for (const auto& x : bo) if (std::fabs(x.m - m) < 1e-6) return &x; return nullptr; };
    std::map<std::pair<int, int>, TCRes> tcs;
    tcs[{1700, 1}] = ReadTC(17.0, 1);
    const std::vector<FC> fc = ReadFC(17.0, 1);
    std::map<int, std::vector<MCr>> mcs;
    mcs[1700] = ReadMC(17.0);
    const Obs* o17 = obsAt(17.0);

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
        for (const auto& x : obs) { const int k = gr->GetN(); gr->SetPoint(k, x.m, 1e5 * FBR * x.rhat[0]); gr->SetPointError(k, 0, 1e5 * FBR * x.rerr[0]); }
        TH1D* fr2 = new TH1D("x3fr2", ";m_{X} [MeV/c^{2}];fitted BR_{X}BR_{ee} [10^{-5}]", 10, 9.5, 40.5);
        fr2->SetMinimum(-6.5); fr2->SetMaximum(6.5); StyleAxis(fr2, 0.075, 0.07); fr2->GetYaxis()->SetTitleOffset(0.7); fr2->GetYaxis()->SetNdivisions(505);
        fr2->Draw("AXIS");
        gr->SetFillColorAlpha(kAzure - 4, 0.5); gr->Draw("3 SAME");
        TGraph* gl = new TGraph(); for (const auto& x : obs) gl->SetPoint(gl->GetN(), x.m, 1e5 * FBR * x.rhat[0]);
        gl->SetLineWidth(2); gl->Draw("L SAME");
        TLine* z = new TLine(9.5, 0, 40.5, 0); z->SetLineStyle(2); z->Draw();
        TLatex t2; t2.SetNDC(); t2.SetTextSize(0.065);
        t2.DrawLatex(0.16, 0.86, "normalisation (1): #hat{r} #times 2.76#times10^{-4} (PDG 2026) #pm Hessian error (band)");
        cv.SaveAs(Form("%s/x17d_X3_p0_pdg26.pdf", kFigD)); cv.SaveAs(Form("%s/x17d_X3_p0_pdg26.png", kFigD));
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
            for (const auto& p : R.pts) { const int k = gT->GetN(); gT->SetPoint(k, 1e5 * F1o17 * p.br, p.cls); gT->SetPointError(k, 0, p.err); gA->SetPoint(k, 1e5 * F1o17 * p.br, p.clsA); }
            TH1D* fr = new TH1D("x4a", "(a) m_{X} = 17 MeV/c^{2}, normalisation (1);BR_{X}BR_{ee} [10^{-5}];CL_{s} on the data", 10, 0, gT->GetN() ? 1e5 * F1o17 * R.pts.back().br * 1.05 : 6);
            fr->SetMinimum(0.01); fr->SetMaximum(1.0); StyleAxis(fr); fr->Draw("AXIS");
            { TLatex tt; tt.SetNDC(); tt.SetTextSize(0.045); tt.DrawLatex(0.15, 0.94, "(a) toy CL_{s} on the data, 17 MeV/c^{2}, (1)"); }
            gA->SetLineColor(kRed + 1); gA->SetLineWidth(2); gA->Draw("L SAME");
            gT->SetMarkerStyle(20); gT->Draw("P SAME");
            TLine* l = new TLine(0, 0.1, fr->GetXaxis()->GetXmax(), 0.1); l->SetLineStyle(2); l->Draw();
            TLegend* lg = new TLegend(0.35, 0.74, 0.97, 0.90); lg->SetBorderSize(0); lg->SetTextSize(0.040); lg->SetFillStyle(0);
            lg->AddEntry(gT, "toys (1000 per point)", "pe"); lg->AddEntry(gA, "asymptotic", "l"); lg->Draw();
            TLatex tx; tx.SetNDC(); tx.SetTextSize(0.042);
            if (R.ok) tx.DrawLatex(0.20, 0.24, Form("UL: asymptotic %.2f, toys %.2f #times10^{-5}", 1e5 * ul26_17, 1e5 * F1o17 * R.ulT));
            if (R.ok) tx.DrawLatex(0.20, 0.18, Form("toys / asymptotic = %.3f", F1o17 * R.ulT / ul26_17));
        }
        // (b) FC construction at 17 MeV/c2
        cv.cd(2); gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.13); gPad->SetRightMargin(0.03);
        {
            TGraph* go = new TGraph(), *g90 = new TGraph(), *g68 = new TGraph();
            for (const auto& p : fc) { go->SetPoint(go->GetN(), 1e5 * F1o17 * p.br, p.tobs); g90->SetPoint(g90->GetN(), 1e5 * F1o17 * p.br, p.t90); g68->SetPoint(g68->GetN(), 1e5 * F1o17 * p.br, p.t68); }
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
            lo *= F1o17; hi *= F1o17; ulFClo = lo; ulFCup = hi;
            TH1D* fr = new TH1D("x4b", "(b) Feldman-Cousins, 17 MeV/c^{2}, (1);BR_{X}BR_{ee} [10^{-5}];t_{r} = -2 ln #lambda(r) / qScale", 10, 0,
                                fc.empty() ? 5 : 1e5 * F1o17 * fc.back().br * 1.03);
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
            if (T) for (double v : T->ul[0]) h->Fill(F1a17 * v);
            StyleAxis(h); h->SetLineWidth(2); h->SetLineColor(kGray + 3); h->SetFillColor(kGray); h->SetMaximum(h->GetMaximum() * 1.3);
            h->Draw("HIST");
            if (T && o17) { int nb2 = 0; for (double v : T->ul[0]) if (F1a17 * v <= ul26_17) ++nb2; pctUL = (double)nb2 / T->ul[0].size(); }
            if (o17) { TLine* l = new TLine(ul26_17, 0, ul26_17, h->GetMaximum()); l->SetLineColor(kRed + 1); l->SetLineWidth(3); l->Draw(); }
            TLatex tx; tx.SetNDC(); tx.SetTextSize(0.042);
            if (T) tx.DrawLatex(0.18, 0.84, Form("%zu toys; median %.2f#times10^{-5}", T->ul[0].size(), 1e5 * F1a17 * x17d::Quantile(T->ul[0], 0.5)));
            tx.SetTextColor(kRed + 1);
            if (o17) tx.DrawLatex(0.18, 0.78, Form("observed %.2f#times10^{-5}: %.0f%% of toys lower", 1e5 * ul26_17, 100 * pctUL));
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
            for (const auto& r : v) { if (r.kind == 0) nom = F1o17 * r.ul[0]; else if (r.kind == 1) a.push_back(F1o17 * r.ul[0]); else b.push_back(F1o17 * r.ul[0]); }
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
        cv.SaveAs(Form("%s/x17d_X4_valid_pdg26.pdf", kFigD)); cv.SaveAs(Form("%s/x17d_X4_valid_pdg26.png", kFigD));
    }

    printf("x17p26_figs: X4 FC interval [%.4e, %.4e], observed percentile in the frozen toys %.3f, empirical p %.3f\n", ulFClo, ulFCup, pctUL, pEmp);

    // ------------------------------------------------------------ x17c figures (blocks of x17c_summary.C)
    Frozen z;
    if (!ReadConfig(z)) return;
    const double FA[3] = {FBR, F2, F3};
    std::vector<x17cs::ToySet> bo2;
    const std::vector<int> f4 = {0, 500, 1000, 1500}, f1 = {0};
    for (double mm : {16.7, 16.85, 17.0}) { x17cs::ToySet T; if (x17cs::ReadToys(mm, 0., f4, T)) bo2.push_back(T); }
    for (double mm : {10.0, 25.0, 40.0}) { x17cs::ToySet T; if (x17cs::ReadToys(mm, 0., f1, T)) bo2.push_back(T); }
    std::vector<x17cs::ToySet> inj;
    for (double b : {1e-5, 3e-5, 6e-5}) { x17cs::ToySet T; if (x17cs::ReadToys(17.0, b, f1, T)) inj.push_back(T); }
    // toy quantiles of the expected limit, PDG 2026: (1) x 2.76/2.68 x k_Asimov(m), (2) x F2, (3) x F3
    struct Q { double m; int norm; double t[5]; };
    std::vector<Q> qs;
    for (const auto& T : bo2) for (int a = 0; a < 3; ++a) {
        const double f = a == 0 ? FBR * kAsi[(int)std::lround(100 * T.m)] : (a == 1 ? F2 : F3);
        Q q; q.m = T.m; q.norm = a + 1;
        for (int k = 0; k < 5; ++k) q.t[k] = f * x17cs::Quant(T.ul[a], x17cs::kQ[k]);
        qs.push_back(q);
        if (a == 0) printf("x17p26_figs: expected toys (1) m %.2f: median %.4e [%.4e, %.4e]\n", T.m, q.t[2], q.t[1], q.t[3]);
    }
    std::vector<double> ix, iy[3], iyE[3];
    for (const auto& T : inj) for (int a = 0; a < 3; ++a) {
        const auto& r = T.rhat[a];
        double mn = 0; for (double v : r) mn += v; mn /= std::max<size_t>(1, r.size());
        double rms = 0; for (double v : r) rms += (v - mn) * (v - mn); rms = std::sqrt(rms / std::max<size_t>(1, r.size() - 1));
        if (a == 0) ix.push_back(T.inj);
        iy[a].push_back(mn); iyE[a].push_back(rms / std::sqrt((double)std::max<size_t>(1, r.size())));
    }
    // ---------------------------------------------------------------- figure: expected limit vs m_X
    {
        TCanvas c("ce", "", 950, 750);
        gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.12); gPad->SetLogy(); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.04);
        const int n = (int)am.size();
        TGraphAsymmErrors b2(n), b1(n);
        TGraph med(n), med2(n), med3(n), proj;
        for (int k = 0; k < n; ++k) {
            b2.SetPoint(k, am[k], ae[0][2][k]); b2.SetPointError(k, 0, 0, ae[0][2][k] - ae[0][0][k], ae[0][4][k] - ae[0][2][k]);
            b1.SetPoint(k, am[k], ae[0][2][k]); b1.SetPointError(k, 0, 0, ae[0][2][k] - ae[0][1][k], ae[0][3][k] - ae[0][2][k]);
            med.SetPoint(k, am[k], ae[0][2][k]); med2.SetPoint(k, am[k], ae[1][2][k]); med3.SetPoint(k, am[k], ae[2][2][k]);
        }
        FILE* fp = fopen(Form("%s/x17ab_proj_mass.txt", kTab), "r");
        if (fp) {
            char line[512];
            while (fgets(line, sizeof line, fp)) { double mm, rx, a0, a1, a2; if (line[0] != '#' && sscanf(line, "%lf %lf %lf %lf %lf", &mm, &rx, &a0, &a1, &a2) == 5) proj.SetPoint(proj.GetN(), mm, FBR * a2); }
            fclose(fp);
        }
        b2.SetFillColor(kOrange); b1.SetFillColor(kGreen + 1); b2.SetLineColor(kOrange); b1.SetLineColor(kGreen + 1);
        b2.SetTitle(";m_{X} [MeV/c^{2}];90% CL_{s} UL on BR(#eta#rightarrow#pi^{+}#pi^{-}X) BR(X#rightarrowe^{+}e^{-})");
        b2.GetXaxis()->SetTitleSize(0.045); b2.GetYaxis()->SetTitleSize(0.042); b2.GetXaxis()->SetLabelSize(0.04); b2.GetYaxis()->SetLabelSize(0.04);
        b2.GetYaxis()->SetTitleOffset(1.6);
        b2.SetMinimum(3e-6); b2.SetMaximum(2e-3); b2.GetXaxis()->SetLimits(9.5, 40.5);
        b2.Draw("A3"); b1.Draw("3 SAME");
        med.SetLineStyle(2); med.SetLineWidth(3); med.Draw("L SAME");
        med2.SetLineColor(kBlue + 1); med2.SetLineWidth(2); med2.SetLineStyle(7); med2.Draw("L SAME");
        med3.SetLineColor(kMagenta + 2); med3.SetLineWidth(2); med3.SetLineStyle(3); med3.Draw("L SAME");
        if (proj.GetN() > 0) { proj.SetLineColor(kGray + 2); proj.SetLineWidth(2); proj.SetLineStyle(4); proj.Draw("L SAME"); }
        // toy quantiles (normalisation 1)
        TGraphAsymmErrors tq;
        for (const auto& q : qs) if (q.norm == 1) { const int k = tq.GetN(); tq.SetPoint(k, q.m, q.t[2]); tq.SetPointError(k, 0, 0, q.t[2] - q.t[1], q.t[3] - q.t[2]); }
        tq.SetMarkerStyle(20); tq.SetMarkerColor(kBlack); tq.SetLineColor(kBlack); tq.SetLineWidth(2); tq.Draw("P SAME");
        TLegend lg(0.30, 0.68, 0.96, 0.95); lg.SetBorderSize(0); lg.SetTextSize(0.027); lg.SetFillStyle(0); lg.SetNColumns(1);
        lg.AddEntry(&med, "(1) N5 fitted: Asimov median (PRIMARY)", "l");
        lg.AddEntry(&b1, "(1) #pm 1#sigma", "f"); lg.AddEntry(&b2, "(1) #pm 2#sigma", "f");
        lg.AddEntry(&tq, "(1) toys: median, 16-84 % (2000 / 500 per mass)", "pe");
        lg.AddEntry(&med2, Form("(2) #eta yield from a_{#eta} = %.2f #pm %.2f: median", z.a0 / F2, z.aSig / F2), "l");
        lg.AddEntry(&med3, Form("(3) #eta yield from the peak, a_{P} = %.2f #pm %.2f: median", z.aP / F3, z.aPSig / F3), "l");
        if (proj.GetN() > 0) lg.AddEntry(&proj, "pre-data projection (x17ab, simulation Asimov)", "l");
        lg.Draw();
        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.028);
        tx.DrawLatex(0.17, 0.16, "EXPECTED (no data in M(e^{+}e^{-}) 12-28 MeV/c^{2} used); S0, isBest #geq 0, X2D");
        tx.DrawLatex(0.17, 0.205, "background from the data side bands; L = 6079.2 nb^{-1}; 2D likelihood, q~_{r}; PDG 2026");
        c.SaveAs(Form("%s/x17c_expected_pdg26.pdf", kFig)); c.SaveAs(Form("%s/x17c_expected_pdg26.png", kFig));
    }
    // ---------------------------------------------------------------- figure: toy UL distribution at 17 MeV/c2
    for (const auto& T : bo2) {
        if (std::fabs(T.m - 17.0) > 1e-6) continue;
        TCanvas c("cu", "", 950, 700);
        gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13); gPad->SetLogx();
        const int nb = 60;
        std::vector<double> be(nb + 1);
        for (int k = 0; k <= nb; ++k) be[k] = 5e-6 * std::pow(40., (double)k / nb);
        TH1D h1("hu1", ";90% CL_{s} UL on BR_{X} BR_{ee} (background-only toys, m_{X} = 17 MeV/c^{2});toys per bin", nb, be.data());
        TH1D h2("hu2", "", nb, be.data()), h3("hu3", "", nb, be.data());
        for (double v : T.ul[0]) h1.Fill(F1a17 * v);
        for (double v : T.ul[1]) h2.Fill(F2 * v);
        for (double v : T.ul[2]) h3.Fill(F3 * v);
        h1.SetLineColor(kBlack); h1.SetLineWidth(3); h2.SetLineColor(kBlue + 1); h2.SetLineWidth(2); h2.SetLineStyle(7);
        h3.SetLineColor(kMagenta + 2); h3.SetLineWidth(2); h3.SetLineStyle(3);
        h1.SetMaximum(1.35 * std::max(h1.GetMaximum(), std::max(h2.GetMaximum(), h3.GetMaximum())));
        h1.GetXaxis()->SetTitleSize(0.042); h1.GetYaxis()->SetTitleSize(0.045);
        h1.Draw("HIST"); h2.Draw("HIST SAME"); h3.Draw("HIST SAME");
        std::vector<TLine*> ls;
        for (int k = 0; k < 5; ++k) {
            const double a = asimovAt(17.0, 0, k);
            TLine* l = new TLine(a, 0, a, h1.GetMaximum() * 0.8);
            l->SetLineColor(k == 2 ? kRed + 1 : kOrange + 7); l->SetLineStyle(k == 2 ? 1 : 2); l->SetLineWidth(2); l->Draw(); ls.push_back(l);
        }
        gPad->SetTopMargin(0.05); gPad->SetRightMargin(0.04);
        TLegend lg(0.56, 0.66, 0.95, 0.93); lg.SetBorderSize(0); lg.SetTextSize(0.03); lg.SetFillStyle(0);
        lg.AddEntry(&h1, Form("(1) N5 fitted, %zu toys", T.ul[0].size()), "l");
        lg.AddEntry(&h2, "(2) a_{#eta} normalisation", "l");
        lg.AddEntry(&h3, "(3) peak normalisation", "l");
        lg.AddEntry(ls[2], "(1) Asimov median", "l"); lg.AddEntry(ls[0], "(1) Asimov #pm1#sigma, #pm2#sigma", "l");
        lg.Draw();
        c.SaveAs(Form("%s/x17c_toys_ul_pdg26.pdf", kFig)); c.SaveAs(Form("%s/x17c_toys_ul_pdg26.png", kFig));
        for (auto* l : ls) delete l;
    }
    // ---------------------------------------------------------------- figure: injection
    if (!ix.empty()) {
        TCanvas c("ci", "", 900, 700);
        gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
        TGraphErrors g[3];
        double b0[3] = {0, 0, 0};
        for (const auto& T : bo2) if (std::fabs(T.m - 17.0) < 1e-6) for (int a = 0; a < 3; ++a) { double s = 0; for (double v : T.rhat[a]) s += v; b0[a] = s / std::max<size_t>(1, T.rhat[a].size()); }
        const int col[3] = {kBlack, kBlue + 1, kMagenta + 2};
        for (int a = 0; a < 3; ++a) {
            g[a].SetPoint(0, 0., FA[a] * b0[a] * 1e5);
            for (size_t k = 0; k < ix.size(); ++k) { g[a].SetPoint(k + 1, FBR * ix[k] * 1e5 * (1 + 0.01 * a), FA[a] * iy[a][k] * 1e5); g[a].SetPointError(k + 1, 0, FA[a] * iyE[a][k] * 1e5); }
            g[a].SetMarkerStyle(20 + a); g[a].SetMarkerColor(col[a]); g[a].SetLineColor(col[a]);
        }
        g[0].SetTitle(";injected BR_{X} BR_{ee} [10^{-5}];mean fitted BR_{X} BR_{ee} [10^{-5}]");
        g[0].GetXaxis()->SetTitleSize(0.045); g[0].GetYaxis()->SetTitleSize(0.045);
        g[0].GetXaxis()->SetLimits(-0.3, 6.6); g[0].SetMinimum(-0.5); g[0].SetMaximum(8);
        g[0].Draw("AP"); g[1].Draw("P SAME"); g[2].Draw("P SAME");
        TLine d(0, 0, 6.5, 6.5); d.SetLineStyle(2); d.Draw();
        TLegend lg(0.16, 0.70, 0.62, 0.89); lg.SetBorderSize(0); lg.SetTextSize(0.032);
        lg.AddEntry(&g[0], "(1) N5 fitted", "pe"); lg.AddEntry(&g[1], "(2) a_{#eta} normalisation", "pe"); lg.AddEntry(&g[2], "(3) peak normalisation", "pe");
        lg.Draw();
        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.028);
        tx.DrawLatex(0.30, 0.17, "m_{X} = 17 MeV/c^{2}; truth = data side-band fit (r = 0) + X17; 500 toys per point"); tx.DrawLatex(0.42, 0.27, "BR in PDG 2026 units (BR_{dir} = 2.76#times10^{-4})");
        tx.DrawLatex(0.30, 0.22, "dashed: fitted = injected; (2), (3): slope N5/(a N5_{ref}) expected");
        c.SaveAs(Form("%s/x17c_inject_pdg26.pdf", kFig)); c.SaveAs(Form("%s/x17c_inject_pdg26.png", kFig));
    }
    fs->Close(); f26->Close(); c.f->Close();
    printf("x17p26_figs: done\n");
}
