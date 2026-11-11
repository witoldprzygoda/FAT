// mass_ee_signal_vs_sim_channels.C
//
// Final e+e- spectrum: the data signal against the simulation broken into channels.
//
// Data (PT3, vertex cut done by FAT): S = N+- - k(M) 2 sqrt(N++ N--), all combinations (isBest >= 0), k from event mixing in
// the data (plots/output/k_factor_data.txt, written by docs/cb_study_sim_pl/make_assets.py; k = 1 below 0.05 GeV/c2).
// Simulation (PT3, eVertReco_z > -500): truth-matched genuine pairs only, never the simulated combinatorial background:
// both legs true e+ and e-, same creation vertex, parent and process. isBest >= 0 is the default selection; isBest == 1 is
// kept for comparison. Decay channel: parent (pi0 = 7, eta = 17, photon = 1 with process 6, then the grandparent of the
// photon); production: sim_geninfo of the leg (707 = primary pi0, 1717 = primary eta). geninfo is a Float_t in PAT, so
// only the small codes are compared exactly.
//
// The file holds no number of generated events, so the simulation is normalised to the data signal in [norm_lo, norm_hi]
// (pi0 Dalitz region): this is a comparison of shapes.
//
// Output: <out>.{root,pdf,png} and a printed table.
// Run: root -l -b -q 'plots/mass_ee_signal_vs_sim_channels.C+'

#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TPad.h>
#include <TStyle.h>
#include <TTree.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <string>
#include <vector>

namespace msc {
constexpr int kNb = 58;
const double* edges() {
    static double e[kNb + 1]; static bool b = false;
    if (!b) { int k = 0; for (int i = 0; i <= 20; ++i) e[k++] = i * 0.0025;
        e[k++] = 0.08; e[k++] = 0.11; e[k++] = 0.14;
        for (int i = 1; i <= 33; ++i) e[k++] = 0.14 + i * 0.02;
        e[k++] = 0.95; e[k++] = 1.40; b = true; }
    return e;
}
enum Ch { kPi0Prim, kPi0Dec, kEtaPrim, kEtaDec, kConvPi0, kConvEta, kConvOth, kOther, kNCh };
const char* kChName[kNCh] = {"pi0Dalitz_primary", "pi0Dalitz_fromDecay", "etaDalitz_primary", "etaDalitz_fromDecay", "conv_gammaFromPi0", "conv_gammaFromEta", "conv_gammaOther", "other"};
struct Leg { Float_t p, th, ph, id, par, gp, proc, vx, vy, vz, gi, gi1; };
constexpr double kMe = 0.51099895;
inline double mass(const Leg& a, const Leg& b) {
    const double d = TMath::DegToRad();
    const double ax = std::sin(a.th * d) * std::cos(a.ph * d), ay = std::sin(a.th * d) * std::sin(a.ph * d), az = std::cos(a.th * d);
    const double bx = std::sin(b.th * d) * std::cos(b.ph * d), by = std::sin(b.th * d) * std::sin(b.ph * d), bz = std::cos(b.th * d);
    const double Ea = std::sqrt(double(a.p) * a.p + kMe * kMe), Eb = std::sqrt(double(b.p) * b.p + kMe * kMe);
    const double m2 = 2 * kMe * kMe + 2 * (Ea * Eb - double(a.p) * b.p * (ax * bx + ay * by + az * bz));
    return m2 > 0 ? std::sqrt(m2) / 1000. : 0.;
}
// genuine pair: opposite true charges, same creation point, parent and process; -1 otherwise
inline int channel(const Leg& a, const Leg& b) {
    if (!((a.id == 2 && b.id == 3) || (a.id == 3 && b.id == 2))) return -1;
    if (a.vx != b.vx || a.vy != b.vy || a.vz != b.vz || a.par != b.par || a.proc != b.proc) return -1;
    if (a.par == 7) return a.gi == 707.f ? kPi0Prim : kPi0Dec;
    if (a.par == 17) return a.gi == 1717.f ? kEtaPrim : kEtaDec;
    if (a.par == 1 && a.proc == 6) return a.gp == 7 ? kConvPi0 : (a.gp == 17 ? kConvEta : kConvOth);
    return kOther;
}
void styleAxis(TH1* f, double ts, double ls, double yoff) {
    f->GetXaxis()->SetTitleSize(ts); f->GetYaxis()->SetTitleSize(ts); f->GetXaxis()->SetLabelSize(ls); f->GetYaxis()->SetLabelSize(ls);
    f->GetYaxis()->SetTitleOffset(yoff); f->GetXaxis()->SetTitleOffset(1.05);
}
}

void mass_ee_signal_vs_sim_channels(const char* sim_file = "/hdd3/przygoda/HADES/files/sim/smash_std_sim_epem.root",
                                    const char* data_file = "plots/output/isbest_fraction_betanew.root",
                                    const char* k_file = "plots/output/k_factor_data.txt",
                                    const char* out = "plots/output/mass_ee_signal_vs_sim_channels",
                                    double norm_lo = 0.05, double norm_hi = 0.11, double vz_min = -500.)
{
    using namespace msc;
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    // ---- simulation: truth-matched pairs per channel
    TFile fs(sim_file); TTree* t = (TTree*)fs.Get("EpEm_ID"); if (!t) { printf("no EpEm_ID in %s\n", sim_file); return; }
    t->SetBranchStatus("*", 0);
    auto sb = [&](const std::string& n, Float_t* a) { if (!t->GetBranch(n.c_str())) { printf("missing branch %s\n", n.c_str()); return; } t->SetBranchStatus(n.c_str(), 1); t->SetBranchAddress(n.c_str(), a); };
    Float_t vz, trig, best; Leg A, B;
    sb("eVertReco_z", &vz); sb("trigbit", &trig); sb("isBest", &best);
    auto leg = [&](const std::string& p, Leg& l) { sb(p + "p", &l.p); sb(p + "theta", &l.th); sb(p + "phi", &l.ph); sb(p + "sim_id", &l.id); sb(p + "sim_parentid", &l.par); sb(p + "sim_grandparentid", &l.gp); sb(p + "sim_processid", &l.proc); sb(p + "sim_vertexx", &l.vx); sb(p + "sim_vertexy", &l.vy); sb(p + "sim_vertexz", &l.vz); sb(p + "sim_geninfo", &l.gi); sb(p + "sim_geninfo1", &l.gi1); };
    leg("ep_", A); leg("em_", B);
    TFile fo(Form("%s.root", out), "RECREATE");
    TH1D* hS[kNCh][2]; const char* vn[2] = {"all", "best"};
    for (int c = 0; c < kNCh; ++c) for (int v = 0; v < 2; ++v) { hS[c][v] = new TH1D(Form("sim_%s_%s", kChName[c], vn[v]), ";M_{ee} [GeV/c^{2}];pairs", kNb, edges()); hS[c][v]->Sumw2(); }
    std::map<int, long> code1[kNCh];
    const Long64_t N = t->GetEntries();
    for (Long64_t i = 0; i < N; ++i) {
        t->GetEntry(i);
        if (vz <= vz_min || best < -0.5f || !(int(trig) & 8192)) continue;
        const int c = channel(A, B); if (c < 0) continue;
        const double m = mass(A, B); hS[c][0]->Fill(m); if (best > 0.5f) hS[c][1]->Fill(m);
        ++code1[c][int(A.gi1)];
    }
    // ---- data
    TFile fd(data_file); TH1* dN = (TH1*)fd.Get("h_all_pm_PT3"); TH1* dP = (TH1*)fd.Get("h_all_pp_PT3"); TH1* dM = (TH1*)fd.Get("h_all_mm_PT3");
    if (!dN || !dP || !dM) { printf("missing data histograms in %s\n", data_file); return; }
    std::vector<double> kv(kNb, 1.), kb(kNb, 0.); { std::ifstream in(k_file); if (!in) { printf("no %s\n", k_file); return; } double lo, hi, k, b; int i = 0; while (in >> lo >> hi >> k >> b && i < kNb) { kv[i] = k; kb[i] = b; ++i; } if (i != kNb) { printf("%s: %d rows instead of %d\n", k_file, i, kNb); return; } }
    fo.cd();
    TH1D* hD = new TH1D("data_signal", ";M_{ee} [GeV/c^{2}];pairs", kNb, edges()), *hD0 = new TH1D("data_signal_k1", ";M_{ee} [GeV/c^{2}];pairs", kNb, edges());
    for (int b = 1; b <= kNb; ++b) {
        const double n = dN->GetBinContent(b), p = dP->GetBinContent(b), m = dM->GetBinContent(b), cb = (p > 0 && m > 0) ? 2 * std::sqrt(p * m) : 0., ecb = (p > 0 && m > 0) ? cb * 0.5 * std::sqrt(1 / p + 1 / m) : 0.;
        hD->SetBinContent(b, n - kv[b - 1] * cb); hD->SetBinError(b, std::sqrt(n + std::pow(kv[b - 1] * ecb, 2) + std::pow(kb[b - 1] * cb, 2)));
        hD0->SetBinContent(b, n - cb); hD0->SetBinError(b, std::sqrt(n + ecb * ecb));
    }
    // ---- totals and normalisation
    TH1D* hT[2]; for (int v = 0; v < 2; ++v) { hT[v] = (TH1D*)hS[0][v]->Clone(Form("sim_total_%s", vn[v])); for (int c = 1; c < kNCh; ++c) hT[v]->Add(hS[c][v]); }
    const int bn1 = hD->FindBin(norm_lo + 1e-6), bn2 = hD->FindBin(norm_hi - 1e-6);
    const double fnorm = hD->Integral(bn1, bn2) / hT[0]->Integral(bn1, bn2);
    TH1D* hPi = (TH1D*)hS[kPi0Prim][0]->Clone("sim_pi0Dalitz"); hPi->Add(hS[kPi0Dec][0]);
    TH1D* hEta = (TH1D*)hS[kEtaPrim][0]->Clone("sim_etaDalitz"); hEta->Add(hS[kEtaDec][0]);
    TH1D* hCon = (TH1D*)hS[kConvPi0][0]->Clone("sim_conversions"); hCon->Add(hS[kConvEta][0]); hCon->Add(hS[kConvOth][0]);
    auto perMass = [&](TH1D* h, double scale) { TH1D* o = (TH1D*)h->Clone(Form("%s_dNdM", h->GetName())); for (int b = 1; b <= kNb; ++b) { const double w = h->GetBinWidth(b); o->SetBinContent(b, h->GetBinContent(b) * scale / w); o->SetBinError(b, h->GetBinError(b) * scale / w); } return o; };
    TH1D* pD = perMass(hD, 1.), *pT = perMass(hT[0], fnorm), *pB = perMass(hT[1], fnorm), *pPi = perMass(hPi, fnorm), *pEta = perMass(hEta, fnorm), *pCon = perMass(hCon, fnorm);
    // ---- table
    auto integ = [](TH1* h, double a, double b) { return h->Integral(h->FindBin(a + 1e-6), h->FindBin(b - 1e-6)); };
    printf("normalisation (data signal / sim genuine, %.2f-%.2f): %.3f\n", norm_lo, norm_hi, fnorm);
    for (int c = 0; c < kNCh; ++c) { printf("  %-22s isBest>=0 %9.0f  isBest==1 kept %.4f  geninfo1:", kChName[c], hS[c][0]->Integral(), hS[c][0]->Integral() > 0 ? hS[c][1]->Integral() / hS[c][0]->Integral() : 0.); std::vector<std::pair<long, int>> v; for (auto& x : code1[c]) v.push_back({x.second, x.first}); std::sort(v.rbegin(), v.rend()); for (size_t i = 0; i < std::min<size_t>(v.size(), 3); ++i) printf(" %d (%ld)", v[i].second, v[i].first); printf("\n"); }
    printf(" window     | data S(k) | sim total x f | data/sim | isBest==1 / isBest>=0 (sim)\n");
    const double W[6][2] = {{0.05, 0.11}, {0.11, 0.14}, {0.14, 0.20}, {0.20, 0.30}, {0.30, 0.42}, {0.42, 0.56}};
    for (auto& w : W) { const double d = integ(hD, w[0], w[1]), s = integ(hT[0], w[0], w[1]) * fnorm; printf(" %.2f-%.2f | %9.0f | %9.0f | %.3f | %.4f\n", w[0], w[1], d, s, d / s, integ(hT[1], w[0], w[1]) / integ(hT[0], w[0], w[1])); }
    // ---- figure
    TCanvas c("c", "", 1000, 950); TPad p1("p1", "", 0, 0.30, 1, 1), p2("p2", "", 0, 0, 1, 0.30);
    p1.SetBottomMargin(0.02); p1.SetLeftMargin(0.14); p1.SetRightMargin(0.03); p1.SetTopMargin(0.03); p1.SetLogy(); p1.SetTicks(1, 1);
    p2.SetTopMargin(0.02); p2.SetBottomMargin(0.30); p2.SetLeftMargin(0.14); p2.SetRightMargin(0.03); p2.SetLogy(); p2.SetTicks(1, 1); p1.Draw(); p2.Draw();
    const double xmax = 0.80; p1.cd(); TH1D* fr = (TH1D*)pD->Clone("frame"); fr->Reset(); fr->GetXaxis()->SetRangeUser(0, xmax); fr->SetMinimum(2e3); fr->SetMaximum(3e10); fr->GetYaxis()->SetTitle("dN/dM_{ee} [(GeV/c^{2})^{-1}]"); fr->GetXaxis()->SetLabelSize(0); styleAxis(fr, 0.052, 0.046, 1.3); fr->GetXaxis()->SetLabelSize(0); fr->Draw("AXIS");
    auto line = [](TH1* h, int col, int st, int w) { h->SetLineColor(col); h->SetLineStyle(st); h->SetLineWidth(w); h->SetMarkerSize(0); h->Draw("HIST SAME"); };
    line(pPi, kOrange + 7, 1, 2); line(pEta, kViolet + 1, 1, 2); line(pCon, kGreen + 2, 1, 2); line(pT, kBlack, 1, 3); line(pB, kGray + 2, 2, 2);
    pD->SetMarkerStyle(20); pD->SetMarkerSize(1.0); pD->SetMarkerColor(kBlue + 1); pD->SetLineColor(kBlue + 1); pD->Draw("E1 X0 SAME");
    TLegend lg(0.50, 0.60, 0.95, 0.95); lg.SetBorderSize(0); lg.SetFillStyle(0); lg.SetTextSize(0.040);
    lg.AddEntry(pD, "data: signal = N_{+-} - k 2#sqrt{N_{++}N_{--}}", "pe"); lg.AddEntry(pT, "sim: genuine pairs, isBest #geq 0", "l"); lg.AddEntry(pB, "sim: genuine pairs, isBest == 1", "l"); lg.AddEntry(pPi, "sim: #pi^{0} Dalitz", "l"); lg.AddEntry(pEta, "sim: #eta Dalitz", "l"); lg.AddEntry(pCon, "sim: #gamma conversions", "l"); lg.Draw();
    p2.cd(); TH1D* r = (TH1D*)pD->Clone("ratio_data_sim"); r->Divide(pT); r->GetXaxis()->SetRangeUser(0, xmax); r->SetMinimum(0.1); r->SetMaximum(50); r->GetXaxis()->SetTitle("M_{ee} [GeV/c^{2}]"); r->GetYaxis()->SetTitle("data / sim"); styleAxis(r, 0.115, 0.10, 0.58); r->GetYaxis()->SetNdivisions(503); r->Draw("E1 X0"); TLine l1(0, 1, xmax, 1); l1.SetLineStyle(2); l1.Draw();
    c.SaveAs(Form("%s.pdf", out)); c.SaveAs(Form("%s.png", out));
    fo.cd(); fo.Write(); fo.Close();
    printf("written %s.{root,pdf,png}\n", out);
}
