// cb_closure_sim.C
//
// Closure test of the like-sign combinatorial-background estimator in
// simulation, for all e+e- combinations and for the best-pair selections.
// Input: a SMASH PAT file with EpEm_ID, EpEp_ID and EmEm_ID (the only one so
// far: /hdd2/przygoda/hades/pp45/049/lepton_sim_pid.root, standard SMASH,
// GEANT Dalitz decays and conversions, no weights, trigbit = 0 -> all events).
//
// For every e+e- combination the origin is classified with the particle history:
//   genuine Dalitz pair      : both legs true leptons (sim_id 2/3) created at the
//                              same GEANT vertex by the same parent decay
//                              (equal parentid > 1, processid, bitwise vertex)
//   genuine conversion pair  : same, parent photon (parentid 1, pair production)
//   genuine primary pair     : both primary (parentid -1), same geninfo and
//                              geninfo2 (decays inside the event generator)
//   combinatorial (CB)       : both true leptons from different sources
//                              (Dalitz + conversion, conversion of gamma1 +
//                              conversion of gamma2 even from the same pi0, ...)
//   fake                     : at least one leg is not a true lepton (hadron misID)
// Repeated hypotheses (same two reconstructed legs written twice by PAT) are
// dropped within an event and counted separately.
//
// Three variants on the same events (event cuts: isBest >= 0, eVertReco_z > -500):
//   all  : every unique combination
//   prod : the stored isBest == 1 (sim production: tracking-chi2 ranking)
//   time : argmin of the reduced relative-time chi2 per event and pattern
//          (= isBest of the exp betanew production)
// Like-sign pairs (EpEp_ID, EmEm_ID) go through the same variants (rankings
// within their own pattern, as in PAT).
//
// Per variant and mass bin: N+-, N++, N--, B_truth (CB + fake), S_truth
// (genuine), CB_hat = 2 sqrt(N++ N--), and the comparison CB_hat vs B_truth,
// S_hat = N+- - CB_hat vs S_truth. Mass from plain p (FAT m_ee).
//
// Run: root -l -b -q 'plots/cb_closure_sim.C+("file.root","tag")'

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TMath.h>
#include <TStyle.h>
#include <TTree.h>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace cbcs {
constexpr int kNb = 58;
const double* edges() {
    static double e[kNb + 1]; static bool b = false;
    if (!b) { int k = 0; for (int i = 0; i <= 20; ++i) e[k++] = i * 0.0025;
        e[k++] = 0.08; e[k++] = 0.11; e[k++] = 0.14;
        for (int i = 1; i <= 33; ++i) e[k++] = 0.14 + i * 0.02;
        e[k++] = 0.95; e[k++] = 1.40; b = true; }
    return e;
}
struct Leg { Float_t p, th, ph, id, par, proc, vx, vy, vz, g, g2, ghost; };
struct Hyp { Leg a, b; float chi2, best, m; int cls; };   // cls: 0 dalitz,1 conv,2 primary,3 CB,4 fake (like-sign: 3/4, 5 same particle)
constexpr int kNV = 3;
const char* kV[kNV] = {"all", "prod", "time"};
enum { kDal, kConv, kPrim, kCB, kFake, kSamePart, kNCls };
const char* kClsName[kNCls] = {"genuine Dalitz", "genuine conversion", "genuine primary", "combinatorial (true leptons)", "fake leg", "same particle twice"};

double mass(const Leg& a, const Leg& b) {
    const double me = 0.51099895;
    const double d2r = TMath::DegToRad();
    const double ax = a.p * std::sin(a.th * d2r) * std::cos(a.ph * d2r), ay = a.p * std::sin(a.th * d2r) * std::sin(a.ph * d2r), az = a.p * std::cos(a.th * d2r);
    const double bx = b.p * std::sin(b.th * d2r) * std::cos(b.ph * d2r), by = b.p * std::sin(b.th * d2r) * std::sin(b.ph * d2r), bz = b.p * std::cos(b.th * d2r);
    const double ea = std::sqrt(a.p * a.p + me * me), eb = std::sqrt(b.p * b.p + me * me);
    const double m2 = 2 * me * me + 2 * (ea * eb - (ax * bx + ay * by + az * bz));
    return m2 > 0 ? std::sqrt(m2) / 1000. : 0.;
}
bool sameLeg(const Leg& x, const Leg& y) { return x.p == y.p && x.th == y.th && x.ph == y.ph; }
bool trueLep(const Leg& x) { return x.id == 2 || x.id == 3; }
int classify(const Leg& a, const Leg& b, bool unlike) {
    if (!trueLep(a) || !trueLep(b)) return kFake;
    const bool sameVtx = a.vx == b.vx && a.vy == b.vy && a.vz == b.vz;
    if (!unlike) return (sameVtx && a.par == b.par && a.proc == b.proc && a.g == b.g) ? kSamePart : kCB;
    if (a.par == -1 && b.par == -1) return (a.g == b.g && a.g2 == b.g2) ? kPrim : kCB;
    if (a.par == b.par && a.par > 0 && a.proc == b.proc && sameVtx) return a.par == 1 ? kConv : kDal;
    return kCB;
}
} // namespace cbcs

void cb_closure_sim(const char* file = "/hdd2/przygoda/hades/pp45/049/lepton_sim_pid.root",
                    const char* tag = "jan26",
                    const char* out_dir = "plots/output",
                    double vertex_z_min = -500.)
{
    using namespace cbcs;
    gStyle->SetOptStat(0);
    TFile* f = TFile::Open(file);
    if (!f || f->IsZombie()) { printf("cannot open %s\n", file); return; }

    // h[variant][pattern 0=+-,1=++,2=--][class]
    TH1D* h[kNV][3][kNCls];
    const char* pn[3] = {"pm", "pp", "mm"};
    for (int v = 0; v < kNV; ++v) for (int p = 0; p < 3; ++p) for (int c = 0; c < kNCls; ++c)
        h[v][p][c] = new TH1D(Form("h_%s_%s_c%d", kV[v], pn[p], c), ";M_{ee} [GeV/c^{2}];pairs", kNb, edges());
    long nDup[3] = {0, 0, 0}, nEvt[3] = {0, 0, 0}, nRows[3] = {0, 0, 0}, nNoBest[3] = {0, 0, 0};

    const char* trees[3] = {"EpEm_ID", "EpEp_ID", "EmEm_ID"};
    const char* pre1[3] = {"ep_", "ep1_", "em1_"};
    const char* pre2[3] = {"em_", "ep2_", "em2_"};
    for (int p = 0; p < 3; ++p) {
        TTree* t = (TTree*)f->Get(trees[p]);
        if (!t) { printf("missing %s\n", trees[p]); return; }
        t->SetBranchStatus("*", 0);
        Float_t run, ev, best, vz, chi2; Leg A, B;
        auto sb = [&](const std::string& n, Float_t* a) { if (!t->GetBranch(n.c_str())) { printf("missing branch %s\n", n.c_str()); throw 1; } t->SetBranchStatus(n.c_str(), 1); t->SetBranchAddress(n.c_str(), a); };
        auto setLeg = [&](const std::string& pr, Leg& L) {
            sb(pr + "p", &L.p); sb(pr + "theta", &L.th); sb(pr + "phi", &L.ph); sb(pr + "sim_id", &L.id);
            sb(pr + "sim_parentid", &L.par); sb(pr + "sim_processid", &L.proc);
            sb(pr + "sim_vertexx", &L.vx); sb(pr + "sim_vertexy", &L.vy); sb(pr + "sim_vertexz", &L.vz);
            sb(pr + "sim_geninfo", &L.g); sb(pr + "sim_geninfo2", &L.g2); sb(pr + "sim_isghost", &L.ghost); };
        try {
            sb("runnumber", &run); sb("event", &ev); sb("isBest", &best); sb("eVertReco_z", &vz);
            sb(std::string(pre1[p]) + "hyp_best_chi2_red", &chi2);
            setLeg(pre1[p], A); setLeg(pre2[p], B);
        } catch (...) { return; }
        std::vector<Hyp> g; float cr = -1, ce = -1; bool have = false;
        auto flush = [&]() {
            if (g.empty()) return;
            // drop repeated hypotheses (like-sign: unordered)
            std::vector<Hyp> u;
            for (auto& x : g) {
                bool seen = false;
                for (auto& y : u) {
                    if (sameLeg(x.a, y.a) && sameLeg(x.b, y.b)) seen = true;
                    if (p > 0 && sameLeg(x.a, y.b) && sameLeg(x.b, y.a)) seen = true;
                }
                if (seen) ++nDup[p]; else u.push_back(x);
            }
            ++nEvt[p];
            int ib = -1, it = 0;
            for (size_t i = 0; i < u.size(); ++i) { if (u[i].best > 0.5f && ib < 0) ib = i; if (u[i].chi2 < u[it].chi2) it = i; }
            if (ib < 0) ++nNoBest[p];
            for (size_t i = 0; i < u.size(); ++i) {
                h[0][p][u[i].cls]->Fill(u[i].m);
                if ((int)i == ib) h[1][p][u[i].cls]->Fill(u[i].m);
                if ((int)i == it) h[2][p][u[i].cls]->Fill(u[i].m);
            }
            g.clear();
        };
        const Long64_t N = t->GetEntries();
        for (Long64_t i = 0; i < N; ++i) {
            t->GetEntry(i);
            if (!have || run != cr || ev != ce) { flush(); cr = run; ce = ev; have = true; }
            if (best < 0.f || vz <= vertex_z_min) continue;
            ++nRows[p];
            Hyp x; x.a = A; x.b = B; x.chi2 = chi2; x.best = best; x.m = mass(A, B); x.cls = classify(A, B, p == 0);
            g.push_back(x);
        }
        flush();
        printf("%s: rows used %ld, events %ld, repeated hypotheses dropped %ld, events without isBest==1 %ld\n", trees[p], nRows[p], nEvt[p], nDup[p], nNoBest[p]);
    }

    // ---- tables ------------------------------------------------------------
    const double rng[4][2] = {{0.0, 0.14}, {0.14, 0.80}, {0.80, 1.40}, {0.0, 1.40}};
    auto I = [&](TH1D* hh, double lo, double hi) { return hh->Integral(hh->FindBin(lo + 1e-6), hh->FindBin(hi - 1e-6)); };
    printf("\nComposition of e+e- pairs (all combinations):\n");
    for (auto& r : rng) {
        printf("  M %.2f-%.2f:", r[0], r[1]);
        for (int c = 0; c < kNCls - 1; ++c) printf("  %s %.0f", kClsName[c], I(h[0][0][c], r[0], r[1]));
        printf("\n");
    }
    printf("Like-sign pairs (all combinations): ++ CB %.0f fake %.0f same-particle %.0f | -- CB %.0f fake %.0f same-particle %.0f\n",
           I(h[0][1][kCB], 0, 1.4), I(h[0][1][kFake], 0, 1.4), I(h[0][1][kSamePart], 0, 1.4), I(h[0][2][kCB], 0, 1.4), I(h[0][2][kFake], 0, 1.4), I(h[0][2][kSamePart], 0, 1.4));

    TH1D *Npm[kNV], *Npp[kNV], *Nmm[kNV], *Btr[kNV], *Str[kNV], *CBh[kNV], *Sh[kNV];
    for (int v = 0; v < kNV; ++v) {
        Npm[v] = (TH1D*)h[v][0][0]->Clone(Form("Npm_%s", kV[v])); Npm[v]->Reset();
        Npp[v] = (TH1D*)Npm[v]->Clone(Form("Npp_%s", kV[v])); Nmm[v] = (TH1D*)Npm[v]->Clone(Form("Nmm_%s", kV[v]));
        Btr[v] = (TH1D*)Npm[v]->Clone(Form("Btruth_%s", kV[v])); Str[v] = (TH1D*)Npm[v]->Clone(Form("Struth_%s", kV[v]));
        for (int c = 0; c < kNCls; ++c) { Npm[v]->Add(h[v][0][c]); Npp[v]->Add(h[v][1][c]); Nmm[v]->Add(h[v][2][c]); }
        Str[v]->Add(h[v][0][kDal]); Str[v]->Add(h[v][0][kConv]); Str[v]->Add(h[v][0][kPrim]);
        Btr[v]->Add(h[v][0][kCB]); Btr[v]->Add(h[v][0][kFake]);
        CBh[v] = (TH1D*)Npm[v]->Clone(Form("CBhat_%s", kV[v])); CBh[v]->Reset();
        for (int b = 1; b <= kNb; ++b) {
            const double a = Npp[v]->GetBinContent(b), c = Nmm[v]->GetBinContent(b);
            if (a > 0 && c > 0) { const double cb = 2 * std::sqrt(a * c); CBh[v]->SetBinContent(b, cb); CBh[v]->SetBinError(b, 0.5 * cb * std::sqrt(1 / a + 1 / c)); }
        }
        Sh[v] = (TH1D*)Npm[v]->Clone(Form("Shat_%s", kV[v])); Sh[v]->Add(CBh[v], -1.);
    }
    printf("\nLike-sign estimate vs true background, per variant (B_truth = combinatorial + fake e+e- pairs):\n");
    printf("  variant | M range    |     N+-     N++     N--  | CB_hat = 2sqrt(N++N--) | B_truth | CB_hat/B_truth      | S_truth | S_hat = N+- - CB_hat | S_hat/S_truth\n");
    for (int v = 0; v < kNV; ++v) for (auto& r : rng) {
        const double npm = I(Npm[v], r[0], r[1]), npp = I(Npp[v], r[0], r[1]), nmm = I(Nmm[v], r[0], r[1]);
        double cb = 0, cbe2 = 0;   // sum of per-bin CB_hat
        for (int b = Npm[v]->FindBin(r[0] + 1e-6); b <= Npm[v]->FindBin(r[1] - 1e-6); ++b) { cb += CBh[v]->GetBinContent(b); cbe2 += std::pow(CBh[v]->GetBinError(b), 2); }
        const double bt = I(Btr[v], r[0], r[1]), st = I(Str[v], r[0], r[1]);
        const double rat = bt > 0 ? cb / bt : 0., erat = (bt > 0 && cb > 0) ? rat * std::sqrt(cbe2 / (cb * cb) + 1. / bt) : 0.;
        printf("  %-7s | %.2f-%.2f | %8.0f %7.0f %7.0f | %10.0f +- %6.0f    | %7.0f | %.3f +- %.3f       | %7.0f | %10.0f           | %.3f\n",
               kV[v], r[0], r[1], npm, npp, nmm, cb, std::sqrt(cbe2), bt, rat, erat, st, npm - cb, st > 0 ? (npm - cb) / st : 0.);
    }
    printf("\nWhat the rankings keep (fraction of the all-combinations count): variant | genuine e+e- | combinatorial e+e- | fake e+e- | e+e+ | e-e-\n");
    for (int v = 1; v < kNV; ++v) for (auto& r : rng) {
        auto fr = [&](TH1D* x, TH1D* y) { const double d = I(y, r[0], r[1]); return d > 0 ? I(x, r[0], r[1]) / d : 0.; };
        TH1D* comb = (TH1D*)h[v][0][kCB]; TH1D* comb0 = (TH1D*)h[0][0][kCB];
        printf("  %-5s M %.2f-%.2f | %.3f | %.3f | %.3f | %.3f | %.3f\n", kV[v], r[0], r[1], fr(Str[v], Str[0]), fr(comb, comb0),
               fr(h[v][0][kFake], h[0][0][kFake]), fr(Npp[v], Npp[0]), fr(Nmm[v], Nmm[0]));
    }

    // ---- figure ------------------------------------------------------------
    TCanvas* c = new TCanvas("c", "", 1800, 1000);
    c->Divide(3, 2);
    const char* vt[kNV] = {"all combinations", "isBest = 1 (tracking #chi^{2}, sim)", "best by time #chi^{2} (as exp isBest)"};
    for (int v = 0; v < kNV; ++v) {
        c->cd(v + 1); gPad->SetLogy(); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13); gPad->SetRightMargin(0.03);
        auto st = [](TH1D* x, int col, int mk) { x->SetMarkerStyle(mk); x->SetMarkerColor(col); x->SetLineColor(col); x->SetMarkerSize(0.9); x->GetXaxis()->SetRangeUser(0, 1.0);
            x->GetXaxis()->SetTitleSize(0.05); x->GetYaxis()->SetTitleSize(0.05); x->GetXaxis()->SetLabelSize(0.045); x->GetYaxis()->SetLabelSize(0.045); };
        st(Npm[v], kBlack, 20); st(CBh[v], kRed + 1, 21); st(Btr[v], kRed + 1, 25); st(Sh[v], kBlue + 1, 22); st(Str[v], kBlue + 1, 26);
        Npm[v]->SetMinimum(0.5); Npm[v]->SetMaximum(Npm[0]->GetMaximum() * 5);
        Npm[v]->Draw("E"); CBh[v]->Draw("E same"); Btr[v]->Draw("E same"); Str[v]->Draw("E same"); Sh[v]->Draw("E same");
        TLegend* lg = new TLegend(0.45, 0.62, 0.97, 0.88); lg->SetBorderSize(0); lg->SetTextSize(0.04); lg->SetFillStyle(0);
        lg->AddEntry(Npm[v], "all: N_{+-}", "p"); lg->AddEntry(CBh[v], "CB = 2#sqrt{N_{++}N_{--}}", "p"); lg->AddEntry(Btr[v], "CB truth", "p");
        lg->AddEntry(Sh[v], "signal = all - CB", "p"); lg->AddEntry(Str[v], "signal truth", "p"); lg->Draw();
        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.045); tx.DrawLatex(0.15, 0.92, vt[v]);
        c->cd(v + 4); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13); gPad->SetRightMargin(0.03); gPad->SetGridy();
        TH1D* r = (TH1D*)CBh[v]->Clone(Form("ratio_%s", kV[v])); r->Divide(Btr[v]);
        for (int b = 1; b <= kNb; ++b) if (Btr[v]->GetBinContent(b) < 20) { r->SetBinContent(b, 0); r->SetBinError(b, 0); }
        r->SetTitle(";M_{ee} [GeV/c^{2}];CB (like-sign) / CB truth"); r->SetMinimum(0.); r->SetMaximum(3.); r->Draw("E");
        TLatex t2; t2.SetNDC(); t2.SetTextSize(0.045); t2.DrawLatex(0.15, 0.92, vt[v]);
    }
    c->SaveAs(Form("%s/cb_closure_sim_%s.pdf", out_dir, tag));
    c->SaveAs(Form("%s/cb_closure_sim_%s.png", out_dir, tag));
    TFile fo(Form("%s/cb_closure_sim_%s.root", out_dir, tag), "RECREATE");
    for (int v = 0; v < kNV; ++v) { Npm[v]->Write(); Npp[v]->Write(); Nmm[v]->Write(); Btr[v]->Write(); Str[v]->Write(); CBh[v]->Write(); Sh[v]->Write();
        for (int p = 0; p < 3; ++p) for (int cc = 0; cc < kNCls; ++cc) h[v][p][cc]->Write(); }
    fo.Close();
}
