// cb_origin_sim.C
//
// Origin breakdown of e+e-, e+e+ and e-e- pairs in a SMASH simulation PAT file
// with EpEm_ID, EpEp_ID and EmEm_ID (the January-2026 standard file
// /hdd2/przygoda/hades/pp45/049/lepton_sim_pid.root), for every unique pair
// (isBest >= 0) and for the best pair (isBest == 1), plus the like-sign
// estimate CB = 2 sqrt(N++ N--) against the true combinatorial background.
// Event cut eVertReco_z > -500. Repeated hypotheses (same two reconstructed
// legs; like-sign also (a,b)/(b,a)) are dropped, the first copy (with isBest) kept.
//
// Leg source (MC truth): pi0 Dalitz (parent pi0), eta Dalitz (parent eta),
// gamma conversion (parent photon, pair production), Compton / delta electron
// (parent photon with Compton process, or knock-on electron from an e+-),
// other true lepton, fake (not a true electron/positron: proton, pion, muon...).
// Pair: genuine = same creation vertex, parent and process of both legs
// (one Dalitz decay, one converted photon); otherwise combinatorial, classified
// by the two leg sources; pairs with a fake leg separately.
//
// Run: root -l -b -q 'plots/cb_origin_sim.C+("file.root")'

#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>
#include <TMath.h>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace cbo {
struct Leg { Float_t p, th, ph, id, par, proc, vx, vy, vz; };
enum Src { kPi0D, kEtaD, kConv, kCompDelta, kOtherLep, kFakeLeg };
int source(const Leg& l) {
    if (l.id != 2 && l.id != 3) return kFakeLeg;
    if (l.par == 7) return kPi0D;
    if (l.par == 17) return kEtaD;
    if (l.par == 1 && l.proc == 6) return kConv;
    if ((l.par == 1 && l.proc == 7) || ((l.par == 2 || l.par == 3) && l.proc == 10)) return kCompDelta;
    return kOtherLep;
}
// pair categories
enum Cat { kGenPi0, kGenEta, kGenConv, kGenOther,
           kDalDal, kDalConv, kConvConv, kWithCompDelta, kWithOther,
           kOneFake, kTwoFake, kSamePart, kNCat };
const char* kCatName[kNCat] = {
    "genuine: pi0 Dalitz", "genuine: eta Dalitz", "genuine: one photon conversion", "genuine: other",
    "CB: Dalitz + Dalitz (two decays)", "CB: Dalitz + conversion", "CB: conversion + conversion (two photons)",
    "CB: with Compton/delta electron", "CB: with other lepton",
    "fake: one leg", "fake: both legs", "same particle twice"};
bool isGenuineCat(int c) { return c <= kGenOther; }
bool isCBCat(int c) { return c >= kDalDal && c <= kTwoFake; }
int category(const Leg& a, const Leg& b) {
    const int sa = source(a), sb = source(b);
    const int nf = (sa == kFakeLeg) + (sb == kFakeLeg);
    if (nf == 2) return kTwoFake;
    if (nf == 1) return kOneFake;
    const bool sameOrigin = a.vx == b.vx && a.vy == b.vy && a.vz == b.vz && a.par == b.par && a.proc == b.proc;
    if (sameOrigin) {
        if (a.id == b.id) return kSamePart;               // like-sign: one particle reconstructed twice
        if (sa == kPi0D) return kGenPi0;
        if (sa == kEtaD) return kGenEta;
        if (sa == kConv) return kGenConv;
        return kGenOther;
    }
    if (sa == kCompDelta || sb == kCompDelta) return kWithCompDelta;
    if (sa == kOtherLep || sb == kOtherLep) return kWithOther;
    const bool da = (sa == kPi0D || sa == kEtaD), db = (sb == kPi0D || sb == kEtaD);
    if (da && db) return kDalDal;
    if (da || db) return kDalConv;
    return kConvConv;
}
double mass(const Leg& a, const Leg& b) {
    const double me = 0.51099895, d = TMath::DegToRad();
    auto v = [&](const Leg& l, double* x) { x[0] = l.p * std::sin(l.th * d) * std::cos(l.ph * d); x[1] = l.p * std::sin(l.th * d) * std::sin(l.ph * d); x[2] = l.p * std::cos(l.th * d); };
    double x[3], y[3]; v(a, x); v(b, y);
    const double m2 = 2 * me * me + 2 * (std::sqrt(a.p * a.p + me * me) * std::sqrt(b.p * b.p + me * me) - (x[0] * y[0] + x[1] * y[1] + x[2] * y[2]));
    return m2 > 0 ? std::sqrt(m2) / 1000. : 0.;
}
double openingAngle(const Leg& a, const Leg& b) {
    const double d = TMath::DegToRad();
    const double c = std::sin(a.th * d) * std::sin(b.th * d) * std::cos((a.ph - b.ph) * d) + std::cos(a.th * d) * std::cos(b.th * d);
    return std::acos(std::max(-1., std::min(1., c))) / d;
}
bool same(const Leg& x, const Leg& y) { return x.p == y.p && x.th == y.th && x.ph == y.ph; }
constexpr int kNR = 3;
const double kR[kNR][2] = {{0.0, 0.14}, {0.14, 0.80}, {0.0, 1.40}};
} // namespace cbo

// oa_min [deg]: pair opening-angle cut applied to every pair (e+e-, e+e+, e-e-)
// AFTER the PAT ranking, i.e. isBest as stored (no re-ranking), as in the data
// analysis where OA is a column; oa_min <= 0 disables it.
void cb_origin_sim(const char* file = "/hdd2/przygoda/hades/pp45/049/lepton_sim_pid.root", double vertex_z_min = -500.,
                   double oa_min = 0.)
{
    using namespace cbo;
    TFile* f = TFile::Open(file);
    // n[pattern][variant 0 all,1 best][range][category]
    static double n[3][2][kNR][kNCat];
    // per-bin spectra on the analysis grid (58 bins) for CB = sum_bins 2 sqrt(N++ N--)
    const int kNb = 58; double ed[kNb + 1]; { int k = 0; for (int i = 0; i <= 20; ++i) ed[k++] = i * 0.0025;
        ed[k++] = 0.08; ed[k++] = 0.11; ed[k++] = 0.14; for (int i = 1; i <= 33; ++i) ed[k++] = 0.14 + i * 0.02; ed[k++] = 0.95; ed[k++] = 1.40; }
    TH1D* hN[3][2]; TH1D* hCBt[2]; TH1D* hGen[2];
    for (int pp = 0; pp < 3; ++pp) for (int v = 0; v < 2; ++v) { hN[pp][v] = new TH1D(Form("hN_%d_%d", pp, v), "", kNb, ed); }
    for (int v = 0; v < 2; ++v) { hCBt[v] = new TH1D(Form("hCBt_%d", v), "", kNb, ed); hGen[v] = new TH1D(Form("hGen_%d", v), "", kNb, ed); }
    long dup[3] = {0, 0, 0}, rows[3] = {0, 0, 0}, evts[3] = {0, 0, 0}, oaCut[3] = {0, 0, 0};
    const char* trees[3] = {"EpEm_ID", "EpEp_ID", "EmEm_ID"};
    const char* p1[3] = {"ep_", "ep1_", "em1_"};
    const char* p2[3] = {"em_", "ep2_", "em2_"};
    for (int p = 0; p < 3; ++p) {
        TTree* t = (TTree*)f->Get(trees[p]);
        t->SetBranchStatus("*", 0);
        Float_t run, ev, best, vz; Leg A, B;
        auto sb = [&](const std::string& nm, Float_t* a) { t->SetBranchStatus(nm.c_str(), 1); t->SetBranchAddress(nm.c_str(), a); };
        auto sl = [&](const std::string& pr, Leg& L) { sb(pr + "p", &L.p); sb(pr + "theta", &L.th); sb(pr + "phi", &L.ph); sb(pr + "sim_id", &L.id);
            sb(pr + "sim_parentid", &L.par); sb(pr + "sim_processid", &L.proc); sb(pr + "sim_vertexx", &L.vx); sb(pr + "sim_vertexy", &L.vy); sb(pr + "sim_vertexz", &L.vz); };
        sb("runnumber", &run); sb("event", &ev); sb("isBest", &best); sb("eVertReco_z", &vz); sl(p1[p], A); sl(p2[p], B);
        struct H { Leg a, b; float best; };
        std::vector<H> g; float cr = -1, ce = -1; bool have = false;
        auto flush = [&]() {
            if (g.empty()) return;
            std::vector<H> u;
            for (auto& x : g) { bool s = false; for (auto& y : u) if ((same(x.a, y.a) && same(x.b, y.b)) || (p > 0 && same(x.a, y.b) && same(x.b, y.a))) s = true;
                if (s) ++dup[p]; else u.push_back(x); }
            ++evts[p];
            for (auto& x : u) {
                if (oa_min > 0 && openingAngle(x.a, x.b) <= oa_min) { ++oaCut[p]; continue; }
                const double m = mass(x.a, x.b); const int c = category(x.a, x.b);
                for (int v = 0; v < 2; ++v) {
                    if (v == 1 && x.best < 0.5f) continue;
                    hN[p][v]->Fill(m);
                    if (p == 0 && isCBCat(c)) hCBt[v]->Fill(m);
                    if (p == 0 && isGenuineCat(c)) hGen[v]->Fill(m);
                }
                for (int r = 0; r < kNR; ++r) if (m >= kR[r][0] && m < kR[r][1]) { n[p][0][r][c] += 1; if (x.best > 0.5f) n[p][1][r][c] += 1; }
            }
            g.clear();
        };
        for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            t->GetEntry(i);
            if (!have || run != cr || ev != ce) { flush(); cr = run; ce = ev; have = true; }
            if (best < 0.f || vz <= vertex_z_min) continue;
            ++rows[p]; g.push_back({A, B, best});
        }
        flush();
        printf("%s: rows %ld, events %ld, repeated hypotheses dropped %ld, pairs removed by OA <= %.1f deg: %ld\n", trees[p], rows[p], evts[p], dup[p], oa_min, oaCut[p]);
    }
    // ---- signal table: true signal, all - true CB, all - CB(formula), per variant ----
    {
        const double sr[10][2] = {{0.0, 0.02}, {0.02, 0.05}, {0.05, 0.14}, {0.14, 0.30}, {0.30, 0.50}, {0.50, 0.80}, {0.80, 1.40},
                                  {0.0, 0.14}, {0.14, 0.80}, {0.0, 1.40}};
        printf("\nSIGNAL TABLE (OA cut %.1f deg). Per variant: N+- | true signal S | N+- - true CB | N+- - CB(formula) | ratios to S\n", oa_min);
        printf("  M range     | variant   |      N+- | true S   | all-trueCB (ratio)  | all-CBformula (ratio) | true CB | CB formula\n");
        for (auto& r : sr) for (int v = 0; v < 2; ++v) {
            const int b1 = hN[0][v]->FindBin(r[0] + 1e-6), b2 = hN[0][v]->FindBin(r[1] - 1e-6);
            double npm = 0, st = 0, bt = 0, cbf = 0;
            for (int b = b1; b <= b2; ++b) {
                npm += hN[0][v]->GetBinContent(b); st += hGen[v]->GetBinContent(b); bt += hCBt[v]->GetBinContent(b);
                const double a = hN[1][v]->GetBinContent(b), m2 = hN[2][v]->GetBinContent(b);
                if (a > 0 && m2 > 0) cbf += 2 * std::sqrt(a * m2);
            }
            printf("  %.2f-%.2f | %-9s | %8.0f | %8.0f | %8.0f (%.3f)    | %8.0f (%.3f)      | %7.0f | %7.0f\n", r[0], r[1], v == 0 ? "isBest>=0" : "isBest==1",
                   npm, st, npm - bt, st > 0 ? (npm - bt) / st : 0., npm - cbf, st > 0 ? (npm - cbf) / st : 0., bt, cbf);
        }
    }
    const char* pn[3] = {"e+e-", "e+e+", "e-e-"};
    for (int r = 0; r < kNR; ++r) {
        printf("\n================ M %.2f - %.2f GeV/c2 ================\n", kR[r][0], kR[r][1]);
        for (int p = 0; p < 3; ++p) {
            double tA = 0, tB = 0; for (int c = 0; c < kNCat; ++c) { tA += n[p][0][r][c]; tB += n[p][1][r][c]; }
            printf("%s: all (isBest>=0) %.0f | best (isBest==1) %.0f | kept %.3f\n", pn[p], tA, tB, tA > 0 ? tB / tA : 0.);
            for (int c = 0; c < kNCat; ++c) {
                if (n[p][0][r][c] < 1) continue;
                printf("   %-44s %9.0f (%5.1f%%) -> %9.0f (%5.1f%%)  kept %.3f\n", kCatName[c], n[p][0][r][c], 100 * n[p][0][r][c] / tA,
                       n[p][1][r][c], tB > 0 ? 100 * n[p][1][r][c] / tB : 0., n[p][1][r][c] / n[p][0][r][c]);
            }
        }
        for (int v = 0; v < 2; ++v) {
            double S = 0, CB = 0, Npm = 0, Npp = 0, Nmm = 0, cbLep = 0;
            for (int c = 0; c < kNCat; ++c) { Npm += n[0][v][r][c]; Npp += n[1][v][r][c]; Nmm += n[2][v][r][c];
                if (isGenuineCat(c)) S += n[0][v][r][c]; if (isCBCat(c)) CB += n[0][v][r][c]; if (c >= kDalDal && c <= kWithOther) cbLep += n[0][v][r][c]; }
            // CB as in the analysis: 2 sqrt(N++ N--) per mass bin, summed over the range
            double est = 0, est_e2 = 0;
            for (int b = hN[0][v]->FindBin(kR[r][0] + 1e-6); b <= hN[0][v]->FindBin(kR[r][1] - 1e-6); ++b) {
                const double a = hN[1][v]->GetBinContent(b), m2 = hN[2][v]->GetBinContent(b);
                if (a > 0 && m2 > 0) { const double cb = 2 * std::sqrt(a * m2); est += cb; est_e2 += 0.25 * cb * cb * (1 / a + 1 / m2); }
            }
            printf("%s: N+- %.0f = genuine %.0f + true CB %.0f (true leptons %.0f, with fake leg %.0f) | N++ %.0f N-- %.0f -> CB = sum over bins 2sqrt(N++N--) = %.0f, ratio to true CB %.3f +- %.3f | signal N+- - CB = %.0f vs genuine %.0f (ratio %.3f)\n",
                   v == 0 ? "isBest>=0" : "isBest==1", Npm, S, CB, cbLep, CB - cbLep, Npp, Nmm, est, CB > 0 ? est / CB : 0.,
                   (CB > 0 && est > 0) ? est / CB * std::sqrt(est_e2 / (est * est) + 1. / CB) : 0., Npm - est, S, S > 0 ? (Npm - est) / S : 0.);
            printf("            (for comparison, 2 sqrt of the range sums: %.0f, ratio %.3f)\n", 2 * std::sqrt(Npp * Nmm), CB > 0 ? 2 * std::sqrt(Npp * Nmm) / CB : 0.);
        }
    }
}
