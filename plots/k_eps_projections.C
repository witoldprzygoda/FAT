// k_eps_projections.C
//
// How do the k factor and the isBest win probability behave in observables other than the pair mass?
// Data, PT3 (vertex cut done by FAT), all combinations, two mass windows (0.14-0.30 and 0.30-0.80 GeV/c2), four observables:
// missing mass of the e+e- pair in p+p at T = 4.53 GeV, polar angle of the pair in the CMS, |cos theta| of the e+ in the
// pair rest frame with respect to the pair direction, and pair pT.
//   k(X) = N+-_mix(X) / 2 sqrt(N++_mix(X) N--_mix(X)), event mixing with a rolling pool, pools of all events and of events
//          with n >= 3 leptons (n = identified e+ and e- of the event);
//   eps_n(X) = pairs with isBest == 1 / all pairs, in events with n = 3 and n = 4 leptons. These classes are mostly
//          combinatorial pairs, so this is a direct, truth-free look at how the ranking depends on X.
//
// Output: plots/output/k_eps_projections.root
// Run: root -l -b -q 'plots/k_eps_projections.C+'

#include <TFile.h>
#include <TH1D.h>
#include <TMath.h>
#include <TStopwatch.h>
#include <TString.h>
#include <TTree.h>
#include <cmath>
#include <cstdio>
#include <vector>

namespace kep {
constexpr double kMe = 0.51099895, kMp = 938.27208, kT = 4530.;
const double kEb = kT + kMp, kPb = std::sqrt(kEb * kEb - kMp * kMp), kEin = kEb + kMp, kBeta = kPb / kEin, kGamma = 1. / std::sqrt(1. - kBeta * kBeta);
constexpr int kNW = 2, kNV = 4;
const double kWlo[kNW] = {0.14, 0.30}, kWhi[kNW] = {0.30, 0.80};
const char* kVn[kNV] = {"mm", "thcms", "cosh", "pt"};
const int kNb[kNV] = {38, 18, 20, 25};
const double kLo[kNV] = {1600., 0., 0., 0.}, kHi[kNV] = {3500., 180., 1., 1000.};
struct MLeg { double ux, uy, uz, e, p; };
struct MEv { std::vector<MLeg> ep, em; };
inline MLeg mk(Float_t p, Float_t th, Float_t ph) { const double d = TMath::DegToRad(), t = th * d, f = ph * d; return {std::sin(t) * std::cos(f), std::sin(t) * std::sin(f), std::cos(t), std::sqrt(double(p) * p + kMe * kMe), double(p)}; }
// returns the mass window (-1 outside) and the four observables
inline int vars(const MLeg& a, const MLeg& b, double* X) {
    const double ax = a.p * a.ux, ay = a.p * a.uy, az = a.p * a.uz, bx = b.p * b.ux, by = b.p * b.uy, bz = b.p * b.uz;
    const double E = a.e + b.e, px = ax + bx, py = ay + by, pz = az + bz, p2 = px * px + py * py + pz * pz, m2 = E * E - p2;
    if (m2 <= 0) return -1;
    const double m = std::sqrt(m2), mg = m / 1000.; int w = -1; for (int i = 0; i < kNW; ++i) if (mg >= kWlo[i] && mg < kWhi[i]) w = i; if (w < 0) return -1;
    const double dE = kEin - E, dz = kPb - pz, mm2 = dE * dE - px * px - py * py - dz * dz; X[0] = mm2 > 0 ? std::sqrt(mm2) : 0.;
    const double pzc = kGamma * (pz - kBeta * E), pc = std::sqrt(px * px + py * py + pzc * pzc); X[1] = pc > 0 ? std::acos(std::max(-1., std::min(1., pzc / pc))) * TMath::RadToDeg() : 0.;
    const double g = E / m, bvx = px / E, bvy = py / E, bvz = pz / E, b2 = bvx * bvx + bvy * bvy + bvz * bvz, bp = bvx * ax + bvy * ay + bvz * az, f = ((g - 1.) / b2) * bp - g * a.e;
    const double sx = ax + f * bvx, sy = ay + f * bvy, sz = az + f * bvz, sn = std::sqrt(sx * sx + sy * sy + sz * sz), pn = std::sqrt(p2); X[2] = (sn > 0 && pn > 0) ? std::fabs(sx * px + sy * py + sz * pz) / (sn * pn) : 0.;
    X[3] = std::hypot(px, py);
    return w;
}
struct Row { Float_t pp, tp, fp, pm, tm, fm, best; };
}

void k_eps_projections(const char* file = "output_epem_exp.root", const char* out = "plots/output/k_eps_projections.root", int pool_depth = 3)
{
    using namespace kep;
    TStopwatch sw; sw.Start();
    TFile f(file); TTree* t = (TTree*)f.Get("dilepton_nt"); if (!t) { printf("no dilepton_nt\n"); return; }
    t->SetBranchStatus("*", 0);
    auto sb = [&](const char* n, Float_t* a) { if (!t->GetBranch(n)) { printf("missing %s\n", n); return; } t->SetBranchStatus(n, 1); t->SetBranchAddress(n, a); };
    Float_t run, ev, trig; Row r; sb("runnumber", &run); sb("event", &ev); sb("trigbit", &trig); sb("isBest", &r.best); sb("ep_p_rec", &r.pp); sb("ep_theta", &r.tp); sb("ep_phi", &r.fp); sb("em_p_rec", &r.pm); sb("em_theta", &r.tm); sb("em_phi", &r.fm);
    TFile fo(out, "RECREATE");
    const char* pn[3] = {"pm", "pp", "mm"}; const char* po[2] = {"all", "nge3"}; const char* cn[2] = {"n3", "n4"};
    TH1D* hM[kNW][kNV][2][3]; TH1D* hE[kNW][kNV][2][2];
    for (int w = 0; w < kNW; ++w) for (int v = 0; v < kNV; ++v) {
        for (int p = 0; p < 2; ++p) for (int q = 0; q < 3; ++q) hM[w][v][p][q] = new TH1D(Form("mix_w%d_%s_%s_%s", w, kVn[v], po[p], pn[q]), "", kNb[v], kLo[v], kHi[v]);
        for (int c = 0; c < 2; ++c) { hE[w][v][c][0] = new TH1D(Form("eps_w%d_%s_%s_all", w, kVn[v], cn[c]), "", kNb[v], kLo[v], kHi[v]); hE[w][v][c][1] = new TH1D(Form("eps_w%d_%s_%s_best", w, kVn[v], cn[c]), "", kNb[v], kLo[v], kHi[v]); }
    }
    struct Pool { std::vector<MEv> buf; size_t head = 0, filled = 0; } P[2]; P[0].buf.resize(pool_depth); P[1].buf.resize(pool_depth);
    auto fillM = [&](int p, int q, const MLeg& a, const MLeg& b) { double X[kNV]; const int w = vars(a, b, X); if (w < 0) return; for (int v = 0; v < kNV; ++v) hM[w][v][p][q]->Fill(X[v]); };
    auto mix = [&](int p, const MEv& E) { Pool& Q = P[p]; for (size_t ib = 0; ib < Q.filled; ++ib) { const MEv& B = Q.buf[ib];
            for (auto& a : E.ep) for (auto& c : B.em) fillM(p, 0, a, c);
            for (auto& c : E.em) for (auto& a : B.ep) fillM(p, 0, a, c);
            for (auto& a : E.ep) for (auto& c : B.ep) fillM(p, 1, a, c);
            for (auto& a : E.em) for (auto& c : B.em) fillM(p, 2, a, c); }
        Q.buf[Q.head] = E; Q.head = (Q.head + 1) % Q.buf.size(); if (Q.filled < Q.buf.size()) ++Q.filled; };
    std::vector<Row> g; bool have = false; Float_t c1 = 0, c2 = 0, trig0 = 0;
    auto flush = [&]() {
        if (g.empty() || trig0 != 8192.f) { g.clear(); return; }
        std::vector<const Row*> Pl, Mi;
        for (auto& x : g) { bool a = false, b = false; for (auto* l : Pl) if (l->pp == x.pp && l->tp == x.tp && l->fp == x.fp) a = true; for (auto* l : Mi) if (l->pm == x.pm && l->tm == x.tm && l->fm == x.fm) b = true; if (!a) Pl.push_back(&x); if (!b) Mi.push_back(&x); }
        MEv E; for (auto* l : Pl) E.ep.push_back(mk(l->pp, l->tp, l->fp)); for (auto* l : Mi) E.em.push_back(mk(l->pm, l->tm, l->fm));
        const int n = (int)Pl.size() + (int)Mi.size();
        mix(0, E); if (n >= 3) mix(1, E);
        if (n == 3 || n == 4) { const int c = n - 3; for (auto& x : g) { double X[kNV]; const int w = vars(mk(x.pp, x.tp, x.fp), mk(x.pm, x.tm, x.fm), X); if (w < 0) continue; for (int v = 0; v < kNV; ++v) { hE[w][v][c][0]->Fill(X[v]); if (x.best > 0.5f) hE[w][v][c][1]->Fill(X[v]); } } }
        g.clear();
    };
    const Long64_t N = t->GetEntries();
    for (Long64_t i = 0; i < N; ++i) { t->GetEntry(i);
        if (!have || ev != c1 || run != c2) { flush(); c1 = ev; c2 = run; have = true; trig0 = trig; }
        g.push_back(r);
        if (i % 10000000 == 0 && i) { printf("  row %lld / %lld (%.0f s)\n", i, N, sw.RealTime()); sw.Continue(); } }
    flush();
    fo.cd(); fo.Write(); fo.Close();
    printf("written %s (%.0f s)\n", out, sw.RealTime());
}
