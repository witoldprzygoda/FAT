// isbest_choice_sim.C
//
// Why the isBest choice between e+e- hypotheses of one event is close to a coin toss (std SMASH simulation).
// PAT3 (user/htimecut.cc, core/hpidpool.cc): per ntuple the hypothesis with the lowest reduced time chi2 gets
// isBest = 1; for a pair chi2 = (D/2)^2 (1/s1^2 + 1/s2^2), D = (tof1 - tof2) - (tof_mom1 - tof_mom2), with the
// nominal resolutions s = 0.25 ns (system 0) and 0.15 ns (otherwise).
//
// The macro (1) checks that formula on every stored pair, (2) measures D for genuine pairs and for background made
// of true leptons, (3) in events with exactly 3 leptons (two hypotheses sharing one leg, one genuine and one background)
// counts how often the genuine pair wins, in total and by the detector of the two non-shared legs.
//
// Run: root -l -b -q 'plots/isbest_choice_sim.C+'

#include <TFile.h>
#include <TH1D.h>
#include <TMath.h>
#include <TStopwatch.h>
#include <TTree.h>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace ics {
struct Leg { Float_t p, th, ph, id, par, proc, vx, vy, vz, tof, tofm, res, sys; };
struct Row { Leg a, b; Float_t best, chi2; };
constexpr double kMe = 0.51099895;
inline bool lepton(const Leg& l) { return l.id == 2 || l.id == 3; }
inline bool genuine(const Row& r) { return lepton(r.a) && lepton(r.b) && r.a.id != r.b.id && r.a.vx == r.b.vx && r.a.vy == r.b.vy && r.a.vz == r.b.vz && r.a.par == r.b.par && r.a.proc == r.b.proc; }
inline bool sameP(const Leg& a, const Leg& b) { return a.p == b.p && a.th == b.th && a.ph == b.ph; }
inline double Dt(const Row& r) { return (double(r.a.tof) - r.b.tof) - (double(r.a.tofm) - r.b.tofm); }
inline double mass(const Row& r) { const double d = TMath::DegToRad(); const double ax = std::sin(r.a.th * d) * std::cos(r.a.ph * d), ay = std::sin(r.a.th * d) * std::sin(r.a.ph * d), az = std::cos(r.a.th * d), bx = std::sin(r.b.th * d) * std::cos(r.b.ph * d), by = std::sin(r.b.th * d) * std::sin(r.b.ph * d), bz = std::cos(r.b.th * d), Ea = std::sqrt(double(r.a.p) * r.a.p + kMe * kMe), Eb = std::sqrt(double(r.b.p) * r.b.p + kMe * kMe); const double m2 = 2 * kMe * kMe + 2 * (Ea * Eb - double(r.a.p) * r.b.p * (ax * bx + ay * by + az * bz)); return m2 > 0 ? std::sqrt(m2) / 1000. : 0.; }
}

void isbest_choice_sim(const char* file = "/hdd3/przygoda/HADES/files/sim/smash_std_sim_epem.root", const char* out = "plots/output/isbest_choice_sim.root", double vertex_z_min = -500.)
{
    using namespace ics;
    TStopwatch sw; sw.Start();
    TFile f(file); TTree* t = (TTree*)f.Get("EpEm_ID"); t->SetBranchStatus("*", 0);
    auto sb = [&](const std::string& n, Float_t* a) { if (!t->GetBranch(n.c_str())) { printf("missing %s\n", n.c_str()); return; } t->SetBranchStatus(n.c_str(), 1); t->SetBranchAddress(n.c_str(), a); };
    Float_t ev, vx, vy, vz; Row r;
    sb("event", &ev); sb("eVertReco_x", &vx); sb("eVertReco_y", &vy); sb("eVertReco_z", &vz); sb("isBest", &r.best); sb("ep_hyp_best_chi2_red", &r.chi2);
    auto leg = [&](const std::string& p, Leg& l) { sb(p + "p", &l.p); sb(p + "theta", &l.th); sb(p + "phi", &l.ph); sb(p + "sim_id", &l.id); sb(p + "sim_parentid", &l.par); sb(p + "sim_processid", &l.proc); sb(p + "sim_vertexx", &l.vx); sb(p + "sim_vertexy", &l.vy); sb(p + "sim_vertexz", &l.vz); sb(p + "tof", &l.tof); sb(p + "tof_mom", &l.tofm); sb(p + "resolution", &l.res); sb(p + "system", &l.sys); };
    leg("ep_", r.a); leg("em_", r.b);

    TFile fo(out, "RECREATE");
    TH1D hD[3][2]; const char* tn[3] = {"gen", "lep", "fake"}; const char* cn[2] = {"n2", "n3"};
    for (int i = 0; i < 3; ++i) for (int c = 0; c < 2; ++c) { hD[i][c] = TH1D(Form("hD_%s_%s", tn[i], cn[c]), ";D [ns];pairs", 240, -6, 6); }
    long nChk = 0, nBad = 0;
    // counts: [mass window 0 = M<0.14, 1 = 0.14-0.80, 2 = all] of (win genuine, total) with detector split of the non-shared legs
    double win[3][4] = {{0}}, tot[3][4] = {{0}};   // 0 = all, 1 = same detector, 2 = genuine leg in RPC (system 0) vs rival not, 3 = genuine leg not RPC vs rival RPC
    double rpcw = 0, rpct = 0, sameDw = 0, sameDt = 0;
    std::vector<Row> g; bool have = false; Float_t c1 = 0, c2 = 0, c3 = 0, c4 = 0;
    auto flush = [&]() {
        if (g.empty()) return;
        std::vector<const Leg*> P, M; for (auto& x : g) { bool a = false, b = false; for (auto* l : P) if (sameP(*l, x.a)) a = true; for (auto* l : M) if (sameP(*l, x.b)) b = true; if (!a) P.push_back(&x.a); if (!b) M.push_back(&x.b); }
        const int n = (int)P.size() + (int)M.size();
        for (auto& x : g) { const double D = Dt(x), s = 0.25 * D * D * (1. / (x.a.res * x.a.res) + 1. / (x.b.res * x.b.res)); ++nChk; if (std::fabs(s - x.chi2) > 1e-3 * std::max(1.f, x.chi2)) ++nBad;
            if (n <= 3) { const int c = n == 2 ? 0 : 1; const int k = genuine(x) ? 0 : ((lepton(x.a) && lepton(x.b)) ? 1 : 2); hD[k][c].Fill(D); } }
        if (n == 3 && g.size() == 2) {
            const bool shareE = sameP(g[0].b, g[1].b);
            const Leg& x0 = shareE ? g[0].a : g[0].b; const Leg& x1 = shareE ? g[1].a : g[1].b;
            if (((x0.sys == 0) != (x1.sys == 0))) { ++rpct; if ((x0.sys == 0 && g[0].best > 0.5f) || (x1.sys == 0 && g[1].best > 0.5f)) ++rpcw; } else { ++sameDt; if (g[0].best > 0.5f) ++sameDw; }
            const bool g0 = genuine(g[0]), g1 = genuine(g[1]);
            if (g0 != g1) { const Row& G = g0 ? g[0] : g[1]; const Leg& xG = g0 ? x0 : x1; const Leg& xB = g0 ? x1 : x0; const double m = mass(G); const int w = m < 0.14 ? 0 : (m < 0.80 ? 1 : -1); const int d = ((xG.sys == 0) == (xB.sys == 0)) ? 1 : (xG.sys == 0 ? 2 : 3); const bool won = G.best > 0.5f;
                for (int wi : {w, 2}) { if (wi < 0) continue; tot[wi][0] += 1; tot[wi][d] += 1; if (won) { win[wi][0] += 1; win[wi][d] += 1; } } }
        }
        g.clear();
    };
    const Long64_t N = t->GetEntries();
    for (Long64_t i = 0; i < N; ++i) { t->GetEntry(i);
        if (!have || ev != c1 || vx != c2 || vy != c3 || vz != c4) { flush(); c1 = ev; c2 = vx; c3 = vy; c4 = vz; have = true; }
        if (vz <= vertex_z_min) continue; g.push_back(r); }
    flush();
    printf("formula check: %ld pairs, |chi2(formula) - chi2(stored)| > 0.1 %%: %ld\n", nChk, nBad);
    const char* wn[3] = {"M(genuine) < 0.14", "0.14-0.80", "all masses"};
    const char* dn[4] = {"all", "same detector", "genuine leg RPC, rival not", "genuine leg not RPC, rival RPC"};
    for (int w = 0; w < 3; ++w) { printf("n=3 events with one genuine and one background hypothesis, %s:\n", wn[w]); for (int d = 0; d < 4; ++d) { const double T = tot[w][d], p = T > 0 ? win[w][d] / T : 0; printf("   %-32s P(genuine wins) = %.4f +- %.4f  (%.0f events)\n", dn[d], p, T > 0 ? std::sqrt(p * (1 - p) / T) : 0., T); } }
    printf("n=3, non-shared legs in different systems: RPC-leg hypothesis wins %.0f / %.0f = %.4f; same system: first row wins %.0f / %.0f = %.4f\n", rpcw, rpct, rpcw / std::max(1., rpct), sameDw, sameDt, sameDw / std::max(1., sameDt));
    for (int i = 0; i < 3; ++i) for (int c = 0; c < 2; ++c) printf("D [ns] %-5s %s: pairs %.0f mean %.4f rms %.4f\n", tn[i], cn[c], hD[i][c].GetEntries(), hD[i][c].GetMean(), hD[i][c].GetRMS());
    // counts for the plots: bin 1 + 4 w + d, w = mass window (0: M<0.14, 1: 0.14-0.80, 2: all), d = detector split as above
    TH1D hW("h_win", "genuine pair wins", 12, 0, 12), hT("h_tot", "events", 12, 0, 12); for (int w = 0; w < 3; ++w) for (int d = 0; d < 4; ++d) { hW.SetBinContent(1 + 4 * w + d, win[w][d]); hT.SetBinContent(1 + 4 * w + d, tot[w][d]); }
    TH1D hR("h_rpc", "RPC-leg wins, RPC-leg total, first row wins, same-system total", 4, 0, 4); hR.SetBinContent(1, rpcw); hR.SetBinContent(2, rpct); hR.SetBinContent(3, sameDw); hR.SetBinContent(4, sameDt);
    fo.cd(); hW.Write(); hT.Write(); hR.Write(); for (int i = 0; i < 3; ++i) for (int c = 0; c < 2; ++c) hD[i][c].Write();
    fo.Close(); printf("written %s (%.0f s)\n", out, sw.RealTime());
}
