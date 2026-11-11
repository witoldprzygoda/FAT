// isbest_choice_data.C
//
// Truth-free test of the isBest choice in events with exactly 3 identified leptons (two e+e- hypotheses sharing one leg).
// A hypothesis with M < 0.03 GeV/c2 is taken as a tagged true pair (conversion or pi0 Dalitz); the test asks how often it
// beats a rival with 0.14 < M < 0.80, in total and by the detector of the two non-shared legs (theta < 45 deg = RPC).
// A ranking that recognised the true pair would give well above 0.5.
//
// Run: root -l -b -q 'plots/isbest_choice_data.C+("data","output_epem_exp.root")'
//      root -l -b -q 'plots/isbest_choice_data.C+("sim","/hdd3/przygoda/HADES/files/sim/smash_std_sim_epem.root")'
#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>
#include <string>
#include <cmath>
#include <cstdio>
#include <vector>
namespace icd { struct R { Float_t pp, tp, fp, pm, tm, fm, best, m; }; constexpr double kMe = 0.51099895;
inline double mass(const R& r) { const double d = 3.14159265358979 / 180.; const double a[3] = {std::sin(r.tp * d) * std::cos(r.fp * d), std::sin(r.tp * d) * std::sin(r.fp * d), std::cos(r.tp * d)}, b[3] = {std::sin(r.tm * d) * std::cos(r.fm * d), std::sin(r.tm * d) * std::sin(r.fm * d), std::cos(r.tm * d)}; const double E1 = std::sqrt(double(r.pp) * r.pp + kMe * kMe), E2 = std::sqrt(double(r.pm) * r.pm + kMe * kMe), m2 = 2 * kMe * kMe + 2 * (E1 * E2 - double(r.pp) * r.pm * (a[0] * b[0] + a[1] * b[1] + a[2] * b[2])); return m2 > 0 ? std::sqrt(m2) / 1000. : 0.; } }
void isbest_choice_data(const char* mode = "data", const char* file = "output_epem_exp.root", double vz_min = -500., const char* out_pattern = "plots/output/isbest_choice_tag_%s.root")
{
    using namespace icd; const bool sim = std::string(mode) == "sim";
    TFile f(file); TTree* t = (TTree*)f.Get(sim ? "EpEm_ID" : "dilepton_nt"); t->SetBranchStatus("*", 0);
    auto sb = [&](const char* n, Float_t* a) { t->SetBranchStatus(n, 1); t->SetBranchAddress(n, a); };
    Float_t ev = 0, run = 0, vx = 0, vy = 0, vz = 0, trig = 0; R r; sb("event", &ev); sb("trigbit", &trig); sb("isBest", &r.best); sb("ep_theta", &r.tp); sb("ep_phi", &r.fp); sb("em_theta", &r.tm); sb("em_phi", &r.fm);
    if (sim) { sb("eVertReco_x", &vx); sb("eVertReco_y", &vy); sb("eVertReco_z", &vz); sb("ep_p", &r.pp); sb("em_p", &r.pm); } else { sb("runnumber", &run); sb("ep_p_rec", &r.pp); sb("em_p_rec", &r.pm); }
    double win[4] = {0, 0, 0, 0}, tot[4] = {0, 0, 0, 0}; std::vector<R> g; bool have = false; Float_t c1 = 0, c2 = 0, c3 = 0, c4 = 0; bool useEv = false;
    auto flush = [&]() { if (g.size() == 2 && useEv) {
            const bool shareM = g[0].pm == g[1].pm && g[0].tm == g[1].tm && g[0].fm == g[1].fm, shareP = g[0].pp == g[1].pp && g[0].tp == g[1].tp && g[0].fp == g[1].fp;
            if (shareM != shareP) { const double m0 = mass(g[0]), m1 = mass(g[1]); int lo = -1; if (m0 < 0.03 && m1 > 0.14 && m1 < 0.80) lo = 0; else if (m1 < 0.03 && m0 > 0.14 && m0 < 0.80) lo = 1;
                if (lo >= 0) { const int hi = 1 - lo; const double tl = shareM ? g[lo].tp : g[lo].tm, th = shareM ? g[hi].tp : g[hi].tm; const bool rl = tl < 45., rh = th < 45.; const int d = (rl == rh) ? 1 : (rl ? 2 : 3); const bool w = g[lo].best > 0.5f; tot[0] += 1; tot[d] += 1; if (w) { win[0] += 1; win[d] += 1; } } } }
        g.clear(); };
    for (Long64_t i = 0, n = t->GetEntries(); i < n; ++i) { t->GetEntry(i); const Float_t a = sim ? vx : run, b = sim ? vy : 0.f, c = sim ? vz : 0.f;
        if (!have || ev != c1 || a != c2 || b != c3 || c != c4) { flush(); c1 = ev; c2 = a; c3 = b; c4 = c; have = true; useEv = sim ? (int(trig) & 8192) != 0 : trig == 8192.f; }
        if (sim && vz <= vz_min) continue; g.push_back(r); }
    flush();
    { TFile fo(Form(out_pattern, mode), "RECREATE"); TH1D hw("h_win", "tagged pair wins: all, same detector, tagged leg RPC, tagged leg TOF", 4, 0, 4), ht("h_tot", "events", 4, 0, 4); for (int d = 0; d < 4; ++d) { hw.SetBinContent(d + 1, win[d]); ht.SetBinContent(d + 1, tot[d]); } hw.Write(); ht.Write(); }
    const char* dn[4] = {"all", "same detector", "tagged pair leg in RPC, rival in TOF", "tagged pair leg in TOF, rival in RPC"};
    for (int d = 0; d < 4; ++d) { const double p = tot[d] > 0 ? win[d] / tot[d] : 0; printf("%s %-40s P(tagged pair M<0.03 wins) = %.4f +- %.4f (%.0f events)\n", mode, dn[d], p, tot[d] > 0 ? std::sqrt(p * (1 - p) / tot[d]) : 0., tot[d]); }
}
