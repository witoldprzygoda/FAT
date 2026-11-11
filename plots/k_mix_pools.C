// k_mix_pools.C
//
// Event-mixing k factor, k_mix(M) = N+-_mix / 2 sqrt(N++_mix N--_mix), with the SAME code on data and on
// simulation, for two truth-free pools of leptons:
//   all  : every event with at least one e+e- hypothesis,
//   nge3 : only events with at least three identified leptons (n+ + n- >= 3); their legs are mostly unpaired,
//          i.e. they resemble the legs that really make the combinatorial background.
// Rolling buffer of the last pool_depth events of the pool, no multiplicity classes, legs unique by float
// equality of (p, theta, phi), events = adjacent rows of one (run, event) or (event, reconstructed vertex).
// A second, truth-free measurement on the same events: in events with exactly 3 leptons the two e+e- hypotheses
// share one leg; if the two other legs are in different detectors (theta < 45 deg taken as RPC), how often does
// the hypothesis with the RPC leg get isBest == 1?
//
// Input : data = output_epem_exp.root, tree dilepton_nt (columns p_rec, PT3 = trigbit 8192; vertex cut done by FAT)
//         sim  = smash_std_sim_epem.root, tree EpEm_ID (columns p; vertex cut eVertReco_z > -500 here)
// Output: root file with h_<pool>_<set>_{pm,pp,mm}, set = PT3, PT2 (data; trigbit 8192, 4096) or all, PT3, even, odd
//         (simulation; PT3 = bit 13, even/odd = parity of the event number, for an independent closure test), and a
//         printed table. Simulation only: truth spectra h_truth_<set>_{B,S,N,pp,mm} (no pool, all combinations).
//         Pools per lepton-multiplicity class: h_c2_, h_c3_, h_c4_ (n = 2, 3, >= 4). PT3 events only, in (M, pT of the pair):
//         h2_<pool>_{pm,pp,mm} for the pools all and nge3 and the same-event h2_truth_{B,S,N,pp,mm} (data: N, pp, mm);
//         hf_<pool>_{pm,pp,mm}: the same pools in 5 MeV bins. oa_min > 0 keeps only mixed pairs with a larger opening angle.
//
// Run: root -l -b -q 'plots/k_mix_pools.C+("data","output_epem_exp.root","plots/output/k_mix_pools_data.root")'
//      root -l -b -q 'plots/k_mix_pools.C+("sim","/hdd3/przygoda/HADES/files/sim/smash_std_sim_epem.root","plots/output/k_mix_pools_sim.root")'

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TStopwatch.h>
#include <TString.h>
#include <TTree.h>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace kmp {
constexpr int kNb = 58;
const double* edges() {
    static double e[kNb + 1]; static bool b = false;
    if (!b) { int k = 0; for (int i = 0; i <= 20; ++i) e[k++] = i * 0.0025;
        e[k++] = 0.08; e[k++] = 0.11; e[k++] = 0.14;
        for (int i = 1; i <= 33; ++i) e[k++] = 0.14 + i * 0.02;
        e[k++] = 0.95; e[k++] = 1.40; b = true; }
    return e;
}
int massBin(double m) {
    const double* e = edges();
    if (m < e[0] || m >= e[kNb]) return 0;
    int lo = 0, hi = kNb;
    while (hi - lo > 1) { const int mid = (lo + hi) / 2; if (m >= e[mid]) lo = mid; else hi = mid; }
    return lo + 1;
}
constexpr double kMe = 0.51099895;
struct MLeg { double ux, uy, uz, e, p; };
struct MEv { std::vector<MLeg> ep, em; };
inline MLeg mk(Float_t p, Float_t th, Float_t ph) {
    const double d = TMath::DegToRad(), t = th * d, f = ph * d;
    return {std::sin(t) * std::cos(f), std::sin(t) * std::sin(f), std::cos(t), std::sqrt(double(p) * p + kMe * kMe), double(p)};
}
constexpr int kNM2 = 7, kNP2 = 6;
const double kM2[kNM2 + 1] = {0.05, 0.14, 0.16, 0.20, 0.26, 0.34, 0.56, 0.80};   // GeV/c2
const double kP2[kNP2 + 1] = {0., 100., 200., 300., 400., 600., 3000.};          // MeV/c
inline double pairPt(const MLeg& a, const MLeg& b) { return std::hypot(a.p * a.ux + b.p * b.ux, a.p * a.uy + b.p * b.uy); }
TH2D* make2(const TString& n) { return new TH2D(n, ";M_{ee} [GeV/c^{2}];p_{T} [MeV/c]", kNM2, kM2, kNP2, kP2); }
struct Pool {
    TH1D* h[3]; TH2D* h2[3] = {nullptr, nullptr, nullptr}; TH1D* hf[3] = {nullptr, nullptr, nullptr}; double cosOA = 2.; std::vector<MEv> buf; size_t head = 0, filled = 0; long nev = 0;
    void init(const TString& name, int depth) { const char* pn[3] = {"pm", "pp", "mm"}; for (int p = 0; p < 3; ++p) h[p] = new TH1D(Form("h_%s_%s", name.Data(), pn[p]), ";M_{ee} [GeV/c^{2}];mixed pairs", kNb, edges()); buf.resize(depth); }
    void initFine(const TString& name) { const char* pn[3] = {"pm", "pp", "mm"}; for (int p = 0; p < 3; ++p) hf[p] = new TH1D(Form("hf_%s_%s", name.Data(), pn[p]), ";M_{ee} [GeV/c^{2}];mixed pairs", 280, 0., 1.4); }
    void init2(const TString& name) { const char* pn[3] = {"pm", "pp", "mm"}; for (int p = 0; p < 3; ++p) h2[p] = make2(Form("h2_%s_%s", name.Data(), pn[p])); }
    inline void add(int p, const MLeg& a, const MLeg& b) { const double cs = a.ux * b.ux + a.uy * b.uy + a.uz * b.uz; if (cs >= cosOA) return; const double m2 = 2 * kMe * kMe + 2 * (a.e * b.e - a.p * b.p * cs); const double m = m2 > 0 ? std::sqrt(m2) / 1000. : 0.; const int bin = massBin(m); if (bin > 0) h[p]->AddBinContent(bin); if (h2[p]) h2[p]->Fill(m, pairPt(a, b)); if (hf[p]) hf[p]->Fill(m); }
    void mix(const MEv& E) {
        for (size_t ib = 0; ib < filled; ++ib) { const MEv& B = buf[ib];
            for (auto& a : E.ep) for (auto& c : B.em) add(0, a, c);
            for (auto& a : E.em) for (auto& c : B.ep) add(0, a, c);
            for (auto& a : E.ep) for (auto& c : B.ep) add(1, a, c);
            for (auto& a : E.em) for (auto& c : B.em) add(2, a, c); }
        buf[head] = E; head = (head + 1) % buf.size(); if (filled < buf.size()) ++filled; ++nev;
    }
};
struct Row { Float_t pp, tp, fp, pm, tm, fm, best; };
}

void k_mix_pools(const char* mode = "data", const char* file = "output_epem_exp.root", const char* out = "plots/output/k_mix_pools_data.root", int pool_depth = 10, double vertex_z_min = -500., double oa_min = 0.)
{
    using namespace kmp;
    const bool sim = TString(mode) == "sim";
    TStopwatch sw; sw.Start();
    TFile f(file); TTree* t = (TTree*)f.Get(sim ? "EpEm_ID" : "dilepton_nt"); if (!t) { printf("no tree in %s\n", file); return; }
    t->SetBranchStatus("*", 0);
    auto sb = [&](const char* n, Float_t* a) { if (!t->GetBranch(n)) { printf("missing branch %s\n", n); return; } t->SetBranchStatus(n, 1); t->SetBranchAddress(n, a); };
    Float_t run = 0, ev = 0, vx = 0, vy = 0, vz = 0, trig = 0; Row r;
    Float_t ida = 0, idb = 0, para = 0, parb = 0, proa = 0, prob = 0, sxa = 0, sya = 0, sza = 0, sxb = 0, syb = 0, szb = 0;
    sb("event", &ev); sb("trigbit", &trig); sb("isBest", &r.best);
    if (sim) { sb("eVertReco_x", &vx); sb("eVertReco_y", &vy); sb("eVertReco_z", &vz); sb("ep_p", &r.pp); sb("em_p", &r.pm);
        sb("ep_sim_id", &ida); sb("em_sim_id", &idb); sb("ep_sim_parentid", &para); sb("em_sim_parentid", &parb); sb("ep_sim_processid", &proa); sb("em_sim_processid", &prob);
        sb("ep_sim_vertexx", &sxa); sb("ep_sim_vertexy", &sya); sb("ep_sim_vertexz", &sza); sb("em_sim_vertexx", &sxb); sb("em_sim_vertexy", &syb); sb("em_sim_vertexz", &szb); }
    else { sb("runnumber", &run); sb("ep_p_rec", &r.pp); sb("em_p_rec", &r.pm); }
    sb("ep_theta", &r.tp); sb("ep_phi", &r.fp); sb("em_theta", &r.tm); sb("em_phi", &r.fm);

    TFile fo(out, "RECREATE");
    const int nSets = sim ? 4 : 2; const char* setName[4] = {sim ? "all" : "PT3", sim ? "PT3" : "PT2", "even", "odd"};
    auto inSet = [&](int s, Float_t tr, Float_t evn) { const int b = int(tr);
        if (!sim) return s == 0 ? b == 8192 : b == 4096;
        if (s == 0) return true; if (s == 1) return (b & 8192) != 0; return ((long(evn) & 1) == 0) == (s == 2); };
    TH1D* hT[4][5]; if (sim) for (int s = 0; s < nSets; ++s) { const char* q[5] = {"B", "S", "N", "pp", "mm"}; for (int i = 0; i < 5; ++i) hT[s][i] = new TH1D(Form("h_truth_%s_%s", setName[s], q[i]), ";M_{ee} [GeV/c^{2}];pairs", kNb, edges()); }
    auto fillT = [&](int kind, Float_t tr, Float_t evn, double m) { for (int s = 0; s < nSets; ++s) if (inSet(s, tr, evn)) hT[s][kind]->Fill(m); };
    auto pairM = [](const MLeg& a, const MLeg& b) { const double cs = a.ux * b.ux + a.uy * b.uy + a.uz * b.uz; const double m2 = 2 * kMe * kMe + 2 * (a.e * b.e - a.p * b.p * cs); return m2 > 0 ? std::sqrt(m2) / 1000. : 0.; };
    Pool pAll[4], pN3[4], pC[4][3];
    for (int s = 0; s < nSets; ++s) { pAll[s].init(Form("all_%s", setName[s]), pool_depth); pN3[s].init(Form("nge3_%s", setName[s]), pool_depth); for (int c = 0; c < 3; ++c) pC[s][c].init(Form("c%d_%s", c + 2, setName[s]), pool_depth); }
    const int sPT3 = sim ? 1 : 0; pAll[sPT3].init2("all"); pN3[sPT3].init2("nge3"); pAll[sPT3].initFine("all"); pN3[sPT3].initFine("nge3");
    if (oa_min > 0) { const double c0 = std::cos(oa_min * TMath::DegToRad()); for (int s = 0; s < nSets; ++s) { pAll[s].cosOA = c0; pN3[s].cosOA = c0; for (int c = 0; c < 3; ++c) pC[s][c].cosOA = c0; } }   // mixed pairs only
    TH2D* hT2[5]; { const char* q[5] = {"B", "S", "N", "pp", "mm"}; for (int i = 0; i < 5; ++i) hT2[i] = make2(Form("h2_truth_%s", q[i])); }
    auto fill2 = [&](int kind, const MLeg& a, const MLeg& b) { hT2[kind]->Fill(pairM(a, b), pairPt(a, b)); };
    long nEv[4] = {0, 0, 0, 0}, nN3[4] = {0, 0, 0, 0}, rpcWin[4] = {0, 0, 0, 0}, rpcTot[4] = {0, 0, 0, 0}, sameDet[4] = {0, 0, 0, 0}, sameTot[4] = {0, 0, 0, 0};

    std::vector<Row> g; bool have = false; Float_t c1 = 0, c2 = 0, c3 = 0, c4 = 0, evTrig = 0, evNum = 0;
    auto sameLegP = [](const Row& a, const Row& b) { return a.pp == b.pp && a.tp == b.tp && a.fp == b.fp; };
    auto sameLegM = [](const Row& a, const Row& b) { return a.pm == b.pm && a.tm == b.tm && a.fm == b.fm; };
    auto flush = [&]() {
        if (g.empty()) return;
        std::vector<const Row*> P, M;
        for (auto& x : g) {
            bool fp = false, fm = false;
            for (auto* l : P) if (sameLegP(*l, x)) fp = true;
            for (auto* l : M) if (sameLegM(*l, x)) fm = true;
            if (!fp) P.push_back(&x);
            if (!fm) M.push_back(&x);
        }
        MEv E; for (auto* l : P) E.ep.push_back(mk(l->pp, l->tp, l->fp)); for (auto* l : M) E.em.push_back(mk(l->pm, l->tm, l->fm));
        const int n = (int)P.size() + (int)M.size();
        for (int s = 0; s < nSets; ++s) { if (!inSet(s, evTrig, evNum)) continue;
            ++nEv[s]; pAll[s].mix(E); if (n >= 3) { ++nN3[s]; pN3[s].mix(E); } pC[s][n <= 2 ? 0 : (n == 3 ? 1 : 2)].mix(E);
            if (n == 3 && g.size() == 2) {   // two hypotheses sharing one leg: which of the two non-shared legs is in RPC?
                const bool shareE = sameLegM(g[0], g[1]);    // e- shared: compare the e+ legs; otherwise compare the e- legs
                const double th0 = shareE ? g[0].tp : g[0].tm, th1 = shareE ? g[1].tp : g[1].tm;
                const bool r0 = th0 < 45., r1 = th1 < 45.;
                if (r0 != r1) { ++rpcTot[s]; if ((r0 && g[0].best > 0.5f) || (r1 && g[1].best > 0.5f)) ++rpcWin[s]; }
                else { ++sameTot[s]; if (g[0].best > 0.5f) ++sameDet[s]; }
            }
        }
        g.clear();
    };
    const Long64_t N = t->GetEntries();
    for (Long64_t i = 0; i < N; ++i) {
        t->GetEntry(i);
        const Float_t a1 = sim ? vx : run, a2 = sim ? vy : 0.f, a3 = sim ? vz : 0.f;
        if (!have || ev != c1 || a1 != c2 || a2 != c3 || a3 != c4) { flush(); c1 = ev; c2 = a1; c3 = a2; c4 = a3; have = true; evTrig = trig; evNum = ev; }
        if (sim && vz <= vertex_z_min) continue;
        g.push_back(r);
        if (inSet(sPT3, trig, ev)) { const MLeg a = mk(r.pp, r.tp, r.fp), b = mk(r.pm, r.tm, r.fm); fill2(2, a, b);
            if (sim) { const bool la = ida == 2 || ida == 3, lb = idb == 2 || idb == 3, same = sxa == sxb && sya == syb && sza == szb && para == parb && proa == prob;
                if (la && lb && same && ida != idb) fill2(1, a, b); else if (!(la && lb && same)) fill2(0, a, b); } }
        if (sim) {   // truth: background = neither genuine nor one particle twice
            const double m = pairM(mk(r.pp, r.tp, r.fp), mk(r.pm, r.tm, r.fm));
            const bool la = ida == 2 || ida == 3, lb = idb == 2 || idb == 3, same = sxa == sxb && sya == syb && sza == szb && para == parb && proa == prob;
            const bool gen = la && lb && same && ida != idb, twice = la && lb && same && ida == idb;
            fillT(2, trig, ev, m); if (gen) fillT(1, trig, ev, m); else if (!twice) fillT(0, trig, ev, m);
        }
        if (i % 10000000 == 0 && i) { printf("  row %lld / %lld (%.0f s)\n", i, N, sw.RealTime()); sw.Continue(); }
    }
    flush();
    if (sim) for (int q = 0; q < 2; ++q) {   // like-sign rows of the same sets
        TTree* tl = (TTree*)f.Get(q == 0 ? "EpEp_ID" : "EmEm_ID"); tl->SetBranchStatus("*", 0); Float_t e2, z2, tr2, p1, t1, f1, p2, t2, f2;
        auto sl = [&](const std::string& n, Float_t* a) { tl->SetBranchStatus(n.c_str(), 1); tl->SetBranchAddress(n.c_str(), a); }; const std::string a = q == 0 ? "ep1_" : "em1_", b = q == 0 ? "ep2_" : "em2_";
        sl("event", &e2); sl("eVertReco_z", &z2); sl("trigbit", &tr2); sl(a + "p", &p1); sl(a + "theta", &t1); sl(a + "phi", &f1); sl(b + "p", &p2); sl(b + "theta", &t2); sl(b + "phi", &f2);
        for (Long64_t i = 0; i < tl->GetEntries(); ++i) { tl->GetEntry(i); if (z2 <= vertex_z_min) continue; const MLeg a2 = mk(p1, t1, f1), b2 = mk(p2, t2, f2); fillT(3 + q, tr2, e2, pairM(a2, b2)); if (inSet(sPT3, tr2, e2)) fill2(3 + q, a2, b2); }
    }
    if (!sim) for (int q = 0; q < 2; ++q) {   // like-sign rows of the data, PT3, for the (M, pT) formula; in these files ep_ and em_ are the two legs of the same charge
        TString fn(file); fn.ReplaceAll("output_epem_exp", q == 0 ? "output_epep_exp" : "output_emem_exp"); TFile fl(fn); TTree* tl = (TTree*)fl.Get("dilepton_nt"); if (!tl) { printf("no like-sign tree in %s\n", fn.Data()); continue; }
        tl->SetBranchStatus("*", 0); Float_t tr2, p1, t1, f1, p2, t2, f2; auto sl = [&](const char* n, Float_t* a) { tl->SetBranchStatus(n, 1); tl->SetBranchAddress(n, a); };
        sl("trigbit", &tr2); sl("ep_p_rec", &p1); sl("ep_theta", &t1); sl("ep_phi", &f1); sl("em_p_rec", &p2); sl("em_theta", &t2); sl("em_phi", &f2);
        for (Long64_t i = 0; i < tl->GetEntries(); ++i) { tl->GetEntry(i); if (tr2 != 8192.f) continue; fill2(3 + q, mk(p1, t1, f1), mk(p2, t2, f2)); }
    }
    for (int s = 0; s < nSets; ++s) {
        printf("\nset %s: events %ld, n>=3 events %ld (%.1f %%); pool depth %d\n", setName[s], nEv[s], nN3[s], 100. * nN3[s] / std::max(1L, nEv[s]), pool_depth);
        printf("n=3, non-shared legs in different detectors: hypothesis with the RPC leg wins %ld / %ld = %.4f +- %.4f; same detector: first row wins %ld / %ld = %.4f\n",
               rpcWin[s], rpcTot[s], double(rpcWin[s]) / std::max(1L, rpcTot[s]), std::sqrt(0.25 / std::max(1L, rpcTot[s])), sameDet[s], sameTot[s], double(sameDet[s]) / std::max(1L, sameTot[s]));
        const double W[7][2] = {{0, 0.05}, {0.05, 0.14}, {0.14, 0.20}, {0.20, 0.30}, {0.30, 0.56}, {0.56, 0.80}, {0.14, 0.80}};
        for (Pool* p : {&pAll[s], &pN3[s]}) {
            printf(" pool %s\n", p->h[0]->GetName());
            for (auto& w : W) { const int b1 = p->h[0]->FindBin(w[0] + 1e-6), b2 = p->h[0]->FindBin(w[1] - 1e-6); double a = 0, g2 = 0; for (int b = b1; b <= b2; ++b) { const double x = p->h[1]->GetBinContent(b), y = p->h[2]->GetBinContent(b); a += p->h[0]->GetBinContent(b); if (x > 0 && y > 0) g2 += 2 * std::sqrt(x * y); } printf("   %.2f-%.2f  k_mix = %.4f  (pm %.0f)\n", w[0], w[1], g2 > 0 ? a / g2 : 0., a); }
        }
    }
    fo.cd();   // opening the like-sign files moved the current directory
    for (int s = 0; s < nSets; ++s) { TH1D* hr = new TH1D(Form("h_rpc_%s", setName[s]), "n=3, theta<45 taken as RPC: RPC-leg wins, RPC-leg total, first row wins, same-detector total", 4, 0, 4); hr->SetBinContent(1, rpcWin[s]); hr->SetBinContent(2, rpcTot[s]); hr->SetBinContent(3, sameDet[s]); hr->SetBinContent(4, sameTot[s]); }
    fo.Write(); fo.Close();
    printf("\nwritten %s (%.0f s)\n", out, sw.RealTime());
}
