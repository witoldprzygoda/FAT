// isbest_truth_sim.C
//
// How well does the PAT best-pair choice pick the GENUINE e+e- pair?
// Simulation (pp45_epem_sim output: all e+e- combinations kept, isBest stored),
// where the source of every lepton is known from MC truth.
//
// Two rankings are evaluated on the SAME events:
//   "prod" : isBest as stored. In the June-2026 SMASH productions it is the
//            argmin of the tracking chi2 (hyp_chi2_combined == hyp_chi2rk_rms,
//            RMS of the two legs' Runge-Kutta chi2), checked 100%.
//   "time" : argmin of the reduced relative-time chi2 (hyp_chi2 =
//            ep_hyp_best_chi2_red) among the active hypotheses of the event,
//            first row on ties (PAT strict-minimum rule). This is the ranking of
//            the exp betanew production (isBest == argmin time chi2, 100%).
//   "rand" : uniform random choice among the hypotheses (reference).
//
// Truth classification of each e+e- combination:
//   genuine       : both legs true leptons (sim_id 2 / 3) from the same source
//                   - both primary (parentid -1): same geninfo AND same geninfo2
//                     (decay index in the event; lepton-enriched files)
//                   - both secondary: same parentid and processid and
//                     bitwise-equal creation vertex (GEANT decay, conversion)
//   combinatorial : everything else.
// Events = runs of ADJACENT rows with the same event number within one tree
// (runnumber is 0 in sim, the event counter restarts within a file).
// Duplicate rows: the June-2026 SMASH PAT productions write, in events with
// 2 e+ and 2 e-, only the two disjoint pairs, each twice (rows A,B,A,B; isBest
// on the first copy). Within an event, a row whose (e+, e-) legs repeat an
// earlier row is dropped (first copy kept). The betanew exp production writes
// every combination exactly once.
// Class = number of distinct e+ and e- legs among the event's rows (2,3,4,>=5).
//
// Printed (and histograms in <out_dir>/isbest_truth_sim_<tag>.root):
//  (A) events with >= 2 hypotheses and exactly ONE genuine one: probability that
//      the ranking picks it, per class and per mass of the genuine pair, with
//      the random reference <1/n_hyp>;
//  (B) eps_S(M) / eps_B(M) = fraction of genuine / combinatorial pairs chosen,
//      per class, for M < 0.14, 0.14-0.80, 0.80-1.40 GeV/c^2 (compare with the
//      data emulation: eps_B,c = 0.50 / 0.27 / 0.17 for n = 3 / 4 / >=5, PT3 MID);
//  (C) the data test topology: two hypotheses sharing one leg, one genuine,
//      split by the META system of the non-shared legs (0 = RPC, else TOF) and
//      by the masses (data: genuine low-mass vs combinatorial MID 0.51,
//      RPC/TOF 0.70, TOF/RPC 0.38).
// Counts are unweighted (the question is the ranking, not the rate).
//
// Run: root -l -b -q 'plots/isbest_truth_sim.C+("f1.root,f2.root","pt3==1","tag")'

#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TRandom3.h>
#include <TString.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <cmath>
#include <cstdio>
#include <vector>

namespace itsim {
constexpr int kNb = 58;
const double* edges() {
    static double e[kNb + 1]; static bool b = false;
    if (!b) { int k = 0; for (int i = 0; i <= 20; ++i) e[k++] = i * 0.0025;
        e[k++] = 0.08; e[k++] = 0.11; e[k++] = 0.14;
        for (int i = 1; i <= 33; ++i) e[k++] = 0.14 + i * 0.02;
        e[k++] = 0.95; e[k++] = 1.40; b = true; }
    return e;
}
struct Row { float m, best, chi2t, sys_p, sys_m; float lp[3], lm[3]; bool genuine; };
constexpr int kNR = 3;                       // rankings
const char* kRName[kNR] = {"prod", "time", "rand"};
} // namespace itsim

void isbest_truth_sim(const char* files = "/hdd3/przygoda/HADES/FAT_out/sim/output_epem_sim_std_jun26.root",
                      const char* trig_sel = "pt3==1",
                      const char* tag = "std_pt3",
                      const char* out_dir = "plots/output",
                      int seed = 4357)
{
    using namespace itsim;
    TRandom3 rnd(seed);
    TChain ch("dilepton_nt");
    TObjArray* fl = TString(files).Tokenize(",");
    for (int i = 0; i < fl->GetEntries(); ++i) ch.Add(((TObjString*)fl->At(i))->GetString());
    const bool selPT3 = TString(trig_sel) == "pt3==1", selPT2 = TString(trig_sel) == "pt2==1";

    Float_t m, best, chi2t, pt3, pt2, ev, sp, sm, pp, pth, pph, mp, mth, mph;
    Float_t idp, idm, parp, parm, prp, prm, gp, gm, g2p, g2m, vxp, vyp, vzp, vxm, vym, vzm;
    struct B { const char* n; Float_t* a; } br[] = {
        {"m_ee", &m}, {"isBest", &best}, {"hyp_chi2", &chi2t}, {"pt3", &pt3}, {"pt2", &pt2}, {"event", &ev},
        {"ep_system", &sp}, {"em_system", &sm}, {"ep_p_rec", &pp}, {"ep_theta", &pth}, {"ep_phi", &pph},
        {"em_p_rec", &mp}, {"em_theta", &mth}, {"em_phi", &mph},
        {"ep_sim_id", &idp}, {"em_sim_id", &idm}, {"ep_sim_parentid", &parp}, {"em_sim_parentid", &parm},
        {"ep_sim_processid", &prp}, {"em_sim_processid", &prm}, {"ep_sim_geninfo", &gp}, {"em_sim_geninfo", &gm},
        {"ep_sim_geninfo2", &g2p}, {"em_sim_geninfo2", &g2m},
        {"ep_sim_vertexx", &vxp}, {"ep_sim_vertexy", &vyp}, {"ep_sim_vertexz", &vzp},
        {"em_sim_vertexx", &vxm}, {"em_sim_vertexy", &vym}, {"em_sim_vertexz", &vzm}};
    ch.SetBranchStatus("*", 0);
    for (auto& b : br) { ch.SetBranchStatus(b.n, 1); ch.SetBranchAddress(b.n, b.a); }

    // (B) spectra: [class][genuine/comb] all, and chosen per ranking
    TH1D *hAll[4][2], *hSel[kNR][4][2];
    const char* tn[2] = {"gen", "comb"};
    for (int c = 0; c < 4; ++c) for (int g = 0; g < 2; ++g) {
        hAll[c][g] = new TH1D(Form("h_all_%s_c%d", tn[g], c + 2), ";M_{ee} [GeV/c^{2}];pairs", kNb, edges());
        for (int r = 0; r < kNR; ++r)
            hSel[r][c][g] = new TH1D(Form("h_%s_%s_c%d", kRName[r], tn[g], c + 2), ";M_{ee} [GeV/c^{2}];pairs", kNb, edges());
    }
    // (A) events with >=2 hypotheses, exactly one genuine: [class][mass bin of genuine: <0.14, 0.14-0.80, >=0.80]
    double aN[4][3] = {}, aW[kNR][4][3] = {}, aInv[4][3] = {};
    // events with >=2 hypotheses and >=1 genuine (any number): P(chosen is genuine)
    double gN[4] = {}, gW[kNR][4] = {};
    // (C) data-test topology: [sys gen][sys comb] and [mass gen][mass comb]
    double cN[2][2] = {}, cW[kNR][2][2] = {}, mN[2][2] = {}, mW[kNR][2][2] = {};
    long nEv = 0, nEvMulti = 0, nProdBad = 0, nTie = 0, nAgree = 0, nDupRows = 0, nDupEv = 0;

    auto same3 = [](const float* a, const float* b) { return a[0] == b[0] && a[1] == b[1] && a[2] == b[2]; };
    std::vector<Row> g; float curEv = -1; int curTree = -1;
    auto flush = [&]() {
        if (g.empty()) return;
        {   // drop repeated (e+, e-) hypotheses, keep the first copy
            std::vector<Row> u; bool dup = false;
            for (auto& r : g) {
                bool seen = false;
                for (auto& x : u) if (same3(x.lp, r.lp) && same3(x.lm, r.lm)) { seen = true; break; }
                if (seen) { ++nDupRows; dup = true; } else u.push_back(r);
            }
            if (dup) ++nDupEv;
            g.swap(u);
        }
        ++nEv;
        const int nh = (int)g.size();
        std::vector<const float*> P, M;
        for (auto& r : g) {
            bool fp = false, fm = false;
            for (auto* x : P) if (same3(x, r.lp)) fp = true;
            for (auto* x : M) if (same3(x, r.lm)) fm = true;
            if (!fp) P.push_back(r.lp);
            if (!fm) M.push_back(r.lm);
        }
        const int n = (int)(P.size() + M.size());
        const int c = n >= 5 ? 3 : std::max(0, n - 2);
        // chosen index per ranking
        int sel[kNR] = {-1, 0, 0};
        int nb = 0; for (int i = 0; i < nh; ++i) if (g[i].best > 0.5f) { sel[0] = i; ++nb; }
        if (nb != 1) { ++nProdBad; g.clear(); return; }
        for (int i = 1; i < nh; ++i) if (g[i].chi2t < g[sel[1]].chi2t) sel[1] = i;
        { int nmin = 0; for (int i = 0; i < nh; ++i) if (g[i].chi2t == g[sel[1]].chi2t) ++nmin; if (nmin > 1 && nh > 1) ++nTie; }
        sel[2] = nh > 1 ? int(rnd.Integer(nh)) : 0;
        if (nh > 1) { ++nEvMulti; if (sel[0] == sel[1]) ++nAgree; }
        int ngen = 0, igen = -1; for (int i = 0; i < nh; ++i) if (g[i].genuine) { ++ngen; igen = i; }
        for (int i = 0; i < nh; ++i) {
            const int k = g[i].genuine ? 0 : 1;
            hAll[c][k]->Fill(g[i].m);
            for (int r = 0; r < kNR; ++r) if (sel[r] == i) hSel[r][c][k]->Fill(g[i].m);
        }
        if (nh >= 2 && ngen >= 1) {
            gN[c] += 1; for (int r = 0; r < kNR; ++r) gW[r][c] += g[sel[r]].genuine ? 1 : 0;
        }
        if (nh >= 2 && ngen == 1) {
            const int mb = g[igen].m < 0.14 ? 0 : (g[igen].m < 0.80 ? 1 : 2);
            aN[c][mb] += 1; aInv[c][mb] += 1.0 / nh;
            for (int r = 0; r < kNR; ++r) aW[r][c][mb] += (sel[r] == igen);
        }
        if (nh == 2 && ngen == 1) {
            const bool shp = same3(g[0].lp, g[1].lp), shm = same3(g[0].lm, g[1].lm);
            if (shp != shm) {
                const int t = igen, f = 1 - t;
                const float st = shm ? g[t].sys_p : g[t].sys_m;
                const float sf = shm ? g[f].sys_p : g[f].sys_m;
                const int a = st == 0 ? 0 : 1, b = sf == 0 ? 0 : 1;
                const int ma = g[t].m < 0.14 ? 0 : 1, mb = g[f].m < 0.14 ? 0 : 1;
                cN[a][b] += 1; mN[ma][mb] += 1;
                for (int r = 0; r < kNR; ++r) { cW[r][a][b] += (sel[r] == t); mW[r][ma][mb] += (sel[r] == t); }
            }
        }
        g.clear();
    };

    const Long64_t N = ch.GetEntries();
    for (Long64_t i = 0; i < N; ++i) {
        ch.GetEntry(i);
        if (ev != curEv || ch.GetTreeNumber() != curTree) { flush(); curEv = ev; curTree = ch.GetTreeNumber(); }
        if (selPT3 && pt3 < 0.5f) continue;
        if (selPT2 && pt2 < 0.5f) continue;
        Row r; r.m = m; r.best = best; r.chi2t = chi2t; r.sys_p = sp; r.sys_m = sm;
        r.lp[0] = pp; r.lp[1] = pth; r.lp[2] = pph; r.lm[0] = mp; r.lm[1] = mth; r.lm[2] = mph;
        bool gen = false;
        if (idp == 2 && idm == 3) {
            if (parp == -1 && parm == -1) gen = (gp == gm && g2p == g2m);
            else if (parp == parm && parp > 0 && prp == prm && vxp == vxm && vyp == vym && vzp == vzm) gen = true;
        }
        r.genuine = gen;
        g.push_back(r);
        if (i % 50000000 == 0 && i) printf("  %lld / %lld rows\n", i, N);
    }
    flush();

    printf("\n=== best-pair choice vs MC truth [%s], %s: %lld rows, %ld events, %ld with >=2 hypotheses ===\n", tag, trig_sel, N, nEv, nEvMulti);
    printf("duplicate rows dropped: %ld in %ld events | events without exactly one isBest==1 row (skipped): %ld | time-chi2 ties at the minimum: %ld | prod and time pick the same pair: %.3f\n",
           nDupRows, nDupEv, nProdBad, nTie, nEvMulti ? double(nAgree) / nEvMulti : 0.);

    const char* cl[4] = {"n=2 ", "n=3 ", "n=4 ", "n>=5"};
    const char* mbn[3] = {"M<0.14", "0.14-0.80", "M>=0.80"};
    printf("\n(A) events with >=2 hypotheses and exactly one genuine: P(genuine chosen) for prod / time / random  [<1/n_hyp>] (events)\n");
    for (int c = 0; c < 4; ++c) for (int mb = 0; mb < 3; ++mb) {
        if (aN[c][mb] < 1) continue;
        printf("   %s  genuine %-9s : %.3f / %.3f / %.3f   [%.3f]  (%.0f)\n", cl[c], mbn[mb],
               aW[0][c][mb] / aN[c][mb], aW[1][c][mb] / aN[c][mb], aW[2][c][mb] / aN[c][mb], aInv[c][mb] / aN[c][mb], aN[c][mb]);
    }
    {   double n = 0, w0 = 0, w1 = 0, w2 = 0, iv = 0;
        for (int c = 0; c < 4; ++c) for (int mb = 0; mb < 3; ++mb) { n += aN[c][mb]; w0 += aW[0][c][mb]; w1 += aW[1][c][mb]; w2 += aW[2][c][mb]; iv += aInv[c][mb]; }
        if (n > 0) printf("   all                     : %.3f / %.3f / %.3f   [%.3f]  (%.0f)\n", w0 / n, w1 / n, w2 / n, iv / n, n); }
    printf("    events with >=2 hypotheses and >=1 genuine: P(chosen pair is genuine) prod / time / random\n");
    for (int c = 0; c < 4; ++c) if (gN[c] > 0) printf("   %s : %.3f / %.3f / %.3f  (%.0f)\n", cl[c], gW[0][c] / gN[c], gW[1][c] / gN[c], gW[2][c] / gN[c], gN[c]);

    const double rng[3][2] = {{0, 0.14}, {0.14, 0.80}, {0.80, 1.40}};
    printf("\n(B) fraction of genuine (eps_S) and combinatorial (eps_B) pairs chosen: prod / time / random\n");
    for (int ri = 0; ri < 3; ++ri) {
        printf("  M %.2f-%.2f GeV/c2\n", rng[ri][0], rng[ri][1]);
        double tA[2] = {0, 0}, tS[kNR][2] = {};
        for (int c = 0; c < 4; ++c) {
            const int b1 = hAll[c][0]->FindBin(rng[ri][0] + 1e-6), b2 = hAll[c][0]->FindBin(rng[ri][1] - 1e-6);
            double a[2], s[kNR][2];
            for (int k = 0; k < 2; ++k) { a[k] = hAll[c][k]->Integral(b1, b2); tA[k] += a[k];
                for (int r = 0; r < kNR; ++r) { s[r][k] = hSel[r][c][k]->Integral(b1, b2); tS[r][k] += s[r][k]; } }
            if (a[0] + a[1] < 1) continue;
            printf("   %s  genuine %9.0f: eps_S %.3f / %.3f / %.3f | combinatorial %9.0f: eps_B %.3f / %.3f / %.3f | genuine share all %.3f -> prod %.3f, time %.3f\n",
                   cl[c], a[0], a[0] ? s[0][0] / a[0] : 0., a[0] ? s[1][0] / a[0] : 0., a[0] ? s[2][0] / a[0] : 0.,
                   a[1], a[1] ? s[0][1] / a[1] : 0., a[1] ? s[1][1] / a[1] : 0., a[1] ? s[2][1] / a[1] : 0.,
                   a[0] / (a[0] + a[1]), (s[0][0] + s[0][1]) ? s[0][0] / (s[0][0] + s[0][1]) : 0., (s[1][0] + s[1][1]) ? s[1][0] / (s[1][0] + s[1][1]) : 0.);
        }
        printf("   all   genuine %9.0f: eps_S %.3f / %.3f / %.3f | combinatorial %9.0f: eps_B %.3f / %.3f / %.3f\n",
               tA[0], tA[0] ? tS[0][0] / tA[0] : 0., tA[0] ? tS[1][0] / tA[0] : 0., tA[0] ? tS[2][0] / tA[0] : 0.,
               tA[1], tA[1] ? tS[0][1] / tA[1] : 0., tA[1] ? tS[1][1] / tA[1] : 0., tA[1] ? tS[2][1] / tA[1] : 0.);
    }

    printf("\n(C) two hypotheses sharing one leg, one genuine: P(genuine chosen) prod / time / random (events)\n");
    const char* S[2] = {"RPC", "TOF"};
    for (int a = 0; a < 2; ++a) for (int b = 0; b < 2; ++b) if (cN[a][b] > 0)
        printf("   non-shared leg: genuine %s, combinatorial %s : %.3f / %.3f / %.3f  (%.0f)\n", S[a], S[b],
               cW[0][a][b] / cN[a][b], cW[1][a][b] / cN[a][b], cW[2][a][b] / cN[a][b], cN[a][b]);
    const char* L[2] = {"M<0.14", "M>=0.14"};
    for (int a = 0; a < 2; ++a) for (int b = 0; b < 2; ++b) if (mN[a][b] > 0)
        printf("   genuine %-7s vs combinatorial %-7s : %.3f / %.3f / %.3f  (%.0f)\n", L[a], L[b],
               mW[0][a][b] / mN[a][b], mW[1][a][b] / mN[a][b], mW[2][a][b] / mN[a][b], mN[a][b]);

    TFile fo(Form("%s/isbest_truth_sim_%s.root", out_dir, tag), "RECREATE");
    for (int c = 0; c < 4; ++c) for (int k = 0; k < 2; ++k) { hAll[c][k]->Write(); for (int r = 0; r < kNR; ++r) hSel[r][c][k]->Write(); }
    fo.Close();
    printf("\nwritten %s/isbest_truth_sim_%s.root\n", out_dir, tag);
}
