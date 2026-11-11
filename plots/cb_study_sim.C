// cb_study_sim.C
//
// Combinatorial-background study on a full SMASH simulation with the PAT trees
// EpEm_ID, EpEp_ID and EmEm_ID (betanew-like production: time-chi2 ranking,
// every combination written once, no isBest = -1 rows, no generator weights).
//
// Every pair is classified with the MC truth:
//   genuine : both legs true leptons of opposite truth charge with the same
//             creation vertex (bitwise), parent and process: pi0 Dalitz, eta
//             Dalitz, photon conversion (photon from pi0 / eta / other), other;
//   CB      : true leptons from different sources (Dalitz + Dalitz, Dalitz +
//             conversion, conversion + conversion, with a secondary electron:
//             Compton, delta, conversion of a bremsstrahlung photon; other);
//   fake    : one leg not a true lepton (sharing a RICH-ring GEANT track with
//             the partner or not) or both legs;
//   same    : one particle reconstructed twice (same origin and same truth id).
//
// For sets {all events, PT3}, ranking {isBest >= 0, isBest == 1} and pair cut
// OA > {0, 2, 4, 9} deg (event cut eVertReco_z > -500) it gives on the 58-bin
// analysis grid: N+-, N++, N--, the true background B_true (CB + fake), the
// like-sign estimate CB_hat = 2 sqrt(N++ N--) per bin, signals, the truth k
// factor, the k factor from event mixing done on this same simulation
// (k_mix,sim, pool without multiplicity classes, legs of written hypotheses),
// and the class-weighted best-pair background sum_c eps_B,c CB_c with the
// truth eps_B,c. Multiplicity class = identified e+ and e- legs of the event
// (like-sign events: from EpEm_ID via the key (event, reconstructed vertex),
// own legs if the event has no e+e- hypothesis).
//
// Run: root -l -b -q 'plots/cb_study_sim.C+("file.root","tag")'

#include <TFile.h>
#include <TH1D.h>
#include <TMath.h>
#include <TString.h>
#include <TStopwatch.h>
#include <TTree.h>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace cbst {
constexpr int kNb = 58;
const double* edges() {
    static double e[kNb + 1]; static bool b = false;
    if (!b) { int k = 0; for (int i = 0; i <= 20; ++i) e[k++] = i * 0.0025;
        e[k++] = 0.08; e[k++] = 0.11; e[k++] = 0.14;
        for (int i = 1; i <= 33; ++i) e[k++] = 0.14 + i * 0.02;
        e[k++] = 0.95; e[k++] = 1.40; b = true; }
    return e;
}
int massBin(double m) {   // 1..kNb, 0 outside the grid
    const double* e = edges();
    if (m < e[0] || m >= e[kNb]) return 0;
    int lo = 0, hi = kNb;
    while (hi - lo > 1) { const int mid = (lo + hi) / 2; if (m >= e[mid]) lo = mid; else hi = mid; }
    return lo + 1;
}
constexpr double kMe = 0.51099895;
constexpr int kNT = 2, kNV = 2, kNOA = 4, kNC = 4;
const char* kTN[kNT] = {"all", "PT3"};
const char* kVN[kNV] = {"isBest>=0", "isBest==1"};
const char* kVS[kNV] = {"all", "best"};
const double kOA[kNOA] = {0., 2., 4., 9.};
const char* kPN[3] = {"pm", "pp", "mm"};
const char* kPT[3] = {"e+e-", "e+e+", "e-e-"};

enum Cat { kGenPi0, kGenEta, kGenConvPi0, kGenConvEta, kGenConvOth, kGenOth,
           kDalDal, kDalConv, kConvConv, kSecondary, kOtherLep,
           kFakeShare, kFakeNoShare, kTwoFake, kSamePart, kNCat };
const char* kCatName[kNCat] = {
    "genuine pi0 Dalitz", "genuine eta Dalitz", "genuine conversion, gamma from pi0", "genuine conversion, gamma from eta",
    "genuine conversion, other gamma", "genuine other",
    "CB Dalitz + Dalitz", "CB Dalitz + conversion", "CB conversion + conversion", "CB with secondary e (Compton/delta/brems)", "CB other leptons",
    "fake leg, ring shared with partner", "fake leg, no shared ring", "both legs fake", "one particle twice"};
inline bool isGen(int c) { return c <= kGenOth; }
inline bool isCBLep(int c) { return c >= kDalDal && c <= kOtherLep; }
inline bool isFake(int c) { return c >= kFakeShare && c <= kTwoFake; }
inline bool isB(int c) { return isCBLep(c) || isFake(c); }
inline bool isDalitz(int c) { return c == kGenPi0 || c == kGenEta; }
inline bool isConv(int c) { return c >= kGenConvPi0 && c <= kGenConvOth; }

struct Leg { Float_t p, th, ph, id, par, gp, proc, vx, vy, vz, r1, r2, r3; double ux, uy, uz, e; };
void derive(Leg& l) {
    const double d = TMath::DegToRad();
    l.ux = std::sin(l.th * d) * std::cos(l.ph * d); l.uy = std::sin(l.th * d) * std::sin(l.ph * d); l.uz = std::cos(l.th * d);
    l.e = std::sqrt(double(l.p) * l.p + kMe * kMe);
}
enum Src { kSPi0, kSEta, kSConv, kSSec, kSOther, kSFake };
int source(const Leg& l) {
    if (l.id != 2 && l.id != 3) return kSFake;
    if (l.par == 7) return kSPi0;
    if (l.par == 17) return kSEta;
    if (l.par == 1 && l.proc == 6) return (l.gp == 2 || l.gp == 3) ? kSSec : kSConv;
    if ((l.par == 1 && l.proc == 7) || ((l.par == 2 || l.par == 3) && l.proc == 10)) return kSSec;
    return kSOther;
}
bool shareRing(const Leg& a, const Leg& b) {
    const Float_t ra[3] = {a.r1, a.r2, a.r3}, rb[3] = {b.r1, b.r2, b.r3};
    for (Float_t x : ra) if (x > 0) for (Float_t y : rb) if (x == y) return true;
    return false;
}
int category(const Leg& a, const Leg& b) {
    const int sa = source(a), sb = source(b);
    const int nf = (sa == kSFake) + (sb == kSFake);
    if (nf == 2) return kTwoFake;
    if (nf == 1) return shareRing(a, b) ? kFakeShare : kFakeNoShare;
    if (a.vx == b.vx && a.vy == b.vy && a.vz == b.vz && a.par == b.par && a.proc == b.proc) {
        if (a.id == b.id) return kSamePart;
        if (sa == kSPi0) return kGenPi0;
        if (sa == kSEta) return kGenEta;
        if (a.par == 1 && a.proc == 6) return a.gp == 7 ? kGenConvPi0 : (a.gp == 17 ? kGenConvEta : kGenConvOth);
        return kGenOth;
    }
    if (sa == kSSec || sb == kSSec) return kSecondary;
    if (sa == kSOther || sb == kSOther) return kOtherLep;
    const bool da = (sa == kSPi0 || sa == kSEta), db = (sb == kSPi0 || sb == kSEta);
    if (da && db) return kDalDal;
    return (da || db) ? kDalConv : kConvConv;
}
inline void pairKin(const Leg& a, const Leg& b, double& cs, double& m) {
    cs = a.ux * b.ux + a.uy * b.uy + a.uz * b.uz;
    const double m2 = 2 * kMe * kMe + 2 * (a.e * b.e - double(a.p) * b.p * cs);
    m = m2 > 0 ? std::sqrt(m2) / 1000. : 0.;
}
inline bool sameLeg(const Leg& x, const Leg& y) { return x.p == y.p && x.th == y.th && x.ph == y.ph; }
unsigned long long evKey(Float_t ev, Float_t vx, Float_t vy, Float_t vz) {
    unsigned long long h = 0xcbf29ce484222325ULL;
    auto mix = [&](Float_t x) { unsigned int u; std::memcpy(&u, &x, 4); h ^= u; h *= 0x100000001b3ULL; h ^= (h >> 29); };
    mix(ev); mix(vx); mix(vy); mix(vz);
    return h;
}
struct Row { Leg a, b; Float_t best, chi2; };
struct LSEvent { unsigned long long key; bool pt3; std::vector<Row> rows; };
struct MLeg { double ux, uy, uz, e, p; };
struct MEv { std::vector<MLeg> ep, em; };
int nearestClass(int n) { return n >= 5 ? 3 : (n <= 2 ? 0 : n - 2); }
} // namespace cbst

void cb_study_sim(const char* file = "/hdd3/przygoda/HADES/files/sim/smash_std_sim_epem.root",
                  const char* tag = "std_full",
                  const char* out_dir = "plots/output",
                  int pool_depth = 10,
                  double vertex_z_min = -500.)
{
    using namespace cbst;
    TStopwatch sw; sw.Start();
    TFile* f = TFile::Open(file);
    if (!f || f->IsZombie()) { printf("cannot open %s\n", file); return; }

    // ---- histograms --------------------------------------------------------
    TFile* fo = TFile::Open(Form("%s/cb_study_sim_%s.root", out_dir, tag), "RECREATE");
    TH1D* h[kNT][kNV][kNOA][3][kNCat];
    TH1D* hc[kNT][kNV][kNOA][3][kNC];          // all pairs per multiplicity class
    TH1D* hcB[kNT][kNV][kNOA][kNC];            // true background (CB + fake) of e+e- per class
    TH1D* hmix[kNT][kNOA][3];
    TH1D* hchi[2][kNCat];                      // log10 time chi2, all events, isBest>=0, OA>2: [M<0.14, 0.14-0.80]
    for (int t = 0; t < kNT; ++t) for (int oa = 0; oa < kNOA; ++oa) {
        for (int p = 0; p < 3; ++p) hmix[t][oa][p] = new TH1D(Form("hmix_%s_oa%g_%s", kTN[t], kOA[oa], kPN[p]), ";M_{ee} [GeV/c^{2}];mixed pairs", kNb, edges());
        for (int v = 0; v < kNV; ++v) {
            for (int p = 0; p < 3; ++p) {
                for (int c = 0; c < kNCat; ++c) h[t][v][oa][p][c] = new TH1D(Form("h_%s_%s_oa%g_%s_cat%d", kTN[t], kVS[v], kOA[oa], kPN[p], c), ";M_{ee} [GeV/c^{2}];pairs", kNb, edges());
                for (int c = 0; c < kNC; ++c) hc[t][v][oa][p][c] = new TH1D(Form("hc_%s_%s_oa%g_%s_c%d", kTN[t], kVS[v], kOA[oa], kPN[p], c + 2), ";M_{ee} [GeV/c^{2}];pairs", kNb, edges());
            }
            for (int c = 0; c < kNC; ++c) hcB[t][v][oa][c] = new TH1D(Form("hcB_%s_%s_oa%g_c%d", kTN[t], kVS[v], kOA[oa], c + 2), ";M_{ee} [GeV/c^{2}];true background pairs", kNb, edges());
        }
    }
    for (int r = 0; r < 2; ++r) for (int c = 0; c < kNCat; ++c) hchi[r][c] = new TH1D(Form("hchi_r%d_cat%d", r, c), ";log_{10}#chi^{2}_{time};pairs", 120, -6., 6.);

    double cosOA[kNOA]; for (int i = 0; i < kNOA; ++i) cosOA[i] = std::cos(kOA[i] * TMath::DegToRad());

    auto fillRows = [&](int p, bool pt3, std::vector<Row>& g, int cl) {
        const int ic = nearestClass(cl);
        for (auto& x : g) {
            double cs, m; pairKin(x.a, x.b, cs, m);
            const int cat = category(x.a, x.b);
            for (int t = 0; t < kNT; ++t) { if (t == 1 && !pt3) continue;
                for (int oa = 0; oa < kNOA; ++oa) { if (oa > 0 && cs >= cosOA[oa]) continue;
                    for (int v = 0; v < kNV; ++v) { if (v == 1 && x.best < 0.5f) continue;
                        h[t][v][oa][p][cat]->Fill(m); hc[t][v][oa][p][ic]->Fill(m);
                        if (p == 0 && isB(cat)) hcB[t][v][oa][ic]->Fill(m);
                    }
                }
            }
            if (p == 0 && cs < cosOA[1]) {
                const int r = m < 0.14 ? 0 : (m < 0.80 ? 1 : -1);
                if (r >= 0) hchi[r][cat]->Fill(x.chi2 > 0 ? std::log10(x.chi2) : -5.99);
            }
        }
    };

    struct Br { Float_t ev, vx, vy, vz, best, trig, chi2; Leg A, B; } b;
    auto attach = [&](TTree* t, const std::string& p1, const std::string& p2, bool truth) {
        t->SetBranchStatus("*", 0);
        auto sb = [&](const std::string& n, Float_t* a) { t->SetBranchStatus(n.c_str(), 1); t->SetBranchAddress(n.c_str(), a); };
        sb("event", &b.ev); sb("eVertReco_x", &b.vx); sb("eVertReco_y", &b.vy); sb("eVertReco_z", &b.vz); sb("trigbit", &b.trig);
        auto leg = [&](const std::string& pr, Leg& L) { sb(pr + "p", &L.p); sb(pr + "theta", &L.th); sb(pr + "phi", &L.ph);
            if (truth) { sb(pr + "sim_id", &L.id); sb(pr + "sim_parentid", &L.par); sb(pr + "sim_grandparentid", &L.gp); sb(pr + "sim_processid", &L.proc);
                sb(pr + "sim_vertexx", &L.vx); sb(pr + "sim_vertexy", &L.vy); sb(pr + "sim_vertexz", &L.vz);
                sb(pr + "rich_sim_geanttrack", &L.r1); sb(pr + "rich_sim_geanttrack2", &L.r2); sb(pr + "rich_sim_geanttrack3", &L.r3); } };
        leg(p1, b.A); leg(p2, b.B);
        if (truth) { sb("isBest", &b.best); sb(p1 + "hyp_best_chi2_red", &b.chi2); }
    };

    // ---- (a) like-sign events into memory ----------------------------------
    std::vector<LSEvent> ls[2];
    std::unordered_set<unsigned long long> lsKeys;
    const char* lsTree[2] = {"EpEp_ID", "EmEm_ID"}; const char* lsP1[2] = {"ep1_", "em1_"}; const char* lsP2[2] = {"ep2_", "em2_"};
    long nRows[3] = {0, 0, 0}, nEvents[3] = {0, 0, 0}, nNoBest[3] = {0, 0, 0};
    for (int k = 0; k < 2; ++k) {
        TTree* t = (TTree*)f->Get(lsTree[k]); attach(t, lsP1[k], lsP2[k], true);
        LSEvent cur; bool have = false; Float_t ce = -1, cx = 0, cy = 0, cz = 0;
        auto push = [&]() { if (have && !cur.rows.empty()) { lsKeys.insert(cur.key); ls[k].push_back(std::move(cur)); } cur = LSEvent(); };
        for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            t->GetEntry(i);
            if (!have || b.ev != ce || b.vx != cx || b.vy != cy || b.vz != cz) { push(); ce = b.ev; cx = b.vx; cy = b.vy; cz = b.vz; have = true; cur.key = evKey(b.ev, b.vx, b.vy, b.vz); cur.pt3 = (int(b.trig) & 8192) != 0; }
            if (b.vz <= vertex_z_min) continue;
            Row r; r.a = b.A; r.b = b.B; r.best = b.best; r.chi2 = b.chi2; derive(r.a); derive(r.b); cur.rows.push_back(r); ++nRows[k + 1];
        }
        push();
        nEvents[k + 1] = (long)ls[k].size();
    }
    printf("[cb_study] like-sign events read: e+e+ %zu, e-e- %zu (%.0f s)\n", ls[0].size(), ls[1].size(), sw.RealTime()); sw.Continue();

    // ---- (b) e+e- pass: pairs, classes, lookup for like-sign events ---------
    std::unordered_map<unsigned long long, std::pair<int, int> > lookup; lookup.reserve(lsKeys.size() * 2);
    {
        TTree* t = (TTree*)f->Get("EpEm_ID"); attach(t, "ep_", "em_", true);
        std::vector<Row> g; Float_t ce = -1, cx = 0, cy = 0, cz = 0; bool have = false, pt3 = false;
        auto flush = [&]() {
            if (g.empty()) return;
            std::vector<const Leg*> P, M;
            for (auto& x : g) {
                bool fp = false, fm = false;
                for (auto* l : P) if (sameLeg(*l, x.a)) fp = true;
                for (auto* l : M) if (sameLeg(*l, x.b)) fm = true;
                if (!fp) P.push_back(&x.a);
                if (!fm) M.push_back(&x.b);
            }
            const int np = (int)P.size(), nm = (int)M.size();
            const unsigned long long key = evKey(ce, cx, cy, cz);
            if (lsKeys.count(key)) lookup[key] = std::make_pair(np, nm);
            int nb = 0; for (auto& x : g) nb += (x.best > 0.5f);
            if (nb != 1) ++nNoBest[0];
            ++nEvents[0];
            fillRows(0, pt3, g, np + nm);
            g.clear();
        };
        for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            t->GetEntry(i);
            if (!have || b.ev != ce || b.vx != cx || b.vy != cy || b.vz != cz) { flush(); ce = b.ev; cx = b.vx; cy = b.vy; cz = b.vz; have = true; pt3 = (int(b.trig) & 8192) != 0; }
            if (b.vz <= vertex_z_min) continue;
            Row r; r.a = b.A; r.b = b.B; r.best = b.best; r.chi2 = b.chi2; derive(r.a); derive(r.b); g.push_back(r); ++nRows[0];
            if (i % 2000000 == 0 && i) printf("   e+e- rows %lld (%.0f s)\n", i, sw.RealTime()), sw.Continue();
        }
        flush();
    }
    printf("[cb_study] e+e- done (%.0f s)\n", sw.RealTime()); sw.Continue();

    // ---- (c) like-sign pairs ---------------------------------------------
    std::vector<MEv> lsOnly[2];   // [all, PT3 flag kept via separate lists]
    std::vector<char> lsOnlyPT3;
    std::vector<MEv> lsOnlyAll;
    long nLsMatched[2] = {0, 0};
    for (int k = 0; k < 2; ++k) for (auto& e : ls[k]) {
        int nb = 0; for (auto& x : e.rows) nb += (x.best > 0.5f);
        if (nb != 1) ++nNoBest[k + 1];
        std::vector<const Leg*> S;
        for (auto& x : e.rows) { bool fa = false, fb = false; for (auto* l : S) { if (sameLeg(*l, x.a)) fa = true; if (sameLeg(*l, x.b)) fb = true; } if (!fa) S.push_back(&x.a); if (!fb && !sameLeg(x.a, x.b)) S.push_back(&x.b); }
        int cl = (int)S.size();
        auto it = lookup.find(e.key);
        if (it != lookup.end()) { cl = it->second.first + it->second.second; ++nLsMatched[k]; }
        else {
            MEv me; for (auto* l : S) { MLeg ml{l->ux, l->uy, l->uz, l->e, l->p}; (k == 0 ? me.ep : me.em).push_back(ml); }
            lsOnlyAll.push_back(std::move(me)); lsOnlyPT3.push_back(e.pt3 ? 1 : 0);
        }
        fillRows(k + 1, e.pt3, e.rows, cl);
    }
    printf("[cb_study] like-sign filled; events found in EpEm_ID: e+e+ %ld, e-e- %ld; like-sign-only events %zu (%.0f s)\n", nLsMatched[0], nLsMatched[1], lsOnlyAll.size(), sw.RealTime()); sw.Continue();

    // ---- (d) event mixing on the simulation --------------------------------
    {
        TTree* t = (TTree*)f->Get("EpEm_ID"); attach(t, "ep_", "em_", false);
        const Long64_t N = t->GetEntries();
        std::vector<MEv> pool[kNT]; size_t head[kNT] = {0, 0}; for (int tt = 0; tt < kNT; ++tt) pool[tt].resize(pool_depth);
        size_t filled[kNT] = {0, 0};
        auto addPair = [&](int tt, int p, const MLeg& a, const MLeg& b) {
            const double cs = a.ux * b.ux + a.uy * b.uy + a.uz * b.uz;
            const double m2 = 2 * kMe * kMe + 2 * (a.e * b.e - a.p * b.p * cs);
            const int bin = massBin(m2 > 0 ? std::sqrt(m2) / 1000. : 0.);
            if (bin <= 0) return;
            for (int oa = 0; oa < kNOA; ++oa) { if (oa > 0 && cs >= cosOA[oa]) continue; hmix[tt][oa][p]->AddBinContent(bin); }
        };
        auto mixEvent = [&](int tt, const MEv& E) {
            for (size_t ib = 0; ib < filled[tt]; ++ib) {
                const MEv& B = pool[tt][ib];
                for (auto& a : E.ep) for (auto& c : B.em) addPair(tt, 0, a, c);
                for (auto& a : E.em) for (auto& c : B.ep) addPair(tt, 0, a, c);
                for (auto& a : E.ep) for (auto& c : B.ep) addPair(tt, 1, a, c);
                for (auto& a : E.em) for (auto& c : B.em) addPair(tt, 2, a, c);
            }
            pool[tt][head[tt]] = E; head[tt] = (head[tt] + 1) % pool_depth; if (filled[tt] < (size_t)pool_depth) ++filled[tt];
        };
        MEv cur; Float_t ce = -1, cx = 0, cy = 0, cz = 0; bool have = false, ok = false, pt3 = false;
        std::vector<Leg> P, M; size_t lsIdx = 0; const size_t nLs = lsOnlyAll.size(); long iEv = 0; const double evTot = std::max<double>(1., nEvents[0]);
        auto flush = [&]() {
            if (!have || !ok || P.empty()) { P.clear(); M.clear(); return; }
            cur.ep.clear(); cur.em.clear();
            for (auto& l : P) cur.ep.push_back({l.ux, l.uy, l.uz, l.e, double(l.p)});
            for (auto& l : M) cur.em.push_back({l.ux, l.uy, l.uz, l.e, double(l.p)});
            mixEvent(0, cur); if (pt3) mixEvent(1, cur);
            ++iEv;
            while (lsIdx < nLs && double(lsIdx) / nLs <= iEv / evTot) { mixEvent(0, lsOnlyAll[lsIdx]); if (lsOnlyPT3[lsIdx]) mixEvent(1, lsOnlyAll[lsIdx]); ++lsIdx; }
            P.clear(); M.clear();
        };
        for (Long64_t i = 0; i < N; ++i) {
            t->GetEntry(i);
            if (!have || b.ev != ce || b.vx != cx || b.vy != cy || b.vz != cz) { flush(); ce = b.ev; cx = b.vx; cy = b.vy; cz = b.vz; have = true; pt3 = (int(b.trig) & 8192) != 0; ok = b.vz > vertex_z_min; }
            if (!ok) continue;
            Leg a = b.A, c = b.B; derive(a); derive(c);
            bool fa = false, fc = false; for (auto& l : P) if (sameLeg(l, a)) fa = true; for (auto& l : M) if (sameLeg(l, c)) fc = true;
            if (!fa) P.push_back(a); if (!fc) M.push_back(c);
        }
        flush();
        printf("[cb_study] mixing done, like-sign-only events injected %zu/%zu (%.0f s)\n", lsIdx, nLs, sw.RealTime()); sw.Continue();
    }

    // ---- derived quantities ------------------------------------------------
    auto zero = [&](const char* n) { TH1D* x = new TH1D(n, ";M_{ee} [GeV/c^{2}];pairs", kNb, edges()); x->SetDirectory(fo); return x; };
    TH1D* Npm[kNT][kNV][kNOA], *Npp[kNT][kNV][kNOA], *Nmm[kNT][kNV][kNOA], *Btrue[kNT][kNV][kNOA], *Blep[kNT][kNV][kNOA], *Bfake[kNT][kNV][kNOA], *Sgen[kNT][kNV][kNOA], *Sdal[kNT][kNV][kNOA], *Sconv[kNT][kNV][kNOA], *CBhat[kNT][kNV][kNOA], *Same[kNT][kNV][kNOA];
    TH1D* NppL[kNT][kNV][kNOA], *NmmL[kNT][kNV][kNOA];
    for (int t = 0; t < kNT; ++t) for (int v = 0; v < kNV; ++v) for (int oa = 0; oa < kNOA; ++oa) {
        TString s = Form("%s_%s_oa%g", kTN[t], kVS[v], kOA[oa]);
        Npm[t][v][oa] = zero("Npm_" + s); Npp[t][v][oa] = zero("Npp_" + s); Nmm[t][v][oa] = zero("Nmm_" + s);
        NppL[t][v][oa] = zero("NppLep_" + s); NmmL[t][v][oa] = zero("NmmLep_" + s);
        Btrue[t][v][oa] = zero("Btrue_" + s); Blep[t][v][oa] = zero("Blep_" + s); Bfake[t][v][oa] = zero("Bfake_" + s);
        Sgen[t][v][oa] = zero("Sgen_" + s); Sdal[t][v][oa] = zero("Sdalitz_" + s); Sconv[t][v][oa] = zero("Sconv_" + s); CBhat[t][v][oa] = zero("CBhat_" + s); Same[t][v][oa] = zero("Same_" + s);
        for (int c = 0; c < kNCat; ++c) {
            Npm[t][v][oa]->Add(h[t][v][oa][0][c]); Npp[t][v][oa]->Add(h[t][v][oa][1][c]); Nmm[t][v][oa]->Add(h[t][v][oa][2][c]);
            if (isB(c)) Btrue[t][v][oa]->Add(h[t][v][oa][0][c]);
            if (isCBLep(c)) { Blep[t][v][oa]->Add(h[t][v][oa][0][c]); NppL[t][v][oa]->Add(h[t][v][oa][1][c]); NmmL[t][v][oa]->Add(h[t][v][oa][2][c]); }
            if (isFake(c)) Bfake[t][v][oa]->Add(h[t][v][oa][0][c]);
            if (isGen(c)) Sgen[t][v][oa]->Add(h[t][v][oa][0][c]);
            if (isDalitz(c)) Sdal[t][v][oa]->Add(h[t][v][oa][0][c]);
            if (isConv(c)) Sconv[t][v][oa]->Add(h[t][v][oa][0][c]);
            if (c == kSamePart) Same[t][v][oa]->Add(h[t][v][oa][0][c]);
        }
        for (int bn = 1; bn <= kNb; ++bn) { const double a = Npp[t][v][oa]->GetBinContent(bn), m = Nmm[t][v][oa]->GetBinContent(bn);
            if (a > 0 && m > 0) { const double cb = 2 * std::sqrt(a * m); CBhat[t][v][oa]->SetBinContent(bn, cb); CBhat[t][v][oa]->SetBinError(bn, 0.5 * cb * std::sqrt(1 / a + 1 / m)); } }
    }

    const int kNR = 9; const double R[kNR][2] = {{0, 0.02}, {0.02, 0.05}, {0.05, 0.14}, {0.14, 0.30}, {0.30, 0.56}, {0.56, 0.80}, {0.80, 1.40}, {0, 0.14}, {0.14, 0.80}};
    auto b1 = [&](double lo) { return massBin(lo + 1e-6); };
    auto b2 = [&](double hi) { return massBin(hi - 1e-6); };
    auto I = [&](TH1D* x, double lo, double hi) { double s = 0; for (int bn = b1(lo); bn <= b2(hi); ++bn) s += x->GetBinContent(bn); return s; };

    for (int t = 0; t < kNT; ++t) {
        printf("\n######## SET: %s events ########\n", kTN[t]);
        printf("rows after vertex cut: e+e- %ld, e+e+ %ld, e-e- %ld | events: %ld, %ld, %ld | events without isBest==1: %ld %ld %ld\n", nRows[0], nRows[1], nRows[2], nEvents[0], nEvents[1], nEvents[2], nNoBest[0], nNoBest[1], nNoBest[2]);
        printf("\n(1) Selection scan: M range | OA | ranking | N+- | S_true (Dalitz, conversion) | B_true (leptons, fake) | CB_hat | CB_hat/B_true | S_hat | S_hat/S_true | S_true/B_true | S_hat/sigma\n");
        for (int r = 0; r < kNR; ++r) for (int oa = 0; oa < kNOA; ++oa) for (int v = 0; v < kNV; ++v) {
            const double lo = R[r][0], hi = R[r][1];
            const double N = I(Npm[t][v][oa], lo, hi), S = I(Sgen[t][v][oa], lo, hi), Sd = I(Sdal[t][v][oa], lo, hi), Sc = I(Sconv[t][v][oa], lo, hi), B = I(Btrue[t][v][oa], lo, hi), Bl = I(Blep[t][v][oa], lo, hi), Bf = I(Bfake[t][v][oa], lo, hi), X = I(Same[t][v][oa], lo, hi);
            double cb = 0, e2 = 0; for (int bn = b1(lo); bn <= b2(hi); ++bn) { cb += CBhat[t][v][oa]->GetBinContent(bn); e2 += std::pow(CBhat[t][v][oa]->GetBinError(bn), 2); }
            const double Sh = N - cb, sig = std::sqrt(N + e2);
            const double rat = B > 0 ? cb / B : 0., erat = (B > 0 && cb > 0) ? rat * std::sqrt(e2 / (cb * cb) + 1. / B) : 0.;
            printf("  %.2f-%.2f | %2g | %-9s | %9.0f | %9.0f (%9.0f, %8.0f) | %8.0f (%8.0f, %7.0f) same %5.0f | %8.0f | %.3f+-%.3f | %9.0f | %.4f | %8.2f | %7.1f\n", lo, hi, kOA[oa], kVN[v], N, S, Sd, Sc, B, Bl, Bf, X, cb, rat, erat, Sh, S > 0 ? Sh / S : 0., B > 0 ? S / B : 0., sig > 0 ? Sh / sig : 0.);
        }
        printf("\n(2) Composition by origin (pairs): for each pattern, OA and ranking, M ranges <0.14 | 0.14-0.80\n");
        for (int p = 0; p < 3; ++p) for (int oa : {0, 1}) {
            printf(" -- %s, OA > %g deg: category | M<0.14: isBest>=0, isBest==1, kept | 0.14-0.80: isBest>=0, isBest==1, kept\n", kPT[p], kOA[oa]);
            for (int c = 0; c < kNCat; ++c) { const double a0 = I(h[t][0][oa][p][c], 0, 0.14), b0 = I(h[t][1][oa][p][c], 0, 0.14), a1 = I(h[t][0][oa][p][c], 0.14, 0.80), b1v = I(h[t][1][oa][p][c], 0.14, 0.80);
                if (a0 + a1 < 1) continue;
                printf("   %-44s | %9.0f %9.0f %.3f | %9.0f %9.0f %.3f\n", kCatName[c], a0, b0, a0 > 0 ? b0 / a0 : 0., a1, b1v, a1 > 0 ? b1v / a1 : 0.); }
        }
        printf("\n(3) k factors, all combinations: M range | OA | k_true(all) | k_true(true leptons) | k_mix,sim | r(all) | r(leptons) | closure B_true vs k_mix*CB_hat\n");
        for (int r = 0; r < kNR; ++r) for (int oa = 0; oa < kNOA; ++oa) {
            const double lo = R[r][0], hi = R[r][1];
            double cb = 0, cbL = 0, mixN = 0, mixF = 0, kcb = 0, e2 = 0, e2L = 0;
            for (int bn = b1(lo); bn <= b2(hi); ++bn) {
                const double c0 = CBhat[t][0][oa]->GetBinContent(bn); cb += c0; e2 += std::pow(CBhat[t][0][oa]->GetBinError(bn), 2);
                const double a = NppL[t][0][oa]->GetBinContent(bn), m = NmmL[t][0][oa]->GetBinContent(bn); if (a > 0 && m > 0) { const double g = 2 * std::sqrt(a * m); cbL += g; e2L += 0.25 * g * g * (1 / a + 1 / m); }
                const double pm = hmix[t][oa][0]->GetBinContent(bn), pp = hmix[t][oa][1]->GetBinContent(bn), mm = hmix[t][oa][2]->GetBinContent(bn);
                if (pp > 0 && mm > 0) { const double g = 2 * std::sqrt(pp * mm); mixN += pm; mixF += g; kcb += (pm / g) * c0; } else kcb += c0;
            }
            const double B = I(Btrue[t][0][oa], lo, hi), Bl = I(Blep[t][0][oa], lo, hi);
            const double kt = cb > 0 ? B / cb : 0., ekt = (B > 0 && cb > 0) ? kt * std::sqrt(1 / B + e2 / (cb * cb)) : 0.;
            const double ktL = cbL > 0 ? Bl / cbL : 0., ektL = (Bl > 0 && cbL > 0) ? ktL * std::sqrt(1 / Bl + e2L / (cbL * cbL)) : 0.;
            const double km = mixF > 0 ? mixN / mixF : 0.;
            printf("  %.2f-%.2f | %2g | %.3f+-%.3f | %.3f+-%.3f | %.4f | %.3f | %.3f | %.3f\n", lo, hi, kOA[oa], kt, ekt, ktL, ektL, km, km > 0 ? kt / km : 0., km > 0 ? ktL / km : 0., B > 0 ? kcb / B : 0.);
        }
        printf("\n(4) Best-pair background: true vs naive vs class-weighted with truth eps_B,c (per bin): M range | OA | B_true(best) | naive CB_hat(best) | ratio | sum_c eps^truth_B,c CB_c(all pairs) | ratio | all-pairs per-class sum_c CB_c / B_true(all)\n");
        for (int r = 7; r < kNR; ++r) for (int oa : {0, 1}) {
            const double lo = R[r][0], hi = R[r][1];
            double naive = 0, corr = 0, cls = 0;
            for (int bn = b1(lo); bn <= b2(hi); ++bn) { naive += CBhat[t][1][oa]->GetBinContent(bn);
                for (int c = 0; c < kNC; ++c) { const double a = hc[t][0][oa][1][c]->GetBinContent(bn), m = hc[t][0][oa][2][c]->GetBinContent(bn); if (a <= 0 || m <= 0) continue; const double g = 2 * std::sqrt(a * m); cls += g;
                    const double ba = hcB[t][0][oa][c]->GetBinContent(bn), bb = hcB[t][1][oa][c]->GetBinContent(bn);
                    double eps = ba > 0 ? bb / ba : (c == 0 ? 1. : hcB[t][1][oa][c]->Integral() / std::max(1., hcB[t][0][oa][c]->Integral())); corr += eps * g; } }
            const double Bb = I(Btrue[t][1][oa], lo, hi), Ba = I(Btrue[t][0][oa], lo, hi);
            printf("  %.2f-%.2f | %g | %.0f | %.0f | %.3f | %.0f | %.3f | %.3f\n", lo, hi, kOA[oa], Bb, naive, Bb > 0 ? naive / Bb : 0., corr, Bb > 0 ? corr / Bb : 0., Ba > 0 ? cls / Ba : 0.);
        }
        printf("\n(5) Multiplicity classes (true background and genuine pairs, OA>0): class | M range | eps_B,c (truth) | eps_S,c (truth)\n");
        for (int c = 0; c < kNC; ++c) for (int r : {7, 8}) { const double lo = R[r][0], hi = R[r][1];
            double ba = I(hcB[t][0][0][c], lo, hi), bb = I(hcB[t][1][0][c], lo, hi), na = I(hc[t][0][0][0][c], lo, hi), nb = I(hc[t][1][0][0][c], lo, hi);
            printf("  n=%d%s | %.2f-%.2f | B: %.0f -> %.0f (%.4f) | genuine+same: %.0f -> %.0f (%.4f)\n", c + 2, c == 3 ? "+" : " ", lo, hi, ba, bb, ba > 0 ? bb / ba : 0., na - ba, nb - bb, (na - ba) > 0 ? (nb - bb) / (na - ba) : 0.); }
    }
    printf("\n(6) log10(time chi2) categories (all events, isBest>=0, OA>2): fraction of pairs with chi2 < 1, < 4, < 9\n");
    for (int r = 0; r < 2; ++r) for (int c = 0; c < kNCat; ++c) { if (hchi[r][c]->GetEntries() < 50) continue; const double tot = hchi[r][c]->Integral(0, 121);
        auto fr = [&](double th) { return hchi[r][c]->Integral(0, hchi[r][c]->FindBin(std::log10(th) - 1e-9)) / tot; };
        printf("  %s %-44s n=%8.0f  %.3f %.3f %.3f\n", r == 0 ? "M<0.14   " : "0.14-0.80", kCatName[c], tot, fr(1), fr(4), fr(9)); }

    fo->Write();
    fo->Close();
    printf("\nwritten %s/cb_study_sim_%s.root (%.0f s)\n", out_dir, tag, sw.RealTime());
}
