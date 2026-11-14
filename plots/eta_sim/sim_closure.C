// ============================================================================
// sim_closure.C  --  weights of the SMASH simulation rechecked on the full production
// ============================================================================
// HADES p+p, T = 4.53 GeV (feb22), pi+ pi- e+ e-. Reads the six FAT outputs
//   /hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/output_{std,lep}_{pippimepem,pippimepep,pippimemem}.root
// (tree scan_nt, every row, the generator codes of the four slots included) and
// answers five questions:
//   1. lep/std closure of the rows the lepton set provides in mode A: pi0 Dalitz
//      pairs, eta Dalitz pairs, two different Dalitz decays (w1 w2), conversions
//      of Dalitz photons, and all rows (std rows dropped in mode A against lep
//      rows owned in mode A); pions from the common SMASH final state (no lepton
//      and no product of a rho0/omega/phi/eta decay in a pion slot); after "pre"
//      and "pre+OA>2"; totals and spectra of M(ee), OA and M4
//   2. partition mode B against mode A, and w_db against w, per M4 bin and per
//      window, for every template (physical background in particular); the
//      difference B - A split by lepton-pair class
//   3. the open rules-v1 gaps marked in sim_flags (bits 2, 12, 13, 14) and the
//      missing natural-decay factors of w_db: rows and weight per window and
//      template, and the change of the template if the gap were closed
//   4. ghost-lepton rows in the signal classes: M4 against the clean rows,
//      peak position and width, window cross sections
//   5. the float-code dictionary over every stored code value of the outputs
//
// Selection "pre" = isBest >= 0, vz > -500, vz > -160 mm, cor_ok == 1, PT3
// ((trigbit & 8192) != 0), the same as simeng::PassPre; "pre+OA>2" adds oa > 2 deg.
// Statistics: per bin sum of weights and variance sum over event groups (adjacent
// rows with equal (event, eVertReco_z) in one tree) of the squared group sum; all
// sums in double precision (TVectorD, TH1D, double accumulators).
//
// Stages (run from /home/przygoda/HADES/FAT/FAT, build directory outside the repository):
//   stage 1, one tree k (0 std_us, 1 std_pp, 2 std_mm, 3 lep_us, 4 lep_pp, 5 lep_mm; -1 = all):
//     root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/aclic_closure", kTRUE)' 'plots/eta_sim/sim_closure.C+(1,3)'
//     writes /hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/closure_pass_<k>.root
//   stage 2 (tables and figures from the six pass files):
//     root -l -b -q -e 'gSystem->SetBuildDir(".../scan/aclic_closure", kTRUE)' 'plots/eta_sim/sim_closure.C+(2)'
//   stage 0: stage 1 for all trees, then stage 2
// Tables: plots/eta_sim/tables/closure_*.txt; figures: plots/eta_sim/figs/closure_*.{png,pdf}
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "TBox.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TObjString.h"
#include "TPad.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TVectorD.h"

#include "SimEngine.h"

namespace scl {

using simtruth::Trk;

const char* const kPassFmt = "/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/closure_pass_%d.root";
const char* const kInName[6] = {"output_std_pippimepem.root", "output_std_pippimepep.root", "output_std_pippimemem.root",
                                "output_lep_pippimepem.root", "output_lep_pippimepep.root", "output_lep_pippimemem.root"};

// ---------------------------------------------------------------------------
// constants of the closure weights (PDG values of src/sim_truth.h)
// ---------------------------------------------------------------------------
// conversions of Dalitz photons, standard set: every photon of a primary pi0 (eta)
// counted, scaled by the share of Dalitz photons under PDG branching ratios:
// per pi0 HGeant gives 2 (1 - 0.01198) + 0.01198 photons, per eta 2 x 0.3756 +
// 0.0604 + 0.0478 = 0.8593 photons directly from the eta (rules v2, HGeant ETA card).
const double kNgamPi0HG = 2.0 * (1.0 - simtruth::kHG_pi0D) + simtruth::kHG_pi0D;     // 1.98802
const double kNgamEtaHG = 2.0 * simtruth::kHG_etaGG + simtruth::kHG_etaD + simtruth::kHG_etaPPG;  // 0.8593
const double kFconvPi0 = simtruth::kPDG_pi0D / kNgamPi0HG;   // std weight factor, pi0 photon conversions
const double kFconvEta = simtruth::kPDG_etaD / kNgamEtaHG;   // std weight factor, eta photon conversions
// windows in M4 [GeV/c2] (as the engine's data/sim comparison)
const int kNWin = 7;
const double kWinEdge[kNWin + 1] = {0.30, 0.455, 0.5225, 0.5336, 0.5542, 0.58, 0.70, 1.10};
const char* const kWinName[kNWin] = {"below bump", "bump", "bump-peak gap", "peak", "peak-0.58", "above peak", "high mass"};

// ---------------------------------------------------------------------------
// binnings
// ---------------------------------------------------------------------------
struct Bins {
    std::vector<double> e;
    bool all = false;   // one bin holding every row (any x, NaN included)
    int n() const { return all ? 1 : (int)e.size() - 1; }
    int find(double x) const {
        if (all) return 0;
        if (!(x >= e.front() && x < e.back())) return -1;
        return (int)(std::upper_bound(e.begin(), e.end(), x) - e.begin()) - 1;
    }
    double lo(int i) const { return all ? 0 : e[i]; }
    double hi(int i) const { return all ? 1 : e[i + 1]; }
};
Bins Uniform(int n, double lo, double hi) {
    Bins b;
    for (int i = 0; i <= n; ++i) b.e.push_back(lo + (hi - lo) * i / n);
    return b;
}
Bins Edges(std::vector<double> e) { Bins b; b.e = e; return b; }

Bins bTot, bMee, bOA, bM4c, bM4, bWin, bM4f, bRes;
void InitBins() {
    bTot.all = true;
    bMee = Edges({0, 0.005, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13,
                  0.14, 0.16, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.60});
    bOA = Edges({0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10, 12, 15, 20, 25, 30, 40, 50, 60, 80, 100, 180});
    {
        std::vector<double> e;
        for (int i = 0; i <= 40; ++i) e.push_back(0.25 + 0.025 * i);
        e.push_back(1.5); e.push_back(2.5);
        bM4c = Edges(e);
    }
    bM4 = Uniform(80, 0.30, 1.10);
    bWin = Edges(std::vector<double>(kWinEdge, kWinEdge + kNWin + 1));
    bM4f = Uniform(110, 0.40, 0.62);
    bRes = Uniform(100, -0.06, 0.04);
}

// ---------------------------------------------------------------------------
// accumulators with event-group variance
// ---------------------------------------------------------------------------
// GAcc: s[i] sum w, v[i] sum over groups of (group sum in bin i)^2, n[i] rows;
// index nb (last) = every filled row (any x). PAcc: two weights a, b per row
// (the same rows under two variants) with the group sums of a^2, b^2, a b, so
// that the ratio B/A gets the error of a ratio estimator with group correlation:
// V(R) = sum_g (b_g - R a_g)^2 / A^2.
struct GAcc {
    const Bins* b = nullptr;
    std::vector<double> s, v, n;
    std::vector<std::pair<int, double>> buf;
    bool dirty = false;
};
struct PAcc {
    const Bins* b = nullptr;
    std::vector<double> a, bb, aa, b2, ab, n;
    std::vector<std::tuple<int, double, double>> buf;
    bool dirty = false;
};
std::map<std::string, GAcc*> gReg;
std::map<std::string, PAcc*> pReg;
std::vector<GAcc*> gDirty;
std::vector<PAcc*> pDirty;

GAcc* G(const std::string& name, const Bins& b) {
    auto it = gReg.find(name);
    if (it != gReg.end()) return it->second;
    GAcc* a = new GAcc;
    a->b = &b;
    const int nb = b.n() + 1;
    a->s.assign(nb, 0.); a->v.assign(nb, 0.); a->n.assign(nb, 0.);
    gReg[name] = a;
    return a;
}
PAcc* P(const std::string& name, const Bins& b) {
    auto it = pReg.find(name);
    if (it != pReg.end()) return it->second;
    PAcc* a = new PAcc;
    a->b = &b;
    const int nb = b.n() + 1;
    a->a.assign(nb, 0.); a->bb.assign(nb, 0.); a->aa.assign(nb, 0.); a->b2.assign(nb, 0.); a->ab.assign(nb, 0.);
    a->n.assign(nb, 0.);
    pReg[name] = a;
    return a;
}
inline void Fill(GAcc* g, double x, double w) {
    const int nb = g->b->n();
    const int ib = g->b->find(x);
    g->s[nb] += w; g->n[nb] += 1;
    if (ib >= 0) { g->s[ib] += w; g->n[ib] += 1; }
    g->buf.emplace_back(ib, w);
    if (!g->dirty) { g->dirty = true; gDirty.push_back(g); }
}
inline void Fill2(PAcc* p, double x, double wa, double wb) {
    const int nb = p->b->n();
    const int ib = p->b->find(x);
    p->a[nb] += wa; p->bb[nb] += wb; p->n[nb] += 1;
    if (ib >= 0) { p->a[ib] += wa; p->bb[ib] += wb; p->n[ib] += 1; }
    p->buf.emplace_back(ib, wa, wb);
    if (!p->dirty) { p->dirty = true; pDirty.push_back(p); }
}
void FlushGroup() {
    for (GAcc* g : gDirty) {
        const int nb = g->b->n();
        auto& e = g->buf;
        std::sort(e.begin(), e.end());
        double all = 0;
        size_t i = 0;
        while (i < e.size()) {
            const int k = e[i].first;
            double acc = 0;
            while (i < e.size() && e[i].first == k) { acc += e[i].second; ++i; }
            all += acc;
            if (k >= 0) g->v[k] += acc * acc;
        }
        g->v[nb] += all * all;
        e.clear();
        g->dirty = false;
    }
    gDirty.clear();
    for (PAcc* p : pDirty) {
        const int nb = p->b->n();
        auto& e = p->buf;
        std::sort(e.begin(), e.end());
        double ta = 0, tb = 0;
        size_t i = 0;
        while (i < e.size()) {
            const int k = std::get<0>(e[i]);
            double a = 0, b = 0;
            while (i < e.size() && std::get<0>(e[i]) == k) { a += std::get<1>(e[i]); b += std::get<2>(e[i]); ++i; }
            ta += a; tb += b;
            if (k >= 0) { p->aa[k] += a * a; p->b2[k] += b * b; p->ab[k] += a * b; }
        }
        p->aa[nb] += ta * ta; p->b2[nb] += tb * tb; p->ab[nb] += ta * tb;
        e.clear();
        p->dirty = false;
    }
    pDirty.clear();
}

// counters (summed over the pass files)
std::map<std::string, double> ctr;
inline void C(const std::string& k, double v = 1.) { ctr[k] += v; }

// ---------------------------------------------------------------------------
// track helpers
// ---------------------------------------------------------------------------
inline long long Ic(float g) { return std::llround((double)g); }
// SMASH decay codes of mothers that the lepton set forces into dileptons, and eta
// products: their pions do not exist in the lepton set (the audit's missingCode)
inline bool MissingCode(float g) {
    const long long c = Ic(g);
    if (c == 52080908LL || c == 55080908LL || c == 52080907LL || c == 55080907LL) return true;
    if (c >= 100000 && c < 1000000) { const long long h = c / 10000; if (h == 41 || h == 52 || h == 55) return true; }
    if (c == 1717) return true;
    return false;
}
// pions from the common SMASH final state: both pion slots hold a hadron (or muon),
// no e+-, nothing from a rho0/omega/phi/eta decay, no product of a forced decay
inline bool CommonPions(const Trk t[4]) {
    for (int i = 2; i < 4; ++i) {
        const Trk& x = t[i];
        if (simtruth::isE(x)) return false;
        if (MissingCode(x.g0) || x.par == 17) return false;
        if (simtruth::lepForcedCode(x.g0)) return false;
    }
    return true;
}
inline bool IsPi0EtaForced(float g) { return g == 75101.f || g == 175101.f; }
// eta descendant of phi -> eta gamma (code 550117, not the phi's own photon or its conversion)
inline bool PhiEtaProduct(const Trk& x) {
    if (x.g0 != 550117.f) return false;
    if (x.par == -1) return false;                 // direct phi product
    if (x.par == 1 && x.gp == -1) return false;    // conversion of the phi's photon
    return true;
}

// ---------------------------------------------------------------------------
// closure groups
// ---------------------------------------------------------------------------
enum CGroup { kPi0D = 0, kEtaD, kDD, kConvPi0, kConvEta, kModeA, kNCG };
const char* const kCGName[kNCG] = {"pi0 Dalitz pair", "eta Dalitz pair", "two Dalitz decays", "conv. pi0 Dalitz photon",
                                   "conv. eta Dalitz photon", "all rows (std dropped / lep owned, mode A)"};
const char* const kCGKey[kNCG] = {"pi0D", "etaD", "DD", "convPi0", "convEta", "modeA"};
enum CObs { kOTot = 0, kOMee, kOOA, kOM4, kOTotNoGhost, kNCO };

// ---------------------------------------------------------------------------
// input row
// ---------------------------------------------------------------------------
#define SCL_SLOT_COLS(X, T)                                                                         \
    X(T##_sim_id) X(T##_sim_parentid) X(T##_sim_grandparentid) X(T##_sim_processid)                 \
    X(T##_sim_geninfo) X(T##_sim_geninfo1) X(T##_sim_geninfo2) X(T##_sim_genweight)                 \
    X(T##_sim_vertexx) X(T##_sim_vertexy) X(T##_sim_vertexz) X(T##_sim_px) X(T##_sim_py) X(T##_sim_pz) \
    X(T##_sim_isghost)
#define SCL_COLS(X)                                                                                 \
    X(event) X(eVertReco_z) X(isBest) X(trigbit) X(cor_ok) X(oa) X(m_pippimepem_cor) X(m_ee_cor)    \
    X(sim_set) X(sim_cls) X(sim_lcls) X(sim_pcls) X(sim_tmpl) X(sim_w) X(sim_w_db) X(sim_own_A)      \
    X(sim_own_B) X(sim_flags) X(m4_true)                                                           \
    SCL_SLOT_COLS(X, ep) SCL_SLOT_COLS(X, em) SCL_SLOT_COLS(X, pip) SCL_SLOT_COLS(X, pim)
struct In {
#define X(n) Float_t n;
    SCL_COLS(X)
#undef X
};
#define SCL_TRK(T)                                                                                  \
    simtruth::Trk{(int)std::lround(in.T##_sim_id), (int)std::lround(in.T##_sim_parentid),           \
                  (int)std::lround(in.T##_sim_grandparentid), (int)std::lround(in.T##_sim_processid), \
                  in.T##_sim_geninfo, in.T##_sim_geninfo1, in.T##_sim_geninfo2, in.T##_sim_genweight, \
                  in.T##_sim_vertexx, in.T##_sim_vertexy, in.T##_sim_vertexz, in.T##_sim_px,         \
                  in.T##_sim_py, in.T##_sim_pz, (int)std::lround(in.T##_sim_isghost)}

// ---------------------------------------------------------------------------
// stage 1: one pass over one output
// ---------------------------------------------------------------------------
void Pass(int k) {
    TStopwatch sw;
    gReg.clear(); pReg.clear(); gDirty.clear(); pDirty.clear(); ctr.clear();
    const int st = k < 3 ? 0 : 1;            // 0 std, 1 lep
    const int h = k % 3;                     // 0 us, 1 pp, 2 mm
    const simtruth::Set set = st == 0 ? simtruth::kStd : simtruth::kLep;
    const std::string fn = std::string(simeng::kProdDir) + "/" + kInName[k];
    TFile* fin = TFile::Open(fn.c_str(), "READ");
    if (!fin || fin->IsZombie()) { printf("cannot open %s\n", fn.c_str()); return; }
    TTree* t = (TTree*)fin->Get("scan_nt");
    In in;
    t->SetBranchStatus("*", 0);
#define X(n) t->SetBranchStatus(#n, 1); t->SetBranchAddress(#n, &in.n);
    SCL_COLS(X)
#undef X
    t->SetCacheSize(400 * 1024 * 1024);

    // ---- accumulator pointers
    GAcc* cl[kNCG][2][kNCO];
    for (int g = 0; g < kNCG; ++g)
        for (int s = 0; s < 2; ++s)
            for (int o = 0; o < kNCO; ++o) {
                const Bins& b = o == kOMee ? bMee : o == kOOA ? bOA : o == kOM4 ? bM4c : bTot;
                cl[g][s][o] = G(Form("cl_%s_s%d_o%d_t%d_h%d", kCGKey[g], s, o, st, h), b);
            }
    // partition: tc 0..3 = T1..T4, 4 = T1+T2; vr 0 A(w) vs B(w), 1 A(w) vs A(w_db), 2 A(w) vs B(w_db),
    // 3 A(w) vs A(w_db x NF), 4 A(w) vs B(w_db x NF) (NF = natural-decay factors of sim_truth.h); bn 0 M4, 1 windows
    PAcc* pa[5][2][5][2];
    for (int tc = 0; tc < 5; ++tc)
        for (int s = 0; s < 2; ++s)
            for (int vr = 0; vr < 5; ++vr)
                for (int bn = 0; bn < 2; ++bn)
                    pa[tc][s][vr][bn] = P(Form("pa_tc%d_s%d_v%d_b%d_t%d_h%d", tc, s, vr, bn, st, h), bn ? bWin : bM4);
    // B - A of the physical background: part 0 std rows owned in B only, part 1 lep rows owned in A only; lcls group
    GAcc* pd[2][4][2];
    for (int s = 0; s < 2; ++s)
        for (int l = 0; l < 4; ++l)
            for (int bn = 0; bn < 2; ++bn)
                pd[s][l][bn] = G(Form("pd_s%d_l%d_b%d_t%d_h%d", s, l, bn, st, h), bn ? bWin : bM4);
    // flags: bit, template, ownership (0 owned A, 1 owned B only, 2 owned by neither) in windows
    GAcc* fl[2][15][4][3] = {};
    GAcc* tt[2][4][3];
    for (int s = 0; s < 2; ++s)
        for (int tp = 0; tp < 4; ++tp)
            for (int o = 0; o < 3; ++o) tt[s][tp][o] = G(Form("tt_s%d_tp%d_o%d_t%d_h%d", s, tp, o, st, h), bWin);
    // gap effects: k 0 bit 13 not owned, 1 bit 14 std not owned in A, 2 phi->eta gamma: W' - W (owned A),
    // 3 bit 2 (owned A), 4 natural factors in w_db: w_db (NF - 1) (owned A), 5 bit 12 (owned A)
    GAcc* ef[2][6][4][2];
    for (int s = 0; s < 2; ++s)
        for (int e = 0; e < 6; ++e)
            for (int tp = 0; tp < 4; ++tp)
                for (int bn = 0; bn < 2; ++bn)
                    ef[s][e][tp][bn] = G(Form("ef_s%d_e%d_tp%d_b%d_t%d_h%d", s, e, tp, bn, st, h), bn ? bWin : bM4);
    // w_db sums of the owned-A rows per template (denominator of the NF gap)
    GAcc* wdb[2][4][2];
    for (int s = 0; s < 2; ++s)
        for (int tp = 0; tp < 4; ++tp)
            for (int bn = 0; bn < 2; ++bn)
                wdb[s][tp][bn] = G(Form("wdb_s%d_tp%d_b%d_t%d_h%d", s, tp, bn, st, h), bn ? bWin : bM4);

    // ---- code dictionary: per column and slot type (0 lepton slots, 1 pion slots), stored value -> rows
    std::unordered_map<uint32_t, long long> codes[3][2];
    std::unordered_map<uint32_t, int> cand[3];   // cache: number of dictionary candidates
    auto ncand = [&](int col, float v) {
        uint32_t key; std::memcpy(&key, &v, 4);
        auto it = cand[col].find(key);
        if (it != cand[col].end()) return it->second;
        simtruth::Code c[8];
        const int n = simtruth::CodeCandidates(col, v, c, 8);
        cand[col][key] = n;
        return n;
    };

    // ---- signal rows (template 1) for the ghost study
    TFile* fout = TFile::Open(Form(kPassFmt, k), "RECREATE");
    TTree* sig = new TTree("sig", "template-1 rows after pre");
    Int_t o_cls, o_ghost, o_gl, o_oa2, o_grp, o_set, o_ownA;
    Float_t o_m4, o_m4t, o_mee, o_oa;
    Double_t o_w, o_wdb;
    sig->Branch("cls", &o_cls); sig->Branch("ghost", &o_ghost); sig->Branch("gl", &o_gl); sig->Branch("oa2", &o_oa2);
    sig->Branch("grp", &o_grp); sig->Branch("set", &o_set); sig->Branch("ownA", &o_ownA);
    sig->Branch("m4", &o_m4); sig->Branch("m4t", &o_m4t); sig->Branch("mee", &o_mee); sig->Branch("oa", &o_oa);
    sig->Branch("w", &o_w); sig->Branch("wdb", &o_wdb);
    sig->SetDirectory(fout);

    simtruth::Opt optPdg;   // pdg = true
    const Long64_t n = t->GetEntries();
    float keyE = -1, keyZ = -1e9;
    Long64_t grp = -1;
    for (Long64_t i = 0; i < n; ++i) {
        t->GetEntry(i);
        if (i == 0 || in.event != keyE || in.eVertReco_z != keyZ) {
            FlushGroup();
            ++grp;
            keyE = in.event; keyZ = in.eVertReco_z;
        }
        const Trk tr[4] = {SCL_TRK(ep), SCL_TRK(em), SCL_TRK(pip), SCL_TRK(pim)};
        const int flags = (int)std::lround(in.sim_flags);
        C("rows");
        // ---- code dictionary (every row)
        {
            bool unk = false, amb = false;
            for (int j = 0; j < 4; ++j) {
                const float gv[3] = {tr[j].g0, tr[j].g1, tr[j].g2};
                for (int c = 0; c < 3; ++c) {
                    uint32_t key; std::memcpy(&key, &gv[c], 4);
                    codes[c][j < 2 ? 0 : 1][key]++;
                    const int nc = ncand(c, gv[c]);
                    if (nc == 0) unk = true;
                    if (nc > 1) amb = true;
                }
            }
            if (unk) C("code_rows_unknown");
            if (amb) C("code_rows_ambiguous");
            if (unk != ((flags & simtruth::kFlagUnknownCode) != 0)) C("code_bit3_mismatch");
            if (amb != ((flags & simtruth::kFlagAmbiguousCode) != 0)) C("code_bit4_mismatch");
            for (int b = 0; b < 15; ++b) if (flags & (1 << b)) C(Form("allrows_bit%d", b));
        }
        const int trig = (int)std::lround(in.trigbit);
        const bool pre = in.isBest >= 0 && in.eVertReco_z > -500.f && in.eVertReco_z > -160.f &&
                         (int)std::lround(in.cor_ok) == 1 && (trig & 8192) != 0;
        if (!pre) continue;
        C("rows_pre");
        const bool oa2 = in.oa > 2.f;
        const int tmpl = (int)std::lround(in.sim_tmpl);
        const int lc = (int)std::lround(in.sim_lcls);
        const int cls = (int)std::lround(in.sim_cls);
        const bool ownA = std::lround(in.sim_own_A) == 1, ownB = std::lround(in.sim_own_B) == 1;
        const double w = in.sim_w, wdbv = in.sim_w_db;
        const double m4 = in.m_pippimepem_cor, mee = in.m_ee_cor, oa = in.oa;
        const bool ghostL = (flags & simtruth::kFlagGhostLepton) != 0;
        const bool clone = (flags & simtruth::kFlagSameParticle) != 0;
        const int tp = tmpl - 1;
        if (tp < 0 || tp > 3) { C("bad_tmpl"); continue; }

        // ---- closure group membership and weight
        double cw[kNCG];
        bool in_g[kNCG] = {false, false, false, false, false, false};
        const bool common = CommonPions(tr);
        const Trk& a = tr[0];
        const Trk& b = tr[1];
        if (h == 0 && common && lc == 1 && simtruth::FD(set, a) && simtruth::FD(set, b) &&
            (st == 0 ? a.par == 7 : a.g0 == 75101.f)) { in_g[kPi0D] = true; cw[kPi0D] = w; }
        if (h == 0 && common && lc == 2 && simtruth::FD(set, a) && simtruth::FD(set, b) &&
            (st == 0 ? a.par == 17 : a.g0 == 175101.f)) { in_g[kEtaD] = true; cw[kEtaD] = w; }
        if (common && !clone && lc == 7 && simtruth::FD(set, a) && simtruth::FD(set, b) &&
            (st == 0 || (IsPi0EtaForced(a.g0) && IsPi0EtaForced(b.g0)))) { in_g[kDD] = true; cw[kDD] = w; }
        if (h == 0 && common && lc == 3) {
            if (st == 0) {
                if (a.par == 1 && b.par == 1 && a.gp == 7 && b.gp == 7 && simtruth::forciblePi0Code(a.g0) &&
                    simtruth::forciblePi0Code(b.g0)) { in_g[kConvPi0] = true; cw[kConvPi0] = w * kFconvPi0; }
                if (a.par == 1 && b.par == 1 && a.gp == 17 && b.gp == 17 && a.g0 == 1717.f && b.g0 == 1717.f) {
                    in_g[kConvEta] = true; cw[kConvEta] = w / simtruth::fEtaGam * kFconvEta;
                }
            } else {
                if (a.g0 == 75101.f && a.g1 == 7.f && b.g0 == 75101.f && b.g1 == 7.f) { in_g[kConvPi0] = true; cw[kConvPi0] = w * simtruth::kPDG_pi0D; }
                if (a.g0 == 175101.f && a.g1 == 17.f && b.g0 == 175101.f && b.g1 == 17.f) { in_g[kConvEta] = true; cw[kConvEta] = w * simtruth::kPDG_etaD; }
            }
        }
        if (st == 0 ? !ownA : ownA) { in_g[kModeA] = true; cw[kModeA] = w; }
        for (int g = 0; g < kNCG; ++g) {
            if (!in_g[g]) continue;
            C(Form("cl_%s_own%d%d", kCGKey[g], ownA ? 1 : 0, ownB ? 1 : 0));
            for (int s = 0; s < 2; ++s) {
                if (s == 1 && !oa2) continue;
                Fill(cl[g][s][kOTot], 0., cw[g]);
                Fill(cl[g][s][kOMee], mee, cw[g]);
                Fill(cl[g][s][kOOA], oa, cw[g]);
                Fill(cl[g][s][kOM4], m4, cw[g]);
                if (!ghostL) Fill(cl[g][s][kOTotNoGhost], 0., cw[g]);
            }
        }

        // ---- natural factors (gap of w_db) and phi -> eta gamma (bit 14)
        const double nf = simtruth::NaturalFactor(tr, optPdg);
        double dPhiEta = 0;
        bool phiEta = false;
        for (int j = 0; j < 4; ++j) if (PhiEtaProduct(tr[j])) phiEta = true;
        if (phiEta) {
            Trk t2[4];
            for (int j = 0; j < 4; ++j) { t2[j] = tr[j]; if (PhiEtaProduct(t2[j])) t2[j].g0 = 1717.f; }
            const double nf2 = simtruth::NaturalFactor(t2, optPdg);
            dPhiEta = nf > 0 ? w * (nf2 / nf - 1.) : 0.;
            C(Form("phieta_rows_tp%d_own%d", tp, ownA ? 1 : 0));
            if (std::fabs(nf2 / nf - 1.) > 1e-9) C(Form("phieta_changed_tp%d_own%d", tp, ownA ? 1 : 0));
        }
        // ---- bit 2: ordering of the counters of the omega lepton (52051) and the second-decay product (7051)
        if ((flags & simtruth::kFlagChainNotApplied) && st == 1) {
            std::set<float> cOm, cSec;
            for (int j = 0; j < 4; ++j) {
                if (!simtruth::FD(set, tr[j]) || tr[j].g0 != 525107.f) continue;
                if (tr[j].g1 == 52051.f) cOm.insert(tr[j].g2);
                if (tr[j].g1 == 7051.f) cSec.insert(tr[j].g2);
            }
            bool sameKey = false;
            for (float c : cOm) if (cSec.count(c)) sameKey = true;
            C(Form("bit2_h%d_own%d", h, ownA ? 1 : 0));
            if (sameKey) C(Form("bit2_samekey_h%d", h));
            if (cOm.size() > 1) C(Form("bit2_twoOmegaCounters_h%d", h));
        }

        for (int s = 0; s < 2; ++s) {
            if (s == 1 && !oa2) continue;
            // ---- partition variants (every template)
            if (ownA || ownB) {
                const double a0 = ownA ? w : 0.;
                const double bv[5] = {ownB ? w : 0., ownA ? wdbv : 0., ownB ? wdbv : 0., ownA ? wdbv * nf : 0.,
                                      ownB ? wdbv * nf : 0.};
                for (int bn = 0; bn < 2; ++bn)
                    for (int vr = 0; vr < 5; ++vr) {
                        Fill2(pa[tp][s][vr][bn], m4, a0, bv[vr]);
                        if (tmpl <= 2) Fill2(pa[4][s][vr][bn], m4, a0, bv[vr]);
                    }
            }
            if (tmpl == 2) {
                const int lg = lc == 1 ? 0 : lc == 2 ? 1 : lc == 7 ? 2 : 3;
                if ((st == 0 && ownB && !ownA) || (st == 1 && ownA && !ownB))
                    for (int bn = 0; bn < 2; ++bn) Fill(pd[s][lg][bn], m4, w);
            }
            // ---- flags and templates per ownership
            const int oc = ownA ? 0 : (ownB ? 1 : 2);
            Fill(tt[s][tp][oc], m4, w);
            for (int bit = 0; bit < 15; ++bit) {
                if (!(flags & (1 << bit))) continue;
                if (!fl[s][bit][tp][oc])
                    fl[s][bit][tp][oc] = G(Form("fl_s%d_b%d_tp%d_o%d_t%d_h%d", s, bit, tp, oc, st, h), bWin);
                Fill(fl[s][bit][tp][oc], m4, w);
            }
            // ---- gap effects
            for (int bn = 0; bn < 2; ++bn) {
                if ((flags & simtruth::kFlagVmWithFC) && oc == 2) Fill(ef[s][0][tp][bn], m4, w);
                if ((flags & simtruth::kFlagV1CodeGap) && st == 0 && !ownA) Fill(ef[s][1][tp][bn], m4, w);
                if (phiEta && ownA) Fill(ef[s][2][tp][bn], m4, dPhiEta);
                if ((flags & simtruth::kFlagChainNotApplied) && ownA) Fill(ef[s][3][tp][bn], m4, w);
                if (ownA) Fill(ef[s][4][tp][bn], m4, wdbv * (nf - 1.));
                if (clone && ownA) Fill(ef[s][5][tp][bn], m4, w);
                if (ownA) Fill(wdb[s][tp][bn], m4, wdbv);
            }
        }

        // ---- template-1 rows
        if (tmpl == 1) {
            o_cls = cls; o_gl = ghostL ? 1 : 0;
            o_ghost = (tr[0].ghost == 1 ? 1 : 0) | (tr[1].ghost == 1 ? 2 : 0) | (tr[2].ghost == 1 ? 4 : 0) |
                      (tr[3].ghost == 1 ? 8 : 0);
            o_oa2 = oa2 ? 1 : 0; o_grp = (int)grp; o_set = st + 1; o_ownA = ownA ? 1 : 0;
            o_m4 = (float)m4; o_m4t = in.m4_true; o_mee = (float)mee; o_oa = (float)oa;
            o_w = w; o_wdb = wdbv;
            sig->Fill();
        }
        if (i % 1000000 == 0 && i > 0) { printf("  [%d] %lld / %lld rows, %.0f s\n", k, i, n, sw.RealTime()); sw.Continue(); fflush(stdout); }
    }
    FlushGroup();
    C("groups", (double)(grp + 1));

    // ---- write
    fout->cd();
    sig->Write();
    for (auto& kv : gReg) {
        TVectorD vs((int)kv.second->s.size(), kv.second->s.data());
        TVectorD vv((int)kv.second->v.size(), kv.second->v.data());
        TVectorD vn((int)kv.second->n.size(), kv.second->n.data());
        fout->WriteObject(&vs, (kv.first + "__s").c_str());
        fout->WriteObject(&vv, (kv.first + "__v").c_str());
        fout->WriteObject(&vn, (kv.first + "__n").c_str());
    }
    for (auto& kv : pReg) {
        const PAcc* p = kv.second;
        const std::vector<double>* arr[6] = {&p->a, &p->bb, &p->aa, &p->b2, &p->ab, &p->n};
        const char* suf[6] = {"__a", "__b", "__aa", "__bb", "__ab", "__n"};
        for (int j = 0; j < 6; ++j) {
            TVectorD v((int)arr[j]->size(), arr[j]->data());
            fout->WriteObject(&v, (kv.first + suf[j]).c_str());
        }
    }
    {
        TTree* ct = new TTree("codes", "stored code values: column, slot type (0 lepton, 1 pion), value, rows");
        Int_t c_col, c_styp, c_tree = k;
        Float_t c_val;
        Long64_t c_cnt;
        ct->Branch("tree", &c_tree); ct->Branch("col", &c_col); ct->Branch("styp", &c_styp);
        ct->Branch("val", &c_val); ct->Branch("cnt", &c_cnt);
        for (int c = 0; c < 3; ++c)
            for (int s = 0; s < 2; ++s)
                for (auto& kv : codes[c][s]) {
                    c_col = c; c_styp = s; std::memcpy(&c_val, &kv.first, 4); c_cnt = kv.second;
                    ct->Fill();
                }
        ct->Write();
    }
    {
        std::string txt;
        for (auto& kv : ctr) txt += Form("%s\t%.17g\n", kv.first.c_str(), kv.second);
        TObjString os(txt.c_str());
        fout->WriteObject(&os, "counters");
    }
    fout->Close();
    delete fout;
    fin->Close();
    delete fin;
    printf("[pass %d] %s: %lld rows, %lld groups, %.0f s, %zu + %zu accumulators\n", k, kInName[k], n, grp + 1,
           sw.RealTime(), gReg.size(), pReg.size());
    for (auto& kv : gReg) delete kv.second;
    for (auto& kv : pReg) delete kv.second;
    gReg.clear(); pReg.clear();
}

}  // namespace scl

#include "sim_closure_analysis.h"

void sim_closure(int stage = 0, int tree = -1) {
    scl::InitBins();
    if (stage == 0 || stage == 1) {
        for (int k = 0; k < 6; ++k) if (tree < 0 || tree == k) scl::Pass(k);
    }
    if (stage == 0 || stage == 2) scl::Analyse();
}

