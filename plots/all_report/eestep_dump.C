// ============================================================================
// eestep_dump.C  --  package eestep (origin of the step at u = 0 in the like-sign e-e- rows): data row dump with
//                    track-level PAT quantities and the links between the e+e-, e+e+ and e-e- rows of one event
// ============================================================================
// Written 2026-11-15 (prefix eestep). Read-only on every input. Inputs:
//   - the data compact rows (scan_v2, ScanEngine.h), selection exactly as edgeu_dump.C: cor_ok == 1, "pre",
//     0.26 <= M(pipi) < 0.46, 0.30 <= M4 < 0.85;
//   - the PAT files /hdd3/przygoda/HADES/files/exp/gen4/lep_<day>.root (the inputs of the production), trees
//     PipPimEpEm_ID (prefixes ep, em), PipPimEpEp_ID (ep1, ep2), PipPimEmEm_ID (em1, em2); a compact row is the PAT
//     entry (file_id = day, pat_entry) of its own tree (main.cc: entry index inside the file).
// Event identity across the three PAT trees: (runnumber, event, start_time, eVertClust_z, totalmult). The event counter
// alone is not unique (it restarts per sub-file); the two reconstructed floats make the key unique (checked below:
// every key must be ONE block of adjacent entries in each tree, otherwise the macro reports it).
// Lepton / pion identity across trees: equal (theta, phi, p) of the PAT candidate (same float values).
//
// Output tables (flat float32 column files, <out>/<table>_<col>.f32, column list in <out>/eestep_layout.txt):
//   ls   every selected like-sign row (k = 1 e+e+, 2 e-e-), all days;
//   usr  the selected e+e- rows with 0.30 <= M(pipi) < 0.40 and -0.10 < u < 0.10 (u = M4 - alldefs::Edge(M(pipi), M(ee)));
// columns: compact quantities, PAT details of the two lepton slots (l1 = ep / ep1 / em1, l2 = em / ep2 / em2) and of the
// event, and the link quantities (see the comments of Link below).
//
// Run (from /home/przygoda/HADES/FAT/FAT; about 6 GB, 10-20 min, one process):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_eestep", kTRUE)' \
//        'plots/all_report/eestep_dump.C+O'
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#include "TFile.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TTree.h"

#include "bedge_common.h"

namespace eestep {
using namespace bedge;
const char* const kOutE = "/hdd3/przygoda/HADES/FAT_out/all_report/eestep";
const char* const kPat = "/hdd3/przygoda/HADES/files/exp/gen4/lep_%03d.root";
FILE* gLog = nullptr;
#define PL(...) do { printf(__VA_ARGS__); if (gLog) fprintf(gLog, __VA_ARGS__); fflush(stdout); } while (0)

// ---------------------------------------------------------------- PAT columns
const char* const kTV[] = {"q", "sector", "system", "p", "p_org", "theta", "phi", "mdcinnerchi2", "mdcouterchi2", "rkchi2",
    "metamatchquality", "metamatchradius", "rkmeta_dx", "rkmeta_dy", "rich_ind", "rich_padnum", "rich_amp", "rich_centr",
    "rich_patmat", "rich_houtra", "richmatchquality", "richmatchqualitynorm", "ringchi2", "delta_theta", "delta_phi",
    "beta_new", "mass2", "tof", "track_length", "ismdcedge", "isBT", "kIsUsed", "kIsLepton", "resolution", "z", "r",
    "dedx_mdc", "isring", "isringtrack", "emc_energy", "isemc", "oa_lept", "oa_hadr", "isUsedVertex", "ismdc3edge",
    "ismdc4edge", "isringnomatch", "isrichmatch_richrk", "isnewrich", "rich_radius", "metaselected"};
const int kNTV = sizeof(kTV) / sizeof(kTV[0]);
const char* const kEV[] = {"lpos_mult", "lneg_mult", "hpos_mult", "hneg_mult", "evtLepMult", "totalmult", "evtVertCand",
                           "trigbit", "isBest", "evtStart"};
const int kNEV = sizeof(kEV) / sizeof(kEV[0]);

struct Trk { float theta, phi, p, ring, in, out, tof, sector, q; };
struct LRow {   // link quantities of one PAT entry
    float run, ev, st, vzc, tm;
    Trk l[2];
    float pip[3], pim[3];   // theta, phi, p
    float isBest;
};
inline bool SameT(const Trk& a, const Trk& b) { return a.theta == b.theta && a.phi == b.phi && a.p == b.p; }
inline bool Same3(const float* a, const float* b) { return a[0] == b[0] && a[1] == b[1] && a[2] == b[2]; }
inline double Ang(const Trk& a, const Trk& b) {   // degrees
    const double d2r = M_PI / 180.;
    const double ta = a.theta * d2r, pa = a.phi * d2r, tb = b.theta * d2r, pb = b.phi * d2r;
    const double c = std::sin(ta) * std::sin(tb) * std::cos(pa - pb) + std::cos(ta) * std::cos(tb);
    return std::acos(std::max(-1., std::min(1., c))) / d2r;
}
std::string KeyOf(const LRow& r) {
    char b[64];
    uint32_t s, z;
    std::memcpy(&s, &r.st, 4); std::memcpy(&z, &r.vzc, 4);
    snprintf(b, sizeof(b), "%d_%d_%08x_%08x_%d", (int)r.run, (int)r.ev, s, z, (int)r.tm);
    return b;
}

// read the link quantities of every entry of one PAT tree
bool ReadLinks(TTree* t, const char* p1, const char* p2, std::vector<LRow>& out) {
    t->SetBranchStatus("*", 0);
    LRow r;
    auto B = [&](const char* n, float* a) {
        if (!t->GetBranch(n)) { printf("missing branch %s\n", n); return false; }
        t->SetBranchStatus(n, 1); t->SetBranchAddress(n, a); return true;
    };
    bool ok = B("runnumber", &r.run) && B("event", &r.ev) && B("start_time", &r.st) && B("eVertClust_z", &r.vzc) &&
              B("totalmult", &r.tm) && B("isBest", &r.isBest);
    const char* pre[2] = {p1, p2};
    for (int s = 0; s < 2; ++s) {
        ok = ok && B(Form("%s_theta", pre[s]), &r.l[s].theta) && B(Form("%s_phi", pre[s]), &r.l[s].phi) &&
             B(Form("%s_p", pre[s]), &r.l[s].p) && B(Form("%s_rich_ind", pre[s]), &r.l[s].ring) &&
             B(Form("%s_mdcinnerchi2", pre[s]), &r.l[s].in) && B(Form("%s_mdcouterchi2", pre[s]), &r.l[s].out) &&
             B(Form("%s_tof", pre[s]), &r.l[s].tof) && B(Form("%s_sector", pre[s]), &r.l[s].sector) &&
             B(Form("%s_q", pre[s]), &r.l[s].q);
    }
    ok = ok && B("pip_theta", &r.pip[0]) && B("pip_phi", &r.pip[1]) && B("pip_p", &r.pip[2]) &&
         B("pim_theta", &r.pim[0]) && B("pim_phi", &r.pim[1]) && B("pim_p", &r.pim[2]);
    if (!ok) return false;
    const Long64_t n = t->GetEntries();
    out.resize(n);
    for (Long64_t i = 0; i < n; ++i) { t->GetEntry(i); out[i] = r; }
    t->ResetBranchAddresses();
    return true;
}
// event blocks: key -> entries; also check that each key is one adjacent block
void Index(const std::vector<LRow>& v, std::unordered_map<std::string, std::vector<int>>& m, long& nSplit) {
    m.clear(); m.reserve(v.size());
    std::string last;
    for (size_t i = 0; i < v.size(); ++i) {
        const std::string k = KeyOf(v[i]);
        auto it = m.find(k);
        if (it != m.end() && k != last) ++nSplit;   // the key re-appears after another key: not adjacent
        m[k].push_back((int)i);
        last = k;
    }
}

// ---------------------------------------------------------------- table
struct Tab {
    std::string name;
    std::vector<std::string> cols;
    std::map<std::string, int> idx;
    std::vector<std::vector<float>> v;
    std::vector<float> row;
    size_t nrows = 0;
    // a column created after nrows rows starts with nrows NaN (no late-column shift, cf. DynamicHNtuple)
    int C(const std::string& c) {
        auto it = idx.find(c);
        if (it != idx.end()) return it->second;
        const int i = (int)cols.size(); cols.push_back(c); idx[c] = i; v.emplace_back(nrows, NAN); row.push_back(NAN); return i;
    }
    void Set(int i, double x) { row[i] = (float)x; }
    void Push() {
        for (size_t i = 0; i < cols.size(); ++i) v[i].push_back(row[i]);
        std::fill(row.begin(), row.end(), NAN);
        ++nrows;
    }
    size_t N() const { return nrows; }
    void Write(FILE* lay) {
        for (size_t i = 0; i < cols.size(); ++i) {
            FILE* fp = fopen(Form("%s/%s_%s.f32", kOutE, name.c_str(), cols[i].c_str()), "wb");
            if (!fp) { PL("cannot write %s\n", cols[i].c_str()); continue; }
            fwrite(v[i].data(), 4, v[i].size(), fp); fclose(fp);
        }
        fprintf(lay, "%s %zu %zu\n", name.c_str(), N(), cols.size());
        for (auto& c : cols) fprintf(lay, "  %s\n", c.c_str());
    }
};
}  // namespace eestep

void eestep_dump() {
    using namespace eestep;
    TStopwatch total;
    gSystem->mkdir(kOutE, kTRUE);
    gLog = fopen(Form("%s/eestep_dump.log", kOutE), "w");
    PL("# eestep_dump.C, 2026-11-15: like-sign and edge-region e+e- rows with PAT track details and event links\n");
    scan::Data d;
    if (!scan::Load(d)) return;
    // selected compact rows per day and sample: (index in d.s[k])
    std::map<int, std::vector<std::pair<int, int>>> perDay;   // day -> (k, i)
    long nk[3] = {0, 0, 0};
    for (int k = 0; k < 3; ++k)
        for (size_t i = 0; i < d.s[k].size(); ++i) {
            const auto& r = d.s[k][i];
            if (r.cor_ok != 1) continue;
            if (!PassSel(kSelPre, r, true, r.dz_pi, r.mm4)) continue;
            if (!(r.mpipi >= 0.26f && r.mpipi < 0.46f)) continue;
            if (!(r.m4 >= 0.30f && r.m4 < 0.85f)) continue;
            perDay[r.file_id].push_back({k, (int)i});
            ++nk[k];
        }
    PL("selected compact rows (as edgeu_dump.C): us %ld pp %ld mm %ld\n", nk[0], nk[1], nk[2]);

    Tab T[2];
    T[0].name = "ls"; T[1].name = "usr";
    const char* trees[3] = {"PipPimEpEm_ID", "PipPimEpEp_ID", "PipPimEmEm_ID"};
    const char* pfx[3][2] = {{"ep", "em"}, {"ep1", "ep2"}, {"em1", "em2"}};
    long nSplitTot[3] = {0, 0, 0}, nChk = 0, nBad = 0;
    double maxdp = 0;
    for (auto& dd : perDay) {
        const int day = dd.first;
        TStopwatch sw;
        TFile* f = TFile::Open(Form(kPat, day), "READ");
        if (!f || f->IsZombie()) { PL("cannot open day %d\n", day); continue; }
        std::vector<LRow> L[3];
        std::unordered_map<std::string, std::vector<int>> M[3];
        for (int k = 0; k < 3; ++k) {
            TTree* t = (TTree*)f->Get(trees[k]);
            if (!t || !ReadLinks(t, pfx[k][0], pfx[k][1], L[k])) { PL("day %d: tree %s unreadable\n", day, trees[k]); return; }
            long ns = 0;
            Index(L[k], M[k], ns);
            nSplitTot[k] += ns;
        }
        // selected compact row of each PAT entry (us, pp, mm)
        std::vector<int> sel[3];
        for (int k = 0; k < 3; ++k) sel[k].assign(L[k].size(), -1);
        for (auto& ki : dd.second) {
            const auto& r = d.s[ki.first][ki.second];
            if (r.pat_entry >= 0 && r.pat_entry < (int)sel[ki.first].size()) sel[ki.first][r.pat_entry] = ki.second;
        }
        // rows that get a detailed entry: all like-sign rows; e+e- rows in the edge region
        std::vector<std::pair<int, int>> det[3];   // (pat_entry, compact index)
        for (auto& ki : dd.second) {
            const auto& r = d.s[ki.first][ki.second];
            if (ki.first == 0) {
                if (!(r.mpipi >= 0.30f && r.mpipi < 0.40f) || !(r.mee < scan::kMpi0)) continue;
                const double e = alldefs::Edge(r.mpipi, r.mee);
                if (e <= 0) continue;
                const double u = r.m4 - e;
                if (!(u > -0.10 && u < 0.10)) continue;
            }
            det[ki.first].push_back({r.pat_entry, ki.second});
        }
        for (int k = 0; k < 3; ++k) {
            std::sort(det[k].begin(), det[k].end());
            if (det[k].empty()) continue;
            TTree* t = (TTree*)f->Get(trees[k]);
            t->SetBranchStatus("*", 0);
            float tv[2][kNTV], ev[kNEV], pp[2][4];
            for (int s = 0; s < 2; ++s)
                for (int j = 0; j < kNTV; ++j) {
                    const char* bn = Form("%s_%s", pfx[k][s], kTV[j]);
                    t->SetBranchStatus(bn, 1); t->SetBranchAddress(bn, &tv[s][j]);
                }
            for (int j = 0; j < kNEV; ++j) { t->SetBranchStatus(kEV[j], 1); t->SetBranchAddress(kEV[j], &ev[j]); }
            const char* pn[2] = {"pip", "pim"};
            const char* pv[4] = {"sector", "q", "kIsLepton", "isring"};
            for (int s = 0; s < 2; ++s)
                for (int j = 0; j < 4; ++j) {
                    const char* bn = Form("%s_%s", pn[s], pv[j]);
                    t->SetBranchStatus(bn, 1); t->SetBranchAddress(bn, &pp[s][j]);
                }
            Tab& O = T[k == 0 ? 1 : 0];
            for (auto& e : det[k]) {
                t->GetEntry(e.first);
                const auto& r = d.s[k][e.second];
                const LRow& me = L[k][e.first];
                // ---- compact quantities
                O.Set(O.C("k"), k); O.Set(O.C("day"), day); O.Set(O.C("pat_entry"), r.pat_entry); O.Set(O.C("run"), r.run);
                O.Set(O.C("grp"), r.grp); O.Set(O.C("isBest"), r.isBest); O.Set(O.C("nrow"), r.nrow);
                O.Set(O.C("n_lep1"), r.n_lep1); O.Set(O.C("n_lep2"), r.n_lep2); O.Set(O.C("n_lep"), r.n_lep);
                O.Set(O.C("n_pip"), r.n_pip); O.Set(O.C("n_pim"), r.n_pim);
                O.Set(O.C("s0"), PassSel(kSelS0, r, true, r.dz_pi, r.mm4) ? 1 : 0);
                O.Set(O.C("m4"), r.m4); O.Set(O.C("mee"), r.mee); O.Set(O.C("mpipi"), r.mpipi); O.Set(O.C("oa"), r.oa);
                O.Set(O.C("oa_lab"), r.oa_lab); O.Set(O.C("mm4"), r.mm4); O.Set(O.C("pt4"), r.pt4); O.Set(O.C("y4"), r.y4);
                O.Set(O.C("vz"), r.vz); O.Set(O.C("vchi2"), r.vchi2); O.Set(O.C("dz_ee"), r.dz_ee); O.Set(O.C("r_ee"), r.r_ee);
                O.Set(O.C("doca_ee"), r.doca_ee); O.Set(O.C("shared_ring"), r.shared_ring); O.Set(O.C("hyp_chi2"), r.hyp_chi2);
                O.Set(O.C("l1_pcor"), r.ep_p); O.Set(O.C("l2_pcor"), r.em_p); O.Set(O.C("l1_prec"), r.ep_p_rec); O.Set(O.C("l2_prec"), r.em_p_rec);
                O.Set(O.C("pip_p"), r.pip_p); O.Set(O.C("pim_p"), r.pim_p); O.Set(O.C("pip_theta"), r.pip_theta); O.Set(O.C("pim_theta"), r.pim_theta);
                O.Set(O.C("pip_phi"), r.pip_phi); O.Set(O.C("pim_phi"), r.pim_phi); O.Set(O.C("pip_mass2"), r.pip_mass2); O.Set(O.C("pim_mass2"), r.pim_mass2);
                O.Set(O.C("n_gam"), r.n_gam); O.Set(O.C("g_e"), r.g_e); O.Set(O.C("neutr_mult"), r.neutr_mult);
                // ---- PAT details
                for (int s = 0; s < 2; ++s)
                    for (int j = 0; j < kNTV; ++j) O.Set(O.C(Form("l%d_%s", s + 1, kTV[j])), tv[s][j]);
                for (int j = 0; j < kNEV; ++j) O.Set(O.C(Form("ev_%s", kEV[j])), ev[j]);
                for (int s = 0; s < 2; ++s)
                    for (int j = 0; j < 4; ++j) O.Set(O.C(Form("%s_%s", pn[s], pv[j])), pp[s][j]);
                // consistency: compact X_p_rec == PAT X_p
                ++nChk;
                const double dp = std::max(std::fabs(r.ep_p_rec - tv[0][3]), std::fabs(r.em_p_rec - tv[1][3]));
                maxdp = std::max(maxdp, dp);
                if (dp > 1e-3) ++nBad;
                // ---- links: the PAT entries of the same event in the three trees
                const std::string key = KeyOf(me);
                int nev[3], nsel[3];
                for (int q = 0; q < 3; ++q) {
                    auto it = M[q].find(key);
                    nev[q] = it == M[q].end() ? 0 : (int)it->second.size();
                    nsel[q] = 0;
                    if (it != M[q].end()) for (int j : it->second) if (sel[q][j] >= 0) ++nsel[q];
                }
                O.Set(O.C("ev_n_us"), nev[0]); O.Set(O.C("ev_n_pp"), nev[1]); O.Set(O.C("ev_n_mm"), nev[2]);
                O.Set(O.C("ev_nsel_us"), nsel[0]); O.Set(O.C("ev_nsel_pp"), nsel[1]); O.Set(O.C("ev_nsel_mm"), nsel[2]);
                // own-tree entries of the event: how many distinct rows, own row rank
                // Link to the other-charge trees. For a like-sign row (k = 1 or 2) the partner tree is e+e- (q = 0):
                //   lk_type 0: no e+e- entry in the event; 1: e+e- entries, none with the same pi+, pi- and one shared lepton;
                //   2: an e+e- entry with the same pi+, pi- and one lepton slot equal ("twin"); lk_slot = 1 / 2 = which slot of
                //   the like-sign row is shared; the other slot ("odd lepton") is compared with the twin's opposite-charge
                //   lepton: lk_ang (deg), lk_prat (p_odd / p_partner), lk_ring / lk_in / lk_out / lk_tof (equal RICH ring index,
                //   inner-segment chi2, outer-segment chi2, META time), lk_m4 / lk_mee / lk_mpipi of the twin if selected.
                // For an e+e- row (k = 0) the partner trees are e-e- (q = 2, "mm") and e+e+ (q = 1, "pp"): the same quantities
                //   with the prefixes lkmm_ and lkpp_ (the twin shares pi+, pi- and the e- (mm) / e+ (pp) of the e+e- row).
                auto link = [&](int q, const char* tag, int keepSlotOfOwnUS) {
                    // keepSlotOfOwnUS: for k = 0 rows, which own slot must be shared (0 = ep, 1 = em); -1 for like-sign rows
                    auto it = M[q].find(key);
                    int type = it == M[q].end() ? 0 : 1, slot = 0;
                    double bestAng = 1e9; int bestJ = -1, bestOdd = -1, bestPartner = -1;
                    double minAng[2] = {1e9, 1e9}; int ringAny[2] = {0, 0};
                    if (it != M[q].end()) {
                        for (int j : it->second) {
                            const LRow& o = L[q][j];
                            if (k != 0) {   // like-sign own row vs e+e- entries: opposite-charge lepton of o is l[0] (ep) for mm, l[1] (em) for pp
                                const int opp = (k == 2) ? 0 : 1, same = 1 - opp;
                                for (int s = 0; s < 2; ++s) {
                                    const double a = Ang(me.l[s], o.l[opp]);
                                    minAng[s] = std::min(minAng[s], a);
                                    if (me.l[s].ring == o.l[opp].ring && me.l[s].ring >= 0) ringAny[s] = 1;
                                }
                                if (!Same3(me.pip, o.pip) || !Same3(me.pim, o.pim)) continue;
                                for (int s = 0; s < 2; ++s) {
                                    if (!SameT(me.l[s], o.l[same])) continue;
                                    const double a = Ang(me.l[1 - s], o.l[opp]);
                                    if (a < bestAng) { bestAng = a; bestJ = j; bestOdd = 1 - s; bestPartner = opp; slot = s + 1; }
                                }
                            } else {        // e+e- own row vs like-sign entries (q = 1 pp: shares the e+; q = 2 mm: shares the e-)
                                const int shared = (q == 2) ? 1 : 0;   // own slot that must appear in o
                                const int oddOwn = 1 - shared;         // own lepton that would be "flipped"
                                for (int s = 0; s < 2; ++s) {
                                    const double a = Ang(me.l[oddOwn], o.l[s]);
                                    minAng[0] = std::min(minAng[0], a);
                                    if (me.l[oddOwn].ring == o.l[s].ring && me.l[oddOwn].ring >= 0) ringAny[0] = 1;
                                }
                                if (!Same3(me.pip, o.pip) || !Same3(me.pim, o.pim)) continue;
                                for (int s = 0; s < 2; ++s) {
                                    if (!SameT(me.l[shared], o.l[s])) continue;
                                    const double a = Ang(me.l[oddOwn], o.l[1 - s]);
                                    if (a < bestAng) { bestAng = a; bestJ = j; bestOdd = 1 - s; bestPartner = oddOwn; slot = s + 1; }
                                }
                            }
                        }
                    }
                    (void)keepSlotOfOwnUS;
                    if (bestJ >= 0) type = 2;
                    O.Set(O.C(Form("%s_type", tag)), type); O.Set(O.C(Form("%s_slot", tag)), slot);
                    O.Set(O.C(Form("%s_minang1", tag)), minAng[0] < 1e8 ? minAng[0] : NAN);
                    O.Set(O.C(Form("%s_minang2", tag)), minAng[1] < 1e8 ? minAng[1] : NAN);
                    O.Set(O.C(Form("%s_ringany1", tag)), ringAny[0]); O.Set(O.C(Form("%s_ringany2", tag)), ringAny[1]);
                    if (bestJ >= 0) {
                        const LRow& o = L[q][bestJ];
                        // odd: for like-sign rows = own slot bestOdd vs o.l[bestPartner]; for e+e- rows = o.l[bestOdd] vs own l[bestPartner]
                        const Trk& A = (k != 0) ? me.l[bestOdd] : o.l[bestOdd];
                        const Trk& Bt = (k != 0) ? o.l[bestPartner] : me.l[bestPartner];
                        O.Set(O.C(Form("%s_ang", tag)), bestAng);
                        O.Set(O.C(Form("%s_prat", tag)), A.p / Bt.p);
                        O.Set(O.C(Form("%s_ring", tag)), A.ring == Bt.ring ? 1 : 0);
                        O.Set(O.C(Form("%s_in", tag)), A.in == Bt.in ? 1 : 0);
                        O.Set(O.C(Form("%s_out", tag)), A.out == Bt.out ? 1 : 0);
                        O.Set(O.C(Form("%s_tof", tag)), A.tof == Bt.tof ? 1 : 0);
                        O.Set(O.C(Form("%s_secsame", tag)), A.sector == Bt.sector ? 1 : 0);
                        O.Set(O.C(Form("%s_isBest", tag)), o.isBest);
                        const int cs = sel[q][bestJ];
                        if (cs >= 0) {
                            const auto& rr = d.s[q][cs];
                            O.Set(O.C(Form("%s_m4", tag)), rr.m4); O.Set(O.C(Form("%s_mee", tag)), rr.mee); O.Set(O.C(Form("%s_mpipi", tag)), rr.mpipi);
                        }
                    } else {
                        for (const char* c : {"ang", "prat", "ring", "in", "out", "tof", "secsame", "isBest", "m4", "mee", "mpipi"})
                            O.C(Form("%s_%s", tag, c));
                    }
                };
                if (k != 0) link(0, "lk", -1);
                else { link(2, "lkmm", 1); link(1, "lkpp", 0); }
                O.Push();
            }
            t->ResetBranchAddresses();
        }
        f->Close(); delete f;
        PL("day %03d: PAT entries us %zu pp %zu mm %zu; detailed rows us %zu pp %zu mm %zu; %.0f s\n", day, L[0].size(), L[1].size(),
           L[2].size(), det[0].size(), det[1].size(), det[2].size(), sw.RealTime());
    }
    PL("event key re-appearing after another key (not one adjacent block): us %ld pp %ld mm %ld\n", nSplitTot[0], nSplitTot[1], nSplitTot[2]);
    PL("compact X_p_rec vs PAT X_p: %ld rows checked, %ld differ by > 1e-3, max |diff| %.3g\n", nChk, nBad, maxdp);
    FILE* lay = fopen(Form("%s/eestep_layout.txt", kOutE), "w");
    for (int i = 0; i < 2; ++i) { T[i].Write(lay); PL("table %s: %zu rows, %zu columns\n", T[i].name.c_str(), T[i].N(), T[i].cols.size()); }
    fclose(lay);
    PL("done in %.0f s\n", total.RealTime());
    if (gLog) fclose(gLog);
}
