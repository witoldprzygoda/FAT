// mixing_build_betanew.C
// =========================================================================
// Event mixing for the e+e- combinatorial-background study (betanew, pp45).
//
// Input : output_epem_exp.root, output_epep_exp.root, output_emem_exp.root
//         (TNtuple "dilepton_nt", Float_t columns; ALL pair combinations,
//         no isBest cut; PT3 = trigbit 8192, PT2 = trigbit 4096).
//         Column prefix note: in the like-sign files ep_* is leg 1 and
//         em_* is leg 2 of the SAME charge (epep: both e+, emem: both e-).
//
// Data layout facts (checked on the full files, 2026-10-31):
//   * runnumber is non-decreasing in all three files (0 violations) and
//     the run order is the same in all files; rows of one event are
//     adjacent.
//   * Inside one runnumber there are (up to) 10 blocks (sub-files); the
//     event number restarts in every block, so (runnumber, event) is NOT a
//     unique key (0.54% of e+e- events share the key with another event of
//     the same run). The like-sign files are sparse, so block boundaries
//     cannot be recovered there reliably.
//   * In the all-combinations scheme every leg of an event with n+>=1 and
//     n->=1 appears in e+e- rows (every e+ is paired with every e-). The
//     like-sign files therefore add only (a) the isBest flag of the
//     like-sign hypothesis and (b) events with n+=0 or n-=0 (like-sign
//     only), whose mixing classes are disjoint from those of e+e- events.
//
// Method (streaming, run by run; memory = events of one run + pools):
//   1. All e+e- events of a run are loaded (file order = time order).
//      Legs e+ / e- are rebuilt from the pair rows, unique by float equality
//      of (p_rec, theta, phi); each leg carries best_pm (member of the
//      isBest==1 e+e- pair) and best_ls (member of the isBest==1 like-sign
//      pair). n+ / n- = numbers of unique legs. Events with a single lepton
//      of one sign and none of the other do not exist (pairs only).
//   2. Like-sign events of the run (adjacent rows with equal (run, event))
//      are matched to e+e- events with the same event number by leg
//      kinematics (float equality); the candidate with the largest number
//      of matching legs wins (>= 1 required). Matched: best_ls flags are
//      transferred, unmatched legs added (counted). Unmatched events are
//      like-sign-only events (n-=0 or n+=0). Every event gets a block
//      index (e+e- events: block boundary = decrease of the event number,
//      reliable at ~270 events per block; like-sign events: own decreases,
//      re-synchronised at every matched event) and all events of the run
//      are processed in (block, event) order, i.e. in time order.
//   3. Same-event spectra from m_ee rows, identical to the recipe macro
//      mass_ee_pt_ratio_betanew.C: h_se_{pm,pp,mm}[_best]_{PT3,PT2} on the
//      58-bin recipe grid and *_fine on 280 x 5 MeV [0, 1.4].
//   4. Mixing: class = trigger x n+ class {0,1,2,>=3} x n- class {0,1,2,>=3}
//      with a rolling buffer of the last pool_depth events of the class
//      (no day class: the buffer keeps time locality by construction).
//      New event E vs each buffered B: +- = (e+ E, e- B) and (e- E, e+ B),
//      ++ = (e+ E, e+ B), -- = (e- E, e- B). Mass as PParticle (REC):
//      m^2 = 2 m_e^2 + 2 (E1 E2 - p1.p2), m_e = 0.51099895 MeV.
//      Variant "best legs": only legs with best_pm (for +-) / best_ls
//      (for ++, --) -> h_mix_*_best_*; same class definition (all legs).
//      Second pool WITHOUT multiplicity classes (trigger only, same
//      pool_depth) -> h_mixnc_*: inside a fixed (n+, n-) class k = 1 by
//      construction and the class-summed k only encodes the (n+, n-)
//      population (like-sign-only events), not acceptance; the no-class
//      pool measures the acceptance/shape effect alone.
//   5. All spectra with Sumw2 (weight 1 -> error sqrt(N)) + TParameters.
//
// max_rows >= 0 : read only the first max_rows rows of the e+e- file and
//      keep only COMPLETE runs (the run cut by max_rows is dropped); the
//      like-sign files are read up to the same last run. The same-event
//      spectra then equal TTree::Draw with "runnumber <= last_run".
//
// Usage (ACLiC):
//   root -l -b -q 'plots/mixing_build_betanew.C+(-1,"plots/output/mixing_betanew.root",10,true)'
//
// @author Witold Przygoda (witold.przygoda@uj.edu.pl)
// @date 2026

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TAxis.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TSystem.h>
#include <TStopwatch.h>
#include <TMath.h>
#include <TString.h>

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdio>
#include <deque>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

namespace mixb {

constexpr const char* kFileEm = "output_epem_exp.root";
constexpr const char* kFilePp = "output_epep_exp.root";
constexpr const char* kFileMm = "output_emem_exp.root";
constexpr const char* kTree   = "dilepton_nt";

constexpr double kMe      = 0.51099895;   // MeV
constexpr int    kNb58    = 58;
constexpr int    kNbFine  = 280;
constexpr double kXmin    = 0.0;
constexpr double kXmax    = 1.4;
constexpr int    kNTrig   = 2;             // 0 = PT3 (8192), 1 = PT2 (4096)
constexpr int    kNMultCl = 4;             // 0, 1, 2, >=3
constexpr long   kDayGap  = 50000;         // runnumber jump defining a "day"
constexpr long   kDaySec  = 86400;         // alternative: runnumber / 86400

const char* trigName(int t) { return t == 0 ? "PT3" : "PT2"; }

// Same construction as plots/mass_ee_pt_ratio_betanew.C::getEdges().
const double* getEdges58() {
    static double edges[kNb58 + 1];
    static bool built = false;
    if (!built) {
        int k = 0;
        for (int i = 0; i <= 20; ++i) edges[k++] = i * 0.0025;
        edges[k++] = 0.0800;
        edges[k++] = 0.1100;
        edges[k++] = 0.1400;
        for (int i = 1; i <= 33; ++i) edges[k++] = 0.140 + i * 0.020;
        edges[k++] = 0.950;
        edges[k++] = 1.400;
        built = true;
    }
    return edges;
}

// ---------------------------------------------------------------- spectra
// One mass spectrum on both grids. Counting is done in plain arrays
// (bin index via TAxis::FindBin, identical to TH1::Fill) and transferred to
// the TH1D at the end (content = N, error = sqrt(N), as Sumw2 with w = 1).
struct Spectrum {
    TH1D* h58   = nullptr;
    TH1D* hfine = nullptr;
    std::vector<double> c58, cfine;
    long long n = 0;

    void init(const std::string& name, const std::string& title) {
        h58 = new TH1D(name.c_str(), title.c_str(), kNb58, getEdges58());
        h58->Sumw2();
        h58->SetDirectory(nullptr);
        hfine = new TH1D((name + "_fine").c_str(), title.c_str(),
                         kNbFine, kXmin, kXmax);
        hfine->Sumw2();
        hfine->SetDirectory(nullptr);
        c58.assign(kNb58 + 2, 0.0);
        cfine.assign(kNbFine + 2, 0.0);
    }
    inline void fill(double m) {
        c58[h58->GetXaxis()->FindBin(m)]     += 1.0;
        cfine[hfine->GetXaxis()->FindBin(m)] += 1.0;
        ++n;
    }
    void finalize() {
        for (int b = 0; b <= kNb58 + 1; ++b) {
            h58->SetBinContent(b, c58[b]);
            h58->SetBinError(b, std::sqrt(c58[b]));
        }
        for (int b = 0; b <= kNbFine + 1; ++b) {
            hfine->SetBinContent(b, cfine[b]);
            hfine->SetBinError(b, std::sqrt(cfine[b]));
        }
        h58->SetEntries(static_cast<double>(n));
        hfine->SetEntries(static_cast<double>(n));
    }
    void write() { h58->Write(); hfine->Write(); }
};

// ---------------------------------------------------------------- legs
struct Leg {
    float  p, th, ph;          // identity (REC momentum MeV/c, deg, deg)
    double px, py, pz, E;      // MeV
    bool   best_pm = false;    // belongs to the isBest==1 e+e- pair
    bool   best_ls = false;    // belongs to the isBest==1 like-sign pair
};

inline void setKinematics(Leg& l) {
    const double p  = l.p;
    const double th = l.th * TMath::DegToRad();
    const double ph = l.ph * TMath::DegToRad();
    const double st = std::sin(th);
    l.px = p * st * std::cos(ph);
    l.py = p * st * std::sin(ph);
    l.pz = p * std::cos(th);
    l.E  = std::sqrt(p * p + kMe * kMe);
}

// Invariant mass [GeV/c^2] of two legs, as PParticle (REC kinematics).
inline double massGeV(const Leg& a, const Leg& b) {
    const double m2 = 2.0 * kMe * kMe +
        2.0 * (a.E * b.E - (a.px * b.px + a.py * b.py + a.pz * b.pz));
    return (m2 > 0.0 ? std::sqrt(m2) : 0.0) * 1e-3;
}

inline int findLeg(const std::vector<Leg>& legs, float p, float th, float ph) {
    for (size_t i = 0; i < legs.size(); ++i)
        if (legs[i].p == p && legs[i].th == th && legs[i].ph == ph)
            return static_cast<int>(i);
    return -1;
}

// Insert leg into a per-event list, unique by float equality; returns index.
int addLeg(std::vector<Leg>& legs, float p, float th, float ph) {
    const int k = findLeg(legs, p, th, ph);
    if (k >= 0) return k;
    Leg l;
    l.p = p; l.th = th; l.ph = ph;
    setKinematics(l);
    legs.push_back(l);
    return static_cast<int>(legs.size() - 1);
}

struct Row {
    float p1, t1, f1, p2, t2, f2, m, best, tb;
};

struct Event {
    long run = 0, ev = 0;
    int  blk = 0;               // sub-file block inside the run (event number restarts)
    int  trig = -1;             // 0 PT3, 1 PT2, -1 other
    std::vector<Leg> pos, neg;
    std::vector<Row> rows_pm;   // e+e- rows (for same-event spectra)
    int  n_rows_pp = 0, n_rows_mm = 0;
};

// ---------------------------------------------------------------- reader
struct Key {
    long run, ev;
    bool operator<(const Key& o) const {
        return run < o.run || (run == o.run && ev < o.ev);
    }
    bool operator==(const Key& o) const { return run == o.run && ev == o.ev; }
};

struct Reader {
    std::string name;
    TFile* f = nullptr;
    TTree* t = nullptr;
    Long64_t i = 0, n = 0, limit = 0;
    Float_t run = 0, ev = 0, tb = 0, best = 0, m = 0;
    Float_t p1 = 0, t1 = 0, f1 = 0, p2 = 0, t2 = 0, f2 = 0;
    bool   has = false;          // current event valid
    Key    key{0, 0};
    std::vector<Row> rows;       // rows of the current event
    bool   rowLoaded = false;    // row i already in the buffers
    long   run_violations = 0;
    long   trig_mismatch = 0;
    long   last_run = -1;

    bool open(const char* fn, Long64_t max_rows) {
        name = fn;
        f = TFile::Open(fn, "READ");
        if (!f || f->IsZombie()) { std::cerr << "cannot open " << fn << "\n"; return false; }
        t = dynamic_cast<TTree*>(f->Get(kTree));
        if (!t) { std::cerr << "no tree " << kTree << " in " << fn << "\n"; return false; }
        n = t->GetEntries();
        limit = (max_rows >= 0 && max_rows < n) ? max_rows : n;
        t->SetBranchStatus("*", 0);
        const char* br[] = {"runnumber", "event", "trigbit", "isBest", "m_ee",
                            "ep_p_rec", "ep_theta", "ep_phi",
                            "em_p_rec", "em_theta", "em_phi"};
        for (const char* b : br) t->SetBranchStatus(b, 1);
        t->SetBranchAddress("runnumber", &run);
        t->SetBranchAddress("event", &ev);
        t->SetBranchAddress("trigbit", &tb);
        t->SetBranchAddress("isBest", &best);
        t->SetBranchAddress("m_ee", &m);
        t->SetBranchAddress("ep_p_rec", &p1);
        t->SetBranchAddress("ep_theta", &t1);
        t->SetBranchAddress("ep_phi", &f1);
        t->SetBranchAddress("em_p_rec", &p2);
        t->SetBranchAddress("em_theta", &t2);
        t->SetBranchAddress("em_phi", &f2);
        return true;
    }

    // Load the next event (all consecutive rows with equal key, rows < limit).
    bool nextEvent() {
        rows.clear();
        has = false;
        if (i >= limit) return false;
        if (!rowLoaded) { t->GetEntry(i); rowLoaded = true; }
        key = Key{static_cast<long>(run), static_cast<long>(ev)};
        if (last_run >= 0 && key.run < last_run) ++run_violations;
        last_run = key.run;
        const float tb0 = tb;
        while (true) {
            rows.push_back(Row{p1, t1, f1, p2, t2, f2, m, best, tb});
            if (tb != tb0) ++trig_mismatch;
            ++i;
            rowLoaded = false;
            if (i >= limit) break;
            t->GetEntry(i);
            rowLoaded = true;
            if (static_cast<long>(run) != key.run || static_cast<long>(ev) != key.ev) break;
        }
        has = true;
        return true;
    }

    // True if the file was cut by max_rows inside run r (row `limit` still
    // belongs to run r).
    bool truncatedInRun(long r) {
        if (limit >= n) return false;
        t->GetEntry(limit);
        rowLoaded = false;
        return static_cast<long>(run) == r;
    }
};

// ---------------------------------------------------------------- pools
struct PoolEvent {
    std::vector<Leg> pos, neg;
};

inline int multClass(int n) { return n >= 3 ? 3 : n; }
inline int trigClass(float tb) {
    if (tb == 8192.f) return 0;
    if (tb == 4096.f) return 1;
    return -1;
}

}  // namespace mixb

using namespace mixb;

void mixing_build_betanew(Long64_t max_rows = -1,
                          const char* out_file = "plots/output/mixing_betanew.root",
                          int pool_depth = 10,
                          bool verbose = true) {
    TStopwatch sw;
    sw.Start();

    // ---------- readers ----------
    Reader r_em, r_pp, r_mm;
    if (!r_em.open(kFileEm, max_rows)) return;
    if (!r_pp.open(kFilePp, -1)) return;
    if (!r_mm.open(kFileMm, -1)) return;
    std::printf("Input rows: epem=%lld (limit %lld)  epep=%lld  emem=%lld\n",
                r_em.n, r_em.limit, r_pp.n, r_mm.n);
    std::printf("pool_depth=%d  max_rows=%lld\n", pool_depth, max_rows);

    // ---------- spectra ----------
    // index: [trig][scheme 0=all 1=best][kind 0=se 1=mix 2=mixnc][charge 0=pm 1=pp 2=mm]
    Spectrum spec[kNTrig][2][3][3];
    const char* kindName[3] = {"se", "mix", "mixnc"};
    const char* kindTitle[3] = {"same-event ", "mixed (classes) ", "mixed (no classes) "};
    const char* chName[3]  = {"pm", "pp", "mm"};
    const char* chTitle[3] = {"e^{+}e^{-}", "e^{+}e^{+}", "e^{-}e^{-}"};
    for (int t = 0; t < kNTrig; ++t)
        for (int s = 0; s < 2; ++s)
            for (int k = 0; k < 3; ++k)
                for (int c = 0; c < 3; ++c) {
                    std::string nm = std::string("h_") + kindName[k] +
                        "_" + chName[c] + (s == 1 ? "_best" : "") + "_" + trigName(t);
                    std::string ti = std::string(kindTitle[k]) +
                        chTitle[c] + (s == 1 ? " (best)" : " (all)") + " " + trigName(t) +
                        ";M_{ee} [GeV/c^{2}];pairs";
                    spec[t][s][k][c].init(nm, ti);
                }

    // ---------- pools and counters ----------
    std::deque<PoolEvent> pools[kNTrig][kNMultCl][kNMultCl];
    std::deque<PoolEvent> pools_nc[kNTrig];
    long long nev_class[kNTrig][kNMultCl][kNMultCl] = {};
    long long nmix_class[kNTrig][kNMultCl][kNMultCl][3] = {};
    long long nmix_best[kNTrig][3] = {};
    long long nmixnc[kNTrig][3] = {}, nmixnc_best[kNTrig][3] = {};
    long long rows_pm_class[kNTrig][kNMultCl][kNMultCl] = {};
    long long n_events = 0, n_events_trig[kNTrig] = {0, 0}, n_events_othertrig = 0;
    long long n_events_epem = 0, n_events_lsonly[2] = {0, 0};
    long long n_events_11 = 0, n_events_11_bad = 0;
    long long n_rows_pm_ne_npnm = 0;         // e+e- rows != n+ * n-
    long long n_rows_pp_ne_comb = 0;         // e+e+ rows != C(n+,2) (events with n+>=2, n->=1)
    long long n_rows_mm_ne_comb = 0;
    long long n_ls_matched[2] = {0, 0}, n_ls_partial[2] = {0, 0}, n_ls_ambig[2] = {0, 0};
    long long n_ls_legs_added[2] = {0, 0};
    long long n_lsonly_key_clash = 0;        // emem-only event with same ev as an epep-only event
    long long n_blk_resync[2] = {0, 0};      // like-sign block index corrected at a matched event
    long long n_legs_bothcharge = 0, n_trig_mismatch_files = 0;
    long long n_rows_used[3] = {0, 0, 0};
    double    max_dm = 0.0, sum_dm2 = 0.0;
    long long n_dm = 0;
    long      n_runs = 0, n_days = 0, n_days86400 = 0, prev_run = -1, first_run = -1, last_run = -1;
    bool      truncated = false;

    // ---- mixing of one event (legs already complete) ----
    auto mixEvent = [&](Event& E) {
        ++n_events;
        const int np = static_cast<int>(E.pos.size());
        const int nm = static_cast<int>(E.neg.size());
        if (np == 1 && nm == 1) { ++n_events_11; if (E.rows_pm.size() != 1) ++n_events_11_bad; }
        for (const Leg& a : E.pos)
            for (const Leg& b : E.neg)
                if (a.p == b.p && a.th == b.th && a.ph == b.ph) ++n_legs_bothcharge;
        if (E.trig < 0) { ++n_events_othertrig; return; }
        ++n_events_trig[E.trig];
        const int cp = multClass(np), cm = multClass(nm);
        ++nev_class[E.trig][cp][cm];
        rows_pm_class[E.trig][cp][cm] += static_cast<long long>(E.rows_pm.size());
        // pairs of E with every event of a pool; counters n[3] (all), nb[3] (best legs)
        auto mixWith = [&](std::deque<PoolEvent>& pool, Spectrum* sp_all, Spectrum* sp_best,
                           long long* n, long long* nb) {
            for (const PoolEvent& B : pool) {
                for (const Leg& a : E.pos)
                    for (const Leg& b : B.neg) {
                        const double m = massGeV(a, b);
                        sp_all[0].fill(m); ++n[0];
                        if (a.best_pm && b.best_pm) { sp_best[0].fill(m); ++nb[0]; }
                    }
                for (const Leg& a : E.neg)
                    for (const Leg& b : B.pos) {
                        const double m = massGeV(a, b);
                        sp_all[0].fill(m); ++n[0];
                        if (a.best_pm && b.best_pm) { sp_best[0].fill(m); ++nb[0]; }
                    }
                for (const Leg& a : E.pos)
                    for (const Leg& b : B.pos) {
                        const double m = massGeV(a, b);
                        sp_all[1].fill(m); ++n[1];
                        if (a.best_ls && b.best_ls) { sp_best[1].fill(m); ++nb[1]; }
                    }
                for (const Leg& a : E.neg)
                    for (const Leg& b : B.neg) {
                        const double m = massGeV(a, b);
                        sp_all[2].fill(m); ++n[2];
                        if (a.best_ls && b.best_ls) { sp_best[2].fill(m); ++nb[2]; }
                    }
            }
            PoolEvent pe;
            pe.pos = E.pos; pe.neg = E.neg;
            pool.push_back(std::move(pe));
            if (static_cast<int>(pool.size()) > pool_depth) pool.pop_front();
        };
        mixWith(pools[E.trig][cp][cm], spec[E.trig][0][1], spec[E.trig][1][1],
                nmix_class[E.trig][cp][cm], nmix_best[E.trig]);
        mixWith(pools_nc[E.trig], spec[E.trig][0][2], spec[E.trig][1][2],
                nmixnc[E.trig], nmixnc_best[E.trig]);
    };

    // ---- same-event spectra from rows of one file (c = 0 pm, 1 pp, 2 mm) ----
    auto fillSE = [&](const std::vector<Row>& rows, int c) {
        for (const Row& w : rows) {
            const int t = trigClass(w.tb);
            if (t < 0) continue;
            spec[t][0][0][c].fill(w.m);
            if (w.best == 1.f) spec[t][1][0][c].fill(w.m);
        }
    };

    // ---------- run-by-run merge ----------
    r_em.nextEvent();
    r_pp.nextEvent();
    r_mm.nextEvent();
    std::vector<Event> evs, lsonly;
    evs.reserve(8192);
    lsonly.reserve(2048);
    std::vector<int> order;
    order.reserve(8192);
    int max_blocks = 0;
    std::unordered_map<long, std::vector<int>> byEv;
    std::unordered_map<long, int> lsonly_pp_ev;   // ev -> count of epep-only events in this run
    long long n_report = 0;

    while (r_em.has || r_pp.has || r_mm.has) {
        long cur = LONG_MAX;
        if (r_em.has) cur = std::min(cur, r_em.key.run);
        if (r_pp.has) cur = std::min(cur, r_pp.key.run);
        if (r_mm.has) cur = std::min(cur, r_mm.key.run);

        // 1. e+e- events of this run
        evs.clear();
        lsonly.clear();
        byEv.clear();
        lsonly_pp_ev.clear();
        int  em_blk = 0;
        long em_prev_ev = -1;
        while (r_em.has && r_em.key.run == cur) {
            evs.emplace_back();
            Event& E = evs.back();
            E.run = r_em.key.run; E.ev = r_em.key.ev;
            if (em_prev_ev >= 0 && E.ev < em_prev_ev) ++em_blk;
            em_prev_ev = E.ev;
            E.blk = em_blk;
            E.trig = trigClass(r_em.rows.front().tb);
            E.rows_pm = r_em.rows;
            for (const Row& w : r_em.rows) {
                const int ia = addLeg(E.pos, w.p1, w.t1, w.f1);
                const int ib = addLeg(E.neg, w.p2, w.t2, w.f2);
                if (w.best == 1.f) { E.pos[ia].best_pm = true; E.neg[ib].best_pm = true; }
                const double dm = massGeV(E.pos[ia], E.neg[ib]) - static_cast<double>(w.m);
                if (std::fabs(dm) > max_dm) max_dm = std::fabs(dm);
                sum_dm2 += dm * dm; ++n_dm;
            }
            byEv[E.ev].push_back(static_cast<int>(evs.size()) - 1);
            r_em.nextEvent();
        }
        if (!r_em.has && r_em.truncatedInRun(cur)) {
            // run cut by max_rows: drop it and stop (like-sign files too)
            truncated = true;
            break;
        }
        if (max_rows >= 0 && !r_em.has && evs.empty()) break;   // epem exhausted at a run boundary

        // 2. like-sign events of this run: match to e+e- events or stand alone
        for (int j = 0; j < 2; ++j) {
            Reader& r = (j == 0) ? r_pp : r_mm;
            const int c = (j == 0) ? 1 : 2;
            int  ls_blk = 0;
            long ls_prev_ev = -1;
            while (r.has && r.key.run == cur) {
                n_rows_used[c] += static_cast<long long>(r.rows.size());
                fillSE(r.rows, c);
                if (ls_prev_ev >= 0 && r.key.ev < ls_prev_ev) ++ls_blk;
                ls_prev_ev = r.key.ev;
                std::vector<Leg> L;
                for (const Row& w : r.rows) {
                    const int ia = addLeg(L, w.p1, w.t1, w.f1);
                    const int ib = addLeg(L, w.p2, w.t2, w.f2);
                    if (w.best == 1.f) { L[ia].best_ls = true; L[ib].best_ls = true; }
                }
                const int ltrig = trigClass(r.rows.front().tb);
                // candidates with the same event number
                int best_idx = -1, best_nmatch = 0, n_best = 0;
                auto it = byEv.find(r.key.ev);
                if (it != byEv.end()) {
                    for (int idx : it->second) {
                        const std::vector<Leg>& ref = (j == 0) ? evs[idx].pos : evs[idx].neg;
                        int nmatch = 0;
                        for (const Leg& l : L) if (findLeg(ref, l.p, l.th, l.ph) >= 0) ++nmatch;
                        if (nmatch > best_nmatch) { best_nmatch = nmatch; best_idx = idx; n_best = 1; }
                        else if (nmatch > 0 && nmatch == best_nmatch) ++n_best;
                    }
                }
                if (best_idx >= 0) {
                    ++n_ls_matched[j];
                    if (n_best > 1) ++n_ls_ambig[j];
                    if (best_nmatch < static_cast<int>(L.size())) ++n_ls_partial[j];
                    Event& E = evs[best_idx];
                    std::vector<Leg>& ref = (j == 0) ? E.pos : E.neg;
                    for (const Leg& l : L) {
                        const int k = findLeg(ref, l.p, l.th, l.ph);
                        if (k >= 0) { if (l.best_ls) ref[k].best_ls = true; }
                        else { ref.push_back(l); ++n_ls_legs_added[j]; }
                    }
                    if (j == 0) E.n_rows_pp += static_cast<int>(r.rows.size());
                    else        E.n_rows_mm += static_cast<int>(r.rows.size());
                    if (E.trig != ltrig) ++n_trig_mismatch_files;
                    if (E.blk != ls_blk) ++n_blk_resync[j];
                    ls_blk = E.blk;                       // re-synchronise block index
                } else {
                    // like-sign-only event (n- = 0 for epep, n+ = 0 for emem)
                    ++n_events_lsonly[j];
                    lsonly.emplace_back();
                    Event& X = lsonly.back();
                    X.run = cur; X.ev = r.key.ev; X.trig = ltrig; X.blk = ls_blk;
                    if (j == 0) { X.pos = L; X.n_rows_pp = static_cast<int>(r.rows.size()); ++lsonly_pp_ev[X.ev]; }
                    else        { X.neg = L; X.n_rows_mm = static_cast<int>(r.rows.size());
                                  if (lsonly_pp_ev.count(X.ev)) ++n_lsonly_key_clash; }
                }
                r.nextEvent();
            }
        }

        // 3. all events of the run in (block, event) order: e+e- events
        //    (same-event spectra, consistency, mixing) and like-sign-only events
        order.clear();
        for (size_t i = 0; i < evs.size(); ++i)    order.push_back(static_cast<int>(i));
        for (size_t i = 0; i < lsonly.size(); ++i) order.push_back(-1 - static_cast<int>(i));
        std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
            const Event& A = (a >= 0) ? evs[a] : lsonly[-1 - a];
            const Event& B = (b >= 0) ? evs[b] : lsonly[-1 - b];
            return A.blk < B.blk || (A.blk == B.blk && A.ev < B.ev);
        });
        for (int idx : order) {
            if (idx < 0) { mixEvent(lsonly[-1 - idx]); continue; }
            Event& E = evs[idx];
            ++n_events_epem;
            n_rows_used[0] += static_cast<long long>(E.rows_pm.size());
            fillSE(E.rows_pm, 0);
            const int np = static_cast<int>(E.pos.size());
            const int nm = static_cast<int>(E.neg.size());
            if (static_cast<int>(E.rows_pm.size()) != np * nm) ++n_rows_pm_ne_npnm;
            if (E.n_rows_pp != np * (np - 1) / 2) ++n_rows_pp_ne_comb;
            if (E.n_rows_mm != nm * (nm - 1) / 2) ++n_rows_mm_ne_comb;
            mixEvent(E);
        }
        if (em_blk + 1 > max_blocks) max_blocks = em_blk + 1;

        // run / day bookkeeping
        if (prev_run < 0) { first_run = cur; n_days = 1; n_days86400 = 1; }
        else {
            if (cur - prev_run > kDayGap) ++n_days;
            if (cur / kDaySec != prev_run / kDaySec) ++n_days86400;
        }
        ++n_runs;
        prev_run = cur;
        last_run = cur;

        if (verbose && (++n_report % 1000) == 0) {
            std::printf("  ... run %ld (#%ld), %lld events, epem row %lld/%lld, %.0f s\n",
                        cur, n_runs, n_events, r_em.i, r_em.limit, sw.RealTime());
            sw.Continue();
        }
    }

    sw.Stop();
    const double t_real = sw.RealTime(), t_cpu = sw.CpuTime();

    // ---------- finalize ----------
    long long nmix_tot[kNTrig][3] = {};
    for (int t = 0; t < kNTrig; ++t)
        for (int a = 0; a < kNMultCl; ++a)
            for (int b = 0; b < kNMultCl; ++b)
                for (int c = 0; c < 3; ++c) nmix_tot[t][c] += nmix_class[t][a][b][c];

    std::printf("\n========== SUMMARY ==========\n");
    std::printf("rows used: epem=%lld epep=%lld emem=%lld%s\n",
                n_rows_used[0], n_rows_used[1], n_rows_used[2],
                truncated ? "  (last run cut by max_rows dropped)" : "");
    std::printf("runs: %ld  (first %ld, last %ld)  days by gap>%ld: %ld   days by run/86400: %ld\n",
                n_runs, first_run, last_run, kDayGap, n_days, n_days86400);
    std::printf("events: %lld  (e+e-: %lld, e+e+ only: %lld, e-e- only: %lld)\n",
                n_events, n_events_epem, n_events_lsonly[0], n_events_lsonly[1]);
    std::printf("events by trigger: PT3 %lld, PT2 %lld, other %lld\n",
                n_events_trig[0], n_events_trig[1], n_events_othertrig);
    std::printf("run-order violations: epem=%ld epep=%ld emem=%ld\n",
                r_em.run_violations, r_pp.run_violations, r_mm.run_violations);
    std::printf("trigbit mismatch within event rows: epem=%ld epep=%ld emem=%ld; across files=%lld\n",
                r_em.trig_mismatch, r_pp.trig_mismatch, r_mm.trig_mismatch, n_trig_mismatch_files);
    std::printf("like-sign events matched to e+e- events: epep %lld (partial %lld, ambiguous %lld, legs added %lld), "
                "emem %lld (partial %lld, ambiguous %lld, legs added %lld)\n",
                n_ls_matched[0], n_ls_partial[0], n_ls_ambig[0], n_ls_legs_added[0],
                n_ls_matched[1], n_ls_partial[1], n_ls_ambig[1], n_ls_legs_added[1]);
    std::printf("emem-only events sharing (run,event) with an epep-only event: %lld\n", n_lsonly_key_clash);
    std::printf("blocks per run (max): %d; like-sign block index re-synchronised at matched events: epep %lld, emem %lld\n",
                max_blocks, n_blk_resync[0], n_blk_resync[1]);
    std::printf("sanity(1): events with n+=n-=1: %lld, of which e+e- rows != 1: %lld\n",
                n_events_11, n_events_11_bad);
    std::printf("           e+e- events with rows != n+*n-: %lld; e+e+ rows != C(n+,2): %lld; e-e- rows != C(n-,2): %lld\n",
                n_rows_pm_ne_npnm, n_rows_pp_ne_comb, n_rows_mm_ne_comb);
    std::printf("sanity(2): mass formula vs m_ee on %lld e+e- rows: max|dm|=%.3e GeV, rms=%.3e GeV\n",
                n_dm, max_dm, n_dm > 0 ? std::sqrt(sum_dm2 / n_dm) : 0.0);
    std::printf("legs with identical (p,theta,phi) as both e+ and e-: %lld\n", n_legs_bothcharge);
    for (int t = 0; t < kNTrig; ++t) {
        std::printf("\n--- %s: events per class n+ (rows) x n- (cols) [e+e- rows in class] / mixed pairs +-, ++, --\n",
                    trigName(t));
        for (int a = 0; a < kNMultCl; ++a) {
            for (int b = 0; b < kNMultCl; ++b)
                std::printf("  n+=%s n-=%s: %10lld ev [%10lld rows]  mix +-=%11lld ++=%10lld --=%10lld\n",
                            a == 3 ? ">=3" : Form("%d", a), b == 3 ? ">=3" : Form("%d", b),
                            nev_class[t][a][b], rows_pm_class[t][a][b],
                            nmix_class[t][a][b][0], nmix_class[t][a][b][1], nmix_class[t][a][b][2]);
        }
        std::printf("  same-event pairs %s: all +-=%lld ++=%lld --=%lld | best +-=%lld ++=%lld --=%lld\n",
                    trigName(t), spec[t][0][0][0].n, spec[t][0][0][1].n, spec[t][0][0][2].n,
                    spec[t][1][0][0].n, spec[t][1][0][1].n, spec[t][1][0][2].n);
        std::printf("  mixed pairs %s (classes): all +-=%lld ++=%lld --=%lld | best legs +-=%lld ++=%lld --=%lld\n",
                    trigName(t), nmix_tot[t][0], nmix_tot[t][1], nmix_tot[t][2],
                    nmix_best[t][0], nmix_best[t][1], nmix_best[t][2]);
        std::printf("  mixed pairs %s (no classes): all +-=%lld ++=%lld --=%lld | best legs +-=%lld ++=%lld --=%lld\n",
                    trigName(t), nmixnc[t][0], nmixnc[t][1], nmixnc[t][2],
                    nmixnc_best[t][0], nmixnc_best[t][1], nmixnc_best[t][2]);
    }
    std::printf("\ntime: real %.1f s, cpu %.1f s\n", t_real, t_cpu);

    // ---------- write ----------
    const TString dir = gSystem->DirName(out_file);
    gSystem->mkdir(dir, true);
    TFile* fo = TFile::Open(out_file, "RECREATE");
    if (!fo || fo->IsZombie()) { std::cerr << "cannot create " << out_file << "\n"; return; }
    for (int t = 0; t < kNTrig; ++t)
        for (int s = 0; s < 2; ++s)
            for (int k = 0; k < 3; ++k)
                for (int c = 0; c < 3; ++c) {
                    spec[t][s][k][c].finalize();
                    spec[t][s][k][c].write();
                }
    auto wl = [&](const char* nm, long long v) { TParameter<Long64_t>(nm, v).Write(); };
    auto wd = [&](const char* nm, double v) { TParameter<double>(nm, v).Write(); };
    wl("pool_depth", pool_depth);
    wl("max_rows", max_rows);
    wl("truncated", truncated ? 1 : 0);
    wl("n_events", n_events);
    wl("n_events_epem", n_events_epem);
    wl("n_events_epep_only", n_events_lsonly[0]);
    wl("n_events_emem_only", n_events_lsonly[1]);
    wl("n_events_PT3", n_events_trig[0]);
    wl("n_events_PT2", n_events_trig[1]);
    wl("n_events_othertrig", n_events_othertrig);
    wl("n_runs", n_runs);
    wl("n_days_gap", n_days);
    wl("n_days_86400", n_days86400);
    wl("first_run", first_run);
    wl("last_run", last_run);
    wl("rows_epem", n_rows_used[0]);
    wl("rows_epep", n_rows_used[1]);
    wl("rows_emem", n_rows_used[2]);
    wl("run_order_violations", r_em.run_violations + r_pp.run_violations + r_mm.run_violations);
    wl("n_events_11", n_events_11);
    wl("n_events_11_bad", n_events_11_bad);
    wl("n_rows_pm_ne_npnm", n_rows_pm_ne_npnm);
    wl("n_rows_pp_ne_comb", n_rows_pp_ne_comb);
    wl("n_rows_mm_ne_comb", n_rows_mm_ne_comb);
    wl("n_ls_matched_pp", n_ls_matched[0]);
    wl("n_ls_matched_mm", n_ls_matched[1]);
    wl("n_ls_partial_pp", n_ls_partial[0]);
    wl("n_ls_partial_mm", n_ls_partial[1]);
    wl("n_ls_ambig_pp", n_ls_ambig[0]);
    wl("n_ls_ambig_mm", n_ls_ambig[1]);
    wl("n_ls_legs_added_pp", n_ls_legs_added[0]);
    wl("n_ls_legs_added_mm", n_ls_legs_added[1]);
    wl("n_lsonly_key_clash", n_lsonly_key_clash);
    wl("n_blk_resync_pp", n_blk_resync[0]);
    wl("n_blk_resync_mm", n_blk_resync[1]);
    wl("max_blocks_per_run", max_blocks);
    wl("n_legs_bothcharge", n_legs_bothcharge);
    wl("n_trig_mismatch_files", n_trig_mismatch_files);
    wd("max_dm_GeV", max_dm);
    wd("time_real_s", t_real);
    for (int t = 0; t < kNTrig; ++t) {
        for (int c = 0; c < 3; ++c) {
            wl(Form("nmix_%s_%s", chName[c], trigName(t)), nmix_tot[t][c]);
            wl(Form("nmix_%s_best_%s", chName[c], trigName(t)), nmix_best[t][c]);
            wl(Form("nmixnc_%s_%s", chName[c], trigName(t)), nmixnc[t][c]);
            wl(Form("nmixnc_%s_best_%s", chName[c], trigName(t)), nmixnc_best[t][c]);
        }
        for (int a = 0; a < kNMultCl; ++a)
            for (int b = 0; b < kNMultCl; ++b) {
                wl(Form("nev_%s_np%d_nm%d", trigName(t), a, b), nev_class[t][a][b]);
                for (int c = 0; c < 3; ++c)
                    wl(Form("nmix_%s_%s_np%d_nm%d", chName[c], trigName(t), a, b),
                       nmix_class[t][a][b][c]);
            }
    }
    TNamed("info",
           "mixing_build_betanew: run-by-run merge of the three files; like-sign events matched to e+e- "
           "events by event number + leg kinematics; classes trigger x n+ x n- (0,1,2,>=3); rolling pool "
           "of pool_depth events per class; legs unique by (p_rec,theta,phi); mass = PParticle REC formula; "
           "best legs: best_pm for +-, best_ls for ++/--; h_mixnc_*: pool per trigger only (no multiplicity classes)").Write();
    fo->Close();
    std::printf("Written: %s\n", out_file);
}
