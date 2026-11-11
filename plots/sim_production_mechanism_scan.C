// sim_production_mechanism_scan.C
//
// Scan output_epem_sim.root / dilepton_nt and tabulate unique combinations
// of the production / decay / vertex sim fields.
//
// Phase-2 of the sim-mechanism research:
//   - Joint tabulation across the 12 ID/info fields per pair
//   - Marginal distributions: ep_sim_id, em_sim_id, ep parent lineage,
//     ep_sim_geninfo1, ep_parent == em_parent breakdown
//   - Same-parent parent ID histogram (likely π0/η/ω/φ decays)
//
// Run:  root -l -b -q 'plots/sim_production_mechanism_scan.C'
//
// NOTE: The dilepton_nt does NOT store separate ep_sim_genweight /
//       em_sim_genweight branches — only a single sim_genweight.
//       The "ep_sim_geninfo" / "em_sim_geninfo" base fields don't exist
//       either — only geninfo1 and geninfo2.  Phase-1 schema was approximate;
//       this macro uses the actual branches and adds processid (which IS
//       present) to the joint key for completeness.

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TStopwatch.h>

#include <cmath>
#include <cstdio>
#include <cstdint>
#include <algorithm>
#include <map>
#include <vector>
#include <string>

namespace {

// Key for the joint tabulation.  All ID-like fields stored as int (cast from
// the Float_t branches; sentinel -1 preserved).
struct JointKey {
    int ep_id, ep_pid, ep_gpid, ep_procid;
    int em_id, em_pid, em_gpid, em_procid;
    int ep_gi1, ep_gi2;
    int em_gi1, em_gi2;

    bool operator<(const JointKey &o) const {
        if (ep_id    != o.ep_id)    return ep_id    < o.ep_id;
        if (ep_pid   != o.ep_pid)   return ep_pid   < o.ep_pid;
        if (ep_gpid  != o.ep_gpid)  return ep_gpid  < o.ep_gpid;
        if (ep_procid!= o.ep_procid)return ep_procid< o.ep_procid;
        if (em_id    != o.em_id)    return em_id    < o.em_id;
        if (em_pid   != o.em_pid)   return em_pid   < o.em_pid;
        if (em_gpid  != o.em_gpid)  return em_gpid  < o.em_gpid;
        if (em_procid!= o.em_procid)return em_procid< o.em_procid;
        if (ep_gi1   != o.ep_gi1)   return ep_gi1   < o.ep_gi1;
        if (ep_gi2   != o.ep_gi2)   return ep_gi2   < o.ep_gi2;
        if (em_gi1   != o.em_gi1)   return em_gi1   < o.em_gi1;
        return em_gi2 < o.em_gi2;
    }
};

struct Counter {
    uint64_t count = 0;
    double   wsum  = 0.0;
};

static inline int f2i(float x) { return (int)std::lround((double)x); }

}  // namespace

void sim_production_mechanism_scan() {
    const char *fname = "output_epem_sim.root";
    const Long64_t MAX_SCAN = 5000000;  // cap per instructions

    TStopwatch sw; sw.Start();

    TFile *f = TFile::Open(fname, "READ");
    if (!f || f->IsZombie()) {
        printf("ERROR: cannot open %s\n", fname);
        return;
    }
    TTree *t = (TTree*)f->Get("dilepton_nt");
    if (!t) {
        printf("ERROR: no dilepton_nt in %s\n", fname);
        return;
    }
    const Long64_t total = t->GetEntries();
    const Long64_t toScan = std::min<Long64_t>(total, MAX_SCAN);

    printf("Opened %s : dilepton_nt entries = %lld ; scanning %lld\n",
           fname, (long long)total, (long long)toScan);

    // Use SetBranchAddress for speed; disable all branches first then enable
    // the 13 we need.
    t->SetBranchStatus("*", 0);
    const char *branches[] = {
        "ep_sim_id", "ep_sim_parentid", "ep_sim_grandparentid",
        "ep_sim_processid", "ep_sim_geninfo1", "ep_sim_geninfo2",
        "em_sim_id", "em_sim_parentid", "em_sim_grandparentid",
        "em_sim_processid", "em_sim_geninfo1", "em_sim_geninfo2",
        "sim_genweight",
    };
    for (const char *b : branches) t->SetBranchStatus(b, 1);

    Float_t ep_id, ep_pid, ep_gpid, ep_procid, ep_gi1, ep_gi2;
    Float_t em_id, em_pid, em_gpid, em_procid, em_gi1, em_gi2;
    Float_t w;

    t->SetBranchAddress("ep_sim_id",            &ep_id);
    t->SetBranchAddress("ep_sim_parentid",      &ep_pid);
    t->SetBranchAddress("ep_sim_grandparentid", &ep_gpid);
    t->SetBranchAddress("ep_sim_processid",     &ep_procid);
    t->SetBranchAddress("ep_sim_geninfo1",      &ep_gi1);
    t->SetBranchAddress("ep_sim_geninfo2",      &ep_gi2);
    t->SetBranchAddress("em_sim_id",            &em_id);
    t->SetBranchAddress("em_sim_parentid",      &em_pid);
    t->SetBranchAddress("em_sim_grandparentid", &em_gpid);
    t->SetBranchAddress("em_sim_processid",     &em_procid);
    t->SetBranchAddress("em_sim_geninfo1",      &em_gi1);
    t->SetBranchAddress("em_sim_geninfo2",      &em_gi2);
    t->SetBranchAddress("sim_genweight",        &w);

    // Aggregators
    std::map<JointKey, Counter> joint;

    std::map<int, Counter> mep_id;              // (A) ep_sim_id marginal
    std::map<int, Counter> mem_id;              // (A) em_sim_id marginal
    std::map<std::pair<int,int>, Counter> mep_lineage;  // (B) (ep_id, ep_pid)
    std::map<std::pair<int,int>, Counter> mem_lineage;  // (B) (em_id, em_pid)
    std::map<int, Counter> mep_gi1;             // (C) ep_geninfo1 marginal
    std::map<int, Counter> mem_gi1;             // (C) em_geninfo1 marginal

    uint64_t sameParent_n = 0;
    uint64_t diffParent_n = 0;
    double   sameParent_w = 0.0;
    double   diffParent_w = 0.0;

    std::map<int, Counter> sameParentID;        // (E) when same parent, what parent ID

    // Loop
    for (Long64_t i = 0; i < toScan; ++i) {
        t->GetEntry(i);

        const int e_id    = f2i(ep_id);
        const int e_pid   = f2i(ep_pid);
        const int e_gpid  = f2i(ep_gpid);
        const int e_proc  = f2i(ep_procid);
        const int e_gi1   = f2i(ep_gi1);
        const int e_gi2   = f2i(ep_gi2);
        const int m_id    = f2i(em_id);
        const int m_pid   = f2i(em_pid);
        const int m_gpid  = f2i(em_gpid);
        const int m_proc  = f2i(em_procid);
        const int m_gi1   = f2i(em_gi1);
        const int m_gi2   = f2i(em_gi2);
        const double ww   = (double)w;

        JointKey k{e_id,e_pid,e_gpid,e_proc, m_id,m_pid,m_gpid,m_proc,
                   e_gi1,e_gi2, m_gi1,m_gi2};
        auto &c = joint[k];
        c.count++; c.wsum += ww;

        auto &a1 = mep_id[e_id]; a1.count++; a1.wsum += ww;
        auto &a2 = mem_id[m_id]; a2.count++; a2.wsum += ww;
        auto &b1 = mep_lineage[{e_id,e_pid}]; b1.count++; b1.wsum += ww;
        auto &b2 = mem_lineage[{m_id,m_pid}]; b2.count++; b2.wsum += ww;
        auto &g1 = mep_gi1[e_gi1]; g1.count++; g1.wsum += ww;
        auto &g2 = mem_gi1[m_gi1]; g2.count++; g2.wsum += ww;

        if (e_pid == m_pid) {
            sameParent_n++; sameParent_w += ww;
            auto &p = sameParentID[e_pid]; p.count++; p.wsum += ww;
        } else {
            diffParent_n++; diffParent_w += ww;
        }

        if ((i+1) % 1000000 == 0) {
            printf("  ... scanned %lld / %lld  (unique joint = %zu)\n",
                   (long long)(i+1), (long long)toScan, joint.size());
        }
    }

    sw.Stop();

    // ===============================================================
    // REPORT
    // ===============================================================
    const uint64_t N = (uint64_t)toScan;

    auto pct = [&](uint64_t c) { return 100.0 * (double)c / (double)N; };

    printf("\n");
    printf("========================================================================\n");
    printf("SIM PRODUCTION MECHANISM SCAN — output_epem_sim.root / dilepton_nt\n");
    printf("========================================================================\n");
    printf("Total entries in tree : %lld\n", (long long)total);
    printf("Entries scanned       : %lld  (%.2f%% of tree)\n",
           (long long)N, 100.0 * (double)N / (double)total);
    printf("Wall time             : %.1f s real / %.1f s cpu\n",
           sw.RealTime(), sw.CpuTime());
    printf("Unique joint combos   : %zu\n", joint.size());
    printf("\n");

    // ---- JOINT TABLE (sorted by descending count) ----------------------
    std::vector<std::pair<JointKey, Counter>> sortedJoint(joint.begin(), joint.end());
    std::sort(sortedJoint.begin(), sortedJoint.end(),
              [](const auto &a, const auto &b){ return a.second.count > b.second.count; });

    printf("------------------------------------------------------------------------\n");
    printf("JOINT TOP COMBINATIONS  (12 fields, sorted by entry count)\n");
    printf("Columns: rank | count | %%cum | %%row | <wbar> | ep(id,pid,gpid,proc,gi1,gi2) | em(id,pid,gpid,proc,gi1,gi2)\n");
    printf("------------------------------------------------------------------------\n");

    uint64_t cum = 0;
    int rank = 0;
    bool reached95 = false;
    const int MAX_PRINT = 50;
    for (auto &kv : sortedJoint) {
        rank++;
        cum += kv.second.count;
        const double cumPct = 100.0 * (double)cum / (double)N;
        const double rowPct = pct(kv.second.count);
        const double wbar = (kv.second.count > 0) ? kv.second.wsum / (double)kv.second.count : 0.0;
        const JointKey &k = kv.first;
        printf("  %3d | %10llu | %6.2f | %6.2f | %.5f | ep(%d,%d,%d,%d,%d,%d) | em(%d,%d,%d,%d,%d,%d)\n",
               rank, (unsigned long long)kv.second.count, cumPct, rowPct, wbar,
               k.ep_id,k.ep_pid,k.ep_gpid,k.ep_procid,k.ep_gi1,k.ep_gi2,
               k.em_id,k.em_pid,k.em_gpid,k.em_procid,k.em_gi1,k.em_gi2);

        if (!reached95 && cumPct >= 95.0) {
            reached95 = true;
            printf("  --- cumulative coverage reached 95.0%% at rank %d ---\n", rank);
            if (rank >= MAX_PRINT) break;
        }
        if (rank >= MAX_PRINT) break;
    }
    if (!reached95) {
        printf("  (95%% coverage not reached within top %d rows)\n", MAX_PRINT);
    }
    printf("\n");

    // ---- (A) Marginal ep_sim_id, em_sim_id ----------------------------
    auto dumpIntMarginal = [&](const char *label, const std::map<int, Counter> &m) {
        std::vector<std::pair<int,Counter>> v(m.begin(), m.end());
        std::sort(v.begin(), v.end(),
                  [](const auto &a, const auto &b){ return a.second.count > b.second.count; });
        printf("------------------------------------------------------------------------\n");
        printf("MARGINAL: %s   (unique values = %zu)\n", label, v.size());
        printf("  value | count | %%row | <wbar>\n");
        for (auto &kv : v) {
            const double wbar = (kv.second.count > 0) ? kv.second.wsum / (double)kv.second.count : 0.0;
            printf("  %6d | %10llu | %6.2f | %.5f\n",
                   kv.first, (unsigned long long)kv.second.count,
                   pct(kv.second.count), wbar);
        }
        printf("\n");
    };
    dumpIntMarginal("ep_sim_id", mep_id);
    dumpIntMarginal("em_sim_id", mem_id);

    // ---- (B) Lineage marginals ----------------------------------------
    auto dumpLineage = [&](const char *label, const std::map<std::pair<int,int>, Counter> &m) {
        std::vector<std::pair<std::pair<int,int>,Counter>> v(m.begin(), m.end());
        std::sort(v.begin(), v.end(),
                  [](const auto &a, const auto &b){ return a.second.count > b.second.count; });
        printf("------------------------------------------------------------------------\n");
        printf("MARGINAL LINEAGE: %s   (unique (id,pid) pairs = %zu)\n", label, v.size());
        printf("  rank | (id,pid) | count | %%row | <wbar>\n");
        const int LIMIT = 30;
        int r = 0;
        for (auto &kv : v) {
            r++;
            const double wbar = (kv.second.count > 0) ? kv.second.wsum / (double)kv.second.count : 0.0;
            printf("  %4d | (%4d,%4d) | %10llu | %6.2f | %.5f\n",
                   r, kv.first.first, kv.first.second,
                   (unsigned long long)kv.second.count,
                   pct(kv.second.count), wbar);
            if (r >= LIMIT) { printf("  ... (truncated at %d)\n", LIMIT); break; }
        }
        printf("\n");
    };
    dumpLineage("(ep_sim_id, ep_sim_parentid)", mep_lineage);
    dumpLineage("(em_sim_id, em_sim_parentid)", mem_lineage);

    // ---- (C) ep_sim_geninfo1 / em_sim_geninfo1 marginals --------------
    dumpIntMarginal("ep_sim_geninfo1", mep_gi1);
    dumpIntMarginal("em_sim_geninfo1", mem_gi1);

    // ---- (D) same parent vs different parent -------------------------
    printf("------------------------------------------------------------------------\n");
    printf("DICHOTOMY: ep_sim_parentid == em_sim_parentid ?\n");
    printf("  SAME parent : %12llu  (%6.2f%%)   sum_w = %.4f\n",
           (unsigned long long)sameParent_n, pct(sameParent_n), sameParent_w);
    printf("  DIFF parent : %12llu  (%6.2f%%)   sum_w = %.4f\n",
           (unsigned long long)diffParent_n, pct(diffParent_n), diffParent_w);
    printf("\n");

    // ---- (E) parent ID histogram for same-parent subset ---------------
    std::vector<std::pair<int,Counter>> spv(sameParentID.begin(), sameParentID.end());
    std::sort(spv.begin(), spv.end(),
              [](const auto &a, const auto &b){ return a.second.count > b.second.count; });
    printf("------------------------------------------------------------------------\n");
    printf("SAME-PARENT parent ID distribution   (unique parent IDs = %zu)\n", spv.size());
    printf("  rank | parent_id | count | %%of-same-parent | %%of-total | <wbar>\n");
    const int LIMIT_E = 30;
    int r = 0;
    for (auto &kv : spv) {
        r++;
        const double pctSP = (sameParent_n > 0)
            ? 100.0 * (double)kv.second.count / (double)sameParent_n : 0.0;
        const double wbar = (kv.second.count > 0) ? kv.second.wsum / (double)kv.second.count : 0.0;
        printf("  %4d | %9d | %10llu | %14.2f | %9.2f | %.5f\n",
               r, kv.first,
               (unsigned long long)kv.second.count,
               pctSP, pct(kv.second.count), wbar);
        if (r >= LIMIT_E) { printf("  ... (truncated at %d)\n", LIMIT_E); break; }
    }
    printf("\n");

    printf("========================================================================\n");
    printf("DONE.  Unique joint combinations = %zu\n", joint.size());
    printf("========================================================================\n");

    f->Close();
    delete f;
}
