// wp2_patcheck.C -- WP2: check that the PAT binary built from the tag prod-feb22-betanew-2026-10-28 (used for the
// empty-target runs) reproduces the production rows: every PipPimEpEm_ID / EpEm_ID row of one full-target DST
// (be2203405213001) processed with that binary must exist, with identical values, in the production file lep_034.root.
// Match key: runnumber, event, isBest, trigbit and the momenta of the tracks (float equality).
#include <cstdio>
#include <set>
#include <tuple>
#include <vector>
#include "TFile.h"
#include "TTree.h"
typedef std::tuple<int, int, float, float, float, float, float, float> Key;
static std::vector<Key> Read(const char* fn, const char* tn, bool fourBody, int runSel) {
    std::vector<Key> v;
    TFile f(fn);
    TTree* t = (TTree*)f.Get(tn);
    if (!t) { printf("no %s in %s\n", tn, fn); return v; }
    Float_t run, ev, isb, tb, vz, a, b, c = 0, d = 0;
    t->SetBranchStatus("*", 0);
    const char* br[] = {"runnumber", "event", "isBest", "trigbit", "eVertReco_z"};
    for (auto x : br) t->SetBranchStatus(x, 1);
    t->SetBranchAddress("runnumber", &run); t->SetBranchAddress("event", &ev); t->SetBranchAddress("isBest", &isb);
    t->SetBranchAddress("trigbit", &tb); t->SetBranchAddress("eVertReco_z", &vz);
    const char* pa = "ep_p"; const char* pb = "em_p";
    t->SetBranchStatus(pa, 1); t->SetBranchStatus(pb, 1); t->SetBranchAddress(pa, &a); t->SetBranchAddress(pb, &b);
    if (fourBody) { t->SetBranchStatus("pip_p", 1); t->SetBranchStatus("pim_p", 1); t->SetBranchAddress("pip_p", &c); t->SetBranchAddress("pim_p", &d); }
    const Long64_t n = t->GetEntries();
    for (Long64_t i = 0; i < n; ++i) {
        t->GetEntry(i);
        if (runSel > 0 && (int)run != runSel) continue;
        v.emplace_back((int)run, (int)ev, isb * 10 + tb / 8192.f, vz, a, b, c, d);
    }
    return v;
}
void wp2_patcheck() {
    const char* mine = "/hdd3/przygoda/HADES/FAT_out/all_report/wp2/full_test/be2203405213001.hld_dst_feb22_out.root";
    const char* prod = "/hdd3/przygoda/HADES/files/exp/gen4/lep_034.root";
    const char* trees[2] = {"PipPimEpEm_ID", "EpEm_ID"};
    for (int k = 0; k < 2; ++k) {
        std::vector<Key> m = Read(mine, trees[k], k == 0, -1);
        if (m.empty()) continue;
        const int run = std::get<0>(m[0]);
        std::vector<Key> p = Read(prod, trees[k], k == 0, run);
        std::multiset<Key> ps(p.begin(), p.end());
        long found = 0;
        for (auto& x : m) { auto it = ps.find(x); if (it != ps.end()) { ++found; ps.erase(it); } }
        printf("%-14s run %d: rows of the test DST %zu, production rows of the run %zu, test rows found identical %ld (%.4f)\n",
               trees[k], run, m.size(), p.size(), found, m.empty() ? 0. : (double)found / m.size());
    }
}
