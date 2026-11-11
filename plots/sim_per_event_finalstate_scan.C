// sim_per_event_finalstate_scan.C
//
// Scan SMASH "default" simulation file and tabulate, per event:
//   1) Final-state topology (from candidate ntuples grouped by event index)
//   2) Resonance decay channels (decoded from sim_geninfo / sim_geninfo1)
//   3) Particle multiplicity histograms
//
// The input file (smash_00.root) contains 3 filtered candidate TNtuples:
//      EpEm_ID         (e+e-)
//      EpEmEpEm_ID     (e+e-e+e-)
//      PipPimEpEm_ID   (pi+ pi- e+ e-)
//
// Each ROW = one reconstructed candidate, with per-particle branches
// (em_sim_id, ep_sim_id, pip_sim_id, pim_sim_id, em_sim_geninfo, ...).
// Rows are grouped per simulated event via the `event` branch.
//
// Generator-truth decay encoding (HADES Pluto / HGeantKine convention):
//     sim_id        = HGeant particle ID of the reconstructed track
//                     (1=gamma, 2=e+, 3=e-, 7=pi0, 8=pi+, 9=pi-, 14=p, 17=eta,
//                      41=rho0, 52=omega, 55=phi, ...)
//     sim_parentid  = -1  -> primary final-state particle
//     sim_processid = 0   -> from generator (not detector secondary)
//     sim_geninfo1  = parent encoding: <parentHGeantID>*1000 + <decayCode>
//                     e.g. 7051   -> pi0 (7)   Dalitz (051)
//                          17051  -> eta (17)  Dalitz
//                          52051  -> omega(52) Dalitz
//                          7,17,41,52,55 (no decay-suffix) -> direct 2-body
//     sim_geninfo2  = generation depth (1 = primary chain step)
//
// Usage:  root -l -b -q 'plots/sim_per_event_finalstate_scan.C'
//
// Tunables: kMaxEvents (events to scan), kInputFile.

#include <TFile.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TString.h>
#include <TLeaf.h>
#include <TROOT.h>

#include <algorithm>
#include <cstdio>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace {

// ---------------- HGeant particle ID -> short name --------------------------
const char* hgeantName(int id) {
    switch (id) {
        case 1:  return "gamma";
        case 2:  return "ep";
        case 3:  return "em";
        case 4:  return "nu";
        case 5:  return "mup";
        case 6:  return "mum";
        case 7:  return "pi0";
        case 8:  return "pip";
        case 9:  return "pim";
        case 10: return "K0L";
        case 11: return "Kp";
        case 12: return "Km";
        case 13: return "n";
        case 14: return "p";
        case 15: return "pbar";
        case 16: return "K0S";
        case 17: return "eta";
        case 18: return "Lambda";
        case 19: return "SigmaP";
        case 20: return "Sigma0";
        case 21: return "SigmaM";
        case 22: return "Xi0";
        case 23: return "XiM";
        case 24: return "Omega";
        case 25: return "nbar";
        case 26: return "LambdaBar";
        case 27: return "SigmaPbar";
        case 28: return "Sigma0bar";
        case 29: return "SigmaMbar";
        case 30: return "Xi0bar";
        case 31: return "XiMbar";
        case 32: return "Omegabar";
        case 41: return "rho0";
        case 42: return "rhoP";
        case 43: return "rhoM";
        case 44: return "Delta++";
        case 45: return "Delta+";
        case 46: return "Delta-";
        case 47: return "Lambda(1405)";
        case 48: return "JPsi";
        case 49: return "psi(2S)";
        case 50: return "Delta+_aux";
        case 51: return "Delta0";
        case 52: return "omega";
        case 53: return "eta_prime";
        case 54: return "etaC";
        case 55: return "phi";
        case 56: return "N1440";
        case 57: return "N1520";
        case 58: return "N1535";
        default: return nullptr;
    }
}
std::string nameOrTag(int id) {
    const char* n = hgeantName(id);
    if (n) return std::string(n);
    char buf[32];
    std::snprintf(buf, sizeof(buf), "[id=%d]", id);
    return std::string(buf);
}

// Canonical sort key for a particle name in a topology string
int sortPriority(const std::string& n) {
    static const std::map<std::string,int> P = {
        {"p",0}, {"n",1}, {"pbar",2}, {"nbar",3},
        {"pip",10}, {"pim",11}, {"pi0",12},
        {"Kp",20}, {"Km",21}, {"K0L",22}, {"K0S",23},
        {"eta",30}, {"eta_prime",31},
        {"rho0",40}, {"rhoP",41}, {"rhoM",42},
        {"omega",50}, {"phi",51},
        {"Lambda",60}, {"SigmaP",61}, {"Sigma0",62}, {"SigmaM",63},
        {"Delta++",70}, {"Delta+",71}, {"Delta0",72}, {"Delta-",73},
        {"ep",80}, {"em",81}, {"mup",82}, {"mum",83}, {"gamma",84},
    };
    auto it = P.find(n);
    if (it != P.end()) return it->second;
    return 999;
}
bool nameLess(const std::string& a, const std::string& b) {
    int pa = sortPriority(a), pb = sortPriority(b);
    if (pa != pb) return pa < pb;
    return a < b;
}

std::string canonicalTopology(std::vector<std::string> names) {
    std::sort(names.begin(), names.end(), nameLess);
    std::string out;
    for (size_t i = 0; i < names.size(); ++i) {
        if (i) out += ' ';
        out += names[i];
    }
    return out;
}

// Decay-code (last 3 digits of geninfo / geninfo1) -> human label.
// Pluto/HADES decay codes. These are the most common ones for dilepton sources.
std::string decayCodeLabel(int parentId, int code) {
    // Generic dilepton decay codes
    if (code == 51 || code == 51) return "e+ e- gamma";  // Dalitz
    if (code == 41) return "e+ e-";                       // direct dielectron
    // For specific resonances, sometimes the 2-body decay is implicit (no suffix)
    switch (parentId) {
        case 7:  // pi0
            if (code == 1)  return "gamma gamma";
            if (code == 51) return "e+ e- gamma";   // Dalitz
            if (code == 7)  return "e+ e-";         // direct (rare)
            break;
        case 17: // eta
            if (code == 1)  return "gamma gamma";
            if (code == 2)  return "pi0 pi0 pi0";
            if (code == 3)  return "pi+ pi- pi0";
            if (code == 5)  return "pi+ pi- gamma";
            if (code == 51) return "e+ e- gamma";   // Dalitz
            break;
        case 41: // rho0
            if (code == 0 || code == 41) return "e+ e-";
            if (code == 1) return "pi+ pi-";
            break;
        case 52: // omega
            if (code == 1)  return "pi+ pi- pi0";
            if (code == 2)  return "pi0 gamma";
            if (code == 3)  return "pi+ pi-";
            if (code == 4)  return "eta gamma";
            if (code == 51) return "e+ e- pi0";     // omega Dalitz
            if (code == 7)  return "e+ e-";         // direct
            break;
        case 55: // phi
            if (code == 1) return "Kp Km";
            if (code == 2) return "K0L K0S";
            if (code == 3) return "rho pi / 3pi";
            if (code == 7) return "e+ e-";
            break;
        case 53: // eta'
            if (code == 1) return "pi+ pi- eta";
            if (code == 2) return "rho0 gamma";
            break;
        case 47: // Lambda(1405)
            if (code == 1) return "Sigma pi";
            break;
        default: break;
    }
    char buf[64];
    std::snprintf(buf, sizeof(buf), "code=%d", code);
    return std::string(buf);
}

// Decode geninfo1 = parent*1000 + decayCode, or geninfo1 = parent (no code)
struct ParentDecay { int parentId; int decayCode; bool hasCode; };
ParentDecay decodeGenInfo1(int g1) {
    // Filter SMASH ID flag bit
    if (g1 >= 1000000) {
        // SMASH-encoded parent id (1000000 + smashPid). Map a few known.
        int smashId = g1 - 1000000;
        // 1000223 = omega, 1000333 = phi, 1000113 = rho0, 1000221 = eta, 1000331 = eta',
        // 1000213 = rho+, ...  Just route them to HGeant equivalents if obvious.
        int hg = -smashId; // fall-through marker
        if (smashId == 1) hg = 7;       // hack: not strictly real
        return {hg, -1, false};
    }
    if (g1 == 0) return {0, -1, false};
    if (g1 < 1000) {
        // Just the parent HGeant id, no decay-code suffix
        return {g1, -1, false};
    }
    int code   = g1 % 1000;
    int parent = g1 / 1000;
    return {parent, code, true};
}

// ---------------- Per-particle "primary?" tag ----------------------------
bool isPrimary(int parentid, int processid) {
    return parentid == -1 && processid == 0;
}

// ---------------- Generic per-ntuple processor ------------------------------
//
// For each row, we collect the candidate particle list (em, ep, pip, pim) and
// their sim_geninfo1 -> parent-decay info. Rows are grouped by `event`.
//
struct ParticleHit {
    int   pid;           // HGeant id
    int   genParentId;   // parent HGeant id (from geninfo1)
    int   genDecayCode;  // decay code (from geninfo1), -1 if absent
    int   gen2;          // generation depth (geninfo2)
};

struct PerEventBucket {
    std::vector<ParticleHit> particles;
};

void addParticle(PerEventBucket& bucket,
                 int pid, int parentid, int processid,
                 int geninfo1, int geninfo2) {
    if (!isPrimary(parentid, processid)) return;
    auto d = decodeGenInfo1(geninfo1);
    ParticleHit h;
    h.pid          = pid;
    h.genParentId  = d.parentId;
    h.genDecayCode = d.decayCode;
    h.gen2         = geninfo2;
    bucket.particles.push_back(h);
}

// Get value from a TNtuple row by branch name; returns 0 if branch missing.
inline float getF(TNtuple* t, const char* name, std::map<std::string, TLeaf*>& cache) {
    auto it = cache.find(name);
    TLeaf* lf = nullptr;
    if (it == cache.end()) {
        lf = t->GetLeaf(name);
        cache[name] = lf;
    } else lf = it->second;
    if (!lf) return 0.f;
    return lf->GetValue();
}

}  // anonymous namespace

void sim_per_event_finalstate_scan() {
    // ---------------------------- Config ----------------------------
    const char* kInputFile  = "/hdd1/HADES/pp45/sim/smash/default/smash_00.root";
    const Long64_t kMaxRowsPerNtuple = 2000000;   // limit per ntuple for speed
    const int   kTopTopologies = 40;
    const int   kTopChannels   = 12;

    // ---------------------------- Open file ----------------------------
    TFile* f = TFile::Open(kInputFile);
    if (!f || f->IsZombie()) {
        std::printf("ERROR: cannot open %s\n", kInputFile);
        return;
    }
    std::printf("File: %s\n", kInputFile);

    // ---------------------------- Per-event aggregator ----------------------------
    // Group particles by event index (from `event` branch) across all 3 ntuples.
    // Note: this is "all candidates per event"; same event can appear with multiple
    // (em,ep) combinations because each row = one reconstructed pair candidate.
    // We DEDUPLICATE per-event particles by (pid, geninfo1, geninfo2) to avoid
    // counting the same generator-level particle multiple times.
    //
    std::map<int, std::set<std::tuple<int,int,int,int>>> evDedup; // event -> set of (pid, genParent, decay, gen2)
    std::map<int, PerEventBucket> events;

    auto pushIfNew = [&](int evid, const ParticleHit& h) {
        auto key = std::make_tuple(h.pid, h.genParentId, h.genDecayCode, h.gen2);
        auto& s = evDedup[evid];
        if (s.insert(key).second) events[evid].particles.push_back(h);
    };

    // ---------- Process EpEm_ID ----------
    {
        TNtuple* t = (TNtuple*) f->Get("EpEm_ID");
        if (!t) {
            std::printf("WARN: EpEm_ID not found\n");
        } else {
            Long64_t N = std::min<Long64_t>(t->GetEntries(), kMaxRowsPerNtuple);
            std::printf("Processing EpEm_ID:        %lld / %lld rows\n",
                        N, t->GetEntries());
            std::map<std::string, TLeaf*> cache;
            for (Long64_t i = 0; i < N; ++i) {
                t->GetEntry(i);
                int evid = (int) getF(t, "event", cache);
                int em_id   = (int) getF(t, "em_sim_id", cache);
                int em_par  = (int) getF(t, "em_sim_parentid", cache);
                int em_proc = (int) getF(t, "em_sim_processid", cache);
                int em_g1   = (int) getF(t, "em_sim_geninfo1", cache);
                int em_g2   = (int) getF(t, "em_sim_geninfo2", cache);
                int ep_id   = (int) getF(t, "ep_sim_id", cache);
                int ep_par  = (int) getF(t, "ep_sim_parentid", cache);
                int ep_proc = (int) getF(t, "ep_sim_processid", cache);
                int ep_g1   = (int) getF(t, "ep_sim_geninfo1", cache);
                int ep_g2   = (int) getF(t, "ep_sim_geninfo2", cache);
                if (isPrimary(em_par, em_proc)) {
                    auto d = decodeGenInfo1(em_g1);
                    pushIfNew(evid, {em_id, d.parentId, d.decayCode, em_g2});
                }
                if (isPrimary(ep_par, ep_proc)) {
                    auto d = decodeGenInfo1(ep_g1);
                    pushIfNew(evid, {ep_id, d.parentId, d.decayCode, ep_g2});
                }
            }
        }
    }

    // ---------- Process PipPimEpEm_ID ----------
    {
        TNtuple* t = (TNtuple*) f->Get("PipPimEpEm_ID");
        if (!t) {
            std::printf("WARN: PipPimEpEm_ID not found\n");
        } else {
            Long64_t N = std::min<Long64_t>(t->GetEntries(), kMaxRowsPerNtuple);
            std::printf("Processing PipPimEpEm_ID:  %lld / %lld rows\n",
                        N, t->GetEntries());
            std::map<std::string, TLeaf*> cache;
            const char* prefixes[4] = {"em","ep","pip","pim"};
            for (Long64_t i = 0; i < N; ++i) {
                t->GetEntry(i);
                int evid = (int) getF(t, "event", cache);
                for (auto* pr : prefixes) {
                    std::string id_n   = std::string(pr) + "_sim_id";
                    std::string par_n  = std::string(pr) + "_sim_parentid";
                    std::string proc_n = std::string(pr) + "_sim_processid";
                    std::string g1_n   = std::string(pr) + "_sim_geninfo1";
                    std::string g2_n   = std::string(pr) + "_sim_geninfo2";
                    int pid   = (int) getF(t, id_n.c_str(), cache);
                    int par   = (int) getF(t, par_n.c_str(), cache);
                    int proc  = (int) getF(t, proc_n.c_str(), cache);
                    int g1    = (int) getF(t, g1_n.c_str(), cache);
                    int g2    = (int) getF(t, g2_n.c_str(), cache);
                    if (isPrimary(par, proc)) {
                        auto d = decodeGenInfo1(g1);
                        pushIfNew(evid, {pid, d.parentId, d.decayCode, g2});
                    }
                }
            }
        }
    }

    // ---------- Process EpEmEpEm_ID (rare four-lepton) ----------
    {
        TNtuple* t = (TNtuple*) f->Get("EpEmEpEm_ID");
        if (!t) {
            std::printf("WARN: EpEmEpEm_ID not found\n");
        } else {
            Long64_t N = std::min<Long64_t>(t->GetEntries(), kMaxRowsPerNtuple);
            std::printf("Processing EpEmEpEm_ID:    %lld / %lld rows\n",
                        N, t->GetEntries());
            std::map<std::string, TLeaf*> cache;
            // The branch naming convention for EpEmEpEm is em/ep/em2/ep2 in some
            // formats; probe what's present. We look for any branch matching
            // <pref>_sim_id.
            std::vector<std::string> prefs;
            TObjArray* br = t->GetListOfBranches();
            for (int k = 0; k < br->GetEntries(); ++k) {
                TString s = br->At(k)->GetName();
                if (s.EndsWith("_sim_id") && !s.Contains("_sim_id_")) {
                    prefs.push_back(std::string(s.Data(), s.Length() - 7));
                }
            }
            for (Long64_t i = 0; i < N; ++i) {
                t->GetEntry(i);
                int evid = (int) getF(t, "event", cache);
                for (auto& pr : prefs) {
                    int pid  = (int) getF(t, (pr + "_sim_id").c_str(), cache);
                    int par  = (int) getF(t, (pr + "_sim_parentid").c_str(), cache);
                    int proc = (int) getF(t, (pr + "_sim_processid").c_str(), cache);
                    int g1   = (int) getF(t, (pr + "_sim_geninfo1").c_str(), cache);
                    int g2   = (int) getF(t, (pr + "_sim_geninfo2").c_str(), cache);
                    if (isPrimary(par, proc)) {
                        auto d = decodeGenInfo1(g1);
                        pushIfNew(evid, {pid, d.parentId, d.decayCode, g2});
                    }
                }
            }
        }
    }

    std::printf("\nTotal distinct events seen: %zu\n", events.size());

    // ---------------------------- Build statistics ----------------------------
    std::map<std::string, Long64_t> topologyCounts;          // canon topology -> count
    std::map<int, std::map<std::string, Long64_t>> decayChannels; // parentId -> (label -> count)
    std::map<int, Long64_t>           parentCounts;
    std::map<std::string, std::map<int, Long64_t>> multHist;  // name -> mult -> count

    for (auto& kv : events) {
        const PerEventBucket& bkt = kv.second;

        // -- Topology --
        std::vector<std::string> names;
        std::map<std::string,int> mult;
        for (auto& h : bkt.particles) {
            std::string nm = nameOrTag(h.pid);
            names.push_back(nm);
            mult[nm]++;
        }
        if (!names.empty()) {
            topologyCounts[canonicalTopology(names)]++;
        }

        // -- Multiplicities --
        for (auto& m : mult) multHist[m.first][m.second]++;
        // also track absence (0) for common species
        for (auto* nm : {"pi0","pip","pim","eta","gamma","p"}) {
            if (mult.find(nm) == mult.end()) multHist[nm][0]++;
        }

        // -- Decay channels: group same-event particles by (parentId, decayCode)
        // so two e from same eta -> one "eta -> e+e-gamma" entry
        std::map<std::tuple<int,int>, std::vector<int>> parentToDaughters;
        for (auto& h : bkt.particles) {
            if (h.genParentId <= 0) continue;
            parentToDaughters[{h.genParentId, h.genDecayCode}].push_back(h.pid);
        }
        for (auto& pd : parentToDaughters) {
            int parent = std::get<0>(pd.first);
            int code   = std::get<1>(pd.first);
            parentCounts[parent]++;
            std::string label;
            if (code >= 0) {
                label = decayCodeLabel(parent, code);
            } else {
                // unknown / not encoded
                std::vector<std::string> dnames;
                for (int d : pd.second) dnames.push_back(nameOrTag(d));
                std::sort(dnames.begin(), dnames.end(), nameLess);
                label = "obs:";
                for (auto& dn : dnames) { label += ' '; label += dn; }
            }
            decayChannels[parent][label]++;
        }
    }

    // ---------------------------- TABLE 1: topology counts ----------------------------
    std::printf("\n");
    std::printf("=================================================================\n");
    std::printf("TABLE 1  Per-event final-state topology (candidate-level)\n");
    std::printf("         Note: built from filtered candidate ntuples — only events\n");
    std::printf("         with at least one reconstructed e+/e- (and optional pi+/pi-)\n");
    std::printf("         candidate appear; full pp final state NOT available in this file.\n");
    std::printf("=================================================================\n");
    std::vector<std::pair<std::string,Long64_t>> topoSorted(topologyCounts.begin(), topologyCounts.end());
    std::sort(topoSorted.begin(), topoSorted.end(),
              [](const auto& a, const auto& b){ return a.second > b.second; });
    Long64_t totalEvents = 0;
    for (auto& kv : topoSorted) totalEvents += kv.second;
    std::printf("  %-50s %12s  %7s\n", "topology", "count", "%");
    std::printf("  ------------------------------------------------------------------------\n");
    int shown = 0;
    for (auto& kv : topoSorted) {
        if (shown >= kTopTopologies) break;
        std::printf("  %-50s %12lld  %6.2f%%\n",
                    kv.first.c_str(), kv.second,
                    100.0 * kv.second / std::max<Long64_t>(1, totalEvents));
        ++shown;
    }
    std::printf("  ------------------------------------------------------------------------\n");
    std::printf("  Total unique topologies: %zu\n", topoSorted.size());
    std::printf("  Total events tabulated : %lld\n", totalEvents);

    // ---------------------------- TABLE 2: decay channels per resonance ----------------------------
    std::printf("\n");
    std::printf("=================================================================\n");
    std::printf("TABLE 2  Resonance decay channels (decoded from sim_geninfo1)\n");
    std::printf("=================================================================\n");

    // Sort parents by count desc
    std::vector<std::pair<int,Long64_t>> parents(parentCounts.begin(), parentCounts.end());
    std::sort(parents.begin(), parents.end(),
              [](const auto& a, const auto& b){ return a.second > b.second; });

    for (auto& p : parents) {
        int parent = p.first;
        Long64_t pn = p.second;
        std::string pname = nameOrTag(parent);
        std::printf("\n  %s  (parentId=%d, total decays seen: %lld)\n",
                    pname.c_str(), parent, pn);
        std::printf("    %-30s %12s  %7s\n", "channel", "count", "%");
        std::printf("    -----------------------------------------------------------\n");
        auto& cm = decayChannels[parent];
        std::vector<std::pair<std::string,Long64_t>> chans(cm.begin(), cm.end());
        std::sort(chans.begin(), chans.end(),
                  [](const auto& a, const auto& b){ return a.second > b.second; });
        int csh = 0;
        for (auto& c : chans) {
            if (csh >= kTopChannels) break;
            std::printf("    %-30s %12lld  %6.2f%%\n",
                        c.first.c_str(), c.second,
                        100.0 * c.second / std::max<Long64_t>(1, pn));
            ++csh;
        }
    }

    // ---------------------------- TABLE 3: multiplicity hists ----------------------------
    std::printf("\n");
    std::printf("=================================================================\n");
    std::printf("TABLE 3  Per-event particle multiplicities (candidates only)\n");
    std::printf("=================================================================\n");
    std::vector<std::string> mkeys = {"ep","em","pip","pim","pi0","eta","gamma","p"};
    for (auto& nm : mkeys) {
        auto it = multHist.find(nm);
        if (it == multHist.end()) continue;
        std::printf("\n  %s multiplicity:\n", nm.c_str());
        for (auto& mb : it->second) {
            std::printf("    n=%d : %lld\n", mb.first, mb.second);
        }
    }

    f->Close();
    std::printf("\nDone.\n");
}
