// ============================================================================
// sim_eta_checks.C  --  checks of the SMASH eta set (sim_set = 3) on scan_nt
// ============================================================================
// Reads scan_nt of the FAT outputs directly (not the compact file):
//   /hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/output_eta_{pippimepem,pippimepep,pippimemem}.root
//   /hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/output_std_pippimepem.root   (closure only)
//
// part bit 1  inventory and recomputation (tables/eta_validate.txt):
//             rows, event groups (adjacent rows with equal (event, eVertReco_z)),
//             isBest == 1 per group, vz > -160, owned rows and sum w (double)
//             per four-particle class and template; every row recomputed from
//             its stored truth columns with src/sim_truth.h EvaluateRow (set 3,
//             PDG BR, BR(eta -> pi+pi-e+e-) = 2.68e-4, BR(X17) = 0, as the
//             production config) and compared with the stored sim_* columns;
//             X17 weight: class-6 rows recomputed with BR(X17) = 1e-5
// part bit 2  float-code check (tables/eta_codes.txt): every stored geninfo,
//             geninfo1, geninfo2 value of the four slots of every row mapped
//             back with simtruth::CodeCandidates; unknown and ambiguous values
//             listed with their rows; flag bits 3 / 4 compared row by row
// part bit 4  normalisation closure and eta production (tables/eta_closure.txt):
//             (i)  eta-set rows with no track of a forced eta decay (natural
//                  rows, not owned) against std rows with no eta product, after
//                  "pre", per M4 region; rows and sum w; event-group errors
//             (ii) std rows with both pions from one eta (pion class 2) / BR
//                  (0.2114 + 0.0478, HGeant) against eta-set rows with both
//                  pions from the forced eta (pion class 9) / 1
//
// Run from /home/przygoda/HADES/FAT/FAT:
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/aclic_eta_checks", kTRUE)' 'plots/eta_sim/sim_eta_checks.C+(7)'
// ============================================================================
#include <cmath>
#include <cstdio>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include "TFile.h"
#include "TStopwatch.h"
#include "TTree.h"

#include "../../src/sim_truth.h"

namespace etachk {

const char* const kDir = "/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew";
const char* const kTab = "/home/przygoda/HADES/FAT/FAT/plots/eta_sim/tables";

#define ECHK_SLOT_COLS(X, T)                                                        \
    X(T##_sim_id) X(T##_sim_parentid) X(T##_sim_grandparentid) X(T##_sim_processid) \
    X(T##_sim_geninfo) X(T##_sim_geninfo1) X(T##_sim_geninfo2) X(T##_sim_genweight) \
    X(T##_sim_vertexx) X(T##_sim_vertexy) X(T##_sim_vertexz)                       \
    X(T##_sim_px) X(T##_sim_py) X(T##_sim_pz) X(T##_sim_isghost)

#define ECHK_COLS(X)                                                                 \
    X(event) X(eVertReco_z) X(isBest) X(trigbit) X(cor_ok) X(m_pippimepem_cor)       \
    X(sim_set) X(sim_cls) X(sim_lcls) X(sim_pcls) X(sim_tmpl) X(sim_w) X(sim_w_db)   \
    X(sim_own_A) X(sim_own_B) X(sim_flags) X(sim_rules_ver)                          \
    X(m4_true) X(mee_true) X(mpipi_true)                                             \
    ECHK_SLOT_COLS(X, ep) ECHK_SLOT_COLS(X, em) ECHK_SLOT_COLS(X, pip) ECHK_SLOT_COLS(X, pim)

struct In {
#define X(n) Float_t n;
    ECHK_COLS(X)
#undef X
};

inline TTree* Open(const std::string& fn, TFile*& f, In& in) {
    f = TFile::Open(fn.c_str(), "READ");
    if (!f || f->IsZombie()) { printf("cannot open %s\n", fn.c_str()); return nullptr; }
    TTree* t = (TTree*)f->Get("scan_nt");
    if (!t) { printf("scan_nt missing in %s\n", fn.c_str()); return nullptr; }
    t->SetBranchStatus("*", 0);
#define X(n) t->SetBranchStatus(#n, 1); t->SetBranchAddress(#n, &in.n);
    ECHK_COLS(X)
#undef X
    t->SetCacheSize(200 * 1024 * 1024);
    return t;
}

inline void Tracks(const In& in, simtruth::Trk t[4]) {
#define TRK(i, T)                                                                                       \
    t[i] = simtruth::Trk{static_cast<int>(in.T##_sim_id), static_cast<int>(in.T##_sim_parentid),        \
                         static_cast<int>(in.T##_sim_grandparentid), static_cast<int>(in.T##_sim_processid), \
                         in.T##_sim_geninfo, in.T##_sim_geninfo1, in.T##_sim_geninfo2, in.T##_sim_genweight, \
                         in.T##_sim_vertexx, in.T##_sim_vertexy, in.T##_sim_vertexz,                     \
                         in.T##_sim_px, in.T##_sim_py, in.T##_sim_pz, static_cast<int>(in.T##_sim_isghost)};
    TRK(0, ep) TRK(1, em) TRK(2, pip) TRK(3, pim)
#undef TRK
}

// "pre" on scan_nt rows: isBest >= 0, z > -500, PT3 (bit test), vz > -160, corrected momenta
inline bool Pre(const In& in) {
    return in.isBest >= 0 && in.eVertReco_z > -500.f && in.eVertReco_z > -160.f && in.cor_ok == 1.f &&
           ((int)in.trigbit & 8192) != 0;
}
inline bool ForcedEta(const simtruth::Trk& x) { return x.g0 == simtruth::kEtaPPee || x.g0 == simtruth::kEtaX17; }
// a track related to an eta of the generator: product of a SMASH final-state eta (geninfo 1717, inherited by
// every descendant), parent or grandparent eta (GEANT id 17), any product of phi -> gamma eta (550117), or a
// track of a forced eta decay (eta set)
inline bool EtaRelated(const simtruth::Trk& x) {
    return x.g0 == 1717.f || x.g0 == 550117.f || x.par == 17 || x.gp == 17 || ForcedEta(x);
}

// sums with event-group variance (rows counted with weight 1 and with sim_w)
struct GSum {
    double n = 0, vn = 0, w = 0, vw = 0;   // totals and group variances
    double gn = 0, gw = 0;                 // current group
    void add(double ww) { gn += 1; gw += ww; }
    void flush() { n += gn; vn += gn * gn; w += gw; vw += gw * gw; gn = gw = 0; }
};
inline void Ratio(double a, double va, double b, double vb, double& r, double& e) {
    r = b > 0 ? a / b : 0;
    e = (a > 0 && b > 0) ? r * std::sqrt(va / (a * a) + vb / (b * b)) : 0;
}

}  // namespace etachk

void sim_eta_checks(int parts = 7) {
    using namespace etachk;
    TStopwatch swAll;
    const char* hypTag[3] = {"pippimepem", "pippimepep", "pippimemem"};
    const int q1[3] = {2, 2, 3}, q2[3] = {3, 2, 3};
    simtruth::Opt opt;          // PDG BR, BR(eta -> pi+pi-e+e-) 2.68e-4, BR(X17) 0: the production config
    simtruth::Opt optX = opt;   // X17 illustration
    optX.brX17 = 1e-5;

    // ======================= parts 1 and 2: one pass over the three eta outputs =======================
    if (parts & 3) {
        FILE* fv = (parts & 1) ? fopen((std::string(kTab) + "/eta_validate.txt").c_str(), "w") : nullptr;
        FILE* fc = (parts & 2) ? fopen((std::string(kTab) + "/eta_codes.txt").c_str(), "w") : nullptr;
        auto PV = [&](const std::string& s) { fputs(s.c_str(), stdout); if (fv) fputs(s.c_str(), fv); };
        auto PC = [&](const std::string& s) { fputs(s.c_str(), stdout); if (fc) fputs(s.c_str(), fc); };
        if (parts & 1) {
            PV("# eta_validate: SMASH eta set (sim_set = 3), scan_nt of the three FAT outputs\n");
            PV("# recomputation: src/sim_truth.h EvaluateRow(kEta, tracks, q1, q2, Opt{pdg, BR(eta->pi+pi-ee) 2.68e-4, BR(X17) 0})\n");
            PV("# 'pre' = isBest >= 0, z > -500, PT3 bit, vz > -160, cor_ok == 1; owned = sim_own_A (== sim_own_B in the eta set)\n");
        }
        // code statistics over all rows of all three files: column -> stored value -> number of track slots
        std::map<float, long long> codeCount[3];
        std::map<float, long long> codeRows[3];   // rows having that value in any slot
        std::map<float, long long> wForced, wOther;
        long long rowsUnk = 0, rowsAmb = 0, rowsFlagUnk = 0, rowsFlagAmb = 0, mismatchUnk = 0, mismatchAmb = 0;
        for (int h = 0; h < 3; ++h) {
            TStopwatch sw;
            TFile* f = nullptr;
            In in;
            TTree* t = Open(std::string(kDir) + "/output_eta_" + hypTag[h] + ".root", f, in);
            if (!t) return;
            const Long64_t n = t->GetEntries();
            long long groups = 0, best1 = 0, grpBad = 0, nBest = 0, rowsVz = 0, rowsPre = 0, set3 = 0, corBad = 0;
            long long ownA = 0, ownAB = 0, nanW = 0, ownVz = 0, notOwnW0 = 0, notOwnForced = 0, ownNoForced = 0;
            double sumOwnVz = 0, sumOwnPre = 0, sumWdbOwnVz = 0;
            std::map<int, std::pair<long long, double>> byCls, byTmpl, byClsTmpl;   // owned, vz > -160
            long long mis[12] = {0};   // mismatches: cls lcls pcls tmpl flags ownA ownB w wdb m4 mee mpipi
            double maxRelW = 0;
            long long x17rows = 0, x17bad = 0, x17other = 0;
            double x17otherW = 0;
            float keyE = -1, keyZ = -1e9;
            for (Long64_t i = 0; i < n; ++i) {
                t->GetEntry(i);
                if (i == 0 || in.event != keyE || in.eVertReco_z != keyZ) {
                    if (i > 0 && nBest != 1) grpBad++;
                    ++groups; nBest = 0; keyE = in.event; keyZ = in.eVertReco_z;
                }
                if (in.isBest == 1.f) { ++best1; ++nBest; }
                if (in.sim_set == 3.f) ++set3;
                if (in.cor_ok != 1.f) ++corBad;
                if (in.sim_w != in.sim_w || in.sim_w_db != in.sim_w_db) ++nanW;
                simtruth::Trk tr[4];
                Tracks(in, tr);
                bool anyForced = false;
                for (int s = 0; s < 4; ++s) if (ForcedEta(tr[s])) anyForced = true;
                const bool own = in.sim_own_A == 1.f;
                if (own) ++ownA;
                if (own && in.sim_own_B == 1.f) ++ownAB;
                if (!own && anyForced) ++notOwnForced;
                if (own && !anyForced) ++ownNoForced;
                if (!own && in.sim_w == 0.f) ++notOwnW0;
                if (in.eVertReco_z > -160.f) {
                    ++rowsVz;
                    if (own) {
                        ++ownVz; sumOwnVz += (double)in.sim_w; sumWdbOwnVz += (double)in.sim_w_db;
                        auto& a = byCls[(int)in.sim_cls]; a.first++; a.second += (double)in.sim_w;
                        auto& b = byTmpl[(int)in.sim_tmpl]; b.first++; b.second += (double)in.sim_w;
                        auto& c = byClsTmpl[10 * (int)in.sim_cls + (int)in.sim_tmpl]; c.first++; c.second += (double)in.sim_w;
                    }
                }
                if (Pre(in)) { ++rowsPre; if (own) sumOwnPre += (double)in.sim_w; }

                // ---- recomputation ----
                if (parts & 1) {
                    const simtruth::RowTruth r = simtruth::EvaluateRow(simtruth::kEta, tr, q1[h], q2[h], opt);
                    if (r.cls != (int)in.sim_cls) mis[0]++;
                    if (r.lcls != (int)in.sim_lcls) mis[1]++;
                    if (r.pcls != (int)in.sim_pcls) mis[2]++;
                    if (r.tmpl != (int)in.sim_tmpl) mis[3]++;
                    if (r.flags != (int)in.sim_flags) mis[4]++;
                    if ((r.ownA ? 1.f : 0.f) != in.sim_own_A) mis[5]++;
                    if ((r.ownB ? 1.f : 0.f) != in.sim_own_B) mis[6]++;
                    const double dw = std::fabs(r.w - (double)in.sim_w), dwdb = std::fabs(r.wdb - (double)in.sim_w_db);
                    if (dw > 1e-6 * std::fabs(r.w) + 1e-30) mis[7]++;
                    if (dwdb > 1e-6 * std::fabs(r.wdb) + 1e-30) mis[8]++;
                    if (r.w > 0) maxRelW = std::max(maxRelW, dw / r.w);
                    if (std::fabs(r.m4 - in.m4_true) > 1e-5) mis[9]++;
                    if (std::fabs(r.mee - in.mee_true) > 1e-5) mis[10]++;
                    if (std::fabs(r.mpipi - in.mpipi_true) > 1e-5) mis[11]++;
                    // X17: class 6 rows with BR(X17) = 1e-5 must have W = 2e-5 exactly (one decay, all four tracks forced)
                    const double wx = simtruth::RowWeight(simtruth::kEta, tr, optX);
                    if (r.cls == 6 && r.ownA) {
                        ++x17rows;
                        if (std::fabs(wx - 2e-5) > 1e-12) ++x17bad;
                    } else if (r.ownA && in.eVertReco_z > -160.f && wx != r.w) {
                        ++x17other; x17otherW += wx - r.w;   // rows that would gain weight from X17 outside class 6
                    }
                }
                // ---- float codes ----
                if (parts & 2) {
                    bool unk = false, amb = false;
                    simtruth::Code cc[4];
                    for (int s = 0; s < 4; ++s) {
                        const float g[3] = {tr[s].g0, tr[s].g1, tr[s].g2};
                        for (int k = 0; k < 3; ++k) {
                            codeCount[k][g[k]]++;
                            const int nc = simtruth::CodeCandidates(k, g[k], cc, 4);
                            if (nc == 0) unk = true;
                            if (nc > 1) amb = true;
                        }
                        if (ForcedEta(tr[s])) wForced[tr[s].w]++; else wOther[tr[s].w]++;
                    }
                    std::set<float> seen[3];
                    for (int s = 0; s < 4; ++s) { seen[0].insert(tr[s].g0); seen[1].insert(tr[s].g1); seen[2].insert(tr[s].g2); }
                    for (int k = 0; k < 3; ++k) for (float v : seen[k]) codeRows[k][v]++;
                    const int fl = (int)in.sim_flags;
                    const bool fu = (fl & simtruth::kFlagUnknownCode) != 0, fa = (fl & simtruth::kFlagAmbiguousCode) != 0;
                    if (unk) ++rowsUnk;
                    if (amb) ++rowsAmb;
                    if (fu) ++rowsFlagUnk;
                    if (fa) ++rowsFlagAmb;
                    if (unk != fu) ++mismatchUnk;
                    if (amb != fa) ++mismatchAmb;
                }
            }
            if (n > 0 && nBest != 1) grpBad++;
            f->Close();
            delete f;
            if (parts & 1) {
                PV(Form("\n[eta %s] %lld rows (sim_set == 3: %lld), %lld event groups, isBest==1 rows %lld, groups without exactly one isBest==1 row %lld, cor_ok != 1 rows %lld, NaN w %lld\n",
                        hypTag[h], n, set3, groups, best1, grpBad, corBad, nanW));
                PV(Form("  vz > -160: %lld rows; owned (own_A) %lld (all vz), own_A and own_B %lld; owned with vz > -160 %lld, sum w %.8g, sum w_db %.8g\n",
                        rowsVz, ownA, ownAB, ownVz, sumOwnVz, sumWdbOwnVz));
                PV(Form("  'pre': %lld rows, sum w of owned rows %.8g\n", rowsPre, sumOwnPre));
                PV(Form("  ownership against the codes: owned without a forced-eta track %lld; not owned with a forced-eta track %lld; not owned with w == 0 %lld (of %lld not owned)\n",
                        ownNoForced, notOwnForced, notOwnW0, n - ownA));
                PV("  owned, vz > -160, per four-particle class: rows, sum w\n");
                for (auto& x : byCls) PV(Form("    class %d: %10lld  %.8g\n", x.first, x.second.first, x.second.second));
                PV("  owned, vz > -160, per template: rows, sum w\n");
                for (auto& x : byTmpl) PV(Form("    T%d: %10lld  %.8g\n", x.first, x.second.first, x.second.second));
                PV("  owned, vz > -160, per (class, template): rows, sum w\n");
                for (auto& x : byClsTmpl) PV(Form("    class %d T%d: %10lld  %.8g\n", x.first / 10, x.first % 10, x.second.first, x.second.second));
                PV(Form("  recomputation mismatches (cls lcls pcls tmpl flags ownA ownB w w_db m4 mee mpipi): %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld ; max rel |w - stored| %.2g\n",
                        mis[0], mis[1], mis[2], mis[3], mis[4], mis[5], mis[6], mis[7], mis[8], mis[9], mis[10], mis[11], maxRelW));
                PV(Form("  X17: owned class-6 rows %lld, of them with W(BR 1e-5) != 2e-5: %lld ; owned rows (vz > -160) outside class 6 whose W changes with BR(X17) = 1e-5: %lld, sum of the change %.6g\n",
                        x17rows, x17bad, x17other, x17otherW));
                PV(Form("  %.1f s\n", sw.RealTime()));
            }
        }
        if (parts & 2) {
            PC("# eta_codes: stored generator codes of the eta set (all rows of the three outputs, four slots, three columns)\n");
            PC("# mapping by src/sim_truth.h CodeCandidates: every integer that rounds to the stored float is tested against the dictionary\n");
            const char* cn[3] = {"geninfo", "geninfo1", "geninfo2"};
            for (int k = 0; k < 3; ++k) {
                long long nUni = 0, nUnk = 0, nAmb = 0, sUni = 0, sUnk = 0, sAmb = 0;
                std::vector<std::string> unkList, ambList, bigList;
                for (auto& x : codeCount[k]) {
                    simtruth::Code cc[8];
                    const int nc = simtruth::CodeCandidates(k, x.first, cc, 8);
                    if (nc == 0) { nUnk++; sUnk += x.second; unkList.push_back(Form("%.10g (%lld slots, %lld rows)", (double)x.first, x.second, codeRows[k][x.first])); }
                    else if (nc > 1) {
                        nAmb++; sAmb += x.second;
                        std::string s = Form("%.10g (%lld slots, %lld rows) ->", (double)x.first, x.second, codeRows[k][x.first]);
                        for (int j = 0; j < nc && j < 8; ++j) s += Form(" %lld", cc[j].code);
                        if (nc > 8) s += " ...";
                        ambList.push_back(s);
                    } else {
                        nUni++; sUni += x.second;
                        if (std::fabs(x.first) >= 16777216.f || cc[0].code != (long long)x.first)
                            bigList.push_back(Form("%.10g -> %lld [%s] (%lld slots, %lld rows)", (double)x.first, cc[0].code,
                                                   simtruth::CodeFamilyName(cc[0].family), x.second, codeRows[k][x.first]));
                    }
                }
                PC(Form("\n## column %s: %zu distinct stored values; unique %lld (%lld slots), unknown %lld (%lld slots), ambiguous %lld (%lld slots)\n",
                        cn[k], codeCount[k].size(), nUni, sUni, nUnk, sUnk, nAmb, sAmb));
                PC("  values >= 2^24 with a unique candidate (stored float -> exact code):\n");
                for (auto& s : bigList) PC("    " + s + "\n");
                PC("  unknown values:\n");
                for (size_t j = 0; j < unkList.size() && j < 40; ++j) PC("    " + unkList[j] + "\n");
                if (unkList.size() > 40) PC(Form("    ... %zu more\n", unkList.size() - 40));
                PC("  ambiguous values:\n");
                for (size_t j = 0; j < ambList.size() && j < 40; ++j) PC("    " + ambList[j] + "\n");
                if (ambList.size() > 40) PC(Form("    ... %zu more\n", ambList.size() - 40));
            }
            {
                simtruth::Code cc[8];
                const int n1 = simtruth::CodeCandidates(0, simtruth::kEtaPPee, cc, 8);
                PC(Form("\nforced codes: stored 1702030848 -> %d candidate(s):", n1));
                for (int j = 0; j < n1 && j < 8; ++j) PC(Form(" %lld [%s]", cc[j].code, simtruth::CodeFamilyName(cc[j].family)));
                const int n2 = simtruth::CodeCandidates(0, simtruth::kEtaX17, cc, 8);
                PC(Form("\nforced codes: stored 17700808 -> %d candidate(s):", n2));
                for (int j = 0; j < n2 && j < 8; ++j) PC(Form(" %lld [%s]", cc[j].code, simtruth::CodeFamilyName(cc[j].family)));
                PC("\n");
            }
            PC(Form("rows with an unknown value %lld (flag bit 3 set: %lld; rows where they differ %lld)\n", rowsUnk, rowsFlagUnk, mismatchUnk));
            PC(Form("rows with an ambiguous value %lld (flag bit 4 set: %lld; rows where they differ %lld)\n", rowsAmb, rowsFlagAmb, mismatchAmb));
            PC("stored genweight of tracks with a forced-eta code (value: slots):");
            for (auto& x : wForced) PC(Form(" %g:%lld", (double)x.first, x.second));
            PC("\nstored genweight of the other tracks (value: slots):");
            int c = 0;
            for (auto& x : wOther) { PC(Form(" %g:%lld", (double)x.first, x.second)); if (++c > 12) { PC(" ..."); break; } }
            PC("\n");
        }
        if (fv) fclose(fv);
        if (fc) fclose(fc);
    }

    // ======================= part 4: closure and eta production =======================
    if (parts & 4) {
        FILE* fo = fopen((std::string(kTab) + "/eta_closure.txt").c_str(), "w");
        auto P = [&](const std::string& s) { fputs(s.c_str(), stdout); if (fo) fputs(s.c_str(), fo); };
        const int nreg = 8;
        const double rlo[nreg] = {0.0, 0.30, 0.455, 0.5336, 0.56, 0.70, 1.00, -1e9};
        const double rhi[nreg] = {0.30, 0.455, 0.5225, 0.5542, 0.70, 1.00, 1e9, 1e9};
        const char* rname[nreg] = {"M4 < 0.30", "0.30-0.455", "bump 0.455-0.5225", "peak 0.5336-0.5542", "0.56-0.70",
                                   "0.70-1.00", "M4 > 1.00", "all"};
        // [set: 0 std, 1 eta][region]
        GSum nat[2][nreg];
        // production: [set][0 all rows, 1 leptons not from the same eta, 2 neither lepton related to any eta]
        GSum prod[2][3];
        // identity of natural rows: (event, M4 bits, vz bits) of std natural rows; matches in the eta set
        std::set<std::tuple<float, float, float>> stdKeys;
        long long natMatch = 0;
        long long natEtaRel = 0, natEtaRel1717 = 0, natTot = 0, stdEtaRel = 0, stdTot = 0;
        for (int s = 0; s < 2; ++s) {
            TStopwatch sw;
            TFile* f = nullptr;
            In in;
            const std::string fn = std::string(kDir) + (s == 0 ? "/output_std_pippimepem.root" : "/output_eta_pippimepem.root");
            TTree* t = Open(fn, f, in);
            if (!t) return;
            const Long64_t n = t->GetEntries();
            float keyE = -1, keyZ = -1e9;
            auto flushAll = [&]() {
                for (int r = 0; r < nreg; ++r) nat[s][r].flush();
                prod[s][0].flush(); prod[s][1].flush(); prod[s][2].flush();
            };
            for (Long64_t i = 0; i < n; ++i) {
                t->GetEntry(i);
                if (i == 0 || in.event != keyE || in.eVertReco_z != keyZ) {
                    flushAll();
                    keyE = in.event; keyZ = in.eVertReco_z;
                }
                if (!Pre(in)) continue;
                simtruth::Trk tr[4];
                Tracks(in, tr);
                bool forced = false, rel = false, r1717 = false;
                for (int k = 0; k < 4; ++k) {
                    if (ForcedEta(tr[k])) forced = true;
                    if (EtaRelated(tr[k])) rel = true;
                    if (tr[k].g0 == 1717.f || tr[k].par == 17 || tr[k].gp == 17) r1717 = true;
                }
                const double m4 = in.m_pippimepem_cor;
                const double w = in.sim_w;
                // (i) natural rows
                bool take = false;
                if (s == 0) { stdTot++; if (rel) stdEtaRel++; take = !rel; }
                else if (!forced) {
                    natTot++;
                    if (rel) natEtaRel++;
                    if (r1717) natEtaRel1717++;
                    take = !rel;
                }
                if (take) {
                    for (int r = 0; r < nreg; ++r)
                        if (m4 >= rlo[r] && m4 < rhi[r]) nat[s][r].add(w);
                    const std::tuple<float, float, float> key(in.event, in.m_pippimepem_cor, in.eVertReco_z);
                    if (s == 0) stdKeys.insert(key);
                    else if (stdKeys.count(key)) natMatch++;
                }
                // (ii) eta production
                const int pc = (int)in.sim_pcls, cl = (int)in.sim_cls;
                if ((s == 0 && pc == 2) || (s == 1 && pc == 9)) {
                    prod[s][0].add(w);
                    const bool sameEtaLeptons = s == 0 ? (cl >= 2 && cl <= 4) : (cl == 5 || cl == 6);
                    if (!sameEtaLeptons) prod[s][1].add(w);
                    if (!EtaRelated(tr[0]) && !EtaRelated(tr[1])) prod[s][2].add(w);
                }
            }
            flushAll();
            f->Close();
            delete f;
            printf("closure pass %s: %.1f s\n", s == 0 ? "std" : "eta", sw.RealTime());
        }
        P("# eta_closure: does the eta set cover the same 10^10 SMASH events as the standard set? (e+e- hypothesis, 'pre')\n");
        P("# 'pre' = isBest >= 0, z > -500, PT3 bit, vz > -160, cor_ok == 1. Errors: event groups; the two sets are treated as\n");
        P("# independent (their GEANT simulations are separate; the shared SMASH events make the true error smaller)\n");
        P("# eta-related track: geninfo 1717 (SMASH eta and descendants), 550117 (phi -> gamma eta), parent or grandparent\n");
        P("#   GEANT id 17, or a forced eta code (1702030848, 17700808)\n");
        P(Form("# std rows after pre: %lld, with an eta-related track %lld (%.2f %%)\n", stdTot, stdEtaRel, 100. * stdEtaRel / std::max(1LL, stdTot)));
        P(Form("# eta-set natural rows (no forced-eta track) after pre: %lld; with another eta-related track %lld (of them 1717 / parent 17: %lld)\n",
               natTot, natEtaRel, natEtaRel1717));
        P("\n## (i) natural rows: eta set (no forced-eta track, no eta-related track) / std (no eta-related track)\n");
        P(Form("%-20s %10s %10s %16s %14s %14s %16s\n", "M4 region", "rows std", "rows eta", "ratio rows", "sum w std", "sum w eta", "ratio sum w"));
        for (int r = 0; r < nreg; ++r) {
            double rr, er, rw, ew;
            Ratio(nat[1][r].n, nat[1][r].vn, nat[0][r].n, nat[0][r].vn, rr, er);
            Ratio(nat[1][r].w, nat[1][r].vw, nat[0][r].w, nat[0][r].vw, rw, ew);
            P(Form("%-20s %10.0f %10.0f %7.4f +- %.4f %14.6g %14.6g %7.4f +- %.4f\n", rname[r], nat[0][r].n, nat[1][r].n, rr, er,
                   nat[0][r].w, nat[1][r].w, rw, ew));
        }
        P(Form("# identical rows: %lld of the %.0f eta-set natural rows have the same (event, M4, vz) floats as a std natural row (%zu keys)\n",
               natMatch, nat[1][nreg - 1].n, stdKeys.size()));
        P("\n## (ii) eta production: rows with both pions from one eta, divided by the probability of a pi+pi- decay of that eta\n");
        P("#   std: pion class 2 / (BR_HGeant(eta -> pi+pi-pi0) + BR_HGeant(eta -> pi+pi-gamma)) = / (0.2114 + 0.0478)\n");
        P("#   eta: pion class 9 / 1 (every eta forced into pi+pi-e+e- or X17 pi+pi-, both with a pi+pi- pair)\n");
        P("#   'all' = every lepton pair; 'other leptons' = the lepton pair is not from the same eta (std: four-class not 2-4,\n");
        P("#   eta: four-class not 5-6), i.e. the pions of the eta combined with leptons from elsewhere;\n");
        P("#   'unrelated leptons' = neither lepton track is eta-related (no lepton of this or any other eta)\n");
        const double brStd = simtruth::kHG_etaPPP0 + simtruth::kHG_etaPPG;
        const char* vname[3] = {"all", "other leptons", "unrelated leptons"};
        for (int v = 0; v < 3; ++v) {
            const double a = prod[1][v].n, va = prod[1][v].vn, b = prod[0][v].n / brStd, vb = prod[0][v].vn / (brStd * brStd);
            double r, e;
            Ratio(a, va, b, vb, r, e);
            P(Form("%-18s std rows %8.0f (/%.4f = %10.1f)   eta rows %8.0f   ratio eta/std %.4f +- %.4f\n",
                   vname[v], prod[0][v].n, brStd, b, a, r, e));
        }
        if (fo) fclose(fo);
    }
    printf("sim_eta_checks: total %.1f s\n", swAll.RealTime());
}
