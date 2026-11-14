// ============================================================================
// sim_validate.C  --  validation of the simulation production (pi+pi-e+e-)
// ============================================================================
// Reads scan_nt of the six FAT outputs
//   /hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/output_{std,lep}_{pippimepem,pippimepep,pippimemem}.root
// and checks, per output:
//   1. rows after the FAT selection and after vz > -160 against the expected numbers
//   2. the partition tables (ownership, lepton-pair class, four-particle class,
//      rows the lepton set provides, diagnostics) built from the STORED columns,
//      written in the format of the reference macro partition.C and compared
//      field by field with the full-file reference tables part_*.txt
//      (all sums in double precision)
//   3. one isBest == 1 row per event group (adjacent rows with equal
//      (event, eVertReco_z); runnumber is 0 in simulation)
//   4. no NaN in sim_w, sim_w_db and the corrected masses; cor_ok, trigger bit
//   5. every row re-evaluated with src/sim_truth.h from the stored raw truth
//      columns: classes, ownership, weights, flags, template, true masses must
//      equal the stored columns (tests the alignment of all sim columns)
//   6. M(pi+pi-e+e-) recomputed from the corrected momenta = stored value
//   7. columns: same 233 names in all six outputs, every branch has as many
//      entries as the tree (a late branch would be shifted), the other ntuples
//
// Output: plots/eta_sim/tables/engine_validate_summary.txt and
//         plots/eta_sim/tables/engine_validate_part_{std,lep}_{us,pp,mm}.txt
//
// Run from /home/przygoda/HADES/FAT/FAT:
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/aclic_validate", kTRUE)' 'plots/eta_sim/sim_validate.C+("all")'
// The argument selects the outputs: "all" or a tag such as "std_us", "lep_pp".
// ============================================================================
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "TBranch.h"
#include "TFile.h"
#include "TKey.h"
#include "TLorentzVector.h"
#include "TObjArray.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TTree.h"
#include "TVector3.h"

#include "../../src/sim_truth.h"

namespace sval {

using namespace simtruth;

const char* const kOutDir = "/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew";
const char* const kRefDir =
    "/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/ppsim/spec";
const char* const kTabDir = "/home/przygoda/HADES/FAT/FAT/plots/eta_sim/tables";

struct Def {
    const char* tag;     // std_us ...
    const char* file;    // output file name
    const char* ref;     // reference table
    int set;             // 1 std, 2 lep
    int q1, q2;          // GEANT ids expected in the two lepton slots
    Long64_t expAll;     // rows after the FAT selection (std: implement.md; lep: production.log only)
    Long64_t expVz;      // rows after vz > -160 (part_*.txt)
};
const Def kDefs[6] = {
    {"std_us", "output_std_pippimepem.root", "part_std_us.txt", 1, 2, 3, 230469, 230115},
    {"std_pp", "output_std_pippimepep.root", "part_std_pp.txt", 1, 2, 2, 1133, 1132},
    {"std_mm", "output_std_pippimemem.root", "part_std_mm.txt", 1, 3, 3, 4791, 4777},
    {"lep_us", "output_lep_pippimepem.root", "part_lep_us.txt", 2, 2, 3, 7441490, 7432170},
    {"lep_pp", "output_lep_pippimepep.root", "part_lep_pp.txt", 2, 2, 2, 205527, 205288},
    {"lep_mm", "output_lep_pippimemem.root", "part_lep_mm.txt", 2, 3, 3, 495711, 495101},
};

struct Acc {
    double n = 0, w = 0, w2 = 0, wdb = 0;
    void add(double x, double xdb) { n++; w += x; w2 += x * x; wdb += xdb; }
};

const double kD2R = 1.74532925199432955e-02;
inline TLorentzVector P4(double p, double th, double ph, double m) {
    const double t = th * kD2R, f = ph * kD2R;
    TLorentzVector v;
    v.SetVectM(TVector3(p * std::sin(t) * std::cos(f), p * std::sin(t) * std::sin(f), p * std::cos(t)), m);
    return v;
}

// ---- parse a partition table into key -> fields -----------------------------
typedef std::map<std::string, std::vector<std::string> > Table;
inline std::string Trim(const std::string& s) {
    size_t a = s.find_first_not_of(' '), b = s.find_last_not_of(' ');
    return a == std::string::npos ? "" : s.substr(a, b - a + 1);
}
inline Table Parse(const std::string& text, std::string& groupLine) {
    Table t;
    std::istringstream in(text);
    std::string line, sec;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        if (line.compare(0, 3, "## ") == 0) { sec = line.substr(3); continue; }
        if (line.compare(0, 11, "owned rows:") == 0) { groupLine = line; continue; }
        const bool diag = sec == "diagnostics";
        const size_t lab = diag ? 80 : 70;
        if (line.size() < lab) continue;
        std::string key = sec + " :: " + Trim(line.substr(0, lab));
        std::istringstream f(line.substr(lab));
        std::vector<std::string> v;
        std::string x;
        while (f >> x) v.push_back(x);
        t[key] = v;
    }
    return t;
}

struct CmpRes { int nRef = 0, nMine = 0, missing = 0, extra = 0, exact = 0, rounding = 0, bad = 0; double maxRel = 0; std::string firstBad; };
inline CmpRes Compare(const Table& ref, const Table& mine) {
    CmpRes r;
    r.nRef = (int)ref.size(); r.nMine = (int)mine.size();
    for (Table::const_iterator it = ref.begin(); it != ref.end(); ++it) {
        Table::const_iterator jt = mine.find(it->first);
        if (jt == mine.end()) { r.missing++; if (r.firstBad.empty()) r.firstBad = "missing: " + it->first; continue; }
        if (it->second == jt->second) { r.exact++; continue; }
        bool ok = it->second.size() == jt->second.size();
        for (size_t k = 0; ok && k < it->second.size(); ++k) {
            const double a = atof(it->second[k].c_str()), b = atof(jt->second[k].c_str());
            if (k == 0 && a != b) { ok = false; break; }                      // row counts must be exact
            const double rel = std::fabs(a - b) / std::max(std::fabs(a), 1e-30);
            // 5 significant digits printed (%.5g), N_eff with one decimal
            if (rel > 1.01e-4 && std::fabs(a - b) > 0.1001) { ok = false; break; }
            if (rel > r.maxRel && std::fabs(a - b) > 0) r.maxRel = rel;
        }
        if (ok) r.rounding++;
        else { r.bad++; if (r.firstBad.empty()) r.firstBad = "differs: " + it->first; }
    }
    for (Table::const_iterator jt = mine.begin(); jt != mine.end(); ++jt)
        if (!ref.count(jt->first)) { r.extra++; if (r.firstBad.empty()) r.firstBad = "extra: " + jt->first; }
    return r;
}

inline std::string ReadFile(const std::string& fn) {
    std::ifstream f(fn.c_str());
    std::stringstream s; s << f.rdbuf();
    return s.str();
}

} // namespace sval

void sim_validate(const char* which = "all") {
    using namespace sval;
    TStopwatch swAll;
    const std::string want = which;
    FILE* sum = fopen((std::string(kTabDir) + "/engine_validate_summary" + (want == "all" ? "" : "_" + want) + ".txt").c_str(), "w");
    if (!sum) { printf("cannot write summary\n"); return; }
    auto P = [&](const char* fmt, auto... a) { printf(fmt, a...); fprintf(sum, fmt, a...); };
    P("sim_validate: production %s\n", kOutDir);

    const char* LPn[12] = {"L0 misid", "L1 same pi0 Dalitz", "L2 same eta Dalitz", "L3 same conversion",
                           "L4 rho0/omega/phi->ee", "L5 omega->pi0ee", "L6 2nd-decay pi0 Dalitz",
                           "L7 two Dalitz (DD)", "L8 Dalitz+conversion (DC)", "L9 other different",
                           "L10 eta->ppee", "L11 X17"};
    const char* C4n[9] = {"fake misid", "fake pi->mu", "bump pi0 Dalitz", "bump conversion", "peak conversion",
                          "peak eta->ppee", "X17", "multi-pion", "CB-like"};
    const char* slot[4] = {"ep", "em", "pip", "pim"};
    const char* simc[18] = {"sim_id", "sim_parentid", "sim_grandparentid", "sim_processid", "sim_geninfo",
                            "sim_geninfo1", "sim_geninfo2", "sim_genweight", "sim_vertexx", "sim_vertexy",
                            "sim_vertexz", "sim_px", "sim_py", "sim_pz", "sim_isghost", "sim_mediumid",
                            "sim_p", "sim_iscommon"};
    std::set<std::string> refCols;
    bool allOk = true;

    for (int k = 0; k < 6; ++k) {
        const Def& D = kDefs[k];
        if (want != "all" && want != D.tag) continue;
        TStopwatch sw;
        const std::string fn = std::string(kOutDir) + "/" + D.file;
        TFile* f = TFile::Open(fn.c_str(), "READ");
        if (!f || f->IsZombie()) { P("[%s] cannot open %s\n", D.tag, fn.c_str()); allOk = false; continue; }
        TTree* t = (TTree*)f->Get("scan_nt");
        if (!t) { P("[%s] scan_nt missing\n", D.tag); allOk = false; continue; }
        const Long64_t n = t->GetEntries();

        // ---- 7. columns ----------------------------------------------------
        std::set<std::string> cols;
        Long64_t shortBranches = 0;
        TObjArray* bl = t->GetListOfBranches();
        for (int i = 0; i < bl->GetEntries(); ++i) {
            TBranch* b = (TBranch*)bl->At(i);
            cols.insert(b->GetName());
            if (b->GetEntries() != n) { shortBranches++; P("[%s] branch %s has %lld entries, tree %lld\n", D.tag, b->GetName(), b->GetEntries(), n); }
        }
        int missingSim = 0;
        for (int s = 0; s < 4; ++s)
            for (int c = 0; c < 18; ++c)
                if (!cols.count(Form("%s_%s", slot[s], simc[c]))) missingSim++;
        for (const char* c : {"ep_rich_sim_geanttrack", "em_rich_sim_geanttrack", "sim_set", "sim_cls", "sim_lcls",
                              "sim_pcls", "sim_tmpl", "sim_w", "sim_w_db", "sim_own_A", "sim_own_B", "sim_flags",
                              "sim_rules_ver", "m4_true", "mee_true", "mpipi_true"})
            if (!cols.count(c)) missingSim++;
        int colDiff = 0;
        if (refCols.empty()) refCols = cols;
        else { for (const auto& c : cols) if (!refCols.count(c)) colDiff++; for (const auto& c : refCols) if (!cols.count(c)) colDiff++; }
        // other trees in the file
        std::string others;
        TIter nk(f->GetListOfKeys());
        while (TKey* key = (TKey*)nk()) {
            if (strcmp(key->GetClassName(), "TTree") && strcmp(key->GetClassName(), "TNtuple")) continue;
            TTree* o = (TTree*)f->Get(key->GetName());
            if (!o) continue;
            others += Form(" %s;%d %lld x %d", key->GetName(), key->GetCycle(), o->GetEntries(), o->GetListOfBranches()->GetEntries());
        }

        // ---- branches --------------------------------------------------------
        float ev, best, trig, vz, corok, setv, cls, lcls, pcls, tmpl, w, wdb, ownA, ownB, flags, rver, m4t, meet, mppt;
        float m4c, meec, mppc, oac, oalc, oaec;
        float pc[4], th[4], ph[4];
        float v[4][18];
        t->SetBranchStatus("*", 0);
        auto bind = [&](const char* b, float* a) { t->SetBranchStatus(b, 1); t->SetBranchAddress(b, a); };
        bind("event", &ev); bind("isBest", &best); bind("trigbit", &trig); bind("eVertReco_z", &vz);
        bind("cor_ok", &corok); bind("sim_set", &setv); bind("sim_cls", &cls); bind("sim_lcls", &lcls);
        bind("sim_pcls", &pcls); bind("sim_tmpl", &tmpl); bind("sim_w", &w); bind("sim_w_db", &wdb);
        bind("sim_own_A", &ownA); bind("sim_own_B", &ownB); bind("sim_flags", &flags); bind("sim_rules_ver", &rver);
        bind("m4_true", &m4t); bind("mee_true", &meet); bind("mpipi_true", &mppt);
        bind("m_pippimepem_cor", &m4c); bind("m_ee_cor", &meec); bind("m_pippim_cor", &mppc); bind("oa", &oac);
        bind("oa_pippim_epem_lab_cor", &oalc); bind("oa_pippim_epem_eta_rest_cor", &oaec);
        for (int s = 0; s < 4; ++s) {
            bind(Form("%s_p_cor", slot[s]), &pc[s]); bind(Form("%s_theta", slot[s]), &th[s]); bind(Form("%s_phi", slot[s]), &ph[s]);
            for (int c = 0; c < 15; ++c) bind(Form("%s_%s", slot[s], simc[c]), &v[s][c]);
        }
        t->SetCacheSize(300 * 1024 * 1024);

        // ---- accumulators ----------------------------------------------------
        std::map<std::string, Acc> mOwn, mOwnLP, mOwnC4, mKey;
        std::map<std::string, long> diag;
        Long64_t nVz = 0, nTrigFail = 0, nBestNeg = 0, nVz500 = 0, nCorBad = 0;
        Long64_t nanW = 0, nanWdb = 0, nanM4 = 0, nanMee = 0, nanMpp = 0, nanOa = 0, nanOal = 0, nanOae = 0, negW = 0;
        Long64_t mis[16] = {0};
        const char* misName[16] = {"cls", "lcls", "pcls", "tmpl", "w", "w_db", "own_A", "own_B", "flags", "m4_true",
                                   "mee_true", "mpipi_true", "sim_set", "rules_ver", "m4_cor recomputed", ""};
        double maxDm4 = 0;
        Long64_t groups = 0, groupsVz = 0, badGroups = 0, badGroupsVz = 0, best1 = 0, best1Vz = 0;
        int nb = 0; bool grpVz = false;
        float pe = -1, pz = -1e9;
        double gsum = 0, gV = 0;
        Opt opt;
        const Set S = (Set)D.set;
        auto closeGroup = [&]() {
            if (groups == 0) return;
            if (nb != 1) { badGroups++; if (grpVz) badGroupsVz++; }
            gV += gsum * gsum; gsum = 0; nb = 0;
        };

        for (Long64_t r = 0; r < n; ++r) {
            t->GetEntry(r);
            if (r == 0 || ev != pe || vz != pz) {
                closeGroup();
                groups++;
                grpVz = vz > -160.f;
                if (grpVz) groupsVz++;
            }
            pe = ev; pz = vz;
            if (best == 1) { nb++; best1++; if (vz > -160.f) best1Vz++; }
            if (!((((int)trig) & 8192) != 0)) nTrigFail++;
            if (best < 0) nBestNeg++;
            if (!(vz > -500.f)) nVz500++;
            if (corok != 1.f) nCorBad++;
            if (w != w) nanW++;
            if (wdb != wdb) nanWdb++;
            if (w < 0) negW++;
            if (m4c != m4c) nanM4++;
            if (meec != meec) nanMee++;
            if (mppc != mppc) nanMpp++;
            if (oac != oac) nanOa++;
            if (oalc != oalc) nanOal++;
            if (oaec != oaec) nanOae++;

            // ---- 5. re-evaluation with sim_truth.h ----
            Trk T[4];
            for (int s = 0; s < 4; ++s) {
                const float* a = v[s];
                T[s] = Trk{(int)a[0], (int)a[1], (int)a[2], (int)a[3], a[4], a[5], a[6], a[7], a[8], a[9], a[10],
                           a[11], a[12], a[13], (int)a[14]};
            }
            const RowTruth R = EvaluateRow(S, T, D.q1, D.q2, opt);
            if (R.cls != (int)cls) mis[0]++;
            if (R.lcls != (int)lcls) mis[1]++;
            if (R.pcls != (int)pcls) mis[2]++;
            if (R.tmpl != (int)tmpl) mis[3]++;
            if ((float)R.w != w) mis[4]++;
            if ((float)R.wdb != wdb) mis[5]++;
            if ((int)R.ownA != (int)ownA) mis[6]++;
            if ((int)R.ownB != (int)ownB) mis[7]++;
            if (R.flags != (int)flags) mis[8]++;
            auto near = [](double a, float b) { return std::fabs(a - b) <= 2e-6 * std::max(1.0, std::fabs(a)); };
            if (!near(R.m4, m4t)) mis[9]++;
            if (!near(R.mee, meet)) mis[10]++;
            if (!near(R.mpipi, mppt)) mis[11]++;
            if ((int)setv != D.set) mis[12]++;
            if ((int)rver != kRulesVersion) mis[13]++;
            // ---- 6. M4 from the corrected momenta ----
            if (corok == 1.f) {
                const TLorentzVector a = P4(pc[0], th[0], ph[0], 0.51099895), b = P4(pc[1], th[1], ph[1], 0.51099895);
                const TLorentzVector c = P4(pc[2], th[2], ph[2], 139.56995), d = P4(pc[3], th[3], ph[3], 139.56995);
                const double m4 = (a + b + c + d).M() / 1000.0;
                const double dm = std::fabs(m4 - m4c);
                if (dm > maxDm4) maxDm4 = dm;
                if (dm > 1e-5) mis[14]++;
            }

            if (!(vz > -160.f)) continue;
            ++nVz;
            // ---- 2. tables from the stored columns (as partition.C) ----
            const bool own = ownA == 1.f;
            const double W = w, Wdb = wdb;
            const std::string o = own ? "owned  " : "dropped";
            mOwn[o].add(W, Wdb);
            if (ownB == 1.f) mOwn["mode B owned"].add(W, 0);
            if (own) gsum += W;
            const int lp = (int)lcls, c4 = (int)cls, pp = (int)pcls;
            mOwnLP[o + " | " + LPn[lp]].add(W, Wdb);
            mOwnC4[o + " | " + C4n[c4]].add(W, Wdb);
            const bool prov = (S == kStd) ? !own : own;
            if (prov) {
                const std::string pk = pp <= 1 ? "fake pions" : pp == 7 ? "both primary" : pp == 6 ? "K0S" : "mixed/secondary";
                mKey[std::string(LPn[lp]) + " + " + pk].add(W, Wdb);
                mKey["ALL provided by lep"].add(W, Wdb);
            }
            const int fl = (int)flags;
            if (fl & 1) diag["flag: stored genweight differs from table"]++;
            if (fl & 2) diag["flag: unknown forced code"]++;
            if (fl & 4) diag["flag: 2nd-decay lepton with omega lepton in row"]++;
            for (int i = 0; i < 4; ++i) {
                const Trk& x = T[i];
                if (FD(S, x)) diag[Form("FD %s g0 %.0f g1 %.0f", i < 2 ? "lepton slot" : "pion slot", x.g0, x.g1)]++;
                if (S != kStd && FC(S, x)) diag[Form("FC %s g1 %.0f proc %d", i < 2 ? "lepton slot" : "pion slot", x.g1, x.proc)]++;
                if (S == kStd && isE(x) && x.proc == 5 && x.par == 7 && x.gp == -1 && !FD(S, x)) diag[Form("std pi0-Dalitz e (gp -1) not FD, g0 %.0f", x.g0)]++;
                if (S == kStd && isE(x) && x.proc == 5 && x.par == 17 && !FD(S, x)) diag[Form("std eta-Dalitz e not FD, g0 %.0f gp %d", x.g0, x.gp)]++;
                if (S == kStd && LFM(x)) diag[Form("LFM %s g0 %.0f id %d par %d", i < 2 ? "lepton slot" : "pion slot", x.g0, x.id, x.par)]++;
                if (S != kStd && x.w > 0 && !forcedCode(S, x.g0)) diag[Form("w>0 with non-forced code %.0f", x.g0)]++;
            }
            if (S == kLep) {
                float cw = -1, cs = -1;
                for (int i = 0; i < 4; ++i) {
                    if (T[i].g0 == 525107.f && T[i].g1 == 52051.f) cw = T[i].g2;
                    if (T[i].g0 == 525107.f && T[i].g1 == 7051.f) cs = T[i].g2;
                }
                if (cw > 0 && cs > 0) diag[Form("omega counter %.0f with 2nd-decay counter %.0f", cw, cs)]++;
            }
        }
        closeGroup();

        // ---- write my table in the format of partition.C ----
        std::string mine;
        {
            char buf[512];
            snprintf(buf, sizeof buf, "# %s scan_nt set %d, vz > -160 of the FAT output: %lld of %lld rows\n", fn.c_str(), D.set, nVz, n);
            mine += buf;
            mine += "# columns: N, sum W (PDG), sum W^2, N_eff=(sum W)^2/sum W^2, sum W (database / note weights)\n";
            auto pr = [&](const char* title, std::map<std::string, Acc>& m) {
                snprintf(buf, sizeof buf, "\n## %s\n", title); mine += buf;
                for (auto& q : m) {
                    snprintf(buf, sizeof buf, "%-70s %10.0f %12.5g %12.5g %12.1f %12.5g\n", q.first.c_str(), q.second.n,
                             q.second.w, q.second.w2, q.second.w2 > 0 ? q.second.w * q.second.w / q.second.w2 : 0., q.second.wdb);
                    mine += buf;
                }
            };
            pr("ownership", mOwn);
            snprintf(buf, sizeof buf, "owned rows: event-group variance sum_g (sum W)^2 = %.6g, N_eff(group) = %.1f\n", gV,
                     gV > 0 ? mOwn["owned  "].w * mOwn["owned  "].w / gV : 0.);
            mine += buf;
            pr("ownership x lepton-pair class", mOwnLP);
            pr("ownership x four-particle class", mOwnC4);
            pr("rows the lepton set provides (std: dropped, lep: owned), by lepton pair + pions", mKey);
            mine += "\n## diagnostics\n";
            for (auto& q : diag) { snprintf(buf, sizeof buf, "%-80s %ld\n", q.first.c_str(), q.second); mine += buf; }
        }
        {
            FILE* o = fopen((std::string(kTabDir) + "/engine_validate_part_" + D.tag + ".txt").c_str(), "w");
            if (o) { fputs(mine.c_str(), o); fclose(o); }
        }
        std::string gRef, gMine;
        const Table tr = Parse(ReadFile(std::string(kRefDir) + "/" + D.ref), gRef);
        const Table tm = Parse(mine, gMine);
        const CmpRes cr = Compare(tr, tm);

        // ---- report ----
        const bool rowsOk = nVz == D.expVz && n == D.expAll;
        int nmis = 0; for (int i = 0; i < 15; ++i) nmis += mis[i] > 0;
        const bool ok = rowsOk && cr.missing == 0 && cr.extra == 0 && cr.bad == 0 && badGroups == 0 && nanW == 0 &&
                        nanWdb == 0 && nanM4 == 0 && nmis == 0 && shortBranches == 0 && missingSim == 0 && colDiff == 0 &&
                        nTrigFail == 0 && nBestNeg == 0 && nVz500 == 0 && gRef == gMine;
        if (!ok) allOk = false;
        P("\n[%s] %s  (%.1f s)  %s\n", D.tag, D.file, sw.RealTime(), ok ? "OK" : "PROBLEM");
        P("  rows %lld (expected %lld)  after vz > -160: %lld (expected %lld)\n", n, D.expAll, nVz, D.expVz);
        P("  owned (mode A) after vz > -160: %.0f rows, sum sim_w %.6g ; dropped %.0f rows, %.6g ; mode B owned %.0f rows, %.6g\n",
          mOwn["owned  "].n, mOwn["owned  "].w, mOwn["dropped"].n, mOwn["dropped"].w, mOwn["mode B owned"].n, mOwn["mode B owned"].w);
        P("  %s\n  reference: %s\n", gMine.c_str(), gRef.c_str());
        P("  table vs %s: %d reference lines, %d mine; identical %d, equal within print rounding %d (max rel %.2g), differ %d, missing %d, extra %d %s\n",
          D.ref, cr.nRef, cr.nMine, cr.exact, cr.rounding, cr.maxRel, cr.bad, cr.missing, cr.extra, cr.firstBad.c_str());
        P("  event groups %lld (vz > -160: %lld); groups without exactly one isBest==1 row: %lld (vz > -160: %lld); isBest==1 rows %lld (%lld)\n",
          groups, groupsVz, badGroups, badGroupsVz, best1, best1Vz);
        P("  trigger bit 8192 missing %lld ; isBest < 0 %lld ; vz <= -500 %lld ; cor_ok != 1 %lld\n", nTrigFail, nBestNeg, nVz500, nCorBad);
        P("  NaN: sim_w %lld, sim_w_db %lld, m4_cor %lld, mee_cor %lld, mpipi_cor %lld, oa %lld, oa_lab %lld, oa_eta %lld ; sim_w < 0 %lld\n",
          nanW, nanWdb, nanM4, nanMee, nanMpp, nanOa, nanOal, nanOae, negW);
        P("  re-evaluation with sim_truth.h, mismatching rows:");
        for (int i = 0; i < 15; ++i) P(" %s %lld,", misName[i], mis[i]);
        P("\n  max |M4 recomputed - m_pippimepem_cor| = %.3g GeV/c2\n", maxDm4);
        P("  columns %zu (sim columns missing %d, differing from the first output %d, branches with fewer entries %lld)\n",
          cols.size(), missingSim, colDiff, shortBranches);
        P("  trees:%s\n", others.c_str());
        f->Close();
        delete f;
    }
    P("\nsim_validate: %s, total %.1f s\n", allOk ? "ALL CHECKS PASSED" : "SOME CHECK FAILED (see above)", swAll.RealTime());
    fclose(sum);
}
