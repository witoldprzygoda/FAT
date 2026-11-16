// ============================================================================
// eestep_simdump.C  --  package eestep: simulated rows with lepton truth and the links between the e+e-, e+e+ and e-e-
//                       rows of one simulated event (charge-flip / wrong-sign-clone rate of the simulation)
// ============================================================================
// Written 2026-11-15 (prefix eestep). Read-only on the simulation compact files (SimEngine.h). The simulation has no
// combinatorial background (rule of the simulation study): these rows are used ONLY to measure how often a true lepton
// appears with the wrong charge, never for k or CB numbers.
// Rows: owned (mode A or B), S0 or pre (PassSel with cor_ok == 1, dz = 0, as edgeu_dump.C), 0.26 <= 0.9952 M(pipi) < 0.46,
// 0.29 <= M4 < 0.86, all three sets (std, lep, eta) and all nine trees.
// Event identity within a set: (event, vx, vy, vz) (generator event number + reconstructed vertex floats); every key must
// be one block of adjacent rows (checked). Track identity across trees: equal (theta, phi, p_rec).
// Tables (flat float32 column files <out>/<table>_<col>.f32; list in <out>/eestep_simlayout.txt):
//   sus  e+e- rows: group g (edgeu_dump.C GroupOf: 0 C2 .. 8 GCnom, -1 other), truth, and the like-sign twins
//        (twmm: an e-e- row with the same pi+, pi-, the same e- and another negative candidate; twpp: e+e+ row with the
//        same pi+, pi-, e+ and another positive candidate); for the twin: the true GEANT id of the other candidate
//        (odd), whether it is the same true particle as the e+e- row's opposite lepton (equal true momentum), angle, p ratio
//   sls  like-sign rows (e+e+ k = 1, e-e- k = 2) with the truth of both slots and the e+e- twin of the same event
// Run (from /home/przygoda/HADES/FAT/FAT; about 10 GB, a few minutes):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_eestep", kTRUE)' \
//        'plots/all_report/eestep_simdump.C+O'
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

#include "TStopwatch.h"
#include "TSystem.h"

#include "sim_eta_physics.h"
#include "bedge_common.h"

namespace eestepsim {
using namespace bedge;
const char* const kOutE = "/hdd3/przygoda/HADES/FAT_out/all_report/eestep";
FILE* gLog = nullptr;
#define PL(...) do { printf(__VA_ARGS__); if (gLog) fprintf(gLog, __VA_ARGS__); fflush(stdout); } while (0)

struct Tab {
    std::string name;
    std::vector<std::string> cols;
    std::map<std::string, int> idx;
    std::vector<std::vector<float>> v;
    std::vector<float> row;
    size_t nrows = 0;
    int C(const std::string& c) {
        auto it = idx.find(c);
        if (it != idx.end()) return it->second;
        const int i = (int)cols.size(); cols.push_back(c); idx[c] = i; v.emplace_back(nrows, NAN); row.push_back(NAN); return i;
    }
    void Set(const std::string& c, double x) { row[C(c)] = (float)x; }
    void Push() {
        for (size_t i = 0; i < cols.size(); ++i) v[i].push_back(row[i]);
        std::fill(row.begin(), row.end(), NAN);
        ++nrows;
    }
    void Write(FILE* lay) {
        for (size_t i = 0; i < cols.size(); ++i) {
            FILE* fp = fopen(Form("%s/%s_%s.f32", kOutE, name.c_str(), cols[i].c_str()), "wb");
            if (!fp) continue;
            fwrite(v[i].data(), 4, v[i].size(), fp); fclose(fp);
        }
        fprintf(lay, "%s %zu %zu\n", name.c_str(), nrows, cols.size());
        for (auto& c : cols) fprintf(lay, "  %s\n", c.c_str());
    }
};
inline int GroupOf(const int* cn, int ncn, bool& isT12) {
    isT12 = false; int g = -1;
    for (int q = 0; q < ncn; ++q) {
        const int c = cn[q];
        if (c == kT12) isT12 = true;
        if (g >= 0) continue;
        switch (c) {
            case kC2: g = gC2; break;
            case kC3: g = gC3; break;
            case kC4: g = gC4; break;
            case kC5: g = gC5; break;
            case kGeta: g = gGeta; break;
            case kGpi0Dom: case kGvecom: case kGothom: g = gGDom; break;
            case kGpi0Dnom: case kGvecnom: case kGothnom: g = gGDnom; break;
            case kGconvom: g = gGCom; break;
            case kGconvnom: g = gGCnom; break;
        }
    }
    return g;
}
std::string KeyOf(const simeng::Row& r) {
    uint32_t a, b, c;
    std::memcpy(&a, &r.vx, 4); std::memcpy(&b, &r.vy, 4); std::memcpy(&c, &r.vz, 4);
    char s[64];
    snprintf(s, sizeof(s), "%d_%d_%08x_%08x_%08x", r.set, r.event, a, b, c);
    return s;
}
struct T3 { float th, ph, p; };
inline T3 Lep(const simeng::Row& r, int s) { return s == 0 ? T3{r.ep_theta, r.ep_phi, r.ep_p_rec} : T3{r.em_theta, r.em_phi, r.em_p_rec}; }
inline bool Same(const T3& a, const T3& b) { return a.th == b.th && a.ph == b.ph && a.p == b.p; }
inline bool SamePions(const simeng::Row& a, const simeng::Row& b) {
    return a.pip_theta == b.pip_theta && a.pip_phi == b.pip_phi && a.pip_p_rec == b.pip_p_rec &&
           a.pim_theta == b.pim_theta && a.pim_phi == b.pim_phi && a.pim_p_rec == b.pim_p_rec;
}
inline double Ang(const T3& a, const T3& b) {
    const double d = M_PI / 180.;
    const double c = std::sin(a.th * d) * std::sin(b.th * d) * std::cos((a.ph - b.ph) * d) + std::cos(a.th * d) * std::cos(b.th * d);
    return std::acos(std::max(-1., std::min(1., c))) / d;
}
inline int SimId(const simeng::Row& r, int s) { return s == 0 ? r.ep_sim_id : r.em_sim_id; }
inline float SimP(const simeng::Row& r, int s) { return s == 0 ? r.ep_sim_p : r.em_sim_p; }
}  // namespace eestepsim

void eestep_simdump() {
    using namespace eestepsim;
    TStopwatch total;
    gSystem->mkdir(kOutE, kTRUE);
    gLog = fopen(Form("%s/eestep_simdump.log", kOutE), "w");
    PL("# eestep_simdump.C, 2026-11-15: simulated rows with lepton truth and like-sign twins\n");
    simphys::Init();
    simphys::Install();
    simphys::LoadAngles();
    simeng::SimData* D = new simeng::SimData();
    if (!simeng::Load(*D, simeng::kCompact, simeng::kMaskAll | simeng::kMaskEta, -1e9, 1e9, 1)) return;
    Tab TU, TL;
    TU.name = "sus"; TL.name = "sls";
    const int setTrees[3][3] = {{simeng::kStdUS, simeng::kStdPP, simeng::kStdMM}, {simeng::kLepUS, simeng::kLepPP, simeng::kLepMM},
                                {simeng::kEtaUS, simeng::kEtaPP, simeng::kEtaMM}};
    for (int st = 0; st < 3; ++st) {
        // selection and key maps of the three trees of this set
        std::vector<char> ok[3];
        std::unordered_map<std::string, std::vector<int>> M[3];
        long nsplit[3] = {0, 0, 0};
        for (int k = 0; k < 3; ++k) {
            const auto& v = D->t[setTrees[st][k]];
            ok[k].assign(v.size(), 0);
            std::string last;
            M[k].reserve(v.size());
            for (size_t i = 0; i < v.size(); ++i) {
                const auto& r = v[i];
                const std::string key = KeyOf(r);
                auto it = M[k].find(key);
                if (it != M[k].end() && key != last) ++nsplit[k];
                M[k][key].push_back((int)i);
                last = key;
                if (r.ownA != 1 && r.ownB != 1) continue;
                const bool s0 = PassSel(kSelS0, r, r.cor_ok == 1, 0.f, r.mm4);
                const bool pre = PassSel(kSelPre, r, r.cor_ok == 1, 0.f, r.mm4);
                if (!s0 && !pre) continue;
                const double mpp = kSPiPi * r.mpipi;
                if (!(mpp >= 0.26 && mpp < 0.46)) continue;
                if (!(r.m4 >= 0.29f && r.m4 < 0.86f)) continue;
                ok[k][i] = 1 | (s0 ? 2 : 0) | (pre ? 4 : 0);
            }
        }
        PL("set %d: rows us %zu pp %zu mm %zu; keys re-appearing after another key: %ld %ld %ld\n", st, D->t[setTrees[st][0]].size(),
           D->t[setTrees[st][1]].size(), D->t[setTrees[st][2]].size(), nsplit[0], nsplit[1], nsplit[2]);
        // ---- e+e- rows
        const auto& U = D->t[setTrees[st][0]];
        long nU = 0, nTw[3] = {0, 0, 0};
        for (size_t i = 0; i < U.size(); ++i) {
            if (!ok[0][i]) continue;
            const auto& r = U[i];
            int cn[6]; bool isT12 = false;
            const int ncn = CompsOf(r, false, cn);
            const int g = ncn > 0 ? GroupOf(cn, ncn, isT12) : -1;
            if (st > 0 && g != gC2 && g != gC3 && g != gC4 && g != gC5) {
                // lep / eta sets: keep only a 1/10 subsample of the other rows (rate denominators), weight x 10 recorded
                if ((r.grp % 10) != 0) continue;
            }
            const double sub = (st > 0 && g != gC2 && g != gC3 && g != gC4 && g != gC5) ? 10. : 1.;
            TU.Set("set", st); TU.Set("tree", setTrees[st][0]); TU.Set("grp", r.grp); TU.Set("g", g); TU.Set("cls", r.cls);
            TU.Set("lcls", r.lcls); TU.Set("tmpl", r.tmpl); TU.Set("flags", r.flags); TU.Set("ghost", r.ghost);
            TU.Set("bits", ok[0][i] | (r.isBest == 1 ? 8 : 0) | (r.ownA == 1 ? 16 : 0) | (r.ownB == 1 ? 32 : 0));
            TU.Set("m4", r.m4); TU.Set("m4t", r.m4_true); TU.Set("mpipi", r.mpipi); TU.Set("mpipit", r.mpipi_true);
            TU.Set("mee", r.mee); TU.Set("meet", r.mee_true); TU.Set("oa", r.oa);
            TU.Set("w", r.w * kSimToRaw * sub); TU.Set("f108", simphys::RowFactor(r, kPhysNom));
            for (int s = 0; s < 2; ++s) {
                const char* L = s == 0 ? "l1" : "l2";
                TU.Set(Form("%s_id", L), SimId(r, s)); TU.Set(Form("%s_par", L), s == 0 ? r.ep_sim_parentid : r.em_sim_parentid);
                TU.Set(Form("%s_proc", L), s == 0 ? r.ep_sim_processid : r.em_sim_processid);
                TU.Set(Form("%s_p", L), s == 0 ? r.ep_p : r.em_p); TU.Set(Form("%s_simp", L), SimP(r, s));
                TU.Set(Form("%s_theta", L), s == 0 ? r.ep_theta : r.em_theta); TU.Set(Form("%s_phi", L), s == 0 ? r.ep_phi : r.em_phi);
                TU.Set(Form("%s_simtheta", L), s == 0 ? r.ep_sim_theta : r.em_sim_theta);
                TU.Set(Form("%s_system", L), s == 0 ? r.ep_system : r.em_system);
            }
            // twins: q = 2 (mm: shares the e-, odd = negative candidate vs our e+), q = 1 (pp: shares the e+, odd vs our e-)
            for (int q = 1; q <= 2; ++q) {
                const char* tg = q == 2 ? "twmm" : "twpp";
                const int shared = q == 2 ? 1 : 0, oddOwn = 1 - shared;
                const auto& V = D->t[setTrees[st][q]];
                auto it = M[q].find(KeyOf(r));
                int type = it == M[q].end() ? 0 : 1, bj = -1, bodd = -1; double ba = 1e9;
                if (it != M[q].end())
                    for (int j : it->second) {
                        const auto& o = V[j];
                        if (!SamePions(r, o)) continue;
                        for (int s = 0; s < 2; ++s) {
                            if (!Same(Lep(r, shared), Lep(o, s))) continue;
                            const double a = Ang(Lep(r, oddOwn), Lep(o, 1 - s));
                            if (a < ba) { ba = a; bj = j; bodd = 1 - s; }
                        }
                    }
                if (bj >= 0) { type = 2; ++nTw[q]; }
                TU.Set(Form("%s_type", tg), type);
                if (bj >= 0) {
                    const auto& o = V[bj];
                    TU.Set(Form("%s_ang", tg), ba);
                    TU.Set(Form("%s_prat", tg), Lep(o, bodd).p / Lep(r, oddOwn).p);
                    TU.Set(Form("%s_oddid", tg), SimId(o, bodd));
                    TU.Set(Form("%s_samepart", tg), (SimP(o, bodd) == SimP(r, oddOwn) && SimId(o, bodd) == SimId(r, oddOwn)) ? 1 : 0);
                    TU.Set(Form("%s_sel", tg), ok[q][bj]);
                    TU.Set(Form("%s_m4", tg), o.m4); TU.Set(Form("%s_mee", tg), o.mee); TU.Set(Form("%s_mpipi", tg), o.mpipi);
                    TU.Set(Form("%s_oddp", tg), Lep(o, bodd).p); TU.Set(Form("%s_ghost", tg), o.ghost);
                } else {
                    for (const char* c : {"ang", "prat", "oddid", "samepart", "sel", "m4", "mee", "mpipi", "oddp", "ghost"}) TU.C(Form("%s_%s", tg, c));
                }
            }
            TU.Push();
            ++nU;
        }
        // ---- like-sign rows
        for (int k = 1; k <= 2; ++k) {
            const auto& V = D->t[setTrees[st][k]];
            for (size_t i = 0; i < V.size(); ++i) {
                if (!ok[k][i]) continue;
                const auto& r = V[i];
                TL.Set("set", st); TL.Set("k", k); TL.Set("tree", setTrees[st][k]); TL.Set("grp", r.grp); TL.Set("cls", r.cls);
                TL.Set("lcls", r.lcls); TL.Set("tmpl", r.tmpl); TL.Set("flags", r.flags); TL.Set("ghost", r.ghost);
                TL.Set("bits", ok[k][i] | (r.isBest == 1 ? 8 : 0) | (r.ownA == 1 ? 16 : 0) | (r.ownB == 1 ? 32 : 0));
                TL.Set("m4", r.m4); TL.Set("m4t", r.m4_true); TL.Set("mpipi", r.mpipi); TL.Set("mpipit", r.mpipi_true);
                TL.Set("mee", r.mee); TL.Set("meet", r.mee_true); TL.Set("oa", r.oa);
                TL.Set("w", r.w * kSimToRaw); TL.Set("f108", simphys::RowFactor(r, kPhysNom));
                for (int s = 0; s < 2; ++s) {
                    const char* L = s == 0 ? "l1" : "l2";
                    TL.Set(Form("%s_id", L), SimId(r, s)); TL.Set(Form("%s_par", L), s == 0 ? r.ep_sim_parentid : r.em_sim_parentid);
                    TL.Set(Form("%s_gpar", L), s == 0 ? r.ep_sim_grandparentid : r.em_sim_grandparentid);
                    TL.Set(Form("%s_proc", L), s == 0 ? r.ep_sim_processid : r.em_sim_processid);
                    TL.Set(Form("%s_p", L), s == 0 ? r.ep_p : r.em_p); TL.Set(Form("%s_simp", L), SimP(r, s));
                    TL.Set(Form("%s_theta", L), s == 0 ? r.ep_theta : r.em_theta); TL.Set(Form("%s_phi", L), s == 0 ? r.ep_phi : r.em_phi);
                    TL.Set(Form("%s_simtheta", L), s == 0 ? r.ep_sim_theta : r.em_sim_theta);
                    TL.Set(Form("%s_system", L), s == 0 ? r.ep_system : r.em_system);
                    TL.Set(Form("%s_rkchi2", L), s == 0 ? r.ep_rkchi2 : r.em_rkchi2);
                }
                // e+e- twin: same pions and one shared lepton; odd = our other slot vs the twin's opposite-charge lepton
                const auto& Uv = D->t[setTrees[st][0]];
                auto it = M[0].find(KeyOf(r));
                const int opp = (k == 2) ? 0 : 1, same = 1 - opp;
                int type = it == M[0].end() ? 0 : 1, bj = -1, bodd = -1; double ba = 1e9;
                if (it != M[0].end())
                    for (int j : it->second) {
                        const auto& o = Uv[j];
                        if (!SamePions(r, o)) continue;
                        for (int s = 0; s < 2; ++s) {
                            if (!Same(Lep(r, s), Lep(o, same))) continue;
                            const double a = Ang(Lep(r, 1 - s), Lep(o, opp));
                            if (a < ba) { ba = a; bj = j; bodd = 1 - s; }
                        }
                    }
                if (bj >= 0) type = 2;
                TL.Set("lk_type", type);
                if (bj >= 0) {
                    const auto& o = Uv[bj];
                    int cn[6]; bool isT12 = false;
                    const int ncn = CompsOf(o, false, cn);
                    TL.Set("lk_odd", bodd + 1); TL.Set("lk_ang", ba); TL.Set("lk_prat", Lep(r, bodd).p / Lep(o, opp).p);
                    TL.Set("lk_samepart", (SimP(r, bodd) == SimP(o, opp) && SimId(r, bodd) == SimId(o, opp)) ? 1 : 0);
                    TL.Set("lk_g", ncn > 0 ? GroupOf(cn, ncn, isT12) : -1); TL.Set("lk_sel", ok[0][bj]);
                    TL.Set("lk_m4", o.m4); TL.Set("lk_mee", o.mee); TL.Set("lk_mpipi", o.mpipi);
                } else {
                    for (const char* c : {"lk_odd", "lk_ang", "lk_prat", "lk_samepart", "lk_g", "lk_sel", "lk_m4", "lk_mee", "lk_mpipi"}) TL.C(c);
                }
                TL.Push();
            }
        }
        PL("set %d: e+e- rows written %ld; with an e+e+ twin %ld, with an e-e- twin %ld; table sls now %zu rows\n", st, nU, nTw[1], nTw[2], TL.nrows);
    }
    delete D;
    FILE* lay = fopen(Form("%s/eestep_simlayout.txt", kOutE), "w");
    TU.Write(lay); TL.Write(lay);
    fclose(lay);
    PL("tables: sus %zu rows %zu cols; sls %zu rows %zu cols; done in %.0f s\n", TU.nrows, TU.cols.size(), TL.nrows, TL.cols.size(), total.RealTime());
    if (gLog) fclose(gLog);
}
