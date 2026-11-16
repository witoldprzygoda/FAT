// ============================================================================
// veestep_dump.C  --  INDEPENDENT verification of the eestep finding (prefix veestep): data row dump with event links
// ============================================================================
// Written 2026-11-15 by the independent verifier of eestep (the e-e- like-sign step at u = 0). No eestep_* code, output or
// number is used. Inputs: the data compact rows (scan_v2, ScanEngine.h) with the frozen selection of defs.h, and the PAT
// files /hdd3/przygoda/HADES/files/exp/gen4/lep_<day>.root (read only).
//
// Part A: rows of the edgeu range (cor_ok == 1, pre, 0.26 <= M(pipi) < 0.46, 0.30 <= M4 < 0.85), all three samples, with
//         the columns k, m4, mpipi, mee, oa, bits (1 S0, 2 pre, 4 isBest == 1), day, pat_entry, grp, nrow, n_lep, lepton and
//         pion momenta and angles, systems, shared_ring, vz.  (check: same row count and same m4 column as edgeu_dump.C)
// Part B: for every LIKE-SIGN row: its PAT entry (lpos_mult, lneg_mult, hpos_mult, hneg_mult, totalmult, isBest, both lepton
//         angles for a join check) and the e+e- content of the same event in PipPimEpEm_ID. The event is found by a key of
//         reconstructed event values (runnumber, event, start_time, eVertReco_z, totalmult); the event counter alone is not
//         unique (it restarts per sub-file). Every key of the e+e- tree is required to be ONE block of adjacent entries:
//         a key seen again after other keys is counted (nonadj) and its rows are not used for links.
//         Per like-sign row: n_us = e+e- entries of the event, n_twin = e+e- entries with the same pi+, pi- (theta, phi)
//         and one lepton equal (theta, phi) to one of the two like-sign leptons.
//         For the e+e- rows of the edge region (|u| < 0.1): n_ls = e-e- entries of the same event, n_flip = e-e- entries with
//         the same pi+, pi-, the same e- and a second e- within 3 deg of the e+ (a negative copy of the e+).
// Run (from /home/przygoda/HADES/FAT/FAT):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_veestep", kTRUE)' \
//        'plots/all_report/veestep_dump.C+O'
// ============================================================================
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#include "TFile.h"
#include "TNtuple.h"
#include "TStopwatch.h"
#include "TSystem.h"

#include "defs.h"

namespace veed {
const char* const kOut = "/hdd3/przygoda/HADES/FAT_out/all_report/veestep";
const char* const kPat = "/hdd3/przygoda/HADES/files/exp/gen4";
FILE* gLog = nullptr;
#define VL(...) do { printf(__VA_ARGS__); if (gLog) fprintf(gLog, __VA_ARGS__); fflush(stdout); } while (0)
template <class T> void W(const std::string& n, const std::vector<T>& v) {
    FILE* fp = fopen((std::string(kOut) + "/" + n).c_str(), "wb");
    if (!fp) { VL("cannot write %s\n", n.c_str()); return; }
    fwrite(v.data(), sizeof(T), v.size(), fp);
    fclose(fp);
}
struct Key {
    uint32_t a[5];
    bool operator==(const Key& o) const { return std::memcmp(a, o.a, sizeof(a)) == 0; }
};
struct KH {
    size_t operator()(const Key& k) const {
        uint64_t h = 1469598103934665603ULL;
        for (int i = 0; i < 5; ++i) { h ^= k.a[i]; h *= 1099511628211ULL; }
        return (size_t)h;
    }
};
inline uint32_t FB(float x) { uint32_t u; std::memcpy(&u, &x, 4); return u; }
inline Key MakeKey(float run, float ev, float st, float vz, float tm) { return Key{{FB(run), FB(ev), FB(st), FB(vz), FB(tm)}}; }
inline double Dang(double t1, double p1, double t2, double p2) {   // angle between two directions [deg]
    const double d = M_PI / 180.;
    const double c = std::sin(t1 * d) * std::sin(t2 * d) * std::cos((p1 - p2) * d) + std::cos(t1 * d) * std::cos(t2 * d);
    return std::acos(std::max(-1.0, std::min(1.0, c))) / d;
}
inline bool Same(float t1, float p1, float t2, float p2) { return t1 == t2 && p1 == p2; }
// compact (phi in -180..180) against PAT (phi in 0..360)
inline bool SameC(float tc, float pc, float tp, float pp) {
    return std::fabs(tc - tp) < 1e-3f && std::fabs(std::remainder((double)pc - (double)pp, 360.0)) < 1e-3;
}
}  // namespace veed

void veestep_dump() {
    using namespace veed;
    TStopwatch total;
    gSystem->mkdir(kOut, kTRUE);
    gLog = fopen(Form("%s/veestep_dump.log", kOut), "w");
    VL("# veestep_dump.C (independent verifier of eestep), 2026-11-15\n# data %s\n", scan::kCompact);
    scan::Data d;
    if (!scan::Load(d)) return;
    // ------------------------------------------------------------------ Part A
    std::vector<int8_t> ck; std::vector<uint8_t> cbits; std::vector<int16_t> cday;
    std::vector<int32_t> cent, cgrp, cnrow, cnlep, cs1, cs2, csr;
    std::vector<float> cm4, cpp, cee, coa, cp1, cp2, ct1, ct2, cf1, cf2, cpip, cpim, ctpip, ctpim, cfpip, cfpim, cvz;
    long nk[3] = {0, 0, 0};
    for (int k = 0; k < 3; ++k) {
        for (const auto& r : d.s[k]) {
            if (r.cor_ok != 1) continue;
            if (!alldefs::PassSel(alldefs::kSelPre, r, true, r.dz_pi, r.mm4)) continue;
            const bool s0 = alldefs::PassSel(alldefs::kSelS0, r, true, r.dz_pi, r.mm4);
            if (!(r.mpipi >= 0.26f && r.mpipi < 0.46f)) continue;
            if (!(r.m4 >= 0.30f && r.m4 < 0.85f)) continue;
            ck.push_back((int8_t)k); cbits.push_back((uint8_t)((s0 ? 1 : 0) | 2 | (r.isBest == 1 ? 4 : 0)));
            cday.push_back((int16_t)r.file_id); cent.push_back(r.pat_entry); cgrp.push_back(r.grp); cnrow.push_back(r.nrow);
            cnlep.push_back(r.n_lep); cs1.push_back(r.ep_system); cs2.push_back(r.em_system); csr.push_back(r.shared_ring);
            cm4.push_back(r.m4); cpp.push_back(r.mpipi); cee.push_back(r.mee); coa.push_back(r.oa);
            cp1.push_back(r.ep_p); cp2.push_back(r.em_p); ct1.push_back(r.ep_theta); ct2.push_back(r.em_theta);
            cf1.push_back(r.ep_phi); cf2.push_back(r.em_phi); cpip.push_back(r.pip_p); cpim.push_back(r.pim_p);
            ctpip.push_back(r.pip_theta); ctpim.push_back(r.pim_theta); cfpip.push_back(r.pip_phi); cfpim.push_back(r.pim_phi);
            cvz.push_back(r.vz);
            ++nk[k];
        }
    }
    const size_t n = ck.size();
    VL("A. rows kept: us %ld pp %ld mm %ld (total %zu)\n", nk[0], nk[1], nk[2], n);
    d.s[0].clear(); d.s[1].clear(); d.s[2].clear();
    d.s[0].shrink_to_fit(); d.s[1].shrink_to_fit(); d.s[2].shrink_to_fit();

    // ------------------------------------------------------------------ Part B
    std::vector<int16_t> lpos(n, -1), lneg(n, -1), hpos(n, -1), hneg(n, -1), tmult(n, -1), nus(n, -1), ntwin(n, -1), nls(n, -1), nflip(n, -1);
    std::vector<uint8_t> joinok(n, 0);
    std::map<int, std::vector<size_t>> byday;
    for (size_t i = 0; i < n; ++i) byday[cday[i]].push_back(i);
    long nJoinBad = 0, nJoin = 0, nNonAdj = 0, nNoKey = 0;
    for (auto& kv : byday) {
        TStopwatch sw;
        const int day = kv.first;
        TFile* f = TFile::Open(Form("%s/lep_%03d.root", kPat, day), "READ");
        if (!f || f->IsZombie()) { VL("cannot open day %d\n", day); continue; }
        // --- e+e- tree: keys and the track directions of every entry
        TNtuple* tu = (TNtuple*)f->Get("PipPimEpEm_ID");
        tu->SetBranchStatus("*", 0);
        const char* ub[13] = {"runnumber", "event", "start_time", "eVertReco_z", "totalmult", "ep_theta", "ep_phi", "em_theta",
                              "em_phi", "pip_theta", "pip_phi", "pim_theta", "pim_phi"};
        float uv[13];
        for (int j = 0; j < 13; ++j) { tu->SetBranchStatus(ub[j], 1); tu->SetBranchAddress(ub[j], &uv[j]); }
        const Long64_t nu = tu->GetEntries();
        std::vector<float> U(nu * 8);
        std::unordered_map<Key, std::pair<Long64_t, int>, KH> umap;
        std::unordered_map<Key, int, KH> ubad;
        Key prev{{0, 0, 0, 0, 0}}; bool havePrev = false;
        for (Long64_t e = 0; e < nu; ++e) {
            tu->GetEntry(e);
            for (int j = 0; j < 8; ++j) U[e * 8 + j] = uv[5 + j];
            const Key kk = MakeKey(uv[0], uv[1], uv[2], uv[3], uv[4]);
            if (havePrev && kk == prev) { umap[kk].second++; continue; }
            auto it = umap.find(kk);
            if (it != umap.end()) { ubad[kk]++; ++nNonAdj; }
            else umap[kk] = std::make_pair(e, 1);
            prev = kk; havePrev = true;
        }
        // --- e-e- and e+e+ trees: the like-sign rows and the event keys of the e-e- tree (for the flip copies)
        for (int k = 1; k <= 2; ++k) {
            TNtuple* t = (TNtuple*)f->Get(k == 1 ? "PipPimEpEp_ID" : "PipPimEmEm_ID");
            const char* s1 = k == 1 ? "ep1" : "em1";
            const char* s2 = k == 1 ? "ep2" : "em2";
            t->SetBranchStatus("*", 0);
            std::vector<std::string> bn = {"runnumber", "event", "start_time", "eVertReco_z", "totalmult", "lpos_mult", "lneg_mult",
                                           "hpos_mult", "hneg_mult", Form("%s_theta", s1), Form("%s_phi", s1), Form("%s_theta", s2),
                                           Form("%s_phi", s2), "pip_theta", "pip_phi", "pim_theta", "pim_phi"};
            float v[17];
            for (size_t j = 0; j < bn.size(); ++j) { t->SetBranchStatus(bn[j].c_str(), 1); t->SetBranchAddress(bn[j].c_str(), &v[j]); }
            for (size_t i : kv.second) {
                if (ck[i] != k) continue;
                t->GetEntry(cent[i]);
                ++nJoin;
                // join check: the compact lepton angles equal the PAT ones (slot 1 = ep, slot 2 = em of the compact row)
                const bool ok = SameC(ct1[i], cf1[i], v[9], v[10]) && SameC(ct2[i], cf2[i], v[11], v[12]) &&
                                SameC(ctpip[i], cfpip[i], v[13], v[14]) && SameC(ctpim[i], cfpim[i], v[15], v[16]);
                if (!ok) { ++nJoinBad; continue; }
                joinok[i] = 1;
                lpos[i] = (int16_t)v[5]; lneg[i] = (int16_t)v[6]; hpos[i] = (int16_t)v[7]; hneg[i] = (int16_t)v[8]; tmult[i] = (int16_t)v[4];
                const Key kk = MakeKey(v[0], v[1], v[2], v[3], v[4]);
                if (ubad.count(kk)) { nus[i] = -2; continue; }
                auto it = umap.find(kk);
                if (it == umap.end()) { nus[i] = 0; ntwin[i] = 0; ++nNoKey; continue; }
                int nt = 0;
                for (Long64_t e = it->second.first; e < it->second.first + it->second.second; ++e) {
                    const float* x = &U[e * 8];
                    if (!(Same(x[4], x[5], v[13], v[14]) && Same(x[6], x[7], v[15], v[16]))) continue;
                    // the lepton of the like-sign charge in the e+e- entry: e+ (x[0], x[1]) for e+e+, e- (x[2], x[3]) for e-e-
                    const float tl = k == 1 ? x[0] : x[2], pl = k == 1 ? x[1] : x[3];
                    if (Same(tl, pl, v[9], v[10]) || Same(tl, pl, v[11], v[12])) ++nt;
                }
                nus[i] = (int16_t)it->second.second; ntwin[i] = (int16_t)nt;
            }
        }
        // --- e+e- rows of the edge region: e-e- entries of the same event and negative copies of the e+
        {
            TNtuple* tm = (TNtuple*)f->Get("PipPimEmEm_ID");
            tm->SetBranchStatus("*", 0);
            const char* mb[13] = {"runnumber", "event", "start_time", "eVertReco_z", "totalmult", "em1_theta", "em1_phi", "em2_theta",
                                  "em2_phi", "pip_theta", "pip_phi", "pim_theta", "pim_phi"};
            float mv[13];
            for (int j = 0; j < 13; ++j) { tm->SetBranchStatus(mb[j], 1); tm->SetBranchAddress(mb[j], &mv[j]); }
            const Long64_t nm = tm->GetEntries();
            std::vector<float> M(nm * 8);
            std::unordered_map<Key, std::pair<Long64_t, int>, KH> mmap;
            std::unordered_map<Key, int, KH> mbad;
            Key pv{{0, 0, 0, 0, 0}}; bool hp = false;
            for (Long64_t e = 0; e < nm; ++e) {
                tm->GetEntry(e);
                for (int j = 0; j < 8; ++j) M[e * 8 + j] = mv[5 + j];
                const Key kk = MakeKey(mv[0], mv[1], mv[2], mv[3], mv[4]);
                if (hp && kk == pv) { mmap[kk].second++; continue; }
                if (mmap.count(kk)) { mbad[kk]++; ++nNonAdj; }
                else mmap[kk] = std::make_pair(e, 1);
                pv = kk; hp = true;
            }
            // e+e- PAT entries of the selected us rows of the edge region
            for (size_t i : kv.second) {
                if (ck[i] != 0) continue;
                const double u = cm4[i] - alldefs::Edge(cpp[i], cee[i]);
                if (!(cee[i] < scan::kMpi0) || alldefs::Edge(cpp[i], cee[i]) < 0 || std::fabs(u) >= 0.100) continue;
                const Long64_t e = cent[i];
                if (e < 0 || e >= nu) { ++nJoinBad; continue; }
                tu->GetEntry(e);
                ++nJoin;
                if (!(SameC(ct1[i], cf1[i], uv[5], uv[6]) && SameC(ct2[i], cf2[i], uv[7], uv[8]))) { ++nJoinBad; continue; }
                joinok[i] = 1;
                const Key kk = MakeKey(uv[0], uv[1], uv[2], uv[3], uv[4]);
                if (mbad.count(kk)) { nls[i] = -2; continue; }
                auto it = mmap.find(kk);
                if (it == mmap.end()) { nls[i] = 0; nflip[i] = 0; continue; }
                int nf = 0;
                for (Long64_t q = it->second.first; q < it->second.first + it->second.second; ++q) {
                    const float* x = &M[q * 8];
                    if (!(Same(x[4], x[5], uv[9], uv[10]) && Same(x[6], x[7], uv[11], uv[12]))) continue;
                    // same e- in one slot, the other slot within 3 deg of the e+
                    if (Same(x[0], x[1], uv[7], uv[8]) && Dang(x[2], x[3], uv[5], uv[6]) < 3.0) ++nf;
                    else if (Same(x[2], x[3], uv[7], uv[8]) && Dang(x[0], x[1], uv[5], uv[6]) < 3.0) ++nf;
                }
                nls[i] = (int16_t)it->second.second; nflip[i] = (int16_t)nf;
            }
        }
        VL("B. day %03d: us entries %lld, rows %zu, %.0f s\n", day, nu, kv.second.size(), sw.RealTime());
        f->Close();
        delete f;
    }
    VL("B. joins %ld, failed join checks %ld, keys seen in two non-adjacent blocks %ld, like-sign rows without e+e- entry %ld\n",
       nJoin, nJoinBad, nNonAdj, nNoKey);
    W("d_k.i8", ck); W("d_bits.u8", cbits); W("d_day.i16", cday); W("d_ent.i32", cent); W("d_grp.i32", cgrp); W("d_nrow.i32", cnrow);
    W("d_nlep.i32", cnlep); W("d_s1.i32", cs1); W("d_s2.i32", cs2); W("d_sr.i32", csr); W("d_m4.f32", cm4); W("d_mpipi.f32", cpp);
    W("d_mee.f32", cee); W("d_oa.f32", coa); W("d_p1.f32", cp1); W("d_p2.f32", cp2); W("d_t1.f32", ct1); W("d_t2.f32", ct2);
    W("d_f1.f32", cf1); W("d_f2.f32", cf2); W("d_ppip.f32", cpip); W("d_ppim.f32", cpim); W("d_vz.f32", cvz);
    W("d_lpos.i16", lpos); W("d_lneg.i16", lneg); W("d_hpos.i16", hpos); W("d_hneg.i16", hneg); W("d_tmult.i16", tmult);
    W("d_nus.i16", nus); W("d_ntwin.i16", ntwin); W("d_nls.i16", nls); W("d_nflip.i16", nflip); W("d_joinok.u8", joinok);
    FILE* fl = fopen(Form("%s/veestep_layout.txt", kOut), "w");
    if (fl) { fprintf(fl, "data %zu\n", n); fclose(fl); }
    VL("done in %.0f s\n", total.RealTime());
    if (gLog) fclose(gLog);
}
