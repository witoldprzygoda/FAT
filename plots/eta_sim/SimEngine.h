// ============================================================================
// SimEngine.h  --  shared basis of the simulation study of pi+ pi- e+ e-
// ============================================================================
// HADES p+p, T = 4.53 GeV (feb22). SMASH simulation, standard (std), lepton
// (lep) and eta (eta) sets, processed by FAT (branch pp45_pippimepem_betanew_sim).
// Header only. Used by plots/eta_sim/sim_prepare.C (writes the compact files)
// and by every later macro of plots/eta_sim/ (read them, fill weighted spectra,
// compare with the data through ../ScanEngine.h).
//
// Eta set (added 2026-11-12, sim_set = 3): every eta of the same 10^10 SMASH
// events is forced into pi+pi-e+e- or X17 pi+pi- (q = 1/2 each). Its owned rows
// (>= 1 track of a forced eta decay, the same in modes A and B) carry
// W = 2 BR(eta -> pi+pi-e+e-) = 5.36e-4 per forced decay (times the natural
// factors of the other tracks); X17 rows carry W = 0 (BR(X17) = 0 in the
// production) and are recomputed only for the X17 illustration (CompX17).
// Rows without a forced-eta track are not owned (they duplicate std). The eta
// trees live in a second file, sim_compact_eta.root (trees eta_us, eta_pp,
// eta_mm, same columns), written by sim_prepare.C(1). Default: the full
// simulated spectrum = std (owned) + lep (owned) + eta (owned) with one
// L_sim (SimOpt::sets = 7); a macro that loads only the std and lep trees
// gets exactly the numbers of the first study.
//
// Contents
//   A. constants: luminosities, files, kinematics shared with the data
//   B. Row (one stored particle combination of the simulation) and the loader
//      of sim_compact.root (trees std_us, std_pp, std_mm, lep_us, lep_pp, lep_mm)
//      and sim_compact_eta.root (trees eta_us, eta_pp, eta_mm)
//   C. kinematics shared by sim_prepare.C and the data side (pair quantities,
//      crossing points; copies of the definitions in ../scan_prepare.C)
//   D. selections: the named presets for simulation AND data rows
//   E. weighted spectra: sum of weights, event-group variance, N_eff; template
//      components (class filters, optional weight override); std + lep + eta
//      owned rows in mode A or B
//   F. absolute normalisation (nb per bin, nb/(GeV/c2), expected data rows)
//   G. data spectra in the same binning and units (N, CB = 2 sqrt(N++ N--),
//      N - CB) through ScanEngine.h, with L_data as a parameter
//   H. weighted efficiency with event-group errors
//   I. plotting helpers (style, labels)
//
// Rules (binding for this study)
//   - no combinatorial background in simulation: the simulated spectrum used
//     for comparisons = template 1 (signal) + template 2 (physical background);
//     templates 3 (fake lepton) and 4 (leptons from different sources) are
//     shown only as separate, labelled curves
//   - default: rows owned in partition mode A (ownA), std + lep added, weight w
//     (PDG branching ratios); variants: mode B (ownB), weight w_db
//   - sigma per bin = sum w / L_sim, L_sim = 221.05 nb^-1 per set (std, lep and
//     eta describe the same SMASH events: owned rows are added with the same L_sim)
//   - X17 (four-particle class 6 of the eta set) is information only: its
//     shape is shown for an illustrative BR, never added to the simulated sum
//   - data: (N - CB) / L_data / bin width, L_data = 6080.5 nb^-1 (producers'
//     note, UNCONFIRMED: every figure says so)
//   - statistics: per bin S = sum w, V = sum over event groups of (sum w of the
//     group's rows in the bin)^2, N_eff = S^2 / V; std and lep variances add
//     (covariance from the shared SMASH events neglected, spec f)
//   - all weighted sums in double precision (TH1D, double accumulators): never
//     TH1F or TTree::Draw into its default TH1F (3.4 % loss measured)
//   - physics weights of the eta decays (added 2026-11-12, sim_eta_physics.h):
//     SimOpt::phys != 0 multiplies w (w_db) by simphys::RowFactor (direct decay
//     to QED x form factors, eta -> 3pi Dalitz plot, eta -> pi pi gamma form
//     factor); default 0 = OFF, the code path of the earlier steps is unchanged
//
// Conventions as in ScanEngine.h: masses in GeV/c2, momenta in MeV/c, angles
// in deg, lengths in mm. Simulation rows all have cor_ok == 1 (validated), so
// the compact file has no *_fb columns. Rapidities are cms rapidities.
//
// Compile macros that include this header with ACLiC and keep the build
// directory outside the repository, e.g.
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/aclic_x", kTRUE)' \
//        'plots/eta_sim/sim_reference.C+'
// ============================================================================
#ifndef SIM_ENGINE_H
#define SIM_ENGINE_H

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "TFile.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TTree.h"
#include "TVector3.h"

#include "../ScanEngine.h"
#include "../../src/sim_truth.h"

namespace simeng {

// ============================================================================
// A. Constants
// ============================================================================
const double kLsim  = 221.05;    // nb^-1, one simulated set (10^10 SMASH events), src/sim_truth.h kLsimNbInv
const double kLdata = 6080.5;    // nb^-1, data, producers' note, UNCONFIRMED (parameter everywhere)
const char* const kLdataNote = "L_{data} = 6080.5 nb^{-1}, unconfirmed";

const char* const kProdDir = "/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew";
const char* const kScanDir = "/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan";
const char* const kCompact = "/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/sim_compact.root";
const char* const kCompactEta = "/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/sim_compact_eta.root";
const char* const kFigDir  = "/home/przygoda/HADES/FAT/FAT/plots/eta_sim/figs";
const char* const kTabDir  = "/home/przygoda/HADES/FAT/FAT/plots/eta_sim/tables";

// beam: T = 4530 MeV, as scan::kTbeam
const double kMeMeV = 0.51099895, kMpiMeV = 139.56995, kMpMeV = 938.27231, kMetaMeV = 547.862;

// trees of the compact files: index = 3 * (set - 1) + hypothesis.
// kNTree = 6 = the trees of sim_compact.root (std, lep), unchanged for the
// macros of the first study; kNTreeAll = 9 adds the eta trees of
// sim_compact_eta.root.
enum Hyp { kUS = 0, kPP = 1, kMM = 2 };
enum TreeIdx { kStdUS = 0, kStdPP = 1, kStdMM = 2, kLepUS = 3, kLepPP = 4, kLepMM = 5, kNTree = 6,
               kEtaUS = 6, kEtaPP = 7, kEtaMM = 8, kNTreeAll = 9 };
inline int TreeIndex(int set, int hyp) { return 3 * (set - 1) + hyp; }
inline const char* TreeName(int i) {
    static const char* n[kNTreeAll] = {"std_us", "std_pp", "std_mm", "lep_us", "lep_pp", "lep_mm",
                                       "eta_us", "eta_pp", "eta_mm"};
    return (i >= 0 && i < kNTreeAll) ? n[i] : "?";
}
// tree masks for Load (bit i = tree i)
const int kMaskStdLep = 63;                                        // the six trees of sim_compact.root
const int kMaskEta    = (1 << kEtaUS) | (1 << kEtaPP) | (1 << kEtaMM);   // 448
const int kMaskAll    = kMaskStdLep | kMaskEta;                    // 511
const int kMaskUS     = (1 << kStdUS) | (1 << kLepUS) | (1 << kEtaUS);   // e+e- trees of all three sets
// set bits of SimOpt::sets: set s (1 std, 2 lep, 3 eta) -> bit 1 << (s - 1)
enum SetBits { kSetStd = 1, kSetLep = 2, kSetEta = 4, kSetsAll = 7 };
inline int SetBit(int set) { return 1 << (set - 1); }
inline const char* SetName(int set) { return set == 1 ? "std" : set == 2 ? "lep" : set == 3 ? "eta" : "?"; }

// template classes (sim_tmpl) and four-particle classes (sim_cls), src/sim_truth.h
enum Tmpl { kTmplSignal = 1, kTmplPhysBkg = 2, kTmplFake = 3, kTmplDiff = 4 };
inline const char* TmplName(int t) {
    switch (t) {
        case 1: return "signal";
        case 2: return "physical background";
        case 3: return "fake lepton";
        case 4: return "different sources";
    }
    return "?";
}
inline const char* ClsName(int c) {
    static const char* n[9] = {"fake (misidentified)", "fake (pi -> mu)", "bump, pi0 Dalitz", "bump, conversion",
                               "peak, conversion", "peak, eta -> pi+pi-e+e-", "X17", "multi-pion", "CB-like"};
    return (c >= 0 && c < 9) ? n[c] : "?";
}
inline const char* LclsName(int c) {
    static const char* n[12] = {"L0 misidentified", "L1 same pi0 Dalitz", "L2 same eta Dalitz", "L3 same conversion",
                                "L4 rho0/omega/phi -> ee", "L5 omega -> pi0 ee", "L6 second-decay pi0 Dalitz",
                                "L7 two Dalitz", "L8 Dalitz + conversion", "L9 other different", "L10 eta -> pi+pi-ee",
                                "L11 X17"};
    return (c >= 0 && c < 12) ? n[c] : "?";
}
inline const char* PclsName(int c) {
    static const char* n[10] = {"P0 misidentified", "P1 pi -> mu", "P2 same eta", "P3 same omega", "P4 same rho0",
                                "P5 same phi", "P6 same K0S", "P7 both primary", "P8 mixed/secondary", "P9 eta set"};
    return (c >= 0 && c < 10) ? n[c] : "?";
}

// ============================================================================
// B. Row and loader
// ============================================================================
// Names equal to the data compact file (scan::Row) wherever both exist, with the
// same definitions (copied code of ../scan_prepare.C). Added here:
//   keys      set (1 std, 2 lep, 3 eta), event (generator event number), trigbit
//   tracks    X_z, X_r (track point of the straight line, as in scan_nt)
//   pairs     p_pipi, th_pipi, pt_pipi, y_pipi; p_ee, th_ee, pt_ee, y_ee
//             (corrected momenta; lab momentum [MeV/c], lab polar angle [deg],
//             transverse momentum [MeV/c], cms rapidity); mm_pipi, mm_ee as data
//   truth     cls, lcls, pcls, tmpl, ownA, ownB, flags (sim_truth.h kFlag*),
//             ghost (bit 0 ep, 1 em, 2 pip, 3 pim: <x>_sim_isghost),
//             X_sim_id, X_sim_parentid, X_sim_grandparentid, X_sim_processid,
//             X_rich_sim_geanttrack (leptons); w, w_db (Double_t, exact copies
//             of the float sim_w, sim_w_db); m4_true, mee_true, mpipi_true
//             (GeV/c2, -1 if misidentified), oa_true (true e+e- opening angle,
//             deg, -1 if a lepton slot is misidentified), X_sim_p (true momentum),
//             X_sim_theta (true polar angle), X_sim_vz, X_sim_vr (creation
//             vertex z and radius, mm)
#define SIM_INT_FIELDS(X)                                                          \
    X(set) X(event) X(pat_entry) X(trigbit) X(grp) X(fold) X(isBest) X(nrow)       \
    X(n_pip) X(n_pim) X(n_lep1) X(n_lep2) X(n_lep) X(cor_ok)                       \
    X(ep_system) X(em_system) X(pip_system) X(pim_system)                          \
    X(ep_isUsedVertex) X(em_isUsedVertex) X(pip_isUsedVertex) X(pim_isUsedVertex)  \
    X(ep_kIsLepton) X(em_kIsLepton) X(shared_ring)                                 \
    X(pip_isring) X(pim_isring) X(evtStart)                                        \
    X(n_gam) X(g_idx) X(neutr_mult) X(fwdet_mult)                                  \
    X(cls) X(lcls) X(pcls) X(tmpl) X(ownA) X(ownB) X(flags) X(ghost)               \
    X(ep_sim_id) X(em_sim_id) X(pip_sim_id) X(pim_sim_id)                          \
    X(ep_sim_parentid) X(em_sim_parentid) X(pip_sim_parentid) X(pim_sim_parentid)  \
    X(ep_sim_grandparentid) X(em_sim_grandparentid)                                \
    X(pip_sim_grandparentid) X(pim_sim_grandparentid)                              \
    X(ep_sim_processid) X(em_sim_processid) X(pip_sim_processid) X(pim_sim_processid) \
    X(ep_rich_sim_geanttrack) X(em_rich_sim_geanttrack)

#define SIM_FLT_FIELDS(X)                                                          \
    X(m4) X(mee) X(mpipi) X(oa) X(oa_lab) X(oa_eta) X(mm4) X(mm_pipi) X(pt4) X(y4)  \
    X(mm_ee)                                                                       \
    X(p_pipi) X(th_pipi) X(pt_pipi) X(y_pipi) X(p_ee) X(th_ee) X(pt_ee) X(y_ee)    \
    X(vx) X(vy) X(vz) X(vchi2)                                                     \
    X(ep_p) X(em_p) X(pip_p) X(pim_p)                                              \
    X(ep_p_rec) X(em_p_rec) X(pip_p_rec) X(pim_p_rec)                              \
    X(ep_theta) X(em_theta) X(pip_theta) X(pim_theta)                              \
    X(ep_phi) X(em_phi) X(pip_phi) X(pim_phi)                                      \
    X(ep_z) X(em_z) X(pip_z) X(pim_z) X(ep_r) X(em_r) X(pip_r) X(pim_r)            \
    X(ep_rkchi2) X(em_rkchi2) X(pip_rkchi2) X(pim_rkchi2)                          \
    X(ep_richmatchquality) X(em_richmatchquality) X(ep_ringchi2) X(em_ringchi2)    \
    X(ep_oa_lept) X(em_oa_lept) X(ep_oa_hadr) X(em_oa_hadr)                        \
    X(pip_mass2) X(pim_mass2) X(pip_dedx_mdc) X(pim_dedx_mdc)                      \
    X(hyp_chi2)                                                                    \
    X(dz_ee) X(r_ee) X(rv_ee) X(doca_ee)                                           \
    X(dz_pipi) X(r_pipi) X(rv_pipi) X(doca_pipi)                                   \
    X(g_e) X(g_theta) X(g_phi) X(m_eeg) X(m5) X(mm5)                               \
    X(T_pi0) X(m4max)                                                              \
    X(fw1_p) X(fw1_theta)                                                          \
    X(m4_true) X(mee_true) X(mpipi_true) X(oa_true)                                \
    X(ep_sim_p) X(em_sim_p) X(pip_sim_p) X(pim_sim_p)                              \
    X(ep_sim_theta) X(em_sim_theta) X(pip_sim_theta) X(pim_sim_theta)              \
    X(ep_sim_vz) X(em_sim_vz) X(pip_sim_vz) X(pim_sim_vz)                          \
    X(ep_sim_vr) X(em_sim_vr) X(pip_sim_vr) X(pim_sim_vr)

#define SIM_DBL_FIELDS(X) X(w) X(w_db)

struct Row {
#define X(n) Int_t n;
    SIM_INT_FIELDS(X)
#undef X
#define X(n) Float_t n;
    SIM_FLT_FIELDS(X)
#undef X
#define X(n) Double_t n;
    SIM_DBL_FIELDS(X)
#undef X
};

inline void BindRow(TTree* t, Row& r, bool write) {
#define X(n) if (write) t->Branch(#n, &r.n, #n "/I"); else t->SetBranchAddress(#n, &r.n);
    SIM_INT_FIELDS(X)
#undef X
#define X(n) if (write) t->Branch(#n, &r.n, #n "/F"); else t->SetBranchAddress(#n, &r.n);
    SIM_FLT_FIELDS(X)
#undef X
#define X(n) if (write) t->Branch(#n, &r.n, #n "/D"); else t->SetBranchAddress(#n, &r.n);
    SIM_DBL_FIELDS(X)
#undef X
}

// The trees in memory (rows in file order: event groups are adjacent).
// t[0..5] std and lep (sim_compact.root), t[6..8] eta (sim_compact_eta.root).
struct SimData {
    std::vector<Row> t[kNTreeAll];
    int loaded = 0;        // mask of the trees read by Load
    std::vector<Row>& tree(int set, int hyp) { return t[TreeIndex(set, hyp)]; }
    const std::vector<Row>& tree(int set, int hyp) const { return t[TreeIndex(set, hyp)]; }
    // set bits (kSetStd | kSetLep | kSetEta) of the loaded trees
    int sets() const {
        int s = 0;
        for (int k = 0; k < kNTreeAll; ++k) if (loaded & (1 << k)) s |= SetBit(k / 3 + 1);
        return s;
    }
};
// set bits of everything loaded in this process (used by OptLabel so that a
// macro which never loads the eta trees keeps its "std + lep" labels)
inline int& LoadedSetBits() { static int m = 0; return m; }

// Load the compact files.
//   treeMask : bit i = tree i (1 std_us, 2 std_pp, 4 std_mm, 8 lep_us, 16 lep_pp, 32 lep_mm,
//              64 eta_us, 128 eta_pp, 256 eta_mm); kMaskStdLep = 63, kMaskEta = 448,
//              kMaskAll = 511 (default), kMaskUS = 73 (e+e- trees of all sets)
//   path     : file of the std and lep trees; pathEta: file of the eta trees (opened only
//              when the mask asks for an eta tree)
//   mlo, mhi : keep only rows with mlo <= m4 < mhi (whole groups are NOT kept
//              together by this cut; the group variance stays correct because
//              rows outside the range never enter a spectrum of m4 in it; do not
//              use a mass range when histogramming another variable)
//   owned    : 0 keep every row, 1 keep rows owned in mode A or B (the others
//              never enter a spectrum; in the eta set the rows without a
//              forced-eta track are not owned and are dropped)
// Memory: about 650 bytes per row (8.37 M rows in the six std/lep trees, 5.4 GB;
// 3.43 M rows in the eta trees, of them 3.25 M owned, about 2.1 GB).
inline bool Load(SimData& d, const char* path = kCompact, int treeMask = kMaskAll, double mlo = -1e9,
                 double mhi = 1e9, int owned = 1, const char* pathEta = kCompactEta) {
    TStopwatch sw;
    d.loaded = 0;
    for (int part = 0; part < 2; ++part) {
        const int k0 = part == 0 ? 0 : kNTree, k1 = part == 0 ? kNTree : kNTreeAll;
        for (int k = k0; k < k1; ++k) d.t[k].clear();
        const int partMask = part == 0 ? kMaskStdLep : kMaskEta;
        if (!(treeMask & partMask)) continue;
        const char* fn = part == 0 ? path : pathEta;
        TFile* f = TFile::Open(fn, "READ");
        if (!f || f->IsZombie()) { printf("simeng::Load: cannot open %s\n", fn); return false; }
        for (int k = k0; k < k1; ++k) {
            if (!(treeMask & (1 << k))) continue;
            TTree* t = (TTree*)f->Get(TreeName(k));
            if (!t) { printf("simeng::Load: tree %s missing in %s\n", TreeName(k), fn); return false; }
            Row r;
            BindRow(t, r, false);
            const Long64_t n = t->GetEntries();
            d.t[k].reserve(n);
            for (Long64_t i = 0; i < n; ++i) {
                t->GetEntry(i);
                if ((mlo > -1e8 || mhi < 1e8) && !(r.m4 >= mlo && r.m4 < mhi)) continue;   // no range: keep NaN rows too
                if (owned && !r.ownA && !r.ownB) continue;
                d.t[k].push_back(r);
            }
            t->ResetBranchAddresses();
            d.loaded |= 1 << k;
        }
        f->Close();
        delete f;
    }
    LoadedSetBits() |= d.sets();
    printf("simeng::Load:");
    for (int k = 0; k < kNTreeAll; ++k) if (d.loaded & (1 << k)) printf(" %s %zu", TreeName(k), d.t[k].size());
    printf(" rows in %.1f s\n", sw.RealTime());
    return true;
}

// ============================================================================
// C. Kinematics (shared by sim_prepare.C and the data side)
// ============================================================================
const double kD2R = 1.74532925199432955e-02;
const double kR2D = 57.2957795130823229;

inline TLorentzVector P4(double p, double theta, double phi, double mMeV) {
    const double th = theta * kD2R, ph = phi * kD2R;
    TLorentzVector v;
    v.SetVectM(TVector3(p * std::sin(th) * std::cos(ph), p * std::sin(th) * std::sin(ph), p * std::cos(th)), mMeV);
    return v;
}
// initial state and cms rapidity, T = scan::kTbeam
inline TLorentzVector Initial() {
    const double Ebeam = scan::kTbeam + kMpMeV;
    return TLorentzVector(0, 0, std::sqrt(Ebeam * Ebeam - kMpMeV * kMpMeV), Ebeam + kMpMeV);
}
inline double YCM() { return Initial().Rapidity(); }

// Pair quantities of the two subsystems from the corrected track momenta
// (p [MeV/c], theta, phi [deg]); the same function serves simulation rows (in
// sim_prepare.C) and data rows (PairKinOf below), so the definitions agree.
struct PairKin { float p_pipi, th_pipi, pt_pipi, y_pipi, p_ee, th_ee, pt_ee, y_ee; };
inline PairKin ComputePairKin(double pEp, double thEp, double phEp, double pEm, double thEm, double phEm,
                              double pPip, double thPip, double phPip, double pPim, double thPim, double phPim) {
    const double yCM = YCM();
    const TLorentzVector ee = P4(pEp, thEp, phEp, kMeMeV) + P4(pEm, thEm, phEm, kMeMeV);
    const TLorentzVector pp = P4(pPip, thPip, phPip, kMpiMeV) + P4(pPim, thPim, phPim, kMpiMeV);
    PairKin k;
    k.p_pipi = (float)pp.P(); k.th_pipi = (float)(pp.Theta() * kR2D); k.pt_pipi = (float)pp.Pt();
    k.y_pipi = (float)(pp.Rapidity() - yCM);
    k.p_ee = (float)ee.P(); k.th_ee = (float)(ee.Theta() * kR2D); k.pt_ee = (float)ee.Pt();
    k.y_ee = (float)(ee.Rapidity() - yCM);
    return k;
}
// for any row type with ep_p, ep_theta, ep_phi, ... (data: scan::Row, NaN when cor_ok == 0)
template <class R> inline PairKin PairKinOf(const R& r) {
    return ComputePairKin(r.ep_p, r.ep_theta, r.ep_phi, r.em_p, r.em_theta, r.em_phi, r.pip_p, r.pip_theta,
                          r.pip_phi, r.pim_p, r.pim_theta, r.pim_phi);
}

// Point of closest approach of two track lines (copy of prep::Crossing of
// ../scan_prepare.C; HYDRA convention base = (-r sin phi, r cos phi, z)).
struct CrossPt { float dz, r, rv, doca; };
inline CrossPt Crossing(double r1, double z1, double th1, double ph1, double r2, double z2, double th2, double ph2,
                        double vx, double vy, double vz) {
    const double t1 = th1 * kD2R, p1 = ph1 * kD2R, t2 = th2 * kD2R, p2 = ph2 * kD2R;
    const TVector3 b1(-r1 * std::sin(p1), r1 * std::cos(p1), z1);
    const TVector3 b2(-r2 * std::sin(p2), r2 * std::cos(p2), z2);
    const TVector3 d1(std::sin(t1) * std::cos(p1), std::sin(t1) * std::sin(p1), std::cos(t1));
    const TVector3 d2(std::sin(t2) * std::cos(p2), std::sin(t2) * std::sin(p2), std::cos(t2));
    const TVector3 w = b1 - b2;
    const double b = d1.Dot(d2), d = d1.Dot(w), e = d2.Dot(w);
    const double den = 1.0 - b * b;
    CrossPt c;
    if (den < 1e-12) {
        c.dz = -9999.f; c.r = -1.f; c.rv = -1.f;
        c.doca = (float)(w - d * d1).Mag();
        return c;
    }
    const double s = (b * e - d) / den, t = (e - b * d) / den;
    const TVector3 q1 = b1 + s * d1, q2 = b2 + t * d2;
    const TVector3 q = 0.5 * (q1 + q2);
    c.dz = (float)(q.Z() - vz);
    c.r = (float)std::hypot(q.X(), q.Y());
    c.rv = (float)std::hypot(q.X() - vx, q.Y() - vy);
    c.doca = (float)(q1 - q2).Mag();
    return c;
}

// ============================================================================
// D. Selections (named presets, for simulation rows AND data rows)
// ============================================================================
// "pre" = isBest >= 0, vertex z > -500 mm, PT3, vertex z > -160 mm, corrected
// momenta (cor_ok == 1). Data: the production keeps only isBest >= 0,
// z > -500 and trigbit == 8192 rows; cor_ok == 1 is applied here. Simulation:
// PT3 = (trigbit & 8192) != 0 (sim PT3 rows carry 12288); cor_ok == 1 in every
// simulated row. Every other preset is "pre" plus its cuts.
enum Preset { kPre = 0, kPreOA2, kPrePiVtx, kSlides, kOldCode, kPeakNom, kBumpNom, kNPreset };
inline const char* PresetName(int p) {
    static const char* n[kNPreset] = {"pre", "pre+OA>2", "pre+pions in vertex fit", "slides", "old code",
                                      "peak nominal", "bump nominal"};
    return (p >= 0 && p < kNPreset) ? n[p] : "?";
}
// short label for figures (ROOT TLatex)
inline const char* PresetLabel(int p) {
    static const char* n[kNPreset] = {
        "pre: isBest #geq 0, PT3, v_{z} > -160 mm",
        "pre + OA > 2#circ",
        "pre + both pions in the vertex fit",
        "slides: pre + #alpha_{lab} < 50#circ, M_{#pi#pi} < 0.46, #alpha_{#eta} > 140#circ, OA > 2#circ",
        "old code: pre + #alpha_{lab} < 50#circ, M_{#pi#pi} < 0.42, #alpha_{#eta} > 140#circ",
        "peak nominal: pre + #alpha_{lab} < 55#circ, M_{#pi#pi} > 0.30",
        "bump nominal: pre + #alpha_{lab} < 40#circ, 0.30 < M_{#pi#pi} < 0.41"};
    return (p >= 0 && p < kNPreset) ? n[p] : "?";
}
inline bool PT3(const Row& r) { return (r.trigbit & 8192) != 0; }
inline bool PT3(const scan::Row&) { return true; }          // data production: trigbit == 8192 only
template <class R> inline bool PassPre(const R& r) {
    return r.isBest >= 0 && r.vz > -500.f && r.vz > -160.f && r.cor_ok == 1 && PT3(r);
}
template <class R> inline bool PassPreset(int p, const R& r) {
    if (!PassPre(r)) return false;
    switch (p) {
        case kPre:      return true;
        case kPreOA2:   return r.oa > 2.f;
        case kPrePiVtx: return r.pip_isUsedVertex == 1 && r.pim_isUsedVertex == 1;
        case kSlides:   return r.oa_lab < 50.f && r.mpipi < 0.46f && r.oa_eta > 140.f && r.oa > 2.f;
        case kOldCode:  return r.oa_lab < 50.f && r.mpipi < 0.42f && r.oa_eta > 140.f;
        case kPeakNom:  return r.oa_lab < 55.f && r.mpipi > 0.30f;
        case kBumpNom:  return r.oa_lab < 40.f && r.mpipi > 0.30f && r.mpipi < 0.41f;
    }
    return false;
}

typedef std::function<bool(const Row&)> Sel;
typedef std::function<double(const Row&)> Var;
inline Sel SelPreset(int p) { return [p](const Row& r) { return PassPreset(p, r); }; }
inline scan::Sel DataSelPreset(int p) { return [p](const scan::Row& r) { return PassPreset(p, r); }; }
inline Sel SelAnd(const Sel& a, const Sel& b) { return [a, b](const Row& r) { return a(r) && b(r); }; }
inline Sel SelBest(const Sel& s) { return [s](const Row& r) { return r.isBest == 1 && s(r); }; }

// ============================================================================
// E. Weighted spectra
// ============================================================================
// A component = a class filter on the truth columns (applied after the
// ownership of the chosen mode), with an optional weight override wf (used
// instead of w / w_db when set; only the X17 illustration uses it).
struct Comp {
    std::string name;
    std::function<bool(const Row&)> f;
    std::function<double(const Row&)> wf;
    Comp() {}
    Comp(const std::string& n, const std::function<bool(const Row&)>& ff,
         const std::function<double(const Row&)>& w = std::function<double(const Row&)>())
        : name(n), f(ff), wf(w) {}
};
inline bool GhostLepton(const Row& r) { return (r.flags & simtruth::kFlagGhostLepton) != 0; }
// X17 signal rows (eta set, four-particle class 6): W = 0 in the production
// (BR(X17) = 0); information only, never part of the simulated spectrum
inline bool X17Row(const Row& r) { return r.cls == 6; }
inline Comp CompAny()        { return {"all rows", [](const Row&) { return true; }}; }
// the simulated spectrum used for comparisons: signal + physical background
// (X17 rows excluded explicitly; they carry W = 0 anyway, so only the row count changes)
inline Comp CompSim()        { return {"signal + physical background", [](const Row& r) { return (r.tmpl == 1 && !X17Row(r)) || r.tmpl == 2; }}; }
inline Comp CompTmpl(int t)  { return {TmplName(t), [t](const Row& r) { return r.tmpl == t; }}; }
inline Comp CompSignal()     { return CompTmpl(kTmplSignal); }
inline Comp CompPhysBkg()    { return CompTmpl(kTmplPhysBkg); }
inline Comp CompFake()       { return CompTmpl(kTmplFake); }
inline Comp CompDiff()       { return CompTmpl(kTmplDiff); }
// signal channels (template 1 split by the four-particle class); noGhost: without ghost-lepton rows
inline Comp CompSigCls(int c, bool noGhost = false) {
    return {std::string(ClsName(c)) + (noGhost ? ", no ghost lepton" : ""),
            [c, noGhost](const Row& r) { return r.tmpl == 1 && r.cls == c && !(noGhost && GhostLepton(r)); }};
}
inline Comp CompBumpDalitz() { return CompSigCls(2); }
inline Comp CompBumpConv()   { return CompSigCls(3); }
inline Comp CompPeakConv()   { return CompSigCls(4); }
// physical background split by the pion-pair class
inline Comp CompBkgPcls(int p) {
    return {std::string("phys. bkg, ") + PclsName(p), [p](const Row& r) { return r.tmpl == 2 && r.pcls == p; }};
}
// any four-particle class (all templates)
inline Comp CompCls(int c) { return {ClsName(c), [c](const Row& r) { return r.cls == c; }}; }

// ---- eta-set components ----
// peak, direct decay eta -> pi+pi-e+e- (four-particle class 5, eta set, W = 2 x 2.68e-4)
inline Comp CompPeakDirect() { return CompSigCls(5); }
// peak total = conversion (eta -> pi+pi-gamma, std) + direct (eta -> pi+pi-e+e-, eta set)
inline Comp CompPeakTotal() {
    return {"peak, conversion + direct", [](const Row& r) { return r.tmpl == 1 && (r.cls == 4 || r.cls == 5); }};
}
// X17 signal (four-particle class 6) with W = 2 BR(eta -> pi+pi- X17, X17 -> e+e-): exact for class 6, where all
// four tracks come from one forced decay (no natural factor). INFORMATION ONLY: illustrative BR (default 1e-5, the
// order of the BESIII limit 1.5e-5); never add it to the simulated spectrum.
inline Comp CompX17(double br = 1e-5) {
    Comp c{Form("X17 signal, BR = %.2g (illustration)", br), [](const Row& r) { return r.tmpl == 1 && X17Row(r); }};
    c.wf = [br](const Row&) { return 2. * br; };
    return c;
}
// restrict a component to the sets of a mask (kSetStd | kSetLep | kSetEta)
inline Comp CompInSets(const Comp& c, int setMask) {
    Comp o{c.name + (setMask == kSetEta ? " (eta set)" : setMask == kSetStd ? " (std set)" : setMask == kSetLep ? " (lep set)"
                     : setMask == (kSetStd | kSetLep) ? " (std + lep)" : ""),
           [c, setMask](const Row& r) { return (SetBit(r.set) & setMask) != 0 && c.f(r); }};
    o.wf = c.wf;
    return o;
}
// physical background (template 2) of the sets in the mask
inline Comp CompPhysBkgSets(int setMask) { return CompInSets(CompPhysBkg(), setMask); }
// physical background of the eta set: template 2 of set 3 (four-particle class 7, and the classes 0 / 1 with a true
// pair and a misidentified or decayed pion); rows with W = 0 (X17 products) carry no weight
inline Comp CompEtaBkg() { return CompPhysBkgSets(kSetEta); }

struct SimOpt {
    int mode = 0;          // 0: partition mode A (ownA), 1: mode B (ownB); the eta set owns the same rows in both
    bool db = false;       // false: w (PDG BR); true: w_db (converter-database BR, no natural-decay factors)
    int sets = kSetsAll;   // bit 0 std, bit 1 lep, bit 2 eta (owned rows of all three added; a set whose trees are
                           // not loaded contributes nothing, so std + lep macros are unchanged)
    int hyp = kUS;         // kUS, kPP, kMM
    int phys = 0;          // physics weights of the eta decays (sim_eta_physics.h, simphys::Code / kPhysMain):
                           // 0 = OFF (default; every result of the earlier steps unchanged); != 0 multiplies w
                           // (or w_db) by simphys::RowFactor(row, phys); never applied to a weight override (Comp::wf)
};
// Physics-weight hook (added 2026-11-12): set by simphys::Install() of sim_eta_physics.h. Used only when
// SimOpt::phys != 0; a macro that never sets phys runs exactly the code of the earlier steps.
typedef double (*PhysFn)(const Row&, int);
inline PhysFn& PhysHook() { static PhysFn f = nullptr; return f; }
inline double PhysFactor(const Row& r, int phys) {
    if (!PhysHook()) {
        static bool warned = false;
        if (!warned) { printf("simeng: SimOpt::phys = %d but no physics weights installed (simphys::Init + Install): factor 1\n", phys); warned = true; }
        return 1.;
    }
    return PhysHook()(r, phys);
}
inline bool Owned(const Row& r, const SimOpt& o) { return o.mode == 0 ? r.ownA == 1 : r.ownB == 1; }
inline double Weight(const Row& r, const SimOpt& o) {
    if (!Owned(r, o)) return 0.;
    const double w = o.db ? r.w_db : r.w;
    return o.phys ? w * PhysFactor(r, o.phys) : w;
}
// label of the summed sets: the sets of o.sets that were loaded in this process
inline std::string OptLabel(const SimOpt& o) {
    const int s = LoadedSetBits() ? (o.sets & LoadedSetBits()) : (o.sets & kSetsAll);
    std::string n;
    for (int k = 1; k <= 3; ++k)
        if (s & SetBit(k)) n += std::string(n.empty() ? "" : " + ") + SetName(k);
    if (n.empty()) n = "none";
    return n + ", mode " + (o.mode ? "B" : "A") + (o.db ? ", w_{db}" : ", w (PDG)") + (o.phys ? Form(", phys %d", o.phys) : "");
}

// binning of the histogrammed variable
struct Binning {
    int nbin = 160;
    double lo = 0.3, hi = 1.1;
    Var x;                 // empty = Row::m4
    double bw() const { return (hi - lo) / nbin; }
};

// Weighted spectrum of one component.
//   s[i]  sum of weights in bin i;  v[i] event-group variance;  n[i] rows (owned, w > 0 or not)
//   S, V, N: inside [lo, hi);  Sall, Vall, Nall: all selected rows (any x, NaN x included)
struct WSpec {
    std::string name;
    int nbin = 0;
    double lo = 0, hi = 0, bw = 0;
    std::vector<double> s, v, n;
    double S = 0, V = 0, N = 0, Sall = 0, Vall = 0, Nall = 0;
    double x(int i) const { return lo + (i + 0.5) * bw; }
    int bin(double m) const { return (m >= lo && m < hi) ? std::min(nbin - 1, (int)((m - lo) / bw)) : -1; }
    double Neff(int i) const { return v[i] > 0 ? s[i] * s[i] / v[i] : 0.; }
    double NeffIn() const { return V > 0 ? S * S / V : 0.; }
    double NeffAll() const { return Vall > 0 ? Sall * Sall / Vall : 0.; }
    // sum and variance in [a, b) (whole bins whose centre is inside; group correlations
    // between bins are NOT included: use FillMany with the range as binning for exact sums)
    std::pair<double, double> Window(double a, double b) const {
        double ss = 0, vv = 0;
        for (int i = 0; i < nbin; ++i) if (x(i) >= a && x(i) < b) { ss += s[i]; vv += v[i]; }
        return std::make_pair(ss, vv);
    }
    void Add(const WSpec& o) {   // independent samples (std and lep trees): sums and variances add
        for (int i = 0; i < nbin; ++i) { s[i] += o.s[i]; v[i] += o.v[i]; n[i] += o.n[i]; }
        S += o.S; V += o.V; N += o.N; Sall += o.Sall; Vall += o.Vall; Nall += o.Nall;
    }
    void Scale(double f) {
        for (int i = 0; i < nbin; ++i) { s[i] *= f; v[i] *= f * f; }
        S *= f; V *= f * f; Sall *= f; Vall *= f * f;
    }
    // unit: 0 sum of weights; 1 sigma per bin [nb] = S / L; 2 dsigma/dx [nb / unit of x] = S / L / bw;
    //       3 expected data rows per bin = S * Ldata / L
    TH1D* Hist(const char* hname, int unit = 2, double L = kLsim, double Ldata = kLdata) const {
        TH1D* h = new TH1D(hname, "", nbin, lo, hi);
        h->SetDirectory(nullptr);
        h->Sumw2();
        const double f = unit == 0 ? 1. : unit == 1 ? 1. / L : unit == 2 ? 1. / (L * bw) : Ldata / L;
        for (int i = 0; i < nbin; ++i) {
            h->SetBinContent(i + 1, s[i] * f);
            h->SetBinError(i + 1, std::sqrt(v[i]) * f);
        }
        return h;
    }
};

namespace detail {
// add the squared group sums of one group to the variance (entries: bin, weight; bin -1 = outside)
inline void FlushGroup(std::vector<std::pair<int, double>>& e, WSpec& sp) {
    if (e.empty()) return;
    std::sort(e.begin(), e.end());
    double all = 0, in = 0;
    size_t i = 0;
    while (i < e.size()) {
        const int b = e[i].first; double acc = 0;
        while (i < e.size() && e[i].first == b) { acc += e[i].second; ++i; }
        all += acc;
        if (b >= 0) { sp.v[b] += acc * acc; in += acc; }
    }
    sp.Vall += all * all;
    sp.V += in * in;
    e.clear();
}
}  // namespace detail

// Fill the weighted spectra of several components for one selection, in ONE
// pass over the trees of opt.sets and opt.hyp. Rows: selection true, owned in
// opt.mode; weight w (or w_db, or the component's wf). Variance per event
// group (grp of the tree). Sets 1 std, 2 lep, 3 eta (bits 1, 2, 4 of opt.sets).
inline std::vector<WSpec> FillMany(const SimData& d, const Sel& sel, const std::vector<Comp>& comps,
                                   const SimOpt& opt = SimOpt(), const Binning& b = Binning()) {
    const int nc = (int)comps.size();
    std::vector<WSpec> out(nc);
    for (int c = 0; c < nc; ++c) {
        WSpec& sp = out[c];
        sp.name = comps[c].name; sp.nbin = b.nbin; sp.lo = b.lo; sp.hi = b.hi; sp.bw = b.bw();
        sp.s.assign(b.nbin, 0.); sp.v.assign(b.nbin, 0.); sp.n.assign(b.nbin, 0.);
    }
    const bool useX = (bool)b.x;
    std::vector<std::vector<std::pair<int, double>>> buf(nc);
    for (int set = 1; set <= 3; ++set) {
        if (!(opt.sets & SetBit(set))) continue;
        const std::vector<Row>& v = d.tree(set, opt.hyp);
        int cur = -1;
        for (size_t i = 0; i < v.size(); ++i) {
            const Row& r = v[i];
            if (r.grp != cur) {
                for (int c = 0; c < nc; ++c) detail::FlushGroup(buf[c], out[c]);
                cur = r.grp;
            }
            if (!Owned(r, opt)) continue;
            if (!sel(r)) continue;
            double w0 = opt.db ? r.w_db : r.w;
            if (opt.phys) w0 *= PhysFactor(r, opt.phys);
            const double xv = useX ? b.x(r) : (double)r.m4;
            const int ib = (xv == xv) ? out[0].bin(xv) : -1;
            for (int c = 0; c < nc; ++c) {
                if (!comps[c].f(r)) continue;
                const double w = comps[c].wf ? comps[c].wf(r) : w0;
                WSpec& sp = out[c];
                sp.Sall += w; sp.Nall += 1;
                if (ib >= 0) { sp.s[ib] += w; sp.n[ib] += 1; sp.S += w; sp.N += 1; }
                buf[c].push_back(std::make_pair(ib, w));
            }
        }
        for (int c = 0; c < nc; ++c) detail::FlushGroup(buf[c], out[c]);
    }
    return out;
}
inline WSpec Fill(const SimData& d, const Sel& sel, const Comp& comp, const SimOpt& opt = SimOpt(),
                  const Binning& b = Binning()) {
    return FillMany(d, sel, std::vector<Comp>(1, comp), opt, b)[0];
}

// Totals (sum w, group variance, rows) of several components, without binning.
struct WSum {
    double S = 0, V = 0, N = 0;
    double Neff() const { return V > 0 ? S * S / V : 0.; }
    double Xsec(double L = kLsim) const { return S / L; }            // nb
    double XsecErr(double L = kLsim) const { return std::sqrt(V) / L; }
    double DataRows(double Ldata = kLdata, double L = kLsim) const { return S * Ldata / L; }
};
inline std::vector<WSum> SumMany(const SimData& d, const Sel& sel, const std::vector<Comp>& comps,
                                 const SimOpt& opt = SimOpt()) {
    Binning b; b.nbin = 1; b.lo = -1e30; b.hi = 1e30;
    b.x = [](const Row&) { return 0.; };
    const std::vector<WSpec> sp = FillMany(d, sel, comps, opt, b);
    std::vector<WSum> out(sp.size());
    for (size_t i = 0; i < sp.size(); ++i) { out[i].S = sp[i].Sall; out[i].V = sp[i].Vall; out[i].N = sp[i].Nall; }
    return out;
}

// ============================================================================
// F. Absolute normalisation
// ============================================================================
inline double XsecPerBin(double sumw, double L = kLsim) { return sumw / L; }                       // nb
inline double DsDx(double sumw, double bw, double L = kLsim) { return sumw / (L * bw); }            // nb / unit of x
inline double ExpectedDataRows(double sumw, double Ldata = kLdata, double L = kLsim) { return sumw * Ldata / L; }

// ============================================================================
// G. Data spectra in the same binning and units (through ScanEngine.h)
// ============================================================================
// N = unlike-sign rows, CB = 2 sqrt(N++ N--) (in wide bins of reb fine bins,
// interpolated to the fine bins exactly as ScanEngine does), N - CB. Errors:
// sqrt(N) for N, the like-sign Poisson propagation for CB (ScanEngine). Rows
// with cor_ok == 0 never enter (the presets require cor_ok == 1 anyway).
struct DSpec {
    int nbin = 0;
    double lo = 0, hi = 0, bw = 0, Ldata = kLdata;
    int reb = 1;
    scan::Spectra sp;      // the ScanEngine spectra (rows)
    double x(int i) const { return lo + (i + 0.5) * bw; }
    // what: "us" (N), "pp", "mm", "cb", "sig" (N - CB);
    // unit: 0 rows; 2 dsigma/dx [nb / unit of x] = rows / Ldata / bw; 1 sigma per bin [nb] = rows / Ldata
    TH1D* Hist(const char* what, const char* hname, int unit = 2) const {
        TH1D* h = sp.Hist(what, hname);
        h->Sumw2();
        if (unit != 0) {
            const double f = unit == 1 ? 1. / Ldata : 1. / (Ldata * bw);
            for (int i = 1; i <= h->GetNbinsX(); ++i) {
                h->SetBinContent(i, h->GetBinContent(i) * f);
                h->SetBinError(i, h->GetBinError(i) * f);
            }
        }
        return h;
    }
};
// largest divisor of nbin whose wide bin is at most 10 MeV (M4) ; 1 for other variables
inline int AutoReb(int nbin, double bw, bool isM4) {
    if (!isM4) return 1;
    int best = 1;
    for (int r = 1; r <= nbin; ++r) if (nbin % r == 0 && r * bw <= 0.0100001) best = r;
    return best;
}
inline DSpec DataSpectrum(const scan::Data& d, const scan::Sel& sel, int nbin, double lo, double hi,
                          scan::Var x = scan::Var(), double Ldata = kLdata, int reb = -1) {
    DSpec o;
    o.nbin = nbin; o.lo = lo; o.hi = hi; o.bw = (hi - lo) / nbin; o.Ldata = Ldata;
    o.reb = reb > 0 ? reb : AutoReb(nbin, o.bw, !(bool)x);
    if (nbin % o.reb != 0) o.reb = 1;            // ScanEngine needs nbin divisible by reb
    scan::FillOpt fo;
    fo.nbin = nbin; fo.lo = lo; fo.hi = hi; fo.reb = o.reb; fo.x = x;
    o.sp = scan::Fill(d, sel, fo);
    return o;
}
inline DSpec DataSpectrumPreset(const scan::Data& d, int preset, int nbin, double lo, double hi,
                                double Ldata = kLdata, int reb = -1) {
    return DataSpectrum(d, DataSelPreset(preset), nbin, lo, hi, scan::Var(), Ldata, reb);
}

// ============================================================================
// H. Weighted efficiency with event-group errors
// ============================================================================
// eps = S_pass / S_all over the rows of `base` (owned, component comp);
// variance from the event groups (pass and fail rows of one group correlated):
//   V(eps) = sum_g [(1 - eps) a_g - eps b_g]^2 / S_all^2, a_g (b_g) = weight of the
// group's rows that pass (fail). Without group correlation this is the formula
// of the spec, [(1-eps)^2 V_pass + eps^2 V_fail] / S_all^2 (EffFromSums).
struct EffRes {
    double eps = 0, err = 0, Spass = 0, Sall = 0, Npass = 0, Nall = 0, NeffAll = 0;
};
inline EffRes EffFromSums(double Sp, double Vp, double Sf, double Vf) {
    EffRes e; e.Spass = Sp; e.Sall = Sp + Sf;
    if (e.Sall > 0) {
        e.eps = Sp / e.Sall;
        e.err = std::sqrt((1 - e.eps) * (1 - e.eps) * Vp + e.eps * e.eps * Vf) / e.Sall;
    }
    return e;
}
// per bin of b.x (or m4) and in total (last element = all rows of base, any x)
inline std::vector<EffRes> EfficiencyVs(const SimData& d, const Sel& base, const Sel& pass, const Comp& comp,
                                        const SimOpt& opt = SimOpt(), const Binning& b = Binning()) {
    const int nb = b.nbin;
    std::vector<EffRes> out(nb + 1);
    std::vector<double> sp(nb + 1, 0.), sa(nb + 1, 0.), np(nb + 1, 0.), na(nb + 1, 0.), va(nb + 1, 0.);
    // pass 1: sums
    for (int set = 1; set <= 3; ++set) {
        if (!(opt.sets & SetBit(set))) continue;
        for (const Row& r : d.tree(set, opt.hyp)) {
            if (!Owned(r, opt) || !base(r) || !comp.f(r)) continue;
            const double w = comp.wf ? comp.wf(r) : (opt.db ? r.w_db : r.w) * (opt.phys ? PhysFactor(r, opt.phys) : 1.);
            const double xv = b.x ? b.x(r) : (double)r.m4;
            const int ib = (xv >= b.lo && xv < b.hi) ? std::min(nb - 1, (int)((xv - b.lo) / b.bw())) : -1;
            const bool ps = pass(r);
            for (int k : {ib, nb}) {
                if (k < 0) continue;
                sa[k] += w; na[k] += 1;
                if (ps) { sp[k] += w; np[k] += 1; }
            }
        }
    }
    std::vector<double> eps(nb + 1, 0.);
    for (int k = 0; k <= nb; ++k) eps[k] = sa[k] > 0 ? sp[k] / sa[k] : 0.;
    // pass 2: group variance of (1 - eps) a - eps b, and of the total weight
    std::vector<double> V(nb + 1, 0.);
    for (int set = 1; set <= 3; ++set) {
        if (!(opt.sets & SetBit(set))) continue;
        const std::vector<Row>& v = d.tree(set, opt.hyp);
        std::vector<std::pair<int, std::pair<double, double>>> g;   // bin, (pass w, fail w)
        auto flush = [&]() {
            if (g.empty()) return;
            std::sort(g.begin(), g.end());
            double ta = 0, tb = 0;
            size_t i = 0;
            while (i < g.size()) {
                const int k = g[i].first; double a = 0, bb = 0;
                while (i < g.size() && g[i].first == k) { a += g[i].second.first; bb += g[i].second.second; ++i; }
                ta += a; tb += bb;
                if (k >= 0) {
                    const double q = (1 - eps[k]) * a - eps[k] * bb;
                    V[k] += q * q; va[k] += (a + bb) * (a + bb);
                }
            }
            const double q = (1 - eps[nb]) * ta - eps[nb] * tb;
            V[nb] += q * q; va[nb] += (ta + tb) * (ta + tb);
            g.clear();
        };
        int cur = -1;
        for (const Row& r : v) {
            if (r.grp != cur) { flush(); cur = r.grp; }
            if (!Owned(r, opt) || !base(r) || !comp.f(r)) continue;
            const double w = comp.wf ? comp.wf(r) : (opt.db ? r.w_db : r.w) * (opt.phys ? PhysFactor(r, opt.phys) : 1.);
            const double xv = b.x ? b.x(r) : (double)r.m4;
            const int ib = (xv >= b.lo && xv < b.hi) ? std::min(nb - 1, (int)((xv - b.lo) / b.bw())) : -1;
            const bool ps = pass(r);
            g.push_back(std::make_pair(ib, std::make_pair(ps ? w : 0., ps ? 0. : w)));
        }
        flush();
    }
    for (int k = 0; k <= nb; ++k) {
        EffRes& e = out[k];
        e.Spass = sp[k]; e.Sall = sa[k]; e.Npass = np[k]; e.Nall = na[k]; e.eps = eps[k];
        e.err = sa[k] > 0 ? std::sqrt(V[k]) / sa[k] : 0.;
        e.NeffAll = va[k] > 0 ? sa[k] * sa[k] / va[k] : 0.;
    }
    return out;
}
inline EffRes Efficiency(const SimData& d, const Sel& base, const Sel& pass, const Comp& comp,
                         const SimOpt& opt = SimOpt()) {
    Binning b; b.nbin = 1; b.lo = -1e30; b.hi = 1e30; b.x = [](const Row&) { return 0.; };
    return EfficiencyVs(d, base, pass, comp, opt, b).back();
}

// ============================================================================
// I. Plot helpers
// ============================================================================
// colours of the simulated components (data: black N, red CB, blue N - CB)
const int kColBumpDalitz = 807;   // kOrange+7
const int kColBumpConv   = 617;   // kMagenta+1
const int kColPeakConv   = 418;   // kGreen+2
const int kColPhysBkg    = 922;   // kGray+2 (line); fill kGray = 920
const int kColPhysFill   = 920;   // kGray
const int kColSim        = 843;   // kTeal+3: signal + physical background (sum)
const int kColFake       = 861;   // kAzure+1, dashed
const int kColDiff       = 403;   // kYellow+3, dashed
const int kColPeakDirect = 601;   // kBlue+1: peak, eta -> pi+pi-e+e- (eta set)
const int kColPeakTotal  = 435;   // kCyan+3: peak, conversion + direct
const int kColX17        = 632;   // kRed: X17 illustration, dashed, never stacked
const int kColEtaBkg     = 882;   // kViolet+2: physical background of the eta set
const int kColEtaFill    = 877;   // kViolet-3
// standard caption lines
// (first study, std + lep only; macros that load the eta set use NoteWithEta)
inline const char* NoteNoDirect() { return "#eta #rightarrow #pi^{+}#pi^{-}e^{+}e^{-} (direct) NOT in the simulation (eta set missing)"; }
inline const char* NoteWithEta() {
    return "#eta #rightarrow #pi^{+}#pi^{-}e^{+}e^{-} (direct) from the eta set, BR = 2.68#times10^{-4} (PDG)";
}
inline std::string NoteX17(double br = 1e-5) {
    return Form("X17: shape for an illustrative BR(#eta #rightarrow #pi^{+}#pi^{-}X17, X17 #rightarrow e^{+}e^{-}) = %.1g "
                "(BESIII limit 1.5#times10^{-5}); information only, not in the sum", br);
}
inline std::string NoteSimNorm() {
    return Form("SMASH std + lep, owned rows; #sigma = #Sigma w / L_{sim}, L_{sim} = %.2f nb^{-1}; no efficiency correction", kLsim);
}
inline std::string NoteSimNormEta() {
    return Form("SMASH std + lep + eta, owned rows (mode A); #sigma = #Sigma w / L_{sim}, L_{sim} = %.2f nb^{-1} per set; "
                "no efficiency correction", kLsim);
}

}  // namespace simeng

#endif  // SIM_ENGINE_H
