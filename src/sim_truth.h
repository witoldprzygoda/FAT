/**
 * @file sim_truth.h
 * @brief Truth classes, set ownership and physical weights of the SMASH simulation rows
 *
 * Rules version 1 (sim_weights_spec.md of 2026-11-12 and its reference code
 * simrules.h). The classes, ownership and weights are numerically identical to
 * simrules.h; this file adds a generator-code dictionary (with recovery of the
 * float-rounded large codes), a flags word per row, the template class of the
 * final specification (sim_spec_final.md: signal / physical background / fake
 * lepton / different sources), the true masses and one call that evaluates a
 * whole row. The additions never change a class, an ownership or a weight.
 *
 * Header only, no ROOT, no structured bindings: used by FAT (main.cc) and by
 * ROOT macros. Codes are compared as the Float_t values stored by PAT.
 *
 * Slots: t[0], t[1] = the two lepton slots (ep/em, ep1/ep2 or em1/em2),
 *        t[2] = pi+ slot, t[3] = pi- slot.
 * Sets:  1 = standard (no forced decays, all weights 0),
 *        2 = lepton (PLUTO forces pi0, eta, rho0, omega, phi into dileptons),
 *        3 = eta (every eta forced into pi+pi-e+e- or X17 pi+pi-).
 *
 * @author Witold Przygoda (witold.przygoda@uj.edu.pl)
 * @date 2026
 */

#ifndef SIM_TRUTH_H
#define SIM_TRUTH_H

#include <cmath>
#include <set>
#include <map>
#include <vector>
#include <utility>
#include <limits>

namespace simtruth {

const int kRulesVersion = 2;   // 2 (2026-11-13): exact HGeant branching ratios from the ETA card; 1: measured values
enum Set { kStd = 1, kLep = 2, kEta = 3 };

// Integrated luminosity of one simulated set (10^10 SMASH events, b_max = 1.2 fm):
// N_exp = (L_data / kLsimNbInv) * sum W. L_data is a parameter of the analysis.
const double kLsimNbInv = 221.05;   // nb^-1

// One track of a row: GEANT id, parent id, grandparent id, creation process,
// generator codes (geninfo, geninfo1, geninfo2), stored generator weight,
// creation vertex [mm] and true momentum [MeV/c]. All as stored by PAT.
// ghost = <x>_sim_isghost (1: PAT marks the track as a ghost). It is not used by
// the rules-v1 classes and weights, only by RowFlags and TemplateClass; an
// initialiser with the first 14 members (as in simrules.h) leaves it 0.
struct Trk { int id, par, gp, proc; float g0, g1, g2, w, vx, vy, vz, px, py, pz; int ghost; };

struct Opt {
  bool pdg = true;              // true: physical PDG branching ratios; false: converter database values (the note's weights)
  double brEtaPPee = 2.68e-4;   // eta -> pi+ pi- e+ e- (PDG)
  double brX17 = 0.;            // eta -> X17 pi+ pi- hypothesis; 0 = X17 rows get W = 0
  int mode = 0;                 // 0 = A (lepton set provides pi0/eta Dalitz rows), 1 = B (standard set keeps them;
                                //     lepton set provides only rows with a rho0/omega/phi forced product)
};

// ---------------- constants ----------------
// HGeant branching ratios of the pi0 and eta decays performed by HGeant (natural decays), EXACT values from the
// HGeant configuration of this production (rules v2, 2026-11-13; rules v1 used the producers' measurement on 1e6
// events: pi0 Dalitz 0.01202, eta 0.3769 / 0.3014 / 0.2074 / 0.0497 / 0.0640):
//  - every HGeant input file of the feb22 p4500p SMASH productions (standard, lepton, eta), e.g.
//    /lustre/hades/user/hadesdst/sub/feb22/p4500p/hgeant/p4500p/smash_standard/hgeant/p4500p.dat, contains the card
//        ETA  37.30 30.20 21.00 4.75 6.0 7.1e-2
//    = eta branchings in percent, in the order (hgeant2-7.4 geant/hgeakey.F: ffkey('ETA',breta,6); geant/uginit.F:
//    dmode1 = 101, 70707, 80907, 80901, 30201, 10107): gg, 3pi0, pi+pi-pi0, pi+pi-gamma, e+e-gamma, gg pi0;
//    HGeant normalises them to their sum (normb), 99.321. The HGeant default is 39.30 32.20 23.00 4.75 0.60 0.071:
//    the card raises e+e-gamma tenfold (6.0 instead of 0.60) and lowers the first three by 2.00 each.
//  - no PI0 card: HGeant default bratio0 = 98.798, 1.198 (gg, e+e-gamma), normalised to 99.996.
// The production ran hgeant2-7.5; its source was not found, the order above is that of 7.4 and agrees with the
// producers' measurement. PDG values as quoted in the note (RPP 2020-2024).
const double kHG_pi0Sum = 98.798 + 1.198;
const double kHG_etaSum = 37.30 + 30.20 + 21.00 + 4.75 + 6.0 + 0.071;            // 99.321
const double kHG_pi0D = 1.198 / kHG_pi0Sum, kPDG_pi0D = 0.01174;                // 0.011980
const double kHG_etaGG = 37.30 / kHG_etaSum, kPDG_etaGG = 0.3936;               // 0.37555
const double kHG_eta3pi0 = 30.20 / kHG_etaSum, kPDG_eta3pi0 = 0.3257;           // 0.30406
const double kHG_etaPPP0 = 21.00 / kHG_etaSum, kPDG_etaPPP0 = 0.2292;           // 0.21144
const double kHG_etaD = 6.0 / kHG_etaSum, kPDG_etaD = 0.0069;                   // 0.060410
const double kHG_etaPPG = 4.75 / kHG_etaSum, kPDG_etaPPG = 0.0422;              // 0.047825
// natural-decay correction factors PDG/HGeant (applied to weight-0 tracks in every set when Opt::pdg)
// (values in brackets: rules v1)
const double fPi0D   = kPDG_pi0D / kHG_pi0D;                                   // 0.9799 (0.9767) per distinct pi0 Dalitz decay
const double fEtaD   = kPDG_etaD / kHG_etaD;                                   // 0.1142 (0.1078) per distinct eta Dalitz decay
const double fEtaPPP0 = kPDG_etaPPP0 / kHG_etaPPP0;                            // 1.0840 (1.105)  eta -> pi+pi-pi0 identified (bump)
const double fEtaPPG = kPDG_etaPPG / kHG_etaPPG;                               // 0.8824 (0.849)  eta -> pi+pi-gamma identified (peak)
const double fEtaPi  = (kPDG_etaPPP0 + kPDG_etaPPG) / (kHG_etaPPP0 + kHG_etaPPG); // 1.0468 (1.056) pions of an eta, mode unknown
const double fEtaPi0 = (3 * kPDG_eta3pi0 + kPDG_etaPPP0) / (3 * kHG_eta3pi0 + kHG_etaPPP0); // 1.0736 (1.085) pi0 of an eta, mode unknown
const double fEtaGam = (2 * kPDG_etaGG + kPDG_etaD + kPDG_etaPPG) / (2 * kHG_etaGG + kHG_etaD + kHG_etaPPG); // 0.9732 (0.964) photon of an eta

// forced decays of the lepton set: probability of the forced channel per mother = BR / q, q = 1 / n_active
struct Forced { float code; double brDB, brPDG, q, wStored; };
const Forced kLepForced[] = {
  {75101.f,  1.198e-2, 1.174e-2, 1.0, 0.01198},    // pi0 -> e+e- gamma
  {175101.f, 0.6e-2,   0.69e-2,  1.0, 0.006},      // eta -> e+e- gamma
  {410203.f, 4.48e-5,  4.72e-5,  1.0, 4.5e-5},     // rho0 -> e+e-
  {520203.f, 7.15e-5,  7.38e-5,  0.5, 1.43e-4},    // omega -> e+e- (2 active channels)
  {525107.f, 5.9e-4,   7.7e-4,   0.5, 1.18e-3},    // omega -> pi0 e+e- (2 active channels); its e+- have g1 = 52051
  {550203.f, 2.99e-4,  2.979e-4, 1.0, 2.99e-4},    // phi -> e+e-
};
const float kEtaPPee = 1702030848.f, kEtaX17 = 17700808.f;  // eta set, as stored (rounded floats of 1702030809, 17700809)

// True masses: nominal masses [MeV/c^2], the values of FAT (src/pparticle.h)
const double kMassElectron = 0.51099895, kMassPion = 139.56995;

inline bool lepForcedCode(float g) { for (const Forced& f : kLepForced) if (g == f.code) return true; return false; }
inline bool forcedCode(Set s, float g) { return s == kLep ? lepForcedCode(g) : s == kEta ? (g == kEtaPPee || g == kEtaX17) : false; }
inline bool isE(const Trk& t) { return t.id == 2 || t.id == 3; }
// codes of SMASH hadronic decays of mothers that the lepton set forces into dileptons (absent there)
inline bool forcibleHadronicCode(float g) {
  return g == 410809.f || g == 520809.f || g == 52080908.f || g == 55080908.f || g == 520107.f || g == 551112.f || g == 551616.f || g == 550117.f;
}
inline bool forciblePi0Code(float g) { return g == 707.f || g == 420807.f || g == 430907.f; }

// FD: direct product of a forced-mode decay of a forcible mother (SMASH final-state pi0 / eta / rho0 / omega / phi).
inline bool FD(Set s, const Trk& t) {
  if (s == kLep) return forcedCode(s, t.g0) && t.proc == 0 && t.par == -1 && t.g1 != 7.f && t.g1 != 17.f;  // 7/17 = Dalitz photon
  if (s == kEta) return forcedCode(s, t.g0) && t.proc == 0 && t.par == -1;
  if (!isE(t) || t.proc != 5 || t.gp != -1) return false;
  if (t.par == 7) return forciblePi0Code(t.g0);
  if (t.par == 17) return t.g0 == 1717.f;
  return false;
}
// FC: any other descendant of a forced decay (conversion of its photon, delta / Compton / bremsstrahlung products).
inline bool FC(Set s, const Trk& t) { return s != kStd && forcedCode(s, t.g0) && !FD(s, t); }
// LFM (standard set only): a track linked to a forcible mother through a non-forced or unidentifiable decay mode.
inline bool LFM(const Trk& t) {
  if (FD(kStd, t)) return false;
  if (forcibleHadronicCode(t.g0) || t.g0 == 1717.f) return true;
  if (forciblePi0Code(t.g0) && isE(t) && (t.par == 1 || t.par == 2 || t.par == 3 || t.par == 7)) return true;
  return false;
}

// Ownership: exactly one set provides each row type. Rows of a set that it does not own are dropped.
// Mode A: std owns all rows except those with >= 1 FD and no LFM track; lep owns rows with >= 1 FD and no FC track.
// Mode B: std owns every row; lep owns rows with >= 1 FD track of rho0/omega/phi (410203, 520203, 525107, 550203) and no FC.
// Eta set (both modes): owns every row with >= 1 track of a forced eta decay (FD or FC); std has no eta -> pi+pi-e+e-.
inline bool Owns(Set s, const Trk t[4], const Opt& opt = Opt()) {
  int nFD = 0, nX = 0, nVM = 0;
  for (int i = 0; i < 4; ++i) {
    if (FD(s, t[i])) { ++nFD; if (s == kLep && t[i].g0 != 75101.f && t[i].g0 != 175101.f) ++nVM; }
    if (s == kStd ? LFM(t[i]) : FC(s, t[i])) ++nX;
  }
  if (s == kEta) return nFD + nX > 0;
  if (s == kStd) return opt.mode == 1 ? true : !(nFD > 0 && nX == 0);
  return (opt.mode == 1 ? nVM > 0 : nFD > 0) && nX == 0;
}

inline bool sameVtx(const Trk& a, const Trk& b) { return a.vx == b.vx && a.vy == b.vy && a.vz == b.vz; }
inline bool sameCodes(const Trk& a, const Trk& b) { return a.g0 == b.g0 && a.g1 == b.g1 && a.g2 == b.g2; }

// Natural-decay factors (weight-0 tracks; only identified modes are corrected). Returns 1 when !opt.pdg.
inline double NaturalFactor(const Trk t[4], const Opt& opt) {
  if (!opt.pdg) return 1.;
  double f = 1.;
  std::set<std::vector<float> > pi0D;                      // distinct pi0 Dalitz decays (creation vertex)
  std::map<std::pair<float,float>, int> eta;                // eta key (g1,g2) of g0 == 1717 tracks -> bit mask
  for (int i = 0; i < 4; ++i) {
    const Trk& x = t[i]; if (x.w > 0) continue;
    if (isE(x) && x.proc == 5 && x.par == 7) pi0D.insert({x.vx, x.vy, x.vz});
    if (x.g0 != 1717.f) continue;
    int& m = eta[std::make_pair(x.g1, x.g2)];
    if (isE(x) && x.proc == 5 && x.par == 17) m |= 1;                  // eta Dalitz
    else if (!isE(x) && x.proc == 5 && x.par == 17) m |= 2;            // pion (any hadron) of the eta
    else if (isE(x) && ((x.par == 7 && x.gp == 17) || (x.proc == 6 && x.gp == 7))) m |= 4;  // e+- from a pi0 of the eta (Dalitz or conversion)
    else if (isE(x) && x.proc == 6 && x.gp == 17) m |= 8;              // conversion of a photon of the eta
  }
  for (size_t k = 0; k < pi0D.size(); ++k) f *= fPi0D;
  for (std::map<std::pair<float,float>, int>::const_iterator it = eta.begin(); it != eta.end(); ++it) {
    int m = it->second;
    if (m & 1) { f *= fEtaD; continue; }                                // a Dalitz eta has no other products
    if ((m & 2) && (m & 4)) f *= fEtaPPP0;
    else if ((m & 2) && (m & 8)) f *= fEtaPPG;
    else if (m & 2) f *= fEtaPi;
    else { if (m & 4) f *= fEtaPi0; if (m & 8) f *= fEtaGam; }
  }
  return f;
}

// Physical weight W of a row of set s (defined for owned rows only). flags: bit0 stored genweight differs from the table,
// bit1 unknown forced code (cannot happen: FD of the lepton set already requires a code of the table; kept as in simrules.h),
// bit2 second-decay product together with a lepton of ANY omega -> pi0 e+e- decay in the row (its counter is not compared;
// the chain factor p(omega -> pi0 ee) is then not applied). (The same bits as in RowFlags.)
// NaturalFactor is 1 when !opt.pdg: the database variant keeps the HGeant branching ratios of the natural decays
// (e.g. a standard-set eta Dalitz pair has W_db = 1, not 0.108).
inline double RowWeight(Set s, const Trk t[4], const Opt& opt, int* flags = 0) {
  int fl = 0; double W = 1.;
  if (s == kLep) {
    std::map<std::pair<float,float>, double> dec; bool second = false, omegaLep = false;
    for (int i = 0; i < 4; ++i) {
      const Trk& x = t[i]; if (!FD(s, x)) continue;
      const Forced* F = 0; for (const Forced& f : kLepForced) if (x.g0 == f.code) F = &f;
      if (!F) { fl |= 2; continue; }
      bool sec = x.g0 == 525107.f && x.g1 == 7051.f;               // pi0 of omega -> pi0 e+e-, decayed again (Dalitz)
      double p;
      if (sec) { second = true; p = opt.pdg ? kPDG_pi0D : kLepForced[0].brDB; if (std::fabs(x.w - 1.44e-4f) > 2e-6f) fl |= 1; }
      else { p = (opt.pdg ? F->brPDG : F->brDB) / F->q; if (std::fabs(x.w - F->wStored) > 1e-3 * F->wStored + 1e-9) fl |= 1; }
      if (x.g0 == 525107.f && x.g1 == 52051.f) omegaLep = true;
      dec[std::make_pair(x.g0, x.g2)] = p;
    }
    for (std::map<std::pair<float,float>, double>::const_iterator it = dec.begin(); it != dec.end(); ++it) W *= it->second;
    if (second) { if (omegaLep) fl |= 4; else W *= (opt.pdg ? kLepForced[4].brPDG : kLepForced[4].brDB) / kLepForced[4].q; }
  } else if (s == kEta) {
    std::set<std::pair<float,float> > dec;
    for (int i = 0; i < 4; ++i) if (FD(s, t[i]) || FC(s, t[i])) dec.insert(std::make_pair(t[i].g0, t[i].g2));
    for (std::set<std::pair<float,float> >::const_iterator it = dec.begin(); it != dec.end(); ++it)
      W *= 2. * (it->first == kEtaPPee ? opt.brEtaPPee : opt.brX17);   // q = 1/2 (two forced channels)
  }
  W *= NaturalFactor(t, opt);
  if (flags) *flags = fl;
  return W;
}

// ---------------- classes ----------------
// lepton pair: 0 misidentified, 1 same pi0 Dalitz, 2 same eta Dalitz, 3 same conversion, 4 rho0/omega/phi -> ee,
// 5 omega -> pi0 ee, 6 second-decay pi0 Dalitz, 7 two different Dalitz decays (DD), 8 Dalitz + conversion (DC),
// 9 other different sources, 10 eta-set eta -> pi+pi-ee leptons, 11 X17 leptons
inline int LeptonPairClass(Set s, const Trk t[4], int q1, int q2) {   // q1, q2: expected GEANT ids of the slots
  const Trk& a = t[0]; const Trk& b = t[1];
  if (a.id != q1 || b.id != q2) return 0;
  bool fa = FD(s, a), fb = FD(s, b);
  if (s != kStd && fa && fb && a.g0 == b.g0 && a.g2 == b.g2 && a.g1 == b.g1) {
    if (a.g0 == 75101.f) return 1;
    if (a.g0 == 175101.f) return 2;
    if (a.g0 == 525107.f) return a.g1 == 52051.f ? 5 : 6;
    if (a.g0 == kEtaPPee) return 10;
    if (a.g0 == kEtaX17) return 11;
    return 4;
  }
  bool da = (a.proc == 5 && (a.par == 7 || a.par == 17)) || fa, db = (b.proc == 5 && (b.par == 7 || b.par == 17)) || fb;
  if (s == kStd || (!fa && !fb)) {
    if (a.proc == 5 && b.proc == 5 && a.par == b.par && (a.par == 7 || a.par == 17) && sameVtx(a, b) && a.gp == b.gp) return a.par == 7 ? 1 : 2;
  }
  if (a.proc == 6 && b.proc == 6 && sameVtx(a, b)) return 3;
  bool ca = a.proc == 6, cb = b.proc == 6;
  if (da && db) return 7;
  if ((da && cb) || (db && ca)) return 8;
  return 9;
}
// pion pair: 0 misidentified, 1 pion decayed to muon only, 2 same eta, 3 same omega, 4 same rho0, 5 same phi, 6 same K0S,
// 7 both primary (SMASH final state, different or no resonance), 8 mixed / secondary, 9 eta-set eta -> pi+pi-ee pions
inline bool sameRes(const Trk& a, const Trk& b) { return a.g0 == b.g0 && a.g1 == b.g1 && std::fabs(a.g2 - b.g2) < 1.5f && sameVtx(a, b); }
inline int PionPairClass(Set s, const Trk t[4]) {
  const Trk& A = t[2]; const Trk& B = t[3];
  bool m0 = A.id != 8, m1 = B.id != 9;
  if (m0 || m1) return ((!m0 || A.id == 5) && (!m1 || B.id == 6)) ? 1 : 0;
  if (s == kEta && FD(s, A) && FD(s, B) && A.g0 == B.g0 && A.g2 == B.g2) return 9;
  if (A.par == 17 && B.par == 17 && A.proc == 5 && B.proc == 5 && sameVtx(A, B) && sameCodes(A, B)) return 2;
  if (A.proc == 0 && B.proc == 0 && sameRes(A, B)) {
    if (A.g0 == 52080908.f || A.g0 == 520809.f) return 3;
    if (A.g0 == 410809.f) return 4;
    if (A.g0 == 55080908.f) return 5;
  }
  if (A.par == 16 && B.par == 16 && sameVtx(A, B)) return 6;
  if (A.proc == 0 && B.proc == 0) return 7;
  return 8;
}
// four-particle class: 0 fake (misidentified), 1 fake (pion decayed to muon), 2 bump via pi0 Dalitz, 3 bump via conversion,
// 4 peak via conversion (eta -> pi+pi-gamma), 5 peak eta -> pi+pi-e+e- (eta set), 6 X17 signal (eta set),
// 7 multi-pion with correlated e+e-, 8 CB-like (leptons from different sources)
inline int FourClass(Set s, const Trk t[4], int lp, int pp) {
  if (lp == 0 || pp == 0) return 0;
  if (pp == 1) return 1;
  const Trk& e = t[0]; const Trk& m = t[1]; const Trk& A = t[2];
  if (pp == 2 && sameCodes(e, A) && sameCodes(m, A)) {
    double d = std::sqrt(std::pow(e.vx - A.vx, 2) + std::pow(e.vy - A.vy, 2) + std::pow(e.vz - A.vz, 2));
    if (lp == 1 && e.gp == 17 && d < 0.01) return 2;
    if (lp == 3 && e.gp == 7) return 3;
    if (lp == 3 && e.gp == 17) return 4;
  }
  if (pp == 9 && (lp == 10 || lp == 11) && e.g2 == A.g2 && e.g0 == A.g0) return lp == 10 ? 5 : 6;
  if (lp >= 1 && lp <= 6) return 7;
  if (lp == 10 || lp == 11) return 7;
  return 8;
}

// ---------------- generator-code dictionary ----------------
// PAT stores the codes as Float_t: values above 2^24 are rounded (52080907 -> 52080908, 55080907 -> 55080908,
// 17700809 -> 17700808, 1702030809 -> 1702030848). CodeCandidates() maps a stored value of one code column
// back to the exact integer code(s) of the dictionary whose float value it is: every integer that rounds to the
// stored float is tested against every family (one integer below 2^24, up to 2^(e-23)+1 integers above).
//   column 0 (geninfo):  SMASH particle 100 g + g = 101 g (g = GEANT id 1..99) and the K0 / K0bar codes 1016, 1610;
//                        SMASH resonance-decay codes, PLUTO forced decays (lepton set), forced eta decays (eta set)
//                        (explicit list kDictCodes)
//   column 1 (geninfo1): 1000000 + SMASH process type (0..999), PLUTO / eta-set parent codes
//   column 2 (geninfo2): 1000000 + SMASH process number (any; a stored value >= 2^24 stands for several process
//                        numbers and is therefore ambiguous), PLUTO / eta-set decay counter (1..9999)
enum CodeColumn { kGeninfo = 0, kGeninfo1 = 1, kGeninfo2 = 2 };
enum CodeFamily { kFamSmashParticle = 1, kFamSmashDecay = 2, kFamPluto = 3, kFamEtaSet = 4,
                  kFamParent = 5, kFamSmashProcess = 6, kFamCounter = 7 };
struct Code { long long code; int family; };
struct DictCode { long long code; int column; int family; const char* text; };
const DictCode kDictCodes[] = {
  {52080907LL,   kGeninfo,  kFamSmashDecay, "omega -> pi+ pi- pi0"},
  {520809LL,     kGeninfo,  kFamSmashDecay, "omega -> pi+ pi-"},
  {410809LL,     kGeninfo,  kFamSmashDecay, "rho0 -> pi+ pi-"},
  {420807LL,     kGeninfo,  kFamSmashDecay, "rho+ -> pi+ pi0"},
  {430907LL,     kGeninfo,  kFamSmashDecay, "rho- -> pi- pi0"},
  {55080907LL,   kGeninfo,  kFamSmashDecay, "phi -> pi+ pi- pi0"},
  {520107LL,     kGeninfo,  kFamSmashDecay, "omega -> gamma pi0"},
  {551112LL,     kGeninfo,  kFamSmashDecay, "phi -> K+ K-"},
  {551010LL,     kGeninfo,  kFamSmashDecay, "phi -> K0 K0bar (GEANT K0L K0L)"},
  {551016LL,     kGeninfo,  kFamSmashDecay, "phi -> K0 K0bar (GEANT K0L K0S)"},
  {551610LL,     kGeninfo,  kFamSmashDecay, "phi -> K0 K0bar (GEANT K0S K0L)"},
  {551616LL,     kGeninfo,  kFamSmashDecay, "phi -> K0 K0bar (GEANT K0S K0S)"},
  {550117LL,     kGeninfo,  kFamSmashDecay, "phi -> gamma eta"},
  {75101LL,      kGeninfo,  kFamPluto,      "pi0 -> e+e- gamma (lepton set)"},
  {175101LL,     kGeninfo,  kFamPluto,      "eta -> e+e- gamma (lepton set)"},
  {410203LL,     kGeninfo,  kFamPluto,      "rho0 -> e+e- (lepton set)"},
  {520203LL,     kGeninfo,  kFamPluto,      "omega -> e+e- (lepton set)"},
  {525107LL,     kGeninfo,  kFamPluto,      "omega -> pi0 e+e- (lepton set)"},
  {550203LL,     kGeninfo,  kFamPluto,      "phi -> e+e- (lepton set)"},
  {1702030809LL, kGeninfo,  kFamEtaSet,     "eta -> pi+ pi- e+ e- (eta set)"},
  {17700809LL,   kGeninfo,  kFamEtaSet,     "eta -> X17 pi+ pi- (eta set)"},
  {7051LL,       kGeninfo1, kFamParent,     "e+- of a pi0 Dalitz decay (PLUTO)"},
  {17051LL,      kGeninfo1, kFamParent,     "e+- of an eta Dalitz decay (PLUTO)"},
  {52051LL,      kGeninfo1, kFamParent,     "e+- of omega -> pi0 e+e- (PLUTO)"},
  {41LL,         kGeninfo1, kFamParent,     "e+- of rho0 -> e+e- (PLUTO)"},
  {52LL,         kGeninfo1, kFamParent,     "e+- of omega -> e+e- (PLUTO)"},
  {55LL,         kGeninfo1, kFamParent,     "e+- of phi -> e+e- (PLUTO)"},
  {7LL,          kGeninfo1, kFamParent,     "photon of a pi0 Dalitz decay (PLUTO)"},
  {17LL,         kGeninfo1, kFamParent,     "photon of an eta Dalitz decay (PLUTO) / product of the eta (eta set)"},
  {70LL,         kGeninfo1, kFamParent,     "e+- of the X17 (eta set)"},
};

inline const char* CodeFamilyName(int family) {
  switch (family) {
    case kFamSmashParticle: return "SMASH particle (100 a + b)";
    case kFamSmashDecay:    return "SMASH resonance decay";
    case kFamPluto:         return "PLUTO forced decay (lepton set)";
    case kFamEtaSet:        return "forced eta decay (eta set)";
    case kFamParent:        return "PLUTO / eta-set parent code";
    case kFamSmashProcess:  return "1000000 + SMASH process";
    case kFamCounter:       return "PLUTO / eta-set decay counter";
  }
  return "unknown";
}

// Family of the exact integer code c of column col (0: not in the dictionary).
inline int CodeFamilyOf(int col, long long c) {
  for (const DictCode& d : kDictCodes) if (d.column == col && d.code == c) return d.family;
  if (col == kGeninfo && ((c >= 101 && c <= 9999 && c % 101 == 0) || c == 1016 || c == 1610)) return kFamSmashParticle;
  if (col == kGeninfo1 && c >= 1000000 && c <= 1000999) return kFamSmashProcess;
  if (col == kGeninfo2 && c >= 1000000) return kFamSmashProcess;
  if (col == kGeninfo2 && c >= 1 && c <= 9999) return kFamCounter;
  return 0;
}

// Candidates of the stored value v of code column col: the integers c with (float) c == v that belong to a
// family of the dictionary. Returns their number (0 = unknown value, > 1 = ambiguous); at most maxOut are
// written to out.
inline int CodeCandidates(int col, float v, Code out[], int maxOut) {
  if (!(v == std::floor(v))) return 0;                               // not an integer (or NaN): unknown
  const double dv = v;
  long long lo = static_cast<long long>(dv), hi = lo;
  if (std::fabs(dv) >= 16777216.) {                                  // float spacing > 1: several integers round to v
    const double up = static_cast<double>(std::nextafter(v, std::numeric_limits<float>::infinity())) - dv;
    const double dn = dv - static_cast<double>(std::nextafter(v, -std::numeric_limits<float>::infinity()));
    lo = static_cast<long long>(std::ceil(dv - dn / 2.));
    hi = static_cast<long long>(std::floor(dv + up / 2.));
  }
  int n = 0;
  for (long long c = lo; c <= hi; ++c) {
    if (static_cast<float>(c) != v) continue;                        // halfway cases: round-to-even decides
    const int fam = CodeFamilyOf(col, c);
    if (fam == 0) continue;
    if (n < maxOut) { out[n].code = c; out[n].family = fam; }
    ++n;
  }
  return n;
}

// ---------------- flags word per row ----------------
// Bits 0-2 come from RowWeight (as in simrules.h), bits 3-14 from RowFlags. None of them changes a class, an
// ownership or a weight of rules v1; they mark rows to inspect or to leave out.
const int kFlagWeightMismatch    = 1 << 0;  // stored genweight of a forced-decay product differs from the table
const int kFlagUnknownForced     = 1 << 1;  // reserved: never set (FD of the lepton set already requires a code of the table)
const int kFlagChainNotApplied   = 1 << 2;  // second-decay product with a lepton of ANY omega -> pi0 e+e- decay in the row
                                            // (counter not compared; chain factor not applied)
const int kFlagUnknownCode       = 1 << 3;  // a geninfo / geninfo1 / geninfo2 value of the four tracks is not in the dictionary
const int kFlagAmbiguousCode     = 1 << 4;  // a stored (rounded) code value has more than one dictionary candidate
const int kFlagElectronInPion    = 1 << 5;  // an e+ or e- (GEANT id 2, 3) sits in a pion slot
const int kFlagSecondDecay       = 1 << 6;  // a product of the second decay (pi0 of omega -> pi0 e+e-, 525107 / 7051) is in the row
const int kFlagSeveralDecays     = 1 << 7;  // more than one distinct forced decay in the row (W is a product)
const int kFlagWeightNotForced   = 1 << 8;  // a track has a stored weight > 0 but no forced code of the configured set
const int kFlagOtherSetCode      = 1 << 9;  // a track carries a forced code of another set (wrong simulation.set?)
const int kFlagGhostLepton       = 1 << 10; // a lepton slot holds a ghost track (<x>_sim_isghost == 1)
const int kFlagGhostPion         = 1 << 11; // a pion slot holds a ghost track
const int kFlagSameParticle      = 1 << 12; // both lepton slots carry the same GEANT particle (equal id, true momentum
                                            // and creation vertex): one particle counted twice (like-sign trees)
const int kFlagVmWithFC          = 1 << 13; // lepton set: a rho0 / omega / phi forced product together with an FC track.
                                            // Rules v1 give such a row to no set (lep drops it, std has no such decay)
const int kFlagV1CodeGap         = 1 << 14; // a code rules v1 do not handle: phi -> K0 K0bar variant 551010 / 551016 /
                                            // 551610 (missing from the forcible list: a std row with an FD lepton is
                                            // dropped in Mode A though no set provides it), or a product of the eta
                                            // of phi -> eta gamma (550117, no eta factor)

inline bool sameParticle(const Trk& a, const Trk& b) {
  return a.id > 0 && a.id == b.id && a.px == b.px && a.py == b.py && a.pz == b.pz && sameVtx(a, b);
}
inline bool vectorMesonForcedCode(float g) { return g == 410203.f || g == 520203.f || g == 525107.f || g == 550203.f; }

inline int RowFlags(Set s, const Trk t[4], const Opt& opt = Opt()) {
  int fl = 0;
  RowWeight(s, t, opt, &fl);                                       // bits 0-2
  std::set<std::pair<float,float> > dec;
  Code c[4];
  for (int i = 0; i < 4; ++i) {
    const Trk& x = t[i];
    const float g[3] = {x.g0, x.g1, x.g2};
    for (int k = 0; k < 3; ++k) {
      const int n = CodeCandidates(k, g[k], c, 4);
      if (n == 0) fl |= kFlagUnknownCode;
      if (n > 1) fl |= kFlagAmbiguousCode;
    }
    if (i >= 2 && isE(x)) fl |= kFlagElectronInPion;
    if (FD(s, x) && x.g0 == 525107.f && x.g1 == 7051.f) fl |= kFlagSecondDecay;
    if (s != kStd && (FD(s, x) || (s == kEta && FC(s, x)))) dec.insert(std::make_pair(x.g0, x.g2));
    if (x.w > 0 && !forcedCode(s, x.g0)) fl |= kFlagWeightNotForced;
    if ((s != kLep && lepForcedCode(x.g0)) || (s != kEta && (x.g0 == kEtaPPee || x.g0 == kEtaX17))) fl |= kFlagOtherSetCode;
  }
  if (dec.size() > 1) fl |= kFlagSeveralDecays;
  bool vm = false, fc = false;
  for (int i = 0; i < 4; ++i) {
    const Trk& x = t[i];
    if (x.ghost == 1) fl |= (i < 2 ? kFlagGhostLepton : kFlagGhostPion);
    if (s == kLep && FD(s, x) && vectorMesonForcedCode(x.g0)) vm = true;
    if (s == kLep && FC(s, x)) fc = true;
    if (x.g0 == 551010.f || x.g0 == 551016.f || x.g0 == 551610.f || (x.g0 == 550117.f && (x.par == 17 || x.gp == 17)))
      fl |= kFlagV1CodeGap;
  }
  if (vm && fc) fl |= kFlagVmWithFC;
  if (sameParticle(t[0], t[1])) fl |= kFlagSameParticle;
  return fl;
}

// ---------------- template class (sim_spec_final.md, section 1) ----------------
//   1 signal: four-particle classes 2-6 (bump Dalitz, bump conversion, peak conversion; eta set: peak
//     eta -> pi+pi-e+e-, X17), also with a ghost lepton (bit 10 tells their share)
//   2 physical background: a true e+e- pair from one source (L1-L6, eta set also L10/L11 outside the signal)
//     with any pion pair (misidentified pions and pi -> mu included), no fake lepton
//   3 fake lepton: misidentified lepton (L0), a ghost track in a lepton slot, or the same particle in both
//     lepton slots (this one has priority over the signal)
//   4 leptons from two different sources (L7-L9): no combinatorial background in the simulation analysis,
//     so these rows are neither signal nor background
// sim_cls 0 / 1 / 7 are NOT these templates: sim_cls 0 and 1 also hold rows with a true e+e- pair and a
// misidentified / decayed pion, and sim_cls 7 keeps rows with a ghost lepton.
inline int TemplateClass(int lcls, int cls, int flags) {
  if (flags & kFlagSameParticle) return 3;
  if (cls >= 2 && cls <= 6) return 1;
  if (lcls == 0 || (flags & kFlagGhostLepton)) return 3;
  if ((lcls >= 1 && lcls <= 6) || lcls == 10 || lcls == 11) return 2;
  return 4;
}

// ---------------- true masses ----------------
// Invariant mass [GeV/c^2] of the tracks i0..i1 (energies E [MeV], momenta P [MeV/c]).
inline double SumMass(const double E[4], const double P[4][3], int i0, int i1) {
  double e = 0., x = 0., y = 0., z = 0.;
  for (int i = i0; i <= i1; ++i) { e += E[i]; x += P[i][0]; y += P[i][1]; z += P[i][2]; }
  const double m2 = e * e - x * x - y * y - z * z;
  return (m2 > 0. ? std::sqrt(m2) : 0.) / 1000.;
}
// M(pi+pi-e+e-), M(ee), M(pipi) [GeV/c^2] from the true momenta [MeV/c] with the nominal masses;
// -1 when a track of the combination is misidentified (lepton slots: GEANT ids q1, q2; pion slots: 8, 9).
inline void TrueMasses(const Trk t[4], int q1, int q2, double& m4, double& mee, double& mpipi) {
  const int idOk[4] = {q1, q2, 8, 9};
  double E[4], P[4][3];
  bool ok[4];
  for (int i = 0; i < 4; ++i) {
    const double m = i < 2 ? kMassElectron : kMassPion;
    P[i][0] = t[i].px; P[i][1] = t[i].py; P[i][2] = t[i].pz;
    E[i] = std::sqrt(P[i][0] * P[i][0] + P[i][1] * P[i][1] + P[i][2] * P[i][2] + m * m);
    ok[i] = t[i].id == idOk[i];
  }
  mee   = (ok[0] && ok[1]) ? SumMass(E, P, 0, 1) : -1.;
  mpipi = (ok[2] && ok[3]) ? SumMass(E, P, 2, 3) : -1.;
  m4    = (ok[0] && ok[1] && ok[2] && ok[3]) ? SumMass(E, P, 0, 3) : -1.;
}

// ---------------- one call per row ----------------
struct RowTruth {
  int lcls, pcls, cls;   // lepton-pair, pion-pair and four-particle class
  double w, wdb;         // weight with opt (PDG BR by default) and with the converter-database BR
  bool ownA, ownB;       // the set owns the row in partition mode A / B
  int flags;             // RowFlags
  int tmpl;              // TemplateClass
  double m4, mee, mpipi; // true masses [GeV/c^2], -1 if misidentified
};

inline RowTruth EvaluateRow(Set s, const Trk t[4], int q1, int q2, const Opt& opt = Opt()) {
  RowTruth r;
  Opt oA = opt; oA.mode = 0;
  Opt oB = opt; oB.mode = 1;
  Opt oDB = opt; oDB.pdg = false;
  r.ownA = Owns(s, t, oA);
  r.ownB = Owns(s, t, oB);
  r.w = RowWeight(s, t, opt);
  r.wdb = RowWeight(s, t, oDB);
  r.lcls = LeptonPairClass(s, t, q1, q2);
  r.pcls = PionPairClass(s, t);
  r.cls = FourClass(s, t, r.lcls, r.pcls);
  r.flags = RowFlags(s, t, opt);
  r.tmpl = TemplateClass(r.lcls, r.cls, r.flags);
  TrueMasses(t, q1, q2, r.m4, r.mee, r.mpipi);
  return r;
}

}  // namespace simtruth

#endif // SIM_TRUTH_H
