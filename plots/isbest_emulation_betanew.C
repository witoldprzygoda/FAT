// isbest_emulation_betanew.C
// =========================================================================
// Data-only test of the "one best pair per event" (isBest) like-sign
// background estimator: emulation of the PAT best-hypothesis ranking on
// mixed (pure combinatorial) pairs, pattern by pattern.
//
// Reference: PAT3 master (git.gsi.de/przygoda/pat3, 73d8881, 2026-10-29),
// user/htimecut.cc (alpha = 0), user/hgraphcut.cc, core/hpidpool.cc:
//   * every hypothesis (pair) of a pattern gets the reduced TOF chi2 of the
//     relative-time reconstruction; for a two-track hypothesis
//         D    = (tof_1 - tof_2) - (tof_mom_1 - tof_mom_2)
//         chi2 = (D/2)^2 (1/res_1^2 + 1/res_2^2)          (ndf = 1)
//         tof_new_1 = tof_mom_1 + D/2,  tof_new_2 = tof_mom_2 - D/2
//         beta_new  = tracklength / (c tof_new)
//     with res = 0.250 ns (system 0) and 0.150 ns (system 1, 2), tof_mom the
//     flight time under the electron hypothesis.  Checked on lep_049.root:
//     chi2, tof_new and beta_new are reproduced exactly (relative
//     differences < 1e-3 only for chi2 < 1e-5, float rounding).
//   * PID: TCutG ep_cut_5sig / em_cut_5sig in the plane (p, beta_new) for
//     both legs (file pp45_pid_cuts_exp_betanew_v1.root; the copy on /hdd2
//     and the one in the PAT repository have identical polygons); a
//     hypothesis with a failing leg is inactive (isBest = -1, not written).
//   * isBest = 1 for the active hypothesis with the lowest chi2 of the
//     pattern in the event (checked: 0 violations on 2.7 M events).
//
// Method (one PAT tree = one pattern; patterns are ranked independently):
//   1. Events are rebuilt from adjacent rows with equal (runnumber, event);
//      legs are unique by float equality of (p, theta, phi).  The native
//      hypotheses are exactly the stored rows (stored chi2, isBest).
//   2. Time frame of the event: only differences of tof - tof_mom between
//      tracks of the same event carry information (the START time is not
//      usable between events: rms of tof - tof_mom is 1.5 ns, while the
//      in-event difference D has rms 0.24 ns).  Per event the residuals
//      delta_l = tof_l - tof_mom_l are expressed relative to the
//      inverse-variance weighted mean (weights 1/sigma_c^2 with sigma_c
//      calibrated in pass 0 from the D distribution per leg category
//      (system x momentum bin): Var(D_cc') = sigma_c^2 + sigma_c'^2).
//      Relative to that frame the unknown common mode has variance 1/W,
//      W = sum 1/sigma_c^2; for Gaussian legs the joint statistics of the
//      residuals are then reproduced exactly when the common mode is
//      sampled from N(0, 1/W).
//   3. Leg replacement: every leg a_k of event A is replaced in turn by a
//      leg b of the same charge from a pooled earlier event B (same trigger,
//      rolling pool of pool_depth events; R donors per leg).  The cross
//      pairs (a_i, b), i != k, are pure combinatorial pairs in the real
//      event environment with unchanged multiplicity; their D is
//      ehat_i - ehat_b + c with c ~ N(0, 1/W_A + 1/W_B).  For every cross
//      pair chi2, tof_new, beta_new and the PID are recomputed as in PAT;
//      the ranking is done among the native hypotheses not involving a_k
//      (stored chi2) and the active cross pairs.
//      Per trigger (PT3 = trigbit 8192, PT2 = 4096), mass bin (58-bin recipe
//      grid + 280 x 5 MeV) and multiplicity class c = total number of
//      identified leptons in the event n = n+ + n- (c2: n = 2, c3, c4,
//      c5: n >= 5; for the like-sign trees the other-charge count comes from
//      a lookup of the same (runnumber, event) in EpEm_ID of the file):
//        h_all      native pairs             h_best      native isBest == 1
//        h_cross    active cross pairs       h_crosswin  cross pairs that win
//        h_nopid    cross pairs failing PID
//      eps_data = best/all, eps_mix = crosswin/cross, per class.  For the
//      like-sign patterns both are measurable -> closure test of the
//      emulation; for e+e- eps_mix is the win probability of a combinatorial
//      pair.  The classes matter because the cross pairs of an event are
//      generated whether or not its native pairs are combinatorial: a
//      two-lepton event with a true pair yields cross pairs but no real
//      combinatorial pair, so the class mix of the cross pairs is not the
//      class mix of the real background.  The analysis macro therefore
//      combines per class: CB_best(m) = sum_c eps_B,c(m) CB_c(m) with
//      CB_c = 2 sqrt(N++_c N--_c) (exact inside a class of fixed total
//      multiplicity for independent legs with any charge asymmetry; for a
//      correlated low-mass pair plus one extra leg, which dominates n = 3,
//      the arithmetic mean N++ + N-- is exact instead: 1.7% difference,
//      quoted as a systematic by the analysis macro), and compares with the
//      naive 2 sqrt(N++_best N--_best).
//   4. Diagnostics: D and log10(chi2) of native and cross pairs in
//      [0.14, 0.80] GeV, counters (h_counters), calibrated sigma table.
//      Assumptions of the leg replacement (documented after review 2026-10-31):
//      the TOF chi2 is blind to the origin of a pair (a combinatorial pair
//      beats a true pair with probability 1/2), legs are exchangeable within
//      a class, and the per-leg residuals are Gaussian (a misidentified
//      hadron pulls the event frame; cross pairs then have a slightly
//      depleted chi2 tail).  Pairs with opening angle < 9 deg have narrower
//      D (0.18 vs 0.24 ns, shared detector elements), which no mixing
//      reproduces: the results are meant for M > 0.14 GeV.
//   5. Event selection: vertex_z_min applies the FAT event cut
//      eVertReco_z > vertex_z_min (FAT main.cc / setup_cuts.h: -500) to the
//      events (as analysed and as donors), so that h_all / h_best equal the
//      FAT dilepton_nt spectra; vertex_z_min = -1e9 disables it.  NOTE: on
//      the betanew files eVertReco_z = -1000 (no track-candidate vertex) in
//      61% of the e+e- rows (day 049), with a strong opening-angle / mass
//      dependence and PT3 39% vs PT2 49% kept.
//
// Usage (ACLiC), one PAT file and one pattern per call:
//   root -l -b -q 'plots/isbest_emulation_betanew.C+("/hdd3/.../lep_049.root","pm","plots/output/isbest_emu/pm_049.root")'
//   patterns: "pm" (EpEm_ID), "pp" (EpEp_ID), "mm" (EmEm_ID)
//   then hadd per pattern over the days (all histograms are sums).
//
// @author Witold Przygoda (witold.przygoda@uj.edu.pl)
// @date 2026

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCutG.h>
#include <TRandom3.h>
#include <TNamed.h>
#include <TStopwatch.h>
#include <TMath.h>
#include <TString.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <deque>
#include <array>
#include <unordered_map>
#include <cstring>
#include <string>
#include <vector>

namespace ibe {

constexpr double kMe     = 0.51099895;   // MeV
constexpr double kC      = 299.792458;   // mm/ns
constexpr int    kNb58   = 58;
constexpr int    kNbFine = 280;
constexpr int    kNTrig  = 2;            // 0 = PT3 (8192), 1 = PT2 (4096)
constexpr int    kNCl    = 4;            // multiplicity classes: n = 2, 3, 4, >= 5
constexpr int    kNCat   = 12;           // system (0,1,2) x momentum bin (4)
constexpr double kMidLo  = 0.14, kMidHi = 0.80;

const char* trigName(int t) { return t == 0 ? "PT3" : "PT2"; }
const char* clName(int c) { static const char* n[kNCl] = {"c2", "c3", "c4", "c5"}; return n[c]; }
int clIndex(int ntot) { return ntot <= 2 ? 0 : (ntot >= 5 ? 3 : ntot - 2); }

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

int momBin(double p) { return p < 150. ? 0 : p < 300. ? 1 : p < 500. ? 2 : 3; }
int category(int sys, double p) {
    int s = sys < 0 ? 0 : (sys > 2 ? 2 : sys);
    return s * 4 + momBin(p);
}

struct Leg {
    int    q = 0;                 // charge
    float  p = 0, th = 0, ph = 0; // identity key (float equality); p = PID momentum
    double pc = 0;                // momentum used for the pair mass = p (FAT m_ee, RECONSTRUCTED kinematics from ep_p/em_p)
    double tof = 0, tofmom = 0, L = 0, res = 0, delta = 0;
    int    sys = 0, cat = 0;
    double sig2 = 0;              // calibrated residual variance (ns^2)
    double ehat = 0;              // residual relative to the event frame
};
struct Hyp { int i = 0, j = 0; double chi2 = 0; int best = 0; double mass = 0; };
struct Event {
    int trig = -1;
    Float_t run = 0, ev = 0;
    Float_t vx = 0, vy = 0, vz = 0;   // reconstructed event vertex (lookup key in sim, see vertex_key)
    bool vtxOk = true;
    int nOther = 0;               // identified leptons of the charge absent from this tree (like-sign trees, from the EpEm_ID lookup)
    int nSameLookup = -1;         // leptons of this tree's charge in the EpEm_ID view (-1: event absent from EpEm_ID / pm tree)
    double W = 0;
    std::vector<Leg> legs;
    std::vector<Hyp> hyps;
};

struct Spec {
    TH1D* h58 = nullptr;
    TH1D* hf  = nullptr;
    void init(const std::string& name, const std::string& title) {
        h58 = new TH1D(name.c_str(), (title + ";M_{ee} [GeV/c^{2}];pairs").c_str(), kNb58, getEdges58());
        hf  = new TH1D((name + "_fine").c_str(), (title + " (5 MeV);M_{ee} [GeV/c^{2}];pairs").c_str(), kNbFine, 0., 1.4);
        h58->Sumw2(); hf->Sumw2();
    }
    void fill(double m) { h58->Fill(m); hf->Fill(m); }
};

double pairMass(const Leg& a, const Leg& b) {
    const double ta = a.th * TMath::DegToRad(), pa = a.ph * TMath::DegToRad();
    const double tb = b.th * TMath::DegToRad(), pb = b.ph * TMath::DegToRad();
    const double ax = a.pc * std::sin(ta) * std::cos(pa), ay = a.pc * std::sin(ta) * std::sin(pa), az = a.pc * std::cos(ta);
    const double bx = b.pc * std::sin(tb) * std::cos(pb), by = b.pc * std::sin(tb) * std::sin(pb), bz = b.pc * std::cos(tb);
    const double Ea = std::sqrt(a.pc * a.pc + kMe * kMe), Eb = std::sqrt(b.pc * b.pc + kMe * kMe);
    const double m2 = 2. * kMe * kMe + 2. * (Ea * Eb - (ax * bx + ay * by + az * bz));
    return m2 > 0. ? std::sqrt(m2) / 1000. : 0.;
}

// PAT HTimeCut, two-track hypothesis (alpha = 0): chi2 and the new times.
inline double patChi2(double D, double res1, double res2) {
    return (D / 2.) * (D / 2.) * (1. / (res1 * res1) + 1. / (res2 * res2));
}

enum Counter { kRows = 1, kEvents, kEvPT3, kEvPT2, kEvOther, kMulti, kChi2Mismatch, kBetaMismatch,
               kBestNotMin, kBestCountBad, kCrossActive, kCrossWin, kNoPid, kNoDonor, kLegDupInconsistent, kEvVertexFail, kLookupMiss, kLookupMismatch, kLegsBeyondPm, kHypRepeated, kNCounters };

} // namespace ibe

void isbest_emulation_betanew(const char* pat_file,
                              const char* pattern   = "pm",
                              const char* out_file  = "plots/output/isbest_emu/pm_test.root",
                              int   R_pt3           = 2,
                              int   R_pt2           = 10,
                              int   pool_depth      = 20,
                              Long64_t max_rows     = -1,
                              int   seed            = 1,
                              const char* cut_file  = "/hdd2/przygoda/hades/pp45/exp/pp45_pid_cuts_exp_betanew_v1.root",
                              Long64_t calib_rows   = 3000000,
                              double vertex_z_min   = -500.,
                              bool   sim_trigger_bits = false,
                              bool   vertex_key       = false)
{
    using namespace ibe;
    TStopwatch sw; sw.Start();
    const std::string pat = pattern;
    std::string tree, pre1, pre2;
    int q1 = 0, q2 = 0;
    if (pat == "pm")      { tree = "EpEm_ID"; pre1 = "ep_";  pre2 = "em_";  q1 = +1; q2 = -1; }
    else if (pat == "pp") { tree = "EpEp_ID"; pre1 = "ep1_"; pre2 = "ep2_"; q1 = +1; q2 = +1; }
    else if (pat == "mm") { tree = "EmEm_ID"; pre1 = "em1_"; pre2 = "em2_"; q1 = -1; q2 = -1; }
    else { printf("ERROR: pattern must be pm, pp or mm\n"); return; }

    // ---- PID cuts (as PAT HGraphCut) ----------------------------------------
    TFile* fcut = TFile::Open(cut_file);
    if (!fcut || fcut->IsZombie()) { printf("ERROR: cannot open cut file %s\n", cut_file); return; }
    TCutG* cut_ep = dynamic_cast<TCutG*>(fcut->Get("ep_cut_5sig"));
    TCutG* cut_em = dynamic_cast<TCutG*>(fcut->Get("em_cut_5sig"));
    if (!cut_ep || !cut_em) { printf("ERROR: ep_cut_5sig / em_cut_5sig missing in %s\n", cut_file); return; }
    auto pidOk = [&](int q, double p, double beta_new) {
        return q > 0 ? cut_ep->IsInside(p, beta_new) : cut_em->IsInside(p, beta_new);
    };

    // ---- input --------------------------------------------------------------
    TFile* fin = TFile::Open(pat_file);
    if (!fin || fin->IsZombie()) { printf("ERROR: cannot open %s\n", pat_file); return; }
    TTree* t = dynamic_cast<TTree*>(fin->Get(tree.c_str()));
    if (!t) { printf("ERROR: tree %s missing in %s\n", tree.c_str(), pat_file); return; }
    t->SetBranchStatus("*", 0);
    Float_t v_run, v_ev, v_trig, v_best, v_chi2, v_vz, v_vx = 0, v_vy = 0;
    struct LegBr { Float_t p, th, ph, tof, tofmom, L, res, sys, bnew, tnew; } b1, b2;
    auto setb = [&](const std::string& n, Float_t* a) {
        if (!t->GetBranch(n.c_str())) { printf("ERROR: branch %s missing\n", n.c_str()); throw 1; }
        t->SetBranchStatus(n.c_str(), 1); t->SetBranchAddress(n.c_str(), a);
    };
    auto setLeg = [&](const std::string& pre, LegBr& b) {
        setb(pre + "p", &b.p); setb(pre + "theta", &b.th); setb(pre + "phi", &b.ph);
        setb(pre + "tof", &b.tof); setb(pre + "tof_mom", &b.tofmom); setb(pre + "tracklength", &b.L);
        setb(pre + "resolution", &b.res); setb(pre + "system", &b.sys); setb(pre + "beta_new", &b.bnew); setb(pre + "tof_new", &b.tnew);
    };
    try {
        setb("runnumber", &v_run); setb("event", &v_ev); setb("trigbit", &v_trig); setb("isBest", &v_best); setb("eVertReco_z", &v_vz);
        if (vertex_key) { setb("eVertReco_x", &v_vx); setb("eVertReco_y", &v_vy); }
        setb(pre1 + "hyp_best_chi2_red", &v_chi2);
        setLeg(pre1, b1); setLeg(pre2, b2);
    } catch (...) { return; }
    const Long64_t nEntries = t->GetEntries();
    const Long64_t nRows = (max_rows >= 0 && max_rows < nEntries) ? max_rows : nEntries;
    printf("[isbest_emu] %s  tree %s  rows %lld (of %lld)  pattern %s  R(PT3)=%d R(PT2)=%d pool=%d seed=%d\n",
           pat_file, tree.c_str(), nRows, nEntries, pattern, R_pt3, R_pt2, pool_depth, seed);

    // ---- pass 0: residual widths per leg category ---------------------------
    // Always from EpEm_ID of the same file (largest statistics), pairs with
    // opening angle > 9 deg only (small-angle pairs have correlated residuals:
    // their D is narrower, e.g. 0.095 vs 0.128 ns for two system-0 legs
    // above 500 MeV/c).  sigma_c = sigma68(D_cc)/sqrt(2) from the diagonal
    // combination; categories without a diagonal sample use the most
    // populated partner with a known sigma: sigma_c^2 = V_cd - sigma_d^2.
    std::vector<double> sig2(kNCat, 0.04);
    {
        // separate file handle: for the pm pattern the analysed tree is EpEm_ID
        // itself and its branch addresses must stay untouched
        TFile* fcal = TFile::Open(pat_file);
        TTree* tc = fcal ? dynamic_cast<TTree*>(fcal->Get("EpEm_ID")) : nullptr;
        if (!tc) { printf("ERROR: EpEm_ID (calibration) missing in %s\n", pat_file); return; }
        tc->SetBranchStatus("*", 0);
        Float_t c_tof1, c_tof2, c_tm1, c_tm2, c_sys1, c_sys2, c_p1, c_p2, c_th1, c_th2, c_ph1, c_ph2;
        auto setc = [&](const char* n, Float_t* a) { tc->SetBranchStatus(n, 1); tc->SetBranchAddress(n, a); };
        setc("ep_tof", &c_tof1); setc("em_tof", &c_tof2); setc("ep_tof_mom", &c_tm1); setc("em_tof_mom", &c_tm2);
        setc("ep_system", &c_sys1); setc("em_system", &c_sys2); setc("ep_p", &c_p1); setc("em_p", &c_p2);
        setc("ep_theta", &c_th1); setc("em_theta", &c_th2); setc("ep_phi", &c_ph1); setc("em_phi", &c_ph2);
        std::vector<TH1D*> hD(kNCat * kNCat, nullptr);
        for (int a = 0; a < kNCat; ++a)
            for (int b = a; b < kNCat; ++b)
                hD[a * kNCat + b] = new TH1D(Form("hDcal_%d_%d", a, b), "", 4000, -5., 5.);
        const Long64_t nCal = std::min(tc->GetEntries(), calib_rows);
        const double d2r = TMath::DegToRad();
        for (Long64_t i = 0; i < nCal; ++i) {
            tc->GetEntry(i);
            const double cosoa = std::sin(c_th1 * d2r) * std::sin(c_th2 * d2r) * std::cos((c_ph1 - c_ph2) * d2r) + std::cos(c_th1 * d2r) * std::cos(c_th2 * d2r);
            if (cosoa > std::cos(9. * d2r)) continue;
            int a = category(int(c_sys1), c_p1), b = category(int(c_sys2), c_p2);
            if (a > b) std::swap(a, b);
            hD[a * kNCat + b]->Fill((c_tof1 - c_tof2) - (c_tm1 - c_tm2));
        }
        std::vector<double> V(kNCat * kNCat, -1.), Nv(kNCat * kNCat, 0.);
        for (int a = 0; a < kNCat; ++a)
            for (int b = a; b < kNCat; ++b) {
                TH1D* h = hD[a * kNCat + b];
                if (h->GetEntries() < 2000) continue;
                double q[3], pr[3] = {0.1587, 0.5, 0.8413};
                h->GetQuantiles(3, q, pr);
                const double s68 = (q[2] - q[0]) / 2.;
                V[a * kNCat + b] = s68 * s68;
                Nv[a * kNCat + b] = h->GetEntries();
            }
        std::vector<int> how(kNCat, 0);   // 0 default, 1 diagonal, 2 partner
        for (int c = 0; c < kNCat; ++c)
            if (V[c * kNCat + c] >= 0.) { sig2[c] = V[c * kNCat + c] / 2.; how[c] = 1; }
        for (int c = 0; c < kNCat; ++c) {
            if (how[c] == 1) continue;
            int bestd = -1; double bestn = 0.;
            for (int d = 0; d < kNCat; ++d) {
                if (d == c || how[d] != 1) continue;
                const int idx = c < d ? c * kNCat + d : d * kNCat + c;
                if (V[idx] >= 0. && Nv[idx] > bestn) { bestn = Nv[idx]; bestd = d; }
            }
            if (bestd >= 0) {
                const int idx = c < bestd ? c * kNCat + bestd : bestd * kNCat + c;
                sig2[c] = std::max(1e-4, V[idx] - sig2[bestd]); how[c] = 2;
            }
        }
        printf("[isbest_emu] residual sigma per category (system 0,1,2 x p-bin <150,<300,<500,>=500 MeV/c), from EpEm_ID pairs with OA > 9 deg (%lld rows):\n", nCal);
        for (int c = 0; c < kNCat; ++c)
            printf("   cat %2d (sys %d, pbin %d): sigma = %.3f ns  [%s]\n", c, c / 4, c % 4, std::sqrt(sig2[c]),
                   how[c] == 1 ? "diagonal" : (how[c] == 2 ? "partner" : "default"));
        double sumdev = 0.; int ndev = 0;
        for (int a = 0; a < kNCat; ++a)
            for (int b = a + 1; b < kNCat; ++b) {
                if (V[a * kNCat + b] < 0. || how[a] != 1 || how[b] != 1) continue;
                const double pred = std::sqrt(sig2[a] + sig2[b]), meas = std::sqrt(V[a * kNCat + b]);
                sumdev += (meas - pred) / pred; ++ndev;
            }
        if (ndev > 0) printf("[isbest_emu] off-diagonal check: mean (measured - predicted)/predicted sigma68(D) = %+.3f over %d combinations\n", sumdev / ndev, ndev);
        for (auto* h : hD) if (h) delete h;
        fcal->Close();
    }
    // ---- multiplicity lookup: (runnumber, event) -> (n+, n-) from EpEm_ID ---
    // Every identified e+ of an event with n+ >= 1 and n- >= 1 is paired with
    // every e- (all combinations written), so the unique legs of the EpEm_ID
    // rows give n+ and n-.  Events absent from EpEm_ID have n+ = 0 or n- = 0.
    // (runnumber, event) repeats between the sub-file blocks of a run for
    // 0.5% of the events (see mixing_build_betanew.C); such a collision gives
    // a wrong class for a like-sign event, which is neglected.
    // Key: (runnumber, event); with vertex_key also the reconstructed event
    // vertex (sim files: runnumber = 0 and the event counter restarts every
    // few dozen events, so (runnumber, event) is far from unique there).
    auto makeKey = [&](Float_t run, Float_t evn, Float_t vx, Float_t vy, Float_t vz) -> unsigned long long {
        unsigned long long k = (unsigned long long)(unsigned int)(long long)run << 32 | (unsigned int)(long long)evn;
        if (vertex_key) {
            auto mix = [](unsigned long long h, Float_t x) { unsigned int u; std::memcpy(&u, &x, 4); h ^= u + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2); return h; };
            k = mix(mix(mix(k * 0x100000001b3ULL, vx), vy), vz);
        }
        return k;
    };
    std::unordered_map<unsigned long long, std::pair<unsigned char, unsigned char> > lookup;
    {
        TFile* flk = TFile::Open(pat_file);
        TTree* tl = flk ? dynamic_cast<TTree*>(flk->Get("EpEm_ID")) : nullptr;
        if (!tl) { printf("ERROR: EpEm_ID (lookup) missing in %s\n", pat_file); return; }
        tl->SetBranchStatus("*", 0);
        Float_t l_run, l_ev, l_p1, l_th1, l_ph1, l_p2, l_th2, l_ph2, l_best = 1, l_vx = 0, l_vy = 0, l_vz = 0;
        auto setl = [&](const char* n, Float_t* a) { tl->SetBranchStatus(n, 1); tl->SetBranchAddress(n, a); };
        setl("runnumber", &l_run); setl("event", &l_ev); setl("ep_p", &l_p1); setl("ep_theta", &l_th1); setl("ep_phi", &l_ph1);
        setl("em_p", &l_p2); setl("em_theta", &l_th2); setl("em_phi", &l_ph2); setl("isBest", &l_best);
        if (vertex_key) { setl("eVertReco_x", &l_vx); setl("eVertReco_y", &l_vy); setl("eVertReco_z", &l_vz); }
        lookup.reserve(size_t(tl->GetEntries() / 1.05) + 16);
        std::vector<std::array<float, 3> > lp, lm;
        Float_t cr = -1.f, ce = -1.f, cvx = 0, cvy = 0, cvz = 0; bool hv = false;
        auto flush = [&]() {
            if (!hv) return;
            const unsigned long long key = makeKey(cr, ce, cvx, cvy, cvz);
            lookup[key] = std::make_pair((unsigned char)std::min<size_t>(lp.size(), 255), (unsigned char)std::min<size_t>(lm.size(), 255));
            lp.clear(); lm.clear();
        };
        auto addu = [](std::vector<std::array<float, 3> >& v, float p, float th, float ph) {
            for (const auto& x : v) if (x[0] == p && x[1] == th && x[2] == ph) return;
            v.push_back({p, th, ph});
        };
        const Long64_t nl = tl->GetEntries();
        for (Long64_t i = 0; i < nl; ++i) {
            tl->GetEntry(i);
            if (l_best < 0.f) continue;   // inactive hypotheses (sim productions only)
            if (!hv || l_run != cr || l_ev != ce) { flush(); cr = l_run; ce = l_ev; cvx = l_vx; cvy = l_vy; cvz = l_vz; hv = true; }
            addu(lp, l_p1, l_th1, l_ph1); addu(lm, l_p2, l_th2, l_ph2);
        }
        flush();
        flk->Close();
        printf("[isbest_emu] lookup: %zu e+e- events (%lld rows, %.0f s)\n", lookup.size(), nl, sw.RealTime()); sw.Continue();
    }
    auto lookupOther = [&](Float_t run, Float_t evn, Float_t vx, Float_t vy, Float_t vz, int q) -> int {   // leptons of charge q in the e+e- event, -1 if absent
        const unsigned long long key = makeKey(run, evn, vx, vy, vz);
        auto it = lookup.find(key);
        if (it == lookup.end()) return -1;
        return q > 0 ? int(it->second.first) : int(it->second.second);
    };
    // Class of a like-sign event: the e+e- pattern defines its class from the
    // legs active in e+e- hypotheses (n+ + n- of EpEm_ID), and eps_B,c is
    // measured on that population.  The PID is hypothesis dependent (beta_new
    // depends on the partner), so a leg can be active in a like-sign hypothesis
    // but in no e+e- hypothesis; for events present in EpEm_ID the class is
    // therefore taken from the EpEm_ID view (n+ + n- of the lookup), which is
    // the definition consistent with eps_B,c, and the like-sign tree's own leg
    // count is used only for like-sign-only events (absent from EpEm_ID).  The
    // counter "like-sign legs beyond EpEm_ID view" counts the events where the
    // own count of the tree's charge exceeds the EpEm_ID count (class-definition
    // systematic, see the analysis macro).

    TH1D* h_sigma = new TH1D(Form("h_sigma_cat_%s", pattern), "calibrated residual sigma per leg category (system*4 + momentum bin);category;#sigma [ns]", kNCat, -0.5, kNCat - 0.5);
    for (int c = 0; c < kNCat; ++c) h_sigma->SetBinContent(c + 1, std::sqrt(sig2[c]));
    h_sigma->SetBit(TH1::kIsAverage);   // hadd averages instead of summing the per-file tables

    // ---- histograms ---------------------------------------------------------
    const char* kinds[5] = {"all", "best", "cross", "crosswin", "nopid"};
    const char* kindTitle[5] = {"native pairs", "native isBest==1", "active cross pairs", "winning cross pairs", "cross pairs failing PID"};
    Spec spec[5][kNTrig][kNCl];
    for (int k = 0; k < 5; ++k)
        for (int tr = 0; tr < kNTrig; ++tr)
            for (int c = 0; c < kNCl; ++c)
                spec[k][tr][c].init(Form("h_%s_%s_%s_%s", kinds[k], pattern, trigName(tr), clName(c)),
                                    Form("%s %s %s %s", kindTitle[k], pattern, trigName(tr), clName(c)));
    TH1D* h_D_nat[kNTrig];  TH1D* h_D_cross[kNTrig];  TH1D* h_chi2_nat[kNTrig];  TH1D* h_chi2_cross[kNTrig];
    TH1D* h_chi2_best[kNTrig];
    for (int tr = 0; tr < kNTrig; ++tr) {
        h_D_nat[tr]      = new TH1D(Form("h_D_native_%s_%s", pattern, trigName(tr)),     "D native, 0.14<M<0.80;D [ns]",            400, -2., 2.);
        h_D_cross[tr]    = new TH1D(Form("h_D_cross_%s_%s", pattern, trigName(tr)),      "D cross, 0.14<M<0.80;D [ns]",             400, -2., 2.);
        h_chi2_nat[tr]   = new TH1D(Form("h_chi2_native_%s_%s", pattern, trigName(tr)),  "log10 chi2 native, 0.14<M<0.80;log_{10}#chi^{2}", 200, -8., 6.);
        h_chi2_best[tr]  = new TH1D(Form("h_chi2_best_%s_%s", pattern, trigName(tr)),    "log10 chi2 native isBest==1, 0.14<M<0.80;log_{10}#chi^{2}", 200, -8., 6.);
        h_chi2_cross[tr] = new TH1D(Form("h_chi2_cross_%s_%s", pattern, trigName(tr)),   "log10 chi2 cross (active), 0.14<M<0.80;log_{10}#chi^{2}", 200, -8., 6.);
        for (TH1D* h : {h_D_nat[tr], h_D_cross[tr], h_chi2_nat[tr], h_chi2_best[tr], h_chi2_cross[tr]}) h->Sumw2();
    }
    TH1D* h_nlegs[kNTrig]; TH1D* h_ntot[kNTrig];
    for (int tr = 0; tr < kNTrig; ++tr) {
        h_nlegs[tr] = new TH1D(Form("h_nlegs_%s_%s", pattern, trigName(tr)), "legs of this tree per event;n legs", 12, -0.5, 11.5);
        h_ntot[tr]  = new TH1D(Form("h_ntot_%s_%s", pattern, trigName(tr)), "identified leptons per event (n+ + n-);n", 12, -0.5, 11.5);
    }
    TH1D* h_cnt = new TH1D(Form("h_counters_%s", pattern), "counters", kNCounters, 0.5, kNCounters + 0.5);
    const char* cntName[kNCounters] = {"", "rows", "events", "events PT3", "events PT2", "events other trigger", "events >=2 hypotheses",
        "native chi2 mismatch", "native beta_new mismatch", "isBest not min chi2", "#isBest!=1 per event", "cross pairs active",
        "cross pairs winning", "cross pairs failing PID", "no donor", "duplicate leg inconsistent", "events failing vertex cut", "lookup miss (like-sign event absent from EpEm_ID)", "lookup mismatch (pm self-check)", "like-sign legs beyond EpEm_ID view", "repeated hypotheses dropped"};
    for (int c = 1; c < kNCounters; ++c) h_cnt->GetXaxis()->SetBinLabel(c, cntName[c]);
    auto cnt = [&](int c, double w = 1.) { h_cnt->AddBinContent(c, w); };

    TRandom3 rnd(seed);
    std::deque<Event> pool[kNTrig];

    // ---- event processing ---------------------------------------------------
    auto processEvent = [&](Event& ev) {
        const int n = int(ev.legs.size());
        if (n < 2 || ev.trig < 0) return;
        if (!ev.vtxOk) { cnt(kEvVertexFail); return; }   // same event selection as FAT main.cc (vertex_z > -500 on eVertReco_z)
        const int tr = ev.trig;
        int ntot = n + std::max(0, ev.nOther);
        if (pat != "pm" && ev.nSameLookup >= 0) {           // event present in EpEm_ID: class from the EpEm_ID view
            if (n > ev.nSameLookup) cnt(kLegsBeyondPm);
            ntot = ev.nSameLookup + std::max(0, ev.nOther);
        }
        const int tp = clIndex(ntot);
        h_nlegs[tr]->Fill(n); h_ntot[tr]->Fill(ntot);
        if (ev.hyps.size() > 1) cnt(kMulti);
        // frame
        double W = 0., num = 0.;
        for (auto& l : ev.legs) { W += 1. / l.sig2; num += l.delta / l.sig2; }
        ev.W = W;
        const double frame = num / W;
        for (auto& l : ev.legs) l.ehat = l.delta - frame;
        // native statistics and checks
        if (pat == "pm") {   // self-check of the lookup on the analysed events
            int np = 0, nm = 0; for (const auto& l : ev.legs) (l.q > 0 ? np : nm)++;
            if (lookupOther(ev.run, ev.ev, ev.vx, ev.vy, ev.vz, +1) != np || lookupOther(ev.run, ev.ev, ev.vx, ev.vy, ev.vz, -1) != nm) cnt(kLookupMismatch);
        }
        int nbest = 0; double minchi2 = 1e300, bestchi2 = -1.;
        for (auto& h : ev.hyps) {
            spec[0][tr][tp].fill(h.mass);
            if (h.best == 1) { spec[1][tr][tp].fill(h.mass); ++nbest; bestchi2 = h.chi2; }
            if (h.chi2 < minchi2) minchi2 = h.chi2;
            if (h.mass > kMidLo && h.mass < kMidHi) {
                const Leg& a = ev.legs[h.i]; const Leg& b = ev.legs[h.j];
                h_D_nat[tr]->Fill(a.delta - b.delta);
                h_chi2_nat[tr]->Fill(std::log10(std::max(h.chi2, 1e-9)));
                if (h.best == 1) h_chi2_best[tr]->Fill(std::log10(std::max(h.chi2, 1e-9)));
            }
        }
        if (nbest != 1) cnt(kBestCountBad);   // ~1 per 2M events: PAT hpidpool.cc marks the best hypothesis by chi2 -> -chi2 and tests < 0, so chi2 == 0 (D == 0 in float) is never flagged
        else if (std::fabs(bestchi2 - minchi2) > 1e-6 * std::max(1., minchi2)) cnt(kBestNotMin);
        // emulation
        const int R = tr == 0 ? R_pt3 : R_pt2;
        if (!pool[tr].empty() && R > 0) {
            for (int k = 0; k < n; ++k) {
                const int q = ev.legs[k].q;
                for (int r = 0; r < R; ++r) {
                    const Event& B = pool[tr][rnd.Integer(pool[tr].size())];
                    int nq = 0;
                    for (const auto& l : B.legs) if (l.q == q) ++nq;
                    if (nq == 0) { cnt(kNoDonor); continue; }
                    int pick = rnd.Integer(nq), bi = -1;
                    for (int l = 0; l < int(B.legs.size()); ++l) if (B.legs[l].q == q && pick-- == 0) { bi = l; break; }
                    const Leg& b = B.legs[bi];
                    // frame of A from the remaining legs (a_k is in no cross pair): for
                    // inverse-variance weights the centred residuals of the legs that
                    // enter the mean are uncorrelated with the mean, so this is exact
                    // for Gaussian legs and removes the pull of a replaced outlier leg.
                    // The donor leg b must stay in the frame of B (a leg centred on a
                    // mean that excludes it is anti-correlated with that mean).
                    double Wk = 0., numk = 0.;
                    for (int i = 0; i < n; ++i) { if (i == k) continue; Wk += 1. / ev.legs[i].sig2; numk += ev.legs[i].delta / ev.legs[i].sig2; }
                    const double framek = numk / Wk;
                    const double c = rnd.Gaus(0., std::sqrt(1. / Wk + 1. / B.W));
                    double bestC = 1e300; bool bestIsCross = false; double bestMass = 0.;
                    for (const auto& h : ev.hyps) {
                        if (h.i == k || h.j == k) continue;
                        if (h.chi2 < bestC) { bestC = h.chi2; bestIsCross = false; }
                    }
                    for (int i = 0; i < n; ++i) {
                        if (i == k) continue;
                        const Leg& a = ev.legs[i];
                        if (pat == "pm" && a.q == q) continue;     // e+e-: partners of the opposite charge only
                        const double D   = (a.delta - framek) - b.ehat + c;
                        const double chi = patChi2(D, a.res, b.res);
                        const double tna = a.tofmom + D / 2., tnb = b.tofmom - D / 2.;
                        const double bna = a.L / (kC * tna), bnb = b.L / (kC * tnb);
                        const double m   = pairMass(a, b);
                        if (!(pidOk(a.q, a.p, bna) && pidOk(b.q, b.p, bnb))) { spec[4][tr][tp].fill(m); cnt(kNoPid); continue; }
                        spec[2][tr][tp].fill(m); cnt(kCrossActive);
                        if (m > kMidLo && m < kMidHi) { h_D_cross[tr]->Fill(D); h_chi2_cross[tr]->Fill(std::log10(std::max(chi, 1e-9))); }
                        if (chi < bestC) { bestC = chi; bestIsCross = true; bestMass = m; }
                    }
                    if (bestIsCross) { spec[3][tr][tp].fill(bestMass); cnt(kCrossWin); }
                }
            }
        }
        // pool
        Event keep; keep.trig = tr; keep.W = ev.W; keep.legs = ev.legs;
        pool[tr].push_back(std::move(keep));
        if (int(pool[tr].size()) > pool_depth) pool[tr].pop_front();
    };

    // ---- pass 1: stream rows, group events ----------------------------------
    Event cur; Float_t curRun = -1.f, curEv = -1.f; bool have = false;
    auto addLeg = [&](Event& ev, const LegBr& b, int q) -> int {
        for (int l = 0; l < int(ev.legs.size()); ++l) {
            const Leg& x = ev.legs[l];
            if (x.q == q && x.p == b.p && x.th == b.th && x.ph == b.ph) {
                if (std::fabs(x.tof - b.tof) > 1e-4 || std::fabs(x.tofmom - b.tofmom) > 1e-4) cnt(kLegDupInconsistent);
                return l;
            }
        }
        Leg l; l.q = q; l.p = b.p; l.th = b.th; l.ph = b.ph; l.pc = b.p; l.tof = b.tof; l.tofmom = b.tofmom; l.L = b.L;   // mass from p (FAT m_ee, KinematicType::RECONSTRUCTED)
        l.res = b.res; l.sys = int(b.sys); l.delta = l.tof - l.tofmom; l.cat = category(l.sys, l.p); l.sig2 = sig2[l.cat];
        ev.legs.push_back(l);
        return int(ev.legs.size()) - 1;
    };
    Long64_t nPrint = std::max<Long64_t>(1, nRows / 10);
    for (Long64_t i = 0; i < nRows; ++i) {
        t->GetEntry(i);
        // Inactive hypotheses (isBest = -1) are not written in the betanew
        // productions; the SMASH sim productions still contain them and FAT
        // drops them, so skip them here too (no-op for betanew data).
        if (v_best < 0.f) continue;
        cnt(kRows);
        if (!have || v_run != curRun || v_ev != curEv) {
            if (have) { processEvent(cur); }
            cur = Event(); curRun = v_run; curEv = v_ev; have = true; cur.run = v_run; cur.ev = v_ev; cur.vx = v_vx; cur.vy = v_vy; cur.vz = v_vz;
            if (sim_trigger_bits) {
                // SMASH sim: bit 13 = PT3, bit 12 = PT2, both set for PT3 events
                // (no downscale). Class 0 = PT3, class 1 = PT2 without PT3.
                const int tb = int(v_trig);
                // trigbit == 0: production without trigger simulation -> all events in slot 0
                cur.trig = (tb == 0 || (tb & 8192)) ? 0 : ((tb & 4096) ? 1 : -1);
            } else {
                cur.trig = (v_trig == 8192.f) ? 0 : (v_trig == 4096.f ? 1 : -1);
            }
            cur.vtxOk = (v_vz > vertex_z_min);
            if (pat == "pm") { cur.nOther = 0; cur.nSameLookup = -1; }
            else {
                cur.nOther = lookupOther(v_run, v_ev, v_vx, v_vy, v_vz, -q1);
                cur.nSameLookup = lookupOther(v_run, v_ev, v_vx, v_vy, v_vz, q1);
                if (cur.nOther < 0) { cur.nOther = 0; cur.nSameLookup = -1; cnt(kLookupMiss); }
            }
            cnt(kEvents);
            cnt(cur.trig == 0 ? kEvPT3 : (cur.trig == 1 ? kEvPT2 : kEvOther));
        }
        // native checks (PAT formulas on the stored legs)
        {
            const double D = (b1.tof - b2.tof) - (b1.tofmom - b2.tofmom);
            const double chi = patChi2(D, b1.res, b2.res);
            if (v_chi2 > 1e-5 && std::fabs(chi - v_chi2) > 1e-3 * v_chi2) cnt(kChi2Mismatch);
            const double bn1 = b1.L / (kC * (b1.tofmom + D / 2.)), bn2 = b2.L / (kC * (b2.tofmom - D / 2.));
            if (std::fabs(bn1 - b1.bnew) > 2e-4 || std::fabs(bn2 - b2.bnew) > 2e-4) cnt(kBetaMismatch);
        }
        Hyp h;
        h.i = addLeg(cur, b1, q1); h.j = addLeg(cur, b2, q2);
        h.chi2 = v_chi2; h.best = int(v_best);
        h.mass = pairMass(cur.legs[h.i], cur.legs[h.j]);
        // Repeated hypotheses (sim productions: (2,2) pairs written twice,
        // like-sign pairs as (a,b) and (b,a)): keep the first copy only.
        // The betanew exp production writes every combination once (no-op).
        bool repeated = false;
        for (const Hyp& o : cur.hyps)
            if ((o.i == h.i && o.j == h.j) || (o.i == h.j && o.j == h.i)) { repeated = true; break; }
        if (repeated) cnt(kHypRepeated);
        else cur.hyps.push_back(h);
        if ((i + 1) % nPrint == 0) printf("   ... %lld rows (%.0f s)\n", i + 1, sw.RealTime()), sw.Continue();
    }
    if (have) processEvent(cur);

    // ---- output -------------------------------------------------------------
    TFile* fout = new TFile(out_file, "recreate");
    for (int k = 0; k < 5; ++k) for (int tr = 0; tr < kNTrig; ++tr) for (int c = 0; c < kNCl; ++c) { spec[k][tr][c].h58->Write(); spec[k][tr][c].hf->Write(); }
    for (int tr = 0; tr < kNTrig; ++tr) { h_D_nat[tr]->Write(); h_D_cross[tr]->Write(); h_chi2_nat[tr]->Write(); h_chi2_best[tr]->Write(); h_chi2_cross[tr]->Write(); h_nlegs[tr]->Write(); h_ntot[tr]->Write(); }
    h_sigma->Write(); h_cnt->Write();
    TNamed info("info", Form("isbest_emulation_betanew: file=%s pattern=%s R_pt3=%d R_pt2=%d pool=%d seed=%d max_rows=%lld cuts=%s vertex_z_min=%g sim_trigger_bits=%d vertex_key=%d",
                             pat_file, pattern, R_pt3, R_pt2, pool_depth, seed, max_rows, cut_file, vertex_z_min, int(sim_trigger_bits), int(vertex_key)));
    info.Write();
    fout->Close();
    printf("[isbest_emu] counters:");
    for (int c = 1; c < kNCounters; ++c) printf("  %s=%.0f", cntName[c], h_cnt->GetBinContent(c));
    printf("\n[isbest_emu] written %s  (%.0f s)\n", out_file, sw.RealTime());
}
