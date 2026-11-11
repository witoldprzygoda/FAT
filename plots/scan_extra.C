// ============================================================================
// scan_extra.C  --  additional cuts for the two eta structures in M(pi+pi-e+e-)
// ============================================================================
// HADES p+p, T = 4.53 GeV, channel pi+ pi- e+ e-. Uses scan_compact.root
// (plots/scan_prepare.C) through plots/ScanEngine.h, plus a few columns read
// directly from scan_nt of the production (extra variables of section N).
//
// Start selections (all rows are particle combinations, isBest >= 0, corrected
// momenta, cor_ok == 1) = the kinematic nominal candidates of the two main scans:
//   peak  : oa_lab < 55 deg and M(pi+pi-) > 0.30 GeV/c2            (scan_peak.C)
//   bump  : oa_lab < 40 deg and 0.30 < M(pi+pi-) < 0.41 GeV/c2     (scan_bump.C)
// (First version: oa_lab < 50 with M(pipi) < 0.46 for the peak and < 0.42 for the
// bump. M(pipi) < 0.46 is flagged as sculpting by scan_peak.C and was removed
// after the review; starting from the nominal candidates also measures the
// interplay of the extra cuts with the kinematic ones.)
//
// Corrections after the review, in short:
//   - the PROPOSED extra cuts are the two vertex cuts (vz > -160 mm, both pions
//     used in the vertex fit); they are confirmed by the photon-tagged
//     eta -> pi+pi-pi0 control sample. vz > -160 tightens a binding decision
//     (z > -500) and needs an explicit decision; the efficiency of isUsedVertex
//     must come from simulation; |z_pi - vz| < 40 mm is evaluated as the
//     geometric alternative.
//   - the pi0-photon veto is OPTIONAL and unverified (no control sample).
//   - the two-fold check now also re-chooses the cut FAMILIES on one fold.
//   - S/sigma_S (F) is statistical only with a frozen shape: it compares
//     selections, it is not a significance. CB = 2 sqrt(N++ N--), no k factor:
//     CB shares are estimates.
//   - "eta -> pi+pi-pi0" here always means: with an e+e- pair from the pi0,
//     Dalitz decay OR converted photon (the two cannot be told apart here).
//
// For every candidate cut (list a..n of the task) the macro measures, by
// pass/fail fits on top of the start selection:
//   peak  : eps_S = S_pass / (S_pass + S_fail) from the frozen-shape fit P_main,
//           eps_B = background (polynomial + CB) kept in the peak window,
//           G = eps_S / sqrt(eps_B) (a cut helps only if G > 1),
//           F = S / sigma_S, S/B and the CB share before and after
//   bump  : eps_tag = efficiency for eta -> pi+pi-pi0 tagged with an ECAL photon
//           (eta peak in the pi0-constrained M(pi+pi-e+e-gamma), see M5c),
//           the fraction of rows and of N - CB kept in the bump window,
//           the summed edge-step height of three M(pi+pi-) slices
//   both  : pass fraction against M4 (sculpting)
// Thresholds are never taken at the maximum of a scan: the working point is the
// centre of the grid region whose figure of merit is within 1 unit of the best.
// Two folds (parity of the event group): thresholds chosen on one fold are
// measured on the other.
//
// Output
//   tables   <prod>/scan/extra_peak_candidates.txt, extra_bump_candidates.txt (every candidate),
//            extra_sculpt.txt (pass fraction against M4), extra_ranking.txt, extra_vertex.txt,
//            extra_tagged.txt (pi0 tag), extra_sbsub.txt (sideband comparison of all variables),
//            extra_combination.txt (best cuts together, model spread, two-fold check)
//   histograms <prod>/scan/extra_scan.root
//   figures  plots/output/eta_scan/extra_*.png (.pdf for the main ones)
//
// Stages: Prepare (data sets) -> Evaluate (all candidates, pass/fail, folds) -> StageTables ->
//   StageN (sideband comparison) -> StageRank -> StageVertex -> StageTagged -> StageCombine ->
//   figures (StageBeforeAfter, StageSummary, StageScans, StageSculptFig, StageSbFig)
//
// Run (ACLiC, build directory outside the repository):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/scan/aclic_extra", kTRUE)' 'plots/scan_extra.C+'
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <functional>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TPad.h"
#include "TPave.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include "ScanEngine.h"

namespace sx {

using scan::Data;
using scan::PeakCfg;
using scan::PeakFit;
using scan::Row;
using scan::Sel;
using scan::Spectra;
using scan::StepCfg;
using scan::StepFit;
typedef const Row& R;

// ============================================================================
// constants of this study
// ============================================================================
const float kOaLabPeak = 55.f;             // peak start: angle (pi+pi-) to (e+e-) in the lab [deg]
const float kMppLoPeak = 0.30f;            // peak start: lower M(pi+pi-) [GeV/c2]
const float kOaLabBump = 40.f;             // bump start
const float kMppLoBump = 0.30f, kMppHiBump = 0.41f;
const char* const kStartPeakTxt = "oa_lab < 55, M(pipi) > 0.30 (kinematic nominal of scan_peak.C)";
const char* const kStartBumpTxt = "oa_lab < 40, 0.30 < M(pipi) < 0.41 (kinematic nominal of scan_bump.C)";
const char* const kStartPeakTex = "oa_{lab} < 55#circ, M(#pi#pi) > 0.30";
const char* const kStartBumpTex = "oa_{lab} < 40#circ, 0.30 < M(#pi#pi) < 0.41";
// proposed common lower cut on the vertex z [mm]. The binding preselection is z > -500:
// tightening it is a NEW decision of the person in charge.
const double kVzCommon = -160.;
const double kSlice[4] = {0.30, 0.34, 0.37, 0.39};   // M(pi+pi-) slices of the edge-step fit

// pi0 tag: photon whose M(e+e-gamma) is closest to 0.140 (column m_eeg). In the
// bump start selection the pi0 sits at 0.1428, sigma 13.3 MeV (Gaussian + pol3).
const float kPi0Lo = 0.110f, kPi0Hi = 0.170f;        // +-2.2 sigma, rounded
// eta peak in the pi0-constrained five-particle mass M5c (see M5c()):
// free Poisson fits, bump start selection, pol2, 0.49-0.62: 0.54637 +- 0.00022,
// 6.16 +- 0.25 MeV (pol3: 0.54625, 6.20). Frozen:
const double kM5Mean = 0.5463, kM5Sigma = 0.0061;
const double kM5Lo = 0.49, kM5Hi = 0.62;
const int kM5Order = 2;
const int kM5Nbin = 160;                   // 2.5 MeV bins, 0.40-0.80
const double kM5HLo = 0.40, kM5HHi = 0.80;
const int kM5Reb = 8;                      // CB of the tagged sample in 20 MeV bins
// beam spot (medians of the vertex x, y of the start selection) [mm]
const double kVx0 = 1.19, kVy0 = -2.39;
const double kPionM2 = 19480.;             // m_pi^2 [MeV^2]

const char* const kTabDir = scan::kScanDir;
const char* kXt = "M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]";
const char* kYt = "rows / 2.5 MeV/c^{2}";

// ---- printing: screen and the current table file -----------------------------
FILE* gOut = nullptr;
void P(const char* fmt, ...) {
    va_list a, b;
    va_start(a, fmt);
    va_copy(b, a);
    vprintf(fmt, a);
    if (gOut) vfprintf(gOut, fmt, b);
    va_end(a);
    va_end(b);
}
void OpenTable(const char* name) {
    if (gOut) fclose(gOut);
    gOut = fopen(Form("%s/%s", kTabDir, name), "w");
    P("# %s  --  written by plots/scan_extra.C\n", name);
    P("# data: %s ; rows = particle combinations (isBest >= 0), corrected momenta, cor_ok == 1\n", scan::kCompact);
}
void CloseTable() { if (gOut) fclose(gOut); gOut = nullptr; }

// ============================================================================
// pi0-constrained five-particle mass
// ============================================================================
// The photon energy is rescaled by lambda so that M(e+e-gamma) = m_pi0:
//   lambda = (m_pi0^2 - m_ee^2) / (M_eeg^2 - m_ee^2)
//   M5c^2  = M4^2 + lambda (M5^2 - M4^2)
// (exact for a massless photon: all terms with the photon are linear in its
// energy). This removes the ECAL energy resolution from the eta peak: width
// 6.1 MeV instead of 9-11 MeV, and the position no longer depends on the pi0
// window.
inline double M5c(R r) {
    if (!(r.m_eeg > 0.f) || !(r.m5 > 0.f)) return -1.;
    const double den = (double)r.m_eeg * r.m_eeg - (double)r.mee * r.mee;
    if (!(den > 0.)) return -1.;
    const double lam = (scan::kMpi0 * scan::kMpi0 - (double)r.mee * r.mee) / den;
    if (!(lam > 0.)) return -1.;
    const double m2 = (double)r.m4 * r.m4 + lam * ((double)r.m5 * r.m5 - (double)r.m4 * r.m4);
    return m2 > 0. ? std::sqrt(m2) : -1.;
}
inline bool Pi0Tag(R r) { return r.m_eeg >= kPi0Lo && r.m_eeg < kPi0Hi; }

PeakCfg Cfg5() {
    PeakCfg c;
    c.mean = kM5Mean; c.sigma = kM5Sigma; c.lo = kM5Lo; c.hi = kM5Hi; c.order = kM5Order;
    c.winLo = kM5Mean - scan::kWinNsig * kM5Sigma; c.winHi = kM5Mean + scan::kWinNsig * kM5Sigma;
    return c;
}
scan::FillOpt Opt5() {
    scan::FillOpt o;
    o.nbin = kM5Nbin; o.lo = kM5HLo; o.hi = kM5HHi; o.reb = kM5Reb;
    o.x = [](R r) { return M5c(r); };
    const PeakCfg c = Cfg5();
    o.windows.push_back(std::make_pair(c.winLo, c.winHi));
    return o;
}

// ============================================================================
// extra columns read from scan_nt (not in the compact file)
// ============================================================================
#define SX_EXT_FIELDS(X)                                                          \
    X(trigbit) X(start_iteration)                                                 \
    X(ep_sector) X(em_sector) X(pip_sector) X(pim_sector)                         \
    X(ep_r) X(em_r) X(pip_r) X(pim_r) X(ep_z) X(em_z) X(pip_z) X(pim_z)           \
    X(ep_beta_new) X(em_beta_new) X(pip_beta_new) X(pim_beta_new)                 \
    X(ep_tof_new) X(em_tof_new) X(pip_tof_new) X(pim_tof_new)                     \
    X(ep_tof_mom) X(em_tof_mom) X(pip_tof_mom) X(pim_tof_mom)                     \
    X(ep_theta_rich) X(em_theta_rich) X(ep_phi_rich) X(em_phi_rich)               \
    X(ep_rich_padnum) X(em_rich_padnum) X(g1_e) X(g2_e) X(fw1_phi)
struct Ext {
#define X(n) Float_t n;
    SX_EXT_FIELDS(X)
#undef X
};

// ---- data sets ---------------------------------------------------------------
// index f: 0 = all rows, 1 = fold 0, 2 = fold 1 (fold = parity of the event group)
Data gAll;                 // every row of the compact file
Data gP[3];                // peak start selection
Data gB[3][4];             // bump start selection: all of it, then the three M(pi+pi-) slices
Data gTagAll;              // every cor_ok row with a pi0 photon (tag sample without kinematic cuts)
Data gT[3];                // bump start selection and a pi0 photon
std::vector<Ext> gExt[3];  // extra columns of the rows of gP[0], per sample
// In the reduced sets the columns file_id and pat_entry are OVERWRITTEN in
// memory with (sample, index into gExt): this addresses the extra columns.
inline const Ext& E(R r) { return gExt[r.file_id][r.pat_entry]; }

inline bool StartPeak(R r) { return r.cor_ok && r.oa_lab < kOaLabPeak && r.mpipi > kMppLoPeak; }
// the bump start is a subset of the peak start
inline bool StartBump(R r) { return StartPeak(r) && r.oa_lab < kOaLabBump && r.mpipi > kMppLoBump && r.mpipi < kMppHiBump; }

bool Prepare() {
    TStopwatch sw;
    if (!scan::Load(gAll)) return false;
    const char* inName[3] = {"output_pippimepem_exp.root", "output_pippimepep_exp.root",
                             "output_pippimemem_exp.root"};
    Long64_t bad = 0;
    for (int k = 0; k < 3; ++k) {
        std::vector<Long64_t> orig;
        for (size_t i = 0; i < gAll.s[k].size(); ++i) {
            const Row& r = gAll.s[k][i];
            if (!StartPeak(r)) continue;
            gP[0].s[k].push_back(r);
            orig.push_back((Long64_t)i);
        }
        // extra columns: row i of the compact tree = entry i of scan_nt (checked with the keys)
        TFile* f = TFile::Open(Form("%s/%s", scan::kProdDir, inName[k]), "READ");
        if (!f || f->IsZombie()) { printf("cannot open %s\n", inName[k]); return false; }
        TTree* t = (TTree*)f->Get("scan_nt");
        if (!t) return false;
        Ext e;
        Float_t pat = 0, fid = 0;
        t->SetBranchStatus("*", 0);
#define X(n) t->SetBranchStatus(#n, 1); t->SetBranchAddress(#n, &e.n);
        SX_EXT_FIELDS(X)
#undef X
        t->SetBranchStatus("pat_entry", 1); t->SetBranchAddress("pat_entry", &pat);
        t->SetBranchStatus("file_id", 1); t->SetBranchAddress("file_id", &fid);
        t->SetCacheSize(200 * 1024 * 1024);
        gExt[k].clear();
        gExt[k].reserve(orig.size());
        for (size_t j = 0; j < orig.size(); ++j) {
            t->GetEntry(orig[j]);
            Row& r = gP[0].s[k][j];
            if ((int)std::lround(pat) != r.pat_entry || (int)std::lround(fid) != r.file_id) bad++;
            gExt[k].push_back(e);
            r.file_id = k;
            r.pat_entry = (int)j;
        }
        f->Close();
        delete f;
    }
    if (bad) { printf("Prepare: %lld rows do not match scan_nt -- stop\n", bad); return false; }
    for (int k = 0; k < 3; ++k) {
        for (const auto& r : gP[0].s[k]) {
            gP[1 + r.fold].s[k].push_back(r);
            if (StartBump(r)) {
                int sl = 0;
                for (int s = 0; s < 3; ++s) if (r.mpipi >= kSlice[s] && r.mpipi < kSlice[s + 1]) sl = s + 1;
                for (int f : {0, 1 + r.fold}) {
                    gB[f][0].s[k].push_back(r);
                    if (sl) gB[f][sl].s[k].push_back(r);
                    if (Pi0Tag(r)) gT[f].s[k].push_back(r);
                }
            }
        }
    }
    for (int k = 0; k < 3; ++k)
        for (const auto& r : gAll.s[k]) if (r.cor_ok && Pi0Tag(r)) gTagAll.s[k].push_back(r);
    printf("Prepare: peak start rows %zu / %zu / %zu, bump start %zu / %zu / %zu, tagged %zu / %zu / %zu (us/pp/mm), %.1f s\n",
           gP[0].s[0].size(), gP[0].s[1].size(), gP[0].s[2].size(),
           gB[0][0].s[0].size(), gB[0][0].s[1].size(), gB[0][0].s[2].size(),
           gT[0].s[0].size(), gT[0].s[1].size(), gT[0].s[2].size(), sw.RealTime());
    return true;
}

// ============================================================================
// candidates
// ============================================================================
inline float Mx(float a, float b) { return a > b ? a : b; }
inline float Mn(float a, float b) { return a < b ? a : b; }
inline double RT(R r) { return std::hypot(r.vx - kVx0, r.vy - kVy0); }

struct Cand {
    std::string grp;     // letter of the task list (a..n)
    std::string fam;     // scan family ("" = single cut)
    std::string name;    // printed label
    double thr = 0;      // threshold (scan families)
    Sel cut;
};
std::vector<Cand> gC;
void Add(const char* grp, const char* name, const Sel& s) {
    Cand c; c.grp = grp; c.name = name; c.cut = s; gC.push_back(c);
}
// thresholds are given from the loosest to the tightest
void AddScan(const char* grp, const char* fam, const char* fmt, const std::vector<double>& thr,
             const std::function<Sel(double)>& mk) {
    for (double t : thr) {
        Cand c; c.grp = grp; c.fam = fam; c.name = Form(fmt, t); c.thr = t; c.cut = mk(t);
        gC.push_back(c);
    }
}
// lepton detector class: 0 both RPC, 1 both TOF, 2 one RPC + one TOF, 3 any lepton with system 2
inline int LepSys(R r) {
    if (r.ep_system == 2 || r.em_system == 2) return 3;
    if (r.ep_system == 0 && r.em_system == 0) return 0;
    if (r.ep_system == 1 && r.em_system == 1) return 1;
    return 2;
}

void BuildCandidates() {
    gC.clear();
    Add("-", "no cut", [](R) { return true; });
    // ---- a. missing mass of the four particles ----
    AddScan("a", "MM4<", "MM4 < %.2f", {2.95, 2.90, 2.85, 2.80, 2.75, 2.70, 2.65, 2.60},
            [](double t) { const float c = t; return Sel([c](R r) { return r.mm4 < c; }); });
    AddScan("a", "MM4>", "MM4 > %.3f", {1.800, 1.876, 1.950, 2.000, 2.100, 2.200},
            [](double t) { const float c = t; return Sel([c](R r) { return r.mm4 > c; }); });
    // ---- b. vertex ----
    AddScan("b", "vz>", "vz > %.0f", {-200, -190, -180, -170, -165, -160, -155, -150, -140},
            [](double t) { const float c = t; return Sel([c](R r) { return r.vz > c; }); });
    AddScan("b", "vz<", "vz < %.0f", {0, -40, -60, -70, -80, -90, -100, -110},
            [](double t) { const float c = t; return Sel([c](R r) { return r.vz < c; }); });
    AddScan("b", "rT<", "rT < %.1f", {40, 30, 20, 15, 10, 7.5, 5},
            [](double t) { return Sel([t](R r) { return RT(r) < t; }); });
    AddScan("b", "vchi2<", "vchi2 < %.0f", {1000, 500, 300, 200, 150, 100, 75, 50},
            [](double t) { const float c = t; return Sel([c](R r) { return r.vchi2 < c; }); });
    // ---- c. pion with a RICH ring: veto ----
    Add("c", "no pion with ring", [](R r) { return r.pip_isring == 0 && r.pim_isring == 0; });
    Add("c", "pi+ without ring", [](R r) { return r.pip_isring == 0; });
    Add("c", "pi- without ring", [](R r) { return r.pim_isring == 0; });
    // ---- d. shared RICH ring of the two leptons ----
    Add("d", "shared ring required", [](R r) { return r.shared_ring == 1; });
    Add("d", "shared ring veto", [](R r) { return r.shared_ring == 0; });
    // ---- e. crossing points ----
    AddScan("e", "dz_ee<", "dz_ee < %.0f", {800, 600, 400, 300, 200, 150, 100, 50},
            [](double t) { const float c = t; return Sel([c](R r) { return r.dz_ee < c; }); });
    AddScan("e", "dz_ee>", "dz_ee > %.0f", {-400, -200, -100, -50, -20},
            [](double t) { const float c = t; return Sel([c](R r) { return r.dz_ee > c; }); });
    AddScan("e", "rv_ee<", "rv_ee < %.0f", {600, 400, 300, 200, 100, 50, 25},
            [](double t) { const float c = t; return Sel([c](R r) { return r.rv_ee < c; }); });
    AddScan("e", "doca_ee<", "doca_ee < %.1f", {40, 30, 20, 15, 10, 7.5, 5},
            [](double t) { const float c = t; return Sel([c](R r) { return r.doca_ee < c; }); });
    AddScan("e", "|dz_pipi|<", "|dz_pipi| < %.0f", {200, 100, 50, 30, 20, 10},
            [](double t) { const float c = t; return Sel([c](R r) { return std::fabs(r.dz_pipi) < c; }); });
    AddScan("e", "rv_pipi<", "rv_pipi < %.0f", {150, 60, 40, 25, 15, 10},
            [](double t) { const float c = t; return Sel([c](R r) { return r.rv_pipi < c; }); });
    AddScan("e", "doca_pipi<", "doca_pipi < %.0f", {100, 50, 30, 20, 15, 10, 5},
            [](double t) { const float c = t; return Sel([c](R r) { return r.doca_pipi < c; }); });
    // ---- f. track-vertex association ----
    Add("f", "isUsedVertex all four", [](R r) {
        return r.ep_isUsedVertex && r.em_isUsedVertex && r.pip_isUsedVertex && r.pim_isUsedVertex; });
    Add("f", "isUsedVertex both leptons", [](R r) { return r.ep_isUsedVertex && r.em_isUsedVertex; });
    Add("f", "isUsedVertex both pions", [](R r) { return r.pip_isUsedVertex && r.pim_isUsedVertex; });
    // ---- g. lepton and pion quality ----
    Add("g", "kIsLepton both", [](R r) { return r.ep_kIsLepton && r.em_kIsLepton; });
    Add("g", "kIsLepton at least one", [](R r) { return r.ep_kIsLepton || r.em_kIsLepton; });
    AddScan("g", "rmq<", "richmatchq (max) < %.1f", {5, 4, 3, 2.5, 2, 1.5, 1},
            [](double t) { const float c = t; return Sel([c](R r) {
                return Mx(r.ep_richmatchquality, r.em_richmatchquality) < c; }); });
    AddScan("g", "ringchi2<", "ringchi2 (max) < %.0f", {600, 500, 400, 350, 300, 250, 200},
            [](double t) { const float c = t; return Sel([c](R r) { return Mx(r.ep_ringchi2, r.em_ringchi2) < c; }); });
    AddScan("g", "|m2pi|<", "|mass2 - m_pi^2| (both pi) < %.0f", {60000, 30000, 20000, 15000, 10000, 7500, 5000},
            [](double t) { const float c = t; return Sel([c](R r) {
                return std::fabs(r.pip_mass2 - kPionM2) < c && std::fabs(r.pim_mass2 - kPionM2) < c; }); });
    AddScan("g", "dedx<", "dE/dx MDC (max pi) < %.1f", {6, 5, 4, 3.5, 3, 2.5, 2.2},
            [](double t) { const float c = t; return Sel([c](R r) { return Mx(r.pip_dedx_mdc, r.pim_dedx_mdc) < c; }); });
    AddScan("g", "dedx>", "dE/dx MDC (min pi) > %.1f", {1.2, 1.3, 1.4, 1.5},
            [](double t) { const float c = t; return Sel([c](R r) { return Mn(r.pip_dedx_mdc, r.pim_dedx_mdc) > c; }); });
    // ---- h. START ----
    Add("h", "evtStart == 1", [](R r) { return r.evtStart == 1; });
    // ---- i. transverse momentum and rapidity of the four particles ----
    AddScan("i", "pt4>", "pt4 > %.0f", {50, 100, 150, 200, 250, 300, 350, 400, 450},
            [](double t) { const float c = t; return Sel([c](R r) { return r.pt4 > c; }); });
    AddScan("i", "y4>", "y4 > %.1f", {-0.8, -0.7, -0.6, -0.5, -0.4, -0.3},
            [](double t) { const float c = t; return Sel([c](R r) { return r.y4 > c; }); });
    AddScan("i", "y4<", "y4 < %.1f", {0.5, 0.4, 0.3, 0.2, 0.1, 0.0},
            [](double t) { const float c = t; return Sel([c](R r) { return r.y4 < c; }); });
    // ---- j. event classes from the adjacency groups ----
    Add("j", "single-row event", [](R r) { return r.nrow == 1; });
    Add("j", "multi-row event", [](R r) { return r.nrow > 1; });
    Add("j", "nrow <= 2", [](R r) { return r.nrow <= 2; });
    Add("j", "one pi+ and one pi-", [](R r) { return r.n_pip == 1 && r.n_pim == 1; });
    Add("j", "two leptons only", [](R r) { return r.n_lep == 2; });
    // ---- k. photons ----
    Add("k", "n_gam == 0", [](R r) { return r.n_gam == 0; });
    Add("k", "n_gam >= 1", [](R r) { return r.n_gam >= 1; });
    Add("k", "no pi0 photon", [](R r) { return !Pi0Tag(r); });
    Add("k", "neutr_mult == 0", [](R r) { return r.neutr_mult == 0; });
    // ---- l. ranking chi2 ----
    AddScan("l", "hyp_chi2<", "hyp_chi2 < %.2f", {8, 5, 3, 2, 1.5, 1.0, 0.75, 0.5, 0.3},
            [](double t) { const float c = t; return Sel([c](R r) { return r.hyp_chi2 < c; }); });
    // ---- m. lepton detector system, momentum, polar angle ----
    Add("m", "leptons both RPC", [](R r) { return LepSys(r) == 0; });
    Add("m", "leptons both TOF", [](R r) { return LepSys(r) == 1; });
    Add("m", "leptons RPC + TOF", [](R r) { return LepSys(r) == 2; });
    Add("m", "no lepton with system 2", [](R r) { return LepSys(r) != 3; });
    AddScan("m", "p_lep>", "p lepton (min) > %.0f", {50, 60, 70, 80, 100, 120, 150},
            [](double t) { const float c = t; return Sel([c](R r) { return Mn(r.ep_p, r.em_p) > c; }); });
    AddScan("m", "p_lep<", "p lepton (max) < %.0f", {1000, 800, 600, 500, 400, 300},
            [](double t) { const float c = t; return Sel([c](R r) { return Mx(r.ep_p, r.em_p) < c; }); });
    AddScan("m", "th_lep>", "theta lepton (min) > %.0f", {12, 15, 18, 20, 22, 25},
            [](double t) { const float c = t; return Sel([c](R r) { return Mn(r.ep_theta, r.em_theta) > c; }); });
    AddScan("m", "th_lep<", "theta lepton (max) < %.0f", {80, 75, 70, 65, 60, 55, 50},
            [](double t) { const float c = t; return Sel([c](R r) { return Mx(r.ep_theta, r.em_theta) < c; }); });
    // ---- n. the three variables with the largest separation in the sideband comparison (stage N,
    //         extra_sbsub.txt): isUsedVertex of the pions (= group f), momentum of the forward-detector
    //         candidate, distance in z between the pion track points and the vertex ----
    Add("n", "isUsedVertex at least one pion", [](R r) { return r.pip_isUsedVertex || r.pim_isUsedVertex; });
    Add("n", "forward candidate (fwdet_mult >= 1)", [](R r) { return r.fwdet_mult >= 1; });
    Add("n", "no forward candidate", [](R r) { return r.fwdet_mult == 0; });
    Add("n", "forward candidate p > 1900", [](R r) { return r.fw1_p > 1900.f; });
    AddScan("n", "|z-vz|pi<", "|z - vz| (both pi) < %.0f", {150, 100, 60, 40, 30, 20, 15, 10},
            [](double t) { const float c = t; return Sel([c](R r) {
                return std::fabs(E(r).pip_z - r.vz) < c && std::fabs(E(r).pim_z - r.vz) < c; }); });
}
int Idx(const std::string& name) {
    for (size_t i = 0; i < gC.size(); ++i) if (gC[i].name == name) return (int)i;
    printf("Idx: candidate '%s' not found\n", name.c_str());
    return 0;
}

// ============================================================================
// evaluation of a list of candidates (pass and fail samples, three fold sets)
// ============================================================================
struct Half {
    // peak, P_main on the peak start selection
    double S = 0, eS = 0, chi2 = 0, prob = 0; int ndf = 0;
    double N = 0, CB = 0, B = 0, Swin = 0, SoB = 0;
    // bump window 0.4550-0.5225, bump start selection
    double bN = 0, bVarN = 0, bCB = 0, bVarCB = 0;
    // edge-step heights in the three M(pi+pi-) slices [rows / 2.5 MeV]
    double h[3] = {0, 0, 0}, eh[3] = {0, 0, 0}; int hok[3] = {0, 0, 0};
    // tagged eta (Poisson fit of the unsubtracted M5c spectrum)
    double S5 = 0, eS5 = 0, N5 = 0, CB5 = 0, SoB5 = 0, chi5 = 0; int ndf5 = 0;
    double F() const { return eS > 0 ? S / eS : 0; }
    double hSum() const { return h[0] + h[1] + h[2]; }
    double hErr() const { return std::sqrt(eh[0] * eh[0] + eh[1] * eh[1] + eh[2] * eh[2]); }
    double Z() const { return hErr() > 0 ? hSum() / hErr() : 0; }     // summed step height / error
    double F5() const { return eS5 > 0 ? S5 / eS5 : 0; }
};
struct CR {
    Half h[3][2];                  // [0 all, 1 fold 0, 2 fold 1][0 pass, 1 fail]
    std::vector<double> pP, pF;    // unlike-sign rows per 10 MeV bin (0.35-0.95), peak start: pass, fail
    std::vector<double> bP, bF;    // the same for the bump start
};

void FillPeak(Half& H, const Spectra& sp) {
    const PeakFit f = scan::FitPeakMain(sp);
    H.S = f.S; H.eS = f.eS; H.chi2 = f.chi2; H.ndf = f.ndf; H.prob = f.prob;
    H.N = f.N; H.CB = f.CB; H.B = f.B; H.Swin = f.Swin; H.SoB = f.SoB;
}
void FillBumpWin(Half& H, const Spectra& sp) {
    const scan::WinSum& w = sp.win[1];
    H.bN = w.n[0]; H.bVarN = w.var[0]; H.bCB = w.cb; H.bVarCB = w.cbVar;
}
void FillStep(Half& H, const Spectra& sp, int s) {
    StepCfg c; c.mpLo = kSlice[s]; c.mpHi = kSlice[s + 1];
    const StepFit f = scan::FitEdgeStep(sp, c);
    H.hok[s] = f.ok; H.h[s] = f.ok ? f.h : 0; H.eh[s] = f.ok ? f.eh : 0;
}
void FillTag(Half& H, const Spectra& sp) {
    const PeakFit f = scan::FitPeakPois(sp, Cfg5());
    H.S5 = f.ok ? f.S : 0; H.eS5 = f.ok ? f.eS : 0; H.N5 = f.N; H.CB5 = f.CB; H.SoB5 = f.SoB;
    H.chi5 = f.chi2; H.ndf5 = f.ndf;
}

std::vector<CR> Evaluate(const std::vector<Cand>& c, bool folds = true) {
    TStopwatch sw;
    std::vector<CR> out(c.size());
    const size_t chunk = 40;
    const scan::FillOpt o5 = Opt5();
    const int nf = folds ? 3 : 1;
    for (size_t i0 = 0; i0 < c.size(); i0 += chunk) {
        const size_t i1 = std::min(c.size(), i0 + chunk);
        std::vector<Sel> sels;
        for (size_t i = i0; i < i1; ++i) {
            const Sel cut = c[i].cut;
            sels.push_back(cut);
            sels.push_back([cut](R r) { return !cut(r); });
        }
        for (int f = 0; f < nf; ++f) {
            {
                const std::vector<Spectra> sp = scan::FillMany(gP[f], sels);
                for (size_t i = i0; i < i1; ++i)
                    for (int q = 0; q < 2; ++q) {
                        const Spectra& s = sp[2 * (i - i0) + q];
                        FillPeak(out[i].h[f][q], s);
                        if (f == 0) (q == 0 ? out[i].pP : out[i].pF) = s.w[0];
                    }
            }
            {
                const std::vector<Spectra> sp = scan::FillMany(gB[f][0], sels);
                for (size_t i = i0; i < i1; ++i)
                    for (int q = 0; q < 2; ++q) {
                        const Spectra& s = sp[2 * (i - i0) + q];
                        FillBumpWin(out[i].h[f][q], s);
                        if (f == 0) (q == 0 ? out[i].bP : out[i].bF) = s.w[0];
                    }
            }
            for (int sl = 0; sl < 3; ++sl) {
                const std::vector<Spectra> sp = scan::FillMany(gB[f][1 + sl], sels);
                for (size_t i = i0; i < i1; ++i)
                    for (int q = 0; q < 2; ++q) FillStep(out[i].h[f][q], sp[2 * (i - i0) + q], sl);
            }
            {
                const std::vector<Spectra> sp = scan::FillMany(gT[f], sels, o5);
                for (size_t i = i0; i < i1; ++i)
                    for (int q = 0; q < 2; ++q) FillTag(out[i].h[f][q], sp[2 * (i - i0) + q]);
            }
        }
    }
    printf("Evaluate: %zu candidates in %.1f s\n", c.size(), sw.RealTime());
    return out;
}

// ---- derived numbers -----------------------------------------------------------
struct Eff { double e = 0, err = 0; };
// a / (a + b) for independent a, b
Eff Ratio(double a, double ea, double b, double eb) {
    Eff r; const double t = a + b;
    if (!(t > 0)) return r;
    r.e = a / t; r.err = std::sqrt(b * b * ea * ea + a * a * eb * eb) / (t * t);
    return r;
}
struct PeakNum {
    Eff epsS; double epsB = 0, G = 0, eG = 0, pull = 0;     // G = eps_S / sqrt(eps_B)
    double F = 0, SoB = 0, cb = 0, closure = 0;             // after the cut; closure = S_pass + S_fail - S_all
};
PeakNum PeakOf(const CR& r, int f, const Half& base) {
    PeakNum p;
    const Half& a = r.h[f][0]; const Half& b = r.h[f][1];
    p.epsS = Ratio(a.S, a.eS, b.S, b.eS);
    const double ba = a.N - a.Swin, bb = b.N - b.Swin;
    p.epsB = (ba + bb) > 0 ? ba / (ba + bb) : 0;
    if (p.epsB > 0) { p.G = p.epsS.e / std::sqrt(p.epsB); p.eG = p.epsS.err / std::sqrt(p.epsB); }
    p.pull = p.eG > 0 ? (p.G - 1.) / p.eG : 0;
    p.F = a.F(); p.SoB = a.SoB; p.cb = a.N > 0 ? a.CB / a.N : 0;
    p.closure = a.S + b.S - base.S;
    return p;
}
struct BumpNum {
    Eff epsTag; double epsN = 0, epsNCB = 0, Gb = 0, eGb = 0, pull = 0;   // Gb = eps_tag / sqrt(eps_N)
    double cb = 0;                 // CB share in the bump window after the cut
    Eff epsH; double hP = 0, ehP = 0, Z = 0, F5 = 0, SoB5 = 0;
};
BumpNum BumpOf(const CR& r, int f, const Half& base) {
    BumpNum p;
    const Half& a = r.h[f][0]; const Half& b = r.h[f][1];
    p.epsTag = Ratio(a.S5, a.eS5, b.S5, b.eS5);
    p.epsN = (a.bN + b.bN) > 0 ? a.bN / (a.bN + b.bN) : 0;
    p.epsNCB = (base.bN - base.bCB) > 0 ? (a.bN - a.bCB) / (base.bN - base.bCB) : 0;
    if (p.epsN > 0) { p.Gb = p.epsTag.e / std::sqrt(p.epsN); p.eGb = p.epsTag.err / std::sqrt(p.epsN); }
    p.pull = p.eGb > 0 ? (p.Gb - 1.) / p.eGb : 0;
    p.cb = a.bN > 0 ? a.bCB / a.bN : 0;
    p.epsH = Ratio(a.hSum(), a.hErr(), b.hSum(), b.hErr());
    p.hP = a.hSum(); p.ehP = a.hErr(); p.Z = a.Z(); p.F5 = a.F5(); p.SoB5 = a.SoB5;
    return p;
}

// Working point of a threshold scan. fom[0] = no cut, fom[1..n] = thresholds from
// the loosest to the tightest. Region = contiguous points around the best one
// with fom >= best - 1; working point = its centre (the looser of the two middle
// points for an even number). Returns the index (0 = no cut); lo, hi = region.
int WorkPoint(const std::vector<double>& fom, int& lo, int& hi) {
    int jm = 0;
    for (size_t j = 1; j < fom.size(); ++j) if (fom[j] > fom[jm]) jm = (int)j;
    lo = hi = jm;
    while (lo > 0 && fom[lo - 1] >= fom[jm] - 1.) --lo;
    while (hi + 1 < (int)fom.size() && fom[hi + 1] >= fom[jm] - 1.) ++hi;
    return (lo + hi) / 2;
}

// scan families: name -> candidate indices in threshold order
std::vector<std::pair<std::string, std::vector<int>>> Families(const std::vector<Cand>& c) {
    std::vector<std::pair<std::string, std::vector<int>>> fam;
    for (size_t i = 0; i < c.size(); ++i) {
        if (c[i].fam.empty()) continue;
        if (fam.empty() || fam.back().first != c[i].fam) fam.push_back(std::make_pair(c[i].fam, std::vector<int>()));
        fam.back().second.push_back((int)i);
    }
    return fam;
}
// working point of a family on fold set f; kind 0: peak F = S/sigma_S, 1: bump Z (summed step / error)
// returns the candidate index, or -1 for "no cut"
int FamilyWP(const std::vector<CR>& res, const std::vector<int>& idx, const Half& base, int f, int kind,
             int* nreg = nullptr) {
    std::vector<double> fom;
    fom.push_back(kind == 0 ? base.F() : base.Z());
    for (int i : idx) fom.push_back(kind == 0 ? res[i].h[f][0].F() : res[i].h[f][0].Z());
    int lo, hi;
    const int j = WorkPoint(fom, lo, hi);
    if (nreg) *nreg = hi - lo + 1;
    return j == 0 ? -1 : idx[j - 1];
}

// ---- sculpting: pass fraction against M4 in 10 MeV bins -----------------------
struct Sculpt {
    double reg[8] = {0, 0, 0, 0, 0, 0, 0, 0};   // pass fraction in the regions kReg
    double slope = 0;      // relative change of the pass fraction per 100 MeV (line fitted on 0.52-0.62 without 0.53-0.56)
    double chi2 = 0; int ndf = 0;
    double dPeak = 0, eDPeak = 0;   // pass fraction in 0.53-0.56 minus the line
    double range = 0;      // (max - min) / mean of the pass fraction over 0.42-0.56
};
const double kReg[9] = {0.42, 0.46, 0.49, 0.52, 0.53, 0.56, 0.58, 0.61, 0.70};
Sculpt SculptOf(const std::vector<double>& p, const std::vector<double>& q) {
    Sculpt s;
    const double lo = scan::kMlo, bw = scan::kBinW * scan::kCBreb;
    auto ib = [&](double m) { return (int)std::floor((m - lo) / bw + 0.5); };
    for (int k = 0; k < 8; ++k) {
        double a = 0, b = 0;
        for (int i = ib(kReg[k]); i < ib(kReg[k + 1]); ++i) { a += p[i]; b += q[i]; }
        s.reg[k] = (a + b) > 0 ? a / (a + b) : 0;
    }
    // straight line through the side regions
    double S0 = 0, Sx = 0, Sxx = 0, Sy = 0, Sxy = 0;
    std::vector<double> xs, ys, es;
    for (int i = ib(0.52); i < ib(0.62); ++i) {
        const double m = lo + (i + 0.5) * bw;
        if (m > 0.53 && m < 0.56) continue;
        const double n = p[i] + q[i];
        if (!(n > 0)) continue;
        const double f = p[i] / n, e = std::sqrt(std::max(f * (1 - f), 1e-6) / n);
        xs.push_back(m - 0.545); ys.push_back(f); es.push_back(e);
        const double w = 1. / (e * e);
        S0 += w; Sx += w * xs.back(); Sxx += w * xs.back() * xs.back(); Sy += w * f; Sxy += w * xs.back() * f;
    }
    const double det = S0 * Sxx - Sx * Sx;
    if (det > 0 && xs.size() > 2) {
        const double a = (Sxx * Sy - Sx * Sxy) / det, b = (S0 * Sxy - Sx * Sy) / det;
        for (size_t i = 0; i < xs.size(); ++i) s.chi2 += std::pow((ys[i] - a - b * xs[i]) / es[i], 2);
        s.ndf = (int)xs.size() - 2;
        s.slope = a > 0 ? 0.1 * b / a : 0;
        double pa = 0, pb = 0;
        for (int i = ib(0.53); i < ib(0.56); ++i) { pa += p[i]; pb += q[i]; }
        const double n = pa + pb, f = n > 0 ? pa / n : 0;
        s.dPeak = f - a;
        s.eDPeak = n > 0 ? std::sqrt(std::max(f * (1 - f), 1e-6) / n + Sxx / det) : 0;
    }
    double mn = 1e9, mx = -1e9, sum = 0; int nn = 0;
    for (int i = ib(0.42); i < ib(0.56); ++i) {
        const double n = p[i] + q[i];
        if (!(n > 0)) continue;
        const double f = p[i] / n;
        mn = std::min(mn, f); mx = std::max(mx, f); sum += f; nn++;
    }
    s.range = (nn > 0 && sum > 0) ? (mx - mn) / (sum / nn) : 0;
    return s;
}

// ============================================================================
// plotting helpers
// ============================================================================
void Style() {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetPadLeftMargin(0.17);
    gStyle->SetPadRightMargin(0.04);
    gStyle->SetPadBottomMargin(0.14);
    gStyle->SetPadTopMargin(0.08);
    gStyle->SetLabelSize(0.05, "xy");
    gStyle->SetTitleSize(0.055, "xy");
    gStyle->SetTitleOffset(1.50, "y");
    gStyle->SetTitleOffset(1.10, "x");
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);
    gStyle->SetEndErrorSize(0);
}
void Deco(TH1* h, int col, int mk, double ms = 0.7, const char* xt = kXt, const char* yt = kYt) {
    h->SetLineColor(col); h->SetMarkerColor(col); h->SetMarkerStyle(mk); h->SetMarkerSize(ms);
    h->GetXaxis()->SetTitle(xt); h->GetYaxis()->SetTitle(yt);
    h->GetXaxis()->SetNdivisions(506);
}
void Head(const char* txt, double x = 0.17, double y = 0.94, double size = 0.042) {
    TLatex l; l.SetNDC(); l.SetTextSize(size); l.SetTextFont(42);
    l.DrawLatex(x, y, txt);
}
void Txt(double x, double y, const char* txt, double size = 0.04, int col = kBlack) {
    TLatex l; l.SetNDC(); l.SetTextSize(size); l.SetTextFont(42); l.SetTextColor(col);
    l.DrawLatex(x, y, txt);
}
void VLine(double x, double y1, double y2, int col, int sty, int wid = 2) {
    TLine* l = new TLine(x, y1, x, y2);
    l->SetLineColor(col); l->SetLineStyle(sty); l->SetLineWidth(wid);
    l->Draw();
}
void HLine(double x1, double x2, double y, int col, int sty, int wid = 1) {
    TLine* l = new TLine(x1, y, x2, y);
    l->SetLineColor(col); l->SetLineStyle(sty); l->SetLineWidth(wid);
    l->Draw();
}
TLegend* Leg(double x1, double y1, double x2, double y2, double size = 0.038) {
    TLegend* l = new TLegend(x1, y1, x2, y2);
    l->SetTextSize(size); l->SetTextFont(42); l->SetFillStyle(0);
    return l;
}
// largest bin content (plus error) inside [a, b]
double MaxIn(const TH1* h, double a, double b) {
    double m = 0;
    const TAxis* ax = h->GetXaxis();
    for (int i = ax->FindFixBin(a + 1e-9); i <= ax->FindFixBin(b - 1e-9); ++i) m = std::max(m, h->GetBinContent(i) + h->GetBinError(i));
    return m;
}
void Save(TCanvas* c, const char* name, bool pdf) {
    c->SaveAs(Form("%s/%s.png", scan::kFigDir, name));
    if (pdf) c->SaveAs(Form("%s/%s.pdf", scan::kFigDir, name));
}
TGraph* Curve(const PeakFit& f, bool total, int col, int sty) {
    TGraph* g = new TGraph();
    const int n = 400;
    for (int i = 0; i <= n; ++i) {
        const double m = f.cfg.lo + (f.cfg.hi - f.cfg.lo) * i / n;
        g->SetPoint(i, m, total ? f.Total(m) : f.Bkg(m));
    }
    g->SetLineColor(col); g->SetLineStyle(sty); g->SetLineWidth(2);
    return g;
}
TH1D* Resid(const Spectra& sp, const std::vector<double>& y, const std::vector<double>& var,
            const PeakFit& f, const char* name) {
    TH1D* h = new TH1D(name, "", sp.nbin, sp.lo, sp.hi);
    h->SetDirectory(nullptr);
    for (int i = 0; i < sp.nbin; ++i) {
        h->SetBinContent(i + 1, y[i] - f.Bkg(sp.x(i)));
        h->SetBinError(i + 1, std::sqrt(var[i]));
    }
    return h;
}

} // namespace sx

namespace sx {

std::vector<CR> gR;          // results of gC
Half gBase[3];               // no cut, per fold set

// ============================================================================
// stage 1: tables of all candidates
// ============================================================================
void StageTables() {
    const std::vector<std::pair<std::string, std::vector<int>>> fam = Families(gC);
    // working points of the scan families (all rows): peak by F, bump by Z
    std::vector<int> wpP(gC.size(), 0), wpB(gC.size(), 0);
    for (const auto& f : fam) {
        const int a = FamilyWP(gR, f.second, gBase[0], 0, 0);
        const int b = FamilyWP(gR, f.second, gBase[0], 0, 1);
        if (a >= 0) wpP[a] = 1;
        if (b >= 0) wpB[b] = 1;
    }

    // ------------------------------------------------------------------ peak
    OpenTable("extra_peak_candidates.txt");
    {
        const Half& b = gBase[0];
        P("#\n# PEAK. Start selection: %s. Fit P_main (frozen Gaussian %.4f / %.1f MeV + pol2 on N - CB, %.4f-%.4f).\n",
          kStartPeakTxt, scan::kPeakMean, 1e3 * scan::kPeakSigma, scan::kFitLo, scan::kFitHi);
        P("# F = S/sigma_S is statistical only (frozen peak shape): it compares cuts, it is NOT a significance. CB/N: CB = 2 sqrt(N++ N--) without a k factor, an estimate.\n");
        P("# The cuts were picked from this list on the full sample: with about 55 entries, pulls of G up to about 2 are what neutral cuts would give (see extra_ranking.txt\n");
        P("# and the two-fold checks in extra_combination.txt). Group e (crossing points): the resolution of the crossing point was never measured, so a failing dz_ee or\n");
        P("# rv_ee cut says nothing about converted photons against Dalitz pairs.\n");
        P("# start: N = %.0f rows in the peak window, CB = %.0f (share %.3f), S = %.0f +- %.0f, F = S/err = %.2f, S/B = %.4f, chi2/ndf = %.1f/%d\n",
          b.N, b.CB, b.CB / b.N, b.S, b.eS, b.F(), b.SoB, b.chi2, b.ndf);
        for (int f = 1; f < 3; ++f)
            P("# start, fold %d: N = %.0f, S = %.0f +- %.0f, F = %.2f, S/B = %.4f\n", f - 1, gBase[f].N, gBase[f].S, gBase[f].eS, gBase[f].F(), gBase[f].SoB);
        P("# epsS = S_pass/(S_pass+S_fail) ; epsB = kept part of (N - S in window) ; G = epsS/sqrt(epsB) (cut helps if G > 1) ; pull = (G-1)/err\n");
        P("# F, S/B, CB/N: after the cut ; clos = S_pass + S_fail - S_start ; G0, G1 = G on fold 0, fold 1 ; wp = working point of a scan (centre of the region F >= Fmax - 1), '-' = no cut\n#\n");
        P("# epsTag = efficiency of the same cut for the photon-tagged eta -> pi+pi-pi0 of the bump start selection (control sample, error 0.5-1.5 %%;\n");
        P("#          meaningful for cuts that do not act on the lepton pair or on photons), Gt = epsTag/sqrt(epsB) ; -1 = not applicable (photon cuts)\n#\n");
        P("%-2s %-36s | %6s %6s | %6s %6s | %6s %6s %5s | %6s %6s | %6s %6s | %6s %5s | %6s %5s | %5s | %9s | %6s %6s | %6s %6s %6s | %2s\n",
          "gr", "cut", "epsS", "err", "epsB", "sqrtB", "G", "err", "pull", "F", "dF", "S/B", "CB/N", "S_pass", "err",
          "S_fail", "err", "clos", "chi2/ndf", "G0", "G1", "epsTag", "err", "Gt", "wp");
        for (size_t i = 0; i < gC.size(); ++i) {
            const PeakNum p = PeakOf(gR[i], 0, gBase[0]);
            const PeakNum p0 = PeakOf(gR[i], 1, gBase[1]), p1 = PeakOf(gR[i], 2, gBase[2]);
            const Half& a = gR[i].h[0][0]; const Half& q = gR[i].h[0][1];
            const BumpNum bt = BumpOf(gR[i], 0, gBase[0]);
            const bool tagOk = gC[i].grp != "k";
            P("%-2s %-36s | %6.3f %6.3f | %6.3f %6.3f | %6.3f %6.3f %5.1f | %6.2f %+6.2f | %6.4f %6.3f | %6.0f %5.0f | %6.0f %5.0f | %5.0f | %5.1f/%-3d | %6.3f %6.3f | %6.3f %6.3f %6.3f | %2s\n",
              gC[i].grp.c_str(), gC[i].name.c_str(), p.epsS.e, p.epsS.err, p.epsB, std::sqrt(p.epsB), p.G, p.eG, p.pull,
              p.F, p.F - gBase[0].F(), p.SoB, p.cb, a.S, a.eS, q.S, q.eS, p.closure, a.chi2, a.ndf, p0.G, p1.G,
              tagOk ? bt.epsTag.e : -1., tagOk ? bt.epsTag.err : -1., (tagOk && p.epsB > 0) ? bt.epsTag.e / std::sqrt(p.epsB) : -1.,
              gC[i].fam.empty() ? "" : (wpP[i] ? "*" : ""));
        }
        P("#\n# working points of the scan families (peak, rule: centre of the grid region with F >= Fmax - 1; 'none' = no cut):\n");
        for (const auto& f : fam) {
            int nreg = 0;
            const int a = FamilyWP(gR, f.second, gBase[0], 0, 0, &nreg);
            if (a < 0) { P("#   %-12s none (region %d grid points incl. no cut)\n", f.first.c_str(), nreg); continue; }
            const PeakNum p = PeakOf(gR[a], 0, gBase[0]);
            P("#   %-12s %-36s region %d points | epsS %.3f +- %.3f epsB %.3f G %.3f +- %.3f F %.2f S/B %.4f\n",
              f.first.c_str(), gC[a].name.c_str(), nreg, p.epsS.e, p.epsS.err, p.epsB, p.G, p.eG, p.F, p.SoB);
        }
    }
    // ------------------------------------------------------------------ bump
    OpenTable("extra_bump_candidates.txt");
    {
        const Half& b = gBase[0];
        P("#\n# BUMP. Start selection: %s. Bump window %.4f-%.4f.\n", kStartBumpTxt, scan::kBumpLo, scan::kBumpHi);
        P("# One definition of the gain in all macros: Gb = (signal kept) / sqrt(rows kept in the bump window) -- the noise is sqrt(N), not sqrt(N - CB).\n");
        P("# CB = 2 sqrt(N++ N--) without a k factor (N--/N++ is about 6 in the bump window): CB/N and N - CB are estimates.\n");
        P("# start: N = %.0f rows in the bump window, CB = %.0f (share %.3f), N - CB = %.0f +- %.0f\n",
          b.bN, b.bCB, b.bCB / b.bN, b.bN - b.bCB, std::sqrt(b.bVarN + b.bVarCB));
        P("# start: edge steps h = %.0f +- %.0f, %.0f +- %.0f, %.0f +- %.0f rows / 2.5 MeV (M(pipi) %.2f-%.2f, %.2f-%.2f, %.2f-%.2f), sum %.0f +- %.0f, Z = %.1f\n",
          b.h[0], b.eh[0], b.h[1], b.eh[1], b.h[2], b.eh[2], kSlice[0], kSlice[1], kSlice[1], kSlice[2], kSlice[2], kSlice[3],
          b.hSum(), b.hErr(), b.Z());
        P("# start: tagged eta (pi0 photon %.3f < M(eeg) < %.3f; eta peak in the pi0-constrained M(pipieeg), frozen %.4f / %.1f MeV, pol%d, %.2f-%.2f, Poisson fit):\n",
          kPi0Lo, kPi0Hi, kM5Mean, 1e3 * kM5Sigma, kM5Order, kM5Lo, kM5Hi);
        P("#        S5 = %.0f +- %.0f, rows in +-1.9 sigma = %.0f, CB = %.0f, S/B = %.3f, chi2/ndf = %.1f/%d\n",
          b.S5, b.eS5, b.N5, b.CB5, b.SoB5, b.chi5, b.ndf5);
        P("# epsTag = S5_pass/(S5_pass+S5_fail) ; epsN = rows kept in the bump window ; epsNCB = (N - CB) kept ; Gb = epsTag/sqrt(epsN) ; pull = (Gb-1)/err\n");
        P("# epsH = h_pass/(h_pass+h_fail), h = summed step height of the three slices ; Z = h/err after the cut ; F5 = S5/err after the cut\n");
        P("# Gb0, Gb1 = Gb on fold 0, 1 ; wp = working point of a scan (centre of the region Z >= Zmax - 1)\n#\n");
        P("%-2s %-36s | %6s %6s | %6s %6s %6s | %6s %6s %5s | %6s %6s | %6s %5s %6s %6s | %6s %6s | %6s %6s | %2s\n",
          "gr", "cut", "epsTag", "err", "epsN", "sqrtN", "epsNCB", "Gb", "err", "pull", "CB/N", "S/B5", "h_pass", "err",
          "epsH", "err", "Z", "dZ", "Gb0", "Gb1", "wp");
        for (size_t i = 0; i < gC.size(); ++i) {
            const BumpNum p = BumpOf(gR[i], 0, gBase[0]);
            const BumpNum p0 = BumpOf(gR[i], 1, gBase[1]), p1 = BumpOf(gR[i], 2, gBase[2]);
            if (gC[i].grp == "k") {   // the tag itself needs a photon: no tagged efficiency for photon cuts
                P("%-2s %-36s | %6s %6s | %6.3f %6.3f %6.3f | %6s %6s %5s | %6.3f %6s | %6.0f %5.0f %6.3f %6.3f | %6.2f %+6.2f | %6s %6s | %2s\n",
                  gC[i].grp.c_str(), gC[i].name.c_str(), "n/a", "", p.epsN, std::sqrt(p.epsN), p.epsNCB, "n/a", "", "", p.cb, "",
                  p.hP, p.ehP, p.epsH.e, p.epsH.err, p.Z, p.Z - gBase[0].Z(), "", "", "");
                continue;
            }
            P("%-2s %-36s | %6.3f %6.3f | %6.3f %6.3f %6.3f | %6.3f %6.3f %5.1f | %6.3f %6.3f | %6.0f %5.0f %6.3f %6.3f | %6.2f %+6.2f | %6.3f %6.3f | %2s\n",
              gC[i].grp.c_str(), gC[i].name.c_str(), p.epsTag.e, p.epsTag.err, p.epsN, std::sqrt(p.epsN), p.epsNCB,
              p.Gb, p.eGb, p.pull, p.cb, p.SoB5, p.hP, p.ehP, p.epsH.e, p.epsH.err, p.Z, p.Z - gBase[0].Z(), p0.Gb, p1.Gb,
              gC[i].fam.empty() ? "" : (wpB[i] ? "*" : ""));
        }
        P("#\n# working points of the scan families (bump, rule: centre of the grid region with Z >= Zmax - 1; 'none' = no cut):\n");
        for (const auto& f : fam) {
            int nreg = 0;
            const int a = FamilyWP(gR, f.second, gBase[0], 0, 1, &nreg);
            if (a < 0) { P("#   %-12s none (region %d grid points incl. no cut)\n", f.first.c_str(), nreg); continue; }
            const BumpNum p = BumpOf(gR[a], 0, gBase[0]);
            P("#   %-12s %-36s region %d points | epsTag %.3f +- %.3f epsN %.3f Gb %.3f +- %.3f epsH %.3f +- %.3f Z %.2f\n",
              f.first.c_str(), gC[a].name.c_str(), nreg, p.epsTag.e, p.epsTag.err, p.epsN, p.Gb, p.eGb, p.epsH.e, p.epsH.err, p.Z);
        }
    }
    // ------------------------------------------------------------------ sculpting
    OpenTable("extra_sculpt.txt");
    {
        P("#\n# Pass fraction of unlike-sign rows against M4 (10 MeV bins). Regions [GeV/c2]:\n#  ");
        for (int k = 0; k < 8; ++k) P(" r%d = %.2f-%.2f", k, kReg[k], kReg[k + 1]);
        P("\n# slope = relative change of the pass fraction per 100 MeV (line fitted on 0.52-0.62 without 0.53-0.56, i.e. the side regions of the peak fit), chi2/ndf of that line,\n");
        P("# dPk = pass fraction in 0.53-0.56 minus the line [%%] with error, range = (max - min)/mean over 0.42-0.56.\n");
        P("# A cut sculpts when the chi2 of the line is bad or dPk is far from (epsS - epsB) * S/N (about 0.06 * (epsS - epsB)).\n#\n");
        for (int w = 0; w < 2; ++w) {
            P("# ---- %s start selection ----\n", w == 0 ? "PEAK" : "BUMP");
            P("%-2s %-36s | %5s %5s %5s %5s %5s %5s %5s %5s | %7s %9s | %6s %5s | %5s\n", "gr", "cut", "r0", "r1", "r2", "r3",
              "r4", "r5", "r6", "r7", "slope", "chi2/ndf", "dPk[%]", "err", "range");
            for (size_t i = 1; i < gC.size(); ++i) {
                const Sculpt s = w == 0 ? SculptOf(gR[i].pP, gR[i].pF) : SculptOf(gR[i].bP, gR[i].bF);
                P("%-2s %-36s | %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f | %+7.3f %5.1f/%-3d | %+6.2f %5.2f | %5.3f\n",
                  gC[i].grp.c_str(), gC[i].name.c_str(), s.reg[0], s.reg[1], s.reg[2], s.reg[3], s.reg[4], s.reg[5], s.reg[6],
                  s.reg[7], s.slope, s.chi2, s.ndf, 100 * s.dPeak, 100 * s.eDPeak, s.range);
            }
        }
    }
    CloseTable();
}


// ============================================================================
// stage N: sideband-subtracted comparison of all available variables (peak)
// ============================================================================
// For every variable the distribution of the unlike-sign rows in the peak window
// is compared with the mean shape of two sidebands of 11.1 MeV on both sides
// (lower 0.5225-0.5336, upper 0.5542-0.5653). Bins: deciles of the sideband
// distribution, or the values themselves for variables with at most 12 values.
//   chi2  = sum (n_i - N p_i)^2 / var_i, p = sideband shape; tests "signal and
//           background have the same distribution"; z = (chi2 - ndf)/sqrt(2 ndf)
//   The same test in two CONTROL windows without a signal (centres 0.5760 and
//   0.6080) shows what a variable does only because it is tied to M4.
//   signal shape s_i = (n_i - B p_i) / S_win (B = N - S_win, S_win from the fit)
//   best one-sided cut on a bin edge: eps_S, eps_B, G = eps_S / sqrt(eps_B)
struct VarDef {
    std::string name;
    std::function<double(R)> f;
    int tied;      // 1 = tied to M4 by kinematics or a parameter of the main scan
};
inline void Dir(double th, double ph, double d[3]) {
    const double t = th * TMath::DegToRad(), p = ph * TMath::DegToRad();
    d[0] = std::sin(t) * std::cos(p); d[1] = std::sin(t) * std::sin(p); d[2] = std::cos(t);
}
inline double AngleDeg(double t1, double p1, double t2, double p2) {
    double a[3], b[3]; Dir(t1, p1, a); Dir(t2, p2, b);
    double c = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    c = std::max(-1., std::min(1., c));
    return std::acos(c) * TMath::RadToDeg();
}
// momentum vector of a pair [MeV/c]
inline void PairP(double p1, double t1, double f1, double p2, double t2, double f2, double v[3]) {
    double a[3], b[3]; Dir(t1, f1, a); Dir(t2, f2, b);
    for (int i = 0; i < 3; ++i) v[i] = p1 * a[i] + p2 * b[i];
}
inline double DPhi(double a, double b) {       // |a - b| folded to 0..180
    double d = std::fabs(a - b);
    while (d > 360.) d -= 360.;
    return d > 180. ? 360. - d : d;
}
inline double BetaPi(double p) { return p / std::sqrt(p * p + 139.57 * 139.57); }

std::vector<VarDef> BuildVars() {
    std::vector<VarDef> v;
    auto add = [&](const char* n, std::function<double(R)> f, int tied = 0) {
        VarDef d; d.name = n; d.f = f; d.tied = tied; v.push_back(d);
    };
    // event
    add("vz", [](R r) { return (double)r.vz; });
    add("rT (vertex)", [](R r) { return RT(r); });
    add("vchi2", [](R r) { return (double)r.vchi2; });
    add("evtStart", [](R r) { return (double)r.evtStart; });
    add("trigbit", [](R r) { return (double)E(r).trigbit; });
    add("start_iteration", [](R r) { return (double)E(r).start_iteration; });
    add("hyp_chi2", [](R r) { return (double)r.hyp_chi2; });
    add("isBest", [](R r) { return (double)r.isBest; });
    add("nrow", [](R r) { return (double)std::min(r.nrow, 4); });
    add("n_pip", [](R r) { return (double)std::min(r.n_pip, 3); });
    add("n_pim", [](R r) { return (double)std::min(r.n_pim, 3); });
    add("n_lep", [](R r) { return (double)std::min(r.n_lep, 4); });
    add("n_gam", [](R r) { return (double)std::min(r.n_gam, 3); });
    add("neutr_mult", [](R r) { return (double)std::min(r.neutr_mult, 4); });
    add("fwdet_mult", [](R r) { return (double)std::min(r.fwdet_mult, 2); });
    add("fw1_theta", [](R r) { return (double)r.fw1_theta; });
    add("fw1_p", [](R r) { return (double)r.fw1_p; });
    add("g1_e", [](R r) { return (double)E(r).g1_e; });
    add("M(eeg) best photon", [](R r) { return (double)r.m_eeg; });
    // pair kinematics
    add("m_ee", [](R r) { return (double)r.mee; }, 1);
    add("oa (e+e-)", [](R r) { return (double)r.oa; }, 1);
    add("M(pipi)", [](R r) { return (double)r.mpipi; }, 1);
    add("oa_lab", [](R r) { return (double)r.oa_lab; }, 1);
    add("MM4", [](R r) { return (double)r.mm4; }, 1);
    add("MM(pipi)", [](R r) { return (double)r.mm_pipi; }, 1);
    add("MM(ee)", [](R r) { return (double)r.mm_ee; }, 1);
    add("pt4", [](R r) { return (double)r.pt4; }, 1);
    add("y4", [](R r) { return (double)r.y4; });
    add("oa (pi+pi-)", [](R r) { return AngleDeg(r.pip_theta, r.pip_phi, r.pim_theta, r.pim_phi); }, 1);
    add("p(e+e-)", [](R r) { double q[3]; PairP(r.ep_p, r.ep_theta, r.ep_phi, r.em_p, r.em_theta, r.em_phi, q);
                             return std::sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]); }, 1);
    add("p(pi+pi-)", [](R r) { double q[3]; PairP(r.pip_p, r.pip_theta, r.pip_phi, r.pim_p, r.pim_theta, r.pim_phi, q);
                               return std::sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]); }, 1);
    add("theta(e+e-)", [](R r) { double q[3]; PairP(r.ep_p, r.ep_theta, r.ep_phi, r.em_p, r.em_theta, r.em_phi, q);
                                 return std::atan2(std::hypot(q[0], q[1]), q[2]) * TMath::RadToDeg(); });
    add("theta(pi+pi-)", [](R r) { double q[3]; PairP(r.pip_p, r.pip_theta, r.pip_phi, r.pim_p, r.pim_theta, r.pim_phi, q);
                                   return std::atan2(std::hypot(q[0], q[1]), q[2]) * TMath::RadToDeg(); });
    add("dphi(pairs)", [](R r) { double a[3], b[3];
                                 PairP(r.ep_p, r.ep_theta, r.ep_phi, r.em_p, r.em_theta, r.em_phi, a);
                                 PairP(r.pip_p, r.pip_theta, r.pip_phi, r.pim_p, r.pim_theta, r.pim_phi, b);
                                 return DPhi(std::atan2(a[1], a[0]) * TMath::RadToDeg(), std::atan2(b[1], b[0]) * TMath::RadToDeg()); }, 1);
    add("p asymmetry leptons", [](R r) { return std::fabs(r.ep_p - r.em_p) / (double)(r.ep_p + r.em_p); });
    add("p asymmetry pions", [](R r) { return std::fabs(r.pip_p - r.pim_p) / (double)(r.pip_p + r.pim_p); });
    // leptons (symmetric)
    add("p lepton min", [](R r) { return (double)Mn(r.ep_p, r.em_p); });
    add("p lepton max", [](R r) { return (double)Mx(r.ep_p, r.em_p); }, 1);
    add("theta lepton min", [](R r) { return (double)Mn(r.ep_theta, r.em_theta); });
    add("theta lepton max", [](R r) { return (double)Mx(r.ep_theta, r.em_theta); });
    add("rkchi2 lepton max", [](R r) { return (double)Mx(r.ep_rkchi2, r.em_rkchi2); });
    add("richmatchq max", [](R r) { return (double)Mx(r.ep_richmatchquality, r.em_richmatchquality); });
    add("ringchi2 max", [](R r) { return (double)Mx(r.ep_ringchi2, r.em_ringchi2); });
    add("ringchi2 min", [](R r) { return (double)Mn(r.ep_ringchi2, r.em_ringchi2); });
    add("kIsLepton sum", [](R r) { return (double)(r.ep_kIsLepton + r.em_kIsLepton); });
    add("isUsedVertex leptons", [](R r) { return (double)(r.ep_isUsedVertex + r.em_isUsedVertex); });
    add("shared ring", [](R r) { return (double)r.shared_ring; });
    add("lepton system class", [](R r) { return (double)LepSys(r); });
    add("oa_lept min |.|", [](R r) { return (double)Mn(std::fabs(r.ep_oa_lept), std::fabs(r.em_oa_lept)); });
    add("oa_lept max |.|", [](R r) { return (double)Mx(std::fabs(r.ep_oa_lept), std::fabs(r.em_oa_lept)); });
    add("oa_hadr min", [](R r) { return (double)Mn(r.ep_oa_hadr, r.em_oa_hadr); });
    add("oa_hadr max", [](R r) { return (double)Mx(r.ep_oa_hadr, r.em_oa_hadr); });
    add("|beta-1| lepton max", [](R r) { return (double)Mx(std::fabs(E(r).ep_beta_new - 1.f), std::fabs(E(r).em_beta_new - 1.f)); });
    add("|tof_new-tof_mom| lepton max", [](R r) { return (double)Mx(std::fabs(E(r).ep_tof_new - E(r).ep_tof_mom),
                                                                   std::fabs(E(r).em_tof_new - E(r).em_tof_mom)); });
    add("|theta_rich-theta| max", [](R r) { return (double)Mx(std::fabs(E(r).ep_theta_rich - r.ep_theta),
                                                             std::fabs(E(r).em_theta_rich - r.em_theta)); });
    add("|phi_rich-phi| max", [](R r) { return std::max(DPhi(E(r).ep_phi_rich, r.ep_phi), DPhi(E(r).em_phi_rich, r.em_phi)); });
    add("rich_padnum min", [](R r) { return (double)Mn(E(r).ep_rich_padnum, E(r).em_rich_padnum); });
    add("rich_padnum max", [](R r) { return (double)Mx(E(r).ep_rich_padnum, E(r).em_rich_padnum); });
    add("leptons same sector", [](R r) { return E(r).ep_sector == E(r).em_sector ? 1. : 0.; });
    add("|r| lepton max", [](R r) { return (double)Mx(std::fabs(E(r).ep_r), std::fabs(E(r).em_r)); });
    add("|z-vz| lepton max", [](R r) { return (double)Mx(std::fabs(E(r).ep_z - r.vz), std::fabs(E(r).em_z - r.vz)); });
    // pions
    add("p pion min", [](R r) { return (double)Mn(r.pip_p, r.pim_p); }, 1);
    add("p pion max", [](R r) { return (double)Mx(r.pip_p, r.pim_p); }, 1);
    add("theta pion min", [](R r) { return (double)Mn(r.pip_theta, r.pim_theta); });
    add("theta pion max", [](R r) { return (double)Mx(r.pip_theta, r.pim_theta); });
    add("rkchi2 pion max", [](R r) { return (double)Mx(r.pip_rkchi2, r.pim_rkchi2); });
    add("|mass2-m_pi^2| max", [](R r) { return (double)Mx(std::fabs(r.pip_mass2 - kPionM2), std::fabs(r.pim_mass2 - kPionM2)); });
    add("dE/dx MDC pion max", [](R r) { return (double)Mx(r.pip_dedx_mdc, r.pim_dedx_mdc); });
    add("dE/dx MDC pion min", [](R r) { return (double)Mn(r.pip_dedx_mdc, r.pim_dedx_mdc); });
    add("pions with ring", [](R r) { return (double)((r.pip_isring != 0) + (r.pim_isring != 0)); });
    add("isUsedVertex pions", [](R r) { return (double)(r.pip_isUsedVertex + r.pim_isUsedVertex); });
    add("pip_system", [](R r) { return (double)r.pip_system; });
    add("pim_system", [](R r) { return (double)r.pim_system; });
    add("|beta-beta(p)| pion max", [](R r) { return std::max(std::fabs(E(r).pip_beta_new - BetaPi(r.pip_p_rec)),
                                                            std::fabs(E(r).pim_beta_new - BetaPi(r.pim_p_rec))); });
    add("|tof_new-tof_mom| pion max", [](R r) { return (double)Mx(std::fabs(E(r).pip_tof_new - E(r).pip_tof_mom),
                                                                 std::fabs(E(r).pim_tof_new - E(r).pim_tof_mom)); });
    add("pions same sector", [](R r) { return E(r).pip_sector == E(r).pim_sector ? 1. : 0.; });
    add("|r| pion max", [](R r) { return (double)Mx(std::fabs(E(r).pip_r), std::fabs(E(r).pim_r)); });
    add("|z-vz| pion max", [](R r) { return (double)Mx(std::fabs(E(r).pip_z - r.vz), std::fabs(E(r).pim_z - r.vz)); });
    add("sectors hit by the four", [](R r) {
        const float q[4] = {E(r).ep_sector, E(r).em_sector, E(r).pip_sector, E(r).pim_sector};
        int n = 0;
        for (int i = 0; i < 4; ++i) { bool nw = true; for (int j = 0; j < i; ++j) if (q[j] == q[i]) nw = false; n += nw; }
        return (double)n; });
    // crossing points
    add("dz_ee", [](R r) { return (double)r.dz_ee; });
    add("rv_ee", [](R r) { return (double)r.rv_ee; });
    add("doca_ee", [](R r) { return (double)r.doca_ee; });
    add("dz_pipi", [](R r) { return (double)r.dz_pipi; });
    add("|dz_pipi|", [](R r) { return (double)std::fabs(r.dz_pipi); });
    add("rv_pipi", [](R r) { return (double)r.rv_pipi; });
    add("doca_pipi", [](R r) { return (double)r.doca_pipi; });
    return v;
}

struct SbRes {
    int iv = 0;
    int nb = 0, ndf = 0; bool discrete = false;
    double chi2[3] = {0, 0, 0}, z[3] = {0, 0, 0};
    double sep = 0;
    double bG = 0, eG = 0, bEpsS = 0, eEpsS = 0, bEpsB = 0, bEdge = 0; int bDir = 0;   // dir +1: x > edge, -1: x < edge, 0: veto value
    std::vector<double> edges;     // inner bin edges
};
const double kSbCen[3] = {scan::kPeakMean, 0.5760, 0.6080};
const double kSbHalf = scan::kWinNsig * scan::kPeakSigma;    // 10.26 MeV
const double kSbWid = scan::kPeakLo - scan::kFitLo;          // 11.14 MeV
// region of a row for window w: 1 = window, 2 = lower sideband, 3 = upper sideband, 0 = none
inline int SbRegion(double m, int w) {
    const double d = m - kSbCen[w];
    if (d >= -kSbHalf && d < kSbHalf) return 1;
    if (d >= -kSbHalf - kSbWid && d < -kSbHalf) return 2;
    if (d >= kSbHalf && d < kSbHalf + kSbWid) return 3;
    return 0;
}

std::vector<VarDef> gV;
std::vector<SbRes> gSb;

SbRes SbAnalyse(int iv, double Swin) {
    const VarDef& vd = gV[iv];
    SbRes res; res.iv = iv;
    std::vector<double> val[3][4];
    for (const auto& r : gP[0].s[0]) {
        if (!(r.m4 >= scan::kFitLo && r.m4 < 0.632)) continue;
        double x = 0; bool have = false;
        for (int w = 0; w < 3; ++w) {
            const int reg = SbRegion(r.m4, w);
            if (!reg) continue;
            if (!have) { x = vd.f(r); have = true; }
            val[w][reg].push_back(x);
        }
    }
    // bin edges from the sidebands of the peak window
    std::vector<double> sb = val[0][2];
    sb.insert(sb.end(), val[0][3].begin(), val[0][3].end());
    std::sort(sb.begin(), sb.end());
    std::vector<double> dist;
    for (size_t i = 0; i < sb.size() && dist.size() <= 12; ++i) if (dist.empty() || sb[i] != dist.back()) dist.push_back(sb[i]);
    std::vector<double>& ed = res.edges;
    if (dist.size() <= 12) {
        res.discrete = true;
        for (size_t i = 0; i + 1 < dist.size(); ++i) ed.push_back(0.5 * (dist[i] + dist[i + 1]));
    } else {
        for (int q = 1; q < 10; ++q) {
            const double e = sb[(size_t)(q * 0.1 * sb.size())];
            if (ed.empty() || e > ed.back()) ed.push_back(e);
        }
        // an edge equal to the smallest value would leave the first bin empty
        if (!ed.empty() && ed.front() <= sb.front()) ed.erase(ed.begin());
    }
    const int nb = (int)ed.size() + 1;
    res.nb = nb; res.ndf = nb - 1;
    if (nb < 2) return res;
    auto bin = [&](double x) { return (int)(std::upper_bound(ed.begin(), ed.end(), x) - ed.begin()); };
    for (int w = 0; w < 3; ++w) {
        std::vector<double> n(nb, 0.), a(nb, 0.), b(nb, 0.);
        for (double x : val[w][1]) n[bin(x)] += 1;
        for (double x : val[w][2]) a[bin(x)] += 1;
        for (double x : val[w][3]) b[bin(x)] += 1;
        const double N = val[w][1].size(), A = val[w][2].size(), B = val[w][3].size();
        if (!(N > 0 && A > 0 && B > 0)) continue;
        std::vector<double> p(nb), vp(nb);
        double c2 = 0;
        for (int i = 0; i < nb; ++i) {
            p[i] = 0.5 * (a[i] / A + b[i] / B);
            vp[i] = 0.25 * (a[i] / (A * A) + b[i] / (B * B));
            const double var = N * p[i] + N * N * vp[i];
            if (var > 0) c2 += (n[i] - N * p[i]) * (n[i] - N * p[i]) / var;
        }
        res.chi2[w] = c2;
        res.z[w] = (c2 - res.ndf) / std::sqrt(2. * res.ndf);
        if (w != 0) continue;
        // signal shape, separation, best one-sided cut
        const double Bk = N - Swin;
        std::vector<double> sg(nb);
        for (int i = 0; i < nb; ++i) {
            sg[i] = (n[i] - Bk * p[i]) / Swin;
            const double den = std::fabs(sg[i]) + p[i];
            if (den > 0) res.sep += 0.5 * (sg[i] - p[i]) * (sg[i] - p[i]) / den;
        }
        auto test = [&](const std::vector<int>& keep, double edge, int dir) {
            double nK = 0, pK = 0, aK = 0, bK = 0;
            for (int i : keep) { nK += n[i]; pK += p[i]; aK += a[i]; bK += b[i]; }
            if (pK < 0.3 || pK > 0.995) return;
            const double q = nK / N, pa = aK / A, pb = bK / B;
            const double eS = (nK - Bk * pK) / Swin;
            const double vS = (N * q * (1 - q) + Bk * Bk * 0.25 * (pa * (1 - pa) / A + pb * (1 - pb) / B)) / (Swin * Swin);
            const double G = eS / std::sqrt(pK);
            if (G > res.bG) {
                res.bG = G; res.eG = std::sqrt(vS) / std::sqrt(pK); res.bEpsS = eS; res.eEpsS = std::sqrt(vS);
                res.bEpsB = pK; res.bEdge = edge; res.bDir = dir;
            }
        };
        for (int j = 1; j < nb; ++j) {
            std::vector<int> lo, hi;
            for (int i = 0; i < j; ++i) lo.push_back(i);
            for (int i = j; i < nb; ++i) hi.push_back(i);
            test(lo, ed[j - 1], -1);
            test(hi, ed[j - 1], +1);
        }
        if (res.discrete && nb > 2)
            for (int j = 0; j < nb; ++j) {
                std::vector<int> k;
                for (int i = 0; i < nb; ++i) if (i != j) k.push_back(i);
                test(k, dist[j], 0);
            }
    }
    return res;
}

// sideband-subtracted (signal) and sideband (background) distributions, unit area
void SbHists(int iv, double Swin, int nbin, double lo, double hi, TH1D*& hs, TH1D*& hb) {
    const VarDef& vd = gV[iv];
    hs = new TH1D(Form("sbS_%d", iv), "", nbin, lo, hi); hs->SetDirectory(nullptr);
    hb = new TH1D(Form("sbB_%d", iv), "", nbin, lo, hi); hb->SetDirectory(nullptr);
    std::vector<double> n(nbin, 0.), a(nbin, 0.), b(nbin, 0.);
    double N = 0, A = 0, B = 0;
    for (const auto& r : gP[0].s[0]) {
        const int reg = SbRegion(r.m4, 0);
        if (!reg) continue;
        const double x = std::max(lo, std::min(hi - 1e-9 * (hi - lo), vd.f(r)));   // under/overflow into the end bins
        const int i = std::min(nbin - 1, (int)((x - lo) / (hi - lo) * nbin));
        if (reg == 1) { n[i] += 1; N += 1; } else if (reg == 2) { a[i] += 1; A += 1; } else { b[i] += 1; B += 1; }
    }
    const double Bk = N - Swin;
    for (int i = 0; i < nbin; ++i) {
        const double p = 0.5 * (a[i] / A + b[i] / B), vp = 0.25 * (a[i] / (A * A) + b[i] / (B * B));
        hb->SetBinContent(i + 1, p); hb->SetBinError(i + 1, std::sqrt(vp));
        hs->SetBinContent(i + 1, (n[i] - Bk * p) / Swin);
        hs->SetBinError(i + 1, std::sqrt(n[i] + Bk * Bk * vp) / Swin);
    }
}

std::vector<int> gTop3;      // indices into gSb of the three variables with the largest separation
void StageN() {
    TStopwatch sw;
    gV = BuildVars();
    const double Swin = gBase[0].Swin;
    gSb.clear();
    for (size_t i = 0; i < gV.size(); ++i) gSb.push_back(SbAnalyse((int)i, Swin));
    std::vector<int> ord(gSb.size());
    for (size_t i = 0; i < ord.size(); ++i) ord[i] = (int)i;
    std::sort(ord.begin(), ord.end(), [&](int a, int b) { return gSb[a].z[0] > gSb[b].z[0]; });
    OpenTable("extra_sbsub.txt");
    P("#\n# N. Peak rows against sideband rows (start selection %s; unlike-sign rows).\n", kStartPeakTxt);
    P("# peak window %.5f-%.5f, sidebands %.4f-%.5f and %.5f-%.5f ; S in the window (fit) = %.0f of %.0f rows.\n",
      scan::kPeakLo, scan::kPeakHi, scan::kFitLo, scan::kPeakLo, scan::kPeakHi, scan::kPeakHi + kSbWid, Swin, gBase[0].N);
    P("# chi2/ndf and z = (chi2 - ndf)/sqrt(2 ndf): window against the mean sideband shape (deciles of the sidebands, or the values).\n");
    P("# z_c1, z_c2: the same in control windows centred at %.4f and %.4f (no signal): a large value means the variable is tied to M4.\n", kSbCen[1], kSbCen[2]);
    P("# sep = 1/2 sum (s - b)^2 / (|s| + b) of the sideband-subtracted signal shape s (inflated by noise: read z, not sep).\n");
    P("# best one-sided cut on a bin edge: epsS, epsB, G = epsS/sqrt(epsB) (maximum of up to 18 tries: biased upward).\n");
    P("# tied = 1: kinematically tied to M4 or a parameter of the main scan (not ranked).\n#\n");
    P("%-30s %4s %3s | %9s %6s | %6s %6s | %6s | %-22s %6s %6s %6s %6s %6s\n", "variable", "tied", "nb", "chi2/ndf", "z",
      "z_c1", "z_c2", "sep", "best cut", "epsS", "err", "epsB", "G", "err");
    for (int i : ord) {
        const SbRes& s = gSb[i];
        const std::string cut = s.bDir == 0 ? Form("veto value %.4g", s.bEdge) : Form("%s %.5g", s.bDir > 0 ? ">" : "<", s.bEdge);
        P("%-30s %4d %3d | %6.1f/%-2d %6.1f | %6.1f %6.1f | %6.4f | %-22s %6.3f %6.3f %6.3f %6.3f %6.3f\n",
          gV[i].name.c_str(), gV[i].tied, s.nb, s.chi2[0], s.ndf, s.z[0], s.z[1], s.z[2], s.sep, cut.c_str(), s.bEpsS, s.eEpsS,
          s.bEpsB, s.bG, s.eG);
    }
    // the three variables with the largest separation: not tied, controls quiet (|z| < 3)
    gTop3.clear();
    for (int i : ord) {
        if (gV[i].tied) continue;
        if (std::fabs(gSb[i].z[1]) > 3. || std::fabs(gSb[i].z[2]) > 3.) continue;
        gTop3.push_back(i);
        if (gTop3.size() == 3) break;
    }
    P("#\n# three variables with the largest separation (not tied to M4, both control windows with |z| < 3):\n");
    for (int i : gTop3)
        P("#   %-30s z = %.1f (controls %.1f, %.1f), best one-sided cut %s %.5g: epsS %.3f +- %.3f, epsB %.3f, G %.3f +- %.3f\n",
          gV[i].name.c_str(), gSb[i].z[0], gSb[i].z[1], gSb[i].z[2], gSb[i].bDir > 0 ? ">" : (gSb[i].bDir < 0 ? "<" : "veto"),
          gSb[i].bEdge, gSb[i].bEpsS, gSb[i].eEpsS, gSb[i].bEpsB, gSb[i].bG, gSb[i].eG);
    CloseTable();
    printf("StageN: %zu variables in %.1f s\n", gV.size(), sw.RealTime());
}


// ============================================================================
// stage R: ranking of the candidates, separately for the peak and for the bump
// ============================================================================
// One entry per single cut and one per scan family (its working point).
struct Entry { int idx; std::string fam; };
std::vector<Entry> Entries(int f, int kind) {
    std::vector<Entry> e;
    for (size_t i = 1; i < gC.size(); ++i) if (gC[i].fam.empty()) { Entry x; x.idx = (int)i; e.push_back(x); }
    for (const auto& fm : Families(gC)) {
        const int a = FamilyWP(gR, fm.second, gBase[f], f, kind);
        if (a >= 0) { Entry x; x.idx = a; x.fam = fm.first; e.push_back(x); }
    }
    return e;
}
// sculpting flags of a candidate for the peak: T = pass fraction tilted by more than 8 % per 100 MeV across
// the fit range, K = pass fraction in the peak region off the side-region line by more than 3 errors after
// allowing for the signal, C = pass and fail yields do not add up (|closure| > 0.15 sigma_S)
std::string PeakFlags(int i) {
    const Sculpt sc = SculptOf(gR[i].pP, gR[i].pF);
    const PeakNum p = PeakOf(gR[i], 0, gBase[0]);
    std::string fl;
    if (std::fabs(sc.slope) > 0.08) fl += "T";
    const double expect = (p.epsS.e - p.epsB) * gBase[0].Swin / gBase[0].N * (20.5 / 30.);
    if (sc.eDPeak > 0 && std::fabs(sc.dPeak - expect) > 3. * sc.eDPeak) fl += "K";
    if (std::fabs(p.closure) > 0.15 * gBase[0].eS) fl += "C";
    return fl.empty() ? "-" : fl;
}
void StageRank() {
    OpenTable("extra_ranking.txt");
    // ---- peak ----
    {
        std::vector<Entry> e = Entries(0, 0);
        std::sort(e.begin(), e.end(), [&](const Entry& a, const Entry& b) {
            return PeakOf(gR[a.idx], 0, gBase[0]).G > PeakOf(gR[b.idx], 0, gBase[0]).G; });
        P("#\n# PEAK ranking by G = epsS/sqrt(epsB) (pass/fail fits; scan families enter with their working point, rule F >= Fmax - 1).\n");
        P("# flags: T = pass fraction tilted > 8 %% per 100 MeV in the fit range, K = kink in the peak region, C = pass + fail does not close, - = none.\n");
        P("# epsTag, Gt: the same cut on the photon-tagged eta -> pi+pi-pi0 (control sample; -1 = not applicable, photon cuts).\n");
        P("# CAUTION: the entries were ranked on the full sample. With %zu entries the largest pulls of G (about 2) are what neutral cuts would give:\n", e.size());
        P("# the ranking alone establishes nothing. A cut counts as established only where the control sample confirms it (Gt = epsTag/sqrt(epsB) > 1\n");
        P("# with the small error of epsTag): the vertex cuts. See the two-fold checks in extra_combination.txt.\n#\n");
        P("%-3s %-2s %-36s | %6s %6s %6s | %6s %6s %5s | %6s %6s %6s | %6s %6s | %6s %6s | %s\n", "#", "gr", "cut", "epsS", "err",
          "epsB", "G", "err", "pull", "dF", "S/B", "CB/N", "G0", "G1", "epsTag", "Gt", "flags");
        int n = 0;
        for (const Entry& x : e) {
            const PeakNum p = PeakOf(gR[x.idx], 0, gBase[0]);
            const PeakNum p0 = PeakOf(gR[x.idx], 1, gBase[1]), p1 = PeakOf(gR[x.idx], 2, gBase[2]);
            const BumpNum bt = BumpOf(gR[x.idx], 0, gBase[0]);
            const bool tagOk = gC[x.idx].grp != "k";
            P("%-3d %-2s %-36s | %6.3f %6.3f %6.3f | %6.3f %6.3f %5.1f | %+6.2f %6.4f %6.3f | %6.3f %6.3f | %6.3f %6.3f | %s\n", ++n,
              gC[x.idx].grp.c_str(), gC[x.idx].name.c_str(), p.epsS.e, p.epsS.err, p.epsB, p.G, p.eG, p.pull, p.F - gBase[0].F(),
              p.SoB, p.cb, p0.G, p1.G, tagOk ? bt.epsTag.e : -1., tagOk ? bt.epsTag.e / std::sqrt(p.epsB) : -1., PeakFlags(x.idx).c_str());
        }
    }
    // ---- bump ----
    {
        std::vector<Entry> e = Entries(0, 1);
        std::vector<Entry> e2;
        for (const Entry& x : e) if (gC[x.idx].grp != "k") e2.push_back(x);
        std::sort(e2.begin(), e2.end(), [&](const Entry& a, const Entry& b) {
            return BumpOf(gR[a.idx], 0, gBase[0]).Gb > BumpOf(gR[b.idx], 0, gBase[0]).Gb; });
        P("#\n# BUMP ranking by Gb = epsTag/sqrt(epsN) (tagged-eta efficiency over the root of the rows kept in the bump window;\n");
        P("# scan families enter with their working point, rule Z >= Zmax - 1, Z = summed edge-step height / error). Photon cuts: see below.\n#\n");
        P("%-3s %-2s %-36s | %6s %6s %6s %6s | %6s %6s %5s | %6s %6s | %6s %6s | %6s %6s\n", "#", "gr", "cut", "epsTag", "err",
          "epsN", "epsNCB", "Gb", "err", "pull", "epsH", "err", "dZ", "CB/N", "Gb0", "Gb1");
        int n = 0;
        for (const Entry& x : e2) {
            const BumpNum p = BumpOf(gR[x.idx], 0, gBase[0]);
            const BumpNum p0 = BumpOf(gR[x.idx], 1, gBase[1]), p1 = BumpOf(gR[x.idx], 2, gBase[2]);
            P("%-3d %-2s %-36s | %6.3f %6.3f %6.3f %6.3f | %6.3f %6.3f %5.1f | %6.3f %6.3f | %+6.2f %6.3f | %6.3f %6.3f\n", ++n,
              gC[x.idx].grp.c_str(), gC[x.idx].name.c_str(), p.epsTag.e, p.epsTag.err, p.epsN, p.epsNCB, p.Gb, p.eGb, p.pull,
              p.epsH.e, p.epsH.err, p.Z - gBase[0].Z(), p.cb, p0.Gb, p1.Gb);
        }
        P("#\n# photon cuts for the bump (no tagged efficiency; epsH = edge-step height kept, GH = epsH/sqrt(epsN)):\n");
        for (size_t i = 1; i < gC.size(); ++i) {
            if (gC[i].grp != "k") continue;
            const BumpNum p = BumpOf(gR[i], 0, gBase[0]);
            P("#   %-28s epsN %.3f epsNCB %.3f epsH %.3f +- %.3f GH %.3f +- %.3f dZ %+.2f\n", gC[i].name.c_str(), p.epsN, p.epsNCB,
              p.epsH.e, p.epsH.err, p.epsH.e / std::sqrt(p.epsN), p.epsH.err / std::sqrt(p.epsN), p.Z - gBase[0].Z());
        }
    }
    // ---- is the ranking a fluctuation? choose on one fold, measure on the other ----
    P("#\n# CROSS-CHECK OF THE RANKING. The five best entries are chosen on fold A (working points also from fold A)\n");
    P("# and measured on fold B, where the choice cannot bias them.\n");
    for (int kind = 0; kind < 2; ++kind)
        for (int A = 0; A < 2; ++A) {
            const int fa = 1 + A, fb = 2 - A;
            std::vector<Entry> e = Entries(fa, kind);
            std::vector<Entry> e2;
            for (const Entry& x : e) if (!(kind == 1 && gC[x.idx].grp == "k")) e2.push_back(x);
            auto g = [&](int i, int f) { return kind == 0 ? PeakOf(gR[i], f, gBase[f]).G : BumpOf(gR[i], f, gBase[f]).Gb; };
            auto ge = [&](int i, int f) { return kind == 0 ? PeakOf(gR[i], f, gBase[f]).eG : BumpOf(gR[i], f, gBase[f]).eGb; };
            std::sort(e2.begin(), e2.end(), [&](const Entry& a, const Entry& b) { return g(a.idx, fa) > g(b.idx, fa); });
            P("#   %s, chosen on fold %d:\n", kind == 0 ? "PEAK" : "BUMP", A);
            for (int j = 0; j < 5 && j < (int)e2.size(); ++j)
                P("#     %-36s G on fold %d = %.3f +- %.3f ; on fold %d = %.3f +- %.3f\n", gC[e2[j].idx].name.c_str(), A,
                  g(e2[j].idx, fa), ge(e2[j].idx, fa), 1 - A, g(e2[j].idx, fb), ge(e2[j].idx, fb));
        }
    CloseTable();
}

// ============================================================================
// stage B: vertex study
// ============================================================================
TGraphErrors* GE(int col, int mk, double ms = 1.0) {
    TGraphErrors* g = new TGraphErrors();
    g->SetLineColor(col); g->SetMarkerColor(col); g->SetMarkerStyle(mk); g->SetMarkerSize(ms); g->SetLineWidth(2);
    return g;
}
void AddPt(TGraphErrors* g, double x, double y, double ey) {
    const int n = g->GetN(); g->SetPoint(n, x, y); g->SetPointError(n, 0, ey);
}
TH1D* Frame(const char* name, double x1, double x2, double y1, double y2, const char* xt, const char* yt) {
    TH1D* h = new TH1D(name, "", 100, x1, x2);
    h->SetDirectory(nullptr);
    h->SetMinimum(y1); h->SetMaximum(y2);
    h->GetXaxis()->SetTitle(xt); h->GetYaxis()->SetTitle(yt);
    h->GetXaxis()->SetNdivisions(506);
    h->Draw("AXIS");
    return h;
}
void StageVertex(TFile* fout) {
    OpenTable("extra_vertex.txt");
    const int nz = 50; const double zlo = -250, zhi = 0;       // 5 mm bins
    // [region: 0 all M4, 1 peak window, 2 bump window (bump start selection)][sample]
    std::vector<double> h[3][3];
    for (int a = 0; a < 3; ++a) for (int k = 0; k < 3; ++k) h[a][k].assign(nz, 0.);
    double out[3] = {0, 0, 0}, tot[3] = {0, 0, 0};
    for (int k = 0; k < 3; ++k)
        for (const auto& r : gP[0].s[k]) {
            tot[k] += 1;
            if (!(r.vz >= zlo && r.vz < zhi)) { out[k] += 1; continue; }
            const int i = (int)((r.vz - zlo) / (zhi - zlo) * nz);
            h[0][k][i] += 1;
            if (r.m4 >= scan::kPeakLo && r.m4 < scan::kPeakHi) h[1][k][i] += 1;
            if (StartBump(r) && r.m4 >= scan::kBumpLo && r.m4 < scan::kBumpHi) h[2][k][i] += 1;
        }
    P("#\n# B. VERTEX. Rows of the peak start selection (%s). z in mm, 10 mm bins.\n", kStartPeakTxt);
    P("# rows outside %.0f..%.0f mm: %.0f of %.0f (unlike-sign), %.0f of %.0f (e+e+), %.0f of %.0f (e-e-)\n", zlo, zhi, out[0],
      tot[0], out[1], tot[1], out[2], tot[2]);
    P("# CB = 2 sqrt(N++ N--) per z bin (no k factor). 'peak' = rows in the peak window, 'bump' = rows of the bump start selection in the bump window.\n#\n");
    P("%8s | %8s %7s %7s %6s | %7s %6s %6s %6s | %7s %6s %6s %6s\n", "z from", "all N", "N++", "N--", "CB/N", "peak N", "N++", "N--",
      "CB/N", "bump N", "N++", "N--", "CB/N");
    for (int i = 0; i < nz; i += 2) {
        double v[3][3];
        for (int a = 0; a < 3; ++a) for (int k = 0; k < 3; ++k) v[a][k] = h[a][k][i] + h[a][k][i + 1];
        P("%8.0f | %8.0f %7.0f %7.0f %6.3f | %7.0f %6.0f %6.0f %6.3f | %7.0f %6.0f %6.0f %6.3f\n", zlo + i * (zhi - zlo) / nz,
          v[0][0], v[0][1], v[0][2], v[0][0] > 0 ? 2 * std::sqrt(v[0][1] * v[0][2]) / v[0][0] : 0.,
          v[1][0], v[1][1], v[1][2], v[1][0] > 0 ? 2 * std::sqrt(v[1][1] * v[1][2]) / v[1][0] : 0.,
          v[2][0], v[2][1], v[2][2], v[2][0] > 0 ? 2 * std::sqrt(v[2][1] * v[2][2]) / v[2][0] : 0.);
    }
    // signal per z slice: peak (P_main) and tagged eta
    const std::vector<double> ze = {-250, -170, -160, -150, -140, -130, -120, -110, -100, -90, -50};
    std::vector<Sel> zs;
    for (size_t i = 0; i + 1 < ze.size(); ++i) {
        const float a = ze[i], b = ze[i + 1];
        zs.push_back([a, b](R r) { return r.vz >= a && r.vz < b; });
    }
    const std::vector<Spectra> spz = scan::FillMany(gP[0], zs);
    const std::vector<Spectra> spt = scan::FillMany(gT[0], zs, Opt5());
    const std::vector<Spectra> spb = scan::FillMany(gB[0][0], zs);
    P("#\n# signal per z slice: peak S (P_main, start selection) and tagged eta S5 (bump start selection); N = rows in the peak window,\n");
    P("# Nb = rows in the bump window. Ratios to the whole sample: (S/N)/(S/N)_all and (S5/Nb)/(S5/Nb)_all.\n#\n");
    P("%8s %8s | %7s %6s %5s %7s %6s %6s | %7s %6s %5s %8s %6s %6s\n", "z from", "to", "N", "S", "err", "S/N", "ratio", "err", "Nb",
      "S5", "err", "S5/Nb", "ratio", "err");
    TGraphErrors* gS = GE(kBlue + 1, 20, 1.3); TGraphErrors* gS5 = GE(kBlue + 1, 24, 1.3);
    const double r0 = gBase[0].S / gBase[0].N, r5 = gBase[0].S5 / gBase[0].bN;
    for (size_t i = 0; i + 1 < ze.size(); ++i) {
        const PeakFit f = scan::FitPeakMain(spz[i]);
        const PeakFit f5 = scan::FitPeakPois(spt[i], Cfg5());
        const double nb = spb[i].win[1].n[0];
        P("%8.0f %8.0f | %7.0f %6.0f %5.0f %7.4f %6.2f %6.2f | %7.0f %6.0f %5.0f %8.5f %6.2f %6.2f\n", ze[i], ze[i + 1], f.N, f.S, f.eS,
          f.N > 0 ? f.S / f.N : 0., f.N > 0 ? f.S / f.N / r0 : 0., f.N > 0 ? f.eS / f.N / r0 : 0., nb, f5.S, f5.eS,
          nb > 0 ? f5.S / nb : 0., nb > 0 ? f5.S / nb / r5 : 0., nb > 0 ? f5.eS / nb / r5 : 0.);
        const double zc = 0.5 * (std::max(ze[i], -180.) + std::min(ze[i + 1], -80.));
        if (f.N > 0) AddPt(gS, zc - 0.8, f.S / f.N / r0, f.eS / f.N / r0);
        if (nb > 0) AddPt(gS5, zc + 0.8, f5.S / nb / r5, f5.eS / nb / r5);
    }
    // ---- the rows below the proposed cut ----
    {
        const float zc = (float)kVzCommon;
        double a[4][2] = {{0, 0}, {0, 0}, {0, 0}, {0, 0}};   // [all cor_ok rows, peak start, peak window, bump window][below, all]
        double used[2][2] = {{0, 0}, {0, 0}}, multi[2] = {0, 0}, nrw[2] = {0, 0};   // peak window: [below, above]
        std::vector<float> zw;
        for (const auto& r : gAll.s[0]) { if (!r.cor_ok) continue; a[0][1] += 1; if (r.vz < zc) a[0][0] += 1; }
        for (const auto& r : gP[0].s[0]) {
            const int lo = r.vz < zc ? 0 : 1;
            a[1][1] += 1; a[1][0] += (lo == 0);
            if (r.m4 >= scan::kPeakLo && r.m4 < scan::kPeakHi) {
                a[2][1] += 1; a[2][0] += (lo == 0); zw.push_back(r.vz);
                nrw[lo] += 1; multi[lo] += (r.nrow > 1);
                used[lo][0] += (r.pip_isUsedVertex && r.pim_isUsedVertex); used[lo][1] += (r.ep_isUsedVertex && r.em_isUsedVertex);
            }
            if (StartBump(r) && r.m4 >= scan::kBumpLo && r.m4 < scan::kBumpHi) { a[3][1] += 1; a[3][0] += (lo == 0); }
        }
        std::sort(zw.begin(), zw.end());
        auto pc = [&](double q) { return zw.empty() ? 0.f : zw[(size_t)(q * (zw.size() - 1))]; };
        const std::vector<Sel> two = {[zc](R r) { return r.vz < zc; }, [zc](R r) { return r.vz >= zc; }};
        const std::vector<Spectra> s2 = scan::FillMany(gP[0], two);
        const PeakFit f1 = scan::FitPeakMain(s2[0]);
        const Spectra s5 = scan::Fill(gT[0], two[0], Opt5());
        PeakFit f5 = scan::FitPeakPois(s5, Cfg5());
        if (!f5.ok) f5 = scan::FitPeakMain(s5, Cfg5());
        P("#\n# THE ROWS BELOW THE PROPOSED CUT vz > %.0f mm (unlike-sign rows). The binding preselection is vz > -500: this cut is a NEW decision.\n", kVzCommon);
        P("#   share of rows with vz < %.0f: all cor_ok rows of the production %.1f %% ; peak start selection, all M4 %.1f %% ; peak window %.1f %% ; bump window %.1f %%\n",
          kVzCommon, 100 * a[0][0] / a[0][1], 100 * a[1][0] / a[1][1], 100 * a[2][0] / a[2][1], 100 * a[3][0] / a[3][1]);
        P("#   vz of the rows in the peak window: 1 %% point %.1f, 10 %% %.1f, 50 %% %.1f, 90 %% %.1f, 99 %% %.1f mm\n", pc(0.01), pc(0.10), pc(0.50), pc(0.90), pc(0.99));
        P("#   peak window, below / above the cut: CB share %.3f / %.3f ; both pions used in the vertex fit %.1f / %.1f %% ; both leptons used %.1f / %.1f %% ;\n",
          s2[0].win[0].cbShare(), s2[1].win[0].cbShare(), 100 * used[0][0] / std::max(1., nrw[0]), 100 * used[1][0] / std::max(1., nrw[1]),
          100 * used[0][1] / std::max(1., nrw[0]), 100 * used[1][1] / std::max(1., nrw[1]));
        P("#                                        rows of events with several rows %.1f / %.1f %%\n", 100 * multi[0] / std::max(1., nrw[0]), 100 * multi[1] / std::max(1., nrw[1]));
        P("#   signal below the cut: eta peak S = %.0f +- %.0f (of %.0f +- %.0f) ; tagged eta S5 = %.0f +- %.0f (of %.0f +- %.0f)\n", f1.S, f1.eS, gBase[0].S, gBase[0].eS,
          f5.S, f5.eS, gBase[0].S5, gBase[0].eS5);
        P("#   NOT identified: what these rows are (no simulation; the like-sign share points to a source without eta signal upstream of the target).\n");
        P("#   isUsedVertex: the pions fail it in %.0f %% of the background rows of the peak window but in about 5 %% of the tagged eta events (see the candidate tables);\n",
          100 * (1 - (used[0][0] + used[1][0]) / std::max(1., nrw[0] + nrw[1])));
        P("#   it is a reconstruction flag, so its efficiency for a cross section must come from simulation. Geometric alternative: |z_pi - vz| < 40 mm.\n");
    }
    P("#\n# threshold scans of the vertex cuts (vz >, vz <, rT <, vchi2 <): see extra_peak_candidates.txt and extra_bump_candidates.txt.\n");
    P("# beam spot used for rT: x = %.2f mm, y = %.2f mm (medians).\n", kVx0, kVy0);

    // ---- figure ----
    TCanvas* c = new TCanvas("c_vertex", "", 1800, 1300);
    c->Divide(2, 2);
    // (a) z distribution
    c->cd(1); gPad->SetLogy();
    TH1D* hz[3];
    const char* nm[3] = {"vz_us", "vz_cb", "vz_sig"};
    for (int j = 0; j < 3; ++j) { hz[j] = new TH1D(nm[j], "", nz, zlo, zhi); hz[j]->SetDirectory(nullptr); }
    for (int i = 0; i < nz; ++i) {
        const double n = h[0][0][i], cb = 2 * std::sqrt(h[0][1][i] * h[0][2][i]);
        hz[0]->SetBinContent(i + 1, n); hz[0]->SetBinError(i + 1, std::sqrt(n));
        hz[1]->SetBinContent(i + 1, cb); hz[1]->SetBinError(i + 1, std::sqrt(h[0][1][i] + h[0][2][i]));
        hz[2]->SetBinContent(i + 1, n - cb); hz[2]->SetBinError(i + 1, std::sqrt(n + h[0][1][i] + h[0][2][i]));
    }
    Deco(hz[0], kBlack, 20, 0.8, "vertex z [mm]", "rows / 5 mm"); Deco(hz[1], kRed + 1, 21, 0.8); Deco(hz[2], kBlue + 1, 24, 0.8);
    hz[0]->SetMinimum(30); hz[0]->SetMaximum(4 * hz[0]->GetMaximum());
    hz[0]->Draw("E1"); hz[1]->Draw("E1 SAME"); hz[2]->Draw("E1 SAME");
    VLine(-160, 30, hz[0]->GetMaximum() / 4, kGray + 2, 2);
    Head(Form("Vertex z, start selection of the peak (%s), all M4", kStartPeakTex));
    { TLegend* l = Leg(0.62, 0.70, 0.94, 0.90); l->AddEntry(hz[0], "all pairs (e^{+}e^{-})", "pe");
      l->AddEntry(hz[1], "CB = 2#sqrt{N_{++}N_{--}}", "pe"); l->AddEntry(hz[2], "all - CB", "pe"); l->Draw(); }
    Txt(0.20, 0.85, "dashed line: z = -160 mm", 0.038, kGray + 2);
    // (b) CB share against z
    c->cd(2);
    Frame("fr_cbz", -220, -60, 0, 0.62, "vertex z [mm]", "CB share = CB / all pairs");
    TGraphErrors* gc[3] = {GE(kBlack, 20, 1.0), GE(kRed + 1, 21, 1.1), GE(kMagenta + 2, 25, 1.1)};
    for (int a = 0; a < 3; ++a)
        for (int i = 0; i + 1 < nz; i += 2) {
            const double n = h[a][0][i] + h[a][0][i + 1], pp = h[a][1][i] + h[a][1][i + 1], mm = h[a][2][i] + h[a][2][i + 1];
            const double z = zlo + (i + 1) * (zhi - zlo) / nz;
            if (z < -220 || z > -60 || n < 50 || pp < 1 || mm < 1) continue;
            const double cb = 2 * std::sqrt(pp * mm);
            AddPt(gc[a], z + (a - 1) * 1.2, cb / n, cb / n * std::sqrt(0.25 / pp + 0.25 / mm + 1. / n));
        }
    for (int a = 0; a < 3; ++a) gc[a]->Draw("P SAME");
    VLine(-160, 0, 0.62, kGray + 2, 2);
    Head("CB share against vertex z (10 mm bins)");
    { TLegend* l = Leg(0.40, 0.68, 0.94, 0.90); l->AddEntry(gc[0], "all M4", "pe");
      l->AddEntry(gc[1], "rows in the peak window", "pe"); l->AddEntry(gc[2], "rows in the bump window (bump start selection)", "pe"); l->Draw(); }
    // (c) signal per row against z
    c->cd(3);
    Frame("fr_sz", -185, -75, -0.6, 2.6, "vertex z [mm]", "signal per row / average");
    HLine(-185, -75, 1, kGray + 1, 2); HLine(-185, -75, 0, kGray + 1, 1);
    gS->Draw("P SAME"); gS5->Draw("P SAME");
    VLine(-160, -0.6, 2.6, kGray + 2, 2);
    Head("Signal content against vertex z (first slice: z < -170, last: z > -90)");
    { TLegend* l = Leg(0.20, 0.72, 0.80, 0.90); l->AddEntry(gS, "#eta peak: S / rows in the peak window", "pe");
      l->AddEntry(gS5, "tagged #eta#rightarrow#pi^{+}#pi^{-}#pi^{0}: S_{5} / rows in the bump window", "pe"); l->Draw(); }
    // (d) scan of the lower edge
    c->cd(4);
    Frame("fr_vzscan", -205, -135, 0.78, 1.10, "lower cut on vertex z [mm]", "fraction kept");
    TGraphErrors* g1 = GE(kBlue + 1, 20, 1.3); TGraphErrors* g2 = GE(kBlue + 1, 24, 1.3);
    TGraphErrors* g3 = GE(kBlack, 21, 1.0); TGraphErrors* g4 = GE(kBlack, 25, 1.0);
    for (size_t i = 0; i < gC.size(); ++i) {
        if (gC[i].fam != "vz>") continue;
        const PeakNum p = PeakOf(gR[i], 0, gBase[0]); const BumpNum b = BumpOf(gR[i], 0, gBase[0]);
        AddPt(g1, gC[i].thr - 0.5, p.epsS.e, p.epsS.err); AddPt(g2, gC[i].thr + 0.5, b.epsTag.e, b.epsTag.err);
        AddPt(g3, gC[i].thr, std::sqrt(p.epsB), 0); AddPt(g4, gC[i].thr, std::sqrt(b.epsN), 0);
    }
    HLine(-205, -135, 1, kGray + 1, 2);
    g3->Draw("LP SAME"); g4->Draw("LP SAME"); g1->Draw("P SAME"); g2->Draw("P SAME");
    Head("Lower cut on vertex z: signal kept against #sqrt{background kept}");
    { TLegend* l = Leg(0.20, 0.17, 0.80, 0.42); l->AddEntry(g1, "#varepsilon_{S}, #eta peak (pass/fail fits)", "pe");
      l->AddEntry(g2, "#varepsilon_{tag}, tagged #eta#rightarrow#pi^{+}#pi^{-}#pi^{0}", "pe");
      l->AddEntry(g3, "#sqrt{#varepsilon_{B}}, peak window", "lp"); l->AddEntry(g4, "#sqrt{#varepsilon_{N}}, bump window", "lp"); l->Draw(); }
    Save(c, "extra_vertex_z", true);
    if (fout) { fout->cd(); for (int j = 0; j < 3; ++j) hz[j]->Write(); }
    CloseTable();
}

// ============================================================================
// stage K: the pi0 tag as a selection for eta -> pi+pi-pi0
// ============================================================================
// model variations of the tagged-eta yield: ranges x polynomial order x width +-0.5 MeV; Poisson fits
struct Spread5 { double half = 0, smin = 0, smax = 0; int nfit = 0, npass = 0; };
Spread5 TagSpread(const Spectra& sp) {
    Spread5 r;
    const double lo[4] = {0.47, 0.49, 0.50, 0.51}, hi[4] = {0.60, 0.62, 0.64, 0.66};
    for (int ir = 0; ir < 4; ++ir)
        for (int ord = 1; ord <= 3; ++ord)
            for (int iw = -1; iw <= 1; ++iw) {
                PeakCfg c = Cfg5(); c.lo = lo[ir]; c.hi = hi[ir]; c.order = ord; c.sigma = kM5Sigma + iw * 0.0005;
                const PeakFit f = scan::FitPeakPois(sp, c);
                r.nfit++;
                if (!f.ok || !(f.prob > scan::kVarPmin)) continue;
                if (r.npass == 0) r.smin = r.smax = f.S;
                r.smin = std::min(r.smin, f.S); r.smax = std::max(r.smax, f.S); r.npass++;
            }
    r.half = 0.5 * (r.smax - r.smin);
    return r;
}
void StageTagged(TFile* fout) {
    OpenTable("extra_tagged.txt");
    const scan::FillOpt o5 = Opt5();
    const PeakCfg c5 = Cfg5();
    const Spectra sp = scan::Fill(gT[0], scan::SelNone(), o5);
    P("#\n# K. PI0 TAG. Rows of the bump start selection (%s) with a photon giving %.3f < M(e+e-gamma) < %.3f\n", kStartBumpTxt, kPi0Lo, kPi0Hi);
    P("# (photon of g1..g3 closest to 0.140). Rows: %.0f unlike-sign, %.0f e+e+, %.0f e-e- (of %zu / %zu / %zu in the start selection).\n",
      sp.nsel[0], sp.nsel[1], sp.nsel[2], gB[0][0].s[0].size(), gB[0][0].s[1].size(), gB[0][0].s[2].size());
    P("# M5c = M(pi+pi-e+e-gamma) with the photon energy rescaled so that M(e+e-gamma) = m_pi0 (see M5c() in the macro).\n#\n");
    // pi0 in M(e+e-gamma): Gaussian + pol3, plain chi2 fit with ROOT
    TH1D* hg = new TH1D("meeg_bumpstart", "", 120, 0., 0.36); hg->SetDirectory(nullptr);
    for (const auto& r : gB[0][0].s[0]) if (r.g_idx > 0) hg->Fill(r.m_eeg);
    TF1* fg = new TF1("fg_pi0", "gausn(0)*0.003+pol3(3)", 0.09, 0.20);
    fg->SetParameters(90000, 0.142, 0.013, hg->GetBinContent(47), 0, 0, 0);
    fg->SetParLimits(1, 0.125, 0.155); fg->SetParLimits(2, 0.004, 0.025);
    hg->Fit(fg, "RQ0");
    P("# pi0 in M(e+e-gamma), bump start selection: position %.4f +- %.4f, width %.1f +- %.1f MeV, content %.0f +- %.0f, chi2/ndf %.1f/%d\n",
      fg->GetParameter(1), fg->GetParError(1), 1e3 * fg->GetParameter(2), 1e3 * fg->GetParError(2), fg->GetParameter(0),
      fg->GetParError(0), fg->GetChisquare(), fg->GetNDF());
    // eta in M5c: free and frozen
    const scan::FreeFit ff = scan::FitPeakFree(sp, c5, true);
    const PeakFit fz = scan::FitPeakPois(sp, c5);
    const PeakFit fm = scan::FitPeakMain(sp, c5);
    const Spread5 spr = TagSpread(sp);
    P("# eta in M5c, free Poisson fit (pol%d, %.2f-%.2f): position %.5f +- %.5f, width %.2f +- %.2f MeV, S5 = %.0f +- %.0f, chi2/ndf %.1f/%d\n",
      kM5Order, kM5Lo, kM5Hi, ff.mean, ff.eMean, 1e3 * ff.sigma, 1e3 * ff.eSigma, ff.S, ff.eS, ff.chi2, ff.ndf);
    P("# eta in M5c, frozen shape %.4f / %.1f MeV: S5 = %.0f +- %.0f (stat) +- %.0f (model spread: %d of %d fits kept, %.0f to %.0f)\n",
      kM5Mean, 1e3 * kM5Sigma, fz.S, fz.eS, spr.half, spr.npass, spr.nfit, spr.smin, spr.smax);
    P("#      window +-1.9 sigma: N = %.0f rows, CB = %.0f (share %.3f), S/B = %.3f, S5/err = %.1f, chi2/ndf %.1f/%d ; chi2 fit of N - CB: %.0f +- %.0f\n",
      fz.N, fz.CB, fz.N > 0 ? fz.CB / fz.N : 0., fz.SoB, fz.F, fz.chi2, fz.ndf, fm.S, fm.eS);
    // the same without the constraint, and with the narrow old window
    {
        scan::FillOpt o = o5; o.x = [](R r) { return (double)r.m5; }; o.windows.clear();
        PeakCfg c = c5; c.mean = 0.5540; c.sigma = 0.0110; c.lo = 0.47; c.hi = 0.64; c.order = 3; c.winLo = c.mean - 1.9 * c.sigma; c.winHi = c.mean + 1.9 * c.sigma;
        const Spectra s0 = scan::Fill(gT[0], scan::SelNone(), o);
        const scan::FreeFit f0 = scan::FitPeakFree(s0, c, true);
        PeakCfg cz = c; if (f0.ok) { cz.mean = f0.mean; cz.sigma = f0.sigma; cz.winLo = cz.mean - 1.9 * cz.sigma; cz.winHi = cz.mean + 1.9 * cz.sigma; }
        const PeakFit p0 = scan::FitPeakPois(s0, cz);
        if (f0.ok)
            P("# for comparison, plain M(pi+pi-e+e-gamma), same rows (pol3, 0.47-0.64, free shape): position %.4f, width %.1f +- %.1f MeV, S = %.0f +- %.0f, S/B = %.3f, chi2/ndf %.1f/%d\n",
              f0.mean, 1e3 * f0.sigma, 1e3 * f0.eSigma, f0.S, f0.eS, p0.SoB, f0.chi2, f0.ndf);
        else
            P("# for comparison, plain M(pi+pi-e+e-gamma), same rows (pol3, 0.47-0.64): free shape did not converge; fixed 0.5540 / 11.0 MeV: S = %.0f +- %.0f, S/B = %.3f, chi2/ndf %.1f/%d\n",
              p0.S, p0.eS, p0.SoB, p0.chi2, p0.ndf);
        const Spectra s1 = scan::Fill(gT[0], [](R r) { return r.m_eeg >= 0.125f && r.m_eeg < 0.145f; }, o);
        PeakCfg c1 = c; c1.mean = 0.5492; c1.sigma = 0.0088; c1.lo = 0.48; c1.hi = 0.64; c1.order = 3;
        const scan::FreeFit f1 = scan::FitPeakFree(s1, c1, true);
        PeakCfg c1z = c1; if (f1.ok) { c1z.mean = f1.mean; c1z.sigma = f1.sigma; } c1z.winLo = c1z.mean - 1.9 * c1z.sigma; c1z.winHi = c1z.mean + 1.9 * c1z.sigma;
        const PeakFit p1 = scan::FitPeakPois(s1, c1z);
        P("# old narrow window 0.125 < M(eeg) < 0.145, plain mass (pol3, 0.48-0.64, free shape): ok %d, position %.4f, width %.1f +- %.1f MeV, S = %.0f +- %.0f, S/B = %.3f\n",
          (int)f1.ok, f1.mean, 1e3 * f1.sigma, 1e3 * f1.eSigma, f1.S, f1.eS, p1.SoB);
        const Spectra s2 = scan::Fill(gT[0], [](R r) { return r.m_eeg >= 0.125f && r.m_eeg < 0.145f; }, o5);
        const PeakFit p2 = scan::FitPeakPois(s2, c5);
        P("# old narrow window, constrained mass, frozen shape: S5 = %.0f +- %.0f, S/B = %.3f\n", p2.S, p2.eS, p2.SoB);
    }
    // ---- the tagged yield by class of the lepton pair: signs only ----
    {
        const double eb[5] = {0., 0.008, 0.020, 0.030, 10.};
        const double kw[4] = {0.500, 0.232, 0.096, 0.172};     // pi0 -> e+e-gamma before acceptance (Kroll-Wada, form factor 1)
        std::vector<Sel> cs;
        for (int i = 0; i < 4; ++i) { const float a = eb[i], b = eb[i + 1]; cs.push_back([a, b](R r) { return r.mee >= a && r.mee < b; }); }
        cs.push_back([](R r) { return r.oa < 2.f; });
        cs.push_back([](R r) { return r.shared_ring == 1; });
        std::vector<Sel> pf;
        for (const Sel& q : cs) { pf.push_back(q); pf.push_back([q](R r) { return !q(r); }); }
        const std::vector<Spectra> sc = scan::FillMany(gT[0], pf, o5);
        const char* nm[6] = {"M(ee) < 8 MeV", "M(ee) 8-20 MeV", "M(ee) 20-30 MeV", "M(ee) > 30 MeV", "OA < 2 deg", "shared RICH ring"};
        P("#\n# tagged eta yield by class of the lepton pair (pass/fail Poisson fits of M5c, frozen shape); share = S5(class) / (S5(class) + S5(rest)).\n");
        for (int i = 0; i < 6; ++i) {
            PeakFit a = scan::FitPeakPois(sc[2 * i], c5), b = scan::FitPeakPois(sc[2 * i + 1], c5);
            if (!a.ok) a = scan::FitPeakMain(sc[2 * i], c5);
            if (!b.ok) b = scan::FitPeakMain(sc[2 * i + 1], c5);
            const Eff e = Ratio(a.S, a.eS, b.S, b.eS);
            P("#   %-18s S5 = %5.0f +- %3.0f, share %.3f +- %.3f%s\n", nm[i], a.S, a.eS, e.e, e.err,
              i < 4 ? Form(" ; pi0 Dalitz decay before acceptance: %.3f", kw[i]) : "");
        }
        P("# READING. These shares are SIGNS only. A pi0 -> gamma gamma with one photon converted to e+e- gives the same M(e+e-gamma) and the same kinematic\n");
        P("# edge in M4 as a Dalitz decay pi0 -> e+e-gamma; the acceptance against M(ee) is unknown without simulation. The shares neither prove nor exclude\n");
        P("# a conversion part. In this study 'eta -> pi+pi-pi0' and 'the bump' therefore mean: with an e+e- pair from the pi0, Dalitz decay or converted photon.\n");
    }
    // where the tagged eta lies in M4, and relative to the predicted edge.
    // M5c > M4, so in an M4 slice the spectrum starts at the slice; the eta peak is fitted locally
    // (frozen Gaussian + straight line, +-4.9 sigma, starting at least 2.5 MeV above the slice).
    P("#\n# tagged eta yield against M4: local Poisson fits of M5c (frozen Gaussian + straight line, 0.5175-0.5775 or from 2.5 MeV above the slice):\n");
    P("%8s %8s | %7s %6s %5s %9s\n", "M4 from", "to", "rows", "S5", "err", "chi2/ndf");
    std::vector<double> me;
    for (double m = 0.36; m < 0.5201; m += 0.01) me.push_back(m);
    std::vector<Sel> ms;
    for (size_t i = 0; i + 1 < me.size(); ++i) {
        const float a = me[i], b = me[i + 1];
        ms.push_back([a, b](R r) { return r.m4 >= a && r.m4 < b; });
    }
    ms.push_back([](R r) { return r.m4 < 0.52f; });
    ms.push_back([](R r) { return r.m4 >= (float)scan::kBumpLo && r.m4 < 0.52f; });
    ms.push_back([](R r) { const double e = scan::EdgeMeasured(r.mpipi); return e > 0 && r.m4 > e + 2. * scan::kEdgeSigma && r.m4 < 0.52f; });
    ms.push_back([](R r) { const double e = scan::EdgeMeasured(r.mpipi); return e > 0 && r.m4 <= e + 2. * scan::kEdgeSigma && r.m4 < 0.52f; });
    const std::vector<Spectra> spm = scan::FillMany(gT[0], ms, o5);
    TGraphErrors* gM = GE(kBlue + 1, 20, 1.2);
    // Poisson fit of the unsubtracted spectrum; if it does not converge (few rows), the chi2 fit of N - CB
    auto fitLocal = [&](const Spectra& q, double m4hi) {
        PeakCfg cc = c5; cc.lo = std::max(0.5175, m4hi + 0.0025); cc.hi = 0.5775; cc.order = 1;
        PeakFit f = scan::FitPeakPois(q, cc);
        if (!f.ok) f = scan::FitPeakMain(q, cc);
        return f;
    };
    double sSum = 0, vSum = 0;
    const size_t nsl = me.size() - 1;
    for (size_t i = 0; i < ms.size(); ++i) {
        const PeakFit f = fitLocal(spm[i], i < nsl ? me[i + 1] : 0.52);
        if (i < nsl) {
            P("%8.4f %8.4f | %7.0f %6.0f %5.0f %6.1f/%-2d%s\n", me[i], me[i + 1], spm[i].nsel[0], f.S, f.eS, f.chi2, f.ndf,
              !f.ok ? "  (fit failed)" : (f.method == 0 ? "  (chi2 fit of N - CB)" : ""));
            if (f.ok) { AddPt(gM, 0.5 * (me[i] + me[i + 1]), f.S / 4., f.eS / 4.); sSum += f.S; vSum += f.eS * f.eS; }
        } else if (i == nsl) {
            P("# all rows with M4 < 0.52, same local fit: S5 = %.0f +- %.0f (sum of the slices: %.0f +- %.0f ; fit of the whole sample: %.0f +- %.0f)\n",
              f.S, f.eS, sSum, std::sqrt(vSum), fz.S, fz.eS);
        } else if (i == nsl + 1) {
            P("# inside the bump window %.4f-0.5200: S5 = %.0f +- %.0f = %.2f of the tagged eta with M4 < 0.52 ; N - CB of the whole bump window = %.0f\n",
              scan::kBumpLo, f.S, f.eS, sSum > 0 ? f.S / sSum : 0., gBase[0].bN - gBase[0].bCB);
        } else if (i == nsl + 2) {
            // few rows and a background that is not a straight line (threshold M5c > M4): the fit is only indicative, so the raw counts are given too
            const double sg = kM5Sigma;
            const scan::WinSum wW = spm[i].WindowApprox(kM5Mean - 1.9 * sg, kM5Mean + 1.9 * sg);
            const scan::WinSum wL = spm[i].WindowApprox(kM5Mean - 4.9 * sg, kM5Mean - 1.9 * sg), wU = spm[i].WindowApprox(kM5Mean + 1.9 * sg, kM5Mean + 4.9 * sg);
            P("# rows more than 2 sigma (%.1f MeV) ABOVE the predicted kinematic edge of their M(pipi), M4 < 0.52: %.0f rows, of which %.0f in the eta window\n",
              2e3 * scan::kEdgeSigma, spm[i].nsel[0], wW.n[0]);
            P("#   of M5c (+-1.9 sigma), %.0f and %.0f in the 3 sigma wide windows below and above it: the spectrum only starts there (M5c > M4), no peak.\n",
              wL.n[0], wU.n[0]);
            P("#   So at most %.0f rows of tagged eta lie above the edge (the local fit is not usable here: %.0f +- %.0f).\n", wW.n[0], f.S, f.eS);
        } else {
            P("# rows below that line, M4 < 0.52: S5 = %.0f +- %.0f, rows %.0f\n", f.S, f.eS, spm[i].nsel[0]);
        }
    }
    P("#\n# comparison with the bump in M4 (start selection): edge steps sum %.0f +- %.0f rows / 2.5 MeV (Z = %.1f), no model-free yield;\n",
      gBase[0].hSum(), gBase[0].hErr(), gBase[0].Z());
    P("# the tag gives a countable yield S5 = %.0f +- %.0f +- %.0f with S/B = %.2f, but only for events whose photon reaches the ECAL.\n",
      fz.S, fz.eS, spr.half, fz.SoB);
    P("# The photon is the one of the pi0 that also gave the e+e- pair (Dalitz decay) or the partner of the converted photon: the tag does not tell which.\n");

    // ---- figure ----
    TCanvas* c = new TCanvas("c_tag", "", 1800, 1300);
    c->Divide(2, 2);
    c->cd(1);
    Deco(hg, kBlack, 20, 0.6, "M(e^{+}e^{-}#gamma) [GeV/c^{2}]", "rows / 3 MeV/c^{2}");
    hg->SetMinimum(0); hg->Draw("E1");
    fg->SetLineColor(kBlue + 1); fg->SetLineWidth(2); fg->Draw("SAME");
    VLine(kPi0Lo, 0, 0.8 * hg->GetMaximum(), kGray + 2, 2); VLine(kPi0Hi, 0, 0.8 * hg->GetMaximum(), kGray + 2, 2);
    Head("Photon closest to the #pi^{0}: M(e^{+}e^{-}#gamma), bump start selection, all pairs");
    Txt(0.52, 0.84, Form("#pi^{0}: %.4f, #sigma = %.1f MeV", fg->GetParameter(1), 1e3 * fg->GetParameter(2)), 0.04, kBlue + 1);
    Txt(0.52, 0.78, Form("tag window %.3f - %.3f", kPi0Lo, kPi0Hi), 0.04, kGray + 2);
    c->cd(2);
    TH1D* h5 = sp.Hist("us", "m5c_us"); TH1D* h5cb = sp.Hist("cb", "m5c_cb"); TH1D* h5s = sp.Hist("sig", "m5c_sig");
    const char* x5 = "M(#pi^{+}#pi^{-}e^{+}e^{-}#gamma), #pi^{0} constrained [GeV/c^{2}]";
    Deco(h5, kBlack, 20, 0.6, x5, kYt); Deco(h5cb, kRed + 1, 21, 0.6); Deco(h5s, kBlue + 1, 24, 0.6);
    const double mx5 = MaxIn(h5, 0.44, 0.70);
    h5->GetXaxis()->SetRangeUser(0.44, 0.70); h5->SetMinimum(0); h5->SetMaximum(1.35 * mx5);
    h5->Draw("E1"); h5cb->Draw("E1 SAME"); h5s->Draw("E1 SAME");
    Curve(fz, true, kBlack, 1)->Draw("L SAME"); Curve(fz, false, kBlack, 2)->Draw("L SAME");
    Head("Tagged #eta#rightarrow#pi^{+}#pi^{-}#pi^{0}: bump start selection and a #pi^{0} photon");
    { TLegend* l = Leg(0.20, 0.68, 0.55, 0.90); l->AddEntry(h5, "all pairs", "pe"); l->AddEntry(h5cb, "CB", "pe");
      l->AddEntry(h5s, "all - CB", "pe"); l->Draw(); }
    Txt(0.58, 0.40, Form("S_{5} = %.0f #pm %.0f #pm %.0f", fz.S, fz.eS, spr.half), 0.042);
    Txt(0.58, 0.34, Form("S/B = %.2f, #sigma = %.1f MeV", fz.SoB, 1e3 * kM5Sigma), 0.042);
    Txt(0.58, 0.28, "fit: Poisson, all pairs", 0.036);
    c->cd(3);
    {
        scan::FillOpt o = o5; o.x = [](R r) { return (double)r.m5; }; o.windows.clear();
        const Spectra s0 = scan::Fill(gT[0], scan::SelNone(), o);
        TH1D* ha = s0.Hist("us", "m5plain_us");
        Deco(ha, kGray + 2, 24, 0.6, "M(#pi^{+}#pi^{-}e^{+}e^{-}#gamma) [GeV/c^{2}]", kYt);
        ha->GetXaxis()->SetRangeUser(0.44, 0.70); ha->SetMinimum(0); ha->SetMaximum(1.35 * mx5);
        TH1D* hb = (TH1D*)h5->Clone("m5c_us2"); hb->GetXaxis()->SetRangeUser(0.44, 0.70);
        ha->Draw("E1"); hb->Draw("E1 SAME");
        Head("Effect of the #pi^{0} constraint (all pairs, same rows)");
        TLegend* l = Leg(0.20, 0.74, 0.62, 0.90); l->AddEntry(ha, "measured photon energy", "pe");
        l->AddEntry(hb, "photon energy from M(e^{+}e^{-}#gamma) = m_{#pi^{0}}", "pe"); l->Draw();
        if (fout) { fout->cd(); ha->Write(); }
    }
    c->cd(4);
    {
        const Spectra sb = scan::Fill(gB[0][0], scan::SelNone());
        TH1D* hb = sb.Hist("sig", "bump_sig_start");
        Deco(hb, kGray + 2, 24, 0.6, kXt, "rows / 2.5 MeV/c^{2}");
        double mx = 0; for (int i = 0; i < gM->GetN(); ++i) mx = std::max(mx, gM->GetY()[i]);
        const double sc = 1.3 * mx / hb->GetMaximum();
        hb->Scale(sc); hb->GetXaxis()->SetRangeUser(0.30, 0.60); hb->SetMinimum(-0.05 * mx); hb->SetMaximum(1.6 * mx);
        hb->Draw("HIST"); gM->Draw("P SAME");
        VLine(scan::kBumpLo, 0, 1.3 * mx, kGray + 2, 2); VLine(scan::kBumpHi, 0, 1.3 * mx, kGray + 2, 2);
        Head("Where the tagged #eta lies in M4 (#eta peak fitted in M4 slices)");
        TLegend* l = Leg(0.20, 0.74, 0.70, 0.90); l->AddEntry(gM, "tagged #eta yield S_{5} per 2.5 MeV (local fits, M4 < 0.52)", "pe");
        l->AddEntry(hb, Form("all - CB of the bump start selection #times %.4f", sc), "l"); l->Draw();
        Txt(0.20, 0.68, "dashed lines: bump window", 0.036, kGray + 2);
    }
    Save(c, "extra_tagged_eta", true);
    if (fout) { fout->cd(); hg->Write(); h5->Write(); h5cb->Write(); h5s->Write(); }
    CloseTable();
}


// ============================================================================
// stage C: the best cuts applied together; model spread; two-fold check
// ============================================================================
Cand AndOf(const std::vector<int>& idx) {
    Cand c; c.grp = "+";
    std::vector<Sel> v;
    for (int i : idx) { v.push_back(gC[i].cut); c.name += (c.name.empty() ? "" : " & ") + gC[i].name; }
    c.cut = [v](R r) { for (const auto& q : v) if (!q(r)) return false; return true; };
    return c;
}
const std::vector<int>* FamIdx(const std::vector<std::pair<std::string, std::vector<int>>>& fam, const char* name) {
    for (const auto& f : fam) if (f.first == name) return &f.second;
    return nullptr;
}
// spread of the ratio S(selection) / S(start) over the model variations kept in both
double RatioSpread(const scan::Spread& a, const scan::Spread& b) {
    double mn = 1e9, mx = -1e9; int n = 0;
    for (size_t i = 0; i < a.fits.size() && i < b.fits.size(); ++i) {
        const PeakFit& x = a.fits[i]; const PeakFit& y = b.fits[i];
        if (!x.ok || !y.ok || !(x.prob > scan::kVarPmin) || !(y.prob > scan::kVarPmin) || !(y.S > 0)) continue;
        mn = std::min(mn, x.S / y.S); mx = std::max(mx, x.S / y.S); n++;
    }
    return n > 1 ? 0.5 * (mx - mn) : 0.;
}
struct BeforeAfter { Sel cutPeak, cutPeakVeto, cutBump; std::string namePeak, nameBump; };
BeforeAfter gBA;

// ---- pass/fail numbers of one cut on an arbitrary data set (the sample without kinematic cuts
//      is not covered by Evaluate()) ------------------------------------------------------------
struct PRow {
    PeakFit f, fp, l1, l2; scan::Spread sv;
    Eff epsS; double epsB = 1, G = 1, eG = 0;
    Eff epsTag; double S5 = 0, eS5 = 0; bool hasTag = false;
    double pullMean = 0, pullRms = 0, pullMean8 = 0, pullRms8 = 0, pullMin = 0;   // displaced-peak scan
};
PRow PeakRowOf(const Data& d, const Data* dTag, const Sel& cut) {
    PRow r;
    const std::vector<Sel> v = {cut, [cut](R q) { return !cut(q); }};
    const std::vector<Spectra> sp = scan::FillMany(d, v);
    r.f = scan::FitPeakMain(sp[0]); r.fp = scan::FitPeakPois(sp[0]);
    r.l1 = scan::FitPeakLikeSign(sp[0], 1); r.l2 = scan::FitPeakLikeSign(sp[0], 2);
    r.sv = scan::PeakModelSpread(sp[0], PeakCfg(), true);
    {   // displaced-peak scan as in scan_peak.C: 13 positions 0.575 ... 0.695 without a known peak; range and window move with the peak
        double s1 = 0, s2 = 0, t1 = 0, t2 = 0; r.pullMin = 1e9;
        for (int i = 0; i < 13; ++i) {
            const double pos = 0.575 + 0.010 * i, sh = pos - scan::kPeakMean;
            PeakCfg c; c.mean = pos; c.lo = scan::kFitLo + sh; c.hi = scan::kFitHi + sh; c.winLo = scan::kPeakLo + sh; c.winHi = scan::kPeakHi + sh;
            const PeakFit q = scan::FitPeakMain(sp[0], c);
            const double pull = q.eS > 0 ? q.S / q.eS : 0;
            s1 += pull; s2 += pull * pull; r.pullMin = std::min(r.pullMin, pull);
            if (i < 8) { t1 += pull; t2 += pull * pull; }
        }
        r.pullMean = s1 / 13; r.pullRms = std::sqrt(s2 / 13); r.pullMean8 = t1 / 8; r.pullRms8 = std::sqrt(t2 / 8);
    }
    r.epsS.e = 1; r.epsS.err = 0;
    if (sp[1].nsel[0] > 0) {
        const PeakFit ff = scan::FitPeakMain(sp[1]);
        r.epsS = Ratio(r.f.S, r.f.eS, ff.S, ff.eS);
        const double ba = r.f.N - r.f.Swin, bb = ff.N - ff.Swin;
        r.epsB = (ba + bb) > 0 ? ba / (ba + bb) : 0;
        if (r.epsB > 0) { r.G = r.epsS.e / std::sqrt(r.epsB); r.eG = r.epsS.err / std::sqrt(r.epsB); }
    }
    if (dTag) {
        const std::vector<Spectra> st = scan::FillMany(*dTag, v, Opt5());
        const PeakFit a = scan::FitPeakPois(st[0], Cfg5());
        r.S5 = a.ok ? a.S : 0; r.eS5 = a.ok ? a.eS : 0;
        double sb = 0, eb = 0;
        if (st[1].nsel[0] > 0) { const PeakFit b = scan::FitPeakPois(st[1], Cfg5()); if (b.ok) { sb = b.S; eb = b.eS; } }
        r.epsTag = Ratio(r.S5, r.eS5, sb, eb);
        r.hasTag = true;
    }
    return r;
}
// sculpting flags that use the pass fraction only (T: tilt, K: kink), for the choice of cut families on one fold
bool SculptTK(int i) {
    const std::string f = PeakFlags(i);
    return f.find('T') != std::string::npos || f.find('K') != std::string::npos;
}
// The fixed procedure "choose the cut families": the best entry (single cut or working point of a scan,
// both taken on fold set fa) of every task group by G; groups with pull > 1 ranked by G; the three best.
std::vector<int> ChooseFamilies(int fa, int kind) {
    const std::vector<Entry> e = Entries(fa, kind);
    auto g = [&](int i) { return kind == 0 ? PeakOf(gR[i], fa, gBase[fa]).G : BumpOf(gR[i], fa, gBase[fa]).Gb; };
    auto pl = [&](int i) { return kind == 0 ? PeakOf(gR[i], fa, gBase[fa]).pull : BumpOf(gR[i], fa, gBase[fa]).pull; };
    std::map<std::string, int> best;
    for (const Entry& x : e) {
        const std::string& grp = gC[x.idx].grp;
        if (kind == 1 && grp == "k") continue;          // photon cuts have no tagged efficiency
        if (kind == 0 && SculptTK(x.idx)) continue;
        const auto it = best.find(grp);
        if (it == best.end() || g(x.idx) > g(it->second)) best[grp] = x.idx;
    }
    std::vector<int> v;
    for (const auto& kv : best) if (pl(kv.second) > 1.) v.push_back(kv.second);
    std::sort(v.begin(), v.end(), [&](int a, int b) { return g(a) > g(b); });
    if (v.size() > 3) v.resize(3);
    return v;
}

void StageCombine() {
    const std::vector<std::pair<std::string, std::vector<int>>> fam = Families(gC);
    OpenTable("extra_combination.txt");
    const std::vector<int>* fvz = FamIdx(fam, "vz>");
    const std::vector<int>* frv = FamIdx(fam, "rv_pipi<");
    const int iU = Idx("isUsedVertex both pions"), iG = Idx("no pi0 photon");
    const int iV = Idx(Form("vz > %.0f", kVzCommon)), iZ = Idx("|z - vz| (both pi) < 40");
    P("#\n# STATUS OF THE CUTS (after the review)\n");
    P("#   PROPOSED : vz > %.0f mm and both pions used in the vertex fit (isUsedVertex). Both are confirmed by the photon-tagged eta -> pi+pi-pi0\n", kVzCommon);
    P("#              control sample (epsTag), which does not depend on the peak fit. Two points need a decision of the person in charge:\n");
    P("#              (1) the binding vertex cut is z > -500; vz > %.0f is a NEW, tighter cut; (2) isUsedVertex is a reconstruction flag: its efficiency\n", kVzCommon);
    P("#              must come from simulation. The geometric alternative |z_pi - vz| < 40 mm (both pions) is listed below; it is easier to model.\n");
    P("#   OPTIONAL : no pi0 photon (veto of rows whose stored photon gives %.3f <= M(e+e-gamma) < %.3f). It has no control sample; alone its G is\n", kPi0Lo, kPi0Hi);
    P("#              about 1.5 sigma above 1. It raises S/B; a gain in S/sigma_S is NOT shown.\n");
    P("#   F = S/sigma_S is statistical only (frozen peak shape) and serves to compare selections; every yield is S +- stat +- model spread.\n");
    P("#   CB/N: CB = 2 sqrt(N++ N--) without a k factor, an estimate.\n");
    {
        int nreg = 0;
        const int a = FamilyWP(gR, *fvz, gBase[0], 0, 0, &nreg), b = FamilyWP(gR, *fvz, gBase[0], 0, 1);
        P("#   vz threshold by the rule (centre of the region within 1 unit of the best figure of merit): peak %s (%d grid points), bump %s; one common value %.0f is used.\n",
          a >= 0 ? gC[a].name.c_str() : "none", nreg, b >= 0 ? gC[b].name.c_str() : "none", kVzCommon);
    }
    // ------------------------------------------------------------------ peak
    const Sel cV = gC[iV].cut, cU = gC[iU].cut, cG = gC[iG].cut, cZ = gC[iZ].cut;
    const Sel cVU = [cV, cU](R r) { return cV(r) && cU(r); };
    const Sel cVZ = [cV, cZ](R r) { return cV(r) && cZ(r); };
    const Sel cVUG = [cV, cU, cG](R r) { return cV(r) && cU(r) && cG(r); };
    const Sel cVZG = [cV, cZ, cG](R r) { return cV(r) && cZ(r) && cG(r); };
    gBA.cutPeak = cVU; gBA.cutPeakVeto = cVUG; gBA.namePeak = gC[iV].name + " & " + gC[iU].name;
    {
        const char* head = "%-70s | %6s %6s | %5s %4s %6s | %5s %4s | %5s %6s | %6s %5s %5s %6s | %6s %5s | %6s %5s %6s | %9s %9s | %8s\n";
        auto line = [&](const char* name, const PRow& r, const PRow& s0, bool isStart, bool tagOk) {
            P("%-70s | %6.0f %6.3f | %5.0f %4.0f %6.0f | %5.0f %4.0f | %5.2f %6.4f | ", name, r.f.N, r.f.CB / r.f.N, r.f.S, r.f.eS, r.sv.half, r.fp.S, r.fp.eS, r.f.F, r.f.SoB);
            if (isStart) P("%6s %5s %5s %6s | %6s %5s | %6s %5s %6s | ", "1", "", "", "1", "1", "", "1", "", "1");
            else {
                P("%6.3f %5.3f %5.3f %6.3f | %6.3f %5.3f | ", r.epsS.e, r.epsS.err, RatioSpread(r.sv, s0.sv), r.epsB, r.G, r.eG);
                if (tagOk && r.hasTag) P("%6.3f %5.3f %6.3f | ", r.epsTag.e, r.epsTag.err, r.epsB > 0 ? r.epsTag.e / std::sqrt(r.epsB) : 0.);
                else P("%6s %5s %6s | ", "n/a", "", "n/a");
            }
            P("%4.0f+-%-3.0f %4.0f+-%-3.0f | %5.1f/%d\n", r.l1.S, r.l1.eS, r.l2.S, r.l2.eS, r.f.chi2, r.f.ndf);
        };
        P("#\n# PEAK. Columns: N, CB/N in the peak window; S = P_main, stat = fit error, spread = half (max - min) over the model variations with p > 0.05;\n");
        P("# Spois = Poisson fit of N; F = S/stat; epsS = S_pass/(S_pass + S_fail), rSpr = model spread of S/S_start, epsB = background kept, G = epsS/sqrt(epsB);\n");
        P("# epsTag = efficiency of the same cut for the photon-tagged eta -> pi+pi-pi0 (control sample), Gtag = epsTag/sqrt(epsB) = the gain predicted\n");
        P("# from the control sample alone (n/a for photon cuts); null++ / null-- = the peak fit on the like-sign spectra.\n");
        // ---- block A: on top of the kinematic nominal ----
        P("#\n# A. On top of the kinematic nominal of the peak: %s. Tag sample: bump start selection.\n", kStartPeakTxt);
        P(head, "selection", "N", "CB/N", "S", "stat", "spread", "Spois", "stat", "F", "S/B", "epsS", "err", "rSpr", "epsB", "G", "err", "epsTag", "err", "Gtag", "null++", "null--", "chi2/ndf");
        const PRow a0 = PeakRowOf(gP[0], &gT[0], scan::SelNone());
        line("nominal (start)", a0, a0, true, true);
        line(gC[iV].name.c_str(), PeakRowOf(gP[0], &gT[0], cV), a0, false, true);
        line(gC[iU].name.c_str(), PeakRowOf(gP[0], &gT[0], cU), a0, false, true);
        const PRow aVU = PeakRowOf(gP[0], &gT[0], cVU);
        line(Form("PROPOSED: %s & %s", gC[iV].name.c_str(), gC[iU].name.c_str()), aVU, a0, false, true);
        line(Form("alternative: %s & %s", gC[iV].name.c_str(), gC[iZ].name.c_str()), PeakRowOf(gP[0], &gT[0], cVZ), a0, false, true);
        line(gC[iG].name.c_str(), PeakRowOf(gP[0], &gT[0], cG), a0, false, false);
        const PRow aVUG = PeakRowOf(gP[0], &gT[0], cVUG);
        line("proposed + optional veto (no pi0 photon)", aVUG, a0, false, false);
        line("alternative + optional veto (no pi0 photon)", PeakRowOf(gP[0], &gT[0], cVZG), a0, false, false);
        {
            // what the veto removes: the rows of the proposed selection that FAIL it
            const Spectra sf = scan::Fill(gP[0], [cVU, cG](R r) { return cVU(r) && !cG(r); });
            const PeakFit ff = scan::FitPeakMain(sf);
            const double blind = aVU.f.N > 0 ? aVU.f.S * ff.N / aVU.f.N : 0;
            P("#   the veto on top of the proposed selection: the %.0f vetoed rows of the peak window hold S = %.0f +- %.0f; a veto blind to the signal would remove %.0f.\n",
              ff.N, ff.S, ff.eS, blind);
            P("#   So the veto enriches the peak, but F changes by %+.2f only (%.2f -> %.2f), with an error of about %.2f: a gain is not shown.\n",
              aVUG.f.F - aVU.f.F, aVU.f.F, aVUG.f.F, std::sqrt(std::max(0., 2 * (1 - aVUG.f.eS / aVU.f.eS))));
        }
        // ---- block B: no kinematic cut ----
        P("#\n# B. Without any kinematic cut (all cor_ok rows). Tag sample: all rows with a pi0 photon.\n");
        P(head, "selection", "N", "CB/N", "S", "stat", "spread", "Spois", "stat", "F", "S/B", "epsS", "err", "rSpr", "epsB", "G", "err", "epsTag", "err", "Gtag", "null++", "null--", "chi2/ndf");
        const PRow b0 = PeakRowOf(gAll, &gTagAll, scan::SelNone());
        line("no kinematic cut", b0, b0, true, true);
        const PRow bVU = PeakRowOf(gAll, &gTagAll, cVU);
        line(Form("%s & %s", gC[iV].name.c_str(), gC[iU].name.c_str()), bVU, b0, false, true);
        line("the same + optional veto (no pi0 photon)", PeakRowOf(gAll, &gTagAll, cVUG), b0, false, false);
        P("#\n# displaced-peak scan (the same fit at 13 positions 0.575-0.695 without a known peak; pull = S/sigma_S; expected mean 0, rms 1; first 8 = range ends below 0.72):\n");
        {
            const PRow* rr[5] = {&b0, &a0, &aVU, &aVUG, &bVU};
            const char* nn[5] = {"no kinematic cut", "nominal (start)", "PROPOSED (nominal + vertex cuts)", "proposed + optional veto", "vertex cuts, no kinematic cut"};
            for (int i = 0; i < 5; ++i)
                P("#   %-34s all 13: mean %+.2f, rms %.2f, most negative %+.2f ; first 8: mean %+.2f, rms %.2f\n", nn[i], rr[i]->pullMean, rr[i]->pullRms,
                  rr[i]->pullMin, rr[i]->pullMean8, rr[i]->pullRms8);
            P("#   A mean away from 0 or an rms above 1 shows that pol2 does not describe the background everywhere: F is not a significance.\n");
        }
        P("#   The gain comes from the vertex cuts, not from the kinematic cuts: F %.2f (no cut) -> %.2f (two vertex cuts alone) ; %.2f (kinematic nominal) -> %.2f (with the two vertex cuts).\n",
          b0.f.F, bVU.f.F, a0.f.F, aVU.f.F);
        P("#   Error of a difference of F between nested samples: sqrt(2 (1 - sigma_pass/sigma_all)), here %.2f.\n", std::sqrt(std::max(0., 2 * (1 - aVU.f.eS / a0.f.eS))));
        // ---- isBest == 1 ----
        P("#\n# explicit comparison, rows with isBest == 1 only (kinematic nominal):\n");
        for (int q = 0; q < 3; ++q) {
            const Sel cut = q == 0 ? scan::SelNone() : (q == 1 ? cVU : cVUG);
            const Spectra sp = scan::Fill(gP[0], [cut](R r) { return r.isBest == 1 && cut(r); });
            const scan::Spread sv = scan::PeakModelSpread(sp);
            const PeakFit f = scan::FitPeakMain(sp), fp = scan::FitPeakPois(sp);
            P("%-70s | %6.0f %6.3f | %5.0f %4.0f %6.0f | %5.0f %4.0f | %5.2f %6.4f |\n",
              q == 0 ? "isBest == 1: nominal (start)" : (q == 1 ? "isBest == 1: proposed" : "isBest == 1: proposed + optional veto"), f.N,
              f.CB / f.N, f.S, f.eS, sv.half, fp.S, fp.eS, f.F, f.SoB);
        }
        // ---- two-fold checks ----
        P("#\n# PEAK, two-fold check 1 (threshold only). The vz threshold is chosen on fold A by the rule (centre of the region F >= Fmax - 1);\n");
        P("# vz and isUsedVertex of both pions are then measured on fold B (pass/fail fits on fold B only).\n");
        double sumS0 = 0, sumS1 = 0, v0 = 0, v1 = 0;
        for (int A = 0; A < 2; ++A) {
            const int fa = 1 + A, fb = 2 - A;
            const int jV = FamilyWP(gR, *fvz, gBase[fa], fa, 0);
            std::vector<int> t2;
            if (jV >= 0) t2.push_back(jV);
            t2.push_back(iU);
            std::vector<Cand> c1(1, AndOf(t2));
            const std::vector<CR> r1 = Evaluate(c1);
            const PeakNum p = PeakOf(r1[0], fb, gBase[fb]);
            const Half& a = r1[0].h[fb][0];
            P("#   chosen on fold %d: %s\n#      on fold %d: start S = %.0f +- %.0f, F = %.2f, S/B = %.4f ; after S = %.0f +- %.0f, F = %.2f, S/B = %.4f, CB/N = %.3f ; epsS = %.3f +- %.3f, epsB = %.3f, G = %.3f +- %.3f\n",
              A, c1[0].name.c_str(), 1 - A, gBase[fb].S, gBase[fb].eS, gBase[fb].F(), gBase[fb].SoB, a.S, a.eS, a.F(), a.SoB, p.cb,
              p.epsS.e, p.epsS.err, p.epsB, p.G, p.eG);
            sumS0 += gBase[fb].S; v0 += gBase[fb].eS * gBase[fb].eS; sumS1 += a.S; v1 += a.eS * a.eS;
        }
        P("#   sum of the two cross-measured halves: start S = %.0f +- %.0f, after S = %.0f +- %.0f ; ratio of F: %.3f\n", sumS0, std::sqrt(v0),
          sumS1, std::sqrt(v1), (sumS1 / std::sqrt(v1)) / (sumS0 / std::sqrt(v0)));
        P("#   This check does NOT cover the choice of the two cuts among the %zu entries of the ranking: that is check 2.\n", Entries(0, 0).size());
        P("#\n# PEAK, two-fold check 2 (choice of the cut FAMILIES). Fixed procedure applied to fold A alone: best entry of every task group by G (single cut, or working\n");
        P("# point of a scan from fold A; entries with a tilted or kinked pass fraction excluded), groups with pull > 1 ranked by G, the three best combined.\n");
        P("# The combination is measured on fold B, which took no part in the choice.\n");
        for (int A = 0; A < 2; ++A) {
            const int fa = 1 + A, fb = 2 - A;
            const std::vector<int> ch = ChooseFamilies(fa, 0);
            if (ch.empty()) { P("#   chosen on fold %d: nothing with pull > 1\n", A); continue; }
            std::vector<Cand> c1(1, AndOf(ch));
            const std::vector<CR> r1 = Evaluate(c1);
            const PeakNum pa = PeakOf(r1[0], fa, gBase[fa]), pb = PeakOf(r1[0], fb, gBase[fb]);
            P("#   chosen on fold %d: %s\n#      G on fold %d (choosing, biased) = %.3f +- %.3f, F %.2f -> %.2f ; G on fold %d (independent) = %.3f +- %.3f, F %.2f -> %.2f, epsS = %.3f +- %.3f, epsB = %.3f\n",
              A, c1[0].name.c_str(), A, pa.G, pa.eG, gBase[fa].F(), r1[0].h[fa][0].F(), 1 - A, pb.G, pb.eG, gBase[fb].F(), r1[0].h[fb][0].F(),
              pb.epsS.e, pb.epsS.err, pb.epsB);
        }
    }
    // ------------------------------------------------------------------ bump
    // both pions used in the vertex fit and the common vertex cut; the pion crossing point (family rv_pipi <,
    // working point by Z) is listed for information only (it adds about one error)
    {
        const int iR = FamilyWP(gR, *frv, gBase[0], 0, 1);
        std::vector<int> trio;
        trio.push_back(iU); trio.push_back(iV);
        if (iR >= 0) trio.push_back(iR);
        const int n = (int)trio.size();
        std::vector<Cand> cc;
        for (int m = 1; m < (1 << n); ++m) {
            std::vector<int> idx;
            for (int j = 0; j < n; ++j) if (m & (1 << j)) idx.push_back(trio[j]);
            cc.push_back(AndOf(idx));
        }
        { std::vector<int> t3 = {iV, iZ}; cc.push_back(AndOf(t3)); }       // geometric alternative
        const std::vector<CR> rc = Evaluate(cc);
        const size_t iProp = 2;          // mask 3 = iU & iV
        gBA.cutBump = cc[iProp].cut; gBA.nameBump = cc[iProp].name;
        const Half& b = gBase[0];
        P("#\n# BUMP: on top of the kinematic nominal of the bump: %s. PROPOSED = %s.\n", kStartBumpTxt, cc[iProp].name.c_str());
        P("# N, CB/N, N-CB: bump window ; epsN, epsNCB: rows and N - CB kept ; epsTag, S5, S/B5: tagged eta (S5 with stat error and model spread) ;\n");
        P("# Gb = epsTag/sqrt(epsN) ; h1, h2, h3: edge-step height [rows / 2.5 MeV] in M(pipi) %.2f-%.2f, %.2f-%.2f, %.2f-%.2f ; Z = sum/error ; epsH = step height kept ;\n",
          kSlice[0], kSlice[1], kSlice[1], kSlice[2], kSlice[2], kSlice[3]);
        P("# GH = epsH/sqrt(epsN): the same gain from the edge step (it has a model spread of 2-5 %% on top of the statistical error).\n#\n");
        P("%-58s | %6s %5s %6s | %5s %6s | %6s %5s | %6s %5s | %5s %4s %4s %5s | %4s %3s %4s %3s %4s %3s | %5s | %6s %5s | %6s %5s\n", "selection", "N", "CB/N",
          "N-CB", "epsN", "epsNCB", "epsTag", "err", "Gb", "err", "S5", "stat", "sprd", "S/B5", "h1", "err", "h2", "err", "h3", "err", "Z", "epsH", "err", "GH", "err");
        {
            const Spread5 t0 = TagSpread(scan::Fill(gT[0], scan::SelNone(), Opt5()));
            P("%-58s | %6.0f %5.3f %6.0f | %5s %6s | %6s %5s | %6s %5s | %5.0f %4.0f %4.0f %5.3f | %4.0f %3.0f %4.0f %3.0f %4.0f %3.0f | %5.2f | %6s %5s | %6s %5s\n",
              "nominal (start)", b.bN, b.bCB / b.bN, b.bN - b.bCB, "1", "1", "1", "", "1", "", b.S5, b.eS5, t0.half, b.SoB5, b.h[0], b.eh[0], b.h[1],
              b.eh[1], b.h[2], b.eh[2], b.Z(), "1", "", "1", "");
        }
        for (size_t i = 0; i < cc.size(); ++i) {
            const BumpNum p = BumpOf(rc[i], 0, gBase[0]);
            const Half& a = rc[i].h[0][0];
            const Spread5 t1 = TagSpread(scan::Fill(gT[0], cc[i].cut, Opt5()));
            const std::string nm = (i == iProp ? "PROPOSED: " : (i + 1 == cc.size() ? "alternative: " : "")) + cc[i].name;
            P("%-58s | %6.0f %5.3f %6.0f | %5.3f %6.3f | %6.3f %5.3f | %6.3f %5.3f | %5.0f %4.0f %4.0f %5.3f | %4.0f %3.0f %4.0f %3.0f %4.0f %3.0f | %5.2f | %6.3f %5.3f | %6.3f %5.3f\n",
              nm.c_str(), a.bN, p.cb, a.bN - a.bCB, p.epsN, p.epsNCB, p.epsTag.e, p.epsTag.err, p.Gb, p.eGb, a.S5, a.eS5, t1.half,
              a.SoB5, a.h[0], a.eh[0], a.h[1], a.eh[1], a.h[2], a.eh[2], a.Z(), p.epsH.e, p.epsH.err,
              p.epsN > 0 ? p.epsH.e / std::sqrt(p.epsN) : 0., p.epsN > 0 ? p.epsH.err / std::sqrt(p.epsN) : 0.);
        }
        // per M(pipi) slice: edge step kept against tagged eta kept (proposed selection)
        {
            const Half& a = rc[iProp].h[0][0];
            P("#\n# proposed selection per M(pipi) slice: edge step kept (h after / h before) against the tagged eta kept (pass/fail fits of the tagged rows of the slice)\n");
            std::vector<Sel> ts;
            for (int s = 0; s < 3; ++s) {
                const float lo = kSlice[s], hi = kSlice[s + 1];
                const Sel cut = cc[iProp].cut;
                ts.push_back([lo, hi, cut](R r) { return r.mpipi >= lo && r.mpipi < hi && cut(r); });
                ts.push_back([lo, hi, cut](R r) { return r.mpipi >= lo && r.mpipi < hi && !cut(r); });
            }
            const std::vector<Spectra> st = scan::FillMany(gT[0], ts, Opt5());
            for (int s = 0; s < 3; ++s) {
                PeakFit fa = scan::FitPeakPois(st[2 * s], Cfg5()), ff = scan::FitPeakPois(st[2 * s + 1], Cfg5());
                if (!fa.ok) fa = scan::FitPeakMain(st[2 * s], Cfg5());
                if (!ff.ok) ff = scan::FitPeakMain(st[2 * s + 1], Cfg5());
                const Eff et = Ratio(fa.S, fa.eS, ff.S, ff.eS);
                P("#   M(pipi) %.2f-%.2f: step %.0f +- %.0f -> %.0f +- %.0f, kept %.2f +- %.2f ; tagged eta kept %.3f +- %.3f\n", kSlice[s], kSlice[s + 1],
                  b.h[s], b.eh[s], a.h[s], a.eh[s], b.h[s] > 0 ? a.h[s] / b.h[s] : 0.,
                  b.h[s] > 0 ? std::sqrt(std::max(0., (1 - 2 * a.h[s] / b.h[s]) * a.eh[s] * a.eh[s] + a.h[s] * a.h[s] / (b.h[s] * b.h[s]) * b.eh[s] * b.eh[s])) / b.h[s] : 0.,
                  et.e, et.err);
            }
        }
        P("#\n# BUMP, two-fold check 1 (threshold only). The vz threshold is chosen on fold A by the rule (centre of the region Z >= Zmax - 1);\n");
        P("# vz and isUsedVertex of both pions are measured on fold B.\n");
        for (int A = 0; A < 2; ++A) {
            const int fa = 1 + A, fb = 2 - A;
            const int jV = FamilyWP(gR, *fvz, gBase[fa], fa, 1);
            std::vector<int> t2(1, iU);
            if (jV >= 0) t2.push_back(jV);
            std::vector<Cand> c1(1, AndOf(t2));
            const std::vector<CR> r1 = Evaluate(c1);
            const BumpNum p = BumpOf(r1[0], fb, gBase[fb]);
            const Half& a = r1[0].h[fb][0];
            P("#   chosen on fold %d: %s\n#      on fold %d: epsTag = %.3f +- %.3f, epsN = %.3f, epsNCB = %.3f, Gb = %.3f +- %.3f ; step sum %.0f +- %.0f -> %.0f +- %.0f, Z %.2f -> %.2f, epsH = %.3f +- %.3f ; CB/N %.3f -> %.3f\n",
              A, c1[0].name.c_str(), 1 - A, p.epsTag.e, p.epsTag.err, p.epsN, p.epsNCB, p.Gb, p.eGb, gBase[fb].hSum(), gBase[fb].hErr(),
              a.hSum(), a.hErr(), gBase[fb].Z(), a.Z(), p.epsH.e, p.epsH.err, gBase[fb].bCB / gBase[fb].bN, p.cb);
        }
        P("#\n# BUMP, two-fold check 2 (choice of the cut FAMILIES; the same fixed procedure as for the peak, figure of merit Gb from the tagged eta):\n");
        for (int A = 0; A < 2; ++A) {
            const int fa = 1 + A, fb = 2 - A;
            const std::vector<int> ch = ChooseFamilies(fa, 1);
            if (ch.empty()) { P("#   chosen on fold %d: nothing with pull > 1\n", A); continue; }
            std::vector<Cand> c1(1, AndOf(ch));
            const std::vector<CR> r1 = Evaluate(c1);
            const BumpNum pa = BumpOf(r1[0], fa, gBase[fa]), pb = BumpOf(r1[0], fb, gBase[fb]);
            P("#   chosen on fold %d: %s\n#      Gb on fold %d (choosing, biased) = %.3f +- %.3f ; Gb on fold %d (independent) = %.3f +- %.3f, epsTag = %.3f +- %.3f, epsN = %.3f, GH = %.3f +- %.3f\n",
              A, c1[0].name.c_str(), A, pa.Gb, pa.eGb, 1 - A, pb.Gb, pb.eGb, pb.epsTag.e, pb.epsTag.err, pb.epsN,
              pb.epsN > 0 ? pb.epsH.e / std::sqrt(pb.epsN) : 0., pb.epsN > 0 ? pb.epsH.err / std::sqrt(pb.epsN) : 0.);
        }
    }
    CloseTable();
}

// ---- fitted spectra before and after the best combination ----------------------
void DrawPeakColumn(TCanvas* c, int col, int ncol, const Spectra& sp, const char* title, TFile* fout, const char* tag) {
    const PeakFit f = scan::FitPeakMain(sp);
    const scan::Spread sv = scan::PeakModelSpread(sp);
    TH1D* hN = sp.Hist("us", Form("pk_%s_us", tag)); TH1D* hC = sp.Hist("cb", Form("pk_%s_cb", tag)); TH1D* hS = sp.Hist("sig", Form("pk_%s_sig", tag));
    Deco(hN, kBlack, 20, 0.6); Deco(hC, kRed + 1, 21, 0.6); Deco(hS, kBlue + 1, 24, 0.6);
    // overview, log scale
    c->cd(col); gPad->SetLogy();
    TH1D* a = (TH1D*)hN->Clone(Form("pk_%s_usA", tag));
    a->GetXaxis()->SetRangeUser(0.40, 0.70); a->SetMinimum(0.3 * std::max(1., hC->GetBinContent(hC->FindBin(0.41)))); a->SetMaximum(6 * MaxIn(hN, 0.40, 0.70));
    a->Draw("E1"); hC->Draw("E1 SAME"); hS->Draw("E1 SAME");
    VLine(scan::kPeakLo, a->GetMinimum(), MaxIn(hN, 0.40, 0.70), kGray + 2, 2, 1); VLine(scan::kPeakHi, a->GetMinimum(), MaxIn(hN, 0.40, 0.70), kGray + 2, 2, 1);
    Head(title, 0.17, 0.94, 0.040);
    { TLegend* l = Leg(0.50, 0.17, 0.94, 0.36, 0.042); l->AddEntry(hN, "all pairs", "pe"); l->AddEntry(hC, "CB = 2#sqrt{N_{++}N_{--}}", "pe");
      l->AddEntry(hS, "all - CB", "pe"); l->Draw(); }
    // linear zoom on N - CB with the fit
    c->cd(col + ncol);
    TH1D* z = (TH1D*)hS->Clone(Form("pk_%s_sigZ", tag));
    z->GetXaxis()->SetRangeUser(0.50, 0.62);
    double mn = 1e30, mx = 0;
    for (int i = z->FindBin(0.5226); i <= z->FindBin(0.6099); ++i) { mn = std::min(mn, z->GetBinContent(i)); mx = std::max(mx, z->GetBinContent(i)); }
    z->SetMinimum(mn - 0.25 * (mx - mn)); z->SetMaximum(mx + 0.65 * (mx - mn));
    z->GetYaxis()->SetTitle("all - CB, rows / 2.5 MeV/c^{2}");
    z->Draw("E1");
    Curve(f, true, kBlue + 1, 1)->Draw("L SAME"); Curve(f, false, kBlue + 1, 2)->Draw("L SAME");
    Head("linear zoom: all - CB with the fit (Gaussian + pol2, frozen shape)", 0.17, 0.94, 0.038);
    Txt(0.20, 0.85, Form("S = %.0f #pm %.0f (stat) #pm %.0f (model)", f.S, f.eS, sv.half), 0.042, kBlue + 1);
    Txt(0.20, 0.79, Form("S/#sigma_{S} = %.1f (stat. only, not a significance)", f.F), 0.042, kBlue + 1);
    Txt(0.20, 0.73, Form("S/B = %.3f,  CB share = %.3f (no k factor)", f.SoB, f.CB / f.N), 0.042, kBlue + 1);
    Txt(0.20, 0.67, Form("#chi^{2}/ndf = %.1f/%d", f.chi2, f.ndf), 0.042, kBlue + 1);
    // polynomial subtracted
    c->cd(col + 2 * ncol);
    TH1D* r = Resid(sp, sp.sig, sp.sigVar, f, Form("pk_%s_res", tag));
    Deco(r, kBlue + 1, 20, 0.7);
    r->GetXaxis()->SetRangeUser(0.50, 0.62); r->GetYaxis()->SetTitle("all - CB - polynomial");
    r->SetMinimum(-450); r->SetMaximum(950);
    r->Draw("E1");
    TGraph* g = new TGraph();
    for (int i = 0; i <= 300; ++i) { const double m = f.cfg.lo + (f.cfg.hi - f.cfg.lo) * i / 300.; g->SetPoint(i, m, f.Total(m) - f.Bkg(m)); }
    g->SetLineColor(kBlue + 1); g->SetLineWidth(2); g->Draw("L SAME");
    HLine(0.50, 0.62, 0, kGray + 1, 2);
    Head("the same after subtracting the fitted polynomial", 0.17, 0.94, 0.038);
    if (fout) { fout->cd(); hN->Write(); hC->Write(); hS->Write(); }
}
void DrawBumpRow(TCanvas* c, int row, const Sel& cut, const char* title, TFile* fout, const char* tag) {
    // column 1: M4 of the bump start selection
    const Spectra sp = scan::Fill(gB[0][0], cut);
    TH1D* hN = sp.Hist("us", Form("bp_%s_us", tag)); TH1D* hC = sp.Hist("cb", Form("bp_%s_cb", tag)); TH1D* hS = sp.Hist("sig", Form("bp_%s_sig", tag));
    Deco(hN, kBlack, 20, 0.6); Deco(hC, kRed + 1, 21, 0.6); Deco(hS, kBlue + 1, 24, 0.6);
    c->cd(row * 5 + 1); gPad->SetLeftMargin(0.21);
    hN->GetYaxis()->SetTitleOffset(1.9);
    hN->GetXaxis()->SetRangeUser(0.40, 0.62); hN->SetMinimum(0); hN->SetMaximum(2.0 * MaxIn(hN, 0.40, 0.62));
    hN->Draw("E1"); hC->Draw("E1 SAME"); hS->Draw("E1 SAME");
    VLine(scan::kBumpLo, 0, MaxIn(hN, 0.40, 0.62), kGray + 2, 2, 1); VLine(scan::kBumpHi, 0, MaxIn(hN, 0.40, 0.62), kGray + 2, 2, 1);
    Head(title, 0.21, 0.94, 0.040);
    { TLegend* l = Leg(0.24, 0.72, 0.50, 0.90, 0.040); l->AddEntry(hN, "all pairs", "pe"); l->AddEntry(hC, "CB", "pe"); l->AddEntry(hS, "all - CB", "pe"); l->Draw(); }
    Txt(0.52, 0.86, "bump window (dashed):", 0.036);
    Txt(0.52, 0.81, Form("N - CB = %.0f", sp.win[1].sig()), 0.036);
    Txt(0.52, 0.76, Form("CB share %.3f (no k)", sp.win[1].cbShare()), 0.036);
    Txt(0.52, 0.69, "fall above 0.52: made by the", 0.032, kRed + 1);
    Txt(0.52, 0.65, "M(#pi#pi) < 0.41 cut, not the bump", 0.032, kRed + 1);
    // columns 2-4: slices with the edge-step fit
    for (int s = 0; s < 3; ++s) {
        const Spectra ss = scan::Fill(gB[0][1 + s], cut);
        StepCfg sc; sc.mpLo = kSlice[s]; sc.mpHi = kSlice[s + 1];
        const StepFit f = scan::FitEdgeStep(ss, sc);
        TH1D* h = ss.Hist("sig", Form("bp_%s_sl%d", tag, s));
        Deco(h, kBlue + 1, 20, 0.6);
        c->cd(row * 5 + 2 + s); gPad->SetLeftMargin(0.21);
        h->GetYaxis()->SetTitleOffset(1.9);
        h->GetXaxis()->SetRangeUser(0.45, 0.58); h->GetYaxis()->SetTitle("all - CB, rows / 2.5 MeV/c^{2}");
        const double mx = MaxIn(h, 0.45, 0.58);
        h->SetMinimum(0.35 * mx); h->SetMaximum(1.30 * mx);
        h->Draw("E1");
        if (f.ok) {
            TGraph* g1 = new TGraph(); TGraph* g0 = new TGraph();
            for (size_t i = 0; i < f.tx.size(); ++i) { g1->SetPoint((int)i, f.tx[i], f.fy[i]); g0->SetPoint((int)i, f.tx[i], f.fy0[i]); }
            g1->SetLineColor(kBlue + 1); g1->SetLineWidth(2); g0->SetLineColor(kBlue + 1); g0->SetLineStyle(2); g0->SetLineWidth(2);
            g1->Draw("L SAME"); g0->Draw("L SAME");
            VLine(f.edge, 0.35 * mx, 1.12 * mx, kGray + 2, 2, 1);
            Txt(0.24, 0.85, Form("step h = %.0f #pm %.0f rows / 2.5 MeV", f.h, f.eh), 0.042, kBlue + 1);
            Txt(0.24, 0.79, Form("#chi^{2}/ndf = %.1f/%d,  #Delta#chi^{2} = %.0f", f.chi2, f.ndf, f.dchi2), 0.042, kBlue + 1);
        }
        Head(Form("%.2f < M(#pi#pi) < %.2f: all - CB, edge-step fit", kSlice[s], kSlice[s + 1]), 0.21, 0.94, 0.040);
        if (fout) { fout->cd(); h->Write(); }
    }
    // column 5: tagged eta
    const Spectra st = scan::Fill(gT[0], cut, Opt5());
    const PeakFit f5 = scan::FitPeakPois(st, Cfg5());
    TH1D* tN = st.Hist("us", Form("tg_%s_us", tag)); TH1D* tC = st.Hist("cb", Form("tg_%s_cb", tag));
    Deco(tN, kBlack, 20, 0.6, "M(#pi^{+}#pi^{-}e^{+}e^{-}#gamma), #pi^{0} constrained [GeV/c^{2}]", kYt); Deco(tC, kRed + 1, 21, 0.6);
    c->cd(row * 5 + 5); gPad->SetLeftMargin(0.21);
    tN->GetYaxis()->SetTitleOffset(1.9);
    tN->GetXaxis()->SetRangeUser(0.46, 0.66); tN->SetMinimum(0); tN->SetMaximum(1.35 * MaxIn(tN, 0.46, 0.66));
    tN->Draw("E1"); tC->Draw("E1 SAME");
    Curve(f5, true, kBlack, 1)->Draw("L SAME"); Curve(f5, false, kBlack, 2)->Draw("L SAME");
    Head("rows with a #pi^{0} photon: tagged #eta#rightarrow#pi^{+}#pi^{-}#pi^{0}", 0.21, 0.94, 0.040);
    Txt(0.24, 0.85, Form("S_{5} = %.0f #pm %.0f,  S/B = %.2f", f5.S, f5.eS, f5.SoB), 0.042);
    { TLegend* l = Leg(0.62, 0.20, 0.94, 0.34, 0.042); l->AddEntry(tN, "all pairs", "pe"); l->AddEntry(tC, "CB", "pe"); l->Draw(); }
    if (fout) { fout->cd(); hN->Write(); hC->Write(); hS->Write(); tN->Write(); tC->Write(); }
}
void StageBeforeAfter(TFile* fout) {
    {
        TCanvas* c = new TCanvas("c_pba", "", 2700, 1900);
        c->Divide(3, 3);
        DrawPeakColumn(c, 1, 3, scan::Fill(gP[0], scan::SelNone()), Form("PEAK, kinematic nominal: %s", kStartPeakTex), fout, "start");
        DrawPeakColumn(c, 2, 3, scan::Fill(gP[0], gBA.cutPeak), Form("PROPOSED: + z_{vtx} > %.0f mm + both #pi in the vertex fit", kVzCommon), fout, "comb");
        DrawPeakColumn(c, 3, 3, scan::Fill(gP[0], gBA.cutPeakVeto), "OPTIONAL, not verified: + no #pi^{0} photon", fout, "veto");
        Save(c, "extra_peak_before_after", true);
    }
    {
        TCanvas* c = new TCanvas("c_bba", "", 3600, 1400);
        c->Divide(5, 2);
        DrawBumpRow(c, 0, scan::SelNone(), Form("BUMP, nominal: %s", kStartBumpTex), fout, "start");
        DrawBumpRow(c, 1, gBA.cutBump, Form("PROPOSED: + z_{vtx} > %.0f mm, both #pi in vtx fit", kVzCommon), fout, "comb");
        Save(c, "extra_bump_before_after", true);
    }
}

// ============================================================================
// stage F: summary figure, threshold scans, sculpting, sideband-subtracted distributions
// ============================================================================
int GrpColor(const std::string& g) {
    static const int col[14] = {kBlack, kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 2,
                                kViolet + 5, kYellow + 3, kPink + 7, kAzure + 7, kSpring + 5, kGray + 2, kTeal + 3};
    const int i = g.empty() ? 0 : (g[0] - 'a');
    return (i >= 0 && i < 14) ? col[i] : kBlack;
}
void StageSummary() {
    const char* grpName[14] = {"a missing mass", "b vertex", "c pion with ring", "d shared ring", "e crossing points",
                               "f isUsedVertex", "g lepton, pion quality", "h evtStart", "i p_{T}, rapidity", "j event class",
                               "k photons", "l hyp_chi2", "m lepton system, p, #theta", "n other (sideband study)"};
    const std::vector<std::pair<std::string, std::vector<int>>> fam = Families(gC);
    TCanvas* c = new TCanvas("c_sum", "", 2200, 2200);
    c->Divide(2, 2);
    for (int pad = 0; pad < 4; ++pad) {
        const int kind = pad % 2; const bool zoom = pad >= 2;
        c->cd(pad + 1);
        const double x1 = zoom ? 0.60 : 0.28, x2 = zoom ? 1.02 : 1.03, y1 = zoom ? 0.78 : 0.28, y2 = zoom ? 1.07 : 1.12;
        Frame(Form("fr_sum%d", pad), x1, x2, y1, y2,
              kind == 0 ? "#varepsilon_{B}: background kept in the peak window" : "#varepsilon_{N}: rows kept in the bump window",
              kind == 0 ? "#varepsilon_{S}: #eta peak content kept" : "#varepsilon_{tag}: tagged #eta#rightarrow#pi^{+}#pi^{-}#pi^{0} kept");
        TGraph* gs = new TGraph(); TGraph* gd = new TGraph();
        for (int i = 0; i <= 200; ++i) { const double x = x1 + (1.0 - x1) * i / 200.; gs->SetPoint(i, x, std::sqrt(x)); if (x > y1) gd->SetPoint(gd->GetN(), x, x); }
        gs->SetLineWidth(3); gs->Draw("L SAME"); gd->SetLineStyle(2); gd->SetLineColor(kGray + 1); gd->Draw("L SAME");
        auto xy = [&](int i, double& x, double& y, double& ey) {
            if (kind == 0) { const PeakNum p = PeakOf(gR[i], 0, gBase[0]); x = p.epsB; y = p.epsS.e; ey = p.epsS.err; }
            else { const BumpNum p = BumpOf(gR[i], 0, gBase[0]); x = p.epsN; y = p.epsTag.e; ey = p.epsTag.err; }
        };
        auto inside = [&](double x, double y) { return x > x1 && x < x2 && y > y1 && y < y2; };
        // scan families: thin line through the thresholds; working point = full circle with error bar
        for (const auto& f : fam) {
            const int col = GrpColor(gC[f.second[0]].grp);
            TGraph* g = new TGraph();
            g->SetPoint(0, 1., 1.);
            for (int i : f.second) { double x, y, ey; xy(i, x, y, ey); if (inside(x, y)) g->SetPoint(g->GetN(), x, y); else break; }
            g->SetLineColor(col); g->SetMarkerColor(col); g->SetMarkerStyle(24); g->SetMarkerSize(0.5); g->Draw("LP SAME");
            const int w = FamilyWP(gR, f.second, gBase[0], 0, kind);
            if (w >= 0) {
                double x, y, ey; xy(w, x, y, ey);
                if (inside(x, y)) { TGraphErrors* gw = GE(col, 20, 1.3); AddPt(gw, x, y, zoom ? ey : 0.); gw->Draw("P SAME"); }
            }
        }
        // single cuts: square with the letter of the group
        for (size_t i = 1; i < gC.size(); ++i) {
            if (!gC[i].fam.empty()) continue;
            if (kind == 1 && gC[i].grp == "k") continue;
            double x, y, ey; xy((int)i, x, y, ey);
            if (!inside(x, y)) continue;
            const int col = GrpColor(gC[i].grp);
            TGraphErrors* g = GE(col, 21, 1.1); AddPt(g, x, y, ey); g->Draw("P SAME");
            TLatex l; l.SetTextSize(0.028); l.SetTextFont(62); l.SetTextColor(col);
            l.DrawLatex(x + 0.004 * (x2 - x1) / 0.42, y - 0.035 * (y2 - y1), gC[i].grp.c_str());
        }
        Head(Form("%s%s: candidate cuts on top of %s", kind == 0 ? "PEAK" : "BUMP", zoom ? ", zoom" : "", kind == 0 ? kStartPeakTex : kStartBumpTex), 0.17, 0.94, 0.036);
        if (!zoom) {
            Txt(0.20, 0.87, "above the thick curve #varepsilon_{signal} = #sqrt{#varepsilon_{background}}: the cut helps", 0.030);
            Txt(0.20, 0.835, "dashed: signal and background reduced alike", 0.030, kGray + 2);
            Txt(0.20, 0.80, "squares + letter: single cuts; lines: threshold scans;", 0.030);
            Txt(0.20, 0.765, "full circle: working point of a scan", 0.030);
            for (int g = 0; g < 14; ++g) {
                const std::string q(1, (char)('a' + g));
                if (kind == 1 && g == 10) continue;
                Txt(0.62, 0.56 - 0.028 * g, grpName[g], 0.027, GrpColor(q));
            }
        } else {
            // the six best entries by G, numbered
            std::vector<Entry> e = Entries(0, kind);
            std::vector<Entry> e2;
            for (const Entry& q : e) if (!(kind == 1 && gC[q.idx].grp == "k")) e2.push_back(q);
            auto g = [&](int i) { return kind == 0 ? PeakOf(gR[i], 0, gBase[0]).G : BumpOf(gR[i], 0, gBase[0]).Gb; };
            auto ge = [&](int i) { return kind == 0 ? PeakOf(gR[i], 0, gBase[0]).eG : BumpOf(gR[i], 0, gBase[0]).eGb; };
            std::sort(e2.begin(), e2.end(), [&](const Entry& a, const Entry& b) { return g(a.idx) > g(b.idx); });
            TPave* box = new TPave(0.40, 0.17, 0.955, 0.48, 1, "NDC");
            box->SetFillColor(kWhite); box->SetLineColor(kGray + 1); box->Draw();
            Txt(0.42, 0.44, "best by G = #varepsilon_{signal} / #sqrt{#varepsilon_{background}} (picked on the full sample):", 0.026);
            for (int j = 0; j < 6 && j < (int)e2.size(); ++j) {
                double x, y, ey; xy(e2[j].idx, x, y, ey);
                if (inside(x, y)) { TLatex l; l.SetTextSize(0.030); l.SetTextFont(62); l.DrawLatex(x - 0.016, y + 0.006, Form("%d", j + 1)); }
                Txt(0.42, 0.40 - 0.034 * j, Form("%d  %s:  %.3f #pm %.3f", j + 1, gC[e2[j].idx].name.c_str(), g(e2[j].idx), ge(e2[j].idx)), 0.028);
            }
        }
    }
    Save(c, "extra_summary_eps", true);
}
void StageScans() {
    const std::vector<std::pair<std::string, std::vector<int>>> fam = Families(gC);
    const char* famP[6] = {"vz>", "rT<", "|z-vz|pi<", "|dz_pipi|<", "MM4<", "pt4>"};
    const char* famB[6] = {"vz>", "rv_pipi<", "|z-vz|pi<", "|dz_pipi|<", "rT<", "MM4<"};
    const char* xt[2][6] = {{"lower cut on vertex z [mm]", "upper cut on r_{T} of the vertex [mm]", "upper cut on |z_{#pi} - z_{vertex}| [mm]",
                             "upper cut on |dz_{#pi#pi}| (crossing point) [mm]", "upper cut on the missing mass [GeV/c^{2}]", "lower cut on p_{T}(4) [MeV/c]"},
                            {"lower cut on vertex z [mm]", "upper cut on r of the #pi#pi crossing [mm]", "upper cut on |z_{#pi} - z_{vertex}| [mm]",
                             "upper cut on |dz_{#pi#pi}| (crossing point) [mm]", "upper cut on r_{T} of the vertex [mm]", "upper cut on the missing mass [GeV/c^{2}]"}};
    for (int kind = 0; kind < 2; ++kind) {
        TCanvas* c = new TCanvas(Form("c_scan%d", kind), "", 2400, 2400);
        c->Divide(3, 4);
        for (int j = 0; j < 6; ++j) {
            const std::vector<int>* idx = FamIdx(fam, kind == 0 ? famP[j] : famB[j]);
            if (!idx) continue;
            double x1 = 1e30, x2 = -1e30;
            for (int i : *idx) { x1 = std::min(x1, gC[i].thr); x2 = std::max(x2, gC[i].thr); }
            const double pad = 0.08 * (x2 - x1); x1 -= pad; x2 += pad;
            TGraphErrors* gSig = GE(kBlue + 1, 20, 1.2); TGraphErrors* gCtl = GE(kBlue + 1, 24, 1.2); TGraphErrors* gBkg = GE(kBlack, 21, 0.9);
            TGraphErrors* gF = GE(kBlack, 20, 1.1);
            double ymin = 1, fmax = kind == 0 ? gBase[0].F() : gBase[0].Z();
            const double off = 0.012 * (x2 - x1);
            for (int i : *idx) {
                const PeakNum p = PeakOf(gR[i], 0, gBase[0]); const BumpNum b = BumpOf(gR[i], 0, gBase[0]);
                if (kind == 0) {
                    AddPt(gSig, gC[i].thr - off, p.epsS.e, p.epsS.err); AddPt(gCtl, gC[i].thr + off, b.epsTag.e, b.epsTag.err);
                    AddPt(gBkg, gC[i].thr, std::sqrt(p.epsB), 0); AddPt(gF, gC[i].thr, p.F, 0);
                    ymin = std::min(ymin, std::min(p.epsS.e - p.epsS.err, std::sqrt(p.epsB))); fmax = std::max(fmax, p.F);
                } else {
                    AddPt(gSig, gC[i].thr - off, b.epsTag.e, b.epsTag.err); AddPt(gCtl, gC[i].thr + off, b.epsH.e, b.epsH.err);
                    AddPt(gBkg, gC[i].thr, std::sqrt(b.epsN), 0); AddPt(gF, gC[i].thr, b.Z, 0);
                    ymin = std::min(ymin, std::min(b.epsH.e - b.epsH.err, std::sqrt(b.epsN))); fmax = std::max(fmax, b.Z);
                }
            }
            const int w = FamilyWP(gR, *idx, gBase[0], 0, kind);
            const int padE = (j < 3) ? j + 1 : j + 4, padF = padE + 3;
            c->cd(padE);
            Frame(Form("fr_sc%d_%d", kind, j), x1, x2, std::max(0.2, ymin - 0.05), 1.16, xt[kind][j], "fraction kept");
            HLine(x1, x2, 1, kGray + 1, 2);
            gBkg->Draw("LP SAME"); gCtl->Draw("P SAME"); gSig->Draw("P SAME");
            if (w >= 0) VLine(gC[w].thr, std::max(0.2, ymin - 0.05), 1.05, kGreen + 2, 1, 2);
            Head(Form("%s: %s", kind == 0 ? "PEAK" : "BUMP", w >= 0 ? Form("working point %s", gC[w].name.c_str()) : "working point: no cut"), 0.17, 0.94, 0.042);
            if (j == 0) {
                TLegend* l = Leg(0.22, 0.17, 0.80, 0.42, 0.045);
                if (kind == 0) { l->AddEntry(gSig, "#varepsilon_{S}, #eta peak (pass/fail fits)", "pe"); l->AddEntry(gCtl, "#varepsilon_{tag}, tagged #eta#rightarrow#pi^{+}#pi^{-}#pi^{0}", "pe");
                                 l->AddEntry(gBkg, "#sqrt{#varepsilon_{B}}, peak window", "lp"); }
                else { l->AddEntry(gSig, "#varepsilon_{tag}, tagged #eta#rightarrow#pi^{+}#pi^{-}#pi^{0}", "pe"); l->AddEntry(gCtl, "#varepsilon_{H}, edge-step height", "pe");
                       l->AddEntry(gBkg, "#sqrt{#varepsilon_{N}}, bump window", "lp"); }
                l->Draw();
            }
            c->cd(padF);
            const double f0 = kind == 0 ? gBase[0].F() : gBase[0].Z();
            double fmin = f0; for (int q = 0; q < gF->GetN(); ++q) fmin = std::min(fmin, gF->GetY()[q]);
            Frame(Form("fr_sf%d_%d", kind, j), x1, x2, std::max(fmin, fmax - 4.) - 0.4, fmax + 0.8, xt[kind][j],
                  kind == 0 ? "F = S/#sigma_{S} (stat. only) after the cut" : "Z = summed edge step / error");
            HLine(x1, x2, f0, kGray + 2, 2, 2); HLine(x1, x2, fmax, kGreen + 2, 1, 1); HLine(x1, x2, fmax - 1, kGreen + 2, 1, 1);
            gF->Draw("LP SAME");
            if (w >= 0) VLine(gC[w].thr, std::max(fmin, fmax - 4.) - 0.4, fmax + 0.8, kGreen + 2, 1, 2);
            Head("figure of merit; dashed: no cut; green: best, best - 1, working point", 0.17, 0.94, 0.040);
        }
        Save(c, kind == 0 ? "extra_scans_peak" : "extra_scans_bump", true);
    }
}
void StageSculptFig() {
    const char* names[6] = {"vz > -160", "isUsedVertex both pions", "no pi0 photon", "|z - vz| (both pi) < 40", "MM4 < 2.80", "pt4 > 300"};
    TCanvas* c = new TCanvas("c_sculpt", "", 2400, 1300);
    c->Divide(3, 2);
    const double lo = scan::kMlo, bw = scan::kBinW * scan::kCBreb;
    for (int j = 0; j < 6; ++j) {
        const int i = Idx(names[j]);
        TGraphErrors* g = GE(kBlack, 20, 1.0);
        double mn = 1, mx = 0;
        for (size_t b = 0; b < gR[i].pP.size(); ++b) {
            const double m = lo + (b + 0.5) * bw, n = gR[i].pP[b] + gR[i].pF[b];
            if (m < 0.40 || m > 0.70 || !(n > 0)) continue;
            const double f = gR[i].pP[b] / n;
            AddPt(g, m, f, std::sqrt(std::max(f * (1 - f), 1e-6) / n));
            mn = std::min(mn, f); mx = std::max(mx, f);
        }
        c->cd(j + 1);
        const double y1 = mn - 0.35 * (mx - mn) - 0.01, y2 = mx + 0.45 * (mx - mn) + 0.01;
        Frame(Form("fr_scu%d", j), 0.40, 0.70, y1, y2, kXt, "pass fraction of all pairs");
        g->Draw("P SAME");
        VLine(scan::kPeakLo, y1, mx, kBlue + 1, 2, 1); VLine(scan::kPeakHi, y1, mx, kBlue + 1, 2, 1);
        VLine(scan::kBumpLo, y1, mx, kGray + 2, 3, 1); VLine(scan::kBumpHi, y1, mx, kGray + 2, 3, 1);
        VLine(scan::kFitHi, y1, mx, kGray + 2, 3, 1);
        const PeakNum p = PeakOf(gR[i], 0, gBase[0]);
        Head(Form("cut: %s  (on top of the peak nominal)", names[j]), 0.17, 0.94, 0.042);
        Txt(0.20, 0.86, Form("#varepsilon_{S} = %.3f #pm %.3f, #varepsilon_{B} = %.3f", p.epsS.e, p.epsS.err, p.epsB), 0.040, kBlue + 1);
        if (j == 0) { Txt(0.20, 0.80, "blue dashed: peak window", 0.034, kBlue + 1);
                      Txt(0.20, 0.75, "grey dotted: bump window, end of the fit range", 0.034, kGray + 2); }
    }
    Save(c, "extra_sculpting", false);
}
void StageSbFig() {
    if (gTop3.size() < 3) return;
    TCanvas* c = new TCanvas("c_sb", "", 2400, 800);
    c->Divide(3, 1);
    for (int j = 0; j < 3; ++j) {
        const int iv = gTop3[j];
        const SbRes& s = gSb[iv];
        int nb = 14; double lo = 0, hi = 1;
        if (s.discrete) { nb = (int)s.edges.size() + 1; lo = s.edges.front() - 1.; hi = s.edges.back() + 1.; }   // values one unit apart
        else {
            std::vector<double> v;
            for (const auto& r : gP[0].s[0]) if (SbRegion(r.m4, 0) > 1) v.push_back(gV[iv].f(r));
            std::sort(v.begin(), v.end());
            lo = v[(size_t)(0.01 * v.size())]; hi = v[(size_t)(0.97 * v.size())];
            if (!(hi > lo)) hi = lo + 1;
        }
        TH1D *hs = nullptr, *hb = nullptr;
        SbHists(iv, gBase[0].Swin, nb, lo, hi, hs, hb);
        c->cd(j + 1);
        std::string xt = gV[iv].name;
        if (xt == "isUsedVertex pions") xt = "number of pions used in the vertex fit";
        if (xt == "fw1_p") xt = "p of the forward candidate [MeV/c] (first bin: none)";
        if (xt == "|z-vz| pion max") xt = "larger |z_{#pi} - z_{vertex}| of the two pions [mm]";
        Deco(hb, kBlack, 1, 0, xt.c_str(), "fraction of rows per bin"); Deco(hs, kBlue + 1, 20, 1.1);
        hb->SetLineWidth(2);
        double mx = 0, mn = 0;
        for (int b = 1; b <= nb; ++b) { mx = std::max(mx, std::max(hb->GetBinContent(b), hs->GetBinContent(b) + hs->GetBinError(b)));
                                        mn = std::min(mn, hs->GetBinContent(b) - hs->GetBinError(b)); }
        hb->SetMinimum(mn - 0.02); hb->SetMaximum(1.45 * mx);
        hb->Draw("HIST"); hs->Draw("E1 SAME");
        HLine(lo, hi, 0, kGray + 1, 2);
        Head(Form("%s: z = %.1f (controls %.1f, %.1f)", gV[iv].name.c_str(), s.z[0], s.z[1], s.z[2]), 0.17, 0.94, 0.042);
        TLegend* l = Leg(0.32, 0.76, 0.94, 0.90, 0.036);
        l->AddEntry(hs, "#eta peak: window - sidebands", "pe"); l->AddEntry(hb, "background: sidebands", "l"); l->Draw();
        if (!s.discrete) Txt(0.32, 0.71, "end bins hold under- and overflow", 0.032, kGray + 2);
    }
    Save(c, "extra_sbsub_top3", false);
}

} // namespace sx

void scan_extra() {
    using namespace sx;
    TStopwatch swAll;
    Style();
    gSystem->mkdir(scan::kFigDir, kTRUE);
    if (!Prepare()) return;
    BuildCandidates();
    gR = Evaluate(gC);
    for (int f = 0; f < 3; ++f) gBase[f] = gR[0].h[f][0];
    StageTables();
    StageN();
    StageRank();
    TFile* fout = TFile::Open(Form("%s/extra_scan.root", kTabDir), "RECREATE");
    StageVertex(fout);
    StageTagged(fout);
    StageCombine();
    StageBeforeAfter(fout);
    StageSummary();
    StageScans();
    StageSculptFig();
    StageSbFig();
    if (fout) { fout->Close(); delete fout; }
    CloseTable();
    printf("scan_extra: total %.1f s\n", swAll.RealTime());
}
