// ============================================================================
// s5_blind_prep.C  --  blind test of the signal extraction: halves, templates, pseudo-data and a sealed key
// ============================================================================
// HADES p+p, T = 4.53 GeV, pi+ pi- e+ e- (e+e- hypothesis). SIMULATION ONLY: no measured data are read (SimEngine.h
// includes ../ScanEngine.h for its types; no data function is called). No combinatorial background exists in the
// simulation and none is added.
//
// Simulation: SMASH std + lep + eta owned rows, partition mode A (SimOpt::sets = 7), weight w (PDG) x physics weights
// of the eta decays, main variant (simphys::Init(); simphys::Install(); SimOpt::phys = simphys::kPhysMain = 108 = V3).
// Simulated sum = template 1 (signal classes 2-5) + template 2 (physical background). Templates 3 (fake lepton) and
// 4 (different sources) only as variant templates; X17 (class 6, W = 2 BR per row) never in the sum.
// Base "pre" = SimEngine kPre. The frozen selections are read from tables/s5_opt_eta_frozen.txt ([PEAK], [DIRECT],
// [BUMP]) and tables/s5_opt_x17_frozen.txt and used unchanged (comparisons as in the optimisation macros: float
// thresholds for the eta cuts as in s5_opt_eta.C, float value against double threshold for X17 as in s5_opt_x17.C).
//
// 1. HALVES. Event group = maximal run of adjacent rows with equal (event, eVertReco_z) in one tree (column grp of the
//    compact files, numbered per tree in file order by sim_prepare.C). Half A = even grp, half B = odd grp (= column
//    'fold'; the same split as the split-half checks of s5_opt_eta.C and s5_opt_x17.C). All rows of a group go to the
//    same half. The three sets come from the same SMASH events; their cross-set correlation is neglected (as in
//    SimEngine).
// 2. TEMPLATES (truth allowed) from half A with the weights x 2, i.e. at the full luminosity of one set
//    (L_sim = 221.05 nb^-1): content = expected rows at L_sim, errors = event-group errors x 2.
//    -> blind/templates/templates.root and blind/templates/README.md (every histogram listed).
// 3. ETA PSEUDO-DATA at L_sim / 2 from half B: Poisson bootstrap per event group. Candidate rows: half B, owned (mode A),
//    pre, template 1 classes 2-5 or template 2, w x physics factor > 0, 0.30 <= M4 < 1.20. Group weight
//    w_g = largest row weight of the group's candidate rows; the group enters n ~ Poisson(w_g) times; in each copy one
//    uniform u is drawn and every row with w_r / w_g > u is kept ("coherent thinning": the rows of largest weight are
//    always kept, rows of the same weight always together, rows of smaller weight in a nested subset). Expected
//    multiplicity of every row = w_r exactly. Each copy becomes one event; the events of a toy are shuffled and
//    numbered. Written: ONLY toy, evt, m4, mee, mm4, dz_pi, pivtx (no truth, no weight, no class, no set).
//    -> blind/pseudo/eta_toys.root
// 4. X17 PSEUDO-DATA at L_data = 6080.5 nb^-1 (UNCONFIRMED), parametric: expectation = half B x 2 x L_data / L_sim in
//    M4 (0.4830-0.6080, 2.5 MeV/c2) x M(ee) (0-100 MeV/c2, 0.5 MeV/c2) after pre + the frozen X17 cuts except the M(ee)
//    window: background T1 (no X17) + T2, plus X17 (class 6, W = 2 BR) at BR in {0, 1e-6, 3e-6, 1e-5, 3e-5}; per cell
//    n = Poisson(mu_bkg) + Poisson(BR s1); 300 toys per BR, shuffled, random six-digit ids.
//    -> blind/pseudo/x17_toys.root
// 5. SEALED KEY (blind/sealed/ only; never printed to stdout, to the log or to the public tables): secret seed, true
//    counts per class (group) of every eta toy in every M4 region of every frozen selection, half-B expectations and
//    the A/B comparison, the injected BR and the expected and sampled counts of every X17 toy.
//
// Public outputs: figs/s5_blind_prep_{templates_eta,templates_x17,toys_check}.png/.pdf,
//   tables/s5_blind_prep_templates.txt (half-A template integrals), tables/s5_blind_prep_checks.txt.
// The branch check of blind/pseudo is s5_blind_prep_check.C.
//
// Run from /home/przygoda/HADES/FAT/FAT (single thread; about 10 GB, a few minutes):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/aclic_s5_blind_prep", kTRUE)' \
//        'plots/eta_sim/s5_blind_prep.C+(0, 200, 300)'
// arguments: seed (0 = a new secret seed from std::random_device, written only to blind/sealed/seed.txt; 1 = read the seed
// from blind/sealed/seed.txt and reproduce the same toys; any other value = that seed), number of eta toys (200), X17
// toys per BR value (300)
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <limits>
#include <random>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "TBox.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TNamed.h"
#include "TPad.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include "SimEngine.h"
#include "sim_eta_physics.h"

namespace s5b {
using namespace simeng;

// ============================================================================
// paths and constants
// ============================================================================
const char* const kTag = "s5_blind_prep";
const char* const kBlind = "/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/blind";
const char* const kEtaFrozen = "/home/przygoda/HADES/FAT/FAT/plots/eta_sim/tables/s5_opt_eta_frozen.txt";
const char* const kX17Frozen = "/home/przygoda/HADES/FAT/FAT/plots/eta_sim/tables/s5_opt_x17_frozen.txt";
std::string BDir(const char* sub) { return std::string(kBlind) + "/" + sub; }
const double kK = kLdata / kLsim;              // sum w -> expected rows at L_data (UNCONFIRMED)
const char* const kLdNote = "L_{data} = 6080.5 nb^{-1} (UNCONFIRMED)";
// M4 windows of this study (fixed for all agents; checked against the frozen file)
const double kPkLo = 0.5333, kPkHi = 0.5587, kBuLo = 0.4576, kBuHi = 0.5255;
double gXLo = 0.5330, gXHi = 0.5555;           // X17 frozen M4 window (read from the frozen file)
// templates and eta pseudo-data
const double kM4Lo = 0.30, kM4Hi = 1.20;
const int kNM4 = 900;                           // 1 MeV/c2
// X17 histograms
const double kXm4Lo = 0.4830, kXm4W = 0.0025;
const int kXm4N = 50;                           // 0.4830-0.6080, 2.5 MeV/c2
const double kXmeeW = 0.5;
const int kXmeeN = 200;                         // 0-100 MeV/c2, 0.5 MeV/c2
const int NBR = 5;
const double kBR[NBR] = {0., 1e-6, 3e-6, 1e-5, 3e-5};

// ============================================================================
// output helpers: public table (stdout + file) and sealed text (file only)
// ============================================================================
FILE* gT = nullptr;   // public table
FILE* gK = nullptr;   // sealed text: NEVER echoed to stdout
void PT(const char* fmt, ...) {
    char buf[65536];
    va_list ap; va_start(ap, fmt); vsnprintf(buf, sizeof buf, fmt, ap); va_end(ap);
    fputs(buf, stdout);
    if (gT) fputs(buf, gT);
}
void PK(const char* fmt, ...) {
    if (!gK) return;
    va_list ap; va_start(ap, fmt); vfprintf(gK, fmt, ap); va_end(ap);
}
int gCnt = 0;
std::string Nm(const char* s) { return Form("%s_%d", s, gCnt++); }

// ============================================================================
// selections
// ============================================================================
struct Cuts {
    std::string name;
    double oa_min = NAN, oa_max = NAN, lab_max = NAN, mpp_min = NAN, mpp_max = NAN, eta_min = NAN, dz_max = NAN,
           mm4_max = NAN, pt4_min = NAN, mee_lo = NAN, mee_hi = NAN;   // mee in MeV/c2
    bool vtx = false;
    bool fcmp = true;   // true: float thresholds (s5_opt_eta.C); false: float value vs double threshold (s5_opt_x17.C)
    // information from the frozen file (cross-check)
    double fom = NAN, purity = NAN, eff = NAN;
};
inline bool Has(double x) { return x == x; }
inline bool Lower(float v, double c, bool f) { return f ? (v > (float)c) : ((double)v > c); }
inline bool Upper(float v, double c, bool f) { return f ? (v < (float)c) : ((double)v < c); }
inline float DzPi(const Row& r) { return std::max(std::fabs(r.pip_z - r.vz), std::fabs(r.pim_z - r.vz)); }
inline bool PassCuts(const Row& r, const Cuts& c) {
    if (!PassPre(r)) return false;
    if (Has(c.oa_min) && !Lower(r.oa, c.oa_min, c.fcmp)) return false;
    if (Has(c.oa_max) && !Upper(r.oa, c.oa_max, c.fcmp)) return false;
    if (Has(c.lab_max) && !Upper(r.oa_lab, c.lab_max, c.fcmp)) return false;
    if (Has(c.mpp_min) && !Lower(r.mpipi, c.mpp_min, c.fcmp)) return false;
    if (Has(c.mpp_max) && !Upper(r.mpipi, c.mpp_max, c.fcmp)) return false;
    if (Has(c.eta_min) && !Lower(r.oa_eta, c.eta_min, c.fcmp)) return false;
    if (c.vtx && !(r.pip_isUsedVertex == 1 && r.pim_isUsedVertex == 1)) return false;
    if (Has(c.dz_max) && !Upper(DzPi(r), c.dz_max, c.fcmp)) return false;
    if (Has(c.mm4_max) && !Upper(r.mm4, c.mm4_max, c.fcmp)) return false;
    if (Has(c.pt4_min) && !Lower(r.pt4, c.pt4_min, c.fcmp)) return false;
    if (Has(c.mee_lo) || Has(c.mee_hi)) {
        const float m = r.mee * 1000.f;
        if (Has(c.mee_lo) && !((double)m >= c.mee_lo)) return false;
        if (Has(c.mee_hi) && !((double)m < c.mee_hi)) return false;
    }
    return true;
}
// compact TLatex label of a selection (figures)
std::string LatexCuts(const Cuts& c) {
    std::string s = "pre";
    auto add = [&](bool on, const std::string& t) { if (on) s += " + " + t; };
    add(Has(c.oa_min), Form("OA > %g#circ", c.oa_min));
    add(Has(c.oa_max), Form("OA < %g#circ", c.oa_max));
    add(Has(c.lab_max), Form("#alpha_{lab} < %g#circ", c.lab_max));
    add(Has(c.mpp_min), Form("M_{#pi#pi} > %g", c.mpp_min));
    add(Has(c.mpp_max), Form("M_{#pi#pi} < %g", c.mpp_max));
    add(Has(c.eta_min), Form("#alpha_{#eta} > %g#circ", c.eta_min));
    add(c.vtx, "#pi#pi in vertex fit");
    add(Has(c.dz_max), Form("|#Deltaz_{#pi}| < %g mm", c.dz_max));
    add(Has(c.mm4_max), Form("MM_{4} < %g GeV/c^{2}", c.mm4_max));
    add(Has(c.pt4_min), Form("p_{T,4} > %g MeV/c", c.pt4_min));
    if (Has(c.mee_lo) || Has(c.mee_hi)) s += Form(" + %g #leq M_{ee} < %g MeV/c^{2}", c.mee_lo, c.mee_hi);
    return s;
}
std::string CutStr(const Cuts& c) {
    std::string s = "pre";
    auto add = [&](bool on, const std::string& t) { if (on) s += " + " + t; };
    add(Has(c.oa_min), Form("OA > %g deg", c.oa_min));
    add(Has(c.oa_max), Form("OA < %g deg", c.oa_max));
    add(Has(c.lab_max), Form("alpha_lab < %g deg", c.lab_max));
    add(Has(c.mpp_min), Form("M(pipi) > %g GeV/c2", c.mpp_min));
    add(Has(c.mpp_max), Form("M(pipi) < %g GeV/c2", c.mpp_max));
    add(Has(c.eta_min), Form("alpha_eta > %g deg", c.eta_min));
    add(c.vtx, "both pions in the vertex fit");
    add(Has(c.dz_max), Form("max|z_pi - z_vtx| < %g mm", c.dz_max));
    add(Has(c.mm4_max), Form("MM4 < %g GeV/c2", c.mm4_max));
    add(Has(c.pt4_min), Form("pT4 > %g MeV/c", c.pt4_min));
    add(Has(c.mee_lo), Form("M(ee) >= %g MeV/c2", c.mee_lo));
    add(Has(c.mee_hi), Form("M(ee) < %g MeV/c2", c.mee_hi));
    return s;
}

// selections of this macro: the three frozen eta selections, the frozen X17 selection (with its M(ee) window; the M4
// window is applied separately), and "pre" (efficiency denominators only)
enum { sPEAK = 0, sDIRECT, sBUMP, sX17, sPRE, NS };
const int NSW = 4;   // selections written (templates, pseudo-data key)
const char* kSN[NS] = {"PEAK", "DIRECT", "BUMP", "X17SEL", "PRE"};
Cuts gCut[NS];
Cuts gCutX17nomee;   // X17 frozen cuts without the M(ee) window (M(ee) templates and the X17 pseudo-data)

// ---------------------------------------------------------------- frozen files
struct KV { std::string sec, key, val; };
std::string Trim(const std::string& s) {
    const size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    const size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}
std::vector<KV> ReadKV(const char* path, bool& ok) {
    std::vector<KV> out;
    std::ifstream in(path);
    ok = in.good();
    std::string line, sec;
    while (std::getline(in, line)) {
        const size_t h = line.find('#');
        if (h != std::string::npos) line = line.substr(0, h);
        line = Trim(line);
        if (line.empty()) continue;
        if (line[0] == '[') { sec = line.substr(1, line.find(']') - 1); continue; }
        const size_t e = line.find('=');
        if (e == std::string::npos) continue;
        out.push_back({sec, Trim(line.substr(0, e)), Trim(line.substr(e + 1))});
    }
    return out;
}
double Num(const std::string& v) { return (v == "none" || v.empty()) ? NAN : atof(v.c_str()); }
bool Flag(const std::string& v) { return v == "1" || v == "yes" || v == "on" || v == "true"; }
bool EtaCuts(const std::vector<KV>& kv, const std::string& sec, Cuts& c) {
    c = Cuts(); c.name = sec; c.fcmp = true;
    int n = 0;
    for (const KV& x : kv) {
        if (x.sec != sec) continue;
        ++n;
        if (x.key == "oa_min") c.oa_min = Num(x.val);
        else if (x.key == "alpha_lab_max") c.lab_max = Num(x.val);
        else if (x.key == "mpipi_min") c.mpp_min = Num(x.val);
        else if (x.key == "mpipi_max") c.mpp_max = Num(x.val);
        else if (x.key == "alpha_eta_min") c.eta_min = Num(x.val);
        else if (x.key == "pions_in_vertex_fit") c.vtx = Flag(x.val);
        else if (x.key == "dz_pi_max") c.dz_max = Num(x.val);
        else if (x.key == "mm4_max") c.mm4_max = Num(x.val);
        else if (x.key == "pt4_min") c.pt4_min = Num(x.val);
        else if (x.key == "fom") c.fom = Num(x.val);
        else if (x.key == "purity") c.purity = Num(x.val);
        else if (x.key == "efficiency") c.eff = Num(x.val);
        else if (x.key == "target" || x.key == "fom_over_best" || x.key == "best_point" ||
                 x.key == "optimisation_bias_frozen_rule" || x.key == "justification") {}
        else { printf("s5b: unknown key '%s' in [%s] of %s: stop\n", x.key.c_str(), sec.c_str(), kEtaFrozen); return false; }
    }
    if (n == 0) { printf("s5b: section [%s] missing in %s\n", sec.c_str(), kEtaFrozen); return false; }
    return true;
}
double gX17sigFrozen = NAN, gX17bkgFrozen = NAN, gX17bkgFrozenErr = NAN;
bool X17Cuts(const std::vector<KV>& kv, Cuts& c) {
    c = Cuts(); c.name = "X17SEL"; c.fcmp = false;
    for (const KV& x : kv) {
        if (x.key == "m4_lo") gXLo = Num(x.val);
        else if (x.key == "m4_hi") gXHi = Num(x.val);
        else if (x.key == "mee_lo_MeV") c.mee_lo = Num(x.val);
        else if (x.key == "mee_hi_MeV") c.mee_hi = Num(x.val);
        else if (x.key == "oa_lo_deg") c.oa_min = Num(x.val);
        else if (x.key == "oa_hi_deg") c.oa_max = Num(x.val);
        else if (x.key == "alpha_lab_hi_deg") c.lab_max = Num(x.val);
        else if (x.key == "mpipi_lo") c.mpp_min = Num(x.val);
        else if (x.key == "mpipi_hi") c.mpp_max = Num(x.val);
        else if (x.key == "alpha_eta_lo_deg") c.eta_min = Num(x.val);
        else if (x.key == "pions_in_vertex_fit") c.vtx = Flag(x.val);
        else if (x.key == "dz_pi_hi_mm") c.dz_max = Num(x.val);
        else if (x.key == "mm4_hi") c.mm4_max = Num(x.val);
        else if (x.key == "pt4_lo") c.pt4_min = Num(x.val);
        else if (x.key == "x17_sigma_per_unit_BR_nb") gX17sigFrozen = Num(x.val);
        else if (x.key == "background_sigma_nb") {
            gX17bkgFrozen = Num(x.val);
            const size_t p = x.val.find("+-");
            if (p != std::string::npos) gX17bkgFrozenErr = atof(x.val.substr(p + 2).c_str());
        }
        else if (x.key == "base" || x.key == "x17_rows_per_unit_BR" || x.key == "background_rows" ||
                 x.key.rfind("expected_", 0) == 0 || x.key.rfind("Z_", 0) == 0) {}
        else { printf("s5b: unknown key '%s' in %s: stop\n", x.key.c_str(), kX17Frozen); return false; }
    }
    return Has(gXLo) && Has(gXHi);
}

// ============================================================================
// components and groups
// ============================================================================
enum { cPeak45 = 0, cDir5, cConv4, cBump23, cBump2, cBump3, cBkg, cSum, cGeta, cGpi0D, cGconv, cGvec, cGoth, cFake, cDiff,
       cX17, cOther, NC };
const char* kCK[NC] = {"peak45", "direct5", "conv4", "bump23", "bump2", "bump3", "bkgT2", "sumT12", "Geta", "Gpi0D",
                       "Gconv", "Gvec", "Goth", "fakeT3", "diffT4", "x17perBR", "otherX17"};
const char* kCD[NC] = {
    "eta peak = signal classes 4 + 5 (conversion peak, std set + direct decay, eta set)",
    "direct eta -> pi+pi-e+e- (class 5, eta set; physics weight V3)",
    "conversion peak eta -> pi+pi-gamma, gamma -> e+e- (class 4, std set)",
    "bump = classes 2 + 3 (eta -> pi+pi-pi0; pi0 Dalitz or pi0 photon converted; photon lost)",
    "bump, class 2 (pi0 -> e+e-gamma)",
    "bump, class 3 (photon of the pi0 converted)",
    "physical background = template 2 (true e+e- pair from one source, not a signal combination); with peak45 and bump23 it "
    "makes the simulated sum",
    "simulated sum = template 1 (classes 2-5) + template 2 = peak45 + bump23 + bkgT2",
    "G_eta: template-2 rows with both pions from one eta (pion class P2) or an eta Dalitz pair (L2)",
    "G_pi0D: remaining template-2 rows with a pi0 Dalitz pair (L1)",
    "G_conv: remaining template-2 rows with a conversion pair (L3)",
    "G_vec: remaining template-2 rows with rho0/omega/phi -> e+e-, omega -> pi0 e+e- or a second-decay pair (L4-L6)",
    "G_oth: all other template-2 rows",
    "template 3, fake lepton (VARIANT only; not in sumT12, not in the pseudo-data)",
    "template 4, leptons from two different sources (VARIANT only; not in sumT12, not in the pseudo-data)",
    "X17 (class 6, eta set, eta -> pi+pi- X17, X17 -> e+e-, m = 17 MeV/c2) PER UNIT BR: W = 2 per row; multiply by BR; never in sumT12",
    "X17 background 'other' = sumT12 - direct5 - conv4 = bump23 + bkgT2"};
const int NG = 9;   // key groups: Gsig2, Gsig3, Gsig4, Gsig5, G_eta, G_pi0D, G_conv, G_vec, G_oth
const char* kGN[NG] = {"Gsig2", "Gsig3", "Gsig4", "Gsig5", "G_eta", "G_pi0D", "G_conv", "G_vec", "G_oth"};
inline int BkgGroup(const Row& r) {   // template 2 rows
    if (r.pcls == 2 || r.lcls == 2) return 4;
    if (r.lcls == 1) return 5;
    if (r.lcls == 3) return 6;
    if (r.lcls >= 4 && r.lcls <= 6) return 7;
    return 8;
}
inline int GroupOf(const Row& r) {
    if (r.tmpl == 1) return (r.cls >= 2 && r.cls <= 5) ? r.cls - 2 : -1;
    if (r.tmpl == 2) return BkgGroup(r);
    return -1;
}
std::vector<Comp> MakeComps() {
    std::vector<Comp> v(NC);
    v[cPeak45] = Comp(kCK[cPeak45], [](const Row& r) { return r.tmpl == 1 && (r.cls == 4 || r.cls == 5); });
    v[cDir5] = Comp(kCK[cDir5], [](const Row& r) { return r.tmpl == 1 && r.cls == 5; });
    v[cConv4] = Comp(kCK[cConv4], [](const Row& r) { return r.tmpl == 1 && r.cls == 4; });
    v[cBump23] = Comp(kCK[cBump23], [](const Row& r) { return r.tmpl == 1 && (r.cls == 2 || r.cls == 3); });
    v[cBump2] = Comp(kCK[cBump2], [](const Row& r) { return r.tmpl == 1 && r.cls == 2; });
    v[cBump3] = Comp(kCK[cBump3], [](const Row& r) { return r.tmpl == 1 && r.cls == 3; });
    v[cBkg] = Comp(kCK[cBkg], [](const Row& r) { return r.tmpl == 2; });
    v[cSum] = Comp(kCK[cSum], [](const Row& r) { return (r.tmpl == 1 && r.cls >= 2 && r.cls <= 5) || r.tmpl == 2; });
    v[cGeta] = Comp(kCK[cGeta], [](const Row& r) { return r.tmpl == 2 && BkgGroup(r) == 4; });
    v[cGpi0D] = Comp(kCK[cGpi0D], [](const Row& r) { return r.tmpl == 2 && BkgGroup(r) == 5; });
    v[cGconv] = Comp(kCK[cGconv], [](const Row& r) { return r.tmpl == 2 && BkgGroup(r) == 6; });
    v[cGvec] = Comp(kCK[cGvec], [](const Row& r) { return r.tmpl == 2 && BkgGroup(r) == 7; });
    v[cGoth] = Comp(kCK[cGoth], [](const Row& r) { return r.tmpl == 2 && BkgGroup(r) == 8; });
    v[cFake] = Comp(kCK[cFake], [](const Row& r) { return r.tmpl == 3; });
    v[cDiff] = Comp(kCK[cDiff], [](const Row& r) { return r.tmpl == 4; });
    v[cX17] = Comp(kCK[cX17], [](const Row& r) { return r.tmpl == 1 && r.cls == 6; }, [](const Row&) { return 2.; });
    v[cOther] = Comp(kCK[cOther], [](const Row& r) { return (r.tmpl == 1 && (r.cls == 2 || r.cls == 3)) || r.tmpl == 2; });
    return v;
}

// ============================================================================
// regions and binnings
// ============================================================================
const int NR = 7;            // disjoint M4 regions; index 7 = X17 M4 window (overlaps 2 and 3)
const int NRX = NR + 1;
const double kRegEdge[NR + 1] = {0.30, kBuLo, kBuHi, kPkLo, kPkHi, 0.70, 1.00, 1.20};
std::string RegName(int k) {
    static const char* n[NR] = {"0.30-0.4576 (below bump)", "bump window 0.4576-0.5255", "gap 0.5255-0.5333",
                                "peak window 0.5333-0.5587", "0.5587-0.70", "0.70-1.00", "1.00-1.20 (above 1.00)"};
    if (k < NR) return n[k];
    return Form("X17 M4 window %.4f-%.4f", gXLo, gXHi);
}
std::string RegShort(int k) {
    static const char* n[NRX] = {"below", "bump", "gap", "peak", "056_070", "070_100", "100_120", "x17win"};
    return n[k];
}
inline int RegionOf(double m4) {
    for (int i = 0; i < NR; ++i) if (m4 >= kRegEdge[i] && m4 < kRegEdge[i + 1]) return i;
    return -1;
}
inline bool InX17Win(double m4) { return m4 >= gXLo && m4 < gXHi; }
inline int XM4Bin(double m4) {   // 0-based bin of the X17 M4 axis, -1 outside
    const double u = (m4 - kXm4Lo) / kXm4W;
    if (!(u >= 0) || u >= kXm4N) return -1;
    return std::min(kXm4N - 1, (int)u);
}
inline int XMeeBin(float meeMeV) {
    const double u = (double)meeMeV / kXmeeW;
    if (!(u >= 0) || u >= kXmeeN) return -1;
    return std::min(kXmeeN - 1, (int)u);
}
int gXwBin0 = -1, gXwBin1 = -1;   // 0-based M4 bins of the X17 window (inclusive)
int gMwBin0 = -1, gMwBin1 = -1;   // 0-based M(ee) bins of the frozen M(ee) window (inclusive)

Binning BM4() { Binning b; b.nbin = kNM4; b.lo = kM4Lo; b.hi = kM4Hi; return b; }
Binning BReg() {
    Binning b; b.nbin = NR; b.lo = 0; b.hi = NR;
    b.x = [](const Row& r) { const int k = RegionOf(r.m4); return k < 0 ? -1. : k + 0.5; };
    return b;
}
Binning BWin() {
    Binning b; b.nbin = 1; b.lo = 0; b.hi = 1;
    b.x = [](const Row& r) { return InX17Win(r.m4) ? 0.5 : -1.; };
    return b;
}
Binning BMee() {   // M(ee) [MeV/c2] of the rows in the X17 M4 window
    Binning b; b.nbin = kXmeeN; b.lo = 0; b.hi = kXmeeN * kXmeeW;
    b.x = [](const Row& r) { return InX17Win(r.m4) ? (double)(r.mee * 1000.f) : -1.; };
    return b;
}
Binning B2D() {    // flattened (M4 bin, M(ee) bin)
    Binning b; b.nbin = kXm4N * kXmeeN; b.lo = 0; b.hi = b.nbin;
    b.x = [](const Row& r) {
        const int ix = XM4Bin(r.m4), iy = XMeeBin(r.mee * 1000.f);
        return (ix < 0 || iy < 0) ? -1. : ix * kXmeeN + iy + 0.5;
    };
    return b;
}
// half: 0 = A (even grp), 1 = B (odd grp), -1 = all
inline bool InHalf(const Row& r, int half) { return half < 0 || (r.grp & 1) == half; }
Sel SelOf(const Cuts& c, int half) { return [c, half](const Row& r) { return InHalf(r, half) && PassCuts(r, c); }; }

// ============================================================================
// style and figures helpers
// ============================================================================
void Style() {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetEndErrorSize(0);
    gStyle->SetTextFont(42);
    gStyle->SetLabelFont(42, "xyz");
    gStyle->SetTitleFont(42, "xyz");
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetLineWidth(2);
    gStyle->SetHistLineWidth(3);
}
// vector PDF with a tight page (EPS -> ps2pdf -dEPSCrop) and a PNG rasterised from that PDF
void Save(TCanvas* c, const char* name) {
    const std::string base = std::string(kFigDir) + "/" + name;
    c->SaveAs((base + ".eps").c_str());
    gSystem->Exec(Form("cd %s && ps2pdf -dEPSCrop %s.eps %s.pdf && rm -f %s.eps && pdftoppm -png -singlefile -scale-to-x %d "
                       "-scale-to-y -1 %s.pdf %s",
                       kFigDir, name, name, name, c->GetWw(), name, name));
}
void Txt(double x, double y, double size, const char* s, int align = 11, int col = 1) {
    TLatex l;
    l.SetNDC();
    l.SetTextFont(42);
    l.SetTextSize(size);
    l.SetTextAlign(align);
    l.SetTextColor(col);
    l.DrawLatex(x, y, s);
}
void AxisStyle(TH1* h, double ls, double ts, double xo = 1.0, double yo = 1.0) {
    for (TAxis* a : {h->GetXaxis(), h->GetYaxis()}) {
        a->SetLabelSize(ls); a->SetTitleSize(ts); a->SetLabelFont(42); a->SetTitleFont(42);
    }
    h->GetXaxis()->SetTitleOffset(xo);
    h->GetYaxis()->SetTitleOffset(yo);
}

// ============================================================================
// results kept between the parts
// ============================================================================
SimData* gD = nullptr;
SimOpt gO;
std::vector<Comp> gComps;
std::vector<WSpec> gReg[3][NS], gWin[3][NS];   // [0 A, 1 B, 2 all][selection][component]: regions and X17 window
std::vector<WSpec> gM4A[NSW];                  // M4 templates, half A (x 2)
std::vector<WSpec> gMeeH[2], g2DH[2];          // X17: M(ee) in the window and 2D, half A / half B (x 2)
// written template histograms (README)
std::vector<std::pair<std::string, std::string>> gReadme;   // name, description
std::vector<TH1*> gTplWrite;
// eta toys: observable totals per toy (all groups) for the checks
struct EtaTot { std::vector<double> n[NSW][NRX]; };
EtaTot gEtaTot;
double gMuTot[NSW][NRX], gVcoh[NSW][NRX], gVind[NSW][NRX];
long long gEtaRowsTotal = 0, gEtaEvtsTotal = 0;
int gNEta = 0, gNX17 = 0;
std::vector<double> gX17obsWin;   // observed rows in the frozen window per X17 toy (observable)

// ============================================================================
// part 1: frozen selections, regions, templates
// ============================================================================
bool SetupCuts() {
    bool ok1 = false, ok2 = false;
    const std::vector<KV> ke = ReadKV(kEtaFrozen, ok1);
    const std::vector<KV> kx = ReadKV(kX17Frozen, ok2);
    if (!ok1 || !ok2) { printf("s5b: cannot read the frozen files\n"); return false; }
    // windows of the frozen eta file must be the windows of this study
    for (const KV& x : ke) {
        if (x.sec != "windows") continue;
        const double v = Num(x.val);
        const double want = x.key == "peak_m4_min" ? kPkLo : x.key == "peak_m4_max" ? kPkHi : x.key == "bump_m4_min" ? kBuLo
                            : x.key == "bump_m4_max" ? kBuHi : NAN;
        if (!Has(want)) continue;
        if (!(std::fabs(v - want) < 1e-9)) { printf("s5b: window %s = %g differs from %g: stop\n", x.key.c_str(), v, want); return false; }
    }
    if (!EtaCuts(ke, "PEAK", gCut[sPEAK]) || !EtaCuts(ke, "DIRECT", gCut[sDIRECT]) || !EtaCuts(ke, "BUMP", gCut[sBUMP])) return false;
    if (!X17Cuts(kx, gCut[sX17])) return false;
    gCut[sPRE] = Cuts(); gCut[sPRE].name = "PRE";
    gCutX17nomee = gCut[sX17]; gCutX17nomee.mee_lo = NAN; gCutX17nomee.mee_hi = NAN; gCutX17nomee.name = "X17 without M(ee) window";
    // window bins of the X17 histograms
    gXwBin0 = (int)std::lround((gXLo - kXm4Lo) / kXm4W);
    gXwBin1 = (int)std::lround((gXHi - kXm4Lo) / kXm4W) - 1;
    gMwBin0 = (int)std::lround(gCut[sX17].mee_lo / kXmeeW);
    gMwBin1 = (int)std::lround(gCut[sX17].mee_hi / kXmeeW) - 1;
    if (std::fabs(kXm4Lo + gXwBin0 * kXm4W - gXLo) > 1e-9 || std::fabs(kXm4Lo + (gXwBin1 + 1) * kXm4W - gXHi) > 1e-9 ||
        std::fabs(gMwBin0 * kXmeeW - gCut[sX17].mee_lo) > 1e-9 || std::fabs((gMwBin1 + 1) * kXmeeW - gCut[sX17].mee_hi) > 1e-9) {
        printf("s5b: the X17 windows are not on the bin edges of the X17 histograms: stop\n");
        return false;
    }
    return true;
}

void PartRegions() {
    TStopwatch sw;
    for (int h = 0; h < 3; ++h)
        for (int s = 0; s < NS; ++s) {
            const int half = h == 2 ? -1 : h;
            gReg[h][s] = FillMany(*gD, SelOf(gCut[s], half), gComps, gO, BReg());
            gWin[h][s] = FillMany(*gD, SelOf(gCut[s], half), gComps, gO, BWin());
            if (h < 2) for (WSpec& w : gReg[h][s]) w.Scale(2.);
            if (h < 2) for (WSpec& w : gWin[h][s]) w.Scale(2.);
        }
    printf("regions filled in %.0f s\n", sw.RealTime());
}
// value, error of a region (k < NR) or the X17 window (k == NR)
inline double RV(int h, int s, int c, int k) { return k < NR ? gReg[h][s][c].s[k] : gWin[h][s][c].s[0]; }
inline double RE(int h, int s, int c, int k) { return std::sqrt(k < NR ? gReg[h][s][c].v[k] : gWin[h][s][c].v[0]); }

TH1D* RegHist(int h, int s, int c, const char* name) {
    TH1D* x = new TH1D(name, "", NRX, 0, NRX);
    x->SetDirectory(nullptr);
    x->Sumw2();
    for (int k = 0; k < NRX; ++k) {
        x->SetBinContent(k + 1, RV(h, s, c, k));
        x->SetBinError(k + 1, RE(h, s, c, k));
        x->GetXaxis()->SetBinLabel(k + 1, RegShort(k).c_str());
    }
    x->GetYaxis()->SetTitle("expected rows at L_{sim} = 221.05 nb^{-1}");
    return x;
}
void AddTpl(TH1* h, const std::string& desc) {
    h->SetTitle(desc.c_str());
    gTplWrite.push_back(h);
    gReadme.push_back(std::make_pair(std::string(h->GetName()), desc));
}

void PartTemplates() {
    TStopwatch sw;
    // M4 templates of the written selections (half A, x 2)
    for (int s = 0; s < NSW; ++s) {
        gM4A[s] = FillMany(*gD, SelOf(gCut[s], 0), gComps, gO, BM4());
        for (WSpec& w : gM4A[s]) w.Scale(2.);
    }
    // X17: M(ee) in the window and 2D, halves A and B
    for (int h = 0; h < 2; ++h) {
        gMeeH[h] = FillMany(*gD, SelOf(gCutX17nomee, h), gComps, gO, BMee());
        g2DH[h] = FillMany(*gD, SelOf(gCutX17nomee, h), gComps, gO, B2D());
        for (WSpec& w : gMeeH[h]) w.Scale(2.);
        for (WSpec& w : g2DH[h]) w.Scale(2.);
    }
    printf("templates filled in %.0f s\n", sw.RealTime());

    // ---------------------------------------------------------------- histograms
    for (int s = 0; s < NSW; ++s) {
        for (int c = 0; c < NC; ++c) {
            TH1D* h = gM4A[s][c].Hist(Form("m4_%s_%s", kSN[s], kCK[c]), 0);
            h->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]");
            h->GetYaxis()->SetTitle("expected rows / 1 MeV/c^{2} at L_{sim} = 221.05 nb^{-1}");
            AddTpl(h, std::string("M4 template, selection ") + kSN[s] + " (" + CutStr(gCut[s]) + "), half A x 2: " + kCD[c]);
        }
        for (int c = 0; c < NC; ++c) {
            TH1D* h = RegHist(0, s, c, Form("reg_%s_%s", kSN[s], kCK[c]));
            AddTpl(h, std::string("exact sums per M4 region (bins: below, bump, gap, peak, 056_070, 070_100, 100_120, x17win), "
                                  "selection ") + kSN[s] + ", half A x 2: " + kCK[c]);
        }
    }
    for (int c = 0; c < NC; ++c) {
        TH1D* h = gMeeH[0][c].Hist(Form("mee_X17_%s", kCK[c]), 0);
        h->GetXaxis()->SetTitle("M(e^{+}e^{-}) [MeV/c^{2}]");
        h->GetYaxis()->SetTitle("expected rows / 0.5 MeV/c^{2} at L_{sim} = 221.05 nb^{-1}");
        AddTpl(h, std::string(Form("M(ee) template in the X17 M4 window %.4f-%.4f GeV/c2, ", gXLo, gXHi)) + CutStr(gCutX17nomee) +
                      ", half A x 2: " + kCD[c]);
    }
    for (int c = 0; c < NC; ++c) {
        TH2D* h = new TH2D(Form("m4mee_X17_%s", kCK[c]), "", kXm4N, kXm4Lo, kXm4Lo + kXm4N * kXm4W, kXmeeN, 0., kXmeeN * kXmeeW);
        h->SetDirectory(nullptr);
        h->Sumw2();
        for (int ix = 0; ix < kXm4N; ++ix)
            for (int iy = 0; iy < kXmeeN; ++iy) {
                const int k = ix * kXmeeN + iy;
                h->SetBinContent(ix + 1, iy + 1, g2DH[0][c].s[k]);
                h->SetBinError(ix + 1, iy + 1, std::sqrt(g2DH[0][c].v[k]));
            }
        h->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]");
        h->GetYaxis()->SetTitle("M(e^{+}e^{-}) [MeV/c^{2}]");
        h->GetZaxis()->SetTitle("expected rows per cell at L_{sim} = 221.05 nb^{-1}");
        AddTpl(h, std::string("M4 x M(ee) template (2.5 MeV/c2 x 0.5 MeV/c2), ") + CutStr(gCutX17nomee) + ", half A x 2: " + kCD[c]);
    }
    // X17 line shapes (unit area)
    {
        TH1D* h = gMeeH[0][cX17].Hist("shape_X17_mee", 0);
        const double a = h->Integral();
        if (a > 0) h->Scale(1. / a);
        h->GetXaxis()->SetTitle("M(e^{+}e^{-}) [MeV/c^{2}]");
        h->GetYaxis()->SetTitle("fraction / 0.5 MeV/c^{2}");
        AddTpl(h, Form("X17 line shape in M(ee), unit area over 0-100 MeV/c2, X17 M4 window, %s, half A", CutStr(gCutX17nomee).c_str()));
        TH1D* g = gM4A[sX17][cX17].Hist("shape_X17_m4", 0);
        const double b = g->Integral();
        if (b > 0) g->Scale(1. / b);
        g->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]");
        g->GetYaxis()->SetTitle("fraction / 1 MeV/c^{2}");
        AddTpl(g, Form("X17 line shape in M4, unit area over 0.30-1.20 GeV/c2, selection X17SEL (%s), half A", CutStr(gCut[sX17]).c_str()));
    }
    // ---------------------------------------------------------------- write the file
    gSystem->mkdir(BDir("templates").c_str(), kTRUE);
    const std::string fn = BDir("templates") + "/templates.root";
    TFile* f = TFile::Open(fn.c_str(), "RECREATE");
    for (TH1* h : gTplWrite) h->Write();
    std::string info = std::string("Templates of the blind test (s5_blind_prep.C). Half A = event groups with even grp; weights x 2 = "
                                   "full luminosity of one set, L_sim = 221.05 nb^-1. Content = expected rows at L_sim (sum of w x "
                                   "physics weight V3, std + lep + eta owned rows, mode A), errors = event-group errors x 2. ") +
                       "Selections: PEAK = " + CutStr(gCut[sPEAK]) + "; DIRECT = " + CutStr(gCut[sDIRECT]) + "; BUMP = " +
                       CutStr(gCut[sBUMP]) + "; X17SEL = " + CutStr(gCut[sX17]) + Form("; X17 M4 window %.4f-%.4f GeV/c2.", gXLo, gXHi);
    TNamed("info", info.c_str()).Write();
    f->Close();
    delete f;
    printf("wrote %s (%zu histograms)\n", fn.c_str(), gTplWrite.size());

    // ---------------------------------------------------------------- README
    const std::string rn = BDir("templates") + "/README.md";
    FILE* r = fopen(rn.c_str(), "w");
    fprintf(r, "# Templates of the blind test (`templates.root`)\n\n");
    fprintf(r, "Written by `plots/eta_sim/s5_blind_prep.C` (simulation only, %s).\n\n", "SMASH std + lep + eta sets");
    fprintf(r, "**Source.** Half A of the simulation = the event groups with an EVEN group number `grp` in each e+e- tree (an event "
               "group is a run of adjacent rows with equal (event, eVertReco_z)). The pseudo-data come from the other half "
               "(odd `grp`), so templates and pseudo-data are statistically independent (apart from the shared SMASH events "
               "of the three sets, neglected as everywhere in this study).\n\n");
    fprintf(r, "**Normalisation.** Weights of half A x 2 = the full luminosity of one set, L_sim = 221.05 nb^-1. Every content is "
               "an EXPECTED NUMBER OF ROWS AT L_sim (sum of w (PDG) x physics weight of the eta decays, main variant V3, rows owned in "
               "partition mode A). Cross section = content / 221.05 nb^-1. For the eta pseudo-data (L_sim / 2) multiply by 0.5; "
               "for the X17 pseudo-data (L_data = 6080.5 nb^-1, UNCONFIRMED) multiply by 6080.5 / 221.05 = %.4f. Errors = "
               "event-group errors of half A x 2 (rows of one event correlated). Double precision (TH1D/TH2D).\n\n", kK);
    fprintf(r, "**Rows.** Base selection pre = isBest >= 0, vertex z > -500 mm, PT3 bit, vertex z > -160 mm, cor_ok == 1; all cut "
               "variables from corrected momenta. Simulated sum = template 1 (signal classes 2-5) + template 2 (physical "
               "background). Templates 3 and 4 are given as variants only (not in the sum, not in the pseudo-data). X17 never in the sum.\n\n");
    fprintf(r, "**Selections** (frozen files `plots/eta_sim/tables/s5_opt_eta_frozen.txt` and `s5_opt_x17_frozen.txt`, unchanged):\n\n");
    for (int s = 0; s < NSW; ++s) fprintf(r, "- `%s` = %s\n", kSN[s], CutStr(gCut[s]).c_str());
    fprintf(r, "- X17 M4 window %.4f <= M4 < %.4f GeV/c2 (M4 bins %d-%d of the 2D histograms, 1-based); M(ee) window "
               "%.1f <= M(ee) < %.1f MeV/c2 (M(ee) bins %d-%d, 1-based)\n", gXLo, gXHi, gXwBin0 + 1, gXwBin1 + 1, gCut[sX17].mee_lo,
            gCut[sX17].mee_hi, gMwBin0 + 1, gMwBin1 + 1);
    fprintf(r, "- M4 windows of the study: bump %.4f-%.4f, peak %.4f-%.4f GeV/c2 (lower edge included, upper excluded)\n\n", kBuLo,
            kBuHi, kPkLo, kPkHi);
    fprintf(r, "The cut values are compared as in the optimisation macros: eta cuts as float thresholds (`mm4 < 2.85f`, "
               "`dz_pi < 40.f`), X17 cuts as float value against double threshold (`mee*1000.f >= 13.0`).\n\n");
    fprintf(r, "**Components** (suffix of the names):\n\n| suffix | meaning |\n|---|---|\n");
    for (int c = 0; c < NC; ++c) fprintf(r, "| `%s` | %s |\n", kCK[c], kCD[c]);
    fprintf(r, "\n**Binning.** `m4_*`: M(pi+pi-e+e-) 0.30-1.20 GeV/c2, 900 bins of 1 MeV/c2 (rebin by 2 for 2 MeV/c2); the window "
               "edges are not on bin edges: use the `reg_*` histograms for exact window sums. `reg_*`: 8 bins = M4 regions "
               "below (0.30-0.4576), bump (0.4576-0.5255), gap (0.5255-0.5333), peak (0.5333-0.5587), 056_070 (0.5587-0.70), "
               "070_100 (0.70-1.00), 100_120 (1.00-1.20), x17win (%.4f-%.4f, overlaps gap and peak); exact sums of the rows, "
               "group errors. `mee_X17_*`: M(ee) 0-100 MeV/c2, 200 bins of 0.5 MeV/c2, rows in the X17 M4 window. `m4mee_X17_*`: "
               "M4 0.4830-0.6080 GeV/c2 (50 bins of 2.5 MeV/c2) x M(ee) 0-100 MeV/c2 (200 bins of 0.5 MeV/c2); the same binning "
               "as the X17 pseudo-data. `shape_*`: unit-area line shapes.\n\n", gXLo, gXHi);
    fprintf(r, "## Every histogram\n\n| name | meaning |\n|---|---|\n");
    for (const auto& x : gReadme) fprintf(r, "| `%s` | %s |\n", x.first.c_str(), x.second.c_str());
    fprintf(r, "| `info` (TNamed) | short description of the file |\n");
    fclose(r);
    printf("wrote %s\n", rn.c_str());
}

// public table of the template integrals and the frozen cross-check
void TablesPart1() {
    gT = fopen((std::string(kTabDir) + "/s5_blind_prep_templates.txt").c_str(), "w");
    PT("# s5_blind_prep_templates: integrals of the blind-test TEMPLATES (half A x 2 = full luminosity of one set), s5_blind_prep.C\n");
    PT("# simulation only; std + lep + eta owned rows, mode A, w (PDG) x physics weights V3 (kPhysMain); expected rows at L_sim = 221.05 nb^-1\n");
    PT("# (= sum w x 2 of half A); error = event-group error x 2; N_eff = S^2 / V. X17 row: per unit BR (W = 2 per row).\n");
    PT("# Rows at L_data = 6080.5 nb^-1 (UNCONFIRMED): multiply by %.4f; at L_sim / 2 (eta pseudo-data): multiply by 0.5.\n", kK);
    PT("# The half-B (pseudo-data) numbers are SEALED and not in this table.\n\n");
    for (int s = 0; s < NSW; ++s) {
        PT("## selection %s = %s\n", kSN[s], CutStr(gCut[s]).c_str());
        PT("%-10s", "component");
        for (int k = 0; k < NRX; ++k) PT(" %24s", RegShort(k).c_str());
        PT("\n");
        for (int c = 0; c < NC; ++c) {
            PT("%-10s", kCK[c]);
            for (int k = 0; k < NRX; ++k) {
                const double v = RV(0, s, c, k), e = RE(0, s, c, k);
                const double ne = e > 0 ? v * v / (e * e) : 0.;
                PT(" %11.4g+-%-8.3g(%5.0f)", v, e, ne);
            }
            PT("\n");
        }
        PT("\n");
    }
    PT("## X17 M(ee) templates in the X17 M4 window %.4f-%.4f, %s (half A x 2, rows at L_sim; X17 per unit BR)\n", gXLo, gXHi,
       CutStr(gCutX17nomee).c_str());
    PT("%-10s %24s %24s %24s\n", "component", "all M(ee) 0-100", Form("M(ee) %.1f-%.1f", gCut[sX17].mee_lo, gCut[sX17].mee_hi), "rows at L_data, window");
    for (int c = 0; c < NC; ++c) {
        const WSpec& w = gMeeH[0][c];
        double sa = 0, va = 0, sw = 0, vw = 0;
        for (int i = 0; i < w.nbin; ++i) {
            sa += w.s[i]; va += w.v[i];
            if (i >= gMwBin0 && i <= gMwBin1) { sw += w.s[i]; vw += w.v[i]; }
        }
        PT("%-10s %13.5g+-%-10.3g %13.5g+-%-10.3g %13.5g+-%-10.3g\n", kCK[c], sa, std::sqrt(va), sw, std::sqrt(vw), sw * kK, std::sqrt(vw) * kK);
    }
    PT("(bin sums: the bin errors are added in quadrature, the correlation of one event's rows across M(ee) bins is not included;\n"
       " the x17win column of X17SEL above has the exact group errors)\n");
    fclose(gT);
    gT = nullptr;
}

void ChecksPart1() {
    PT("## 1. frozen selections used (read from the frozen files, unchanged)\n");
    for (int s = 0; s < NS; ++s) PT("%-8s %s\n", kSN[s], CutStr(gCut[s]).c_str());
    PT("X17 M4 window %.4f-%.4f GeV/c2 = 2D M4 bins %d-%d (1-based); M(ee) window %.1f-%.1f MeV/c2 = bins %d-%d (1-based)\n\n", gXLo,
       gXHi, gXwBin0 + 1, gXwBin1 + 1, gCut[sX17].mee_lo, gCut[sX17].mee_hi, gMwBin0 + 1, gMwBin1 + 1);
    PT("## 2. cross-check: the frozen numbers recomputed on the FULL sample (both halves) with this macro's cut code\n");
    struct T { int s, reg; std::vector<int> sig; const char* what; };
    const T tg[3] = {{sPEAK, 3, {cPeak45}, "PEAK: classes 4+5, peak window"},
                     {sDIRECT, 3, {cDir5}, "DIRECT: class 5, peak window"},
                     {sBUMP, 1, {cBump23}, "BUMP: classes 2+3, bump window"}};
    PT("%-34s %10s %10s %10s %10s %10s %10s\n", "target", "FOM", "frozen", "purity", "frozen", "eff", "frozen");
    for (const T& t : tg) {
        double S = 0, Spre = 0;
        for (int c : t.sig) { S += RV(2, t.s, c, t.reg); Spre += RV(2, sPRE, c, t.reg); }
        const double tot = RV(2, t.s, cSum, t.reg);
        const double fom = S * kK / std::sqrt(tot * kK);
        PT("%-34s %10.3f %10.3f %10.4f %10.4f %10.4f %10.4f\n", t.what, fom, gCut[t.s].fom, S / tot, gCut[t.s].purity, S / Spre, gCut[t.s].eff);
    }
    const double xs = RV(2, sX17, cX17, NR) / kLsim, xb = RV(2, sX17, cSum, NR) / kLsim, xbe = RE(2, sX17, cSum, NR) / kLsim;
    PT("X17SEL in the X17 M4 window: X17 per unit BR %.5f nb (frozen %.5f); background T1 (no X17) + T2 %.5f +- %.5f nb (frozen %.5f +- %.5f)\n\n",
       xs, gX17sigFrozen, xb, xbe, gX17bkgFrozen, gX17bkgFrozenErr);
    // consistency of the window definitions with the 2D binning (counts of rows, all halves, X17 cuts without M(ee))
    long long nIn = 0, nMis = 0;
    for (int set = 1; set <= 3; ++set)
        for (const Row& r : gD->tree(set, kUS)) {
            if (!r.ownA || !PassCuts(r, gCutX17nomee)) continue;
            const int ix = XM4Bin(r.m4);
            const bool a = InX17Win(r.m4), b = ix >= gXwBin0 && ix <= gXwBin1;
            if (a) ++nIn;
            if (a != b) ++nMis;
            const int iy = XMeeBin(r.mee * 1000.f);
            const float m = r.mee * 1000.f;
            const bool c1 = (double)m >= gCut[sX17].mee_lo && (double)m < gCut[sX17].mee_hi, c2 = iy >= gMwBin0 && iy <= gMwBin1;
            if (c1 != c2) ++nMis;
        }
    PT("window bins of the X17 histograms against the direct window cuts: %lld rows in the M4 window, %lld mismatches (M4 or M(ee))\n",
       nIn, nMis);
    // X17 window sum of the 2D template against the exact window sum (half A)
    {
        double s2 = 0, s1 = 0;
        for (int ix = gXwBin0; ix <= gXwBin1; ++ix)
            for (int iy = gMwBin0; iy <= gMwBin1; ++iy) s2 += g2DH[0][cSum].s[ix * kXmeeN + iy];
        s1 = RV(0, sX17, cSum, NR);
        PT("half A, X17SEL in the X17 window: sum of the 2D window cells %.6f, exact window sum %.6f (expected rows at L_sim)\n\n", s2, s1);
    }
}

// ============================================================================
// part 2: eta pseudo-data (Poisson bootstrap of half B per event group)
// ============================================================================
struct Cand {
    float m4, mee, mm4, dz, lp1, lp2;   // lp1, lp2: true momenta of the lepton slots (lepton-pair identity, variance study)
    unsigned char pivtx, g, bits, xwin;
    signed char reg;
    double w, ratio;
};
struct CGrp { int first, n, set; double wg, L; };
std::vector<Cand> gC;
std::vector<CGrp> gG;
inline double U01(std::mt19937_64& g) { return (double)(g() >> 11) * (1.0 / 9007199254740992.0); }
inline int Pois(std::mt19937_64& g, double mean, double L) {
    if (mean > 30.) return std::poisson_distribution<int>(mean)(g);
    int n = 0;
    double p = U01(g);
    while (p > L) { ++n; p *= U01(g); }
    return n;
}
inline double SumMin(std::vector<double>& w) {   // sum over ordered pairs (r, s) of min(w_r, w_s)
    std::sort(w.begin(), w.end());
    const size_t n = w.size();
    double s = 0;
    for (size_t i = 0; i < n; ++i) s += w[i] * (1. + 2. * (double)(n - 1 - i));
    return s;
}

void BuildCandidates() {
    gC.clear(); gG.clear();
    long long nRow[4] = {0, 0, 0, 0};
    for (int set = 1; set <= 3; ++set) {
        const std::vector<Row>& v = gD->tree(set, kUS);
        int cur = -1;
        CGrp G{0, 0, set, 0., 1.};
        auto close = [&]() {
            if (G.n > 0) {
                for (int i = G.first; i < G.first + G.n; ++i) gC[i].ratio = gC[i].w / G.wg;
                G.L = std::exp(-G.wg);
                gG.push_back(G);
            }
        };
        for (const Row& r : v) {
            if (r.grp != cur) {
                close();
                cur = r.grp;
                G = CGrp{(int)gC.size(), 0, set, 0., 1.};
            }
            if ((r.grp & 1) != 1) continue;                // half B
            if (r.ownA != 1 || !PassPre(r)) continue;
            const int g = GroupOf(r);
            if (g < 0) continue;                           // T1 classes 2-5 and T2 only (X17, T3, T4 excluded)
            if (!(r.m4 >= kM4Lo && r.m4 < kM4Hi)) continue;
            const double w = Weight(r, gO);
            if (!(w > 0)) continue;
            Cand c;
            c.m4 = r.m4; c.mee = r.mee; c.mm4 = r.mm4; c.dz = DzPi(r);
            c.lp1 = r.ep_sim_p; c.lp2 = r.em_sim_p;
            c.pivtx = (r.pip_isUsedVertex == 1 && r.pim_isUsedVertex == 1) ? 1 : 0;
            c.g = (unsigned char)g;
            c.bits = 0;
            for (int s = 0; s < NSW; ++s) if (PassCuts(r, gCut[s])) c.bits |= (unsigned char)(1 << s);
            c.reg = (signed char)RegionOf(r.m4);
            c.xwin = InX17Win(r.m4) ? 1 : 0;
            c.w = w; c.ratio = 1.;
            gC.push_back(c);
            ++G.n;
            G.wg = std::max(G.wg, w);
            ++nRow[set];
        }
        close();
    }
    printf("half B candidates: %zu rows in %zu event groups (std %lld, lep %lld, eta %lld rows)\n", gC.size(), gG.size(), nRow[1], nRow[2],
           nRow[3]);
}

// expectations and variance predictions (half B): mu per selection x region x group (sealed), totals (public)
double gMuG[NSW][NRX][NG];
void Expectations() {
    for (int s = 0; s < NSW; ++s)
        for (int k = 0; k < NRX; ++k) {
            gMuTot[s][k] = gVcoh[s][k] = gVind[s][k] = 0;
            for (int g = 0; g < NG; ++g) gMuG[s][k][g] = 0;
        }
    std::vector<double> ws;
    std::vector<std::pair<std::pair<float, float>, double>> wp;
    for (const CGrp& G : gG) {
        for (int s = 0; s < NSW; ++s)
            for (int k = 0; k < NRX; ++k) {
                ws.clear(); wp.clear();
                for (int i = G.first; i < G.first + G.n; ++i) {
                    const Cand& c = gC[i];
                    if (!(c.bits & (1 << s))) continue;
                    if (k < NR ? c.reg != k : !c.xwin) continue;
                    ws.push_back(c.w);
                    wp.push_back(std::make_pair(std::make_pair(c.lp1, c.lp2), c.w));
                    gMuG[s][k][c.g] += c.w;
                    gMuTot[s][k] += c.w;
                }
                if (ws.empty()) continue;
                gVcoh[s][k] += SumMin(ws);
                std::sort(wp.begin(), wp.end());
                size_t a = 0;
                while (a < wp.size()) {
                    size_t b = a;
                    std::vector<double> sub;
                    while (b < wp.size() && wp[b].first == wp[a].first) { sub.push_back(wp[b].second); ++b; }
                    gVind[s][k] += SumMin(sub);
                    a = b;
                }
            }
    }
}

void PartEtaToys(ULong64_t seed, int nToy) {
    TStopwatch sw;
    BuildCandidates();
    Expectations();
    gNEta = nToy;
    std::mt19937_64 rng(seed);
    gSystem->mkdir(BDir("pseudo").c_str(), kTRUE);
    gSystem->mkdir(BDir("sealed").c_str(), kTRUE);
    const std::string fn = BDir("pseudo") + "/eta_toys.root";
    TFile* f = TFile::Open(fn.c_str(), "RECREATE");
    TTree* t = new TTree("eta_toys", "pi+pi-e+e- pseudo-data rows (pre applied), M4 0.30-1.20 GeV/c2");
    Int_t toy, evt, pivtx;
    Float_t m4, mee, mm4, dz;
    t->Branch("toy", &toy, "toy/I");
    t->Branch("evt", &evt, "evt/I");
    t->Branch("m4", &m4, "m4/F");
    t->Branch("mee", &mee, "mee/F");
    t->Branch("mm4", &mm4, "mm4/F");
    t->Branch("dz_pi", &dz, "dz_pi/F");
    t->Branch("pivtx", &pivtx, "pivtx/I");
    // sealed key tree
    const std::string kn = BDir("sealed") + "/eta_key.root";
    TFile* fk = TFile::Open(kn.c_str(), "RECREATE");
    TTree* tk = new TTree("key", "SEALED: true rows per group of every eta toy, selection and M4 region");
    Int_t ktoy, ksel, kreg, kn9[NG], ktot;
    tk->Branch("toy", &ktoy, "toy/I");
    tk->Branch("sel", &ksel, "sel/I");
    tk->Branch("reg", &kreg, "reg/I");
    tk->Branch("n", kn9, Form("n[%d]/I", NG));
    tk->Branch("ntot", &ktot, "ntot/I");
    std::vector<int> K((size_t)nToy * NSW * NRX * NG, 0);
    auto KI = [&](int tt, int s, int k, int g) -> int& { return K[(((size_t)tt * NSW + s) * NRX + k) * NG + g]; };
    for (int s = 0; s < NSW; ++s) for (int k = 0; k < NRX; ++k) gEtaTot.n[s][k].assign(nToy, 0.);
    struct Ev { int g; double u; };
    std::vector<Ev> ev;
    gEtaRowsTotal = gEtaEvtsTotal = 0;
    for (int tt = 0; tt < nToy; ++tt) {
        ev.clear();
        for (size_t gi = 0; gi < gG.size(); ++gi) {
            const CGrp& G = gG[gi];
            const int n = Pois(rng, G.wg, G.L);
            for (int k = 0; k < n; ++k) ev.push_back(Ev{(int)gi, U01(rng)});
        }
        std::shuffle(ev.begin(), ev.end(), rng);
        f->cd();
        for (size_t e = 0; e < ev.size(); ++e) {
            const CGrp& G = gG[ev[e].g];
            const double u = ev[e].u;
            for (int i = G.first; i < G.first + G.n; ++i) {
                const Cand& c = gC[i];
                if (!(c.ratio > u)) continue;
                toy = tt; evt = (Int_t)e; m4 = c.m4; mee = c.mee; mm4 = c.mm4; dz = c.dz; pivtx = c.pivtx;
                t->Fill();
                ++gEtaRowsTotal;
                for (int s = 0; s < NSW; ++s) {
                    if (!(c.bits & (1 << s))) continue;
                    if (c.reg >= 0) { KI(tt, s, c.reg, c.g)++; gEtaTot.n[s][c.reg][tt] += 1; }
                    if (c.xwin) { KI(tt, s, NR, c.g)++; gEtaTot.n[s][NR][tt] += 1; }
                }
            }
        }
        gEtaEvtsTotal += (long long)ev.size();
        if (tt % 20 == 0) printf("eta toy %d done (%.0f s)\n", tt, sw.RealTime()), sw.Continue();
    }
    f->cd();
    t->Write("", TObject::kOverwrite);
    TNamed("info", "Eta pseudo-data of the blind test: rows passing pre (isBest >= 0, vz > -500 mm, PT3 bit, vz > -160 mm, cor_ok), "
                   "0.30 <= M4 < 1.20 GeV/c2, at an integrated luminosity of 110.525 nb^-1. Columns: toy (toy index), evt (event number in the "
                   "toy; rows of one event share it), m4 [GeV/c2], mee [GeV/c2], mm4 [GeV/c2], dz_pi = max(|z_pi+ - z_vtx|, "
                   "|z_pi- - z_vtx|) [mm], pivtx (1 = both pions used in the vertex fit). Unweighted: one entry = one row.")
        .Write();
    f->Close();
    delete f;
    // sealed key
    fk->cd();
    for (int tt = 0; tt < nToy; ++tt)
        for (int s = 0; s < NSW; ++s)
            for (int k = 0; k < NRX; ++k) {
                ktoy = tt; ksel = s; kreg = k; ktot = 0;
                for (int g = 0; g < NG; ++g) { kn9[g] = KI(tt, s, k, g); ktot += kn9[g]; }
                tk->Fill();
            }
    tk->Write();
    TTree* te = new TTree("expect", "SEALED: half-B expectation per toy (sum of w), selection, M4 region, group");
    Int_t esel, ereg;
    Double_t emu[NG], emut, evc, evi;
    te->Branch("sel", &esel, "sel/I");
    te->Branch("reg", &ereg, "reg/I");
    te->Branch("mu", emu, Form("mu[%d]/D", NG));
    te->Branch("mutot", &emut, "mutot/D");
    te->Branch("vcoh", &evc, "vcoh/D");
    te->Branch("vind", &evi, "vind/D");
    for (int s = 0; s < NSW; ++s)
        for (int k = 0; k < NRX; ++k) {
            esel = s; ereg = k; emut = gMuTot[s][k]; evc = gVcoh[s][k]; evi = gVind[s][k];
            for (int g = 0; g < NG; ++g) emu[g] = gMuG[s][k][g];
            te->Fill();
        }
    te->Write();
    TNamed("info", "SEALED key of the eta pseudo-data. sel: 0 PEAK, 1 DIRECT, 2 BUMP, 3 X17SEL; reg: 0 below 0.4576, 1 bump window "
                   "0.4576-0.5255, 2 gap, 3 peak window 0.5333-0.5587, 4 0.5587-0.70, 5 0.70-1.00, 6 1.00-1.20, 7 X17 M4 window; "
                   "n[9]: rows of groups Gsig2, Gsig3, Gsig4, Gsig5, G_eta, G_pi0D, G_conv, G_vec, G_oth.")
        .Write();
    fk->Close();
    delete fk;

    // sealed text
    PK("\n## ETA PSEUDO-DATA (half B, L_sim / 2 = 110.525 nb^-1): %d toys\n", nToy);
    PK("# true rows per toy, derived from the sampled rows. Targets: PEAK = Gsig4 + Gsig5 in the peak window (PEAK selection);\n");
    PK("# DIRECT = Gsig5 in the peak window (DIRECT selection); BUMP = Gsig2 + Gsig3 in the bump window (BUMP selection)\n");
    PK("# expectation per toy (half-B sum of w): PEAK %.3f, DIRECT %.3f, BUMP %.3f\n", gMuG[sPEAK][3][2] + gMuG[sPEAK][3][3],
       gMuG[sDIRECT][3][3], gMuG[sBUMP][1][0] + gMuG[sBUMP][1][1]);
    PK("%5s %9s %9s %9s | %8s %8s %8s %8s %8s | %8s %8s | %8s %8s %8s\n", "toy", "PEAK", "DIRECT", "BUMP", "PK:Gs4", "PK:Gs5", "PK:tot",
       "DI:Gs5", "DI:tot", "BU:Gs2", "BU:Gs3", "BU:tot", "X17:Gs5", "X17:tot");
    for (int tt = 0; tt < nToy; ++tt) {
        const int pk = KI(tt, sPEAK, 3, 2) + KI(tt, sPEAK, 3, 3), di = KI(tt, sDIRECT, 3, 3), bu = KI(tt, sBUMP, 1, 0) + KI(tt, sBUMP, 1, 1);
        PK("%5d %9d %9d %9d | %8d %8d %8.0f %8d %8.0f | %8d %8d | %8.0f %8d %8.0f\n", tt, pk, di, bu, KI(tt, sPEAK, 3, 2), KI(tt, sPEAK, 3, 3),
           gEtaTot.n[sPEAK][3][tt], KI(tt, sDIRECT, 3, 3), gEtaTot.n[sDIRECT][3][tt], KI(tt, sBUMP, 1, 0), KI(tt, sBUMP, 1, 1),
           gEtaTot.n[sBUMP][1][tt], KI(tt, sX17, NR, 3), gEtaTot.n[sX17][NR][tt]);
    }
    // closure of the bootstrap per group: toy mean against the expectation
    PK("\n# bootstrap closure (sealed): toy mean of the true rows against the half-B expectation, per selection x region x group\n");
    PK("%-7s %-8s %-7s %12s %12s %12s %8s\n", "sel", "region", "group", "expect", "toy mean", "toy rms", "pull");
    int nBad = 0, nTest = 0;
    for (int s = 0; s < NSW; ++s)
        for (int k = 0; k < NRX; ++k)
            for (int g = 0; g < NG; ++g) {
                double m = 0, m2 = 0;
                for (int tt = 0; tt < nToy; ++tt) { const double x = KI(tt, s, k, g); m += x; m2 += x * x; }
                m /= nToy;
                const double rms = std::sqrt(std::max(0., m2 / nToy - m * m));
                const double mu = gMuG[s][k][g];
                const double pull = rms > 0 ? (m - mu) / (rms / std::sqrt((double)nToy)) : 0.;
                if (mu > 0.5 || m > 0.5) {
                    ++nTest;
                    if (std::fabs(pull) > 3.5) ++nBad;
                    PK("%-7s %-8s %-7s %12.3f %12.3f %12.3f %8.2f\n", kSN[s], RegShort(k).c_str(), kGN[g], mu, m, rms, pull);
                }
            }
    PK("bootstrap closure: %d of %d entries with |pull| > 3.5\n", nBad, nTest);
    printf("eta toys: %d toys, %lld rows, %lld events written to %s (%.0f s); key written to the sealed directory\n", nToy, gEtaRowsTotal,
           gEtaEvtsTotal, fn.c_str(), sw.RealTime());
}

// ============================================================================
// part 4: X17 pseudo-data (parametric, L_data UNCONFIRMED)
// ============================================================================
void PartX17Toys(ULong64_t seed, int nPer) {
    TStopwatch sw;
    std::mt19937_64 rng(seed);
    const int ncell = kXm4N * kXmeeN;
    // expectations from half B (x 2 already) scaled to L_data
    std::vector<double> mub(ncell), s1(ncell);
    for (int k = 0; k < ncell; ++k) { mub[k] = g2DH[1][cSum].s[k] * kK; s1[k] = g2DH[1][cX17].s[k] * kK; }
    double mubW = 0, s1W = 0, mubA = 0, s1A = 0;
    for (int ix = 0; ix < kXm4N; ++ix)
        for (int iy = 0; iy < kXmeeN; ++iy) {
            const int k = ix * kXmeeN + iy;
            mubA += mub[k]; s1A += s1[k];
            if (ix >= gXwBin0 && ix <= gXwBin1 && iy >= gMwBin0 && iy <= gMwBin1) { mubW += mub[k]; s1W += s1[k]; }
        }
    // toy list: (BR index, replica), shuffled; random unique six-digit ids
    struct XT { int ibr, rep, id; };
    std::vector<XT> L;
    for (int b = 0; b < NBR; ++b) for (int i = 0; i < nPer; ++i) L.push_back(XT{b, i, 0});
    std::shuffle(L.begin(), L.end(), rng);
    std::set<int> used;
    std::uniform_int_distribution<int> did(100000, 999999);
    for (XT& x : L) { int id; do { id = did(rng); } while (used.count(id)); used.insert(id); x.id = id; }
    gNX17 = (int)L.size();
    gSystem->mkdir(BDir("pseudo").c_str(), kTRUE);
    const std::string fn = BDir("pseudo") + "/x17_toys.root";
    TFile* f = TFile::Open(fn.c_str(), "RECREATE");
    TTree* tl = new TTree("toys", "list of the X17 pseudo-data toys in file order");
    Int_t tid, tidx;
    tl->Branch("id", &tid, "id/I");
    tl->Branch("index", &tidx, "index/I");
    // sealed key
    const std::string kn = BDir("sealed") + "/x17_key.root";
    TFile* fk = TFile::Open(kn.c_str(), "RECREATE");
    TTree* tk = new TTree("key", "SEALED: injected BR and counts of every X17 toy");
    Int_t kid, kidx, kibr, ns, nb, nobs, nsA, nobsA;
    Double_t kbr, kmus, kmub, kmusA, kmubA;
    tk->Branch("id", &kid, "id/I");
    tk->Branch("index", &kidx, "index/I");
    tk->Branch("ibr", &kibr, "ibr/I");
    tk->Branch("br", &kbr, "br/D");
    tk->Branch("mu_s_win", &kmus, "mu_s_win/D");
    tk->Branch("mu_b_win", &kmub, "mu_b_win/D");
    tk->Branch("n_s_win", &ns, "n_s_win/I");
    tk->Branch("n_b_win", &nb, "n_b_win/I");
    tk->Branch("n_obs_win", &nobs, "n_obs_win/I");
    tk->Branch("mu_s_all", &kmusA, "mu_s_all/D");
    tk->Branch("mu_b_all", &kmubA, "mu_b_all/D");
    tk->Branch("n_s_all", &nsA, "n_s_all/I");
    tk->Branch("n_obs_all", &nobsA, "n_obs_all/I");
    gX17obsWin.assign(L.size(), 0.);
    PK("\n## X17 PSEUDO-DATA at L_data = 6080.5 nb^-1 (UNCONFIRMED): %d toys (%d per BR)\n", gNX17, nPer);
    PK("# expectation from half B (x 2 x L_data / L_sim), %s, M4 0.4830-0.6080 x M(ee) 0-100 MeV/c2:\n", CutStr(gCutX17nomee).c_str());
    PK("#   background T1 (no X17) + T2: frozen window %.3f rows, whole 2D range %.3f rows\n", mubW, mubA);
    PK("#   X17 per unit BR: frozen window %.6g rows, whole 2D range %.6g rows\n", s1W, s1A);
    PK("%6s %8s %5s %10s %10s %10s %8s %8s %9s %8s %9s\n", "index", "id", "ibr", "BR", "mu_s_win", "mu_b_win", "n_s_win", "n_b_win",
       "n_obs_win", "n_s_all", "n_obs_all");
    for (size_t i = 0; i < L.size(); ++i) {
        const double br = kBR[L[i].ibr];
        TH2D* h = new TH2D(Form("toy_%06d_m4mee", L[i].id), Form("pseudo-data toy %06d: M4 x M(ee)", L[i].id), kXm4N, kXm4Lo,
                           kXm4Lo + kXm4N * kXm4W, kXmeeN, 0., kXmeeN * kXmeeW);
        h->SetDirectory(nullptr);
        h->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]");
        h->GetYaxis()->SetTitle("M(e^{+}e^{-}) [MeV/c^{2}]");
        h->GetZaxis()->SetTitle("rows");
        TH1D* hm = new TH1D(Form("toy_%06d_mee", L[i].id), Form("pseudo-data toy %06d: M(ee) in the X17 M4 window", L[i].id), kXmeeN, 0.,
                            kXmeeN * kXmeeW);
        hm->SetDirectory(nullptr);
        hm->GetXaxis()->SetTitle("M(e^{+}e^{-}) [MeV/c^{2}]");
        hm->GetYaxis()->SetTitle("rows / 0.5 MeV/c^{2}");
        ns = nb = nobs = nsA = nobsA = 0;
        for (int ix = 0; ix < kXm4N; ++ix)
            for (int iy = 0; iy < kXmeeN; ++iy) {
                const int k = ix * kXmeeN + iy;
                const int b = mub[k] > 0 ? Pois(rng, mub[k], std::exp(-mub[k])) : 0;
                const double ms = br * s1[k];
                const int s = ms > 0 ? Pois(rng, ms, std::exp(-ms)) : 0;
                const int n = b + s;
                h->SetBinContent(ix + 1, iy + 1, n);
                if (ix >= gXwBin0 && ix <= gXwBin1) hm->SetBinContent(iy + 1, hm->GetBinContent(iy + 1) + n);
                nsA += s; nobsA += n;
                if (ix >= gXwBin0 && ix <= gXwBin1 && iy >= gMwBin0 && iy <= gMwBin1) { ns += s; nb += b; nobs += n; }
            }
        h->SetEntries(nobsA);
        hm->SetEntries(hm->Integral());
        f->cd();
        h->Write();
        hm->Write();
        delete h;
        delete hm;
        tid = L[i].id; tidx = (Int_t)i;
        tl->Fill();
        kid = L[i].id; kidx = (Int_t)i; kibr = L[i].ibr; kbr = br; kmus = br * s1W; kmub = mubW; kmusA = br * s1A; kmubA = mubA;
        tk->Fill();
        gX17obsWin[i] = nobs;
        PK("%6zu %08d %5d %10.2g %10.3f %10.3f %8d %8d %9d %8d %9d\n", i, L[i].id, L[i].ibr, br, kmus, kmub, ns, nb, nobs, nsA, nobsA);
    }
    f->cd();
    tl->Write();
    TNamed("info", Form("X17 pseudo-data of the blind test at L_data = 6080.5 nb^-1 (UNCONFIRMED): %d toys in random order with random "
                        "six-digit ids (tree 'toys'). Per toy: toy_<id>_m4mee = observed rows in M4 0.4830-0.6080 GeV/c2 (50 bins of "
                        "2.5 MeV/c2) x M(ee) 0-100 MeV/c2 (200 bins of 0.5 MeV/c2), after %s; toy_<id>_mee = its M(ee) projection "
                        "over the X17 M4 window %.4f-%.4f GeV/c2 (M4 bins %d-%d). Observed rows only.",
                        gNX17, CutStr(gCutX17nomee).c_str(), gXLo, gXHi, gXwBin0 + 1, gXwBin1 + 1))
        .Write();
    f->Close();
    delete f;
    fk->cd();
    tk->Write();
    TNamed("info", "SEALED key of the X17 pseudo-data: per toy id the injected BR (eta -> pi+pi- X17, X17 -> e+e-), the expected "
                   "signal and background rows in the frozen window (M4 and M(ee)) and in the whole 2D range, and the sampled "
                   "signal / background / observed counts.")
        .Write();
    fk->Close();
    delete fk;
    printf("X17 toys: %d toys written to %s (%.0f s); key written to the sealed directory\n", gNX17, fn.c_str(), sw.RealTime());
}

// ============================================================================
// sealed comparison of the halves (templates against the pseudo-data source)
// ============================================================================
void SealedHalves() {
    PK("\n## HALF A (templates) against HALF B (pseudo-data source), both x 2 (expected rows at L_sim); difference / sqrt(VA + VB)\n");
    for (int s = 0; s < NSW; ++s) {
        PK("# selection %s = %s\n", kSN[s], CutStr(gCut[s]).c_str());
        PK("%-10s %-8s %14s %14s %8s\n", "component", "region", "half A", "half B", "pull");
        for (int c = 0; c < NC; ++c)
            for (int k = 0; k < NRX; ++k) {
                const double a = RV(0, s, c, k), b = RV(1, s, c, k), ea = RE(0, s, c, k), eb = RE(1, s, c, k);
                if (a == 0 && b == 0) continue;
                const double pl = (ea * ea + eb * eb) > 0 ? (b - a) / std::sqrt(ea * ea + eb * eb) : 0.;
                PK("%-10s %-8s %8.4g+-%-6.2g %8.4g+-%-6.2g %8.2f\n", kCK[c], RegShort(k).c_str(), a, ea, b, eb, pl);
            }
    }
    PK("\n# X17 window counts at L_data (x %.4f): half A / half B\n", kK);
    for (int c : {cSum, cDir5, cConv4, cOther, cX17}) {
        double a = 0, b = 0, va = 0, vb = 0;
        for (int ix = gXwBin0; ix <= gXwBin1; ++ix)
            for (int iy = gMwBin0; iy <= gMwBin1; ++iy) {
                const int k = ix * kXmeeN + iy;
                a += g2DH[0][c].s[k]; b += g2DH[1][c].s[k]; va += g2DH[0][c].v[k]; vb += g2DH[1][c].v[k];
            }
        PK("%-10s %12.5g +- %-10.3g %12.5g +- %-10.3g (cell sums, per unit BR for x17perBR)\n", kCK[c], a * kK, std::sqrt(va) * kK, b * kK,
           std::sqrt(vb) * kK);
    }
}

// ============================================================================
// public checks after the toys (totals only, no composition)
// ============================================================================
void ChecksToys() {
    PT("## 3. split into halves (public part; the half-B numbers per set and per component are sealed)\n");
    {
        long long rows[4][2] = {{0, 0}, {0, 0}, {0, 0}, {0, 0}};
        double sw[4][2] = {{0, 0}, {0, 0}, {0, 0}, {0, 0}};
        std::set<int> grps[4][2];
        for (int set = 1; set <= 3; ++set)
            for (const Row& r : gD->tree(set, kUS)) {
                if (r.ownA != 1 || !PassPre(r) || GroupOf(r) < 0) continue;
                const double w = Weight(r, gO);
                if (!(w > 0)) continue;
                const int h = r.grp & 1;
                ++rows[set][h]; sw[set][h] += w; grps[set][h].insert(r.grp);
            }
        PT("rule: half A = event groups with even grp, half B = odd grp, per e+e- tree (grp = event-group index in file order)\n");
        PT("rows: pre, owned (mode A), T1 classes 2-5 + T2, w x phys > 0, 0.30 <= M4 < 1.20\n");
        PT("%-5s %12s %12s %16s\n", "set", "groups A", "rows A", "sum w A (x 1)");
        for (int set = 1; set <= 3; ++set)
            PT("%-5s %12zu %12lld %16.4f\n", SetName(set), grps[set][0].size(), rows[set][0], sw[set][0]);
        PT("half B in total: %zu groups, %lld rows (= the candidates of the eta pseudo-data)\n\n",
           grps[1][1].size() + grps[2][1].size() + grps[3][1].size(), rows[1][1] + rows[2][1] + rows[3][1]);
    }
    if (gNEta > 0) {
        PT("## 4. eta pseudo-data: %d toys at L_sim / 2 = 110.525 nb^-1, %lld rows and %lld events in total (%.0f rows per toy)\n", gNEta,
           gEtaRowsTotal, gEtaEvtsTotal, (double)gEtaRowsTotal / gNEta);
        PT("# totals of ALL rows (signal + background together; observable) per selection and M4 region across the toys:\n");
        PT("# mean and variance of the toys against the bootstrap prediction (half-B sum of w; variance of the coherent-thinning model\n");
        PT("# V_coh = sum over groups of sum_rs min(w_r, w_s); V_ind = the same with independent lepton pairs; plain Poisson = mean)\n");
        PT("%-7s %-8s %10s %10s %8s %10s %10s %10s %10s %7s\n", "sel", "region", "toy mean", "predicted", "pull", "toy var", "V_coh",
           "V_ind", "Poisson", "var/Vc");
        for (int s = 0; s < NSW; ++s)
            for (int k = 0; k < NRX; ++k) {
                double m = 0, m2 = 0;
                for (int tt = 0; tt < gNEta; ++tt) { const double x = gEtaTot.n[s][k][tt]; m += x; m2 += x * x; }
                m /= gNEta;
                const double var = std::max(0., (m2 - gNEta * m * m) / (gNEta - 1));
                const double pull = gVcoh[s][k] > 0 ? (m - gMuTot[s][k]) / std::sqrt(gVcoh[s][k] / gNEta) : 0.;
                PT("%-7s %-8s %10.2f %10.2f %8.2f %10.1f %10.1f %10.1f %10.1f %7.3f\n", kSN[s], RegShort(k).c_str(), m, gMuTot[s][k], pull,
                   var, gVcoh[s][k], gVind[s][k], gMuTot[s][k], gVcoh[s][k] > 0 ? var / gVcoh[s][k] : 0.);
            }
        PT("(expected spread of var/V_coh for %d toys: about +-%.3f)\n\n", gNEta, std::sqrt(2. / (gNEta - 1)));
    }
    if (gNX17 > 0) {
        double m = 0, mn = 1e30, mx = -1e30;
        for (double x : gX17obsWin) { m += x; mn = std::min(mn, x); mx = std::max(mx, x); }
        m /= gX17obsWin.size();
        PT("## 5. X17 pseudo-data: %d toys (five BR values mixed, shuffled, random ids); observed rows in the frozen window: mean %.1f, "
           "min %.0f, max %.0f\n\n", gNX17, m, mn, mx);
    }
}

// ============================================================================
// READMEs of the pseudo-data and sealed directories
// ============================================================================
void WriteReadmes(int nEta, int nX17per) {
    FILE* r = fopen((BDir("pseudo") + "/README.md").c_str(), "w");
    fprintf(r, "# Pseudo-data of the blind test\n\n");
    fprintf(r, "Written by `plots/eta_sim/s5_blind_prep.C`. SIMULATION ONLY (SMASH std + lep + eta sets, rows owned in partition "
               "mode A, PDG weights and the physics weights of the eta decays, main variant V3). There is no combinatorial "
               "background in the simulation and none is added: nothing has to be subtracted. The pseudo-data contain the rows of "
               "template 1 (signal classes 2-5) and template 2 (physical background); fake-lepton and different-source rows "
               "(templates 3, 4) are not included. They come from half B of the simulation (event groups with an odd group "
               "number); the templates (`../templates/templates.root`) come from the independent half A.\n\n");
    fprintf(r, "Do not open `../sealed/` before the extraction procedure and its results are frozen.\n\n");
    fprintf(r, "## `eta_toys.root`: eta pseudo-data, %d toys at L = L_sim / 2 = 110.525 nb^-1\n\n", nEta);
    fprintf(r, "Tree `eta_toys`, one entry per row (combination), unweighted. Every row already passes the base selection pre "
               "(isBest >= 0, vertex z > -500 mm, PT3 bit, vertex z > -160 mm, cor_ok == 1); 0.30 <= M4 < 1.20 GeV/c2. Each toy is "
               "a Poisson bootstrap of half B per event group at the luminosity of half B (statistics of a real measurement with "
               "L = 110.5 nb^-1), so the toys differ by honest Poisson fluctuations but share the same parent sample.\n\n");
    fprintf(r, "| branch | type | meaning |\n|---|---|---|\n");
    fprintf(r, "| `toy` | Int_t | toy index 0..%d |\n", nEta - 1);
    fprintf(r, "| `evt` | Int_t | event number inside the toy (events in random order); rows of one event share it. isBest >= 0 keeps "
               "several rows per event: use `evt` for event-correlated errors |\n");
    fprintf(r, "| `m4` | Float_t | M(pi+pi-e+e-) [GeV/c2], corrected momenta |\n");
    fprintf(r, "| `mee` | Float_t | M(e+e-) [GeV/c2], corrected momenta (the X17 cuts use `mee*1000.f` in MeV/c2) |\n");
    fprintf(r, "| `mm4` | Float_t | missing mass of pi+pi-e+e- [GeV/c2], corrected momenta |\n");
    fprintf(r, "| `dz_pi` | Float_t | max(abs(z_pi+ - z_vertex), abs(z_pi- - z_vertex)) [mm] |\n");
    fprintf(r, "| `pivtx` | Int_t | 1 = both pions used in the vertex fit |\n\n");
    fprintf(r, "Also `info` (TNamed): a short description. Frozen selections (apply on top of pre, which is already applied):\n\n");
    fprintf(r, "- PEAK: `dz_pi < 40.f && mm4 < 2.85f`; target window 0.5333 <= m4 < 0.5587\n");
    fprintf(r, "- DIRECT: `dz_pi < 40.f`; target window 0.5333 <= m4 < 0.5587\n");
    fprintf(r, "- BUMP: no further cut; target window 0.4576 <= m4 < 0.5255\n");
    fprintf(r, "- X17SEL: `pivtx == 1 && mee*1000.f >= %.1f && mee*1000.f < %.1f`; M4 window %.4f <= m4 < %.4f (there is no X17 in the "
               "eta pseudo-data)\n\n", gCut[sX17].mee_lo, gCut[sX17].mee_hi, gXLo, gXHi);
    fprintf(r, "(the exact definitions are in `../templates/README.md`; lower window edges included, upper excluded)\n\n");
    fprintf(r, "## `x17_toys.root`: X17 pseudo-data, %d toys at L_data = 6080.5 nb^-1 (UNCONFIRMED)\n\n", NBR * nX17per);
    fprintf(r, "Parametric toys: Poisson counts per cell around an expectation taken from half B, with an X17 contribution "
               "(eta -> pi+pi- X17, X17 -> e+e-) of unknown branching ratio injected in some of them. The toys are in random order "
               "and carry random six-digit ids; the injected values are only in the sealed key.\n\n");
    fprintf(r, "- tree `toys`: `id` (Int_t, the toy id) and `index` (Int_t, position in the file)\n");
    fprintf(r, "- `toy_<id>_m4mee` (TH2D): observed rows, x = M(pi+pi-e+e-) 0.4830-0.6080 GeV/c2 (50 bins of 2.5 MeV/c2), y = M(e+e-) "
               "0-100 MeV/c2 (200 bins of 0.5 MeV/c2), after %s\n", CutStr(gCutX17nomee).c_str());
    fprintf(r, "- `toy_<id>_mee` (TH1D): its M(e+e-) projection over the X17 M4 window %.4f-%.4f GeV/c2 (x bins %d-%d, 1-based)\n",
            gXLo, gXHi, gXwBin0 + 1, gXwBin1 + 1);
    fprintf(r, "- frozen M(e+e-) window %.1f-%.1f MeV/c2 = y bins %d-%d (1-based)\n", gCut[sX17].mee_lo, gCut[sX17].mee_hi, gMwBin0 + 1,
            gMwBin1 + 1);
    fprintf(r, "- `info` (TNamed): a short description\n\n");
    fprintf(r, "Templates for the X17 analysis: `../templates/templates.root`, histograms `mee_X17_*` and `m4mee_X17_*` (same binning; "
               "they are at L_sim = 221.05 nb^-1: multiply by 6080.5 / 221.05 = %.4f for these toys).\n", kK);
    fclose(r);

    FILE* k = fopen((BDir("sealed") + "/README.md").c_str(), "w");
    fprintf(k, "# SEALED key of the blind test (s5_blind_prep.C)\n\n");
    fprintf(k, "Not for the analyst before the extraction is frozen.\n\n");
    fprintf(k, "| file | content |\n|---|---|\n");
    fprintf(k, "| `seed.txt` | master seed of the run (reproduces every toy with `s5_blind_prep.C+(1, %d, %d)`, which reads this file) |\n", nEta, nX17per);
    fprintf(k, "| `key_summary.txt` | half A against half B per component, selection and M4 region; true rows of every eta toy for the "
               "three targets (PEAK = Gsig4 + Gsig5 in the peak window, DIRECT = Gsig5 in the peak window, BUMP = Gsig2 + Gsig3 in "
               "the bump window) with the expectation; bootstrap closure per group; X17 expectations and every X17 toy (id, BR, "
               "expected and sampled signal / background in the frozen window) |\n");
    fprintf(k, "| `eta_key.root` | tree `key` (toy, sel 0 PEAK / 1 DIRECT / 2 BUMP / 3 X17SEL, reg 0-7, n[9] rows of Gsig2, Gsig3, "
               "Gsig4, Gsig5, G_eta, G_pi0D, G_conv, G_vec, G_oth, ntot); tree `expect` (sel, reg, mu[9] = half-B sum of w per toy, "
               "mutot, vcoh, vind); regions 0 below 0.4576, 1 bump window, 2 gap, 3 peak window, 4 0.5587-0.70, 5 0.70-1.00, "
               "6 1.00-1.20, 7 X17 M4 window |\n");
    fprintf(k, "| `x17_key.root` | tree `key` (id, index, ibr, br, mu_s_win, mu_b_win, n_s_win, n_b_win, n_obs_win, mu_s_all, mu_b_all, "
               "n_s_all, n_obs_all) |\n");
    fclose(k);
}

// ============================================================================
// figures
// ============================================================================
TH1D* Reb(TH1D* h, int n, const char* nm) {
    TH1D* r = (TH1D*)h->Rebin(n, nm);
    r->SetDirectory(nullptr);
    return r;
}
TH1D* Dummy(int fillCol, int lineCol, int lineStyle, int lineWidth) {
    TH1D* h = new TH1D(Nm("dum").c_str(), "", 1, 0, 1);
    h->SetDirectory(nullptr);
    if (fillCol >= 0) { h->SetFillColor(fillCol); h->SetFillStyle(1001); } else h->SetFillStyle(0);
    h->SetLineColor(lineCol); h->SetLineStyle(lineStyle); h->SetLineWidth(lineWidth);
    return h;
}
void FigTemplatesEta() {
    TCanvas* c = new TCanvas("cTplEta", "", 2400, 2850);
    c->SetFillColor(0);
    const int comps[5] = {cBkg, cBump2, cBump3, cConv4, cDir5};
    const int fill[5] = {kGray, kOrange + 1, kOrange + 7, kRed + 1, kPink + 7};
    const char* lab[5] = {"physical background (T2)", "bump, #pi^{0} Dalitz (class 2)", "bump, #pi^{0} #gamma conv. (class 3)",
                          "peak, #eta#rightarrow#pi^{+}#pi^{-}#gamma, #gamma conv. (class 4)",
                          "peak, direct #eta#rightarrow#pi^{+}#pi^{-}e^{+}e^{-} (class 5)"};
    TPad* hd = new TPad("hd", "", 0, 0.88, 1, 1);
    hd->Draw();
    hd->cd();
    Txt(0.01, 0.84, 0.125, "Blind-test templates from half A (weights #times 2 = full luminosity of one set, L_{sim} = 221.05 nb^{-1})");
    Txt(0.01, 0.68, 0.125, "SMASH std + lep + eta, owned rows (mode A), w (PDG) #times physics weights V3; simulation only, no CB");
    {
        TLegend* lg = new TLegend(0.01, 0.02, 0.99, 0.62);
        lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.115); lg->SetNColumns(3);
        for (int i = 4; i >= 0; --i) lg->AddEntry(Dummy(fill[i], i == 0 ? kGray + 2 : fill[i], 1, 1), lab[i], "f");
        lg->AddEntry(Dummy(-1, kAzure + 1, 2, 4), "fake lepton (T3, variant, not summed)", "l");
        lg->AddEntry(Dummy(-1, kBlue + 3, 7, 4), "different sources (T4, variant, not summed)", "l");
        lg->Draw();
    }
    c->cd();
    for (int s = 0; s < 3; ++s) {
        for (int col = 0; col < 2; ++col) {
            const double x0 = col == 0 ? 0 : 0.5, x1 = col == 0 ? 0.5 : 1.0;
            const double y1 = 0.88 - s * 0.2933, y0 = y1 - 0.2933;
            c->cd();
            TPad* p = new TPad(Nm("p").c_str(), "", x0, y0, x1, y1);
            p->SetLeftMargin(0.16); p->SetRightMargin(0.03); p->SetBottomMargin(0.16); p->SetTopMargin(0.10);
            if (col == 0) p->SetLogy();
            p->Draw();
            p->cd();
            const int rb = col == 0 ? 5 : 2;
            std::vector<TH1D*> hs;   // cumulative: hs[i] = components 0..i
            for (int i = 0; i < 5; ++i) {
                TH1D* h0 = gM4A[s][comps[i]].Hist(Nm("h").c_str(), 0);
                TH1D* h = Reb(h0, rb, Nm("hr").c_str());
                if (i) h->Add(hs[i - 1]);
                h->SetFillColor(fill[i]); h->SetFillStyle(1001); h->SetLineColor(i == 0 ? kGray + 2 : fill[i]); h->SetLineWidth(1);
                hs.push_back(h);
            }
            TH1D* h3 = Reb(gM4A[s][cFake].Hist(Nm("h").c_str(), 0), rb, Nm("hr").c_str());
            TH1D* h4 = Reb(gM4A[s][cDiff].Hist(Nm("h").c_str(), 0), rb, Nm("hr").c_str());
            h3->SetLineColor(kAzure + 1); h3->SetLineStyle(2); h3->SetLineWidth(3); h3->SetFillStyle(0);
            h4->SetLineColor(kBlue + 3); h4->SetLineStyle(7); h4->SetLineWidth(3); h4->SetFillStyle(0);
            const double xlo = col == 0 ? 0.30 : 0.42, xhi = col == 0 ? 1.20 : 0.62;
            TH1D* fr = new TH1D(Nm("fr").c_str(), "", 100, xlo, xhi);
            fr->SetDirectory(nullptr);
            double ymax = 0;
            for (int b = 1; b <= hs[4]->GetNbinsX(); ++b)
                if (hs[4]->GetBinCenter(b) > xlo && hs[4]->GetBinCenter(b) < xhi) ymax = std::max(ymax, hs[4]->GetBinContent(b));
            if (col == 0) { fr->SetMinimum(0.05); fr->SetMaximum(ymax * 8); }
            else { fr->SetMinimum(0); fr->SetMaximum(ymax * 1.30); }
            fr->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]");
            fr->GetYaxis()->SetTitle(Form("expected rows / %d MeV/c^{2} at L_{sim}", rb));
            AxisStyle(fr, 0.055, 0.060, 1.05, 1.25);
            fr->GetXaxis()->SetNdivisions(col == 0 ? 509 : 505);
            fr->Draw("axis");
            for (int i = 4; i >= 0; --i) hs[i]->Draw("hist same");
            h3->Draw("hist same");
            h4->Draw("hist same");
            gPad->RedrawAxis();
            if (col == 1) {
                const double yl = ymax * 1.30;
                for (double xw : {kBuLo, kBuHi, kPkLo, kPkHi}) {
                    TLine* l = new TLine(xw, 0, xw, yl * 0.84);
                    l->SetLineStyle(2); l->SetLineWidth(3); l->SetLineColor(kBlack);
                    l->Draw();
                }
                Txt(0.33, 0.80, 0.050, "bump window");
                Txt(0.62, 0.80, 0.050, "peak");
            }
            Txt(0.17, 0.92, 0.055, Form("%s: %s", kSN[s], LatexCuts(gCut[s]).c_str()));
        }
    }
    Save(c, "s5_blind_prep_templates_eta");
    delete c;
}
void FigTemplatesX17() {
    TCanvas* c = new TCanvas("cTplX17", "", 2400, 1500);
    const int comps[3] = {cOther, cConv4, cDir5};
    const int fill[3] = {kGray, kAzure - 4, kBlue + 1};
    const char* lab[3] = {"other T1 + T2 (bump, physical background)", "#eta#rightarrow#pi^{+}#pi^{-}#gamma, #gamma converted (class 4)",
                          "direct #eta#rightarrow#pi^{+}#pi^{-}e^{+}e^{-} (class 5)"};
    TPad* hd = new TPad("hd", "", 0, 0.80, 1, 1);
    hd->Draw();
    hd->cd();
    Txt(0.01, 0.80, 0.150, Form("X17 templates from half A (#times 2): %s, %.4f #leq M_{4} < %.4f GeV/c^{2}", LatexCuts(gCutX17nomee).c_str(),
                                gXLo, gXHi));
    Txt(0.01, 0.60, 0.150, Form("rows at %s; dashed lines: frozen M(ee) window %.1f-%.1f MeV/c^{2}", kLdNote, gCut[sX17].mee_lo,
                                gCut[sX17].mee_hi));
    {
        TLegend* lg = new TLegend(0.01, 0.02, 0.99, 0.56);
        lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.140); lg->SetNColumns(2);
        for (int i = 2; i >= 0; --i) lg->AddEntry(Dummy(fill[i], fill[i], 1, 1), lab[i], "f");
        lg->AddEntry(Dummy(-1, kRed + 1, 2, 5), "X17 at BR = 10^{-5} (hypothesis, never in the sum)", "l");
        lg->Draw();
    }
    for (int col = 0; col < 2; ++col) {
        c->cd();
        TPad* p = new TPad(Nm("p").c_str(), "", col * 0.5, 0, col * 0.5 + 0.5, 0.80);
        p->SetLeftMargin(col == 0 ? 0.16 : 0.19); p->SetRightMargin(0.03); p->SetBottomMargin(0.15); p->SetTopMargin(0.04);
        if (col == 0) p->SetLogy();
        p->Draw();
        p->cd();
        if (col == 0) {
            std::vector<TH1D*> hs;   // cumulative
            double ymax = 0;
            for (int i = 0; i < 3; ++i) {
                TH1D* h = gMeeH[0][comps[i]].Hist(Nm("h").c_str(), 3, kLsim, kLdata);   // rows at L_data
                if (i) h->Add(hs[i - 1]);
                h->SetFillColor(fill[i]); h->SetFillStyle(1001); h->SetLineColor(fill[i]); h->SetLineWidth(1);
                hs.push_back(h);
            }
            TH1D* hx = gMeeH[0][cX17].Hist(Nm("h").c_str(), 3, kLsim, kLdata);
            hx->Scale(1e-5);
            hx->SetLineColor(kRed + 1); hx->SetLineStyle(2); hx->SetLineWidth(5); hx->SetFillStyle(0);
            for (int b = 1; b <= 120; ++b) ymax = std::max(ymax, hs[2]->GetBinContent(b));
            TH1D* fr = new TH1D(Nm("fr").c_str(), "", 100, 0, 60);
            fr->SetDirectory(nullptr);
            fr->SetMinimum(0.5); fr->SetMaximum(ymax * 3);
            fr->GetXaxis()->SetTitle("M(e^{+}e^{-}) [MeV/c^{2}]");
            fr->GetYaxis()->SetTitle("expected rows / 0.5 MeV/c^{2} at L_{data}");
            AxisStyle(fr, 0.050, 0.055, 1.05, 1.30);
            fr->Draw("axis");
            for (int i = 2; i >= 0; --i) hs[i]->Draw("hist same");
            hx->Draw("hist same");
            gPad->RedrawAxis();
            for (double xw : {gCut[sX17].mee_lo, gCut[sX17].mee_hi}) {
                TLine* l = new TLine(xw, 0.5, xw, ymax * 2);
                l->SetLineStyle(2); l->SetLineWidth(3);
                l->Draw();
            }
        } else {
            TH1D* h = gMeeH[0][cX17].Hist(Nm("h").c_str(), 0);
            const double a = h->Integral();
            if (a > 0) h->Scale(1. / a);
            h->SetLineColor(kRed + 1); h->SetLineWidth(3); h->SetMarkerColor(kRed + 1); h->SetMarkerStyle(20); h->SetMarkerSize(1.4);
            TH1D* fr = new TH1D(Nm("fr").c_str(), "", 100, 0, 60);
            fr->SetDirectory(nullptr);
            fr->SetMinimum(0); fr->SetMaximum(h->GetMaximum() * 1.25);
            fr->GetXaxis()->SetTitle("M(e^{+}e^{-}) [MeV/c^{2}]");
            fr->GetYaxis()->SetTitle("fraction / 0.5 MeV/c^{2}");
            AxisStyle(fr, 0.050, 0.055, 1.05, 1.65);
            fr->Draw("axis");
            h->Draw("e1 same");
            double fw = 0;
            for (int b = gMwBin0 + 1; b <= gMwBin1 + 1; ++b) fw += h->GetBinContent(b);
            for (double xw : {gCut[sX17].mee_lo, gCut[sX17].mee_hi}) {
                TLine* l = new TLine(xw, 0, xw, h->GetMaximum() * 1.1);
                l->SetLineStyle(2); l->SetLineWidth(3);
                l->Draw();
            }
            Txt(0.52, 0.88, 0.048, "X17 line, unit area");
            Txt(0.52, 0.82, 0.048, Form("in window: %.3f", fw));
            Txt(0.52, 0.76, 0.048, "m_{X17} = 17 MeV/c^{2}");
        }
    }
    Save(c, "s5_blind_prep_templates_x17");
    delete c;
}
void FigToys() {
    if (gNEta <= 0 && gNX17 <= 0) return;
    TCanvas* c = new TCanvas("cToys", "", 2400, 1850);
    TPad* hd = new TPad("hd", "", 0, 0.91, 1, 1);
    hd->Draw();
    hd->cd();
    Txt(0.01, 0.60, 0.27, Form("Pseudo-data checks, observable totals only: %d eta toys at L = L_{sim}/2 = 110.5 nb^{-1} (half B),", gNEta));
    Txt(0.01, 0.17, 0.27, Form("%d X17 toys at %s; blue: prediction of the Poisson bootstrap", gNX17, kLdNote));
    c->cd();
    struct Pn { int s, k; const char* t; };
    const Pn pn[3] = {{sPEAK, 3, "eta toys: PEAK selection, peak window"}, {sDIRECT, 3, "eta toys: DIRECT selection, peak window"},
                      {sBUMP, 1, "eta toys: BUMP selection, bump window"}};
    for (int i = 0; i < 4; ++i) {
        c->cd();
        const double x0 = (i % 2) * 0.5, y1 = 0.91 - (i / 2) * 0.455;
        TPad* p = new TPad(Nm("p").c_str(), "", x0, y1 - 0.455, x0 + 0.5, y1);
        p->SetLeftMargin(0.15); p->SetRightMargin(0.03); p->SetBottomMargin(0.16); p->SetTopMargin(0.09);
        p->Draw();
        p->cd();
        std::vector<double> v;
        const char* t = "";
        double pred = 0, vc = 0;
        if (i < 3) {
            if (gNEta <= 0) continue;
            v = gEtaTot.n[pn[i].s][pn[i].k];
            t = pn[i].t;
            pred = gMuTot[pn[i].s][pn[i].k]; vc = gVcoh[pn[i].s][pn[i].k];
        } else {
            if (gNX17 <= 0) continue;
            v = gX17obsWin;
            t = "X17 toys: rows in the frozen M_{4} #times M(ee) window";
        }
        const double mn = *std::min_element(v.begin(), v.end()), mx = *std::max_element(v.begin(), v.end());
        double lo = mn, hi = mx;
        if (i < 3) { lo = std::min(lo, pred - 4 * std::sqrt(vc)); hi = std::max(hi, pred + 4 * std::sqrt(vc)); }
        const double pad = 0.05 * (hi - lo) + 1;
        TH1D* h = new TH1D(Nm("h").c_str(), "", i < 3 ? 25 : 40, lo - pad, hi + pad);
        h->SetDirectory(nullptr);
        for (double x : v) h->Fill(x);
        h->SetLineColor(kBlack); h->SetLineWidth(3); h->SetFillColor(kGray); h->SetFillStyle(1001);
        h->SetMaximum(h->GetMaximum() * 1.75);
        h->GetXaxis()->SetTitle("rows per toy");
        h->GetYaxis()->SetTitle("toys");
        AxisStyle(h, 0.050, 0.055, 1.10, 1.10);
        h->GetXaxis()->SetNdivisions(506);
        h->Draw("hist");
        Txt(0.15, 0.93, 0.052, t);
        if (i < 3) {
            double m = 0, m2 = 0;
            for (double x : v) { m += x; m2 += x * x; }
            m /= v.size();
            const double var = v.size() > 1 ? (m2 - v.size() * m * m) / (v.size() - 1) : 0.;
            TF1* g = new TF1(Nm("g").c_str(), "gaus", lo - pad, hi + pad);
            g->SetNpx(400);
            g->SetParameters(v.size() * h->GetBinWidth(1) / std::sqrt(2 * M_PI * vc), pred, std::sqrt(vc));
            g->SetLineColor(kAzure + 2); g->SetLineWidth(4);
            g->Draw("same");
            Txt(0.18, 0.84, 0.046, Form("toys: mean %.1f, rms %.1f", m, std::sqrt(var)));
            Txt(0.18, 0.78, 0.046, Form("prediction: %.1f, #sqrt{V_{coh}} = %.1f", pred, std::sqrt(vc)), 11, kAzure + 2);
            Txt(0.18, 0.72, 0.046, Form("variance / V_{coh} = %.2f", var / vc));
        } else {
            Txt(0.18, 0.84, 0.046, "five injected BR values mixed;");
            Txt(0.18, 0.78, 0.046, "the mapping is sealed");
        }
    }
    Save(c, "s5_blind_prep_toys_check");
    delete c;
}

// ============================================================================
// main
// ============================================================================
ULong64_t SplitMix(ULong64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}
void Run(ULong64_t seedIn, int nEta, int nX17) {
    TStopwatch sw;
    Style();
    if (!SetupCuts()) return;
    gSystem->mkdir(BDir("templates").c_str(), kTRUE);
    gSystem->mkdir(BDir("pseudo").c_str(), kTRUE);
    gSystem->mkdir(BDir("sealed").c_str(), kTRUE);
    // secret seed (sealed only)
    ULong64_t seed = seedIn;
    if (seed == 1) {   // reproduce: read the master seed from the sealed file (so that it never appears in a command line or log)
        std::ifstream in(BDir("sealed") + "/seed.txt");
        std::string line;
        seed = 0;
        while (std::getline(in, line)) {
            const size_t p = line.find("seed = ");
            if (p != std::string::npos && line[0] != '#') { seed = std::strtoull(line.c_str() + p + 7, nullptr, 10); break; }
        }
        if (seed <= 1) { printf("s5b: no valid seed in %s/seed.txt: stop\n", BDir("sealed").c_str()); return; }
        printf("s5b: master seed read from the sealed directory (reproduction run)\n");
    }
    if (!seed) { std::random_device rd; seed = ((ULong64_t)rd() << 32) ^ (ULong64_t)rd(); if (seed <= 1) seed = 2; }
    gK = fopen((BDir("sealed") + "/key_summary.txt").c_str(), "w");
    {
        FILE* fs = fopen((BDir("sealed") + "/seed.txt").c_str(), "w");
        const time_t now = time(nullptr);
        std::string tm = ctime(&now);
        while (!tm.empty() && (tm.back() == '\n' || tm.back() == '\r')) tm.pop_back();
        fprintf(fs, "# SEALED: master seed of s5_blind_prep.C (run %s)\n", tm.c_str());
        fprintf(fs, "seed = %llu\n# reproduce: s5_blind_prep.C+(1, %d, %d) reads this file (or pass the value as the first argument)\n",
                (unsigned long long)seed, nEta, nX17);
        fclose(fs);
    }
    PK("# SEALED KEY of the blind test (s5_blind_prep.C). Do not show to the analyst before the extraction is frozen.\n");
    PK("# eta toys: %d at L_sim / 2; X17 toys: %d per BR at L_data = 6080.5 nb^-1 (UNCONFIRMED)\n", nEta, nX17);
    // physics weights and rows
    if (!simphys::Init()) { printf("simphys::Init failed\n"); return; }
    simphys::Install();
    if (!simphys::LoadAngles()) { printf("simphys::LoadAngles failed\n"); return; }
    gO = SimOpt();
    gO.phys = simphys::kPhysMain;   // mode A, w (PDG), sets 7, e+e-
    gComps = MakeComps();
    gD = new SimData;
    if (!Load(*gD, kCompact, kMaskUS, kM4Lo, kM4Hi, 1)) return;
    printf("loaded in %.0f s\n", sw.RealTime());
    sw.Continue();
    // part 1
    PartRegions();
    PartTemplates();
    TablesPart1();
    gT = fopen((std::string(kTabDir) + "/s5_blind_prep_checks.txt").c_str(), "w");
    PT("# s5_blind_prep_checks: checks of the blind-test preparation (s5_blind_prep.C); simulation only; public part (no sealed numbers)\n");
    PT("# rows: std + lep + eta owned rows, mode A, w (PDG) x physics weights V3; pre; T1 classes 2-5 + T2 (X17, T3, T4 not in the pseudo-data)\n\n");
    ChecksPart1();
    SealedHalves();
    // part 2 and 4 (independent random streams derived from the master seed)
    PartEtaToys(SplitMix(seed ^ 0x1111111111111111ULL), nEta);
    printf("angles: hit %lld miss %lld\n", simphys::AngHit(), simphys::AngMiss());
    PartX17Toys(SplitMix(seed ^ 0x2222222222222222ULL), nX17);
    ChecksToys();
    PT("## 6. files\n");
    for (const char* fnm : {"templates/templates.root", "templates/README.md", "pseudo/eta_toys.root", "pseudo/x17_toys.root"}) {
        Long_t id, flags, mt; Long64_t size = 0;
        gSystem->GetPathInfo((std::string(kBlind) + "/" + fnm).c_str(), &id, &size, &flags, &mt);
        PT("%s/%s  %.1f MB\n", kBlind, fnm, size / 1048576.);
    }
    PT("sealed directory: %s/sealed (key_summary.txt, eta_key.root, x17_key.root, seed.txt)\n", kBlind);
    PT("physics-weight angles (V3): hit %lld, miss %lld\n", simphys::AngHit(), simphys::AngMiss());
    fclose(gT);
    gT = nullptr;
    fclose(gK);
    gK = nullptr;
    WriteReadmes(nEta, nX17);
    // figures
    FigTemplatesEta();
    FigTemplatesX17();
    FigToys();
    printf("s5_blind_prep done in %.0f s\n", sw.RealTime());
}

}  // namespace s5b

void s5_blind_prep(ULong64_t seed = 0, int nEta = 200, int nX17 = 300) { s5b::Run(seed, nEta, nX17); }
