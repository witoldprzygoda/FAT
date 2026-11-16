// ============================================================================
// defs.h  --  frozen definitions of SPEC report_all (Sect. 2), namespace alldefs
// ============================================================================
// HADES p+p, T = 4.5 GeV (feb22), p p -> pi+ pi- e+ e- X, M4 = M(pi+pi-e+e-).
// Written by WP0 (2026-11-14). Every agent of the report reads its numbers from
// here and its histograms from hists_v1.root (naming scheme: section N below and
// plots/all_report/README.md). Never change a value here without a new version
// of hists_v1.root.
//
// Units: masses GeV/c2 (M(ee) also GeV/c2 in code; MeV only in names/labels),
// momenta MeV/c, angles deg, lengths mm, cross sections nb, L nb^-1.
//
// Contents
//   A. normalisation (L, L_sim, 64/63, live time, old-table factor, conversions)
//   B. beam energies, periods (days -> I / IIA / IIB / III)
//   C. mass transfer of the simulation (s, r, M(pipi) scale) and the edge E(m, q)
//   D. selections pre, S0 and variants (data and simulation rows)
//   E. layouts F10, F5g, X2D, G3, BEFORE, fine
//   F. windows and control regions
//   G. CB construction parameters
//   H. simulation components (classes C2-C6, groups G_*, T3, T4) and variants
//   N. names in hists_v1.root
// ============================================================================
#ifndef ALLREPORT_DEFS_H
#define ALLREPORT_DEFS_H

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "TString.h"

#include "SimEngine.h"      // simeng:: (simulation rows) and, through it, ../ScanEngine.h (scan::, data rows)
#include "sim_ghostvar.h"   // simeng::GhostPair (ghost variant (a))

namespace alldefs {

// ============================================================================
// A. Normalisation (spec D7, D16, D17, WP2)
// ============================================================================
const double kL        = 6079.2;        // nb^-1, analysed runs (17 557 runs present in the data; START-free chain "N";
                                        // LH2 only, i.e. yields need the empty-target subtraction f_win of WP2).
                                        // [v0 verification 2026-11-14: wording made explicit] Sector-2 hardware trigger
                                        // deficit of days 034-044: the elastic normalisation removes it (factor S2fac
                                        // of the per-day CSV), so L itself is unbiased; the same deficit still reduces
                                        // the analysis' own PT3 yields on those days and must be corrected by C_sec
                                        // (spec D16, WP2; A3 note: "NOT in L: handle it in the analysis' own trigger
                                        // efficiency"). Nothing in this file applies C_sec.
const double kLrelErr  = 0.069;         // total relative uncertainty of L (HADES syst 5.56 %, sigma_int 1.03 %, Argonne 4.0 %)
const double kLsim     = 221.05;        // nb^-1 per simulated set (std, lep, eta describe the same 1e10 SMASH events)
const double kPT3      = 64.0 / 63.0;   // data yields x 64/63 (1/64 of the PT3 events carry the word 4096; FAT selects 8192)
const double kLiveUp   = 0.073;         // live time: one-sided +7.3 % / -0 % on cross sections, NOT applied (D17)
const double kOldTab   = 1.0161;        // old tables (L 6080.5, no 64/63) x 1.0161 = this normalisation
// conversions
const double kNbToRaw  = kL / kPT3;             // sigma [nb] -> expected RAW data rows (x L x 63/64) = 5984.20
const double kSimToRaw = kL / kLsim / kPT3;     // sum w      -> expected RAW data rows = 27.0718
const double kRawToNb  = kPT3 / kL;             // raw data rows (N - CB) -> nb  (x 64/63 / L)
// PDG / SMASH numbers used by later work packages (A2 note, sim report)
const double kBRdir    = 2.68e-4;       // BR(eta -> pi+pi-e+e-), PDG (+- 0.11e-4)
const double kSigEtaSMASH = 1.0665;     // mb, SMASH sigma(pp -> eta X)

// ============================================================================
// B. Beam and periods
// ============================================================================
const double kTdata = 4530.;            // MeV, data analysis (FAT config)
const double kTsim  = 4500.;            // MeV, SMASH generation (sim MM4 and y4 are recomputed with this value)
// periods of the feb22 beam time by day (= file_id of the data compact rows)
//   I = 034-046 (13 days), IIA = 047-048, IIB = 049-053, III = 059-067 (days 054-058 absent)
enum Period { kPerI = 0, kPerIIA, kPerIIB, kPerIII, kNPer };
inline int PeriodOfDay(int day) {
    if (day >= 34 && day <= 46) return kPerI;
    if (day >= 47 && day <= 48) return kPerIIA;
    if (day >= 49 && day <= 53) return kPerIIB;
    if (day >= 59 && day <= 67) return kPerIII;
    return -1;
}
inline const char* PeriodName(int p) { static const char* n[kNPer] = {"I", "IIA", "IIB", "III"}; return (p >= 0 && p < kNPer) ? n[p] : "?"; }
// L6_SF per period (sf_luminosity_perperiod.csv, all runs; the analysed-run L is 6079.2 of 6080.52 in total)
const double kLper[kNPer] = {2810.669, 372.610, 1078.134, 1819.110};

// ============================================================================
// C. Mass transfer (D8) and the bump edge E(m, q) (Sect. 2)
// ============================================================================
const double kS      = 0.99615;   // sim M4' = s [m4_true + r (m4 - m4_true)]
const double kR      = 1.0;
const double kSPiPi  = 0.9952;    // sim M(pipi) x 0.9952 (K0S: sim 1.0018 vs data 0.9970); M(ee) not scaled
const double kSK0    = 0.99521;   // variant s = K0S-based value
const double kSerr   = 0.0009;    // variant s +- 0.0009
const double kRlo = 0.89, kRhi = 1.10;
// rows without a true four-particle mass (m4_true <= 0: a misidentified track) get M4' = s m4 for every r
inline double M4sim(double m4, double m4true, double s = kS, double r = kR) {
    if (!(m4true > 0.)) return s * m4;
    return s * (m4true + r * (m4 - m4true));
}
// edge of the bump on the DATA scale for m = M(pipi) (data scale) and q = M(ee) [GeV/c2]:
//   E(m, q) = 0.99274 sqrt(m_eta^2 - m_eta (E_pi0 - p_pi0)(1 - q^2/m_pi0^2)),
//   E_pi0 = (m_eta^2 + m_pi0^2 - m'^2)/(2 m_eta), p_pi0 = sqrt(E_pi0^2 - m_pi0^2), m' = m/0.9970
// returns -1 above the limit (m > 0.4116: no bump)
const double kEdgeScale = 0.99274, kEdgePiPi = 0.9970;
inline double Edge(double m, double q = 0.) {
    const double me = scan::kMeta, mp0 = scan::kMpi0;
    const double mpr = m / kEdgePiPi;
    const double E0 = (me * me + mp0 * mp0 - mpr * mpr) / (2. * me);
    const double p2 = E0 * E0 - mp0 * mp0;
    if (!(p2 >= 0.)) return -1.;
    const double a = me * me - me * (E0 - std::sqrt(p2)) * (1. - q * q / (mp0 * mp0));
    return a > 0. ? kEdgeScale * std::sqrt(a) : -1.;
}
const double kEdgeLimit = 0.4116;   // data-scale M(pipi) above which there is no bump

// ============================================================================
// D. Selections (Sect. 2, D1, WP0 variants)
// ============================================================================
// pre = isBest >= 0 && vz > -500 && PT3 && vz > -160 && cor_ok == 1   (simeng::PassPre; data PT3 = production word 8192,
//       sim PT3 = (trigbit & 8192) != 0)
// S0  = pre && pip_isUsedVertex == 1 && pim_isUsedVertex == 1          (NOMINAL for bump, peak and X17)
// variants: dz40 = pre && dz_pi < 40 mm (DIRECT); PEAK = pre && dz_pi < 40 && MM4 < 2.85 GeV/c2 (sim MM4 at 4500 MeV);
//           vz150 / vz170 = S0 with vz > -150 / vz > -170 instead of -160; none = pre without vz > -160 (data check);
//           isBest == 1 (data only); cor_ok fallback (data only: cor_ok == 0 rows with the reconstructed pi-, *_fb values)
// dz_pi = max(|pip_z - vz|, |pim_z - vz|) (float arithmetic, = sim study DzPi; data column dz_pi of scan_v2)
enum SelId { kSelS0 = 0, kSelPre, kSelDz40, kSelPeak, kSelVz150, kSelVz170, kSelNone, kNSel };
inline const char* SelTag(int s) {
    static const char* n[kNSel] = {"S0", "pre", "dz40", "PEAK", "vz150", "vz170", "none"};
    return (s >= 0 && s < kNSel) ? n[s] : "?";
}
inline const char* SelText(int s) {
    static const char* n[kNSel] = {
        "S0 = pre + both pions in the vertex fit",
        "pre = isBest >= 0, vz > -500, PT3, vz > -160 mm, cor_ok",
        "pre + dz_pi < 40 mm (DIRECT)",
        "PEAK = pre + dz_pi < 40 mm + MM4 < 2.85 GeV/c2 (sim MM4 at 4500 MeV)",
        "S0 with vz > -150 mm",
        "S0 with vz > -170 mm",
        "none = isBest >= 0, vz > -500, PT3, cor_ok (no vz > -160, no pion-vertex cut)"};
    return (s >= 0 && s < kNSel) ? n[s] : "?";
}
inline float DzPi(float pipz, float pimz, float vz) { return std::max(std::fabs(pipz - vz), std::fabs(pimz - vz)); }
// base without the vz > -160 cut: isBest >= 0, vz > -500, PT3, cor_ok (corOk = cor_ok == 1 or the fallback)
template <class R> inline bool PassBase(const R& r, bool corOk) {
    return r.isBest >= 0 && r.vz > -500.f && corOk && simeng::PT3(r);
}
template <class R> inline bool PiVtx(const R& r) { return r.pip_isUsedVertex == 1 && r.pim_isUsedVertex == 1; }
// selection s for a row; dzpi and mm4 are passed explicitly (sim: dz_pi computed, MM4 at 4500 MeV)
template <class R> inline bool PassSel(int s, const R& r, bool corOk, float dzpi, float mm4) {
    if (!PassBase(r, corOk)) return false;
    switch (s) {
        case kSelS0:    return r.vz > -160.f && PiVtx(r);
        case kSelPre:   return r.vz > -160.f;
        case kSelDz40:  return r.vz > -160.f && dzpi < 40.f;
        case kSelPeak:  return r.vz > -160.f && dzpi < 40.f && mm4 < 2.85f;
        case kSelVz150: return r.vz > -150.f && PiVtx(r);
        case kSelVz170: return r.vz > -170.f && PiVtx(r);
        case kSelNone:  return true;
    }
    return false;
}

// ============================================================================
// E. Layouts
// ============================================================================
// M(pipi) slices (data scale): S1 < 0.34, S2 0.34-0.37, S3 0.37-0.39, S4 0.39-0.43, S5 >= 0.43
const int kNSlice = 5;
const double kSliceEdge[kNSlice - 1] = {0.34, 0.37, 0.39, 0.43};
inline int Slice(double mpipi) {
    if (!(mpipi == mpipi)) return -1;
    int j = 0;
    while (j < kNSlice - 1 && mpipi >= kSliceEdge[j]) ++j;
    return j;
}
// M(ee) classes: E1 < 10 MeV, E2 >= 10 MeV
const double kEsplit = 0.010;
inline int Eclass(double mee) { return !(mee == mee) ? -1 : (mee < kEsplit ? 0 : 1); }
const double kF5gMee = 0.020;   // F5g: M(ee) >= 20 MeV only
// F10 / F5g M4 binning: 5 MeV in 0.400-0.520, 2.5 MeV in 0.520-0.570, 5 MeV in 0.570-0.700 (70 bins)
inline std::vector<double> F10Edges() {
    std::vector<double> e;
    for (int i = 0; i <= 24; ++i) e.push_back(0.400 + 0.005 * i);          // 0.400 ... 0.520
    for (int i = 1; i <= 20; ++i) e.push_back(0.520 + 0.0025 * i);         // ... 0.570
    for (int i = 1; i <= 26; ++i) e.push_back(0.570 + 0.005 * i);          // ... 0.700
    for (auto& x : e) x = std::round(x * 1e6) / 1e6;
    return e;
}
const int kNF10 = 70;
const double kF10Lo = 0.400, kF10Hi = 0.700;
// X2D: M4 0.5225-0.6100 in 2.5 MeV (35) x M(ee) 2-62 MeV in 2 MeV (30)
const int kX2Dnx = 35, kX2Dny = 30;
const double kX2Dxlo = 0.5225, kX2Dxhi = 0.6100, kX2Dylo = 0.002, kX2Dyhi = 0.062;
// G3 (superset of F10/F5g for the binning, range, slice-edge, M(ee)-split, CR-K and K0S-veto variants of WP3/WP5/WP7):
//   M4 0.400-0.850 in 2.5 MeV (180) x M(pipi) 27 bins x M(ee) 5 classes
//   M(pipi) edges: 0, then 5 MeV from 0.330 to 0.445, then 0.485, 0.507 (K0S veto 0.485-0.507), 2.0
//   M(ee) edges [GeV/c2]: 0, 0.008, 0.010, 0.012, 0.020, 2.0
const int kG3nx = 180;
const double kG3xlo = 0.400, kG3xhi = 0.850;
inline std::vector<double> G3PiPiEdges() {
    std::vector<double> e{0.0};
    for (int i = 0; i <= 23; ++i) e.push_back(std::round((0.330 + 0.005 * i) * 1e6) / 1e6);   // 0.330 ... 0.445
    e.push_back(0.485); e.push_back(0.507); e.push_back(2.0);
    return e;
}
inline std::vector<double> G3MeeEdges() { return {0.0, 0.008, 0.010, 0.012, 0.020, 2.0}; }
// fine master M4 grid: 0.30-1.10 in 2.5 MeV (320 bins); BEFORE M4: 0.30-1.10 in 5 MeV (160)
const int kNfine = 320;
const double kFineLo = 0.30, kFineHi = 1.10, kFineW = 0.0025;
inline int FineBin(double m4) { return (m4 >= kFineLo && m4 < kFineHi) ? std::min(kNfine - 1, (int)((m4 - kFineLo) / kFineW)) : -1; }
const int kNbefM4 = 160;
// BEFORE variables in the 8 M4 regions (Table B1 regions, below)
enum BVar { kBmpipi = 0, kBmee, kBoa, kBvz, kBmm4, kBpt4, kBy4, kNBVar };
inline const char* BVarTag(int v) { static const char* n[kNBVar] = {"mpipi", "mee", "oa", "vz", "mm4", "pt4", "y4"}; return n[v]; }
inline const char* BVarTitle(int v) {
    static const char* n[kNBVar] = {"M(#pi^{+}#pi^{-}) [GeV/c^{2}]", "M(e^{+}e^{-}) [GeV/c^{2}]", "OA(e^{+}e^{-}) [deg]",
                                    "v_{z} [mm]", "MM(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]", "p_{T}(#pi^{+}#pi^{-}e^{+}e^{-}) [MeV/c]",
                                    "y_{cm}(#pi^{+}#pi^{-}e^{+}e^{-})"};
    return n[v];
}
const int    kBVnb[kNBVar] = {320, 250, 180, 200, 180, 60, 80};
const double kBVlo[kNBVar] = {0.25, 0.0, 0.0, -320., 1.6, 0., -2.};
const double kBVhi[kNBVar] = {1.05, 0.5, 180., 80., 3.4, 1200., 2.};
// BEFORE map: M4 0.30-1.10 in 10 MeV (80) x M(pipi) 0.25-1.05 in 10 MeV (80)
const int kMapNx = 80, kMapNy = 80;
const double kMapXlo = 0.30, kMapXhi = 1.10, kMapYlo = 0.25, kMapYhi = 1.05;

// ============================================================================
// F. Windows and control regions (data scale; simulation rows use M4')
// ============================================================================
// bump B 0.4550-0.5225; peak P 0.53364-0.55416 (= 0.5439 +- 1.9 x 5.4 MeV); side bands 0.5250-0.5336, 0.5542-0.5628;
// X17 counting cross-check M4 0.5312-0.5566 (= sim 0.5333-0.5587 x 0.99615) and M(ee) 15.5-24.5 MeV;
// CR-U 0.565 <= M4 < 0.700 all slices; CR-K S5 and 0.400 <= M4 < 0.5225 (spec: "< 0.522"; on the F10 grid the
// nearest bin edge is 0.5225, see the WP0 note); CR-L 0.5225-0.5300 (test only);
// Table B1 regions R0 < 0.455, R1 = B, R2 0.5225-0.53364, R3 = P, R4 0.55416-0.565, R5 0.565-0.70, R6 0.70-0.90,
// R7 0.90-1.10; jump J = R(0.560-0.595) / R(0.4875-0.5225); X17P = P and M(ee) 15.5-24.5 MeV (the count 5258 of D19 uses this
// window, not X17c).
const double kBlo = 0.4550, kBhi = 0.5225;
const double kPlo = 0.53364, kPhi = 0.55416;
struct Window { const char* tag; double lo, hi; double meeLo, meeHi; int slice; const char* text; };
// meeLo/meeHi < 0: no M(ee) condition; slice < 0: all slices
enum WinId { kWB = 0, kWP, kWSBlo, kWSBhi, kWX17c, kWCRU, kWCRK, kWCRL, kWR0, kWR1, kWR2, kWR3, kWR4, kWR5, kWR6, kWR7,
             kWJlo, kWJhi, kWFit, kWBef, kWX17P, kNWin };
const Window kWin[kNWin] = {
    {"B", 0.4550, 0.5225, -1, -1, -1, "bump window"},
    {"P", 0.53364, 0.55416, -1, -1, -1, "peak window"},
    {"SBlo", 0.5250, 0.5336, -1, -1, -1, "peak lower side band"},
    {"SBhi", 0.5542, 0.5628, -1, -1, -1, "peak upper side band"},
    {"X17c", 0.5312, 0.5566, 0.0155, 0.0245, -1, "X17 counting window, M(ee) 15.5-24.5 MeV"},
    {"CRU", 0.565, 0.700, -1, -1, -1, "CR-U"},
    {"CRK", 0.400, 0.5225, -1, -1, 4, "CR-K (S5)"},
    {"CRL", 0.5225, 0.5300, -1, -1, -1, "CR-L (test only)"},
    {"R0", 0.0, 0.455, -1, -1, -1, "M4 < 0.455"},
    {"R1", 0.4550, 0.5225, -1, -1, -1, "B"},
    {"R2", 0.5225, 0.53364, -1, -1, -1, "0.5225-0.5336"},
    {"R3", 0.53364, 0.55416, -1, -1, -1, "P"},
    {"R4", 0.55416, 0.565, -1, -1, -1, "0.5542-0.565"},
    {"R5", 0.565, 0.700, -1, -1, -1, "0.565-0.70"},
    {"R6", 0.700, 0.900, -1, -1, -1, "0.70-0.90"},
    {"R7", 0.900, 1.100, -1, -1, -1, "0.90-1.10"},
    {"Jlo", 0.4875, 0.5225, -1, -1, -1, "jump denominator band"},
    {"Jhi", 0.560, 0.595, -1, -1, -1, "jump numerator band"},
    {"Fit", 0.400, 0.700, -1, -1, -1, "F10 range"},
    {"Bef", 0.300, 1.100, -1, -1, -1, "BEFORE range"},
    {"X17P", 0.53364, 0.55416, 0.0155, 0.0245, -1, "peak window P and M(ee) 15.5-24.5 MeV (the window of the count 5258 quoted in D19)"}};
const int kNReg = 8;   // the Table B1 regions R0..R7 = windows kWR0..kWR7 (BEFORE variables are filled per region)
inline int Region(double m4) {
    for (int k = 0; k < kNReg; ++k) if (m4 >= kWin[kWR0 + k].lo && m4 < kWin[kWR0 + k].hi) return k;
    return -1;
}
inline bool InWindow(int w, double m4, double mee, double mpipi) {
    const Window& x = kWin[w];
    if (!(m4 >= x.lo && m4 < x.hi)) return false;
    if (x.meeLo >= 0 && !(mee >= x.meeLo && mee < x.meeHi)) return false;
    if (x.slice >= 0 && Slice(mpipi) != x.slice) return false;
    return true;
}
// control-region membership of an F10 bin (bin given by its lower and upper M4 edge and its slice index)
inline bool InCRU(double lo, double hi, int) { return lo >= 0.565 - 1e-9 && hi <= 0.700 + 1e-9; }
inline bool InCRK(double lo, double hi, int slice) { return slice == 4 && lo >= 0.400 - 1e-9 && hi <= 0.5225 + 1e-9; }
inline bool InCRL(double lo, double hi, int) { return lo >= 0.5225 - 1e-9 && hi <= 0.5300 + 1e-9; }
// tilt variable of the background model (D11)
inline double TiltX(double m4) { return (m4 - 0.55) / 0.1; }

// ============================================================================
// G. CB construction (D5, Sect. 2)
// ============================================================================
// CB = k 2 sqrt(N++ N--) in 10 MeV M4 cells (4 fine bins of 2.5 MeV on the grid from 0.30), linearly interpolated to
// the centres of the fine bins exactly as ScanEngine.h (constant beyond the first/last cell centre), k = 1. The CB of a
// cell of any layout (M(pipi) slice x M(ee) class, M(ee) bin, ...) = CB(fine bin) x f(cell | 50 MeV band), f = the
// fraction of the like-sign rows (N++ + N--) of that 50 MeV M4 band (bands 0.30-0.35, 0.35-0.40, ...) that fall in the
// cell. A window or region not aligned to the fine grid takes fractional fine bins (as ScanEngine WindowApprox).
// BEFORE variables of a region: CB(region) x f(variable bin | region). Error: 200 Poisson event-group replicas of the
// like-sign rows (each like-sign event group gets a Poisson(1) weight), which propagate both the 10 MeV cell counts and f.
// Variant CB20: the same with 20 MeV cells.
const int kCBreb = 4;          // fine bins per CB cell (10 MeV)
const int kCBreb20 = 8;        // 20 MeV variant
const double kBandW = 0.050;   // GeV/c2
const int kNBand = 16;         // 0.30-1.10
inline int BandOfFine(int i) { return i < 0 ? -1 : std::min(kNBand - 1, i / 20); }   // 20 fine bins per band
inline int Band(double m4) { return BandOfFine(FineBin(m4)); }
const int kNRepCB = 200;
const int kNRepSim = 100;

// ============================================================================
// H. Simulation components and variants (D6, WP0)
// ============================================================================
// classes: C2 bump pi0 Dalitz (T1, cls 2), C3 bump conversion (cls 3), C4 peak conversion eta -> pi+pi-gamma (cls 4),
// C5 direct eta -> pi+pi-e+e- (cls 5, eta set, V3), C6 X17 (cls 6, PER UNIT BR: weight 2 per row, never in sums);
// groups of T2 (first match): G_eta (pion class P2 or lepton class L2), G_pi0D (L1), G_conv (L3), G_vec (L4-L6),
// G_oth (rest); "_om" = both pions from one omega or rho0 (pcls 3 or 4), "_nom" = the rest; G_D = G_pi0D + G_vec + G_oth,
// G_C = G_conv; T12 = C2-C5 + T2 (the simulated sum, never X17); T3 fake lepton, T4 different sources (shown apart).
// Ghost variant (a): ghost-flagged true pairs (sim_ghostvar.h GhostPair) move from T3 into T2 (and its groups).
enum Comp {
    kC2 = 0, kC3, kC4, kC5, kC6, kGeta, kGpi0D, kGpi0Dom, kGpi0Dnom, kGconv, kGconvom, kGconvnom, kGvec, kGvecom, kGvecnom,
    kGoth, kGothom, kGothnom, kGD, kT12, kT3, kT4, kNComp };
inline const char* CompTag(int c) {
    static const char* n[kNComp] = {"C2", "C3", "C4", "C5", "C6", "Geta", "Gpi0D", "Gpi0D_om", "Gpi0D_nom", "Gconv",
                                    "Gconv_om", "Gconv_nom", "Gvec", "Gvec_om", "Gvec_nom", "Goth", "Goth_om", "Goth_nom",
                                    "GD", "T12", "T3", "T4"};
    return (c >= 0 && c < kNComp) ? n[c] : "?";
}
inline const char* CompText(int c) {
    static const char* n[kNComp] = {
        "C2 bump, pi0 Dalitz (cls 2)", "C3 bump, conversion (cls 3)", "C4 peak, eta -> pi+pi-gamma conversion (cls 4)",
        "C5 direct eta -> pi+pi-e+e- (cls 5)", "C6 X17 per unit BR (cls 6, weight 2)", "G_eta (T2, P2 or L2)",
        "G_pi0D (T2, L1)", "G_pi0D, omega/rho0 pions", "G_pi0D, other pions", "G_conv (T2, L3)", "G_conv, omega/rho0 pions",
        "G_conv, other pions", "G_vec (T2, L4-L6)", "G_vec, omega/rho0 pions", "G_vec, other pions", "G_oth (T2, rest)",
        "G_oth, omega/rho0 pions", "G_oth, other pions", "G_D = G_pi0D + G_vec + G_oth", "T12 = C2-C5 + T2 (simulated sum)",
        "T3 fake lepton", "T4 different sources"};
    return (c >= 0 && c < kNComp) ? n[c] : "?";
}
// group of a template-2 row (s5_blind_prep.C BkgGroup): 4 G_eta, 5 G_pi0D, 6 G_conv, 7 G_vec, 8 G_oth
inline int BkgGroup(const simeng::Row& r) {
    if (r.pcls == 2 || r.lcls == 2) return 4;
    if (r.lcls == 1) return 5;
    if (r.lcls == 3) return 6;
    if (r.lcls >= 4 && r.lcls <= 6) return 7;
    return 8;
}
inline bool OmegaPions(const simeng::Row& r) { return r.pcls == 3 || r.pcls == 4; }
// components of a row (ghostA: variant (a)); returns the number written to out (at most 6)
inline int CompsOf(const simeng::Row& r, bool ghostA, int* out) {
    int n = 0;
    int t = r.tmpl;
    if (ghostA && simeng::GhostPair(r)) t = 2;
    if (t == 1) {
        if (r.cls >= 2 && r.cls <= 5) { out[n++] = kC2 + (r.cls - 2); out[n++] = kT12; }
        else if (r.cls == 6) out[n++] = kC6;
        return n;
    }
    if (t == 2) {
        const int g = BkgGroup(r);
        const bool om = OmegaPions(r);
        switch (g) {
            case 4: out[n++] = kGeta; break;
            case 5: out[n++] = kGpi0D; out[n++] = om ? kGpi0Dom : kGpi0Dnom; out[n++] = kGD; break;
            case 6: out[n++] = kGconv; out[n++] = om ? kGconvom : kGconvnom; break;
            case 7: out[n++] = kGvec; out[n++] = om ? kGvecom : kGvecnom; out[n++] = kGD; break;
            default: out[n++] = kGoth; out[n++] = om ? kGothom : kGothnom; out[n++] = kGD; break;
        }
        out[n++] = kT12;
        return n;
    }
    if (t == 3) { out[n++] = kT3; return n; }
    if (t == 4) { out[n++] = kT4; return n; }
    return n;
}
// physics-weight codes (sim_eta_physics.h): 108 nominal V3 + KLOE + conv F_V + bkg rules
const int kPhysNom = 108, kPhysV2 = 107, kPhysBump1aY = 116, kPhysBumpPS = 100, kPhysConvOff = 76, kPhysGS = 236, kPhysOff = 0;
// simulation variants (directory suffix)
enum SimVar { kVnom = 0, kVV2, kVbump1aY, kVbumpPS, kVconvOff, kVGS, kVwtsOff, kVmodeB, kVwdb, kVghostA, kVsK0, kVsUp, kVsDn,
              kVr089, kVr110, kNSimVar };
inline const char* SimVarTag(int v) {
    static const char* n[kNSimVar] = {"nom", "V2", "bump1aY", "bumpPS", "convOff", "GS", "wtsOff", "modeB", "wdb", "ghostA",
                                      "sK0", "sUp", "sDn", "r089", "r110"};
    return (v >= 0 && v < kNSimVar) ? n[v] : "?";
}
inline const char* SimVarText(int v) {
    static const char* n[kNSimVar] = {
        "nominal: mode A, w x V3 (108), s = 0.99615, r = 1", "direct V2 (107)", "bump 1+aY (116)", "bump phase space (100)",
        "conversion F_V off (76)", "F_V Gounaris-Sakurai (236)", "physics weights off (0)", "ownership mode B", "weight w_db",
        "ghost variant (a): ghost-flagged true pairs in T2", "s = 0.99521", "s = 0.99705", "s = 0.99525", "r = 0.89", "r = 1.10"};
    return (v >= 0 && v < kNSimVar) ? n[v] : "?";
}
inline int SimVarPhys(int v) {
    switch (v) {
        case kVV2: return kPhysV2; case kVbump1aY: return kPhysBump1aY; case kVbumpPS: return kPhysBumpPS;
        case kVconvOff: return kPhysConvOff; case kVGS: return kPhysGS; case kVwtsOff: return kPhysOff;
    }
    return kPhysNom;
}
inline double SimVarS(int v) { return v == kVsK0 ? kSK0 : v == kVsUp ? kS + kSerr : v == kVsDn ? kS - kSerr : kS; }
inline double SimVarR(int v) { return v == kVr089 ? kRlo : v == kVr110 ? kRhi : kR; }

// ============================================================================
// N. Names in hists_v1.root (see README.md)
// ============================================================================
const char* const kHistFile   = "/hdd3/przygoda/HADES/FAT_out/all_report/hists_v1.root";
const char* const kSealedFile = "/hdd3/przygoda/HADES/FAT_out/all_report/sealed/hists_v1_x17_sealed.root";
// layout histogram names (suffix after "<what>_" or "<comp>_")
inline std::string F10Name(int slice, int ecl) { return Form("F10_S%dE%d", slice + 1, ecl + 1); }
inline std::string F5gName(int slice) { return Form("F5g_S%d", slice + 1); }
inline std::string BVarName(int v, int reg) { return Form("B_%s_R%d", BVarTag(v), reg); }
// data subsets
enum DataSub { kSubAll = 0, kSubFold0, kSubFold1, kSubPerI, kSubPerIIA, kSubPerIIB, kSubPerIII, kNSub };
inline const char* SubTag(int s) {
    static const char* n[kNSub] = {"all", "fold0", "fold1", "perI", "perIIA", "perIIB", "perIII"};
    return (s >= 0 && s < kNSub) ? n[s] : "?";
}

}  // namespace alldefs

#endif  // ALLREPORT_DEFS_H
