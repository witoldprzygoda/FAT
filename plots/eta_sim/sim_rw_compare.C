// ============================================================================
// sim_rw_compare.C  --  eta peak and eta-enriched comparisons with the physics-corrected templates
// ============================================================================
// HADES p+p, T = 4.53 GeV, pi+ pi- e+ e- (e+e- hypothesis). The eta decays of the SMASH simulation
// are reweighted to physics models (sim_eta_physics.h, built by sim_rw_build.C; nothing is produced
// again): direct eta -> pi+pi-e+e- (eta set, PLUTO plain four-body, M(ee) not QED) to Kroll-Wada x
// p*^3 lambda^{3/2} x |F_V|^2 x |F_VMD|^2 x e+ and pi+ decay angles (main variant V3 since the review of
// 2026-11-12; V2 without the angles was the main variant before), the bump eta -> pi+pi-pi0 (HGeant, phase
// space ASSUMED) to the KLOE-2 2016 Dalitz plot, the conversion eta -> pi+pi-gamma (HGeant, phase
// space ASSUMED) with |F_V|^2 p*^2 k^2, and the background rules of the header. This macro redoes the
// eta-peak and eta-enriched comparisons of sim_eta_peak.C, sim_eta_enriched.C and sim_eta_summary.C
// BEFORE (phys = 0, the earlier numbers) and AFTER (phys = simphys::kPhysMain = 108, V3) the weights, and
// for the variants V0, V1, V2 (direct decay), V4 (plane term, check), 1 + aY (bump, the X17 slides) and
// Gounaris-Sakurai F_V. Review of 2026-11-12: main V2 -> V3; the two-template M(pipi) fit uses the variance of
// the rest (no double counting of the bump template); absolute data/sim ratios labelled as scale- and
// method-dependent; scale-free results (class shares, efficiencies) listed first.
//
// Part 1 (bit 1)  the eta peak: composition (direct / conversion) in pre, pre+OA>2, E1 (slides), peak
//                 nominal, E2; peak content per OA and M(ee) class, data (P_main per class, stat (+)
//                 fit-model errors) against simulation (shapes and absolute), one common factor and the
//                 free conversion / direct factors; efficiencies of seven single cuts (data pass/fail fits
//                 against simulation); peak position and width (Gaussian +-2 sigma, 1 MeV bins)
// Part 2 (bit 2)  the bump: simulated bump in the eta-enriched selections before / after the Dalitz weight;
//                 M(pipi) in the bump window, data N - CB against simulation; M(pipi) > 0.34 efficiency of
//                 the photon-tagged bump (data: ECAL-tagged eta of the data study)
// Part 3 (bit 4)  eta-enriched selections (pre for reference, E1, E2, E3, E3oa, E4, E4oa; M4 < 1 GeV/c2):
//                 data N, CB, N - CB against the corrected simulation in absolute scale with ratio pads;
//                 region ratios and the jump at the bump edge; template fits M1 (bump, peak, bkg) and M3
//                 (bkg split eta-related / multi-pion), errors x sqrt(chi2/ndf); the peak data/sim ratio by
//                 the methods A, Aw, B, C, D, E of reconcile.md (pre, pre+OA>2, E1, E2)
// Part 4 (bit 8)  variants V0, V1, V2, V3 (main), V4, bump 1 + aY, F_V Gounaris-Sakurai: headline numbers
//
// Rules: simulation = std (owned, mode A) + lep (owned, mode A) + eta (owned), w (PDG) x physics weight,
// T1 + T2, no CB; variant (a) adds the ghost-flagged true pairs; sigma = sum w / L_sim, L_sim = 221.05
// nb^-1 per set; X17 never summed. Data: (N - CB) / L_data, CB = 2 sqrt(N++ N--) (10 MeV), L_data =
// 6080.5 nb^-1 UNCONFIRMED. pre = isBest >= 0, PT3, vz > -160 mm, cor_ok. All sums in double (TH1D).
// Reuse: #includes sim_eta_summary.C (which includes sim_eta_enriched.C and sim_eta_peak.C) for the
// selections, data peak fits (P_main + 48-variant model spread), template fit and plot helpers; none of
// them is changed and none of their table-writing functions is called. Checks of the 'off' numbers
// against the earlier tables: tables/rw_compare_checks.txt.
// Output: figs/rw_compare_*.png (+ pdf for the key ones), tables/rw_compare_*.txt
// Run from /home/przygoda/HADES/FAT/FAT (about 5-10 min, about 10 GB):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/aclic_rw_compare", kTRUE)' 'plots/eta_sim/sim_rw_compare.C+(15)'
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TPad.h"
#include "TStopwatch.h"
#include "TStyle.h"

#include "sim_eta_summary.C"
#include "sim_eta_physics.h"

namespace src {
using namespace simeng;
typedef std::pair<double, double> VE;
const double kNaN = std::numeric_limits<double>::quiet_NaN();
const double kToRows = kLdata / kLsim;
const int kColData = kBlack, kColCB = kRed + 1, kColSig = kBlue + 1;
const int kColBefore = kGray + 2, kColAfter = kColSim;

// ---------------------------------------------------------------- output
FILE* gT = nullptr;   // current table
FILE* gK = nullptr;   // tables/rw_compare_checks.txt
void T(const char* fmt, ...) {
    char b[32768];
    va_list ap; va_start(ap, fmt); vsnprintf(b, sizeof b, fmt, ap); va_end(ap);
    fputs(b, stdout);
    if (gT) fputs(b, gT);
}
void K(const char* fmt, ...) {
    char b[8192];
    va_list ap; va_start(ap, fmt); vsnprintf(b, sizeof b, fmt, ap); va_end(ap);
    fputs(b, stdout);
    if (gK) fputs(b, gK);
}
void OpenT(const char* n) {
    gT = fopen(Form("%s/%s", kTabDir, n), "w");
    if (!gT) printf("cannot write %s\n", n);
}
void CloseT() { if (gT) fclose(gT); gT = nullptr; }
int gCnt = 0;
std::string Nm(const char* p) { return Form("src_%s_%d", p, gCnt++); }
void Save(TCanvas* c, const char* name, bool pdf) {
    c->SaveAs(Form("%s/%s.png", kFigDir, name));
    if (pdf) c->SaveAs(Form("%s/%s.pdf", kFigDir, name));
}
int gNchk = 0, gNbad = 0;
void Check(const char* what, double mine, double ref, double tol) {
    ++gNchk;
    if (!(ref == ref)) { K("  %-78s this %11.5g  ref: none\n", what, mine); return; }
    const double d = mine - ref;
    const bool ok = std::fabs(d) <= tol;
    if (!ok) ++gNbad;
    K("  %-78s this %11.5g  ref %11.5g  diff %+9.2g  %s\n", what, mine, ref, d, ok ? "OK" : "DIFFERS");
}
VE Ratio(double a, double ea, double b, double eb) {
    const double r = a / b;
    return VE(r, std::fabs(r) * std::hypot(ea / a, eb / b));
}
double Pull(double s, double es, double d, double ed) { return (s - d) / std::sqrt(es * es + ed * ed); }
// ---------------------------------------------------------------- references read from the regenerated study tables
// (replace the rules-v1 numbers that were written into the checks; rules v2 since 2026-11-13)
// first line after the header "## <sec>" (followed by the end of the line, a blank or ':') whose collapsed text starts
// with <start> (followed by the end or a blank); empty if not found
std::string TabLine(const char* file, const std::string& sec, const std::string& start) {
    const std::string h = "## " + sec, st = ses::Collapse(start);
    bool in = false;
    for (const std::string& l : ses::Lines(file)) {
        if (l.compare(0, 3, "## ") == 0) {
            in = l.compare(0, h.size(), h) == 0 && (l.size() == h.size() || l[h.size()] == ' ' || l[h.size()] == ':');
            continue;
        }
        if (!in) continue;
        const std::string c = ses::Collapse(l);
        if (c.compare(0, st.size(), st) == 0 && (c.size() == st.size() || c[st.size()] == ' ')) return l;
    }
    printf("sim_rw_compare: no line '%s' in section '%s' of tables/%s\n", start.c_str(), sec.c_str(), file);
    return "";
}
// tables/rw_weights_xsec.txt (sim_rw_build.C): sigma [nb] 'off' (idx 0) or 'main' (idx 1) of a component in a window
double XsecTab(const char* sec, const char* comp, const char* win, int idx) {
    const std::vector<std::vector<VE>> g = ses::Groups(TabLine("rw_weights_xsec.txt", sec, std::string(comp) + " " + win));
    return g.size() > 1 && idx < (int)g[1].size() ? g[1][idx].first : kNaN;
}
// the same line: ratio with / without weights of a variant (k = 0 V0, 1 V1, 2 V2, 3 V4, 4 1+aY, 5 GS; after the second '|')
double XsecRatioTab(const char* sec, const char* comp, const char* win, int k) {
    const std::string l = TabLine("rw_weights_xsec.txt", sec, std::string(comp) + " " + win);
    const size_t b1 = l.find('|'), b2 = b1 == std::string::npos ? b1 : l.find('|', b1 + 1);
    if (b2 == std::string::npos) return kNaN;
    double x[6];
    const int n = sscanf(l.c_str() + b2 + 1, "%lf %lf %lf %lf %lf %lf", &x[0], &x[1], &x[2], &x[3], &x[4], &x[5]);
    return k < n ? x[k] : kNaN;
}
// tables/rw_weights_peak.txt (sim_rw_build.C): one-factor chi2 of a fit line ("direct V2", "direct V3 main") of a class
// section ("pre, OA5", "pre, MEE5", "E2, OA5"); printed with one decimal (the 'off' fits of that table differ from the ones
// of this macro by up to 0.1 in E2, so the 'before' checks keep the values of tables/eta_summary.txt)
double Chi2Tab(const char* sec, const char* fit) {
    const std::string l = TabLine("rw_weights_peak.txt", sec, fit);
    const size_t p = l.find("chi2/ndf ");
    return p == std::string::npos ? kNaN : atof(l.c_str() + p + 9);
}

// ---------------------------------------------------------------- variants of the physics weights
struct VarDef { std::string tag, lab, tex; int code; };
std::vector<VarDef> gV;
// main = V3 since the review of 2026-11-12 (vV2 = the earlier main variant); vV4 = V3 x plane term (check)
enum { vOff = 0, vV0, vV1, vMain, vV2, vSl, vGS, vV4, kNV };
void InitVars() {
    using namespace simphys;
    gV = {{"off", "no physics weights (as produced)", "no weights", 0},
          {"V0", "direct V0 (KW x p*^3 lambda^3/2) + bump KLOE + conv F_V", "V0", Code(kV0, kBumpKLOE, 1, 1)},
          {"V1", "direct V1 (V0 x |F_V|^2) + bump KLOE + conv F_V", "V1", Code(kV1, kBumpKLOE, 1, 1)},
          {"main", "MAIN: direct V3 (V1 x |F_VMD|^2 x e+ and pi+ decay angles) + bump KLOE + conv F_V", "V3 (main)", kPhysMain},
          {"V2", "direct V2 (V3 without the decay angles; main before the review) + bump KLOE + conv F_V", "V2", kPhysV2},
          {"1+aY", "direct V3 + bump 1 + aY (X17 slides) + conv F_V", "bump 1+aY", Code(kV3, kBumpSlides, 1, 1)},
          {"GS", "main with F_V Gounaris-Sakurai", "F_{V} GS", kPhysMain | (1 << 7)},
          {"V4", "direct V4 (V3 x plane term R_phi, check) + bump KLOE + conv F_V", "V4 (planes)", Code(kV4, kBumpKLOE, 1, 1)}};
}
SimOpt Opt(int iv) { SimOpt o; o.phys = gV[iv].code; return o; }

// ---------------------------------------------------------------- shared state
SimData* gD = nullptr;        // e+e- trees of std, lep, eta; owned rows; M4 < 1.1 (seen::gD points here)
SimData* gS = nullptr;        // reduced copy: template 1, classes 2-5
const scan::Data* gDD = nullptr;

// ---------------------------------------------------------------- components
Comp CDir() { return CompPeakDirect(); }
Comp CConv() { return CompPeakConv(); }
Comp CTot() { return CompPeakTotal(); }
Comp CBump() { return Comp("bump D+C", [](const Row& r) { return r.tmpl == 1 && (r.cls == 2 || r.cls == 3); }); }
Comp CTag() {
    return Comp("bump D+C, photon-tagged", [](const Row& r) {
        return r.tmpl == 1 && (r.cls == 2 || r.cls == 3) && r.g_idx > 0 && r.m_eeg >= 0.119f && r.m_eeg < 0.173f;
    });
}
Sel SWin(const Sel& s, double lo, double hi) { return [s, lo, hi](const Row& r) { return r.m4 >= lo && r.m4 < hi && s(r); }; }
scan::Sel DWin(const scan::Sel& s, double lo, double hi) {
    return [s, lo, hi](const scan::Row& r) { return r.m4 >= lo && r.m4 < hi && s(r); };
}

// ---------------------------------------------------------------- notes for figures
std::string SimNote() {
    return Form("simulation: SMASH std + lep + eta owned rows (mode A), w (PDG) #times physics weight, T1 + T2, #Sigma w / L_{sim}, L_{sim} = %.2f nb^{-1}; "
                "no CB in simulation; no efficiency correction", kLsim);
}
std::string WNote() {
    return "physics weights (main): #eta#rightarrow#pi^{+}#pi^{-}e^{+}e^{-} Kroll-Wada #times p*^{3}#lambda^{3/2}|F_{V}|^{2}|F_{VMD}|^{2} #times decay angles (V3); "
           "#eta#rightarrow#pi^{+}#pi^{-}#pi^{0} KLOE-2 Dalitz plot; #eta#rightarrow#pi^{+}#pi^{-}#gamma |F_{V}|^{2}; HGeant phase space ASSUMED";
}
std::string DNote() {
    return std::string("data: N - CB, CB = 2#sqrt{N_{++}N_{--}} (10 MeV); ") + seen::kLdTex + "; absolute data/sim ratios scale with 6080.5/L_{data}: not physics results";
}

// =====================================================================================================
// Part 1: the eta peak
// =====================================================================================================
const int kNPS = 5;
const char* kPST[kNPS] = {"pre", "pre+OA>2", "E1 (slides)", "peak nominal", "E2"};
const char* kPSs[kNPS] = {"pre", "pre+OA>2", "E1", "peak nom.", "E2"};
Sel PSim(int i) {
    if (i == 0) return SelPreset(kPre);
    if (i == 1) return SelPreset(kPreOA2);
    if (i == 2) return SelPreset(kSlides);
    if (i == 3) return SelPreset(kPeakNom);
    return seen::SimSel(seen::sE2);
}
scan::Sel PDat(int i) {
    if (i == 0) return DataSelPreset(kPre);
    if (i == 1) return DataSelPreset(kPreOA2);
    if (i == 2) return DataSelPreset(kSlides);
    if (i == 3) return DataSelPreset(kPeakNom);
    return seen::DatSel(seen::sE2);
}
const char* kCST[2] = {"OA5", "MEE5"};
const char* kCSx[2] = {"OA(e^{+}e^{-}) class [deg]", "M(e^{+}e^{-}) class [MeV/c^{2}]"};
inline int CIdx(const sep::ClassSet& c, double v) {
    for (size_t j = 0; j < c.lo.size(); ++j) if (v >= c.lo[j] && v < c.hi[j]) return (int)j;
    return -1;
}
struct PeakData { sep::DFit tot; sep::DFit cls[2][5]; scan::FreeFit pos; };
PeakData gPD[kNPS];
struct PeakSim { WSum all[3], win[3]; double c[2][5][2], cv[2][5][2]; };   // [dir, conv, tot]; classes [cs][j][dir, conv]
PeakSim gPS[kNV][kNPS];
bool gHavePS[kNV] = {};

void PeakDataFits() {
    TStopwatch sw;
    std::vector<scan::Sel> sels;
    for (int i = 0; i < kNPS; ++i) {
        const scan::Sel base = PDat(i);
        sels.push_back(base);
        for (int k = 0; k < 2; ++k) {
            const sep::ClassSet cs = sep::FindCS(kCST[k]);
            const int var = cs.var;
            for (size_t j = 0; j < cs.lo.size(); ++j) {
                const double lo = cs.lo[j], hi = cs.hi[j];
                sels.push_back([base, var, lo, hi](const scan::Row& r) {
                    if (!base(r)) return false;
                    const double v = sep::CVar(var, r);
                    return v >= lo && v < hi;
                });
            }
        }
    }
    const std::vector<sep::DFit> f = sep::FitDataSels(*gDD, sels);
    size_t n = 0;
    for (int i = 0; i < kNPS; ++i) {
        gPD[i].tot = f[n++];
        for (int k = 0; k < 2; ++k)
            for (int j = 0; j < 5; ++j) gPD[i].cls[k][j] = f[n++];
        const scan::Spectra sp = scan::Fill(*gDD, PDat(i));
        gPD[i].pos = scan::FitPeakFree(sp);
    }
    printf("PeakDataFits: %zu fits, %.1f s\n", f.size(), sw.RealTime());
}
void PeakSimFill(int iv) {
    const SimOpt o = Opt(iv);
    const std::vector<Comp> cc = {CDir(), CConv(), CTot()};
    for (int i = 0; i < kNPS; ++i) {
        PeakSim& p = gPS[iv][i];
        const Sel s = PSim(i);
        const std::vector<WSum> a = SumMany(*gS, s, cc, o);
        const std::vector<WSum> w = SumMany(*gS, SWin(s, scan::kPeakLo, scan::kPeakHi), cc, o);
        for (int c = 0; c < 3; ++c) { p.all[c] = a[c]; p.win[c] = w[c]; }
        for (int k = 0; k < 2; ++k) {
            const sep::ClassSet cs = sep::FindCS(kCST[k]);
            Binning b; b.nbin = 5; b.lo = 0; b.hi = 5;
            b.x = [cs](const Row& r) { const int j = CIdx(cs, sep::CVar(cs.var, r)); return j < 0 ? -1. : j + 0.5; };
            const std::vector<WSpec> ws = FillMany(*gS, s, {CDir(), CConv()}, o, b);
            for (int j = 0; j < 5; ++j)
                for (int c = 0; c < 2; ++c) { p.c[k][j][c] = ws[c].s[j] * kToRows; p.cv[k][j][c] = ws[c].v[j] * kToRows * kToRows; }
        }
    }
    gHavePS[iv] = true;
}
sep::MixRes Mix(int iv, int ip, int k) {
    std::vector<double> D, eD, C, eC, X, eX;
    for (int j = 0; j < 5; ++j) {
        // a class that is empty by construction (OA < 2 with OA > 2 in the selection: no simulated content) is left out of
        // the fits (error 0 = skipped by sep::MixFit), so that ndf counts only the populated classes
        const bool empty = !(gPS[iv][ip].c[k][j][0] + gPS[iv][ip].c[k][j][1] > 0);
        D.push_back(gPD[ip].cls[k][j].S); eD.push_back(empty ? 0. : gPD[ip].cls[k][j].eT());
        C.push_back(gPS[iv][ip].c[k][j][1]); eC.push_back(std::sqrt(gPS[iv][ip].cv[k][j][1]));
        X.push_back(gPS[iv][ip].c[k][j][0]); eX.push_back(std::sqrt(gPS[iv][ip].cv[k][j][0]));
    }
    return sep::MixFit(D, eD, C, eC, X, eX);
}

// ---- single-cut efficiencies
struct ECut { const char* lab; const char* txt; int kind; double val; bool bstart; const char* keyPeak; const char* keyBump; };
const int kNEC = 7;
const ECut kEC[kNEC] = {
    {"OA > 2#circ", "OA > 2", sep::kOaLo, 2, false, "oa_lo 2 pre", "oa_lo 2 no cut"},
    {"OA < 4#circ", "OA < 4", sep::kOaHi, 4, false, "oa_hi 4 pre", "oa_hi 4 no cut"},
    {"M_{#pi#pi} > 0.34", "M(pipi) > 0.34", sep::kMppLo, 0.34, false, "mpipi_lo 0.34 pre", "mpipi_lo 0.34 no cut"},
    {"#alpha_{#eta} > 160#circ", "eta_rest > 160", sep::kEta, 160, false, "eta_rest 160 pre", "eta_rest 160 no cut"},
    {"p_{T,4} > 400 MeV/c", "pT4 > 400", sep::kPt4, 400, true, "pt4_lo 400 pre", "pt4_lo 400 start"},
    {"#alpha_{lab} < 50#circ", "oa_lab < 50", sep::kLab, 50, false, "oa_lab 50 pre", "oa_lab 50 no cut"},
    {"both #pi in vertex fit", "both pions in vertex fit", sep::kVtx, 2, true, "vtx both pions pre", "vtx both pions start"}};
enum { eDir = 0, eConv, eTot, eBump, eTag, eNE };
struct EffOut { EffRes s[kNV][eNE]; sep::Eff d; VE dTab, dTag, sTagTab; bool have[kNV]; };
EffOut gE[kNEC];
void Efficiencies() {
    TStopwatch sw;
    const sep::Cut base = sep::CPre();
    sep::Cut bstart; bstart.vzLo = -1e9; bstart.labHi = 40; bstart.mppLo = 0.30; bstart.mppHi = 0.41;   // sim_cuts.C bump start
    std::vector<scan::Sel> dsel;
    for (int k = 0; k < kNEC; ++k) {
        sep::Cut cp = base; sep::SetCut(cp, kEC[k].kind, kEC[k].val);
        dsel.push_back(sep::DatSel(cp));
        dsel.push_back([base, cp](const scan::Row& r) { return sep::PassBase(r) && sep::PassCut(r, base) && !sep::PassCut(r, cp); });
    }
    const std::vector<sep::DFit> df = sep::FitDataSels(*gDD, dsel);
    const std::vector<Comp> cc = {CDir(), CConv(), CTot(), CBump(), CTag()};
    for (int k = 0; k < kNEC; ++k) {
        EffOut& o = gE[k];
        o.d = sep::PassFail(df[2 * k], df[2 * k + 1]);
        o.dTab = ses::GetVE("eta_peak_eff.txt", "", kEC[k].keyPeak, 3, 1);
        o.dTag = ses::GetVE("cuts_compare_data.txt", "", kEC[k].keyBump, 2, 0);
        o.sTagTab = ses::GetVE("cuts_compare_data.txt", "", kEC[k].keyBump, 2, 1);
        sep::Cut cp = base; sep::SetCut(cp, kEC[k].kind, kEC[k].val);
        sep::Cut cb = bstart; sep::SetCut(cb, kEC[k].kind, kEC[k].val);
        for (int iv = 0; iv < kNV; ++iv) {
            const SimOpt opt = Opt(iv);
            for (int c = 0; c < eNE; ++c) {
                const bool bs = kEC[k].bstart && (c == eBump || c == eTag);
                o.s[iv][c] = bs ? Efficiency(*gS, sep::SimSel(bstart), sep::SimSel(cb), cc[c], opt)
                                : Efficiency(*gS, sep::SimSel(base), sep::SimSel(cp), cc[c], opt);
            }
            o.have[iv] = true;
        }
    }
    printf("Efficiencies: %.1f s\n", sw.RealTime());
}

// ---- position and width
struct PosOut { sep::GFit g[3]; double tm[3], tmE[3], trms[3]; };
PosOut gPos[2][kNPS];
TH1D* gPosH[2][3];
void PeakPositions() {
    const std::vector<Comp> cc = {CDir(), CConv(), CTot()};
    for (int m = 0; m < 2; ++m) {
        const int iv = m ? vMain : vOff;
        const SimOpt o = Opt(iv);
        for (int ip = 0; ip < kNPS; ++ip) {
            Binning b; b.lo = 0.49; b.hi = 0.60; b.nbin = 110;
            const std::vector<WSpec> w = FillMany(*gS, PSim(ip), cc, o, b);
            for (int c = 0; c < 3; ++c) {
                TH1D* h = w[c].Hist(Nm("pk").c_str(), 0);
                gPos[m][ip].g[c] = sep::GaussIter(h, 2.0);
                if (ip == 0) gPosH[m][c] = h;
            }
            // truncated mean and rms within +-12 MeV of the weighted median
            std::vector<std::pair<double, double>> v[3];
            const Sel s = PSim(ip);
            for (int set = 1; set <= 3; ++set)
                for (const Row& r : gS->tree(set, kUS)) {
                    if (!s(r)) continue;
                    const double w0 = Weight(r, o);
                    if (!(w0 != 0) || !(r.m4 >= 0.49f && r.m4 < 0.60f)) continue;
                    const int c = r.cls == 5 ? 0 : (r.cls == 4 ? 1 : -1);
                    if (c < 0) continue;
                    v[c].push_back({r.m4, w0}); v[2].push_back({r.m4, w0});
                }
            for (int c = 0; c < 3; ++c) {
                const double med = sep::WQuant(v[c], 0.5);
                double s0 = 0, s1 = 0, s2 = 0, sww = 0;
                for (const auto& p : v[c])
                    if (std::fabs(p.first - med) < 0.012) { s0 += p.second; s1 += p.second * p.first; s2 += p.second * p.first * p.first; sww += p.second * p.second; }
                const double mean = s1 / s0, rms = std::sqrt(std::max(0., s2 / s0 - mean * mean)), neff = s0 * s0 / sww;
                gPos[m][ip].tm[c] = mean; gPos[m][ip].tmE[c] = rms / std::sqrt(neff); gPos[m][ip].trms[c] = rms;
            }
        }
    }
}

void PeakTables() {
    OpenT("rw_compare_peak.txt");
    T("# rw_compare_peak: the eta peak BEFORE (no physics weights) and AFTER (main weights: direct V3 = V2 x decay angles, bump KLOE, conversion F_V)\n");
    T("#   - sim_rw_compare.C; V2 (no decay angles, the main variant before the review of 2026-11-12) is listed where it matters\n");
    T("# READ FIRST: scale-free results (class shares, section 0; efficiencies, section 4) are the results; absolute contents and data/sim\n");
    // rules v2 (2026-11-13): methods A-E after the weights 0.279-1.025, ratio 3.67 (rules v1: 0.284-1.032, 3.63)
    T("#   ratios (sections 1-3) scale with 6080.5 / L_data (UNCONFIRMED) and the peak ratio is method-dependent (A-E differ by a factor of about 3.7)\n");
    T("# simulation: std + lep + eta owned rows (mode A), w (PDG) x physics weight; direct = eta -> pi+pi-e+e- (class 5, eta set), conversion = eta -> pi+pi-gamma\n");
    T("#   (class 4, std); expected data rows = sum w x L_data / L_sim, L_data = 6080.5 nb^-1 UNCONFIRMED (every data/sim scales with 6080.5 / L_data)\n");
    T("# data: P_main of N - CB (Gaussian frozen 0.5439 / 5.4 MeV, pol2, 0.5225-0.61) per selection and per class; errors stat (+) fit-model spread\n");
    T("#   (48 variants of the data study); selections: pre = isBest >= 0, PT3, vz > -160, cor_ok; E1 = slides; E2 = E1 without OA > 2\n\n");
    T("## 0. scale-free: share of the peak at high M(ee) and high OA (data P_main per class, stat (+) model, classes uncorrelated)\n");
    T("# the one-factor class fits have little power (data errors 30-40 %% per class; before the weights p = 0.05 / 0.03 in pre already);\n");
    T("# these shares are the clearest test: the generator M(ee) of the direct decay put almost nothing at M(ee) > 20 MeV\n");
    T("%-13s %-12s | %-15s | %-8s %-8s %-8s | %-22s\n", "selection", "share", "data", "sim bef", "sim V2", "sim V3", "pull bef / V2 / V3");
    for (int ip = 0; ip < kNPS; ++ip)
        for (int k = 0; k < 2; ++k) {
            const int j0 = k == 1 ? 2 : 3;   // M(ee) > 20 MeV: classes 2-4; OA > 6 deg: classes 3-4
            double dH = 0, dL = 0, vH = 0, vL = 0;
            for (int j = 0; j < 5; ++j) {
                const double x = gPD[ip].cls[k][j].S, e = gPD[ip].cls[k][j].eT();
                if (j >= j0) { dH += x; vH += e * e; } else { dL += x; vL += e * e; }
            }
            const double fD = dH / (dH + dL), eD = std::sqrt((1 - fD) * (1 - fD) * vH + fD * fD * vL) / (dH + dL);
            double fs[3];
            const int ivs[3] = {vOff, vV2, vMain};
            for (int q = 0; q < 3; ++q) {
                double h = 0, t = 0;
                for (int j = 0; j < 5; ++j) { const double x = gPS[ivs[q]][ip].c[k][j][0] + gPS[ivs[q]][ip].c[k][j][1]; t += x; if (j >= j0) h += x; }
                fs[q] = h / t;
            }
            T("%-13s %-12s | %5.3f +- %-6.3f | %8.3f %8.3f %8.3f | %+5.1f / %+5.1f / %+5.1f\n", kPST[ip], k == 1 ? "M(ee) > 20" : "OA > 6 deg", fD, eD, fs[0], fs[1], fs[2],
              (fs[0] - fD) / eD, (fs[1] - fD) / eD, (fs[2] - fD) / eD);
        }
    T("# pulls with the data error only (simulation statistical errors are below 0.01)\n\n");
    T("## 1. composition [nb] (all M4 | peak window 0.5336-0.5542), before -> after\n");
    T("%-13s | %-34s | %-34s | %-30s | %-17s\n", "selection", "direct all M4 (window)", "conversion all M4 (window)", "total all M4 (window)", "direct share");
    for (int ip = 0; ip < kNPS; ++ip) {
        const PeakSim& a = gPS[vOff][ip]; const PeakSim& b = gPS[vMain][ip];
        T("%-13s | %6.4f -> %6.4f (%6.4f -> %6.4f) | %6.4f -> %6.4f (%6.4f -> %6.4f) | %6.3f -> %6.3f (%5.3f -> %5.3f) | %5.3f -> %5.3f\n", kPST[ip],
          a.all[0].Xsec(), b.all[0].Xsec(), a.win[0].Xsec(), b.win[0].Xsec(), a.all[1].Xsec(), b.all[1].Xsec(), a.win[1].Xsec(), b.win[1].Xsec(),
          a.all[2].Xsec(), b.all[2].Xsec(), a.win[2].Xsec(), b.win[2].Xsec(), a.all[0].S / a.all[2].S, b.all[0].S / b.all[2].S);
    }
    T("# errors (event groups, after): ");
    for (int ip = 0; ip < kNPS; ++ip) T("%s direct %.4f conv %.4f total %.4f nb; ", kPSs[ip], gPS[vMain][ip].all[0].XsecErr(), gPS[vMain][ip].all[1].XsecErr(), gPS[vMain][ip].all[2].XsecErr());
    T("\n\n## 2. peak content: data P_main against the simulated peak over all M4 (method A), expected data rows\n");
    T("# NOT a physics result: scales with 6080.5 / L_data (UNCONFIRMED) and depends on the method (A-E differ by a factor of about 3.7 after the weights,\n");
    T("#   rw_compare_enriched.txt section 4)\n");
    T("%-13s %24s %22s %22s %18s %18s\n", "selection", "data P_main (stat, model)", "sim total before", "sim total after", "data/sim before", "data/sim after");
    for (int ip = 0; ip < kNPS; ++ip) {
        const sep::DFit& d = gPD[ip].tot;
        const WSum& a = gPS[vOff][ip].all[2]; const WSum& b = gPS[vMain][ip].all[2];
        const VE ra = Ratio(d.S, d.eT(), a.S * kToRows, std::sqrt(a.V) * kToRows), rb = Ratio(d.S, d.eT(), b.S * kToRows, std::sqrt(b.V) * kToRows);
        T("%-13s %7.0f +- %4.0f +- %-5.0f %10.0f +- %-6.0f %10.0f +- %-6.0f %8.3f +- %-6.3f %8.3f +- %-6.3f\n", kPST[ip], d.S, d.eS, d.eM, a.S * kToRows,
          std::sqrt(a.V) * kToRows, b.S * kToRows, std::sqrt(b.V) * kToRows, ra.first, ra.second, rb.first, rb.second);
    }
    T("# errors of data/sim: stat (+) model of the data and simulation statistics; NOT the absolute scale (about 20 %%, L_data unconfirmed)\n");
    T("\n## 3. peak content per class (rows at L_data): data against direct / conversion / total, before -> after; fits\n");
    T("# one factor: data_i = a1 x (sim total_i); two factors: data_i = a x conv_i + b x direct_i (sep::MixFit; sim errors in the variance)\n");
    for (int ip = 0; ip < kNPS; ++ip)
        for (int k = 0; k < 2; ++k) {
            const sep::ClassSet cs = sep::FindCS(kCST[k]);
            T("\n### %s, %s classes (%s)\n", kPST[ip], kCST[k], k == 0 ? "OA [deg]" : "M(ee) [GeV/c2]");
            T("%-8s | %18s %13s | %8s %8s %8s | %8s %8s %8s | %-15s %-15s | %-13s %6s %6s\n", "class", "data (stat+model)", "frac data", "dir bef", "conv bef",
              "tot bef", "dir aft", "conv aft", "tot aft", "data/sim bef", "data/sim aft", "frac: data", "bef", "aft");
            std::vector<double> dv, de;
            double tb = 0, ta = 0;
            for (int j = 0; j < 5; ++j) {
                dv.push_back(gPD[ip].cls[k][j].S); de.push_back(gPD[ip].cls[k][j].eT());
                tb += gPS[vOff][ip].c[k][j][0] + gPS[vOff][ip].c[k][j][1]; ta += gPS[vMain][ip].c[k][j][0] + gPS[vMain][ip].c[k][j][1];
            }
            for (int j = 0; j < 5; ++j) {
                const double* cb = gPS[vOff][ip].c[k][j]; const double* ca = gPS[vMain][ip].c[k][j];
                const VE fr = ses::Frac(dv, de, j);
                const double sb = cb[0] + cb[1], sa = ca[0] + ca[1];
                const double eb = std::sqrt(gPS[vOff][ip].cv[k][j][0] + gPS[vOff][ip].cv[k][j][1]), ea = std::sqrt(gPS[vMain][ip].cv[k][j][0] + gPS[vMain][ip].cv[k][j][1]);
                const VE rb = Ratio(dv[j], de[j], sb, eb), ra = Ratio(dv[j], de[j], sa, ea);
                T("%-8s | %7.0f +- %-7.0f %5.3f+-%5.3f | %8.0f %8.0f %8.0f | %8.0f %8.0f %8.0f | %6.3f +- %-6.3f %6.3f +- %-6.3f | %5.3f+-%5.3f %6.3f %6.3f\n",
                  cs.name[j].c_str(), dv[j], de[j], fr.first, fr.second, cb[0], cb[1], sb, ca[0], ca[1], sa, rb.first, rb.second, ra.first, ra.second,
                  fr.first, fr.second, sb / tb, sa / ta);
            }
            for (int m = 0; m < 2; ++m) {
                const sep::MixRes x = Mix(m ? vMain : vOff, ip, k);
                T("  %-6s one factor a1 = %.3f +- %.3f, chi2/ndf %.1f/%d (p %.3f) | two factors: conv %+.3f +- %.3f, direct %+.3f +- %.3f (rho %+.2f), chi2/ndf %.1f/%d;"
                  " direct share of the sim peak %.3f, implied in data %.2f +- %.2f\n",
                  m ? "after" : "before", x.a1, x.ea1, x.chi21, x.ndf1, TMath::Prob(x.chi21, x.ndf1), x.a, x.ea, x.b, x.eb, x.rho, x.chi2, x.ndf, x.fdirSim,
                  x.fdir, x.efdir);
            }
        }
    T("\n## 4. single-cut efficiencies, base pre, all M4 (simulation: event-group errors; data peak: pass/fail P_main fits, stat (+) model)\n");
    T("# bump: D+C all M4 and photon-tagged (0.119 <= M(ee gamma) < 0.173, emulates the ECAL tag, UNVERIFIED); pT4 and vertex: bump at the bump start base\n");
    T("#   (pre without vz + oa_lab < 40 + 0.30 < M(pipi) < 0.41), as the data study; data bump = ECAL-tagged eta of the data study (tables/cuts_compare_data.txt)\n");
    T("# total V2 = the peak total with the direct decay in V2 (no decay angles; the main variant before the review)\n");
    T("%-26s | %-29s %-29s %-29s %-13s | %-15s | %-8s %-8s %-8s | %-34s | %-15s\n", "cut", "direct bef -> aft", "conversion bef -> aft", "total bef -> aft", "total V2",
      "data peak", "pull bef", "pull V2", "pull aft", "bump tagged bef / KLOE / 1+aY", "data tagged");
    for (int k = 0; k < kNEC; ++k) {
        const EffOut& o = gE[k];
        const EffRes* b = o.s[vOff]; const EffRes* a = o.s[vMain];
        const EffRes* c2 = o.s[vV2];
        T("%-26s | %5.3f+-%5.3f -> %5.3f+-%5.3f %5.3f+-%5.3f -> %5.3f+-%5.3f %5.3f+-%5.3f -> %5.3f+-%5.3f %5.3f+-%5.3f | %5.3f +- %5.3f | %+6.1f   %+6.1f   %+6.1f   | "
          "%5.3f+-%5.3f / %5.3f / %5.3f       | %5.3f +- %5.3f\n",
          kEC[k].txt, b[eDir].eps, b[eDir].err, a[eDir].eps, a[eDir].err, b[eConv].eps, b[eConv].err, a[eConv].eps, a[eConv].err, b[eTot].eps, b[eTot].err,
          a[eTot].eps, a[eTot].err, c2[eTot].eps, c2[eTot].err, o.d.e, o.d.err, Pull(b[eTot].eps, b[eTot].err, o.d.e, o.d.err),
          Pull(c2[eTot].eps, c2[eTot].err, o.d.e, o.d.err), Pull(a[eTot].eps, a[eTot].err, o.d.e, o.d.err),
          b[eTag].eps, b[eTag].err, a[eTag].eps, o.s[vSl][eTag].eps, o.dTag.first, o.dTag.second);
    }
    T("# pulls = (sim total - data) / sqrt(err_sim^2 + err_data^2); data peak of M(pipi) cuts and MM4 are sculpted (data study); eta_rest > 160 data:\n");
    T("#   fail-sample fits have p < 0.05 in all variants (review)\n");
    T("# bump all D+C: ");
    for (int k = 0; k < kNEC; ++k) T("%s %.3f -> %.3f (1+aY %.3f); ", kEC[k].txt, gE[k].s[vOff][eBump].eps, gE[k].s[vMain][eBump].eps, gE[k].s[vSl][eBump].eps);
    T("\n\n## 5. position and width of the simulated peak: Gaussian +-2 sigma iterated, 1 MeV bins, reconstructed corrected M4; data FitPeakFree (same selection)\n");
    T("%-13s %-11s | %-26s %-22s | %-26s %-22s | %-20s %-12s\n", "selection", "sample", "mean before [GeV/c2]", "sigma before [MeV]", "mean after [GeV/c2]",
      "sigma after [MeV]", "trunc. mean bef/aft", "trunc rms bef/aft [MeV]");
    const char* cn[3] = {"direct", "conversion", "total"};
    for (int ip = 0; ip < kNPS; ++ip) {
        for (int c = 0; c < 3; ++c) {
            const sep::GFit& b = gPos[0][ip].g[c]; const sep::GFit& a = gPos[1][ip].g[c];
            T("%-13s %-11s | %10.5f +- %-10.5f %7.2f +- %-8.2f | %10.5f +- %-10.5f %7.2f +- %-8.2f | %8.5f / %8.5f  %5.2f / %5.2f\n", kPST[ip], cn[c], b.m, b.em,
              b.s * 1e3, b.es * 1e3, a.m, a.em, a.s * 1e3, a.es * 1e3, gPos[0][ip].tm[c], gPos[1][ip].tm[c], gPos[0][ip].trms[c] * 1e3, gPos[1][ip].trms[c] * 1e3);
        }
        const scan::FreeFit& f = gPD[ip].pos;
        T("%-13s %-11s | %10.5f +- %-10.5f %7.2f +- %-8.2f | (data, statistical errors; fit-choice spread in tables/eta_peak_position.txt)\n", kPST[ip], "data",
          f.mean, f.eMean, f.sigma * 1e3, f.eSigma * 1e3);
    }
    T("# K0S scale sim/data 1.0048 +- 0.0007: a simulated position of 0.5459 corresponds to 0.5433 in data units\n");
    CloseT();
    // checks of the 'off' numbers
    K("\n## peak: 'off' against the earlier tables\n");
    Check("pre direct all M4 [nb] (eta_engine_reference 1.385)", gPS[vOff][0].all[0].Xsec(), 1.3852, 6e-4);
    // references of the weight-rules-dependent checks: read from the regenerated tables (rules v2, 2026-11-13); the rules-v1
    // values were 1.1792, 1.2537, 1.8471, 1.3832, 2.0038, -1.11, chi2 1.6 / 3.9 / 1.2 / 0.7 / 3.3 and OA > 2 0.591 / 0.700 / 0.711
    Check("pre conversion all M4 [nb] (rw_weights_xsec)", gPS[vOff][0].all[1].Xsec(), XsecTab("pre", "peak conversion (4)", "all M4", 0), 6e-4);
    // main = V3 (rw_weights_xsec.txt of the corrected sim_rw_build.C); V2 = the earlier main (earlier rw_weights_xsec / rw_compare values)
    Check("pre direct all M4 after V3 [nb] (rw_weights_xsec 1.4168)", gPS[vMain][0].all[0].Xsec(), 1.4168, 6e-4);
    Check("E1 peak total window after V3 [nb] (rw_weights_xsec)", gPS[vMain][2].win[2].Xsec(), XsecTab("E1", "peak total (4+5)", "peak win", 1), 6e-4);
    Check("E2 peak total window after V3 [nb] (rw_weights_xsec)", gPS[vMain][4].win[2].Xsec(), XsecTab("E2", "peak total (4+5)", "peak win", 1), 6e-4);
    Check("pre direct all M4 V2 [nb] (earlier main 1.6442)", gPS[vV2][0].all[0].Xsec(), 1.6442, 6e-4);
    Check("E1 peak total window V2 [nb] (rw_weights_xsec off x ratio V2)", gPS[vV2][2].win[2].Xsec(),
          XsecTab("E1", "peak total (4+5)", "peak win", 0) * XsecRatioTab("E1", "peak total (4+5)", "peak win", 2), 6e-4);
    Check("E2 peak total window V2 [nb] (rw_weights_xsec off x ratio V2)", gPS[vV2][4].win[2].Xsec(),
          XsecTab("E2", "peak total (4+5)", "peak win", 0) * XsecRatioTab("E2", "peak total (4+5)", "peak win", 2), 6e-4);
    Check("pre data P_main [rows] (eta_summary_peak_methods 4671)", gPD[0].tot.S, 4671, 0.6);
    Check("pre data model spread [rows] (958)", gPD[0].tot.eM, 958, 0.6);
    Check("E2 data P_main [rows] (4047)", gPD[4].tot.S, 4047, 0.6);
    Check("E1 data P_main [rows] (3187)", gPD[2].tot.S, 3187, 0.6);
    Check("pre OA <2 data [rows] (rw_weights_peak 1039)", gPD[0].cls[0][0].S, 1039, 0.6);
    Check("pre MEE 20-50 data [rows] (751)", gPD[0].cls[1][2].S, 751, 0.6);
    Check("E2 OA >10 data [rows] (865)", gPD[4].cls[0][4].S, 865, 0.6);
    Check("pre OA one factor chi2 before (9.5)", Mix(vOff, 0, 0).chi21, 9.5, 0.051);
    Check("pre MEE one factor chi2 before (10.8)", Mix(vOff, 0, 1).chi21, 10.8, 0.051);
    Check("E2 OA one factor chi2 before (14.8)", Mix(vOff, 4, 0).chi21, 14.8, 0.051);
    Check("pre OA two factors conv before (eta_summary)", Mix(vOff, 0, 0).a,
          ses::GetVE("eta_summary.txt", "4.", "free conversion / direct factors:", -1, 0).first, 0.006);
    Check("pre OA one factor chi2 after V3 (rw_weights_peak)", Mix(vMain, 0, 0).chi21, Chi2Tab("pre, OA5", "direct V3 main"), 0.051);
    Check("pre MEE one factor chi2 after V3 (1.0)", Mix(vMain, 0, 1).chi21, 1.0, 0.051);
    Check("E2 OA one factor chi2 after V3 (rw_weights_peak)", Mix(vMain, 4, 0).chi21, Chi2Tab("E2, OA5", "direct V3 main"), 0.051);
    Check("pre OA one factor chi2 V2 (rw_weights_peak)", Mix(vV2, 0, 0).chi21, Chi2Tab("pre, OA5", "direct V2"), 0.051);
    Check("pre MEE one factor chi2 V2 (rw_weights_peak)", Mix(vV2, 0, 1).chi21, Chi2Tab("pre, MEE5", "direct V2"), 0.051);
    Check("E2 OA one factor chi2 V2 (rw_weights_peak)", Mix(vV2, 4, 0).chi21, Chi2Tab("E2, OA5", "direct V2"), 0.051);
    for (int k = 0; k < kNEC; ++k) {
        Check(Form("%s: data peak eff (eta_peak_eff)", kEC[k].txt), gE[k].d.e, gE[k].dTab.first, 6e-4);
        Check(Form("%s: sim bump tagged eff before (cuts_compare_data)", kEC[k].txt), gE[k].s[vOff][eTag].eps, gE[k].sTagTab.first, 6e-4);
    }
    Check("OA > 2: sim total before (rw_weights_peak)", gE[0].s[vOff][eTot].eps, ses::GetVE("rw_weights_peak.txt", "single-cut", "OA > 2", 3, 0).first, 6e-4);
    Check("OA > 2: sim total after V3 (rw_weights_peak)", gE[0].s[vMain][eTot].eps, ses::GetVE("rw_weights_peak.txt", "single-cut", "OA > 2", 3, 1).first, 6e-4);
    // no table holds the V2 total: rules v2 (2026-11-13) value of this macro (rules v1: 0.711)
    Check("OA > 2: sim total V2 (0.708)", gE[0].s[vV2][eTot].eps, 0.708, 6e-4);
    Check("M(pipi) > 0.34: bump all after KLOE (rw_weights_peak 0.681)", gE[2].s[vMain][eBump].eps, 0.681, 6e-4);
    Check("pre sim total Gaussian mean before (eta_peak_position 0.54605 +- 0.00024 fit choice)", gPos[0][0].g[2].m, 0.54605, 3e-4);
}

// ---- figures of part 1
void FigPeakClasses() {
    TCanvas* cv = new TCanvas(Nm("cPC").c_str(), "", 3200, 1900);
    cv->cd();
    seen::MkPad(Nm("h").c_str(), 0, 0.90, 1, 1, false, 0, 0, 0, 0);
    seen::Txt(0.01, 0.72, "#eta peak per class of OA(e^{+}e^{-}) and M(e^{+}e^{-}): data against the simulation BEFORE (grey, open) and AFTER (teal, full) the physics weights", 0.21);
    seen::Txt(0.01, 0.44, (WNote() + "; direct = #eta set, class 5; conversion = std, class 4").c_str(), 0.13);
    seen::Txt(0.01, 0.16, (DNote() + "; data: P_{main} per class, stat #oplus fit-model spread (48 variants); pre = isBest #geq 0, PT3, v_{z} > -160 mm; E2 = slides without OA > 2#circ").c_str(), 0.13);
    const int ips[2] = {0, 4};
    for (int row = 0; row < 2; ++row) {
        const int ip = ips[row];
        const double y1 = row == 0 ? 0.46 : 0.0, y2 = row == 0 ? 0.90 : 0.44;
        for (int k = 0; k < 2; ++k) {
            cv->cd();
            seen::MkPad(Nm("p").c_str(), k * 0.25, y1, (k + 1) * 0.25, y2, false, 0.15, 0.02, 0.09, 0.13);
            const sep::ClassSet cs = sep::FindCS(kCST[k]);
            TH1D* fr = ses::CatFrame(Nm("f").c_str(), cs.name, "fraction of the peak content", 0, 1.05, 0.05, 0.06, 1.35);
            fr->GetXaxis()->SetTitle(kCSx[k]); fr->GetXaxis()->SetTitleSize(0.05); fr->GetXaxis()->SetTitleOffset(1.1);
            fr->Draw("axis");
            std::vector<double> dv, de, x0, x1, x2;
            std::vector<VE> fd, fb, fa, fdir, fconv;
            double tb = 0, ta = 0, tda = 0, tca = 0;
            for (int j = 0; j < 5; ++j) {
                dv.push_back(gPD[ip].cls[k][j].S); de.push_back(gPD[ip].cls[k][j].eT());
                tb += gPS[vOff][ip].c[k][j][0] + gPS[vOff][ip].c[k][j][1]; ta += gPS[vMain][ip].c[k][j][0] + gPS[vMain][ip].c[k][j][1];
                tda += gPS[vMain][ip].c[k][j][0]; tca += gPS[vMain][ip].c[k][j][1];
            }
            for (int j = 0; j < 5; ++j) {
                fd.push_back(ses::Frac(dv, de, j));
                fb.push_back(VE((gPS[vOff][ip].c[k][j][0] + gPS[vOff][ip].c[k][j][1]) / tb, 0.));
                fa.push_back(VE((gPS[vMain][ip].c[k][j][0] + gPS[vMain][ip].c[k][j][1]) / ta, 0.));
                fdir.push_back(VE(gPS[vMain][ip].c[k][j][0] / tda, 0.)); fconv.push_back(VE(gPS[vMain][ip].c[k][j][1] / tca, 0.));
                x0.push_back(j - 0.18); x1.push_back(j); x2.push_back(j + 0.18);
            }
            TGraphErrors* gdir = ses::Pts(x1, fdir, kColPeakDirect, 0, 0, 3); gdir->SetLineStyle(2); gdir->Draw("L same");
            TGraphErrors* gcon = ses::Pts(x1, fconv, kColPeakConv, 0, 0, 3); gcon->SetLineStyle(3); gcon->Draw("L same");
            TGraphErrors* gb = ses::Pts(x0, fb, kColBefore, 25, 2.0); gb->Draw("P same");
            TGraphErrors* ga = ses::Pts(x2, fa, kColAfter, 21, 2.0); ga->Draw("P same");
            TGraphErrors* gd = ses::Pts(x1, fd, kColData, 20, 2.0, 3); gd->Draw("P same");
            const sep::MixRes mb = Mix(vOff, ip, k), ma = Mix(vMain, ip, k);
            seen::Txt(0.18, 0.93, Form("#bf{%s}, %s", kPST[ip], kCST[k][0] == 'O' ? "OA classes" : "M_{ee} classes"), 0.055);
            seen::Txt(0.17, 0.84, Form("one factor #chi^{2}/ndf: %.1f/%d #rightarrow %.1f/%d", mb.chi21, mb.ndf1, ma.chi21, ma.ndf1), 0.040);
            seen::Txt(0.17, 0.79, Form("2 factors after: a_{conv} %.2f#pm%.2f, a_{dir} %.2f#pm%.2f", ma.a, ma.ea, ma.b, ma.eb), 0.036);
            seen::Txt(0.17, 0.74, Form("before: a_{conv} %.2f#pm%.2f, a_{dir} %.2f#pm%.2f", mb.a, mb.ea, mb.b, mb.eb), 0.036, kColBefore);
            {   // scale-free share at high M(ee) / high OA (the clearest test; the class chi2 tests have little power)
                const int j0 = k == 1 ? 2 : 3;
                double dH = 0, dL = 0, vH = 0, vL = 0, sH[2] = {0, 0}, sT[2] = {0, 0};
                for (int j = 0; j < 5; ++j) {
                    const double x = gPD[ip].cls[k][j].S, e = gPD[ip].cls[k][j].eT();
                    if (j >= j0) { dH += x; vH += e * e; } else { dL += x; vL += e * e; }
                    for (int m = 0; m < 2; ++m) {
                        const int iv = m ? vMain : vOff;
                        const double y = gPS[iv][ip].c[k][j][0] + gPS[iv][ip].c[k][j][1];
                        sT[m] += y; if (j >= j0) sH[m] += y;
                    }
                }
                const double fD = dH / (dH + dL), eD = std::sqrt((1 - fD) * (1 - fD) * vH + fD * fD * vL) / (dH + dL);
                seen::Txt(0.17, 0.69, Form("%s: data %.2f#pm%.2f, sim %.3f#rightarrow%.2f", k == 1 ? "M_{ee} > 20 MeV" : "OA > 6#circ", fD, eD, sH[0] / sT[0], sH[1] / sT[1]), 0.033);
            }
            if (row == 0 && k == 0) {
                TLegend* lg = new TLegend(0.45, 0.42, 0.97, 0.64); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.034);
                lg->AddEntry(gd, "data (P_{main}, stat #oplus model)", "pe");
                lg->AddEntry(gb, "sim peak total, before", "p");
                lg->AddEntry(ga, "sim peak total, after", "p");
                lg->AddEntry(gdir, "sim direct, after", "l");
                lg->AddEntry(gcon, "sim conversion, after", "l");
                lg->Draw();
            }
        }
    }
    // column 3: composition (rows at L_data) per selection; column 4: fit factors
    cv->cd();
    seen::MkPad(Nm("p").c_str(), 0.50, 0.46, 0.75, 0.90, true, 0.15, 0.02, 0.09, 0.13);
    {
        std::vector<std::string> lab(kPSs, kPSs + kNPS);
        TH1D* fr = ses::CatFrame(Nm("f").c_str(), lab, "peak content [rows at L_{data}]", 1500, 300000, 0.05, 0.05, 1.35);
        fr->Draw("axis");
        std::vector<double> x(kNPS), xb(kNPS), xa(kNPS), xd(kNPS), xc(kNPS);
        std::vector<VE> d, b, a, dd, cc, db, cb;
        for (int ip = 0; ip < kNPS; ++ip) {
            x[ip] = ip; xb[ip] = ip - 0.2; xa[ip] = ip + 0.2; xd[ip] = ip + 0.1; xc[ip] = ip - 0.1;
            d.push_back(VE(gPD[ip].tot.S, gPD[ip].tot.eT()));
            b.push_back(VE(gPS[vOff][ip].all[2].S * kToRows, std::sqrt(gPS[vOff][ip].all[2].V) * kToRows));
            a.push_back(VE(gPS[vMain][ip].all[2].S * kToRows, std::sqrt(gPS[vMain][ip].all[2].V) * kToRows));
            dd.push_back(VE(gPS[vMain][ip].all[0].S * kToRows, 0)); cc.push_back(VE(gPS[vMain][ip].all[1].S * kToRows, 0));
            db.push_back(VE(gPS[vOff][ip].all[0].S * kToRows, 0)); cb.push_back(VE(gPS[vOff][ip].all[1].S * kToRows, 0));
        }
        ses::Pts(xb, b, kColBefore, 25, 2.0)->Draw("P same");
        ses::Pts(xa, a, kColAfter, 21, 2.0)->Draw("P same");
        ses::Pts(xd, dd, kColPeakDirect, 22, 1.8)->Draw("P same");
        ses::Pts(xc, cc, kColPeakConv, 23, 1.8)->Draw("P same");
        ses::Pts(xd, db, kColPeakDirect, 26, 1.6)->Draw("P same");
        ses::Pts(xc, cb, kColPeakConv, 32, 1.6)->Draw("P same");
        ses::Pts(x, d, kColData, 20, 2.0, 3)->Draw("P same");
        seen::Txt(0.18, 0.93, "#bf{peak content (all M_{4})}", 0.055);
        TLegend* lg = new TLegend(0.17, 0.66, 0.97, 0.88); lg->SetNColumns(2); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.036);
        lg->AddEntry(seen::Dummy(Nm("d").c_str(), kColData, 2, 1, 0, 20, kColData), "data P_{main}", "p");
        lg->AddEntry(seen::Dummy(Nm("d").c_str(), kColAfter, 2, 1, 0, 21, kColAfter), "sim total after", "p");
        lg->AddEntry(seen::Dummy(Nm("d").c_str(), kColBefore, 2, 1, 0, 25, kColBefore), "sim total before", "p");
        lg->AddEntry(seen::Dummy(Nm("d").c_str(), kColPeakDirect, 2, 1, 0, 22, kColPeakDirect), "direct after (open: before)", "p");
        lg->AddEntry(seen::Dummy(Nm("d").c_str(), kColPeakConv, 2, 1, 0, 23, kColPeakConv), "conversion after (open: before)", "p");
        lg->Draw();
    }
    cv->cd();
    seen::MkPad(Nm("p").c_str(), 0.50, 0.0, 0.75, 0.44, false, 0.15, 0.02, 0.09, 0.13);
    {
        std::vector<std::string> lab(kPSs, kPSs + kNPS);
        TH1D* fr = ses::CatFrame(Nm("f").c_str(), lab, "one-factor fit: #chi^{2}/ndf", 0, 5.0, 0.05, 0.05, 1.35);
        fr->Draw("axis");
        ses::HLine(-0.5, kNPS - 0.5, 1., kBlack, 2, 1);
        const int mk[2][2] = {{24, 25}, {20, 21}};
        for (int m = 0; m < 2; ++m)
            for (int k = 0; k < 2; ++k) {
                std::vector<double> x; std::vector<VE> y;
                for (int ip = 0; ip < kNPS; ++ip) {
                    const sep::MixRes r = Mix(m ? vMain : vOff, ip, k);
                    x.push_back(ip + (k ? 0.12 : -0.12)); y.push_back(VE(r.chi21 / r.ndf1, 0));
                }
                ses::Pts(x, y, m ? kColAfter : kColBefore, mk[m][k], 2.0)->Draw("P same");
            }
        seen::Txt(0.18, 0.93, "#bf{one common factor on the classes}", 0.055);
        seen::Txt(0.18, 0.85, "circles: OA, squares: M_{ee} classes; open grey: before, teal: after", 0.034);
        seen::Txt(0.18, 0.80, "low power: data errors 30-40 % per class", 0.034);
        seen::Txt(0.18, 0.75, "(before: p = 0.05 / 0.03 in pre)", 0.034);
    }
    cv->cd();
    seen::MkPad(Nm("p").c_str(), 0.75, 0.0, 1.0, 0.90, false, 0.15, 0.02, 0.06, 0.08);
    {
        std::vector<std::string> lab;
        for (int ip = 0; ip < kNPS; ++ip) lab.push_back(kPSs[ip]);
        TH1D* fr = ses::CatFrame(Nm("f").c_str(), lab, "free factors: data = a_{conv} #times conv + a_{dir} #times direct", -3.0, 3.5, 0.035, 0.035, 1.9);
        fr->Draw("axis");
        ses::HLine(-0.5, kNPS - 0.5, 0., kBlack, 2, 1);
        const int col[2] = {kColPeakConv, kColPeakDirect};
        for (int m = 0; m < 2; ++m)
            for (int k = 0; k < 2; ++k)
                for (int f = 0; f < 2; ++f) {
                    std::vector<double> x; std::vector<VE> y;
                    for (int ip = 0; ip < kNPS; ++ip) {
                        const sep::MixRes r = Mix(m ? vMain : vOff, ip, k);
                        x.push_back(ip - 0.3 + 0.2 * (2 * f + k) + 0.04 * m);
                        y.push_back(f == 0 ? VE(r.a, r.ea) : VE(r.b, r.eb));
                    }
                    ses::Pts(x, y, m ? col[f] : kColBefore, m ? (k ? 21 : 20) : (k ? 25 : 24), 1.6)->Draw("P same");
                }
        seen::Txt(0.18, 0.95, "#bf{two free factors (conversion, direct)}", 0.038);
        seen::Txt(0.18, 0.91, "green: a_{conv} after, blue: a_{dir} after; grey open: before", 0.030);
        seen::Txt(0.18, 0.88, "circles: OA classes, squares: M_{ee} classes", 0.030);
        seen::Txt(0.18, 0.12, "before: a_{conv} < 0 (unphysical); after: both > 0", 0.032);
        seen::Txt(0.18, 0.09, "factors include the absolute scale (L_{data} unconfirmed)", 0.030);
    }
    Save(cv, "rw_compare_peak_classes", true);
}
void FigPeakEff() {
    TCanvas* cv = new TCanvas(Nm("cPE").c_str(), "", 2800, 1800);
    cv->cd();
    seen::MkPad(Nm("h").c_str(), 0, 0.90, 1, 1, false, 0, 0, 0, 0);
    seen::Txt(0.01, 0.70, "Single-cut efficiencies of the #eta peak (base pre, all M_{4}): data (pass/fail P_{main}) against the simulated peak before and after the physics weights", 0.20);
    seen::Txt(0.01, 0.42, WNote().c_str(), 0.13);
    seen::Txt(0.01, 0.14, Form("data errors: stat #oplus fit-model spread; sim: event-group errors; M_{#pi#pi} cuts sculpted in data; #alpha_{#eta} > 160#circ data: fail fits p < 0.05; "
                               "OA > 2#circ pull %+.1f #rightarrow %+.1f (V3; V2 %+.1f): reduced, not removed",
                               Pull(gE[0].s[vOff][eTot].eps, gE[0].s[vOff][eTot].err, gE[0].d.e, gE[0].d.err), Pull(gE[0].s[vMain][eTot].eps, gE[0].s[vMain][eTot].err, gE[0].d.e, gE[0].d.err),
                               Pull(gE[0].s[vV2][eTot].eps, gE[0].s[vV2][eTot].err, gE[0].d.e, gE[0].d.err)), 0.13);
    std::vector<std::string> lab;
    for (int k = 0; k < kNEC; ++k) lab.push_back(kEC[k].lab);
    cv->cd();
    seen::MkPad(Nm("p").c_str(), 0, 0.32, 1, 0.90, false, 0.07, 0.01, 0.03, 0.10);
    TH1D* f1 = ses::CatFrame(Nm("f").c_str(), lab, "efficiency", 0.1, 1.08, 0.055, 0.050, 0.55);
    f1->Draw("axis");
    std::vector<double> xd, xb, xa, xdir, xconv;
    std::vector<VE> d, b, a, dir, conv, dirb;
    for (int k = 0; k < kNEC; ++k) {
        xd.push_back(k); xb.push_back(k - 0.22); xa.push_back(k + 0.22); xdir.push_back(k + 0.32); xconv.push_back(k + 0.12);
        d.push_back(VE(gE[k].d.e, gE[k].d.err));
        b.push_back(VE(gE[k].s[vOff][eTot].eps, gE[k].s[vOff][eTot].err)); a.push_back(VE(gE[k].s[vMain][eTot].eps, gE[k].s[vMain][eTot].err));
        dir.push_back(VE(gE[k].s[vMain][eDir].eps, gE[k].s[vMain][eDir].err)); conv.push_back(VE(gE[k].s[vMain][eConv].eps, gE[k].s[vMain][eConv].err));
        dirb.push_back(VE(gE[k].s[vOff][eDir].eps, gE[k].s[vOff][eDir].err));
    }
    TGraphErrors* gb = ses::Pts(xb, b, kColBefore, 25, 2.2); gb->Draw("P same");
    TGraphErrors* ga = ses::Pts(xa, a, kColAfter, 21, 2.2); ga->Draw("P same");
    std::vector<double> xdb; for (double x : xdir) xdb.push_back(x - 0.44);
    TGraphErrors* gdb = ses::Pts(xdb, dirb, kColPeakDirect, 26, 1.6); gdb->Draw("P same");
    TGraphErrors* gdir = ses::Pts(xdir, dir, kColPeakDirect, 22, 1.6); gdir->Draw("P same");
    TGraphErrors* gconv = ses::Pts(xconv, conv, kColPeakConv, 23, 1.6); gconv->Draw("P same");
    TGraphErrors* gd = ses::Pts(xd, d, kColData, 20, 2.4, 3); gd->Draw("P same");
    for (int k = 1; k < kNEC; ++k) { TLine* l = new TLine(k - 0.5, 0.1, k - 0.5, 1.08); l->SetLineColor(kGray); l->SetLineStyle(3); l->Draw(); }
    TLegend* lg = new TLegend(0.30, 0.13, 0.99, 0.30); lg->SetNColumns(3); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.035);
    lg->AddEntry(gd, "data peak", "pe"); lg->AddEntry(gb, "sim peak total, before", "pe"); lg->AddEntry(ga, "sim peak total, after", "pe");
    lg->AddEntry(gdb, "sim direct, before", "p"); lg->AddEntry(gdir, "sim direct, after", "p"); lg->AddEntry(gconv, "sim conversion, after", "p");
    lg->Draw();
    cv->cd();
    seen::MkPad(Nm("p").c_str(), 0, 0, 1, 0.32, false, 0.07, 0.01, 0.02, 0.22);
    TH1D* f2 = ses::CatFrame(Nm("f").c_str(), lab, "(sim - data)/#sigma", -7, 7, 0.10, 0.09, 0.30);
    f2->Draw("axis");
    ses::HLine(-0.5, kNEC - 0.5, 0, kBlack, 2, 1); ses::HLine(-0.5, kNEC - 0.5, 2, kGray + 1, 3, 1); ses::HLine(-0.5, kNEC - 0.5, -2, kGray + 1, 3, 1);
    std::vector<VE> pb, pa;
    for (int k = 0; k < kNEC; ++k) {
        pb.push_back(VE(Pull(b[k].first, b[k].second, d[k].first, d[k].second), 0));
        pa.push_back(VE(Pull(a[k].first, a[k].second, d[k].first, d[k].second), 0));
    }
    ses::Pts(xb, pb, kColBefore, 25, 2.2)->Draw("P same");
    ses::Pts(xa, pa, kColAfter, 21, 2.2)->Draw("P same");
    Save(cv, "rw_compare_peak_eff", true);
}
void FigPeakPos() {
    TCanvas* cv = new TCanvas(Nm("cPP").c_str(), "", 2700, 1000);
    cv->cd();
    seen::MkPad(Nm("h").c_str(), 0, 0.88, 1, 1, false, 0, 0, 0, 0);
    seen::Txt(0.01, 0.55, "Position and width of the simulated #eta peak before and after the physics weights (Gaussian #pm2#sigma, 1 MeV bins), with the data (FitPeakFree)", 0.30);
    seen::Txt(0.01, 0.15, "simulated M_{4} not rescaled (K0S scale sim/data 1.0048: 0.5459 in simulation = 0.5433 in data units); data errors statistical only", 0.22);
    cv->cd();
    seen::MkPad(Nm("p").c_str(), 0, 0, 0.34, 0.88, false, 0.15, 0.02, 0.06, 0.13);
    TH1D* fr = ses::Frame(Nm("f").c_str(), 100, 0.51, 0.58, "M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}]", "fraction per MeV", 0, 0.085, 0.05, 1.4);
    fr->Draw("axis");
    const int col[3] = {kColPeakDirect, kColPeakConv, kColPeakTotal};
    for (int m = 0; m < 2; ++m)
        for (int c = 0; c < 3; ++c) {
            TH1D* h = (TH1D*)gPosH[m][c]->Clone(Nm("h").c_str());
            h->Scale(1. / h->Integral());
            seen::StyleH(h, col[c], 3, m ? 1 : 2);
            h->Draw("HIST same");
        }
    seen::Txt(0.18, 0.88, "#bf{pre}: dashed before, solid after", 0.05);
    seen::Txt(0.18, 0.82, "#color[601]{direct}, #color[418]{conversion}, #color[435]{total}", 0.05);
    for (int q = 0; q < 2; ++q) {
        cv->cd();
        seen::MkPad(Nm("p").c_str(), 0.34 + q * 0.33, 0, 0.67 + q * 0.33, 0.88, false, 0.17, 0.02, 0.06, 0.13);
        std::vector<std::string> lab(kPSs, kPSs + kNPS);
        TH1D* f = q == 0 ? ses::CatFrame(Nm("f").c_str(), lab, "peak position [GeV/c^{2}]", 0.540, 0.549, 0.05, 0.05, 1.6)
                         : ses::CatFrame(Nm("f").c_str(), lab, "peak width #sigma [MeV/c^{2}]", 3, 9, 0.05, 0.05, 1.4);
        f->Draw("axis");
        std::vector<double> xd, xb, xa;
        std::vector<VE> d, b, a;
        for (int ip = 0; ip < kNPS; ++ip) {
            xd.push_back(ip); xb.push_back(ip - 0.2); xa.push_back(ip + 0.2);
            const scan::FreeFit& ff = gPD[ip].pos; const sep::GFit& gb = gPos[0][ip].g[2]; const sep::GFit& ga = gPos[1][ip].g[2];
            if (q == 0) { d.push_back(VE(ff.mean, ff.eMean)); b.push_back(VE(gb.m, gb.em)); a.push_back(VE(ga.m, ga.em)); }
            else { d.push_back(VE(ff.sigma * 1e3, ff.eSigma * 1e3)); b.push_back(VE(gb.s * 1e3, gb.es * 1e3)); a.push_back(VE(ga.s * 1e3, ga.es * 1e3)); }
        }
        ses::Pts(xb, b, kColBefore, 25, 2.0)->Draw("P same");
        ses::Pts(xa, a, kColAfter, 21, 2.0)->Draw("P same");
        ses::Pts(xd, d, kColData, 20, 2.0, 3)->Draw("P same");
        seen::Txt(0.20, 0.88, q == 0 ? "black: data; grey open: sim total before; teal: after" : "width of the simulated total and of the data peak", 0.042);
    }
    Save(cv, "rw_compare_peak_position", false);
}

// =====================================================================================================
// Part 2 and 3 share the selections: 0 pre (reference), 1 E1, 2 E2, 3 E3, 4 E3oa, 5 E4, 6 E4oa
// =====================================================================================================
const int kNES = 7;
const char* kEST[kNES] = {"pre", "E1", "E2", "E3", "E3oa", "E4", "E4oa"};
Sel ESim(int i) { return i == 0 ? SelPreset(kPre) : seen::SimSel(i - 1); }
scan::Sel EDat(int i) { return i == 0 ? DataSelPreset(kPre) : seen::DatSel(i - 1); }
std::string ELatex(int i) { return i == 0 ? std::string("pre: isBest #geq 0, PT3, v_{z} > -160 mm, cor_ok") : seen::SelLatex(i - 1); }

// ---- part 2: bump
struct BumpSim { WSum all[3], win[3]; };   // D, C, D+C; M4 < 1 and bump window
BumpSim gBS[kNV][kNES];
const int kNMP = 3;
const int kMPSel[kNMP] = {1, 2, 0};      // E1, E2, pre
std::vector<seen::DH> gMppD[kNMP];
std::vector<WSpec> gMppS[kNV][kNMP];
std::vector<WSpec> gBumpM4[kNV][kNES];
Binning MppBinning() { Binning b; b.nbin = 20; b.lo = 0.27; b.hi = 0.47; b.x = [](const Row& r) { return (double)r.mpipi; }; return b; }
void BumpPart() {
    TStopwatch sw;
    const std::vector<Comp> cb = {CompBumpDalitz(), CompBumpConv(), CBump()};
    for (int iv : {vOff, vMain, vSl}) {
        const SimOpt o = Opt(iv);
        for (int is = 0; is < kNES; ++is) {
            const std::vector<WSum> a = SumMany(*gS, SWin(ESim(is), -1e9, 1.0), cb, o);
            const std::vector<WSum> w = SumMany(*gS, SWin(ESim(is), scan::kBumpLo, scan::kBumpHi), cb, o);
            for (int c = 0; c < 3; ++c) { gBS[iv][is].all[c] = a[c]; gBS[iv][is].win[c] = w[c]; }
            Binning b; b.nbin = 80; b.lo = 0.40; b.hi = 0.60;
            gBumpM4[iv][is] = FillMany(*gS, ESim(is), cb, o, b);
        }
        for (int j = 0; j < kNMP; ++j) {
            const int is = kMPSel[j];
            // components: 0 T1 + T2, 1 bump D+C, 2 variant (a), 3 rest = T1 + T2 without the bump (its own event-group variance)
            const Comp cs = CompSim(), cb = CBump();
            const Comp rest("rest", [cs, cb](const Row& r) { return cs.f(r) && !cb.f(r); });
            gMppS[iv][j] = FillMany(*gD, SWin(ESim(is), scan::kBumpLo, scan::kBumpHi), {cs, cb, CompSimVar(1), rest}, o, MppBinning());
        }
    }
    for (int j = 0; j < kNMP; ++j) {
        const Binning b = MppBinning();
        seen::DVar xm = [](const scan::Row& r) { return (double)r.mpipi; };
        gMppD[j] = seen::FillDH(*gDD, DWin(EDat(kMPSel[j]), scan::kBumpLo, scan::kBumpHi), {xm}, {b});
    }
    printf("BumpPart: %.1f s\n", sw.RealTime());
}
// M(pipi) in the bump window: fraction above 0.34, mean, shape chi2 (simulation normalised to the data), absolute ratio
struct MppStat { double fD, fDe, fS, meanD, meanS, chi2; int ndf; double ratio, ratioE; };
MppStat MppCompare(int j, int iv, int comp) {
    MppStat o;
    const seen::DH& d = gMppD[j][0];
    const WSpec& s = gMppS[iv][j][comp];
    double D = 0, DV = 0, Da = 0, DaV = 0, S = 0, SV = 0, Sa = 0, mD = 0, mS = 0;
    for (int i = 0; i < d.nb; ++i) {
        const double x = d.lo + (i + 0.5) * d.bw;
        D += d.sig(i); DV += d.sigVar(i); S += s.s[i]; SV += s.v[i];
        mD += x * d.sig(i); mS += x * s.s[i];
        if (x > 0.34) { Da += d.sig(i); DaV += d.sigVar(i); Sa += s.s[i]; }
    }
    o.fD = Da / D; o.fDe = std::sqrt((1 - o.fD) * (1 - o.fD) * DaV + o.fD * o.fD * (DV - DaV)) / D;
    o.fS = Sa / S; o.meanD = mD / D; o.meanS = mS / S;
    const double k = D / (S * kToRows);
    o.chi2 = 0; o.ndf = -1;
    for (int i = 0; i < d.nb; ++i) {
        const double v = d.sigVar(i) + k * k * s.v[i] * kToRows * kToRows;
        if (!(v > 0) || (d.n[i] < 5 && s.s[i] <= 0)) continue;
        o.chi2 += std::pow(d.sig(i) - k * s.s[i] * kToRows, 2) / v; ++o.ndf;
    }
    o.ratio = D / (S * kToRows); o.ratioE = o.ratio * std::sqrt(DV / (D * D) + SV / (S * S));
    return o;
}

// two-template fit of the data M(pipi) in the bump window: data_i = a x bump_i + c x rest_i (rest = T1 + T2 - bump); variance =
// data (N + CB) + a^2 var(bump) + c^2 var(rest), iterated. Review of 2026-11-12: the earlier version used var(T1 + T2) for the rest,
// which already contains var(bump): the bump template variance was counted twice and the chi2 came out too small.
struct MppFitRes { double a = 0, ea = 0, c = 0, ec = 0, rho = 0, chi2 = 0; int ndf = 0; };
MppFitRes MppFit(int j, int iv) {
    MppFitRes o;
    const seen::DH& d = gMppD[j][0];
    const WSpec& B = gMppS[iv][j][1]; const WSpec& R = gMppS[iv][j][3];
    double a = 0.5, c = 1.0, Saa = 0, Sac = 0, Scc = 0;
    for (int it = 0; it < 50; ++it) {
        double Sya = 0, Syc = 0; Saa = Sac = Scc = 0;
        for (int i = 0; i < d.nb; ++i) {
            const double b = B.s[i] * kToRows, r = R.s[i] * kToRows;
            const double v = d.sigVar(i) + (a * a * B.v[i] + c * c * R.v[i]) * kToRows * kToRows;
            if (!(v > 0) || (b <= 0 && r <= 0)) continue;
            Saa += b * b / v; Sac += b * r / v; Scc += r * r / v; Sya += d.sig(i) * b / v; Syc += d.sig(i) * r / v;
        }
        const double det = Saa * Scc - Sac * Sac;
        if (!(det > 0)) return o;
        a = (Sya * Scc - Syc * Sac) / det; c = (Saa * Syc - Sac * Sya) / det;
    }
    const double det = Saa * Scc - Sac * Sac;
    o.a = a; o.c = c; o.ea = std::sqrt(Scc / det); o.ec = std::sqrt(Saa / det); o.rho = -Sac / det / (o.ea * o.ec);
    int n = 0;
    for (int i = 0; i < d.nb; ++i) {
        const double b = B.s[i] * kToRows, r = R.s[i] * kToRows;
        const double v = d.sigVar(i) + (a * a * B.v[i] + c * c * R.v[i]) * kToRows * kToRows;
        if (!(v > 0) || (b <= 0 && r <= 0)) continue;
        o.chi2 += std::pow(d.sig(i) - a * b - c * r, 2) / v; ++n;
    }
    o.ndf = n - 2;
    return o;
}

// ---- part 3: eta-enriched spectra, regions, fits, peak methods
std::vector<WSpec> eLog[2][kNES], eLin[2][kNES], eReg[2][kNES];
DSpec dLog[kNES], dLin[kNES];
std::vector<scan::Spectra> dReg;           // win[2 + k]: region k (0..7); win[10], win[11]: edge bands
const double kBand[2][2] = {{0.4875, 0.5225}, {0.56, 0.595}};
WSum gBand[kNV][kNES][2][2];               // [variant][sel][band][T1+T2, (a)]
bool gHaveBand[kNV] = {};
void EnrichedSpectra() {
    TStopwatch sw;
    const std::vector<Comp> comps = seen::Comps();
    Binning bLog; bLog.nbin = 130; bLog.lo = 0.35; bLog.hi = 1.00;
    Binning bLin; bLin.nbin = 120; bLin.lo = 0.40; bLin.hi = 0.70;
    const Binning bReg = seen::RegBinning();
    for (int m = 0; m < 2; ++m) {
        const SimOpt o = Opt(m ? vMain : vOff);
        for (int is = 0; is < kNES; ++is) {
            eLog[m][is] = FillMany(*gD, ESim(is), comps, o, bLog);
            eLin[m][is] = FillMany(*gD, ESim(is), comps, o, bLin);
            eReg[m][is] = FillMany(*gD, ESim(is), comps, o, bReg);
        }
    }
    std::vector<scan::Sel> dsels;
    for (int is = 0; is < kNES; ++is) {
        dLog[is] = DataSpectrum(*gDD, EDat(is), bLog.nbin, bLog.lo, bLog.hi);
        dLin[is] = DataSpectrum(*gDD, EDat(is), bLin.nbin, bLin.lo, bLin.hi);
        dsels.push_back(EDat(is));
    }
    scan::FillOpt fo; fo.nbin = 260; fo.lo = 0.35; fo.hi = 1.00; fo.reb = 4;
    for (int k = 0; k <= seen::kNReg; ++k) fo.windows.push_back(std::make_pair(seen::kReg[k].lo, seen::kReg[k].hi));
    for (int b = 0; b < 2; ++b) fo.windows.push_back(std::make_pair(kBand[b][0], kBand[b][1]));
    dReg = scan::FillMany(*gDD, dsels, fo);
    printf("EnrichedSpectra: %.1f s\n", sw.RealTime());
}
void BandSums(int iv) {
    const SimOpt o = Opt(iv);
    for (int is = 0; is < kNES; ++is)
        for (int b = 0; b < 2; ++b) {
            const std::vector<WSum> w = SumMany(*gD, SWin(ESim(is), kBand[b][0], kBand[b][1]), {CompSim(), CompSimVar(1)}, o);
            gBand[iv][is][b][0] = w[0]; gBand[iv][is][b][1] = w[1];
        }
    gHaveBand[iv] = true;
}
// data/sim in a region k (0..7 regions of seen::kReg) for spectra set m, ghost variant g
VE RegRatio(int m, int is, int k, int g) {
    const scan::WinSum& w = dReg[is].win[2 + k];
    const int c = g ? seen::cSumA : seen::cSum;
    const WSpec& s = eReg[m][is][c];
    const double ss = k < seen::kNReg ? s.s[k] : s.S, sv = k < seen::kNReg ? s.v[k] : s.V;
    return Ratio(w.sig() / kLdata, w.sigErr() / kLdata, ss / kLsim, std::sqrt(sv) / kLsim);
}
VE BandRatio(int iv, int is, int b, int g) {
    const scan::WinSum& w = dReg[is].win[2 + seen::kNReg + 1 + b];
    const WSum& s = gBand[iv][is][b][g];
    return Ratio(w.sig() / kLdata, w.sigErr() / kLdata, s.Xsec(), s.XsecErr());
}
VE Jump(int iv, int is, int g) {
    const VE l = BandRatio(iv, is, 0, g), r = BandRatio(iv, is, 1, g);
    return VE(r.first / l.first, r.first / l.first * std::hypot(l.second / l.first, r.second / r.first));
}

// ---- template fits
struct FitSet { seen::FitRes m1[2][2][2]; seen::FitRes m3[2]; bool ok = false; };   // M1 [range][scale][ghost]; M3 [ghost] (range 0, no scale)
FitSet gFit[kNV][kNES];
DSpec gFitDS[kNES][2];
struct Rng { double lo, hi; int nb; };
const Rng kRng[2] = {{0.40, 0.75, 70}, {0.43, 0.70, 54}};
seen::FitRes TFit(const DSpec& DS, const std::vector<WSpec>& TS, const std::vector<int>& grp, int nb) {
    std::vector<int> bins;
    for (int i = 0; i < nb; ++i) bins.push_back(i);
    std::vector<std::vector<double>> TT, TV;
    for (int k : grp) {
        std::vector<double> t(nb), tv(nb);
        for (int i = 0; i < nb; ++i) { t[i] = TS[k].s[i] * kToRows; tv[i] = TS[k].v[i] * kToRows * kToRows; }
        TT.push_back(t); TV.push_back(tv);
    }
    return seen::TemplateFit(DS.sp, bins, TT, TV);
}
void FitDataPrep() {
    for (int is = 0; is < kNES; ++is)
        for (int ir = 0; ir < 2; ++ir) gFitDS[is][ir] = DataSpectrum(*gDD, EDat(is), kRng[ir].nb, kRng[ir].lo, kRng[ir].hi);
}
void Fits(int iv, int isOnly = -1) {
    const SimOpt o = Opt(iv);
    const std::vector<Comp> fc = seen::FitComps();
    for (int is = 0; is < kNES; ++is) {
        if (isOnly >= 0 && is != isOnly) continue;
        FitSet& F = gFit[iv][is];
        for (int ir = 0; ir < 2; ++ir)
            for (int sc = 0; sc < 2; ++sc) {
                Binning b; b.nbin = kRng[ir].nb; b.lo = kRng[ir].lo; b.hi = kRng[ir].hi;
                const double f = sc ? seen::kScaleK0S : 1.;
                b.x = [f](const Row& r) { return (double)r.m4 * f; };
                const std::vector<WSpec> TS = FillMany(*gD, ESim(is), fc, o, b);
                for (int g = 0; g < 2; ++g) {
                    F.m1[ir][sc][g] = TFit(gFitDS[is][ir], TS, {seen::fBump, seen::fPeak, g ? seen::fKA : seen::fK}, b.nbin);
                    if (ir == 0 && sc == 0)
                        F.m3[g] = TFit(gFitDS[is][ir], TS, {seen::fBump, seen::fPeak, g ? seen::fKEA : seen::fKE, g ? seen::fKMA : seen::fKM}, b.nbin);
                }
            }
        F.ok = true;
    }
}
VE Scaled(const seen::FitRes& r, int k) { return VE(r.a[k], r.ea[k] * std::sqrt(std::max(1., r.chi2 / std::max(1, r.ndf)))); }
void M1Range(const FitSet& F, int k, double& lo, double& hi) {
    lo = 1e9; hi = -1e9;
    for (int ir = 0; ir < 2; ++ir)
        for (int sc = 0; sc < 2; ++sc)
            for (int g = 0; g < 2; ++g) {
                if (!F.m1[ir][sc][g].ok) continue;
                lo = std::min(lo, F.m1[ir][sc][g].a[k]); hi = std::max(hi, F.m1[ir][sc][g].a[k]);
            }
}

// ---- peak data/sim by method (copy of ses::PartF with the physics option; reconcile.md / eta_summary_peak_methods.txt)
const int kNMS = 4;
const char* kMST[kNMS] = {"pre", "pre+OA>2", "E1", "E2"};
Sel MSim(int i) { return i == 0 ? SelPreset(kPre) : i == 1 ? SelPreset(kPreOA2) : seen::SimSel(i == 2 ? seen::sE1 : seen::sE2); }
scan::Sel MDat(int i) { return i == 0 ? DataSelPreset(kPre) : i == 1 ? DataSelPreset(kPreOA2) : seen::DatSel(i == 2 ? seen::sE1 : seen::sE2); }
struct MethData { sep::DFit df; double D[3], DV[3]; scan::Spectra sp; DSpec DS[2]; bool ok = false; };
MethData gMD[kNMS];
struct MethRes {
    VE A, Aw, B, C, Ca, D, E, Ea; double mA = 0, mAw = 0, mB = 0, mC = 0, mCa = 0;
    double Dlo = 0, Dhi = 0, Dchi = 0; int Dndf = 0; double Dscaled = 0;
    double Elo = 0, Ehi = 0, rL = 0, rR = 0, rW = 0;
    double simAll = 0, simAllE = 0, simWin = 0, pkOnly = 0, simTot = 0, simTotA = 0;
    bool ok = false;
};
MethRes gM[kNV][kNMS];
void MethDataPrep() {
    for (int is = 0; is < kNMS; ++is) {
        MethData& md = gMD[is];
        scan::FillOpt fo;
        fo.windows = {{0.5225, scan::kPeakLo}, {scan::kPeakLo, scan::kPeakHi}, {scan::kPeakHi, 0.5650}};
        md.sp = scan::FillMany(*gDD, {MDat(is)}, fo)[0];
        md.df = sep::DoFit(md.sp);
        for (int b = 0; b < 3; ++b) { md.D[b] = md.sp.win[2 + b].sig(); md.DV[b] = std::pow(md.sp.win[2 + b].sigErr(), 2); }
        for (int ir = 0; ir < 2; ++ir) md.DS[ir] = DataSpectrum(*gDD, MDat(is), kRng[ir].nb, kRng[ir].lo, kRng[ir].hi);
        md.ok = true;
    }
}
MethRes Methods(int is, int iv) {
    MethRes m;
    const MethData& md = gMD[is];
    const SimData& d = *gD;
    const SimOpt opt = Opt(iv);
    const double sc = seen::kScaleK0S;
    static const std::vector<scan::PeakCfg> vc = ses::VarCfgs();
    const Sel ssel = MSim(is);
    const double dS = md.df.S, dSe = md.df.eS, dSm = md.df.eM;
    const Comp cPk = CompPeakTotal();
    const Comp cT = CompSim(), cTa = CompSimVar(1);
    const Comp cB{"sim without peak", [](const Row& r) { return ((r.tmpl == 1 && !X17Row(r)) || r.tmpl == 2) && !(r.tmpl == 1 && (r.cls == 4 || r.cls == 5)); }};
    const Comp cBa{"(a) without peak", [cTa](const Row& r) { return cTa.f(r) && !(r.tmpl == 1 && (r.cls == 4 || r.cls == 5)); }};
    const std::vector<Comp> cc = {cPk, cT, cTa, cB, cBa};
    enum { kPk = 0, kT, kTa, kB, kBa };
    const std::vector<WSum> all = SumMany(d, ssel, {cPk}, opt);
    m.simAll = all[0].S * kToRows; m.simAllE = std::sqrt(all[0].V) * kToRows;
    m.A = Ratio(dS, dSe, m.simAll, m.simAllE); m.mA = m.A.first * dSm / dS;
    Binning bP; bP.nbin = md.sp.nbin; bP.lo = md.sp.lo; bP.hi = md.sp.hi;
    bP.x = [sc](const Row& r) { return r.m4 * sc; };
    const std::vector<WSpec> Tm = FillMany(d, ssel, cc, opt, bP);
    std::vector<std::vector<double>> y(cc.size(), std::vector<double>(bP.nbin)), v = y;
    for (size_t c = 0; c < cc.size(); ++c)
        for (int i = 0; i < bP.nbin; ++i) { y[c][i] = Tm[c].s[i] * kToRows; v[c][i] = Tm[c].v[i] * kToRows * kToRows; }
    {
        const double fw = scan::GaussFrac(scan::kPeakLo, scan::kPeakHi, scan::kPeakMean, scan::kPeakSigma);
        const Sel inW = SelAnd(ssel, [sc](const Row& r) { const double x = r.m4 * sc; return x >= scan::kPeakLo && x < scan::kPeakHi; });
        const std::vector<WSum> w = SumMany(d, inW, {cPk}, opt);
        m.simWin = w[0].S * kToRows;
        m.Aw = Ratio(dS * fw, dSe * fw, w[0].S * kToRows, std::sqrt(w[0].V) * kToRows); m.mAw = m.Aw.first * dSm / dS;
    }
    auto pmRatio = [&](int comp, VE& out, double& model, double* simS) {
        double e = 0;
        const double s0 = ses::PmainCfg(y[comp], v[comp], bP.lo, bP.bw(), scan::PeakCfg(), &e);
        if (simS) *simS = s0;
        out = Ratio(dS, dSe, s0, e);
        double lo = 1e30, hi = -1e30; int n = 0;
        for (size_t k = 0; k < vc.size() && k < md.df.vS.size(); ++k) {
            if (!md.df.vOk[k]) continue;
            const double sv = ses::PmainCfg(y[comp], v[comp], bP.lo, bP.bw(), vc[k]);
            if (!(sv > 0)) continue;
            const double r = md.df.vS[k] / sv; lo = std::min(lo, r); hi = std::max(hi, r); ++n;
        }
        model = n >= 2 ? 0.5 * (hi - lo) : out.first * dSm / dS;
    };
    pmRatio(kPk, m.B, m.mB, &m.pkOnly);
    pmRatio(kT, m.C, m.mC, &m.simTot);
    pmRatio(kTa, m.Ca, m.mCa, &m.simTotA);
    // method E: side bands (K0S-scaled simulation)
    {
        const double ed[4] = {0.5225, scan::kPeakLo, scan::kPeakHi, 0.5650};
        double Bk[2][3], BkV[2][3], Pk[3], PkV[3];
        for (int b = 0; b < 3; ++b) {
            const double lo = ed[b], hi = ed[b + 1];
            const Sel sb = SelAnd(ssel, [sc, lo, hi](const Row& r) { const double x = r.m4 * sc; return x >= lo && x < hi; });
            const std::vector<WSum> w = SumMany(d, sb, {cB, cBa, cPk}, opt);
            for (int g = 0; g < 2; ++g) { Bk[g][b] = w[g].S * kToRows; BkV[g][b] = w[g].V * kToRows * kToRows; }
            Pk[b] = w[2].S * kToRows; PkV[b] = w[2].V * kToRows * kToRows;
        }
        const double cL = 0.5 * (ed[0] + ed[1]), cW = 0.5 * (ed[1] + ed[2]), cR = 0.5 * (ed[2] + ed[3]);
        const double a = (cR - cW) / (cR - cL);
        const double* D = md.D; const double* DV = md.DV;
        auto solve = [&](int g, int mode, double* rl, double* rr, double* rw) {
            double k = 0.5, rL = 1, rR = 1, rW = 1;
            for (int it = 0; it < 200; ++it) {
                rL = (D[0] - k * Pk[0]) / Bk[g][0]; rR = (D[2] - k * Pk[2]) / Bk[g][2];
                rW = mode == 0 ? a * rL + (1 - a) * rR : (mode == 1 ? rL : rR);
                k = (D[1] - rW * Bk[g][1]) / Pk[1];
            }
            if (rl) *rl = rL;
            if (rr) *rr = rR;
            if (rw) *rw = rW;
            const double vL = (DV[0] + rL * rL * BkV[g][0]) / (Bk[g][0] * Bk[g][0]), vR = (DV[2] + rR * rR * BkV[g][2]) / (Bk[g][2] * Bk[g][2]);
            const double vW = mode == 0 ? a * a * vL + (1 - a) * (1 - a) * vR : (mode == 1 ? vL : vR);
            const double ek = std::sqrt(DV[1] + rW * rW * BkV[g][1] + Bk[g][1] * Bk[g][1] * vW + k * k * PkV[1]) / Pk[1];
            return VE(k, ek);
        };
        m.E = solve(0, 0, &m.rL, &m.rR, &m.rW);
        const double k1 = solve(0, 1, nullptr, nullptr, nullptr).first, k2 = solve(0, 2, nullptr, nullptr, nullptr).first;
        m.Elo = std::min(k1, k2); m.Ehi = std::max(k1, k2);
        m.Ea = solve(1, 0, nullptr, nullptr, nullptr);
    }
    // method D: template fit M1 (2 ranges x 2 scales x 2 ghost variants)
    {
        const std::vector<Comp> fc = seen::FitComps();
        m.Dlo = 1e9; m.Dhi = -1e9;
        for (int ir = 0; ir < 2; ++ir)
            for (int scl = 0; scl < 2; ++scl) {
                Binning b; b.nbin = kRng[ir].nb; b.lo = kRng[ir].lo; b.hi = kRng[ir].hi;
                const double f = scl ? sc : 1.;
                b.x = [f](const Row& r) { return (double)r.m4 * f; };
                const std::vector<WSpec> TS = FillMany(d, ssel, fc, opt, b);
                for (int gv = 0; gv < 2; ++gv) {
                    const seen::FitRes r = TFit(md.DS[ir], TS, {seen::fBump, seen::fPeak, gv ? seen::fKA : seen::fK}, b.nbin);
                    if (!r.ok) continue;
                    m.Dlo = std::min(m.Dlo, r.a[1]); m.Dhi = std::max(m.Dhi, r.a[1]);
                    if (ir == 0 && scl == 0 && gv == 0) {
                        m.D = VE(r.a[1], r.ea[1]); m.Dchi = r.chi2; m.Dndf = r.ndf;
                        m.Dscaled = r.ea[1] * std::sqrt(std::max(1., r.chi2 / r.ndf));
                    }
                }
            }
    }
    m.ok = true;
    return m;
}

void EnrichedTables() {
    OpenT("rw_compare_enriched.txt");
    T("# rw_compare_enriched: data N - CB against the simulation in the eta-enriched selections (M4 < 1 GeV/c2), BEFORE and AFTER the physics weights\n");
    T("# simulation: std + lep + eta owned rows (mode A), w (PDG) x physics weight (main: direct V3 = V2 x decay angles, bump KLOE, conversion F_V,\n");
    T("#   background rules);\n");
    T("# NOT physics results: every absolute data/sim ratio and alpha below scales with 6080.5 / L_data (UNCONFIRMED); the eta peak data/sim depends on\n");
    T("#   the method (section 4: A-E span a factor of about 3.7 after the weights); the template fits are poor (chi2/ndf 2.4-3.6)\n");
    T("#   T1 + T2 (default) and (a) = + ghost-flagged true pairs; sigma = sum w / L_sim (221.05 nb^-1); no CB in simulation\n");
    T("# data: (N - CB) / L_data, CB = 2 sqrt(N++ N--) in 10 MeV bins; L_data = 6080.5 nb^-1 UNCONFIRMED (every ratio and alpha scales with 6080.5 / L_data)\n");
    T("# pre = isBest >= 0, PT3, vz > -160 mm, cor_ok (reference, not enriched); E1 = slides (oa_lab < 50, eta_rest > 140, M(pipi) < 0.46, OA > 2);\n");
    T("#   E2 = E1 without OA > 2; E3 = old code (M(pipi) < 0.42, no OA); E3oa = E3 + OA > 2; E4 = eta band (no M(pipi) cut); E4oa = E4 + OA > 2\n");
    T("# regions above 0.56 are 80-90 %% multi-pion background, which SMASH does not describe: listed, not used as validation\n\n");
    T("## 1. data / simulation per M4 region (T1 + T2; (a) in brackets), before -> after\n");
    T("%-5s %-27s %16s %16s | %-26s %-26s\n", "sel", "region", "data N-CB [nb]", "sim bef->aft [nb]", "data/sim before (a)", "data/sim after (a)");
    for (int is = 0; is < kNES; ++is) {
        for (int k = 0; k <= seen::kNReg; ++k) {
            const VE rb = RegRatio(0, is, k, 0), ra = RegRatio(1, is, k, 0), rba = RegRatio(0, is, k, 1), raa = RegRatio(1, is, k, 1);
            const scan::WinSum& w = dReg[is].win[2 + k];
            const WSpec& sb = eReg[0][is][seen::cSum]; const WSpec& sa = eReg[1][is][seen::cSum];
            const double xb = (k < seen::kNReg ? sb.s[k] : sb.S) / kLsim, xa = (k < seen::kNReg ? sa.s[k] : sa.S) / kLsim;
            T("%-5s %-27s %8.4f+-%-7.4f %7.3f -> %-7.3f | %6.3f +- %-5.3f (%5.3f)      %6.3f +- %-5.3f (%5.3f)\n", kEST[is], seen::kReg[k].name, w.sig() / kLdata,
              w.sigErr() / kLdata, xb, xa, rb.first, rb.second, rba.first, ra.first, ra.second, raa.first);
        }
        for (int b = 0; b < 2; ++b) {
            const VE rb = BandRatio(vOff, is, b, 0), ra = BandRatio(vMain, is, b, 0);
            T("%-5s %-27s %16s %16s | %6.3f +- %-5.3f (%5.3f)      %6.3f +- %-5.3f (%5.3f)\n", kEST[is], Form("band %.4f-%.4f", kBand[b][0], kBand[b][1]), "", "",
              rb.first, rb.second, BandRatio(vOff, is, b, 1).first, ra.first, ra.second, BandRatio(vMain, is, b, 1).first);
        }
        const VE jb = Jump(vOff, is, 0), ja = Jump(vMain, is, 0);
        T("%-5s %-27s %16s %16s | %6.3f +- %-5.3f (%5.3f)      %6.3f +- %-5.3f (%5.3f)\n\n", kEST[is], "JUMP = band 0.56 / band 0.49", "", "", jb.first, jb.second,
          Jump(vOff, is, 1).first, ja.first, ja.second, Jump(vMain, is, 1).first);
    }
    T("## 2. simulated component shares of T1 + T2 per region, after (multi-pion, eta-related bkg, bump Dalitz, bump conv., peak conv., peak direct)\n");
    for (int is = 1; is < kNES; ++is) {
        T("%-5s", kEST[is]);
        for (int k : {0, 1, 3, 5}) {
            const WSpec& s = eReg[1][is][seen::cSum];
            T(" | %s:", seen::kReg[k].tag);
            for (int c : seen::kStack) T(" %.3f", eReg[1][is][c].s[k] / s.s[k]);
        }
        T("\n");
    }
    T("\n## 3. template fits of data N - CB (M1: bump, peak (direct + conv., simulated ratio), physical bkg; M3: bkg split eta-related / multi-pion)\n");
    T("# nominal: 0.40-0.75, 5 MeV bins, no mass scale, default ghost treatment; errors x sqrt(chi2/ndf); range = over 2 fit ranges x 2 mass scales\n");
    T("#   (none, K0S 0.9952) x 2 ghost variants. alpha = 1: SMASH (+ weights) at L_data = 6080.5 nb^-1 (UNCONFIRMED). These factors compare shapes, NOT yields\n");
    T("%-5s %-6s | %-17s %-17s %-17s %-12s | %-13s %-13s %-13s | %-17s %-17s %-17s %-17s %-12s\n", "sel", "", "M1 a_bump", "M1 a_peak", "M1 a_bkg", "chi2/ndf",
      "range a_bump", "range a_peak", "range a_bkg", "M3 a_bump", "M3 a_peak", "M3 a_bkgEta", "M3 a_bkgMP", "M3 chi2/ndf");
    for (int is = 0; is < kNES; ++is)
        for (int m = 0; m < 2; ++m) {
            const FitSet& F = gFit[m ? vMain : vOff][is];
            const seen::FitRes& r = F.m1[0][0][0]; const seen::FitRes& r3 = F.m3[0];
            double lo[3], hi[3];
            for (int k = 0; k < 3; ++k) M1Range(F, k, lo[k], hi[k]);
            T("%-5s %-6s | %6.3f +- %-6.3f %6.3f +- %-6.3f %6.3f +- %-6.3f %6.1f/%-5d | %5.3f-%-7.3f %5.3f-%-7.3f %5.3f-%-7.3f | %6.3f +- %-6.3f %6.3f +- %-6.3f %6.3f +- %-6.3f "
              "%6.3f +- %-6.3f %6.1f/%-5d\n",
              kEST[is], m ? "after" : "before", Scaled(r, 0).first, Scaled(r, 0).second, Scaled(r, 1).first, Scaled(r, 1).second, Scaled(r, 2).first,
              Scaled(r, 2).second, r.chi2, r.ndf, lo[0], hi[0], lo[1], hi[1], lo[2], hi[2], Scaled(r3, 0).first, Scaled(r3, 0).second, Scaled(r3, 1).first,
              Scaled(r3, 1).second, Scaled(r3, 2).first, Scaled(r3, 2).second, Scaled(r3, 3).first, Scaled(r3, 3).second, r3.chi2, r3.ndf);
        }
    T("\n## 4. eta peak, data / simulation by method (reconcile.md; as tables/eta_summary_peak_methods.txt), before -> after\n");
    T("# A  data P_main / sim peak all M4; Aw  inside the window; B  P_main on data and on the sim peak template; C  P_main on data and on the sim\n");
    T("#   T1 + T2 (C(a): variant (a)); D  template fit M1 alpha_peak (error x sqrt(chi2/ndf), range over ranges / scales / ghost); E  side bands\n");
    T("#   (range: left or right band alone; E(a): variant (a)). Errors: stat +- model (48 data fit variants); add 20 %% absolute scale (L_data unconfirmed)\n");
    T("%-9s %-5s | %-28s | %-28s | %s\n", "sel", "meth", "before", "after", "notes (after)");
    for (int is = 0; is < kNMS; ++is) {
        const MethRes& b = gM[vOff][is]; const MethRes& a = gM[vMain][is];
        T("%-9s %-5s | %6.3f +- %-5.3f +- %-5.3f      | %6.3f +- %-5.3f +- %-5.3f      | sim peak all M4 %.0f -> %.0f rows (data %.0f +- %.0f +- %.0f)\n", kMST[is], "A",
          b.A.first, b.A.second, b.mA, a.A.first, a.A.second, a.mA, b.simAll, a.simAll, gMD[is].df.S, gMD[is].df.eS, gMD[is].df.eM);
        T("%-9s %-5s | %6.3f +- %-5.3f +- %-5.3f      | %6.3f +- %-5.3f +- %-5.3f      | sim peak in window %.0f -> %.0f\n", kMST[is], "Aw", b.Aw.first, b.Aw.second, b.mAw,
          a.Aw.first, a.Aw.second, a.mAw, b.simWin, a.simWin);
        T("%-9s %-5s | %6.3f +- %-5.3f +- %-5.3f      | %6.3f +- %-5.3f +- %-5.3f      | P_main on the sim peak template %.0f -> %.0f\n", kMST[is], "B", b.B.first, b.B.second,
          b.mB, a.B.first, a.B.second, a.mB, b.pkOnly, a.pkOnly);
        T("%-9s %-5s | %6.3f +- %-5.3f +- %-5.3f      | %6.3f +- %-5.3f +- %-5.3f      | P_main on T1+T2 %.0f -> %.0f; C(a) %.3f -> %.3f\n", kMST[is], "C", b.C.first, b.C.second,
          b.mC, a.C.first, a.C.second, a.mC, b.simTot, a.simTot, b.Ca.first, a.Ca.first);
        T("%-9s %-5s | %6.3f +- %-5.3f [%.2f-%.2f]  | %6.3f +- %-5.3f [%.2f-%.2f]  | chi2/ndf %.0f/%d -> %.0f/%d\n", kMST[is], "D", b.D.first, b.Dscaled, b.Dlo, b.Dhi,
          a.D.first, a.Dscaled, a.Dlo, a.Dhi, b.Dchi, b.Dndf, a.Dchi, a.Dndf);
        T("%-9s %-5s | %6.3f +- %-5.3f [%.2f-%.2f]  | %6.3f +- %-5.3f [%.2f-%.2f]  | E(a) %.3f -> %.3f; bkg ratio left/right %.3f/%.3f -> %.3f/%.3f\n", kMST[is], "E",
          b.E.first, b.E.second, b.Elo, b.Ehi, a.E.first, a.E.second, a.Elo, a.Ehi, b.Ea.first, a.Ea.first, b.rL, b.rR, a.rL, a.rR);
    }
    double pl[2] = {1e9, 1e9}, ph[2] = {-1e9, -1e9}, ol[2] = {1e9, 1e9}, oh[2] = {-1e9, -1e9};
    for (int m = 0; m < 2; ++m)
        for (int is = 0; is < kNMS; ++is) {
            const MethRes& x = gM[m ? vMain : vOff][is];
            for (double v : {x.B.first, x.C.first, x.Ca.first}) { pl[m] = std::min(pl[m], v); ph[m] = std::max(ph[m], v); }
            for (double v : {x.D.first, x.E.first, x.Ea.first}) { ol[m] = std::min(ol[m], v); oh[m] = std::max(oh[m], v); }
        }
    T("# P_main on both sides (B, C, C(a)): before %.2f-%.2f, after %.2f-%.2f; with the simulated background shape (D, E, E(a)): before %.2f-%.2f, after %.2f-%.2f\n",
      pl[0], ph[0], pl[1], ph[1], ol[0], oh[0], ol[1], oh[1]);
    CloseT();
    // checks
    K("\n## enriched: 'off' against the earlier tables\n");
    // references: tables/eta_summary.txt section 2 (rules v2 since 2026-11-13; rules v1: 0.944, 1.266, 0.893, 1.158)
    Check("E1 bump window data/sim before (eta_summary)", RegRatio(0, 1, 1, 0).first, ses::GetVE("eta_summary.txt", "2.", "E1 bump window", 0, 3).first, 6e-4);
    Check("E1 peak window data/sim before (eta_summary)", RegRatio(0, 1, 3, 0).first, ses::GetVE("eta_summary.txt", "2.", "E1 peak window", 0, 3).first, 6e-4);
    Check("E2 bump window data/sim before (eta_summary)", RegRatio(0, 2, 1, 0).first, ses::GetVE("eta_summary.txt", "2.", "E2 bump window", 0, 3).first, 6e-4);
    Check("E2 peak window data/sim before (eta_summary)", RegRatio(0, 2, 3, 0).first, ses::GetVE("eta_summary.txt", "2.", "E2 peak window", 0, 3).first, 6e-4);
    Check("E1 0.56-0.70 data/sim before (1.361)", RegRatio(0, 1, 5, 0).first, 1.361, 6e-4);
    const std::vector<ses::FitRow> fits = ses::ParseFits();
    for (int is = 1; is < kNES; ++is) {
        const ses::FitRow* f = ses::FindFit(fits, kEST[is], 0, false, false, "M1");
        const ses::FitRow* f3 = ses::FindFit(fits, kEST[is], 0, false, false, "M3");
        const seen::FitRes& r = gFit[vOff][is].m1[0][0][0];
        Check(Form("%s M1 a_bump before (eta_enriched_fit)", kEST[is]), r.a[0], f ? f->Par("a_bump").first : kNaN, 6e-4);
        Check(Form("%s M1 a_peak before", kEST[is]), r.a[1], f ? f->Par("a_peak").first : kNaN, 6e-4);
        Check(Form("%s M1 a_bkg before", kEST[is]), r.a[2], f ? f->Par("a_bkg").first : kNaN, 6e-4);
        Check(Form("%s M1 chi2 before", kEST[is]), r.chi2, f ? f->chi2 : kNaN, 0.051);
        Check(Form("%s M3 a_bkgEta before", kEST[is]), gFit[vOff][is].m3[0].a[2], f3 ? f3->Par("a_bkgEta").first : kNaN, 6e-4);
    }
    const char* mt[4] = {"pre", "pre+OA>2", "E1", "E2"};
    for (int is = 0; is < kNMS; ++is) {
        const MethRes& b = gM[vOff][is];
        Check(Form("%s method A before (eta_summary_peak_methods)", mt[is]), b.A.first, ses::GetVE("eta_summary_peak_methods.txt", "", Form("%s A", mt[is]), -1, 0).first, 6e-4);
        Check(Form("%s method B before", mt[is]), b.B.first, ses::GetVE("eta_summary_peak_methods.txt", "", Form("%s B", mt[is]), -1, 0).first, 6e-4);
        Check(Form("%s method C before", mt[is]), b.C.first, ses::GetVE("eta_summary_peak_methods.txt", "", Form("%s C", mt[is]), -1, 0).first, 6e-4);
        Check(Form("%s method D before", mt[is]), b.D.first, ses::GetVE("eta_summary_peak_methods.txt", "", Form("%s D", mt[is]), -1, 0).first, 6e-4);
        Check(Form("%s method E before", mt[is]), b.E.first, ses::GetVE("eta_summary_peak_methods.txt", "", Form("%s E", mt[is]), -1, 0).first, 6e-4);
    }
}

void BumpTables() {
    OpenT("rw_compare_bump.txt");
    T("# rw_compare_bump: the simulated eta bump (eta -> pi+pi-pi0, e+e- from the pi0: Dalitz decay D or photon conversion C) BEFORE and AFTER the\n");
    T("#   Dalitz-plot weight (KLOE-2 2016, X^2 -> Xmax(Y)^2/3; variant 1 + aY of the X17 slides); HGeant eta -> 3 pi ASSUMED phase space\n");
    T("# simulation: std + lep + eta owned rows, mode A, w (PDG) x physics weight; sigma = sum w / L_sim; data (N - CB) / L_data, L_data = 6080.5 nb^-1 UNCONFIRMED\n\n");
    T("## 1. simulated bump cross section [nb]: M4 < 1.0 | bump window 0.455-0.5225; before -> KLOE (main) / 1 + aY\n");
    T("%-5s | %-30s %-30s %-38s | %-38s\n", "sel", "Dalitz (M4<1)", "conversion (M4<1)", "D+C (M4<1)", "D+C bump window");
    for (int is = 0; is < kNES; ++is) {
        const BumpSim& b = gBS[vOff][is]; const BumpSim& a = gBS[vMain][is]; const BumpSim& s = gBS[vSl][is];
        T("%-5s | %6.3f -> %6.3f / %-6.3f       %6.3f -> %6.3f / %-6.3f       %6.3f +- %5.3f -> %6.3f / %-6.3f (x%.3f) | %6.3f +- %5.3f -> %6.3f / %-6.3f (x%.3f)\n",
          kEST[is], b.all[0].Xsec(), a.all[0].Xsec(), s.all[0].Xsec(), b.all[1].Xsec(), a.all[1].Xsec(), s.all[1].Xsec(), b.all[2].Xsec(), b.all[2].XsecErr(),
          a.all[2].Xsec(), s.all[2].Xsec(), a.all[2].S / b.all[2].S, b.win[2].Xsec(), b.win[2].XsecErr(), a.win[2].Xsec(), s.win[2].Xsec(), a.win[2].S / b.win[2].S);
    }
    T("\n## 2. data / simulation in the bump window and the jump at the bump edge (T1 + T2; (a) in brackets), before -> after\n");
    T("# jump = data/sim in %.4f-%.4f divided by data/sim in %.4f-%.4f (the upper band is mostly multi-pion background)\n", kBand[1][0], kBand[1][1], kBand[0][0], kBand[0][1]);
    T("%-5s | %-30s %-30s | %-30s %-30s | %-30s\n", "sel", "bump window bef", "bump window aft", "band 0.49-0.52 bef -> aft", "band 0.56-0.595 bef -> aft", "jump bef -> aft");
    for (int is = 0; is < kNES; ++is) {
        const VE rb = RegRatio(0, is, 1, 0), ra = RegRatio(1, is, 1, 0);
        const VE lb = BandRatio(vOff, is, 0, 0), la = BandRatio(vMain, is, 0, 0), ub = BandRatio(vOff, is, 1, 0), ua = BandRatio(vMain, is, 1, 0);
        const VE jb = Jump(vOff, is, 0), ja = Jump(vMain, is, 0);
        T("%-5s | %6.3f +- %-5.3f (%5.3f)        %6.3f +- %-5.3f (%5.3f)        | %6.3f -> %6.3f +- %-5.3f       %6.3f -> %6.3f +- %-5.3f       | %6.3f -> %6.3f +- %-5.3f (1+aY %.3f)\n",
          kEST[is], rb.first, rb.second, RegRatio(0, is, 1, 1).first, ra.first, ra.second, RegRatio(1, is, 1, 1).first, lb.first, la.first, la.second, ub.first,
          ua.first, ua.second, jb.first, ja.first, ja.second, gHaveBand[vSl] ? Jump(vSl, is, 0).first : kNaN);
    }
    T("\n## 3. M(pipi) in the bump window (M4 0.455-0.5225): data N - CB (CB per 10 MeV M(pipi) bin) against the simulated T1 + T2\n");
    T("# shape chi2: simulation scaled to the data integral (composition fixed); abs = data / sim (absolute, L_data unconfirmed);\n");
    T("# fit: data = a_bump x bump + a_rest x (T1+T2 - bump), variance data + a^2 var(bump) + c^2 var(rest) (corrected after the review of 2026-11-12:\n");
    T("#   the earlier variance used var(T1+T2) for the rest and counted the bump template twice)\n");
    T("%-5s %-7s | %-16s %-10s | %-10s %-10s %-12s %-18s\n", "sel", "sim", "frac M>0.34", "mean", "chi2/ndf", "", "", "abs data/sim");
    for (int j = 0; j < kNMP; ++j) {
        const MppStat d0 = MppCompare(j, vOff, 0);
        T("%-5s %-7s | %6.3f +- %-6.3f %8.4f   |\n", kEST[kMPSel[j]], "data", d0.fD, d0.fDe, d0.meanD);
        for (int iv : {vOff, vMain, vSl}) {
            const MppStat s = MppCompare(j, iv, 0), sb = MppCompare(j, iv, 1);
            const MppFitRes f = MppFit(j, iv);
            T("%-5s %-7s | %6.3f            %8.4f   | %6.1f/%-4d  bump D+C alone: frac %.3f mean %.4f | %6.3f +- %-6.3f | fit a_bump %.3f +- %.3f, a_rest %.3f +- %.3f"
              " (rho %+.2f), chi2/ndf %.1f/%d\n", kEST[kMPSel[j]], gV[iv].tag.c_str(), s.fS, s.meanS, s.chi2, s.ndf, sb.fS, sb.meanS, s.ratio, s.ratioE,
              f.a, f.ea, f.c, f.ec, f.rho, f.chi2, f.ndf);
        }
        T("  per bin data/sim (before | after):");
        const seen::DH& d = gMppD[j][0];
        for (int i = 0; i < d.nb; ++i) {
            const double x = d.lo + (i + 0.5) * d.bw;
            const double sb = gMppS[vOff][j][0].s[i] * kToRows, sa = gMppS[vMain][j][0].s[i] * kToRows;
            if (sb > 0 && sa > 0 && d.n[i] > 0) T(" %.3f: %.2f|%.2f", x, d.sig(i) / sb, d.sig(i) / sa);
        }
        T("\n");
    }
    T("\n## 4. M(pipi) > 0.34 efficiency of the bump (base pre, all M4; sim tagged = photon 0.119 <= M(ee gamma) < 0.173, emulated ECAL tag, UNVERIFIED)\n");
    T("%-26s %-18s %-18s %-18s | %-16s\n", "sample", "before", "KLOE (main)", "1 + aY", "data");
    const EffOut& e = gE[2];
    T("%-26s %6.3f +- %-8.3f %6.3f +- %-8.3f %6.3f +- %-8.3f | %6.3f +- %-6.3f (ECAL-tagged eta, data study, no vz cut)\n", "bump D+C photon-tagged",
      e.s[vOff][eTag].eps, e.s[vOff][eTag].err, e.s[vMain][eTag].eps, e.s[vMain][eTag].err, e.s[vSl][eTag].eps, e.s[vSl][eTag].err, e.dTag.first, e.dTag.second);
    T("%-26s %6.3f +- %-8.3f %6.3f +- %-8.3f %6.3f +- %-8.3f\n", "bump D+C all", e.s[vOff][eBump].eps, e.s[vOff][eBump].err, e.s[vMain][eBump].eps, e.s[vMain][eBump].err,
      e.s[vSl][eBump].eps, e.s[vSl][eBump].err);
    T("# pulls (tagged): before %+.1f, KLOE %+.1f, 1 + aY %+.1f\n", Pull(e.s[vOff][eTag].eps, e.s[vOff][eTag].err, e.dTag.first, e.dTag.second),
      Pull(e.s[vMain][eTag].eps, e.s[vMain][eTag].err, e.dTag.first, e.dTag.second), Pull(e.s[vSl][eTag].eps, e.s[vSl][eTag].err, e.dTag.first, e.dTag.second));
    T("\n## 5. the other single cuts on the bump (tagged; pT4 and vertex at the bump start base), before -> KLOE / 1 + aY, data tagged\n");
    for (int k = 0; k < kNEC; ++k)
        T("%-26s %6.3f +- %5.3f -> %6.3f / %6.3f | data %s\n", kEC[k].txt, gE[k].s[vOff][eTag].eps, gE[k].s[vOff][eTag].err, gE[k].s[vMain][eTag].eps,
          gE[k].s[vSl][eTag].eps, gE[k].dTag.first == gE[k].dTag.first ? Form("%.3f +- %.3f", gE[k].dTag.first, gE[k].dTag.second) : "-");
    T("\n## 6. what the KLOE weight does and does not fix (review of 2026-11-12)\n");
    T("# fixed: the M(pipi) > 0.34 efficiency of the photon-tagged bump (section 4).\n");
    T("# not fixed: with the composition of the simulation fixed, the M(pipi) shape in the bump window gets WORSE (section 3), the data fraction above\n");
    T("#   0.34 is overshot, and the weight opens a deficit of the simulation below the bump; the corrected bump plus background describe neither\n");
    T("#   M(pipi) nor M4 unless the bump is scaled by a factor that depends on the observable (below). The bump problem is not shown to be\n");
    T("#   normalisation only. Factors and ratios are absolute (L_data UNCONFIRMED) and come from poor fits: comparisons, not yields.\n");
    T("%-5s | %-24s %-24s | %-18s %-18s | %-16s %-16s | %-26s\n", "sel", "Mpipi shape chi2 bef->aft", "Mpipi frac>0.34 d / s", "a_bump (Mpipi fit)", "chi2/ndf",
      "alpha_bump M1 bef", "alpha_bump M1 aft", "data/sim below bump bef->aft");
    for (int j = 0; j < kNMP; ++j) {
        const int is = kMPSel[j];
        const MppStat s0 = MppCompare(j, vOff, 0), s1 = MppCompare(j, vMain, 0);
        const MppFitRes f0 = MppFit(j, vOff), f1 = MppFit(j, vMain);
        const bool hf = gFit[vOff][is].ok && gFit[vMain][is].ok;
        const VE ab = hf ? Scaled(gFit[vOff][is].m1[0][0][0], 0) : VE(kNaN, kNaN), aa = hf ? Scaled(gFit[vMain][is].m1[0][0][0], 0) : VE(kNaN, kNaN);
        const VE bb = RegRatio(0, is, 0, 0), ba = RegRatio(1, is, 0, 0);
        T("%-5s | %6.1f -> %6.1f / %-6d %5.3f / %5.3f -> %5.3f  | %5.3f +- %-5.3f -> %5.3f +- %-5.3f | %4.1f -> %4.1f / %-3d | %5.3f +- %-6.3f %5.3f +- %-6.3f | %5.3f +- %5.3f -> %5.3f +- %5.3f\n",
          kEST[is], s0.chi2, s1.chi2, s1.ndf, s0.fD, s0.fS, s1.fS, f0.a, f0.ea, f1.a, f1.ea, f0.chi2, f1.chi2, f1.ndf, ab.first, ab.second, aa.first, aa.second,
          bb.first, bb.second, ba.first, ba.second);
    }
    T("# alpha_bump M1: template fit 0.40-0.75 (errors x sqrt(chi2/ndf)); pre is the reference selection, not eta-enriched; below bump = M4 0.35-0.455\n");
    CloseT();
    K("\n## bump: checks\n");
    Check("M(pipi) > 0.34 tagged before (0.444)", gE[2].s[vOff][eTag].eps, 0.444, 6e-4);
    // references: tables/rw_weights_xsec.txt (rules v2 since 2026-11-13; rules v1: 13.9862, 10.0757, 5.7325)
    Check("pre bump D+C all M4 before [nb] (rw_weights_xsec)", gBS[vOff][0].all[2].Xsec(), XsecTab("pre", "bump total (2+3)", "all M4", 0), 2e-3);
    Check("pre bump D+C bump window after [nb] (rw_weights_xsec)", gBS[vMain][0].win[2].Xsec(), XsecTab("pre", "bump total (2+3)", "bump win", 1), 6e-4);
    Check("E1 bump D+C bump window after [nb] (rw_weights_xsec)", gBS[vMain][1].win[2].Xsec(), XsecTab("E1", "bump total (2+3)", "bump win", 1), 6e-4);
}

// ---- figures of parts 2 and 3
TH1D* StackSum(const std::vector<WSpec>& w, int upTo, const char* name) {
    TH1D* h = w[seen::kStack[0]].Hist(name, 2);
    for (int k = 1; k <= upTo; ++k) { TH1D* x = w[seen::kStack[k]].Hist(Nm("x").c_str(), 2); h->Add(x); delete x; }
    return h;
}
void DrawM4Column(TCanvas* cv, int is, double x1, double x2, double yTop, bool withLin) {
    // log panel 0.35-1.00 with ratio pad, optional linear panel 0.40-0.70 with ratio pad
    const double yLogBot = withLin ? yTop - 0.46 : 0.0;
    const double yLogR = yLogBot + (yTop - yLogBot) * 0.30;
    cv->cd();
    seen::MkPad(Nm("p").c_str(), x1, yLogR, x2, yTop, true, 0.15, 0.02, 0.07, 0.02);
    TH1D* hN = dLog[is].Hist("us", Nm("n").c_str(), 2); TH1D* hC = dLog[is].Hist("cb", Nm("c").c_str(), 2); TH1D* hS = dLog[is].Hist("sig", Nm("s").c_str(), 2);
    TH1D* sb = eLog[0][is][seen::cSum].Hist(Nm("sb").c_str(), 2); TH1D* sa = eLog[1][is][seen::cSum].Hist(Nm("sa").c_str(), 2);
    TH1D* saa = eLog[1][is][seen::cSumA].Hist(Nm("saa").c_str(), 2);
    const double ymax = seen::HMax(hN) * 60, ymin = std::max(1e-1, seen::HMinPos(hS) * 0.3);
    TH1D* fr = ses::Frame(Nm("f").c_str(), 100, 0.35, 1.00, "", seen::kYM4, ymin, ymax, 0.06, 1.15);
    fr->GetXaxis()->SetLabelSize(0); fr->Draw("axis");
    const int lc[6] = {kColPhysBkg, kColEtaBkg, kColBumpDalitz, kColBumpConv, kColPeakConv, kColPeakDirect};
    for (int k = 0; k < 6; ++k) {
        TH1D* h = eLog[1][is][seen::kStack[k]].Hist(Nm("k").c_str(), 2);
        seen::StyleH(h, lc[k], 2, 1); h->Draw("HIST same");
    }
    seen::StyleH(sb, kColBefore, 3, 2); sb->Draw("HIST same");
    seen::StyleH(saa, kColAfter, 2, 3); saa->Draw("HIST same");
    seen::StyleH(sa, kColAfter, 3, 1); sa->Draw("HIST same");
    seen::StyleM(hN, kColData, 20, 0.7); hN->Draw("E1 same");
    seen::StyleM(hC, kColCB, 21, 0.6); hC->Draw("E1 same");
    seen::StyleM(hS, kColSig, 20, 0.7); hS->Draw("E1 same");
    seen::VLines({scan::kBumpLo, scan::kBumpHi, scan::kPeakLo, scan::kPeakHi, 0.56}, ymin, ymax);
    seen::Txt(0.17, 0.95, Form("#bf{%s}", ELatex(is).c_str()), 0.040);
    const VE rb = RegRatio(0, is, 1, 0), ra = RegRatio(1, is, 1, 0), pb = RegRatio(0, is, 3, 0), pa = RegRatio(1, is, 3, 0);
    const double tsz = withLin ? 0.040 : 0.036;
    seen::Txt(0.17, 0.87, Form("data/sim, bump window: %.3f #rightarrow %.3f", rb.first, ra.first), tsz);
    seen::Txt(0.17, 0.82, Form("data/sim, peak window: %.3f #rightarrow %.3f", pb.first, pa.first), tsz);
    seen::Txt(0.17, 0.77, Form("jump (0.56-0.595 / 0.49-0.52): %.3f #rightarrow %.3f", Jump(vOff, is, 0).first, Jump(vMain, is, 0).first), tsz);
    cv->cd();
    seen::MkPad(Nm("p").c_str(), x1, yLogBot, x2, yLogR, false, 0.15, 0.02, 0.02, 0.30);
    TH1D* fr2 = ses::Frame(Nm("f").c_str(), 100, 0.35, 1.00, seen::kXM4, "data/sim", 0.4, 1.8, 0.12, 0.55);
    fr2->GetYaxis()->SetNdivisions(505); fr2->Draw("axis");
    ses::HLine(0.35, 1.00, 1., kBlack, 2, 1);
    TH1D* rbh = seen::RatioH(hS, sb, Nm("rb").c_str()); seen::StyleM(rbh, kColBefore, 24, 0.7); rbh->Draw("E1 same");
    TH1D* rah = seen::RatioH(hS, sa, Nm("ra").c_str()); seen::StyleM(rah, kColSig, 20, 0.7); rah->Draw("E1 same");
    seen::VLines({scan::kBumpLo, scan::kBumpHi, scan::kPeakLo, scan::kPeakHi, 0.56}, 0.4, 1.8);
    if (!withLin) return;
    const double yLinTop = yLogBot - 0.02, yLinR = 0.0 + (yLinTop - 0.0) * 0.30;
    cv->cd();
    seen::MkPad(Nm("p").c_str(), x1, yLinR, x2, yLinTop, false, 0.15, 0.02, 0.05, 0.02);
    TH1D* lS = dLin[is].Hist("sig", Nm("ls").c_str(), 2);
    TH1D* lsb = eLin[0][is][seen::cSum].Hist(Nm("lsb").c_str(), 2);
    TH1D* stk[6];
    for (int k = 0; k < 6; ++k) stk[k] = StackSum(eLin[1][is], k, Nm("st").c_str());
    const double ym = std::max(seen::HMax(lS), seen::HMax(lsb)) * 1.15;
    TH1D* fr3 = ses::Frame(Nm("f").c_str(), 100, 0.40, 0.70, "", seen::kYM4, 0, ym, 0.06, 1.15);
    fr3->GetXaxis()->SetLabelSize(0); fr3->Draw("axis");
    const int fc[6] = {kColPhysFill, kColEtaFill, kColBumpDalitz, kColBumpConv, kColPeakConv, kColPeakDirect};
    for (int k = 5; k >= 0; --k) { seen::StyleH(stk[k], lc[k], 1, 1, fc[k]); stk[k]->Draw("HIST same"); }
    seen::StyleH(lsb, kColBefore, 3, 2); lsb->Draw("HIST same");
    seen::StyleM(lS, kColSig, 20, 0.8); lS->Draw("E1 same");
    seen::VLines({scan::kBumpLo, scan::kBumpHi, scan::kPeakLo, scan::kPeakHi, 0.56}, 0, ym);
    gPad->RedrawAxis();
    seen::Txt(0.17, 0.88, "stack: after the weights; grey dashed: sum before", 0.045);
    cv->cd();
    seen::MkPad(Nm("p").c_str(), x1, 0.0, x2, yLinR, false, 0.15, 0.02, 0.02, 0.30);
    TH1D* fr4 = ses::Frame(Nm("f").c_str(), 100, 0.40, 0.70, seen::kXM4, "data/sim", 0.4, 1.8, 0.12, 0.55);
    fr4->GetYaxis()->SetNdivisions(505); fr4->Draw("axis");
    ses::HLine(0.40, 0.70, 1., kBlack, 2, 1);
    TH1D* lra = seen::RatioH(lS, stk[5], Nm("lra").c_str()); seen::StyleM(lra, kColSig, 20, 0.7);
    TH1D* lrb = seen::RatioH(lS, lsb, Nm("lrb").c_str()); seen::StyleM(lrb, kColBefore, 24, 0.7);
    lrb->Draw("E1 same"); lra->Draw("E1 same");
    seen::VLines({scan::kBumpLo, scan::kBumpHi, scan::kPeakLo, scan::kPeakHi, 0.56}, 0.4, 1.8);
}
void M4Legend(double x1, double y1, double x2, double y2, int ncol, double ts) {
    TLegend* lg = new TLegend(x1, y1, x2, y2);
    lg->SetNColumns(ncol); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(ts);
    lg->AddEntry(seen::Dummy(Nm("l").c_str(), kColData, 2, 1, 0, 20, kColData), "data N", "p");
    lg->AddEntry(seen::Dummy(Nm("l").c_str(), kColCB, 2, 1, 0, 21, kColCB), "data CB", "p");
    lg->AddEntry(seen::Dummy(Nm("l").c_str(), kColSig, 2, 1, 0, 20, kColSig), "data N - CB (ratio: / sim after)", "p");
    lg->AddEntry(seen::Dummy(Nm("l").c_str(), kColBefore, 2, 1, 0, 24, kColBefore), "ratio: data / sim before", "p");
    lg->AddEntry(seen::Dummy(Nm("l").c_str(), kColAfter, 3, 1, 0), "sim T1 + T2 after the weights", "l");
    lg->AddEntry(seen::Dummy(Nm("l").c_str(), kColBefore, 3, 2, 0), "sim T1 + T2 before", "l");
    lg->AddEntry(seen::Dummy(Nm("l").c_str(), kColAfter, 2, 3, 0), "after, variant (a) + ghost pairs", "l");
    const int lc[6] = {kColPhysBkg, kColEtaBkg, kColBumpDalitz, kColBumpConv, kColPeakConv, kColPeakDirect};
    const int fc[6] = {kColPhysFill, kColEtaFill, kColBumpDalitz, kColBumpConv, kColPeakConv, kColPeakDirect};
    const char* ln[6] = {"bkg multi-pion", "bkg eta-related", "bump #pi^{0} Dalitz (KLOE)", "bump conversion (KLOE)", "peak #pi^{+}#pi^{-}#gamma conv. (F_{V})",
                         "peak #pi^{+}#pi^{-}e^{+}e^{-} direct (V3)"};
    for (int k = 5; k >= 0; --k) lg->AddEntry(seen::Dummy(Nm("l").c_str(), lc[k], 2, 1, fc[k]), ln[k], "f");
    lg->Draw();
}
void FigM4() {
    {
        TCanvas* cv = new TCanvas(Nm("cM12").c_str(), "", 2800, 2600);
        cv->cd();
        seen::MkPad(Nm("h").c_str(), 0, 0.90, 1, 1, false, 0, 0, 0, 0);
        seen::Txt(0.01, 0.82, "Data against the physics-corrected simulation in E1 (slides) and E2 (slides without OA > 2#circ), absolute scale; top: log 0.35-1.00 GeV/c^{2} (simulated components as lines, after), bottom: stacked linear zoom", 0.115);
        seen::Txt(0.01, 0.66, (SimNote() + "; " + DNote()).c_str(), 0.085);
        seen::Txt(0.01, 0.53, WNote().c_str(), 0.085);
        M4Legend(0.0, 0.0, 1.0, 0.50, 6, 0.085);
        DrawM4Column(cv, 1, 0.0, 0.5, 0.90, true);
        DrawM4Column(cv, 2, 0.5, 1.0, 0.90, true);
        Save(cv, "rw_compare_m4_E12", true);
    }
    {
        TCanvas* cv = new TCanvas(Nm("cM34").c_str(), "", 3600, 1500);
        cv->cd();
        seen::MkPad(Nm("h").c_str(), 0, 0.84, 1, 1, false, 0, 0, 0, 0);
        seen::Txt(0.01, 0.84, "Data against the physics-corrected simulation in E3 (old code) and E4 (#eta band), each without and with OA > 2#circ, absolute scale, M_{4} 0.35-1.00 GeV/c^{2}", 0.09);
        seen::Txt(0.01, 0.70, (SimNote() + "; " + DNote()).c_str(), 0.065);
        M4Legend(0.0, 0.0, 1.0, 0.64, 7, 0.065);
        for (int c = 0; c < 4; ++c) DrawM4Column(cv, 3 + c, c * 0.25, (c + 1) * 0.25, 0.84, false);
        Save(cv, "rw_compare_m4_E34", true);
    }
}
void FigRegions() {
    TCanvas* cv = new TCanvas(Nm("cR").c_str(), "", 3000, 1700);
    cv->cd();
    seen::MkPad(Nm("h").c_str(), 0, 0.90, 1, 1, false, 0, 0, 0, 0);
    seen::Txt(0.01, 0.68, "Data (N - CB) / simulation per M_{4} region in the #eta-enriched selections, before (grey open) and after (blue full) the physics weights; teal: after, variant (a)", 0.20);
    seen::Txt(0.01, 0.30, "regions above 0.56 GeV/c^{2} are 80-90 % multi-pion background (not judged); last point: jump = data/sim(0.56-0.595) / data/sim(0.4875-0.5225); L_{data} = 6080.5 nb^{-1} UNCONFIRMED", 0.15);
    const std::vector<std::string> lab = {"<0.455", "bump", "gap", "peak", "0.554-0.56", "0.56-0.70", "0.70-1.0", "all", "jump"};
    for (int p = 0; p < 6; ++p) {
        const int is = p + 1;
        cv->cd();
        const double x1 = (p % 3) / 3., y1 = p < 3 ? 0.45 : 0.0;
        seen::MkPad(Nm("p").c_str(), x1, y1, x1 + 1. / 3, y1 + 0.45, false, 0.12, 0.02, 0.08, 0.10);
        TH1D* fr = ses::CatFrame(Nm("f").c_str(), lab, "data / simulation", 0.5, 1.6, 0.05, 0.045, 1.1);
        fr->Draw("axis");
        ses::HLine(-0.5, 8.5, 1., kBlack, 2, 1);
        std::vector<double> xb, xa, xaa;
        std::vector<VE> b, a, aa;
        for (int k = 0; k <= seen::kNReg; ++k) {
            xb.push_back(k - 0.2); xa.push_back(k + 0.05); xaa.push_back(k + 0.25);
            b.push_back(RegRatio(0, is, k, 0)); a.push_back(RegRatio(1, is, k, 0)); aa.push_back(RegRatio(1, is, k, 1));
        }
        xb.push_back(8 - 0.2); xa.push_back(8 + 0.05); xaa.push_back(8 + 0.25);
        b.push_back(Jump(vOff, is, 0)); a.push_back(Jump(vMain, is, 0)); aa.push_back(Jump(vMain, is, 1));
        ses::Pts(xb, b, kColBefore, 24, 1.8)->Draw("P same");
        ses::Pts(xaa, aa, kColAfter, 22, 1.6)->Draw("P same");
        ses::Pts(xa, a, kColSig, 20, 1.9)->Draw("P same");
        seen::Txt(0.14, 0.93, Form("#bf{%s}", seen::SelLatex(is - 1).c_str()), 0.040);
    }
    Save(cv, "rw_compare_regions", false);
}
void FigFits() {
    TCanvas* cv = new TCanvas(Nm("cF").c_str(), "", 3000, 1800);
    cv->cd();
    seen::MkPad(Nm("h").c_str(), 0, 0.90, 1, 1, false, 0, 0, 0, 0);
    seen::Txt(0.01, 0.70, "Template fits of data N - CB (M_{4} 0.40-0.75) with the simulated templates before (grey open) and after (teal full) the physics weights; errors #times #sqrt{#chi^{2}/ndf}", 0.19);
    seen::Txt(0.01, 0.40, "M1: bump, peak (direct + conversion, simulated ratio), physical background; M3: background split into #eta-related and multi-pion; dashed bars: range over 2 fit ranges #times 2 mass scales #times 2 ghost variants (after)", 0.14);
    seen::Txt(0.01, 0.12, "#alpha = 1: the simulation at L_{data} = 6080.5 nb^{-1} (UNCONFIRMED); the factors compare shapes in poor fits, they are NOT yields", 0.14);
    std::vector<std::string> lab(kEST, kEST + kNES);
    const char* yt[6] = {"M1 #alpha_{bump}", "M1 #alpha_{peak}", "M1 #alpha_{bkg}", "M1 #chi^{2}/ndf", "M3 #alpha_{bump}", "M3 #alpha_{bkg}: #eta-related (circles), multi-pion (squares)"};
    const double ylo[6] = {0, 0, 0.8, 0, -0.3, 0}, yhi[6] = {1.4, 2.2, 1.6, 7.0, 1.4, 3.0};
    for (int p = 0; p < 6; ++p) {
        cv->cd();
        const double x1 = (p % 3) / 3., y1 = p < 3 ? 0.45 : 0.0;
        seen::MkPad(Nm("p").c_str(), x1, y1, x1 + 1. / 3, y1 + 0.45, false, 0.14, 0.02, 0.06, 0.10);
        TH1D* fr = ses::CatFrame(Nm("f").c_str(), lab, yt[p], ylo[p], yhi[p], 0.05, 0.055, 1.25);
        fr->Draw("axis");
        if (p != 3) ses::HLine(-0.5, kNES - 0.5, 1., kBlack, 2, 1);
        for (int m = 0; m < 2; ++m) {
            const int iv = m ? vMain : vOff;
            std::vector<double> x, x2; std::vector<VE> y, y2;
            for (int is = 0; is < kNES; ++is) {
                const FitSet& F = gFit[iv][is];
                const seen::FitRes& r = F.m1[0][0][0]; const seen::FitRes& r3 = F.m3[0];
                x.push_back(is + (m ? 0.15 : -0.15));
                if (p < 3) y.push_back(Scaled(r, p));
                else if (p == 3) y.push_back(VE(r.chi2 / r.ndf, 0));
                else if (p == 4) y.push_back(Scaled(r3, 0));
                else { y.push_back(Scaled(r3, 2)); x2.push_back(is + (m ? 0.25 : -0.05)); y2.push_back(Scaled(r3, 3)); }
                if (m && p < 3) {
                    double lo, hi; M1Range(F, p, lo, hi);
                    TLine* l = new TLine(is + 0.25, lo, is + 0.25, hi); l->SetLineColor(kColAfter); l->SetLineStyle(2); l->SetLineWidth(3); l->Draw();
                }
            }
            ses::Pts(x, y, m ? kColAfter : kColBefore, m ? 20 : 24, 1.9)->Draw("P same");
            if (p == 5) ses::Pts(x2, y2, m ? kColAfter : kColBefore, m ? 21 : 25, 1.9)->Draw("P same");
        }
        if (p == 3) for (int is = 0; is < kNES; ++is) {
            const seen::FitRes& rb = gFit[vOff][is].m3[0]; const seen::FitRes& ra = gFit[vMain][is].m3[0];
            std::vector<double> x = {is - 0.15 + 0.08, is + 0.15 + 0.08}; std::vector<VE> y = {VE(rb.chi2 / rb.ndf, 0), VE(ra.chi2 / ra.ndf, 0)};
            ses::Pts(x, y, kColEtaBkg, 27, 1.6)->Draw("P same");
        }
        if (p == 3) seen::Txt(0.17, 0.86, "violet diamonds: M3 (before, after)", 0.045);
    }
    Save(cv, "rw_compare_fits", true);
}
void FigMethods() {
    TCanvas* cv = new TCanvas(Nm("cMe").c_str(), "", 3000, 1500);
    cv->cd();
    seen::MkPad(Nm("h").c_str(), 0, 0.86, 1, 1, false, 0, 0, 0, 0);
    seen::Txt(0.01, 0.72, "#eta peak: data / simulation by method (reconcile.md), before (grey open) and after (full) the physics weights", 0.20);
    seen::Txt(0.01, 0.42, "A: P_{main} / sim content all M_{4}; Aw: in the window; B: P_{main} on data and on the sim peak template; C: P_{main} on data and on the sim T1 + T2; D: template fit M1 #alpha_{peak} "
                          "(dashed: range over ranges, scales, ghost); E: side bands (dashed: left or right band alone)", 0.12);
    seen::Txt(0.01, 0.14, "thick bars: stat #oplus fit-model spread (48 variants); thin: also the absolute scale incl. L_{data} (#pm20 %, L_{data} = 6080.5 nb^{-1} UNCONFIRMED)", 0.12);
    const std::vector<std::string> lab = {"A", "Aw", "B", "C", "D", "E"};
    const int col[kNMS] = {kBlack, kGray + 2, kRed + 1, kBlue + 1};
    for (int is = 0; is < kNMS; ++is) {
        cv->cd();
        seen::MkPad(Nm("p").c_str(), is * 0.25, 0, (is + 1) * 0.25, 0.86, false, 0.15, 0.02, 0.06, 0.08);
        TH1D* fr = ses::CatFrame(Nm("f").c_str(), lab, "#eta peak: data / simulation", 0.0, 2.0, 0.05, 0.065, 1.3);
        fr->Draw("axis");
        ses::HLine(-0.5, 5.5, 1., kBlack, 2, 1);
        for (int m = 0; m < 2; ++m) {
            const MethRes& r = gM[m ? vMain : vOff][is];
            const double off = m ? 0.15 : -0.15;
            const VE val[6] = {r.A, r.Aw, r.B, r.C, VE(r.D.first, r.Dscaled), r.E};
            const double mod[6] = {r.mA, r.mAw, r.mB, r.mC, 0, 0};
            for (int k = 0; k < 6; ++k) {
                const double inner = std::hypot(val[k].second, mod[k]), outer = std::sqrt(inner * inner + std::pow(0.2 * val[k].first, 2));
                TGraphErrors* go = new TGraphErrors(1); go->SetPoint(0, k + off, val[k].first); go->SetPointError(0, 0, outer);
                go->SetLineColor(m ? col[is] : kColBefore); go->SetLineWidth(1); go->SetMarkerSize(0); go->Draw("P same");
                TGraphErrors* gi = new TGraphErrors(1); gi->SetPoint(0, k + off, val[k].first); gi->SetPointError(0, 0, inner);
                gi->SetLineColor(m ? col[is] : kColBefore); gi->SetLineWidth(4); gi->SetMarkerColor(m ? col[is] : kColBefore);
                gi->SetMarkerStyle(m ? 20 : 24); gi->SetMarkerSize(2.0); gi->Draw("P same");
            }
            for (int k : {4, 5}) {
                const double lo = k == 4 ? r.Dlo : r.Elo, hi = k == 4 ? r.Dhi : r.Ehi;
                TLine* l = new TLine(k + off + 0.05, lo, k + off + 0.05, hi); l->SetLineColor(m ? col[is] : kColBefore); l->SetLineStyle(2); l->SetLineWidth(3); l->Draw();
            }
        }
        seen::Txt(0.18, 0.95, Form("#bf{%s}", kMST[is]), 0.06);
    }
    Save(cv, "rw_compare_methods", true);
}
void FigBump() {
    TCanvas* cv = new TCanvas(Nm("cB").c_str(), "", 3000, 1900);
    cv->cd();
    seen::MkPad(Nm("h").c_str(), 0, 0.90, 1, 1, false, 0, 0, 0, 0);
    seen::Txt(0.01, 0.70, "The #eta bump (#eta#rightarrow#pi^{+}#pi^{-}#pi^{0}, e^{+}e^{-} from the #pi^{0}) before and after the Dalitz-plot weight (KLOE-2 2016; dotted: 1 + aY of the X17 slides)", 0.20);
    seen::Txt(0.01, 0.40, (SimNote() + "; HGeant #eta#rightarrow3#pi ASSUMED phase space").c_str(), 0.13);
    seen::Txt(0.01, 0.12, (DNote() + "; bottom: M_{#pi#pi} of the rows in the bump window M_{4} 0.455-0.5225 GeV/c^{2}").c_str(), 0.13);
    const int js[2] = {0, 1};   // E1, E2 in kMPSel
    for (int c = 0; c < 2; ++c) {
        const int is = kMPSel[js[c]];
        cv->cd();
        seen::MkPad(Nm("p").c_str(), c / 3., 0.46, (c + 1) / 3., 0.90, false, 0.15, 0.02, 0.08, 0.12);
        TH1D* hb = gBumpM4[vOff][is][2].Hist(Nm("b").c_str(), 2); TH1D* ha = gBumpM4[vMain][is][2].Hist(Nm("a").c_str(), 2);
        TH1D* hs = gBumpM4[vSl][is][2].Hist(Nm("s").c_str(), 2);
        TH1D* hdb = gBumpM4[vOff][is][0].Hist(Nm("db").c_str(), 2); TH1D* hda = gBumpM4[vMain][is][0].Hist(Nm("da").c_str(), 2);
        TH1D* hcb = gBumpM4[vOff][is][1].Hist(Nm("cb").c_str(), 2); TH1D* hca = gBumpM4[vMain][is][1].Hist(Nm("ca").c_str(), 2);
        TH1D* fr = ses::Frame(Nm("f").c_str(), 100, 0.40, 0.60, seen::kXM4, seen::kYM4, 0, seen::HMax(hb) * 1.2, 0.05, 1.4);
        fr->Draw("axis");
        seen::StyleH(hb, kColBefore, 3, 2); hb->Draw("HIST same");
        seen::StyleH(hs, kOrange + 3, 3, 3); hs->Draw("HIST same");
        seen::StyleH(ha, kBlack, 3, 1); ha->Draw("HIST same");
        seen::StyleH(hdb, kColBumpDalitz, 2, 2); hdb->Draw("HIST same"); seen::StyleH(hda, kColBumpDalitz, 2, 1); hda->Draw("HIST same");
        seen::StyleH(hcb, kColBumpConv, 2, 2); hcb->Draw("HIST same"); seen::StyleH(hca, kColBumpConv, 2, 1); hca->Draw("HIST same");
        seen::VLines({scan::kBumpLo, scan::kBumpHi}, 0, seen::HMax(hb) * 1.2);
        seen::Txt(0.18, 0.93, Form("#bf{%s}: sim bump D+C (black), Dalitz (orange), conv. (magenta)", kEST[is]), 0.036);
        seen::Txt(0.18, 0.86, Form("bump window: %.2f #rightarrow %.2f nb (#times%.3f); 1+aY %.2f nb", gBS[vOff][is].win[2].Xsec(), gBS[vMain][is].win[2].Xsec(),
                                   gBS[vMain][is].win[2].S / gBS[vOff][is].win[2].S, gBS[vSl][is].win[2].Xsec()), 0.036);
        seen::Txt(0.18, 0.80, "dashed: before; solid: after (KLOE); dotted: 1 + aY", 0.036);
        // M(pipi) panel
        cv->cd();
        seen::MkPad(Nm("p").c_str(), c / 3., 0.0, (c + 1) / 3., 0.44, false, 0.15, 0.02, 0.08, 0.13);
        const seen::DH& d = gMppD[js[c]][0];
        const double sc = 1. / (kLdata * d.bw);
        TH1D* hd = d.Hist(Nm("d").c_str(), 2, sc);
        TH1D* sb = gMppS[vOff][js[c]][0].Hist(Nm("sb").c_str(), 2); TH1D* sa = gMppS[vMain][js[c]][0].Hist(Nm("sa").c_str(), 2);
        TH1D* bb = gMppS[vOff][js[c]][1].Hist(Nm("bb").c_str(), 2); TH1D* ba = gMppS[vMain][js[c]][1].Hist(Nm("ba").c_str(), 2);
        TH1D* fr2 = ses::Frame(Nm("f").c_str(), 100, 0.27, 0.47, "M_{#pi^{+}#pi^{-}} [GeV/c^{2}]", "d#sigma/dM_{#pi#pi} [nb/(GeV/c^{2})]", 0,
                               std::max(seen::HMax(hd), seen::HMax(sb)) * 1.65, 0.05, 1.4);
        fr2->Draw("axis");
        seen::StyleH(sb, kColBefore, 3, 2); sb->Draw("HIST same"); seen::StyleH(sa, kColAfter, 3, 1); sa->Draw("HIST same");
        seen::StyleH(bb, kColBumpDalitz, 2, 2); bb->Draw("HIST same"); seen::StyleH(ba, kColBumpDalitz, 2, 1); ba->Draw("HIST same");
        seen::StyleM(hd, kColSig, 20, 1.1); hd->Draw("E1 same");
        const MppStat s0 = MppCompare(js[c], vOff, 0), s1 = MppCompare(js[c], vMain, 0);
        const MppFitRes f0 = MppFit(js[c], vOff), f1 = MppFit(js[c], vMain);
        seen::Txt(0.18, 0.93, Form("#bf{%s}, bump window: data N-CB (blue), sim T1+T2 (dashed before, teal after), bump (orange)", kEST[is]), 0.032);
        seen::Txt(0.18, 0.87, Form("fraction M_{#pi#pi} > 0.34: data %.3f #pm %.3f; sim %.3f #rightarrow %.3f", s0.fD, s0.fDe, s0.fS, s1.fS), 0.036);
        seen::Txt(0.18, 0.82, Form("shape #chi^{2}/ndf %.0f/%d #rightarrow %.0f/%d", s0.chi2, s0.ndf, s1.chi2, s1.ndf), 0.036);
        seen::Txt(0.18, 0.77, Form("fit a_{bump}B + a_{rest}R: a_{bump} %.2f#pm%.2f #rightarrow %.2f#pm%.2f, a_{rest} %.2f #rightarrow %.2f", f0.a, f0.ea, f1.a, f1.ea, f0.c, f1.c), 0.032);
        seen::Txt(0.18, 0.72, Form("   #chi^{2}/ndf %.0f/%d #rightarrow %.0f/%d%s", f0.chi2, f0.ndf, f1.chi2, f1.ndf,
                                   gFit[vMain][is].ok ? Form("; M1 #alpha_{bump} %.2f #rightarrow %.2f", gFit[vOff][is].m1[0][0][0].a[0], gFit[vMain][is].m1[0][0][0].a[0]) : ""), 0.032);
        seen::Txt(0.18, 0.67, "fixed composition: shape worse; bump factor differs per observable", 0.030);
    }
    // efficiency panel
    cv->cd();
    seen::MkPad(Nm("p").c_str(), 2 / 3., 0.0, 1.0, 0.90, false, 0.15, 0.02, 0.06, 0.10);
    {
        const std::vector<std::string> lab = {"OA > 2", "OA < 4", "M_{#pi#pi} > 0.34", "p_{T,4} > 400", "#alpha_{lab} < 50", "#pi vertex"};
        const int ks[6] = {0, 1, 2, 4, 5, 6};
        TH1D* fr = ses::CatFrame(Nm("f").c_str(), lab, "efficiency of the photon-tagged bump", 0.2, 1.12, 0.045, 0.040, 1.5);
        fr->Draw("axis");
        std::vector<double> xd, xb, xa, xs; std::vector<VE> d, b, a, s;
        for (int q = 0; q < 6; ++q) {
            const EffOut& e = gE[ks[q]];
            xd.push_back(q); xb.push_back(q - 0.2); xa.push_back(q + 0.15); xs.push_back(q + 0.3);
            d.push_back(e.dTag); b.push_back(VE(e.s[vOff][eTag].eps, e.s[vOff][eTag].err)); a.push_back(VE(e.s[vMain][eTag].eps, e.s[vMain][eTag].err));
            s.push_back(VE(e.s[vSl][eTag].eps, e.s[vSl][eTag].err));
        }
        TGraphErrors* gb = ses::Pts(xb, b, kColBefore, 25, 2.0); gb->Draw("P same");
        TGraphErrors* gs = ses::Pts(xs, s, kOrange + 3, 26, 1.8); gs->Draw("P same");
        TGraphErrors* ga = ses::Pts(xa, a, kColBumpDalitz, 21, 2.0); ga->Draw("P same");
        TGraphErrors* gd = ses::Pts(xd, d, kColData, 20, 2.2, 3); gd->Draw("P same");
        TLegend* lg = new TLegend(0.40, 0.11, 0.98, 0.30); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.030);
        lg->AddEntry(gd, "data: ECAL-tagged #eta (data study)", "pe");
        lg->AddEntry(gb, "sim tagged bump, before", "pe");
        lg->AddEntry(ga, "sim tagged bump, KLOE weight", "pe");
        lg->AddEntry(gs, "sim tagged bump, 1 + aY", "pe");
        lg->Draw();
        seen::Txt(0.17, 0.96, "#bf{bump: single-cut efficiencies (base pre; p_{T,4}, vertex: bump start)}", 0.028);
        const EffOut& e = gE[2];
        seen::Txt(0.17, 0.915, Form("M_{#pi#pi} > 0.34: data %.3f #pm %.3f; sim %.3f #pm %.3f #rightarrow %.3f #pm %.3f (1+aY %.3f)", e.dTag.first, e.dTag.second,
                                   e.s[vOff][eTag].eps, e.s[vOff][eTag].err, e.s[vMain][eTag].eps, e.s[vMain][eTag].err, e.s[vSl][eTag].eps), 0.023);
        seen::Txt(0.17, 0.885, "sim tag: photon with 0.119 #leq M_{ee#gamma} < 0.173 (emulated, UNVERIFIED)", 0.026);
    }
    Save(cv, "rw_compare_bump", true);
}

// =====================================================================================================
// Part 4: variants
// =====================================================================================================
struct VarHead {
    double pkPre = 0, pkPreE = 0, dirPre = 0, dirShare = 0;   // rows all M4 (pre)
    double effOA2 = 0, effOA2E = 0, effOA4 = 0;
    sep::MixRes mOA, mMEE, mOAE2;
    double hiMee = 0;            // share of the simulated peak (pre) at M(ee) > 20 MeV
    double tag034 = 0, tag034E = 0;
    VE bumpWin[2], peakWin[2];   // E1, E2 region ratios
    VE jump[2];
    seen::FitRes m1[2];
    MethRes meth[2];             // E1, E2 (methods A, B, C, D, E)
};
VarHead gVH[kNV];
void Variants() {
    TStopwatch sw;
    for (int iv = 0; iv < kNV; ++iv) {
        VarHead& h = gVH[iv];
        if (!gHavePS[iv]) PeakSimFill(iv);
        h.pkPre = gPS[iv][0].all[2].S * kToRows; h.pkPreE = std::sqrt(gPS[iv][0].all[2].V) * kToRows;
        h.dirPre = gPS[iv][0].all[0].S * kToRows; h.dirShare = gPS[iv][0].all[0].S / gPS[iv][0].all[2].S;
        h.effOA2 = gE[0].s[iv][eTot].eps; h.effOA2E = gE[0].s[iv][eTot].err; h.effOA4 = gE[1].s[iv][eTot].eps;
        h.mOA = Mix(iv, 0, 0); h.mMEE = Mix(iv, 0, 1); h.mOAE2 = Mix(iv, 4, 0);
        { double hh = 0, tt = 0; for (int j = 0; j < 5; ++j) { const double y = gPS[iv][0].c[1][j][0] + gPS[iv][0].c[1][j][1]; tt += y; if (j >= 2) hh += y; } h.hiMee = hh / tt; }
        h.tag034 = gE[2].s[iv][eTag].eps; h.tag034E = gE[2].s[iv][eTag].err;
        if (!gHaveBand[iv]) BandSums(iv);
        for (int j = 0; j < 2; ++j) {
            const int is = 1 + j;   // E1, E2
            const SimOpt o = Opt(iv);
            const std::vector<WSum> wb = SumMany(*gD, SWin(ESim(is), scan::kBumpLo, scan::kBumpHi), {CompSim()}, o);
            const std::vector<WSum> wp = SumMany(*gD, SWin(ESim(is), scan::kPeakLo, scan::kPeakHi), {CompSim()}, o);
            const scan::WinSum& db = dReg[is].win[2 + 1]; const scan::WinSum& dp = dReg[is].win[2 + 3];
            h.bumpWin[j] = Ratio(db.sig() / kLdata, db.sigErr() / kLdata, wb[0].Xsec(), wb[0].XsecErr());
            h.peakWin[j] = Ratio(dp.sig() / kLdata, dp.sigErr() / kLdata, wp[0].Xsec(), wp[0].XsecErr());
            h.jump[j] = Jump(iv, is, 0);
            if (!gFit[iv][is].ok) Fits(iv, is);
            h.m1[j] = gFit[iv][is].m1[0][0][0];
            if (!gM[iv][2 + j].ok) gM[iv][2 + j] = Methods(2 + j, iv);
            h.meth[j] = gM[iv][2 + j];
        }
        printf("Variants: %s done (%.0f s)\n", gV[iv].tag.c_str(), sw.RealTime()); sw.Continue();
    }
    OpenT("rw_compare_variants.txt");
    T("# rw_compare_variants: headline numbers for the variants of the physics weights (sim_eta_physics.h codes)\n");
    for (int iv = 0; iv < kNV; ++iv) T("#   %-5s code %3d: %s\n", gV[iv].tag.c_str(), gV[iv].code, gV[iv].lab.c_str());
    T("# simulation as in rw_compare_peak.txt / rw_compare_enriched.txt; L_data = 6080.5 nb^-1 UNCONFIRMED; V3 (main) and V4 use the true decay angles\n");
    T("# scale-free rows first (shares, efficiencies, chi2); absolute rows (contents, data/sim, alpha, methods) scale with 6080.5 / L_data and the peak\n");
    T("#   data/sim depends on the method: NOT physics results\n\n");
    T("%-48s", "quantity");
    for (int iv = 0; iv < kNV; ++iv) T(" %13s", gV[iv].tag.c_str());
    T("\n");
    auto row = [&](const char* name, std::function<double(const VarHead&)> f, const char* fmt) {
        T("%-48s", name);
        for (int iv = 0; iv < kNV; ++iv) T(fmt, f(gVH[iv]));
        T("\n");
    };
    row("pre: sim peak total all M4 [rows at L_data]", [](const VarHead& h) { return h.pkPre; }, " %13.0f");
    row("pre: sim direct all M4 [rows]", [](const VarHead& h) { return h.dirPre; }, " %13.0f");
    row("pre: direct share of the sim peak", [](const VarHead& h) { return h.dirShare; }, " %13.3f");
    row("pre: data/sim method A (data 4671 rows)", [](const VarHead& h) { return gPD[0].tot.S / h.pkPre; }, " %13.3f");
    row("pre: eff OA > 2, sim peak total (data 0.778)", [](const VarHead& h) { return h.effOA2; }, " %13.3f");
    row("pre: pull OA > 2 (sim - data) / sigma", [](const VarHead& h) { return Pull(h.effOA2, h.effOA2E, gE[0].d.e, gE[0].d.err); }, " %13.2f");
    row("pre: share of the sim peak at M(ee) > 20 MeV", [](const VarHead& h) { return h.hiMee; }, " %13.3f");
    row("pre: eff OA < 4, sim peak total (data 0.586)", [](const VarHead& h) { return h.effOA4; }, " %13.3f");
    row("pre OA classes: one factor chi2 (ndf 4)", [](const VarHead& h) { return h.mOA.chi21; }, " %13.2f");
    row("pre M(ee) classes: one factor chi2 (ndf 4)", [](const VarHead& h) { return h.mMEE.chi21; }, " %13.2f");
    row("E2 OA classes: one factor chi2 (ndf 4)", [](const VarHead& h) { return h.mOAE2.chi21; }, " %13.2f");
    row("pre OA classes: one factor a1", [](const VarHead& h) { return h.mOA.a1; }, " %13.3f");
    row("pre OA classes: two factors, a_conv", [](const VarHead& h) { return h.mOA.a; }, " %13.3f");
    row("pre OA classes: two factors, a_conv error", [](const VarHead& h) { return h.mOA.ea; }, " %13.3f");
    row("pre OA classes: two factors, a_dir", [](const VarHead& h) { return h.mOA.b; }, " %13.3f");
    row("pre OA classes: two factors, a_dir error", [](const VarHead& h) { return h.mOA.eb; }, " %13.3f");
    row("pre M(ee) classes: two factors, a_conv", [](const VarHead& h) { return h.mMEE.a; }, " %13.3f");
    row("pre M(ee) classes: two factors, a_dir", [](const VarHead& h) { return h.mMEE.b; }, " %13.3f");
    row("bump tagged: eff M(pipi) > 0.34 (data 0.705)", [](const VarHead& h) { return h.tag034; }, " %13.3f");
    for (int j = 0; j < 2; ++j) {
        const char* s = j ? "E2" : "E1";
        T("%-48s", Form("%s: data/sim bump window", s)); for (int iv = 0; iv < kNV; ++iv) T(" %13.3f", gVH[iv].bumpWin[j].first); T("\n");
        T("%-48s", Form("%s: data/sim peak window", s)); for (int iv = 0; iv < kNV; ++iv) T(" %13.3f", gVH[iv].peakWin[j].first); T("\n");
        T("%-48s", Form("%s: jump at the bump edge", s)); for (int iv = 0; iv < kNV; ++iv) T(" %13.3f", gVH[iv].jump[j].first); T("\n");
        T("%-48s", Form("%s: M1 alpha_bump", s)); for (int iv = 0; iv < kNV; ++iv) T(" %13.3f", gVH[iv].m1[j].a[0]); T("\n");
        T("%-48s", Form("%s: M1 alpha_peak", s)); for (int iv = 0; iv < kNV; ++iv) T(" %13.3f", gVH[iv].m1[j].a[1]); T("\n");
        T("%-48s", Form("%s: M1 alpha_bkg", s)); for (int iv = 0; iv < kNV; ++iv) T(" %13.3f", gVH[iv].m1[j].a[2]); T("\n");
        T("%-48s", Form("%s: M1 chi2 (ndf 67)", s)); for (int iv = 0; iv < kNV; ++iv) T(" %13.1f", gVH[iv].m1[j].chi2); T("\n");
        const char* mn[5] = {"A", "B", "C", "D", "E"};
        for (int q = 0; q < 5; ++q) {
            T("%-48s", Form("%s: peak data/sim, method %s", s, mn[q]));
            for (int iv = 0; iv < kNV; ++iv) {
                const MethRes& m = gVH[iv].meth[j];
                const double v = q == 0 ? m.A.first : q == 1 ? m.B.first : q == 2 ? m.C.first : q == 3 ? m.D.first : m.E.first;
                T(" %13.3f", v);
            }
            T("\n");
        }
    }
    T("\n# spread of the variants around the main one V3 (V0, V1, V2, V4, 1+aY, GS): largest |x - main| / main per quantity\n");
    auto spread = [&](const char* name, std::function<double(const VarHead&)> f) {
        double mx = 0; std::string who;
        for (int iv : {vV0, vV1, vV2, vV4, vSl, vGS}) {
            const double r = std::fabs(f(gVH[iv]) - f(gVH[vMain])) / std::fabs(f(gVH[vMain]));
            if (r > mx) { mx = r; who = gV[iv].tag; }
        }
        T("%-48s %6.1f %% (%s)\n", name, 100 * mx, who.c_str());
    };
    spread("pre sim peak total", [](const VarHead& h) { return h.pkPre; });
    spread("pre eff OA > 2", [](const VarHead& h) { return h.effOA2; });
    spread("E1 data/sim peak window", [](const VarHead& h) { return h.peakWin[0].first; });
    spread("E2 data/sim bump window", [](const VarHead& h) { return h.bumpWin[1].first; });
    spread("E2 M1 alpha_bump", [](const VarHead& h) { return h.m1[1].a[0]; });
    spread("E2 M1 alpha_peak", [](const VarHead& h) { return h.m1[1].a[1]; });
    spread("E2 method B", [](const VarHead& h) { return h.meth[1].B.first; });
    spread("E2 method C", [](const VarHead& h) { return h.meth[1].C.first; });
    spread("bump tagged eff M(pipi) > 0.34", [](const VarHead& h) { return h.tag034; });
    CloseT();
    printf("Variants: %.1f s\n", sw.RealTime());
}
void FigVariants() {
    TCanvas* cv = new TCanvas(Nm("cV").c_str(), "", 3000, 1700);
    cv->cd();
    seen::MkPad(Nm("h").c_str(), 0, 0.90, 1, 1, false, 0, 0, 0, 0);
    seen::Txt(0.01, 0.65, "Headline numbers for the variants of the physics weights: direct decay V0-V4 (V3 = main), bump 1 + aY (X17 slides), F_{V} Gounaris-Sakurai; 'no weights' = as produced", 0.20);
    seen::Txt(0.01, 0.25, "black lines / bands: data; L_{data} = 6080.5 nb^{-1} UNCONFIRMED; V3 = V2 #times e^{+}, #pi^{+} decay angles; V4 = V3 #times decay-plane term (checked: below 0.5 %)", 0.15);
    std::vector<std::string> lab;
    for (int iv = 0; iv < kNV; ++iv) lab.push_back(gV[iv].tex);
    struct Pan { const char* yt; double lo, hi; };
    const Pan pans[8] = {{"pre: sim peak total [rows]", 0, 30000}, {"pre: eff OA > 2 (peak total)", 0.5, 0.9}, {"pre: one-factor #chi^{2} (ndf 4)", 0, 16},
                         {"pre OA: a_{conv} (circles), a_{dir} (squares)", -2.5, 2.5}, {"bump tagged: eff M_{#pi#pi} > 0.34", 0.3, 0.85},
                         {"data/sim: bump (circles), peak (squares)", 0.7, 1.4}, {"E2, M1: #alpha_{bump} (circles), #alpha_{peak} (squares)", 0, 1.6},
                         {"E2 peak data/sim: B (circ.), C (sq.), E (diam.)", 0, 1.6}};
    for (int p = 0; p < 8; ++p) {
        cv->cd();
        const double x1 = (p % 4) / 4., y1 = p < 4 ? 0.45 : 0.0;
        seen::MkPad(Nm("p").c_str(), x1, y1, x1 + 0.25, y1 + 0.45, false, 0.16, 0.02, 0.07, 0.12);
        TH1D* fr = ses::CatFrame(Nm("f").c_str(), lab, pans[p].yt, pans[p].lo, pans[p].hi, 0.05, 0.045, 1.4);
        fr->Draw("axis");
        std::vector<double> x0, x1v, x2;
        std::vector<VE> a, b, c;
        for (int iv = 0; iv < kNV; ++iv) {
            const VarHead& h = gVH[iv];
            x0.push_back(iv - 0.15); x1v.push_back(iv + 0.15); x2.push_back(iv);
            switch (p) {
                case 0: a.push_back(VE(h.pkPre, h.pkPreE)); b.push_back(VE(h.dirPre, 0)); break;
                case 1: a.push_back(VE(h.effOA2, h.effOA2E)); break;
                case 2: a.push_back(VE(h.mOA.chi21, 0)); b.push_back(VE(h.mMEE.chi21, 0)); break;
                case 3: a.push_back(VE(h.mOA.a, h.mOA.ea)); b.push_back(VE(h.mOA.b, h.mOA.eb)); break;
                case 4: a.push_back(VE(h.tag034, h.tag034E)); break;
                case 5: a.push_back(h.bumpWin[1]); b.push_back(h.peakWin[1]); c.push_back(h.bumpWin[0]); break;
                case 6: a.push_back(Scaled(h.m1[1], 0)); b.push_back(Scaled(h.m1[1], 1)); break;
                case 7: a.push_back(h.meth[1].B); b.push_back(h.meth[1].C); c.push_back(h.meth[1].E); break;
            }
        }
        if (p == 0) { ses::HLine(-0.5, kNV - 0.5, gPD[0].tot.S, kBlack, 1, 3); seen::Txt(0.2, 0.86, Form("line: data %.0f #pm %.0f; squares: direct", gPD[0].tot.S, gPD[0].tot.eT()), 0.042); }
        if (p == 1) { ses::HLine(-0.5, kNV - 0.5, gE[0].d.e, kBlack, 1, 3); ses::HLine(-0.5, kNV - 0.5, gE[0].d.e - gE[0].d.err, kBlack, 2, 1); ses::HLine(-0.5, kNV - 0.5, gE[0].d.e + gE[0].d.err, kBlack, 2, 1); }
        if (p == 2) { ses::HLine(-0.5, kNV - 0.5, 4, kBlack, 2, 1); seen::Txt(0.2, 0.86, "circles: OA, squares: M_{ee} classes (pre)", 0.045); }
        if (p == 3 || p == 6 || p == 7 || p == 5) ses::HLine(-0.5, kNV - 0.5, p == 3 ? 0. : 1., kBlack, 2, 1);
        if (p == 4) { ses::HLine(-0.5, kNV - 0.5, gE[2].dTag.first, kBlack, 1, 3); ses::HLine(-0.5, kNV - 0.5, gE[2].dTag.first - gE[2].dTag.second, kBlack, 2, 1); ses::HLine(-0.5, kNV - 0.5, gE[2].dTag.first + gE[2].dTag.second, kBlack, 2, 1); }
        if (p == 5) seen::Txt(0.2, 0.86, "E2; diamonds: E1 bump window", 0.045);
        ses::Pts(p == 0 || p == 1 || p == 4 ? x2 : x0, a, kColAfter, 20, 1.9)->Draw("P same");
        if (!b.empty()) ses::Pts(x1v, b, p == 0 ? kColPeakDirect : kColSig, 21, 1.9)->Draw("P same");
        if (!c.empty()) ses::Pts(x2, c, kColEtaBkg, 33, 2.2)->Draw("P same");
    }
    Save(cv, "rw_compare_variants", false);
}

}  // namespace src

// =====================================================================================================
void sim_rw_compare(int parts = 15) {
    using namespace simeng;
    using namespace src;
    TStopwatch swAll;
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0); gStyle->SetEndErrorSize(0);
    TH1::AddDirectory(kFALSE);
    InitVars();
    gK = fopen(Form("%s/rw_compare_checks.txt", kTabDir), "w");
    K("# rw_compare_checks: sim_rw_compare.C, checks of the recomputed 'off' (and some 'after') numbers against the earlier tables; inputs\n");
    // physics weights
    if (!simphys::Init()) { printf("simphys::Init failed\n"); return; }
    simphys::Install();
    if (!simphys::LoadAngles()) { printf("simphys::LoadAngles failed\n"); return; }
    for (int iv = 0; iv < kNV; ++iv) K("variant %-5s code %3d = %s (%s)\n", gV[iv].tag.c_str(), gV[iv].code, simphys::PhysLabel(gV[iv].code).c_str(), gV[iv].lab.c_str());
    // simulation
    SimData* d = new SimData;
    if (!Load(*d, kCompact, kMaskUS, seen::kLoadLo, seen::kLoadHi, 1)) return;
    gD = d; seen::gD = d;
    if (!seen::AttachER(*d, gK)) { printf("eta-related flags not aligned: STOP\n"); return; }
    SimData* s = new SimData;
    size_t nS = 0;
    for (int k = 0; k < kNTreeAll; ++k) {
        for (const Row& r : d->t[k]) if (r.tmpl == 1 && r.cls >= 2 && r.cls <= 5) s->t[k].push_back(r);
        nS += s->t[k].size();
    }
    s->loaded = d->loaded;
    gS = s;
    K("simulation: rows loaded %zu (std_us %zu, lep_us %zu, eta_us %zu); reduced signal classes 2-5: %zu\n", d->t[kStdUS].size() + d->t[kLepUS].size() + d->t[kEtaUS].size(),
      d->t[kStdUS].size(), d->t[kLepUS].size(), d->t[kEtaUS].size(), nS);
    // data
    scan::Data* dd = new scan::Data;
    if (!scan::Load(*dd, scan::kCompact, seen::kLoadLo, seen::kLoadHi)) return;
    gDD = dd;
    // angle coverage of V3
    {
        simphys::AngHit() = 0; simphys::AngMiss() = 0;
        const SimOpt o = Opt(vMain);
        double sw = 0;
        for (const Row& r : d->t[kEtaUS]) if (r.tmpl == 1 && r.cls == 5) sw += Weight(r, o);
        K("V3 (main): class-5 rows of eta_us (M4 < 1.1, owned): angles found %lld, missing %lld; sum w %.6g\n", simphys::AngHit(), simphys::AngMiss(), sw);
    }

    // shared: simulated peak contents for off and main, data peak fits, efficiencies
    PeakSimFill(vOff); PeakSimFill(vMain); PeakSimFill(vV2);
    PeakDataFits();
    Efficiencies();
    if (parts & 1) {
        PeakPositions();
        PeakTables();
        FigPeakClasses();
        FigPeakEff();
        FigPeakPos();
        printf("part 1 done: %.0f s\n", swAll.RealTime()); swAll.Continue();
    }
    if (parts & (2 | 4 | 8)) {
        EnrichedSpectra();
        BandSums(vOff); BandSums(vMain); BandSums(vSl);
    }
    if (parts & 2) {
        BumpPart();
        if (!(parts & 4)) BumpTables();   // with part 3, written after the template fits (section 6 quotes alpha_bump)
        printf("part 2 done: %.0f s\n", swAll.RealTime()); swAll.Continue();
    }
    if (parts & 4) {
        FitDataPrep();
        Fits(vOff); Fits(vMain);
        MethDataPrep();
        for (int is = 0; is < kNMS; ++is) { gM[vOff][is] = Methods(is, vOff); gM[vMain][is] = Methods(is, vMain); }
        EnrichedTables();
        if (parts & 2) BumpTables();
        FigM4();
        FigRegions();
        FigFits();
        FigMethods();
        if (parts & 2) FigBump();
        printf("part 3 done: %.0f s\n", swAll.RealTime()); swAll.Continue();
    }
    if (parts & 8) {
        if (!(parts & 4)) { FitDataPrep(); MethDataPrep(); }
        Variants();
        FigVariants();
    }
    K("\n%d checks, %d outside the tolerance\n", gNchk, gNbad);
    fclose(gK); gK = nullptr;
    printf("sim_rw_compare: done in %.0f s; checks %d, outside tolerance %d\n", swAll.RealTime(), gNchk, gNbad);
}
