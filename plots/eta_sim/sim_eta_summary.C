// ============================================================================
// sim_eta_summary.C  --  key figures and headline numbers of the eta-set step
// ============================================================================
// HADES p+p, T = 4.53 GeV, pi+ pi- e+ e- (e+e- hypothesis). The SMASH eta set
// (sim_set = 3: the direct decay eta -> pi+pi-e+e-) is added to the std and lep
// sets. User guidance (2026-11-12): SMASH describes the multi-pion background
// above the eta structures badly, so the simulation is judged only in
// eta-enriched selections (mainly p p eta, some p p eta X) with M4 < 1 GeV/c2;
// a mismatch is a result for discussion, not a failure.
//
// One consistent set of key figures (figs/, png + pdf):
//   eta_summary_a_sim_m4   simulated dsigma/dM4 per component, including the direct decay, for E1 and E2
//   eta_summary_b_datasim  data N, CB, N - CB against the simulation in E1 and E2, absolute, ratio pads,
//                          M4 0.35-1.00 (log) and the linear zoom 0.40-0.70
//   eta_summary_c_alpha    fitted factors alpha (bump, peak, background) across selections, fit ranges,
//                          mass scale and ghost variant; the peak ratio data/sim by the data method (P_main)
//   eta_summary_d_classes  peak content per OA and M(ee) class, data against simulation (pre and E2)
//   eta_summary_e_eff      cut efficiencies of the direct decay, the conversion peak and the bump, next to the
//                          data-derived ones, with pulls
//   eta_summary_f_peak_methods  the eta peak, data / simulation by method (added after the review of 2026-11-12)
// Review of 2026-11-12: peak ratio by method (f, tables/eta_summary_peak_methods.txt), fit-model spread on data peak
// contents and efficiencies, direct decay labelled 'generator M(ee) model not QED', QED variant illustration only,
// fit errors scaled by sqrt(chi2/ndf), regions above 0.56 not used as validation (central vertex slice added).
// Tables: tables/eta_summary.txt (headline numbers), tables/eta_summary_checks.txt (cross-checks of every
// recomputed number against the study tables).
//
// Reuse: this file #includes sim_eta_enriched.C (namespace seen: selections E1-E4, eta-related flags and
// their cache, components, template fit, plot helpers) and sim_eta_peak.C (namespace sep: single cuts,
// class sets, data peak fits, the Kroll-Wada (QED) variant of the direct decay, the mixture fit). Neither
// file is changed and their functions sim_eta_enriched() / sim_eta_peak() are not called.
// Recomputed here: figures a, b, d, e (simulation and data), and the M1 fits and P_main values of E1 and E2
// as a check of figure c. Read from tables: the 320 template fits and the P_main values of
// tables/eta_enriched_fit.txt (figure c), the data-derived bump efficiencies (ECAL-tagged eta, edge step)
// and the Dalitz-plot-weighted simulated tagged efficiency of tables/cuts_compare_data.txt (figure e).
//
// Rules: simulation = std (owned, mode A) + lep (owned, mode A) + eta (owned), w (PDG), T1 + T2; variant (a)
// adds the ghost-flagged true pairs (sim_ghostvar.h); sigma = sum w / L_sim, L_sim = 221.05 nb^-1 per set;
// no CB in the simulation; X17 shape for BR = 1e-5 information only, never summed. Data: (N - CB) / L_data,
// CB = 2 sqrt(N++ N--) (10 MeV), L_data = 6080.5 nb^-1 UNCONFIRMED. pre = isBest >= 0, PT3, vz > -160 mm,
// corrected momenta (cor_ok). All weighted sums in double (WSpec, TH1D).
//
// Parts (bits, default 31): 1 figures a and b, 2 figure c, 4 figure d, 8 figure e, 16 figure f; the tables collect what
// the parts computed. Run from /home/przygoda/HADES/FAT/FAT (about 2-3 min):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/aclic_eta_summary", kTRUE)' 'plots/eta_sim/sim_eta_summary.C+'
// ============================================================================
#include <cctype>
#include <cstring>
#include <cstdarg>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>

#include "TBox.h"
#include "TGraphErrors.h"
#include "TMath.h"

#include "sim_eta_enriched.C"
#include "sim_eta_peak.C"

namespace ses {
using namespace simeng;
using seen::Axes;
using seen::Dummy;
using seen::HMax;
using seen::HMinPos;
using seen::MkPad;
using seen::RatioH;
using seen::StyleH;
using seen::StyleM;
using seen::Txt;
using seen::VLines;

// ---------------------------------------------------------------- output
FILE* gT = nullptr;   // tables/eta_summary.txt
FILE* gC = nullptr;   // tables/eta_summary_checks.txt
void T(const char* fmt, ...) {
    char b[32768];
    va_list ap; va_start(ap, fmt); vsnprintf(b, sizeof b, fmt, ap); va_end(ap);
    fputs(b, stdout);
    if (gT) fputs(b, gT);
}
void C(const char* fmt, ...) {
    char b[32768];
    va_list ap; va_start(ap, fmt); vsnprintf(b, sizeof b, fmt, ap); va_end(ap);
    fputs(b, stdout);
    if (gC) fputs(b, gC);
}
int gCnt = 0;
std::string Nm(const char* p) { return Form("ses_%s_%d", p, gCnt++); }
void Save(TCanvas* c, const char* name) {
    c->SaveAs(Form("%s/%s.png", kFigDir, name));
    c->SaveAs(Form("%s/%s.pdf", kFigDir, name));
}
const double kToRows = kLdata / kLsim;
const int kColData = kBlack, kColCB = kRed + 1, kColSig = kBlue + 1;
const double kNaN = std::numeric_limits<double>::quiet_NaN();
typedef std::pair<double, double> VE;
// worst relative / absolute difference of the cross-checks
double gWorstRel = 0, gWorstPull = 0;
void Check(const char* what, double mine, double ref, double tolAbs) {
    if (!(ref == ref)) { C("  %-62s this %12.5g  table: no value\n", what, mine); return; }
    const double d = mine - ref;
    const bool ok = std::fabs(d) <= tolAbs || (!(ref == ref) && !(mine == mine));
    C("  %-62s this %12.5g  table %12.5g  diff %+10.3g  %s\n", what, mine, ref, d, ok ? "OK" : "DIFFERS");
    if (ref == ref && std::fabs(ref) > 0) gWorstRel = std::max(gWorstRel, std::fabs(d / ref));
    if (!ok) gWorstPull += 1;
}

// ---------------------------------------------------------------- reading the study tables
std::string Collapse(const std::string& s) {
    std::string o; bool sp = false;
    for (char ch : s) {
        if (ch == ' ' || ch == '\t') { sp = true; continue; }
        if (sp && !o.empty()) o += ' ';
        sp = false; o += ch;
    }
    return o;
}
const std::vector<std::string>& Lines(const char* file) {   // cached; const char* avoids a temporary argument
    static std::map<std::string, std::vector<std::string>> cache;
    auto it = cache.find(file);
    if (it != cache.end()) return it->second;
    std::vector<std::string>& v = cache[file];
    std::ifstream in((std::string(kTabDir) + "/" + file).c_str());
    std::string l;
    while (std::getline(in, l)) v.push_back(l);
    if (v.empty()) printf("ses: table %s missing or empty\n", file);
    return v;
}
// first line of section sec ("## <sec>:" or "## <sec> ..."; empty = whole file) whose collapsed text starts
// with the collapsed `start` followed by a blank or the end of the line
bool FindLine(const std::string& file, const std::string& sec, const std::string& start, std::string& out) {
    bool in = sec.empty();
    const std::string st = Collapse(start);
    for (const std::string& l : Lines(file.c_str())) {
        if (l.compare(0, 3, "## ") == 0) {
            if (!sec.empty())
                in = l.compare(0, 4 + sec.size(), "## " + sec + ":") == 0 || l.compare(0, 4 + sec.size(), "## " + sec + " ") == 0;
            continue;
        }
        if (!in || l.empty() || l[0] == '#') continue;
        const std::string c = Collapse(l);
        if (c.compare(0, st.size(), st) == 0 && (c.size() == st.size() || c[st.size()] == ' ')) { out = l; return true; }
    }
    return false;
}
// value +- error pairs of a line per '|' group; a lone '-' is a missing value (NaN pair)
std::vector<std::vector<VE>> Groups(const std::string& line) {
    std::vector<std::string> tok;
    std::string cur;
    for (char ch : line) {
        if (ch == '|') { if (!cur.empty()) tok.push_back(cur); cur.clear(); tok.push_back("|"); }
        else if (ch == ' ' || ch == '\t') { if (!cur.empty()) tok.push_back(cur); cur.clear(); }
        else cur += ch;
    }
    if (!cur.empty()) tok.push_back(cur);
    std::vector<std::vector<VE>> g(1);
    for (size_t i = 0; i < tok.size(); ++i) {
        const std::string& t = tok[i];
        if (t == "|") { g.push_back(std::vector<VE>()); continue; }
        if (t == "-") { g.back().push_back(VE(kNaN, kNaN)); continue; }
        const size_t p = t.find("+-");
        if (p != std::string::npos && p > 0 && p + 2 < t.size()) {
            g.back().push_back(VE(atof(t.substr(0, p).c_str()), atof(t.substr(p + 2).c_str())));
            continue;
        }
        if (i + 2 < tok.size() && tok[i + 1] == "+-") {
            g.back().push_back(VE(atof(t.c_str()), atof(tok[i + 2].c_str())));
            i += 2;
        }
    }
    return g;
}
VE GetVE(const std::string& file, const std::string& sec, const std::string& start, int group, int idx) {
    std::string l;
    if (!FindLine(file, sec, start, l)) return VE(kNaN, kNaN);
    const std::vector<std::vector<VE>> g = Groups(l);
    if (group < 0) {   // all pairs of the line, in order
        std::vector<VE> all;
        for (const auto& x : g) all.insert(all.end(), x.begin(), x.end());
        return idx < (int)all.size() ? all[idx] : VE(kNaN, kNaN);
    }
    if (group >= (int)g.size() || idx >= (int)g[group].size()) return VE(kNaN, kNaN);
    return g[group][idx];
}

// template fits of tables/eta_enriched_fit.txt
struct FitRow {
    std::string sel; int rng = 0; bool k0s = false, a = false; std::string model;
    std::vector<std::string> pn; std::vector<VE> p; double chi2 = 0; int ndf = 0; VE pk, bm;
    VE Par(const char* n) const {
        for (size_t i = 0; i < pn.size(); ++i) if (pn[i] == n) return p[i];
        return VE(kNaN, kNaN);
    }
};
std::vector<FitRow> ParseFits() {
    std::vector<FitRow> out;
    std::string sec;
    for (const std::string& l : Lines("eta_enriched_fit.txt")) {
        if (l.compare(0, 3, "## ") == 0) {
            std::istringstream is(l.substr(3));
            is >> sec;
            if (!sec.empty() && sec[sec.size() - 1] == ':') sec.erase(sec.size() - 1); else sec.clear();
            continue;
        }
        if (sec.empty() || l.empty() || !isdigit((unsigned char)l[0])) continue;
        std::istringstream is(l);
        FitRow f; f.sel = sec;
        std::string range, scale, var;
        is >> f.rng >> range >> scale >> var >> f.model;
        f.k0s = scale == "K0S"; f.a = var == "(a)";
        std::string t;
        while (is >> t) {
            if (t.find('/') != std::string::npos) { sscanf(t.c_str(), "%lf/%d", &f.chi2, &f.ndf); break; }
            std::string v; is >> v;
            double a = 0, e = 0;
            if (sscanf(v.c_str(), "%lf+-%lf", &a, &e) == 2) { f.pn.push_back(t); f.p.push_back(VE(a, e)); }
        }
        std::string pv, pm; double x = kNaN, ex = kNaN, y = kNaN, ey = kNaN;
        is >> pv >> x >> pm >> ex >> y >> pm >> ey;
        f.pk = VE(x, ex); f.bm = VE(y, ey);
        out.push_back(f);
    }
    return out;
}
const FitRow* FindFit(const std::vector<FitRow>& v, const std::string& sel, int rng, bool k0s, bool a, const char* model) {
    for (const FitRow& f : v) if (f.sel == sel && f.rng == rng && f.k0s == k0s && f.a == a && f.model == model) return &f;
    return nullptr;
}

// ---------------------------------------------------------------- shared state
SimData* gD = nullptr;        // simulation, e+e- trees of the three sets, owned rows, M4 < 1.1
SimData* gDS = nullptr;       // reduced copy: signal classes 2-6 (template 1) only
const scan::Data* gDD = nullptr;
SimOpt gOpt;                  // std + lep + eta, mode A, w, e+e-
const int kJ[2] = {seen::sE1, seen::sE2};   // the two main selections
const char* kJTag[2] = {"E1", "E2"};
const char* kJName[2] = {"E1 (slides)", "E2 (slides without OA)"};

// figures a and b
std::vector<Comp> gComps;
std::vector<WSpec> sLog[2], sLin[2], sSg[2], sReg[2];
std::vector<WSum> sTot[2];
DSpec dLog[2], dLin[2];
std::vector<scan::Spectra> dReg;
double rat[2][8], ratE[2][8], ratA[2][8], ratAE[2][8], dNb[2][8], dNbE[2][8], sNb[2][8], sNbE[2][8], aNb[2][8], aNbE[2][8];
double pur[2], purE[2], purA[2], purAE[2];
// central target slice -140 < vz < -100 mm (review finding 6): data/sim per region, T1 + T2 and variant (a)
double ratC[2][8], ratCE[2][8], ratCA[2][8];
std::vector<EffRes> purReg[2];
bool gHaveAB = false;

// figure c
std::vector<FitRow> gFits;
struct PMRow { VE data, dataHand, sim54, sim62, simA, pkOnly, simUns, r54, r62, rA, rUns; bool ok = false; };
PMRow gPM[seen::kNSel];
bool gHaveC = false;

// figure d
struct ClassRes {
    std::string sel, cs; std::vector<std::string> name;
    std::vector<double> sd, sde, sc, sce, st, ste, sk, ske, sx, dS, dSe, dSs, dSm;   // dSe = stat (+) model (review finding 4)
    sep::MixRes mix, mixKW, mixStat;
    double dTot = 0, dTotE = 0;     // data P_main of the selection (no class)
};
std::vector<ClassRes> gCls;
bool gHaveD = false;

// figure e
struct CutDef { const char* lab; const char* txt; int kind; double val; const char* keyPeak; const char* keyBump; bool bumpStart; bool sculpt; };
const int kNCut = 11;
const CutDef kCut[kNCut] = {
    {"OA > 2#circ", "OA > 2", sep::kOaLo, 2, "oa_lo 2 pre", "oa_lo 2 no cut", false, false},
    {"OA < 4#circ", "OA < 4", sep::kOaHi, 4, "oa_hi 4 pre", "oa_hi 4 no cut", false, false},
    {"#alpha_{lab} < 50#circ", "oa_lab < 50", sep::kLab, 50, "oa_lab 50 pre", "oa_lab 50 no cut", false, false},
    {"M_{#pi#pi} > 0.34", "M(pipi) > 0.34", sep::kMppLo, 0.34, "mpipi_lo 0.34 pre", "mpipi_lo 0.34 no cut", false, false},
    {"M_{#pi#pi} < 0.40", "M(pipi) < 0.40", sep::kMppHi, 0.40, "mpipi_hi 0.40 pre", "mpipi_hi 0.40 no cut", false, true},
    {"M_{#pi#pi} < 0.46", "M(pipi) < 0.46", sep::kMppHi, 0.46, "mpipi_hi 0.46 pre", "mpipi_hi 0.46 no cut", false, true},
    {"#alpha_{#eta} > 140#circ", "eta_rest > 140", sep::kEta, 140, "eta_rest 140 pre", "eta_rest 140 no cut", false, false},
    {"#alpha_{#eta} > 160#circ", "eta_rest > 160", sep::kEta, 160, "eta_rest 160 pre", "eta_rest 160 no cut", false, false},
    {"MM_{4} < 2.80", "MM4 < 2.80", sep::kMm4, 2.80, "mm4_hi 2.80 pre", "mm4_hi 2.80 no cut", false, true},
    {"p_{T,4} > 400", "pT4 > 400 MeV/c", sep::kPt4, 400, "pt4_lo 400 pre", "pt4_lo 400 start", true, false},
    {"both #pi in vertex", "both pions in vertex fit", sep::kVtx, 2, "vtx both pions pre", "vtx both pions start", true, false}};
enum { eDir = 0, eConv, eTot, eTotKW, eBump, eTag, eSbu, eNE };
struct EffOut {
    EffRes s[eNE];                  // simulation at base pre
    EffRes sBumpStart, sTagStart;   // bump all and tagged at the bump start base (pT4, vertex)
    VE dPeak, dPeakTab, dTag, dEdge, tTag, tSbu, tBump, tTagDal, tDir, tConv, tTot, tTotKW;
    bool degTag = false, degEdge = false;   // data value with zero error in the study table (degenerate fit): not used
};
EffOut gEff[kNCut];
bool gHaveE = false;

// ---------------------------------------------------------------- small helpers
std::string SelDefLine() {
    return "E1 (slides) = pre + #alpha_{lab} < 50#circ, #alpha_{#eta} > 140#circ, M_{#pi#pi} < 0.46 GeV/c^{2}, OA > 2#circ;   "
           "E2 = E1 without OA > 2#circ;   pre = isBest #geq 0, PT3, v_{z} > -160 mm, corrected momenta, e^{+}e^{-}";
}
std::string SimLine() {
    return Form("simulation: SMASH std + lep + eta, owned rows (mode A), w (PDG), #Sigma w / L_{sim} / bin width, L_{sim} = %.2f nb^{-1} per set; "
                "T1 + T2, no CB in simulation; no efficiency correction", kLsim);
}
std::string DataLine() { return std::string("data: rows / L_{data} / bin width, ") + seen::kLdTex + "; CB = 2#sqrt{N_{++}N_{--}} in 10 MeV bins"; }
TH1D* Frame(const char* name, int n, double lo, double hi, const char* xt, const char* yt, double ylo, double yhi, double ts = 0.05,
            double yoff = 1.2) {
    TH1D* f = new TH1D(name, "", n, lo, hi); f->SetDirectory(nullptr);
    Axes(f, xt, yt, ts, yoff);
    f->SetMinimum(ylo); f->SetMaximum(yhi);
    return f;
}
TH1D* CatFrame(const char* name, const std::vector<std::string>& labels, const char* yt, double ylo, double yhi, double ts, double lsz,
               double yoff = 1.0) {
    const int n = (int)labels.size();
    TH1D* f = new TH1D(name, "", n, -0.5, n - 0.5); f->SetDirectory(nullptr);
    for (int i = 0; i < n; ++i) f->GetXaxis()->SetBinLabel(i + 1, labels[i].c_str());
    f->GetXaxis()->SetLabelSize(lsz); f->GetXaxis()->LabelsOption("h");
    f->GetYaxis()->SetTitle(yt); f->GetYaxis()->SetTitleSize(ts); f->GetYaxis()->SetLabelSize(ts * 0.85);
    f->GetYaxis()->SetTitleOffset(yoff); f->GetYaxis()->SetNdivisions(505);
    f->SetMinimum(ylo); f->SetMaximum(yhi);
    return f;
}
TGraphErrors* Pts(const std::vector<double>& x, const std::vector<VE>& y, int col, int mk, double ms, int lw = 2) {
    TGraphErrors* g = new TGraphErrors();
    for (size_t i = 0; i < x.size(); ++i) {
        if (!(y[i].first == y[i].first)) continue;
        const int n = g->GetN();
        g->SetPoint(n, x[i], y[i].first);
        g->SetPointError(n, 0., y[i].second == y[i].second ? y[i].second : 0.);
    }
    g->SetMarkerColor(col); g->SetLineColor(col); g->SetMarkerStyle(mk); g->SetMarkerSize(ms); g->SetLineWidth(lw);
    return g;
}
void HLine(double x1, double x2, double y, int col, int st, int w) {
    TLine* l = new TLine(x1, y, x2, y); l->SetLineColor(col); l->SetLineStyle(st); l->SetLineWidth(w); l->Draw();
}
std::string F2(VE v, int p = 2) { return Form("%.*f #pm %.*f", p, v.first, p, v.second); }
// fraction f_i = x_i / sum x of independent class contents: var f_i = [(1 - f_i)^2 e_i^2 + f_i^2 sum_{j != i} e_j^2] / T^2
VE Frac(const std::vector<double>& x, const std::vector<double>& e, size_t i) {
    double t = 0, v = 0;
    for (size_t j = 0; j < x.size(); ++j) { t += x[j]; v += e[j] * e[j]; }
    const double f = x[i] / t;
    return VE(f, std::sqrt((1 - f) * (1 - f) * e[i] * e[i] + f * f * (v - e[i] * e[i])) / t);
}

// ============================================================================
// Shared fills of E1 and E2 (figures a and b)
// ============================================================================
void FillAB() {
    TStopwatch sw;
    const SimData& d = *gD;
    gComps = seen::Comps();
    Binning bLog; bLog.nbin = 130; bLog.lo = 0.35; bLog.hi = 1.00;
    Binning bLin; bLin.nbin = 120; bLin.lo = 0.40; bLin.hi = 0.70;
    Binning bSg; bSg.nbin = 90; bSg.lo = 0.42; bSg.hi = 0.60;
    const Binning bReg = seen::RegBinning();
    for (int j = 0; j < 2; ++j) {
        const Sel sl = seen::SimSel(kJ[j]);
        sLog[j] = FillMany(d, sl, gComps, gOpt, bLog);
        sLin[j] = FillMany(d, sl, gComps, gOpt, bLin);
        sSg[j] = FillMany(d, sl, gComps, gOpt, bSg);
        sReg[j] = FillMany(d, sl, gComps, gOpt, bReg);
        sTot[j] = SumMany(d, seen::SimSelM(kJ[j], -1e9, 1.0), gComps, gOpt);
        dLog[j] = DataSpectrum(*gDD, seen::DatSel(kJ[j]), bLog.nbin, bLog.lo, bLog.hi);
        dLin[j] = DataSpectrum(*gDD, seen::DatSel(kJ[j]), bLin.nbin, bLin.lo, bLin.hi);
        // purity of p p eta (X) (seen definition)
        const Sel base = seen::SimSelM(kJ[j], -1e9, 1.0);
        const Sel pass = [](const Row& r) { return seen::EtaX(r); };
        const EffRes e = Efficiency(d, base, pass, CompSim(), gOpt);
        const EffRes ea = Efficiency(d, base, pass, CompSimVar(1), gOpt);
        pur[j] = e.eps; purE[j] = e.err; purA[j] = ea.eps; purAE[j] = ea.err;
        purReg[j] = EfficiencyVs(d, base, pass, CompSim(), gOpt, bReg);
    }
    scan::FillOpt fo; fo.nbin = 260; fo.lo = 0.35; fo.hi = 1.00; fo.reb = 4;
    for (int k = 0; k <= seen::kNReg; ++k) fo.windows.push_back(std::make_pair(seen::kReg[k].lo, seen::kReg[k].hi));
    dReg = scan::FillMany(*gDD, {seen::DatSel(kJ[0]), seen::DatSel(kJ[1])}, fo);
    for (int j = 0; j < 2; ++j)
        for (int k = 0; k <= seen::kNReg; ++k) {
            const scan::WinSum& w = dReg[j].win[2 + k];
            const double dn = w.sig() / kLdata, de = w.sigErr() / kLdata;
            double ss, sv, sa, sav;
            if (k < seen::kNReg) { ss = sReg[j][seen::cSum].s[k]; sv = sReg[j][seen::cSum].v[k]; sa = sReg[j][seen::cSumA].s[k]; sav = sReg[j][seen::cSumA].v[k]; }
            else { ss = sReg[j][seen::cSum].S; sv = sReg[j][seen::cSum].V; sa = sReg[j][seen::cSumA].S; sav = sReg[j][seen::cSumA].V; }
            ss /= kLsim; sa /= kLsim; sv /= kLsim * kLsim; sav /= kLsim * kLsim;
            dNb[j][k] = dn; dNbE[j][k] = de; sNb[j][k] = ss; sNbE[j][k] = std::sqrt(sv); aNb[j][k] = sa; aNbE[j][k] = std::sqrt(sav);
            rat[j][k] = dn / ss; ratE[j][k] = rat[j][k] * std::sqrt(de * de / (dn * dn) + sv / (ss * ss));
            ratA[j][k] = dn / sa; ratAE[j][k] = ratA[j][k] * std::sqrt(de * de / (dn * dn) + sav / (sa * sa));
        }
    // ---- central target slice (review finding 6): the region ratios depend on the vertex slice
    {
        const std::vector<Comp> cc = {CompSim(), CompSimVar(1)};
        for (int j = 0; j < 2; ++j) {
            const int is = kJ[j];
            const Sel sl = SelAnd(seen::SimSel(is), [](const Row& r) { return r.vz > -140.f && r.vz < -100.f; });
            const std::vector<WSpec> sc = FillMany(d, sl, cc, gOpt, bReg);
            const std::vector<WSum> st = SumMany(d, SelAnd(sl, [](const Row& r) { return r.m4 >= 0.35f && r.m4 < 1.0f; }), cc, gOpt);
            const scan::Spectra sp = scan::FillMany(*gDD, {[is](const scan::Row& r) { return seen::PassSel(is, r) && r.vz > -140.f && r.vz < -100.f; }}, fo)[0];
            for (int k = 0; k <= seen::kNReg; ++k) {
                const scan::WinSum& w = sp.win[2 + k];
                const double dn = w.sig(), de = w.sigErr();
                const double ss = (k < seen::kNReg ? sc[0].s[k] : st[0].S) * kToRows, sv = (k < seen::kNReg ? sc[0].v[k] : st[0].V) * kToRows * kToRows;
                const double sa = (k < seen::kNReg ? sc[1].s[k] : st[1].S) * kToRows;
                ratC[j][k] = dn / ss; ratCE[j][k] = ratC[j][k] * std::sqrt(de * de / (dn * dn) + sv / (ss * ss)); ratCA[j][k] = dn / sa;
            }
        }
    }
    gHaveAB = true;
    printf("FillAB: %.1f s\n", sw.RealTime());

    // ---- checks against the enriched study tables
    C("\n## figures a, b: E1 / E2 against tables/eta_enriched_regions.txt and eta_enriched_purity.txt\n");
    for (int j = 0; j < 2; ++j) {
        for (int k = 0; k <= seen::kNReg; ++k) {
            std::string reg = seen::kReg[k].name;
            Check(Form("%s %s: data N-CB [nb]", kJTag[j], reg.c_str()), dNb[j][k], GetVE("eta_enriched_regions.txt", kJTag[j], reg, 0, 0).first, 6e-5);
            Check(Form("%s %s: sim T1+T2 [nb]", kJTag[j], reg.c_str()), sNb[j][k], GetVE("eta_enriched_regions.txt", kJTag[j], reg, 0, 1).first, 6e-5);
            Check(Form("%s %s: data/sim", kJTag[j], reg.c_str()), rat[j][k], GetVE("eta_enriched_regions.txt", kJTag[j], reg, 0, 3).first, 6e-5);
            Check(Form("%s %s: data/sim (a)", kJTag[j], reg.c_str()), ratA[j][k], GetVE("eta_enriched_regions.txt", kJTag[j], reg, 0, 4).first, 6e-5);
        }
        Check(Form("%s purity M4 < 1.0", kJTag[j]), pur[j], GetVE("eta_enriched_purity.txt", "", kJTag[j], -1, 0).first, 6e-5);
        Check(Form("%s purity M4 < 1.0 (a)", kJTag[j]), purA[j], GetVE("eta_enriched_purity.txt", "", kJTag[j], -1, 1).first, 6e-5);
        const int cidx[6] = {seen::cBD, seen::cBC, seen::cPC, seen::cPD, seen::cKE, seen::cKM};
        for (int q = 0; q < 6; ++q) {
            // composition block: the second line of the purity table that starts with the tag
            std::string l; bool first = true, found = false;
            for (const std::string& x : Lines("eta_enriched_purity.txt")) {
                const std::string c = Collapse(x);
                if (c.compare(0, strlen(kJTag[j]) + 1, std::string(kJTag[j]) + " ") == 0) {
                    if (first) { first = false; continue; }
                    l = x; found = true; break;
                }
            }
            const std::vector<std::vector<VE>> g = found ? Groups(l) : std::vector<std::vector<VE>>(1);
            const double ref = (found && q < (int)g[0].size()) ? g[0][q].first : kNaN;
            Check(Form("%s composition M4 < 1.0: %s [nb]", kJTag[j], seen::CompShort(cidx[q])), sTot[j][cidx[q]].Xsec(), ref, 6e-5);
        }
    }
    // engine reference: slides (= E1), peak window, read from tables/eta_engine_reference.txt (sim_eta_reference.C)
    // (rules v1 hard-coded 0.9642 and 5.243; rules v2 (2026-11-13) 0.9803 and 5.253)
    const double pt = (sReg[0][seen::cPC].s[3] + sReg[0][seen::cPD].s[3]) / kLsim;
    const std::string erSec = "selection \"slides\", peak";
    Check("E1 peak window: peak total [nb] (eta_engine_reference, slides)", pt, GetVE("eta_engine_reference.txt", erSec, "peak total (C4 + C5)", 0, 0).first, 6e-4);
    Check("E1 peak window: T1+T2 [nb] (eta_engine_reference, slides)", sReg[0][seen::cSum].s[3] / kLsim,
          GetVE("eta_engine_reference.txt", erSec, "T1+T2 full simulated spectrum", 0, 0).first, 6e-4);
}

// ============================================================================
// Figure a: simulated M4 per component
// ============================================================================
void FigA() {
    using namespace seen;
    TCanvas* cv = new TCanvas(Nm("cA").c_str(), "", 2700, 2100);
    cv->cd();
    MkPad(Nm("hA").c_str(), 0, 0.87, 1, 1, false, 0, 0, 0, 0);
    Txt(0.01, 0.84, "Simulated d#sigma/dM_{#pi^{+}#pi^{-}e^{+}e^{-}} per component, with the direct decay #eta #rightarrow #pi^{+}#pi^{-}e^{+}e^{-} (eta set; generator M_{ee} model "
                    "not QED), in the #eta-enriched selections E1 and E2", 0.120);
    Txt(0.01, 0.69, SimLine().c_str(), 0.08);
    Txt(0.01, 0.57, SelDefLine().c_str(), 0.08);
    Txt(0.01, 0.45, NoteX17(1e-5).c_str(), 0.08, kColX17);
    TLegend* lg = new TLegend(0.0, 0.0, 1.0, 0.40);
    lg->SetNColumns(5); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.078);
    for (int k = 5; k >= 0; --k) lg->AddEntry(Dummy(Nm("lg").c_str(), kStackLine[k], 2, 1, kStackFill[k]), kStackLeg[k], "f");
    lg->AddEntry(Dummy(Nm("lg").c_str(), kColPeakTotal, 4, 1, 0), "peak total = conversion + direct", "l");
    lg->AddEntry(Dummy(Nm("lg").c_str(), kColSim, 4, 1, 0), "simulation T1 + T2 (sum)", "l");
    lg->AddEntry(Dummy(Nm("lg").c_str(), kColSim, 3, 2, 0), "variant (a): + ghost-flagged true pairs", "l");
    lg->AddEntry(Dummy(Nm("lg").c_str(), kColX17, 3, 2, 0), "X17 shape, BR 10^{-5} (information, not in the sum)", "l");
    lg->Draw();

    // common y ranges of the two rows
    double yLogHi = 0, yLogLo = 1e30, yLinHi = 0, ySgHi = 0;
    for (int j = 0; j < 2; ++j) {
        yLogHi = std::max(yLogHi, HMax(sLog[j][cSumA].Hist(Nm("t").c_str(), 2)));
        yLogLo = std::min(yLogLo, HMinPos(sLog[j][cPD].Hist(Nm("t").c_str(), 2)));
        yLinHi = std::max(yLinHi, HMax(sLin[j][cSumA].Hist(Nm("t").c_str(), 2)));
        for (int c : {cBD, cBC, cPC, cPD}) ySgHi = std::max(ySgHi, HMax(sSg[j][c].Hist(Nm("t").c_str(), 2)));
        TH1D* pt = sSg[j][cPC].Hist(Nm("t").c_str(), 2); pt->Add(sSg[j][cPD].Hist(Nm("t").c_str(), 2));
        ySgHi = std::max(ySgHi, HMax(pt));
    }
    yLogLo = std::max(yLogLo * 0.5, yLogHi * 2e-5);
    const char letters[2][3] = {{'a', 'b', 'c'}, {'d', 'e', 'f'}};
    for (int j = 0; j < 2; ++j) {
        const double y1 = j == 0 ? 0.435 : 0.0, y2 = j == 0 ? 0.87 : 0.435;
        // (1) log, all components
        cv->cd(); MkPad(Nm("pa").c_str(), 0.0, y1, 0.345, y2, true, 0.14, 0.04, 0.09, 0.12);
        TH1D* fr = Frame(Nm("fr").c_str(), 130, 0.35, 1.00, kXM4, kYM4, yLogLo, yLogHi * 40, 0.05, 1.3);
        fr->Draw("axis");
        TH1D* hk = sLog[j][cKM].Hist(Nm("h").c_str(), 2); StyleH(hk, kColPhysBkg, 1, 1, kColPhysFill); hk->Draw("hist same");
        for (int k = 1; k < 6; ++k) {
            TH1D* h = sLog[j][kStack[k]].Hist(Nm("h").c_str(), 2); StyleH(h, kStackLine[k], 3, 1); h->Draw("hist same");
        }
        TH1D* hx = sLog[j][cX17].Hist(Nm("h").c_str(), 2); StyleH(hx, kColX17, 3, 2); hx->Draw("hist same");
        TH1D* hs = sLog[j][cSum].Hist(Nm("h").c_str(), 2); StyleH(hs, kColSim, 4, 1); hs->Draw("hist same");
        TH1D* ha = sLog[j][cSumA].Hist(Nm("h").c_str(), 2); StyleH(ha, kColSim, 3, 2); ha->Draw("hist same");
        VLines({scan::kBumpLo, scan::kBumpHi, scan::kPeakLo, scan::kPeakHi}, yLogLo, yLogHi * 40);
        fr->Draw("axis same");
        Txt(0.15, 0.93, Form("(%c) %s: log scale, 5 MeV bins", letters[j][0], kJTag[j]), 0.046);
        Txt(0.17, 0.85, Form("M_{4} < 1: T1 + T2 = %.1f #pm %.1f nb", sTot[j][cSum].Xsec(), sTot[j][cSum].XsecErr()), 0.036);
        Txt(0.17, 0.80, Form("purity p p #eta (X) = %.3f #pm %.3f", pur[j], purE[j]), 0.036, kColPeakTotal);
        Txt(0.17, 0.75, "dotted: bump and peak windows", 0.032, kGray + 2);
        // (2) linear zoom, stacked
        cv->cd(); MkPad(Nm("pb").c_str(), 0.345, y1, 0.67, y2, false, 0.13, 0.04, 0.09, 0.12);
        TH1D* fb = Frame(Nm("fr").c_str(), 120, 0.40, 0.70, kXM4, kYM4, 0, yLinHi * 1.12, 0.05, 1.25);
        fb->Draw("axis");
        TH1D* st[6];
        for (int k = 0; k < 6; ++k) {
            st[k] = sLin[j][kStack[k]].Hist(Nm("h").c_str(), 2);
            if (k > 0) st[k]->Add(st[k - 1]);
            StyleH(st[k], kStackLine[k], 1, 1, kStackFill[k]);
        }
        for (int k = 5; k >= 0; --k) st[k]->Draw("hist same");
        TH1D* hs2 = sLin[j][cSum].Hist(Nm("h").c_str(), 2); StyleH(hs2, kColSim, 4, 1); hs2->Draw("hist same");
        TH1D* ha2 = sLin[j][cSumA].Hist(Nm("h").c_str(), 2); StyleH(ha2, kColSim, 3, 2); ha2->Draw("hist same");
        TH1D* hx2 = sLin[j][cX17].Hist(Nm("h").c_str(), 2); StyleH(hx2, kColX17, 3, 2); hx2->Draw("hist same");
        VLines({scan::kBumpLo, scan::kBumpHi, scan::kPeakLo, scan::kPeakHi}, 0, yLinHi * 1.12);
        fb->Draw("axis same");
        Txt(0.14, 0.93, Form("(%c) %s: linear, 2.5 MeV bins, stacked", letters[j][1], kJTag[j]), 0.046);
        // (3) signals alone
        cv->cd(); MkPad(Nm("pc").c_str(), 0.67, y1, 1.0, y2, false, 0.13, 0.04, 0.09, 0.12);
        TH1D* fc = Frame(Nm("fr").c_str(), 90, 0.42, 0.60, kXM4, kYM4, 0, ySgHi * 1.55, 0.05, 1.25);
        fc->Draw("axis");
        const int sc[4] = {cBD, cBC, cPC, cPD};
        const int scol[4] = {kColBumpDalitz, kColBumpConv, kColPeakConv, kColPeakDirect};
        for (int q = 0; q < 4; ++q) {
            TH1D* h = sSg[j][sc[q]].Hist(Nm("h").c_str(), 2); StyleH(h, scol[q], 3, 1); h->Draw("hist same");
            TH1D* he = (TH1D*)h->Clone(Nm("he").c_str()); he->SetLineWidth(1); he->Draw("E1 same");
        }
        TH1D* ht = sSg[j][cPC].Hist(Nm("h").c_str(), 2); ht->Add(sSg[j][cPD].Hist(Nm("h").c_str(), 2));
        StyleH(ht, kColPeakTotal, 4, 1); ht->Draw("hist same");
        TH1D* hx3 = sSg[j][cX17].Hist(Nm("h").c_str(), 2); StyleH(hx3, kColX17, 3, 2); hx3->Draw("hist same");
        VLines({scan::kBumpLo, scan::kBumpHi, scan::kPeakLo, scan::kPeakHi}, 0, ySgHi * 1.55);
        fc->Draw("axis same");
        Txt(0.14, 0.93, Form("(%c) %s: the signals alone, 2 MeV bins", letters[j][2], kJTag[j]), 0.046);
        auto win = [&](int c, int k, double& v, double& e) { v = sReg[j][c].s[k] / kLsim; e = std::sqrt(sReg[j][c].v[k]) / kLsim; };
        double a, ea, b, eb, c1, ec, dd1, ed;
        win(cBD, 1, a, ea); win(cBC, 1, b, eb);
        Txt(0.15, 0.86, Form("bump window [nb]: #pi^{0} Dalitz %.2f#pm%.2f, conv. %.2f#pm%.2f", a, ea, b, eb), 0.031);
        win(cPC, 3, c1, ec); win(cPD, 3, dd1, ed);
        Txt(0.15, 0.81, Form("peak window [nb]: conv. %.3f#pm%.3f, direct %.3f#pm%.3f", c1, ec, dd1, ed), 0.031);
        Txt(0.15, 0.76, Form("peak total %.3f#pm%.3f nb = %.0f#pm%.0f data rows", c1 + dd1, std::sqrt(ec * ec + ed * ed),
                             (c1 + dd1) * kLdata, std::sqrt(ec * ec + ed * ed) * kLdata), 0.031, kColPeakTotal);
        Txt(0.15, 0.71, Form("X17 (BR 10^{-5}, information): %.3f nb in the peak window", sReg[j][cX17].s[3] / kLsim), 0.031, kColX17);
    }
    Save(cv, "eta_summary_a_sim_m4");
}

// ============================================================================
// Figure b: data against simulation, absolute, with ratio pads
// ============================================================================
void FigB() {
    using namespace seen;
    TCanvas* cv = new TCanvas(Nm("cB").c_str(), "", 2700, 2400);
    cv->cd();
    MkPad(Nm("hB").c_str(), 0, 0.885, 1, 1, false, 0, 0, 0, 0);
    Txt(0.01, 0.84, "d#sigma/dM_{#pi^{+}#pi^{-}e^{+}e^{-}}: data against the full simulation (with the direct decay) in the #eta-enriched selections E1 and E2, absolute scale",
        0.13);
    Txt(0.01, 0.69, DataLine().c_str(), 0.08);
    Txt(0.01, 0.58, SimLine().c_str(), 0.08);
    Txt(0.01, 0.47, (SelDefLine() + ";   ratio pads: data (N - CB) / simulation, full blue T1 + T2, open teal variant (a); horizontal bars: region averages").c_str(), 0.08);
    SharedLegend(0.0, 0.0, 1.0, 0.42, 6, 0.075);

    double yLogHi = 0, yLogLo = 1e30, yLinHi = 0;
    for (int j = 0; j < 2; ++j) {
        yLogHi = std::max(yLogHi, HMax(dLog[j].Hist("us", Nm("t").c_str(), 2)));
        yLogLo = std::min(yLogLo, HMinPos(dLog[j].Hist("cb", Nm("t").c_str(), 2)));
        yLinHi = std::max(yLinHi, std::max(HMax(dLin[j].Hist("sig", Nm("t").c_str(), 2)), HMax(sLin[j][cSumA].Hist(Nm("t").c_str(), 2))));
    }
    yLogLo = std::max(yLogLo * 0.3, yLogHi * 1e-4);
    const char* qn[7] = {"0.35-0.455", "bump", "gap", "peak", "0.554-0.56", "0.56-0.70", "0.70-1.00"};
    for (int j = 0; j < 2; ++j) {
        const double y1 = j == 0 ? 0.4425 : 0.0, y2 = j == 0 ? 0.885 : 0.4425, ys = y1 + 0.30 * (y2 - y1);
        for (int v = 0; v < 2; ++v) {
            const bool lin = v == 1;
            const double x1 = lin ? 0.5 : 0.0, x2 = lin ? 1.0 : 0.5;
            const std::vector<WSpec>& S = lin ? sLin[j] : sLog[j];
            const DSpec& D = lin ? dLin[j] : dLog[j];
            const double xlo = lin ? 0.40 : 0.35, xhi = lin ? 0.70 : 1.00;
            cv->cd(); MkPad(Nm("pm").c_str(), x1, ys, x2, y2, !lin, 0.12, 0.02, 0.10, 0.0);
            TH1D* hN = D.Hist("us", Nm("hN").c_str(), 2);
            TH1D* hC = D.Hist("cb", Nm("hC").c_str(), 2);
            TH1D* hS = D.Hist("sig", Nm("hS").c_str(), 2);
            TH1D* hSum = S[cSum].Hist(Nm("hSum").c_str(), 2);
            TH1D* hSumA = S[cSumA].Hist(Nm("hSumA").c_str(), 2);
            TH1D* fr = Frame(Nm("fr").c_str(), S[cSum].nbin, xlo, xhi, kXM4, kYM4, lin ? 0 : yLogLo, lin ? yLinHi * 1.15 : yLogHi * 60, 0.055, 1.05);
            fr->Draw("axis");
            if (lin) {
                TH1D* hs[6];
                for (int k = 0; k < 6; ++k) {
                    hs[k] = S[kStack[k]].Hist(Nm("hst").c_str(), 2);
                    if (k > 0) hs[k]->Add(hs[k - 1]);
                    StyleH(hs[k], kStackLine[k], 1, 1, kStackFill[k]);
                }
                for (int k = 5; k >= 0; --k) hs[k]->Draw("hist same");
            } else {
                TH1D* hk = S[cKM].Hist(Nm("hk").c_str(), 2); StyleH(hk, kColPhysBkg, 1, 1, kColPhysFill); hk->Draw("hist same");
                for (int k = 1; k < 6; ++k) { TH1D* h = S[kStack[k]].Hist(Nm("hl").c_str(), 2); StyleH(h, kStackLine[k], 3, 1); h->Draw("hist same"); }
                StyleM(hN, kColData, 20, 0.9); hN->Draw("E1 same");
                StyleM(hC, kColCB, 21, 0.9); hC->Draw("E1 same");
            }
            TH1D* hx = S[cX17].Hist(Nm("hx").c_str(), 2); StyleH(hx, kColX17, 3, 2); hx->Draw("hist same");
            StyleH(hSum, kColSim, 4, 1); hSum->Draw("hist same");
            StyleH(hSumA, kColSim, 3, 2); hSumA->Draw("hist same");
            StyleM(hS, kColSig, 20, lin ? 1.1 : 0.9); hS->Draw("E1 same");
            VLines({scan::kBumpLo, scan::kBumpHi, scan::kPeakLo, scan::kPeakHi}, fr->GetMinimum(), fr->GetMaximum());
            fr->Draw("axis same");
            Txt(0.13, 0.92, Form("(%c) %s; %s", 'a' + 2 * j + v, kJName[j], lin ? "linear, 2.5 MeV, stacked; data N - CB only" : "log, 5 MeV"), 0.045);
            if (!lin) {
                Txt(0.15, 0.84, "data/sim per region (T1+T2 | variant (a)):", 0.036);
                const int qs[6] = {0, 1, 3, 5, 6, 7};
                const char* qq[6] = {"0.35-0.455", "bump", "peak", "0.56-0.70", "0.70-1.00", "0.35-1.00"};
                for (int q = 0; q < 6; ++q)
                    Txt(0.15 + 0.27 * (q % 3), 0.785 - 0.05 * (q / 3), Form("%s: %.2f#pm%.2f | %.2f", qq[q], rat[j][qs[q]], ratE[j][qs[q]], ratA[j][qs[q]]), 0.033);
                Txt(0.15, 0.68, Form("purity p p #eta (X), M_{4} < 1: %.3f #pm %.3f (sim.); bump window %.2f, peak window %.2f", pur[j], purE[j],
                                     purReg[j][1].eps, purReg[j][3].eps), 0.033, kColPeakTotal);
                Txt(0.15, 0.63, Form("#color[2]{above 0.56: multi-pion bkg %.0f%% (0.56-0.70), %.0f%% (0.70-1.00) of the sim: SMASH not judged there}",
                                     100 * sReg[j][cKM].s[5] / sReg[j][cSum].s[5], 100 * sReg[j][cKM].s[6] / sReg[j][cSum].s[6]), 0.031);
            } else {
                double pc = sReg[j][cPC].s[3] + sReg[j][cPD].s[3];
                Txt(0.15, 0.84, Form("peak window: sim. peak total %.3f nb of %.2f nb (T1+T2); data N - CB %.2f nb", pc / kLsim, sNb[j][3], dNb[j][3]), 0.034);
                Txt(0.15, 0.79, "X17 (BR 10^{-5}) dashed red, not stacked", 0.034, kColX17);
            }
            // ratio pad
            cv->cd(); MkPad(Nm("pr").c_str(), x1, y1, x2, ys, false, 0.12, 0.02, 0.0, 0.30);
            TH1D* r1 = RatioH(hS, hSum, Nm("r1").c_str());
            TH1D* r2 = RatioH(hS, hSumA, Nm("r2").c_str());
            TH1D* frr = Frame(Nm("frr").c_str(), S[cSum].nbin, xlo, xhi, kXM4, "data/sim", 0.3, 2.3, 0.12, 0.45);
            frr->GetYaxis()->SetNdivisions(505);
            frr->Draw("axis");
            HLine(xlo, xhi, 1., kBlack, 2, 1);
            StyleM(r1, kColSig, 20, 0.8); r1->Draw("E1 same");
            StyleM(r2, kColSim, 24, 0.8); r2->Draw("E1 same");
            for (int k = 0; k < seen::kNReg; ++k) {
                const double a = std::max(xlo, seen::kReg[k].lo), b = std::min(xhi, seen::kReg[k].hi);
                if (!(a < b)) continue;
                HLine(a, b, rat[j][k], kBlue + 3, 1, 5);
                HLine(a, b, ratA[j][k], kColSim, 2, 4);
            }
            VLines({scan::kBumpLo, scan::kBumpHi, scan::kPeakLo, scan::kPeakHi}, 0.3, 2.3);
            (void)qn;
        }
    }
    Save(cv, "eta_summary_b_datasim");
}

// ============================================================================
// Figure c: fit factors across selections (tables/eta_enriched_fit.txt) + check of E1/E2 M1 and P_main
// ============================================================================
void PartC() {
    TStopwatch sw;
    gFits = ParseFits();
    printf("PartC: %zu fit rows parsed\n", gFits.size());
    for (int is = 0; is < seen::kNSel; ++is) {
        std::string l;
        PMRow& p = gPM[is];
        if (!FindLine("eta_enriched_fit.txt", "peak", seen::kSel[is].tag, l)) continue;
        const std::vector<std::vector<VE>> g = Groups(l);
        if (g.size() < 3 || g[0].size() < 2 || g[1].size() < 5 || g[2].size() < 4) continue;
        p.data = g[0][0]; p.dataHand = g[0][1];
        p.sim54 = g[1][0]; p.sim62 = g[1][1]; p.simA = g[1][2]; p.pkOnly = g[1][3]; p.simUns = g[1][4];
        p.r54 = g[2][0]; p.r62 = g[2][1]; p.rA = g[2][2]; p.rUns = g[2][3];
        p.ok = true;
    }
    // ---- check: recompute M1 (0.40-0.75, no scale, T12 and (a)) and the data P_main of E1 and E2
    C("\n## figure c: M1 fits (0.40-0.75, no scale) and data P_main of E1 / E2 recomputed against tables/eta_enriched_fit.txt\n");
    const std::vector<Comp> fc = seen::FitComps();
    for (int j = 0; j < 2; ++j) {
        const int is = kJ[j];
        const DSpec D = DataSpectrum(*gDD, seen::DatSel(is), 70, 0.40, 0.75);
        std::vector<int> bins; for (int i = 0; i < D.nbin; ++i) bins.push_back(i);
        Binning b; b.nbin = 70; b.lo = 0.40; b.hi = 0.75; b.x = [](const Row& r) { return (double)r.m4; };
        const std::vector<WSpec> TS = FillMany(*gD, seen::SimSel(is), fc, gOpt, b);
        for (int gv = 0; gv < 2; ++gv) {
            const int grp[3] = {seen::fBump, seen::fPeak, gv ? seen::fKA : seen::fK};
            std::vector<std::vector<double>> TT, TV;
            for (int k : grp) {
                std::vector<double> t(b.nbin), tv(b.nbin);
                for (int i = 0; i < b.nbin; ++i) { t[i] = TS[k].s[i] * kToRows; tv[i] = TS[k].v[i] * kToRows * kToRows; }
                TT.push_back(t); TV.push_back(tv);
            }
            const seen::FitRes r = seen::TemplateFit(D.sp, bins, TT, TV);
            const FitRow* f = FindFit(gFits, kJTag[j], 0, false, gv == 1, "M1");
            const char* pn[3] = {"a_bump", "a_peak", "a_bkg"};
            for (int k = 0; k < 3; ++k)
                Check(Form("%s M1 %s %s", kJTag[j], gv ? "(a)" : "T12", pn[k]), r.a[k], f ? f->Par(pn[k]).first : kNaN, 6e-4);
            Check(Form("%s M1 %s chi2", kJTag[j], gv ? "(a)" : "T12"), r.chi2, f ? f->chi2 : kNaN, 0.06);
        }
        const scan::PeakFit pm = scan::FitPeakMain(scan::FillMany(*gDD, {seen::DatSel(is)}, scan::FillOpt())[0]);
        Check(Form("%s data P_main [rows]", kJTag[j]), pm.S, gPM[is].data.first, 0.6);
    }
    gHaveC = true;
    printf("PartC: %.1f s\n", sw.RealTime());
}
void FigC() {
    using namespace seen;
    TCanvas* cv = new TCanvas(Nm("cC").c_str(), "", 2700, 2050);
    cv->cd();
    MkPad(Nm("hC").c_str(), 0, 0.855, 1, 1, false, 0, 0, 0, 0);
    Txt(0.01, 0.87, "Template fits of data N - CB with free factors #alpha on the SMASH components, per #eta-enriched selection (a comparison for discussion, not a measurement)", 0.13);
    Txt(0.01, 0.75, "#alpha = 1: the SMASH prediction for L_{data} = 6080.5 nb^{-1} (UNCONFIRMED; every #alpha scales with 6080.5 / L_{data}).  "
                    "M1: bump (#pi^{0} Dalitz + conversion), peak (direct + conversion, SMASH ratio), physical background.  M3: background split into #eta-related and multi-pion", 0.06);
    Txt(0.01, 0.65, "fits: 5 MeV bins, data covariance (N + interpolated CB) + #alpha^{2} #times template variance; numbers from tables/eta_enriched_fit.txt (sim_eta_enriched.C), E1 / E2 M1 recomputed here. "
                    "P_{main}: data-study peak fit (Gaussian 0.5439 / 5.4 MeV, pol2) on data and on the simulation (M_{4} #times 0.9952)", 0.06);
    Txt(0.01, 0.55, "selections: E1 slides, E2 slides without OA, E3 old code (M_{#pi#pi} < 0.42), E4 #eta band (no M_{#pi#pi} cut), 'oa' + OA > 2#circ, 'mm' + MM_{4} < 2.85, 'v' + both pions in the vertex fit; all pre.  "
                    "#color[2]{Poor fits (#chi^{2}/ndf up to 3.3): error bars scaled by #sqrt{#chi^{2}/ndf}; the factors are NOT yields}", 0.06);
    // variant markers: 4 range/scale combinations x ghost variant (filled T1+T2, open (a))
    const int vcol[4] = {kAzure + 2, kOrange + 7, kGreen + 2, kMagenta + 1};
    const int vmk[4][2] = {{20, 24}, {21, 25}, {22, 26}, {23, 32}};
    const char* vname[4] = {"0.40-0.75", "0.43-0.70", "0.40-0.75, M_{4}^{sim} #times 0.9952 (K0S)", "0.43-0.70, #times 0.9952"};
    TLegend* lg = new TLegend(0.0, 0.0, 1.0, 0.45);
    lg->SetNColumns(5); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.062);
    for (int c = 0; c < 4; ++c) {
        TH1D* h = Dummy(Nm("lg").c_str(), vcol[c], 2, 1, 0, vmk[c][0], vcol[c]);
        lg->AddEntry(h, Form("fit range %s", vname[c]), "p");
    }
    lg->AddEntry(Dummy(Nm("lg").c_str(), kGray + 2, 2, 1, 0, 24, kGray + 2), "open markers: variant (a), + ghost-flagged true pairs", "p");
    lg->AddEntry(Dummy(Nm("lg").c_str(), kBlack, 2, 1, 0, 29, kBlack), "P_{main} data / sim (5.4 MeV)", "p");
    lg->AddEntry(Dummy(Nm("lg").c_str(), kBlack, 2, 1, 0, 30, kBlack), "P_{main} data / sim, variant (a)", "p");
    lg->AddEntry(Dummy(Nm("lg").c_str(), kColEtaBkg, 2, 1, 0, 20, kColEtaBkg), "M3: #alpha_{bkg}, #eta-related", "p");
    lg->AddEntry(Dummy(Nm("lg").c_str(), kColPhysBkg, 2, 1, 0, 21, kColPhysBkg), "M3: #alpha_{bkg}, multi-pion", "p");
    lg->Draw();

    std::vector<std::string> labs;
    for (int is = 0; is < seen::kNSel; ++is) labs.push_back(seen::kSel[is].tag);
    auto offs = [](int c, int a) { return -0.38 + 0.085 * (2 * c + a); };
    struct PanelDef { const char* title; const char* yt; const char* model; const char* par; double ylo, yhi; };
    const PanelDef pd[4] = {{"(a) #alpha_{bump}, model M1", "#alpha_{bump}", "M1", "a_bump", 0.0, 0.8},
                            {"(b) #alpha_{peak} (M1); stars: P_{main} data/sim", "#alpha_{peak}, data/sim", "M1", "a_peak", 0.0, 1.9},
                            {"(c) #alpha_{bkg}, model M1", "#alpha_{bkg}", "M1", "a_bkg", 0.9, 1.5},
                            {"(d) #alpha_{bump}, model M3 (background split)", "#alpha_{bump}", "M3", "a_bump", 0.0, 0.8}};
    for (int p = 0; p < 6; ++p) {
        const double x1 = (p % 3) / 3.0, x2 = (p % 3 + 1) / 3.0, y1 = p < 3 ? 0.4275 : 0.0, y2 = p < 3 ? 0.855 : 0.4275;
        cv->cd(); MkPad(Nm("pc").c_str(), x1, y1, x2, y2, false, 0.13, 0.02, 0.10, 0.10);
        const char* yt = p < 4 ? pd[p].yt : (p == 4 ? "#alpha_{bkg} (M3)" : "#chi^{2}/ndf");
        const double ylo = p < 4 ? pd[p].ylo : (p == 4 ? 0.4 : 0.0), yhi = p < 4 ? pd[p].yhi : (p == 4 ? 2.95 : 4.5);
        TH1D* fr = CatFrame(Nm("fr").c_str(), labs, yt, ylo, yhi, 0.05, 0.045, 1.2);
        fr->Draw("axis");
        if (ylo < 1 && yhi > 1) HLine(-0.5, seen::kNSel - 0.5, 1., kBlack, 2, 2);
        for (int is = 1; is < seen::kNSel; ++is) { TLine* l = new TLine(is - 0.5, ylo, is - 0.5, yhi); l->SetLineColor(kGray); l->SetLineStyle(3); l->Draw(); }
        double vmin = 1e9, vmax = -1e9;
        if (p < 4 || p == 5) {
            for (int c = 0; c < 4; ++c)
                for (int a = 0; a < 2; ++a) {
                    std::vector<double> x; std::vector<VE> y;
                    for (int is = 0; is < seen::kNSel; ++is) {
                        const FitRow* f = FindFit(gFits, seen::kSel[is].tag, c % 2, c >= 2, a == 1, p == 5 ? "M1" : pd[p].model);
                        if (!f) continue;
                        VE v = p == 5 ? VE(f->chi2 / f->ndf, 0.) : f->Par(pd[p].par);
                        if (p != 5) v.second *= std::sqrt(std::max(1., f->chi2 / f->ndf));   // review finding 5: scaled errors
                        x.push_back(is + offs(c, a)); y.push_back(v);
                        if (v.first == v.first) { vmin = std::min(vmin, v.first); vmax = std::max(vmax, v.first); }
                    }
                    TGraphErrors* g = Pts(x, y, vcol[c], vmk[c][a], 1.5);
                    g->Draw("P same");
                }
        }
        if (p == 1) {   // P_main data/sim
            std::vector<double> x0, x1v; std::vector<VE> y0, y1v;
            for (int is = 0; is < seen::kNSel; ++is) {
                if (!gPM[is].ok) continue;
                x0.push_back(is + 0.33); y0.push_back(gPM[is].r54);
                x1v.push_back(is + 0.43); y1v.push_back(gPM[is].rA);
            }
            Pts(x0, y0, kBlack, 29, 2.2)->Draw("P same");
            Pts(x1v, y1v, kBlack, 30, 2.2)->Draw("P same");
        }
        if (p == 4) {   // M3 background factors, 0.40-0.75, no scale
            const char* pn[2] = {"a_bkgEta", "a_bkgMP"};
            const int pc[2] = {kColEtaBkg, kColPhysBkg}, pm[2][2] = {{20, 24}, {21, 25}};
            for (int q = 0; q < 2; ++q)
                for (int a = 0; a < 2; ++a) {
                    std::vector<double> x; std::vector<VE> y;
                    for (int is = 0; is < seen::kNSel; ++is) {
                        const FitRow* f = FindFit(gFits, seen::kSel[is].tag, 0, false, a == 1, "M3");
                        if (!f) continue;
                        VE v = f->Par(pn[q]); v.second *= std::sqrt(std::max(1., f->chi2 / f->ndf));
                        x.push_back(is - 0.27 + 0.18 * (2 * q + a)); y.push_back(v);
                    }
                    Pts(x, y, pc[q], pm[q][a], 1.6)->Draw("P same");
                }
        }
        fr->Draw("axis same");
        const char* title = p < 4 ? pd[p].title : (p == 4 ? "(e) #alpha_{bkg} of M3: #eta-related, multi-pion" : "(f) #chi^{2}/ndf of M1 (ndf 67 or 51)");
        Txt(0.14, 0.92, title, 0.047);
        if (p < 4 && vmax > vmin) Txt(0.16, 0.84, Form("all selections and variants: %.2f - %.2f", vmin, vmax), 0.04);
        if (p == 0 || p == 3) Txt(0.16, 0.79, "poor fits, not a yield: sim bump too large vs sim bkg at the edge", 0.032);
        if (p == 4) Txt(0.16, 0.79, "#eta-related bkg mostly > 1 while #alpha_{bump} < 1: contradictory", 0.034);
        if (p == 1) Txt(0.16, 0.79, "stars: P_{main} on data / on the simulation (sim M_{4} #times 0.9952)", 0.034);
        if (p == 4) Txt(0.16, 0.84, "fit range 0.40-0.75, no scale; filled T1+T2, open (a)", 0.04);
        if (p == 5) Txt(0.16, 0.84, Form("range of M1 #chi^{2}/ndf: %.2f - %.2f", vmin, vmax), 0.04);
    }
    Save(cv, "eta_summary_c_alpha");
}

// ============================================================================
// Figure d: peak content per OA and M(ee) class, data against simulation
// ============================================================================
void PartD() {
    TStopwatch sw;
    const SimData& ds = *gDS;
    const sep::KWTab kw = sep::MakeKW(ds);
    const std::vector<Comp> cc = {CompPeakDirect(), CompPeakConv(), CompPeakTotal(), sep::CompDirectKW(kw), CompX17(1e-5)};
    struct Job { std::string sel, cs; Sel ssel; scan::Sel dsel; };
    std::vector<Job> jobs;
    for (int s = 0; s < 2; ++s)
        for (const char* csn : {"OA5", "MEE5"}) {
            Job jb; jb.sel = s == 0 ? "pre" : "E2"; jb.cs = csn;
            jb.ssel = s == 0 ? SelPreset(kPre) : seen::SimSel(seen::sE2);
            jb.dsel = s == 0 ? DataSelPreset(kPre) : seen::DatSel(seen::sE2);
            jobs.push_back(jb);
        }
    // data fits: per class, and per selection without class
    std::vector<scan::Sel> dsels;
    std::vector<std::vector<int>> idx(jobs.size());
    for (size_t k = 0; k < jobs.size(); ++k) {
        const sep::ClassSet cs = sep::FindCS(jobs[k].cs);
        for (size_t i = 0; i < cs.lo.size(); ++i) {
            const double lo = cs.lo[i], hi = cs.hi[i]; const int var = cs.var; const scan::Sel b = jobs[k].dsel;
            idx[k].push_back((int)dsels.size());
            dsels.push_back([b, lo, hi, var](const scan::Row& r) { if (!b(r)) return false; const double x = sep::CVar(var, r); return x >= lo && x < hi; });
        }
    }
    const int iTot0 = (int)dsels.size();
    dsels.push_back(DataSelPreset(kPre));
    dsels.push_back(seen::DatSel(seen::sE2));
    const std::vector<sep::DFit> df = sep::FitDataSels(*gDD, dsels);
    for (size_t k = 0; k < jobs.size(); ++k) {
        const sep::ClassSet cs = sep::FindCS(jobs[k].cs);
        ClassRes o; o.sel = jobs[k].sel; o.cs = jobs[k].cs; o.name = cs.name;
        for (size_t i = 0; i < cs.lo.size(); ++i) {
            const double lo = cs.lo[i], hi = cs.hi[i]; const int var = cs.var;
            const std::vector<WSum> s = SumMany(ds, SelAnd(jobs[k].ssel, [lo, hi, var](const Row& r) { const double x = sep::CVar(var, r); return x >= lo && x < hi; }), cc, gOpt);
            o.sd.push_back(s[0].S * kToRows); o.sde.push_back(std::sqrt(s[0].V) * kToRows);
            o.sc.push_back(s[1].S * kToRows); o.sce.push_back(std::sqrt(s[1].V) * kToRows);
            o.st.push_back(s[2].S * kToRows); o.ste.push_back(std::sqrt(s[2].V) * kToRows);
            o.sk.push_back(s[3].S * kToRows); o.ske.push_back(std::sqrt(s[3].V) * kToRows);
            o.sx.push_back(s[4].S * kToRows);
            o.dS.push_back(df[idx[k][i]].S); o.dSe.push_back(df[idx[k][i]].eT());
            o.dSs.push_back(df[idx[k][i]].eS); o.dSm.push_back(df[idx[k][i]].eM);
        }
        o.mix = sep::MixFit(o.dS, o.dSe, o.sc, o.sce, o.sd, o.sde);
        o.mixKW = sep::MixFit(o.dS, o.dSe, o.sc, o.sce, o.sk, o.ske);
        o.mixStat = sep::MixFit(o.dS, o.dSs, o.sc, o.sce, o.sd, o.sde);
        const sep::DFit& ft = df[iTot0 + (o.sel == "pre" ? 0 : 1)];
        o.dTot = ft.S; o.dTotE = ft.eT();
        gCls.push_back(o);
    }
    // checks against tables/eta_peak_classes.txt (pre)
    C("\n## figure d: peak classes (pre) against tables/eta_peak_classes.txt\n");
    for (const ClassRes& o : gCls) {
        if (o.sel != "pre") continue;
        const std::string sec = "pre, " + o.cs;
        for (size_t i = 0; i < o.name.size(); ++i) {
            Check(Form("pre %s %s: sim total [rows]", o.cs.c_str(), o.name[i].c_str()), o.st[i], GetVE("eta_peak_classes.txt", sec, o.name[i], 1, 0).first, 0.6);
            Check(Form("pre %s %s: data S [rows]", o.cs.c_str(), o.name[i].c_str()), o.dS[i], GetVE("eta_peak_classes.txt", sec, o.name[i], 4, 0).first, 0.6);
        }
        Check(Form("pre %s: one-factor a1", o.cs.c_str()), o.mix.a1, GetVE("eta_peak_classes.txt", sec, "mixture fit:", -1, 2).first, 6e-4);
        Check(Form("pre %s: mixture a(conv)", o.cs.c_str()), o.mix.a, GetVE("eta_peak_classes.txt", sec, "mixture fit:", -1, 0).first, 6e-4);
        Check(Form("pre %s: KW mixture b(direct) (illustration)", o.cs.c_str()), o.mixKW.b, GetVE("eta_peak_classes.txt", sec, "ILLUSTRATION ONLY", -1, 1).first, 6e-4);
    }
    gHaveD = true;
    printf("PartD: %.1f s\n", sw.RealTime());
}
void FigD() {
    TCanvas* cv = new TCanvas(Nm("cD").c_str(), "", 2800, 2050);
    cv->cd();
    MkPad(Nm("hD").c_str(), 0, 0.855, 1, 1, false, 0, 0, 0, 0);
    Txt(0.01, 0.87, "#eta peak content per class of OA(e^{+}e^{-}) (top) and M_{ee} (bottom): simulation (direct decay + conversion) against data, selections pre and E2", 0.13);
    Txt(0.01, 0.75, "simulation: std + lep + eta owned rows, w (PDG), all M_{4}, expected data rows = #Sigma w #times L_{data} / L_{sim}; data: peak content = P_{main} fit of N - CB per class "
                    "(Gaussian frozen 0.5439 / 5.4 MeV, pol2, 0.5225-0.61); L_{data} = 6080.5 nb^{-1} UNCONFIRMED", 0.062);
    Txt(0.01, 0.65, "#color[2]{sim direct decay: generator M_{ee} model NOT QED (accepted density largest at threshold); its class shapes are not for physics.} "
                    "Dashed: QED-reshaped direct = ILLUSTRATION of shapes only (flat acceptance and fixed content assumed), no fit result quoted", 0.062);
    Txt(0.01, 0.55, "pre = isBest #geq 0, PT3, v_{z} > -160 mm, corrected momenta; E2 = pre + #alpha_{lab} < 50#circ, #alpha_{#eta} > 140#circ, M_{#pi#pi} < 0.46 (no OA cut); "
                    "data errors: stat #oplus fit-model spread per class (48 variants); fraction errors include the correlation with the total", 0.062);
    TLegend* lg = new TLegend(0.0, 0.0, 1.0, 0.45);
    lg->SetNColumns(4); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.064);
    lg->AddEntry(Dummy(Nm("lg").c_str(), kColPeakDirect, 3, 1, 0), "sim direct #eta #rightarrow #pi^{+}#pi^{-}e^{+}e^{-}, gen. M_{ee} not QED (shape)", "l");
    lg->AddEntry(Dummy(Nm("lg").c_str(), kColPeakConv, 3, 1, 0), "sim conversion #eta #rightarrow #pi^{+}#pi^{-}#gamma (shape)", "l");
    lg->AddEntry(Dummy(Nm("lg").c_str(), kColPeakTotal, 4, 1, 0), "sim total (shape)", "l");
    lg->AddEntry(Dummy(Nm("lg").c_str(), kColPeakDirect, 3, 2, 0), "direct, QED M_{ee}: illustration only", "l");
    lg->AddEntry(Dummy(Nm("lg").c_str(), kColPeakDirect, 1, 1, kColPeakDirect), "sim direct (absolute, on top)", "f");
    lg->AddEntry(Dummy(Nm("lg").c_str(), kColPeakConv, 1, 1, kColPeakConv), "sim conversion (absolute)", "f");
    lg->AddEntry(Dummy(Nm("lg").c_str(), kColData, 2, 1, 0, 20, kColData), "data: P_{main} fit of N - CB per class (stat #oplus model)", "p");
    lg->AddEntry(Dummy(Nm("lg").c_str(), kColX17, 3, 2, 0), "X17, BR 10^{-5} (information, not in the sum)", "l");
    lg->Draw();
    const char* xt[2] = {"OA(e^{+}e^{-}) class [deg]", "M_{ee} class [MeV/c^{2}]"};
    for (int row = 0; row < 2; ++row) {
        const double y1 = row == 0 ? 0.4275 : 0.0, y2 = row == 0 ? 0.855 : 0.4275;
        for (int s = 0; s < 2; ++s) {
            const ClassRes* o = nullptr;
            for (const ClassRes& x : gCls) if (x.sel == (s == 0 ? "pre" : "E2") && x.cs == (row == 0 ? "OA5" : "MEE5")) o = &x;
            if (!o) continue;
            const int n = (int)o->name.size();
            double tD = 0, tS = 0, tC = 0, tT = 0, tK = 0;
            for (int i = 0; i < n; ++i) { tD += o->dS[i]; tS += o->sd[i]; tC += o->sc[i]; tT += o->st[i]; tK += o->sk[i]; }
            // shapes
            const double xs1 = s * 0.5, xs2 = s * 0.5 + 0.22;
            cv->cd(); MkPad(Nm("ps").c_str(), xs1, y1, xs2, y2, false, 0.17, 0.02, 0.10, 0.13);
            TH1D* fr = CatFrame(Nm("fr").c_str(), o->name, "fraction of the peak content", 0, 1.0, 0.055, 0.075, 1.35);
            fr->GetXaxis()->SetTitle(xt[row]); fr->GetXaxis()->SetTitleSize(0.055); fr->GetXaxis()->SetTitleOffset(1.0);
            fr->Draw("axis");
            TH1D* h[4];
            const std::vector<double>* src[4] = {&o->sd, &o->sc, &o->st, &o->sk};
            const double tt[4] = {tS, tC, tT, tK};
            const int hc[4] = {kColPeakDirect, kColPeakConv, kColPeakTotal, kColPeakDirect}, hl[4] = {3, 3, 4, 3}, hs[4] = {1, 1, 1, 2};
            for (int q = 0; q < 4; ++q) {
                h[q] = new TH1D(Nm("hf").c_str(), "", n, -0.5, n - 0.5); h[q]->SetDirectory(nullptr);
                for (int i = 0; i < n; ++i) h[q]->SetBinContent(i + 1, (*src[q])[i] / tt[q]);
                StyleH(h[q], hc[q], hl[q], hs[q]); h[q]->Draw("hist same");
            }
            std::vector<double> xv; std::vector<VE> yv;
            for (int i = 0; i < n; ++i) { xv.push_back(i); yv.push_back(Frac(o->dS, o->dSe, i)); }
            Pts(xv, yv, kColData, 20, 1.8, 3)->Draw("P same");
            fr->Draw("axis same");
            Txt(0.18, 0.92, Form("(%c) %s: shapes (each normalised to 1)", 'a' + 4 * row + 2 * s, o->sel.c_str()), 0.05);
            Txt(0.20, 0.84, Form("sim direct share %.2f", tS / tT), 0.045, kColPeakTotal);
            // absolute + ratio
            const double xa1 = s * 0.5 + 0.22, xa2 = s * 0.5 + 0.5, ysp = y1 + 0.32 * (y2 - y1);
            cv->cd(); MkPad(Nm("pab").c_str(), xa1, ysp, xa2, y2, false, 0.16, 0.03, 0.10, 0.0);
            double ymax = 0;
            for (int i = 0; i < n; ++i) ymax = std::max(ymax, std::max(o->st[i] + o->ste[i], o->dS[i] + o->dSe[i]));
            TH1D* fa = CatFrame(Nm("fa").c_str(), o->name, "peak content [rows at L_{data}]", 0, ymax * 1.75, 0.06, 0.0, 1.25);
            fa->Draw("axis");
            TH1D* hcv = new TH1D(Nm("hc").c_str(), "", n, -0.5, n - 0.5); hcv->SetDirectory(nullptr);
            TH1D* htt = new TH1D(Nm("ht").c_str(), "", n, -0.5, n - 0.5); htt->SetDirectory(nullptr);
            TH1D* hx = new TH1D(Nm("hx").c_str(), "", n, -0.5, n - 0.5); hx->SetDirectory(nullptr);
            for (int i = 0; i < n; ++i) {
                hcv->SetBinContent(i + 1, o->sc[i]); htt->SetBinContent(i + 1, o->st[i]); htt->SetBinError(i + 1, o->ste[i]);
                hx->SetBinContent(i + 1, o->sx[i]);
            }
            StyleH(htt, kColPeakDirect, 1, 1, kColPeakDirect); htt->Draw("hist same");
            StyleH(hcv, kColPeakConv, 1, 1, kColPeakConv); hcv->Draw("hist same");
            TH1D* hte = (TH1D*)htt->Clone(Nm("hte").c_str()); hte->SetFillStyle(0); hte->SetLineColor(kBlack); hte->SetLineWidth(1); hte->Draw("E1 same");
            StyleH(hx, kColX17, 3, 2); hx->Draw("hist same");
            std::vector<VE> ya; for (int i = 0; i < n; ++i) ya.push_back(VE(o->dS[i], o->dSe[i]));
            Pts(xv, ya, kColData, 20, 1.8, 3)->Draw("P same");
            fa->Draw("axis same");
            Txt(0.17, 0.92, Form("(%c) %s: absolute = method A (fig. f), L_{data} UNCONF.", 'a' + 4 * row + 2 * s + 1, o->sel.c_str()), 0.045);
            Txt(0.19, 0.84, Form("sum: data %.0f#pm%.0f, sim %.0f#pm%.0f rows", tD,
                                 [&]() { double v = 0; for (double e : o->dSe) v += e * e; return std::sqrt(v); }(), tT,
                                 [&]() { double v = 0; for (double e : o->ste) v += e * e; return std::sqrt(v); }()), 0.037);
            Txt(0.19, 0.785, Form("one factor %.3f#pm%.3f, #chi^{2}/ndf %.1f/%d (stat only: %.0f/%d)", o->mix.a1, o->mix.ea1, o->mix.chi21, o->mix.ndf1,
                                  o->mixStat.chi21, o->mixStat.ndf1), 0.037, kColPeakTotal);
            Txt(0.19, 0.73, Form("free: conv. %.2f#pm%.2f, direct %.2f#pm%.2f (%.1f/%d)", o->mix.a, o->mix.ea, o->mix.b, o->mix.eb,
                                 o->mix.chi2, o->mix.ndf), 0.037);
            if (o->mix.a < 0) Txt(0.19, 0.675, "   negative conversion factor: unphysical", 0.037);
            Txt(0.19, 0.62, "data: more peak at OA > 10#circ, M_{ee} > 20 MeV than either sim component", 0.034);
            cv->cd(); MkPad(Nm("prt").c_str(), xa1, y1, xa2, ysp, false, 0.16, 0.03, 0.0, 0.30);
            double rmax = 0;
            std::vector<VE> rr;
            for (int i = 0; i < n; ++i) {
                const double r = o->dS[i] / o->st[i];
                const double e = std::fabs(r) * std::sqrt(std::pow(o->dSe[i] / o->dS[i], 2) + std::pow(o->ste[i] / o->st[i], 2));
                rr.push_back(VE(r, e)); rmax = std::max(rmax, r + e);
            }
            TH1D* fq = CatFrame(Nm("fq").c_str(), o->name, "data / sim", 0, std::min(6.0, std::max(1.6, rmax * 1.15)), 0.12, 0.14, 0.55);
            fq->GetXaxis()->SetTitle(xt[row]); fq->GetXaxis()->SetTitleSize(0.12); fq->GetXaxis()->SetTitleOffset(1.0);
            fq->Draw("axis");
            HLine(-0.5, n - 0.5, 1., kBlack, 2, 1);
            HLine(-0.5, n - 0.5, o->mix.a1, kColPeakTotal, 2, 3);
            Pts(xv, rr, kColData, 20, 1.6, 3)->Draw("P same");
        }
    }
    Save(cv, "eta_summary_d_classes");
}

// ============================================================================
// Figure e: cut efficiencies
// ============================================================================
void PartE() {
    TStopwatch sw;
    const SimData& ds = *gDS;
    const sep::KWTab kw = sep::MakeKW(ds);
    std::vector<Comp> cc(eNE);
    cc[eDir] = CompPeakDirect(); cc[eConv] = CompPeakConv(); cc[eTot] = CompPeakTotal(); cc[eTotKW] = sep::CompTotalKW(kw);
    cc[eBump] = Comp("bump D+C", [](const Row& r) { return r.tmpl == 1 && (r.cls == 2 || r.cls == 3); });
    cc[eTag] = Comp("bump D+C, photon-tagged", [](const Row& r) {
        return r.tmpl == 1 && (r.cls == 2 || r.cls == 3) && r.g_idx > 0 && r.m_eeg >= 0.119f && r.m_eeg < 0.173f;
    });
    cc[eSbu] = Comp("S_bump", [](const Row& r) { return r.tmpl == 1 && (r.cls == 2 || r.cls == 3) && r.m4 >= scan::kBumpLo && r.m4 < scan::kBumpHi; });
    const sep::Cut base = sep::CPre();
    sep::Cut bstart; bstart.vzLo = -1e9; bstart.labHi = 40; bstart.mppLo = 0.30; bstart.mppHi = 0.41;   // sim_cuts.C bump start
    std::vector<scan::Sel> dsel;
    for (int k = 0; k < kNCut; ++k) {
        sep::Cut cp = base; sep::SetCut(cp, kCut[k].kind, kCut[k].val);
        const Sel bs = sep::SimSel(base), ps = sep::SimSel(cp);
        for (int c = 0; c < eNE; ++c) gEff[k].s[c] = Efficiency(ds, bs, ps, cc[c], gOpt);
        if (kCut[k].bumpStart) {
            sep::Cut cb = bstart; sep::SetCut(cb, kCut[k].kind, kCut[k].val);
            gEff[k].sBumpStart = Efficiency(ds, sep::SimSel(bstart), sep::SimSel(cb), cc[eBump], gOpt);
            gEff[k].sTagStart = Efficiency(ds, sep::SimSel(bstart), sep::SimSel(cb), cc[eTag], gOpt);
        }
        dsel.push_back(sep::DatSel(cp));
        dsel.push_back([base, cp](const scan::Row& r) { return sep::PassBase(r) && sep::PassCut(r, base) && !sep::PassCut(r, cp); });
    }
    const std::vector<sep::DFit> df = sep::FitDataSels(*gDD, dsel);
    for (int k = 0; k < kNCut; ++k) {
        const sep::Eff e = sep::PassFail(df[2 * k], df[2 * k + 1]);
        EffOut& o = gEff[k];
        o.dPeak = e.ok ? VE(e.e, e.err) : VE(kNaN, kNaN);
        o.dPeakTab = GetVE("eta_peak_eff.txt", "", kCut[k].keyPeak, 3, 1);
        o.tDir = GetVE("eta_peak_eff.txt", "", kCut[k].keyPeak, 1, 0);
        o.tConv = GetVE("eta_peak_eff.txt", "", kCut[k].keyPeak, 1, 1);
        o.tTot = GetVE("eta_peak_eff.txt", "", kCut[k].keyPeak, 1, 2);
        o.tTotKW = GetVE("eta_peak_eff.txt", "", kCut[k].keyPeak, 2, 1);
        o.dTag = GetVE("cuts_compare_data.txt", "", kCut[k].keyBump, 2, 0);
        o.tTag = GetVE("cuts_compare_data.txt", "", kCut[k].keyBump, 2, 1);
        o.dEdge = GetVE("cuts_compare_data.txt", "", kCut[k].keyBump, 3, 0);
        o.tSbu = GetVE("cuts_compare_data.txt", "", kCut[k].keyBump, 3, 1);
        o.tBump = GetVE("cuts_compare_data.txt", "", kCut[k].keyBump, 4, 0);
        o.tTagDal = GetVE("cuts_compare_data.txt", "", kCut[k].keyBump, 5, 0);
        if (o.dTag.first == o.dTag.first && !(o.dTag.second > 0)) { o.degTag = true; o.dTag = VE(kNaN, kNaN); }
        if (o.dEdge.first == o.dEdge.first && !(o.dEdge.second > 0)) { o.degEdge = true; o.dEdge = VE(kNaN, kNaN); }
    }
    C("\n## figure e: efficiencies against tables/eta_peak_eff.txt (peak, base pre) and tables/cuts_compare_data.txt (bump; sim at pre, or at the bump start for pT4 / vertex)\n");
    for (int k = 0; k < kNCut; ++k) {
        const EffOut& o = gEff[k];
        Check(Form("%s: sim direct", kCut[k].txt), o.s[eDir].eps, o.tDir.first, 6e-4);
        Check(Form("%s: sim conversion", kCut[k].txt), o.s[eConv].eps, o.tConv.first, 6e-4);
        Check(Form("%s: sim total", kCut[k].txt), o.s[eTot].eps, o.tTot.first, 6e-4);
        Check(Form("%s: sim total, QED direct", kCut[k].txt), o.s[eTotKW].eps, o.tTotKW.first, 6e-4);
        Check(Form("%s: data peak (pass/fail fits, pre)", kCut[k].txt), o.dPeak.first, o.dPeakTab.first, 6e-4);
        const double stag = kCut[k].bumpStart ? o.sTagStart.eps : o.s[eTag].eps;
        const double sbum = kCut[k].bumpStart ? o.sBumpStart.eps : o.s[eBump].eps;
        Check(Form("%s: sim bump tagged", kCut[k].txt), stag, o.tTag.first, 6e-4);
        Check(Form("%s: sim bump all", kCut[k].txt), sbum, o.tBump.first, 6e-4);
        if (!kCut[k].bumpStart) Check(Form("%s: sim S_bump", kCut[k].txt), o.s[eSbu].eps, o.tSbu.first, 6e-4);
    }
    gHaveE = true;
    printf("PartE: %.1f s\n", sw.RealTime());
}
void FigE() {
    TCanvas* cv = new TCanvas(Nm("cE").c_str(), "", 2800, 2100);
    cv->cd();
    MkPad(Nm("hE").c_str(), 0, 0.85, 1, 1, false, 0, 0, 0, 0);
    Txt(0.01, 0.87, "Efficiency of single cuts for the simulated #eta peak (direct decay, conversion) and #eta bump, next to the data-derived efficiencies", 0.13);
    Txt(0.01, 0.75, "simulation: std + lep + eta owned rows, w (PDG), all M_{4}, base pre (isBest #geq 0, PT3, v_{z} > -160 mm, corrected momenta); #varepsilon = #Sigma w(pass) / #Sigma w(base), event-group errors. "
                    "Data peak: pass/fail P_{main} fits of N - CB at the same base", 0.058);
    Txt(0.01, 0.65, "data bump (tables/cuts_compare_data.txt, data study): ECAL-tagged #eta (M_{ee#gamma} in the #pi^{0} window) and the edge-step sum E, base without the v_{z} cut; "
                    "p_{T,4} and vertex: bump start base (#alpha_{lab} < 40#circ, 0.30 < M_{#pi#pi} < 0.41) in data and simulation", 0.058);
    Txt(0.01, 0.55, "grey bands: data peak efficiency biased by sculpting (data study); data peak errors: stat #oplus fit-model spread (efficiency per model variant); sim direct: generator M_{ee} model NOT QED; "
                    "open circles: QED-reshaped direct = illustration only (flat acceptance, fixed content assumed), no pulls; edge step at M_{#pi#pi} < 0.40 degenerate, not shown", 0.058);
    TLegend* lg = new TLegend(0.0, 0.0, 1.0, 0.46);
    lg->SetNColumns(4); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.062);
    lg->AddEntry(Dummy(Nm("lg").c_str(), kColPeakDirect, 2, 1, 0, 22, kColPeakDirect), "sim peak: direct (gen. M_{ee} not QED)", "p");
    lg->AddEntry(Dummy(Nm("lg").c_str(), kColPeakConv, 2, 1, 0, 21, kColPeakConv), "sim peak: conversion #pi^{+}#pi^{-}#gamma", "p");
    lg->AddEntry(Dummy(Nm("lg").c_str(), kColPeakTotal, 2, 1, 0, 20, kColPeakTotal), "sim peak: total", "p");
    lg->AddEntry(Dummy(Nm("lg").c_str(), kColPeakTotal, 2, 1, 0, 24, kColPeakTotal), "total, QED-M_{ee} direct (illustration)", "p");
    lg->AddEntry(Dummy(Nm("lg").c_str(), kColData, 2, 1, 0, 29, kColData), "data peak: pass/fail fits (stat #oplus model)", "p");
    lg->AddEntry(Dummy(Nm("lg").c_str(), kColBumpDalitz, 2, 1, 0, 20, kColBumpDalitz), "sim bump (#pi^{0} Dalitz + conversion)", "p");
    lg->AddEntry(Dummy(Nm("lg").c_str(), kOrange + 3, 2, 1, 0, 21, kOrange + 3), "sim bump, photon-tagged", "p");
    lg->AddEntry(Dummy(Nm("lg").c_str(), kOrange + 3, 2, 1, 0, 25, kOrange + 3), "sim tagged, #eta#rightarrow3#pi Dalitz-plot weight (UNVERIFIED)", "p");
    lg->AddEntry(Dummy(Nm("lg").c_str(), kColBumpConv, 2, 1, 0, 22, kColBumpConv), "sim bump in the bump window (S_{bump})", "p");
    lg->AddEntry(Dummy(Nm("lg").c_str(), kColData, 2, 1, 0, 21, kColData), "data: ECAL-tagged #eta", "p");
    lg->AddEntry(Dummy(Nm("lg").c_str(), kColData, 2, 1, 0, 26, kColData), "data: bump edge step E", "p");
    lg->Draw();

    std::vector<std::string> labs;
    for (int k = 0; k < kNCut; ++k) labs.push_back(kCut[k].lab);
    // (a) peak
    cv->cd(); MkPad(Nm("pe1").c_str(), 0, 0.49, 1, 0.85, false, 0.06, 0.01, 0.09, 0.10);
    TH1D* f1 = CatFrame(Nm("f").c_str(), labs, "efficiency", 0.4, 1.08, 0.06, 0.062, 0.45);
    f1->Draw("axis");
    for (int k = 0; k < kNCut; ++k)
        if (kCut[k].sculpt) { TBox* b = new TBox(k - 0.48, 0.4, k + 0.48, 1.08); b->SetFillColor(kGray); b->SetFillStyle(1001); b->Draw(); }
    for (int k = 1; k < kNCut; ++k) { TLine* l = new TLine(k - 0.5, 0.4, k - 0.5, 1.08); l->SetLineColor(kGray + 1); l->SetLineStyle(3); l->Draw(); }
    HLine(-0.5, kNCut - 0.5, 1., kBlack, 2, 1);
    auto series = [&](int which, double off, int col, int mk, double ms) {
        std::vector<double> x; std::vector<VE> y;
        for (int k = 0; k < kNCut; ++k) {
            VE v(kNaN, kNaN);
            const EffOut& o = gEff[k];
            switch (which) {
                case 0: v = VE(o.s[eDir].eps, o.s[eDir].err); break;
                case 1: v = VE(o.s[eConv].eps, o.s[eConv].err); break;
                case 2: v = VE(o.s[eTot].eps, o.s[eTot].err); break;
                case 3: v = VE(o.s[eTotKW].eps, o.s[eTotKW].err); break;
                case 4: v = o.dPeak; break;
                case 5: v = kCut[k].bumpStart ? VE(o.sBumpStart.eps, o.sBumpStart.err) : VE(o.s[eBump].eps, o.s[eBump].err); break;
                case 6: v = kCut[k].bumpStart ? VE(o.sTagStart.eps, o.sTagStart.err) : VE(o.s[eTag].eps, o.s[eTag].err); break;
                case 7: v = o.tTagDal; break;
                case 8: v = kCut[k].bumpStart ? VE(kNaN, kNaN) : VE(o.s[eSbu].eps, o.s[eSbu].err); break;
                case 9: v = o.dTag; break;
                case 10: v = o.dEdge; break;
            }
            x.push_back(k + off); y.push_back(v);
        }
        Pts(x, y, col, mk, ms, 3)->Draw("P same");
    };
    series(0, -0.30, kColPeakDirect, 22, 2.0);
    series(1, -0.15, kColPeakConv, 21, 1.8);
    series(2, 0.0, kColPeakTotal, 20, 2.0);
    series(3, 0.13, kColPeakTotal, 24, 2.0);
    series(4, 0.30, kColData, 29, 3.4);
    f1->Draw("axis same");
    Txt(0.065, 0.93, "(a) #eta peak: simulated direct decay, conversion and total against the data peak (pass/fail fits)", 0.055);
    // (b) bump
    cv->cd(); MkPad(Nm("pe2").c_str(), 0, 0.20, 1, 0.49, false, 0.06, 0.01, 0.10, 0.12);
    TH1D* f2 = CatFrame(Nm("f").c_str(), labs, "efficiency", 0.2, 1.08, 0.075, 0.075, 0.36);
    f2->Draw("axis");
    for (int k = 1; k < kNCut; ++k) { TLine* l = new TLine(k - 0.5, 0.2, k - 0.5, 1.08); l->SetLineColor(kGray + 1); l->SetLineStyle(3); l->Draw(); }
    HLine(-0.5, kNCut - 0.5, 1., kBlack, 2, 1);
    series(5, -0.32, kColBumpDalitz, 20, 2.0);
    series(6, -0.18, kOrange + 3, 21, 1.8);
    series(7, -0.06, kOrange + 3, 25, 1.8);
    series(8, 0.07, kColBumpConv, 22, 2.0);
    series(9, 0.21, kColData, 21, 2.0);
    series(10, 0.33, kColData, 26, 2.2);
    f2->Draw("axis same");
    Txt(0.065, 0.92, "(b) #eta bump (#eta #rightarrow #pi^{+}#pi^{-}#pi^{0}): simulated classes against the data-study tagged-#eta and edge-step efficiencies", 0.068);
    // (c) pulls
    cv->cd(); MkPad(Nm("pe3").c_str(), 0, 0.0, 1, 0.20, false, 0.06, 0.01, 0.10, 0.20);
    TH1D* f3 = CatFrame(Nm("f").c_str(), labs, "(sim - data)/#sigma", -9, 9, 0.11, 0.12, 0.25);
    f3->Draw("axis");
    for (int k = 1; k < kNCut; ++k) { TLine* l = new TLine(k - 0.5, -9, k - 0.5, 9); l->SetLineColor(kGray + 1); l->SetLineStyle(3); l->Draw(); }
    HLine(-0.5, kNCut - 0.5, 0., kBlack, 1, 1);
    HLine(-0.5, kNCut - 0.5, 2., kGray + 2, 2, 1);
    HLine(-0.5, kNCut - 0.5, -2., kGray + 2, 2, 1);
    auto pull = [](VE s, VE d) {
        if (!(s.first == s.first) || !(d.first == d.first)) return VE(kNaN, kNaN);
        const double e = std::sqrt(s.second * s.second + d.second * d.second);
        return e > 0 ? VE(std::max(-8.9, std::min(8.9, (s.first - d.first) / e)), 0.) : VE(kNaN, kNaN);
    };
    for (int q = 0; q < 4; ++q) {
        if (q == 1) continue;   // no pulls for the QED illustration (review finding 3)
        std::vector<double> x; std::vector<VE> y;
        for (int k = 0; k < kNCut; ++k) {
            const EffOut& o = gEff[k];
            VE v;
            if (q == 0) v = pull(VE(o.s[eTot].eps, o.s[eTot].err), o.dPeak);
            if (q == 2) v = pull(kCut[k].bumpStart ? VE(o.sTagStart.eps, o.sTagStart.err) : VE(o.s[eTag].eps, o.s[eTag].err), o.dTag);
            if (q == 3) v = kCut[k].bumpStart ? VE(kNaN, kNaN) : pull(VE(o.s[eSbu].eps, o.s[eSbu].err), o.dEdge);
            x.push_back(k - 0.27 + 0.18 * q); y.push_back(v);
        }
        const int col[4] = {kColPeakTotal, kColPeakTotal, kOrange + 3, kColBumpConv}, mk[4] = {20, 24, 21, 22};
        Pts(x, y, col[q], mk[q], 1.8)->Draw("P same");
    }
    f3->Draw("axis same");
    Txt(0.065, 0.91, "(c) pulls: peak total vs data peak (data error stat #oplus model); tagged bump (square), bump window (triangle) vs data; |pull| < 9", 0.095);
    Save(cv, "eta_summary_e_eff");
}

// ============================================================================
// Figure f: eta peak normalisation, data / simulation by method (review finding 1)
// ============================================================================
// The peak ratio data/sim depends on how the peak content is defined. Methods (all on the same rows):
//   A    data P_main (whole Gaussian, frozen 0.5439 / 5.4 MeV) / simulated peak content over ALL M4
//   Aw   data P_main inside the window (x 0.943) / simulated peak inside the window (sim M4 x 0.9952)
//   B    data P_main / P_main applied to the simulated PEAK TEMPLATE alone (sim M4 x 0.9952, same fit)
//   C    data P_main / P_main applied to the simulated TOTAL T1 + T2 (sim M4 x 0.9952, same fit); C(a): variant (a)
//   D    template fit M1 (bump, peak, background free): alpha_peak; 0.40-0.75, no scale, T1 + T2;
//        its error scaled by sqrt(chi2/ndf); the range over 2 fit ranges x 2 mass scales x 2 ghost variants
//   E    side bands: data (N - CB) / simulated background without the peak in 0.5225-0.5336 and 0.5542-0.5650,
//        linearly interpolated to the window; implied data peak in the window / simulated peak in the window
//        (iterated for the peak tails in the bands); range: background ratio of the left or of the right band alone
// Errors: stat (data and template statistics); model = half range over the 48 fit-model variants of the data study
// (pol 1-4, 4 ranges, width -+ 0.5 MeV; methods B, C: the same variant applied to data and simulation, half range of
// the ratio); absolute scale incl. L_data: 20 % (reconcile.md: 15-20 %); ghost: A, Aw, B also against the peak without
// rows with a ghost-flagged lepton ('b'), C and E also with variant (a).
const int kNMs = 4;
const char* kMsTag[kNMs] = {"pre", "pre+OA>2", "E1", "E2"};
const double kAbsScale = 0.20;
struct MethRes {
    VE A, Aw, B, C, Ca, D, E, Ea; double mA = 0, mAw = 0, mB = 0, mC = 0, mCa = 0;   // ratio and stat error; model half range
    double Ab = 0, Awb = 0, Bb = 0;      // ghost variant (b): against the peak without ghost-lepton rows
    double Dlo = 0, Dhi = 0, Dchi = 0; int Dndf = 0; double Dscaled = 0;
    double Elo = 0, Ehi = 0, Eunsc = 0, rL = 0, rR = 0, rW = 0;
    double dS = 0, dSe = 0, dSm = 0, simAll = 0, simAllE = 0, simWin = 0, pkOnly = 0, simTot = 0, simTotA = 0;
};
MethRes gMeth[kNMs];
bool gHaveF = false;
// P_main with a general configuration on a simulated spectrum (weighted linear least squares, template variances)
double PmainCfg(const std::vector<double>& y, const std::vector<double>& var, double lo, double bw, const scan::PeakCfg& c, double* err = nullptr) {
    const double x0 = 0.5 * (c.lo + c.hi), dx = 0.5 * (c.hi - c.lo);
    double vmax = 0;
    for (double v : var) vmax = std::max(vmax, v);
    std::vector<double> yy, vv, F;
    const int np = 2 + c.order;
    for (size_t i = 0; i < y.size(); ++i) {
        const double xc = lo + (i + 0.5) * bw;
        if (!(xc >= c.lo && xc < c.hi)) continue;
        yy.push_back(y[i]); vv.push_back(std::max(var[i], 1e-6 * vmax + 1e-12));
        const double t = (xc - x0) / dx;
        F.push_back(scan::GaussFrac(xc - 0.5 * bw, xc + 0.5 * bw, c.mean, c.sigma));
        double tp = 1; for (int j = 0; j <= c.order; ++j) { F.push_back(tp); tp *= t; }
    }
    const scan::LinRes r = scan::LinFit((int)yy.size(), np, yy.data(), vv.data(), F.data());
    if (err) *err = std::sqrt(std::max(0., r.cov[0][0]));
    return r.ok ? r.par[0] : kNaN;
}
// the 48 configurations of scan::PeakModelSpread, in its order
std::vector<scan::PeakCfg> VarCfgs() {
    std::vector<scan::PeakCfg> v;
    for (int ir = 0; ir < scan::kNVarRange; ++ir)
        for (int ord = scan::kVarOrdMin; ord <= scan::kVarOrdMax; ++ord)
            for (int iw = -1; iw <= 1; ++iw) {
                scan::PeakCfg c; c.lo = scan::kVarLo[ir]; c.hi = scan::kVarHi[ir]; c.order = ord; c.sigma = scan::kPeakSigma + iw * scan::kVarDSigma;
                v.push_back(c);
            }
    return v;
}
void PartF() {
    TStopwatch sw;
    const SimData& d = *gD;
    const double sc = seen::kScaleK0S;
    const std::vector<scan::PeakCfg> vc = VarCfgs();
    for (int is = 0; is < kNMs; ++is) {
        MethRes& m = gMeth[is];
        Sel ssel; scan::Sel dsel;
        if (is == 0) { ssel = SelPreset(kPre); dsel = DataSelPreset(kPre); }
        if (is == 1) { ssel = SelPreset(kPreOA2); dsel = DataSelPreset(kPreOA2); }
        if (is == 2) { ssel = seen::SimSel(seen::sE1); dsel = seen::DatSel(seen::sE1); }
        if (is == 3) { ssel = seen::SimSel(seen::sE2); dsel = seen::DatSel(seen::sE2); }
        // ---- data: P_main and its model variants; N - CB in the side bands
        scan::FillOpt fo;
        fo.windows = {{0.5225, scan::kPeakLo}, {scan::kPeakLo, scan::kPeakHi}, {scan::kPeakHi, 0.5650}};
        const scan::Spectra sp = scan::FillMany(*gDD, {dsel}, fo)[0];
        const sep::DFit df = sep::DoFit(sp);
        m.dS = df.S; m.dSe = df.eS; m.dSm = df.eM;
        double D[3], DV[3];
        for (int b = 0; b < 3; ++b) { D[b] = sp.win[2 + b].sig(); DV[b] = std::pow(sp.win[2 + b].sigErr(), 2); }
        // ---- simulation
        const Comp cPk = CompPeakTotal();
        const Comp cPkNG{"peak, no ghost lepton", [](const Row& r) { return r.tmpl == 1 && (r.cls == 4 || r.cls == 5) && !GhostLepton(r); }};
        const Comp cT = CompSim(), cTa = CompSimVar(1);
        const Comp cB{"sim without peak", [](const Row& r) { return ((r.tmpl == 1 && !X17Row(r)) || r.tmpl == 2) && !(r.tmpl == 1 && (r.cls == 4 || r.cls == 5)); }};
        const Comp cBa{"(a) without peak", [cTa](const Row& r) { return cTa.f(r) && !(r.tmpl == 1 && (r.cls == 4 || r.cls == 5)); }};
        const std::vector<Comp> cc = {cPk, cPkNG, cT, cTa, cB, cBa};
        enum { kPk = 0, kPkNG, kT, kTa, kB, kBa };
        const std::vector<WSum> all = SumMany(d, ssel, cc, gOpt);
        m.simAll = all[kPk].S * kToRows; m.simAllE = std::sqrt(all[kPk].V) * kToRows;
        const double simAllNG = all[kPkNG].S * kToRows, simAllNGE = std::sqrt(all[kPkNG].V) * kToRows;
        // method A
        auto ratio = [](double a, double ea, double b, double eb) { const double r = a / b; return VE(r, std::fabs(r) * std::hypot(ea / a, eb / b)); };
        m.A = ratio(m.dS, m.dSe, m.simAll, m.simAllE); m.mA = m.A.first * m.dSm / m.dS;
        m.Ab = m.dS / simAllNG;
        // binned simulation, M4 x 0.9952, the data binning
        Binning bP; bP.nbin = sp.nbin; bP.lo = sp.lo; bP.hi = sp.hi;
        bP.x = [sc](const Row& r) { return r.m4 * sc; };
        const std::vector<WSpec> T = FillMany(d, ssel, cc, gOpt, bP);
        std::vector<std::vector<double>> y(cc.size(), std::vector<double>(bP.nbin)), v = y;
        for (size_t c = 0; c < cc.size(); ++c)
            for (int i = 0; i < bP.nbin; ++i) { y[c][i] = T[c].s[i] * kToRows; v[c][i] = T[c].v[i] * kToRows * kToRows; }
        // method Aw: window contents (data: whole Gaussian x window fraction; sim: K0S-scaled M4 inside the window)
        {
            const double fw = scan::GaussFrac(scan::kPeakLo, scan::kPeakHi, scan::kPeakMean, scan::kPeakSigma);
            const Sel inW = SelAnd(ssel, [sc](const Row& r) { const double x = r.m4 * sc; return x >= scan::kPeakLo && x < scan::kPeakHi; });
            const std::vector<WSum> w = SumMany(d, inW, {cPk, cPkNG}, gOpt);
            m.simWin = w[0].S * kToRows;
            m.Aw = ratio(m.dS * fw, m.dSe * fw, w[0].S * kToRows, std::sqrt(w[0].V) * kToRows); m.mAw = m.Aw.first * m.dSm / m.dS;
            m.Awb = m.dS * fw / (w[1].S * kToRows);
        }
        // methods B and C: nominal and per variant (the same configuration on data and simulation)
        auto pmRatio = [&](int comp, VE& out, double& model, double* simS) {
            double e = 0;
            const double s0 = PmainCfg(y[comp], v[comp], bP.lo, bP.bw(), scan::PeakCfg(), &e);
            if (simS) *simS = s0;
            out = ratio(m.dS, m.dSe, s0, e);
            double lo = 1e30, hi = -1e30; int n = 0;
            for (size_t k = 0; k < vc.size() && k < df.vS.size(); ++k) {
                if (!df.vOk[k]) continue;
                const double sv = PmainCfg(y[comp], v[comp], bP.lo, bP.bw(), vc[k]);
                if (!(sv > 0)) continue;
                const double r = df.vS[k] / sv; lo = std::min(lo, r); hi = std::max(hi, r); ++n;
            }
            model = n >= 2 ? 0.5 * (hi - lo) : out.first * m.dSm / m.dS;
        };
        double tmp = 0;
        pmRatio(kPk, m.B, m.mB, &m.pkOnly);
        { const double sNG = PmainCfg(y[kPkNG], v[kPkNG], bP.lo, bP.bw(), scan::PeakCfg()); m.Bb = m.dS / sNG; }
        pmRatio(kT, m.C, m.mC, &m.simTot);
        pmRatio(kTa, m.Ca, m.mCa, &m.simTotA);
        (void)tmp;
        // method E: side bands (scaled and unscaled simulation)
        for (int us = 0; us < 2; ++us) {
            const double s = us ? 1. : sc;
            const double ed[4] = {0.5225, scan::kPeakLo, scan::kPeakHi, 0.5650};
            double Bk[2][3], BkV[2][3], Pk[3], PkV[3];
            for (int b = 0; b < 3; ++b) {
                const double lo = ed[b], hi = ed[b + 1];
                const Sel sb = SelAnd(ssel, [s, lo, hi](const Row& r) { const double x = r.m4 * s; return x >= lo && x < hi; });
                const std::vector<WSum> w = SumMany(d, sb, {cB, cBa, cPk}, gOpt);
                for (int g = 0; g < 2; ++g) { Bk[g][b] = w[g].S * kToRows; BkV[g][b] = w[g].V * kToRows * kToRows; }
                Pk[b] = w[2].S * kToRows; PkV[b] = w[2].V * kToRows * kToRows;
            }
            const double cL = 0.5 * (ed[0] + ed[1]), cW = 0.5 * (ed[1] + ed[2]), cR = 0.5 * (ed[2] + ed[3]);
            const double a = (cR - cW) / (cR - cL);
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
            if (us == 0) {
                m.E = solve(0, 0, &m.rL, &m.rR, &m.rW);
                const double k1 = solve(0, 1, nullptr, nullptr, nullptr).first, k2 = solve(0, 2, nullptr, nullptr, nullptr).first;
                m.Elo = std::min(k1, k2); m.Ehi = std::max(k1, k2);
                m.Ea = solve(1, 0, nullptr, nullptr, nullptr);
            } else
                m.Eunsc = solve(0, 0, nullptr, nullptr, nullptr).first;
        }
        // method D: template fit M1 (2 ranges x 2 scales x 2 ghost variants)
        {
            const std::vector<Comp> fc = seen::FitComps();
            struct Rng { double lo, hi; int nb; };
            const Rng rng[2] = {{0.40, 0.75, 70}, {0.43, 0.70, 54}};
            m.Dlo = 1e9; m.Dhi = -1e9;
            for (int ir = 0; ir < 2; ++ir) {
                const DSpec DS = DataSpectrum(*gDD, dsel, rng[ir].nb, rng[ir].lo, rng[ir].hi);
                std::vector<int> bins; for (int i = 0; i < DS.nbin; ++i) bins.push_back(i);
                for (int scl = 0; scl < 2; ++scl) {
                    Binning b; b.nbin = rng[ir].nb; b.lo = rng[ir].lo; b.hi = rng[ir].hi;
                    const double f = scl ? sc : 1.;
                    b.x = [f](const Row& r) { return (double)r.m4 * f; };
                    const std::vector<WSpec> TS = FillMany(d, ssel, fc, gOpt, b);
                    for (int gv = 0; gv < 2; ++gv) {
                        const int grp[3] = {seen::fBump, seen::fPeak, gv ? seen::fKA : seen::fK};
                        std::vector<std::vector<double>> TT, TV;
                        for (int k : grp) {
                            std::vector<double> t(b.nbin), tv(b.nbin);
                            for (int i = 0; i < b.nbin; ++i) { t[i] = TS[k].s[i] * kToRows; tv[i] = TS[k].v[i] * kToRows * kToRows; }
                            TT.push_back(t); TV.push_back(tv);
                        }
                        const seen::FitRes r = seen::TemplateFit(DS.sp, bins, TT, TV);
                        if (!r.ok) continue;
                        m.Dlo = std::min(m.Dlo, r.a[1]); m.Dhi = std::max(m.Dhi, r.a[1]);
                        if (ir == 0 && scl == 0 && gv == 0) {
                            m.D = VE(r.a[1], r.ea[1]); m.Dchi = r.chi2; m.Dndf = r.ndf;
                            m.Dscaled = r.ea[1] * std::sqrt(std::max(1., r.chi2 / r.ndf));
                        }
                    }
                }
            }
        }
        printf("PartF %s: A %.3f  Aw %.3f  B %.3f  C %.3f  C(a) %.3f  D %.3f  E %.3f [%.3f, %.3f]\n", kMsTag[is], m.A.first, m.Aw.first, m.B.first,
               m.C.first, m.Ca.first, m.D.first, m.E.first, m.Elo, m.Ehi);
    }
    // ---- checks against the study tables
    C("\n## figure f: peak methods against tables/eta_peak_composition.txt (method A) and tables/eta_enriched_fit.txt (methods B, C, D)\n");
    Check("pre: method A data/sim (eta_peak_composition)", gMeth[0].A.first, GetVE("eta_peak_composition.txt", "", "pre", -1, 8).first, 6e-4);
    for (int j = 0; j < 2; ++j) {
        const int is = kJ[j];
        if (!gPM[is].ok) continue;
        Check(Form("%s: data P_main [rows]", kJTag[j]), gMeth[2 + j].dS, gPM[is].data.first, 0.6);
        Check(Form("%s: P_main on the peak template [rows]", kJTag[j]), gMeth[2 + j].pkOnly, gPM[is].pkOnly.first, 0.6);
        Check(Form("%s: P_main on T1+T2 [rows]", kJTag[j]), gMeth[2 + j].simTot, gPM[is].sim54.first, 0.6);
        Check(Form("%s: P_main on (a) [rows]", kJTag[j]), gMeth[2 + j].simTotA, gPM[is].simA.first, 0.6);
        const FitRow* f = FindFit(gFits, kJTag[j], 0, false, false, "M1");
        Check(Form("%s: M1 alpha_peak", kJTag[j]), gMeth[2 + j].D.first, f ? f->Par("a_peak").first : kNaN, 6e-4);
    }
    // ---- table
    FILE* f = fopen(Form("%s/eta_summary_peak_methods.txt", kTabDir), "w");
    auto W = [&](const char* fmt, ...) {
        char b[8192]; va_list ap; va_start(ap, fmt); vsnprintf(b, sizeof b, fmt, ap); va_end(ap); fputs(b, stdout); if (f) fputs(b, f);
    };
    W("# eta_summary_peak_methods: the eta peak, data / simulation, by method (review finding 1: the ratio depends on the method)\n");
    W("# simulation: std + lep + eta owned rows, mode A, w (PDG); peak = direct eta -> pi+pi-e+e- (generator M(ee) model not QED) + conversion\n");
    W("# data: rows of scan_compact (vz > -160 mm, isBest >= 0, PT3, cor_ok); L_data = 6080.5 nb^-1 UNCONFIRMED: every ratio scales with 6080.5 / L_data\n");
    W("# A  : data P_main (whole Gaussian 0.5439 / 5.4 MeV, pol2, 0.5225-0.61) / simulated peak content over ALL M4\n");
    W("# Aw : data P_main x 0.943 (window share) / simulated peak inside the window 0.5336-0.5542 (sim M4 x 0.9952)\n");
    W("# B  : data P_main / P_main on the simulated PEAK TEMPLATE alone (sim M4 x 0.9952)\n");
    W("# C  : data P_main / P_main on the simulated TOTAL T1 + T2 (sim M4 x 0.9952); C(a): on variant (a)\n");
    W("# D  : template fit M1 alpha_peak (0.40-0.75, no scale, T1 + T2), error x sqrt(chi2/ndf); range over 2 ranges x 2 scales x 2 ghost variants\n");
    W("# E  : side bands: data / simulated background without the peak in 0.5225-0.5336 (left) and 0.5542-0.5650 (right), linear interpolation to\n");
    W("#      the window; (data window - r x sim background) / sim peak in the window; range = left or right band ratio alone; E(a): variant (a)\n");
    W("# errors: ratio +- stat +- model (fit-model spread: 48 variants, the same variant on data and simulation for B and C); 'tot' adds the\n");
    W("#   absolute scale incl. L_data (%.0f %%); ghost: (b) = against the simulated peak WITHOUT rows with a ghost-flagged lepton\n\n", kAbsScale * 100);
    W("%-9s data P_main [rows] (stat, model)    sim peak all M4   sim peak in window   P_main(peak tmpl)   P_main(T1+T2)   P_main((a))\n", "sel");
    for (int is = 0; is < kNMs; ++is) {
        const MethRes& m = gMeth[is];
        W("%-9s %6.0f +- %-4.0f +- %-5.0f (%.0f %%)    %6.0f +- %-4.0f    %6.0f              %6.0f              %6.0f          %6.0f\n", kMsTag[is], m.dS, m.dSe, m.dSm,
          100 * m.dSm / m.dS, m.simAll, m.simAllE, m.simWin, m.pkOnly, m.simTot, m.simTotA);
    }
    W("\n%-9s %-6s %-36s %-8s %-10s\n", "sel", "method", "data/sim +- stat +- model", "tot", "ghost var.");
    for (int is = 0; is < kNMs; ++is) {
        const MethRes& m = gMeth[is];
        auto line = [&](const char* tag, VE r, double mod, double gvar, const char* gname) {
            const double tot = std::sqrt(r.second * r.second + mod * mod + std::pow(kAbsScale * r.first, 2));
            W("%-9s %-6s %6.3f +- %-6.3f +- %-6.3f          %6.3f   %s %.3f\n", kMsTag[is], tag, r.first, r.second, mod, tot, gname, gvar);
        };
        line("A", m.A, m.mA, m.Ab, "(b)");
        line("Aw", m.Aw, m.mAw, m.Awb, "(b)");
        line("B", m.B, m.mB, m.Bb, "(b)");
        line("C", m.C, m.mC, m.Ca.first, "(a)");
        {
            const VE r(m.D.first, m.Dscaled);
            const double tot = std::sqrt(r.second * r.second + std::pow(kAbsScale * r.first, 2));
            W("%-9s %-6s %6.3f +- %-6.3f (x sqrt(chi2/ndf), chi2/ndf %.0f/%d)  %6.3f   range %.3f-%.3f\n", kMsTag[is], "D", r.first, r.second, m.Dchi, m.Dndf, tot, m.Dlo, m.Dhi);
        }
        {
            const double tot = std::sqrt(m.E.second * m.E.second + std::pow(kAbsScale * m.E.first, 2));
            W("%-9s %-6s %6.3f +- %-6.3f   range %.3f-%.3f (bkg ratio left %.3f, right %.3f, window %.3f)  %6.3f   (a) %.3f; unscaled sim M4 %.3f\n", kMsTag[is], "E",
              m.E.first, m.E.second, m.Elo, m.Ehi, m.rL, m.rR, m.rW, tot, m.Ea.first, m.Eunsc);
        }
    }
    // overall ranges
    double pl = 1e9, ph = -1e9, ol = 1e9, oh = -1e9;
    for (int is = 0; is < kNMs; ++is) {
        const MethRes& m = gMeth[is];
        for (double x : {m.B.first, m.C.first, m.Ca.first}) { pl = std::min(pl, x); ph = std::max(ph, x); }
        for (double x : {m.D.first, m.E.first, m.Ea.first}) { ol = std::min(ol, x); oh = std::max(oh, x); }
    }
    W("\n# P_main on both sides (B, C, C(a)): %.2f-%.2f; methods using the simulated background shape (D, E, E(a)): %.2f-%.2f; all-M4 content (A): %.2f-%.2f\n", pl, ph, ol, oh,
      std::min({gMeth[0].A.first, gMeth[1].A.first, gMeth[2].A.first, gMeth[3].A.first}), std::max({gMeth[0].A.first, gMeth[1].A.first, gMeth[2].A.first, gMeth[3].A.first}));
    W("# reading: the data peak is not larger than the simulated one; by the data-study method (P_main) it is %.2f-%.2f of the simulation,\n", pl, ph);
    W("#   methods that use the simulated background shape give %.2f-%.2f; the peak normalisation is known only to a factor of about 2\n", ol, oh);
    W("#   (a 10-15 %% change of the background level under the peak moves data/sim of the peak by a factor 2-3: method E range)\n");
    if (f) fclose(f);
    gHaveF = true;
    printf("PartF: %.1f s\n", sw.RealTime());
}
void FigF() {
    TCanvas* cv = new TCanvas(Nm("cF").c_str(), "", 2600, 1500);
    cv->cd();
    MkPad(Nm("hF").c_str(), 0, 0.84, 1, 1, false, 0, 0, 0, 0);
    Txt(0.01, 0.80, "#eta peak: data / simulation by method - the ratio depends on how the peak content is defined (not a measurement of the #eta yield)", 0.17);
    Txt(0.01, 0.60, "simulation: SMASH std + lep + eta owned rows, w (PDG), peak = direct #eta#rightarrow#pi^{+}#pi^{-}e^{+}e^{-} (generator M_{ee} model not QED) + conversion #eta#rightarrow#pi^{+}#pi^{-}#gamma; "
                    "data: N - CB, v_{z} > -160 mm, isBest #geq 0, PT3, cor_ok; L_{data} = 6080.5 nb^{-1} UNCONFIRMED", 0.10);
    Txt(0.01, 0.42, "A: data P_{main} / sim content, all M_{4};  Aw: inside the window;  B: P_{main} on data and on the sim peak template;  C: P_{main} on data and on the sim total "
                    "(open: variant (a));  D: template fit M1, #alpha_{peak} (dashed bar: range over fit ranges,", 0.095);
    Txt(0.01, 0.24, "mass scale, ghost variant);  E: side bands, background ratio interpolated (dashed bar: left or right band alone; open: variant (a));  error bars: thick = stat #oplus fit-model spread "
                    "(48 variants), thin = also the absolute scale incl. L_{data} (#pm20%)", 0.095);
    Txt(0.01, 0.06, "pre = isBest #geq 0, PT3, v_{z} > -160 mm, corrected momenta;  E1 = pre + #alpha_{lab} < 50#circ, #alpha_{#eta} > 140#circ, M_{#pi#pi} < 0.46, OA > 2#circ;  E2 = E1 without OA > 2#circ;  "
                    "open circles on A/Aw/B: sim peak without ghost-lepton rows", 0.095);
    cv->cd(); MkPad(Nm("pF").c_str(), 0, 0, 0.74, 0.84, false, 0.07, 0.01, 0.03, 0.13);
    const std::vector<std::string> lab = {"A", "Aw", "B", "C", "D", "E"};
    TH1D* fr = CatFrame(Nm("f").c_str(), lab, "#eta peak: data / simulation", 0.0, 2.2, 0.05, 0.065, 0.65);
    fr->GetXaxis()->SetTitle("method"); fr->GetXaxis()->SetTitleSize(0.05);
    fr->Draw("axis");
    for (int k = 1; k < 6; ++k) { TLine* l = new TLine(k - 0.5, 0.0, k - 0.5, 2.2); l->SetLineColor(kGray + 1); l->SetLineStyle(3); l->Draw(); }
    HLine(-0.5, 5.5, 1., kBlack, 2, 2);
    const int col[kNMs] = {kBlack, kGray + 2, kRed + 1, kBlue + 1}, mk[kNMs] = {20, 21, 33, 22};
    TLegend* lg = new TLegend(0.08, 0.80, 0.60, 0.96); lg->SetNColumns(4); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.040);
    for (int is = 0; is < kNMs; ++is) {
        const MethRes& m = gMeth[is];
        const double off = -0.27 + 0.18 * is;
        auto pt = [&](double x, VE r, double mod, int mst, bool scaleOnly) {
            const double inner = std::hypot(r.second, mod), outer = std::sqrt(inner * inner + std::pow(kAbsScale * r.first, 2));
            TGraphErrors* go = new TGraphErrors(1); go->SetPoint(0, x, r.first); go->SetPointError(0, 0, outer);
            go->SetLineColor(col[is]); go->SetLineWidth(1); go->SetMarkerSize(0); go->Draw("P same");
            TGraphErrors* gi = new TGraphErrors(1); gi->SetPoint(0, x, r.first); gi->SetPointError(0, 0, inner);
            gi->SetLineColor(col[is]); gi->SetLineWidth(4); gi->SetMarkerColor(col[is]); gi->SetMarkerStyle(mst); gi->SetMarkerSize(2.2); gi->Draw("P same");
            (void)scaleOnly;
            return gi;
        };
        auto open = [&](double x, double y) {
            TGraph* g = new TGraph(1); g->SetPoint(0, x, y); g->SetMarkerStyle(24); g->SetMarkerSize(1.8); g->SetMarkerColor(col[is]); g->Draw("P same");
        };
        auto bar = [&](double x, double lo, double hi) {
            TLine* l = new TLine(x + 0.035, lo, x + 0.035, hi); l->SetLineColor(col[is]); l->SetLineStyle(2); l->SetLineWidth(3); l->Draw();
        };
        TGraphErrors* g0 = pt(0 + off, m.A, m.mA, mk[is], false); open(0 + off, m.Ab);
        pt(1 + off, m.Aw, m.mAw, mk[is], false); open(1 + off, m.Awb);
        pt(2 + off, m.B, m.mB, mk[is], false); open(2 + off, m.Bb);
        pt(3 + off, m.C, m.mC, mk[is], false); open(3 + off, m.Ca.first);
        pt(4 + off, VE(m.D.first, m.Dscaled), 0., mk[is], false); bar(4 + off, m.Dlo, m.Dhi);
        pt(5 + off, m.E, 0., mk[is], false); bar(5 + off, m.Elo, m.Ehi); open(5 + off, m.Ea.first);
        lg->AddEntry(g0, kMsTag[is], "p");
    }
    lg->Draw();
    fr->Draw("axis same");
    // right: reading
    cv->cd(); MkPad(Nm("pT").c_str(), 0.74, 0, 1, 0.84, false, 0.03, 0.02, 0.03, 0.10);
    double pl = 1e9, ph = -1e9, ol = 1e9, oh = -1e9;
    for (int is = 0; is < kNMs; ++is) {
        const MethRes& m = gMeth[is];
        for (double x : {m.B.first, m.C.first, m.Ca.first}) { pl = std::min(pl, x); ph = std::max(ph, x); }
        for (double x : {m.D.first, m.E.first, m.Ea.first}) { ol = std::min(ol, x); oh = std::max(oh, x); }
    }
    Txt(0.04, 0.94, "#bf{reading}", 0.075);
    Txt(0.04, 0.86, "The data peak is not larger", 0.060);
    Txt(0.04, 0.80, "than the simulated one.", 0.060);
    Txt(0.04, 0.72, Form("P_{main} on both sides (B, C):"), 0.060);
    Txt(0.04, 0.66, Form("   %.2f - %.2f", pl, ph), 0.060);
    Txt(0.04, 0.58, "with the simulated background", 0.060);
    Txt(0.04, 0.52, Form("shape (D, E): %.2f - %.2f", ol, oh), 0.060);
    Txt(0.04, 0.44, "The peak normalisation is known", 0.060);
    Txt(0.04, 0.38, "only to a factor of about 2.", 0.060);
    Txt(0.04, 0.28, "A (all-M_{4} content) adds the", 0.055, kGray + 2);
    Txt(0.04, 0.22, "23% of the simulated peak outside", 0.055, kGray + 2);
    Txt(0.04, 0.16, "the window (sim peak 6.7 MeV wide).", 0.055, kGray + 2);
    Save(cv, "eta_summary_f_peak_methods");
}

// ============================================================================
// Table eta_summary.txt
// ============================================================================
void WriteSummary() {
    gT = fopen(Form("%s/eta_summary.txt", kTabDir), "w");
    T("# eta_summary: headline numbers of the eta-set step (sim_eta_summary.C): SMASH std + lep + eta against data in eta-enriched selections, M4 < 1 GeV/c2\n");
    T("# simulation: std (owned, mode A) + lep (owned, mode A) + eta (owned), w (PDG); sigma = sum w / L_sim, L_sim = 221.05 nb^-1 per set; T1 + T2; no CB\n");
    T("#   variant (a): T1 + T2 + ghost-flagged true pairs (sim_ghostvar.h); X17 for BR = 1e-5 is information only, never in a sum\n");
    T("# data: (N - CB) / L_data, CB = 2 sqrt(N++ N--) in 10 MeV bins; L_data = 6080.5 nb^-1 UNCONFIRMED (every absolute ratio and alpha scales with 6080.5 / L_data)\n");
    T("# selections: pre = isBest >= 0, PT3, vz > -160 mm, corrected momenta; E1 (slides) = pre + oa_lab < 50, eta_rest > 140, M(pipi) < 0.46, OA > 2; E2 = E1 without OA > 2\n");
    T("# errors: statistical (event groups for the simulation, Poisson + CB propagation for data) unless stated; data PEAK contents and peak\n");
    T("#   efficiencies also carry the fit-model spread of the data study (48 variants: pol 1-4, 4 ranges, width -+ 0.5 MeV); std-lep covariance neglected\n");
    T("# the simulated direct decay eta -> pi+pi-e+e- has a NON-QED generator M(ee) (accepted density largest at threshold, close to 1/M^2,\n");
    T("#   tables/eta_peak_mee_true.txt): every simulated direct-decay number is labelled 'generator M(ee) model not QED' and is not for physics\n");
    T("# [table] = read from a study table (source given), everything else recomputed by this macro; cross-checks in tables/eta_summary_checks.txt\n");
    T("# revised after the adversarial review of 2026-11-12 (six findings): see section 0\n");
    if (gHaveF && gHaveAB && gHaveC) {
        double pl = 1e9, ph = -1e9, ol = 1e9, oh = -1e9, al = 1e9, ah = -1e9;
        for (int is = 0; is < kNMs; ++is) {
            const MethRes& m = gMeth[is];
            for (double x : {m.B.first, m.C.first, m.Ca.first}) { pl = std::min(pl, x); ph = std::max(ph, x); }
            for (double x : {m.D.first, m.E.first, m.Ea.first}) { ol = std::min(ol, x); oh = std::max(oh, x); }
            al = std::min(al, m.A.first); ah = std::max(ah, m.A.first);
        }
        double bl = 1e9, bh = -1e9;
        for (const FitRow& f : gFits) if (f.model == "M1" || f.model == "M3") { const VE v = f.Par("a_bump"); if (v.first == v.first) { bl = std::min(bl, v.first); bh = std::max(bh, v.first); } }
        double rl = 1e9, rh = -1e9;
        for (int j = 0; j < 2; ++j) for (int k : {5, 6}) for (double x : {rat[j][k], ratA[j][k], ratC[j][k], ratCA[j][k]}) { rl = std::min(rl, x); rh = std::max(rh, x); }
        T("\n## 0. Headline statements (corrected after the review)\n");
        T("1. eta PEAK: the data peak is not larger than the simulated one. By the data-study method (P_main on data and simulation) it is %.2f-%.2f\n", pl, ph);
        T("   of the simulation; methods that use the simulated background shape (template fit, side bands) give %.2f-%.2f. The peak normalisation is\n", ol, oh);
        T("   known only to a factor of about 2 (section 6). WITHDRAWN: 'about 3.3 times larger' (method A, all-M4 content: %.2f-%.2f) and 'about 2 times'.\n", al, ah);
        T("2. The simulated direct decay eta -> pi+pi-e+e- has a non-QED generator M(ee) (largest at threshold, close to 1/M^2): its yield, the\n");
        T("   direct/conversion ratio, the class shapes and the OA efficiencies are not for physics. To report to the producers (PLUTO database_decays_eta.dat).\n");
        T("3. The QED-reshaped direct decay is an ILLUSTRATION of shapes only (flat acceptance and fixed content assumed); no fit result of it is quoted.\n");
        T("4. Qualitatively the data peak has far more content at OA > 10 deg and M(ee) > 20 MeV than either simulated component (section 4; data errors\n");
        T("   now include the fit-model spread, so the chi2 values are smaller than before).\n");
        T("5. eta BUMP: the simulated bump is too large relative to the simulated background at the bump edge; the fit factor lies between about %.2f\n", bl);
        double el = 1e9, eh = -1e9;
        for (const FitRow& f : gFits) if (f.model == "M3") { const VE v = f.Par("a_bkgEta"); if (v.first == v.first) { el = std::min(el, v.first); eh = std::max(eh, v.first); } }
        T("   and %.2f depending on model, range, mass scale and ghost variant (poor fits, not a yield). The eta-related background gets %.2f-%.2f (mostly > 1)\n", bh, el, eh);
        T("   in the same fits, which contradicts a bump made too large by the same etas: never combine bump and peak factors into an eta production factor.\n");
        T("6. Above the eta structures the spectrum is dominated by multi-pion background (sim share E1/E2: %.0f/%.0f %% in 0.56-0.70, %.0f/%.0f %% in 0.70-1.00), which\n",
          100 * sReg[0][seen::cKM].s[5] / sReg[0][seen::cSum].s[5], 100 * sReg[1][seen::cKM].s[5] / sReg[1][seen::cSum].s[5],
          100 * sReg[0][seen::cKM].s[6] / sReg[0][seen::cSum].s[6], 100 * sReg[1][seen::cKM].s[6] / sReg[1][seen::cSum].s[6]);
        T("   SMASH does not describe. Its level there is %.2f-%.2f times SMASH (E1, E2; T1+T2 or variant (a); vz > -160 or the central slice\n", rl, rh);
        T("   -140 < vz < -100 mm), with L_data unconfirmed. No selection reaches eta dominance outside the bump window (sim purity E1/E2: bump window\n");
        T("   %.2f/%.2f, peak window %.2f/%.2f, 0.56-0.70 %.2f/%.2f). data/sim close to 1 below the bump is a coincidence of components, not a validation.\n",
          purReg[0][1].eps, purReg[1][1].eps, purReg[0][3].eps, purReg[1][3].eps, purReg[0][5].eps, purReg[1][5].eps);
    }
    if (gHaveAB) {
        T("\n## 1. Simulated cross sections in E1 and E2 [nb] (M4 < 1.0 / bump window 0.455-0.5225 / peak window 0.5336-0.5542)\n");
        T("%-34s | %-50s | %-50s\n", "component", "E1: M4 < 1.0      bump window      peak window", "E2: M4 < 1.0      bump window      peak window");
        const int cs[9] = {seen::cBD, seen::cBC, seen::cPC, seen::cPD, -1, seen::cKE, seen::cKM, seen::cSum, seen::cSumA};
        const char* cn[9] = {"bump: pi0 Dalitz", "bump: conversion", "peak: conversion (pi+pi-gamma)", "peak: direct (gen. M(ee) not QED)",
                             "peak: total", "phys. bkg, eta-related", "phys. bkg, multi-pion", "simulation T1 + T2", "variant (a)"};
        for (int q = 0; q < 9; ++q) {
            T("%-34s |", cn[q]);
            for (int j = 0; j < 2; ++j) {
                double v[3], e[3];
                if (cs[q] >= 0) {
                    v[0] = sTot[j][cs[q]].Xsec(); e[0] = sTot[j][cs[q]].XsecErr();
                    for (int w = 0; w < 2; ++w) { const int k = w == 0 ? 1 : 3; v[w + 1] = sReg[j][cs[q]].s[k] / kLsim; e[w + 1] = std::sqrt(sReg[j][cs[q]].v[k]) / kLsim; }
                } else {
                    v[0] = sTot[j][seen::cPC].Xsec() + sTot[j][seen::cPD].Xsec();
                    e[0] = std::hypot(sTot[j][seen::cPC].XsecErr(), sTot[j][seen::cPD].XsecErr());
                    for (int w = 0; w < 2; ++w) {
                        const int k = w == 0 ? 1 : 3;
                        v[w + 1] = (sReg[j][seen::cPC].s[k] + sReg[j][seen::cPD].s[k]) / kLsim;
                        e[w + 1] = std::sqrt(sReg[j][seen::cPC].v[k] + sReg[j][seen::cPD].v[k]) / kLsim;
                    }
                }
                T(" %7.3f +- %-6.3f %7.3f +- %-6.3f %7.3f +- %-6.3f |", v[0], e[0], v[1], e[1], v[2], e[2]);
            }
            T("\n");
        }
        T("%-34s |", "X17, BR 1e-5 (information)");
        for (int j = 0; j < 2; ++j) T(" %7.4f           %7.4f           %7.4f          |", sTot[j][seen::cX17].Xsec(), sReg[j][seen::cX17].s[1] / kLsim, sReg[j][seen::cX17].s[3] / kLsim);
        T("\n%-34s |", "purity p p eta (X), sim");
        for (int j = 0; j < 2; ++j) T(" %7.3f +- %-6.3f %7.3f +- %-6.3f %7.3f +- %-6.3f |", pur[j], purE[j], purReg[j][1].eps, purReg[j][1].err, purReg[j][3].eps, purReg[j][3].err);
        T("\n# purity = (signal classes 2-5 + physical-background rows with a track from an eta) / (T1 + T2)\n");
        T("# expected data rows of the simulated peak total in the peak window (x L_data / L_sim): E1 %.0f +- %.0f, E2 %.0f +- %.0f\n",
          (sReg[0][seen::cPC].s[3] + sReg[0][seen::cPD].s[3]) * kToRows, std::sqrt(sReg[0][seen::cPC].v[3] + sReg[0][seen::cPD].v[3]) * kToRows,
          (sReg[1][seen::cPC].s[3] + sReg[1][seen::cPD].s[3]) * kToRows, std::sqrt(sReg[1][seen::cPC].v[3] + sReg[1][seen::cPD].v[3]) * kToRows);

        T("\n## 2. Data (N - CB) against the simulation per M4 region (absolute; nb)\n");
        T("%-4s %-28s %20s %20s %20s %18s %18s\n", "sel", "region", "data N-CB [nb]", "sim T1+T2 [nb]", "sim (a) [nb]", "data/sim", "data/sim (a)");
        for (int j = 0; j < 2; ++j)
            for (int k = 0; k <= seen::kNReg; ++k)
                T("%-4s %-28s %9.3f +- %-7.3f %9.3f +- %-7.3f %9.3f +- %-7.3f %7.3f +- %-7.3f %7.3f +- %-7.3f\n", kJTag[j], seen::kReg[k].name, dNb[j][k], dNbE[j][k],
                  sNb[j][k], sNbE[j][k], aNb[j][k], aNbE[j][k], rat[j][k], ratE[j][k], ratA[j][k], ratAE[j][k]);
        T("\n# the same ratios in the central target slice -140 < vz < -100 mm (data and simulation), T1+T2 | variant (a); and the simulated shares\n");
        T("%-4s %-28s %20s %12s | %10s %10s\n", "sel", "region", "data/sim central", "(a) central", "purity", "multi-pion");
        for (int j = 0; j < 2; ++j)
            for (int k = 0; k <= seen::kNReg; ++k) {
                const double mp = k < seen::kNReg ? sReg[j][seen::cKM].s[k] / sReg[j][seen::cSum].s[k] : sTot[j][seen::cKM].S / sTot[j][seen::cSum].S;
                const double pu = k < seen::kNReg ? purReg[j][k].eps : pur[j];
                T("%-4s %-28s %9.3f +- %-8.3f %12.3f | %10.3f %10.3f\n", kJTag[j], seen::kReg[k].name, ratC[j][k], ratCE[j][k], ratCA[j][k], pu, mp);
            }
        T("# reading (review finding 6): 0.56-0.70 and 0.70-1.00 are multi-pion background, which SMASH does not describe; their ratios are listed,\n");
        T("#   not used as a validation; the ratios move with the ghost treatment, the vertex slice and L_data (unconfirmed). Below the bump the ratio\n");
        T("#   near 1 is a coincidence of a too large simulated bump tail and a too small eta-related background, not a check of the eta yield.\n");
    }
    if (gHaveC) {
        T("\n## 3. Fit factors alpha (template fits of data N - CB, M4 0.40-0.75 or 0.43-0.70; [table] tables/eta_enriched_fit.txt; E1/E2 M1 recomputed)\n");
        T("# alpha = 1: SMASH at L_data = 6080.5 nb^-1 (UNCONFIRMED). M1 = bump, peak (direct + conversion, SMASH ratio), physical background.\n");
        T("# POOR FITS (review finding 5): errors below are scaled by sqrt(chi2/ndf); the model spread (next lines) is far larger; the factors are\n");
        T("#   comparisons of shapes, NOT yields; P_main data/sim: see section 6 for all methods\n");
        T("%-6s | %-26s %-26s %-26s %-12s %-9s | %-26s %-26s\n", "sel", "M1 a_bump", "M1 a_peak", "M1 a_bkg", "chi2/ndf", "p", "M1 (a) a_bkg", "P_main data/sim (5.4 MeV)");
        for (int is = 0; is < seen::kNSel; ++is) {
            const FitRow* f = FindFit(gFits, seen::kSel[is].tag, 0, false, false, "M1");
            const FitRow* fa = FindFit(gFits, seen::kSel[is].tag, 0, false, true, "M1");
            if (!f || !fa) continue;
            const double sf = std::sqrt(std::max(1., f->chi2 / f->ndf)), sfa = std::sqrt(std::max(1., fa->chi2 / fa->ndf));
            T("%-6s | %8.3f +- %-14.3f %8.3f +- %-14.3f %8.3f +- %-14.3f %6.1f/%-5d %-9.2g | %8.3f +- %-14.3f %8.3f +- %-14.3f\n", seen::kSel[is].tag, f->Par("a_bump").first,
              f->Par("a_bump").second * sf, f->Par("a_peak").first, f->Par("a_peak").second * sf, f->Par("a_bkg").first, f->Par("a_bkg").second * sf, f->chi2, f->ndf,
              TMath::Prob(f->chi2, f->ndf), fa->Par("a_bkg").first, fa->Par("a_bkg").second * sfa, gPM[is].r54.first, gPM[is].r54.second);
        }
        // ranges over everything
        const char* pn[3] = {"a_bump", "a_peak", "a_bkg"};
        for (const char* model : {"M1", "M3"}) {
            for (int k = 0; k < 3; ++k) {
                if (std::string(model) == "M3" && k == 2) continue;
                double lo = 1e9, hi = -1e9;
                for (const FitRow& f : gFits) {
                    if (f.model != model) continue;
                    const VE v = f.Par(pn[k]);
                    if (v.first == v.first) { lo = std::min(lo, v.first); hi = std::max(hi, v.first); }
                }
                T("range of %s %s over 10 selections x 2 ranges x 2 scales x 2 ghost variants: %.3f - %.3f\n", model, pn[k], lo, hi);
            }
        }
        double lo = 1e9, hi = -1e9, plo = 1e9, phi = -1e9;
        for (const FitRow& f : gFits) if (f.model == "M1") {
            lo = std::min(lo, f.chi2 / f.ndf); hi = std::max(hi, f.chi2 / f.ndf);
            if (f.rng == 0) { plo = std::min(plo, TMath::Prob(f.chi2, f.ndf)); phi = std::max(phi, TMath::Prob(f.chi2, f.ndf)); }
        }
        T("range of M1 chi2/ndf: %.2f - %.2f; fit probability in 0.40-0.75: %.1e - %.1e\n", lo, hi, plo, phi);
        {
            double bl = 1e9, bh = -1e9;
            for (const FitRow& f : gFits) if (f.model == "M1" || f.model == "M3") { const VE v = f.Par("a_bump"); if (v.first == v.first) { bl = std::min(bl, v.first); bh = std::max(bh, v.first); } }
            T("reading: the simulated bump is too large relative to the simulated background at the bump edge; the fit factor lies between about %.2f and %.2f\n", bl, bh);
            T("  depending on model (M1, M3), fit range, mass scale and ghost variant (poor fits, not a yield; the earlier word 'stable' is withdrawn).\n");
            T("  M3 gives the eta-related background mostly > 1 (range below) while the bump (same etas) < 1: contradictory; never combine bump and peak\n");
            T("  factors into an eta factor. The P_main data/sim below carries statistical errors only; with the fit-model spread: section 6.\n");
        }
        for (const char* pq : {"a_bkgEta", "a_bkgMP"}) {
            double l2 = 1e9, h2 = -1e9;
            for (const FitRow& f : gFits) if (f.model == "M3") { const VE v = f.Par(pq); l2 = std::min(l2, v.first); h2 = std::max(h2, v.first); }
            T("range of M3 %s: %.3f - %.3f\n", pq, l2, h2);
        }
        double pl = 1e9, ph = -1e9;
        for (int is = 0; is < seen::kNSel; ++is) if (gPM[is].ok) { pl = std::min(pl, gPM[is].r54.first); ph = std::max(ph, gPM[is].r54.first); }
        T("range of the P_main peak ratio data/sim (5.4 MeV, sim M4 x 0.9952): %.3f - %.3f; without OA cut (E2, E3, E4, E2mm, E2v) and with OA > 2 (E1, E3oa, E4oa, E1mm, E1v) differ:\n", pl, ph);
        for (int is = 0; is < seen::kNSel; ++is)
            if (gPM[is].ok) T("  %-5s data P_main %5.0f +- %-4.0f  sim %5.0f +- %-5.0f  data/sim %.3f +- %.3f  (a) %.3f +- %.3f\n", seen::kSel[is].tag, gPM[is].data.first,
                              gPM[is].data.second, gPM[is].sim54.first, gPM[is].sim54.second, gPM[is].r54.first, gPM[is].r54.second, gPM[is].rA.first, gPM[is].rA.second);
    }
    if (gHaveD) {
        T("\n## 4. Peak content per class (P_main fits of data N - CB per class; simulation in expected data rows; all M4)\n");
        for (const ClassRes& o : gCls) {
            T("\n# %s, %s classes\n", o.sel.c_str(), o.cs == "OA5" ? "OA [deg]" : "M(ee) [MeV/c2]");
            T("%-8s %16s %16s %16s %10s %16s | %15s %15s %15s %15s\n", "class", "sim direct", "sim conversion", "sim total", "QED dir", "data", "frac data",
              "frac sim tot", "frac sim dir", "frac QED dir");
            double tD = 0, tT = 0, tS = 0, tK = 0, vD = 0, vT = 0;
            for (size_t i = 0; i < o.name.size(); ++i) { tD += o.dS[i]; tT += o.st[i]; tS += o.sd[i]; tK += o.sk[i]; vD += o.dSe[i] * o.dSe[i]; vT += o.ste[i] * o.ste[i]; }
            for (size_t i = 0; i < o.name.size(); ++i)
                T("%-8s %7.0f +- %-6.0f %7.0f +- %-6.0f %7.0f +- %-6.0f %10.0f %7.0f +- %-6.0f | %6.3f +- %-5.3f %15.3f %15.3f %15.3f\n", o.name[i].c_str(), o.sd[i], o.sde[i],
                  o.sc[i], o.sce[i], o.st[i], o.ste[i], o.sk[i], o.dS[i], o.dSe[i], Frac(o.dS, o.dSe, i).first, Frac(o.dS, o.dSe, i).second, o.st[i] / tT, o.sd[i] / tS, o.sk[i] / tK);
            T("sum      %16s %16s %7.0f +- %-6.0f %10s %7.0f +- %-6.0f   (data P_main without classes: %.0f +- %.0f); data/sim %.3f +- %.3f\n", "", "", tT, std::sqrt(vT), "", tD,
              std::sqrt(vD), o.dTot, o.dTotE, tD / tT, (tD / tT) * std::sqrt(vD / (tD * tD) + vT / (tT * tT)));
            T("data errors: stat (+) model; model / stat per class:");
            for (size_t i = 0; i < o.name.size(); ++i) T(" %s %.0f/%.0f", o.name[i].c_str(), o.dSm[i], o.dSs[i]);
            T("\n");
            T("one factor on the simulated total: %.3f +- %.3f, chi2/ndf %.1f/%d (p %.2g); stat-only errors (before the review): chi2/ndf %.1f/%d\n", o.mix.a1, o.mix.ea1,
              o.mix.chi21, o.mix.ndf1, TMath::Prob(o.mix.chi21, o.mix.ndf1), o.mixStat.chi21, o.mixStat.ndf1);
            T("free conversion / direct factors: %.3f +- %.3f / %.3f +- %.3f, chi2/ndf %.1f/%d%s; stat-only: %.3f +- %.3f / %.3f +- %.3f\n", o.mix.a, o.mix.ea, o.mix.b,
              o.mix.eb, o.mix.chi2, o.mix.ndf, o.mix.a < 0 ? "  (negative conversion factor: unphysical)" : "", o.mixStat.a, o.mixStat.ea, o.mixStat.b, o.mixStat.eb);
            T("(QED-reshaped direct decay: illustration of shapes only, column 'QED dir'; no fit result quoted)\n");
        }
    }
    if (gHaveE) {
        T("\n## 5. Efficiency of single cuts (simulation: base pre, all M4; data peak: pass/fail P_main fits at base pre; data bump [table] cuts_compare_data.txt)\n");
        T("# bump data base: no vz cut (1D scans); pT4 and vertex: bump start base (pre without vz cut + oa_lab < 40 + 0.30 < M(pipi) < 0.41), sim at the same base\n");
        T("# (s) = data peak value biased by sculpting (data study); data peak errors = stat (+) model (per-variant efficiency, half range)\n");
        T("# sim direct: generator M(ee) model not QED; 'total QED' = ILLUSTRATION only (flat acceptance and fixed content assumed), not a result\n");
        T("%-26s | %13s %13s %13s %13s %13s | %13s %13s %13s %13s %13s\n", "cut", "sim direct", "sim conv.", "sim total", "total QED", "data peak", "sim bump",
          "sim tagged", "sim S_bump", "data tagged", "data edge");
        for (int k = 0; k < kNCut; ++k) {
            const EffOut& o = gEff[k];
            const EffRes& sb = kCut[k].bumpStart ? o.sBumpStart : o.s[eBump];
            const EffRes& stg = kCut[k].bumpStart ? o.sTagStart : o.s[eTag];
            auto f = [](VE v) { return v.first == v.first ? std::string(Form("%5.3f+-%5.3f", v.first, v.second)) : std::string("      -      "); };
            T("%-26s | %13s %13s %13s %13s %13s | %13s %13s %13s %13s %13s\n", (std::string(kCut[k].txt) + (kCut[k].sculpt ? " (s)" : "")).c_str(),
              f(VE(o.s[eDir].eps, o.s[eDir].err)).c_str(), f(VE(o.s[eConv].eps, o.s[eConv].err)).c_str(), f(VE(o.s[eTot].eps, o.s[eTot].err)).c_str(),
              f(VE(o.s[eTotKW].eps, o.s[eTotKW].err)).c_str(), f(o.dPeak).c_str(), f(VE(sb.eps, sb.err)).c_str(), f(VE(stg.eps, stg.err)).c_str(),
              f(kCut[k].bumpStart ? VE(kNaN, kNaN) : VE(o.s[eSbu].eps, o.s[eSbu].err)).c_str(), f(o.dTag).c_str(), f(o.dEdge).c_str());
        }
    }
    if (gHaveE) {
        for (int k = 0; k < kNCut; ++k) {
            if (gEff[k].degTag) T("# %s: data tagged value has zero error in cuts_compare_data.txt (degenerate fit): not used\n", kCut[k].txt);
            if (gEff[k].degEdge) T("# %s: data edge-step value has zero error in cuts_compare_data.txt (degenerate fit): not used\n", kCut[k].txt);
        }
    }
    if (gHaveF) {
        T("\n## 6. eta peak, data / simulation by method (review finding 1; tables/eta_summary_peak_methods.txt, figure f)\n");
        T("# A all-M4 content; Aw window content; B P_main on data and on the sim peak template; C P_main on data and on the sim total (C(a): variant (a));\n");
        T("# D template fit M1 alpha_peak (error x sqrt(chi2/ndf); range over ranges, scales, ghost variants); E side bands (range: left / right band alone)\n");
        T("# errors: +- stat +- model; 'tot' adds the absolute scale incl. L_data (20 %%); ghost: (b) = sim peak without ghost-lepton rows, (a) = variant (a)\n");
        T("%-9s | %-26s | %-26s | %-26s | %-26s | %-24s | %-24s\n", "sel", "A", "Aw", "B", "C", "D [range]", "E [range]");
        for (int is = 0; is < kNMs; ++is) {
            const MethRes& m = gMeth[is];
            T("%-9s | %5.2f +- %4.2f +- %4.2f (%4.2f) | %5.2f +- %4.2f +- %4.2f (%4.2f) | %5.2f +- %4.2f +- %4.2f (%4.2f) | %5.2f +- %4.2f +- %4.2f (%4.2f) | %4.2f +- %4.2f [%4.2f-%4.2f] | %4.2f +- %4.2f [%4.2f-%4.2f]\n",
              kMsTag[is], m.A.first, m.A.second, m.mA, m.Ab, m.Aw.first, m.Aw.second, m.mAw, m.Awb, m.B.first, m.B.second, m.mB, m.Bb, m.C.first, m.C.second, m.mC,
              m.Ca.first, m.D.first, m.Dscaled, m.Dlo, m.Dhi, m.E.first, m.E.second, m.Elo, m.Ehi);
        }
        T("# in brackets after A, Aw, B: ghost variant (b); after C: variant (a). Data P_main [rows] (stat, model):");
        for (int is = 0; is < kNMs; ++is) T("  %s %.0f +- %.0f +- %.0f", kMsTag[is], gMeth[is].dS, gMeth[is].dSe, gMeth[is].dSm);
        T("\n");
    }
    T("\n## UNVERIFIED / caveats\n");
    T("- L_data = 6080.5 nb^-1 (producers' note): every absolute ratio and alpha scales with 6080.5 / L_data; absolute scale known to about 15-20 %%.\n");
    T("- The eta-set generator M(ee) model of eta -> pi+pi-e+e- (non-QED on accepted rows; the PLUTO model and generator-level rows are needed).\n");
    T("- The QED (Kroll-Wada) reshaping assumes an acceptance flat in the true M(ee) and a fixed content; both arbitrary (illustration only).\n");
    T("- The eta-related / multi-pion split uses geninfo and two generations of parent ids (sim_eta_enriched.C).\n");
    T("- The K0S mass scale (sim/data 1.0048) is transferred from M(pipi) to M4 in the scaled fits and in P_main.\n");
    {
        double l0 = 1e9, h0 = -1e9, l1 = 1e9, h1 = -1e9;
        for (const FitRow& f : gFits) {
            if (f.model != "M1") continue;
            const double c = f.chi2 / f.ndf;
            if (f.rng == 0) { l0 = std::min(l0, c); h0 = std::max(h0, c); } else { l1 = std::min(l1, c); h1 = std::max(h1, c); }
        }
        T("- Template fits neglect bin-to-bin template correlations; chi2/ndf of M1 is %.2f-%.2f in 0.40-0.75 and %.2f-%.2f in 0.43-0.70: the fitted\n"
          "  templates do not describe the data shape, the alphas are a comparison, not measured scale factors.\n", l0, h0, l1, h1);
    }
    T("- Data bump efficiencies come from the data study (tagged eta, edge step) at a base without the vz cut.\n");
    T("- Ghost-flagged true pairs: which treatment rule 1 means is for the user to confirm (variant (a) shown everywhere).\n");
    fclose(gT); gT = nullptr;
}

}  // namespace ses

// ============================================================================
void sim_eta_summary(int parts = 31) {
    using namespace simeng;
    using namespace ses;
    TStopwatch swAll;
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0); gStyle->SetEndErrorSize(0);
    gStyle->SetTextFont(42); gStyle->SetLabelFont(42, "xyz"); gStyle->SetTitleFont(42, "xyz"); gStyle->SetLegendFont(42);
    gStyle->SetFrameLineWidth(1); gStyle->SetLineWidth(1);
    TH1::AddDirectory(kFALSE);
    gC = fopen(Form("%s/eta_summary_checks.txt", kTabDir), "w");
    C("# eta_summary_checks: numbers recomputed by sim_eta_summary.C against the tables of the study macros\n");
    C("# 'this' = this macro, 'table' = the study table; OK when |diff| is below the printing precision of the table\n");

    // MM4 cut of the enriched study (only E1mm / E2mm use it; read for completeness)
    {
        std::string l;
        for (const std::string& x : Lines("eta_enriched_mm4.txt"))
            if (x.compare(0, 14, "CHOSEN MM4 cut") == 0) l = x;
        double v = 2.85;
        if (!l.empty()) sscanf(l.c_str(), "CHOSEN MM4 cut: MM4 < %lf", &v);
        seen::gMMcut = v;
        C("MM4 cut of tables/eta_enriched_mm4.txt: %.2f\n", v);
    }
    // ---- load
    SimData* d = new SimData;
    if (!Load(*d, kCompact, kMaskUS, seen::kLoadLo, seen::kLoadHi, 1)) return;
    gD = d;
    seen::gD = d;
    {
        TFile* f = TFile::Open(seen::kCache, "READ");
        const bool have = f && !f->IsZombie() && f->Get("er_std") && f->Get("er_lep") && f->Get("er_eta");
        if (f) { f->Close(); delete f; }
        if (!have) { printf("eta-related cache missing: building it\n"); if (!seen::BuildCache(gC)) { printf("cache build FAILED\n"); return; } }
    }
    if (!seen::AttachER(*d, gC)) { printf("eta-related flags not aligned: STOP\n"); return; }
    scan::Data* dd = new scan::Data;
    if (!scan::Load(*dd, scan::kCompact, seen::kLoadLo, seen::kLoadHi)) return;
    gDD = dd;
    // reduced copy with the signal classes (template 1, classes 2-6) for parts d and e
    SimData* ds = new SimData;
    for (int k = 0; k < kNTreeAll; ++k)
        for (const Row& r : d->t[k]) if (r.tmpl == 1 && r.cls >= 2 && r.cls <= 6) ds->t[k].push_back(r);
    ds->loaded = d->loaded;
    gDS = ds;
    {
        size_t n = 0; for (int k = 0; k < kNTreeAll; ++k) n += ds->t[k].size();
        C("reduced copy (template 1, classes 2-6): %zu rows\n", n);
        // the eta set holds no bump rows (no eta -> pi+pi-pi0 there): owned rows of classes 2/3 per set
        double w[4] = {0, 0, 0, 0};
        for (int s = 1; s <= 3; ++s)
            for (const Row& r : ds->tree(s, kUS)) if ((r.cls == 2 || r.cls == 3) && Owned(r, gOpt)) w[s] += r.w;
        C("bump classes 2+3, owned sum w per set: std %.6g, lep %.6g, eta %.6g\n", w[1], w[2], w[3]);
    }

    if (parts & 1) { FillAB(); FigA(); FigB(); }
    if (parts & 2) { PartC(); FigC(); }
    if (parts & 4) { PartD(); FigD(); }
    if (parts & 8) { PartE(); FigE(); }
    if (parts & 16) { PartF(); FigF(); }
    WriteSummary();
    C("\nworst relative difference over all checks: %.2e; checks outside the tolerance: %.0f\n", gWorstRel, gWorstPull);
    fclose(gC); gC = nullptr;
    printf("sim_eta_summary(%d): %.0f s\n", parts, swAll.RealTime());
}
