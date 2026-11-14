// ============================================================================
// sim_review.C  --  corrections after the adversarial review of the study
// ============================================================================
// HADES p+p, T = 4.53 GeV, pi+ pi- e+ e-. Simulation: SMASH std + lep, owned
// rows (mode A default), w (PDG), sigma = sum w / L_sim (221.05 nb^-1), event-
// group variance, double precision. Data: (N - CB) / L_data, CB = 2 sqrt(N++
// N--) in 10 MeV M4 bins (no k), L_data = 6080.5 nb^-1 UNCONFIRMED.
// eta -> pi+pi-e+e- (direct decay) is NOT in the simulation.
//
// Review findings treated here (see README, section "Review corrections"):
//   1  ghost-flagged true e+e- pairs: in T1 when signal, in T3 when background.
//      Variants (sim_ghostvar.h): (a) T1 + T2 + ghost-flagged true pairs,
//      (b) T1 without ghost rows + T2. Decomposition of T3, S/B per variant.
//   2  R (data / simulation) as a table of variants: M4 window x vertex-z
//      slice x ghost treatment x mode / weight (L_data scales all by one factor)
//   3  vertex z distributions of data and simulation; every ratio also in the
//      central target slice -140 < vz < -100 and the edge slices
//   5  OA(e+e-) and M(e+e-) classes with the ghost variants, mode B, central vz;
//      cumulative OA > x
//   6  template-fit variants: ghost-pair background, no-ghost signal, central vz
//   8  position and width of the simulated conversion peak against the fit
//      choices (binning, range, ghost rows)
// Parts (bit mask): 1 ghost + R + regions + S/B + M4 figures, 2 vertex z,
//   4 OA / M(ee) classes, 8 peak fits, 16 template-fit variants.
//
// Figures: figs/review_*.png (.pdf); tables: tables/review_*.txt
// Run from /home/przygoda/HADES/FAT/FAT:
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/aclic_review", kTRUE)' 'plots/eta_sim/sim_review.C+(31)'
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TStopwatch.h"
#include "TStyle.h"

#include "SimEngine.h"
#include "sim_ghostvar.h"
#include "sim_datacomp.C"      // dcmp: data filler with group variance, ratio, chi2, plot helpers (not modified)
#include "sim_templatefit.C"   // stf: templates, data covariance, fit, implied contents

namespace srev {
using namespace simeng;
using dcmp::Nm;
using dcmp::P2;
using dcmp::Txt;

const char* kXM4 = "M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}]";
const std::string kNoteL = "L_{data} = 6080.5 nb^{-1}, unconfirmed: every data value scales with 1/L_{data}";
const std::string kNoteS = "simulation: SMASH std + lep, owned rows (mode A unless stated), w (PDG); #sigma = #Sigma w / L_{sim}, L_{sim} = 221.05 nb^{-1}; no efficiency correction";
const std::string kNoteD = "#eta #rightarrow #pi^{+}#pi^{-}e^{+}e^{-} (direct decay) is NOT in the simulation (eta set not available)";
const std::string kNoteC = "data: CB = 2#sqrt{N_{++}N_{--}} in 10 MeV M_{4} bins (no k factor); data N black, CB red, N - CB blue";

// ---------------------------------------------------------------- windows of M4
struct Win { const char* name; double a, b; };
const std::vector<Win> kRegs = {{"below the bump", 0.30, 0.455},      {"bump window", scan::kBumpLo, scan::kBumpHi},
                                {"gap bump-peak", scan::kBumpHi, scan::kPeakLo}, {"peak window", scan::kPeakLo, scan::kPeakHi},
                                {"0.554-0.58", scan::kPeakHi, 0.58},  {"0.56-0.70", 0.56, 0.70},
                                {"0.70-0.90", 0.70, 0.90},            {"0.90-1.10", 0.90, 1.10},
                                {"R 0.60-0.80", 0.60, 0.80},          {"R 0.56-0.62", 0.56, 0.62},
                                {"R 0.65-0.85", 0.65, 0.85},          {"all 0.30-1.10", 0.30, 1.10}};
const int kIR = 8;        // index of the default R window
const int kNRwin = 4;     // R windows: 8, 9, 10, 11
inline std::vector<std::pair<double, double>> WinList(const std::vector<Win>& w) {
    std::vector<std::pair<double, double>> o;
    for (const Win& x : w) o.push_back(std::make_pair(x.a, x.b));
    return o;
}

// ---------------------------------------------------------------- simulated window sums
// sums of several components in several M4 windows (row values in [a, b)), exact event-group variance
struct WT {
    int nc = 0, nw = 0;
    std::vector<double> S, V, N;
    double s(int c, int w) const { return S[c * nw + w]; }
    double v(int c, int w) const { return V[c * nw + w]; }
    double xs(int c, int w) const { return S[c * nw + w] / kLsim; }
    double xe(int c, int w) const { return std::sqrt(V[c * nw + w]) / kLsim; }
    double neff(int c, int w) const { return V[c * nw + w] > 0 ? S[c * nw + w] * S[c * nw + w] / V[c * nw + w] : 0; }
};
WT SimWin(const SimData& d, const Sel& sel, const std::vector<Comp>& comps, const std::vector<std::pair<double, double>>& wins,
          const SimOpt& o = SimOpt(), const Var& x = Var()) {
    WT t; t.nc = (int)comps.size(); t.nw = (int)wins.size();
    t.S.assign(t.nc * t.nw, 0.); t.V.assign(t.nc * t.nw, 0.); t.N.assign(t.nc * t.nw, 0.);
    std::vector<double> g(t.nc * t.nw, 0.);
    std::vector<int> touched;
    auto flush = [&]() {
        for (int i : touched) if (g[i] != 0) { t.V[i] += g[i] * g[i]; g[i] = 0; }
        touched.clear();
    };
    for (int set = 1; set <= 2; ++set) {
        if (!(o.sets & set)) continue;
        const std::vector<Row>& v = d.tree(set, o.hyp);
        int cur = -1;
        for (const Row& r : v) {
            if (r.grp != cur) { flush(); cur = r.grp; }
            if (!Owned(r, o) || !sel(r)) continue;
            const double w = o.db ? r.w_db : r.w;
            const double m = x ? x(r) : (double)r.m4;
            if (!(m == m)) continue;
            for (int c = 0; c < t.nc; ++c) {
                if (!comps[c].f(r)) continue;
                for (int k = 0; k < t.nw; ++k) {
                    if (m < wins[k].first || m >= wins[k].second) continue;
                    const int i = c * t.nw + k;
                    g[i] += w; t.S[i] += w; t.N[i] += 1; touched.push_back(i);
                }
            }
        }
        flush();
    }
    return t;
}
struct RV { double r = 0, e = 0; };
inline RV Div(double a, double ea, double b, double eb) {
    RV q;
    if (b > 0 && a != 0) { q.r = a / b; q.e = std::fabs(q.r) * std::sqrt(ea * ea / (a * a) + eb * eb / (b * b)); }
    return q;
}
inline RV DataSim(const scan::WinSum& w, const WT& t, int c, int k) {
    return Div(w.sig() / kLdata, w.sigErr() / kLdata, t.xs(c, k), t.xe(c, k));
}

// selections: vertex-z slice iv (sim_ghostvar.h) and an extra cut (0 none, 1 OA > 2)
inline Sel SSel(int iv, int e = 0) { return [iv, e](const Row& r) { return PassVz(r, iv) && (e == 0 || r.oa > 2.f); }; }
inline scan::Sel DSel(int iv, int e = 0) { return [iv, e](const scan::Row& r) { return PassVz(r, iv) && (e == 0 || r.oa > 2.f); }; }
inline SimOpt OptV(int i) { SimOpt o; if (i == 1) o.mode = 1; if (i == 2) o.db = true; return o; }
const char* kOptName[3] = {"A, w", "B, w", "A, w_db"};

inline TCanvas* Canvas(const char* name, int w, int h) { TCanvas* c = new TCanvas(name, "", w, h); c->SetFillColor(0); return c; }
inline TPad* Pad(double x1, double y1, double x2, double y2, bool logy, double lm, double rm, double tm, double bm) {
    TPad* p = new TPad(Nm("pd"), "", x1, y1, x2, y2);
    p->SetLeftMargin(lm); p->SetRightMargin(rm); p->SetTopMargin(tm); p->SetBottomMargin(bm);
    p->SetLogy(logy); p->SetTicks(1, 1); p->SetFillStyle(0); p->Draw(); p->cd();
    return p;
}
inline TH1D* Fr(double xlo, double xhi, double ylo, double yhi, const char* xt, const char* yt, double lab, double tit, double yoff = 1.2,
                double xoff = 1.0) {
    TH1D* f = new TH1D(Nm("fr"), "", 100, xlo, xhi);
    f->SetDirectory(nullptr); f->SetMinimum(ylo); f->SetMaximum(yhi); f->SetStats(0);
    f->GetXaxis()->SetTitle(xt); f->GetYaxis()->SetTitle(yt);
    f->GetXaxis()->SetLabelSize(lab); f->GetYaxis()->SetLabelSize(lab);
    f->GetXaxis()->SetTitleSize(tit); f->GetYaxis()->SetTitleSize(tit);
    f->GetYaxis()->SetTitleOffset(yoff); f->GetXaxis()->SetTitleOffset(xoff);
    f->GetYaxis()->SetNdivisions(505); f->GetYaxis()->SetMaxDigits(4);
    f->Draw("AXIS");
    return f;
}
inline void Hdr(TCanvas* c, const std::vector<std::pair<std::string, int>>& lines, double top = 0.985, double dy = 0.022, double sz = 0.016) {
    c->cd();
    double y = top;
    for (size_t i = 0; i < lines.size(); ++i) {
        TLatex t; t.SetNDC(); t.SetTextFont(i == 0 ? 62 : 42); t.SetTextSize(i == 0 ? sz * 1.15 : sz); t.SetTextColor(lines[i].second);
        t.SetTextAlign(13); t.DrawLatex(0.01, y, lines[i].first.c_str());
        y -= (i == 0 ? 1.2 : 1.0) * dy;
    }
}
inline void Save(TCanvas* c, const char* base) {
    c->SaveAs(Form("%s/%s.png", kFigDir, base));
    c->SaveAs(Form("%s/%s.pdf", kFigDir, base));
}
// histogram of a data DH (dcmp) in nb/unit; what: 0 N, 3 CB, 4 N - CB
inline TH1D* DHist(const dcmp::DH& h, int what) { return h.Hist(what, Nm("dh"), 2); }

// ============================================================================
// Part 1: ghost-flagged true pairs, R variants, regions, S/B, M4 figures
// ============================================================================
void Part1(const SimData& sd, const scan::Data& dd) {
    TStopwatch sw;
    // ---------------- (1) decomposition of T3 and the ghost rows, pre, mode A, w ----------------
    {
        FILE* f = dcmp::OpenTab("review_ghost.txt");
        P2(f, "# Review finding 1: ghost-flagged true e+e- pairs. Selection pre, std + lep owned rows, mode A, w (PDG); sigma [nb] = sum w / %.2f\n", kLsim);
        P2(f, "# src/sim_truth.h TemplateClass: cls 2-6 -> T1 before the ghost check; any other row with a ghost-flagged lepton -> T3 'fake lepton'.\n");
        P2(f, "# 'ghost-flagged true pair' = T3, ghost lepton flag, lepton-pair class L1-L6 (one true source), not the same particle twice.\n");
        const std::vector<Comp> dc = {
            CompTmpl(1),
            {"T1 rows with a ghost lepton", [](const Row& r) { return r.tmpl == 1 && GhostLepton(r); }},
            CompTmpl(2),
            CompTmpl(3),
            CompGhostPair(),
            {"T3 ghost + L1-L6 incl. same particle", [](const Row& r) { return r.tmpl == 3 && GhostLepton(r) && r.lcls >= 1 && r.lcls <= 6; }},
            {"T3 same particle twice", [](const Row& r) { return r.tmpl == 3 && (r.flags & simtruth::kFlagSameParticle); }},
            {"T3 misidentified (L0), not same", [](const Row& r) { return r.tmpl == 3 && r.lcls == 0 && !(r.flags & simtruth::kFlagSameParticle); }},
            {"T3 other (ghost, L7-L9)", [](const Row& r) { return r.tmpl == 3 && !GhostPair(r) && !(r.flags & simtruth::kFlagSameParticle) && r.lcls != 0; }},
            CompTmpl(4),
            CompSimVar(0), CompSimVar(1), CompSimVar(2)};
        const std::vector<Win> w = {{"all M4", -1e9, 1e9}, {"0.30-1.10", 0.30, 1.10}, {"0.60-0.80", 0.60, 0.80},
                                    {"bump window", scan::kBumpLo, scan::kBumpHi}, {"peak window", scan::kPeakLo, scan::kPeakHi}};
        const WT t = SimWin(sd, SelPreset(kPre), dc, WinList(w));
        P2(f, "%-40s", "component");
        for (const Win& x : w) P2(f, " %24s", x.name);
        P2(f, "\n");
        for (int c = 0; c < t.nc; ++c) {
            P2(f, "%-40s", dc[c].name.c_str());
            for (int k = 0; k < t.nw; ++k) P2(f, " %9.4g +- %-6.2g (%5.0f)", t.xs(c, k), t.xe(c, k), t.neff(c, k));
            P2(f, "\n");
        }
        P2(f, "# (numbers: sigma +- error [nb] (N_eff))\n");
        P2(f, "share of T3 that is ghost-flagged true pairs: all M4 %.3f, 0.60-0.80 %.3f\n", t.s(4, 0) / t.s(3, 0), t.s(4, 2) / t.s(3, 2));
        P2(f, "share of T1 with a ghost lepton: all M4 %.3f\n", t.s(1, 0) / t.s(0, 0));
        // ghost share of one-source pairs per set
        const Comp one = {"one-source pairs", [](const Row& r) {
                              return r.tmpl <= 3 && r.lcls >= 1 && r.lcls <= 6 && !(r.flags & simtruth::kFlagSameParticle); }};
        const Comp oneG = {"one-source pairs with ghost", [one](const Row& r) { return one.f(r) && GhostLepton(r); }};
        for (int s = 1; s <= 2; ++s) {
            SimOpt o; o.sets = s;
            const WT q = SimWin(sd, SelPreset(kPre), {one, oneG}, {{-1e9, 1e9}, {0.60, 0.80}}, o);
            P2(f, "ghost share of true one-source pairs (T1-T3, L1-L6), %s-owned rows: all M4 %.4f, 0.60-0.80 %.4f\n", s == 1 ? "std" : "lep",
               q.s(1, 0) / q.s(0, 0), q.s(1, 1) / q.s(0, 1));
        }
        // groups with three or more lepton tracks among close pairs
        {
            const Sel s = [](const Row& r) { return PassPre(r) && r.oa < 2.f && r.m4 >= 0.30f && r.m4 < 1.10f; };
            const std::vector<Comp> cc = {CompGhostPair(), {"gp3", [](const Row& r) { return GhostPair(r) && r.n_lep >= 3; }},
                                          CompSimVar(0), {"t12_3", [](const Row& r) { return InSimVar(r, 0) && r.n_lep >= 3; }}};
            const WT q = SimWin(sd, s, cc, {{-1e9, 1e9}});
            double nA = 0, n3 = 0;
            for (const scan::Row& r : dd.s[0]) {
                if (!r.cor_ok || !PassPre(r) || !(r.oa < 2.f) || !(r.m4 >= 0.30f && r.m4 < 1.10f)) continue;
                nA += 1; if (r.n_lep >= 3) n3 += 1;
            }
            P2(f, "rows with >= 3 lepton tracks in the event, OA < 2 deg, M4 0.30-1.10, pre: ghost-flagged true pairs %.4f, T1+T2 %.4f, data (all e+e- rows) %.4f\n",
               q.s(1, 0) / q.s(0, 0), q.s(3, 0) / q.s(2, 0), n3 / nA);
        }
        if (f) fclose(f);
        printf("Part1 (1) %.0f s\n", sw.RealTime()); sw.Continue();
    }

    // ---------------- (2) R and regions: vz slice x extra cut x option x ghost variant ----------------
    const std::vector<std::pair<double, double>> wl = WinList(kRegs);
    std::vector<scan::Sel> ds;
    for (int e = 0; e < 2; ++e) for (int iv = 0; iv < 5; ++iv) ds.push_back(DSel(iv, e));
    scan::FillOpt fo; fo.nbin = 320; fo.lo = 0.30; fo.hi = 1.10; fo.reb = 4; fo.windows = wl;
    const std::vector<scan::Spectra> dsp = scan::FillMany(dd, ds, fo);
    auto DW = [&](int e, int iv, int k) -> const scan::WinSum& { return dsp[e * 5 + iv].win[2 + k]; };
    // components: 0-2 sim variants; 3 S_bump(def), 4 S_bump(b); 5 S_peak(def), 6 S_peak(b); 7 B(def), 8 B(a); 9 ghost pairs
    const std::vector<Comp> comps = {CompSimVar(0), CompSimVar(1), CompSimVar(2), CompSigVar(0, 2, 3), CompSigVar(2, 2, 3),
                                     CompSigVar(0, 4, 4), CompSigVar(2, 4, 4), CompBkgVar(0), CompBkgVar(1), CompGhostPair()};
    std::vector<WT> st(2 * 5 * 3);   // [e][iv][opt]
    for (int e = 0; e < 2; ++e)
        for (int iv = 0; iv < 5; ++iv)
            for (int io = 0; io < 3; ++io) st[(e * 5 + iv) * 3 + io] = SimWin(sd, SSel(iv, e), comps, wl, OptV(io));
    auto ST = [&](int e, int iv, int io) -> const WT& { return st[(e * 5 + iv) * 3 + io]; };
    printf("Part1 (2) sums %.0f s\n", sw.RealTime()); sw.Continue();

    {
        FILE* f = dcmp::OpenTab("review_R_variants.txt");
        P2(f, "# Review finding 2: R = data (N - CB) / simulation in an M4 window, for every variant. L_data = %.1f nb^-1 UNCONFIRMED: all R scale with 1/L_data.\n", kLdata);
        P2(f, "# data: CB = 2 sqrt(N++ N--) in 10 MeV M4 bins, no k; event-group variances (ScanEngine windows). sim: std + lep owned rows, L_sim = %.2f nb^-1.\n", kLsim);
        P2(f, "# ghost variants: T12 = default T1 + T2; (a) = T1 + T2 + ghost-flagged true pairs; (b) = T1 without ghost rows + T2 (sim_ghostvar.h)\n");
        P2(f, "# vz slices [mm]: pre = vz > -160; central = -140..-100; edge = -160..-140; down = vz > -100; novz = z > -500 only. 'B, w_db' is not physical (closure) and omitted.\n");
        P2(f, "# errors: statistical (data and simulation); the std-lep covariance is neglected\n");
        for (int e = 0; e < 2; ++e) {
            P2(f, "\n## selection %s\n", e ? "pre + OA > 2 (same vz slices)" : "pre (isBest >= 0, PT3, cor_ok, vz slice)");
            P2(f, "%-8s %-14s %16s | %15s %15s %15s | %15s %15s %15s | %15s %15s %15s\n", "vz", "window", "data N-CB [nb]", "A,w: T12", "A,w: (a)",
               "A,w: (b)", "B,w: T12", "B,w: (a)", "B,w: (b)", "A,wdb: T12", "A,wdb: (a)", "A,wdb: (b)");
            for (int iv = 0; iv < 5; ++iv)
                for (int k : {8, 9, 10, 6, 11}) {
                    const scan::WinSum& w = DW(e, iv, k);
                    P2(f, "%-8s %-14s %8.4g +- %-5.2g |", kVzSl[iv].tag, kRegs[k].name, w.sig() / kLdata, w.sigErr() / kLdata);
                    for (int io = 0; io < 3; ++io) {
                        for (int v = 0; v < 3; ++v) { const RV q = DataSim(w, ST(e, iv, io), v, k); P2(f, " %6.3f +- %-5.3f", q.r, q.e); }
                        P2(f, io < 2 ? " |" : "\n");
                    }
                }
        }
        // overall spread
        double rmin = 1e9, rmax = -1e9;
        for (int iv = 0; iv < 4; ++iv) for (int k : {8, 9, 10}) for (int io = 0; io < 3; ++io) for (int v = 0; v < 3; ++v) {
            const RV q = DataSim(DW(0, iv, k), ST(0, iv, io), v, k);
            rmin = std::min(rmin, q.r); rmax = std::max(rmax, q.r);
        }
        P2(f, "\n# spread of R over pre/central/edge/down x windows 0.60-0.80, 0.56-0.62, 0.65-0.85 x options x ghost variants: %.3f - %.3f\n", rmin, rmax);
        {
            double a = 1e9, b = -1e9;
            for (int iv : {0, 1}) for (int io = 0; io < 3; ++io) for (int v = 0; v < 3; ++v) {
                const RV q = DataSim(DW(0, iv, 8), ST(0, iv, io), v, 8);
                a = std::min(a, q.r); b = std::max(b, q.r);
            }
            P2(f, "# spread of R in 0.60-0.80 over pre and central vz x options x ghost variants: %.3f - %.3f\n", a, b);
        }
        if (f) fclose(f);
    }
    {
        FILE* f = dcmp::OpenTab("review_regions.txt");
        P2(f, "# Review findings 1-4: data (N - CB) / simulation per M4 region, per vz slice and ghost variant; mode A, w (PDG)\n");
        P2(f, "# L_data = %.1f nb^-1 UNCONFIRMED; 'after R' = ratio / R(0.60-0.80) of the same slice and variant (error of R not included)\n", kLdata);
        P2(f, "# eta -> pi+pi-e+e- (direct) is NOT in the simulation; the peak-window ratio follows the background trend of the neighbouring regions\n");
        for (int e = 0; e < 2; ++e)
            for (int iv = 0; iv < 5; ++iv) {
                const WT& t = ST(e, iv, 0);
                const RV R0 = DataSim(DW(e, iv, kIR), t, 0, kIR), Ra = DataSim(DW(e, iv, kIR), t, 1, kIR), Rb = DataSim(DW(e, iv, kIR), t, 2, kIR);
                P2(f, "\n## %s, vz slice %s: R(0.60-0.80) T12 %.3f +- %.3f, (a) %.3f +- %.3f, (b) %.3f +- %.3f\n", e ? "pre+OA>2" : "pre", kVzSl[iv].tag,
                   R0.r, R0.e, Ra.r, Ra.e, Rb.r, Rb.e);
                P2(f, "%-15s %16s %16s %8s %15s %15s %15s %9s %9s %9s %18s\n", "region", "data N-CB [nb]", "sim T12 [nb]", "N_eff", "data/T12",
                   "data/(a)", "data/(b)", "T12 aft R", "(a) aft R", "(b) aft R", "ghost pairs [nb]");
                for (int k = 0; k < (int)kRegs.size(); ++k) {
                    const scan::WinSum& w = DW(e, iv, k);
                    const RV q0 = DataSim(w, t, 0, k), qa = DataSim(w, t, 1, k), qb = DataSim(w, t, 2, k);
                    P2(f, "%-15s %8.4g +- %-5.2g %8.4g +- %-5.2g %8.0f %6.3f +- %-5.3f %6.3f +- %-5.3f %6.3f +- %-5.3f %9.3f %9.3f %9.3f %9.4g +- %-5.2g\n",
                       kRegs[k].name, w.sig() / kLdata, w.sigErr() / kLdata, t.xs(0, k), t.xe(0, k), t.neff(0, k), q0.r, q0.e, qa.r, qa.e, qb.r, qb.e,
                       q0.r / R0.r, qa.r / Ra.r, qb.r / Rb.r, t.xs(9, k), t.xe(9, k));
                }
            }
        if (f) fclose(f);
    }
    {
        FILE* f = dcmp::OpenTab("review_sb.txt");
        P2(f, "# Review finding 1: simulated S/B in the data windows for the three ghost treatments; mode A, w (PDG), sigma [nb]\n");
        P2(f, "# S bump = signal classes 2+3 (bump window), S peak = class 4 (peak window 0.5336-0.5542); B = T2 (default, (b)) or T2 + ghost-flagged true pairs ((a))\n");
        P2(f, "%-9s %-6s %16s %16s %16s %16s | %16s %16s %16s\n", "vz slice", "window", "S (T1, all)", "S (b) no ghost", "B = T2", "B (a)", "S/B as used",
           "S/B (b) no ghosts", "S/B (a) ghosts both");
        for (int iv : {0, 1, 2}) {
            const WT& t = ST(0, iv, 0);
            for (int win : {1, 3}) {
                const int cs = win == 1 ? 3 : 5;
                const RV u = Div(t.xs(cs, win), t.xe(cs, win), t.xs(7, win), t.xe(7, win));
                const RV b = Div(t.xs(cs + 1, win), t.xe(cs + 1, win), t.xs(7, win), t.xe(7, win));
                const RV a = Div(t.xs(cs, win), t.xe(cs, win), t.xs(8, win), t.xe(8, win));
                P2(f, "%-9s %-6s %7.4g +- %-6.2g %7.4g +- %-6.2g %7.4g +- %-6.2g %7.4g +- %-6.2g | %7.3f +- %-6.3f %7.3f +- %-6.3f %7.3f +- %-6.3f\n",
                   kVzSl[iv].tag, win == 1 ? "bump" : "peak", t.xs(cs, win), t.xe(cs, win), t.xs(cs + 1, win), t.xe(cs + 1, win), t.xs(7, win),
                   t.xe(7, win), t.xs(8, win), t.xe(8, win), u.r, u.e, b.r, b.e, a.r, a.e);
            }
        }
        P2(f, "# (S/B errors neglect the S-B correlation; S and B are different rows, so it is small)\n");
        if (f) fclose(f);
    }

    // ---------------- figure: R variants ----------------
    {
        TCanvas* c = Canvas("cRv", 2200, 1150);
        Hdr(c, {{"R = data (N - CB) / simulation: variants of the M_{4} window, vertex-z slice, ghost treatment and mode / weight; selection pre", 1},
                {kNoteL + "; R is not determined to better than the spread shown here", kBlue + 1},
                {kNoteS, 1},
                {"ghost treatment: circles T1 + T2 (default), squares (a) T1 + T2 + ghost-flagged true pairs, triangles (b) T1 without ghost rows + T2", 1},
                {"left: R in M_{4} 0.60 - 0.80 per v_{z} slice (filled mode A, w; open mode B, w; crosses mode A, w_{db}); right: R per M_{4} window [GeV/c^{2}], mode A, w (filled pre, open central v_{z}, diamonds edge v_{z})", 1}},
            0.985, 0.026, 0.016);
        // left: R(0.60-0.80) vs vz slice; filled A,w; open B,w; crosses A,w_db
        c->cd();
        Pad(0.0, 0.0, 0.5, 0.87, false, 0.12, 0.02, 0.06, 0.10);
        TH1D* fA = Fr(-0.5, 4.5, 0.75, 3.0, "", "R = data (N - CB) / sim, M_{4} 0.60 - 0.80", 0.040, 0.044, 1.25);
        fA->GetXaxis()->SetLabelSize(0); fA->GetXaxis()->SetTickLength(0);
        TLatex lb; lb.SetTextFont(42); lb.SetTextSize(0.031); lb.SetTextAlign(23);
        const char* vzShort[5] = {"v_{z} > -160", "-140..-100", "-160..-140", "v_{z} > -100", "no v_{z} cut"};
        for (int iv = 0; iv < 5; ++iv) lb.DrawLatex(iv, 0.72, vzShort[iv]);
        const int vcol[3] = {kBlue + 1, kGreen + 2, kOrange + 7}, vst[3] = {20, 21, 22};
        TLegend* lg = new TLegend(0.14, 0.77, 0.97, 0.93); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.030); lg->SetNColumns(3);
        for (int v = 0; v < 3; ++v)
            for (int io = 0; io < 3; ++io) {
                TGraphErrors* g = new TGraphErrors();
                for (int iv = 0; iv < 5; ++iv) {
                    const RV q = DataSim(DW(0, iv, kIR), ST(0, iv, io), v, kIR);
                    g->SetPoint(iv, iv - 0.27 + 0.18 * v + 0.05 * io, q.r); g->SetPointError(iv, 0, q.e);
                }
                const int stl = io == 0 ? vst[v] : io == 1 ? vst[v] + 4 : (v == 0 ? 5 : v == 1 ? 25 : 26);
                g->SetMarkerStyle(io == 2 ? (v == 0 ? 34 : v == 1 ? 28 : 32) : stl); g->SetMarkerColor(vcol[v]); g->SetLineColor(vcol[v]);
                g->SetMarkerSize(1.8); g->SetLineWidth(2);
                g->Draw("P SAME");
                lg->AddEntry(g, Form("%s, %s", v == 0 ? "T1+T2" : v == 1 ? "(a)" : "(b)", kOptName[io]), "p");
            }
        lg->Draw();
        dcmp::HLine(-0.5, 4.5, 1.0, kGray + 2, 2, 1);
        // right: R per window for pre and central, mode A, w
        c->cd();
        Pad(0.5, 0.0, 1.0, 0.87, false, 0.12, 0.02, 0.06, 0.10);
        TH1D* fB = Fr(-0.5, 3.5, 0.75, 3.0, "", "R = data (N - CB) / sim (mode A, w)", 0.040, 0.044, 1.25);
        fB->GetXaxis()->SetLabelSize(0); fB->GetXaxis()->SetTickLength(0);
        const int kw[4] = {9, 8, 10, 11};
        for (int i = 0; i < 4; ++i) lb.DrawLatex(i, 0.72, Form("M_{4} %.2f - %.2f", kRegs[kw[i]].a, kRegs[kw[i]].b));
        TLegend* lg2 = new TLegend(0.14, 0.77, 0.97, 0.93); lg2->SetBorderSize(0); lg2->SetFillStyle(0); lg2->SetTextSize(0.030); lg2->SetNColumns(2);
        for (int iv : {0, 1, 2})
            for (int v = 0; v < 3; ++v) {
                TGraphErrors* g = new TGraphErrors();
                for (int i = 0; i < 4; ++i) {
                    const RV q = DataSim(DW(0, iv, kw[i]), ST(0, iv, 0), v, kw[i]);
                    g->SetPoint(i, i - 0.3 + 0.2 * iv + 0.06 * v, q.r); g->SetPointError(i, 0, q.e);
                }
                g->SetMarkerStyle(iv == 0 ? vst[v] : iv == 1 ? vst[v] + 4 : 33); g->SetMarkerColor(vcol[v]); g->SetLineColor(vcol[v]);
                g->SetMarkerSize(1.8); g->SetLineWidth(2); g->Draw("P SAME");
                lg2->AddEntry(g, Form("%s, %s", kVzSl[iv].name, v == 0 ? "T1+T2" : v == 1 ? "(a)" : "(b)"), "p");
            }
        lg2->Draw();
        dcmp::HLine(-0.5, 3.5, 1.0, kGray + 2, 2, 1);
        Save(c, "review_R_variants");
        delete c;
    }

    // ---------------- figure: data/sim against M4 for the variants ----------------
    {
        Binning b; b.nbin = 80; b.lo = 0.30; b.hi = 1.10;
        struct Cfg { int iv, v, col, st; const char* lab; };
        const std::vector<Cfg> cfg = {{0, 0, kBlue + 1, 20, "pre, T1 + T2 (default)"},
                                      {0, 1, kGreen + 2, 21, "pre, (a) + ghost-flagged true pairs"},
                                      {0, 2, kOrange + 7, 22, "pre, (b) no ghost rows"},
                                      {1, 0, kBlue + 1, 24, "central v_{z}, T1 + T2"},
                                      {1, 1, kGreen + 2, 25, "central v_{z}, (a)"},
                                      {2, 1, kMagenta + 1, 33, "edge v_{z} (-160..-140), (a)"}};
        std::vector<TH1D*> rAbs, rShp;
        for (const Cfg& q : cfg) {
            const WSpec w = Fill(sd, SSel(q.iv), CompSimVar(q.v), SimOpt(), b);
            const dcmp::DH h = dcmp::FillData(dd, DSel(q.iv), 80, 0.30, 1.10, 1);
            TH1D* hd = DHist(h, 4); TH1D* hs = w.Hist(Nm("hs"), 2);
            TH1D* r = dcmp::Ratio(hd, hs, Nm("rA"));
            const RV R = DataSim(DW(0, q.iv, kIR), ST(0, q.iv, 0), q.v, kIR);
            TH1D* r2 = (TH1D*)r->Clone(Nm("rS")); r2->Scale(1.0 / R.r);
            for (TH1D* x : {r, r2}) { x->SetMarkerStyle(q.st); x->SetMarkerColor(q.col); x->SetLineColor(q.col); x->SetMarkerSize(1.3); x->SetLineWidth(2); }
            rAbs.push_back(r); rShp.push_back(r2);
        }
        TCanvas* c = Canvas("cMv", 2000, 1700);
        Hdr(c, {{"Data (N - CB) / simulation against M_{4} for the ghost treatments and vertex-z slices (10 MeV bins)", 1},
                {kNoteL + "; " + kNoteC, kBlue + 1},
                {kNoteS, 1},
                {kNoteD + "; dotted: bump and peak windows of the data", kRed + 1}},
            0.985, 0.022, 0.0145);
        for (int pnl = 0; pnl < 2; ++pnl) {
            c->cd();
            Pad(0.0, pnl ? 0.0 : 0.455, 1.0, pnl ? 0.455 : 0.91, false, 0.08, 0.02, 0.05, 0.12);
            Fr(0.30, 1.10, pnl ? 0.3 : 0.4, pnl ? 1.7 : 2.9, kXM4, pnl ? "(data/sim) / R_{0.60-0.80}" : "data (N - CB) / simulation", 0.045, 0.050, 0.75);
            dcmp::HLine(0.30, 1.10, 1.0, kGray + 2, 2, 1);
            for (double x : {scan::kBumpLo, scan::kBumpHi, scan::kPeakLo, scan::kPeakHi}) dcmp::VLine(x, pnl ? 0.3 : 0.4, pnl ? 1.7 : 2.9, kGray + 1, 3);
            std::vector<TH1D*>& v = pnl ? rShp : rAbs;
            for (TH1D* h : v) h->Draw("E1 P SAME");
            if (pnl == 0) {
                TLegend* lg = new TLegend(0.09, 0.72, 0.97, 0.94); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.036); lg->SetNColumns(3);
                lg->SetHeader("absolute scale (R = data/sim in 0.60 - 0.80):");
                for (size_t i = 0; i < cfg.size(); ++i) {
                    const RV R = DataSim(DW(0, cfg[i].iv, kIR), ST(0, cfg[i].iv, 0), cfg[i].v, kIR);
                    lg->AddEntry(v[i], Form("%s (R = %.3f)", cfg[i].lab, R.r), "p");
                }
                lg->Draw();

            } else {
                Txt(0.97, 0.88, "shape: each curve divided by its own R in 0.60 - 0.80", 0.045, 1, 31);
            }
        }
        Save(c, "review_m4_variants");
        delete c;
    }

    // ---------------- figure: M4 spectra in the central and the edge vz slice ----------------
    {
        TCanvas* c = Canvas("cMc", 2200, 2000);
        Hdr(c, {{"M_{4} spectra of data and simulation in the central (-140 < v_{z} < -100 mm) and the edge (-160 < v_{z} < -140 mm) target slices", 1},
                {kNoteL + "; " + kNoteC, kBlue + 1},
                {kNoteS + "; fake lepton and different sources not shown", 1},
                {"solid teal: T1 + T2 (default); dotted teal: (a) with ghost-flagged true pairs (dotted azure: those pairs alone); ratio pads: blue N - CB / (T1 + T2), open green N - CB / (a); zoom: (a) stacked", 1},
                {kNoteD, kRed + 1}},
            0.985, 0.020, 0.0135);
        Binning b; b.nbin = 160; b.lo = 0.30; b.hi = 1.10;
        const std::vector<Comp> cc = {CompSimVar(0), CompSimVar(1), CompPhysBkg(), CompGhostPair(), CompBumpDalitz(), CompBumpConv(), CompPeakConv()};
        int col = 0;
        for (int iv : {1, 2}) {
            const std::vector<WSpec> w = FillMany(sd, SSel(iv), cc, SimOpt(), b);
            const dcmp::DH h = dcmp::FillData(dd, DSel(iv), 160, 0.30, 1.10, 2);
            TH1D* hN = DHist(h, 0); TH1D* hC = DHist(h, 3); TH1D* hD = DHist(h, 4);
            dcmp::StyleM(hN, 1, 20, 0.9); dcmp::StyleM(hC, 2, 21, 0.8); dcmp::StyleM(hD, 4, 20, 0.9);
            std::vector<TH1D*> hs(cc.size());
            for (size_t k = 0; k < cc.size(); ++k) hs[k] = w[k].Hist(Nm("hm"), 2);
            dcmp::StyleL(hs[0], kColSim, 3, 1); dcmp::StyleL(hs[1], kColSim, 3, 2); dcmp::StyleL(hs[2], kColPhysBkg, 1, 1, kColPhysFill);
            dcmp::StyleL(hs[3], kColFake, 2, 3); dcmp::StyleL(hs[4], kColBumpDalitz, 2, 1); dcmp::StyleL(hs[5], kColBumpConv, 2, 1);
            dcmp::StyleL(hs[6], kColPeakConv, 2, 1);
            const double x0 = 0.5 * col, x1 = x0 + 0.5;
            // log pad
            c->cd();
            Pad(x0, 0.56, x1, 0.885, true, 0.16, 0.02, 0.10, 0.0);
            TH1D* f = Fr(0.30, 1.10, 0.03, 3e7, "", "d#sigma/dM_{4} [nb/(GeV/c^{2})]", 0.045, 0.050, 1.25);
            f->GetXaxis()->SetLabelSize(0);
            hs[2]->Draw("HIST SAME");
            for (int k : {3, 4, 5, 6, 1, 0}) hs[k]->Draw("HIST SAME");
            hN->Draw("P SAME"); hC->Draw("P SAME"); hD->Draw("P SAME");
            gPad->RedrawAxis();
            Txt(0.17, 0.995, Form("%s mm (absolute scale)", kVzSl[iv].name), 0.050, 1, 13);
            TLegend* lg = new TLegend(0.18, 0.62, 0.48, 0.88); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.040);
            lg->AddEntry(hN, "data: N", "p"); lg->AddEntry(hC, "data: CB", "p"); lg->AddEntry(hD, "data: N - CB", "p"); lg->Draw();
            TLegend* lg2 = new TLegend(0.49, 0.56, 0.98, 0.89); lg2->SetBorderSize(0); lg2->SetFillStyle(0); lg2->SetTextSize(0.038);
            lg2->AddEntry(hs[0], "sim T1 + T2", "l"); lg2->AddEntry(hs[1], "sim (a) T1 + T2 + ghost pairs", "l");
            lg2->AddEntry(hs[2], "sim physical background", "f"); lg2->AddEntry(hs[3], "sim ghost-flagged true pairs", "l");
            lg2->AddEntry(hs[4], "sim bump, #pi^{0} Dalitz", "l"); lg2->AddEntry(hs[5], "sim bump, conversion", "l");
            lg2->AddEntry(hs[6], "sim peak, conversion", "l"); lg2->SetTextSize(0.034); lg2->Draw();
            // ratio pad
            c->cd();
            Pad(x0, 0.43, x1, 0.56, false, 0.16, 0.02, 0.0, 0.30);
            TH1D* fr = Fr(0.30, 1.10, 0.3, 2.7, kXM4, "data/sim", 0.10, 0.11, 0.55, 1.0);
            (void)fr;
            dcmp::HLine(0.30, 1.10, 1.0, kGray + 2, 2, 1);
            TH1D* r0 = dcmp::Ratio(hD, hs[0], Nm("r0")); dcmp::StyleM(r0, 4, 20, 0.8);
            TH1D* ra = dcmp::Ratio(hD, hs[1], Nm("ra")); dcmp::StyleM(ra, kGreen + 2, 24, 0.8);
            r0->Draw("E1 P SAME"); ra->Draw("E1 P SAME");
            // zoom
            c->cd();
            Pad(x0, 0.0, x1, 0.40, false, 0.16, 0.02, 0.08, 0.12);
            TH1D* bg = (TH1D*)hs[2]->Clone(Nm("zb")); bg->Add(hs[3]);
            TH1D* s1 = (TH1D*)bg->Clone(Nm("z1")); s1->Add(hs[4]);
            TH1D* s2 = (TH1D*)s1->Clone(Nm("z2")); s2->Add(hs[5]);
            TH1D* s3 = (TH1D*)s2->Clone(Nm("z3")); s3->Add(hs[6]);
            double ym = 0;
            for (int i = 1; i <= hD->GetNbinsX(); ++i) if (hD->GetBinCenter(i) > 0.40 && hD->GetBinCenter(i) < 0.70) ym = std::max({ym, hD->GetBinContent(i) + hD->GetBinError(i), s3->GetBinContent(i)});
            Fr(0.40, 0.70, 0, 1.5 * ym, kXM4, "d#sigma/dM_{4} [nb/(GeV/c^{2})]", 0.040, 0.045, 1.45);
            dcmp::StyleL(s3, kColPeakConv, 1, 1, kColPeakConv); dcmp::StyleL(s2, kColBumpConv, 1, 1, kColBumpConv);
            dcmp::StyleL(s1, kColBumpDalitz, 1, 1, kColBumpDalitz); dcmp::StyleL(bg, kColPhysBkg, 1, 1, kColPhysFill);
            s3->Draw("HIST SAME"); s2->Draw("HIST SAME"); s1->Draw("HIST SAME"); bg->Draw("HIST SAME");
            TH1D* hs0 = (TH1D*)hs[0]->Clone(Nm("z0")); hs0->Draw("HIST SAME");
            TH1D* hD2 = (TH1D*)hD->Clone(Nm("zD")); hD2->Draw("E1 P SAME");
            TH1D* hC2 = (TH1D*)hC->Clone(Nm("zC")); hC2->Draw("E1 P SAME");
            gPad->RedrawAxis();
            TLegend* lz = new TLegend(0.17, 0.55, 0.62, 0.93); lz->SetBorderSize(0); lz->SetFillStyle(0); lz->SetTextSize(0.036);
            lz->AddEntry(hD2, "data: N - CB", "p"); lz->AddEntry(hC2, "data: CB", "p");
            lz->AddEntry(s3, "sim (a): + peak, conversion", "f"); lz->AddEntry(s2, "sim (a): + bump, conversion", "f");
            lz->AddEntry(s1, "sim (a): + bump, #pi^{0} Dalitz", "f"); lz->AddEntry(bg, "sim (a): phys. bkg + ghost pairs", "f");
            lz->AddEntry(hs0, "sim T1 + T2 (default)", "l"); lz->Draw();
            Txt(0.17, 0.995, Form("%s mm: linear zoom, absolute scale, 5 MeV bins", kVzSl[iv].name), 0.040, 1, 13);
            ++col;
        }
        Save(c, "review_m4_vzslices");
        delete c;
    }
    printf("Part1 done %.0f s\n", sw.RealTime());
}

// ============================================================================
// Part 2: vertex z of data and simulation
// ============================================================================
void Part2(const SimData& sd, const scan::Data& dd) {
    TStopwatch sw;
    const int nb = 46;
    const double lo = -260, hi = -30;
    const std::vector<Win> regs = {{"0.30-0.455 (below the bump)", 0.30, 0.455}, {"bump window", scan::kBumpLo, scan::kBumpHi},
                                   {"peak window", scan::kPeakLo, scan::kPeakHi}, {"0.56-0.70", 0.56, 0.70},
                                   {"0.60-0.80 (R window)", 0.60, 0.80},          {"0.90-1.10", 0.90, 1.10}};
    Binning b; b.nbin = nb; b.lo = lo; b.hi = hi; b.x = [](const Row& r) { return (double)r.vz; };
    FILE* f = dcmp::OpenTab("review_vz.txt");
    P2(f, "# Review finding 3: vertex z of data and simulation, selection pre WITHOUT the vz > -160 cut (isBest >= 0, PT3, cor_ok, z > -500)\n");
    P2(f, "# data N - CB per 5 mm vz bin (CB = 2 sqrt(N++ N--) of the same vz bin), L_data = %.1f nb^-1 UNCONFIRMED; sim std + lep, mode A, w, L_sim = %.2f\n", kLdata, kLsim);
    P2(f, "# sim T12 = T1 + T2 (default); (a) = T1 + T2 + ghost-flagged true pairs. sigma per bin [nb]. eta -> pi+pi-e+e- not in the simulation.\n");
    TCanvas* c = Canvas("cVz", 2400, 2000);
    Hdr(c, {{"Vertex z of data and simulation in six M_{4} regions; selection pre without the v_{z} > -160 mm cut (dashed lines: -160, -140, -100 mm)", 1},
            {kNoteL + "; data: N black, CB red (per 5 mm v_{z} bin), N - CB blue", kBlue + 1},
            {kNoteS + "; solid teal T1 + T2, dotted green (a) with ghost-flagged true pairs", 1},
            {kNoteD, kRed + 1}},
        0.985, 0.020, 0.0135);
    for (size_t ir = 0; ir < regs.size(); ++ir) {
        const double a = regs[ir].a, bb = regs[ir].b;
        const Sel s = [a, bb](const Row& r) { return PassBaseNoVz(r) && r.m4 >= a && r.m4 < bb; };
        const scan::Sel dsel = [a, bb](const scan::Row& r) { return PassBaseNoVz(r) && r.m4 >= a && r.m4 < bb; };
        const std::vector<WSpec> w = FillMany(sd, s, {CompSimVar(0), CompSimVar(1)}, SimOpt(), b);
        const dcmp::DH h = dcmp::FillData(dd, dsel, nb, lo, hi, 1, dcmp::DVar([](const scan::Row& r) { return (double)r.vz; }));
        // unit: nb per mm
        TH1D* hN = h.Hist(0, Nm("vN"), 2); TH1D* hC = h.Hist(3, Nm("vC"), 2); TH1D* hD = h.Hist(4, Nm("vD"), 2);
        TH1D* s0 = w[0].Hist(Nm("v0"), 2); TH1D* s1 = w[1].Hist(Nm("v1"), 2);
        P2(f, "\n## M4 %s\n%-14s %8s %8s %18s %18s %18s %15s %15s\n", regs[ir].name, "vz [mm]", "N", "CB", "data N-CB [nb]", "sim T12 [nb]",
           "sim (a) [nb]", "data/T12", "data/(a)");
        for (int i = 0; i < nb; ++i) {
            const double ds = h.sig[i] / kLdata, de = std::sqrt(std::max(0., h.sigVar[i])) / kLdata;
            const double x0 = w[0].s[i] / kLsim, e0 = std::sqrt(w[0].v[i]) / kLsim, x1 = w[1].s[i] / kLsim, e1 = std::sqrt(w[1].v[i]) / kLsim;
            const RV q0 = Div(ds, de, x0, e0), q1 = Div(ds, de, x1, e1);
            P2(f, "%6.0f..%-6.0f %8.0f %8.1f %8.4g +- %-7.2g %8.4g +- %-7.2g %8.4g +- %-7.2g %6.3f +- %-6.3f %6.3f +- %-6.3f\n", lo + i * 5, lo + (i + 1) * 5,
               h.n[0][i], h.cb[i], ds, de, x0, e0, x1, e1, q0.r, q0.e, q1.r, q1.e);
        }
        // panel
        const int ix = ir % 3, iy = ir / 3;
        const double x1p = ix / 3.0, x2p = (ix + 1) / 3.0, ytop = iy == 0 ? 0.90 : 0.45, ybot = iy == 0 ? 0.45 : 0.0;
        dcmp::PP P = dcmp::MkPair(c, x1p, ybot, x2p, ytop, 0.33, true, 0.17, 0.08);
        P.up->cd();
        dcmp::StyleM(hN, 1, 20, 0.8); dcmp::StyleM(hC, 2, 21, 0.7); dcmp::StyleM(hD, 4, 20, 0.8);
        dcmp::StyleL(s0, kColSim, 3, 1); dcmp::StyleL(s1, kGreen + 2, 3, 2);
        const double ymax = dcmp::MaxOf({hN, s1}, false);
        TH1D* fu = dcmp::Frame(hN, ymax * 2e-4, ymax * 5, "v_{z} [mm]", "d#sigma/dv_{z} [nb/mm]", 0.017, P.hu, false, 1.45);
        fu->Draw("AXIS");
        s0->Draw("HIST SAME"); s1->Draw("HIST SAME");
        hN->Draw("P SAME"); hC->Draw("P SAME"); hD->Draw("P SAME");
        for (double x : {-160., -140., -100.}) dcmp::VLine(x, ymax * 2e-4, ymax * 5, kGray + 2, 2, 2);
        Txt(0.19, 0.84, Form("M_{4} %s", regs[ir].name), 0.050, 1, 11);
        if (ir == 0) {
            TLegend* lg = new TLegend(0.19, 0.50, 0.62, 0.80); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.042);
            lg->AddEntry(hN, "data N", "p"); lg->AddEntry(hC, "data CB", "p"); lg->AddEntry(hD, "data N - CB", "p");
            lg->AddEntry(s0, "sim T1 + T2", "l"); lg->AddEntry(s1, "sim (a)", "l"); lg->Draw();
        }
        P.lo->cd();
        TH1D* r0 = dcmp::Ratio(hD, s0, Nm("vr0")); dcmp::StyleM(r0, 4, 20, 0.8);
        TH1D* r1 = dcmp::Ratio(hD, s1, Nm("vr1")); dcmp::StyleM(r1, kGreen + 2, 24, 0.8);
        TH1D* fl = dcmp::Frame(hN, 0.0, 3.0, "v_{z} [mm]", "data/sim", 0.017, P.hl, true, 1.45);
        fl->GetYaxis()->SetTitleOffset(1.45 * P.hl / P.hu); fl->Draw("AXIS");
        dcmp::HLine(lo, hi, 1.0, kGray + 2, 2, 1);
        for (double x : {-160., -140., -100.}) dcmp::VLine(x, 0, 3, kGray + 2, 2, 2);
        r0->Draw("E1 P SAME"); r1->Draw("E1 P SAME");
    }
    if (f) fclose(f);
    Save(c, "review_vz");
    delete c;
    printf("Part2 done %.0f s\n", sw.RealTime());
}

// ============================================================================
// Part 4: OA(e+e-) and M(e+e-) classes with the variants
// ============================================================================
void Part4(const SimData& sd, const scan::Data& dd) {
    TStopwatch sw;
    const std::vector<Win> regs = {{"M4 0.56-0.80", 0.56, 0.80}, {"bump window", scan::kBumpLo, scan::kBumpHi}, {"M4 0.30-0.45", 0.30, 0.45}};
    const std::vector<std::pair<double, double>> wl = WinList(regs);
    // data: per vz (pre, central): all, 5 OA classes, 5 M(ee) classes, 7 cumulative OA > x
    const double cumX[7] = {0, 1, 2, 3, 4, 6, 10};
    std::vector<scan::Sel> ds;
    for (int iv : {0, 1}) {
        const scan::Sel b = DSel(iv);
        ds.push_back(b);
        for (int k = 0; k < 5; ++k) ds.push_back([b, k](const scan::Row& r) { return b(r) && dcmp::OAcls(r.oa) == k; });
        for (int k = 0; k < 5; ++k) ds.push_back([b, k](const scan::Row& r) { return b(r) && dcmp::MEEcls(r.mee) == k; });
        for (int k = 0; k < 7; ++k) { const float x = cumX[k]; ds.push_back([b, x](const scan::Row& r) { return b(r) && r.oa > x; }); }
    }
    scan::FillOpt fo; fo.nbin = 320; fo.lo = 0.30; fo.hi = 1.10; fo.reb = 4; fo.windows = wl;
    fo.windows.push_back(std::make_pair(0.60, 0.80));
    const std::vector<scan::Spectra> dsp = scan::FillMany(dd, ds, fo);
    auto DW = [&](int ivi, int s, int r) -> const scan::WinSum& { return dsp[ivi * 18 + s].win[2 + r]; };
    // simulation: configurations (vz index into {pre, central}, mode)
    const std::vector<Comp> comps = {CompSimVar(0), CompSimVar(1), CompSimVar(2)};
    Binning bOA; bOA.nbin = 5; bOA.lo = 0; bOA.hi = 5; bOA.x = [](const Row& r) { return dcmp::OAcls(r.oa) + 0.5; };
    Binning bME; bME.nbin = 5; bME.lo = 0; bME.hi = 5; bME.x = [](const Row& r) { return dcmp::MEEcls(r.mee) + 0.5; };
    // [ivi][mode][kind][region] -> WSpec per comp
    std::vector<WSpec> S[2][2][2][3];
    WT tot[2][2];
    for (int ivi = 0; ivi < 2; ++ivi)
        for (int md = 0; md < 2; ++md) {
            SimOpt o; o.mode = md;
            for (int kind = 0; kind < 2; ++kind)
                for (int r = 0; r < 3; ++r) {
                    const double a = regs[r].a, b = regs[r].b;
                    const Sel s = SelAnd(SSel(ivi), [a, b](const Row& x) { return x.m4 >= a && x.m4 < b; });
                    S[ivi][md][kind][r] = FillMany(sd, s, comps, o, kind ? bME : bOA);
                }
            tot[ivi][md] = SimWin(sd, SSel(ivi), comps, wl, o);
        }
    // cumulative: R(OA > x) in 0.60-0.80, mode A and B
    WT cum[2][2][7];
    for (int ivi = 0; ivi < 2; ++ivi)
        for (int md = 0; md < 2; ++md)
            for (int k = 0; k < 7; ++k) {
                SimOpt o; o.mode = md;
                const float x = cumX[k];
                cum[ivi][md][k] = SimWin(sd, SelAnd(SSel(ivi), [x](const Row& r) { return r.oa > x; }), comps, {{0.60, 0.80}}, o);
            }
    printf("Part4 sums %.0f s\n", sw.RealTime()); sw.Continue();

    FILE* f = dcmp::OpenTab("review_oa_classes.txt");
    P2(f, "# Review finding 5: data (N - CB) / simulation per OA(e+e-) and M(e+e-) class; CB per class in 10 MeV M4 bins (no k); L_data = %.1f UNCONFIRMED\n", kLdata);
    P2(f, "# variants: vz slice pre (vz > -160) or central (-140..-100); mode A or B; ghost treatment T12 (default), (a) + ghost-flagged true pairs, (b) no ghost rows\n");
    P2(f, "# 'rel' = ratio of the class / ratio of all classes of the same region and variant (shape across the classes)\n");
    const char* cfgName[4] = {"pre A", "pre B", "central A", "central B"};
    double relv[2][3][2][2][3][5];   // [kind][region][ivi][mode][var][class]
    double rele[2][3][2][2][3][5];
    for (int kind = 0; kind < 2; ++kind)
        for (int r = 0; r < 3; ++r) {
            P2(f, "\n## %s classes, %s\n%-8s %-10s %16s", kind ? "M(e+e-) [MeV/c2]" : "OA(e+e-) [deg]", regs[r].name, "class", "config", "data N-CB [nb]");
            for (int v = 0; v < 3; ++v) P2(f, " %15s %8s", Form("abs %s", v == 0 ? "T12" : v == 1 ? "(a)" : "(b)"), "rel");
            P2(f, "\n");
            for (int k = 0; k <= 5; ++k)
                for (int ivi = 0; ivi < 2; ++ivi)
                    for (int md = 0; md < 2; ++md) {
                        const scan::WinSum& w = k < 5 ? DW(ivi, 1 + 5 * kind + k, r) : DW(ivi, 0, r);
                        const scan::WinSum& wa = DW(ivi, 0, r);
                        P2(f, "%-8s %-10s %8.4g +- %-5.2g", k < 5 ? (kind ? dcmp::kMEEname[k] : dcmp::kOAname[k]) : "all", cfgName[ivi * 2 + md],
                           w.sig() / kLdata, w.sigErr() / kLdata);
                        for (int v = 0; v < 3; ++v) {
                            const WSpec& sp = S[ivi][md][kind][r][v];
                            const double xs = k < 5 ? sp.s[k] / kLsim : sp.S / kLsim, xe = k < 5 ? std::sqrt(sp.v[k]) / kLsim : std::sqrt(sp.V) / kLsim;
                            const RV q = Div(w.sig() / kLdata, w.sigErr() / kLdata, xs, xe);
                            const RV qa = Div(wa.sig() / kLdata, wa.sigErr() / kLdata, sp.S / kLsim, std::sqrt(sp.V) / kLsim);
                            const double rel = qa.r > 0 ? q.r / qa.r : 0;
                            if (k < 5) { relv[kind][r][ivi][md][v][k] = rel; rele[kind][r][ivi][md][v][k] = qa.r > 0 ? q.e / qa.r : 0; }
                            P2(f, " %6.3f +- %-6.3f %8.3f", q.r, q.e, rel);
                        }
                        P2(f, "\n");
                    }
        }
    P2(f, "\n## cumulative OA > x: R = data / sim in M4 0.60-0.80 (as the cuts study frames it); rel = R(OA > x) / R(OA > 0)\n");
    P2(f, "%-6s %-10s %16s", "OA > x", "config", "data N-CB [nb]");
    for (int v = 0; v < 3; ++v) P2(f, " %15s %7s", Form("R %s", v == 0 ? "T12" : v == 1 ? "(a)" : "(b)"), "rel");
    P2(f, "\n");
    double cumv[2][2][3][7], cume[2][2][3][7];
    for (int k = 0; k < 7; ++k)
        for (int ivi = 0; ivi < 2; ++ivi)
            for (int md = 0; md < 2; ++md) {
                const scan::WinSum& w = DW(ivi, 11 + k, 3);
                P2(f, "%-6.0f %-10s %8.4g +- %-5.2g", cumX[k], cfgName[ivi * 2 + md], w.sig() / kLdata, w.sigErr() / kLdata);
                for (int v = 0; v < 3; ++v) {
                    const RV q = DataSim(w, cum[ivi][md][k], v, 0);
                    const RV q0 = DataSim(DW(ivi, 11, 3), cum[ivi][md][0], v, 0);
                    cumv[ivi][md][v][k] = q.r / q0.r; cume[ivi][md][v][k] = q.e / q0.r;
                    P2(f, " %6.3f +- %-6.3f %7.3f", q.r, q.e, q.r / q0.r);
                }
                P2(f, "\n");
            }
    if (f) fclose(f);

    // figure
    TCanvas* c = Canvas("cOA", 2400, 1900);
    Hdr(c, {{"Data (N - CB) / simulation per e^{+}e^{-} class, relative to all classes of the same variant (shape across the classes)", 1},
            {kNoteL + "; CB of each class in 10 MeV M_{4} bins (no k); errors: data and simulation statistics of the class", kBlue + 1},
            {kNoteS, 1},
            {"variants: blue T1 + T2 (default), green (a) + ghost-flagged true pairs, orange (b) no ghost rows; filled pre (v_{z} > -160), open central v_{z}; circles mode A, squares mode B", 1},
            {kNoteD, kRed + 1}},
        0.985, 0.020, 0.0135);
    const int vcol[3] = {kBlue + 1, kGreen + 2, kOrange + 7};
    struct PanelDef { int kind, r; const char* title; };
    const PanelDef pd[3] = {{0, 0, "OA(e^{+}e^{-}) classes, M_{4} 0.56 - 0.80"}, {0, 1, "OA(e^{+}e^{-}) classes, bump window"},
                            {1, 0, "M(e^{+}e^{-}) classes, M_{4} 0.56 - 0.80"}};
    for (int p = 0; p < 4; ++p) {
        c->cd();
        const double x1 = (p % 2) * 0.5, y2 = p < 2 ? 0.88 : 0.44;
        Pad(x1, y2 - 0.44, x1 + 0.5, y2, false, 0.12, 0.02, 0.08, 0.12);
        if (p < 3) {
            TH1D* fr = Fr(0, 5, 0.7, 1.35, "", "(data/sim of the class) / (data/sim of all)", 0.042, 0.042, 1.3);
            for (int k = 0; k < 5; ++k) fr->GetXaxis()->SetBinLabel(k * 20 + 10, pd[p].kind ? dcmp::kMEEname[k] : dcmp::kOAname[k]);
            fr->GetXaxis()->SetLabelSize(0.055); fr->GetXaxis()->LabelsOption("h");
            dcmp::HLine(0, 5, 1.0, kGray + 2, 2, 1);
            for (int ivi = 0; ivi < 2; ++ivi)
                for (int md = 0; md < 2; ++md)
                    for (int v = 0; v < 3; ++v) {
                        if (md == 1 && v == 2) continue;
                        TGraphErrors* g = new TGraphErrors();
                        for (int k = 0; k < 5; ++k) {
                            g->SetPoint(k, k + 0.5 - 0.33 + 0.22 * v + 0.11 * ivi + 0.05 * md, relv[pd[p].kind][pd[p].r][ivi][md][v][k]);
                            g->SetPointError(k, 0, rele[pd[p].kind][pd[p].r][ivi][md][v][k]);
                        }
                        g->SetMarkerStyle(md == 0 ? (ivi == 0 ? 20 : 24) : (ivi == 0 ? 21 : 25)); g->SetMarkerColor(vcol[v]); g->SetLineColor(vcol[v]);
                        g->SetMarkerSize(1.5); g->SetLineWidth(2); g->Draw("P SAME");
                    }
            Txt(0.14, 0.86, pd[p].title, 0.050, 1, 11);
            if (pd[p].kind == 1) Txt(0.14, 0.80, "M(e^{+}e^{-}) classes in MeV/c^{2}", 0.040, 1, 11);
        } else {
            Fr(-0.5, 10.5, 0.7, 1.15, "OA(e^{+}e^{-}) > x [deg]", "R(OA > x) / R(OA > 0), M_{4} 0.60 - 0.80", 0.042, 0.042, 1.3);
            dcmp::HLine(-0.5, 10.5, 1.0, kGray + 2, 2, 1);
            TLegend* lg = new TLegend(0.14, 0.14, 0.75, 0.42); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.034); lg->SetNColumns(2);
            for (int ivi = 0; ivi < 2; ++ivi)
                for (int md = 0; md < 2; ++md)
                    for (int v = 0; v < 3; ++v) {
                        if (md == 1 && v == 2) continue;
                        TGraphErrors* g = new TGraphErrors();
                        for (int k = 0; k < 7; ++k) { g->SetPoint(k, cumX[k] + 0.08 * v + 0.04 * ivi, cumv[ivi][md][v][k]); g->SetPointError(k, 0, cume[ivi][md][v][k]); }
                        g->SetMarkerStyle(md == 0 ? (ivi == 0 ? 20 : 24) : (ivi == 0 ? 21 : 25)); g->SetMarkerColor(vcol[v]); g->SetLineColor(vcol[v]);
                        g->SetMarkerSize(1.5); g->SetLineWidth(2); g->Draw("PL SAME");
                        if (md == 0) lg->AddEntry(g, Form("%s, %s", ivi ? "central" : "pre", v == 0 ? "T1+T2" : v == 1 ? "(a)" : "(b)"), "p");
                    }
            lg->Draw();
            Txt(0.14, 0.86, "cumulative cut OA > x (framing of the cuts study)", 0.050, 1, 11);
        }
    }
    Save(c, "review_oa_classes");
    delete c;
    printf("Part4 done %.0f s\n", sw.RealTime());
}

// ============================================================================
// Part 8: position and width of the simulated conversion peak
// ============================================================================
void Part8(const SimData& sd) {
    FILE* f = dcmp::OpenTab("review_peak.txt");
    P2(f, "# Review finding 8: Gaussian fits of the simulated peak-conversion class (eta -> pi+pi-gamma, photon converted), pre, mode A, w (PDG)\n");
    P2(f, "# weighted histogram (event-group errors), chi2 fit with the bin errors (empty bins skipped), range mean +- n sigma iterated 6 times\n");
    P2(f, "# data peak (all decays incl. eta -> pi+pi-e+e-, data study, other selection: no vz cut): 0.5439 GeV/c2, sigma 5.4 MeV\n");
    P2(f, "%-10s %-8s %-6s %22s %20s %12s\n", "rows", "bin", "range", "mean [GeV/c2]", "sigma [MeV]", "chi2/ndf");
    struct Res { int g, ib, ir; double m, em, s, es; };
    std::vector<Res> res;
    const double bws[3] = {0.0005, 0.001, 0.002};
    const double nsg[4] = {1.5, 2.0, 2.5, 3.0};
    TH1D* hShow[2] = {nullptr, nullptr};
    TF1* fShow[2][2] = {{nullptr, nullptr}, {nullptr, nullptr}};
    for (int g = 0; g < 2; ++g)
        for (int ib = 0; ib < 3; ++ib) {
            Binning b; b.lo = 0.49; b.hi = 0.60; b.nbin = (int)std::lround((b.hi - b.lo) / bws[ib]);
            const WSpec w = Fill(sd, SelPreset(kPre), CompSigCls(4, g == 1), SimOpt(), b);
            TH1D* h = w.Hist(Nm("pk"), 0);
            if (ib == 1) hShow[g] = h;
            for (int ir = 0; ir < 4; ++ir) {
                TF1* fg = new TF1(Nm("fg"), "gaus", 0.535, 0.556);
                fg->SetParameters(h->GetMaximum(), 0.5458, 0.0065);
                h->Fit(fg, "QNR0");
                for (int it = 0; it < 6; ++it) {
                    const double m = fg->GetParameter(1), s = std::fabs(fg->GetParameter(2));
                    fg->SetRange(m - nsg[ir] * s, m + nsg[ir] * s);
                    h->Fit(fg, "QNR0");
                }
                Res r{g, ib, ir, fg->GetParameter(1), fg->GetParError(1), std::fabs(fg->GetParameter(2)) * 1e3, fg->GetParError(2) * 1e3};
                res.push_back(r);
                P2(f, "%-10s %-8s %-6s %10.5f +- %-8.5f %8.2f +- %-7.2f %6.1f/%-4d\n", g ? "no ghost" : "all", Form("%.1f MeV", bws[ib] * 1e3),
                   Form("+-%.1f", nsg[ir]), r.m, r.em, r.s, r.es, fg->GetChisquare(), fg->GetNDF());
                if (ib == 1 && (ir == 0 || ir == 3)) fShow[g][ir == 0 ? 0 : 1] = fg;
            }
        }
    // summary
    for (int gg = 0; gg < 4; ++gg) {
        const int g = gg % 2;
        const bool coarse = gg >= 2;   // 1 and 2 MeV bins only (0.5 MeV bins: many empty bins, skipped by the chi2 fit)
        std::vector<double> m, s, em, es;
        for (const Res& r : res) if (r.g == g && (!coarse || r.ib > 0)) { m.push_back(r.m); s.push_back(r.s); em.push_back(r.em); es.push_back(r.es); }
        if (coarse) P2(f, "\n(1 and 2 MeV bins only)");
        auto med = [](std::vector<double> v) { std::sort(v.begin(), v.end()); const size_t n = v.size(); return n % 2 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]); };
        auto rms = [](const std::vector<double>& v) { double a = 0, q = 0; for (double x : v) { a += x; q += x * x; } a /= v.size(); return std::sqrt(std::max(0., q / v.size() - a * a)); };
        P2(f, "\n%s rows: position median %.5f, range %.5f-%.5f, rms over fit choices %.5f, median stat error %.5f GeV/c2\n", g ? "no-ghost" : "all",
           med(m), *std::min_element(m.begin(), m.end()), *std::max_element(m.begin(), m.end()), rms(m), med(em));
        P2(f, "%s rows: width median %.2f, range %.2f-%.2f, rms over fit choices %.2f, median stat error %.2f MeV\n", g ? "no-ghost" : "all", med(s),
           *std::min_element(s.begin(), s.end()), *std::max_element(s.begin(), s.end()), rms(s), med(es));
    }
    P2(f, "# quote: position = median (stat = median error) (fit = rms over binning and range); width likewise. Width minus data 5.4 MeV is not significant\n");
    P2(f, "# when stat and fit-choice errors are combined. M(pi+pi-) scale from K0S (spectra_eta.txt): sim 1.0018 +- 0.0005, data 0.9970 +- 0.0004 -> sim/data 1.0048 +- 0.0007\n");
    if (f) fclose(f);
    // figure
    TCanvas* c = Canvas("cPk", 2200, 1100);
    Hdr(c, {{"Simulated #eta #rightarrow #pi^{+}#pi^{-}#gamma (converted photon) peak: Gaussian position and width against the fit choices; selection pre, mode A, w", 1},
            {kNoteS, 1},
            {"the data peak (0.5439 GeV/c^{2}, 5.4 MeV) contains all eta decays incl. the direct #eta #rightarrow #pi^{+}#pi^{-}e^{+}e^{-} and comes from another selection (data study, no v_{z} cut)", kRed + 1}},
        0.985, 0.028, 0.017);
    c->cd();
    Pad(0.0, 0.0, 0.42, 0.90, false, 0.17, 0.02, 0.05, 0.11);
    double ym = 0;
    for (int i = 1; i <= hShow[0]->GetNbinsX(); ++i) ym = std::max(ym, hShow[0]->GetBinContent(i) + hShow[0]->GetBinError(i));
    for (int g = 0; g < 2; ++g) { hShow[g]->Scale(1.0 / (kLsim * 0.001)); }
    for (int g = 0; g < 2; ++g) for (int k = 0; k < 2; ++k) if (fShow[g][k]) fShow[g][k]->SetParameter(0, fShow[g][k]->GetParameter(0) / (kLsim * 0.001));
    Fr(0.51, 0.58, 0, 1.9 * ym / (kLsim * 0.001), kXM4, "d#sigma/dM_{4} [nb/(GeV/c^{2})]", 0.035, 0.040, 1.9);
    dcmp::StyleM(hShow[0], kColPeakConv, 20, 1.0); dcmp::StyleM(hShow[1], kGray + 2, 24, 0.9);
    hShow[0]->Draw("E1 P SAME"); hShow[1]->Draw("E1 P SAME");
    for (int g = 0; g < 2; ++g)
        for (int k = 0; k < 2; ++k)
            if (fShow[g][k]) { fShow[g][k]->SetLineColor(g ? kGray + 2 : kColPeakConv); fShow[g][k]->SetLineStyle(k ? 2 : 1); fShow[g][k]->SetLineWidth(2); fShow[g][k]->Draw("SAME"); }
    dcmp::VLine(0.5439, 0, 1.9 * ym / (kLsim * 0.001), kBlue + 1, 2, 2);
    TLegend* lg = new TLegend(0.19, 0.70, 0.97, 0.94); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.032);
    lg->AddEntry(hShow[0], "sim peak conversion, all rows (1 MeV)", "p");
    lg->AddEntry(hShow[1], "sim peak conversion, no ghost rows", "p");
    lg->AddEntry(fShow[0][0], "Gaussian, #pm1.5#sigma range", "l");
    lg->AddEntry(fShow[0][1], "Gaussian, #pm3#sigma range", "l");
    lg->Draw();
    Txt(0.19, 0.66, "dashed blue: data peak position 0.5439", 0.030, kBlue + 1);
    for (int q = 0; q < 2; ++q) {
        c->cd();
        Pad(0.42, q ? 0.0 : 0.45, 1.0, q ? 0.45 : 0.90, false, 0.13, 0.02, 0.05, 0.14);
        TH1D* fr = Fr(-0.5, 23.5, q ? 4.0 : 0.5435, q ? 10.0 : 0.5480, "fit choice (bin width / range)", q ? "Gaussian #sigma [MeV]" : "Gaussian mean [GeV/c^{2}]",
                      0.045, 0.050, q ? 1.0 : 1.25);
        (void)fr;
        TGraphErrors* g[2] = {new TGraphErrors(), new TGraphErrors()};
        int n[2] = {0, 0};
        for (const Res& r : res) {
            const double x = r.ib * 4 + r.ir + r.g * 12;
            g[r.g]->SetPoint(n[r.g], x, q ? r.s : r.m); g[r.g]->SetPointError(n[r.g], 0, q ? r.es : r.em); ++n[r.g];
        }
        for (int k = 0; k < 2; ++k) { g[k]->SetMarkerStyle(k ? 24 : 20); g[k]->SetMarkerColor(k ? kGray + 2 : kColPeakConv); g[k]->SetLineColor(k ? kGray + 2 : kColPeakConv); g[k]->SetMarkerSize(1.4); g[k]->SetLineWidth(2); g[k]->Draw("P SAME"); }
        dcmp::HLine(-0.5, 23.5, q ? 5.4 : 0.5439, kBlue + 1, 2, 2);
        Txt(0.12, 0.86, q ? "width (dashed: data 5.4 MeV)" : "position (dashed: data 0.5439)", 0.055, 1, 11);
        if (q) Txt(0.12, 0.20, "x: 0-11 all rows, 12-23 no ghost rows; groups of 4: 0.5, 1, 2 MeV bins; within a group #pm1.5, 2, 2.5, 3#sigma", 0.040, 1, 11);
    }
    Save(c, "review_peak_fits");
    delete c;
}

// ============================================================================
// Part 16: template-fit variants
// ============================================================================
void Part16(SimData& sd, const scan::Data& dd) {
    TStopwatch sw;
    const std::vector<stf::Model> M = stf::Models();
    const int use[4] = {0, 3, 5, 6};
    const double gm = scan::kPeakMean, gs = 0.00674;
    struct Cfg { const char* name; int iv, v; };
    const std::vector<Cfg> cfg = {{"pre, default T1+T2", 0, 0},
                                  {"pre, (a) ghost-flagged true pairs in the background", 0, 1},
                                  {"pre, (b) signal without ghost rows", 0, 2},
                                  {"central vz, default", 1, 0},
                                  {"central vz, (a)", 1, 1},
                                  {"edge vz, (a)", 2, 1}};
    FILE* f = dcmp::OpenTab("review_fit_variants.txt");
    P2(f, "# Review finding 6: template-fit variants (sim_templatefit.C machinery), 0.43-0.75, 2.5 MeV bins, mode A, w; PRELIMINARY\n");
    P2(f, "# eta -> pi+pi-e+e- (direct) is in no template; L_data = %.1f nb^-1 UNCONFIRMED (alphas scale with 1/L_data). G at %.4f, sigma %.5f\n", kLdata, gm, gs);
    P2(f, "# rel. errors per 2.5 MeV point (median over the fit range): data sqrt(V_data)/|y|, background template sqrt(C_bkg)/T_bkg\n");
    P2(f, "%-52s %-4s %16s %16s %16s %18s %10s %8s %22s %18s\n", "variant", "mod", "a_bump (a_eta)", "a_peak", "a_bkg (a_rest)", "extra", "chi2/ndf",
       "tvShare", "S implied (par, tpl)", "E implied");
    struct Out { int c, m; double ab, eab, ap, eap, ak, eak, chi; int ndf; };
    std::vector<Out> outs;
    for (size_t ic = 0; ic < cfg.size(); ++ic) {
        TmplSave save;
        RelabelTmpl(sd, cfg[ic].v, save);
        stf::TplSet T;
        stf::FillTpl(sd, SSel(cfg[ic].iv), stf::OptOf(0), 1.0, T);
        RestoreTmpl(sd, save);
        stf::DSet D;
        stf::FillData(dd, DSel(cfg[ic].iv), D);
        const stf::FitIn in = stf::MakeIn(D, T, 1, 0.43, 0.75, gm, gs);
        // relative errors
        std::vector<double> rd, rt;
        for (int q = 0; q < in.npts; ++q) {
            const double tb = in.T[stf::cVM][q] + in.T[stf::cRest][q];
            const double cb = in.CT[stf::cVM][stf::cVM][q * in.npts + q] + in.CT[stf::cRest][stf::cRest][q * in.npts + q] +
                              2 * in.CT[stf::cVM][stf::cRest][q * in.npts + q];
            if (std::fabs(in.y[q]) > 0) rd.push_back(std::sqrt(in.Vd[q * in.npts + q]) / std::fabs(in.y[q]));
            if (tb > 0) rt.push_back(std::sqrt(cb) / tb);
        }
        std::sort(rd.begin(), rd.end()); std::sort(rt.begin(), rt.end());
        P2(f, "\n## %s: data S(P_main) %.0f +- %.0f, E %.0f +- %.0f; median rel. error per point: data %.4f, background template %.4f\n", cfg[ic].name,
           D.pf.S, D.pf.eS, D.E, D.eE, rd[rd.size() / 2], rt[rt.size() / 2]);
        for (int j : use) {
            stf::FitOut o = stf::DoFit(in, M[j]);
            stf::Implied(o, M[j], T, D, gm, gs);
            std::string extra;
            if (j == 3) extra = Form("a_VM %.3f+-%.3f", o.p[2], o.e[2]);
            if (j == 5) extra = Form("b %.3f+-%.3f", o.p[3], o.e[3]);
            if (j == 6) extra = Form("b %.3f+-%.3f A_G %.0f+-%.0f", o.p[2], o.e[2], o.p[3], o.e[3]);
            // a_bump, a_peak, a_bkg per model
            double ab = o.p[0], eab = o.e[0], ap = 0, eap = 0, ak = 0, eak = 0;
            if (j == 0 || j == 5) { ap = o.p[1]; eap = o.e[1]; ak = o.p[2]; eak = o.e[2]; }
            if (j == 3) { ap = o.p[1]; eap = o.e[1]; ak = o.p[3]; eak = o.e[3]; }
            if (j == 6) { ap = o.p[0]; eap = o.e[0]; ak = o.p[1]; eak = o.e[1]; }
            outs.push_back({(int)ic, j, ab, eab, ap, eap, ak, eak, o.chi2, o.ndf});
            P2(f, "%-52s %-4s %7.3f +- %-6.3f %7.3f +- %-6.3f %7.3f +- %-6.3f %18s %5.0f/%-4d %8.3f %7.0f +- %-5.0f +- %-5.0f %7.0f +- %-5.0f\n", cfg[ic].name,
               M[j].tag.c_str(), ab, eab, ap, eap, ak, eak, extra.c_str(), o.chi2, o.ndf, o.tvShare, o.S, o.eS, o.tS, o.E, o.eE);
        }
        printf("fit variant %s done (%.0f s)\n", cfg[ic].name, sw.RealTime()); sw.Continue();
    }
    if (f) fclose(f);
    // figure: a_bump, a_peak, a_bkg, chi2/ndf per variant and model
    TCanvas* c = Canvas("cFv", 2200, 1500);
    Hdr(c, {{"PRELIMINARY template fit of data N - CB (0.43 - 0.75, 2.5 MeV): fitted factors for the ghost treatments and vertex-z slices", 1},
            {kNoteL + "; #alpha = 1 means agreement with SMASH", kBlue + 1},
            {kNoteS, 1},
            {kNoteD + "; models: M1 three templates, M4 split background, M6 background tilt, M7 tilt + Gaussian", kRed + 1}},
        0.985, 0.024, 0.0155);
    const int mcol[4] = {kBlack, kAzure + 1, kRed + 1, kGreen + 2}, mst[4] = {20, 21, 22, 33};
    const char* yt[4] = {"#alpha_{bump} (#alpha_{#eta} in M7)", "#alpha_{peak} (= #alpha_{#eta} in M7)", "#alpha_{bkg} (#alpha_{rest} in M4)", "#chi^{2}/ndf"};
    const double ylo[4] = {-0.3, -1.0, 0.9, 0.5}, yhi[4] = {1.0, 4.5, 2.4, 5.5};
    for (int p = 0; p < 4; ++p) {
        c->cd();
        const double x1 = (p % 2) * 0.5, y2 = p < 2 ? 0.90 : 0.45;
        Pad(x1, y2 - 0.45, x1 + 0.5, y2, false, 0.12, 0.02, 0.05, 0.20);
        TH1D* fr = Fr(-0.5, cfg.size() - 0.5, ylo[p], yhi[p], "", yt[p], 0.040, 0.045, 1.2);
        for (size_t ic = 0; ic < cfg.size(); ++ic) {
            const int bin = fr->GetXaxis()->FindBin((double)ic);
            fr->GetXaxis()->SetBinLabel(bin, Form("%s", ic == 0 ? "pre" : ic == 1 ? "pre (a)" : ic == 2 ? "pre (b)" : ic == 3 ? "central" : ic == 4 ? "central (a)" : "edge (a)"));
        }
        fr->GetXaxis()->SetLabelSize(0.060); fr->GetXaxis()->LabelsOption("h");
        dcmp::HLine(-0.5, cfg.size() - 0.5, p == 3 ? 1.0 : (p == 2 ? 1.0 : 1.0), kGray + 2, 2, 1);
        if (p < 2) dcmp::HLine(-0.5, cfg.size() - 0.5, 0.0, kGray + 1, 3, 1);
        TLegend* lg = (p == 0) ? new TLegend(0.14, 0.72, 0.97, 0.92) : nullptr;
        if (lg) { lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.040); lg->SetNColumns(4); }
        for (int im = 0; im < 4; ++im) {
            TGraphErrors* g = new TGraphErrors();
            int n = 0;
            for (const Out& o : outs) {
                if (o.m != use[im]) continue;
                const double x = o.c - 0.24 + 0.16 * im;
                const double y = p == 0 ? o.ab : p == 1 ? o.ap : p == 2 ? o.ak : o.chi / o.ndf;
                const double e = p == 0 ? o.eab : p == 1 ? o.eap : p == 2 ? o.eak : 0;
                g->SetPoint(n, x, y); g->SetPointError(n, 0, e); ++n;
            }
            g->SetMarkerStyle(mst[im]); g->SetMarkerColor(mcol[im]); g->SetLineColor(mcol[im]); g->SetMarkerSize(1.6); g->SetLineWidth(2);
            g->Draw("P SAME");
            if (lg) lg->AddEntry(g, M[use[im]].tag.c_str(), "p");
        }
        if (lg) lg->Draw();
    }
    Save(c, "review_fit_variants");
    delete c;
    printf("Part16 done %.0f s\n", sw.RealTime());
}

}  // namespace srev

// ============================================================================
void sim_review(int parts = 31) {
    using namespace srev;
    TStopwatch sw;
    TH1::AddDirectory(kFALSE);
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0); gStyle->SetEndErrorSize(0); gStyle->SetLineScalePS(1.5);
    if (parts < 0) { printf("sim_review: compiled only\n"); return; }
    SimData sd;
    if (!Load(sd, kCompact, (1 << kStdUS) | (1 << kLepUS))) return;
    scan::Data dd;
    if (!scan::Load(dd, scan::kCompact)) return;
    if (parts & 1) Part1(sd, dd);
    if (parts & 2) Part2(sd, dd);
    if (parts & 4) Part4(sd, dd);
    if (parts & 8) Part8(sd);
    if (parts & 16) Part16(sd, dd);
    printf("sim_review(%d): total %.0f s\n", parts, sw.RealTime());
}
