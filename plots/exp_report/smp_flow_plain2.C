// smp_flow_plain2.C -- COPY of plots/exp_report/smp_flow_cnt.C made by the short-report figure agent, redo shortfig_r (integ6, 2026-11-16):
// internal catalogue labels in the figure text replaced by plain names (glossary_v2.md); outputs redirected:
// figures of the scope -> docs/latex/<dir>/<base>_plain*.{pdf,png}; everything else -> /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ6/shortfig_r/run/smp
// Built by /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ6/shortfig_r/mk_plain_r.py; the physics code is unchanged.
// smp_flow_cnt.C -- COPY of plots/exp_report/smp_flow.C made by figure agent A (integ4, 2026-11-15) with
// the figure strings 'rows' -> 'counts' (analyst's terminology) and the outputs of the copy redirected:
// figures of the scope -> docs/latex/figs_exp/<base>_cnt.{pdf,png}; everything else -> /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figA/work/smp_flow
// Built by /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figA/work/mk_cnt.py; the physics code is unchanged.
// ============================================================================
// smp_flow.C  --  cut flow after the FAT production, rates per day and per
//                 sector, eta-peak content per period (data)
// ============================================================================
// HADES p+p, T = 4.53 GeV (feb22), channel pi+ pi- e+ e-. EXPERIMENTAL DATA ONLY.
// One row = one stored particle combination; event group = adjacent rows with
// equal (file_id, run, event) (scan_prepare.C, column grp).
//
// What it does
//   1. reads all rows of the compact file (trees us = e+e-, pp = e+e+, mm = e-e-),
//      and, row by row in the same order, the track sectors (ep/em/pip/pim_sector)
//      from scan_nt of the three FAT outputs (file_id and pat_entry are compared
//      row by row);
//   2. cut flow after FAT (cumulative steps, then the named selections, each
//      defined from pre), for isBest >= 0 and isBest == 1, without and with
//      OA(e+e-) > 2 deg: rows, event groups, rows in the data peak window
//      0.53364-0.55416, the data bump window 0.4550-0.5225, the frozen
//      simulation windows 0.5333-0.5587 and 0.4576-0.5255 GeV/c2, and for e+e- the
//      CB = 2 sqrt(N++ N--) in these windows (10 MeV bins interpolated, k = 1);
//   3. rows per day of every selection (variance = sum over event groups of k^2,
//      k = selected rows of the group) -> rows per nb^-1 x 64/63 per day, with the
//      luminosity of the analysed runs of the day (smp_common.h), luminosity-weighted
//      mean <R> = sum rows x 64/63 / sum L, chi2 test of
//      constancy, mean of days 034-044 (sector-2 trigger deficit) against 045-067;
//   4. share of the tracks per sector (0-5), leptons and pions, per day and per
//      period, after pre;
//   5. eta peak per period: P_main (chi2 fit of N - CB, frozen Gaussian 0.5439
//      GeV/c2 / 5.4 MeV + pol2 in 0.5225-0.6100) with the model spread (half range
//      over the 48 variations with p > 0.05) for pre and peak nominal + vertex
//      (also with OA > 2 and isBest == 1); S x 64/63 / L per period; chi2 test of
//      constancy over the four periods; free position per period (mass-scale
//      check); bump-window N - CB per nb^-1 (indicator: bump + physical background);
//   6. the four runs with no_hf = 1 (L = 0): rows and the effect of removing them.
// Selections (corrected momenta; already in the FAT production: isBest >= 0,
// vertex z > -500 mm, PT3 trigbit == 8192)
//   pre                  = cor_ok == 1 && isBest >= 0 && vz > -160 mm   (= frozen BUMP)
//   vertex               = pip_isUsedVertex == 1 && pim_isUsedVertex == 1
//   peak nominal         = pre + oa_lab < 55 deg + M(pi+pi-) > 0.30 GeV/c2
//   bump nominal         = pre + oa_lab < 40 deg + 0.30 < M(pi+pi-) < 0.41 GeV/c2
//   DIRECT (frozen, sim) = pre + dz_pi < 40 mm
//   PEAK   (frozen, sim) = pre + dz_pi < 40 mm + MM4 < 2.85 GeV/c2
//   X17 region           = pre + 0.5333 <= M4 < 0.5587 + 15.5 <= M(ee) < 24.5 MeV/c2 + vertex
//                          (also with the data peak window 0.53364-0.55416)
//
// Inputs : /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/scan_v2/scan_compact.root
//          /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/output_pippim{epem,epep,emem}_exp.root (scan_nt, sectors)
//          .../report_exp/smp/pat/runs_<day>.txt (smp_pat_counts.C) and .../report_exp/smp/lumi_inputs/*.csv
// Outputs: .../report_exp/smp/flow.txt (all numbers; machine-readable lines CF, DAY, SEC, PK, NOHF)
//          docs/latex/figs_exp/smp_rates_perday, smp_peak_periods, smp_peak_spectra_pre,
//          smp_peak_spectra_nomvtx, smp_sectors (.pdf + .png)
//          docs/latex/tables_exp/smp_rates_perday.txt, smp_peak_periods.txt/.tex, smp_sectors.txt/.tex
//
// Run (from the repository root /home/przygoda/HADES/FAT/FAT; about 3 GB of memory, 1-2 min):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/aclic_smp", kTRUE)' 'plots/exp_report/smp_flow.C+'
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <functional>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "TBox.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TDatime.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TPad.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include "../ScanEngine.h"
#include "smp_common_plain2.h"

namespace smpf {

using scan::PeakCfg;
using scan::PeakFit;
using scan::Row;
using scan::Sel;
using scan::Spectra;

FILE* gOut = nullptr;
void P(const char* fmt, ...) {
    va_list a; va_start(a, fmt); vprintf(fmt, a); va_end(a);
    if (gOut) { va_start(a, fmt); vfprintf(gOut, fmt, a); va_end(a); }
}

// ---- selection pieces ------------------------------------------------------------
inline bool fPre(const Row& r)     { return r.cor_ok == 1 && r.isBest >= 0 && r.vz > -160.f; }
inline bool fVtx(const Row& r)     { return r.pip_isUsedVertex == 1 && r.pim_isUsedVertex == 1; }
inline bool fDz(const Row& r)      { return r.dz_pi < 40.f; }
inline bool fMM(const Row& r)      { return r.mm4 < 2.85f; }
inline bool fPeakKin(const Row& r) { return r.oa_lab < 55.f && r.mpipi > 0.30f; }
inline bool fBumpKin(const Row& r) { return r.oa_lab < 40.f && r.mpipi > 0.30f && r.mpipi < 0.41f; }
inline bool fX17ee(const Row& r)   { const float m = r.mee * 1000.f; return m >= 15.5f && m < 24.5f; }
const double kSimPkLo = 0.5333, kSimPkHi = 0.5587;     // frozen simulation windows
const double kSimBuLo = 0.4576, kSimBuHi = 0.5255;
inline bool fWinSim(const Row& r)  { const double m = r.m4; return m >= kSimPkLo && m < kSimPkHi; }
inline bool fWinData(const Row& r) { const double m = r.m4; return m >= scan::kPeakLo && m < scan::kPeakHi; }

struct Step { std::string key, name; Sel f; };
std::vector<Step> Steps() {
    std::vector<Step> v;
    v.push_back({"fat",    "FAT rows (isBest>=0, vz>-500, PT3)",      [](const Row&) { return true; }});
    v.push_back({"cor",    "+ cor_ok == 1",                           [](const Row& r) { return r.cor_ok == 1; }});
    v.push_back({"pre",    "+ vz > -160 mm  (= pre = BUMP)",          [](const Row& r) { return fPre(r); }});
    v.push_back({"vtx",    "+ both pions in the vertex fit",          [](const Row& r) { return fPre(r) && fVtx(r); }});
    v.push_back({"vtxdz",  "+ dz_pi < 40 mm",                         [](const Row& r) { return fPre(r) && fVtx(r) && fDz(r); }});
    v.push_back({"direct", "DIRECT = pre + dz_pi < 40",               [](const Row& r) { return fPre(r) && fDz(r); }});
    v.push_back({"peakf",  "PEAK = pre + dz_pi < 40 + MM4 < 2.85",    [](const Row& r) { return fPre(r) && fDz(r) && fMM(r); }});
    v.push_back({"pnom",   "peak nominal",                            [](const Row& r) { return fPre(r) && fPeakKin(r); }});
    v.push_back({"pnomv",  "peak nominal + vertex",                   [](const Row& r) { return fPre(r) && fPeakKin(r) && fVtx(r); }});
    v.push_back({"bnom",   "bump nominal",                            [](const Row& r) { return fPre(r) && fBumpKin(r); }});
    v.push_back({"bnomv",  "bump nominal + vertex",                   [](const Row& r) { return fPre(r) && fBumpKin(r) && fVtx(r); }});
    v.push_back({"x17s",   "X17 region (sim window 0.5333-0.5587)",   [](const Row& r) { return fPre(r) && fWinSim(r) && fX17ee(r) && fVtx(r); }});
    v.push_back({"x17d",   "X17 region (data window 0.53364-0.55416)",[](const Row& r) { return fPre(r) && fWinData(r) && fX17ee(r) && fVtx(r); }});
    return v;
}
const int kNVar = 4;
const char* kVarName[kNVar] = {"isBest>=0", "isBest==1", "isBest>=0, OA>2", "isBest==1, OA>2"};
inline Sel WithVar(const Sel& s, int v) {
    switch (v) {
        case 1: return [s](const Row& r) { return r.isBest == 1 && s(r); };
        case 2: return [s](const Row& r) { return r.oa > 2.f && s(r); };
        case 3: return [s](const Row& r) { return r.isBest == 1 && r.oa > 2.f && s(r); };
        default: return s;
    }
}

// ---- period groups of the peak fits ------------------------------------------------
const int kNG = 8;
const char* kGName[kNG] = {"I", "IIA", "IIB", "III", "034-044", "045-046", "all", "045-067"};
const char* kGLong[kNG] = {"I (034-046)", "IIA (047-048)", "IIB (049-053)", "III (059-067)",
                           "034-044 (sector-2 deficit)", "045-046", "all 29 days", "045-067"};
inline bool InGroup(int g, int day) {
    if (g < 4) return smp::PeriodOf(day) == g;
    if (g == 4) return day >= 34 && day <= 44;
    if (g == 5) return day >= 45 && day <= 46;
    if (g == 6) return true;
    return day >= 45;
}

// ---- no_hf runs (L = 0) that are in the data ------------------------------------------
const std::set<int> kNoHF = {1244739, 3274681, 2109357, 2109532};

// ---- loading all rows + sectors ---------------------------------------------------------
// sector code per row: 3 bits per track (ep, em, pip, pim), 7 = invalid
bool LoadAll(scan::Data& d, std::vector<unsigned short> sec[3]) {
    TStopwatch sw;
    TFile* f = TFile::Open(scan::kCompact, "READ");
    if (!f || f->IsZombie()) return false;
    const char* fat[3] = {"output_pippimepem_exp.root", "output_pippimepep_exp.root", "output_pippimemem_exp.root"};
    for (int k = 0; k < 3; ++k) {
        TTree* t = (TTree*)f->Get(scan::SampleName(k));
        Row r;
        if (scan::BindRow(t, r, false) > 0) { P("compact file without the pion z columns\n"); return false; }
        TFile* g = TFile::Open(Form("%s/%s", scan::kProdDir, fat[k]), "READ");
        TTree* s = (TTree*)g->Get("scan_nt");
        const Long64_t n = t->GetEntries();
        if (s->GetEntries() != n) { P("row count mismatch %s\n", fat[k]); return false; }
        s->SetBranchStatus("*", 0);
        Float_t sv[4], fid = 0, pe = 0;
        const char* sb[4] = {"ep_sector", "em_sector", "pip_sector", "pim_sector"};
        for (int q = 0; q < 4; ++q) { s->SetBranchStatus(sb[q], 1); s->SetBranchAddress(sb[q], &sv[q]); }
        s->SetBranchStatus("file_id", 1); s->SetBranchAddress("file_id", &fid);
        s->SetBranchStatus("pat_entry", 1); s->SetBranchAddress("pat_entry", &pe);
        d.s[k].clear(); d.s[k].reserve(n);
        sec[k].clear(); sec[k].reserve(n);
        Long64_t mism = 0, bad = 0;
        for (Long64_t i = 0; i < n; ++i) {
            t->GetEntry(i); s->GetEntry(i);
            if ((int)std::lround(fid) != r.file_id || (int)std::lround(pe) != r.pat_entry) mism++;
            unsigned short code = 0;
            for (int q = 0; q < 4; ++q) {
                int sc = (int)std::lround(sv[q]);
                if (sc < 0 || sc > 5) { sc = 7; bad++; }
                code |= (unsigned short)(sc << (3 * q));
            }
            d.s[k].push_back(r);
            sec[k].push_back(code);
        }
        t->ResetBranchAddresses();
        g->Close(); delete g;
        P("load %s: %lld rows, (file_id, pat_entry) mismatches with scan_nt %lld, invalid sectors %lld\n",
          scan::SampleName(k), n, mism, bad);
        if (mism) return false;
    }
    f->Close(); delete f;
    d.fallback = false;
    P("load: %.1f s\n", sw.RealTime());
    return true;
}

// ---- per-day accumulators -----------------------------------------------------------
const int kMaxDay = 68;
struct DayAcc { double rows[kMaxDay], var[kMaxDay], ev[kMaxDay]; DayAcc() { for (int i = 0; i < kMaxDay; ++i) rows[i] = var[i] = ev[i] = 0; } };

struct Cnt { double rows = 0, ev = 0, w[4] = {0, 0, 0, 0}; };

// weighted mean and chi2 of constancy
struct Const { double mean = 0, emean = 0, chi2 = 0; int ndf = 0; double prob = 0, rmsRel = 0; };
Const Constancy(const std::vector<double>& x, const std::vector<double>& e) {
    Const c; double sw = 0, swx = 0; int n = 0;
    for (size_t i = 0; i < x.size(); ++i) if (e[i] > 0) { sw += 1 / (e[i] * e[i]); swx += x[i] / (e[i] * e[i]); n++; }
    if (n == 0) return c;
    c.mean = swx / sw; c.emean = 1 / std::sqrt(sw);
    double s2 = 0;
    for (size_t i = 0; i < x.size(); ++i) if (e[i] > 0) {
        c.chi2 += (x[i] - c.mean) * (x[i] - c.mean) / (e[i] * e[i]);
        s2 += (x[i] / c.mean - 1) * (x[i] / c.mean - 1);
    }
    c.ndf = n - 1; c.prob = c.ndf > 0 ? TMath::Prob(c.chi2, c.ndf) : 1; c.rmsRel = std::sqrt(s2 / n);
    return c;
}
// Weighted mean with external weights w (verification fix 2026-11-14): for the per-day rates
// R_i = N_i x 64/63 / L_i the weights w_i = L_i give mean = sum N x 64/63 / sum L, the maximum-
// likelihood constant rate. The stat-weighted mean (Constancy) gives a weight L_i/R_i to each day
// and is pulled towards the low days when the scatter is far above the statistics (chi2/ndf
// 40-250 here): -0.3 % for e+e-, -2 % for the like-sign rates. chi2 and rms are taken about
// this mean, with the statistical errors e.
Const ConstancyW(const std::vector<double>& x, const std::vector<double>& e, const std::vector<double>& w) {
    Const c; double sw = 0, swx = 0, sw2e2 = 0; int n = 0;
    for (size_t i = 0; i < x.size(); ++i) if (e[i] > 0 && w[i] > 0) { sw += w[i]; swx += w[i] * x[i]; sw2e2 += w[i] * w[i] * e[i] * e[i]; n++; }
    if (n == 0) return c;
    c.mean = swx / sw; c.emean = std::sqrt(sw2e2) / sw;
    double s2 = 0;
    for (size_t i = 0; i < x.size(); ++i) if (e[i] > 0 && w[i] > 0) {
        c.chi2 += (x[i] - c.mean) * (x[i] - c.mean) / (e[i] * e[i]);
        s2 += (x[i] / c.mean - 1) * (x[i] / c.mean - 1);
    }
    c.ndf = n - 1; c.prob = c.ndf > 0 ? TMath::Prob(c.chi2, c.ndf) : 1; c.rmsRel = std::sqrt(s2 / n);
    return c;
}

// draw helpers
void PeriodLines(TPad* p, double ylo, double yhi, double ylab, double tsize) {
    p->cd();
    TBox* b = new TBox(33.5, ylo, 44.5, yhi);
    b->SetFillColor(TColor::GetColor("#e3e3e3")); b->SetFillStyle(1001); b->SetLineWidth(0);
    b->Draw();
    const double xs[3] = {46.5, 48.5, 56.0};
    for (double x : xs) {
        TLine* l = new TLine(x, ylo, x, yhi); l->SetLineStyle(2); l->SetLineColor(kGray + 2); l->Draw();
    }
    if (tsize > 0) {
        TLatex t; t.SetTextSize(tsize); t.SetTextAlign(21); t.SetTextColor(kGray + 3);
        t.DrawLatex(40.0, ylab, "days 34#minus46"); t.SetTextSize(0.58 * tsize); t.DrawLatex(47.5, ylab, "47#minus48"); t.SetTextSize(tsize); t.DrawLatex(51.0, ylab, "49#minus53"); t.DrawLatex(63.0, ylab, "59#minus67");
    }
}

} // namespace smpf

void smp_flow_plain2()
{
    using namespace smpf;
    TStopwatch swAll;
    smp::Style();
    gSystem->mkdir(smp::kWork, kTRUE); gSystem->mkdir(smp::kWorkOut, kTRUE); gSystem->mkdir(smp::kPlainScratch, kTRUE);
    gSystem->mkdir(smp::kFigOut, kTRUE);
    gSystem->mkdir(smp::kTabOut, kTRUE);
    gOut = fopen(Form("%s/flow.txt", smp::kWorkOut), "w");
    P("# smp_flow.C  (%s)\n", TDatime().AsString());
    P("# data: %s\n", scan::kCompact);

    // ---- luminosity of the analysed runs per day ----
    std::map<int, double> Lday = smp::AnalysedLumiPerDay(false);
    std::map<int, smp::DayLumi> dl = smp::ReadDayLumi();
    double Lsum = 0; for (auto& p : Lday) Lsum += p.second;
    P("# L6_SF of the analysed runs (EpEm_ID runs present): %.3f nb^-1\n", Lsum);
    double Lg[kNG];
    for (int g = 0; g < kNG; ++g) { Lg[g] = 0; for (auto& p : Lday) if (InGroup(g, p.first)) Lg[g] += p.second; }

    // ---- load ----
    scan::Data d;
    std::vector<unsigned short> sec[3];
    if (!LoadAll(d, sec)) { P("load failed\n"); fclose(gOut); return; }

    // ---- selections: 13 steps x 4 variants ----
    const std::vector<Step> st = Steps();
    const int nst = (int)st.size();
    const int ns = nst * kNVar;
    std::vector<Sel> sels(ns);
    for (int v = 0; v < kNVar; ++v) for (int i = 0; i < nst; ++i) sels[v * nst + i] = WithVar(st[i].f, v);

    std::vector<Cnt> cnt[3];
    std::vector<DayAcc> day[3];
    // sector counts after pre (isBest >= 0): [k][day][0 leptons / 1 pions][sector]
    std::vector<double> secC[3];
    // rows after pre with at least one track / one lepton / one pion in sector 2: [k][day][0 rows, 1 any, 2 lepton, 3 pion]
    std::vector<double> secRow[3];
    for (int k = 0; k < 3; ++k) { cnt[k].assign(ns, Cnt()); day[k].assign(ns, DayAcc()); secC[k].assign(kMaxDay * 2 * 6, 0.); secRow[k].assign(kMaxDay * 4, 0.); }
    // no_hf runs
    double nohfRows[3][2] = {{0, 0}, {0, 0}, {0, 0}};       // FAT rows, pre rows
    std::map<int, double> nohfPerRun;                       // e+e- FAT rows per no_hf run

    const int iPre = 2;
    TStopwatch swC;
    for (int k = 0; k < 3; ++k) {
        const std::vector<Row>& rows = d.s[k];
        int cur = -1, curDay = 0;
        std::vector<int> gcnt(ns, 0);
        auto flush = [&]() {
            if (cur < 0) return;
            for (int s = 0; s < ns; ++s) if (gcnt[s]) {
                cnt[k][s].ev += 1;
                day[k][s].ev[curDay] += 1;
                day[k][s].var[curDay] += (double)gcnt[s] * gcnt[s];
                gcnt[s] = 0;
            }
        };
        for (size_t i = 0; i < rows.size(); ++i) {
            const Row& r = rows[i];
            if (r.grp != cur) { flush(); cur = r.grp; curDay = r.file_id; }
            const double m = r.m4;
            const bool w0 = m >= scan::kPeakLo && m < scan::kPeakHi, w1 = m >= scan::kBumpLo && m < scan::kBumpHi;
            const bool w2 = m >= kSimPkLo && m < kSimPkHi, w3 = m >= kSimBuLo && m < kSimBuHi;
            for (int s = 0; s < ns; ++s) {
                if (!sels[s](r)) continue;
                Cnt& c = cnt[k][s];
                c.rows += 1; if (w0) c.w[0] += 1; if (w1) c.w[1] += 1; if (w2) c.w[2] += 1; if (w3) c.w[3] += 1;
                day[k][s].rows[r.file_id] += 1;
                gcnt[s]++;
            }
            if (fPre(r)) {
                const unsigned short code = sec[k][i];
                bool s2l = false, s2p = false;
                for (int q = 0; q < 4; ++q) {
                    const int sc = (code >> (3 * q)) & 7;
                    if (sc > 5) continue;
                    secC[k][(r.file_id * 2 + (q < 2 ? 0 : 1)) * 6 + sc] += 1;
                    if (sc == 2) { if (q < 2) s2l = true; else s2p = true; }
                }
                double* sr = &secRow[k][r.file_id * 4];
                sr[0] += 1; if (s2l || s2p) sr[1] += 1; if (s2l) sr[2] += 1; if (s2p) sr[3] += 1;
            }
            if (kNoHF.count(r.run)) {
                nohfRows[k][0] += 1;
                if (fPre(r)) nohfRows[k][1] += 1;
                if (k == 0) nohfPerRun[r.run] += 1;
            }
        }
        flush();
    }
    P("# counting: %.1f s\n", swC.RealTime());

    // ---- CB in the windows (FillMany; rows with cor_ok == 1) ----
    scan::FillOpt fo;
    fo.windows.push_back(std::make_pair(kSimPkLo, kSimPkHi));
    fo.windows.push_back(std::make_pair(kSimBuLo, kSimBuHi));
    // The X17 selections contain an M4 window; the like-sign 10 MeV bins at the window edges would be
    // only partly filled and the interpolated CB too low (verification fix 2026-11-14: CB/N 0.018 -> 0.019
    // in the data window). The CB of the X17 rows is therefore formed from the same selection WITHOUT
    // the M4 cut and integrated over the window (valid for the windows inside the M4 cut of the row).
    std::vector<Sel> selsCB = sels;
    {
        const Sel x17noM4 = [](const Row& r) { return fPre(r) && fX17ee(r) && fVtx(r); };
        for (int v = 0; v < kNVar; ++v)
            for (int i = 0; i < nst; ++i)
                if (st[i].key == "x17s" || st[i].key == "x17d") selsCB[v * nst + i] = WithVar(x17noM4, v);
    }
    std::vector<Spectra> spCF = scan::FillMany(d, selsCB, fo);

    // ---- print the cut flow ----
    P("\n==== A. CUT FLOW AFTER FAT (rows = stored combinations; events = event groups) ====\n");
    P("# windows: pk = data peak 0.53364-0.55416, bu = data bump 0.4550-0.5225, pkS/buS = frozen sim windows 0.5333-0.5587 / 0.4576-0.5255 GeV/c2\n");
    P("# CB = 2 sqrt(N++ N--) per 10 MeV bin, interpolated (k = 1); CB/N given for e+e- rows; '-' = not defined (rows with cor_ok == 0 included)\n");
    const char* hn[3] = {"e+e-", "e+e+", "e-e-"};
    for (int v = 0; v < kNVar; ++v) {
        P("\n-- variant: %s --\n", kVarName[v]);
        P("%-42s | %9s %9s %8s %8s %8s %8s | %7s %7s | %8s %8s | %8s %8s\n", "step", "rows e+e-", "ev e+e-", "N pk", "N bu", "N pkS", "N buS",
          "CB/N pk", "CB/N bu", "rows e+e+", "ev e+e+", "rows e-e-", "ev e-e-");
        for (int i = 0; i < nst; ++i) {
            const int s = v * nst + i;
            const Cnt& u = cnt[0][s];
            const Spectra& sp = spCF[s];
            char c0[16], c1[16];
            const bool x17 = (st[i].key == "x17s" || st[i].key == "x17d");
            if (i == 0) { snprintf(c0, 16, "-"); snprintf(c1, 16, "-"); }
            else {
                // CB/N = CB / rows of the step in the window (for X17 the CB comes from the selection without the M4 cut)
                snprintf(c0, 16, "%.4f", u.w[0] > 0 ? sp.win[0].cb / u.w[0] : 0.);
                if (x17) snprintf(c1, 16, "-");
                else snprintf(c1, 16, "%.4f", u.w[1] > 0 ? sp.win[1].cb / u.w[1] : 0.);
            }
            P("%-42s | %9.0f %9.0f %8.0f %8.0f %8.0f %8.0f | %7s %7s | %8.0f %8.0f | %8.0f %8.0f\n", st[i].name.c_str(), u.rows, u.ev,
              u.w[0], u.w[1], u.w[2], u.w[3], c0, c1, cnt[1][s].rows, cnt[1][s].ev, cnt[2][s].rows, cnt[2][s].ev);
        }
    }
    // machine-readable
    for (int v = 0; v < kNVar; ++v)
        for (int i = 0; i < nst; ++i) {
            const int s = v * nst + i;
            for (int k = 0; k < 3; ++k) {
                const Cnt& c = cnt[k][s];
                const Spectra& sp = spCF[s];
                P("CF %d %d %s %d %.0f %.0f %.0f %.0f %.0f %.0f %.3f %.3f %.3f %.3f\n", v, i, st[i].key.c_str(), k, c.rows, c.ev,
                  c.w[0], c.w[1], c.w[2], c.w[3],
                  (i > 0 && k == 0) ? sp.win[0].cb : -1., (i > 0 && k == 0) ? sp.win[1].cb : -1.,
                  (i > 0 && k == 0) ? sp.win[2].cb : -1., (i > 0 && k == 0) ? sp.win[3].cb : -1.);
            }
        }

    // ---- B. rows per nb^-1 per day ----
    P("\n==== B. ROWS PER nb^-1 PER DAY (x 64/63), selection pre; L = L6_SF of the analysed runs of the day ====\n");
    const std::vector<int> days = smp::Days();
    struct RateSeries { std::vector<double> r, e, rb, eb; };
    RateSeries rs[3];
    const int sPre0 = 0 * nst + iPre, sPre1 = 1 * nst + iPre;
    FILE* tr = fopen(Form("%s/smp_rates_perday.txt", smp::kTabOut), "w");
    fprintf(tr, "# smp_rates_perday.txt -- written by plots/exp_report/smp_flow.C (%s)\n", TDatime().AsString());
    fprintf(tr, "# rows per day of the selection pre (cor_ok, isBest >= 0, vz > -160 mm; FAT: vz > -500, PT3) and the visible rate\n");
    fprintf(tr, "# R = rows x 64/63 / L [rows per nb^-1]; L = L6_SF summed over the analysed runs of the day (START-free, revision N)\n");
    fprintf(tr, "# error of R: rows with the event-group correlation (sqrt of sum k^2) (+) L stat; isBest==1 rows = events\n");
    fprintf(tr, "# sec2 = 1 for days 034-044 (sector-2 hardware trigger deficit, not in L)\n");
    fprintf(tr, "%-4s %-4s %4s %9s | %8s %9s %8s %7s | %7s %8s %7s %6s | %7s %8s %7s %6s\n", "day", "per", "sec2", "L[nb-1]",
            "us rows", "R_us", "err", "R_us(b1)", "pp rows", "R_pp", "err", "R(b1)", "mm rows", "R_mm", "err", "R(b1)");
    for (int dd : days) {
        const double L = Lday[dd];
        const double lstat = dl.count(dd) ? dl[dd].L6stat / std::max(dl[dd].L6, 1e-9) : 0;
        fprintf(tr, "%-4d %-4s %4d %9.3f |", dd, smp::PeriodName(smp::PeriodOf(dd)), smp::Sec2Deficit(dd) ? 1 : 0, L);
        for (int k = 0; k < 3; ++k) {
            const double n0 = day[k][sPre0].rows[dd], v0 = day[k][sPre0].var[dd], n1 = day[k][sPre1].rows[dd];
            const double R = L > 0 ? n0 * smp::kPT3fac / L : 0;
            const double eR = L > 0 ? std::hypot(std::sqrt(v0) * smp::kPT3fac / L, R * lstat) : 0;
            const double Rb = L > 0 ? n1 * smp::kPT3fac / L : 0;
            const double eRb = L > 0 ? std::hypot(std::sqrt(n1) * smp::kPT3fac / L, Rb * lstat) : 0;
            rs[k].r.push_back(R); rs[k].e.push_back(eR); rs[k].rb.push_back(Rb); rs[k].eb.push_back(eRb);
            fprintf(tr, " %8.0f %9.3f %8.4f %7.3f |", n0, R, eR, Rb);
            P("DAY %d %d %.0f %.1f %.0f %.4f %.5f %.4f %.5f\n", dd, k, n0, v0, n1, R, eR, Rb, eRb);
        }
        fprintf(tr, "\n");
    }
    // constancy (luminosity-weighted mean <R> = sum rows x 64/63 / sum L; see ConstancyW)
    fprintf(tr, "\n# constancy over the 29 days: <R> = luminosity-weighted mean = sum rows x 64/63 / sum L; chi2/ndf about <R> with the stat errors;\n");
    fprintf(tr, "# rms of R/<R> over the days; (stat-weighted mean in brackets, for reference only: biased low when the scatter exceeds the statistics)\n");
    P("\n# constancy of R over days (luminosity-weighted mean)\n");
    Const cR[3], cRb[3];
    std::vector<double> Lv;
    for (int dd : days) Lv.push_back(Lday[dd]);
    for (int k = 0; k < 3; ++k) {
        cR[k] = ConstancyW(rs[k].r, rs[k].e, Lv); cRb[k] = ConstancyW(rs[k].rb, rs[k].eb, Lv);
        const Const cS = Constancy(rs[k].r, rs[k].e);
        std::vector<double> ra, ea, la, rb, eb, lb;
        for (size_t i = 0; i < days.size(); ++i) {
            if (smp::Sec2Deficit(days[i])) { ra.push_back(rs[k].r[i]); ea.push_back(rs[k].e[i]); la.push_back(Lv[i]); }
            else { rb.push_back(rs[k].r[i]); eb.push_back(rs[k].e[i]); lb.push_back(Lv[i]); }
        }
        const Const ca = ConstancyW(ra, ea, la), cb = ConstancyW(rb, eb, lb);
        const char* line = Form("%s: <R> = %.3f +- %.4f rows/nb^-1 (stat-weighted %.3f), chi2/ndf = %.0f/%d (p %.1e), rms(R/<R>) = %.3f ; "
                                "isBest==1: <R> = %.3f, chi2/ndf %.0f/%d, rms %.3f ; days 034-044: %.3f (chi2 %.0f/%d, rms %.3f), "
                                "045-067: %.3f (chi2 %.0f/%d, rms %.3f), ratio %.4f",
                                hn[k], cR[k].mean, cR[k].emean, cS.mean, cR[k].chi2, cR[k].ndf, cR[k].prob, cR[k].rmsRel,
                                cRb[k].mean, cRb[k].chi2, cRb[k].ndf, cRb[k].rmsRel, ca.mean, ca.chi2, ca.ndf, ca.rmsRel,
                                cb.mean, cb.chi2, cb.ndf, cb.rmsRel, ca.mean / cb.mean);
        fprintf(tr, "# %s\n", line);
        P("%s\n", line);
    }
    // like-sign / unlike-sign per day; mean = ratio of the sums of rows (weights = e+e- rows x 64/63)
    {
        std::vector<double> x, e, w;
        for (size_t i = 0; i < days.size(); ++i) {
            const double a = rs[1].r[i] + rs[2].r[i], b = rs[0].r[i];
            x.push_back(a / b);
            e.push_back(a / b * std::hypot(std::hypot(rs[1].e[i], rs[2].e[i]) / a, rs[0].e[i] / b));
            w.push_back(b * Lv[i]);
        }
        const Const c = ConstancyW(x, e, w);
        fprintf(tr, "# (e+e+ + e-e-)/e+e- rows per day: mean (ratio of the sums) %.4f, chi2/ndf %.0f/%d, rms %.3f\n", c.mean, c.chi2, c.ndf, c.rmsRel);
        P("# (e+e+ + e-e-)/e+e- rows per day: mean (ratio of the sums) %.4f, chi2/ndf %.0f/%d, rms %.3f\n", c.mean, c.chi2, c.ndf, c.rmsRel);
    }
    fclose(tr);

    // ---- C. sectors ----
    P("\n==== C. SHARE OF TRACKS PER SECTOR, selection pre (isBest >= 0) ====\n");
    FILE* ts = fopen(Form("%s/smp_sectors.txt", smp::kTabOut), "w");
    fprintf(ts, "# smp_sectors.txt -- written by plots/exp_report/smp_flow.C (%s)\n", TDatime().AsString());
    fprintf(ts, "# share of the tracks per sector (0-5) after pre; leptons = both lepton slots, pions = pi+ and pi- slots;\n");
    fprintf(ts, "# rows = stored combinations (a track shared by several rows counts in each); binomial errors on the track counts\n");
    const int kNSG = 6;
    const char* sgName[kNSG] = {"034-044", "045-046", "IIA", "IIB", "III", "all"};
    auto inSG = [](int g, int dd) {
        if (g == 0) return dd >= 34 && dd <= 44;
        if (g == 1) return dd >= 45 && dd <= 46;
        if (g == 2) return smp::PeriodOf(dd) == 1;
        if (g == 3) return smp::PeriodOf(dd) == 2;
        if (g == 4) return smp::PeriodOf(dd) == 3;
        return true;
    };
    double shareSG[3][kNSG][2][6], eshareSG[3][kNSG][2][6];
    for (int k = 0; k < 3; ++k) {
        fprintf(ts, "\n# %s rows\n%-8s %-7s %8s | %15s %15s %15s %15s %15s %15s\n", hn[k], "days", "tracks", "N", "sector 0", "sector 1",
                "sector 2", "sector 3", "sector 4", "sector 5");
        for (int g = 0; g < kNSG; ++g)
            for (int t = 0; t < 2; ++t) {
                double c[6] = {0, 0, 0, 0, 0, 0}, tot = 0;
                for (int dd : days) if (inSG(g, dd)) for (int q = 0; q < 6; ++q) c[q] += secC[k][(dd * 2 + t) * 6 + q];
                for (int q = 0; q < 6; ++q) tot += c[q];
                fprintf(ts, "%-8s %-7s %8.0f |", sgName[g], t == 0 ? "leptons" : "pions", tot);
                for (int q = 0; q < 6; ++q) {
                    const double p = tot > 0 ? c[q] / tot : 0, e = tot > 0 ? std::sqrt(p * (1 - p) / tot) : 0;
                    shareSG[k][g][t][q] = p; eshareSG[k][g][t][q] = e;
                    fprintf(ts, " %8.4f+-%.4f", p, e);
                }
                fprintf(ts, "\n");
            }
    }
    // per day, e+e-
    fprintf(ts, "\n# per day, e+e- rows: share of sector 0..5, leptons | pions\n");
    std::vector<double> shDay[2][6];
    for (int dd : days) {
        fprintf(ts, "%03d", dd);
        for (int t = 0; t < 2; ++t) {
            double tot = 0; for (int q = 0; q < 6; ++q) tot += secC[0][(dd * 2 + t) * 6 + q];
            fprintf(ts, " |");
            for (int q = 0; q < 6; ++q) {
                const double p = tot > 0 ? secC[0][(dd * 2 + t) * 6 + q] / tot : 0;
                shDay[t][q].push_back(p);
                fprintf(ts, " %.4f", p);
                P("SEC %d %d %d %.0f %.0f\n", dd, t, q, secC[0][(dd * 2 + t) * 6 + q], tot);
            }
        }
        fprintf(ts, "\n");
    }
    // rows with a track in sector 2 and the implied loss on the deficit days
    fprintf(ts, "\n# rows (pre) with at least one track in sector 2: fraction per day (any track | a lepton | a pion), e+e-, e+e+, e-e-\n");
    for (int dd : days) {
        fprintf(ts, "%03d", dd);
        for (int k = 0; k < 3; ++k) {
            const double* sr = &secRow[k][dd * 4];
            fprintf(ts, " | %.4f %.4f %.4f", sr[1] / sr[0], sr[2] / sr[0], sr[3] / sr[0]);
            P("SECROW %d %d %.0f %.0f %.0f %.0f\n", k, dd, sr[0], sr[1], sr[2], sr[3]);
        }
        fprintf(ts, "\n");
    }
    // loss estimate: rows without a sector-2 track are assumed unaffected; the fraction of rows with a
    // sector-2 track expected on the deficit days is taken from the reference days
    fprintf(ts, "\n# loss estimate on the deficit days: expected sector-2 rows = f_ref/(1-f_ref) x (rows without sector-2 track);\n");
    fprintf(ts, "# loss = (expected - observed)/(rows + expected - observed); assumes rows without a sector-2 track are unaffected\n");
    const int kNL = 2;
    const int defLo[kNL] = {34, 34}, defHi[kNL] = {44, 45}, refLo[kNL] = {45, 46};
    for (int k = 0; k < 3; ++k)
        for (int L = 0; L < kNL; ++L) {
            double a[4] = {0, 0, 0, 0}, rf[4] = {0, 0, 0, 0};
            for (int dd : days) {
                const double* sr = &secRow[k][dd * 4];
                for (int q = 0; q < 4; ++q) { if (dd >= defLo[L] && dd <= defHi[L]) a[q] += sr[q]; else if (dd >= refLo[L]) rf[q] += sr[q]; }
            }
            const double fref = rf[1] / rf[0], fobs = a[1] / a[0];
            const double expd = fref / (1 - fref) * (a[0] - a[1]);
            const double loss = (expd - a[1]) / (a[0] + expd - a[1]);
            const double frl = rf[2] / rf[0], fol = a[2] / a[0];
            const double expl = frl / (1 - frl) * (a[0] - a[2]);
            const double lossl = (expl - a[2]) / (a[0] + expl - a[2]);
            const char* line = Form("%s: deficit days %03d-%03d vs reference %03d-067: rows with a sector-2 track %.4f vs %.4f -> row loss %.4f "
                                    "(lepton in sector 2: %.4f vs %.4f -> %.4f)", hn[k], defLo[L], defHi[L], refLo[L], fobs, fref, loss, fol, frl, lossl);
            fprintf(ts, "# %s\n", line);
            P("%s\n", line);
        }
    fclose(ts);
    for (int k = 0; k < 3; ++k)
        P("%s: sector-2 share leptons 034-044 %.4f vs 045-067 (IIA,IIB,III) %.4f / %.4f / %.4f; pions %.4f vs %.4f / %.4f / %.4f\n", hn[k],
          shareSG[k][0][0][2], shareSG[k][2][0][2], shareSG[k][3][0][2], shareSG[k][4][0][2],
          shareSG[k][0][1][2], shareSG[k][2][1][2], shareSG[k][3][1][2], shareSG[k][4][1][2]);
    {   // LaTeX: e+e- shares per group
        FILE* o = fopen(Form("%s/smp_sectors.tex", smp::kTabOut), "w");
        fprintf(o, "%% smp_sectors.tex -- written by plots/exp_report/smp_flow.C (numbers: tables_exp/smp_sectors.txt)\n");
        fprintf(o, "\\begin{table}[H]\n\\centering\\small\n\\setlength{\\tabcolsep}{4pt}\n");
        fprintf(o, "\\caption{Share of the tracks per sector (in \\%%) in the $e^+e^-$ rows of the selection pre (all stored combinations, "
                   "isBest $\\geq 0$, corrected momenta, PT3, $v_z>-160$~mm), for the two lepton slots and the two pion slots. "
                   "Days 034--044 carry the sector-2 hardware trigger deficit, which is not contained in the luminosity. "
                   "Binomial errors on the track counts are 0.02--0.1 percentage points.}\n");
        fprintf(o, "\\label{tab:smp_sectors}\n\\begin{tabular}{llrrrrrr}\n\\toprule\n");
        fprintf(o, "Days & Tracks & S0 & S1 & S2 & S3 & S4 & S5 \\\\\n\\midrule\n");
        for (int t = 0; t < 2; ++t) {
            for (int g = 0; g < kNSG; ++g) {
                fprintf(o, "%s & %s", sgName[g], t == 0 ? "leptons" : "pions");
                for (int q = 0; q < 6; ++q) fprintf(o, " & %.2f", 100 * shareSG[0][g][t][q]);
                fprintf(o, " \\\\\n");
            }
            if (t == 0) fprintf(o, "\\midrule\n");
        }
        fprintf(o, "\\bottomrule\n\\end{tabular}\n\\end{table}\n");
        fclose(o);
    }

    // ---- D. eta peak per period ----
    P("\n==== D. ETA PEAK PER PERIOD (P_main, frozen 0.5439 GeV/c2 / 5.4 MeV + pol2, 0.5225-0.6100, N - CB) ====\n");
    const int kNB = 6;
    const char* bName[kNB] = {"pre", "peak nominal + vertex", "pre + OA>2", "peak nominal + vertex + OA>2",
                              "pre, isBest==1", "peak nominal + vertex, isBest==1"};
    const char* bKey[kNB] = {"pre", "pnomv", "pre_oa2", "pnomv_oa2", "pre_b1", "pnomv_b1"};
    std::vector<Sel> bsel(kNB);
    bsel[0] = st[2].f; bsel[1] = st[8].f;
    bsel[2] = WithVar(st[2].f, 2); bsel[3] = WithVar(st[8].f, 2);
    bsel[4] = WithVar(st[2].f, 1); bsel[5] = WithVar(st[8].f, 1);
    std::vector<Sel> ps;
    for (int b = 0; b < kNB; ++b)
        for (int g = 0; g < kNG; ++g) {
            const Sel base = bsel[b];
            ps.push_back([base, g](const Row& r) { return InGroup(g, r.file_id) && base(r); });
        }
    // no_hf removed (all days): pre and peak nominal + vertex
    for (int b = 0; b < 2; ++b) {
        const Sel base = bsel[b];
        ps.push_back([base](const Row& r) { return !kNoHF.count(r.run) && base(r); });
    }
    std::vector<Spectra> spP = scan::FillMany(d, ps);
    PeakFit pf[kNB][kNG];
    scan::Spread spr[kNB][kNG];
    scan::FreeFit ff[kNB][kNG];
    for (int b = 0; b < kNB; ++b)
        for (int g = 0; g < kNG; ++g) {
            const Spectra& sp = spP[b * kNG + g];
            pf[b][g] = scan::FitPeakMain(sp);
            spr[b][g] = scan::PeakModelSpread(sp);
            if (b < 2) ff[b][g] = scan::FitPeakFree(sp);
        }
    FILE* tp = fopen(Form("%s/smp_peak_periods.txt", smp::kTabOut), "w");
    fprintf(tp, "# smp_peak_periods.txt -- written by plots/exp_report/smp_flow.C (%s)\n", TDatime().AsString());
    fprintf(tp, "# eta peak content per period, P_main: chi2 fit of N - CB (CB = 2 sqrt(N++ N--), 10 MeV bins interpolated, k = 1),\n");
    fprintf(tp, "# frozen Gaussian 0.5439 GeV/c2 / 5.4 MeV + pol2 in 0.5225-0.6100 GeV/c2; S = whole Gaussian; stat = fit error;\n");
    fprintf(tp, "# model = half range of S over the 48 variations with p > 0.05 (n = number kept); s = S x 64/63 / L [rows per nb^-1],\n");
    fprintf(tp, "# L = L6_SF of the analysed runs of the group (stat error of L negligible, 0.1 %%); visible yield, NOT corrected for\n");
    fprintf(tp, "# acceptance or efficiency. N_pk, CB_pk = rows and CB in the data peak window 0.53364-0.55416. free: position and\n");
    fprintf(tp, "# width fitted (mass-scale check). W_bu = N - CB in the data bump window 0.4550-0.5225 per nb^-1 (indicator only:\n");
    fprintf(tp, "# eta -> pi+pi-pi0 bump plus physical background).\n");
    struct PkLine { double s, es, ms; };
    PkLine pl[kNB][kNG];
    for (int b = 0; b < kNB; ++b) {
        fprintf(tp, "\n# selection: %s\n", bName[b]);
        fprintf(tp, "%-27s %8s | %7s %6s %6s %3s %6s %5s | %7s %7s %7s | %7s %6s | %10s %9s | %7s %6s\n", "group", "L[nb-1]", "S", "stat", "model",
                "n", "chi2", "ndf", "s/nb", "stat", "model", "N_pk", "CB/N", "free mean", "sigma", "W_bu/nb", "err");
        for (int g = 0; g < kNG; ++g) {
            const PeakFit& f = pf[b][g];
            const scan::Spread& s = spr[b][g];
            const Spectra& sp = spP[b * kNG + g];
            const double L = Lg[g];
            pl[b][g].s = f.S * smp::kPT3fac / L; pl[b][g].es = f.eS * smp::kPT3fac / L; pl[b][g].ms = s.half * smp::kPT3fac / L;
            const double W = sp.win[1].sig() * smp::kPT3fac / L, eW = sp.win[1].sigErr() * smp::kPT3fac / L;
            char fm[64] = "-";
            if (b < 2 && ff[b][g].ok) snprintf(fm, 64, "%.5f+-%.5f %5.2f+-%.2f", ff[b][g].mean, ff[b][g].eMean, 1000 * ff[b][g].sigma, 1000 * ff[b][g].eSigma);
            fprintf(tp, "%-27s %8.2f | %7.0f %6.0f %6.0f %3d %6.1f %5d | %7.4f %7.4f %7.4f | %7.0f %6.3f | %20s | %7.3f %6.3f\n", kGLong[g], L,
                    f.S, f.eS, s.half, s.npass, f.chi2, f.ndf, pl[b][g].s, pl[b][g].es, pl[b][g].ms, f.N, f.N > 0 ? f.CB / f.N : 0., fm, W, eW);
            P("PK %s %d %.3f %.2f %.2f %.2f %d %.2f %d %.5f %.5f %.5f %.0f %.1f %.4f %.4f\n", bKey[b], g, L, f.S, f.eS, s.half, s.npass, f.chi2, f.ndf,
              pl[b][g].s, pl[b][g].es, pl[b][g].ms, f.N, f.CB, W, eW);
        }
        // constancy over the four periods (stat) and 034-044 vs 045-067
        std::vector<double> x, e, eL;
        for (int g = 0; g < 4; ++g) { x.push_back(pl[b][g].s); e.push_back(pl[b][g].es); }
        const Const c = Constancy(x, e);
        const double dif = pl[b][4].s - pl[b][7].s, edif = std::hypot(pl[b][4].es, pl[b][7].es);
        fprintf(tp, "# constancy over I, IIA, IIB, III (stat errors): mean %.4f +- %.4f rows/nb^-1, chi2/ndf = %.2f/%d, p = %.3f\n",
                c.mean, c.emean, c.chi2, c.ndf, c.prob);
        fprintf(tp, "# 034-044 minus 045-067: %.4f +- %.4f (%.2f sigma); ratio %.3f +- %.3f\n", dif, edif, dif / edif,
                pl[b][4].s / pl[b][7].s, pl[b][4].s / pl[b][7].s * std::hypot(pl[b][4].es / pl[b][4].s, pl[b][7].es / pl[b][7].s));
        P("# %s: constancy over periods chi2/ndf = %.2f/%d p = %.3f, mean %.4f +- %.4f; 034-044 - 045-067 = %.4f +- %.4f\n", bName[b],
          c.chi2, c.ndf, c.prob, c.mean, c.emean, dif, edif);
    }
    // no_hf
    {
        const Spectra& a0 = spP[kNB * kNG + 0];
        const Spectra& a1 = spP[kNB * kNG + 1];
        const PeakFit q0 = scan::FitPeakMain(a0), q1 = scan::FitPeakMain(a1);
        fprintf(tp, "\n# the four no_hf runs (L = 0) removed (all days): pre S = %.1f +- %.1f (with them %.1f), peak nominal + vertex S = %.1f +- %.1f (with them %.1f)\n",
                q0.S, q0.eS, pf[0][6].S, q1.S, q1.eS, pf[1][6].S);
        P("NOHF fit pre %.2f %.2f %.2f pnomv %.2f %.2f %.2f\n", q0.S, q0.eS, pf[0][6].S, q1.S, q1.eS, pf[1][6].S);
        P("NOHF rows FAT us %.0f pp %.0f mm %.0f ; pre us %.0f pp %.0f mm %.0f\n", nohfRows[0][0], nohfRows[1][0], nohfRows[2][0],
          nohfRows[0][1], nohfRows[1][1], nohfRows[2][1]);
        for (auto& p : nohfPerRun) P("NOHFRUN %d %.0f\n", p.first, p.second);
        fprintf(tp, "# rows of the no_hf runs: FAT e+e- %.0f, e+e+ %.0f, e-e- %.0f; pre e+e- %.0f, e+e+ %.0f, e-e- %.0f\n",
                nohfRows[0][0], nohfRows[1][0], nohfRows[2][0], nohfRows[0][1], nohfRows[1][1], nohfRows[2][1]);
    }
    fclose(tp);

    // LaTeX of the peak per period (pre and peak nominal + vertex)
    {
        FILE* o = fopen(Form("%s/smp_peak_periods.tex", smp::kTabOut), "w");
        fprintf(o, "%% smp_peak_periods.tex -- written by plots/exp_report/smp_flow.C (numbers: tables_exp/smp_peak_periods.txt)\n");
        fprintf(o, "\\begin{table}[H]\n\\centering\\small\n\\setlength{\\tabcolsep}{4pt}\n");
        fprintf(o, "\\caption{$\\eta$-peak content per period. $P_\\mathrm{main}$: $\\chi^2$ fit of $N-\\mathrm{CB}$ "
                   "(CB $=2\\sqrt{N_{++}N_{--}}$, $k=1$) with the frozen Gaussian (0.5439~GeV/$c^2$, $\\sigma=5.4$~MeV) and a "
                   "second-order polynomial in 0.5225--0.6100~GeV/$c^2$; $S$ = whole Gaussian $\\pm$ stat $\\pm$ model spread (old envelope). "
                   "All stored combinations (isBest $\\geq 0$), corrected momenta, PT3, $v_z>-160$~mm. "
                   "$s=Y_S=S\\times\\frac{64}{63}/L$ in visible rows per nb$^{-1}$, $L$ = START-free luminosity of the analysed runs; "
                   "visible yield, not corrected for acceptance or efficiency. Constancy over I, IIA, IIB, III with the statistical errors.}\n");
        fprintf(o, "\\label{tab:smp_peak_periods}\n\\begin{tabular}{lr rr rr}\n\\toprule\n");
        fprintf(o, " & & \\multicolumn{2}{c}{pre} & \\multicolumn{2}{c}{peak nominal + vertex} \\\\\n");
        fprintf(o, "Period (days) & $L$ [nb$^{-1}$] & $S$ & $s$ [rows/nb$^{-1}$] & $S$ & $s$ [rows/nb$^{-1}$] \\\\\n\\midrule\n");
        const int order[7] = {0, 1, 2, 3, 4, 7, 6};
        for (int q = 0; q < 7; ++q) {
            const int g = order[q];
            if (q == 4 || q == 6) fprintf(o, "\\midrule\n");
            fprintf(o, "%s & %.1f", kGLong[g], Lg[g]);
            for (int b = 0; b < 2; ++b)
                fprintf(o, " & $%.0f\\pm%.0f\\pm%.0f$ & $%.3f\\pm%.3f$", pf[b][g].S, pf[b][g].eS, spr[b][g].half, pl[b][g].s, pl[b][g].es);
            fprintf(o, " \\\\\n");
        }
        fprintf(o, "\\midrule\n$\\chi^2/\\mathrm{ndf}$ (I, IIA, IIB, III) & ");
        for (int b = 0; b < 2; ++b) {
            std::vector<double> x, e;
            for (int g = 0; g < 4; ++g) { x.push_back(pl[b][g].s); e.push_back(pl[b][g].es); }
            const Const c = Constancy(x, e);
            fprintf(o, " & \\multicolumn{2}{c}{%.1f/%d ($p=%.2f$)}", c.chi2, c.ndf, c.prob);
        }
        fprintf(o, " \\\\\n\\bottomrule\n\\end{tabular}\n\\end{table}\n");
        fclose(o);
    }

    // ================= FIGURES =================
    // ---- 1. rates per day ----
    {
        TCanvas* c = new TCanvas("cR", "", 900, 1000);
        const double yl = 0.68, yh = 1.40;
        const double fr[4] = {1.0, 0.685, 0.385, 0.0};
        const char* lab[3] = {"e^{+}e^{-} counts", "e^{+}e^{+} counts", "e^{-}e^{-} counts"};
        for (int k = 0; k < 3; ++k) {
            c->cd();
            TPad* p = new TPad(Form("pR%d", k), "", 0, fr[k + 1], 1, fr[k]);
            p->SetLeftMargin(0.12); p->SetRightMargin(0.03);
            p->SetTopMargin(k == 0 ? 0.05 : 0.02); p->SetBottomMargin(k == 2 ? 0.22 : 0.02);
            p->Draw(); p->cd();
            // text sizes in pad units: equal absolute size in all three pads
            const double h = fr[k] - fr[k + 1], ts = 0.024 / h;
            TH1D* fr0 = new TH1D(Form("frR%d", k), "", 34, 33.5, 67.5);
            fr0->SetMinimum(yl); fr0->SetMaximum(yh);
            fr0->GetYaxis()->SetTitle("R / #LTR#GT");
            fr0->GetXaxis()->SetTitle("day of 2022");
            fr0->GetYaxis()->SetNdivisions(507);
            fr0->GetYaxis()->SetLabelSize(ts); fr0->GetYaxis()->SetTitleSize(1.1 * ts);
            fr0->GetYaxis()->SetTitleOffset(0.50);
            fr0->GetXaxis()->SetLabelSize(k == 2 ? ts : 0); fr0->GetXaxis()->SetTitleSize(k == 2 ? 1.1 * ts : 0);
            fr0->GetXaxis()->SetTitleOffset(1.0);
            fr0->Draw("axis");
            PeriodLines(p, yl, yh, yl + 0.035 * (yh - yl), k == 2 ? 1.1 * ts : 0.);
            TLine* one = new TLine(33.5, 1, 67.5, 1); one->SetLineStyle(3); one->Draw();
            TGraphErrors* g0 = new TGraphErrors(), *g1 = new TGraphErrors();
            for (size_t i = 0; i < days.size(); ++i) {
                g0->SetPoint(i, days[i] - 0.15, rs[k].r[i] / cR[k].mean); g0->SetPointError(i, 0, rs[k].e[i] / cR[k].mean);
                g1->SetPoint(i, days[i] + 0.15, rs[k].rb[i] / cRb[k].mean); g1->SetPointError(i, 0, rs[k].eb[i] / cRb[k].mean);
            }
            const int col = k == 0 ? kBlack : kRed + 1;
            g0->SetMarkerStyle(20); g0->SetMarkerSize(1.0); g0->SetMarkerColor(col); g0->SetLineColor(col);
            g1->SetMarkerStyle(24); g1->SetMarkerSize(1.0); g1->SetMarkerColor(col); g1->SetLineColor(col);
            g0->Draw("P same"); g1->Draw("P same");
            fr0->Draw("axis same");
            // label above the data (top of the frame)
            TLatex t; t.SetTextSize(ts); t.SetTextAlign(13);
            t.DrawLatex(34.0, yh - 0.03 * (yh - yl),
                        Form("%s:  #LTR#GT = %.1f (isBest #geq 0), %.1f (isBest = 1) counts per nb^{-1}", lab[k], cR[k].mean, cRb[k].mean));
            if (k == 0) {
                TLegend* lg = new TLegend(0.70, 0.05, 0.97, 0.27);
                lg->SetTextSize(ts); lg->SetFillStyle(0);
                lg->AddEntry(g0, "isBest #geq 0", "pe"); lg->AddEntry(g1, "isBest = 1", "pe");
                lg->Draw();
            }
        }
        smp::Save(c, std::string(smp::kFigOut) + "/smp_rates_perday");
    }

    // ---- 2. peak content per period (summary) ----
    {
        TCanvas* c = new TCanvas("cP", "", 900, 560);
        c->SetLeftMargin(0.13); c->SetBottomMargin(0.16);
        const int order[7] = {0, 1, 2, 3, 4, 7, 6};
        TH1D* fr0 = new TH1D("frP", "", 7, 0, 7);
        for (int q = 0; q < 7; ++q) fr0->GetXaxis()->SetBinLabel(q + 1, kGName[order[q]]);
        fr0->GetXaxis()->SetLabelSize(0.065);
        double ymax = 0;
        for (int b = 0; b < 2; ++b) for (int q = 0; q < 7; ++q) ymax = std::max(ymax, pl[b][order[q]].s + std::hypot(pl[b][order[q]].es, pl[b][order[q]].ms));
        fr0->SetMinimum(0); fr0->SetMaximum(1.35 * ymax);
        fr0->GetYaxis()->SetTitle("S #times 64/63 / L  [counts per nb^{-1}]");
        fr0->GetYaxis()->SetTitleOffset(1.05);
        fr0->Draw("axis");
        TLine* sep = new TLine(4, 0, 4, 1.35 * ymax); sep->SetLineStyle(2); sep->SetLineColor(kGray + 1); sep->Draw();
        TLine* sep2 = new TLine(6, 0, 6, 1.35 * ymax); sep2->SetLineStyle(2); sep2->SetLineColor(kGray + 1); sep2->Draw();
        TGraphErrors* gs[2]; TGraphErrors* gm[2];
        const int mk[2] = {20, 25}; const double dx[2] = {-0.12, 0.12};
        for (int b = 0; b < 2; ++b) {
            gs[b] = new TGraphErrors(); gm[b] = new TGraphErrors();
            std::vector<double> x, e;
            for (int g = 0; g < 4; ++g) { x.push_back(pl[b][g].s); e.push_back(pl[b][g].es); }
            const Const cc = Constancy(x, e);
            TBox* band = new TBox(0, cc.mean - cc.emean, 4, cc.mean + cc.emean);
            band->SetFillColorAlpha(b == 0 ? kBlue - 9 : kAzure - 9, 0.35); band->SetLineWidth(0); band->Draw();
            TLine* ml = new TLine(0, cc.mean, 4, cc.mean); ml->SetLineColor(b == 0 ? kBlue + 1 : kAzure + 2); ml->SetLineStyle(b == 0 ? 1 : 7); ml->Draw();
            for (int q = 0; q < 7; ++q) {
                const int g = order[q];
                gs[b]->SetPoint(q, q + 0.5 + dx[b], pl[b][g].s); gs[b]->SetPointError(q, 0, pl[b][g].es);
                gm[b]->SetPoint(q, q + 0.5 + dx[b], pl[b][g].s); gm[b]->SetPointError(q, 0.05, std::hypot(pl[b][g].es, pl[b][g].ms));
            }
            const int col = b == 0 ? kBlue + 1 : kAzure + 2;
            gm[b]->SetFillStyle(0); gm[b]->SetLineColor(col); gm[b]->SetLineWidth(1);
            gs[b]->SetMarkerStyle(mk[b]); gs[b]->SetMarkerSize(1.4); gs[b]->SetMarkerColor(col); gs[b]->SetLineColor(col); gs[b]->SetLineWidth(2);
            gm[b]->Draw("5 same"); gs[b]->Draw("P same");
        }
        fr0->Draw("axis same");
        TLegend* lg = new TLegend(0.15, 0.73, 0.75, 0.91);
        lg->SetTextSize(0.042); lg->SetFillStyle(0);
        for (int b = 0; b < 2; ++b) {
            std::vector<double> x, e;
            for (int g = 0; g < 4; ++g) { x.push_back(pl[b][g].s); e.push_back(pl[b][g].es); }
            const Const cc = Constancy(x, e);
            lg->AddEntry(gs[b], Form("%s: #chi^{2}/ndf (I-III) = %.1f/%d", b == 0 ? "pre" : "peak nominal + vertex", cc.chi2, cc.ndf), "pe");
        }
        lg->Draw();
        TLatex t; t.SetNDC(); t.SetTextSize(0.036); t.SetTextColor(kGray + 2);
        t.DrawLatex(0.15, 0.20, "bars: stat;  boxes: stat #oplus model spread");
        smp::Save(c, std::string(smp::kFigOut) + "/smp_peak_periods");
    }

    // ---- 3. N - CB spectra per period (pre and peak nominal + vertex) ----
    for (int b = 0; b < 2; ++b) {
        const int gl[6] = {0, 1, 2, 3, 4, 7};
        TCanvas* c = new TCanvas(Form("cS%d", b), "", 1200, 760);
        c->Divide(3, 2, 0.002, 0.002);
        // common y range over the panels (per nb^-1, 5 MeV bins)
        double ymin = 0, ymax = 0;
        std::vector<TH1D*> hs(6);
        for (int q = 0; q < 6; ++q) {
            const int g = gl[q];
            const Spectra& sp = spP[b * kNG + g];
            TH1D* h = new TH1D(Form("hS%d_%d", b, q), "", sp.nbin / 2, sp.lo, sp.hi);
            for (int i = 0; i < sp.nbin / 2; ++i) {
                const double v = (sp.sig[2 * i] + sp.sig[2 * i + 1]) * smp::kPT3fac / Lg[g];
                const double e = std::sqrt(sp.sigVar[2 * i] + sp.sigVar[2 * i + 1]) * smp::kPT3fac / Lg[g];
                h->SetBinContent(i + 1, v); h->SetBinError(i + 1, e);
                const double x = h->GetBinCenter(i + 1);
                if (x > 0.47 && x < 0.62) { ymax = std::max(ymax, v + e); ymin = std::min(ymin, v - e); }
            }
            hs[q] = h;
        }
        ymax *= 1.30; ymin = std::min(0., ymin * 1.1);
        for (int q = 0; q < 6; ++q) {
            const int g = gl[q];
            TPad* p = (TPad*)c->cd(q + 1);
            p->SetLeftMargin(0.19); p->SetRightMargin(0.03); p->SetBottomMargin(0.15); p->SetTopMargin(0.06);
            TH1D* h = hs[q];
            h->GetXaxis()->SetRangeUser(0.47, 0.62);
            h->SetMinimum(ymin); h->SetMaximum(ymax);
            h->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]");
            h->GetYaxis()->SetTitle("(N - CB) #times 64/63 / L  [per 5 MeV/c^{2} per nb^{-1}]");
            h->GetXaxis()->SetNdivisions(505);
            h->GetXaxis()->SetLabelSize(0.055); h->GetXaxis()->SetTitleSize(0.058);
            h->GetYaxis()->SetLabelSize(0.055); h->GetYaxis()->SetTitleSize(0.050); h->GetYaxis()->SetTitleOffset(1.70);
            h->SetMarkerStyle(20); h->SetMarkerSize(0.8); h->SetMarkerColor(kBlue + 1); h->SetLineColor(kBlue + 1);
            h->Draw("E1");
            TLine* z = new TLine(0.47, 0, 0.62, 0); z->SetLineStyle(3); z->Draw();
            // fitted curves, per 5 MeV, inside the fit range
            const PeakFit& f = pf[b][g];
            TGraph* gt = new TGraph(); TGraph* gb = new TGraph();
            int np = 0;
            for (double m = scan::kFitLo; m <= scan::kFitHi + 1e-9; m += 0.0005) {
                gt->SetPoint(np, m, 2 * f.Total(m) * smp::kPT3fac / Lg[g]);
                gb->SetPoint(np, m, 2 * f.Bkg(m) * smp::kPT3fac / Lg[g]);
                np++;
            }
            gt->SetLineColor(kBlue + 2); gt->SetLineWidth(2); gt->Draw("L same");
            gb->SetLineColor(kGray + 2); gb->SetLineWidth(2); gb->SetLineStyle(2); gb->Draw("L same");
            TLatex t; t.SetNDC(); t.SetTextSize(0.056);
            t.DrawLatex(0.23, 0.86, kGLong[g]);
            t.SetTextSize(0.050);
            t.DrawLatex(0.23, 0.79, Form("L = %.0f nb^{-1}", Lg[g]));
            t.DrawLatex(0.23, 0.73, Form("S = %.0f #pm %.0f #pm %.0f", f.S, f.eS, spr[b][g].half));
            t.DrawLatex(0.23, 0.67, Form("s = %.3f #pm %.3f counts/nb^{-1}", pl[b][g].s, pl[b][g].es));
        }
        smp::Save(c, std::string(smp::kFigOut) + (b == 0 ? "/smp_peak_spectra_pre" : "/smp_peak_spectra_nomvtx"));
    }

    // ---- 4. sector shares per day ----
    {
        TCanvas* c = new TCanvas("cSec", "", 1100, 520);
        c->Divide(2, 1, 0.002, 0.002);
        const int colS[6] = {kGreen + 2, kOrange + 1, kRed + 1, kAzure + 1, kViolet + 1, kGray + 2};
        const int mkS[6] = {21, 22, 20, 23, 33, 34};
        double ylo = 1, yhi = 0;
        for (int t = 0; t < 2; ++t) for (int q = 0; q < 6; ++q) for (double v : shDay[t][q]) { ylo = std::min(ylo, v); yhi = std::max(yhi, v); }
        ylo = std::floor(ylo * 100 - 1.5) / 100; yhi = std::ceil(yhi * 100 + 2.5) / 100;
        for (int t = 0; t < 2; ++t) {
            TPad* p = (TPad*)c->cd(t + 1);
            p->SetLeftMargin(0.15); p->SetRightMargin(0.03);
            TH1D* f0 = new TH1D(Form("frSec%d", t), "", 34, 33.5, 67.5);
            f0->SetMinimum(ylo); f0->SetMaximum(yhi);
            f0->GetXaxis()->SetTitle("day of 2022");
            f0->GetYaxis()->SetTitle(t == 0 ? "share of lepton tracks per sector" : "share of pion tracks per sector");
            f0->Draw("axis");
            PeriodLines(p, ylo, yhi, ylo + 0.03 * (yhi - ylo), 0.050);
            TLegend* lg = new TLegend(0.17, 0.82, 0.96, 0.91);
            lg->SetNColumns(6); lg->SetTextSize(0.046); lg->SetFillStyle(1001); lg->SetFillColor(kWhite);
            for (int q = 0; q < 6; ++q) {
                TGraph* g = new TGraph();
                for (size_t i = 0; i < days.size(); ++i) g->SetPoint(i, days[i], shDay[t][q][i]);
                g->SetMarkerStyle(mkS[q]); g->SetMarkerColor(colS[q]); g->SetLineColor(colS[q]);
                g->SetMarkerSize(q == 2 ? 1.3 : 1.0);
                g->Draw("P same");
                lg->AddEntry(g, Form("S%d", q), "p");
            }
            lg->Draw();
            f0->Draw("axis same");
        }
        smp::Save(c, std::string(smp::kFigOut) + "/smp_sectors");
    }

    P("\n# total run time %.1f s\n", swAll.RealTime());
    fclose(gOut); gOut = nullptr;
}
