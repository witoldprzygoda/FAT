// ============================================================================
// cb_checks.C  --  combinatorial background (CB) of pi+pi-e+e-: data checks and
//                  sensitivity of the eta results to the CB
// ============================================================================
// HADES p+p, T = 4.53 GeV (feb22), channel pi+ pi- e+ e-. EXPERIMENTAL DATA ONLY.
// CB = 2 sqrt(N++ N--) per 10 MeV bin of M4 = M(pi+pi-e+e-), from the like-sign
// hypotheses pi+pi-e+e+ (N++) and pi+pi-e-e- (N--), linearly interpolated to the
// 2.5 MeV bins (plots/ScanEngine.h), k = 1.
//
// What it does
//   A. like-sign spectra N++, N-- and CB against M4, M(e+e-) and OA(e+e-);
//      ratio N--/N++ and CB/N; per selection and M4 region: N, N++, N--,
//      N--/N++, CB, CB/N, (N++ + N--)/CB, N - CB with event-group errors;
//   B. lepton-multiplicity classes of the like-sign events (join of the trees
//      by file, run and event vertex): share of the like-sign formula per class
//      n = n(e+) + n(e-) and the expected isBest == 1 bias of the plain formula
//      (one row per event and tree; equal-probability choice among the lepton
//      hypotheses of the e+e- tree);
//   C. null tests: frozen eta-peak fit on the like-sign spectra (Poisson
//      likelihood, Gaussian 0.5439 GeV/c2 / 5.4 MeV + pol2), at the eta position
//      and at 46 positions 0.40-0.85 GeV/c2; P_main on N - CB at 13 displaced
//      positions 0.575-0.695; edge-step fit (bump) applied to the like-sign
//      spectra in the four M(pi+pi-) slices;
//   D. sensitivity: eta peak S (P_main), bump edge-step sum E and the bump-
//      window N - CB (W) when the CB is scaled (0.90-1.10; quoted: 0.97, 1.03),
//      for CB bins of 2.5, 5, 10 (nominal), 20 MeV, step or cubic instead of
//      linear interpolation, the CB multiplied by the event-mixing k(M4) of
//      plots/exp_report/cb_kmix.C (read from its ROOT file; run cb_kmix.C
//      first, otherwise this variant is skipped), the like-sign SUM shape
//      normalised to the integrated geometric mean, and the plain like-sign sum
//      (extreme); the same for isBest == 1, and with the CB of isBest == 1
//      scaled by 0.82 (correction found in the e+e- study, for illustration);
//   E. figures and tables.
//
// Selections (rows = all stored combinations, isBest >= 0; corrected momenta,
// cor_ok == 1 (no fallback); already in the FAT production: event vertex
// z > -500 mm and PT3 = trigbit 8192)
//   FAT only     = cor_ok && isBest >= 0                    (no vz > -160 cut)
//   pre          = cor_ok && isBest >= 0 && vz > -160 mm
//   +vertex      = + pip_isUsedVertex == 1 && pim_isUsedVertex == 1
//   +OA>2        = + OA(e+e-) > 2 deg   (not a fixed cut: numbers with and without)
//   peak nominal = pre + oa_lab < 55 deg + M(pi+pi-) > 0.30 GeV/c2
//   bump nominal = pre + oa_lab < 40 deg + 0.30 < M(pi+pi-) < 0.41 GeV/c2
//   frozen simulation selections (sim study, section 2 of A2):
//     DIRECT = pre + dz_pi < 40 mm ; PEAK = DIRECT + MM4 < 2.85 GeV/c2 ; BUMP = pre ;
//     X17 = pre + 15.5 <= M(ee) < 24.5 MeV/c2 + both pions in the vertex fit
//     (the X17 M4 window 0.5333-0.5587 is the sim-centred PEAK window)
//   variant isBest == 1 of the main selections (comparison only)
// Windows: data peak window 0.53364-0.55416, data bump window 0.4550-0.5225
// (the frozen sim windows 0.5333-0.5587 / 0.4576-0.5255 are simulation-centred;
// they appear only as extra regions in the tables).
//
// Inputs : /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/scan_v2/scan_compact.root
//          (scan::kCompact) through ../ScanEngine.h and cb_common.h;
//          /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/cb/cb_kmix.root (k_mix, optional)
// Outputs: /home/przygoda/HADES/FAT/FAT/docs/latex/figs_exp/cb_spectra.pdf, cb_share.pdf,
//          cb_null.pdf, cb_sensitivity.pdf (+ .png; PDF via EPS + ps2pdf -dEPSCrop)
//          /home/przygoda/HADES/FAT/FAT/docs/latex/tables_exp/cb_checks.txt (all numbers),
//          cb_share.tex, cb_null.tex, cb_sens.tex
//          /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/cb/cb_checks.root
//
// Run (from /home/przygoda/HADES/FAT/FAT) after cb_kmix.C, about 1 min, 3 GB of memory:
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/aclic_cb", kTRUE)' 'plots/exp_report/cb_checks.C+'
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TSystem.h"

#include "cb_common.h"

namespace cbc {

// ---- M4 regions: index into Spectra::win -------------------------------------------
// 0 peak window, 1 bump window (engine), then the extra windows below
struct Reg { const char* name; double lo, hi; };
const Reg kExtra[] = {
    {"0.350-0.455", 0.350, 0.455},
    {"lowSB 0.5225-0.5336", 0.5225, scan::kPeakLo},
    {"0.5542-0.600", scan::kPeakHi, 0.600},
    {"0.600-0.700", 0.600, 0.700},
    {"0.700-0.950", 0.700, 0.950},
    {"simPEAK 0.5333-0.5587", 0.5333, 0.5587},
    {"simBUMP 0.4576-0.5255", 0.4576, 0.5255},
    {"all 0.35-0.95", 0.350, 0.950},
    {"fit 0.5225-0.6100", scan::kFitLo, scan::kFitHi},
};
const int kNExtra = sizeof(kExtra) / sizeof(kExtra[0]);
// order of the regions in the tables (indices into win)
const int kRegOrder[] = {2, 1, 3, 0, 4, 5, 6, 9};
const char* const kRegName[] = {"0.350-0.455", "bump 0.4550-0.5225", "0.5225-0.5336", "peak 0.5336-0.5542",
                                "0.5542-0.600", "0.600-0.700", "0.700-0.950", "all 0.35-0.95"};
const int kNReg = 8;

scan::FillOpt M4Opt(int reb = scan::kCBreb) {
    scan::FillOpt o; o.reb = reb;
    for (int i = 0; i < kNExtra; ++i) o.windows.push_back(std::make_pair(kExtra[i].lo, kExtra[i].hi));
    return o;
}

// ---- event key for the join of the three trees ---------------------------------------
struct EvKey {
    int file, run; uint32_t z, x, y;
    bool operator==(const EvKey& o) const { return file == o.file && run == o.run && z == o.z && x == o.x && y == o.y; }
};
struct EvKeyHash {
    size_t operator()(const EvKey& k) const {
        uint64_t h = (uint64_t)k.file * 1000003ULL ^ ((uint64_t)k.run << 20);
        h ^= (uint64_t)k.z * 0x9E3779B97F4A7C15ULL; h ^= ((uint64_t)k.x << 32) ^ ((uint64_t)k.y * 0xC2B2AE3D27D4EB4FULL);
        return (size_t)(h ^ (h >> 29));
    }
};
inline uint32_t Bits(float f) { uint32_t u; std::memcpy(&u, &f, 4); return u; }
inline EvKey KeyOf(const Row& r) { EvKey k; k.file = r.file_id; k.run = r.run; k.z = Bits(r.vz); k.x = Bits(r.vx); k.y = Bits(r.vy); return k; }

// ratio histogram a/b per bin, a and b independent (variances va, vb; default Poisson)
TH1D* RatioH(const char* name, const std::vector<double>& a, const std::vector<double>& b, int n, double lo, double hi,
             const std::vector<double>* va = nullptr, const std::vector<double>* vb = nullptr) {
    TH1D* h = new TH1D(name, "", n, lo, hi); h->SetDirectory(nullptr);
    for (int i = 0; i < n; ++i) {
        if (!(a[i] > 0 && b[i] > 0)) continue;
        const double r = a[i] / b[i];
        const double ra = (va ? (*va)[i] : a[i]) / (a[i] * a[i]), rb = (vb ? (*vb)[i] : b[i]) / (b[i] * b[i]);
        h->SetBinContent(i + 1, r); h->SetBinError(i + 1, r * std::sqrt(ra + rb));
    }
    return h;
}

// k of the event mixing (plots/exp_report/cb_kmix.C) per 10 MeV bin; empty if absent
std::vector<double> ReadK(TFile* f, const std::string& name) {
    std::vector<double> k;
    if (!f) return k;
    TH1D* h = (TH1D*)f->Get(name.c_str());
    if (!h) return k;
    for (int j = 1; j <= h->GetNbinsX(); ++j) k.push_back(h->GetBinContent(j));
    return k;
}

} // namespace cbc

// ============================================================================
void cb_checks() {
    using namespace cbc;
    TStopwatch swAll;
    gSystem->mkdir(kWork, kTRUE);
    gSystem->mkdir(kFigOut, kTRUE);
    gSystem->mkdir(kTabOut, kTRUE);
    OutFile() = fopen(Form("%s/cb_checks.txt", kTabOut), "w");
    P("cb_checks.txt -- written by plots/exp_report/cb_checks.C (experimental data only)\n");
    P("input %s\n", scan::kCompact);
    P("rows = stored combinations; cor_ok == 1 (no fallback); FAT: vertex z > -500 mm, PT3 (trigbit 8192).\n");
    P("CB = 2 sqrt(N++ N--) in 10 MeV bins of M4, linearly interpolated to 2.5 MeV (k = 1) unless stated.\n");
    P("Windows: data peak 0.53364-0.55416, data bump 0.4550-0.5225 GeV/c2. Errors of window sums: event groups.\n\n");

    Data d;
    if (!scan::Load(d)) return;
    const std::vector<SD> sels = Selections();
    const int nS = (int)sels.size();
    std::vector<Sel> sv;
    for (const auto& s : sels) sv.push_back(s.sel);

    // ======================================================================
    // A. M4 spectra of all selections, nominal CB
    // ======================================================================
    std::vector<Spectra> sp = scan::FillMany(d, sv, M4Opt());
    P("A. CB SHARE PER SELECTION AND M4 REGION [GeV/c2]: CB/N (N = unlike-sign rows)\n");
    P("%-34s", "selection");
    for (int q = 0; q < kNReg; ++q) P(" | %-11.11s", kRegName[q]);
    P("\n");
    for (int s = 0; s < nS; ++s) {
        P("%-34s", sels[s].name.c_str());
        for (int q = 0; q < kNReg; ++q) {
            const WinSum& w = sp[s].win[kRegOrder[q]];
            P(" | %5.3f%6s", w.cbShare(), "");
        }
        P("\n");
    }
    P("\nA2. DETAILS in the peak and bump windows (data windows). r = N--/N++, AM/GM = (N++ + N--)/(2 sqrt(N++ N--)) summed in the window,\n");
    P("    disp = (sum over event groups of k^2)/N for the unlike-sign rows (over-dispersion of the window count)\n");
    P("%-34s | %-6s | %8s %7s %7s %5s %8s %6s %6s %6s | %9s %6s %5s\n", "selection", "window", "N", "N++", "N--", "r", "CB", "eCB", "CB/N",
      "AM/GM", "N-CB", "err", "disp");
    for (int s = 0; s < nS; ++s)
        for (int w = 0; w < 2; ++w) {
            const WinSum& x = sp[s].win[w];
            const double am = x.n[1] + x.n[2], gm = 2 * std::sqrt(x.n[1] * x.n[2]);
            P("%-34s | %-6s | %8.0f %7.0f %7.0f %5.2f %8.1f %6.1f %6.4f %6.3f | %9.0f %6.0f %5.3f\n", sels[s].name.c_str(), w == 0 ? "peak" : "bump",
              x.n[0], x.n[1], x.n[2], x.n[1] > 0 ? x.n[2] / x.n[1] : 0., x.cb, std::sqrt(x.cbVar), x.cbShare(), gm > 0 ? am / gm : 0., x.sig(),
              x.sigErr(), x.n[0] > 0 ? x.var[0] / x.n[0] : 0.);
        }
    P("\nA3. LIKE-SIGN COUNTS, all M4 0.35-0.95 and whole sample: rows N++, N--, r = N--/N++, rows per event group (pp, mm, us)\n");
    {
        long ng[3] = {0, 0, 0}, nr[3] = {0, 0, 0};
        for (int k = 0; k < 3; ++k) {
            int cur = -1;
            for (const auto& r : d.s[k]) {
                if (!r.cor_ok || !Pre(r)) continue;
                nr[k]++;
                if (r.grp != cur) { ng[k]++; cur = r.grp; }
            }
        }
        P("   pre, any M4: rows us %ld, pp %ld, mm %ld; groups us %ld, pp %ld, mm %ld; rows/group %.3f %.3f %.3f; N--/N++ = %.3f\n",
          nr[0], nr[1], nr[2], ng[0], ng[1], ng[2], (double)nr[0] / ng[0], (double)nr[1] / ng[1], (double)nr[2] / ng[2], (double)nr[2] / nr[1]);
    }
    const int iFat = Find(sels, "fat"), iPre = Find(sels, "pre"), iPreOa = Find(sels, "preoa"), iPreVtx = Find(sels, "prevtx");
    const int iDir = Find(sels, "direct"), iPkSim = Find(sels, "peaksim"), iPk = Find(sels, "peaknom"), iPkV = Find(sels, "peakvtx");
    const int iPkOa = Find(sels, "peakoa"), iBm = Find(sels, "bumpnom"), iBmV = Find(sels, "bumpvtx"), iBmOa = Find(sels, "bumpoa");
    const int iX17 = Find(sels, "x17"), iPreB = Find(sels, "pre_b1"), iPkB = Find(sels, "peaknom_b1"), iPkVB = Find(sels, "peakvtx_b1");
    const int iBmB = Find(sels, "bumpnom_b1"), iBmVB = Find(sels, "bumpvtx_b1"), iPreVB = Find(sels, "prevtx_b1");
    const int iPreOaB = Find(sels, "preoa_b1"), iDirB = Find(sels, "direct_b1");

    // ---- isBest: kept fractions ----
    P("\nA4. isBest == 1 AGAINST isBest >= 0 (same selection): kept fraction of N, N++, N--, CB, and the CB/N of both, per window\n");
    P("%-24s | %-6s | %6s %6s %6s %6s | %6s %6s\n", "selection", "window", "N", "N++", "N--", "CB", "CB/N>=0", "CB/N==1");
    {
        const int pa[][2] = {{iPre, iPreB}, {iPreOa, iPreOaB}, {iPreVtx, iPreVB}, {iDir, iDirB}, {iPk, iPkB}, {iPkV, iPkVB}, {iBm, iBmB}, {iBmV, iBmVB}};
        for (const auto& p : pa)
            for (int w : {2, 1, 0, 5, 6}) {
                const WinSum &a = sp[p[0]].win[w], &b = sp[p[1]].win[w];
                const char* wn = w == 0 ? "peak" : (w == 1 ? "bump" : (w == 2 ? "<0.455" : (w == 5 ? "0.6-0.7" : "0.7-.95")));
                P("%-24s | %-6s | %6.3f %6.3f %6.3f %6.3f | %6.4f %6.4f\n", sels[p[0]].name.c_str(), wn, b.n[0] / a.n[0], b.n[1] / a.n[1], b.n[2] / a.n[2],
                  b.cb / a.cb, a.cbShare(), b.cbShare());
            }
    }

    // ======================================================================
    // B. lepton-multiplicity classes (join of the three trees)
    // ======================================================================
    P("\nB. LEPTON-MULTIPLICITY CLASSES. Event key = (file_id, run, event vertex x, y, z as stored floats); the e+e- tree gives\n");
    P("   n(e+) = distinct e+ tracks, n(e-) = distinct e- tracks of the event; like-sign events without an e+e- row have 0 leptons of the other sign.\n");
    std::unordered_map<EvKey, std::pair<int, int>, EvKeyHash> usMap;
    long nDupKey = 0;
    {
        int cur = -1;
        for (const auto& r : d.s[0]) {
            if (r.grp == cur) continue;
            cur = r.grp;
            const EvKey k = KeyOf(r);
            if (usMap.count(k)) { nDupKey++; continue; }
            usMap[k] = std::make_pair(r.n_lep1, r.n_lep2);
        }
    }
    // class per row of the like-sign trees: 2, 3, 4 (= >= 4)
    std::vector<int> lsCls[3];
    long nMatch[3] = {0, 0, 0}, nGrp[3] = {0, 0, 0}, nAgree[3] = {0, 0, 0};
    for (int k = 1; k <= 2; ++k) {
        lsCls[k].assign(d.s[k].size(), 2);
        int cur = -1, cls = 2;
        for (size_t i = 0; i < d.s[k].size(); ++i) {
            const Row& r = d.s[k][i];
            if (r.grp != cur) {
                cur = r.grp; nGrp[k]++;
                int other = 0;
                auto it = usMap.find(KeyOf(r));
                if (it != usMap.end()) {
                    nMatch[k]++;
                    other = (k == 1) ? it->second.second : it->second.first;
                    const int same = (k == 1) ? it->second.first : it->second.second;
                    if (same == r.n_lep) nAgree[k]++;
                }
                cls = std::min(4, r.n_lep + other);
            }
            lsCls[k][i] = cls;
        }
    }
    lsCls[0].assign(d.s[0].size(), 2);
    for (size_t i = 0; i < d.s[0].size(); ++i) lsCls[0][i] = std::min(4, d.s[0][i].n_lep1 + d.s[0][i].n_lep2);
    P("   e+e- groups %zu (duplicate keys %ld); e+e+ groups %ld, matched to an e+e- event %ld (%.3f), same n(e+) in both trees %ld;\n",
      usMap.size(), nDupKey, nGrp[1], nMatch[1], (double)nMatch[1] / nGrp[1], nAgree[1]);
    P("   e-e- groups %ld, matched %ld (%.3f), same n(e-) %ld\n", nGrp[2], nMatch[2], (double)nMatch[2] / nGrp[2], nAgree[2]);
    // per class sums for pre and isBest == 1, in regions
    {
        const double rlo[4] = {scan::kBumpLo, scan::kPeakLo, 0.60, 0.35}, rhi[4] = {scan::kBumpHi, scan::kPeakHi, 0.70, 0.95};
        const char* rn[4] = {"bump window", "peak window", "0.60-0.70", "0.35-0.95"};
        double n[2][3][3][4];   // [best][tree][class 2,3,4][region]
        std::memset(n, 0, sizeof(n));
        for (int k = 0; k < 3; ++k)
            for (size_t i = 0; i < d.s[k].size(); ++i) {
                const Row& r = d.s[k][i];
                if (!r.cor_ok || !Pre(r)) continue;
                const int c = lsCls[k][i] - 2;
                for (int q = 0; q < 4; ++q) {
                    if (!(r.m4 >= rlo[q] && r.m4 < rhi[q])) continue;
                    n[0][k][c][q] += 1;
                    if (r.isBest == 1) n[1][k][c][q] += 1;
                }
            }
        P("   selection pre. Per class n = n(e+) + n(e-) (2, 3, >= 4): rows N(us), N++, N--, GM_c = 2 sqrt(N++_c N--_c);\n");
        P("   sum_c GM_c / GM(all) (class sum against the plain formula); equal-probability expectation of the naive isBest == 1 excess:\n");
        P("   X = sum_c GM_c / sum_c eps_c GM_c with eps_2 = 1, eps_3 = 1/2, eps_4 = 1/4 (e+e- hypotheses of the e+e- tree that compete;\n");
        P("   the like-sign tree of such events has one lepton pair, and the pion combinations are assumed to compete equally in both trees).\n");
        for (int q = 0; q < 4; ++q) {
            double gmAll = 2 * std::sqrt((n[0][1][0][q] + n[0][1][1][q] + n[0][1][2][q]) * (n[0][2][0][q] + n[0][2][1][q] + n[0][2][2][q]));
            double gmSum = 0, gmEps = 0, gmc[3];
            const double eps[3] = {1.0, 0.5, 0.25};
            for (int c = 0; c < 3; ++c) { gmc[c] = 2 * std::sqrt(n[0][1][c][q] * n[0][2][c][q]); gmSum += gmc[c]; gmEps += eps[c] * gmc[c]; }
            P("   %-12s | us %8.0f %7.0f %6.0f | ++ %6.0f %5.0f %4.0f | -- %6.0f %5.0f %4.0f | GM_c %7.1f %6.1f %5.1f | share %5.3f %5.3f %5.3f | sum/GM %5.3f | X %5.3f\n",
              rn[q], n[0][0][0][q], n[0][0][1][q], n[0][0][2][q], n[0][1][0][q], n[0][1][1][q], n[0][1][2][q], n[0][2][0][q], n[0][2][1][q], n[0][2][2][q],
              gmc[0], gmc[1], gmc[2], gmc[0] / gmSum, gmc[1] / gmSum, gmc[2] / gmSum, gmSum / gmAll, gmSum / gmEps);
        }
        P("   isBest == 1 per class: rows kept (us), naive formula on isBest == 1 rows GM_c(best) = 2 sqrt(N++_c,best N--_c,best),\n");
        P("   kept fractions eps_us = N_c,best/N_c and eps_ls = GM_c(best)/GM_c, and N_c,best - GM_c(best) (signal + true CB minus naive CB).\n");
        P("   Equal-probability expectation for CB rows: eps_us(CB) = eps_ls x 1 (class 2), x 1/2 (class 3), x 1/4 (class >= 4).\n");
        for (int q = 0; q < 4; ++q) {
            P("   %-12s |", rn[q]);
            for (int c = 0; c < 3; ++c) {
                const double gb = 2 * std::sqrt(n[1][1][c][q] * n[1][2][c][q]), ga = 2 * std::sqrt(n[0][1][c][q] * n[0][2][c][q]);
                P(" class %d: N_best %7.0f GM_best %6.0f eps_us %.3f eps_ls %.3f N-GM %7.0f |", c + 2, n[1][0][c][q], gb,
                  n[0][0][c][q] > 0 ? n[1][0][c][q] / n[0][0][c][q] : 0., ga > 0 ? gb / ga : 0., n[1][0][c][q] - gb);
            }
            P("\n");
        }
    }

    // ======================================================================
    // C. null tests
    // ======================================================================
    P("\nC1. FROZEN PEAK FIT ON THE LIKE-SIGN SPECTRA (Poisson likelihood, Gaussian 0.5439 / 5.4 MeV + pol2, 0.5225-0.6100).\n");
    P("    Expected 0. S(sum) = fit of N++ + N--. Also P_main on N - CB (the measurement) for comparison.\n");
    P("%-34s | %7s %5s %6s | %6s %5s | %6s %5s | %6s %5s %5s | pulls 0.40-0.85 (46 positions): ++ mean rms, -- mean rms, sum mean rms | displaced N-CB 0.575-0.695: mean rms\n",
      "selection", "S", "stat", "spread", "S(++)", "err", "S(--)", "err", "S(sum)", "err", "pull");
    const int nLsPos = 46;
    std::vector<std::vector<double>> lsPull(nS, std::vector<double>(nLsPos, 0.)), dispPull(nS, std::vector<double>(13, 0.));
    std::vector<PeakFit> pmain(nS), lsP(nS), lsM(nS), lsS(nS);
    std::vector<scan::Spread> spread(nS);
    for (int s = 0; s < nS; ++s) {
        pmain[s] = scan::FitPeakMain(sp[s]);
        spread[s] = scan::PeakModelSpread(sp[s]);
        lsP[s] = scan::FitPeakLikeSign(sp[s], 1);
        lsM[s] = scan::FitPeakLikeSign(sp[s], 2);
        const Spectra sum = LikeSignAsSignal(sp[s], 3);
        lsS[s] = scan::FitPeakLikeSign(sum, 1);
        double m[3] = {0, 0, 0}, q2[3] = {0, 0, 0};
        for (int i = 0; i < nLsPos; ++i) {
            const double pos = 0.40 + 0.010 * i, sh = pos - scan::kPeakMean;
            PeakCfg c; c.mean = pos; c.lo = scan::kFitLo + sh; c.hi = scan::kFitHi + sh; c.winLo = scan::kPeakLo + sh; c.winHi = scan::kPeakHi + sh;
            const PeakFit a = scan::FitPeakLikeSign(sp[s], 1, c), b = scan::FitPeakLikeSign(sp[s], 2, c), e = scan::FitPeakLikeSign(sum, 1, c);
            const double pa = a.eS > 0 ? a.S / a.eS : 0, pb = b.eS > 0 ? b.S / b.eS : 0, pe = e.eS > 0 ? e.S / e.eS : 0;
            m[0] += pa; q2[0] += pa * pa; m[1] += pb; q2[1] += pb * pb; m[2] += pe; q2[2] += pe * pe;
            lsPull[s][i] = pe;
        }
        double dm = 0, dq = 0;
        for (int i = 0; i < 13; ++i) {
            const double pos = 0.575 + 0.010 * i, sh = pos - scan::kPeakMean;
            PeakCfg c; c.mean = pos; c.lo = scan::kFitLo + sh; c.hi = scan::kFitHi + sh; c.winLo = scan::kPeakLo + sh; c.winHi = scan::kPeakHi + sh;
            const PeakFit f = scan::FitPeakMain(sp[s], c);
            const double p = f.eS > 0 ? f.S / f.eS : 0;
            dispPull[s][i] = p; dm += p; dq += p * p;
        }
        auto ms = [](double a, double b, int n, double& mean, double& rms) { mean = a / n; rms = std::sqrt(b / n); };
        double mm[3], rr[3], dmean, drms;
        for (int k = 0; k < 3; ++k) ms(m[k], q2[k], nLsPos, mm[k], rr[k]);
        ms(dm, dq, 13, dmean, drms);
        P("%-34s | %7.0f %5.0f %6.0f | %6.0f %5.0f | %6.0f %5.0f | %6.0f %5.0f %5.2f | %5.2f %4.2f, %5.2f %4.2f, %5.2f %4.2f | %5.2f %4.2f\n",
          sels[s].name.c_str(), pmain[s].S, pmain[s].eS, spread[s].half, lsP[s].S, lsP[s].eS, lsM[s].S, lsM[s].eS, lsS[s].S, lsS[s].eS,
          lsS[s].eS > 0 ? lsS[s].S / lsS[s].eS : 0., mm[0], rr[0], mm[1], rr[1], mm[2], rr[2], dmean, drms);
    }
    P("    pulls = S/sigma_S with Poisson errors; neighbouring positions share most of their fit range (87.5 MeV wide, 10 MeV steps),\n");
    P("    so the pulls are correlated; event groups with several rows make the Poisson errors slightly too small (see disp in A2).\n");

    // ---- bump: edge steps on N - CB and on the like-sign spectra ----
    P("\nC2. BUMP EDGE-STEP FIT (fixed edge, pol2, slices of M(pi+pi-)) ON N - CB AND ON THE LIKE-SIGN SPECTRA (null test: h = 0 expected)\n");
    P("    h in rows per 2.5 MeV; E = sum over the four slices 0.30-0.34, 0.34-0.37, 0.37-0.39, 0.39-0.40.\n");
    const char* bumpKeys[] = {"fat", "pre", "preoa", "prevtx", "bumpnom", "bumpvtx", "bumpoa", "peakvtx", "pre_b1", "bumpnom_b1", "bumpvtx_b1"};
    const int nB = sizeof(bumpKeys) / sizeof(bumpKeys[0]);
    std::vector<int> bIdx;
    for (int b = 0; b < nB; ++b) bIdx.push_back(Find(sels, bumpKeys[b]));
    std::vector<std::vector<Spectra>> slSp(nB);   // [selection][slice], nominal binning
    {
        std::vector<Sel> ss;
        for (int b = 0; b < nB; ++b)
            for (int i = 0; i < kNSig; ++i) {
                const Sel base = sels[bIdx[b]].sel; const float lo = kSlLo[i], hi = kSlHi[i];
                ss.push_back([=](const Row& r) { return r.mpipi >= lo && r.mpipi < hi && base(r); });
            }
        const std::vector<Spectra> all = scan::FillMany(d, ss, M4Opt());
        for (int b = 0; b < nB; ++b) for (int i = 0; i < kNSig; ++i) slSp[b].push_back(all[b * kNSig + i]);
    }
    std::vector<EVal> eNom(nB), eLs(nB), eLsP(nB), eLsM(nB);
    P("%-30s | %7s %5s | %-35s | %-45s | %7s %5s %5s | ++: E err | --: E err\n", "selection", "E", "stat", "h per slice (N - CB)", "h per slice (N++ + N--)", "E(ls)", "err", "pull");
    for (int b = 0; b < nB; ++b) {
        eNom[b] = EdgeSum(slSp[b]);
        std::vector<Spectra> ls, lp, lm;
        for (int i = 0; i < kNSig; ++i) { ls.push_back(LikeSignAsSignal(slSp[b][i], 3)); lp.push_back(LikeSignAsSignal(slSp[b][i], 1)); lm.push_back(LikeSignAsSignal(slSp[b][i], 2)); }
        eLs[b] = EdgeSum(ls); eLsP[b] = EdgeSum(lp); eLsM[b] = EdgeSum(lm);
        P("%-30s | %7.1f %5.1f | %6.1f %6.1f %6.1f %6.1f%s | %6.1f+-%4.1f %6.1f+-%4.1f %6.1f+-%4.1f %6.1f+-%4.1f | %7.1f %5.1f %5.2f | %6.1f %5.1f | %6.1f %5.1f\n",
          sels[bIdx[b]].name.c_str(), eNom[b].E, eNom[b].eE, eNom[b].h[0], eNom[b].h[1], eNom[b].h[2], eNom[b].h[3], eNom[b].ok ? "  " : " !",
          eLs[b].h[0], eLs[b].eh[0], eLs[b].h[1], eLs[b].eh[1], eLs[b].h[2], eLs[b].eh[2], eLs[b].h[3], eLs[b].eh[3], eLs[b].E, eLs[b].eE,
          eLs[b].eE > 0 ? eLs[b].E / eLs[b].eE : 0., eLsP[b].E, eLsP[b].eE, eLsM[b].E, eLsM[b].eE);
    }
    P("    The like-sign E is in like-sign rows; the CB is 2 sqrt(N++ N--) = about %.2f of N++ + N-- (pre, bump window), so a step in the\n",
      2 * std::sqrt(sp[iPre].win[1].n[1] * sp[iPre].win[1].n[2]) / (sp[iPre].win[1].n[1] + sp[iPre].win[1].n[2]));
    P("    CB would be that fraction of E(ls).\n");

    // ======================================================================
    // D. sensitivity
    // ======================================================================
    P("\nD1. SENSITIVITY OF THE ETA PEAK S (P_main) AND OF THE BUMP-WINDOW W = N - CB TO THE CB. dS = S(variant) - S(nominal).\n");
    P("    variants: CB x 0.97, x 1.03; CB bins 2.5, 5, 20 MeV (nominal 10, refilled); step and cubic interpolation (10 MeV bins);\n");
    P("    'kmix' = CB x k(M4) of the event mixing (cb_kmix.C, variant A, the selection's own k; 'n/a' where none);\n");
    P("    'LS-sum shape' = (N++ + N--) x integrated 2 sqrt(N++ N--)/(N++ + N--); 'AM' = N++ + N-- (extreme, not an estimator).\n");
    TFile* fk = TFile::Open(kKmixFile, "READ");
    if (!fk || fk->IsZombie()) { P("    (no %s: run cb_kmix.C first; the kmix variant is skipped)\n", kKmixFile); fk = nullptr; }
    auto kOf = [&](const std::string& key) { return ReadK(fk, "k_A_" + key); };
    const char* sensKeys[] = {"fat", "pre", "preoa", "prevtx", "direct", "peaksim", "peaknom", "peakvtx", "peakoa", "bumpnom", "bumpvtx", "x17",
                              "pre_b1", "preoa_b1", "prevtx_b1", "peaknom_b1", "peakvtx_b1"};
    const int nSens = sizeof(sensKeys) / sizeof(sensKeys[0]);
    std::vector<int> sIdx;
    for (int k = 0; k < nSens; ++k) sIdx.push_back(Find(sels, sensKeys[k]));
    // refilled binning variants
    const int rebs[3] = {1, 2, 8};
    std::vector<std::vector<Spectra>> spReb(3);
    for (int q = 0; q < 3; ++q) {
        std::vector<Sel> ss;
        for (int k = 0; k < nSens; ++k) ss.push_back(sels[sIdx[k]].sel);
        spReb[q] = scan::FillMany(d, ss, M4Opt(rebs[q]));
    }
    const int nVar = 10;
    const char* varName[nVar] = {"x0.97", "x1.03", "2.5MeV", "5MeV", "20MeV", "step", "cubic", "kmix", "LSsum", "AM"};
    const int iKmix = 7;
    std::vector<std::vector<double>> dSv(nSens, std::vector<double>(nVar, 0.)), dWv(nSens, std::vector<double>(nVar, 0.));
    std::vector<bool> hasK(nSens, false);
    P("%-34s | %7s %5s %6s %6s | %s | W = N-CB bump: nominal, then dW/W for the same variants\n", "selection", "S", "stat", "spread", "CB/N",
      "dS for x0.97 x1.03 2.5MeV 5MeV 20MeV step cubic kmix LSsum AM");
    for (int k = 0; k < nSens; ++k) {
        const Spectra& s0 = sp[sIdx[k]];
        const PeakFit f0 = pmain[sIdx[k]];
        const double w0 = s0.win[1].sig();
        std::vector<Spectra> v;
        v.push_back(Modify(s0, 0.97)); v.push_back(Modify(s0, 1.03));
        for (int q = 0; q < 3; ++q) v.push_back(spReb[q][k]);
        v.push_back(Modify(s0, 1.0, 1)); v.push_back(Modify(s0, 1.0, 2));
        const std::vector<double> kk = sels[sIdx[k]].best ? std::vector<double>() : kOf(sels[sIdx[k]].key);
        hasK[k] = ((int)kk.size() == s0.nw);
        v.push_back(hasK[k] ? Modify(s0, 1.0, 0, &kk) : s0);
        v.push_back(Modify(s0, 1.0, 3)); v.push_back(Modify(s0, 1.0, 4));
        P("%-34s | %7.0f %5.0f %6.0f %6.4f |", sels[sIdx[k]].name.c_str(), f0.S, f0.eS, spread[sIdx[k]].half, f0.N > 0 ? f0.CB / f0.N : 0.);
        for (int q = 0; q < nVar; ++q) {
            const PeakFit f = scan::FitPeakMain(v[q]);
            dSv[k][q] = f.ok ? f.S - f0.S : 0;
            dWv[k][q] = w0 != 0 ? (v[q].win[1].sig() - w0) / w0 : 0;
            if (q == iKmix && !hasK[k]) { P(" %6s", "n/a"); continue; }
            P(" %6.0f", dSv[k][q]);
        }
        P(" | %8.0f", w0);
        for (int q = 0; q < nVar; ++q) { if (q == iKmix && !hasK[k]) P(" %6s", "n/a"); else P(" %6.3f", dWv[k][q]); }
        P("\n");
    }
    // relative: dS/S
    P("\n    the same as relative changes dS/S:\n");
    for (int k = 0; k < nSens; ++k) {
        P("%-34s |", sels[sIdx[k]].name.c_str());
        for (int q = 0; q < nVar; ++q) { if (q == iKmix && !hasK[k]) P(" %7s", "n/a"); else P(" %7.4f", pmain[sIdx[k]].S != 0 ? dSv[k][q] / pmain[sIdx[k]].S : 0.); }
        P("   (stat/S %.3f, spread/S %.3f)\n", pmain[sIdx[k]].eS / pmain[sIdx[k]].S, spread[sIdx[k]].half / pmain[sIdx[k]].S);
    }

    // ---- isBest == 1 with the CB scaled by 0.82 ----
    P("\nD2. isBest == 1: the plain formula on isBest == 1 rows over-estimates the CB (e+e- study: x 1.21 simulation, x 1.22 data emulation).\n");
    P("    This channel (section B, equal-probability choice among the competing lepton hypotheses, class shares of the data):\n");
    P("    naive / class-weighted = 1.14-1.16, i.e. a factor 0.86-0.88. Illustration with both factors, 0.82 and 0.87: S, W and CB/N:\n");
    for (double fsc : {0.82, 0.87}) {
        P("    factor %.2f:\n", fsc);
        const int pb[] = {iPreB, iPreOaB, iPreVB, iPkB, iPkVB, iBmB, iBmVB};
        for (int k : pb) {
            const Spectra m = Modify(sp[k], fsc);
            const PeakFit a = pmain[k], b = scan::FitPeakMain(m);
            P("%-40s | S %6.0f +- %4.0f -> %6.0f (%+.3f) | W %7.0f -> %7.0f (%+.3f) | CB/N peak %.4f -> %.4f, bump %.4f -> %.4f\n", sels[k].name.c_str(),
              a.S, a.eS, b.S, (b.S - a.S) / a.S, sp[k].win[1].sig(), m.win[1].sig(), (m.win[1].sig() - sp[k].win[1].sig()) / sp[k].win[1].sig(),
              sp[k].win[0].cbShare(), m.win[0].cbShare(), sp[k].win[1].cbShare(), m.win[1].cbShare());
        }
    }
    P("    S(isBest==1)/S(isBest>=0) at nominal CB: pre %.3f, peak nominal %.3f, peak nominal + vertex %.3f\n", pmain[iPreB].S / pmain[iPre].S,
      pmain[iPkB].S / pmain[iPk].S, pmain[iPkVB].S / pmain[iPkV].S);

    // ---- bump E sensitivity ----
    P("\nD3. SENSITIVITY OF THE BUMP EDGE-STEP SUM E (rows per 2.5 MeV) TO THE CB, same variants (binning variants refilled).\n");
    std::vector<std::vector<std::vector<Spectra>>> slReb(3, std::vector<std::vector<Spectra>>(nB));
    for (int q = 0; q < 3; ++q) {
        std::vector<Sel> ss;
        for (int b = 0; b < nB; ++b)
            for (int i = 0; i < kNSig; ++i) {
                const Sel base = sels[bIdx[b]].sel; const float lo = kSlLo[i], hi = kSlHi[i];
                ss.push_back([=](const Row& r) { return r.mpipi >= lo && r.mpipi < hi && base(r); });
            }
        const std::vector<Spectra> all = scan::FillMany(d, ss, M4Opt(rebs[q]));
        for (int b = 0; b < nB; ++b) for (int i = 0; i < kNSig; ++i) slReb[q][b].push_back(all[b * kNSig + i]);
    }
    std::vector<std::vector<double>> dEv(nB, std::vector<double>(nVar, 0.));
    P("%-30s | %7s %5s | dE for x0.97 x1.03 2.5MeV 5MeV 20MeV step cubic kmix LSsum AM | dE/E\n", "selection", "E", "stat");
    P("    (kmix for E: k of the selection and M(pi+pi-) slice where cb_kmix.C made one, else the k of the selection)\n");
    std::vector<bool> hasKb(nB, false);
    for (int b = 0; b < nB; ++b) {
        std::vector<EVal> ev;
        ev.push_back(EdgeSum(ModifyAll(slSp[b], 0.97, 0))); ev.push_back(EdgeSum(ModifyAll(slSp[b], 1.03, 0)));
        for (int q = 0; q < 3; ++q) ev.push_back(EdgeSum(slReb[q][b]));
        ev.push_back(EdgeSum(ModifyAll(slSp[b], 1.0, 1))); ev.push_back(EdgeSum(ModifyAll(slSp[b], 1.0, 2)));
        {
            std::vector<std::vector<double>> ks(kNSig);
            bool ok = !sels[bIdx[b]].best;
            for (int i = 0; i < kNSig && ok; ++i) {
                ks[i] = ReadK(fk, Form("k_A_%s_sl%d", sels[bIdx[b]].key.c_str(), i));
                if ((int)ks[i].size() != slSp[b][i].nw) ks[i] = kOf(sels[bIdx[b]].key);
                ok = ((int)ks[i].size() == slSp[b][i].nw);
            }
            hasKb[b] = ok;
            ev.push_back(ok ? EdgeSum(ModifyAll(slSp[b], 1.0, 0, &ks)) : eNom[b]);
        }
        ev.push_back(EdgeSum(ModifyAll(slSp[b], 1.0, 3))); ev.push_back(EdgeSum(ModifyAll(slSp[b], 1.0, 4)));
        P("%-30s | %7.1f %5.1f |", sels[bIdx[b]].name.c_str(), eNom[b].E, eNom[b].eE);
        for (int q = 0; q < nVar; ++q) { dEv[b][q] = ev[q].E - eNom[b].E; if (q == iKmix && !hasKb[b]) P(" %6s", "n/a"); else P(" %6.1f", dEv[b][q]); }
        P(" |");
        for (int q = 0; q < nVar; ++q) { if (q == iKmix && !hasKb[b]) P(" %6s", "n/a"); else P(" %6.3f", eNom[b].E != 0 ? dEv[b][q] / eNom[b].E : 0.); }
        P("\n");
    }

    // ---- scale curves for the figure ----
    const int nF = 21;
    const char* curveKeys[] = {"pre", "preoa", "peaknom", "peakvtx", "bumpnom", "bumpvtx", "pre_b1"};
    const int nC = sizeof(curveKeys) / sizeof(curveKeys[0]);
    std::vector<TGraph*> gS(nC), gW(nC), gE(nC);
    P("\nD4. S(f)/S(1), W(f)/W(1), E(f)/E(1) for the CB scaled by f (numbers only):\n");
    for (int c = 0; c < nC; ++c) {
        const int k = Find(sels, curveKeys[c]);
        int b = -1;
        for (int q = 0; q < nB; ++q) if (std::string(bumpKeys[q]) == curveKeys[c]) b = q;
        gS[c] = new TGraph(); gW[c] = new TGraph(); gE[c] = new TGraph();
        P("%-30s |", sels[k].name.c_str());
        for (int i = 0; i < nF; ++i) {
            const double f = 0.90 + 0.01 * i;
            const Spectra m = Modify(sp[k], f);
            const PeakFit p = scan::FitPeakMain(m);
            gS[c]->SetPoint(i, f, p.S / pmain[k].S);
            gW[c]->SetPoint(i, f, m.win[1].sig() / sp[k].win[1].sig());
            if (b >= 0) gE[c]->SetPoint(gE[c]->GetN(), f, EdgeSum(ModifyAll(slSp[b], f, 0)).E / eNom[b].E);
            if (i % 5 == 0) P(" f=%.2f: %.4f %.4f%s", f, p.S / pmain[k].S, m.win[1].sig() / sp[k].win[1].sig(),
                              b >= 0 ? Form(" %.4f", EdgeSum(ModifyAll(slSp[b], f, 0)).E / eNom[b].E) : "");
        }
        P("\n");
    }

    // ======================================================================
    // E. spectra in M(ee) and OA (selection pre) for the figure
    // ======================================================================
    scan::FillOpt oMee; oMee.nbin = 100; oMee.lo = 0.0; oMee.hi = 0.50; oMee.reb = 1; oMee.x = [](const Row& r) { return (double)r.mee; };
    scan::FillOpt oOa;  oOa.nbin = 60;  oOa.lo = 0.0;  oOa.hi = 120.; oOa.reb = 1; oOa.x = [](const Row& r) { return (double)r.oa; };
    scan::FillOpt oM4;  oM4.nbin = 60;  oM4.lo = 0.35; oM4.hi = 0.95; oM4.reb = 1;
    std::vector<Sel> s1(1, sels[iPre].sel);
    const Spectra spMee = scan::FillMany(d, s1, oMee)[0];
    const Spectra spOa  = scan::FillMany(d, s1, oOa)[0];
    const Spectra spM4w = scan::FillMany(d, s1, oM4)[0];   // 10 MeV bins, CB per bin (no interpolation)
    P("\nE. pre: N--/N++ and CB/N in M(ee) and OA classes\n");
    {
        const double eLo[] = {0.0, 0.005, 0.010, 0.020, 0.050, 0.100, 0.140, 0.200, 0.300};
        const double eHi[] = {0.005, 0.010, 0.020, 0.050, 0.100, 0.140, 0.200, 0.300, 0.500};
        for (int q = 0; q < 9; ++q) {
            const WinSum w = spMee.WindowApprox(eLo[q], eHi[q]);
            P("   M(ee) %5.3f-%5.3f: N %8.0f N++ %6.0f N-- %6.0f r %5.2f CB/N %6.4f\n", eLo[q], eHi[q], w.n[0], w.n[1], w.n[2], w.n[1] > 0 ? w.n[2] / w.n[1] : 0., w.cbShare());
        }
        const double oLo[] = {0, 2, 4, 6, 10, 15, 20, 30, 50, 80};
        const double oHi[] = {2, 4, 6, 10, 15, 20, 30, 50, 80, 120};
        for (int q = 0; q < 10; ++q) {
            const WinSum w = spOa.WindowApprox(oLo[q], oHi[q]);
            P("   OA %3.0f-%3.0f deg: N %8.0f N++ %6.0f N-- %6.0f r %5.2f CB/N %6.4f\n", oLo[q], oHi[q], w.n[0], w.n[1], w.n[2], w.n[1] > 0 ? w.n[2] / w.n[1] : 0., w.cbShare());
        }
    }

    // ======================================================================
    // figures
    // ======================================================================
    Style();
    TFile fr(Form("%s/cb_checks.root", kWork), "RECREATE");
    // ---- Figure 1: like-sign spectra ----
    {
        TCanvas* c = new TCanvas("cb_spectra", "", 1500, 1250);
        const Spectra* S3[3] = {&spM4w, &spMee, &spOa};
        const char* xt[3] = {"M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]", "M(e^{+}e^{-}) [GeV/c^{2}]", "OA(e^{+}e^{-}) [deg]"};
        const char* yt[3] = {"rows / 10 MeV/c^{2}", "rows / 5 MeV/c^{2}", "rows / 2 deg"};
        const double hgt[3] = {0.46, 0.27, 0.27};
        double ytop = 1.0;
        for (int row = 0; row < 3; ++row) {
            const double ybot = ytop - hgt[row];
            for (int col = 0; col < 3; ++col) {
                c->cd();
                TPad* p = new TPad(Form("p%d%d", row, col), "", col / 3.0, ybot, (col + 1) / 3.0, ytop);
                p->SetLeftMargin(0.17); p->SetRightMargin(0.03);
                p->SetTopMargin(row == 0 ? 0.08 : 0.03); p->SetBottomMargin(row == 2 ? 0.25 : 0.05);
                p->Draw(); p->cd();
                const Spectra& s = *S3[col];
                const double tsc = 0.46 / hgt[row];
                if (row == 0) {
                    p->SetLogy();
                    TH1D* hN = s.Hist("us", Form("hN_%d", col)); TH1D* hP = s.Hist("pp", Form("hP_%d", col));
                    TH1D* hM = s.Hist("mm", Form("hM_%d", col)); TH1D* hC = s.Hist("cb", Form("hC_%d", col));
                    Deco(hN, kBlack, 20); Deco(hP, kGreen + 2, 25); Deco(hM, kOrange + 7, 26); Deco(hC, kRed, 20, 0.6);
                    double mx = hN->GetMaximum(), mn = 1e9;
                    for (int i = 1; i <= hP->GetNbinsX(); ++i) if (hP->GetBinContent(i) > 0) mn = std::min(mn, hP->GetBinContent(i));
                    hN->SetMaximum(mx * 60); hN->SetMinimum(std::max(0.5, mn * 0.5));
                    hN->GetYaxis()->SetTitle(yt[col]); hN->GetXaxis()->SetLabelSize(0);
                    hN->GetYaxis()->SetLabelSize(0.050); hN->GetYaxis()->SetTitleSize(0.055); hN->GetYaxis()->SetTitleOffset(1.45);
                    hN->Draw("E1"); hC->Draw("E1 same"); hP->Draw("E1 same"); hM->Draw("E1 same"); hN->Draw("E1 same");
                    hN->Write(); hP->Write(); hM->Write(); hC->Write();
                    if (col == 0) {
                        TLegend* l = new TLegend(0.20, 0.70, 0.95, 0.90); l->SetNColumns(2);
                        l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(0.048);
                        l->AddEntry(hN, "N (e^{+}e^{-})", "pe"); l->AddEntry(hC, "CB = 2#sqrt{N_{++}N_{--}}", "pe");
                        l->AddEntry(hP, "N_{++} (e^{+}e^{+})", "pe"); l->AddEntry(hM, "N_{--} (e^{-}e^{-})", "pe");
                        l->Draw();
                    }
                    Txt(0.19, 0.935, Form("(%c) data, pre, isBest #geq 0", 'a' + col), 0.050);
                } else {
                    TH1D* h = nullptr;
                    if (row == 1) {
                        h = RatioH(Form("hR_%d", col), s.n[2], s.n[1], s.nbin, s.lo, s.hi);
                        Deco(h, kBlack, 24); h->SetMinimum(0.0); h->SetMaximum(11.9); h->GetYaxis()->SetTitle("N_{--} / N_{++}");
                    } else {
                        std::vector<double> cbv(s.cb), cbvv(s.cbVar);
                        h = RatioH(Form("hS_%d", col), cbv, s.n[0], s.nbin, s.lo, s.hi, &cbvv, nullptr);
                        Deco(h, kRed, 20, 0.6); h->SetMinimum(0.0); h->SetMaximum(0.85); h->GetYaxis()->SetTitle("CB / N");
                        h->GetXaxis()->SetTitle(xt[col]);
                    }
                    h->GetYaxis()->SetLabelSize(0.050 * tsc); h->GetYaxis()->SetTitleSize(0.055 * tsc); h->GetYaxis()->SetTitleOffset(1.45 / tsc);
                    h->GetYaxis()->SetNdivisions(505);
                    h->GetXaxis()->SetLabelSize(row == 2 ? 0.050 * tsc * 0.85 : 0); h->GetXaxis()->SetTitleSize(0.055 * tsc * 0.85);
                    h->GetXaxis()->SetTitleOffset(1.0);
                    h->Draw("E1");
                    h->Write();
                    if (row == 1) { TLine* l = new TLine(s.lo, 1, s.hi, 1); l->SetLineStyle(2); l->SetLineColor(kGray + 1); l->Draw(); }
                }
            }
            ytop = ybot;
        }
        Save(c, "cb_spectra");
        c->Write();
    }
    // ---- Figure 2: CB/N against M4 for several selections ----
    {
        TCanvas* c = new TCanvas("cb_share", "", 1400, 560);
        c->Divide(2, 1, 0.002, 0.002);
        const int grpA[] = {iFat, iPre, iPreVtx, iDir, iPkSim};
        const int grpB[] = {iPre, iPreOa, iPk, iBm, iPreB, iPkV};
        const int colA[] = {kGray + 2, kBlack, kBlue + 1, kGreen + 2, kOrange + 7};
        // [review fix] one selection key with spec_likesign: peak nom. + vertex violet, isBest = 1 cyan, bump nominal spring green
        const int colB[] = {kBlack, kGray + 1, kAzure + 7, kSpring + 5, kCyan + 2, kViolet + 1};
        const int mkA[] = {24, 20, 21, 22, 23};
        const int mkB[] = {20, 27, 25, 26, 30, 33};   // pre + OA>2: open diamond (open circle = "FAT only" in panel a)
        for (int pnl = 0; pnl < 2; ++pnl) {
            c->cd(pnl + 1);
            gPad->SetLeftMargin(0.13); gPad->SetRightMargin(0.03);
            TH1D* fr0 = new TH1D(Form("frs%d", pnl), "", 10, 0.35, 0.95);
            fr0->SetMinimum(0.0); fr0->SetMaximum(0.27);
            fr0->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]"); fr0->GetYaxis()->SetTitle("CB / N  (10 MeV/c^{2} bins)");
            fr0->GetYaxis()->SetTitleOffset(1.10);
            fr0->Draw("AXIS");
            TLegend* l = new TLegend(0.15, 0.68, 0.97, 0.90); l->SetNColumns(2);
            l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(0.040);
            const int n = pnl == 0 ? 5 : 6;
            for (int q = 0; q < n; ++q) {
                const int k = pnl == 0 ? grpA[q] : grpB[q];
                const Spectra& s = sp[k];
                TGraphErrors* g = new TGraphErrors();
                for (int j = 0; j < s.nw; ++j) {
                    if (!(s.w[0][j] > 0)) continue;
                    const double x = s.lo + (j + 0.5) * s.reb * s.bw + (q - 2) * 0.0008;
                    const double r = s.cbW[j] / s.w[0][j];
                    const double e = std::sqrt(s.cbWVar[j]) / s.w[0][j];
                    g->SetPoint(g->GetN(), x, r); g->SetPointError(g->GetN() - 1, 0, e);
                }
                Deco(g, pnl == 0 ? colA[q] : colB[q], pnl == 0 ? mkA[q] : mkB[q], 0.8);
                g->SetLineWidth(1);
                g->Draw("P same");
                g->Write(Form("g_share_%s", sels[k].key.c_str()));
                l->AddEntry(g, sels[k].name.c_str(), "pe");
            }
            l->Draw();
            Txt(0.14, 0.935, pnl == 0 ? "(a) data, isBest #geq 0: vertex cuts" : "(b) data: kinematic cuts, OA > 2#circ, isBest = 1 (pre)", 0.045);
            for (double x : {scan::kBumpLo, scan::kBumpHi, scan::kPeakLo, scan::kPeakHi}) {
                TLine* ln = new TLine(x, 0, x, 0.17); ln->SetLineStyle(3); ln->SetLineColor(kGray + 1); ln->Draw();
            }
        }
        Save(c, "cb_share");
        c->Write();
    }
    // ---- Figure 3: null tests ----
    {
        TCanvas* c = new TCanvas("cb_null", "", 1500, 520);
        c->Divide(3, 1, 0.002, 0.002);
        // (a) like-sign sum near the eta with the null fit
        c->cd(1); gPad->SetLeftMargin(0.19);
        {
            const Spectra sum = LikeSignAsSignal(sp[iPre], 3);
            TH1D* h = sum.Hist("pp", "h_lssum");
            h->GetXaxis()->SetRangeUser(0.50, 0.65);
            h->GetXaxis()->SetNdivisions(505);
            Deco(h, kRed, 20, 0.7);
            h->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-}e^{#pm}e^{#pm}) [GeV/c^{2}]"); h->GetYaxis()->SetTitle("like-sign rows / 2.5 MeV/c^{2}");
            h->GetYaxis()->SetTitleOffset(1.65);
            double mx = 0; for (int i = h->FindBin(0.50); i <= h->FindBin(0.65); ++i) mx = std::max(mx, h->GetBinContent(i));
            h->SetMinimum(0); h->SetMaximum(mx * 1.45);
            h->Draw("E1");
            const PeakFit f = lsS[iPre];
            TGraph* gt = new TGraph(); TGraph* gb = new TGraph();
            for (int i = 0; i <= 200; ++i) {
                const double x = scan::kFitLo + (scan::kFitHi - scan::kFitLo) * i / 200.;
                gt->SetPoint(i, x, f.Total(x)); gb->SetPoint(i, x, f.Bkg(x));
            }
            gt->SetLineColor(kBlue + 1); gt->SetLineWidth(2); gb->SetLineColor(kGray + 2); gb->SetLineStyle(2); gb->SetLineWidth(2);
            gb->Draw("L same"); gt->Draw("L same");
            TLegend* l = new TLegend(0.22, 0.68, 0.97, 0.90); l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(0.048);
            l->AddEntry(h, "N_{++} + N_{--}, pre", "pe");
            l->AddEntry(gt, Form("frozen #eta fit: S = %.0f #pm %.0f", f.S, f.eS), "l");
            l->AddEntry(gb, "its pol2 part", "l");
            l->Draw();
            Txt(0.20, 0.935, "(a) like-sign rows near m_{#eta}", 0.054);
        }
        // (b) pulls against position
        c->cd(2); gPad->SetLeftMargin(0.15);
        {
            TH1D* fr0 = new TH1D("frpull", "", 10, 0.38, 0.87);
            fr0->SetMinimum(-5.0); fr0->SetMaximum(9.0);
            fr0->GetXaxis()->SetNdivisions(505);
            fr0->GetXaxis()->SetTitle("position of the test Gaussian [GeV/c^{2}]"); fr0->GetYaxis()->SetTitle("S / #sigma_{S}");
            fr0->Draw("AXIS");
            const int ks[3] = {iPre, iPkV, iBm};
            const int cl[3] = {kRed, kOrange + 7, kMagenta + 1};
            const int mk[3] = {20, 25, 26};
            TLegend* l = new TLegend(0.17, 0.62, 0.97, 0.90); l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(0.046);
            for (int q = 0; q < 3; ++q) {
                TGraph* g = new TGraph();
                for (int i = 0; i < nLsPos; ++i) g->SetPoint(i, 0.40 + 0.01 * i + (q - 1) * 0.002, lsPull[ks[q]][i]);
                Deco(g, cl[q], mk[q], 0.8); g->Draw("P same");
                l->AddEntry(g, Form("N_{++}+N_{--}, %s", sels[ks[q]].name.c_str()), "p");
            }
            TGraph* gd = new TGraph();
            for (int i = 0; i < 13; ++i) gd->SetPoint(i, 0.575 + 0.01 * i, dispPull[iPre][i]);
            Deco(gd, kBlue + 1, 21, 0.9); gd->Draw("P same");
            // (verifier fix 2026-11-14) the measured eta point, S/sigma_S = 12.7 for pre, lies far above the
            // frame (-5..9); it was drawn as an invisible star. Its value is now given in the legend only.
            l->AddEntry(gd, "N - CB, pre, displaced positions", "p");
            l->AddEntry((TObject*)nullptr, Form("N - CB at m_{#eta}: %.1f #sigma (off scale)",
                                                pmain[iPre].S / pmain[iPre].eS), "");
            l->Draw();
            for (double y : {-2.0, 0.0, 2.0}) { TLine* ln = new TLine(0.38, y, 0.87, y); ln->SetLineStyle(y == 0 ? 1 : 3); ln->SetLineColor(kGray + 1); ln->Draw(); }
            Txt(0.16, 0.935, "(b) frozen peak fit at other positions", 0.054);
        }
        // (c) edge-step heights per slice: N - CB against like-sign
        c->cd(3); gPad->SetLeftMargin(0.15);
        {
            TH1D* fr0 = new TH1D("fredge", "", 4, 0.5, 4.5);
            fr0->SetMinimum(-100); fr0->SetMaximum(780);
            fr0->GetYaxis()->SetTitle("step height h [rows / 2.5 MeV/c^{2}]");
            const char* nm[4] = {"0.30-0.34", "0.34-0.37", "0.37-0.39", "0.39-0.40"};
            for (int i = 0; i < 4; ++i) fr0->GetXaxis()->SetBinLabel(i + 1, nm[i]);
            fr0->GetXaxis()->SetLabelSize(0.060); fr0->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-}) slice [GeV/c^{2}]");
            fr0->Draw("AXIS");
            const int bb[2] = {1, 4};   // pre, bump nominal
            const int cl[2] = {kBlue + 1, kAzure + 7};
            TLegend* l = new TLegend(0.13, 0.60, 0.99, 0.90); l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(0.036);
            for (int q = 0; q < 2; ++q) {
                TGraphErrors* g = new TGraphErrors(); TGraphErrors* gl = new TGraphErrors();
                for (int i = 0; i < 4; ++i) {
                    g->SetPoint(i, i + 1 - 0.12 + 0.08 * q, eNom[bb[q]].h[i]); g->SetPointError(i, 0, eNom[bb[q]].eh[i]);
                    gl->SetPoint(i, i + 1 + 0.08 + 0.08 * q, eLs[bb[q]].h[i]); gl->SetPointError(i, 0, eLs[bb[q]].eh[i]);
                }
                Deco(g, cl[q], q == 0 ? 20 : 24, 1.0); Deco(gl, q == 0 ? kRed : kRed - 7, q == 0 ? 21 : 25, 1.0);
                g->Draw("P same"); gl->Draw("P same");
                l->AddEntry(g, Form("N - CB, %s: E = %.0f #pm %.0f", sels[bIdx[bb[q]]].name.c_str(), eNom[bb[q]].E, eNom[bb[q]].eE), "pe");
                l->AddEntry(gl, Form("N_{++}+N_{--}, %s: E = %.0f #pm %.0f", sels[bIdx[bb[q]]].name.c_str(), eLs[bb[q]].E, eLs[bb[q]].eE), "pe");
            }
            l->Draw();
            TLine* ln = new TLine(0.5, 0, 4.5, 0); ln->SetLineStyle(2); ln->SetLineColor(kGray + 1); ln->Draw();
            Txt(0.16, 0.935, "(c) bump edge steps: e^{+}e^{-} vs like-sign", 0.054);
        }
        Save(c, "cb_null");
        c->Write();
    }
    // ---- Figure 4: sensitivity to the CB variants ----
    {
        TCanvas* c = new TCanvas("cb_sensitivity", "", 1500, 560);
        c->Divide(3, 1, 0.002, 0.002);
        const char* vlab[nVar] = {"#times0.97", "#times1.03", "2.5 MeV", "5 MeV", "20 MeV", "step", "cubic", "k_{mix}", "LS shape", "N_{++}+N_{--}"};
        const char* yt[3] = {"#DeltaS / S  [%]   (#eta peak, P_{main})", "#DeltaE / E  [%]   (bump edge steps)", "#DeltaW / W  [%]   (N - CB, bump window)"};
        const char* tt[3] = {"(a) #eta-peak content S", "(b) bump edge-step sum E", "(c) N - CB in the bump window"};
        const char* fk3[] = {"pre", "peakvtx", "bumpvtx", "pre_b1"};
        const int fcol[] = {kBlack, kBlue + 1, kGreen + 2, kMagenta + 1};
        const int fmk[] = {20, 21, 22, 30};
        for (int pnl = 0; pnl < 3; ++pnl) {
            c->cd(pnl + 1); gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.20);
            TH1D* fr0 = new TH1D(Form("frsen%d", pnl), "", nVar, -0.5, nVar - 0.5);
            fr0->SetMinimum(-7.0); fr0->SetMaximum(7.0);
            for (int q = 0; q < nVar; ++q) fr0->GetXaxis()->SetBinLabel(q + 1, vlab[q]);
            fr0->GetXaxis()->SetLabelSize(0.058); fr0->GetXaxis()->LabelsOption("v");
            fr0->GetYaxis()->SetTitle(yt[pnl]); fr0->GetYaxis()->SetTitleOffset(1.15);
            fr0->Draw("AXIS");
            for (double y : {-3.0, 0.0, 3.0}) { TLine* ln = new TLine(-0.5, y, nVar - 0.5, y); ln->SetLineStyle(y == 0 ? 1 : 3); ln->SetLineColor(kGray + 1); ln->Draw(); }
            TLegend* l = new TLegend(0.17, 0.69, 0.97, 0.91); l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(0.037);
            const char* shortName[4] = {"pre", "peak nominal + vertex", "bump nominal + vertex", "pre, isBest = 1"};
            for (int q4 = 0; q4 < 4; ++q4) {
                TGraph* g = new TGraph();
                int ks = -1, kb = -1;
                for (int k = 0; k < nSens; ++k) if (std::string(sensKeys[k]) == fk3[q4]) ks = k;
                for (int b = 0; b < nB; ++b) if (std::string(bumpKeys[b]) == fk3[q4]) kb = b;
                for (int q = 0; q < nVar; ++q) {
                    double y = 0; bool ok = true;
                    if (pnl == 0) { ok = ks >= 0 && (q != iKmix || hasK[ks]); if (ok) y = 100. * dSv[ks][q] / pmain[sIdx[ks]].S; }
                    else if (pnl == 1) { ok = kb >= 0 && (q != iKmix || hasKb[kb]); if (ok) y = 100. * dEv[kb][q] / eNom[kb].E; }
                    else { ok = ks >= 0 && (q != iKmix || hasK[ks]); if (ok) y = 100. * dWv[ks][q]; }
                    if (ok) g->SetPoint(g->GetN(), q + (q4 - 1.5) * 0.12, y);
                }
                if (g->GetN() == 0) continue;
                Deco(g, fcol[q4], fmk[q4], 1.0); g->Draw("P same");
                const int si = Find(sels, fk3[q4]);
                double rel = 0;
                if (pnl == 1 && kb >= 0) rel = 100. * eNom[kb].eE / eNom[kb].E;
                else if (ks >= 0 && pnl == 0) rel = 100. * pmain[sIdx[ks]].eS / pmain[sIdx[ks]].S;
                else if (ks >= 0) rel = 100. * sp[sIdx[ks]].win[1].sigErr() / sp[sIdx[ks]].win[1].sig();
                (void)si;
                l->AddEntry(g, Form("%s  (stat. %.1f %%)", shortName[q4], rel), "p");
            }
            l->Draw();
            Txt(0.16, 0.935, tt[pnl], 0.048);
        }
        Save(c, "cb_sensitivity");
        c->Write();
    }
    fr.Close();
    if (fk) fk->Close();

    // ======================================================================
    // LaTeX tables
    // ======================================================================
    {
        FILE* t = fopen(Form("%s/cb_share.tex", kTabOut), "w");
        fprintf(t, "%% written by plots/exp_report/cb_checks.C\n");
        fprintf(t, "\\begin{table}[H]\\centering\\footnotesize\\setlength{\\tabcolsep}{3pt}\n");
        fprintf(t, "\\caption{Combinatorial background share CB/$N$ in regions of $M(\\pi^+\\pi^-e^+e^-)$ [GeV/$c^2$], experimental data. "
                   "CB $=2\\sqrt{N_{++}N_{--}}$ in 10~MeV/$c^2$ bins, linearly interpolated, $k=1$. Rows are all stored combinations "
                   "(isBest $\\geq 0$) unless marked isBest $=1$; corrected momenta (cor\\_ok $=1$), PT3, vertex $z>-500$~mm; ``pre'' adds $v_z>-160$~mm; "
                   "``vertex'' = both pions in the event-vertex fit. Peak window 0.53364--0.55416, bump window 0.4550--0.5225 (data windows). "
                   "Last column: $N_{--}/N_{++}$ in the peak window. Not normalised to luminosity.}\n");
        fprintf(t, "\\label{tab:cb_share}\n\\resizebox{\\textwidth}{!}{%%\n\\begin{tabular}{lrrrrrrrr}\n\\toprule\n");
        fprintf(t, "selection & 0.35--0.455 & bump & 0.5225--0.5336 & peak & 0.5542--0.60 & 0.60--0.70 & 0.70--0.95 & $N_{--}/N_{++}$\\\\\n\\midrule\n");
        const int rowsT[] = {iFat, iPre, iPreOa, iPreVtx, iDir, iPkSim, iPk, iPkV, iPkOa, iBm, iBmV, iBmOa, iX17, iPreB, iPkB, iPkVB, iBmB};
        for (int k : rowsT) {
            fprintf(t, "%s", TexName(sels[k].name).c_str());
            for (int q = 0; q < 7; ++q) fprintf(t, " & %.3f", sp[k].win[kRegOrder[q]].cbShare());
            const WinSum& w = sp[k].win[0];
            fprintf(t, " & %.2f\\\\\n", w.n[1] > 0 ? w.n[2] / w.n[1] : 0.);
        }
        fprintf(t, "\\bottomrule\n\\end{tabular}}\n\\end{table}\n");
        fclose(t);
    }
    {
        FILE* t = fopen(Form("%s/cb_null.tex", kTabOut), "w");
        fprintf(t, "%% written by plots/exp_report/cb_checks.C\n");
        fprintf(t, "\\begin{table}[H]\\centering\\footnotesize\\setlength{\\tabcolsep}{3pt}\n");
        fprintf(t, "\\caption{Null tests of the like-sign spectra, experimental data, isBest $\\geq 0$ unless marked. $S$ = $\\eta$-peak content "
                   "from $P_\\text{main}$ ($N-\\text{CB}$, frozen Gaussian 0.5439~GeV/$c^2$, $\\sigma=5.4$~MeV, pol2, 0.5225--0.6100) with stat.\\ error. "
                   "$S(++)$, $S(--)$, $S(\\Sigma)$: the same frozen fit (Poisson likelihood) on $N_{++}$, $N_{--}$ and $N_{++}+N_{--}$. "
                   "Pulls: mean and rms of $S/\\sigma_S$ of $N_{++}+N_{--}$ for the test Gaussian at 46 positions 0.40--0.85~GeV/$c^2$. "
                   "$E$: bump edge-step sum on $N-\\text{CB}$; $E_\\text{ls}$: the same fit on $N_{++}+N_{--}$ (rows per 2.5~MeV/$c^2$). "
                   "Not normalised to luminosity.}\n");
        fprintf(t, "\\label{tab:cb_null}\n\\resizebox{\\textwidth}{!}{%%\n\\begin{tabular}{lrrrrrr}\n\\toprule\n");
        fprintf(t, "selection & $S$ & $S(++)$ & $S(--)$ & $S(\\Sigma)$ & pulls mean/rms & $E$ / $E_\\text{ls}$\\\\\n\\midrule\n");
        for (int k : {iFat, iPre, iPreOa, iPreVtx, iDir, iPkSim, iPk, iPkV, iBm, iBmV, iPreB}) {
            double m = 0, q2 = 0;
            for (int i = 0; i < nLsPos; ++i) { m += lsPull[k][i]; q2 += lsPull[k][i] * lsPull[k][i]; }
            int b = -1; for (int q = 0; q < nB; ++q) if (bIdx[q] == k) b = q;
            fprintf(t, "%s & $%.0f\\pm%.0f$ & $%.0f\\pm%.0f$ & $%.0f\\pm%.0f$ & $%.0f\\pm%.0f$ & $%.2f$ / $%.2f$ & %s\\\\\n", TexName(sels[k].name).c_str(),
                    pmain[k].S, pmain[k].eS, lsP[k].S, lsP[k].eS, lsM[k].S, lsM[k].eS, lsS[k].S, lsS[k].eS, m / nLsPos, std::sqrt(q2 / nLsPos),
                    b >= 0 ? Form("$%.0f\\pm%.0f$ / $%.0f\\pm%.0f$", eNom[b].E, eNom[b].eE, eLs[b].E, eLs[b].eE) : "--");
        }
        fprintf(t, "\\bottomrule\n\\end{tabular}}\n\\end{table}\n");
        fclose(t);
    }
    {
        FILE* t = fopen(Form("%s/cb_sens.tex", kTabOut), "w");
        fprintf(t, "%% written by plots/exp_report/cb_checks.C\n");
        fprintf(t, "\\begin{table}[H]\\centering\\footnotesize\\setlength{\\tabcolsep}{2.2pt}\n");
        fprintf(t, "\\caption{Change of the $\\eta$-peak content $S$ ($P_\\text{main}$) and of the bump edge-step sum $E$, in per cent, when the CB is "
                   "changed; experimental data, isBest $\\geq 0$ unless marked. Columns: CB scaled by 0.97 and 1.03; CB bins of 2.5, 5 and 20~MeV/$c^2$ "
                   "instead of 10; step and cubic instead of linear interpolation; CB multiplied by the event-mixing $k(M)$ of the selection "
                   "(Table~\\ref{tab:cb_kmix}; n/a: none made); like-sign sum shape normalised to the integrated geometric mean; plain like-sign "
                   "sum $N_{++}+N_{--}$ (extreme, not an estimator). Last columns: statistical error and model spread of $S$ (of $E$: stat.\\ only), "
                   "in per cent. Not normalised to luminosity.}\n");
        fprintf(t, "\\label{tab:cb_sens}\n\\resizebox{\\textwidth}{!}{%%\n\\begin{tabular}{lrrrrrrrrrrrr}\n\\toprule\n");
        fprintf(t, "selection & $\\times0.97$ & $\\times1.03$ & 2.5 & 5 & 20 & step & cubic & $k_\\text{mix}$ & LS shape & $N_{++}+N_{--}$ & stat & spread\\\\\n\\midrule\n");
        fprintf(t, "\\multicolumn{13}{l}{$\\eta$-peak content $S$:}\\\\\n");
        for (int k = 0; k < nSens; ++k) {
            const int s = sIdx[k];
            fprintf(t, "%s", TexName(sels[s].name).c_str());
            for (int q = 0; q < nVar; ++q) {
                if (q == iKmix && !hasK[k]) { fprintf(t, " & n/a"); continue; }
                fprintf(t, " & %+.1f", 100.0 * dSv[k][q] / pmain[s].S);
            }
            fprintf(t, " & %.1f & %.1f\\\\\n", 100.0 * pmain[s].eS / pmain[s].S, 100.0 * spread[s].half / pmain[s].S);
        }
        fprintf(t, "\\midrule\n\\multicolumn{13}{l}{bump edge-step sum $E$:}\\\\\n");
        for (int b = 0; b < nB; ++b) {
            fprintf(t, "%s", TexName(sels[bIdx[b]].name).c_str());
            for (int q = 0; q < nVar; ++q) {
                if (q == iKmix && !hasKb[b]) { fprintf(t, " & n/a"); continue; }
                fprintf(t, " & %+.1f", 100.0 * dEv[b][q] / eNom[b].E);
            }
            fprintf(t, " & %.1f & --\\\\\n", 100.0 * eNom[b].eE / eNom[b].E);
        }
        fprintf(t, "\\midrule\n\\multicolumn{13}{l}{$N-\\text{CB}$ in the bump window ($W$):}\\\\\n");
        for (int k = 0; k < nSens; ++k) {
            const int s = sIdx[k];
            fprintf(t, "%s", TexName(sels[s].name).c_str());
            for (int q = 0; q < nVar; ++q) {
                if (q == iKmix && !hasK[k]) { fprintf(t, " & n/a"); continue; }
                fprintf(t, " & %+.1f", 100.0 * dWv[k][q]);
            }
            fprintf(t, " & %.1f & --\\\\\n", 100.0 * sp[s].win[1].sigErr() / sp[s].win[1].sig());
        }
        fprintf(t, "\\bottomrule\n\\end{tabular}}\n\\end{table}\n");
        fclose(t);
    }
    P("\nrun time %.1f s\n", swAll.RealTime());
    fclose(OutFile()); OutFile() = nullptr;
}
