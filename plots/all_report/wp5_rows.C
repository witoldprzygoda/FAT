// ============================================================================
// wp5_rows.C  --  WP5 (main extraction, stage B): event-group covariances of the data N
// ============================================================================
// SPEC report_all, WP5. Written by WP5 (prefix wp5), 2026-11-14.
//
// The stage-B fits use a sandwich error with the event-group covariance of the observed counts
// (one event group can give several rows = several cells). wp3_rows.C made it for the F10 cells of
// S0 and pre (all data). WP5 also needs it for
//   * the F5g cells (M(ee) >= 20 MeV; F-gold): S0, pre, S0_best;
//   * the subsets of S0: fold0, fold1 (event-group parity), periods I, IIA, IIB, III;
//   * the data with isBest == 1 (S0_best: all, fold0, fold1);
//   * pre: fold0, fold1.
//   Sigma_ab = sum over event groups g of n_ga n_gb (n_ga = selected rows of group g in cell a).
//   Rows are grouped by ADJACENCY of grp (never by a global map), exactly as wp3_rows.C.
// Check: the cell counts and the diagonal equal N and err^2 of hists_v1.root, and the S0/pre F10
// matrices equal those of wp3_rows.root.
//
// Output: /hdd3/przygoda/HADES/FAT_out/all_report/wp5/wp5_egcov.root (TMatrixDSym egcov_<sel>_<sub>_<F10|F5g>)
// Run (from /home/przygoda/HADES/FAT/FAT; 2.7 GB, about 1 min):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_wp5", kTRUE)' \
//        'plots/all_report/wp5_rows.C+O'
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1D.h"
#include "TMatrixDSym.h"
#include "TSystem.h"

#include "hists_io.h"

namespace wp5r {
using namespace alldefs;

const char* const kOutDir = "/hdd3/przygoda/HADES/FAT_out/all_report/wp5";

std::vector<double> gF10e;
// M4 bin (0..69) of a row on the F10 grid (identical to wp0_fill.C / wp3_rows.C: fine bin first, then its centre)
inline int F10x(const scan::Row& r) {
    const int fi = FineBin(r.m4);
    if (fi < 0) return -1;
    const double c = kFineLo + (fi + 0.5) * kFineW;
    if (c < gF10e.front() || c >= gF10e.back()) return -1;
    return (int)(std::upper_bound(gF10e.begin(), gF10e.end(), c) - gF10e.begin()) - 1;
}
inline int F10Cell(const scan::Row& r) {
    const int x = F10x(r);
    if (x < 0) return -1;
    const int s = Slice(r.mpipi), e = Eclass(r.mee);
    if (s < 0 || e < 0) return -1;
    return (2 * s + e) * kNF10 + x;
}
inline int F5gCell(const scan::Row& r) {
    const int x = F10x(r);
    if (x < 0) return -1;
    const int s = Slice(r.mpipi);
    if (s < 0 || !(r.mee >= kF5gMee)) return -1;
    return s * kNF10 + x;
}

struct Acc {
    std::string sel, sub, lay;
    int selId = 0;        // alldefs selection
    bool best = false;    // isBest == 1 in addition
    int subId = 0;        // kSubAll, kSubFold0, ...
    bool f5g = false;
    int n = 0;
    std::vector<double> S, cnt;
    std::vector<int> cells;
    long ngroups = 0, nmulti = 0;
    void Init() { n = f5g ? 5 * kNF10 : 10 * kNF10; S.assign((size_t)n * n, 0.); cnt.assign(n, 0.); }
    void Flush() {
        if (cells.empty()) return;
        ++ngroups;
        if (cells.size() > 1) ++nmulti;
        for (size_t a = 0; a < cells.size(); ++a)
            for (size_t b = 0; b < cells.size(); ++b) S[(size_t)cells[a] * n + cells[b]] += 1.0;
        cells.clear();
    }
};

}  // namespace wp5r

void wp5_rows() {
    using namespace alldefs;
    using namespace wp5r;
    gSystem->mkdir(kOutDir, kTRUE);
    gF10e = F10Edges();
    FILE* out = fopen(Form("%s/wp5_rows.txt", kOutDir), "w");
    scan::Data d;
    if (!scan::Load(d)) return;

    std::vector<Acc> acc;
    auto add = [&](const char* sel, int selId, bool best, int subId, bool f5g) {
        Acc a; a.sel = sel; a.selId = selId; a.best = best; a.subId = subId; a.f5g = f5g; a.sub = SubTag(subId);
        a.lay = f5g ? "F5g" : "F10"; a.Init(); acc.push_back(a);
    };
    // S0: every subset, F10 and F5g
    for (int s = 0; s < kNSub; ++s) { add("S0", kSelS0, false, s, false); add("S0", kSelS0, false, s, true); }
    // S0_best: all and folds
    for (int s : {(int)kSubAll, (int)kSubFold0, (int)kSubFold1}) { add("S0_best", kSelS0, true, s, false); add("S0_best", kSelS0, true, s, true); }
    // pre: all (F10 as a check against wp3, F5g) and folds
    for (int s : {(int)kSubAll, (int)kSubFold0, (int)kSubFold1}) { add("pre", kSelPre, false, s, false); add("pre", kSelPre, false, s, true); }

    int cur = -1;
    for (const auto& r : d.us()) {
        if (r.grp != cur) { for (auto& a : acc) a.Flush(); cur = r.grp; }
        const bool pS0 = PassSel(kSelS0, r, r.cor_ok == 1, r.dz_pi, r.mm4);
        const bool pPre = PassSel(kSelPre, r, r.cor_ok == 1, r.dz_pi, r.mm4);
        if (!pS0 && !pPre) continue;
        const int c10 = F10Cell(r), c5 = F5gCell(r);
        if (c10 < 0 && c5 < 0) continue;
        const int fold = kSubFold0 + (r.fold & 1);
        const int per = PeriodOfDay(r.file_id);
        const int perSub = per >= 0 ? kSubPerI + per : -1;
        for (auto& a : acc) {
            const bool ps = a.selId == kSelS0 ? pS0 : pPre;
            if (!ps) continue;
            if (a.best && r.isBest != 1) continue;
            if (a.subId != kSubAll && a.subId != fold && a.subId != perSub) continue;
            const int c = a.f5g ? c5 : c10;
            if (c < 0) continue;
            a.cells.push_back(c);
            a.cnt[c] += 1;
        }
    }
    for (auto& a : acc) a.Flush();

    // checks against hists_v1.root and wp3_rows.root; write
    TFile* fh = TFile::Open(kHistFile);
    TFile* f3 = TFile::Open("/hdd3/przygoda/HADES/FAT_out/all_report/wp3/wp3_rows.root");
    TFile* fo = TFile::Open(Form("%s/wp5_egcov.root", kOutDir), "RECREATE");
    for (auto& a : acc) {
        double maxd = 0, sumN = 0, diag = 0, off = 0;
        const int nh = a.f5g ? 5 : 10;
        for (int h = 0; h < nh; ++h) {
            const std::string lay = a.f5g ? F5gName(h) : F10Name(h / 2, h % 2);
            TH1* N = allio::Data(fh, a.sel.c_str(), a.sub.c_str(), "N", lay);
            if (!N) continue;
            for (int x = 0; x < kNF10; ++x) {
                const int c = h * kNF10 + x;
                maxd = std::max(maxd, std::fabs(N->GetBinContent(x + 1) - a.cnt[c]));
                maxd = std::max(maxd, std::fabs(N->GetBinError(x + 1) * N->GetBinError(x + 1) - a.S[(size_t)c * a.n + c]));
                sumN += a.cnt[c];
            }
            delete N;
        }
        for (int i = 0; i < a.n; ++i) for (int j = 0; j < a.n; ++j) { if (i == j) diag += a.S[(size_t)i * a.n + j]; else off += a.S[(size_t)i * a.n + j]; }
        double d3 = -1;
        if (a.sub == "all" && !a.f5g && !a.best) {
            TMatrixDSym* M = (TMatrixDSym*)f3->Get(Form("egcov_%s", a.sel.c_str()));
            if (M) { d3 = 0; for (int i = 0; i < a.n; ++i) for (int j = 0; j < a.n; ++j) d3 = std::max(d3, std::fabs((*M)(i, j) - a.S[(size_t)i * a.n + j])); }
        }
        TMatrixDSym S(a.n);
        for (int i = 0; i < a.n; ++i) for (int j = 0; j < a.n; ++j) S(i, j) = a.S[(size_t)i * a.n + j];
        fo->cd();
        S.Write(Form("egcov_%s_%s_%s", a.sel.c_str(), a.sub.c_str(), a.lay.c_str()));
        const char* line = Form("egcov_%s_%s_%s: rows %.0f, groups %ld (%ld with >1 cell), max|N - hists_v1| and |diag - err^2| %.3g, "
                                "var(total)/N %.4f%s\n", a.sel.c_str(), a.sub.c_str(), a.lay.c_str(), sumN, a.ngroups, a.nmulti, maxd,
                                sumN > 0 ? (diag + off) / sumN : 0., d3 >= 0 ? Form(", max|diff to wp3_rows| %.3g", d3) : "");
        printf("%s", line); fprintf(out, "%s", line);
    }
    fo->Close();
    fclose(out);
}
