// ============================================================================
// x17c_disp.C  --  WP9 X-c: event-group covariance of the data rows in the X2D SIDE-BAND cells (S0)
// ============================================================================
// SPEC report_all, WP9 X-c (prefix x17c). Written 2026-11-14.
//
// Why: one event can give several rows (all stored combinations, isBest >= 0), so the cell counts are not independent
// Poisson numbers. The X17 likelihood is a Poisson likelihood of the rows; its test statistic must be divided by the
// variance factor phi = w^T S w / w^T diag(N) w of the linear combination the X17 fit uses (w = X17 template / model).
// This macro builds S_ab = sum over event groups g of n_ga n_gb for the X2D cells of the M(ee) SIDE BANDS only.
//
// BLINDING (D19): rows with 12 <= M(ee) < 28 MeV/c2 (the 2 MeV bins that overlap 13-27 MeV) are dropped as the very
// first operation on every row, before the selection and before any counter; nothing of them is stored or printed.
//
// Cells: X2D M4 0.5225-0.6100 (35 bins of 2.5 MeV, a row's bin = the X2D bin of the centre of its 2.5 MeV fine bin,
// exactly as wp0_fill.C) x side-band M(ee) rows 2-12 MeV (5 rows) and 28-62 MeV (17 rows) = 22 x 35 = 770 cells.
// Rows grouped by ADJACENCY of grp (never a global map), as wp3_rows.C / wp5_rows.C.
// Check: the cell counts and the diagonal equal N and err^2 of hists_v1.root data/S0/all/N_X2D (side-band rows).
//
// Output: /hdd3/przygoda/HADES/FAT_out/all_report/x17c/x17c_egcov.root (TMatrixDSym egcov_sb, TVectorD cnt_sb, index map)
// Run (from /home/przygoda/HADES/FAT/FAT; < 3 GB, about 1 min):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_x17c", kTRUE)' \
//        'plots/all_report/x17c_disp.C+O'
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

#include "TFile.h"
#include "TH2.h"
#include "TMatrixDSym.h"
#include "TSystem.h"
#include "TVectorD.h"

#include "hists_io.h"

namespace x17cd {
using namespace alldefs;
const char* const kOutDir = "/hdd3/przygoda/HADES/FAT_out/all_report/x17c";
const double kBlindLo = 0.012, kBlindHi = 0.028;   // GeV/c2: never read
// side-band row index (0..21) of an M(ee) value, -1 outside; X2D row j (0..29) returned in jx2d
inline int SbRow(double mee, int& jx2d) {
    jx2d = -1;
    if (!(mee >= kX2Dylo && mee < kX2Dyhi)) return -1;
    const int j = std::min(kX2Dny - 1, (int)((mee - kX2Dylo) / 0.002));
    jx2d = j;
    if (j <= 4) return j;            // 2-12 MeV
    if (j >= 13) return j - 8;       // 28-62 MeV -> 5..21
    return -1;
}
inline int X2Dx(double m4) {
    const int fi = FineBin(m4);
    if (fi < 0) return -1;
    const double c = kFineLo + (fi + 0.5) * kFineW;
    if (c < kX2Dxlo || c >= kX2Dxhi) return -1;
    return std::min(kX2Dnx - 1, (int)((c - kX2Dxlo) / 0.0025));
}
}  // namespace x17cd

void x17c_disp() {
    using namespace alldefs;
    using namespace x17cd;
    gSystem->mkdir(kOutDir, kTRUE);
    scan::Data d;
    if (!scan::Load(d, scan::kCompact, 0.40, 0.75, 1)) return;   // us rows only, generous M4 range on m4_fb
    const int nrow = 22, nc = nrow * kX2Dnx;
    std::vector<double> S((size_t)nc * nc, 0.), cnt(nc, 0.);
    std::vector<int> cells;
    long ngroups = 0, nmulti = 0;
    long nBlindSkipped = 0;   // only a counter of skipped rows, no content
    auto flush = [&]() {
        if (cells.empty()) return;
        ++ngroups;
        if (cells.size() > 1) ++nmulti;
        for (size_t a = 0; a < cells.size(); ++a)
            for (size_t b = 0; b < cells.size(); ++b) S[(size_t)cells[a] * nc + cells[b]] += 1.0;
        cells.clear();
    };
    int cur = -1;
    for (const auto& r : d.us()) {
        if (r.grp != cur) { flush(); cur = r.grp; }
        // BLINDING: first operation on the row
        if (r.mee >= kBlindLo && r.mee < kBlindHi) { ++nBlindSkipped; continue; }
        if (!PassSel(kSelS0, r, r.cor_ok == 1, r.dz_pi, r.mm4)) continue;
        int jx;
        const int jr = SbRow(r.mee, jx);
        if (jr < 0) continue;
        const int ix = X2Dx(r.m4);
        if (ix < 0) continue;
        const int c = jr * kX2Dnx + ix;
        cells.push_back(c);
        cnt[c] += 1;
    }
    flush();
    // check against hists_v1.root (side-band rows only)
    TFile* fh = TFile::Open(kHistFile, "READ");
    TH2* N = (TH2*)fh->Get("data/S0/all/N_X2D");
    double maxd = 0, maxe = 0, tot = 0, diag = 0, all = 0;
    for (int jr = 0; jr < nrow; ++jr) {
        const int j = jr <= 4 ? jr : jr + 8;   // X2D row (0-based)
        for (int i = 0; i < kX2Dnx; ++i) {
            const int c = jr * kX2Dnx + i;
            maxd = std::max(maxd, std::fabs(N->GetBinContent(i + 1, j + 1) - cnt[c]));
            maxe = std::max(maxe, std::fabs(std::pow(N->GetBinError(i + 1, j + 1), 2) - S[(size_t)c * nc + c]));
            tot += cnt[c];
        }
    }
    for (int a = 0; a < nc; ++a) for (int b = 0; b < nc; ++b) { all += S[(size_t)a * nc + b]; if (a == b) diag += S[(size_t)a * nc + b]; }
    printf("x17c_disp: side-band cells %d, rows %.0f, groups %ld (%ld with > 1 side-band row); max|N - hists_v1| %.3g, max|diag - err^2| %.3g\n",
           nc, tot, ngroups, nmulti, maxd, maxe);
    printf("x17c_disp: var(total)/N = %.4f (within-cell part %.4f); rows skipped by the blinding filter (count only): %ld\n", all / tot, diag / tot,
           nBlindSkipped);
    TFile* fo = TFile::Open(Form("%s/x17c_egcov.root", kOutDir), "RECREATE");
    TMatrixDSym M(nc);
    for (int a = 0; a < nc; ++a) for (int b = 0; b < nc; ++b) M(a, b) = S[(size_t)a * nc + b];
    TVectorD v(nc);
    for (int a = 0; a < nc; ++a) v(a) = cnt[a];
    M.Write("egcov_sb");
    v.Write("cnt_sb");
    TVectorD rows(nrow);   // X2D row (0-based) of every side-band row index
    for (int jr = 0; jr < nrow; ++jr) rows(jr) = jr <= 4 ? jr : jr + 8;
    rows.Write("x2d_row_of_sb");
    fo->Close();
    FILE* t = fopen(Form("%s/x17c_disp.txt", kOutDir), "w");
    fprintf(t, "x17c_disp: S0, X2D side-band cells (M(ee) 2-12 and 28-62 MeV/c2), %d cells, rows %.0f, event groups %ld (%ld with > 1 row)\n", nc, tot, ngroups,
            nmulti);
    fprintf(t, "check vs hists_v1 N_X2D: max|N diff| %.3g, max|err^2 diff| %.3g\n", maxd, maxe);
    fprintf(t, "var(total)/N %.4f, within-cell %.4f\n", all / tot, diag / tot);
    fclose(t);
    fh->Close();
}
