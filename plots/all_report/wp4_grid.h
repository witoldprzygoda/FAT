// ============================================================================
// wp4_grid.h  --  cell grid of the WP4 stress closure (SPEC report_all, prefix wp4)
// ============================================================================
// The toys are generated on "refined" cells: the 70 F10 M4 bins (0.400-0.700) split
// at every window edge that is not an F10 edge (P 0.53364/0.55416, side bands
// 0.5336/0.5542/0.5628, X17c 0.5312/0.5566), times 15 (M(pipi) slice x M(ee) class)
// cells, class e3: 0 = M(ee) < 10 MeV (E1), 1 = 10-20 MeV, 2 = >= 20 MeV (F5g).
// F10 E2 = e3 1 + 2, F5g = e3 2. Every window of the estimators (B, P, SB, CR-U,
// CR-K) is a union of refined cells, so window counts are exact.
// ============================================================================
#ifndef WP4_GRID_H
#define WP4_GRID_H

#include <algorithm>
#include <cmath>
#include <vector>

#include "defs.h"

namespace wp4 {
using namespace alldefs;

const int kNj = kNSlice;   // 5 slices
const int kNe3 = 3;        // M(ee) classes of the toy cells
const int kNh15 = kNj * kNe3;
inline int H15(int j, int e3) { return j * kNe3 + e3; }
inline int E3class(double mee) { return !(mee == mee) ? -1 : (mee < 0.010 ? 0 : (mee < 0.020 ? 1 : 2)); }
inline int E2of(int e3) { return e3 == 0 ? 0 : 1; }   // F10 class of a toy class

// refined M4 edges
inline std::vector<double> RefEdges() {
    std::vector<double> e = F10Edges();
    const double extra[7] = {0.5312, 0.5336, 0.53364, 0.55416, 0.5542, 0.5566, 0.5628};
    for (double x : extra) e.push_back(x);
    std::sort(e.begin(), e.end());
    return e;
}
// 0 underflow, 1..n, n + 1 overflow (NaN -> overflow)
inline int FindBinE(const std::vector<double>& e, double x) {
    const int n = (int)e.size() - 1;
    if (!(x == x)) return n + 1;
    if (x < e[0]) return 0;
    if (x >= e[n]) return n + 1;
    return (int)(std::upper_bound(e.begin(), e.end(), x) - e.begin());
}

struct Grid {
    std::vector<double> ref, f10;   // edges
    int nR = 0, nF = 0;
    std::vector<int> refToF10;       // refined bin (1..nR) -> F10 bin (1..nF)
    std::vector<double> xc;          // tilt variable x at the refined-cell centre (index 1..nR)
    std::vector<double> xF;          // tilt variable at the F10 bin centre (1..nF)
    void Init() {
        ref = RefEdges(); f10 = F10Edges();
        nR = (int)ref.size() - 1; nF = (int)f10.size() - 1;
        refToF10.assign(nR + 2, -1); xc.assign(nR + 2, 0.); xF.assign(nF + 2, 0.);
        for (int r = 1; r <= nR; ++r) {
            const double c = 0.5 * (ref[r - 1] + ref[r]);
            refToF10[r] = FindBinE(f10, c);
            xc[r] = TiltX(c);
        }
        for (int b = 1; b <= nF; ++b) xF[b] = TiltX(0.5 * (f10[b - 1] + f10[b]));
    }
    // refined cells inside [lo, hi) (edges must be refined edges)
    std::vector<int> RefIn(double lo, double hi) const {
        std::vector<int> v;
        for (int r = 1; r <= nR; ++r)
            if (ref[r - 1] >= lo - 1e-9 && ref[r] <= hi + 1e-9) v.push_back(r);
        return v;
    }
};

// components of the toy templates
enum TComp { kTC2 = 0, kTC3, kTC4, kTC5, kTGeta, kTGDom, kTGDnom, kTGCom, kTGCnom, kNTComp };
inline const char* TCompTag(int c) {
    static const char* n[kNTComp] = {"C2", "C3", "C4", "C5", "Geta", "GDom", "GDnom", "GCom", "GCnom"};
    return n[c];
}
// mass-transfer variants filled by wp4_fill.C
enum MVar { kMnom = 0, kMf, kMsUp, kMr089, kNMVar };
inline const char* MVarTag(int v) { static const char* n[kNMVar] = {"nom", "f", "sUp", "r089"}; return n[v]; }
inline double MVarS(int v) { return (v == kMf || v == kMsUp) ? kS + kSerr : kS; }
inline double MVarR(int v) { return (v == kMf || v == kMr089) ? kRlo : kR; }

const char* const kOutDir = "/hdd3/przygoda/HADES/FAT_out/all_report/wp4";
const char* const kTmplFile = "/hdd3/przygoda/HADES/FAT_out/all_report/wp4/wp4_templates.root";

}  // namespace wp4

#endif
