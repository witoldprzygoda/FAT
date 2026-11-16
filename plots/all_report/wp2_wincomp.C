// ============================================================================
// wp2_wincomp.C  --  WP2: composition of the empty-target rows that pass S0, and the bounds on b
// ============================================================================
// Plain words: the S0 rows of the empty-target runs (vz > -160) come from three places: the entrance window
// (vz -160..-145 mm), the inside of the cell (-145..-105 mm) and the exit window and beyond (> -105 mm). Rows on the
// beam axis inside the cell could come from H2 gas left in the "empty" cell; in the full-target runs that gas is
// replaced by the liquid, so subtracting it would be wrong. The radial distance r of the vertex from the beam spot
// separates gas (same r profile as LH2: beam-like) from walls and foils at larger r (r profile of the material in
// front of the cell, -250..-170 mm). The on-axis share g of the in-cell empty rows, from the fraction with r < 4 mm,
//     g = (p_cell - p_up) / (p_LH2 - p_up),
// gives an upper limit of the possible gas over-subtraction g x (in-cell part of f_win).
// Bounds on the flux ratio b: b <= F/E in the entrance-window bins [-170,-160) (the full rows there are at least
// the material); the fit wp2_vzfit.C gives the lower variant.
// Output: docs/latex/tables_all/wp2_wincomp.txt
// Run: root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_wp2", kTRUE)' \
//           'plots/all_report/wp2_wincomp.C+O'
// ============================================================================
#include <cmath>
#include <cstdio>

#include "defs.h"

namespace wp2c {
const char* const kEmpty = "/hdd3/przygoda/HADES/FAT_out/all_report/wp2/empty/scan2/scan_compact_empty.root";
const int kDays[3] = {34, 47, 65};
const double kLper[3] = {2809.453, 372.610 + 1078.022, 1819.110};
inline int DI(int d) { return d == 34 ? 0 : d == 47 ? 1 : d == 65 ? 2 : -1; }
enum { zEnt = 0, zCell, zExit, kNZ };
const char* kZN[kNZ] = {"entrance [-160,-145)", "cell [-145,-105)", "exit [-105,...)"};
enum { wAll = 0, wB, wP, wBP, kNW };
const char* kWN[kNW] = {"all M4", "B", "P", "B+P"};
}  // namespace wp2c

void wp2_wincomp() {
    using namespace wp2c;
    scan::Data F, E;
    if (!scan::Load(F) || !scan::Load(E, kEmpty)) return;
    // beam spot from the full-target LH2 rows of the three days
    double sx = 0, sy = 0, sn = 0;
    for (const auto& r : F.s[0]) {
        if (DI(r.file_id) < 0 || !alldefs::PassBase(r, r.cor_ok == 1) || !alldefs::PiVtx(r)) continue;
        if (r.vz > -140 && r.vz < -110) { sx += r.vx; sy += r.vy; sn += 1; }
    }
    sx /= sn; sy /= sn;
    // accumulators [target][sample][day][z region][window]; on-axis counts [target][day][region 0 up,1 cell,2 LH2-centre]
    double n[2][3][3][kNZ][kNW] = {}, U[2][3] = {}, Fe[2][3] = {}, ax[2][3][3] = {}, axn[2][3][3] = {};
    for (int t = 0; t < 2; ++t) {
        const scan::Data& D = t ? E : F;
        for (int s = 0; s < 3; ++s)
            for (const auto& r : D.s[s]) {
                const int d = DI(r.file_id);
                if (d < 0 || !alldefs::PassBase(r, r.cor_ok == 1)) continue;
                if (r.vz >= -200.f && r.vz < -170.f) U[t][d] += 1;   // nominal normaliser ch[-200,-170)
                if (!alldefs::PiVtx(r)) continue;
                if (s == 0 && r.vz >= -170.f && r.vz < -160.f) Fe[t][d] += 1;
                if (s == 0) {
                    const double rr = std::hypot(r.vx - sx, r.vy - sy);
                    const int reg = (r.vz >= -250 && r.vz < -170) ? 0 : (r.vz >= -145 && r.vz < -105) ? 1 : -1;
                    if (reg >= 0) { axn[t][d][reg] += 1; if (rr < 4) ax[t][d][reg] += 1; }
                    if (r.vz >= -140 && r.vz < -110) { axn[t][d][2] += 1; if (rr < 4) ax[t][d][2] += 1; }
                }
                if (r.vz <= -160.f) continue;
                const int z = r.vz < -145.f ? zEnt : r.vz < -105.f ? zCell : zExit;
                n[t][s][d][z][wAll] += 1;
                if (r.m4 >= alldefs::kBlo && r.m4 < alldefs::kBhi) n[t][s][d][z][wB] += 1;
                if (r.m4 >= alldefs::kPlo && r.m4 < alldefs::kPhi) n[t][s][d][z][wP] += 1;
                if ((r.m4 >= alldefs::kBlo && r.m4 < alldefs::kBhi) || (r.m4 >= alldefs::kPlo && r.m4 < alldefs::kPhi)) n[t][s][d][z][wBP] += 1;
            }
    }
    FILE* o = fopen("/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all/wp2_wincomp.txt", "w");
    fprintf(o, "# wp2_wincomp.C: composition of f_win (S0, N - CB) by vertex region, gas bound, bounds on b; beam spot x %.2f y %.2f mm\n", sx, sy);
    fprintf(o, "\n[b bounds] day: b_ch = U_F/U_E (ch[-200,-170)) ; b_max = F/E of the unlike-sign S0-without-vz rows in [-170,-160)\n");
    double b[3];
    for (int d = 0; d < 3; ++d) {
        b[d] = U[0][d] / U[1][d];
        const double bm = Fe[0][d] / Fe[1][d];
        fprintf(o, "  %03d  b_ch %.3f +- %.3f   b_max %.3f +- %.3f   b_max/b_ch %.3f\n", kDays[d], b[d],
                b[d] * std::sqrt(1 / U[0][d] + 1 / U[1][d]), bm, bm * std::sqrt(1 / Fe[0][d] + 1 / Fe[1][d]), bm / b[d]);
    }
    const double Ls = kLper[0] + kLper[1] + kLper[2];
    fprintf(o, "\n[on-axis share] rows with r < 4 mm (S0 without vz cut, unlike-sign), all three days\n");
    double p[2][3];
    for (int t = 0; t < 2; ++t)
        for (int reg = 0; reg < 3; ++reg) {
            double a = 0, nn = 0;
            for (int d = 0; d < 3; ++d) { a += ax[t][d][reg]; nn += axn[t][d][reg]; }
            p[t][reg] = a / nn;
            fprintf(o, "  %-5s %-28s p(r<4) = %.3f (n %.0f)\n", t ? "empty" : "full", reg == 0 ? "front material [-250,-170)" : reg == 1 ? "cell [-145,-105)" : "cell centre [-140,-110)",
                    p[t][reg], nn);
        }
    const double g = (p[1][1] - p[1][0]) / (p[0][2] - p[1][0]);
    fprintf(o, "  on-axis (gas-like) share of the in-cell empty rows g = (p_cell,E - p_front,E)/(p_LH2,F - p_front,E) = %.3f\n", g);
    fprintf(o, "\n[composition] f_win [%%] (N - CB, S0, nominal b) split by vertex region, L-weighted over the three days; gas bound = g x cell part\n");
    for (int w = 0; w < kNW; ++w) {
        double fz[kNZ] = {}, ftot = 0;
        for (int d = 0; d < 3; ++d) {
            double YF = 0;
            double ye[kNZ];
            for (int z = 0; z < kNZ; ++z) {
                auto nmc = [&](int t) {
                    const double nu = n[t][0][d][z][w], np = n[t][1][d][z][w], nm = n[t][2][d][z][w];
                    return nu - ((np > 0 && nm > 0) ? 2 * std::sqrt(np * nm) : 0.);
                };
                ye[z] = nmc(1);
                // full yield: N - CB over all regions (CB of the whole vz > -160 range)
            }
            double nu = 0, np = 0, nm = 0;
            for (int z = 0; z < kNZ; ++z) { nu += n[0][0][d][z][w]; np += n[0][1][d][z][w]; nm += n[0][2][d][z][w]; }
            YF = nu - 2 * std::sqrt(np * nm);
            for (int z = 0; z < kNZ; ++z) fz[z] += kLper[d] / Ls * b[d] * ye[z] / YF;
        }
        for (int z = 0; z < kNZ; ++z) ftot += fz[z];
        fprintf(o, "  %-7s total %6.2f :", kWN[w], 100 * ftot);
        for (int z = 0; z < kNZ; ++z) fprintf(o, "  %s %6.2f", kZN[z], 100 * fz[z]);
        fprintf(o, "   gas bound %5.2f\n", 100 * g * fz[zCell]);
        fprintf(o, "wc|%s,%.4f,%.4f,%.4f,%.4f,%.4f\n", kWN[w], 100 * ftot, 100 * fz[zEnt], 100 * fz[zCell], 100 * fz[zExit], 100 * g * fz[zCell]);
    }
    fprintf(o, "# note: the like-sign CB of the empty sample is taken per vertex region (small counts; the per-region split is indicative)\n");
    fclose(o);
    FILE* in = fopen("/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all/wp2_wincomp.txt", "r");
    char line[512];
    while (fgets(line, sizeof line, in)) printf("%s", line);
    fclose(in);
}
