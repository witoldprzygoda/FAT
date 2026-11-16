// ============================================================================
// c2recon_bedge.C  --  package c2recon: the C2 gamma* factor that the per-slice edge estimator (bedge) actually needs
// ============================================================================
// Written 2026-11-15 (prefix c2recon). New file. It #includes bedge_ana.C to use the bedge fit code as it is
// (namespace bedgea: Evaluate, FitStep with the final data covariance, Quant), reads the bedge arrays and writes
// NOTHING into the bedge directory (bedge_ana() itself is never called).
//
//   1. nominal bedge evaluation (reproduces E_data 1141.86, E_sim 2915.70, the linear C2 / C3 steps per class);
//   2. the C2 template of the nominal configuration rebuilt from c2recon_rows.txt (S0 rows, W = w x RowFactor(108) x
//      27.0718, the bedge bins) -> must equal sim_S.bin [nom][C2] bin by bin;
//   3. the same template with every row x its QED weight (c2recon_weights.txt; density variants A, B, A coarse, flat,
//      A x KLOE-X) substituted into sim_S.bin [nom][C2], and the full bedge evaluation repeated:
//        - the linear C2 step per class (data covariance) -> the factor c_edge = h(C2 QED) / h(C2 iso);
//        - E_sim with the bedge (WP6) own weights -> f = E_data / E_sim(QED) directly;
//        - the bedge "recommended" formula with the WP8 class means 0.879 / 0.882, for comparison;
//   4. the linear functional of every slice x class fit (h of a unit spectrum in each M4 bin, data covariance), for
//      the row-level statistics in c2recon_combine.py.
//
// Output: /hdd3/przygoda/HADES/FAT_out/all_report/c2recon/c2recon_bedge.txt (KEY lines), c2recon_bedge_lin.txt
// Run (from /home/przygoda/HADES/FAT/FAT; about 2 GB, under a minute):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/c2recon/aclic", kTRUE)' \
//        'plots/all_report/c2recon_bedge.C+O'
// ============================================================================
#include "bedge_ana.C"

#include <fstream>
#include <map>
#include <sstream>

namespace c2rb {
using namespace bedgea;
const char* const kDir = "/hdd3/przygoda/HADES/FAT_out/all_report/c2recon";
FILE* gO = nullptr;
#define PO(...) do { printf(__VA_ARGS__); if (gO) fprintf(gO, __VA_ARGS__); fflush(stdout); } while (0)

struct RowB { int tree; long idx; int s0; double W; int e, p, ib; double q[5]; };

bool ReadRows(std::vector<RowB>& v) {
    std::ifstream in(Form("%s/c2recon_rows.txt", kDir));
    if (!in) { PO("cannot open c2recon_rows.txt\n"); return false; }
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream s(line);
        RowB r; double x[22];
        s >> r.tree >> r.idx;
        for (int k = 0; k < 22; ++k) s >> x[k];   // set pat_entry grp s0 best cor_ok W108 W0 W116 W100 m4 m4_true mpipi mpipi_true mee mee_true e p ib tag m5c g_e
        r.s0 = (int)x[3]; r.W = x[6]; r.e = (int)x[16]; r.p = (int)x[17]; r.ib = (int)x[18];
        for (int k = 0; k < 5; ++k) r.q[k] = 1.0;
        v.push_back(r);
    }
    std::ifstream iw(Form("%s/c2recon_weights.txt", kDir));
    if (!iw) { PO("cannot open c2recon_weights.txt\n"); return false; }
    std::map<std::pair<int, long>, size_t> key;
    for (size_t i = 0; i < v.size(); ++i) key[{v[i].tree, v[i].idx}] = i;
    long nj = 0;
    while (std::getline(iw, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream s(line);
        int t; long ix; double q[5];
        s >> t >> ix >> q[0] >> q[1] >> q[2] >> q[3] >> q[4];
        auto it = key.find({t, ix});
        if (it == key.end()) continue;
        for (int k = 0; k < 5; ++k) v[it->second].q[k] = q[k];
        ++nj;
    }
    PO("rows %zu, weights joined %ld\n", v.size(), nj);
    return nj == (long)v.size();
}

void Run() {
    gO = fopen(Form("%s/c2recon_bedge.txt", kDir), "w");
    PO("# c2recon_bedge.C, 2026-11-15: the C2 gamma* factor of the per-slice edge estimator (bedge code)\n");
    if (!ReadArr(Form("%s/data_N.bin", kOut), DN, SizeDN()) || !ReadArr(Form("%s/data_R.bin", kOut), DR, SizeDR()) ||
        !ReadArr(Form("%s/sim_S.bin", kOut), SS, SizeSS()) || !ReadArr(Form("%s/sim_V.bin", kOut), SV, SizeSS()) ||
        !ReadArr(Form("%s/sim_R.bin", kOut), SR, SizeSR()) || !ReadArr(Form("%s/sim_W.bin", kOut), SW, SizeSW()) ||
        !ReadArr(Form("%s/sim_WV.bin", kOut), SWV, SizeSW()) || !ReadArr(Form("%s/sim_WR.bin", kOut), SWR, SizeSWR()) ||
        !ReadArr(Form("%s/sim_hq.bin", kOut), HQ, SizeHq())) return;
    InitSets();
    for (int k = 0; k < 3; ++k) CF[2][k] = 0.5 * (CF[0][k] + CF[1][k]);
    std::vector<RowB> rows;
    if (!ReadRows(rows)) return;

    // 1. nominal
    Var nom; nom.tag = "nominal";
    Res R0; Evaluate(nom, R0);
    double q0[kNQ]; Quant(MakeInp(R0), q0);
    PO("1. nominal: E_data %.4f +- %.4f, E_sim %.4f, f %.5f; f_rec (C2 x 0.879 / 0.882) %.5f  [bedge: 1141.8584, 2915.6983, 0.39162, 0.4184]\n",
       R0.c[0].Ed, std::sqrt(R0.c[0].vEd), R0.c[0].Es, q0[0], q0[5]);
    for (int k = 0; k < 3; ++k)
        PO("   class %d: Ed %.4f Es %.4f lin T12 %.4f C2 %.4f C3 %.4f\n", k, R0.c[k].Ed, R0.c[k].Es, R0.c[k].lin[0], R0.c[k].lin[1], R0.c[k].lin[2]);
    PO("KEY nom Ed %.6f %.6f %.6f Es %.6f %.6f %.6f linC2 %.6f %.6f %.6f linC3 %.6f %.6f %.6f linT12 %.6f %.6f %.6f f %.6f frec %.6f\n",
       R0.c[0].Ed, R0.c[1].Ed, R0.c[2].Ed, R0.c[0].Es, R0.c[1].Es, R0.c[2].Es, R0.c[0].lin[1], R0.c[1].lin[1], R0.c[2].lin[1],
       R0.c[0].lin[2], R0.c[1].lin[2], R0.c[2].lin[2], R0.c[0].lin[0], R0.c[1].lin[0], R0.c[2].lin[0], q0[0], q0[5]);

    // 2. the C2 template rebuilt from the rows
    const size_t n1 = (size_t)kNCls * kNPP * kNM4;
    // qv < 0: isotropic; 0-4: the QED weight variants; 5: constant class factors 0.879 / 0.882 (bedge "rec");
    // 6: constant 0.8855; 7: constant class factors of the bedge edge toy 0.8983 / 0.8628
    auto fill = [&](int qv) {
        std::vector<double> h(n1, 0.);
        for (const RowB& r : rows) {
            if (r.s0 != 1 || r.p < 0 || r.ib < 0 || r.e < 0) continue;
            double q = 1.0;
            if (qv >= 0 && qv < 5) q = r.q[qv];
            else if (qv == 5) q = r.e == 0 ? CF[2][1] : CF[2][2];
            else if (qv == 6) q = CF[2][0];
            else if (qv == 7) q = r.e == 0 ? 0.8983 : 0.8628;
            h[((size_t)r.e * kNPP + r.p) * kNM4 + r.ib] += r.W * q;
        }
        return h;
    };
    std::vector<double> iso = fill(-1);
    std::vector<double> orig(n1);
    for (size_t j = 0; j < n1; ++j) orig[j] = SS[IdxSS(sNom, gC2, 0, 0, 0) + j];
    double maxd = 0, sumA = 0, sumB = 0;
    for (size_t j = 0; j < n1; ++j) { maxd = std::max(maxd, std::fabs(iso[j] - orig[j])); sumA += iso[j]; sumB += orig[j]; }
    PO("2. C2 template from c2recon_rows.txt vs sim_S.bin [nom][C2]: max |diff| %.3e, sums %.6f / %.6f\n", maxd, sumA, sumB);

    // 3. QED-weighted templates through the full evaluation
    const char* qn[8] = {"A", "B", "Ac", "flat", "A_kloe", "const_rec", "const_8855", "const_edgetoy"};
    for (int qv = 0; qv < 8; ++qv) {
        std::vector<double> hq = fill(qv);
        for (size_t j = 0; j < n1; ++j) SS[IdxSS(sNom, gC2, 0, 0, 0) + j] = hq[j];
        Res R; Evaluate(nom, R);
        const double c1 = R.c[1].lin[1] / R0.c[1].lin[1], c2 = R.c[2].lin[1] / R0.c[2].lin[1], ca = R.c[0].lin[1] / R0.c[0].lin[1];
        const double ce = (R.c[1].lin[1] + R.c[2].lin[1]) / (R0.c[1].lin[1] + R0.c[2].lin[1]);
        // bedge formula with the class factors c1, c2 (as fc_rec does)
        const double red = ((1 - c1) * R0.c[1].lin[1] + (1 - c2) * R0.c[2].lin[1]) / (R0.c[1].lin[0] + R0.c[2].lin[0]);
        const double fLin = R0.c[0].Ed / (R0.c[0].Es * (1 - red));
        const double fDir = R.c[0].Ed / R.c[0].Es;
        PO("3. weights %-13s: factor (linear C2 step) all %.4f, E1 %.4f, E2 %.4f, E1+E2 %.4f; f = E_data / E_sim(QED, bedge own weights) %.5f; "
           "bedge formula with these factors %.5f; E_sim %.4f (E1 %.4f, E2 %.4f)\n", qn[qv], ca, c1, c2, ce, fDir, fLin, R.c[0].Es, R.c[1].Es, R.c[2].Es);
        PO("KEY qed %s cAll %.6f cE1 %.6f cE2 %.6f cE12 %.6f fDir %.6f fLin %.6f Es %.6f %.6f %.6f linC2 %.6f %.6f %.6f fE1 %.6f fE2 %.6f\n", qn[qv], ca, c1, c2, ce, fDir, fLin,
           R.c[0].Es, R.c[1].Es, R.c[2].Es, R.c[0].lin[1], R.c[1].lin[1], R.c[2].lin[1], R.c[1].Ed / R.c[1].Es, R.c[2].Ed / R.c[2].Es);
    }
    for (size_t j = 0; j < n1; ++j) SS[IdxSS(sNom, gC2, 0, 0, 0) + j] = orig[j];

    // 4. linear functionals of the slice x class fits (data covariance), M4 bin by bin
    {
        FILE* fl = fopen(Form("%s/c2recon_bedge_lin.txt", kDir), "w");
        fprintf(fl, "# class k (0 all, 1 E1, 2 E2) slice i pa pb | coefficient of each M4 bin ib (0..%d) in h (data covariance)\n", kNM4 - 1);
        const std::vector<ESlice>& set = gSets[0];
        double chk[3] = {0, 0, 0};
        for (int k = 0; k < 3; ++k) {
            for (size_t i = 0; i < set.size(); ++i) {
                const int pa = PB(set[i].lo), pb = PB(set[i].hi);
                Tmpl* tS = nullptr;
                std::vector<double> a(kNM4, 0.);
                std::vector<double> y(kNM4, 0.);
                for (int b = 0; b < kNM4; ++b) {
                    std::fill(y.begin(), y.end(), 0.); y[b] = 1.0;
                    const scan::StepFit fl1 = FitStep(SimSpec(y.data(), kNM4, scan::kMlo, kBwD), R0.c[k].cS[i], tS, &R0.c[k].Vd[i]);
                    a[b] = fl1.ok ? fl1.h : 0.;
                }
                // zero spectrum offset (should be 0) and check against the C2 step of this slice
                std::fill(y.begin(), y.end(), 0.);
                const scan::StepFit f0 = FitStep(SimSpec(y.data(), kNM4, scan::kMlo, kBwD), R0.c[k].cS[i], tS, &R0.c[k].Vd[i]);
                const std::vector<double> yc2 = SY(sNom, kMC2, k - 1, pa, pb);
                double h = 0; for (int b = 0; b < kNM4; ++b) h += a[b] * yc2[b];
                chk[k] += h;
                fprintf(fl, "%d %zu %d %d", k, i, pa, pb);
                for (int b = 0; b < kNM4; ++b) fprintf(fl, " %.10g", a[b] - (f0.ok ? f0.h : 0.));
                fprintf(fl, "\n");
            }
        }
        fclose(fl);
        PO("4. linear functionals written; sum_b a_b C2_b per class: %.4f %.4f %.4f (lin C2: %.4f %.4f %.4f)\n", chk[0], chk[1], chk[2],
           R0.c[0].lin[1], R0.c[1].lin[1], R0.c[2].lin[1]);
    }
    fclose(gO);
}
}  // namespace c2rb

void c2recon_bedge() { c2rb::Run(); }
