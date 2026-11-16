// ============================================================================
// peaktail_bedge_v3.C  --  rulings D38 + D40 (15 November 2026): the slope of the per-slice edge f_eta(edge) with the size t
//                          of the low-mass tail of the simulated eta peak, evaluated with the D38 class factors of the
//                          per-row QED weights (E1 0.868178, E2 0.853877, from c2recon_bedge.C) instead of the WP8 means.
// ============================================================================
// Written 2026-11-15 (numbers agent of the fourth integration). New file; bedge_ana.C is #included unchanged (namespace
// bedgea: Evaluate, Quant); the C4 / C5 core and tail arrays of package peaktail (bedge_C45_core.bin, bedge_C45_tail.bin) are
// read, the tail is scaled exactly as in peaktail_bedge.C (core x c(t) + t x tail, rows kept per class over 0.30-0.40).
// Nothing is written into the bedge or peaktail directories.
// Output: /hdd3/przygoda/HADES/FAT_out/all_report/bedge_v3/peaktail_bedge_v3.txt ("PT3" lines)
// Run (from /home/przygoda/HADES/FAT/FAT; about 2 GB, under a minute):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/bedge_v3/aclic", kTRUE)' \
//        'plots/all_report/peaktail_bedge_v3.C+O'
// ============================================================================
#include "bedge_ana.C"

namespace ptb3 {
using namespace bedgea;
const char* const kOutP = "/hdd3/przygoda/HADES/FAT_out/all_report/peaktail";
const char* const kOut3 = "/hdd3/przygoda/HADES/FAT_out/all_report/bedge_v3";
FILE* gO = nullptr;
#define PO(...) do { printf(__VA_ARGS__); if (gO) fprintf(gO, __VA_ARGS__); fflush(stdout); } while (0)
const int kNA = 2 * kNCls * kNPP * kNM4;
inline size_t IdxA(int k, int e, int p, int i) { return (((size_t)k * kNCls + e) * kNPP + p) * kNM4 + i; }
std::vector<double> AC, AT, SS0;
double NN[2], NT[2];
void SetPeak(double t) {
    const int gg[2] = {gC4, gC5};
    for (int k = 0; k < 2; ++k) {
        const double c = (NN[k] - t * NT[k]) / (NN[k] - NT[k]);
        for (int e = 0; e < kNCls; ++e)
            for (int p = 0; p < kNPP; ++p)
                for (int i = 0; i < kNM4; ++i) SS[IdxSS(sNom, gg[k], e, p, i)] = AC[IdxA(k, e, p, i)] * c + t * AT[IdxA(k, e, p, i)];
    }
}
void Restore() { for (int g : {gC4, gC5}) for (int e = 0; e < kNCls; ++e) for (int p = 0; p < kNPP; ++p) for (int i = 0; i < kNM4; ++i) SS[IdxSS(sNom, g, e, p, i)] = SS0[IdxSS(sNom, g, e, p, i)]; }
}  // namespace ptb3

void peaktail_bedge_v3() {
    using namespace ptb3;
    gO = fopen(Form("%s/peaktail_bedge_v3.txt", kOut3), "w");
    gL = fopen(Form("%s/peaktail_bedge_v3_bedgelog.txt", kOut3), "w");
    PO("# peaktail_bedge_v3.C: per-slice edge against the peak-tail size t, C2 factors of D38 and of D35\n");
    if (!ReadArr(Form("%s/data_N.bin", kOut), DN, SizeDN()) || !ReadArr(Form("%s/sim_S.bin", kOut), SS, SizeSS()) ||
        !ReadArr(Form("%s/sim_hq.bin", kOut), HQ, SizeHq())) return;
    if (!ReadArr(Form("%s/bedge_C45_core.bin", kOutP), AC, kNA) || !ReadArr(Form("%s/bedge_C45_tail.bin", kOutP), AT, kNA)) return;
    SS0 = SS;
    InitSets();
    const int gg[2] = {gC4, gC5};
    double maxd = 0;
    for (int k = 0; k < 2; ++k) {
        NN[k] = NT[k] = 0;
        for (int e = 0; e < kNCls; ++e)
            for (int p = 0; p < kNPP; ++p)
                for (int i = 0; i < kNM4; ++i) {
                    const double c = AC[IdxA(k, e, p, i)], t = AT[IdxA(k, e, p, i)];
                    maxd = std::max(maxd, std::fabs(c + t - SS[IdxSS(sNom, gg[k], e, p, i)]));
                    if (p >= 6 && p < 26) { NN[k] += c + t; NT[k] += t; }
                }
    }
    PO("   core + tail vs sim_S: max |diff| %.3e; tail share C4 %.4f, C5 %.4f\n", maxd, NT[0] / NN[0], NT[1] / NN[1]);
    const double tv[] = {0.0, 0.5, 1.0, 1.5, 2.0};
    const char* lab[2] = {"D35", "D38"};
    for (int mode = 0; mode < 2; ++mode) {
        if (mode == 0) { for (int k = 0; k < 3; ++k) CF[2][k] = 0.5 * (CF[0][k] + CF[1][k]); }
        else { CF[2][0] = 0.859863; CF[2][1] = 0.868178; CF[2][2] = 0.853877; }
        Var nom; nom.tag = "nominal";
        double q1[kNQ];
        SetPeak(1.0); { Res R; Evaluate(nom, R); Quant(MakeInp(R), q1); } Restore();
        PO("   %s (CF E1 %.6f, E2 %.6f): t = 1: f_all %.6f fc_rec %.6f a_k48 %.6f a_k48_unc %.6f\n", lab[mode], CF[2][1], CF[2][2], q1[0], q1[5], q1[12], q1[13]);
        for (double t : tv) {
            SetPeak(t);
            Res R; Evaluate(nom, R);
            double q[kNQ]; Quant(MakeInp(R), q);
            Restore();
            PO("PT3 %s t %.2f Es %.4f f_all %.6f (%+.4f) fc_rec %.6f (%+.4f) a_k48 %.6f (%+.4f) a_k48_unc %.6f (%+.4f)\n", lab[mode], t, R.c[0].Es,
               q[0], q[0] / q1[0] - 1, q[5], q[5] / q1[5] - 1, q[12], q[12] / q1[12] - 1, q[13], q[13] / q1[13] - 1);
        }
    }
    fclose(gO); fclose(gL);
}
