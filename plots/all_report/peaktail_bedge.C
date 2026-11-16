// ============================================================================
// peaktail_bedge.C  --  package peaktail, task 3: the low-mass tail of the eta peak propagated into the per-slice edge
//                       f_eta(edge) of package bedge, with bedge's OWN code (bedge_ana.C included unchanged)
// ============================================================================
// Written 2026-11-15 (prefix peaktail). bedge fits every M(pi+pi-) slice in M4 with pol2 + a free Gaussian eta term (no tail)
// + step, in the data and (same fit) in the simulated T12 spectrum; f = E_data / E_sim (C2 corrected: Quant, "fc_rec").
// The simulated peak (C4, C5) enters E_sim with its SMASH tail. Here the C4 and C5 spectra of the nominal configuration
// (sim_S.bin, configuration 0) are replaced in memory by
//     core x c(t) + t x tail   (c(t) = (N - t N_tail)/N_core per class over 0.30-0.40; rows kept)
// with the core/tail split of peaktail_prop.py (bedge_C45_core.bin, bedge_C45_tail.bin: tail = simulated rows whose leptons'
// momentum deficit shifts M4 by more than 5 MeV/c2), or by the QED radiation toy of C5 (bedge_C45_rad.bin), and the
// nominal bedge evaluation (Evaluate, Quant) is rerun. Nothing is written into the bedge directory.
// Outputs: /hdd3/przygoda/HADES/FAT_out/all_report/peaktail/peaktail_bedge.txt (log and machine-readable lines "PT ...")
// Run (from /home/przygoda/HADES/FAT/FAT; about 2 GB, a minute):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/peaktail/aclic", kTRUE)' \
//        'plots/all_report/peaktail_bedge.C+O'
// ============================================================================
#include "bedge_ana.C"

namespace peaktailb {
using namespace bedgea;
const char* const kOutP = "/hdd3/przygoda/HADES/FAT_out/all_report/peaktail";
FILE* gO = nullptr;
#define PO(...) do { printf(__VA_ARGS__); if (gO) fprintf(gO, __VA_ARGS__); fflush(stdout); } while (0)
const int kNA = 2 * kNCls * kNPP * kNM4;
inline size_t IdxA(int k, int e, int p, int i) { return (((size_t)k * kNCls + e) * kNPP + p) * kNM4 + i; }
std::vector<double> AC, AT, AR, SS0;
double NN[2], NT[2];
void SetPeak(double t4, double t5, int useRad, double s4 = 1.0, double s5 = 1.0) {
    const int gg[2] = {gC4, gC5};
    const double tt[2] = {t4, t5}, ss[2] = {s4, s5};
    for (int k = 0; k < 2; ++k) {
        const double c = (NN[k] - tt[k] * NT[k]) / (NN[k] - NT[k]);
        for (int e = 0; e < kNCls; ++e)
            for (int p = 0; p < kNPP; ++p)
                for (int i = 0; i < kNM4; ++i) {
                    double v;
                    if (useRad && k == 1) v = AR[IdxA(k, e, p, i)];
                    else v = AC[IdxA(k, e, p, i)] * c + tt[k] * AT[IdxA(k, e, p, i)];
                    SS[IdxSS(sNom, gg[k], e, p, i)] = ss[k] * v;
                }
    }
}
void Restore() { for (int g : {gC4, gC5}) for (int e = 0; e < kNCls; ++e) for (int p = 0; p < kNPP; ++p) for (int i = 0; i < kNM4; ++i) SS[IdxSS(sNom, g, e, p, i)] = SS0[IdxSS(sNom, g, e, p, i)]; }
void Report(const char* tag, const char* text, const Var& v, double qref[kNQ]) {
    Res R; Evaluate(v, R);
    double q[kNQ]; Quant(MakeInp(R), q);
    // bump-only convention (h_sim from C2 + C3 alone, C2 corrected): E_data / (c2 lin_C2 + lin_C3)
    const double c2 = (CF[2][1] * R.c[1].lin[1] + CF[2][2] * R.c[2].lin[1]) / (R.c[1].lin[1] + R.c[2].lin[1]);
    const double fb = R.c[0].Ed / (c2 * R.c[0].lin[1] + R.c[0].lin[2]);
    const double rest = R.c[0].lin[0] - R.c[0].lin[1] - R.c[0].lin[2];
    PO("   %-12s E_sim %8.2f (T12 lin %8.2f, rest %7.2f) f_all %.5f  fc_rec %.5f (%+6.2f %%)  E1 %.5f E2 %.5f  bump-only %.5f   %s\n", tag, R.c[0].Es,
       R.c[0].lin[0], rest, q[0], q[5], qref ? 100 * (q[5] / qref[5] - 1) : 0.0, q[10], q[11], fb, text);
    if (gO) fprintf(gO, "PT %s Es %.6f lin %.6f rest %.6f f_all %.6f fc_rec %.6f E1 %.6f E2 %.6f bumponly %.6f\n", tag, R.c[0].Es, R.c[0].lin[0], rest, q[0], q[5], q[10], q[11], fb);
    if (!qref) return;
}
}  // namespace peaktailb

void peaktail_bedge() {
    using namespace peaktailb;
    gO = fopen(Form("%s/peaktail_bedge.txt", kOutP), "w");
    gL = fopen(Form("%s/peaktail_bedge_bedgelog.txt", kOutP), "w");   // bedge's own log goes into my directory
    PO("# peaktail_bedge.C: the peak tail in the per-slice edge (bedge code)\n");
    if (!ReadArr(Form("%s/data_N.bin", kOut), DN, SizeDN()) || !ReadArr(Form("%s/sim_S.bin", kOut), SS, SizeSS()) ||
        !ReadArr(Form("%s/sim_hq.bin", kOut), HQ, SizeHq())) return;
    if (!ReadArr(Form("%s/bedge_C45_core.bin", kOutP), AC, kNA) || !ReadArr(Form("%s/bedge_C45_tail.bin", kOutP), AT, kNA) ||
        !ReadArr(Form("%s/bedge_C45_rad.bin", kOutP), AR, kNA)) return;
    SS0 = SS;
    InitSets();
    for (int k = 0; k < 3; ++k) CF[2][k] = 0.5 * (CF[0][k] + CF[1][k]);
    // normalisation per class over 0.30-0.40 (bins 6..25) and the check of core + tail against sim_S
    double maxd = 0;
    const int gg[2] = {gC4, gC5};
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
    PO("   core + tail vs bedge sim_S (C4, C5): max |diff| = %.3e; tail share 0.30-0.40: C4 %.4f of %.1f, C5 %.4f of %.1f\n", maxd, NT[0] / NN[0], NN[0], NT[1] / NN[1], NN[1]);
    Var nom; nom.tag = "nominal";
    double q0[kNQ];
    { Res R; Evaluate(nom, R); Quant(MakeInp(R), q0); }
    PO("   nominal reproduced: f_all %.5f (WP6 0.39162), fc_rec %.5f (bedge 0.418)\n", q0[0], q0[5]);
    // read t_data, s_data, t_qed from the check file of peaktail_prop.py
    double td = 0.72, sdv = 0.38, tq = 1.19;
    {
        FILE* fc = fopen(Form("%s/bedge_C45_check.txt", kOutP), "r");
        char line[512];
        while (fc && fgets(line, sizeof line, fc)) { if (!strncmp(line, "t_data", 6)) sscanf(line, "t_data %lf s_data %lf t_qed_C5 %lf", &td, &sdv, &tq); }
        if (fc) fclose(fc);
    }
    PO("   t_data %.3f +- %.3f, t_qed(C5) %.3f\n", td, sdv, tq);
    PO("\n   tag          E_sim (WP6 weights; T12 with the data covariance; its rest = peak + backgrounds) f_all fc_rec (change) E1 E2 bump-only\n");
    struct C_ { const char* tag; double t4, t5; int rad; double s4, s5; const char* text; };
    char b1[64], b2[64], b3[64], b4[64];
    snprintf(b1, 64, "t = t_data = %.2f", td); snprintf(b2, 64, "t = %.2f", td + sdv); snprintf(b3, 64, "t = %.2f", td - sdv); snprintf(b4, 64, "QED tail of C5, t_C5 = %.2f", tq);
    const C_ cs[] = {{"nominal", 1, 1, 0, 1, 1, "SMASH tail"}, {"t_data", td, td, 0, 1, 1, b1}, {"t_data+1s", td + sdv, td + sdv, 0, 1, 1, b2},
                     {"t_data-1s", td - sdv, td - sdv, 0, 1, 1, b3}, {"t0", 0, 0, 0, 1, 1, "no tail (t = 0)"}, {"t2", 2, 2, 0, 1, 1, "doubled tail"},
                     {"qed", 1, tq, 0, 1, 1, b4}, {"qed_rows", 1, 1, 1, 1, 1, "QED toy rows of C5"},
                     {"comp", 1, 1, 0, 0.237 / 0.418249, 0.497 / 0.418249, "peak in E_sim with the Part III composition (C4 x 0.237/a_P, C5 x 0.497/a_P)"},
                     {"comp_t0", 0, 0, 0, 0.237 / 0.418249, 0.497 / 0.418249, "composition, no tail"},
                     {"nopeak", 0, 0, 0, 1e-12, 1e-12, "no simulated peak in E_sim at all"}};
    for (const C_& c : cs) {
        SetPeak(c.t4, c.t5, c.rad, c.s4, c.s5);
        Report(c.tag, c.text, nom, q0);
        Restore();
    }
    // the same with the simulation fitted with the final data covariance (bedge variant wgtData)
    Var vw = nom; vw.wgt = 1;
    double qw[kNQ]; { Res R; Evaluate(vw, R); Quant(MakeInp(R), qw); }
    PO("\n   with the simulation fitted with the data covariance (variant wgtData): nominal fc_rec %.5f\n", qw[5]);
    for (const C_& c : cs) {
        if (strcmp(c.tag, "t0") && strcmp(c.tag, "t2") && strcmp(c.tag, "t_data")) continue;
        SetPeak(c.t4, c.t5, c.rad, c.s4, c.s5);
        Report(c.tag, c.text, vw, qw);
        Restore();
    }
    fclose(gO); fclose(gL);
}
