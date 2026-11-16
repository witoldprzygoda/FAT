// ============================================================================
// wp4_mcsplit.C  --  WP4: which template MC noise causes the bias and the spread (prefix wp4)
// ============================================================================
// Asimov data of scenario (a) (truth = full simulated sample: smoothed non-eta background, raw eta classes, data CB),
// fitted 100 times with templates in which only one part comes from a bootstrap replica of the full sample
// (hists_v1 _rep, WP0):
//   "signal"     : C2, C3, C4, C5, G_eta from the replica, smoothed non-eta background from the full sample
//   "background" : signal templates from the full sample, smoothed G_D and G_C from the replica
//   "both"       : everything from the replica (= the RF configuration of wp4_toys.C)
// Printed: mean shift (replica mean - true) and RMS over the replicas, relative to the true value.
// Also the same split for the half samples (truth half h, signal templates half s, background templates half b).
// Output: docs/latex/tables_all/wp4_mcsplit.txt
//
// Run (from /home/przygoda/HADES/FAT/FAT):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_wp4", kTRUE)' 'plots/all_report/wp4_mcsplit.C+O'
// ============================================================================
#include <cmath>
#include <cstdio>

#include "wp4_lib.h"

void wp4_mcsplit() {
    using namespace wp4;
    Grid G; G.Init();
    TFile* fh = TFile::Open(kHistFile);
    DataIn Din;
    if (!LoadData(fh, G, Din)) return;
    TFile* ft = TFile::Open(kTmplFile);
    Tmpl* Tf = new Tmpl();
    LoadTmplFull(ft, kMnom, G, *Tf);
    FitTmpl* F = new FitTmpl();
    MakeFitTmpl(*Tf, G, *F);
    TruthSmooth* TS = new TruthSmooth();
    MakeTruthSmooth(*Tf, G, *TS);
    std::vector<FitTmpl*> Fr;
    std::vector<Tmpl*> Tr;
    if (!LoadReplicas(fh, G, *Tf, 100, Fr, Tr)) return;
    Truth U;
    MakeTruth(kScen[0], *Tf, *TS, Din, G, U);
    Obs O;
    MakeToy(U, Din, G, nullptr, O);
    const int nq = 7;
    const int qs[nq] = {qFMa, qFMYb, qFMYp, qFMrho, qFGBRG, qB4Yb, qPSBYp};
    FILE* out = fopen("/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all/wp4_mcsplit.txt", "w");
    fprintf(out, "# wp4_mcsplit.txt -- plots/all_report/wp4_mcsplit.C; scenario (a), Asimov data (truth = full sample), S0, L = 6079.2 nb^-1\n");
    fprintf(out, "# full-sample templates with ONE part replaced by a bootstrap replica (100 replicas of WP0): mean shift and RMS [%% of true]\n");
    fprintf(out, "%-26s", "replaced part");
    for (int i = 0; i < nq; ++i) fprintf(out, " %20s", QName(qs[i]));
    fprintf(out, "\n");
    const char* nm[3] = {"signal (C2-C5, G_eta)", "background (G_D, G_C)", "both"};
    for (int mode = 0; mode < 3; ++mode) {
        double s[nq] = {0}, s2[nq] = {0};
        for (int r = 0; r < 100; ++r) {
            FitTmpl X = *F;
            Tmpl* T = nullptr;
            if (mode == 0) {
                T = new Tmpl(*Tf);
                for (int c : {kTC2, kTC3, kTC4, kTC5, kTGeta}) {
                    for (int k = 0; k < 10; ++k) T->F10[c][k] = Tr[r]->F10[c][k];
                    for (int j = 0; j < kNj; ++j) T->F15[c][H15(j, 2)] = Tr[r]->F15[c][H15(j, 2)];
                }
                X.T = T;
                for (int t = 0; t < kNSig; ++t) for (int w = 0; w < 4; ++w) X.W[t][w] = Fr[r]->W[t][w];
            } else if (mode == 1) {
                for (int k = 0; k < 10; ++k) { X.gd10[k] = Fr[r]->gd10[k]; X.gc10[k] = Fr[r]->gc10[k]; X.gdR[k] = Fr[r]->gdR[k]; X.gcR[k] = Fr[r]->gcR[k]; }
                for (int j = 0; j < kNj; ++j) { X.gd5[j] = Fr[r]->gd5[j]; X.gc5[j] = Fr[r]->gc5[j]; }
            } else {
                X = *Fr[r];
            }
            Result R;
            RunEstimators(O, Din, X, G, R);
            for (int i = 0; i < nq; ++i) { s[i] += R.est[qs[i]]; s2[i] += R.est[qs[i]] * R.est[qs[i]]; }
            delete T;
        }
        fprintf(out, "%-26s", nm[mode]);
        for (int i = 0; i < nq; ++i) {
            const double m = s[i] / 100, rm = std::sqrt(std::max(0., s2[i] / 100 - m * m)), t = U.val[qs[i]];
            fprintf(out, "   %+6.1f (rms %5.1f)", 100 * (m - t) / t, 100 * rm / t);
        }
        fprintf(out, "\n");
    }
    // half samples
    fprintf(out, "\n# half samples (x2): truth half / signal templates half / background templates half -> bias [%% of true]\n");
    Tmpl* T2[2]; FitTmpl* F2[2]; TruthSmooth* S2[2];
    for (int h = 0; h < 2; ++h) {
        T2[h] = new Tmpl(); LoadTmpl(ft, kMnom, h, G, *T2[h]);
        F2[h] = new FitTmpl(); MakeFitTmpl(*T2[h], G, *F2[h]);
        S2[h] = new TruthSmooth(); MakeTruthSmooth(*T2[h], G, *S2[h]);
    }
    for (int h = 0; h < 2; ++h)
        for (int sg = 0; sg < 2; ++sg)
            for (int bk = 0; bk < 2; ++bk) {
                Truth V; MakeTruth(kScen[0], *T2[h], *S2[h], Din, G, V);
                Obs Oh; MakeToy(V, Din, G, nullptr, Oh);
                FitTmpl X = *F2[bk];
                Tmpl* Tm = new Tmpl(*T2[bk]);
                for (int c : {kTC2, kTC3, kTC4, kTC5, kTGeta}) {
                    for (int k = 0; k < 10; ++k) { Tm->F10[c][k] = T2[sg]->F10[c][k]; Tm->F10v[c][k] = T2[sg]->F10v[c][k]; }
                    for (int j = 0; j < kNj; ++j) { Tm->F15[c][H15(j, 2)] = T2[sg]->F15[c][H15(j, 2)]; Tm->F15v[c][H15(j, 2)] = T2[sg]->F15v[c][H15(j, 2)]; }
                }
                X.T = Tm;
                for (int t = 0; t < kNSig; ++t) for (int w = 0; w < 4; ++w) X.W[t][w] = F2[sg]->W[t][w];
                Result R;
                RunEstimators(Oh, Din, X, G, R);
                fprintf(out, "truth %c signal %c background %c:", "AB"[h], "AB"[sg], "AB"[bk]);
                for (int i = 0; i < nq; ++i) fprintf(out, "  %s %+6.1f", QName(qs[i]), 100 * (R.est[qs[i]] - V.val[qs[i]]) / V.val[qs[i]]);
                fprintf(out, "\n");
                delete Tm;
            }
    fclose(out);
    printf("wp4_mcsplit: wrote docs/latex/tables_all/wp4_mcsplit.txt\n");
}
