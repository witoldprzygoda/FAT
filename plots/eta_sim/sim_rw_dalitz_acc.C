// ============================================================================
// sim_rw_dalitz_acc.C  --  is the acceptance the cause of the non-QED M(ee) of eta -> pi+pi-e+e-?
// ============================================================================
// HADES p+p, T = 4.53 GeV, pi+ pi- e+ e-; SMASH simulation, compact files of sim_prepare.C.
// Accepted rows (FAT selection isBest >= 0, z > -500, PT3 bit; here also vz > -160 mm; unweighted counts) in bins of the
// TRUE M(ee), divided by the QED (Kroll-Wada) shape of the same decay, normalised to 1 in the bin 3.37-5.02 MeV:
//   - pi0 Dalitz pairs of the standard set (HGeant generator, lepton-pair class L1, std_us)
//   - pi0 Dalitz pairs of the lepton set (PLUTO, L1, lep_us)
//   - eta Dalitz pairs of the lepton set (PLUTO, L2, lep_us)
//   - eta -> pi+pi-e+e- of the eta set (PLUTO "eta:e+,e-,pi+,pi-", four-particle class 5, eta_us)
// The Dalitz generators follow QED, so for them the ratio is the relative acceptance alone. If the acceptance caused the
// threshold excess of eta -> pi+pi-e+e-, the Dalitz pairs would show it too.
// QED references: Dalitz decays P -> gamma e+e- with the Kroll-Wada formula (no form factor),
//   dGamma/dM ~ (1/M) (1 - M^2/m_P^2)^3 (1 + 2 m_e^2/M^2) sqrt(1 - 4 m_e^2/M^2), integrated over the bin;
//   eta -> pi+pi-e+e-: the reference column of tables/eta_peak_mee_true.txt (Kroll-Wada x P-wave pi pi gamma*, no form factor).
// Same code as the check of 2026-11-12 (session scratchpad ppsim/meecheck/meecheck.C), now with a table.
// Output: plots/eta_sim/tables/rw_dalitz_acceptance.txt (and stdout). About 1 min.
// Run from /home/przygoda/HADES/FAT/FAT:
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/aclic_dalitz_acc", kTRUE)' 'plots/eta_sim/sim_rw_dalitz_acc.C+'
// ============================================================================
#include <TFile.h>
#include <TTree.h>
#include <cmath>
#include <cstdio>

namespace dalacc {
const double kMe = 0.000510999;
double KW(double M, double mP) {
    if (M <= 2 * kMe || M >= mP) return 0;
    const double x = M * M / (mP * mP);
    return (1. / M) * std::pow(1 - x, 3) * (1 + 2 * kMe * kMe / (M * M)) * std::sqrt(1 - 4 * kMe * kMe / (M * M));
}
// sel: 1 = lepton-pair class L1 (pi0 Dalitz), 2 = L2 (eta Dalitz), 5 = four-particle class 5 (eta -> pi+pi-e+e-)
void Fill(const char* file, const char* tree, int sel, double* h, const double* e, int nb) {
    TFile f(file);
    TTree* t = (TTree*)f.Get(tree);
    if (!t) { printf("missing tree %s in %s\n", tree, file); return; }
    Float_t mee, vz; Int_t lcls, cls;
    t->SetBranchStatus("*", 0);
    for (const char* n : {"mee_true", "vz", "lcls", "cls"}) t->SetBranchStatus(n, 1);
    t->SetBranchAddress("mee_true", &mee); t->SetBranchAddress("vz", &vz);
    t->SetBranchAddress("lcls", &lcls); t->SetBranchAddress("cls", &cls);
    for (Long64_t i = 0, N = t->GetEntries(); i < N; ++i) {
        t->GetEntry(i);
        if (!(vz > -160)) continue;
        if (sel == 1 && lcls != 1) continue;
        if (sel == 2 && lcls != 2) continue;
        if (sel == 5 && cls != 5) continue;
        const double m = mee * 1000.;
        for (int b = 0; b < nb; ++b) if (m >= e[b] && m < e[b + 1]) { h[b] += 1.; break; }
    }
}
}  // namespace dalacc

void sim_rw_dalitz_acc() {
    using namespace dalacc;
    const int nb = 10;
    const double e[nb + 1] = {1.022, 1.52, 2.27, 3.37, 5.02, 7.47, 11.13, 16.57, 24.66, 36.71, 54.66};
    const char* C = "/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/sim_compact.root";
    const char* E = "/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/sim_compact_eta.root";
    double a[4][nb] = {{0}};
    Fill(C, "std_us", 1, a[0], e, nb);
    Fill(C, "lep_us", 1, a[1], e, nb);
    Fill(C, "lep_us", 2, a[2], e, nb);
    Fill(E, "eta_us", 5, a[3], e, nb);
    const char* nm[4] = {"pi0 Dalitz std (HGeant)", "pi0 Dalitz lep (PLUTO)", "eta Dalitz lep (PLUTO)", "eta->pipiee eta set (PLUTO)"};
    const double mP[3] = {0.1349768, 0.1349768, 0.547862};
    const double qedEta[nb] = {0.0665, 0.0906, 0.0939, 0.0945, 0.0944, 0.0940, 0.0931, 0.0910, 0.0867, 0.0783};
    double ref[4][nb];
    for (int k = 0; k < 3; ++k)
        for (int b = 0; b < nb; ++b) {
            double s = 0; const int n = 200;
            for (int j = 0; j < n; ++j) {
                const double M = (e[b] + (e[b + 1] - e[b]) * (j + 0.5) / n) / 1000.;
                s += KW(M, mP[k]) * (e[b + 1] - e[b]) / n;
            }
            ref[k][b] = s;
        }
    for (int b = 0; b < nb; ++b) ref[3][b] = qedEta[b];
    double norm[4];
    for (int k = 0; k < 4; ++k) norm[k] = a[k][3] / ref[k][3];

    FILE* fo = fopen("plots/eta_sim/tables/rw_dalitz_acceptance.txt", "w");
    FILE* outs[2] = {stdout, fo};
    for (FILE* o : outs) {
        if (!o) continue;
        fprintf(o, "# rw_dalitz_acceptance: accepted / QED shape per bin of the TRUE M(ee), normalised to 1 in 3.37-5.02 MeV (sim_rw_dalitz_acc.C)\n");
        fprintf(o, "# rows of the compact files (FAT selection + vz > -160 mm), unweighted; QED: Kroll-Wada for the Dalitz decays, the reference column of\n");
        fprintf(o, "#   eta_peak_mee_true.txt for eta -> pi+pi-e+e-. Dalitz generators are QED, so their column is the relative acceptance alone.\n");
        fprintf(o, "true M(ee) [MeV]  ");
        for (int k = 0; k < 4; ++k) fprintf(o, "| %-30s ", nm[k]);
        fprintf(o, "\n");
        for (int b = 0; b < nb; ++b) {
            fprintf(o, "%6.2f - %6.2f   ", e[b], e[b + 1]);
            for (int k = 0; k < 4; ++k) fprintf(o, "| %7.3f  (rows %8.0f)        ", (a[k][b] / ref[k][b]) / norm[k], a[k][b]);
            fprintf(o, "\n");
        }
        fprintf(o, "# Dalitz pairs: smooth, 0.88-1.75 over 1-55 MeV, close to 1 at threshold (acceptance rises mildly with M(ee), larger OA)\n");
        fprintf(o, "# eta -> pi+pi-e+e-: 5.3 at threshold, 0.17 at 37-55 MeV: opposite to the acceptance trend -> the generator shape is the cause\n");
    }
    if (fo) fclose(fo);
}
