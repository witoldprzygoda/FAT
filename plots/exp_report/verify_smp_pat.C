// ============================================================================
// verify_smp_pat.C  --  independent PAT-level counts of one day file for the
//                       verification of the smp chapter (cut-flow steps 0-2)
// ============================================================================
// HADES p+p, T = 4.53 GeV (feb22), channel pi+ pi- e+ e-. EXPERIMENTAL DATA ONLY.
// Verification macro (written by the verifier of the smp chapter). Input:
// /hdd3/przygoda/HADES/files/exp/gen4/lep_<day>.root. With TTree::GetEntries
// (selection strings, not the loop of smp_pat_counts.C) it prints for the trees
// PipPimEpEm_ID, PipPimEpEp_ID, PipPimEmEm_ID: entries, rows with isBest == 1
// (= events, one best row per event), rows with eVertReco_z > -500, rows with
// eVertReco_z > -500 && trigbit == 8192, and the entries of EpEm_ID.
// Output: one line per tree to
//   /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/verify_smp/pat_<day>.txt
// Run (from /home/przygoda/HADES/FAT/FAT), all days 6 in parallel:
//   for d in 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 59 60 61 62 63 64 65 66 67; do echo $d; done | \
//   xargs -P 6 -I{} root -l -b -q 'plots/exp_report/verify_smp_pat.C({})'
// ============================================================================
#include <cstdio>
#include "TFile.h"
#include "TTree.h"

void verify_smp_pat(int day)
{
    TFile* f = TFile::Open(Form("/hdd3/przygoda/HADES/files/exp/gen4/lep_%03d.root", day), "READ");
    FILE* o = fopen(Form("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/verify_smp/pat_%03d.txt", day), "w");
    const char* tr[3] = {"PipPimEpEm_ID", "PipPimEpEp_ID", "PipPimEmEm_ID"};
    for (int k = 0; k < 3; ++k) {
        TTree* t = (TTree*)f->Get(tr[k]);
        const Long64_t n = t->GetEntries();
        const Long64_t b1 = t->GetEntries("isBest == 1");
        const Long64_t v = t->GetEntries("eVertReco_z > -500");
        const Long64_t vb = t->GetEntries("eVertReco_z > -500 && isBest == 1");
        const Long64_t vt = t->GetEntries("eVertReco_z > -500 && trigbit == 8192");
        fprintf(o, "%03d %s %lld %lld %lld %lld %lld\n", day, tr[k], n, b1, v, vb, vt);
    }
    TTree* e = (TTree*)f->Get("EpEm_ID");
    fprintf(o, "%03d EpEm_ID %lld\n", day, e->GetEntries());
    fclose(o);
    f->Close();
}
