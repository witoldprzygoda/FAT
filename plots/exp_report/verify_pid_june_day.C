// verify_pid_june_day.C -- INDEPENDENT VERIFICATION (PID chapter, Table pid_samples, row "L: all slots"): exact polygon
// test of the four slots of PipPimEpEm_ID in ONE day file of the June 2026 lepton production (no PID cut).
//
// What it does
//   All rows with eVertReco_z > -500 (all isBest values = no PID cut). For each slot (pip, pim, ep, em) it counts the
//   track entries and how many lie inside the slot's own polygon at the PAT level (p, pi: 3 sigma; e: 5 sigma) and at
//   1 sigma, with TCutG::IsInside on (<slot>_p, <slot>_beta_new). Prints one line per file:
//   "JUNE <day> rows <n> | <slot> <entries> <inside PAT> <inside 1sig> ...". Sum the lines over the 29 days to compare
//   with the chapter's 29-day numbers (bin centres of the banana histogram).
// Inputs : one file /hdd1/HADES/pp45/exp/gen4/default/lep_<day>.root (read-only); data polygons (vae25 copy).
// Output : stdout (redirect to a log file).
// Run (all 29 days, 6 in parallel):
//   ls /hdd1/HADES/pp45/exp/gen4/default/lep_0*.root | xargs -P 6 -I{} root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/aclic_pid_verify", kTRUE)' \
//        '/home/przygoda/HADES/FAT/FAT/plots/exp_report/verify_pid_june_day.C+("{}")' | grep ^JUNE
//   (compile once before the parallel run so that the processes do not build the library at the same time)
#include <TFile.h>
#include <TTree.h>
#include <TCutG.h>
#include <TString.h>
#include <TSystem.h>
#include <cstdio>

void verify_pid_june_day(const char* file)
{
   TFile fc("/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/remote/pid_cuts/pp45_pid_cuts_exp_betanew_v1.root", "READ");
   const char* sl[4] = {"pip", "pim", "ep", "em"};
   TCutG *cP[4], *c1[4];
   for (int k = 0; k < 4; ++k) {
      cP[k] = (TCutG*)fc.Get(Form("%s_cut_%s", sl[k], k < 2 ? "3sig" : "5sig"));
      c1[k] = (TCutG*)fc.Get(Form("%s_cut_1sig", sl[k]));
   }
   TFile* f = TFile::Open(file, "READ");
   TTree* t = (TTree*)f->Get("PipPimEpEm_ID");
   t->SetBranchStatus("*", 0);
   float vz, p[4], b[4];
   t->SetBranchStatus("eVertReco_z", 1); t->SetBranchAddress("eVertReco_z", &vz);
   for (int k = 0; k < 4; ++k) {
      t->SetBranchStatus(Form("%s_p", sl[k]), 1); t->SetBranchAddress(Form("%s_p", sl[k]), &p[k]);
      t->SetBranchStatus(Form("%s_beta_new", sl[k]), 1); t->SetBranchAddress(Form("%s_beta_new", sl[k]), &b[k]);
   }
   double n = 0, inP[4] = {0, 0, 0, 0}, in1[4] = {0, 0, 0, 0};
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (!(vz > -500.f)) continue;
      n += 1;
      for (int k = 0; k < 4; ++k) { inP[k] += cP[k]->IsInside(p[k], b[k]); in1[k] += c1[k]->IsInside(p[k], b[k]); }
   }
   TString day = TString(gSystem->BaseName(file))(4, 3);
   printf("JUNE %s rows %.0f", day.Data(), n);
   for (int k = 0; k < 4; ++k) printf(" | %s %.0f %.0f", sl[k], inP[k], in1[k]);
   printf("\n");
}
