// pid_sim_startcheck.C -- PID chapter of the experimental-data report: cross-check that the SIMULATION has no START
// problem in the DST track sorter (the only use of simulation in this chapter).
//
// What it does
//   Reads (read-only) the simulated PAT outputs of the betanew simulation production, tree PipPimEpEm_ID, rows with
//   isBest >= 0 and eVertReco_z > -500 (all trigger words), and counts for the e+ and e- slots: DST beta < 0.9
//   (lepton pass of HParticleTrackSorter failed), kIsLepton == 0, the minimum DST beta and tof, and the share of rows
//   with evtStart == 1. The same numbers for the data are in tables_exp/pid_dstbeta.txt.
// Inputs : /hdd3/przygoda/HADES/files/sim/ppee_std_all.root, /hdd3/przygoda/HADES/files/sim/ppee_lep_all.root
//          (the inputs of the simulation branch, lists sim_std.list / sim_lep.list of FAT_simwt)
// Output : /home/przygoda/HADES/FAT/FAT/docs/latex/tables_exp/pid_sim_startcheck.txt
// Run:
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/aclic_pid", kTRUE)' \
//        '/home/przygoda/HADES/FAT/FAT/plots/exp_report/pid_sim_startcheck.C+'
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <cstdio>
#include <algorithm>

void pid_sim_startcheck()
{
   const char* files[2] = {"/hdd3/przygoda/HADES/files/sim/ppee_std_all.root", "/hdd3/przygoda/HADES/files/sim/ppee_lep_all.root"};
   FILE* out = fopen("/home/przygoda/HADES/FAT/FAT/docs/latex/tables_exp/pid_sim_startcheck.txt", "w");
   fprintf(out, "# Simulation (SMASH, betanew sim production), PipPimEpEm_ID rows with isBest >= 0 && eVertReco_z > -500, all trigger words\n");
   fprintf(out, "# per lepton slot: share with DST beta < 0.9, share with kIsLepton == 0, minimum DST beta, minimum DST tof [ns]; share of rows with evtStart == 1\n");
   for (const char* fn : files) {
      TFile* f = TFile::Open(fn, "READ");
      if (!f || f->IsZombie()) { fprintf(out, "%s: cannot open\n", fn); continue; }
      TTree* t = (TTree*)f->Get("PipPimEpEm_ID");
      t->SetEstimate(t->GetEntries() + 1);
      const TString base = "isBest>=0&&eVertReco_z>-500";
      const double n = t->GetEntries(base);
      fprintf(out, "%s: rows %.0f\n", fn, n);
      for (const char* s : {"ep", "em"}) {
         const double low = t->GetEntries(base + Form("&&%s_beta<0.9", s)) / n;
         const double kl0 = t->GetEntries(base + Form("&&%s_kIsLepton==0", s)) / n;
         t->Draw(Form("%s_beta:%s_tof", s, s), base, "goff");
         double bmin = 1e9, tmin = 1e9;
         for (Long64_t i = 0; i < t->GetSelectedRows(); ++i) { bmin = std::min(bmin, t->GetV1()[i]); tmin = std::min(tmin, t->GetV2()[i]); }
         fprintf(out, "  %s: DST beta < 0.9 %.5f ; kIsLepton == 0 %.5f ; min DST beta %.3f ; min DST tof %.3f ns\n", s, low, kl0, bmin, tmin);
      }
      fprintf(out, "  evtStart == 1: %.4f of the rows\n", t->GetEntries(base + "&&evtStart==1") / n);
      f->Close();
   }
   fclose(out);
   printf("written /home/przygoda/HADES/FAT/FAT/docs/latex/tables_exp/pid_sim_startcheck.txt\n");
}
