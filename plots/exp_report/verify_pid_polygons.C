// verify_pid_polygons.C -- INDEPENDENT VERIFICATION of Table pid_polygons (PID chapter): centre and half width of the
// 1 sigma and PAT-level polygons of the data (betanew_v1) and simulation (betanew_v6) PID cut files.
//
// What it does
//   For each species (p, pi+, pi-, e+, e-) and p = 200, 500, 1000, 1500 MeV/c, scans beta_new from 0 to 1.3 in steps of
//   1e-5 with TCutG::IsInside (a different method from the chapter, which intersected the polygon edges with a
//   vertical line) and prints the lowest and highest beta inside, the centre and the half width; also the momentum
//   range (min, max x of the polygon points) and whether the e- polygon equals the e+ polygon point by point.
// Inputs : data cuts /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/remote/pid_cuts/pp45_pid_cuts_exp_betanew_v1.root
//          simulation cuts /hdd3/przygoda/HADES/files/sim/pp45_pid_cuts_sim_betanew_v6.root
// Output : /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/pid/verify/verify_pid_polygons.txt
// Run:
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/aclic_pid_verify", kTRUE)' \
//        '/home/przygoda/HADES/FAT/FAT/plots/exp_report/verify_pid_polygons.C+'
#include <TFile.h>
#include <TCutG.h>
#include <TMath.h>
#include <cstdio>
#include <cmath>

void verify_pid_polygons()
{
   const char* fn[2] = {"/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/remote/pid_cuts/pp45_pid_cuts_exp_betanew_v1.root",
                        "/hdd3/przygoda/HADES/files/sim/pp45_pid_cuts_sim_betanew_v6.root"};
   const char* sp[5] = {"p", "pip", "pim", "ep", "em"};
   const char* lv[2][5] = {{"1sig", "1sig", "1sig", "1sig", "1sig"}, {"3sig", "3sig", "3sig", "5sig", "5sig"}};
   const double pp[4] = {200, 500, 1000, 1500};
   FILE* fo = fopen("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/pid/verify/verify_pid_polygons.txt", "w");
   TFile* f[2] = {TFile::Open(fn[0]), TFile::Open(fn[1])};
   for (int s = 0; s < 5; ++s) {
      for (int v = 0; v < 2; ++v) {
         TCutG* c = (TCutG*)f[v]->Get(Form("%s_cut_%s", sp[s], lv[1][s]));
         const double xmin = TMath::MinElement(c->GetN(), c->GetX()), xmax = TMath::MaxElement(c->GetN(), c->GetX());
         fprintf(fo, "%-3s %s PAT-level polygon: %d points, x range %.0f - %.0f MeV/c\n", sp[s], v == 0 ? "data v1" : "sim  v6", c->GetN(), xmin, xmax);
      }
      for (int j = 0; j < 4; ++j) {
         fprintf(fo, "    p = %4.0f:", pp[j]);
         for (int l = 0; l < 2; ++l)
            for (int v = 0; v < 2; ++v) {
               TCutG* c = (TCutG*)f[v]->Get(Form("%s_cut_%s", sp[s], lv[l][s]));
               double lo = -1, hi = -1;
               for (int i = 0; i <= 130000; ++i) {
                  const double b = i * 1e-5;
                  if (c->IsInside(pp[j], b)) { if (lo < 0) lo = b; hi = b; }
               }
               fprintf(fo, "  %s %s centre %.4f half %.4f |", lv[l][s], v == 0 ? "v1" : "v6", 0.5 * (lo + hi), 0.5 * (hi - lo));
            }
         fprintf(fo, "\n");
      }
   }
   // e- = e+ ?
   for (int v = 0; v < 2; ++v) {
      int ndiff = 0;
      for (const char* l : {"1sig", "2sig", "25sig", "3sig", "35sig", "5sig"}) {
         TCutG* a = (TCutG*)f[v]->Get(Form("ep_cut_%s", l));
         TCutG* b = (TCutG*)f[v]->Get(Form("em_cut_%s", l));
         if (a->GetN() != b->GetN()) { ++ndiff; continue; }
         for (int i = 0; i < a->GetN(); ++i) if (a->GetX()[i] != b->GetX()[i] || a->GetY()[i] != b->GetY()[i]) { ++ndiff; break; }
      }
      fprintf(fo, "%s: e- polygons differing from e+ at %d of 6 levels\n", v == 0 ? "data v1" : "sim v6", ndiff);
      fprintf(fo, "%s: d_cut present: %s\n", v == 0 ? "data v1" : "sim v6", f[v]->Get("d_cut") ? "yes" : "no");
   }
   fclose(fo);
}
