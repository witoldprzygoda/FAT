// verify_pid_samples.C -- INDEPENDENT VERIFICATION of Table pid_samples (PID samples without PID cut) and of
// pid_sim_startcheck.txt (PID chapter).
//
// What it does (exact TCutG::IsInside tests on the track values, i.e. not from the binned banana histograms):
//   1. June 2026 lepton production (no PID cut), ONE day (036): EpEm_ID and PipPimEpEm_ID, all rows (all isBest values),
//      eVertReco_z > -500: share of each slot inside its own polygon at 1 sigma, 2 sigma and the PAT level
//      (p, pi: 3 sigma; e: 5 sigma) and, for the lepton slots, inside the pion polygon of the same charge.
//   2. June 2026 hadron production had_034, tree PipPim_ID, first 5 M rows, eVertReco_z > -500, tracks WITHOUT a RICH
//      ring (isring == 0): pi+ and pi- inside their own polygons and pi+ inside the proton and e+ polygons.
//      The chapter used 4 days x 40 M rows (and bin centres), so only approximate agreement is expected.
//   3. Simulation START check: ppee_std_all.root and ppee_lep_all.root (read-only), PipPimEpEm_ID, isBest >= 0 &&
//      eVertReco_z > -500: share of e+ / e- with DST beta < 0.9, minimum DST beta, share of rows with evtStart == 1.
// Inputs : /hdd1/HADES/pp45/exp/gen4/default/lep_036.root, /hdd2/przygoda/HADES/pp45/exp/had_034.root,
//          /hdd3/przygoda/HADES/files/sim/ppee_{std,lep}_all.root, data polygons (vae25 copy, md5 b8f0cc2f...)
// Output : /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/pid/verify/verify_pid_samples.txt
// Run:
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/aclic_pid_verify", kTRUE)' \
//        '/home/przygoda/HADES/FAT/FAT/plots/exp_report/verify_pid_samples.C+'
#include <TFile.h>
#include <TTree.h>
#include <TCutG.h>
#include <TString.h>
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>
#include <map>

namespace vps {
const char* kCuts = "/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/remote/pid_cuts/pp45_pid_cuts_exp_betanew_v1.root";
FILE* gOut = nullptr;
TFile* gC = nullptr;
// polygons read once and cached (TFile::Get would read a new copy at every call)
std::map<std::string, TCutG*> gCache;
TCutG* Cut(const char* sp, const char* lv)
{
   const std::string k = Form("%s_cut_%s", sp, lv);
   auto it = gCache.find(k);
   if (it != gCache.end()) return it->second;
   TCutG* c = (TCutG*)gC->Get(k.c_str());
   gCache[k] = c;
   return c;
}

// slots: name, polygon species, sign (+1/-1); requireNoRing: only tracks with isring == 0
void Sample(const char* file, const char* tree, std::vector<std::string> slots, std::vector<std::string> species, long long maxrows, bool requireNoRing)
{
   TFile* f = TFile::Open(file, "READ");
   TTree* t = (TTree*)f->Get(tree);
   t->SetBranchStatus("*", 0);
   const int ns = slots.size();
   float vz;
   std::vector<float> p(ns), b(ns), ring(ns);
   t->SetBranchStatus("eVertReco_z", 1); t->SetBranchAddress("eVertReco_z", &vz);
   for (int i = 0; i < ns; ++i) {
      for (const char* x : {"_p", "_beta_new", "_isring"}) t->SetBranchStatus((slots[i] + x).c_str(), 1);
      t->SetBranchAddress((slots[i] + "_p").c_str(), &p[i]);
      t->SetBranchAddress((slots[i] + "_beta_new").c_str(), &b[i]);
      t->SetBranchAddress((slots[i] + "_isring").c_str(), &ring[i]);
   }
   std::vector<double> n(ns, 0), i1(ns, 0), i2(ns, 0), iP(ns, 0), iPi(ns, 0), iProt(ns, 0), iE(ns, 0);
   Long64_t nent = t->GetEntries();
   if (maxrows > 0 && nent > maxrows) nent = maxrows;
   for (Long64_t e = 0; e < nent; ++e) {
      t->GetEntry(e);
      if (!(vz > -500.f)) continue;
      for (int i = 0; i < ns; ++i) {
         if (requireNoRing && ring[i] > 0.5f) continue;
         const std::string& s = species[i];
         const bool pos = (s == "p" || s == "pip" || s == "ep");
         const char* pat = (s == "ep" || s == "em") ? "5sig" : "3sig";
         n[i] += 1;
         i1[i] += Cut(s.c_str(), "1sig")->IsInside(p[i], b[i]);
         i2[i] += Cut(s.c_str(), "2sig")->IsInside(p[i], b[i]);
         iP[i] += Cut(s.c_str(), pat)->IsInside(p[i], b[i]);
         iPi[i] += Cut(pos ? "pip" : "pim", "3sig")->IsInside(p[i], b[i]);
         iE[i] += Cut(pos ? "ep" : "em", "5sig")->IsInside(p[i], b[i]);
         if (pos) iProt[i] += Cut("p", "3sig")->IsInside(p[i], b[i]);
      }
   }
   fprintf(gOut, "%s %s (rows read %lld%s):\n", file, tree, (long long)nent, requireNoRing ? ", tracks without ring" : "");
   for (int i = 0; i < ns; ++i)
      fprintf(gOut, "  slot %-4s entries %10.0f  own 1sig %.3f  2sig %.3f  PAT level %.3f | inside pi(same q) 3sig %.3f  e(same q) 5sig %.3f  p 3sig %.3f\n",
              slots[i].c_str(), n[i], i1[i] / n[i], i2[i] / n[i], iP[i] / n[i], iPi[i] / n[i], iE[i] / n[i], iProt[i] / n[i]);
   fflush(gOut);
   f->Close();
}
}

void verify_pid_samples()
{
   using namespace vps;
   gOut = fopen("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/pid/verify/verify_pid_samples.txt", "w");
   gC = TFile::Open(kCuts, "READ");
   Sample("/hdd1/HADES/pp45/exp/gen4/default/lep_036.root", "EpEm_ID", {"ep", "em"}, {"ep", "em"}, 0, false);
   Sample("/hdd1/HADES/pp45/exp/gen4/default/lep_036.root", "PipPimEpEm_ID", {"pip", "pim", "ep", "em"}, {"pip", "pim", "ep", "em"}, 0, false);
   Sample("/hdd2/przygoda/HADES/pp45/exp/had_034.root", "PipPim_ID", {"pip", "pim"}, {"pip", "pim"}, 5000000, true);

   // ---- simulation START check
   for (const char* fn : {"/hdd3/przygoda/HADES/files/sim/ppee_std_all.root", "/hdd3/przygoda/HADES/files/sim/ppee_lep_all.root"}) {
      TFile* f = TFile::Open(fn, "READ");
      TTree* t = (TTree*)f->Get("PipPimEpEm_ID");
      t->SetBranchStatus("*", 0);
      float isb, vz, bp, bm, est;
      for (const char* x : {"isBest", "eVertReco_z", "ep_beta", "em_beta", "evtStart"}) t->SetBranchStatus(x, 1);
      t->SetBranchAddress("isBest", &isb); t->SetBranchAddress("eVertReco_z", &vz);
      t->SetBranchAddress("ep_beta", &bp); t->SetBranchAddress("em_beta", &bm); t->SetBranchAddress("evtStart", &est);
      double n = 0, lp = 0, lm = 0, s1 = 0, minp = 9, minm = 9;
      for (Long64_t e = 0; e < t->GetEntries(); ++e) {
         t->GetEntry(e);
         if (!(isb > -0.5f && vz > -500.f)) continue;
         n += 1; lp += bp < 0.9f; lm += bm < 0.9f; s1 += std::lround(est) == 1;
         minp = std::min<double>(minp, bp); minm = std::min<double>(minm, bm);
      }
      fprintf(gOut, "%s: rows %.0f  e+ DST beta < 0.9 %.5f  e- %.5f  min beta e+ %.3f e- %.3f  evtStart==1 %.5f\n", fn, n, lp / n, lm / n, minp, minm, s1 / n);
      fflush(gOut);
      f->Close();
   }
   fclose(gOut);
}
