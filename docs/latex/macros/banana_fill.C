// banana_fill.C -- run on vae25 (read-only on the samples), fills beta_new versus q*p
// for the slots of one hypothesis tree of a PAT sample written WITHOUT PID cut.
//   pattern : file pattern (TChain::Add)
//   tree    : hypothesis tree, e.g. PipPimEpEm (lepton sample LEP_merged) or PPipPim (hadron sample HH_elossfix)
//   slots   : comma-separated "<slot>:<sign>", e.g. "pip:+,ep:+,pim:-,em:-"
//   out     : output ROOT file
//   maxfiles: number of files to use (0 = all)
// Selection: isBest==1 && eVertReco_z>-500 (no trigger, no PID cut in the samples).
// Histograms (q*p in MeV/c, 10 MeV/c x 0.0025): hAll, hSlot_<slot>, hTrue_<species> (true GEANT id of the track:
// 14 p, 8 pi+, 9 pi-, 2 e+, 3 e-, other) and hSlotTrue_<slot>_<species>.
#include <TChain.h>
#include <TFile.h>
#include <TH2D.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TSystem.h>
#include <TString.h>
#include <TRegexp.h>
#include <algorithm>
#include <vector>
#include <cstdio>
#include <cmath>
void banana_fill(const char* pattern, const char* tree, const char* slots, const char* out, int maxfiles = 0)
{
   TChain ch(tree);
   TString pat(pattern);
   TString dir = gSystem->DirName(pattern), base = gSystem->BaseName(pattern);
   std::vector<TString> files;
   void* d = gSystem->OpenDirectory(dir);
   const char* e;
   TRegexp re(base, kTRUE);
   while ((e = gSystem->GetDirEntry(d))) { TString s(e); if (s.Index(re) == 0 && s.EndsWith(".root")) files.push_back(dir + "/" + s); }
   gSystem->FreeDirectory(d);
   std::sort(files.begin(), files.end());
   int nuse = (maxfiles > 0 && maxfiles < (int)files.size()) ? maxfiles : (int)files.size();
   for (int i = 0; i < nuse; ++i) ch.Add(files[i]);
   printf("FILES %d of %zu ENTRIES %lld\n", nuse, files.size(), ch.GetEntries());

   TObjArray* tok = TString(slots).Tokenize(",");
   const int ns = tok->GetEntries();
   std::vector<TString> sl(ns); std::vector<int> sg(ns);
   for (int i = 0; i < ns; ++i) { TString t = ((TObjString*)tok->At(i))->GetString(); sl[i] = t(0, t.Index(":")); sg[i] = t.EndsWith("-") ? -1 : +1; }
   std::vector<float> p(ns), b(ns), id(ns); float best, vz;
   ch.SetBranchStatus("*", 0);
   for (int i = 0; i < ns; ++i) {
      ch.SetBranchStatus(sl[i] + "_p", 1);        ch.SetBranchAddress(sl[i] + "_p", &p[i]);
      ch.SetBranchStatus(sl[i] + "_beta_new", 1); ch.SetBranchAddress(sl[i] + "_beta_new", &b[i]);
      ch.SetBranchStatus(sl[i] + "_sim_id", 1);   ch.SetBranchAddress(sl[i] + "_sim_id", &id[i]);
   }
   ch.SetBranchStatus("isBest", 1);      ch.SetBranchAddress("isBest", &best);
   ch.SetBranchStatus("eVertReco_z", 1); ch.SetBranchAddress("eVertReco_z", &vz);

   TFile fo(out, "RECREATE");
   const int nx = 600, ny = 520; const double x0 = -2500, x1 = 3500, y0 = 0., y1 = 1.3;
   const char* sp[6] = {"p", "pip", "pim", "ep", "em", "other"};
   TH2D* hAll = new TH2D("hAll", ";q p [MeV/c];beta_new", nx, x0, x1, ny, y0, y1);
   TH2D* hT[6]; for (int k = 0; k < 6; ++k) hT[k] = new TH2D(Form("hTrue_%s", sp[k]), "", nx, x0, x1, ny, y0, y1);
   std::vector<TH2D*> hS(ns); std::vector<std::vector<TH2D*>> hST(ns, std::vector<TH2D*>(6));
   for (int i = 0; i < ns; ++i) {
      hS[i] = new TH2D("hSlot_" + sl[i], "", nx, x0, x1, ny, y0, y1);
      for (int k = 0; k < 6; ++k) hST[i][k] = new TH2D(Form("hSlotTrue_%s_%s", sl[i].Data(), sp[k]), "", nx, x0, x1, ny, y0, y1);
   }
   Long64_t nb = 0, n = ch.GetEntries();
   for (Long64_t ev = 0; ev < n; ++ev) {
      ch.GetEntry(ev);
      if (best != 1 || vz <= -500) continue;
      ++nb;
      for (int i = 0; i < ns; ++i) {
         double x = sg[i] * p[i];
         int g = (int)lround(id[i]), k = 5;
         if (g == 14) k = 0; else if (g == 8) k = 1; else if (g == 9) k = 2; else if (g == 2) k = 3; else if (g == 3) k = 4;
         hAll->Fill(x, b[i]); hT[k]->Fill(x, b[i]); hS[i]->Fill(x, b[i]); hST[i][k]->Fill(x, b[i]);
      }
   }
   printf("BESTROWS %lld\n", nb);
   for (int i = 0; i < ns; ++i) { printf("SLOT %-4s", sl[i].Data()); for (int k = 0; k < 6; ++k) printf(" %s=%.0f", sp[k], hST[i][k]->GetEntries()); printf("\n"); }
   fo.Write(); fo.Close();
}
