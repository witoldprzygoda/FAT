// rep_stats.C -- statistics of the simulation samples for the report (tables/stat_*.txt)
//
// For each set (std, lep, eta) and channel (pi+pi-e+e-, pi+pi-e+e+, pi+pi-e-e-):
//   PAT input rows (highest key cycle of the hypothesis tree), FAT rows (scan_nt: isBest >= 0, eVertReco_z > -500,
//   PT3 bit), rows with eVertReco_z > -160, owned rows (partition A) and their sum of sim_w (double),
//   sigma = sum(sim_w) / L_sim. For pi+pi-e+e- also the signal classes (sim_tmpl == 1, sim_cls 2..6).
// Weights: rules v2 of src/sim_truth.h (since 2026-11-13; the stored sim_w), PDG branching ratios, BR(eta -> pi+pi-e+e-) = 2.68e-4, BR(X17) = 0
// (the physics reweighting of the eta generator is NOT applied here).
//
// Run: root -l -b -q 'docs/latex/macros/rep_stats.C+'
#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <TROOT.h>
#include <ROOT/RDataFrame.hxx>
#include <cstdio>
#include <vector>
#include <string>

namespace rs {
const char* kIn  = "/hdd3/przygoda/HADES/files/sim/ppee_%s_all.root";
const char* kOut = "/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/output_%s_%s.root";
const double kLsim = 221.05; // nb^-1

Long64_t patRows(const char* set, const char* tree, int& cycle, int& ncycles)
{
   TFile f(Form(kIn, set));
   Long64_t n = -1;
   cycle = -1;
   ncycles = 0;
   for (auto* o : *f.GetListOfKeys()) {
      TKey* k = (TKey*)o;
      if (std::string(k->GetName()) != tree) continue;
      ++ncycles;
      if (k->GetCycle() > cycle) {
         cycle = k->GetCycle();
         n = ((TTree*)k->ReadObj())->GetEntries();
      }
   }
   return n;
}
} // namespace rs

void rep_stats()
{
   using namespace rs;
   ROOT::EnableImplicitMT(12);
   const char* sets[3] = {"std", "lep", "eta"};
   const char* ch[3] = {"pippimepem", "pippimepep", "pippimemem"};
   const char* tr[3] = {"PipPimEpEm_ID", "PipPimEpEp_ID", "PipPimEmEm_ID"};
   FILE* f = fopen("/home/przygoda/HADES/FAT/FAT/docs/latex/tables/stat_samples.txt", "w");
   fprintf(f, "# set channel PAT_rows(cycle/ncycles) FAT_rows vz>-160 owned_A(vz>-160) sumW_owned(vz>-160) sigma_nb  (L_sim = %.2f nb^-1)\n", kLsim);
   FILE* g = fopen("/home/przygoda/HADES/FAT/FAT/docs/latex/tables/stat_classes.txt", "w");
   fprintf(g, "# pi+pi-e+e-, owned A, vz > -160: rows and sum of sim_w per template / signal class; sigma = sumW / L_sim [nb]\n");
   for (int s = 0; s < 3; ++s)
      for (int c = 0; c < 3; ++c) {
         int cyc, ncyc;
         Long64_t npat = patRows(sets[s], tr[c], cyc, ncyc);
         ROOT::RDataFrame df("scan_nt", Form(kOut, sets[s], ch[c]));
         auto all = df.Count();
         auto dv = df.Filter("eVertReco_z > -160");
         auto nv = dv.Count();
         auto own = dv.Filter("sim_own_A == 1").Define("w", "(double)sim_w");
         auto no = own.Count();
         auto sw = own.Sum<double>("w");
         // per template and class (only pi+pi-e+e-)
         std::vector<ROOT::RDF::RResultPtr<ULong64_t>> nt, nc;
         std::vector<ROOT::RDF::RResultPtr<double>> wt, wc;
         if (c == 0) {
            for (int t = 1; t <= 4; ++t) {
               auto ft = own.Filter(Form("sim_tmpl == %d", t));
               nt.push_back(ft.Count());
               wt.push_back(ft.Sum<double>("w"));
            }
            for (int k = 2; k <= 6; ++k) {
               auto fk = own.Filter(Form("sim_tmpl == 1 && sim_cls == %d", k));
               nc.push_back(fk.Count());
               wc.push_back(fk.Sum<double>("w"));
            }
         }
         fprintf(f, "%-4s %-11s %10lld(;%d/%d) %10llu %10llu %10llu %14.4f %12.4f\n", sets[s], ch[c], npat, cyc, ncyc,
                 *all, *nv, *no, *sw, *sw / kLsim);
         printf("%-4s %-11s %10lld(;%d/%d) %10llu %10llu %10llu %14.4f %12.4f\n", sets[s], ch[c], npat, cyc, ncyc,
                *all, *nv, *no, *sw, *sw / kLsim);
         if (c == 0) {
            const char* tn[4] = {"signal", "phys_bkg", "fake_lepton", "diff_sources"};
            for (int t = 0; t < 4; ++t)
               fprintf(g, "%-4s tmpl %d %-13s %10llu %14.4f %12.5f\n", sets[s], t + 1, tn[t], *nt[t], *wt[t], *wt[t] / kLsim);
            const char* cn[5] = {"bump_pi0Dalitz", "bump_conv", "peak_conv", "peak_direct", "X17"};
            for (int k = 0; k < 5; ++k)
               fprintf(g, "%-4s cls %d %-15s %10llu %14.4f %12.5f\n", sets[s], k + 2, cn[k], *nc[k], *wc[k], *wc[k] / kLsim);
         }
      }
   fclose(f);
   fclose(g);
}
