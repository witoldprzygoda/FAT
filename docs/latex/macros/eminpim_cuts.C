#include <ROOT/RDataFrame.hxx>
#include <TROOT.h>
#include <cstdio>
#include <string>
#include <vector>
void eminpim_cuts() {
  ROOT::EnableImplicitMT(6);
  const char* dir = "/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/";
  struct C { const char* name; const char* expr; };
  std::vector<C> cuts = {
    {"(none)", "true"},
    {"dEdx_MDC < 2.2", "pim_dedx_mdc < 2.2"},
    {"dEdx_MDC < 2.0", "pim_dedx_mdc < 2.0"},
    {"dEdx_MDC < 1.9", "pim_dedx_mdc < 1.9"},
    {"rkchi2 < 50", "pim_rkchi2 < 50"},
    {"not (0 <= oa_hadr < 2 deg)", "!(pim_oa_hadr >= 0 && pim_oa_hadr < 2)"},
    {"dEdx < 2.0 && rkchi2 < 50", "pim_dedx_mdc < 2.0 && pim_rkchi2 < 50"},
    {"dEdx < 2.0 && !(0<=oa_hadr<2)", "pim_dedx_mdc < 2.0 && !(pim_oa_hadr >= 0 && pim_oa_hadr < 2)"},
  };
  const char* pb[4] = {"pim_p_cor < 300", "pim_p_cor >= 300 && pim_p_cor < 600", "pim_p_cor >= 600", "true"};
  const char* pn[4] = {"p<300", "300-600", "p>=600", "all"};
  for (const char* set : {"std", "lep"}) {
    ROOT::RDataFrame df("scan_nt", std::string(dir) + "output_" + set + "_pippimepem.root");
    auto base = df.Filter("eVertReco_z > -160");
    auto pi = base.Filter("(int)pim_sim_id == 9");
    auto el = base.Filter("(int)pim_sim_id == 3");
    std::vector<std::vector<ROOT::RDF::RResultPtr<ULong64_t>>> np(4), ne(4);
    for (int b = 0; b < 4; ++b) for (auto& c : cuts) { np[b].push_back(pi.Filter(pb[b]).Filter(c.expr).Count()); ne[b].push_back(el.Filter(pb[b]).Filter(c.expr).Count()); }
    printf("\n=== %s set, pi- slot after PID (vertex z > -160): kept fraction pi- / e-  per momentum bin [MeV/c]; last column: e- share of the slot after the cut\n", set);
    printf("%-30s", "cut"); for (int b = 0; b < 4; ++b) printf(" %17s", pn[b]); printf(" %10s\n", "e share");
    printf("%-30s", "  (counts pi- / e-)"); for (int b = 0; b < 4; ++b) printf(" %8llu/%-8llu", *np[b][0], *ne[b][0]); printf("\n");
    for (size_t k = 0; k < cuts.size(); ++k) {
      printf("%-30s", cuts[k].name);
      for (int b = 0; b < 4; ++b) printf("    %5.3f / %5.3f", *np[b][0] ? (double)*np[b][k] / *np[b][0] : 0., *ne[b][0] ? (double)*ne[b][k] / *ne[b][0] : 0.);
      printf(" %9.2f%%\n", 100. * *ne[3][k] / (double)(*np[3][k] + *ne[3][k]));
    }
  }
}
