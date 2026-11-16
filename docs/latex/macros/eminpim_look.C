#include <ROOT/RDataFrame.hxx>
#include <TROOT.h>
#include <cstdio>
#include <string>
#include <vector>
#include <algorithm>
void eminpim_look() {
  ROOT::EnableImplicitMT(6);
  const char* dir = "/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/";
  for (const char* set : {"std", "lep"}) {
    ROOT::RDataFrame df("scan_nt", std::string(dir) + "output_" + set + "_pippimepem.root");
    auto base = df.Filter("eVertReco_z > -160").Define("vr", "sqrt(pim_sim_vertexx*pim_sim_vertexx+pim_sim_vertexy*pim_sim_vertexy)");
    auto el = base.Filter("(int)pim_sim_id == 3");
    auto pi = base.Filter("(int)pim_sim_id == 9 && pim_p_cor >= 200");
    auto nE = el.Count();
    auto tgt = el.Filter("pim_sim_vertexz > -200 && pim_sim_vertexz < -60 && vr < 20").Count();
    auto dn = el.Filter("pim_sim_vertexz > 100").Count();
    auto conv = el.Filter("(int)pim_sim_processid == 6").Count();
    auto comp = el.Filter("(int)pim_sim_processid == 7").Count();
    auto dec = el.Filter("(int)pim_sim_processid == 5").Count();
    auto gh = el.Filter("pim_sim_isghost > 0").Count();
    auto ring = el.Filter("pim_isring > 0").Count();
    auto pring = pi.Filter("pim_isring > 0").Count();
    auto nP = pi.Count();
    printf("\n=== %s: e- in pi- slot after PID: %llu ; created in target region %.1f %%, at z > 100 mm %.1f %% ; conversion %.1f %%, Compton %.1f %%, decay %.1f %% ; ghost %.1f %% ; isring %.1f %% (pi- p>=200: isring %.1f %%)\n",
      set, *nE, 100.**tgt / *nE, 100.**dn / *nE, 100.**conv / *nE, 100.**comp / *nE, 100.**dec / *nE, 100.**gh / *nE, 100.**ring / *nE, 100.**pring / *nP);
    for (const char* v : {"pim_dedx_mdc", "pim_oa_hadr", "pim_oa_lept", "pim_rkchi2", "pim_p_cor", "pim_theta"}) {
      auto ve = el.Take<float>(v); auto vp = pi.Take<float>(v);
      std::vector<float> a(ve->begin(), ve->end()), b(vp->begin(), vp->end());
      std::sort(a.begin(), a.end()); std::sort(b.begin(), b.end());
      auto q = [](std::vector<float>& x, double f) { return x.empty() ? -999.f : x[(size_t)(f * (x.size() - 1))]; };
      printf("  %-14s e-: q10 %8.2f q50 %8.2f q90 %8.2f   | pi- (p>=200): q10 %8.2f q50 %8.2f q90 %8.2f\n", v, q(a, .1), q(a, .5), q(a, .9), q(b, .1), q(b, .5), q(b, .9));
    }
  }
}
