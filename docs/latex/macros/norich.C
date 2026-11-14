// e- in the pi- slot (no RICH ring correlation) versus e- in the e- slot (with ring), PipPimEpEm without PID cut,
// isBest==1, eVertReco_z>-500: creation process, place, parent, momentum, theta, RICH variables.
#include <TChain.h>
#include <TSystem.h>
#include <map>
#include <string>
#include <cmath>
#include <cstdio>
void norich(const char* pattern, int nfiles) {
  TChain ch("PipPimEpEm");
  void* d = gSystem->OpenDirectory(gSystem->DirName(pattern)); const char* e; std::vector<TString> fs;
  while ((e = gSystem->GetDirEntry(d))) { TString s(e); if (s.BeginsWith("lep_sim_") && s.EndsWith(".root")) fs.push_back(TString(gSystem->DirName(pattern)) + "/" + s); }
  std::sort(fs.begin(), fs.end()); for (int i = 0; i < nfiles && i < (int)fs.size(); ++i) ch.Add(fs[i]);
  const char* sl[2] = {"pim", "em"};
  float id[2], proc[2], par[2], vx[2], vy[2], vz[2], p[2], th[2], rind[2], rmq[2], isr[2], ghost[2], med[2], kl[2];
  float best, ez;
  ch.SetBranchStatus("*", 0);
  auto B = [&](const char* n, float* a) { ch.SetBranchStatus(n, 1); ch.SetBranchAddress(n, a); };
  for (int s = 0; s < 2; ++s) {
    std::string x = sl[s];
    B((x + "_sim_id").c_str(), &id[s]); B((x + "_sim_processid").c_str(), &proc[s]); B((x + "_sim_parentid").c_str(), &par[s]);
    B((x + "_sim_vertexx").c_str(), &vx[s]); B((x + "_sim_vertexy").c_str(), &vy[s]); B((x + "_sim_vertexz").c_str(), &vz[s]);
    B((x + "_sim_p").c_str(), &p[s]); B((x + "_theta").c_str(), &th[s]); B((x + "_rich_ind").c_str(), &rind[s]);
    B((x + "_richmatchquality").c_str(), &rmq[s]); B((x + "_isring").c_str(), &isr[s]); B((x + "_sim_isghost").c_str(), &ghost[s]);
    B((x + "_sim_mediumid").c_str(), &med[s]); B((x + "_kIsLepton").c_str(), &kl[s]);
  }
  B("isBest", &best); B("eVertReco_z", &ez);
  std::map<std::string, long> c[2]; long n[2] = {0, 0};
  auto procName = [](int q) { switch (q) { case 0: return "primary"; case 5: return "decay"; case 6: return "pair (conversion)"; case 7: return "Compton";
    case 8: return "photoeffect"; case 9: return "bremsstrahlung"; case 10: return "delta ray"; case 11: return "annihilation"; case 12: return "hadronic"; default: return "other"; } };
  for (Long64_t i = 0; i < ch.GetEntries(); ++i) {
    ch.GetEntry(i);
    if (best != 1 || ez <= -500) continue;
    for (int s = 0; s < 2; ++s) {
      if ((int)lround(id[s]) != 3) continue;
      ++n[s];
      c[s][std::string("proc: ") + procName((int)lround(proc[s]))]++;
      int pp = (int)lround(par[s]);
      c[s][std::string("parent: ") + (pp == 1 ? "photon" : pp == 7 ? "pi0" : pp == 17 ? "eta" : pp == 2 || pp == 3 ? "e+-" : pp == 0 ? "none" : (pp == 8 || pp == 9) ? "pi+-" : pp == 5 || pp == 6 ? "mu" : "other")]++;
      double r = std::hypot(vx[s], vy[s]);
      const char* where = (vz[s] > -200 && vz[s] < -60 && r < 20) ? "target region" : (vz[s] < -200 ? "upstream z<-200" : (r < 20 ? "beam axis, z>-60" : "off axis r>20 mm"));
      c[s][std::string("vertex: ") + where]++;
      double zb = vz[s];
      c[s][std::string("vertex z: ") + (zb < -200 ? "<-200" : zb < -60 ? "-200..-60" : zb < 200 ? "-60..200" : zb < 500 ? "200..500" : zb < 1000 ? "500..1000" : ">1000")]++;
      double pm = p[s];
      c[s][std::string("p_sim: ") + (pm < 50 ? "<50" : pm < 100 ? "50-100" : pm < 200 ? "100-200" : pm < 400 ? "200-400" : ">400")]++;
      double t = th[s];
      c[s][std::string("theta: ") + (t < 18 ? "<18" : t < 30 ? "18-30" : t < 60 ? "30-60" : t < 85 ? "60-85" : ">85")]++;
      c[s][std::string("rich_ind: ") + (rind[s] >= 0 ? ">=0 (ring assigned)" : "<0")]++;
      c[s][std::string("isring: ") + (isr[s] > 0 ? "1" : "0")]++;
      c[s][std::string("kIsLepton: ") + (kl[s] > 0 ? "1" : "0")]++;
      c[s][std::string("isghost: ") + (ghost[s] > 0 ? "1" : "0")]++;
      c[s][std::string("medium: ") + std::to_string((int)lround(med[s]))]++;
    }
  }
  printf("RESULT e- in pim slot (no ring): %ld   e- in em slot (ring): %ld   files %d\n", n[0], n[1], nfiles);
  std::map<std::string, int> keys; for (int s = 0; s < 2; ++s) for (auto& kv : c[s]) keys[kv.first] = 1;
  for (auto& k : keys) { double a = n[0] ? 100. * c[0][k.first] / n[0] : 0, b = n[1] ? 100. * c[1][k.first] / n[1] : 0; if (a < 0.5 && b < 0.5) continue; printf("RESULT %-34s  pim-slot %6.1f %%   em-slot %6.1f %%\n", k.first.c_str(), a, b); }
}
