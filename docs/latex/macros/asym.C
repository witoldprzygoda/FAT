// Ring-less e+ and e- tracks in the hadron PID sample without PID cut (HH_elossfix, tree PipPim, ALL combinations):
// every distinct track of an event counted once. e+ / e- versus momentum, theta and creation place.
#include <TChain.h>
#include <TSystem.h>
#include <TTree.h>
#include <set>
#include <tuple>
#include <cmath>
#include <cstdio>
#include <vector>
#include <algorithm>
void asym(const char* dir, int nfiles) {
  TChain ch("PipPim");
  void* d = gSystem->OpenDirectory(dir); const char* e; std::vector<TString> fs;
  while ((e = gSystem->GetDirEntry(d))) { TString s(e); if (s.EndsWith("_hadron_out.root")) fs.push_back(TString(dir) + "/" + s); }
  std::sort(fs.begin(), fs.end()); for (int i = 0; i < nfiles && i < (int)fs.size(); ++i) ch.Add(fs[i]);
  float id[2], p[2], th[2], ph[2], sp[2], vz[2], vx[2], vy[2], proc[2], ez, ex, ey, best;
  const char* sl[2] = {"pip", "pim"};
  ch.SetBranchStatus("*", 0);
  auto B = [&](const char* n, float* a) { ch.SetBranchStatus(n, 1); ch.SetBranchAddress(n, a); };
  for (int s = 0; s < 2; ++s) { std::string x = sl[s];
    B((x + "_sim_id").c_str(), &id[s]); B((x + "_p").c_str(), &p[s]); B((x + "_theta").c_str(), &th[s]); B((x + "_phi").c_str(), &ph[s]);
    B((x + "_sim_p").c_str(), &sp[s]); B((x + "_sim_vertexz").c_str(), &vz[s]); B((x + "_sim_vertexx").c_str(), &vx[s]); B((x + "_sim_vertexy").c_str(), &vy[s]);
    B((x + "_sim_processid").c_str(), &proc[s]); }
  B("eVertReco_z", &ez); B("eVertReco_x", &ex); B("eVertReco_y", &ey); B("isBest", &best);
  // counters: [charge 0=e+ (pip slot),1=e- (pim slot)][p bin][theta bin], downstream-created and target-created
  const int NP = 5, NT = 4; double pe[NP + 1] = {0, 100, 200, 300, 500, 1e9}; double te[NT + 1] = {0, 20, 30, 50, 90};
  long cnt[2][2][NP][NT] = {}; long npi[2] = {0, 0};
  std::set<std::tuple<int, int, int>> seen[2]; float lastz = -1e9, lastx = -1e9, lasty = -1e9; long events = 0;
  for (Long64_t i = 0; i < ch.GetEntries(); ++i) {
    ch.GetEntry(i);
    if (ez != lastz || ex != lastx || ey != lasty) { seen[0].clear(); seen[1].clear(); lastz = ez; lastx = ex; lasty = ey; ++events; }
    if (ez <= -500) continue;
    for (int s = 0; s < 2; ++s) {
      auto key = std::make_tuple((int)lround(p[s] * 10), (int)lround(th[s] * 100), (int)lround(ph[s] * 100));
      if (!seen[s].insert(key).second) continue;            // track already counted in this event
      int g = (int)lround(id[s]);
      if (g == (s == 0 ? 8 : 9)) { npi[s]++; continue; }
      if (g != (s == 0 ? 2 : 3)) continue;
      double r = std::hypot(vx[s], vy[s]);
      int w = (vz[s] > -200 && vz[s] < -60 && r < 20) ? 0 : 1;   // 0 target region, 1 elsewhere
      int ip = 0; while (ip < NP - 1 && sp[s] >= pe[ip + 1]) ++ip;
      int it = 0; while (it < NT - 1 && th[s] >= te[it + 1]) ++it;
      cnt[s][w][ip][it]++;
    }
  }
  printf("RESULT events %ld ; distinct ring-less tracks: pi+ %ld, pi- %ld\n", events, npi[0], npi[1]);
  for (int w = 0; w < 2; ++w) {
    long te_ = 0, tm = 0; for (int a = 0; a < NP; ++a) for (int b = 0; b < NT; ++b) { te_ += cnt[0][w][a][b]; tm += cnt[1][w][a][b]; }
    printf("RESULT %s: e+ %ld  e- %ld  ratio e+/e- %.3f\n", w ? "created outside the target region" : "created in the target region", te_, tm, tm ? (double)te_ / tm : 0.);
    for (int a = 0; a < NP; ++a) { long x = 0, y = 0; for (int b = 0; b < NT; ++b) { x += cnt[0][w][a][b]; y += cnt[1][w][a][b]; }
      printf("RESULT    true p %4.0f-%-5.0f MeV/c: e+ %6ld e- %6ld ratio %.3f\n", pe[a], a == NP - 1 ? 9999. : pe[a + 1], x, y, y ? (double)x / y : 0.); }
    for (int b = 0; b < NT; ++b) { long x = 0, y = 0; for (int a = 0; a < NP; ++a) { x += cnt[0][w][a][b]; y += cnt[1][w][a][b]; }
      printf("RESULT    theta %2.0f-%2.0f deg:         e+ %6ld e- %6ld ratio %.3f\n", te[b], te[b + 1], x, y, y ? (double)x / y : 0.); }
  }
}
