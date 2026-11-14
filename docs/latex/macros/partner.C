// Conversion pairs among the ring-less tracks of the hadron PID sample (PipPim, all combinations): for every distinct
// ring-less e- / e+ from a photon conversion, is the partner of the other charge (same GEANT creation vertex) also
// reconstructed as a ring-less track in the event? Ghost fractions by charge.
#include <TChain.h>
#include <TSystem.h>
#include <map>
#include <set>
#include <tuple>
#include <cmath>
#include <cstdio>
#include <vector>
#include <algorithm>
struct T { int q; float vx, vy, vz, th; int ghost; };
void partner(const char* dir, int nfiles) {
  TChain ch("PipPim");
  void* d = gSystem->OpenDirectory(dir); const char* e; std::vector<TString> fs;
  while ((e = gSystem->GetDirEntry(d))) { TString s(e); if (s.EndsWith("_hadron_out.root")) fs.push_back(TString(dir) + "/" + s); }
  std::sort(fs.begin(), fs.end()); for (int i = 0; i < nfiles && i < (int)fs.size(); ++i) ch.Add(fs[i]);
  float id[2], p[2], th[2], ph[2], vz[2], vx[2], vy[2], proc[2], gh[2], ez, ex, ey;
  const char* sl[2] = {"pip", "pim"};
  ch.SetBranchStatus("*", 0);
  auto B = [&](const char* n, float* a) { ch.SetBranchStatus(n, 1); ch.SetBranchAddress(n, a); };
  for (int s = 0; s < 2; ++s) { std::string x = sl[s];
    B((x + "_sim_id").c_str(), &id[s]); B((x + "_p").c_str(), &p[s]); B((x + "_theta").c_str(), &th[s]); B((x + "_phi").c_str(), &ph[s]);
    B((x + "_sim_vertexz").c_str(), &vz[s]); B((x + "_sim_vertexx").c_str(), &vx[s]); B((x + "_sim_vertexy").c_str(), &vy[s]);
    B((x + "_sim_processid").c_str(), &proc[s]); B((x + "_sim_isghost").c_str(), &gh[s]); }
  B("eVertReco_z", &ez); B("eVertReco_x", &ex); B("eVertReco_y", &ey);
  std::set<std::tuple<int, int, int, int>> seen; std::vector<T> trk;
  long nconv[2] = {0, 0}, withPartner[2] = {0, 0}, ghost[2] = {0, 0}, partnerGhost[2] = {0, 0}, nall[2] = {0, 0}, ghall[2] = {0, 0};
  auto flush = [&]() {
    for (auto& a : trk) {
      ++nall[a.q]; if (a.ghost) ++ghall[a.q];
      ++nconv[a.q]; if (a.ghost) ++ghost[a.q];
      for (auto& b : trk) if (b.q != a.q && std::fabs(b.vz - a.vz) < 0.1 && std::fabs(b.vx - a.vx) < 0.1 && std::fabs(b.vy - a.vy) < 0.1) { ++withPartner[a.q]; if (b.ghost) ++partnerGhost[a.q]; break; }
    }
    trk.clear(); seen.clear(); };
  float lz = -1e9, lx = -1e9, ly = -1e9;
  for (Long64_t i = 0; i < ch.GetEntries(); ++i) {
    ch.GetEntry(i);
    if (ez != lz || ex != lx || ey != ly) { flush(); lz = ez; lx = ex; ly = ey; }
    if (ez <= -500) continue;
    for (int s = 0; s < 2; ++s) {
      int g = (int)lround(id[s]); int q = (g == 2) ? 0 : (g == 3 ? 1 : -1);
      if (q < 0 || (int)lround(proc[s]) != 6) continue;
      auto key = std::make_tuple(s, (int)lround(p[s] * 10), (int)lround(th[s] * 100), (int)lround(ph[s] * 100));
      if (!seen.insert(key).second) continue;
      trk.push_back({q, vx[s], vy[s], vz[s], th[s], gh[s] > 0});
    }
  }
  flush();
  for (int q = 0; q < 2; ++q)
    printf("RESULT %s from conversions (ring-less): %ld ; ghost %.1f %% ; partner of the other charge also reconstructed: %.1f %% (partner ghost in %.1f %% of those)\n",
      q ? "e-" : "e+", nconv[q], 100. * ghost[q] / nconv[q], 100. * withPartner[q] / nconv[q], withPartner[q] ? 100. * partnerGhost[q] / withPartner[q] : 0.);
}
