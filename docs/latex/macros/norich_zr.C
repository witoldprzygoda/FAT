#include <TChain.h>
#include <TSystem.h>
#include <TH2D.h>
#include <cmath>
#include <cstdio>
#include <vector>
#include <algorithm>
void norich_zr(const char* dir, int nfiles) {
  TChain ch("PipPimEpEm");
  void* d = gSystem->OpenDirectory(dir); const char* e; std::vector<TString> fs;
  while ((e = gSystem->GetDirEntry(d))) { TString s(e); if (s.BeginsWith("lep_sim_") && s.EndsWith(".root")) fs.push_back(TString(dir) + "/" + s); }
  std::sort(fs.begin(), fs.end()); for (int i = 0; i < nfiles && i < (int)fs.size(); ++i) ch.Add(fs[i]);
  float id, isr, vx, vy, vz, best, ez, med;
  ch.SetBranchStatus("*", 0);
  auto B = [&](const char* n, float* a) { ch.SetBranchStatus(n, 1); ch.SetBranchAddress(n, a); };
  B("pim_sim_id", &id); B("pim_isring", &isr); B("pim_sim_vertexx", &vx); B("pim_sim_vertexy", &vy); B("pim_sim_vertexz", &vz);
  B("pim_sim_mediumid", &med); B("isBest", &best); B("eVertReco_z", &ez);
  TH2D h("h", "", 40, -200, 800, 16, 0, 800);
  std::vector<double> zc; for (int k = 0; k < 41; ++k) zc.push_back(-1000 + 25 * k);
  std::vector<long> nsph(zc.size(), 0); long n = 0;
  long nmed[64] = {0}; double zm[64] = {0}, rm[64] = {0};
  for (Long64_t i = 0; i < ch.GetEntries(); ++i) { ch.GetEntry(i);
    if (best != 1 || ez <= -500 || (int)lround(id) != 3 || isr > 0) continue;
    double r = std::hypot(vx, vy); ++n; h.Fill(vz, r);
    int m = (int)lround(med); if (m >= 0 && m < 64) { nmed[m]++; zm[m] += vz; rm[m] += r; }
    for (size_t k = 0; k < zc.size(); ++k) { double dd = std::hypot(r, vz - zc[k]); if (std::fabs(dd - 872.) < 6.) nsph[k]++; } }
  printf("RESULT e- without ring in the pi- slot: %ld\n", n);
  size_t kb = 0; for (size_t k = 0; k < zc.size(); ++k) if (nsph[k] > nsph[kb]) kb = k;
  printf("RESULT sphere of radius 872 mm (+-6): best centre z = %.0f mm with %.1f %% of the vertices (next: %.0f -> %.1f %%, %.0f -> %.1f %%)\n", zc[kb], 100. * nsph[kb] / n, zc[kb > 0 ? kb - 1 : kb], 100. * nsph[kb > 0 ? kb - 1 : kb] / n, zc[kb + 1 < zc.size() ? kb + 1 : kb], 100. * nsph[kb + 1 < zc.size() ? kb + 1 : kb] / n);
  for (int m = 0; m < 64; ++m) if (nmed[m] > 0.01 * n) printf("RESULT medium %2d: %5.1f %%  mean z %6.0f mm  mean r %5.0f mm\n", m, 100. * nmed[m] / n, zm[m] / nmed[m], rm[m] / nmed[m]);
  std::vector<std::pair<double, int>> cells; for (int ix = 1; ix <= 40; ++ix) for (int iy = 1; iy <= 16; ++iy) cells.push_back({h.GetBinContent(ix, iy), h.GetBin(ix, iy)});
  std::sort(cells.rbegin(), cells.rend());
  for (int k = 0; k < 12; ++k) { int ix, iy, iz; h.GetBinXYZ(cells[k].second, ix, iy, iz); printf("RESULT cell z %4.0f..%4.0f  r %3.0f..%3.0f : %5.1f %%\n", h.GetXaxis()->GetBinLowEdge(ix), h.GetXaxis()->GetBinUpEdge(ix), h.GetYaxis()->GetBinLowEdge(iy), h.GetYaxis()->GetBinUpEdge(iy), 100. * cells[k].first / n); }
}
